#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <cmath>

void setInlets(float sigma, float u_lid, int it, float uTop[][100], float uBot[][100], float uLeft[][100],
               float uRight[][100], int height, int lenght)
{
    const float ret = 1.0 - std::exp(-static_cast<double>(it * it) / (2.0 * sigma * sigma));
    for (int j = 0; j < lenght; j++)
    {
        uTop[0][j] = u_lid * ret;
        uBot[0][j] = 0.0;
    }
    for (int i = 0; i < height; i++)
    {
        uLeft[1][i] = 0.0;
        uRight[1][i] = 0.0;
    }
}
// u [coordinata x o y][posizione della cella delle x][posizione della cella delle y]
void computeEquilibrium(float u[][100][100], float rho[][100], float fEq[][100][100], int height, int lenght, int number_of_directions, float velocities[][2], float weigths[])
{
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < lenght; j++)
        {
            for (int k = 0; k < number_of_directions; k++)
            {
                float temp1 = 1.5 * (u[0][i][j] * u[0][i][j] + u[1][i][j] * u[1][i][j]);
                float temp2 = 3.0 * (u[0][i][j] * velocities[k][0] + u[1][i][j] * velocities[k][1]);
                fEq[k][i][j] = weigths[k] * rho[i][j] * (1.0 + temp2 + 0.5 * temp2 * temp2 - temp1);
            }
        }
    }
}

void computeMacroscopic(float f[][100][100], float u[][100][100], float rho[][100], int height, int lenght, int number_of_directions, float velocities[][2])
{
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < lenght; j++)
        {
            rho[i][j] = 0;
            u[0][i][j] = 0;
            u[1][i][j] = 0;
            for (int k = 0; k < number_of_directions; k++)
            {
                rho[i][j] += f[k][i][j];
                u[0][i][j] += velocities[k][0] * f[k][i][j];
                u[1][i][j] += velocities[k][1] * f[k][i][j];
            }
            u[0][i][j] /= rho[i][j];
            u[1][i][j] /= rho[i][j];
        }
    }
}
void collision_and_streaming(float om_p, float om_m, int height, int lenght, int number_of_directions, const float fEq[][100][100], float f[][100][100], float fNew[][100][100], int opposites[])
{
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < lenght; j++)
        {
            fNew[0][i][j] = (1.0 - om_p) * f[0][i][j] + om_p * fEq[0][i][j];
            f[0][i][j] = fNew[0][i][j];
        }
    }

    // collide for other indices
    for (int k = 0; k < number_of_directions; k++)
    {
        auto kOpposite = opposites[k];
        for (int i = 0; i < height; ++i)
        {
            for (int j = 0; j < lenght; ++j)
            {
                fNew[k][i][j] = (1.0 - 0.5 * (om_p + om_m)) * f[k][i][j] - 0.5 * (om_p - om_m) * f[kOpposite][i][j] +
                                0.5 * (om_p + om_m) * fEq[k][i][j] + 0.5 * (om_p - om_m) * fEq[kOpposite][i][j];
            }
        }
    }
    // stream
    for (int i = 0; i < height - 1; i++)
    {
        for (int j = 0; j < lenght - 1; j++)
        {
            f[1][i][j + 1] = fNew[1][i][j];
            f[2][i][j + 1] = fNew[2][i][j];
            f[3][i][j] = fNew[3][i + 1][j];
            f[4][i + 1][j] = fNew[4][i][j];
            f[5][i][j + 1] = fNew[5][i + 1][j];
            f[6][i][j] = fNew[6][i + 1][j + 1];
            f[7][i + 1][j] = fNew[7][i][j + 1];
            f[8][i + 1][j + 1] = fNew[8][i][j];
        }
    }
    for (int i = 0; i < height; i++)
    {
        f[1][i][lenght-1] = fNew[1][i][lenght-1];     //! TO CHECK
        f[3][i][lenght-1] = fNew[3][i][lenght-1]; //! TO CHECK
    }
    for (int j = 0; j < lenght; j++) //! TO CHECK
    {
        f[2][height-1][j] = fNew[2][height-1][j]; //! TO CHECK
        f[4][height-1][j] = fNew[4][height-1][j]; //! TO CHECK
    }
}

void zou_he_bottom_left_corner_velocity(float f[][100][100], float u[][100][100], float rho[][100], int height)
{
    u[0][height - 1][0] = u[0][height - 2][0]; // forse questa va messa come la uWall
    u[1][height - 1][0] = u[1][height - 2][0]; // idem

    rho[height - 1][0] = rho[height - 2][0];

    f[1][height - 1][0] = f[3][height - 1][0] + (2.0 / 3.0) * rho[height - 1][0] * u[0][height - 1][0];

    f[2][height - 1][0] = f[4][height - 1][0] + (2.0 / 3.0) * rho[height - 1][0] * u[1][height - 1][0];

    f[5][height - 1][0] = f[7][height - 1][0] + (1.0 / 6.0) * rho[height - 1][0] * u[0][height - 1][0] + (1.0 / 6.0) * rho[height - 1][0] * u[1][height - 1][0];

    f[6][height - 1][0] = 0.0;
    f[8][height - 1][0] = 0.0;

    f[0][height - 1][0] = rho[height - 1][0] - f[1][height - 1][0] - f[2][height - 1][0] - f[3][height - 1][0] - f[4][height - 1][0] - f[5][height - 1][0] - f[6][height - 1][0] -
                          f[7][height - 1][0] - f[8][height - 1][0];
}

void zou_he_bottom_right_corner_velocity(float f[][100][100], float u[][100][100], float rho[][100], int height, int lenght)
{
    u[0][height - 1][lenght - 1] = u[0][height - 2][lenght - 2]; // forse questa va messa come la uWall
    u[1][height - 1][lenght - 1] = u[1][height - 2][lenght - 2]; // idem

    rho[height - 1][lenght - 1] = rho[height - 2][lenght - 2];

    f[3][height - 1][lenght - 1] = f[1][height - 1][lenght - 1] - (2.0 / 3.0) * rho[height - 1][lenght - 1] * u[0][height - 1][lenght - 1];

    f[2][height - 1][lenght - 1] = f[4][height - 1][lenght - 1] + (2.0 / 3.0) * rho[height - 1][lenght - 1] * u[1][height - 1][lenght - 1];

    f[6][height - 1][lenght - 1] = f[8][height - 1][lenght - 1] - (1.0 / 6.0) * rho[height - 1][lenght - 1] * u[0][height - 1][lenght - 1] + (1.0 / 6.0) * rho[height - 1][lenght - 1] * u[1][height - 1][lenght - 1];

    f[5][height - 1][lenght - 1] = 0.0;
    f[7][height - 1][lenght - 1] = 0.0;

    f[0][height - 1][lenght - 1] = rho[height - 1][lenght - 1] - f[1][height - 1][lenght - 1] - f[2][height - 1][lenght - 1] - f[3][height - 1][lenght - 1] - f[4][height - 1][lenght - 1] - f[5][height - 1][lenght - 1] - f[6][height - 1][lenght - 1] -
                                   f[7][height - 1][lenght - 1] - f[8][height - 1][lenght - 1];
}

void zou_he_top_left_corner_velocity(float f[][100][100], float u[][100][100], float rho[][100])
{
    u[0][0][0] = u[0][1][0]; // forse questa va messa come la uWall
    u[1][0][0] = u[1][1][0]; // idem

    rho[0][0] = rho[1][0];

    f[1][0][0] = f[3][0][0] + (2.0 / 3.0) * rho[0][0] * u[0][0][0];

    f[4][0][0] = f[2][0][0] - (2.0 / 3.0) * rho[0][0] * u[1][0][0];

    f[8][0][0] = f[6][0][0] - (1.0 / 6.0) * rho[0][0] * u[0][0][0] + (1.0 / 6.0) * rho[0][0] * u[1][0][0];

    f[5][0][0] = 0.0;
    f[7][0][0] = 0.0;

    f[0][0][0] = rho[0][0] - f[1][0][0] - f[2][0][0] - f[3][0][0] - f[4][0][0] - f[5][0][0] - f[6][0][0] - f[7][0][0] -
                 f[8][0][0];
}

void zou_he_top_right_corner_velocity(float f[][100][100], float u[][100][100], float rho[][100], int lenght)
{

    u[0][0][lenght - 1] = u[0][0][lenght - 2]; // forse questa va messa come la uWall
    u[1][0][lenght - 1] = u[1][0][lenght - 2]; // idem

    rho[0][lenght - 1] = rho[0][lenght - 2];

    f[3][0][lenght - 1] = f[1][0][lenght - 1] - (2.0 / 3.0) * rho[0][lenght - 1] * u[0][0][lenght - 1];
    f[4][0][lenght - 1] = f[2][0][lenght - 1] - (2.0 / 3.0) * rho[0][lenght - 1] * u[1][0][lenght - 1];
    f[7][0][lenght - 1] = f[5][0][lenght - 1] - (1.0 / 6.0) * rho[0][lenght - 1] * u[0][0][lenght - 1] - (1.0 / 6.0) * rho[0][lenght - 1] * u[1][0][lenght - 1];

    f[6][0][lenght - 1] = 0.0;
    f[8][0][lenght - 1] = 0.0;

    f[0][0][lenght - 1] = rho[0][lenght - 1] - f[1][0][lenght - 1] - f[2][0][lenght - 1] - f[3][0][lenght - 1] - f[4][0][lenght - 1] - f[5][0][lenght - 1] - f[6][0][lenght - 1] -
                          f[7][0][lenght - 1] - f[8][0][lenght - 1];
}

void zou_he_left_wall_velocity(float f[][100][100], float u[][100][100], float rho[][100], float uLeft[][100], int height)
{
    for (int i = 0; i < height; i++)
    {
        u[0][i][0] = uLeft[0][i];
        u[1][i][0] = uLeft[1][i];

        rho[i][0] = (f[0][i][0] + f[2][i][0] + f[4][i][0] + 2.0 * (f[3][i][0] + f[6][i][0] + f[7][i][0])) / (1.0 - u[0][i][0]);

        f[1][i][0] = f[3][i][0] - 2.0 / 3.0 * rho[i][0] * u[0][i][0];

        f[5][i][0] = f[7][i][0] - 0.5 * (f[2][i][0] - f[4][i][0]) + 1.0 / 6.0 * rho[i][0] * u[0][i][0] + 0.5 * rho[i][0] * u[1][i][0];

        f[8][i][0] = f[6][i][0] + 0.5 * (f[2][i][0] - f[4][i][0]) + 1.0 / 6.0 * rho[i][0] * u[0][i][0] - 0.5 * rho[i][0] * u[1][i][0];
    }
}

void zou_he_right_wall_velocity(float f[][100][100], float u[][100][100], float rho[][100], float uRight[][100], int height, int lenght)
{
    for (int i = 0; i < height; i++)
    {
        u[0][i][lenght - 1] = uRight[0][i];
        u[1][i][lenght - 1] = uRight[1][i];

        rho[i][lenght - 1] = (f[0][i][lenght - 1] + f[2][i][lenght - 1] + f[4][i][lenght - 1] + 2.0 * (f[1][i][lenght - 1] + f[5][i][lenght - 1] + f[8][i][lenght - 1])) /
                             (1.0 + u[0][i][lenght - 1]);

        f[3][i][lenght - 1] = f[1][i][lenght - 1] - 2.0 / 3.0 * rho[i][lenght - 1] * u[0][i][lenght - 1];

        f[6][i][lenght - 1] = f[8][i][lenght - 1] - 0.5 * (f[2][i][lenght - 1] - f[4][i][lenght - 1]) - 1.0 / 6.0 * rho[i][lenght - 1] * u[0][i][lenght - 1] +
                              0.5 * rho[i][lenght - 1] * u[1][i][lenght - 1];

        f[7][i][lenght - 1] = f[5][i][lenght - 1] + 0.5 * (f[2][i][lenght - 1] - f[4][i][lenght - 1]) - 1.0 / 6.0 * rho[i][lenght - 1] * u[0][i][lenght - 1] - 0.5 * rho[i][lenght - 1] * u[1][i][lenght - 1];
    }
}
void zou_he_top_wall_velocity(float f[][100][100], float u[][100][100], float rho[][100], float uTop[][100],
                              int lenght)
{
    for (int j = 0; j < lenght; j++)
    {
        u[0][0][j] = uTop[0][j];
        u[1][0][j] = uTop[1][j];

        rho[0][j] = (f[0][0][j] + f[1][0][j] + f[3][0][j] + 2.0 * (f[2][0][j] + f[5][0][j] + f[6][0][j])) / (1.0 + u[1][0][j]);

        f[4][0][j] = f[2][0][j] - 2.0 / 3.0 * rho[0][j] * u[1][0][j];

        f[7][0][j] = f[5][0][j] + 0.5 * (f[1][0][j] - f[3][0][j]) - 1.0 / 6.0 * rho[0][j] * u[1][0][j] - 0.5 * rho[0][j] * u[0][0][j];

        f[8][0][j] = f[6][0][j] - 0.5 * (f[1][0][j] - f[3][0][j]) - 1.0 / 6.0 * rho[0][j] * u[1][0][j] + 0.5 * rho[0][j] * u[0][0][j];
    }
}

void zou_he_bottom_wall_velocity(float f[][100][100], float u[][100][100], float rho[][100], float uBot[][100],
                                 int height, int lenght)
{
    for (int j = 0; j < lenght; j++)
    {
        u[0][height - 1][j] = uBot[0][j];
        u[1][height - 1][j] = uBot[1][j];

        rho[height - 1][j] = (f[0][height - 1][j] + f[1][height - 1][j] + f[3][height - 1][j] + 2.0 * (f[4][height - 1][j] + f[7][height - 1][j] + f[8][height - 1][j])) /
                             (1.0 - u[1][height - 1][j]);

        f[2][height - 1][j] = f[4][height - 1][j] + 2.0 / 3.0 * rho[height - 1][j] * u[1][height - 1][j];

        f[5][height - 1][j] = f[7][height - 1][j] - 0.5 * (f[1][height - 1][j] - f[3][height - 1][j]) + 1.0 / 6.0 * rho[height - 1][j] * u[1][height - 1][j] +
                              0.5 * rho[height - 1][j] * u[0][height - 1][j];

        f[6][height - 1][j] = f[8][height - 1][j] + 0.5 * (f[1][height - 1][j] - f[3][height - 1][j]) + 1.0 / 6.0 * rho[height - 1][j] * u[1][height - 1][j] -
                              0.5 * rho[height - 1][j] * u[0][height - 1][j];
    }
}

int main()
{

    int lenght = 100;
    int height = 100;
    int npts = 100;

    float re_lbm = 100.0f; // reynolds number

    float rho_lbm = 1.0; // densità
    float t_max = 20.0;
    float x_min = 0.0;
    float x_max = 1.0;
    float y_min = 0.0;
    float y_max = 1.0;
    float c_s = 1.0f / std::sqrt(3.0f);
    int number_of_directions = 9; // 9 direzioni di velocità
                                  // # D2Q9 Velocities

                                  
    float u_lbm = 0.2; // Define the variable "u_lbm"
    float N = 100.0f;  // Define the variable "N"
    
    int ny = 100;

    float nu_lbm = u_lbm * N / re_lbm; // 0.2

    
    float tau_lbm = 0.5 + nu_lbm / (c_s * c_s); // 1.1
    float dt = re_lbm * nu_lbm / (npts * npts); // 0.002

      
    int nx = (int)std::round(ny * (x_max - x_min) / (y_max - y_min)); // 100: lattice width

    int it_max = (int)std::round(t_max / dt); // 10'000
    float sigma = 10.0 * nx;                  // 1000

    float tau_p_lbm = tau_lbm;
    float lambda_trt = 1.0 / 4.0; // Best for stability
    float tau_m_lbm = lambda_trt / (tau_p_lbm - 0.5) + 0.5;
    float om_p = 1.0 / tau_p_lbm;
    float om_m = 1.0 / tau_m_lbm;
    float u[2][100][100];



    float velocities[number_of_directions][2] = {
        // le nostre velocities
        {0., 0.},
        {1., 0.},
        {0., -1.},
        {-1., 0.},
        {0., 1.},
        {1., -1.},
        {-1., -1.},
        {-1., 1.},
        {1., 1.},
    };
    // # Weights
    float weigths[number_of_directions] = {
        4. / 9.,
        1. / 9.,
        1. / 9.,
        1. / 9.,
        1. / 9.,
        1. / 36.,
        1. / 36.,
        1. / 36.,
        1. / 36.,
    };
    // # Array for bounce-back
    int opposites[number_of_directions] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

  
    for (int a = 0; a < 2; a++)
    {
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < lenght; j++)
            {
                u[a][i][j] = 0;
            }
        }
    }

    float f[number_of_directions][100][100];
    float fEq[number_of_directions][100 ][100];
    float fNew[number_of_directions][100 ][100];
    float rho[100][100];
    float uTop[2][100];
    float uBot[2][100];
    float uLeft[2][100];
    float uRight[2][100];


    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < lenght; j++)
        {
            for (int k = 0; k < number_of_directions; k++)
            {
                fNew[k][i][j] = 0;
                fEq[k][i][j] = 0;
                f[k][i][j] = rho_lbm;
                rho[i][j] = 1.0;
            }
        }
    }

    for (int a = 0; a < 2; a++)
    {
        for (int j = 0; j < height; j++)
        {
            uTop[a][j] = 0;
            uBot[a][j] = 0;
            uLeft[a][j] = 0;
            uRight[a][j] = 0;
        }
    }

    // Initialize and compute first equilibrium
    setInlets(sigma, u_lbm, 0, uTop, uBot, uLeft, uRight, height, lenght); // 0 perchè è la prima iterazione (in realtà inizializzazione)
    computeEquilibrium(u, rho, fEq, height, lenght, number_of_directions, velocities, weigths);

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < height; j++)
        {
            for (int k = 0; k < number_of_directions; k++)
            {
                f[k][i][j] = fEq[k][i][j];
            }
        }
    }
    // open output file
    std::ofstream fout("output.txt", std::ios::out);

    // write shape
    fout << nx << " " << ny << '\n';

    for (int it = 0; it < it_max + 1; ++it)
    {
        if (it % 100 == 0)
        {
            // save data
            fout << it << '\n';
            for (int j = 0; j < lenght; j++)
            {
                for (int i = 0; i < height; i++)
                {
                    fout << rho[i][j] << ' ';
                }
            }
            fout << '\n';
            for (int j = 0; j < lenght; j++)
            {
                for (int i = 0; i < height; i++)
                {
                    fout << u[0][i][j] << ' ';
                }
            }
            fout << '\n';
            for (int j = 0; j < lenght; j++)
            {
                for (int i = 0; i < height; i++)
                {
                    fout << u[1][i][j] << ' ';
                }
            }
            fout << '\n';
            std::cout << "it = " << it << '\n';
        }

        // 1. Set inlets
        setInlets(sigma, u_lbm, it, uTop, uBot, uLeft, uRight, height, lenght);
        // 2. Compute macroscopic fields
        computeMacroscopic(f, u, rho, height, lenght, number_of_directions, velocities);
        computeEquilibrium(u, rho, fEq, height, lenght, number_of_directions, velocities, weigths);
        // 5. Streaming
        collision_and_streaming(om_p, om_m, height, lenght, number_of_directions, fEq, f, fNew, opposites);
        // 6. Boundary conditions
        zou_he_bottom_wall_velocity(f, u, rho, uBot, height, lenght);
        zou_he_left_wall_velocity(f, u, rho, uLeft, height);
        zou_he_right_wall_velocity(f, u, rho, uRight, height, lenght);
        zou_he_top_wall_velocity(f, u, rho, uTop, lenght);
        zou_he_bottom_left_corner_velocity(f, u, rho, height);
        zou_he_top_left_corner_velocity(f, u, rho);
        zou_he_top_right_corner_velocity(f, u, rho, lenght);
        zou_he_bottom_right_corner_velocity(f, u, rho, height, lenght);
        // TBD: Compute observables (drag, lift, etc)
    }

    // close output file
    fout.close();

    return 0;
}
