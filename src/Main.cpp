#include <bits/stdc++.h>
#define DIM 50

void setInlets(float u_lid, float u[][DIM][DIM], int height, int lenght)
{
    // top wall
    for (int i = 0; i < lenght; i++)
    {
        u[0][0][i] = u_lid;
        u[1][0][i] = 0;
    }
    // bottom wall
    for (int i = 1; i < lenght - 1; i++)
    {
        u[0][height - 1][i] = 0;
        u[1][height - 1][i] = 0;
    }
    // left wall
    for (int i = 1; i < height - 1; i++)
    {
        u[0][i][0] = 0;
        u[1][i][0] = 0;
    }
    // right wall
    for (int i = 1; i < height - 1; i++)
    {
        u[0][i][lenght - 1] = 0;
        u[1][i][lenght - 1] = 0;
    }
}

void computeEquilibrium(float u[][DIM][DIM], float rho[][DIM], float fEq[][DIM][DIM], int height, int lenght,
                        int number_of_directions, float velocities[][2], float weigths[])
{
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < lenght; j++)
        {
            float temp1 = 1.5 * (u[0][i][j] * u[0][i][j] + u[1][i][j] * u[1][i][j]);
            for (int k = 0; k < number_of_directions; k++)
            {
                float temp2 = 3.0 * (u[0][i][j] * velocities[k][0] + u[1][i][j] * velocities[k][1]);
                fEq[k][i][j] = weigths[k] * rho[i][j] * (1.0 + temp2 + 0.5 * temp2 * temp2 - temp1);
            }
        }
    }
}

void computeMacroscopic(float f[][DIM][DIM], float u[][DIM][DIM], float rho[][DIM], int height, int lenght,
                        int number_of_directions, float velocities[][2])
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

void collision_and_streaming(float om_p, float om_m, int height, int lenght, int number_of_directions,
                             const float fEq[][DIM][DIM], float f[][DIM][DIM], float fNew[][DIM][DIM], int opposites[])
{
    // collide and stream for index 0
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < lenght; j++)
        {
            fNew[0][i][j] = (1.0 - om_p) * f[0][i][j] + om_p * fEq[0][i][j];
            f[0][i][j] = fNew[0][i][j];
        }
    }

    // collide for other indices
    for (int k = 1; k < number_of_directions; k++)
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
    // stream for other indices
    for (int i = 0; i < height - 1; i++)
    {
        for (int j = 0; j < lenght - 1; j++)
        {
            f[1][i][j + 1] = fNew[1][i][j];
            f[2][i][j] = fNew[2][i + 1][j];
            f[3][i][j] = fNew[3][i][j + 1];
            f[4][i + 1][j] = fNew[4][i][j];
            f[5][i][j + 1] = fNew[5][i + 1][j];
            f[6][i][j] = fNew[6][i + 1][j + 1];
            f[7][i + 1][j] = fNew[7][i][j + 1];
            f[8][i + 1][j + 1] = fNew[8][i][j];
        }
    }
}

void zou_he_bottom_left_corner_velocity(float f[][DIM][DIM], float u[][DIM][DIM], float rho[][DIM], int height)
{

    f[1][height - 1][0] = f[3][height - 1][0] + (2.0 / 3.0) * rho[height - 1][0] * u[0][height - 1][0];
    f[2][height - 1][0] = f[4][height - 1][0] + (2.0 / 3.0) * rho[height - 1][0] * u[1][height - 1][0];
    f[5][height - 1][0] = f[7][height - 1][0] + (1.0 / 6.0) * rho[height - 1][0] * u[0][height - 1][0] + (1.0 / 6.0) * rho[height - 1][0] * u[1][height - 1][0];

    f[6][height - 1][0] = 0.0;
    f[8][height - 1][0] = 0.0;

    f[0][height - 1][0] = rho[height - 1][0] - f[1][height - 1][0] - f[2][height - 1][0] - f[3][height - 1][0] -
                          f[4][height - 1][0] - f[5][height - 1][0] - f[6][height - 1][0] - f[7][height - 1][0] -
                          f[8][height - 1][0];
}

void zou_he_bottom_right_corner_velocity(float f[][DIM][DIM], float u[][DIM][DIM], float rho[][DIM], int height,
                                         int lenght)
{

    f[3][height - 1][lenght - 1] = f[1][height - 1][lenght - 1] - (2.0 / 3.0) * rho[height - 1][lenght - 1] * u[0][height - 1][lenght - 1];
    f[2][height - 1][lenght - 1] = f[4][height - 1][lenght - 1] + (2.0 / 3.0) * rho[height - 1][lenght - 1] * u[1][height - 1][lenght - 1];
    f[6][height - 1][lenght - 1] = f[8][height - 1][lenght - 1] - (1.0 / 6.0) * rho[height - 1][lenght - 1] * u[0][height - 1][lenght - 1] + (1.0 / 6.0) * rho[height - 1][lenght - 1] * u[1][height - 1][lenght - 1];

    f[5][height - 1][lenght - 1] = 0.0;
    f[7][height - 1][lenght - 1] = 0.0;

    f[0][height - 1][lenght - 1] =
        rho[height - 1][lenght - 1] - f[1][height - 1][lenght - 1] - f[2][height - 1][lenght - 1] -
        f[3][height - 1][lenght - 1] - f[4][height - 1][lenght - 1] - f[5][height - 1][lenght - 1] -
        f[6][height - 1][lenght - 1] - f[7][height - 1][lenght - 1] - f[8][height - 1][lenght - 1];
}

void zou_he_top_left_corner_velocity(float f[][DIM][DIM], float u[][DIM][DIM], float rho[][DIM])
{

    f[1][0][0] = f[3][0][0] + (2.0 / 3.0) * rho[0][0] * u[0][0][0];
    f[4][0][0] = f[2][0][0] - (2.0 / 3.0) * rho[0][0] * u[1][0][0];
    f[8][0][0] = f[6][0][0] - (1.0 / 6.0) * rho[0][0] * u[0][0][0] + (1.0 / 6.0) * rho[0][0] * u[1][0][0];

    f[5][0][0] = 0.0;
    f[7][0][0] = 0.0;

    f[0][0][0] = rho[0][0] - f[1][0][0] - f[2][0][0] - f[3][0][0] - f[4][0][0] - f[5][0][0] - f[6][0][0] - f[7][0][0] -
                 f[8][0][0];
}

void zou_he_top_right_corner_velocity(float f[][DIM][DIM], float u[][DIM][DIM], float rho[][DIM], int lenght)
{

    f[3][0][lenght - 1] = f[1][0][lenght - 1] - (2.0 / 3.0) * rho[0][lenght - 1] * u[0][0][lenght - 1];
    f[4][0][lenght - 1] = f[2][0][lenght - 1] - (2.0 / 3.0) * rho[0][lenght - 1] * u[1][0][lenght - 1];
    f[7][0][lenght - 1] = f[5][0][lenght - 1] - (1.0 / 6.0) * rho[0][lenght - 1] * u[0][0][lenght - 1] - (1.0 / 6.0) * rho[0][lenght - 1] * u[1][0][lenght - 1];

    f[6][0][lenght - 1] = 0.0;
    f[8][0][lenght - 1] = 0.0;

    f[0][0][lenght - 1] = rho[0][lenght - 1] - f[1][0][lenght - 1] - f[2][0][lenght - 1] - f[3][0][lenght - 1] -
                          f[4][0][lenght - 1] - f[5][0][lenght - 1] - f[6][0][lenght - 1] - f[7][0][lenght - 1] -
                          f[8][0][lenght - 1];
}

void zou_he_left_wall_velocity(float f[][DIM][DIM], float u[][DIM][DIM], float rho[][DIM], int height)
{
    for (int i = 1; i < height - 1; i++)
    {

        rho[i][0] = (f[0][i][0] + f[2][i][0] + f[4][i][0] + 2.0 * (f[3][i][0] + f[6][i][0] + f[7][i][0])) / (1.0 - u[0][i][0]);

        f[1][i][0] = f[3][i][0] - 2.0 / 3.0 * rho[i][0] * u[0][i][0];

        f[5][i][0] = f[7][i][0] - 0.5 * (f[2][i][0] - f[4][i][0]) + 1.0 / 6.0 * rho[i][0] * u[0][i][0] + 0.5 * rho[i][0] * u[1][i][0];

        f[8][i][0] = f[6][i][0] + 0.5 * (f[2][i][0] - f[4][i][0]) + 1.0 / 6.0 * rho[i][0] * u[0][i][0] - 0.5 * rho[i][0] * u[1][i][0];
    }
}

void zou_he_right_wall_velocity(float f[][DIM][DIM], float u[][DIM][DIM], float rho[][DIM], int height, int lenght)
{
    for (int i = 1; i < height - 1; i++)
    {
        rho[i][lenght - 1] = (f[0][i][lenght - 1] + f[2][i][lenght - 1] + f[4][i][lenght - 1] + 2.0 * (f[1][i][lenght - 1] + f[5][i][lenght - 1] + f[8][i][lenght - 1])) / (1.0 + u[0][i][lenght - 1]);

        f[3][i][lenght - 1] = f[1][i][lenght - 1] - 2.0 / 3.0 * rho[i][lenght - 1] * u[0][i][lenght - 1];

        f[6][i][lenght - 1] = f[8][i][lenght - 1] - 0.5 * (f[2][i][lenght - 1] - f[4][i][lenght - 1]) - 1.0 / 6.0 * rho[i][lenght - 1] * u[0][i][lenght - 1] + 0.5 * rho[i][lenght - 1] * u[1][i][lenght - 1];

        f[7][i][lenght - 1] = f[5][i][lenght - 1] + 0.5 * (f[2][i][lenght - 1] - f[4][i][lenght - 1]) -  1.0 / 6.0 * rho[i][lenght - 1] * u[0][i][lenght - 1] - 0.5 * rho[i][lenght - 1] * u[1][i][lenght - 1];
    }
}

void zou_he_top_wall_velocity(float f[][DIM][DIM], float u[][DIM][DIM], float rho[][DIM], int lenght)
{
    for (int j = 1; j < lenght - 1; j++)
    {

        rho[0][j] = (f[0][0][j] + f[1][0][j] + f[3][0][j] + 2.0 * (f[2][0][j] + f[5][0][j] + f[6][0][j])) / (1.0 + u[1][0][j]);

        f[4][0][j] = f[2][0][j] - 2.0 / 3.0 * rho[0][j] * u[1][0][j];

        f[7][0][j] = f[5][0][j] + 0.5 * (f[1][0][j] - f[3][0][j]) - 1.0 / 6.0 * rho[0][j] * u[1][0][j] - 0.5 * rho[0][j] * u[0][0][j];

        f[8][0][j] = f[6][0][j] - 0.5 * (f[1][0][j] - f[3][0][j]) - 1.0 / 6.0 * rho[0][j] * u[1][0][j] + 0.5 * rho[0][j] * u[0][0][j];
    }
}

void zou_he_bottom_wall_velocity(float f[][DIM][DIM], float u[][DIM][DIM], float rho[][DIM], int height, int lenght)
{
    for (int j = 1; j < lenght - 1; j++)
    {

        rho[height - 1][j] = (f[0][height - 1][j] + f[1][height - 1][j] + f[3][height - 1][j] + 2.0 * (f[4][height - 1][j] + f[7][height - 1][j] + f[8][height - 1][j])) / (1.0 - u[1][height - 1][j]);

        f[2][height - 1][j] = f[4][height - 1][j] + 2.0 / 3.0 * rho[height - 1][j] * u[1][height - 1][j];

        f[5][height - 1][j] = f[7][height - 1][j] - 0.5 * (f[1][height - 1][j] - f[3][height - 1][j]) + 1.0 / 6.0 * rho[height - 1][j] * u[1][height - 1][j] + 0.5 * rho[height - 1][j] * u[0][height - 1][j];

        f[6][height - 1][j] = f[8][height - 1][j] + 0.5 * (f[1][height - 1][j] - f[3][height - 1][j]) + 1.0 / 6.0 * rho[height - 1][j] * u[1][height - 1][j] - 0.5 * rho[height - 1][j] * u[0][height - 1][j];
    }
}

int main()
{
    int length = DIM;
    int height = DIM;
    float t_max = 10.0;
    float re_lbm = 200.0;
    float u_lbm = 0.2;

    int it_max = (int)std::round(t_max * length / u_lbm);
    float sigma = 10.0 * length;
    int nx = length;
    int ny = height;
    float om_p = 1.0 / (0.5 + 3.0 * u_lbm * length / re_lbm);
    float om_m = 1.0 / (1.0 / (12.0 * u_lbm * length / re_lbm) + 0.5);

    int number_of_directions = 9;

    float u[2][DIM][DIM];
    float velocities[number_of_directions][2] = {
        // le nostre velocities
        {0., 0.}, {1., 0.}, {0., -1.}, {-1., 0.}, {0., 1.}, {1., -1.}, {-1., -1.}, {-1., 1.}, {1., 1.},
    };
    // # Weights
    float weigths[number_of_directions] = {
        4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36.,
    };
    // # Array for bounce-back
    int opposites[number_of_directions] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    for (int a = 0; a < 2; a++)
    {
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < length; j++)
            {
                u[a][i][j] = 0;
            }
        }
    }

    float f[number_of_directions][DIM][DIM];
    float fEq[number_of_directions][DIM][DIM];
    float fNew[number_of_directions][DIM][DIM];
    float rho[DIM][DIM];

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < length; j++)
        {
            for (int k = 0; k < number_of_directions; k++)
            {
                f[k][i][j] = 1.0;
            }
        }
    }

    // open output file
    std::ofstream fout("output.txt", std::ios::out);

    // write shape
    fout << nx << " " << ny << '\n';

    for (int it = 0; it <= it_max; ++it)
    {
        if (it % 100 == 0)
        {
            // save data
            fout << it << '\n';
            for (int j = 0; j < length; j++)
            {
                for (int i = 0; i < height; i++)
                {
                    fout << u[0][j][i] << ' ';
                }
            }
            fout << '\n';
            for (int j = 0; j < length; j++)
            {
                for (int i = 0; i < height; i++)
                {
                    fout << u[1][j][i] << ' ';
                }
            }
            fout << '\n';
            std::cout << "it = " << it << '\n';
        }

        computeMacroscopic(f, u, rho, height, length, number_of_directions, velocities);
        computeEquilibrium(u, rho, fEq, height, length, number_of_directions, velocities, weigths);
        collision_and_streaming(om_p, om_m, height, length, number_of_directions, fEq, f, fNew, opposites);
        const float uLidNow = u_lbm * (1.0 - std::exp(-static_cast<double>(it * it) / (2.0 * sigma * sigma)));
        setInlets(uLidNow, u, height, length);
        zou_he_bottom_wall_velocity(f, u, rho, height, length);
        zou_he_left_wall_velocity(f, u, rho, height);
        zou_he_right_wall_velocity(f, u, rho, height, length);
        zou_he_top_wall_velocity(f, u, rho, length);
        zou_he_bottom_left_corner_velocity(f, u, rho, height);
        zou_he_top_left_corner_velocity(f, u, rho);
        zou_he_top_right_corner_velocity(f, u, rho, length);
        zou_he_bottom_right_corner_velocity(f, u, rho, height, length);
        // TBD: Compute observables (drag, lift, etc)
    }

    // close output file
    fout.close();

    return 0;
}
