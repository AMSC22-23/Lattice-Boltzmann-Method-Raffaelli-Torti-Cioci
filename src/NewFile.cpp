#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unordered_map>

#define Nx = 10
#define Ny = 10

void save_velocity(int it, const float u[][100][100])
{
    char buffer[128];
    std::sprintf(buffer, "u.%.6intd.txt", it);
    std::ofstream fout(buffer, std::ios::out | std::ios::binary);
    fout.write(reinterpret_cast<const char *>(u[0][0]), sizeof(float) * 100 * 100);
}

static std::unordered_map<std::string, size_t> profiler;
#define profile(code, profiler_entry)                                                                                  \
    do                                                                                                                 \
    {                                                                                                                  \
        using namespace std::chrono;                                                                                   \
        const auto start_time = high_resolution_clock::now();                                                          \
        code;                                                                                                          \
        const auto end_time = high_resolution_clock::now();                                                            \
        profiler[profiler_entry] += duration_cast<nanoseconds>(end_time - start_time).count();                         \
    } while (0)

constexpr int number_of_directions = 9; //9 direzioni di velocità
// # D2Q9 Velocities
// constexpr double C[2][Q] = {{0., 1., -1., 0., 0., 1., -1., -1., 1.},
//                            {0., 0., 0., 1., -1., 1., -1, 1., -1.}};
constexpr float velocities[number_of_directions][2] = { // le nostre velocities
    {0., 0.}, {1., 0.}, {0., -1.}, {-1., 0.}, {0., 1.}, {1., -1.}, {-1., -1.}, {-1., 1.}, {1., 1.},
};
// # Weights
constexpr float weigths[number_of_directions] = {
    4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36.,
};
// # Array for bounce-back
constexpr int opposites[number_of_directions] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

void setVelocitiesAtWalls (double sigma, double u_lid, int it, float uTop[][100], float uBot[][100], float uLeft[][100], float uRight[][100] )
{
    const float ret = 1.0 - std::exp(-static_cast<double>(it*it) / ( 2.0 * sigma * sigma));
    for (int i = 0; i < 100; i++)
    {
        uTop[0][i] = u_lid * ret;
        uBot[0][i] = 0.0;
        uLeft[1][i] = 0.0;
        uRight[1][i] = 0.0;
    }
}
// u [coordinata x o y][posizione della cella delle x][posizione della cella delle y]
void computeEquilibrium(const float u[][100][100], const float rho[][100],  float fEq[][100][100])
{
    for (int i = 0; i < 9; i++)
    {
        for(int j=0; j<100; j++)
        {
            for(int k=0; k<100; k++)
            {
                const float temp1 = 1.5 * ( u[0][j][k] * u[0][j][k] + u[1][j][k] * u[1][j][k]);
                const float temp2 = 3 * ( velocities[i][0] * u[0][j][k] + velocities[i][1] * u[1][j][k]);
                fEq[i][j][k] = weigths[i] * rho[j][k] * (1.0 + temp2 + 0.5 * temp2 * temp2 - temp1);
            }
        }
    }
}

void computeMacroscopic ( float f[][100][100], float u[][100][100], float rho[][100])
{
    for(int i = 0; i<100; i++)
    {
        for (int j = 0; j<100; j++)
        {
            rho[i][j] = 0;
            u[0][i][j] = 0;
            u[1][i][j] = 0;
            for(int k = 0; k<9; k++)
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
void collision_and_streaming(double om_p, double om_m, int lx, int ly, const float fEq[][100][100], float f[][100][100], float fNew[][100][100])
{
    for (int i = 0; i < lx + 1; ++i)
    {
        for(int j = 0; j < ly + 1; ++j)
        {
            fNew[0][i][j] = (1.0 - om_p) * f[0][i][j] + om_p * fEq[0][i][j];
            f[0][i][j] = fNew[0][i][j];
        }
    }

    //collide for other indices
    for(int i = 0; i < 9; ++i)
    {
        const auto qb = opposites[i];
        for(int j = 0; j < lx + 1; ++j)
        {
            for(int k = 0; k < ly + 1; ++k)
            {
                fNew[i][j][k] = (1.0 - 0.5 * (om_p+ om_m)) * f[i][j][k] - 0.5 * (om_p - om_m) * f[qb][j][k] + 0.5 * (om_m + om_p) * fEq[i][j][k] + 0.5 * (om_m - om_p) * fEq[qb][j][k];
            }
        }
    }

    //stream
    for(int i = 0; i < lx; i++)
    {
        for(int j = 0; j < ly; j++)
        {
            f[1][i+1][j] = fNew[1][i][j];
            f[2][i][j] = fNew[2][i][j+1];
            f[3][i][j] = fNew[3][i+1][j];
            f[4][i][j+1] = fNew[4][i][j];
            f[5][i+1][j] = fNew[5][i][j+1];
            f[6][i][j] = fNew[6][i+1][j+1];
            f[7][i][j+1] = fNew[7][i+1][j];
            f[8][i+1][j+1] = fNew[8][i][j];
        }
    }
    for(int i = 0; i < lx; i++)
    {
        f[1][i+1][ly] = fNew[1][i][ly];
        f[3][i][ly] = fNew[3][i+1][ly];
    }
    for(int j = 0; j < ly; j++)
    {
        f[2][lx][j] = fNew[2][lx][j+1];
        f[4][lx][j+1] = fNew[4][lx][j];
    }
}

void zou_he_bottom_left_corner_velocity(float f[][100][100], float u[][100][100], float rho[][100], int ly)
{
    u[0][0][ly] = u[0][1][ly]; //forse questa va messa come la uWall
    u[1][0][ly] = u[1][1][ly]; //idem

    rho[0][ly] = rho[1][ly];

    f[1][0][ly] = f[3][0][ly] + (2.0 / 3.0) * rho[0][ly] * u[0][0][ly];

    f[2][0][ly] = f[4][0][ly] + (2.0 / 3.0) * rho[0][ly] * u[1][0][ly];

    f[5][0][ly] = f[7][0][ly] + (1.0 / 6.0) * rho[0][ly] * u[0][0][ly] + (1.0 / 6.0) * rho[0][ly] * u[1][0][ly];

    f[6][0][ly] = 0.0;
    f[8][0][ly] = 0.0;

    f[0][0][ly] = rho[0][ly] - f[1][0][ly] - f[2][0][ly] - f[3][0][ly] - f[4][0][ly] - f[5][0][ly] - f[6][0][ly] - f[7][0][ly] -
                 f[8][0][ly];
}

void zou_he_bottom_right_corner_velocity(float f[][100][100], float u[][100][100], float rho[][100], int lx)
{
    u[0][lx][0] = u[0][lx-1][0]; //forse questa va messa come la uWall
    u[1][lx][0] = u[1][lx-1][0]; //idem

    rho[lx][0] = rho[lx-1][0];

    f[3][lx][0] = f[1][lx][0] - (2.0 / 3.0) * rho[9][0] * u[0][lx][0];

    f[2][lx][0] = f[4][lx][0] + (2.0 / 3.0) * rho[9][0] * u[1][lx][0];

    f[6][lx][0] = f[8][lx][0] - (1.0 / 6.0) * rho[9][0] * u[0][lx][0] + (1.0 / 6.0) * rho[9][0] * u[1][lx][0];

    f[5][lx][0] = 0.0;
    f[7][lx][0] = 0.0;

    f[0][lx][0] = rho[9][0] - f[1][lx][0] - f[2][lx][0] - f[3][lx][0] - f[4][lx][0] - f[5][lx][0] - f[6][lx][0] - f[7][lx][0] -
                 f[8][lx][0];
}

void zou_he_top_left_corner_velocity( float f[][100][100], float u[][100][100], float rho[][100])
{
    u[0][0][0] = u[0][1][0]; //forse questa va messa come la uWall
    u[1][0][0] = u[1][1][0]; //idem

    rho[0][0] = rho[1][0];

    f[1][0][0] = f[3][0][0] + (2.0 / 3.0) * rho[0][0] * u[0][0][0];

    f[4][0][0] = f[2][0][0] - (2.0 / 3.0) * rho[0][0] * u[1][0][0];

    f[8][0][0] = f[6][0][0] - (1.0 / 6.0) * rho[0][0] * u[0][0][0] + (1.0 / 6.0) * rho[0][0] * u[1][0][0];

    f[5][0][0] = 0.0;
    f[7][0][0] = 0.0;

    f[0][0][0] = rho[0][0] - f[1][0][0] - f[2][0][0] - f[3][0][0] - f[4][0][0] - f[5][0][0] - f[6][0][0] - f[7][0][0] -
                 f[8][0][0];
}

void zou_he_top_right_corner_velocity( float f[][100][100], float u[][100][100], float rho[][100], int lx)
{
    u[0][lx][0] = u[0][lx-1][0]; //forse questa va messa come la uWall
    u[1][lx][0] = u[1][lx-1][0]; //idem

    rho[lx][0] = rho[lx-1][0];

    f[3][lx][0] = f[1][lx][0] - (2.0 / 3.0) * rho[lx][0] * u[0][lx][0];

    f[4][lx][0] = f[2][lx][0] - (2.0 / 3.0) * rho[lx][0] * u[1][lx][0];

    f[7][lx][0] = f[5][lx][0] - (1.0 / 6.0) * rho[lx][0] * u[0][lx][0] - (1.0 / 6.0) * rho[lx][0] * u[1][lx][0];

    f[6][lx][0] = 0.0;
    f[8][lx][0] = 0.0;

    f[0][lx][0] = rho[lx][0] - f[1][lx][0] - f[2][lx][0] - f[3][lx][0] - f[4][lx][0] - f[5][lx][0] - f[6][lx][0] - f[7][lx][0] -
                 f[8][lx][0];
}

void zou_he_left_wall_velocity( float f[][100][100], float u[][100][100], float rho[][100], int ly, float uLeft[][100])
{
    for(int i = 0; i < ly ; i++)
    {
        u[0][0][i] = uLeft[0][i];
        u[1][0][i] = uLeft[1][i];

        rho[0][i] = (f[0][0][i] + f[2][0][i] + f[4][0][i] + 2.0 * (f[3][0][i] + f[6][0][i] + f[7][0][i])) / (1.0 - u[0][0][i]);

        f[1][0][i] = f[3][0][i] + 2.0/3.0 * rho[0][i] * u[0][0][i];
        f[5][0][i] = f[7][0][i] - 0.5 * (f[2][0][i] - f[4][0][i]) + 1.0/6.0 * rho[0][i] * u[0][0][i] + 0.5 * rho[0][i] * u[1][0][i];
        f[8][0][i] = f[6][0][i] + 0.5 * (f[2][0][i] - f[4][0][i]) + 1.0/6.0 * rho[0][i] * u[0][0][i] - 0.5 * rho[0][i] * u[1][0][i];
    }
}
void zou_he_right_wall_velocity(float f[][100][100], float u[][100][100], float rho[][100], int ly, float uRight[][100], int lx)
{
    for (int i = 0; i<ly; i++)
    {
        u[0][lx][i] = uRight[0][i];
        u[1][lx][i] = uRight[1][i];

        rho[lx][i] = (f[0][lx][i] + f[2][lx][i] + f[4][lx][i] + 2.0 * (f[1][lx][i] + f[5][lx][i] + f[8][lx][i])) / (1.0 + u[0][lx][i]);

        f[3][lx][i] = f[1][lx][i] - 2.0 / 3.0 * rho[lx][i] * u[0][lx][i];
        f[6][lx][i] = f[8][lx][i] - 0.5 * (f[2][lx][i] - f[4][lx][i]) - 1.0/6.0 * rho[lx][i] * u[0][lx][i] - 0.5 * rho[lx][i] * u[1][lx][i];
        f[7][lx][i] = f[5][lx][i] + 0.5 * (f[2][lx][i] - f[4][lx][i]) - 1.0/6.0 * rho[lx][i] * u[0][lx][i] - 0.5 * rho[lx][i] * u[1][lx][i];
    }
}
void zou_he_top_wall_velocity(float f[][100][100], float u[][100][100], float rho[][100], int ly, float uTop[][100], int lx)
{
    for(int i = 0; i < lx; i++)
    {
        u[0][i][0] = uTop[0][i];
        u[1][i][0] = uTop[1][i];

        rho[i][0] = (f[0][i][0] + f[1][i][0] + f[3][i][0] + 2.0 * (f[2][i][0] + f[5][i][0] + f[6][i][0])) / (1.0 + u[1][i][0]);

        f[4][i][0] = f[2][i][0] - 2.0/3.0 * rho[i][0] * u[1][i][0];
        f[8][i][0] = f[6][i][0] - 0.5 * (f[1][i][0] - f[3][i][0]) - 1.0/6.0 * rho[i][0] * u[1][i][0] + 0.5 * rho[i][0] * u[0][i][0];
        f[7][i][0] = f[5][i][0] + 0.5 * (f[1][i][0] - f[3][i][0]) - 1.0/6.0 * rho[i][0] * u[1][i][0] - 0.5 * rho[i][0] * u[0][i][0];
    }
}

void zou_he_bottom_wall_velocity(float f[][100][100], float u[][100][100], float rho[][100], int ly, float uBot[][100], int lx)
{
    for(int i = 0; i < lx; i++)
    {
        u[0][i][ly] = uBot[0][i];
        u[1][i][ly] = uBot[1][i];

        rho[i][ly] = (f[0][i][ly] + f[1][i][ly] + f[3][i][ly] + 2.0 * (f[4][i][ly] + f[7][i][ly] + f[8][i][ly])) / (1.0 - u[1][i][ly]);

        f[2][i][ly] = f[4][i][ly] + 2.0 / 3.0 * rho[i][ly] * u[1][i][ly];
        f[5][i][ly] = f[7][i][ly] - 0.5 * (f[1][i][ly] - f[3][i][ly]) + 1.0/6.0 * rho[i][ly] * u[1][i][ly] + 0.5 * rho[i][ly] * u[0][i][ly];
        f[6][i][ly] = f[8][i][ly] + 0.5 * (f[1][i][ly] - f[3][i][ly]) + 1.0/6.0 * rho[i][ly] * u[1][i][ly] - 0.5 * rho[i][ly] * u[0][i][ly];
    }
}


int main(){

    using namespace std::chrono;
    const auto t0 = high_resolution_clock::now();

    const double re_lbm = 100.0; // reynolds number

    const double u_lbm = 0.2; // velocità della lid
    const double rho_lbm = 1.0; //densità 
    const double t_max = 20.0; 
    const double x_min = 0.0; 
    const double x_max = 1.0;
    const double y_min = 0.0;
    const double y_max = 1.0;
    const double c_s = 1.0 / std::sqrt(3.0); 

    const int npts = 100; 
    const int ny = npts; // 100: lattice height
    const double nu_lbm = u_lbm * npts / re_lbm; // 0.2
    const double tau_lbm = 0.5 + nu_lbm / (c_s * c_s); // 1.1
    const double dt = re_lbm * nu_lbm / (npts * npts); // 0.002
    // const double dy = (y_max - y_min) / ny;
    // const double dx = dy;
    const int nx = (int)std::round(ny * (x_max - x_min) / (y_max - y_min)); // 100: lattice width

    const int it_max = (int)std::round(t_max / dt); // 10'000
    const double sigma = 10.0 * nx; // 1000

    const int lx = nx - 1; // 99
    const int ly = ny - 1; // 99

    const double tau_p = tau_lbm;
    const double lambda_trt = 1.0 / 4.0; // Best for stability
    const double tau_m = lambda_trt / (tau_p - 0.5) + 0.5;
    const double om_p = 1.0 / tau_p;
    const double om_m = 1.0 / tau_m;
    float u[2][100][100];
    for(int i = 0; i< 2; i++)
    {
        for(int j = 0; j<100; j++)
        {
            for(int k = 0; k<10; k++)
            {
                u[i][j][k] = 0;
            }
        }
    }
    float fNew[9][100][100];
    float f[9][100][100];
    float fEq[9][100][100];
    float rho [100][100];
    for(int i = 0; i< 9; i++)
    {
        for(int j = 0; j<100; j++)
        {
            for(int k = 0; k<100; k++)
            {
                fEq[i][j][k] = 0;
                f[i][j][k] = rho_lbm;
                rho[j][k] = 1.0;
            }
        }
    }
    float uTop [2][100];
    float uBot [2][100];
    float uLeft [2][100];
    float uRight [2][100];
    for (int i = 0; i < 2; i++)
    {
        for(int j=0; j<100; j++)
        {
            uTop[i][j] = 0;
            uBot[i][j] = 0;
            uLeft[i][j] = 0;
            uRight[i][j] = 0;
        }
    }
    
    // Initialize and compute first equilibrium
    setVelocitiesAtWalls(sigma, u_lbm, 0, uTop, uBot, uLeft, uRight); //0 perchè è la prima iterazione (in realtà inizializzazione)
    computeEquilibrium(u, rho, fEq);
    for(int i = 0; i< 9; i++)
    {
        for(int j = 0; j<100; j++)
        {
            for(int k = 0; k<100; k++)
            {
                f[i][j][k] = fEq[i][j][k];
            }
        }
    }

    for (int it = 0; it < it_max + 1; ++it)
    {
        profile(
            if (it % 100 == 0) {
                std::cout << "Iteration: " << it << " / " << it_max << std::endl;
                save_velocity(it, u);
            },
            "IO");
        // 1. Set inlets
        profile(setVelocitiesAtWalls(sigma, u_lbm, it, uTop, uBot, uLeft, uRight);, "set_inlets");
        // 2. Compute macroscopic fields
        profile(computeMacroscopic(f, u, rho);, "compute_macroscopic");
        // 4. Compute equilibrium state
        profile(computeEquilibrium(u, rho, fEq);, "compute_equilibrium");
        // 5. Streaming
        profile(collision_and_streaming(om_p, om_m, lx, ly, fEq, f, fNew);, "collision_and_streaming");
        // 6. Boundary conditions
        profile(zou_he_bottom_wall_velocity(f, u, rho, ly, uBot, lx); zou_he_left_wall_velocity(f, u, rho, ly, uLeft);
                zou_he_right_wall_velocity(f, u, rho, ly, uRight, lx); zou_he_top_wall_velocity(f, u, rho, ly, uTop, lx);
                zou_he_bottom_left_corner_velocity(f, u, rho, ly); zou_he_top_left_corner_velocity(f, u, rho);
                zou_he_top_right_corner_velocity(f, u, rho, lx); zou_he_bottom_right_corner_velocity(f, u, rho, lx);
                , "boundary_conditions");
        // TBD: Compute observables (drag, lift, etc)
    }

  

    return 0;
}


