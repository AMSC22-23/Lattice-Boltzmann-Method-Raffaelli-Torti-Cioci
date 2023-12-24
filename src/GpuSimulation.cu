#include "GpuSimulation.cuh"

__device__ void zouHe(float *f, float *new_f, float *rho, float *ux, float *uy, int *boundary)
{
    return;
}

__device__ void updateMacro(float *f, float *new_f, float *rho, float *ux, float *uy, int *boundary,
                            const int *velocitiesX, const int *velocitiesY)
{
    *rho = 0;
    *ux = 0;
    *uy = 0;
    for (int i = 0; i < 9; i++)
    {
        *rho += f[i];
        *ux += f[i] * velocitiesX[i];
        *uy += f[i] * velocitiesY[i];
    }
}

__device__ void setInlets(float *ux, float *uy, const int problem_type, const float u_lid, const int row, const int col,
                          const int nx, const int ny)
{
    if (problem_type == 1)
    {
        // if i'm a any boundary set ux uy to 0
        if (row == 0 || row == ny - 1 || col == 0 || col == nx - 1)
        {
            *ux = 0;
            *uy = 0;
        }
        // if i'm on the lid set ux to u_lid
        if (row == 0)
        {
            *ux = u_lid;
        }
    }
}

__device__ void equilibriumCollision(float *f, float *new_f, float *ux, float *uy, const int *velocitiesX,
                                     const int *velocitiesY, const float *weights, const float om_p, const float om_m)
{
    return;
}

__global__ void GpuSimulation::step1(const int nx, const int ny, const int it, const int problem_type,
                                     const float u_lid, const float om_p, const float om_m, float *f, float *new_f,
                                     float *rho, float *ux, float *uy, int *boundary, bool *obstacle)
{
    const int row = blockIdx.y * blockDim.y + threadIdx.y;
    const int col = blockIdx.x * blockDim.x + threadIdx.x;

    // return if out of bounds or obstacle
    if (row >= ny || col >= nx || obstacle[row * nx + col])
        return;

    const int velocitiesX[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    const int velocitiesY[9] = {0, 0, -1, 0, 1, -1, -1, 1, 1};
    const float weights[9] = {4.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0, 1.0 / 9.0,
                              1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
    const int index = row * nx + col;
    const int index9 = index * 9;

    // zouhe if not it 0
    if (it != 0)
    {
        zouHe(&f[index9], &new_f[index9], &rho[index], &ux[index], &uy[index], &boundary[index]);
    }

    // update macro
    updateMacro(&f[index9], &new_f[index9], &rho[index], &ux[index], &uy[index], &boundary[index], velocitiesX,
                velocitiesY);

    // set inlets
    setInlets(&ux[index], &uy[index], problem_type, u_lid, row, col, nx, ny);

    // equilibrium and collision
}

__global__ void GpuSimulation::step2(const int nx, const int ny, float *f, float *new_f, int *boundary, bool *obstacle)
{
    const int row = blockIdx.y * blockDim.y + threadIdx.y;
    const int col = blockIdx.x * blockDim.x + threadIdx.x;

    // return if out of bounds or obstacle
    if (row >= ny || col >= nx || obstacle[row * nx + col])
        return;

    const int velocitiesX[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    const int velocitiesY[9] = {0, 0, -1, 0, 1, -1, -1, 1, 1};
}
