#include "GpuSimulation.cuh"

__device__ void step1dev(const int nx, const int ny, const int it, const int problem_type, const float u_lid,
                         const float om_p, const float halfOmpOmmSum, const float halfOmpOmmSub, const int row,
                         const int col, float *f, float *new_f, float &rho, float &ux, float &uy, int &boundary)
{
    const int velocitiesX[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    const int velocitiesY[9] = {0, 0, -1, 0, 1, -1, -1, 1, 1};
    const float weights[9] = {4.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0, 1.0 / 9.0,
                              1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
    const int opposite[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    // zouHe if it != 0
    if (it != 0)
    {
        // top wall
        if (boundary == 2)
        {
            rho = (f[0] + f[1] + f[3] + 2.0 * (f[2] + f[5] + f[6])) / (1.0 + uy);
            f[4] = f[2] - 2.0 / 3.0 * rho * uy;
            f[7] = f[5] + 0.5 * (f[1] - f[3]) - 0.5 * rho * ux - 1.0 / 6.0 * rho * uy;
            f[8] = f[6] - 0.5 * (f[1] - f[3]) + 0.5 * rho * ux - 1.0 / 6.0 * rho * uy;
        }
        // right wall
        else if (boundary == 1)
        {
            rho = (f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8])) / (1.0 + ux);
            f[3] = f[1] - 2.0 / 3.0 * rho * ux;
            f[6] = f[8] - 0.5 * (f[2] - f[4]) - 1.0 / 6.0 * rho * ux + 0.5 * rho * uy;
            f[7] = f[5] + 0.5 * (f[2] - f[4]) - 1.0 / 6.0 * rho * ux - 0.5 * rho * uy;
        }
        // bottom wall
        else if (boundary == 4)
        {
            rho = (f[0] + f[1] + f[3] + 2.0 * (f[4] + f[7] + f[8])) / (1.0 - uy);
            f[2] = f[4] + 2.0 / 3.0 * rho * uy;
            f[5] = f[7] - 0.5 * (f[1] - f[3]) + 0.5 * rho * ux + 1.0 / 6.0 * rho * uy;
            f[6] = f[8] + 0.5 * (f[1] - f[3]) - 0.5 * rho * ux + 1.0 / 6.0 * rho * uy;
        }
        // left wall
        else if (boundary == 3)
        {
            rho = (f[0] + f[2] + f[4] + 2.0 * (f[3] + f[6] + f[7])) / (1.0 - ux);
            f[1] = f[3] - 2.0 / 3.0 * rho * ux;
            f[5] = f[7] - 0.5 * (f[2] - f[4]) + 1.0 / 6.0 * rho * ux + 0.5 * rho * uy;
            f[8] = f[6] + 0.5 * (f[2] - f[4]) + 1.0 / 6.0 * rho * ux - 0.5 * rho * uy;
        }
        // top right corner
        else if (boundary == 5)
        {
            f[3] = f[1] - 2.0 / 3.0 * rho * ux;
            f[4] = f[2] - 2.0 / 3.0 * rho * uy;
            f[7] = f[5] - 1.0 / 6.0 * rho * ux - 1.0 / 6.0 * rho * uy;
            f[8] = 0;
            f[6] = 0;
            f[0] = rho - f[1] - f[2] - f[3] - f[4] - f[5] - f[7];
        }
        // bottom right corner
        else if (boundary == 8)
        {
            f[3] = f[1] - 2.0 / 3.0 * rho * ux;
            f[2] = f[4] + 2.0 / 3.0 * rho * uy;
            f[6] = f[8] + 1.0 / 6.0 * rho * uy - 1.0 / 6.0 * rho * ux;
            f[7] = 0;
            f[5] = 0;
            f[0] = rho - f[1] - f[2] - f[3] - f[4] - f[6] - f[8];
        }
        // bottom left corner
        else if (boundary == 7)
        {
            f[1] = f[3] + 2.0 / 3.0 * rho * ux;
            f[2] = f[4] + 2.0 / 3.0 * rho * uy;
            f[5] = f[7] + 1.0 / 6.0 * rho * ux + 1.0 / 6.0 * rho * uy;
            f[6] = 0;
            f[8] = 0;
            f[0] = rho - f[1] - f[2] - f[3] - f[4] - f[5] - f[7];
        }
        // top left corner
        else if (boundary == 6)
        {
            f[1] = f[3] + 2.0 / 3.0 * rho * ux;
            f[4] = f[2] - 2.0 / 3.0 * rho * uy;
            f[8] = f[6] - 1.0 / 6.0 * rho * ux + 1.0 / 6.0 * rho * uy;
            f[7] = 0;
            f[5] = 0;
            f[0] = rho - f[1] - f[2] - f[3] - f[4] - f[6] - f[8];
        }
    }

    // update macro
    rho = 0;
    ux = 0;
    uy = 0;
    for (int i = 0; i < 9; i++)
    {
        rho += f[i];
        ux += f[i] * velocitiesX[i];
        uy += f[i] * velocitiesY[i];
    }
    ux /= rho;
    uy /= rho;

    // set inlets
    if (problem_type == 1)
    {
        // if i'm a any boundary set ux uy to 0
        if (row == 0 || row == ny - 1 || col == 0 || col == nx - 1)
        {
            ux = 0;
            uy = 0;
        }
        // if i'm on the lid set ux to u_lid
        if (row == 0)
        {
            ux = u_lid;
        }
    }

    // equilibrium
    float feq[9];
    const float temp1 = 1.5 * (ux * ux + uy * uy);
    for (int i = 0; i < 9; i++)
    {
        const float temp2 = 3.0 * (velocitiesX[i] * ux + velocitiesY[i] * uy);
        feq[i] = weights[i] * rho * (1.0 + temp2 + 0.5 * temp2 * temp2 - temp1);
    }

    // collision for index 0
    new_f[0] = (1.0 - om_p) * f[0] + om_p * feq[0];

    // collision for other indices
    for (int i = 1; i < 9; i++)
    {
        new_f[i] = (1.0 - halfOmpOmmSum) * f[i] - halfOmpOmmSub * f[opposite[i]] + halfOmpOmmSum * feq[i] +
                   halfOmpOmmSub * feq[opposite[i]];
    }
}

__global__ void GpuSimulation::step1(const int nx, const int ny, const int it, const int problem_type,
                                     const float u_lid, const float om_p, const float halfOmpOmmSum,
                                     const float halfOmpOmmSub, float *f, float *new_f, float *rho, float *ux,
                                     float *uy, int *boundary, bool *obstacle)
{
    const int row = blockIdx.y * blockDim.y + threadIdx.y;
    const int col = blockIdx.x * blockDim.x + threadIdx.x;

    // return if out of bounds or obstacle
    if (row >= ny || col >= nx || obstacle[row * nx + col])
        return;

    const int index = row * nx + col;
    const int index9 = index * 9;

    step1dev(nx, ny, it, problem_type, u_lid, om_p, halfOmpOmmSum, halfOmpOmmSub, row, col, &f[index9], &new_f[index9],
             rho[index], ux[index], uy[index], boundary[index]);
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

    const int index = row * nx + col;
    const int index9 = index * 9;

    // stream for index 0
    f[index9] = new_f[index9];

    const int boundary_here = boundary[index];
    // stream for other indices
    for (int i = 1; i < 9; i++)
    {
        // check if there's a boundary in the way
        if ((velocitiesX[boundary_here] != velocitiesX[i] || velocitiesX[boundary_here] == 0) &&
            (velocitiesY[boundary_here] != velocitiesY[i] || velocitiesY[boundary_here] == 0))
        {
            // obtain new cell coordinates
            const int new_index9 = (row + velocitiesY[i]) * nx * 9 + (col + velocitiesX[i]) * 9;
            // stream
            f[new_index9 + i] = new_f[index9 + i];
        }
    }
}
