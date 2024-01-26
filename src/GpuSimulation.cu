#include "GpuSimulation.cuh"
#include "Lattice.hpp"
#include <iostream>

__device__ void step1dev(const int nx, const int ny, const int it, const int problem_type, const float u_lid,
                         const float om_p, const float halfOmpOmmSum, const float halfOmpOmmSub, const int row,
                         const int col, float *f, float *new_f, float *rho, float &ux, float &uy, const int *boundary)
{
    const int velocitiesX[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    const int velocitiesY[9] = {0, 0, -1, 0, 1, -1, -1, 1, 1};
    const float weights[9] = {4.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0, 1.0 / 9.0,
                              1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
    const int opposite[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    const int index = row * nx + col;
    float &rho_here = rho[index];

    // if i'm a any boundary set u to 0
    if (row == 0 || row == ny - 1 || col == 0 || col == nx - 1)
    {
        ux = 0;
        uy = 0;
    }

    // set lid inlet
    if (problem_type == 1 && row == 0)
    {
        ux = u_lid;
    }
    // set parabolic profile inlet
    else if (problem_type == 2 && col == 0)
    {
        const float halfDim = static_cast<float>(ny - 1) / 2.0;
        const float temp = static_cast<float>(row / halfDim) - 1.0;
        const float mul = 1.0 - temp * temp;
        ux = u_lid * mul;
    }

    // zouHe

    // top wall
    if (row == 0 && col != 0 && col != nx - 1)
    {
        rho_here = (f[0] + f[1] + f[3] + 2.0 * (f[2] + f[5] + f[6])) / (1.0 + uy);
        f[4] = f[2] - 2.0 / 3.0 * rho_here * uy;
        f[7] = f[5] + 0.5 * (f[1] - f[3]) - 0.5 * rho_here * ux - 1.0 / 6.0 * rho_here * uy;
        f[8] = f[6] - 0.5 * (f[1] - f[3]) + 0.5 * rho_here * ux - 1.0 / 6.0 * rho_here * uy;
    }
    // right wall
    else if (col == nx - 1 && row != 0 && row != ny - 1)
    {
        if (problem_type == 1)
        {
            rho_here = (f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8])) / (1.0 + ux);
            f[3] = f[1] - 2.0 / 3.0 * rho_here * ux;
            f[6] = f[8] - 0.5 * (f[2] - f[4]) - 1.0 / 6.0 * rho_here * ux + 0.5 * rho_here * uy;
            f[7] = f[5] + 0.5 * (f[2] - f[4]) - 1.0 / 6.0 * rho_here * ux - 0.5 * rho_here * uy;
        }
        else if (problem_type == 2)
        {
            rho_here = 1;
            ux = f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8]) - 1.0;
            f[3] = f[1] - 2.0 / 3.0 * ux;
            f[6] = f[8] - 0.5 * (f[2] - f[4]) - 1.0 / 6.0 * ux;
            f[7] = f[5] + 0.5 * (f[2] - f[4]) - 1.0 / 6.0 * ux;
        }
    }
    // bottom wall
    else if (row == ny - 1 && col != 0 && col != nx - 1)
    {
        rho_here = (f[0] + f[1] + f[3] + 2.0 * (f[4] + f[7] + f[8])) / (1.0 - uy);
        f[2] = f[4] + 2.0 / 3.0 * rho_here * uy;
        f[5] = f[7] - 0.5 * (f[1] - f[3]) + 0.5 * rho_here * ux + 1.0 / 6.0 * rho_here * uy;
        f[6] = f[8] + 0.5 * (f[1] - f[3]) - 0.5 * rho_here * ux + 1.0 / 6.0 * rho_here * uy;
    }
    // left wall
    else if (col == 0 && row != 0 && row != ny - 1)
    {
        rho_here = (f[0] + f[2] + f[4] + 2.0 * (f[3] + f[7] + f[6])) / (1.0 - ux);
        f[1] = f[3] + 2.0 / 3.0 * rho_here * ux;
        f[5] = f[7] - 0.5 * (f[2] - f[4]) + 1.0 / 6.0 * rho_here * ux + 0.5 * rho_here * uy;
        f[8] = f[6] + 0.5 * (f[2] - f[4]) + 1.0 / 6.0 * rho_here * ux - 0.5 * rho_here * uy;
    }
    // top right corner
    else if (row == 0 && col == nx - 1)
    {
        rho_here = rho[index - 1];
        f[3] = f[1] - 2.0 / 3.0 * rho_here * ux;
        f[4] = f[2] - 2.0 / 3.0 * rho_here * uy;
        f[7] = f[5] - 1.0 / 6.0 * rho_here * ux - 1.0 / 6.0 * rho_here * uy;
        f[8] = 0;
        f[6] = 0;
        f[0] = rho_here - f[1] - f[2] - f[3] - f[4] - f[5] - f[7];
    }
    // bottom right corner
    else if (row == ny - 1 && col == nx - 1)
    {
        rho_here = rho[index - 1];
        f[3] = f[1] - 2.0 / 3.0 * rho_here * ux;
        f[2] = f[4] + 2.0 / 3.0 * rho_here * uy;
        f[6] = f[8] + 1.0 / 6.0 * rho_here * uy - 1.0 / 6.0 * rho_here * ux;
        f[7] = 0;
        f[5] = 0;
        f[0] = rho_here - f[1] - f[2] - f[3] - f[4] - f[6] - f[8];
    }
    // bottom left corner
    else if (row == ny - 1 && col == 0)
    {
        rho_here = rho[index + 1];
        f[1] = f[3] + 2.0 / 3.0 * rho_here * ux;
        f[2] = f[4] + 2.0 / 3.0 * rho_here * uy;
        f[5] = f[7] + 1.0 / 6.0 * rho_here * ux + 1.0 / 6.0 * rho_here * uy;
        f[6] = 0;
        f[8] = 0;
        f[0] = rho_here - f[1] - f[2] - f[3] - f[4] - f[5] - f[7];
    }
    // top left corner
    else if (row == 0 && col == 0)
    {
        rho_here = rho[index + 1];
        f[1] = f[3] + 2.0 / 3.0 * rho_here * ux;
        f[4] = f[2] - 2.0 / 3.0 * rho_here * uy;
        f[8] = f[6] + 1.0 / 6.0 * rho_here * ux - 1.0 / 6.0 * rho_here * uy;
        f[7] = 0;
        f[5] = 0;
        f[0] = rho_here - f[1] - f[2] - f[3] - f[4] - f[6] - f[8];
    }

    // update macro
    rho_here = 0;
    ux = 0;
    uy = 0;
    for (int i = 0; i < 9; i++)
    {
        rho_here += f[i];
        ux += f[i] * velocitiesX[i];
        uy += f[i] * velocitiesY[i];
    }
    ux /= rho_here;
    uy /= rho_here;

    // equilibrium
    float feq[9];
    const float temp1 = 1.5 * (ux * ux + uy * uy);
    for (int i = 0; i < 9; i++)
    {
        const float temp2 = 3.0 * (velocitiesX[i] * ux + velocitiesY[i] * uy);
        feq[i] = weights[i] * rho_here * (1.0 + temp2 + 0.5 * temp2 * temp2 - temp1);
    }

    // collision for index 0
    new_f[0] = (1.0 - om_p) * f[0] + om_p * feq[0];

    // collision for other indices
    for (int i = 1; i < 9; i++)
    {
        new_f[i] = (1.0 - halfOmpOmmSum) * f[i] - halfOmpOmmSub * f[opposite[i]] + halfOmpOmmSum * feq[i] +
                   halfOmpOmmSub * feq[opposite[i]];
    }

    if (problem_type == 2)
    {
        // regular bounce back
        if (boundary[0] == 1)
        {
            f[3] = new_f[1];
        }
        else if (boundary[0] == -1)
        {
            f[1] = new_f[3];
        }
        if (boundary[1] == 1)
        {
            f[2] = new_f[4];
        }
        else if (boundary[1] == -1)
        {
            f[4] = new_f[2];
        }
        if (boundary[2] == 1)
        {
            f[6] = new_f[8];
        }
        else if (boundary[2] == -1)
        {
            f[8] = new_f[6];
        }
        if (boundary[3] == 1)
        {
            f[5] = new_f[7];
        }
        else if (boundary[3] == -1)
        {
            f[7] = new_f[5];
        }
    }
}

__global__ void step1(const int nx, const int ny, const int it, const int problem_type, const float u_lid,
                      const float om_p, const float halfOmpOmmSum, const float halfOmpOmmSub, float *f, float *new_f,
                      float *rho, float *ux, float *uy, const int *boundary, const bool *obstacle)
{
    const int row = blockIdx.y * blockDim.y + threadIdx.y;
    const int col = blockIdx.x * blockDim.x + threadIdx.x;

    // return if out of bounds or obstacle
    if (row >= ny || col >= nx || obstacle[row * nx + col])
        return;

    const int index = row * nx + col;
    const int index9 = index * 9;
    const int index4 = index * 4;

    step1dev(nx, ny, it, problem_type, u_lid, om_p, halfOmpOmmSum, halfOmpOmmSub, row, col, &f[index9], &new_f[index9],
             rho, ux[index], uy[index], &boundary[index4]);
}

__global__ void step2(const int nx, const int ny, float *f, const float *new_f, const int *boundary,
                      const bool *obstacle)
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

    // stream for other indices
    for (int i = 1; i < 9; i++)
    {
        // obtain new indices
        const int new_row = row + velocitiesY[i];
        const int new_col = col + velocitiesX[i];
        const int new_index = new_row * nx + new_col;
        // stream if new index is not out of bounds or obstacle
        if (new_row >= 0 && new_row < ny && new_col >= 0 && new_col < nx && !obstacle[new_index])
        {
            const int new_index9 = new_index * 9;
            f[new_index9 + i] = new_f[index9 + i];
        }
    }
}

__global__ void calculateLiftAndDragKernel(float *lift, float *drag, const float *f, const float *new_f,
                                           const int *boundary, const bool *obstacle, const int nx, const int ny)
{
    const int row = blockIdx.y * blockDim.y + threadIdx.y;
    const int col = blockIdx.x * blockDim.x + threadIdx.x;

    // return if out of bounds or obstacle
    if (row >= ny || col >= nx || obstacle[row * nx + col])
        return;

    const int index = row * nx + col;
    const int index4 = index * 4;
    const int index9 = index * 9;
    const int *boundary_here = &boundary[index4];
    const float *f_here = &f[index9];
    const float *new_f_here = &new_f[index9];

    // Calculate lift and drag for this thread's data
    const float adj = -2.0 / 0.04;
    float localLift = 0;
    float localDrag = 0;

    if (boundary_here[0] == 1)
    {
        localDrag += new_f_here[1] + f_here[3];
    }
    else if (boundary_here[0] == -1)
    {
        localDrag -= new_f_here[3] + f_here[1];
    }
    if (boundary_here[1] == 1)
    {
        localLift += new_f_here[4] + f_here[2];
    }
    else if (boundary_here[1] == -1)
    {
        localLift -= new_f_here[2] + f_here[4];
    }
    if (boundary_here[2] == 1)
    {
        localDrag += new_f_here[8] + f_here[6];
        localLift += new_f_here[8] + f_here[6];
    }
    else if (boundary_here[2] == -1)
    {
        localDrag -= new_f_here[6] + f_here[8];
        localLift -= new_f_here[6] + f_here[8];
    }
    if (boundary_here[3] == 1)
    {
        localDrag -= new_f_here[7] + f_here[5];
        localLift += new_f_here[7] + f_here[5];
    }
    else if (boundary_here[3] == -1)
    {
        localDrag += new_f_here[5] + f_here[7];
        localLift -= new_f_here[5] + f_here[7];
    }

    // Atomic add to global lift and drag
    atomicAdd(lift, localLift * adj);
    atomicAdd(drag, localDrag * adj);
}

void GpuSimulation::cudaCaller(const NDimensionalMatrix<Cell> &cells, const float sigma, const float omP,
                               const float omM, const int maxIt, const float uLid, const int problemType,
                               const int plotSteps, std::ofstream &velocity_out, std::ofstream &lift_drag_out)
{
    const int nx = cells.getShape().at(0);
    const int ny = cells.getShape().at(1);
    const int totalSize = nx * ny;
    int timeInstant = 0;

    // device allocations

    float *host_f, *host_new_f, *host_rho, *host_ux, *host_uy, *dev_f, *dev_new_f, *dev_rho, *dev_ux, *dev_uy;
    int *host_boundary, *dev_boundary;
    bool *host_obstacle, *dev_obstacle;

    // allocate memory on device
    cudaMalloc((void **)&dev_f, totalSize * 9 * sizeof(float));
    cudaMalloc((void **)&dev_new_f, totalSize * 9 * sizeof(float));
    cudaMalloc((void **)&dev_rho, totalSize * sizeof(float));
    cudaMalloc((void **)&dev_ux, totalSize * sizeof(float));
    cudaMalloc((void **)&dev_uy, totalSize * sizeof(float));
    cudaMalloc((void **)&dev_boundary, totalSize * 4 * sizeof(int));
    cudaMalloc((void **)&dev_obstacle, totalSize * sizeof(bool));

    // allocate memory on host
    cudaMallocHost((void **)&host_f, totalSize * 9 * sizeof(float));
    cudaMallocHost((void **)&host_new_f, totalSize * 9 * sizeof(float));
    cudaMallocHost((void **)&host_rho, totalSize * sizeof(float));
    cudaMallocHost((void **)&host_ux, totalSize * sizeof(float));
    cudaMallocHost((void **)&host_uy, totalSize * sizeof(float));
    cudaMallocHost((void **)&host_boundary, totalSize * 4 * sizeof(int));
    cudaMallocHost((void **)&host_obstacle, totalSize * sizeof(bool));

    // set host data
    for (int i = 0; i < totalSize; ++i)
    {
        const Cell &cell = cells.getElementAtFlatIndex(i);
        const std::vector<float> &f = cell.getF();
        const std::vector<float> &new_f = cell.getNewF();
        const std::vector<float> &macroU = cell.getMacroU();
        const std::vector<int> &boundary = cell.getBoundary();
        const bool &obstacle = cell.isObstacle();

        for (int j = 0; j < 9; ++j)
        {
            host_f[i * 9 + j] = f.at(j);
            host_new_f[i * 9 + j] = new_f.at(j);
        }
        host_rho[i] = cell.getRho();
        host_ux[i] = macroU.at(0);
        host_uy[i] = macroU.at(1);
        for (int j = 0; j < 4; ++j)
        {
            host_boundary[i * 4 + j] = boundary.at(j);
        }
        host_obstacle[i] = obstacle;
    }

    // copy host data to device
    cudaMemcpy(dev_f, host_f, totalSize * 9 * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_new_f, host_new_f, totalSize * 9 * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_rho, host_rho, totalSize * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_ux, host_ux, totalSize * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_uy, host_uy, totalSize * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_boundary, host_boundary, totalSize * 4 * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_obstacle, host_obstacle, totalSize * sizeof(bool), cudaMemcpyHostToDevice);

    // free host memory
    cudaFreeHost(host_f);
    cudaFreeHost(host_new_f);
    cudaFreeHost(host_rho);
    cudaFreeHost(host_boundary);
    cudaFreeHost(host_obstacle);

    // variables for lift and drag
    float drag, lift;
    float *dev_drag, *dev_lift;
    cudaMalloc((void **)&dev_drag, sizeof(float));
    cudaMalloc((void **)&dev_lift, sizeof(float));

    const dim3 threadsPerBlock(24, 24);
    const dim3 numBlocks(ceil(nx / 24.0), ceil(ny / 24.0));
    // loop
    const float temp = 2.0 * sigma * sigma;
    const float halfOmpOmmSub = 0.5 * (omP - omM);
    const float halfOmpOmmSum = 0.5 * (omP + omM);
    while (timeInstant <= maxIt)
    {
        const float uLidNow = uLid * (1.0 - std::exp(-static_cast<double>(timeInstant * timeInstant) / temp));
        step1<<<numBlocks, threadsPerBlock>>>(nx, ny, timeInstant, problemType, uLidNow, omP, halfOmpOmmSum,
                                              halfOmpOmmSub, dev_f, dev_new_f, dev_rho, dev_ux, dev_uy, dev_boundary,
                                              dev_obstacle);
        step2<<<numBlocks, threadsPerBlock>>>(nx, ny, dev_f, dev_new_f, dev_boundary, dev_obstacle);

        if (problemType == 2 && timeInstant % (maxIt / plotSteps) == 0)
        {
            // clear lift and drag
            cudaMemset(dev_lift, 0, sizeof(float));
            cudaMemset(dev_drag, 0, sizeof(float));
            // Launch kernel
            calculateLiftAndDragKernel<<<numBlocks, threadsPerBlock>>>(dev_lift, dev_drag, dev_f, dev_new_f,
                                                                       dev_boundary, dev_obstacle, nx, ny);

            // Copy results back to host
            cudaMemcpy(&lift, dev_lift, sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(&drag, dev_drag, sizeof(float), cudaMemcpyDeviceToHost);
        }

        // write to file every maxIt/plotSteps time steps
        if (timeInstant % (maxIt / plotSteps) == 0)
        {
            // copy ux and uy to host
            cudaMemcpy(host_ux, dev_ux, totalSize * sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(host_uy, dev_uy, totalSize * sizeof(float), cudaMemcpyDeviceToHost);

            // write to files time instant
            velocity_out << timeInstant << '\n';
            lift_drag_out << timeInstant << '\n';

            // write ux
            for (int i = 0; i < totalSize; ++i)
            {
                velocity_out << host_ux[i] << ' ';
            }
            velocity_out << '\n';
            // write uy
            for (int i = 0; i < totalSize; ++i)
            {
                velocity_out << host_uy[i] << ' ';
            }
            velocity_out << '\n';

            // lift and drag
            if (problemType == 2)
            {
                lift_drag_out << drag << ' ' << lift << '\n';
            }

            // print to console
            std::cout << "Time step: " << timeInstant << '\n';
        }

        // advance time
        timeInstant++;
    }
    cudaFreeHost(host_ux);
    cudaFreeHost(host_uy);
}
