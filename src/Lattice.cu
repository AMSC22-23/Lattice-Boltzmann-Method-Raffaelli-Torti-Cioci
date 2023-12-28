#include "GpuSimulation.cuh"
#include "Lattice.hpp"
#include <iostream>
#include <omp.h>

int __calculateBoundary(const std::vector<int> &inputBoundary)
{
    // calculate host boundary with streaming conventions.
    // 0 is no boundary, 1 is right, 2 is up, 3 is left, 4 is down, 5 is up-right, 6 is up-left, 7 is down-left, 8 is
    // down-right

    if (inputBoundary.at(0) == 0 && inputBoundary.at(1) == 0)
        return 0;
    else if (inputBoundary.at(0) == 1 && inputBoundary.at(1) == 0)
        return 1;
    else if (inputBoundary.at(0) == 0 && inputBoundary.at(1) == -1)
        return 2;
    else if (inputBoundary.at(0) == -1 && inputBoundary.at(1) == 0)
        return 3;
    else if (inputBoundary.at(0) == 0 && inputBoundary.at(1) == 1)
        return 4;
    else if (inputBoundary.at(0) == 1 && inputBoundary.at(1) == -1)
        return 5;
    else if (inputBoundary.at(0) == -1 && inputBoundary.at(1) == -1)
        return 6;
    else if (inputBoundary.at(0) == -1 && inputBoundary.at(1) == 1)
        return 7;
    else if (inputBoundary.at(0) == 1 && inputBoundary.at(1) == 1)
        return 8;
    else
        return -1;
}

Lattice::Lattice(const std::string &filename)
{
    // Read the lattice from the file
    std::ifstream file;
    file.open(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file");
    }

    // read type of problem
    file >> problemType;
    file.get(); // Skip the newline

    // Read the number of cells in each dimension until newline
    std::vector<int> shape;
    int dimensions = 0;
    while (file.peek() != '\n')
    {
        int numCells;
        file >> numCells;
        shape.push_back(numCells);
        ++dimensions;
    }
    if (dimensions == 2)
    {
        structure = Structure::D2Q9;
    }
    /*
    else if (dimensions == 3)
    {
        structure = Structure::D3Q27;
    }
    */
    else
    {
        throw std::runtime_error("Invalid number of dimensions");
    }
    file.get(); // Skip the newline

    // Read Reynolds number
    float reynolds;
    file >> reynolds;

    // Read the simulation time
    float simulationTime;
    file >> simulationTime;

    if (problemType == 1)
    {
        file >> uLid;
    }
    file.get(); // Skip the newline

    // calculate simulation parameters
    sigma = 10.0 * shape.at(0);
    omP = 1.0 / (0.5 + 3.0 * uLid * shape.at(0) / reynolds);
    omM = 1.0 / (1.0 / (12.0 * uLid * shape.at(0) / reynolds) + 0.5);
    maxIt = (int)std::round(simulationTime * shape.at(0) / uLid);

    // Initialize the cells
    cells = NDimensionalMatrix<Cell>(shape);

    // Read the obstacles : for each newline, read the coordinates of the obstacle
    NDimensionalMatrix<bool> obstacles(shape);
    for (int i = 0; i < obstacles.getTotalSize(); ++i)
    {
        obstacles.setElementAtFlatIndex(i, false);
    }
    while (file.peek() != EOF)
    {
        std::vector<int> indices;
        for (int i = 0; i < dimensions; ++i)
        {
            int index;
            file >> index;
            indices.push_back(index);
        }
        obstacles.setElement(indices, true);
        file.get(); // Skip the newline
    }

    //  Initialize the cells one by one
    std::vector<float> f;
    std::vector<int> boundary;
    std::vector<int> indices;
    for (int i = 0; i < cells.getTotalSize(); ++i)
    {
        indices = cells.getIndicesAtFlatIndex(i);
        const bool &obstacle = obstacles.getElementAtFlatIndex(i);

        boundary.clear();
        f.clear();

        for (int k = 0; k < dimensions; ++k)
        {
            const int indexOfCurrDimension = indices.at(k);
            const int lenghtOfCurrDimension = shape.at(k);
            if (indexOfCurrDimension == 0)
            {
                boundary.push_back(-1);
            }
            else if (indexOfCurrDimension == lenghtOfCurrDimension - 1)
            {
                boundary.push_back(1);
            }
            else
            {
                boundary.push_back(0);
            }

            /*
            // check if there is an obstacle in the adjacent cell
            if (indexOfCurrDimension > 0)
            {
                std::vector<int> adjacentIndices = indices;
                adjacentIndices.at(i) -= 1;
                if (obstacles.getElement(adjacentIndices))
                {
                    boundary.at(i) = -1;
                }
            }

            if (indexOfCurrDimension < lenghtOfCurrDimension - 1)
            {
                std::vector<int> adjacentIndices = indices;
                adjacentIndices.at(i) += 1;
                if (obstacles.getElement(adjacentIndices))
                {
                    boundary.at(i) = 1;
                }
            }
            */
        }

        // f is 1
        for (int j = 0; j < structure.velocity_directions; ++j)
        {
            f.push_back(1);
        }

        cells.setElementAtFlatIndex(i, Cell(structure, boundary, obstacle, f));
    }

    // Close the file
    file.close();
}

void Lattice::simulate(std::ofstream &file)
{
    const float temp = 2.0 * sigma * sigma;
    const float halfOmpOmmSub = 0.5 * (omP - omM);
    const float halfOmpOmmSum = 0.5 * (omP + omM);
    while (timeInstant <= maxIt)
    {
        const float uLidNow = uLid * (1.0 - std::exp(-static_cast<double>(timeInstant * timeInstant) / temp));
// update cells
#pragma omp parallel
        {
#pragma omp for
            for (int j = 0; j < cells.getTotalSize(); ++j)
            {
                if (timeInstant != 0)
                {
                    cells.getElementAtFlatIndex(j).zouHe();
                }
                cells.getElementAtFlatIndex(j).updateMacro(structure);
                cells.getElementAtFlatIndex(j).setInlets(structure, uLidNow, problemType);
                cells.getElementAtFlatIndex(j).equilibriumCollision(structure, omP, halfOmpOmmSum, halfOmpOmmSub);
            }
#pragma omp for
            for (int j = 0; j < cells.getTotalSize(); ++j)
            {
                cells.getElementAtFlatIndex(j).streaming(*this, cells.getIndicesAtFlatIndex(j));
            }
        }

        // write to file every maxIt/100 time steps
        if (timeInstant % (maxIt / 100) == 0)
        {
            // write to file time instant
            file << timeInstant << '\n';

            // loop dimensions
            for (int i = 0; i < structure.dimensions; ++i)
            {
                // write to file macroU
                for (int j = 0; j < cells.getTotalSize(); ++j)
                {
                    file << cells.getElementAtFlatIndex(j).getMacroU().at(i) << ' ';
                }
                file << '\n';
            }
            // print to console every 100 time steps
            std::cout << "Time step: " << timeInstant << '\n';
        }

        // advance time
        timeInstant++;
    }
}

/// @brief supports only 2D lattice
void Lattice::simulateGpu(std::ofstream &file)
{
    const int nx = cells.getShape().at(0);
    const int ny = cells.getShape().at(1);

    // device allocations

    float *host_f, *host_new_f, *host_rho, *host_ux, *host_uy, *dev_f, *dev_new_f, *dev_rho, *dev_ux, *dev_uy;
    int *host_boundary, *dev_boundary;
    bool *host_obstacle, *dev_obstacle;

    // allocate memory on device
    cudaMalloc((void **)&dev_f, cells.getTotalSize() * 9 * sizeof(float));
    cudaMalloc((void **)&dev_new_f, cells.getTotalSize() * 9 * sizeof(float));
    cudaMalloc((void **)&dev_rho, cells.getTotalSize() * sizeof(float));
    cudaMalloc((void **)&dev_ux, cells.getTotalSize() * sizeof(float));
    cudaMalloc((void **)&dev_uy, cells.getTotalSize() * sizeof(float));
    cudaMalloc((void **)&dev_boundary, cells.getTotalSize() * sizeof(int));
    cudaMalloc((void **)&dev_obstacle, cells.getTotalSize() * sizeof(bool));

    // allocate memory on host
    cudaMallocHost((void **)&host_f, cells.getTotalSize() * 9 * sizeof(float));
    cudaMallocHost((void **)&host_new_f, cells.getTotalSize() * 9 * sizeof(float));
    cudaMallocHost((void **)&host_rho, cells.getTotalSize() * sizeof(float));
    cudaMallocHost((void **)&host_ux, cells.getTotalSize() * sizeof(float));
    cudaMallocHost((void **)&host_uy, cells.getTotalSize() * sizeof(float));
    cudaMallocHost((void **)&host_boundary, cells.getTotalSize() * sizeof(int));
    cudaMallocHost((void **)&host_obstacle, cells.getTotalSize() * sizeof(bool));

#pragma omp parallel for
    // set host data
    for (int i = 0; i < cells.getTotalSize(); ++i)
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
        host_boundary[i] = __calculateBoundary(boundary);
        host_obstacle[i] = obstacle;
    }

    // copy host data to device
    cudaMemcpy(dev_f, host_f, cells.getTotalSize() * 9 * sizeof(float),
               cudaMemcpyHostToDevice);
    cudaMemcpy(dev_new_f, host_new_f, cells.getTotalSize() * 9 * sizeof(float),
               cudaMemcpyHostToDevice);
    cudaMemcpy(dev_rho, host_rho, cells.getTotalSize() * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_ux, host_ux, cells.getTotalSize() * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_uy, host_uy, cells.getTotalSize() * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_boundary, host_boundary, cells.getTotalSize() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_obstacle, host_obstacle, cells.getTotalSize() * sizeof(bool), cudaMemcpyHostToDevice);

    // free host memory
    cudaFreeHost(host_f);
    cudaFreeHost(host_new_f);
    cudaFreeHost(host_rho);
    cudaFreeHost(host_boundary);
    cudaFreeHost(host_obstacle);

    const dim3 threadsPerBlock(24, 24);
    const dim3 numBlocks(ceil(cells.getShape().at(0) / 24.0), ceil(cells.getShape().at(1) / 24.0));
    // loop
    const float temp = 2.0 * sigma * sigma;
    const float halfOmpOmmSub = 0.5 * (omP - omM);
    const float halfOmpOmmSum = 0.5 * (omP + omM);
    while (timeInstant <= maxIt)
    {
        const float uLidNow = uLid * (1.0 - std::exp(-static_cast<double>(timeInstant * timeInstant) / temp));
        GpuSimulation::step1<<<numBlocks, threadsPerBlock>>>(nx, ny, timeInstant, problemType, uLidNow, omP, halfOmpOmmSum, halfOmpOmmSub, dev_f,
                                                             dev_new_f, dev_rho, dev_ux, dev_uy, dev_boundary,
                                                             dev_obstacle);
        GpuSimulation::step2<<<numBlocks, threadsPerBlock>>>(nx, ny, dev_f, dev_new_f, dev_boundary, dev_obstacle);

        // write to file every maxIt/100 time steps
        if (timeInstant % (maxIt / 100) == 0)
        {
            // copy ux and uy to host
            cudaMemcpy(host_ux, dev_ux, cells.getTotalSize() * sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(host_uy, dev_uy, cells.getTotalSize() * sizeof(float), cudaMemcpyDeviceToHost);
            // write to file time instant
            file << timeInstant << '\n';

            // write ux
                            for (int i = 0; i < cells.getTotalSize(); ++i)
                {
                    file << host_ux[i] << ' ';
                            }
            file << '\n';
            // write uy
                            for (int i = 0; i < cells.getTotalSize(); ++i)
                {
                    file << host_uy[i] << ' ';
                            }
            file << '\n';
            // print to console
            std::cout << "Time step: " << timeInstant << '\n';
        }

        // advance time
        timeInstant++;
    }
    cudaFreeHost(host_ux);
    cudaFreeHost(host_uy);
}

Cell &Lattice::getCellAtIndices(const std::vector<int> &indices)
{
    return cells.getElement(indices);
}

Cell &Lattice::getCellAtIndices(const int x, const int y)
{
    return cells.getElement(x, y);
}

Cell &Lattice::getCellAtIndices(const int *indices)
{
    if (structure.dimensions == 2)
    {
        return cells.getElement(indices[0], indices[1]);
    }
    else if (structure.dimensions == 3)
    {
        return cells.getElement(indices[0], indices[1], indices[2]);
    }
    else
    {
        throw std::runtime_error("Invalid number of dimensions");
    }
}

Cell &Lattice::getCellAtIndices(const int x, const int y, const int z)
{
    return cells.getElement(x, y, z);
}

const std::vector<int> Lattice::getShape() const
{
    return cells.getShape();
}

bool Lattice::isLid() const
{
    return problemType == 1;
}

const Structure &Lattice::getStructure() const
{
    return structure;
}