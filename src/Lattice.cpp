#include "Lattice.hpp"
#include <iostream>
#include <omp.h>

std::vector<int> boundaryCalculator2d(const std::vector<int> &indices, const NDimensionalMatrix<bool> &obstacleMatrix)
{
    const int &col = indices.at(0);
    const int &row = indices.at(1);
    const int &nx = obstacleMatrix.getShape().at(0);
    const int &ny = obstacleMatrix.getShape().at(1);
    std::vector<int> boundary_here(4);

    // determine boundary based on neighboring obstacles

    // horizontal
    if (col > 0 && obstacleMatrix.getConstCopy(col - 1, row))
    {
        boundary_here[0] = -1;
    }
    else if (col < nx - 1 && obstacleMatrix.getConstCopy(col + 1, row))
    {
        boundary_here[0] = 1;
    }
    else
    {
        boundary_here[0] = 0;
    }
    // vertical
    if (row > 0 && obstacleMatrix.getConstCopy(col, row - 1))
    {
        boundary_here[1] = -1;
    }
    else if (row < ny - 1 && obstacleMatrix.getConstCopy(col, row + 1))
    {
        boundary_here[1] = 1;
    }
    else
    {
        boundary_here[1] = 0;
    }
    // main diagonal
    if (col > 0 && row > 0 && obstacleMatrix.getConstCopy(col - 1, row - 1))
    {
        boundary_here[2] = -1;
    }
    else if (col < nx - 1 && row < ny - 1 && obstacleMatrix.getConstCopy(col + 1, row + 1))
    {
        boundary_here[2] = 1;
    }
    else
    {
        boundary_here[2] = 0;
    }
    // secondary diagonal
    if (col > 0 && row < ny - 1 && obstacleMatrix.getConstCopy(col - 1, row + 1))
    {
        boundary_here[3] = 1;
    }
    else if (col < nx - 1 && row > 0 && obstacleMatrix.getConstCopy(col + 1, row - 1))
    {
        boundary_here[3] = -1;
    }
    else
    {
        boundary_here[3] = 0;
    }

    return boundary_here;
}

Lattice::Lattice(std::ifstream &file_in, const int plotSteps) : plotSteps(plotSteps)
{
    // read type of problem
    file_in >> problemType;
    file_in.get(); // Skip the newline

    // Read the number of cells in each dimension until newline
    std::vector<int> shape;
    int dimensions = 0;
    while (file_in.peek() != '\n')
    {
        int numCells;
        file_in >> numCells;
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
    file_in.get(); // Skip the newline

    // Read Reynolds number
    float reynolds;
    file_in >> reynolds;

    // Read the simulation steps
    file_in >> maxIt;

    // set input U
    file_in >> uLid;

    file_in.get(); // Skip the newline

    // calculate simulation parameters
    float nu = uLid * shape.at(1) / reynolds;
    if (problemType == 2 || problemType == 3)
    {
        nu *= 2.0 / 3.0;
    }

    const float tau = 0.5 + nu * 3.0;
    // const float dt = reynolds * nu / (shape.at(1) * shape.at(1));
    sigma = 10.0 * shape.at(1);

    const float lambda_trt = 1.0 / 4.0;
    const float tau_minus = lambda_trt / (tau - 0.5) + 0.5;
    omP = 1.0 / tau;
    omM = 1.0 / tau_minus;

    // Initialize the cells
    cells = NDimensionalMatrix<Cell>(shape);

    // Read the obstacles : for each newline, read the coordinates of the obstacle
    NDimensionalMatrix<bool> obstacles(shape);
    for (int i = 0; i < obstacles.getTotalSize(); ++i)
    {
        obstacles.setElementAtFlatIndex(i, false);
    }
    std::vector<int> indices;
    while (file_in.peek() != EOF)
    {
        indices.clear();
        for (int i = 0; i < dimensions; ++i)
        {
            int index;
            file_in >> index;
            indices.push_back(index);
        }
        obstacles.setElement(indices, true);
        file_in.get(); // Skip the newline
    }

    //  Initialize the cells one by one
    std::vector<int> boundary;
    bool obstacle;
    for (int i = 0; i < cells.getTotalSize(); ++i)
    {
        indices = cells.getIndicesAtFlatIndex(i);

        if (problemType == 2)
        {
            obstacle = obstacles.getConstCopy(indices.at(0), indices.at(1));
            boundary = boundaryCalculator2d(indices, obstacles);
        }
        else
        {
            obstacle = false;
            boundary = std::vector<int>(4, 0);
        }

        cells.setElementAtFlatIndex(i, Cell(structure, boundary, obstacle, indices));
        cells.getElementAtFlatIndex(i).initEq(structure, problemType);
        cells.getElementAtFlatIndex(i).updateMacro(structure);
    }
}

void Lattice::simulate(std::ofstream &velocity_out, std::ofstream &lift_drag_out)
{
    const float temp = 2.0 * sigma * sigma;
    const float halfOmpOmmSub = 0.5 * (omP - omM);
    const float halfOmpOmmSum = 0.5 * (omP + omM);
    float drag, lift;

    while (timeInstant <= maxIt)
    {
        const float uLidNow = uLid * (1.0 - std::exp(-static_cast<double>(timeInstant * timeInstant) / temp));

#pragma omp parallel
        {
#pragma omp for
            // inlets, zouhe, bb, macro, eq coll
            for (int j = 0; j < cells.getTotalSize(); ++j)
            {
                cells.getElementAtFlatIndex(j).setInlets(*this, uLidNow);
                cells.getElementAtFlatIndex(j).zouHe(*this);
                cells.getElementAtFlatIndex(j).updateMacro(structure);
                cells.getElementAtFlatIndex(j).equilibriumCollision(structure, omP, halfOmpOmmSum, halfOmpOmmSub);
                cells.getElementAtFlatIndex(j).bounce_back_obstacle();
            }

#pragma omp for
            // streaming
            for (int j = 0; j < cells.getTotalSize(); ++j)
            {
                cells.getElementAtFlatIndex(j).streaming(*this);
            }

            // drag and lift if problemType == 2 and about to print
            if (problemType == 2 && timeInstant % (maxIt / plotSteps) == 0)
            {
#pragma omp single
                {
                    drag = 0;
                    lift = 0;
                }

#pragma omp for
                for (int j = 0; j < cells.getTotalSize(); ++j)
                {
                    cells.getElementAtFlatIndex(j).dragAndLift(drag, lift);
                }
            }
        }

        // write to files every maxIt/plotSteps time steps
        if (timeInstant % (maxIt / plotSteps) == 0)
        {
            // write to files time instant
            velocity_out << timeInstant << '\n';
            lift_drag_out << timeInstant << '\n';

            // loop dimensions
            for (int i = 0; i < structure.dimensions; ++i)
            {
                // write to file macroU
                for (int j = 0; j < cells.getTotalSize(); ++j)
                {
                    velocity_out << cells.getElementAtFlatIndex(j).getMacroU().at(i) << ' ';
                }
                velocity_out << '\n';
            }

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
}

#ifdef USE_CUDA
#include "GpuSimulation.cuh"
/// @brief supports only 2D lattice
void Lattice::simulateGpu(std::ofstream &velocity_out, std::ofstream &lift_drag_out)
{
    GpuSimulation::cudaCaller(cells, sigma, omP, omM, maxIt, uLid, problemType, plotSteps, velocity_out, lift_drag_out);
}
#endif

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
        return cells.getElement(indices[0], indices[1]);

    else if (structure.dimensions == 3)
        return cells.getElement(indices[0], indices[1], indices[2]);

    else
        throw std::runtime_error("Invalid number of dimensions");
}

Cell &Lattice::getCellAtIndices(const int x, const int y, const int z)
{
    return cells.getElement(x, y, z);
}

const Cell &Lattice::getCellAtIndex(const int index) const
{
    return cells.getElementAtFlatIndex(index);
}

const std::vector<int> Lattice::getShape() const
{
    return cells.getShape();
}

const Structure &Lattice::getStructure() const
{
    return structure;
}

float Lattice::getCloseRho(const std::vector<int> &indices)
{
    const int xLen = getShape().at(0) - 1;
    const int yLen = getShape().at(1) - 1;
    if (indices.at(0) == 0 && indices.at(1) == yLen) // bottom left corner
        return getCellAtIndices(indices.at(0) + 1, indices.at(1)).getRho();

    if (indices.at(0) == xLen && indices.at(1) == yLen) // bottom right corner
        return getCellAtIndices(indices.at(0) - 1, indices.at(1)).getRho();

    if (indices.at(0) == 0 && indices.at(1) == 0) // top left corner
        return getCellAtIndices(indices.at(0) + 1, indices.at(1)).getRho();

    if (indices.at(0) == xLen && indices.at(1) == 0) // top right corner
        return getCellAtIndices(indices.at(0) - 1, indices.at(1)).getRho();

    else
        throw std::runtime_error("Invalid indices");
}

int Lattice::getProblemType() const
{
    return problemType;
}
