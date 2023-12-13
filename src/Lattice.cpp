#include "Lattice.hpp"
#include <random>
#include <iostream>

Lattice::Lattice(std::string filename)
{
    // Read the lattice from the file
    std::ifstream file;
    file.open(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file");
    }

    // read type of problem
    int problemType;
    file >> problemType;
    if (problemType == 1)
    {
        lid = true;
    }
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
    auto reynoldsNumber = 0.0f;
    file >> reynoldsNumber;

    // Read macroscopic velocity
    std::vector<float> macroU;
    for (int i = 0; i < dimensions; ++i)
    {
        float velocity;
        file >> velocity;
        macroU.push_back(velocity);
    }

    // Read viscosity
    auto mu = 0.0f;
    file >> mu;

    file.get(); // Skip the newline

    // Initialize the lattice
    cells = NDimensionalMatrix<Cell>(shape);

    NDimensionalMatrix<bool> obstacles(shape);
    // Read the obstacles : for each newline, read the coordinates of the obstacle
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

    // Set up a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.99, 1.01);

    //  Initialize the cells one by one
    std::vector<float> f;
    std::vector<int> boundary;
    std::vector<int> indices;

    //  Initialize the cells one by one
    for (int i = 0; i < cells.getTotalSize(); ++i)
    {
        indices = cells.getIndicesAtFlatIndex(i);
        const bool &obstacle = obstacles.getElementAtFlatIndex(i);

        boundary.clear();
        f.clear();

        for (int i = 0; i < dimensions; ++i)
        {
            int indexOfCurrDimension = indices.at(i);
            int lenghtOfCurrDimension = shape.at(i);
            if (indexOfCurrDimension == 0)
            {
                switch (i)
                {
                case 0:
                    boundary.push_back(-1);
                    break;
                case 1:
                case 2:
                    boundary.push_back(1);
                    break;
                default:
                    throw std::runtime_error("Invalid dimension");
                }
            }
            else if (indexOfCurrDimension == lenghtOfCurrDimension - 1)
            {
                switch (i)
                {
                case 0:
                    boundary.push_back(1);
                    break;
                case 1:
                case 2:
                    boundary.push_back(-1);
                    break;
                default:
                    throw std::runtime_error("Invalid dimension");
                }
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

        // f is 1 with 0.01 random noise
        for (int i = 0; i < structure.velocity_directions; ++i)
        {
            f.push_back(dis(gen));
        }

        // ! old constructor
        // cells.setElementAtFlatIndex(i, Cell(structure, boundary, obstacle, reynoldsNumber, shape.at(0), mu));
        cells.setElementAtFlatIndex(i, Cell(structure, boundary, obstacle, f));
    }

    // Close the file
    file.close();
}

void Lattice::update(const float deltaTime, std::ofstream &file)
{
    std::vector<int> indices;
    // update all cells part one
    for (int i = 0; i < cells.getTotalSize(); ++i)
    {
        if (!cells.getElementAtFlatIndex(i).isObstacle())
        {
            indices = cells.getIndicesAtFlatIndex(i);
            cells.getElementAtFlatIndex(i).update1(deltaTime, *this, indices);
        }
    }
    // update all cells part two
    for (int i = 0; i < cells.getTotalSize(); ++i)
    {
        if (!cells.getElementAtFlatIndex(i).isObstacle())
        {
            indices = cells.getIndicesAtFlatIndex(i);
            cells.getElementAtFlatIndex(i).update2(getShape(), indices, structure);
        }
    }
    // write to file time instant
    file << timeInstant << '\n';
    // write to file rho
    for (int i = 0; i < cells.getTotalSize(); ++i)
    {
        // TODO ignoring obstacles for now
        file << cells.getElementAtFlatIndex(i).getRho() << ' ';
    }
    file << '\n';

    // loop dimensions
    for (int i = 0; i < structure.dimensions; ++i)
    {
        // write to file macroU
        for (int j = 0; j < cells.getTotalSize(); ++j)
        {
            // TODO ignoring obstacles for now
            file << cells.getElementAtFlatIndex(j).getMacroU().at(i) << ' ';
        }
        file << '\n';
    }
    // advance time
    timeInstant++;
    // print every 50 time steps
    if (timeInstant % 50 == 0)
    {
        std::cout << "Time step: " << timeInstant << '\n';
    }
}

/*
const Cell &Lattice::getCellAtIndices(std::vector<int> indices) const
{
    return cells.getElement(indices);
}
*/

Cell &Lattice::getCellAtIndices(std::vector<int> indices)
{
    return cells.getElement(indices);
}

const std::vector<int> Lattice::getShape()
{
    return cells.getShape();
}

bool Lattice::isLid()
{
    return lid;
}

const Structure &Lattice::getStructure() const
{
    return structure;
}