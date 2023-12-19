#include "Lattice.hpp"
#include <iostream>

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

        for (int i = 0; i < dimensions; ++i)
        {
            const int indexOfCurrDimension = indices.at(i);
            const int lenghtOfCurrDimension = shape.at(i);
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
        for (int i = 0; i < structure.velocity_directions; ++i)
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
    std::vector<int> indices;

    // loop
    while (timeInstant <= maxIt)
    {
        const float uLidNow =
            uLid * (1.0 - std::exp(-static_cast<double>(timeInstant * timeInstant) / (2.0 * sigma * sigma)));
        // update cells
        for (int j = 0; j < cells.getTotalSize(); ++j)
        {
            cells.getElementAtFlatIndex(j).updateMacro(structure);
        }
        for (int j = 0; j < cells.getTotalSize(); ++j)
        {
            cells.getElementAtFlatIndex(j).updateFeq(structure);
        }
        for (int j = 0; j < cells.getTotalSize(); ++j)
        {
            indices = cells.getIndicesAtFlatIndex(j);
            cells.getElementAtFlatIndex(j).collisionStreaming(*this, indices, omP, omM);
        }
        for (int j = 0; j < cells.getTotalSize(); ++j)
        {
            cells.getElementAtFlatIndex(j).setInlets(structure, uLidNow, problemType);
        }
        for (int j = 0; j < cells.getTotalSize(); ++j)
        {
            cells.getElementAtFlatIndex(j).zouHe();
        }

        // write to file every 100 time steps
        if (timeInstant % 100 == 0)
        {
            // write to file time instant
            file << timeInstant << '\n';

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
            // print to console every 100 time steps
            std::cout << "Time step: " << timeInstant << '\n';
        }

        // advance time
        timeInstant++;
    }
}

Cell &Lattice::getCellAtIndices(std::vector<int> indices)
{
    return cells.getElement(indices);
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