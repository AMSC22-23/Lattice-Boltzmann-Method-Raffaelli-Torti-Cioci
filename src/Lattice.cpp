#include "Lattice.hpp"

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
    if (dimensions != 2 && dimensions != 3)
    {
        throw std::runtime_error("Invalid number of dimensions");
    }
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

    // TODO Initialize the cells
    for (int i = 0; i < cells.getTotalSize(); ++i)
    {
        // TODO cells.setElementAtFlatIndex(i, Cell(argomentiDiCostruttoreDiCell));
        std::vector<int> indices = cells.getIndicesAtFlatIndex(i);
        bool obstacle = obstacles.getElement(indices);

        // calculate boundary based on edges of the lattice and adjacent obstacles
        std::vector<int> boundary;

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
                    boundary.push_back(1);
                    break;
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
                    boundary.push_back(-1);
                    break;
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
    }

    // Close the file
    file.close();
}

void Lattice::update(const float deltaTime)
{
    for (int i = 0; i < cells.getShape().at(0); i++)
    {
        for (int j = 0; j < cells.getShape().at(1); j++)
        {
            Cell cell = cells.getElement({i, j});
            if (!cell.isObstacle())
            {
                cell.update(deltaTime, *this, {i, j});
            }
        }
    }
}

const Cell Lattice::getCellAtIndex(std::vector<int> index)
{
    return cells.getElement(index);
}

std::vector<int> Lattice::getShape()
{
    return cells.getShape();
}

bool Lattice::isLid()
{
    return lid;
}