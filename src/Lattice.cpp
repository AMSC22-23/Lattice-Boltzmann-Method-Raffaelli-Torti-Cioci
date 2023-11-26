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
    for (bool obstacle : obstacles)
    {
        obstacle = false;
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
        obstacles.getElement(indices) = true;
        file.get(); // Skip the newline
    }

    // Initialize the cells
    for (auto it = cells.begin(); it != cells.end(); ++it)
    {
        std::vector<int> indices = it.getIndices();
        bool obstacle = obstacles.getElement(indices);

        // calculate boundary based on edges of the lattice and adjacent obstacles
        std::vector<int> boundary;
        for (int i = 0; i < dimensions; ++i)
        {
            int index = indices.at(i);
            int numCells = shape.at(i);
            if (index == 0)
            {
                boundary.push_back(-1);
            }
            else if (index == numCells - 1)
            {
                boundary.push_back(1);
            }
            else
            {
                boundary.push_back(0);
            }

            // check if there is an obstacle in the adjacent cell
            if (index > 0)
            {
                std::vector<int> adjacentIndices = indices;
                adjacentIndices.at(i) -= 1;
                if (obstacles.getElement(adjacentIndices))
                {
                    boundary.at(i) = -1;
                }
            }

            if (index < numCells - 1)
            {
                std::vector<int> adjacentIndices = indices;
                adjacentIndices.at(i) += 1;
                if (obstacles.getElement(adjacentIndices))
                {
                    boundary.at(i) = 1;
                }
            }
        }

        it.emplace(indices, Cell(indices, boundary, obstacle));
    }

    // Close the file
    file.close();
}

Lattice::~Lattice()
{
}

void Lattice::update(const float deltaTime)
{
    for (int i = 0; i < cells.getShape().at(0); i++)
    {
        for (int j = 0; j < cells.getShape().at(1); j++)
        {
            Cell &cell = cells.getElement({i, j});
            if (!cell.isObstacle())
            {
                cell.update(deltaTime, *this, {i, j});
            }
        }
    }
}

Cell &Lattice::getCellAtIndex(std::vector<int> index)
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