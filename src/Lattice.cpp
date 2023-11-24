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

    // Read the number of cells in each dimension until newline
    std::vector<unsigned int> shape;
    unsigned int dimensions = 0;
    while (file.peek() != '\n')
    {
        unsigned int numCells;
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

    // Read the obstacles : for each newline, read the coordinates of the obstacle
    while (file.peek() != EOF)
    {
        std::vector<unsigned int> indices;
        for (unsigned int i = 0; i < dimensions; ++i)
        {
            unsigned int index;
            file >> index;
            indices.push_back(index);
        }
        file.get(); // Skip the newline
        cells.getElement(indices).setObstacle();
    }

    // Close the file
    file.close();
}

Lattice::~Lattice()
{
}

void Lattice::update(const float deltaTime)
{
    for (unsigned int i = 0; i < cells.getShape().at(0); i++)
    {
        for (unsigned int j = 0; j < cells.getShape().at(1); j++)
        {
            Cell &cell = cells.getElement({i, j});
            if (!cell.isObstacle())
            {
                cell.update(deltaTime, *this);
            }
        }
    }
}

Cell &Lattice::getCellAtIndex(std::vector<unsigned int> index)
{
    return cells.getElement(index);
}

std::vector<unsigned int> Lattice::getShape()
{
    return cells.getShape();
}