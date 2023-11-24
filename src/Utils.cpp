#ifndef UTILS_HPP
#define UTILS_HPP

#include <stdexcept>
#include <vector>

// D2Q9 lattice
static struct
{
    unsigned int length = 9;
    float weights[9] = {4.0f / 9.0f,  1.0f / 9.0f,  1.0f / 9.0f,  1.0f / 9.0f, 1.0f / 9.0f,
                        1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f};
    int velocities[9][2] = {{0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
} D2Q9;

// D3Q27 lattice
static struct
{
    unsigned int length = 27;
    float weights[27] = {8.0f / 27.0f,  2.0f / 27.0f,  2.0f / 27.0f, 2.0f / 27.0f, 2.0f / 27.0f,  2.0f / 27.0f,
                         2.0f / 27.0f,  1.0f / 54.0f,  1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f,  1.0f / 54.0f,
                         1.0f / 54.0f,  1.0f / 54.0f,  1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f,  1.0f / 54.0f,
                         1.0f / 54.0f,  1.0f / 54.0f,  1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 216.0f, 1.0f / 216.0f,
                         1.0f / 216.0f, 1.0f / 216.0f, 1.0f / 216.0f};
    int velocities[27][3] = {{0, 0, 0},   {1, 0, 0},  {0, 1, 0},   {0, 0, 1},    {-1, 0, 0},  {0, -1, 0}, {0, 0, -1},
                             {1, 1, 0},   {1, -1, 0}, {-1, 1, 0},  {-1, -1, 0},  {1, 0, 1},   {1, 0, -1}, {-1, 0, 1},
                             {-1, 0, -1}, {0, 1, 1},  {0, 1, -1},  {0, -1, 1},   {0, -1, -1}, {1, 1, 1},  {1, 1, -1},
                             {1, -1, 1},  {-1, 1, 1}, {-1, -1, 1}, {-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}};
} D3Q27;

template <typename T> class NDimensionalMatrix
{
  private:
    std::vector<T> data;
    std::vector<unsigned int> dimensions;

  public:
    NDimensionalMatrix(const std::vector<unsigned int> &dims)
    {
        this->dimensions = dims;

        // Calculate the total size of the matrix
        unsigned int totalSize = 1;
        for (unsigned int dim : dims)
        {
            totalSize *= dim;
        }

        // Initialize the data vector with the calculated size
        data.resize(totalSize);
    }

    T &getElement(const std::vector<unsigned int> &indices)
    {
        // Validate the number of indices
        if (indices.size() != dimensions.size())
        {
            throw std::runtime_error("Invalid number of indices");
        }

        // Calculate the flat index using a formula
        unsigned int flatIndex = 0;
        unsigned int multiplier = 1;
        for (unsigned int i = 0; i < dimensions.size(); ++i)
        {
            flatIndex += indices[i] * multiplier;
            multiplier *= dimensions[i];
        }

        // Return the reference to the element
        return data[flatIndex];
    }

    const std::vector<unsigned int> &getShape()
    {
        return dimensions;
    }

    NDimensionalMatrix() = default; // Default constructor
};

#endif // UTILS_HPP
