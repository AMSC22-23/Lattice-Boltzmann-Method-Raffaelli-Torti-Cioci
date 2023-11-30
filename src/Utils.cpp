#ifndef UTILS_HPP
#define UTILS_HPP

#include "Cell.hpp"
#include <algorithm>
#include <execution>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <cmath>

/*
D2Q9 lattice
* 6 2 5
* 3 0 1
* 7 4 8
const struct
{
    int dimensions = 2;
    int velocity_directions = 9;
    std::unordered_map<int, float> weights = {{0, 4.0f / 9.0f},  {1, 1.0f / 9.0f},  {2, 1.0f / 9.0f},
                                              {3, 1.0f / 9.0f},  {4, 1.0f / 9.0f},  {5, 1.0f / 36.0f},
                                              {6, 1.0f / 36.0f}, {7, 1.0f / 36.0f}, {8, 1.0f / 36.0f}};
    // map of velocities for each direction
    std::unordered_map<int, std::vector<int>> velocities = {{0, {0, 0}},  {1, {1, 0}},   {2, {0, 1}},
                                                            {3, {-1, 0}}, {4, {0, -1}},  {5, {1, 1}},
                                                            {6, {-1, 1}}, {7, {-1, -1}}, {8, {1, -1}}};
    // map of opposite velocities for bounce back respect to x axis symmetry
    std::unordered_map<int, int> oppositeX = {{6, 7}, {7, 6}, {2, 4}, {4, 2}, {5, 8}, {8, 5}};
    // map of opposite velocities for bounce back respect to y axis symmetry
    std::unordered_map<int, int> oppositeY = {{5, 6}, {6, 5}, {1, 3}, {3, 1}, {7, 8}, {8, 7}};
    // map of opposite velocities for bounce back respect to both axis symmetry
    std::unordered_map<int, int> opposite = {{1, 3}, {3, 1}, {2, 4}, {4, 2}, {5, 7}, {7, 5}, {6, 8}, {8, 6}};

} D2Q9;

// D3Q27 lattice
const struct
{
    int dimensions = 3;

    int velocity_directions = 27;

    std::unordered_map<int, float> weights = {
        {0, 8.0f / 27.0f},   {1, 2.0f / 27.0f},  {2, 2.0f / 27.0f},  {3, 2.0f / 27.0f},  {4, 2.0f / 27.0f},
        {5, 2.0f / 27.0f},   {6, 2.0f / 27.0f},  {7, 1.0f / 54.0f},  {8, 1.0f / 54.0f},  {9, 1.0f / 54.0f},
        {10, 1.0f / 54.0f},  {11, 1.0f / 54.0f}, {12, 1.0f / 54.0f}, {13, 1.0f / 54.0f}, {14, 1.0f / 54.0f},
        {15, 1.0f / 54.0f},  {16, 1.0f / 54.0f}, {17, 1.0f / 54.0f}, {18, 1.0f / 54.0f}, {19, 1.0f / 54.0f},
        {20, 1.0f / 54.0f},  {21, 1.0f / 54.0f}, {22, 1.0f / 54.0f}, {23, 1.0f / 54.0f}, {24, 1.0f / 216.0f},
        {25, 1.0f / 216.0f}, {26, 1.0f / 216.0f}};

    std::unordered_map<int, std::vector<int>> velocities = {
        {0, {0, 0, 0}},    {1, {1, 0, 0}},   {2, {0, 1, 0}},   {3, {0, 0, 1}},    {4, {-1, 0, 0}},
        {5, {0, -1, 0}},   {6, {0, 0, -1}},  {7, {1, 1, 0}},   {8, {1, -1, 0}},   {9, {-1, 1, 0}},
        {10, {-1, -1, 0}}, {11, {1, 0, 1}},  {12, {1, 0, -1}}, {13, {-1, 0, 1}},  {14, {-1, 0, -1}},
        {15, {0, 1, 1}},   {16, {0, 1, -1}}, {17, {0, -1, 1}}, {18, {0, -1, -1}}, {19, {1, 1, 1}},
        {20, {1, 1, -1}},  {21, {1, -1, 1}}, {22, {-1, 1, 1}}, {23, {-1, -1, 1}}, {24, {-1, -1, -1}},
        {25, {1, -1, -1}}, {26, {-1, 1, -1}}};

    // map of opposite velocities for bounce back respect to both axis symmetry
    std::unordered_map<int, int> opposite = {{1, 4},   {4, 1},   {2, 5},   {5, 2},   {3, 6},   {6, 3},   {7, 10},
                                             {10, 7},  {8, 11},  {11, 8},  {9, 12},  {12, 9},  {13, 16}, {16, 13},
                                             {14, 17}, {17, 14}, {15, 18}, {18, 15}, {19, 22}, {22, 19}, {20, 23},
                                             {23, 20}, {21, 24}, {24, 21}, {25, 26}};
} D3Q27;
*/

template <class T> class NDimensionalMatrix
{
  private:
    std::vector<T> data;
    std::vector<int> dimensions;

  public:
    NDimensionalMatrix(const std::vector<int> &dims)
    {
        this->dimensions = dims;

        // Calculate the total size of the matrix
        int totalSize = 1;
        for (int dim : dims)
        {
            totalSize *= dim;
        }

        // Initialize the data vector with the calculated size
        data.resize(totalSize);
    }

    const T &getElement(const std::vector<int> &indices) const
    {
        // Validate the number of indices
        if (indices.size() != dimensions.size())
        {
            throw std::runtime_error("Invalid number of indices");
        }

        // Calculate the flat index using a formula
        int flatIndex = 0;
        int multiplier = 1;
        for (int i = 0; i < dimensions.size(); ++i)
        {
            flatIndex += indices[i] * multiplier;
            multiplier *= dimensions[i];
        }

        // Return the reference to the element
        return data.at(flatIndex);
    }

    T &getMutableElement(const std::vector<int> &indices)
    {
        // Validate the number of indices
        if (indices.size() != dimensions.size())
        {
            throw std::runtime_error("Invalid number of indices");
        }

        // Calculate the flat index using a formula
        int flatIndex = 0;
        int multiplier = 1;
        for (int i = 0; i < dimensions.size(); ++i)
        {
            flatIndex += indices[i] * multiplier;
            multiplier *= dimensions[i];
        }

        // Return the reference to the element
        return data.at(flatIndex);
    }

    const std::vector<int> &getShape() const
    {
        return dimensions;
    }

    const int getTotalSize() const
    {
        return data.size();
    }

    const std::vector<int> getIndicesAtFlatIndex(int flatIndex) const
    {
        // Validate the flat index
        if (flatIndex >= data.size())
        {
            throw std::runtime_error("Invalid flat index");
        }

        // Calculate the indices using a formula
        std::vector<int> indices(dimensions.size());
        for (int i = dimensions.size() - 1; i >= 0; --i)
        {
            indices[i] = flatIndex % dimensions[i];
            flatIndex /= dimensions[i];
        }

        // Return the indices
        return indices;
    }

    T &getMutableElementAtFlatIndex(int index)
    {
        return data.at(index);
    }

    const T &getElementAtFlatIndex(int index) const
    {
        return data.at(index);
    }

    void setElement(const std::vector<int> &indices, const T &element)
    {
        // Validate the number of indices
        if (indices.size() != dimensions.size())
        {
            throw std::runtime_error("Invalid number of indices");
        }

        // Calculate the flat index using a formula
        int flatIndex = 0;
        int multiplier = 1;
        for (int i = 0; i < dimensions.size(); ++i)
        {
            flatIndex += indices[i] * multiplier;
            multiplier *= dimensions[i];
        }

        // Set the element
        data.at(flatIndex) = element;
    }

    void setElementAtFlatIndex(int index, const T &element)
    {
        data.at(index) = element;
    }

    NDimensionalMatrix() = default; // Default constructor
};

template <class T> float scalar_product_parallel(const std::vector<std::vector<T>> &vectors)
{
    std::vector<T> product_vector = vectors[0];

    for (int i = 1; i < vectors.size(); i++)
    {
        std::transform(std::execution::par, product_vector.begin(), product_vector.end(), vectors[i].begin(),
                       product_vector.begin(), [](T value1, T value2) { return value1 * value2; });
    }

    return std::accumulate(product_vector.begin(), product_vector.end(), 0.0f);
}

template <class T> std::vector<T> scalar_vector_product_parallel(const T &scalar, const std::vector<T> &vector)
{
    std::vector<T> product(vector.size());
    std::transform(std::execution::par, vector.begin(), vector.end(), product.begin(),
                   [&scalar](const T &value) { return scalar * value; });
    return product;
}

template <class T> float vector_modulus_parallel(const std::vector<T> &vector)
{
    T modulus_squared = std::transform_reduce(std::execution::par, vector.begin(), vector.end(), T(0), std::plus<>(),
                                              [](T value) { return value * value; });

    return std::sqrt(modulus_squared);
}

#endif // UTILS_HPP
