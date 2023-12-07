#ifndef UTILS_HPP
#define UTILS_HPP

#include "Cell.hpp"
#include <algorithm>
#include <cmath>
#include <execution>
#include <stdexcept>
#include <unordered_map>
#include <vector>

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

        const auto &returnValue = data.at(flatIndex);
        // Return the reference to the element
        return returnValue;
    }

    T &getElement(const std::vector<int> &indices)
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
        return const_cast<T &>(const_cast<const NDimensionalMatrix *>(this)->getElement(indices));
    }

    /*T &getMutableElement(const std::vector<int> &indices)
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
    }*/

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

    T &getElementAtFlatIndex(const int index)
    {
        return const_cast<T &>(const_cast<const NDimensionalMatrix *>(this)->getElementAtFlatIndex(index));
    }

    const T &getElementAtFlatIndex(const int index) const
    {
        const auto &returnValue = data.at(index);
        return returnValue;
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

    void setElementAtFlatIndex(const int index, const T &element)
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
