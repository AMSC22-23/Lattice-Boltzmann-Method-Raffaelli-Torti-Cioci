#ifndef UTILS_HPP
#define UTILS_HPP

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
    int calculateFlatIndex(const std::vector<int> &indices) const
    {
        // Validate the number of indices
        if (indices.size() != dimensions.size())
        {
            throw std::runtime_error("Invalid number of indices");
        }

        // Calculate the flat index using a formula
        int flatIndex = 0;
        for (int i = indices.size() - 1; i >= 0; --i)
        {
            flatIndex = flatIndex * dimensions[i] + indices[i];
        }

        // Return the flat index
        return flatIndex;
    }

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

    T &getElement(const std::vector<int> &indices)
    {
        return data.at(calculateFlatIndex(indices));
    }

    const T getConstCopy(const int x, const int y) const 
    {
        return data.at(x + y * dimensions[0]);
    }

    T &getElement(const int x, const int y)
    {
        return data.at(x + y * dimensions[0]);
    }

    T &getElement(const int x, const int y, const int z)
    {
        return data.at(x + y * dimensions[0] + z * dimensions[0] * dimensions[1]);
    }

    const std::vector<int> &getShape() const
    {
        return dimensions;
    }

    int getTotalSize() const
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
        for (int i = 0; i < dimensions.size(); ++i)
        {
            indices[i] = flatIndex % dimensions[i]; //
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
        const T &returnValue = data.at(index);
        return returnValue;
    }

    void setElement(const std::vector<int> &indices, const T &element)
    {
        data.at(calculateFlatIndex(indices)) = element;
    }

    void setElementAtFlatIndex(const int index, const T &element)
    {
        data.at(index) = element;
    }

    NDimensionalMatrix() = default; // Default constructor
};

#endif // UTILS_HPP
