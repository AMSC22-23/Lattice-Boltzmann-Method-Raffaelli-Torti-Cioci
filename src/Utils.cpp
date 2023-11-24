#include <vector>
#include "Utils.hpp"

template <typename T>
NDimensionalMatrix<T>::NDimensionalMatrix(const std::vector<unsigned int>& dims) {
    this->dimensions = dims;

    // Calculate the total size of the matrix
    unsigned int totalSize = 1;
    for (unsigned int dim : dims) {
        totalSize *= dim;
    }

    // Initialize the data vector with the calculated size
    data.resize(totalSize);
}

template <typename T>
T& NDimensionalMatrix<T>::getElement(const std::vector<unsigned int>& indices) {
    // Validate the number of indices
    if (indices.size() != dimensions.size()) {
        throw std::runtime_error("Invalid number of indices");
    }

    // Calculate the flat index using a formula
    unsigned int flatIndex = 0;
    unsigned int multiplier = 1;
    for (unsigned int i = 0; i < dimensions.size(); ++i) {
        flatIndex += indices[i] * multiplier;
        multiplier *= dimensions[i];
    }

    // Return the reference to the element
    return data[flatIndex];
}

template <typename T>
std::vector<unsigned int>& NDimensionalMatrix<T>::getShape() {
    return dimensions;
}
