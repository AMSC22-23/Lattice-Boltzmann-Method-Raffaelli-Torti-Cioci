#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>



// D1Q3 lattice
typedef struct {
    static constexpr unsigned int length = 3;
    static constexpr float weights[3] = {2.0f/3.0f, 1.0f/6.0f, 1.0f/6.0f};
    static constexpr int velocities[3] = {0, 1, -1};
} D1Q3;

// D2Q9 lattice
typedef struct {
    static constexpr unsigned int length = 9;
    static constexpr float weights[9] = {4.0f/9.0f, 1.0f/9.0f, 1.0f/9.0f, 1.0f/9.0f,
                                         1.0f/9.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f};
    static constexpr int velocities[9][2] = {{0, 0}, {1, 0}, {0, 1}, {-1, 0},
                                             {0, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
} D2Q9;

template <typename T>
class NDimensionalMatrix {
private:
    std::vector<T> data;
    std::vector<unsigned int> dimensions;

public:
    NDimensionalMatrix(const std::vector<unsigned int>& dims);
    T& getElement(const std::vector<unsigned int>& indices);
    unsigned int getDimensionLength(const unsigned int dim);
};

#endif // UTILS_HPP
