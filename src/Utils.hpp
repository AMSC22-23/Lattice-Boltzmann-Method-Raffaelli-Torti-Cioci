#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>

// D1Q3 lattice
static struct {
    unsigned int length = 3;
    float weights[3] = {2.0f/3.0f, 1.0f/6.0f, 1.0f/6.0f};
    int velocities[3] = {0, 1, -1};
} D1Q3;

// D2Q9 lattice
static struct {
    unsigned int length = 9;
    float weights[9] = {4.0f/9.0f, 1.0f/9.0f, 1.0f/9.0f, 1.0f/9.0f,
                                         1.0f/9.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f};
    int velocities[9][2] = {{0, 0}, {1, 0}, {0, 1}, {-1, 0},
                                             {0, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
} D2Q9;

// D3Q27 lattice
static struct {
    unsigned int length = 27;
    float weights[27] = {8.0f/27.0f, 2.0f/27.0f, 2.0f/27.0f, 2.0f/27.0f,
                        2.0f/27.0f, 2.0f/27.0f, 2.0f/27.0f, 1.0f/54.0f, 1.0f/54.0f,
                        1.0f/54.0f, 1.0f/54.0f, 1.0f/54.0f, 1.0f/54.0f, 1.0f/54.0f,
                        1.0f/54.0f, 1.0f/54.0f, 1.0f/54.0f, 1.0f/54.0f, 1.0f/54.0f,
                        1.0f/54.0f, 1.0f/54.0f, 1.0f/54.0f, 1.0f/216.0f, 1.0f/216.0f,
                        1.0f/216.0f, 1.0f/216.0f, 1.0f/216.0f};
    int velocities[27][3] = {
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
        {-1, 0, 0}, {0, -1, 0}, {0, 0, -1}, {1, 1, 0}, {1, -1, 0},
        {-1, 1, 0}, {-1, -1, 0}, {1, 0, 1}, {1, 0, -1}, {-1, 0, 1},
        {-1, 0, -1}, {0, 1, 1}, {0, 1, -1}, {0, -1, 1}, {0, -1, -1},
        {1, 1, 1}, {1, 1, -1}, {1, -1, 1}, {-1, 1, 1}, {-1, -1, 1},
        {-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}};
} D3Q27;

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
