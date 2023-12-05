#include "Structure.hpp"

// Definition of the copy assignment operator
Structure &Structure::operator=(const Structure &other)
{
    if (this != &other)
    {
        dimensions = other.dimensions;
        velocity_directions = other.velocity_directions;
        weights = other.weights;
        velocities = other.velocities;
        opposite = other.opposite;
    }
    return *this;
}

Structure::Structure(int dimensions, int velocity_directions, std::vector<float> weights,
                     std::vector<std::vector<int>> velocities, std::vector<int> opposite)
    : dimensions(dimensions), velocity_directions(velocity_directions), weights(weights), velocities(velocities),
      opposite(opposite)
{
}

// Initialization of D2Q9
Structure Structure::D2Q9 = Structure(2, 9,
                                      {{4.0 / 9.0},
                                       {1.0 / 9.0},
                                       {1.0 / 9.0},
                                       {1.0 / 9.0},
                                       {1.0 / 9.0},
                                       {1.0 / 36.0},
                                       {1.0 / 36.0},
                                       {1.0 / 36.0},
                                       {1.0 / 36.0}},
                                      {{0, 0},
                                       {1, 0},
                                       {0, 1},
                                       {-1, 0},
                                       {0, -1},
                                       {1, 1},
                                       {-1, 1},
                                       {-1, -1},
                                       {1, -1}},
                                      {0, 3, 4, 1, 2, 7, 8, 5, 6});

// TODO 3D
/*
// Initialization of D3Q27
Structure Structure::D3Q27 = Structure(
    3, 27, {{0, 8.0f / 27.0f},   {1, 2.0f / 27.0f},  {2, 2.0f / 27.0f},  {3, 2.0f / 27.0f},  {4, 2.0f / 27.0f},
            {5, 2.0f / 27.0f},   {6, 2.0f / 27.0f},  {7, 1.0f / 54.0f},  {8, 1.0f / 54.0f},  {9, 1.0f / 54.0f},
            {10, 1.0f / 54.0f},  {11, 1.0f / 54.0f}, {12, 1.0f / 54.0f}, {13, 1.0f / 54.0f}, {14, 1.0f / 54.0f},
            {15, 1.0f / 54.0f},  {16, 1.0f / 54.0f}, {17, 1.0f / 54.0f}, {18, 1.0f / 54.0f}, {19, 1.0f / 54.0f},
            {20, 1.0f / 54.0f},  {21, 1.0f / 54.0f}, {22, 1.0f / 54.0f}, {23, 1.0f / 54.0f}, {24, 1.0f / 216.0f},
            {25, 1.0f / 216.0f}, {26, 1.0f / 216.0f}},
    {{0, {0, 0, 0}},     {1, {1, 0, 0}},    {2, {0, 1, 0}},    {3, {0, 0, 1}},   {4, {-1, 0, 0}},   {5, {0, -1, 0}},
     {6, {0, 0, -1}},    {7, {1, 1, 0}},    {8, {1, -1, 0}},   {9, {-1, 1, 0}},  {10, {-1, -1, 0}}, {11, {1, 0, 1}},
     {12, {1, 0, -1}},   {13, {-1, 0, 1}},  {14, {-1, 0, -1}}, {15, {0, 1, 1}},  {16, {0, 1, -1}},  {17, {0, -1, 1}},
     {18, {0, -1, -1}},  {19, {1, 1, 1}},   {20, {1, 1, -1}},  {21, {1, -1, 1}}, {22, {-1, 1, 1}},  {23, {-1, -1, 1}},
     {24, {-1, -1, -1}}, {25, {1, -1, -1}}, {26, {-1, 1, -1}}},
    {{0, 0},   {1, 4},   {4, 1},   {2, 5},   {5, 2},   {3, 6},   {6, 3},   {7, 10},  {10, 7},
     {8, 11},  {11, 8},  {9, 12},  {12, 9},  {13, 16}, {16, 13}, {14, 17}, {17, 14}, {15, 18},
     {18, 15}, {19, 22}, {22, 19}, {20, 23}, {23, 20}, {21, 24}, {24, 21}, {25, 26}});
*/