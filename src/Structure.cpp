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


// Initialization of D3Q27
Structure Structure::D3Q27 = Structure(
    3, 27, {{8.0f / 27.0f},   {2.0f / 27.0f},  {2.0f / 27.0f},  {2.0f / 27.0f},  {2.0f / 27.0f},
            {2.0f / 27.0f},   {2.0f / 27.0f},  {1.0f / 54.0f},  {1.0f / 54.0f},  {1.0f / 54.0f},
            {1.0f / 54.0f},  {1.0f / 54.0f}, {1.0f / 54.0f}, {1.0f / 54.0f}, {1.0f / 54.0f},
            {1.0f / 54.0f},  {1.0f / 54.0f}, {1.0f / 54.0f}, {1.0f / 54.0f}, {1.0f / 54.0f},
            {1.0f / 54.0f},  {1.0f / 54.0f}, {1.0f / 54.0f}, {1.0f / 54.0f}, {1.0f / 216.0f},
            {1.0f / 216.0f}, {1.0f / 216.0f}},
    {{0, 0, 0},    {1, 0, 0},    {0, 1, 0},   {0, 0, 1},   {-1, 0, 0},  {0, -1, 0},
     {0, 0, -1},   {1, 1, 0},    {1, -1, 0},  {-1, 1, 0},  {-1, -1, 0}, {1, 0, 1},
     {1, 0, -1},   {-1, 0, 1},  {-1, 0, -1}, {0, 1, 1},  {0, 1, -1}, {0, -1, 1},
     {0, -1, -1},  {1, 1, 1},   {1, 1, -1},  {1, -1, 1}, {-1, 1, 1}, {-1, -1, 1},
     {-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}},
    {0, 4, 5, 6, 1, 2, 3, 10, 11, 8, 9, 7, 14, 15, 12, 13, 16, 17, 18, 19, 20, 21, 22, 23, 24});