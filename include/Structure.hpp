#ifndef STRUCTURE_HPP
#define STRUCTURE_HPP

#include <unordered_map>
#include <vector>

class Structure
{
  public:
    int dimensions;
    int velocity_directions;
    std::vector<float> weights;
    std::vector<std::vector<float>> velocities_by_direction;
    std::vector<std::vector<float>> velocities_by_dimension;
    std::vector<int> opposite;

    Structure &operator=(const Structure &other);

    Structure() = default;

    Structure(int dimensions, int velocity_directions, std::vector<float> weights,
              std::vector<std::vector<float>> velocities_by_direction,
              std::vector<std::vector<float>> velocities_by_dimension, std::vector<int> opposite);

    static Structure D2Q9;
    static Structure D3Q27;
};

#endif // STRUCTURE_HPP