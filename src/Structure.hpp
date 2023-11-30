#pragma once

#include <unordered_map>
#include <vector>

class Structure
{
  public:
    int dimensions;
    int velocity_directions;
    std::unordered_map<int, float> weights;
    std::unordered_map<int, std::vector<int>> velocities;
    std::unordered_map<int, int> opposite;

    Structure &operator=(const Structure &other);

    Structure() = default;

    Structure(int dim, int vel_dir, const std::unordered_map<int, float> &w,
              const std::unordered_map<int, std::vector<int>> &v, const std::unordered_map<int, int> &opp);

    static Structure D2Q9;
    static Structure D3Q27;
};
