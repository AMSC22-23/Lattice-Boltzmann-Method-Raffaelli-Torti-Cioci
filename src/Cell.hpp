#ifndef CELL_HPP
#define CELL_HPP

#include "Structure.hpp"
#include <vector>
class Lattice;

class Cell
{
  public:
    Cell(const Structure &structure, const std::vector<int> &boundary, const bool &obstacle,
         const std::vector<float> &macroUInput, const float &reynoldsNumber, const float &length, const float &mu);
    void update(const float deltaTime, Lattice &lattice, const std::vector<int> &cellPosition);
    void setFAtIndex(const int index, const float value);
    void setNewFAtIndex(const int index, const float value);
    void updatePartTwo(const Structure &structure);
    const float &getRho() const;
    const std::vector<float> &getMacroU() const;
    bool isObstacle() const;
    Cell() = default;
    // copy operator
    Cell &operator=(const Cell &other);

  private:
    void updateFeq(const Structure &structure);
    void collision(const Structure &structure, const float deltaTime);
    void streaming(Lattice &lattice, const std::vector<int> &position);
    std::vector<float> f;            // Distribution function (length == Qx)
    std::vector<float> newF;         // Updated distribution function (length == Qx)
    std::vector<float> feq;          // Equilibrium Distribution function (length == Qx)
    std::vector<float> macroU;       // Macroscopic velocity (length == Dx)
    std::vector<float> marcoRhoU;    // Momentum density (rho * u) (length == Dx)
    std::vector<int> boundary; // boundary conditions (length == Dx)
    bool obstacle = {false};         // Is this cell an obstacle?
    float rho;                       // Macroscopic density
};
// boundary is an array of two elements; eache element can be 0, 1 or -1.
// if we have a boundary {-1, 1} it means that we can't go in the "opposite direction of the x axis" (not on left)
// and we can't go in the direction of the y axis (not up).
// sum up: the first element of the vector tells us which direction on x we can't follow in that cell and
// the second element tells us which direction on y we can't follow. If one of the element is 0 it means that
// we have no boundary in this direction

#endif // CELL_HPP