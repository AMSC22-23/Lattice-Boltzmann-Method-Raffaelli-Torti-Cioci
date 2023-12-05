#ifndef CELL_HPP
#define CELL_HPP

#include "Structure.hpp"
#include <vector>
class Lattice;

class Cell
{
  public:
    // constructor
    Cell(const Structure &structure, const std::vector<int> &_boundary, const bool &_obstacle,
         const std::vector<float> &_f);

    // collide and stream
    void update1(const float deltaTime, Lattice &lattice, const std::vector<int> &cellPosition);
    // update F and macro variables
    void update2(const Structure &structure);

    // getters and setters
    const float &getRho() const;
    const std::vector<float> &getMacroU() const;
    bool isObstacle() const;
    void setNewFAtIndex(const int index, const float value);

    // other
    Cell() = default;
    Cell &operator=(const Cell &other);

  private:
    void collision(const Structure &structure, const float deltaTime);
    void collision_fast(const Structure &structure, const float deltaTime);
    void streaming(Lattice &lattice, const std::vector<int> &position);
    void updateFeq(const Structure &structure, const bool isLidCell); // updates feq
    void updateMacro(const Structure &structure); // updates rho and macroU
    void updateF(const Structure &structure); // updates f using newF

    std::vector<float> f;    // Distribution  (length == Qx)
    std::vector<float> newF; // Distribution streamed from neighboring cells (length == Qx)
    std::vector<float> feq;  // Equilibrium distribution  (length == Qx)

    std::vector<float> macroU; // Macroscopic velocity (length == Dx)
    float rho;                 // Macroscopic density

    std::vector<int> boundary; // boundary conditions (length == Dx)
    bool obstacle = {false};   // Is this cell an obstacle?
};
// boundary is an array of two elements; eache element can be 0, 1 or -1.
// if we have a boundary {-1, 1} it means that we can't go in the "opposite direction of the x axis" (not on left)
// and we can't go in the direction of the y axis (not up).
// sum up: the first element of the vector tells us which direction on x we can't follow in that cell and
// the second element tells us which direction on y we can't follow. If one of the element is 0 it means that
// we have no boundary in this direction

#endif // CELL_HPP