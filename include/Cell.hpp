#ifndef CELL_HPP
#define CELL_HPP

#include "Structure.hpp"
#include <vector>
class Lattice;

class Cell
{
  public:
    Cell(const Structure &structure, const std::vector<int> &_boundary, const bool _obstacle,
         const std::vector<int> &position);

    // update
    void updateMacro(const Structure &structure);
    void equilibriumCollision(const Structure &structure, const float omP, const float halfOmpOmmSum,
                              const float halfOmpOmmSub);
    void initEq(const Structure &structure, const int problemType);
    void streaming(Lattice &lattice);
    void setInlets(Lattice &lattice, const float uLidNow);
    void zouHe(Lattice &lattice);
    void bounceBackObstacle();
    void dragAndLift(float &drag, float &lift);

    // getters and setters
    const float &getRho() const;
    const std::vector<float> &getMacroU() const;
    const std::vector<float> &getF() const;
    const std::vector<float> &getNewF() const;
    const std::vector<int> &getBoundary() const;
    bool isObstacle() const;
    void setFAtIndex(const int index, const float &value);

    // other
    Cell() = default;

  private:
    std::vector<float> f;    // Distribution  (length == Qx)
    std::vector<float> newF; // Distribution streamed from neighboring cells (length == Qx)

    std::vector<float> macroU; // Macroscopic velocity (length == Dx)
    float rho;                 // Macroscopic density

    std::vector<int> boundary; // boundary conditions (length == Dx)
    bool obstacle = {false};   // Is this cell an obstacle?
    std::vector<int> position; // position of the cell in the lattice
};

#endif // CELL_HPP