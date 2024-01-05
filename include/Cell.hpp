#ifndef CELL_HPP
#define CELL_HPP

#include "Structure.hpp"
#include <vector>
class Lattice;

class Cell
{
  public:
    Cell(const Structure &structure, const std::vector<int> &_boundary, const bool _obstacle,
         const std::vector<float> &_f);

    // update
    void updateMacro(const Structure &structure);
    void equilibriumCollision(const Structure &structure, const float omP, const float halfOmpOmmSum,
                              const float halfOmpOmmSub);
    void streaming(Lattice &lattice, const std::vector<int> &position);
    void setInlets(const Structure &structure, const float &uLid, const std::vector<int> &position, const int &problemType, const int &dim); //this works
    void zouHe(const int problemType, const std::vector<float> &closeU, const float &uLidNow, const std::vector<int> &position, const float &closeRho);

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
};
/*

*/

#endif // CELL_HPP