#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "Cell.hpp"
#include "Utils.cpp"
#include "Structure.hpp"
#include <fstream>
#include <iterator>
#include <string>

class Lattice
{
  public:
    Lattice(std::ifstream &file_in, const int plotSteps);
    void simulate(std::ofstream &velocity_out, std::ofstream &lift_drag_out);
    void simulateGpu(std::ofstream &velocity_out, std::ofstream &lift_drag_out);
    Cell &getCellAtIndices(const std::vector<int> &indices);
    Cell &getCellAtIndices(const int x, const int y);
    Cell &getCellAtIndices(const int x, const int y, const int z);
    Cell &getCellAtIndices(const int *indices);
    const std::vector<float> &getCloseU(const std::vector<int> &indices);
    float getCloseRho(const std::vector<int> &indices);
    const Cell &getCellAtIndex(const int index) const;
    const std::vector<int> getShape() const;
    const Structure &getStructure() const;
    int getProblemType() const;

  private:
    NDimensionalMatrix<Cell> cells;
    int problemType;
    Structure structure;
    int timeInstant = 0;
    float sigma;
    float uLid;
    float omP;
    float omM;
    int maxIt;
    int plotSteps;
};

#endif // LATTICE_HPP