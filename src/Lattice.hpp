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
    Lattice(std::string filename);
    void update(const float deltaTime, std::ofstream &file);
    const Cell &getCellAtIndices(std::vector<int> index) const;
    Cell &getMutableCellAtIndices(std::vector<int> index);
    const std::vector<int> getShape();
    bool isLid();
    const Structure &getStructure() const;

  private:
    NDimensionalMatrix<Cell> cells;
    bool lid = false;
    Structure structure;
    int timeInstant = 0;
};

#endif // LATTICE_HPP