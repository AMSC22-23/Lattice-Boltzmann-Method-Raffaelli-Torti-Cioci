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
    Lattice(const std::string &filename);
    void simulate(std::ofstream &file_out);
    Cell &getCellAtIndices(std::vector<int> indices);
    const std::vector<int> getShape() const;
    bool isLid() const;
    const Structure &getStructure() const;

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
};

#endif // LATTICE_HPP