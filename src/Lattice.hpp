#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "Cell.hpp"
#include "Utils.cpp"
#include <fstream>
#include <iterator>
#include <string>

class Lattice
{
  public:
    Lattice(std::string filename);
    ~Lattice();
    void update(const float deltaTime);
    Cell &getCellAtIndex(std::vector<int> index);
    std::vector<int> getShape();
    bool isLid();

  private:
    NDimensionalMatrix<Cell> cells;
    bool lid = false;
};

#endif // LATTICE_HPP