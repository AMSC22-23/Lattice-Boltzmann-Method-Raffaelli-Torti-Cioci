#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "Utils.hpp"
#include "Cell.hpp"
#include <string>

class Lattice {
    public:
        Lattice(std::string filename);
        ~Lattice();
        void update(const float deltaTime);
        Cell& getCellAtIndex(std::vector<unsigned int> index);

    private:
        NDimensionalMatrix<Cell> cells;
};

#endif // LATTICE_HPP