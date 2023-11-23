#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "Utils.hpp"
#include "Cell.hpp"

class Lattice {
    public:
        Lattice();
        ~Lattice();
        void update(const float deltaTime);
        Cell& getCellAtIndex(const std::vector<unsigned int>& indices);

    private:
        NDimensionalMatrix<Cell> cells;
};

#endif // LATTICE_HPP