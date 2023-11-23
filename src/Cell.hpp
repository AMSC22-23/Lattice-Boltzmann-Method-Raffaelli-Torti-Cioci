#ifndef CELL_HPP
#define CELL_HPP
#include "Lattice.hpp"


class Cell {
    public:
        Cell();
        ~Cell();
        void update(const float deltaTime, const Lattice& lattice);

    private:
       uint f[9]; // Distribution functions
       bool obstacle; // Is this cell an obstacle?
       float ux, uy; // Macroscopic velocity
       uint rho; // Macroscopic density
};

        
#endif // CELL_HPP