#ifndef CELL_HPP
#define CELL_HPP
#include "Lattice.hpp"


class Cell {
    public:
        Cell();
        ~Cell();
        void update(const float deltaTime, const Lattice& lattice);

    private:
       void updateFeq(std::vector<float> feq, const float& ux, const float&  uy, const float& rho);
       
       float  f[9]; // Distribution function
       float feq[9]; // Equilibrium Distribution function
       float Omega [9]; // Collision operator 
       bool obstacle; // Is this cell an obstacle?
       float ux, uy; // Macroscopic velocity
       float momx, momy; // Momentum density
       float  rho; // Macroscopic density
};

        
#endif // CELL_HPP