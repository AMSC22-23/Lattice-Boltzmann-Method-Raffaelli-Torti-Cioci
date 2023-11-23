#ifndef CELL_HPP
#define CELL_HPP
#include "Lattice.hpp"


class Cell {
    public:
        Cell();
        ~Cell();
        void update(const float deltaTime, const Lattice& lattice);

    private:
       void updateFeq(float feq[], const float& ux, const float&  uy, const float& rho);
       void collision(const std::vector<float> feq, std::vector<float> f, const float dt  );
       void streaming();
       float f[9]; // Distribution function
       float feq[9]; // Equilibrium Distribution function
       float Omega [9]; // Collision operator 
       bool obstacle; // Is this cell an obstacle?
       float uX, uY; // Macroscopic velocity
       float momX, momY; // Momentum density
       float  rho; // Macroscopic density
       const int posX, posY; //cell position in the lattice

};

        
#endif // CELL_HPP