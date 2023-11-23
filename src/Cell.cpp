#include <cmath>
#include "Cell.hpp"
#include "Utils.hpp"
#define TAU 0.5f
#define R 8.214f
#define T 300f
#define CS std::sqrt(1/3)




void Cell::update(const float deltaTime, const Lattice& lattice){
    //momentum density, f, feq, omega
    rho = 0;
    for (int i = 0; i<9; i++)    //update rho
    {
        rho += f[i];
    };
    
    ux = ( f[1] + f[5] + f[8] - ( f[3] + f[6] + f[7] )) / rho;
    uy = ( f[2] + f[5] + f[6] - ( f[4] + f[7] + f[8] )) / rho;

}

void Cell::updateFeq(std::vector<float> feq, const float& ux, const float&  uy, const float& rho) {
    float uprod = ux * ux + uy * uy;
    for (int i=0; i<9; i++) {
        float temp = ux * D2Q9.velocities[i][0] + uy * D2Q9.velocities[i][1];
        f[i] = D2Q9.weights[i] * rho * (1 + temp / std::pow(CS, 2) + std::pow(temp, 2) / (2 * std::pow(CS, 4)) - uprod / (2 * std::pow(CS, 2)));
    }
}
