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
    
    uX = ( f[1] + f[5] + f[8] - ( f[3] + f[6] + f[7] )) / rho;
    uY = ( f[2] + f[5] + f[6] - ( f[4] + f[7] + f[8] )) / rho;
    
    updateFeq(feq, uX, uY, rho);
    collison(feq, f, dt);
    streaming(f);

}
// TODO : controllare la questione della & per feq: io voglio modificare tutti gli elementi del vettore
void Cell::updateFeq(float feq[], const float& uX, const float&  uY, const float& rho) { 
    float uprod = uX * uX + uY * uY;
    for (int i=0; i<9; i++) {
        float temp = uX * D2Q9.velocities[i][0] + uY * D2Q9.velocities[i][1];
        feq.[i] = D2Q9.weights[i] * rho * (1 + temp / std::pow(CS, 2) + std::pow(temp, 2) / (2 * std::pow(CS, 4)) - uprod / (2 * std::pow(CS, 2)));
    }
}

void collision(const float feq[], float f[], const float dt  )
{
    for (int i=0; i<9; i++) {
        f[i] = f[i] * (1 - dt / TAU) + feq[i] * dt / TAU;

    }
}

void streaming(const float fstar[])
{
    for (int i = 0; i < 9; i++)
    {
        Cell& newCell = getCellAtIndex(posX + velocities[i][0], posY + velocities[i][1]);
        newcell.f[i] = fstar[i];
    }
}

