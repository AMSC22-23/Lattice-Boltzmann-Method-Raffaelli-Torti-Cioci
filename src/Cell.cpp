#include "Cell.hpp"
#include "Lattice.hpp"
#include "Utils.cpp"
#include <cmath>

#define TAU 0.5f
#define R 8.214f
#define T 300f
#define CS std::sqrt(1 / 3)

void Cell::update(const float deltaTime, Lattice &lattice)
{
    // TODO fibo, marti: capire omega
    rho = 0;
    for (int i = 0; i < 9; i++) // update rho
    {
        rho += f[i];
    };

    uX = (f[1] + f[5] + f[8] - (f[3] + f[6] + f[7])) / rho;
    uY = (f[2] + f[5] + f[6] - (f[4] + f[7] + f[8])) / rho;

    updateFeq(feq, uX, uY, rho);
    collision(feq, f, deltaTime);
    streaming(f, lattice);
}

void Cell::updateFeq(std::vector<float> &feq, const float &uX, const float &uY, const float &rho)
{
    float uprod = uX * uX + uY * uY;
    for (int i = 0; i < 9; i++)
    {
        float temp = uX * D2Q9.velocities[i][0] + uY * D2Q9.velocities[i][1];
        feq.at(i) =
            D2Q9.weights[i] * rho *
            (1 + temp / std::pow(CS, 2) + std::pow(temp, 2) / (2 * std::pow(CS, 4)) - uprod / (2 * std::pow(CS, 2)));
    }
}

void Cell::collision(const std::vector<float> &feq, std::vector<float> &f, const float dt)
{
    for (int i = 0; i < 9; i++)
    {
        f.at(i) = f.at(i) * (1 - dt / TAU) + feq.at(i) * dt / TAU;
    }
}

void Cell::streaming(const std::vector<float> fstar, Lattice &lattice)
{
    for (int i = 0; i < 9; i++)
    {
        Cell &newCell =
            lattice.getCellAtIndex({position.at(0) + D2Q9.velocities[i][0], position.at(1) + D2Q9.velocities[i][1]});
        newCell.setFAtIndex(i, fstar.at(i));
    }
}

void Cell::setFAtIndex(const unsigned int index, const float value)
{
    f.at(index) = value;
}

void Cell::setObstacle()
{
    obstacle = true;
}

bool Cell::isObstacle()
{
    return obstacle;
}