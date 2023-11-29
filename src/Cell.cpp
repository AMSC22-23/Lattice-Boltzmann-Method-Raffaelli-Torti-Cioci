#include "Cell.hpp"
#include "Lattice.hpp"
#include "Utils.cpp"
#include <cmath>

#define TAU 0.5f
#define R 8.214f
#define T 300f
#define CS std::sqrt(1 / 3)

// TODO stabilire parametri di inizializzazione
Cell::Cell(const std::vector<int> &boundary, const bool &obstacle)
    : boundary(boundary), obstacle(obstacle), f(D2Q9.velocity_directions, 0), newF(D2Q9.velocity_directions, 0),
      feq(D2Q9.velocity_directions, 0)
{
    // TODO Initialize the distribution function (e inizializzare anche il resto...)
    for (int i = 0; i < D2Q9.velocity_directions; i++)
    {
        /* non si sa da dove arrivi
        f.at(i) = D2Q9.weights.at(i) *
                  (1 + 3 * (D2Q9.velocities.at(i).at(0) * uX + D2Q9.velocities.at(i).at(1) * uY) / std::pow(CS, 2) +
                   9 * (D2Q9.velocities.at(i).at(0) * uX + D2Q9.velocities.at(i).at(1) * uY) *
                       (D2Q9.velocities.at(i).at(0) * uX + D2Q9.velocities.at(i).at(1) * uY) / (2 * std::pow(CS, 4)) -
                   3 * (uX * uX + uY * uY) / (2 * std::pow(CS, 2)));
                   */
    }
}

void Cell::update(const float deltaTime, Lattice &lattice, const std::vector<int> &cellPosition)
{
    rho = 0;
    for (int i = 0; i < D2Q9.velocity_directions; i++) // update rho
    {
        rho += f.at(i);
    }

    /*
    // TODO rendere 2 e 3 dimensionale
    uX = (f.at(1) + f.at(5) + f.at(8) - (f.at(3) + f.at(6) + f.at(7))) / rho;
    uY = (f.at(2) + f.at(5) + f.at(6) - (f.at(4) + f.at(7) + f.at(8))) / rho;
    */

    updateFeq();
    collision(deltaTime);                // TODO may return fstar
    streaming(f, lattice, cellPosition); // TODO be careful of fstar

    // Copy fnew to f
    for (int i = 0; i < D2Q9.velocity_directions; i++)
    {
        f.at(i) = newF.at(i);
    }
}

void Cell::updateFeq()
{
    // TODO riscrivere per u n-dimensionale
    float uProd = scalar_product(macroU, macroU);
    /*
    float uprod = uX * uX + uY * uY;
    for (int i = 0; i < D2Q9.velocity_directions; i++)
    {
        float temp = uX * D2Q9.velocities.at(i).at(0) + uY * D2Q9.velocities.at(i).at(1);
        feq.at(i) =
            D2Q9.weights.at(i) * rho *
            (1 + temp / std::pow(CS, 2) + std::pow(temp, 2) / (2 * std::pow(CS, 4)) - uprod / (2 * std::pow(CS, 2)));
    }
    */
}

void Cell::collision(const float deltaTime)
{
    for (int i = 0; i < D2Q9.velocity_directions; i++)
    {
        f.at(i) = f.at(i) * (1 - deltaTime / TAU) + feq.at(i) * deltaTime / TAU;
    }
}

// TODO read fstar again
void Cell::streaming(const std::vector<float> fstar, Lattice &lattice, const std::vector<int> &position)
{
    // consider one velocity at a time
    for (int i = 0; i < D2Q9.velocity_directions; i++)
    {
        // if there is no boundary in the direction of the velocity, we stream
        if (boundary.at(0) != D2Q9.velocities.at(i).at(0) && boundary.at(1) != D2Q9.velocities.at(i).at(1))
        {
            Cell newCell = lattice.getCellAtIndex(
                {position.at(0) + D2Q9.velocities.at(i).at(0), position.at(1) + D2Q9.velocities.at(i).at(1)});

            newCell.setFAtIndex(i, fstar.at(i));
        }
        else // otherwise we bounce back (potentially in both directions)
        {

            if (lattice.isLid() && position.at(0) == 1) //  we are going against the moving wall
            {
                const float lid_velocity = 0.1f; //  only x component of the velocity
                newF.at(D2Q9.opposite.at(i)) = fstar.at(i) - 2 * D2Q9.weights.at(i) * rho *
                                                                 (D2Q9.velocities.at(i).at(0) * lid_velocity) /
                                                                 std::pow(CS, 2);
            }
            else
            {
                // we have to bounce BACK in the direction of the velocity we are considerint
                // we don't have to bounce FORWARD (we are in a no-slip condition between the fluid and the resting
                // wall)

                newF.at(D2Q9.opposite.at(i)) = fstar.at(i);
            }
        }
    }
}

void Cell::setFAtIndex(const int index, const float value)
{
    newF.at(index) = value;
}

bool Cell::isObstacle()
{
    return obstacle;
}
