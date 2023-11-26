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
    rho = 0;
    for (int i = 0; i < D2Q9.velocity_directions; i++) // update rho
    {
        rho += f.at(i);
    }

    uX = (f.at(1) + f.at(5) + f.at(8) - (f.at(3) + f.at(6) + f.at(7))) / rho;
    uY = (f.at(2) + f.at(5) + f.at(6) - (f.at(4) + f.at(7) + f.at(8))) / rho;

    updateFeq(feq, uX, uY, rho);
    collision(feq, f, deltaTime);
    streaming(f, lattice);

    // Copy fnew to f
    for (int i = 0; i < D2Q9.velocity_directions; i++)
    {
        f.at(i) = newF.at(i);
    }
}

void Cell::updateFeq(std::vector<float> &feq, const float &uX, const float &uY, const float &rho)
{
    float uprod = uX * uX + uY * uY;
    for (int i = 0; i < D2Q9.velocity_directions; i++)
    {
        float temp = uX * D2Q9.velocities.at(i).at(0) + uY * D2Q9.velocities.at(i).at(1);
        feq.at(i) =
            D2Q9.weights.at(i) * rho *
            (1 + temp / std::pow(CS, 2) + std::pow(temp, 2) / (2 * std::pow(CS, 4)) - uprod / (2 * std::pow(CS, 2)));
    }
}

void Cell::collision(const std::vector<float> &feq, std::vector<float> &f, const float dt)
{
    for (int i = 0; i < D2Q9.velocity_directions; i++)
    {
        f.at(i) = f.at(i) * (1 - dt / TAU) + feq.at(i) * dt / TAU;
    }
}

void Cell::streaming(const std::vector<float> fstar, Lattice &lattice)
{
    // consider one velocity at a time
    for (int i = 0; i < D2Q9.velocity_directions; i++)
    {
        // if there is no boundary in the direction of the velocity, we stream
        if (boundary.at(0) != D2Q9.velocities.at(i).at(0) && boundary.at(1) != D2Q9.velocities.at(i).at(1))
        {
            Cell &newCell = lattice.getCellAtIndex(
                {position.at(0) + D2Q9.velocities.at(i).at(0), position.at(1) + D2Q9.velocities.at(i).at(1)});

            newCell.setFAtIndex(i, fstar.at(i));
        }
        else // otherwise we bounce back (potentially in both directions)
        {
            

            
            if (lattice.isLid() && position.at(0) == 1 ) //  we are going against the moving wall
            {
                const float lid_velocity = 0.1f; //  only x component of the velocity
                newF.at(D2Q9.opposite.at(i)) = fstar.at(i) - 2 * D2Q9.weights.at(i) * rho * (D2Q9.velocities.at(i).at(0) * lid_velocity) / std::pow(CS, 2);
            }
            else{
                //  we have to bounce BACK in the direction of the velocity we are considerint
                // we don't have to bounce FORWARD (we are in a no-slip condition between the fluid and the resting wall)
                
                newF.at(D2Q9.opposite.at(i)) = fstar.at(i);

            }
        
   
        }
    }
}

void Cell::setFAtIndex(const int index, const float value)
{
    newF.at(index) = value;
}

void Cell::setObstacle()
{
    obstacle = true;
}

bool Cell::isObstacle()
{
    return obstacle;
}
