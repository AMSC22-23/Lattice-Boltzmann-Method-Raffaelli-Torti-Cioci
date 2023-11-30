#include "Cell.hpp"
#include "Lattice.hpp"
#include "Utils.cpp"
#include <cmath>
#include <numeric>

#define TAU 0.5f
#define R 8.214f
#define T 300f
#define CS std::sqrt(1 / 3)

Cell::Cell(const std::vector<int> &boundary, const bool &obstacle, const std::vector<float> &macroUInput, const float &rhoInput)
    : boundary(boundary), obstacle(obstacle), f(D2Q9.velocity_directions, 0), newF(D2Q9.velocity_directions, 0),
      feq(D2Q9.velocity_directions, 0), rho(rhoInput), macroU(macroUInput) 

// to give in input macroU and rho, coerent with Reynolds number
{
    for (int i = 0; i < D2Q9.velocity_directions; i++)
    {
        updateFeq();
        //inizializzare/modificare: feq
        f.at(i) = feq.at(i);
        
    }
    
    for (int i = 0; i < D2Q9.dimensions; i++) //updating macroRhoU
    {
        marcoRhoU.at(i) = rho * macroU.at(i);
    }
}

void Cell::update(const float deltaTime, Lattice &lattice, const std::vector<int> &cellPosition)
{
    rho = 0;
    for (int i = 0; i < D2Q9.velocity_directions; i++) // update rho
    {
        rho += f.at(i);
    }

    

    updateFeq();
    collision(deltaTime);               
    streaming(lattice, cellPosition); 


    // Copy fnew to f
    for (int i = 0; i < D2Q9.velocity_directions; i++)
    {
        f.at(i) = newF.at(i);
    }
}

void Cell::updateFeq()
{
 
    float uProd = std::inner_product(macroU.begin(), macroU.end(), macroU.begin(), 0.0f);
    
    for (int i = 0; i < D2Q9.velocity_directions; i++)
    { 
        float temp = std::inner_product(macroU.begin(), macroU.end(), D2Q9.velocities.at(i).begin(), 0.0f);
        feq.at(i) =
            D2Q9.weights.at(i) * rho *
            (1 + temp / std::pow(CS, 2) + std::pow(temp, 2) / (2 * std::pow(CS, 4)) - uProd / (2 * std::pow(CS, 2)));
    }
    
}

void Cell::collision(const float deltaTime)
{
    for (int i = 0; i < D2Q9.velocity_directions; i++)
    {
        f.at(i) = f.at(i) * (1 - deltaTime / TAU) + feq.at(i) * deltaTime / TAU; // in the algorithm it is called fstar
    }
}

void Cell::streaming(Lattice &lattice, const std::vector<int> &position)
{
    // consider one velocity at a time
    for (int i = 0; i < D2Q9.velocity_directions; i++)
    {
        // if there is no boundary in the direction of the velocity, we stream
        if (boundary.at(0) != D2Q9.velocities.at(i).at(0) && boundary.at(1) != D2Q9.velocities.at(i).at(1))
        {
            Cell newCell = lattice.getCellAtIndex(
                {position.at(0) + D2Q9.velocities.at(i).at(0), position.at(1) + D2Q9.velocities.at(i).at(1)});

            newCell.setNewFAtIndex(i, f.at(i));
        }
        else // otherwise we bounce back (potentially in both directions)
        {

            if (lattice.isLid() && position.at(0) == 1) //  we are going against the moving wall
            {
                const float lid_velocity = 0.1f; //  only x component of the velocity
                newF.at(D2Q9.opposite.at(i)) = f.at(i) - 2 * D2Q9.weights.at(i) * rho *
                                                                 (D2Q9.velocities.at(i).at(0) * lid_velocity) /
                                                                 std::pow(CS, 2);
            }
            else
            {
                // we have to bounce BACK in the direction of the velocity we are considerint
                // we don't have to bounce FORWARD (we are in a no-slip condition between the fluid and the resting
                // wall)

                newF.at(D2Q9.opposite.at(i)) = f.at(i);
            }
        }
    }
}

void Cell::setFAtIndex(const int index, const float value)
{
    f.at(index) = value;
}

void Cell::setNewFAtIndex(const int index, const float value)
{
    newF.at(index) = value;
}

bool Cell::isObstacle()
{
    return obstacle;
}
