#include "Cell.hpp"
#include "Lattice.hpp"
#include "Utils.cpp"
#include <cmath>
#include <numeric>

#define TAU 0.5f
#define CS std::sqrt(1 / 3)

Cell::Cell(const Structure &structure, const std::vector<int> &boundary, const bool &obstacle,
           const std::vector<float> &macroUInput, const float &reynoldsNumber, const float &length, const float &mu)
    : boundary(boundary), obstacle(obstacle), f(structure.velocity_directions, 0),
      newF(structure.velocity_directions, 0), feq(structure.velocity_directions, 0), macroU(macroUInput)

// to give in input macroU and rho, coerent with Reynolds number
{
    rho = reynoldsNumber * mu / (length * vector_modulus_parallel(macroU));

    for (int i = 0; i < structure.velocity_directions; i++)
    {
        updateFeq(structure);
        f.at(i) = feq.at(i);
    }

    marcoRhoU = scalar_vector_product_parallel(rho, macroU);
}

void Cell::update(const float deltaTime, Lattice &lattice, const std::vector<int> &cellPosition)
{
    const Structure &structure = lattice.getStructure();
    rho = 0;
    for (int i = 0; i < structure.velocity_directions; i++) // update rho
    {
        rho += f.at(i);
    }

    updateFeq(structure);
    collision(structure, deltaTime);
    streaming(lattice, cellPosition);

    // Copy fnew to f
    for (int i = 0; i < structure.velocity_directions; i++)
    {
        f.at(i) = newF.at(i);
    }

    // Update macroscopic velocity
    for (int i = 0; i < structure.dimensions; i++)
    {
        macroU.at(i) = 0;
        /*
        for (int j = 0; j < D2Q9.velocity_directions; j++)
        {
            macroU.at(i) += D2Q9.weights.at(j) * D2Q9.velocities.at(j).at(i) * f.at(j);
        }
        */
        // !
        std::vector<float> weights;
        std::vector<float> velocities;
        for (int j = 0; j < structure.velocity_directions; j++)
        {
            weights.push_back(structure.weights.at(j));
            velocities.push_back(structure.velocities.at(j).at(i));
        }
        macroU.at(i) += scalar_product_parallel<float>({weights, velocities, f});
        // !
    }
}

void Cell::updateFeq(const Structure &structure)
{

    auto uProd = std::inner_product(macroU.begin(), macroU.end(), macroU.begin(), 0.0f);

    for (int i = 0; i < structure.velocity_directions; i++)
    {
        // ! auto temp = std::inner_product(macroU.begin(), macroU.end(), D2Q9.velocities.at(i).begin(), 0.0f);
        std::vector<float> velocities;
        for (int j = 0; j < structure.dimensions; j++)
        {
            velocities.push_back(structure.velocities.at(i).at(j));
        }
        auto temp = scalar_product_parallel<float>({macroU, velocities});
        // !
        feq.at(i) =
            structure.weights.at(i) * rho *
            (1 + temp / std::pow(CS, 2) + std::pow(temp, 2) / (2 * std::pow(CS, 4)) - uProd / (2 * std::pow(CS, 2)));
    }
}

void Cell::collision(const Structure &structure, const float deltaTime)
{
    for (int i = 0; i < structure.velocity_directions; i++)
    {
        f.at(i) = f.at(i) * (1 - deltaTime / TAU) + feq.at(i) * deltaTime / TAU; // in the algorithm it is called fstar
    }
}

void Cell::streaming(Lattice &lattice, const std::vector<int> &position)
{
    const Structure &structure = lattice.getStructure();
    // consider one velocity at a time
    for (int i = 0; i < structure.velocity_directions; i++)
    {
        // if there is no boundary in the direction of the velocity, we stream
        if (boundary.at(0) != structure.velocities.at(i).at(0) && boundary.at(1) != structure.velocities.at(i).at(1))
        {
            Cell newCell = lattice.getCellAtIndices(
                {position.at(0) + structure.velocities.at(i).at(0), position.at(1) + structure.velocities.at(i).at(1)});

            newCell.setNewFAtIndex(i, f.at(i));
        }
        else // otherwise we bounce back (potentially in both directions)
        {

            if (lattice.isLid() && position.at(0) == 1) //  we are going against the moving wall
            {
                const float lid_velocity = 0.1f; //  only x component of the velocity
                newF.at(structure.opposite.at(i)) = f.at(i) - 2 * structure.weights.at(i) * rho *
                                                                  (structure.velocities.at(i).at(0) * lid_velocity) /
                                                                  std::pow(CS, 2);
            }
            else
            {
                // we have to bounce BACK in the direction of the velocity we are considerint
                // we don't have to bounce FORWARD (we are in a no-slip condition between the fluid and the resting
                // wall)

                newF.at(structure.opposite.at(i)) = f.at(i);
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

bool Cell::isObstacle() const
{
    return obstacle;
}

Cell &Cell::operator=(const Cell &other)
{
    f = other.f;
    newF = other.newF;
    feq = other.feq;
    macroU = other.macroU;
    marcoRhoU = other.marcoRhoU;
    rho = other.rho;
    return *this;
}