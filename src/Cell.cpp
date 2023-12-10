#include "Cell.hpp"
#include "Lattice.hpp"
#include "Utils.cpp"
#include <cmath>
#include <numeric>

static constexpr float TAU = 0.5;
static constexpr float CS = 0.57735;

/*
* old constructor
Cell::Cell(const Structure &structure, const std::vector<int> &_boundary, const bool &_obstacle,
           const float &reynoldsNumber, const float &length, const float &mu)
    : boundary(_boundary), obstacle(_obstacle), f(structure.velocity_directions, 0),
      newF(structure.velocity_directions, 0), feq(structure.velocity_directions, 0), macroU(structure.dimensions, 0)
{
    float ulid = reynoldsNumber * mu / (rho * length);

    for (int i = 0; i < structure.velocity_directions; i++)
    {
        updateFeq(structure);
        f.at(i) = feq.at(i);
    }

    marcoRhoU = scalar_vector_product_parallel(rho, macroU);
}
*/

Cell::Cell(const Structure &structure, const std::vector<int> &_boundary, const bool &_obstacle,
           const std::vector<float> &_f)
    : f(_f), newF(_f), boundary(_boundary), obstacle(_obstacle), macroU(structure.dimensions, 0), feq(_f)
{
    updateMacro(structure);
}

void Cell::update1(const float deltaTime, Lattice &lattice, const std::vector<int> &cellPosition)
{
    const Structure &structure = lattice.getStructure();
    bool isLidCell = false;
    if (lattice.isLid() && cellPosition.at(1) == 0)
    {
        isLidCell = true;
    }

    updateFeq(structure);
    collision(structure, deltaTime);
    streaming(lattice, cellPosition, isLidCell);
}

void Cell::update2(const Structure &structure)
{
    updateF(structure);
    updateMacro(structure);
}

void Cell::updateF(const Structure &structure)
{
    // Copy fnew to f
    for (int i = 0; i < structure.velocity_directions; i++)
    {
        f.at(i) = newF.at(i);
    }
}

void Cell::updateMacro(const Structure &structure)
{
    // update macroscopic density (rho)
    rho = std::accumulate(f.begin(), f.end(), 0.0f);
    // Update macroscopic velocity
    for (int i = 0; i < structure.dimensions; i++)
    {
        // ! without weights
        macroU.at(i) = scalar_product_parallel<float>({structure.velocities_by_dimension.at(i), f}) / rho;
    }
}

void Cell::updateFeq(const Structure &structure)
{
    // square of modulus of macroU
    const float uProd = scalar_product_parallel<float>({macroU, macroU});

    for (int i = 0; i < structure.velocity_directions; i++)
    {
        float temp = scalar_product_parallel<float>({macroU, structure.velocities_by_direction.at(i)});
        float temp2 = std::pow(CS, 2);
        feq.at(i) = structure.weights.at(i) * rho *
                    (1.0 + temp / temp2 + std::pow(temp, 2) / (2.0 * temp2 * temp2) - uProd / (2.0 * temp2));
        feq.at(i) = std::max(feq.at(i), 0.0f);
    }
}

void Cell::collision(const Structure &structure, const float deltaTime)
{
    for (int i = 0; i < structure.velocity_directions; i++)
    {
        f.at(i) = f.at(i) * (1 - deltaTime / TAU) + feq.at(i) * deltaTime / TAU; // in the algorithm it is called fstar
    }
}

void Cell::collision_fast(const Structure &structure, const float deltaTime)
{
    for (int i = 0; i < structure.velocity_directions; i++)
    {
        f.at(i) = feq.at(i);
    }
}

void Cell::streaming(Lattice &lattice, const std::vector<int> &position, const bool &isLidCell)
{
    const Structure &structure = lattice.getStructure();
    // consider one velocity at a time
    for (int i = 0; i < structure.velocity_directions; i++)
    {
        // if there is no boundary in the direction of the velocity, we stream
        if (boundary.at(0) != structure.velocities_by_direction.at(i).at(0) &&
            boundary.at(1) != structure.velocities_by_direction.at(i).at(1))
        {
            Cell newCell =
                lattice.getCellAtIndices({position.at(0) + structure.velocities_by_direction_int.at(i).at(0),
                                          position.at(1) - structure.velocities_by_direction_int.at(i).at(1)});

            newCell.setNewFAtIndex(i, f.at(i));
        }
        else // otherwise we bounce back (potentially in both directions)
        {

            if (isLidCell) //  we are going against the moving wall
            {
                // ! set lid Velocity based on reynolds number...
                const float lid_velocity = 0.1f;
                rho = f.at(0) + f.at(1) + f.at(3) + 2 * (f.at(2) + f.at(5) + f.at(6)); // NEBB system
                newF.at(4) = f.at(2);
                newF.at(7) = f.at(5) + 0.5 * (f.at(1) - f.at(3) - lid_velocity * rho);
                newF.at(8) = f.at(6) - 0.5 * (f.at(1) - f.at(3) - lid_velocity * rho);
            }
            else
            {
                newF.at(structure.opposite.at(i)) = f.at(i);
            }
        }
    }
}

void Cell::setNewFAtIndex(const int index, const float value)
{
    newF.at(index) = value;
}

const std::vector<float> &Cell::getMacroU() const
{
    return macroU;
}

const float &Cell::getRho() const
{
    return rho;
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
    boundary = other.boundary;
    obstacle = other.obstacle;
    rho = other.rho;
    return *this;
}