#include "Cell.hpp"
#include "Lattice.hpp"
#include "Utils.cpp"
#include <algorithm>
#include <cmath>
#include <numeric>

static constexpr float TAU = 0.5;
static constexpr float CS = 0.57735;

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

    updateFeq(structure, isLidCell);
    collision_fast(structure, deltaTime);
    streaming(lattice, cellPosition);
}

void Cell::update2(const Structure &structure)
{
    updateF(structure);
    updateMacro(structure);
}

void Cell::updateFeq(const Structure &structure, const bool isLidCell)
{
    if (isLidCell) // apply zou-he boundary condition
    {
        // TODO as of right now it will be fixed
        feq.at(0) = feq.at(2) = feq.at(4) = 1;
        feq.at(1) = feq.at(5) = feq.at(8) = 4;
        feq.at(3) = feq.at(6) = feq.at(7) = 0.25;

        return;
    }

    // square of modulus of macroU
    float uProd = scalar_product_parallel<float>({macroU, macroU});

    for (int i = 0; i < structure.velocity_directions; i++)
    {
        float temp = scalar_product_parallel<float>({macroU, structure.velocities_by_direction.at(i)});
        float temp2 = std::pow(CS, 2);
        feq.at(i) = structure.weights.at(i) * rho *
                    (1.0 + temp / temp2 + std::pow(temp, 2) / (2.0 * temp2 * temp2) - uProd / (2.0 * temp2));
        feq.at(i) = std::max(feq.at(i), 0.0f);
    }
}

void Cell::updateMacro(const Structure &structure)
{
    // Update macroscopic velocity
    for (int i = 0; i < structure.dimensions; i++)
    {
        macroU.at(i) = 0;
        macroU.at(i) += scalar_product_parallel<float>({structure.weights, structure.velocities_by_dimension.at(i), f});
    }
    // update macroscopic density (rho)
    rho = std::accumulate(f.begin(), f.end(), 0.0f);
}

void Cell::updateF(const Structure &structure)
{
    // Copy fnew to f
    for (int i = 0; i < structure.velocity_directions; i++)
    {
        f.at(i) = newF.at(i);
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

void Cell::streaming(Lattice &lattice, const std::vector<int> &position)
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
        // otherwise we bounce back
        else
        {
            newF.at(structure.opposite.at(i)) = f.at(i);
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