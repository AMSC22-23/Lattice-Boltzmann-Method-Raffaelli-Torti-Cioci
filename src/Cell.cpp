#include "Cell.hpp"
#include "Lattice.hpp"
#include "Utils.cpp"
#include <cmath>
#include <numeric>

static constexpr float TAU = 1.1;
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

    updateFeq(structure, isLidCell);  // needs macroU. writes feq
    collision_fast(structure, deltaTime);  // needs feq, f. writes f
    streaming(lattice, cellPosition); // needs f. writes newF in neighbors
}

void Cell::update2(const std::vector<int> &shape, const std::vector<int> &position, const Structure &structure)
{
    // if cell is at boundary
    if (boundary.at(0) != 0 || boundary.at(1) != 0)
    {
        // ! need to force macroU at boundaries here.
        zouHe(shape, position); // needs newF. writes newF for boundary cells
    }
    else
    {
    }
    // ! move this up
    updateF(structure);     // needs newF. writes f
    updateMacro(structure); // needs f, writes rho, macroU
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
        macroU.at(i) = scalar_product_parallel<float>({structure.velocities_by_dimension.at(i), f}) / rho;
    }
}

void Cell::updateFeq(const Structure &structure, const bool &isLidCell)
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
            Cell &newCell =
                lattice.getCellAtIndices({position.at(0) + structure.velocities_by_direction_int.at(i).at(0),
                                          position.at(1) - structure.velocities_by_direction_int.at(i).at(1)});

            newCell.setNewFAtIndex(i, f.at(i));
        }
    }
}

void Cell::zouHe(const std::vector<int> &shape, const std::vector<int> &position)
{
    if (position.at(1) == 0 && position.at(0) != 0 && position.at(0) != shape.at(0) - 1) //  top wall (lid), not corner
    {
        // unknowns: 4, 7, 8
        const float lid_velocity = 0.2;
        rho = newF.at(0) + newF.at(1) + newF.at(3) + 2 * (newF.at(2) + newF.at(5) + newF.at(6)); // NEBB system
        const float temp = 0.5 * (newF.at(1) - newF.at(3) - lid_velocity * rho);
        newF.at(4) = newF.at(2);
        newF.at(7) = newF.at(5) + temp;
        newF.at(8) = newF.at(6) - temp;
    }
    else if (position.at(0) == shape.at(0) - 1 && position.at(1) == 0) // top right corner
    {
        // unknowns: 0, 6, 3, 7, 4, 8
        const float lid_velocity = 0.2;
        newF.at(3) = newF.at(1) - 2.0 / 3.0 * lid_velocity * rho;
        newF.at(7) = newF.at(5) - 1.0 / 6.0 * lid_velocity * rho;
        newF.at(4) = newF.at(2);
        newF.at(6) = -1.0 / 12.0 * lid_velocity * rho;
        newF.at(8) = 1.0 / 12.0 * lid_velocity * rho;
        // !
        newF.at(0) = rho - (newF.at(1) + newF.at(2) + newF.at(3) + newF.at(4) + newF.at(5) + newF.at(6) + newF.at(7) +
                            newF.at(8));
    }
    else if (position.at(1) == 0 && position.at(0) == 0) // top left corner
    {
        // unknowns: 0, 5, 1, 8, 4, 7
        const float lid_velocity = 0.2;
        newF.at(1) = newF.at(3) + 2.0 / 3.0 * lid_velocity * rho;
        newF.at(8) = newF.at(6) + 1.0 / 6.0 * lid_velocity * rho;
        newF.at(4) = newF.at(2);
        newF.at(5) = 1.0 / 12.0 * lid_velocity * rho;
        newF.at(7) = -1.0 / 12.0 * lid_velocity * rho;
        // !
        newF.at(0) = rho - (newF.at(1) + newF.at(2) + newF.at(3) + newF.at(4) + newF.at(5) + newF.at(6) + newF.at(7) +
                            newF.at(8));
    }
    else if (position.at(1) == shape.at(1) - 1 && position.at(0) == 0) // bottom left corner
    {
        // unknowns: 0, 1, 2, 5, 6, 8
        newF.at(1) = newF.at(3);
        newF.at(5) = newF.at(7);
        newF.at(2) = newF.at(4);
        newF.at(6) = 0;
        newF.at(8) = 0;
        // !
        newF.at(0) = rho - (newF.at(1) + newF.at(2) + newF.at(3) + newF.at(4) + newF.at(5) + newF.at(6) + newF.at(7) +
                            newF.at(8));
    }
    else if (position.at(1) == shape.at(1) - 1 && position.at(0) == shape.at(0) - 1) // bottom right corner
    {
        // unknowns: 0, 2, 3, 5, 7, 6
        newF.at(2) = newF.at(4);
        newF.at(6) = newF.at(8);
        newF.at(3) = newF.at(1);
        newF.at(7) = 0;
        newF.at(5) = 0;
        // !
        newF.at(0) = rho - (newF.at(1) + newF.at(2) + newF.at(3) + newF.at(4) + newF.at(5) + newF.at(6) + newF.at(7) +
                            newF.at(8));
    }
    else if (position.at(1) == shape.at(1) - 1 && position.at(0) != 0 &&
             position.at(0) != shape.at(0) - 1) // bottom wall
    {
        // unknowns: 2, 5, 6
        rho = newF.at(0) + newF.at(1) + newF.at(3) + 2 * (newF.at(4) + newF.at(7) + newF.at(8)); // NEBB system
        const float temp = 0.5 * (newF.at(1) - newF.at(3));
        newF.at(2) = newF.at(4);
        newF.at(5) = newF.at(7) + temp;
        newF.at(6) = newF.at(8) - temp;
    }
    else if (position.at(0) == shape.at(0) - 1 && position.at(1) != 0 &&
             position.at(1) != shape.at(1) - 1) // right wall
    {
        // unknowns: 3, 6, 7
        rho = newF.at(0) + newF.at(2) + newF.at(4) + 2 * (newF.at(5) + newF.at(1) + newF.at(8)); // NEBB system
        const float temp = 0.5 * (newF.at(2) - newF.at(4));
        newF.at(3) = newF.at(1);
        newF.at(7) = newF.at(5) + temp;
        newF.at(6) = newF.at(8) - temp;
    }
    else if (position.at(0) == 0 && position.at(1) != 0 && position.at(1) != shape.at(1) - 1) // left wall
    {
        // unknowns: 1, 5, 8
        rho = newF.at(0) + newF.at(2) + newF.at(4) + 2 * (newF.at(3) + newF.at(6) + newF.at(7)); // NEBB system
        const float temp = 0.5 * (newF.at(2) - newF.at(4));
        newF.at(1) = newF.at(3);
        newF.at(8) = newF.at(6) + temp;
        newF.at(5) = newF.at(7) - temp;
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