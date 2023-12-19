#include "Cell.hpp"
#include "Lattice.hpp"
#include "Utils.cpp"

Cell::Cell(const Structure &structure, const std::vector<int> &_boundary, const bool &_obstacle,
           const std::vector<float> &_f)
{
    boundary = _boundary;
    obstacle = _obstacle;
    f = _f;
    feq = std::vector<float>(structure.velocity_directions, 0.0);

    macroU = std::vector<float>(structure.dimensions, 0.0);
    newF = std::vector<float>(structure.velocity_directions, 0.0);
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

void Cell::updateFeq(const Structure &structure)
{
    // update equilibrium distribution
    const float temp1 = 1.5 * scalar_product_parallel<float>({macroU, macroU});
    for (int i = 0; i < structure.velocity_directions; i++)
    {
        const float temp2 = 3.0 * scalar_product_parallel<float>({structure.velocities_by_direction.at(i), macroU});
        feq.at(i) = structure.weights.at(i) * rho * (1.0 + temp2 + 0.5 * temp2 * temp2 - temp1);
    }
}

void Cell::collisionStreaming(Lattice &lattice, const std::vector<int> &position, const float &omP, const float &omM)
{
    const Structure &structure = lattice.getStructure();

    // collide and stream for index 0
    newF.at(0) = (1.0 - omP) * f.at(0) + omP * feq.at(0);
    f.at(0) = newF.at(0);

    // collide for other indices
    for (int i = 1; i < structure.velocity_directions; i++)
    {
        newF.at(i) = (1.0 - 0.5 * (omP + omM)) * f.at(i) - 0.5 * (omP - omM) * f.at(structure.opposite.at(i)) +
                     0.5 * (omP + omM) * feq.at(i) + 0.5 * (omP - omM) * feq.at(structure.opposite.at(i));
    }

    // stream for other indices
    for (int i = 1; i < structure.velocity_directions; i++)
    {
        // if there is no boundary in the direction of the velocity, we stream
        if (boundary.at(0) != structure.velocities_by_direction.at(i).at(0) &&
            boundary.at(1) != structure.velocities_by_direction.at(i).at(1))
        {
            Cell &newCell =
                lattice.getCellAtIndices({position.at(0) + structure.velocities_by_direction_int.at(i).at(0),
                                          position.at(1) + structure.velocities_by_direction_int.at(i).at(1)});

            newCell.setFAtIndex(i, newF.at(i));
        }
    }
}

void Cell::setInlets(const Structure &structure, const float &uLid, const int &problemType)
{
    // if i'm at any boundary set macroU to 0
    if (boundary.at(0) == 1 || boundary.at(1) == 1 || boundary.at(0) == -1 || boundary.at(1) == -1)
    {
        macroU.at(0) = 0;
        macroU.at(1) = 0;
    }
    // if i'm at top wall set macroU.x to uLid
    if (boundary.at(1) == -1)
    {
        macroU.at(0) = uLid;
    }
}

void Cell::zouHe()
{
    if (boundary.at(0) == 0 && boundary.at(1) == -1) //  top wall (lid), not corner
    {
        rho = f.at(0) + f.at(1) + f.at(3) + 2 * (f.at(2) + f.at(5) + f.at(6)) / (1 + macroU.at(1));
        f.at(4) = f.at(2) - 2.0 / 3.0 * rho * macroU.at(1);
        f.at(7) = f.at(5) + 0.5 * (f.at(1) - f.at(3)) - 0.5 * rho * macroU.at(0) - 1.0 / 6.0 * rho * macroU.at(1);
        f.at(8) = f.at(6) - 0.5 * (f.at(1) - f.at(3)) + 0.5 * rho * macroU.at(0) - 1.0 / 6.0 * rho * macroU.at(1);
    }
    else if (boundary.at(0) == 1 && boundary.at(1) == 0) // right wall
    {
        rho = f.at(0) + f.at(2) + f.at(4) + 2 * (f.at(1) + f.at(5) + f.at(8)) / (1 + macroU.at(0));
        f.at(3) = f.at(1) - 2.0 / 3.0 * rho * macroU.at(0);
        f.at(7) = f.at(5) + 0.5 * (f.at(2) - f.at(4)) - 1.0 / 6.0 * rho * macroU.at(0) - 0.5 * rho * macroU.at(1);
        f.at(6) = f.at(8) - 0.5 * (f.at(2) - f.at(4)) + 1.0 / 6.0 * rho * macroU.at(0) - 0.5 * rho * macroU.at(1);
    }
    else if (boundary.at(0) == 0 && boundary.at(1) == 1) // bottom wall
    {
        rho = f.at(0) + f.at(1) + f.at(3) + 2 * (f.at(4) + f.at(7) + f.at(8)) / (1 - macroU.at(1));
        f.at(2) = f.at(4) + 2.0 / 3.0 * rho * macroU.at(1);
        f.at(5) = f.at(7) + 0.5 * (f.at(3) - f.at(1)) + 0.5 * rho * macroU.at(0) + 1.0 / 6.0 * rho * macroU.at(1);
        f.at(6) = f.at(8) - 0.5 * (f.at(3) - f.at(1)) - 0.5 * rho * macroU.at(0) + 1.0 / 6.0 * rho * macroU.at(1);
    }
    else if (boundary.at(0) == -1 && boundary.at(1) == 0) // left wall
    {
        rho = f.at(0) + f.at(2) + f.at(4) + 2 * (f.at(3) + f.at(6) + f.at(7)) / (1 - macroU.at(0));
        f.at(1) = f.at(3) + 2.0 / 3.0 * rho * macroU.at(0);
        f.at(5) = f.at(7) + 0.5 * (f.at(4) - f.at(2)) + 1.0 / 6.0 * rho * macroU.at(0) + 0.5 * rho * macroU.at(1);
        f.at(8) = f.at(6) - 0.5 * (f.at(4) - f.at(2)) - 1.0 / 6.0 * rho * macroU.at(0) + 0.5 * rho * macroU.at(1);
    }
    else if (boundary.at(0) == 1 && boundary.at(1) == -1) // top right corner
    {
        f.at(3) = f.at(1) - 2.0 / 3.0 * rho * macroU.at(0);
        f.at(4) = f.at(2) - 2.0 / 3.0 * rho * macroU.at(1);
        f.at(7) = f.at(5) - 1.0 / 6.0 * rho * macroU.at(0) - 1.0 / 6.0 * rho * macroU.at(1);
        f.at(8) = 0;
        f.at(6) = 0;
        f.at(0) = rho - f.at(1) - f.at(2) - f.at(3) - f.at(4) - f.at(5) - f.at(7);
    }
    else if (boundary.at(0) == 1 && boundary.at(1) == 1) // bottom right corner
    {
        f.at(3) = f.at(1) - 2.0 / 3.0 * rho * macroU.at(0);
        f.at(2) = f.at(4) + 2.0 / 3.0 * rho * macroU.at(1);
        f.at(6) = f.at(8) - 1.0 / 6.0 * rho * macroU.at(0) + 1.0 / 6.0 * rho * macroU.at(1);
        f.at(5) = 0;
        f.at(7) = 0;
        f.at(0) = rho - f.at(1) - f.at(2) - f.at(3) - f.at(4) - f.at(6) - f.at(8);
    }
    else if (boundary.at(0) == -1 && boundary.at(1) == 1) // bottom left corner
    {
        f.at(1) = f.at(3) + 2.0 / 3.0 * rho * macroU.at(0);
        f.at(2) = f.at(4) + 2.0 / 3.0 * rho * macroU.at(1);
        f.at(5) = f.at(7) + 1.0 / 6.0 * rho * macroU.at(0) + 1.0 / 6.0 * rho * macroU.at(1);
        f.at(6) = 0;
        f.at(8) = 0;
        f.at(0) = rho - f.at(1) - f.at(2) - f.at(3) - f.at(4) - f.at(5) - f.at(7);
    }
    else if (boundary.at(0) == -1 && boundary.at(1) == -1) // top left corner
    {
        f.at(1) = f.at(3) + 2.0 / 3.0 * rho * macroU.at(0);
        f.at(4) = f.at(2) - 2.0 / 3.0 * rho * macroU.at(1);
        f.at(8) = f.at(6) - 1.0 / 6.0 * rho * macroU.at(0) + 1.0 / 6.0 * rho * macroU.at(1);
        f.at(5) = 0;
        f.at(7) = 0;
        f.at(0) = rho - f.at(1) - f.at(2) - f.at(3) - f.at(4) - f.at(6) - f.at(8);
    }
}

void Cell::setFAtIndex(const int index, const float &value)
{
    f.at(index) = value;
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
