#include "Cell.hpp"
#include "Lattice.hpp"
#include "Utils.cpp"

float scalarProduct(const std::vector<float> &a, const std::vector<float> &b)
{
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0f);
}

float scalarProduct(const std::vector<float> &a, const std::vector<float> &b, const std::vector<float> &c)
{
    float res = 0;
    auto a_it = a.begin();
    auto b_it = b.begin();
    auto c_it = c.begin();

    for (; a_it != a.end() && b_it != b.end() && c_it != c.end(); ++a_it, ++b_it, ++c_it)
    {
        res += (*a_it) * (*b_it) * (*c_it);
    }

    return res;
}

Cell::Cell(const Structure &structure, const std::vector<int> &_boundary, const bool _obstacle,
           const std::vector<float> &_f)
{
    boundary = _boundary;
    obstacle = _obstacle;
    f = _f;

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
        macroU.at(i) = 0;
        for (int j = 0; j < structure.velocity_directions; j++)
        {
            macroU.at(i) += structure.velocities_by_dimension.at(i).at(j) * f.at(j);
        }
        macroU.at(i) /= rho;
    }
}

void Cell::equilibriumCollision(const Structure &structure, const float omP, const float halfOmpOmmSum,
                                const float halfOmpOmmSub)
{
    // equilibrium
    float feq[structure.velocity_directions] = {0.0};

    const float macroUSquareProd = scalarProduct(macroU, macroU);
    const float temp1 = 1.5 * macroUSquareProd;
    for (int i = 0; i < structure.velocity_directions; i++)
    {
        const float temp2 = 3.0 * scalarProduct(structure.velocities_by_direction.at(i), macroU);
        feq[i] = structure.weights.at(i) * rho * (1.0 + temp2 + 0.5 * temp2 * temp2 - temp1);
    }

    // collision for index 0
    newF.at(0) = (1.0 - omP) * f.at(0) + omP * feq[0];

    // collision for other indices
    for (int i = 1; i < structure.velocity_directions; i++)
    {
        newF.at(i) = (1.0 - halfOmpOmmSum) * f.at(i) - halfOmpOmmSub * f.at(structure.opposite.at(i)) +
                     halfOmpOmmSum * feq[i] + halfOmpOmmSub * feq[structure.opposite.at(i)];
    }
}

void Cell::streaming(Lattice &lattice, const std::vector<int> &position)
{
    const Structure &structure = lattice.getStructure();

    // stream for index 0
    f.at(0) = newF.at(0);
    bool stream;
    int new_position[structure.dimensions];

    // stream for other indices
    for (int i = 1; i < structure.velocity_directions; i++)
    {
        // if there is no boundary in the direction of the velocity...
        stream = true;
        for (int j = 0; j < structure.dimensions; j++)
        {
            if (boundary.at(j) == (int)structure.velocities_by_direction.at(i).at(j) && boundary.at(j) != 0)
            {
                stream = false;
                break;
            }
        }
        // ...stream
        if (stream)
        {
            for (int j = 0; j < structure.dimensions; j++)
            {
                new_position[j] = position.at(j) + (int)structure.velocities_by_direction.at(i).at(j);
            }
            lattice.getCellAtIndices(new_position).setFAtIndex(i, newF.at(i));
        }
    }
}

void Cell::setInlets(const Structure &structure, const float uLid, const int problemType)
{
    if (problemType == 1)
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
}

void Cell::zouHe()
{
    if (boundary.at(0) == 0 && boundary.at(1) == -1) // top wall
    {
        rho = (f.at(0) + f.at(1) + f.at(3) + 2.0 * (f.at(2) + f.at(5) + f.at(6))) / (1.0 + macroU.at(1));
        f.at(4) = f.at(2) - 2.0 / 3.0 * rho * macroU.at(1);
        f.at(7) = f.at(5) + 0.5 * (f.at(1) - f.at(3)) - 0.5 * rho * macroU.at(0) - 1.0 / 6.0 * rho * macroU.at(1);
        f.at(8) = f.at(6) - 0.5 * (f.at(1) - f.at(3)) + 0.5 * rho * macroU.at(0) - 1.0 / 6.0 * rho * macroU.at(1);
    }
    else if (boundary.at(0) == 1 && boundary.at(1) == 0) // right wall
    {
        rho = (f.at(0) + f.at(2) + f.at(4) + 2.0 * (f.at(1) + f.at(5) + f.at(8))) / (1.0 + macroU.at(0));
        f.at(3) = f.at(1) - 2.0 / 3.0 * rho * macroU.at(0);
        f.at(6) = f.at(8) - 0.5 * (f.at(2) - f.at(4)) - 1.0 / 6.0 * rho * macroU.at(0) + 0.5 * rho * macroU.at(1);
        f.at(7) = f.at(5) + 0.5 * (f.at(2) - f.at(4)) - 1.0 / 6.0 * rho * macroU.at(0) - 0.5 * rho * macroU.at(1);
    }
    else if (boundary.at(0) == 0 && boundary.at(1) == 1) // bottom wall
    {
        rho = (f.at(0) + f.at(1) + f.at(3) + 2.0 * (f.at(4) + f.at(7) + f.at(8))) / (1.0 - macroU.at(1));
        f.at(2) = f.at(4) + 2.0 / 3.0 * rho * macroU.at(1);
        f.at(5) = f.at(7) - 0.5 * (f.at(1) - f.at(3)) + 0.5 * rho * macroU.at(0) + 1.0 / 6.0 * rho * macroU.at(1);
        f.at(6) = f.at(8) + 0.5 * (f.at(1) - f.at(3)) - 0.5 * rho * macroU.at(0) + 1.0 / 6.0 * rho * macroU.at(1);
    }
    else if (boundary.at(0) == -1 && boundary.at(1) == 0) // left wall
    {
        rho = (f.at(0) + f.at(2) + f.at(4) + 2.0 * (f.at(3) + f.at(6) + f.at(7))) / (1.0 - macroU.at(0));
        f.at(1) = f.at(3) - 2.0 / 3.0 * rho * macroU.at(0);
        f.at(5) = f.at(7) - 0.5 * (f.at(2) - f.at(4)) + 1.0 / 6.0 * rho * macroU.at(0) + 0.5 * rho * macroU.at(1);
        f.at(8) = f.at(6) + 0.5 * (f.at(2) - f.at(4)) + 1.0 / 6.0 * rho * macroU.at(0) - 0.5 * rho * macroU.at(1);
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
        f.at(6) = f.at(8) + 1.0 / 6.0 * rho * macroU.at(1) - 1.0 / 6.0 * rho * macroU.at(0);
        f.at(7) = 0;
        f.at(5) = 0;
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
        f.at(7) = 0;
        f.at(5) = 0;
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

const std::vector<float> &Cell::getF() const
{
    return f;
}

const std::vector<float> &Cell::getNewF() const
{
    return newF;
}

const std::vector<int> &Cell::getBoundary() const
{
    return boundary;
}