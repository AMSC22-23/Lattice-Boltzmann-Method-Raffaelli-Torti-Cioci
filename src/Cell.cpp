#include "Cell.hpp"
#include "Lattice.hpp"
#include "Utils.cpp"

#include <iostream>

auto scalarProduct(const std::vector<float> &a, const std::vector<float> &b) -> float
{
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0f);
}

auto scalarProduct(const std::vector<float> &a, const std::vector<float> &b, const std::vector<float> &c) -> float
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
           const std::vector<int> &_position)
    : boundary(_boundary), obstacle(_obstacle), position(_position)
{
    rho = 1;// @note why not initialize in the initializer list?
    macroU = std::vector<float>(structure.dimensions, 0.0);//@note A strange way to resize a vector to a given size
    // macroU.resize(structure.dimensions, 0.0); // @note A better way to resize a vector to a given size
    f = std::vector<float>(structure.velocity_directions, 0.0);
    newF = std::vector<float>(structure.velocity_directions, 0.0);
}

/// @brief sets rho to sum of F and macroU to weighted sum of F normalized by rho
/// @param structure
void Cell::updateMacro(const Structure &structure)
{
    if (obstacle)
        return;

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

/// @brief collides F and stores the result in newF
/// @param structure
/// @param omP
/// @param halfOmpOmmSum
/// @param halfOmpOmmSub
void Cell::equilibriumCollision(const Structure &structure, const float omP, const float halfOmpOmmSum,
                                const float halfOmpOmmSub)
{

    if (obstacle)
        return;

    // equilibrium
    float feq[structure.velocity_directions] = {0.0}; // @note This is a C99 feature, not C++11
    // And also an error: a variable size array cannot be initialized
    // @note In C++11, you can use std::vector<float> feq(tructure.velocity_directions.0);
    const float macroUSquareProd = scalarProduct(macroU, macroU);
    const float temp1 = 1.5 * macroUSquareProd;
    for (int i = 0; i < structure.velocity_directions; i++)
    {
        const float temp2 = 3.0 * scalarProduct(structure.velocities_by_direction.at(i), macroU);//@note why float?
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

void Cell::initEq(const Structure &structure, const int problemType)
{
    if (obstacle && problemType == 2)
    {
        macroU.at(0) = std::numeric_limits<double>::quiet_NaN();// Better use optional to indicate missing values
        macroU.at(1) = std::numeric_limits<double>::quiet_NaN();
        return;
    }

    // equilibrium
    const float macroUSquareProd = scalarProduct(macroU, macroU);
    const float temp1 = 1.5 * macroUSquareProd;
    for (int i = 0; i < structure.velocity_directions; i++)
    {
        const float temp2 = 3.0 * scalarProduct(structure.velocities_by_direction.at(i), macroU);
        f[i] = structure.weights.at(i) * rho * (1.0 + temp2 + 0.5 * temp2 * temp2 - temp1);
    }
}

/// @brief streams newF to adjacent cells' F
/// @param lattice
void Cell::streaming(Lattice &lattice)
{
    if (obstacle)
        return;

    const Structure &structure = lattice.getStructure();

    // stream for index 0
    f.at(0) = newF.at(0);
    int new_position[structure.dimensions];

    // stream for other indices
    for (int i = 1; i < structure.velocity_directions; i++)
    {
        // calculate new position
        for (int j = 0; j < structure.dimensions; j++)
        {
            new_position[j] = position.at(j) + (int)structure.velocities_by_direction.at(i).at(j);
        }

        // stream if new position is inside the lattice and is not an obstacle
        if (new_position[0] >= 0 && new_position[0] < lattice.getShape().at(0) && new_position[1] >= 0 &&
            new_position[1] < lattice.getShape().at(1) && !lattice.getCellAtIndices(new_position).isObstacle())
        {
            lattice.getCellAtIndices(new_position).setFAtIndex(i, newF.at(i));// @note Implicity decay ann array to a pointer. Not nice
        }
    }
}

/// @brief sets macroU at walls depending on the problem type
/// @param lattice
/// @param uLidNow
void Cell::setInlets(Lattice &lattice, const float uLidNow)
{
    if (obstacle)
        return;

    // 2d only
    const int xLen = lattice.getShape().at(0);
    const int yLen = lattice.getShape().at(1);
    const int problemType = lattice.getProblemType();

    // if I'm at any wall set macroU to 0
    if (position.at(0) == 0 || position.at(1) == 0 || position.at(0) == xLen - 1 || position.at(1) == yLen - 1)
    {
        macroU.at(0) = 0;
        macroU.at(1) = 0;
    }

    switch (problemType)
    {
    case 1:
        // if i'm at top wall set macroU.x to uLidNow
        if (position.at(1) == 0)
            macroU.at(0) = uLidNow;

        break;
    // ! Checked by marti
    case 2:
        // if I'm at the left wall, calculate parabolic profile and set macroU.x to it
        if (position.at(0) == 0)
        {
            const float halfDim = static_cast<float>(yLen) / 2.0;
            const float temp = (static_cast<float>(position.at(1)) / halfDim) - 1.0;
            const float mul = 1.0 - temp * temp;
            macroU.at(0) = uLidNow * mul;
        }
        break;
    default:
        break;
    }
}

void Cell::zouHe(Lattice &lattice)
{
    if (obstacle)
        return;

    // 2d only
    const int xLen = lattice.getShape().at(0);
    const int yLen = lattice.getShape().at(1);

    const int problemType = lattice.getProblemType();

    // top wall
    if (position.at(0) != 0 && position.at(0) != xLen - 1 && position.at(1) == 0)
    {
        rho = (f.at(0) + f.at(1) + f.at(3) + 2.0 * (f.at(2) + f.at(5) + f.at(6))) / (1.0 + macroU.at(1));
        f.at(4) = f.at(2) - 2.0 / 3.0 * rho * macroU.at(1);
        f.at(7) = f.at(5) + 0.5 * (f.at(1) - f.at(3)) - 0.5 * rho * macroU.at(0) - 1.0 / 6.0 * rho * macroU.at(1);
        f.at(8) = f.at(6) - 0.5 * (f.at(1) - f.at(3)) + 0.5 * rho * macroU.at(0) - 1.0 / 6.0 * rho * macroU.at(1);
    }
    // right wall
    else if (position.at(0) == xLen - 1 && position.at(1) != 0 && position.at(1) != yLen - 1)
    {
        switch (problemType)
        {
        case 1:
            rho = (f.at(0) + f.at(2) + f.at(4) + 2.0 * (f.at(1) + f.at(5) + f.at(8))) / (1.0 + macroU.at(0));
            f.at(3) = f.at(1) - 2.0 / 3.0 * rho * macroU.at(0);
            f.at(6) = f.at(8) - 0.5 * (f.at(2) - f.at(4)) - 1.0 / 6.0 * rho * macroU.at(0) + 0.5 * rho * macroU.at(1);
            f.at(7) = f.at(5) + 0.5 * (f.at(2) - f.at(4)) - 1.0 / 6.0 * rho * macroU.at(0) - 0.5 * rho * macroU.at(1);
            break;
        case 2:
            rho = 1;
            macroU.at(0) = f.at(0) + f.at(2) + f.at(4) + 2.0 * (f.at(1) + f.at(5) + f.at(8)) - 1.0;
            f.at(3) = f.at(1) - 2.0 / 3.0 * macroU.at(0);
            f.at(6) = f.at(8) - 0.5 * (f.at(2) - f.at(4)) - 1.0 / 6.0 * macroU.at(0);
            f.at(7) = f.at(5) + 0.5 * (f.at(2) - f.at(4)) - 1.0 / 6.0 * macroU.at(0);
            break;
        default:
            break;
        }
    }
    // bottom wall
    else if (position.at(0) != 0 && position.at(0) != xLen - 1 && position.at(1) == yLen - 1)
    {
        rho = (f.at(0) + f.at(1) + f.at(3) + 2.0 * (f.at(4) + f.at(7) + f.at(8))) / (1.0 - macroU.at(1));
        f.at(2) = f.at(4) + 2.0 / 3.0 * rho * macroU.at(1);
        f.at(5) = f.at(7) - 0.5 * (f.at(1) - f.at(3)) + 0.5 * rho * macroU.at(0) + 1.0 / 6.0 * rho * macroU.at(1);
        f.at(6) = f.at(8) + 0.5 * (f.at(1) - f.at(3)) - 0.5 * rho * macroU.at(0) + 1.0 / 6.0 * rho * macroU.at(1);
    }
    // left wall
    else if (position.at(0) == 0 && position.at(1) != 0 && position.at(1) != yLen - 1)
    {
        rho = (f.at(0) + f.at(2) + f.at(4) + 2.0 * (f.at(3) + f.at(6) + f.at(7))) / (1.0 - macroU.at(0));
        f.at(1) = f.at(3) + 2.0 / 3.0 * rho * macroU.at(0);
        f.at(5) = f.at(7) - 0.5 * (f.at(2) - f.at(4)) + 1.0 / 6.0 * rho * macroU.at(0) + 0.5 * rho * macroU.at(1);
        f.at(8) = f.at(6) + 0.5 * (f.at(2) - f.at(4)) + 1.0 / 6.0 * rho * macroU.at(0) - 0.5 * rho * macroU.at(1);
    }

    // top right corner
    else if (position.at(0) == xLen - 1 && position.at(1) == 0)
    {
        rho = lattice.getCloseRho(position);

        f.at(3) = f.at(1) - 2.0 / 3.0 * rho * macroU.at(0);
        f.at(4) = f.at(2) - 2.0 / 3.0 * rho * macroU.at(1);
        f.at(7) = f.at(5) - 1.0 / 6.0 * rho * macroU.at(0) - 1.0 / 6.0 * rho * macroU.at(1);
        f.at(8) = 0;
        f.at(6) = 0;
        f.at(0) = rho - f.at(1) - f.at(2) - f.at(3) - f.at(4) - f.at(5) - f.at(7);
    }
    // bottom right corner
    else if (position.at(0) == xLen - 1 && position.at(1) == yLen - 1)
    {
        rho = lattice.getCloseRho(position);

        f.at(3) = f.at(1) - 2.0 / 3.0 * rho * macroU.at(0);
        f.at(2) = f.at(4) + 2.0 / 3.0 * rho * macroU.at(1);
        f.at(6) = f.at(8) + 1.0 / 6.0 * rho * macroU.at(1) - 1.0 / 6.0 * rho * macroU.at(0);
        f.at(7) = 0;
        f.at(5) = 0;
        f.at(0) = rho - f.at(1) - f.at(2) - f.at(3) - f.at(4) - f.at(6) - f.at(8);
    }
    // bottom left corner
    else if (position.at(0) == 0 && position.at(1) == yLen - 1)
    {
        rho = lattice.getCloseRho(position);

        f.at(1) = f.at(3) + 2.0 / 3.0 * rho * macroU.at(0);
        f.at(2) = f.at(4) + 2.0 / 3.0 * rho * macroU.at(1);
        f.at(5) = f.at(7) + 1.0 / 6.0 * rho * macroU.at(0) + 1.0 / 6.0 * rho * macroU.at(1);
        f.at(6) = 0;
        f.at(8) = 0;
        f.at(0) = rho - f.at(1) - f.at(2) - f.at(3) - f.at(4) - f.at(5) - f.at(7);
    }
    // top left corner
    else if (position.at(0) == 0 && position.at(1) == 0)
    {
        rho = lattice.getCloseRho(position);

        f.at(1) = f.at(3) + 2.0 / 3.0 * rho * macroU.at(0);
        f.at(4) = f.at(2) - 2.0 / 3.0 * rho * macroU.at(1);
        f.at(8) = f.at(6) + 1.0 / 6.0 * rho * macroU.at(0) - 1.0 / 6.0 * rho * macroU.at(1);
        f.at(7) = 0;
        f.at(5) = 0;
        f.at(0) = rho - f.at(1) - f.at(2) - f.at(3) - f.at(4) - f.at(6) - f.at(8);
    }
}

void Cell::bounceBackObstacle()
{
    if (obstacle)
        return;
    // regular bounce back
    if (boundary.at(0) == 1)
    {
        f.at(3) = newF.at(1);
    }
    if (boundary.at(0) == -1)
    {
        f.at(1) = newF.at(3);
    }
    if (boundary.at(1) == 1)
    {
        f.at(2) = newF.at(4);
    }
    if (boundary.at(1) == -1)
    {
        f.at(4) = newF.at(2);
    }
    if (boundary.at(2) == 1)
    {
        f.at(6) = newF.at(8);
    }
    if (boundary.at(2) == -1)
    {
        f.at(8) = newF.at(6);
    }
    if (boundary.at(3) == 1)
    {
        f.at(5) = newF.at(7);
    }
    if (boundary.at(3) == -1)
    {
        f.at(7) = newF.at(5);
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

void Cell::dragAndLift(float &drag, float &lift)
{
    const float adj = -2.0 / 0.04;
    float localDrag = 0;
    float localLift = 0;

    if (obstacle || (boundary.at(0) == 0 && boundary.at(1) == 0 && boundary.at(2) == 0 && boundary.at(3) == 0))
        return;

    if (boundary.at(0) == 1)
    {
        localDrag += newF.at(1) + f.at(3);
    }
    else if (boundary.at(0) == -1)
    {
        localDrag -= newF.at(3) + f.at(1);
    }
    if (boundary.at(1) == 1)
    {
        localLift += newF.at(4) + f.at(2);
    }
    else if (boundary.at(1) == -1)
    {
        localLift -= newF.at(2) + f.at(4);
    }
    if (boundary.at(2) == 1)
    {
        const float t = newF.at(8) + f.at(6);
        localDrag += t;
        localLift += t;
    }
    else if (boundary.at(2) == -1)
    {
        const float t = newF.at(6) + f.at(8);
        localDrag -= t;
        localLift -= t;
    }
    if (boundary.at(3) == 1)
    {
        const float t = newF.at(7) + f.at(5);
        localDrag -= t;
        localLift += t;
    }
    else if (boundary.at(3) == -1)
    {
        const float t = newF.at(5) + f.at(7);
        localDrag += t;
        localLift -= t;
    }

#pragma omp atomic
    drag += localDrag * adj;
#pragma omp atomic
    lift += localLift * adj;
}
