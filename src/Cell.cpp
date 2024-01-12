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
           const std::vector<float> &_f, const std::vector<int> &_position)
    : boundary(_boundary), obstacle(_obstacle), f(_f), position(_position)
{
    rho = 1.0;
    macroU = std::vector<float>(structure.dimensions, 0.0);
    newF = std::vector<float>(structure.velocity_directions, 0.0);
    if (obstacle)
    {
        rho = 0;
        for(int i = 0; i < 9; i++)
        {
            f.at(i) = 0;
        }
    }
}

void Cell::updateMacro(const Structure &structure)
{
    if (obstacle)
    {
        macroU.at(0) = 0;
        macroU.at(1) = 0;
        rho = 0;
        return;
    }
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

    if (obstacle)
        return;

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

void Cell::streaming(Lattice &lattice)
{
    if (obstacle)
        return;

    const Structure &structure = lattice.getStructure();

    // stream for index 0
    f.at(0) = newF.at(0);
    bool stream;
    int new_position[structure.dimensions];

    // stream for other indices
    for (int i = 1; i < structure.velocity_directions; i++)
    {

        // calculate new position
        for (int j = 0; j < structure.dimensions; j++)
        {
            new_position[j] = position.at(j) + (int)structure.velocities_by_direction.at(i).at(j);
        }

        // check if new position is inside the lattice and is not an obstacle
        if (new_position[0] >= 0 && new_position[0] < lattice.getShape().at(0) && new_position[1] >= 0 &&
            new_position[1] < lattice.getShape().at(1) && !lattice.getCellAtIndices(new_position).isObstacle())
        {
            lattice.getCellAtIndices(new_position).setFAtIndex(i, newF.at(i));
        }
    }
}

void Cell::setInlets(Lattice &lattice, const float uLidNow)
{
    if (obstacle)
    {
        macroU.at(0) = 0; //! magari farli diventare nan
        macroU.at(1) = 0;
        return;
    }
    // 2d only
    const int xLen = lattice.getShape().at(0);
    const int yLen = lattice.getShape().at(1);
    const Structure &structure = lattice.getStructure();
    const int problemType = lattice.getProblemType();
    switch (problemType)
    {
    case 1:
        // if i'm at any boundary set macroU to 0
        if (position.at(0) == 0 || position.at(1) == 0 || position.at(0) == xLen - 1 || position.at(1) == yLen - 1)
        {
            macroU.at(0) = 0;
            macroU.at(1) = 0;
        }
        // if i'm at top wall set macroU.x to uLidNow
        if (position.at(1) == 0)
            macroU.at(0) = uLidNow;

        break;
    case 2:
    case 3:
        if ((position.at(0) != 0 || position.at(0) != xLen) && (position.at(1) != 0 || position.at(1) != yLen))
            break;
        if (position.at(0) == 0)
        {
            const float halfDim = static_cast<float>(yLen) / 2.0;
            const float temp = (static_cast<float>(position.at(1)) / halfDim) - 1.0;
            const float mul = 1.0 - temp * temp;
            macroU.at(0) = uLidNow * 2.5 * mul;
            macroU.at(1) = 0;
        }

        if (position.at(0) == xLen)
            macroU.at(1) = 0;

        if (position.at(1) == 0 || position.at(1) == yLen)
            if (position.at(0) != 0)
                macroU.at(0) = 0;
    
        if(boundary.at(0) != 0)//no slip condition
            macroU.at(0) =0;

        if(boundary.at(1) != 0)//no slip condition
            macroU.at(1) =0;

        break;
    default:
        break;
    }
}

void Cell::zouHe(Lattice &lattice, const float uLidNow)
{
    if (obstacle)
        return;
    //! ora controllo: se sono una cella centrale esco 
    if((position.at(0) != 0 && position.at(0) != lattice.getShape().at(0) - 1) && (position.at(1) != 0 && position.at(1) != lattice.getShape().at(1) - 1))
        return;

    // only 2D...
    const int xLen = lattice.getShape().at(0) - 1;
    const int yLen = lattice.getShape().at(1) - 1;
    const int problemType = lattice.getProblemType();

    if (position.at(0) != 0 && position.at(0) != xLen && position.at(1) == 0) // top wall
    {
        if (problemType == 2 || problemType == 3)
            macroU.at(0) = 0; //!
            
        rho = (f.at(0) + f.at(1) + f.at(3) + 2.0 * (f.at(2) + f.at(5) + f.at(6))) / (1.0 + macroU.at(1));
        f.at(4) = f.at(2) - 2.0 / 3.0 * rho * macroU.at(1);
        f.at(7) = f.at(5) + 0.5 * (f.at(1) - f.at(3)) - 0.5 * rho * macroU.at(0) - 1.0 / 6.0 * rho * macroU.at(1);
        f.at(8) = f.at(6) - 0.5 * (f.at(1) - f.at(3)) + 0.5 * rho * macroU.at(0) - 1.0 / 6.0 * rho * macroU.at(1);
    }
    else if (position.at(0) == xLen && position.at(1) != 0 && position.at(1) != yLen) // right wall
    {
        switch (problemType)
        {
        case 1:
            rho = (f.at(0) + f.at(2) + f.at(4) + 2.0 * (f.at(1) + f.at(5) + f.at(8))) / (1.0 + macroU.at(0));
            break;
        case 2:
        case 3:
            //rho = lattice.getCloseRho(position);
            macroU.at(1) = 0; //! potrebbe essere superfluo
            //macroU.at(0) = lattice.getCloseU(position).at(0);
            break;
        default:
            // Handle other cases or throw an error
            break;
        }
        f.at(3) = f.at(1) - 2.0 / 3.0 * rho * macroU.at(0);
        f.at(6) = f.at(8) - 0.5 * (f.at(2) - f.at(4)) - 1.0 / 6.0 * rho * macroU.at(0) + 0.5 * rho * macroU.at(1);
        f.at(7) = f.at(5) + 0.5 * (f.at(2) - f.at(4)) - 1.0 / 6.0 * rho * macroU.at(0) - 0.5 * rho * macroU.at(1);
    }
    else if (position.at(0) != 0 && position.at(0) != xLen && position.at(1) == yLen) // bottom wall
    {
        if (problemType == 2 || problemType == 3)
            macroU.at(0) = 0; //!
        
        rho = (f.at(0) + f.at(1) + f.at(3) + 2.0 * (f.at(4) + f.at(7) + f.at(8))) / (1.0 - macroU.at(1));
        f.at(2) = f.at(4) + 2.0 / 3.0 * rho * macroU.at(1);
        f.at(5) = f.at(7) - 0.5 * (f.at(1) - f.at(3)) + 0.5 * rho * macroU.at(0) + 1.0 / 6.0 * rho * macroU.at(1);
        f.at(6) = f.at(8) + 0.5 * (f.at(1) - f.at(3)) - 0.5 * rho * macroU.at(0) + 1.0 / 6.0 * rho * macroU.at(1);
    }
    else if (position.at(0) == 0 && position.at(1) != 0 && position.at(1) != yLen) // left wall
    {

        if (problemType == 2 || problemType == 3)
            macroU.at(1) = 0;
        
        rho = (f.at(0) + f.at(2) + f.at(4) + 2.0 * (f.at(3) + f.at(6) + f.at(7))) / (1.0 - macroU.at(0));
        f.at(1) = f.at(3) - 2.0 / 3.0 * rho * macroU.at(0);
        f.at(5) = f.at(7) - 0.5 * (f.at(2) - f.at(4)) + 1.0 / 6.0 * rho * macroU.at(0) + 0.5 * rho * macroU.at(1);
        f.at(8) = f.at(6) + 0.5 * (f.at(2) - f.at(4)) + 1.0 / 6.0 * rho * macroU.at(0) - 0.5 * rho * macroU.at(1);
    }
    else if (position.at(0) == xLen && position.at(1) == 0) // top right corner
    {
        if (problemType == 2 || problemType == 3)
        {
            macroU.at(0) = lattice.getCloseU(position).at(0);
            macroU.at(1) = 0;
            rho = lattice.getCloseRho(position);
        }
        f.at(3) = f.at(1) - 2.0 / 3.0 * rho * macroU.at(0);
        f.at(4) = f.at(2) - 2.0 / 3.0 * rho * macroU.at(1);
        f.at(7) = f.at(5) - 1.0 / 6.0 * rho * macroU.at(0) - 1.0 / 6.0 * rho * macroU.at(1);
        f.at(8) = 0;
        f.at(6) = 0;
        f.at(0) = rho - f.at(1) - f.at(2) - f.at(3) - f.at(4) - f.at(5) - f.at(7);
    }
    else if (position.at(0) == xLen && position.at(1) == yLen) // bottom right corner
    {
        if (problemType == 2 || problemType == 3)
        {
            macroU.at(0) = lattice.getCloseU(position).at(0);
            macroU.at(1) = 0;
            rho = lattice.getCloseRho(position);
        }
        f.at(3) = f.at(1) - 2.0 / 3.0 * rho * macroU.at(0);
        f.at(2) = f.at(4) + 2.0 / 3.0 * rho * macroU.at(1);
        f.at(6) = f.at(8) + 1.0 / 6.0 * rho * macroU.at(1) - 1.0 / 6.0 * rho * macroU.at(0);
        f.at(7) = 0;
        f.at(5) = 0;
        f.at(0) = rho - f.at(1) - f.at(2) - f.at(3) - f.at(4) - f.at(6) - f.at(8);
    }
    else if (position.at(0) == 0 && position.at(1) == yLen) // bottom left corner
    {
        if (problemType == 2 || problemType == 3)
        {
            macroU.at(0) = lattice.getCloseU(position).at(0);
            macroU.at(1) = 0;
            rho = lattice.getCloseRho(position);
        }
        f.at(1) = f.at(3) + 2.0 / 3.0 * rho * macroU.at(0);
        f.at(2) = f.at(4) + 2.0 / 3.0 * rho * macroU.at(1);
        f.at(5) = f.at(7) + 1.0 / 6.0 * rho * macroU.at(0) + 1.0 / 6.0 * rho * macroU.at(1);
        f.at(6) = 0;
        f.at(8) = 0;
        f.at(0) = rho - f.at(1) - f.at(2) - f.at(3) - f.at(4) - f.at(5) - f.at(7);
    }
    else if (position.at(0) == 0 && position.at(1) == 0) // top left corner
    {
        if (problemType == 2 || problemType == 3)
        {
            macroU.at(0) = lattice.getCloseU(position).at(0);
            macroU.at(1) = 0;
            rho = lattice.getCloseRho(position);
        }
        f.at(1) = f.at(3) + 2.0 / 3.0 * rho * macroU.at(0);
        f.at(4) = f.at(2) - 2.0 / 3.0 * rho * macroU.at(1);
        f.at(8) = f.at(6) - 1.0 / 6.0 * rho * macroU.at(0) + 1.0 / 6.0 * rho * macroU.at(1);
        f.at(7) = 0;
        f.at(5) = 0;
        f.at(0) = rho - f.at(1) - f.at(2) - f.at(3) - f.at(4) - f.at(6) - f.at(8);
    }
}

void Cell::bounce_back_obstacle() //! chiamare solo per le celle che non stanno al bordo (per quelle si usa ZouHe, inoltre chiamare zouHe solo per quelle al bordo)
{
    if (obstacle)
        return;
    // regular bounce back
    if (boundary.at(0) == 0 && boundary.at(1) == 0 && boundary.at(2) == 0 && boundary.at(3) == 0)
        return;
    // f.at(0) is already okay
    // first check the x and y directions
    if (boundary.at(0) == 1)
    {
        macroU.at(0) = 0; //! no slip condition
        f.at(6) = newF.at(8); //! ho cambiato tutte le f con newf percè lui usava quelle
        f.at(3) = newF.at(1);
        f.at(7) = newF.at(5);
        f.at(8) = 0;
        f.at(1) = 0;
        f.at(5) = 0;
    }
    if (boundary.at(0) == -1)
    {
        macroU.at(0) = 0;
        f.at(8) = newF.at(6);
        f.at(1) = newF.at(3);
        f.at(5) = newF.at(7);
        f.at(6) = 0;
        f.at(3) = 0;
        f.at(7) = 0;
    }
    if (boundary.at(1) == 1)
    {
        macroU.at(1) = 0;
        f.at(6) = newF.at(8);
        f.at(2) = newF.at(4);
        f.at(5) = newF.at(7);
        f.at(8) = 0;
        f.at(4) = 0;
        f.at(7) = 0;
    }
    if (boundary.at(1) == -1)
    {
        macroU.at(1) = 0;
        f.at(8) = newF.at(6);
        f.at(4) = newF.at(2);
        f.at(7) = newF.at(5);
        f.at(6) = 0;
        f.at(2) = 0;
        f.at(5) = 0;
    }
    if (boundary.at(2) != 0 && boundary.at(0) == 0 &&
        boundary.at(1) == 0) // spigolo in fuori, solo una velocità rimbalza
    {
        if (boundary.at(2) == 1)
        {
            f.at(6) = newF.at(8);
            f.at(8) = 0;
        }
        else
        {
            f.at(8) = newF.at(6);
            f.at(6) = 0;
        }
    }
    if (boundary.at(3) != 0 && boundary.at(0) == 0 &&
        boundary.at(1) == 0) // spigolo in fuori, solo una velocità rimbalza
    {
        if (boundary.at(3) == 1)
        {
            f.at(5) = newF.at(7);
            f.at(7) = 0;
        }
        else
        {
            f.at(7) = newF.at(5);
            f.at(5) = 0;
        }
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