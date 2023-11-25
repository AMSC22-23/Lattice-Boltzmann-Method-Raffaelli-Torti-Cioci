#ifndef CELL_HPP
#define CELL_HPP

#include <vector>
class Lattice;

class Cell
{
  public:
    Cell() = default;
    ~Cell() = default;
    void update(const float deltaTime, Lattice &lattice);
    void setFAtIndex(const unsigned int index, const float value);
    void setObstacle();
    bool isObstacle();

  private:
    void updateFeq(std::vector<float> &feq, const float &ux, const float &uy, const float &rho);
    void collision(const std::vector<float> &feq, std::vector<float> &f, const float dt);
    void streaming(const std::vector<float> fstar, Lattice &lattice);
    std::vector<float> f;                     // Distribution function
    std::vector<float> feq;                   // Equilibrium Distribution function
    std::vector<float> Omega;                 // Collision operator
    bool obstacle;                            // Is this cell an obstacle?
    float uX, uY;                             // Macroscopic velocity
    float momX, momY;                         // Momentum density
    float rho;                                // Macroscopic density
    const std::vector<unsigned int> position; // cell position in the lattice
    const int boundary[2];                    // array of two elements that describes the boundary of the cell
};
//boundary is an array of two elements; eache element can be 0, 1 or -1.
//if we have a boundary {-1, 1} it means that we can't go in the "opposite direction of the x axis" (not on left)
//and we can't go in the direction of the y axis (not up).
//sum up: the first element of the vector tells us which direction on x we can't follow in that cell and
//the second element tells us which direction on y we can't follow. If one of the element is 0 it means that
//we have no boundary in this direction

#endif // CELL_HPP