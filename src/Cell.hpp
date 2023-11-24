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
};

#endif // CELL_HPP