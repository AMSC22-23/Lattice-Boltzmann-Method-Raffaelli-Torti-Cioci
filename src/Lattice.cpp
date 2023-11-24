#include "Lattice.hpp"

Lattice::Lattice(std::string filename) {
    
}

Lattice::~Lattice() {}

void Lattice::update(const float deltaTime) {
    
}

Cell& Lattice::getCellAtIndex(std::vector<unsigned int> index) {
    return cells.getElement(index);
}