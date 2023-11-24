#include "Lattice.hpp"
#include <iostream>

using namespace std;

int main()
{
    std::cout << "Hello, World!" << std::endl;

    Lattice lattice("data/lid-driven-cavity.txt");

    cout << "Lattice shape: ";
    for (auto i : lattice.getShape())
    {
        cout << i << " ";
    }
    cout << endl;

    return 0;
}