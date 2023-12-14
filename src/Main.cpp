#include "Lattice.hpp"
#include <iostream>

using namespace std;

int main(int argc, char const *argv[])
{
    // --fast flag sets f to feq IGNORED FOR NOW
    // may be automatic if deltaTime is big enough (5 times TAU)

    // needs at least 3 arguments: input file, delta time, number of time steps
    if (argc < 4)
    {
        cout << "Usage: " << argv[0] << " <input file> <delta time> <number of time steps>" << endl;
        return 1;
    }

    // set cmd line arguments
    string filename = argv[1];
    float deltaTime = stof(argv[2]);
    int timeSteps = stoi(argv[3]);


    return 0;
}

// !  add dimensions
