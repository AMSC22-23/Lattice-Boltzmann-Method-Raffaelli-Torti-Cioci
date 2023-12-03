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

    // create lattice
    Lattice lattice(filename);

    // open output file
    ofstream file;
    file.open("output.txt");
    if (!file.is_open())
    {
        throw runtime_error("Could not open file");
    }

    // write dimensions in first line
    for (int i = 0; i < lattice.getShape().size(); i++)
    {
        file << lattice.getShape().at(i) << ' ';
    }
    file << '\n';

    // update lattice for timeSteps
    for (int i = 0; i < timeSteps; i++)
    {
        lattice.update(deltaTime, file);
    }

    // close output file
    file.close();

    return 0;
}

// !  add dimensions
