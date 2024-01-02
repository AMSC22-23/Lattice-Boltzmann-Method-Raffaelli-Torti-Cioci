#include "Lattice.hpp"
#include <iostream>

using namespace std;

int main(int argc, char const *argv[])
{
    // needs input file argument
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <input file> [-gpu]" << endl;
        return 1;
    }

    // set cmd line arguments
    string filename_in = argv[1];

    // check for -gpu flag
    bool gpuFlag = false;
    for(int i = 2; i < argc; i++) {
        if(string(argv[i]) == "-gpu") {
            gpuFlag = true;
            break;
        }
    }

    // create lattice
    Lattice lattice(filename_in);

    // open output file
    ofstream file_out;
    file_out.open("output.txt");
    if (!file_out.is_open())
    {
        throw runtime_error("Could not open file");
    }

    // write dimensions in first line
    for (int i = 0; i < lattice.getShape().size(); i++)
    {
        file_out << lattice.getShape().at(i) << ' ';
    }
    file_out << '\n';

    #ifdef USE_CUDA
    if (gpuFlag)
    {
        lattice.simulateGpu(file_out);
    }
    else
    {
        lattice.simulate(file_out);
    }
    #endif

    #ifndef USE_CUDA
    lattice.simulate(file_out);
    #endif

    // close output file
    file_out.close();

    return 0;
}
