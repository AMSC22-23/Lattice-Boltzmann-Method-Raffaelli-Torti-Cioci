#include "Lattice.hpp"
#include <iostream>

using namespace std;

int main(int argc, char const *argv[])
{
    // needs input file argument
    if (argc < 3)
    {
        cout << "Usage: " << argv[0] << " <input file> <plot frames> [-gpu]" << endl;
        return 1;
    }

    // set cmd line arguments
    const string filename_in = argv[1];
    const int plotFrames = atoi(argv[2]);

    // check for -gpu flag
    bool gpuFlag = false;
    for (int i = 3; i < argc; i++)
    {
        if (string(argv[i]) == "-gpu")
        {
            gpuFlag = true;
            break;
        }
    }

    // open input file
    ifstream file_in(filename_in);
    if (!file_in.is_open())
    {
        throw runtime_error("Could not open input file");
    }

    // create lattice
    Lattice lattice(file_in, plotFrames);

    // close input file
    file_in.close();

    // open output files
    ofstream velocity_out("outputs/velocity_out.txt");
    ofstream lift_drag_out("outputs/lift_drag_out.txt");
    if (!velocity_out.is_open() || !lift_drag_out.is_open())
    {
        throw runtime_error("Could not create output files");
    }

    // write dimensions in first line
    for (int i = 0; i < (int)lattice.getShape().size(); i++)
    {
        velocity_out << lattice.getShape().at(i) << ' ';
    }
    velocity_out << '\n';

#ifdef USE_CUDA
    if (gpuFlag)
    {
        lattice.simulateGpu(velocity_out, lift_drag_out);
    }
    else
    {
        lattice.simulate(velocity_out, lift_drag_out);
    }
#endif

#ifndef USE_CUDA
    lattice.simulate(velocity_out, lift_drag_out);
#endif

    // close output files
    velocity_out.close();
    lift_drag_out.close();

    return 0;
}
