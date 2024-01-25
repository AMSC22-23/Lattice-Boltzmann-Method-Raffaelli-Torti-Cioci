#ifndef GPUSIMULATION_CUH
#define GPUSIMULATION_CUH

#include "Cell.hpp"
#include "Structure.hpp"
#include "Utils.cpp"
#include <fstream>

namespace GpuSimulation
{
void cudaCaller(const NDimensionalMatrix<Cell> &cells, const float sigma, const float omP, const float omM,
                const int maxIt, const float uLid, const int problemType, const Structure &structure,
                const int plotSteps, std::ofstream &velocity_out, std::ofstream &lift_drag_out);
} // namespace GpuSimulation

#endif // GPUSIMULATION_CUH