#ifndef GPUSIMULATION_CUH
#define GPUSIMULATION_CUH

#include <fstream>
#include "Cell.hpp"
#include "Utils.cpp"
#include "Structure.hpp"

namespace GpuSimulation
{
void cudaCaller(const NDimensionalMatrix<Cell> &cells, const float sigma, const float omP, const float omM,
                const int maxIt, const float uLid, const int problemType, const Structure &structure,
                std::ofstream &file);
} // namespace GpuSimulation

#endif // GPUSIMULATION_CUH