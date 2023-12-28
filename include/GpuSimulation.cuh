#ifndef GPUSIMULATION_CUH
#define GPUSIMULATION_CUH

namespace GpuSimulation
{
__global__ void step1(const int nx, const int ny, const int it, const int problem_type, const float u_lid,
                      const float om_p, const float halfOmpOmmSum, const float halfOmpOmmSub, float *f, float *new_f, float *rho, float *ux, float *uy,
                      int *boundary, bool *obstacle);
__global__ void step2(const int nx, const int ny, float *f, float *new_f, int *boundary, bool *obstacle);
} // namespace GpuSimulation

#endif // GPUSIMULATION_CUH