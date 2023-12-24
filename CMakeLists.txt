cmake_minimum_required(VERSION 3.20)

set(CUDAToolkit_ROOT "/usr/local/cuda")
set(CMAKE_CUDA_COMPILER "${CUDAToolkit_ROOT}/bin/nvcc")
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CUDA_STANDARD 23)
set(CMAKE_CUDA_ARCHITECTURES 80)

project(lbm CUDA CXX)

include_directories(include)

# find openmp package
find_package(OpenMP REQUIRED)
# find cuda package
find_package(CUDA)

# Add the executable target
add_executable(lbm src/Main.cpp src/Lattice.cu src/Cell.cpp src/Structure.cpp include/Utils.cpp src/GpuSimulation.cu)
target_link_libraries(lbm PRIVATE OpenMP::OpenMP_CXX)
