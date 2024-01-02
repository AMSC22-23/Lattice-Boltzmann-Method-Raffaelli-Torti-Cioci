cmake_minimum_required(VERSION 3.20)


set(CMAKE_CXX_STANDARD 23)


project(lbm CXX)

include_directories(include)

# find openmp package
find_package(OpenMP REQUIRED)
# find cuda package
find_package(CUDA)
# find tbb
find_package(TBB REQUIRED)

if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -G")
endif ()

# Add the library target for CUDA-specific code
if(CUDA_FOUND)
    set(CMAKE_CUDA_STANDARD 23)    
    set(CUDAToolkit_ROOT "/usr/local/cuda")
    set(CMAKE_CUDA_COMPILER "${CUDAToolkit_ROOT}/bin/nvcc")
    enable_language(CUDA)
    set(CMAKE_CUDA_ARCHITECTURES 80)
    include_directories(${CUDA_INCLUDE_DIRS})
    add_library(lbm_cuda_lib src/GpuSimulation.cu)
    target_link_libraries(lbm_cuda_lib PRIVATE ${CUDA_LIBRARIES} OpenMP::OpenMP_CXX TBB::tbb)
    set_target_properties(lbm_cuda_lib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_compile_definitions(lbm_cuda_lib PRIVATE -DUSE_CUDA)
endif()
    
# Add the main executable target
add_executable(lbm src/Main.cpp src/Lattice.cpp src/Cell.cpp src/Structure.cpp include/Utils.cpp)
    
# Link the main executable with the CUDA-specific library if CUDA is found
if(CUDA_FOUND)
    target_link_libraries(lbm PRIVATE lbm_cuda_lib)
    target_compile_definitions(lbm PRIVATE -DUSE_CUDA)
endif()

target_link_libraries(lbm PRIVATE OpenMP::OpenMP_CXX TBB::tbb)
