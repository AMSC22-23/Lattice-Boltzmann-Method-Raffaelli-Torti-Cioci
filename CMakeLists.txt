cmake_minimum_required(VERSION 3.20)

project(lbm)

set(CMAKE_CXX_STANDARD 23)

include_directories(include)

# find openmp package
find_package(OpenMP REQUIRED)
# find cuda package
find_package(CUDA)


# Add the executable target
add_executable(lbm src/Main.cpp src/Lattice.cpp src/Cell.cpp src/Structure.cpp include/Utils.cpp)
target_link_libraries(lbm PRIVATE OpenMP::OpenMP_CXX)
