# FILEPATH: /home/andre/Lattice-Boltzmann-Method-Raffaelli-Torti-Cioci/CMakeLists.txt

cmake_minimum_required(VERSION 3.20)

project(lbm)

set(CMAKE_CXX_STANDARD 23)

find_package(TBB REQUIRED)
# find openmp package
find_package(OpenMP REQUIRED)

# Add the executable target
add_executable(lbm src/Main.cpp src/Lattice.cpp src/Cell.cpp src/Utils.cpp src/Structure.cpp)
target_link_libraries(lbm PRIVATE TBB::tbb)
target_link_libraries(lbm PRIVATE OpenMP::OpenMP_CXX)
