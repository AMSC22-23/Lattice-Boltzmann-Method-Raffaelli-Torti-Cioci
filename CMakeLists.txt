# FILEPATH: /home/andre/Lattice-Boltzmann-Method-Raffaelli-Torti-Cioci/CMakeLists.txt

cmake_minimum_required(VERSION 3.22)

project(lbm)

set(CMAKE_CXX_STANDARD 23)

# Add the executable target
add_executable(lbm src/Main.cpp src/Lattice.cpp src/Cell.cpp src/Utils.cpp)
