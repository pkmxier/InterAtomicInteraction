# see https://cmake.org/cmake/help/latest/module/FindBoost.html
find_package(Boost REQUIRED)

include_directories(${Boost_INCLUDE_DIR})

cmake_minimum_required(VERSION 3.13)  # CMake version check
project(atomicEnergy)               # Create project "simple_example"
set(CMAKE_CXX_STANDARD 14)            # Enable c++14 standard

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

# Add main.cpp file of project root directory as source file
set(SOURCE_FILES main.cpp atom.cpp vector3d.cpp Solver.cpp Solver.h Lattice.cpp Lattice.h Optimizer.cpp Optimizer.h Timer.cpp Timer.h Parameters.cpp Parameters.h)

# Add executable target with source files listed in SOURCE_FILES variable
add_executable(atomicEnergy ${SOURCE_FILES})