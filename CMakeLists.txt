cmake_minimum_required(VERSION 4.0)
project(VascularCircuit0D LANGUAGES CXX)
add_executable(solver_serial.exe solver_serial.cpp)

# Copy input files to build directory for easier execution
configure_file(${CMAKE_SOURCE_DIR}/input/input.inp ${CMAKE_BINARY_DIR}/input/input.inp COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/input/AAo.csv ${CMAKE_BINARY_DIR}/input/AAo.csv COPYONLY)