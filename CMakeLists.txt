cmake_minimum_required(VERSION 3.28)
project(4islaki)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_executable(LU LU_Decomposition.cpp)

add_executable(run_throw run_throw.cpp)

add_executable(Zeidel Zeidel.cpp)