cmake_minimum_required(VERSION 3.28)
project(4islaki)

set(CMAKE_CXX_STANDARD 20)


file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS src/*.cpp)
add_library(matrix ${SRC_FILES})
target_include_directories(matrix PUBLIC include)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_link_libraries(matrix PRIVATE asan)
endif()

add_executable(LU LU_Decomposition.cpp)
add_executable(run_throw run_throw.cpp)