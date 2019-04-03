cmake_minimum_required(VERSION 3.13)
project(LinearAlgebra)

set(CMAKE_CXX_STANDARD 17)

add_library(LinearAlgebra SHARED linearalgebra.cpp linearalgebra.h)
add_executable(TEST test.cpp)
target_link_libraries(TEST LinearAlgebra)
