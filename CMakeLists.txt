cmake_minimum_required(VERSION 3.10)
project(MatrixAlgebra C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

include_directories(inc)

add_executable(matrix_alg
    src/algebra.c
    src/main.c
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)