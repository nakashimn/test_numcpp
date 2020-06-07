cmake_minimum_required(VERSION 3.8.0)
project(blas VERSION 0.1.0)

add_executable(blas main.h main.cpp)


target_include_directories(blas
    PUBLIC include
)
