cmake_minimum_required(VERSION 3.8.0)
project(Test_NumCpp VERSION 0.1.0)

add_executable(Test_NumCpp main.h main.cpp)

target_include_directories(Test_NumCpp
    PUBLIC include
)

set_property(DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}
    PROPERTY
    VS_STARTUP_PROJECT
    Test_NumCpp
)
