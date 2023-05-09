cmake_minimum_required( VERSION 3.10 )

project( expansion )

include_directories( ./src )

add_executable( expansion ./src/main.cpp ./src/solver.cpp )