cmake_minimum_required( VERSION 3.10 )

project( project )

include_directories( ./src )

add_executable( project ./src/main.cpp ./src/solver.cpp )