cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Example 1
project(ex1 VERSION 1.0)

add_executable(ex1 ${CMAKE_CURRENT_SOURCE_DIR}/src/newton.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/example1.cpp)

target_include_directories(ex1 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Example 2
project(ex2 VERSION 1.0)

add_executable(ex2 ${CMAKE_CURRENT_SOURCE_DIR}/src/newton.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/example2.cpp)

target_include_directories(ex2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Example 3
project(ex3 VERSION 1.0)

add_executable(ex3 ${CMAKE_CURRENT_SOURCE_DIR}/src/newton.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/example3.cpp)

target_include_directories(ex3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)