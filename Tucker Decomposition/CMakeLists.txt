cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(tucker VERSION 1.0)

file(GLOB_RECURSE SRC_FILES src/*.cpp)
add_executable(tucker ${SRC_FILES})

target_include_directories(tucker PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)