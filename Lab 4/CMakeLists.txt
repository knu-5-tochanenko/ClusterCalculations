cmake_minimum_required(VERSION 3.20)
project(Lab_4)

set(CMAKE_CXX_STANDARD 14)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(Lab_4 main.cpp)

target_link_libraries(Lab_4 ${MPI_C_LIBRARIES})
