cmake_minimum_required(VERSION 3.20)
project(Lab_3)

set(CMAKE_CXX_STANDARD 14)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(Lab_3 main.cpp)

target_link_libraries(Lab_3 ${MPI_C_LIBRARIES})
