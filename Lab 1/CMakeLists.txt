cmake_minimum_required(VERSION 3.20)
project (Lab_1)

set(CMAKE_CXX_STANDARD 14)

find_package(MPI REQUIRED)

add_executable(Lab_1 main.cpp)
target_link_libraries(Lab_1 PUBLIC MPI::MPI_CXX)

