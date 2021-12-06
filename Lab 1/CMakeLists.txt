cmake_minimum_required(VERSION 3.20)
project(Lab_1)

include_directories( ./src )
set(CMAKE_CXX_STANDARD 14)
include_directories(/MPI)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(Lab_1 main.cpp BlockSchemaAlgorithm.cpp BlockSchemaAlgorithm.h HorisontalStripesAlgorithm.cpp HorisontalStripesAlgorithm.h SerialAlgorithm.cpp SerialAlgorithm.h VerticalStripesAlgorithm.cpp VerticalStripesAlgorithm.h)

target_link_libraries(Lab_1 ${MPI_C_LIBRARIES})
