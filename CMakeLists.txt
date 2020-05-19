cmake_minimum_required(VERSION 3.16)
project(CalculationMethodsLab2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(SRC_LIST
        Sources/Main.cpp
        Sources/Matrix.cpp
        Sources/Utils.cpp
        Sources/NonlinearSolvers.cpp
        Sources/QRalgorithm.cpp Sources/QRalgorithm.h)

include_directories(Sources)

add_executable(${PROJECT_NAME} ${SRC_LIST})