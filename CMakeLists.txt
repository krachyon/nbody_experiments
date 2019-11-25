cmake_minimum_required(VERSION 3.15)
project(experiment)

set(PYBIND11_CPP_STANDARD -std=c++17)
set(CMAKE_CXX_STANDARD 17)

set(PYTHON_EXECUTABLE "/usr/bin/python")
set(PYTHON_LIBRARY "/usr/lib/libpython.so")
find_package(Python COMPONENTS Development)
find_package(pybind11 REQUIRED)
find_package(Eigen3 REQUIRED)





find_package(OpenMP)
find_package(Boost)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Python_INCLUDE_DIRS})

pybind11_add_module(nbody nbody.h module.cpp)
add_executable(runNbody nbody.h executable.cpp)

target_link_libraries(nbody PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(runNbody PUBLIC OpenMP::OpenMP_CXX)




set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
