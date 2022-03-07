cmake_minimum_required(VERSION 3.14)
project(H5hut VERSION 2.0.0)

option(HDF5_PREFER_PARALLEL "Prefer parallel HDF5" ON)
option(USE_FORTRAN "Build Fortran library" ON)
option(USE_PYTHON "Build Python interface" OFF)

find_package(MPI REQUIRED COMPONENTS C)
find_package(HDF5 REQUIRED)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include_directories("${PROJECT_SOURCE_DIR}/src/include")
add_subdirectory(src/h5core)

if(USE_FORTRAN)
  enable_language(Fortran)
  add_subdirectory(src/Fortran)
endif(USE_FORTRAN)

if(USE_PYTHON)
  add_subdirectory(src/Python)
endif(USE_PYTHON)
