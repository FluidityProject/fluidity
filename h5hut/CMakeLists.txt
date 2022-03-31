cmake_minimum_required(VERSION 3.14)
project(H5hut VERSION 2.0.0)

option(H5HUT_USE_FORTRAN "Build Fortran library" ON)
option(H5HUT_USE_PYTHON "Build Python interface" OFF)

add_subdirectory(src/h5core)

if(H5HUT_USE_FORTRAN)
  enable_language(Fortran)
  add_subdirectory(src/Fortran)
endif()

if(H5HUT_USE_PYTHON)
  add_subdirectory(src/Python)
endif()
