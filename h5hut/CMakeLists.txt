cmake_minimum_required(VERSION 3.12)
project(H5hut)

# build configuration
option(USE_FORTRAN "Build Fortran library" ON)
option(USE_PYTHON "Build Python interface" OFF)
option(USE_PARALLEL "Build parallel IO version of library" OFF)

# tell CMake we're using Fortran
if (USE_FORTRAN)
  enable_language(Fortran)
endif (USE_FORTRAN)

# find HDF5 and ensure it matches our parallel setting
set(HDF5_PREFER_PARALLEL ${USE_PARALLEL})
find_package(HDF5 REQUIRED)
if (USE_PARALLEL AND NOT HDF5_IS_PARALLEL)
  message(FATAL_ERROR "Requested parallel build, but no parallel HDF5 library available")
endif (USE_PARALLEL AND NOT HDF5_IS_PARALLEL)

# additionally build with MPI if parallel is requested
if (USE_PARALLEL)
  find_package(MPI REQUIRED)
  # deal with bug in Xenial OpenMPI
  string(STRIP ${MPI_C_LINK_FLAGS} MPI_C_LINK_FLAGS)
endif (USE_PARALLEL)

# required for building Python library
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# core h5hut library
include_directories("${PROJECT_SOURCE_DIR}/src/include")
add_subdirectory(src/h5core)

# fortran library and module
if (USE_FORTRAN)
  add_subdirectory(src/Fortran)
endif (USE_FORTRAN)

# python module
if (USE_PYTHON)
  cmake_minimum_required(VERSION 3.14)
  add_subdirectory(src/Python)
endif (USE_PYTHON)
