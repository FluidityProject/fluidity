cmake_minimum_required(VERSION 3.19)
project(fluidity
  VERSION 4.1.19
  HOMEPAGE_URL "http://fluidityproject.github.io/"
  LANGUAGES C CXX Fortran
)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Set build type; \
typical values include Debug, Release, RelWithDebInfo and MinSizeRel." FORCE
  )
endif()

set(CMAKE_PREFIX_PATH "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

set(CMAKE_CXX_STANDARD 17)

option(HDF5_PREFER_PARALLEL "Look for parallel HDF5 first" ON)

include(ExternalProject)
include(FetchContent)

find_package(PkgConfig REQUIRED)

# Check for VTK
find_package(VTK REQUIRED)

# Check for Python 3
find_package(Python3 3.6...<3.10 REQUIRED Interpreter Development NumPy)
# # Check for assess Python package
# execute_process(
#   COMMAND ${Python3_EXECUTABLE} -c "import assess"
#   RESULT_VARIABLE Python3_ASSESS
# )

# Check for CGAL
find_package(CGAL REQUIRED)

# Check for PETSc
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
message(STATUS "Found PETSc ${PETSC_VERSION}")

# Check for NetCDF
pkg_check_modules(NETCDF REQUIRED netcdf netcdf-fortran)

# Check for Zoltan
find_path(ZOLTAN_INCLUDE_DIR zoltan.mod REQUIRED)
find_library(ZOLTAN_LIBRARY zoltan REQUIRED)

# Check for UDUNITS
find_library(UDUNITS_LIBRARY udunits2 REQUIRED)

# Check for ExodusII
find_library(EXODUSII_LIBRARY exoIIv2c REQUIRED)

# Add libspatialindex
FetchContent_Declare(libspatialindex
  GIT_REPOSITORY https://github.com/libspatialindex/libspatialindex
  GIT_TAG origin/master
  LOG_DOWNLOAD ON
  LOG_CONFIGURE ON
  LOG_BUILD ON
  LOG_INSTALL ON
)
FetchContent_MakeAvailable(libspatialindex)

# Add Judy
ExternalProject_Add(libjudy
  PREFIX libjudy
  URL https://sourceforge.net/projects/judy/files/judy/Judy-1.0.5/Judy-1.0.5.tar.gz/download
  CONFIGURE_COMMAND ./configure --enable-64-bit --prefix=<SOURCE_DIR>
  BUILD_COMMAND make
  BUILD_IN_SOURCE TRUE
  INSTALL_COMMAND make install exec_prefix=<SOURCE_DIR>
  TEST_COMMAND make check
  TEST_BEFORE_INSTALL TRUE
  LOG_DOWNLOAD ON
  LOG_CONFIGURE ON
  LOG_BUILD ON
  LOG_INSTALL ON
  LOG_TEST ON
)
ExternalProject_Add_Step(libjudy byproducts
  DEPENDEES install
  BYPRODUCTS <SOURCE_DIR>/lib/libJudy.a
  LOG ON
)
ExternalProject_Get_property(libjudy SOURCE_DIR)
set(libjudy_SOURCE_DIR ${SOURCE_DIR})
add_library(Judy STATIC IMPORTED)
set_target_properties(Judy
  PROPERTIES IMPORTED_LOCATION ${libjudy_SOURCE_DIR}/lib/libJudy.a
)

# Generate some source files
add_custom_command(
  OUTPUT ${CMAKE_SOURCE_DIR}/preprocessor/register_diagnostics.F90
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/make_register_diagnostics.py
  MAIN_DEPENDENCY ${CMAKE_SOURCE_DIR}/tools/make_register_diagnostics.py
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  VERBATIM
)
add_custom_command(
  OUTPUT ${CMAKE_SOURCE_DIR}/preprocessor/check_options.F90
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/make_check_options.py
  MAIN_DEPENDENCY ${CMAKE_SOURCE_DIR}/tools/make_check_options.py
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  VERBATIM
)
add_custom_command(
  OUTPUT ${CMAKE_SOURCE_DIR}/diagnostics/Diagnostic_Fields_New.F90
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/make_diagnostic_fields.py
  MAIN_DEPENDENCY ${CMAKE_SOURCE_DIR}/tools/make_diagnostic_fields.py
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/diagnostics
  VERBATIM
)

file(WRITE "${CMAKE_BINARY_DIR}/include/confdefs.h"
"#define F77_FUNC(name,NAME) name ## _
#define F77_FUNC_(name,NAME) name ## _"
)
configure_file(include/confdefs.h libvtkfortran/include COPYONLY)
configure_file(include/confdefs.h libadaptivity/include COPYONLY)

file(WRITE "${CMAKE_BINARY_DIR}/include/vtk.h"
"#ifndef VTK_H
#define VTK_H

#ifdef HAVE_VTK
#include <vtkBMPWriter.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkCellDerivatives.h>
#include <vtkCell.h>
#include <vtkClipDataSet.h>
#include <vtkContourGrid.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkDataSetReader.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkGenericCell.h>
#include <vtkHexahedron.h>
#include <vtkIdList.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkMPIController.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkPolyData.h>
#include <vtkShortArray.h>
#include <vtkStructuredGrid.h>
#include <vtkTetra.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkZLibDataCompressor.h>

#ifndef vtkFloatingPointType
#define vtkFloatingPointType vtkFloatingPointType
typedef float vtkFloatingPointType;
#endif

#endif
#endif"
)
configure_file(include/vtk.h libvtkfortran/include COPYONLY)
configure_file(include/vtk.h libadaptivity/include COPYONLY)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
  set(CMAKE_Fortran_FLAGS "-fbounds-check \
-Wall -Wno-surprising -Wno-unused-function -Wno-unused-dummy-argument \
-ffast-math \
-frecord-marker=4 \
-ffree-line-length-none -ffixed-line-length-none \
-fdefault-real-8 -fdefault-double-8 \
${CMAKE_Fortran_FLAGS}"
  )
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-f77rtl -O0 -g")
  set(CMAKE_Fortran_FLAGS "-ffastmath \
-frecord-marker=4 \
-extend_source \
-r8 \
${CMAKE_Fortran_FLAGS}"
  )
else()
  message("CMAKE_Fortran_COMPILER full path: ${CMAKE_Fortran_COMPILER}")
  message("Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
  message("No optimised Fortran compiler flags are known, we just try -O2.")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
endif()

configure_file(tools/version-info.in bin/version-info)
execute_process(COMMAND
  ${CMAKE_BINARY_DIR}/bin/version-info ${CMAKE_BINARY_DIR}/include/version.h
)

# Build libfluidity.a
add_library(fluidity STATIC main/mainfl.F90 main/Fluids.F90 main/Usage.cpp)
add_dependencies(fluidity
  adaptivity mba2d spatialindex spud vtkfortran
  H5hutf Judy
)
target_compile_definitions(fluidity
  PUBLIC
  DOUBLEP
  HAVE_ADAPTIVITY
  HAVE_BLAS
  HAVE_CGAL
  HAVE_EXODUSII
  HAVE_HYPRE
  HAVE_LIBNETCDF
  HAVE_LIBNETCDFF
  HAVE_MEMORY_STATS
  HAVE_MBA_2D
  HAVE_MPI
  HAVE_NUMPY
  HAVE_PETSC
  HAVE_PYTHON
  HAVE_VTK
  HAVE_ZOLTAN
  PACKAGE_NAME=fluidity
  PACKAGE_TARNAME=fluidity
  SIGNAL
  STREAM_IO
)
target_include_directories(fluidity
  PUBLIC
  "${VTK_INCLUDE_DIRS}"
  "${Python3_INCLUDE_DIRS}"
  "${Python3_NumPy_INCLUDE_DIRS}"

  "${PETSC_INCLUDE_DIRS}"

  ${ZOLTAN_INCLUDE_DIR}
  "${NETCDF_INCLUDE_DIRS}"

  ${libspatialindex_SOURCE_DIR}/include
  ${libjudy_SOURCE_DIR}/src

  include
  libwm
)
target_link_libraries(fluidity
  adaptivity
  mba2d
  spatialindex
  spud
  vtkfortran

  H5hutf
  Judy

  "${VTK_LIBRARIES}"
  "${Python3_LIBRARIES}"
  "${CGAL_LIBRARIES}"

  PkgConfig::PETSC
  ${PETSC_STATIC_LDFLAGS}

  ${ZOLTAN_LIBRARY}
  ${NETCDF_STATIC_LDFLAGS}
  ${UDUNITS_LIBRARY}
  ${EXODUSII_LIBRARY}
)
vtk_module_autoinit(TARGETS fluidity MODULES ${VTK_LIBRARIES})

add_subdirectory(assemble)
add_subdirectory(bathymetry)
add_subdirectory(climatology)
add_subdirectory(debug)
add_subdirectory(diagnostics)
add_subdirectory(error_measures)
add_subdirectory(femtools)
add_subdirectory(fldecomp)
add_subdirectory(forward_interfaces)
add_subdirectory(h5hut)
add_subdirectory(horizontal_adaptivity)
add_subdirectory(libadaptivity)
add_subdirectory(libmba2d)
add_subdirectory(libspud)
add_subdirectory(libvtkfortran)
add_subdirectory(libwm)
add_subdirectory(ocean_forcing)
add_subdirectory(parameterisation)
add_subdirectory(population_balance)
add_subdirectory(preprocessor)
add_subdirectory(python)
add_subdirectory(sediments)
# add_subdirectory(tests)
add_subdirectory(tools)

# Build the fluidity binary
add_executable(fluidity-bin main.cpp)
set_target_properties(fluidity-bin PROPERTIES OUTPUT_NAME fluidity)
target_link_libraries(fluidity-bin fluidity)

# # Add unittests
# add_custom_target(unittest
#   "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_BINARY_DIR}/tests/data" "${CMAKE_BINARY_DIR}/bin/tests/data"
#   COMMAND bin/unittestharness bin/tests
# )
# add_dependencies(unittest ${unittests} unittestharness)
# foreach(mesh IN LISTS unittest_meshes_3d)
#   add_dependencies(unittest ${mesh})
# endforeach()
#
# add_custom_target(fltools DEPENDS ${tools} fldecomp)
#
# add_custom_target(shorttest
#   COMMAND bin/testharness -l short
# )
# add_dependencies(shorttest fluidity-bin fltools)
#
# add_custom_target(mediumtest
#   COMMAND ${CMAKE_BINARY_DIR}/bin/testharness -l medium
# )
# add_dependencies(mediumtest fluidity-bin fltools)
