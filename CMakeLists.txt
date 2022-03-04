cmake_minimum_required(VERSION 3.14)
project(fluidity
  VERSION 4.1.19
  HOMEPAGE_URL "http://fluidityproject.github.io/"
  LANGUAGES C CXX Fortran
)

include(ExternalProject)
include(FetchContent)

find_package(PkgConfig REQUIRED)

macro(add_library_sources)
  file(RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
  foreach(_src ${ARGN})
    if(_relPath)
      list(APPEND LIB_SRCS "${_relPath}/${_src}")
    else()
      list(APPEND LIB_SRCS "${_src}")
    endif()
  endforeach()
  if(_relPath)
    # propagate LIB_SRCS to parent directory
    set(LIB_SRCS ${LIB_SRCS} PARENT_SCOPE)
  endif()
endmacro()

# Ensure that libraries and binaries end up in the right place
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

# Check for MPI
find_package(MPI REQUIRED)

# Check for VTK
find_package(VTK REQUIRED)

# Check for BLAS
find_package(BLAS REQUIRED)

# Check for HDF5
find_package(HDF5 REQUIRED)

# Check for Python 3
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

# Check for PETSc
pkg_search_module(PETSC PETSc REQUIRED)
message(STATUS "Found PETSc " ${PETSC_VERSION})

# Check for NetCDF
pkg_check_modules(NETCDF netcdf netcdf-fortran REQUIRED)
find_library(NETCDF_LIB netcdf)
find_library(NETCDFF_LIB netcdff)
find_path(NETCDF_INCLUDE_DIR netcdf.inc)

# Check for libspud
find_library(SPUD_LIB spud REQUIRED)
find_path(SPUD_INCLUDE_DIR spud.mod)

# Check for UDUNITS
find_library(UDUNITS_LIB udunits2 REQUIRED)

# # Check for exodusii
# if (USE_EXODUSII)
#    message(STATUS "link to exodusii")
#    find_library(EXODUSII_LIB NAMES exoIIv2c exodus)
#       set(USE_EXODUSII ON CACHE BOOL
#       "include exodusii functionality")
#    add_definitions(-DHAVE_EXODUSII=yes -DHAVE_LIBEXOIIV2C=1)
# else()
#    message(STATUS "No exodusii support")
#    set(USE_EXODUSII OFF CACHE BOOL
#       "Build and link exodusii")
# endif()

# # Check for CGAL
# if (USE_CGAL)
#    message(STATUS "Link in CGAL support")
#    find_package(CGAL REQUIRED)
#    include(${CGAL_USE_FILE})
#    set(USE_CGAL ON CACHE BOOL
#       "Link to CGAL")
#    add_definitions(-DHAVE_CGAL=yes)
#    set(CGAL_LIBS ${CGAL_LIB} )
# else()
#    message(STATUS "No CGAL support")
#    set(USE_CGAL OFF CACHE BOOL
#       "Link to CGAL")
# endif()

# Add Zoltan
ExternalProject_Add(libzoltan
  PREFIX libzoltan
  GIT_REPOSITORY https://github.com/sandialabs/Zoltan
  GIT_TAG origin/main
  CONFIGURE_COMMAND ../libzoltan/configure --prefix=<BINARY_DIR> --enable-f90interface --with-gnumake --with-scotch --with-scotch-incdir=$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include --with-scotch-libdir=$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib --with-parmetis --with-parmetis-incdir=$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include --with-parmetis-libdir=$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
  BUILD_COMMAND make everything
  INSTALL_COMMAND make install
  LOG_DOWNLOAD ON
  LOG_CONFIGURE ON
  LOG_BUILD ON
  LOG_INSTALL ON
)
ExternalProject_Add_Step(libzoltan byproducts
  COMMAND ""
  DEPENDEES install
  BYPRODUCTS <BINARY_DIR>/lib/libzoltan.a
  WORKING_DIRECTORY <BINARY_DIR>
  LOG ON
)
ExternalProject_Get_Property(libzoltan BINARY_DIR)
set(libzoltan_BINARY_DIR ${BINARY_DIR})
add_library(Zoltan STATIC IMPORTED)
set_target_properties(Zoltan
  PROPERTIES IMPORTED_LOCATION ${libzoltan_BINARY_DIR}/lib/libzoltan.a
)

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
  COMMAND ""
  DEPENDEES install
  BYPRODUCTS <SOURCE_DIR>/lib/libJudy.a
  WORKING_DIRECTORY <SOURCE_DIR>
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
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  DEPENDS ${CMAKE_SOURCE_DIR}/tools/make_register_diagnostics.py
  VERBATIM
)
add_custom_command(
  OUTPUT ${CMAKE_SOURCE_DIR}/preprocessor/check_options.F90
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/make_check_options.py
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  DEPENDS ${CMAKE_SOURCE_DIR}/tools/make_check_options.py
  VERBATIM
)
add_custom_command(
  OUTPUT ${CMAKE_SOURCE_DIR}/diagnostics/Diagnostic_Fields_New.F90
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/make_diagnostic_fields.py
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/diagnostics
  DEPENDS ${CMAKE_SOURCE_DIR}/tools/make_diagnostic_fields.py
  VERBATIM
)

# FFLAGS depend on the compiler
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  # gfortran
  set(CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops  -O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
  set(CMAKE_Fortran_FLAGS
  	"${CMAKE_Fortran_FLAGS} -frecord-marker=4 -ffree-line-length-none -ffixed-line-length-none -fdefault-real-8 -fdefault-double-8"
  )
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  # ifort (untested)
  set(CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-f77rtl -O0 -g")
  set(CMAKE_Fortran_FLAGS
  	"${CMAKE_Fortran_FLAGS} -ffastmath -frecord-marker=4 -extend_source -r8"
  )
else (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  message("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message("Fortran compiler: " ${CMAKE_Fortran_COMPILER_ID})
  message("No optimized Fortran compiler flags are known, we just try -O2.")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

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
add_subdirectory(libvtkfortran)
add_subdirectory(libwm)
add_subdirectory(ocean_forcing)
add_subdirectory(parameterisation)
add_subdirectory(population_balance)
add_subdirectory(preprocessor)
add_subdirectory(sediments)
add_subdirectory(tests)
add_subdirectory(tools)

add_definitions(
  -DPACKAGE_NAME="fluidity"
  -DPACKAGE_TARNAME="fluidity"
  -DDOUBLEP
  -DHAVE_ADAPTIVITY
  -DHAVE_MEMORY_STATS=1
  -DHAVE_MPI
  -DHAVE_NUMPY
  -DSIGNAL
)

# Build libfluidity.a
add_library(fluidity STATIC
  "${LIB_SRCS}"
  main/mainfl.F90
  main/Fluids.F90
  main/Usage.cpp
)

target_include_directories(fluidity
  PUBLIC
  "${MPI_Fortran_INCLUDE_DIRS}"
  "${VTK_INCLUDE_DIRS}"
  "${BLAS_LINKER_FLAGS}"
  "${HDF5_Fortran_INCLUDE_DIRS}"
  "${Python3_INCLUDE_DIRS}"

  "${PETSC_INCLUDE_DIRS}"

  ${NETCDF_INCLUDE_DIR}
  ${SPUD_INCLUDE_DIR}

  ${libspatialindex_SOURCE_DIR}/include

  ${CMAKE_SOURCE_DIR}/include

  ${CMAKE_BINARY_DIR}/libadaptivity/include
  ${CMAKE_BINARY_DIR}/libvtkfortran/include

  ${CMAKE_SOURCE_DIR}/h5hut/src/include
  ${CMAKE_SOURCE_DIR}/libwm
)
target_link_libraries(fluidity
  Zoltan

  "${MPI_Fortran_LIBRARIES}"
  "${VTK_LIBRARIES}"
  "${BLAS_LIBRARIES}"
  "${HDF5_Fortran_LIBRARIES}"

  "${PETSC_LINK_LIBRARIES}"
)
add_dependencies(fluidity
  adaptivity
  mba2d
  vtkfortran

  H5hutf
  Zoltan
)

# Build the fluidity binary
add_executable(fluidity-bin main.cpp)
set_target_properties(fluidity-bin
  PROPERTIES
  OUTPUT_NAME fluidity
)
target_link_libraries(fluidity-bin
  adaptivity
  fluidity
  mba2d
  spatialindex
  vtkfortran

  H5hutf
  Judy
  Zoltan

  "${MPI_Fortran_LIBRARIES}"
  "${VTK_LIBRARIES}"
  "${BLAS_LIBRARIES}"
  "${HDF5_Fortran_LIBRARIES}"
  "${Python3_LIBRARIES}"

  -L${PETSC_LIBRARY_DIRS}

  ${NETCDF_LIB}
  ${NETCDFF_LIB}

  ${SPUD_LIB}
  ${UDUNITS_LIB}
)

# Execute unittests
add_custom_target(unittest
  "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_BINARY_DIR}/tests/data" "${CMAKE_BINARY_DIR}/bin/tests/data"
  COMMAND bin/unittestharness bin/tests
)
add_dependencies(unittest ${unittests} unittestharness)
foreach(mesh IN LISTS unittest_meshes_3d)
  add_dependencies(unittest ${mesh})
endforeach()

add_custom_target(fltools DEPENDS ${tools} fldecomp)
file(
  WRITE "${CMAKE_BINARY_DIR}/include/confdefs.h"
"#define F77_FUNC(name,NAME) name ## _
#define F77_FUNC_(name,NAME) name ## _"
)
file(
  WRITE "${CMAKE_BINARY_DIR}/libvtkfortran/include/confdefs.h"
"#define F77_FUNC(name,NAME) name ## _
#define F77_FUNC_(name,NAME) name ## _"
)
file(
  WRITE "${CMAKE_BINARY_DIR}/include/vtk.h"
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

file(
  COPY "${CMAKE_BINARY_DIR}/include/vtk.h"
  DESTINATION "${CMAKE_BINARY_DIR}/libvtkfortran/include/"
)

configure_file("tools/version-info.in" "bin/version-info")

execute_process(
  COMMAND
  "${CMAKE_BINARY_DIR}/bin/version-info" "${CMAKE_BINARY_DIR}/include/version.h"
)

add_custom_target(shorttest
  COMMAND bin/testharness -l short
)
add_dependencies(shorttest fluidity-bin fltools)

add_custom_target(mediumtest
  COMMAND ${CMAKE_BINARY_DIR}/bin/testharness -l medium
)
add_dependencies(mediumtest fluidity-bin fltools)
