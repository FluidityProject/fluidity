cmake_minimum_required(VERSION 3.19)
project(fluidity
  VERSION 4.1.19
  HOMEPAGE_URL "http://fluidityproject.github.io/"
  LANGUAGES C CXX Fortran
)

function(build_unittest_executable unittest_sources)
  foreach(unittest_source ${unittest_sources})
    get_filename_component(unittest ${unittest_source} NAME_WE)

    add_executable(${unittest} ${unittest_source} test_main.cpp)
    set_target_properties(${unittest}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests
    )
    target_compile_definitions(${unittest} PUBLIC TESTNAME=${unittest}_)
    target_link_libraries(${unittest} fluidity)
  endforeach()
endfunction()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Set build type; \
typical values include Debug, Release, RelWithDebInfo and MinSizeRel." FORCE
  )
elseif(CMAKE_BUILD_TYPE STREQUAL Debug)
  set(DDEBUG DDEBUG)
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

set(CMAKE_CXX_STANDARD 17)

option(HDF5_PREFER_PARALLEL "Look for parallel HDF5 first" ON)

include(ExternalProject)
include(FetchContent)

find_program(SED NAMES sed gsed REQUIRED)

find_package(PkgConfig REQUIRED)

# Check for MPI
find_package(MPI REQUIRED C CXX Fortran)

# Check for HDF5
find_package(HDF5 REQUIRED Fortran)

# Check for VTK
find_package(VTK REQUIRED
  CommonCore
  CommonDataModel
  FiltersCore
  FiltersGeneral
  FiltersVerdict
  IOLegacy
  IOCore
  IOImage
  IOParallelXML
  IOXML
  ParallelMPI
)

# Check for Python 3
find_package(Python3 REQUIRED Interpreter Development NumPy)
option(PREPEND_FLUIDITY_PYTHONPATH "Whether to prepend CMAKE_SOURCE_PATH to PYTHONPATH in all commands" ON)
if(PREPEND_FLUIDITY_PYTHONPATH)
  set(FLUIDITY_PYTHON_PATH "${CMAKE_SOURCE_DIR}/python")
  configure_file(tools/prepend-python-path.in bin/prepend-python-path @ONLY)
  set(PYENV prepend-python-path)
  Python3_EXECUTABLE = "prepend-python-path ${PYTHON3_EXECUTABLE}"
endif()

# # Check for OpenMP
# find_package(OpenMP QUIET COMPONENTS Fortran)

# Check for CGAL
find_package(CGAL QUIET)
if(CGAL_FOUND)
  message(STATUS "Found CGAL ${CGAL_VERSION}")
endif()

# Check for PETSc
if(DEFINED ENV{PETSC_ARCH})
  set(CMAKE_PREFIX_PATH $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
else()
  set(CMAKE_PREFIX_PATH $ENV{PETSC_DIR})
endif()
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
message(STATUS "Found PETSc ${PETSC_VERSION}")

# # Check for ARPACK
# pkg_search_module(ARPACK REQUIRED arpack)
# message(STATUS "Found ARPACK ${ARPACK_VERSION}")

# Check for NetCDF
pkg_search_module(NETCDF REQUIRED netcdf)
message(STATUS "Found NetCDF ${NETCDF_VERSION}")

# Check for TCMalloc
pkg_search_module(TCMALLOC QUIET libtcmalloc)
if(TCMALLOC_FOUND)
  message(STATUS "Found TCMalloc ${TCMALLOC_VERSION}")
endif()

# Check for UDUNITS
pkg_search_module(UDUNITS QUIET udunits)
if(UDUNITS_FOUND)
  message(STATUS "Found UDUNITS ${UDUNITS_VERSION}")
endif()

# Check for libsupermesh
pkg_search_module(LIBSUPERMESH QUIET libsupermesh)
if(LIBSUPERMESH_FOUND)
  message(STATUS "Found libsupermesh ${LIBSUPERMESH_VERSION}")
endif()

# Check for ExodusII
find_library(EXODUSII_LIBRARY exoIIv2c)
find_path(EXODUSII_INCLUDE_DIR exodusII.h)
if(EXODUSII_LIBRARY)
  message(STATUS "Found ExodusII")
endif()

find_library(SPATIALINDEX_LIBRARY spatialindex REQUIRED)
find_path(SPATIALINDEX_INCLUDE_DIR SpatialIndex.h PATH_SUFFIXES spatialindex REQUIRED)
if(SPATIALINDEX_LIBRARY)
  message(STATUS "Found SpatialIndex")
endif()

find_library(JUDY_LIBRARY Judy REQUIRED)
find_path(JUDY_INCLUDE_DIR Judy.h REQUIRED)
if(JUDY_LIBRARY)
  message(STATUS "Found Judy")
endif()

find_library(PARMETIS_LIBRARY parmetis REQUIRED)
find_path(PARMETIS_INCLUDE_DIR parmetis.h REQUIRED)
if(PARMETIS_LIBRARY)
  message(STATUS "Found Parmetis")
endif()
find_library(PARMETIS_LIBRARY metis REQUIRED)
# Add Zoltan
find_package(Zoltan REQUIRED)



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
configure_file("${CMAKE_BINARY_DIR}/include/confdefs.h"
  libvtkfortran/include/confdefs.h
  COPYONLY
)
configure_file("${CMAKE_BINARY_DIR}/include/confdefs.h"
  libadaptivity/include/confdefs.h
  COPYONLY
)

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
#include <vtkMPICommunicator.h>
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
configure_file("${CMAKE_BINARY_DIR}/include/vtk.h"
  libvtkfortran/include/vtk.h
  COPYONLY
)
configure_file("${CMAKE_BINARY_DIR}/include/vtk.h"
  libadaptivity/include/vtk.h
  COPYONLY
)

set(CMAKE_C_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS_DEBUG "-O0")
set(CMAKE_C_FLAGS "-fsanitize=bounds \
-fno-unsafe-math-optimizations -frounding-math -fsignaling-nans \
-g"
)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-O0")
set(CMAKE_CXX_FLAGS "-fsanitize=bounds \
-fno-unsafe-math-optimizations -frounding-math -fsignaling-nans \
-g"
)
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0")
set(CMAKE_Fortran_FLAGS "-fcheck=bounds \
-fdefault-real-8 -fdefault-double-8 \
-ffixed-line-length-none -ffree-line-length-none \
-frecord-marker=4 \
-fno-unsafe-math-optimizations -frounding-math -fsignaling-nans \
-g"
)
if(${CMAKE_C_COMPILER_VERSION} GREATER_EQUAL 10)
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch ${CMAKE_Fortran_FLAGS}")
endif()
if(DEFINED GITHUB_ACTIONS)
  string(APPEND CMAKE_C_FLAGS " -march=sandybridge")
  string(APPEND CMAKE_CXX_FLAGS " -march=sandybridge")
  string(APPEND CMAKE_Fortran_FLAGS " -march=sandybridge")
else()
  string(APPEND CMAKE_C_FLAGS " -march=native")
  string(APPEND CMAKE_CXX_FLAGS " -march=native")
  string(APPEND CMAKE_Fortran_FLAGS " -march=native")
endif()

configure_file(tools/version-info.in bin/version-info @ONLY)
execute_process(COMMAND bin/version-info include/version.h OUTPUT_QUIET)

# Build the fluidity static library
add_library(fluidity STATIC $<TARGET_OBJECTS:spud>)
target_compile_definitions(fluidity PUBLIC
  ${DDEBUG}
  DOUBLEP
  HAVE_ADAPTIVITY
  HAVE_BLAS
  # HAVE_HYPERLIGHT
  HAVE_HYPRE
  # HAVE_LIBARPACK
  $<$<BOOL:${CGAL_FOUND}>:HAVE_LIBCGAL>
  $<$<BOOL:${EXODUSII_LIBRARY}>:HAVE_LIBEXOIIV2C>
  HAVE_LIBNETCDF
  $<$<BOOL:${LIBSUPERMESH_FOUND}>:HAVE_LIBSUPERMESH>
  $<$<BOOL:${TCMALLOC_FOUND}>:HAVE_LIBTCMALLOC>
  $<$<BOOL:${UDUNITS_FOUND}>:HAVE_LIBUDUNITS>
  HAVE_MEMORY_STATS
  HAVE_MBA_2D
  # HAVE_MBA_3D
  HAVE_MPI
  HAVE_NUMPY
  HAVE_PETSC
  HAVE_PYTHON
  HAVE_VTK
  HAVE_ZOLTAN
  # IO_ADVANCE_BUG
  # PACKAGE_NAME=fluidity
  # PACKAGE_TARNAME=fluidity
  SED=${SED}
  SIGNAL
  # SIGNAL_HAVE_FLAG
  STREAM_IO
  # USING_FEMDEM
  # _OPENMP
)
target_include_directories(fluidity PUBLIC
  ${HDF5_Fortran_INCLUDE_DIRS}
  ${Python3_INCLUDE_DIRS}
  ${Python3_NumPy_INCLUDE_DIRS}
  # ${OpenMP_Fortran_INCLUDE_DIRS}

  ${PETSC_INCLUDE_DIRS}
  ${NETCDF_INCLUDE_DIRS}

  $<$<BOOL:${EXODUSII_INCLUDE_DIR}>:${EXODUSII_INCLUDE_DIR}>

  ${libspatialindex_SOURCE_DIR}/include
  ${LIBJUDY_INCLUDE_DIR}
  ${ZOLTAN_INCLUDE_DIR}

  ${CMAKE_BINARY_DIR}/include
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/libspud/include
  ${CMAKE_SOURCE_DIR}/libwm
)
target_link_libraries(fluidity
  adaptivity mba2d spatialindex vtkfortran
  H5hutf Judy

  ${ZOLTAN_LIBRARIES}
  ${HDF5_Fortran_LIBRARIES}
  ${VTK_LIBRARIES}
  ${Python3_LIBRARIES}
  # ${OpenMP_Fortran_LIBRARIES}
  $<TARGET_NAME_IF_EXISTS:CGAL::CGAL>

  PkgConfig::PETSC
  ${PETSC_STATIC_LDFLAGS}
  # ${ARPACK_STATIC_LDFLAGS}
  ${NETCDF_STATIC_LDFLAGS}
  $<$<BOOL:${TCMALLOC_STATIC_LDFLAGS}>:${TCMALLOC_STATIC_LDFLAGS}>
  $<$<BOOL:${UDUNITS_STATIC_LDFLAGS}>:${UDUNITS_STATIC_LDFLAGS}>
  $<$<BOOL:${LIBSUPERMESH_STATIC_LDFLAGS}>:${LIBSUPERMESH_STATIC_LDFLAGS}>

  $<$<BOOL:${EXODUSII_LIBRARY}>:${EXODUSII_LIBRARY}>
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
add_subdirectory(main)
add_subdirectory(ocean_forcing)
add_subdirectory(parameterisation)
add_subdirectory(population_balance)
add_subdirectory(preprocessor)
add_subdirectory(python)
add_subdirectory(sediments)
add_subdirectory(tests)
add_subdirectory(tools)

# Build the fluidity binary
add_executable(fluidity_binary main.cpp)
set_target_properties(fluidity_binary PROPERTIES OUTPUT_NAME fluidity)
target_link_libraries(fluidity_binary fluidity)

# Install Python packages included in Fluidity's repository
add_custom_target(install_python_packages
  COMMAND ${CMAKE_MAKE_PROGRAM} fluidity_python spud_diamond spud_dxdiff spud_python
  VERBATIM
)

# Uninstall Python packages included in Fluidity's repository
add_custom_target(uninstall_python_packages
  COMMAND ${Python3_EXECUTABLE} -m pip uninstall -y diamond dxdiff fluidity libspud
  VERBATIM
)

# Build the fluidity manual
add_custom_target(manual
  COMMAND BUILD_TOOL=cmake make
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/manual
  VERBATIM
)

# Add unittests
execute_process(
  COMMAND ${CMAKE_COMMAND} -E make_directory bin/tests
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
execute_process(
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/tests/data bin/tests/data
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
add_custom_target(unittest COMMAND bin/unittestharness -x unit-tests.xml VERBATIM)

# Add tests
add_custom_target(vshorttest
  COMMAND ${Python3_EXECUTABLE} bin/testharness -x vshort-tests.xml -l vshort
  VERBATIM
)
add_custom_target(shorttest
  COMMAND ${Python3_EXECUTABLE} bin/testharness -x short-tests.xml -l short
  VERBATIM
)
add_custom_target(mediumtest
  COMMAND ${Python3_EXECUTABLE} bin/testharness -x medium-tests.xml -l medium
  VERBATIM
)
add_custom_target(longtest
  COMMAND ${Python3_EXECUTABLE} bin/testharness -x long-tests.xml -l long
  VERBATIM
)
add_custom_target(vlongtest
  COMMAND ${Python3_EXECUTABLE} bin/testharness -x vlong-tests.xml -l vlong
  VERBATIM
)
