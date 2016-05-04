cmake_minimum_required(VERSION 2.8)
enable_language (Fortran)
project(fluidity)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")

#check for mpi
find_package(MPI REQUIRED)
#check for vtk
find_package(VTK REQUIRED)
# check for PETSc
set(PETSC_DIR "/usr/lib/petscdir/3.6.3" CACHE STRING
  "Base directory for the PETsc installation")
find_program(petscarch petscarch)
if(petscarch)
  execute_process(COMMAND petscarch
    OUTPUT_VARIABLE PETSC_ARCH)
endif()
set(PETSC_ARCH ${PETSC_ARCH} CACHE STRING
  "Architecture for the PETSc installation")
find_package(PETSc REQUIRED
  HINTS ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf
    /usr/lib/petsc/conf)
find_library (PETSC_LIB petsc HINTS ${PETSc_DIR}/../../)
add_definitions(-DHAVE_PETSC=1
	        -DHAVE_PETSC_MODULES=1)
get_filename_component(PETSC_LIB_DIR ${PETSC_LIB} DIRECTORY)
find_path(PETSC_EXTRAS_DIR petscversion.h
          HINTS ${PETSC_LIB_DIR}/../../include)
# seperate check for hypre
find_library(HYPRE_LIB HYPRE
             HINTS "${PETSC_LIB_DIR}")
message(STATUS "using PETSc:${PETSC_LIB}")
message(STATUS "using HYPRE:${HYPRE_LIB}")
add_definitions(-DHAVE_HYPRE=1)
#check for zoltan
find_path(ZOLTAN_INCLUDE_DIR zoltan.h
  HINTS ${PETSC_PACKAGE_INCLUDES} ${PETSC_EXTRAS_DIR})
find_library(ZOLTAN_LIB zoltan 
  HINTS ${PETSC_LIB_DIR})
if(ZOLTAN_LIB)
   message(STATUS "Using Zoltan library:" ${ZOLTAN_LIB})
   add_definitions(-DHAVE_ZOLTAN)
endif()
# netcdf searching
find_package(PkgConfig)
pkg_check_modules(NETCDF QUIET netcdf)
pkg_check_modules(NETCDFF QUIET netcdf-fortran)
find_library(NETCDF_LIB netcdf
             HINTS ${NETCDF_LIBDIR} ${NETCDF_LIBRARY_DIRS} )
find_library(NETCDFF_LIB netcdff
             HINTS ${NETCDFF_LIBDIR} ${NETCDFF_LIBRARY_DIRS} )
if (NETCDF_LIB)
  add_definitions(-DHAVE_LIBNETCDF)
endif()

#check for spud
find_library(SPUD_LIB spud)
if(NOT SPUD_LIB)
  add_subdirectory(libspud)
  set(SPUD_LIB spud CACHE STRING
      "libspud library to link fluidity against")
endif()
get_filename_component(SPUD_LIB_DIR ${SPUD_LIB} DIRECTORY)
find_path(SPUD_INCLUDE_DIR spud.mod
          HINTS ${SPUD_LIB_DIR} ${SPUD_LIB_DIR}/../include)
add_definitions(-DHAVE_SPUD=1)
#check for judy
find_library(JUDY_LIB Judy)
if(NOT JUDY_LIB)
  add_subdirectory(libjudy)
  set(JUDY_LIB judy CACHE STRING
      "libjudy library to link fluidity against")
endif()
add_definitions(-DHAVE_JUDY=1)
#check for udunits
find_library(UDUNITS_LIB udunits2 udunits)
add_definitions(-DHAVE_LIBUDUNITS=1)
#check for spatialindex
find_library(spatialindex spatialindex)
if(NOT spatialindex)
   add_subdirectory(spatialindex-1.8.0)
   set(spatialindex spatialindex CACHE STRING
     "Copy of spatial index to link against")
endif()

include(${VTK_USE_FILE})
add_definitions(-DHAVE_VTK=1)

include_directories(${MPI_INCLUDE_PATH} ${PETSC_PACKAGE_INCLUDES} ${ZOLTAN_INCLUDE_DIR} ./libwm ./libvtkfortran/include ${NETCDF_INCLUDE_DIRS} ${PETSC_EXTRAS_DIR} ${SPUD_INCLUDE_DIR} ${CMAKE_BINARY_DIR}/include)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(SIDX_INCLUDE_DIR "../include")
set(SIDX_BIN_DIR "../lib")
set(SIDX_LIB_DIR "../lib")
#add_subdirectory(spatialindex-1.8.0 EXCLUDE_FROM_ALL)
add_subdirectory(libvtkfortran)
add_subdirectory(libadaptivity)
add_definitions(-DHAVE_ADAPTIVITY)

find_package(PythonLibs REQUIRED)
add_definitions(-DHAVE_PYTHON)
execute_process ( COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib(True)" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)

find_path(PYTHON_NUMPY_FOUND numpy/version.py HINTS ${PYTHON_SITE_PACKAGES})
if(PYTHON_NUMPY_FOUND)
   message(STATUS "Found Numpy")
   include_directories(${PYTHON_SITE_PACKAGES}/numpy/core/include)
   add_definitions(-DHAVE_NUMPY)
endif()

file(WRITE "${CMAKE_BINARY_DIR}/teststream.F90" 
     "program teststream
      open(unit = 0, access = 'stream')
      end program teststream")
try_compile(RESULT_VAR ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}/teststream.F90 OUTPUT_VARIABLE STREAM)
if (STREAM)
   message(STATUS "Streaming IO works")
   add_definitions(-DSTREAM_IO)
else(STREAM)
   message(STATUS "No streaming IO")
endif()
file(REMOVE "${CMAKE_BINARY_DIR}/teststream.F90")

include_directories(./include ${MPI_INCLUDE_PATH} ${PETSC_PACKAGE_INCLUDES} ${ZOLTAN_INCLUDE_DIR} ./libwm ./libvtkfortran/include ${PYTHON_INCLUDE_DIRS})

execute_process(COMMAND ${CMAKE_COMMAND} "-E" chdir "${CMAKE_SOURCE_DIR}" python tools/make_register_diagnostics.py)
execute_process(COMMAND ${CMAKE_COMMAND} "-E" chdir "${CMAKE_SOURCE_DIR}" python tools/make_check_options.py)

macro (add_Library_sources)
    file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND LIB_SRCS "${_relPath}/${_src}")
        else()
            list (APPEND LIB_SRCS "${_src}")
        endif()
    endforeach()
    if (_relPath)
        # propagate LIB_SRCS to parent directory
        set (LIB_SRCS ${LIB_SRCS} PARENT_SCOPE)
    endif()
endmacro()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops  -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -frecord-marker=4 -ffree-line-length-none -ffixed-line-length-none -fdefault-real-8 -fdefault-double-8")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffastmath -frecord-marker=4 -extend_source -r8")
else (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${CMAKE_Fortran_COMPILER_ID})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

add_library_sources()
add_subdirectory(debug)
add_subdirectory(libwm)
add_subdirectory(femtools)
add_subdirectory(bathymetry)
add_subdirectory(ocean_forcing)
add_subdirectory(sediments)
add_subdirectory(parameterisation)
add_subdirectory(forward_interfaces)
add_subdirectory(horizontal_adaptivity)
add_subdirectory(preprocessor)
add_subdirectory(error_measures)
add_subdirectory(assemble)
add_subdirectory(diagnostics)
add_subdirectory(reduced_modelling)
add_subdirectory(tests)
add_subdirectory(fldecomp)
add_subdirectory(climatology)


#execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "copy_directory" ${CMAKE_SOURCE_DIR}/python ${CMAKE_BINARY_DIR}/python)

if(${USE_LIBMBA} MATCHES ON)
   message(STATUS "link to libmba2d")
   add_subdirectory(libmba2d)
   set(MBA2D_LIBS mba2d)
   add_definitions(-DHAVE_MBA_2D=1)
else()
   message(STATUS "Not using libmba2d")
   set(USE_LIBMBA OFF CACHE STRING
      "Build and link libmba2d")
endif()

add_library(fluidity STATIC ${LIB_SRCS} main/mainfl.F90 main/Fluids.F90 main/Usage.cpp preprocessor/check_options.F90)
add_dependencies(fluidity vtkfortran)
add_executable(fluidity-bin main.cpp)

SET_TARGET_PROPERTIES(fluidity-bin PROPERTIES OUTPUT_NAME fluidity)

if(MPI_COMPILE_FLAGS)
	set_target_properties(fluidity PROPERTIES
      	COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(fluidity PROPERTIES
      LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

if(APPLE)
   find_library(COCOA_LIBRARY Cocoa)
   find_library(EXPAT_LIB expat)
   MARK_AS_ADVANCED (COCOA_LIBRARY EXPAT_LIB)
   set(EXTRA_LIBS ${COCOA_LIBRARY} ${EXPAT_LIB})
endif()

target_link_libraries(fluidity-bin fluidity vtkfortran ${spatialindex} adaptivity ${VTK_LIBRARIES} ${MPI_LIBRARIES} ${ZOLTAN_LIB} ${PETSC_PACKAGE_LIBS} ${PETSC_LIB} ${HYPRE_LIB} ${SPUD_LIB} ${JUDY_LIB} ${UDUNITS_LIB} ${MBA2D_LIBS} ${PYTHON_LIBRARIES} ${NETCDF_LIB} ${NETCDFF_LIB} ${EXTRA_LIBS})

add_definitions(-DPACKAGE_NAME="fluidity"
	        -DPACKAGE_TARNAME="fluidity"
		-DSIGNAL
                -DHAVE_MPI
                -DDOUBLEP
		-DHAVE_MEMORY_STATS=1)

add_subdirectory(tools)

add_custom_target(unittest
        "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_SOURCE_DIR}/tests/data"  "${CMAKE_BINARY_DIR}/bin/tests/data"
	COMMAND bin/unittestharness bin/tests)
add_dependencies(unittest ${unittests})
foreach(mesh IN LISTS unittest_meshes_3d)
    add_dependencies(unittest ${mesh})
endforeach()

add_custom_target(fltools DEPENDS ${tools} fldecomp)
file(WRITE "${CMAKE_SOURCE_DIR}/include/confdefs.h"
"#define F77_FUNC(name,NAME) name ## _
#define F77_FUNC_(name,NAME) name ## _")

file(WRITE "${CMAKE_SOURCE_DIR}/libvtkfortran/include/confdefs.h"
"#define F77_FUNC(name,NAME) name ## _
#define F77_FUNC_(name,NAME) name ## _")

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
#endif")

file(COPY "${CMAKE_BINARY_DIR}/include/vtk.h" DESTINATION "${CMAKE_BINARY_DIR}/libvtkfortran/include/")

configure_file("tools/version-info.in" "bin/version-info")

execute_process(COMMAND "${CMAKE_BINARY_DIR}/bin/version-info" "${CMAKE_BINARY_DIR}/include/version.h")

add_custom_target(test  
  COMMAND ${CMAKE_COMMAND} -"-E" "chdir" "${CMAKE_BINARY_DIR}/tests" "../bin/testharness -l short"
  DEPENDS fluidity-bin fltools)

add_custom_target(mediumtest 
  COMMAND ${CMAKE_COMMAND} -"-E" "chdir" "${CMAKE_BINARY_DIR}/tests" "../bin/testharness -l medium"
  DEPENDS fluidity-bin fltools)