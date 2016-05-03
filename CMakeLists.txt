cmake_minimum_required(VERSION 2.8)
enable_language (Fortran C CXX)
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
find_package(PETSc REQUIRED)
add_definitions(-DHAVE_PETSC=1
	        -DHAVE_PETSC_MODULES=1)
# seperate check for hypre
find_library(HYBRE_LIB HYPRE)
add_definitions(-DHAVE_HYPRE=1)
#check for zoltan
find_path(ZOLTAN_INCLUDE_DIR zoltan.h)
find_library(ZOLTAN_LIB zoltan)
if(ZOLTAN_LIB)
   message(STATUS "Zoltan library found!")
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
add_definitions(-DHAVE_LIBNETCDF)

#check for spud
find_library(SPUD_LIB spud)
add_definitions(-DHAVE_SPUD=1)
#check for judy
find_library(JUDY_LIB Judy)
add_definitions(-DHAVE_JUDY=1)
#check for udunits
find_library(UDUNITS_LIB udunits2 udunits)
add_definitions(-DHAVE_LIBUDUNITS=1)
#check for spatialindex
find_library(spatialindex spatialindex)

include(${VTK_USE_FILE})
add_definitions(-DHAVE_VTK=1)

find_library (PETSC_LIB petsc HINTS ${PETSc_DIR}/../../)

include_directories(${MPI_INCLUDE_PATH} ${PETSC_PACKAGE_INCLUDES} ${ZOLTAN_INCLUDE_DIR} ./libwm ./libvtkfortran/include ${NETCDF_INCLUDE_DIRS})

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
   message(STATUS "Numpy found!")
   include_directories(${PYTHON_SITE_PACKAGES}/numpy/core/include)
   add_definitions(-DHAVE_NUMPY)
endif()

file(WRITE "${CMAKE_BINARY_DIR}/teststream.F90" 
     "program teststream
      open(unit = 0, access = 'stream')
      end program teststream")
try_compile(RESULT_VAR ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}/teststream.F90 OUTPUT_VARIABLE STREAM)
if (STREAM)
   message(STATUS "Streaming io works")
   add_definitions(-DSTREAM_IO)
else(STREAM)
   message(STATUS "No streaming io")
endif()
file(REMOVE "${CMAKE_BINARY_DIR}/teststream.F90")

include_directories(./include ${MPI_INCLUDE_PATH} ${PETSC_PACKAGE_INCLUDES} ${ZOLTAN_INCLUDE_DIR} ./libwm ./libvtkfortran/include ${PYTHON_INCLUDE_DIRS})


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

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops  -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -frecord-marker=4 -ffree-line-length-none -ffixed-line-length-none -fdefault-real-8 -fdefault-double-8")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

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

execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "copy_if_different" ${CMAKE_SOURCE_DIR}/python ${CMAKE_BINARY_DIR}/python)
execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "copy_directory" ${CMAKE_SOURCE_DIR}/python ${CMAKE_BINARY_DIR}/python)

if(${USING_TWOD_ADAPTIVITY} MATCHES ON)
   message(STATUS "link to libmba2d")
   add_subdirectory(libmba2d)
   set(MBA2D_LIBS mba2d)
   add_definitions(-DHAVE_MBA_2D=1)
else()
   message(STATUS "Not using libmba2d")
   set(USING_TWOD_ADAPTIVITY OFF)
endif()

add_library(fluidity ${LIB_SRCS} main/mainfl.F90 main/Fluids.F90 main/Usage.cpp preprocessor/check_options.F90)
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

target_link_libraries(fluidity-bin fluidity vtkfortran ${spatialindex} adaptivity ${VTK_LIBRARIES} ${MPI_LIBRARIES} ${PETSC_PACKAGE_LIBS} ${PETSC_LIB} ${HYPRE_LIB} ${ZOLTAN_LIB} ${SPUD_LIB} ${JUDY_LIB} ${UDUNITS_LIB} ${MBA2D_LIBS} ${PYTHON_LIBRARIES} ${NETCDF_LIB} ${EXTRA_LIBS})

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
