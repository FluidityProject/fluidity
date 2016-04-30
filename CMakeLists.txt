cmake_minimum_required(VERSION 2.6)
enable_language (Fortran CXX)
project(fluidity)

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")

#check for mpi
find_package(MPI REQUIRED)
#check for vtk
find_package(VTK REQUIRED)
# check for PETSc
find_package(PETSc REQUIRED)
add_definitions(-DHAVE_PETSC=1
	        -DHAVE_PETSC_MODULES=1)
#check for zoltan
find_path(ZOLTAN_INCLUDE_DIR zoltan.h)
find_library(ZOLTAN_LIB zoltan)
if(ZOLTAN_LIB)
   message(STATUS "Zoltan library found!")
   add_definitions(-DHAVE_ZOLTAN)
endif()
#check for spud
find_library(SPUD_LIB spud)
add_definitions(-DHAVE_SPUD=1)
#check for judy
find_library(JUDY_LIB Judy)
add_definitions(-DHAVE_JUDY=1)
#check for udunits
find_library(UDUNITS_LIB udunits)
add_definitions(-DHAVE_UDUNITS=1)
#check for spatialindex
find_library(spatialindex spatialindex)

include(${VTK_USE_FILE})
add_definitions(-DHAVE_VTK=1)

find_library (PETSC_LIB petsc HINTS ${PETSc_DIR}/../../)

include_directories(${MPI_INCLUDE_PATH} ${PETSC_PACKAGE_INCLUDES} ${ZOLTAN_INCLUDE_DIR} ./libwm ./libvtkfortran/include)

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

include_directories(./include ${MPI_INCLUDE_PATH} ${PETSC_PACKAGE_INCLUDES} ${ZOLTAN_INCLUDE_DIR} ./libwm ./libvtkfortran/include)


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
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -frecord-marker=4 -ffree-line-length-none -ffixed-line-length-none -fdefault-real-8 -fdefault-double-8")

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

if(MBA2D)
   add_subdirectory(libmba2d)
   set(MBA2D_LIBS mba2d)
   add_definitions(-DHAVE_MBA_2D=1)
endif()

add_library(fluidity ${LIB_SRCS})
add_executable(fluidity-bin main.cpp main/mainfl.F90 main/Fluids.F90 main/Usage.cpp preprocessor/check_options.F90)

SET_TARGET_PROPERTIES(fluidity-bin PROPERTIES OUTPUT_NAME fluidity)

if(MPI_COMPILE_FLAGS)
	set_target_properties(fluidity PROPERTIES
      	COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(fluidity PROPERTIES
      LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

target_link_libraries(fluidity-bin fluidity vtkfortran ${spatialindex} adaptivity ${VTK_LIBRARIES} ${MPI_LIBRARIES} ${PETSC_PACKAGE_LIBS} ${PETSC_LIB} ${ZOLTAN_LIB} ${SPUD_LIB} ${JUDY_LIB} ${UDUNITS_LIB} ${MBA2D_LIBS} ${PYTHON_LIBRARIES})

add_definitions(-DPACKAGE_NAME="fluidity"
	        -DPACKAGE_TARNAME="fluidity"
		-DSIGNAL
		-DHAVE_MEMORY_STATS=1)
