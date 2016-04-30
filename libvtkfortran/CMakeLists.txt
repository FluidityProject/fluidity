cmake_minimum_required(VERSION 2.6)
enable_language (Fortran CXX)
project(vtkfortran)
find_package(MPI REQUIRED)
find_package(VTK REQUIRED)
include_directories(./include)
include(${VTK_USE_FILE})

include_directories(${MPI_INCLUDE_PATH})

add_library(vtkfortran vtkfortran.cpp vtkmeshio.cpp fvtkfortran.F90)
target_link_libraries(vtkfortran ${MPI_LIBRARIES} ${VTK_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(vtkfortran PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS} -DHAVE_VTK")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(vtkfortran PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_VTK")

