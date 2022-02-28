cmake_minimum_required(VERSION 3.10)
project(vtkfortran C CXX Fortran)

find_package(MPI REQUIRED)
find_package(VTK REQUIRED)

include_directories(./include)
include_directories(${MPI_INCLUDE_PATH} ${VTK_INCLUDE_DIRS})

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

add_library(vtkfortran STATIC vtkfortran.cpp vtkmeshio.cpp fvtkfortran.F90)
target_link_libraries(vtkfortran ${MPI_LIBRARIES} ${VTK_LIBRARIES})


add_definitions("-DDOUBLEP")

set_target_properties(
	vtkfortran
	PROPERTIES
	COMPILE_OPTIONS "${MPI_Fortran_COMPILE_OPTIONS}"
	COMPILE_DEFINITIONS "${MPI_Fortran_COMPILE_DEFINITIONS}"
	LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_VTK")
