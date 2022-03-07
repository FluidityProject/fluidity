cmake_minimum_required(VERSION 3.14)
project(vtkfortran C CXX Fortran)

find_package(MPI REQUIRED)
find_package(VTK REQUIRED)

set(CMAKE_Fortran_MODULE_DIRECTORY
	${CMAKE_BINARY_DIR}/include
)

add_library(vtkfortran STATIC
	fvtkfortran.F90
	vtkfortran.cpp
	vtkmeshio.cpp
)
target_include_directories(vtkfortran
  PUBLIC
	include
	"${MPI_Fortran_INCLUDE_DIRS}"
	"${VTK_INCLUDE_DIRS}"
)
target_link_libraries(vtkfortran
	"${MPI_Fortran_LIBRARIES}"
  "${VTK_LIBRARIES}"
)

add_definitions(-DDOUBLEP)

set_target_properties(
	vtkfortran
	PROPERTIES
	COMPILE_OPTIONS "${MPI_Fortran_COMPILE_OPTIONS}"
	COMPILE_DEFINITIONS "${MPI_Fortran_COMPILE_DEFINITIONS}"
	LINK_OPTIONS "${MPI_Fortran_LINK_FLAGS}"
)

set(CMAKE_CXX_FLAGS
	"${CMAKE_CXX_FLAGS} -DHAVE_VTK"
)
