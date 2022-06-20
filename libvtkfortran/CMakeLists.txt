cmake_minimum_required(VERSION 3.21)
project(vtkfortran LANGUAGES CXX Fortran)

add_library(vtkfortran STATIC fvtkfortran.F90 vtkfortran.cpp vtkmeshio.cpp)
target_compile_definitions(vtkfortran PUBLIC DOUBLEP HAVE_VTK)
target_include_directories(vtkfortran PUBLIC
  include ${CMAKE_CURRENT_BINARY_DIR}/include
)
target_link_libraries(vtkfortran "${VTK_LIBRARIES}")
