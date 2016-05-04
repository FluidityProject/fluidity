set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ../lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)

enable_language (Fortran)

#check for mpi
find_package(MPI REQUIRED)

if(MPI_COMPILE_FLAGS)
	set_target_properties(adaptivity PROPERTIES
      	COMPILE_FLAGS "${MPI_COMPILE_FLAGS} -DHAVE_VTK")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(adaptivity PROPERTIES
      LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

file(GLOB lib_srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} adapt3d/src/*.F adapt3d/src/*.F90 load_balance/src/*.cpp load_balance/src/*.F)

list(REMOVE_ITEM lib_srcs load_balance/src/calculate_submeshes.cpp load_balance/src/metis.cpp load_balance/src/unn_t.cpp load_balance/src/vtkio.cpp)

set(CMAKE_Fortran_MODULE_DIRECTORY include)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -ffixed-line-length-none -fdefault-real-8")

add_definitions(-DHAVE_MPI -DUSING_DOUBLE_PRECISION=1)

add_library(adaptivity STATIC ${lib_srcs})
target_include_directories(adaptivity BEFORE PUBLIC ${MPI_INCLUDE_PATH} include adapt3d/include load_balance/include)
