set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ../lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)

enable_language (Fortran)

# Check for MPI
find_package(MPI REQUIRED)
file(GLOB lib_srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} adapt3d/src/*.F adapt3d/src/*.F90 load_balance/src/*.cpp load_balance/src/*.F)

list(REMOVE_ITEM lib_srcs load_balance/src/calculate_submeshes.cpp load_balance/src/metis.cpp load_balance/src/unn_t.cpp load_balance/src/vtkio.cpp)

set(CMAKE_Fortran_MODULE_DIRECTORY include)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -ffixed-line-length-none -fdefault-real-8")

add_definitions(-DHAVE_MPI -DUSING_DOUBLE_PRECISION=1)

add_library(adaptivity STATIC ${lib_srcs})
target_include_directories(adaptivity BEFORE PUBLIC ${MPI_Fortran_INCLUDE_DIRS} include adapt3d/include load_balance/include)

set_target_properties(
	adaptivity
	PROPERTIES
	COMPILE_OPTIONS "${MPI_Fortran_COMPILE_OPTIONS}"
	COMPILE_DEFINITIONS "${MPI_Fortran_COMPILE_DEFINITIONS}"
	LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
)
