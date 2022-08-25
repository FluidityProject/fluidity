find_package(SWIG REQUIRED)
find_package(Python COMPONENTS Development NumPy REQUIRED)
include(UseSWIG)

# add dependency (not sure if CMake pick this up?)
set_property(SOURCE H5hut.i PROPERTY DEPENDS numpy.i)

# define SWIG library
swig_add_library(H5hutpy TYPE SHARED
  LANGUAGE Python
  SOURCES H5hut.i)
swig_link_libraries(H5hutpy H5hut)

set_property(TARGET H5hutpy PROPERTY SWIG_GENERATED_INCLUDE_DIRECTORIES ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})
get_property(support_files TARGET H5hutpy PROPERTY SWIG_SUPPORT_FILES)

install(TARGETS H5hutpy
  DESTINATION python)
install(FILES ${support_files}
  DESTINATION python)
