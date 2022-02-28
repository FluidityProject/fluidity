add_executable(fldecomp EXCLUDE_FROM_ALL fldmain.cpp fldgmsh.cpp fldtriangle.cpp partition.cpp)

target_link_libraries(fldecomp fluidity vtkfortran ${spatialindex} adaptivity ${VTK_LIBRARIES} ${MPI_LIBRARIES} -L${PETSC_LIB_DIR} ${ZOLTAN_LIB} ${PETSC_PACKAGE_LIBS} ${PETSC_LIB} ${SPUD_LIB})

add_definitions(-DHAVE_MPI)
