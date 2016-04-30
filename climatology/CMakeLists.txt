add_executable(create_climatology_atlas EXCLUDE_FROM_ALL create_climatology_atlas.cpp)
target_link_libraries(create_climatology_atlas fluidity vtkfortran ${spatialindex} adaptivity ${VTK_LIBRARIES} ${MPI_LIBRARIES} ${PETSC_PACKAGE_LIBS} ${PETSC_LIB} ${ZOLTAN_LIB} ${SPUD_LIB} ${JUDY_LIB} ${UDUNITS_LIB} ${MBA2D_LIBS} ${PYTHON_LIBRARIES})
