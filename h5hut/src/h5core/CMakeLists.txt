configure_file(config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

add_library(H5hut
  h5_attachments.c
  h5_attribs.c
  h5b_attribs.c
  h5b_io.c
  h5b_model.c
  h5_err.c
  h5_file.c
  h5_init.c
  h5_log.c
  h5_model.c
  h5_syscall.c
  h5t_adjacencies.c
  h5t_io.c
  h5t_map.c
  h5t_model.c
  h5t_octree.c
  h5t_retrieve.c
  h5t_store.c
  h5t_tags.c
  h5u_io.c
  h5u_model.c
  private/h5_fcmp.c
  private/h5_hdf5.c
  private/h5_hsearch.c
  private/h5_io.c
  private/h5_lustre.c
  private/h5_maps.c
  private/h5_qsort.c
  private/h5_qsort_r.c
  private/h5t_access_tetm.c
  private/h5t_access_trim.c
  private/h5t_adjacencies_tetm.c
  private/h5t_adjacencies_trim.c
  private/h5t_core.c
  private/h5t_core_tetm.c
  private/h5t_core_trim.c
  private/h5t_io_tetm.c
  private/h5t_io_trim.c
  private/h5t_model_tetm.c
  private/h5t_model_trim.c
  private/h5t_ref_elements.c
  private/h5t_retrieve_tetm.c
  private/h5t_retrieve_trim.c
  private/h5t_store_tetm.c
  private/h5t_store_trim.c
)
target_include_directories(H5hut
  PUBLIC "${MPI_C_INCLUDE_DIRS}" "${HDF5_Fortran_INCLUDE_DIRS}"
  . ../include
)
install(TARGETS H5hut ARCHIVE)
