#include "detector_parallelc.h"

// Functions for a shared MPI counterfor creating unique detector IDs in parallel

struct mpi_counter_t *id_counter;

/* Initialise the shared ID counter for generating detector IDs in parallel */
extern "C" void init_id_counter_c(){
  int i;
  id_counter = (struct mpi_counter_t *)malloc(sizeof(struct mpi_counter_t));

  MPI_Comm_rank(MPI_COMM_WORLD, &(id_counter->rank));
  MPI_Comm_size(MPI_COMM_WORLD, &(id_counter->size));

  // Rank 0 is always the host, where we initialise the data array
  if(id_counter->rank == 0){
     MPI_Alloc_mem(id_counter->size * sizeof(int), MPI_INFO_NULL, &(id_counter->data));
     for(i=0; i<id_counter->size; i++){
        id_counter->data[i] = 0;
     }
     MPI_Win_create(id_counter->data, id_counter->size * sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &(id_counter->win));
  }else{
     id_counter->data = NULL;
     MPI_Win_create(id_counter->data, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &(id_counter->win));
  }
  id_counter->myval = 0;
}

/* Get the next detector ID from the shared MPI counter */
extern "C" void get_next_detector_id_c(int *next_id){
  int *vals = (int *)malloc( id_counter->size * sizeof(int) );
  int i, val;
  int increment=1;

  MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, id_counter->win);

  for(i=0; i<id_counter->size; i++) {
     if(i == id_counter->rank){
        MPI_Accumulate(&increment, 1, MPI_INT, 0, i, 1, MPI_INT, MPI_SUM, id_counter->win);
     }else{
        MPI_Get(&vals[i], 1, MPI_INT, 0, i, 1, MPI_INT, id_counter->win);
     }
  }

  MPI_Win_unlock(0, id_counter->win);

  id_counter->myval += increment;
  vals[id_counter->rank] = id_counter->myval;
  *next_id = 0;
  for(i=0; i<id_counter->size; i++){
      *next_id += vals[i];
  }

  free(vals);
}

/* Delete the shared ID counter */
extern "C" void delete_id_counter_c() {
    // Rank 0 is always the data host
    if(id_counter->rank == 0){
        MPI_Free_mem(id_counter->data);
    }
    MPI_Win_free(&(id_counter->win));
    free(id_counter);
}
