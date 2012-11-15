#include "confdefs.h"
#include <stdlib.h>
#include <mpi.h>

struct mpi_counter_t {
    MPI_Win win;
    int  myval;
    int *data;
    int rank, size;
};

/* Initialise the shared ID counter for generating detector IDs in parallel */
extern "C" void init_id_counter_c();

/* Get the next detector ID from the shared MPI counter */
extern "C" void get_next_detector_id_c(int *next_id);

/* Delete the shared ID counter */
extern "C" void delete_id_counter_c();
