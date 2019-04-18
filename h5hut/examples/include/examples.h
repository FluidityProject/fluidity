#ifndef __EXAMPLES_H
#define __EXAMPLES_H

#if !defined (H5_HAVE_PARALLEL)

#define MPI_COMM_WORLD (0)
#define MPI_Init(argc, argv)
#define MPI_Comm_size(comm,nprocs) {comm = 0; *nprocs = 1;}
#define MPI_Comm_rank(comm,rank)  {comm = 0; (void)(comm); *rank = 0;}
#define MPI_Finalize()

#endif
#endif
