!
!  Copyright (c) 2006-2015, The Regents of the University of California,
!  through Lawrence Berkeley National Laboratory (subject to receipt of any
!  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
!  Institut (Switzerland).  All rights reserved.!
!
!  License: see file COPYING in top level of source distribution.
!
include 'H5hut.f90'

program openclose

  use H5hut

  implicit none

#if defined(PARALLEL_IO)
  include 'mpif.h'
#endif
  integer :: comm, rank
  integer*8 :: file_id, status
  integer*8 :: props
  
#if defined(PARALLEL_IO)
  integer :: ierr
  comm = MPI_COMM_WORLD
  call mpi_init(ierr)
  call mpi_comm_rank(comm, rank, ierr)
#else
  comm = 0
  rank = 1
#endif
  
  props = h5_createprop_file ()
#if defined(PARALLEL_IO)
  status = h5_setprop_file_mpio_collective (props, comm)
#endif
  file_id = h5_openfile ("testfile.h5", H5_O_WRONLY, props)
  status = h5_closeprop (props)
  status = h5_closefile (file_id);

#if defined(PARALLEL_IO)
  call mpi_finalize(ierr)
#endif
  
end program openclose
