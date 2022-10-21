!
!  Copyright (c) 2006-2015, The Regents of the University of California,
!  through Lawrence Berkeley National Laboratory (subject to receipt of any
!  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
!  Institut (Switzerland).  All rights reserved.!
!
!  License: see file COPYING in top level of source distribution.
!
include 'H5hut.f90'

program write_setview
  use H5hut
  implicit none

#if defined(PARALLEL_IO)
  include 'mpif.h'
#endif
  
  ! name of output file
  character (len=*), parameter :: fname = "example_setview.h5"

  ! H5hut verbosity level
  integer*8, parameter :: h5_verbosity = H5_VERBOSE_DEFAULT

  ! we are going to write multiple consecutive blocks
  integer*8, parameter :: num_blocks = 4;
  integer*8, parameter :: num_particles_per_block = 32

  integer   :: comm_rank = 0
  integer*8 :: file, h5_ierror
  integer*8 :: i, j, offset
  integer*4, allocatable :: data(:)

  ! initialize MPI & H5hut
#if defined(PARALLEL_IO)
  integer   :: comm, mpi_ierror
  comm = MPI_COMM_WORLD
  call mpi_init (mpi_ierror)
  call mpi_comm_rank (comm, comm_rank, mpi_ierror)
#endif
  
  call h5_abort_on_error ()
  call h5_set_verbosity_level (h5_verbosity)

  ! open file and create first step
  file = h5_openfile (FNAME, H5_O_WRONLY, H5_PROP_DEFAULT)
  h5_ierror = h5_setstep(file, 1_8)

  ! If we want to write consecutive blocks, the 'view' can be defined
  ! with H5PartSetview(). Otherwise we have to define the total number
  ! of particles with H5PartSetNumParticles().
  offset = comm_rank * num_blocks * num_particles_per_block+1
  h5_ierror = h5pt_setview (file, offset, offset + num_blocks*num_particles_per_block - 1)

  !  write multiple consecutive blocks
  allocate (data (num_particles_per_block))
  do i = 1, num_blocks
     ! create fake data
     do j = 1, num_particles_per_block
        data (j) = int((j-1) + (i-1)*num_particles_per_block + offset - 1)
     end do
     h5_ierror = h5pt_setview (file, offset + (i-1)*num_particles_per_block, offset - 1 + i*num_particles_per_block)
     ! write data
     h5_ierror = h5pt_writedata_i4 (file, "data", data)
  end do

  ! cleanup
  deallocate (data)
  h5_ierror = h5_closefile (file)

#if defined(PARALLEL_IO)
  call mpi_finalize (mpi_ierror)
#endif
  
end program write_setview
