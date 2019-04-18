!
!  Copyright (c) 2006-2015, The Regents of the University of California,
!  through Lawrence Berkeley National Laboratory (subject to receipt of any
!  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
!  Institut (Switzerland).  All rights reserved.!
!
!  License: see file COPYING in top level of source distribution.
!
include 'H5hut.f90'

program read_core_vfd
  use H5hut
  implicit none

#if defined(PARALLEL_IO)
  include 'mpif.h'
#endif

  ! name of input file
  character (len=*), parameter :: fname = "example_core_vfd.h5"

  ! H5hut verbosity level
  integer*8, parameter :: h5_verbosity = H5_VERBOSE_DEFAULT

  integer*8 :: file, h5_ierror
  integer*8 :: prop
  integer*8 :: num_particles
  integer*4 :: i
  integer*4, allocatable :: data(:)

  ! initialize MPI & H5hut
  integer   :: comm_rank = 0
#if defined(PARALLEL_IO)
  integer   :: comm, mpi_ierror
  comm = MPI_COMM_WORLD
  call mpi_init (mpi_ierror)
  call mpi_comm_rank (comm, comm_rank, mpi_ierror)
#endif
  call h5_abort_on_error ()
  call h5_set_verbosity_level (h5_verbosity)

  ! open file and create first step
  prop = h5_createprop_file ()
  h5_ierror = h5_setprop_file_corevfd (prop);
  file = h5_openfile (fname, H5_O_RDONLY, prop)
  h5_ierror = h5_setstep(file, 1_8)

  ! with core cfd we read the hole file on all cores!
  num_particles = h5pt_getnpoints (file)
  write (*, "('[proc ', i4, ']: particles in view: ', i4)") &
       comm_rank, num_particles
  
  ! read and print data
  allocate (data (num_particles))
  h5_ierror = h5pt_readdata_i4 (file, "data", data)
  do i = 1, int (num_particles)
     write (*, "('[proc ', i4, ']: local index = ', i4, ', value = ', i4)") &
          comm_rank, i, data(i)
  end do

  ! cleanup
  deallocate (data)
  h5_ierror = h5_closefile (file)

#if defined(PARALLEL_IO)
  call mpi_finalize (mpi_ierror)
#endif

end program read_core_vfd
