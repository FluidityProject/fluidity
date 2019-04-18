!
!  Copyright (c) 2006-2015, The Regents of the University of California,
!  through Lawrence Berkeley National Laboratory (subject to receipt of any
!  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
!  Institut (Switzerland).  All rights reserved.!
!
!  License: see file COPYING in top level of source distribution.
!
include 'H5hut.f90'

program write_stridedf
  use H5hut
  implicit none

#if defined(PARALLEL_IO)
  include 'mpif.h'
#endif
  
  ! the file name we want to read
  character (len=*), parameter :: FNAME =       "example_strided.h5"
  integer*8, parameter :: NPOINTS =             99

  integer :: rank = 0
  integer*8 :: file, status
  integer*4 :: i
  real*8, allocatable :: particles(:)
  integer*8, allocatable :: id(:)

  ! init MPI & H5hut
#if defined(PARALLEL_IO)
  integer :: comm, ierr
  comm = MPI_COMM_WORLD
  call mpi_init(ierr)
  call mpi_comm_rank(comm, rank, ierr)
#endif

  call h5_abort_on_error ()

  ! create fake data
  allocate(particles(6*NPOINTS), id(NPOINTS))
  do i = 0, NPOINTS-1
    particles (6*i + 1) = 0.0 + real (i + NPOINTS*rank)
    particles (6*i + 2) = 0.1 + real (i + NPOINTS*rank)
    particles (6*i + 3) = 0.2 + real (i + NPOINTS*rank)
    particles (6*i + 4) = 0.3 + real (i + NPOINTS*rank)
    particles (6*i + 5) = 0.4 + real (i + NPOINTS*rank)
    particles (6*i + 6) = 0.5 + real (i + NPOINTS*rank)
    id(i+1) = i+NPOINTS*rank
  enddo

  ! open the a file for parallel writing and ceate step #0
  file = h5_openfile (FNAME, H5_O_WRONLY, H5_PROP_DEFAULT)
  status = h5_setstep(file, 1_8)

  ! define number of items this processor will write and set the
  ! in-memory striding
  status = h5pt_setnpoints_strided (file, NPOINTS, 6_8)

  ! write strided data
  status = h5pt_writedata_r8 (file, "x",  particles(1))
  status = h5pt_writedata_r8 (file, "y",  particles(2))
  status = h5pt_writedata_r8 (file, "z",  particles(3))
  status = h5pt_writedata_r8 (file, "px", particles(4))
  status = h5pt_writedata_r8 (file, "py", particles(5))
  status = h5pt_writedata_r8 (file, "pz", particles(6))

  ! disable striding to write the ID's
  status = h5pt_setnpoints(file, NPOINTS)
  status = h5pt_writedata_i8(file, "id", id)

  ! cleanup
  status = h5_closefile (file)

  deallocate(particles, id)

#if defined(PARALLEL_IO)
  call mpi_finalize(ierr)
#endif

end program write_stridedf
