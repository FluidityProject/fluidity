!
!  Copyright (c) 2006-2015, The Regents of the University of California,
!  through Lawrence Berkeley National Laboratory (subject to receipt of any
!  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
!  Institut (Switzerland).  All rights reserved.!
!
!  License: see file COPYING in top level of source distribution.
!
include 'H5hut.f90'

program read_stridedf
  use H5hut
  implicit none
#if defined(PARALLEL_IO)
  include 'mpif.h'
#endif
  
  ! name of input file
  character (len=*), parameter :: fname = "example_strided.h5"

  ! H5hut verbosity level
  integer*8, parameter :: h5_verbosity = H5_VERBOSE_DEFAULT
  
  integer*8 :: file, h5_ierror
  integer*8 :: num_particles, num_particles_total
  real*8, allocatable :: data(:)
  integer*8 :: i, start
  
  ! initialize MPI & H5hut
#if defined(PARALLEL_IO)
  integer   :: comm, comm_size, comm_rank, mpi_ierror
  comm = MPI_COMM_WORLD
  call mpi_init (mpi_error)
  call mpi_comm_size (comm, comm_size, mpi_ierror)
  call mpi_comm_rank (comm, comm_rank, mpi_ierror)
#else
  integer   :: comm_size = 1
  integer   :: comm_rank = 0
#endif
  
  call h5_abort_on_error ()
  call h5_set_verbosity_level (h5_verbosity)

  ! open file and go to first step
  file = h5_openfile (fname, H5_O_RDONLY, H5_PROP_DEFAULT)
  h5_ierror = h5_setstep(file, 1_8)

  ! compute number of particles this process has to read
  num_particles_total = h5pt_getnpoints (file)
  num_particles = num_particles_total / comm_size
  if (comm_rank+1 == comm_size) then
     num_particles = num_particles + mod (num_particles_total, comm_size)
  end if

  write (*, "('Total number of particles: ', i8)") num_particles_total
  write (*, "('Number of particles on this core: ', i8)") num_particles

  ! set number of particeles and memory stride
  h5_ierror = h5pt_setnpoints_strided (file, num_particles, 6_8)

  ! read data
  allocate (data (6*num_particles))
  h5_ierror = h5pt_readdata_r8 (file, "x",  data(1:))
  h5_ierror = h5pt_readdata_r8 (file, "y",  data(2:))
  h5_ierror = h5pt_readdata_r8 (file, "z",  data(3:))
  h5_ierror = h5pt_readdata_r8 (file, "px", data(4:))
  h5_ierror = h5pt_readdata_r8 (file, "py", data(5:))
  h5_ierror = h5pt_readdata_r8 (file, "pz", data(6:))

  ! print dataset "x"
  start = 1
  do i = start, num_particles*6, 6
     write (*, "('[proc ', i4, ']: global index = ', i4, '; local index = ', i4, ', value = ', f10.2)") &
          comm_rank, start+i-2, i, data(i)
  end do
  
  ! cleanup
  deallocate (data)
  h5_ierror = h5_closefile (file)

#if defined(PARALLEL_IO)
  call mpi_finalize (mpi_ierror)
#endif
  
end program read_stridedf
