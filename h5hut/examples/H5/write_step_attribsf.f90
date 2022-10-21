  !
  !  Copyright (c) 2006-2015, The Regents of the University of California,
  !  through Lawrence Berkeley National Laboratory (subject to receipt of any
  !  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  !  Institut (Switzerland).  All rights reserved.!
  !
  !  License: see file COPYING in top level of source distribution.
  !
  include 'H5hut.f90'

  program write_step_attribs
    use H5hut
    implicit none

#if defined(PARALLEL_IO)
    include 'mpif.h'
#endif
    
    integer*8, parameter :: verbosity_level =          1
    character (len=*), parameter :: FNAME =            "example_step_attribs.h5"

    character (len=*), parameter :: ATTR_STRING =      "StepAttrString"
    character (len=*), parameter :: ATTR_I4 =          "StepAttrInt32"
    character (len=*), parameter :: ATTR_I8 =          "StepAttrInt64"
    character (len=*), parameter :: ATTR_R4 =          "StepAttrFloat32"
    character (len=*), parameter :: ATTR_R8 =          "StepAttrFloat64"

    character (len=*),parameter :: string_value =      "This is a string attribute bound to this step."

    integer*4, parameter, dimension(*) :: i4_value =   (/0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144/)
    integer*8, parameter, dimension(*) :: i8_value =   (/42, 43, 44, 45/)
    real*4, parameter, dimension(*) ::    r4_value =   (/2.71828/)
    real*8, parameter, dimension(*) ::    r8_value =   (/3.141592653589793238462643383279502884197169/)

    integer*8 :: file_id, status

#if defined(PARALLEL_IO)
    integer :: ierr
    call mpi_init(ierr)
#endif
    
    ! abort program on any H5hut error
    call h5_abort_on_error()

    call h5_set_verbosity_level (verbosity_level)

    !  MPI_COMM_WORLD is used, if file is opened with default properties
    file_id = h5_openfile (FNAME, H5_O_WRONLY, H5_PROP_DEFAULT)

    ! open step 1
    status = h5_setstep (file_id, int8(1))

    ! write attributes
    status = h5_writestepattrib_string (file_id, ATTR_STRING, string_value)
    status = h5_writestepattrib_i4     (file_id, ATTR_I4,     i4_value, int8(size(i4_value, 1)))
    status = h5_writestepattrib_i8     (file_id, ATTR_I8,     i8_value, int8(size(i8_value, 1)))
    status = h5_writestepattrib_r4     (file_id, ATTR_R4,     r4_value, int8(size(r4_value, 1)))
    status = h5_writestepattrib_r8     (file_id, ATTR_R8,     r8_value, int8(size(r8_value, 1)))

    ! cleanup
    status = h5_closefile (file_id)

#if defined(PARALLEL_IO)
    call mpi_finalize(ierr)
#endif
  end program write_step_attribs
