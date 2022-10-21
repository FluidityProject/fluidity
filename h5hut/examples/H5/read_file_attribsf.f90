  !
  !  Copyright (c) 2006-2015, The Regents of the University of California,
  !  through Lawrence Berkeley National Laboratory (subject to receipt of any
  !  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  !  Institut (Switzerland).  All rights reserved.!
  !
  !  License: see file COPYING in top level of source distribution.
  !
  include 'H5hut.f90'

  program read_file_attribs
    use H5hut
    implicit none
#if defined(PARALLEL_IO)    
    include 'mpif.h'
#endif

    ! the file name we want to read
    character (len=*), parameter :: FNAME =               "example_file_attribs.h5"

    ! verbosity level: set it to
    ! - 1 to see error messages, if something goes wrong
    ! - 0 to get no output
    ! - a power of 2 minus one to get lot of output
    ! zeror zero
    integer*8, parameter :: verbosity_level =             1

    ! we know the attribute names!
    character (len=*), parameter :: ATTR_STRING =         "FileAttrString"
    character (len=*), parameter :: ATTR_I4 =             "FileAttrInt32"
    character (len=*), parameter :: ATTR_I8 =             "FileAttrInt64"
    character (len=*), parameter :: ATTR_R4 =             "FileAttrFloat32"
    character (len=*), parameter :: ATTR_R8 =             "FileAttrFloat64"

    ! for formated output
    character (len=128) ::    fmt

    ! attribute values. Note: allocatable strings aren't supported in Fortran90
    character (len=256) ::    string_value
    integer*4, allocatable :: i4_value (:)
    integer*8, allocatable :: i8_value (:)
    real*4, allocatable ::    r4_value (:)
    real*8, allocatable ::    r8_value (:)

    ! H5hut file id
    integer*8 :: file_id

    ! H5hut API status return 
    integer*8 status

    ! type of attribute
    integer*8 type

    ! len of attribute
    integer*8 len

    ! loop index
    integer*8 i

#if defined(PARALLEL_IO)
    ! used for mpi error return
    integer :: ierr

    call mpi_init (ierr)
#endif

    ! abort program on any H5hut error
    call h5_abort_on_error ()

    call h5_set_verbosity_level (verbosity_level)

    !  MPI_COMM_WORLD is used, if file is opened with default properties
    file_id = h5_openfile (FNAME, H5_O_RDONLY, H5_PROP_DEFAULT)

    ! read and output string attribute
    status = h5_getfileattribinfo_by_name (file_id, ATTR_STRING, type, len)
    status = h5_readfileattrib_string (file_id, ATTR_STRING, string_value)
    write (fmt, "(a, i0, a)")  '(a', len, ')'
    write (*, "(a, ' = ')", advance='no') ATTR_STRING
    write (*, fmt) string_value

    ! read and output 32bit integer attribute
    status = h5_getfileattribinfo_by_name (file_id, ATTR_I4, type, len)
    allocate (i4_value(len))
    status = h5_readfileattrib_i4 (file_id, ATTR_I4, i4_value)
    write (fmt, "(a, i0, a)")  '(', len, 'i4)'
    write (*, "(a, ' =')", advance='no') ATTR_I4
    write (*, fmt) (i4_value(i), i = 1, len)

    ! read and output 64bit integer attribute
    status = h5_getfileattribinfo_by_name (file_id, ATTR_I8, type, len)
    allocate (i8_value(len))
    status = h5_readfileattrib_i8 (file_id, ATTR_I8, i8_value)
    write (fmt, "(a, i0, a)")  '(', len, 'i4)'
    write (*, "(a, ' =')", advance='no') ATTR_I8
    write (*, fmt) (i8_value(i), i = 1, len)

    ! read and output 32bit floating point attribute
    status = h5_getfileattribinfo_by_name (file_id, ATTR_R4, type, len)
    allocate (r4_value(len))
    status = h5_readfileattrib_r4 (file_id, ATTR_R4, r4_value)
    write (fmt, "(a, i0, a)")  '(', len, 'f10.5)'
    write (*, "(a, ' =')", advance='no') ATTR_R4
    write (*, fmt) (r4_value(i), i = 1, len)

    ! read and output 64bit floating point attribute
    status = h5_getfileattribinfo_by_name (file_id, ATTR_R8, type, len)
    allocate (r8_value(len))
    status = h5_readfileattrib_r8 (file_id, ATTR_R8, r8_value)
    write (fmt, "(a, i0, a)")  '(', len, 'f10.5)'
    write (*, "(a, ' =')", advance='no') ATTR_R8
    write (*, fmt) (r8_value(i), i = 1, len)

    ! cleanup
    status = h5_closefile (file_id)
#if defined(PARALLEL_IO)
    call mpi_finalize(ierr)
#endif
  end program read_file_attribs
