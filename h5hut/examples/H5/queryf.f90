  !
  !  Copyright (c) 2006-2017, The Regents of the University of California,
  !  through Lawrence Berkeley National Laboratory (subject to receipt of any
  !  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  !  Institut (Switzerland).  All rights reserved.!
  !
  !  License: see file COPYING in top level of source distribution.
  !
  include 'H5hut.f90'
  
  program query
    use H5hut
    implicit none
#if defined(PARALLEL_IO)
    include 'mpif.h'
#endif

    ! the file name we want to read
    character (len=*), parameter :: FNAME1 = "example_file_attribs.h5"
    character (len=*), parameter :: FNAME2 = "example_step_attribs.h5"

    ! verbosity level: set it to a power of 2 minus one or zero
    integer*8, parameter :: verbosity_level =             1

    ! used for mpi error return
    integer :: ierr

#if defined(PARALLEL_IO)
    call mpi_init (ierr)
#endif

    ! abort program on any H5hut error
    call h5_abort_on_error ()

    call h5_set_verbosity_level (verbosity_level)

    call query_file (FNAME1);
    call query_file (FNAME2);

#if defined(PARALLEL_IO)
    call mpi_finalize(ierr)
#endif
    call exit (ierr)

  contains

    subroutine query_file (fname)
      character(len=*), intent(in):: fname

      integer*8 file_id
      integer*8 i, n, status

      write (*, '("File: ", a)') fname

      ! if file properties is set to default, MPI_COMM_WORLD will be used
      file_id = h5_openfile (fname, H5_O_RDONLY, H5_PROP_DEFAULT)

      ! query and output file attribs
      call query_file_attribs (file_id);

      ! query # of steps, if > 0: go to first step, query and output step attribs
      n = h5_getnsteps (file_id);
      write (*, '(T8, "Number of steps: ", I0)') n

      if (n > 0) then
         ! go to first step 
         i = 0;
         do
            if (h5_hasstep (file_id, i)) exit
            i = i+1
         end do
         call query_step_attribs (file_id, i);
      end if
      status = h5_closefile (file_id)
    end subroutine query_file

    ! print header for attribute metadata table
    subroutine print_header (n)
      integer*8, intent(in):: n

      character(len=6),  parameter :: idx =  "idx"
      character(len=30), parameter :: name = "name"
      character(len=15), parameter :: type = "type"
      character(len=10), parameter :: dim =  "dim"

      if (n > 0) then
         write (*, '(T8, A, 1X, A, 1X, A, 1X, A)') idx, name, type, dim
      end if

    end subroutine print_header

    ! output attribute metadata
    subroutine print_query_result (i, name, type, dim)
      integer*8, intent(in):: i
      character(len=*), intent(in):: name
      integer*8, intent(in):: type
      integer*8, intent(in):: dim

      character(len=30) name_
      character(len=15) type_, type_char
      character(len=6)  i_
      character(len=10) dim_

      if (type == H5_FLOAT64_T) then
         type_char = "H5_FLOAT64_T"
      else if (type == H5_FLOAT32_T) then
         type_char = "H5_FLOAT32_T"
      else if (type == H5_INT64_T) then
         type_char = "H5_INT64_T"
      else if (type == H5_INT32_T) then
         type_char = "H5_INT32_T"
      else if (type == H5_STRING_T) then
         type_char = "H5_STRING_T"
      else
         type_char = "unknown type"
      end if

      write (i_,    '(I6)')  i
      write (name_, '(A30)') name
      write (type_, '(A15)') type_char
      write (dim_,  '(I10)') dim

      write (*, '(T8, A6, 1X, A30, 1X, A15, 1X, A10)')  adjustl(i_), adjustl(name_), adjustl(type_), adjustl(dim_)
    end subroutine print_query_result

    subroutine query_file_attribs (file_id)
      integer*8, intent(in):: file_id

      integer*8 status
      integer*8 i, n
      character(len=H5_MAX_NAME_LEN) name
      integer*8 type, dim

      ! query # of file attributes
      n = h5_getnfileattribs (file_id);
      write  (*, '(T8, "Number of file attributes: ", I0)') n

      ! output name and type of all file attribute
      call print_header (n);
      do i = 1, n
         status = h5_getfileattribinfo (file_id, i, name, type, dim);
         call print_query_result (i, name, type, dim);
      end do
      write (*,*)

    end subroutine query_file_attribs

    subroutine query_step_attribs (file_id, stepno)
      integer*8, intent(in):: file_id
      integer*8, intent(in):: stepno

      integer*8 status
      integer*8 i, n
      character(len=H5_MAX_NAME_LEN) name
      integer*8 type, dim

      ! Go to step #1
      status = h5_setstep (file_id, stepno);

      ! query # of step attributes
      n = h5_getnstepattribs (file_id)
      write  (*, '(T8, "Number of step attributes: ", i0)') n

      ! output name and type of all step attribute
      call print_header (n)
      do i = 1, n
         status = h5_getstepattribinfo (file_id, i, name, type, dim)
         call print_query_result (i, name, type, dim)
      end do

    end subroutine query_step_attribs
  end program query
