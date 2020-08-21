!
!  Copyright (c) 2006-2015, The Regents of the University of California,
!  through Lawrence Berkeley National Laboratory (subject to receipt of any
!  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
!  Institut (Switzerland).  All rights reserved.!
!
!  License: see file COPYING in top level of source distribution.
!
include 'H5hut.f90'

program read_write_scalar_field
  use H5hut
  implicit   none

#if defined(PARALLEL_IO)
  include 'mpif.h'

  integer :: comm = 0
  integer :: mpi_err
#endif
  
  integer :: nargs = 0
  integer :: comm_rank = 0
  integer :: comm_size = 1
  integer*8 :: h5_err
  integer :: i
  character(len=32) :: arg_str
  integer :: opt_read = 0
  integer :: opt_write = 0
  integer :: opt_with_ghosts = 0
  character(len=128) :: fname
  integer*8 :: layout   (6)
  integer*8 :: layout1  (6,1)
  integer*8 :: layout8  (6,8)
  integer*8 :: layout8g (6,8)
  integer*8 :: layout16 (6,16)
  integer*8 :: layout16g(6,16)
  integer*8 :: layout32 (6,32)
  integer*8 :: layout32g(6,32)
  
  data layout1  / 1,64,  1,64,   1,512 /

  data layout8  / 1,64,  1,64,   1, 64, &
                  1,64,  1,64,  65,128, &
                  1,64,  1,64, 129,192, &
                  1,64,  1,64, 193,256, &
                  1,64,  1,64, 257,320, &
                  1,64,  1,64, 321,384, &
                  1,64,  1,64, 385,448, &
                  1,64,  1,64, 449,512 /

  data layout8g / 1,64,  1,64,   1, 65, &
                  1,64,  1,64,  64,129, &
                  1,64,  1,64, 128,193, &
                  1,64,  1,64, 192,257, &
                  1,64,  1,64, 256,321, &
                  1,64,  1,64, 320,385, &
                  1,64,  1,64, 384,449, &
                  1,64,  1,64, 448,512  /

  data layout16 / 1,64,  1,32,   1, 64, &
                  1,64, 33,64,   1, 64, &
                  1,64,  1,32,  65,128, &
                  1,64, 33,64,  65,128, &
                  1,64,  1,32, 129,192, &
                  1,64, 33,64, 129,192, &
                  1,64,  1,32, 193,256, &
                  1,64, 33,64, 193,256, &
                  1,64,  1,32, 257,320, &
                  1,64, 33,64, 257,320, &
                  1,64,  1,32, 321,384, &
                  1,64, 33,64, 321,384, &
                  1,64,  1,32, 385,448, &
                  1,64, 33,64, 385,448, &
                  1,64,  1,32, 449,512, &
                  1,64, 33,64, 449,512 /

  data layout16g/ 1,64,  1,33,   1, 65, &
                  1,64, 32,64,   1, 65, &
                  1,64,  1,33,  64,129, &
                  1,64, 32,64,  64,129, &
                  1,64,  1,33, 128,193, &
                  1,64, 32,64, 128,193, &
                  1,64,  1,33, 192,257, &
                  1,64, 32,64, 192,257, &
                  1,64,  1,33, 256,321, &
                  1,64, 32,64, 256,321, &
                  1,64,  1,33, 320,385, &
                  1,64, 32,64, 320,385, &
                  1,64,  1,33, 384,449, &
                  1,64, 32,64, 384,449, &
                  1,64,  1,33, 448,512, &
                  1,64, 32,64, 448,512  /

  data layout32 / 1,32,  1,32,   1, 64, &
                  1,32, 33,64,   1, 64, &
                 33,64,  1,32,   1, 64, &
                 33,64, 33,64,   1, 64, &
                  1,32,  1,32,  65,128, &
                  1,32, 33,64,  65,128, &
                 33,64,  1,32,  65,128, &
                 33,64, 33,64,  65,128, &
                  1,32,  1,32, 129,192, &
                  1,32, 33,64, 129,192, &
                 33,64,  1,32, 129,192, &
                 33,64, 33,64, 129,192, &
                  1,32,  1,32, 193,256, &
                  1,32, 33,64, 193,256, &
                 33,64,  1,32, 193,256, &
                 33,64, 33,64, 193,256, &
                  1,32,  1,32, 257,320, &
                  1,32, 33,64, 257,320, &
                 33,64,  1,32, 257,320, &
                 33,64, 33,64, 257,320, &
                  1,32,  1,32, 321,384, &
                  1,32, 33,64, 321,384, &
                 33,64,  1,32, 321,384, &
                 33,64, 33,64, 321,384, &
                  1,32,  1,32, 385,448, &
                  1,32, 33,64, 385,448, &
                 33,64,  1,32, 385,448, &
                 33,64, 33,64, 385,448, &
                  1,32,  1,32, 449,512, &
                  1,32, 33,64, 449,512, &
                 33,64,  1,32, 449,512, &
                 33,64, 33,64, 449,512  /

  data layout32g/ 1,33,  1,33,   1, 65, &
                  1,33, 32,64,   1, 65, &
                 32,64,  1,33,   1, 65, &
                 32,64, 32,64,   1, 65, &
                  1,33,  1,33,  64,129, &
                  1,33, 32,64,  64,129, &
                 32,64,  1,33,  64,129, &
                 32,64, 32,64,  64,129, &
                  1,33,  1,33, 128,193, &
                  1,33, 32,64, 128,193, &
                 32,64,  1,33, 128,193, &
                 32,64, 32,64, 128,193, &
                  1,33,  1,33, 192,257, &
                  1,33, 32,64, 192,257, &
                 32,64,  1,33, 192,257, &
                 32,64, 32,64, 192,257, &
                  1,33,  1,33, 256,321, &
                  1,33, 32,64, 256,321, &
                 32,64,  1,33, 256,321, &
                 32,64, 32,64, 256,321, &
                  1,33,  1,33, 320,385, &
                  1,33, 32,64, 320,385, &
                 32,64,  1,33, 320,385, &
                 32,64, 32,64, 320,385, &
                  1,33,  1,33, 384,449, &
                  1,33, 32,64, 384,449, &
                 32,64,  1,33, 384,449, &
                 32,64, 32,64, 384,449, &
                  1,33,  1,33, 448,512, &
                  1,33, 32,64, 448,512, &
                 32,64,  1,33, 448,512, &
                 32,64, 32,64, 448,512  /
  nargs = iargc ()
  if (nargs == 0) then
     print *, "usage: read_write_scalarfield -w | -r [-g]"      
     call exit (1)
  end if
  do i = 1, nargs
     call getarg (i, arg_str)
     if (arg_str == "-r") then
        opt_read = 1
     else if (arg_str == "-w") then
        opt_write = 1
     else if (arg_str == "-g") then
        opt_with_ghosts = 1
     else
        print *, "Illegal option ", arg_str, "\n"
        print *, "Usage: read_write_scalarfield -w | -r [-g]"
        call exit (1)
     end if
  end do
  
  ! init MPI & H5hut
#if defined(PARALLEL_IO)
  comm = MPI_COMM_WORLD
  call mpi_init(mpi_err)
  call mpi_comm_rank(comm, comm_rank, mpi_err)
  call mpi_comm_size (comm, comm_size, mpi_err)
#else
  comm_size = 1
  comm_rank = 0
#endif
  call h5_abort_on_error ()
  call h5_set_verbosity_level (511_8)

  selectcase (comm_size)
  case (1)
     fname = "blockfile1.h5"
     layout = layout1 (:, comm_rank+1)

  case (8)
     if (opt_with_ghosts == 1) then
        fname = "blockfile8g.h5"
        layout = layout8g (:, comm_rank+1)
     else
        fname = "blockfile8.h5"
        layout = layout8 (:, comm_rank+1)
     end if
     
  case (16)
     if (opt_with_ghosts == 1) then
        fname = "blockfile16g.h5"
        layout = layout16g (:, comm_rank+1)
     else
        fname = "blockfile16.h5"
        layout = layout16 (:, comm_rank+1)
     end if
     
  case (32)
     if (opt_with_ghosts == 1) then
        fname = "blockfile32g.h5"
        layout = layout32g (:, comm_rank+1)
     else
        fname = "blockfile32.h5"
        layout = layout32 (:, comm_rank+1)
     end if

  case default
     print *, "Run this test on 1, 8, 16 or 32 cores!"
#if defined(PARALLEL_IO)
     call mpi_finalize
#endif
     call exit (1)
  end select

  if (opt_write == 1) then
     h5_err = write_file (fname, comm_rank, layout)
     if (h5_err < 0) then
        print "('[proc ', I3, ']: Faild to write file ', A, '!')", comm_rank, fname
     end if
     
  else if (opt_read == 1) then
     h5_err = read_file (fname, comm_rank, layout)
     if (h5_err < 0) then
        print "('[proc ', I3, ']: Faild to read file ', A, '!')", comm_rank, fname
     end if
     
  endif
  print "('[proc ', I3, ']: Cleanup.')", comm_rank
#if defined(PARALLEL_IO)
  call mpi_finalize
#endif
  print "('[proc ', I3, ']: Done.')", comm_rank
  call exit (0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
   contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   integer*8 function write_file (fname, comm_rank, layout)
     use H5hut
     implicit   none

     character(len=*), intent(in) :: fname
     integer, intent(in) ::          comm_rank
     integer*8, intent(in) ::        layout(6)

     integer*8 :: file
     integer*8 :: timestep = 1

     print "('[proc ', I3, ']: Open file for writing ...')", comm_rank
     file = h5_openfile (fname, H5_O_WRONLY, H5_PROP_DEFAULT)
     h5_err = h5_setstep (file, timestep)
     h5_err = write_field (file, comm_rank, layout)
     h5_err = write_attributes (file)
     h5_err = h5_closefile (file)

     write_file = H5_SUCCESS
   end function write_file

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   integer*8 function write_field (file, comm_rank, layout)
     use H5hut
     implicit none

     integer*8, intent(in) :: file
     integer, intent(in) ::   comm_rank
     integer*8, intent(in) :: layout(6)

     integer*8 :: i, j, k
     integer*8 :: i_start
     integer*8 :: i_end
     integer*8 :: j_start
     integer*8 :: j_end
     integer*8 :: k_start
     integer*8 :: k_end
     integer*8 :: i_dims
     integer*8 :: j_dims
     integer*8 :: k_dims
     real*8 :: value
  
     real*8, dimension(:,:,:), allocatable :: data

     i_start = layout(1)
     i_end   = layout(2)
     j_start = layout(3)
     j_end   = layout(4)
     k_start = layout(5)
     k_end   = layout(6)
     i_dims  = i_end - i_start + 1
     j_dims  = j_end - j_start + 1
     k_dims  = k_end - k_start + 1

     allocate ( data (i_dims,j_dims, k_dims) )

     print "('[proc ', I3, ']: Defining layout for writing ...')", comm_rank
     print "('[proc ', I3, ']: ', I3, ':', I3, ', ', I3, ':', I3,', ', I3, ':', I3)", &
          comm_rank, &
          i_start, i_end, &
          j_start, j_end, &
          k_start, k_end

     h5_err = h5bl_3d_setview (file, i_start, i_end, j_start, j_end, k_start, k_end)

     do i = 1, i_dims
        do j = 1, j_dims
           do k = 1, k_dims
              value = (k-1) + 1000*(j-1) + 100000*(i-1) + 10000000*comm_rank
              data(i,j,k) = value
           end do
        end do
     end do

     print "('[proc ', I3, ']: Writing field ...')", comm_rank
     h5_err = h5bl_3d_write_scalar_field_r8 ( file, "TestField", data )

     write_field = 0
   END FUNCTION write_field

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   integer*8 function write_attributes (file)
     use H5hut
     implicit none
     integer*8, intent(in) :: file
     integer*8 :: h5_err = 0
     character(len=128) :: s_val
     integer*8 :: i8_val(1)
     integer*4 :: i4_val(1)
     real*8 :: r8_val(1)
     real*4 :: r4_val(1)

     print "('[proc ', I3, ']: Writing string attribute ...')", comm_rank
     s_val = "42"
     h5_err = h5bl_writefieldattrib_string ( file, "TestField", "TestString", s_val ) 

     print "('[proc ', I3, ']: Writing int64 attribute ...')", comm_rank
     i8_val(1) = 42
     h5_err = h5bl_writefieldattrib_i8 ( file, "TestField", "TestInt64", i8_val, 1_8 )

     print "('[proc ', I3, ']: Writing int32 attribute ...')", comm_rank
     i4_val(1) = 42
     h5_err = h5bl_writefieldattrib_i4 ( file, "TestField", "TestInt32", i4_val, 1_8 )

     print "('[proc ', I3, ']: Writing float64 attribute ...')", comm_rank
     r8_val(1) = 42.0
     h5_err = h5bl_writefieldattrib_r8 ( file, "TestField", "TestFloat64", r8_val, 1_8 )

     print "('[proc ', I3, ']: Writing float32 attribute ...')", comm_rank
     r4_val(1) = 42.0
     h5_err = h5bl_writefieldattrib_r4 ( file, "TestField", "TestFloat32", r4_val, 1_8 )
 
     write_attributes = 0
   end function write_attributes

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   integer*8 function read_file (fname, comm_rank, layout)
     use H5hut
     implicit none

     character(len=*), intent(in) :: fname
     integer, intent(in) ::          comm_rank
     integer*8, intent(in) ::        layout(6)
     integer*8 :: file
     integer*8 :: timestep = 1

     print "('[proc ', I3, ']: Open file for reading ...')", comm_rank
     file = h5_openfile (FNAME, H5_O_RDONLY, H5_PROP_DEFAULT)
     h5_err = h5_setstep (file, timestep)

     h5_err =  read_field (file, comm_rank, layout)
     h5_err = read_attributes (file)
     h5_err = h5_closefile (file)

     read_file = 0
   end function read_file

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   integer*8 function read_field (file, comm_rank, layout)
     use H5hut
     implicit none

     integer*8, intent(in) :: file
     integer, intent(in) ::          comm_rank
     integer*8, intent(in) ::        layout(6)

     integer*8 :: i, j, k
     integer*8 :: i_start
     integer*8 :: i_end
     integer*8 :: j_start
     integer*8 :: j_end
     integer*8 :: k_start
     integer*8 :: k_end
     integer*8 :: i_dims
     integer*8 :: j_dims
     integer*8 :: k_dims
     real*8 :: value
  
     real*8, dimension(:,:,:), allocatable :: data

     i_start = layout(1)
     i_end   = layout(2)
     j_start = layout(3)
     j_end   = layout(4)
     k_start = layout(5)
     k_end   = layout(6)
     i_dims  = i_end - i_start + 1
     j_dims  = j_end - j_start + 1
     k_dims  = k_end - k_start + 1

     allocate ( data (i_dims, j_dims, k_dims) )

     print "('[proc ', I3, ']: Defining layout for reading ...')", comm_rank
     print "('[proc ', I3, ']: ', I3, ':', I3, ', ', I3, ':', I3,', ', I3, ':', I3)", &
          comm_rank, &
          i_start, i_end, &
          j_start, j_end, &
          k_start, k_end

     h5_err = h5bl_3d_setview ( file, i_start, i_end, j_start, j_end, k_start, k_end )

     print "('[proc ', I3, ']: Reading field ...')", comm_rank
     h5_err = h5bl_3d_read_scalar_field_r8 ( file, "TestField", data )

     do i = 1, i_dims
        do j = 1, j_dims
           do k = 1, k_dims
              value = (k-1) + 1000*(j-1) + 100000*(i-1) + 10000000*comm_rank
              if (data(i,j,k) /= value) then
                 print "('[proc ', I3, ']: error: data(',I4,',',I4,',',I4,') = ',F10.2,' /= ',F10.2)", &
                      i, j, k, data(i,j,k), value
                 read_field = -2
                 return
              end if
           end do
        end do
     end do
     read_field = 0
   end function read_field

   integer*8 function read_attributes (file)
     use H5hut
     implicit none
     integer*8, intent(in) :: file

     integer*8 :: h5_err = 0
     character(len=128) :: s_val
     integer*8 :: i8_val(1)
     integer*4 :: i4_val(1)
     real*8 :: r8_val(1)
     real*4 :: r4_val(1)

     print "('[proc ', I3, ']: Reading string attribute ...')", comm_rank
     h5_err = h5bl_readfieldattrib_string ( file, "TestField", "TestString", s_val ) 
     IF ( s_val /= "42" ) THEN
        print "('[proc ', I3, ']: Error reading string attribute: Value is ', A, ' but should be 42')", &
             comm_rank, s_val
     end if

     print "('[proc ', I3, ']: Reading int64 attribute ...')", comm_rank
     h5_err = h5bl_readfieldattrib_i8 ( file, "TestField", "TestInt64", i8_val )
     if ( i8_val(1) /= 42 ) then
        print "('[proc ', I3, ']: Error reading int64 attribute: Value is ', I8, ' but should be 42')", &
             comm_rank, i8_val(1)
     end if

     print "('[proc ', I3, ']: Reading int32 attribute ...')", comm_rank
     h5_err = h5bl_readfieldattrib_i4 ( file, "TestField", "TestInt32", i4_val )
     if ( i4_val(1) /= 42 ) then
        print "('[proc ', I3, ']: Error reading int32 attribute: Value is ', I8, ' but should be 42')", &
             comm_rank, i4_val(1)
     end if

     print "('[proc ', I3, ']: Reading float64 attribute ...')", comm_rank
     h5_err = h5bl_readfieldattrib_r8 ( file, "TestField", "TestFloat64", r8_val )
     if ( r8_val(1) /= 42.0 ) then
        print "('[proc ', I3, ']: Error reading float64 attribute: Value is ', F10.2, ' but should be 42.0')", &
             comm_rank, r8_val(1)
     end if

     print "('[proc ', I3, ']: Reading float32 attribute ...')", comm_rank
     h5_err = h5bl_readfieldattrib_r4 ( file, "TestField", "TestFloat32", r4_val )
     if ( r4_val(1) /= 42.0 ) then
        print "('[proc ', I3, ']: Error reading float32 attribute: Value is ', F10.2, ' but should be 42.0')", &
             comm_rank, r4_val(1)
     end if

     read_attributes = h5_err
   end function read_attributes

 end program
