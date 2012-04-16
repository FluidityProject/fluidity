! This is part of the netCDF package.
! Copyright 2006 University Corporation for Atmospheric Research/Unidata.
! See COPYRIGHT file for conditions of use.
      
! This is a simple example which reads a small dummy array, from a
! netCDF data file created by the companion program simple_xy_wr.f90.
      
! This is intended to illustrate the use of the netCDF fortran 77
! API. This example program is part of the netCDF tutorial, which can
! be found at:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial
      
! Full documentation of the netCDF Fortran 90 API can be found at:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

! $Id: simple_xy_rd.f90,v 1.7 2006/12/09 18:44:58 russ Exp $
module read_netcdf_file

contains

subroutine simple_xy_rd(FILE_NAME, NX, NY, data_in)
  use netcdf
  implicit none

  ! This is the name of the data file we will read. 
  character (len = *), intent(in) :: FILE_NAME

  ! We are reading 2D data, a NX x NY grid. 
  integer, intent(in) :: NX, NY
  real, dimension(:,:), allocatable, intent(out) :: data_in

  ! This will be the netCDF ID for the file and data variable.
  integer :: ncid, varid

  ! Loop indexes, and error handling.
  integer :: x, y

  allocate(data_in(NY, NX))

  ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
  ! the file.
  call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )

  ! Get the varid of the data variable, based on its name.
  call check( nf90_inq_varid(ncid, "data", varid) )

  ! Read the data.
  call check( nf90_get_var(ncid, varid, data_in) )

  ! Check the data.
  !do x = 1, NX
  !   do y = 1, NY
  !      if (data_in(y, x) /= (x - 1) * NY + (y - 1)) then
  !         print *, "data_in(", y, ", ", x, ") = ", data_in(y, x)
  !         stop "Stopped"
  !      end if
  !   end do
  !end do

  ! Close the file, freeing all resources.
  call check( nf90_close(ncid) )

  print *,"*** SUCCESS reading example file ", FILE_NAME, "! "

end subroutine simple_xy_rd

  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      !print *, trim(nf90_strerror(status))
      !print *, nf90_strerror(status)
      stop "Stopped"
    end if
  end subroutine check  


end module read_netcdf_file
