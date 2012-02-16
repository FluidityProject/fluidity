! This module loads NEMO data into the states specified in the flml
#include "fdebug.h"
module load_netcdf_module

use global_parameters
use spud
use fields
use coordinates
use Field_Options

implicit none

logical :: on_sphere

contains

subroutine set_scalar_field_from_netcdf(field,path,position)

  type(scalar_field), intent(inout) :: field
  character(len=*), intent(in) :: path
  type(vector_field), intent(in) :: position
  character(len=OPTION_PATH_LEN) :: field_name, format
  real :: gravity_magnitude
  integer :: stat

  ! Are we getting data on a cartesian or lon-lat grid?
  on_sphere = have_option('/geometry/spherical_earth/')

  call load_netcdf_values(field,path,position)

  call get_option(trim(path)//"/from_netcdf/format", format)

  select case (format)
    ! It is possible to manipulate imported netCDF data here.
    ! Note that the 'Free-surface height' format option is no longer supported.
    ! In this case the 'free_surface' node under pressure initial conditions should be used.
    
    case ("raw")
      ewrite(3,*) "The data used to initialise field " // trim(field%name) // &
        ", has been treated as the raw values and not post-processed after import."

    case default
      ewrite(-1,*) "The format " // trim(format) // ", is not a recognised method to handle netCDF files."
      ewrite(-1,*) "Please specify a valid format for initialiasing the field " // trim(field%name) // &
        ", at " // trim(path) // "/from_netcdf/format"
      FLAbort("Fatal error initialising a field from a netCDF file.")
    
  end select

end subroutine

subroutine load_netcdf_values(field,path,position)

  type(scalar_field), intent(inout) :: field
  character(len=*), intent(in) :: path
  type(vector_field), intent(in) :: position
  real, dimension(position%dim,node_count(position)) :: temp_pos
  character(len=FIELD_NAME_LEN) :: filename
  real, dimension(:), allocatable :: X, Y, Z
  integer :: NNodes, i

  assert(node_count(field)==node_count(position))

  call get_option(trim(path)//"/from_netcdf/file_name", filename)

  ewrite(1,*) "Populating field " // trim(field%name) // &
    ", with data from the netCDF file: " // trim(filename)

  NNodes = node_count(field)

  allocate(X(NNodes), Y(NNodes), Z(NNodes))

  if (on_sphere) then
    ! Convert x,y,z coords to lon,lat
    do i=1,NNodes
      temp_pos(:,i)=node_val(position,i)
    end do
    call LongitudeLatitude(temp_pos, X, Y)
  else
    if (position%dim==2) then
      do i=1,NNodes
        temp_pos(:,i)=node_val(position,i)
        X(i)=temp_pos(1,i)
        Y(i)=temp_pos(2,i)
      end do
    else if (position%dim==3) then
      ! Tri-linear interpolation is not currently supported but
      ! this is here for when someone needs it.
      do i=1,NNodes
        temp_pos(:,i)=node_val(position,i)
        X(i)=temp_pos(1,i)
        Y(i)=temp_pos(2,i)
      end do
    else
      FLExit("This dimension is currently not supported")
    end if
  end if

  call get_field_values(trim(filename)//char(0), X, Y, Z, NNodes)

  do i=1,NNodes
     call set(field,i,Z(i))
  enddo

  deallocate(X, Y, Z)

end subroutine

end module load_netcdf_module
