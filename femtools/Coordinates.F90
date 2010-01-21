!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    C.Pain@Imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module Coordinates
  use FLDebug
  use vector_tools
  use fields
  use global_parameters
  implicit none
  
  private
  
  logical::initialised=.false.
  real, parameter:: earth_radius = 6378000
  real, parameter:: rad_to_deg = 180.0/pi
  real, parameter:: deg_to_rad = pi/180.0
  
  public:: &
       LongitudeLatitude,  &
       cart2spher, spher2cart, ll2r3_rotate, rotate2ll, &
       earth_radius, higher_order_sphere_projection
       
  interface LongitudeLatitude
     module procedure LongitudeLatitude_single, LongitudeLatitude_multiple
  end interface

contains
    
  subroutine LongitudeLatitude_single(xyz, longitude, latitude, height)
    real, dimension(:), intent(in):: xyz
    real, intent(out):: longitude, latitude
    real, intent(out), optional:: height
    real r
    
    assert( size(xyz)==3 )
    r = sqrt(sum(xyz**2))
    if(r<1.0) then
       write(0, *) "XYZ = ", xyz
       FLAbort("coordinate doesn't appear to be on earth surface")
    end if

    if(present(height)) then
       height = r - earth_radius
    end if
    longitude = rad_to_deg*atan2(xyz(2), xyz(1))
    latitude = 90.0 - rad_to_deg*acos(xyz(3)/r)
    
  end subroutine LongitudeLatitude_single
  
  subroutine LongitudeLatitude_multiple(xyz, longitude, latitude, height)
    real, dimension(:,:), intent(in):: xyz
    real, dimension(:), intent(out):: longitude, latitude
    real, dimension(:), intent(out), optional::height
    
    integer i
    
    if (present(height)) then
       do i=1, size(xyz,2)
          call LongitudeLatitude_single( xyz(:,i), &
              longitude(i), latitude(i), height(i))
       end do
    else
       do i=1, size(xyz,2)
          call LongitudeLatitude_single( xyz(:,i), &
              longitude(i), latitude(i))
       end do
    end if
    
  end subroutine LongitudeLatitude_multiple
    
  elemental subroutine ll2r3_rotate(longitude, latitude, u, v, r3u, r3v, r3w)
    real, intent(in)::longitude, latitude, u, v
    real, intent(out)::r3u, r3v, r3w
    real t
    
    r3w = v*cos(deg_to_rad*latitude)
    t = v*sin(deg_to_rad*latitude)

    r3v = u*cos(deg_to_rad*longitude) - t*sin(deg_to_rad*longitude)
    r3u = -(u*sin(deg_to_rad*longitude) + t*cos(deg_to_rad*longitude))
    
  end subroutine ll2r3_rotate

  ! rotates vector in cartesian to align with lat/long
  elemental subroutine rotate2ll(longitude, latitude, r3u, r3v, r3w, u, v)
    real, intent(in)  :: longitude, latitude, r3u, r3v, r3w
    real, intent(out) :: u, v
    real lat
    real long
    lat = deg_to_rad*latitude
    long = deg_to_rad*longitude 
    
    u = -(r3u*sin(long)) + r3w*cos(long)
    v = r3u*cos(long)*sin(lat) + r3v*sin(long)*sin(lat) - r3w*cos(lat)

  end subroutine rotate2ll

  subroutine spher2cart(x,y,lat,long,prime_meridian,horiz_rescale)
    real, intent(in)::lat,long,prime_meridian,horiz_rescale
    real, intent(out)::x,y
    
    !  lat/long in degrees assumed at the moment            
    x = earth_radius*cos(lat/rad_to_deg)*(long - prime_meridian)
    y = earth_radius*lat      
    
    ! scale from real units to something that model is using      
    x = x/horiz_rescale
    y = y/horiz_rescale
    
  end subroutine spher2cart
  
  pure subroutine cart2spher(x,y,lat,long,prime_meridian,horiz_rescale,rad)
    real, intent(in)::x,y,prime_meridian,horiz_rescale
    logical, intent(in)::rad
    real, intent(out)::lat,long
    
    real xtmp,ytmp

    xtmp = x*horiz_rescale
    ytmp = y*horiz_rescale

    !     ewrite(3,*) horiz_rescale,x,y

    lat  = ytmp/earth_radius
    long = xtmp/(earth_radius*cos(lat))


    if(rad) then
       long = long + prime_meridian/rad_to_deg
       if(long.lt.0.0) long = 2*pi + long
    else
       ! convert to degrees
       lat  = lat*rad_to_deg
       long = long*rad_to_deg + prime_meridian
       if(long.lt.0.0) long = 360.0 + long
    endif

  end subroutine cart2spher
  
  subroutine higher_order_sphere_projection(positions, s_positions)
    !!< Given a P1 'positions' field and a Pn 's_positions' field, bends the 
    !!< elements of the 's_positions' field onto the sphere
    type(vector_field), intent(inout):: positions
    type(vector_field), intent(inout):: s_positions
    
    real rold, rnew
    integer i
  
    type(scalar_field):: radius, s_radius
    real, dimension(positions%dim):: xyz

    ewrite(1,*), 'In higher_order_sphere_projection'
    
    call allocate(s_radius, s_positions%mesh, "HigherOrderRadius")
    radius=magnitude(positions)
    call remap_field(radius, s_radius)
    
    ! then bend by adjusting to the linearly interpolated radius
    do i=1, node_count(s_positions)
       xyz=node_val(s_positions, i)
       rold=sqrt(sum(xyz**2))
       rnew=node_val(s_radius, i)
       call set(s_positions, i, xyz*rnew/rold)
    end do
    
    call deallocate(s_radius)
    call deallocate(radius)    
  
  end subroutine higher_order_sphere_projection

end module Coordinates
