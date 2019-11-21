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
!    amcgsoftware@imperial.ac.uk
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
  use fldebug
  use iso_c_binding
  use global_parameters
  use futils, only: int2str
  use vector_tools
  use spud
  use parallel_tools, only: isparallel
  use halos_base
  use sparse_tools
  use parallel_fields, only: zero_non_owned
  use fields
  use sparse_tools_petsc
  use state_module
  use halos

  implicit none
  
  private
  
  logical::initialised=.false.
  real, parameter:: rad_to_deg = 180.0/pi
  real, parameter:: deg_to_rad = pi/180.0
  
  public:: &
       LongitudeLatitude,  &
       spherical_polar_2_cartesian, cartesian_2_spherical_polar, &
       spherical_polar_2_cartesian_c, cartesian_2_spherical_polar_c, &
       ll2r3_rotate, &
       lon_lat_height_2_spherical_polar, spherical_polar_2_lon_lat_height, &
       lon_lat_height_2_cartesian, cartesian_2_lon_lat_height, &
       lon_lat_height_2_cartesian_c, cartesian_2_lon_lat_height_c, &
       vector_spherical_polar_2_cartesian, vector_cartesian_2_spherical_polar, &
       vector_lon_lat_height_2_cartesian, vector_cartesian_2_lon_lat_height, &
       vector_lon_lat_height_2_cartesian_c, vector_cartesian_2_lon_lat_height_c, &
       tensor_spherical_polar_2_cartesian, &
       higher_order_sphere_projection, &
       radial_inward_normal_at_quad_ele, radial_inward_normal_at_quad_face, &
       rotate_diagonal_to_sphere_gi, rotate_diagonal_to_sphere_face, &
       rotate_ct_m_sphere, rotate_momentum_to_sphere, &
       rotate_velocity_sphere, rotate_velocity_back_sphere, &
       Coordinates_check_options

  interface LongitudeLatitude
     module procedure LongitudeLatitude_single, LongitudeLatitude_multiple
  end interface

  interface spherical_polar_2_cartesian
     module procedure spherical_polar_2_cartesian, &
                      spherical_polar_2_cartesian_field
  end interface

  interface cartesian_2_spherical_polar
     module procedure cartesian_2_spherical_polar, &
                      cartesian_2_spherical_polar_field
  end interface

  interface vector_spherical_polar_2_cartesian
     module procedure vector_spherical_polar_2_cartesian, &
                      vector_spherical_polar_2_cartesian_field
  end interface

  interface vector_cartesian_2_spherical_polar
     module procedure vector_cartesian_2_spherical_polar, &
                      vector_cartesian_2_spherical_polar_field
  end interface

contains
    
  subroutine LongitudeLatitude_single(xyz, longitude, latitude)
    real, dimension(:), intent(in):: xyz
    real, intent(out):: longitude, latitude
    real r
    
    assert( size(xyz)==3 )
    r = sqrt(sum(xyz**2))
    if(r<1.0) then
       ! May need to include a tolerance here
       write(0, *) "XYZ = ", xyz
       ewrite(-1,*) "Unit vector r on Earth's surface is of size, ", r
       FLAbort("Coordinate doesn't appear to be on the Earth's surface")
    end if

    longitude = rad_to_deg*atan2(xyz(2), xyz(1))
    latitude = 90.0 - rad_to_deg*acos(xyz(3)/r)
    
  end subroutine LongitudeLatitude_single
  
  subroutine LongitudeLatitude_multiple(xyz, longitude, latitude)
    real, dimension(:,:), intent(in):: xyz
    real, dimension(:), intent(out):: longitude, latitude
    
    integer i
    
     do i=1, size(xyz,2)
        call LongitudeLatitude_single( xyz(:,i), &
            longitude(i), latitude(i))
     end do
  
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

  subroutine spherical_polar_2_cartesian(radius,theta,phi,x,y,z)
    !Subroutine for calculation of Cartesian coordinates from spherical-polar
    ! coordinates.
    implicit none

    real, intent(in) :: radius  !Distance from centre of sphere
    real, intent(in) :: theta   !Polar angle, in radians
    real, intent(in) :: phi     !Azimuthal angle, in radians
    real, intent(out) :: x,y,z  !Cartesian coordinates
    
    x = radius*sin(theta)*cos(phi)
    y = radius*sin(theta)*sin(phi)      
    z = radius*cos(theta)
    
  end subroutine spherical_polar_2_cartesian
  
  subroutine spherical_polar_2_cartesian_c(radius,theta,phi,x,y,z) bind(c)
    !C-inter-operable subroutine for calculation of Cartesian coordinates
    ! from spherical-polar coordinates.
    implicit none
    
    real(kind=c_double) :: radius  !Distance from centre of sphere
    real(kind=c_double) :: theta   !Polar angle, in radians
    real(kind=c_double) :: phi     !Azimuthal angle, in radians
    real(kind=c_double) :: x,y,z   !Cartesian coordinates

    real :: radius_f
    real :: theta_f
    real :: phi_f
    real :: x_f,y_f,z_f

    !Cast input variables to Fortran intrinsic types.
    radius_f = real(radius)
    theta_f = real(theta)
    phi_f = real(phi)

    !Convert coordinates
    call spherical_polar_2_cartesian(radius_f,theta_f,phi_f,x_f,y_f,z_f)

    !Cast output variables to C-inter-operable types.
    x = real(x_f, kind=c_double)
    y = real(y_f, kind=c_double)
    z = real(z_f, kind=c_double)

  end subroutine spherical_polar_2_cartesian_c

  subroutine cartesian_2_spherical_polar(x,y,z,radius,theta,phi)
    !Subroutine for calculation of spherical-polar coordinates from cartesian.
    implicit none

    real, intent(in) :: x,y,z   !cartesian coordinates
    real, intent(out) :: radius !Distance from centre of sphere
    real, intent(out) :: theta  !Polar angle, in radians
    real, intent(out) :: phi    !Azimuthal angle, in radians

    radius = sqrt(x**2 + y**2 + z**2)
    theta = acos(z/radius)
    phi = atan2(y,x)

  end subroutine cartesian_2_spherical_polar

  subroutine cartesian_2_spherical_polar_c(x, y, z, radius, theta, phi) bind(c)
    !C-inter-operable subroutine for calculation of spherical-polar coordinates
    ! from Cartesian coordinates.
    implicit none
    
    real(kind=c_double) :: x,y,z   !cartesian coordinates
    real(kind=c_double) :: radius  !Distance from centre of sphere
    real(kind=c_double) :: theta   !Polar angle, in radians
    real(kind=c_double) :: phi     !Azimuthal angle, in radians

    real :: x_f,y_f,z_f
    real :: radius_f
    real :: theta_f
    real :: phi_f

    !Cast input variables to fortran intrinsic types.
    x_f = real(x)
    y_f = real(y)
    z_f = real(z)

    !Convert coordinates
    call cartesian_2_spherical_polar(x_f, y_f, z_f, radius_f, theta_f, phi_f)

    !Cast output variables to C-inter-operable types.
    radius = real(radius_f, kind=c_double)
    theta = real(theta_f, kind=c_double)
    phi = real(phi_f, kind=c_double)

  end subroutine cartesian_2_spherical_polar_c
  
  subroutine spherical_polar_2_cartesian_field(spherical_polar_coordinate_field, &
                                               cartesian_coordinate_field)
    !Subroutine for conversion of a spherical-polar coordinate field into a cartesian
    ! coordinate field.
    implicit none

    type(vector_field) :: spherical_polar_coordinate_field
    type(vector_field) :: cartesian_coordinate_field
    integer :: node
    real, dimension(3) :: XYZ, RTP !arrays containing a single node's position vector
                                   ! in cartesian & spherical-polar bases 

    do node=1,node_count(spherical_polar_coordinate_field)
      RTP = node_val(spherical_polar_coordinate_field, node)
      call spherical_polar_2_cartesian(RTP(1), RTP(2), RTP(3), XYZ(1), XYZ(2), XYZ(3))
      call set(cartesian_coordinate_field, node, XYZ)
    enddo

  end subroutine spherical_polar_2_cartesian_field

  subroutine cartesian_2_spherical_polar_field(cartesian_coordinate_field, &
                                               spherical_polar_coordinate_field)
    !Subroutine for conversion of a cartesian coordinate field into a spherical-polar
    ! coordinate field.
    implicit none

    type(vector_field) :: cartesian_coordinate_field
    type(vector_field) :: spherical_polar_coordinate_field
    integer :: node
    real, dimension(3) :: XYZ, RTP !arrays containing a single node's position vector
                                   ! components in cartesian & spherical-polar bases 

    do node=1,node_count(cartesian_coordinate_field)
      XYZ = node_val(cartesian_coordinate_field, node)
      call cartesian_2_spherical_polar(XYZ(1), XYZ(2), XYZ(3), RTP(1), RTP(2), RTP(3))
      call set(spherical_polar_coordinate_field, node, RTP)
    enddo

  end subroutine cartesian_2_spherical_polar_field

  subroutine lon_lat_height_2_spherical_polar(longitude, latitude, height, &
                                              radius, theta, phi, &
                                              referenceRadius)
    !Subroutine for conversion of longitude-latitude-height coordinates on a 
    !  sphere to spherical-polar coordinates. Longitude and latitude must be
    !  in degrees, polar coordinates are returned into radians
    implicit none

    real, intent(in) :: longitude !in degrees
    real, intent(in) :: latitude  !in degrees
    real, intent(in) :: height
    real, intent(out) :: radius !Distance from centre of sphere
    real, intent(out) :: theta  !Polar angle, in radians
    real, intent(out) :: phi    !Azimuthal angle, in radians
    real, intent(in), optional :: referenceRadius !Distance form the centre of
                                                  ! the sphere to its surface
    real :: pi

    pi=4*atan(1.0)

    !Convert longitude to azimuthal angle and latitude in polar angle; in radians.
    phi = longitude*pi/180.
    theta = (90.- latitude)*pi/180.

    !Convert height to distance from origin
    ! Check if referenceRadius is present. If not use default value
    ! of surface radius, available in global_parameters module
    if(present(referenceRadius)) then
      radius = height + referenceRadius
    else
      radius = height + surface_radius
    endif

  end subroutine lon_lat_height_2_spherical_polar

  subroutine spherical_polar_2_lon_lat_height(radius, theta, phi, &
                                              longitude, latitude, height, &
                                              referenceRadius)
    !Subroutine for conversion of spherical-polar coordinates to
    !  longitude-latitude-height coordinates. The polar coordinates must
    !  be given in radians. Longitude and latitude are returned in
    !  degrees. If referenceRadius is specified, height is measured as the
    !  radial distance relative to that radius, ie it is the distance relative to the
    !  surface of the sphere. if referenceRadius is absent height is the distance
    !  from the center of the sphere.
    implicit none

    real, intent(in) :: radius !Distance from centre of sphere
    real, intent(in) :: theta  !Polar angle, in radians
    real, intent(in) :: phi    !Azimuthal angle, in radians
    real, intent(out) :: longitude !in degrees
    real, intent(out) :: latitude  !in degrees
    real, intent(out) :: height
    real, intent(in), optional :: referenceRadius !distance form the centre of
                                                  ! the sphere to its surface
    real :: pi

    pi=4*atan(1.0)

    longitude = phi*180.0/pi
    latitude = (pi/2 - theta)*180.0/pi

    !If referenceRadius is present, subtract it from the radial distance
    if(present(referenceRadius)) then
      height = radius - referenceRadius
    else
      height = radius - surface_radius
    endif

  end subroutine spherical_polar_2_lon_lat_height

  subroutine lon_lat_height_2_cartesian(longitude, latitude, height, &
                                        x, y, z, &
                                        referenceRadius)
    !Subroutine for conversion of longitude-latitude-height coordinates into 
    ! Cartesian coordinates. If referenceRadius is specified, height is measured
    ! as the radial distance relative to that radius, i.e. it is the distance
    ! relative to the surface of the sphere.
    implicit none

    real, intent(in) :: longitude !in degrees
    real, intent(in) :: latitude  !in degrees
    real, intent(in) :: height
    real, intent(out) :: x,y,z   !Cartesian coordinates
    real, intent(in), optional :: referenceRadius

    real :: radius !Distance from centre of sphere
    real :: theta  !Polar angle, in radians
    real :: phi    !Azimuthal angle, in radians

    !Convert longitude-latitude-height into spherical-polar coordinates.
    ! Check if referenceRadius is present. If not use default value
    ! of surface radius, available in global_parameters module
    if(present(referenceRadius)) then
      call lon_lat_height_2_spherical_polar(longitude, latitude, height, &
                                            radius, theta, phi, &
                                            referenceRadius)
    else
      call lon_lat_height_2_spherical_polar(longitude, latitude, height, &
                                            radius, theta, phi, &
                                            surface_radius)
    endif


    !convert spherical-polar coordinates to Cartesian
    call spherical_polar_2_cartesian(radius,theta,phi,x,y,z)

  end subroutine lon_lat_height_2_cartesian

  subroutine lon_lat_height_2_cartesian_c(longitude, latitude, height, &
                                          x, y, z, &
                                          referenceRadius) bind(c)
    !C-inter-operable subroutine for conversion of longitude-latitude-height into
    ! spherical-polar coordinates. referenceRadius must be specified, i.e. height
    ! is always measured as the radial distance relative to that radius and denotes
    ! the distance from the surface of the sphere.
    implicit none
    
    real(kind=c_double) :: longitude        !Longitude, in radians.
    real(kind=c_double) :: latitude         !Latitude, in radians.
    real(kind=c_double) :: height           !Distance from surface of sphere.
    real(kind=c_double) :: x,y,z            !Cartesian coordinates.
    real(kind=c_double) :: referenceRadius  !Sphere radius.

    real :: longitude_f
    real :: latitude_f
    real :: height_f
    real :: x_f,y_f,z_f
    real :: referenceRadius_f

    !Cast input variables to Fortran intrinsic types.
    longitude_f = real(longitude)
    latitude_f = real(latitude)
    height_f = real(height)
    referenceRadius_f = real(referenceRadius)

    !Convert coordinates
    call lon_lat_height_2_cartesian(longitude_f, latitude_f, height_f, &
                                    x_f, y_f, z_f, &
                                    referenceRadius_f)

    !Cast output variables to C-inter-operable types.
    x = real(x_f, kind=c_double)
    y = real(y_f, kind=c_double)
    z = real(z_f, kind=c_double)

  end subroutine lon_lat_height_2_cartesian_c

  subroutine cartesian_2_lon_lat_height(x, y, z, longitude, latitude, height, &
                                        referenceRadius)
    !Subroutine for conversion of Cartesian coordinates into longitude-latitude-height
    ! If referenceRadius is specified, height is measures as the radial distance relative
    ! to that radius.
    implicit none

    real, intent(in) :: x,y,z   !Cartesian coordinates
    real, intent(out) :: longitude !in degrees
    real, intent(out) :: latitude  !in degrees
    real, intent(out) :: height
    real, intent(in), optional :: referenceRadius
    real :: radius !Distance from centre of sphere
    real :: theta  !Polar angle, in radians
    real :: phi    !Azimuthal angle, in radians

    !convert Cartesian coordinates to spherical-polar
    call cartesian_2_spherical_polar(x,y,z,radius,theta,phi)

    !Convert polar angle into latitude and azimuthal angle into longitude; in radians.
    if(present(referenceRadius)) then
      call spherical_polar_2_lon_lat_height(radius, theta, phi, &
                                            longitude, latitude, height, &
                                            referenceRadius)
    else
      call spherical_polar_2_lon_lat_height(radius, theta, phi, &
                                            longitude, latitude, height)
    endif


  end subroutine cartesian_2_lon_lat_height

  subroutine cartesian_2_lon_lat_height_c(x, y, z, longitude, latitude, height, &
                                          referenceRadius) bind(c)
    !C-inter-operable subroutine for conversion of Cartesian coordinates into
    ! longitude-latitude-height.
    implicit none

    real(kind=c_double) :: x,y,z   !Cartesian coordinates
    real(kind=c_double) :: longitude !in degrees
    real(kind=c_double) :: latitude  !in degrees
    real(kind=c_double) :: height
    real(kind=c_double) :: referenceRadius

    real :: x_f,y_f,z_f
    real :: longitude_f
    real :: latitude_f
    real :: height_f
    real :: referenceRadius_f

    !Cast input variables to Fortran intrinsic types.
    x_f = real(x)
    y_f = real(y)
    z_f = real(z)

    referenceRadius_f = real(referenceRadius)

    !Convert coordinates
    call cartesian_2_lon_lat_height(x_f, y_f, z_f, longitude_f, latitude_f, height_f, &
                                        referenceRadius_f)

    !Cast output variables to C-inter-operable types.
    longitude = real(longitude_f, kind=c_double)
    latitude = real(latitude_f, kind=c_double)
    height = real(height_f, kind=c_double)

  end subroutine cartesian_2_lon_lat_height_c

  subroutine transformation_matrix_cartesian_2_spherical_polar(xCoord, yCoord, zCoord, R, RT)
    !Subroutine calculating transformation matrix for spherical-polar to/from Cartesian
    ! tensor transformations. The routine also returns the transposed transformation matrix
    implicit none

    real, intent(in) :: xCoord  !x-component of position vector
    real, intent(in) :: yCoord  !y-component of position vector
    real, intent(in) :: zCoord  !z-component of position vector
    real, dimension(3,3), intent(out) :: R   !Transformation matrix
    real, dimension(3,3), intent(out) :: RT  !Transposed transformation matrix

    real :: radius        !Distance from centre of sphere
    real :: theta         !Polar angle, in radians
    real :: phi           !Azimuthal angle, in radians

    !Calculate position-vector components in spherical-polar basis
    call cartesian_2_spherical_polar(xCoord, yCoord, zCoord, radius, theta, phi)

    R(1,1)=sin(theta)*cos(phi)
    R(1,2)=sin(theta)*sin(phi)
    R(1,3)=cos(theta)
    R(2,1)=cos(theta)*cos(phi)
    R(2,2)=cos(theta)*sin(phi)
    R(2,3)=-sin(theta)
    R(3,1)=-sin(phi)
    R(3,2)=cos(phi)
    R(3,3)=0.0

    RT = TRANSPOSE(R)

  end subroutine transformation_matrix_cartesian_2_spherical_polar

  subroutine vector_spherical_polar_2_cartesian(radial, polar, azimuthal, &
                                                radius, theta, phi, &
                                                xComp, yComp, zComp, &
                                                xCoord, yCoord, zCoord)
    !Subroutine for vector change of basis: from spherical-polar to cartesian. The
    ! coordinates of the position vector are also transformed
    implicit none

    real, intent(in) :: radial        !Radial component of vector
    real, intent(in) :: polar         !Polar component of vector
    real, intent(in) :: azimuthal     !Azimuthal  component of vector
    real, intent(in) :: radius        !Distance from centre of sphere
    real, intent(in) :: theta         !Polar angle, in radians
    real, intent(in) :: phi           !Azimuthal angle, in radians
    real, intent(out) :: xComp        !1st vector component in cartesian basis
    real, intent(out) :: yComp        !2nd vector component in cartesian basis
    real, intent(out) :: zComp        !3rd vector component in cartesian basis
    real, intent(out) :: xCoord       !1st vector component of position vector in cartesian basis
    real, intent(out) :: yCoord       !2nd vector component of position vector in cartesian basis
    real, intent(out) :: zCoord       !3rd vector component of position vector in cartesian basis

    real, dimension(3) :: cartesianComponents
    real, dimension(3,3) :: R   !Transformation matrix
    real, dimension(3,3) :: RT  !Transposed transformation matrix

    !Calculate position-vector components in cartesian system
    call spherical_polar_2_cartesian(radius, theta, phi, xCoord, yCoord, zCoord)

    !Calculate transformation matrix
    call transformation_matrix_cartesian_2_spherical_polar(xCoord, yCoord, zCoord, R, RT)

    !Evaluate vector components in Cartesian basis
    cartesianComponents = matmul(RT,(/radial, polar, azimuthal/))
    xComp = cartesianComponents(1)
    yComp = cartesianComponents(2)
    zComp = cartesianComponents(3)

  end subroutine vector_spherical_polar_2_cartesian

  subroutine vector_cartesian_2_spherical_polar(xComp, yComp, zComp, &
                                                xCoord, yCoord, zCoord, &
                                                radial, polar, azimuthal, &
                                                radius, theta, phi)
    !Subroutine for vector change of basis: from Cartesian to spherical-polar. The
    ! coordinates of the position vector are also transformed
    implicit none

    real, intent(in) :: xComp         !1st vector component in cartesian basis
    real, intent(in) :: yComp         !2nd vector component in cartesian basis
    real, intent(in) :: zComp         !3rd vector component in cartesian basis
    real, intent(in) :: xCoord        !1st vector component of position vector in cartesian basis
    real, intent(in) :: yCoord        !2nd vector component of position vector in cartesian basis
    real, intent(in) :: zCoord        !3rd vector component of position vector in cartesian basis
    real, intent(out) :: radial       !Radial component of vector
    real, intent(out) :: polar        !Polar component of vector
    real, intent(out) :: azimuthal    !Azimuthal  component of vector
    real, intent(out) :: radius       !Distance from centre of sphere
    real, intent(out) :: theta        !Polar angle, in radians
    real, intent(out) :: phi          !Azimuthal angle, in radians

    real, dimension(3) :: sphericalPolarComponents
    real, dimension(3,3) :: R   !Transformation matrix
    real, dimension(3,3) :: RT  !Transposed transformation matrix

    !Calculate position-vector components in spherical-polar system
    call cartesian_2_spherical_polar(xCoord, yCoord, zCoord, radius, theta, phi)

    !Calculate transformation matrix
    call transformation_matrix_cartesian_2_spherical_polar(xCoord, yCoord, zCoord, R, RT)

    !Evaluate vector components in spherical-polar basis
    sphericalPolarComponents = matmul(R,(/xComp, yComp, zComp/))
    radial = sphericalPolarComponents(1)
    polar = sphericalPolarComponents(2)
    azimuthal = sphericalPolarComponents(3)

  end subroutine vector_cartesian_2_spherical_polar

  subroutine vector_lon_lat_height_2_cartesian(zonalComponent,&
                                               meridionalComponent,&
                                               verticalComponent, &
                                               longitude, &
                                               latitude, &
                                               height, &
                                               xComp, yComp, zComp, &
                                               xCoord, yCoord, zCoord, &
                                               referenceRadius)
    !Subroutine for change of basis of a vector from meridional-zonal-vertical
    !  components to cartesian components.
    implicit none

    real, intent(in) :: zonalComponent      !Vector component tangential to parallel
    real, intent(in) :: meridionalComponent !Vector component tangential to meridian
    real, intent(in) :: verticalComponent   !Vecor component in the vertical (radial)
    real, intent(in) :: longitude 
    real, intent(in) :: latitude
    real, intent(in) :: height
    real, intent(out) :: xComp          !1st vector component in cartesian basis
    real, intent(out) :: yComp          !2nd vector component in cartesian basis
    real, intent(out) :: zComp          !3rd vector component in cartesian basis
    real, intent(out) :: xCoord         !1st vector component of position vector
                                        ! in Cartesian basis
    real, intent(out) :: yCoord         !2nd vector component of position vector
                                        ! in Cartesian basis
    real, intent(out) :: zCoord         !3rd vector component of position vector
                                        ! in Cartesian basis
    real, intent(in), optional :: referenceRadius
    real :: radial       !Radial component of vector
    real :: polar        !Polar component of vector
    real :: azimuthal    !Azimuthal  component of vector
    real :: radius       !Distance from centre of sphere
    real :: theta        !Polar angle, in radians
    real :: phi          !Azimuthal angle, in radians

    !Convert zonal-meridional-vertical components to spherical-polar
    azimuthal = zonalComponent
    polar = -meridionalComponent
    radial = verticalComponent
    !Convert longitude-latitude-height to spherical-polar.
    ! If referenceRadius is present then pass that to coordinate conversion routine,
    ! height then is the radial distance of a point from the sphere with radius=
    ! referenceRadius. Otherwise height is simply the distance from the Cartesian
    ! coordinate origin.
    if(present(referenceRadius)) then
      call lon_lat_height_2_spherical_polar(longitude, latitude, height, &
                                            radius, theta, phi, &
                                            referenceRadius)
    else
      call lon_lat_height_2_spherical_polar(longitude, latitude, height, &
                                            radius, theta, phi)
    endif
    !convert spherical-polar components to cartesian.
    call vector_spherical_polar_2_cartesian(radial, polar, azimuthal, &
                                            radius, theta, phi, &
                                            xComp, yComp, zComp, &
                                            xCoord, yCoord, zCoord)

  end subroutine vector_lon_lat_height_2_cartesian

  subroutine vector_cartesian_2_lon_lat_height(xComp, yComp, zComp, &
                                               xCoord, yCoord, zCoord, &
                                               zonalComponent,&
                                               meridionalComponent,&
                                               verticalComponent, &
                                               longitude, &
                                               latitude, &
                                               height, &
                                               referenceRadius)
    !Subroutine for change of basis of a vector from cartesian to
    !  meridional-zonal-vertical.
    implicit none

    real, intent(in) :: xComp          !1st vector component in cartesian basis
    real, intent(in) :: yComp          !2nd vector component in cartesian basis
    real, intent(in) :: zComp          !3rd vector component in cartesian basis
    real, intent(in) :: xCoord         !1st vector component of position vector
                                       ! in Cartesian basis
    real, intent(in) :: yCoord         !2nd vector component of position vector
                                       ! in Cartesian basis
    real, intent(in) :: zCoord         !3rd vector component of position vector
                                       ! in Cartesian basis
    real, intent(out) :: zonalComponent      !Vector component tangential to parallel
    real, intent(out) :: meridionalComponent !Vector component tangential to meridian
    real, intent(out) :: verticalComponent   !Vector component in the vertical (radial)
    real, intent(out) :: longitude 
    real, intent(out) :: latitude
    real, intent(out) :: height
    real, intent(in), optional :: referenceRadius
    real :: radial       !Radial component of vector
    real :: polar        !Polar component of vector
    real :: azimuthal    !Azimuthal  component of vector
    real :: radius       !Distance from centre of sphere
    real :: theta        !Polar angle, in radians
    real :: phi          !Azimuthal angle, in radians

    !Convert cartesian components to spherical-polar
    call vector_cartesian_2_spherical_polar(xComp, yComp, zComp, &
                                            xCoord, yCoord, zCoord, &
                                            radial, polar, azimuthal, &
                                            radius, theta, phi)
    !Convert cartesian coordinates to longitude-latitude-radius
    if(present(referenceRadius)) then
      call cartesian_2_lon_lat_height(xCoord, yCoord, zCoord, &
                                      longitude, latitude, height, &
                                      referenceRadius)
    else
      call cartesian_2_lon_lat_height(xCoord, yCoord, zCoord, &
                                      longitude, latitude, height)
    endif
    !Convert spherical-polar components to zonal-meridional-vertical
    zonalComponent = azimuthal
    meridionalComponent = -polar
    verticalComponent = radial

  end subroutine vector_cartesian_2_lon_lat_height

  subroutine vector_lon_lat_height_2_cartesian_c(zonalComponent,&
                                                 meridionalComponent,&
                                                 verticalComponent, &
                                                 longitude, &
                                                 latitude, &
                                                 height, &
                                                 xComp, yComp, zComp, &
                                                 xCoord, yCoord, zCoord, &
                                                 referenceRadius) bind(c)
    !C-interoperable subroutine for change of basis of a vector from
    !  meridional-zonal-vertical components to cartesian components. Note that
    !  unlike the Fortran version of the present routine, referenceRadius is
    !  a mandatory argument.
    implicit none

    real(kind=c_double), intent(in) :: zonalComponent      !Vector component tangential
                                                           ! to parallel
    real(kind=c_double), intent(in) :: meridionalComponent !Vector component tangential
                                                           ! to meridian
    real(kind=c_double), intent(in) :: verticalComponent   !Vecor component in the
                                                           ! vertical (radial)
    real(kind=c_double), intent(in) :: longitude 
    real(kind=c_double), intent(in) :: latitude
    real(kind=c_double), intent(in) :: height
    real(kind=c_double), intent(out) :: xComp      !1st vector component in
                                                   ! cartesian basis
    real(kind=c_double), intent(out) :: yComp      !2nd vector component in
                                                   ! cartesian basis
    real(kind=c_double), intent(out) :: zComp      !3rd vector component in
                                                   ! cartesian basis
    real(kind=c_double), intent(out) :: xCoord     !1st vector component of
                                                   ! position vector in cartesian basis
    real(kind=c_double), intent(out) :: yCoord     !2nd vector component of
                                                   ! position vector in cartesian basis
    real(kind=c_double), intent(out) :: zCoord     !3rd vector component of
                                                   ! position vector in cartesian basis
    real(kind=c_double), intent(in) :: referenceRadius

    real :: zonalComponent_f      !Vector component tangential to parallel
    real :: meridionalComponent_f !Vector component tangential to meridian
    real :: verticalComponent_f   !Vecor component in the vertical (radial)
    real :: longitude_f 
    real :: latitude_f
    real :: height_f
    real :: xComp_f          !1st vector component in cartesian basis
    real :: yComp_f          !2nd vector component in cartesian basis
    real :: zComp_f          !3rd vector component in cartesian basis
    real :: xCoord_f         !1st vector component of position vector in cartesian basis
    real :: yCoord_f         !2nd vector component of position vector in cartesian basis
    real :: zCoord_f         !3rd vector component of position vector in cartesian basis
    real :: referenceRadius_f

    !Convert C-types in to Fortran intrinsic types.
    zonalComponent_f = real(zonalComponent)
    meridionalComponent_f = real(meridionalComponent)
    verticalComponent_f = real(verticalComponent)
    longitude_f = real(longitude)
    latitude_f = real(latitude)
    height_f = real(height)
    referenceRadius_f = real(referenceRadius)

    !Convert coordinates and components.
    call vector_lon_lat_height_2_cartesian(zonalComponent_f,&
                                           meridionalComponent_f,&
                                           verticalComponent_f, &
                                           longitude_f, &
                                           latitude_f, &
                                           height_f, &
                                           xComp_f, yComp_f, zComp_f, &
                                           xCoord_f, yCoord_f, zCoord_f, &
                                           referenceRadius_f)

    !Convert Fortran intrinsic types to C-types.
    xComp = real(xComp_f, kind=c_double)
    yComp = real(yComp_f, kind=c_double)
    zComp = real(zComp_f, kind=c_double)
    xCoord = real(xCoord_f, kind=c_double)
    yCoord = real(yCoord_f, kind=c_double)
    zCoord = real(zCoord_f, kind=c_double)

  end subroutine vector_lon_lat_height_2_cartesian_c

  subroutine vector_cartesian_2_lon_lat_height_c(xComp, yComp, zComp, &
                                                 xCoord, yCoord, zCoord, &
                                                 zonalComponent,&
                                                 meridionalComponent,&
                                                 verticalComponent, &
                                                 longitude, &
                                                 latitude, &
                                                 height, &
                                                 referenceRadius) bind (c)
    !C inter-operable subroutine for change of basis of a vector from Cartesian to
    !  meridional-zonal-vertical. Note that
    !  unlike the Fortran version of the present routine, referenceRadius is
    !  a mandatory argument.
    implicit none

    real(kind=c_double), intent(in) :: xComp          !1st vector component in
                                                      ! cartesian basis
    real(kind=c_double), intent(in) :: yComp          !2nd vector component in
                                                      ! cartesian basis
    real(kind=c_double), intent(in) :: zComp          !3rd vector component in
                                                      ! cartesian basis
    real(kind=c_double), intent(in) :: xCoord         !1st vector component of position
                                                      ! vector in Cartesian basis
    real(kind=c_double), intent(in) :: yCoord         !2nd vector component of position
                                                      ! vector in Cartesian basis
    real(kind=c_double), intent(in) :: zCoord         !3rd vector component of position
                                                      ! vector in Cartesian basis
    real(kind=c_double), intent(out) :: zonalComponent      !Vector component tangential
                                                            ! to parallel
    real(kind=c_double), intent(out) :: meridionalComponent !Vector component tangential
                                                            ! to meridian
    real(kind=c_double), intent(out) :: verticalComponent   !Vector component in the
                                                            ! vertical (radial)
    real(kind=c_double), intent(out) :: longitude
    real(kind=c_double), intent(out) :: latitude
    real(kind=c_double), intent(out) :: height
    real(kind=c_double), intent(in) :: referenceRadius

    real :: xComp_f          !1st vector component in cartesian basis
    real :: yComp_f          !2nd vector component in cartesian basis
    real :: zComp_f          !3rd vector component in cartesian basis
    real :: xCoord_f         !1st vector component of position vector
                                         ! in Cartesian basis
    real :: yCoord_f         !2nd vector component of position vector
                                         ! in Cartesian basis
    real :: zCoord_f         !3rd vector component of position vector
                                         ! in Cartesian basis
    real :: zonalComponent_f      !Vector component tangential to parallel
    real :: meridionalComponent_f !Vector component tangential to meridian
    real :: verticalComponent_f   !Vector component in the vertical (radial)
    real :: longitude_f
    real :: latitude_f
    real :: height_f
    real :: referenceRadius_f

    !Convert C-types in to Fortran intrinsic types.
    xComp_f = real(xComp)
    yComp_f = real(yComp)
    zComp_f = real(zComp)
    xCoord_f = real(xCoord)
    yCoord_f = real(yCoord)
    zCoord_f = real(zCoord)
    referenceRadius_f = real(referenceRadius)

    !Convert coordinates and components.
    call vector_cartesian_2_lon_lat_height(xComp_f, yComp_f, zComp_f, &
                                           xCoord_f, yCoord_f, zCoord_f, &
                                           zonalComponent_f, &
                                           meridionalComponent_f, &
                                           verticalComponent_f, &
                                           longitude_f, &
                                           latitude_f, &
                                           height_f, &
                                           referenceRadius_f)

    !Convert Fortran intrinsic types to C-types.
    zonalComponent = real(zonalComponent_f, kind=c_double)
    meridionalComponent = real(meridionalComponent_f, kind=c_double)
    verticalComponent = real(verticalComponent_f, kind=c_double)
    longitude = real(longitude_f, kind=c_double)
    latitude = real(latitude_f, kind=c_double)
    height = real(height_f, kind=c_double)

  end subroutine vector_cartesian_2_lon_lat_height_c

  subroutine vector_spherical_polar_2_cartesian_field(spherical_polar_vector_field, &
                                                      spherical_polar_coordinate_field, &
                                                      cartesian_vector_field, &
                                                      cartesian_coordinate_field)
    !Subroutine for change of basis of a cartesian vector field into a spherical-polar
    ! vector field. This routine also converts and returns the position vector component
    ! fields
    implicit none

    type(vector_field) :: spherical_polar_vector_field
    type(vector_field) :: spherical_polar_coordinate_field
    type(vector_field) :: cartesian_vector_field
    type(vector_field) :: cartesian_coordinate_field
    integer :: node
    real, dimension(3) :: XYZ, RTP !arrays containing a signel node's position vector
                                   ! in cartesian & spherical-polar bases 
    real, dimension(3) :: cartesianComponents, sphericalPolarComponents

    assert(node_count(spherical_polar_coordinate_field) == node_count(cartesian_coordinate_field))

    do node=1,node_count(spherical_polar_coordinate_field)
      RTP = node_val(spherical_polar_coordinate_field, node)
      sphericalPolarComponents = node_val(spherical_polar_vector_field, node)
      call vector_spherical_polar_2_cartesian(sphericalPolarComponents(1), &
                                              sphericalPolarComponents(2), &
                                              sphericalPolarComponents(3), &
                                              RTP(1), RTP(2), RTP(3), &
                                              cartesianComponents(1), &
                                              cartesianComponents(2), &
                                              cartesianComponents(3), &
                                              XYZ(1), XYZ(2), XYZ(3))
      call set(cartesian_coordinate_field, node, XYZ)
      call set(cartesian_vector_field, node, cartesianComponents)
    enddo
  end subroutine vector_spherical_polar_2_cartesian_field

  subroutine vector_cartesian_2_spherical_polar_field(cartesian_vector_field, &
                                                      cartesian_coordinate_field, &
                                                      spherical_polar_vector_field, &
                                                      spherical_polar_coordinate_field)
    !Subroutine for change of basis of a cartesian vector field into a spherical-polar
    ! vector field. This routine also converts and returns the position vector component
    ! fields
    implicit none

    type(vector_field) :: cartesian_vector_field
    type(vector_field) :: cartesian_coordinate_field
    type(vector_field) :: spherical_polar_vector_field
    type(vector_field) :: spherical_polar_coordinate_field
    integer :: node
    real, dimension(3) :: XYZ, RTP !arrays containing a signel node's position vector
                                   ! in cartesian & spherical-polar bases 
    real, dimension(3) :: cartesianComponents, sphericalPolarComponents

    assert(node_count(spherical_polar_coordinate_field) == node_count(cartesian_coordinate_field) )

    do node=1,node_count(spherical_polar_coordinate_field)
      XYZ = node_val(cartesian_coordinate_field, node)
      cartesianComponents = node_val(cartesian_vector_field, node)
      call vector_cartesian_2_spherical_polar(cartesianComponents(1), &
                                              cartesianComponents(2), &
                                              cartesianComponents(3), &
                                              XYZ(1), XYZ(2), XYZ(3), & 
                                              sphericalPolarComponents(1), &
                                              sphericalPolarComponents(2), &
                                              sphericalPolarComponents(3), &
                                              RTP(1), RTP(2), RTP(3))
      call set(spherical_polar_coordinate_field, node, RTP)
      call set(spherical_polar_vector_field, node, sphericalPolarComponents)
    enddo
  end subroutine vector_cartesian_2_spherical_polar_field

  subroutine tensor_spherical_polar_2_cartesian(sphericalPolarComponents, &
                                                radius, theta, phi, &
                                                cartesianComponents, &
                                                xCoord, yCoord, zCoord)
    !Subroutine for tensor change of basis: From spherical-polar to cartesian. The
    ! coordinates of the position vector are also transformed. The tensor must
    ! be a 3x3 tensor.
    implicit none

    real, intent(in), dimension(3,3) :: sphericalPolarComponents   !Tensor 
                                      ! components in spherical-polar basis
    real, intent(in) :: radius        !Distance from centre of sphere
    real, intent(in) :: theta         !Polar angle, in radians
    real, intent(in) :: phi           !Azimuthal angle, in radians
    real, intent(out), dimension(3,3) :: cartesianComponents       !Tensor
                                      ! components in Cartesian bisis
    real, intent(out) :: xCoord       !1st vector component of position vector
                                      ! in cartesian basis
    real, intent(out) :: yCoord       !2nd vector component of position vector
                                      ! in cartesian basis
    real, intent(out) :: zCoord       !3rd vector component of position vector
                                      ! in cartesian basis

    real, dimension(3,3) :: R   !Transformation matrix
    real, dimension(3,3) :: RT  !Transposed transformation matrix

    !Calculate position-vector components in cartesian system
    call spherical_polar_2_cartesian(radius, theta, phi, xCoord, yCoord, zCoord)

    !Calculate transformation matrix
    call transformation_matrix_cartesian_2_spherical_polar(xCoord, yCoord, zCoord, R, RT)

    !Evaluate vector components in Cartesian basis
    cartesianComponents = matmul(matmul(RT, sphericalPolarComponents), R)

  end subroutine tensor_spherical_polar_2_cartesian

  subroutine higher_order_sphere_projection(positions, s_positions)
    !!< Given a P1 'positions' field and a Pn 's_positions' field, bends the 
    !!< elements of the 's_positions' field onto the sphere
    type(vector_field), intent(inout):: positions
    type(vector_field), intent(inout):: s_positions
    
    real rold, rnew
    integer i
  
    type(scalar_field):: radius, s_radius
    real, dimension(positions%dim):: xyz

    ewrite(1,*) 'In higher_order_sphere_projection'
    
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

  function radial_inward_normal_at_quad_ele(positions, ele_number) result(quad_val)
    ! Return the direction of gravity at the quadrature points of and element.
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: ele_number
    real, dimension(positions%dim,ele_ngi(positions,ele_number)) :: X_quad, quad_val
    integer :: i,j

    X_quad=ele_val_at_quad(positions, ele_number)

    do j=1,ele_ngi(positions,ele_number)
      do i=1,positions%dim
        quad_val(i,j)=-X_quad(i,j)/sqrt(sum(X_quad(:,j)**2))
      end do
    end do

  end function radial_inward_normal_at_quad_ele

  function radial_inward_normal_at_quad_face(positions, face_number) result(quad_val)
    ! Return the direction of gravity at the quadrature points of and element.
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: face_number
    real, dimension(positions%dim,face_ngi(positions,face_number)) :: X_quad, quad_val
    integer :: i,j

    X_quad=face_val_at_quad(positions, face_number)

    do j=1,face_ngi(positions,face_number)
      do i=1,positions%dim
        quad_val(i,j)=-X_quad(i,j)/sqrt(sum(X_quad(:,j)**2))
      end do
    end do

  end function radial_inward_normal_at_quad_face

  function rotate_diagonal_to_sphere_gi(positions, ele_number, diagonal) result(quad_val)
    ! Given the diagonal of a tensor in cartesian coordinates, this function
    ! transforms the tensor components to a spherical-polar basis. This result
    ! is given by R(diagonal)R^T where R is the matrix of Eigen vectors of the
    ! spherical-polar basis, expressed in the cartesian basis.
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: ele_number
    real, dimension(positions%dim,ele_ngi(positions,ele_number)), intent(in) :: diagonal
    real, dimension(positions%dim,ele_ngi(positions,ele_number)) :: X_quad
    real, dimension(positions%dim,positions%dim) :: R, RT
    real, dimension(positions%dim,positions%dim,ele_ngi(positions,ele_number)) :: diagonal_T, quad_val
    real :: radius, theta, phi !distance form origin, polar angle, azimuthal angle
    integer :: i

    assert(positions%dim==3)

    X_quad=ele_val_at_quad(positions, ele_number)

    diagonal_T=0.0
    do i=1,positions%dim
      diagonal_T(i,i,:)=diagonal(i,:)
    end do

    do i=1,ele_ngi(positions,ele_number)
      ! Calculate the spherical-polar coordinates of the point
      call cartesian_2_spherical_polar(X_quad(1,i), X_quad(2,i), X_quad(3,i), radius, theta, phi)

      R(1,1)=-sin(phi)
      R(1,2)=cos(theta)*cos(phi)
      R(1,3)=sin(theta)*cos(phi)
      R(2,1)=cos(phi)
      R(2,2)=cos(theta)*sin(phi)
      R(2,3)=sin(theta)*sin(phi)
      R(3,1)=0
      R(3,2)=-sin(theta)
      R(3,3)=cos(theta)

      RT=R
      call invert(RT)
      quad_val(:,:,i)=matmul((matmul(R,diagonal_T(:,:,i))),RT)

    end do

  end function rotate_diagonal_to_sphere_gi

  function rotate_diagonal_to_sphere_face(positions, face_number, diagonal) result(quad_val)
    ! Given the diagonal of a tensor in cartesian coordinates, this function
    ! transforms the tensor components to a spherical-polar basis. This result
    ! is given by R(diagonal)R^T ! where R is the matrix of Eigen vectors of the
    ! spherical-polar basis, expressed in the cartesian basis.
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: face_number
    real, dimension(positions%dim,face_ngi(positions,face_number)), intent(in) :: diagonal
    real, dimension(positions%dim,face_ngi(positions,face_number)) :: X_quad
    real, dimension(positions%dim,positions%dim) :: R, RT
    real, dimension(positions%dim,positions%dim,face_ngi(positions,face_number)) :: diagonal_T, quad_val
    real :: radius, theta, phi !distance form origin, polar angle, azimuthal angle
    integer :: i

    assert(positions%dim==3)

    X_quad=face_val_at_quad(positions, face_number)

    diagonal_T=0.0
    do i=1,positions%dim
      diagonal_T(i,i,:)=diagonal(i,:)
    end do

    do i=1,face_ngi(positions,face_number)
      ! Calculate the spherical-polar coordinates of the point
      call cartesian_2_spherical_polar(X_quad(1,i), X_quad(2,i), X_quad(3,i), radius, theta, phi)

      R(1,1)=-sin(phi)
      R(1,2)=cos(theta)*cos(phi)
      R(1,3)=sin(theta)*cos(phi)
      R(2,1)=cos(phi)
      R(2,2)=cos(theta)*sin(phi)
      R(2,3)=sin(theta)*sin(phi)
      R(3,1)=0
      R(3,2)=-sin(theta)
      R(3,3)=cos(theta)

      RT=R
      call invert(RT)
      quad_val(:,:,i)=matmul((matmul(R,diagonal_T(:,:,i))),RT)

    end do

  end function rotate_diagonal_to_sphere_face

  subroutine rotate_ct_m_sphere(state, ct_m, u)

    type(block_csr_matrix), intent(inout):: ct_m
    type(vector_field), intent(in) :: u

    type(vector_field) :: sphere_normal, sphere_tangent1, sphere_tangent2
    integer, dimension(:), pointer:: rowcol
    real, dimension(u%dim, u%dim):: local_rotation
    real, dimension(u%dim):: ct_xyz, ct_rot
    real, dimension(:), pointer:: rowval
    integer:: node, i, j, k, rotated_node
    
    type(state_type), intent(in) :: state
    type(vector_field), pointer :: position
    type(vector_field) :: u_position
    real, dimension(u%dim) :: x, node_normal, node_tangent1, node_tangent2
    real :: radius, theta, phi !distance form origin, polar angle, azimuthal angle

    ewrite(1,*) "Inside rotate_ct_m_sphere"

    assert( all(blocks(ct_m) == (/ 1, u%dim /)) )

    position => extract_vector_field(state, "Coordinate")
    call allocate(u_position, u%dim, u%mesh, name="VelocityCoordinate")
    call remap_field(position, u_position)

    if (associated(u%mesh%halos)) then
      call halo_update(u_position)
    end if

    assert(u%dim==3)

    call allocate(sphere_normal, u%dim, u%mesh, name="sphere_normal")
    call allocate(sphere_tangent1, u%dim, u%mesh, name="sphere_tangent1")
    call allocate(sphere_tangent2, u%dim, u%mesh, name="sphere_tangent2")

    do node=1, node_count(u)

      !Extract the cartesian coordinates of the node.
      x=node_val(u_position, node)

      !Calculate spherical-polar coordinates.
      call cartesian_2_spherical_polar(x(1),x(2),x(3),radius,theta,phi)

      node_normal=(/sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)/)
      node_tangent1=(/-sin(phi),cos(phi),0.0/)
      node_tangent2=(/cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)/)

      call set(sphere_normal, node, node_normal)
      call set(sphere_tangent1, node, node_tangent1)
      call set(sphere_tangent2, node, node_tangent2)

    end do

    if (associated(u%mesh%halos)) then
      call halo_update(sphere_normal)
      call halo_update(sphere_tangent1)
      call halo_update(sphere_tangent2)
    end if

    do i=1, size(ct_m, 1)
      rowcol => row_m_ptr(ct_m, i)
      do j=1, size(rowcol)
        rotated_node=rowcol(j)
        ! construct local rotation matrix
        local_rotation(1,:)=node_val(sphere_tangent1, rotated_node)
        local_rotation(2,:)=node_val(sphere_tangent2, rotated_node)
        local_rotation(3,:)=node_val(sphere_normal, rotated_node)

        ! look up ct_m values of row i, column rowcol(j) in xyz orientation
        do k=1, blocks(ct_m,2)
          rowval => row_val_ptr(ct_m, 1, k, i)
          ct_xyz(k)=rowval(j)
        end do
        ! rotate to tangent1, tangent2, normal orientation
        ct_rot=matmul( local_rotation, ct_xyz)
        ! put back in the matrix
        do k=1, blocks(ct_m,2)
          rowval => row_val_ptr(ct_m, 1, k, i)
          rowval(j)=ct_rot(k)
        end do
      end do
    end do

    call deallocate(u_position)
    call deallocate(sphere_normal)
    call deallocate(sphere_tangent1)
    call deallocate(sphere_tangent2)

  end subroutine rotate_ct_m_sphere

  subroutine rotate_momentum_to_sphere(big_m, rhs, u, state, dg)

    type(petsc_csr_matrix), intent(inout):: big_m
    type(vector_field), intent(inout):: rhs
    type(vector_field), intent(inout):: u
    type(state_type), intent(inout):: state
    logical, intent(in) :: dg

    type(petsc_csr_matrix), pointer:: rotation_sphere
    type(petsc_csr_matrix):: rotated_big_m
    type(vector_field):: result
    integer :: stat

    ewrite(1,*) "Inside rotate_momentum_to_sphere"

    rotation_sphere => extract_petsc_csr_matrix(state, "RotationMatrixSphere", stat=stat)

    if (stat/=0) then
      allocate(rotation_sphere)
      call create_rotation_matrix_sphere(rotation_sphere, u, state)
      call insert(state, rotation_sphere, "RotationMatrixSphere")
    end if

    ! rotate big_m:
    call ptap(rotated_big_m, big_m, rotation_sphere)

    ! rotate rhs:
    ! need to have separate copy of the field, because of intent(out) and intent(in)
    ! of mult_T call, as result%val points at the same space as rhs%val, this directly
    ! puts the result in rhs as well 
    result=rhs 
    call mult_T(result, rotation_sphere, rhs)
    if (dg) then
      ! We have just poluted the halo rows of the rhs. This is incorrect
      ! in the dg case due to the non-local assembly system employed.
      call zero_non_owned(rhs)
    end if
    ! rotate u:
    if (dg) then
      call zero_non_owned(u)
    end if
    result=u ! same story
    call mult_T(result, rotation_sphere, u)

    ! throw out unrotated big_m and replace with rotated:
    call deallocate(big_m)
    big_m=rotated_big_m

    if (stat/=0) then
      call deallocate(rotation_sphere)
      deallocate(rotation_sphere)
    end if

  end subroutine rotate_momentum_to_sphere

  subroutine create_rotation_matrix_sphere(rotation_sphere, u, state)

    type(petsc_csr_matrix), intent(out):: rotation_sphere
    type(vector_field), intent(in):: u
    type(state_type), intent(in) :: state

    type(halo_type), pointer:: halo
    type(vector_field) :: sphere_normal, sphere_tangent1, sphere_tangent2
    real, dimension(u%dim) :: x, node_normal, node_tangent1, node_tangent2
    real :: radius, theta, phi !distance form origin, polar angle, azimuthal angle
    real, dimension(u%dim, u%dim):: local_rotation
    integer, dimension(:), allocatable:: dnnz, onnz
    integer:: node, nodes, mynodes
    logical:: parallel

    type(vector_field), pointer :: position
    type(vector_field) :: u_position    

    ewrite(1,*) "Inside create_rotation_matrix_sphere"

    nodes=node_count(u)
    if (associated(u%mesh%halos)) then
       halo => u%mesh%halos(1)
       mynodes=halo_nowned_nodes(halo)
    else
       nullify(halo)
       mynodes=nodes
    end if
    parallel=IsParallel()

    allocate(dnnz(1:mynodes*u%dim), onnz(1:mynodes*u%dim))
    onnz=0
    ! default is just a 1.0 on the diagonal (no rotation)
    dnnz=1

    do node=1, mynodes
      if (any(dnnz(node:node+(u%dim-1)*mynodes:mynodes)>1)) then
        FLExit("Two rotated specifications for the same node.")
      end if
      dnnz( node:node+(u%dim-1)*mynodes:mynodes)=u%dim
    end do

    call allocate(rotation_sphere, nodes, nodes, &
         dnnz, onnz, (/ u%dim, u%dim /), "RotationMatrixSphere", halo=halo)

    position => extract_vector_field(state, "Coordinate")
    call allocate(u_position, u%dim, u%mesh, name="VelocityCoordinate")
    call remap_field(position, u_position)

    call allocate(sphere_normal, u%dim, u%mesh, name="sphere_normal")
    call allocate(sphere_tangent1, u%dim, u%mesh, name="sphere_tangent1")
    call allocate(sphere_tangent2, u%dim, u%mesh, name="sphere_tangent2")

    do node=1, mynodes

      !Extract the cartesian coordinates of the node.
      x=node_val(u_position, node)

      !Calculate spherical-polar coordinates.
      call cartesian_2_spherical_polar(x(1),x(2),x(3),radius,theta,phi)

      node_normal=(/sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)/)
      node_tangent1=(/-sin(phi),cos(phi),0.0/)
      node_tangent2=(/cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)/)

      call set(sphere_normal, node, node_normal)
      call set(sphere_tangent1, node, node_tangent1)
      call set(sphere_tangent2, node, node_tangent2)

    end do

    do node=1, mynodes
      local_rotation(:,1)=node_val(sphere_tangent1, node)
      local_rotation(:,2)=node_val(sphere_tangent2, node)
      local_rotation(:,3)=node_val(sphere_normal, node)

      call addto(rotation_sphere, node, node, local_rotation)
    end do

    call assemble(rotation_sphere)

    call deallocate(u_position)
    call deallocate(sphere_normal)
    call deallocate(sphere_tangent1)
    call deallocate(sphere_tangent2)

  end subroutine create_rotation_matrix_sphere

  subroutine rotate_velocity_sphere(vfield, state)

    type(vector_field), intent(inout):: vfield
    type(state_type), intent(inout):: state
    
    type(vector_field), pointer:: u
    type(vector_field):: result
    type(petsc_csr_matrix), pointer:: rotation_sphere
    integer :: stat
    
    rotation_sphere => extract_petsc_csr_matrix(state, "RotationMatrixSphere", stat=stat)
    if (stat/=0) then
      allocate(rotation_sphere)
      u => extract_vector_field(state, "Velocity")
      call create_rotation_matrix_sphere(rotation_sphere, u, state)
      call insert(state, rotation_sphere, "RotationMatrixSphere")
    end if
    
    result=vfield ! see note in rotate_momentum_equation
    call mult_T(result, rotation_sphere, vfield)

    if (stat/=0) then
      call deallocate(rotation_sphere)
      deallocate(rotation_sphere)
    end if

  end subroutine rotate_velocity_sphere
  
  subroutine rotate_velocity_back_sphere(vfield, state)

    type(vector_field), intent(inout):: vfield
    type(state_type), intent(inout):: state
    
    type(vector_field), pointer:: u
    type(vector_field):: result
    type(petsc_csr_matrix), pointer:: rotation_sphere
    integer :: stat
    
    rotation_sphere => extract_petsc_csr_matrix(state, "RotationMatrixSphere", stat=stat)
    if (stat/=0) then
      allocate(rotation_sphere)
      u => extract_vector_field(state, "Velocity")
      call create_rotation_matrix_sphere(rotation_sphere, u, state)
      call insert(state, rotation_sphere, "RotationMatrixSphere")
    end if
    
    result=vfield ! see note in rotate_momentum_equation
    call mult(result, rotation_sphere, vfield)

    if (stat/=0) then
      call deallocate(rotation_sphere)
      deallocate(rotation_sphere)
    end if

  end subroutine rotate_velocity_back_sphere

  ! Coordinates options checking
  subroutine Coordinates_check_options

    integer :: nmat, m

    ! Pressure stabilisation does not currently work with a p2 or higher
    ! coordinate fields. Check that this term is not enabled and, if it is,
    ! exit.
    nmat = option_count("/material_phase")
    do m = 0, nmat-1
      if (have_option('/geometry/spherical_earth/superparametric_mapping/') &
          .and. have_option("/material_phase["//int2str(m)//"]/scalar_field::Pressure/"// &
                            "prognostic/spatial_discretisation/continuous_galerkin") &
          .and. .not.have_option("/material_phase["//int2str(m)// &
                                 "]/scalar_field::Pressure/prognostic/spatial_discretisation/"// &
                                 "continuous_galerkin/remove_stabilisation_term")) then
        ewrite(-1,*) "Pressure stabilisation does not currently work with 2nd order or higher coordinate meshes. Please enable"
        ewrite(-1,*) "remove_stabilisation_term under the spatial discretisation tab of your pressure field. Things should work"
        ewrite(-1,*) "nicely then. Thanks!"
        FLExit("Pressure stabilisation is not currently compatible with coordinate fields of order >1.")
      end if
    end do

  end subroutine Coordinates_check_options

end module Coordinates
