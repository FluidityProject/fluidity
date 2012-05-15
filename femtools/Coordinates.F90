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
  use FLDebug
  use vector_tools
  use fields
  use global_parameters
  use spud
  use halos
  use halos_base
  use sparse_tools_petsc
  use state_module

  implicit none
  
  private
  
  logical::initialised=.false.
  real, parameter:: rad_to_deg = 180.0/pi
  real, parameter:: deg_to_rad = pi/180.0
  
  public:: &
       LongitudeLatitude,  &
       ll2r3_rotate, rotate2ll, &
       higher_order_sphere_projection, &
       sphere_inward_normal_at_quad_ele, sphere_inward_normal_at_quad_face, &
       rotate_diagonal_to_sphere_gi, rotate_diagonal_to_sphere_face, &
       rotate_ct_m_sphere, rotate_momentum_to_sphere, &
       rotate_velocity_sphere, rotate_velocity_back_sphere, &
       Coordinates_check_options

  interface LongitudeLatitude
     module procedure LongitudeLatitude_single, LongitudeLatitude_multiple
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

  ! rotates vector in cartesian to align with lat/long
  elemental subroutine rotate2ll(longitude, latitude, r3u, r3v, r3w, u, v)
    real, intent(in)  :: longitude, latitude, r3u, r3v, r3w
    real, intent(out) :: u, v
    real lat
    real long
    lat = deg_to_rad*latitude
    long = deg_to_rad*longitude 
    
    u = -(r3u*sin(long)) + r3v*cos(long)
    v = -r3u*cos(long)*sin(lat) - r3v*sin(long)*sin(lat) + r3w*cos(lat)

  end subroutine rotate2ll

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

  function sphere_inward_normal_at_quad_ele(positions, ele_number) result(quad_val)
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

  end function sphere_inward_normal_at_quad_ele

  function sphere_inward_normal_at_quad_face(positions, face_number) result(quad_val)
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

  end function sphere_inward_normal_at_quad_face

  function rotate_diagonal_to_sphere_gi(positions, ele_number, diagonal) result(quad_val)
    ! Given the diagonal of a tensor, this function rotates it to a spherical coordinate system. 
    ! This result is given by R(diagonal)R^T where R is the matrix of Eigen vectors of the 
    ! spherical coordinate system.
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: ele_number
    real, dimension(positions%dim,ele_ngi(positions,ele_number)), intent(in) :: diagonal
    real, dimension(positions%dim,ele_ngi(positions,ele_number)) :: X_quad
    real, dimension(positions%dim,positions%dim) :: R, RT
    real, dimension(positions%dim,positions%dim,ele_ngi(positions,ele_number)) :: diagonal_T, quad_val
    real :: rad, phi, theta
    integer :: i

    assert(positions%dim==3)

    X_quad=ele_val_at_quad(positions, ele_number)

    diagonal_T=0.0
    do i=1,positions%dim
      diagonal_T(i,i,:)=diagonal(i,:)
    end do

    do i=1,ele_ngi(positions,ele_number)
      rad=sqrt(sum(X_quad(:,i)**2))
      phi=atan2(X_quad(2,i),X_quad(1,i))
      theta=acos(X_quad(3,i)/rad)

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
    ! Given the diagonal of a tensor, this function rotates it to a spherical coordinate system. 
    ! This result is given by R(diagonal)R^T where R is the matrix of Eigen vectors of the 
    ! spherical coordinate system.
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: face_number
    real, dimension(positions%dim,face_ngi(positions,face_number)), intent(in) :: diagonal
    real, dimension(positions%dim,face_ngi(positions,face_number)) :: X_quad
    real, dimension(positions%dim,positions%dim) :: R, RT
    real, dimension(positions%dim,positions%dim,face_ngi(positions,face_number)) :: diagonal_T, quad_val
    real :: rad, phi, theta
    integer :: i

    assert(positions%dim==3)

    X_quad=face_val_at_quad(positions, face_number)

    diagonal_T=0.0
    do i=1,positions%dim
      diagonal_T(i,i,:)=diagonal(i,:)
    end do

    do i=1,face_ngi(positions,face_number)
      rad=sqrt(sum(X_quad(:,i)**2))
      phi=atan2(X_quad(2,i),X_quad(1,i))
      theta=acos(X_quad(3,i)/rad)

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
    real :: phi, theta, rad

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

      x=node_val(u_position, node)

      rad=sqrt(sum(x(:)**2))
      phi=atan2(x(2),x(1))
      theta=acos(x(3)/rad)

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
    real :: rad, phi, theta
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

      x=node_val(u_position, node)

      rad=sqrt(sum(x(:)**2))
      phi=atan2(x(2),x(1))
      theta=acos(x(3)/rad)

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
      if (have_option('/geometry/spherical_earth/superparametric_mapping/').and. &
        (.not.have_option("/material_phase["//int2str(m)// &
        "]/scalar_field::Pressure/prognostic/spatial_discretisation/continuous_galerkin/remove_stabilisation_term"))) then
        ewrite(-1,*) "Pressure stabilisation does not currently work with 2nd order or higher coordinate meshes. Please enable"
        ewrite(-1,*) "remove_stabilisation_term under the spatial discretisation tab of your pressure field. Things should work"
        ewrite(-1,*) "nicely then. Thanks!"
        FLExit("Pressure stabilisation is not currently compatible with coordinate fields of order >1.")
      end if
    end do

  end subroutine Coordinates_check_options

end module Coordinates
