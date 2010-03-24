!    Copyright (C) 2007 Imperial College London and others.
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


#include "confdefs.h"
#include "fdebug.h"

module vertical_extrapolation_module
!!< Module containing routines for vertical extrapolation on
!!< fully unstructured 3D meshes. Also contains routines for updating
!!< distance to top and bottom fields.
use fldebug
use elements
use sparse_tools
use fields
use state_module
use transform_elements
use boundary_conditions
use parallel_tools
use global_parameters, only: new_options, real_4, real_8
use spud
use dynamic_bin_sort_module
use pickers
use coordinates, only: earth_radius
use boundary_conditions
use integer_set_module
use vtk_interfaces
implicit none

interface VerticalExtrapolation
  module procedure VerticalExtrapolationScalar, &
    VerticalExtrapolationMultiple, VerticalExtrapolationVector
end interface VerticalExtrapolation

interface vertical_integration
  module procedure vertical_integration_scalar, &
    vertical_integration_multiple, vertical_integration_vector
end interface vertical_integration

interface verticalshellmapper_find
  module procedure verticalshellmapper_find_sp

  subroutine verticalshellmapper_find(x, y, z, nnodes, &
    & senlist, ntri, fortran_zero_index, &
    & flat_earth, tri_ids, shape_fxn)
    use global_parameters, only: real_8
    implicit none
    integer, intent(in) :: nnodes
    integer, intent(in) :: ntri
    real(kind = real_8), dimension(nnodes), intent(in) :: x
    real(kind = real_8), dimension(nnodes), intent(in) :: y
    real(kind = real_8), dimension(nnodes), intent(in) :: z
    integer, dimension(ntri), intent(in) :: senlist
    integer, intent(in) :: fortran_zero_index
    integer, intent(in) :: flat_earth
    integer, dimension(nnodes), intent(out) :: tri_ids
    real(kind = real_8), dimension(*), intent(out) :: shape_fxn
  end subroutine verticalshellmapper_find
end interface verticalshellmapper_find

!! element i is considered to be above j iff the inner product of the face
!! normal point from i to j and the gravity normal is bigger than this
real, parameter:: VERTICAL_INTEGRATION_EPS=1.0e-8

! the two hemispheres are both projected to the unit circle, but the projected southern
! hemisphere will have its origin shifted. The projected hemispheres are slightly bigger
! than the unit circle, as a bit of overlap with the opposite hemisphere is included around
! the equator. So this shift needs to be big enough to not make both projected hemispheres
! overlap
real, parameter:: HEMISPHERE_SHIFT=10.0
! elements in the overlap that are still included in the projected northern hemisphere are
! elements for which any of its nodes has a z-coordinate bigger than -HEMISPHERE_OVERLAP
! (and z<+HEMISPHERE_OVERLAP for the southern hemisphere)
real, parameter:: HEMISPHERE_OVERLAP=0.01*earth_radius
  
private

public CalculateTopBottomDistance
public VerticalExtrapolation, vertical_integration
public vertical_element_ordering
public VerticalProlongationOperator
public vertical_extrapolation_module_check_options

public :: calculate_hydrostatic_pressure, &
  & subtract_hydrostatic_pressure_gradient
character(len = *), parameter, public :: hp_name = "HydrostaticPressure"

contains

subroutine verticalshellmapper_find_sp(x, y, z, nnodes, &
  & senlist, ntri, fortran_zero_index, &
  & flat_earth, tri_ids, shape_fxn)
  integer, intent(in) :: nnodes
  integer, intent(in) :: ntri
  real(kind = real_4), dimension(nnodes), intent(in) :: x
  real(kind = real_4), dimension(nnodes), intent(in) :: y
  real(kind = real_4), dimension(nnodes), intent(in) :: z
  integer, dimension(ntri), intent(in) :: senlist
  integer, intent(in) :: fortran_zero_index
  integer, intent(in) :: flat_earth
  integer, dimension(nnodes), intent(out) :: tri_ids
  real(kind = real_4), dimension(:), intent(out) :: shape_fxn
  
  real(kind = real_8), dimension(size(shape_fxn)) :: lshape_fxn
  
  call verticalshellmapper_find(real(x, kind = real_8), real(y, kind = real_8), real(z, kind = real_8), nnodes, &
    & senlist, ntri, fortran_zero_index, &
    & flat_earth, tri_ids, lshape_fxn)
  shape_fxn = lshape_fxn
  
end subroutine verticalshellmapper_find_sp
  
subroutine UpdateDistanceField(state, name, vertical_coordinate)
  ! This sub calculates the vertical distance to the free surface
  ! and bottom of the ocean to all nodes. The results are stored
  ! in the 'DistanceToBottom/FreeSurface' fields from state.
  type(state_type), intent(inout):: state
  character(len=*), intent(in):: name
  type(scalar_field), intent(in):: vertical_coordinate

  ! Local variables
  type(vector_field), pointer:: positions, vertical_normal
  type(scalar_field), pointer:: distance
  
  integer, pointer, dimension(:):: surface_element_list

  ! the distance field to compute:
  distance => extract_scalar_field(state, name)
  ! its first boundary condition is on the related top or bottom mesh
  call get_boundary_condition(distance, 1, &
    surface_element_list=surface_element_list)
  positions => extract_vector_field(state, "Coordinate")
  vertical_normal => extract_vector_field(state, "GravityDirection")
  
  ! in each node of the mesh, set "distance" to the vertical coordinate 
  ! of this node projected to the above/below surface mesh
  call VerticalExtrapolation(vertical_coordinate, distance, positions, &
    vertical_normal, surface_element_list, surface_name=name)
    
  ! the distance is then calculated by subtracting its own vertical coordinate
  call addto(distance, vertical_coordinate, scale=-1.0)
  
  if (name=="DistanceToBottom") then
    ! make distance to bottom positive
    call scale(distance, -1.0)
  end if
  
end subroutine UpdateDistanceField
  
subroutine CalculateTopBottomDistance(state)
  !! This sub calculates the vertical distance to the free surface
  !! and bottom of the ocean to all nodes. The results are stored
  !! in the 'DistanceToBottom/Top' fields from state.
  type(state_type), intent(inout):: state
    
  type(mesh_type), pointer:: xmesh
  type(scalar_field):: vertical_coordinate
  
  xmesh => extract_mesh(state, "CoordinateMesh")
  call allocate(vertical_coordinate, xmesh, "VerticalCoordinate")
  call calculate_vertical_coordinate(state, vertical_coordinate)
  call UpdateDistanceField(state, "DistanceToTop", vertical_coordinate)
  call UpdateDistanceField(state, "DistanceToBottom", vertical_coordinate)
  call deallocate(vertical_coordinate)

end subroutine CalculateTopBottomDistance
  
subroutine calculate_vertical_coordinate(state, vertical_coordinate)
  !! Computes a vertical coordinate, i.e. a scalar field such that
  !! for each 2 nodes above each other, the difference of the field
  !! in these nodes gives the distance between them.
  type(state_type), intent(inout):: state
  type(scalar_field), intent(inout):: vertical_coordinate
    
  type(vector_field), pointer:: positions, gravity_normal
  type(scalar_field):: positions_magnitude
  
  positions => extract_vector_field(state, "Coordinate")
  if(have_option('/geometry/spherical_earth')) then
    ! use the radius as vertical coordinate
    ! that is, the l2-norm of the coordinate field
    positions_magnitude=magnitude(positions)
    call set(vertical_coordinate, positions_magnitude)
    call deallocate(positions_magnitude)
  else
    gravity_normal => extract_vector_field(state, "GravityDirection")
    assert(gravity_normal%field_type==FIELD_TYPE_CONSTANT)
    call inner_product(vertical_coordinate, gravity_normal, positions)
    ! gravity points down, we want a vertical coordinate that increases upward
    call scale(vertical_coordinate, -1.0)
  end if
  
end subroutine calculate_vertical_coordinate

subroutine VerticalExtrapolationScalar(from_field, to_field, &
  positions, vertical_normal, surface_element_list, surface_name)
  !!< This sub extrapolates the values on a horizontal 2D surface
  !!< in the vertical direction to 3D fields
  !! The from_fields should be 3D fields of which only the values on the
  !! 2D horizontal surface are used.
  type(scalar_field), intent(in):: from_field
  !! Resulting extrapolated field. May be the same field or a field on
  !! a different mesh (different degree).
  type(scalar_field), intent(inout):: to_field
  !! positions, and upward normal vector on the whole domain
  type(vector_field), target, intent(inout):: positions
  type(vector_field), target, intent(in):: vertical_normal
  !! the surface elements (faces numbers) that make up the surface
  integer, dimension(:), intent(in):: surface_element_list
  !! If provided the projected surface mesh onto horizontal coordinates 
  !! and its associated rtree/pickers are cached under this name and 
  !! attached to the 'positions'. In this case when called again with
  !! the same 'positions' and the same surface_name,
  !! the same surface_element_list should again be provided.
  character(len=*), optional, intent(in):: surface_name

  type(scalar_field), dimension(1):: to_fields
    
  to_fields=(/ to_field /)
  call VerticalExtrapolationMultiple( (/ from_field /) , to_fields, &
      positions, vertical_normal, surface_element_list, &
      surface_name=surface_name)
  
end subroutine VerticalExtrapolationScalar

subroutine VerticalExtrapolationVector(from_field, to_field, &
  positions, vertical_normal, surface_element_list, surface_name)
  !!< This sub extrapolates the values on a horizontal 2D surface
  !!< in the vertical direction to 3D fields
  !! The from_fields should be 3D fields of which only the values on the
  !! 2D horizontal surface are used.
  type(vector_field), intent(in):: from_field
  !! Resulting extrapolated field. May be the same field or a field on
  !! a different mesh (different degree).
  type(vector_field), intent(inout):: to_field
  !! positions, and upward normal vector on the whole domain
  type(vector_field), target, intent(inout):: positions
  type(vector_field), target, intent(in):: vertical_normal
  !! the surface elements (faces numbers) that make up the surface
  integer, dimension(:), intent(in):: surface_element_list
  !! If provided the projected surface mesh onto horizontal coordinates 
  !! and its associated rtree/pickers are cached under this name and 
  !! attached to the 'positions'. In this case when called again with
  !! the same 'positions' and the same surface_name,
  !! the same surface_element_list should again be provided.
  character(len=*), optional, intent(in):: surface_name
  
  type(scalar_field), dimension(from_field%dim):: from_field_components, to_field_components
  integer i
  
  assert(from_field%dim==to_field%dim)
  
  do i=1, from_field%dim
     from_field_components(i)=extract_scalar_field(from_field, i)
     to_field_components(i)=extract_scalar_field(to_field, i)
  end do
    
  call VerticalExtrapolationMultiple( from_field_components, to_field_components, &
        positions, vertical_normal, surface_element_list, &
        surface_name=surface_name)
  
end subroutine VerticalExtrapolationVector

subroutine VerticalExtrapolationMultiple(from_fields, to_fields, &
  positions, vertical_normal, surface_element_list, surface_name)
  !!< This sub extrapolates the values on a horizontal 2D surface
  !!< in the vertical direction to 3D fields
  !! The from_fields should be 3D fields of which only the values on the
  !! 2D horizontal surface are used.
  !! This version takes multiple from_fields at the same time and extrapolates
  !! to to_fields, such that the surface search only has to be done once. This 
  !! will only work if all the from_fields are on the same mesh, and on all the 
  !! to_field are on the same (possibly a different) mesh.
  type(scalar_field), dimension(:), intent(in):: from_fields
  !! Resulting extrapolated field. May be the same field or a field on
  !! a different mesh (different degree).
  type(scalar_field), dimension(:), intent(inout):: to_fields
  !! positions, and upward normal vector on the whole domain
  type(vector_field), target, intent(inout):: positions
  type(vector_field), target, intent(in):: vertical_normal

  !! the surface elements (faces numbers) that make up the surface
  integer, dimension(:), intent(in):: surface_element_list
  !! If provided the projected surface mesh onto horizontal coordinates 
  !! and its associated rtree/pickers are cached under this name and 
  !! attached to the 'positions'. In this case when called again with
  !! the same 'positions' and the same surface_name,
  !! the same surface_element_list should again be provided.
  character(len=*), optional, intent(in):: surface_name

  character(len=FIELD_NAME_LEN):: lsurface_name
  real, dimension(:,:), allocatable:: loc_coords
  integer, dimension(:), allocatable:: seles
  integer i, j, to_nodes, face

  assert(size(from_fields)==size(to_fields))
  assert(element_count(from_fields(1))==element_count(to_fields(1)))
  do i=2, size(to_fields)
     assert(to_fields(1)%mesh==to_fields(i)%mesh)
     assert(from_fields(1)%mesh==from_fields(i)%mesh)
  end do
  
  to_nodes=node_count(to_fields(1))
  ! local coordinates is one more than horizontal coordinate dim
  allocate( seles(to_nodes), loc_coords(1:positions%dim, 1:to_nodes) )
  
  if (present(surface_name)) then
    lsurface_name=surface_name
  else
    lsurface_name="TempSurfaceName"
  end if
  
  ! project the positions of to_fields(1)%mesh into the horizontal plane
  ! and returns face numbers 'seles' and loc_coords to tell where
  ! these projected nodes are found in the surface mesh
  call horizontal_picker(to_fields(1)%mesh, positions, vertical_normal, &
      surface_element_list, lsurface_name, &
      seles, loc_coords)
  
  ! interpolate using the returned faces and loc_coords
  do i=1, size(to_fields)
    do j=1, to_nodes
      face=seles(j)
      call set(to_fields(i), j, &
           dot_product( eval_shape( face_shape(from_fields(i), face), loc_coords(:,j) ), &
             face_val( from_fields(i), face ) ))
     end do
  end do
  
  if (.not. present(surface_name)) then
    call remove_boundary_condition(positions, "TempSurfaceName")
  end if
  
end subroutine VerticalExtrapolationMultiple
  
subroutine horizontal_picker(mesh, positions, vertical_normal, &
  surface_element_list, surface_name, &
  seles, loc_coords)
  !! Searches the nodes of 'mesh' in the surface mesh above.
  !! Returns the surface elements 'seles' that each node lies under
  !! and the loc_coords in this element of this node projected
  !! upward (radially on the sphere) onto the surface mesh
  
  !! mesh
  type(mesh_type), intent(in):: mesh
  !! a valid positions field for the whole domain, not necessarily on 'mesh'
  !! for instance in a periodic domain, mesh is periodic and positions should not be
  type(vector_field), target, intent(inout):: positions
  !! upward normal vector on the whole domain
  type(vector_field), target, intent(in):: vertical_normal
  !! the surface elements (faces numbers) that make up the surface
  integer, dimension(:), intent(in):: surface_element_list
  !! The projected surface mesh onto horizontal coordinates 
  !! and its associated rtree/pickers are cached under this name and 
  !! attached to the 'positions'. When called again with
  !! the same 'positions' and the same surface_name,
  !! the same surface_element_list should again be provided.
  character(len=*), intent(in):: surface_name
  
  !! returned surface elements (face numbers in 'mesh')
  !! and loc coords each node has been found in
  integer, dimension(node_count(mesh)), intent(out):: seles
  real, dimension(positions%dim,node_count(mesh)), intent(out):: loc_coords
  
  type(vector_field):: mesh_positions
  type(vector_field), pointer:: horizontal_positions
  integer, dimension(:), pointer:: horizontal_mesh_list
  real, dimension(:,:), allocatable:: horizontal_coordinate
  real, dimension(vertical_normal%dim):: normal_vector
  real, dimension(positions%dim):: xyz
  integer:: i, stat
  
  assert(.not. mesh_periodic(positions))
  
  if (mesh==positions%mesh) then
    mesh_positions=positions
    ! make mesh_positions indep. ref. of the field, so we can deallocate it 
    ! safely without destroying positions
    call incref(mesh_positions)
  else
    call allocate(mesh_positions, positions%dim, mesh, &
      name='ToPositions_VerticalExtrapolation')
    call remap_field(positions, mesh_positions, stat)
    if (stat/=0 .and. stat/=REMAP_ERR_HIGHER_LOWER_CONTINUOUS .and. &
      stat/=REMAP_ERR_UNPERIODIC_PERIODIC) then
      ! Mapping from higher order to lower order is allowed for coordinates
      ! (well depends on how the higher order is derived from the lower order)
      !
      ! Mapping to periodic coordinates is ok in this case as we only need
      ! locations for the nodes individually (i.e. we don't care about elements
      ! in 'mesh') - the created horizontal_positions will be non-periodic
      ! So that we should be able to find nodes on the periodic boundary on
      ! either side. Using this to interpolate is consistent as long as
      ! the interpolated from field is indeed periodic.
      FLAbort("Unknown error in remmaping positions in horizontal_picker.")
    end if
  end if
  
  ! create an array of the horizontally projected coordinates of the 'mesh' nodes
  allocate( horizontal_coordinate(1:positions%dim-1, 1:node_count(mesh)) )
  if (have_option('/geometry/spherical_earth')) then
    assert(mesh_positions%dim==3)
    do i=1, node_count(mesh)
      xyz=node_val(mesh_positions, i)
      if (xyz(3)>0.0) then
        horizontal_coordinate(:,i) = &
          map2horizontal_sphere( node_val(mesh_positions, i), +1.0)
      else
        horizontal_coordinate(:,i) = &
          map2horizontal_sphere( node_val(mesh_positions, i), -1.0)
        horizontal_coordinate(1,i)=horizontal_coordinate(1,i)+HEMISPHERE_SHIFT
      end if
    end do
  else
    assert( vertical_normal%field_type==FIELD_TYPE_CONSTANT )
    normal_vector=node_val(vertical_normal,1)

    do i=1, node_count(mesh)
      horizontal_coordinate(:,i) = &
          map2horizontal( node_val(mesh_positions, i), normal_vector )
    end do
  end if
    
  call get_horizontal_positions(positions, surface_element_list, &
      vertical_normal, surface_name, &
      horizontal_positions, horizontal_mesh_list)
    
  call picker_inquire( horizontal_positions, horizontal_coordinate, &
    seles, loc_coords, global=.false. )
    
  ! in the spherical case some of the surface elements may be duplicated
  ! within horizontal positions, the returned seles refer to elements
  ! of this horizontal mesh and horiozontal_mesh_list is a map from these
  ! to face numbers in the full mesh
  
  ! we return face numbers however
  do i=1, size(seles)
    if (seles(i)>0) then
      seles(i)=horizontal_mesh_list(seles(i))
    else
      ewrite(0,*) "For node with coordinate", node_val(mesh_positions, i)
      ewrite(0,*) "no top surface node was found."
      FLAbort("Something wrong with the geometry.")
    end if
  end do
    
  call deallocate(mesh_positions)

end subroutine horizontal_picker
  
subroutine get_horizontal_positions(positions, surface_element_list, vertical_normal, surface_name, &
  horizontal_positions, horizontal_mesh_list)
! returns a horizontal positions field over the surface mesh indicated by
! 'surface_element_list'. This field will be created and cached on 'positions'
! as a surface field attached to a dummy boundary condition under the name surface_name
  type(vector_field), intent(inout):: positions
  type(vector_field), intent(in):: vertical_normal
  integer, dimension(:), intent(in):: surface_element_list
  character(len=*), intent(in):: surface_name
  
  ! returns the horizontal positions, and a mapping between elements
  ! in this horizontally projected mesh and facets in 'positions'
  ! in the case of the sphere, this is not the same as surface_element_list
  ! as some of the facets have been duplicated (around the equator)
  type(vector_field), pointer:: horizontal_positions
  integer, dimension(:), pointer:: horizontal_mesh_list
  
  if (.not. has_boundary_condition_name(positions, surface_name)) then
    if (have_option('/geometry/spherical_earth')) then
      call create_horizontal_positions_sphere(positions, &
        surface_element_list, surface_name)
    else
      call create_horizontal_positions_flat(positions, &
        surface_element_list, vertical_normal, surface_name)
    end if
  end if
    
  horizontal_positions => extract_surface_field(positions, surface_name, &
    trim(surface_name)//"HorizontalCoordinate")
    
!   call vtk_write_fields('horizontal_mesh', 0,      horizontal_positions, &
!     horizontal_positions%mesh)    
    
  call get_boundary_condition(positions, surface_name, &
    surface_element_list=horizontal_mesh_list)
  
end subroutine get_horizontal_positions
  
subroutine create_horizontal_positions_flat(positions, surface_element_list, vertical_normal, surface_name)
! adds a "boundary condition" to 'positions' with an associated vector surface field containing a dim-1 
! horizontal coordinate field that can be used to map from the surface mesh specified 
! by 'surface_element_list'. This "boundary condition" will be stored under the name 'surface_name'
! The horizontal coordinates are created by projecting out the component
! in the direction of 'vertical_normal' and then throwing out the x, y or z 
! coordinate that is most aligned with 'vertical_normal'
  type(vector_field), intent(inout):: positions
  type(vector_field), intent(in):: vertical_normal
  integer, dimension(:), intent(in):: surface_element_list
  character(len=*), intent(in):: surface_name
    
  type(mesh_type), pointer:: surface_mesh
  type(vector_field):: horizontal_positions
  integer, dimension(:), pointer:: surface_node_list
  real, dimension(vertical_normal%dim):: normal_vector
  integer:: i, node
  
  assert(vertical_normal%field_type==FIELD_TYPE_CONSTANT)
  normal_vector=node_val(vertical_normal, 1)
  
  call add_boundary_condition_surface_elements(positions, &
    name=surface_name, type="verticalextrapolation", &
    surface_element_list=surface_element_list)
  ! now get back the created surface mesh
  ! and surface_node_list a mapping between node nos in the projected mesh and node nos in the original 'positions' mesh
  call get_boundary_condition(positions, name=surface_name, &
    surface_mesh=surface_mesh, surface_node_list=surface_node_list)
  call allocate( horizontal_positions, dim=positions%dim-1, mesh=surface_mesh, &
    name=trim(surface_name)//"HorizontalCoordinate" )
    
  call insert_surface_field(positions, name=surface_name, &
    surface_field=horizontal_positions)
    
  do i=1, size(surface_node_list)
    node=surface_node_list(i)
    call set(horizontal_positions, i, &
        map2horizontal(node_val(positions, node), normal_vector))
  end do
    
  call deallocate( horizontal_positions )

end subroutine create_horizontal_positions_flat

subroutine create_horizontal_positions_sphere(positions, surface_element_list, surface_name)
! adds a "boundary condition" to 'positions' with an associated vector surface field containing a dim-1 
! horizontal coordinate field that can be used to map from the surface mesh specified 
! by 'surface_element_list'. This "boundary condition" will be stored under the name 'surface_name'
! The horizontal coordinates are created by a stereographic projection from the unit sphere to the
! xy-plane. The northern hemisphere is projected to the unit-circle (including some extra surface
! elements on the southern hemisphere around the equator). The southern hemisphere (including some
! northern surface elements near the equator) is also projected to a unit circle but translated 
! away from the origin to not overlap with the projected northern hemisphere. Thus the created 
! projected horizontal positions field consists of two disjoint areas in the plane corresponding 
! to both hemispheres, and elements in the original surface mesh (near the equator) may appear 
! twice in the projected field.
  type(vector_field), intent(inout):: positions
  integer, dimension(:), intent(in):: surface_element_list
  character(len=*), intent(in):: surface_name
  
  type(vector_field), pointer:: horizontal_positions_north, horizontal_positions_south
  type(vector_field):: horizontal_positions
  type(mesh_type), pointer:: surface_mesh_north, surface_mesh_south
  type(mesh_type):: surface_mesh
  integer, dimension(:), pointer:: surface_element_list_north, surface_element_list_south
  integer, dimension(:), allocatable:: surface_element_list_combined
  integer:: nodes_north, elements_south, elements_north
  integer:: i
  
  ! first create 2 separate horizontal coordinate fields
  ! for each hemisphere
  
  call create_horizontal_positions_hemisphere(positions, &
    surface_element_list, trim(surface_name)//'North', +1.0)
    
  call create_horizontal_positions_hemisphere(positions, &
    surface_element_list, trim(surface_name)//'South', -1.0)
  
  ! these are stored as bcs under the positions
  ! retrieve this information  back:
  
  call get_boundary_condition(positions, name=trim(surface_name)//'North', &
    surface_mesh=surface_mesh_north, &
    surface_element_list=surface_element_list_north)
  call get_boundary_condition(positions, name=trim(surface_name)//'South', &
    surface_mesh=surface_mesh_south, &
    surface_element_list=surface_element_list_south)    
    
  horizontal_positions_north => extract_surface_field(positions, &
      trim(surface_name)//'North', trim(surface_name)//"NorthHorizontalCoordinate")
  horizontal_positions_south => extract_surface_field(positions, &
      trim(surface_name)//'South', trim(surface_name)//"SouthHorizontalCoordinate")    
    
  ! merge these 2 meshes
  surface_mesh=merge_meshes( (/ surface_mesh_north, surface_mesh_south /) )
  
  ! and merge the positions field
  call allocate( horizontal_positions, positions%dim-1, surface_mesh, &
    trim(surface_name)//"HorizontalCoordinate" )
    
  nodes_north=node_count(surface_mesh_north)
  do i=1, nodes_north
    call set( horizontal_positions, i, &
      node_val(horizontal_positions_north, i))
  end do
  ! but translate the projected southern hemisphere to the left
  ! to not overlap it with the northern hemisphere
  do i=1, node_count(surface_mesh_south)
    call set( horizontal_positions, nodes_north+i, &
      node_val(horizontal_positions_south, i)+(/ HEMISPHERE_SHIFT, 0.0 /) )
  end do
    
  ! merge the surface element lists
  ! (note that this is different than the incoming surface_element_list
  ! as it will have some duplicate equatorial elements)
  elements_north=size(surface_element_list_north)
  elements_south=size(surface_element_list_south)
  allocate(surface_element_list_combined(1:elements_north+elements_south))
  surface_element_list_combined(1:elements_north)=surface_element_list_north
  surface_element_list_combined(elements_north+1:)=surface_element_list_south
  
  ! finally the bc for the combined surface mesh:
  ! (note that this creates a new surface mesh different than
  ! the merged mesh, that we won't be using)
  call add_boundary_condition_surface_elements( &
    positions, name=surface_name, type="verticalextrapolation", &
    surface_element_list=surface_element_list_combined)
  
  ! insert the horizontal positions under this bc
  call insert_surface_field( positions, name=surface_name, &
    surface_field=horizontal_positions)
  
  ! everything is safely stored, so we can deallocate our references
  call deallocate(horizontal_positions)
  call deallocate(surface_mesh)
  deallocate(surface_element_list_combined)
  
  ! also we won't need the 2 hemisphere bcs anymore
  call remove_boundary_condition( positions, trim(surface_name)//'North')
  call remove_boundary_condition( positions, trim(surface_name)//'South')
    
end subroutine create_horizontal_positions_sphere

subroutine create_horizontal_positions_hemisphere(positions, &
  surface_element_list, surface_name, hemi_sign)
  type(vector_field), intent(inout):: positions
  integer, dimension(:), intent(in):: surface_element_list
  character(len=*), intent(in):: surface_name
  real, intent(in):: hemi_sign

  type(vector_field):: horizontal_positions
  type(mesh_type), pointer:: surface_mesh
  real, dimension(2):: xy
  integer, dimension(:), pointer:: nodes, surface_node_list
  integer:: i, j, sele, node
  
  type(integer_set):: surface_element_set
  
  call allocate(surface_element_set)
  
  do i=1, size(surface_element_list)
    sele=surface_element_list(i)
    if (any(hemi_sign*face_val(positions, 3, sele)>-HEMISPHERE_OVERLAP)) then
      call insert(surface_element_set, sele)
    end if
  end do
    
  call add_boundary_condition_surface_elements(positions, &
    name=surface_name, type="verticalextrapolation", &
    surface_element_list=set2vector(surface_element_set))
  call deallocate(surface_element_set)
  
  call get_boundary_condition(positions, name=surface_name, &
    surface_mesh=surface_mesh, surface_node_list=surface_node_list)
    
  call allocate( horizontal_positions, dim=2, mesh=surface_mesh, &
    name=trim(surface_name)//"HorizontalCoordinate" )
    
  call insert_surface_field(positions, name=surface_name, &
    surface_field=horizontal_positions)
    
  do i=1, element_count(horizontal_positions)
    nodes => ele_nodes(horizontal_positions, i)
    do j=1, size(nodes)
      node=surface_node_list( nodes(j) )
      xy=map2horizontal_sphere(node_val(positions, node), hemi_sign)
      assert(abs(xy(1))<HEMISPHERE_SHIFT/2.0)
      call set( horizontal_positions, nodes(j), xy)
    end do
  end do
    
  call deallocate(horizontal_positions)
    
end subroutine create_horizontal_positions_hemisphere

function map2horizontal(xyz, normal_vector)
real, dimension(:), intent(in):: xyz, normal_vector
real, dimension(size(xyz)-1):: map2horizontal

  real, dimension(size(xyz)):: hxyz
  integer:: i, c, takeout
  
  ! first subtract of the vertical component
  hxyz=xyz-dot_product(xyz, normal_vector)*normal_vector
  
  ! then leave out the "most vertical" coordinate
  takeout=maxloc(abs(normal_vector), dim=1)
  
  c=1
  do i=1, size(xyz)
    if (i==takeout) cycle
    map2horizontal(c)=hxyz(i)
    c=c+1
  end do

end function map2horizontal

function map2horizontal_sphere(xyz, hemi_sign)
real, dimension(3), intent(in):: xyz
real, intent(in):: hemi_sign
real, dimension(2):: map2horizontal_sphere

  real:: r
  
  r=sqrt(sum(xyz**2))
  map2horizontal_sphere=xyz(1:2)/(r+hemi_sign*xyz(3))
    
end function map2horizontal_sphere

function VerticalProlongationOperator(positions, &
  surface_positions, reduce_columns, owned_nodes)
  !! creates a prolongation operator that prolongates values on
  !! a surface mesh to a full mesh below using interpolation
  !! with VerticalShellMapper. The transpose of this prolongation
  !! operator can be used as a restriction/clustering operator 
  !! from the full mesh to the surface.
  type(csr_matrix) :: VerticalProlongationOperator
  !! the locations of the nodes to which to prolongate, this can be
  !! on any mesh, the nodes are only used a set of loose points
  type(vector_field), target, intent(inout):: positions
  integer::flat_earth_int
  !! positions of the surface mesh - its surface node numbering will 
  !! be used as the numbering of the columns, i.e. of the 
  !! aggregates/clusters. This has to be on a non-periodic mesh to ensure
  !! nodes on the periodic boundary are found
  type(vector_field), intent(in):: surface_positions
  !! if present and true, columns of the prolongation operator
  !! that are not used, i.e. not interpolated from are taken out
  !! and the columns renumbered. This means the column indexing does
  !! no longer correspond to the ordering of surface_positions
  logical, intent(in), optional:: reduce_columns
  !! if present only construct prolongation from first 1:owned_nodes nodes
  !! this can be used in parallel (assuming owned nodes come first) to construct
  !! a local vertical prolongation operator provided the mesh is columnar and decomposed
  !! along columns (2d decomposition). The surface_positions don't necessarily need to be
  !! ordered such that owned nodes are first. In that case you do however need to use
  !! reduce_columns==.true. to remove the inbetween non-owned surface nodes
  integer, intent(in), optional:: owned_nodes

  type(csr_sparsity):: sparsity
  real, dimension(:), allocatable:: x, y, z
  integer, dimension(:), pointer:: sele_nodes
  integer, dimension(:), allocatable:: senlist
  integer, dimension(:), allocatable:: tri_ids
  real, dimension(:), allocatable:: shape_fxn
  integer, dimension(:), allocatable:: surface_node2all
  integer, dimension(:), allocatable:: surface_node2used
  integer, dimension(:), allocatable:: findrm, colm
  real, dimension(:), allocatable:: mat
  real coef
  integer count
  integer i, j, k, nod, snod, entries, nnodes
  logical lreduce_columns
  ! number of vertices in surface element, has to be 3, a triangle
  integer, parameter:: surface_vertices=3
  real, dimension(1:surface_vertices):: loc_coords
  integer, dimension(1:surface_vertices):: flv
  ! coefficient have to be at least this otherwise they're negligable in an interpolation
  real, parameter :: COEF_EPS=1d-10
  
  if(have_option('/geometry/spherical_earth')) then
     flat_earth_int = 0
  else
     flat_earth_int = 1
  end if

  assert(size(local_vertices(ele_shape(surface_positions,1)))==surface_vertices)
  ! local node number of the vertices of a face (surface element)
  !ATTENTION doesn't work with meshes with varying element types
  flv=local_vertices(ele_shape(surface_positions,1))
  
  if (present(reduce_columns)) then
     lreduce_columns=reduce_columns
  else
     lreduce_columns=.false.
  end if
  
  ! the nodes we feed to VerticalShellMapper() are the nodes given
  ! by positions and added to those the vertices of the surface mesh
  ! map between surface node number and this total list of nodes
  allocate( surface_node2all(1:node_count(surface_positions)) )
  surface_node2all=0
  
  ! first in this list are the given positions
  if (present(owned_nodes)) then
    nnodes=owned_nodes
  else
    nnodes=node_count(positions)
  end if
  count=nnodes
  
  ! make linear, triangular surface mesh, regardless of surface mesh degree 
  ! the node numbering refers to all the positions given to VerticalShellMapper
  allocate( senlist(element_count(surface_positions) * surface_vertices) )
 
  entries=1
  do i=1, element_count(surface_positions)
    ! the global nodes of each surface element
    sele_nodes => ele_nodes(surface_positions, i)
    ! copy the vertices thereof in senlist
    do j=1, surface_vertices
      snod=sele_nodes(flv(j))
      if (surface_node2all(snod)==0) then
        ! this vertex is new and added to the list
        count=count+1
        surface_node2all(snod)=count
      end if
      senlist( entries )=surface_node2all(snod)
      entries=entries+1
    end do
  end do
  
  ! Allocate some tempory space. This could be cached.
  allocate(tri_ids(1:count))
  allocate(shape_fxn(1:count*surface_vertices))
  allocate(x(1:count), y(1:count), z(1:count))
  ! first the given nodes:
  x(1:nnodes)=positions%val(X_)%ptr(1:nnodes)
  y(1:nnodes)=positions%val(Y_)%ptr(1:nnodes)
  z(1:nnodes)=positions%val(Z_)%ptr(1:nnodes)
  ! then surface vertices
  do i=1, node_count(surface_positions)
    nod=surface_node2all(i)
    if (nod/=0) then
      x(nod)=surface_positions%val(X_)%ptr(i)
      y(nod)=surface_positions%val(Y_)%ptr(i)
      z(nod)=surface_positions%val(Z_)%ptr(i)
    end if
  end do
  
  ! Perform spatial search
  call VerticalShellMapper_find(x, y, z, count, senlist, &
       element_count(surface_positions), 1, flat_earth_int, tri_ids, shape_fxn)
  deallocate( senlist, surface_node2all )
  
  ! count upper estimate for n/o entries for sparsity
  entries=0
  do nod=1, nnodes
     entries=entries+ele_loc(surface_positions, tri_ids(nod))
  end do
  ! preliminate matrix:  
  allocate( mat(1:entries), findrm(1:nnodes+1), &
    colm(1:entries) )
  
  if (lreduce_columns) then
     ! not all surface nodes may be used (i.e. interpolated from)
     ! unused nodes must be removed and the surface node numbering adjusted
     allocate(surface_node2used(1:node_count(surface_positions)))
     surface_node2used=0
     count=0 ! counts number of used surface nodes
  end if

  entries=0 ! this time only count nonzero entries
  do nod=1, nnodes

     ! beginning of each row in mat
     findrm(nod)=entries+1
     ! local coordinates of the node projected on the surface element
     loc_coords=shape_fxn( (nod-1)*surface_vertices+1 : nod*surface_vertices )
      
     sele_nodes => ele_nodes(surface_positions, tri_ids(nod))
     
     do j=1, size(sele_nodes)
       coef=eval_shape(ele_shape(surface_positions, tri_ids(nod)), j, loc_coords)
       if (abs(coef)>COEF_EPS) then
         snod=sele_nodes(j)
         if (lreduce_columns) then
           if (surface_node2used(snod)==0) then
             ! as of yet unused surface node
             count=count+1
             surface_node2used(snod)=count
           end if
           ! this is the column index we're gonna use instead
           snod=surface_node2used(snod)
         end if
         entries=entries+1
         colm(entries)=snod
         mat(entries)=coef
       end if
     end do
       
  end do
  findrm(nod)=entries+1
  
  if (.not. lreduce_columns) then
    ! we haven't counted used surface nodes, instead we're including all
    count=node_count(surface_positions)
  end if
  
  call allocate(sparsity, nnodes, count, &
     entries, diag=.false., name="VerticalProlongationSparsity")
  sparsity%findrm=findrm
  sparsity%colm=colm(1:entries)
  ! for lots of applications it's good to have sorted rows
  call sparsity_sort(sparsity)
    
  call allocate(VerticalProlongationOperator, sparsity, &
    name="VerticalProlongationOperator")
  call deallocate(sparsity)

  ! as the sparsity has been sorted the ordering of mat(:) no longer
  ! matches that of sparsity%colm, however it still matches the original 
  ! unsorted colm(:)
  do i=1, nnodes
     do k=findrm(i), findrm(i+1)-1
       j=colm(k)
       call set(VerticalProlongationOperator, i, j, mat(k))
     end do
  end do
  
  if (lreduce_columns) then
    deallocate( surface_node2used )
  end if
  deallocate( findrm, colm, mat )
  deallocate( tri_ids, shape_fxn )
  
end function VerticalProlongationOperator
  
subroutine vertical_element_ordering(ordered_elements, face_normal_gravity, optimal_ordering)
!!< Calculates an element ordering such that each element is
!!< is preceded by all elements above it.
integer, dimension(:), intent(out):: ordered_elements
!! need to supply face_normal_gravity matrix,
!! created by compute_face_normal_gravity() subroutine below
type(csr_matrix), intent(in):: face_normal_gravity
!! returns .true. if an optimal ordering is found, i.e there are no 
!! cycles, i.o.w. elements that are (indirectly) above and below each other 
!! at the same time (deck of cards problem).
logical, optional, intent(out):: optimal_ordering
  
  type(dynamic_bin_type) dbin
  real, dimension(:), pointer:: inn
  integer, dimension(:), pointer:: neigh
  integer, dimension(:), allocatable:: bin_list
  integer i, j, elm, bin_no
  logical warning
  
  assert( size(ordered_elements)==size(face_normal_gravity,1) )
  
  ! create binlist, i.e. assign each element to a bin, according to
  ! the number of elements above it
  allocate(bin_list(1:size(ordered_elements)))
  do i=1, size(ordered_elements)
    neigh => row_m_ptr(face_normal_gravity, i)
    inn => row_val_ptr(face_normal_gravity, i)
    ! elements with no element above it go in bin 1
    ! elements with n elements above it go in bin n+1
    ! neigh>0 so we don't count exterior boundary faces
    bin_list(i)=count( inn<-VERTICAL_INTEGRATION_EPS .and. neigh>0 )+1
  end do
    
  call allocate(dbin, bin_list)
  
  warning=.false.
  do i=1, size(ordered_elements)
    ! pull an element from the first non-empty bin
    ! (hopefully an element with no unprocessed elements above)
    call pull_element(dbin, elm, bin_no)
    ordered_elements(i)=elm
    ! if this is bin one then it is indeed an element with no unprocessed 
    !  elements above, otherwise issue a warning
    if (bin_no>1) warning=.true.
    
    ! update elements below:
    
    ! adjacent elements:
    neigh => row_m_ptr(face_normal_gravity, elm)
    inn => row_val_ptr(face_normal_gravity, elm)
    do j=1, size(neigh)
       if (inn(j)>VERTICAL_INTEGRATION_EPS .and. neigh(j)>0) then
         ! element neigh(j) is below element i, therefore now has one
         ! less unprocessed element above it, so can be moved to
         ! lower bin.
         if (.not. element_pulled(dbin, neigh(j))) then
           ! but only if neigh(j) itself hasn't been selected yet
           ! (which might happen for imperfect vertical orderings)
           call move_element(dbin, neigh(j), bin_list(neigh(j))-1)
         end if
       end if
    end do
  end do
  
  if (warning) then
    ! this warning may be reduced (in verbosity level) if it occurs frequently:
    ewrite(0,*) "Warning: vertical_element_ordering has detected a cycle."
    ewrite(0,*) "(deck of cards problem). This may reduce the efficiency"
    ewrite(0,*) "of your vertically sweeping solve."
  end if
  
  if (present(optimal_ordering)) then
    optimal_ordering=.not. warning
  end if
  
  call deallocate(dbin)
  
end subroutine vertical_element_ordering
  
subroutine compute_face_normal_gravity(face_normal_gravity, &
  positions, vertical_normal)
!!< Returns a matrix where A_ij is the inner product of the face normal
!!< and the gravity normal vector of the face between element i and j.
type(csr_matrix), intent(out):: face_normal_gravity
type(vector_field), target, intent(in):: positions, vertical_normal

  type(mesh_type), pointer:: mesh
  real, dimension(:), pointer:: face_normal_gravity_val
  real, dimension(:), allocatable:: detwei_f
  real, dimension(:,:), allocatable:: face_normal, gravity_normal
  integer, dimension(:), pointer:: neigh, faces
  real inn, area
  integer sngi, nloc, i, k
    
  mesh => positions%mesh
  call allocate(face_normal_gravity, mesh%faces%face_list%sparsity)
  call zero(face_normal_gravity)
  
  sngi=face_ngi(mesh, 1)
  nloc=ele_loc(mesh,1)
  allocate( detwei_f(1:sngi), &
    face_normal(1:positions%dim, 1:sngi), &
    gravity_normal(1:positions%dim, 1:sngi))
  
  do i=1, element_count(mesh)
     ! elements adjacent to element i
     ! this is a row (column indeces) in the mesh%faces%face_list matrix
     neigh => ele_neigh(mesh, i)
     ! the surrounding faces
     ! this is a row (integer values) in the mesh%faces%face_list matrix
     faces => ele_faces(mesh, i)
     do k=1, size(neigh)
        if (neigh(k)>i .or. neigh(k)<=0) then
           ! only handling neigh(k)>i to ensure anti-symmetry of the matrix
           ! (and more efficient of course)
           call transform_facet_to_physical(positions, faces(k), &
                detwei_f=detwei_f, &
                normal=face_normal)
           gravity_normal=face_val_at_quad(vertical_normal, faces(k))
           area=sum(detwei_f)
           ! inner product of face normal and vertical normal
           ! integrated over face
           inn=sum(matmul(face_normal*gravity_normal, detwei_f))/area
           if (neigh(k)>0) then
              call set(face_normal_gravity, i, neigh(k), inn)
              call set(face_normal_gravity, neigh(k), i, -inn)
           else
              ! exterior surface: matrix entry does not have valid
              ! column index, still want to store its value, so we
              ! use a pointer
              face_normal_gravity_val => row_val_ptr(face_normal_gravity, i)
              face_normal_gravity_val(k)=inn
           end if
        end if
     end do
        
  end do

end subroutine compute_face_normal_gravity
  
subroutine vertical_integration_scalar(from_field, to_field, &
    positions, vertical_normal, surface_element_list, rhs)
!!< See description vertical_integration_multiple
type(scalar_field), intent(in):: from_field
type(scalar_field), intent(in):: to_field
type(vector_field), intent(in):: positions, vertical_normal
integer, dimension(:), intent(in):: surface_element_list
type(scalar_field), optional, intent(in):: rhs

  type(scalar_field) to_fields(1)
    
  to_fields=(/ to_field /)
  if (present(rhs)) then     
     call vertical_integration_multiple( (/ from_field /), to_fields, &
        positions, vertical_normal, surface_element_list, rhs=(/ rhs /) )
  else  
     call vertical_integration_multiple( (/ from_field /), to_fields, &
        positions, vertical_normal, surface_element_list)
  end if

end subroutine vertical_integration_scalar

subroutine vertical_integration_vector(from_field, to_field, &
    positions, vertical_normal, surface_element_list, rhs)
!!< See description vertical_integration_multiple
type(vector_field), intent(in):: from_field
type(vector_field), intent(in):: to_field
type(vector_field), intent(in):: positions, vertical_normal
integer, dimension(:), intent(in):: surface_element_list
type(vector_field), optional, intent(in):: rhs
  
  type(scalar_field), dimension(from_field%dim):: from_field_components, &
     to_field_components, rhs_components
  integer i
  
  assert(from_field%dim==to_field%dim)
  
  do i=1, from_field%dim
     from_field_components(i)=extract_scalar_field(from_field, i)
     to_field_components(i)=extract_scalar_field(to_field, i)
     if (present(rhs)) then
        rhs_components(i)=extract_scalar_field(rhs, i)
     end if
  end do
  
  if (present(rhs)) then
     call vertical_integration_multiple( from_field_components, &
        to_field_components, positions, vertical_normal, &
        surface_element_list, rhs=rhs_components)
  else
     call vertical_integration_multiple( from_field_components, &
        to_field_components, positions, vertical_normal, &
        surface_element_list)
  end if
  
end subroutine vertical_integration_vector

subroutine calculate_hydrostatic_pressure(state)
  type(state_type), intent(in) :: state
  
  integer :: stat
  integer, dimension(:), pointer :: surface_element_list
  real :: gravity_magnitude
  type(mesh_type) :: from_hp_mesh
  type(mesh_type), pointer :: surface_mesh
  type(scalar_field) :: lbuoyancy, from_hp
  type(scalar_field), pointer :: buoyancy, hp, topdis
  type(vector_field), pointer :: positions, gravity
  
  hp => extract_scalar_field(state, hp_name, stat = stat)
  if(stat /= 0) return
  if(.not. continuity(hp) == -1) then
    FLExit("HydrostaticPressure requires a discontinuous mesh")
  end if
  
  positions => extract_vector_field(state, "Coordinate")
  
  call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
  buoyancy => extract_scalar_field(state, "VelocityBuoyancyDensity")
  assert(ele_count(buoyancy) == ele_count(hp))
  ewrite_minmax(buoyancy%val)
  call allocate(lbuoyancy, buoyancy%mesh, "Buoyancy")
  call set(lbuoyancy, buoyancy)
  call scale(lbuoyancy, gravity_magnitude)
        
  gravity => extract_vector_field(state, "GravityDirection")
  assert(gravity%dim == mesh_dim(hp))
  assert(ele_count(gravity) == ele_count(hp))
  
  topdis => extract_scalar_field(state, "DistanceToTop")
  call get_boundary_condition(topdis, 1, surface_mesh = surface_mesh, surface_element_list = surface_element_list) 
  from_hp_mesh = make_mesh(surface_mesh, shape = face_shape(hp, 1), continuity = -1)
  call allocate(from_hp, from_hp_mesh, hp_name // "BoundaryCondition")
  call deallocate(from_hp_mesh)
  call zero(from_hp)
  
  call vertical_integration(from_hp, hp, positions, gravity, surface_element_list, lbuoyancy)
  
  call deallocate(from_hp)
  call deallocate(lbuoyancy)
  
  ewrite_minmax(hp%val)

end subroutine calculate_hydrostatic_pressure

subroutine subtract_hydrostatic_pressure_gradient(mom_rhs, state)
  !!< Subtract the HydrostaticPressure gradient from the momentum equation
  !!< RHS
  
  type(vector_field), intent(inout) :: mom_rhs
  type(state_type), intent(inout) :: state
  
  integer :: i
  type(vector_field), pointer :: positions
  type(scalar_field), pointer :: hp
  
  ewrite(1, *) "In subtract_hydrostatic_pressure_gradient"
          
  hp => extract_scalar_field(state, hp_name)
          
  ! Apply to momentum equation
  assert(ele_count(hp) == ele_count(mom_rhs))
  
  positions => extract_vector_field(state, "Coordinate")
  assert(positions%dim == mom_rhs%dim)
  assert(ele_count(positions) == ele_count(mom_rhs))

  do i = 1, mom_rhs%dim
    ewrite_minmax(mom_rhs%val(i)%ptr)
  end do
  
  do i = 1, ele_count(mom_rhs)
    call subtract_given_hydrostatic_pressure_gradient_element(i, positions,hp, mom_rhs)
  end do
  
  do i = 1, mom_rhs%dim
    ewrite_minmax(mom_rhs%val(i)%ptr)
  end do
  
  ewrite(1, *) "Exiting subtract_hydrostatic_pressure_gradient"

end subroutine subtract_hydrostatic_pressure_gradient

subroutine subtract_given_hydrostatic_pressure_gradient_element(ele, positions, hp, mom_rhs)
  !!< Subtract the element-wise contribution of the HydrostaticPressure
  !!< gradient from the momentum equation RHS

  integer, intent(in) :: ele
  type(vector_field), intent(in) :: positions
  type(scalar_field), intent(in) :: hp
  type(vector_field), intent(inout) :: mom_rhs
  
  real, dimension(ele_ngi(positions, ele)) :: detwei
  real, dimension(ele_loc(hp, ele), ele_ngi(hp, ele), mom_rhs%dim) :: dn_t
      
  assert(ele_ngi(positions, ele) == ele_ngi(mom_rhs, ele))
  assert(ele_ngi(hp, ele) == ele_ngi(mom_rhs, ele))
      
  call transform_to_physical(positions, ele, ele_shape(hp, ele), &
    & dshape = dn_t, detwei = detwei)
    
  ! /
  ! | -N_A grad gp dV
  ! /
  call addto(mom_rhs, ele_nodes(mom_rhs, ele), -shape_vector_rhs(ele_shape(mom_rhs, ele), transpose(ele_grad_at_quad(hp, ele, dn_t)), detwei))

end subroutine subtract_given_hydrostatic_pressure_gradient_element

subroutine vertical_integration_multiple(from_fields, to_fields, &
    positions, vertical_normal, surface_element_list, rhs)
!!< This subroutine solves: dP/dz=rhs using DG
!!< It can be used for vertical integration downwards (dP/dz=0) as a drop
!!< in replacement of VerticalExtrapolation hence its similar interface.
!!< The field P is the to_field. A boundary condition is given by the 
!!< from_field. Again completely similar to VerticalExtrapolation, it may
!!< be defined as a surface field on surface elements given by 
!!< surface_element_list or it may be a field on the complete mesh
!!< in which case only its values on these surface elements are used.
!!< vertical_normal specifies the direction in which to integrate (usually downwards)
!!< If not specified rhs is assumed zero.
!!<
!!< This version accepts multiple from_fields, to_fields and rhs
type(scalar_field), dimension(:), intent(in):: from_fields
type(scalar_field), dimension(:), intent(inout):: to_fields
type(vector_field), intent(in):: positions, vertical_normal
integer, dimension(:), intent(in):: surface_element_list
type(scalar_field), dimension(:), optional, intent(in):: rhs
  
  type(csr_matrix) face_normal_gravity
  type(element_type), pointer:: ele_shp, face_shp, x_face_shp
  real, dimension(:), pointer:: inn
  real, dimension(:,:,:), allocatable:: surface_rhs, dele_shp
  real, dimension(:,:), allocatable:: ele_mat, face_mat, ele_rhs
  real, dimension(:), allocatable:: detwei, detwei_f
  integer, dimension(:), pointer:: neigh, ele_nds, faces, face_lnds
  integer, dimension(:), allocatable:: ordered_elements, face_nds, face_nds2
  integer nloc, snloc, ngi, sngi
  integer i, j, k, f, f2, elm, it, noit
  logical optimal_ordering, from_surface_fields
  
  assert( size(from_fields)==size(to_fields) )
  
  ! computes inner product of face normal and gravity (see above)
  call compute_face_normal_gravity(face_normal_gravity, &
     positions, vertical_normal)
  
  ! determine an ordering for the elements based on this
  allocate( ordered_elements(1:element_count(positions)) )
  call vertical_element_ordering(ordered_elements, face_normal_gravity, &
     optimal_ordering)
     
  ! General initalisation
  !-----------------------
  ! various grid numbers
  nloc=ele_loc(to_fields(1), 1)
  snloc=face_loc(to_fields(1), 1)
  ngi=ele_ngi(positions, 1)
  sngi=face_ngi(positions, 1)
  ! shape functions
  ele_shp => ele_shape(to_fields(1), 1)
  face_shp => face_shape(to_fields(1), 1)
  x_face_shp => face_shape(positions, 1)
  ! various allocations:
  allocate( &
     surface_rhs(1:snloc, 1:size(to_fields), 1:surface_element_count(positions)), &
     ele_mat(1:nloc, 1:nloc), ele_rhs(1:nloc,1:size(to_fields)), &
     dele_shp(1:nloc, 1:ngi, 1:positions%dim), &
     face_mat(1:snloc, 1:snloc), face_nds(1:snloc), face_nds2(1:snloc), &
     detwei(1:ngi), detwei_f(1:sngi))
     
  if (element_count(from_fields(1))==size(surface_element_list)) then
     ! from_fields are fields over the surface mesh only
     ! so we're using all of its values:
     from_surface_fields=.true.
     ! check the other fields as well:
     do k=2, size(from_fields)
       assert( element_count(from_fields(k))==size(surface_element_list) )
     end do
  else
     ! from_fields are on the full mesh and we only extract its values
     ! at the specified surface_elements
     
     from_surface_fields=.false.     
  end if     
     
  surface_rhs=0
  ! Compute contribution of exterior surface integral (boundary condition) to rhs
  !-----------------------
  do i=1, size(surface_element_list)
     f=surface_element_list(i)
     call transform_facet_to_physical(positions, f, detwei_f)
     face_mat=-shape_shape(face_shp, face_shp, detwei_f)
     do k=1, size(from_fields)
        if (from_surface_fields) then
           ! we need to use ele_val, where i is the element number in the surface_mesh
           surface_rhs(:,k,f)=surface_rhs(:,k,f)+ &
              matmul(face_mat, ele_val(from_fields(k), i))
        else
           ! we can simply use face_val with face number f
           surface_rhs(:,k,f)=surface_rhs(:,k,f)+ &
              matmul(face_mat, face_val(from_fields(k), f))
        end if
     end do
  end do
  
  ! Solution loop
  !-----------------------
  if (optimal_ordering) then
    noit=1
  else
    noit=10
  end if
  
  do it=1, noit
     do i=1, element_count(positions)
       
        elm=ordered_elements(i)
        
        ! construct diagonal matrix block for this element
        call transform_to_physical(positions, elm, &
           shape=ele_shp, dshape=dele_shp, detwei=detwei)
        ele_mat=shape_vector_dot_dshape(ele_shp, &
           ele_val_at_quad(vertical_normal,elm), &
           dele_shp, detwei)

        ! initialise rhs
        if (present(rhs)) then
           do k=1, size(to_fields)
              ele_rhs(:,k)=shape_rhs(ele_shp, detwei*ele_val_at_quad(rhs(k), elm))
           end do
        else
           ele_rhs=0.0
        end if

        
        ! then add contribution of surface integrals of incoming
        ! faces to the rhs and matrix
        neigh => row_m_ptr(face_normal_gravity, elm)
        inn => row_val_ptr(face_normal_gravity, elm)
        faces => ele_faces(positions, elm)
        do j=1, size(neigh)
          if (inn(j)<-VERTICAL_INTEGRATION_EPS) then
            call transform_facet_to_physical(positions, faces(j), &
                 detwei_f)
            face_mat=-shape_shape(face_shp, face_shp, detwei_f)*inn(j)
            
            face_nds=face_global_nodes(to_fields(1), faces(j))
            face_lnds => face_local_nodes(to_fields(1)%mesh, faces(j))
            ele_mat(face_lnds,face_lnds)=ele_mat(face_lnds,face_lnds)+face_mat
            
            if (neigh(j)>0) then
               ! face of element neigh(j), facing elm:
               f2=ele_face(positions, neigh(j), elm)
               face_nds2=face_global_nodes(to_fields(1), f2)
               
               do k=1, size(to_fields)
                  ele_rhs(face_lnds,k)=ele_rhs(face_lnds,k)+ &
                     matmul(face_mat, node_val(to_fields(k), face_nds2))
               end do
            else
               ! note that we've already multiplied with face_mat above, but not with inn(j)
               ele_rhs(face_lnds,:)=ele_rhs(face_lnds,:)+surface_rhs(:,:,faces(j))*inn(j)
            end if
          end if
        end do
        
        call invert(ele_mat)
        
        ! compute values for the to_fields:
        ele_nds => ele_nodes(to_fields(1), elm)
        do k=1, size(to_fields)
           call set( to_fields(k), ele_nds, matmul(ele_mat, ele_rhs(:,k)) )
        end do
     end do
  end do
    
  call deallocate(face_normal_gravity)
  
end subroutine vertical_integration_multiple
  
subroutine vertical_extrapolation_module_check_options
  
  if (have_option("/geometry/ocean_boundaries")) then
    if (.not. have_option("/physical_parameters/gravity")) then
      ewrite(0,*) "If you select /geometry/ocean_boundaries, you also need to &
        &set /physical_parameters/gravity"
      FLAbort("Missing gravity!")
    end if
  end if
  
end subroutine vertical_extrapolation_module_check_options

end module vertical_extrapolation_module
