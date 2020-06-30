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


#include "confdefs.h"
#include "fdebug.h"

module vertical_extrapolation_module
!!< Module containing routines for vertical extrapolation on
!!< fully unstructured 3D meshes. Also contains routines for updating
!!< distance to top and bottom fields.
use fldebug
use global_parameters, only: FIELD_NAME_LEN
use elements
use parallel_tools
use spud
use data_structures
use sparse_tools
use transform_elements
use fetools
use parallel_fields
use fields
use state_module
use boundary_conditions
use dynamic_bin_sort_module
use pickers
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

!! element i is considered to be above j iff the inner product of the face
!! normal point from i to j and the gravity normal is bigger than this
real, parameter:: VERTICAL_INTEGRATION_EPS=1.0e-8

! the two hemispheres are both projected to the unit circle, but the projected southern
! hemisphere will have its origin shifted. The projected hemispheres are slightly bigger
! than the unit circle, as a bit of overlap with the opposite hemisphere is included around
! the equator. So this shift needs to be big enough to not make both projected hemispheres
! overlap
real, parameter:: HEMISPHERE_SHIFT=10.0
  
private

public CalculateTopBottomDistance
public VerticalExtrapolation, vertical_integration
public vertical_element_ordering
public VerticalProlongationOperator
public vertical_extrapolation_module_check_options
public compute_face_normal_gravity

! with /geometry/spherical_earth the vertical extrapolation is done using a gnomonic
! projection of each of the 3*dim sides of the cubed sphere. The gnomonic projection projects
! each of the sides to the [-1,1]x[-1,1] square (dim=3 case) - but we allow some overlap so that
! facets on the surface mesh are represented in at least one (and possibly multiple) sides of the
! cubed sphere projection. This is done by including any facet that maps entirely
! within [-w,w]x[-w,w] where w=GNOMONIC_BOX_WIDTH. The 3*dim images of the projected sides of the
! cubed sphere are then shifted (shifting the origin with multiples of 3*w) - so they no
! longer overlap. To search any arbitrary other position (not necessarily on the surface),
! we choose the most suitable of the 3*dim projections by looking at its maxloc(abs(xyz)),
! project and apply the same shift. This way we should find a triangle within one of the 3*dim
! projected sides of the cubed sphere surface mesh.
! Choosing GNOMONIC_BOX_WIDTH too small, we may end up with large surface triangles that do not
! map entirely within any of the projected sides. In 2D, the minimal angle of such a facet (line
! segment) would have to be 2*atan(w)-90 deg (for w=10 that's 78.5 deg)
real, parameter :: GNOMONIC_BOX_WIDTH = 10.0

contains
  
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
  
  to_nodes=nowned_nodes(to_fields(1))
  ! local coordinates is one more than horizontal coordinate dim
  allocate( seles(to_nodes), loc_coords(1:positions%dim, 1:to_nodes) )
  
  if (present(surface_name)) then
    lsurface_name=surface_name
  else
    lsurface_name="TempSurfaceName"
  end if
  
  ! project the positions of to_fields(1)%mesh into the horizontal plane
  ! and returns 'seles' (facet numbers) and loc_coords
  ! to tell where these projected nodes are found in the surface mesh
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
    
  if (IsParallel()) then
    do i=1, size(to_fields)
      call halo_update(to_fields(i))
    end do
  end if
  
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
  
  !! returned surface elements (facet numbers in 'mesh')
  !! and loc coords each node has been found in
  !! size(seles)==size(loc_coords,2)==nowned_nodes(mesh)
  integer, dimension(:), intent(out):: seles
  real, dimension(:,:), intent(out):: loc_coords
  
  type(vector_field):: mesh_positions
  type(vector_field), pointer:: horizontal_positions
  integer, dimension(:), pointer:: horizontal_mesh_list
  real, dimension(:,:), allocatable:: horizontal_coordinate
  real, dimension(vertical_normal%dim):: normal_vector
  integer:: i, stat, nodes
  
  assert(.not. mesh_periodic(positions))
  
  ! search only the first owned nodes
  nodes=nowned_nodes(mesh)
  
  assert( size(seles)==nodes )
  assert(size(loc_coords,1)==positions%dim)
  assert( size(loc_coords,2)==nodes )
  
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
  allocate( horizontal_coordinate(1:positions%dim-1, 1:nodes) )
  if (have_option('/geometry/spherical_earth')) then
    do i=1, nodes
      horizontal_coordinate(:,i) = map2horizontal_sphere(node_val(mesh_positions, i))
    end do
  else
    assert( vertical_normal%field_type==FIELD_TYPE_CONSTANT )
    normal_vector=node_val(vertical_normal,1)

    do i=1, nodes
      horizontal_coordinate(:,i) = map2horizontal(node_val(mesh_positions, i), normal_vector)
    end do
  end if
    
  call get_horizontal_positions(positions, surface_element_list, &
      vertical_normal, surface_name, &
      horizontal_positions, horizontal_mesh_list)
    
  call picker_inquire(horizontal_positions, horizontal_coordinate, &
    seles, loc_coords, global=.false. )
    
  ! in the spherical case some of the surface elements may be duplicated
  ! within horizontal positions, the returned seles should refer to entries
  ! in surface_element_list however - also check for nodes not found
  do i=1, size(seles)
    if (seles(i)<0) then
      ewrite(-1,*) "For node with coordinate", node_val(mesh_positions, i)
      ewrite(-1,*) "no top surface node was found."
      FLAbort("Something wrong with the geometry.")
    else
      seles(i) = horizontal_mesh_list(seles(i))
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
  
  ! Returns the horizontal positions on a horizontal mesh, this mesh
  ! may have some of the facets in surface_element_list duplicated -
  ! this is only in the spherical case where the horizontal mesh
  ! consists of multiple disjoint regions representing the sides of a cubed sphere
  ! and surface elements near the boundaries of these regions may be represented multiple times.
  ! Therefore we also return a map between elements in the horizontal positions mesh
  ! and facet numbers in the full mesh (in the planar case this is the same as surface_element_list)
  ! horizontal_positions is a borrowed reference, don't deallocate
  ! also horizontal_mesh_list should not be deallocated
  type(vector_field), pointer:: horizontal_positions
  integer, dimension(:), pointer:: horizontal_mesh_list
  
  if (.not. has_boundary_condition_name(positions, surface_name)) then
    call create_horizontal_positions(positions, &
        surface_element_list, vertical_normal, surface_name)
  end if
    
  horizontal_positions => extract_surface_field(positions, surface_name, &
    trim(surface_name)//"HorizontalCoordinate")

  call get_boundary_condition(positions, surface_name, &
      surface_element_list=horizontal_mesh_list)

end subroutine get_horizontal_positions
  
subroutine create_horizontal_positions(positions, surface_element_list, vertical_normal, surface_name)
! adds a "boundary condition" to 'positions' with an associated vector surface field containing a dim-1 
! horizontal coordinate field that can be used to map from the surface mesh specified 
! by 'surface_element_list'. This "boundary condition" will be stored under the name 'surface_name'
  type(vector_field), intent(inout):: positions
  type(vector_field), intent(in):: vertical_normal
  integer, dimension(:), intent(in):: surface_element_list
  character(len=*), intent(in):: surface_name
    
  type(mesh_type), pointer:: surface_mesh
  type(vector_field):: horizontal_positions, surface_positions
  type(vector_boundary_condition), pointer :: last_bc
  integer, dimension(:), pointer:: surface_node_list, surface_element_map
  real, dimension(vertical_normal%dim):: normal_vector
  integer:: i, node, nobcs
  
  call add_boundary_condition_surface_elements(positions, &
    name=surface_name, type="verticalextrapolation", &
    surface_element_list=surface_element_list)
  ! now get back the created surface mesh
  ! and surface_node_list a mapping between nodes in surface_mesh and nodes in positions%mesh
  call get_boundary_condition(positions, name=surface_name, &
    surface_mesh=surface_mesh, surface_node_list=surface_node_list)
    
  if (have_option('/geometry/spherical_earth')) then
    call allocate(surface_positions, dim=positions%dim, mesh=surface_mesh, &
      name="TempSurfacePositions")
    call remap_field_to_surface(positions, surface_positions, surface_element_list)
    call create_horizontal_positions_sphere(surface_positions, &
      horizontal_positions, surface_element_map, surface_name)

    call deallocate(surface_mesh)
    ! surface_mesh is a pointer, so this replaces the original surface_mesh with the new surface_mesh
    ! with duplicated facets for facets that occur in multiple coordinate domains
    surface_mesh = horizontal_positions%mesh
    call incref(surface_mesh)

    ! now we need to fix its surface_element_list
    ! we have surface_element_map which maps from elements in surface_positions to entries
    ! in the original surface_element_list, which in turn maps to facets in the positions%mesh
    ! the new surface_element_list is therefore simply a composition of these two maps

    nobcs = size(positions%bc%boundary_condition)
    last_bc => positions%bc%boundary_condition(nobcs)
    deallocate(last_bc%surface_element_list)
    allocate(last_bc%surface_element_list(1:size(surface_element_map)))
    last_bc%surface_element_list = surface_element_list(surface_element_map)
    deallocate(surface_element_map)

    ! similarly, reconstruct the surface_node_list (although I'm not sure we use this)
    deallocate(last_bc%surface_node_list)
    allocate(last_bc%surface_node_list(1:node_count(surface_positions)))
    do i=1, element_count(surface_positions)
      last_bc%surface_node_list(ele_nodes(surface_positions, i)) = face_global_nodes(positions, surface_element_list(i))
    end do

    call deallocate(surface_positions)
    
  else

    ! flat case

    assert(vertical_normal%field_type==FIELD_TYPE_CONSTANT)
    normal_vector=node_val(vertical_normal, 1)

    call allocate(horizontal_positions, dim=positions%dim-1, mesh=surface_mesh, &
      name=trim(surface_name)//"HorizontalCoordinate" )
    do i=1, size(surface_node_list)
      node=surface_node_list(i)
      call set(horizontal_positions, i, &
          map2horizontal(node_val(positions, node), normal_vector))
    end do
  end if

  call insert_surface_field(positions, name=surface_name, &
    surface_field=horizontal_positions)
    
  call deallocate(horizontal_positions)

end subroutine create_horizontal_positions

subroutine create_horizontal_positions_sphere(surface_positions, &
    horizontal_positions, element_map, mesh_name)
  type(vector_field), intent(in) :: surface_positions
  type(vector_field), intent(out) :: horizontal_positions
  integer, dimension(:), pointer :: element_map
  character(len=*), intent(in):: mesh_name

  type(mesh_type) :: horizontal_mesh
  real, dimension(1:surface_positions%dim, ele_loc(surface_positions,1)) :: xyz
  real, dimension(:,:), allocatable :: hxy
  integer, dimension(:), allocatable :: ele_start
  real shift
  integer ele, i, nele, dim, nloc

  dim = surface_positions%dim
  nloc = size(xyz, 2)

  allocate(hxy(dim-1, nloc*element_count(surface_positions)*3))
  allocate(ele_start(1:element_count(surface_positions)+1))

  nele = 0
  do ele = 1, element_count(surface_positions)
    xyz = ele_val(surface_positions, ele)
    ele_start(ele) = nele + 1
    do i = 1, dim
      if (all(abs(xyz(:,:)/spread(xyz(i,:),1,dim))<GNOMONIC_BOX_WIDTH)) then
        if (.not.(all(xyz(i,:)>0.) .or. all(xyz(i,:)<0.))) cycle
        shift = sign(3*i*GNOMONIC_BOX_WIDTH, xyz(i,1))
        nele = nele + 1
        call add_horizontal_ele(hxy(:, (nele-1)*nloc+1:nele*nloc), i, shift)
      end if
    end do
    if (nele<ele_start(ele)) then
      ewrite(-1,*) "This probably means the surface mesh is too coarse"
      FLAbort("Failed to construct horizontally projected surface mesh")
    end if
  end do
  ele_start(ele) = nele + 1

  call allocate(horizontal_mesh, nele*nloc, nele, surface_positions%mesh%shape, mesh_name)
  call allocate(horizontal_positions, dim-1, horizontal_mesh, trim(mesh_name)//"HorizontalCoordinate")

  horizontal_mesh%continuity = -1
  horizontal_mesh%ndglno = (/ (i, i=1, nele*nloc) /)
  horizontal_positions%val =  hxy(:, 1:nele*nloc)

  allocate(element_map(1:nele))
  do ele=1, element_count(surface_positions)
    element_map(ele_start(ele):ele_start(ele+1)-1) = ele
  end do

  call deallocate(horizontal_mesh)
  deallocate(hxy, ele_start)

  contains

  subroutine add_horizontal_ele(hxy, i, shift)
    real, dimension(:,:), intent(out) :: hxy
    integer, intent(in) :: i
    real, intent(in) :: shift

    integer k

    do k = 1, i-1
      hxy(k, :) = xyz(k, :)/xyz(i, :)
    end do
    do k = i+1, dim
      hxy(k-1, :) = xyz(k, :)/xyz(i, :)
    end do
    hxy(1,:) = hxy(1,:) + shift

  end subroutine add_horizontal_ele

end subroutine create_horizontal_positions_sphere

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

function map2horizontal_sphere(xyz) result (xy)
real, dimension(:), intent(in):: xyz
real, dimension(size(xyz)-1):: xy

  integer i

  i = maxloc(abs(xyz), dim=1)
  xy(1:i-1) = xyz(1:i-1)/xyz(i)
  xy(i:) = xyz(i+1:)/xyz(i)
  xy(1) = xy(1) + sign(3*i*GNOMONIC_BOX_WIDTH, xyz(i))
    
end function map2horizontal_sphere

function VerticalProlongationOperator(mesh, positions, vertical_normal, &
  surface_element_list, surface_mesh)
  !! creates a prolongation operator that prolongates values on
  !! a surface mesh to a full mesh below using the same interpolation
  !! as the vertical extrapolation code above. The transpose of this prolongation
  !! operator can be used as a restriction/clustering operator 
  !! from the full mesh to the surface.
  type(csr_matrix) :: VerticalProlongationOperator
  !! the mesh to which to prolongate, its nodes are only considered as
  !! a set of loose points
  type(mesh_type), target, intent(in):: mesh
  !! positions on the whole domain (doesn't have to be the same mesh)
  type(vector_field), intent(inout):: positions
  !! upward normal vector on the whole domain
  type(vector_field), target, intent(in):: vertical_normal
  !! the face numbers of the surface mesh
  integer, dimension(:), intent(in):: surface_element_list
  !! Optionally a surface mesh may be provided that represents the nodes
  !! /from/ which to interpolate (may be of different signature than 'mesh')
  !! Each column in the prolongator will correspond to a node in this surface_mesh.
  !! If not provided, the "from mesh" is considered to consist of faces
  !! given by surface_element_list (and same cont. and order as 'mesh').
  !! In this case however empty columns, i.e. surface nodes from which no
  !! value in the full mesh is interpolated, are removed and there is no
  !! necessary relation between column numbering and surface node numbering.
  type(mesh_type), intent(in), target, optional:: surface_mesh

  type(csr_sparsity):: sparsity
  type(mesh_type), pointer:: lsurface_mesh
  type(integer_hash_table):: facet2sele
  real, dimension(:,:), allocatable:: loc_coords
  real, dimension(:), allocatable:: mat
  real, dimension(:), allocatable:: coefs
  integer, dimension(:), pointer:: snodes
  integer, dimension(:), allocatable:: seles, colm, findrm, snod2used_snod
  integer i, j, k, rows, entries, count, snod, sele
  ! coefficient have to be at least this otherwise they're negligable in an interpolation
  real, parameter :: COEF_EPS=1d-10
  
  ! only assemble the rows associted with nodes we own
  rows=nowned_nodes(mesh)
  
  ! local coordinates is one more than horizontal coordinate dim
  allocate( seles(rows), loc_coords(1:positions%dim, 1:rows) )
  
  ! project the positions of to_fields(1)%mesh into the horizontal plane
  ! and returns 'seles' (indices in surface_element_list) and loc_coords 
  ! to tell where these projected nodes are found in the surface mesh
  call horizontal_picker(mesh, positions, vertical_normal, &
      surface_element_list, "TempSurfaceName", &
      seles, loc_coords)

  ! count upper estimate for n/o entries for sparsity
  entries=0
  do i=1, rows
     entries=entries+face_loc(mesh, seles(i))
  end do
  ! preliminary matrix:  
  allocate( mat(1:entries), findrm(1:rows+1), &
    colm(1:entries) )
  
  if (.not. present(surface_mesh)) then
     ! We use the entire surface mesh of 'mesh'
     lsurface_mesh => mesh%faces%surface_mesh
     ! Not all surface nodes may be used (i.e. interpolated from) - even
     ! within surface elements that /are/ in surface_element-list.
     ! We need a map between global surface node numbering and 
     ! a consecutive numbering of used surface nodes.
     ! (this will be the column numbering)
     allocate(snod2used_snod(1:node_count(lsurface_mesh)))
     snod2used_snod=0
     count=0 ! counts number of used surface nodes
  else
     lsurface_mesh => surface_mesh
     ! seles are facet element numbers in 'mesh'
     ! need to invert surface_element_list to convert these to element numbers in surface_mesh
     call invert_set(surface_element_list, facet2sele)
     do i=1, size(seles)
       seles(i) = fetch(facet2sele, seles(i))
     end do
  end if

  allocate(coefs(1:lsurface_mesh%shape%loc))
  entries=0 ! this time only count nonzero entries
  do i=1, rows

     ! beginning of each row in mat
     findrm(i)=entries+1

     sele = seles(i)
     snodes => ele_nodes(lsurface_mesh, sele)
     coefs = eval_shape(ele_shape(lsurface_mesh, sele), loc_coords(:,i))
     
     do j=1, size(snodes)
       snod=snodes(j)
       if (abs(coefs(j))>COEF_EPS) then
         if (.not. present(surface_mesh)) then
           if (snod2used_snod(snod)==0) then
             ! as of yet unused surface node
             count=count+1
             snod2used_snod(snod)=count
           end if
           ! this is the column index we're gonna use instead
           snod=snod2used_snod(snod)
         end if
         entries=entries+1
         colm(entries)=snod
         mat(entries)=coefs(j)
       end if
     end do
       
  end do
  findrm(i)=entries+1
  
  if (present(surface_mesh)) then
    ! we haven't counted used surface nodes, instead we're using all 
    ! nodes of surface mesh as columns
    count=node_count(surface_mesh)
  end if
  
  call allocate(sparsity, rows, count, &
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
  do i=1, rows
     do k=findrm(i), findrm(i+1)-1
       j=colm(k)
       call set(VerticalProlongationOperator, i, j, mat(k))
     end do
  end do
  
  if (.not. present(surface_mesh)) then
    deallocate( snod2used_snod )
  else
    call deallocate(facet2sele)
  end if
  deallocate( findrm, colm, mat, coefs )
  deallocate( seles, loc_coords )
  
  call remove_boundary_condition(positions, "TempSurfaceName")
  
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
    ewrite(-1,*) "Warning: vertical_element_ordering has detected a cycle."
    ewrite(-1,*) "(deck of cards problem). This may reduce the efficiency"
    ewrite(-1,*) "of your vertically sweeping solve."
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
     ! this is a row (column indices) in the mesh%faces%face_list matrix
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
      ewrite(-1,*) "If you select /geometry/ocean_boundaries, you also need to "//&
        &"set /physical_parameters/gravity"
      FLExit("Missing gravity!")
    end if
  end if
  
end subroutine vertical_extrapolation_module_check_options

end module vertical_extrapolation_module
