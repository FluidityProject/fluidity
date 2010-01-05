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
use sparse_tools
use state_module
use fields
use transform_elements
use boundary_conditions
use parallel_tools
use global_parameters, only: new_options, real_4, real_8
use spud
use dynamic_bin_sort_module
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
  
private

public CalculateTopBottomDistance
public VerticalExtrapolation, vertical_integration
public vertical_element_ordering
public VerticalProlongationOperator

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
  type(vector_field), pointer:: positions
  type(scalar_field), pointer:: distance
  
  integer, pointer, dimension(:):: surface_element_list

  ! the distance field to compute:
  distance => extract_scalar_field(state, name)
  ! its first boundary condition is on the related top or bottom mesh
  call get_boundary_condition(distance, 1, &
    surface_element_list=surface_element_list)
  positions => extract_vector_field(state, "Coordinate")
  
  ! in each node of the mesh, set "distance" to the vertical coordinate 
  ! of this node projected to the above/below surface mesh
  call VerticalExtrapolation(vertical_coordinate, distance, positions, &
    surface_element_list)
    
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
    ! at the moment only 3d and n_g=(0,0,-1) works
    assert(gravity_normal%val(1)%ptr(1)==0.0)
    assert(gravity_normal%val(2)%ptr(1)==0.0)
    assert(gravity_normal%val(3)%ptr(1)==-1.0)
    call inner_product(vertical_coordinate, gravity_normal, positions)
    ! gravity points down, we want a vertical coordinate that increases upward
    call scale(vertical_coordinate, -1.0)
  end if
  
end subroutine calculate_vertical_coordinate

subroutine VerticalExtrapolationScalar(from_field, to_field, &
  positions, surface_element_list)
  !!< This sub extrapolates the values on a horizontal 2D surface
  !!< in the vertical direction to a 3D field
  !! The from_field may be a 3D field of which only the values on the
  !! 2D horizontal surface are used, or it is a 2D field defined
  !! on the surface mesh.
  type(scalar_field), intent(in):: from_field
  !! Resulting extrapolated field. May be the same field or a field on
  !! a different mesh (different degree).
  type(scalar_field), intent(in):: to_field
  !! Needed for extrapolation (not necessarily same degree as to_field)
  type(vector_field), target, intent(inout):: positions
  !! the surface elements (faces numbers) that make up the surface
  integer, dimension(:), intent(in):: surface_element_list
  
  call VerticalExtrapolationMultiple( (/ from_field /) , (/ to_field /), &
      positions, surface_element_list)
  
end subroutine VerticalExtrapolationScalar

subroutine VerticalExtrapolationVector(from_field, to_field, &
  positions, surface_element_list)
  !!< This sub extrapolates the values on a horizontal 2D surface
  !!< in the vertical direction to a 3D field
  !! The from_field may be a 3D field of which only the values on the
  !! 2D horizontal surface are used, or it is a 2D field defined
  !! on the surface mesh.
  type(vector_field), intent(in):: from_field
  !! Resulting extrapolated field. May be the same field or a field on
  !! a different mesh (different degree).
  type(vector_field), intent(in):: to_field
  !! Needed for extrapolation (not necessarily same degree as to_field)
  type(vector_field), target, intent(inout):: positions
  !! the surface elements (faces numbers) that make up the surface
  integer, dimension(:), intent(in):: surface_element_list
  
  type(scalar_field), dimension(from_field%dim):: from_field_components, to_field_components
  integer i
  
  assert(from_field%dim==to_field%dim)
  
  do i=1, from_field%dim
     from_field_components(i)=extract_scalar_field(from_field, i)
     to_field_components(i)=extract_scalar_field(to_field, i)
  end do
    
  call VerticalExtrapolationMultiple( from_field_components, to_field_components, &
        positions, surface_element_list)
  
end subroutine VerticalExtrapolationVector

subroutine VerticalExtrapolationMultiple(from_fields, to_fields, &
  positions, surface_element_list)
  !!< This sub extrapolates the values on a horizontal 2D surface
  !!< in the vertical direction to 3D fields
  !! The from_fields may be 3D fields of which only the values on the
  !! 2D horizontal surface are used, or it is 2D fields defined
  !! on the surface mesh.
  !! This version takes multiple from_fields at the same time and extrapolates
  !! to to_fields, such that the surface search only has to be done once. This 
  !! will only work if all the from_fields are on the same mesh, and on all the 
  !! to_field are on the same (possibly a different) mesh.
  type(scalar_field), dimension(:), intent(in):: from_fields
  !! Resulting extrapolated field. May be the same field or a field on
  !! a different mesh (different degree).
  type(scalar_field), dimension(:), intent(in):: to_fields
  !! Needed for extrapolation (not necessarily same degree as to_field)
  type(vector_field), target, intent(inout):: positions
  integer::flat_earth_int

  !! the surface elements (faces numbers) that make up the surface
  integer, dimension(:), intent(in):: surface_element_list

  type(vector_field) to_positions
  type(scalar_field), dimension(size(to_fields)):: to_fields_non_periodic
  type(scalar_field) to_field_copy
  type(mesh_type), pointer:: x_mesh
  type(mesh_type):: to_mesh
  type(element_type), pointer:: from_shape, to_shape
  real, dimension(:), pointer:: x, y, z
  real, dimension(:), allocatable:: from_val, loc_coords
  real value
  integer, dimension(:), allocatable:: flv, senlist, sele_nodes
  integer, dimension(:), allocatable:: tri_ids
  real, dimension(:), allocatable:: shape_fxn
  integer xnonod, from_vertices, from_loc
  integer i, j, sele, nod

  if(have_option('/geometry/spherical_earth')) then
     flat_earth_int = 0
  else
     flat_earth_int = 1
  end if

  assert(size(from_fields)==size(to_fields))
  do i=2, size(to_fields)
     assert(to_fields(1)%mesh==to_fields(i)%mesh)
     assert(from_fields(1)%mesh==from_fields(i)%mesh)
  end do
  
  x_mesh => positions%mesh
  if (to_fields(1)%mesh%periodic) then
     ! we first extrapolate to a non-periodic field
     to_shape => ele_shape(to_fields(1), 1)
     if (to_shape%degree==x_mesh%shape%degree) then
        ! can put it on the x_mesh
        call allocate(to_fields_non_periodic(1), x_mesh, &
             & name='ToFieldNonPeriodic_VerticalExtrapolation')
     else
        ! need to make our own non-periodic mesh
        to_mesh=make_mesh(x_mesh, to_shape, continuity=continuity(to_fields(1)),&
             name="NonPeriodicMesh_VerticalExtrapolation") 
        call allocate(to_fields_non_periodic(1), to_mesh, &
             & name='ToFieldNonPeriodic_VerticalExtrapolation')
     end if
     ! we only need to allocate one field and then reuse it for the rest
     to_fields_non_periodic(2:)=to_fields_non_periodic(1)
  else
     ! directly put result in to_fields
     to_fields_non_periodic=to_fields
  end if
     
  
  if (to_fields_non_periodic(1)%mesh==x_mesh) then
    to_positions=positions
    ! make to_positions indep. ref. of the field, so we can deallocate it 
    ! safely without destroying positions
    call incref(to_positions)
  else
    call allocate(to_positions, positions%dim, to_fields_non_periodic(1)%mesh, &
      name='ToPositions_VerticalExtrapolation')
    call remap_field(positions, to_positions)
  end if
  x => to_positions%val(1)%ptr
  y => to_positions%val(2)%ptr
  z => to_positions%val(3)%ptr
  xnonod=size(x) ! note that this is the same number of nodes as in to_fields_non_periodic
  
  ! number of vertices in surface element, has to be 3, a triangle
  from_vertices=3
  assert(size(local_vertices(face_shape(to_positions,1)))==from_vertices)
  ! local node number of the vertices of a face (surface element)
  !ATTENTION doesn't work with meshes with varying element types
  allocate(flv(1:from_vertices))
  flv=local_vertices(face_shape(to_positions,1))
  
  allocate( senlist(size(surface_element_list) * from_vertices) )
  allocate( sele_nodes(1:face_loc(to_positions,1)) )
  
  ! make linear, triangular surface mesh, regardless of mesh degree of
  ! to_field/to_positions, but using its node numbering
  do i=1, size(surface_element_list)
    sele=surface_element_list(i)
    ! the global nodes of each surface element
    sele_nodes=face_global_nodes(to_positions, sele)
    ! copy the vertices thereof in senlist
    do j=1, from_vertices
      senlist((i-1)*from_vertices+j)=sele_nodes(flv(j))
    end do
  end do
  deallocate( sele_nodes, flv)
  
  ! Allocate some tempory space. This could be cached.
  allocate(tri_ids(xnonod))
  allocate(shape_fxn(xnonod*from_vertices))
  
  ! Perform spatial search
  call VerticalShellMapper_find(x, y, z, xnonod, senlist, &
       size(surface_element_list), 1, flat_earth_int, tri_ids, shape_fxn)
  deallocate( senlist )
  call deallocate(to_positions)
  
  if (element_count(from_fields(1))==size(surface_element_list)) then
     ! assume the from field is on the surface mesh only
     ! elements of this field are therefore the surface elements of to_fields
     from_loc=ele_loc(from_fields(1), 1)
     from_shape => ele_shape(from_fields(1), 1)
  else
     ! the from field is on the full mesh and we're only interested
     ! in its values at its faces that form the surface mesh
     from_loc=face_loc(from_fields(1),1)
     from_shape => face_shape(from_fields(1),1)
  end if
  allocate( from_val(1:from_loc), loc_coords(1:from_vertices) )
  
  do i=1, size(to_fields)
    do nod=1, xnonod
      if (element_count(from_fields(i))==size(surface_element_list)) then
         from_val=ele_val( from_fields(i), tri_ids(nod) )
      else
         ! surface element/triangle above this node
         sele=surface_element_list(tri_ids(nod))
         from_val=face_val( from_fields(i), sele )
      end if
      value=0.0
      ! local coordinates of the node projected on the surface element
      loc_coords=shape_fxn( (nod-1)*from_vertices+1 : nod*from_vertices )
      
      
      do j=1, from_loc
        value=value + eval_shape(from_shape, j, loc_coords)* &
           from_val(j)
      end do
        
      call set(to_fields_non_periodic(i), nod, value)
      
    end do
  end do
    
  if (to_fields(1)%mesh%periodic) then
     ewrite(1,*) 'Mapping from non-periodic to periodic in VerticalExtrapolation'
     do i=1, size(to_fields)
        ! this is a bit of a hack: we want to_fields to be intent(in)
        ! so we remap the field to a copy of to_fields(i) that points to the same value space
        to_field_copy=to_fields(i)
        call remap_field(to_fields_non_periodic(i), to_field_copy)
     end do
     ! we only allocated one, so only need to deallocate one
     call deallocate(to_fields_non_periodic(1))
     if(to_shape%degree/=1) then
        call deallocate(to_mesh)
     end if
  endif
  
  deallocate( tri_ids, shape_fxn, from_val, loc_coords )
  
end subroutine VerticalExtrapolationMultiple

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

end module vertical_extrapolation_module
