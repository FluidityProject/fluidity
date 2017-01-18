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
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
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

module spherical_adaptivity
  use spud
  use fldebug
  use global_parameters, only : OPTION_PATH_LEN
  use mpi_interfaces
  use futils, only: int2str
  use parallel_tools
  use fields
  use parallel_fields
  use state_module
  use vtk_interfaces
  use vertical_extrapolation_module
  use field_options
  use boundary_conditions
  use pickers
  use reserve_state_module
  use populate_state_module
  use halos
  implicit none

  private

  public prepare_spherical_adaptivity, spherical_adaptivity_pop_in, spherical_adaptivity_pop_out

  contains

  subroutine prepare_spherical_adaptivity(states, current_positions)
    type(state_type), dimension(:), intent(inout) :: states
    type(vector_field), intent(inout) :: current_positions ! FIXME: intent(in)

    type(vector_field), pointer :: external_positions
    type(vector_field) :: horizontal_positions, top_positions, bottom_positions
    type(mesh_type), pointer :: external_mesh
    integer, dimension(:), pointer :: element_map

    if (.not. has_vector_field(states(1), "BaseGeometryMeshHorizontalCoordinate")) then
      external_mesh => get_external_mesh(states)
      external_positions => get_external_coordinate_field(states(1), external_mesh)
      assert(external_positions%mesh==external_mesh)

      call create_horizontal_positions_sphere(external_positions, &
        horizontal_positions, element_map, "BaseGeometryMesh")
      top_positions = construct_base_geometry_surface_positions(horizontal_positions, &
        external_mesh, element_map, current_positions, 'top', 'BaseGeometryTopPositions')
      bottom_positions = construct_base_geometry_surface_positions(horizontal_positions, &
        external_mesh, element_map, current_positions, 'bottom', 'BaseGeometryBottomPositions')

      if (isparallel()) then
        call base_geometry_allgather(external_mesh, element_map, horizontal_positions, top_positions, bottom_positions)
      end if

      !call insert(reserve_state, horizontal_positions%mesh, horizontal_positions%mesh%name)
      call insert(reserve_state, horizontal_positions, horizontal_positions%name)
      call insert(reserve_state, top_positions, top_positions%name)
      call insert(reserve_state, bottom_positions, bottom_positions%name)

      call insert(states, horizontal_positions, horizontal_positions%name)
      call insert(states, top_positions, top_positions%name)
      call insert(states, bottom_positions, bottom_positions%name)

      call deallocate(horizontal_positions)
      call deallocate(top_positions)
      call deallocate(bottom_positions)
      ! we can throw away element_map, as we only care about looking up top/bottom_positions
      ! which is on the same mesh as horizontal_positions - so that we never have to map back to
      ! the element numbering of external_mesh
      deallocate(element_map)
    end if

  end subroutine prepare_spherical_adaptivity

  subroutine base_geometry_allgather(external_mesh, element_map, horizontal_positions, top_positions, bottom_positions)
    type(mesh_type), intent(in) :: external_mesh
    integer, dimension(:), intent(in) :: element_map
    type(vector_field), intent(inout) :: horizontal_positions, top_positions, bottom_positions

    type(mesh_type) :: universal_mesh
    integer, dimension(:), allocatable :: owned_per_proc, displacements
    integer owned_elements, universal_elements
    integer rank, nprocs, communicator, ierr
    integer i, j, dim

    assert(associated(external_mesh%element_halos))
    assert(external_mesh%element_halos(1)%ordering_scheme==HALO_ORDER_TRAILING_RECEIVES)
    do i=1, element_count(horizontal_positions)
      if (.not. element_owned(external_mesh, element_map(i))) exit
    end do
    owned_elements = max(i-1, 0)
    dim = top_positions%dim
    assert(all(horizontal_positions%mesh%ndglno(1:owned_elements*dim)==(/ (i, i=1, owned_elements*dim) /)))

    communicator = halo_communicator(external_mesh%element_halos(1))
    rank = getrank(communicator)
    nprocs = getnprocs(communicator)

    ! work out n/o element owned by each proc
    allocate(owned_per_proc(1:nprocs), displacements(1:nprocs))
    call mpi_allgather(owned_elements, 1, getpinteger(), owned_per_proc, 1, getpinteger(), communicator, ierr)
    assert(ierr == MPI_SUCCESS)
    universal_elements = sum(owned_per_proc)

    ! now change it to store n/o owned nodes per proc
    owned_per_proc = owned_per_proc * dim
    ! and compute starting index in buffer for each when node data is gathered (0-indexed)
    j = 0
    do i = 1, size(displacements)
      displacements(i) = j
      j = j + owned_per_proc(i)
    end do

    ! setup trivial universal DG mesh
    call allocate(universal_mesh, universal_elements*dim, universal_elements, &
      horizontal_positions%mesh%shape, horizontal_positions%mesh%name)
    universal_mesh%continuity = -1
    do i=1, element_count(universal_mesh)
      call set_ele_nodes(universal_mesh, i, (/ (j, j=(i-1)*dim+1, i*dim) /))
    end do

    call allgather_positions(horizontal_positions, universal_mesh, owned_elements*dim, owned_per_proc, displacements, communicator)
    call allgather_positions(top_positions, universal_mesh, owned_elements*dim, owned_per_proc, displacements, communicator)
    call allgather_positions(bottom_positions, universal_mesh, owned_elements*dim, owned_per_proc, displacements, communicator)

  end subroutine base_geometry_allgather

  subroutine allgather_positions(positions, universal_mesh, owned_nodes, owned_per_proc, displacements, communicator)
    type(vector_field), intent(inout) :: positions
    type(mesh_type), intent(in) :: universal_mesh
    integer, intent(in) :: owned_nodes
    integer, dimension(:), intent(in) :: owned_per_proc, displacements
    integer, intent(in) :: communicator

    type(vector_field) :: universal_positions
    integer gdim, ierr

    gdim = positions%dim
    call allocate(universal_positions, gdim, universal_mesh, positions%name)
    call mpi_allgatherv(positions%val, owned_nodes*gdim, getpreal(), &
      universal_positions%val, owned_per_proc*gdim, displacements*gdim, getpreal(), &
      communicator, ierr)
    assert(ierr == MPI_SUCCESS)
    call deallocate(positions)
    positions = universal_positions

  end subroutine allgather_positions

  function construct_base_geometry_surface_positions(external_horizontal_positions, external_mesh, element_map, positions, &
    surface_name, positions_name) result (base_geometry_surface_positions)
    type(vector_field), intent(inout) :: external_horizontal_positions
    type(mesh_type), intent(in) :: external_mesh
    integer, dimension(:), intent(in) :: element_map
    type(vector_field), intent(inout) :: positions ! FIXME: should be intent(in)
    character(len=*), intent(in) :: surface_name, positions_name
    type(vector_field) :: base_geometry_surface_positions

    type(vector_field) :: surface_positions, surface_horizontal_positions, external_surface_positions
    type(mesh_type), pointer :: surface_mesh
    character(len=OPTION_PATH_LEN) :: extrude_region_option_path, surface_id_option_path
    real, dimension(:,:), allocatable :: loc_coords
    real, dimension(:), allocatable :: max_loc_coords
    integer, dimension(:), pointer :: surface_element_list, nodes
    integer, dimension(:), allocatable :: eles, surface_ids
    integer i, j


    extrude_region_option_path = trim(positions%mesh%option_path) // '/from_mesh/extrude/regions'
    allocate(surface_ids(option_count(extrude_region_option_path)))
    do i=1, size(surface_ids)
      surface_id_option_path = trim(extrude_region_option_path) // '[' // int2str(i-1) // ']/' // trim(surface_name) // '_surface_id'
      call get_option(surface_id_option_path, surface_ids(i))
    end do

    ! FIXME: adding a temp. bc. only to create a surface mesh from boundary ids
    call add_boundary_condition(positions, name=surface_name, type='temporary', boundary_ids=surface_ids)
    deallocate(surface_ids)

    call get_boundary_condition(positions, name=surface_name, surface_mesh=surface_mesh, &
     surface_element_list=surface_element_list)

    call allocate(surface_positions, positions%dim, surface_mesh, "TempSurfacePositions")
    call remap_field_to_surface(positions, surface_positions, surface_element_list)

    call remove_boundary_condition(positions, name=surface_name)

    call allocate(surface_horizontal_positions, positions%dim-1, surface_mesh, "TempSurfaceHorizontalPositions")
    do i=1, node_count(surface_positions)
      call set(surface_horizontal_positions, i, map2horizontal_sphere(node_val(surface_positions, i)))
    end do

    allocate(eles(1:node_count(surface_positions)), loc_coords(1:positions%dim, 1:node_count(surface_positions)))
    call picker_inquire(external_horizontal_positions, surface_horizontal_positions, &
      eles, loc_coords, global=.false.)
    call deallocate(surface_horizontal_positions)

    ! first we calculate the top or bottom positions on the external mesh
    call allocate(external_surface_positions, positions%dim, external_mesh, "TempExternalSurfacePositions")
    allocate(max_loc_coords(1:node_count(external_mesh)))
    max_loc_coords = 0.
    do i=1, size(eles)
      nodes => ele_nodes(external_mesh, element_map(eles(i)))
      do j=1, size(nodes)
        if (loc_coords(j, i)>max_loc_coords(nodes(j))) then
         max_loc_coords(nodes(j)) = loc_coords(j, i)
         call set(external_surface_positions, nodes(j), node_val(surface_positions, i))
       end if
     end do
    end do
    call deallocate(surface_positions)

    do i=1, size(max_loc_coords)
     if (max_loc_coords(i)-1.<-1e-8) then
       ewrite(-1, *) "Node in base geometry "  //  trim(external_mesh%name) // &
         " not present in " // trim(surface_name) // "-surface of mesh " // trim(positions%mesh%name)
       FLExit("Base geometry mesh not consistent with spherically adapted mesh")
     end if
    end do
    deallocate(max_loc_coords, eles, loc_coords)

    ! then we map it onto the external_horizontal_positions%mesh
    call allocate(base_geometry_surface_positions, positions%dim, external_horizontal_positions%mesh, &
        positions_name)
    do i=1, element_count(base_geometry_surface_positions)
      call set(base_geometry_surface_positions, ele_nodes(base_geometry_surface_positions, i), &
        ele_val(external_surface_positions, element_map(i)))
    end do
    call deallocate(external_surface_positions)

  end function construct_base_geometry_surface_positions

  subroutine spherical_adaptivity_pop_in(states, positions)
    type(state_type), dimension(:), intent(in) :: states
    type(vector_field), intent(inout) :: positions

    type(vector_field), pointer :: horizontal_base_geometry, top_positions, bottom_positions
    type(vector_field) :: horizontal_positions
    real, dimension(:, :), allocatable :: loc_coords
    real, dimension(positions%dim) :: xyz_top, xyz_bottom
    real r, r_top, r_bottom
    integer, dimension(:), allocatable :: eles
    integer i

    horizontal_base_geometry => extract_vector_field(reserve_state, "BaseGeometryMeshHorizontalCoordinate")
    call allocate(horizontal_positions, positions%dim-1, positions%mesh, "HorizontalPositions")
    do i = 1, node_count(positions)
      call set(horizontal_positions, i, map2horizontal_sphere(node_val(positions, i)))
    end do
    allocate(eles(1:node_count(positions)), loc_coords(1:positions%dim, 1:node_count(positions)))
    call picker_inquire(horizontal_base_geometry, horizontal_positions, eles, loc_coords, global=.false.)

    top_positions => extract_vector_field(reserve_state, "BaseGeometryTopPositions")
    bottom_positions => extract_vector_field(reserve_state, "BaseGeometryBottomPositions")
    do i = 1, node_count(positions)
      call ray_plane_interpolate_xyz_radius(ele_val(top_positions, eles(i)), node_val(positions, i), xyz_top, r_top)
      call ray_plane_interpolate_xyz_radius(ele_val(bottom_positions, eles(i)), node_val(positions, i), xyz_bottom, r_bottom)
      r = sqrt(sum(node_val(positions, i)**2))
      call set(positions, i, ((r-r_bottom) * xyz_top + (r_top-r) * xyz_bottom)/(r_top-r_bottom))
    end do

    deallocate(eles, loc_coords)

  end subroutine spherical_adaptivity_pop_in

  subroutine spherical_adaptivity_pop_out(states, positions)
    type(state_type), dimension(:), intent(in) :: states
    type(vector_field), intent(inout) :: positions

    type(vector_field), pointer :: horizontal_base_geometry, top_positions, bottom_positions
    type(vector_field) :: horizontal_positions
    real, dimension(:, :), allocatable :: loc_coords
    real, dimension(positions%dim) :: xyz_top, xyz_bottom
    real r_top, r_bottom, r_current, r_new, xyz_j
    integer, dimension(:), allocatable :: eles
    integer i, j

    horizontal_base_geometry => extract_vector_field(states, "BaseGeometryMeshHorizontalCoordinate")
    call allocate(horizontal_positions, positions%dim-1, positions%mesh, "HorizontalPositions")
    do i = 1, node_count(positions)
      call set(horizontal_positions, i, map2horizontal_sphere(node_val(positions, i)))
    end do
    allocate(eles(1:node_count(positions)), loc_coords(1:positions%dim, 1:node_count(positions)))
    call picker_inquire(horizontal_base_geometry, horizontal_positions, eles, loc_coords, global=.false.)

    top_positions => extract_vector_field(states, "BaseGeometryTopPositions")
    bottom_positions => extract_vector_field(states, "BaseGeometryBottomPositions")
    do i = 1, node_count(positions)
      call ray_plane_interpolate_xyz_radius(ele_val(top_positions, eles(i)), node_val(positions, i), xyz_top, r_top)
      call ray_plane_interpolate_xyz_radius(ele_val(bottom_positions, eles(i)), node_val(positions, i), xyz_bottom, r_bottom)
      j = maxloc(abs(xyz_top), dim=1)
      xyz_j = node_val(positions, j, i)
      r_current = sqrt(sum(node_val(positions, i)**2))
      r_new = ((xyz_j-xyz_bottom(j)) * r_top + (xyz_top(j)-xyz_j) * r_bottom)/(xyz_top(j)-xyz_bottom(j))
      call set(positions, i, r_new * node_val(positions, i) / r_current)
    end do

    deallocate(eles, loc_coords)

  end subroutine spherical_adaptivity_pop_out

  subroutine ray_plane_interpolate_xyz_radius(xyz_plane, xyz_ray, xyz, r)
    real, dimension(:,:), intent(in) :: xyz_plane
    real, dimension(:), intent(in) :: xyz_ray
    real, dimension(:), intent(out) :: xyz
    real, intent(out) :: r

    real, dimension(size(xyz)) :: tuv

    call ray_plane_intersection(xyz_plane, xyz_ray, tuv)
    ! tuv(1) is coefficient on ray
    xyz = tuv(1) * xyz_ray
    ! overwrite tuv(1) with additional plane coordinate for easy interpolation
    tuv(1) = 1.0 - sum(tuv(2:))
    r = dot_product(tuv, sqrt(sum(xyz_plane**2, dim=1)))

  end subroutine ray_plane_interpolate_xyz_radius

  subroutine ray_plane_intersection(xyz_plane, xyz_ray, tuv)
    real, dimension(:,:), intent(in) :: xyz_plane
    real, dimension(:), intent(in) :: xyz_ray
    real, dimension(:), intent(out) :: tuv

    real, dimension(size(tuv), size(tuv)) :: m
    integer i

    m(:,1) = xyz_ray
    do i=2, size(m,1)
      m(:,i) = xyz_plane(:,1) - xyz_plane(:,i)
    end do
    call invert(m)
    tuv = matmul(m, xyz_plane(:,1))

  end subroutine ray_plane_intersection

end module spherical_adaptivity
