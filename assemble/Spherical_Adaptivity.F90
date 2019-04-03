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
  use iso_c_binding
  use spud
  use fldebug
  use global_parameters, only : OPTION_PATH_LEN, FIELD_NAME_LEN, no_active_processes
  use vector_tools
  use mpi_interfaces
  use futils, only: int2str
  use parallel_tools
  use quadrature
  use elements
  use parallel_fields
  use fields
  use state_module
  use mesh_files
  use vtk_interfaces
  use boundary_conditions
  use pickers
  use vertical_extrapolation_module
  use field_options
  use reserve_state_module
  use halos
  use populate_state_module
  use adapt_integration, only: pain_functional
  implicit none

  private

  public prepare_spherical_adaptivity, spherical_adaptivity_pop_in, spherical_adaptivity_pop_out

  public check_inverted_elements, transform_metric

  integer, private :: MPI_TYPE_LOC_PHYS_COORDS, MPI_TYPE_LOC_PHYS_COORDS_LEN

  contains

  subroutine prepare_spherical_adaptivity(states, current_positions)
    type(state_type), dimension(:), intent(inout) :: states
    type(vector_field), intent(inout) :: current_positions ! FIXME: intent(in)

    type(vector_field) :: base_geometry_horizontal_positions, top_positions, bottom_positions, external_base_geometry
    integer, dimension(:), pointer :: element_map

    ewrite(1,*) "Inside prepare_spherical_adaptivity"

    if (.not. has_vector_field(states(1), "BaseGeometryMeshHorizontalCoordinate")) then
      external_base_geometry = read_base_geometry()

      call create_horizontal_positions_sphere(external_base_geometry, &
          base_geometry_horizontal_positions, element_map, "BaseGeometryMesh")
      top_positions = construct_base_geometry_surface_positions(base_geometry_horizontal_positions, &
          external_base_geometry%mesh, element_map, current_positions, &
          "/geometry/spherical_earth/base_geometry/top_surface_ids", 'BaseGeometryTopPositions')
      bottom_positions = construct_base_geometry_surface_positions(base_geometry_horizontal_positions, &
          external_base_geometry%mesh, element_map, current_positions, &
          "/geometry/spherical_earth/base_geometry/bottom_surface_ids", 'BaseGeometryBottomPositions')

      call insert(reserve_state, base_geometry_horizontal_positions, base_geometry_horizontal_positions%name)
      call insert(reserve_state, top_positions, top_positions%name)
      call insert(reserve_state, bottom_positions, bottom_positions%name)
      ! also need to insert the mesh separately, as this avoid printing of left-over ("tagged") references
      ! during adaptivity
      call insert(reserve_state, base_geometry_horizontal_positions%mesh, base_geometry_horizontal_positions%mesh%name)

      call insert(states, base_geometry_horizontal_positions, base_geometry_horizontal_positions%name)
      call insert(states, top_positions, top_positions%name)
      call insert(states, bottom_positions, bottom_positions%name)

      call deallocate(base_geometry_horizontal_positions)
      call deallocate(top_positions)
      call deallocate(bottom_positions)
      ! we can throw away element_map and external_base_geometry, as we only care about looking up top/bottom_positions
      ! which is on the same mesh as base_geometry_horizontal_positions
      deallocate(element_map)
      call deallocate(external_base_geometry)
    end if

  end subroutine prepare_spherical_adaptivity

  function read_base_geometry() result (position)
    type(vector_field) :: position

    character(len=OPTION_PATH_LEN) :: option_path, mesh_file_name
    character(len=FIELD_NAME_LEN) :: mesh_file_format
    integer :: quad_degree
    type(mesh_type) :: mesh
    type(element_type) :: shape
    type(quadrature_type) :: quad
    integer, dimension(:, :), allocatable :: ndglno
    integer :: dim, rank, ele, nelements, nnodes, stat
    integer :: quad_family, save_no_active_processes

    option_path = '/geometry/spherical_earth/base_geometry'
    if (.not. have_option(option_path)) then
      FLExit("Missing option /geometry/spherical_earth/base_geometry needed for spherical adaptivity.")
    end if
    call get_option(trim(option_path) // "/file_name", mesh_file_name)
    call get_option(trim(option_path) // "/format[0]/name", mesh_file_format)

    call get_option("/geometry/dimension", dim)
    call get_option("/geometry/quadrature/degree", quad_degree)
    quad_family = get_quad_family()

    rank = getrank()
    save_no_active_processes = no_active_processes
    no_active_processes = 1
    if (rank==0) then
      position = read_mesh_files(trim(mesh_file_name), &
          quad_degree=quad_degree, &
          quad_family=quad_family, &
          format=mesh_file_format)
      position%name = "ExternalBaseGeometryMeshCoordinate"
      position%mesh%name = "ExternalBaseGeometryMesh"

      if (.not. (position%dim == dim .or. mesh_dim(position) == dim-1)) then
        FLExit("Wrong dimension for base_geometry mesh.")
      end if
      if (associated(position%mesh%columns) .or. position%mesh%faces%has_discontinuous_internal_boundaries) then
        FLExit("Base geometry mesh is not supposed to have column ids or internal boundaries")
      end if

      nelements = element_count(position)
      nnodes = node_count(position)
      call MPI_bcast(nelements, 1, getpinteger(), 0, MPI_COMM_FEMTOOLS, stat)
      call MPI_bcast(position%mesh%ndglno, dim*nelements, getpinteger(), 0, MPI_COMM_FEMTOOLS, stat)
      call MPI_bcast(position%val, dim*nnodes, getpreal(), 0, MPI_COMM_FEMTOOLS, stat)

    else

      quad = make_quadrature(dim, dim-1, degree=quad_degree, family=quad_family)
      shape = make_element_shape(dim, dim-1, 1, quad)

      call MPI_bcast(nelements, 1, getpinteger(), 0, MPI_COMM_FEMTOOLS, stat)
      allocate(ndglno(1:dim, nelements))
      call MPI_bcast(ndglno, dim*nelements, getpinteger(), 0, MPI_COMM_FEMTOOLS, stat)
      nnodes = maxval(ndglno)
      call allocate(mesh, nodes=nnodes, elements=nelements, shape=shape, name="ExternalBaseGeometryMesh")
      call allocate(position, dim, mesh,  "ExternalBaseGeometryMeshCoordinate")
      do ele=1, nelements
        call set_ele_nodes(mesh, ele, ndglno(:, ele))
      end do
      call add_faces(position%mesh)

      call MPI_bcast(position%val, dim*nnodes, getpreal(), 0, MPI_COMM_FEMTOOLS, stat)

      ! Reference counting cleanups.
      call deallocate(mesh)
      call deallocate(quad)
      call deallocate(shape)

    end if
    no_active_processes = save_no_active_processes

  end function read_base_geometry

  subroutine base_geometry_allgather(external_mesh, element_map, base_geometry_horizontal_positions, top_positions, bottom_positions)
    type(mesh_type), intent(in) :: external_mesh
    integer, dimension(:), intent(in) :: element_map
    type(vector_field), intent(inout) :: base_geometry_horizontal_positions, top_positions, bottom_positions

    type(mesh_type) :: universal_mesh
    integer, dimension(:), allocatable :: owned_per_proc, displacements
    integer owned_elements, universal_elements
    integer rank, nprocs, communicator, ierr
    integer i, j, dim

    assert(associated(external_mesh%element_halos))
    assert(external_mesh%element_halos(1)%ordering_scheme==HALO_ORDER_TRAILING_RECEIVES)
    do i=1, element_count(base_geometry_horizontal_positions)
      if (.not. element_owned(external_mesh, element_map(i))) exit
    end do
    owned_elements = max(i-1, 0)
    dim = top_positions%dim
    assert(all(base_geometry_horizontal_positions%mesh%ndglno(1:owned_elements*dim)==(/ (i, i=1, owned_elements*dim) /)))

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
      base_geometry_horizontal_positions%mesh%shape, base_geometry_horizontal_positions%mesh%name)
    universal_mesh%continuity = -1
    do i=1, element_count(universal_mesh)
      call set_ele_nodes(universal_mesh, i, (/ (j, j=(i-1)*dim+1, i*dim) /))
    end do

    call allgather_positions(base_geometry_horizontal_positions, universal_mesh, owned_elements*dim, owned_per_proc, displacements, communicator)
    call allgather_positions(top_positions, universal_mesh, owned_elements*dim, owned_per_proc, displacements, communicator)
    call allgather_positions(bottom_positions, universal_mesh, owned_elements*dim, owned_per_proc, displacements, communicator)

    call deallocate(universal_mesh)
    deallocate(owned_per_proc, displacements)

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

  subroutine max_first_argument(invec, inoutvec, len, type)
    ! subroutine used as mpi reduction operator which
    ! reduces an array of type MPI_TYPE_LOC_PHYS_COORDS
    ! where each entry exists of MPI_TYPE_LOC_PHYS_COORDS_LEN reals (nreals below)
    ! the reduction takes the maximum value of the first real, and the other nreals-1
    ! reals are taken from the same entry whose first real has the maximum value
    real, dimension(*) :: invec, inoutvec
    integer :: len, type

    integer :: size, ierr, i
    integer :: nreals

    assert(type==MPI_TYPE_LOC_PHYS_COORDS)
    nreals = MPI_TYPE_LOC_PHYS_COORDS_LEN
    call mpi_type_size(type, size, ierr)
    assert(c_sizeof(1.0)*nreals==size)

    do i=1, len
      if (invec((i-1)*nreals+1)>inoutvec((i-1)*nreals+1)) then
        inoutvec((i-1)*nreals+1:i*nreals) = invec((i-1)*nreals+1:i*nreals)
      end if
    end do

  end subroutine max_first_argument

  function construct_base_geometry_surface_positions(base_geometry_horizontal_positions, external_base_mesh, element_map, positions, &
    surface_ids_option_path, positions_name) result (base_geometry_surface_positions)
    type(vector_field), intent(inout) :: base_geometry_horizontal_positions
    type(mesh_type), intent(in) :: external_base_mesh
    integer, dimension(:), intent(in) :: element_map
    type(vector_field), intent(inout) :: positions ! FIXME: should be intent(in)
    character(len=*), intent(in) :: surface_ids_option_path, positions_name
    type(vector_field) :: base_geometry_surface_positions

    type(vector_field) :: surface_positions, surface_horizontal_positions, external_surface_positions
    type(mesh_type), pointer :: surface_mesh
    real, dimension(:,:), allocatable :: loc_coords
    real, dimension(:,:), allocatable :: max_loc_phys_coords, max_loc_phys_coords_global
    integer, dimension(:), pointer :: surface_element_list, nodes
    integer, dimension(:), allocatable :: eles, surface_ids
    integer shape_option(2)
    integer mpi_max_loc_phys_coords_op
    integer i, j, ierr

    shape_option = option_shape(surface_ids_option_path)
    allocate(surface_ids(1:shape_option(1)))
    call get_option(surface_ids_option_path, surface_ids)

    ! Step 1: create surface_positions and surface_horizontal_positions on the (local) surface_mesh
    ! defined by surface_ids

    ! FIXME: adding a temp. bc. only to create a surface mesh from boundary ids
    call add_boundary_condition(positions, name="TemporarySurface", type='temporary', boundary_ids=surface_ids)
    deallocate(surface_ids)

    call get_boundary_condition(positions, name="TemporarySurface", surface_mesh=surface_mesh, &
     surface_element_list=surface_element_list)

    call allocate(surface_positions, positions%dim, surface_mesh, "TempSurfacePositions")
    call remap_field_to_surface(positions, surface_positions, surface_element_list)

    call allocate(surface_horizontal_positions, positions%dim-1, surface_mesh, "TempSurfaceHorizontalPositions")
    do i=1, node_count(surface_positions)
      call set(surface_horizontal_positions, i, map2horizontal_sphere(node_val(surface_positions, i)))
    end do

    call remove_boundary_condition(positions, name="TemporarySurface")

    ! Step 2: search for nodes of the local surface_mesh in the base_geometry mesh (using horizontal coordinates)
    ! The physical (3D) coordinates of those surface mesh nodes that have the same horizontal coordinate
    ! as a node in the base geometry mesh will be used to set the physical (3D) coordinates of the base geometry node
    allocate(eles(1:node_count(surface_positions)), loc_coords(1:positions%dim, 1:node_count(surface_positions)))
    call picker_inquire(base_geometry_horizontal_positions, surface_horizontal_positions, &
      eles, loc_coords, global=.false.)
    call deallocate(surface_horizontal_positions)

    ! Step 3: for each node of the external_base_mesh, copy the physical coordinates (in surface_positions)
    ! of the node in the surface_mesh that has the maximum local coordinate wrt to the external_base_mesh node
    !
    ! the maximum local coordinate of each node in surface_mesh with respect to nodes in external_base_mesh
    ! a surface_mesh node that coincides (in horizontal coordinates) with a node in external_base_mesh
    ! will have a local coordinate of 1.0

    ! max_loc_coords(1,:) will contain this maximum local coordinate, and max_loc_coords(2:,:)
    ! the corresponding physical coordinate from surface_positions
    allocate(max_loc_phys_coords(1+positions%dim, 1:node_count(external_base_mesh)))
    max_loc_phys_coords = 0.
    do i=1, size(eles)
      nodes => ele_nodes(external_base_mesh, element_map(eles(i)))
      do j=1, size(nodes)
        if (loc_coords(j, i)>max_loc_phys_coords(1, nodes(j))) then
         max_loc_phys_coords(1, nodes(j)) = loc_coords(j, i)
         max_loc_phys_coords(2:, nodes(j)) = node_val(surface_positions, i)
       end if
     end do
    end do
    call deallocate(surface_positions)

    ! Step 4: now take the maximum of max_loc_coords(1, :) over all processes, and copy the corresponding physical coordinate

    ! define a mpi type storing one entry max_loc_coords(:,i)
    call mpi_type_contiguous(1+positions%dim, getpreal(), MPI_TYPE_LOC_PHYS_COORDS, ierr)
    call mpi_type_commit(MPI_TYPE_LOC_PHYS_COORDS, ierr)
    MPI_TYPE_LOC_PHYS_COORDS_LEN = 1+positions%dim
    ! define the reduction operator that takes the maximum local coordinates and copies the physical coordinates
    call mpi_op_create(max_first_argument, 1, mpi_max_loc_phys_coords_op, ierr)
    ! applies this reduction over all processes and store the result in max_loc_phys_coords_global
    allocate(max_loc_phys_coords_global(1+positions%dim, 1:node_count(external_base_mesh)))
    call mpi_allreduce(max_loc_phys_coords, max_loc_phys_coords_global, size(max_loc_phys_coords, 2), &
      MPI_TYPE_LOC_PHYS_COORDS, mpi_max_loc_phys_coords_op, MPI_COMM_FEMTOOLS, ierr)
    call mpi_type_free(MPI_TYPE_LOC_PHYS_COORDS, ierr)
    call mpi_op_free(mpi_max_loc_phys_coords_op, ierr)

    ! Step 5: if any max_local_coords_global(1,:)/=1 it means the corresponding node in the external_base_mesh does not
    ! have a node with the same horizontal coordinates in the surface_mesh
    do i=1, size(max_loc_phys_coords_global, 2)
     if (max_loc_phys_coords_global(1, i)-1.<-1e-8) then
       ewrite(-1, *) "Constructing " // positions_name
       ewrite(-1, *) "Node in base geometry "  //  trim(external_base_mesh%name) // &
         " not present in top or bottom surface of mesh " // trim(positions%mesh%name)
       FLExit("Base geometry mesh not consistent with spherically adapted mesh")
     end if
    end do

    ! Step 6: copy the result to a vector field on the external_base_mesh
    call allocate(external_surface_positions, positions%dim, external_base_mesh, "TempExternalSurfacePositions")
    call set_all(external_surface_positions, max_loc_phys_coords_global(2:,:))
    deallocate(max_loc_phys_coords, max_loc_phys_coords_global, eles, loc_coords)

    ! Step 7: then we map it onto the base_geometry_horizontal_positions%mesh (which contains discontinuties)
    call allocate(base_geometry_surface_positions, positions%dim, base_geometry_horizontal_positions%mesh, &
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

    ewrite(1,*) "Inside spherical_adaptivity_pop_in"

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

    call deallocate(horizontal_positions)
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

    ewrite(1,*) "Inside spherical_adaptivity_pop_out"

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
      j = maxloc(abs(xyz_top), dim=1)
      xyz_j = node_val(positions, j, i)
      r_current = sqrt(sum(node_val(positions, i)**2))
      r_new = ((xyz_j-xyz_bottom(j)) * r_top + (xyz_top(j)-xyz_j) * r_bottom)/(xyz_top(j)-xyz_bottom(j))
      call set(positions, i, r_new * node_val(positions, i) / r_current)
    end do

    call deallocate(horizontal_positions)
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

  subroutine ray_plane_intersection_gradient(xyz_plane, xyz_ray, t, grad_tuv)
    real, dimension(:,:), intent(in) :: xyz_plane
    real, dimension(:), intent(in) :: xyz_ray
    real, intent(in) :: t
    real, dimension(:, :), intent(out) :: grad_tuv

    real, dimension(size(xyz_ray), size(xyz_ray)) :: m
    real, dimension(size(xyz_ray)) :: rhs
    integer i, k

    ! in the previous subroutine we have solved:
    !     x1_i = t x_i + u (x1-x2)_i + v (x1-x3)_i
    ! taking the gradient:
    !      0_ki = grad_k(t) x_i + t I_ki + grad_k(u) (x1-x2)_i + grad_k(v) (x1-x3)_i
    ! which gives a linear s/o eqns for each k to solve grad_k(t,u,v)
    ! note that for each k, the lhs is the same, only the rhs, -t I_ki, changes

    m(:,1) = xyz_ray
    do i=2, size(m,1)
      m(:,i) = xyz_plane(:,1) - xyz_plane(:,i)
    end do
    call invert(m)

    do k=1, size(grad_tuv,1)
      rhs = 0.
      rhs(k) = -t
      grad_tuv(k,:) = matmul(m, rhs)
    end do

  end subroutine ray_plane_intersection_gradient

  subroutine check_inverted_elements(positions, base_geometry, metric)
    ! routine that checks whether elements have been inverted during
    ! popping out to the spherical geometry after and adapt
    ! It will print loads of debugging information about any inverted
    ! element it finds before producing an error message (if any are found).
    type(vector_field), intent(in) :: positions, base_geometry
    type(tensor_field), intent(in) :: metric

    real, dimension(positions%dim, positions%dim+1) :: tet_xyz
    real vol1, vol2
    integer ele, j, ele2, k, dim
    integer, dimension(:), pointer:: neigh, facets, nodes, nodes2
    logical :: any_inverted

    ewrite(1, *) "Inside check_inverted_elements"

    dim = positions%dim
    if (dim<2 .or. dim>3 .or. ele_loc(positions, 1)/=dim+1) then
      FLAbort("This routine should only be used for linear meshes in 2d or 3d")
    endif

    any_inverted = .false.
    do ele=1, element_count(positions)
      neigh => ele_neigh(positions, ele)
      facets => ele_faces(positions, ele)
      nodes => ele_nodes(positions, ele)
      do j=1, size(neigh)
        if (neigh(j)>0) then
          tet_xyz(:, 1:dim) = face_val(positions, facets(j))
          tet_xyz(:, dim+1) = node_val(positions, nodes(j))
          if (dim==3) then
            vol1 = tetvol(tet_xyz(1,:), tet_xyz(2,:), tet_xyz(3,:))
          else
            vol1 = cross_product2(tet_xyz(:,2)-tet_xyz(:,1), tet_xyz(:,3)-tet_xyz(:,1))
          end if

          ele2 = neigh(j)
          nodes2 => ele_nodes(positions, ele2)
          k = local_face_number(positions, face_neigh(positions, facets(j)))
          tet_xyz(:, positions%dim+1) = node_val(positions, nodes2(k))
          if (dim==3) then
            vol1 = tetvol(tet_xyz(1,:), tet_xyz(2,:), tet_xyz(3,:))
          else
            vol1 = cross_product2(tet_xyz(:,2)-tet_xyz(:,1), tet_xyz(:,3)-tet_xyz(:,1))
          end if

          if (vol1*vol2>0.0) then
            ewrite(2,*) "INVERTED ELEMENT!!!"
            ewrite(2,*) "Element 1 ***"
            call ele_info(ele)
            ewrite(2,*) "Element 2 ***"
            call ele_info(ele2)
            any_inverted = .true.
          end if
        end if
      end do
    end do

    if (any_inverted) then
      ewrite(0,*) "Inverted elements were found after popping out to the new adapted spherical geometry."
      ewrite(0,*) "This typically means that either the base geometry (the original input mesh) is too coarse,"
      ewrite(0,*) "or the aspect ratio bound needs to be lowered."
      FLExit("Inverted elements in spherical adaptivity.")
    end if

    contains

    subroutine ele_info(ele)
      integer, intent(in) :: ele
      integer, dimension(:), pointer :: nodes
      integer :: i
      real :: f

      ewrite(2, *) "Element number: ", ele
      ewrite(2, *) "Popped out positions:", ele_val(positions, ele)
      ewrite(2, *) "Base geometry positions:", ele_val(base_geometry, ele)
      f = pain_functional(ele, base_geometry, metric, verbose=.true.)
      ewrite(2,*) "Pain functional:", f
      ewrite(2,*) "Metric at nodes:"
      nodes => ele_nodes(metric, ele)
      do i=1, size(nodes)
        ewrite(2,*) node_val(metric, nodes(i))
      end do

    end subroutine ele_info

  end subroutine check_inverted_elements

  subroutine transform_metric(states, positions, metric)
    type(state_type), dimension(:), intent(in) :: states
    type(vector_field), intent(inout) :: positions
    type(tensor_field), intent(inout) :: metric

    type(vector_field), pointer :: horizontal_base_geometry, top_positions, bottom_positions
    type(vector_field) :: horizontal_positions
    real, dimension(:, :), allocatable :: loc_coords
    real, dimension(positions%dim) :: xyz_top, xyz_bottom
    real, dimension(positions%dim, positions%dim) :: Jac
    real r_top, r_bottom, r_current, r_new, xyz_j
    integer, dimension(:), allocatable :: eles
    integer i, j

    assert(positions%mesh == metric%mesh)
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
      j = maxloc(abs(xyz_top), dim=1)
      xyz_j = node_val(positions, j, i)
      r_current = sqrt(sum(node_val(positions, i)**2))
      r_new = ((xyz_j-xyz_bottom(j)) * r_top + (xyz_top(j)-xyz_j) * r_bottom)/(xyz_top(j)-xyz_bottom(j))
      Jac = pop_out_jacobian(node_val(positions, i), r_current, r_new, &
              ele_val(top_positions, eles(i)), ele_val(bottom_positions, eles(i)), &
              r_top, r_bottom)
      call set(metric, i, matmul(transpose(Jac), matmul(node_val(metric, i), Jac)))
    end do

    call deallocate(horizontal_positions)
    deallocate(eles, loc_coords)

  end subroutine transform_metric

  function pop_out_jacobian(xyz, r_current, r_new, xyz_top, xyz_bottom, r_top, r_bottom) result (J)
    real, dimension(:), intent(in) :: xyz
    real, intent(in) :: r_current, r_new, r_top, r_bottom
    real, dimension(:,:), intent(in) :: xyz_top, xyz_bottom
    real, dimension(size(xyz), size(xyz)) :: J

    real, dimension(size(xyz), size(xyz)) :: grad_tuv_t, grad_tuv_b
    real, dimension(size(xyz)) :: grad_r_new, r_t_nodes, r_b_nodes
    real t_t, t_b
    integer i, k

    ! taking the gradient of phi(x) = (x/r_current) * r_new

    ! first grad(x/r_current) * r_new = ( I/r_current - x\otimes x/r_current**3 ) * r_new
    do i=1, size(xyz)
      do k=1, size(xyz)
        J(i,k) = -xyz(i)*xyz(k)/r_current**3
      end do
      J(i,i) = J(i,i) + 1./r_current
    end do

    J = J*r_new

    ! then add (x/r_current) \otimes grad(r_new)

    ! we have (using a linear interpolation along the ray t*x for t_b<t<t_r at t=1)
    !     r_new = [(t_t-1)*r_b + (1-t_b)*r_t]/(t_t-t_b)
    ! thus:
    !     grad r_new = [(grad t_t)*r_b-(grad t_b)*r_t]/(t_t-t_b) - r_new/(t_t-t_b) * grad(t_t-t_b)
    !                  + [(t_t-1)*grad r_b + (1-t_b)*grad r_t]/(t_t-t_b)
    !                = [(r_b-r_new)*grad t_t + (t_t-1)*grad r_b + (r_new-r_t)*grad t_b + (1-t_b)*grad r_t]/(t_t-t_b)
    ! with:
    !     r_b = (1-u_b-v_b) r_b1 + u_b r_b2 + v_b r_b3
    !     grad r_b = (r_b2-r_b1)*grad u_b + (r_b3-r_b1)*grad v_b
    ! (m.m. for r_t)

    t_b = r_bottom/r_current
    t_t = r_top/r_current
    call ray_plane_intersection_gradient(xyz_bottom, xyz, t_b, grad_tuv_b)
    call ray_plane_intersection_gradient(xyz_top, xyz, t_t, grad_tuv_t)
    r_b_nodes = sqrt(sum(xyz_bottom**2, dim=1))
    r_t_nodes = sqrt(sum(xyz_top**2, dim=1))
    grad_r_new = (r_bottom-r_new)*grad_tuv_t(:,1) + (r_new-r_top)*grad_tuv_b(:,1)
    do i=2, size(xyz)
      grad_r_new = grad_r_new + (t_t-1.0)*(r_b_nodes(i)-r_b_nodes(1))*grad_tuv_b(:,i) &
                            & + (1.0-t_b)*(r_t_nodes(i)-r_t_nodes(1))*grad_tuv_t(:,i)
    end do
    grad_r_new = grad_r_new/(t_t-t_b)

    do i=1, size(xyz)
      do k=1, size(xyz)
        J(i,k) = J(i,k) + xyz(i)/r_current * grad_r_new(k)
      end do
    end do

  end function pop_out_jacobian

end module spherical_adaptivity
