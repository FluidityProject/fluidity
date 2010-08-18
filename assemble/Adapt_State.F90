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

module adapt_state_module
! these 5 need to be on top and in this order, so as not to confuse silly old intel compiler 
  use quadrature
  use elements
  use sparse_tools
  use fields
  use state_module
!
  use adaptivity_1d
  use adapt_integration, only : adapt_mesh_3d => adapt_mesh
  use adaptive_timestepping
  use checkpoint
  use diagnostic_fields_wrapper
  use discrete_properties_module
  use edge_length_module
  use eventcounter
  use field_options
  use global_parameters, only : OPTION_PATH_LEN, periodic_boundary_option_path, adaptivity_mesh_name, domain_bbox, topology_mesh_name
  use hadapt_extrude
  use hadapt_metric_based_extrude
  use halos
  use interpolation_manager
  use interpolation_module
  use mba_adapt_module
  use mba2d_integration
  use mba3d_integration
  use metric_assemble
  use anisotropic_gradation, only: use_anisotropic_gradation
  use parallel_tools
  use boundary_conditions
  use boundary_conditions_from_options
  use populate_state_module
  use project_metric_to_surface_module
  use reserve_state_module
  use sam_integration
  use tictoc  
  use timeloop_utilities
  use fields_halos
  use data_structures
  use detector_data_types
  use diagnostic_variables
  use intersection_finder_module
  use diagnostic_variables
  use pickers
#ifdef HAVE_ZOLTAN
  use zoltan_integration
#endif

  use mesh_files
  
  implicit none
  
  private
  
  public :: adapt_mesh, adapt_state, adapt_state_first_timestep
  public :: insert_metric_for_interpolation, extract_and_remove_metric, sam_options
  public :: adapt_state_module_check_options
  
  interface adapt_state
    module procedure adapt_state_single, adapt_state_multiple
  end interface adapt_state
  
contains

  subroutine adapt_mesh_simple(old_positions, metric, new_positions, node_ownership, force_preserve_regions, &
      lock_faces, allow_boundary_elements)
    type(vector_field), intent(in) :: old_positions
    type(tensor_field), intent(inout) :: metric
    type(vector_field), intent(out) :: new_positions
    integer, dimension(:), pointer, optional :: node_ownership
    logical, intent(in), optional :: force_preserve_regions
    type(integer_set), intent(in), optional :: lock_faces
    logical, intent(in), optional :: allow_boundary_elements

    assert(.not. mesh_periodic(old_positions))

    select case(old_positions%dim)
      case(1)
        call adapt_mesh_1d(old_positions, metric, new_positions, &
          & node_ownership = node_ownership, force_preserve_regions = force_preserve_regions)
      case(2)
        call adapt_mesh_mba2d(old_positions, metric, new_positions, &
          & force_preserve_regions=force_preserve_regions, lock_faces=lock_faces, &
          & allow_boundary_elements=allow_boundary_elements)
      case(3)
        if(have_option("/mesh_adaptivity/hr_adaptivity/adaptivity_library/libmba3d")) then
          assert(.not. present(lock_faces))
          call adapt_mesh_mba3d(old_positions, metric, new_positions, &
                             force_preserve_regions=force_preserve_regions)
        else
          call adapt_mesh_3d(old_positions, metric, new_positions, node_ownership = node_ownership, &
                             force_preserve_regions=force_preserve_regions, lock_faces=lock_faces)
        end if
      case default
        FLAbort("Mesh adaptivity requires a 1D, 2D or 3D mesh")
    end select
  end subroutine adapt_mesh_simple

  subroutine adapt_mesh_periodic(old_positions, metric, new_positions, force_preserve_regions)
    type(vector_field), intent(in) :: old_positions
    type(tensor_field), intent(inout) :: metric
    type(vector_field), intent(out) :: new_positions
    logical, intent(in), optional :: force_preserve_regions

    ! Periodic adaptivity variables
    integer :: no_bcs, bc, i, j, k, l
    type(integer_set) :: lock_faces, surface_ids
    type(vector_field) :: unwrapped_positions_A, unwrapped_positions_B, intermediate_positions
    integer, dimension(2) :: shape_option
    integer, dimension(:), allocatable :: surface_id
    type(tensor_field) :: unwrapped_metric_A, unwrapped_metric_B, intermediate_metric
    integer :: stat
    type(csr_sparsity), pointer :: eelist, nelist, periodic_eelist
    type(scalar_field) :: front_field
    integer :: ele, ele2, node, face
    integer, dimension(:), pointer :: neighbours, eles, periodic_neighbours, neighbours2
    integer :: new_physical_colour, new_aliased_colour
    integer :: face_count, existing_face_count
    integer, dimension(:), allocatable :: sndgln, boundary_ids, element_owners, fnodes
    integer :: floc
    integer, dimension(:), allocatable :: physical_colours, aliased_colours
    type(integer_set) :: nodes_to_move, periodic_nodes_to_move, extra_nodes_to_move
    real, dimension(:,:), allocatable:: aliased_positions, physical_positions
    character(len=OPTION_PATH_LEN) :: periodic_mapping_python
    type(integer_hash_table) :: aliased_to_new_node_number
    integer, dimension(:), pointer :: nodes, faces
    real, dimension(:, :), allocatable :: tmp_bbox
    type(integer_set) :: new_aliased_faces, new_physical_faces, old_physical_nodes
    type(integer_set) :: front_contained_nodes, front_face_nodes
    type(integer_set) :: other_surface_ids
    integer :: dim, sid

    integer, save :: delete_me = 1

    assert(mesh_periodic(metric))
    assert(metric%dim == old_positions%dim)
    dim = metric%dim

    no_bcs = option_count(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions')

    intermediate_positions = old_positions
    call incref(intermediate_positions)

    intermediate_metric = metric
    call incref(metric)

    ! As written, this is quadratic in the number of boundary conditions. I'm not too stressed
    ! about that

    do bc=0,no_bcs-1

      shape_option = option_shape(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(bc)//']/physical_boundary_ids')
      allocate(physical_colours(shape_option(1)))
      call get_option(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(bc)//']/physical_boundary_ids', physical_colours)
      shape_option = option_shape(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(bc)//']/aliased_boundary_ids')
      allocate(aliased_colours(shape_option(1)))
      call get_option(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(bc)//']/aliased_boundary_ids', aliased_colours)

      ! Step a). Unwrap the periodic input.
      unwrapped_positions_A = make_mesh_unperiodic_from_options(intermediate_positions, trim(periodic_boundary_option_path(dim)))
      call allocate(unwrapped_metric_A, unwrapped_positions_A%mesh, trim(metric%name))
      call remap_field(intermediate_metric, unwrapped_metric_A)

      ! We don't need the periodic mesh anymore
      call deallocate(intermediate_positions)
      call deallocate(intermediate_metric)

      ! Step b). Collect all the faces that need to be locked through the adapt
      call allocate(lock_faces)
      call allocate(surface_ids)
      call allocate(other_surface_ids)

      ! Collect all the relevant surface labels
      do j=0,no_bcs-1
        ! Physical ...
        shape_option = option_shape(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(j)//']/physical_boundary_ids')
        allocate(surface_id(shape_option(1)))
        call get_option(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(j)//']/physical_boundary_ids', surface_id)
        call insert(surface_ids, surface_id)
        if (j /= bc) then
          call insert(other_surface_ids, surface_id)
        end if
        deallocate(surface_id)

        ! and aliased
        shape_option = option_shape(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(j)//']/aliased_boundary_ids')
        allocate(surface_id(shape_option(1)))
        call get_option(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(j)//']/aliased_boundary_ids', surface_id)
        call insert(surface_ids, surface_id)
        if (j /= bc) then
          call insert(other_surface_ids, surface_id)
        end if
        deallocate(surface_id)
      end do

      ! With the relevant surface labels, loop through the mesh and fetch the information from them
      do j=1,surface_element_count(unwrapped_positions_A)
        if (has_value(surface_ids, surface_element_id(unwrapped_positions_A, j))) then
          call insert(lock_faces, j)
        end if
      end do
      call deallocate(surface_ids)

      ! Step c). Adapt the mesh, locking appropriately, and interpolate the metric
      ! mesh_0: incoming unwrapped mesh
      call vtk_write_fields("mesh", 0, position=unwrapped_positions_A, model=unwrapped_positions_A%mesh)
      call vtk_write_surface_mesh("surface", 0, unwrapped_positions_A)
      call adapt_mesh_simple(unwrapped_positions_A, unwrapped_metric_A, unwrapped_positions_B, & 
                & force_preserve_regions=force_preserve_regions, &
                & lock_faces=lock_faces, allow_boundary_elements=.true.)
      ! mesh_1: first adapted mesh
      call vtk_write_fields("mesh", 1, position=unwrapped_positions_B, model=unwrapped_positions_B%mesh)
      call vtk_write_surface_mesh("surface", 1, unwrapped_positions_B)
      call allocate(unwrapped_metric_B, unwrapped_positions_B%mesh, trim(metric%name))
      call linear_interpolation(unwrapped_metric_A, unwrapped_positions_A, unwrapped_metric_B, unwrapped_positions_B)
      call deallocate(lock_faces)
      call deallocate(unwrapped_positions_A)
      call deallocate(unwrapped_metric_A)

      ! Step d). Reperiodise
      intermediate_positions = make_mesh_periodic_from_options(unwrapped_positions_B, periodic_boundary_option_path(dim))
      intermediate_positions%mesh%name = "TmpMesh"
      intermediate_positions%mesh%option_path = periodic_boundary_option_path(dim)

      call allocate(intermediate_metric, intermediate_positions%mesh, trim(metric%name))
      call remap_field(unwrapped_metric_B, intermediate_metric, stat=stat)
      assert(stat /= REMAP_ERR_DISCONTINUOUS_CONTINUOUS)
      assert(stat /= REMAP_ERR_HIGHER_LOWER_CONTINUOUS)
      ! mesh_2: first adapted mesh, periodised
      call vtk_write_fields("mesh", 2, position=intermediate_positions, model=intermediate_positions%mesh)
      call vtk_write_surface_mesh("surface", 2, intermediate_positions)

      ! Step e). Advance a front in the new mesh using the unwrapped nelist from the aliased boundary
      ! until the front contains no nodes on the boundary; this forms the new cut
      nelist => extract_nelist(unwrapped_positions_B)
      shape_option = option_shape(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(bc)//']/aliased_boundary_ids')
      allocate(surface_id(shape_option(1)))
      call get_option(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(bc)//']/aliased_boundary_ids', surface_id)

      call allocate(nodes_to_move)
      call allocate(periodic_nodes_to_move)
      ! start from the nodes on the aliased boundary
      do j=1,surface_element_count(unwrapped_positions_B)
        if (any(surface_id == surface_element_id(unwrapped_positions_B, j))) then
          call insert(nodes_to_move, face_global_nodes(unwrapped_positions_B, j))
          call insert(periodic_nodes_to_move, face_global_nodes(intermediate_positions, j))
        end if
      end do
      deallocate(surface_id)
        
      ! if any nodes on the other periodic are in the elements directly
      ! adjacent move them too
      if (key_count(other_surface_ids)>0) then
        call allocate(extra_nodes_to_move)
        do i=1, key_count(nodes_to_move)
          node = fetch(nodes_to_move, i)
          neighbours => row_m_ptr(nelist, node)
          do j=1, size(neighbours)
            ele = neighbours(j)
            faces => ele_faces(unwrapped_positions_B, ele)
            do k=1, size(faces)
              face = faces(k)
              if (face>surface_element_count(unwrapped_positions_B)) cycle
              sid = surface_element_id(unwrapped_positions_B, face)
              if (has_value(other_surface_ids, sid)) then
                call insert(nodes_to_move, face_global_nodes(unwrapped_positions_B, face))
                ! work out face on the other side
                periodic_neighbours => ele_neigh(intermediate_positions, ele)
                face = ele_face(intermediate_positions, periodic_neighbours(k), ele)
                call insert(extra_nodes_to_move, face_global_nodes(unwrapped_positions_B, face))
              end if
            end do
          end do
        end do
        call insert(nodes_to_move, set2vector(extra_nodes_to_move))
        call deallocate(extra_nodes_to_move)
      end if
      
      ! the new cut now consists of the faces in the elements adjacent
      ! to the nodes to move that are not in "nodes_to_move" itself
      ! the nodes on the cut will in fact not be moved as they will retain
      ! their aliased position in the periodic position field
      call allocate(new_physical_faces)
      call allocate(new_aliased_faces)
      
      do i=1, key_count(nodes_to_move)
        node = fetch(nodes_to_move, i)
        neighbours => row_m_ptr(nelist, node)
        do j=1, size(neighbours)
           ele = neighbours(j)
           call insert(periodic_nodes_to_move, ele_nodes(intermediate_positions, ele))
           faces => ele_faces(unwrapped_positions_B, ele)
           neighbours2 => ele_neigh(unwrapped_positions_B, ele)
           do k=1, size(faces)
             face = faces(k)
             ele2 = neighbours2(k)
             if (ele2<=0) cycle
             if (.not. any(has_value(nodes_to_move,face_global_nodes(unwrapped_positions_B, face))) &
               .and. .not. any(has_value(nodes_to_move,ele_nodes(unwrapped_positions_B, ele2))) &
               ) then
               call insert(new_physical_faces, face)               
               ! opposite face becomes aliased
               face = ele_face(unwrapped_positions_B, ele2, ele)
               call insert(new_aliased_faces, face)
             end if
           end do
        end do
      end do
      call deallocate(nodes_to_move)
      
      ! don't want to move nodes on the new cut, as they'll retain
      ! their original before-aliasing position
      floc = face_loc(intermediate_positions, 1)
      allocate(fnodes(1:floc))
      do j=1, key_count(new_physical_faces)
        face = fetch(new_physical_faces, j)
        fnodes = face_global_nodes(intermediate_positions, face)
        do k=1, floc
          node = fnodes(k)
          if (has_value(periodic_nodes_to_move, node)) then
            call remove(periodic_nodes_to_move, node)
          end if
        end do
      end do
      deallocate(fnodes)
        
      ! now move the nodes in the periodic positions field
      allocate(aliased_positions(intermediate_positions%dim, key_count(periodic_nodes_to_move)))
      allocate(physical_positions(intermediate_positions%dim, key_count(periodic_nodes_to_move)))

      do j=1,key_count(periodic_nodes_to_move)
        node = fetch(periodic_nodes_to_move, j)
        aliased_positions(:, j) = node_val(intermediate_positions, node)
      end do

      call get_option(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(bc)//']/coordinate_map', periodic_mapping_python)
      call set_from_python_function(physical_positions, periodic_mapping_python, aliased_positions, time=0.0)

      do j=1,key_count(periodic_nodes_to_move)
        node = fetch(periodic_nodes_to_move, j)
        call set(intermediate_positions, node, physical_positions(:, j))
      end do

      deallocate(physical_positions)
      deallocate(aliased_positions)
      call deallocate(periodic_nodes_to_move)
      
      ! Step f). Colour those faces on either side of the new cut

      ! Choose a new colour that isn't used
      new_physical_colour = maxval(physical_colours)
      new_aliased_colour = maxval(aliased_colours)
      ! the only forward-compatible way of doing this is to fetch the information from
      ! the current faces into the primitive data structures, and then re-call add_faces
      ! fixme: for mixed meshes
      ! note well: we /lose/ the old faces, as we don't need that cut anymore
      ! (and can't retain it through the adapt anyway)

      ! first things first: find out how many faces we have
      face_count = key_count(new_physical_faces) + key_count(new_aliased_faces)
      do j=1,surface_element_count(intermediate_positions)
        if (.not. (any(surface_element_id(intermediate_positions, j) == physical_colours) .or. &
          & any(surface_element_id(intermediate_positions, j) == aliased_colours))) then
          face_count = face_count + 1
        end if
      end do
      
      allocate(boundary_ids(face_count))
      allocate(element_owners(face_count))
      allocate(sndgln(face_count * floc))

      ! fetch the existing information
      l = 0
      do j=1,surface_element_count(intermediate_positions)
        if (.not. (any(surface_element_id(intermediate_positions, j) == physical_colours) .or. &
          & any(surface_element_id(intermediate_positions, j) == aliased_colours))) then
        
          l = l + 1
          boundary_ids(l) = surface_element_id(intermediate_positions, j)
          element_owners(l) = face_ele(intermediate_positions, j)
          sndgln( (l-1)*floc + 1:l*floc ) = face_global_nodes(intermediate_positions, j)
        end if
      end do
        
      ! and now fetch the information for the faces we are adding
      
      do j=1, key_count(new_physical_faces)
        l = l + 1
        face = fetch(new_physical_faces, j)
        boundary_ids(l) = new_physical_colour
        element_owners(l) = face_ele(unwrapped_positions_B, face)
        sndgln( (l-1)*floc + 1:l*floc ) = face_global_nodes(intermediate_positions, face)
      end do
      call deallocate(new_physical_faces)
      
      do j=1, key_count(new_aliased_faces)
        l = l + 1
        face = fetch(new_aliased_faces, j)
        boundary_ids(l) = new_aliased_colour
        element_owners(l) = face_ele(unwrapped_positions_B, face)
        sndgln( (l-1)*floc + 1:l*floc ) = face_global_nodes(intermediate_positions, face)
      end do
      call deallocate(new_aliased_faces)

      assert(l == face_count)
      
      ! deallocate the old faces, and rebuild
      ! mesh_3: first adapted mesh, periodised with moving front nodes moved over
      call vtk_write_fields("mesh", 3, position=intermediate_positions, model=intermediate_positions%mesh)
      call vtk_write_surface_mesh("surface", 3, intermediate_positions)
      call deallocate_faces(intermediate_positions%mesh)
      call add_faces(intermediate_positions%mesh, sndgln=sndgln, element_owner=element_owners, boundary_ids=boundary_ids)
      intermediate_metric%mesh = intermediate_positions%mesh
      ! mesh_4: first adapted mesh, periodised with moving front nodes moved over and surface mesh relabeled
      call vtk_write_fields("mesh", 4, position=intermediate_positions, model=intermediate_positions%mesh)
      call vtk_write_surface_mesh("surface", 4, intermediate_positions)
      deallocate(sndgln)
      deallocate(element_owners)
      deallocate(boundary_ids)

      ! Step g). Unwrap again
      ! We need to fiddle with the options tree to mark the aliased and physical surface IDs appropriately

      unwrapped_positions_A = make_mesh_unperiodic_from_options(intermediate_positions, trim(periodic_boundary_option_path(dim)), & 
                                aliased_to_new_node_number=aliased_to_new_node_number, stat=stat)

      ! mesh_5: adapted once, moved up, ready to adapt again
      call vtk_write_fields("mesh", 5, position=unwrapped_positions_A, model=unwrapped_positions_A%mesh)
      call vtk_write_surface_mesh("surface", 5, unwrapped_positions_A)
      
      call allocate(unwrapped_metric_A, unwrapped_positions_A%mesh, trim(metric%name))
      call remap_field(intermediate_metric, unwrapped_metric_A)

      call deallocate(intermediate_positions)
      call deallocate(intermediate_metric)

      assert(has_faces(unwrapped_positions_A%mesh))

      call deallocate(unwrapped_positions_B)
      call deallocate(unwrapped_metric_B)

      ! Step h). Adapt again

      call allocate(lock_faces)
      call allocate(surface_ids)

      ! Collect all the relevant surface labels
      do j=0,no_bcs-1
        ! Physical ...
        shape_option = option_shape(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(j)//']/physical_boundary_ids')
        allocate(surface_id(shape_option(1)))
        call get_option(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(j)//']/physical_boundary_ids', surface_id)
        call insert(surface_ids, surface_id)
        deallocate(surface_id)

        ! and aliased
        shape_option = option_shape(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(j)//']/aliased_boundary_ids')
        allocate(surface_id(shape_option(1)))
        call get_option(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(j)//']/aliased_boundary_ids', surface_id)
        call insert(surface_ids, surface_id)
        deallocate(surface_id)
      end do

      ! With the relevant surface labels, loop through the mesh and fetch the information from them
      do j=1,surface_element_count(unwrapped_positions_A)
        if (has_value(surface_ids, surface_element_id(unwrapped_positions_A, j))) then
          call insert(lock_faces, j)
        end if
      end do
      call deallocate(surface_ids)

      call adapt_mesh_simple(unwrapped_positions_A, unwrapped_metric_A, unwrapped_positions_B, &
                  & force_preserve_regions=force_preserve_regions, &
                  & lock_faces=lock_faces, allow_boundary_elements=.true.)
      call deallocate(lock_faces)
      ! mesh_7: after second adapt
      call vtk_write_fields("mesh", 7, position=unwrapped_positions_B, model=unwrapped_positions_B%mesh)
      call vtk_write_surface_mesh("surface", 7, unwrapped_positions_B)

      call allocate(unwrapped_metric_B, unwrapped_positions_B%mesh, trim(metric%name))
      call linear_interpolation(unwrapped_metric_A, unwrapped_positions_A, unwrapped_metric_B, unwrapped_positions_B)
      call deallocate(unwrapped_positions_A)
      call deallocate(unwrapped_metric_A)

      ! Step i). Reperiodise for the next go around!
      intermediate_positions = make_mesh_periodic_from_options(unwrapped_positions_B, periodic_boundary_option_path(dim))
      intermediate_positions%mesh%option_path = periodic_boundary_option_path(dim)
      call allocate(intermediate_metric, intermediate_positions%mesh, trim(metric%name))
      call remap_field(unwrapped_metric_B, intermediate_metric, stat=stat)
      assert(stat /= REMAP_ERR_DISCONTINUOUS_CONTINUOUS)
      assert(stat /= REMAP_ERR_HIGHER_LOWER_CONTINUOUS)

      ! Step j). If the user has specified an inverse coordinate map, then let's loop through the nodes in the mesh
      ! and map them back to inside the bounding box of the domain
      if (have_option(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(bc)//']/inverse_coordinate_map')) then

        ! nodes_to_move stores the potential nodes to map.
        ! The rule is: we map any element which contains any node such that:
        ! the node's location is not inside the original bounding box, and
        ! the image of the node under the inverse mapping is inside the original bounding box

        ! So first let's find the nodes outside the bounding box
        call allocate(nodes_to_move)
        allocate(tmp_bbox(unwrapped_positions_B%dim, 2))
        do j=1, node_count(unwrapped_positions_B)
          tmp_bbox(:, 1) = node_val(unwrapped_positions_B, j)
          tmp_bbox(:, 2) = node_val(unwrapped_positions_B, j)
          if (.not. bbox_predicate(domain_bbox(1:unwrapped_positions_B%dim, :), tmp_bbox)) then
            call insert(nodes_to_move, j)
          end if
        end do

        allocate(aliased_positions(intermediate_positions%dim, key_count(nodes_to_move)))
        allocate(physical_positions(intermediate_positions%dim, key_count(nodes_to_move)))

        do j=1,key_count(nodes_to_move)
          physical_positions(:, j) = node_val(unwrapped_positions_B, fetch(nodes_to_move, j))
        end do

        ! and map those nodes
        call get_option(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(bc)//']/inverse_coordinate_map', periodic_mapping_python)
        call set_from_python_function(aliased_positions, periodic_mapping_python, physical_positions, time=0.0)

        ! now let's loop through those nodes, and if the image is inside the bounding box, mark
        ! the elements to map
        front_field = piecewise_constant_field(intermediate_positions%mesh, "AdvancingFront")
        call zero(front_field)
        nelist => extract_nelist(unwrapped_positions_B)
        eelist => extract_eelist(unwrapped_positions_B)
        periodic_eelist => extract_eelist(intermediate_positions)

        do j=1,key_count(nodes_to_move)
          tmp_bbox(:, 1) = aliased_positions(:, j)
          tmp_bbox(:, 2) = aliased_positions(:, j)
          if (bbox_predicate(domain_bbox(1:unwrapped_positions_B%dim, :), tmp_bbox)) then
            eles => row_m_ptr(nelist, fetch(nodes_to_move, j))
            do k=1,size(eles)
              call set(front_field, eles(k), 1.0)
            end do
          end if
        end do

        deallocate(physical_positions)
        deallocate(aliased_positions)
        call deallocate(nodes_to_move)
        deallocate(tmp_bbox)

        ! For doubly periodic, we need to make sure that the front is "periodic" in some sense.
        ! Otherwise, bad things can happen where the node on one side of the other BC wants to be
        ! mapped, but the other node doesn't, and there's no consistent solution.
        do ele=1,ele_count(front_field)
          if (node_val(front_field, ele) == 1.0) then
            neighbours => ele_neigh(intermediate_positions, ele)
            do j=1,size(neighbours)
              face = ele_face(intermediate_positions, ele, neighbours(j))
              if (face > surface_element_count(intermediate_positions)) cycle
              if (has_value(other_surface_ids, surface_element_id(intermediate_positions, face))) then
                call set(front_field, neighbours(j), 1.0)
              end if
            end do
          end if
        end do

        ! mesh_8: after second adapt, showing elements to be moved back into bounding box
        call vtk_write_fields("mesh", 8, position=unwrapped_positions_B, model=unwrapped_positions_B%mesh, sfields=(/front_field/))
        call vtk_write_surface_mesh("surface", 8, unwrapped_positions_B)

        ! mesh_8: same thing periodised
        call vtk_write_fields("mesh", 9, position=intermediate_positions, model=intermediate_positions%mesh, sfields=(/front_field/))
        call vtk_write_surface_mesh("surface", 9, intermediate_positions)

        ! OK. Now we know which elements we are mapping, it is very similar to the shuffling
        ! around we did earlier. The two main subtasks are to
        ! a) update the positions of the periodic nodes appropriately
        ! b) change the faces of the periodic mesh

        ! First thing: let's build two sets
        ! that store the current lists of physical and aliased faces.

        call allocate(new_aliased_faces)
        call allocate(new_physical_faces)
        call allocate(old_physical_nodes)
        existing_face_count = 0
        do j=1,surface_element_count(intermediate_positions)
          if (surface_element_id(intermediate_positions, j) == new_physical_colour) then
            call insert(new_physical_faces, j)
            call insert(old_physical_nodes, face_global_nodes(intermediate_positions, j))
          else if (surface_element_id(intermediate_positions, j) == new_aliased_colour) then
            call insert(new_aliased_faces, j)
          else
            existing_face_count = existing_face_count + 1
          end if
        end do

        call allocate(front_contained_nodes)
        call allocate(front_face_nodes)
        do ele=1,ele_count(front_field)
          if (node_val(front_field, ele) == 1.0) then
            call insert(front_contained_nodes, ele_nodes(intermediate_positions, ele))
            neighbours => row_m_ptr(eelist, ele)
            periodic_neighbours => row_m_ptr(periodic_eelist, ele)
            faces => ele_faces(intermediate_positions, ele)
            neighbourloop: do k=1,size(neighbours)
              j = neighbours(k)
              face = faces(k)

              if (has_value(new_physical_faces, face)) then
                call insert(front_face_nodes, face_global_nodes(intermediate_positions, face))
                call remove(new_physical_faces, face)
                call remove(new_aliased_faces, ele_face(intermediate_positions, periodic_neighbours(k), ele))
              end if

              if (j > 0) then
                if (node_val(front_field, j) /= 1.0) then
                  face = ele_face(intermediate_positions, ele, j)
                  call insert(new_aliased_faces, face)

                  face = ele_face(intermediate_positions, j, ele)
                  call insert(new_physical_faces, face)
                end if
              end if

            end do neighbourloop

            nodes => ele_nodes(intermediate_positions, ele)
            do k=1,size(nodes)
              if (has_value(old_physical_nodes, nodes(k))) then
                call insert(front_face_nodes, nodes(k))
              end if
            end do

          end if
        end do

        call deallocate(old_physical_nodes)

        ! Now pack into the primitive data structures

        allocate(boundary_ids(existing_face_count + 2 * key_count(new_physical_faces)))
        allocate(element_owners(existing_face_count + 2 * key_count(new_physical_faces)))
        allocate(sndgln(floc * (existing_face_count + 2 * key_count(new_physical_faces))))

        l = 1
        do j=1,surface_element_count(intermediate_positions)
          if (surface_element_id(intermediate_positions, j) == new_physical_colour) then
            cycle
          else if (surface_element_id(intermediate_positions, j) == new_aliased_colour) then
            cycle
          else
            boundary_ids(l) = surface_element_id(intermediate_positions, j)
            element_owners(l) = face_ele(intermediate_positions, j)
            sndgln( (l-1)*floc + 1:l*floc ) = face_global_nodes(intermediate_positions, j)
            l = l + 1
          end if
        end do
        assert(l == existing_face_count + 1)

        do j=1,key_count(new_physical_faces)
          face = fetch(new_physical_faces, j)
          boundary_ids(l) = new_physical_colour
          element_owners(l) = face_ele(intermediate_positions, face)
          sndgln( (l-1)*floc + 1:l*floc ) = face_global_nodes(intermediate_positions, face)
          l = l + 1

          face = fetch(new_aliased_faces, j)
          boundary_ids(l) = new_aliased_colour
          element_owners(l) = face_ele(intermediate_positions, face)
          sndgln( (l-1)*floc + 1:l*floc ) = face_global_nodes(intermediate_positions, face)
          l = l + 1
        end do
        assert(l == size(boundary_ids) + 1)
        !call vtk_write_internal_face_mesh("surface", 11, unwrapped_positions_B, face_sets=(/new_physical_faces, new_aliased_faces/))
        assert(key_count(new_physical_faces) == key_count(new_aliased_faces))

        call deallocate(new_physical_faces)
        call deallocate(new_aliased_faces)


        ! deallocate the old faces, and rebuild
        call deallocate_faces(intermediate_positions%mesh)
        call add_faces(intermediate_positions%mesh, sndgln=sndgln, element_owner=element_owners, boundary_ids=boundary_ids)
        call vtk_write_surface_mesh("surface", 12, intermediate_positions)
        intermediate_metric%mesh = intermediate_positions%mesh
        deallocate(sndgln)
        deallocate(element_owners)
        deallocate(boundary_ids)

        call set_minus(nodes_to_move, front_contained_nodes, front_face_nodes)
        call deallocate(front_contained_nodes)
        call deallocate(front_face_nodes)

        allocate(aliased_positions(intermediate_positions%dim, key_count(nodes_to_move)))
        allocate(physical_positions(intermediate_positions%dim, key_count(nodes_to_move)))

        do j=1,key_count(nodes_to_move)
          physical_positions(:, j) = node_val(intermediate_positions, fetch(nodes_to_move, j))
        end do

        call set_from_python_function(aliased_positions, periodic_mapping_python, physical_positions, time=0.0)

        allocate(tmp_bbox(intermediate_positions%dim, 2))
        do j=1,key_count(nodes_to_move)
          tmp_bbox(:, 1) = aliased_positions(:, j)
          tmp_bbox(:, 2) = aliased_positions(:, j)
          call set(intermediate_positions, fetch(nodes_to_move, j), aliased_positions(:, j))
        end do
        deallocate(tmp_bbox)

        deallocate(physical_positions)
        deallocate(aliased_positions)
        call deallocate(nodes_to_move)

        call deallocate(front_field)
      end if

      call deallocate(unwrapped_positions_B)
      call deallocate(unwrapped_metric_B)

      deallocate(physical_colours)
      deallocate(aliased_colours)

      ! mesh_10: final mesh
      call vtk_write_fields("mesh", 10, position=intermediate_positions, model=intermediate_positions%mesh)
      call vtk_write_surface_mesh("surface", 10, intermediate_positions)

      call deallocate(other_surface_ids)
    end do

    new_positions = intermediate_positions
    new_positions%option_path = old_positions%option_path
    new_positions%mesh%option_path = old_positions%mesh%option_path
    new_positions%name = old_positions%name
    new_positions%mesh%name = old_positions%mesh%name

    call deallocate(intermediate_metric)

    call vtk_write_fields("adapted_mesh", delete_me, position=new_positions, model=new_positions%mesh)

    unwrapped_positions_A = make_mesh_unperiodic_from_options(intermediate_positions, trim(periodic_boundary_option_path(dim)))
    call vtk_write_fields("adapted_mesh_unwrapped", delete_me, position=unwrapped_positions_A, model=unwrapped_positions_A%mesh)
    call vtk_write_surface_mesh("adapted_surface_unwrapped", delete_me, unwrapped_positions_A)
    call deallocate(unwrapped_positions_A)
    delete_me = delete_me + 1

  end subroutine adapt_mesh_periodic
  
  subroutine adapt_mesh(old_positions, metric, new_positions, node_ownership, force_preserve_regions)
    !!< A wrapper to select the appropriate adapt_mesh routine.
    !!< If the input is periodic, then apply the algorithm for adapting periodic meshes.
    
    type(vector_field), intent(in) :: old_positions
    type(tensor_field), intent(inout) :: metric
    type(vector_field), intent(out) :: new_positions
    integer, dimension(:), pointer, optional :: node_ownership
    logical, intent(in), optional :: force_preserve_regions

#ifdef DDEBUG
    if(present(node_ownership)) then
      assert(.not. associated(node_ownership))
    end if
#endif

    ! Periodic case
    if (mesh_periodic(old_positions)) then
      call adapt_mesh_periodic(old_positions, metric, new_positions, force_preserve_regions=force_preserve_regions)
    ! Nonperiodic case
    else
      call adapt_mesh_simple(old_positions, metric, new_positions, node_ownership=node_ownership, force_preserve_regions=force_preserve_regions)
    end if
  end subroutine adapt_mesh

  subroutine adapt_state_single(state, metric, initialise_fields)
   
    type(state_type), intent(inout) :: state
    type(tensor_field), intent(inout) :: metric
    !! If present and .true., initialise fields rather than interpolate them
    logical, optional, intent(in) :: initialise_fields
    
    type(state_type), dimension(1) :: states
    
    states = (/state/)
    call adapt_state(states, metric, initialise_fields = initialise_fields)
    state = states(1)
    
  end subroutine adapt_state_single

  subroutine adapt_state_multiple(states, metric, initialise_fields)

    type(state_type), dimension(:), intent(inout) :: states
    type(tensor_field), intent(inout) :: metric
    !! If present and .true., initialise fields rather than interpolate them
    logical, optional, intent(in) :: initialise_fields
    
    call tictoc_clear(TICTOC_ID_SERIAL_ADAPT)
    call tictoc_clear(TICTOC_ID_DATA_MIGRATION)
    call tictoc_clear(TICTOC_ID_DATA_REMAP)
    call tictoc_clear(TICTOC_ID_ADAPT)
    
    call tic(TICTOC_ID_ADAPT)
    
    call adapt_state_internal(states, metric, initialise_fields = initialise_fields)
    
    call toc(TICTOC_ID_ADAPT)
       
    call tictoc_report(2, TICTOC_ID_SERIAL_ADAPT)
    call tictoc_report(2, TICTOC_ID_DATA_MIGRATION)
    call tictoc_report(2, TICTOC_ID_DATA_REMAP)
    call tictoc_report(2, TICTOC_ID_ADAPT)
    
  end subroutine adapt_state_multiple
  
  subroutine adapt_state_first_timestep(states)
    !!< Subroutine to adapt the supplied states at the simulation start
    
    type(state_type), dimension(:), intent(inout) :: states
    
    character(len = *), parameter :: base_path = "/mesh_adaptivity/hr_adaptivity/adapt_at_first_timestep"
    integer :: adapt_iterations, i
    type(mesh_type), pointer :: old_mesh
    type(tensor_field) :: metric
    type(vector_field), pointer :: output_positions
    real :: dt
    
    ewrite(1, *) "In adapt_state_first_timestep"
    
    call get_option(trim(base_path) // "/number_of_adapts", adapt_iterations)
    
    do i = 1, adapt_iterations
      ewrite(2, "(a,i0,a,i0)") "Performing first timestep adapt ", i, " of ", adapt_iterations
      
      ! Recalculate diagnostics, as error metric formulations may need them
      call allocate_and_insert_auxilliary_fields(states)
      call copy_to_stored_values(states,"Old")
      call copy_to_stored_values(states,"Iterated")
      call relax_to_nonlinear(states)
      
      call calculate_diagnostic_variables(states)
    
      call enforce_discrete_properties(states)
      if(have_option("/timestepping/adaptive_timestep/at_first_timestep")) then
        ! doing this here helps metric advection get the right amount of advection
        call get_option("/timestepping/timestep", dt)
        call calc_cflnumber_field_based_dt(states, dt, force_calculation = .true.)
        call set_option("/timestepping/timestep", dt)
      end if
      
      ! Form the new metric
      old_mesh => extract_mesh(states(1), "CoordinateMesh")
      call allocate(metric, old_mesh, "ErrorMetric")
      call assemble_metric(states, metric)
      
      ! Adapt state, initialising fields from the options tree rather than
      ! interpolating them
      call adapt_state(states, metric, initialise_fields = .true.)
    end do
    
    if(have_option(trim(base_path) // "/output_adapted_mesh")) then
      output_positions => extract_vector_field(states(1), "Coordinate")

      if(isparallel()) then
        call write_mesh_files(parallel_filename("first_timestep_adapted_mesh"), output_positions)
        call write_halos("first_timestep_adapted_mesh", output_positions%mesh)
      else
        call write_mesh_files("first_timestep_adapted_mesh", output_positions)
      end if

    end if
    
    ewrite(1, *) "Exiting adapt_state_first_timestep"
  
  end subroutine adapt_state_first_timestep
      
  subroutine adapt_state_internal(states, metric, initialise_fields)
    !!< Adapt the supplied states according to the supplied metric. In parallel,
    !!< additionally re-load-balance with libsam. metric is deallocated by this
    !!< routine. Based on adapt_state_2d.

    type(state_type), dimension(:), intent(inout) :: states
    type(tensor_field), intent(inout) :: metric
    !! If present and .true., re-initialise fields with their initial condition.
    !! This means that the fields are not interpolated but rather reinitialise
    !! according to the specified initial condition in the options tree, except
    !! if these fields are initialised from_file (checkpointed).
    logical, optional, intent(in) :: initialise_fields
    
    character(len = FIELD_NAME_LEN) :: metric_name
    integer :: i, j, max_adapt_iteration
    integer, dimension(:), pointer :: node_ownership
    type(state_type), dimension(size(states)) :: interpolate_states
    type(mesh_type), pointer :: old_linear_mesh
    type(vector_field) :: old_positions, new_positions
    logical :: vertical_only

    ! Vertically structured adaptivity stuff
    type(vector_field) :: extruded_positions
    type(tensor_field) :: full_metric
    logical :: vertically_structured_adaptivity

    ewrite(1, *) "In adapt_state_internal"
    
    nullify(node_ownership)
    
    max_adapt_iteration = adapt_iterations()

    vertically_structured_adaptivity = have_option( &
     &  "/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity")
    vertical_only = have_option(&
        & "/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity/inhomogenous_vertical_resolution/adapt_in_vertical_only")

    ! Don't need to strip the level 2 halo with Zoltan .. in fact, we don't want to
#ifndef HAVE_ZOLTAN
    if(isparallel()) then
      ! In parallel, strip off the level 2 halo (expected by libsam). The level
      ! 2 halo is restored on the final adapt iteration by libsam.
      call strip_level_2_halo(states, metric, initialise_fields=initialise_fields)
    end if
#endif

    do i = 1, max_adapt_iteration
      if(max_adapt_iteration > 1) then
        ewrite(2, "(a,i0)") "Performing adapt ", i
      end if
      
      ! Select mesh to adapt. Has to be linear and continuous.
      ! For vertically_structured_adaptivity, this is the horizontal mesh!
      call find_mesh_to_adapt(states(1), old_linear_mesh)
      ewrite(2, *) "External mesh to be adapted: " // trim(old_linear_mesh%name)
      if (mesh_periodic(old_linear_mesh)) then
        old_positions = extract_vector_field(states(1), trim(old_linear_mesh%name) // "Coordinate")
        call incref(old_positions)
      else
        ! Extract the mesh field to be adapted (takes a reference)
        old_positions = get_coordinate_field(states(1), old_linear_mesh)
      end if
      ewrite(2, *) "Mesh field to be adapted: " // trim(old_positions%name)
      
      call prepare_vertically_structured_adaptivity(states, metric, full_metric, &
                                                    extruded_positions, old_positions)
      
      call initialise_boundcount(old_linear_mesh, old_positions)
     
      do j = 1, size(states)
        ! Reference fields to be interpolated in interpolate_states
        ! (if initialise_fields then leave out those fields that can be reinitialised)
        call select_fields_to_interpolate(states(j), interpolate_states(j), &
          & first_time_step = initialise_fields)
      end do

      do j = 1, size(states)
        call deallocate(states(j))
      end do

      if(isparallel()) then
        ! Update the fields to be interpolated, just in case
        call halo_update(interpolate_states)
      end if 
      
      ! Before we start allocating any new objects we tag all references to
      ! current objects before the adapt so we can later on check they have all
      ! been deallocated
      call tag_references()
      
      ! Generate a new mesh field based on the current mesh field and the input
      ! metric
      if (.not. vertical_only) then
        call adapt_mesh(old_positions, metric, new_positions, node_ownership = node_ownership, &
          & force_preserve_regions=initialise_fields)
      else
        call allocate(new_positions,old_positions%dim,old_positions%mesh,name=trim(old_positions%name))
        call set(new_positions,old_positions)
      end if
      ! We're done with old_positions, so we may deallocate it
      call deallocate(old_positions)

      ! Insert the new mesh field and linear mesh into all states
      call insert(states, new_positions%mesh, name = new_positions%mesh%name)
      call insert(states, new_positions, name = new_positions%name)
      
      call perform_vertically_inhomogenous_step(states, new_positions, full_metric, extruded_positions)

      ! Insert meshes from reserve states
      call restore_reserved_meshes(states)
      ! Next we recreate all derived meshes
      call insert_derived_meshes(states)
      
      if(vertically_structured_adaptivity) then
        call deallocate(metric)
        call deallocate(new_positions)
        
        metric = full_metric
        new_positions = get_coordinate_field(states(1), extract_mesh(states(1), topology_mesh_name))
      end if

      if(isparallel()) then
        ! Update the detector element ownership data
        call search_for_detectors(detector_list, new_positions)
      end if

      ! Then reallocate all fields 
      call allocate_and_insert_fields(states)
      ! Insert fields from reserve states
      call restore_reserved_fields(states)
      ! Add on the boundary conditions again
      call populate_boundary_conditions(states)
      ! Set their values
      call set_boundary_conditions_values(states)
      
      if(i < max_adapt_iteration .or. isparallel()) then
        ! If there are remaining adapt iterations, or we will be calling
        ! sam_drive or zoltan_drive, insert the old metric into interpolate_states(1) and a
        ! new metric into states(1), for interpolation
        call insert_metric_for_interpolation(metric, new_positions%mesh, &
                                             interpolate_states(1), states(1), &
                                             metric_name = metric_name)
      end if
      ! We're done with the old metric, so we may deallocate it / drop our
      ! reference
      call deallocate(metric)
      ! We're done with the new_positions, so we may drop our reference
      call deallocate(new_positions) 
      
      
      ! Interpolate fields
      if(associated(node_ownership)) then
        call interpolate(interpolate_states, states, map = node_ownership)
      else
        call interpolate(interpolate_states, states)
      end if
      
      ! Deallocate the old fields used for interpolation, referenced in
      ! interpolate_states
      do j = 1, size(states)
        call deallocate(interpolate_states(j))
      end do
      if(associated(node_ownership)) then
        ! Deallocate the node ownership mapping
        deallocate(node_ownership)
        nullify(node_ownership)
      end if
      
      if(i < max_adapt_iteration .or. isparallel()) then
        ! If there are remaining adapt iterations, extract the new metric for
        ! the next adapt iteration. If we will be calling sam_drive, always
        ! extract the new metric.
        metric = extract_and_remove_metric(states(1), metric_name)
      end if
      
      if(present_and_true(initialise_fields)) then
        ! Reinitialise the prognostic fields (where possible)
        call initialise_prognostic_fields(states)
        ! Prescribed fields are recalculated
        ! NOTE: we don't have exclude_interpolated, as the only prescribed
        ! fields that are interpolated are from_file which will be skipped
        ! anyway because initial_mesh = .false., and the routine  doesn't know
        ! we're not interpolating other prescribed fields with interpolation
        ! options
        call set_prescribed_field_values(states)
      else
        ! Prescribed fields are recalculated (except those with interpolation 
        ! options)
        call set_prescribed_field_values(states, exclude_interpolated = .true.)
      end if
      
      ! If strong bc or weak that overwrite then enforce the bc on the fields
      call set_dirichlet_consistent(states)
      ! Insert aliased fields in state
      call alias_fields(states)      

      zoltan_drive_call=.false.

      if(isparallel()) then
#ifdef HAVE_ZOLTAN
        ! Re-load-balance using zoltan
        call zoltan_drive(states, i, metric = metric)
        zoltan_drive_call=.true.         
#else
        ! Re-load-balance using libsam
        call sam_drive(states, sam_options(i, max_adapt_iteration), metric = metric)
#endif
        if(i == max_adapt_iteration) then
          ! On the last adapt iteration the metric was interpolated
          ! only for re-load-balancing, hence it must be deallocated
          call deallocate(metric)
        end if
      end if

!     call deallocate(new_positions)
      
      if(vertical_only) then
        ewrite(2,*) "Using vertical_only adaptivity, so skipping the printing of references"
      else if (no_reserved_meshes()) then
        ewrite(2, *) "Tagged references remaining:"
        call print_tagged_references(0)
      else
        ewrite(2, *) "There are reserved meshes, so skipping printing of references."
      end if
 
      call write_adapt_state_debug_output(states, i, max_adapt_iteration, &
        & initialise_fields = initialise_fields)   
      
      call incrementeventcounter(EVENT_ADAPTIVITY)
      call incrementeventcounter(EVENT_MESH_MOVEMENT)
    end do

    ewrite(1, *) "Exiting adapt_state_internal"

  end subroutine adapt_state_internal
  
  subroutine insert_metric_for_interpolation(metric, new_mesh, old_state, new_state, metric_name)
    !!< Insert the old metric into old_states and a new metric into new_states, 
    !!< for interpolation
    
    type(tensor_field), intent(in) :: metric
    type(mesh_type), intent(in) :: new_mesh
    type(state_type), intent(inout) :: old_state
    type(state_type), intent(inout) :: new_state
    character(len = *), optional, intent(out) :: metric_name
    
    type(tensor_field) :: new_metric
    
    assert(.not. has_tensor_field(old_state, metric%name))
    call insert(old_state, metric, metric%name)

    call allocate(new_metric, new_mesh, metric%name)
    assert(.not. has_tensor_field(new_state, new_metric%name))
    call insert(new_state, new_metric, new_metric%name)

    if(present(metric_name)) metric_name = new_metric%name
    
    call deallocate(new_metric)
    
  end subroutine insert_metric_for_interpolation
  
  function extract_and_remove_metric(state, metric_name) result(metric)
    !!< Extract and remove the metric from the supplied state. metric takes
    !!< a reference in this routine.
    
    type(state_type), intent(inout) :: state
    character(len = *), intent(in) :: metric_name
    
    type(tensor_field) :: metric
    
    type(tensor_field), pointer :: metric_ptr
    
    ! Extract the metric
    metric_ptr => extract_tensor_field(state, metric_name)
    metric = metric_ptr
#ifdef DDEBUG
    ! Check the metric
    call check_metric(metric)
#endif
    ! Take a reference to the metric
    call incref(metric)
    ! and remove it from state
    call remove_tensor_field(state, metric%name)
    
  end function extract_and_remove_metric
  
  function adapt_iterations()
    !!< Return the number of adapt / re-load-balance iterations
    
    integer :: adapt_iterations
    
    if(isparallel()) then
      adapt_iterations = 3
    else
      adapt_iterations = 1
    end if
    
  end function adapt_iterations
  
  pure function sam_options(adapt_iteration, max_adapt_iteration)
    !!< Return sam options array
    
    integer, intent(in) :: adapt_iteration
    integer, intent(in) :: max_adapt_iteration
    
    integer, dimension(10) :: sam_options
    
    sam_options = 0
    
    ! Target number of partitions - 0 indicates size of MPI_COMM_WORLD
    sam_options(1) = 0
    
    ! Graph partitioning options:    
    !sam_options(2) = 1  ! Clean partitioning to optimise the length of the 
                         ! interface boundary.    
    if(adapt_iteration < max_adapt_iteration) then
      ! Diffusive method -- fast partitioning, small partition movement
      ! thus edges are weighed to avoid areas of high activity. 
      ! sam_options(2) = 2  ! Local diffusion
      sam_options(2) = 3  ! Directed diffusion
    else
      ! Clean partitioning to optimise the length of the interface boundary.
      ! This partitioning is then remapped  onto the original partitioning to
      ! maximise overlap and therefore the volume of data migration.
      sam_options(2) = 4
    end if
    
    ! Heterogerious options (disabled)
    sam_options(3) = 1

    ! Node and edge weight options
    if(adapt_iteration < max_adapt_iteration) then
      ! Node weights are based on an estimate of the density of nodes in the
      ! region of a node after adaption
      sam_options(4) = 2
      ! Calculate edge weights as being the maximum length in metric space of
      ! any element surrounding the edge. This should give high weights to
      ! elements that are likely to be involved in adaption.
      sam_options(5) = 2
      ! Mixed formulation options
      sam_options(6) = 1 ! Disabled
                         ! Do not restore the level 2 halo
    else
      ! No node weights
      sam_options(4) = 1
      ! No edge weights
      sam_options(5) = 1
      ! Mixed formulation options
      sam_options(6) = 2 ! Enabled
                         ! Restore the level 2 halo
    end if

  end function sam_options

  subroutine prepare_vertically_structured_adaptivity(states, metric, full_metric, &
                                                      extruded_positions, old_positions)
    type(state_type), dimension(:), intent(inout) :: states
    ! the metric will be collapsed, and the uncollapsed full_metric stored in full_metric
    type(tensor_field), intent(inout) :: metric, full_metric
    type(vector_field), intent(inout) :: extruded_positions
    ! old positions of the horizontal mesh
    type(vector_field), intent(in) :: old_positions

    integer, save:: adaptcnt=0
    logical :: vertically_structured_adaptivity
    logical :: vertically_inhomogenous_adaptivity
    logical :: include_bottom_metric
    logical :: split_gradation

    type(scalar_field) :: edge_lengths
      
    vertically_structured_adaptivity = have_option( &
     &  "/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity")
    vertically_inhomogenous_adaptivity = have_option( &
     &  "/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity/inhomogenous_vertical_resolution")
    include_bottom_metric = have_option( &
     &  "/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity/include_bottom_metric")
    split_gradation = have_option( &
     &  "/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity/split_gradation")

    if (vertically_structured_adaptivity) then
      ! project full mesh metric to horizontal surface mesh metric
      full_metric=metric
      call project_metric_to_surface(full_metric, old_positions, metric)
      
      ! include the components of the full_metric that are tangential to the bathymetry 
      ! in the surface metric
      ! state only required for DistanceToBottom and Coordinate so states(1) is fine
      if(include_bottom_metric) call incorporate_bathymetric_metric(states(1), full_metric, old_positions, metric)
      
      if(split_gradation) then
        ! apply gradation just to the horizontal metric for now
        call apply_horizontal_gradation(states(1), metric, full_metric, old_positions)
      end if
      
      ! apply limiting to enforce maximum number of nodes
      call limit_metric(old_positions, metric)
      if (have_option('/mesh_adaptivity/hr_adaptivity/debug/write_metric_stages')) then
        call allocate(edge_lengths, metric%mesh, "EdgeLengths")
        call get_edge_lengths(metric, edge_lengths)
        call vtk_write_fields('horizontal_metric', adaptcnt, &
          old_positions, old_positions%mesh, &
          sfields=(/ edge_lengths /), tfields=(/ metric /) )
        adaptcnt=adaptcnt+1
        call deallocate(edge_lengths)
      end if
      
      if (vertically_inhomogenous_adaptivity) then
         ! we need the position field later on for vertical adaptivity
         ! this takes a reference so that it's prevented from the big deallocate in adapt_state
         extruded_positions = get_coordinate_field(states(1), full_metric%mesh)
      end if
    end if
    
  end subroutine prepare_vertically_structured_adaptivity

  subroutine perform_vertically_inhomogenous_step(states, new_positions, full_metric, extruded_positions)
    type(state_type), intent(inout), dimension(:) :: states
    type(vector_field), intent(inout) :: new_positions
    type(tensor_field), intent(inout) :: full_metric
    type(vector_field), intent(inout) :: extruded_positions

    logical :: vertically_inhomogenous_adaptivity

    type(vector_field) :: background_positions
    type(tensor_field) :: background_full_metric

    logical :: split_gradation

    vertically_inhomogenous_adaptivity = have_option( &
     &  "/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity/inhomogenous_vertical_resolution")
    split_gradation = have_option( &
     &  "/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity/split_gradation")

    if (vertically_inhomogenous_adaptivity) then
       ! first we create a background mesh: this is an extrusion of 
       ! the new horizontal mesh using the layer depths specified
       ! for the initial extruded mesh
       call extrude(new_positions, full_metric%mesh%option_path, &
          background_positions)
       
       ! now map the old full metric on this background mesh
       call allocate(background_full_metric, &
          background_positions%mesh, name="BackgroundFullMetric")
       ! get old extruded positions (takes reference)
       ! temp. fix: old and new metric need same name for linear_interpolation -ask Patrick
       background_full_metric%name=full_metric%name
       call linear_interpolation( full_metric, extruded_positions, &
          & background_full_metric, background_positions)
       background_full_metric%name="BackgroundFullMetric"
       ! we can do away with the old positions now
       call deallocate(extruded_positions)
       
       if(split_gradation) then
          call apply_vertical_gradation(states(1), background_full_metric, background_positions, new_positions)
       end if
       
       ! extrude with adaptivity, computes new extruded_positions
       call metric_based_extrude(new_positions, background_positions, &
          background_full_metric, extruded_positions)
       ! and we're done with the background stuff
       call deallocate(background_full_metric)
       call deallocate(background_positions)
       
       ! insert the new positions in state:
       ! give it a generic temporary name, so that it'll be picked up and
       ! adjusted by insert_derived meshes later on:
       extruded_positions%name="AdaptedExtrudedPositions"
       call insert(states, extruded_positions, name="AdaptedExtrudedPositions")
       ! and drop our reference:
       call deallocate(extruded_positions)
    end if
  end subroutine perform_vertically_inhomogenous_step
  
  subroutine write_adapt_state_debug_output(states, adapt_iteration, max_adapt_iteration, initialise_fields)
    !!< Diagnostic output for mesh adaptivity
  
    type(state_type), dimension(:), intent(in) :: states
    !! The current iteration the adapt-re-load-balance loop
    integer, intent(in) :: adapt_iteration
    !! The total number of iterations to be performed in the
    !! adapt-re-load-balance loop
    integer, intent(in) :: max_adapt_iteration
    !! If present and .true., initialise fields rather than interpolate them
    logical, optional, intent(in) :: initialise_fields
    
    character(len = FIELD_NAME_LEN) :: file_name
    character(len = *), parameter :: base_path = "/mesh_adaptivity/hr_adaptivity/debug"
    integer :: max_output, stat
    type(mesh_type), pointer :: mesh
    type(vector_field) :: positions
    
    integer, save :: cp_no = 0, mesh_dump_no = 0, state_dump_no = 0
    
    if(.not. have_option(base_path)) then
      ! No debug output options
      return
    end if
    
    if(have_option(base_path // "/write_adapted_mesh")) then
      ! Debug mesh output. These are output on every adapt iteration.

      file_name = adapt_state_debug_file_name("adapted_mesh", mesh_dump_no, adapt_iteration, max_adapt_iteration)
      call find_mesh_to_adapt(states(1), mesh)
      positions = get_coordinate_field(states(1), mesh)
      call write_mesh_files(file_name, positions)
      if(isparallel()) then
        file_name = adapt_state_debug_file_name("adapted_mesh", mesh_dump_no, adapt_iteration, max_adapt_iteration, &
                                                add_parallel = .false.)  ! parallel extension is added by write_halos
        call write_halos(file_name, positions%mesh)
      end if
      call deallocate(positions)
      
      mesh_dump_no = mesh_dump_no + 1
    end if
    
    if(have_option(base_path // "/write_adapted_state")) then   
      ! Debug vtu output. These are output on every adapt iteration.

      file_name = adapt_state_debug_file_name("adapted_state", state_dump_no, adapt_iteration, max_adapt_iteration, add_parallel = .false.)   
      call vtk_write_state(file_name, state = states)
      
      state_dump_no = state_dump_no + 1
    end if
    
    if(adapt_iteration == max_adapt_iteration .and. have_option(base_path // "/checkpoint")) then
      ! Debug checkpointing. These are only output on the final adapt iteration.
    
      if(present_and_true(initialise_fields)) then
        ! If we're adapting with field initialisation rather than interpolation
        ! then we probably don't want to overwrite the field initialisation
        ! options by checkpointing, as any subsequent adapt with field
        ! initialisation will read (and consistently interpolate) the debug
        ! checkpoint. Applies to first timestep adapts.        
        ewrite(1, *) "Adapt checkpoint skipped, as adapt performed with field initialisation"
      else
        ewrite(1, "(a,i0)") "Performing adapt checkpoint ", cp_no
        
        call checkpoint_simulation(states, postfix = "adapt_checkpoint", cp_no = cp_no)
        
        cp_no = cp_no + 1
        
        call get_option(base_path // "/checkpoint/max_checkpoint_count", max_output, stat = stat)
        if(stat == SPUD_NO_ERROR) cp_no = modulo(cp_no, max_output)
      end if

    end if
    
  contains
  
    function adapt_state_debug_file_name(base_name, dump_no, adapt_iteration, max_adapt_iteration, add_parallel) result(file_name)
      !!< Form an adapt diagnostic output filename

      !! Filename base
      character(len = *), intent(in) :: base_name
      integer, intent(in) :: dump_no
      !! The current iteration the adapt-re-load-balance loop
      integer, intent(in) :: adapt_iteration
      !! The total number of iterations to be performed in the
      !! adapt-re-load-balance loop
      integer, intent(in) :: max_adapt_iteration
      !! If present and .false., do not convert into a parallel file_name
      logical, optional, intent(in) :: add_parallel
      
      character(len = len_trim(base_name) + 1 + int2str_len(dump_no) + 1 + int2str_len(adapt_iteration) + parallel_filename_len("")) :: file_name
      
      file_name = trim(base_name) // "_" // int2str(dump_no)
      if(max_adapt_iteration > 1) file_name = trim(file_name) // "_" // int2str(adapt_iteration)
      if(.not. present_and_false(add_parallel) .and. isparallel()) file_name = parallel_filename(file_name)
      
    end function adapt_state_debug_file_name
        
  end subroutine write_adapt_state_debug_output
  
  subroutine adapt_state_module_check_options
  
    integer :: max_output, stat
    
    call get_option("/mesh_adaptivity/hr_adaptivity/debug/checkpoint/max_checkpoint_count", max_output, stat = stat)
    if(stat == SPUD_NO_ERROR) then
      if(max_output <= 0) then
        FLExit("Max adaptivity debug checkpoint count must be positive")
      end if
    end if
  
  end subroutine adapt_state_module_check_options

end module adapt_state_module
