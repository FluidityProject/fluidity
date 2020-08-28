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

module sam_integration

  use fldebug
  use global_parameters, only : OPTION_PATH_LEN, FIELD_NAME_LEN
  use futils
  use reference_counting, only: print_tagged_references
  use quadrature
  use elements
  use spud
  use mpi_interfaces
  use parallel_tools
  use memory_diagnostics
  use data_structures
  use ieee_arithmetic
  use metric_tools
  use fields
  use state_module
  use field_options
  use halos
  use surfacelabels
  use node_boundary
  use boundary_conditions
  use tictoc
  use detector_data_types
  use boundary_conditions_from_options
  use reserve_state_module
  use pickers
  use detector_tools
  use diagnostic_variables
  use populate_state_module
  use surface_id_interleaving

  implicit none
  
  interface 
    !subroutine flstriph2(nnodes, nprivatenodes, nprocs, &
    !  & volumeenlist, nvolumeelems, nloc, &
    !  & surfaceenlist, surfaceids, nsurfaceelems, snloc, &
    !  & x, y, z, &
    !  & fields, nfields, fstride, &
    !  & metric, &
    !  & scatter, nscatter)
    !  implicit none
    !  integer, intent(inout) :: nnodes
    !  integer, intent(in) :: nprivatenodes
    !  integer, intent(in) :: nprocs
    !  integer, intent(inout) :: nvolumeelems
    !  integer, intent(in) :: nloc
    !  integer, dimension(nvolumeelems * nloc), intent(inout) :: volumeenlist
    !  integer, intent(inout) :: nsurfaceelems
    !  integer, intent(in) :: snloc
    !  integer, dimension(nsurfaceelems * snloc), intent(inout) :: surfaceenlist
    !  integer, dimension(nsurfaceelems), intent(inout) :: surfaceids
    !  real, dimension(nnodes), intent(inout) :: x
    !  real, dimension(nnodes), intent(inout) :: y
    !  real, dimension(nnodes), intent(inout) :: z
    !  integer, intent(inout) :: nfields
    !  integer, intent(inout) :: fstride
    !  real, dimension(nnodes * nfields * fstride), intent(inout) :: fields
    !  real, dimension(nnodes * 9), intent(inout) :: metric
    !  integer, intent(inout) :: nscatter
    !  integer, dimension(nscatter), intent(inout) :: scatter
    !end subroutine flstriph2
  
    subroutine sam_init_c(dim, nonods, totele, stotel, &
                       & gather, atosen, &
                       & scater, atorec, &
                       & ncolga, nscate, nprocs, &
                       & NDGLNO, nloc, &
                       & SNDGLN, SURFID, snloc, &
                       & X, Y, Z, &
                       & metric, FIELDS, NFIELDS, &
                       & options, MESTP1)
       implicit none
       integer, intent(in) :: dim
       integer, intent(in) :: nonods, totele, stotel
       integer, intent(in) :: ncolga
       integer, intent(in) :: nscate
       integer, intent(in) :: nprocs
       integer, intent(in) :: nloc, snloc
       integer, dimension(ncolga), intent(in) :: gather
       integer, dimension(nprocs + 1), intent(in) :: atosen
       integer, dimension(nscate), intent(in) :: scater
       integer, dimension(nprocs + 1), intent(in) :: atorec
       integer, intent(in), dimension(stotel * snloc) :: SNDGLN
       integer, intent(in), dimension(stotel) :: SURFID
       integer, intent(in), dimension(totele * nloc) :: NDGLNO
       real, dimension(nonods), intent(in) :: X, Y, Z
       real, dimension(dim ** 2 * nonods), intent(in) :: metric
       integer, intent(in) :: NFIELDS
       real, dimension(NFIELDS * NONODS), intent(in) :: FIELDS
       integer, dimension(10), intent(in) :: options
       real, intent(in) :: MESTP1
     end subroutine sam_init_c
   end interface 

   interface sam_migrate
     subroutine sam_migrate_c
     end subroutine sam_migrate_c
   end interface sam_migrate

   interface sam_add_field
     module procedure sam_add_field_scalar, sam_add_field_vector, sam_add_field_tensor
   end interface

   interface sam_query
     subroutine sam_query_c(NONODS, TOTELE, STOTEL, ncolga, nscate, pncolga, pnscate)
       implicit none
       integer, intent(out) :: NONODS, TOTELE, STOTEL
       integer, intent(out) :: ncolga
       integer, intent(out) :: nscate
       integer, intent(out) :: pncolga
       integer, intent(out) :: pnscate
     end subroutine sam_query_c
   end interface sam_query

   interface sam_cleanup
     subroutine sam_cleanup_c
     end subroutine sam_cleanup_c
   end interface sam_cleanup
   
   interface sam_export_mesh
     subroutine sam_export_mesh_c(nonods, totele, stotel, nloc, snloc, nodx, nody, nodz, enlist, senlist, surfid)
       implicit none
       integer, intent(in) :: nonods, totele, stotel, nloc, snloc
       real, dimension(nonods), intent(out) :: nodx
       real, dimension(nonods), intent(out) :: nody
       real, dimension(nonods), intent(out) :: nodz
       integer, dimension(totele * nloc), intent(out) :: enlist
       integer, dimension(stotel * snloc), intent(out) :: senlist
       integer, dimension(stotel), intent(out) :: surfid
     end subroutine sam_export_mesh_c
   end interface sam_export_mesh
   
   interface sam_export_halo
     subroutine sam_export_halo_c(colgat, atosen, scater, atorec, ncolga, nscate, nprocs, pnodes, nnodes)
       implicit none
       integer, intent(in) :: ncolga
       integer, intent(in) :: nscate
       integer, intent(in) :: nprocs
       integer, dimension(ncolga), intent(out) :: colgat
       integer, dimension(nprocs + 1), intent(out) :: atosen
       integer, dimension(nscate), intent(out) :: scater
       integer, dimension(nprocs + 1), intent(out) :: atorec
       integer, intent(out) :: pnodes
       integer, intent(out) :: nnodes
     end subroutine sam_export_halo_c
   end interface sam_export_halo

   interface sam_export_phalo
     subroutine sam_export_phalo_c(pcolgat, patosen, pscater, patorec, pncolga, pnscate, nprocs, ppnodes, pnnodes)
       implicit none
       integer, intent(in) :: pncolga
       integer, intent(in) :: pnscate
       integer, intent(in) :: nprocs
       integer, dimension(pncolga), intent(out) :: pcolgat
       integer, dimension(nprocs + 1), intent(out) :: patosen
       integer, dimension(pnscate), intent(out) :: pscater
       integer, dimension(nprocs + 1), intent(out) :: patorec
       integer, intent(out) :: ppnodes
       integer, intent(out) :: pnnodes
     end subroutine sam_export_phalo_c
   end interface sam_export_phalo
   
   interface sam_add_field
     subroutine sam_add_field_c(field_data, nnodes)
       implicit none
       integer, intent(in) :: nnodes
       real, dimension(nnodes), intent(in) :: field_data
     end subroutine sam_add_field_c
   end interface sam_add_field

   interface sam_pop_field
     subroutine sam_pop_field_c(field_data, nnodes)
       implicit none
       integer, intent(in) :: nnodes
       real, dimension(nnodes), intent(out) :: field_data
     end subroutine sam_pop_field_c
   end interface sam_pop_field
   
   interface sam_export_node_ownership
     subroutine sam_export_node_ownership_c(node_ownership, nnodes)
       implicit none
       integer, intent(in) :: nnodes
       integer, dimension(nnodes), intent(out) :: node_ownership
     end subroutine sam_export_node_ownership_c
   end interface sam_export_node_ownership
   
   private
   
   public :: sam_drive, strip_level_2_halo, sam_integration_check_options

   contains
   
     subroutine strip_level_2_halo(states, metric, external_mesh_name, initialise_fields)
       !!< Strip the level 2 halo from the supplied states and error metric.
       !!< Replaces flstriph2.
       
       type(state_type), dimension(:), intent(inout) :: states
       type(tensor_field), optional, intent(inout) :: metric
       character(len=FIELD_NAME_LEN), optional, intent(in) :: external_mesh_name
       logical, optional, intent(in) :: initialise_fields
       
       character(len = FIELD_NAME_LEN) :: linear_coordinate_field_name
       integer :: i, j, nlocal_dets, stat
       integer, dimension(:), allocatable :: renumber
       logical, dimension(:), allocatable :: keep
       type(halo_type), pointer :: level_1_halo, level_2_halo
       type(mesh_type) :: new_linear_mesh, old_linear_mesh
       type(mesh_type), pointer :: old_linear_mesh_ptr
       type(scalar_field), pointer :: new_s_field
       type(vector_field) :: new_positions, old_positions
       type(vector_field), pointer :: new_v_field
       type(state_type), dimension(:), allocatable :: interpolate_states
       type(tensor_field), pointer :: new_t_field
       type(tensor_field) :: new_metric
       
       ewrite(1, *) "In strip_level_2_halo"
       
       ! Find the external mesh. Must be linear and continuous.
       old_linear_mesh_ptr => get_external_mesh(states, external_mesh_name=external_mesh_name)

       old_linear_mesh = old_linear_mesh_ptr
       old_linear_mesh_ptr => null()
       call incref(old_linear_mesh)
       ewrite(2, *) "External mesh: " // trim(old_linear_mesh%name)
       if(trim(old_linear_mesh%name) == "CoordinateMesh") then
          linear_coordinate_field_name = "Coordinate"
       else
          linear_coordinate_field_name = trim(old_linear_mesh%name) // "Coordinate"
       end if
       ewrite(2, *) "Mesh field: " // trim(linear_coordinate_field_name)
       
       ! Extract the mesh field
       old_positions = extract_vector_field(states(1), linear_coordinate_field_name)
       call incref(old_positions)
       assert(old_positions%mesh == old_linear_mesh)
       
       call initialise_boundcount(old_linear_mesh, old_positions)
       
       ! Use select_fields_to_interpolate to reference all non-recoverable
       ! information in interpolate_states
       allocate(interpolate_states(size(states)))
       do i = 1, size(states)
         call select_fields_to_interpolate(states(i), interpolate_states(i), first_time_step=initialise_fields)
         ! If the old mesh field is referenced in interpolate_states(i), remove
         ! it (it will be dealt with seperately)
         call remove_vector_field(interpolate_states(i), old_positions%name, stat)
       end do
       
       ! Extract the level 1 and level 2 halos
       assert(associated(old_positions%mesh%halos))
       assert(size(old_positions%mesh%halos) >= 2)
       level_1_halo => old_positions%mesh%halos(1)
       call incref(level_1_halo)
       level_2_halo => old_positions%mesh%halos(2)
       call incref(level_2_halo)       

       ! Deallocate all recoverable information
       do i = 1, size(states)
         call deallocate(states(i))
       end do
       
       ! Find the nodes to keep
       allocate(keep(node_count(old_positions)))
       call find_nodes_to_keep(keep, level_1_halo, level_2_halo)
       call deallocate(level_2_halo)
       
       ! Generate the renumbering map
       allocate(renumber(size(keep)))
       call create_renumbering_map(renumber, keep)
       
       ewrite(2, *) "Stripping level 2 halo from the external mesh"
       call generate_stripped_linear_mesh(old_linear_mesh, new_linear_mesh, level_1_halo, keep, renumber)
       call insert(states, new_linear_mesh, new_linear_mesh%name)
       
       ewrite(2, *) "Stripping level 2 halo from the mesh field"
       call allocate(new_positions, mesh_dim(new_linear_mesh), new_linear_mesh, old_positions%name)
       call generate_stripped_vector_field(old_positions, old_linear_mesh, new_positions, new_linear_mesh, keep)
       call insert(states, new_positions, new_positions%name)
       call deallocate(old_positions)
       
       nlocal_dets = default_stat%detector_list%length
       call allsum(nlocal_dets)
       if(nlocal_dets > 0) call halo_transfer_detectors(old_linear_mesh, new_positions)
       call deallocate(new_positions)
       
       ! Insert meshes from reserve states
       call restore_reserved_meshes(states)
       ! Next we recreate all derived meshes
       call insert_derived_meshes(states)
       ! Then reallocate all fields 
       call allocate_and_insert_fields(states)
       ! Insert fields from reserve states
       call restore_reserved_fields(states)
       
       ! Strip the level 2 halo from all fields in states
       do i = 1, size(interpolate_states)
         do j = 1, scalar_field_count(interpolate_states(i))
           assert(associated(interpolate_states(i)%scalar_fields(j)%ptr))
           ewrite(2, *) "Stripping level 2 halo from field " // trim(interpolate_states(i)%scalar_fields(j)%ptr%name) // " in state " // trim(states(i)%name)
           new_s_field => extract_scalar_field(states(i), interpolate_states(i)%scalar_fields(j)%ptr%name)
           call generate_stripped_scalar_field(interpolate_states(i)%scalar_fields(j)%ptr, old_linear_mesh, new_s_field, new_linear_mesh, keep)
         end do
         
         do j = 1, vector_field_count(interpolate_states(i))
           assert(associated(interpolate_states(i)%vector_fields(j)%ptr))
           if(trim(interpolate_states(i)%vector_fields(j)%ptr%name) == trim(linear_coordinate_field_name)) cycle
           ewrite(2, *) "Stripping level 2 halo from field " // trim(interpolate_states(i)%vector_fields(j)%ptr%name) // " in state " // trim(states(i)%name)
           new_v_field => extract_vector_field(states(i), interpolate_states(i)%vector_fields(j)%ptr%name)
           call generate_stripped_vector_field(interpolate_states(i)%vector_fields(j)%ptr, old_linear_mesh, new_v_field, new_linear_mesh, keep)
         end do
         
         do j = 1, tensor_field_count(interpolate_states(i))
           assert(associated(interpolate_states(i)%tensor_fields(j)%ptr))
           ewrite(2, *) "Stripping level 2 halo from field " // trim(interpolate_states(i)%tensor_fields(j)%ptr%name) // " in state " // trim(states(i)%name)
           new_t_field => extract_tensor_field(states(i), interpolate_states(i)%tensor_fields(j)%ptr%name)
           call generate_stripped_tensor_field(interpolate_states(i)%tensor_fields(j)%ptr, old_linear_mesh, new_t_field, new_linear_mesh, keep)
         end do
       end do
       
       do i = 1, size(interpolate_states)
         call deallocate(interpolate_states(i))
       end do
       deallocate(interpolate_states)
       
       ewrite(2, *) "Renumbering level 1 halo"
       call renumber_halo(level_1_halo, renumber)
#ifdef DDEBUG
       if(isparallel()) then
         assert(halo_verifies(level_1_halo, extract_vector_field(states(1), linear_coordinate_field_name)))
       end if
#endif
       call deallocate(level_1_halo)
       
       if(present(metric)) then
         ewrite(2, *) "Stripping level 2 halo from metric " // trim(metric%name)
         call allocate(new_metric, new_linear_mesh, metric%name)
         call generate_stripped_tensor_field(metric, old_linear_mesh, new_metric, new_linear_mesh, keep)
         call deallocate(metric)
         call allocate(metric, new_linear_mesh, metric%name)
         call set(metric, new_metric)
         call deallocate(new_metric)
#ifdef DDEBUG
         call check_metric(metric)
#endif
       end if
       
       call deallocate(old_linear_mesh)
       call deallocate(new_linear_mesh)

       deallocate(keep)
       deallocate(renumber)
       
       ! The following is the same as the tail of populate_state:    
       ! Prescribed fields are recalculated
       call set_prescribed_field_values(states, exclude_interpolated=.true.)
       ! Add on the boundary conditions again
       call populate_boundary_conditions(states)
       ! Set their values
       call set_boundary_conditions_values(states)
       ! if strong bc or weak that overwrite then enforce the bc on the fields
       call set_dirichlet_consistent(states)
       ! Insert aliased fields in state
       call alias_fields(states)
       
       if(no_reserved_meshes()) then
          ewrite(2, *) "Tagged references remaining:"
          call print_tagged_references(0)
       else
          ewrite(2, *) "There are reserved meshes, so skipping printing of references."
       end if
              
       ewrite(1, *) "Exiting strip_level_2_halo"
       
     end subroutine strip_level_2_halo

     subroutine find_nodes_to_keep(keep, level_1_halo, level_2_halo)
       !!< Set the keep array, deciding whether a node should be kept
       
       logical, dimension(:), intent(out) :: keep
       type(halo_type), intent(in) :: level_1_halo
       type(halo_type), intent(in) :: level_2_halo
       
       integer :: i, j
       
       assert(size(keep) >= max_halo_node(level_1_halo))
       assert(size(keep) >= max_halo_node(level_2_halo))
       assert(halo_proc_count(level_1_halo) == halo_proc_count(level_2_halo))
       
       keep = .true.
       do i = 1, halo_proc_count(level_2_halo)
         do j = 1, halo_receive_count(level_2_halo, i)
           keep(halo_receive(level_2_halo, i, j)) = .false.
         end do
       end do
       do i = 1, halo_proc_count(level_1_halo)
         do j = 1, halo_receive_count(level_1_halo, i)
           keep(halo_receive(level_1_halo, i, j)) = .true.
         end do
       end do
       
     end subroutine find_nodes_to_keep
     
     subroutine create_renumbering_map(renumber, keep)
       !!< Generate the map used for node renumbering when stripping the halo.
       !!< renumber is negative if the node is to be stripped, and forms
       !!< a consecutive list elsewhere.
       
       integer, dimension(:), intent(out) :: renumber
       logical, dimension(size(renumber)), intent(in) :: keep
       
       integer :: i, index
       
       index = 0
       renumber = -1
       do i = 1, size(keep)
         if(keep(i)) then
           index = index + 1
           renumber(i) = index
         end if
       end do
       
     end subroutine create_renumbering_map
     
     subroutine generate_stripped_linear_mesh(input_linear_mesh, output_linear_mesh, level_1_halo, keep, renumber)
       !!< Generate a new mesh based on the input mesh, with nodes stripped
       !!< as specified by keep and renumber. output_linear_mesh is allocated
       !!< by this routine.
       
       type(mesh_type), intent(in) :: input_linear_mesh
       type(mesh_type), intent(out) :: output_linear_mesh
       type(halo_type), intent(in) :: level_1_halo
       logical, dimension(:), intent(in) :: keep
       integer, dimension(size(keep)), intent(in) :: renumber
       
       integer :: i, j, new_boundary_ids_size, new_coplanar_ids_size, new_nelms, new_ndglno_size, new_nnodes, new_region_ids_size, new_sndgln_size, nowned_nodes
       integer, dimension(:), allocatable :: sloc, new_boundary_ids, new_coplanar_ids, new_ndglno, new_region_ids, new_sndgln
       integer, dimension(:), pointer :: loc
       logical :: keep_element
       type(element_type) :: output_shape
       type(quadrature_type) :: output_quad
       
       assert(size(keep) == node_count(input_linear_mesh))
       assert(trailing_receives_consistent(level_1_halo))

       nowned_nodes = halo_nowned_nodes(level_1_halo)

       ! Count number of nodes in the stripped mesh
       if(size(renumber) > 0) then
         new_nnodes = max(maxval(renumber), 0)
       else
         new_nnodes = 0
       end if
       
       ! Strip volume elements
       new_nelms = 0
       allocate(new_ndglno(size(input_linear_mesh%ndglno)))
       new_ndglno_size = 0
       allocate(new_region_ids(ele_count(input_linear_mesh)))
       new_region_ids_size = 0
       volume_element_loop: do i = 1, ele_count(input_linear_mesh)
         assert(ele_loc(input_linear_mesh, i) == ele_loc(input_linear_mesh, 1))
         loc => ele_nodes(input_linear_mesh, i)
         keep_element = .true.
         element_loc_node_loop: do j = 1, size(loc)
           assert(loc(j) >= lbound(keep, 1) .and. loc(j) <= ubound(keep, 1))
           if(.not. keep(loc(j))) then
             keep_element = .false.
             exit element_loc_node_loop
           end if
         end do element_loc_node_loop
         ! This really shouldn't be necessary, but some elements have all
         ! nodes in the level 1 halo and no owned nodes!
         if(keep_element .and. .not. any(loc <= nowned_nodes)) then
           keep_element = .false.
           !ewrite(2, *) "Warning: Element found that has all nodes in the level 1 halo and no owned nodes"
         end if
         keep_element = keep_element .and. any(loc <= nowned_nodes)
         if(keep_element) then
           new_nelms = new_nelms + 1
           new_ndglno_size = new_ndglno_size + size(loc)
           ! Renumber nodes in the new volume element list
           add_volume_element_loop: do j = 1, size(loc)
             assert(loc(j) >= lbound(renumber, 1) .and. loc(j) <= ubound(renumber, 1))
             assert(renumber(loc(j)) > 0)
             new_ndglno(new_ndglno_size - size(loc) + j) = renumber(loc(j))
           end do add_volume_element_loop
           if(associated(input_linear_mesh%region_ids)) then
             new_region_ids_size = new_region_ids_size + 1
             new_region_ids(new_region_ids_size) = input_linear_mesh%region_ids(i)
           end if
         end if
       end do volume_element_loop
       
       ! Strip surface elements
       if(surface_element_count(input_linear_mesh) > 0) then
         allocate(new_sndgln(surface_element_count(input_linear_mesh) * face_loc(input_linear_mesh, 1)))
       else
         allocate(new_sndgln(0))
       end if
       new_sndgln_size = 0
       allocate(new_boundary_ids(surface_element_count(input_linear_mesh)))
       new_boundary_ids_size = 0
       allocate(new_coplanar_ids(surface_element_count(input_linear_mesh)))
       new_coplanar_ids_size = 0
       surface_element_loop: do i = 1, surface_element_count(input_linear_mesh)
         assert(face_loc(input_linear_mesh, i) == face_loc(input_linear_mesh, 1))
         loc => ele_nodes(input_linear_mesh, face_ele(input_linear_mesh, i))
         keep_element = .true.
         face_loc_node_loop: do j = 1, size(loc)
           assert(loc(j) >= lbound(keep, 1) .and. loc(j) <= ubound(keep, 1))
           if(.not. keep(loc(j))) then
             keep_element = .false.
             exit face_loc_node_loop
           end if
         end do face_loc_node_loop
         ! This really shouldn't be necessary, but some elements have all
         ! nodes in the level 1 halo and no owned nodes!
         if(keep_element .and. .not. any(loc <= nowned_nodes)) then
           keep_element = .false.
           !ewrite(2, *) "Warning: Surface element found attached to element that has all nodes in the level 1 halo and no owned nodes"
         end if
         if(keep_element) then
           allocate(sloc(face_loc(input_linear_mesh, i)))
           sloc = face_global_nodes(input_linear_mesh, i)
           new_sndgln_size = new_sndgln_size + size(sloc)
           ! Renumber nodes in the new surface element list
           add_surface_element_loop: do j = 1, size(sloc)
             assert(sloc(j) >= lbound(renumber, 1) .and. sloc(j) <= ubound(renumber, 1))
             assert(renumber(sloc(j)) > 0)
             new_sndgln(new_sndgln_size - size(sloc) + j) = renumber(sloc(j))
           end do add_surface_element_loop
           if(associated(input_linear_mesh%faces%boundary_ids)) then
             new_boundary_ids_size = new_boundary_ids_size + 1
             new_boundary_ids(new_boundary_ids_size) = input_linear_mesh%faces%boundary_ids(i)
           end if
           if(associated(input_linear_mesh%faces%coplanar_ids)) then
             new_coplanar_ids_size = new_coplanar_ids_size + 1
             new_coplanar_ids(new_coplanar_ids_size) = input_linear_mesh%faces%coplanar_ids(i)
           end if
           deallocate(sloc)
         end if
       end do surface_element_loop
       
       ! Construct the new mesh
       output_quad = make_quadrature(ele_loc(input_linear_mesh, 1), mesh_dim(input_linear_mesh), degree = input_linear_mesh%shape%quadrature%degree)
       output_shape = make_element_shape(ele_loc(input_linear_mesh, 1), mesh_dim(input_linear_mesh), input_linear_mesh%shape%degree, output_quad)
       call allocate(output_linear_mesh, new_nnodes, new_nelms, output_shape, name = input_linear_mesh%name)
       call deallocate(output_quad)
       call deallocate(output_shape)
       output_linear_mesh%ndglno = new_ndglno(1:new_ndglno_size)
       if(isparallel()) then
         if(associated(input_linear_mesh%faces%boundary_ids)) then
           call add_faces(output_linear_mesh, sndgln = new_sndgln(1:new_sndgln_size), boundary_ids = new_boundary_ids(1:new_boundary_ids_size))
         else
           call add_faces(output_linear_mesh, sndgln = new_sndgln(1:new_sndgln_size))
         end if
       else
         if(associated(input_linear_mesh%faces%boundary_ids)) then
           call add_faces(output_linear_mesh, sndgln = new_sndgln(1:new_sndgln_size), boundary_ids = new_boundary_ids(1:new_boundary_ids_size))
         else
           call add_faces(output_linear_mesh, sndgln = new_sndgln(1:new_sndgln_size))
         end if
       end if
       if(associated(input_linear_mesh%faces%coplanar_ids)) then
         allocate(output_linear_mesh%faces%coplanar_ids(new_coplanar_ids_size))
         output_linear_mesh%faces%coplanar_ids = new_coplanar_ids(1:new_coplanar_ids_size)
       end if
       if(associated(input_linear_mesh%region_ids)) then
         allocate(output_linear_mesh%region_ids(new_region_ids_size))
         output_linear_mesh%region_ids = new_region_ids(1:new_region_ids_size)
       end if
       output_linear_mesh%option_path = input_linear_mesh%option_path
       
       allocate(output_linear_mesh%halos(1))
       output_linear_mesh%halos(1)=level_1_halo
       call incref(level_1_halo)
       if(.not. serial_storage_halo(level_1_halo)) then  ! Cannot derive halos in serial
         allocate(output_linear_mesh%element_halos(1))
         call derive_element_halo_from_node_halo(output_linear_mesh, &
           & ordering_scheme = HALO_ORDER_TRAILING_RECEIVES)
       end if

       deallocate(new_ndglno)
       deallocate(new_boundary_ids)
       deallocate(new_coplanar_ids)
       deallocate(new_sndgln)
       deallocate(new_region_ids)
       
     end subroutine generate_stripped_linear_mesh
    
     subroutine generate_stripped_scalar_field(input_field, input_linear_mesh, output_field, output_linear_mesh, keep)
       !!< Generate a new field based on the input field, with nodes stripped
       !!< as specified by keep.
       
       type(scalar_field), intent(in) :: input_field
       type(mesh_type), intent(inout) :: input_linear_mesh
       type(scalar_field), intent(inout) :: output_field
       type(mesh_type), intent(inout) :: output_linear_mesh
       logical, dimension(:), intent(in) :: keep
      
       integer :: i, index
       type(scalar_field) :: input_linear_field, output_linear_field
      
       assert(mesh_dim(input_field) == mesh_dim(output_field))
       assert(size(keep) == node_count(input_linear_mesh))
       
       call allocate(input_linear_field, input_linear_mesh, input_field%name)
       call remap_field(input_field, input_linear_field)
       
       call allocate(output_linear_field, output_linear_mesh, input_field%name)
       
       index = 0
       do i = 1, node_count(input_linear_field)
         if(keep(i)) then
           index = index + 1
           assert(index <= node_count(output_linear_field))
           call set(output_linear_field, index, node_val(input_linear_field, i))
         end if
       end do
       
       call remap_field(output_linear_field, output_field)
       
       call deallocate(input_linear_field)
       call deallocate(output_linear_field)
       
     end subroutine generate_stripped_scalar_field
    
     subroutine generate_stripped_vector_field(input_field, input_linear_mesh, output_field, output_linear_mesh, keep)
       !!< Generate a new field based on the input field, with nodes stripped
       !!< as specified by keep.
        
       type(vector_field), intent(in) :: input_field
       type(mesh_type), intent(in) :: input_linear_mesh
       type(vector_field), target, intent(inout) :: output_field
       type(mesh_type), intent(in) :: output_linear_mesh
       logical, dimension(:), intent(in) :: keep
      
       integer :: i, index
       type(vector_field) :: input_linear_field, output_linear_field
      
       assert(mesh_dim(input_field) == mesh_dim(output_field))
       assert(size(keep) == node_count(input_linear_mesh))
       
       call allocate(input_linear_field, mesh_dim(input_linear_mesh), input_linear_mesh, input_field%name)
       call remap_field(input_field, input_linear_field)
       
       call allocate(output_linear_field, mesh_dim(output_linear_mesh), output_linear_mesh, input_field%name)
       
       index = 0
       do i = 1, node_count(input_linear_field)
         if(keep(i)) then
           index = index + 1
           assert(index <= node_count(output_linear_field))
           call set(output_linear_field, index, node_val(input_linear_field, i))
         end if
       end do
       
       call remap_field(output_linear_field, output_field)
       
       call deallocate(input_linear_field)
       call deallocate(output_linear_field)
       
     end subroutine generate_stripped_vector_field    
    
     subroutine generate_stripped_tensor_field(input_field, input_linear_mesh, output_field, output_linear_mesh, keep)
       !!< Generate a new field based on the input field, with nodes stripped
       !!< as specified by keep. output_field is allocated by this routine.
       
       type(tensor_field), intent(in) :: input_field
       type(mesh_type), intent(in) :: input_linear_mesh
       type(tensor_field), intent(inout) :: output_field
       type(mesh_type), intent(in) :: output_linear_mesh
       logical, dimension(:), intent(in) :: keep
      
       integer :: i, index
       type(tensor_field) :: input_linear_field, output_linear_field
      
       assert(mesh_dim(input_field) == mesh_dim(output_field))
       assert(size(keep) == node_count(input_linear_mesh))
       
       call allocate(input_linear_field, input_linear_mesh, input_field%name)
       call remap_field(input_field, input_linear_field)
       
       call allocate(output_linear_field, output_linear_mesh, input_field%name)
       
       index = 0
       do i = 1, node_count(input_linear_field)
         if(keep(i)) then
           index = index + 1
           assert(index <= node_count(output_linear_field))
           call set(output_linear_field, index, node_val(input_linear_field, i))
         end if
       end do
       
       call remap_field(output_linear_field, output_field)
       
       call deallocate(input_linear_field)
       call deallocate(output_linear_field)
       
     end subroutine generate_stripped_tensor_field  
    
     subroutine renumber_halo(halo, renumber)
       !!< Renumber the supplied halo according to the specified renumbering.
      
       type(halo_type), intent(inout) :: halo
       integer, dimension(:), intent(in) :: renumber
      
       integer :: i, j, receive

       do i = 1, halo_proc_count(halo)
         do j = 1, halo_receive_count(halo, i)
           receive = halo_receive(halo, i, j)
           assert(receive >= lbound(renumber, 1) .and. receive <= ubound(renumber, 1))
           assert(renumber(receive) > 0)
           call set_halo_receive(halo, i, j, renumber(receive))
         end do
       end do
       
       assert(halo_valid_for_communication(halo))
      
     end subroutine renumber_halo

     subroutine sam_drive(states, options, metric, external_mesh_name, initialise_fields)
       type(state_type), dimension(:), intent(inout) :: states
       integer, dimension(10), intent(in) :: options
       type(tensor_field), optional, intent(inout) :: metric
       character(len=FIELD_NAME_LEN), intent(in), optional :: external_mesh_name
       ! if present and true: don't bother redistributing fields that can be reinitialised
       logical, intent(in), optional :: initialise_fields
       
       type(state_type), dimension(size(states)) :: interpolate_states

       integer :: field, state

       type(scalar_field) :: linear_s
       type(vector_field) :: linear_v
       type(tensor_field) :: linear_t

       character(len=FIELD_NAME_LEN), dimension(:, :), allocatable, target :: namelist_s, namelist_v, namelist_t
       type(scalar_field), pointer :: field_s
       type(vector_field), pointer :: field_v
       type(tensor_field), pointer :: field_t

       integer, dimension(size(states)) :: scount, vcount, tcount

       type(element_type) :: linear_shape

       integer :: max_coplanar_id, nlocal_dets
       integer, dimension(:), allocatable :: boundary_ids, coplanar_ids
       
       integer :: NNODP, NONODS, TOTELE, STOTEL, ncolga, nscate, pncolga, pnscate
       type(mesh_type) :: old_linear_mesh
       type(mesh_type), pointer :: linear_mesh

       integer, dimension(:), allocatable :: ATOSEN, ATOREC
       integer, dimension(:), allocatable :: COLGAT, SCATER
       type(vector_field), target :: new_positions
       type(vector_field), pointer :: old_positions
       integer :: dim, snloc, nloc
       integer, dimension(:), allocatable :: senlist, surface_ids
       character(len=FIELD_NAME_LEN) :: linear_mesh_name, linear_coordinate_field_name, metric_name
       character(len=OPTION_PATH_LEN) :: linear_mesh_option_path
       integer :: component, component_i, component_j

       real, dimension(:,:), allocatable :: xyz
       real, dimension(:), allocatable :: value
   
       integer :: stat


       ewrite(1, *) "In sam_drive"
       call tic(TICTOC_ID_DATA_REMAP)

       ! Step 1. Initialise sam.
       if(present(metric)) then
         call sam_init(states, options, max_coplanar_id, metric = metric, external_mesh_name=external_mesh_name)
       else
         call sam_init(states, options, max_coplanar_id, external_mesh_name=external_mesh_name)
       end if

       ! Step 2. Supply sam with all the fields it needs to migrate.
       old_linear_mesh = get_external_mesh(states, external_mesh_name=external_mesh_name)

       nlocal_dets = default_stat%detector_list%length
       call allsum(nlocal_dets)
       if(nlocal_dets > 0) then
         ! Detector communication required. Take a reference to the old mesh.
         ! If no detector communication is required, avoid taking a reference to
         ! save memory.
         call incref(old_linear_mesh)
       end if
       
       linear_mesh_name = old_linear_mesh%name
       linear_mesh_option_path = old_linear_mesh%option_path
       if(trim(linear_mesh_name) == "CoordinateMesh") then
          linear_coordinate_field_name="Coordinate"
       else
          linear_coordinate_field_name=trim(linear_mesh_name)//"Coordinate"
       end if

       old_positions => extract_vector_field(states(1), trim(linear_coordinate_field_name))
       dim = old_positions%dim
       linear_shape = ele_shape(old_linear_mesh, 1)

       nloc = old_linear_mesh%shape%loc
       snloc = old_linear_mesh%faces%surface_mesh%shape%loc
       call incref(linear_shape)

       call allocate(linear_s, old_linear_mesh, "LinearScalarField")
       call allocate(linear_v, dim, old_linear_mesh, "LinearVectorField")
       call allocate(linear_t, old_linear_mesh, "LinearTensorField")

       ! Get the sizes so I can allocate the right amount of headers
       ! Multiple states is a REAL pain

       ! Select fields to interpolate
       do state=1,size(states)
         call select_fields_to_interpolate(states(state), interpolate_states(state), &
           no_positions=.true., first_time_step=initialise_fields)
       end do

       ! Record the headers of the interpolated fields
       scount = 0
       vcount = 0
       tcount = 0
       do state=1,size(states)
         scount(state) = scalar_field_count(interpolate_states(state))
         vcount(state) = vector_field_count(interpolate_states(state))
         tcount(state) = tensor_field_count(interpolate_states(state))
       end do

       allocate(namelist_s(size(states), maxval(scount)))
       allocate(namelist_v(size(states), maxval(vcount)))
       allocate(namelist_t(size(states), maxval(tcount)))

       do state=1,size(states)
         do field=1,scalar_field_count(interpolate_states(state))
           field_s => extract_scalar_field(interpolate_states(state), field)
           namelist_s(state, field) = field_s%name
         end do

         do field=1,vector_field_count(interpolate_states(state))
           field_v => extract_vector_field(interpolate_states(state), field)
           namelist_v(state, field) = field_v%name
         end do

         do field=1,tensor_field_count(interpolate_states(state))
           field_t => extract_tensor_field(interpolate_states(state), field)
           namelist_t(state, field) = field_t%name
         end do
       end do

       do state=1,size(states)
         do field=1,scount(state)
           field_s => extract_scalar_field(interpolate_states(state), trim(namelist_s(state, field)))
           call remap_field(field_s, linear_s,stat)
           call check_sam_linear_remap_validity(stat, trim(namelist_s(state, field)))
           call sam_add_field(linear_s)
           call remove_scalar_field(states(state), trim(namelist_s(state, field)))
           call remove_scalar_field(interpolate_states(state), trim(namelist_s(state, field)))
         end do

         do field=1,vcount(state)
           field_v => extract_vector_field(interpolate_states(state), trim(namelist_v(state, field)))
           call remap_field(field_v, linear_v,stat)
           call check_sam_linear_remap_validity(stat, trim(namelist_v(state, field)))
           call sam_add_field(linear_v)
           call remove_vector_field(states(state), trim(namelist_v(state, field)))
           call remove_vector_field(interpolate_states(state), trim(namelist_v(state, field)))
         end do

         do field=1,tcount(state)
           field_t => extract_tensor_field(interpolate_states(state), trim(namelist_t(state, field)))
           call remap_field(field_t, linear_t,stat)
           call check_sam_linear_remap_validity(stat, trim(namelist_t(state, field)))
           call sam_add_field(linear_t)
           call remove_tensor_field(states(state), trim(namelist_t(state, field)))
           call remove_tensor_field(interpolate_states(state), trim(namelist_t(state, field)))
         end do
       end do
       
       if(present(metric)) then
         ! Add the metric
         call remap_field(metric, linear_t)
         call sam_add_field(linear_t)
       end if

       call deallocate(linear_s)
       call deallocate(linear_v)
       call deallocate(linear_t)
       
       if(present(metric)) metric_name = metric%name

       ! Step 3. Deallocate.
       
       ! Deallocate the states
       do state=1,size(states)
         call deallocate(states(state))
         call deallocate(interpolate_states(state))
       end do
       if(present(metric)) then
         ! Deallocate the metric
         call deallocate(metric)
       end if
       
       ! Step 4. Migrate.
       ewrite(1, *) "Calling sam_migrate from sam_drive"
       call tic(TICTOC_ID_DATA_MIGRATION)
       call sam_migrate
       call toc(TICTOC_ID_DATA_MIGRATION) 
       ewrite(1, *) "Exited sam_migrate"

       ! Step 5. Now, we need to reconstruct.
       
       ! Query the statistics of the new mesh.
       ewrite(1, *) "Calling sam_query from sam_drive"
       call sam_query(nonods, totele, stotel, ncolga, nscate, pncolga, pnscate)
       ewrite(1, *) "Exited sam_query"
       
       !!! sanity check

       if (nonods==0) then
          FLExit("Libsam has produced an empty partition for your problem. Please consider reconfiguring with Zoltan")
       end if

       ! Export mesh data from sam
       allocate(linear_mesh)
       call allocate(linear_mesh, nonods, totele, linear_shape, linear_mesh_name)
       call deallocate(linear_shape)
       call allocate(new_positions, dim, linear_mesh, linear_coordinate_field_name)         
       call deallocate(linear_mesh)
       deallocate(linear_mesh)
       linear_mesh => new_positions%mesh            
       allocate(surface_ids(stotel))
       allocate(senlist(stotel * snloc))
       ewrite(1, *) "Calling sam_export_mesh from sam_drive"
       allocate(xyz(1:nonods, 1:3))
       call sam_export_mesh(nonods, totele, stotel, nloc, snloc, &
                          & xyz(:,1), xyz(:,2), xyz(:,3), &
                          & linear_mesh%ndglno, senlist, surface_ids)
       new_positions%val=transpose(xyz(:,1:new_positions%dim))
       deallocate(xyz)
       ewrite(1, *) "Exited sam_export_mesh"
       linear_mesh%option_path = linear_mesh_option_path
         
       if(nlocal_dets > 0) then
         ! Communicate the local detectors
         call sam_transfer_detectors(old_linear_mesh, new_positions)
         call deallocate(old_linear_mesh)
       end if
       
       ! Add the surface mesh data
       allocate(boundary_ids(stotel))
       allocate(coplanar_ids(stotel))
       call deinterleave_surface_ids(surface_ids, max_coplanar_id, boundary_ids, coplanar_ids)
       call add_faces(linear_mesh, sndgln=senlist, boundary_ids=boundary_ids)
       deallocate(boundary_ids)
       allocate(linear_mesh%faces%coplanar_ids(stotel))
       linear_mesh%faces%coplanar_ids = coplanar_ids
       deallocate(coplanar_ids)
       deallocate(surface_ids)
       deallocate(senlist)

       ! Check that the level 2 halo is around
       if(pncolga >= 0) then
         allocate(linear_mesh%halos(2))
       else
         allocate(linear_mesh%halos(1))
       end if

       ! Export the level 1 halo
       allocate(colgat(ncolga))
       allocate(scater(nscate))
       allocate(atosen(getnprocs() + 1))
       allocate(atorec(getnprocs() + 1))
       ewrite(1, *) "Calling sam_export_halo from sam_drive"
       call sam_export_halo(colgat, atosen, scater, atorec, ncolga, nscate, getnprocs(), nnodp, nonods)
       ewrite(1, *) "Exited sam_export_halo"
       assert(nonods == node_count(linear_mesh))
       ! Form a halo from the primitive data structures
       call form_halo_from_raw_data(linear_mesh%halos(1), getnprocs(), colgat, atosen, scater, atorec, nowned_nodes = nnodp)
       ! Deallocate the primitive data structures
       deallocate(colgat)
       deallocate(atosen)
       deallocate(scater)
       deallocate(atorec)
#ifdef DDEBUG
       if(isparallel()) then
         ! Check the new halo
         assert(trailing_receives_consistent(linear_mesh%halos(1)))
         assert(halo_valid_for_communication(linear_mesh%halos(1)))
         assert(halo_verifies(linear_mesh%halos(1), new_positions))
       end if
#endif
      
       if(pncolga >= 0) then
         ! In "mixed formulation" export the level 2 halo
         assert(pnscate >= 0)
         allocate(colgat(pncolga))
         allocate(scater(pnscate))
         allocate(atosen(getnprocs() + 1))
         allocate(atorec(getnprocs() + 1))
         ewrite(1, *) "Calling sam_export_phalo from sam_drive"
         call sam_export_phalo(colgat, atosen, scater, atorec, pncolga, pnscate, getnprocs(), nnodp, nonods)
         ewrite(1, *) "Exited sam_export_phalo"
         assert(nnodp == halo_nowned_nodes(linear_mesh%halos(1)))
         assert(nonods == node_count(linear_mesh))
         ! Form a halo from the primitive data structures
         call form_halo_from_raw_data(linear_mesh%halos(2), getnprocs(), colgat, atosen, scater, atorec, nowned_nodes = nnodp)
         ! Deallocate the primitive data structures
         deallocate(colgat)
         deallocate(atosen)
         deallocate(scater)
         deallocate(atorec)
         ! Check the new halo
         assert(trailing_receives_consistent(linear_mesh%halos(2)))
         assert(halo_valid_for_communication(linear_mesh%halos(2)))
         assert(halo_verifies(linear_mesh%halos(2), new_positions))

         ! Derive the elements halo
         allocate(linear_mesh%element_halos(2))
         call derive_element_halo_from_node_halo(linear_mesh, &
           & ordering_scheme = HALO_ORDER_TRAILING_RECEIVES)  
       else
         if(.not. serial_storage_halo(linear_mesh%halos(1))) then  ! Cannot derive halos in serial
           allocate(linear_mesh%element_halos(1))
           call derive_element_halo_from_node_halo(linear_mesh, &
             & ordering_scheme = HALO_ORDER_TRAILING_RECEIVES)
         else
           allocate(linear_mesh%element_halos(0)) 
         end if
       end if

       ! Insert the positions and linear mesh into all states
       call insert(states, linear_mesh, trim(linear_mesh_name))
       call insert(states, new_positions, trim(linear_coordinate_field_name))
       
       ! Insert meshes from reserve state
       call restore_reserved_meshes(states)
       ! Next we recreate all derived meshes
       call insert_derived_meshes(states)
       ! Then reallocate all fields 
       call allocate_and_insert_fields(states, dont_allocate_prognostic_value_spaces = .true.)
       
       ! Now extract all fields in the reverse order we put them in
       call allocate(linear_s, linear_mesh, "LinearScalarField")
       call allocate(linear_v, dim, linear_mesh, "LinearVectorField")
       call allocate(linear_t, linear_mesh, "LinearTensorField")
       
       if(present(metric)) then
         do component_i=dim,1,-1
           do component_j=dim,1,-1
             call sam_pop_field(linear_t%val(component_i, component_j, :), node_count(linear_mesh))
           end do
         end do
         call allocate(metric, linear_mesh, name = metric_name)
         call remap_field(linear_t, metric)
#ifdef DDEBUG
         call check_metric(metric)
#endif
       end if
       
       allocate(value(1:node_count(linear_mesh)))
       ! Extract field data from sam
       do state=size(states),1,-1
         do field=tcount(state),1,-1
           do component_i=dim,1,-1
             do component_j=dim,1,-1
               call sam_pop_field(linear_t%val(component_i, component_j, :), node_count(linear_mesh))
             end do
           end do

           field_t => extract_tensor_field(states(state), trim(namelist_t(state, field)))
           deallocate(field_t%val)
           allocate(field_t%val(dim, dim, node_count(field_t%mesh)))
#ifdef HAVE_MEMORY_STATS
           call register_allocation("tensor_field", "real", node_count(field_t%mesh)*mesh_dim(field_t%mesh)**2, name=trim(field_t%name))
#endif
           field_t%field_type = FIELD_TYPE_NORMAL

           call remap_field(linear_t, field_t)
         end do

         do field=vcount(state),1,-1
           do component=dim,1,-1
             call sam_pop_field(value, node_count(linear_mesh))
             call set_all(linear_v, component, value)
           end do

           field_v => extract_vector_field(states(state), trim(namelist_v(state, field)))
           deallocate(field_v%val)
           allocate(field_v%val(dim,node_count(field_v%mesh)))
#ifdef HAVE_MEMORY_STATS
           call register_allocation("vector_field", "real", node_count(field_v%mesh)*mesh_dim(field_v%mesh), name=trim(field_v%name))
#endif
           field_v%field_type = FIELD_TYPE_NORMAL

           call remap_field(linear_v, field_v)
         end do

         do field=scount(state),1,-1
           call sam_pop_field(linear_s%val, node_count(linear_mesh))

           field_s => extract_scalar_field(states(state), trim(namelist_s(state, field)))
           deallocate(field_s%val)
           allocate(field_s%val(node_count(field_s%mesh)))
#ifdef HAVE_MEMORY_STATS
           call register_allocation("scalar_field", "real", size(field_s%val), name=trim(field_s%name))
#endif
           field_s%field_type = FIELD_TYPE_NORMAL

           call remap_field(linear_s, field_s)
         end do
       end do

       call deallocate(linear_s)
       call deallocate(linear_v)
       call deallocate(linear_t)
       
       deallocate(namelist_s)
       deallocate(namelist_v)
       deallocate(namelist_t)
       
       ! We're done with the new positions now so we may drop our reference
       call deallocate(new_positions)

       ! Make sure all fields have their value space allocated, even those that
       ! aren't interpolated
       call allocate_remaining_fields(states)

       ! The following is the same as the tail of populate_state:    
       ! Prescribed fields are recalculated
       call set_prescribed_field_values(states, exclude_interpolated=.true.)
       ! Add on the boundary conditions again
       call populate_boundary_conditions(states)
       ! Set their values
       call set_boundary_conditions_values(states)
       ! if strong bc or weak that overwrite then enforce the bc on the fields
       call set_dirichlet_consistent(states)
       ! Insert aliased fields in state
       call alias_fields(states)
       
       ! Step 6. Cleanup
       
       ewrite(1, *) "Calling sam_cleanup from sam_drive"
       call sam_cleanup
       ewrite(1, *) "Exited sam_cleanup"
       
       call ewrite_load_imbalance(2, "Owned nodes:", nnodp)
       call ewrite_load_imbalance(2, "Total nodes:", nonods)
       call ewrite_load_imbalance(2, "Total elements:", totele)

       call toc(TICTOC_ID_DATA_REMAP)       
       ewrite(1, *) "Exiting sam_drive"
     end subroutine sam_drive

     subroutine allocate_remaining_fields(states)
       !!< Allocate all fields that have not had their value spaces allocated,
       !!< but which are non-constant

       type(state_type), dimension(:), intent(inout) :: states

       integer :: i, j, k
       type(scalar_field), pointer :: s_field
       type(tensor_field), pointer :: t_field
       type(vector_field), pointer :: v_field

       do i = 1, size(states)
         do j = 1, scalar_field_count(states(i))
           s_field => extract_scalar_field(states(i), j)
           if(s_field%field_type == FIELD_TYPE_DEFERRED) then
             deallocate(s_field%val)
             allocate(s_field%val(node_count(s_field%mesh)))
#ifdef HAVE_MEMORY_STATS
             call register_allocation("scalar_field", "real", size(s_field%val), name=trim(s_field%name))
#endif
             s_field%field_type = FIELD_TYPE_NORMAL
             call zero(s_field)
           end if
         end do
         do j = 1, vector_field_count(states(i))
           v_field => extract_vector_field(states(i), j)
           if(v_field%field_type == FIELD_TYPE_DEFERRED) then
             deallocate(v_field%val)
             allocate(v_field%val(mesh_dim(v_field%mesh),node_count(v_field%mesh)))
#ifdef HAVE_MEMORY_STATS
             call register_allocation("vector_field", "real", node_count(v_field%mesh)*mesh_dim(v_field%mesh), name=trim(v_field%name))
#endif
             v_field%field_type = FIELD_TYPE_NORMAL
             call zero(v_field)
           end if
         end do
         do j = 1, tensor_field_count(states(i))
           t_field => extract_tensor_field(states(i), j)
           if(t_field%field_type == FIELD_TYPE_DEFERRED) then
             deallocate(t_field%val)
             allocate(t_field%val(mesh_dim(t_field%mesh), mesh_dim(t_field%mesh), node_count(t_field%mesh)))
#ifdef HAVE_MEMORY_STATS
             call register_allocation("tensor_field", "real", node_count(t_field%mesh)*mesh_dim(t_field%mesh)**2, name=trim(t_field%name))
#endif
             t_field%field_type = FIELD_TYPE_NORMAL
             call zero(t_field)
           end if
         end do
       end do

     end subroutine allocate_remaining_fields

     subroutine sam_init(states, options, max_coplanar_id, metric, external_mesh_name)
       !!< Initialise sam, with the external mesh in the supplied states.

       type(state_type), dimension(:), intent(in) :: states
       integer, intent(out) :: max_coplanar_id
       type(tensor_field), optional, intent(in) :: metric
       character(len=FIELD_NAME_LEN), optional, intent(in) :: external_mesh_name

       ! sam_init_c variables
       integer :: nonods, totele, stotel
       integer, dimension(:), allocatable :: scater, atorec, gather, atosen
       integer :: nscate
       integer, dimension(:), pointer :: ndglno
       integer, dimension(:), allocatable :: surfid, sndgln
       integer :: nloc, snloc
       real, dimension(:), pointer :: x, y, z
       real, dimension(:), allocatable :: metric_handle
       integer :: nfields
       real, dimension(:), pointer :: fields
       real, dimension(1), target :: dummy
       integer, dimension(10), intent(in) :: options
       real :: mestp1
       real, dimension(:,:), allocatable :: xyz

       integer :: dim, i, j, nprocs
       type(halo_type) :: halo
       type(mesh_type), pointer :: mesh
       type(vector_field), pointer :: positions
       
       ewrite(1, *) "In sam_init"

       ! Extract the external mesh and mesh field
       mesh => get_external_mesh(states, external_mesh_name=external_mesh_name)
       positions => extract_vector_field(states(1), "Coordinate")
       dim = mesh_dim(mesh)
       nonods = node_count(mesh)
       totele = ele_count(mesh)
       stotel = unique_surface_element_count(mesh)
       nloc = mesh%shape%loc
       snloc = mesh%faces%surface_mesh%shape%loc
       ndglno => mesh%ndglno
       allocate(sndgln(stotel * snloc))
       call getsndgln(mesh, sndgln)
       allocate(surfid(unique_surface_element_count(mesh)))
       call interleave_surface_ids(mesh, surfid, max_coplanar_id)
       
       ! Extract the level 1 halo
       nprocs = getnprocs()       
       if(halo_count(mesh) > 0) then
         halo = mesh%halos(1)
         assert(trailing_receives_consistent(halo))
         assert(halo_valid_for_communication(halo))
         
         ! Copy the halo data into primitive data structures
         allocate(gather(halo_all_sends_count(halo)))
         allocate(atosen(halo_proc_count(halo) + 1))
         nscate = halo_all_receives_count(halo)
         allocate(scater(nscate))
         allocate(atorec(halo_proc_count(halo) + 1))
         call extract_raw_halo_data(halo, gather, atosen, scater, atorec)
       else
         if(isparallel()) then
           ewrite(-1, *) "Warning: sam_init called in parallel with no level one halo"
         end if
         allocate(gather(0))
         allocate(atosen(nprocs + 1))
         nscate = 0
         allocate(scater(nscate))
         allocate(atorec(nprocs + 1))
         atosen = 0
         atorec = 0
       end if
       
       ! Form the metric
       allocate(metric_handle(dim * dim * nonods))
       if(present(metric)) then
         metric_handle = reshape(metric%val, (/nonods * dim ** 2/))
       else
         ! Gerard: Do we really want to allocate 9 * nonods
         ! just for the identity? For a man worried about memory
         ! this is incredibly wasteful.
         metric_handle = 0.0
         forall(i = 0:nonods - 1, j = 0:dim - 1)
           metric_handle(i * dim ** 2 + 1 + j * dim) = 1.0
         end forall
       end if

       ! The field data is taken care of later
       nfields = 0
       dummy = 0
       fields => dummy
       
       allocate(xyz(1:nonods,1:3))
       xyz(:,1:positions%dim)=transpose(positions%val)
       xyz(:,positions%dim+1:)=0.0

       call get_option('/mesh_adaptivity/hr_adaptivity/functional_tolerance', mestp1, default = 0.0)

       ewrite(1, *) "Calling sam_init_c from sam_init"
       call sam_init_c(dim, nonods, totele, stotel, &
                          & gather, atosen, &
                          & scater(1:nscate), atorec, &
                          & size(gather), nscate, nprocs, &
                          & ndglno(1:totele * nloc), nloc, &
                          & sndgln(1:stotel * snloc), surfid(1:stotel), snloc, &
                          & xyz(:,1), xyz(:,2), xyz(:,3), &
                          & metric_handle(1:nonods * dim ** 2), fields, nfields, &
                          & options, mestp1)
       deallocate(xyz)
       ewrite(1, *) "Exited sam_init_c"
       
       deallocate(sndgln)
       deallocate(surfid)
       deallocate(gather)
       deallocate(atosen)
       deallocate(scater)
       deallocate(atorec)
       deallocate(metric_handle)
       
       ewrite(1, *) "Exiting sam_init"

     end subroutine sam_init

     subroutine sam_add_field_scalar(field)
       type(scalar_field), intent(in) :: field
       
       call sam_add_field(field%val, node_count(field))
       
     end subroutine sam_add_field_scalar

     subroutine sam_add_field_vector(field)
       type(vector_field), intent(in) :: field
       
       real, dimension(:), allocatable:: value
       integer :: i
       
       allocate( value(1:node_count(field)) )
       do i=1,field%dim
         value=field%val(i,:)
         call sam_add_field(value, node_count(field))
       end do
       
     end subroutine sam_add_field_vector

     subroutine sam_add_field_tensor(field)
       type(tensor_field), intent(in) :: field
       
       integer :: i, j
       
       do i=1,mesh_dim(field%mesh)
         do j=1,mesh_dim(field%mesh)
           call sam_add_field(field%val(i, j, :), node_count(field))
         end do
       end do
       
     end subroutine sam_add_field_tensor
     
  subroutine sam_transfer_detectors(old_mesh, new_positions)
    type(mesh_type), intent(in) :: old_mesh
    type(vector_field), intent(inout) :: new_positions
    
    integer :: nnodes
    integer, dimension(:), allocatable :: node_ownership
    
    nnodes = node_count(old_mesh)
    allocate(node_ownership(nnodes))
    call sam_export_node_ownership(node_ownership, nnodes)
    node_ownership = node_ownership + 1
    
    call transfer_detectors(old_mesh, new_positions, node_ownership)
    
    deallocate(node_ownership)
  
  end subroutine sam_transfer_detectors
  
  subroutine halo_transfer_detectors(old_mesh, new_positions)
    type(mesh_type), intent(in) :: old_mesh
    type(vector_field), intent(inout) :: new_positions
  
    integer :: nhalos, nnodes
    integer, dimension(:), allocatable :: node_ownership
    type(halo_type), pointer :: halo
    
    nhalos = halo_count(old_mesh)
    if(nhalos == 0) return
    halo => old_mesh%halos(nhalos)
    
    nnodes = node_count(old_mesh)
    allocate(node_ownership(nnodes))
    call get_node_owners(halo, node_ownership)
    
    call transfer_detectors(old_mesh, new_positions, node_ownership)
    
    deallocate(node_ownership)
    
  end subroutine halo_transfer_detectors
     
  subroutine transfer_detectors(old_mesh, new_positions, node_ownership)
    type(mesh_type), intent(in) :: old_mesh
    type(vector_field), intent(inout) :: new_positions
    integer, dimension(node_count(old_mesh)) :: node_ownership
    
    integer :: communicator, i, j, nhalos,  nprocs, &
      & owner, procno
    type(detector_type), pointer :: next_node, node
    type(halo_type), pointer :: halo
    
    integer, parameter :: idata_size = 2
    integer :: rdata_size
    
    integer, dimension(:), allocatable :: nsends, data_index
    type(integer_vector), dimension(:), allocatable :: isend_data
    type(real_vector), dimension(:), allocatable :: rsend_data
    
    integer, dimension(:), allocatable :: nreceives
    type(integer_vector), dimension(:), allocatable :: ireceive_data
    type(real_vector), dimension(:), allocatable :: rreceive_data
    
    integer :: ierr, tag
    integer, dimension(:), allocatable :: requests, statuses
    
    nhalos = halo_count(old_mesh)
    if(nhalos == 0) return
    halo => old_mesh%halos(nhalos)
    communicator = halo_communicator(halo)
    procno = getprocno(communicator = communicator)
    nprocs = halo_proc_count(halo)
    
    rdata_size = new_positions%dim
    
    allocate(nsends(nprocs))
    nsends = 0
    node => default_stat%detector_list%first
    do while(associated(node))
      if(node%element > 0) then
        owner = minval(node_ownership(ele_nodes(old_mesh, node%element)))
        if(owner /= procno) then
          nsends(owner) = nsends(owner) + 1
        end if
      end if      
    
      node => node%next
    end do
    
    allocate(isend_data(nprocs))
    allocate(rsend_data(nprocs))    
    do i = 1, nprocs
      allocate(isend_data(i)%ptr(nsends(i) * idata_size))
      allocate(rsend_data(i)%ptr(nsends(i) * rdata_size))
    end do
    
    allocate(data_index(nprocs))
    data_index = 0
    node => default_stat%detector_list%first
    do while(associated(node))
      next_node => node%next
    
      if(node%element > 0) then
        owner = minval(node_ownership(ele_nodes(old_mesh, node%element)))
        if(owner /= procno) then
          ! Pack this node for sending
          
          ! Integer data
          isend_data(owner)%ptr(data_index(owner) * idata_size + 1) = node%id_number
          ! Real data
          rsend_data(owner)%ptr(data_index(owner) * rdata_size + 1:data_index(owner) * rdata_size + new_positions%dim) = node%position
                    
          data_index(owner) = data_index(owner) + 1
                    
          ! Delete this node from the detector list
          call delete(node, default_stat%detector_list)
        end if
      end if      
    
      node => next_node
    end do
    deallocate(data_index)
    
    ewrite(2, *) "Detectors to be sent: ", sum(nsends)
    
    allocate(nreceives(nprocs))
    nreceives = invert_comms_sizes(nsends, communicator = communicator)
    
    ewrite(2, *) "Detectors to be received: ", sum(nreceives)
    
    allocate(ireceive_data(nprocs))
    allocate(rreceive_data(nprocs))
    do i = 1, nprocs
      allocate(ireceive_data(i)%ptr(nreceives(i) * idata_size))
      allocate(rreceive_data(i)%ptr(nreceives(i) * rdata_size))
    end do
    
    ! Set up non-blocking communications
    allocate(requests(nprocs * 4))
    requests = MPI_REQUEST_NULL
    tag = next_mpi_tag()

    do i = 1, nprocs
      ! Non-blocking sends
      if(nsends(i) > 0) then
        call mpi_isend(isend_data(i)%ptr, nsends(i) * idata_size, getpinteger(), i - 1, tag, communicator, requests(i), ierr)
        assert(ierr == MPI_SUCCESS)
        call mpi_isend(rsend_data(i)%ptr, nsends(i) * rdata_size, getpreal(), i - 1, tag, communicator, requests(i + nprocs), ierr)
        assert(ierr == MPI_SUCCESS)
      end if
      
      ! Non-blocking receives
      if(nreceives(i) > 0) then
        call mpi_irecv(ireceive_data(i)%ptr, nreceives(i) * idata_size, getpinteger(), i - 1, tag, communicator, requests(i + 2 * nprocs), ierr)
        assert(ierr == MPI_SUCCESS)
        call mpi_irecv(rreceive_data(i)%ptr, nreceives(i) * rdata_size, getpreal(), i - 1, tag, communicator, requests(i + 3 * nprocs), ierr)
        assert(ierr == MPI_SUCCESS)
      end if
    end do   
        
    ! Wait for all non-blocking communications to complete
    allocate(statuses(MPI_STATUS_SIZE * size(requests)))
    call mpi_waitall(size(requests), requests, statuses, ierr)
    assert(ierr == MPI_SUCCESS)
    deallocate(statuses)
    deallocate(requests)
    
    deallocate(nsends)
    do i = 1, nprocs
      deallocate(isend_data(i)%ptr)
      deallocate(rsend_data(i)%ptr)
    end do
    deallocate(isend_data)
    deallocate(rsend_data)
    
    do i = 1, nprocs
      do j = 1, nreceives(i)
        ! Unpack the node
        allocate(node)
        
        ! Integer data
        node%id_number = ireceive_data(i)%ptr((j - 1) * idata_size + 1)
                
        ! Real data
        allocate(node%position(new_positions%dim))
        node%position = rreceive_data(i)%ptr((j - 1) * rdata_size + 1:(j - 1) * rdata_size + new_positions%dim)
        
        ! Recoverable data, not communicated
        node%name = default_stat%detector_list%detector_names(node%id_number)
        allocate(node%local_coords(new_positions%dim + 1))
        
        call insert(node, default_stat%detector_list)
      end do
      
      deallocate(ireceive_data(i)%ptr)
      deallocate(rreceive_data(i)%ptr)
    end do      
    
    deallocate(nreceives)
    deallocate(ireceive_data)
    deallocate(rreceive_data)
    
    ! Update the detector element ownership data
    call search_for_detectors(default_stat%detector_list, new_positions)
  
  end subroutine transfer_detectors

  function load_imbalance(count)
    !!< Calculates the load imbalance metric:
    !!<   (Max in a domain - Mean in a domain)
    !!<   (----------------------------------)
    !!<   (     Mean nodes in a domain       )

    integer, intent(in) :: count

    real :: load_imbalance
    
#ifdef HAVE_MPI
    integer :: i, max_count, min_count, ierr
    integer, dimension(:), allocatable :: node_counts
    real :: mean
    
    allocate(node_counts(getnprocs()))
    
    call mpi_gather(count, 1, getpinteger(), node_counts, 1, getpinteger(), 0, MPI_COMM_FEMTOOLS, ierr)
    assert(ierr == MPI_SUCCESS)
    
    if(getprocno() == 1) then
      max_count = node_counts(1)
      min_count = node_counts(1)
      mean = node_counts(1)
      do i = 2, size(node_counts)
        max_count = max(max_count, node_counts(i))
        min_count = min(min_count, node_counts(i))
        mean = mean + node_counts(i)
      end do
      mean = mean / size(node_counts)
      load_imbalance = (max_count - mean) / mean
    end if
    
    call mpi_bcast(load_imbalance, 1, getpreal(), 0, MPI_COMM_FEMTOOLS, ierr)
    assert(ierr == MPI_SUCCESS)
    
    deallocate(node_counts)
#else    
    load_imbalance = 0.0
#endif
    
  end function load_imbalance
  
  subroutine ewrite_load_imbalance(debug_level, prefix, count)
    !!< ewrite the load imbalance at the supplied debug level, based upon the
    !!< supplied count for this process, adding a prefix.
    
    integer, intent(in) :: debug_level
    character(len = *), intent(in) :: prefix
    integer, intent(in) :: count
    
    integer :: max_count, min_count
    real :: mean, imbalance
    
    if(debug_level > current_debug_level) return
    
    max_count = count
    call allmax(max_count)
    
    min_count = count
    call allmin(min_count)
    
    mean = count
    call allmean(mean)

    imbalance = load_imbalance(count)         
    
    if(getprocno() == 1) then
      ewrite(debug_level, *) prefix
      ewrite(debug_level, *) "Mean = ", mean
      ewrite(debug_level, *) "Min. = ", min_count
      ewrite(debug_level, *) "Max. = ", max_count
      ewrite(debug_level, *) "Imbalance = ", imbalance
    end if
  
  end subroutine ewrite_load_imbalance

  subroutine check_sam_linear_remap_validity(stat,name)
    integer, intent(in) :: stat
    character(len = * ) :: name


    !! short circuit for the trivial case
    if (stat==0) return
    
    ewrite(0,*) "For field ", trim(name)

    select case(stat)
    case (REMAP_ERR_DISCONTINUOUS_CONTINUOUS)
       ewrite(0,*) "Unable to redistribute discontinuous field."
    case(REMAP_ERR_HIGHER_LOWER_CONTINUOUS)
       ewrite(0,*) "Unable to redistribute higher order field."
    case(REMAP_ERR_UNPERIODIC_PERIODIC)
       ewrite(0,*) "Unable to redistribute field on periodic mesh."
    case(REMAP_ERR_BUBBLE_LAGRANGE)
       ewrite(0,*) "Unable to redistribute field on finite element mesh with bubble function."
    case default
       ewrite(0,*) "Failure to remap. This probably means a discretisation type that is not handled by sam"
    end select
    FLExit("This discretisation is not supported in parallel with libsam. Please consider reconfiguring with Zoltan")

  end subroutine check_sam_linear_remap_validity
  
  subroutine sam_integration_check_options

    integer :: i
    character (len=OPTION_PATH_LEN) :: continuity_var

    !!< Check libsam integration related options
   
    if(.not. isparallel()) then
      ! Nothing to check
      return
    end if       
    
#ifndef HAVE_ZOLTAN
    ewrite(2, *) "Checking libsam integration related options"

    if( have_option("/flredecomp") ) then
       FLExit("Specification of flredecomp parameters in the options tree is not supported with libsam. Please remove or reconfigure with Zoltan")
    end if
    
    if(have_option("/mesh_adaptivity/hr_adaptivity/preserve_mesh_regions")) then
      FLExit("Preserving of mesh regions through adapts is not supported in parallel with libsam. Please reconfigure with Zoltan")
    end if

    if(option_count('/geometry/mesh/from_mesh/extrude')>0) then
       FLExit("Mesh extrusion is not supported in parallel with libsam. Please reconfigure with Zoltan")
    end if

    if(option_count('/geometry/mesh/from_mesh/mesh_shape')+option_count('/geometry/mesh/from_mesh/mesh_continuity')>0 .and. &
       have_option('/mesh_adaptivity/hr_adaptivity')) then
      ! there are meshes that change the mesh_shape or continuity
      ! from this we assume: 1) there are non P1 meshes, 2) there are fields on these meshes that need to be distributed by SAM
      ! neither of those are necessarily always true, but it will be the case in 98% of the cases

      ! sam does not handle such fields

      ! we allow this for non-adaptive cases - sam is then only going to be invoked during an flredecomp - if the fields
      ! on the non-P1 can be represcribed they don't actually need to be handled by sam. This means flredecomp on a non-checkpoint
      ! .flml will usually still work.
      ewrite(0,*) "It appears you have non P1 meshes (mesh that are not linear and continuous) and are using adaptivity."
      ewrite(0,*) "For this to work you need to reconfigure with Zoltan."
      FLExit("Non supported discretisation for sam with parallel adaptivity.")
    end if

    ewrite(2, *) "Finished checking libsam integration related options"
#endif

  end subroutine sam_integration_check_options

end module sam_integration
