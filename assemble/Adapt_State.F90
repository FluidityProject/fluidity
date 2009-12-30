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
  use global_parameters, only : OPTION_PATH_LEN
  use hadapt_extrude
  use hadapt_metric_based_extrude
  use halos
  use interpolation_manager
  use interpolation_module
  use mba_adapt_module
  use mba2d_integration
  use mba3d_integration
  use metric_assemble
  use parallel_tools
  use boundary_conditions_from_options
  use populate_state_module
  use project_metric_to_surface_module
  use reserve_state_module
  use sam_integration
  use tictoc  
  use timeloop_utilities
  use write_triangle
#ifdef HAVE_ZOLTAN
  use zoltan_integration
#endif
  
  implicit none
  
  private
  
  public :: adapt_mesh, adapt_state, adapt_state_first_timestep
  public :: insert_metric_for_interpolation, extract_and_remove_metric, sam_options
  public :: adapt_state_module_check_options
  
  interface adapt_state
    module procedure adapt_state_single, adapt_state_multiple
  end interface adapt_state
  
contains
  
  subroutine adapt_mesh(old_positions, metric, new_positions, node_ownership, force_preserve_regions)
    !!< A simple wrapper to select the appropriate adapt_mesh routine
    
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
    
    select case(old_positions%dim)
      case(1)
        call adapt_mesh_1d(old_positions, metric, new_positions, &
          & node_ownership = node_ownership, force_preserve_regions = force_preserve_regions)
      case(2)
        call adapt_mesh_mba2d(old_positions, metric, new_positions, &
          & force_preserve_regions=force_preserve_regions)
      case(3)
        if(have_option("/mesh_adaptivity/hr_adaptivity/adaptivity_library/libmba3d")) then
          call adapt_mesh_mba3d(old_positions, metric, new_positions, &
                             force_preserve_regions=force_preserve_regions)
        else
          call adapt_mesh_3d(old_positions, metric, new_positions, node_ownership = node_ownership, &
                             force_preserve_regions=force_preserve_regions)
        end if
      case default
        FLAbort("Mesh adaptivity requires a 1D, 2D or 3D mesh")
    end select
    
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
        call write_triangle_files(parallel_filename("first_timestep_adapted_mesh"), output_positions)
        call write_halos("first_timestep_adapted_mesh", output_positions%mesh)
      else
        call write_triangle_files("first_timestep_adapted_mesh", output_positions)
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

    ! Vertically structured adaptivity stuff
    type(vector_field) :: extruded_positions
    type(tensor_field) :: full_metric

    ewrite(1, *) "In adapt_state_internal"
    
    nullify(node_ownership)
    
    max_adapt_iteration = adapt_iterations()

    ! Don't need to strip the level 2 halo with Zoltan .. in fact, we don't want to
#ifndef HAVE_ZOLTAN
    if(isparallel()) then
      ! In parallel, strip off the level 2 halo (expected by libsam). The level
      ! 2 halo is restored on the final adapt iteration by libsam.
      call strip_level_2_halo(states, metric)
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
      ! Extract the mesh field to be adapted (takes a reference)
      old_positions = get_coordinate_field(states(1), old_linear_mesh)
      ewrite(2, *) "Mesh field to be adapted: " // trim(old_positions%name)
      
      call prepare_vertically_structured_adaptivity(states, metric, full_metric, extruded_positions, old_positions)
      
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
        call halo_update(interpolate_states, level = 1)
      end if 
      
      ! Before we start allocating any new objects we tag all references to
      ! current objects before the adapt so we can later on check they have all
      ! been deallocated
      call tag_references()
      
      ! Generate a new mesh field based on the current mesh field and the input
      ! metric
      call adapt_mesh(old_positions, metric, new_positions, node_ownership = node_ownership, &
        & force_preserve_regions=initialise_fields)
      
      ! Insert the new mesh field and linear mesh into all states
      call insert(states, new_positions%mesh, name = new_positions%mesh%name)
      call insert(states, new_positions, name = new_positions%name)
      
      call perform_vertically_inhomogenous_step(states, new_positions, full_metric, extruded_positions)

      ! We're done with old_positions, so we may deallocate it
      call deallocate(old_positions)
      
      ! Insert meshes from reserve states
      call restore_reserved_meshes(states)
      ! Next we recreate all derived meshes
      call insert_derived_meshes(states)
      ! Then reallocate all fields 
      call allocate_and_insert_fields(states)
      ! Insert fields from reserve states
      call restore_reserved_fields(states)
      
      if(i < max_adapt_iteration .or. isparallel()) then
        ! If there are remaining adapt iterations, or we will be calling
        ! sam_drive, insert the old metric into interpolate_states(1) and a
        ! new metric into states(1), for interpolation
        call insert_metric_for_interpolation(metric, new_positions%mesh, interpolate_states(1), states(1), metric_name = metric_name)
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
      
      ! The following is the same as the tail of populate_state:
      ! Add on the boundary conditions again
      call populate_boundary_conditions(states)
      ! Set their values
      call set_boundary_conditions_values(states)
      ! If strong bc or weak that overwrite then enforce the bc on the fields
      call set_dirichlet_consistent(states)
      ! Insert aliased fields in state
      call alias_fields(states)      

      if(isparallel()) then
#ifdef HAVE_ZOLTAN
        call zoltan_drive(states, i, metric=metric)
        if (i == max_adapt_iteration) then
          call deallocate(metric)
        end if
#else
        ! Re-load-balance using libsam
        call sam_drive(states, sam_options(i, max_adapt_iteration), metric = metric)
        if(i == max_adapt_iteration) then
          ! On the last adapt iteration the metric was interpolated
          ! only for sam_drive, hence it must be deallocated
          call deallocate(metric)
        end if
#endif
      end if
      
      if(no_reserved_meshes()) then
       ewrite(2, *) "Tagged references remaining:"
       call print_tagged_references(0)
      else
       ewrite(2, *) "There are reserved meshes, so skipping printing of references."
      end if
      
      call write_adapt_state_debug_output(states, i, max_adapt_iteration)      
      
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

  subroutine prepare_vertically_structured_adaptivity(states, metric, full_metric, extruded_positions, old_positions)
    type(state_type), dimension(:), intent(inout) :: states
    ! the metric will be collapsed, and the uncollapsed full_metric stored in full_metric
    type(tensor_field), intent(inout) :: metric, full_metric
    type(vector_field), intent(inout) :: extruded_positions
    ! old positions of the horizontal mesh
    type(vector_field), intent(in) :: old_positions

    integer, save:: adaptcnt=0
    logical :: vertically_structured_adaptivity
    logical :: vertically_inhomogenous_adaptivity

    type(scalar_field):: edge_lengths
      
    vertically_structured_adaptivity = have_option( &
     &  "/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity")
    vertically_inhomogenous_adaptivity = have_option( &
     &  "/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity/inhomogenous_vertical_resolution")


    if (vertically_structured_adaptivity) then
      ! project full mesh metric to horizontal surface mesh metric
      full_metric=metric
      call project_metric_to_surface(full_metric, old_positions, metric)
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
         ! we need the full_metric and its position field later on for vertical adaptivity
         ! this takes a reference so that it's prevented from the big deallocate in adapt_state
         extruded_positions = get_coordinate_field(states(1), full_metric%mesh)
      else
         ! otherwise we're done with it:
         call deallocate(full_metric)
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

    vertically_inhomogenous_adaptivity = have_option( &
     &  "/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity/inhomogenous_vertical_resolution")

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
       ! we can do away with the old metric now
       call deallocate(full_metric)
       call deallocate(extruded_positions)
       
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
  
  subroutine write_adapt_state_debug_output(states, adapt_iteration, max_adapt_iteration)
    !!< Diagnostic output for mesh adaptivity
  
    type(state_type), dimension(:), intent(in) :: states
    integer, intent(in) :: adapt_iteration
    integer, intent(in) :: max_adapt_iteration
    
    character(len = FIELD_NAME_LEN) :: file_name
    character(len = *), parameter :: base_path = "/mesh_adaptivity/hr_adaptivity/debug"
    integer :: max_output, stat
    type(mesh_type), pointer :: mesh
    type(vector_field) :: positions
    
    integer, save :: cp_no = 0, mesh_dump_no = 0, state_dump_no = 0
    
    if(.not. have_option(base_path)) return
    
    if(have_option(base_path // "/write_adapted_mesh")) then
      file_name = adapt_state_debug_file_name("adapted_mesh", mesh_dump_no, adapt_iteration, max_adapt_iteration)
      call find_mesh_to_adapt(states(1), mesh)
      positions = get_coordinate_field(states(1), mesh)
      call write_triangle_files(file_name, positions)
      if(isparallel()) then
        file_name = adapt_state_debug_file_name("adapted_mesh", mesh_dump_no, adapt_iteration, max_adapt_iteration, &
                                                add_parallel = .false.)  ! parallel extension is added by write_halos
        call write_halos(file_name, positions%mesh)
      end if
      call deallocate(positions)
      
      mesh_dump_no = mesh_dump_no + 1
    end if
    
    if(have_option(base_path // "/write_adapted_state")) then   
      file_name = adapt_state_debug_file_name("adapted_state", state_dump_no, adapt_iteration, max_adapt_iteration, add_parallel = .false.)   
      call vtk_write_state(file_name, state = states)
      
      state_dump_no = state_dump_no + 1
    end if
    
    if(adapt_iteration == max_adapt_iteration .and. have_option(base_path // "/checkpoint")) then
      call checkpoint_simulation(states, postfix = "adapt_checkpoint", cp_no = cp_no)
      
      cp_no = cp_no + 1
      
      call get_option(base_path // "/checkpoint/max_checkpoint_count", max_output, stat = stat)
      if(stat == SPUD_NO_ERROR) cp_no = modulo(cp_no, max_output)
    end if
    
  contains
  
    function adapt_state_debug_file_name(base_name, dump_no, adapt_iteration, max_adapt_iteration, add_parallel) result(file_name)
      character(len = *), intent(in) :: base_name
      integer, intent(in) :: dump_no
      integer, intent(in) :: adapt_iteration
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
