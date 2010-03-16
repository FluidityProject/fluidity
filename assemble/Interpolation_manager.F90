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

module interpolation_manager

  use interpolation_module
  use conservative_interpolation_module
  use dg_interpolation_module
  use state_module
  use fields
  use spud
  use global_parameters, only : OPTION_PATH_LEN
  use supermesh_construction
  use field_options
  use intersection_finder_module
  use linked_lists
  use boundary_conditions_from_options
  use tictoc
  use geostrophic_pressure
  use quicksort
  
  implicit none
  
  private
  
  public :: interpolate, interpolation_manager_check_options
  
  interface collect_fields_to_interpolate
    module procedure collect_fields_to_interpolate_single_state, collect_fields_to_interpolate_multiple_states
  end interface

contains

  subroutine interpolate(states_old, states_new, map)
    !! OK! We need to figure out what algorithm to use when and where.
    type(state_type), dimension(:), intent(inout) :: states_old, states_new
    !! Map from new nodes to old elements
    integer, dimension(:), optional, intent(in) :: map

    ! The fields organised by mesh:
    type(state_type), dimension(:), allocatable :: meshes_old, meshes_new
    ! Within each mesh, the fields organised by algorithm.
    type(state_type), dimension(:), allocatable :: alg_old, alg_new
    
    integer :: state, state_cnt, mesh
    integer :: mesh_cnt
    character(len=FIELD_NAME_LEN) :: mesh_name
    
    integer :: unique_mesh_index
    integer, dimension(:), allocatable :: mesh_ele_counts, mesh_indices, unique_mesh_indices
    logical, dimension(:), allocatable :: duplicate_mesh

    type(scalar_field), pointer :: field_s
    type(vector_field), pointer :: field_v, field_v_2
    type(tensor_field), pointer :: field_t
    integer :: field

    character(len=255), dimension(4), parameter :: algorithms = (/&
                       & "consistent_interpolation ", &
                       & "interpolation_galerkin   ", &
                       & "geostrophic_interpolation", &
                       & "grandy_interpolation     " /)
    integer :: alg_cnt, alg

    type(mesh_type), pointer :: old_mesh, new_mesh
    type(vector_field) :: old_pos, new_pos
    type(ilist), dimension(:), allocatable :: map_BA
    integer :: ele

    logical :: all_consistent_interpolation, all_linear_meshes
    type(element_type), pointer :: field_shape

    integer :: no_fields
    
    integer :: stat

    ewrite(1, *) "In interpolate"
    
    assert(size(states_old) == size(states_new))
    state_cnt = size(states_old)

    ! First thing: check the usual case. If all field request consistent
    ! interpolation, and all fields are on linear meshes, the just pass over to
    ! linear interpolation

    all_consistent_interpolation = .true.
    all_linear_meshes = .true.

    consistent_linear_state_loop: do state=1,state_cnt
      do field=1,scalar_field_count(states_old(state))
        field_s => extract_scalar_field(states_old(state), field)
        ! If the field has no option path, assume consistent interpolation
        if(len_trim(field_s%option_path) /= 0) then
          assert(have_option(trim(field_s%option_path)//"/prognostic").or.have_option(trim(field_s%option_path)//"/prescribed"))

          if (.not. have_option(trim(complete_field_path(field_s%option_path, stat)) // &
                   & "/consistent_interpolation")) then
            all_consistent_interpolation = .false.
          end if
          
          field_shape => ele_shape(field_s, 1)
          if(field_shape%degree /= 1) then
            all_linear_meshes = .false.
          end if
        end if
      end do

      do field=1,vector_field_count(states_old(state))
        field_v => extract_vector_field(states_old(state), field)
        ! If the field has no option path, assume consistent interpolation
        if(len_trim(field_v%option_path) /= 0) then
          if (len_trim(field_v%name) >= len_trim("Coordinate")) then
            if (field_v%name(len_trim(field_v%name)-10+1:) == trim("Coordinate")) then
              cycle
            end if
          end if
          assert(have_option(trim(field_v%option_path)//"/prognostic").or.have_option(trim(field_v%option_path)//"/prescribed"))

          if (.not. have_option(trim(complete_field_path(field_v%option_path, stat)) // &
                   & "/consistent_interpolation")) then
            all_consistent_interpolation = .false.
          end if
          
          field_shape => ele_shape(field_v, 1)
          if(field_shape%degree /= 1) then
            all_linear_meshes = .false.
          end if
        end if
      end do

      do field=1,tensor_field_count(states_old(state))
        field_t => extract_tensor_field(states_old(state), field)
        ! If the field has no option path, assume consistent interpolation
        if(len_trim(field_t%option_path) /= 0) then
          assert(have_option(trim(field_t%option_path)//"/prognostic").or.have_option(trim(field_t%option_path)//"/prescribed"))

          if (.not. have_option(trim(complete_field_path(field_t%option_path, stat)) // &
                   & "/consistent_interpolation")) then
            all_consistent_interpolation = .false.
          end if
          
          field_shape => ele_shape(field_t, 1)
          if(field_shape%degree /= 1) then
            all_linear_meshes = .false.
          end if
        end if
      end do

      if(.not. all_consistent_interpolation .and. .not. all_linear_meshes) exit consistent_linear_state_loop
    end do consistent_linear_state_loop
    
    if(all_consistent_interpolation .and. all_linear_meshes) then
      ewrite(2, *) "All fields are on linear meshes and use consistent interpolation"

      call tictoc_clear(TICTOC_ID_INTERPOLATION)
      call tic(TICTOC_ID_INTERPOLATION)

      ! Assuming here that "map" is for the linear mesh
      if(present(map)) then
        call linear_interpolate_states(states_old, states_new, map = map)
      else
        call linear_interpolate_states(states_old, states_new)
      end if

      call toc(TICTOC_ID_INTERPOLATION)
      call tictoc_report(2, TICTOC_ID_INTERPOLATION)

      ewrite(1, *) "Exiting interpolate"
      
      return
    end if

    ewrite(2, *) "Not all fields are on linear meshes and use consistent interpolation"
    ewrite(2, *) "Gathering fields for more general interpolation"

    ! OK! So we have some work to do.
    ! First, let's organise the fields according to what mesh
    ! they are on.

    ! The following is a temporary hack - it identifies the unique meshes in
    ! states_old(1)
    
    ! 1. Sort the meshes by their element count to reduce the number of
    ! mesh comparisons
    mesh_cnt = mesh_count(states_old(1))
    allocate(mesh_ele_counts(mesh_cnt))
    do mesh = 1, mesh_cnt
      old_mesh => extract_mesh(states_old(1), mesh)
      mesh_ele_counts(mesh) = ele_count(old_mesh)
    end do
    allocate(mesh_indices(size(mesh_ele_counts)))
    call qsort(mesh_ele_counts, mesh_indices)
    call apply_permutation(mesh_ele_counts, mesh_indices)
        
    ! 2. Count up the unique meshes
    allocate(unique_mesh_indices(mesh_cnt))
    allocate(duplicate_mesh(mesh_cnt))
    duplicate_mesh = .false.
    unique_mesh_indices(mesh_indices(1)) = 1
    mesh_cnt = 1  
    do mesh = 2, size(mesh_ele_counts)
      if(mesh_ele_counts(mesh) == mesh_ele_counts(mesh - 1)) then
        old_mesh => extract_mesh(states_old(1), mesh_indices(mesh))
        if(old_mesh == extract_mesh(states_old(1), mesh_indices(mesh - 1))) then
          ewrite(2, *) "Duplicate mesh: " // trim(old_mesh%name)
          unique_mesh_indices(mesh_indices(mesh)) = unique_mesh_indices(mesh_indices(mesh - 1))
          duplicate_mesh(mesh_indices(mesh)) = .true.
          cycle
        end if
      end if

      mesh_cnt = mesh_cnt + 1
      unique_mesh_indices(mesh_indices(mesh)) = mesh_cnt
    end do
    deallocate(mesh_indices)
    deallocate(mesh_ele_counts)
    ewrite(2, *) "Unique meshes: ", mesh_cnt
    ewrite(2, *) "Duplicate meshes: ", mesh_count(states_old(1)) - mesh_cnt
        
    allocate(meshes_old(mesh_cnt))
    allocate(meshes_new(mesh_cnt))
    allocate(alg_old(mesh_cnt))
    allocate(alg_new(mesh_cnt))
    alg_cnt = size(algorithms)
    
    do mesh = 1, mesh_count(states_old(1))
      old_mesh => extract_mesh(states_old(1), mesh)
      mesh_name = old_mesh%name
      unique_mesh_index = unique_mesh_indices(mesh)
      assert(unique_mesh_index >= 1 .and. unique_mesh_index <= mesh_cnt)
      if(.not. duplicate_mesh(mesh)) then
        call insert(meshes_old(unique_mesh_index), old_mesh, "Mesh")
        call insert(meshes_new(unique_mesh_index), extract_mesh(states_new(1), mesh_name), "Mesh")
      end if
      
      do state=1,state_cnt
        do field=1,scalar_field_count(states_old(state))
          field_s => extract_scalar_field(states_old(state), field)
          if (trim(field_s%mesh%name) == trim(mesh_name)) then
            ! we need to append the state name here to make this safe for
            ! multi-material_phase/state... let's just hope you aren't going
            ! to try to pull this out of state by its name!
            call insert(meshes_old(unique_mesh_index), field_s, trim(states_new(state)%name)//"::"//trim(field_s%name))
            field_s => extract_scalar_field(states_new(state), trim(field_s%name))
            ! we need to append the state name here to make this safe for
            ! multi-material_phase/state... let's just hope you aren't going
            ! to try to pull this out of state by its name!
            call insert(meshes_new(unique_mesh_index), field_s, trim(states_new(state)%name)//"::"//trim(field_s%name))
          end if
        end do

        do field=1,vector_field_count(states_old(state))
          field_v => extract_vector_field(states_old(state), field)
          if (trim(field_v%mesh%name) == trim(mesh_name)) then
            ! we need to append the state name here to make this safe for
            ! multi-material_phase/state... let's just hope you aren't going
            ! to try to pull this out of state by its name!
            call insert(meshes_old(unique_mesh_index), field_v, trim(states_new(state)%name)//"::"//trim(field_v%name))
            field_v => extract_vector_field(states_new(state), trim(field_v%name))
            ! we need to append the state name here to make this safe for
            ! multi-material_phase/state... let's just hope you aren't going
            ! to try to pull this out of state by its name!
            call insert(meshes_new(unique_mesh_index), field_v, trim(states_new(state)%name)//"::"//trim(field_v%name))
          end if
        end do

        do field=1,tensor_field_count(states_old(state))
          field_t => extract_tensor_field(states_old(state), field)
          if (trim(field_t%mesh%name) == trim(mesh_name)) then
            ! we need to append the state name here to make this safe for
            ! multi-material_phase/state... let's just hope you aren't going
            ! to try to pull this out of state by its name!
            call insert(meshes_old(unique_mesh_index), field_t, trim(states_new(state)%name)//"::"//trim(field_t%name))
            field_t => extract_tensor_field(states_new(state), trim(field_t%name))
            ! we need to append the state name here to make this safe for
            ! multi-material_phase/state... let's just hope you aren't going
            ! to try to pull this out of state by its name!
            call insert(meshes_new(unique_mesh_index), field_t, trim(states_new(state)%name)//"::"//trim(field_t%name))
          end if
        end do

      end do
    end do    
    deallocate(unique_mesh_indices)
    deallocate(duplicate_mesh)

    ! Great! Now let's loop over the fields associated with each mesh
    ! and group them by algorithm.

    old_pos = extract_vector_field(states_old(1), "Coordinate")
    new_pos = extract_vector_field(states_new(1), "Coordinate")
   
    ! Do we need an element intersection map?
    if(.not. all_consistent_interpolation) then
      allocate(map_BA(ele_count(new_pos)))
      map_BA = intersection_finder(new_pos, old_pos)
    end if

    alg_loop: do alg = 1, alg_cnt
      ewrite(2, *) "  Considering algorithm " // trim(algorithms(alg))
      
      call tictoc_clear(TICTOC_ID_INTERPOLATION)
      call tic(TICTOC_ID_INTERPOLATION)

      select case(trim(algorithms(alg)))
        case("consistent_interpolation")
          do mesh = 1, mesh_cnt
            old_mesh => extract_mesh(meshes_old(mesh), "Mesh")
            new_mesh => extract_mesh(meshes_new(mesh), "Mesh")
            
            call insert(alg_old(mesh), old_mesh, "Mesh")
            call insert(alg_new(mesh), new_mesh, "Mesh")
            call insert(alg_old(mesh), old_pos, "Coordinate")
            call insert(alg_new(mesh), new_pos, "Coordinate")
            
            call collect_fields_to_interpolate(interpolate_field_consistent, meshes_new(mesh), meshes_old(mesh), alg_new(mesh), alg_old(mesh))
          
            if(field_count(alg_old(mesh)) > 1) then           
              if(present(map)) then
                ! Cannot assume here that "map" applies for new_mesh (as
                ! new_mesh may have any degree)
                if(size(map) == node_count(new_mesh)) then
                  call linear_interpolation(alg_old(mesh), alg_new(mesh), map = map)
                else
                  call linear_interpolation(alg_old(mesh), alg_new(mesh))
                end if
              else
                call linear_interpolation(alg_old(mesh), alg_new(mesh))
              end if
            end if
          
            call deallocate(alg_old(mesh))
            call deallocate(alg_new(mesh))
          end do
        case("interpolation_galerkin")
          do mesh = 1, mesh_cnt
            old_mesh => extract_mesh(meshes_old(mesh), "Mesh")
            new_mesh => extract_mesh(meshes_new(mesh), "Mesh")
            
            call insert(alg_old(mesh), old_mesh, "Mesh")
            call insert(alg_new(mesh), new_mesh, "Mesh")
            call insert(alg_old(mesh), old_pos, "Coordinate")
            call insert(alg_new(mesh), new_pos, "Coordinate")
          end do
          
          call collect_fields_to_interpolate(interpolate_field_galerkin_projection, meshes_new, meshes_old, alg_new, alg_old)
        
          no_fields = 0
          do mesh = 1, mesh_cnt
            no_fields = no_fields + field_count(alg_old(mesh))
          end do
          if(no_fields > mesh_cnt) then ! there will always be a Coordinate per mesh
            assert(allocated(map_BA))
            call interpolation_galerkin(alg_old, alg_new, map_BA = map_BA)
          end if
        
          do mesh=1,mesh_cnt
            call deallocate(alg_old(mesh))
            call deallocate(alg_new(mesh))
          end do
        case("geostrophic_interpolation")
          do state = 1, state_cnt
            do field = 1, vector_field_count(states_old(state))
              field_v_2 => extract_vector_field(states_old(state), field)
              if(.not. interpolate_field_geostrophic(field_v_2%option_path)) cycle
              if(field_v_2%name == "Coordinate") cycle
              field_v => extract_vector_field(states_new(state), field_v_2%name)
              call geostrophic_interpolation(states_old(state), field_v_2, &
                & states_new(state), field_v, &
                & map_BA = map_BA)
            end do
          end do
        case("grandy_interpolation")
          do mesh = 1, mesh_cnt
            old_mesh => extract_mesh(meshes_old(mesh), "Mesh")
            new_mesh => extract_mesh(meshes_new(mesh), "Mesh")
            
            call insert(alg_old(mesh), old_mesh, "Mesh")
            call insert(alg_new(mesh), new_mesh, "Mesh")
            call insert(alg_old(mesh), old_pos, "Coordinate")
            call insert(alg_new(mesh), new_pos, "Coordinate")
          end do
          
          call collect_fields_to_interpolate(interpolate_field_grandy_interpolation, meshes_new, meshes_old, alg_new, alg_old)
        
          no_fields = 0
          do mesh = 1, mesh_cnt
            no_fields = no_fields + field_count(alg_old(mesh))
          end do
          if(no_fields > mesh_cnt) then ! there will always be a Coordinate per mesh
            assert(allocated(map_BA))
            call grandy_projection(alg_old, alg_new, map_BA = map_BA)
          end if
        
          do mesh=1,mesh_cnt
            call deallocate(alg_old(mesh))
            call deallocate(alg_new(mesh))
          end do
          
        case("no_interpolation")
          ! nothing to be done obviously
        case default
          FLAbort("Unrecognised interpolation algorithm")
      end select

      call toc(TICTOC_ID_INTERPOLATION)
      call tictoc_report(2, TICTOC_ID_INTERPOLATION)

    end do alg_loop
    
    do mesh=1,mesh_cnt
      call deallocate(meshes_old(mesh))
      call deallocate(meshes_new(mesh))
    end do

    deallocate(meshes_old)
    deallocate(meshes_new)
    
    deallocate(alg_old)
    deallocate(alg_new)

    if(allocated(map_BA)) then
      do ele=1,size(map_BA)
        call deallocate(map_BA(ele))
      end do
      deallocate(map_BA)
    end if
    
    ewrite(1, *) "Exiting interpolate"
    
  end subroutine interpolate
  
  function interpolate_field_consistent(option_path) result(interpolate)
    character(len = *), intent(in) :: option_path

    logical :: interpolate
  
    interpolate = .false.
    if(len_trim(option_path) == 0) then
      interpolate = .true.
    else if(have_option(trim(complete_field_path(option_path)) // "/consistent_interpolation")) then
      interpolate = .true.
    end if
    
  end function interpolate_field_consistent
  
  function interpolate_field_galerkin_projection(option_path) result(interpolate)
    character(len = *), intent(in) :: option_path

    logical :: interpolate
    
    character(len = OPTION_PATH_LEN) :: base_path
   
    interpolate = .false.
    if(len_trim(option_path) == 0) return

    base_path = trim(complete_field_path(option_path))
    
    interpolate = have_option(trim(base_path) // "/galerkin_projection") &
      & .and. .not. have_option(trim(base_path) // "/galerkin_projection/supermesh_free")
    
  end function interpolate_field_galerkin_projection

  function interpolate_field_grandy_interpolation(option_path) result(interpolate)
    character(len = *), intent(in) :: option_path

    logical :: interpolate
    
    character(len = OPTION_PATH_LEN) :: base_path
   
    interpolate = .false.
    if(len_trim(option_path) == 0) return

    base_path = trim(complete_field_path(option_path))
    
    interpolate = have_option(trim(base_path) // "/grandy_interpolation")
  end function interpolate_field_grandy_interpolation
  
  function interpolate_field_geostrophic(option_path) result(interpolate)
    character(len = *), intent(in) :: option_path

    logical :: interpolate
    
    character(len = OPTION_PATH_LEN) :: base_path
   
    interpolate = .false.
    if(len_trim(option_path) == 0) return

    base_path = trim(complete_field_path(option_path))
    
    interpolate = have_option(trim(base_path) // "/geostrophic_interpolation")
    
  end function interpolate_field_geostrophic
  
  function interpolate_field_galerkin_projection_cg_supermesh_free(option_path) result(interpolate)
    character(len = *), intent(in) :: option_path

    logical :: interpolate
    
    character(len = OPTION_PATH_LEN) :: base_path
    
    interpolate = .false.
    if(len_trim(option_path) == 0) return

    base_path = trim(complete_field_path(option_path))
    
    interpolate = have_option(trim(base_path) // "/galerkin_projection/continuous") &
      & .and. have_option(trim(base_path) // "/galerkin_projection/supermesh_free")
      
#ifdef DDEBUG
    if(interpolate) then
      assert(.not. have_option(trim(base_path) // "/galerkin_projection/continuous/bounded"))
    end if
#endif
    
  end function interpolate_field_galerkin_projection_cg_supermesh_free
  
  function interpolate_field_galerkin_projection_dg_supermesh_free(option_path) result(interpolate)
    character(len = *), intent(in) :: option_path

    logical :: interpolate
    
    character(len = OPTION_PATH_LEN) :: base_path
    
    interpolate = .false.
    if(len_trim(option_path) == 0) return

    base_path = trim(complete_field_path(option_path))
    
    interpolate = have_option(trim(base_path) // "/galerkin_projection/discontinuous") &
      & .and. have_option(trim(base_path) // "/galerkin_projection/supermesh_free")
    
  end function interpolate_field_galerkin_projection_dg_supermesh_free
  
  subroutine collect_fields_to_interpolate_single_state(test, input_state_new, input_state_old, output_state_new, output_state_old)
    interface
      function test(option_path)
        implicit none
        character(len = *), intent(in) :: option_path
        logical :: test
      end function test
    end interface
    type(state_type), intent(inout) :: input_state_new
    type(state_type), intent(inout) :: input_state_old
    type(state_type), intent(inout) :: output_state_new
    type(state_type), intent(inout) :: output_state_old
    
    type(state_type), dimension(1) :: input_states_new
    type(state_type), dimension(1) :: input_states_old
    type(state_type), dimension(1) :: output_states_new
    type(state_type), dimension(1) :: output_states_old
    
    input_states_new = (/input_state_new/)
    input_states_old = (/input_state_old/)
    output_states_new = (/output_state_new/)
    output_states_old = (/output_state_old/)
    
    call collect_fields_to_interpolate(test, input_states_new, input_states_old, output_states_new, output_states_old)
    
    input_state_new = input_states_new(1)
    input_state_old = input_states_old(1)
    output_state_new = output_states_new(1)
    output_state_old = output_states_old(1)
    
  end subroutine collect_fields_to_interpolate_single_state

  subroutine collect_fields_to_interpolate_multiple_states(test, input_states_new, input_states_old, output_states_new, output_states_old)
    !!< Collect all fields in the supplied input_old and input_new states that
    !!< pass the supplied tests, and insert them into output_state_new and
    !!< output_state_old
    
    interface
      function test(option_path)
        implicit none
        character(len = *), intent(in) :: option_path
        logical :: test
      end function test
    end interface
    type(state_type), dimension(:), intent(in) :: input_states_new
    type(state_type), dimension(:), intent(in) :: input_states_old
    type(state_type), dimension(:), intent(inout) :: output_states_new
    type(state_type), dimension(:), intent(inout) :: output_states_old

    integer :: i, j
    type(scalar_field), pointer :: s_field_new, s_field_old
    type(tensor_field), pointer :: t_field_new, t_field_old
    type(vector_field), pointer :: v_field_new , v_field_old
    
#ifdef DDEBUG
    assert(size(input_states_new)==size(input_states_old))
    assert(size(output_states_new)==size(output_states_old))
    assert(size(input_states_new)==size(output_states_new))
#endif
    
    do j = 1, size(input_states_new)
    
      assert(scalar_field_count(input_states_new(j)) == scalar_field_count(input_states_old(j)))
      
      do i = 1, scalar_field_count(input_states_new(j))
        s_field_new => extract_scalar_field(input_states_new(j), i)
        s_field_old => extract_scalar_field(input_states_old(j), i)
#ifdef DDEBUG
        assert(trim(s_field_new%name) == trim(s_field_old%name))
        assert(trim(s_field_new%option_path) == trim(s_field_old%option_path))
#endif
        if(.not.aliased(s_field_new).and.test(s_field_new%option_path)) then
          ewrite(2, *) "    Found ", trim(s_field_new%name)
          ! make sure we keep the multi-material_phase/state safe names from state
          call insert(output_states_new(j), s_field_new, trim(input_states_new(j)%scalar_names(i)))
          call insert(output_states_old(j), s_field_old, trim(input_states_old(j)%scalar_names(i)))
        end if
      end do
      
      assert(vector_field_count(input_states_new(j)) == vector_field_count(input_states_old(j)))
      
      do i = 1, vector_field_count(input_states_new(j))
        v_field_new => extract_vector_field(input_states_new(j), i)
        v_field_old => extract_vector_field(input_states_old(j), i)
#ifdef DDEBUG
        assert(trim(v_field_new%name) == trim(v_field_old%name))
        assert(trim(v_field_new%option_path) == trim(v_field_old%option_path))
#endif
        if(v_field_new%name == "Coordinate") cycle
        if(.not.aliased(v_field_new).and.test(v_field_new%option_path)) then
          ewrite(2, *) "    Found ", trim(v_field_new%name)
          ! make sure we keep the multi-material_phase/state safe names from state
          call insert(output_states_new(j), v_field_new, trim(input_states_new(j)%vector_names(i)))
          call insert(output_states_old(j), v_field_old, trim(input_states_old(j)%vector_names(i)))
        end if
      end do
      
      assert(tensor_field_count(input_states_new(j)) == tensor_field_count(input_states_old(j)))
      
      do i = 1, tensor_field_count(input_states_new(j))
        t_field_new => extract_tensor_field(input_states_new(j), i)
        t_field_old => extract_tensor_field(input_states_old(j), i)
#ifdef DDEBUG
        assert(trim(t_field_new%name) == trim(t_field_old%name))
        assert(trim(t_field_new%option_path) == trim(t_field_old%option_path))
#endif
        if(.not.aliased(t_field_new).and.test(t_field_new%option_path)) then
          ewrite(2, *) "    Found ", trim(t_field_new%name)
          ! make sure we keep the multi-material_phase/state safe names from state
          call insert(output_states_new(j), t_field_new, trim(input_states_new(j)%tensor_names(i)))
          call insert(output_states_old(j), t_field_old, trim(input_states_old(j)%tensor_names(i)))
        end if
      end do
    end do
    
  end subroutine collect_fields_to_interpolate_multiple_states
  
  subroutine interpolation_manager_check_options()
    !!< Check interpolation algorithm selection options

    character(len = OPTION_PATH_LEN) :: base_path
    integer :: i, j

    ewrite(2, *) "Checking interpolation algorithm selection options"

    do i = 0, option_count("/material_phase") - 1
      do j = 0, option_count("/material_phase[" // int2str(i) // "]/scalar_field") - 1
        base_path = complete_field_path("/material_phase[" // int2str(i) // "]/scalar_field[" // int2str(j) // "]")
        if(have_option(trim(base_path) // "/galerkin_projection/supermesh_free")) then
          FLExit("Supermesh free Galerkin projection is not yet available")
        end if
      end do
      do j = 0, option_count("/material_phase[" // int2str(i) // "]/vector_field") - 1
        base_path = complete_field_path("/material_phase[" // int2str(i) // "]/vector_field[" // int2str(j) // "]")
        if(have_option(trim(base_path) // "/galerkin_projection/supermesh_free")) then
          FLExit("Supermesh free Galerkin projection is not yet available")
        end if
      end do
      do j = 0, option_count("/material_phase[" // int2str(i) // "]/tensor_field") - 1
        base_path = complete_field_path("/material_phase[" // int2str(i) // "]/tensor_field[" // int2str(j) // "]")
        if(have_option(trim(base_path) // "/galerkin_projection/supermesh_free")) then
          FLExit("Supermesh free Galerkin projection is not yet available")
        end if
      end do
    end do

    ewrite(2, *) "Finished checking interpolation algorithm selection options"

  end subroutine interpolation_manager_check_options

end module interpolation_manager
