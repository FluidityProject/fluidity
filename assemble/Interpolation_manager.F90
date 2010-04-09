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
  use global_parameters, only : OPTION_PATH_LEN, periodic_boundary_option_path
  use supermesh_construction
  use field_options
  use intersection_finder_module
  use linked_lists
  use boundary_conditions_from_options
  use tictoc
  use populate_state_module
  use vtk_interfaces
  use geostrophic_pressure
  
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
    ! For periodic interpolation, we need to interpolate the unwrapped fields,
    ! but keep the periodic fields too
    type(state_type), dimension(:), allocatable :: periodic_new
    ! Within each mesh, the fields organised by algorithm.
    type(state_type), dimension(:), allocatable :: alg_old, alg_new
    
    integer :: state, state_cnt, mesh
    integer :: mesh_cnt
    character(len=FIELD_NAME_LEN) :: mesh_name

    type(scalar_field), pointer :: field_s, p_field_s
    type(vector_field), pointer :: field_v, p_field_v
    type(tensor_field), pointer :: field_t, p_field_t
    integer :: field

    character(len=255), dimension(3), parameter :: algorithms = (/&
                       & "consistent_interpolation", &
                       & "interpolation_galerkin  ", &
                       & "grandy_interpolation    " /)
    integer :: alg_cnt, alg

    type(mesh_type), pointer :: old_mesh, new_mesh
    type(vector_field) :: old_pos, new_pos
    type(ilist), dimension(:), allocatable :: map_BA
    integer :: ele

    logical :: all_consistent_interpolation, all_linear_meshes, &
      any_periodic_meshes
    type(element_type), pointer :: field_shape

    integer :: no_fields
    
    integer :: stat

    ewrite(1, *) "In interpolate"

    call initialise_geostrophic_interpolation(states_old, states_new)

    mesh_cnt = option_count("/geometry/mesh")
    allocate(meshes_old(mesh_cnt))
    allocate(meshes_new(mesh_cnt))
    allocate(periodic_new(mesh_cnt))
    allocate(alg_old(mesh_cnt))
    allocate(alg_new(mesh_cnt))
    alg_cnt = size(algorithms)
    state_cnt = size(states_old)
    assert(size(states_old) == size(states_new))

    ! First thing: check the usual case. If all field request consistent
    ! interpolation, and all fields are on linear meshes, the just pass over to
    ! linear interpolation

    all_consistent_interpolation = .true.
    all_linear_meshes = .true.
    any_periodic_meshes = .false.

    consistent_linear_state_loop: do state=1,state_cnt
      do field=1,scalar_field_count(states_old(state))
        field_s => extract_scalar_field(states_old(state), field)
        ! If the field has no option path, assume consistent interpolation
        if(len_trim(field_s%option_path) /= 0) then
          all_consistent_interpolation = all_consistent_interpolation &
            .and. have_option(trim(complete_field_path(  &
               field_s%option_path, stat)) // "/consistent_interpolation")
          
          field_shape => ele_shape(field_s, 1)
          all_linear_meshes = all_linear_meshes .and. (field_shape%degree == 1)
          any_periodic_meshes = any_periodic_meshes .or. mesh_periodic(field_s)
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
          
          all_consistent_interpolation = all_consistent_interpolation &
            .and. have_option(trim(complete_field_path( &
               field_v%option_path, stat)) // "/consistent_interpolation")
          
          field_shape => ele_shape(field_v, 1)
          all_linear_meshes = all_linear_meshes .and. (field_shape%degree == 1)
          any_periodic_meshes = any_periodic_meshes .or. mesh_periodic(field_v)
        end if
      end do

      do field=1,tensor_field_count(states_old(state))
        field_t => extract_tensor_field(states_old(state), field)
        ! If the field has no option path, assume consistent interpolation
        if(len_trim(field_t%option_path) /= 0) then
        
          all_consistent_interpolation = all_consistent_interpolation &
            .and. have_option(trim(complete_field_path( &
               field_t%option_path, stat)) // "/consistent_interpolation")
          
          field_shape => ele_shape(field_t, 1)
          all_linear_meshes = all_linear_meshes .and. (field_shape%degree == 1)
          any_periodic_meshes = any_periodic_meshes .or. mesh_periodic(field_t)
        end if
      end do

    end do consistent_linear_state_loop
    
    if(all_consistent_interpolation &
      .and. all_linear_meshes .and. .not. any_periodic_meshes) then
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
    
    old_pos = extract_vector_field(states_old(1), "Coordinate")
    new_pos = extract_vector_field(states_new(1), "Coordinate")
    
    ! OK! So we have some work to do.
    ! First, let's organise the fields according to what mesh
    ! they are on.
    do mesh=1,mesh_cnt

      call get_option("/geometry/mesh["//int2str(mesh-1)//"]/name", mesh_name)

      do state=1,state_cnt
        do field=1,scalar_field_count(states_old(state))
          field_s => extract_scalar_field(states_old(state), field)
          if (trim(field_s%mesh%name) == trim(mesh_name)) then
            ! we need to append the state name here to make this safe for
            ! multi-material_phase/state... let's just hope you aren't going
            ! to try to pull this out of state by its name!
            call insert(meshes_old(mesh), field_s, trim(states_new(state)%name)//"::"//trim(field_s%name))
            field_s => extract_scalar_field(states_new(state), trim(field_s%name))
            ! we need to append the state name here to make this safe for
            ! multi-material_phase/state... let's just hope you aren't going
            ! to try to pull this out of state by its name!
            call insert(meshes_new(mesh), field_s, trim(states_new(state)%name)//"::"//trim(field_s%name))
          end if
        end do

        do field=1,vector_field_count(states_old(state))
          field_v => extract_vector_field(states_old(state), field)
          if (trim(field_v%mesh%name) == trim(mesh_name)) then
            
            if (field_v%name=="Coordinate" .or. field_v%name==trim(mesh_name)//"Coordinate") cycle
            ! we need to append the state name here to make this safe for
            ! multi-material_phase/state... let's just hope you aren't going
            ! to try to pull this out of state by its name!
            call insert(meshes_old(mesh), field_v, trim(states_new(state)%name)//"::"//trim(field_v%name))
            field_v => extract_vector_field(states_new(state), trim(field_v%name))
            ! we need to append the state name here to make this safe for
            ! multi-material_phase/state... let's just hope you aren't going
            ! to try to pull this out of state by its name!
            call insert(meshes_new(mesh), field_v, trim(states_new(state)%name)//"::"//trim(field_v%name))
          end if
        end do

        do field=1,tensor_field_count(states_old(state))
          field_t => extract_tensor_field(states_old(state), field)
          if (trim(field_t%mesh%name) == trim(mesh_name)) then
            ! we need to append the state name here to make this safe for
            ! multi-material_phase/state... let's just hope you aren't going
            ! to try to pull this out of state by its name!
            call insert(meshes_old(mesh), field_t, trim(states_new(state)%name)//"::"//trim(field_t%name))
            field_t => extract_tensor_field(states_new(state), trim(field_t%name))
            ! we need to append the state name here to make this safe for
            ! multi-material_phase/state... let's just hope you aren't going
            ! to try to pull this out of state by its name!
            call insert(meshes_new(mesh), field_t, trim(states_new(state)%name)//"::"//trim(field_t%name))
          end if
        end do

      end do

      old_mesh => extract_mesh(states_old(1), trim(mesh_name))
      call insert(meshes_old(mesh), old_mesh, "Mesh")

      new_mesh => extract_mesh(states_new(1), trim(mesh_name))
      call insert(meshes_new(mesh), new_mesh, "Mesh")
      
      call insert(meshes_old(mesh), old_pos, "Coordinate")
      call insert(meshes_new(mesh), new_pos, "Coordinate")

      if (mesh_periodic(new_mesh)) then
        
        ! Nice work, soldier! Now if the mesh is periodic we need to expand the mesh in the plane,
        ! as if we've adapted the domain edges won't match up, will they?
        
        ! first we keep a copy of meshes_new
        do field=1,scalar_field_count(meshes_new(mesh))
          field_s => extract_scalar_field(meshes_new(mesh), field)
          call insert(periodic_new(mesh), field_s, trim(field_s%name))
        end do

        do field=1,vector_field_count(meshes_new(mesh))
          field_v => extract_vector_field(meshes_new(mesh), field)
          call insert(periodic_new(mesh), field_v, trim(field_v%name))
        end do

        do field=1,tensor_field_count(meshes_new(mesh))
          field_t => extract_tensor_field(meshes_new(mesh), field)
          call insert(periodic_new(mesh), field_t, trim(field_t%name))
        end do

        ! create an expanded, non-periodic version of meshes_old
        ! and a non-periodic version of meshes_new, and remap all fields
        ! to it
        call prepare_periodic_states_for_interpolation(meshes_old(mesh), meshes_new(mesh))
      end if
      
    end do

    ! Great! Now let's loop over the fields associated with each mesh
    ! and group them by algorithm.
    
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
            call get_option("/geometry/mesh[" // int2str(mesh - 1) // "]/name", mesh_name)
            old_mesh => extract_mesh(states_old(1), trim(mesh_name))
            new_mesh => extract_mesh(states_new(1), trim(mesh_name))
            old_pos = extract_vector_field(meshes_old(mesh), "Coordinate")
            new_pos = extract_vector_field(meshes_new(mesh), "Coordinate")
            
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
            call get_option("/geometry/mesh[" // int2str(mesh - 1) // "]/name", mesh_name)
            old_mesh => extract_mesh(states_old(1), trim(mesh_name))
            new_mesh => extract_mesh(states_new(1), trim(mesh_name))
            old_pos = extract_vector_field(meshes_old(mesh), "Coordinate")
            new_pos = extract_vector_field(meshes_new(mesh), "Coordinate")
            
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
        case("grandy_interpolation")
          do mesh = 1, mesh_cnt
            call get_option("/geometry/mesh[" // int2str(mesh - 1) // "]/name", mesh_name)
            old_mesh => extract_mesh(states_old(1), trim(mesh_name))
            new_mesh => extract_mesh(states_new(1), trim(mesh_name))
            old_pos = extract_vector_field(meshes_old(mesh), "Coordinate")
            new_pos = extract_vector_field(meshes_new(mesh), "Coordinate")
            
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
    
    if (any_periodic_meshes) then
      do mesh=1,mesh_cnt
        call get_option("/geometry/mesh["//int2str(mesh-1)//"]/name", mesh_name)
        old_mesh => extract_mesh(states_old(1), trim(mesh_name))

        if (.not. mesh_periodic(old_mesh)) cycle

        ! If we are periodic, we have interpolated the unwrapped version
        ! So let's remap to the periodic one we actually need
        do field=1,scalar_field_count(meshes_new(mesh))
          field_s => extract_scalar_field(meshes_new(mesh), field)
          p_field_s => extract_scalar_field(periodic_new(mesh), trim(field_s%name))
          assert(trim(field_s%name) == trim(p_field_s%name))
          call remap_field(field_s, p_field_s, stat=stat)
        end do

        do field=1,vector_field_count(meshes_new(mesh))
          field_v => extract_vector_field(meshes_new(mesh), field)
          if (trim(field_v%name) == "Coordinate") cycle
          p_field_v => extract_vector_field(periodic_new(mesh), trim(field_v%name))
          assert(trim(field_v%name) == trim(p_field_v%name))
          call remap_field(field_v, p_field_v, stat=stat)
        end do

        do field=1,tensor_field_count(meshes_new(mesh))
          field_t => extract_tensor_field(meshes_new(mesh), field)
          p_field_t => extract_tensor_field(periodic_new(mesh), trim(field_v%name))
          assert(trim(field_t%name) == trim(p_field_t%name))
          call remap_field(field_t, p_field_t, stat=stat)
        end do

        call deallocate(periodic_new(mesh))
      end do
    end if
    
    do mesh=1, mesh_cnt
      call deallocate(meshes_old(mesh))
      call deallocate(meshes_new(mesh))
    end do

    deallocate(meshes_old)
    deallocate(meshes_new)
    deallocate(periodic_new)
    
    deallocate(alg_old)
    deallocate(alg_new)

    if(allocated(map_BA)) then
      do ele=1,size(map_BA)
        call deallocate(map_BA(ele))
      end do
      deallocate(map_BA)
    end if

    call finalise_geostrophic_interpolation(states_new)
    
    ewrite(1, *) "Exiting interpolate"

  end subroutine interpolate
  
  function interpolate_field_consistent(option_path) result(interpolate)
    character(len = *), intent(in) :: option_path

    logical :: interpolate
    
    integer :: stat
  
    interpolate = .false.
    if(len_trim(option_path) == 0) then
      interpolate = .true.
    else if(have_option(trim(complete_field_path(option_path, stat = stat)) // "/consistent_interpolation")) then
      interpolate = .true.
    end if
    
  end function interpolate_field_consistent
  
  function interpolate_field_galerkin_projection(option_path) result(interpolate)
    character(len = *), intent(in) :: option_path

    logical :: interpolate
    
    character(len = OPTION_PATH_LEN) :: base_path
    integer :: stat
   
    interpolate = .false.
    if(len_trim(option_path) == 0) return

    base_path = trim(complete_field_path(option_path, stat = stat))
    
    interpolate = have_option(trim(base_path) // "/galerkin_projection") &
      & .and. .not. have_option(trim(base_path) // "/galerkin_projection/supermesh_free")
    
  end function interpolate_field_galerkin_projection

  function interpolate_field_grandy_interpolation(option_path) result(interpolate)
    character(len = *), intent(in) :: option_path

    logical :: interpolate
    
    character(len = OPTION_PATH_LEN) :: base_path
    integer :: stat
   
    interpolate = .false.
    if(len_trim(option_path) == 0) return

    base_path = trim(complete_field_path(option_path, stat = stat))
    
    interpolate = have_option(trim(base_path) // "/grandy_interpolation")
  end function interpolate_field_grandy_interpolation
  
  function interpolate_field_galerkin_projection_cg_supermesh_free(option_path) result(interpolate)
    character(len = *), intent(in) :: option_path

    logical :: interpolate
    
    character(len = OPTION_PATH_LEN) :: base_path
    integer :: stat
    
    interpolate = .false.
    if(len_trim(option_path) == 0) return

    base_path = trim(complete_field_path(option_path, stat = stat))
    
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
    integer :: stat
    
    interpolate = .false.
    if(len_trim(option_path) == 0) return

    base_path = trim(complete_field_path(option_path, stat = stat))
    
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
        if(index(trim(v_field_new%name), "Coordinate") /= 0) cycle
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

  subroutine prepare_periodic_states_for_interpolation(mesh_old, mesh_new)
    type(state_type), intent(inout) :: mesh_old, mesh_new

    type(vector_field), pointer :: positions_new
    type(mesh_type), pointer :: new_mesh
    type(mesh_type) :: new_mesh_unperiodic

    type(vector_field) :: expanded_positions
    type(vector_field), pointer :: positions_old
    type(mesh_type), pointer :: old_mesh
    type(mesh_type) :: old_mesh_unperiodic, old_mesh_expanded

    integer :: field
    type(scalar_field), pointer :: p_field_s, u_field_s
    type(scalar_field) :: field_s
    type(vector_field), pointer :: p_field_v, u_field_v
    type(vector_field) :: field_v
    type(tensor_field), pointer :: p_field_t, u_field_t
    type(tensor_field) :: field_t
      
    type(element_type), pointer :: x_shape, mesh_shape

    integer :: j, node, total_multiple

    type(state_type) :: unwrapped_state

    integer :: dim

    ! Easiest first. Unwrap the coordinate of mesh_new, reallocate all the fields,
    ! and remap periodic -> nonperiodic.

    positions_new => extract_vector_field(mesh_new, "Coordinate")
    new_mesh => extract_mesh(mesh_new, "Mesh")
    dim = positions_new%dim
    
    x_shape => positions_new%mesh%shape
    mesh_shape => new_mesh%shape
    if (mesh_shape%degree/=x_shape%degree .or. &
      continuity(new_mesh)/=0) then
      ! make a non-periodic version of the mesh
      new_mesh_unperiodic = make_mesh(positions_new%mesh, mesh_shape, &
        continuity=continuity(new_mesh), name=trim(new_mesh%name)//"Unperiodic")
    else
      new_mesh_unperiodic = positions_new%mesh
      call incref(new_mesh_unperiodic)
    end if
    
    do field=1,scalar_field_count(mesh_new)
      p_field_s => extract_scalar_field(mesh_new, field)
      call allocate(field_s, new_mesh_unperiodic, trim(p_field_s%name))
      field_s%option_path = p_field_s%option_path
      call remap_field(p_field_s, field_s)
      call insert(mesh_new, field_s, trim(mesh_new%scalar_names(field)))
      call deallocate(field_s)
    end do

    do field=1,vector_field_count(mesh_new)
      p_field_v => extract_vector_field(mesh_new, field)
      if (trim(p_field_v%name) == "Coordinate") then
        cycle
      end if
      call allocate(field_v, p_field_v%dim, new_mesh_unperiodic, trim(p_field_v%name))
      field_v%option_path = p_field_v%option_path
      call remap_field(p_field_v, field_v)
      call insert(mesh_new, field_v, trim(mesh_new%vector_names(field)))
      call deallocate(field_v)
    end do

    do field=1,tensor_field_count(mesh_new)
      p_field_t => extract_tensor_field(mesh_new, field)
      call allocate(field_t, new_mesh_unperiodic, trim(p_field_t%name))
      field_t%option_path = p_field_t%option_path
      call remap_field(p_field_t, field_t)
      call insert(mesh_new, field_t, trim(mesh_new%tensor_names(field)))
      call deallocate(field_t)
    end do

    ! replace mesh with non-periodic version:
    call insert(mesh_new, new_mesh_unperiodic, "Mesh")
    call deallocate(new_mesh_unperiodic)

    ! OK. Let's do the same for the old fields, with the additional twist
    ! that we duplicate the mesh on either side of the periodic boundary.
    ! Since the boundaries of the new domain don't match up with the old one,
    ! this is necessary so the new domain can find the relevant information.

    positions_old => extract_vector_field(mesh_old, "Coordinate")
    old_mesh => extract_mesh(mesh_old, "Mesh")
    dim = positions_old%dim
    
    x_shape => positions_old%mesh%shape
    mesh_shape => old_mesh%shape
    if (mesh_shape%degree/=x_shape%degree .or. &
      continuity(old_mesh)/=0) then
      ! make a non-periodic version of the mesh
      old_mesh_unperiodic = make_mesh(positions_old%mesh, mesh_shape, &
        continuity=continuity(old_mesh), name=trim(old_mesh%name)//"Unperiodic")
    else
      old_mesh_unperiodic = positions_old%mesh
      call incref(old_mesh_unperiodic)
    end if
    
    call expand_periodic_mesh(positions_old, old_mesh_unperiodic, &
      expanded_positions, old_mesh_expanded)
    
    total_multiple = element_count(old_mesh_expanded)/element_count(old_mesh_unperiodic)
    assert(element_count(old_mesh_expanded) == total_multiple*element_count(old_mesh_unperiodic))
        
    ! Let's unwrap everything and put it into the unwrapped_state.
    do field=1,scalar_field_count(mesh_old)
      p_field_s => extract_scalar_field(mesh_old, field)
      call allocate(field_s, old_mesh_unperiodic, trim(p_field_s%name))
      call remap_field(p_field_s, field_s)
      call insert(unwrapped_state, field_s, trim(p_field_s%name))
      call deallocate(field_s)
    end do

    do field=1,vector_field_count(mesh_old)
      p_field_v => extract_vector_field(mesh_old, field)
      call allocate(field_v, p_field_v%dim, old_mesh_unperiodic, trim(p_field_v%name))
      call remap_field(p_field_v, field_v)
      call insert(unwrapped_state, field_v, trim(p_field_v%name))
      call deallocate(field_v)
    end do

    do field=1,tensor_field_count(mesh_old)
      p_field_t => extract_tensor_field(mesh_old, field)
      call allocate(field_t, old_mesh_unperiodic, trim(p_field_t%name))
      call remap_field(p_field_t, field_t)
      call insert(unwrapped_state, field_t, trim(p_field_t%name))
      call deallocate(field_t)
    end do

    ! now copy the unwrapped fields into expanded versions, and put them
    ! in mesh_old states
    
    do field=1,scalar_field_count(mesh_old)
      p_field_s => extract_scalar_field(mesh_old, field)
      u_field_s => extract_scalar_field(unwrapped_state, trim(p_field_s%name))
      call allocate(field_s, old_mesh_expanded, trim(p_field_s%name))
      field_s%option_path = p_field_s%option_path
      do node=1,node_count(u_field_s)
        do j=0,total_multiple-1
          call set(field_s, node + j*node_count(u_field_s), node_val(u_field_s, node))
        end do
      end do
      call insert(mesh_old, field_s, trim(mesh_old%scalar_names(field)))
      call deallocate(field_s)
    end do

    do field=1,vector_field_count(mesh_old)
      p_field_v => extract_vector_field(mesh_old, field)
      if (trim(p_field_v%name) == "Coordinate") then
        cycle
      end if
      u_field_v => extract_vector_field(unwrapped_state, trim(p_field_v%name))
      call allocate(field_v, p_field_v%dim, old_mesh_expanded, trim(p_field_v%name))
      field_v%option_path = p_field_v%option_path
      do node=1,node_count(u_field_v)
        do j=0,total_multiple-1
          call set(field_v, node + j*node_count(u_field_v), node_val(u_field_v, node))
        end do
      end do
      call insert(mesh_old, field_v, trim(mesh_old%vector_names(field)))
      call deallocate(field_v)
    end do

    do field=1,tensor_field_count(mesh_old)
      p_field_t => extract_tensor_field(mesh_old, field)
      u_field_t => extract_tensor_field(unwrapped_state, trim(p_field_t%name))
      call allocate(field_t, old_mesh_expanded, trim(p_field_t%name))
      field_t%option_path = p_field_t%option_path
      do node=1,node_count(u_field_t)
        do j=0,total_multiple-1
          call set(field_t, node + j*node_count(u_field_t), node_val(u_field_t, node))
        end do
      end do
      call insert(mesh_old, field_t, trim(mesh_old%tensor_names(field)))
      call deallocate(field_t)
    end do
      
    ! replace mesh and positions with expanded non-periodic version:
    call insert(mesh_old, old_mesh_expanded, "Mesh")
    call insert(mesh_old, expanded_positions, "Coordinate")
    
    ! debugging vtus:
    !call insert(unwrapped_state, old_mesh_unperiodic, "Mesh")
    !call vtk_write_state("unwrapped", 0, model="Mesh", state=(/unwrapped_state/))
    !call vtk_write_state("mesh_new", 0, model="Mesh", state=(/mesh_new/))
    !call vtk_write_state("mesh_old", 0, model="Mesh", state=(/mesh_old/))
      
    call deallocate(old_mesh_unperiodic)
    call deallocate(old_mesh_expanded)
    call deallocate(expanded_positions)
    call deallocate(unwrapped_state)

  end subroutine prepare_periodic_states_for_interpolation
    
  subroutine expand_periodic_mesh(positions_in, mesh_in, &
    expanded_positions, expanded_mesh)
    ! creates an expanded version of the provided positions_in field
    ! and the (possibly different) mesh_in, adding periodic copies 
    ! by applying the periodic (inverse) maps
    type(vector_field), intent(in):: positions_in
    type(mesh_type), intent(in):: mesh_in
    type(vector_field), intent(out):: expanded_positions
    type(mesh_type), intent(out):: expanded_mesh
    
    type(vector_field) :: positions
    type(mesh_type) :: mesh, expanded_x_mesh
    character(len=OPTION_PATH_LEN) :: periodic_mapping_python, inverse_periodic_mapping_python
    real, dimension(:, :), allocatable :: forward_mapped_nodes, inverse_mapped_nodes
    integer :: multiple, j, bc, no_bcs, ele, node, total_multiple
    integer :: dim
    
    dim = positions_in%dim
    
    no_bcs = option_count(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions')
    total_multiple = 1
    
    positions = positions_in
    mesh = mesh_in
    call incref(positions)
    call incref(mesh)
    
    do bc=0,no_bcs-1

      ! figure out how many copies to make
      multiple = 2 ! for the current mesh, as well as the one under the forward periodic mapping
      if (have_option(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(bc)//']/inverse_coordinate_map')) then
        multiple = multiple + 1 ! plus the inverse periodic mapping, if we have it
      end if
      total_multiple = total_multiple * multiple

      ! expand the coordinate mesh
      call allocate(expanded_x_mesh, multiple * node_count(positions), multiple * ele_count(positions), &
                 &  positions%mesh%shape, positions%mesh%name)
      do ele=1,ele_count(positions)
        do j=0,multiple-1
          call set_ele_nodes(expanded_x_mesh, ele + j*ele_count(positions), ele_nodes(positions, ele) + j*node_count(positions))
        end do
      end do
      
      ! do the same for the provided 'mesh_in'
      if (positions_in%mesh==mesh_in) then
        expanded_mesh = expanded_x_mesh
        call incref(expanded_mesh)
      else        
        call allocate(expanded_mesh, multiple * node_count(mesh), multiple * ele_count(mesh), &
                 &  mesh%shape, mesh%name)
        expanded_mesh%continuity=mesh%continuity
        do ele=1,ele_count(mesh)
          do j=0,multiple-1
            call set_ele_nodes(expanded_mesh, ele + j*ele_count(mesh), ele_nodes(mesh, ele) + j*node_count(mesh))
          end do
        end do        
      end if
      
      call allocate(expanded_positions, dim, expanded_x_mesh, positions%name)
      call deallocate(expanded_x_mesh)

      do node=1, node_count(positions)
        call set(expanded_positions, node, node_val(positions, node))
      end do

      call get_option(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(bc)//']/coordinate_map', periodic_mapping_python)
      allocate(forward_mapped_nodes(dim, node_count(positions)))
      call set_from_python_function(forward_mapped_nodes, periodic_mapping_python, positions, time=0.0)
      do node=1, node_count(positions)
        call set(expanded_positions, node + node_count(positions), forward_mapped_nodes(:, node))
      end do
      deallocate(forward_mapped_nodes)

      if (have_option(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(bc)//']/inverse_coordinate_map')) then
        call get_option(trim(periodic_boundary_option_path(dim)) // '/from_mesh/periodic_boundary_conditions['//int2str(bc)//']/inverse_coordinate_map', inverse_periodic_mapping_python)
        allocate(inverse_mapped_nodes(dim, node_count(positions)))
        call set_from_python_function(inverse_mapped_nodes, inverse_periodic_mapping_python, positions, time=0.0)
        do node=1,node_count(positions)
          call set(expanded_positions, node + 2*node_count(positions), inverse_mapped_nodes(:, node))
        end do
        deallocate(inverse_mapped_nodes)
      end if

      call deallocate(mesh)
      call deallocate(positions)
      positions = expanded_positions
      mesh = expanded_mesh
    end do

  end subroutine expand_periodic_mesh

end module interpolation_manager
