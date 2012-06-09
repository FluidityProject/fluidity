#include "fdebug.h"

module lebiology_python
  use fldebug
  use fields
  use detector_data_types
  use detector_tools
  use detector_parallel
  use ieee_arithmetic, only: ieee_is_nan

  implicit none
  
  private
  
  public :: lebiology_init_module, lebiology_set_stage_id, &
            lebiology_add_variables, lebiology_add_envfields, &
            lebiology_add_foods, lebiology_prepare_pyfunc, &
            lebiology_initialise_agent, lebiology_move_agent, &
            lebiology_update_agent, get_new_agent_list

  type(detector_linked_list), target, save :: new_agent_list

  interface

    subroutine lebiology_init_module() bind(c, name='initlebiology')
      use :: iso_c_binding
      implicit none

    end subroutine lebiology_init_module

    subroutine lebiology_add_fg_varname(fg, fglen, var, varlen, stat) &
           bind(c, name='lebiology_add_fg_varname_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: fglen, varlen
      character(kind=c_char), dimension(fglen), intent(in) :: fg
      character(kind=c_char), dimension(varlen), intent(in) :: var
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_add_fg_varname

    subroutine lebiology_add_fg_envname(fg, fglen, env, envlen, stat) &
           bind(c, name='lebiology_add_fg_envname_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: fglen, envlen
      character(kind=c_char), dimension(fglen), intent(in) :: fg
      character(kind=c_char), dimension(envlen), intent(in) :: env
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_add_fg_envname

    subroutine lebiology_add_fg_foodname(fg, fglen, food, foodlen, variety, varietylen, stat) &
           bind(c, name='lebiology_add_fg_foodname_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: fglen, foodlen, varietylen
      character(kind=c_char), dimension(fglen), intent(in) :: fg
      character(kind=c_char), dimension(foodlen), intent(in) :: food
      character(kind=c_char), dimension(varietylen), intent(in) :: variety
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_add_fg_foodname

    subroutine lebiology_add_fg_stage_id(fg, fglen, stage, stagelen, id, stat) &
           bind(c, name='lebiology_add_fg_stage_id_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: fglen, stagelen
      character(kind=c_char), dimension(fglen), intent(in) :: fg
      character(kind=c_char), dimension(stagelen), intent(in) :: stage
      real(c_double), intent(in) :: id
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_add_fg_stage_id

    subroutine lebiology_compile_function(fg, fglen, key, keylen, func, funclen, stat) &
           bind(c, name='lebiology_compile_function_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: fglen, keylen, funclen
      character(kind=c_char), dimension(fglen), intent(in) :: fg
      character(kind=c_char), dimension(keylen), intent(in) :: key
      character(kind=c_char), dimension(funclen), intent(in) :: func
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_compile_function

    subroutine lebiology_agent_init(fg, fglen, key, keylen, vars, n_vars, stat) &
           bind(c, name='lebiology_agent_init_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: fglen, keylen, n_vars
      character(kind=c_char), dimension(fglen), intent(in) :: fg
      character(kind=c_char), dimension(keylen), intent(in) :: key
      real(c_double), dimension(n_vars), intent(inout) :: vars
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_agent_init

    subroutine lebiology_agent_update(fg, fglen, key, keylen, food, foodlen, &
           vars, n_vars, envvals, n_envvals, fvariety, frequest, fthreshold, fingest, n_fvariety, dt, stat) &
           bind(c, name='lebiology_agent_update_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: fglen, keylen, foodlen
      integer(c_int), intent(in), value :: n_vars, n_envvals, n_fvariety
      character(kind=c_char), dimension(fglen), intent(in) :: fg
      character(kind=c_char), dimension(keylen), intent(in) :: key
      character(kind=c_char), dimension(foodlen), intent(in) :: food
      real(c_double), dimension(n_vars), intent(inout) :: vars
      real(c_double), dimension(n_envvals), intent(inout) :: envvals
      real(c_double), dimension(n_fvariety), intent(inout) :: fvariety, frequest, fingest, fthreshold
      real(c_double), intent(in) :: dt
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_agent_update

    subroutine lebiology_agent_move(fg, fglen, key, keylen, pos, n_pos, &
           vars, n_vars, var_inds, dt, vector, stat) &
           bind(c, name='lebiology_agent_move_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: fglen, keylen, n_pos, n_vars
      character(kind=c_char), dimension(fglen), intent(in) :: fg
      character(kind=c_char), dimension(keylen), intent(in) :: key
      real(c_double), dimension(n_pos), intent(inout) :: pos
      real(c_double), dimension(n_vars), intent(inout) :: vars
      integer(c_int), dimension(n_vars), intent(inout) :: var_inds
      real(c_double), intent(in) :: dt
      real(c_double), dimension(n_pos), intent(out) :: vector
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_agent_move

  end interface

contains

  function get_new_agent_list() result(agent_list_ptr)
    type(detector_linked_list), pointer :: agent_list_ptr

    agent_list_ptr => new_agent_list
  end function get_new_agent_list

  subroutine lebiology_add_variables(fgroup)
    type(functional_group), intent(inout) :: fgroup
    integer :: v, stat

    stat=0
    do v=1, size(fgroup%variables)
       call lebiology_add_fg_varname(trim(fgroup%name), len_trim(fgroup%name), &
             trim(fgroup%variables(v)%name), len_trim(fgroup%variables(v)%name), stat)
    end do

    if (stat < 0) then
       ewrite(-1, *) "Error setting variable names for FG::"//trim(fgroup%name)
       FLExit("Python error in LE-Biology")
    end if
  end subroutine lebiology_add_variables

  subroutine lebiology_add_envfields(fgroup)
    type(functional_group), intent(inout) :: fgroup
    integer :: f, stat

    stat=0
    do f=1, size(fgroup%envfield_names)
       call lebiology_add_fg_envname(trim(fgroup%name), len_trim(fgroup%name), &
             trim(fgroup%envfield_names(f)), len_trim(fgroup%envfield_names(f)), stat)
    end do

    if (stat < 0) then
       ewrite(-1, *) "Error setting environment field names for FG::"//trim(fgroup%name)
       FLExit("Python error in LE-Biology")
    end if
  end subroutine lebiology_add_envfields

  subroutine lebiology_add_foods(fgroup)
    type(functional_group), intent(inout) :: fgroup
    type(food_set) :: fset
    integer :: f, s, stat

    stat=0
    do f=1, size(fgroup%food_sets)
       fset = fgroup%food_sets(f)
       do s=1, size(fset%varieties)
          call lebiology_add_fg_foodname(trim(fgroup%name), len_trim(fgroup%name), &
                 trim(fset%name), len_trim(fset%name), trim(fset%varieties(s)%name), &
                 len_trim(fset%varieties(s)%name), stat)
       end do
    end do

    if (stat < 0) then
       ewrite(-1, *) "Error setting food variety names for FG::"//trim(fgroup%name)
       FLExit("Python error in LE-Biology")
    end if
  end subroutine lebiology_add_foods

  subroutine lebiology_set_stage_id(fgroup, stage_name, id)
    type(functional_group), intent(inout) :: fgroup
    character(len=*), intent(in) :: stage_name
    real, intent(in) :: id
    integer :: s, stat

    stat=0
    call lebiology_add_fg_stage_id(trim(fgroup%name), len_trim(fgroup%name), &
             trim(stage_name), len_trim(stage_name), id, stat)

    if (stat < 0) then
       ewrite(-1, *) "Error setting stage ID for FG::"//trim(fgroup%name)
       FLExit("Python error in LE-Biology")
    end if
  end subroutine lebiology_set_stage_id

  subroutine lebiology_prepare_pyfunc(fgroup, key, func)
    type(functional_group), intent(inout) :: fgroup
    character(len=*), intent(in) :: key, func
    integer :: stat

    stat=0
    call lebiology_compile_function(trim(fgroup%name), len_trim(fgroup%name), &
             trim(key), len_trim(key), trim(func), len_trim(func), stat)

    if (stat < 0) then
       ewrite(-1, *) "Error compiling for FG::"//trim(fgroup%name)
       FLExit("Python error in LE-Biology")
    end if
  end subroutine lebiology_prepare_pyfunc

  subroutine lebiology_initialise_agent(fgroup, key, agent)
    type(functional_group), intent(inout) :: fgroup
    character(len=*), intent(in) :: key
    type(detector_type), intent(inout) :: agent
    integer :: stat

    stat=0
    call lebiology_agent_init(trim(fgroup%name), len_trim(fgroup%name), &
             trim(key), len_trim(key), agent%biology, size(agent%biology), stat)

    if (stat < 0) then
       ewrite(-1, *) "Error initialising agent for FG::"//trim(fgroup%name)
       FLExit("Python error in LE-Biology")
    end if    
  end subroutine lebiology_initialise_agent

  subroutine lebiology_move_agent(fgroup, key, agent, dt, vector)
    type(functional_group), intent(inout) :: fgroup
    character(len=*), intent(in) :: key
    type(detector_type), intent(inout) :: agent
    real, intent(in) :: dt
    real, dimension(size(agent%position)), intent(out) :: vector

    real, dimension(size(fgroup%motion_var_inds)) :: state_vars
    integer :: v, stat

    if (size(fgroup%motion_var_inds) > 0) then
       do v=1, size(fgroup%motion_var_inds)
          state_vars(v) = agent%biology( fgroup%motion_var_inds(v) )
       end do
    end if

    stat=0
    call lebiology_agent_move(trim(fgroup%name), len_trim(fgroup%name), &
            trim(key), len_trim(key), agent%update_vector, size(agent%update_vector), &
            state_vars, size(state_vars), fgroup%motion_var_inds, dt, vector, stat) 

    if (size(fgroup%motion_var_inds) > 0) then
       do v=1, size(fgroup%motion_var_inds)
          if (ieee_is_nan(state_vars(v))) then
             FLExit('NaN motion variable detected in '//trim(fgroup%name)//"::"//trim(key))
          end if
          agent%biology( fgroup%motion_var_inds(v) ) = state_vars(v)
       end do
    end if

    if (stat < 0) then
       ewrite(-1, *) "Error moving agent "//int2str(agent%id_number)//" for FG::"//trim(fgroup%name)
       FLExit("Python error in LE-Biology")
    end if
  end subroutine lebiology_move_agent

  subroutine lebiology_update_agent(fgroup, key, foodname, agent, xfield, envfields, foodfields, dt)
    type(functional_group), intent(inout) :: fgroup
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: foodname
    type(detector_type), intent(inout) :: agent
    type(vector_field), pointer, intent(inout) :: xfield
    type(scalar_field_pointer), dimension(:), pointer, intent(inout) :: envfields
    type(scalar_field_pointer), dimension(:), pointer, intent(inout) :: foodfields
    real, intent(in) :: dt

    real, dimension(size(envfields)) :: envfield_vals
    real, dimension(size(foodfields)) :: foodfield_vals
    real, dimension(:), allocatable :: foodval_ele
    real :: path_total
    real, dimension(size(agent%biology)) :: agent_state_copy
    integer :: f, v, e, stat

    agent_state_copy = agent%biology

    ! Sample environment fields
    do f=1, size(envfields)
       ! Evaluate at detector position
       envfield_vals(f) = eval_field(agent%element, envfields(f)%ptr, agent%local_coords)
    end do

    ! Add food concentrations
    do f=1, size(foodfields)

       if (fgroup%food_sets(1)%path_integrate .and. allocated(agent%ele_path)) then

          ! Integrate along the path of the agent
          allocate(foodval_ele(size(agent%ele_path)))
          do e=1, size(agent%ele_path)
             foodval_ele(e) = agent%ele_dist(e) * integral_element(foodfields(f)%ptr, xfield, agent%ele_path(e)) / element_volume(xfield, agent%ele_path(e))
          end do

          foodfield_vals(f) = sum(foodval_ele)
          path_total = sum(agent%ele_dist)
          if (path_total > 0.0) then
             foodfield_vals(f) = foodfield_vals(f) / path_total
          else
             foodfield_vals(f) = integral_element(foodfields(f)%ptr, xfield, agent%element) / element_volume(xfield, agent%element)
          end if
          deallocate(foodval_ele)

       else
          foodfield_vals(f) = integral_element(foodfields(f)%ptr, xfield, agent%element) / element_volume(xfield, agent%element)
       end if
    end do

    stat=0
    call lebiology_agent_update(trim(fgroup%name), len_trim(fgroup%name), &
             trim(key), len_trim(key), trim(foodname), len_trim(foodname), &
             agent%biology, size(agent%biology), envfield_vals, size(envfield_vals), &
             foodfield_vals, agent%food_requests, agent%food_thresholds, agent%food_ingests, &
             size(foodfield_vals), dt, stat)

    do v=1, size(agent%biology)
       if (ieee_is_nan(agent%biology(v))) then
          ewrite(-1,*) "Error updating agent ", agent%id_number, ", agent state was:"
          do f=1, size(agent_state_copy)
             ewrite(-1,*) trim(fgroup%variables(f)%name), " --> ", agent_state_copy(f)
          end do
          FLExit('NaN agent variable detected in '//trim(fgroup%name)//"::"//trim(key))
       end if
    end do

    do f=1, size(foodfields)
       if (ieee_is_nan(agent%food_requests(f))) then
          FLExit('NaN agent food request detected in '//trim(fgroup%name)//"::"//trim(key))
       end if
       if (ieee_is_nan(agent%food_thresholds(f))) then
          FLExit('NaN agent food thresholds detected in '//trim(fgroup%name)//"::"//trim(key))
       end if
       if (ieee_is_nan(agent%food_ingests(f))) then
          FLExit('NaN agent food ingest detected in '//trim(fgroup%name)//"::"//trim(key))
       end if
    end do

    if (stat < 0) then
       ewrite(-1, *) "Error updating agent "//int2str(agent%id_number)//" for FG::"//trim(fgroup%name)
       FLExit("Python error in LE-Biology")
    end if
  end subroutine lebiology_update_agent

  subroutine fl_add_agent(vars, n_vars, pos, n_pos) &
         bind(c, name='fl_add_agent_c')
    use :: iso_c_binding
    implicit none
    integer(c_int), intent(inout) :: n_vars, n_pos
    real(c_double), dimension(n_vars), intent(inout) :: vars
    real(c_double), dimension(n_pos), intent(inout) :: pos

    type(detector_type), pointer :: agent

    allocate(agent)
    allocate(agent%position(n_pos))
    allocate(agent%biology(n_vars))
    allocate(agent%local_coords(n_pos+1))

    call get_next_detector_id(agent%id_number)
    agent%name = trim(int2str(agent%id_number))
    agent%type = LAGRANGIAN_DETECTOR

    agent%position = pos
    agent%biology = vars

    call insert(agent, new_agent_list)
  end subroutine fl_add_agent

end module lebiology_python
