#include "fdebug.h"

module lebiology_python
  use fldebug
  use fields
  use detector_data_types
  use detector_tools
  use detector_parallel
  use ieee_arithmetic, only: ieee_is_nan
  use iso_c_binding
  use element_path_list, only: elepath_list => linked_list, &
                               elepath_list_next => list_next

  implicit none
  
  private
  
  public :: lebiology_init_module, lebiology_set_stage_id, &
            lebiology_load_kernel, lebiology_reload_persistent, &
            lebiology_add_variables, lebiology_add_envfields, &
            lebiology_add_foods, lebiology_prepare_pyfunc, &
            lebiology_initialise_agent, lebiology_move_agent, &
            lebiology_update_agent, get_new_agent_list, &
            lebiology_pre_sample_environment

  type(detector_linked_list), target, save :: new_agent_list

  interface

    subroutine lebiology_init_module(dim) bind(c, name='initlebiology')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: dim
    end subroutine lebiology_init_module

    subroutine lebiology_add_fg_varname(fg, var, stat) &
           bind(c, name='lebiology_add_fg_varname_c')
      use :: iso_c_binding
      implicit none
      character(kind=c_char), dimension(*), intent(in) :: fg
      character(kind=c_char), dimension(*), intent(in) :: var
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_add_fg_varname

    subroutine lebiology_add_fg_envname(fg, env, stat) &
           bind(c, name='lebiology_add_fg_envname_c')
      use :: iso_c_binding
      implicit none
      character(kind=c_char), dimension(*), intent(in) :: fg
      character(kind=c_char), dimension(*), intent(in) :: env
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_add_fg_envname

    subroutine lebiology_add_fg_foodname(fg, food, variety, stat) &
           bind(c, name='lebiology_add_fg_foodname_c')
      use :: iso_c_binding
      implicit none
      character(kind=c_char), dimension(*), intent(in) :: fg
      character(kind=c_char), dimension(*), intent(in) :: food
      character(kind=c_char), dimension(*), intent(in) :: variety
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_add_fg_foodname

    subroutine lebiology_add_fg_stage_id(fg, stage, id, stat) &
           bind(c, name='lebiology_add_fg_stage_id_c')
      use :: iso_c_binding
      implicit none
      character(kind=c_char), dimension(*), intent(in) :: fg
      character(kind=c_char), dimension(*), intent(in) :: stage
      real(c_double), intent(in) :: id
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_add_fg_stage_id

    subroutine lebiology_fg_kernel_load(fg, key, module, kernel, param, stat) &
         bind(c, name='lebiology_fg_kernel_load_c')
      use :: iso_c_binding
      implicit none
      character(kind=c_char), dimension(*), intent(in) :: fg
      character(kind=c_char), dimension(*), intent(in) :: key
      character(kind=c_char), dimension(*), intent(in) :: module
      character(kind=c_char), dimension(*), intent(in) :: kernel
      character(kind=c_char), dimension(*), intent(in) :: param
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_fg_kernel_load

    subroutine lebiology_reload_persistent() &
         bind(c, name='lebiology_reload_persistent_c')
      use :: iso_c_binding
      implicit none
    end subroutine lebiology_reload_persistent

    subroutine lebiology_compile_function(fg, key, func, stat) &
           bind(c, name='lebiology_compile_function_c')
      use :: iso_c_binding
      implicit none
      character(kind=c_char), dimension(*), intent(in) :: fg
      character(kind=c_char), dimension(*), intent(in) :: key
      character(kind=c_char), dimension(*), intent(in) :: func
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_compile_function

    subroutine lebiology_agent_init(fg, key, vars, n_vars, stat) &
           bind(c, name='lebiology_agent_init_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: n_vars
      character(kind=c_char), dimension(*), intent(in) :: fg
      character(kind=c_char), dimension(*), intent(in) :: key
      real(c_double), dimension(n_vars), intent(inout) :: vars
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_agent_init

    subroutine lebiology_parallel_prepare(fg, key, food, vars, n_vars, &
         envvals, n_envvals, fvariety, fingest, n_fvariety, &
         agent_id, dt, persistent, stat) &
         bind(c, name='lebiology_parallel_prepare_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: n_vars, n_envvals, n_fvariety
      character(kind=c_char), dimension(*), intent(in) :: fg
      character(kind=c_char), dimension(*), intent(in) :: key
      character(kind=c_char), dimension(*), intent(in) :: food
      real(c_double), dimension(n_vars), intent(inout) :: vars
      real(c_double), dimension(n_envvals), intent(inout) :: envvals
      real(c_double), dimension(n_fvariety), intent(inout) :: fvariety, fingest
      real(c_double), intent(in) :: dt
      integer(c_int), intent(in), value :: agent_id, persistent
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_parallel_prepare

    subroutine lebiology_parallel_finish(fg, food, vars, n_vars, &
         frequest, fthreshold, n_fvariety, agent_id, stat) &
         bind(c, name='lebiology_parallel_finish_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: n_vars, n_fvariety
      character(kind=c_char), dimension(*), intent(in) :: fg
      character(kind=c_char), dimension(*), intent(in) :: food
      real(c_double), dimension(n_vars), intent(inout) :: vars
      real(c_double), dimension(n_fvariety), intent(inout) :: frequest, fthreshold
      integer(c_int), intent(in), value :: agent_id
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_parallel_finish

    subroutine lebiology_kernel_update(fg, key, food, vars, n_vars, &
         envvals, n_envvals, fvariety, frequest, fthreshold, fingest, &
         n_fvariety, dt, persistent, stat) &
         bind(c, name='lebiology_kernel_update_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: n_vars, n_envvals, n_fvariety
      character(kind=c_char), dimension(*), intent(in) :: fg
      character(kind=c_char), dimension(*), intent(in) :: key
      character(kind=c_char), dimension(*), intent(in) :: food
      real(c_double), dimension(n_vars), intent(inout) :: vars
      real(c_double), dimension(n_envvals), intent(inout) :: envvals
      real(c_double), dimension(n_fvariety), intent(inout) :: fvariety, frequest, fingest, fthreshold
      real(c_double), intent(in) :: dt
      integer(c_int), intent(in), value :: persistent
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_kernel_update

    subroutine lebiology_agent_update(fg, key, food, vars, n_vars, envvals, n_envvals, &
         fvariety, frequest, fthreshold, fingest, n_fvariety, dt, dropout, stat) &
         bind(c, name='lebiology_agent_update_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: n_vars, n_envvals, n_fvariety
      character(kind=c_char), dimension(*), intent(in) :: fg
      character(kind=c_char), dimension(*), intent(in) :: key
      character(kind=c_char), dimension(*), intent(in) :: food
      real(c_double), dimension(n_vars), intent(inout) :: vars
      real(c_double), dimension(n_envvals), intent(inout) :: envvals
      real(c_double), dimension(n_fvariety), intent(inout) :: fvariety, frequest, fingest, fthreshold
      real(c_double), intent(in) :: dt
      integer(c_int), intent(out) :: dropout, stat
    end subroutine lebiology_agent_update

    subroutine lebiology_agent_move(fg, key, pos, n_pos, &
           vars, n_vars, var_inds, dt, vector, stat) &
           bind(c, name='lebiology_agent_move_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value ::n_pos, n_vars
      character(kind=c_char), dimension(*), intent(in) :: fg
      character(kind=c_char), dimension(*), intent(in) :: key
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
       call lebiology_add_fg_varname(trim(fgroup%name)//C_NULL_CHAR, &
             trim(fgroup%variables(v)%name)//C_NULL_CHAR, stat)
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
       call lebiology_add_fg_envname(trim(fgroup%name)//C_NULL_CHAR, &
             trim(fgroup%envfield_names(f))//C_NULL_CHAR, stat)
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
          call lebiology_add_fg_foodname(trim(fgroup%name)//C_NULL_CHAR, &
                 trim(fset%name)//C_NULL_CHAR, trim(fset%varieties(s)%name)//C_NULL_CHAR, stat)
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
    call lebiology_add_fg_stage_id(trim(fgroup%name)//C_NULL_CHAR, &
             trim(stage_name)//C_NULL_CHAR, id, stat)

    if (stat < 0) then
       ewrite(-1, *) "Error setting stage ID for FG::"//trim(fgroup%name)
       FLExit("Python error in LE-Biology")
    end if
  end subroutine lebiology_set_stage_id

  subroutine lebiology_load_kernel(fgroup, key, module, kernel, paramset)
    type(functional_group), intent(inout) :: fgroup
    character(len=*), intent(in) :: key, module, kernel, paramset
    integer :: stat

    stat=0
    call lebiology_fg_kernel_load(trim(fgroup%name)//C_NULL_CHAR, &
             trim(key)//C_NULL_CHAR, trim(module)//C_NULL_CHAR, &
             trim(kernel)//C_NULL_CHAR, trim(paramset)//C_NULL_CHAR, stat)

    if (stat < 0) then
       ewrite(-1, *) "Error loading kernel for FG::"//trim(fgroup%name)
       ewrite(-1, *) "Kernel name: "//trim(kernel)
       FLExit("Python error in LE-Biology")
    end if
  end subroutine lebiology_load_kernel

  subroutine lebiology_prepare_pyfunc(fgroup, key, func)
    type(functional_group), intent(inout) :: fgroup
    character(len=*), intent(in) :: key, func
    integer :: stat

    stat=0
    call lebiology_compile_function(trim(fgroup%name)//C_NULL_CHAR, &
             trim(key)//C_NULL_CHAR, trim(func)//C_NULL_CHAR, stat)

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
    call lebiology_agent_init(trim(fgroup%name)//C_NULL_CHAR, &
             trim(key)//C_NULL_CHAR, agent%biology, size(agent%biology), stat)

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
    call lebiology_agent_move(trim(fgroup%name)//C_NULL_CHAR, trim(key)//C_NULL_CHAR, &
            agent%update_vector, size(agent%update_vector), &
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

  subroutine lebiology_pre_sample_environment(agent, fields, sampling)
    type(detector_type), intent(inout) :: agent
    type(scalar_field_pointer), dimension(:), pointer, intent(inout) :: fields
    ! Sampling method to use for each field
    integer, dimension(size(fields)), intent(inout) :: sampling

    integer :: f, i_env

    ! Sample environment fields
    i_env = 1
    do f=1, size(fields)
       if (sampling(f) == VARSMPL_OLDZ) then
          agent%env_samples(i_env) = eval_field(agent%element, fields(f)%ptr, agent%local_coords)
          i_env = i_env + 1
       end if
    end do
  end subroutine lebiology_pre_sample_environment

  subroutine lebiology_sample_environment(agent, xfield, fields, sampling, fieldvals)
    type(detector_type), intent(inout) :: agent
    type(vector_field), pointer, intent(inout) :: xfield
    type(scalar_field_pointer), dimension(:), pointer, intent(inout) :: fields
    ! Sampling method to use for each field
    integer, dimension(size(fields)), intent(inout) :: sampling
    real, dimension(size(fields)), intent(inout) :: fieldvals

    type(elepath_list), pointer :: path_ele
    real :: path_total, ele_integral, ele_volume
    integer :: f, i_env

    ! Sample environment fields
    i_env = 1
    do f=1, size(fields)
       if (sampling(f) == VARSMPL_INTZ .and. associated(agent%path_elements)) then
          fieldvals(f) = 0.0
          path_total = 0.0
          path_ele => agent%path_elements
          do while( associated(path_ele) )
              ele_integral = integral_element(fields(f)%ptr, xfield, path_ele%data%ele)
              ele_volume = element_volume(xfield, path_ele%data%ele)
              fieldvals(f) = fieldvals(f) + path_ele%data%dist * ele_integral / ele_volume
              path_total = path_total + path_ele%data%dist

              path_ele => elepath_list_next(path_ele)
          end do

          if (path_total > 0.0) then
             fieldvals(f) = fieldvals(f) / path_total
          else
             fieldvals(f) = integral_element(fields(f)%ptr, xfield, agent%element) / element_volume(xfield, agent%element)
          end if
       elseif (sampling(f) == VARSMPL_OLDZ) then
          ! Use buffered field value from beginning of timestep
          fieldvals(f) = agent%env_samples(i_env)
          i_env = i_env + 1
       elseif (sampling(f) == VARSMPL_NEWZ .or. &
            sampling(f) == VARSMPL_INTZ .and. .not. associated(agent%path_elements)) then
          ! Evaluate at current agent position
          fieldvals(f) = eval_field(agent%element, fields(f)%ptr, agent%local_coords)
       else
          FLExit("Sampling error: Have you used the right tracking scheme?")
       end if
    end do
  end subroutine lebiology_sample_environment

  subroutine lebiology_update_agent(fgroup, key, foodname, agent, xfield, envfields, &
       foodfields, dt, use_kernel_func, use_persistent)
    type(functional_group), intent(inout) :: fgroup
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: foodname
    type(detector_type), intent(inout) :: agent
    type(vector_field), pointer, intent(inout) :: xfield
    type(scalar_field_pointer), dimension(:), pointer, intent(inout) :: envfields
    type(scalar_field_pointer), dimension(:), pointer, intent(inout) :: foodfields
    real, intent(in) :: dt
    logical, intent(in) :: use_kernel_func, use_persistent

    real, dimension(size(envfields)) :: envfield_vals
    real, dimension(size(foodfields)) :: foodfield_vals
    integer, dimension(size(foodfields)) :: food_sampling
    real :: path_total, ele_integral, ele_volume
    real, dimension(size(agent%biology)) :: agent_state_copy
    type(elepath_list), pointer :: path_ele
    integer :: f, v, e, persistent, dropout_agent, stat

    agent_state_copy = agent%biology

    call lebiology_sample_environment(agent, xfield, envfields, fgroup%envfield_sampling, envfield_vals)

    if (size(fgroup%food_sets) > 0) then
       food_sampling(:) = fgroup%food_sets(1)%sampling
       call lebiology_sample_environment(agent, xfield, foodfields, food_sampling, foodfield_vals)
    end if

    stat=0
    if (use_kernel_func) then                   
       if (use_persistent) then
          persistent = 1
       else
          persistent = 0
       end if
       call lebiology_kernel_update(trim(fgroup%name)//C_NULL_CHAR, &
             trim(key)//C_NULL_CHAR, trim(foodname)//C_NULL_CHAR, &
             agent%biology, size(agent%biology), envfield_vals, size(envfield_vals), &
             foodfield_vals, agent%food_requests, agent%food_thresholds, agent%food_ingests, &
             size(foodfield_vals), dt, persistent, stat)
    else
       call lebiology_agent_update(trim(fgroup%name)//C_NULL_CHAR, &
             trim(key)//C_NULL_CHAR, trim(foodname)//C_NULL_CHAR, &
             agent%biology, size(agent%biology), envfield_vals, size(envfield_vals), &
             foodfield_vals, agent%food_requests, agent%food_thresholds, agent%food_ingests, &
             size(foodfield_vals), dt, dropout_agent, stat)
    end if

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

    if (dropout_agent > 0) then
       agent%dropout = .true.
    end if

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

    agent%path_elements => null()

    call insert(agent, new_agent_list)
  end subroutine fl_add_agent

end module lebiology_python
