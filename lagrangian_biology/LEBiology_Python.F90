#include "fdebug.h"

module lebiology_python
  use fldebug
  use fields
  use detector_data_types
  use Profiler
  use ieee_arithmetic, only: ieee_is_nan

  implicit none
  
  private
  
  public :: lebiology_init_module, lebiology_set_stage_id, &
            lebiology_add_variables, lebiology_add_envfields, &
            lebiology_prepare_pyfunc, lebiology_initialise_agent, &
            lebiology_update_agent

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

    subroutine lebiology_agent_update(fg, fglen, key, keylen, vars, n_vars, &
           envvals, n_envvals, dt, stat) &
           bind(c, name='lebiology_agent_update_c')
      use :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: fglen, keylen, n_vars, n_envvals
      character(kind=c_char), dimension(fglen), intent(in) :: fg
      character(kind=c_char), dimension(keylen), intent(in) :: key
      real(c_double), dimension(n_vars), intent(inout) :: vars
      real(c_double), dimension(n_envvals), intent(inout) :: envvals
      real(c_double), intent(in) :: dt
      integer(c_int), intent(out) :: stat
    end subroutine lebiology_agent_update

  end interface

contains

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

  subroutine lebiology_update_agent(fgroup, key, agent, envfields, dt)
    type(functional_group), intent(inout) :: fgroup
    character(len=*), intent(in) :: key
    type(detector_type), intent(inout) :: agent
    type(scalar_field_pointer), dimension(:), pointer, intent(inout) :: envfields
    real, intent(in) :: dt

    real, dimension(size(envfields)) :: envfield_vals
    integer :: f, v, stat

    call profiler_tic(trim(fgroup%name)//"::"//trim(key))

    ! Sample environment fields
    do f=1, size(envfields)
       envfield_vals(f) = eval_field(agent%element, envfields(f)%ptr, agent%local_coords)
    end do

    stat=0
    call lebiology_agent_update(trim(fgroup%name), len_trim(fgroup%name), &
             trim(key), len_trim(key), agent%biology, size(agent%biology), &
             envfield_vals, size(envfield_vals), dt, stat)

    do v=1, size(agent%biology)
       if (ieee_is_nan(agent%biology(v))) then
          FLExit('NaN agent variable detected in '//trim(fgroup%name)//"::"//trim(key))
       end if
    end do

    call profiler_toc(trim(fgroup%name)//"::"//trim(key))

    if (stat < 0) then
       ewrite(-1, *) "Error updating agent "//int2str(agent%id_number)//" for FG::"//trim(fgroup%name)
       FLExit("Python error in LE-Biology")
    end if
  end subroutine lebiology_update_agent

end module lebiology_python
