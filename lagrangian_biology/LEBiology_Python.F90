#include "fdebug.h"

module lebiology_python
  use fldebug
  use detector_data_types

  implicit none
  
  private
  
  public :: lebiology_init_module, lebiology_add_variables, &
            lebiology_set_stage_id

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

end module lebiology_python
