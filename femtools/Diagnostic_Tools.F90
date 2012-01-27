!    Copyright (C) 2006 Imperial College London and others.   
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
#include "version.h"

module diagnostic_tools
  use futils
  use global_parameters, only: integer_size, real_size
  use c_interfaces

  implicit none

  private

  public :: field_tag, constant_tag, initialise_constant_diagnostics

contains

  function field_tag(name, column, statistic, material_phase_name, components)
    !!< Create a field tag for the given entries.
    character(len=*), intent(in) :: name
    integer, intent(in) :: column
    character(len=*), intent(in) :: statistic
    character(len=*), intent(in), optional :: material_phase_name 
    integer, intent(in), optional :: components
    character(len=254) :: field_tag

    character(len=254) :: front_buffer, material_buffer, components_buffer, end_buffer

    write(front_buffer,'(a,i0,a)') '<field column="',column,'" name="'&
            &//trim(name)//'" statistic="'//trim(statistic)//'"'

    if (present(material_phase_name)) then
        write(material_buffer,'(a)') ' material_phase="'//&
            trim(material_phase_name)//'"'
    else
        material_buffer = ''
    end if

    if (present(components)) then
        write(components_buffer,'(a,i0,a)') ' components="', components, '"' 
    else
        components_buffer = ''
    end if

    end_buffer = '/>'

    field_tag = trim(front_buffer)//trim(material_buffer)//trim(components_buffer)//trim(end_buffer)

  end function field_tag

  function constant_tag(name, type, value)
    !!< Create a field tag for the given entries.
    character(len=*), intent(in) :: name, type, value
    
    character(len=254) :: constant_tag
    
    constant_tag='<constant name="'//trim(name)&
         &//'" type="'//trim(type)//'" value="'//trim(value)//'" />'

  end function constant_tag

  subroutine initialise_constant_diagnostics(unit, binary_format)
    !!< Output constant values in the header of the stat file.
    integer, intent(in) :: unit
    !! If present and .true., indicates binary output format
    logical, optional, intent(in) :: binary_format
    
    character(len=254) :: buffer, value_buffer

#ifdef __FLUIDITY_VERSION__
    value_buffer = __FLUIDITY_VERSION__
#else
    value_buffer="Unknown"
#endif
    buffer=constant_tag(name="FluidityVersion", type="string", value=trim(value_buffer))
    write(unit, '(a)') trim(buffer)
    
    value_buffer = __DATE__ // " " // __TIME__
    buffer=constant_tag(name="CompileTime", type="string", value=trim(value_buffer))
    write(unit, '(a)') trim(buffer)

    value_buffer=date_and_time_string()
    buffer=constant_tag(name="StartTime", type="string", value=trim(value_buffer))
    write(unit, '(a)') trim(buffer)
    
    call get_environment_variable("HOSTNAME", value_buffer, default = "Unknown")
    buffer=constant_tag(name="HostName", type="string", value=trim(value_buffer))
    write(unit, '(a)') trim(buffer)
    
    ! Constant values
    if(present_and_true(binary_format)) then
      buffer = constant_tag(name = "format", type = "string", value = "binary")
      write(unit, '(a)') trim(buffer)
      buffer = constant_tag(name = "real_size", type = "integer", value = int2str(real_size))
      write(unit, '(a)') trim(buffer)
      buffer = constant_tag(name = "integer_size", type = "integer", value = int2str(integer_size))
      write(unit, '(a)') trim(buffer)
    else
      buffer = constant_tag(name = "format", type = "string", value = "plain_text")
      write(unit, '(a)') trim(buffer)
    end if
    
  contains 
    
    function date_and_time_string() 
      character(len=254) :: date_and_time_string

      character(len=8) :: date
      character(len=10) :: time
      character(len=5) :: zone

      call date_and_time(date, time, zone)
      
      date_and_time_string=date//" "//time//zone

    end function date_and_time_string

  end subroutine initialise_constant_diagnostics

end module diagnostic_tools
