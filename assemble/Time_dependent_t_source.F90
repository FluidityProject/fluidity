module time_dependent_t_source
#include "fdebug.h"
  use fldebug
  use Futils
  implicit none

  logical :: initialised = .false.
  logical :: have_time_dependent_temperature_source = .false.
  real :: time_period,amplification_factor, pulse_length
contains  

  subroutine get_time_dependent_temperature_source()
    integer :: unit, io
    namelist/temp_source/have_time_dependent_temperature_source, time_period, &
         & amplification_factor, pulse_length

    if(.not.initialised) then
       ewrite(3,*) 'Reading in file'
       unit=free_unit()
       open(unit=unit, file="time_dependent_temperature_source.dat", action="read",&
            & status="old", iostat=io)
       if(io == 0) then
          ewrite(3,*) 'Reading data'
          read(unit, nml=temp_source)
          initialised = .true.
       end if
    end if

    ewrite(3,*) 'have_time_dependent_temperature_source', &
         have_time_dependent_temperature_source
    ewrite(3,*) 'time_period',time_period
    ewrite(3,*) 'amplification_factor', amplification_factor
    ewrite(3,*) 'pulse_length', pulse_length

  end subroutine get_time_dependent_temperature_source

  function get_time_dependent_source_factor(T) result (factor)
    real, intent(in) :: T
    real :: factor

    !locals
    real :: local_T 

    if(have_time_dependent_temperature_source) then
    
       local_T = T - time_period*floor(T/time_period)
       ewrite(3,*) 'local_T, T', local_T, T
       if(local_T<pulse_length) then
          factor = amplification_factor
       else
          factor = 0.0
       end if
    else
       factor = 1.0
    end if
    ewrite(3,*) 'factor', factor

  end function get_time_dependent_source_factor

end module time_dependent_t_source
