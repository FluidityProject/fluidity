#include "fdebug.h"

module momentum_source
  !!< This module adds a momentum source specified by a python code.
  use state_module
  use fields
  use global_parameters, only : PYTHON_FUNC_LEN
  use spud
  implicit none
  
  private
  
  public :: SIMPLE_SOURCE_python

contains

  subroutine SIMPLE_SOURCE_python(state)
    use FLDebug
    use AllSorts
    IMPLICIT NONE
    type(state_type), intent(inout) :: state

    type(vector_field) :: position, source, lsource
    real :: time
    logical, save :: first_time=.true., file_exists
    character(len=PYTHON_FUNC_LEN), save :: func
    character(len=PYTHON_FUNC_LEN) :: buffer
    integer :: unit

    if (first_time) then
       first_time=.false.
       inquire(file="velocitysource.py",exist=file_exists)
    
       if (.not.file_exists) return

       unit=free_unit()
       open(unit, file="velocitysource.py", action="read",&
            & status="old")
       read(unit, '(a)', end=42) func
       
       ! Read all the lines of the file and put in newlines between them.
       do
          read(unit, '(a)', end=42) buffer
         
          func=trim(func)//achar(10)//trim(buffer)
          
       end do
       
42     func=trim(func)
       close(unit)

    else
       if (.not.file_exists) return
    end if

    source=extract_vector_field(state, "VelocitySource")
    position=extract_vector_field(state, "Coordinate")
    
    call get_option("/timestepping/current_time", time)

    call allocate(lsource, source%dim, source%mesh)
    
    call set_from_python_function(lsource, func, position, time)
    
    call addto(source, lsource)

    call deallocate(lsource)

    ewrite(2,*) 'end Subroutine simple_source'

  END SUBROUTINE SIMPLE_SOURCE_python

end module momentum_source
