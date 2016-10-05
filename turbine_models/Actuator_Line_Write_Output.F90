#include "fdebug.h"

module actuator_line_write_output
    
    use fldebug
    use spud
    use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, PYTHON_FUNC_LEN, pi
    use futils

    use actuator_line_turbine
    use actuator_line_element

contains

    subroutine actuator_line_element_write_output(act_line,dump_num)

        implicit none
        type(ActuatorLineType),intent(in) :: act_line
        integer,intent(in) :: dump_num
        character(LEN=40) :: Format
        integer :: ielem


        open(2017,File=trim(act_line%name)//'_'//int2str(dump_num)//'.dat')
        write(2017,*) 'ielem, AOA, CD, CDL,CM'
        Format="(I5,A,F10.5,A,F10.5,A,F10.5,A,F10.5)"
        do ielem=1,act_line%NElem
        write(2017,Format) ielem,',',act_line%EAOA_Last(ielem),',',act_line%ECD(ielem),',',act_line%ECL(ielem),',',act_line%ECM(ielem)
        end do
        close(2017)


 
    end subroutine actuator_line_element_write_output

end module actuator_line_write_output
