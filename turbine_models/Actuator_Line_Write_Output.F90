#include"fdebug.h"

module actuator_line_write_output
    
    use fldebug
    use spud
    use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, PYTHON_FUNC_LEN, pi
    use futils

    use actuator_line_turbine
    use actuator_line_element

contains


    subroutine actuator_line_element_write_output(act_line,dir)

        implicit none
        type(ActuatorLineType),intent(in) :: act_line
        character(len=100),intent(in) :: dir
        character(LEN=22) :: Format
        integer :: ielem
        open(2017,File=trim(dir)//'/'//trim(act_line%name)//'.load')
        write(2017,*) 'ielem,X,Y,Z,r/R,pitch,AOA,adot,RE,ur,eps,CL,CD,CM25,Cn,Ct,Fn,Ft,tipFactor'
        Format="(I5,A,18(E14.7,A))"
        do ielem=1,act_line%NElem
write(2017,Format)ielem,',',act_line%PEx(ielem),',',act_line%PEy(ielem),',',act_line%PEz(ielem),',',act_line%ERdist(ielem)/act_line%L,',',act_line%Epitch(ielem)*180/pi,',',act_line%EAOA(ielem)*180/pi,',',act_line%EAOAdot(ielem)*pi/180,',',act_line%ERE(ielem),',',act_line%EUr(ielem),',',act_line%Eepsilon(ielem),',',act_line%ECL(ielem),',',act_line%ECD(ielem),',',act_line%ECM(ielem),',',act_line%ECN(ielem),',',act_line%ECT(ielem),',',act_line%EFn(ielem),',',act_line%EFt(ielem),',',act_line%EEndeffects_factor(ielem)
        end do
        close(2017)
 
    end subroutine actuator_line_element_write_output
    
    subroutine dynamic_stall_write_output(act_line,dir)

        implicit none
        type(ActuatorLineType),intent(in) :: act_line
        character(len=100),intent(in) :: dir
        character(LEN=22) :: Format
        integer :: ielem
        open(2018,File=trim(dir)//'/'//trim(act_line%name)//'.dynstall')
        write(2018,*) 'ielem,rdist/R,pitch,AOA,f'
        Format="(I5,A,15(E14.7,A))"
        do ielem=1,act_line%NElem
write(2018,Format)ielem,',',act_line%ERdist(ielem)/act_line%L,',',act_line%Epitch(ielem)*180/pi,',',act_line%EAOA(ielem)*180/pi,',',act_line%EDynstall(ielem)%F
        end do
        close(2018)
 
    end subroutine dynamic_stall_write_output
    
    subroutine actuator_line_turbine_write_output(turbine,dir)

        implicit none
        type(TurbineType),intent(in) :: turbine
        character(len=100),intent(in) :: dir
        character(LEN=22) :: Format
        
        open(2016,File=trim(dir)//'/'//trim(turbine%name)//'.perf')
        write(2016,*) 'Number of Revs, Urotor_averaged, CFx , CFy , CFz , CT , CTR , CP'
        Format="(8(E14.7,A))"
write(2016,Format) turbine%AzimAngle/(2*pi),',',turbine%Urotor_average,',' ,turbine%CFx,',',turbine%CFy,',',turbine%CFz,',',turbine%CT,',',turbine%CTR,',',turbine%CP
        close(2016)

    end subroutine actuator_line_turbine_write_output

end module actuator_line_write_output
