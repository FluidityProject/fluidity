subroutine test_actuator_line


 use fldebug
 use spud
 use global_parameters, only:pi
 use futils

 use actuator_line_model
 use actuator_line_source
    
 implicit none

 integer:: istep, Ntimesteps,ial
 real   :: current_time, final_time, dt

 write(*,*) '================================================='
 write(*,*) 'This is a reduced model for uALM/uBEM model that '
 write(*,*) 'that works without a dynamic inflow model        '
 write(*,*) '================================================='

 call load_options("test_actuator_line_model.flml")
 
 Ntimesteps=100
 deltaT=0.01
 visc=1.5e-5
 current_time=0.0
 
 call actuator_line_model_init
 
 write(*,*) 'Number of Actuator lines :', Nal
 write(*,*) 'Number of Elements of the AL :', Actuatorline(1:Nal)%Nelem

 do istep=1,Ntimesteps

 !> Set the velocities
    do ial=1, Nal
        Actuatorline(ial)%EVx(:)=1.0
        Actuatorline(ial)%EVy(:)=0.0
        Actuatorline(ial)%EVz(:)=0.0
    end do
    
    call actuator_line_model_compute_forces    

    current_time=current_time+deltaT

    call actuator_line_model_update(current_time,deltaT)
   
    do ial=1, Nal
    write(*,*) Actuatorline(ial)%EAOA_Last(11)
    call actuator_line_model_write_output(istep)
    end do

end do

end subroutine test_actuator_line
