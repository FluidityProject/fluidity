
subroutine test_turbine

 use fldebug
 use spud
 use global_parameters, only:pi
 use futils

 use actuator_line_model
 use actuator_line_source
    
 implicit none

 integer:: istep, Ntimesteps,itur,iblade,nargin,FNLength,status,dump_step
 real   :: current_time, final_time, dt
 character(80) :: InputFN
 
 write(*,*) '================================================='
 write(*,*) 'This is a reduced model for the ALM '
 write(*,*) '================================================='

 call load_options("Test.flml")
 
 Ntimesteps=1000
 deltaT=0.001
 visc=1.5e-5
 current_time=0.0
 dump_step=1
 
 call actuator_line_model_init
 
 write(6,*) 'Number of Turbines :', Ntur
 write(6,*) 'Number of Blades :', Turbine(1:Ntur)%NBlades
    

 write(6,*) 'Starting the Time-loop :' 
 do istep=1,Ntimesteps

 !> Set the velocities
    do itur=1, Ntur
        do iblade=1,Turbine(itur)%NBlades
        Turbine(itur)%Blade(iblade)%EVx(:)=10.0
        Turbine(itur)%Blade(iblade)%EVy(:)=0.0
        Turbine(itur)%Blade(iblade)%EVz(:)=0.0
        enddo
    end do

    call actuator_line_model_compute_forces    
    
    current_time=current_time+deltaT

    call actuator_line_model_update(current_time,deltaT)
  
    if (mod(istep,dump_step)==0) then
    write(6,*) 'Writing output for time step : ', istep*deltaT 
    call actuator_line_model_write_output(istep)
    endif

end do

end subroutine test_turbine

