#include "fdebug.h"

module actuator_line_model

    use fldebug
    use spud
    use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, PYTHON_FUNC_LEN, pi
    use futils

    !use the actuator_line Modules
    use airfoils
    use actuator_line_model_utils
    use actuator_line_element
    use actuator_line_turbine
    use actuator_line_write_output
    use dynstall

    implicit none

    type(ActuatorLineType), allocatable, save :: Actuatorline(:)
    type(TurbineType), allocatable, save :: Turbine(:) ! Turbine 
    integer,save :: Ntur, Nal ! Number of the turbines 
    real,save :: deltaT, Visc, ctime
    logical,save :: actuator_line_model_writeFlag=.true.

    public  actuator_line_model_init, actuator_line_model_compute_forces, actuator_line_model_update 

contains

    subroutine actuator_line_model_init

        implicit none
        integer :: itur,ial 

        ewrite(1,*) 'Entering the actuator_line_model_init '

        !### Specify Turbines
        call get_turbine_options

        if (Ntur>0) then 
            do itur=1,Ntur
            call set_turbine_geometry(Turbine(itur))
            end do
        endif

        !### Speficy Actuator Lines

        call get_actuatorline_options 
        if(Nal>0) then
            do ial=1,Nal
            call set_actuatorline_geometry(actuatorline(ial))
            end do
        endif

        ewrite(1,*) 'Exiting the actuator_line_model_init '

    end subroutine actuator_line_model_init

    subroutine actuator_line_model_write_output(dump_no)

        implicit none
        integer,intent(in) :: dump_no
        integer :: itur,ial,iblade
        character(len=100) :: dir

        call system('mkdir -p ALM/'//adjustl(trim(dirname(dump_no))))
        
        dir='ALM/'//adjustl(trim(dirname(dump_no)))

        if (Ntur>0) then
            do itur=1,Ntur
                call actuator_line_turbine_write_output(turbine(itur),dir)
                do iblade=1,turbine(itur)%NBlades
                 call actuator_line_element_write_output(turbine(itur)%Blade(iblade),dir)
                end do
                if(turbine(itur)%Has_Tower) then 
                call actuator_line_element_write_output(turbine(itur)%tower,dir)
                endif
            end do
        endif

        if (Nal>0) then
            do ial=1,Nal
            call actuator_line_element_write_output(actuatorline(ial),dir)

            if (actuatorline(ial)%do_dynamic_stall) then
               call dynamic_stall_write_output(actuatorline(ial),dir) 
            end if
       
            end do
        endif
    
    end subroutine actuator_line_model_write_output

    subroutine get_turbine_options

        implicit none

        character(len=OPTION_PATH_LEN)::  turbine_name, actuatorline_name
        integer :: i,j,k
        integer, parameter :: MaxReadLine = 1000    
        character(MaxReadLine) :: FN    ! path to geometry input file 
        integer :: NElem, nfoils
        character(MaxReadLine) :: ReadLine
        character(len=OPTION_PATH_LEN) :: section_path
        character(len=OPTION_PATH_LEN), allocatable :: turbine_path(:), actuatorline_path(:)

        ewrite(2,*) 'Entering get_turbine_options'

        Ntur = option_count("/turbine_models/actuator_line_model/turbine")
        ewrite(2,*) 'Number of Turbines : ', Ntur

        ! Allocate Turbines Arrays
        Allocate(Turbine(Ntur))
        Allocate(turbine_path(Ntur))

        ! ==========================================
        ! Get Turbines' options and INITIALIZE THEM
        ! ==========================================
        do i=1, Ntur 

        turbine_path(i)="/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]"
        call get_option("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/name",Turbine(i)%name)
        Turbine(i)%ID=i    
        !###########1 Blade Specs #############################################   
        call get_option("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/location/",Turbine(i)%origin) 
        call get_option("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/Blades/number_of_blades/",Turbine(i)%NBlades)
        call get_option("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/Blades/blade_geometry/file_name",Turbine(i)%blade_geom_file)

        ! Allocate Blades
        Allocate(Turbine(i)%Blade(Turbine(i)%NBlades))
        ! Count how many Airfoil Sections are available
        nfoils = option_count("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/Blades/static_foil_data/foil")
        if(nfoils==0) then
            FLExit("You need to provide at least on static_foils_data entry for the computation of the blade forces")
        end if
        ewrite(2,*) 'Number of Static Foil Data available  for the analysis of the blades: ', nfoils
        ! Allocate the memory of the Airfoils

        do j=1,Turbine(i)%NBlades
        Turbine(i)%Blade(j)%NAirfoilData=nfoils
        Allocate(Turbine(i)%Blade(j)%AirfoilData(nfoils))
        do k=1, Turbine(i)%Blade(j)%NAirfoilData
        call get_option(trim(turbine_path(i))//"/Blades/static_foil_data/foil["//int2str(k-1)//"]/foil_file",Turbine(i)%Blade(j)%AirfoilData(k)%afname)
        ! Read and Store Airfoils
        call airfoil_init_data(Turbine(i)%Blade(j)%AirfoilData(k))
        end do
        end do

        ! ## Tower ?
        if (have_option("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/tower")) then
            Turbine(i)%Has_Tower=.true.
            call get_option("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/tower/offset",Turbine(i)%TowerOffset)
            call get_option("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/tower/tower_geometry/file_name",Turbine(i)%Tower%geom_file)
            call get_option("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/tower/tower_loading/Drag_coeff",Turbine(i)%TowerDrag)
            call get_option("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/tower/tower_loading/Lift_coeff",Turbine(i)%TowerLift)
            call get_option("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/tower/tower_loading/Strouhal_number",Turbine(i)%TowerStrouhal)
        
        endif

        !#############2  Get turbine_specs #################
        ! Check the typ of Turbine (choose between Horizontal and Vertical Axis turbines) 
        if(have_option(trim(turbine_path(i))//"/type/Horizontal_Axis")) then
            Turbine(i)%Type='Horizontal_Axis'
            call get_option(trim(turbine_path(i))//"/type/Horizontal_Axis/hub_tilt_angle",Turbine(i)%hub_tilt_angle)
            call get_option(trim(turbine_path(i))//"/type/Horizontal_Axis/blade_cone_angle",Turbine(i)%blade_cone_angle)
            call get_option(trim(turbine_path(i))//"/type/Horizontal_Axis/yaw_angle",Turbine(i)%yaw_angle)
        elseif(have_option(trim(turbine_path(i))//"/type/Vertical_Axis")) then
            FLExit("At the moment only the Horizontal_Axis Turbine is available")
        else
            FLExit("You should not be here")
        end if

        !##############3 Get Operation Options ######################
        if (have_option(trim(turbine_path(i))//"/operation/constant_rotational_velocity")) then
            Turbine(i)%Is_constant_rotation_operated= .true.
            call get_option("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/operation/constant_rotational_velocity/omega",Turbine(i)%angularVel)
            call get_option("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/operation/constant_rotational_velocity/TSR",Turbine(i)%TSR)
            if(have_option(trim(turbine_path(i))//"/operation/constant_rotational_velocity/rotation_direction/clockwise")) then
                Turbine(i)%IsClockwise=.true.
            elseif(have_option(trim(turbine_path(i))//"/operation/constant_rotational_velocity/rotation_direction/counter_clockwise")) then
                Turbine(i)%IsCounterClockwise=.true.
            else
                FLExit("You should not be here. The options are clockwise and counterclockwise")
            endif 
        else if(have_option(trim("turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]")//"/operation/force_based_rotational_velocity")) then
            Turbine(i)%Is_force_based_operated = .true. 
        else
            FLExit("At the moment only the constant and the force_based rotational velocity models are supported") 
        endif

        !##################4 Get Unsteady Modelling Options ##################
        if(have_option(trim(turbine_path(i))//"/Blades/unsteady_modelling/added_mass")) then
            do j=1,Turbine(i)%NBlades
            Turbine(i)%Blade(j)%do_added_mass=.true.
            end do
        endif

        if(have_option(trim(turbine_path(i))//"/Blades/unsteady_modelling/dynamic_stall")) then
            do j=1,Turbine(i)%NBlades
            Turbine(i)%Blade(j)%do_dynamic_stall=.true.
            end do
        endif

        if(have_option(trim(turbine_path(i))//"/Blades/unsteady_modelling/tip_loss_correction")) then
            Turbine(i)%do_tip_correction=.true.
        endif


        end do

        ewrite(2,*) 'Exiting get_turbine_options'

    end subroutine get_turbine_options 

    subroutine get_actuatorline_options

        implicit none

        character(len=OPTION_PATH_LEN)::  actuatorline_name
        integer :: i,j,k
        integer, parameter :: MaxReadLine = 1000    
        character(MaxReadLine) :: FN    ! path to geometry input file 
        integer :: NElem, iDoF
        character(MaxReadLine) :: ReadLine
        character(len=OPTION_PATH_LEN) :: section_path
        character(len=OPTION_PATH_LEN), allocatable :: actuatorline_path(:)

        ewrite(2,*) 'Entering get_actuatorline_options'

        Nal = option_count("/turbine_models/actuator_line_model/actuatorline")
        ewrite(2,*) 'Number of Actuatorlines : ', Nal

        ! Allocate Turbines Arrays
        allocate(actuatorline(Nal))
        Allocate(actuatorline_path(Nal))

        ! ==========================================
        ! Get Turbines' options and INITIALIZE THEM
        ! ==========================================
        do i=1, Nal
        actuatorline_path(i)="/turbine_models/actuator_line_model/actuatorline["//int2str(i-1)//"]"
        call get_option("/turbine_models/actuator_line_model/actuatorline["//int2str(i-1)//"]/name",Actuatorline(i)%name)
        call get_option("/turbine_models/actuator_line_model/actuatorline["//int2str(i-1)//"]/location/",Actuatorline(i)%COR) 
        call get_option("/turbine_models/actuator_line_model/actuatorline["//int2str(i-1)//"]/geometry_file/file_name",Actuatorline(i)%geom_file)

        ! Count how many Airfoil Sections are available
        Actuatorline(i)%NAirfoilData=option_count("/turbine_models/actuator_line_model/actuatorline["//int2str(i-1)//"]/airfoil_sections/section") 
        ewrite(2,*) 'Number of Airfoils available : ', Actuatorline(i)%NAirfoilData
        ! Allocate the memory of the Airfoils
        Allocate(Actuatorline(i)%AirfoilData(Actuatorline(i)%NAirfoilData))

        do k=1, Actuatorline(i)%NAirfoilData

        call get_option(trim(Actuatorline_path(i))//"/airfoil_sections/section["//int2str(k-1)//"]/airfoil_file",Actuatorline(i)%AirfoilData(k)%afname)

        ! Read and Store Airfoils
        call airfoil_init_data(Actuatorline(i)%AirfoilData(k))
        end do

        !##################4 Get Unsteady Modelling Options ##################
        if(have_option(trim(actuatorline_path(i))//"/unsteady_modelling/added_mass")) then
            Actuatorline%do_added_mass=.true.
        endif

        if(have_option(trim(actuatorline_path(i))//"/unsteady_modelling/added_mass")) then
            Actuatorline%do_dynamic_stall=.true.
        endif

        !##################4 Get Pitching Opions ##################
        if(have_option(trim(actuatorline_path(i))//"/pitch_control")) then
            Actuatorline%pitch_control=.true.
            call get_option(trim(actuatorline_path(i))//"/pitch_control/start_time",actuatorline(i)%pitch_start_time)
            call get_option(trim(actuatorline_path(i))//"/pitch_control/end_time",actuatorline(i)%pitch_end_time)

            if(have_option(trim(actuatorline_path(i))//"/pitch_control/harmonic")) then
                call get_option(trim(actuatorline_path(i))//"/pitch_control/harmonic/initial_pitch_angle",actuatorline(i)%pitch_angle_init)
                call get_option(trim(actuatorline_path(i))//"/pitch_control/harmonic/pitch_amplitude",actuatorline(i)%pitchAmp)
                call get_option(trim(actuatorline_path(i))//"/pitch_control/harmonic/angular_pitching_frequency",actuatorline(i)%angular_pitch_freq)
            endif
        endif

        end do

        ewrite(2,*) 'Exiting get_actuatorline_options'

    end subroutine get_actuatorline_options 

    subroutine actuator_line_model_update(current_time,dt)

        implicit none
        real, intent(in) :: dt, current_time
        integer :: i,j, Nstation
        real :: theta 

        ewrite(1,*) 'Entering the actuator_line_model_update'
        ! This routine updates the location of the actuator lines

        ctime=current_time
        DeltaT=dt

        if (Ntur>0) then
            do i=1,Ntur
            if(Turbine(i)%Is_constant_rotation_operated) then
                theta=Turbine(i)%angularVel*DeltaT
                Turbine(i)%AzimAngle=Turbine(i)%AzimAngle+theta
                call rotate_turbine(turbine(i),theta)
                call Compute_Turbine_RotVel(Turbine(i))  
            endif
            enddo
        endif

        if (Nal>0) then
            do i=1,Nal
            if(ActuatorLine(i)%pitch_control.and.ctime > ActuatorLine(i)%pitch_start_time.and.ctime < ActuatorLine(i)%pitch_end_time) then    
                !> Do harmonic pitch control for all elements of the actuator line
                Nstation=ActuatorLine(i)%NElem+1
                do j=1,Nstation
                ActuatorLine(i)%pitch(j)=ActuatorLine(i)%pitch_angle_init+actuatorline(i)%pitchAmp*sin(actuatorline(i)%angular_pitch_freq*(ctime-ActuatorLine(i)%pitch_start_time))
                end do
                call pitch_actuator_line(actuatorline(i))
            endif
            enddo
        endif
        ewrite(1,*) 'Exiting the actuator_line_model_update'

    end subroutine actuator_line_model_update

    subroutine actuator_line_model_compute_forces

        implicit none

        integer :: i,j,k
        real :: theta, pitchangle, omega
        ! Zero the Source Term at each time step

        ewrite(1,*) 'Entering the actuator_line_model_compute_forces'

        if (Ntur>0) then

            ! Get into each Turbine and Compute the Forces blade by blade and element by element
            do i=1,Ntur 
            ! Blades
            do j=1,Turbine(i)%Nblades
            call Compute_ActuatorLine_Forces(Turbine(i)%Blade(j),visc,deltaT)    
            end do
            call Compute_Turbine_Tip_Correction(Turbine(i))
            call Compute_performance(Turbine(i))

            ! Tower
            if(Turbine(i)%has_tower) then
                call Compute_Tower_Forces(Turbine(i)%Tower,visc,ctime,Turbine(i)%TowerDrag,Turbine(i)%TowerLift,Turbine(i)%TowerStrouhal)
            endif
            
            end do
        end if

        if (Nal>0) then
            do i=1,Nal
            call Compute_ActuatorLine_Forces(ActuatorLine(i),visc,deltaT)
            end do
        end if

        ewrite(1,*) 'Exiting actuator_line_model_compute_forces'
        return

    end subroutine actuator_line_model_compute_forces

end module actuator_line_model
