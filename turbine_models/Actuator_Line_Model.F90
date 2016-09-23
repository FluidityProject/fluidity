!    Copyright (C) 2008 Imperial College London and others.
!    
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
    use dynstall

    implicit none

    type(ActuatorLineType), allocatable, save :: Actuatorline(:)
    type(TurbineType), allocatable, save :: Turbine(:) ! Turbine 
    integer,save :: Ntur, Nal ! Number of the turbines 
    real,save :: deltaT, Visc
    logical,save :: actuator_line_model_writeFlag=.false.

    public  actuator_line_model_init, actuator_line_model_compute_forces, actuator_line_model_update 

contains
    
    subroutine actuator_line_model_init
   
    implicit none
    integer :: itur,ial 
    
    ewrite(1,*) 'Entering the actuator_line_model_init '

    ! Get Global Time 
    call get_option("/timestepping/timestep",deltaT) 

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
    
    subroutine actuator_line_model_init_output
        
        implicit none
        integer :: itur,ial

        if (Ntur>0) then
        call init_turbine_output_file
        endif
        
        if (Nal>0) then
        call init_actuator_line_output_file
        endif

    end subroutine actuator_line_model_init_output
    
    subroutine actuator_line_model_write_output
        
        implicit none
        integer :: itur,ial

        if (Ntur>0) then
        do itur=1,Ntur
        call write_turbine_output_file(Turbine(itur))
        end do
        endif
        
        if (Nal>0) then
        do ial=1,Nal
        call write_actuator_line_output_file(actuatorline(ial))
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

        ! ## HUB ?
        if (have_option("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/hub")) then
        Turbine(i)%Has_Hub=.true.
        call get_option("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/hub/hub_geometry/file_name",Turbine(i)%Hub%geom_file)
        nfoils=option_count("turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/hub/static_foil_data/foil")
        ewrite(2,*) 'Number of Static Foil Data available for the analysis of the hub: ', nfoils
        Allocate(Turbine(i)%hub%AirfoilData(nfoils))
        
        turbine(i)%hub%NAirfoilData=nfoils
        
        do k=1, Turbine(i)%hub%NAirfoilData
            call get_option(trim(turbine_path(i))//"/hub/static_foil_data/foil["//int2str(k-1)//"]/foil_file",Turbine(i)%hub%AirfoilData(k)%afname)   
            ! Read and Store Airfoils
            call airfoil_init_data(Turbine(i)%hub%AirfoilData(k))
        end do

        endif
        
        ! ## Tower ?
        if (have_option("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/tower")) then
        Turbine(i)%Has_Tower=.true.
        call get_option("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/tower/offset",Turbine(i)%TowerOffset)
        call get_option("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/tower/tower_geometry/file_name",Turbine(i)%Tower%geom_file)
        nfoils=option_count("/turbine_models/actuator_line_model/turbine["//int2str(i-1)//"]/tower/static_foil_data/foil")
        ewrite(2,*) 'Number of Static Foil Data available for the analysis of the tower: ', nfoils
        Allocate(Turbine(i)%tower%AirfoilData(nfoils))
       
        turbine(i)%Tower%NAirfoilData=nfoils

        do k=1, Turbine(i)%tower%NAirfoilData
            call get_option(trim(turbine_path(i))//"/tower/static_foil_data/foil["//int2str(k-1)//"]/foil_file",Turbine(i)%tower%AirfoilData(k)%afname)   
            ! Read and Store Airfoils
            call airfoil_init_data(Turbine(i)%tower%AirfoilData(k))
        end do

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
    
  
   end do
    
   ewrite(2,*) 'Exiting get_actuatorline_options'

    end subroutine get_actuatorline_options 

    subroutine actuator_line_model_update

    implicit none
    integer :: i
    real :: theta
    ! This routine updates the location of the actuator lines
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

    end subroutine actuator_line_model_update

    subroutine actuator_line_model_compute_forces
 
    implicit none
    
    integer :: i,j,k
    real :: theta, pitchangle, omega
    ! Zero the Source Term at each time step
        
    ewrite(1,*) 'Entering the actuator_line_model_update'
     
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
            call Compute_ActuatorLine_Forces(Turbine(i)%Tower,visc,deltaT)
        endif
        
        ! Hub
        if(Turbine(i)%has_hub) then
            call Compute_ActuatorLine_Forces(Turbine(i)%hub,visc,deltaT)
        endif

        end do
    end if

    if (Nal>0) then
    do i=1,Nal
    call Compute_ActuatorLine_Forces(ActuatorLine(i),visc,deltaT)
    end do
    end if

    ewrite(1,*) 'Exiting actuator_line_model_update'
    return

    end subroutine actuator_line_model_compute_forces
 
end module actuator_line_model
