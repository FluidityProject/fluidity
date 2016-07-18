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

    subroutine get_turbine_options
    
    implicit none
    
    character(len=OPTION_PATH_LEN)::  turbine_name, actuatorline_name
    integer :: i,j,k
    integer, parameter :: MaxReadLine = 1000    
    character(MaxReadLine) :: FN    ! path to geometry input file 
    integer :: NElem
    character(MaxReadLine) :: ReadLine
    character(len=OPTION_PATH_LEN) :: section_path
    character(len=OPTION_PATH_LEN), allocatable :: turbine_path(:), actuatorline_path(:)
    
    ewrite(2,*) 'Entering get_turbine_options'
    
    Ntur = option_count("/actuator_line_model/turbine")
    ewrite(2,*) 'Number of Turbines : ', Ntur
    
    ! Allocate Turbines Arrays
    Allocate(Turbine(Ntur))
    Allocate(turbine_path(Ntur))
    
    ! ==========================================
    ! Get Turbines' options and INITIALIZE THEM
    ! ==========================================
    do i=1, Ntur 
       turbine_path(i)="/actuator_line_model/turbine["//int2str(i-1)//"]"
       call get_option("/actuator_line_model/turbine["//int2str(i-1)//"]/name",Turbine(i)%name)
       call get_option("/actuator_line_model/turbine["//int2str(i-1)//"]/blade_specs_file/file_name",Turbine(i)%blade_specs_file)
       ! Count how many Airfoil Sections are available
       Turbine(i)%NAirfoilData=option_count("/actuator_line_model/turbine["//int2str(i-1)//"]/airfoil_sections/section") 
       ewrite(2,*) 'Number of Airfoils available : ', Turbine(i)%NAirfoilData
       ! Allocate the memory of the Airfoils
       Allocate(Turbine(i)%AirfoilData(Turbine(i)%NAirfoilData))
       
       
       do k=1, Turbine(i)%NAirfoilData

        call get_option(trim(turbine_path(i))//"/airfoil_sections/section["//int2str(k-1)//"]/airfoil_file",Turbine(i)%AirfoilData(k)%afname)
           
           ! Read and Store Airfoils
           call airfoil_init_data(Turbine(i)%AirfoilData(k))
       end do

   !########## Get turbine_specs #################
   ! Check the typ of Turbine (choose between Horizontal and Vertical Axis turbines) 
   if(have_option(trim(turbine_path(i))//"/turbine_specs/type/Horizontal_Axis")) then
        Turbine(i)%Type='Horizontal_Axis'
        call get_option(trim(turbine_path(i))//"/turbine_specs/type/Horizontal_Axis/Number_of_blades",Turbine(i)%Nblades)
        call get_option(trim(turbine_path(i))//"/turbine_specs/type/Horizontal_Axis/origin",Turbine(i)%origin)
        call get_option(trim(turbine_path(i))//"/turbine_specs/type/Horizontal_Axis/hub_tilt_angle",Turbine(i)%hub_tilt_angle)
        call get_option(trim(turbine_path(i))//"/turbine_specs/type/Horizontal_Axis/blade_cone_angle",Turbine(i)%blade_cone_angle)
        call get_option(trim(turbine_path(i))//"/turbine_specs/type/Horizontal_Axis/yaw_angle",Turbine(i)%yaw_angle)
    elseif(have_option(trim(turbine_path(i))//"/turbine_specs/type/Vertical_Axis")) then
    
        FLExit("At the moment only the Horizontal_Axis Turbine is available")
    else
        FLExit("You should not be here")
    end if
   
   !##############3 Get Operation Options ######################
       if (have_option(trim(turbine_path(i))//"/operation/constant_rotational_velocity")) then
            Turbine(i)%Is_constant_rotation_operated= .true.
            call get_option("/actuator_line_model/turbine["//int2str(i-1)//"]/operation/constant_rotational_velocity/omega",Turbine(i)%angularVel)
            call get_option("/actuator_line_model/turbine["//int2str(i-1)//"]/operation/constant_rotational_velocity/TSR",Turbine(i)%TSR)
            if(have_option(trim(turbine_path(i))//"/operation/constant_rotational_velocity/rotation_direction/clockwise")) then
                Turbine(i)%IsClockwise=.true.
            elseif(have_option(trim(turbine_path(i))//"/operation/constant_rotational_velocity/rotation_direction/counter_clockwise")) then
                Turbine(i)%IsCounterClockwise=.true.
            else
                FLExit("You should not be here. The options are clockwise and counterclockwise")
            endif 
        else if(have_option(trim("/actuator_line_model/turbine["//int2str(i-1)//"]")//"/operation/force_based_rotational_velocity")) then
            Turbine(i)%Is_force_based_operated = .true. 
       else
           FLExit("At the moment only the constant and the force_based rotational velocity models are supported") 
       endif

    !##################4 Get Unsteady Modelling Options ##################
    !if(have_option(trim(turbine_path(i))//"/unsteady_modelling/added_mass")) then

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
    
    Nal = option_count("/actuator_line_model/actuatorline")
    ewrite(2,*) 'Number of Actuatorlines : ', Nal
    
    ! Allocate Turbines Arrays
    allocate(actuatorline(Nal))
    Allocate(actuatorline_path(Nal))
    
    ! ==========================================
    ! Get Turbines' options and INITIALIZE THEM
    ! ==========================================
    do i=1, Nal
       actuatorline_path(i)="/actuator_line_model/actuatorline["//int2str(i-1)//"]"
       call get_option("/actuator_line_model/actuatorline["//int2str(i-1)//"]/name",Actuatorline(i)%name)
       call get_option("/actuator_line_model/actuatorline["//int2str(i-1)//"]/geometry_file/file_name",Actuatorline(i)%geom_file)
       
       ! Count how many Airfoil Sections are available
       Actuatorline(i)%NAirfoilData=option_count("/actuator_line_model/actuatorline["//int2str(i-1)//"]/airfoil_sections/section") 
       ewrite(2,*) 'Number of Airfoils available : ', Actuatorline(i)%NAirfoilData
       ! Allocate the memory of the Airfoils
       Allocate(Actuatorline(i)%AirfoilData(Actuatorline(i)%NAirfoilData))
        
       do k=1, Actuatorline(i)%NAirfoilData
           
        call get_option(trim(Actuatorline_path(i))//"/airfoil_sections/section["//int2str(k-1)//"]/airfoil_file",Actuatorline(i)%AirfoilData(k)%afname)
           
           ! Read and Store Airfoils
           call airfoil_init_data(Actuatorline(i)%AirfoilData(k))
       end do
 
       !############## Get Actuator line DoF ######################
       Actuatorline(i)%NDoF=option_count("actuator_line_model/actuatorline/DoF")       
       
       allocate(actuatorline(i)%RotN(Actuatorline(i)%NDoF,3))
       ! So far only pitch works
       do iDof=1,actuatorline(i)%NDoF
       if (have_option(trim(actuatorline_path(i))//"/DoF/Pitch")) then
           actuatorline(i)%RotN(iDoF,:)=actuatorline(i)%SpanWise
       else if(have_option(trim(actuatorline_path(i))//"/DoF/Roll")) then
           actuatorline(i)%RotN(iDoF,:)=(/1.0,0.0,0.0/)
       else if(have_option(trim(actuatorline_path(i))//"/DoF/Yaw")) then
           actuatorline(i)%RotN(iDoF,:)=(/0.0,0.0,1.0/) 
       else 
           FLExit("Option not available. Options are: Pitch, Roll and Yaw ") 
       endif
       end do
       
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
        do j=1,Turbine(i)%Nblades
            call Compute_ActuatorLine_Tip_Correction(Turbine(i)%Blade(j),Turbine(i)%TSR,Turbine(i)%Nblades,Turbine(i)%Rmax)
            call Compute_ActuatorLine_Forces(Turbine(i)%Blade(j),visc,deltaT)    
            ! Compute Rotor Torque --here
        end do
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
