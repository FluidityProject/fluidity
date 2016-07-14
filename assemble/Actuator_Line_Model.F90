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
  use sparse_tools
  use fetools
  use fields
  use futils
  use elements
  use transform_elements
  use state_module
  use state_fields_module
  use boundary_conditions
  use solvers
  use python_state
  use sparsity_patterns_meshes
  use field_options
  use fefields

  !use the ALTurbine Modules
  use airfoils
  use actuator_line_model_utils
  use dynstall

  implicit none
!GGGGG
! Define the types that will be used

type ActuatorLineType
    integer :: NElem                    ! Number of Elements of the Blade
    character(len=100):: name           ! Actuator line name
    character(len=100):: geom_file      ! Actuator line file name (is not used for the turbines)

    !#######################################################################
    ! Station Values
    logical :: FlipN =.false.           ! Flip Normal
    real, allocatable :: QCx(:)         ! Blade quarter-chord line x coordinates at element ends
    real, allocatable :: QCy(:)         ! Blade quarter-chord line y coordinates at element ends
    real, allocatable :: QCz(:)         ! Blade quarter-chord line z coordinates at element ends
    real, allocatable :: tx(:)          ! Blade unit tangent vector (rearward chord line direction) x-comp at element ends
    real, allocatable :: ty(:)          ! Blade unit tangent vector (rearward chord line direction) y-comp at element ends 
    real, allocatable :: tz(:)          ! Blade unit tangent vector (rearward chord line direction) z-comp at element ends  
    real, allocatable :: C(:)        ! Blade chord length at element ends
    real, allocatable :: thick(:) 
    !######################################################################
    
    !################################################################
    ! Element Element Values
    real, allocatable :: PEx(:)         ! Element centre x coordinates
    real, allocatable :: PEy(:)         ! Element centre y coordinates
    real, allocatable :: PEz(:)         ! Element centre z coordinates
    real, allocatable :: tEx(:)         ! Element unit tangent vector (rearward chord line direction) x-component
    real, allocatable :: tEy(:)         ! Element unit tangent vector (rearward chord line direction) y-component
    real, allocatable :: tEz(:)         ! Element unit tangent vector (rearward chord line direction) z-component
    real, allocatable :: nEx(:)         ! Element unit normal vector x-component
    real, allocatable :: nEy(:)         ! Element unit normal vector y-component
    real, allocatable :: nEz(:)         ! Element unit normal vector z-component
    real, allocatable :: sEx(:)         ! Element unit spanwise vector x-component 
    real, allocatable :: sEy(:)         ! Element unit spanwise vector y-component
    real, allocatable :: sEz(:)         ! Element unit spanwise vector z-component
    real, allocatable :: EC(:)          ! Element chord lenght
    real, allocatable :: EDS(:)         ! Element spanwise distance (length)
    real, allocatable :: EArea(:)       ! Element Area
    real, allocatable :: Eepsilon(:)    ! Element Force Projection Parameter
    real, allocatable :: ERdist(:)      ! Element Distance from the origin 
    real, allocatable :: ETtoC(:)       ! Element thickness to Chord ratio
    real, allocatable :: EAOA(:)        ! Element Last angle of Attack (used in added mass terms)
    real, allocatable :: EUn(:)         ! Element Last normal velocity (used in added mass terms)
    real, allocatable :: EAOA_LAST(:)   ! Element Last angle of Attack (used in added mass terms)
    real, allocatable :: EUn_LAST(:)    ! Element Last normal velocity (used in added mass terms)
    
    ! Velocity of the Fluid at the actuator Line Locations
    real, allocatable :: EVx(:)         ! Element Local fluid Velocity in the global x-direction
    real, allocatable :: EVy(:)         ! Element Local fluid Velocity in the global y-direction
    real, allocatable :: EVz(:)         ! Element Local fluid Velocity in the global z-direction
    
    ! Body Velocity of the Actuator Line 
    real, allocatable :: EVbx(:)        ! Element Local body Velocity in the global x-direction
    real, allocatable :: EVby(:)        ! Element Local body Velocity in the global y-direction
    real, allocatable :: EVbz(:)        ! Element Local body Velocity in the global z-direction

    ! Element Forces in the nts direction
    real, allocatable :: EFn(:)         ! Element Force in the normal direction
    real, allocatable :: EFt(:)         ! Element Force in the tangential direction (rearward chord line direction) 
    real, allocatable :: EFs(:)         ! Element Force in the spanwise direction

    ! Element Forces and Torque in the xyz direction
    real, allocatable :: EFx(:)         ! Element Force in the global x-direction
    real, allocatable :: EFy(:)         ! Element Force in the global y-direction
    real, allocatable :: EFz(:)         ! Element Force in the global z-direction
    real, allocatable :: ETRx(:)     ! Element Torque over the point of rotation 
    real, allocatable :: ETRy(:)     ! Element Torque over the point of rotation 
    real, allocatable :: ETRz(:)     ! Element Torque over the point of rotation 

    
    ! Element Airfoil Data
    type(AirfoilType), allocatable :: EAirfoil(:) ! Element Airfoil 
    integer :: NAirfoilData
    type(AirfoilType), allocatable :: AirfoilData(:) ! Element Airfoil 
    type(LB_Type), allocatable :: E_LB_Model(:)   ! Element Leishman-Beddoes Model
    !##########################################################################################

    !##########################################################################################
    ! Forces and Torques on the ActuatorLine 
    real :: Fn     ! Force in the normal direction
    real :: Ft     ! Force in the tangential direction (rearward chord line direction) 
    real :: Fs     ! Force in the spanwise direction

    real :: Fx     ! Element Force in the global x-direction
    real :: Fy     ! Element Force in the global y-direction
    real :: Fz     ! Element Force in the global z-direction
    real :: TRX    ! Element Torque X over the point of rotation 
    real :: TRY    ! Element Torque Y over the point of rotation 
    real :: TRZ    ! Element Torque Z over the point of rotation 

    real :: Area   ! Effective Airfoil Area

    ! Degrees of Freedom
    
    real :: COR(3)       ! Center of Rotation
    real :: SpanWise(3)   ! Point of Rotation
    integer          :: NDoF          ! Number of Degrees of Freedom
    real,allocatable :: RotN(:,:)     ! Maximum Degrees of Freedom 6 
    
end type ActuatorLineType

type TurbineType
    character(len=100) :: name
    character(len=100) :: blade_specs_file
    character(len=100) :: type
    integer :: NBlades, NAirfoilData
    real, dimension(3) :: RotN, origin ! Rotational vectors in the normal and perpendicular directions
    real :: hub_tilt_angle, blade_cone_angle, yaw_angle 
    real :: Rmax ! Reference radius, velocity, viscosity
    real :: A   ! Rotor area
    real :: angularVel
    real :: AzimAngle=0.0
    logical :: Is_constant_rotation_operated = .false. ! For a constant rotational velocity (in Revolutions Per Minute)
    logical :: Is_force_based_operated = .false. ! For a forced based rotational velocity (Computed during the simulation)
    logical :: IsRotating = .false.
    logical :: IsClockwise = .false.
    logical :: IsCounterClockwise = .false. 
    logical :: Has_Hub, Has_Tower
    type(ActuatorLineType), allocatable :: Blade(:)
    type(AirfoilType), allocatable :: AirfoilData(:)
    type(ActuatorLineType) :: hub, tower
    real :: CP  ! Power coefficient 
    real :: CTR ! Torque coefficient 
    real :: CFx ! Fx coefficient 
    real :: CFy ! Fy coefficient 
    real :: CFz ! Fz coefficient 
    real :: CT  ! Thrust coefficient    
end type TurbineType


    type(ActuatorLineType), allocatable, save :: Actuatorline(:)
    type(TurbineType), allocatable, save :: Turbine(:) ! Turbine 
    real, allocatable :: Sx(:),Sy(:),Sz(:),Se(:),Su(:),Sv(:),Sw(:),Sc(:),Sfx(:),Sfy(:),Sfz(:)
    integer,save :: Ntur, Nal ! Number of the turbines 
    integer,save :: NSource
    real,save :: deltaT, Visc
 
    public  actuator_line_model_init, actuator_line_model_compute_forces 

contains
    
    subroutine actuator_line_model_init
   
        implicit none
        integer :: itur,ial 
    ewrite(1,*) 'Entering the actuator_line_model_init '

    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    !* This a routine from which fluidity reads the information about the turbine model 
    !* We start with a single turbine model 
    ! In the final version of the code we should be able to enable 
    ! multiple turbines. 
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

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

    !############## Specify Source Terms ########################
    call allocate_source_terms

    ewrite(1,*) 'Exiting the actuator_line_model_init '

    end subroutine actuator_line_model_init
    
    subroutine allocate_source_terms

    implicit none
    integer :: counter,itur,iblade,ielem,ial
    
    ewrite(1,*) 'entering allocate_source_terms'
    counter=0
    if (Ntur>0) then
        do itur=1,Ntur
            do iblade=1,Turbine(itur)%Nblades
                do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
                counter=counter+1
                end do
            end do
        end do
    endif
    
    if (Nal>0) then
        do ial=1,Nal
            do ielem=1,actuatorline(ial)%NElem
                counter=counter+1
            end do
        end do
    endif
    NSource=counter
    allocate(Sx(NSource),Sy(NSource),Sz(NSource),Sc(Nsource),Su(NSource),Sv(NSource),Sw(NSource),Se(NSource),Sfx(NSource),Sfy(NSource),Sfz(NSource))

    ewrite(1,*) 'exiting allocate_source_terms'

    end subroutine allocate_source_terms

    subroutine get_locations
    
    implicit none
    integer :: counter,itur,iblade,ielem,ial

    counter=0
    if (Ntur>0) then
        do itur=1,Ntur
            do iblade=1,Turbine(itur)%Nblades
                do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
                counter=counter+1
                Sx(counter)=Turbine(itur)%Blade(iblade)%PEX(ielem)
                Sy(counter)=Turbine(itur)%Blade(iblade)%PEY(ielem)
                Sz(counter)=Turbine(itur)%Blade(iblade)%PEZ(ielem)
                Sc(counter)=Turbine(itur)%Blade(iblade)%EC(ielem)
                end do
            end do
        end do
    endif
    
    if (Nal>0) then
        do ial=1,Nal
            do ielem=1,actuatorline(ial)%NElem
                counter=counter+1
                Sx(counter)=actuatorline(ial)%PEX(ielem)
                Sy(counter)=actuatorline(ial)%PEY(ielem)
                Sz(counter)=actuatorline(ial)%PEZ(ielem)
                Sc(counter)=actuatorline(ial)%EC(ielem)
            end do
        end do
    endif
  
    end subroutine get_locations
    
    subroutine set_vel
    
    implicit none
    integer :: counter,itur,iblade,ielem,ial
    
    ewrite(1,*) 'Entering set_vel'
    counter=0
    if (Ntur>0) then
        do itur=1,Ntur
            do iblade=1,Turbine(itur)%Nblades
                do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
                counter=counter+1
                Turbine(itur)%Blade(iblade)%EVX(ielem)=Su(counter)
                Turbine(itur)%Blade(iblade)%EVY(ielem)=Sv(counter)
                Turbine(itur)%Blade(iblade)%EVZ(ielem)=Sw(counter)
                end do
            end do
        end do
    endif
    
    if (Nal>0) then
        do ial=1,Nal
            do ielem=1,actuatorline(ial)%NElem
                counter=counter+1
                actuatorline(ial)%EVX(ielem)=Su(counter)
                actuatorline(ial)%EVY(ielem)=Sv(counter)
                actuatorline(ial)%EVZ(ielem)=Sw(counter)
            end do
        end do
    endif

    ewrite(1,*) 'Exiting set_vel'
  
    end subroutine set_vel
    
    subroutine get_forces
    
    implicit none
    integer :: counter,itur,iblade,ielem,ial

    counter=0
    if (Ntur>0) then
        do itur=1,Ntur
            do iblade=1,Turbine(itur)%Nblades
                do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
                counter=counter+1
                Sfx(counter)=Turbine(itur)%Blade(iblade)%EFX(ielem)
                Sfy(counter)=Turbine(itur)%Blade(iblade)%EFY(ielem)
                Sfz(counter)=Turbine(itur)%Blade(iblade)%EFZ(ielem)
                end do
            end do
        end do
    endif
    
    if (Nal>0) then
        do ial=1,Nal
            do ielem=1,actuatorline(ial)%NElem
                counter=counter+1
                Sfx(counter)=actuatorline(ial)%EFX(ielem)
                Sfy(counter)=actuatorline(ial)%EFY(ielem)
                Sfz(counter)=actuatorline(ial)%EFZ(ielem)
            end do
        end do
    endif
  
    end subroutine get_forces
    
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
    Allocate(actuatorline(Nal))
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
            call rotate_turbines(theta)
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
    do i=1,Ntur
        do j=1,Turbine(i)%Nblades
            call Compute_ActuatorLine_Forces(Turbine(i)%Blade(j))
        end do
    end do
    end if

    if (Nal>0) then
    do i=1,Nal
    call Compute_ActuatorLine_Forces(ActuatorLine(i))
    ewrite(2,*) 'Forces X,Y,Z', ActuatorLine(i)%FX, ActuatorLine(i)%FY, ActuatorLine(i)%FZ  
    end do
    end if

    ewrite(1,*) 'Exiting actuator_line_model_update'
    return

end subroutine actuator_line_model_compute_forces

subroutine set_turbine_geometry(turbine)

    implicit none
    type(TurbineType),intent(inout) :: turbine
    real, allocatable :: rR(:),ctoR(:),pitch(:),thick(:)
    real :: SVec(3), theta, origin(3)
    integer :: Nstations, iblade, Istation

    ewrite(2,*) 'Entering set_turbine_geometry'

    ewrite(2,*) 'Turbine Name : ', turbine%name 
    ewrite(2,*) '============================='
    ewrite(2,*) 'Number of Blades : ', turbine%Nblades
    ewrite(2,*) 'Origin           : ', turbine%origin
    ewrite(2,*) 'Hub Tilt Angle   : ', turbine%hub_tilt_angle
    ewrite(2,*) 'Blade Cone Angle : ', turbine%blade_cone_angle
    ewrite(2,*) 'Yaw Angle        : ', turbine%yaw_angle

    allocate(turbine%blade(turbine%Nblades))

    call read_actuatorline_geometry(turbine%blade_specs_file,turbine%Rmax,SVec,rR,ctoR,pitch,thick,Nstations)
    ! Make sure that the spanwise is [0 0 1]
    Svec = (/0.0,0.0,1.0/)
    ! Make sure that origin is [0,0,0] : we set everything to origin 0 and then translate the
    ! turbine to the actual origin(this is for simplicity)
    theta=2*pi/turbine%Nblades
    do iblade=1,turbine%Nblades
    call allocate_actuatorline(Turbine%blade(iblade),Nstations)
    turbine%blade(iblade)%name=turbine%name//int2str(iblade)
    
    turbine%blade(iblade)%COR=turbine%origin
    
    do istation=1,Nstations
    turbine%blade(iblade)%QCx(istation)=rR(istation)*turbine%Rmax*Svec(1)+turbine%blade(iblade)%COR(1)
    turbine%blade(iblade)%QCy(istation)=rR(istation)*turbine%Rmax*Svec(2)+turbine%blade(iblade)%COR(2)
    turbine%blade(iblade)%QCz(istation)=rR(istation)*turbine%Rmax*Svec(3)+turbine%blade(iblade)%COR(3)
    if(turbine%IsCounterClockwise) then
        turbine%RotN=(/-1.0,0.0,0.0/)
        turbine%blade(iblade)%tx(istation)=sin(pitch(istation)/180.0*pi)    
        turbine%blade(iblade)%ty(istation)=-cos(pitch(istation)/180.0*pi)    
        turbine%blade(iblade)%tz(istation)= 0.0
        turbine%blade(iblade)%C(istation)=ctoR(istation)*turbine%Rmax
        turbine%blade(iblade)%thick(istation)=thick(istation)
    elseif(turbine%IsClockwise) then
        turbine%RotN=(/1.0,0.0,0.0/)
        turbine%blade(iblade)%tx(istation)=sin(pitch(istation)/180.0*pi)    
        turbine%blade(iblade)%ty(istation)=cos(pitch(istation)/180.0*pi)    
        turbine%blade(iblade)%tz(istation)= 0.0
        turbine%blade(iblade)%C(istation)=ctoR(istation)*turbine%Rmax
        turbine%blade(iblade)%thick(istation)=thick(istation)
        turbine%blade(iblade)%FlipN = .true.
    endif
    end do

    ! Always rotate counterclockwise to assign the turbine blades
    call rotate_actuatorline(turbine%blade(iblade),turbine%blade(iblade)%COR,(/-1.0,0.0,0.0/),(iblade-1)*theta)   
    ! Rotate through incidence (hub tilt) and coning angle
    call make_actuatorline_geometry(turbine%blade(iblade))
    ! Populate element Airfoils 
    call populate_blade_airfoils(turbine%blade(iblade)%NElem,turbine%Blade(iblade)%EAirfoil,turbine%AirfoilData,turbine%Blade(iblade)%ETtoC)
    
    turbine%Blade(iblade)%EAOA_LAST(:)=1.e7
    turbine%Blade(iblade)%EUn_LAST(:)=1.e7
    
    end do
    
    call Compute_Turbine_RotVel

    ewrite(2,*) 'Exiting set_turbine_geometry'

end subroutine set_turbine_geometry

subroutine set_actuatorline_geometry(actuatorline)

    implicit none
    type(ActuatorLineType),intent(inout) :: actuatorline
    real, allocatable :: rR(:),ctoR(:),pitch(:),thick(:)
    real :: SVec(3), theta, origin(3), length 
    integer :: Nstations, iact, Istation, ielem

    ewrite(2,*) 'Entering set_actuatorline_geometry'

    ewrite(2,*) 'Actuatorline Name : ', actuatorline%name 
    ewrite(2,*) '============================='

    call read_actuatorline_geometry(actuatorline%geom_file,length,SVec,rR,ctoR,pitch,thick,Nstations)
    
    call allocate_actuatorline(actuatorline,Nstations)
   
    actuatorline%SpanWise=SVec
    actuatorline%COR=(/0.0, 0.0, 0.0/)

    ! The directions of vectors etc are just hacked ...
    do istation=1,Nstations
    actuatorline%QCx(istation)=rR(istation)*length*Svec(1)
    actuatorline%QCy(istation)=rR(istation)*length*Svec(2)
    actuatorline%QCz(istation)=rR(istation)*length*Svec(3)      
    actuatorline%tx(istation)= cos(pitch(istation)/180.0*pi)    
    actuatorline%ty(istation)= 0.0
    actuatorline%tz(istation)= -sin(pitch(istation)/180.0*pi)     
    actuatorline%C(istation)=ctoR(istation)*length
    actuatorline%thick(istation)=thick(istation)
    actuatorline%FlipN=.true.
    end do
   
    call make_actuatorline_geometry(actuatorline)
    ! Compute the Area

    do ielem=1,actuatorline%Nelem
    actuatorline%Area=actuatorline%Area+actuatorline%EArea(ielem)
    end do
    
    ! Populate element Airfoils 
    call populate_blade_airfoils(actuatorline%NElem,actuatorline%EAirfoil,actuatorline%AirfoilData,actuatorline%ETtoC)
    
    actuatorline%EAOA_LAST(:)=1.e7
    actuatorline%EUn_LAST(:)=1.e7
    
    do ielem=1,actuatorline%Nelem
    call dystl_init_LB(actuatorline%E_LB_Model(ielem))
    end do
    ewrite(2,*) 'Exiting set_actuatorline_geometry'

end subroutine set_actuatorline_geometry
! Turb
subroutine calculate_performance(turbine)

    implicit none
    type(TurbineType), intent(inout) :: turbine
    real :: TRX_i,TRY_i,TRZ_i,FX_i,FY_i,FZ_i, FX,FY,FZ,TRX,TRY,TRZ
    real :: U_ref, R, A
    integer :: iblade, ielem
    ewrite(2,*) 'In calculate_performance'

    !turbine%CFx=FX/(0.5*A*U_ref**2)
    !turbine%CFy=FY/(0.5*A*U_ref**2)
    !turbine%CFz=Fz/(0.5*A*U_ref**2)
    !turbine%CT=sqrt(turbine%CFx**2.0+turbine%CFy**2.0+turbine%CFz**2.0)
    !turbine%CTR=TRX/(0.5*A*R*U_ref**2.0)
    !turbine%CP= turbine%CTR*turbine%TSR
    
    ewrite(2,*) '--------------------------------------------------------'
    ewrite(2,*) 'Calculate performance for Turbine : ',turbine%name
    ewrite(2,*) 'Azimuthal Angle (degrees) : ', turbine%AzimAngle
    ewrite(2,*) 'Thrust Coefficient : ', turbine%CT
    ewrite(2,*) 'Torque Coefficient : ', turbine%CTR
    ewrite(2,*) 'Power Coefficient : ', turbine%CP
    ewrite(2,*) '--------------------------------------------------------'
    
    ewrite(2,*) 'Exiting calculate_performance'

end subroutine calculate_performance

subroutine Compute_Turbine_RotVel
    implicit none
    integer :: iturb,iblade,ielem
    real :: wRotX,wRotY,wRotZ,Rx,Ry,Rz,ublade,vblade,wblade
    real :: RotX,RotY,RotZ 
    ewrite(2,*) 'Entering Compute_Turbine_Local_Vel '
    
    !========================================================
    ! Compute Element local rotational velocity
    !========================================================
    do iturb=1,Ntur
    do iblade=1,Turbine(iturb)%NBlades
    do ielem=1,Turbine(iturb)%Blade(iblade)%Nelem
    
    wRotX=Turbine(iturb)%angularVel*Turbine(iturb)%RotN(1)
    wRotY=Turbine(iturb)%angularVel*Turbine(iturb)%RotN(2)
    wRotZ=Turbine(iturb)%angularVel*Turbine(iturb)%RotN(3)
    
    RotX=Turbine(iturb)%RotN(1)
    RotY=Turbine(iturb)%RotN(2)
    RotZ=Turbine(iturb)%RotN(3)
    

    Rx=-Turbine(iturb)%origin(1)+Turbine(iturb)%Blade(iblade)%PEx(ielem);
    Ry=-Turbine(iturb)%origin(2)+Turbine(iturb)%Blade(iblade)%PEy(ielem);
    Rz=-Turbine(iturb)%origin(3)+Turbine(iturb)%Blade(iblade)%PEz(ielem);
    Turbine(iturb)%Blade(iblade)%ERdist(ielem)=sqrt(Rx**2+Ry**2+Rz**2)


   ! ! Find the cross product Ublade = Omega x R
    call cross(wRotX,wRotY,wRotZ,Rx,Ry,Rz,ublade,vblade,wblade)
    
    Turbine(iturb)%Blade(iblade)%EVbx(ielem)=-ublade
    Turbine(iturb)%Blade(iblade)%EVby(ielem)=-vblade
    Turbine(iturb)%Blade(iblade)%EVbz(ielem)=-wblade
    
    end do
    end do

    end do


    ewrite(2,*) 'Entering Compute_Turbine_Local_Vel '

end subroutine Compute_Turbine_RotVel

subroutine Compute_ActuatorLine_Forces(act_line)
       
    implicit none
    type(ActuatorLineType),intent(inout) :: act_line
    real :: R(3)
    real :: wRotX,wRotY,wRotZ,Rx,Ry,Rz,ub,vb,wb,u,v,w
    real :: xe,ye,ze,nxe,nye,nze,txe,tye,tze,sxe,sye,sze,ElemArea,ElemChord
    real :: urdn,urdc, wP,ur,alpha,Re,alpha5,alpha75,adotnorm, A, B, C, dUnorm
    real :: CL,CD,CN,CT,CLCirc,CM25,MS,FN,FT,FS,FX,Fy,Fz,te, F1, g1
    real :: TRx,TRy,TRz, TRn,TRt, TRs, RotX,RotY,RotZ, dal, wPNorm, relem, ds
    integer :: ielem
  
    ewrite(2,*) 'Entering Compute_Forces '
     
    !===========================================================
    ! Assign global values to local values (to make life easier
    !==========================================================
    do ielem=1,act_line%NElem
    
    xe=act_line%PEX(ielem)
    ye=act_line%PEY(ielem)
    ze=act_line%PEZ(ielem)

    nxe=act_line%nEx(ielem)
    nye=act_line%nEy(ielem)
    nze=act_line%nEz(ielem)
    txe=act_line%tEx(ielem)
    tye=act_line%tEy(ielem)
    tze=act_line%tEz(ielem)
    sxe=act_line%sEx(ielem)
    sye=act_line%sEy(ielem)
    sze=act_line%sEz(ielem)
    ElemArea=act_line%EArea(ielem)
    ElemChord=act_line%EC(ielem) 
    u=act_line%EVx(ielem)
    v=act_line%EVy(ielem)
    w=act_line%EVz(ielem) 

    ub=act_line%EVbx(ielem)
    vb=act_line%EVby(ielem)
    wb=act_line%EVbz(ielem)
    
    !==============================================================
    ! Calculate element normal and tangential velocity components. 
    !==============================================================
    urdn=nxe*(u+ub)+nye*(v+vb)+nze*(w+wb)! Normal 
    urdc=txe*(u+ub)+tye*(v+vb)+tze*(w+wb)! Tangential
    ur=sqrt(urdn**2.0+urdc**2.0)
    alpha=atan2(urdn,urdc)
    Re = ur*ElemChord/Visc
    alpha5=alpha
    alpha75=alpha
    
    !=========================================================
    ! Compute rate of change of Unormal and angle of attack
    !=========================================================
    if(act_line%EAOA_Last(ielem)>1e6) then
    dal=0
    dUnorm=0
    else
    dal=(alpha75-act_line%EAOA_Last(ielem))
    dUnorm=urdn-act_line%EUn_last(ielem)
    endif
    adotnorm=dal/deltaT*ElemChord/(2.0*max(ur,0.001)) ! adot*c/(2*U)
    A = urdn/max(ur,0.001)
    B = ElemChord*dUnorm/(deltaT*max(ur**2,0.001))
    C = urdn*urdc/max(ur**2,0.001)

    !====================================
    ! Compute the Aerofoil Coefficients
    !====================================
    call compute_aeroCoeffs(act_line%EAirfoil(ielem),act_line%E_LB_Model(ielem),alpha75,alpha5,Re,A,B,C,adotnorm,CN,CT,CM25,CL,CLCIrc,CD)
    
    !===================================
    ! Update Dynamic Stall Model 
    !===================================
    ds=2.0*ur*DeltaT/ElemChord
    call LB_UpdateStates(act_line%E_LB_MODEL(ielem),ds)

    !========================================================
    ! Apply Coeffs to calculate tangential and normal Forces
    !========================================================
    FN=0.5*CN*ElemArea*ur**2.0
    FT=0.5*CT*ElemArea*ur**2.0
    FS=0.0 ! Makes sure that there is no spanwise force
    MS=0.5*CM25*ElemChord*ElemArea*ur**2.0

    !===============================================
    ! Compute forces in the X, Y, Z axis and torque  
    !===============================================
    FX=FN*nxe+FT*txe+FS*sxe
    FY=FN*nye+FT*tye+FS*sye
    FZ=FN*nze+FT*tze+FS*sze
    
    call cross(xe-act_line%COR(1),ye-act_line%COR(2),ze-act_line%COR(3),FN,FT,FS,TRn,TRt,TRs)
    
    TRX=TRN*nxe+TRT*txe+(TRS+MS)*sxe
    TRY=TRN*nye+TRT*tye+(TRS+MS)*sye
    TRZ=TRN*nze+TRT*tze+(TRS+MS)*sze

    !==========================================
    ! Assign the derived types
    !==========================================
    ! Local Forces
    act_line%EFN(ielem)=FN
    act_line%EFT(ielem)=FT
    act_line%EFS(ielem)=FS
    
    ! Global Forces and Torques
    act_line%EFX(ielem)=FX
    act_line%EFY(ielem)=FY
    act_line%EFZ(ielem)=FZ
    act_line%ETRX(ielem)=TRX
    act_line%ETRY(ielem)=TRY
    act_line%ETRZ(ielem)=TRZ

    !===============================================
    !! Set the AOA_LAST before exiting the routine
    !===============================================
    act_line%EAOA_LAST(ielem)=alpha75 
    act_line%EUn_last(ielem)=urdn 
    end do

    !=============================================
    ! Compute Total Forces
    !============================================
    act_line%FX=0.0
    act_line%FY=0.0
    act_line%FZ=0.0
    
    do ielem=1,act_line%NElem
        act_line%FX=act_line%FX+act_line%EFX(ielem)
        act_line%FY=act_line%FY+act_line%EFY(ielem)
        act_line%FZ=act_line%FZ+act_line%EFZ(ielem)
    end do

    ewrite(2,*) 'Exiting Compute_Forces'

end subroutine compute_Actuatorline_Forces

subroutine populate_blade_airfoils(NElem,EAirfoil,AirfoilData,ETtoC)

    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    ! This routine initialises the airfoil struct for the blades
    ! by interpolating from the data
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    implicit none
    integer,intent(IN) :: NElem
    type(AirfoilType),dimension(:) :: EAirfoil, AirfoilData
    real,dimension(:) :: ETtoC
    real,allocatable ::  Thicks(:), diffthicks(:)
    integer :: ielem,idata,NData,imin,imax,iint

    ewrite(2,*) 'Entering populate_blade_airfoils'    
    ! We need to interpolate from two or more 
    NData=size(AirfoilData)
    allocate(Thicks(NData),diffthicks(NData))
    

    do ielem=1,NElem
    EAirfoil(ielem)%afname = int2str(ielem)
    call allocate_airfoil(EAirfoil(ielem),MaxAOAVals,MaxReVals) 
    
    Thicks(:)=0.0
    diffthicks(:)=0.0
    imin=0
    imax=0
    iint=0
    
        do idata=1,NData
             Thicks(idata)=AirfoilData(idata)%tc
             diffthicks(idata)=abs(AirfoilData(idata)%tc-ETtoC(ielem))
        end do 
        imin=minloc(Thicks,1)
        imax=maxloc(Thicks,1)
        iint=minloc(diffthicks,1)
        
        if(ETtoC(ielem)>=Thicks(imax)) then
             call copy_airfoil_values(EAirfoil(ielem),AirfoilData(imax))
        elseif(ETtoC(ielem)<=Thicks(imin)) then
             call copy_airfoil_values(EAirfoil(ielem),AirfoilData(imin))
        else
             call copy_airfoil_values(EAirfoil(ielem),AirfoilData(iint))
        endif

    end do 

    ewrite(2,*) 'Exiting populate_blade_airfoils'

end subroutine populate_blade_airfoils
! Turb
subroutine rotate_turbines(theta)
    implicit none

    real :: theta,nrx,nry,nrz,px,py,pz 
    integer :: j,ielem,i
    real :: vrx,vry,vrz,VMag
    real :: xtmp,ytmp,ztmp, txtmp, tytmp, tztmp
    ! Rotates data in blade arrays. Rotate element end geometry and recalculate element geometry.

    ewrite(1,*) 'Entering rotate_turbines'
    do i=1,Ntur
        
        ! Specify the rotation axis and the normal vector of rotation
        
        nrx=Turbine(i)%RotN(1)
        nry=Turbine(i)%RotN(2)
        nrz=Turbine(i)%RotN(3)
        
        px=Turbine(i)%origin(1)
        py=Turbine(i)%origin(2)
        pz=Turbine(i)%origin(3)


        do j=1,Turbine(i)%NBlades
            do ielem=1,Turbine(i)%Blade(j)%Nelem+1
            ! Blade end locations (quarter chord). xBE(MaxSegEnds)
            xtmp=Turbine(i)%Blade(j)%QCx(ielem)
            ytmp=Turbine(i)%Blade(j)%QCy(ielem)
            ztmp=Turbine(i)%Blade(j)%QCz(ielem)
            
            Call QuatRot(xtmp,ytmp,ztmp,theta,nrx,nry,nrz,px,py,pz,vrx,vry,vrz)
            Turbine(i)%Blade(j)%QCx(ielem)=vrx                                       
            Turbine(i)%Blade(j)%QCy(ielem)=vry                                       
            Turbine(i)%Blade(j)%QCz(ielem)=vrz                                  
            
            txtmp=Turbine(i)%Blade(j)%tx(ielem)
            tytmp=Turbine(i)%Blade(j)%ty(ielem)
            tztmp=Turbine(i)%Blade(j)%tz(ielem)
            
            ! Tangent vectors
            Call QuatRot(txtmp,tytmp,tztmp,theta,nrx,nry,nrz,px,py,pz,vrx,vry,vrz)
            VMag=sqrt(vrx**2+vry**2+vrz**2)
            Turbine(i)%Blade(j)%tx(ielem)=vrx/VMag                                      
            Turbine(i)%Blade(j)%ty(ielem)=vry/VMag                                  
            Turbine(i)%Blade(j)%tz(ielem)=vrz/VMag                                       
  
            end do
            
            call make_actuatorline_geometry(Turbine(i)%Blade(j))
        end do 
    end do
    
    ewrite(1,*) 'Exiting rotate_turbines'

end subroutine rotate_turbines  

subroutine rotate_actuatorline(actuatorline,origin,rotN,theta)

    implicit none
    type(ActuatorLineType),intent(inout) :: actuatorline
    real,intent(in) :: origin(3), rotN(3)
    real :: theta,nrx,nry,nrz,px,py,pz 
    integer :: j,ielem,i
    real :: vrx,vry,vrz,VMag
    real :: xtmp,ytmp,ztmp, txtmp, tytmp, tztmp
    ! Rotates data in blade arrays. Rotate element end geometry and recalculate element geometry.

    ewrite(1,*) 'Entering rotate_actuatorline'
        
        ! Specify the rotation axis and the normal vector of rotation
        
        nrx=rotN(1)
        nry=RotN(2)
        nrz=RotN(3)
        
        px=origin(1)
        py=origin(2)
        pz=origin(3)


        do ielem=1,actuatorline%Nelem+1
        ! Blade end locations (quarter chord). xBE(MaxSegEnds)
        xtmp=actuatorline%QCx(ielem)
        ytmp=actuatorline%QCy(ielem)
        ztmp=actuatorline%QCz(ielem)
        
        Call QuatRot(xtmp,ytmp,ztmp,theta,nrx,nry,nrz,px,py,pz,vrx,vry,vrz)
        actuatorline%QCx(ielem)=vrx                                       
        actuatorline%QCy(ielem)=vry                                       
        actuatorline%QCz(ielem)=vrz                                  
        
        txtmp=actuatorline%tx(ielem)
        tytmp=actuatorline%ty(ielem)
        tztmp=actuatorline%tz(ielem)
        
        ! Tangent vectors
        Call QuatRot(txtmp,tytmp,tztmp,theta,nrx,nry,nrz,px,py,pz,vrx,vry,vrz)
        VMag=sqrt(vrx**2+vry**2+vrz**2)
        actuatorline%tx(ielem)=vrx/VMag                                      
        actuatorline%ty(ielem)=vry/VMag                                  
        actuatorline%tz(ielem)=vrz/VMag                                       
  
        end do
    
    ewrite(1,*) 'Exiting rotate_actuatorline'


end subroutine rotate_actuatorline

subroutine allocate_actuatorline(actuatorline,NStations)

    implicit none
    
    type(ActuatorLineType) :: actuatorline
    integer,intent(in) :: Nstations
    integer :: NElem
    
    Nelem=Nstations-1
    actuatorline%Nelem = Nelem

    allocate(actuatorline%QCx(NElem+1))
    allocate(actuatorline%QCy(NElem+1))
    allocate(actuatorline%QCz(NElem+1))
    allocate(actuatorline%tx(NElem+1))
    allocate(actuatorline%ty(NElem+1))
    allocate(actuatorline%tz(NElem+1))
    allocate(actuatorline%C(NElem+1))
    allocate(actuatorline%thick(NElem+1))
    allocate(actuatorline%PEx(NElem))
    allocate(actuatorline%PEy(NElem))
    allocate(actuatorline%PEz(NElem))
    allocate(actuatorline%tEx(NElem))
    allocate(actuatorline%tEy(NElem))
    allocate(actuatorline%tEz(NElem))
    allocate(actuatorline%nEx(NElem))
    allocate(actuatorline%nEy(NElem))
    allocate(actuatorline%nEz(NElem))
    allocate(actuatorline%sEx(NElem))
    allocate(actuatorline%sEy(NElem))
    allocate(actuatorline%sEz(NElem))
    allocate(actuatorline%EC(NElem))
    allocate(actuatorline%EDS(NElem))
    allocate(actuatorline%EArea(NElem))
    allocate(actuatorline%ETtoC(NElem))
    allocate(actuatorline%Eepsilon(NElem))
    allocate(actuatorline%EAirfoil(Nelem))
    allocate(actuatorline%E_LB_Model(Nelem))
    allocate(actuatorline%ERdist(Nelem))
    allocate(actuatorline%EVx(NElem))
    allocate(actuatorline%EVy(NElem))
    allocate(actuatorline%EVz(NElem))
    allocate(actuatorline%EVbx(NElem))
    allocate(actuatorline%EVby(NElem))
    allocate(actuatorline%EVbz(NElem))
    allocate(actuatorline%EAOA(Nelem))
    allocate(actuatorline%EUn(Nelem))
    allocate(actuatorline%EAOA_LAST(Nelem))
    allocate(actuatorline%EUn_LAST(Nelem))
    allocate(actuatorline%EFn(NElem))
    allocate(actuatorline%EFt(NElem))
    allocate(actuatorline%EFs(NElem))
    allocate(actuatorline%EFx(NElem))
    allocate(actuatorline%EFy(NElem))
    allocate(actuatorline%EFz(NElem))
    allocate(actuatorline%ETRx(NElem))
    allocate(actuatorline%ETRy(NElem))
    allocate(actuatorline%ETRz(NElem))
    
end subroutine allocate_actuatorline
    
subroutine read_actuatorline_geometry(FN,Rmax,SpanwiseVec,rR,ctoR,pitch,thick,Nstations)
    
    implicit none
    character(len=100),intent(in)  :: FN ! FileName of the geometry file
    real,dimension(3) :: SpanwiseVec
    real, allocatable,intent(out) :: rR(:),ctoR(:),pitch(:),thick(:)  
    real, intent(out) :: Rmax  
    integer, intent(out) :: Nstations
    integer :: i,j
    character(1000) :: ReadLine
    
    open(15,file=FN)
    ! Read the Number of Blades
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) Rmax 
    
    !Read Spanwise actuator line axis
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) SpanwiseVec(1), SpanwiseVec(2), SpanwiseVec(3)
     
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) Nstations

    allocate(rR(Nstations),ctoR(Nstations),pitch(Nstations),thick(Nstations))
    ! Read the stations specs
    do i=1,NStations
    
    read(15,'(A)') ReadLine ! Blade ....

    read(ReadLine,*) rR(i), ctoR(i), pitch(i), thick(i)

    end do
    
    close(15)


end subroutine read_actuatorline_geometry 

subroutine make_actuatorline_geometry(blade)

    implicit none

    type(ActuatorLineType),intent(INOUT) :: blade ! For simplity I leave it as blade. In fact this is an actuator line 
    integer :: BNum
    integer :: nbe, nei, nej, j
    real :: sEM, tEM, nEM, dx,dy,dz
    real :: PE(3), sE(3), tE(3), normE(3), P1(3), P2(3), P3(3), P4(3), V1(3), V2(3), V3(3), V4(3), A1(3), A2(3)

    ewrite(2,*) 'Entering make_actuatorline_geometry'
    ! Calculates element geometry from element end geometry

    ! JCM: Eventually, should just be able to loop through Blades(BNum) data structure
    ! While data is still held in arrays concatenated across blades, need to replicate
    ! nbe (stored in configr) from Blades(1).NElem
    nbe=blade%NElem

    do j=1,nbe
    nej=1+j

    ! Element center locations
    blade%PEx(nej-1)=(blade%QCx(nej)+blade%QCx(nej-1))/2.0
    blade%PEy(nej-1)=(blade%QCy(nej)+blade%QCy(nej-1))/2.0
    blade%PEz(nej-1)=(blade%QCz(nej)+blade%QCz(nej-1))/2.0

    ! Element length

    ! Set spannwise and tangential vectors
    sE=(/blade%QCx(nej)-blade%QCx(nej-1),blade%QCy(nej)-blade%QCy(nej-1),blade%QCz(nej)-blade%QCz(nej-1)/) ! nominal element spanwise direction set opposite to QC line
    sEM=sqrt(dot_product(sE,sE))

    blade%EDS(nej-1) = sEM

    sE=sE/sEM
    tE=(/blade%tx(nej)+blade%tx(nej-1),blade%ty(nej)+blade%ty(nej-1),blade%tz(nej)+blade%tz(nej-1)/)/2.0
    ! Force tE normal to sE
    tE=tE-dot_product(tE,sE)*sE
    tEM=sqrt(dot_product(tE,tE))
    tE=tE/tEM
    blade%sEx(nej-1)=sE(1)
    blade%sEy(nej-1)=sE(2)
    blade%sEz(nej-1)=sE(3)
    blade%tEx(nej-1)=tE(1)
    blade%tEy(nej-1)=tE(2)
    blade%tEz(nej-1)=tE(3)

    ! Calc normal vector
    Call cross(sE(1),sE(2),sE(3),tE(1),tE(2),tE(3),normE(1),normE(2),normE(3))
    nEM=sqrt(dot_product(normE,normE))
    normE=normE/nEM
    blade%nEx(nej-1)=normE(1)
    blade%nEy(nej-1)=normE(2)
    blade%nEz(nej-1)=normE(3)

    if (blade%FlipN) then
    blade%nEx(nej-1)=-blade%nEx(nej-1)
    blade%nEy(nej-1)=-blade%nEy(nej-1)
    blade%nEz(nej-1)=-blade%nEz(nej-1)
    blade%sEx(nej-1)=-blade%sEx(nej-1)
    blade%sEy(nej-1)=-blade%sEy(nej-1)
    blade%sEz(nej-1)=-blade%sEz(nej-1)
    endif
    ! Calc element area and chord
    P1=(/blade%QCx(nej-1)-0.25*blade%C(nej-1)*blade%tx(nej-1),blade%QCy(nej-1)-0.25*blade%C(nej-1)*blade%ty(nej-1),blade%QCz(nej-1)-0.25*blade%C(nej-1)*blade%tz(nej-1)/)

    P2=(/blade%QCx(nej-1)+0.75*blade%C(nej-1)*blade%tx(nej-1),blade%QCy(nej-1)+0.75*blade%C(nej-1)*blade%ty(nej-1),blade%QCz(nej-1)+0.75*blade%C(nej-1)*blade%tz(nej-1)/)

    P3=(/blade%QCx(nej)+0.75*blade%C(nej)*blade%tx(nej),blade%QCy(nej)+0.75*blade%C(nej)*blade%ty(nej),blade%QCz(nej)+0.75*blade%C(nej)*blade%tz(nej)/)

    P4=(/blade%QCx(nej)-0.25*blade%C(nej)*blade%tx(nej),blade%QCy(nej)-0.25*blade%C(nej)*blade%ty(nej),blade%QCz(nej)-0.25*blade%C(nej)*blade%tz(nej)/)

    V1=P2-P1
    V2=P3-P2
    V3=P4-P3
    V4=P1-P4
    ! Calc quad area from two triangular facets
    Call cross(V1(1),V1(2),V1(3),V2(1),V2(2),V2(3),A1(1),A1(2),A1(3))
    A1=A1/2.0
    Call cross(V3(1),V3(2),V3(3),V4(1),V4(2),V4(3),A2(1),A2(2),A2(3))
    A2=A2/2.0
    blade%EArea(nej-1)=sqrt(dot_product(A1,A1))+sqrt(dot_product(A2,A2))
    ! Calc average element chord from area and span
    blade%EC(nej-1)=blade%EArea(nej-1)/sEM
    blade%ETtoC(nej-1)=0.5*(blade%thick(nej)+blade%thick(nej-1))
    end do
    ewrite(2,*) 'Exiting make_actuatorline_geometry'

End SUBROUTINE make_actuatorline_geometry 

end module actuator_line_model
