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

module alturbine

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
  use alturbine_utils
  use lbdynstall

  implicit none
!GGGGG
! Define the types that will be used

type BladeType
    integer :: NElem                ! Number of Elements of the Blade
    integer :: FlipN                ! 
    real, allocatable :: QCx(:)     ! Blade quarter-chord line x coordinates at element ends
    real, allocatable :: QCy(:)     ! Blade quarter-chord line y coordinates at element ends
    real, allocatable :: QCz(:)     ! Blade quarter-chord line z coordinates at element ends
    real, allocatable :: tx(:)      ! Blade unit tangent vector (rearward chord line direction) x-componenst at element ends
    real, allocatable :: ty(:)      ! Blade unit tangent vector (rearward chord line direction) y-componenst at element ends 
    real, allocatable :: tz(:)      ! Blade unit tangent vector (rearward chord line direction) z-componenst at element ends  
    real, allocatable :: CtoR(:)    ! Blade chord to turbine reference radius ratio at element ends
    real, allocatable :: PEx(:)     ! Element centre x coordinates
    real, allocatable :: PEy(:)     ! Element centre y coordinates
    real, allocatable :: PEz(:)     ! Element centre z coordinates
    real, allocatable :: tEx(:)     ! Element unit tangent vector (rearward chord line direction) x-component
    real, allocatable :: tEy(:)     ! Element unit tangent vector (rearward chord line direction) y-component
    real, allocatable :: tEz(:)     ! Element unit tangent vector (rearward chord line direction) z-component
    real, allocatable :: nEx(:)     ! Element unit normal vector x-component
    real, allocatable :: nEy(:)     ! Element unit normal vector y-component
    real, allocatable :: nEz(:)     ! Element unit normal vector z-component
    real, allocatable :: sEx(:)     ! Element unit spanwise vector x-component 
    real, allocatable :: sEy(:)     ! Element unit spanwise vector y-component
    real, allocatable :: sEz(:)     ! Element unit spanwise vector z-component
    real, allocatable :: EC(:)      ! Element chord lenght
    real, allocatable :: CircSign(:)! Direction of segment circulation on wake
    real, allocatable :: EArea(:)   ! Element Area
    real, allocatable :: ETtoC(:)! Element thickness to Chord ratio
    
    type(AirfoilType), allocatable :: EAirfoil(:) ! Element Airfoil 
    type(LB_Type), allocatable :: E_LB_Model(:)   ! Element Leishman-Beddoes Model

    ! Momentum Sink Forces in the nts direction
    real, allocatable :: CFn(:)     ! Element Force in the normal direction
    real, allocatable :: CFt(:)     ! Element Force in the tangential direction (rearward chord line direction) 
    real, allocatable :: CFs(:)     ! Element Force in the spanwise direction

    ! Momentum Sink Forces in the xyz direction
    real, allocatable :: CFx(:)     ! Element Force in the global x-direction
    real, allocatable :: CFy(:)     ! Element Force in the global y-direction
    real, allocatable :: CFz(:)     ! Element Force in the global z-direction
    
end type BladeType

type TurbineType

    real,dimension(3) :: axis_loc
    character(len=100) :: turb_name
    character(len=100) :: geom_file 
    integer :: NBlades, NAirfoilData
    real, dimension(3) :: RotN, RotP ! Rotational vectors in the normal and perpendicular directions
    real :: at, Rmax
    real :: RPM 
    logical :: Is_constant_rotation_operated = .false. ! For a constant rotational velocity (in Revolutions Per Minute)
    logical :: Is_force_based_operated = .false. ! For a forced based rotational velocity (Computed during the simulation)
    type(BladeType), allocatable :: Blade(:)
    type(AirfoilType), allocatable :: AirfoilData(:)
    real :: CP  ! Power coefficient 
    real :: CTR ! Torque coefficient 
    real :: CFx ! Fx coefficient 
    real :: CFy ! Fy coefficient 
    real :: CFz ! Fz coefficient 
    real :: CT  ! Thrust coefficient

end type TurbineType

    type(TurbineType), allocatable :: Turbine(:) ! Turbine 
    integer :: notur, NBlades, NElem      ! Number of the turbines 
    type(vector_field), pointer :: TurbineSource ! Pointer to the turbine source term
    
    private turbine_geometry_read, allocate_turbine_elements, allocate_turbine_blades 
    public  turbine_init, turbine_timeloop

contains
    
    subroutine turbine_init(state)
    
    implicit none
    
    type(state_type) :: state
    character(len=OPTION_PATH_LEN)::  turbine_name
    integer :: i, j,k
    integer, parameter :: MaxReadLine = 1000    
    character(MaxReadLine) :: FN    ! path to geometry input file 
    integer :: NElem
    character(MaxReadLine) :: ReadLine
    character(len=OPTION_PATH_LEN), allocatable :: turbine_path(:)
    character(len=OPTION_PATH_LEN) :: section_path
    ewrite(1,*) 'Entering the ALTurbine_init '

    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    !* This a routine from which fluidity reads the information about the turbine model 
    !* We start with a single turbine model 
    ! In the final version of the code we should be able to enable 
    ! multiple turbines. 
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

    notur = option_count("/ALM_Turbine/alm_turbine")
    ewrite(2,*) 'Number of Actuator Line Turbines : ', notur
    ! Allocate Turbines Array 
    Allocate(Turbine(notur))
    Allocate(turbine_path(notur))
    
    
    do i=1, notur
       
       turbine_path(i)="/ALM_Turbine/alm_turbine["//int2str(i-1)//"]"
       call get_option("/ALM_Turbine/alm_turbine["//int2str(i-1)//"]/name",Turbine(i)%turb_name)
       call get_option("/ALM_Turbine/alm_turbine["//int2str(i-1)//"]/geometry_file/file_name",Turbine(i)%geom_file)
       
       ! Count how many Airfoil Sections are available
       Turbine(i)%NAirfoilData=option_count("/ALM_Turbine/alm_turbine["//int2str(i-1)//"]/airfoil_sections/section") 
       ewrite(2,*) 'Number of Airfoils available : ', Turbine(i)%NAirfoilData
       ! Allocate the memory of the Airfoils
       Allocate(Turbine(i)%AirfoilData(Turbine(i)%NAirfoilData))
       
       call turbine_geometry_read(i,Turbine(i)%geom_file) 
       ewrite(2,*) 'Turbine ',i,' : ', Turbine(i)%turb_name
       ewrite(2,*) '---------------'
       ewrite(2,*) 'Axis location : ',Turbine(i)%RotP
       ewrite(2,*) 'Geometry file :    ',Turbine(i)%geom_file
       
       do k=1, Turbine(i)%NAirfoilData
           
        call get_option(trim(turbine_path(i))//"/airfoil_sections/section["//int2str(k-1)//"]/airfoil_file",Turbine(i)%AirfoilData(k)%afname)
           
           ! Read and Store Airfoils
           call airfoil_init_data(Turbine(i)%AirfoilData(k))
            

       end do

       do j=1,Turbine(i)%NBlades 
       ! Populate Turbine Airfoil Sections
        ewrite(2,*) 'Populating blade airfoils for turbine ',i,' Blade ', j
        call populate_blade_airfoils(Turbine(i)%Blade(j)%NElem,Turbine(i)%Blade(j)%EAirfoil,Turbine(i)%AirfoilData,Turbine(i)%Blade(j)%ETtoC)
       enddo

       ! Check the type of Turbine Operation
       if (have_option(trim(turbine_path(i))//"/operation/constant_rotational_velocity")) then
            Turbine(i)%Is_constant_rotation_operated= .true.
            call get_option("/ALM_Turbine/alm_turbine["//int2str(i-1)//"]/operation/constant_rotational_velocity/RPM",Turbine(i)%RPM)
       else if(have_option(trim("/ALM_Turbine/alm_turbine["//int2str(i-1)//"]")//"/operation/force_based_rotational_velocity")) then
            Turbine(i)%Is_force_based_operated = .true. 
       else
           FLExit("At the moment only the constant and the force_based rotational velocity models are supported") 
       endif
       
       ! Specify the operation mode of the turbine
       if(Turbine(i)%Is_constant_rotation_operated) then
       ewrite(2,*) 'Constant rotational velocity : ', Turbine(i)%Is_constant_rotation_operated 
       ewrite(2,*) 'RPM : ', Turbine(i)%RPM
       else
       ewrite(2,*) 'Forced-based rotational velocity : ', Turbine(i)%Is_force_based_operated
       endif
       ewrite(2,*) ' '
       
   
   end do
  
   ! Creating the Source Term that will allow to model the Turbine input
   !TurbineSource => extract_vector_field(states,"VelocitySource")
   ! Zeroing the Source Term (This should be done at each time step)
   !call zero(TurbineSource)
   ewrite(1,*) 'Exiting the ALTurbine_init'

end subroutine turbine_init

subroutine turbine_timeloop(state,exclude_nonrecalculated)
 
    implicit none
    
    type(state_type), intent(inout) :: state
    logical, intent(in), optional :: exclude_nonrecalculated
    type(vector_field), pointer :: positions, velocity
    integer :: i,j,k
    real :: dt, theta
    ! Zero the Source Term at each time step
    
    ewrite(1,*) 'Entering the turbine_timeloop'
   
    
    ! First we need to get the dt in order to rotate the turbine
    call get_option("/timestepping/timestep",dt) 
    
    do i=1,notur
    ! Depending on whether the turbine is using a constant or a forced
    ! Based model for its operation: we will have the following options
     if(Turbine(i)%Is_constant_rotation_operated) then
        ewrite(2,*) 'Operating Turbine with a constant rotational Velocity'
        theta=Turbine(i)%RPM*2*pi/60*dt ! 1 revolution/minute = 2 pi tads / 60 s
        call rotate_turbines(theta) 

    elseif(Turbine(i)%Is_force_based_operated) then
        ewrite(2,*) 'Operating Turbine with a force-based approach'
     else
         FLExit("At the moment only constant rotational velocity and the force-based approach are supported") 
     endif
    
    end do

    ewrite(1,*) 'Exiting the turbine_timeloop'
    return

end subroutine turbine_timeloop

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

subroutine rotate_turbines(theta)

    implicit none

    real :: theta,nrx,nry,nrz,px,py,pz 
    integer :: j,ielem,i
    real :: vrx,vry,vrz,VMag
    real :: xtmp,ytmp,ztmp, txtmp, tytmp, tztmp
    ! Rotates data in blade arrays. Rotate element end geometry and recalculate element geometry.

    ewrite(1,*) 'Entering rotate_turbines'
    do i=1,notur
        
        ! Specify the rotation axis and the normal vector of rotation
        
        nrx=Turbine(i)%RotN(1)
        nry=Turbine(i)%RotN(2)
        nrz=Turbine(i)%RotN(3)
        
        px=Turbine(i)%RotP(1)
        py=Turbine(i)%RotP(2)
        pz=Turbine(i)%RotP(3)


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
            
            call set_blade_geometry(Turbine(i)%Blade(j))
        end do 
    end do
    
    ewrite(1,*) 'Exiting rotate_turbines'

end subroutine rotate_turbines
    
subroutine allocate_turbine_blades(ITurbine,NBlades)
    implicit none

    integer :: ITurbine,NBlades

    allocate(Turbine(ITurbine)%Blade(NBlades))
    
end subroutine allocate_turbine_blades

subroutine allocate_turbine_elements(ITurbine,IBlade,NElem)

    implicit none

    integer :: ITurbine,IBlade,NElem
     
    allocate(Turbine(ITurbine)%Blade(IBlade)%QCx(NElem+1))
    allocate(Turbine(ITurbine)%Blade(IBlade)%QCy(NElem+1))
    allocate(Turbine(ITurbine)%Blade(IBlade)%QCz(NElem+1))
    allocate(Turbine(ITurbine)%Blade(IBlade)%tx(NElem+1))
    allocate(Turbine(ITurbine)%Blade(IBlade)%ty(NElem+1))
    allocate(Turbine(ITurbine)%Blade(IBlade)%tz(NElem+1))
    allocate(Turbine(ITurbine)%Blade(IBlade)%CtoR(NElem+1))
    allocate(Turbine(ITurbine)%Blade(IBlade)%PEx(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%PEy(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%PEz(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%tEx(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%tEy(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%tEz(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%nEx(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%nEy(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%nEz(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%sEx(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%sEy(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%sEz(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%EC(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%EArea(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%CircSign(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%ETtoC(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%EAirfoil(Nelem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%E_LB_Model(Nelem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%CFn(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%CFt(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%CFs(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%CFx(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%CFy(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%CFz(NElem))

    
end subroutine allocate_turbine_elements
    
subroutine turbine_geometry_read(i,FN)
    
    implicit none
    character(len=*) :: FN ! FileName of the geometry file
    integer :: i,j
    character(1000) :: ReadLine
    
    if(i>notur) then
        FLExit("turbine index in turbine_geometry_read has exceeded the number of the turbines")
    endif
    
    open(15,file=FN)
    ! Read the Number of Blades
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%NBlades
     
    call allocate_turbine_blades(i,Turbine(i)%NBlades)
    
    ! Read the Turbine rotation axis normal vector (x y z values)
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%RotN(1), Turbine(i)%RotN(2), Turbine(i)%RotN(3)
    
    ! Read the Turbine rotation origin point (x y z values)
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%RotP(1), Turbine(i)%RotP(2), Turbine(i)%RotP(3)
    
    ! Read the perpendicular Axis of Rotation
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%at

    ! Read the perpendicular Axis of Rotation
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Rmax 

    read(15,'(A)') ReadLine ! Defines the Type of the Turbine

    
    do j=1,Turbine(i)%NBlades
    
        read(15,'(A)') ReadLine ! Blade ....
    
        !Read Number of Elements
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%NElem

        call allocate_turbine_elements(i,j,Turbine(i)%Blade(j)%NElem)
        
        !Read FlipN
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%FlipN

        !Read QCx(1:NElem+1)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%QCx(1:Turbine(i)%Blade(j)%Nelem+1)
        
        !Read QCy(1:NElem+1)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%QCy(1:Turbine(i)%Blade(j)%Nelem+1)
        
        !Read QCz(1:NElem+1)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%QCz(1:Turbine(i)%Blade(j)%Nelem+1)
        
        !Read tx(1:NElem+1)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%tx(1:Turbine(i)%Blade(j)%Nelem+1)
        
        !Read ty(1:NElem+1)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%ty(1:Turbine(i)%Blade(j)%Nelem+1)
        
        !Read tz(1:NElem+1)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%tz(1:Turbine(i)%Blade(j)%Nelem+1)
      
        !Read CtoR(1:NElem+1)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%CtoR(1:Turbine(i)%Blade(j)%Nelem+1)
        
        !Read PEx(1:NElem)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%PEx(1:Turbine(i)%Blade(j)%Nelem)

        !Read PEy(1:NElem)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%PEy(1:Turbine(i)%Blade(j)%Nelem)
        
        !Read PEz(1:NElem)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%PEz(1:Turbine(i)%Blade(j)%Nelem)
        
        !Read tEx(1:NElem)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%tEx(1:Turbine(i)%Blade(j)%Nelem)
        
        !Read tEy(1:NElem)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%tEy(1:Turbine(i)%Blade(j)%Nelem)
        
        !Read tEz(1:NElem)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%tEz(1:Turbine(i)%Blade(j)%Nelem)
        
        !Read nEx(1:NElem)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%nEx(1:Turbine(i)%Blade(j)%Nelem)
        
        !Read nEy(1:NElem)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%nEy(1:Turbine(i)%Blade(j)%Nelem)
        
        !Read nEz(1:NElem)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%nEz(1:Turbine(i)%Blade(j)%Nelem)
        
        !Read sEx(1:NElem)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%sEx(1:Turbine(i)%Blade(j)%Nelem)
        
        !Read sEy(1:NElem)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%sEy(1:Turbine(i)%Blade(j)%Nelem)
        
        !Read sEz(1:NElem)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%sEz(1:Turbine(i)%Blade(j)%Nelem)
        
        !Read ECtoR(1:NElem)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%EC(1:Turbine(i)%Blade(j)%Nelem)
        
        !Read EAreaR(1:NElem
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%EArea(1:Turbine(i)%Blade(j)%Nelem)
        
        !Read iSect(1:NElem)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%ETtoC(1:Turbine(i)%Blade(j)%Nelem)
        
    end do
    
    close(15)


end subroutine turbine_geometry_read

    SUBROUTINE set_blade_geometry(blade)

        implicit none
        
        type(BladeType),intent(INOUT) :: blade
        integer :: BNum
        integer :: nbe, nei, FlipN, nej, j
        real :: sEM, tEM, nEM
        real :: PE(3), sE(3), tE(3), normE(3), P1(3), P2(3), P3(3), P4(3), V1(3), V2(3), V3(3), V4(3), A1(3), A2(3)
    
    ewrite(2,*) 'Entering set_blade_geometry'
    ! Calculates element geometry from element end geometry

    ! JCM: Eventually, should just be able to loop through Blades(BNum) data structure
    ! While data is still held in arrays concatenated across blades, need to replicate
    ! nbe (stored in configr) from Blades(1).NElem
        nbe=blade%NElem

        FlipN=blade%FlipN
        
        do j=1,nbe
            nej=1+j

            ! Element center locations
            blade%PEx(nej-1)=(blade%QCx(nej)+blade%QCx(nej-1))/2.0
            blade%PEy(nej-1)=(blade%QCy(nej)+blade%QCy(nej-1))/2.0
            blade%PEz(nej-1)=(blade%QCz(nej)+blade%QCz(nej-1))/2.0

  ! Set spanwise and tangential vectors
   sE=-(/blade%QCx(nej)-blade%QCx(nej-1),blade%QCy(nej)-blade%QCy(nej-1),blade%QCz(nej)-blade%QCz(nej-1)/) ! nominal element spanwise direction set opposite to QC line
		    sEM=sqrt(dot_product(sE,sE))
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

            ! Flip normal direction if requested
		    blade%CircSign(nej)=1.0
		    if (FlipN .eq. 1) then
		        blade%nEx(nej-1)= -blade%nEx(nej-1)
		        blade%nEy(nej-1)= -blade%nEy(nej-1)
		        blade%nEz(nej-1)= -blade%nEz(nej-1)
		        blade%sEx(nej-1)= -blade%sEx(nej-1)
		        blade%sEy(nej-1)= -blade%sEy(nej-1)
		        blade%sEz(nej-1)= -blade%sEz(nej-1)
		        blade%tEx(nej-1)= -blade%tEx(nej-1)
		        blade%tEy(nej-1)= -blade%tEy(nej-1)
		        blade%tEz(nej-1)= -blade%tEz(nej-1)
                blade%CircSign(nej-1)=-1.0
		    end if

		    ! Calc element area and chord
    P1=(/blade%QCx(nej-1)-0.25*blade%CtoR(nej-1)*blade%tx(nej-1),blade%QCy(nej-1)-0.25*blade%CtoR(nej-1)*blade%ty(nej-1),blade%QCz(nej-1)-0.25*blade%CtoR(nej-1)*blade%tz(nej-1)/)

    P2=(/blade%QCx(nej-1)+0.75*blade%CtoR(nej-1)*blade%tx(nej-1),blade%QCy(nej-1)+0.75*blade%CtoR(nej-1)*blade%ty(nej-1),blade%QCz(nej-1)+0.75*blade%CtoR(nej-1)*blade%tz(nej-1)/)

    P3=(/blade%QCx(nej)+0.75*blade%CtoR(nej)*blade%tx(nej),blade%QCy(nej)+0.75*blade%CtoR(nej)*blade%ty(nej),blade%QCz(nej)+0.75*blade%CtoR(nej)*blade%tz(nej)/)

    P4=(/blade%QCx(nej)-0.25*blade%CtoR(nej)*blade%tx(nej),blade%QCy(nej)-0.25*blade%CtoR(nej)*blade%ty(nej),blade%QCz(nej)-0.25*blade%CtoR(nej)*blade%tz(nej)/)

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

        end do
        ewrite(2,*) 'Exiting set_blade_geometry'

    End SUBROUTINE set_blade_geometry 

end module alturbine

