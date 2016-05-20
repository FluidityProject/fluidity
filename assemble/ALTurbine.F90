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

  !use the ALTurbine_Airfoil Module
  use alturbine_airfoils
  use alturbine_utils
  implicit none
!GGGGG
! Define the types that will be used


type BladeType
    integer :: NElem
    integer :: FlipN
    real, allocatable :: QCx(:)
    real, allocatable :: QCy(:)
    real, allocatable :: QCz(:)
    real, allocatable :: tx(:)
    real, allocatable :: ty(:)
    real, allocatable :: tz(:)
    real, allocatable :: CtoR(:)
    real, allocatable :: PEx(:)
    real, allocatable :: PEy(:)
    real, allocatable :: PEz(:)
    real, allocatable :: tEx(:)
    real, allocatable :: tEy(:)
    real, allocatable :: tEz(:)
    real, allocatable :: nEx(:)
    real, allocatable :: nEy(:)
    real, allocatable :: nEz(:)
    real, allocatable :: sEx(:)
    real, allocatable :: sEy(:)
    real, allocatable :: sEz(:)
    real, allocatable :: ECtoR(:)
    real, allocatable :: EAreaR(:)
    integer, allocatable :: iSect(:)
    
    ! Momentum Sink Forces in the nts direction
    real, allocatable :: CFn(:)
    real, allocatable :: CFt(:)
    real, allocatable :: CFs(:)

    ! Momentum Sink Forces in the xyz direction
    real, allocatable :: CFx(:)
    real, allocatable :: CFy(:)
    real, allocatable :: CFz(:)
    
end type BladeType

type TurbineType

    real,dimension(3) :: axis_loc
    character(len=100) :: turb_name
    character(len=100) :: geom_file 
    integer :: NBlades, NAirfoils
    real, dimension(3) :: RotN, RotP ! Rotational vectors in the normal and perpendicular directions
    real :: at, Rmax
    real :: RPM 
    logical :: Is_constant_rotation_operated = .false. ! For a constant rotational velocity (in Revolutions Per Minute)
    logical :: Is_force_based_operated = .false. ! For a forced based rotational velocity (Computed during the simulation)
    type(BladeType), allocatable :: Blade(:)
    type(AirfoilType), allocatable :: Airfoil(:)
    real :: CP  ! Power coefficient 
    real :: CTR ! Torque coefficient 
    real :: CFx ! Fx coefficient 
    real :: CFy ! Fy coefficient 
    real :: CFz ! Fz coefficient 
    real :: CT  ! Thrust coefficient

end type TurbineType

    type(TurbineType), allocatable :: Turbine(:) ! Turbine 
    integer :: notur, NBlades, NElem      ! Number of the turbines 
    
    private turbine_geometry_read, allocate_turbine_elements, allocate_turbine_blades 
    public  turbine_init, turbine_timeloop

contains
    
    subroutine turbine_init(states)
    
    implicit none
    
    type(state_type), dimension(:) :: states
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
       Turbine(i)%NAirfoils=option_count("/ALM_Turbine/alm_turbine["//int2str(i-1)//"]/airfoil_sections/section") 
       ewrite(2,*) 'Number of Airfoils available : ', Turbine(i)%NAirfoils
       ! Allocate the memory of the Airfoils
       Allocate(Turbine(i)%Airfoil(Turbine(i)%NAirfoils))
       
       call turbine_geometry_read(i,Turbine(i)%geom_file) 
       ewrite(2,*) 'Turbine ',i,' : ', Turbine(i)%turb_name
       ewrite(2,*) '---------------'
       ewrite(2,*) 'Axis location : ',Turbine(i)%RotP
       ewrite(2,*) 'Geometry file :    ',Turbine(i)%geom_file
       
       do k=1, Turbine(i)%NAirfoils
           
        call get_option(trim(turbine_path(i))//"/airfoil_sections/section["//int2str(k-1)//"]/airfoil_file",Turbine(i)%Airfoil(k)%afname)
           
           ! Read and Store Airfoils
           call airfoil_init_data(Turbine(i)%Airfoil(k))

       end do

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
    
   ewrite(1,*) 'Exiting the ALTurbine_init'

end subroutine turbine_init

subroutine turbine_timeloop(states,exclude_nonrecalculated)

    use pickers_inquire
    
    implicit none
    
    type(state_type), dimension(:), intent(inout) :: states
    logical, intent(in), optional :: exclude_nonrecalculated
    type(vector_field), pointer :: positions, velocity
    type(vector_field) :: v_field
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

   
    !* Finds the element number end local coordinates of the element where the point of interest
    !* belongs in.
    
    !do j=1,Turbine(i)%NBlade
    !    k=1,Turbine(i)%Blade(j)%NElem+1
    !    SCoords(1)=
    !call picker_inquire(positions,Turbine(i)%Blade(1:Turbine(i)%NBlades)%,ele,local_coord,.false.)
    ! 
    !* Evaluates the velocity at the point of interest 
    !value_vel = eval_field(ele,velocity, local_coord) 
    !
    !!* Temporary Parameters that will be deleted and normaly be loaded
    !!* From a file.
    !
    !radius= 0.25
    !epsilon_par=0.05
    !Area = pi*radius**2 ! Termporary computation of the Area 
    !Force_coeff(1)=-1.1
    !Force_coeff(2)=0

    !do i = 1, node_count(v_field)
    !    ! Compute a molification function in 2D 
    !    Rcoords=node_val(positions,i)
    !        dr2=0
    !        d=0
    !        do j=1, v_field%dim
    !            d(j) = Scoords(j)-Rcoords(j)
    !            dr2=dr2+d(j)**2
    !        end do

    !        do j=1, v_field%dim
    !        
    !        if(v_field%dim==2) then    
    !            kernel(j) = 1/(epsilon_par**2*pi)*exp(-dr2/epsilon_par**2)
    !        else 
    !            FLAbort("3D source not implemented yet")
    !        endif

    !        end do
   
    !    call set(v_field, i, 0.5*Area*(value_vel(1)**2+value_vel(2)**2)*kernel*Force_coeff)
    !end do

    end do

    ewrite(1,*) 'Exiting the turbine_timeloop'

end subroutine turbine_timeloop


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
    allocate(Turbine(ITurbine)%Blade(IBlade)%ECtoR(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%EAreaR(NElem))
    allocate(Turbine(ITurbine)%Blade(IBlade)%iSect(NElem))
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
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%ECtoR(1:Turbine(i)%Blade(j)%Nelem)
        
        !Read EAreaR(1:NElem+1)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%EAreaR(1:Turbine(i)%Blade(j)%Nelem)
        
        !Read iSect(1:NElem+1)
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%Blade(j)%iSect(1:Turbine(i)%Blade(j)%Nelem)
        
    end do
    
    close(15)


end subroutine turbine_geometry_read

    SUBROUTINE set_blade_geometry

        implicit none

        integer :: BNum
        integer :: nbe, nei, FlipN, nej, j
    !    real :: sEM, tEM, nEM
    !    real :: PE(3), sE(3), tE(3), normE(3), P1(3), P2(3), P3(3), P4(3), V1(3), V2(3), V3(3), V4(3), A1(3), A2(3)

    !    ! Calculates element geometry from element end geometry

    !    ! JCM: Eventually, should just be able to loop through Blades(BNum) data structure
    !    ! While data is still held in arrays concatenated across blades, need to replicate
    !    ! nbe (stored in configr) from Blades(1).NElem
    !    nbe=Blades(1)%NElem

    !    FlipN=Blades(BNum)%FlipN

    !    nei=1+(BNum-1)*(nbe+1)

    !    do j=1,nbe
    !        nej=nei+j

    !        ! Element center locations
    !        xBC(nej)=(xBE(nej)+xBE(nej-1))/2.0
    !        yBC(nej)=(yBE(nej)+yBE(nej-1))/2.0
    !        zBC(nej)=(zBE(nej)+zBE(nej-1))/2.0

    !        ! Set spanwise and tangential vectors
	!	    sE=-(/xBE(nej)-xBE(nej-1),yBE(nej)-yBE(nej-1),zBE(nej)-zBE(nej-1)/) ! nominal element spanwise direction set opposite to QC line
	!	    sEM=sqrt(dot_product(sE,sE))
	!	    sE=sE/sEM
	!	    tE=(/txBE(nej)+txBE(nej-1),tyBE(nej)+tyBE(nej-1),tzBE(nej)+tzBE(nej-1)/)/2.0
	!	    ! Force tE normal to sE
	!	    tE=tE-dot_product(tE,sE)*sE
	!	    tEM=sqrt(dot_product(tE,tE))
	!	    tE=tE/tEM
	!	    sxBC(nej)=sE(1)
	!	    syBC(nej)=sE(2)
	!	    szBC(nej)=sE(3)
    !        txBC(nej)=tE(1)
    !        tyBC(nej)=tE(2)
    !        tzBC(nej)=tE(3)

	!	    ! Calc normal vector
	!	    Call cross(sE(1),sE(2),sE(3),tE(1),tE(2),tE(3),normE(1),normE(2),normE(3))
	!	    nEM=sqrt(dot_product(normE,normE))
	!	    normE=normE/nEM
	!	    nxBC(nej)=normE(1)
	!	    nyBC(nej)=normE(2)
	!	    nzBC(nej)=normE(3)

	!	    ! Flip normal direction if requested
	!	    CircSign(nej)=1.0
	!	    if (FlipN .eq. 1) then
	!            nxBC(nej)=-nxBC(nej)
	!            nyBC(nej)=-nyBC(nej)
	!            nzBC(nej)=-nzBC(nej)
    !            sxBC(nej)=-sxBC(nej)
    !            syBC(nej)=-syBC(nej)
    !            szBC(nej)=-szBC(nej)
    !            CircSign(nej)=-1.0
	!	    end if

	!	    ! Calc element area and chord
	!	    P1=(/xBE(nej-1)-0.25*CtoR(nej-1)*txBE(nej-1),yBE(nej-1)-0.25*CtoR(nej-1)*tyBE(nej-1),zBE(nej-1)-0.25*CtoR(nej-1)*tzBE(nej-1)/)
	!	    P2=(/xBE(nej-1)+0.75*CtoR(nej-1)*txBE(nej-1),yBE(nej-1)+0.75*CtoR(nej-1)*tyBE(nej-1),zBE(nej-1)+0.75*CtoR(nej-1)*tzBE(nej-1)/)
	!	    P3=(/xBE(nej)+0.75*CtoR(nej)*txBE(nej),yBE(nej)+0.75*CtoR(nej)*tyBE(nej),zBE(nej)+0.75*CtoR(nej)*tzBE(nej)/)
	!	    P4=(/xBE(nej)-0.25*CtoR(nej)*txBE(nej),yBE(nej)-0.25*CtoR(nej)*tyBE(nej),zBE(nej)-0.25*CtoR(nej)*tzBE(nej)/)
	!	    V1=P2-P1
	!	    V2=P3-P2
	!	    V3=P4-P3
	!	    V4=P1-P4
	!	    ! Calc quad area from two triangular facets
	!	    Call cross(V1(1),V1(2),V1(3),V2(1),V2(2),V2(3),A1(1),A1(2),A1(3))
	!	    A1=A1/2.0
    !        Call cross(V3(1),V3(2),V3(3),V4(1),V4(2),V4(3),A2(1),A2(2),A2(3))
    !        A2=A2/2.0
	!	    eArea(nej)=sqrt(dot_product(A1,A1))+sqrt(dot_product(A2,A2))
	!	    ! Calc average element chord from area and span
	!	    eChord(nej)=eArea(nej)/sEM

    !    end do
    End SUBROUTINE set_blade_geometry 

end module alturbine

