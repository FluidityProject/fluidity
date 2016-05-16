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
  use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, PYTHON_FUNC_LEN
  use sparse_tools
  use fetools
  use fields
  use futils
  use elements
  use transform_elements
  use state_module
  use boundary_conditions
  use solvers
  use python_state
  use sparsity_patterns_meshes
  use field_options
  use fefields

  !use the ALTurbine_Airfoil Module
  use alturbine_airfoils
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

    ! Current total blade output (nomalized by machine scale parameters)
    real :: CP  ! Power coefficient due to this blade
    real :: CTR ! Torque coefficient due to this blade
    real :: CFx ! Fx coefficient due to this blade
    real :: CFy ! Fy coefficient due to this blade
    real :: CFz ! Fz coefficient due to this blade
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
    logical :: Is_forced_based_operated = .false. ! For a forced based rotational velocity (Computed during the simulation)
    type(BladeType), allocatable :: Blade(:)
    type(AirfoilType), allocatable :: Airfoil(:)
    
end type TurbineType

    type(TurbineType), allocatable :: Turbine(:) ! Turbine 
    integer :: notur, NBlades, NElem      ! Number of the turbines 
    
    private turbine_geometry_read, allocate_turbine_elements, allocate_turbine_blades, set_turbine_location
    public  turbine_init

contains
    
    subroutine turbine_init(state)
    
    implicit none
    
    type(state_type), intent(inout) :: state
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
       ewrite(2,*) 'Number of airfoils : ', Turbine(i)%NAirfoils
       ! Allocate the memory of the Airfoils
       Allocate(Turbine(i)%Airfoil(Turbine(i)%NAirfoils))
       
       do k=1, Turbine(i)%NAirfoils
           
        call get_option(trim(turbine_path(i))//"/airfoil_sections/section["//int2str(k-1)//"]/airfoil_file",Turbine(i)%Airfoil(k)%aftitle)
           
           ! Read and Store Airfoils
           call airfoil_init(Turbine(i)%Airfoil(k))

       end do

       ! Check the type of Turbine Operation
       if (have_option(trim(turbine_path(i))//"/operation/constant_rotational_velocity")) then
            Turbine(i)%Is_constant_rotation_operated= .true.
            call get_option("/ALM_Turbine/alm_turbine["//int2str(i-1)//"]/operation/constant_rotational_velocity/RPM",Turbine(i)%RPM)
       else if(have_option(trim("/ALM_Turbine/alm_turbine["//int2str(i-1)//"]")//"/operation/force_based_rotational_velocity")) then
            Turbine(i)%Is_forced_based_operated = .true. 
       else
           FLExit("At the moment only the constant and the force_based rotational velocity models are supported") 
       endif
       
       call turbine_geometry_read(i,Turbine(i)%geom_file) 
       ewrite(2,*) 'Turbine ',i,' : ', Turbine(i)%turb_name
       ewrite(2,*) '---------------'
       ewrite(2,*) 'Axis location : ',Turbine(i)%RotP
       ewrite(2,*) 'Geometry file :    ',Turbine(i)%geom_file
       if(Turbine(i)%Is_constant_rotation_operated) then
       ewrite(2,*) 'Constant rotational velocity : ', Turbine(i)%Is_constant_rotation_operated 
       ewrite(2,*) 'RPM : ', Turbine(i)%RPM
       else
       ewrite(2,*) 'Forced-based rotational velocity : ', Turbine(i)%Is_forced_based_operated
       endif
       ewrite(2,*) ' '
       
       call set_turbine_location(i,Turbine(i)%axis_loc)
   
   end do
    
   ewrite(1,*) 'Exiting the ALTurbine_init'
stop

end subroutine turbine_init

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

    subroutine set_turbine_location(i,axis_loc)
        implicit none
        integer :: i
        real,dimension(3) :: axis_loc
        
        !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
        ! This Routine translates the turbine axis (and all the turbine) from the (0,0,0)
        ! to the given axis_location of 
        !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
             
         

    end subroutine set_turbine_location

end module alturbine

