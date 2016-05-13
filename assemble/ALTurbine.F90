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
    character(len=100) :: name
    character(len=100) :: geom_file 
    integer :: NBlades
    real, dimension(3) :: RotN, RotP ! Rotational vectors in the normal and perpendicular directions
    type(BladeType), allocatable :: Blade(:)
    
end type TurbineType

    type(TurbineType), allocatable :: Turbine(:) ! Turbine 
    integer :: notur, NBlades, NElem      ! Number of the turbines 
    integer, parameter :: NBlades_max = 4 ! Maximum Number of Blades
    integer, parameter :: NElem_max = 50  ! Maximum Number of Elements per blade
    
    private turbine_geometry_read, allocate_turbine_elements
    public  turbine_init

contains
    
    subroutine turbine_init(state)
    
    implicit none
    
    type(state_type), intent(inout) :: state
    character(len=OPTION_PATH_LEN):: turbine_path, turbine_name, bc_name
    integer :: i, j
    integer, parameter :: MaxReadLine = 1000    
    character(MaxReadLine) :: FN    ! path to geometry input file 
    integer :: NElem
    character(MaxReadLine) :: ReadLine

    ewrite(1,*) 'Hi there from the turbine init subroutine'

    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    !* This a routine from which fluidity reads the information about the turbine
    !* We start with a single turb 
    ! In the final version of the code we should be able to enable 
    ! multiple turbines. This is we we f
    ! loop through turbines
    notur = option_count("/ALM_Turbine/alm_turbine")
    ewrite(1,*) 'Number of Actuator Line Turbines : ', notur
    ! Allocate Turbines Array 
    
    call allocate_turbine_elements(notur,NBlades_max,NElem_max)

    do i=1, notur
       !turbine_path="/ALM_Turbine/alm_turbine["//int2str(i)//"]"  
       call get_option("/ALM_Turbine/alm_turbine["//int2str(i-1)//"]/name",Turbine(i)%name)
       call get_option("/ALM_Turbine/alm_turbine["//int2str(i-1)//"]/axis_location",Turbine(i)%axis_loc)
       call get_option("/ALM_Turbine/alm_turbine["//int2str(i-1)//"]/geometry_file/file_name",Turbine(i)%geom_file)
       ewrite(1,*) Turbine(i)%name
       ewrite(1,*) Turbine(i)%axis_loc
       ewrite(1,*) Turbine(i)%geom_file
       call turbine_geometry_read(i,Turbine(i)%geom_file)

   end do
   
   stop

end subroutine turbine_init
    
subroutine allocate_turbine_elements(notur,NBlades_max,NElem_max)

    implicit none

    integer :: notur,NBlades_max,NElem_max
     
    allocate(Turbine(notur)%Blade(NBlades_max)%QCx(NElem_max+1))
    allocate(Turbine(notur)%Blade(NBlades_max)%QCy(NElem_max+1))
    allocate(Turbine(notur)%Blade(NBlades_max)%QCz(NElem_max+1))
    allocate(Turbine(notur)%Blade(NBlades_max)%tx(NElem_max+1))
    allocate(Turbine(notur)%Blade(NBlades_max)%ty(NElem_max+1))
    allocate(Turbine(notur)%Blade(NBlades_max)%tz(NElem_max+1))
    allocate(Turbine(notur)%Blade(NBlades_max)%CtoR(NElem_max+1))
    allocate(Turbine(notur)%Blade(NBlades_max)%PEx(NElem_max))
    allocate(Turbine(notur)%Blade(NBlades_max)%PEy(NElem_max))
    allocate(Turbine(notur)%Blade(NBlades_max)%PEz(NElem_max))
    allocate(Turbine(notur)%Blade(NBlades_max)%tEx(NElem_max))
    allocate(Turbine(notur)%Blade(NBlades_max)%tEy(NElem_max))
    allocate(Turbine(notur)%Blade(NBlades_max)%tEz(NElem_max))
    allocate(Turbine(notur)%Blade(NBlades_max)%nEx(NElem_max))
    allocate(Turbine(notur)%Blade(NBlades_max)%nEy(NElem_max))
    allocate(Turbine(notur)%Blade(NBlades_max)%nEz(NElem_max))
    allocate(Turbine(notur)%Blade(NBlades_max)%sEx(NElem_max))
    allocate(Turbine(notur)%Blade(NBlades_max)%sEy(NElem_max))
    allocate(Turbine(notur)%Blade(NBlades_max)%sEz(NElem_max))
    allocate(Turbine(notur)%Blade(NBlades_max)%ECtoR(NElem_max))
    allocate(Turbine(notur)%Blade(NBlades_max)%EAreaR(NElem_max))
    allocate(Turbine(notur)%Blade(NBlades_max)%iSect(NElem_max))

    
end subroutine allocate_turbine_elements
    
subroutine turbine_geometry_read(i,FN)
    
    implicit none
    character(len=*) :: FN ! FileName of the geometry file
    integer :: i
    character(1000) :: ReadLine
    if(i>notur) then
        FLExit("turbine index in turbine geometry_read has exceeded the number of the turbines")
    endif
    
    open(15,file=FN)
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) Turbine(i)%NBlades
    
    if(NBlades>4) then
        FLExit("The maximum number of blades is 4")
    endif
   
    write(*,*) Turbine(i)%NBlades
    close(15)


end subroutine turbine_geometry_read



end module alturbine

