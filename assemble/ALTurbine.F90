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

    type(BladeType), allocatable :: Blades(:)   ! Input blade geometry

	real, allocatable :: xBE(:)     ! X location for each blade segment end (quarter chord) 	
	real, allocatable :: yBE(:)     ! Y location for each blade segment end (quarter chord)		
	real, allocatable :: zBE(:)     ! Z location for each blade segment end (quarter chord)	

    real, allocatable :: txBE(:)    ! Tangential X for each blade segment end
    real, allocatable :: tyBE(:)    ! Tangential Y for each blade segment end
    real, allocatable :: tzBE(:)    ! Tangential Z for each blade segment end

    real, allocatable :: CtoR(:)    ! Chord to radius ratio for each blade segment end

    real, allocatable :: xBC(:)     ! X location for each blade segment center (quarter chord)         
    real, allocatable :: yBC(:)     ! Y location for each blade segment center (quarter chord)         
    real, allocatable :: zBC(:)     ! Z location for each blade segment center (quarter chord) 	

    real :: dSGeom                  ! Geometry discretization level used in vortex core calculation
    real :: CrRef                   ! Ref chord to radius ratio 	
	real, allocatable :: nxBC(:)    ! Normal X for each blade segment
	real, allocatable :: nyBC(:)    ! Normal Y for each blade segment
	real, allocatable :: nzBC(:)    ! Normal Z for each blade segment 
	real, allocatable :: txBC(:)    ! Tangential X for each blade segment 	
	real, allocatable :: tyBC(:)    ! Tangential Y for each blade segment 	
	real, allocatable :: tzBC(:)    ! Tangential Z for each blade segment 	
    real, allocatable :: sxBC(:)    ! Spanwise X for each blade segment 	
	real, allocatable :: syBC(:)    ! Spanwise Y for each blade segment 	
	real, allocatable :: szBC(:)    ! Spanwise Z for each blade segment 
    ! JCM: note CircSign could be made a function of blade not element with the new geometry spec, or eliminated as it is a function of FlipN
    real, allocatable :: CircSign(:)        ! Direction of segment circulation on wake grid at positive lift
	real, allocatable :: eArea(:)           ! Element area to radius ratio for each element
	real, allocatable :: eChord(:)          ! Element chord to radius ratio for each element
    integer, allocatable :: iSect(:)        ! Array of indicies of the section table to apply to each blade element

    ! Current sum of output over all blades (nomalized by machine scale parameters)
    real :: CP_B  ! Power coefficient due to all blades
    real :: CTR_B ! Torque coefficient due to all blades
    real :: CFx_B ! Fx coefficient due to all blades
    real :: CFy_B ! Fy coefficient due to all blades
    real :: CFz_B ! Fz coefficient due to all blades


contains
    
    subroutine turbine_init(state)
     implicit none
     type(state_type), intent(inout) :: state
     ewrite(1,*) 'Hi there from the turbine init subroutine'
     !* This a routine from which fluidity reads the information about the turbine
     !* We start with a single turb
    
     call read_geometry

    end subroutine turbine_init
    
    subroutine read_geometry
    
    implicit none
    
    character(len=OPTION_PATH_LEN):: turbine_path, turbine_name, bc_name
    integer :: notur, i, j
    integer, parameter :: MaxReadLine = 1000    
    character(MaxReadLine) :: FN    ! path to geometry input file 
    integer :: NElem
    character(MaxReadLine) :: ReadLine


    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    ! In the final version of the code we should be able to enable 
    ! multiple turbines. This is we we f
    ! loop through turbines
    notur = option_count("/ALM_Turbine/turbine")
    do i=0, notur-1
       turbine_path="/ALM_Turbine/turbine["//int2str(i)//"]"
       write(*,*) notur
       stop
    
    end do
     
    
    end subroutine read_geometry


end module alturbine

