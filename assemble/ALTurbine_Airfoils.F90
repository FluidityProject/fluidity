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

module alturbine_airfoils 

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
  
  type AirfoilType
  
  ! Naming and others
  character(len=OPTION_PATH_LEN) :: aftitle ! Title for each airfoil section 
  integer :: camb ! Camber flag for each section
  real    :: tc      ! Thickness to chord ration for each section
  real    :: alzer   ! Zero lift AOA for each section
  
  ! Airfoild section coefficient data
  real, allocatable :: TA(:,:)   ! Table AOA values
  real, allocatable :: TCL(:,:)  ! Table CL values
  real, allocatable :: TCD(:,:)  ! Table CD values
  real, allocatable :: TCM(:,:)  ! Table Cm values
  real, allocatable :: TRE(:)    ! Table Reynolds Number values  
  real, allocatable :: nTBL(:)   ! Number of AOA values for each Re number, in each section data table
  real  :: nRET   ! Number of Re number values in each section data table

  ! Interpolation warning flags
  integer :: ilxtp
  integer :: iuxtp

  ! Airfoil params for Leisham Beddoes dyn stall
  real, allocatable :: CLaData(:)
  real, allocatable :: CLCritPData(:)
  real, allocatable :: CLCritNData(:)
  
  end type AirfoilType

  type(AirfoilType), allocatable :: airfoil(:)

  ! Private subroutines
  private intp
 
  ! Public subroutines
  public airfoils_init

contains
   
    subroutine airfoils_init(airfoil_path)
    
    implicit none
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    ! This routine initialises the airfoil struct by reading 
    ! the existing input files and allocating the memory of the
    ! airfoil array
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    character(len=OPTION_PATH_LEN) :: airfoil_path

    ewrite(1,*) 'In airfoils_init'
     

    ewrite(1,*) 'Exiting airdoil_init'

    end subroutine airfoils_init
    
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    ! These routines have been taken directly out of CACTUS and changed
    ! slightly so they fit the needs of the present coupling
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

    subroutine intp(RE,ALPHA,CL,CD,CM25,airfoil)   
        
        implicit none

        real :: RE, ALPHA, CL, CD, CM25
        integer :: i,j               
        real :: XRE, XA 
        real :: CLA(2),CDA(2),CM25A(2)                                      
        type(AirfoilType) :: airfoil 
        integer :: U1, X1, iUB, iLB, NTB, L1
        logical :: NotDone                                               

        !    INTERPOLATE ON RE NO. AND ANGLE OF ATTACK TO GET AIRFOIL CHARACTERISTICS                                            
        CLA(:)=0.0                                                        
        CDA(:)=0.0  
        CM25A(:)=0.0

        if (RE >= airfoil%TRE(1)) then                                                                                                   

            ! Find Re upper and lower bounds.                                     
            NotDone=.true.    
            iUB=2                                                                 
            do while (NotDone)   

                if (RE <= airfoil%TRE(iUB)) then
                    ! Done
                    NotDone=.false.
                    if (RE == airfoil%TRE(iUB)) then
                        iLB=iUB
                    else
                        iLB=iUB-1                                                           
                        XRE=(RE-airfoil%TRE(iLB)/(airfoil%TRE(iUB)-airfoil%TRE(iLB)))
                    end if
                else
                    if (iUB == airfoil%nRET) then       
                        ! warning: no upper bound in table, take last point and set warning...
                        NotDone=.false.                                                       
                        iLB=iUB                                                           
                        XRE=0.0                                                           
                        airfoil%IUXTP=1
                    else    
                        ! No upper bound, increment and continue                                
                        iUB=iUB+1
                    end if
                end if

            end do

        else        
            ! warning: no lower bound in table, take first point and set warning                                               
            iLB=1                                                             
            iUB=1                                                             
            XRE=0.0                                                                                                 
            airfoil%ILXTP=1
        end if


        ! INTERPOLATE ON THE ANGLE OF ATTACK                               
        i=1                                                               
        do j=iLB,iUB                                                  

            NTB=airfoil%NTBL(j) ! # of alpha values in table for this section                                                
        
            ! Find upper and lower bound indicies on alpha                                                                     

            ! DO INTERVAL HALVING LOOK UP                                      

            U1=NTB                                                                                              
            L1=1                                                              
            X1=NTB/2 
            NotDone=.true. 

            do while (NotDone)                                                        
                if (ALPHA < airfoil%TA(X1,J)) then
                    U1=X1                                                                                                            
                else    
                    L1=X1 
                end if

                if ((U1-L1) == 1) then
                    NotDone=.false.
                else 
                    X1=L1+(U1-L1)/2
                end if
            end do

            ! DO STRAIGHT LINE INTERPOLATION ON ALPHA                          

            XA=(ALPHA-airfoil%TA(L1,J))/(airfoil%TA(U1,J)-airfoil%TA(L1,J))                  
            CLA(I)=airfoil%TCL(L1,J)+XA*(airfoil%TCL(U1,J)-airfoil%TCL(L1,J))                
            CDA(I)=airfoil%TCD(L1,J)+XA*(airfoil%TCD(U1,J)-airfoil%TCD(L1,J))    
            CM25A(I)=airfoil%TCM(L1,J)+XA*(airfoil%TCM(U1,J)-airfoil%TCM(L1,J)) 

            I=I+1           
        end do

        ! DO STRAIGHT LINE INTERPOLATION ON RE NO.                         

        CL=CLA(1)+XRE*(CLA(2)-CLA(1))                                     
        CD=CDA(1)+XRE*(CDA(2)-CDA(1))  
        CM25=CM25A(1)+XRE*(CM25A(2)-CM25A(1))

    END SUBROUTINE intp

end module alturbine_airfoils 
