
!    Copyright (C) 2006 Imperial College London and others.
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
SUBROUTINE ADJOINT(filename, filenamelen)
  use FLDebug
  use signals
  use spud

  IMPLICIT NONE  ! ***************************************************************
  ! THIS PROGRAM SOLVES THE NAVIER STOKES, RADIATION(WITH OR 
  ! WITHOUT EVENT) AND OR ADVECTION-DIFFUSION TYPE OF EQUATIONS. 
  ! ***************************************************************
  integer, intent(in) :: filenamelen
  character(len=filenamelen), intent(in) :: filename


    if(have_option("/model/fluids/reduced/initial")) then
       
       !##############################################################
       !       Reduced Fluidity Model
       !#############################################################
       ewrite(3,*) '/model/fluids/reduced/initial'
       call initialise_signals     
       CALL ADreducedFLUIDS_Initial(filename, filenamelen)
       
    else if(have_option("/model/fluids/reduced/forward")) then
       ewrite(3,*) '/model/fluids/reduced/forward'
       !####################################################################
       !     REDUCED ADJOINT MODEL for inversion of initial condition
       !####################################################################
       
       
       CALL Reducedfluids(filename, filenamelen)
       
    else  if(have_option("/model/fluids/adjoint/sensitivity")) then
       
       ewrite(3,*) '/model/fluids/adjoint'
       !###################################################################
       !     ADJOINT MODEL for sensitivity analysis
       !####################################################################
       
       CALL ADFLUIDS_SENSITIVITY(filename, filenamelen)
       
    else if(have_option(trim("/model/fluids/adjoint/bc"))) then
       ewrite(3,*) '/model/fluids/adjoint/bc'
       
       !################################################
       !     FULL ADJOINT MODEL for boundary condition
       !################################################
       
       CALL ADFLUIDS_BC(filename, filenamelen)
       
    else if(have_option(trim("/model/fluids/adjoint/bctest"))) then
       CALL ADFLUIDS_BC(filename, filenamelen)
!       CALL ADFLUIDS_BC_TEST(filename, filenamelen)
 
    endif

  RETURN
END SUBROUTINE ADJOINT
