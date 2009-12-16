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
!    C.Pain@Imperial.ac.uk
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

module solmom_module
  use FLDebug
  use solvers
  use AllSorts
  use flcomms_module
implicit none
contains    

      SUBROUTINE SOLMOM(VEL,VELGES,IGUESS,F,PHI,DT,&
     &     DMI,BIGM,FINDRM,COLM,CENTRM,NCOLM, &
     &     SAVDIA,&
     &     NONODS,&
     &     GMRESQ,&
     &     PARA,halo_tag,NNODP,&
     &     option_path)
      INTEGER NNODP
      INTEGER NONODS
      INTEGER GMRESQ
      INTEGER NCOLM
!     This sub is only called when COUPLE=0. 
!     This sub solves for velocity VEL with rhs vec F and acc PHI
!     DT=time step. 
!     If IGUESS then have a guess of velocity VEL
      INTEGER MATSTR,ZERO
      PARAMETER(MATSTR=0,ZERO=0)
      REAL VEL(NONODS),F(NONODS)
      REAL VELGES(NONODS)
      LOGICAL IGUESS
      REAL BIGM(NCOLM),PHI(NONODS),DT
!     
!     NR=11*NONODS if bi-conjigate gradient method is used. 
!     NR=10*NONODS if conjigate gradient SQUARED method is used. 
!     NR=5*NONODS if UZAWA=1 and NR=2*NONODS if multipass method is used. 
!     NR=can be unlimited if GMRES is used -it takes as much memory as it can.
      REAL SAVDIA(NONODS)
      REAL DMI(NONODS)
!     R contains temperory vectors for solver. 
!     
      INTEGER FINDRM(NONODS+1),COLM(NCOLM),CENTRM(NONODS)
!     For parallel processing ...
!     NNODP=no of subdomain nodes excluding halo nodes. 
      INTEGER  PARA,halo_tag
!     options tree path to field that we're solving for
!     this might be velocity (momentum) or any other random scalar field
!     people deemed fit to stick in solmom()
      character(len=*), intent(in):: option_path
      INTEGER I,KITS
!     
!     This sub solves a single momentum equation. 
      
      ewrite(3,*)'INSIDE SOLMOM'
      IF(IGUESS) THEN
      do I=1,NONODS
            PHI(I)=(VELGES(I)-VEL(I))/DT
      end do
      ELSE
      do I=1,NONODS
            PHI(I)=0.
      end do
      ENDIF

      do I=1,NONODS
         SAVDIA(I)=BIGM(CENTRM(I))
         BIGM(CENTRM(I))=BIGM(CENTRM(I))+DMI(I)
      end do

      ewrite(3,*) 'gmresq',gmresq
      IF(GMRESQ.NE.0) THEN
         CALL GMRES(PHI,F,NONODS,NONODS,NNODP,.FALSE.,&
     &        BIGM,FINDRM,COLM,NCOLM,NCOLM,&
     & halo_tag,KITS,&
     &        option_path=option_path)
      ENDIF

!     Put diagonal back.
      do I=1,NONODS
       BIGM(CENTRM(I))=SAVDIA(I)
      end do
!     Work out vel from accereration. 
      do I=1,NONODS
       VEL(I)= PHI(I)*DT + VEL(I)
      end do
      IF(PARA.EQ.1) CALL HALGET(VEL,NONODS,NONODS,NNODP,halo_tag)
      
end subroutine solmom
end module solmom_module
