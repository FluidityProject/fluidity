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
!
!
!
!
module limit_gradient
  implicit none

  private
  public :: limgra

contains

  SUBROUTINE LIMGRA(D3,NGI,NLOC,NONODS,TOTELE, &
      &             ELE,GI,&
      &             NDGLNO, UD,VD,WD, U, HOVERQ, &
      &             N,NX,NY,NZ, RESTRI,&
      ! 2-d Inverse of Jacobian...
      &             DETJ, AGI,BGI,&
      &             CGI,DGI,&
      ! 3-d inverse of Jacobian...
      &             A11,A12,A13,&
      &             A21,A22,A23,&
      &             A31,A32,A33,&
      ! THE LIMITING VARIABLE IN[0,1]
      &             UUD,UVD,UWD,LOVERQ )
    REAL TOLER
    PARAMETER(TOLER=1.E-7)
    LOGICAL D3
    INTEGER NGI,NLOC,NONODS,TOTELE
    INTEGER ELE,GI 
    INTEGER NDGLNO(NLOC*TOTELE) 
    REAL UD(NGI),VD(NGI),WD(NGI)  
    REAL U(NONODS)
    REAL HOVERQ 
    REAL N(NLOC,NGI),NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    LOGICAL RESTRI
    ! CENTRE is true when we are to apply diffusion in cross stream 
    ! directions. 
    REAL DETJ,&
        &             AGI,BGI,&
        &             CGI,DGI,&
        ! 3-d inverse of Jacobian...
        &             A11,A12,A13,&
        &             A21,A22,A23,&
        &             A31,A32,A33
    REAL UUD,UVD,UWD,LOVERQ 
    ! LIM1,LIM2,LIM3 are the limiters for x,y & z directions respectively.
    ! LIMITU 
    ! LOCAL VARIABLES...
    !
    REAL UDERI,RN
    !
    REAL HX,HY,HZ
    INTEGER ILOC,NODX,NOD
  
    ! This subroutine limits the gradiant.
    ! Determine UUD, UVD at the Gauss pt GI
  
    IF(D3) THEN
      UUD=0.
      UVD=0.
      UWD=0.
      UDERI=0.
      do ILOC=1,NLOC
          NOD=NDGLNO((ELE-1)*NLOC+ILOC)
          UDERI=UDERI+(UD(GI)*NX(ILOC,GI)+VD(GI)*NY(ILOC,GI)&
              &           +WD(GI)*NZ(ILOC,GI))*U(NOD)
          UUD=UUD+NX(ILOC,GI)*U(NOD)
          UVD=UVD+NY(ILOC,GI)*U(NOD)
          UWD=UWD+NZ(ILOC,GI)*U(NOD)
      END DO
      RN=MAX(UUD**2+UVD**2+UWD**2,TOLER)
      UUD=UDERI*UUD/RN
      UVD=UDERI*UVD/RN
      UWD=UDERI*UWD/RN
      !
      HX=ABS(A11*UUD+A12*UVD+A13*UWD)
      HY=ABS(A21*UUD+A22*UVD+A23*UWD)
      HZ=ABS(A31*UUD+A32*UVD+A33*UWD)
      LOVERQ=2./MAX(HX,HY,HZ,TOLER)
    ELSE
      UUD=0.
      UVD=0.
      UDERI=0.
      do ILOC=1,NLOC
          NOD=NDGLNO((ELE-1)*NLOC+ILOC)
          UDERI=UDERI+(UD(GI)*NX(ILOC,GI)+VD(GI)*NY(ILOC,GI))&
              &           *U(NOD)
          UUD=UUD+NX(ILOC,GI)*U(NOD)
          UVD=UVD+NY(ILOC,GI)*U(NOD)
      END DO
      RN=MAX(UUD**2+UVD**2,TOLER)
      UUD=UDERI*UUD/RN
      UVD=UDERI*UVD/RN
      !
      HX=ABS( DGI*UUD-BGI*UVD)/DETJ
      HY=ABS(-CGI*UUD+AGI*UVD)/DETJ
      LOVERQ=2./MAX(HX,HY,TOLER)
    ENDIF
    ! Amend to take into account streamline upwinding
    IF(RESTRI) LOVERQ=MAX(0.0,LOVERQ-HOVERQ)
    !
    RETURN
  END SUBROUTINE LIMGRA

end module limit_gradient

!
!
