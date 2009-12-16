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

module assemble_boundary_conditions
  use fldebug
  implicit none

  private
  public :: boucon, bouvis

contains

  SUBROUTINE BOUCON(IDIM,NONODS,DM1,NOBCU,BCU1,BCU2,FX,DM2,NOBCV,BCV1,BCV2&
       &,FY,DM3,NOBCW,BCW1,BCW2,FZ)

    !     FORMULATE BOUNDARY CONDITIONS
    !     This subroutine adds in the DIRICHLET boundary conditions. 
    !     IDIM controls how many values to set. 
    !     BCU2,BCV2,BCW2 contain the nodes at which the boundary 
    !     conditions are to be applied. 
    !     BCU1 contains  the value of the specified acceleration of U 
    !     that is BCU1=(Unew - Uold)/DT
    REAL INFINY,INFIN1,INFIN2,NEAR
#ifndef DOUBLEP
    !     single prec...
    PARAMETER(INFIN1=10.0E+10)
    PARAMETER(INFIN2=10.0E+10)
    !     double precision...
#else 
    PARAMETER(INFIN1=10.0E+25)
    PARAMETER(INFIN2=10.0E+25)
#endif
    PARAMETER(NEAR=1./INFIN1)

    INTEGER IDIM,NONODS,NOBCU,NOBCV,NOBCW
    INTEGER BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
    REAL BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
    REAL DM1(NONODS),DM2(NONODS),DM3(NONODS)
    REAL FX(NONODS),FY(NONODS),FZ(NONODS)

    INTEGER ROW,I

    DO I=1,NONODS
       DM1(I)=0.
    END DO
    DO I=1,NOBCU
       ROW=BCU2(I)
       INFINY=INFIN2
       IF(ABS(BCU1(I)).GT.NEAR) THEN 
          INFINY=INFIN1
       ENDIF
       DM1(ROW)=INFINY
       FX(ROW)=INFINY* BCU1(I)
    END DO

    IF(IDIM.GE.2) THEN
       DO I=1,NONODS
          DM2(I)=0.
       END DO
       DO I=1,NOBCV
          ROW=BCV2(I)
          INFINY=INFIN2
          IF(ABS(BCV1(I)).GT.NEAR) THEN
             INFINY=INFIN1
          ENDIF
          DM2(ROW)=INFINY
          FY(ROW)=INFINY* BCV1(I)
       END DO
    ENDIF

    IF(IDIM.GE.3) THEN
       DO I=1,NONODS
          DM3(I)=0.
       END DO
       DO I=1,NOBCW
          ROW=BCW2(I)
          INFINY=INFIN2
          IF(ABS(BCW1(I)).GT.NEAR) THEN
             INFINY=INFIN1
          ENDIF
          DM3(ROW)=INFINY
          FZ(ROW)=INFINY* BCW1(I)
       END DO
    ENDIF

  END SUBROUTINE BOUCON

  SUBROUTINE BOUVIS(IDIM, NONODS, DM1,DM2,DM3,&
       &     BIGM1,NBIGM,NCOLM,CENTRM)
    INTEGER IDIM,NONODS,NBIGM,NCOLM,CENTRM(NONODS)
    REAL BIGM1(NBIGM), DM1(NONODS),DM2(NONODS),DM3(NONODS)
    INTEGER I,ICENT

    IF(IDIM.EQ.3) THEN
       IF(NBIGM.EQ.6*NCOLM) THEN
          do  I=1,NONODS! Was loop 146
             ICENT=CENTRM(I)
             BIGM1(ICENT)        =BIGM1(ICENT)        +DM1(I)
             BIGM1(ICENT+2*NCOLM)=BIGM1(ICENT+2*NCOLM)+DM2(I)
             BIGM1(ICENT+5*NCOLM)=BIGM1(ICENT+5*NCOLM)+DM3(I)
          end do ! Was loop 146
       ELSE
          do  I=1,NONODS! Was loop 2146
             ICENT=CENTRM(I)
             BIGM1(ICENT)        =BIGM1(ICENT)        +DM1(I)
             BIGM1(ICENT+4*NCOLM)=BIGM1(ICENT+4*NCOLM)+DM2(I)
             BIGM1(ICENT+8*NCOLM)=BIGM1(ICENT+8*NCOLM)+DM3(I)
          end do ! Was loop 2146
       ENDIF
    ENDIF
    IF(IDIM.EQ.2) THEN
       IF(NBIGM.EQ.3*NCOLM) THEN
          do  I=1,NONODS! Was loop 246
             ICENT=CENTRM(I)
             BIGM1(ICENT)        =BIGM1(ICENT)        +DM1(I)
             BIGM1(ICENT+2*NCOLM)=BIGM1(ICENT+2*NCOLM)+DM2(I)
          end do ! Was loop 246
       ELSE
          do  I=1,NONODS! Was loop 346
             ICENT=CENTRM(I)
             BIGM1(ICENT)        =BIGM1(ICENT)        +DM1(I)
             BIGM1(ICENT+3*NCOLM)=BIGM1(ICENT+3*NCOLM)+DM2(I)
          end do ! Was loop 346
       ENDIF
    ENDIF
  END SUBROUTINE BOUVIS

end module assemble_boundary_conditions
