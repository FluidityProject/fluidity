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
#include "confdefs.h"
#include "fmangle.h"

INTEGER FUNCTION NID(i,j,STRIDE)
  INTEGER, INTENT(IN)::i,j,STRIDE
  NID = (i-1)*STRIDE + j
END FUNCTION NID

PROGRAM INTERPO

  IMPLICIT NONE
  EXTERNAL NID
  INTEGER NID
  INTEGER,PARAMETER::NOX1=1000,NOY1=500
  INTEGER,PARAMETER::NONODS1=NOX1*NOY1
  REAL X1(NONODS1),Y1(NONODS1),H1(NONODS1)
  INTEGER,PARAMETER::TOTELE1=(NOY1-1)*(NOX1-1)*2,NLOC1=3
  INTEGER NDGLNO1(TOTELE1*NLOC1)
  REAL,PARAMETER::LX1=1.0/(NOX1-1),LY1=1.0/(NOY1-1)

  INTEGER,PARAMETER::NOX2=333,NOY2=333
  INTEGER,PARAMETER::NONODS2=NOX2*NOY2
  REAL X2(NONODS2),Y2(NONODS2),H2(NONODS2),INTERPOL(NONODS2)
  INTEGER,PARAMETER::TOTELE2=(NOY2-1)*(NOX2-1)*2,NLOC2=3
  INTEGER NDGLNO2(TOTELE2*NLOC2)
  REAL,PARAMETER::LX2=1.0/(NOX2-1),LY2=1.0/(NOY2-1) 

  ! LOCAL.......
  INTEGER ::NOD,II,JJ,ELE,ILOC,NOD1,I1,I2,IERROR
  INTEGER,PARAMETER::NFIELDS=1
  INTEGER FIELDS1(NFIELDS),FIELDS2(NFIELDS)
  REAL L2NORM

  FIELDS1(1) = 1
  FIELDS2(1) = 1

  ! Mesh 1 and variable H1
  !-------------

  DO NOD = 1,NONODS1
     X1(NOD) = 0.0
     Y1(NOD) = 0.0
     H1(NOD) = 0.0
  ENDDO

  DO II =1,TOTELE1*NLOC1
     NDGLNO1(II) = -1
  ENDDO

  DO II = 1,NOX1
     DO JJ = 1,NOY1
        NOD = NID(II,JJ,NOY1)
        X1(NOD) = LX1*(II-1)
        Y1(NOD) = LY1*(JJ-1)
        H1(NOD) = SIN( X1(NOD) ) + SIN( Y1(NOD) )
     ENDDO
  ENDDO
  
  ELE=0
  DO II=1,NOX1-1
     DO JJ=1,NOY1-1
        NDGLNO1(ELE*NLOC1+1) = NID(II,  JJ,  NOY1)
        NDGLNO1(ELE*NLOC1+2) = NID(II+1,JJ,  NOY1)
        NDGLNO1(ELE*NLOC1+3) = NID(II,  JJ+1,NOY1)
        ELE = ELE+1

        NDGLNO1(ELE*NLOC1+1) = NID(II+1,JJ,  NOY1)
        NDGLNO1(ELE*NLOC1+2) = NID(II+1,JJ+1,NOY1)
        NDGLNO1(ELE*NLOC1+3) = NID(II,  JJ+1,NOY1)
        ELE = ELE+1
     ENDDO
  ENDDO
  
  ! Mesh 2 and variable H2
  !-------------
  DO NOD = 1,NONODS2
     X2(NOD) = 0.0
     Y2(NOD) = 0.0
     H2(NOD) = 0.0
  ENDDO

  DO II =1,TOTELE2*NLOC2
     NDGLNO2(II)=-1
  ENDDO

  DO II = 1,NOX2
     DO JJ = 1,NOY2
        NOD = NID(II,JJ,NOY2)
        X2(NOD) = LX2*(II-1)
        Y2(NOD) = LY2*(JJ-1)
        H2(NOD) = SIN( X2(NOD) ) + SIN( Y2(NOD) )
        INTERPOL(NOD) = 0
     ENDDO
  ENDDO

  ELE=0
  DO II=1,NOX2-1
     DO JJ=1,NOY2-1
        NDGLNO2(ELE*NLOC2+1) = NID(II,  JJ,  NOY2)
        NDGLNO2(ELE*NLOC2+2) = NID(II+1,JJ,  NOY2)
        NDGLNO2(ELE*NLOC2+3) = NID(II,  JJ+1,NOY2)
        ELE = ELE+1

        NDGLNO2(ELE*NLOC2+1) = NID(II+1,JJ,  NOY2)
        NDGLNO2(ELE*NLOC2+2) = NID(II+1,JJ+1,NOY2)
        NDGLNO2(ELE*NLOC2+3) = NID(II,  JJ+1,NOY2)
        ELE = ELE+1
     ENDDO
  ENDDO

  ewrite(3,*)  "calling fltri3totri3()"
  CALL fltri3totri3(NONODS1, TOTELE1, X1, Y1, &
       NDGLNO1, H1, FIELDS1, NFIELDS,         &
       NONODS2, TOTELE2, X2, Y2,              &
       NDGLNO2, INTERPOL, FIELDS2, IERROR)

  L2NORM = 0
  DO II=1,NONODS2
    L2NORM = L2NORM + (INTERPOL(II)-H2(II))*(INTERPOL(II)-H2(II))
  END DO
  ewrite(3,*)  "L2NORM = ", L2NORM
END PROGRAM INTERPO
