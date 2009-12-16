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
SUBROUTINE FORWARDDATA_ADJUVW(U,V,W,NONODS,NPHASE,D3,ITSTIME,DIRNAME)
   use FLDebug
  IMPLICIT NONE
  INTEGER ::NONODS,NPHASE
  REAL    ::U(NONODS*NPHASE),V(NONODS*NPHASE),W(NONODS*NPHASE)
  LOGICAL ::D3
  INTEGER ::ITSTIME
  CHARACTER(40) DIRNAME
  CHARACTER(40) FILENAME_ADJUVW,FILENAME
  INTEGER I,K1,K2

  WRITE(FILENAME_ADJUVW, 10) ITSTIME
10 FORMAT('/adjuvw_data',I5.5)

  ewrite(3,*) DIRNAME
  ewrite(3,*) FILENAME_ADJUVW
  K1 = index( DIRNAME, ' ' ) - 1
  K2 = index( FILENAME_ADJUVW, ' ' ) - 1
  FILENAME=DIRNAME(1:K1) // FILENAME_ADJUVW(1:K2)
  ewrite(3,*) FILENAME
  OPEN(1,FILE = FILENAME)

  WRITE(1,*) (U(I),I=1,NONODS*NPHASE)

  WRITE(1,*) (V(I),I=1,NONODS*NPHASE)

  IF(D3) THEN
     WRITE(1,*) (W(I),I=1,NONODS*NPHASE)
  ENDIF

  CLOSE(1)

END SUBROUTINE FORWARDDATA_ADJUVW


!===================================================

SUBROUTINE FORWARDDATA_ADJH(HLSTDT,NONODS,NPHASE,D3,ITSTIME,DIRNAME)
   use FLDebug
  IMPLICIT NONE
  INTEGER ::NONODS,NPHASE
  REAL    ::HLSTDT(NONODS*NPHASE)
  LOGICAL ::D3
  INTEGER ::ITSTIME
  CHARACTER(40) DIRNAME
  CHARACTER(40) FILENAME_ADJH,FILENAME
  INTEGER I,K1,K2

  WRITE(FILENAME_ADJH, 10) ITSTIME
10 FORMAT('/adjh_data',I5.5)

  ewrite(3,*) DIRNAME
  ewrite(3,*) FILENAME_ADJH
  K1 = index( DIRNAME, ' ' ) - 1
  K2 = index( FILENAME_ADJH, ' ' ) - 1
  FILENAME=DIRNAME(1:K1) // FILENAME_ADJH(1:K2)
  ewrite(3,*) FILENAME
  OPEN(1,FILE = FILENAME)

  WRITE(1,*) (HLSTDT(I),I=1,NONODS*NPHASE)

  CLOSE(1)

END SUBROUTINE FORWARDDATA_ADJH

SUBROUTINE FORWARDDATA_ADJH1(NOCVA,NONODS,NPHASE,G,IMEM_CVA,MAXREC,MRECORD,ITSTIME,NTIME,NOBCH_F,DIRNAME)
   use FLDebug
  IMPLICIT NONE
  INTEGER ::NOCVA,NONODS,NPHASE,MAXREC
  REAL    ::G(NOCVA)
  INTEGER  MRECORD(MAXREC),IMEM_CVA(NOCVA)
  INTEGER ::NOBCH_F(NTIME)
  LOGICAL ::D3
  INTEGER ::ITSTIME,NTIME
  CHARACTER(40) DIRNAME
  CHARACTER(40) FILENAME_ADJH,FILENAME
  INTEGER I,K1,K2,NOD
  REAL,DIMENSION(:),ALLOCATABLE ::GLOCAL

  ALLOCATE(GLOCAL(NONODS*NPHASE) )
  GLOCAL(1:NONODS*NPHASE) = 0.0

  DO I= MRECORD(ITSTIME),MRECORD(ITSTIME)+NOBCH_F(NTIME-ITSTIME+1)-1
    NOD=IMEM_CVA(I)
    GLOCAL(NOD) = G(I)
  ENDDO

  WRITE(FILENAME_ADJH, 10) NTIME-ITSTIME+1
10 FORMAT('/adjh_data',I5.5)

  ewrite(3,*) DIRNAME
  ewrite(3,*) FILENAME_ADJH
  K1 = index( DIRNAME, ' ' ) - 1
  K2 = index( FILENAME_ADJH, ' ' ) - 1
  FILENAME=DIRNAME(1:K1) // FILENAME_ADJH(1:K2)
  ewrite(3,*) FILENAME

  OPEN(1,FILE = FILENAME)
  WRITE(1,*) (GLOCAL(I),I=1,NONODS*NPHASE)
  CLOSE(1)

  DEALLOCATE(GLOCAL)

END SUBROUTINE FORWARDDATA_ADJH1
