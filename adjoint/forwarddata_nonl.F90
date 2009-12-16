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
SUBROUTINE FORWARDDATA_NONL(NU,NV,NW,NONODS,NPHASE,D3,ITSTIME,DIRNAME)
   use FLDebug
  IMPLICIT NONE
  INTEGER ::NONODS,NPHASE
  REAL    ::NU(NONODS*NPHASE),NV(NONODS*NPHASE),NW(NONODS*NPHASE)
  LOGICAL ::D3
  INTEGER ::ITSTIME
  CHARACTER(40) DIRNAME
  CHARACTER(40) FILENAME_NONL,FILENAME
  INTEGER I,K1,K2

  WRITE(FILENAME_NONL, 10) ITSTIME
10 FORMAT('/nonlinear_data',I5.5)

  ewrite(3,*) DIRNAME
  ewrite(3,*) FILENAME_NONL
  K1 = index( DIRNAME, ' ' ) - 1
  K2 = index( FILENAME_NONL, ' ' ) - 1
  FILENAME=DIRNAME(1:K1) // FILENAME_NONL(1:K2)
  ewrite(3,*) FILENAME
  OPEN(1,FILE = FILENAME)
  ewrite(3,*) NONODS
  ewrite(3,*) NPHASE
  ewrite(3,*) NU(1)
  ewrite(3,*) NU(NONODS*NPHASE)
  WRITE(1,*) (NU(I),I=1,NONODS*NPHASE)

  WRITE(1,*) (NV(I),I=1,NONODS*NPHASE)

  IF(D3) THEN
     WRITE(1,*) (NW(I),I=1,NONODS*NPHASE)
  ENDIF

  CLOSE(1)

END SUBROUTINE FORWARDDATA_NONL
