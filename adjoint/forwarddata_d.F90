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
SUBROUTINE FORWARDDATA_D(DU,BCU2,DV,BCV2,DW,BCW2,   &
     NOBCU,NOBCV,NOBCW,D3,ITSTIME,MAXREC,MRECORD,          &
     NOCVA,D,GOTFRE)
   use FLDebug
  IMPLICIT NONE
  INTEGER ::NOBCU,NOBCV,NOBCW
  INTEGER ::BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
  REAL    ::DU(NOBCU),DV(NOBCV),DW(NOBCW)
  LOGICAL ::D3
  INTEGER ::NOCVA
  REAL	  ::D(NOCVA)
  INTEGER ::ITSTIME,MAXREC
  INTEGER ::MRECORD(MAXREC)

  LOGICAL ::GOTFRE

  INTEGER ::K1,K2,II

  ewrite(3,*) 'NOBCU,NOBCV,NOBCW',NOBCU,NOBCV,NOBCW
  IF(GOTFRE) THEN
    K1 = MRECORD(ITSTIME)
    K2 = MRECORD(ITSTIME)+NOBCU-1
    D(K1:K2) = DU

  ELSE
    K1 = MRECORD(ITSTIME)
    K2 = MRECORD(ITSTIME)+NOBCU-1
    D(K1:K2) = DU

    K1 = K2+1
    K2 = K2+NOBCV
    D(K1:K2) = DV

    IF(D3) THEN
      K1 = K2+1
      K2 = K2+NOBCW
      D(K1:K2) = DW
    ENDIF

  ENDIF

END SUBROUTINE FORWARDDATA_D
