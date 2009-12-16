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

module position_in_matrix

  implicit none

  private

  public :: posinmat

contains

  SUBROUTINE POSINMAT(POSMAT,GLOBI,GLOBJ,&
       &               NONODS,FINDRM,COLM,NCOLM)
    IMPLICIT NONE
    INTEGER NONODS,NCOLM
    INTEGER POSMAT,GLOBI,GLOBJ,FINDRM(NONODS+1),COLM(NCOLM)
    ! Local variables...
    INTEGER INUM,LOWER,UPPER
    ! Find position in matrix POSMAT which has column GLOBJ
    ! Find POSMAT. ***************************************
    LOWER=FINDRM(GLOBI) 
    UPPER=FINDRM(GLOBI+1)-1
7000 CONTINUE
    INUM=LOWER+(UPPER-LOWER+1)/2 
    IF(GLOBJ.GE.COLM(INUM) )  THEN 
       LOWER=INUM
    ELSE
       UPPER=INUM
    ENDIF
    IF(UPPER-LOWER.LE.1) THEN
       IF(GLOBJ.EQ.COLM(LOWER)) THEN
          POSMAT=LOWER
       ELSE
          POSMAT=UPPER
       ENDIF
       GOTO 9000
    ENDIF
    GOTO 7000
9000 CONTINUE
    ! Found POSMAT ***************************************
    RETURN 
  END SUBROUTINE POSINMAT

end module position_in_matrix

