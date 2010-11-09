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
SUBROUTINE FORWARDDATA_CONTROL_VARIABLES(CVA,IMEM_CVA,D,             &
        CVA_X,CVA_Y,CVA_Z,NOCVA,D3)
   use FLDebug
  IMPLICIT NONE
  INTEGER ::NOCVA
  REAL    ::CVA(NOCVA),CVA_X(NOCVA),CVA_Y(NOCVA),CVA_Z(NOCVA)
  REAL    ::D(NOCVA)  
  INTEGER ::IMEM_CVA(NOCVA)
  LOGICAL ::D3
! local......
  INTEGER ::I

  OPEN(1,FILE = 'control_variables.dat',STATUS='REPLACE')
  WRITE(1,*) NOCVA
  WRITE(1,*) (IMEM_CVA(I),I=1,NOCVA)
  WRITE(1,*) (CVA(I),I=1,NOCVA)
!  WRITE(1,*) (CVA_X(I),I=1,NOCVA)
!  WRITE(1,*) (CVA_Y(I),I=1,NOCVA)
!  IF(D3) THEN
!     WRITE(1,*) (CVA_Z(I),I=1,NOCVA)
!  ENDIF
  WRITE(1,*) (D(I),I=1,NOCVA)
  CLOSE(1)

END SUBROUTINE FORWARDDATA_CONTROL_VARIABLES
