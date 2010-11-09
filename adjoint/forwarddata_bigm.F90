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

SUBROUTINE FORWARDDATA_BIGM(BIGM1,NBIGM,D3,ITSTIME,DIRNAME)
   use FLDebug
  IMPLICIT NONE
  LOGICAL ::D3
  INTEGER ::NBIGM
  REAL    ::BIGM1(NBIGM)
  INTEGER ::NOVA
  INTEGER ::ITSTIME
  CHARACTER(40) DIRNAME
  CHARACTER(40) FILENAME_BIGM,FILENAME
  CHARACTER(40) PWD
  INTEGER I,K1,K2
  
  WRITE(FILENAME_BIGM, 10) ITSTIME
10 FORMAT('/bigm_data',I5.5)
   ewrite(3,*) 'DIRNAD = ', DIRNAME
   ewrite(3,*) 'FILENAME_BIGM = ', FILENAME_BIGM
         K1 = index( DIRNAME, ' ' ) - 1
        K2 = index( FILENAME_BIGM, ' ' ) - 1
        FILENAME=DIRNAME(1:K1) // FILENAME_BIGM(1:K2)
        ewrite(3,*)FILENAME
         OPEN(1,FILE = FILENAME)
 !K1 = LEN_TRIM( DIRNAME )
 ! K2 = LEN_TRIM( FILENAME_BIGM )

 ! FILENAME=DIRNAME(1:K1) // FILENAME_BIGM(1:K2)
  ! print*,'FILENAME = ',FILENAME
  OPEN(1,FILE=FILENAME(1:LEN_TRIM(FILENAME)))
  
  NOVA = SIZE(BIGM1)
  WRITE(1,*) (BIGM1(I),I=1,NOVA)
  
  CLOSE(1)
  
END SUBROUTINE FORWARDDATA_BIGM
