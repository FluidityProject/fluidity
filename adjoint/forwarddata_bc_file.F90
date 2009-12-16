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
	SUBROUTINE FORWARDDATA_BC_FILE(BCU1,BCU2,BCV1,BCV2,BCW1,BCW2,   &
	           NOBCU,NOBCV,NOBCW,D3,ITSTIME,DIRNAME)
	   use FLDebug
        IMPLICIT NONE
	INTEGER ::NOBCU,NOBCV,NOBCW
	INTEGER ::BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
	REAL    ::BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
	LOGICAL ::D3
	INTEGER ::NOVA
	INTEGER ::ITSTIME
	CHARACTER(40) ::DIRNAME

	CHARACTER(40) ::FILENAME,FILENAME_BC
	INTEGER ::I,K1,K2

	WRITE(FILENAME_BC, 10) ITSTIME
10 	FORMAT('/bc_data',I5.5)
	ewrite(3,*) DIRNAME
	ewrite(3,*) FILENAME_BC
	K1 = index( DIRNAME, ' ' ) - 1
	K2 = index( FILENAME_BC, ' ' ) - 1
	ewrite(3,*) 'K2=',K2
	FILENAME=DIRNAME(1:K1) // FILENAME_BC(1:K2)
	ewrite(3,*) FILENAME
 	OPEN(1,FILE = FILENAME)
	NOVA = SIZE(BCU1)
	ewrite(3,*) 'NOVA=',NOVA
	ewrite(3,*) 'NOBCU',NOBCU
	WRITE(1,*) (BCU1(I),I=1,NOVA)
	NOVA = SIZE(BCU2)
	WRITE(1,*) (BCU2(I),I=1,NOVA)

	NOVA = SIZE(BCV1)
	ewrite(3,*) 'NOVA=',NOVA
	ewrite(3,*) 'NOBCV',NOBCV
	WRITE(1,*) (BCV1(I),I=1,NOVA)
	NOVA = SIZE(BCV2)
	WRITE(1,*) (BCV2(I),I=1,NOVA)

	IF(D3) THEN
	NOVA = SIZE(BCW1)
	ewrite(3,*) 'NOVA=',NOVA
	ewrite(3,*) 'NOBCW',NOBCW
	WRITE(1,*) (BCW1(I),I=1,NOVA)
	NOVA = SIZE(BCW2)
	WRITE(1,*) (BCW2(I),I=1,NOVA)
	ENDIF

	CLOSE(1)

        END SUBROUTINE FORWARDDATA_BC
