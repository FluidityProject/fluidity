
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
!
#include "fdebug.h"
#include "fmangle.h"

SUBROUTINE CREATDIR(DIRNAME)
  use FLDebug
  IMPLICIT NONE
  CHARACTER(*) DIRNAME
  INTEGER ILEN, IERR
  
  IERR=0
  
  WRITE(DIRNAME, 10)
10 FORMAT('adresults')
  
  ILEN = LEN_TRIM(DIRNAME)

  CALL flmkdir(DIRNAME,ILEN,IERR)
  
  IF(IERR .LT. 0) THEN
     ewrite(3,*)  &
          'WARNING: The dir already exists. Are you sure you want this?'
  END IF
  
END SUBROUTINE CREATDIR

