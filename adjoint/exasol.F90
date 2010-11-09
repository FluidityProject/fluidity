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

  SUBROUTINE EXASOL(NOSNOD,NOSTIM,DT,UEXAC,VEXAC,WEXAC,SX,SY,SZ,STIME,D3,FILNAM_SOURCE,  &
       FSZERO,FSSCALE)
    
    !----------------------------------------------------
    ! This subroutine is used to get the exact solutions
    !----------------------------------------------------
    
    use FLDebug
    IMPLICIT NONE
    INTEGER ::NOSTIM,NOSNOD,ITSTIME
    REAL    ::UEXAC(NOSNOD,NOSTIM),VEXAC(NOSNOD,NOSTIM),WEXAC(NOSNOD,NOSTIM)
    REAL    ::SX(NOSNOD),SY(NOSNOD),SZ(NOSNOD)
    REAL    ::STIME(NOSTIM),ST,DT,FSZERO,FSSCALE
    LOGICAL ::D3
    CHARACTER(40) ::UNAME,VNAME,WNAME,SXNAME,SYNAME,SZNAME,STIMENAME
    CHARACTER*240 ::FILNAM_SOURCE
    ! Local variables.......
    INTEGER ::II,KK,I
    REAL    ::AA
    
    
    STIME(1:NOSTIM) = 0.0
    
    DO KK =2,NOSTIM
       STIME(KK) = STIME(KK-1)+DT
    ENDDO
    UEXAC=0.0
    VEXAC=0.0
    WEXAC=0.0
    
    OPEN(1,FILE=FILNAM_SOURCE,STATUS='OLD')
    
    
    
    ! Read the obervational data on free surface.... ignor the height change for the time being
    
    DO II=1,NOSNOD
       READ(1,*) UNAME
       ewrite(3,*) UNAME
       READ(1,*) SX(II),SY(II)
       DO KK =1,NOSTIM
          READ(1,*) ST,UEXAC(II,KK)
          UEXAC(II,KK)=UEXAC(II,KK)*FSSCALE+FSZERO
       ENDDO
    ENDDO
    
    CLOSE(1)

  END SUBROUTINE EXASOL
  
  
  SUBROUTINE EXASOL_UVW_READ(NOSNOD,NOSTIM,DT,UEXAC,VEXAC,WEXAC,SX,SY,SZ,STIME,D3,FILNAM_SOURCE)
    
    !----------------------------------------------------
    ! This subroutine is used to get the exact solutions
    !----------------------------------------------------
    
    use FLDebug
    IMPLICIT NONE
    INTEGER ::NOSTIM,NOSNOD,ITSTIME
    REAL    ::UEXAC(NOSNOD,NOSTIM),VEXAC(NOSNOD,NOSTIM),WEXAC(NOSNOD,NOSTIM)
    REAL    ::SX(NOSNOD),SY(NOSNOD),SZ(NOSNOD)
    REAL    ::STIME(NOSTIM),ST,DT,FSZERO,FSSCALE
    LOGICAL ::D3
    CHARACTER(40) ::UNAME,VNAME,WNAME,SXNAME,SYNAME,SZNAME,STIMENAME
    CHARACTER*240 ::FILNAM_SOURCE
    ! Local variables.......
    INTEGER ::II,KK,I
    REAL    ::AA,ACCTIM
    
    
    UEXAC=0.0
    VEXAC=0.0
    WEXAC=0.0
    SX=0.0
    SY=0.0
    SZ=0.0
    STIME=0.0
    
    DO KK =2,NOSTIM
       STIME(KK) = STIME(KK-1)+DT
    ENDDO
    
    OPEN(1,FILE='AdjSourceData.dat',STATUS='OLD')
    ! Read the obervational data on free surface.... ignor the height change for the time being
    
    DO KK =1,NOSTIM
       IF( KK.EQ.1 ) THEN
          READ(1,*) UNAME
          READ(1,*) (SX(II),II=1,NOSNOD)
          READ(1,*) UNAME
          READ(1,*) (SY(II),II=1,NOSNOD)
          IF(D3) THEN
             READ(1,*) UNAME
             READ(1,*) (SZ(II),II=1,NOSNOD)
          ENDIF
          
          READ(1,*) UNAME
          READ(1,*) (STIME(II),II=1,NOSTIM)
       ENDIF
       
       READ(1,*)  ACCTIM
       READ(1,*) (UEXAC(II,KK),II=1,NOSNOD)
       READ(1,*) (VEXAC(II,KK),II=1,NOSNOD)
       IF(D3) THEN
          READ(1,*) (WEXAC(II,KK),II=1,NOSNOD)
       ENDIF
       
    ENDDO
    
    CLOSE(1)
    
  END SUBROUTINE EXASOL_UVW_READ
  
  
  
 
 SUBROUTINE EXASOL_SEN(NOSNOD,NOSTIM,DT,UEXAC,VEXAC,WEXAC,SX,SY,SZ,STIME,D3,FILNAM_SOURCE,  &
      FSZERO,FSSCALE)
   
   !----------------------------------------------------
   ! This subroutine is used to get the exact solutions
   !----------------------------------------------------
   
   use FLDebug
   IMPLICIT NONE
   INTEGER ::NOSTIM,NOSNOD,ITSTIME
   REAL    ::UEXAC(NOSNOD,NOSTIM),VEXAC(NOSNOD,NOSTIM),WEXAC(NOSNOD,NOSTIM)
   REAL    ::SX(NOSNOD),SY(NOSNOD),SZ(NOSNOD)
   REAL    ::STIME(NOSTIM),ST,DT,FSZERO,FSSCALE
   LOGICAL ::D3
   CHARACTER(40) ::UNAME,VNAME,WNAME,SXNAME,SYNAME,SZNAME,STIMENAME
   CHARACTER*240 ::FILNAM_SOURCE
   ! Local variables.......
   INTEGER ::II,KK,I
   REAL    ::AA
   
   
   STIME(1:NOSTIM) = 0.0
   
   DO KK =2,NOSTIM
      STIME(KK) = STIME(KK-1)+DT
   ENDDO
   UEXAC=0.0
   VEXAC=0.0
   WEXAC=0.0
   
   OPEN(1,FILE=FILNAM_SOURCE,STATUS='OLD')
   
   
   
   ! Read the obervational data on free surface.... ignor the height change for the time being
   
   DO II=1,NOSNOD
      READ(1,*) UNAME
      ewrite(3,*) UNAME
      READ(1,*) SX(II),SY(II)
      !    DO KK =1,NOSTIM
      !      READ(1,*) ST,UEXAC(II,KK)
      !      UEXAC(II,KK)=UEXAC(II,KK)*FSSCALE+FSZERO
      !    ENDDO
   ENDDO
   
   CLOSE(1)
   
 END SUBROUTINE EXASOL_SEN
 
 
