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
  SUBROUTINE GRADIENT_BC_SIMPLE(UADJ,VADJ,WADJ,                        &
     NONODS,NOBCU_F,NOBCV_F,NOBCW_F,MAXGRAD,D3,                              &
     ITSTIME,MAXREC,MRECORD,NOCVA,G,IMEM_CVA,NTIME)          
   use FLDebug
  IMPLICIT NONE  
  !*************************************************************************
  ! This subroutine is used to calculate the gradient of the Cost function
  !*************************************************************************
  ! Find gradient and functional...
  REAL         ::MAXGRAD
  LOGICAL   ::D3
  INTEGER   ::NONODS,NTIME
  INTEGER   ::NOBCU_F(NTIME),NOBCV_F(NTIME),NOBCW_F(NTIME)
  REAL               ::UADJ(NONODS),VADJ(NONODS),WADJ(NONODS)
  INTEGER   ::NOCVA       ! number of the control variables
  REAL               ::G(NOCVA)    ! the vector of the gradient
  INTEGER   ::IMEM_CVA(NOCVA)
  INTEGER   ::ITSTIME,MAXREC
  INTEGER             ::MRECORD(MAXREC)
  INTEGER             ::LSGRAD
  ! Local....
  INTEGER ::K1,K2                 ! the records for writing to G
  INTEGER ::II,NOD

  IF(ITSTIME.EQ.1) THEN
     MAXGRAD = 0.0
  ENDIF

  ! Local total size of G at each time step
  LSGRAD = NOBCU_F(NTIME-ITSTIME+1)+NOBCV_F(NTIME-ITSTIME+1)
  IF(D3) THEN
     LSGRAD = LSGRAD + NOBCW_F(NTIME-ITSTIME+1)
  ENDIF

  ewrite(3,*) 'NOCVA in SUB_grad',NOCVA
  K1 = MRECORD(ITSTIME)
  K2 = K1+NOBCU_F(NTIME-ITSTIME+1)-1
  ewrite(3,*) 'K1,K2',K1,K2
  DO II=1,NOBCU_F(NTIME-ITSTIME+1)
     NOD=IMEM_CVA(K1+II-1)
!     ewrite(3,*) 'NOD,UADJ(NOD)',NOD,UADJ(NOD)
     if(nod.eq.0) then
     ewrite(3,*) 'NOD,ii',NOD,ii
     stop 3333
     endif
     G(K1+II-1) = UADJ(NOD)
    if(itstime.eq.ntime) then
       G(K1+II-1) = 0.0  ! not adjust the initial condition
    endif
     MAXGRAD = MAX(MAXGRAD,ABS(G(K1+II-1)) )
  ENDDO

  IF(ITSTIME.EQ.1) THEN
     OPEN(4,FILE='gradient_time.dat',STATUS='REPLACE')
  ELSE
     OPEN(4,FILE='gradient_time.dat',POSITION='APPEND')
  ENDIF
  write(4,*) 'ITSTIME=',NTIME-ITSTIME+1
  write(4,*) (G(K1+II-1),II=1,NOBCU_F(NTIME-ITSTIME+1))
  CLOSE(4)

  K1 = K2+1             
  K2 = K1+NOBCV_F(NTIME-ITSTIME+1)-1
  ewrite(3,*) 'K1,K2',K1,K2
  DO II=1,NOBCV_F(NTIME-ITSTIME+1)
     NOD=IMEM_CVA(K1+II-1)
!     if(nod.eq.0) then
!     ewrite(3,*) 'NOD,ii',NOD,ii
!     stop 3333
!     endif
     ewrite(3,*) 'II,NOCVA,nod,nonods',II,NOCVA,nod,nonods
     G(K1+II-1) = VADJ(NOD)
     MAXGRAD = MAX(MAXGRAD,ABS(G(K1+II-1)) )
  ENDDO

  IF(D3) THEN
     K1 = K2+1           
     K2 = K1+NOBCW_F(NTIME-ITSTIME+1)-1
     ewrite(3,*) 'K1,K2',K1,K2
     DO II=1,NOBCW_F(NTIME-ITSTIME+1)
      !             ewrite(3,*) 'II,NOCVA,k1+II-1',II,NOCVA,k1+II-1
        NOD=IMEM_CVA(K1+II-1)
     if(nod.eq.0) then
     ewrite(3,*) 'NOD,ii',NOD,ii
     stop 3333
     endif
!        ewrite(3,*) 'NOD,NONODS',NOD,NONODS
        G(K1+II-1) = WADJ(NOD)
        MAXGRAD = MAX(MAXGRAD,ABS(G(K1+II-1)) )
     ENDDO
  ENDIF
  !20           continue
  ewrite(3,*) 'MAXGRAD inside grad_BC',MAXGRAD
  IF( (ITSTIME+1).LE.NTIME) THEN
    MRECORD(ITSTIME+1) = K2-LSGRAD                 &
     -NOBCU_F(NTIME-ITSTIME)-NOBCV_F(NTIME-ITSTIME)-NOBCW_F(NTIME-ITSTIME)+1
  ENDIF

END SUBROUTINE GRADIENT_BC_SIMPLE

!========================================


  SUBROUTINE GRADIENT_BC_FS1(NONODS,NOBCH_F,HADJ,MAXGRAD,D3,           &
     ITSTIME,MAXREC,MRECORD,NOCVA,G,IMEM_CVA,NTIME,GCOVARI,            &
     GOTFRE)          
!     HADJ,FSNDID,GOTFRE )
   use FLDebug
  IMPLICIT NONE  
  !*************************************************************************
  ! This subroutine is used to calculate the gradient of the Cost function
  !*************************************************************************
  ! Find gradient and functional...
  REAL      ::MAXGRAD
  LOGICAL   ::D3
  INTEGER   ::NONODS,NTIME
  INTEGER   ::NOBCH_F(NTIME)
  REAL      ::HADJ(NONODS),FSNDID(NONODS)
  INTEGER   ::NOCVA       ! number of the control variables
  REAL               ::G(NOCVA)    ! the vector of the gradient
  INTEGER   ::IMEM_CVA(NOCVA)
  INTEGER   ::ITSTIME,MAXREC
  INTEGER   ::MRECORD(MAXREC)
  INTEGER   ::LSGRAD
  LOGICAL  ::GOTFRE
  REAL               ::GCOVARI(NOCVA)

  ! Local....
  INTEGER ::K1,K2                 ! the records for writing to G
  INTEGER ::II,NOD

  IF(ITSTIME.EQ.1) THEN
     MAXGRAD = 0.0
  ENDIF

! Free surface.....
!-----------------
        open(1,file='grd2.dat') 
           write(1,*) NOBCH_F(NTIME-ITSTIME+1)
           write(1,*) (HADJ(II),II=1,NONODS)
        close(1)

  LSGRAD = NOBCH_F(NTIME-ITSTIME+1)

  ewrite(3,*) 'NOCVA in SUB_grad and free surface',NOCVA


  K1 = MRECORD(ITSTIME)
  K2 = K1+NOBCH_F(NTIME-ITSTIME+1)-1
  ewrite(3,*) 'K1,K2',K1,K2

  DO II=1,NOBCH_F(NTIME-ITSTIME+1)
       NOD=IMEM_CVA(K1+II-1)
       ewrite(3,*) 'nod',nod
       G(K1+II-1) = HADJ(NOD)
       MAXGRAD = MAX(MAXGRAD,ABS(G(K1+II-1)) )
    if(itstime.eq.ntime) then
       G(K1+II-1) = 0.0  ! not adjust the initial condition
    endif
  ENDDO
 

  ewrite(3,*) 'MAXGRAD inside grad_BC_FS',MAXGRAD
  IF( (ITSTIME+1).LE.NTIME) THEN
    MRECORD(ITSTIME+1) = K2-LSGRAD                 &
     -NOBCH_F(NTIME-ITSTIME)+1
  ENDIF

END SUBROUTINE GRADIENT_BC_FS1

!========================================


  SUBROUTINE GRADIENT_BC_FS2(NONODS,NONODS_F,NOBCH_F,NTIME,NOBCU,HADJ,U_F,                &
                            BCU1,BCU2,MAXGRAD,D3,                                         &
                            ITSTIME,MAXREC,MRECORD,NOCVA,G,IMEM_CVA,                      &
                            REGU,LAMDA,GCOVARI,GOTFRE)          
!     HADJ,FSNDID,GOTFRE )
   use FLDebug
  IMPLICIT NONE  
  !*************************************************************************
  ! This subroutine is used to calculate the gradient of the Cost function
  !*************************************************************************
  ! Find gradient and functional...
  REAL      ::MAXGRAD
  LOGICAL   ::D3
  INTEGER   ::NONODS,NONODS_F,NTIME,NOBCU
  INTEGER   ::NOBCH_F(NTIME)
  REAL      ::HADJ(NONODS),U_F(NONODS_F)
  REAL      ::BCU1(NOBCU)
  integer   ::BCU2(NOBCU)
  INTEGER   ::NOCVA       ! number of the control variables
  REAL               ::G(NOCVA)    ! the vector of the gradient
  INTEGER   ::IMEM_CVA(NOCVA)
  INTEGER   ::ITSTIME,MAXREC
  INTEGER   ::MRECORD(MAXREC)
  INTEGER   ::LSGRAD
  INTEGER   ::REGU
  REAL      ::LAMDA
  LOGICAL   ::GOTFRE
  REAL               ::GCOVARI(NOCVA)

  ! Local....
  INTEGER ::K1,K2                 ! the records for writing to G
  INTEGER ::II,JJ,NOD_F,NOD

  IF(ITSTIME.EQ.1) THEN
     MAXGRAD = 0.0
  ENDIF

! Free surface.....
!-----------------

  LSGRAD = NOBCH_F(NTIME-ITSTIME+1)

  ewrite(3,*) 'NOCVA in SUB_grad and free surface',NOCVA


  K1 = MRECORD(ITSTIME)
  K2 = K1+NOBCH_F(NTIME-ITSTIME+1)-1
  ewrite(3,*) 'K1,K2',K1,K2

  DO II=1,NOBCH_F(NTIME-ITSTIME+1)
       NOD_F=IMEM_CVA(K1+II-1)
       NOD=BCU2(II) 
     IF(REGU.EQ.1) THEN
       G(K1+II-1) = HADJ(NOD)*U_F(NOD_F)+LAMDA*HADJ(NOD)
     ELSE
       G(K1+II-1) = HADJ(NOD)*U_F(NOD_F)
     ENDIF
       MAXGRAD = MAX(MAXGRAD,ABS(G(K1+II-1)) )

  ENDDO

  ewrite(3,*) 'MAXGRAD inside grad_BC_FS',MAXGRAD
  IF( (ITSTIME+1).LE.NTIME) THEN
    MRECORD(ITSTIME+1) = K2-LSGRAD                 &
     -NOBCH_F(NTIME-ITSTIME)+1
  ENDIF

END SUBROUTINE GRADIENT_BC_FS2
