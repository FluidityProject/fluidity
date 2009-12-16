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
       SUBROUTINE GRADIENT_BC(BIGM1,BCU1,BCU2,BCV1,BCV2,                            &
                 BCW1,BCW2,UADJ,VADJ,WADJ,                                           &
                 FINDRM,COLM,                                                        &
                 NONODS,NBIGM,NCOLM,                                                 &
                 NOBCU,NOBCV,NOBCW,MAXGRAD,D3,ITSTIME,MAXREC,MRECORD,NOCVA,G,        &
                 ML,THETA,DT,ACCTIM,IMEM_CVA)

!*************************************************************************
! This subroutine is used to calculate the gradient of the Cost function
!*************************************************************************
! Find gradient and functional...

     use FLDebug
        IMPLICIT NONE
  REAL      ::MAXGRAD
  LOGICAL   ::D3
  INTEGER   ::NONODS,NBIGM,NCOLM
  INTEGER   ::NOBCU,NOBCV,NOBCW            !passed dwon from forward model
  REAL      ::BIGM1(NBIGM)
  INTEGER    ::FINDRM(NONODS+1),COLM(NCOLM)
  REAL    ::UADJ(NONODS),VADJ(NONODS),WADJ(NONODS)
  REAL      ::BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
  INTEGER   ::BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW) !passed dwon from forward model
  INTEGER   ::NOCVA       ! number of the control variables
  REAL    ::G(NOCVA)    ! the vector of the gradient
  INTEGER   ::IMEM_CVA(NOCVA)
  INTEGER   ::ITSTIME,MAXREC
  INTEGER    ::MRECORD(MAXREC)
  REAL    ::ML(NONODS)
  INTEGER    ::LSGRAD
  REAl    ::THETA,DT,ACCTIM
! Local....
  INTEGER ::K1,K2        ! the records for writing to G
  INTEGER ::II,NOD

  IF(ITSTIME.EQ.1) THEN
  MAXGRAD = 0.0
  ENDIF

! Local total size of G at each time step
  LSGRAD = NOBCU+NOBCV
  IF(D3) THEN
    LSGRAD = LSGRAD + NOBCW
  ENDIF
!  LSGRAD = NOBCU
  
  ewrite(3,*) 'NOCVA in SUB_grad',NOCVA
  K1 = MRECORD(ITSTIME)
  K2 = K1+NOBCU-1
   ewrite(3,*) 'K1,K2',K1,K2
        DO II=1,NOBCU
          NOD=IMEM_CVA(K1+II-1)
  ewrite(3,*) 'NOD,UADJ(NOD)',NOD,UADJ(NOD)
     G(K1+II-1) = UADJ(NOD)
     MAXGRAD = MAX(MAXGRAD,ABS(G(K1+II-1)) )
   ENDDO

  IF(ITSTIME.EQ.1) THEN
     OPEN(1,FILE='gradient_time.dat',STATUS='REPLACE')
  ELSE
     OPEN(1,FILE='gradient_time.dat',POSITION='APPEND')
  ENDIF
  write(1,*) 'ITSTIME=',ITSTIME
  write(1,*) (G(K1+II-1),II=1,NOBCU)
  CLOSE(1)


!  CALL G_BC(BIGM1,BCU1,BCU2,UADJ,VADJ,WADJ,                              &
!             FINDRM,COLM,NONODS,NBIGM,NCOLM,NOBCU,MAXGRAD,ITSTIME,             &
!             MAXREC,MRECORD,NOCVA,G,ML,THETA,DT,LSGRAD,ACCTIM,K1,K2,1)

  K1 = K2+1             
  K2 = K1+NOBCV-1
  ewrite(3,*) 'K1,K2',K1,K2
  ewrite(3,*) 'BCV2(136)',BCV2(136)
        DO II=1,NOBCV
           NOD=IMEM_CVA(K1+II-1)
     G(K1+II-1) = VADJ(NOD)
     MAXGRAD = MAX(MAXGRAD,ABS(G(K1+II-1)) )
   ENDDO
!  CALL G_BC(BIGM1,BCV1,BCV2,UADJ,VADJ,WADJ,                               &
!             FINDRM,COLM,NONODS,NBIGM,NCOLM,NOBCV,MAXGRAD,ITSTIME,              &
!             MAXREC,MRECORD,NOCVA,G,ML,THETA,DT,LSGRAD,ACCTIM,K1,K2,2)

  IF(D3) THEN
    K1 = K2+1           
    K2 = K1+NOBCW-1
    ewrite(3,*) 'K1,K2',K1,K2
        DO II=1,NOBCW
!  ewrite(3,*) 'II,NOCVA,k1+II-1',II,NOCVA,k1+II-1
            NOD=IMEM_CVA(K1+II-1)
     G(K1+II-1) = WADJ(NOD)
     MAXGRAD = MAX(MAXGRAD,ABS(G(K1+II-1)) )
   ENDDO
!    CALL G_BC(BIGM1,BCW1,BCW2,UADJ,VADJ,WADJ,                             &
!             FINDRM,COLM,NONODS,NBIGM,NCOLM,NOBCW,MAXGRAD,ITSTIME,              &
!             MAXREC,MRECORD,NOCVA,G,ML,THETA,DT,LSGRAD,ACCTIM,K1,K2,3)
  ENDIF
!20  continue
  ewrite(3,*) 'MAXGRAD inside grad_BC',MAXGRAD
  MRECORD(ITSTIME+1) = K2-2*LSGRAD+1
!  MRECORD(ITSTIME+2) = MRECORD(ITSTIME+1)+2*(NOBCU+NOBCV+NOBCW)+1

  END SUBROUTINE GRADIENT_BC

  SUBROUTINE G_BC(BIGM1,BCU1,BCU2,UADJ,VADJ,WADJ,                        &
             FINDRM,COLM,NONODS,NBIGM,NCOLM,NOBCU,MAXGRAD,ITSTIME,             &
             MAXREC,MRECORD,NOCVA,G,ML,THETA,DT,LSGRAD,ACCTIM,K1,K2,TU)

     use FLDebug
        IMPLICIT NONE
  REAL    ::MAXGRAD
  LOGICAL ::D3
  INTEGER ::NOBCU,NONODS,NBIGM,NCOLM
  REAL    ::BIGM1(NBIGM)
  INTEGER  ::FINDRM(NONODS+1),COLM(NCOLM)
  REAL  ::UADJ(NONODS),VADJ(NONODS),WADJ(NONODS)
  REAL  ::BCU1(NOBCU)
  INTEGER ::BCU2(NOBCU)
  INTEGER ::NOCVA       ! number of the control variables
  REAL  ::G(NOCVA)    ! the vector of the gradient
  INTEGER ::ITSTIME,MAXREC
  INTEGER  ::MRECORD(MAXREC)
  INTEGER  ::LSGRAD
  REAL  ::ML(NONODS)
  REAl  ::THETA,DT,ACCTIM
  
! Local........
        REAL,DIMENSION(:),ALLOCATABLE:: UGRAD
  REAL  ::DSDM
  INTEGER ::ICOL,COUNT,COUNT2,NOD,NOD2
  INTEGER ::II,K1,K2
  INTEGER ::TU       ! TU=1,2,3 for U,V,W respectively

  ALLOCATE(UGRAD(NOBCU))
        UGRAD(1:NOBCU) = 0.0

! (ML+DT*THETA*A)( U(n+1)-(U(n) )/DT = -A*U(n)+f
! (ML+DT*THETA*A)*U(n+1)/DT = (ML+DT*THETA*A)*U(n)/DT - A*U(n)+f
! BIGM1 = ML+DT*THETA*A; DT*A = (BIGM1-ML)/THETA
! BIGM1*U(n+1)/DT = BIGM1*U(n)/DT - A*U(n)+f
! PARTIAL S/PARTIAL M_l = DSDM = (BIGM1-A*DT)/DT
! PARTIAL S/PARTIAL M_l = DSDM = ( (THETA-1)*BIGM1+ML )/(DT*THETA)


! CALCULATE PARTIAL S/PARTIAL M_l
         DO II=1,NOBCU
           NOD=BCU2(II)

! find coln NOD in BIGM1
! search for the coln in the neighbouring nodes 
! to NOD (ie. next to the boundary)
!  ewrite(3,*)  'size of the FINDRM',size(FINDRM)
!  ewrite(3,*)  'NOD'
!  ewrite(3,*) NOD
!  ewrite(3,*) 'FINDRM(NOD)',FINDRM(NOD)
           DO COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
             NOD2=COLM(COUNT)
!  ewrite(3,*) 'COLM(COUNT)',COLM(COUNT)
! Seach BIGM1 matrix..and find PARTIAL S/PARTIAL M_l: SUGRAD....
! For block1 in BIGM --,mutiplied by UADJ
! locate dignal...
            IF(BIGM1((TU-1)*3*NCOLM+COUNT).LT.1.0E10) THEN
             IF((NOD.EQ.NOD2).AND.(TU.EQ.1)) THEN
                     UGRAD(II)=UGRAD(II)                                  
             ELSE
               IF(ACCTIM.NE.0.0) THEN
! the gradient at the last time step....
!                 ewrite(3,*) 'ACCTIM,II,LSGRAD',ACCTIM,II,LSGRAD
                 DSDM = ( (THETA-1)*BIGM1((TU-1)*3*NCOLM+COUNT)       &
                             )/(DT*THETA)
                       G(K1+II-1+LSGRAD)=G(K1+II-1+LSGRAD)+                &
                                        DSDM*UADJ(NOD2)
               ENDIF                      
                       UGRAD(II)=UGRAD(II)                                 &
                           -BIGM1((TU-1)*3*NCOLM+COUNT)*UADJ(NOD2)/DT                   

             ENDIF
!  ewrite(3,*) 'BIGM1((TU-1)*3*NCOLM+COUNT)',BIGM1((TU-1)*3*NCOLM+COUNT)
            ENDIF
!  goto 10

! For block2 in BIGM --, mutiplied by VADJ
! locate diagal...
            IF(BIGM1((TU-1)*3*NCOLM+NCOLM+COUNT).LT.1.0E10) THEN
             IF((NOD.EQ.NOD2).AND.(TU.EQ.2)) THEN
                     UGRAD(II)=UGRAD(II)                                   
             ELSE
               IF(ACCTIM.NE.0.0) THEN
! the gradient at the last time step....
                 DSDM = ( (THETA-1)*BIGM1((TU-1)*3*NCOLM+NCOLM+COUNT)     &
                             )/(DT*THETA)
                       G(K1+II-1+LSGRAD)=G(K1+II-1+LSGRAD)+                  &
                                        DSDM*VADJ(NOD2)
               ENDIF                      
                       UGRAD(II)=UGRAD(II)                                   &
                            -BIGM1((TU-1)*3*NCOLM+NCOLM+COUNT)*VADJ(NOD2)/DT
             ENDIF
            ENDIF

! For block3 in BIGM --, mutiplied by WADJ
           IF(D3) THEN
! locate dignal...
            IF(BIGM1((TU-1)*3*NCOLM+2*NCOLM+COUNT).LT.1.0E10) THEN
             IF((NOD.EQ.NOD2).AND.(TU.EQ.3)) THEN
                     UGRAD(II)=UGRAD(II)                                   
             ELSE
               IF(ACCTIM.NE.0.0) THEN
! the gradient at the last time step....
                 DSDM = ( (THETA-1)*BIGM1((TU-1)*3*NCOLM+2*NCOLM+COUNT)     &
                             )/(DT*THETA)
                       G(K1+II-1+LSGRAD)=G(K1+II-1+LSGRAD)+                  &
                                        DSDM*WADJ(NOD2)
               ENDIF                      
                       UGRAD(II)=UGRAD(II)                                   &
                          -BIGM1((TU-1)*3*NCOLM+2*NCOLM+COUNT)*WADJ(NOD2)/DT
             ENDIF
            ENDIF

           ENDIF

!               ENDIF
          END DO

        END DO

!10  continue
  DO II=1,NOBCU
     MAXGRAD = MAX(MAXGRAD,ABS(UGRAD(II)))
!  if(ABS(UGRAD(II)).gt.3.0)                                   &
!            ewrite(3,*) 'MAXGRAD in SUB_grad',ABS(UGRAD(II)),MAXGRAD,II
  END DO

! write to the total vector of the gradient.......
  G(K1:K2) = UGRAD

  DEALLOCATE(UGRAD)

       END SUBROUTINE G_BC

