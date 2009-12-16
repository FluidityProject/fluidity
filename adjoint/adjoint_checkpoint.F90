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

module adjoint_checkpointing
  implicit none
  use FLDebug 

contains

  SUBROUTINE FBcheckpoint(T,STEPS,LTIME,NTIME,LX,NX,RNU,RU,&
       &  RNUADJO,RUADJO,RNUSADJ,RUSADJ,RNUGRAD,RUGRAD,RNUEXAC,RUEXAC,&
       &  BNU,BU,NU0,U0,BNUGRAD,BUGRAD,NUGRAD0,UGRAD0,&
       &  FUNCT,VALFUN,WEIGHT,MAXGRAD,TYPEINVE,GLOITS,LINITS)



    ! ....Main vatirables and adjoint variables...
    REAL LX,T,LTIME,FUNCT,VALFUN
    INTEGER NX,STEPS,TYPEINVE,NTIME
    INTEGER GLOITS,LINITS
    REAL MAXGRAD
    REAL RNU(0:NTIME+1,0:NX+1),RU(0:NTIME+1,0:NX+1)
    REAL RNUADJO(0:NTIME+1,0:NX+1),RUADJO(0:NTIME+1,0:NX+1)
    REAL RNUSADJ(0:NTIME+1,0:NX+1),RUSADJ(0:NTIME+1,0:NX+1)
    REAL RNUEXAC(0:STEPS,0:NX+1),RUEXAC(0:STEPS,0:NX+1)
    REAL RNUGRAD(0:NTIME+1,0:NX+1),RUGRAD(0:NTIME+1,0:NX+1)
    REAL WEIGHT(0:STEPS,0:NX+1)
    LOGICAL ADJOINT,NOADJOINT
    PARAMETER(ADJOINT=.TRUE.,NOADJOINT=.FALSE.)

    !     ..
    !     .. Checkpoint  ..
    !     .. Parameters ..        

    REAL TAKSHT,ADVAN,FSTURN,YUTURN
    PARAMETER (TAKSHT=1,ADVAN=2,FSTURN=3,YUTURN=4)
    REAL RESTRE,TRMATE,ERROR
    INTEGER CHECKUP
    PARAMETER (CHECKUP=64)
    REAL CU(0:NX),CUSTOR(-1:CHECKUP,0:NX),REC(0:NX)
    REAL CNU(0:NX),CNUSTOR(-1:CHECKUP,0:NX)
    REAL BNU(0:STEPS,0:1),BU(0:STEPS,0:1)
    REAL BNUGRAD(0:STEPS,0:1),BUGRAD(0:STEPS,0:1)
    REAL NUGRAD0(0:NX+1),UGRAD0(0:NX+1)
    REAL NU0(0:NX+1),U0(0:NX+1)
    PARAMETER (RESTRE=5,TRMATE=6,ERROR=7)
    !     .. Local Scalars ..
    INTEGER CAPO,CHECK,FINE,INFO,SNAPS,WHATDO
    INTEGER OLDCAPO
    REAL X,DT
    PARAMETER (INFO=3)
    ! ...Running time....
    !      INTEGER FSTEPS1,FSTEPS2,FSTEPS
    !      REAL TIMERATIO
    !     ..
    !     .. External Functions ..
    INTEGER REVOLV,I,J,ADJUST
    EXTERNAL REVOLV
    !.....
    SNAPS = ADJUST(STEPS)

    !     ..
    ewrite(3,FMT=*) 'ENTER:  STEPS,SNAPS,INFO'
    !      READ (*,FMT=*) STEPS,SNAPS,INFO

    ewrite(3,FMT=*) STEPS,SNAPS,INFO
    DT = T/STEPS
    ewrite(3,*) 'DT=',DT

    !      ewrite(3,*) 'BNU',BNU
    !      ewrite(3,*) 'NU0',NU0
    ! ------------------------------------------
    ! if (LINITS.NE.1) then
    !    only run the forward model
    ! else
    !    run the forward and backward model 
    ! endif
    ! ------------------------------------------


    IF(LINITS.NE.1) THEN
       do I = 1,NX
          RNU(0,I) = NU0(I)
          RU(0,I) =  U0(I)
       END DO

       ! Solve the wave eqn...
       VALFUN = 0.0
       do J = 0,STEPS-1
          RNU(1,0) = BNU(J+1,0)
          RNU(1,NX+1) = BNU(J+1,1)
          RU(1,0) = BU(J+1,0)
          RU(1,NX+1) = BU(J+1,1)
          IF (J.NE.0) THEN
             do I =1,NX
                RNU(0,I) = CNU(I)
                RU(0,I) =  CU(I)       
             END DO
          ENDIF
          CALL WAVSOL(RNU,RU,  RNUSADJ,RUSADJ,&
               &            LTIME,LX,NTIME,NX,NOADJOINT)   
          do I =1,NX
             CU(I) = RU(1,I)
             CNU(I) = RNU(1,I)
          END DO

          ! Find functional...
          CALL CFINGRA(RNUGRAD,RUGRAD, RNUADJO, RUADJO,&
               &                    RNU,RU, RNUEXAC,RUEXAC,FUNCT,STEPS,&
               &                    BNUGRAD,BUGRAD,NUGRAD0,UGRAD0,&
               &                    J,NX,NTIME,LX,LTIME,TYPEINVE,WEIGHT)
          VALFUN = VALFUN+FUNCT 
       END DO
       ewrite(3,*) 'GLOITS,FUNCT,VALFUN ', GLOITS,FUNCT,VALFUN
       GOTO 30
    ENDIF

    CAPO = 0
    FINE = STEPS + CAPO
    CHECK = -1
    ! !!!!DO I = 0,NX+1
    do I = 1,NX
       ! Initial condition....
       RNU(0,I) = NU0(I)
       RU(0,I) = U0(I)
       ! give the inital values for the first takeshot....
       CNU(I) = NU0(I)
       CU(I) = U0(I)
    END DO
    !**********************************************************
10  CONTINUE
    OLDCAPO =CAPO
    WHATDO = REVOLV(CHECK,CAPO,FINE,SNAPS,INFO)
    !---------------------------------------------------------
    IF ((WHATDO.EQ.TAKSHT) .AND. (INFO.GT.1)) THEN
       ewrite(3,FMT=9000) CAPO
       !C!!!!!!!!DO I =0,NX
       do I =1,NX
          CUSTOR(CHECK,I) = CU(I)
          CNUSTOR(CHECK,I) = CNU(I)
       END DO
    END IF

    !----------------------------------------------------------
    IF ((WHATDO.EQ.ADVAN) .AND. (INFO.GT.2)) THEN
       ewrite(3,FMT=9010) CAPO
       IF (OLDCAPO.EQ.0) THEN
          !  !!!!!!DO I = 0,NX+1
          do I = 1,NX
             RNU(0,I) = NU0(I)
             RU(0,I) =  U0(I)
          END DO
       ELSE
          do I =1,NX
             RNU(0,I) = CNUSTOR(CHECK,I)
             RU(0,I) =  CUSTOR(CHECK,I)
          END DO
       ENDIF
       ! Solve the wave eqn...
       do J = OLDCAPO,CAPO-1
          RNU(1,0) = BNU(J+1,0)
          RNU(1,NX+1) = BNU(J+1,1)
          RU(1,0) = BU(J+1,0)
          RU(1,NX+1) = BU(J+1,1)
          IF (J.NE.OLDCAPO) THEN
             !  !!!!!!     DO I = 0,NX
             do I =1,NX
                RNU(0,I) = CNU(I)
                RU(0,I) =  CU(I)       
             END DO
          ENDIF
          CALL WAVSOL(RNU,RU,  RNUSADJ,RUSADJ,&
               &            LTIME,LX,NTIME,NX,NOADJOINT)   
          do I =1,NX
             CU(I) = RU(1,I)
             CNU(I) = RNU(1,I)
          END DO
       ENDDO
       !        ewrite(3,*) 'RNU',RNU
    END IF

    !----------------------------------------------------------
    IF ((WHATDO.EQ.FSTURN) .AND. (INFO.GT.2)) THEN
       ewrite(3,FMT=9020) CAPO
       ! initial and boundary conditions......
       RNU(1,0) = BNU(CAPO+1,0)
       RNU(1,NX+1) = BNU(CAPO+1,1)
       RU(1,0) = BU(CAPO+1,0)
       RU(1,NX+1) = BU(CAPO+1,1)
       do I =1,NX
          RNU(0,I) = CNU(I)
          RU(0,I) =  CU(I)       
       END DO
       !  running the forward model......
       CALL WAVSOL(RNU,RU,RNUSADJ,RUSADJ,&
            &            LTIME,LX,NTIME,NX,NOADJOINT)   
       !        ewrite(3,*) 'RNU',RNU
       ! Initianise adjoint variables....
       do I =1,NX
          !           RNUADJO(2,I) = RNU(1,I)-CNUSTOR(CHECK,I)
          !            RUADJO(2,I) =  RU(1,I)- CUSTOR(CHECK,I)
          RNUADJO(2,I) = 0.0
          RUADJO(2,I) = 0.0
       END DO
       ! The first backward running
       ! ....Calculate the residuals...
       CALL CRESSOU( RNUSADJ,RUSADJ, RNU,RU,RNUEXAC,RUEXAC,&
            &                CAPO,STEPS,LX,NTIME,NX,TYPEINVE,WEIGHT)

       !      ewrite(3,*) 'WEIGHT',(WEIGHT(CAPO,I),I=0,NX+1)
       !      ewrite(3,*) 'RNU',(RNU(CAPO,I),I=0,NX+1)
       !      ewrite(3,*) 'RNUEXAC',(RNUEXAC(CAPO,I),I=0,NX+1)
       !      ewrite(3,*) 'RNUSADJ',RNUSADJ
       ! ....Solve adjoint eqn...
       CALL WAVSOL(RNUADJO,RUADJO,  RNUSADJ,RUSADJ,&
            &            LTIME,LX,NTIME,NX,ADJOINT)

       do I =1,NX
          !           CU(I) = RU(1,I)
          !          CNU(I) = RNU(1,I)
          RNUADJO(2,I) = RNUADJO(1,I)
          RUADJO(2,I) =  RUADJO(1,I)
       END DO
       ! Find gradient and functional...
       ewrite(3,*) 'fingra'
       !         MAXGRAD = 0.0
       VALFUN = 0.0
       CALL CFINGRA(RNUGRAD,RUGRAD, RNUADJO, RUADJO,&
            &                    RNU,RU, RNUEXAC,RUEXAC,FUNCT,STEPS,&
            &                    BNUGRAD,BUGRAD,NUGRAD0,UGRAD0,&
            &                    CAPO,NX,NTIME,LX,LTIME,TYPEINVE,WEIGHT)
       VALFUN = VALFUN+FUNCT 
       !      IF (GLOITS.EQ.1.AND.LINITS.EQ.1) THEN
       IF (LINITS.EQ.1) THEN
          MAXGRAD = 0.0
          do I = 0,NTIME+1
             do J = 0,NX+1
                MAXGRAD = MAX(ABS(RNUGRAD(I,J)),ABS(RUGRAD(I,J)),MAXGRAD)
             ENDDO
          ENDDO
       ENDIF

    END IF

    !--------------------------------------------------------------
    IF ((WHATDO.EQ.YUTURN) .AND. (INFO.GT.2)) THEN
       ewrite(3,FMT=9030) CAPO
       ! INitial and boundary condition......
       RNU(1,0) = BNU(CAPO+1,0)
       RNU(1,NX+1) = BNU(CAPO+1,1)
       RU(1,0) = BU(CAPO+1,0)
       RU(1,NX+1) = BU(CAPO+1,1)
       do I =1,NX
          RNU(0,I) = CNU(I)
          RU(0,I) =  CU(I)       
       END DO
       ! running the forward model......
       CALL WAVSOL(RNU,RU,  RNUSADJ,RUSADJ,&
            &            LTIME,LX,NTIME,NX,NOADJOINT)  
       ! Backward running
       ! ....Calculate the residuals...
       CALL CRESSOU( RNUSADJ,RUSADJ, RNU,RU,RNUEXAC,RUEXAC,&
            &                CAPO,STEPS,LX,NTIME,NX,TYPEINVE,WEIGHT)
       !        ewrite(3,*) 'RNUSADJ',RNUSADJ
       ! ....Solve adjoint eqn...
       CALL WAVSOL(RNUADJO,RUADJO,  RNUSADJ,RUSADJ,&
            &            LTIME,LX,NTIME,NX,ADJOINT)
       do I =1,NX
          !           CU(I) = RU(1,I)
          !           CNU(I) = RNU(1,I)
          RNUADJO(2,I) = RNUADJO(1,I)
          RUADJO(2,I) =  RUADJO(1,I)
       END DO
       ! Find gradient and functional...
       ewrite(3,*) 'fingra'
       CALL CFINGRA(RNUGRAD,RUGRAD, RNUADJO, RUADJO,&
            &                    RNU,RU, RNUEXAC,RUEXAC,FUNCT,STEPS,&
            &                    BNUGRAD,BUGRAD,NUGRAD0,UGRAD0,&
            &                    CAPO,NX,NTIME,LX,LTIME,TYPEINVE,WEIGHT)
       !      ewrite(3,*) 'RNUGRAD',RNUGRAD
       !      ewrite(3,*) 'BNUGRAD',BNUGRAD
       !      ewrite(3,*) 'NUGRAD0',NUGRAD0
       VALFUN = VALFUN+FUNCT 
       !      IF (GLOITS.EQ.1.AND.LINITS.EQ.1) THEN
       IF (LINITS.EQ.1) THEN
          do I = 0,NTIME+1
             do J = 0,NX+1
                MAXGRAD = MAX(ABS(RNUGRAD(I,J)),ABS(RUGRAD(I,J)),MAXGRAD)
             ENDDO
          ENDDO
       ENDIF
       ewrite(3,*) '---------------'
       ewrite(3,*) 'GLOITS,FUNCT,VALFUN ', GLOITS,FUNCT,VALFUN
       ewrite(3,*) 'GLOITS,MAXGRAD ',GLOITS,MAXGRAD
150    FORMAT(1X,I4,1X,F10.3)

    END IF

    !--------------------------------------------------------
    IF ((WHATDO.EQ.RESTRE) .AND. (INFO.GT.2)) THEN
       ewrite(3,FMT=9040) CAPO
!!!!!!!DO I = 0,NX
       do I = 1,NX
          CNU(I) = CNUSTOR(CHECK,I)
          CU(I) = CUSTOR(CHECK,I)
       END DO
    END IF

    !-----------------------------------------------------------
    IF (WHATDO.EQ.ERROR) THEN
       ewrite(3,FMT=*) ' irregular termination of treeverse'
       IF (INFO.EQ.10) THEN
          ewrite(3,FMT=*)&
               &          ' number of checkpoints stored exceeds CHEKUP,'
          ewrite(3,FMT=*) ' increase constant CHEKUP and recompile'
       END IF
       IF (INFO.EQ.11) THEN
          ewrite(3,FMT=*) ' number of checkpoints stored = ',&
               &          CHECK + 1,' exceeds SNAPS,'
          ewrite(3,FMT=*) ' ensure SNAPS > 0 and ',&
               &          'increase initial FINE'
       END IF

       IF (INFO.EQ.12) THEN
          ewrite(3,FMT=*) ' error occurs in NUMFRW'
       ENDIF
       IF (INFO.EQ.13) THEN
          ewrite(3,FMT=*) ' enhancement of FINE, SNAPS = ',&
               &          CHECK + 1,'checkpoints stored, increase SNAPS'
       END IF
       IF (INFO.EQ.14) THEN
          ewrite(3,FMT=*) ' number of SNAPS = ',SNAPS,&
               &          ' exceeds CHEKUP,'
          ewrite(3,FMT=*) ' increase constant CHEKUP and recompile'
       END IF
       IF (INFO.EQ.15) THEN
          ewrite(3,FMT=*) ' number of reps exceeds REPSUP, '
          ewrite(3,FMT=*) ' increase constant REPSUP and recompile'
       END IF

    END IF

    !------------------------------------------------------------
    IF ((WHATDO.EQ.TRMATE) .OR. (WHATDO.EQ.ERROR)) THEN
       GO TO 20
    ELSE
       GO TO 10
    END IF
    !****************************************************************
20  CONTINUE
    ewrite(3,*) 'STEPS,SNAPS,INFO',STEPS,SNAPS,INFO
    !      ewrite(3,*) 'RNU',RNU

    IF (LINITS.EQ.1) THEN
       MAXGRAD = 0.0
       do J = 0,NX+1
          MAXGRAD = MAX(ABS(NUGRAD0(J)),ABS(UGRAD0(J)),MAXGRAD)
       ENDDO
       do I=0,NTIME+1
          do J=0,1
             MAXGRAD = MAX(ABS(BNUGRAD(I,J)),ABS(BUGRAD(I,J)),MAXGRAD)
          END DO
       END DO
    ENDIF
30  CONTINUE

9000 FORMAT (' takeshot at',I6)
9010 FORMAT (' advance to',I7)
9020 FORMAT (' firsturn at',I6)
9030 FORMAT (' youturn at',I7)
9040 FORMAT (' restore at',I7)
  END subroutine FBcheckpoint

end module adjoint_checkpointing



