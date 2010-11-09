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
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
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

#ifdef TESTING_HALGET

SUBROUTINE ResetHalo(V, len)
  IMPLICIT NONE
  INTEGER, INTENT(IN)::len
  REAL, INTENT(OUT)::V(len)
  INTEGER I
  DO I=1, len
     V(I) = -1.0
  END DO
  RETURN
END SUBROUTINE ResetHalo

LOGICAL FUNCTION HaloOk(halo, haloLen, NProcs, MyRank)
  use fldebug
  IMPLICIT NONE
  INTEGER, INTENT(IN)::haloLen, NProcs, MyRank
  REAL, INTENT(IN)::halo((NProcs-1)*haloLen)
  INTEGER I,J,POS
  HaloOk = .TRUE.
  POS = 1
  DO I=0, NProcs-1
     IF (I.NE.MyRank) THEN
        DO J=1,haloLen
           IF (INT(halo(POS)+0.5).NE.MyRank) THEN
              HaloOk = .FALSE.
              GOTO 10
           END IF
           POS = POS+1
        END DO
     END IF
  END DO

10 CONTINUE

  IF(.NOT.HaloOk) THEN
     DO I=0, NProcs-1
        DO J=1,haloLen
           ewrite(3,*)  MyRank, I, halo(I*haloLen+J)
        END DO
     END DO
  END IF 

  RETURN
END FUNCTION HaloOk

PROGRAM MAIN
  use FLDebug
  use parallel_tools
  IMPLICIT NONE
  EXTERNAL HaloOk
  LOGICAL HaloOk
  REAL, ALLOCATABLE, DIMENSION(:)::Array, Buffer
  INTEGER, ALLOCATABLE, DIMENSION(:)::IMEM
  INTEGER, PARAMETER::blockLen=1, stride=1, fieldCnt=1, NTrials=1000
  INTEGER, ALLOCATABLE, DIMENSION(:)::ATOSEN, Gather, ATOREC, Scatter

  INTEGER SizeHalo
  INTEGER MyRank, NProcs, IERROR
  INTEGER NCOLGA, NSCATE
  INTEGER I, J, N, POS
  REAL*8 Time0, TimeNB1(NTrials), TimeNB2(NTrials)
  REAL*8 TimeNB1_mean, TimeNB2_mean
  REAL*8 TimeNB1_variance, TimeNB2_variance
  REAL*8 TimeNB1_err, TimeNB2_err
  INCLUDE 'mpif.h'

  CALL MPI_INIT(IERROR)

  NProcs = GetNProcs()
  MyRank = GetRank()

  IF(MyRank.EQ.0) THEN
     ewrite(2,FMT='(A9,A19,A9,A19,A9)') "SizeHalo", &
          "Non-Blocking 1(ms)", "error",           &
          "Non-Blocking 2(ms)", "error"
  END IF

  ALLOCATE( ATOSEN(NProcs+1) )
  ALLOCATE( ATOREC(NProcs+1) )

  SizeHalo=64
  DO N=1, 13
     SizeHalo = SizeHalo*2
 
     ATOSEN(1) = 1
     ATOREC(1) = 1
     DO I=1, NProcs
        IF(I.EQ.(MyRank+1)) THEN
           ATOSEN(I+1) = ATOSEN(I)
           ATOREC(I+1) = ATOREC(I)
        ELSE
           ATOSEN(I+1) = ATOSEN(I)+SizeHalo
           ATOREC(I+1) = ATOREC(I)+SizeHalo
        END IF
     END DO

     NCOLGA = ATOSEN(NProcs + 1) - 1
     NSCATE = ATOREC(NProcs + 1) - 1

     IF( ALLOCATED(Array) ) DEALLOCATE(Array)   
     ALLOCATE( Array((2*NProcs-1)*SizeHalo) )
     DO I=0, NProcs-1
        DO J=1,SizeHalo
           Array(I*SizeHalo+J) = REAL(I)
        END DO
     END DO

     IF( ALLOCATED(Buffer) ) DEALLOCATE(Buffer)   
     ALLOCATE( Buffer(NProcs*SizeHalo*100) )

     IF( ALLOCATED(Gather) ) DEALLOCATE(Gather)   
     ALLOCATE( Gather((NProcs-1)*SizeHalo) )

     IF( ALLOCATED(Scatter) ) DEALLOCATE(Scatter)   
     ALLOCATE( Scatter((NProcs-1)*SizeHalo) )
     
     POS = 1
     DO I=0, NProcs-1
        DO J=1, SizeHalo
           IF(I.NE.MyRank) THEN
              Gather(POS)  = I*SizeHalo + J
              Scatter(POS) = NProcs*SizeHalo + POS
              POS = POS+1
           END IF
        END DO
     END DO
     
     DO J=1, NTrials
        CALL ResetHalo(Array(NProcs*SizeHalo+1), (NProcs-1)*SizeHalo)
        
        ! Test 1
        Time0 = MPI_WTIME()        
        CALL HalgetNB(Array, blockLen, stride, fieldCnt, &
             ATOSEN, Gather, ATOREC, Scatter)        
        TimeNB1(J) = MPI_WTIME() - Time0
        
        IF(.NOT.HaloOk(Array(NProcs*SizeHalo+1), SizeHalo, NProcs, MyRank)) THEN
           ewrite(-1,*)  "HalgetNB() is foobar"
           CALL MPI_ABORT(MPI_COMM_WORLD, MPI_ERR_OTHER, IERROR)
        END IF
        CALL ResetHalo(Array(NProcs*SizeHalo+1), (NProcs-1)*SizeHalo)
        
        ! Test 2
        Time0 = MPI_WTIME()
        CALL HalgetNB_simple(Array, 1, ATOSEN, Gather, ATOREC, Scatter)
        TimeNB2(J) = MPI_WTIME() - Time0
        
        IF(.NOT.HaloOk(Array(NProcs*SizeHalo+1), SizeHalo, NProcs, MyRank)) THEN
           ewrite(-1,*)  "HalgetNB_simple() is foobar"
           CALL MPI_ABORT(MPI_COMM_WORLD, MPI_ERR_OTHER, IERROR)
        END IF
     END DO

     TimeNB1_mean = 0.0
     TimeNB2_mean = 0.0
     DO J=1, NTrials
        TimeNB1_mean = TimeNB1_mean + TimeNB1(J)/NTrials
        TimeNB2_mean = TimeNB2_mean + TimeNB2(J)/NTrials
     END DO
     
     TimeNB1_variance = 0.0
     TimeNB2_variance = 0.0
     DO J=1, NTrials
        TimeNB1_variance = TimeNB1_variance + ((TimeNB1_mean - TimeNB1(J))**2)/NTrials
        TimeNB2_variance = TimeNB2_variance + ((TimeNB2_mean - TimeNB2(J))**2)/NTrials
     END DO
     
     TimeNB1_err = SQRT(TimeNB1_variance)
     TimeNB2_err = SQRT(TimeNB2_variance)
     
     IF(MyRank.EQ.0) THEN
        ewrite(2,FMT='(I9,I19,I9,I19,I9)') SizeHalo, &
             INT(1000.0*TimeNB1_mean+0.5), INT(1000.0*TimeNB1_err+0.5), &
             INT(1000.0*TimeNB2_mean+0.5), INT(1000.0*TimeNB2_err+0.5)
     END IF
     
  END DO
  CALL MPI_FINALIZE(IERROR)

END PROGRAM MAIN

#endif

! this may be a duplicate - need to compare

#ifdef TESTING_HALGET

SUBROUTINE ResetHalo(V, len)
  IMPLICIT NONE
  INTEGER, INTENT(IN)::len
  REAL, INTENT(OUT)::V(len)
  INTEGER I
  DO I=1, len
     V(I) = -1.0
  END DO
  RETURN
END SUBROUTINE ResetHalo

LOGICAL FUNCTION HaloOk(halo, haloLen, NProcs, MyRank)
  use fldebug
  IMPLICIT NONE
  INTEGER, INTENT(IN)::haloLen, NProcs, MyRank
  REAL, INTENT(IN)::halo((NProcs-1)*haloLen)
  INTEGER I,J,POS
  HaloOk = .TRUE.
  POS = 1
  DO I=0, NProcs-1
     IF (I.NE.MyRank) THEN
        DO J=1,haloLen
           IF (INT(halo(POS)+0.5).NE.MyRank) THEN
              HaloOk = .FALSE.
              GOTO 10
           END IF
           POS = POS+1
        END DO
     END IF
  END DO

10 CONTINUE

  IF(.NOT.HaloOk) THEN
     DO I=0, NProcs-1
        DO J=1,haloLen
           ewrite(3,*)  MyRank, I, halo(I*haloLen+J)
        END DO
     END DO
  END IF 

  RETURN
END FUNCTION HaloOk

PROGRAM MAIN
  use fldebug
  use parallel_tools
  IMPLICIT NONE
  EXTERNAL HaloOk
  LOGICAL HaloOk
  REAL, ALLOCATABLE, DIMENSION(:)::Array, Buffer
  INTEGER, ALLOCATABLE, DIMENSION(:)::IMEM
  INTEGER, PARAMETER::blockLen=1, stride=1, fieldCnt=1, NTrials=1000
  INTEGER, ALLOCATABLE, DIMENSION(:)::ATOSEN, Gather, ATOREC, Scatter

  INTEGER SizeHalo
  INTEGER MyRank, NProcs, IERROR
  INTEGER NCOLGA, NSCATE
  INTEGER I, J, N, POS
  REAL*8 Time0, TimeNB1(NTrials), TimeNB2(NTrials)
  REAL*8 TimeNB1_mean, TimeNB2_mean
  REAL*8 TimeNB1_variance, TimeNB2_variance
  REAL*8 TimeNB1_err, TimeNB2_err
  INCLUDE 'mpif.h'

  CALL MPI_INIT(IERROR)

  NProcs = GetNProcs()
  MyRank = GetRank()

  IF(MyRank.EQ.0) THEN
     ewrite(3,FMT='(A9,A19,A9,A19,A9)') "SizeHalo", &
          "Non-Blocking 1(ms)", "error",           &
          "Non-Blocking 2(ms)", "error"
  END IF

  ALLOCATE( ATOSEN(NProcs+1) )
  ALLOCATE( ATOREC(NProcs+1) )

  SizeHalo=64
  DO N=1, 13
     SizeHalo = SizeHalo*2
 
     ATOSEN(1) = 1
     ATOREC(1) = 1
     DO I=1, NProcs
        IF(I.EQ.(MyRank+1)) THEN
           ATOSEN(I+1) = ATOSEN(I)
           ATOREC(I+1) = ATOREC(I)
        ELSE
           ATOSEN(I+1) = ATOSEN(I)+SizeHalo
           ATOREC(I+1) = ATOREC(I)+SizeHalo
        END IF
     END DO

     NCOLGA = ATOSEN(NProcs + 1) - 1
     NSCATE = ATOREC(NProcs + 1) - 1

     IF( ALLOCATED(Array) ) DEALLOCATE(Array)   
     ALLOCATE( Array((2*NProcs-1)*SizeHalo) )
     DO I=0, NProcs-1
        DO J=1,SizeHalo
           Array(I*SizeHalo+J) = REAL(I)
        END DO
     END DO

     IF( ALLOCATED(Buffer) ) DEALLOCATE(Buffer)   
     ALLOCATE( Buffer(NProcs*SizeHalo*100) )

     IF( ALLOCATED(Gather) ) DEALLOCATE(Gather)   
     ALLOCATE( Gather((NProcs-1)*SizeHalo) )

     IF( ALLOCATED(Scatter) ) DEALLOCATE(Scatter)   
     ALLOCATE( Scatter((NProcs-1)*SizeHalo) )
     
     POS = 1
     DO I=0, NProcs-1
        DO J=1, SizeHalo
           IF(I.NE.MyRank) THEN
              Gather(POS)  = I*SizeHalo + J
              Scatter(POS) = NProcs*SizeHalo + POS
              POS = POS+1
           END IF
        END DO
     END DO
     
     DO J=1, NTrials
        CALL ResetHalo(Array(NProcs*SizeHalo+1), (NProcs-1)*SizeHalo)
        
        ! Test 1
        Time0 = MPI_WTIME()        
        CALL HalgetNB(Array, blockLen, stride, fieldCnt, &
             ATOSEN, Gather, ATOREC, Scatter)        
        TimeNB1(J) = MPI_WTIME() - Time0
        
        IF(.NOT.HaloOk(Array(NProcs*SizeHalo+1), SizeHalo, NProcs, MyRank)) THEN
           ewrite(3,*)  "HalgetNB() is foobar"
           CALL MPI_ABORT(MPI_COMM_WORLD, MPI_ERR_OTHER, IERROR)
        END IF
        CALL ResetHalo(Array(NProcs*SizeHalo+1), (NProcs-1)*SizeHalo)
        
        ! Test 2
        Time0 = MPI_WTIME()
        CALL HalgetNB_simple(Array, 1, ATOSEN, Gather, ATOREC, Scatter)
        TimeNB2(J) = MPI_WTIME() - Time0
        
        IF(.NOT.HaloOk(Array(NProcs*SizeHalo+1), SizeHalo, NProcs, MyRank)) THEN
           ewrite(3,*)  "HalgetNB_simple() is foobar"
           CALL MPI_ABORT(MPI_COMM_WORLD, MPI_ERR_OTHER, IERROR)
        END IF
     END DO

     TimeNB1_mean = 0.0
     TimeNB2_mean = 0.0
     DO J=1, NTrials
        TimeNB1_mean = TimeNB1_mean + TimeNB1(J)/NTrials
        TimeNB2_mean = TimeNB2_mean + TimeNB2(J)/NTrials
     END DO
     
     TimeNB1_variance = 0.0
     TimeNB2_variance = 0.0
     DO J=1, NTrials
        TimeNB1_variance = TimeNB1_variance + ((TimeNB1_mean - TimeNB1(J))**2)/NTrials
        TimeNB2_variance = TimeNB2_variance + ((TimeNB2_mean - TimeNB2(J))**2)/NTrials
     END DO
     
     TimeNB1_err = SQRT(TimeNB1_variance)
     TimeNB2_err = SQRT(TimeNB2_variance)
     
     IF(MyRank.EQ.0) THEN
        ewrite(3,FMT='(I9,I19,I9,I19,I9)') SizeHalo, &
             INT(1000.0*TimeNB1_mean+0.5), INT(1000.0*TimeNB1_err+0.5), &
             INT(1000.0*TimeNB2_mean+0.5), INT(1000.0*TimeNB2_err+0.5)
     END IF
     
  END DO
  CALL MPI_FINALIZE(IERROR)

END PROGRAM MAIN

#endif
