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

module NONLINCG_module
  use FLDebug
  use signals
  use spud
  implicit none
contains

  SUBROUTINE NONLINCG(GLOITS,LINITS,VALFUN,CVA,G,       &
       GOLDR,LINCONVER,NOCVA,GOLD,D,CVAOLD,FSZERO)
    !-----------------------------------------------------------------------
    ! NONLINEAR CONJUGATE GRADIENT METHOD:
    ! This program is to get the minimum of the cost function
    ! and its gradient using the Nonlinear Conjugate Gradient. 
    ! The main procedures:
    ! 1.Give the search direction D(I)=-G(I)+BATA*D(I-1)
    !   where, the initial search direction: D(1) = -G(1)
    !          BATA is calculated using the Polak-Ribiere
    ! 2.Find the minimum of the search step NUVAL using the line search
    ! 3.If the line search converged, CVA(I) = CVA(I-1)+NUVAL*D(I) and start 
    !     another line search from the new search direction (repeat 1-3)
    !   If the line search hasn't converged, repeat 2-3
    ! The main frame
    ! DO GLOITS = 1,NGLOITS ----nonlinear conjugate gradient loop
    !   Start the line search
    !     if LINITS.eq.1
    !     call the subroutine CALCFG to calculate both the cost function
    !                                  and the gradient
    !     else
    !     call the subroutine CALCFG to calculate only the cost function
    !     call the subroutine NONLINCG
    !     if the line search hasn't converged, LINITS =LINITS+1 and
    !        repeat the above line search
    !     else if the line search has converged, then
    !   End the line search
    !  If the nonlinear CG's criterion is satisfied, then jump out the loop
    !  END LOOP---nonlinear conjugate gradient loop   
    ! VARIABLES which are passed down
    ! GLOITS----The current global iteration number
    ! VALFUN----The cost function
    ! CVA(NOCVA)----The control variables 
    !         which we wish to get to the minimum CVA of the cost function
    ! G(NOCVA)----The gradient of the cost function
    ! NOCVA----The dimensional number for CVA,G,D
    ! GOLDR----The golden rate
    ! VARIABLES which are passed back or calculated in this subroutine
    ! D(NOCVA)----The search direction
    ! CVAOLD(NOCVA)----The previous control variables
    ! GOLD(NOCVA)----The previous gradient of the cost function
    ! BATA ---- the cofficient used to calculate the search direction
    ! LINCONVER----The logical variable 
    !              which is used to adjust whether the line search has converged
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER GLOITS,LINITS,NOCVA
    REAL,DIMENSION(NOCVA)::G,CVA
    REAL VALFUN
    REAL GOLDR
    LOGICAL LINCONVER
    REAL FSZERO

    ! Local variables.....
    INTEGER ::NSOR
    REAL,DIMENSION(NOCVA)::GOLD,D,CVAOLD
    REAL    ::SQG,SQGOLD,SQG1,BATA
    REAL    ::NUVAL,MAXCVA,MAXD,INISTEP,MAXCD
    REAL    ::LINERR,CVAFREE
    !        REAL,PARAMETER ::LINERR0=0.5
    REAL,PARAMETER ::LINERR0=0.5e-10
    REAL,PARAMETER ::RELAXP=0.5
    REAL METHOD ! METHOD =1.0, Polak-Ribiere; 2.0, Positive Polak-Ribiere
    !  SAVE GOLD,D,CVAOLD
    SAVE  ::NUVAL
    REAL min_CVA
    INTEGER ::I,min_nod
    REAL MAXF
    INTEGER MAXNSOR,JJ

    METHOD=2.0
    ewrite(3,*) 'NOCVA=',size(CVA)

    LINERR = LINERR0
    IF(GLOITS.GT.3) THEN
       LINERR = LINERR0/10.0
    ELSE
       IF(GLOITS.GT.6) THEN
          LINERR = LINERR0/100.0
       ENDIF
    ENDIF

    ! Find the new search direction D
    !---------------------------------


    IF (LINITS.EQ.1) THEN
       LINCONVER=.FALSE.
       ! Set the initial search direction
       IF (GLOITS.EQ.1) THEN
          !         IF( (GLOITS.EQ.1).OR.( ABS(MOD(GLOITS,5)-0.0).LT.1.0E-6) ) THEN
          !*****
          if(.true.) then
             !*********
             DO I =1,NOCVA
                D(I) = -G(I)
                CVAOLD(I) = CVA(I)
                GOLD(I) = G(I)
             END DO
             !*****
          else

             open(1,file='d1.dat')
             read(1,*) (d(i),i=1,nocva)
             close(1)
             DO I =1,NOCVA
                CVAOLD(I) = CVA(I)
                G(I) = -D(I)
                GOLD(I) = G(I)
             END DO
          endif
          !********

       ELSE  
          ! SQG --- the sum of G(I)^2
          ! SQGOLD --- the previous SQG: the sum of G(I-1)^2
          ! SQR1 --- the sum of G(I)*G(I-1)
          !         ewrite(3,*) 'GOLD',(GOLD(I),I=1,NOCVA)
          SQG = 0.0
          SQGOLD = 0.0
          SQG1 = 0.0
          ! calculate the BATA using Polak-Ribiere formula......
          DO I =1,NOCVA
             SQG = SQG+G(I)*G(I)
             SQGOLD = SQGOLD+GOLD(I)*GOLD(I)
             SQG1 = SQG1+G(I)*GOLD(I)
          END DO

          IF(METHOD.EQ.1) THEN
             BATA = (SQG-SQG1)/SQGOLD
          ELSE
             BATA = MAX( (SQG-SQG1)/SQGOLD, 0.0)
          ENDIF
          ewrite(3,*) 'SQG,SQGOLD,SQG1',SQG,SQGOLD,SQG1
          ewrite(3,*) 'BATA',BATA
          ! calculate the search direction......
          DO I =1,NOCVA
             D(I) = -G(I)+BATA*D(I)
             if(G(I).NE.0.0) then
                ewrite(3,*) 'I,BATA,G(I),D(I)222',I,BATA,G(I),D(I)
             endif
             CVAOLD(I) = CVA(I)
             GOLD(I) = G(I)
          END DO
       ENDIF ! end IF(GLOITS.EQ.1) 

    ENDIF

    if(.false.) then
       ! decent method
       DO I =1,NOCVA
          D(I) = -G(I)
          !        CVAOLD(I) = CVA(I)
          !        GOLD(I) = G(I)
       END DO
    endif

    ! Start the line search......
    ! Calculate the search step NUVAL...
    !-----------------------------

    ! Calculate the INISTEP which is used to find the initial three points
    MAXCVA = 0.0
    MAXD = 0.0
    MAXCD = 0.0
    DO I = 1,NOCVA
       MAXCVA = MAX(ABS(CVA(I)),MAXCVA)
       MAXD = MAX(ABS(D(I)),MAXD)
       MAXCD =MAX(ABS(CVA(I)/D(I)),MAXCD)
    END DO
    !  IF (GLOITS.EQ.1) THEN
    IF( (GLOITS.EQ.1).OR.( ABS(MOD(GLOITS,5)-0.0).LT.1.0E-6) ) THEN
       !          INISTEP = 0.01
       !**headland    INISTEP = 0.2*MAXCVA/MAXD
       !    INISTEP = -0.1*MAXCVA/MAXD
       INISTEP = -0.4*MAXCVA/MAXD
       !          INISTEP = 0.05*MAXCD          
       !     INISTEP = -3.8
       if(have_option(trim("/model/fluids/adjoint/bctest"))) then
          !         INISTEP = 0.2
          INISTEP = -0.2
       endif
       !    INISTEP = 0.0001 ! for healdland
       !    INISTEP = 0.00000000001
    ELSE
       INISTEP = 0.01*MAXCVA/MAXD
       !**    INISTEP = 0.05*MAXCVA/MAXD
       !**headland    INISTEP = 0.1*MAXCVA/MAXD
    ENDIF

    ! CALL the line search subroutine
    ewrite(3,*) 'LINITS=',LINITS
    CALL LINESEARCH(LINITS,VALFUN,NUVAL,GOLDR,    &
         INISTEP,LINCONVER,LINERR)
    ewrite(3,*) 'LINCONVER=',LINCONVER

    ! Do we start another line search along the new search direction?
    ! if not, continue the linear search
    !**
    !****  IF(linits.eq.4) LINCONVER=.true.
    IF (LINCONVER) THEN
       ! Start line search along the new search direction
       ! ... if the nonlinear CG criterion isn't satisfied.......
       DO I = 1,NOCVA
          CVAOLD(I)=CVA(I) 
       END DO
    ELSE
       ! Continue this line search
       ! Renew CVA --- the control variables
       if(linits.eq.3) then
          open(2,file='nuval.dat')
       else
          open(2,file='nuval.dat',position='append')
       endif
       write(2,*) NUVAL,VALFUN,INISTEP,MAXCVA,MAXD
       !***** for the free surface
       if(.false.) then
10        continue
          min_CVA =1.0e12
          min_nod = 0
          DO I = 1,NOCVA
             !       CVA(I) = CVAOLD(I)+NUVAL*D(I)
             !       CVAFREE = CVAOLD(I)+NUVAL*D(I)
             CVAFREE = CVAOLD(I)+RELAXP*NUVAL*D(I)
             IF(CVAFREE.LT.min_CVA) THEN
                min_CVA = CVAFREE
                min_nod = I
             ENDIF
          END DO
          if ( (min_CVA+FSZERO).lt.0.4*FSZERO) then
             NUVAL = (-0.6*FSZERO - CVAOLD(min_nod))/(RELAXP*D(min_nod))
             write (2,*) 'the previous nuval is too larger, reduced to', NUVAL
             GOTO 10
          endif

          close(2)

       endif
       !************
       DO I = 1,NOCVA
          CVA(I) = CVAOLD(I)+RELAXP*NUVAL*D(I)
       END DO
    ENDIF
    open(1,file='d.dat')
    write(1,*) 'RELAXP,NUVAL',RELAXP,NUVAL
    write(1,*) (d(i),i=1,nocva)
    close(1)
    open(1,file='cvaold.dat')
    write(1,*) 'RELAXP,NUVAL',RELAXP,NUVAL
    write(1,*) (cvaold(i),i=1,nocva)
    close(1)
    open(1,file='cva.dat')
    write(1,*) 'RELAXP,NUVAL',RELAXP,NUVAL
    write(1,*) (cva(i),i=1,nocva)
    close(1)

200 CONTINUE

    IF (LINCONVER) THEN
       ! Start line search along the new search direction
       ! ... if the nonlinear CG criterion isn't satisfied.......
       DO I = 1,NOCVA
          CVAOLD(I)=CVA(I) 
       END DO
       LINCONVER=.true. 
       LINITS =1
       NUVAL = 0.0
    ENDIF

    RETURN
  END SUBROUTINE NONLINCG

  SUBROUTINE LINESEARCH(LINITS,VALFUN,NUVAL,  &
       GOLDR,INISTEP,LINCONVER,LINERR)
    !------------------------------------------------------------------------
    ! This subroutine calculates the search step using the line search
    ! The main procedures
    ! 1. Set the initial three points
    ! 2. Start the line search----
    !        fit a quadratic polynomial to the cost function evaluated 
    !        at the three points, then find the new point where the
    !        cost function is minimum
    ! 3. If the line search dosen't satisfy the criterion, the new point
    !       replaces one of the previous points, then start the new line
    !       search
    ! VARIABLES:
    ! LINITS ---- the current intertion number of the line search
    ! VALFUN ---- the cost function
    ! NUVAL ---- the new position (search step) where the function is minimum
    ! GOLDR ---- the gold rate
    ! INISTEP ---- parameter used to calculate the initial three points
    ! LINCONVER ---- the logical parameter to adjust whether the line search 
    !                has converged
    ! Note: When LINITS =1,2,3, we search the initial three points. The line
    !       search starts from LINITS =3
    !-------------------------------------------------------------------------

    ! INISTEP=initial step size or NUVAL values
    IMPLICIT NONE
    INTEGER LINITS
    REAL GOLDR
    REAL VALFUN
    REAL NUVAL,INISTEP
    REAL MINNU,MAXNU
    PARAMETER (MINNU=-1E20,MAXNU=1E20)
    LOGICAL LINCONVER
    ! Local variables.......
    INTEGER NSOR,I
    PARAMETER(NSOR=3)
    REAL NU(NSOR),F(NSOR)
    REAL GOLDR2,LINERR,VALFUN0
    !  PARAMETER(LINERR=0.002)
    !  PARAMETER(LINERR0=0.5)
    !  PARAMETER(LINERR=0.05)
    LOGICAL GOTNUV
    SAVE F,NU,GOLDR2,VALFUN0

    ! Conversion ?
    LINCONVER = .FALSE.


    ! If LINITS=1,2,3,then choose the initial three points/positions.......

    ! fing the first point.......
    IF (LINITS.EQ.1) THEN
       ! Set NU(1) = 0
       NU(1) = 0.0
       ! Set F(1) =VALFUN
       F(1)=VALFUN
       ! Set LINITS =2
       LINITS = LINITS+1
    ENDIF

    ! Find F(2) at NU(2) = NU(1).....
    IF (LINITS.EQ.2) THEN
       NU(2)=NU(1)+INISTEP
       NUVAL = NU(2)
       GOTO 900
    ENDIF

    IF (LINITS.EQ.3) THEN
       ! set F(2) = VALFUN
       F(2) = VALFUN
       ! Find the third point
       NU(3) = NU(2)+INISTEP
       NUVAL = NU(3)
       GOTO 900 
    ENDIF

    IF (LINITS.EQ.4) THEN
       ! Set F(3) =VALFUN and the initial VALFUN0
       F(3) = VALFUN
       VALFUN0 =MIN(F(1),F(2),F(3))
       ewrite(3,*) 'Finish the choice of the initial three points'
       ewrite(3,*) 'NU1,NU2,NU3',NU(1), NU(2),NU(3)
       ewrite(3,*) 'F1,F2,F3',F(1),F(2),F(3)
       GOTO 10
    ENDIF


    ! Has the last line search converged?
    ewrite(3,*) 'VALFUN0=',VALFUN0    
    IF( (ABS(VALFUN-VALFUN0).LT.LINERR).OR.VALFUN.LT.2.0E-8) THEN 
       LINCONVER=.TRUE.
       ! Re-initilise LINITS and NUVAL for the next line search along the new direction
       LINITS =1
       NUVAL = 0.0
       GOTO 1000  
    ELSE
       ! continue the current line search
       VALFUN0 = VALFUN
    ENDIF
    ewrite(3,*) 'VALFUN=',VALFUN0    

10  CONTINUE
    !....    Obtain new value of NU to explore

    GOTNUV=.FALSE.
    IF(LINITS.GE.5) THEN
       ewrite(3,*) 'NUVAL=, in sub_nlcg',NUVAL

       !....   If VALFUN is better than the other 3 increase golden ratio
       IF (VALFUN .LT. MIN(F(1),F(2),F(3)) ) THEN
          GOLDR2= MIN(2.0*GOLDR2,GOLDR)
          ewrite(3,*) 'GOLDR2=',GOLDR2
       ELSE
          IF (VALFUN .GT. MAX(F(1),F(2),F(3)) )             &
               GOLDR2 = MAX(0.5*GOLDR2,0.5) 
       ENDIF

       !....       WORK OUT What to do with 4th point?
       !....       if interpolating always accept new point replacing one of the 
       !....       points
       IF((NUVAL.GT.NU(1)).AND.(NUVAL.LT.NU(3))) THEN
          IF( (NUVAL.LT.NU(2)) ) THEN
             !....   Replace 3rd point
             NU(3)=NUVAL
             F(3) =VALFUN
          ELSE
             !....   Replace 1ST point
             NU(1)=NUVAL
             F(1) =VALFUN
          ENDIF
       ELSE  
          !....   We are extrapolating...
          IF( VALFUN.LT.MAX(F(1),F(2),F(3)) ) THEN
             !....   Replace furthest point in NU(i) away from NUVAL. 
             IF(NUVAL.LE.NU(1)) THEN
                NU(3)=NUVAL
                F(3)=VALFUN
             ELSE
                NU(1)=NUVAL
                F(1)=VALFUN
             ENDIF
          ELSE
             !....   If new pt is worse than the other 3 then take mid pt of nearest 2 NU's.
             GOTNUV=.TRUE.
             IF(NUVAL.GT.NU(3)) NUVAL=0.5*(NU(2)+NU(3))
             IF(NUVAL.LT.NU(1)) NUVAL=0.5*(NU(1)+NU(2))
          ENDIF
       ENDIF

       !....    endof IF(LINITS.GE.5) THEN...
    ELSE
       GOLDR2=GOLDR
       !....    endof IF(LINITS.GE.5) THEN ELSE...
    ENDIF
    !

    IF(LINITS.GE.4) THEN
       !....   Predict a new minimum NUVAL. Sort in order of increasing NU.
       CALL SORTNU(3,NU,F)
       !     IF(.NOT.GOTNUV) THEN
       !....   Find minimum of quadratic (NUVAL) - with constraints. 
       CALL MINQUA_F(NUVAL,3,NU,F,MINNU,MAXNU,GOLDR2)
       !     ENDIF
    ENDIF
    ewrite(3,*) '--------------------------'
    ewrite(3,*) 'NU', NU
    ewrite(3,*) 'F', F
    ewrite(3,*) '--------------------------'
    ! set the next line search number
900 LINITS = LINITS+1    
1000 CONTINUE
    RETURN
  END SUBROUTINE LINESEARCH
  !  
  !
  !
  !         SUBROUTINE MINQUA(NUVAL,NSOR,NU,F,MINNU,MAXNU,GOLDR,LINCONVER)
  SUBROUTINE MINQUA_F(NUVAL,NSOR,NU,F,MINNU,MAXNU,GOLDR)
    !---------------------------------------------------------------------------------------
    !
    ! This sub works out the new value of NUVAL to sample 
    ! at given the three pairs (NU(1..3),F(1..3)) to maximise F. 
    ! Work out the min of quadratic through the 3 points 
    ! and find minimum subject to constraints of line 
    ! search. 
    ! MINNU,MAXNU are the min and max values that NUVAL can take on. 
    ! Make sure, new pt. NUVAL is not too close to neighbouring pts NU(i)
    ! (controlled by 1/GOLDR).
    ! When extrapolating, make sure we do not extrapolate too far
    ! (controlled by GOLDR).
    ! GOLDR is specified in the setup-file *.inv0
    !
    !---------------------------------------------------------------------------------------


    !.... Declaration of Variables
    IMPLICIT NONE

    REAL GOLDR,GOLINT
    !.... GOLINT is the golden ratio for interpolation 
    !(new point can't be too close to any of old points)
    PARAMETER(GOLINT=5)
    INTEGER NSOR
    REAL NU(NSOR),F(NSOR)
    REAL NUVAL
    REAL MINNU,MAXNU
    !.... Local variables
    REAL Y,X,AA,BB,CC
    REAL A,B,C,D
    REAL A11,A12,A21,A22,DET,DET2
    REAL RMIN,RMAX
    ! The matrix is   A   B
    !                 C   D.
    INTEGER I,J
    REAL EPS 
    PARAMETER(EPS=1.E-30)
    LOGICAL LINCONVER
    !
    !.... 
    ewrite(3,*) '****** just in MINQUA'
    ewrite(3,*) 'NU:',NU
    ewrite(3,*) 'F:',F
    !

    A=NU(2)-NU(1)
    B=NU(2)*NU(2) - NU(1)*NU(1) 
    C=NU(3)-NU(1)
    D=NU(3)*NU(3) - NU(1)*NU(1) 

    ewrite(3,*) 'A,B:',A,B
    ewrite(3,*) 'C,D:',C,D
    !      ewrite(3,*) 'DET=',DET
    !      ewrite(3,*) 'GOING TO SOLVE'

    !.... Inverse...
    DET=A*D-B*C
    !         DET2=1.0
    !         DET2=A*D-B*C
    !         DET =MAX(ABS(A),ABS(B),ABS(C),ABS(D))
    !....
    LINCONVER=.FALSE.
    ewrite(3,*) 'DET=',DET
    !      IF(DET .EQ. 0.) THEN
    IF(DET .LT. EPS) THEN 
       LINCONVER=.TRUE.
       RETURN
    ENDIF
    A11= D/DET
    A12=-B/DET
    A21=-C/DET
    A22= A/DET
    !
    BB=A11*(F(2)-F(1))+A12*(F(3)-F(1))
    CC=A21*(F(2)-F(1))+A22*(F(3)-F(1))
    AA=-(NU(1)*BB + NU(1)*NU(1)*CC) + F(1)
    !
    ewrite(3,*) 'FINISHED SMLIN3 AA,BB,CC:',AA,BB,CC
    ewrite(3,*) 'checking curve:'
    ewrite(3,*) 'from interpolated curve f1,f2,f3:',  &
         AA+BB*NU(1)+CC*NU(1)**2,                    &             
         AA+BB*NU(2)+CC*NU(2)**2,                    &
         AA+BB*NU(3)+CC*NU(3)**2                     
    !
    IF(CC.LE.1.E-7) THEN
       !.... Draw line through points 1 & 3. 
       IF(F(3).LT.F(1)) THEN
          NUVAL=MIN(NU(3)+ GOLDR*(NU(3)-NU(1)),MAXNU) 
       ELSE
          NUVAL=MAX(NU(1)- GOLDR*(NU(3)-NU(1)),MINNU) 
       ENDIF
    ELSE
       NUVAL=-BB/(2.*CC)
       ewrite(3,*) 'before corrections BB,CC,NUVAL=',BB,CC,NUVAL
       NUVAL=MIN(NUVAL,MAXNU) 
       NUVAL=MAX(NUVAL,MINNU) 
    ENDIF
    !
    !.... Make sure NUVAL IS BOUNDED BY GOLDR=golden ratio.
    RMIN=1./GOLINT
    RMAX= GOLDR 
    IF(NUVAL.LE.NU(1)) THEN
       NUVAL=MIN(NUVAL, NU(1)-RMIN*(NU(2)-NU(1)) )
       NUVAL=MAX(NUVAL, NU(1)-RMAX*(NU(2)-NU(1)) )
    ENDIF
    IF ( (NUVAL.GE.NU(1)).AND.(NUVAL.LE.NU(2)) ) THEN
       !.... Make sure that NUVAL is not too close to NU(1) or NU(2)
       NUVAL=MAX(NUVAL,NU(1)+RMIN*(NU(2)-NU(1)))
       NUVAL=MIN(NUVAL,NU(2)-RMIN*(NU(2)-NU(1)))
    ENDIF
    IF (NUVAL.GE.NU(3)) THEN
       NUVAL=MAX(NUVAL, NU(3)+RMIN*(NU(3)-NU(2)) )
       NUVAL=MIN(NUVAL, NU(3)+RMAX*(NU(3)-NU(2)) )
    ENDIF
    IF ( (NUVAL.GE.NU(2)).AND.(NUVAL.LE.NU(3)) ) THEN
       NUVAL=MAX(NUVAL,NU(2)+RMIN*(NU(3)-NU(2)))
       NUVAL=MIN(NUVAL,NU(3)-RMIN*(NU(3)-NU(2)))
    ENDIF
    ewrite(3,*) 'finishing minqua NUVAL=',NUVAL
    !     
    RETURN
  END SUBROUTINE MINQUA_F
  !
  !
end module NONLINCG_module
