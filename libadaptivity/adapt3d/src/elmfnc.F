#include "fdebug_adapt.h"

!#define LIPNIKOV_FUNCTIONAL
#ifdef LIPNIKOV_FUNCTIONAL
      function elmfnc(biglst, nodlst, ielm, nod1, nod2, nod3, nod4,
     &         undlen) result(functional)

      use addebug
      implicit none
      include 'blknew.i'
      include 'blkbig.i'

      integer :: ielm, nod1, nod2, nod3, nod4
      real :: undlen
      integer :: cache 
      integer :: ielm2 
      real :: functional
      real, parameter :: m=2.0*huge(0)-3.0, c= -1.0*huge(0)+2.0 ! mapping between reals and ints for the cache
      real, dimension(9) :: metric
      real, dimension(3) :: tmppos
      real, dimension(3, 3) :: positions
      real :: volume ! volume with respect to the metric
      real :: metric_determinant ! determinant of the metric
      real :: tetvol ! volume of the tet in physical space
      real :: edgelens ! 2-norm of vector of edge lengths in metric space
      real :: F ! from Agouzal et al, 1999
      real :: x
      real, parameter :: factor = 6.0**(2.5) * sqrt(2.0)
      logical :: negate
      integer, dimension(4) :: nods
      integer :: k
      real, parameter :: bst = 12.0 / sqrt(2.0) ! inverse of ideal tet volume
      integer :: intval ! undlen cast to an integer


C - Sort the nodes into ascending order
C
      NODS(1) = NOD1
      NODS(2) = NOD2
      NODS(3) = NOD3
      NODS(4) = NOD4
      NEGATE = .false.
C
      K = NODS(1)
      IF( K .GT. NODS(2) ) THEN
        NODS(1) = NODS(2)
        NODS(2) = K
        NEGATE = .NOT. NEGATE
      ELSE
        K = NODS(2)
      END IF
      IF( K .GT. NODS(3) ) THEN
        NODS(2) = NODS(3)
        NODS(3) = K
        NEGATE = .NOT. NEGATE
      ELSE
        K = NODS(3)
      END IF
      IF( K .GT. NODS(4) ) THEN
        NODS(3) = NODS(4)
        NODS(4) = K
        NEGATE = .NOT. NEGATE
      END IF
      K = NODS(1)
      IF( K .GT. NODS(2) ) THEN
        NODS(1) = NODS(2)
        NODS(2) = K
        NEGATE = .NOT. NEGATE
      ELSE
        K = NODS(2)
      END IF
      IF( K .GT. NODS(3) ) THEN
        NODS(2) = NODS(3)
        NODS(3) = K
        NEGATE = .NOT. NEGATE
      END IF
      K = NODS(1)
      IF( K .GT. NODS(2) ) THEN
        NODS(1) = NODS(2)
        NODS(2) = K
        NEGATE = .NOT. NEGATE
      END IF

      if (ielm > 0) then
        cache = biglst(7, ielm)
        ielm2 = biglst(nxtbig, ielm)
      else
        cache = 0
        ielm2 = -1
      end if

      ! if we have no cached information, then we need to compute it, na ja?
      if (cache == 0) then
        ! first things first, let us fetch the metric and positions
        tmppos = nodlst(1:3, nods(1))
        metric = nodlst(7:15, nods(1))

        positions(:, 1) = nodlst(1:3, nods(2)) - tmppos
        metric = metric + nodlst(7:15, nods(2))

        positions(:, 2) = nodlst(1:3, nods(3)) - tmppos
        metric = metric + nodlst(7:15, nods(3))

        positions(:, 3) = nodlst(1:3, nods(4)) - tmppos
        metric = metric + nodlst(7:15, nods(4))

        metric = metric / 4.0

        metric_determinant = 
     &  metric(1) * (metric(5) * metric(9) - metric(6) * metric(8))
     &- metric(2) * (metric(4) * metric(9) - metric(6) * metric(7))
     &+ metric(3) * (metric(4) * metric(8) - metric(5) * metric(7))

        tetvol = 
     &  positions(1,1) * (positions(2,2) * positions(3,3) - positions(2,3) * positions(3,2))
     &- positions(1,2) * (positions(2,1) * positions(3,3) - positions(2,3) * positions(3,1))
     &+ positions(1,3) * (positions(2,1) * positions(3,2) - positions(2,2) * positions(3,1))
        tetvol = tetvol / 6.0

        undlen = sqrt(metric_determinant) * tetvol
!        if(negate) undlen = - undlen
        volume = abs(undlen)
        undlen = undlen * bst

        !write(0,*) "computed undlen: ", undlen
        tmppos = positions(:, 1)
        edgelens =
     &  tmppos(1) * (metric(1) * tmppos(1) + metric(2) * tmppos(2) + metric(3) * tmppos(3))
     &+ tmppos(2) * (metric(4) * tmppos(1) + metric(5) * tmppos(2) + metric(6) * tmppos(3))
     &+ tmppos(3) * (metric(7) * tmppos(1) + metric(8) * tmppos(2) + metric(9) * tmppos(3))

        tmppos = positions(:, 2)
        edgelens = edgelens +
     &  tmppos(1) * (metric(1) * tmppos(1) + metric(2) * tmppos(2) + metric(3) * tmppos(3))
     &+ tmppos(2) * (metric(4) * tmppos(1) + metric(5) * tmppos(2) + metric(6) * tmppos(3))
     &+ tmppos(3) * (metric(7) * tmppos(1) + metric(8) * tmppos(2) + metric(9) * tmppos(3))

        tmppos = positions(:, 3)
        edgelens = edgelens +
     &  tmppos(1) * (metric(1) * tmppos(1) + metric(2) * tmppos(2) + metric(3) * tmppos(3))
     &+ tmppos(2) * (metric(4) * tmppos(1) + metric(5) * tmppos(2) + metric(6) * tmppos(3))
     &+ tmppos(3) * (metric(7) * tmppos(1) + metric(8) * tmppos(2) + metric(9) * tmppos(3))

        tmppos = positions(:, 1) - positions(:, 2)
        edgelens = edgelens +
     &  tmppos(1) * (metric(1) * tmppos(1) + metric(2) * tmppos(2) + metric(3) * tmppos(3))
     &+ tmppos(2) * (metric(4) * tmppos(1) + metric(5) * tmppos(2) + metric(6) * tmppos(3))
     &+ tmppos(3) * (metric(7) * tmppos(1) + metric(8) * tmppos(2) + metric(9) * tmppos(3))

        tmppos = positions(:, 1) - positions(:, 3)
        edgelens = edgelens +
     &  tmppos(1) * (metric(1) * tmppos(1) + metric(2) * tmppos(2) + metric(3) * tmppos(3))
     &+ tmppos(2) * (metric(4) * tmppos(1) + metric(5) * tmppos(2) + metric(6) * tmppos(3))
     &+ tmppos(3) * (metric(7) * tmppos(1) + metric(8) * tmppos(2) + metric(9) * tmppos(3))

        tmppos = positions(:, 2) - positions(:, 3)
        edgelens = edgelens +
     &  tmppos(1) * (metric(1) * tmppos(1) + metric(2) * tmppos(2) + metric(3) * tmppos(3))
     &+ tmppos(2) * (metric(4) * tmppos(1) + metric(5) * tmppos(2) + metric(6) * tmppos(3))
     &+ tmppos(3) * (metric(7) * tmppos(1) + metric(8) * tmppos(2) + metric(9) * tmppos(3))

        edgelens = sqrt(edgelens)

        x = edgelens / sqrt(6.0); x = min(x, 1.0/x)
        F = (x * (2-x))**3
        functional = factor * (volume / edgelens**3) * F

        undlen = min(max(undlen, -2.01E+3), 2.01E+3)
        intval = int(undlen * 1000000)
        undlen = float(intval) / 1000000
        if (ielm > 0) then
          biglst(7, ielm)  = intval
          biglst(7, ielm2) = int(m * functional + c)
        end if
      else
        undlen = float(cache) / 1000000
        functional = (biglst(7, ielm2) - c) / m
      end if

      if (negate) undlen = -undlen
      !write(0,*) "negate: ", negate, "; returning undlen: ", undlen
      assert(functional >= 0.0)
      assert(functional <= 1.0)

      ! And swap for adaptivity's convention
      functional = 5.0/max(functional, 5.0e-30) - 5.0

      end function elmfnc
#else
C Copyright (C) 2006 Imperial College London and others.
C 
C Please see the AUTHORS file in the main source directory for a full list
C of copyright holders.
C 
C Adrian Umpleby
C Applied Modelling and Computation Group
C Department of Earth Science and Engineering
C Imperial College London
C 
C adrian@Imperial.ac.uk
C 
C This library is free software; you can redistribute it and/or
C modify it under the terms of the GNU Lesser General Public
C License as published by the Free Software Foundation; either
C version 2.1 of the License.
C 
C This library is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C Lesser General Public License for more details.
C 
C You should have received a copy of the GNU Lesser General Public
C License along with this library; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
C USA
#include "ewrite.h"
      REAL FUNCTION ELMFNC( BIGLST, NODLST, IELM,
     :                      NOD1, NOD2, NOD3, NOD4, UNDLEN )
C-----------------------------------------------------------------------
C
C - This subroutine calculates the element functional.
C - The functional is calculated from the four node pointers (NOD1-NOD4)
C - which have to be first arranged in ascending order (this makes sure
C - that the element will have the same functional every time it is
C - calculated, by removing round-off errors because the order of nodes
C - is the same each time).
C - If IELM>0 and element IELM has not already had its functional worked
C - out then xeveral critical values (sum of edge functionals etc) are
C - worked out and stored for use next time the functional is required
C - for this element. This reduces the amount of calculation required.
C
C-----------------------------------------------------------------------
      use write_log
      IMPLICIT NONE
C
      REAL, PARAMETER :: BST = 4.898979486  ! sqrt(24)
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkerr.i'
C
      REAL ANGLEFACT
      COMMON / BLKTMP /  ANGLEFACT
C
      INTEGER IELM, NOD1, NOD2, NOD3, NOD4
C
      REAL UNDLEN
C
      integer i, j, k, nods(4), intval, nxt1, s
C
      REAL X(4), Y(4), Z(4), V1(3), V2(3), V3(3), D1, D2, D3, D, Q,
     :     M(9), AREAS(4), VOL, EDGFNC, a, b, c,xx(4),yy(4),zz(4),
     :     L(6), det, getdtm, angfnc
c
      logical negate, tstdts
c
      tstdts = .false.
c
      negate = .false.
C
      ELMFNC = 1E+30
      UNDLEN = 0.0
C
      IF( IERR .NE. 0 ) RETURN
C
C - Sort the nodes into ascending order
C
      NODS(1) = NOD1
      NODS(2) = NOD2
      NODS(3) = NOD3
      NODS(4) = NOD4
C
      K = NODS(1)
      IF( K .GT. NODS(2) ) THEN
        NODS(1) = NODS(2)
        NODS(2) = K
        NEGATE = .NOT. NEGATE
      ELSE
        K = NODS(2)
      END IF
      IF( K .GT. NODS(3) ) THEN
        NODS(2) = NODS(3)
        NODS(3) = K
        NEGATE = .NOT. NEGATE
      ELSE
        K = NODS(3)
      END IF
      IF( K .GT. NODS(4) ) THEN
        NODS(3) = NODS(4)
        NODS(4) = K
        NEGATE = .NOT. NEGATE
      END IF
      K = NODS(1)
      IF( K .GT. NODS(2) ) THEN
        NODS(1) = NODS(2)
        NODS(2) = K
        NEGATE = .NOT. NEGATE
      ELSE
        K = NODS(2)
      END IF
      IF( K .GT. NODS(3) ) THEN
        NODS(2) = NODS(3)
        NODS(3) = K
        NEGATE = .NOT. NEGATE
      END IF
      K = NODS(1)
      IF( K .GT. NODS(2) ) THEN
        NODS(1) = NODS(2)
        NODS(2) = K
        NEGATE = .NOT. NEGATE
      END IF
C
cc      if( debug ) print*,'   sorted: ',nods
C
C - Read the in-sphere radius stored in BIGLST
C
      INTVAL = 0
      NXT1   = 0
      IF( IELM .GT. 0 ) THEN
         INTVAL = BIGLST(7,IELM)
         NXT1   = BIGLST(NXTBIG,IELM)
      END IF
C
C - If there is no stored in-sphere info, then calculate it
C
      IF( INTVAL .EQ. 0 ) THEN
C
         DO J = 1, 9
            M(J) = 0.0
         END DO
C
         DO I = 1, 4
            if( tstdts .or. debug ) then
               det = getdtm( nodlst(7,nods(i)) )
               if( det .lt. 1e-20 ) then
                  print*,'+++ ELMFNC: Warning, small det:', det,nods(i)
                  print*,(nodlst(j,nods(i)),j=7,9)
                  print*,(nodlst(j,nods(i)),j=10,12)
                  print*,(nodlst(j,nods(i)),j=13,15)
               end if
            end if
            X(I) = NODLST(1,NODS(I))
            Y(I) = NODLST(2,NODS(I))
            Z(I) = NODLST(3,NODS(I))
            DO J = 1, 9
               M(J) = M(J) + NODLST(J+6,NODS(I))*0.25
            END DO
         END DO
C
C - get all the element info with respect to metric M
C - (Volume, Face areas, Edge lengths, In-sphere, Quality)
C
         CALL MTETIN( X, Y, Z, M, VOL, AREAS, L, UNDLEN, Q )
C
cc         if( debug ) print*,'   vol,rad,qly: ',vol,undlen,q
C
C - If we are using the quality of the element (instead of the in-sphere)...
C
         if( useqly ) undlen = 1.0/q/bst
C
         D = UNDLEN*BST
C
C - Don't worry about the maximum in-sphere radius, since such large
C - sizes will not contribute anything significant to the functional
C
            D = MIN(MAX(D,-2.01E+3),2.01E+3)
C
            INTVAL = INT(D*1000000)
            IF( ABS(INTVAL) .LE. 1 ) THEN
              IF( INTVAL .EQ. 0 .AND. IELM .LE. 0 ) THEN
                ELMFNC = 1E+32
                IF( NEGATE ) UNDLEN = -UNDLEN
                RETURN
              END IF
              IF( D .GT. 0.0 ) THEN
                INTVAL = 1
              ELSE IF( D .LT. 0.0 ) THEN
                INTVAL = -1
              ELSE
                write(0,*) 
     :                '*** ELMFNC: GOT ZERO IN-SPHERE!',IELM
                IF( IELM .LE. 0 ) write(0,*)
     :                '  Nodes: ',NOD1,NOD2,NOD3,NOD4
              END IF
              IF( IELM .GT. 0 .AND. INTVAL .EQ. 0 ) THEN
                write(0,*)
     :               '+++ Existing elem with zero in-sphere:',ielm
                write(0,*)
     :               '  d,v,r,q:',d,vol,undlen,q
                write(0,*)
     :               '  node 1: ',(nodlst(j,nod1),j=1,3)
                write(0,*)
     :               '  node 2: ',(nodlst(j,nod2),j=1,3)
                write(0,*)
     :               '  node 3: ',(nodlst(j,nod3),j=1,3)
                write(0,*)
     :               '  node 4: ',(nodlst(j,nod4),j=1,3)
                write(0,*)
     :               '  Element metric:'
                write(0,*)(m(i),i=1,3)
                write(0,*)(m(i),i=4,6)
                write(0,*)(m(i),i=7,9)
                do i = 1, 9
                  m(i) = 0.0
                end do
                m(1) = 1.0
                m(5) = 1.0
                m(9) = 1.0
                write(0,*)
     :               '  Calculating values in flat space...'
                call mtetin( x, y, z, m, vol, areas, l, undlen, q )
                write(0,*)
     :               '    r,v,q: ',undlen,vol,q
                IERR = -676
                RETURN
              END IF
            END IF
            IF( IELM .GT. 0 ) BIGLST(7,IELM) = INTVAL
            D = FLOAT(INTVAL)/1000000
C
      ELSE
C
         D = FLOAT(INTVAL)/1000000
C
      END IF
C
C - work out the *signed* in-sphere radius (or quality) to be returned
C - (this is based on the node numbering originally sent into the call
C - and is relative to the value calculated when the nodes are in ascending
C - order in the calculation).
C
      UNDLEN = D/BST
      IF( NEGATE ) UNDLEN = -UNDLEN
C
C - Calculate the in-sphere contribution to the element functional.
C - This is 'fiddled' so it gets large extremely quickly once the
C - in-sphere becomes very small. Note that the value is continuous,
C - even though different equations are used for different ranges of
C - the in-sphere radius.
C
      D = ABS(D)
      IF( D .GT. 9.9E-6 ) THEN
         D = 1.0/D - 1.0
         D = D*D
      ELSE IF( D .GT. 9.9E-7 ) THEN
         D = D*23713.74
         D = 1.0/D
         D = D*D*D*D
         D = D*D*D*D
      ELSE
         D = 4E+26
      END IF
C
      EDGON = .NOT. EDGON
C
      IF( EDGON ) THEN
C
         INTVAL = 0
C
C - retreive the edge functional info (if NXT1>0)
C
         IF( NXT1 .GT. 0 ) INTVAL = BIGLST(7,NXT1)
C
         IF( INTVAL .EQ. 0 ) THEN
C
C - calculate (and store, if NXT1>0) half the sum of edge functionals
C
c           L(1) = 1.0 - L(1)
c           L(2) = 1.0 - L(2)
c           L(3) = 1.0 - L(3)
c           L(4) = 1.0 - L(4)
c           L(5) = 1.0 - L(5)
c           L(6) = 1.0 - L(6)
c           D1   = L(1)*L(1) + L(2)*L(2) + L(3)*L(3)
c     :          + L(4)*L(4) + L(5)*L(5) + L(6)*L(6)
c
           D1 = 0.5 * ( EDGFNC( BIGLST, NODLST, NOD1, NOD2, D2 )
     :                + EDGFNC( BIGLST, NODLST, NOD1, NOD3, D3 )
     :                + EDGFNC( BIGLST, NODLST, NOD1, NOD4, D2 )
     :                + EDGFNC( BIGLST, NODLST, NOD2, NOD3, D3 )
     :                + EDGFNC( BIGLST, NODLST, NOD2, NOD4, D2 )
     :                + EDGFNC( BIGLST, NODLST, NOD3, NOD4, D3 ) )
C
           IF( ANGLEFACT .GT. 0.0 )
     :         D1 = D1 + ANGFNC( BIGLST, NODLST, NODS )*ANGLEFACT
C
C - large edges will be close to the functional bound, so we can cap them
C
           D1 = MIN(D1,2.01E+6)
           INTVAL = INT(D1*1000)
           IF( INTVAL .EQ. 0 ) INTVAL = 1
           IF( NXT1 .GT. 0 ) BIGLST(7,NXT1) = INTVAL
           D1 = FLOAT(INTVAL)/1000
C
         ELSE
           D1 = FLOAT(INTVAL)/1000
         END IF
C
      ELSE
         D1 = 0.0
      END IF
C
      EDGON = .NOT. EDGON
C
C - This is an attempt to fix a problem where the in-sphere radius
C - contributes too much to the functional for elements which have
C - some kind of constraint to their edges (probably through the
C - geometry), meaning the edges end up becoming far too large (as it
C - tries to increase the in-sphere radius by any means possible).
C - The basic strategy to combat this is to increase the contribution
C - from the edges as they become too large.
C - It's worth remembering that having long edges is actually far worse
C - for interpolation errors than having small/thin elements!
C - Note that the value of the functional contribution is continuous,
C - even though different equations are used for ranges of d1...
C
      if( d1 .gt. 2E+6 ) then
         d1 = 1.6E+25
c         d1 = 3.2E+26
c      else if( d1 .gt. 10 ) then
c         d1 = d1 + d1*d1*d1*d1*d1 - 91102.0
      else if( d1 .gt. 6 ) then
         d1 = d1 + d1*d1*d1*d1 - 1102.0
      else if( d1 .gt. 3 ) then
         d1 = d1 + d1*d1*d1 - 22.0
      else if( d1 .gt. 2 ) then
         d1 = d1 + d1*d1 - 4.0
      end if
C
C - Work out the entire element functional from half the sum of the edge
C - functionals plus the in-sphere functional.
C
      ELMFNC = D1 + D
C
      RETURN
      END
C
C
C
      REAL FUNCTION ANGFNC( BIGLST, NODLST, NODS )
C
      IMPLICIT NONE
C
      INTEGER NODS(4)
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      REAL NX1, NY1, NZ1, NX2, NY2, NZ2, NX3, NY3, NZ3, NX4, NY4, NZ4,
     :     X(6), Y(6), Z(6), V, D
C
      INTEGER I, J, K
C
      angfnc = 1E+30
C
      NX1 = NODLST(1,NODS(1))
      NY1 = NODLST(2,NODS(1))
      NZ1 = NODLST(3,NODS(1))
C
      DO I = 1, 3
         X(I) = NODLST(1,NODS(I+1)) - NX1
         Y(I) = NODLST(2,NODS(I+1)) - NY1
         Z(I) = NODLST(3,NODS(I+1)) - NZ1
      END DO
C
      NX1 = NODLST(1,NODS(2))
      NY1 = NODLST(2,NODS(2))
      NZ1 = NODLST(3,NODS(2))
C
      DO I = 4, 5
         X(I) = NODLST(1,NODS(I-1)) - NX1
         Y(I) = NODLST(2,NODS(I-1)) - NY1
         Z(I) = NODLST(3,NODS(I-1)) - NZ1
      END DO
C
      X(6) = NODLST(1,NODS(4)) - NODLST(1,NODS(3))
      Y(6) = NODLST(2,NODS(4)) - NODLST(2,NODS(3))
      Z(6) = NODLST(3,NODS(4)) - NODLST(3,NODS(3))
C
C - work out normals to faces
C
      NX1 = Y(2)*Z(4) - Y(4)*Z(2)
      NY1 = Z(2)*X(4) - Z(4)*X(2)
      NZ1 = X(2)*Y(4) - X(4)*Y(2)
C
      V = NX1*NX1 + NY1*NY1 + NZ1*NZ1
      if( v .lt. 1e-23 ) then
         return
         print*,'got small normal 1: ',v
         print*,nx1,ny1,nz1
         print*,nods
      else
         V = 1.0/SQRT(V)
         NX1 = NX1*V
         NY1 = NY1*V
         NZ1 = NZ1*V
      end if
C
      NX2 = Y(3)*Z(1) - Y(1)*Z(3)
      NY2 = Z(3)*X(1) - Z(1)*X(3)
      NZ2 = X(3)*Y(1) - X(1)*Y(3)
C
      V = NX2*NX2 + NY2*NY2 + NZ2*NZ2
      if( v .lt. 1e-23 ) then
         return
         print*,'got small normal 2: ',v
         print*,nx2,ny2,nz2
         print*,nods
      else
         V = 1.0/SQRT(V)
         NX2 = NX2*V
         NY2 = NY2*V
         NZ2 = NZ2*V
      end if
C
      NX3 = Y(5)*Z(4) - Y(4)*Z(5)
      NY3 = Z(5)*X(4) - Z(4)*X(5)
      NZ3 = X(5)*Y(4) - X(4)*Y(5)
C
      V = NX3*NX3 + NY3*NY3 + NZ3*NZ3
      if( v .lt. 1e-23 ) then
         return
         print*,'got small normal 3: ',v
         print*,nx3,ny3,nz3
         print*,nods
      else
         V = 1.0/SQRT(V)
         NX3 = NX3*V
         NY3 = NY3*V
         NZ3 = NZ3*V
      end if
C
      NX4 = Y(3)*Z(6) - Y(6)*Z(3)
      NY4 = Z(3)*X(6) - Z(6)*X(3)
      NZ4 = X(3)*Y(6) - X(6)*Y(3)
C
      V = NX4*NX4 + NY4*NY4 + NZ4*NZ4
      if( v .lt. 1e-23 ) then
         return
         print*,'got small normal 4: ',v
         print*,nx4,ny4,nz4
         print*,nods
      else
         V = 1.0/SQRT(V)
         NX4 = NX4*V
         NY4 = NY4*V
         NZ4 = NZ4*V
      end if
C
C - work out the functionals for the angle between each pair of faces
C
      D = NX1*NX2 + NY1*NY2 + NZ1*NZ2 + 1.0
      IF( D .GT. 1E-15 ) THEN
         D = ( 1.0/D - 1.0 )
      ELSE
         return
         print*,'got d = ',d,' for normals 1+2...'
         print*,nx1,ny1,nz1
         print*,nx2,ny2,nz2
         print*,x
         print*,y
         print*,z
         D = 1E+15
      END IF
C
      V = D*D
C
      D = NX1*NX3 + NY1*NY3 + NZ1*NZ3 + 1.0
      IF( D .GT. 1E-15 ) THEN
         D = ( 1.0/D - 1.0 )
      ELSE
         return
         print*,'got d = ',d,' for normals 1+3...'
         print*,nx1,ny1,nz1
         print*,nx3,ny3,nz3
         print*,x
         print*,y
         print*,z
         D = 1E+15
      END IF
C
      V = V + D*D
C
      D = NX1*NX4 + NY1*NY4 + NZ1*NZ4 + 1.0
      IF( D .GT. 1E-15 ) THEN
         D = ( 1.0/D - 1.0 )
      ELSE
         return
         print*,'got d = ',d,' for normals 1+4...'
         print*,nx1,ny1,nz1
         print*,nx4,ny4,nz4
         print*,x
         print*,y
         print*,z
         D = 1E+15
      END IF
C
      V = V + D*D
C
      D = NX2*NX3 + NY2*NY3 + NZ2*NZ3 + 1.0
      IF( D .GT. 1E-15 ) THEN
         D = ( 1.0/D - 1.0 )
      ELSE
         return
         print*,'got d = ',d,' for normals 2+3...'
         print*,nx2,ny2,nz2
         print*,nx3,ny3,nz3
         print*,x
         print*,y
         print*,z
         D = 1E+15
      END IF
C
      V = V + D*D
C
      D = NX2*NX4 + NY2*NY4 + NZ2*NZ4 + 1.0
      IF( D .GT. 1E-15 ) THEN
         D = ( 1.0/D - 1.0 )
      ELSE
         return
         print*,'got d = ',d,' for normals 2+4...'
         print*,nx2,ny2,nz2
         print*,nx4,ny4,nz4
         print*,x
         print*,y
         print*,z
         D = 1E+15
      END IF
C
      V = V + D*D
C
      D = NX3*NX4 + NY3*NY4 + NZ3*NZ4 + 1.0
      IF( D .GT. 1E-15 ) THEN
         D = ( 1.0/D - 1.0 )
      ELSE
         return
         print*,'got d = ',d,' for normals 3+4...'
         print*,nx3,ny3,nz3
         print*,nx4,ny4,nz4
         print*,x
         print*,y
         print*,z
         D = 1E+15
      END IF
C
      ANGFNC = (V + D*D)/6
C
      RETURN
      END
C
#endif
