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
      SUBROUTINE SECTLP( XSTART, YSTART, ZSTART, DIRX, DIRY, DIRZ,
     :                   X, Y, Z, NPNTS, XNORM, YNORM, ZNORM,
     :                   BOUNDS, XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX,
     :                   XSECT, YSECT, ZSECT, DIST, INSIDE )
C
      IMPLICIT NONE
C=====================================================================
C
C - This subroutine takes a line from XSTART,YSTART,ZSTART in *unit* direction
C - DIRX,DIRY,DIRZ and a plane defined by the *unit* normal XNORM,YNORM,ZNORM
C - and the first point in X,Y,Z.
C - It returns in DIST the distance along the line from XSTART,YSTART,ZSTART
C - to the intersection point, which is returned in XSECT,YSECT,ZSECT.
C
C - If INSIDE is sent as .true. then it checks to see if the intersection
C - point is inside the convex closed NPNTS sided polygon defined by the NPNTS
C - points in X,Y,Z. If it is not inside it returns INSIDE as .false.
C - Note that the polygon can be either clockwise or anticlockwise - the
C - orientation is worked out relative to the plane normal.
C - NOTE THAT INSIDE SHOULD *NEVER* BE SENT AS A CONSTANT!!!
C
C - If BOUNDS is .true. it expects the min/max bounds of X,Y,Z to be provided
C - in XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX, otherwise it calculates them itself and
C - returns them (these are used to speed up checking if the intersect point
C - is inside the polygon). Note that you should only supply bounds to this
C - subroutine that have been previously calculated by this subroutine.
C - If you want to re-use them on subsequent calls, you should set BOUNDS to
C - .false. on the first call, store the min/max bounds returned, and set
C - BOUNDS to .true. on subsequent calls.
C
C - NB: There is a possibility that a point very slightly outside the polygon
C - will be returned as inside. However, the code has been written such that
C - the converse should never be true.
C
C=====================================================================
C
      INTEGER NPNTS
C
      REAL XSTART, YSTART, ZSTART, DIRX, DIRY, DIRZ,
     :     X(NPNTS), Y(NPNTS), Z(NPNTS), DIST,
     :     XEPS, YEPS, ZEPS, XNORM, YNORM, ZNORM, XMAX, YMAX, ZMAX,
     :     XSECT, YSECT, ZSECT, XMIN, YMIN, ZMIN
C
      LOGICAL BOUNDS, INSIDE
C
C - local
C
      INTEGER I, J
C
      REAL XX, YY, ZZ, X2, Y2, Z2, DOT, DIRECT, XDIFF, YDIFF, ZDIFF,
     :     DIFF, DX, DY, DZ, X1, Y1, Z1
C
      DOT = -DIRX*XNORM - DIRY*YNORM - DIRZ*ZNORM
C
      IF( ABS( DOT ) .LT. 1E-10 ) THEN
         INSIDE = .FALSE.
         DIST   = -1E+30
         XSECT = 0.0
         YSECT = 0.0
         ZSECT = 0.0
         GOTO 9000
      END IF
C
      DX = XSTART - X(1)
      DY = YSTART - Y(1)
      DZ = ZSTART - Z(1)
C
      DIST = DX*XNORM + DY*YNORM + DZ*ZNORM
      DIST = DIST/DOT
C
      XSECT = XSTART + DIST*DIRX
      YSECT = YSTART + DIST*DIRY
      ZSECT = ZSTART + DIST*DIRZ
C
      XEPS = MAX( ABS(XSTART), ABS(DIST*DIRX) )
      YEPS = MAX( ABS(YSTART), ABS(DIST*DIRY) )
      ZEPS = MAX( ABS(ZSTART), ABS(DIST*DIRZ) )
C
      IF( .NOT. INSIDE ) GOTO 9000
C
      INSIDE = .FALSE.
C
      IF( NPNTS .LT. 3 ) GOTO 9000
C
      IF( .NOT. BOUNDS ) THEN
C
         XMIN = X(1)
         XMAX = X(1)
         YMIN = Y(1)
         YMAX = Y(1)
         ZMIN = Z(1)
         ZMAX = Z(1)
C
         DO I = 2, NPNTS
            XMIN = MIN( XMIN, X(I) )
            XMAX = MAX( XMAX, X(I) )
            YMIN = MIN( YMIN, Y(I) )
            YMAX = MAX( YMAX, Y(I) )
            ZMIN = MIN( ZMIN, Z(I) )
            ZMAX = MAX( ZMAX, Z(I) )
         END DO
C
         XEPS = MAX( XEPS, ABS(XMAX+XMIN) )
         YEPS = MAX( YEPS, ABS(YMAX+YMIN) )
         ZEPS = MAX( ZEPS, ABS(ZMAX+ZMIN) )
C
         XDIFF = MAX( (XMAX-XMIN)*1E-4, XEPS*1E-7 )
         YDIFF = MAX( (YMAX-YMIN)*1E-4, YEPS*1E-7 )
         ZDIFF = MAX( (ZMAX-ZMIN)*1E-4, ZEPS*1E-7 )
C
         DIFF = MAX(MAX(XDIFF,YDIFF),ZDIFF)
C
         XDIFF = MAX( XDIFF, DIFF*1E-4 )
         YDIFF = MAX( YDIFF, DIFF*1E-4 )
         ZDIFF = MAX( ZDIFF, DIFF*1E-4 )
C
         XMIN = XMIN - XDIFF
         XMAX = XMAX + XDIFF
         YMIN = YMIN - YDIFF
         YMAX = YMAX + YDIFF
         ZMIN = ZMIN - ZDIFF
         ZMAX = ZMAX + ZDIFF
C
      ELSE
C
         XDIFF = (XMAX - XMIN)*1E-4
         YDIFF = (YMAX - YMIN)*1E-4
         ZDIFF = (ZMAX - ZMIN)*1E-4
C
         DIFF = MAX(MAX(XDIFF,YDIFF),ZDIFF)
C
C         XDIFF = MAX( XDIFF, DIFF*1E-4 )
C         YDIFF = MAX( YDIFF, DIFF*1E-4 )
C         ZDIFF = MAX( ZDIFF, DIFF*1E-4 )
C
      END IF
C
      IF( XSECT .LT. XMIN ) GOTO 9000
      IF( XSECT .GT. XMAX ) GOTO 9000
      IF( YSECT .LT. YMIN ) GOTO 9000
      IF( YSECT .GT. YMAX ) GOTO 9000
      IF( ZSECT .LT. ZMIN ) GOTO 9000
      IF( ZSECT .GT. ZMAX ) GOTO 9000
C
c      if( zsect.lt.1e-5 ) then
c         print*,'Node vectors are:'
c         do ii = 1, npnts
c            print*,x(ii)-xsect,y(ii)-ysect,z(ii)-zsect
c         end do
c      end if
c
      X1 = X(1) - X(2)
      X2 = X(3) - X(2)
      Y1 = Y(1) - Y(2)
      Y2 = Y(3) - Y(2)
      Z1 = Z(1) - Z(2)
      Z2 = Z(3) - Z(2)
C
      XX = Y1*Z2 - Z1*Y2
      YY = Z1*X2 - X1*Z2
      ZZ = X1*Y2 - Y1*X2
C
      DIRECT = XX*XNORM + YY*YNORM + ZZ*ZNORM
C
      IF( DIRECT .NE. 0.0 ) THEN
         DIRECT = DIRECT/ABS(DIRECT)
      ELSE
         DIRECT = 1.0
      END IF
C
      X1 = X(NPNTS) - XSECT
      Y1 = Y(NPNTS) - YSECT
      Z1 = Z(NPNTS) - ZSECT
C
      DO J = 1, NPNTS
C
         X2 = X(J) - XSECT
         Y2 = Y(J) - YSECT
         Z2 = Z(J) - ZSECT
C
         XX = Y1*Z2 - Y2*Z1
         YY = Z1*X2 - Z2*X1
         ZZ = X1*Y2 - X2*Y1
C
         DOT = XX*XNORM + YY*YNORM + ZZ*ZNORM
c
c         if( zsect.lt.1e-5 ) then
c            print*,'nodes ',j-1,j,dot*direct
c            print*,'xx,yy,zz:',xx,yy,zz
c         end if
C
         IF( DOT*DIRECT .GT. DIFF ) GOTO 9000
C
         X1 = X2
         Y1 = Y2
         Z1 = Z2
C
      END DO
C
      INSIDE = .TRUE.
C
 9000 RETURN
      END
C
