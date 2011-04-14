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
      REAL FUNCTION EDGDST( BIGLST, NODLST, IEDGE, X, Y, Z,
     :                      ONX, ONY, ONZ,
     :                      DIRX, DIRY, DIRZ, LNGTH,
     :                      MNDST, MNALLW )
C-----------------------------------------------------------------------
C
C - This function checks to see if an edge (IEDGE) punctures a face (defined
C - by three points in X,Y,Z (and with normal ONX,ONY,ONZ).
C - It also finds the distance between IEDGE and the edge defined by the first
C - two nodes in X,Y,Z.
C
C - This is calculated in different ways depending on the situation:
C -     i) the distance between the edges if they are parallel
C -    ii) the distance between the edges if both points of closest approach
C -        lie within the edge
C -   iii) the distance between an edge and the closest node of the other edge
C -        if the points of closest approach are within one edge and not the other
C -    iv) the distance between the closest nodes if the points of closest
C -        approach are within neither edge
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER IEDGE
C
      REAL MNDST, ONX, ONY, ONZ, DIRX, DIRY, DIRZ, LNGTH,
     :     X(3), Y(3), Z(3), MNALLW
C
      INTEGER NODS(2), NXT, LST, IFLAG, INEW, NOD1, NOD2
C
      REAL DIRX2, DIRY2, DIRZ2, LNGTH2, XX, YY, ZZ, DST, NX, NY, NZ,
     :     XSECT, YSECT, ZSECT, XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX,
     :     POS1, POS2, DOT, DD, VX, VY, VZ
C
      LOGICAL INSIDE
C
c      print*,' '
c      print*,'Into EDGDST... ',ierr
c      print*,'   X: ',x
c      print*,'   Y: ',y
c      print*,'   Z: ',z
c      print*,'   normal: ',onx,ony,onz
c      print*,'   dir: ',dirx,diry,dirz
c      print*,'   lngth,mndst: ',lngth,mndst
c
      EDGDST = 0.0
C
      CALL EDGINF( BIGLST, NODLST, IEDGE, NODS, INEW, IFLAG, NXT, LST )
C
      IF( IERR .NE. 0 ) RETURN
c
c      print*,'   edge nodes: ',nods(1),nods(2)
C
      XX = NODLST(1,NODS(1))
      YY = NODLST(2,NODS(1))
      ZZ = NODLST(3,NODS(1))
C
      DIRX2 = NODLST(1,NODS(2)) - XX
      DIRY2 = NODLST(2,NODS(2)) - YY
      DIRZ2 = NODLST(3,NODS(2)) - ZZ
C
      LNGTH2 = DIRX2*DIRX2 + DIRY2*DIRY2 + DIRZ2*DIRZ2
C
      LNGTH2 = SQRT(LNGTH2)
C
      DIRX2 = DIRX2/LNGTH2
      DIRY2 = DIRY2/LNGTH2
      DIRZ2 = DIRZ2/LNGTH2
c
c      print*,'     dir2: ',dirx2,diry2,dirz2
c      print*,'   lngth2: ',lngth2
C
      DOT = DIRX*DIRX2 + DIRY*DIRY2 + DIRZ*DIRZ2
C
      IF( ABS(DOT) .GT. 0.99999 ) THEN
C
C - parallel edges
C
         NX = XX - X(2)
         NY = YY - Y(2)
         NZ = ZZ - Z(2)
C
         DOT = NX*DIRX + NY*DIRY + NZ*DIRZ
         DOT = NX*NX + NY*NY + NZ*NZ - DOT*DOT
C
         EDGDST = SQRT(DOT)
C
c         print*,'Leaving EDGDST...'
         dd = sqrt(dot)
         if( dd .lt. mnallw ) then
            print*,'     +++ Parallel edge separation:',dd
            print*,'         edge nodes: ',nods(1),nods(2)
         end if
         RETURN
C
      END IF
C
      DOT = DIRX2*ONX + DIRY2*ONY + DIRZ2*ONZ
C
      IF( ABS(DOT) .LT. 1E-7 ) THEN
C
C - edge lies within plane of face
C
         INSIDE  = .FALSE.
C
      ELSE
C
         INSIDE = .TRUE.
C
C - check intersection of edge and face
C
c         print*,'   calling SECTLP...'
         CALL SECTLP( XX, YY, ZZ, DIRX2, DIRY2, DIRZ2,
     :                X, Y, Z, 3, ONX, ONY, ONZ,
     :                .FALSE., XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX,
     :                XSECT, YSECT, ZSECT, DST, INSIDE )
c         print*,'   inside,dst: ',inside,dst
c         print*,'  sect: ',xsect,ysect,zsect
C
         IF( INSIDE ) THEN
            IF( DST .GE. -1E-5 .AND. DST .LE. LNGTH2+1E-5 ) THEN
                print*,'     +++ Got intersection with edge ',nods,dst
                EDGDST = 0.0
                RETURN
            END IF
         END IF
C
      END IF
C
C - find minimum distance between infinitely extended edges
C
      XX = X(2) - XX
      YY = Y(2) - YY
      ZZ = Z(2) - ZZ
C
      NX = DIRY*DIRZ2 - DIRZ*DIRY2
      NY = DIRZ*DIRX2 - DIRX*DIRZ2
      NZ = DIRX*DIRY2 - DIRY*DIRX2
C
      DOT = NX*NX + NY*NY + NZ*NZ
C
      DOT = 1.0/SQRT(DOT)
C
      DST = ABS(NX*XX + NY*YY + NZ*ZZ)*DOT
C
C - find distances along extended edges for points of closest approach
C
      VX = DIRZ2*NY-NZ*DIRY2
      VY = DIRX2*NZ-NX*DIRZ2
      VZ = DIRY2*NX-NY*DIRX2
C
      DOT = -VX*DIRX - VY*DIRY - VZ*DIRZ
c      print*,'   got a1: ',vx,vy,vz,dot
C
      POS1 = ( XX*VX + YY*VY + ZZ*VZ )/DOT
C
      VX = DIRZ*NY-NZ*DIRY
      VY = DIRX*NZ-NX*DIRZ
      VZ = DIRY*NX-NY*DIRX
C
      DOT = VX*DIRX2 + VY*DIRY2 + VZ*DIRZ2
c      print*,'   got a2: ',vx,vy,vz,dot
C
      POS2 = ( XX*VX + YY*VY + ZZ*VZ )/DOT
C
C - check if the minimum distance is too far away to worry about
C
c      print*,'   Got dst,p1,p2:',dst,pos1,pos2
      IF( DST .GT. MNDST ) THEN
         EDGDST = 1E+30
c         print*,'Leaving EDGDST...'
c         print*,'  Beyond checking'
         RETURN
      END IF
C
C - check if closest approach for other edge is within bounds, otherwise pick best node
C
      IF( POS1 .LT. -1E-5 ) THEN
         NOD1 = 2
      ELSE IF( POS1 .GT. LNGTH+1E-5 ) THEN
         NOD1 = 3
      ELSE
         NOD1 = 0
      END IF
C
C - check if closest approach for IEDGE is within bounds, otherwise pick best node
C
      NOD2 = 0
      IF( POS2 .LT. -1E-5 ) THEN
         NOD2 = 1
      ELSE IF( POS2 .GT. LNGTH2+1E-5 ) THEN
         NOD2 = 2
      ELSE IF( NOD1 .EQ. 0 ) THEN
C
C - use minimum distance from edge IEDGE to face edge
C
         EDGDST = DST
c         print*,'Leaving EDGDST...'
         if( dst .lt. mnallw ) then
            print*,'    +++ Edge to edge:',dst
            print*,'         edge nodes: ',nods(1),nods(2)
c            print*,'         normal: ',onx,ony,onz
c            print*,'         facnd1: ',x(1),y(1),z(1)
c            print*,'         facnd2: ',x(2),y(2),z(2)
c            print*,'         facnd3: ',x(3),y(3),z(3)
c            print*,'         dirn1 : ',dirx,diry,dirz
c            print*,'         dirn2 : ',dirx2,diry2,dirz2
c            print*,'         ndtofc: ',xx,yy,zz
c            print*,'         lengths: ',lngth,lngth2
c            print*,'         p1,p2,dst:',pos1,pos2,dst
c            xx = nodlst(1,nods(1)) + pos2*dirx2
c            yy = nodlst(2,nods(1)) + pos2*diry2
c            zz = nodlst(3,nods(1)) + pos2*dirz2
c            nx = x(2) + pos1*dirx
c            ny = y(2) + pos1*diry
c            nz = z(2) + pos1*dirz
c            print*,'         nearest 1: ',nx,ny,nz
c            print*,'         nearest 2: ',xx,yy,zz
c            xx = xx - nx
c            yy = yy - ny
c            zz = zz - nz
c            dst = xx*xx + yy*yy + zz*zz
c            print*,'         calculated distance: ',dst
         end if
         RETURN
      ELSE
C
C - use minimum distance from face node (NOD1) to edge IEDGE
C
         XX = X(NOD1) - NODLST(1,NODS(1))
         YY = Y(NOD1) - NODLST(2,NODS(1))
         ZZ = Z(NOD1) - NODLST(3,NODS(1))
         DOT = XX*DIRX2 + YY*DIRY2 + ZZ*DIRZ2
         XX = XX - DIRX2*DOT
         YY = YY - DIRY2*DOT
         ZZ = ZZ - DIRZ2*DOT
         DOT = XX*XX + YY*YY + ZZ*ZZ
         EDGDST = SQRT(DOT)
c         print*,'Leaving EDGDST...'
         dd = SQRT(DOT)
         if( dd .lt. mnallw ) then
            print*,'    +++ Face node to edge:',dd,nod1
            print*,'         edge nodes: ',nods(1),nods(2)
            print*,'         dst,p1,p2:',dst,pos1,pos2
         end if
         RETURN
      END IF
C
      IF( NOD1 .EQ. 0 ) THEN
C
C - use minimum distance from edge IEDGE node (NOD2) to face edge
C
         XX = X(1) - NODLST(1,NODS(NOD2))
         YY = Y(1) - NODLST(2,NODS(NOD2))
         ZZ = Z(1) - NODLST(3,NODS(NOD2))
         DOT = XX*DIRX + YY*DIRY + ZZ*DIRZ
         XX = XX + DIRX*DOT
         YY = YY + DIRY*DOT
         ZZ = ZZ + DIRZ*DOT
         DOT = XX*XX + YY*YY + ZZ*ZZ
         EDGDST = SQRT(DOT)
c         print*,'Leaving EDGDST...'
         dd = SQRT(DOT)
         if( dd .lt. mnallw ) then
            print*,'    +++ Face edge to node:',dd,nods(nod2)
            print*,'         edge nodes: ',nods(1),nods(2)
            print*,'         dst,p1,p2:',dst,pos1,pos2
         end if
c         XX = X(1) - NODLST(1,NODS(NOD2))
c         YY = Y(1) - NODLST(2,NODS(NOD2))
c         ZZ = Z(1) - NODLST(3,NODS(NOD2))
c         EDGDST = ABS(XX*ONX + YY*ONY + ZZ*ONZ)
cc         print*,'Leaving EDGDST...'
c         dd = ABS(XX*ONX + YY*ONY + ZZ*ONZ)
c         if( dd .lt. mnallw ) then
c            print*,'    +++ Node to face:',dd,nods(nod2)
c            print*,'         edge nodes: ',nods(1),nods(2)
c         end if
C
      ELSE
C
C - use minimum distance from edge IEDGE node (NOD2) to face edge node (NOD1)
C
         XX = X(NOD1) - NODLST(1,NODS(NOD2))
         YY = Y(NOD1) - NODLST(2,NODS(NOD2))
         ZZ = Z(NOD1) - NODLST(3,NODS(NOD2))
         DOT = XX*XX + YY*YY + ZZ*ZZ
         EDGDST = SQRT(DOT)
c         print*,'Leaving EDGDST...'
         dd = sqrt(dot)
         if( dd .lt. mnallw ) then
            print*,'    +++ Node to node:',dd,nod1,nods(nod2)
            print*,'         edge nodes: ',nods(1),nods(2)
            print*,'         dst,p1,p2:',dst,pos1,pos2
         end if
C
      END IF
C
      RETURN
      END
C
