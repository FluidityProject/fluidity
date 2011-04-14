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
      SUBROUTINE MINANG( BIGLST, NODLST, DX, DY, DZ, INOD, IFAC, DOTN )
C-----------------------------------------------------------------------
C
C - This subroutine finds the minimum angle between a *unit* vector
C - (DX,DY,DZ) and the faces around node INOD. IFAC is used in NODFCS
C - to start the search for faces around node INOD.
C - DOTN returns the cosine of the angle between the unit vector and
C - the normal of the face with the minimum angle.
C - This should be positive if the unit vector points outside the
C - front faces around node INOD.
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkfrt.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER INOD, IFAC
C
      REAL DX, DY, DZ, DOTN
C
      INTEGER FCPTRS(MAXFCS), EDPTRS(MAXFCS), NDPTRS(MAXFCS), NUMFAC,
     :        MAXFAC, MINEDG, I
C
      REAL MAXDOT, DOT, XO, YO, ZO, DXL, DYL ,DZL, DXN, DYN, DZN,
     :     NX, NY, NZ, QX, QY, QZ, VX, VY, VZ, MAXNX, MAXNY, MAXNZ,
     :     DOT1, DOT2, DXA, DYA, DZA, MINDOT
C
      LOGICAL BTWEEN, REVERS
C
c      print*,' '
c      print*,'In MINANG... ',ierr
c      print*,'   INOD,IFAC: ',INOD,IFAC
c      PRINT*,'   DX,DY,DZ: ',DX,DY,DZ
C
      CALL NODFCS( BIGLST, NODLST, INOD, IFAC, FCPTRS,
     :             EDPTRS, NDPTRS, NUMFAC, REVERS )
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( NUMFAC .LE. 0 ) THEN
         IERR = -80
         PRINT*,'*** MINANG: NO FACES AROUND NODE'
         PRINT*,INOD,IFAC,NUMFAC
         RETURN
      END IF
C
      MAXFAC = 0
      MINDOT = -1.0
C
      MINEDG = 0
      MAXDOT = 1E-8
C
      MAXNX  = 0.0
      MAXNY  = 0.0
      MAXNZ  = 0.0
C
      XO = NODLST(1,INOD)
      YO = NODLST(2,INOD)
      ZO = NODLST(3,INOD)
C
      DXL = NODLST(1,NDPTRS(NUMFAC)) - XO
      DYL = NODLST(2,NDPTRS(NUMFAC)) - YO
      DZL = NODLST(3,NDPTRS(NUMFAC)) - ZO
C
      DOT = DXL*DXL + DYL*DYL + DZL*DZL
      DOT = 1.0/SQRT(DOT)
C
      DXL = DXL*DOT
      DYL = DYL*DOT
      DZL = DZL*DOT
C
      DO I = 1, NUMFAC
C
c         print*,'>>> Doing face ',fcptrs(i)
c
         DXN = NODLST(1,NDPTRS(I)) - XO
         DYN = NODLST(2,NDPTRS(I)) - YO
         DZN = NODLST(3,NDPTRS(I)) - ZO
C
         DOT = DXN*DXN + DYN*DYN + DZN*DZN
         DOT = 1.0/SQRT(DOT)
C
         DXN = DXN*DOT
         DYN = DYN*DOT
         DZN = DZN*DOT
C
         DOT = DX*DXN + DY*DYN + DZ*DZN
c         PRINT*,'    edge ang: ',dot,edptrs(i),ndptrs(i)
C
         IF( DOT .GT. MINDOT ) THEN
            MINDOT = DOT
            MINEDG = I
            IF( MINDOT .GT. SQRT(MAXDOT) ) THEN
c               print*,'Found smaller angle with edge ',i
c               print*,'mindot,maxdot: ',mindot,maxdot
               MAXFAC = 0
               MAXDOT = MINDOT*MINDOT
            END IF
         END IF
c
c         print*,'    last line: ',dxl,dyl,dzl
c         print*,'    this line: ',dxn,dyn,dzn
C
         NX = DYL*DZN - DZL*DYN
         NY = DZL*DXN - DXL*DZN
         NZ = DXL*DYN - DYL*DXN
C
         DOT = NX*NX + NY*NY + NZ*NZ
         DOT = 1.0/SQRT(DOT)
C
         IF( REVERS ) DOT = -DOT
C
         NX = NX*DOT
         NY = NY*DOT
         NZ = NZ*DOT
c         print*,'    normal: ',nx,ny,nz
C
         QX = NY*DZ - DY*NZ
         QY = NZ*DX - DZ*NX
         QZ = NX*DY - DX*NY
C
         DOTN = QX*QX + QY*QY + QZ*QZ
c         print*,'    DOTN: ',dotn
C
         IF( DOTN .LT. MAXDOT ) GOTO 900
C
         IF( DOTN .LT. MINDOT*MINDOT ) GOTO 900
C
         DOT = 1.0/SQRT(DOTN)
C
         QX = QX*DOT
         QY = QY*DOT
         QZ = QZ*DOT
C
         VX = NY*QZ - QY*NZ
         VY = NZ*QX - QZ*NX
         VZ = NX*QY - QX*NY
C
         DOT = VX*DX + VY*DY + VZ*DZ
C
         IF( DOT .LT. 0.0 ) THEN
            VX = -VX
            VY = -VY
            VZ = -VZ
         END IF
c         print*,'    Q: ',qx,qy,qz
c         print*,'    V: ',vx,vy,vz
C
         DOT  = DXL*DXN + DYL*DYN + DZL*DZN
         DOT1 = DXL*VX  + DYL*VY  + DZL*VZ
         DOT2 = DXN*VX  + DYN*VY  + DZN*VZ
c         print*,'    dot,dot1,dot2: ',dot,dot1,dot2
C
         DXA = DXL + DXN
         DYA = DYL + DYN
         DZA = DZL + DZN
C
         BTWEEN = .TRUE.
C
         IF( DOT1 .LT. DOT ) THEN
C
            BTWEEN = .FALSE.
C
         ELSE IF( DOT2 .LT. DOT ) THEN
C
            BTWEEN = .FALSE.
C
         ELSE IF( DOT .LT. 0.0 ) THEN
C
            BTWEEN = ( DXA*VX + DYA*VY + DZA*VZ .GT. 0.0 )
C
         END IF
C
         IF( BTWEEN ) THEN
C
            MAXDOT = DOTN
            MAXFAC = I
            MAXNX  = NX
            MAXNY  = NY
            MAXNZ  = NZ
c            print*,'    * Got between * '
C
c         ELSE IF( MAXFAC .LE. 0 ) THEN
C
c            BTWEEN = .TRUE.
c            DOT1 = -DOT1
c            DOT2 = -DOT2
C
c            IF( DOT1 .LT. DOT ) THEN
C
c               BTWEEN = .FALSE.
C
c            ELSE IF( DOT2 .LT. DOT ) THEN
C
c               BTWEEN = .FALSE.
C
c            ELSE IF( DOT .LT. 0.0 ) THEN
C
c               BTWEEN = ( DXA*VX + DYA*VY + DZA*VZ .GT. 0.0 )
C
c            END IF
C
c            IF( BTWEEN ) THEN
C
c               MAXDOT = DOTN
c               MAXFAC = -I
c               MAXNX  = NX
c               MAXNY  = NY
c               MAXNZ  = NZ
c               print*,'    * Got between opposite * '
C
c            END IF
C
         END IF
C
 900     CONTINUE
C
         DXL = DXN
         DYL = DYN
         DZL = DZN
C
      END DO
C
C - catch the case where the vector is nearest to an edge rather than
C - falling between the edges of a face.
C
      IF( MAXFAC .EQ. 0 .AND. MINEDG .GT. 0 ) THEN
C
c         print*,'Using MINEDG: ',edptrs(minedg),minedg
         I = MINEDG-1
         IF( I .EQ. 0 ) I = NUMFAC
C
         DXL = NODLST(1,NDPTRS(I)) - XO
         DYL = NODLST(2,NDPTRS(I)) - YO
         DZL = NODLST(3,NDPTRS(I)) - ZO
C
         DOT = DXL*DXL + DYL*DYL + DZL*DZL
         DOT = 1.0/SQRT(DOT)
C
         DXL = DXL*DOT
         DYL = DYL*DOT
         DZL = DZL*DOT
C
         DXA = NODLST(1,NDPTRS(MINEDG)) - XO
         DYA = NODLST(2,NDPTRS(MINEDG)) - YO
         DZA = NODLST(3,NDPTRS(MINEDG)) - ZO
C
         DOT = DXA*DXA + DYA*DYA + DZA*DZA
         DOT = 1.0/SQRT(DOT)
C
         DXA = DXA*DOT
         DYA = DYA*DOT
         DZA = DZA*DOT
C
         I = MINEDG+1
         IF( I .EQ. NUMFAC+1 ) I = 1
C
         DXN = NODLST(1,NDPTRS(I)) - XO
         DYN = NODLST(2,NDPTRS(I)) - YO
         DZN = NODLST(3,NDPTRS(I)) - ZO
C
         DOT = DXN*DXN + DYN*DYN + DZN*DZN
         DOT = 1.0/SQRT(DOT)
C
         DXN = DXN*DOT
         DYN = DYN*DOT
         DZN = DZN*DOT
C
         DOT = DXL*DXA + DYL*DYA + DZL*DZA
C
         DOT1 = SQRT(1.0-DOT*DOT)
C
         DOT2 = 0.05/DOT1
C
         DOT1 = 0.9987492 - DOT2*DOT
C
         VX  = DOT1*DXA + DOT2*DXL
         VY  = DOT1*DYA + DOT2*DYL
         VZ  = DOT1*DZA + DOT2*DZL
c         print*,'Got V1: ',vx,vy,vz
C
         DOT = DXN*DXA + DYN*DYA + DZN*DZA
C
         DOT1 = SQRT(1.0-DOT*DOT)
C
         DOT2 = 0.05/DOT1
C
         DOT1 = 0.9987492 - DOT2*DOT
C
         QX  = DOT1*DXA + DOT2*DXN
         QY  = DOT1*DYA + DOT2*DYN
         QZ  = DOT1*DZA + DOT2*DZN
c         print*,'Got V2: ',qx,qy,qz
C
         DOT1 = VX*DX + VY*DY + VZ*DZ
         DOT2 = QX*DX + QY*DY + QZ*DZ
C
         IF( DOT1 .GT. DOT2 ) THEN
            MAXFAC = MINEDG
            MAXNX = DYL*DZA - DZL*DYA
            MAXNY = DZL*DXA - DXL*DZA
            MAXNZ = DXL*DYA - DYL*DXA
            DOT = MAXNX*MAXNX + MAXNY*MAXNY + MAXNZ*MAXNZ
            IF( REVERS ) THEN
               DOT = -1.0/SQRT(DOT)
            ELSE
               DOT = 1.0/SQRT(DOT)
            END IF
            MAXNX = MAXNX*DOT
            MAXNY = MAXNY*DOT
            MAXNZ = MAXNZ*DOT
         ELSE
            MAXFAC = MINEDG+1
            IF( MAXFAC .GT. NUMFAC ) MAXFAC = 1
            MAXNX = DYN*DZA - DZN*DYA
            MAXNY = DZN*DXA - DXN*DZA
            MAXNZ = DXN*DYA - DYN*DXA
            DOT = MAXNX*MAXNX + MAXNY*MAXNY + MAXNZ*MAXNZ
            IF( REVERS ) THEN
               DOT = 1.0/SQRT(DOT)
            ELSE
               DOT = -1.0/SQRT(DOT)
            END IF
            MAXNX = MAXNX*DOT
            MAXNY = MAXNY*DOT
            MAXNZ = MAXNZ*DOT
         END IF
c         print*,'Using face ',fcptrs(maxfac),dot1,dot2
C
      END IF
C
      IF( MAXFAC .NE. 0 ) THEN
         DOTN = MAXNX*DX + MAXNY*DY + MAXNZ*DZ
      ELSE
         IERR = -81
         PRINT*,'*** MINANG: FAILED TO FIND MINIMUM ANGLE FACE'
      END IF
C
c      PRINT*,'Leaving MINANG...',ierr
c      print*,'   DOTN: ',DOTN
c      print*,'  (maxfac: ',fcptrs(abs(maxfac)),maxfac,' )'
c      print*,'  (maxnrm: ',maxnx,maxny,maxnz,' )'
C
      RETURN
      END
C
