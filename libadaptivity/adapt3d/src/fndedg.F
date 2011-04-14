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
      SUBROUTINE FNDEDG( I1, I2, ENLBAS, ENLIST, NELBAS, NELIST,
     :                   NNOD, NELM, SZENLS, SZNELS, IELM,
     :                   GREATR, GEOM3D )
C-----------------------------------------------------------------------
C
C - This subroutine looks for elements that have nodes I1 and I2,
C - and returns the first one it finds in IELM - unless IELM was sent in
C - as that element, in which case it returns the next one.
C - If GREATR and GEOM3D are true then it will return an element greater
C - than IELM. If GREATR is false and GEOM3D is true then it will return
C - an element less than IELM.
C - GREATR is ignored in the 2d case (GEOM3D false).
C - If no satisfactory element is found it returns IELM as -1.
C
C - This is supposed to be used in 2d (GEOM3D false) to find connected
C - elements - by supplying the original element in IELM, and the nodes
C - of an edge in I1 and I2, it will return in IELM the element that is
C - connected to the original IELM via the edge with nodes I1 and I2.
C
C - It is also for use in 3d (GEOM3D true) for finding elements around
C - an edge, by setting GREATR to true and sending in successively larger
C - elements until -1 is returned, then using GREATR false to find the
C - smaller ones until -1 is returned.
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER I1, I2, NNOD, NELM, SZENLS, SZNELS, IELM
C
      INTEGER ENLBAS(NELM+1), NELBAS(NNOD+1), ENLIST(SZENLS),
     :        NELIST(SZNELS)
C
      INTEGER I, J, ICONEL
C
      LOGICAL GREATR, GEOM3D
C
c      PRINT*,' '
c      print*,'Into FNDEDG...'
c      PRINT*,'   I1, I2: ',I1,I2
c      PRINT*,'     IELM: ',IELM
C
      IF( GEOM3D .AND. GREATR ) THEN
C
        DO I = NELBAS(I1)+1, NELBAS(I1+1)
C
          ICONEL = NELIST(I)
c         print*,'     checking ',iconel,enlbas(iconel),ENLBAS(ICONEL+1)
C
          IF( ICONEL .GT. IELM ) THEN
C
            DO J = ENLBAS(ICONEL)+1, ENLBAS(ICONEL+1)
               IF( ENLIST(J) .EQ. I2 ) THEN
                  IELM = ICONEL
c      print*,'Leaving FNDEDG... Found IELM:', IELM
                  GOTO 100
               END IF
            END DO
C
          END IF
C
        END DO
C
      ELSE IF( GEOM3D .AND. .NOT. GREATR ) THEN
C
        DO I = NELBAS(I1)+1, NELBAS(I1+1)
C
          ICONEL = NELIST(I)
c         print*,'     checking ',iconel,enlbas(iconel),ENLBAS(ICONEL+1)
C
          IF( ICONEL .LT. IELM ) THEN
C
            DO J = ENLBAS(ICONEL)+1, ENLBAS(ICONEL+1)
               IF( ENLIST(J) .EQ. I2 ) THEN
                  IELM = ICONEL
c      print*,'Leaving FNDEDG... Found IELM:', IELM
                  GOTO 100
               END IF
            END DO
C
          END IF
C
        END DO
C
      ELSE
C
        DO I = NELBAS(I1)+1, NELBAS(I1+1)
C
          ICONEL = NELIST(I)
c         print*,'     checking ',iconel,enlbas(iconel),ENLBAS(ICONEL+1)
C
          IF( ICONEL .NE. IELM ) THEN
C
            DO J = ENLBAS(ICONEL)+1, ENLBAS(ICONEL+1)
               IF( ENLIST(J) .EQ. I2 ) THEN
                  IELM = ICONEL
c      print*,'Leaving FNDEDG... Found IELM:', IELM
                  GOTO 100
               END IF
            END DO
C
          END IF
C
        END DO
C
      END IF
C
c      print*,'FAILED! Leaving FNDEDG... IELM:', IELM
      IELM = -1
C
  100 CONTINUE
C
      RETURN
      END
C
