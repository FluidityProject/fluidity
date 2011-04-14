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
      SUBROUTINE NODFCS( BIGLST, NODLST, IPTR, IFAC, FCPTRS,
     :                   EDPTRS, NDPTRS, CURADD, REVERS )
C-----------------------------------------------------------------------
C
C - This subroutine finds a list of faces which surround a node,
C - and their edges (in the appropriate order for the faces so that
C - the first two edges define the first face), returning the lists
C - in FCPTRS and EDPTRS and the number of them in CURADD.
C - It also returns information about the orientation of the surfaces:
C - the cross product of each adjacent pair of edges gives the normal
C - of each face - but it may point the wrong way if the faces go
C - around the wrong way. If this is the case then REVERS is returned
C - as TRUE.
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
      INTEGER FCPTRS(MAXFCS), EDPTRS(MAXFCS), NDPTRS(MAXFCS),
     :        IPTR, IFAC, CURADD
C
      LOGICAL REVERS
C
      INTEGER EDGLST(3), CONFAC(3), NXT, LST, IELM, INEW, IFLAG,
     :        CURCHK, I, FCHSND, FC1, FC2, FC3
C
      LOGICAL GTFRFL
C
c      print*,' '
c      print*,'In NODFCS... ',ierr
c      print*,'   INOD,IFAC: ',IPTR,IFAC
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( IPTR .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** NODFCS: INVALID NODE POINTER'
         RETURN
      END IF
C
      IF( IFAC .LE. 0 ) THEN
         IERR = -2
         PRINT*,'*** NODFCS: INVALID STARTING FACE POINTER'
         RETURN
      END IF
C
      CURADD = 0
      CURCHK = 1
C
      CALL FACINF( IFAC, CONFAC, EDGLST, IELM, INEW, IFLAG, NXT, LST )
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( FCHSND( BIGLST, NODLST, IFAC, IPTR ) .EQ. 0 ) THEN
         IERR = -93
         PRINT*,'*** NODFCS: STARTING FACE DOES NOT CONTAIN NODE'
         RETURN
      END IF
C
      FC1 = CONFAC(1)
      FC2 = CONFAC(2)
      FC3 = CONFAC(3)
C
C 100  CONTINUE
C
      CURADD = CURADD + 1
      FCPTRS(CURADD) = IFAC
      CALL STFRFL(IFAC,3,.TRUE.)
C
 200  CONTINUE
C
      IF( CONFAC(1) .GT. 0 ) THEN
         IF( FCHSND( BIGLST, NODLST, CONFAC(1), IPTR ) .NE. 0 ) THEN
            IF( .NOT. GTFRFL(CONFAC(1),3) ) THEN
               CURADD = CURADD + 1
               FCPTRS(CURADD) = CONFAC(1)
               EDPTRS(CURADD-1) = EDGLST(1)
               CALL STFRFL(CONFAC(1),3,.TRUE.)
               GOTO 300
            END IF
         END IF
      END IF
C
      IF( CONFAC(2) .GT. 0 ) THEN
         IF( FCHSND( BIGLST, NODLST, CONFAC(2), IPTR ) .NE. 0 ) THEN
            IF( .NOT. GTFRFL(CONFAC(2),3) ) THEN
               CURADD = CURADD + 1
               FCPTRS(CURADD) = CONFAC(2)
               EDPTRS(CURADD-1) = EDGLST(2)
               CALL STFRFL(CONFAC(2),3,.TRUE.)
               GOTO 300
            END IF
         END IF
      END IF
C
      IF( CONFAC(3) .GT. 0 ) THEN
         IF( FCHSND( BIGLST, NODLST, CONFAC(3), IPTR ) .NE. 0 ) THEN
            IF( .NOT. GTFRFL(CONFAC(3),3) ) THEN
               CURADD = CURADD + 1
               FCPTRS(CURADD) = CONFAC(3)
               EDPTRS(CURADD-1) = EDGLST(3)
               CALL STFRFL(CONFAC(3),3,.TRUE.)
               GOTO 300
            END IF
         END IF
      END IF
C
      IF( CONFAC(1) .EQ. IFAC ) THEN
         EDPTRS(CURADD) = EDGLST(1)
      ELSE IF( CONFAC(2) .EQ. IFAC ) THEN
         EDPTRS(CURADD) = EDGLST(2)
      ELSE IF( CONFAC(3) .EQ. IFAC ) THEN
         EDPTRS(CURADD) = EDGLST(3)
      ELSE
         PRINT*,'*** NODFCS: LIST OF FACES NOT CLOSED'
         IERR = -99
         GOTO 400
      END IF
C
      IF( BIGLST(1,EDPTRS(CURADD)) .EQ. IPTR ) THEN
         NDPTRS(CURADD) = BIGLST(2,EDPTRS(CURADD))
      ELSE IF( BIGLST(2,EDPTRS(CURADD)) .EQ. IPTR ) THEN
         NDPTRS(CURADD) = BIGLST(1,EDPTRS(CURADD))
      ELSE
         PRINT*,'*** NODFCS: FIRST EDGE DOES NOT HAVE NODE'
         PRINT*,'EDGE,NODE: ',EDPTRS(CURADD),IPTR
         PRINT*,'NODES:',BIGLST(1,EDPTRS(CURADD)),
     :                   BIGLST(2,EDPTRS(CURADD))
         IERR = -91
      END IF
C
      GOTO 400
C
 300  CONTINUE
C
      IF( BIGLST(1,EDPTRS(CURADD-1)) .EQ. IPTR ) THEN
         NDPTRS(CURADD-1) = BIGLST(2,EDPTRS(CURADD-1))
      ELSE IF( BIGLST(2,EDPTRS(CURADD-1)) .EQ. IPTR ) THEN
         NDPTRS(CURADD-1) = BIGLST(1,EDPTRS(CURADD-1))
      ELSE
         PRINT*,'*** NODFCS: EDGE DOES NOT HAVE NODE'
         PRINT*,'EDGE,NODE: ',EDPTRS(CURADD-1),IPTR
         PRINT*,'NODES:',BIGLST(1,EDPTRS(CURADD-1)),
     :                   BIGLST(2,EDPTRS(CURADD-1))
         IERR = -91
         GOTO 400
      END IF
C
      CALL FACINF( FCPTRS(CURADD), CONFAC, EDGLST, IELM, INEW,
     :             IFLAG, NXT, LST )
C
      GOTO 200
C
 400  CONTINUE
C
      IF( IERR .EQ. 0 ) THEN
         IF( FCPTRS(2) .EQ. FC1 ) THEN
            REVERS = ( FCPTRS(CURADD) .EQ. FC2 )
         ELSE IF( FCPTRS(2) .EQ. FC2 ) THEN
            REVERS = ( FCPTRS(CURADD) .EQ. FC3 )
         ELSE IF( FCPTRS(2) .EQ. FC3 ) THEN
            REVERS = ( FCPTRS(CURADD) .EQ. FC1 )
         ELSE
            IERR = -92
            PRINT*,'*** NODFCS: INCONSISTENCY IN FACE LIST'
         END IF
      END IF
C
      DO CURCHK = 1, CURADD
         CALL STFRFL(FCPTRS(CURCHK),3,.FALSE.)
      END DO
C
c      print*,'Leaving NODFCS... ',ierr
c      print*,'  faces: ',(FCPTRS(CURCHK),CURCHK=1,CURADD)
c      print*,'  edges: ',(EDPTRS(CURCHK),CURCHK=1,CURADD)
c      print*,'  nodes: ',(NDPTRS(CURCHK),CURCHK=1,CURADD)
c      PRINT*,'  REVERS:',REVERS
C
      RETURN
      END
C
