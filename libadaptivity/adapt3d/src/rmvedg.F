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
      SUBROUTINE RMVEDG( BIGLST, NODLST, IPTR )
C-----------------------------------------------------------------------
C
C - This subroutine removes an edge from BIGLST.
C
C - IMPORTANT NOTE: you should ALWAYS make sure that any elements that
C - have this edge are also removed a.s.a.p. when calling this subroutine.
C - This subroutine makes no attempt to correct for the elements that
C - may have this edge. Use EDGELS to find which elements have an edge
C - BEFORE removing it!
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blknxt.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER IPTR
C
      INTEGER NODS(2), INEWED, NXT, LST, IFLAG
C
c      if( iptr .eq. 32454 ) then
c         print*,'>>>>>>>>> removing edge ',iptr,' <<<<<<<<<'
c         print*,'   nodes: ',nods
c         pause
c      end if
c      print*,nextls
c      print*,' '
c      print*,'Into RMVEDG...',ierr
c      print*,'     IPTR:',iptr
C
c      IF( IERR .NE. 0 ) RETURN
C
c      CALL EDGINF( BIGLST, NODLST, IPTR, NODS, INEWED, IFLAG, NXT, LST )
C
      IF( IERR .NE. 0 ) RETURN
C
      LST = BIGLST( LSTBIG, IPTR )
      NXT = BIGLST( NXTBIG, IPTR )
C
cc      IF( BIGLST(4,IPTR) .GT. 0 ) THEN
cc         IERR = -2
cc         PRINT*,'*** RMVEDG: NOT AN EDGE BLOCK'
cc         RETURN
cc      END IF
C
cc      IF( BIGLST(4,IPTR) .EQ. 0 ) THEN
cc         IERR = -1
cc         PRINT*,'*** RMVEDG: BLOCK IS ALREADY EMPTY'
cc         print*,iptr
cc         RETURN
cc      END IF
C
cc      IF( BIGLST(1,IPTR) .LE. 0 ) THEN
cc         IERR = -2
cc         PRINT*,'*** RMVEDG: FOUND INVALID NODE 1'
cc         print*,iptr,biglst(1,iptr),biglst(2,iptr)
cc         RETURN
cc      END IF
C
cc      IF( BIGLST(2,IPTR) .LE. 0 ) THEN
cc         IERR = -2
cc         PRINT*,'*** RMVEDG: FOUND INVALID NODE 2'
cc         print*,iptr,biglst(1,iptr),biglst(2,iptr)
cc         RETURN
cc      END IF
C
      BIGLST( 1, IPTR ) = 0
      BIGLST( 2, IPTR ) = 0
      BIGLST( 3, IPTR ) = 0
      BIGLST( 4, IPTR ) = 0
C
c      IF( INEWED .GT. 0 ) THEN
c         NEWEDG(INEWED) = NEWEDG(NNEWED)
c         NNEWED = NNEWED - 1
c      END IF
C
      DO IFLAG = 1, MXNEXT
         IF( NEXTLS(IFLAG) .EQ. IPTR .AND. ENDBIG .NE. IPTR ) THEN
c            print*,'>>>> RMVELE: REMOVING NEXT EDGE: ',
c     :              IPTR,NXT,NELEMS,NEDGES
            NEXTLS(IFLAG) = NXT
         ELSE IF( NEXTLS(IFLAG) .EQ. EMTBIG ) THEN
c            print*,'>>>> RMVELE: REMOVING LAST EDGE: ',
c     :              IPTR,EMTBIG,NELEMS,NEDGES
            NEXTLS(IFLAG) = IPTR
         END IF
      END DO
C
      IF( IPTR .EQ. STTBIG ) STTBIG = NXT
C
      IF( LST .GT. 0 ) BIGLST( NXTBIG, LST ) = NXT
      IF( NXT .GT. 0 ) BIGLST( LSTBIG, NXT ) = LST
C
      IF( ENDBIG .EQ. IPTR ) ENDBIG = LST
C
      BIGLST( LSTBIG, IPTR ) = ENDBIG
      BIGLST( NXTBIG, IPTR ) = EMTBIG
C
      IF( EMTBIG .GT. 0 ) BIGLST( LSTBIG, EMTBIG ) = IPTR
      IF( ENDBIG .GT. 0 ) BIGLST( NXTBIG, ENDBIG ) = IPTR
C
      BIGLST( 7, IPTR ) = 0
C
      EMTBIG = IPTR
C
      NEDGES = NEDGES - 1
c
c      print*,'>>>>>>>>> removed edge ',iptr,' <<<<<<<<<'
c      print*,nextls
C
      RETURN
      END
C
