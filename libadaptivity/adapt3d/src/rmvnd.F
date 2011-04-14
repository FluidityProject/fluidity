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
      SUBROUTINE RMVND( BIGLST, NODLST, IPTR )
C-----------------------------------------------------------------------
C
C - This subroutine removes a node from NODLST.
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER IPTR
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER NXT, LST, I
C
c      print*,' '
c      print*,'Into RMVND...',ierr
c      print*,'     IPTR:',iptr
C
      IF( IERR .NE. 0 ) RETURN
C
cc      IF( IPTR .LE. 0 ) THEN
cc         IERR = -1
cc         PRINT*,'*** RMVND: INVALID NODE POINTER'
cc         RETURN
cc      END IF
C
      NODLST( 1, IPTR ) = 0.0
      NODLST( 2, IPTR ) = 0.0
      NODLST( 3, IPTR ) = 0.0
C
      NODLST( 6, IPTR ) = 0.0
C
      DO I = 1, 9
         NODLST(I+6,IPTR) = 0.0
      END DO
C
      NODLST( 16, IPTR ) = 0.0
      NODLST( 17, IPTR ) = 0.0
C
      LST = INT(NODLST(LSTNOD,IPTR))
      NXT = INT(NODLST(NXTNOD,IPTR))
C
      IF( IPTR .EQ. STTNOD ) STTNOD = NXT
C
      IF( LST .GT. 0 ) NODLST( NXTNOD, LST ) = FLOAT(NXT)
      IF( NXT .GT. 0 ) NODLST( LSTNOD, NXT ) = FLOAT(LST)
C
      IF( ENDNOD .EQ. IPTR ) ENDNOD = LST
C
      NODLST( LSTNOD, IPTR )   = NODLST( LSTNOD, EMTNOD )
      NODLST( NXTNOD, IPTR )   = FLOAT(EMTNOD)
C
      IF( EMTNOD .GT. 0 ) NODLST( LSTNOD, EMTNOD ) = FLOAT(IPTR)
      IF( ENDNOD .GT. 0 ) NODLST( NXTNOD, ENDNOD ) = FLOAT(IPTR)
C
      EMTNOD = IPTR
C
      NUMNDS = NUMNDS - 1
C
c      print*,'Leaving RMVND'
c
      RETURN
      END
C
