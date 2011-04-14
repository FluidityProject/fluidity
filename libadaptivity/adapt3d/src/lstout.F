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
      SUBROUTINE LSTOUT( BIGLST, NODLST, outfil )
C-----------------------------------------------------------------------
C
C
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkerr.i'
c
      integer outfil
C
      INTEGER IPTR, I, ICNT, ICON
C
C      print*,'-------------------------------------------------------'
C      print*,' '
C      PRINT*,'Number of nodes: ',NUMNDS,STTNOD,ENDNOD
C
      IPTR = STTNOD
      ICNT = 0
      ICON = 0
C
      DO I = 1, NUMNDS
C
C         PRINT*,'NODE ',IPTR
C         PRINT*,'       X: ',NODLST(1,IPTR)
C         PRINT*,'       Y: ',NODLST(2,IPTR)
C         PRINT*,'       Z: ',NODLST(3,IPTR)
C         PRINT*,' NXT,LST: ',INT(NODLST(NXTNOD,IPTR))
C     :                      ,INT(NODLST(LSTNOD,IPTR))
         ICNT = ICNT + 1
C
         IF( NODLST(6,IPTR) .LT. 0.0 ) THEN
            ICON = ICON + 1
            WRITE( outfil, 999 ) ICNT,ICON
         END IF
C
         IPTR = INT(NODLST(NXTNOD,IPTR))
C
         IF( IPTR .LE. 0 ) THEN
            IERR = -1
            PRINT*,'*** FATAL ERROR: NXT =',IPTR
            RETURN
         END IF
C
      END DO
C
 999  FORMAT( 'CONS ',2I5 )
      RETURN
      END
C
