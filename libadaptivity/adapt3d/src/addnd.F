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
      SUBROUTINE ADDND( BIGLST, NODLST,
     :                  ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                  SZENLS, SZNELS, NNOD,   NELM,
     :                  NODX,   NODY,   NODZ,   ORGMTX,
     :                  X, Y, Z, IPTR, INEWND, IELM, SEARCH )
C-----------------------------------------------------------------------
C
C - This subroutine adds a node into nodlst.
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER IPTR, INEWND, IELM
C
      REAL X, Y, Z
C
      LOGICAL SEARCH
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkorg.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER I
c
c      print*,' '
c      print*,'Into ADDND...',ierr
c      print*,'        X:',X
c      print*,'        Y:',Y
c      print*,'        Z:',Z
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** ADDND: GOT ERROR ON ENTRY'
         RETURN
      END IF
C
      IPTR = EMTNOD
      IF( IPTR .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** ADDNOD: FAILED TO ADD NODE TO NODLST - FULL'
         RETURN
      END IF
C
cc      IF( NODLST(6,IPTR) .NE. 0.0 ) THEN
cc         IERR = -2
cc         PRINT*,'*** ADDND: NODE BLOCK IS NOT EMPTY'
cc         RETURN
cc      END IF
C
      CALL FNDELM( ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :             SZENLS, SZNELS, NNOD,   NELM,
     :             NODX,   NODY,   NODZ,   ORGMTX,
     :             X, Y, Z, IELM,  SEARCH, NODLST(7,IPTR) )
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** ADDND: GOT ERROR FROM FNDELM'
         RETURN
      END IF
C
      NODLST(1,IPTR) = X
      NODLST(2,IPTR) = Y
      NODLST(3,IPTR) = Z
C
      NODLST(6,IPTR) = 1.0
      NODLST(16,IPTR) = FLOAT(ABS(IELM))
C      CALL FIXEIG( NODLST(1,IPTR), NODLST(7,IPTR), ANISOT )
C      DO I = 1, 9
C         NODLST(I+6,IPTR) = 0.0
C      END DO
C      NODLST( 7,IPTR) = 44.44444
C      NODLST(11,IPTR) = 44.44444
C      NODLST(15,IPTR) = 400.0
C
      NODLST(LSTNOD,IPTR) = FLOAT(ENDNOD)
      IF( ENDNOD .GT. 0 ) NODLST(NXTNOD,ENDNOD) = FLOAT(IPTR)
      ENDNOD = IPTR
      EMTNOD = INT(NODLST(NXTNOD,IPTR))
C
      NUMNDS = NUMNDS + 1
C
c      if( anisot .gt. 1 ) then
c       print*,'Leaving ADDND...'
c       print*,'     IPTR,IELM:',iptr,ielm
c       print*,'   INEWND:',inewnd
c       print*,'     Metric: '
c       print*,'     ',nodlst(7,iptr),nodlst(10,iptr),nodlst(13,iptr)
c       print*,'     ',nodlst(8,iptr),nodlst(11,iptr),nodlst(14,iptr)
c       print*,'     ',nodlst(9,iptr),nodlst(12,iptr),nodlst(15,iptr)
c      end if
C
      RETURN
      END
C
