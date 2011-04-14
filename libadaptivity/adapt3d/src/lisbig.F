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
      SUBROUTINE LISBIG( BIGLST, NODLST, LISNDS, LISELS, LISEDS )
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
C
      LOGICAL LISNDS, LISELS, LISEDS
C
      INTEGER IPTR, NXT, LST, INEW, I, IFLAG, IREG,
     :        NODS(2), CONELM(4), EDGLST(6)
C
      LOGICAL GETFLG, GTNDFL
C
      IF( LISNDS ) THEN
C
       print*,'-------------------------------------------------------'
       print*,' '
       WRITE(*,900) 'Number of nodes: ',NUMNDS,STTNOD,ENDNOD
C
       IPTR = STTNOD
C
       DO I = 1, NUMNDS
C
         IF( GTNDFL(BIGLST, NODLST,IPTR,1) ) THEN
            WRITE(*,900) 'NODE ',IPTR
            WRITE(*,910) '       X: ',NODLST(1,IPTR)
            WRITE(*,910) '       Y: ',NODLST(2,IPTR)
            WRITE(*,910) '       Z: ',NODLST(3,IPTR)
            WRITE(*,900) '   FLAGS: ',INT(NODLST(6,IPTR))
            WRITE(*,900) ' NXT,LST: ',INT(NODLST(NXTNOD,IPTR)),
     :                                INT(NODLST(LSTNOD,IPTR))
         ELSE
            NXT = NODLST(NXTNOD,IPTR)
            LST = NODLST(LSTNOD,IPTR)
            WRITE(*,900) '*** PREMATURE EMPTY BLOCK ',IPTR
            WRITE(*,900) '      IFLAG: ',NODLST(6,IPTR)
            WRITE(*,900) '    NXT,LST: ',NXT,LST
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
      END IF
C
      IF( LISELS .OR. LISEDS ) THEN
C
       print*,'-------------------------------------------------------'
       print*,' '
       IF( LISEDS ) WRITE(*,900) 'Number of edges: ',NEDGES,STTBIG
       IF( LISELS ) WRITE(*,900) 'Number of elements: ',NELEMS,ENDBIG
C
       IPTR = STTBIG
C
       DO I = 1, NEDGES+NELEMS
C
         IF( .NOT. GETFLG(BIGLST,NODLST,IPTR,1) ) THEN
            NXT = BIGLST(NXTBIG,IPTR)
            LST = BIGLST(LSTBIG,IPTR)
            WRITE(*,900) '*** PREMATURE EMPTY BLOCK ',IPTR,I
            WRITE(*,900) '      IFLAG: ',BIGLST(4,IPTR)
            WRITE(*,900) '    NXT,LST: ',NXT,LST
         ELSE IF( BIGLST(4,IPTR) .LT. 0 ) THEN
            CALL EDGINF( BIGLST, NODLST, IPTR, NODS,
     :                   INEW, IFLAG, NXT, LST )
            IF( LISEDS ) THEN
               WRITE(*,900) 'EDGE ',IPTR,I
               WRITE(*,900) '      NODES: ',NODS
               WRITE(*,900) '  INEW,FLAG: ',INEW,IFLAG
               WRITE(*,900) '    NXT,LST: ',NXT,LST
            END IF
         ELSE IF( BIGLST(4,IPTR) .GT. 0 ) THEN
            CALL ELMINF( BIGLST, NODLST, IPTR, CONELM, EDGLST,
     :                   INEW, IFLAG, IREG, NXT, LST )
            IF( LISELS ) THEN
               WRITE(*,900) 'ELEMENT ',IPTR,I,IREG
               WRITE(*,900) '     CONELM: ',CONELM
               WRITE(*,900) '     EDGLST: ',EDGLST
               WRITE(*,900) '  INEW,FLAG: ',INEW,IFLAG
               WRITE(*,900) '    NXT,LST: ',NXT,LST
            END IF
         END IF
C
         IPTR = NXT
C
         IF( NXT .LE. 0 ) THEN
            IERR = -1
            WRITE(*,900) '*** FATAL ERROR: NXT =',NXT
            RETURN
         END IF
C
       END DO
C
      END IF
C
 900  FORMAT( A, 6I6 )
 910  FORMAT( A, 1PE13.6 )
      RETURN
      END
C
