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
#include "ewrite.h"
      SUBROUTINE CHKGMY( BIGLST, NODLST, SRFGMY )
C-----------------------------------------------------------------------
C
C - SRFGMY should be TRUE if the surface mesh is to be left intact
C
C-----------------------------------------------------------------------
      use write_log
      IMPLICIT NONE
C
      LOGICAL SRFGMY
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER IPTR, CONELM(4), EDGLST(6), IELM, INEW, IFLAG, NXT, LST,
     :        I, ELPTRS(3*MAXEDE), EDPTRS(3*MAXEDE+3), NDPTRS(MAXEDE+1),
     :        NODS(4), SURFAC, NUMEDE, I1, I2, I3, I4, INTRNL, IREG,
     :        NUMEDG, NUMELM, IEDG, J, cntnod, cntedg, cntelm
C
      REAL NX1, NY1, NZ1, NX2, NY2, NZ2, DOT, CUTOFF
C
      LOGICAL GETFLG, GTNDFL, OK
C
      PARAMETER( CUTOFF = 0.999999 )
C
      OK = .TRUE.
C
      ewrite(3,*) 
     :     'Into CHKGMY... ',ierr
C
C - Loop through edges, clearing bit four of the flag
C
      cntedg = 0
      cntelm = 0
      IPTR = STTBIG
C
  10  IF( BIGLST(4,IPTR) .LT. 0 ) THEN
         NXT = BIGLST(NXTBIG,IPTR)
         CALL SETFLG( BIGLST, NODLST, IPTR, 4, .FALSE. )
         cntedg = cntedg + 1
      ELSE IF( BIGLST(4,IPTR) .GT. 0 ) THEN
         CALL ELMINF( BIGLST, NODLST, IPTR, CONELM, EDGLST,
     :                INEW, IFLAG, IREG, NXT, LST )
         cntelm = cntelm + 1
      END IF
C
      IPTR = NXT
C
      IF( GETFLG(BIGLST,NODLST,IPTR,1) ) GOTO 10
c
      ewrite(2,*) 'Edge flags cleared: ',cntedg
      ewrite(2,*) 'Elements counted: ',cntelm
c
      if( ierr .ne. 0 ) goto 900
C
C - Loop through nodes, clearing bit six of the flag
C
      cntnod = 0
      IPTR = STTNOD
C
 15   CALL STNDFL( BIGLST, NODLST, IPTR, 6, .FALSE. )
      NXT = INT(NODLST(NXTNOD,IPTR))
      cntnod = cntnod + 1
C
      IPTR = NXT
C
      IF( GTNDFL(BIGLST,NODLST,IPTR,1) ) GOTO 15
c
      ewrite(2,*) 'Node flags cleared: ',cntnod
c
      if( ierr .ne. 0 ) goto 900
C
C - Loop through elements, picking out edges
C
      IPTR = STTBIG
C
  20  IF( BIGLST(4,IPTR) .LT. 0 ) THEN
C
         NXT = BIGLST(NXTBIG,IPTR)
C
      ELSE IF( BIGLST(4,IPTR) .GT. 0 ) THEN
C
c         print*,'Checking element ',iptr
         CALL ELMINF( BIGLST, NODLST, IPTR, CONELM, EDGLST,
     :                INEW, IFLAG, IREG, NXT, LST )
C
         DO I = 1, 6
C
            IF( .NOT. GETFLG(BIGLST,NODLST,EDGLST(I),4) ) THEN
C
C - for edges not already done, find the surrounding elements
C
               IEDG = EDGLST(I)
c               print*,'  checking edge ',iedg
C
               CALL EDGELS( BIGLST, NODLST, IEDG, IPTR, ELPTRS,
     :                      EDPTRS, NDPTRS, NUMEDE, SURFAC, INTRNL )
C
               IF( IERR .NE. 0 ) RETURN
               I1 = BIGLST(1,IEDG)
               I2 = BIGLST(2,IEDG)
C
C - and work out what sort of edge it is (internal/surface/geometry)
C
c               CALL SETFLG( BIGLST, NODLST, IEDG, 2, .FALSE. )
               IF( SURFAC .GT. 0 ) THEN
c                  CALL SETFLG( BIGLST, NODLST, IEDG, 3, .FALSE. )
                  IF( INTRNL .LT. 0 ) THEN
c                     CALL SETFLG( BIGLST, NODLST, IEDG, 2, .TRUE. )
                     if( .not. getflg(BIGLST,NODLST,iedg,2)
     :                    .or. getflg(BIGLST,NODLST,iedg,3) ) then
                        ewrite(1,*) 
     :                       'Edge ',IEDG,numede,surfac,intrnl
                        ewrite(1,*) 
     :                       'Geometry edge: ',i1,i2,dot,
     :                       getflg(BIGLST,NODLST,iedg,2),
     :                       getflg(BIGLST,NODLST,iedg,3)
                        ok = .false.
                     end if
                  ELSE IF( SRFGMY ) THEN
                     if( .not. getflg(BIGLST,NODLST,iedg,2) ) then
                        ewrite(1,*) 
     :                       'Edge ',IEDG,numede,surfac,intrnl
                        ewrite(1,*) 
     :                       'Geometry edge: ',i1,i2,dot,
     :                       getflg(BIGLST,NODLST,iedg,2),
     :                       getflg(BIGLST,NODLST,iedg,3)
                        ok = .false.
                     end if
                  ELSE
                     I1 = BIGLST(1,IEDG)
                     I2 = BIGLST(2,IEDG)
                     I3 = NDPTRS(NUMEDE)
                     I4 = NDPTRS(NUMEDE+1)
                     CALL FACNRM( BIGLST, NODLST,
     :                            I1, I2, I3, NX1, NY1, NZ1 )
                     CALL FACNRM( BIGLST, NODLST,
     :                            I2, I1, I4, NX2, NY2, NZ2 )
                     DOT = NX1*NX2 + NY1*NY2 + NZ1*NZ2
                     IF( ABS(DOT) .LT. CUTOFF ) THEN
c                        CALL SETFLG( BIGLST, NODLST, IEDG, 2, .TRUE. )
                        if( .not. getflg(BIGLST,NODLST,iedg,2) ) then
                           ewrite(1,*) 
     :                          'Edge ',IEDG,numede,surfac,intrnl
                           ewrite(1,*) 
     :                          'Geometry edge: ',i1,i2,dot,
     :                          getflg(BIGLST,NODLST,iedg,2),
     :                          getflg(BIGLST,NODLST,iedg,3)
                           ok = .false.
                        end if
                     else
                        if( getflg(BIGLST,NODLST,iedg,2) ) then
                           ewrite(1,*) 
     :                          'Edge ',
     :                          IEDG,numede,surfac,intrnl
                           ewrite(1,*) 
     :                          'Surface edge: ',i1,i2,dot,
     :                          getflg(BIGLST,NODLST,iedg,2),
     :                          getflg(BIGLST,NODLST,iedg,3)
                           ok = .false.
                        end if
                     END IF
                  END IF
               ELSE IF( INTRNL .LT. -1 ) THEN
c                  CALL SETFLG( BIGLST, NODLST, IEDG, 3, .FALSE. )
                  I1 = BIGLST(1,IEDG)
                  I2 = BIGLST(2,IEDG)
                  I3 = NDPTRS(-INTRNL-1)
                  I4 = NDPTRS(-SURFAC)
                  CALL FACNRM( BIGLST, NODLST,
     :                         I1, I2, I3, NX1, NY1, NZ1 )
                  CALL FACNRM( BIGLST, NODLST,
     :                         I2, I1, I4, NX2, NY2, NZ2 )
                  DOT = NX1*NX2 + NY1*NY2 + NZ1*NZ2
                  IF( ABS(DOT) .LT. CUTOFF ) THEN
c                     CALL SETFLG( BIGLST, NODLST, IEDG, 2, .TRUE. )
                     if( .not. getflg(BIGLST,NODLST,iedg,2)
     :                    .or. getflg(BIGLST,NODLST,iedg,3) ) then
                        ewrite(1,*) 
     :                       'Edge ',IEDG,numede,surfac,intrnl
                        ewrite(1,*) 
     :                       'Int.geom. edge: ',i1,i2,dot,
     :                       getflg(BIGLST,NODLST,iedg,2),
     :                       getflg(BIGLST,NODLST,iedg,3)
                        ok = .false.
                     end if
                  else
                     if( getflg(BIGLST,NODLST,iedg,2) .or.
     :                   getflg(BIGLST,NODLST,iedg,3) ) then
                        ewrite(1,*)
     :                       'Edge ',IEDG,numede,surfac,intrnl
                        ewrite(1,*) 
     :                       'Int.srf. edge: ',i1,i2,dot,
     :                       getflg(BIGLST,NODLST,iedg,2),
     :                       getflg(BIGLST,NODLST,iedg,3)
                        ok = .false.
                     end if
                  END IF
               ELSE IF( INTRNL .EQ. -1 ) THEN
c                  CALL SETFLG( BIGLST, NODLST, IEDG, 3, .FALSE. )
c                  CALL SETFLG( BIGLST, NODLST, IEDG, 2, .TRUE. )
                  if( .not. getflg(BIGLST,NODLST,iedg,2)
     :                 .or. getflg(BIGLST,NODLST,iedg,3) ) then
                     ewrite(1,*) 
     :                    'Edge ',IEDG,numede,surfac,intrnl
                     ewrite(1,*) 
     :                    'Int.geom. edge: ',i1,i2,
     :                    getflg(BIGLST,NODLST,iedg,2),
     :                    getflg(BIGLST,NODLST,iedg,3)
                     ok = .false.
                  end if
               ELSE
c                  CALL SETFLG( BIGLST, NODLST, IEDG, 3, .TRUE. )
                  if( .not. getflg(BIGLST,NODLST,iedg,3) ) then
                     ewrite(1,*) 
     :                    'Edge ',IEDG,numede,surfac,intrnl
                     ewrite(1,*) 
     :                    'Internal edge: ',i1,i2,
     :                    getflg(BIGLST,NODLST,iedg,2),
     :                    getflg(BIGLST,NODLST,iedg,3)
                     ok = .false.
                  end if
               END IF
C
               CALL SETFLG( BIGLST, NODLST, IEDG, 4, .TRUE. )
C
            END IF
C
         END DO
C
      END IF
C
      IPTR = NXT
C
      IF( GETFLG(BIGLST,NODLST,IPTR,1) ) GOTO 20
C
      IPTR = STTBIG
C
C - Loop through elements, picking out nodes
C
  30  IF( BIGLST(4,IPTR) .LT. 0 ) THEN
C
         NXT = BIGLST(NXTBIG,IPTR)
C
      ELSE IF( BIGLST(4,IPTR) .GT. 0 ) THEN
C
         CALL ELNODS( BIGLST, NODLST, IPTR, NODS, NXT, LST, .false. )
C
         DO I = 1, 4
C
            IF( .NOT. GTNDFL(BIGLST,NODLST,NODS(I),6) ) THEN
C
C - for nodes not already done, find list of surrounding elements/edges/nodes
C
               CALL NDCNCT( BIGLST, NODLST, NODS(I), IPTR, ELPTRS, 
     :                      EDPTRS, NDPTRS, NUMEDG, NUMELM )
C
               IF( IERR .NE. 0 ) RETURN
C
C - and work out what sort of node it is from the types of edges
C
               I2 = 0
               I3 = 0
               I4 = 0
C
               DO I1 = 1, NUMEDG
                  IF( GETFLG(BIGLST,NODLST,EDPTRS(I1),2) ) THEN
                     I2 = I2 + 1
                  ELSE IF( GETFLG(BIGLST,NODLST,EDPTRS(I1),3) ) THEN
                     I3 = I3 + 1
                  END IF
               END DO
C
c               CALL STNDFL(BIGLST,NODLST,NODS(I),2,.FALSE.)
c               CALL STNDFL(BIGLST,NODLST,NODS(I),3,.FALSE.)
c               CALL STNDFL(BIGLST,NODLST,NODS(I),4,.FALSE.)
               IF( I2 .GT. 2 ) THEN
c                  CALL STNDFL(BIGLST,NODLST,NODS(I),2,.TRUE.)
                  if( .not. gtndfl(BIGLST,NODLST,nods(i),2) ) then
                   ewrite(1,*) 
     :                    'Geometry node: ',nods(i),
     :                    (nodlst(j,nods(i)),j=1,3),
     :                    (gtndfl(BIGLST,NODLST,nods(i),j),j=2,4)
                     ok = .false.
                  end if
               ELSE IF( I3 .EQ. NUMEDG ) THEN
                  CALL STNDFL(BIGLST,NODLST,NODS(I),3,.TRUE.)
                  if( .not. gtndfl(BIGLST,NODLST,nods(i),3) ) then
                   ewrite(1,*)
     :                    'Internal node: ',nods(i),
     :                    (nodlst(j,nods(i)),j=1,3),
     :                    (gtndfl(BIGLST,NODLST,nods(i),j),j=2,4)
                     ok = .false.
                  end if
               ELSE IF( I2 .EQ. 2 ) THEN
c                  CALL STNDFL(BIGLST,NODLST,NODS(I),4,.TRUE.)
                  if( .not. gtndfl(BIGLST,NODLST,nods(i),4) ) then
                   ewrite(1,*) 
     :                    'Splitter node: ',nods(i),
     :                    (nodlst(j,nods(i)),j=1,3),
     :                    (gtndfl(BIGLST,NODLST,nods(i),j),j=2,4)
                     ok = .false.
                  end if
               END IF
C
               CALL STNDFL(BIGLST,NODLST,NODS(I),6,.TRUE.)
C
            END IF
C
         END DO
C
      END IF
C
      IPTR = NXT
C
      IF( GETFLG(BIGLST,NODLST,IPTR,1) ) GOTO 30
c
 900  if( ierr .ne. 0 ) ok = .false.
C
c      if( .not. ok ) then
c         call system('beep')
c         pause
c      end if
      ewrite(3,*) 'Leaving CHKGMY'
C
      RETURN
      END
C
