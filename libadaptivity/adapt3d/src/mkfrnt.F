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
      SUBROUTINE MKFRNT( BIGLST, NODLST,
     :                   NODX, NODY, NODZ, MAXNOD, ELMSUR, NDSTSF,
     :                   ENLBAS, ENLIST, NELBAS, NELIST, EELIST, CNCTSF,
     :                   NELM, NNOD, SZENLS, SZNELS, GEOM3D )
C-----------------------------------------------------------------------
C
C - This subroutine creates the adaptive-front structures (edges in BIGLST and
C - faces in FRTLST) required to fill in the domain with a first crude mesh.
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
      INCLUDE 'blknwf.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER NELM, SZENLS, SZNELS, NNOD, MAXNOD, NDSTSF
C
      INTEGER ELMSUR(NELM), ENLBAS(NELM+1), NELBAS(NNOD+1),
     :        ENLIST(SZENLS), NELIST(SZNELS), CNCTSF(NDSTSF),
     :        EELIST(NELM*4)
C
      REAL NODX(MAXNOD), NODY(MAXNOD), NODZ(MAXNOD)
C
      LOGICAL GEOM3D
C
      INTEGER CONFAC(3), EDGLST(3), TMPLST(3), IBAS, IELM, NLN, ILN,
     :        CONELM(3), JBAS, IPTR, INEWFC, INOD1, INOD2, INOD3,
     :        NODS(2), NXT, LST, INEWED, IFLAG, JELM
C
      IS3DMS = GEOM3D
      IS2DFR = GEOM3D
      ADPTNG = .FALSE.
      debug = .false.
C
      NNEWFC = 0
c      NNEWED = 0
c      NNEWEL = 0
      NFACES = 0
      NELEMS = 0
      NEDGES = 0
C
      STTBIG = 1
      ENDBIG = -1
      STTFRT = 1
      ENDFRT = -1
C
      EDGLST(1) = 0
      EDGLST(2) = 0
      EDGLST(3) = 0
C
      PRINT*,'Setting up initial list pointers...'
C
      CALL INITLS( EMTBIG, BIGLST, SIZBIG, NXTBIG, LSTBIG, MAXBIG )
      CALL INITLS( EMTFRT, FRTLST, SIZFRT, NXTFRT, LSTFRT, MAXFRT )
      CALL INIRLS( EMTNOD, NODLST, SIZNOD, NXTNOD, LSTNOD, MXNODS )
C
      DO ILN = 1, NDSTSF
         CNCTSF(ILN) = 0
      END DO
C
      STTNOD = 1
      STFRND = 1
C
      PRINT*,'Creating node list...'
C
      DO ILN = 1, NNOD
         NODLST(1,ILN) = NODX(ILN)
         NODLST(2,ILN) = NODY(ILN)
         NODLST(3,ILN) = NODZ(ILN)
         NODLST(6,ILN) = 3.0
         CALL FIXEIG( NODLST(1,ILN), NODLST(7,ILN), ANISOT )
C         NODLST( 7,ILN) = 44.44444
C         NODLST(11,ILN) = 44.44444
C         NODLST(15,ILN) = 400.0
      END DO
C
      NUMNDS = NNOD
      NFRTND = NNOD
      ENDNOD = NNOD
      EMTNOD = NNOD+1
C
      IF( .NOT. GEOM3D ) THEN
C
         DO IELM = 1, NELM
C
            IBAS = ENLBAS(IELM)
            NLN  = ENLBAS(IELM+1) - IBAS
C
            IF( NLN .NE. 2 ) THEN
               IERR = -12
               PRINT*,'*** MKFRNT: INCORRECT NODES IN 1D SURFACES'
               RETURN
            END IF
C
            DO ILN = 1, NLN
C
               INOD1 = ENLIST(IBAS+ILN)
               JBAS  = NELBAS(INOD1)
               IF( NELIST(JBAS+1) .NE. IELM ) THEN
                  CONELM(ILN) = NELIST(JBAS+1)
               ELSE IF( NELIST(JBAS+2) .NE. IELM ) THEN
                  CONELM(ILN) = NELIST(JBAS+2)
               ELSE
                  IERR = -10
                  PRINT*,'*** MKFRNT: CONNECTION ERROR IN 1D NELIST'
                  RETURN
               END IF
C
               IF( CONELM(ILN) .LE. 0 ) THEN
                  IERR = -11
                  PRINT*,'*** MKFRNT: INVALID 1D SURFACE'
                  RETURN
               END IF
C
            END DO
C
            CONELM(3) = 0
C
            INEWFC = -1
C
            CALL ADDFAC( CONELM, .FALSE., EDGLST, -ELMSUR(IELM),
     :                   IPTR, INEWFC )
C
            IF( IERR .NE. 0 ) RETURN
C
            IF( IPTR .NE. IELM ) THEN
               IERR = -20
               PRINT*,'*** MKFRNT: INCONSISTENCY CREATING FRTLST'
               RETURN
            END IF
C
         END DO
C
      ELSE
C
         PRINT*,'Creating face and edge lists...'
C
         DO IELM = 1, NELM
C
            if( debug ) print*,'>>> Adding face ',ielm
C
            IBAS = ENLBAS(IELM)
            NLN  = ENLBAS(IELM+1) - IBAS
C
            IF( NLN .NE. 3 ) THEN
               IERR = -12
               PRINT*,'*** MKFRNT: INCORRECT NODES IN 2D SURFACES'
               RETURN
            END IF
C
            INOD1 = ENLIST(IBAS+1)
            INOD2 = ENLIST(IBAS+2)
            INOD3 = ENLIST(IBAS+3)
C
c            CONELM(1) = IELM
c            CONELM(2) = IELM
c            CONELM(3) = IELM
C
c            CALL FNDEDG( INOD1, INOD2, ENLBAS, ENLIST,
c     :                   NELBAS, NELIST, NNOD, NELM,
c     :                   SZENLS, SZNELS, CONELM(1) )
C
c            CALL FNDEDG( INOD2, INOD3, ENLBAS, ENLIST,
c     :                   NELBAS, NELIST, NNOD, NELM,
c     :                   SZENLS, SZNELS, CONELM(2) )
C
c            CALL FNDEDG( INOD3, INOD1, ENLBAS, ENLIST,
c     :                   NELBAS, NELIST, NNOD, NELM,
c     :                   SZENLS, SZNELS, CONELM(3) )
C
c            IF( CONELM(1) .NE. EELIST(IELM*3-2) ) THEN
c               PRINT*,'INCONSISTENCY 1: ',(EELIST(IELM*3-ILN),
c     :                                                  ILN=2,0,-1)
c            END IF
C
c            IF( CONELM(2) .NE. EELIST(IELM*3-1) ) THEN
c               PRINT*,'INCONSISTENCY 2: ',(EELIST(IELM*3-ILN),
c     :                                                  ILN=2,0,-1)
c            END IF
C
c            IF( CONELM(3) .NE. EELIST(IELM*3) ) THEN
c               PRINT*,'INCONSISTENCY 3: ',(EELIST(IELM*3-ILN),
c     :                                                  ILN=2,0,-1)
c            END IF
C
            CONELM(1) = EELIST(IELM*3-2)
            CONELM(2) = EELIST(IELM*3-1)
            CONELM(3) = EELIST(IELM*3)
            if( debug ) print*,'    connected faces: ',CONELM
C
            IF( CONELM(1) .GT. IELM ) THEN
               NODS(1) = INOD1
               NODS(2) = INOD2
c               print*,'   Nodes ',inod1,inod2
               INEWED = -1
               CALL ADDEDG( BIGLST, NODLST, NODS, IPTR, INEWED )
               IF( IERR .NE. 0 ) RETURN
               CALL SETFLG( BIGLST, NODLST, IPTR, 2, .TRUE. )
               EDGLST(1) = IPTR
               if( debug ) print*,'   Added edge ',iptr,inod1,inod2
            ELSE IF( CONELM(1) .NE. 0 ) THEN
               IPTR = CONELM(1)
c               print*,'   Nodes ',inod1,inod2,iptr
               CALL FACINF( IPTR, CONFAC, TMPLST,
     :                      JELM, INEWFC, IFLAG, NXT, LST )
c               print*,'   Face info: ',tmplst
               CALL EDGINF( BIGLST, NODLST, TMPLST(1), NODS,
     :                      INEWED, IFLAG, NXT, LST )
               IF( ( NODS(1).EQ.INOD1 .AND. NODS(2).EQ.INOD2 )
     :        .OR. ( NODS(1).EQ.INOD2 .AND. NODS(2).EQ.INOD1 ) ) THEN
                  EDGLST(1) = TMPLST(1)
c                  print*,'   Nodes ',nods(1),nods(2)
                  if( debug ) print*,'   Edge info: ',tmplst(1),nods
                  GOTO 100
               END IF
               CALL EDGINF( BIGLST, NODLST, TMPLST(2), NODS,
     :                      INEWED, IFLAG, NXT, LST )
               IF( ( NODS(1).EQ.INOD1 .AND. NODS(2).EQ.INOD2 )
     :        .OR. ( NODS(1).EQ.INOD2 .AND. NODS(2).EQ.INOD1 ) ) THEN
                  EDGLST(1) = TMPLST(2)
c                  print*,'   Nodes ',nods(1),nods(2)
                  if( debug ) print*,'   Edge info: ',tmplst(2),nods
                  GOTO 100
               END IF
               CALL EDGINF( BIGLST, NODLST, TMPLST(3), NODS,
     :                      INEWED, IFLAG, NXT, LST )
               IF( ( NODS(1).EQ.INOD1 .AND. NODS(2).EQ.INOD2 )
     :        .OR. ( NODS(1).EQ.INOD2 .AND. NODS(2).EQ.INOD1 ) ) THEN
                  EDGLST(1) = TMPLST(3)
c                  print*,'   Nodes ',nods(1),nods(2)
                  if( debug ) print*,'   Edge info: ',tmplst(3),nods
                  GOTO 100
               END IF
               IERR = -21
               PRINT*,'*** MKFRNT: INCONSISTENCY IN LISTS 1'
               RETURN
            ELSE
               IERR = -22
               PRINT*,'*** MKFRNT: FOUND NON-CONNECTED FACE'
               RETURN
            END IF
C
 100        IF( CONELM(2) .GT. IELM ) THEN
               NODS(1) = INOD2
               NODS(2) = INOD3
c               print*,'   Nodes ',inod2,inod3
               INEWED = -1
               CALL ADDEDG( BIGLST, NODLST, NODS, IPTR, INEWED )
               IF( IERR .NE. 0 ) RETURN
               CALL SETFLG( BIGLST, NODLST, IPTR, 2, .TRUE. )
               EDGLST(2) = IPTR
               if( debug ) print*,'   Added edge ',iptr,inod2,inod3
            ELSE IF( CONELM(2) .NE. 0 ) THEN
               IPTR = CONELM(2)
c               print*,'   Nodes ',inod2,inod3,iptr
               CALL FACINF( IPTR, CONFAC, TMPLST,
     :                      JELM, INEWFC, IFLAG, NXT, LST )
               CALL EDGINF( BIGLST, NODLST, TMPLST(1), NODS,
     :                      INEWED, IFLAG, NXT, LST )
c               if( debug ) print*,'   Edge info: ',tmplst(1),nods
               IF( ( NODS(1).EQ.INOD3 .AND. NODS(2).EQ.INOD2 )
     :        .OR. ( NODS(1).EQ.INOD2 .AND. NODS(2).EQ.INOD3 ) ) THEN
                  EDGLST(2) = TMPLST(1)
                  if( debug ) print*,'   Edge info: ',tmplst(1),nods
                  GOTO 200
               END IF
               CALL EDGINF( BIGLST, NODLST, TMPLST(2), NODS,
     :                      INEWED, IFLAG, NXT, LST )
               IF( ( NODS(1).EQ.INOD3 .AND. NODS(2).EQ.INOD2 )
     :        .OR. ( NODS(1).EQ.INOD2 .AND. NODS(2).EQ.INOD3 ) ) THEN
                  EDGLST(2) = TMPLST(2)
c                  print*,'   Nodes ',nods(1),nods(2)
                  if( debug ) print*,'   Edge info: ',tmplst(2),nods
                  GOTO 200
               END IF
               CALL EDGINF( BIGLST, NODLST, TMPLST(3), NODS,
     :                      INEWED, IFLAG, NXT, LST )
               IF( ( NODS(1).EQ.INOD3 .AND. NODS(2).EQ.INOD2 )
     :        .OR. ( NODS(1).EQ.INOD2 .AND. NODS(2).EQ.INOD3 ) ) THEN
                  EDGLST(2) = TMPLST(3)
c                  print*,'   Nodes ',nods(1),nods(2)
                  if( debug ) print*,'   Edge info: ',tmplst(3),nods
                  GOTO 200
               END IF
               IERR = -21
               PRINT*,'*** MKFRNT: INCONSISTENCY IN LISTS 2'
               RETURN
            ELSE
               IERR = -22
               PRINT*,'*** MKFRNT: FOUND NON-CONNECTED FACE'
               RETURN
            END IF
C
 200        IF( CONELM(3) .GT. IELM ) THEN
               NODS(1) = INOD3
               NODS(2) = INOD1
c               print*,'   Nodes ',inod3,inod1
               INEWED = -1
               CALL ADDEDG( BIGLST, NODLST, NODS, IPTR, INEWED )
               IF( IERR .NE. 0 ) RETURN
               CALL SETFLG( BIGLST, NODLST, IPTR, 2, .TRUE. )
               EDGLST(3) = IPTR
               if( debug ) print*,'   Added edge ',iptr,inod3,inod1
            ELSE IF( CONELM(3) .NE. 0 ) THEN
               IPTR = CONELM(3)
c               print*,'   Nodes ',inod3,inod1,iptr
               CALL FACINF( IPTR, CONFAC, TMPLST,
     :                      JELM, INEWFC, IFLAG, NXT, LST )
               CALL EDGINF( BIGLST, NODLST, TMPLST(1), NODS,
     :                      INEWED, IFLAG, NXT, LST )
               IF( ( NODS(1).EQ.INOD1 .AND. NODS(2).EQ.INOD3 )
     :        .OR. ( NODS(1).EQ.INOD3 .AND. NODS(2).EQ.INOD1 ) ) THEN
                  EDGLST(3) = TMPLST(1)
c                  print*,'   Nodes ',nods(1),nods(2)
                  if( debug ) print*,'   Edge info: ',tmplst(1),nods
                  GOTO 300
               END IF
               CALL EDGINF( BIGLST, NODLST, TMPLST(2), NODS,
     :                      INEWED, IFLAG, NXT, LST )
               IF( ( NODS(1).EQ.INOD1 .AND. NODS(2).EQ.INOD3 )
     :        .OR. ( NODS(1).EQ.INOD3 .AND. NODS(2).EQ.INOD1 ) ) THEN
                  EDGLST(3) = TMPLST(2)
c                  print*,'   Nodes ',nods(1),nods(2)
                  if( debug ) print*,'   Edge info: ',tmplst(2),nods
                  GOTO 300
               END IF
               CALL EDGINF( BIGLST, NODLST, TMPLST(3), NODS,
     :                      INEWED, IFLAG, NXT, LST )
               IF( ( NODS(1).EQ.INOD1 .AND. NODS(2).EQ.INOD3 )
     :        .OR. ( NODS(1).EQ.INOD3 .AND. NODS(2).EQ.INOD1 ) ) THEN
                  EDGLST(3) = TMPLST(3)
c                  print*,'   Nodes ',nods(1),nods(2)
                  if( debug ) print*,'   Edge info: ',tmplst(3),nods
                  GOTO 300
               END IF
               IERR = -21
               PRINT*,'*** MKFRNT: INCONSISTENCY IN LISTS 3'
               RETURN
            ELSE
               IERR = -22
               PRINT*,'*** MKFRNT: FOUND NON-CONNECTED FACE'
               RETURN
            END IF
C
 300        INEWFC = -1
C
            CALL ADDFAC( CONELM, .FALSE., EDGLST, -ELMSUR(IELM),
     :                   IPTR, INEWFC )
C
            IF( IERR .NE. 0 ) RETURN
C
            IF( IPTR .NE. IELM ) THEN
               IERR = -20
               PRINT*,'*** MKFRNT: INCONSISTENCY CREATING FRTLST'
               RETURN
            END IF
C
         END DO
C
      END IF
C
      CALL STFRGM( BIGLST, NODLST )
C
      print*,' '
      print*,'++++++++++ Finished making front ++++++++++++'
      print*,'STTFRT,EMTFRT,ENDFRT: ',STTFRT,EMTFRT,ENDFRT
      print*,'NFACES,NNEWFC: ',NFACES,NNEWFC
      print*,' '
c      debug = .false.
      RETURN
      END
C
