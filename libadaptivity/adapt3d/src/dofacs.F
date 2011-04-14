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
      SUBROUTINE DOFACS( BIGLST, NODLST, IFAC1, CONFC1, EDGLS1,
     :                   IELM1, NXT1, E12,
     :                   THSDOT, THSREL, MINDOT, MINREL, DOIT, IREG )
C-----------------------------------------------------------------------
C
C - This subroutine takes two faces (IFAC1 and CONFC1(E12)) and sees what
C - happens if their opposing nodes are joined by an edge, to form two new faces
C - that, along with the original two, form a tetrahedron.
C - If they are already joined (i.e. there is a third face that is common to the
C - other two) then it sees what happens if a new face is added such that a
C - tetrahedron is formed from this and the other three faces.
C
C - The checks that need to be performed on new faces/edges are:
C -   i) New edge must lie inside domain (ANGFAC and MINANG check this)
C -  ii) New faces must not be punctured by another edge (CHKFAC checks this)
C - iii) New edge must not be too close (in a way defined by the size of the
C -      existing faces/edges) to another edge (CHKFAC also checks this)
C
C - MINDOT and MINREL give the minimum allowed cosine and length aspect ratios
C - to form a tetrahedron (i.e. none will be formed if the angle between a
C - new edge and an existing face is too small or negative, or if the length
C - of a new edge or the distance to another edge is too small compared to the
C - lengths of the existing edges). THSDOT and THSREL return the cos of the
C - minimum edge angle and the minimum relative length aspect ratio.
C
C - DOIT is a flag which, if true, makes this subroutine actually add a new
C - acceptable tetrahedron. If false then no tetrahedron will be added, but
C - THSDOT and THSREL are returned as if a new tetrahedron had been formed.
C - This allows a decision to be made as to which tets are best to add first.
C
C - IREG contains the element 'region' number that should be given to a tet
C - that is added.
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
      INTEGER CONFC1(3), EDGLS1(3), IFAC1, E12, IELM1, NXT1, IREG
C
      REAL MINDOT, MINREL, THSDOT, THSREL
C
      LOGICAL DOIT
C
      INTEGER IFAC2, CONFC3(3), EDGLS3(3), IELM2, IELM3, IFAC3,
     :        CONFC2(3), EDGLS2(3), INEW, IFLAG, NXT2, NXT3, LST,
     :        EDG1, EDG2, EDG3, EDG4, F3, EDG13, EDG23, IEDG, IELM,
     :        NODS12(2), NODS13(2), NODS23(2), NODE1, NODE2,
     :        N12, N13, N14, N23, N24, N34, NODS(2), IPTR, CONELM(4),
     :        CONFAC(6), EDGLST(6), IFC1, IFC2, N1, N2, N3, N4, N123,
     :        FAC1, FAC2, FAC3, FAC4, EDG12
C
      REAL DOT1, DOT2, RAD1, RAD2, DST1, DST2, DST12, DST13, DST23,
     :     LENEDG, DOTN, DXE, DYE, DZE, D1, D2, D3
C
      LOGICAL REVERS, GTFRFL, GETFLG
c
      if( debug ) then
       print*,' '
       print*,'Into DOFACS... ',ierr,doit
       print*,'    IFAC1:',ifac1,ielm1
       print*,'   CONFC1:',confc1
       print*,'   EDGLS1:',edgls1
       PRINT*,'     NXT1:',nxt1
       PRINT*,'    EDG12:',e12
      end if
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( E12 .LE. 0 .OR. E12 .GT. 3 ) THEN
         IERR = -40
         PRINT*,'*** DOFACS: INVALID EDGE VALUE'
         RETURN
      END IF
C
      IF( IFAC1 .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** DOFACS: INVALID FACE POINTER'
         RETURN
      END IF
C
      EDG12 = E12
      IFAC2 = CONFC1(EDG12)
C
c      if( doit )
c      print*,'    Checking faces ',ifac1,ifac2
c
      CALL FACINF( IFAC2, CONFC2, EDGLS2, IELM2, INEW, IFLAG,
     :             NXT2, LST )
C
      IF( IERR .NE. 0 ) RETURN
c
      if( debug ) then
       print*,'    IFAC2:',ifac2,ielm2
       print*,'   CONFC2:',confc2
       print*,'   EDGLS2:',edgls2
      end if
C
      IFAC3 = 0
      EDG13 = MOD(EDG12,3) + 1
      EDG4  = MOD(EDG13,3) + 1
C
      IFAC3 = CONFC1(EDG13)
      EDG1  = EDG4
      IF( IFAC3 .EQ. CONFC2(1) ) THEN
         EDG23 = 1
         GOTO 100
      ELSE IF( IFAC3 .EQ. CONFC2(2) ) THEN
         EDG23 = 2
         GOTO 100
      ELSE IF( IFAC3 .EQ. CONFC2(3) ) THEN
         EDG23 = 3
         GOTO 100
      END IF
C
      IFAC3 = CONFC1(EDG4)
      EDG1  = EDG13
      EDG13 = EDG4
      IF( IFAC3 .EQ. CONFC2(1) ) THEN
         EDG23 = 1
         GOTO 100
      ELSE IF( IFAC3 .EQ. CONFC2(2) ) THEN
         EDG23 = 2
         GOTO 100
      ELSE IF( IFAC3 .EQ. CONFC2(3) ) THEN
         EDG23 = 3
         GOTO 100
      END IF
C
      IF( EDG12 .EQ. 1 ) THEN
         EDG1 = EDGLS1(2)
         EDG2 = EDGLS1(3)
         FAC1 = CONFC1(2)
         FAC2 = CONFC1(3)
      ELSE IF( EDG12 .EQ. 2 ) THEN
         EDG1 = EDGLS1(3)
         EDG2 = EDGLS1(1)
         FAC1 = CONFC1(3)
         FAC2 = CONFC1(1)
      ELSE
         EDG1 = EDGLS1(1)
         EDG2 = EDGLS1(2)
         FAC1 = CONFC1(1)
         FAC2 = CONFC1(2)
      END IF
C
      EDG12 = EDGLS1(EDG12)
C
      IF( EDG12 .EQ. EDGLS2(1) ) THEN
         EDG3 = EDGLS2(2)
         EDG4 = EDGLS2(3)
         FAC3 = CONFC2(2)
         FAC4 = CONFC2(3)
      ELSE IF( EDG12 .EQ. EDGLS2(2) ) THEN
         EDG3 = EDGLS2(3)
         EDG4 = EDGLS2(1)
         FAC3 = CONFC2(3)
         FAC4 = CONFC2(1)
      ELSE IF( EDG12 .EQ. EDGLS2(3) ) THEN
         EDG3 = EDGLS2(1)
         EDG4 = EDGLS2(2)
         FAC3 = CONFC2(1)
         FAC4 = CONFC2(2)
      ELSE
         IERR = -41
         PRINT*,'*** DOFACS: FACE2 DOES NOT HAVE THE COMMON EDGE'
         RETURN
      END IF
C
C      CALL CMNNOD( EDG4, EDG1, N14, N4, N1 )
C      CALL CMNNOD( EDG2, EDG3, N23, N2, N3 )
C
      CALL EDG4ND( BIGLST, NODLST, EDG1, EDG2, EDG3, EDG4,
     :             N12, N34, N14, N23 )
C
      DXE = NODLST(1,N34) - NODLST(1,N12)
      DYE = NODLST(2,N34) - NODLST(2,N12)
      DZE = NODLST(3,N34) - NODLST(3,N12)
C
c      print*,'>>> Wondering about joining nodes ',n12,n34
C
      DOT1 = DXE*DXE + DYE*DYE + DZE*DZE
      DOT1 = 1.0/SQRT(DOT1)
C
      DXE = DXE*DOT1
      DYE = DYE*DOT1
      DZE = DZE*DOT1
C
      CALL ANGFAC( BIGLST, NODLST, N12, N23, N14, DXE, DYE, DZE,
     :             THSDOT )
c      print*,'     cos of face ',ifac1,' angle: ',thsdot
C
c      if( ielm1 .lt. 0 .and. ielm2 .lt. 0 ) then
c      if( getflg(BIGLST,NODLST,edg12,2) .and. abs(thsdot).lt.1e-5 ) then
c         call setflg(BIGLST,NODLST,edg12,2,.false.)
c         print*,'--- Info: made edge non-geometry ',edg12,
c     :          biglst(4,edg12)
c         stop
c      else if( .not. getflg(BIGLST,NODLST,edg12,3) .and.
c     :         .not. getflg(BIGLST,NODLST,edg12,2)
c     :                               .and. abs(thsdot).ge.1e-5 ) then
c         call setflg(BIGLST,NODLST,edg12,2,.true.)
c         print*,'--- Info: made edg geometry ',edg12,biglst(4,edg12)
c         print*,'Nodes: ',biglst(1,edg12),biglst(2,edg12),thsdot
c         stop
c      end if
c      end if
c
      THSDOT = -THSDOT
      IF( ( THSDOT .LT. MINDOT ) .AND. DOIT ) THEN
         THSREL = MINREL
         GOTO 900
      END IF
C
      CALL ANGFAC( BIGLST, NODLST, N34, N23, N14, DXE, DYE, DZE, DOTN )
c      print*,'     cos of face ',ifac2,' angle: ',dotn
C
      IF( ( DOTN .GT. -MINDOT ) .AND. ( DOIT .OR. DOTN .GE. 0.0 ) ) THEN
         THSDOT = -DOTN
         THSREL = MINREL
         GOTO 900
      END IF
C
      THSDOT = MIN( THSDOT, -DOTN )
C
      CALL MINANG( BIGLST, NODLST, DXE, DYE, DZE, N12, IFAC1, DOTN )
c      print*,'     cos of node ',n12,' min ang: ',dotn
C
      IF( ( DOTN .GT. -MINDOT ) .AND. ( DOIT .OR. DOTN .GE. 0.0 ) ) THEN
         THSDOT = -DOTN
         THSREL = MINREL
         GOTO 900
      END IF
C
      THSDOT = MIN( THSDOT, -DOTN )
C
      CALL MINANG( BIGLST, NODLST, -DXE, -DYE, -DZE, N34, IFAC2, DOTN )
c      print*,'     cos of node ',n34,' min ang: ',dotn
C
      IF( ( DOTN .GT. -MINDOT ) .AND. DOIT ) THEN
         THSDOT = -DOTN
         THSREL = MINREL
         GOTO 900
      END IF
C
      THSDOT = MIN( THSDOT, -DOTN )
C
      DST12 = LENEDG( BIGLST, NODLST, EDG12 )
C
      D1 = LENEDG( BIGLST, NODLST, EDG4 )
      D2 = LENEDG( BIGLST, NODLST, EDG1 )
C
      DST1 = MAX( D1, D2 )
      DST1 = MAX( DST12, DST1 )
C
      D3 = MIN( D1, D2 )
      D3 = MIN( DST12, D3 )
C
      RAD1 = MIN( MINREL*DST1, D3 )
C
      DST2 = DST1
C
c      print*,'Calling CHKFAC 1:',DST1,DST2,RAD1,RAD2
      CALL CHKFAC( BIGLST, NODLST, N14, N34, N12,
     :             IFAC1, IFAC2, 0, FAC4, FAC1, 0,
     :             RAD1, NODE1, DST1 )
C
c      print*,'1) DST1,DST2,RAD1,RAD2:',DST1,DST2,RAD1,RAD2
      THSREL = DST1/DST2
      IF( ( DST1 .LT. RAD1 ) .AND. DOIT ) GOTO 900
C
      D1 = LENEDG( BIGLST, NODLST, EDG3 )
      D2 = LENEDG( BIGLST, NODLST, EDG2 )
C
      DST2 = MAX( D1, D2 )
      DST2 = MAX( DST12, DST2 )
C
      D3 = MIN( D1, D2 )
      D3 = MIN( DST12, D3 )
C
      RAD2 = MIN( DST2*MINREL, D3 )
C
      DST1 = DST2
C
c      print*,'Calling CHKFAC 2:',DST1,DST2,RAD1,RAD2
      CALL CHKFAC( BIGLST, NODLST, N23, N12, N34,
     :             IFAC1, IFAC2, 0, FAC3, FAC2, 0,
     :             RAD2, NODE2, DST2 )
c      print*,'2) DST1,DST2,RAD1,RAD2:',DST1,DST2,RAD1,RAD2
C
      THSREL = MIN( DST2/DST1, THSREL )
C
      IF( ( DST2 .LT. RAD2 ) .AND. DOIT ) GOTO 900
C
      IF( .NOT. DOIT ) GOTO 900
C
C - add a new edge with nodes N14 and N23
C
      NODS(1) = N12
      NODS(2) = N34
C
      CALL ADDEDG( BIGLST, NODLST, NODS, IEDG, INEW )
C
      CALL SETFLG( BIGLST, NODLST, IEDG, 3, .TRUE. )
C
C - remove the two faces from FRTLST
C
C      IF( NXT1 .EQ. IFAC2 ) NXT1 = NXT2
C
      IF( IFAC1 .EQ. ENDFRT ) THEN
         NXT1 = IFAC1
      END IF
C
      CALL RMVFAC( IFAC2 )
      CALL RMVFAC( IFAC1 )
C
      IF( IERR .NE. 0 ) RETURN
C
      IFAC1 = 0
C
C - add new element into BIGLST
C
      CONELM(1) = IELM1
      CONELM(2) = IELM2
      CONELM(3) = 0
      CONELM(4) = 0
C
      EDGLST(1) = EDG4
      EDGLST(2) = EDG1
      EDGLST(3) = IEDG
      EDGLST(4) = EDG2
      EDGLST(5) = EDG3
      EDGLST(6) = EDG12
C
      CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, IELM, INEW, IREG )
C
      if( ierr .eq. 12345 ) then
         print*,'--- got bad oppedg in dofacs 1'
         ierr = 0
      endif
c
      IF( IERR .NE. 0 ) RETURN
C
      IF( IELM .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** DOFACS: FAILED TO CREATE NEW ELEMENT'
         RETURN
      END IF
C
C - add two new faces into FRTLST and connect them together
C
      CONFAC(1) = FAC4
      CONFAC(2) = FAC1
      CONFAC(3) = 0
      CONFAC(4) = FAC2
      CONFAC(5) = FAC3
C
      CALL ADDFAC( CONFAC, .TRUE., EDGLST, IELM, IFC1, INEW )
C
      IF( IFC1 .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** DOFACS: FAILED TO CREATE NEW FACE 1'
         RETURN
      END IF
C
      CALL ADDFAC( CONFAC(3), .TRUE., EDGLST(3), IELM, IFC2, INEW )
C
      IF( IFC2 .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** DOFACS: FAILED TO CREATE NEW FACE 2'
         RETURN
      END IF
C
      CALL CNCTFC( IFC1, IFC2 )
c
      print*,'>>>>>>> NEW EDGE, ELEMENT & FACES CREATED HERE <<<<<<<'
c      ierr = 1
C
      GOTO 900
C
 100  CONTINUE
C
c      if( doit )
c      print*,'     found common connected face ',ifac3
C
      THSREL = MINREL/2
      THSDOT = MINDOT/2
C
      IF( .NOT. GTFRFL(IFAC3,2) ) GOTO 900
C
      CALL FACINF( IFAC3, CONFC3, EDGLS3, IELM3, INEW, IFLAG,
     :             NXT3, LST )
C
      IF( IERR .NE. 0 ) RETURN
C
C - EDGLS1(EDG12) is the common edge between faces 1 and 2
C - EDGLS1(EDG13) is the common edge between faces 1 and 3
C - EDGLS1(EDG1)  is the remaining edge of face 1
C - EDGLS2(EDG23) is the common edge between faces 2 and 3
C
C - figure out which way round the edges are going,
C - and then perform some consistency checks
C
      IF( MOD(EDG1,3) .EQ. EDG13-1 ) THEN
         REVERS = .TRUE.
         IF( MOD(EDG13,3) .NE. EDG12-1 ) IERR = -45
         IF( MOD(EDG12,3) .NE. EDG1-1  ) IERR = -45
         IF( IERR .NE. 0 ) THEN
            PRINT*,'*** DOFACS: INCONSISTENCY 2 IN EDGES OF FACE1'
c            print*,'edges of face 1: ',edgls1
c            print*,'edges of face 2: ',edgls2
c            print*,'edges of face 3: ',edgls3
c            print*,'edg1,edg12,edg13,edg23: ',edg1,edg12,edg13,edg23
            RETURN
         END IF
      ELSE
         REVERS = .FALSE.
         IF( MOD(EDG13,3) .NE. EDG1-1  ) IERR = -46
         IF( MOD(EDG12,3) .NE. EDG13-1 ) IERR = -46
         IF( MOD(EDG1,3)  .NE. EDG12-1 ) IERR = -46
         IF( IERR .NE. 0 ) THEN
            PRINT*,'*** DOFACS: INCONSISTENCY 1 IN EDGES OF FACE1'
c            print*,'edges of face 1: ',edgls1
c            print*,'edges of face 2: ',edgls2
c            print*,'edges of face 3: ',edgls3
c            print*,'edg1,edg12,edg13,edg23: ',edg1,edg12,edg13,edg23
            RETURN
         END IF
      END IF
c      print*,'   got revers: ',revers
C
C - need to find which edge of FACE2 is EDG2 - we know it's not EDG23
C
      EDG2 = MOD(EDG23,3) + 1
      IF( CONFC2(EDG2).EQ.IFAC1 ) THEN
         EDG2 = MOD(EDG2,3) + 1
      END IF
C
C - check that worked OK...
C
      IF( CONFC2(EDG2).EQ.IFAC1 .OR. CONFC2(EDG2).EQ.IFAC3 ) THEN
         IERR = -47
         PRINT*,'*** DOFACS: ERROR IN EDGES OF FACE 2'
         print*,'edges of face 1: ',edgls1
         print*,'edges of face 2: ',edgls2
         print*,'edges of face 3: ',edgls3
         print*,'edg1,edg2: ',edg1,edg2
         print*,'edg12,edg13,edg23: ',edg12,edg13,edg23
         RETURN
      END IF
C
C - now find which edge of FACE3 is EDG3
C
      IF( EDGLS3(1) .EQ. EDGLS1(EDG13) ) THEN
         EDG3 = 3
      ELSE IF( EDGLS3(1) .EQ. EDGLS2(EDG23) ) THEN
         EDG3 = 2
      ELSE
         EDG3 = 1
      END IF
C
C - if the edges go anticlockwise for FACE1 then we must also
C - swap round the local edge value for EDG3
C
      IF( REVERS ) THEN
         IF( EDG3 .EQ. 3 ) THEN
            EDG3 = 2
         ELSE IF( EDG3 .EQ. 2 ) THEN
            EDG3 = 3
         END IF
      END IF
c
c      print*,'  edges of face 1: ',edgls1
c      print*,'  edges of face 2: ',edgls2
c      print*,'  edges of face 3: ',edgls3
c      print*,'   edg1,edg2,edg3: ',edg1,edg2,edg3
c      print*,'edg12,edg13,edg23: ',edg12,edg13,edg23
C
C - find connected faces for EDG1, EDG2 and EDG3
C
c      print*,'  faces of face 1: ',confc1
c      print*,'  faces of face 2: ',confc2
c      print*,'  faces of face 3: ',confc3
      FAC1 = CONFC1(EDG1)
      FAC2 = CONFC2(EDG2)
      FAC3 = CONFC3(EDG3)
C
C - check that worked OK
C
      IF( FAC3 .eq. IFAC1 .or. FAC3 .eq. IFAC2 ) THEN
         IERR = -48
         PRINT*,'*** DOFACS: ERROR IN EDGES OF FACE 3'
         print*,'  edges of face 1: ',edgls1
         print*,'  edges of face 2: ',edgls2
         print*,'  edges of face 3: ',edgls3
         print*,'   edg1,edg2,edg3: ',edg1,edg2,edg3
         print*,'edg12,edg13,edg23: ',edg12,edg13,edg23
         RETURN
      END IF
C
c      print*,'  got faces: ',fac1,fac2,fac3
C
C - now that all edges have been found, we can replace the local
C - edge values with their BIGLST pointers
C
      EDG1  = EDGLS1(EDG1)
      EDG2  = EDGLS2(EDG2)
      EDG3  = EDGLS3(EDG3)
      EDG12 = EDGLS1(EDG12)
      EDG13 = EDGLS1(EDG13)
      EDG23 = EDGLS2(EDG23)
c      print*,'   Got edges...'
c      print*,'          edg1,edg2,edg3: ',edg1,edg2,edg3
c      print*,'       edg12,edg13,edg23: ',edg12,edg13,edg23
C
      CALL EDG3ND( BIGLST, NODLST, EDG12, EDG13, EDG23, N123,
     :             N12, N13, N23 )
C
      DXE = NODLST(1,N23) - NODLST(1,N13)
      DYE = NODLST(2,N23) - NODLST(2,N13)
      DZE = NODLST(3,N23) - NODLST(3,N13)
C
      DOT1 = DXE*DXE + DYE*DYE + DZE*DZE
      DOT1 = 1.0/SQRT(DOT1)
C
      DXE = DXE*DOT1
      DYE = DYE*DOT1
      DZE = DZE*DOT1
C
      THSREL = MINREL
C
      IF( REVERS ) THEN
         CALL ANGFAC( BIGLST, NODLST, N13, N123, N12, DXE, DYE, DZE,
     :                THSDOT )
         if( ielm1 .lt. 0 .and. ielm2 .lt. 0 ) then
         if( getflg(BIGLST,NODLST,edg12,2) .and.
     :       abs(thsdot).lt.1e-5 ) then
            call setflg(BIGLST,NODLST,edg12,2,.false.)
            print*,'--- Info: made edge non-geometry ',edg12,
     :             biglst(4,edg12)
            stop
         else if( .not. getflg(BIGLST,NODLST,edg12,3) .and.
     :            .not. getflg(BIGLST,NODLST,edg12,2)
     :            .and. abs(thsdot).ge.1e-5 ) then
            call setflg(BIGLST,NODLST,edg12,2,.true.)
            print*,'--- Info: made edg geometry ',edg12,biglst(4,edg12)
            print*,'Nodes: ',biglst(1,edg12),biglst(2,edg12),thsdot
            stop
         end if
         end if
         THSDOT = -THSDOT
         IF( ( THSDOT .LT. MINDOT ) .AND. DOIT ) GOTO 900
         CALL ANGFAC( BIGLST, NODLST, N23, N123, N12, DXE, DYE, DZE,
     :                DOTN )
         IF( ( DOTN .GT. -MINDOT ) .AND. DOIT ) THEN
            THSDOT = -DOTN
            GOTO 900
         END IF
      ELSE
         CALL ANGFAC( BIGLST, NODLST, N13, N12, N123, DXE, DYE, DZE,
     :                THSDOT )
c         print*,'   1st dot: ',thsdot
         if( ielm1 .lt. 0 .and. ielm2 .lt. 0 ) then
         if( getflg(BIGLST,NODLST,edg12,2) .and.
     :       abs(thsdot).lt.1e-5 ) then
            call setflg(BIGLST,NODLST,edg12,2,.false.)
            print*,'--- Info: made edge non-geometry ',edg12,
     :             biglst(4,edg12)
            stop
         else if( .not. getflg(BIGLST,NODLST,edg12,3) .and.
     :            .not. getflg(BIGLST,NODLST,edg12,2)
     :            .and. abs(thsdot).ge.1e-5 ) then
            call setflg(BIGLST,NODLST,edg12,2,.true.)
            print*,'--- Info: made edg geometry ',edg12,biglst(4,edg12)
            print*,'Nodes: ',biglst(1,edg12),biglst(2,edg12),thsdot
            stop
         end if
         end if
         THSDOT = -THSDOT
         IF( ( THSDOT .LT. MINDOT ) .AND. DOIT ) GOTO 900
         CALL ANGFAC( BIGLST, NODLST, N23, N12, N123, DXE, DYE, DZE,
     :                DOTN )
c         print*,'   2nd dot: ',dotn
         IF( ( DOTN .GT. -MINDOT ) .AND. DOIT ) THEN
            THSDOT = -DOTN
            GOTO 900
         END IF
      END IF
C
      IF( IERR .NE. 0 ) RETURN
C
      THSDOT = MIN( THSDOT, -DOTN )
C
C      IF( DOT1 .LT. 0.2 .OR. DOT2 .LT. 0.2 ) GOTO 900
C
      DST12 = LENEDG( BIGLST, NODLST, EDG12 )
      DST13 = LENEDG( BIGLST, NODLST, EDG13 )
      DST23 = LENEDG( BIGLST, NODLST, EDG23 )
C
      DST1 = MAX( DST12, DST13 )
      DST1 = MAX( DST1,  DST23 )
C
      DST2 = MIN( DST12, DST13 )
      DST2 = MIN( DST1,  DST23 )
C
      RAD1 = MIN( MINREL, DST2/DST1 )
C
      THSREL = DST1
C
      IF( REVERS ) THEN
         CALL CHKFAC( BIGLST, NODLST, N12, N13, N23,
     :                IFAC1, IFAC2, IFAC3,
     :                FAC1, FAC2, FAC3, RAD1, NODE1, THSREL )
      ELSE
c         print*,'   before chkfac: ',dst12,dst13,dst23
         CALL CHKFAC( BIGLST, NODLST, N12, N23, N13,
     :                IFAC1, IFAC2, IFAC3,
     :                FAC1, FAC2, FAC3, RAD1, NODE1, THSREL )
c         print*,'   after chkfac: ',thsrel
      END IF
C
      THSREL = THSREL/DST1
C
      IF( ( THSREL .LT. RAD1 ) .AND. DOIT ) GOTO 900
C
      IF( .NOT. DOIT ) GOTO 900
C
C - remove three faces
C
      IF( IFAC1 .EQ. ENDFRT ) THEN
         NXT1 = IFAC1
      END IF
C
C      IF( NXT1 .EQ. IFAC2 ) NXT1 = NXT2
      CALL RMVFAC( IFAC2 )
c      print*,'   removed face ',ifac2
C      IF( NXT1 .EQ. IFAC3 ) NXT1 = NXT3
      CALL RMVFAC( IFAC3 )
c      print*,'   removed face ',ifac3
C
      CALL RMVFAC( IFAC1 )
c      print*,'   removed face ',ifac1
C
C - remove node N123 from the front
C
      CALL RMFRND( BIGLST, NODLST, N123 )
c      print*,'   removed node ',n123
C
      IF( IERR .NE. 0 ) RETURN
C
      IFAC1 = 0
C
C - add new element into BIGLST
C
      CONELM(1) = 0
      CONELM(2) = IELM1
      CONELM(3) = IELM2
      CONELM(4) = IELM3
C
      EDGLST(1) = EDG1
      EDGLST(2) = EDG2
      EDGLST(3) = EDG3
      EDGLST(6) = EDG13
      EDGLST(5) = EDG23
      EDGLST(4) = EDG12
C
      CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, IELM, INEW, IREG )
c
      if( ierr .eq. 12345 ) then
         print*,'--- got bad oppedg in dofacs 2'
         ierr = 0
      endif
C
      IF( IELM .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** DOFACS: FAILED TO CREATE NEW ELEMENT'
         RETURN
      END IF
C
C - add new face into FRTLST
C
      CONFAC(1) = FAC1
      CONFAC(2) = FAC2
      CONFAC(3) = FAC3
      CONFAC(4) = FAC2
      CONFAC(5) = FAC1
C
      EDGLST(1) = EDG1
      EDGLST(2) = EDG2
      EDGLST(3) = EDG3
      EDGLST(4) = EDG2
      EDGLST(5) = EDG1
C
      IF( REVERS ) THEN
         CALL ADDFAC( CONFAC(3), .TRUE., EDGLST(3), IELM, IFC1, INEW )
      ELSE
         CALL ADDFAC( CONFAC, .TRUE., EDGLST, IELM, IFC1, INEW )
      END IF
C
      IF( IFC1 .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** DOFACS: FAILED TO CREATE NEW FACE'
         RETURN
      END IF
c
      print*,'>>>>>>> NEW ELEMENT & FACE CREATED HERE <<<<<<<'
c      ierr = 1
C
 900  CONTINUE
C
c      print*,'Leaving DOFACS'
      if( doit ) print*,'     THSDOT,THSREL: ',THSDOT,THSREL
C
      RETURN
      END
C
