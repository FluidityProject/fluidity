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
      SUBROUTINE GETMSH( MSHFIL, MSHBIN, X, Y, Z,
     :                   GETSRF, GETBAS, GETTYP, GETREG, GETMAT, GETSOR,
     :                   GETDOM, GETNWN, GETCON,
     :                   ENLBAS, ENLIST, NEWNOD, NODDOM,
     :                   ELMBAS, ELMTYP, ELMREG, ELMAT, ELMSOR,
     :                   CONNOD, CONTYP,
     :                   SZENLS, NNOD, NELM, NGEOM, NCON, OK )
C  ---------------------------------------------------------------------
C
C
C
C  ---------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER MSHFIL, NNOD, NELM, SZENLS, NGEOM, NCON,
     :        ELMREG(NELM), ENLIST(SZENLS), NEWNOD(NNOD),
     :        NODDOM(NNOD), ELMAT(NELM), ENLBAS(NELM+1), ELMTYP(NELM),
     :        ELMSOR(NELM), CONNOD(NCON), CONTYP(NCON), ELMBAS(NELM)
C
      REAL X(NNOD), Y(NNOD), Z(NNOD)
C
      LOGICAL MSHBIN, GETSRF, GETBAS, GETTYP, GETREG, GETMAT, GETSOR,
     :        GETDOM, GETNWN, GETCON, OK
C
      INTEGER IELM, INOD, IBAS, KBAS, MAXSIZ, NLN, I, J, JELM, ICON,
     :        IDOM, INEW, ITYP, IREG, IMAT, ISOR, NODES(21), JBAS, ILN
C
c      print*,'In GETMSH ',mshbin
c
      OK = .TRUE.
C
      MAXSIZ = SZENLS
C
      SZENLS = 0
C
      ENLBAS(1) = SZENLS
C
      JELM = 0
C
      DO IELM = 1, NELM
C
         IF( MSHBIN ) THEN
            READ( MSHFIL ) KBAS, ITYP, IREG, IMAT, ISOR, NLN,
     :                     ( NODES(I), I = 1, NLN )
         ELSE
            NLN = 21
            READ( MSHFIL, '(5X,15I5)' ) KBAS, ITYP, IREG, IMAT, ISOR,
     :                     ( NODES(I), I = 1, 10 )
         END IF
C
         IF( NGEOM .LE. 3 .OR.  MOD(ITYP,2)  .EQ. 0 ) GO TO 10
         IF( NGEOM .LE. 5 .AND. MOD(ITYP,10) .LE. 5 ) GO TO 10
         IF( NGEOM .EQ. 6 .AND. MOD(ITYP,10) .LE. 2 ) GO TO 10
         IF( MSHBIN ) GOTO 10
         READ( MSHFIL, '(5X,15I5)' ) ( NODES(I), I = 11, 21 )
C
  10     CONTINUE
C
C - ignore surface elements, but note the nodes that are in them for surfaces
C
         if( mod(ityp,2).eq.0 .and. .not. getsrf ) then
            goto 30
         end if
C
C - count number of nodes and add to element node list
C
         ILN = 0
         DO I = 1, NLN
            IF( NODES(I) .EQ. 0 ) GOTO 20
            INOD   = NODES(I)
            ILN    = ILN    + 1
            SZENLS = SZENLS + 1
            IF( SZENLS .GT. MAXSIZ ) THEN
               WRITE(0,*)  '*** ELEMENT NODE LIST TOO LARGE'
               OK = .FALSE.
               SZENLS = -1
               RETURN
            END IF
            ENLIST(SZENLS) = INOD
         END DO
C
  20     JELM = JELM + 1
         NLN  = ILN
C
         IF( MOD(ITYP,2) .EQ. 0 ) THEN
C
            IF( NGEOM .EQ. 6 ) THEN
              IF( NLN .LT. 3 ) OK = .FALSE.
            ELSE IF( NGEOM .GT. 3 ) THEN
              IF( NLN .LT. 2 ) OK = .FALSE.
            ELSE IF( NLN .LT. 1 ) THEN
              OK = .FALSE.
            END IF
C
         ELSE
C
            IF( NGEOM .EQ. 6 ) THEN
               IF( NLN .LT. 4 ) OK = .FALSE.
            ELSE IF( NGEOM .GT. 3 ) THEN
               IF( NLN .LT. 3 ) OK = .FALSE.
            ELSE IF( NLN .LT. 2 ) THEN
               OK = .FALSE.
            END IF
C
         END IF
C
         IF( .NOT. OK ) THEN
            WRITE(0,*) 'Element has too few nodes: ',jelm,nln
            RETURN
         END IF
C
         IF( GETBAS ) ELMBAS(JELM) = KBAS
         IF( GETTYP ) ELMTYP(JELM) = ITYP
         IF( GETREG ) ELMREG(JELM) = IREG
         IF( GETSOR ) ELMSOR(JELM) = ISOR
         IF( GETMAT ) ELMAT(JELM)  = IMAT
C
         ENLBAS(JELM+1) = SZENLS
c          print*,'elm ',jelm,kbas,ityp,ireg,imat,isor,nln,
c     :                       (nodes(i),i=1,nln)
C
  30     CONTINUE
C
      END DO
c
      WRITE(0,*) 'Elements read: ',jelm,'  out of ',nelm
C
      NELM = JELM
C
      DO INOD = 1, NNOD
C
         IF( MSHBIN ) THEN
            READ( MSHFIL ) X(INOD), Y(INOD), Z(INOD), INEW,
     :                                      KBAS, IDOM, JBAS
         ELSE
            READ( MSHFIL, '(10X,3E15.7,4I5)' ) X(INOD), Y(INOD),
     :               Z(INOD), INEW, KBAS, IDOM, JBAS
         END IF
C
         IF( GETDOM ) NODDOM(INOD) = IDOM
         IF( GETNWN ) NEWNOD(INOD) = INEW
c         print*,'node ',inod,x(inod),y(inod),z(inod),inew,kbas,idom,jbas
C
      END DO
C
      IF( GETCON) THEN
         DO ICON = 1, NCON
            IF( MSHBIN ) THEN
               READ( MSHFIL ) CONNOD(ICON), CONTYP(ICON)
            ELSE
               READ( MSHFIL, '(10X,2I5)' ) CONNOD(ICON), CONTYP(ICON)
            END IF
         END DO
      END IF
C
      RETURN
      END
C
