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
      SUBROUTINE MKESCL( NNOD, NELM, ENLBAS, ENLIST, SZENLS,
     :                   NELBAS, NELIST, SZNELS,
     :                   NSELM, SZSNLS, SNLBAS, SNLIST, SURFID,
     :                   IS3D, EELIST )
C-----------------------------------------------------------------------
C
C - This forms the element connectivity list (EELIST) given the element node
C - list (ENLIST), its base pointer (ENLBAS) and its size (SZENLS), and the
C - info returned by the previous call to MKNELS containing the nodal element
C - list etc., plus info about whether this is a 3D problem (IS3D - logical).
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER NELM, NNOD, NSELM, SZENLS, SZNELS, ENLBAS(NELM+1),
     :        ENLIST(SZENLS), NELBAS(NNOD+1), NELIST(SZNELS),
     :        EELIST(NELM*4), SZSNLS, SNLBAS(NSELM+1), SNLIST(SZSNLS),
     :        SURFID(NSELM)
C
      LOGICAL IS3D
C
      INTEGER I, J, IBAS, ILN, IELM, JELM, NLN, INOD, I1, I2, I3, I4,
     :        J1, J2, J3
C
C - make sure it's all zero
C
      DO I = 1, 4*NELM
         EELIST(I) = 0
      END DO
C
      J = 1
C
C - check for surface connectivity
C
      IF( SZSNLS .GT. 1 ) THEN
C
         DO JELM = 1, NSELM
C
            IBAS = SNLBAS(JELM)
            NLN  = SNLBAS(JELM+1) - IBAS
            IELM = 0
C
            IF( IS3D ) THEN
C
               J1 = SNLIST(IBAS+1)
               J2 = SNLIST(IBAS+2)
               J3 = SNLIST(IBAS+3)
C
               CALL FNDSRF( J1, J2, J3, J3, ENLBAS, ENLIST, NELBAS,
     :                      NELIST, NNOD, NELM, SZENLS, SZNELS, IELM )
C
               IF( IELM .GT. 0 .AND. IELM .LE. NELM ) THEN
C
                  IBAS = ENLBAS(IELM)
                  NLN  = ENLBAS(IELM+1) - IBAS
C
                  DO I = 1, NLN
                     IF( J1 .EQ. ENLIST(IBAS+I) ) I1 = I
                     IF( J2 .EQ. ENLIST(IBAS+I) ) I2 = I
                     IF( J3 .EQ. ENLIST(IBAS+I) ) I3 = I
                  END DO
C
                  I4 = I1 + I2 + I3
                  IF( I4 .GE. 6 .AND. I4 .LE. 9 ) THEN
                     EELIST(IELM*4+I4-9) = -SURFID(JELM)
                  ELSE
                     WRITE(0,*)
     :                     '*** MKESCL: failed to connect volel to srf'
                     WRITE(0,*) jelm,surfid(jelm)
                     WRITE(0,*) j1,j2,j3
                     WRITE(0,*) ielm,ibas,nln
                     WRITE(0,*) (enlist(ibas+i),i=1,nln)
                  END IF
C
               END IF
C
            ELSE
C
               J1 = SNLIST(IBAS+1)
               J2 = SNLIST(IBAS+2)
C
               CALL FNDEDG( J1, J2, ENLBAS, ENLIST, NELBAS, NELIST,
     :                      NNOD, NELM, SZENLS, SZNELS, IELM,
     :                      .FALSE., .FALSE. )
C
            END IF
C
         END DO
C
      END IF
C
      DO IELM = 1, NELM
C
         IBAS = ENLBAS(IELM)
         NLN  = ENLBAS(IELM+1) - IBAS
C
         IF( IS3D ) THEN
C
            I1 = ENLIST(IBAS+1)
            I2 = ENLIST(IBAS+2)
            I3 = ENLIST(IBAS+3)
            I4 = ENLIST(IBAS+4)
C
            JELM = IELM
            CALL FNDSRF( I1, I2, I3, I3, ENLBAS, ENLIST, NELBAS, NELIST,
     :                   NNOD, NELM, SZENLS, SZNELS, JELM )
C
            IF( JELM .GT. 0 .AND. JELM .LE. NELM ) THEN
               IF( EELIST(J) .LT. 0 .AND. EELIST(J) .GE. -NSELM ) THEN
                  WRITE(0,*) '*** MKESCL: Already connected to surface?'
                  WRITE(0,*) ielm,ibas,nln
                  WRITE(0,*) i1,i2,i3,i4
                  WRITE(0,*) jelm,eelist(j),j
               END IF
               EELIST( J ) = JELM
            END IF
            J = J + 1
C
            JELM = IELM
            CALL FNDSRF( I1, I2, I4, I4, ENLBAS, ENLIST, NELBAS, NELIST,
     :                   NNOD, NELM, SZENLS, SZNELS, JELM )
C
            IF( JELM .GT. 0 .AND. JELM .LE. NELM ) THEN
               IF( EELIST(J) .LT. 0 .AND. EELIST(J) .GE. -NSELM ) THEN
                  WRITE(0,*) '*** MKESCL: Already connected to surface?'
                  WRITE(0,*) ielm,ibas,nln
                  WRITE(0,*) i1,i2,i3,i4
                  WRITE(0,*) jelm,eelist(j),j
               END IF
               EELIST( J ) = JELM
            END IF
            J = J + 1
C
            JELM = IELM
            CALL FNDSRF( I1, I3, I4, I4, ENLBAS, ENLIST, NELBAS, NELIST,
     :                   NNOD, NELM, SZENLS, SZNELS, JELM )
C
            IF( JELM .GT. 0 .AND. JELM .LE. NELM ) THEN
               IF( EELIST(J) .LT. 0 .AND. EELIST(J) .GE. -NSELM ) THEN
                  WRITE(0,*) '*** MKESCL: Already connected to surface?'
                  WRITE(0,*) ielm,ibas,nln
                  WRITE(0,*) i1,i2,i3,i4
                  WRITE(0,*) jelm,eelist(j),j
               END IF
               EELIST( J ) = JELM
            END IF
            J = J + 1
C
            JELM = IELM
            CALL FNDSRF( I2, I3, I4, I4, ENLBAS, ENLIST, NELBAS, NELIST,
     :                   NNOD, NELM, SZENLS, SZNELS, JELM )
C
            IF( JELM .GT. 0 .AND. JELM .LE. NELM ) THEN
               IF( EELIST(J) .LT. 0 .AND. EELIST(J) .GE. -NSELM ) THEN
                  WRITE(0,*) '*** MKESCL: Already connected to surface?'
                  WRITE(0,*) ielm,ibas,nln
                  WRITE(0,*) i1,i2,i3,i4
                  WRITE(0,*) jelm,eelist(j),j
               END IF
               EELIST( J ) = JELM
            END IF
            J = J + 1
C
         ELSE
C
            I1 = ENLIST(IBAS+1)
            I2 = ENLIST(IBAS+2)
            I3 = ENLIST(IBAS+3)
C
            JELM = IELM
            CALL FNDEDG( I1, I2, ENLBAS, ENLIST, NELBAS, NELIST,
     :                   NNOD, NELM, SZENLS, SZNELS, JELM,
     :                   .FALSE., .FALSE. )
C
            IF( JELM .GT. 0 .AND. JELM .LE. NELM ) THEN
               IF( EELIST(J) .LT. 0 .AND. EELIST(J) .GE. -NSELM ) THEN
                  WRITE(0,*) '*** MKESCL: Already connected to surface?'
                  WRITE(0,*) ielm,ibas,nln
                  WRITE(0,*) i1,i2,i3,i4
                  WRITE(0,*) jelm,eelist(j),j
               END IF
               EELIST( J ) = JELM
            END IF
            J = J + 1
C
            JELM = IELM
            CALL FNDEDG( I2, I3, ENLBAS, ENLIST, NELBAS, NELIST,
     :                   NNOD, NELM, SZENLS, SZNELS, JELM,
     :                   .FALSE., .FALSE. )
C
            IF( JELM .GT. 0 .AND. JELM .LE. NELM ) THEN
               IF( EELIST(J) .LT. 0 .AND. EELIST(J) .GE. -NSELM ) THEN
                  WRITE(0,*) '*** MKESCL: Already connected to surface?'
                  WRITE(0,*) ielm,ibas,nln
                  WRITE(0,*) i1,i2,i3,i4
                  WRITE(0,*) jelm,eelist(j),j
               END IF
               EELIST( J ) = JELM
            END IF
            J = J + 1
C
            JELM = IELM
            CALL FNDEDG( I1, I3, ENLBAS, ENLIST, NELBAS, NELIST,
     :                   NNOD, NELM, SZENLS, SZNELS, JELM,
     :                   .FALSE., .FALSE. )
C
            IF( JELM .GT. 0 .AND. JELM .LE. NELM ) THEN
               IF( EELIST(J) .LT. 0 .AND. EELIST(J) .GE. -NSELM ) THEN
                  WRITE(0,*) '*** MKESCL: Already connected to surface?'
                  WRITE(0,*) ielm,ibas,nln
                  WRITE(0,*) i1,i2,i3,i4
                  WRITE(0,*) jelm,eelist(j),j
               END IF
               EELIST( J ) = JELM
            END IF
            J = J + 1
C
         END IF
C
      END DO
C
c      IF( IS3D ) THEN
c         J = 4
c      ELSE
c         J = 3
c      END IF
c
c      print*,'Element connectivity list: '
c      print*,eelist
c      do ielm = 1, nelm
c        print*,ielm,(eelist(ielm*J-iln),iln=J-1,0,-1)
c      end do
C
      RETURN
      END
C
C
C
      SUBROUTINE CHKCNSTNCY( NNOD, NELM, ENLBAS, ENLIST, SZENLS,
     :                       IS3D, EELIST, XNOD, YNOD, ZNOD )
C
      IMPLICIT NONE
C
      INTEGER NNOD, NELM, SZENLS
C
      INTEGER ENLBAS(NELM+1), ENLIST(SZENLS), EELIST(NELM*4)
C
      REAL XNOD(NNOD), YNOD(NNOD), ZNOD(NNOD)
C
      LOGICAL IS3D
C
C - local   
C
      INTEGER IELM, JELM, NLN, ILN, I, J, K, IBAS, INOD, NODLST(8), CNT,
     :        L
C
      REAL X(8), Y(8), Z(8), V1, V2, TETVOL
C
C - check if any element is connected to another element twice
C
      I = 0
      DO IELM = 1, NELM
        DO J = 1, 3
          IF( EELIST(I+J) .GT. 0 ) THEN
            DO K = J+1, 4
              IF( EELIST(I+J) .EQ. EELIST(I+K) ) THEN
                WRITE(0,*)
     :              '*** CHKCNSTNCY: FOUND BAD ELEMENT CONNECTIONS:'
                WRITE(0,*)'  element: ',IELM,' ( with ',EELIST(I+J),' )'
                WRITE(0,*) '  conelm: ',(EELIST(I+L),L=1,4)
                IBAS = ENLBAS(IELM)+1
                WRITE(0,*) '  elnds1: ',(ENLIST(L),L=IBAS,IBAS+3)
                IBAS = ENLBAS(EELIST(I+J))+1
                WRITE(0,*) '  elnds2: ',(ENLIST(L),L=IBAS,IBAS+3)
                STOP
              END IF
            END DO
          END IF
        END DO
        I = I + 4
      END DO
C
C - check all element pairs are consistent
C
      I = 4
C
      DO IELM = 2, NELM
C
         IBAS = ENLBAS(IELM)
         NLN  = ENLBAS(IELM+1) - IBAS
C
         IF( NLN .NE. 4 ) THEN
            WRITE(0,*) '*** CHKCNSTNCY: FOUND NON-TETRAHEDRAL ELEMENT!'
            WRITE(0,*) '  element ',IELM,' has ',NLN,' nodes'
            STOP
         END IF
C
         DO ILN = 1, 4
            NODLST(ILN) = ENLIST(IBAS+ILN)
         END DO
C
         DO J = 1, 4
C
            JELM = EELIST(I+J)
C
            IF( JELM .LT. IELM .AND. JELM .GT. 0 ) THEN
C
               IBAS = ENLBAS(JELM)
               NLN  = ENLBAS(JELM+1) - IBAS
               DO ILN = 1, 4
                  NODLST(ILN+4) = ENLIST(IBAS+ILN)
               END DO
C
               CNT  = 0
               DO ILN = 1, 4
                  INOD = NODLST(ILN)
                  IF( NODLST(5) .EQ. INOD ) THEN
                     NODLST(5) = NODLST(ILN+4)
                     NODLST(ILN+4) = INOD
                     CNT = CNT + 1
                  ELSE IF( NODLST(6) .EQ. INOD ) THEN
                     NODLST(6) = NODLST(ILN+4)
                     NODLST(ILN+4) = INOD
                     CNT = CNT + 1
                  ELSE IF( NODLST(7) .EQ. INOD ) THEN
                     NODLST(7) = NODLST(ILN+4)
                     NODLST(ILN+4) = INOD
                     CNT = CNT + 1
                  ELSE IF( NODLST(8) .EQ. INOD ) THEN
                     NODLST(8) = NODLST(ILN+4)
                     NODLST(ILN+4) = INOD
                     CNT = CNT + 1
                  END IF
               END DO
C
               IF( CNT .NE. 3 ) THEN
                 WRITE(0,*) '*** CHKCNSTNCY: ELEMS DO NOT SHARE 3 NODES'
                 WRITE(0,*) '  elems: ',IELM,JELM
                 WRITE(0,*) '  nodls1: ',(NODLST(K),K=1,4)
                 WRITE(0,*) '  nodls2: ',(NODLST(K),K=5,8)
                 STOP
               END IF
C
               DO ILN = 1, 8
                  INOD   = NODLST(ILN)
                  X(ILN) = XNOD(INOD)
                  Y(ILN) = YNOD(INOD)
                  Z(ILN) = ZNOD(INOD)
               END DO
C
               V1 = TETVOL( X(1), Y(1), Z(1) )
               V2 = TETVOL( X(5), Y(5), Z(5) )
C
               IF( V1*V2 .GE. 0.0 ) THEN
                  WRITE(0,*) '*** CHKCNSTNCY: ELEMENTS NOT CONSISTENT:'
                  WRITE(0,*) '  elements: ',IELM,JELM
                  WRITE(0,*) '  vols (should have opp. signs): ',V1,V2
                  WRITE(0,*) '  nodls1: ',(NODLST(K),K=1,4)
                  DO K = 1, 4
                     WRITE(0,*) X(K),Y(K),Z(K)
                  END DO
                  WRITE(0,*) '  nodls2: ',(NODLST(K),K=5,8)
                  DO K = 5, 8
                     WRITE(0,*) X(K),Y(K),Z(K)
                  END DO
                  STOP
               END IF
C
            END IF
C
         END DO
C
         I = I + 4
C
      END DO
C
      RETURN
      END
C
