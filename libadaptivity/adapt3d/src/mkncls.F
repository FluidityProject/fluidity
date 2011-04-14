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
      SUBROUTINE MKNCLS( NNOD, NELM, ENLBAS, ENLIST, SZENLS,
     :                   NELBAS, NELIST, SZNELS,
     :                   NNLBAS, NNLIST, SZNNLS )
C-----------------------------------------------------------------------
C
C - This forms the node connectivity list (NNLIST) given the element node
C - list (ENLIST), its base pointer (ENLBAS) and its size (SZENLS), and the
C - info returned by the call to MKNELS containing the nodal element
C - list etc.
C - The connected nodes are listed in ascending order except for itself
C - which always comes first in the list.
C 
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER NELM, NNOD, SZENLS, SZNELS, SZNNLS, ENLBAS(NELM+1),
     :        ENLIST(SZENLS), NELBAS(NNOD+1), NELIST(SZNELS),
     :        NNLIST(SZNNLS), NNLBAS(NNOD+1)
C
      INTEGER INOD, I, IELM, J, JNOD, CURPOS, K, L
C
      SZNNLS = 0
C
      NNLBAS(1) = SZNNLS
C
      DO INOD = 1, NNOD
C
         SZNNLS = SZNNLS + 1
         CURPOS = SZNNLS
         NNLIST(SZNNLS) = INOD
C
         DO I = NELBAS(INOD)+1, NELBAS(INOD+1)
C
            IELM = NELIST(I)
C
            DO J = ENLBAS(IELM)+1, ENLBAS(IELM+1)
C
               JNOD = ENLIST(J)
C
               IF( JNOD .NE. INOD ) THEN
C
                  K = CURPOS
C
  100             IF( K .LT. SZNNLS ) THEN
                     K = K + 1
                     IF( NNLIST(K) .LT. JNOD ) GOTO 100
                     IF( NNLIST(K) .EQ. JNOD ) GOTO 300
                     SZNNLS = SZNNLS + 1
                     DO L = SZNNLS, K+1, -1
                        NNLIST(L) = NNLIST(L-1)
                     END DO
                     NNLIST(K) = JNOD
                  ELSE
                     SZNNLS = SZNNLS + 1
                     NNLIST(SZNNLS) = JNOD
                  END IF
C
               END IF
C
  300        END DO
C
         END DO
C
         NNLBAS(INOD+1) = SZNNLS
C
      END DO
C
      RETURN
      END
C
