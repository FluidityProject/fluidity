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
      SUBROUTINE MKNELS( NNOD, NELM, ENLBAS, ENLIST, SZENLS,
     :                   NELBAS, NELIST, SZNELS, OFFSET )
C-----------------------------------------------------------------------
C
C - This forms NELIST, the nodal element list (i.e. elems round a node),
C - and its base pointer (NELBAS), given the element node list (ENLIST),
C - its base pointer (ENLBAS), its size (SZENLS) and the number of nodes
C - (NNOD) and elements (NELM). The size of NELIST (SZNELS) is also returned
C - but always turns out to be the same size as the element node list (SZENLS).
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER NELM, NNOD, SZENLS, ENLBAS(NELM+1), ENLIST(SZENLS),
     :        NELBAS(NNOD+1), NELIST(SZENLS), SZNELS, OFFSET
C
      INTEGER I, IBAS, ILN, IELM, NLN, INOD
C
c      print*,'In MKNELS: ',NNOD,NELM,SZENLS,OFFSET
      DO INOD = 1, NNOD+1
         NELBAS(INOD) = 0
      END DO
C
C - count occurences of each node in ENLIST
C
c      DO I = 1, SZENLS
c         INOD = ENLIST(I) + 1
c         NELBAS(INOD) = NELBAS(INOD) + 1
c      END DO
C
      DO IELM = 1, NELM
         IBAS = ENLBAS(IELM)
         NLN  = ENLBAS(IELM+1) - IBAS
         IBAS = IBAS - 1 + OFFSET
         if(ibas+nln.gt.szenls)
     :      WRITE(0,*) 'MKNELS: bad IBAS: ',ielm,ibas,nln,szenls
         DO ILN = 1, NLN
            INOD = ENLIST(IBAS+ILN) + 1
            if(inod.gt.nnod+1.or.inod.le.1)
     :         WRITE(0,*) 'MKNELS: bad INOD: ',inod-1,ielm,ibas,iln
            NELBAS(INOD) = NELBAS(INOD) + 1
         END DO
      END DO
c      print*,'NELBAS as counter: ',(nelbas(i),i=1,21)
C
C - work out node element list base pointer (NELBAS) from count
C
      IF( NELBAS(1) .NE. 0 ) THEN
         WRITE(0,*) '*** MKNELS: MESH INCONSISTENCY'
c         PRINT*,'ENLIST: '
c         PRINT*,ENLIST
         STOP
      END IF
C
      DO I = 1, NNOD
         IF( NELBAS(I+1) .EQ. 0 ) THEN
            WRITE(0,*) '+++ MKNELS: Warning, node not used: ',I
         END IF
         NELBAS(I+1) = NELBAS(I+1) + NELBAS(I)
      END DO
c      print*,'NELBAS before: ',(nelbas(i),i=1,21)
C
C - work out node element list (NELIST) usint ENLIST and NELBAS
C
      DO IELM = 1, NELM
         IBAS = ENLBAS(IELM)
         NLN = ENLBAS(IELM+1) - IBAS
         IBAS = IBAS - 1 + OFFSET
         DO ILN = 1, NLN
            INOD = ENLIST(IBAS+ILN)
            NELBAS(INOD) = NELBAS(INOD) + 1
            NELIST(NELBAS(INOD)) = IELM
         END DO
      END DO
C
c      print*,'NELBAS after: ',(nelbas(i),i=1,21)
c
      DO I = NNOD, 1, -1
         NELBAS(I+1) = NELBAS(I) + 1 - OFFSET
      END DO
      NELBAS(1) = 1 - OFFSET
C
c      print*,'NELBAS added up: ',(nelbas(i),i=1,21)
c
      SZNELS = NELBAS(NNOD+1) - 1 + OFFSET
C
      IF( SZNELS .NE. SZENLS-ENLBAS(1)+1-OFFSET ) THEN
         WRITE(0,*) 'MKNELS: Something went wrong!!',SZENLS,SZNELS
         WRITE(0,*) 'ENLBAS: ',(enlbas(i),i=1,10)
         WRITE(0,*) 'ENLIST: ',(ENLIST(i),i=enlbas(1)+offset,enlbas(11))
         WRITE(0,*) 'NELBAS: ',(nelbas(i),i=1,10)
         WRITE(0,*) 'NELIST: ',(nelist(i),i=nelbas(1)+offset,nelbas(11))
cccc         pause
      end if
C
c      print*,'NELBAS: ',(nelbas(i),i=1,nnod+1)
c      print*,'NELIST: ',(nelist(i),i=1,sznels)
C
      RETURN
      END
C
