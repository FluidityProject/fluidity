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
      SUBROUTINE GETSIZ( BIGLST, NODLST, NNOD, NELM, SZENLS, SZNNLS,
     :                   NSELM, SZSNLS )
C-----------------------------------------------------------------------
C
C - This returns the size of the new adapted mesh; the new number of nodes
C - (NNOD), the new number of elements (NELM) and the new size of the
C - element node list (SZENLS - which is always 4*NELM for 3D), and the
C - size of the node-to-node (node connectivity) list (SZNNLS).
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER NNOD, NELM, SZENLS, SZNNLS, NSELM, SZSNLS
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INTEGER IPTR
C
      NNOD = NUMNDS
      NELM = NELEMS
C
      NSELM = 0
C
      IPTR = STTBIG
C
 10   IF( BIGLST(4,IPTR) .LT. 0 ) THEN
         IPTR = BIGLST(NXTBIG,IPTR)
         GOTO 10
      ELSE IF( BIGLST(4,IPTR) .GT. 0 ) THEN
         IF( BIGLST(1,IPTR) .LE. 0 ) NSELM = NSELM + 1
         IF( BIGLST(2,IPTR) .LE. 0 ) NSELM = NSELM + 1
         IF( BIGLST(3,IPTR) .LE. 0 ) NSELM = NSELM + 1
         IF( BIGLST(4,IPTR)/2 .EQ. 0 ) NSELM = NSELM + 1
         IPTR = BIGLST(NXTBIG,IPTR)
         IPTR = BIGLST(NXTBIG,IPTR)
         IPTR = BIGLST(NXTBIG,IPTR)
         GOTO 10
      END IF
C
      IF( IS3DMS ) THEN
         SZENLS = 4*NELEMS
         SZSNLS = 3*NSELM
      ELSE
         SZENLS = 3*NELEMS
         SZSNLS = 2*NSELM
      END IF
C
      SZNNLS = 2*NEDGES + NNOD
C
      RETURN
      END
C
