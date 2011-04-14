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
      SUBROUTINE ANGEDG( BIGLST, NODLST, N1, N2, DXE, DYE, DZE, DOTE )
C-----------------------------------------------------------------------
C
C - This subroutine takes the nodes supplied and returns in DOTE
C - the cosine of the angle between the *unit* vector given by DXE,
C - DYE, DZE and the edge (from N1 to N2).
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
      INTEGER N1, N2
C
      REAL DOTE, DXE, DYE, DZE
C
      REAL DOT, DX1, DY1, DZ1
c
c      print*,' '
c      print*,'Into ANGEDG... ',ierr
c      print*,'       N1,N1:',N1,N2
c      print*,'    DX,DY,DZ:',DXE,DYE,DZE
C
      IF( IERR .NE. 0 ) RETURN
C
      DX1 = NODLST(1,N2) - NODLST(1,N1)
      DY1 = NODLST(2,N2) - NODLST(2,N1)
      DZ1 = NODLST(3,N2) - NODLST(3,N1)
C
      DOTE = DXE*DX1 + DYE*DY1 + DZE*DZ1
C
      DOT  = DX1*DX1 + DY1*DY1 + DZ1*DZ1
C
      DOTE = DOTE/SQRT(DOT)
C
c      print*,'Leaving ANGEDG... '
c      print*,'       DOTE:',DOTE
C
      RETURN
      END
C
