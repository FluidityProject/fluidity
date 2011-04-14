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
      SUBROUTINE ANGFAC( BIGLST, NODLST, N1, N2, N3, DXE, DYE, DZE,
     :                   DOTN )
C-----------------------------------------------------------------------
C
C - This subroutine takes the nodes supplied and returns in DOTN
C - the cosine of the angle between the *unit* vector given by DXE,
C - DYE, DZE and the face normal.
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
      INTEGER N1, N2, N3
C
      REAL DOTN, DXE, DYE, DZE
C
      REAL ONX1, ONY1, ONZ1, DOT, DX1, DY1, DZ1, DX2, DY2, DZ2
c
c      print*,' '
c      print*,'Into ANGFAC... ',ierr
c      print*,'    N1,N1,N3:',N1,N2,N3
c      print*,'    DX,DY,DZ:',DXE,DYE,DZE
C
      IF( IERR .NE. 0 ) RETURN
C
      DX1 = NODLST(1,N3) - NODLST(1,N1)
      DY1 = NODLST(2,N3) - NODLST(2,N1)
      DZ1 = NODLST(3,N3) - NODLST(3,N1)
C
      DX2 = NODLST(1,N2) - NODLST(1,N1)
      DY2 = NODLST(2,N2) - NODLST(2,N1)
      DZ2 = NODLST(3,N2) - NODLST(3,N1)
C
      ONX1 = DY1*DZ2 - DY2*DZ1
      ONY1 = DZ1*DX2 - DZ2*DX1
      ONZ1 = DX1*DY2 - DX2*DY1
C
      DOT = ONX1*ONX1 + ONY1*ONY1 + ONZ1*ONZ1
      DOT = 1.0/SQRT(DOT)
C
      ONX1 = ONX1*DOT
      ONY1 = ONY1*DOT
      ONZ1 = ONZ1*DOT
C
      DOTN = ONX1*DXE + ONY1*DYE + ONZ1*DZE
C
c      print*,'Leaving ANGFAC... '
c      print*,'       DOTN:',DOTN
C
      RETURN
      END
C
