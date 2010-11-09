!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

SUBROUTINE FORWARDDATA_XYZ0(X,Y,Z,NONODS,NPHASE,TOTELE,NLOC,NDGLNO,  &
     D3,ITSTIME,DIRNAME)
  IMPLICIT NONE
  INTEGER ::NONODS,NPHASE,TOTELE,NLOC
  INTEGER ::NDGLNO(TOTELE*NLOC)
  REAL    ::X(NONODS*NPHASE),Y(NONODS*NPHASE),Z(NONODS*NPHASE)
  LOGICAL ::D3
  INTEGER ::ITSTIME
  CHARACTER(40) DIRNAME
  CHARACTER(40) FILENAME_XYZ,FILENAME
  INTEGER I,K1,K2

  WRITE(FILENAME_XYZ, 10) ITSTIME
10 FORMAT('/xyz_data_ori',I5.5)

  ewrite(3,*) DIRNAME
  ewrite(3,*) FILENAME_XYZ
  K1 = index( DIRNAME, ' ' ) - 1
  K2 = index( FILENAME_XYZ, ' ' ) - 1
  FILENAME=DIRNAME(1:K1) // FILENAME_XYZ(1:K2)
  ewrite(3,*) FILENAME
  OPEN(1,FILE = FILENAME)

  WRITE(1,*) NONODS,NPHASE,TOTELE,NLOC
  WRITE(1,*) (NDGLNO(I),I=1,TOTELE*NLOC)
  WRITE(1,*) (X(I),I=1,NONODS*NPHASE)
  WRITE(1,*) (Y(I),I=1,NONODS*NPHASE)
  IF(D3) THEN
     WRITE(1,*) (Z(I),I=1,NONODS*NPHASE)
  ENDIF

  CLOSE(1)

END SUBROUTINE FORWARDDATA_XYZ
