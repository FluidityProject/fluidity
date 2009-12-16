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
!    C.Pain@Imperial.ac.uk
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

#include "fdebug.h"
module legacy_cv_numbering
  use FLDebug
  implicit none

  private
  public :: SETELM, VOLNEI, SURNEI

contains

      SUBROUTINE SETELM( ELETYP, NGI, NLOC,&
!     - INTEGERS
     &     OPTELM, SELETYP, &
     &     SNGI,   SNLOC,   SVNGI,&
!     - LOGICALS
     &     D3,     DCYL,    REDQUAD   )
!     -----------------------------------------------
!     
!     - this subroutine the necessary parameters for the different
!     - element types: 
!     
!     - i.e.
!     
!     element type               | optelm
!     -----------------------------------
!     - linear        triangles        1
!     - quadratic     triangles        2
!     - bi-linear     quadrilaterals   3
!     - bi-quadratic  quadrilaterals   4
!     - linear        tetrahedra       5
!     - quadratic     tetrahedra       6 
!     - tri-linear    hexahedra        7
!     - tri-quadratic hexahedra        8
!     
!     -------------------------------
!     - date last modified : 04/11/2003
!     -------------------------------
!     
      IMPLICIT NONE
!     
      INTEGER ELETYP,  NGI
      INTEGER NLOC,    OPTELM
      INTEGER SELETYP
      INTEGER SNGI,    SNLOC
      INTEGER SVNGI
!     
      LOGICAL D3, DCYL, REDQUAD
!     
      IF( REDQUAD ) THEN
!     
         IF(OPTELM.EQ.1) THEN
!     
            D3      = .FALSE.
            DCYL    = .FALSE.
            ELETYP  =  2
            NGI     =  4
            NLOC    =  3
            SELETYP =  1
!             SFNGI   = 20
            SNGI    =  2
            SNLOC   =  2
            SVNGI   =  3
!     
         ELSE IF(OPTELM.EQ.2) THEN
!     
            D3      = .FALSE.
            DCYL    = .FALSE.
            ELETYP  =  7
            NGI     =  7
            NLOC    =  6
            SELETYP =  4
!             SFNGI   = 20
            SNGI    =  3
            SNLOC   =  3
            SVNGI   = 18
!     
         ELSE IF(OPTELM.EQ.3) THEN
!     
            D3      = .FALSE.
            DCYL    = .FALSE.
            ELETYP  =  1
            NGI     =  4
            NLOC    =  4
            SELETYP =  1
!             SFNGI   = 20
            SNGI    =  2
            SNLOC   =  2
            SVNGI   =  4
!     
         ELSE IF(OPTELM.EQ.4) THEN
!     
            D3      = .FALSE.
            DCYL    = .FALSE.
            ELETYP  =  5
            NGI     =  9
            NLOC    =  9
            SELETYP =  4
!             SFNGI   = 20
            SNGI    =  3
            SNLOC   =  3
            SVNGI   = 24
!     
         ELSE IF(OPTELM.EQ.5) THEN
!     linear tets
            D3      = .TRUE.
            DCYL    = .FALSE.
            ELETYP  =  4
            NGI     =  4
            NLOC    =  4
            SELETYP =  3
!             SFNGI   = 20
            SNGI    =  4
            SNLOC   =  3
            SVNGI   =  6
!     
         ELSE IF(OPTELM.EQ.6) THEN
!     
            D3      = .TRUE.
            DCYL    = .FALSE.
            ELETYP  =  8
            NGI     = 11
            NLOC    = 10
            SELETYP =  6
!             SFNGI   = 20
            SNGI    =  7 
!     SNGI = 36
            SNLOC   =  6
            SVNGI   = 96
!     
         ELSE IF(OPTELM.EQ.7) THEN
!     trilinear hexes
            D3      = .TRUE.
            DCYL    = .FALSE.
            ELETYP  =  3
            NGI     =  8
            NLOC    =  8
            SELETYP =  2
!             SFNGI   = 20
            SNGI    =  4 
            SNLOC   =  4
            SVNGI   = 12
!     
         ELSE IF(OPTELM.EQ.8) THEN
!     
            D3      = .TRUE.
            DCYL    = .FALSE.
            ELETYP  =   6
            NGI     =  27
            NLOC    =  27
            SELETYP =   5
!             SFNGI   =  20
            SNGI    =   9 
            SNLOC   =   9
            SVNGI   = 216
!     
         END IF
!     
      ELSE
!     
         IF(OPTELM.EQ.1) THEN
!     
            D3      = .FALSE.
            DCYL    = .FALSE.
            ELETYP  =  2
            NGI     = 12
            NLOC    =  3
            SELETYP =  1
!             SFNGI   = 20
            SNGI    =  2
            SNLOC   =  2
            SVNGI   =  3
!     
         ELSE IF(OPTELM.EQ.2) THEN
!     
            D3      = .FALSE.
            DCYL    = .FALSE.
            ELETYP  =  7
            NGI     = 36
            NLOC    =  6
            SELETYP =  4
!             SFNGI   = 20
            SNGI    =  6
            SNLOC   =  3
            SVNGI   = 18
!     
         ELSE IF(OPTELM.EQ.3) THEN
!     
            D3      = .FALSE.
            DCYL    = .FALSE.
            ELETYP  =  1
            NGI     = 16
            NLOC    =  4
            SELETYP =  1
!             SFNGI   = 20
            SNGI    =  2
            SNLOC   =  2
            SVNGI   =  4
!     
         ELSE IF(OPTELM.EQ.4) THEN
!     
            D3      = .FALSE.
            DCYL    = .FALSE.
            ELETYP  =  5
            NGI     = 36
            NLOC    =  9
            SELETYP =  4
!             SFNGI   = 20
            SNGI    =  6
            SNLOC   =  3
            SVNGI   = 24
!     
         ELSE IF(OPTELM.EQ.5) THEN
!     
            D3      = .TRUE.
            DCYL    = .FALSE.
            ELETYP  =  4
            NGI     = 32
            NLOC    =  4
            SELETYP =  3
!             SFNGI   = 20
            SNGI    =  3
            SNLOC   =  3
            SVNGI   =  6
!     
         ELSE IF(OPTELM.EQ.6) THEN
!     
            D3      =  .TRUE.
            DCYL    =  .FALSE.
            ELETYP  =   8
            NGI     = 128
            NLOC    =  10
            SELETYP =   6
!             SFNGI   =  20
            SNGI    =  36 
            SNLOC   =   6
            SVNGI   =  96
!     
         ELSE IF(OPTELM.EQ.7) THEN
!     
            D3      = .TRUE.
            DCYL    = .FALSE.
            ELETYP  =  3
            NGI     = 64
            NLOC    =  8
            SELETYP =  2
!             SFNGI   = 20
            SNGI    =  4 
            SNLOC   =  4
            SVNGI   = 12
!     
         ELSE IF(OPTELM.EQ.8) THEN
!     
            D3      = .TRUE.
            DCYL    = .FALSE.
            ELETYP  =   6
            NGI     = 216
            NLOC    =  27
            SELETYP =   5
!             SFNGI   =  20
            SNGI    =  36 
            SNLOC   =   9
            SVNGI   = 216
!     
         END IF
!     
      END IF
!     
      END SUBROUTINE 

      SUBROUTINE VOLNEI( NEILOC, NLOC, SVNGI, ELETYP )
!     ------------------------------------------------
!     
!     - this subroutine calculates NEILOC which is the 
!     - array containing information given a local node
!     - and an integration point what is the other opposing 
!     - local node. 
!     
!     -------------------------------
!     - date last modified : 11/10/2002
!     -------------------------------
!     
      IMPLICIT NONE
!     
      INTEGER ELETYP, NLOC, SVNGI
!     
      INTEGER NEILOC(NLOC,SVNGI)
!     
!     - local variables
!     
      INTEGER ILOC
!     
      NEILOC = 0
!     
      IF(ELETYP.EQ.1) THEN
      ! quadrilaterals
!     
         ILOC=1
         NEILOC(ILOC,1)=3
         NEILOC(ILOC,3)=2
!     
         ILOC=2
         NEILOC(ILOC,2)=4
         NEILOC(ILOC,3)=1
!     
         ILOC=3
         NEILOC(ILOC,1)=1
         NEILOC(ILOC,4)=4
!     
         ILOC=4
         NEILOC(ILOC,2)=2
         NEILOC(ILOC,4)=3
!     
      ELSE IF(ELETYP.EQ.2) THEN
!     
         ILOC=1
         NEILOC(ILOC,1)=2
         NEILOC(ILOC,3)=3
!     
         ILOC=2
         NEILOC(ILOC,1)=1
         NEILOC(ILOC,2)=3
!     
         ILOC=3
         NEILOC(ILOC,2)=2
         NEILOC(ILOC,3)=1
!     
      ELSE IF(ELETYP.EQ.3) THEN
!     ! hexes

         ILOC=1
         NEILOC(ILOC,1)=3
         NEILOC(ILOC,3)=2
         NEILOC(ILOC,5)=5
!     
         ILOC=2
         NEILOC(ILOC,2)=4
         NEILOC(ILOC,3)=1
         NEILOC(ILOC,6)=6
!     
         ILOC=3
         NEILOC(ILOC,1)=1
         NEILOC(ILOC,4)=4
         NEILOC(ILOC,7)=7
!     
         ILOC=4
         NEILOC(ILOC,2)=2
         NEILOC(ILOC,4)=3
         NEILOC(ILOC,8)=8
!     
         ILOC=5
         NEILOC(ILOC,5)=1
         NEILOC(ILOC,9)=7
         NEILOC(ILOC,11)=6
!     
         ILOC=6
         NEILOC(ILOC,6)=2
         NEILOC(ILOC,10)=8
         NEILOC(ILOC,11)=5
!     
         ILOC=7
         NEILOC(ILOC,7)=3
         NEILOC(ILOC,9)=5
         NEILOC(ILOC,12)=8
!     
         ILOC=8
         NEILOC(ILOC,8)=4
         NEILOC(ILOC,10)=6
         NEILOC(ILOC,12)=7
!     
      ELSE IF(ELETYP.EQ.4) THEN
!     
         ILOC=1
         NEILOC(ILOC,1)=2
         NEILOC(ILOC,3)=3
         NEILOC(ILOC,4)=4
!     
         ILOC=2
         NEILOC(ILOC,1)=1
         NEILOC(ILOC,2)=3
         NEILOC(ILOC,5)=4
!     
         ILOC=3
         NEILOC(ILOC,2)=2
         NEILOC(ILOC,3)=1
         NEILOC(ILOC,6)=4
!     
         ILOC=4
         NEILOC(ILOC,4)=1
         NEILOC(ILOC,5)=2
         NEILOC(ILOC,6)=3
!     
      ELSE IF(ELETYP.EQ.5) THEN
!     
         ILOC=1
         NEILOC(ILOC,1)=2
         NEILOC(ILOC,2)=2
         NEILOC(ILOC,5)=4
         NEILOC(ILOC,6)=4
!     
         ILOC=2
         NEILOC(ILOC,1)=1
         NEILOC(ILOC,2)=1
         NEILOC(ILOC,7)=5
         NEILOC(ILOC,8)=5
         NEILOC(ILOC,3)=3
         NEILOC(ILOC,4)=3      
!     
         ILOC=3
         NEILOC(ILOC,3) =2
         NEILOC(ILOC,4) =2
         NEILOC(ILOC,9) =6
         NEILOC(ILOC,10)=6
!     
         ILOC=4
         NEILOC(ILOC,5) =1
         NEILOC(ILOC,6) =1
         NEILOC(ILOC,11)=5
         NEILOC(ILOC,12)=5
         NEILOC(ILOC,15)=7
         NEILOC(ILOC,16)=7
!     
         ILOC=5
         NEILOC(ILOC,7) =2
         NEILOC(ILOC,8) =2
         NEILOC(ILOC,13)=6
         NEILOC(ILOC,14)=6
         NEILOC(ILOC,18)=8
         NEILOC(ILOC,17)=8
         NEILOC(ILOC,12)=4
         NEILOC(ILOC,11)=4
!     
         ILOC=6
         NEILOC(ILOC,9) =3
         NEILOC(ILOC,10)=3
         NEILOC(ILOC,13)=5
         NEILOC(ILOC,14)=5
         NEILOC(ILOC,19)=9
         NEILOC(ILOC,20)=9
!     
         ILOC=7
         NEILOC(ILOC,15)=4
         NEILOC(ILOC,16)=4
         NEILOC(ILOC,21)=8
         NEILOC(ILOC,22)=8
!     
         ILOC=8
         NEILOC(ILOC,21)=7
         NEILOC(ILOC,22)=7
         NEILOC(ILOC,17)=5
         NEILOC(ILOC,18)=5
         NEILOC(ILOC,23)=9
         NEILOC(ILOC,24)=9
!     
         ILOC=9
         NEILOC(ILOC,23)=8
         NEILOC(ILOC,24)=8
         NEILOC(ILOC,19)=6
         NEILOC(ILOC,20)=6
!     
      ELSE IF(ELETYP.EQ.6) THEN
!     
         ILOC=1
!     
!     - face 1 
!     
         NEILOC(ILOC,1)=2
         NEILOC(ILOC,2)=2
         NEILOC(ILOC,3)=2
         NEILOC(ILOC,4)=2
!     
!     - face 3
!     
         NEILOC(ILOC,9)=4
         NEILOC(ILOC,10)=4
         NEILOC(ILOC,11)=4
         NEILOC(ILOC,12)=4
!     
!     - face 13
!     
         NEILOC(ILOC,49)=10
         NEILOC(ILOC,50)=10
         NEILOC(ILOC,51)=10
         NEILOC(ILOC,52)=10
!     
         ILOC=2
!     
!     - face 1
!     
         NEILOC(ILOC,1)=1
         NEILOC(ILOC,2)=1
         NEILOC(ILOC,3)=1
         NEILOC(ILOC,4)=1
!     
!     - face 4
!     
         NEILOC(ILOC,13)=5
         NEILOC(ILOC,14)=5
         NEILOC(ILOC,15)=5
         NEILOC(ILOC,16)=5
!     
!     - face 2
!     
         NEILOC(ILOC,5)=3
         NEILOC(ILOC,6)=3
         NEILOC(ILOC,7)=3
         NEILOC(ILOC,8)=3
!     
!     - face 14
!     
         NEILOC(ILOC,53)=11
         NEILOC(ILOC,54)=11
         NEILOC(ILOC,55)=11
         NEILOC(ILOC,56)=11
!     
         ILOC=3
!     
!     - face 2
!     
         NEILOC(ILOC,5)=2
         NEILOC(ILOC,6)=2
         NEILOC(ILOC,7)=2
         NEILOC(ILOC,8)=2
!     
!     - face 5
!     
         NEILOC(ILOC,17)=6
         NEILOC(ILOC,18)=6
         NEILOC(ILOC,19)=6
         NEILOC(ILOC,20)=6
!     
!     - face 15
!     
         NEILOC(ILOC,57)=12
         NEILOC(ILOC,58)=12
         NEILOC(ILOC,59)=12
         NEILOC(ILOC,60)=12
!     
         ILOC=4
!     
!     - face 3
!     
         NEILOC(ILOC,9)=1
         NEILOC(ILOC,10)=1
         NEILOC(ILOC,11)=1
         NEILOC(ILOC,12)=1
!     
!     - face 6
!     
         NEILOC(ILOC,21)=5
         NEILOC(ILOC,22)=5
         NEILOC(ILOC,23)=5
         NEILOC(ILOC,24)=5
!     
!     - face 8
!     
         NEILOC(ILOC,29)=7
         NEILOC(ILOC,30)=7
         NEILOC(ILOC,31)=7
         NEILOC(ILOC,32)=7
!     
!     - face 16
!     
         NEILOC(ILOC,61)=13
         NEILOC(ILOC,62)=13
         NEILOC(ILOC,63)=13
         NEILOC(ILOC,64)=13
!     
         ILOC=5
!     
!     - face 4
!     
         NEILOC(ILOC,13)=2
         NEILOC(ILOC,14)=2
         NEILOC(ILOC,15)=2
         NEILOC(ILOC,16)=2
!     
!     - face 7
!     
         NEILOC(ILOC,25)=6
         NEILOC(ILOC,26)=6
         NEILOC(ILOC,27)=6
         NEILOC(ILOC,28)=6
!     
!     - face 9
!     
         NEILOC(ILOC,33)=8
         NEILOC(ILOC,34)=8
         NEILOC(ILOC,35)=8
         NEILOC(ILOC,36)=8
!     
!     - face 6
!     
         NEILOC(ILOC,21)=4
         NEILOC(ILOC,22)=4
         NEILOC(ILOC,23)=4
         NEILOC(ILOC,24)=4
!     
!     - face 17
!     
         NEILOC(ILOC,65)=14
         NEILOC(ILOC,66)=14
         NEILOC(ILOC,67)=14
         NEILOC(ILOC,68)=14
!     
         ILOC=6
!     
!     - face 5
!     
         NEILOC(ILOC,17)=3
         NEILOC(ILOC,18)=3
         NEILOC(ILOC,19)=3
         NEILOC(ILOC,20)=3
!     
!     - face 7
!     
         NEILOC(ILOC,25)=5
         NEILOC(ILOC,26)=5
         NEILOC(ILOC,27)=5
         NEILOC(ILOC,28)=5
!     
!     - face 10
!     
         NEILOC(ILOC,37)=9
         NEILOC(ILOC,38)=9
         NEILOC(ILOC,39)=9
         NEILOC(ILOC,40)=9
!     
!     - face 18
!     
         NEILOC(ILOC,69)=15
         NEILOC(ILOC,70)=15
         NEILOC(ILOC,71)=15
         NEILOC(ILOC,72)=15
!     
         ILOC=7
!     
!     - face 8
!     
         NEILOC(ILOC,29)=4
         NEILOC(ILOC,30)=4
         NEILOC(ILOC,31)=4
         NEILOC(ILOC,32)=4
!     
!     - face 11
!     
         NEILOC(ILOC,41)=8
         NEILOC(ILOC,42)=8
         NEILOC(ILOC,43)=8
         NEILOC(ILOC,44)=8
!     
!     - face 19
!     
         NEILOC(ILOC,73)=16
         NEILOC(ILOC,74)=16
         NEILOC(ILOC,75)=16
         NEILOC(ILOC,76)=16
!     
         ILOC=8
!     
!     - face 11
!     
         NEILOC(ILOC,41)=7
         NEILOC(ILOC,42)=7
         NEILOC(ILOC,43)=7
         NEILOC(ILOC,44)=7
!     
!     - face 9
!     
         NEILOC(ILOC,33)=5
         NEILOC(ILOC,34)=5
         NEILOC(ILOC,35)=5
         NEILOC(ILOC,36)=5
!     
!     - face 12
!     
         NEILOC(ILOC,45)=9
         NEILOC(ILOC,46)=9
         NEILOC(ILOC,47)=9
         NEILOC(ILOC,48)=9
!     
!     - face 20
!     
         NEILOC(ILOC,77)=17
         NEILOC(ILOC,78)=17
         NEILOC(ILOC,79)=17
         NEILOC(ILOC,80)=17
!     
         ILOC=9
!     
!     - face 10
!     
         NEILOC(ILOC,37)=6
         NEILOC(ILOC,38)=6
         NEILOC(ILOC,39)=6
         NEILOC(ILOC,40)=6
!     
!     - face 12
!     
         NEILOC(ILOC,45)=8
         NEILOC(ILOC,46)=8
         NEILOC(ILOC,47)=8
         NEILOC(ILOC,48)=8
!     
!     - face 21
!     
         NEILOC(ILOC,81)=18
         NEILOC(ILOC,82)=18
         NEILOC(ILOC,83)=18
         NEILOC(ILOC,84)=18
!     
         ILOC=10
!     
!     - face 13
!     
         NEILOC(ILOC,49)=1
         NEILOC(ILOC,50)=1
         NEILOC(ILOC,51)=1
         NEILOC(ILOC,52)=1
!     
!     - face 22
!     
         NEILOC(ILOC,85)=11
         NEILOC(ILOC,86)=11
         NEILOC(ILOC,87)=11
         NEILOC(ILOC,88)=11
!     
!     - face 24
!     
         NEILOC(ILOC,93)=13
         NEILOC(ILOC,94)=13
         NEILOC(ILOC,95)=13
         NEILOC(ILOC,96)=13
!     
!     - face 34
!     
         NEILOC(ILOC,133)=19
         NEILOC(ILOC,134)=19
         NEILOC(ILOC,135)=19
         NEILOC(ILOC,136)=19
!     
         ILOC=11
!     
!     - face 14
!     
         NEILOC(ILOC,53)=2
         NEILOC(ILOC,54)=2
         NEILOC(ILOC,55)=2
         NEILOC(ILOC,56)=2
!     
!     - face 22
!     
         NEILOC(ILOC,85)=10
         NEILOC(ILOC,86)=10
         NEILOC(ILOC,87)=10
         NEILOC(ILOC,88)=10
!     
!     - face 25
!     
         NEILOC(ILOC,97)=14
         NEILOC(ILOC,98)=14
         NEILOC(ILOC,99)=14
         NEILOC(ILOC,100)=14
!     
!     - face 23
!     
         NEILOC(ILOC,89)=12
         NEILOC(ILOC,90)=12
         NEILOC(ILOC,91)=12
         NEILOC(ILOC,92)=12
!     
!     - face 35
!     
         NEILOC(ILOC,137)=20
         NEILOC(ILOC,138)=20
         NEILOC(ILOC,139)=20
         NEILOC(ILOC,140)=20
!     
         ILOC=12
!     
!     - face 15
!     
         NEILOC(ILOC,57)=3
         NEILOC(ILOC,58)=3  
         NEILOC(ILOC,59)=3 
         NEILOC(ILOC,60)=3
!     
!     - face 23
!     
         NEILOC(ILOC,89)=11
         NEILOC(ILOC,90)=11
         NEILOC(ILOC,91)=11
         NEILOC(ILOC,92)=11
!     
!     - face 26
!     
         NEILOC(ILOC,101)=15 
         NEILOC(ILOC,102)=15
         NEILOC(ILOC,103)=15
         NEILOC(ILOC,104)=15
!     
!     - face 36
!     
         NEILOC(ILOC,141)=21 
         NEILOC(ILOC,142)=21
         NEILOC(ILOC,143)=21
         NEILOC(ILOC,144)=21
!     
         ILOC=13
!     
!     - face 16
!     
         NEILOC(ILOC,61)=4
         NEILOC(ILOC,62)=4
         NEILOC(ILOC,63)=4
         NEILOC(ILOC,64)=4
!     
!     - face 24
!     
         NEILOC(ILOC,93)=10
         NEILOC(ILOC,94)=10
         NEILOC(ILOC,95)=10
         NEILOC(ILOC,96)=10
!     
!     - face 27
!     
         NEILOC(ILOC,105)=14
         NEILOC(ILOC,106)=14
         NEILOC(ILOC,107)=14
         NEILOC(ILOC,108)=14
!     
!     - face 29
!     
         NEILOC(ILOC,113)=16
         NEILOC(ILOC,114)=16
         NEILOC(ILOC,115)=16
         NEILOC(ILOC,116)=16
!     
!     - face 37
!     
         NEILOC(ILOC,145)=22
         NEILOC(ILOC,146)=22
         NEILOC(ILOC,147)=22
         NEILOC(ILOC,148)=22
!     
         ILOC=14
!     
!     - face 17
!     
         NEILOC(ILOC,65)=5
         NEILOC(ILOC,66)=5
         NEILOC(ILOC,67)=5
         NEILOC(ILOC,68)=5
!     
!     - face 25
!     
         NEILOC(ILOC,97)=11
         NEILOC(ILOC,98)=11
         NEILOC(ILOC,99)=11
         NEILOC(ILOC,100)=11
!     
!     - face 28
!     
         NEILOC(ILOC,109)=15
         NEILOC(ILOC,110)=15
         NEILOC(ILOC,111)=15
         NEILOC(ILOC,112)=15
!     
!     - face 30
!     
         NEILOC(ILOC,117)=17
         NEILOC(ILOC,118)=17
         NEILOC(ILOC,119)=17
         NEILOC(ILOC,120)=17
!     
!     - face 27
!     
         NEILOC(ILOC,105)=13
         NEILOC(ILOC,106)=13
         NEILOC(ILOC,107)=13
         NEILOC(ILOC,108)=13
!     
!     - face 38
!     
         NEILOC(ILOC,149)=23
         NEILOC(ILOC,150)=23
         NEILOC(ILOC,151)=23
         NEILOC(ILOC,152)=23
!     
         ILOC=15
!     
!     - face 18
!     
         NEILOC(ILOC,69)=6
         NEILOC(ILOC,70)=6
         NEILOC(ILOC,71)=6
         NEILOC(ILOC,72)=6
!     
!     - face 26
!     
         NEILOC(ILOC,101)=12
         NEILOC(ILOC,102)=12
         NEILOC(ILOC,103)=12
         NEILOC(ILOC,104)=12
!     
!     - face 28
!     
         NEILOC(ILOC,109)=14
         NEILOC(ILOC,110)=14
         NEILOC(ILOC,111)=14
         NEILOC(ILOC,112)=14
!     
!     - face 31
!     
         NEILOC(ILOC,121)=18
         NEILOC(ILOC,122)=18
         NEILOC(ILOC,123)=18
         NEILOC(ILOC,124)=18
!     
!     - face 39
!     
         NEILOC(ILOC,153)=24
         NEILOC(ILOC,154)=24
         NEILOC(ILOC,155)=24
         NEILOC(ILOC,156)=24
!     
         ILOC=16
!     
!     - face 19
!     
         NEILOC(ILOC,73)=7
         NEILOC(ILOC,74)=7
         NEILOC(ILOC,75)=7
         NEILOC(ILOC,76)=7
!     
!     - face 29
!     
         NEILOC(ILOC,113)=13
         NEILOC(ILOC,114)=13
         NEILOC(ILOC,115)=13
         NEILOC(ILOC,116)=13
!     
!     - face 32
!     
         NEILOC(ILOC,125)=17
         NEILOC(ILOC,126)=17
         NEILOC(ILOC,127)=17
         NEILOC(ILOC,128)=17
!     
!     - face 40
!     
         NEILOC(ILOC,157)=25
         NEILOC(ILOC,158)=25
         NEILOC(ILOC,159)=25
         NEILOC(ILOC,160)=25
!     
         ILOC=17
!     
!     - face 20
!     
         NEILOC(ILOC,77)=8
         NEILOC(ILOC,78)=8
         NEILOC(ILOC,79)=8
         NEILOC(ILOC,80)=8
!     
!     - face 32
!     
         NEILOC(ILOC,125)=16
         NEILOC(ILOC,126)=16
         NEILOC(ILOC,127)=16
         NEILOC(ILOC,128)=16
!     
!     - face 30
!     
         NEILOC(ILOC,117)=14
         NEILOC(ILOC,118)=14
         NEILOC(ILOC,119)=14
         NEILOC(ILOC,120)=14
!     
!     - face 33
!     
         NEILOC(ILOC,129)=18
         NEILOC(ILOC,130)=18
         NEILOC(ILOC,131)=18
         NEILOC(ILOC,132)=18
!     
!     - face 41
!     
         NEILOC(ILOC,161)=26
         NEILOC(ILOC,162)=26
         NEILOC(ILOC,163)=26
         NEILOC(ILOC,164)=26
!     
         ILOC=18
!     
!     - face 21
!     
         NEILOC(ILOC,81)=9
         NEILOC(ILOC,82)=9
         NEILOC(ILOC,83)=9
         NEILOC(ILOC,84)=9
!     
!     - face 31
!     
         NEILOC(ILOC,121)=15
         NEILOC(ILOC,122)=15
         NEILOC(ILOC,123)=15
         NEILOC(ILOC,124)=15
!     
!     - face 33
!     
         NEILOC(ILOC,129)=17
         NEILOC(ILOC,130)=17
         NEILOC(ILOC,131)=17
         NEILOC(ILOC,132)=17
!     
!     - face 42
!     
         NEILOC(ILOC,165)=27
         NEILOC(ILOC,166)=27
         NEILOC(ILOC,167)=27
         NEILOC(ILOC,168)=27
!     
         ILOC=19
!     
!     - face 34
!     
         NEILOC(ILOC,133)=10
         NEILOC(ILOC,134)=10
         NEILOC(ILOC,135)=10
         NEILOC(ILOC,136)=10
!     
!     - face 43
!     
         NEILOC(ILOC,169)=20
         NEILOC(ILOC,170)=20
         NEILOC(ILOC,171)=20
         NEILOC(ILOC,172)=20
!     
!     - face 45
!     
         NEILOC(ILOC,177)=22
         NEILOC(ILOC,178)=22
         NEILOC(ILOC,179)=22
         NEILOC(ILOC,180)=22
!     
         ILOC=20
!     
!     - face 35
!     
         NEILOC(ILOC,137)=11
         NEILOC(ILOC,138)=11
         NEILOC(ILOC,139)=11
         NEILOC(ILOC,140)=11
!     
!     - face 43
!     
         NEILOC(ILOC,169)=19
         NEILOC(ILOC,170)=19
         NEILOC(ILOC,171)=19
         NEILOC(ILOC,172)=19
!     
!     - face 46
!     
         NEILOC(ILOC,181)=23
         NEILOC(ILOC,182)=23
         NEILOC(ILOC,183)=23
         NEILOC(ILOC,184)=23
!     
!     - face 44
!     
         NEILOC(ILOC,173)=21
         NEILOC(ILOC,174)=21
         NEILOC(ILOC,175)=21
         NEILOC(ILOC,176)=21
!     
         ILOC=21
!     
!     - face 36
!     
         NEILOC(ILOC,141)=12
         NEILOC(ILOC,142)=12
         NEILOC(ILOC,143)=12
         NEILOC(ILOC,144)=12
!     
!     - face 44
!     
         NEILOC(ILOC,173)=20
         NEILOC(ILOC,174)=20
         NEILOC(ILOC,175)=20
         NEILOC(ILOC,176)=20
!     
!     - face 47
!     
         NEILOC(ILOC,185)=24
         NEILOC(ILOC,186)=24
         NEILOC(ILOC,187)=24
         NEILOC(ILOC,188)=24
!     
         ILOC=22
!     
!     - face 37
!     
         NEILOC(ILOC,145)=13
         NEILOC(ILOC,146)=13
         NEILOC(ILOC,147)=13
         NEILOC(ILOC,148)=13
!     
!     - face 45
!     
         NEILOC(ILOC,177)=19
         NEILOC(ILOC,178)=19
         NEILOC(ILOC,179)=19
         NEILOC(ILOC,180)=19
!     
!     - face 48
!     
         NEILOC(ILOC,189)=23
         NEILOC(ILOC,190)=23
         NEILOC(ILOC,191)=23
         NEILOC(ILOC,192)=23
!     
!     - face 50
!     
         NEILOC(ILOC,197)=25
         NEILOC(ILOC,198)=25
         NEILOC(ILOC,199)=25
         NEILOC(ILOC,200)=25
!     
         ILOC=23
!     
!     - face 38
!     
         NEILOC(ILOC,149)=14
         NEILOC(ILOC,150)=14
         NEILOC(ILOC,151)=14
         NEILOC(ILOC,152)=14
!     
!     - face 46
!     
         NEILOC(ILOC,181)=20
         NEILOC(ILOC,182)=20
         NEILOC(ILOC,183)=20
         NEILOC(ILOC,184)=20
!     
!     - face 49
!     
         NEILOC(ILOC,193)=24
         NEILOC(ILOC,194)=24
         NEILOC(ILOC,195)=24
         NEILOC(ILOC,196)=24
!     
!     - face 51
!     
         NEILOC(ILOC,201)=26
         NEILOC(ILOC,202)=26
         NEILOC(ILOC,203)=26
         NEILOC(ILOC,204)=26
!     
!     - face 48
!     
         NEILOC(ILOC,189)=22
         NEILOC(ILOC,190)=22
         NEILOC(ILOC,191)=22
         NEILOC(ILOC,192)=22
!     
         ILOC=24
!     
!     - face 39
!     
         NEILOC(ILOC,153)=15
         NEILOC(ILOC,154)=15
         NEILOC(ILOC,155)=15
         NEILOC(ILOC,156)=15
!     
!     - face 47
!     
         NEILOC(ILOC,185)=21
         NEILOC(ILOC,186)=21
         NEILOC(ILOC,187)=21
         NEILOC(ILOC,188)=21
!     
!     - face 49
!     
         NEILOC(ILOC,193)=23
         NEILOC(ILOC,194)=23
         NEILOC(ILOC,195)=23
         NEILOC(ILOC,196)=23
!     
!     - face 52
!     
         NEILOC(ILOC,205)=27
         NEILOC(ILOC,206)=27
         NEILOC(ILOC,207)=27
         NEILOC(ILOC,208)=27
!     
         ILOC=25
!     
!     - face 40
!     
         NEILOC(ILOC,157)=16
         NEILOC(ILOC,158)=16
         NEILOC(ILOC,159)=16
         NEILOC(ILOC,160)=16
!     
!     - face 50
!     
         NEILOC(ILOC,197)=22
         NEILOC(ILOC,198)=22
         NEILOC(ILOC,199)=22
         NEILOC(ILOC,200)=22
!     
!     - face 53
!     
         NEILOC(ILOC,209)=26
         NEILOC(ILOC,210)=26
         NEILOC(ILOC,211)=26
         NEILOC(ILOC,212)=26
!     
         ILOC=26
!     
!     - face 41
!     
         NEILOC(ILOC,161)=17
         NEILOC(ILOC,162)=17
         NEILOC(ILOC,163)=17
         NEILOC(ILOC,164)=17
!     
!     - face 53
!     
         NEILOC(ILOC,209)=25
         NEILOC(ILOC,210)=25
         NEILOC(ILOC,211)=25
         NEILOC(ILOC,212)=25
!     
!     - face 51
!     
         NEILOC(ILOC,201)=23
         NEILOC(ILOC,202)=23
         NEILOC(ILOC,203)=23
         NEILOC(ILOC,204)=23
!     
!     - face 54
!     
         NEILOC(ILOC,213)=27
         NEILOC(ILOC,214)=27
         NEILOC(ILOC,215)=27
         NEILOC(ILOC,216)=27
!     
         ILOC=27
!     
!     - face 42
!     
         NEILOC(ILOC,165)=18
         NEILOC(ILOC,166)=18
         NEILOC(ILOC,167)=18
         NEILOC(ILOC,168)=18
!     
!     - face 52
!     
         NEILOC(ILOC,205)=24
         NEILOC(ILOC,206)=24
         NEILOC(ILOC,207)=24
         NEILOC(ILOC,208)=24
!     
!     - face 54
!     
         NEILOC(ILOC,213)=26
         NEILOC(ILOC,214)=26
         NEILOC(ILOC,215)=26
         NEILOC(ILOC,216)=26
!     
      ELSE IF(ELETYP.EQ.7) THEN
!     
         ILOC=1
!     
!     - face 1
!     
         NEILOC(ILOC,1) = 2
         NEILOC(ILOC,2) = 2
!     
!     - face 2
!     
         NEILOC(ILOC,3) = 6
         NEILOC(ILOC,4) = 6
!     
         ILOC=2
!     
!     - face 1
!     
         NEILOC(ILOC,1) = 1
         NEILOC(ILOC,2) = 1
!     
!     - face 7
!     
         NEILOC(ILOC,13) = 6
         NEILOC(ILOC,14) = 6
!     
!     - face 8
!     
         NEILOC(ILOC,15) = 4
         NEILOC(ILOC,16) = 4
!     
!     - face 3
!     
         NEILOC(ILOC,5) = 3
         NEILOC(ILOC,6) = 3
!     
         ILOC=3
!     
!     - face 3
!     
         NEILOC(ILOC,5) = 2
         NEILOC(ILOC,6) = 2
!     
!     - face 3
!     
         NEILOC(ILOC,7) = 4
         NEILOC(ILOC,8) = 4
!     
         ILOC=4
!     
!     - face 4
!     
         NEILOC(ILOC,7) = 3
         NEILOC(ILOC,8) = 3
!     
!     - face 8
!     
         NEILOC(ILOC,15) = 2
         NEILOC(ILOC,16) = 2
!     
!     - face 9
!     
         NEILOC(ILOC,17) = 6
         NEILOC(ILOC,18) = 6
!     
!     - face 5
!     
         NEILOC(ILOC,9)  = 5
         NEILOC(ILOC,10) = 5
!     
         ILOC=5
!     
!     - face 5
!     
         NEILOC(ILOC,9)  = 4  
         NEILOC(ILOC,10) = 4
!     
!     - face 6
!     
         NEILOC(ILOC,11) = 6 
         NEILOC(ILOC,12) = 6
!     
         ILOC=6
!     
!     - face 6
!     
         NEILOC(ILOC,11) = 5
         NEILOC(ILOC,12) = 5
!     
!     - face 9
!     
         NEILOC(ILOC,17) = 4
         NEILOC(ILOC,18) = 4
!     
!     - face 7
!     
         NEILOC(ILOC,13) = 2 
         NEILOC(ILOC,14) = 2
!     
!     - face 2
!     
         NEILOC(ILOC,3) = 1
         NEILOC(ILOC,4) = 1
!     
      ELSE IF(ELETYP.EQ.8) THEN
!     
         ILOC=1
!     
!     - face 1
!     
         NEILOC(ILOC,1) = 6
         NEILOC(ILOC,2) = 6
         NEILOC(ILOC,3) = 6
         NEILOC(ILOC,4) = 6
!     
!     - face 2
!     
         NEILOC(ILOC,5) = 2
         NEILOC(ILOC,6) = 2
         NEILOC(ILOC,7) = 2
         NEILOC(ILOC,8) = 2
!     
!     - face 3
!     
         NEILOC(ILOC,9)  = 7
         NEILOC(ILOC,10) = 7
         NEILOC(ILOC,11) = 7
         NEILOC(ILOC,12) = 7
!     
         ILOC=2
!     
!     - face 2
!     
         NEILOC(ILOC,5) = 1
         NEILOC(ILOC,6) = 1
         NEILOC(ILOC,7) = 1
         NEILOC(ILOC,8) = 1
!     
!     - face 5
!     
         NEILOC(ILOC,17) = 3
         NEILOC(ILOC,18) = 3
         NEILOC(ILOC,19) = 3
         NEILOC(ILOC,20) = 3
!     
!     - face 20
!     
         NEILOC(ILOC,77) = 8
         NEILOC(ILOC,78) = 8
         NEILOC(ILOC,79) = 8
         NEILOC(ILOC,80) = 8
!     
!     - face 21
!     
         NEILOC(ILOC,81) = 7
         NEILOC(ILOC,82) = 7
         NEILOC(ILOC,83) = 7
         NEILOC(ILOC,84) = 7
!     
!     - face 23
!     
         NEILOC(ILOC,89) = 6
         NEILOC(ILOC,80) = 6
         NEILOC(ILOC,91) = 6
         NEILOC(ILOC,92) = 6
!     
!     - face 24
!     
         NEILOC(ILOC,93) = 4
         NEILOC(ILOC,94) = 4
         NEILOC(ILOC,95) = 4
         NEILOC(ILOC,96) = 4
!     
         ILOC=3
!     
!     - face 4
!     
         NEILOC(ILOC,13) = 4
         NEILOC(ILOC,14) = 4
         NEILOC(ILOC,15) = 4
         NEILOC(ILOC,16) = 4
!     
!     - face 5
!     
         NEILOC(ILOC,17) = 2
         NEILOC(ILOC,18) = 2
         NEILOC(ILOC,19) = 2
         NEILOC(ILOC,20) = 2
!     
!     - face 6
!     
         NEILOC(ILOC,21) = 8
         NEILOC(ILOC,22) = 8
         NEILOC(ILOC,23) = 8
         NEILOC(ILOC,24) = 8
!     
         ILOC=4
!     
!     - face 4 
!     
         NEILOC(ILOC,13) = 3
         NEILOC(ILOC,14) = 3
         NEILOC(ILOC,15) = 3
         NEILOC(ILOC,16) = 3
!     
!     - face 8
!     
         NEILOC(ILOC,29) = 5
         NEILOC(ILOC,30) = 5
         NEILOC(ILOC,31) = 5
         NEILOC(ILOC,32) = 5
!     
!     - face 18
!     
         NEILOC(ILOC,69) = 9
         NEILOC(ILOC,70) = 9
         NEILOC(ILOC,71) = 9
         NEILOC(ILOC,72) = 9
!     
!     - face 19
!     
         NEILOC(ILOC,73) = 8
         NEILOC(ILOC,74) = 8
         NEILOC(ILOC,75) = 8
         NEILOC(ILOC,76) = 8
!     
!     - face 22
!     
         NEILOC(ILOC,85) = 6
         NEILOC(ILOC,86) = 6
         NEILOC(ILOC,87) = 6
         NEILOC(ILOC,88) = 6
!     
!     - face 24
!     
         NEILOC(ILOC,93) = 2
         NEILOC(ILOC,94) = 2
         NEILOC(ILOC,95) = 2
         NEILOC(ILOC,96) = 2
!     
         ILOC=5
!     
!     - face 7
!     
         NEILOC(ILOC,25) = 6
         NEILOC(ILOC,26) = 6
         NEILOC(ILOC,27) = 6
         NEILOC(ILOC,28) = 6
!     
!     - face 8
!     
         NEILOC(ILOC,29) = 4
         NEILOC(ILOC,30) = 4
         NEILOC(ILOC,31) = 4
         NEILOC(ILOC,32) = 4
!     
!     - face 9
!     
         NEILOC(ILOC,33) = 9
         NEILOC(ILOC,34) = 9
         NEILOC(ILOC,35) = 9
         NEILOC(ILOC,36) = 9
!     
         ILOC=6
!     
!     - face 1
!     
         NEILOC(ILOC,1) = 1
         NEILOC(ILOC,2) = 1
         NEILOC(ILOC,3) = 1
         NEILOC(ILOC,4) = 1
!     
!     - face 7
!     
         NEILOC(ILOC,25) = 5
         NEILOC(ILOC,26) = 5
         NEILOC(ILOC,27) = 5
         NEILOC(ILOC,28) = 5
!     
!     - face 16
!     
         NEILOC(ILOC,61) = 7
         NEILOC(ILOC,62) = 7
         NEILOC(ILOC,63) = 7
         NEILOC(ILOC,64) = 7
!     
!     - face 17
!     
         NEILOC(ILOC,65) = 9
         NEILOC(ILOC,66) = 9
         NEILOC(ILOC,67) = 9
         NEILOC(ILOC,68) = 9
!     
!     - face 22
!     
         NEILOC(ILOC,85) = 4
         NEILOC(ILOC,86) = 4
         NEILOC(ILOC,87) = 4
         NEILOC(ILOC,88) = 4
!     
!     - face 23
!     
         NEILOC(ILOC,89) = 2
         NEILOC(ILOC,90) = 2
         NEILOC(ILOC,91) = 2
         NEILOC(ILOC,92) = 2
!     
         ILOC=7
!     
!     - face 3
!     
         NEILOC(ILOC,9)  = 1
         NEILOC(ILOC,10) = 1
         NEILOC(ILOC,11) = 1
         NEILOC(ILOC,12) = 1
!     
!     - face 11
!     
         NEILOC(ILOC,41) = 10
         NEILOC(ILOC,42) = 10
         NEILOC(ILOC,43) = 10
         NEILOC(ILOC,44) = 10
!     
!     - face 13
!     
         NEILOC(ILOC,49) = 9
         NEILOC(ILOC,50) = 9
         NEILOC(ILOC,51) = 9
         NEILOC(ILOC,52) = 9
!     
!     - face 15
!     
         NEILOC(ILOC,57) = 8
         NEILOC(ILOC,58) = 8
         NEILOC(ILOC,59) = 8
         NEILOC(ILOC,60) = 8
!     
!     - face 16
!     
         NEILOC(ILOC,61) = 6
         NEILOC(ILOC,62) = 6
         NEILOC(ILOC,63) = 6
         NEILOC(ILOC,64) = 6
!     
!     - face 21
!     
         NEILOC(ILOC,81) = 2
         NEILOC(ILOC,82) = 2
         NEILOC(ILOC,83) = 2
         NEILOC(ILOC,84) = 2
!     
         ILOC=8
!     
!     - face 6
!     
         NEILOC(ILOC,21) = 3
         NEILOC(ILOC,22) = 3
         NEILOC(ILOC,23) = 3
         NEILOC(ILOC,24) = 3
!     
!     - face 12
!     
         NEILOC(ILOC,45) = 10
         NEILOC(ILOC,46) = 10
         NEILOC(ILOC,47) = 10
         NEILOC(ILOC,48) = 10
!     
!     - face 14
!     
         NEILOC(ILOC,53) = 9
         NEILOC(ILOC,54) = 9
         NEILOC(ILOC,55) = 9
         NEILOC(ILOC,56) = 9
!     
!     - face 15
!     
         NEILOC(ILOC,57) = 7
         NEILOC(ILOC,58) = 7
         NEILOC(ILOC,59) = 7
         NEILOC(ILOC,60) = 7
!     
!     - face 19
!     
         NEILOC(ILOC,73) = 4
         NEILOC(ILOC,74) = 4
         NEILOC(ILOC,75) = 4
         NEILOC(ILOC,76) = 4
!     
!     - face 20
!     
         NEILOC(ILOC,77) = 2
         NEILOC(ILOC,78) = 2
         NEILOC(ILOC,79) = 2
         NEILOC(ILOC,80) = 2
!     
         ILOC=9
!     
!     - face 9
!     
         NEILOC(ILOC,33) = 5
         NEILOC(ILOC,34) = 5
         NEILOC(ILOC,35) = 5
         NEILOC(ILOC,36) = 5
!     
!     - face 10
!     
         NEILOC(ILOC,37) = 10
         NEILOC(ILOC,38) = 10
         NEILOC(ILOC,39) = 10
         NEILOC(ILOC,40) = 10
!     
!     - face 13
!     
         NEILOC(ILOC,49) = 7
         NEILOC(ILOC,50) = 7
         NEILOC(ILOC,51) = 7
         NEILOC(ILOC,52) = 7
!     
!     - face 14
!     
         NEILOC(ILOC,53) = 8
         NEILOC(ILOC,54) = 8
         NEILOC(ILOC,55) = 8
         NEILOC(ILOC,56) = 8
!     
!     - face 17
!     
         NEILOC(ILOC,65) = 6
         NEILOC(ILOC,66) = 6
         NEILOC(ILOC,67) = 6
         NEILOC(ILOC,68) = 6
!     
!     - face 18
!     
         NEILOC(ILOC,69) = 4
         NEILOC(ILOC,70) = 4
         NEILOC(ILOC,71) = 4
         NEILOC(ILOC,72) = 4
!     
         ILOC=10
!     
!     - face 10
!     
         NEILOC(ILOC,37) = 9
         NEILOC(ILOC,38) = 9
         NEILOC(ILOC,39) = 9
         NEILOC(ILOC,40) = 9
!     
!     - face 11
!     
         NEILOC(ILOC,41) = 7
         NEILOC(ILOC,42) = 7
         NEILOC(ILOC,43) = 7
         NEILOC(ILOC,44) = 7
!     
!     - face 12
!     
         NEILOC(ILOC,45) = 8
         NEILOC(ILOC,46) = 8
         NEILOC(ILOC,47) = 8
         NEILOC(ILOC,48) = 8
!     
      ENDIF
!     END SELECT
!     
      END SUBROUTINE

      SUBROUTINE SURNEI( SNEILOC, SNLOC, SNGI, SELETYP )
!     ------------------------------------------------
!     
! - this subroutine calculates SNEILOC which is the array containing
! - information if SGI is in the same CV as SILOC then set
! - SNEILOC(SILOC,SGI) = 1 else SNEILOC(SILOC,SGI) = 0. Note that
! - SNEILOC is related to the FV basis functions M(ILOC,GI) for the
! - 2-D elements but instead of 1.0 and 0.0 we have integer values.
!     
!     -------------------------------
!     - date last modified : 23/10/2006
!     - added by crgw (taken from RADIANT)
!     -------------------------------
!     
      IMPLICIT NONE
!     
      INTEGER SELETYP, SNLOC, SNGI
!     
      INTEGER SNEILOC(SNLOC,SNGI)
!     
!     - local variables
!     
      INTEGER COUNT, INC, SILOC, SGI
!     
      SNEILOC = 0
!     
!     
      IF(SELETYP.EQ.1) THEN
!
! - linear line surface elements
!
          COUNT = 1
!
      do SILOC = 1, SNLOC! Was loop 
!
      do SGI = COUNT, COUNT! Was loop 
!
                SNEILOC(SILOC,SGI) = 1
!
             END DO
!
             COUNT = COUNT + 1
!
          END DO
!
! - This is a little note for the future if we are
! - using 16 point quadrature for the surface elements
! - then change the loop DO SGI = COUNT, COUNT to
! - DO SGI = COUNT, COUNT + 3. Also change the loop
! - increment to COUNT = COUNT + 4. This is similar for
! - triangles.
!     
      ELSE IF(SELETYP.EQ.2) THEN
!
! - linear quadrilateral surface elements
! 
          COUNT = 1
!
      do SILOC = 1, SNLOC! Was loop 
!
      do SGI = COUNT, COUNT! Was loop 
!
                SNEILOC(SILOC,SGI) = 1
!
            END DO
!
            COUNT = COUNT + 1
!
         END DO
!     
      ELSE IF(SELETYP.EQ.3) THEN
!
! - linear triangular surface elements
!
          COUNT = 1
!
      do SILOC = 1, SNLOC! Was loop 
!
      do SGI = COUNT, COUNT! Was loop 
!
                SNEILOC(SILOC,SGI) = 1
!
             END DO
!
             COUNT = COUNT + 1
!
          END DO    
!     
      ELSE IF(SELETYP.EQ.4) THEN
!
! - quadratic line surface elements
!
          COUNT = 1
!
      do SILOC = 1, SNLOC! Was loop 
!
      do SGI = COUNT, COUNT+1! Was loop 
!
                SNEILOC(SILOC,SGI) = 1
!
             END DO
!
             COUNT = COUNT + 2    
!
          END DO     
!     
      ELSE IF(SELETYP.EQ.5) THEN
!
! - quadratic quadrilateral surface elements
!
          COUNT = 1
!
      do SILOC = 1, SNLOC! Was loop 
!
      do SGI = COUNT, COUNT + 3! Was loop 
!
                SNEILOC(SILOC,SGI) = 1
!
             END DO
!
             COUNT = COUNT + 4
!
          END DO
!     
      ELSE IF(SELETYP.EQ.6) THEN
!
! - quadratic triangular surface elements
!
           COUNT = 1
!
      do SILOC = 1, SNLOC! Was loop 
!
              IF ((SILOC.EQ.1).OR.(SILOC.EQ.3).OR.(SILOC.EQ.5)) THEN
!
                  INC = 3
!
              ELSE
!
                  INC = 7
!
              ENDIF
!
      do SGI = COUNT, COUNT+INC! Was loop 
!
                 SNEILOC(SILOC,SGI) = 1
!
              END DO
!
              IF ((SILOC.EQ.1).OR.(SILOC.EQ.3).OR.(SILOC.EQ.5)) THEN
!
                  COUNT = COUNT + 4
!
              ELSE
!
                  COUNT = COUNT + 8
!
              ENDIF
!
           END DO    
!     
! - note that in future that the DO GI = COUNT, COUNT+3
! - will become DO GI = COUNT, COUNT+NGI/NLOC
!
! - note the use of NGI/NLOC which gives the number of Gauss
! - points associated with a CV. We might consider the use
! - of another variable NCVGI which gives the number of Gauss
! - points associated with a CV.
!     
      ENDIF
!     END SELECT
!     
      END SUBROUTINE


end module legacy_cv_numbering
