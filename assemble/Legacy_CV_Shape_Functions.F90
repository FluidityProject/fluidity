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
module legacy_cv_shape_functions
  !!< Generate shape functions for elements of arbitrary polynomial degree.
  use FLDebug
  use legacy_cv_numbering
  implicit none

   private
   public :: SHAPESV, SHAPESE
  
contains


      SUBROUTINE SHAPESV( ELETYP,  NEILOC,    &
!     - INTEGERS
     &     NGI,     NLOC,   &
     &     SVNGI,&
!     - REALS
     &     SVN,&
     &     SVNLX,   SVNLY,  &
     &     SVWEIGH    &
     &     ,M&
     &     )
!     ----------------------------------------------
!     
!     - this subroutine generates the FE basis functions, weights and the
!     - derivatives of the shape functions for a variety of elements. 
!     - The routine also generates the shape functions and derivatives 
!     - associated with the CV surfaces and also the FV basis functions.
!     
!     -------------------------------
!     - date last modified : 24/05/2003
!     -------------------------------
!     
      IMPLICIT NONE
!     
      INTEGER NGI, NLOC, SVNGI
!     
      INTEGER ELETYP
!     
      INTEGER NEILOC(NLOC,SVNGI)
!     
      REAL SVN(NLOC,SVNGI),   SVNLX(NLOC,SVNGI)
      REAL SVNLY(NLOC,SVNGI)
!     
      REAL SVWEIGH(SVNGI)
!     Local variable...
!       REAL M(1000)
      REAL M(NLOC,NGI)
!     
!     
      IF(ELETYP.EQ.1) THEN
!     
         CALL FVQUAD( NGI,    NLOC, SVNGI,&
     &        M,      SVN,  SVNLX,&
     &        SVWEIGH              )
!     
      ELSE IF(ELETYP.EQ.2) THEN
!     
         CALL FVTRI( NGI,    NLOC, SVNGI,&
     &        M,      SVN,  SVNLX,&
     &        SVWEIGH              )
!     
      ELSE IF(ELETYP.EQ.3) THEN
!     
         CALL FVHEX( NGI,    NLOC,    SVNGI,&
     &        M,      SVN,     SVNLX, &
     &        SVNLY,  SVWEIGH         )
!     
      ELSE IF(ELETYP.EQ.4) THEN
!     
         CALL FVTET( NGI,    NLOC,    SVNGI,&
     &        M,      SVN,     SVNLX, &
     &        SVNLY,  SVWEIGH         )
!     
      ELSE IF(ELETYP.EQ.5) THEN
!     
         CALL FVQQUAD(  NGI,    NLOC,    SVNGI,&
     &        M,      SVN,     SVNLX, &
     &        SVWEIGH                 )
!     
      ELSE IF(ELETYP.EQ.6) THEN
!     
         CALL FVQHEX( NGI,   NLOC,   SVNGI,&
     &        M,     SVN,    SVNLX,&
     &        SVNLY, SVWEIGH        )
!     
      ELSE IF(ELETYP.EQ.7) THEN
!     
!     - quadratic triangles
!     
         CALL FVQTRI( NGI,    NLOC, SVNGI,&
     &        M,      SVN,  SVNLX,  &
     &        SVWEIGH              )
!     
      ELSE IF(ELETYP.EQ.8) THEN
!     
         CALL FVQTET( NGI,   NLOC,   SVNGI,&
     &        M,     SVN,    SVNLX,&
     &        SVNLY, SVWEIGH        )
!     
      END IF
!     
!     
      CALL VOLNEI( NEILOC, NLOC, SVNGI, ELETYP )
!     
      RETURN
      END SUBROUTINE

      SUBROUTINE SHAPESE( SELETYP,  SNEILOC,    &
!     - INTEGERS
     &     SNGI,     SNLOC,   &
!     - REALS
     &     SN,&
     &     SNLX,   SNLY,  &
     &     SWEIGH    &
     &     )
!     ----------------------------------------------
!     
!     - this subroutine generates the FE basis functions, weights and the
!     - derivatives of the surface shape functions for a variety of elements. 
!     
!     -------------------------------
!     - date last modified : 24/05/2003
!     -------------------------------
!     
      use fldebug
      IMPLICIT NONE
!     
      INTEGER SNGI, SNLOC
!     
      INTEGER SELETYP
!     
      INTEGER SNEILOC(SNLOC,SNGI)
!     
      REAL SN(SNLOC,SNGI),   SNLX(SNLOC,SNGI)
      REAL SNLY(SNLOC,SNGI)
!     
      REAL SWEIGH(SNGI)
!     
      IF(SELETYP.EQ.1) THEN
!     
          CALL SFELIN( SNGI, SNLOC, SN, SNLX, SWEIGH )
!     
      ELSE IF(SELETYP.EQ.2) THEN
!     
          CALL SFEQUAD( SNGI,   SNLOC,&
     &                  SN,     SNLX,  &
     &                  SNLY,   SWEIGH )
!     
      ELSE IF(SELETYP.EQ.3) THEN
!     
          CALL SFETRI( SNGI,  SNLOC,&
     &                 SN,    SNLX,  SNLY,&
     &                 SWEIGH             )
!     
      ELSE IF(SELETYP.EQ.4) THEN
!     
          CALL SFEQLIN( SNGI, SNLOC, SN, SNLX, SWEIGH )
!     
      ELSE IF(SELETYP.EQ.5) THEN
!     
          CALL SFEQQUAD( SNGI, SNLOC,&
     &                   SN,   SNLX,&
     &                   SNLY, SWEIGH )
!     
      ELSE IF(SELETYP.EQ.6) THEN
!     
          CALL SFEQTRI( SNGI, SNLOC,&
     &                  SN,   SNLX,&
     &                  SNLY, SWEIGH )
!     
      END IF
!     
!     
      CALL SURNEI( SNEILOC, SNLOC, SNGI, SELETYP )
! 

      RETURN
      END SUBROUTINE


! ELETYPS

      SUBROUTINE FVQUAD( NGI,    NLOC,    SVNGI,&
!     - REALS
     &     M,      SVN,     SVNLX, &
     &     SVWEIGH                 )
!     --------------------------------------------
!     
!     - this routine generates the shape functions associated
!     - with the FV's i.e. their surfaces and volume shape 
!     - functions and derivatives. The surface shape functions
!     - are the values of the FE volume shape functions evaluated
!     - on the surfaces of the CV's.
!     
!     -------------------------------
!     - date last modified : 07/11/2002
!     -------------------------------
!     
      IMPLICIT NONE
!     
      INTEGER NGI, NLOC, SVNGI
!     
      REAL M(NLOC,NGI),     SVWEIGH(SVNGI)
!     
      REAL SVN(NLOC,SVNGI), SVNLX(NLOC,SVNGI)
!     
!     - local variables
!     
      INTEGER ILOC, JLOC, IFACE, ICOORD
!     
      INTEGER GI, GJ
!     
!     - note that NCOORD is the number of co-ordinates.
!     - NFACE is the number of surfaces internal to an 
!     - element 
!     
      INTEGER NCOORD, NFACE
      PARAMETER( NCOORD = 2, NFACE = 4 )
!     
      INTEGER ENLOC, FNGI, FNLOC
      PARAMETER( ENLOC = 4, FNGI = 1, FNLOC = 2  )
!     
      REAL POS(NCOORD),    DPDXI(NCOORD)
      REAL DPDETA(NCOORD)
!     
      REAL XIGP(FNGI)
!     
      REAL XIP(FNLOC)
!     
      REAL XI(ENLOC), ETA(ENLOC)
!     
      REAL CORN(NFACE,NCOORD,FNLOC)
!     
!     - FV basis functions for use in calculating the integral of the
!     - shape functions over a subcell of a CV.
!     
      M = 0.0
!     
!     COUNT = 1
!     
!     DO ILOC = 1, NLOC
!     
!     DO GI = COUNT, COUNT+3
!     
!     M(ILOC,GI) = 1.0
!     
!     END DO
!     
!     COUNT = COUNT + 4
!     
!     END DO
!     
!     - Note that CORN(IFACE,ICOORD,JLOC) tabulates the positions
!     - of the corners (vertices) of the faces of a subcell in volume 
!     - co-ordinates. Note also that IFACE ranges from 1-4, ICOORD
!     - ranges from 1-2 and JLOC ranges from 1-2. IFACE signifies 
!     - the particular subcell face, ICOORD signfies the co-ordinates
!     - of the vertices of the face in (XI,ETA) co-ordinates
!     - and JLOC signifies the vertex of the face of the subcell
!     - which are straight lines.
!     
!     - face 1
!     
      CORN(1,1,1) = -1.0
      CORN(1,2,1) =  0.0
!     
      CORN(1,1,2) =  0.0
      CORN(1,2,2) =  0.0
!     
!     - face 2
!     
      CORN(2,1,1) =  0.0
      CORN(2,2,1) =  0.0
!     
      CORN(2,1,2) =  1.0
      CORN(2,2,2) =  0.0
!     
!     - face 3
!     
      CORN(3,1,1) =  0.0
      CORN(3,2,1) = -1.0
!     
      CORN(3,1,2) =  0.0
      CORN(3,2,2) =  0.0
!     
!     - face 4
!     
      CORN(4,1,1) =  0.0
      CORN(4,2,1) =  0.0
!     
      CORN(4,1,2) =  0.0
      CORN(4,2,2) =  1.0
!     
!     - define the Gaussian integration points in XIP space.
!     
      XIGP(1)  = 0.0
!     
!     - define positions of vertices of line element
!     - associated with a subcell face in XIP
!     
      XIP(1) = -1.0
      XIP(2) =  1.0
!     
!     - define positions of vertices of quadrilateral element in 
!     - (XI,ETA) space.
!     
      XI(1)  = -1.0
      ETA(1) = -1.0
!     
      XI(2)  =  1.0
      ETA(2) = -1.0
!     
      XI(3)  = -1.0
      ETA(3) =  1.0
!     
      XI(4)  =  1.0
      ETA(4) =  1.0
!     
!     - generate values of the FE basis functions at the quadrature
!     - points on the faces of the subcells.
!     
      do ILOC = 1, NLOC! Was loop 
!     
      do IFACE = 1, NFACE! Was loop 
!     
      do GJ = 1, FNGI! Was loop 
!     
               GI = IFACE
!     
      do ICOORD = 1, NCOORD! Was loop 
!     
                  POS(ICOORD)    = 0.0
                  DPDXI(ICOORD)  = 0.0
                  DPDETA(ICOORD) = 0.0
!     
      do JLOC = 1, FNLOC! Was loop 
!     
                     POS(ICOORD) = POS(ICOORD)&
     &                    + CORN(IFACE,ICOORD,JLOC)&
     &                    *0.5*(1.+XIP(JLOC)*XIGP(GJ))
!     
                     DPDXI(ICOORD) = DPDXI(ICOORD)&
     &                    + CORN(IFACE,ICOORD,JLOC)&
     &                    *0.5*XIP(JLOC)
!     
                  END DO
!     
               END DO
!     
               SVN(ILOC,GI) = 0.25*( 1.0 + XI(ILOC)*POS(1) )&
     &              *( 1.0 + ETA(ILOC)*POS(2) )
!     
               SVNLX(ILOC,GI) = 0.25*XI(ILOC)*DPDXI(1)&
     &              *( 1.0 + ETA(ILOC)*POS(2) )&
     &              + 0.25*( 1.0 + XI(ILOC)*POS(1) )&
     &              *ETA(ILOC)*DPDXI(2)
!     
            END DO
!     
         END DO
!     
      END DO 
!     
!     - set the weights for the surface integration
!     
      do GI = 1, SVNGI! Was loop 
!     
         SVWEIGH(GI) = 2.0
!     
      END DO
!     
      END SUBROUTINE
!     
!     
!     
!     
      SUBROUTINE FVTRI( NGI,    NLOC, SVNGI,&
     &     M,      SVN,  SVNLX,  &
     &     SVWEIGH              )
!     ----------------------------------------
!     
!     - this routine generates the shape functions associated
!     - with the FV's i.e. their surfaces and volume shape 
!     - functions and derivatives.
!     
!     -------------------------------
!     - date last modified : 06/10/2002
!     -------------------------------
!     
      IMPLICIT NONE
!     
      INTEGER NGI, NLOC, SVNGI
!     
      REAL M(NLOC,NGI),     SVN(NLOC,NGI)
      REAL SVNLX(NLOC,NGI), SVWEIGH(SVNGI)
!     
!     - local variables
!     
      INTEGER ILOC, JLOC, IFACE, ICOORD
!     
      INTEGER GI, GJ
!     
!     - note that NCOORD is the number of co-ordinates.
!     - NFACE is the number of surfaces internal to an 
!     - element 
!     
      INTEGER NCOORD, NFACE
      PARAMETER( NCOORD = 3, NFACE = 3 )
!     
      INTEGER FNGI, FNLOC
      PARAMETER( FNGI = 1, FNLOC = 2 )
!     
      REAL POS(NCOORD), DPDXI(NCOORD)
!     
      REAL XIGP(FNGI)
!     
      REAL XI(FNLOC)
!     
      REAL CORN(NFACE,NCOORD,FNLOC)
!     
!     - FV basis functions for use in calculating the integral of the
!     - shape functions over a subcell of a CV.
!     
      M = 0.0
!     
!     COUNT = 1
!     
!     DO ILOC = 1, NLOC
!     
!     DO GI = COUNT, COUNT+3
!     
!     M(ILOC,GI) = 1.0
!     
!     END DO
!     
!     COUNT = COUNT + 4
!     
!     END DO      
!     
!     - Note that CORN(IFACE,ICOORD,JLOC) tabulates the positions
!     - of the corners (vertices) of the faces of a subcell in volume 
!     - co-ordinates. Note also that IFACE ranges from 1-3, ICOORD
!     - ranges from 1-3 and JLOC ranges from 1-2. IFACE signifies 
!     - the particular subcell face, ICOORD signfies the co-ordinates
!     - of the vertices of the face in volume co-ordinates (L1,L2,L3)
!     - and JLOC signifies the vertex of the face of the subcell
!     - which are straight lines.
!     
!     - face 1
!     
      CORN(1,1,1) = 0.33333333
      CORN(1,2,1) = 0.33333333
      CORN(1,3,1) = 0.33333333
!     
      CORN(1,1,2) = 0.5
      CORN(1,2,2) = 0.5
      CORN(1,3,2) = 0.0
!     
!     - face 2
!     
      CORN(2,1,1) = 0.33333333
      CORN(2,2,1) = 0.33333333
      CORN(2,3,1) = 0.33333333
!     
      CORN(2,1,2) = 0.0
      CORN(2,2,2) = 0.5
      CORN(2,3,2) = 0.5 
!     
!     - face 3
!     
      CORN(3,1,1) = 0.33333333
      CORN(3,2,1) = 0.33333333
      CORN(3,3,1) = 0.33333333
!     
      CORN(3,1,2) = 0.5
      CORN(3,2,2) = 0.0
      CORN(3,3,2) = 0.5 
!     
!     - define the Gaussian integration points in XI space.
!     
      XIGP(1) = 0.0
!     
!     - define positions of end points of line element 
!     - associated with a subcell face in XI space.
!     
      XI(1) = -1.0
      XI(2) =  1.0
!     
!     - generate values of the FE basis functions at the quadrature
!     - points on the faces of the subcells.
!     
      do ILOC = 1, NLOC! Was loop 
!     
      do IFACE = 1, NFACE! Was loop 
!     
      do GJ = 1, FNGI! Was loop 
!     
               GI = IFACE
!     
      do ICOORD = 1, NCOORD! Was loop 
!     
                  POS(ICOORD)    = 0.0
                  DPDXI(ICOORD)  = 0.0
!     
      do JLOC = 1, FNLOC! Was loop 
!     
                     POS(ICOORD) = POS(ICOORD)&
     &                    + CORN(IFACE,ICOORD,JLOC)&
     &                    *0.5*(1.+XI(JLOC)*XIGP(GJ))
!     
                     DPDXI(ICOORD) = DPDXI(ICOORD)&
     &                    + CORN(IFACE,ICOORD,JLOC)&
     &                    *0.5*XI(JLOC)
!     
                  END DO
!     
               END DO
!     
               SVN(ILOC,GI)   = POS(ILOC)
               SVNLX(ILOC,GI) = DPDXI(ILOC)
!     
            END DO
!     
         END DO
!     
      END DO    
!     
!     - set the weights for the surface integration
!     
      do GI = 1, SVNGI! Was loop 
!     
         SVWEIGH(GI) = 2.0
!     
      END DO
!     
      END SUBROUTINE

      SUBROUTINE FVHEX( NGI,    NLOC,    SVNGI,&
!     - REALS
     &     M,      SVN,     SVNLX, &
     &     SVNLY,  SVWEIGH          )
!     --------------------------------------------
!     
!     - this routine generates the shape functions associated
!     - with the FV's i.e. their surfaces and volume shape 
!     - functions and derivatives.
!     
!     -------------------------------
!     - date last modified : 06/10/2002
!     -------------------------------
!     
      IMPLICIT NONE
!     
      INTEGER NGI, NLOC, SVNGI
!     
      REAL M(NLOC,NGI), SVN(NLOC,SVNGI)
!     
      REAL SVNLX(NLOC,SVNGI), SVNLY(NLOC,SVNGI)
      REAL SVWEIGH(SVNGI)
!     
!     - local variables
!     
      INTEGER ILOC, JLOC, IFACE, ICOORD
!     
      INTEGER COUNT, GI, GJ
!     
!     - note that NCOORD is the number of co-ordinates.
!     - NFACE is the number of surfaces internal to an 
!     - element 
!     
      INTEGER NCOORD, NFACE
      PARAMETER( NCOORD = 3, NFACE = 12 )
!     
      INTEGER ENLOC, FNGI, FNLOC
      PARAMETER( ENLOC = 8, FNGI = 1, FNLOC = 4 )
!     
      REAL POS(NCOORD),    DPDXI(NCOORD)
      REAL DPDETA(NCOORD)
!     
      REAL XIGP(FNGI), ETAGP(FNGI)
!     
      REAL XIP(FNLOC), ETAP(FNLOC)
!     
      REAL XI(ENLOC), ETA(ENLOC), ZETA(ENLOC)
!     
      REAL CORN(NFACE,NCOORD,FNLOC)
!     
!     - FV basis functions for use in calculating the integral of the
!     - shape functions over a subcell of a CV.
!     
      M = 0.0
!     
      COUNT = 1
!     
!     DO ILOC = 1, NLOC
!     
!     DO GI = COUNT, COUNT+7
!     
!     M(ILOC,GI) = 1.0
!     
!     END DO
!     
!     COUNT = COUNT + 8
!     
!     END DO    
!     
!     - Note that CORN(IFACE,ICOORD,JLOC) tabulates the positions
!     - of the corners (vertices) of the faces of a subcell in volume 
!     - co-ordinates. Note also that IFACE ranges from 1-12, ICOORD
!     - ranges from 1-3 and JLOC ranges from 1-4. IFACE signifies 
!     - the particular subcell face, ICOORD signfies the co-ordinates
!     - of the vertices of the face in (XI,ETA,ZETA) co-ordinates
!     - and JLOC signifies the vertex of the face of the subcell
!     - which are quadrilaterals.
!     
!     - face 1
!     
      CORN(1,1,1) = -1.0
      CORN(1,2,1) =  0.0
      CORN(1,3,1) = -1.0
!     
      CORN(1,1,2) =  0.0
      CORN(1,2,2) =  0.0
      CORN(1,3,2) = -1.0
!     
      CORN(1,1,3) = -1.0
      CORN(1,2,3) =  0.0
      CORN(1,3,3) =  0.0
!     
      CORN(1,1,4) =  0.0 
      CORN(1,2,4) =  0.0
      CORN(1,3,4) =  0.0
!     
!     - face 2
!     
      CORN(2,1,1) =  0.0
      CORN(2,2,1) =  0.0
      CORN(2,3,1) = -1.0
!     
      CORN(2,1,2) =  1.0
      CORN(2,2,2) =  0.0
      CORN(2,3,2) = -1.0
!     
      CORN(2,1,3) =  0.0
      CORN(2,2,3) =  0.0
      CORN(2,3,3) =  0.0
!     
      CORN(2,1,4) =  1.0
      CORN(2,2,4) =  0.0
      CORN(2,3,4) =  0.0     
!     
!     - face 3
!     
      CORN(3,1,1) =  0.0
      CORN(3,2,1) = -1.0
      CORN(3,3,1) = -1.0
!     
      CORN(3,1,2) =  0.0
      CORN(3,2,2) =  0.0
      CORN(3,3,2) = -1.0
!     
      CORN(3,1,3) =  0.0 
      CORN(3,2,3) = -1.0 
      CORN(3,3,3) =  0.0 
!     
      CORN(3,1,4) =  0.0
      CORN(3,2,4) =  0.0
      CORN(3,3,4) =  0.0
!     
!     - face 4
!     
      CORN(4,1,1) =  0.0
      CORN(4,2,1) =  0.0
      CORN(4,3,1) = -1.0
!     
      CORN(4,1,2) =  0.0
      CORN(4,2,2) =  1.0
      CORN(4,3,2) = -1.0
!     
      CORN(4,1,3) =  0.0
      CORN(4,2,3) =  0.0
      CORN(4,3,3) =  0.0 
!     
      CORN(4,1,4) =  0.0
      CORN(4,2,4) =  1.0
      CORN(4,3,4) =  0.0
!     
!     - face 5
!     
      CORN(5,1,1) = -1.0
      CORN(5,2,1) = -1.0
      CORN(5,3,1) =  0.0
!     
      CORN(5,1,2) =  0.0
      CORN(5,2,2) = -1.0
      CORN(5,3,2) =  0.0
!     
      CORN(5,1,3) = -1.0
      CORN(5,2,3) =  0.0
      CORN(5,3,3) =  0.0
!     
      CORN(5,1,4) =  0.0
      CORN(5,2,4) =  0.0
      CORN(5,3,4) =  0.0
!     
!     - face 6
!     
      CORN(6,1,1) =  0.0
      CORN(6,2,1) = -1.0
      CORN(6,3,1) =  0.0
!     
      CORN(6,1,2) =  1.0
      CORN(6,2,2) = -1.0
      CORN(6,3,2) =  0.0 
!     
      CORN(6,1,3) =  0.0
      CORN(6,2,3) =  0.0
      CORN(6,3,3) =  0.0    
!     
      CORN(6,1,4) =  1.0
      CORN(6,2,4) =  0.0
      CORN(6,3,4) =  0.0   
!     
!     - face 7
!     
      CORN(7,1,1) = -1.0
      CORN(7,2,1) =  0.0
      CORN(7,3,1) =  0.0
!     
      CORN(7,1,2) =  0.0
      CORN(7,2,2) =  0.0
      CORN(7,3,2) =  0.0
!     
      CORN(7,1,3) = -1.0
      CORN(7,2,3) =  1.0
      CORN(7,3,3) =  0.0
!     
      CORN(7,1,4) =  0.0
      CORN(7,2,4) =  1.0
      CORN(7,3,4) =  0.0
!     
!     - face 8
!     
      CORN(8,1,1) = 0.0
      CORN(8,2,1) = 0.0
      CORN(8,3,1) = 0.0
!     
      CORN(8,1,2) = 1.0
      CORN(8,2,2) = 0.0
      CORN(8,3,2) = 0.0
!     
      CORN(8,1,3) = 0.0
      CORN(8,2,3) = 1.0
      CORN(8,3,3) = 0.0
!     
      CORN(8,1,4) = 1.0
      CORN(8,2,4) = 1.0
      CORN(8,3,4) = 0.0
!     
!     - face 9
!     
      CORN(9,1,1) = -1.0
      CORN(9,2,1) =  0.0
      CORN(9,3,1) =  0.0
!     
      CORN(9,1,2) =  0.0
      CORN(9,2,2) =  0.0
      CORN(9,3,2) =  0.0
!     
      CORN(9,1,3) = -1.0
      CORN(9,2,3) =  0.0
      CORN(9,3,3) =  1.0
!     
      CORN(9,1,4) =  0.0
      CORN(9,2,4) =  0.0
      CORN(9,3,4) =  1.0
!     
!     - face 10
!     
      CORN(10,1,1) = 0.0
      CORN(10,2,1) = 0.0
      CORN(10,3,1) = 0.0
!     
      CORN(10,1,2) = 1.0
      CORN(10,2,2) = 0.0
      CORN(10,3,2) = 0.0
!     
      CORN(10,1,3) = 0.0
      CORN(10,2,3) = 0.0
      CORN(10,3,3) = 1.0
!     
      CORN(10,1,4) = 1.0
      CORN(10,2,4) = 0.0
      CORN(10,3,4) = 1.0
!     
!     - face 11
!     
      CORN(11,1,1) =  0.0
      CORN(11,2,1) = -1.0
      CORN(11,3,1) =  0.0
!     
      CORN(11,1,2) =  0.0
      CORN(11,2,2) =  0.0
      CORN(11,3,2) =  0.0
!     
      CORN(11,1,3) =  0.0
      CORN(11,2,3) = -1.0
      CORN(11,3,3) =  1.0
!     
      CORN(11,1,4) =  0.0
      CORN(11,2,4) =  0.0
      CORN(11,3,4) =  1.0
!     
!     - face 12
!     
      CORN(12,1,1) = 0.0 
      CORN(12,2,1) = 0.0
      CORN(12,3,1) = 0.0
!     
      CORN(12,1,2) = 0.0
      CORN(12,2,2) = 1.0
      CORN(12,3,2) = 0.0
!     
      CORN(12,1,3) = 0.0
      CORN(12,2,3) = 0.0
      CORN(12,3,3) = 1.0
!     
      CORN(12,1,4) = 0.0
      CORN(12,2,4) = 1.0
      CORN(12,3,4) = 1.0
!     
!     - define the Gaussian integration points in (XIP,ETAP) space.
!     
      XIGP(1)  = 0.0
      ETAGP(1) = 0.0
!     
!     - define positions of vertices of quadrilateral element 
!     - associated with a subcell face in (XIP,ETAP) space.
!     
      XIP(1)  = -1.0
      ETAP(1) = -1.0
!     
      XIP(2)  =  1.0
      ETAP(2) = -1.0
!     
      XIP(3)  = -1.0
      ETAP(3) =  1.0
!     
      XIP(4)  =  1.0
      ETAP(4) =  1.0
!     
!     - define positions of vertices of hexhedral element in 
!     - (XI,ETA,ZETA) space.
!     
      XI(1)   = -1.0
      ETA(1)  = -1.0
      ZETA(1) = -1.0
!     
      XI(2)   =  1.0
      ETA(2)  = -1.0
      ZETA(2) = -1.0
!     
      XI(3)   = -1.0
      ETA(3)  =  1.0
      ZETA(3) = -1.0
!     
      XI(4)   =  1.0
      ETA(4)  =  1.0
      ZETA(4) = -1.0
!     
      XI(5)   = -1.0
      ETA(5)  = -1.0
      ZETA(5) =  1.0
!     
      XI(6)   =  1.0
      ETA(6)  = -1.0
      ZETA(6) =  1.0
!     
      XI(7)   = -1.0
      ETA(7)  =  1.0
      ZETA(7) =  1.0
!     
      XI(8)   =  1.0
      ETA(8)  =  1.0
      ZETA(8) =  1.0
!     
!     - generate values of the FE basis functions at the quadrature
!     - points on the faces of the subcells.
!     
      do ILOC = 1, NLOC! Was loop 
!     
      do IFACE = 1, NFACE! Was loop 
!     
      do GJ = 1, FNGI! Was loop 
!     
               GI = IFACE
!     
      do ICOORD = 1, NCOORD! Was loop 
!     
                  POS(ICOORD)    = 0.0
                  DPDXI(ICOORD)  = 0.0
                  DPDETA(ICOORD) = 0.0
!     
      do JLOC = 1, FNLOC! Was loop 
!     
                     POS(ICOORD) = POS(ICOORD)&
     &                    + CORN(IFACE,ICOORD,JLOC)&
     &                    *0.25*(1.+XIP(JLOC)*XIGP(GJ))&
     &                    *(1.+ETAP(JLOC)*ETAGP(GJ))
!     
                     DPDXI(ICOORD) = DPDXI(ICOORD)&
     &                    + CORN(IFACE,ICOORD,JLOC)&
     &                    *0.25*XIP(JLOC)&
     &                    *(1.+ETAP(JLOC)*ETAGP(GJ))
!     
                     DPDETA(ICOORD) = DPDETA(ICOORD)&
     &                    + CORN(IFACE,ICOORD,JLOC)&
     &                    *0.25*(1.+XIP(JLOC)*XIGP(GJ))&
     &                    *ETAP(JLOC)
!     
                  END DO
!     
               END DO
!     
               SVN(ILOC,GI) = 0.125*( 1.0 + XI(ILOC)*POS(1) )&
     &              *( 1.0 + ETA(ILOC)*POS(2) )&
     &              *( 1.0 + ZETA(ILOC)*POS(3) )
!     
               SVNLX(ILOC,GI) = 0.125*XI(ILOC)*DPDXI(1)&
     &              *( 1.0 + ETA(ILOC)*POS(2) )&
     &              *( 1.0 + ZETA(ILOC)*POS(3) )&
     &              + 0.125*( 1.0 + XI(ILOC)*POS(1) )&
     &              *ETA(ILOC)*DPDXI(2)&
     &              *( 1.0 + ZETA(ILOC)*POS(3) )&
     &              + 0.125*( 1.0 + XI(ILOC)*POS(1) )&
     &              *( 1.0 + ETA(ILOC)*POS(2) )&
     &              *ZETA(ILOC)*DPDXI(3)
!     
               SVNLY(ILOC,GI) = 0.125*XI(ILOC)*DPDETA(1)&
     &              *( 1.0 + ETA(ILOC)*POS(2) )&
     &              *( 1.0 + ZETA(ILOC)*POS(3) )&
     &              + 0.125*( 1.0 + XI(ILOC)*POS(1) )&
     &              *ETA(ILOC)*DPDETA(2)&
     &              *( 1.0 + ZETA(ILOC)*POS(3) )&
     &              + 0.125*( 1.0 + XI(ILOC)*POS(1) )&
     &              *( 1.0 + ETA(ILOC)*POS(2) )&
     &              *ZETA(ILOC)*DPDETA(3)
!     
            END DO
!     
         END DO
!     
      END DO 
!     
!     - set the weights for the surface integration note that this
!     - is 2*2.0 = 4.0 for the weight because it is a product of
!     - two 1-D weights for integration over a 2-D quadrilateral.
!     
      do GI = 1, SVNGI! Was loop 
!     
         SVWEIGH(GI) = 4.0
!     
      END DO
!     
      END SUBROUTINE
!     
      SUBROUTINE FVTET( NGI,   NLOC, SVNGI,&
     &     M,     SVN,  SVNLX,&
     &     SVNLY, SVWEIGH      )
!     ---------------------------------------
!     
!     - this routine generates the shape functions associated
!     - with the FV's i.e. their surfaces and volume shape 
!     - functions and derivatives. The surface shape functions
!     - are the values of the FE volume shape functions evaluated
!     - on the surfaces of the CV's.
!     
!     -------------------------------
!     - date last modified : 06/10/2002
!     -------------------------------
!     
      use fldebug
      IMPLICIT NONE
!     
      INTEGER NGI, NLOC, SVNGI
!     
      REAL M(NLOC,NGI),     SVWEIGH(SVNGI)
!     
      REAL SVN(NLOC,SVNGI), SVNLX(NLOC,SVNGI), SVNLY(NLOC,SVNGI)
!     
!     - local variables
!     
      INTEGER ILOC, JLOC, IFACE, ICOORD
!     
      INTEGER GI, GJ
!     
!     - note that NCOORD is the number of co-ordinates.
!     - NFACE is the number of surfaces internal to an 
!     - element 
!     
      INTEGER NCOORD, NFACE
      PARAMETER( NCOORD = 4, NFACE = 6 )
!     
      INTEGER FNGI, FNLOC
      PARAMETER( FNGI = 1, FNLOC = 4 )
!     
      REAL POS(NCOORD),    DPDXI(NCOORD)
      REAL DPDETA(NCOORD)
!     
      REAL XIGP(FNGI), ETAGP(FNGI)
!     
      REAL XI(FNLOC), ETA(FNLOC)
!     
      REAL CORN(NFACE,NCOORD,FNLOC)
!     
      INTEGER COUNT
!
!     - FV basis functions for use in calculating the integral of the
!     - shape functions over a subcell of a CV.
!     
      ewrite(2,*) 'SUBROUTINE FVTET()'
      
      M=0.0
      COUNT = 1
      do ILOC = 1, NLOC! Was loop 
      do GI = COUNT, COUNT+7! Was loop 
          M(ILOC,GI) = 1.0
        END DO
        COUNT = COUNT + 8
      END DO

!     
!     - Note that CORN(IFACE,ICOORD,JLOC) tabulates the positions
!     - of the corners (vertices) of the faces of a subcell in volume 
!     - co-ordinates. Note also that IFACE ranges from 1-6, ICOORD
!     - ranges from 1-4 and JLOC ranges from 1-4. IFACE signifies 
!     - the particular subcell face, ICOORD signfies the co-ordinates
!     - of the vertices of the face in volume co-ordinates (L1,L2,L3,L4)
!     - and JLOC signifies the vertex of the face of the subcell
!     - which are quadrilaterals.
!     
!     - face 1
!     
      CORN(1,1,1) = 0.5
      CORN(1,2,1) = 0.5 
      CORN(1,3,1) = 0.0
      CORN(1,4,1) = 0.0
!     
      CORN(1,1,2) = 0.33333333
      CORN(1,2,2) = 0.33333333
      CORN(1,3,2) = 0.33333333
      CORN(1,4,2) = 0.0
!     
      CORN(1,1,3) = 0.33333333
      CORN(1,2,3) = 0.33333333
      CORN(1,3,3) = 0.0
      CORN(1,4,3) = 0.33333333
!     
      CORN(1,1,4) = 0.25
      CORN(1,2,4) = 0.25
      CORN(1,3,4) = 0.25
      CORN(1,4,4) = 0.25
!     
!     - face 2
!     
      CORN(2,1,1) = 0.33333333
      CORN(2,2,1) = 0.33333333
      CORN(2,3,1) = 0.33333333
      CORN(2,4,1) = 0.0
!     
      CORN(2,1,2) = 0.0
      CORN(2,2,2) = 0.5 
      CORN(2,3,2) = 0.5
      CORN(2,4,2) = 0.0
!     
      CORN(2,1,3) = 0.25
      CORN(2,2,3) = 0.25 
      CORN(2,3,3) = 0.25
      CORN(2,4,3) = 0.25
!     
      CORN(2,1,4) = 0.0
      CORN(2,2,4) = 0.33333333
      CORN(2,3,4) = 0.33333333
      CORN(2,4,4) = 0.33333333
!     
!     - face 3
!     
      CORN(3,1,1) = 0.33333333
      CORN(3,2,1) = 0.33333333
      CORN(3,3,1) = 0.33333333
      CORN(3,4,1) = 0.0
!     
      CORN(3,1,2) = 0.5
      CORN(3,2,2) = 0.0 
      CORN(3,3,2) = 0.5
      CORN(3,4,2) = 0.0
!     
      CORN(3,1,3) = 0.25
      CORN(3,2,3) = 0.25 
      CORN(3,3,3) = 0.25
      CORN(3,4,3) = 0.25
!     
      CORN(3,1,4) = 0.33333333
      CORN(3,2,4) = 0.0
      CORN(3,3,4) = 0.33333333
      CORN(3,4,4) = 0.33333333
!     
!     - face 4
!     
      CORN(4,1,1) = 0.5
      CORN(4,2,1) = 0.0
      CORN(4,3,1) = 0.0
      CORN(4,4,1) = 0.5
!     
      CORN(4,1,2) = 0.33333333
      CORN(4,2,2) = 0.33333333
      CORN(4,3,2) = 0.0
      CORN(4,4,2) = 0.33333333
!     
      CORN(4,1,3) = 0.33333333
      CORN(4,2,3) = 0.0
      CORN(4,3,3) = 0.33333333
      CORN(4,4,3) = 0.33333333
!     
      CORN(4,1,4) = 0.25
      CORN(4,2,4) = 0.25
      CORN(4,3,4) = 0.25
      CORN(4,4,4) = 0.25
!     
!     - face 5
!     
      CORN(5,1,1) = 0.33333333
      CORN(5,2,1) = 0.33333333
      CORN(5,3,1) = 0.0
      CORN(5,4,1) = 0.33333333
!     
      CORN(5,1,2) = 0.0
      CORN(5,2,2) = 0.5
      CORN(5,3,2) = 0.0
      CORN(5,4,2) = 0.5
!     
      CORN(5,1,3) = 0.25
      CORN(5,2,3) = 0.25
      CORN(5,3,3) = 0.25
      CORN(5,4,3) = 0.25
!     
      CORN(5,1,4) = 0.0
      CORN(5,2,4) = 0.33333333
      CORN(5,3,4) = 0.33333333
      CORN(5,4,4) = 0.33333333
!     
!     - face 6
!     
      CORN(6,1,1) = 0.25
      CORN(6,2,1) = 0.25
      CORN(6,3,1) = 0.25
      CORN(6,4,1) = 0.25
!     
      CORN(6,1,2) = 0.0
      CORN(6,2,2) = 0.33333333
      CORN(6,3,2) = 0.33333333
      CORN(6,4,2) = 0.33333333
!     
      CORN(6,1,3) = 0.33333333
      CORN(6,2,3) = 0.0
      CORN(6,3,3) = 0.33333333
      CORN(6,4,3) = 0.33333333
!     
      CORN(6,1,4) = 0.0
      CORN(6,2,4) = 0.0
      CORN(6,3,4) = 0.5
      CORN(6,4,4) = 0.5
!     
!     - define the Gaussian integration points in (XI,ETA) space.
!     
      XIGP(1)  = 0.0
      ETAGP(1) = 0.0
!     
!     - define positions of vertices of quadrilateral element 
!     - associated with a subcell face in (XI,ETA) space.
!     
      XI(1)  = -1.0
      ETA(1) = -1.0
!     
      XI(2)  =  1.0
      ETA(2) = -1.0
!     
      XI(3)  = -1.0
      ETA(3) =  1.0
!     
      XI(4)  =  1.0
      ETA(4) =  1.0
!     
!     - generate values of the FE basis functions at the quadrature
!     - points on the faces of the subcells.
!     
      do ILOC = 1, NLOC! Was loop 
!     
      do IFACE = 1, NFACE! Was loop 
!     
      do GJ = 1, FNGI! Was loop 
!     
               GI = IFACE
!     
      do ICOORD = 1, NCOORD! Was loop 
!     
                  POS(ICOORD)    = 0.0
                  DPDXI(ICOORD)  = 0.0
                  DPDETA(ICOORD) = 0.0
!     
      do JLOC = 1, FNLOC! Was loop 
!     
                     POS(ICOORD) = POS(ICOORD)&
     &                    + CORN(IFACE,ICOORD,JLOC)&
     &                    *0.25*(1.+XI(JLOC)*XIGP(GJ))&
     &                    *(1.+ETA(JLOC)*ETAGP(GJ))
!     
                     DPDXI(ICOORD) = DPDXI(ICOORD)&
     &                    + CORN(IFACE,ICOORD,JLOC)&
     &                    *0.25*XI(JLOC)&
     &                    *(1.+ETA(JLOC)*ETAGP(GJ))
!     
                     DPDETA(ICOORD) = DPDETA(ICOORD)&
     &                    + CORN(IFACE,ICOORD,JLOC)&
     &                    *0.25*(1.+XI(JLOC)*XIGP(GJ))&
     &                    *ETA(JLOC)
!     
                  END DO
!     
               END DO
!     
               SVN(ILOC,GI)   = POS(ILOC)
               SVNLX(ILOC,GI) = DPDXI(ILOC)
               SVNLY(ILOC,GI) = DPDETA(ILOC) 
!     
            END DO
!     
         END DO
!     
      END DO         
!     
!     - set the weights for the quadrature points
!     
      do GI = 1, SVNGI! Was loop 
!     
         SVWEIGH(GI) = 4.0
!     
      END DO
!     
      END SUBROUTINE
!     
!     
      SUBROUTINE FVQQUAD( NGI,    NLOC,    SVNGI,&
!     - REALS
     &     M,      SVN,     SVNLX, &
     &     SVWEIGH                 )
!     ---------------------------------------------
!     
!     - this routine generates the shape functions associated
!     - with the FV's i.e. their surfaces and volume shape 
!     - functions and derivatives. The surface shape functions
!     - are the values of the FE volume shape functions evaluated
!     - on the surfaces of the CV's.
!     
      IMPLICIT NONE
!     
      INTEGER NGI, NLOC, SVNGI
!     
      REAL M(NLOC,NGI),     SVWEIGH(SVNGI)
!     
      REAL SVN(NLOC,SVNGI), SVNLX(NLOC,SVNGI)
!     
!     - local variables
!     
      INTEGER ILOC, JLOC, IFACE, ICOORD
!     
      INTEGER COUNT, GI, GJ
!     
!     - note that NCOORD is the number of co-ordinates.
!     - NFACE is the number of surfaces internal to an 
!     - element 
!     
      INTEGER NCOORD, NFACE
      PARAMETER( NCOORD = 2, NFACE = 12 )
!     
      INTEGER FNGI, FNLOC
      PARAMETER( FNGI = 2, FNLOC = 2  )
!     
      REAL POS(NCOORD),    DPDXI(NCOORD)
      REAL DPDETA(NCOORD)
!     
      REAL XI, ETA
!     
      REAL XIGP(FNGI)
!     
      REAL XIP(FNLOC)
!     
      REAL CORN(NFACE,NCOORD,FNLOC)
!     
      REAL LIJXI(3),   LIJETA(3)
      REAL DLIJXIDXI(3), DLIJETADXI(3)    
!     
      INTEGER I, J
!     
!     - FV basis functions for use in calculating the integral of the
!     - shape functions over a subcell of a CV.
!     
      M = 0.0
!     
      COUNT = 1
!     
      do ILOC = 1, NLOC! Was loop 
!     
      do GI = COUNT, COUNT + 3! Was loop 
!     
            M(ILOC,GI) = 1.0
!     
         END DO
!     
         COUNT = COUNT + 4
!     
      END DO
!     
!     - Note that CORN(IFACE,ICOORD,JLOC) tabulates the positions
!     - of the corners (vertices) of the faces of a subcell in volume 
!     - co-ordinates. Note also that IFACE ranges from 1-4, ICOORD
!     - ranges from 1-2 and JLOC ranges from 1-2. IFACE signifies 
!     - the particular subcell face, ICOORD signfies the co-ordinates
!     - of the vertices of the face in (XI,ETA) co-ordinates
!     - and JLOC signifies the vertex of the face of the subcell
!     - which are straight lines.
!     
!     - face 1
!     
      CORN(1,1,1) = -0.5
      CORN(1,2,1) = -1.0
!     
      CORN(1,1,2) = -0.5
      CORN(1,2,2) = -0.5
!     
!     - face 2
!     
      CORN(2,1,1) =  0.5
      CORN(2,2,1) = -1.0
!     
      CORN(2,1,2) =  0.5
      CORN(2,2,2) = -0.5
!     
!     - face 3
!     
      CORN(3,1,1) = -1.0
      CORN(3,2,1) = -0.5
!     
      CORN(3,1,2) = -0.5
      CORN(3,2,2) = -0.5
!     
!     - face 4
!     
      CORN(4,1,1) = -0.5
      CORN(4,2,1) = -0.5
!     
      CORN(4,1,2) =  0.5
      CORN(4,2,2) = -0.5
!     
!     - face 5
!     
      CORN(5,1,1) =  0.5
      CORN(5,2,1) = -0.5
!     
      CORN(5,1,2) =  1.0
      CORN(5,2,2) = -0.5
!     
!     - face 6
!     
      CORN(6,1,1) = -0.5
      CORN(6,2,1) = -0.5
!     
      CORN(6,1,2) = -0.5
      CORN(6,2,2) =  0.5
!     
!     - face 7
!     
      CORN(7,1,1) =  0.5
      CORN(7,2,1) = -0.5
!     
      CORN(7,1,2) =  0.5
      CORN(7,2,2) =  0.5
!     
!     - face 8
!     
      CORN(8,1,1) = -1.0
      CORN(8,2,1) =  0.5
!     
      CORN(8,1,2) = -0.5
      CORN(8,2,2) =  0.5
!     
!     - face 9
!     
      CORN(9,1,1) = -0.5
      CORN(9,2,1) =  0.5
!     
      CORN(9,1,2) =  0.5
      CORN(9,2,2) =  0.5
!     
!     - face 10
!     
      CORN(10,1,1) =  0.5
      CORN(10,2,1) =  0.5
!     
      CORN(10,1,2) =  1.0
      CORN(10,2,2) =  0.5
!     
!     - face 11
!     
      CORN(11,1,1) = -0.5
      CORN(11,2,1) =  0.5
!     
      CORN(11,1,2) = -0.5
      CORN(11,2,2) =  1.0
!     
!     - face 12
!     
      CORN(12,1,1) =  0.5
      CORN(12,2,1) =  0.5
!     
      CORN(12,1,2) =  0.5
      CORN(12,2,2) =  1.0
!     
!     - define the Gaussian integration points in XIP space.
!     
      XIGP(1)  = -1/SQRT(3.0) 
      XIGP(2)  =  1/SQRT(3.0) 
!     
!     - define positions of vertices of line element
!     - associated with a subcell face in XIP
!     
      XIP(1) = -1.0
      XIP(2) =  1.0
!     
!     - generate values of the FE basis functions at the quadrature
!     - points on the faces of the subcells.
!     
      do IFACE = 1, NFACE! Was loop 
!     
      do GJ = 1, FNGI! Was loop 
!     
            GI = (IFACE-1)*FNGI+GJ
!     
      do ICOORD = 1, NCOORD! Was loop 
!     
               POS(ICOORD)    = 0.0
               DPDXI(ICOORD)  = 0.0
               DPDETA(ICOORD) = 0.0
!     
      do JLOC = 1, FNLOC! Was loop 
!     
                  POS(ICOORD) = POS(ICOORD)&
     &                 + CORN(IFACE,ICOORD,JLOC)&
     &                 *0.5*(1.+XIP(JLOC)*XIGP(GJ))
!     
                  DPDXI(ICOORD) = DPDXI(ICOORD)&
     &                 + CORN(IFACE,ICOORD,JLOC)&
     &                 *0.5*XIP(JLOC)
!     
               END DO
!     
            END DO
!     
            XI  = POS(1)
            ETA = POS(2)
!     
            LIJXI(1)   =  0.5*XI*(XI-1.0)
            LIJXI(2)   =  1.0-XI*XI
            LIJXI(3)   =  0.5*XI*(XI+1.0)
!     
            LIJETA(1)  =  0.5*ETA*(ETA-1.0)
            LIJETA(2)  =  1.0-ETA*ETA
            LIJETA(3)  =  0.5*ETA*(ETA+1.0)
!     
            DLIJXIDXI(1) =  0.5*(2.0*XI-1.0)*DPDXI(1)
            DLIJXIDXI(2) = -2.0*XI*DPDXI(1)
            DLIJXIDXI(3) =  0.5*(2.0*XI+1.0)*DPDXI(1)
!     
            DLIJETADXI(1) =  0.5*(2.0*ETA-1.0)*DPDXI(2)
            DLIJETADXI(2) = -2.0*ETA*DPDXI(2)
            DLIJETADXI(3) =  0.5*(2.0*ETA+1.0)*DPDXI(2)
!     
      do I = 1, 3! Was loop 
!     
      do J = 1, 3! Was loop 
!     
                  ILOC = I+(J-1)*3 
!     
                  SVN(ILOC,GI) = LIJXI(I)*LIJETA(J)
!     
                  SVNLX(ILOC,GI) = LIJXI(I)*DLIJETADXI(J)&
     &                 + LIJETA(J)*DLIJXIDXI(I)
!     
               END DO
!     
            END DO
!     
         END DO
!     
      END DO
!     
!     - set the weights for the surface integration
!     
      do GI = 1, SVNGI! Was loop 
!     
         SVWEIGH(GI) = 1.0
!     
      END DO
!     
      END SUBROUTINE

      SUBROUTINE FVQHEX( NGI,   NLOC,   SVNGI,&
!     - REALS            
     &     M,     SVN,    SVNLX,&
     &     SVNLY, SVWEIGH        )
!     ------------------------------------------
!     
!     - this routine generates the shape functions associated
!     - with the FV's i.e. their surfaces and volume shape 
!     - functions and derivatives. The surface shape functions
!     - are the values of the FE volume shape functions evaluated
!     - on the surfaces of the CV's.
!     
!     -------------------------------
!     - date last modified : 07/11/2002
!     -------------------------------
!     
      IMPLICIT NONE
!     
      INTEGER NGI, NLOC, SVNGI
!     
      REAL M(NLOC,NGI),     SVWEIGH(SVNGI)
!     
      REAL SVN(NLOC,SVNGI), SVNLX(NLOC,SVNGI)
!     
      REAL SVNLY(NLOC,SVNGI)
!     
!     - local variables
!     
      INTEGER ILOC, JLOC, IFACE, ICOORD
!     
      INTEGER COUNT, GI, GJ
!     
!     - note that NCOORD is the number of co-ordinates.
!     - NFACE is the number of surfaces internal to an 
!     - element 
!     
      INTEGER NCOORD, NFACE
      PARAMETER( NCOORD = 3, NFACE = 54 )
!     
      INTEGER FNGI, FNLOC
      PARAMETER( FNGI = 4, FNLOC = 4  )
!     
      REAL POS(NCOORD),    DPDXI(NCOORD)
      REAL DPDETA(NCOORD)
!     
      REAL XI, ETA, ZETA
!     
      REAL XIPGP(FNGI), ETAPGP(FNGI)
!     
      REAL XIP(FNLOC), ETAP(FNLOC)
!     
      REAL CORN(NFACE,NCOORD,FNLOC)
!     
      REAL LIJXI(3),      LIJETA(3),      LIJZETA(3)
      REAL DLIJXIDXI(3),  DLIJETADXI(3),  DLIJZETADXI(3)
      REAL DLIJXIDETA(3), DLIJETADETA(3), DLIJZETADETA(3)    
!     
      INTEGER I, J, K 
!     
!     - FV basis functions for use in calculating the integral of the
!     - shape functions over a subcell of a CV.
!     
      M = 0.0
!     
      COUNT = 1
!     
!     DO ILOC = 1, NLOC
!     
!     DO GI = COUNT, COUNT + 7
!     
!     M(ILOC,GI) = 1.0
!     
!     END DO
!     
!     COUNT = COUNT + 8
!     
!     END DO     
!     
!     - Note that CORN(IFACE,ICOORD,JLOC) tabulates the positions
!     - of the corners (vertices) of the faces of a subcell in volume 
!     - co-ordinates. Note also that IFACE ranges from 1-54, ICOORD
!     - ranges from 1-3 and JLOC ranges from 1-4. IFACE signifies 
!     - the particular subcell face, ICOORD signfies the co-ordinates
!     - of the vertices of the face in (XI,ETA) co-ordinates
!     - and JLOC signifies the vertex of the face of the subcell
!     - which are straight lines.
!     
!     - face 1
!     
      CORN(1,1,1) = -0.5
      CORN(1,2,1) = -1.0
      CORN(1,3,1) = -1.0
!     
      CORN(1,1,2) = -0.5
      CORN(1,2,2) = -0.5
      CORN(1,3,2) = -1.0
!     
      CORN(1,1,3) = -0.5
      CORN(1,2,3) = -1.0
      CORN(1,3,3) = -0.5
!     
      CORN(1,1,4) = -0.5
      CORN(1,2,4) = -0.5
      CORN(1,3,4) = -0.5
!     
!     - face 2
!     
      CORN(2,1,1) =  0.5
      CORN(2,2,1) = -1.0
      CORN(2,3,1) = -1.0
!     
      CORN(2,1,2) =  0.5
      CORN(2,2,2) = -0.5
      CORN(2,3,2) = -1.0
!     
      CORN(2,1,3) =  0.5
      CORN(2,2,3) = -1.0
      CORN(2,3,3) = -0.5
!     
      CORN(2,1,4) =  0.5
      CORN(2,2,4) = -0.5
      CORN(2,3,4) = -0.5
!     
!     - face 3
!     
      CORN(3,1,1) = -1.0
      CORN(3,2,1) = -0.5
      CORN(3,3,1) = -1.0
!     
      CORN(3,1,2) = -0.5
      CORN(3,2,2) = -0.5
      CORN(3,3,2) = -1.0
!     
      CORN(3,1,3) = -1.0
      CORN(3,2,3) = -0.5
      CORN(3,3,3) = -0.5
!     
      CORN(3,1,4) = -0.5
      CORN(3,2,4) = -0.5
      CORN(3,3,4) = -0.5
!     
!     - face 4
!     
      CORN(4,1,1) = -0.5
      CORN(4,2,1) = -0.5
      CORN(4,3,1) = -1.0
!     
      CORN(4,1,2) =  0.5
      CORN(4,2,2) = -0.5
      CORN(4,3,2) = -1.0
!     
      CORN(4,1,3) = -0.5
      CORN(4,2,3) = -0.5
      CORN(4,3,3) = -0.5
!     
      CORN(4,1,4) =  0.5
      CORN(4,2,4) = -0.5
      CORN(4,3,4) = -0.5
!     
!     - face 5
!     
      CORN(5,1,1) =  0.5
      CORN(5,2,1) = -0.5
      CORN(5,3,1) = -1.0
!     
      CORN(5,1,2) =  1.0
      CORN(5,2,2) = -0.5
      CORN(5,3,2) = -1.0
!     
      CORN(5,1,3) =  0.5
      CORN(5,2,3) = -0.5
      CORN(5,3,3) = -0.5
!     
      CORN(5,1,4) =  1.0
      CORN(5,2,4) = -0.5
      CORN(5,3,4) = -0.5
!     
!     - face 6
!     
      CORN(6,1,1) = -0.5
      CORN(6,2,1) = -0.5
      CORN(6,3,1) = -1.0
!     
      CORN(6,1,2) = -0.5
      CORN(6,2,2) =  0.5
      CORN(6,3,2) = -1.0
!     
      CORN(6,1,3) = -0.5
      CORN(6,2,3) = -0.5
      CORN(6,3,3) = -0.5
!     
      CORN(6,1,4) = -0.5
      CORN(6,2,4) =  0.5
      CORN(6,3,4) = -0.5
!     
!     - face 7
!     
      CORN(7,1,1) =  0.5
      CORN(7,2,1) = -0.5
      CORN(7,3,1) = -1.0
!     
      CORN(7,1,2) =  0.5
      CORN(7,2,2) =  0.5
      CORN(7,3,2) = -1.0
!     
      CORN(7,1,3) =  0.5
      CORN(7,2,3) = -0.5
      CORN(7,3,3) = -0.5
!     
      CORN(7,1,4) =  0.5
      CORN(7,2,4) =  0.5
      CORN(7,3,4) = -0.5
!     
!     - face 8
!     
      CORN(8,1,1) = -1.0
      CORN(8,2,1) =  0.5
      CORN(8,3,1) = -1.0
!     
      CORN(8,1,2) = -0.5
      CORN(8,2,2) =  0.5
      CORN(8,3,2) = -1.0
!     
      CORN(8,1,3) = -1.0
      CORN(8,2,3) =  0.5
      CORN(8,3,3) = -0.5
!     
      CORN(8,1,4) = -0.5
      CORN(8,2,4) =  0.5
      CORN(8,3,4) = -0.5
!     
!     - face 9
!     
      CORN(9,1,1) = -0.5
      CORN(9,2,1) =  0.5
      CORN(9,3,1) = -1.0
!     
      CORN(9,1,2) =  0.5
      CORN(9,2,2) =  0.5
      CORN(9,3,2) = -1.0
!     
      CORN(9,1,3) = -0.5
      CORN(9,2,3) =  0.5
      CORN(9,3,3) = -0.5
!     
      CORN(9,1,4) =  0.5
      CORN(9,2,4) =  0.5
      CORN(9,3,4) = -0.5
!     
!     - face 10
!     
      CORN(10,1,1) =  0.5
      CORN(10,2,1) =  0.5
      CORN(10,3,1) = -1.0
!     
      CORN(10,1,2) =  1.0
      CORN(10,2,2) =  0.5
      CORN(10,3,2) = -1.0
!     
      CORN(10,1,3) =  0.5
      CORN(10,2,3) =  0.5
      CORN(10,3,3) = -0.5
!     
      CORN(10,1,4) =  1.0
      CORN(10,2,4) =  0.5
      CORN(10,3,4) = -0.5
!     
!     - face 11
!     
      CORN(11,1,1) = -0.5
      CORN(11,2,1) =  0.5
      CORN(11,3,1) = -1.0
!     
      CORN(11,1,2) = -0.5
      CORN(11,2,2) =  1.0
      CORN(11,3,2) = -1.0
!     
      CORN(11,1,3) = -0.5
      CORN(11,2,3) =  0.5
      CORN(11,3,3) = -0.5
!     
      CORN(11,1,4) = -0.5
      CORN(11,2,4) =  1.0
      CORN(11,3,4) = -0.5
!     
!     - face 12
!     
      CORN(12,1,1) =  0.5
      CORN(12,2,1) =  0.5
      CORN(12,3,1) = -1.0
!     
      CORN(12,1,2) =  0.5
      CORN(12,2,2) =  1.0
      CORN(12,3,2) = -1.0
!     
      CORN(12,1,3) =  0.5
      CORN(12,2,3) =  0.5
      CORN(12,3,3) = -0.5
!     
      CORN(12,1,4) =  0.5
      CORN(12,2,4) =  1.0
      CORN(12,3,4) = -0.5 
!     
!     - face 13
!     
      CORN(13,1,1) = -1.0
      CORN(13,2,1) = -1.0
      CORN(13,3,1) = -0.5
!     
      CORN(13,1,2) = -0.5
      CORN(13,2,2) = -1.0
      CORN(13,3,2) = -0.5
!     
      CORN(13,1,3) = -1.0
      CORN(13,2,3) = -0.5
      CORN(13,3,3) = -0.5
!     
      CORN(13,1,4) = -0.5
      CORN(13,2,4) = -0.5
      CORN(13,3,4) = -0.5
!     
!     - face 14
!     
      CORN(14,1,1) = -0.5
      CORN(14,2,1) = -1.0
      CORN(14,3,1) = -0.5
!     
      CORN(14,1,2) =  0.5
      CORN(14,2,2) = -1.0
      CORN(14,3,2) = -0.5
!     
      CORN(14,1,3) = -0.5
      CORN(14,2,3) = -0.5
      CORN(14,3,3) = -0.5
!     
      CORN(14,1,4) =  0.5
      CORN(14,2,4) = -0.5
      CORN(14,3,4) = -0.5
!     
!     - face 15
!     
      CORN(15,1,1) =  0.5
      CORN(15,2,1) = -1.0
      CORN(15,3,1) = -0.5
!     
      CORN(15,1,2) =  1.0
      CORN(15,2,2) = -1.0
      CORN(15,3,2) = -0.5
!     
      CORN(15,1,3) =  0.5
      CORN(15,2,3) = -0.5
      CORN(15,3,3) = -0.5
!     
      CORN(15,1,4) =  1.0
      CORN(15,2,4) = -0.5
      CORN(15,3,4) = -0.5
!     
!     - face 16
!     
      CORN(16,1,1) = -1.0
      CORN(16,2,1) = -0.5
      CORN(16,3,1) = -0.5
!     
      CORN(16,1,2) = -0.5
      CORN(16,2,2) = -0.5
      CORN(16,3,2) = -0.5
!     
      CORN(16,1,3) = -1.0
      CORN(16,2,3) =  0.5
      CORN(16,3,3) = -0.5
!     
      CORN(16,1,4) = -0.5
      CORN(16,2,4) =  0.5
      CORN(16,3,4) = -0.5
!     
!     - face 17
!     
      CORN(17,1,1) = -0.5
      CORN(17,2,1) = -0.5
      CORN(17,3,1) = -0.5
!     
      CORN(17,1,2) =  0.5
      CORN(17,2,2) = -0.5
      CORN(17,3,2) = -0.5
!     
      CORN(17,1,3) = -0.5
      CORN(17,2,3) =  0.5
      CORN(17,3,3) = -0.5
!     
      CORN(17,1,4) =  0.5
      CORN(17,2,4) =  0.5
      CORN(17,3,4) = -0.5
!     
!     - face 18
!     
      CORN(18,1,1) =  0.5
      CORN(18,2,1) = -0.5
      CORN(18,3,1) = -0.5
!     
      CORN(18,1,2) =  1.0
      CORN(18,2,2) = -0.5
      CORN(18,3,2) = -0.5
!     
      CORN(18,1,3) =  0.5
      CORN(18,2,3) =  0.5
      CORN(18,3,3) = -0.5
!     
      CORN(18,1,4) =  1.0
      CORN(18,2,4) =  0.5
      CORN(18,3,4) = -0.5
!     
!     - face 19
!     
      CORN(19,1,1) = -1.0
      CORN(19,2,1) =  0.5
      CORN(19,3,1) = -0.5
!     
      CORN(19,1,2) = -0.5
      CORN(19,2,2) =  0.5
      CORN(19,3,2) = -0.5
!     
      CORN(19,1,3) = -1.0
      CORN(19,2,3) =  1.0
      CORN(19,3,3) = -0.5
!     
      CORN(19,1,4) = -0.5
      CORN(19,2,4) =  1.0
      CORN(19,3,4) = -0.5
!     
!     - face 20
!     
      CORN(20,1,1) = -0.5
      CORN(20,2,1) =  0.5
      CORN(20,3,1) = -0.5
!     
      CORN(20,1,2) =  0.5
      CORN(20,2,2) =  0.5
      CORN(20,3,2) = -0.5
!     
      CORN(20,1,3) = -0.5
      CORN(20,2,3) =  1.0
      CORN(20,3,3) = -0.5
!     
      CORN(20,1,4) =  0.5
      CORN(20,2,4) =  1.0
      CORN(20,3,4) = -0.5
!     
!     - face 21
!     
      CORN(21,1,1) =  0.5
      CORN(21,2,1) =  0.5
      CORN(21,3,1) = -0.5
!     
      CORN(21,1,2) =  1.0
      CORN(21,2,2) =  0.5
      CORN(21,3,2) = -0.5
!     
      CORN(21,1,3) =  0.5
      CORN(21,2,3) =  1.0
      CORN(21,3,3) = -0.5
!     
      CORN(21,1,4) =  1.0
      CORN(21,2,4) =  1.0
      CORN(21,3,4) = -0.5
!     
!     - face 22
!     
      CORN(22,1,1) = -0.5
      CORN(22,2,1) = -1.0
      CORN(22,3,1) = -0.5
!     
      CORN(22,1,2) = -0.5
      CORN(22,2,2) = -0.5
      CORN(22,3,2) = -0.5
!     
      CORN(22,1,3) = -0.5
      CORN(22,2,3) = -1.0
      CORN(22,3,3) =  0.5
!     
      CORN(22,1,4) = -0.5
      CORN(22,2,4) = -0.5
      CORN(22,3,4) =  0.5
!     
!     - face 23
!     
      CORN(23,1,1) =  0.5
      CORN(23,2,1) = -1.0
      CORN(23,3,1) = -0.5
!     
      CORN(23,1,2) =  0.5
      CORN(23,2,2) = -0.5
      CORN(23,3,2) = -0.5
!     
      CORN(23,1,3) =  0.5
      CORN(23,2,3) = -1.0
      CORN(23,3,3) =  0.5
!     
      CORN(23,1,4) =  0.5
      CORN(23,2,4) = -0.5
      CORN(23,3,4) =  0.5
!     
!     - face 24
!     
      CORN(24,1,1) = -1.0
      CORN(24,2,1) = -0.5
      CORN(24,3,1) = -0.5
!     
      CORN(24,1,2) = -0.5
      CORN(24,2,2) = -0.5
      CORN(24,3,2) = -0.5
!     
      CORN(24,1,3) = -1.0
      CORN(24,2,3) = -0.5
      CORN(24,3,3) =  0.5
!     
      CORN(24,1,4) = -0.5
      CORN(24,2,4) = -0.5
      CORN(24,3,4) =  0.5
!     
!     - face 25
!     
      CORN(25,1,1) = -0.5
      CORN(25,2,1) = -0.5
      CORN(25,3,1) = -0.5
!     
      CORN(25,1,2) =  0.5
      CORN(25,2,2) = -0.5
      CORN(25,3,2) = -0.5
!     
      CORN(25,1,3) = -0.5
      CORN(25,2,3) = -0.5
      CORN(25,3,3) =  0.5
!     
      CORN(25,1,4) =  0.5
      CORN(25,2,4) = -0.5
      CORN(25,3,4) =  0.5
!     
!     - face 26
!     
      CORN(26,1,1) =  0.5
      CORN(26,2,1) = -0.5
      CORN(26,3,1) = -0.5
!     
      CORN(26,1,2) =  1.0
      CORN(26,2,2) = -0.5
      CORN(26,3,2) = -0.5
!     
      CORN(26,1,3) =  0.5
      CORN(26,2,3) = -0.5
      CORN(26,3,3) =  0.5
!     
      CORN(26,1,4) =  1.0
      CORN(26,2,4) = -0.5
      CORN(26,3,4) =  0.5
!     
!     - face 27
!     
      CORN(27,1,1) = -0.5
      CORN(27,2,1) = -0.5
      CORN(27,3,1) = -0.5
!     
      CORN(27,1,2) = -0.5
      CORN(27,2,2) =  0.5
      CORN(27,3,2) = -0.5
!     
      CORN(27,1,3) = -0.5
      CORN(27,2,3) = -0.5
      CORN(27,3,3) =  0.5
!     
      CORN(27,1,4) = -0.5
      CORN(27,2,4) =  0.5
      CORN(27,3,4) =  0.5
!     
!     - face 28
!     
      CORN(28,1,1) =  0.5
      CORN(28,2,1) = -0.5
      CORN(28,3,1) = -0.5
!     
      CORN(28,1,2) =  0.5
      CORN(28,2,2) =  0.5
      CORN(28,3,2) = -0.5
!     
      CORN(28,1,3) =  0.5
      CORN(28,2,3) = -0.5
      CORN(28,3,3) =  0.5
!     
      CORN(28,1,4) =  0.5
      CORN(28,2,4) =  0.5
      CORN(28,3,4) =  0.5
!     
!     - face 29
!     
      CORN(29,1,1) = -1.0
      CORN(29,2,1) =  0.5
      CORN(29,3,1) = -0.5
!     
      CORN(29,1,2) = -0.5
      CORN(29,2,2) =  0.5
      CORN(29,3,2) = -0.5
!     
      CORN(29,1,3) = -1.0
      CORN(29,2,3) =  0.5
      CORN(29,3,3) =  0.5
!     
      CORN(29,1,4) = -0.5
      CORN(29,2,4) =  0.5
      CORN(29,3,4) =  0.5
!     
!     - face 30
!     
      CORN(30,1,1) = -0.5
      CORN(30,2,1) =  0.5
      CORN(30,3,1) = -0.5
!     
      CORN(30,1,2) =  0.5
      CORN(30,2,2) =  0.5
      CORN(30,3,2) = -0.5
!     
      CORN(30,1,3) = -0.5
      CORN(30,2,3) =  0.5
      CORN(30,3,3) =  0.5
!     
      CORN(30,1,4) =  0.5
      CORN(30,2,4) =  0.5
      CORN(30,3,4) =  0.5
!     
!     - face 31
!     
      CORN(31,1,1) =  0.5 
      CORN(31,2,1) =  0.5 
      CORN(31,3,1) = -0.5
!     
      CORN(31,1,2) =  1.0
      CORN(31,2,2) =  0.5
      CORN(31,3,2) = -0.5
!     
      CORN(31,1,3) =  0.5
      CORN(31,2,3) =  0.5
      CORN(31,3,3) =  0.5
!     
      CORN(31,1,4) =  1.0
      CORN(31,2,4) =  0.5
      CORN(31,3,4) =  0.5
!     
!     - face 32
!     
      CORN(32,1,1) = -0.5
      CORN(32,2,1) =  0.5
      CORN(32,3,1) = -0.5
!     
      CORN(32,1,2) = -0.5
      CORN(32,2,2) =  1.0
      CORN(32,3,2) = -0.5
!     
      CORN(32,1,3) = -0.5
      CORN(32,2,3) =  0.5
      CORN(32,3,3) =  0.5
!     
      CORN(32,1,4) = -0.5
      CORN(32,2,4) =  1.0
      CORN(32,3,4) =  0.5
!     
!     - face 33
!     
      CORN(33,1,1) =  0.5
      CORN(33,2,1) =  0.5
      CORN(33,3,1) = -0.5
!     
      CORN(33,1,2) =  0.5
      CORN(33,2,2) =  1.0
      CORN(33,3,2) = -0.5
!     
      CORN(33,1,3) =  0.5
      CORN(33,2,3) =  0.5
      CORN(33,3,3) =  0.5
!     
      CORN(33,1,4) =  0.5
      CORN(33,2,4) =  1.0
      CORN(33,3,4) =  0.5
!     
!     - face 34
!     
      CORN(34,1,1) = -1.0
      CORN(34,2,1) = -1.0
      CORN(34,3,1) =  0.5
!     
      CORN(34,1,2) = -0.5
      CORN(34,2,2) = -1.0
      CORN(34,3,2) =  0.5
!     
      CORN(34,1,3) = -1.0
      CORN(34,2,3) = -0.5
      CORN(34,3,3) =  0.5
!     
      CORN(34,1,4) = -0.5
      CORN(34,2,4) = -0.5
      CORN(34,3,4) =  0.5
!     
!     - face 35
!     
      CORN(35,1,1) = -0.5
      CORN(35,2,1) = -1.0
      CORN(35,3,1) =  0.5
!     
      CORN(35,1,2) =  0.5
      CORN(35,2,2) = -1.0
      CORN(35,3,2) =  0.5
!     
      CORN(35,1,3) = -0.5
      CORN(35,2,3) = -0.5
      CORN(35,3,3) =  0.5
!     
      CORN(35,1,4) =  0.5
      CORN(35,2,4) = -0.5
      CORN(35,3,4) =  0.5
!     
!     - face 36
!     
      CORN(36,1,1) =  0.5
      CORN(36,2,1) = -1.0
      CORN(36,3,1) =  0.5
!     
      CORN(36,1,2) =  1.0
      CORN(36,2,2) = -1.0
      CORN(36,3,2) =  0.5
!     
      CORN(36,1,3) =  0.5
      CORN(36,2,3) = -0.5
      CORN(36,3,3) =  0.5
!     
      CORN(36,1,4) =  1.0
      CORN(36,2,4) = -0.5
      CORN(36,3,4) =  0.5
!     
!     - face 37
!     
      CORN(37,1,1) = -1.0
      CORN(37,2,1) = -0.5
      CORN(37,3,1) =  0.5
!     
      CORN(37,1,2) = -0.5
      CORN(37,2,2) = -0.5
      CORN(37,3,2) =  0.5
!     
      CORN(37,1,3) = -1.0
      CORN(37,2,3) =  0.5
      CORN(37,3,3) =  0.5
!     
      CORN(37,1,4) = -0.5
      CORN(37,2,4) =  0.5
      CORN(37,3,4) =  0.5
!     
!     - face 38
!     
      CORN(38,1,1) = -0.5
      CORN(38,2,1) = -0.5
      CORN(38,3,1) =  0.5
!     
      CORN(38,1,2) =  0.5
      CORN(38,2,2) = -0.5
      CORN(38,3,2) =  0.5
!     
      CORN(38,1,3) = -0.5
      CORN(38,2,3) =  0.5
      CORN(38,3,3) =  0.5
!     
      CORN(38,1,4) =  0.5
      CORN(38,2,4) =  0.5
      CORN(38,3,4) =  0.5
!     
!     - face 39
!     
      CORN(39,1,1) =  0.5
      CORN(39,2,1) = -0.5
      CORN(39,3,1) =  0.5
!     
      CORN(39,1,2) =  1.0
      CORN(39,2,2) = -0.5
      CORN(39,3,2) =  0.5
!     
      CORN(39,1,3) =  0.5
      CORN(39,2,3) =  0.5
      CORN(39,3,3) =  0.5
!     
      CORN(39,1,4) =  1.0
      CORN(39,2,4) =  0.5
      CORN(39,3,4) =  0.5
!     
!     - face 40
!     
      CORN(40,1,1) = -1.0
      CORN(40,2,1) =  0.5
      CORN(40,3,1) =  0.5
!     
      CORN(40,1,2) = -0.5
      CORN(40,2,2) =  0.5
      CORN(40,3,2) =  0.5
!     
      CORN(40,1,3) = -1.0
      CORN(40,2,3) =  1.0
      CORN(40,3,3) =  0.5
!     
      CORN(40,1,4) = -0.5
      CORN(40,2,4) =  1.0
      CORN(40,3,4) =  0.5
!     
!     - face 41
!     
      CORN(41,1,1) = -0.5
      CORN(41,2,1) =  0.5
      CORN(41,3,1) =  0.5
!     
      CORN(41,1,2) =  0.5
      CORN(41,2,2) =  0.5
      CORN(41,3,2) =  0.5
!     
      CORN(41,1,3) = -0.5
      CORN(41,2,3) =  1.0
      CORN(41,3,3) =  0.5
!     
      CORN(41,1,4) =  0.5
      CORN(41,2,4) =  1.0
      CORN(41,3,4) =  0.5
!     
!     - face 42
!     
      CORN(42,1,1) =  0.5
      CORN(42,2,1) =  0.5
      CORN(42,3,1) =  0.5
!     
      CORN(42,1,2) =  1.0
      CORN(42,2,2) =  0.5
      CORN(42,3,2) =  0.5
!     
      CORN(42,1,3) =  0.5
      CORN(42,2,3) =  1.0
      CORN(42,3,3) =  0.5
!     
      CORN(42,1,4) =  1.0
      CORN(42,2,4) =  1.0
      CORN(42,3,4) =  0.5
!     
!     - face 43
!     
      CORN(43,1,1) = -0.5
      CORN(43,2,1) = -1.0
      CORN(43,3,1) =  0.5
!     
      CORN(43,1,2) = -0.5
      CORN(43,2,2) = -0.5
      CORN(43,3,2) =  0.5
!     
      CORN(43,1,3) = -0.5
      CORN(43,2,3) = -1.0
      CORN(43,3,3) =  1.0
!     
      CORN(43,1,4) = -0.5
      CORN(43,2,4) = -0.5
      CORN(43,3,4) =  1.0
!     
!     - face 44
!     
      CORN(44,1,1) =  0.5
      CORN(44,2,1) = -1.0
      CORN(44,3,1) =  0.5
!     
      CORN(44,1,2) =  0.5
      CORN(44,2,2) = -0.5
      CORN(44,3,2) =  0.5
!     
      CORN(44,1,3) =  0.5
      CORN(44,2,3) = -1.0
      CORN(44,3,3) =  1.0
!     
      CORN(44,1,4) =  0.5
      CORN(44,2,4) = -0.5
      CORN(44,3,4) =  1.0
!     
!     - face 45
!     
      CORN(45,1,1) = -1.0
      CORN(45,2,1) = -0.5
      CORN(45,3,1) =  0.5
!     
      CORN(45,1,2) = -0.5
      CORN(45,2,2) = -0.5
      CORN(45,3,2) =  0.5
!     
      CORN(45,1,3) = -1.0
      CORN(45,2,3) = -0.5
      CORN(45,3,3) =  1.0
!     
      CORN(45,1,4) = -0.5
      CORN(45,2,4) = -0.5
      CORN(45,3,4) =  1.0
!     
!     - face 46
!     
      CORN(46,1,1) = -0.5
      CORN(46,2,1) = -0.5
      CORN(46,3,1) =  0.5
!     
      CORN(46,1,2) =  0.5
      CORN(46,2,2) = -0.5
      CORN(46,3,2) =  0.5
!     
      CORN(46,1,3) = -0.5
      CORN(46,2,3) = -0.5
      CORN(46,3,3) =  1.0
!     
      CORN(46,1,4) =  0.5
      CORN(46,2,4) = -0.5
      CORN(46,3,4) =  1.0
!     
!     - face 47
!     
      CORN(47,1,1) =  0.5
      CORN(47,2,1) = -0.5
      CORN(47,3,1) =  0.5
!     
      CORN(47,1,2) =  1.0
      CORN(47,2,2) = -0.5
      CORN(47,3,2) =  0.5
!     
      CORN(47,1,3) =  0.5
      CORN(47,2,3) = -0.5
      CORN(47,3,3) =  1.0
!     
      CORN(47,1,4) =  1.0
      CORN(47,2,4) = -0.5
      CORN(47,3,4) =  1.0
!     
!     - face 48
!     
      CORN(48,1,1) = -0.5
      CORN(48,2,1) = -0.5
      CORN(48,3,1) =  0.5
!     
      CORN(48,1,2) = -0.5
      CORN(48,2,2) =  0.5
      CORN(48,3,2) =  0.5
!     
      CORN(48,1,3) = -0.5
      CORN(48,2,3) = -0.5
      CORN(48,3,3) =  1.0
!     
      CORN(48,1,4) = -0.5
      CORN(48,2,4) =  0.5
      CORN(48,3,4) =  1.0
!     
!     - face 49
!     
      CORN(49,1,1) =  0.5
      CORN(49,2,1) = -0.5
      CORN(49,3,1) =  0.5
!     
      CORN(49,1,2) =  0.5
      CORN(49,2,2) =  0.5
      CORN(49,3,2) =  0.5
!     
      CORN(49,1,3) =  0.5
      CORN(49,2,3) = -0.5
      CORN(49,3,3) =  1.0
!     
      CORN(49,1,4) =  0.5
      CORN(49,2,4) =  0.5
      CORN(49,3,4) =  1.0
!     
!     - face 50
!     
      CORN(50,1,1) = -1.0
      CORN(50,2,1) =  0.5
      CORN(50,3,1) =  0.5
!     
      CORN(50,1,2) = -0.5
      CORN(50,2,2) =  0.5
      CORN(50,3,2) =  0.5
!     
      CORN(50,1,3) = -1.0
      CORN(50,2,3) =  0.5
      CORN(50,3,3) =  1.0
!     
      CORN(50,1,4) = -0.5
      CORN(50,2,4) =  0.5
      CORN(50,3,4) =  1.0
!     
!     - face 51
!     
      CORN(51,1,1) = -0.5
      CORN(51,2,1) =  0.5
      CORN(51,3,1) =  0.5
!     
      CORN(51,1,2) =  0.5
      CORN(51,2,2) =  0.5
      CORN(51,3,2) =  0.5
!     
      CORN(51,1,3) = -0.5
      CORN(51,2,3) =  0.5
      CORN(51,3,3) =  1.0
!     
      CORN(51,1,4) =  0.5
      CORN(51,2,4) =  0.5
      CORN(51,3,4) =  1.0
!     
!     - face 52
!     
      CORN(52,1,1) =  0.5
      CORN(52,2,1) =  0.5
      CORN(52,3,1) =  0.5
!     
      CORN(52,1,2) =  1.0
      CORN(52,2,2) =  0.5
      CORN(52,3,2) =  0.5
!     
      CORN(52,1,3) =  0.5
      CORN(52,2,3) =  0.5
      CORN(52,3,3) =  1.0
!     
      CORN(52,1,4) =  1.0
      CORN(52,2,4) =  0.5
      CORN(52,3,4) =  1.0
!     
!     - face 53
!     
      CORN(53,1,1) = -0.5
      CORN(53,2,1) =  0.5
      CORN(53,3,1) =  0.5
!     
      CORN(53,1,2) = -0.5
      CORN(53,2,2) =  1.0
      CORN(53,3,2) =  0.5
!     
      CORN(53,1,3) = -0.5
      CORN(53,2,3) =  0.5
      CORN(53,3,3) =  1.0
!     
      CORN(53,1,4) = -0.5
      CORN(53,2,4) =  1.0
      CORN(53,3,4) =  1.0
!     
!     - face 54
!     
      CORN(54,1,1) =  0.5
      CORN(54,2,1) =  0.5
      CORN(54,3,1) =  0.5
!     
      CORN(54,1,2) =  0.5
      CORN(54,2,2) =  1.0
      CORN(54,3,2) =  0.5
!     
      CORN(54,1,3) =  0.5
      CORN(54,2,3) =  0.5
      CORN(54,3,3) =  1.0
!     
      CORN(54,1,4) =  0.5
      CORN(54,2,4) =  1.0
      CORN(54,3,4) =  1.0
!     
!     - define the Gaussian integration points in XIP space.
!     
      XIPGP(1)  = -1.0/SQRT(3.0)  
      ETAPGP(1) = -1.0/SQRT(3.0)  
!     
      XIPGP(2)  =  1.0/SQRT(3.0)  
      ETAPGP(2) = -1.0/SQRT(3.0)  
!     
      XIPGP(3)  = -1.0/SQRT(3.0)  
      ETAPGP(3) =  1.0/SQRT(3.0)  
!     
      XIPGP(4)  =  1.0/SQRT(3.0)  
      ETAPGP(4) =  1.0/SQRT(3.0)  
!     
!     - define positions of vertices of quadrilateral element
!     - associated with a subcell face in XIP
!     
      XIP(1)  = -1.0
      ETAP(1) = -1.0
!     
      XIP(2)  =  1.0
      ETAP(2) = -1.0
!     
      XIP(3)  = -1.0
      ETAP(3) =  1.0
!     
      XIP(4)  =  1.0
      ETAP(4) =  1.0
!     
!     - generate values of the FE basis functions at the quadrature
!     - points on the faces of the subcells.
!     
      do IFACE = 1, NFACE! Was loop 
!     
      do GJ = 1, FNGI! Was loop 
!     
            GI = (IFACE-1)*FNGI+GJ
!     
      do ICOORD = 1, NCOORD! Was loop 
!     
               POS(ICOORD)    = 0.0
               DPDXI(ICOORD)  = 0.0
               DPDETA(ICOORD) = 0.0
!     
      do JLOC = 1, FNLOC! Was loop 
!     
                  POS(ICOORD) = POS(ICOORD)&
     &                 + CORN(IFACE,ICOORD,JLOC)&
     &                 *0.25*(1.+XIP(JLOC)*XIPGP(GJ))&
     &                 *(1.+ETAP(JLOC)*ETAPGP(GJ))
!     
                  DPDXI(ICOORD) = DPDXI(ICOORD)&
     &                 + CORN(IFACE,ICOORD,JLOC)&
     &                 *0.25*XIP(JLOC)&
     &                 *(1.+ETAP(JLOC)*ETAPGP(GJ))
!     
                  DPDETA(ICOORD) = DPDETA(ICOORD)&
     &                 + CORN(IFACE,ICOORD,JLOC)&
     &                 *0.25*(1.+XIP(JLOC)*XIPGP(GJ))&
     &                 *ETAP(JLOC)
!     
               END DO
!     
            END DO
!     
            XI   = POS(1)
            ETA  = POS(2)
            ZETA = POS(3)
!     
            LIJXI(1)   =  0.5*XI*(XI-1.0)
            LIJXI(2)   =  1.0-XI*XI
            LIJXI(3)   =  0.5*XI*(XI+1.0)
!     
            LIJETA(1)  =  0.5*ETA*(ETA-1.0)
            LIJETA(2)  =  1.0-ETA*ETA
            LIJETA(3)  =  0.5*ETA*(ETA+1.0)
!     
            LIJZETA(1) =  0.5*ZETA*(ZETA-1.0)
            LIJZETA(2) =  1.0-ZETA*ZETA
            LIJZETA(3) =  0.5*ZETA*(ZETA+1.0)
!     
            DLIJXIDXI(1)   =  0.5*(2.0*XI-1.0)*DPDXI(1)
            DLIJXIDXI(2)   = -2.0*XI*DPDXI(1)
            DLIJXIDXI(3)   =  0.5*(2.0*XI+1.0)*DPDXI(1)
!     
            DLIJETADXI(1)  =  0.5*(2.0*ETA-1.0)*DPDXI(2)
            DLIJETADXI(2)  = -2.0*ETA*DPDXI(2)
            DLIJETADXI(3)  =  0.5*(2.0*ETA+1.0)*DPDXI(2)
!     
            DLIJZETADXI(1) =  0.5*(2.0*ZETA-1.0)*DPDXI(3)
            DLIJZETADXI(2) = -2.0*ZETA*DPDXI(3)
            DLIJZETADXI(3) =  0.5*(2.0*ZETA+1.0)*DPDXI(3)
!     
            DLIJXIDETA(1)  =  0.5*(2.0*XI-1.0)*DPDETA(1)
            DLIJXIDETA(2)  = -2.0*XI*DPDETA(1)
            DLIJXIDETA(3)  =  0.5*(2.0*XI+1.0)*DPDETA(1)
!     
            DLIJETADETA(1) =  0.5*(2.0*ETA-1.0)*DPDETA(2)
            DLIJETADETA(2) = -2.0*ETA*DPDETA(2)
            DLIJETADETA(3) =  0.5*(2.0*ETA+1.0)*DPDETA(2)
!     
            DLIJZETADETA(1) =  0.5*(2.0*ZETA-1.0)*DPDETA(3)
            DLIJZETADETA(2) = -2.0*ZETA*DPDETA(3)
            DLIJZETADETA(3) =  0.5*(2.0*ZETA+1.0)*DPDETA(3)
!     
      do I = 1, 3! Was loop 
!     
      do J = 1, 3! Was loop 
!     
      do K = 1, 3! Was loop 
!     
                     ILOC = I+(J-1)*3+(K-1)*9
!     
                     SVN(ILOC,GI)   = LIJXI(I)*LIJETA(J)*LIJZETA(K)
!     
                     SVNLX(ILOC,GI) = LIJXI(I)*LIJETA(J)*DLIJZETADXI(K)&
     &                    + LIJXI(I)*DLIJETADXI(J)*LIJZETA(K)&
     &                    + DLIJXIDXI(I)*LIJETA(J)*LIJZETA(K)
!     
                     SVNLY(ILOC,GI) = LIJXI(I)*LIJETA(J)*DLIJZETADETA(K)&
     &                    + LIJXI(I)*DLIJETADETA(J)*LIJZETA(K)&
     &                    + DLIJXIDETA(I)*LIJETA(J)*LIJZETA(K)
!     
                  END DO
!     
               END DO
!     
            END DO
!     
         END DO
!     
      END DO
!     
!     - set the weights for the surface integration
!     
      do GI = 1, SVNGI! Was loop 
!     
         SVWEIGH(GI) = 1.0
!     
      END DO
!     
!     stop 11
      END SUBROUTINE

      SUBROUTINE FVQTRI( NGI,    NLOC, SVNGI,&
     &     M,      SVN,  SVNLX,  &
     &     SVWEIGH              )
!     ----------------------------------------
!     
!     - this routine generates the shape functions associated
!     - with the FV's i.e. their surfaces and volume shape 
!     - functions and derivatives.
!     
!     -------------------------------
!     - date last modified : 06/10/2002
!     -------------------------------
!     
      IMPLICIT NONE
!     
      INTEGER NGI, NLOC, SVNGI
!     
      REAL M(NLOC,NGI),       SVN(NLOC,SVNGI)
      REAL SVNLX(NLOC,SVNGI), SVWEIGH(SVNGI)
!     
!     - local variables
!     
      INTEGER ILOC, JLOC, IFACE, ICOORD
!     
      INTEGER COUNT, GI, GJ
!     
!     - note that NCOORD is the number of co-ordinates.
!     - NFACE is the number of surfaces internal to an 
!     - element 
!     
      INTEGER NCOORD, NFACE
      PARAMETER( NCOORD = 3, NFACE = 9 )
!     
      INTEGER INC, FNGI, FNLOC
      PARAMETER( FNGI = 2, FNLOC = 2  )
!     
      REAL POS(NCOORD),    DPDXI(NCOORD)
!     
      REAL XIPGP(FNGI)
!     
      REAL XIP(FNLOC)
!     
      REAL CORN(NFACE,NCOORD,FNLOC)
!     
      REAL L1, L2, L3
!     
!     - FV basis functions for use in calculating the integral of the
!     - shape functions over a subcell of a CV.
!     
      M = 0.0
!     
      COUNT = 1
!     
      do ILOC = 1, NLOC! Was loop 
!     
         IF((ILOC.EQ.1).OR.(ILOC.EQ.3).OR.(ILOC.EQ.5)) THEN
            INC = 3
         ELSE
            INC = 7
         END IF
!     
      do GI = COUNT, COUNT+INC! Was loop 
!     
            M(ILOC,GI) = 1.0
!     
         END DO
!     
         IF((ILOC.EQ.1).OR.(ILOC.EQ.3).OR.(ILOC.EQ.5)) THEN
            COUNT = COUNT + 4
         ELSE
            COUNT = COUNT + 8
         END IF
!     
      END DO 
!     
!     - Note that CORN(IFACE,ICOORD,JLOC) tabulates the positions
!     - of the corners (vertices) of the faces of a subcell in volume 
!     - co-ordinates. Note also that IFACE ranges from 1-3, ICOORD
!     - ranges from 1-3 and JLOC ranges from 1-2. IFACE signifies 
!     - the particular subcell face, ICOORD signfies the co-ordinates
!     - of the vertices of the face in volume co-ordinates (L1,L2,L3)
!     - and JLOC signifies the vertex of the face of the subcell
!     - which are straight lines.
!     
!     - face 1
!     
      CORN(1,1,1) = 0.75
      CORN(1,2,1) = 0.25
      CORN(1,3,1) = 0.0
!     
      CORN(1,1,2) = 0.66666666
      CORN(1,2,2) = 0.16666666
      CORN(1,3,2) = 0.16666666
!     
!     - face 2
!     
      CORN(2,1,1) = 0.66666666
      CORN(2,2,1) = 0.16666666
      CORN(2,3,1) = 0.16666666
!     
      CORN(2,1,2) = 0.75
      CORN(2,2,2) = 0.0
      CORN(2,3,2) = 0.25
!     
!     - face 3
!     
      CORN(3,1,1) = 0.25
      CORN(3,2,1) = 0.75
      CORN(3,3,1) = 0.0
!     
      CORN(3,1,2) = 0.16666666
      CORN(3,2,2) = 0.66666666
      CORN(3,3,2) = 0.16666666
!     
!     - face 4
!     
      CORN(4,1,1) = 0.16666666
      CORN(4,2,1) = 0.66666666
      CORN(4,3,1) = 0.16666666
!     
      CORN(4,1,2) = 0.0
      CORN(4,2,2) = 0.75
      CORN(4,3,2) = 0.25
!     
!     - face 5
!     
      CORN(5,1,1) = 0.0
      CORN(5,2,1) = 0.25
      CORN(5,3,1) = 0.75
!     
      CORN(5,1,2) = 0.16666666
      CORN(5,2,2) = 0.16666666
      CORN(5,3,2) = 0.66666666
!     
!     - face 6
!     
      CORN(6,1,1) = 0.25
      CORN(6,2,1) = 0.0
      CORN(6,3,1) = 0.75
!     
      CORN(6,1,2) = 0.16666666
      CORN(6,2,2) = 0.16666666
      CORN(6,3,2) = 0.66666666
!     
!     - face 7
!     
      CORN(7,1,1) = 0.66666666
      CORN(7,2,1) = 0.16666666
      CORN(7,3,1) = 0.16666666
!     
      CORN(7,1,2) = 0.33333333
      CORN(7,2,2) = 0.33333333
      CORN(7,3,2) = 0.33333333
!     
!     - face 8
!     
      CORN(8,1,1) = 0.33333333
      CORN(8,2,1) = 0.33333333
      CORN(8,3,1) = 0.33333333
!     
      CORN(8,1,2) = 0.16666666
      CORN(8,2,2) = 0.66666666
      CORN(8,3,2) = 0.16666666
!     
!     - face 9
!     
      CORN(9,1,1) = 0.33333333
      CORN(9,2,1) = 0.33333333
      CORN(9,3,1) = 0.33333333
!     
      CORN(9,1,2) = 0.16666666
      CORN(9,2,2) = 0.16666666
      CORN(9,3,2) = 0.66666666
!     
!     - define the Gaussian integration points in XIP space.
!     
      XIPGP(1)  = -1/SQRT(3.0) 
      XIPGP(2)  =  1/SQRT(3.0) 
!     
!     - define positions of vertices of line element
!     - associated with a subcell face in XIP
!     
      XIP(1) = -1.0
      XIP(2) =  1.0
!     
!     - generate values of the FE basis functions at the quadrature
!     - points on the faces of the subcells.
!     
      do IFACE = 1, NFACE! Was loop 
!     
      do GJ = 1, FNGI! Was loop 
!     
            GI = (IFACE-1)*FNGI+GJ
!     
      do ICOORD = 1, NCOORD! Was loop 
!     
               POS(ICOORD)   = 0.0
               DPDXI(ICOORD) = 0.0
!     
      do JLOC = 1, FNLOC! Was loop 
!     
                  POS(ICOORD) = POS(ICOORD)&
     &                 + CORN(IFACE,ICOORD,JLOC)&
     &                 *0.5*(1.+XIP(JLOC)*XIPGP(GJ))
!     
                  DPDXI(ICOORD) = DPDXI(ICOORD)&
     &                 + CORN(IFACE,ICOORD,JLOC)&
     &                 *0.5*XIP(JLOC)
!     
               END DO
!     
            END DO
!     
            L1 = POS(1)
            L2 = POS(2)
            L3 = POS(3)
!     
            SVN(1,GI) = (2.0*L1-1.0)*L1
            SVN(2,GI) =  4.0*L1*L2
            SVN(3,GI) = (2.0*L2-1.0)*L2
            SVN(4,GI) =  4.0*L2*L3
            SVN(5,GI) = (2.0*L3-1.0)*L3
            SVN(6,GI) =  4.0*L1*L3
!     
            SVNLX(1,GI) = (4.0*L1-1.0)*DPDXI(1)
            SVNLX(2,GI) =  4.0*DPDXI(1)*L2 + 4.0*L1*DPDXI(2)
            SVNLX(3,GI) = (4.0*L2-1.0)*DPDXI(2)
            SVNLX(4,GI) =  4.0*DPDXI(2)*L3 + 4.0*L2*DPDXI(3)
            SVNLX(5,GI) = (4.0*L3-1.0)*DPDXI(3)
            SVNLX(6,GI) =  4.0*DPDXI(1)*L3 + 4.0*L1*DPDXI(3)
!     
         END DO
!     
      END DO
!     
      do GI = 1, SVNGI! Was loop 
!     
         SVWEIGH(GI) = 1.0
!     
      END DO
!     
      END SUBROUTINE

      SUBROUTINE FVQTET( NGI,   NLOC,   SVNGI,&
!     - REALS            
     &     M,     SVN,    SVNLX,&
     &     SVNLY, SVWEIGH        )
!     ------------------------------------------
!     
!     - this routine generates the shape functions associated
!     - with the FV's i.e. their surfaces and volume shape 
!     - functions and derivatives. The surface shape functions
!     - are the values of the FE volume shape functions evaluated
!     - on the surfaces of the CV's.
!     
!     -------------------------------
!     - date last modified : 07/11/2002
!     -------------------------------
!     
      IMPLICIT NONE
!     
      INTEGER NGI, NLOC, SVNGI
!     
      REAL M(NLOC,NGI),     SVWEIGH(SVNGI)
!     
      REAL SVN(NLOC,SVNGI), SVNLX(NLOC,SVNGI)
!     
      REAL SVNLY(NLOC,SVNGI)
!     
!     - local variables
!     
      INTEGER ILOC, JLOC, IFACE, ICOORD
!     
      INTEGER COUNT, GI, GJ
!     
!     - note that NCOORD is the number of co-ordinates.
!     - NFACE is the number of surfaces internal to an 
!     - element 
!     
      INTEGER NCOORD, NFACE
      PARAMETER( NCOORD = 4, NFACE = 24 )
!     
      INTEGER FNGI, FNLOC
      PARAMETER( FNGI = 4, FNLOC = 4  )
!     
      REAL POS(NCOORD),    DPDXI(NCOORD)
      REAL DPDETA(NCOORD)
!     
      REAL XIPGP(FNGI), ETAPGP(FNGI)
!     
      REAL XIP(FNLOC), ETAP(FNLOC)
!     
      REAL CORN(NFACE,NCOORD,FNLOC)
!     
      REAL L1, L2, L3, L4
!     
      INTEGER INC
!     
!     - FV basis functions for use in calculating the integral of the
!     - shape functions over a subcell of a CV.
!     
      M = 0.0
!     
      COUNT = 1
!     
      do ILOC = 1, NLOC! Was loop 
!     
         IF((ILOC.EQ.1).OR.(ILOC.EQ.3).OR.&
     &        (ILOC.EQ.5).OR.(ILOC.EQ.10)) THEN
            INC = 7
         ELSE
            INC = 15
         END IF
!     
      do GI = COUNT, COUNT + INC! Was loop 
!     
            M(ILOC,GI) = 1.0
!     
         END DO
!     
         IF((ILOC.EQ.1).OR.(ILOC.EQ.3).OR.&
     &        (ILOC.EQ.5).OR.(ILOC.EQ.10)) THEN
            COUNT = COUNT + 8
         ELSE
            COUNT = COUNT + 16
         END IF
!     
      END DO
!     
!     - Note that CORN(IFACE,ICOORD,JLOC) tabulates the positions
!     - of the corners (vertices) of the faces of a subcell in volume 
!     - co-ordinates. Note also that IFACE ranges from 1-24, ICOORD
!     - ranges from 1-4 and JLOC ranges from 1-4. IFACE signifies 
!     - the particular subcell face, ICOORD signfies the co-ordinates
!     - of the vertices of the face in volume co-ordinates (L1,L2,L3,L4)
!     - and JLOC signifies the vertex of the face of the subcell
!     - which are straight lines.
!     
!     - face 1
!     
      CORN(1,1,1) = 0.75
      CORN(1,2,1) = 0.25
      CORN(1,3,1) = 0.0
      CORN(1,4,1) = 0.0
!     
      CORN(1,1,2) = 0.66666666
      CORN(1,2,2) = 0.16666666
      CORN(1,3,2) = 0.16666666
      CORN(1,4,2) = 0.0
!     
      CORN(1,1,3) = 0.66666666
      CORN(1,2,3) = 0.16666666
      CORN(1,3,3) = 0.0
      CORN(1,4,3) = 0.16666666
!     
      CORN(1,1,4) = 0.625
      CORN(1,2,4) = 0.125
      CORN(1,3,4) = 0.125
      CORN(1,4,4) = 0.125
!     
!     - face 2
!     
      CORN(2,1,1) = 0.66666666
      CORN(2,2,1) = 0.16666666
      CORN(2,3,1) = 0.16666666
      CORN(2,4,1) = 0.0
!     
      CORN(2,1,2) = 0.75
      CORN(2,2,2) = 0.0
      CORN(2,3,2) = 0.25
      CORN(2,4,2) = 0.0
!     
      CORN(2,1,3) = 0.625
      CORN(2,2,3) = 0.125
      CORN(2,3,3) = 0.125
      CORN(2,4,3) = 0.125
!     
      CORN(2,1,4) = 0.66666666
      CORN(2,2,4) = 0.0
      CORN(2,3,4) = 0.16666666
      CORN(2,4,4) = 0.16666666
!     
!     - face 3
!     
      CORN(3,1,1) = 0.75
      CORN(3,2,1) = 0.0
      CORN(3,3,1) = 0.0
      CORN(3,4,1) = 0.25
!     
      CORN(3,1,2) = 0.66666666
      CORN(3,2,2) = 0.16666666
      CORN(3,3,2) = 0.0
      CORN(3,4,2) = 0.16666666
!     
      CORN(3,1,3) = 0.66666666
      CORN(3,2,3) = 0.0
      CORN(3,3,3) = 0.16666666
      CORN(3,4,3) = 0.16666666
!     
      CORN(3,1,4) = 0.625
      CORN(3,2,4) = 0.125
      CORN(3,3,4) = 0.125
      CORN(3,4,4) = 0.125
!     
!     - face 4
!     
      CORN(4,1,1) = 0.0
      CORN(4,2,1) = 0.25
      CORN(4,3,1) = 0.75
      CORN(4,4,1) = 0.0
!     
      CORN(4,1,2) = 0.16666666
      CORN(4,2,2) = 0.16666666
      CORN(4,3,2) = 0.66666666
      CORN(4,4,2) = 0.0
!     
      CORN(4,1,3) = 0.0
      CORN(4,2,3) = 0.16666666
      CORN(4,3,3) = 0.66666666
      CORN(4,4,3) = 0.16666666
!     
      CORN(4,1,4) = 0.125
      CORN(4,2,4) = 0.125
      CORN(4,3,4) = 0.625
      CORN(4,4,4) = 0.125
!     
!     - face 5
!     
      CORN(5,1,1) = 0.25
      CORN(5,2,1) = 0.0
      CORN(5,3,1) = 0.75
      CORN(5,4,1) = 0.0
!     
      CORN(5,1,2) = 0.16666666
      CORN(5,2,2) = 0.16666666
      CORN(5,3,2) = 0.66666666
      CORN(5,4,2) = 0.0
!     
      CORN(5,1,3) = 0.16666666
      CORN(5,2,3) = 0.0
      CORN(5,3,3) = 0.66666666
      CORN(5,4,3) = 0.16666666
!     
      CORN(5,1,4) = 0.125
      CORN(5,2,4) = 0.125
      CORN(5,3,4) = 0.625
      CORN(5,4,4) = 0.125
!     
!     - face 6
!     
      CORN(6,1,1) = 0.0
      CORN(6,2,1) = 0.0
      CORN(6,3,1) = 0.75
      CORN(6,4,1) = 0.25
!     
      CORN(6,1,2) = 0.16666666
      CORN(6,2,2) = 0.0
      CORN(6,3,2) = 0.66666666
      CORN(6,4,2) = 0.16666666
!     
      CORN(6,1,3) = 0.0
      CORN(6,2,3) = 0.16666666
      CORN(6,3,3) = 0.66666666
      CORN(6,4,3) = 0.16666666
!     
      CORN(6,1,4) = 0.125
      CORN(6,2,4) = 0.125
      CORN(6,3,4) = 0.625
      CORN(6,4,4) = 0.125
!     
!     - face 7
!     
      CORN(7,1,1) = 0.25
      CORN(7,2,1) = 0.75
      CORN(7,3,1) = 0.0
      CORN(7,4,1) = 0.0
!     
      CORN(7,1,2) = 0.16666666
      CORN(7,2,2) = 0.66666666
      CORN(7,3,2) = 0.16666666
      CORN(7,4,2) = 0.0
!     
      CORN(7,1,3) = 0.16666666
      CORN(7,2,3) = 0.66666666
      CORN(7,3,3) = 0.0
      CORN(7,4,3) = 0.16666666
!     
      CORN(7,1,4) = 0.125
      CORN(7,2,4) = 0.625
      CORN(7,3,4) = 0.125
      CORN(7,4,4) = 0.125
!     
!     - face 8
!     
      CORN(8,1,1) = 0.0
      CORN(8,2,1) = 0.75
      CORN(8,3,1) = 0.25
      CORN(8,4,1) = 0.0
!     
      CORN(8,1,2) = 0.16666666
      CORN(8,2,2) = 0.66666666
      CORN(8,3,2) = 0.16666666
      CORN(8,4,2) = 0.0
!     
      CORN(8,1,3) = 0.0
      CORN(8,2,3) = 0.66666666
      CORN(8,3,3) = 0.16666666
      CORN(8,4,3) = 0.16666666
!     
      CORN(8,1,4) = 0.125
      CORN(8,2,4) = 0.625
      CORN(8,3,4) = 0.125
      CORN(8,4,4) = 0.125
!     
!     - face 9
!     
      CORN(9,1,1) = 0.0
      CORN(9,2,1) = 0.75
      CORN(9,3,1) = 0.0
      CORN(9,4,1) = 0.25
!     
      CORN(9,1,2) = 0.0
      CORN(9,2,2) = 0.66666666
      CORN(9,3,2) = 0.16666666
      CORN(9,4,2) = 0.16666666
!     
      CORN(9,1,3) = 0.16666666
      CORN(9,2,3) = 0.66666666 
      CORN(9,3,3) = 0.0
      CORN(9,4,3) = 0.16666666 
!     
      CORN(9,1,4) = 0.125
      CORN(9,2,4) = 0.625
      CORN(9,3,4) = 0.125
      CORN(9,4,4) = 0.125
!     
!     - face 10
!     
      CORN(10,1,1) = 0.0
      CORN(10,2,1) = 0.25
      CORN(10,3,1) = 0.0
      CORN(10,4,1) = 0.75
!     
      CORN(10,1,2) = 0.16666666 
      CORN(10,2,2) = 0.16666666 
      CORN(10,3,2) = 0.0
      CORN(10,4,2) = 0.66666666 
!     
      CORN(10,1,3) = 0.0
      CORN(10,2,3) = 0.16666666 
      CORN(10,3,3) = 0.16666666 
      CORN(10,4,3) = 0.66666666 
!     
      CORN(10,1,4) = 0.125
      CORN(10,2,4) = 0.125
      CORN(10,3,4) = 0.125
      CORN(10,4,4) = 0.625
!     
!     - face 11
!     
      CORN(11,1,1) = 0.25
      CORN(11,2,1) = 0.0
      CORN(11,3,1) = 0.0
      CORN(11,4,1) = 0.75
!     
      CORN(11,1,2) = 0.16666666 
      CORN(11,2,2) = 0.16666666 
      CORN(11,3,2) = 0.0
      CORN(11,4,2) = 0.66666666 
!     
      CORN(11,1,3) = 0.16666666 
      CORN(11,2,3) = 0.0
      CORN(11,3,3) = 0.16666666 
      CORN(11,4,3) = 0.66666666 
!     
      CORN(11,1,4) = 0.125
      CORN(11,2,4) = 0.125
      CORN(11,3,4) = 0.125
      CORN(11,4,4) = 0.625
!     
!     - face 12
!     
      CORN(12,1,1) = 0.0
      CORN(12,2,1) = 0.0
      CORN(12,3,1) = 0.25
      CORN(12,4,1) = 0.75
!     
      CORN(12,1,2) = 0.0
      CORN(12,2,2) = 0.16666666 
      CORN(12,3,2) = 0.16666666 
      CORN(12,4,2) = 0.66666666 
!     
      CORN(12,1,3) = 0.16666666 
      CORN(12,2,3) = 0.0
      CORN(12,3,3) = 0.16666666 
      CORN(12,4,3) = 0.66666666 
!     
      CORN(12,1,4) = 0.125
      CORN(12,2,4) = 0.125
      CORN(12,3,4) = 0.125
      CORN(12,4,4) = 0.625
!     
!     - face 13
!     
      CORN(13,1,1) = 0.16666666 
      CORN(13,2,1) = 0.16666666 
      CORN(13,3,1) = 0.0
      CORN(13,4,1) = 0.66666666 
!     
      CORN(13,1,2) = 0.33333333
      CORN(13,2,2) = 0.33333333
      CORN(13,3,2) = 0.0
      CORN(13,4,2) = 0.33333333
!     
      CORN(13,1,3) = 0.125
      CORN(13,2,3) = 0.125
      CORN(13,3,3) = 0.125
      CORN(13,4,3) = 0.625
!     
      CORN(13,1,4) = 0.25
      CORN(13,2,4) = 0.25
      CORN(13,3,4) = 0.25
      CORN(13,4,4) = 0.25
!     
!     - face 14
!     
      CORN(14,1,1) = 0.0
      CORN(14,2,1) = 0.16666666 
      CORN(14,3,1) = 0.16666666 
      CORN(14,4,1) = 0.66666666 
!     
      CORN(14,1,2) = 0.0
      CORN(14,2,2) = 0.33333333
      CORN(14,3,2) = 0.33333333
      CORN(14,4,2) = 0.33333333
!     
      CORN(14,1,3) = 0.125
      CORN(14,2,3) = 0.125
      CORN(14,3,3) = 0.125
      CORN(14,4,3) = 0.625
!     
      CORN(14,1,4) = 0.25
      CORN(14,2,4) = 0.25
      CORN(14,3,4) = 0.25
      CORN(14,4,4) = 0.25
!     
!     - face 15
!     
      CORN(15,1,1) = 0.16666666 
      CORN(15,2,1) = 0.0
      CORN(15,3,1) = 0.16666666 
      CORN(15,4,1) = 0.66666666 
!     
      CORN(15,1,2) = 0.33333333
      CORN(15,2,2) = 0.0
      CORN(15,3,2) = 0.33333333
      CORN(15,4,2) = 0.33333333
!     
      CORN(15,1,3) = 0.125
      CORN(15,2,3) = 0.125
      CORN(15,3,3) = 0.125
      CORN(15,4,3) = 0.625
!     
      CORN(15,1,4) = 0.25
      CORN(15,2,4) = 0.25
      CORN(15,3,4) = 0.25
      CORN(15,4,4) = 0.25
!     
!     - face 16
!     
      CORN(16,1,1) = 0.66666666 
      CORN(16,2,1) = 0.16666666 
      CORN(16,3,1) = 0.0
      CORN(16,4,1) = 0.16666666 
!     
      CORN(16,1,2) = 0.33333333
      CORN(16,2,2) = 0.33333333
      CORN(16,3,2) = 0.0
      CORN(16,4,2) = 0.33333333
!     
      CORN(16,1,3) = 0.625
      CORN(16,2,3) = 0.125
      CORN(16,3,3) = 0.125
      CORN(16,4,3) = 0.125
!     
      CORN(16,1,4) = 0.25
      CORN(16,2,4) = 0.25
      CORN(16,3,4) = 0.25
      CORN(16,4,4) = 0.25
!     
!     - face 17
!     
      CORN(17,1,1) = 0.33333333
      CORN(17,2,1) = 0.33333333
      CORN(17,3,1) = 0.0
      CORN(17,4,1) = 0.33333333
!     
      CORN(17,1,2) = 0.16666666 
      CORN(17,2,2) = 0.66666666
      CORN(17,3,2) = 0.0
      CORN(17,4,2) = 0.16666666 
!     
      CORN(17,1,3) = 0.25
      CORN(17,2,3) = 0.25
      CORN(17,3,3) = 0.25
      CORN(17,4,3) = 0.25
!     
      CORN(17,1,4) = 0.125
      CORN(17,2,4) = 0.625
      CORN(17,3,4) = 0.125
      CORN(17,4,4) = 0.125
!     
!     - face 18
!     
      CORN(18,1,1) = 0.0
      CORN(18,2,1) = 0.66666666
      CORN(18,3,1) = 0.16666666 
      CORN(18,4,1) = 0.16666666 
!     
      CORN(18,1,2) = 0.0
      CORN(18,2,2) = 0.33333333
      CORN(18,3,2) = 0.33333333
      CORN(18,4,2) = 0.33333333
!     
      CORN(18,1,3) = 0.125
      CORN(18,2,3) = 0.625
      CORN(18,3,3) = 0.125
      CORN(18,4,3) = 0.125
!     
      CORN(18,1,4) = 0.25
      CORN(18,2,4) = 0.25
      CORN(18,3,4) = 0.25
      CORN(18,4,4) = 0.25
!     
!     - face 19
!     
      CORN(19,1,1) = 0.0
      CORN(19,2,1) = 0.33333333
      CORN(19,3,1) = 0.33333333
      CORN(19,4,1) = 0.33333333
!     
      CORN(19,1,2) = 0.0
      CORN(19,2,2) = 0.16666666
      CORN(19,3,2) = 0.66666666 
      CORN(19,4,2) = 0.16666666 
!     
      CORN(19,1,3) = 0.25
      CORN(19,2,3) = 0.25
      CORN(19,3,3) = 0.25
      CORN(19,4,3) = 0.25
!     
      CORN(19,1,4) = 0.125
      CORN(19,2,4) = 0.125
      CORN(19,3,4) = 0.625
      CORN(19,4,4) = 0.125
!     
!     - face 20
!     
      CORN(20,1,1) = 0.16666666
      CORN(20,2,1) = 0.0
      CORN(20,3,1) = 0.66666666 
      CORN(20,4,1) = 0.16666666
!     
      CORN(20,1,2) = 0.33333333
      CORN(20,2,2) = 0.0
      CORN(20,3,2) = 0.33333333
      CORN(20,4,2) = 0.33333333
!     
      CORN(20,1,3) = 0.125
      CORN(20,2,3) = 0.125
      CORN(20,3,3) = 0.625
      CORN(20,4,3) = 0.125
!     
      CORN(20,1,4) = 0.25
      CORN(20,2,4) = 0.25
      CORN(20,3,4) = 0.25
      CORN(20,4,4) = 0.25
!     
!     - face 21
!     
      CORN(21,1,1) = 0.33333333
      CORN(21,2,1) = 0.0
      CORN(21,3,1) = 0.33333333
      CORN(21,4,1) = 0.33333333
!     
      CORN(21,1,2) = 0.66666666
      CORN(21,2,2) = 0.0
      CORN(21,3,2) = 0.16666666
      CORN(21,4,2) = 0.16666666
!     
      CORN(21,1,3) = 0.25
      CORN(21,2,3) = 0.25
      CORN(21,3,3) = 0.25
      CORN(21,4,3) = 0.25
!     
      CORN(21,1,4) = 0.625
      CORN(21,2,4) = 0.125
      CORN(21,3,4) = 0.125
      CORN(21,4,4) = 0.125
!     
!     - face 22
!     
      CORN(22,1,1) = 0.16666666
      CORN(22,2,1) = 0.66666666
      CORN(22,3,1) = 0.16666666
      CORN(22,4,1) = 0.0
!     
      CORN(22,1,2) = 0.33333333
      CORN(22,2,2) = 0.33333333
      CORN(22,3,2) = 0.33333333
      CORN(22,4,2) = 0.0
!     
      CORN(22,1,3) = 0.125
      CORN(22,2,3) = 0.625
      CORN(22,3,3) = 0.125
      CORN(22,4,3) = 0.125
!     
      CORN(22,1,4) = 0.25
      CORN(22,2,4) = 0.25
      CORN(22,3,4) = 0.25
      CORN(22,4,4) = 0.25
!     
!     - face 23
!     
      CORN(23,1,1) = 0.66666666
      CORN(23,2,1) = 0.16666666
      CORN(23,3,1) = 0.16666666
      CORN(23,4,1) = 0.0
!     
      CORN(23,1,2) = 0.33333333
      CORN(23,2,2) = 0.33333333
      CORN(23,3,2) = 0.33333333
      CORN(23,4,2) = 0.0
!     
      CORN(23,1,3) = 0.625
      CORN(23,2,3) = 0.125
      CORN(23,3,3) = 0.125
      CORN(23,4,3) = 0.125
!     
      CORN(23,1,4) = 0.25
      CORN(23,2,4) = 0.25
      CORN(23,3,4) = 0.25
      CORN(23,4,4) = 0.25
!     
!     - face 24
!     
      CORN(24,1,1) = 0.33333333
      CORN(24,2,1) = 0.33333333
      CORN(24,3,1) = 0.33333333
      CORN(24,4,1) = 0.0
!     
      CORN(24,1,2) = 0.16666666
      CORN(24,2,2) = 0.16666666
      CORN(24,3,2) = 0.66666666
      CORN(24,4,2) = 0.0
!     
      CORN(24,1,3) = 0.25
      CORN(24,2,3) = 0.25
      CORN(24,3,3) = 0.25
      CORN(24,4,3) = 0.25
!     
      CORN(24,1,4) = 0.125
      CORN(24,2,4) = 0.125
      CORN(24,3,4) = 0.625
      CORN(24,4,4) = 0.125
!     
!     - define the Gaussian integration points in XIP space.
!     
      XIPGP(1)  = -1.0/SQRT(3.0)  
      ETAPGP(1) = -1.0/SQRT(3.0)  
!     
      XIPGP(2)  =  1.0/SQRT(3.0)  
      ETAPGP(2) = -1.0/SQRT(3.0)  
!     
      XIPGP(3)  = -1.0/SQRT(3.0)  
      ETAPGP(3) =  1.0/SQRT(3.0)  
!     
      XIPGP(4)  =  1.0/SQRT(3.0)  
      ETAPGP(4) =  1.0/SQRT(3.0)  
!     
!     - define positions of vertices of quadrilateral element
!     - associated with a subcell face in XIP
!     
      XIP(1)  = -1.0
      ETAP(1) = -1.0
!     
      XIP(2)  =  1.0
      ETAP(2) = -1.0
!     
      XIP(3)  = -1.0
      ETAP(3) =  1.0
!     
      XIP(4)  =  1.0
      ETAP(4) =  1.0
!     
!     - generate values of the FE basis functions at the quadrature
!     - points on the faces of the subcells.
!     
      do IFACE = 1, NFACE! Was loop 
!     
      do GJ = 1, FNGI! Was loop 
!     
            GI = (IFACE-1)*FNGI+GJ
!     
      do ICOORD = 1, NCOORD! Was loop 
!     
               POS(ICOORD)    = 0.0
               DPDXI(ICOORD)  = 0.0
               DPDETA(ICOORD) = 0.0
!     
      do JLOC = 1, FNLOC! Was loop 
!     
                  POS(ICOORD) = POS(ICOORD)&
     &                 + CORN(IFACE,ICOORD,JLOC)&
     &                 *0.25*(1.+XIP(JLOC)*XIPGP(GJ))&
     &                 *(1.+ETAP(JLOC)*ETAPGP(GJ))
!     
                  DPDXI(ICOORD) = DPDXI(ICOORD)&
     &                 + CORN(IFACE,ICOORD,JLOC)&
     &                 *0.25*XIP(JLOC)&
     &                 *(1.+ETAP(JLOC)*ETAPGP(GJ))
!     
                  DPDETA(ICOORD) = DPDETA(ICOORD)&
     &                 + CORN(IFACE,ICOORD,JLOC)&
     &                 *0.25*(1.+XIP(JLOC)*XIPGP(GJ))&
     &                 *ETAP(JLOC)
!     
               END DO
!     
            END DO
!     
            L1 = POS(1)
            L2 = POS(2)
            L3 = POS(3)
            L4 = POS(4)
!     
            SVN(1,GI)  = (2.0*L1-1.0)*L1
            SVN(2,GI)  =  4.0*L1*L3
            SVN(3,GI)  = (2.0*L3-1.0)*L3
            SVN(4,GI)  =  4.0*L2*L3
            SVN(5,GI)  = (2.0*L2-1.0)*L2
            SVN(6,GI)  =  4.0*L1*L2
            SVN(7,GI)  =  4.0*L1*L4
            SVN(8,GI)  =  4.0*L3*L4
            SVN(9,GI)  =  4.0*L2*L4
            SVN(10,GI) = (2.0*L4-1.0)*L4
!     
            SVNLX(1,GI)  = (4.0*L1-1.0)*DPDXI(1)
            SVNLX(2,GI)  =  4.0*DPDXI(1)*L3 + 4.0*L1*DPDXI(3)
            SVNLX(3,GI)  = (4.0*L3-1.0)*DPDXI(3)
            SVNLX(4,GI)  =  4.0*DPDXI(2)*L3 + 4.0*L2*DPDXI(3)
            SVNLX(5,GI)  = (4.0*L2-1.0)*DPDXI(2)
            SVNLX(6,GI)  =  4.0*DPDXI(1)*L2 + 4.0*L1*DPDXI(2)
            SVNLX(7,GI)  =  4.0*DPDXI(1)*L4 + 4.0*L1*DPDXI(4)
            SVNLX(8,GI)  =  4.0*DPDXI(3)*L4 + 4.0*L3*DPDXI(4)
            SVNLX(9,GI)  =  4.0*DPDXI(2)*L4 + 4.0*L2*DPDXI(4)
            SVNLX(10,GI) = (4.0*L4-1.0)*DPDXI(4)
!     
            SVNLY(1,GI)  = (4.0*L1-1.0)*DPDETA(1)
            SVNLY(2,GI)  =  4.0*DPDETA(1)*L3 + 4.0*L1*DPDETA(3)
            SVNLY(3,GI)  = (4.0*L3-1.0)*DPDETA(3)
            SVNLY(4,GI)  =  4.0*DPDETA(2)*L3 + 4.0*L2*DPDETA(3)
            SVNLY(5,GI)  = (4.0*L2-1.0)*DPDETA(2)
            SVNLY(6,GI)  =  4.0*DPDETA(1)*L2 + 4.0*L1*DPDETA(2)
            SVNLY(7,GI)  =  4.0*DPDETA(1)*L4 + 4.0*L1*DPDETA(4)
            SVNLY(8,GI)  =  4.0*DPDETA(3)*L4 + 4.0*L3*DPDETA(4)
            SVNLY(9,GI)  =  4.0*DPDETA(2)*L4 + 4.0*L2*DPDETA(4)
            SVNLY(10,GI) = (4.0*L4-1.0)*DPDETA(4)
!     
         END DO
!     
      END DO
!     
!     - set the weights for the surface integration
!     
      do GI = 1, SVNGI! Was loop 
!     
         SVWEIGH(GI) = 1.0
!     
      END DO
!     
      END SUBROUTINE

      SUBROUTINE SFELIN( SNGI, SNLOC, SN, SNLX, SWEIGH )
!     --------------------------------------------------
!
! - this subroutine generates the FE surface
! - shape functions and the associated derivatives
! - and weights.
!
!   -------------------------------
! - date last modified : 14/06/2003
!   -------------------------------
!
      IMPLICIT NONE
!
      INTEGER SNGI, SNLOC
!
      REAL SN(SNLOC,SNGI), SNLX(SNLOC,SNGI)
!
      REAL SWEIGH(SNGI)
!
! - local variables
!
      INTEGER SILOC, SJLOC
      INTEGER SGI,   SGJ
!
      REAL XIGP(1), XI(2)
!
! - positions of Gauss points associated with CV's in (ETAP) space
! 
      XIGP(1) = 0.0
!
! - positions of vertices of line element in (ETA) space
!
      XI(1) = -1.0
      XI(2) =  1.0
!
! - values of FE basis functions at volume integration points
!
      do SILOC = 1, SNLOC! Was loop 
!
      do SJLOC = 1, SNLOC! Was loop 
!
      do SGJ = 1, 1! Was loop 
!
               SGI = (SJLOC-1)*(SNLOC-1)+SGJ
!
               SN(SILOC,SGI) = 0.5*(1.+XI(SILOC)&
     &                        *0.5*(XI(SJLOC)+XIGP(SGJ)))
!
            END DO
!
         END DO
!
      END DO
!
! - derivatives of FE basis functions at the Gauss points
!
      do SILOC = 1, SNLOC! Was loop 
!
      do SGI = 1, SNGI! Was loop 
!
            SNLX(SILOC,SGI) = 0.5*XI(SILOC)
!
         END DO
!
      END DO
!
! - note that the weight of the normal Gauss point is 2.0
! - but this is reduced by 0.5 i.e. 1/2 because of the 
! - second Jacobian
!
      do SGI = 1, SNGI! Was loop 
!
         SWEIGH(SGI) = 1.0
!
      END DO
!
      END SUBROUTINE

      SUBROUTINE SFEQUAD( SNGI, SNLOC, &
! - REALS
     &                    SN,   SNLX, &
     &                    SNLY, SWEIGH )
!     ---------------------------------
!
! - this subroutine generates the FE volume shape 
! - functions and the associated derivatives and
! - weights.
!
!   -------------------------------
! - date last modified : 12/10/2002
!   -------------------------------
!
      IMPLICIT NONE
!
      INTEGER SNGI, SNLOC
!
      REAL SN(SNLOC,SNGI),   SNLX(SNLOC,SNGI)
      REAL SNLY(SNLOC,SNGI), SWEIGH(SNGI)
!
! - local variables
!
      INTEGER SILOC, SJLOC, SKLOC, ICOORD
!
      INTEGER SGI, SGJ
!
! - note that NCOORD is the number of co-ordinates,
! - SBNLOC is the number of local nodes for a subcell
! - and SBNGI is the number of local integration points.
!
      INTEGER NCOORD, SBNLOC, SBNGI
      PARAMETER( NCOORD = 2, SBNLOC = 4, SBNGI = 1 )
!
      REAL POS(NCOORD), DPDXI(NCOORD), DPDETA(NCOORD)
!
      REAL XIGP(SBNGI), ETAGP(SBNGI)
!
      REAL XI(SBNLOC), ETA(SBNLOC)
!
      REAL CORN(SNLOC,NCOORD,SBNLOC)
!
! - positions of the corners (vertices) of the subcells in area
! - co-ordinates. Note that CORN(SJLOC,ICOORD,SKLOC). Note that
! - SJLOC ranges from 1-4, and ICOORD ranges from 1-2 and SKLOC
! - ranges from 1-4. SJLOC signifies the particular subcell
! - ICOORD signifies the co-ordinate of the subcell (XI,ETA)
! - and SKLOC signifies the vertex of the subcell (which are quad)
!
!
! - subcell 1
!
      CORN(1,1,1) = -1.0
      CORN(1,2,1) = -1.0
!
      CORN(1,1,2) =  0.0
      CORN(1,2,2) = -1.0
!
      CORN(1,1,3) = -1.0
      CORN(1,2,3) =  0.0
!
      CORN(1,1,4) =  0.0
      CORN(1,2,4) =  0.0
!
! - subcell 2
!      
      CORN(2,1,1) =  0.0
      CORN(2,2,1) = -1.0
!
      CORN(2,1,2) =  1.0
      CORN(2,2,2) = -1.0
!
      CORN(2,1,3) =  0.0
      CORN(2,2,3) =  0.0
!
      CORN(2,1,4) =  1.0
      CORN(2,2,4) =  0.0
!
! - subcell 3
!      
      CORN(3,1,1) = -1.0
      CORN(3,2,1) =  0.0
!
      CORN(3,1,2) =  0.0
      CORN(3,2,2) =  0.0
!
      CORN(3,1,3) = -1.0
      CORN(3,2,3) =  1.0
!
      CORN(3,1,4) =  0.0
      CORN(3,2,4) =  1.0
!
! - subcell 4
!      
      CORN(4,1,1) =  0.0
      CORN(4,2,1) =  0.0
!
      CORN(4,1,2) =  1.0
      CORN(4,2,2) =  0.0
!
      CORN(4,1,3) =  0.0
      CORN(4,2,3) =  1.0
!
      CORN(4,1,4) =  1.0
      CORN(4,2,4) =  1.0
!
! - define Gaussian integration points in (XIP,ETAP) space.
!
      XIGP(1)   = 0.0  
      ETAGP(1)  = 0.0
!
! - define positions of vertices of quadrilateral element in (XIP,ETAP) space.
!
      XI(1)   = -1.0
      ETA(1)  = -1.0
!
      XI(2)   =  1.0
      ETA(2)  = -1.0
!
      XI(3)   = -1.0
      ETA(3)  =  1.0
!
      XI(4)   =  1.0
      ETA(4)  =  1.0
!
! - generate values of the FE basis functions at Gaussian integration points
!
      do SILOC = 1, SNLOC! Was loop 
!
        do SJLOC = 1, SNLOC! Was loop 
!
            do SGJ = 1, SBNGI! Was loop 
!
               SGI = SJLOC-1 + SGJ
!
               do ICOORD = 1, NCOORD! Was loop 
!
                  POS(ICOORD)     = 0.0
                  DPDXI(ICOORD)   = 0.0
                  DPDETA(ICOORD)  = 0.0
!
                  do SKLOC = 1, SBNLOC! Was loop 
!
                     POS(ICOORD) = POS(ICOORD) &
     &                           + CORN(SJLOC,ICOORD,SKLOC)&
     &                             *0.25*(1.+XI(SKLOC)*XIGP(SGJ))&
     &                             *(1.+ETA(SKLOC)*ETAGP(SGJ))
!
                     DPDXI(ICOORD) = DPDXI(ICOORD) &
     &                             + CORN(SJLOC,ICOORD,SKLOC)&
     &                               *0.25*XI(SKLOC)&
     &                               *(1.+ETA(SKLOC)*ETAGP(SGJ))
!
                     DPDETA(ICOORD) = DPDETA(ICOORD) &
     &                              + CORN(SJLOC,ICOORD,SKLOC)&
     &                                *0.25*(1.+XI(SKLOC)*XIGP(SGJ))&
     &                                *ETA(SKLOC)
!
                  END DO  ! skloc
!
                END DO  ! icoord
!
            SN(SILOC,SGI) = 0.25*( 1.0 + XI(SILOC)*POS(1) )&
     &                        *( 1.0 + ETA(SILOC)*POS(2) )
!
            SNLX(SILOC,SGI) = 0.25*XI(SILOC)*DPDXI(1)&
     &                         *( 1.0 + ETA(SILOC)*POS(2) )&
     &                   + 0.25*( 1.0 + XI(SILOC)*POS(1) )&
     &                         *ETA(SILOC)*DPDXI(2)
!
            SNLY(SILOC,SGI) = 0.25*XI(SILOC)*DPDETA(1)&
     &                         *( 1.0 + ETA(SILOC)*POS(2) )&
     &                   + 0.25*( 1.0 + XI(SILOC)*POS(1) )&
     &                         *ETA(SILOC)*DPDETA(2)
!
            END DO  ! sgj
!
         END DO   ! sjloc
!
      END DO  ! siloc
!
      do SGI = 1, SNGI! Was loop 
!
         SWEIGH(SGI) = 4.0
!
      END DO
!
      END SUBROUTINE

      SUBROUTINE SFETRI( SNGI, SNLOC, &
! - REALS
     &                   SN,   SNLX, &
     &                   SNLY, SWEIGH )
!     ---------------------------------
!
! - this subroutine generates the FE shape functions 
! - and the associated derivatives and weights for 
! - surface elements.
! 
!   ------------------------------- 
! - date last modified : 12/10/2002
!   -------------------------------
!
      IMPLICIT NONE
!
      INTEGER SNGI, SNLOC
!
      REAL SN(SNLOC,SNGI),   SNLX(SNLOC,SNGI)
      REAL SNLY(SNLOC,SNGI), SWEIGH(SNGI) 
!
! - local variables
!
      INTEGER SILOC, SJLOC, SKLOC, ICOORD
!
      INTEGER SGI, SGJ
!
! - note that NCOORD is the number of co-ordinates,
! - SBNLOC is the number of local nodes for a subcell
! - and SBNGI is the number of local integration points.
!
      INTEGER NCOORD, SBNLOC, SBNGI
      PARAMETER( NCOORD = 3, SBNLOC = 4, SBNGI = 1 )
!
      REAL POS(NCOORD), DPDXI(NCOORD), DPDETA(NCOORD)
!
      REAL XIGP(SBNGI), ETAGP(SBNGI) 
!
      REAL XI(SBNLOC),   ETA(SBNLOC)
!
      REAL CORN(SNLOC,NCOORD,SBNLOC)
!
! - positions of the corners (vertices) of the subcells in area
! - co-ordinates. Note that CORN(SJLOC,ICOORD,SKLOC). Note that
! - SJLOC ranges from 1-3, and ICOORD ranges from 1-3 and SKLOC
! - ranges from 1-4. SJLOC signifies the particular subcell
! - SJLOC signifies the co-ordinate of the subcell (L1,L2,L3) 
! - and SKLOC signifies the vertex of the subcell (which are quads)
!
!
! - subcell 1
!
      CORN(1,1,1) = 1.0
      CORN(1,2,1) = 0.0
      CORN(1,3,1) = 0.0
!
      CORN(1,1,2) = 0.5
      CORN(1,2,2) = 0.5
      CORN(1,3,2) = 0.0
!
      CORN(1,1,3) = 0.5
      CORN(1,2,3) = 0.0
      CORN(1,3,3) = 0.5
!
      CORN(1,1,4) = 0.33333333
      CORN(1,2,4) = 0.33333333
      CORN(1,3,4) = 0.33333333
!
! - subcell 2
!
      CORN(2,1,1) = 0.5
      CORN(2,2,1) = 0.5
      CORN(2,3,1) = 0.0
!
      CORN(2,1,2) = 0.0
      CORN(2,2,2) = 1.0
      CORN(2,3,2) = 0.0
!
      CORN(2,1,3) = 0.33333333
      CORN(2,2,3) = 0.33333333
      CORN(2,3,3) = 0.33333333
!
      CORN(2,1,4) = 0.0
      CORN(2,2,4) = 0.5
      CORN(2,3,4) = 0.5
!
! - subcell 3
!
      CORN(3,1,1) = 0.5
      CORN(3,2,1) = 0.0
      CORN(3,3,1) = 0.5
!
      CORN(3,1,2) = 0.33333333
      CORN(3,2,2) = 0.33333333
      CORN(3,3,2) = 0.33333333
!
      CORN(3,1,3) = 0.0
      CORN(3,2,3) = 0.0
      CORN(3,3,3) = 1.0
!
      CORN(3,1,4) = 0.0
      CORN(3,2,4) = 0.5
      CORN(3,3,4) = 0.5
!
! - define Gaussian integration points in (XI,ETA) space.
!
      XIGP(1)  =  0.0
      ETAGP(1) =  0.0
!
! - define postions of quadrilateral element in (XI,ETA) space.
!
      XI(1)  = -1.0
      ETA(1) = -1.0
!
      XI(2)  =  1.0
      ETA(2) = -1.0
!
      XI(3)  = -1.0
      ETA(3) =  1.0
!
      XI(4)  =  1.0
      ETA(4) =  1.0
!
! - generate values of FE basis functions at Gaussian integration points
!
      do SILOC = 1, SNLOC! Was loop 
!
      do SJLOC = 1, SNLOC! Was loop 
!
      do SGJ = 1, SBNGI! Was loop 
!
              SGI = SJLOC-1 + SGJ
!
      do ICOORD = 1, NCOORD! Was loop 
!
                 POS(ICOORD)     = 0.0
                 DPDXI(ICOORD)   = 0.0
                 DPDETA(ICOORD)  = 0.0
!
      do SKLOC = 1, SBNLOC! Was loop 
!
                    POS(ICOORD) = POS(ICOORD) &
     &                          + CORN(SJLOC,ICOORD,SKLOC)*&
     &                            0.25*(1.+XI(SKLOC)*XIGP(SGJ))&
     &                            *(1.+ETA(SKLOC)*ETAGP(SGJ))
!
!
                    DPDXI(ICOORD) = DPDXI(ICOORD) &
     &                            + CORN(SJLOC,ICOORD,SKLOC)&
     &                              *0.25*XI(SKLOC)&
     &                              *(1.+ETA(SKLOC)*ETAGP(SGJ))
!
                    DPDETA(ICOORD) = DPDETA(ICOORD) &
     &                             + CORN(SJLOC,ICOORD,SKLOC)&
     &                               *0.25*(1.+XI(SKLOC)*XIGP(SGJ))&
     &                               *ETA(SKLOC)
!
                  END DO
!
               END DO
!
               SN(SILOC,SGI)   = POS(SILOC)
!
               SNLX(SILOC,SGI) = DPDXI(SILOC)
!
               SNLY(SILOC,SGI) = DPDETA(SILOC)
!
            END DO
!
         END DO
!
      END DO
!
!      do siloc = 1, snloc
!
!         write( *,* ) ( sn(siloc,sgi), sgi = 1, sngi )
!
!      end do
! 
!      stop 11
!
      do SGI = 1, SNGI! Was loop 
!
         SWEIGH(SGI) = 4.0
!
      END DO
!
      END SUBROUTINE

      SUBROUTINE SFEQLIN( SNGI, SNLOC, SN, SNLX, SWEIGH )
!     ---------------------------------------------------
!
! - this subroutine generates the FE surface
! - shape functions and the associated derivatives
! - and weights.
!
!   -------------------------------
! - date last modified : 14/06/2003
!   -------------------------------
!
      IMPLICIT NONE
!
      INTEGER SNGI, SNLOC
!
      REAL SN(SNLOC,SNGI), SNLX(SNLOC,SNGI)
!
      REAL SWEIGH(SNGI)
!
! - local variables
!
      INTEGER SIFACE, SILOC, SJLOC
      INTEGER SGI,    SGJ
!
! - note that SNFACE is the number of surfaces internal to
! - a surface element.
!
      INTEGER SNFACE, SFNGI, SFNLOC
!
      PARAMETER( SNFACE = 3, SFNGI = 2, SFNLOC = 2 ) 
!
      REAL POS, DPDXI, XI
      REAL DLIJXIDXI(3), LIJXI(3)
      REAL XIP(SFNLOC),  XIPGP(SFNGI) 
      REAL CORN(SNFACE,SFNLOC)
!
! - face 1
!
      CORN(1,1) = -1.0 
      CORN(1,2) = -0.5
!
! - face 2
!
      CORN(2,1) = -0.5
      CORN(2,2) =  0.5
!
! - face 3
!
      CORN(3,1) =  0.5
      CORN(3,2) =  1.0
!
! - positions of Gauss points associated with SCVs in ETAP space
!
      XIPGP(1) = -1/SQRT(3.0)
      XIPGP(2) =  1/SQRT(3.0)
!
! - define positions of line element assciated with a SCV
! - in XIP space.
!
      XIP(1) = -1.0
      XIP(2) =  1.0
!
! - generate values of the FE basis functions at the quadrature 
! - points on the faces of the subcells.
!
      do SIFACE = 1, SNFACE! Was loop 
!
      do SGJ = 1, SFNGI! Was loop 
!
               SGI = (SIFACE-1)*SFNGI+SGJ
!
               POS   = 0.0
               DPDXI = 0.0
!
      do SJLOC = 1, SFNLOC! Was loop 
!
                  POS   = POS &
     &                  + CORN(SIFACE,SJLOC)&
     &                    *0.5*(1.+XIP(SJLOC)*XIPGP(SGJ))
!
                  DPDXI = DPDXI&
     &                  + CORN(SIFACE,SJLOC)*0.5*XIP(SJLOC)
!
               END DO
!
               XI = POS
!
               LIJXI(1)   =  0.5*XI*(XI-1.0)
               LIJXI(2)   =  1.0-XI*XI
               LIJXI(3)   =  0.5*XI*(XI+1.0)
!
               DLIJXIDXI(1) =  0.5*(2.0*XI-1.0)*DPDXI
               DLIJXIDXI(2) = -2.0*XI*DPDXI
               DLIJXIDXI(3) =  0.5*(2.0*XI+1.0)*DPDXI
!
      do SILOC = 1, SNLOC! Was loop 
!
                  SN(SILOC,SGI)   = LIJXI(SILOC)
!
                  SNLX(SILOC,SGI) = DLIJXIDXI(SILOC)
!  
               END DO
!
            END DO
!
         END DO
!
! - set the weights for the surface integration
!
      do SGI = 1, SNGI! Was loop 
!
        SWEIGH(SGI) = 1.0
!
      END DO
!
      END SUBROUTINE

      SUBROUTINE SFEQQUAD( SNGI, SNLOC, &
! - REALS
     &                     SN,   SNLX, &
     &                     SNLY, SWEIGH )
!     -----------------------------------
!
! - this subroutine generates the FE surface shape 
! - functions and the associated derivatives and
! - weights.
!
!   -------------------------------
! - date last modified : 12/10/2002
!   -------------------------------
!
      IMPLICIT NONE
!
      INTEGER SNGI, SNLOC
!
      REAL SN(SNLOC,SNGI),   SNLX(SNLOC,SNGI)
      REAL SNLY(SNLOC,SNGI), SWEIGH(SNGI)
!
! - local variables
!
      INTEGER SILOC, SJLOC, SKLOC, ICOORD
!
      INTEGER SGI, SGJ
!
! - note that NCOORD is the number of co-ordinates,
! - SBNLOC is the number of local nodes for a subcell
! - and SBNGI is the number of local integration points.
!
      INTEGER NCOORD, NSCV, SBNLOC, SBNGI
      PARAMETER( NCOORD = 2, NSCV = 9 )
      PARAMETER( SBNLOC = 4, SBNGI = 4 )
!
      REAL POS(NCOORD), DPDXI(NCOORD), DPDETA(NCOORD)
!
      REAL XIGP(SBNGI), ETAGP(SBNGI)
!
      REAL XI, ETA
!
      REAL XIP(SBNLOC), ETAP(SBNLOC)
!
      REAL CORN(NSCV,NCOORD,SBNLOC)
!
! - for derivatives wrt the second parametric space.
!
      REAL LIJXI(3),      LIJETA(3)
      REAL DLIJXIDXI(3),  DLIJXIDETA(3)
      REAL DLIJETADXI(3), DLIJETADETA(3)  
!
      INTEGER I, J
!
! - positions of the corners (vertices) of the subcells in area
! - co-ordinates. Note that CORN(JLOC,ICOORD,KLOC). Note that
! - JLOC ranges from 1-9, and ICOORD ranges from 1-2 and KLOC
! - ranges from 1-4. JLOC signifies the particular subcell
! - ICOORD signifies the co-ordinate of the subcell (XI,ETA)
! - and KLOC signifies the vertex of the subcell (which are quad)
!
!
! - subcell 1
!
      CORN(1,1,1) = -1.0
      CORN(1,2,1) = -1.0
!
      CORN(1,1,2) = -0.5
      CORN(1,2,2) = -1.0
!
      CORN(1,1,3) = -1.0
      CORN(1,2,3) = -0.5
!
      CORN(1,1,4) = -0.5
      CORN(1,2,4) = -0.5
!
! - subcell 2
!
      CORN(2,1,1) = -0.5
      CORN(2,2,1) = -1.0
!
      CORN(2,1,2) =  0.5
      CORN(2,2,2) = -1.0
!
      CORN(2,1,3) = -0.5
      CORN(2,2,3) = -0.5
!
      CORN(2,1,4) =  0.5
      CORN(2,2,4) = -0.5
!
! - subcell 3
! 
      CORN(3,1,1) =  0.5
      CORN(3,2,1) = -1.0
!
      CORN(3,1,2) =  1.0
      CORN(3,2,2) = -1.0
!
      CORN(3,1,3) =  0.5
      CORN(3,2,3) = -0.5
!
      CORN(3,1,4) =  1.0
      CORN(3,2,4) = -0.5
!
! - subcell 4
! 
      CORN(4,1,1) = -1.0
      CORN(4,2,1) = -0.5
!
      CORN(4,1,2) = -0.5
      CORN(4,2,2) = -0.5
!
      CORN(4,1,3) = -1.0
      CORN(4,2,3) =  0.5
!
      CORN(4,1,4) = -0.5
      CORN(4,2,4) =  0.5
!
! - subcell 5
! 
      CORN(5,1,1) = -0.5
      CORN(5,2,1) = -0.5
!
      CORN(5,1,2) =  0.5
      CORN(5,2,2) = -0.5
!
      CORN(5,1,3) = -0.5
      CORN(5,2,3) =  0.5
!
      CORN(5,1,4) =  0.5
      CORN(5,2,4) =  0.5
!
! - subcell 6
! 
      CORN(6,1,1) =  0.5
      CORN(6,2,1) = -0.5
!
      CORN(6,1,2) =  1.0
      CORN(6,2,2) = -0.5
!
      CORN(6,1,3) =  0.5
      CORN(6,2,3) =  0.5
!
      CORN(6,1,4) =  1.0
      CORN(6,2,4) =  0.5
!
! - subcell 7
! 
      CORN(7,1,1) = -1.0
      CORN(7,2,1) =  0.5
!
      CORN(7,1,2) = -0.5
      CORN(7,2,2) =  0.5
!
      CORN(7,1,3) = -1.0
      CORN(7,2,3) =  1.0
!
      CORN(7,1,4) = -0.5
      CORN(7,2,4) =  1.0
!
! - subcell 8
! 
      CORN(8,1,1) = -0.5
      CORN(8,2,1) =  0.5
!
      CORN(8,1,2) =  0.5
      CORN(8,2,2) =  0.5
!
      CORN(8,1,3) = -0.5
      CORN(8,2,3) =  1.0
!
      CORN(8,1,4) =  0.5
      CORN(8,2,4) =  1.0
!
! - subcell 9
! 
      CORN(9,1,1) =  0.5
      CORN(9,2,1) =  0.5
!
      CORN(9,1,2) =  1.0
      CORN(9,2,2) =  0.5
!
      CORN(9,1,3) =  0.5
      CORN(9,2,3) =  1.0
!
      CORN(9,1,4) =  1.0
      CORN(9,2,4) =  1.0
!
! - define Gaussian integration points in (XIP,ETAP) space.
!
      XIGP(1)   = -1/SQRT(3.0)  
      ETAGP(1)  = -1/SQRT(3.0)  
!
      XIGP(2)   =  1/SQRT(3.0) 
      ETAGP(2)  = -1/SQRT(3.0) 
!
      XIGP(3)   = -1/SQRT(3.0)   
      ETAGP(3)  =  1/SQRT(3.0)   
!
      XIGP(4)   =  1/SQRT(3.0)     
      ETAGP(4)  =  1/SQRT(3.0)   
!
! - define positions of vertices of quadrilateral element in (XIP,ETAP) space.
!
      XIP(1)   = -1.0
      ETAP(1)  = -1.0
!
      XIP(2)   =  1.0
      ETAP(2)  = -1.0
!
      XIP(3)   = -1.0
      ETAP(3)  =  1.0
!
      XIP(4)   =  1.0
      ETAP(4)  =  1.0
!
! - generate values of the FE basis functions at Gaussian integration points
!
      do SJLOC = 1, NSCV! Was loop 
!
      do SGJ = 1, SBNGI! Was loop 
!
            SGI = (SJLOC-1)*SBNLOC+SGJ
!
      do ICOORD = 1, NCOORD! Was loop 
!
               POS(ICOORD)     = 0.0
               DPDXI(ICOORD)   = 0.0
               DPDETA(ICOORD)  = 0.0
!
      do SKLOC = 1, SBNLOC! Was loop 
!
                  POS(ICOORD) = POS(ICOORD) &
     &                          + CORN(SJLOC,ICOORD,SKLOC)&
     &                          *0.25*(1.+XIP(SKLOC)*XIGP(SGJ))&
     &                           *(1.+ETAP(SKLOC)*ETAGP(SGJ))
!
                  DPDXI(ICOORD) = DPDXI(ICOORD) &
     &                            + CORN(SJLOC,ICOORD,SKLOC)&
     &                            *0.25*XIP(SKLOC)&
     &                            *(1.+ETAP(SKLOC)*ETAGP(SGJ))
!
                  DPDETA(ICOORD) = DPDETA(ICOORD) &
     &                             + CORN(SJLOC,ICOORD,SKLOC)&
     &                             *0.25*(1.+XIP(SKLOC)*XIGP(SGJ))&
     &                             *ETAP(SKLOC)
!
               END DO
!
            END DO
!
            XI  = POS(1)
            ETA = POS(2)
!
            LIJXI(1)   =  0.5*XI*(XI-1.0)
            LIJXI(2)   =  1.0-XI*XI
            LIJXI(3)   =  0.5*XI*(XI+1.0)
!
            LIJETA(1)  =  0.5*ETA*(ETA-1.0)
            LIJETA(2)  =  1.0-ETA*ETA
            LIJETA(3)  =  0.5*ETA*(ETA+1.0)
!
            DLIJXIDXI(1) =  0.5*(2.0*XI-1.0)*DPDXI(1)
            DLIJXIDXI(2) = -2.0*XI*DPDXI(1)
            DLIJXIDXI(3) =  0.5*(2.0*XI+1.0)*DPDXI(1)
!
            DLIJXIDETA(1) =  0.5*(2.0*XI-1.0)*DPDETA(1)
            DLIJXIDETA(2) = -2.0*XI*DPDETA(1)
            DLIJXIDETA(3) =  0.5*(2.0*XI+1.0)*DPDETA(1)
!
            DLIJETADXI(1) =  0.5*(2.0*ETA-1.0)*DPDXI(2)
            DLIJETADXI(2) = -2.0*ETA*DPDXI(2)
            DLIJETADXI(3) =  0.5*(2.0*ETA+1.0)*DPDXI(2)
!           
            DLIJETADETA(1) =  0.5*(2.0*ETA-1.0)*DPDETA(2)
            DLIJETADETA(2) = -2.0*ETA*DPDETA(2)
            DLIJETADETA(3) =  0.5*(2.0*ETA+1.0)*DPDETA(2)
!
      do I = 1, 3! Was loop 
!
      do J = 1, 3! Was loop 
!
                  SILOC = I+(J-1)*3
!
                  SN(SILOC,SGI)   = LIJXI(I)*LIJETA(J)
!
                  SNLX(SILOC,SGI) = LIJXI(I)*DLIJETADXI(J)&
     &                            + DLIJXIDXI(I)*LIJETA(J)
!
                  SNLY(SILOC,SGI) = LIJXI(I)*DLIJETADETA(J)&
     &                            + DLIJXIDETA(I)*LIJETA(J)
!
               END DO
!
            END DO               
!
         END DO
!
      END DO
!
      do SGI = 1, SNGI! Was loop 
!
         SWEIGH(SGI) = 1.0
!
      END DO
!
      END SUBROUTINE

      SUBROUTINE SFEQTRI( SNGI,   SNLOC, &
! - REALS
     &                    SN,     SNLX, &
     &                    SNLY,   SWEIGH )
!     ------------------------------------
!
! - this subroutine generates the FE surface shape 
! - functions and the associated derivatives and
! - weights.
! 
!   ------------------------------- 
! - date last modified : 12/10/2002
!   -------------------------------
!
      IMPLICIT NONE
!
      INTEGER SNGI, SNLOC
!
      REAL SN(SNLOC,SNGI),   SNLX(SNLOC,SNGI)
      REAL SNLY(SNLOC,SNGI), SWEIGH(SNGI) 
!
! - local variables
!
      INTEGER SILOC, SJLOC, SKLOC, ICOORD
!
      INTEGER SGI, SGJ
!
! - note that NCOORD is the number of co-ordinates,
! - SBNLOC is the number of local nodes for a subcell
! - and SBNGI is the number of local integration points.
!
      INTEGER NCOORD, NSCV, SBNLOC, SBNGI
      PARAMETER( NCOORD = 3, NSCV  = 9 )
      PARAMETER( SBNLOC = 4, SBNGI = 4 )
!
      REAL POS(NCOORD), DPDXI(NCOORD), DPDETA(NCOORD)
!
      REAL XIPGP(SBNGI), ETAPGP(SBNGI) 
!
      REAL XIP(SBNLOC),   ETAP(SBNLOC)
!
      REAL CORN(NSCV,NCOORD,SBNLOC)
!
      REAL L1, L2, L3
!
! - positions of the corners (vertices) of the subcells in area
! - co-ordinates. Note that CORN(JLOC,ICOORD,KLOC). Note that
! - JLOC ranges from 1-9, and ICOORD ranges from 1-3 and KLOC
! - ranges from 1-4. JLOC signifies the particular subcell
! - JLOC signifies the co-ordinate of the subcell (L1,L2,L3) 
! - and KLOC signifies the vertex of the subcell (which are quads)
!
! - Note that for quadratic triangles the mid-side nodes
! - produce pentagonal shaped CVs these are broken into 
! - two quadrilateral shaped CVs. Thus there are a total
! - of 9 quadrilateral shaped SCVs. But really there are
! - only the six SCVs with three SCVs split into two. 
!
! - subcell 1
!
      CORN(1,1,1) = 1.0
      CORN(1,2,1) = 0.0
      CORN(1,3,1) = 0.0
!
      CORN(1,1,2) = 0.75
      CORN(1,2,2) = 0.25
      CORN(1,3,2) = 0.0
!
      CORN(1,1,3) = 0.75
      CORN(1,2,3) = 0.0
      CORN(1,3,3) = 0.25
!                                 
      CORN(1,1,4) = 0.66666666
      CORN(1,2,4) = 0.16666666
      CORN(1,3,4) = 0.16666666
!
! - subcell 2
!
      CORN(2,1,1) = 0.75
      CORN(2,2,1) = 0.25
      CORN(2,3,1) = 0.0
!
      CORN(2,1,2) = 0.5
      CORN(2,2,2) = 0.5
      CORN(2,3,2) = 0.0
!
      CORN(2,1,3) = 0.66666666
      CORN(2,2,3) = 0.16666666
      CORN(2,3,3) = 0.16666666
!                                 
      CORN(2,1,4) = 0.33333333
      CORN(2,2,4) = 0.33333333
      CORN(2,3,4) = 0.33333333
!
! - subcell 3
!     
      CORN(3,1,1) = 0.5
      CORN(3,2,1) = 0.5
      CORN(3,3,1) = 0.0
!
      CORN(3,1,2) = 0.25
      CORN(3,2,2) = 0.75
      CORN(3,3,2) = 0.0
!
      CORN(3,1,3) = 0.33333333
      CORN(3,2,3) = 0.33333333
      CORN(3,3,3) = 0.33333333
!                                 
      CORN(3,1,4) = 0.16666666
      CORN(3,2,4) = 0.66666666
      CORN(3,3,4) = 0.16666666
!
! - subcell 4
! 
      CORN(4,1,1) = 0.25
      CORN(4,2,1) = 0.75
      CORN(4,3,1) = 0.0
!
      CORN(4,1,2) = 0.0
      CORN(4,2,2) = 1.0
      CORN(4,3,2) = 0.0
!
      CORN(4,1,3) = 0.16666666
      CORN(4,2,3) = 0.66666666
      CORN(4,3,3) = 0.16666666
!                                 
      CORN(4,1,4) = 0.0
      CORN(4,2,4) = 0.75
      CORN(4,3,4) = 0.25
!
! - subcell 5
!
      CORN(5,1,1) = 0.16666666
      CORN(5,2,1) = 0.66666666
      CORN(5,3,1) = 0.16666666
!
      CORN(5,1,2) = 0.0
      CORN(5,2,2) = 0.75
      CORN(5,3,2) = 0.25
!
      CORN(5,1,3) = 0.33333333
      CORN(5,2,3) = 0.33333333
      CORN(5,3,3) = 0.33333333
!                                 
      CORN(5,1,4) = 0.0
      CORN(5,2,4) = 0.5
      CORN(5,3,4) = 0.5
!
! - subcell 6
!
      CORN(6,1,1) = 0.33333333
      CORN(6,2,1) = 0.33333333
      CORN(6,3,1) = 0.33333333
!
      CORN(6,1,2) = 0.0
      CORN(6,2,2) = 0.5
      CORN(6,3,2) = 0.5
!
      CORN(6,1,3) = 0.16666666
      CORN(6,2,3) = 0.16666666
      CORN(6,3,3) = 0.66666666
!                                 
      CORN(6,1,4) = 0.0
      CORN(6,2,4) = 0.25
      CORN(6,3,4) = 0.75
!
! - subcell 7
!
      CORN(7,1,1) = 0.16666666
      CORN(7,2,1) = 0.16666666
      CORN(7,3,1) = 0.66666666
!
      CORN(7,1,2) = 0.0
      CORN(7,2,2) = 0.25
      CORN(7,3,2) = 0.75
!
      CORN(7,1,3) = 0.25
      CORN(7,2,3) = 0.0
      CORN(7,3,3) = 0.75
!                                 
      CORN(7,1,4) = 0.0
      CORN(7,2,4) = 0.0
      CORN(7,3,4) = 1.0
!
! - subcell 8
!      
      CORN(8,1,1) = 0.5
      CORN(8,2,1) = 0.0
      CORN(8,3,1) = 0.5
!
      CORN(8,1,2) = 0.33333333
      CORN(8,2,2) = 0.33333333
      CORN(8,3,2) = 0.33333333
!
      CORN(8,1,3) = 0.25
      CORN(8,2,3) = 0.0
      CORN(8,3,3) = 0.75
!                                 
      CORN(8,1,4) = 0.16666666
      CORN(8,2,4) = 0.16666666
      CORN(8,3,4) = 0.66666666
!
! - subcell 9
!  
      CORN(9,1,1) = 0.75
      CORN(9,2,1) = 0.0
      CORN(9,3,1) = 0.25
!
      CORN(9,1,2) = 0.66666666
      CORN(9,2,2) = 0.16666666
      CORN(9,3,2) = 0.16666666
!
      CORN(9,1,3) = 0.5
      CORN(9,2,3) = 0.0
      CORN(9,3,3) = 0.5
!                                 
      CORN(9,1,4) = 0.33333333
      CORN(9,2,4) = 0.33333333
      CORN(9,3,4) = 0.33333333
!
! - define Gaussian integration points in (XIP,ETAP) space.
!
      XIPGP(1)  = -1./SQRT(3.0)
      ETAPGP(1) = -1./SQRT(3.0)
!
      XIPGP(2)  =  1./SQRT(3.0)
      ETAPGP(2) = -1./SQRT(3.0)
!
      XIPGP(3)  = -1./SQRT(3.0)
      ETAPGP(3) =  1./SQRT(3.0)
!
      XIPGP(4)  =  1./SQRT(3.0)
      ETAPGP(4) =  1./SQRT(3.0)
!
! - define postions of quadrilateral element in (XIP,ETAP) space.
!
      XIP(1)  = -1.0
      ETAP(1) = -1.0
!
      XIP(2)  =  1.0
      ETAP(2) = -1.0
!
      XIP(3)  = -1.0
      ETAP(3) =  1.0
!
      XIP(4)  =  1.0
      ETAP(4) =  1.0
!
! - generate values of FE basis function at Gaussian integration
! - points.
!
      do SJLOC = 1, NSCV! Was loop 
!
      do SGJ = 1, SBNGI! Was loop 
!
            SGI = (SJLOC-1)*SBNLOC+SGJ
!
      do ICOORD = 1, NCOORD! Was loop 
!
               POS(ICOORD)     = 0.0
               DPDXI(ICOORD)   = 0.0
               DPDETA(ICOORD)  = 0.0
!
      do SKLOC = 1, SBNLOC! Was loop 
!
                  POS(ICOORD) = POS(ICOORD) &
     &                        + CORN(SJLOC,ICOORD,SKLOC)*&
     &                          0.25*(1.+XIP(SKLOC)*XIPGP(SGJ))&
     &                          *(1.+ETAP(SKLOC)*ETAPGP(SGJ))
!
                  DPDXI(ICOORD) = DPDXI(ICOORD) &
     &                          + CORN(SJLOC,ICOORD,SKLOC)&
     &                            *0.25*XIP(SKLOC)&
     &                            *(1.+ETAP(SKLOC)*ETAPGP(SGJ))
!
                  DPDETA(ICOORD) = DPDETA(ICOORD) &
     &                           + CORN(SJLOC,ICOORD,SKLOC)&
     &                             *0.25*(1.+XIP(SKLOC)*XIPGP(SGJ))&
     &                             *ETAP(SKLOC) 
!
               END DO
!
            END DO
!
            L1 = POS(1)
            L2 = POS(2)
            L3 = POS(3)
!
            SN(1,SGI) = (2.0*L1-1.0)*L1
            SN(2,SGI) =  4.0*L1*L2
            SN(3,SGI) = (2.0*L2-1.0)*L2
            SN(4,SGI) =  4.0*L2*L3
            SN(5,SGI) = (2.0*L3-1.0)*L3
            SN(6,SGI) =  4.0*L1*L3
!
            SNLX(1,SGI) = (4.0*L1-1.0)*DPDXI(1)
            SNLX(2,SGI) =  4.0*DPDXI(1)*L2 + 4.0*L1*DPDXI(2)
            SNLX(3,SGI) = (4.0*L2-1.0)*DPDXI(2)
            SNLX(4,SGI) =  4.0*DPDXI(2)*L3 + 4.0*L2*DPDXI(3)
            SNLX(5,SGI) = (4.0*L3-1.0)*DPDXI(3)
            SNLX(6,SGI) =  4.0*DPDXI(1)*L3 + 4.0*L1*DPDXI(3)
!
            SNLY(1,SGI) = (4.0*L1-1.0)*DPDETA(1)
            SNLY(2,SGI) =  4.0*DPDETA(1)*L2 + 4.0*L1*DPDETA(2)
            SNLY(3,SGI) = (4.0*L2-1.0)*DPDETA(2)
            SNLY(4,SGI) =  4.0*DPDETA(2)*L3 + 4.0*L2*DPDETA(3)
            SNLY(5,SGI) = (4.0*L3-1.0)*DPDETA(3)
            SNLY(6,SGI) =  4.0*DPDETA(1)*L3 + 4.0*L1*DPDETA(3)
!
         END DO
!
      END DO
!
      do SGI = 1, SNGI! Was loop 
!
         SWEIGH(SGI) = 1.0
!
      END DO
!
      END SUBROUTINE

end module legacy_cv_shape_functions
