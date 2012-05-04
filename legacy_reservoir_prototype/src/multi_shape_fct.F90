
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
#include "fdebug.h"


!!!!========================================!!!!
!!!!     SHAPE FUNCTIONS SUBRTS             !!!!
!!!!========================================!!!!


  module shape_functions

    use fldebug
    use shape_functions_Linear_Quadratic
    use shape_functions_NDim

    use state_module
    use spud
    use global_parameters, only: option_path_len

  contains

!!!
!!! SHAPESE AND RELATED SUBRTS & FUNCTIONS

    SUBROUTINE SHAPESE( SELETYP,  SNEILOC,    &
         SNGI, SNLOC,   &
         SN, &
         SNLX, SNLY,  &
         SWEIGH )
      !- this subroutine generates the FE basis functions, weights and the
      !- derivatives of the surface shape functions for a variety of elements. 

      IMPLICIT NONE
      INTEGER, intent( in ) :: SELETYP, SNGI, SNLOC
      INTEGER, DIMENSION( SNLOC, SNGI ), intent( inout ) :: SNEILOC
      REAL, DIMENSION( SNLOC, SNGI ), intent( inout ) :: SN, SNLX, SNLY
      REAL, DIMENSION( SNGI ), intent( inout ) ::  SWEIGH

      ewrite(3,*) 'In SHAPESE'

      SELECT CASE( SELETYP )

      CASE( 1 ) ; CALL SFELIN( SNGI, SNLOC, SN, SNLX, SWEIGH )

      CASE( 2 ) ; CALL SFEQUAD( SNGI, SNLOC, &
           SN, SNLX,  &
           SNLY, SWEIGH )

      CASE( 3 ) ; CALL SFETRI( SNGI, SNLOC, &
           SN, SNLX,  SNLY, &
           SWEIGH )

      CASE( 4 ) ; CALL SFEQLIN( SNGI, SNLOC, SN, SNLX, SWEIGH )

      CASE( 5 ) ; CALL SFEQQUAD( SNGI, SNLOC, &
           SN, SNLX, &
           SNLY, SWEIGH )

      CASE( 6 ) ; CALL SFEQTRI( SNGI, SNLOC, &
           SN, SNLX, &
           SNLY, SWEIGH )

      END SELECT

      CALL SURNEI( SNEILOC, SNLOC, SNGI, SELETYP )

      ewrite(3,*) 'Leaving SHAPESE'


      RETURN
    END SUBROUTINE SHAPESE



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
      INTEGER, intent( in ) :: SELETYP, SNGI, SNLOC
      !     
      INTEGER, DIMENSION( SNLOC, SNGI ), intent( inout )  :: SNEILOC
      !     
      !     - local variables
      !     
      INTEGER :: COUNT, INC, SILOC, SGI

      ewrite(3,*) 'In SURNEI'


      SNEILOC = 0

      SELECT CASE ( SELETYP )

      CASE( 1 : 3 ) ! Linear line (1), quadrilateral (2), triangular (3) surface elements

         COUNT = 1
         DO SILOC = 1, SNLOC

            DO SGI = COUNT, COUNT
               SNEILOC( SILOC, SGI ) = 1
            END DO

            COUNT = COUNT + 1

         END DO
         !
         ! This is a little note for the future if we are: using 16 point quadrature
         ! for the surface elements then change the loop DO SGI = COUNT, COUNT to 
         ! DO SGI = COUNT, COUNT + 3. Also change the loop increment to 
         ! COUNT = COUNT + 4. This is similar for triangles.
         !

      CASE( 4 ) ! Quadratic line surface elements

         COUNT = 1          
         DO SILOC = 1, SNLOC

            DO SGI = COUNT, COUNT + 1
               SNEILOC( SILOC, SGI ) = 1
            END DO

            COUNT = COUNT + 2

         END DO

      CASE( 5 ) ! Quadratic quadrilateral surface elements

         COUNT = 1          
         DO SILOC = 1, SNLOC

            DO SGI = COUNT, COUNT + 3
               SNEILOC( SILOC, SGI ) = 1
            END DO

            COUNT = COUNT + 4

         END DO


      CASE( 6 ) ! Quadratic triangular surface elements

         COUNT = 1

         Loop_SNLOC: DO SILOC = 1, SNLOC 

            IF(( SILOC == 1 ) .OR. ( SILOC == 3 ) .OR. ( SILOC == 5 )) THEN
               INC = 3
            ELSE
               INC = 7
            ENDIF

            DO SGI = COUNT, COUNT + INC
               SNEILOC( SILOC, SGI ) = 1
            END DO

            IF(( SILOC == 1 ) .OR. ( SILOC == 3 ) .OR. ( SILOC == 5 )) THEN
               COUNT = COUNT + 4
            ELSE
               COUNT = COUNT + 8
            ENDIF

         END DO Loop_SNLOC
         !     
         ! - note that in future that the DO GI = COUNT, COUNT+3
         ! - will become DO GI = COUNT, COUNT+NGI/NLOC
         !
         ! - note the use of NGI/NLOC which gives the number of Gauss
         ! - points associated with a CV. We might consider the use
         ! - of another variable NCVGI which gives the number of Gauss
         ! - points associated with a CV.

      CASE DEFAULT
         FLAbort(" Wrong option for surface element type ")

      END SELECT

      ewrite(3,*) 'Leaving SURNEI'

      RETURN

    END SUBROUTINE SURNEI
    !
    !
    !
    !

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
      INTEGER, intent( in ) :: SNGI, SNLOC
      !
      REAL, DIMENSION( SNLOC, SNGI ), intent( inout ) :: SN, SNLX
      !
      REAL, DIMENSION( SNGI ), intent( inout ) :: SWEIGH
      !
      ! - local variables
      !
      INTEGER ::SILOC, SJLOC, SGI, SGJ
      !
      REAL, DIMENSION( : ) , allocatable :: XIGP, XI
      INTEGER, PARAMETER :: ONE = 1, TWO = 2

      ewrite(3,*) 'In SFELIN'

      ALLOCATE( XIGP( ONE ))
      ALLOCATE( XI( TWO ))
      !
      ! - positions of Gauss points associated with CV's in (ETAP) space
      ! 
      XIGP = 0.0
      !
      ! - positions of vertices of line element in (ETA) space
      !
      XI(1) = -1.0
      XI(2) =  1.0
      !
      ! - values of FE basis functions at volume integration points
      !
      do SILOC = 1, SNLOC
         do SJLOC = 1, SNLOC
            do SGJ = 1, 1

               SGI = ( SJLOC - 1 ) * ( SNLOC - 1 ) + SGJ
               SN( SILOC, SGI ) = 0.5 * ( 1. + XI(SILOC) &
                    * 0.5 * ( XI( SJLOC ) + XIGP( SGJ )))
            END DO
         END DO
      END DO
      !
      ! - derivatives of FE basis functions at the Gauss points
      !
      do SILOC = 1, SNLOC
         do SGI = 1, SNGI
            SNLX( SILOC, SGI ) = 0.5 * XI( SILOC )
         END DO
      END DO
      !
      ! - note that the weight of the normal Gauss point is 2.0
      ! - but this is reduced by 0.5 i.e. 1/2 because of the 
      ! - second Jacobian
      !
      do SGI = 1, SNGI
         SWEIGH( SGI ) = 1.0
      END DO

      DEALLOCATE( XIGP )
      DEALLOCATE( XI )

      ewrite(3,*) 'Leaving SFELIN'

      RETURN
      !
    END SUBROUTINE SFELIN

    SUBROUTINE SFEQUAD( SNGI, SNLOC, &
                                ! - REALS
         SN,   SNLX, &
         SNLY, SWEIGH )
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
      INTEGER, PARAMETER :: NCOORD = 2, SBNLOC = 4, SBNGI = 1 
      REAL, DIMENSION( : ), allocatable :: POS, DPDXI, DPDETA, XIGP, ETAGP, &
           XI, ETA
      REAL, DIMENSION( : , : , : ), allocatable :: CORN

      ewrite(3,*) 'In SFEQUAD'

      ALLOCATE( POS( NCOORD ))
      ALLOCATE( DPDXI( NCOORD ))
      ALLOCATE( DPDETA( NCOORD ))
      ALLOCATE( XIGP( SBNGI ))
      ALLOCATE( ETAGP( SBNGI ))
      ALLOCATE( XI( SBNLOC ))
      ALLOCATE( ETA( SBNLOC ))
      ALLOCATE( CORN( SNLOC, NCOORD, SBNLOC ))

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
      Loop_SILOC: DO SILOC = 1, SNLOC

         Loop_SJLOC: do SJLOC = 1, SNLOC

            Loop_SGJ: DO SGJ = 1, SBNGI

               SGI = SJLOC-1 + SGJ

               Loop_ICOORD: DO ICOORD = 1, NCOORD

                  POS(ICOORD)     = 0.0
                  DPDXI(ICOORD)   = 0.0
                  DPDETA(ICOORD)  = 0.0

                  Loop_SKLOC: do SKLOC = 1, SBNLOC

                     POS(ICOORD) = POS(ICOORD) &
                          + CORN(SJLOC,ICOORD,SKLOC)&
                          *0.25*(1.+XI(SKLOC)*XIGP(SGJ))&
                          *(1.+ETA(SKLOC)*ETAGP(SGJ))

                     DPDXI(ICOORD) = DPDXI(ICOORD) &
                          + CORN(SJLOC,ICOORD,SKLOC)&
                          *0.25*XI(SKLOC)&
                          *(1.+ETA(SKLOC)*ETAGP(SGJ))

                     DPDETA(ICOORD) = DPDETA(ICOORD) &
                          + CORN(SJLOC,ICOORD,SKLOC)&
                          *0.25*(1.+XI(SKLOC)*XIGP(SGJ))&
                          *ETA(SKLOC)

                  END DO Loop_SKLOC

               END DO Loop_ICOORD
               !
               SN(SILOC,SGI) = 0.25*( 1.0 + XI(SILOC)*POS(1) )&
                    *( 1.0 + ETA(SILOC)*POS(2) )
               !
               SNLX(SILOC,SGI) = 0.25*XI(SILOC)*DPDXI(1)&
                    *( 1.0 + ETA(SILOC)*POS(2) )&
                    + 0.25*( 1.0 + XI(SILOC)*POS(1) )&
                    *ETA(SILOC)*DPDXI(2)
               !
               SNLY(SILOC,SGI) = 0.25*XI(SILOC)*DPDETA(1)&
                    *( 1.0 + ETA(SILOC)*POS(2) )&
                    + 0.25*( 1.0 + XI(SILOC)*POS(1) )&
                    *ETA(SILOC)*DPDETA(2)
               !
            END DO Loop_SGJ   ! sjloc
            !
         END DO Loop_SJLOC  ! sjloc
         !
      END DO Loop_SILOC  ! siloc

      do SGI = 1, SNGI! Was loop 
         !
         SWEIGH(SGI) = 4.0
         !
      END DO


      DEALLOCATE( POS )
      DEALLOCATE( DPDXI )
      DEALLOCATE( DPDETA )
      DEALLOCATE( XIGP )
      DEALLOCATE( ETAGP )
      DEALLOCATE( XI )
      DEALLOCATE( ETA )
      DEALLOCATE( CORN )

      ewrite(3,*) 'Leaving SFEQUAD'


      RETURN
      !
    END SUBROUTINE SFEQUAD


    SUBROUTINE SFETRI( SNGI, SNLOC, &
         SN,   SNLX, &
         SNLY, SWEIGH )
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
      !
      INTEGER, PARAMETER :: NCOORD = 3, SBNLOC = 4, SBNGI = 1 
      REAL, DIMENSION( : ), allocatable :: POS, DPDXI, DPDETA, XIGP, ETAGP, &
           XI, ETA
      REAL, DIMENSION( : , : , : ), allocatable :: CORN

      ewrite(3,*) 'In SFETRI'

      ALLOCATE( POS( NCOORD ))
      ALLOCATE( DPDXI( NCOORD ))
      ALLOCATE( DPDETA( NCOORD ))
      ALLOCATE( XIGP( SBNGI ))
      ALLOCATE( ETAGP( SBNGI ))
      ALLOCATE( XI( SBNLOC ))
      ALLOCATE( ETA( SBNLOC ))
      ALLOCATE( CORN( SNLOC, NCOORD, SBNLOC ))

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
      Loop_SNLOCI: DO SILOC = 1, SNLOC

         Loop_SNLOCJ: DO SJLOC = 1, SNLOC

            Loop_SG: DO SGJ = 1, SBNGI
               !
               SGI = SJLOC-1 + SGJ
               !
               Loop_COORD: DO ICOORD = 1, NCOORD
                  !
                  POS(ICOORD)     = 0.0
                  DPDXI(ICOORD)   = 0.0
                  DPDETA(ICOORD)  = 0.0
                  !
                  Loop_SNLOCK: DO SKLOC = 1, SBNLOC
                     !
                     POS(ICOORD) = POS(ICOORD) &
                          + CORN(SJLOC,ICOORD,SKLOC)*&
                          0.25*(1.+XI(SKLOC)*XIGP(SGJ))&
                          *(1.+ETA(SKLOC)*ETAGP(SGJ))


                     DPDXI(ICOORD) = DPDXI(ICOORD) &
                          + CORN(SJLOC,ICOORD,SKLOC)&
                          *0.25*XI(SKLOC)&
                          *(1.+ETA(SKLOC)*ETAGP(SGJ))

                     DPDETA(ICOORD) = DPDETA(ICOORD) &
                          + CORN(SJLOC,ICOORD,SKLOC)&
                          *0.25*(1.+XI(SKLOC)*XIGP(SGJ))&
                          *ETA(SKLOC)

                  END DO Loop_SNLOCK

               END DO Loop_COORD

               SN(SILOC,SGI)   = POS(SILOC)

               SNLX(SILOC,SGI) = DPDXI(SILOC)

               SNLY(SILOC,SGI) = DPDETA(SILOC)

            END DO Loop_SG
            !
         END DO Loop_SNLOCJ
         !
      END DO Loop_SNLOCI
      !
      do SGI = 1, SNGI 
         !
         SWEIGH(SGI) = 4.0
         !
      END DO

      DEALLOCATE( POS )
      DEALLOCATE( DPDXI )
      DEALLOCATE( DPDETA )
      DEALLOCATE( XIGP )
      DEALLOCATE( ETAGP )
      DEALLOCATE( XI )
      DEALLOCATE( ETA )
      DEALLOCATE( CORN )

      ewrite(3,*) 'Leaving SFETRI'

      RETURN
      !
    END SUBROUTINE SFETRI


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

      ewrite(3,*) 'In SFEQLIN'
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
      XIPGP(1) = -1./SQRT(3.0)
      XIPGP(2) =  1./SQRT(3.0)
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
      do SIFACE = 1, SNFACE
         !
         do SGJ = 1, SFNGI
            !
            SGI = (SIFACE-1)*SFNGI+SGJ
            !
            POS   = 0.0
            DPDXI = 0.0
            !
            do SJLOC = 1, SFNLOC

               POS   = POS &
                    + CORN(SIFACE,SJLOC)&
                    *0.5*(1.+XIP(SJLOC)*XIPGP(SGJ))

               DPDXI = DPDXI&
                    + CORN(SIFACE,SJLOC)*0.5*XIP(SJLOC)

            END DO

            XI = POS

            LIJXI(1)   =  0.5*XI*(XI-1.0)
            LIJXI(2)   =  1.0-XI*XI
            LIJXI(3)   =  0.5*XI*(XI+1.0)

            DLIJXIDXI(1) =  0.5*(2.0*XI-1.0)*DPDXI
            DLIJXIDXI(2) = -2.0*XI*DPDXI
            DLIJXIDXI(3) =  0.5*(2.0*XI+1.0)*DPDXI

            do SILOC = 1, SNLOC! Was loop 

               SN(SILOC,SGI)   = LIJXI(SILOC)

               SNLX(SILOC,SGI) = DLIJXIDXI(SILOC)

            END DO

         END DO

      END DO
      !
      ! - set the weights for the surface integration
      !
      do SGI = 1, SNGI! Was loop 
         !
         SWEIGH(SGI) = 1.0
         !
      END DO

      ewrite(3,*) 'Leaving SFEQLIN'

      RETURN
      !
    END SUBROUTINE SFEQLIN

    SUBROUTINE SFEQQUAD( SNGI, SNLOC, &
                                ! - REALS
         SN,   SNLX, &
         SNLY, SWEIGH )
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
      INTEGER, PARAMETER :: NCOORD = 2, SBNLOC = 4, SBNGI = 4,  NSCV = 9, THREE = 3
      REAL, DIMENSION( : ), allocatable :: POS, DPDXI, DPDETA, XIGP, ETAGP, XIP, ETAP, &
           LIJXI, LIJETA, DLIJXIDXI, DLIJXIDETA, DLIJETADXI, DLIJETADETA
      REAL, DIMENSION( : , : , : ), allocatable :: CORN
      REAL :: XI, ETA
      INTEGER ::  I, J

      ewrite(3,*) 'In SFEQQUAD'

      ALLOCATE( POS( NCOORD ))
      ALLOCATE( DPDXI( NCOORD ))
      ALLOCATE( DPDETA( NCOORD ))
      ALLOCATE( XIGP( SBNGI ))
      ALLOCATE( ETAGP( SBNGI ))
      ALLOCATE( XIP( SBNLOC ))
      ALLOCATE( ETAP( SBNLOC ))
      ALLOCATE( CORN( NSCV, NCOORD, SBNLOC ))
      ALLOCATE( LIJXI( THREE )) ! for derivatives wrt the second parametric space.
      ALLOCATE( LIJETA( THREE ))
      ALLOCATE( DLIJXIDXI( THREE ))
      ALLOCATE( DLIJXIDETA( THREE ))
      ALLOCATE( DLIJETADXI( THREE ))
      ALLOCATE( DLIJETADETA( THREE ))

      !
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
                       + CORN(SJLOC,ICOORD,SKLOC)&
                       *0.25*(1.+XIP(SKLOC)*XIGP(SGJ))&
                       *(1.+ETAP(SKLOC)*ETAGP(SGJ))
                  !
                  DPDXI(ICOORD) = DPDXI(ICOORD) &
                       + CORN(SJLOC,ICOORD,SKLOC)&
                       *0.25*XIP(SKLOC)&
                       *(1.+ETAP(SKLOC)*ETAGP(SGJ))
                  !
                  DPDETA(ICOORD) = DPDETA(ICOORD) &
                       + CORN(SJLOC,ICOORD,SKLOC)&
                       *0.25*(1.+XIP(SKLOC)*XIGP(SGJ))&
                       *ETAP(SKLOC)
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
                       + DLIJXIDXI(I)*LIJETA(J)
                  !
                  SNLY(SILOC,SGI) = LIJXI(I)*DLIJETADETA(J)&
                       + DLIJXIDETA(I)*LIJETA(J)
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

      DEALLOCATE( POS )
      DEALLOCATE( DPDXI )
      DEALLOCATE( DPDETA )
      DEALLOCATE( XIGP )
      DEALLOCATE( ETAGP )
      DEALLOCATE( XIP )
      DEALLOCATE( ETAP )
      DEALLOCATE( CORN )
      DEALLOCATE( LIJXI ) ! for derivatives wrt the second parametric space.
      DEALLOCATE( LIJETA )
      DEALLOCATE( DLIJXIDXI )
      DEALLOCATE( DLIJXIDETA )
      DEALLOCATE( DLIJETADXI )
      DEALLOCATE( DLIJETADETA )

      ewrite(3,*) 'Leaving SFEQQUAD'

      RETURN
      !
    END SUBROUTINE SFEQQUAD




    SUBROUTINE SFEQTRI( SNGI,   SNLOC, &
                                ! - REALS
         SN,     SNLX, &
         SNLY,   SWEIGH )
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
      INTEGER  SJLOC, SKLOC, ICOORD
      !
      INTEGER SGI, SGJ
      !
      ! - note that NCOORD is the number of co-ordinates,
      ! - SBNLOC is the number of local nodes for a subcell
      ! - and SBNGI is the number of local integration points.
      !
      INTEGER, PARAMETER :: NCOORD = 3, SBNLOC = 4, SBNGI = 4,  NSCV = 9, THREE = 3
      REAL, DIMENSION( : ), allocatable :: POS, DPDXI, DPDETA, XIPGP, ETAPGP, XIP, ETAP
      REAL, DIMENSION( : , : , : ), allocatable :: CORN
      REAL ::  L1, L2, L3

      ewrite(3,*) 'In SFEQTRI'

      ALLOCATE( POS( NCOORD ))
      ALLOCATE( DPDXI( NCOORD ))
      ALLOCATE( DPDETA( NCOORD ))
      ALLOCATE( XIPGP( SBNGI ))
      ALLOCATE( ETAPGP( SBNGI ))
      ALLOCATE( CORN( NSCV, NCOORD, SBNLOC ))
      ALLOCATE( XIP( SBNLOC ))
      ALLOCATE( ETAP( SBNLOC ))
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
                       + CORN(SJLOC,ICOORD,SKLOC)*&
                       0.25*(1.+XIP(SKLOC)*XIPGP(SGJ))&
                       *(1.+ETAP(SKLOC)*ETAPGP(SGJ))
                  !
                  DPDXI(ICOORD) = DPDXI(ICOORD) &
                       + CORN(SJLOC,ICOORD,SKLOC)&
                       *0.25*XIP(SKLOC)&
                       *(1.+ETAP(SKLOC)*ETAPGP(SGJ))
                  !
                  DPDETA(ICOORD) = DPDETA(ICOORD) &
                       + CORN(SJLOC,ICOORD,SKLOC)&
                       *0.25*(1.+XIP(SKLOC)*XIPGP(SGJ))&
                       *ETAP(SKLOC) 
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

      DEALLOCATE( POS )
      DEALLOCATE( DPDXI )
      DEALLOCATE( DPDETA )
      DEALLOCATE( XIPGP )
      DEALLOCATE( ETAPGP )
      DEALLOCATE( CORN )
      DEALLOCATE( XIP )
      DEALLOCATE( ETAP )

      ewrite(3,*) 'Leaving SFEQTRI'

      RETURN
      !
    END SUBROUTINE SFEQTRI
    !

!!!
!!!   SHAPESV AND RELATED SUBRTS & FUNCTIONS
!!!

    subroutine shape_cv_n( ndim, cv_ele_type, &
         cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, &
         n, nlx, nly, nlz, &
         un, unlx, unly, unlz )
      ! Shape functions associated with volume integration using both CV basis 
      ! functions CVN as well as FEM basis functions N (and its derivatives NLX, NLY, NLZ)
      ! also for velocity basis functions UN, UNLX, UNLY, UNLZ
      implicit none
      integer, intent( in ) :: ndim, cv_ele_type, cv_ngi, cv_nloc, u_nloc
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: cvn
      real, dimension( cv_ngi ), intent( inout ) :: cvweigh
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: n, nlx, nly, nlz
      real, dimension( u_nloc, cv_ngi ), intent( inout ) :: un, unlx, unly, unlz

      ewrite(3,*) 'In SHAPE_CV_N'

      Select Case( cv_ele_type )
      case( 1, 2 ) ! 1D
         call quad_1d_shape( cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, n, nlx, un, unlx ) 
         nly = 0. 
         nlz = 0. 
         unly = 0. 
         unlz = 0. 

      case( 5, 6, 9, 10 ) ! Quadrilaterals and Hexahedra
         call quad_nd_shape( ndim, cv_ele_type, cv_ngi, cv_nloc, u_nloc, cvn, cvweigh, &
              n, nlx, nly, nlz, &
              un, unlx, unly, unlz )

      case( 3, 4, 7, 8 ) ! Triangles and Tetrahedra
         call vol_cv_tri_tet_shape( cv_ele_type, ndim, cv_ngi, cv_nloc, u_nloc, cvn, &
              cvweigh, n, nlx, nly, nlz, &
              un, unlx, unly, unlz )
         !stop 12
      case default; FLExit( "Wrong integer for CV_ELE_TYPE" )
      end Select

      ewrite(3,*) 'Leaving SHAPE_CV_N'

    end subroutine shape_cv_n


    subroutine cv_fem_shape_funs( &
         ndim, cv_ele_type, &
         cv_ngi, cv_ngi_short, cv_nloc, u_nloc, cvn, cvn_short, &
                                ! Volume shape functions
         cvweight, cvfen, cvfenlx, cvfenly, cvfenlz, &
         cvweight_short, cvfen_short, cvfenlx_short, cvfenly_short, cvfenlz_short, &
         ufen, ufenlx, ufenly, ufenlz, &
                                ! Surface of each CV shape functions
         scvngi, cv_neiloc, cv_on_face, cvfem_on_face, &
         scvfen, scvfenslx, scvfensly, scvfeweigh, &
         scvfenlx, scvfenly, scvfenlz, &
         sufen, sufenslx, sufensly, &
         sufenlx, sufenly, sufenlz, &
                                ! Surface element shape funcs
         u_on_face, ufem_on_face, nface, &
         sbcvngi, sbcvfen, sbcvfenslx, sbcvfensly, sbcvfeweigh, sbcvfenlx, sbcvfenly, sbcvfenlz, &
         sbufen, sbufenslx, sbufensly, sbufenlx, sbufenly, sbufenlz, &
         cv_sloclist, u_sloclist, cv_snloc, u_snloc, &
                                ! Define the gauss points that lie on the surface of the CV
         findgpts, colgpts, ncolgpts, &
         sele_overlap_scale, QUAD_OVER_WHOLE_ELE ) 
      ! This subrt defines the sub-control volume and FEM shape functions.
      ! Shape functions associated with volume integration using both CV basis 
      ! functions CVN as well as FEM basis functions CVFEN (and its derivatives 
      ! CVFENLX, CVFENLY, CVFENLZ)
      implicit none
      integer, intent( in ) :: ndim, cv_ele_type, cv_ngi, cv_ngi_short, cv_nloc, u_nloc
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: cvn
      real, dimension( cv_nloc, cv_ngi_short ), intent( inout ) :: cvn_short
      real, dimension( cv_ngi ), intent( inout ) :: cvweight
      real, dimension( cv_nloc, cv_ngi ), intent( inout ) :: cvfen, cvfenlx, cvfenly, cvfenlz
      real, dimension( cv_ngi_short ), intent( inout ) :: cvweight_short
      real, dimension( cv_nloc, cv_ngi_short ), intent( inout ) :: cvfen_short, cvfenlx_short, &
           cvfenly_short, cvfenlz_short
      real, dimension( u_nloc, cv_ngi ), intent( inout ) :: ufen, ufenlx, ufenly, ufenlz
      integer, intent( in ) :: scvngi
      integer, dimension( cv_nloc, scvngi ), intent( inout ) :: cv_neiloc
      logical, dimension( cv_nloc, scvngi ), intent( inout ) :: cv_on_face, cvfem_on_face
      real, dimension( cv_nloc, scvngi ), intent( inout ) :: scvfen, scvfenslx, scvfensly
      real, dimension( scvngi ), intent( inout ) :: scvfeweigh
      real, dimension( cv_nloc, scvngi ), intent( inout ) :: scvfenlx, scvfenly, scvfenlz
      real, dimension( u_nloc, scvngi ), intent( inout ) :: sufen, sufenslx, sufensly, sufenlx, &
           sufenly, sufenlz
      logical, dimension( u_nloc, scvngi ), intent( inout ) :: u_on_face, ufem_on_face
      integer, intent( in ) :: nface, sbcvngi
      logical, intent( in ) :: QUAD_OVER_WHOLE_ELE
      ! if QUAD_OVER_WHOLE_ELE then dont divide element into CV's to form quadrature.
      real, dimension( cv_snloc, sbcvngi ), intent( inout ) :: sbcvfen, sbcvfenslx, sbcvfensly
      real, dimension( sbcvngi ), intent( inout ) :: sbcvfeweigh
      real, dimension( cv_snloc, sbcvngi ), intent( inout ) :: sbcvfenlx, sbcvfenly, sbcvfenlz
      integer, intent( in ) :: cv_snloc, u_snloc
      real, dimension( u_snloc, sbcvngi ), intent( inout ) :: sbufen, sbufenslx, sbufensly, &
           sbufenlx, sbufenly, sbufenlz
      integer, dimension( nface, cv_snloc ), intent( inout ) :: cv_sloclist
      integer, dimension( nface, u_snloc ), intent( inout ) :: u_sloclist
      integer, dimension( cv_nloc + 1 ), intent( inout ) :: findgpts
      integer, dimension( cv_nloc * scvngi ), intent( inout ) :: colgpts
      integer, intent( inout ) :: ncolgpts
      real, dimension( cv_nloc ), intent( inout ) :: sele_overlap_scale
      ! Local variables
      logical, dimension( :, : ), allocatable :: u_on_face2, ufem_on_face2
      integer, dimension( :, : ), allocatable :: u_sloclist2
      real, dimension( :, : ), allocatable :: ufen2, ufenlx2, ufenly2, ufenlz2, & 
           sufen2, sufenslx2, sufensly2, sufenlx2, sufenly2, sufenlz2, &
           sbufen2, sbufenslx2, sbufensly2, sbufenlx2, sbufenly2, sbufenlz2 
      character( len = option_path_len ) :: overlapping_path 
      logical :: is_overlapping   
      integer :: u_nloc2, ilev, ilev2, u_snloc2, u_ele_type2, gi, cv_iloc

      ewrite(3,*) 'in  cv_fem_shape_funs subrt'

      is_overlapping = .false.
      call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           overlapping_path )
      if( trim( overlapping_path ) == 'overlapping' ) is_overlapping = .true.

      ! Sele_Overlap_Scale is the scaling needed to convert to overlapping element surfaces
      if( is_overlapping ) then
         Select Case( cv_nloc )
         case( 3 :  )
            sele_overlap_scale = 4. 
         case( 2 )
            sele_overlap_scale = 2.
         case default
            sele_overlap_scale = 1.
         end Select
      else
         sele_overlap_scale = 1.
      end if

      Conditional_OverlappingMethod1: if( is_overlapping ) then
         ! Define basis functions (overlapping) for porous media flow
         u_nloc2 = u_nloc / cv_nloc

         allocate( ufen2( u_nloc2, cv_ngi_short ) )
         allocate( ufenlx2( u_nloc2, cv_ngi_short ) )
         allocate( ufenly2( u_nloc2, cv_ngi_short ) )
         allocate( ufenlz2( u_nloc2, cv_ngi_short ) )

         call shape_cv_n( ndim, cv_ele_type, &
              cv_ngi_short, cv_nloc, u_nloc2, cvn_short, cvweight_short, &
              cvfen_short, cvfenlx_short, cvfenly_short, cvfenlz_short, &
              ufen2, ufenlx2, ufenly2, ufenlz2 )

         ! Defining the base functions for each level
         ufen = 0. ; ufenlx = 0. ; ufenly = 0. ; ufenlz = 0. ; cvn = 0.
         Loop_ILEV1: do ilev = 1, cv_nloc
            ufen( 1 + ( ilev - 1) * u_nloc2 : ilev * u_nloc2, &
                 1 + ( ilev - 1 ) * cv_ngi_short : ilev * cv_ngi_short ) = &
                 ufen2( 1 : u_nloc2, 1 : cv_ngi_short )
            ufenlx( 1 + ( ilev - 1) * u_nloc2 : ilev * u_nloc2, &
                 1 + ( ilev - 1 ) * cv_ngi_short : ilev * cv_ngi_short ) = &
                 ufenlx2( 1 : u_nloc2, 1 : cv_ngi_short )
            ufenly( 1 + ( ilev - 1) * u_nloc2 : ilev * u_nloc2, &
                 1 + ( ilev - 1 ) * cv_ngi_short : ilev * cv_ngi_short ) = &
                 ufenly2( 1 : u_nloc2, 1 : cv_ngi_short )
            ufenlz( 1 + ( ilev - 1) * u_nloc2 : ilev * u_nloc2, &
                 1 + ( ilev - 1 ) * cv_ngi_short : ilev * cv_ngi_short ) = &
                 ufenlz2( 1 : u_nloc2, 1 : cv_ngi_short )

            cvfen( 1 : cv_nloc, &
                 1 + ( ilev - 1 ) * cv_ngi_short : ilev * cv_ngi_short ) = &
                 cvfen_short( 1 : cv_nloc, 1 : cv_ngi_short )
            cvfenlx( 1 : cv_nloc, &
                 1 + ( ilev - 1 ) * cv_ngi_short : ilev * cv_ngi_short ) = &
                 cvfenlx_short( 1 : cv_nloc, 1 : cv_ngi_short )
            cvfenly( 1 : cv_nloc, &
                 1 + ( ilev - 1 ) * cv_ngi_short : ilev * cv_ngi_short ) = &
                 cvfenly_short( 1 : cv_nloc, 1 : cv_ngi_short )
            cvfenlz( 1 : cv_nloc, &
                 1 + ( ilev - 1 ) * cv_ngi_short : ilev * cv_ngi_short ) = &
                 cvfenlz_short( 1 : cv_nloc, 1 : cv_ngi_short )

            cvweight( 1 + ( ilev - 1 ) * cv_ngi_short : &
                 ilev * cv_ngi_short ) = &
                 cvweight_short( 1 : cv_ngi_short )
            cvn( ilev, &
                 1 + ( ilev - 1 ) * cv_ngi_short : ilev * cv_ngi_short ) = 1.

         end do Loop_ILEV1

      else ! if it is not overlapping formulation
         call shape_cv_n( ndim, cv_ele_type, &
              cv_ngi, cv_nloc, u_nloc, cvn, cvweight, &
              cvfen, cvfenlx, cvfenly, cvfenlz, &
              ufen, ufenlx, ufenly, ufenlz )
         cvn_short = cvn
         cvfen_short = cvfen
         cvfenlx_short = cvfenlx
         cvfenly_short = cvfenly
         cvfenlz_short = cvfenlz
         cvweight_short = cvweight       

      end if Conditional_OverlappingMethod1
      ewrite(3,*) 'out of shape_cv_n - CVWEIGHT', CVWEIGHT

      !
      !(a) scvfen( cv_nloc, scvngi ): the shape function evaluated for each node 
      !          at each surface gauss point
      !(b) scvfenslx[y/z]( cv_nloc, scvngi ): the surface derivatives of the shape 
      !          function for each node at those same points, and the derivatives 
      !          of the shape
      !(c) scvfeweigh( scvngi ): the Gauss weights to use when integrating around 
      !          the control volume surface
      !(d) cv_neiloc( cv_nloc, scvngi ): neighbour node for a given node/gauss-point
      !          pair. This also include quadature points around the element. 
      !

      Conditional_OverlappingMethod2: if( is_overlapping ) then
         u_ele_type2 = 1
         u_nloc2 = u_nloc / cv_nloc
         allocate( ufem_on_face2( u_nloc2, scvngi ) )
         allocate( sufen2( u_nloc2, scvngi ) )
         allocate( sufenslx2( u_nloc2, scvngi ) )
         allocate( sufensly2( u_nloc2, scvngi ) )
         allocate( sufenlx2( u_nloc2, scvngi ) )
         allocate( sufenly2( u_nloc2, scvngi ) )
         allocate( sufenlz2( u_nloc2, scvngi ) )
         ewrite(3,*)'cv_nloc, cv_ngi, scvngi:', cv_nloc, cv_ngi, scvngi
         ewrite(3,*)'u_nloc, u_nloc2:', u_nloc, u_nloc2


         call shapesv_fem_plus( scvngi, cv_neiloc, cv_on_face, cvfem_on_face, &
              ufem_on_face2, &
              cv_ele_type, cv_nloc, scvfen, scvfenslx, scvfensly, scvfeweigh, &
              scvfenlx, scvfenly, scvfenlz, &
              u_nloc2, sufen2, sufenslx2, sufensly2, &
              sufenlx2, sufenly2, sufenlz2, &
              ndim )

         !call U_Volnei( cv_ele_type, cv_nloc, u_nloc, scvngi, &
         !     cv_neiloc,   &
         !     u_on_face )

         !         do ilev = 1, u_nloc
         !            ewrite(3,*)'**u_on_face**:', ( u_on_face( ilev, gi ), gi = 1, scvngi )
         !         end do
         !         do ilev = 1, cv_nloc
         !            ewrite(3,*)'**cvfem_on_face**:', ( cvfem_on_face( ilev, gi ), gi = 1, scvngi )
         !         end do

         sufen = 0. ; sufenslx = 0. ; sufensly = 0. 
         sufenlx = 0. ; sufenly = 0. ; sufenlz = 0.
         Loop_ILEV2: do ilev = 1, cv_nloc
            !==  u_on_face( 1 + ( ilev - 1 ) * u_nloc2 : ilev * u_nloc2, &
            !==       1 : scvngi ) = &
            !==       u_on_face2( 1 : u_nloc2, 1 : scvngi )
            !            u_on_face2( 1 : u_nloc2, 1 : scvngi ) = &
            !                 u_on_face( 1 + ( ilev - 1 ) * u_nloc2 : ilev * u_nloc2, &
            !                 1 : scvngi ) 
            !            ufem_on_face2( 1 : u_nloc2, 1 : scvngi ) = &
            !                 ufem_on_face( 1 + ( ilev - 1 ) * u_nloc2 : ilev * u_nloc2, &
            !                 1 : scvngi ) 
            sufen( 1 + ( ilev - 1 ) * u_nloc2 : ilev * u_nloc2, &
                 1 : scvngi ) = &
                 sufen2( 1 : u_nloc2, 1 : scvngi )
            sufenslx( 1 + ( ilev - 1 ) * u_nloc2 : ilev * u_nloc2, &
                 1 : scvngi ) = &
                 sufenslx2( 1 : u_nloc2, 1 : scvngi )
            sufensly( 1 + ( ilev - 1 ) * u_nloc2 : ilev * u_nloc2, &
                 1 : scvngi ) = &
                 sufensly2( 1 : u_nloc2, 1 : scvngi )
            sufenlx( 1 + ( ilev - 1 ) * u_nloc2 : ilev * u_nloc2, &
                 1 : scvngi ) = &
                 sufenlx2( 1 : u_nloc2, 1 : scvngi )
            sufenly( 1 + ( ilev - 1 ) * u_nloc2 : ilev * u_nloc2, &
                 1 : scvngi ) = &
                 sufenly2( 1 : u_nloc2, 1 : scvngi )
            sufenlz( 1 + ( ilev - 1 ) * u_nloc2 : ilev * u_nloc2, &
                 1 : scvngi ) = &
                 sufenlz2( 1 : u_nloc2, 1 : scvngi )
            ufem_on_face( 1 + ( ilev - 1 ) * u_nloc2 : ilev * u_nloc2, &
                 1 : scvngi ) = ufem_on_face2( 1 : u_nloc2, 1 : scvngi )
         end do Loop_ILEV2

      else
         call shapesv_fem_plus( scvngi, cv_neiloc, cv_on_face, cvfem_on_face, &
              ufem_on_face, &
              cv_ele_type, cv_nloc, scvfen, scvfenslx, scvfensly, scvfeweigh, &
              scvfenlx, scvfenly, scvfenlz, &
              u_nloc, sufen, sufenslx, sufensly, &
              sufenlx, sufenly, sufenlz, &
              ndim )
         !call U_Volnei( cv_ele_type, cv_nloc, u_nloc, scvngi, &
         !     cv_neiloc,   &
         !     u_on_face )
      end if Conditional_OverlappingMethod2

      ! Determine the surface element shape functions from those 
      ! calculated in SHAPESV_FEM_PLUS and also CV_SLOCLIST( NFACE,CV_SNLOC )
      !ewrite(3,*)'u_on_face2( u_iloc, i_scvngi):'
      !do ilev = 1, u_nloc2
      !   ewrite(3,*) ilev, ( u_on_face2( ilev, ilev2 ), ilev2 = 1, scvngi )
      !end do

      Conditional_OverlappingMethod3: if( is_overlapping ) then
         u_ele_type2 = 1
         u_nloc2 = u_nloc / cv_nloc
         u_snloc2 = u_snloc / cv_nloc

         ewrite(3,*) 'U_NLOC  , U_SNLOC  ', u_nloc  , u_snloc
         ewrite(3,*) 'U_NLOC2, U_SNLOC2', u_nloc2, u_snloc2
         ewrite(3,*) 'SBCVNGI                   ', sbcvngi
         ewrite(3,*) 'scvfen:', size( scvfen ) 
         do  cv_iloc = 1, cv_nloc
            ewrite(3,*) ( scvfen( cv_iloc, gi ), gi = 1, scvngi )
         end do

         allocate( u_sloclist2( 1 : nface, u_snloc2 ) )
         allocate( sbufen2( u_snloc2, sbcvngi ) )
         allocate( sbufenslx2( u_snloc2, sbcvngi ) )
         allocate( sbufensly2( u_snloc2, sbcvngi ) )
         allocate( sbufenlx2( u_snloc2, sbcvngi ) )
         allocate( sbufenly2( u_snloc2, sbcvngi ) )
         allocate( sbufenlz2( u_snloc2, sbcvngi ) )

         call det_suf_ele_shape( scvngi, nface, &
              cvfem_on_face, &
              cv_nloc, scvfen, scvfenslx, scvfensly, scvfeweigh, &
              scvfenlx, scvfenly, scvfenlz, &
              u_nloc2, sufen2, sufenslx2, sufensly2, &
              sufenlx2, sufenly2, sufenlz2, &
              sbcvngi, sbcvfen, sbcvfenslx, sbcvfensly, sbcvfeweigh, &
              sbcvfenlx, sbcvfenly, sbcvfenlz, &
              sbufen2, sbufenslx2, sbufensly2, &
              sbufenlx2, sbufenly2, sbufenlz2, &
              cv_sloclist, u_sloclist2, cv_snloc, u_snloc2, &
              ndim, cv_ele_type )
         ewrite(3,*) 'sbufen2:', sbufen2
         ewrite(3,*) 'sbcvfen:', sbcvfen
         ewrite(3,*) 'sufen2:', sufen2
         ewrite(3,*) 'scvfen:', scvfen
         stop 2882

         Loop_ILEV3: do ilev = 1, cv_nloc

            sbufen( 1 + ( ilev - 1 ) * u_snloc2 : ilev * u_snloc2, &
                 1 : sbcvngi ) = &
                 sbufen2( 1 : u_snloc2, 1 : sbcvngi )
            sbufenslx( 1 + ( ilev - 1 ) * u_snloc2 : ilev * u_snloc2, &
                 1 : sbcvngi ) = &
                 sbufenslx2( 1 : u_snloc2, 1 : sbcvngi )
            sbufensly( 1 + ( ilev - 1 ) * u_snloc2 : ilev * u_snloc2, &
                 1 : sbcvngi ) = &
                 sbufensly2( 1 : u_snloc2, 1 : sbcvngi )
            sbufenlx( 1 + ( ilev - 1 ) * u_snloc2 : ilev * u_snloc2, &
                 1 : sbcvngi ) = &
                 sbufenlx2( 1 : u_snloc2, 1 : sbcvngi )
            sbufenly( 1 + ( ilev - 1 ) * u_snloc2 : ilev * u_snloc2, &
                 1 : sbcvngi ) = &
                 sbufenly2( 1 : u_snloc2, 1 : sbcvngi )
            sbufenlz( 1 + ( ilev - 1 ) * u_snloc2 : ilev * u_snloc2, &
                 1 : sbcvngi ) = &
                 sbufenlz2( 1 : u_snloc2, 1 : sbcvngi )
            u_sloclist( 1 : nface, &
                 1 + ( ilev - 1 ) * u_snloc2 : ilev * u_snloc2 ) = &
                 u_sloclist2( 1 : nface, 1 : u_snloc2 ) + &
                 ( ilev - 1 ) * u_nloc2
         end do Loop_ILEV3

      else
         call det_suf_ele_shape( scvngi, nface, &
              cvfem_on_face, &
              cv_nloc, scvfen, scvfenslx, scvfensly, scvfeweigh, &
              scvfenlx, scvfenly, scvfenlz, &
              u_nloc, sufen2, sufenslx, sufensly, &
              sufenlx, sufenly, sufenlz, &
              sbcvngi, sbcvfen, sbcvfenslx, sbcvfensly, sbcvfeweigh, &
              sbcvfenlx, sbcvfenly, sbcvfenlz,  &
              sbufen, sbufenslx, sbufensly, &
              sbufenlx, sbufenly, sbufenlz, &
              cv_sloclist, u_sloclist, cv_snloc, u_snloc, &
              ndim, cv_ele_type )
      end if Conditional_OverlappingMethod3

      ! Define the gauss points that lie on the surface of the
      ! control volume surrounding a given local node (iloc)
      ! that is FINDGPTS, COLGPTS, NCOLGPTS
      call gaussiloc( findgpts, colgpts, ncolgpts, &
           cv_neiloc, cv_nloc, scvngi )  


      ewrite(3,*) 'cv_on_face: ', cv_on_face
      ewrite(3,*) 'u_on_face: ', u_on_face
      ewrite(3,*) 'cv_sloclist: ', cv_sloclist
      ewrite(3,*) 'u_sloclist: ', u_sloclist
      ewrite(3,*) 'leaving cv_fem_shape_funs subrt, ncolgpts', ncolgpts

      if( is_overlapping ) then
         deallocate( ufen2 )
         deallocate( ufenlx2 )
         deallocate( ufenly2 )
         deallocate( ufenlz2 )
         deallocate( ufem_on_face2 )
         deallocate( sufen2 )
         deallocate( sufenslx2 )
         deallocate( sufensly2 )
         deallocate( sufenlx2 )
         deallocate( sufenly2 )
         deallocate( sufenlz2 )
         deallocate( sbufen2 )
         deallocate( sbufenslx2 )
         deallocate( sbufensly2 )
         deallocate( sbufenlx2 )
         deallocate( sbufenly2 )
         deallocate( sbufenlz2 )
         deallocate( u_sloclist2 )
      end if


      ! Set to zero anything that should be zero in case it was not pre-defined
      if( ndim < 2 ) then
         cvfenly = 0.0 ; cvfenly_short = 0.0 ; ufenly = 0.0 ; scvfenslx = 0.0 ; & 
              scvfenly = 0.0 ; sufenslx = 0.0 ; sufenly = 0.0 ; sbcvfenslx = 0.0 ;  &
              sbcvfenly = 0.0 ; sbufenslx = 0.0 ; sbufenly = 0.0 

      elseif( ndim < 3 ) then
         cvfenlz = 0.0 ; cvfenlz_short = 0.0 ; ufenlz = 0.0 ; scvfensly = 0.0 ; &
              scvfenlz = 0.0 ; sufensly = 0.0 ; sufenlz = 0.0 ; sbcvfensly = 0.0 ; &
              sbcvfenlz = 0.0 ; sbufensly = 0.0 ;sbufenlz = 0.0

      end if

      return
    end subroutine cv_fem_shape_funs


    SUBROUTINE DET_SUF_ELE_SHAPE( SCVNGI, NFACE, &  
         CVFEM_ON_FACE, &
         CV_NLOC, SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
         SCVFENLX, SCVFENLY, SCVFENLZ,  &
         U_NLOC,  SUFEN, SUFENSLX, SUFENSLY,  &
         SUFENLX, SUFENLY, SUFENLZ,  &
         SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, &
         SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
         SBUFEN, SBUFENSLX, SBUFENSLY, &
         SBUFENLX, SBUFENLY, SBUFENLZ, &
         CV_SLOCLIST, U_SLOCLIST, CV_SNLOC, U_SNLOC, &
         NDIM, CV_ELE_TYPE )
      !     
      !     - this subroutine generates the FE basis functions, weights and the
      !     - derivatives of the shape functions for a variety of elements on the 
      !     - control volume boundaries.  
      !     - The routine also generates the shape functions and derivatives 
      !     - associated with the CV surfaces and also the FV basis functions.
      !     -------------------------------
      !     - date last modified : 21/02/2012
      !     -------------------------------
 
      IMPLICIT NONE

      INTEGER, intent( in ) :: SCVNGI, CV_NLOC, U_NLOC, NFACE, &
           SBCVNGI, CV_SNLOC, U_SNLOC
      LOGICAL, DIMENSION( CV_NLOC, SCVNGI ), intent( in ) :: CVFEM_ON_FACE
      ! CV_ON_FACE(CV_KLOC,GI)=.TRUE. if CV_KLOC is on the face that GI is centred on.
      REAL, DIMENSION( CV_NLOC, SCVNGI ), intent( in ) :: SCVFEN, SCVFENSLX, SCVFENSLY, &       
           SCVFENLX, SCVFENLY, SCVFENLZ
      REAL, DIMENSION( SCVNGI ), intent( inout ) :: SCVFEWEIGH
      REAL, DIMENSION( U_NLOC, SCVNGI ), intent( in ) :: SUFEN, SUFENSLX, SUFENSLY, &       
           SUFENLX, SUFENLY, SUFENLZ
      REAL, DIMENSION( CV_SNLOC, SBCVNGI ), intent( inout ) :: SBCVFEN, SBCVFENSLX, &
           SBCVFENSLY, SBCVFENLX, SBCVFENLY, SBCVFENLZ
      REAL, DIMENSION( SBCVNGI ), intent( inout ) :: SBCVFEWEIGH
      REAL, DIMENSION( U_SNLOC, SBCVNGI ), intent( inout ) :: SBUFEN, SBUFENSLX, SBUFENSLY, &
           SBUFENLX, SBUFENLY, SBUFENLZ
      INTEGER, DIMENSION( NFACE, CV_SNLOC ), intent( inout ) ::  CV_SLOCLIST
      INTEGER, DIMENSION( NFACE, U_SNLOC ), intent( inout ) ::  U_SLOCLIST
      INTEGER, intent( in ) :: NDIM, CV_ELE_TYPE
      ! Local variables
      INTEGER :: CV_KLOC, CV_SKLOC, U_KLOC, U_SKLOC, CV_BSNGI

      ewrite(3,*) 'In DET_SUF_ELE_SHAPE'
    
      ! Obtain SBCVFEN from SCVFEN: 
      print *,'for cv:'
      CALL SCVFEN_2_SBCVFEN( CV_NLOC, CV_SNLOC, SCVNGI, SBCVNGI, &
           CV_NLOC, CV_SNLOC, CVFEM_ON_FACE, &
           SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFENLX, SBCVFENLY, SBCVFENLZ, SBCVFEWEIGH, &
           SCVFEN, SCVFENSLX, SCVFENSLY, SCVFENLX, SCVFENLY, SCVFENLZ, SCVFEWEIGH )

      print *,'U_NLOC, U_SNLOC, SCVNGI, SBCVNGI, CV_NLOC, CV_SNLOC:', &
               U_NLOC, U_SNLOC, SCVNGI, SBCVNGI, CV_NLOC, CV_SNLOC
      print *,'for u:'
      ! Obtain SBUFEN from SUFEN: 
      CALL SCVFEN_2_SBCVFEN( U_NLOC, U_SNLOC, SCVNGI, SBCVNGI, &
           CV_NLOC, CV_SNLOC, CVFEM_ON_FACE, &
           SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, SBUFENLY, SBUFENLZ, SBCVFEWEIGH, &
           SUFEN, SUFENSLX, SUFENSLY, SUFENLX, SUFENLY, SUFENLZ, SCVFEWEIGH )
       print *,'SBUFEN:',SBUFEN
       print *,'sufen:',sufen

      ! Determine CV_SLOCLIST & U_SLOCLIST
      CALL DETERMIN_SLOCLIST( CV_SLOCLIST, CV_NLOC, CV_SNLOC, SCVNGI, NFACE, &
           NDIM, CV_ELE_TYPE )
      IF( U_SNLOC == 1 ) THEN
         U_SLOCLIST( 1, 1 ) = 1
         U_SLOCLIST( 2, 1 ) = U_NLOC
      ELSE
         CALL DETERMIN_SLOCLIST( U_SLOCLIST, U_NLOC, U_SNLOC, SCVNGI, NFACE, &
              NDIM, CV_ELE_TYPE )
      ENDIF
      ewrite(3,*)'CV_SNLOC, U_SNLOC, SCVNGI:', CV_SNLOC, U_SNLOC, SCVNGI
      ewrite(3,*)'CV_SLOCLIST:', CV_SLOCLIST
      ewrite(3,*)'U_SLOCLIST:', U_SLOCLIST
       stop 2982

      RETURN
    END SUBROUTINE DET_SUF_ELE_SHAPE


    subroutine scvfen_2_sbcvfen( cv_nloc, cv_snloc, scvngi, sbcvngi, &
         cv_nloc_cells, cv_snloc_cells, cvfem_on_face, &
         sbcvfen, sbcvfenslx, sbcvfensly, sbcvfenlx, sbcvfenly, sbcvfenlz, sbcvfeweigh, &
         scvfen, scvfenslx, scvfensly, scvfenlx, scvfenly, scvfenlz, scvfeweigh )
      ! Compute SBCVFEN from SCVFEN
      implicit none
      integer, intent( in ) :: cv_nloc, cv_snloc, scvngi, sbcvngi
      integer, intent( in ) :: cv_nloc_cells, cv_snloc_cells
      logical, dimension( cv_nloc_cells, scvngi ), intent( in ) :: cvfem_on_face
      real, dimension( cv_snloc, sbcvngi ), intent( inout ) :: sbcvfen, sbcvfenslx, &
           sbcvfensly, sbcvfenlx, sbcvfenly, sbcvfenlz
      real, dimension( sbcvngi ), intent( inout ) :: sbcvfeweigh
      real, dimension( cv_nloc, scvngi ), intent( in ) ::  scvfen, &
           scvfenslx, scvfensly, scvfenlx, scvfenly, scvfenlz
      real, dimension( scvngi ), intent( in ) :: scvfeweigh
      ! Local variables
      logical, dimension( : ), allocatable :: candidate_gi,candidate_gi2
      integer :: cv_siloc, cv_iloc, cv_bsgi, cv_sgi, cv_iloc_cells
      real :: r_prodt

      ewrite(3,*) ' In scvfen_2_sbcvfen'
      ewrite(3,*) ' cv_nloc, cv_snloc, scvngi, sbcvngi, cv_nloc_cells, cv_snloc_cells:', &
           cv_nloc, cv_snloc, scvngi, sbcvngi, cv_nloc_cells, cv_snloc_cells

      allocate( candidate_gi( scvngi ) )
      allocate( candidate_gi2( scvngi ) )

      ! The CV_SNLOC surface nodes are the only nodes that are candidates. 
      ! They are the 1st nodes in the local list for the volumes
! WE CAN DELETE THE next bit probably...
!      candidate_gi = .false.
!      Loop_SGI1: do cv_sgi = 1, scvngi
!         r_prodt = 1.
!         Loop_NLOC: do cv_iloc = 1, cv_snloc
!            ewrite(3,*)'cv_iloc, cv_snloc, cv_nloc, cv_sgi, scvngi:', cv_iloc, cv_snloc, cv_nloc, cv_sgi, scvngi
!            r_prodt = r_prodt * scvfen( cv_iloc, cv_sgi )
!         end do Loop_NLOC
!         ewrite(3,*) 'cv_sgi, r_prodt:', cv_sgi, r_prodt
!         if( r_prodt > 1.e-5 ) candidate_gi( cv_sgi ) = .true.
!      end do Loop_SGI1


      do cv_sgi = 1, scvngi
         candidate_gi2( cv_sgi ) = .true.
         do cv_iloc_cells = 1, cv_snloc_cells
            if( .not.cvfem_on_face(cv_iloc_cells,cv_sgi) ) candidate_gi2( cv_sgi ) = .false.
            print *,'cv_iloc_cells,cv_sgi, cvfem_on_face(cv_iloc_cells,cv_sgi):', &
                     cv_iloc_cells,cv_sgi, cvfem_on_face(cv_iloc_cells,cv_sgi)
         end do
      end do

      print *,'candidate_gi2:',candidate_gi2

      Loop_SNLOC: do cv_siloc = 1, cv_snloc
         cv_iloc = cv_siloc
         cv_bsgi = 0
         Loop_SGI2: do cv_sgi = 1, scvngi
!            Conditional_1: if( candidate_gi( cv_sgi ) ) then
               Conditional_2: if( candidate_gi2( cv_sgi ) ) then
                  cv_bsgi = cv_bsgi + 1
                  ewrite(3,*) 'cv_siloc, cv_bsgi,cv_iloc, cv_sgi:', &
                       cv_siloc, cv_bsgi,cv_iloc, cv_sgi
                  ewrite(3,*) 'scvfen( cv_iloc, cv_sgi ):', scvfen( cv_iloc, cv_sgi )
                  sbcvfen( cv_siloc, cv_bsgi ) = scvfen( cv_iloc, cv_sgi )
                  sbcvfenslx( cv_siloc, cv_bsgi ) = scvfenslx( cv_iloc, cv_sgi )
                  sbcvfensly( cv_siloc, cv_bsgi ) = scvfensly( cv_iloc, cv_sgi )
                  sbcvfenlx( cv_siloc, cv_bsgi ) = scvfenlx( cv_iloc, cv_sgi )
                  sbcvfenly( cv_siloc, cv_bsgi ) = scvfenly( cv_iloc, cv_sgi )
                  sbcvfenlz( cv_siloc, cv_bsgi ) = scvfenlz( cv_iloc, cv_sgi )
                  sbcvfeweigh( cv_bsgi ) = scvfeweigh( cv_sgi )
               end if Conditional_2
!            end if Conditional_1
         end do Loop_SGI2
      end do Loop_SNLOC

      if(cv_bsgi.ne.sbcvngi) then
         print *,'cv_bsgi,sbcvngi:',cv_bsgi,sbcvngi
         FLAbort("cv_bsgi.ne.sbcvngi")
      endif

      deallocate( candidate_gi )
      deallocate( candidate_gi2 )
!         if(cv_nloc_cells.ne.cv_nloc) stop 29892

      return
    end subroutine scvfen_2_sbcvfen







    subroutine shapesv_fem_plus( scvngi, cv_neiloc, cv_on_face, cvfem_on_face, &
         ufem_on_face, &
         cv_ele_type, cv_nloc, scvfen, scvfenslx, scvfensly, scvfeweigh, &
         scvfenlx, scvfenly, scvfenlz, &
         u_nloc, sufen, sufenslx, sufensly, &
         sufenlx, sufenly, sufenlz, &
         ndim )  
      implicit none
      !-     
      !- This subroutine generates the FE basis functions, weights and the
      !- derivatives of the shape functions for a variety of elements on the 
      !- control volume boundaries.  
      !- The routine also generates the shape functions and derivatives 
      !- associated with the CV surfaces and also the FV basis functions.
      !-
      !- date last modified : 29/11/2011
      !-
      integer, intent( in ) :: scvngi, cv_nloc
      integer, dimension( cv_nloc, scvngi ), intent( inout ) :: cv_neiloc
      logical, dimension( cv_nloc, scvngi ), intent( inout ) :: cv_on_face, cvfem_on_face
      logical, dimension( u_nloc, scvngi ), intent( inout ) :: ufem_on_face
      integer, intent( in ) :: cv_ele_type
      real, dimension( cv_nloc, scvngi ), intent( inout ) :: scvfen, scvfenslx, scvfensly
      real, dimension( scvngi ), intent( inout ) :: scvfeweigh
      real, dimension( cv_nloc, scvngi ), intent( inout ) :: scvfenlx, scvfenly, scvfenlz
      integer, intent( in ) :: u_nloc
      real, dimension( u_nloc, scvngi ), intent( inout ) :: sufen, sufenslx, sufensly, &
           sufenlx, sufenly, sufenlz
      integer, intent( in ) :: ndim
      ! Local variables
      integer :: iloc, gi
      logical :: tri_tet
      integer, dimension( :, : ), allocatable :: cvfem_neiloc, ufem_neiloc
      real, dimension( :, : ), allocatable :: m, mu, cvn_dummy
      real, dimension( : ), allocatable :: cvweigh_dummy

      ewrite(3,*)' In ShapesV_Fem_Plus '

      ! Allocating space
      allocate( m( cv_nloc, scvngi ) )
      allocate( mu( cv_nloc, scvngi ) )
      allocate( cvn_dummy( cv_nloc, scvngi ) )
      allocate( cvweigh_dummy( scvngi ) )
      allocate( cvfem_neiloc( cv_nloc, scvngi ) )
      allocate( ufem_neiloc( u_nloc, scvngi ) )

      tri_tet=.false.

      Cond_ShapeType: Select Case( cv_ele_type )
      case( 1, 2 ) ! 1D
         call fv_1d_quad( scvngi, cv_nloc, scvfen, scvfenslx, scvfensly, scvfeweigh, &
              scvfenlx, scvfenly, scvfenlz ) ! For scalar fields
         call fv_1d_quad( scvngi, u_nloc, sufen, sufenslx, sufensly, scvfeweigh, &
              sufenlx, sufenly, sufenlz ) ! For U fields

      case( 3, 4, 7, 8 ) ! Triangle and Tetrahedra 
         tri_tet=.true.
         call suf_cv_tri_tet_shape( cv_ele_type, ndim, scvngi, cv_nloc, u_nloc, scvfeweigh, &
              scvfen, scvfenlx, scvfenly, scvfenlz, scvfenslx, scvfensly,  &
              sufen, sufenlx, sufenly, sufenlz, sufenslx, sufensly, &
              cv_neiloc, cvfem_neiloc, ufem_neiloc )
      case( 5 ) ! Bi-linear Quadrilateral
         call fvquad( scvngi, cv_nloc, scvngi, &
              m, scvfen, scvfenslx, & 
              scvfeweigh )
         call fvquad( scvngi, u_nloc, scvngi, &
              mu, sufen, sufenslx, & 
              scvfeweigh ) 
         call quad_nd_shape( ndim, cv_ele_type, scvngi, cv_nloc, u_nloc, cvn_dummy, cvweigh_dummy, &
              scvfen, scvfenlx, scvfenly, scvfenlz, &
              sufen, sufenlx, sufenly, sufenlz )

      case( 6 ) ! Tri-linear Quadrilateral
         call fvqquad( scvngi, cv_nloc, scvngi, &
              m, scvfen, scvfenslx, & 
              scvfeweigh )
         call fvqquad( scvngi, u_nloc, scvngi, &
              mu, sufen, sufenslx, & 
              scvfeweigh )
         call quad_nd_shape( ndim, cv_ele_type, scvngi, cv_nloc, u_nloc, cvn_dummy, cvweigh_dummy, &
              scvfen, scvfenlx, scvfenly, scvfenlz, &
              sufen, sufenlx, sufenly, sufenlz )

      case( 9 ) ! Tri-linear Hexahedron
         call fvhex( scvngi, cv_nloc, scvngi, &
              m, scvfen, scvfenslx, & 
              scvfensly, scvfeweigh )
         call fvhex( scvngi, u_nloc, scvngi, &
              mu, sufen, sufenslx, & 
              sufensly, scvfeweigh )
         call quad_nd_shape( ndim, cv_ele_type, scvngi, cv_nloc, u_nloc, cvn_dummy, cvweigh_dummy, &
              scvfen, scvfenlx, scvfenly, scvfenlz, &
              sufen, sufenlx, sufenly, sufenlz )

      case( 10 ) ! Tri-linear Hexahedron
         call fvqhex( scvngi, cv_nloc, scvngi, &
              m, scvfen, scvfenslx, & 
              scvfensly, scvfeweigh )
         call fvqhex( scvngi, u_nloc, scvngi, &
              mu, sufen, sufenslx, & 
              sufensly, scvfeweigh )
         call quad_nd_shape( ndim, cv_ele_type, scvngi, cv_nloc, u_nloc, cvn_dummy, cvweigh_dummy, &
              scvfen, scvfenlx, scvfenly, scvfenlz, &
              sufen, sufenlx, sufenly, sufenlz )

      case default; FLExit( "Wrong integer for CV_ELE_TYPE" )

      end Select Cond_ShapeType


      if(.not.tri_tet) then
         call volnei( cv_neiloc, cvfem_neiloc, cv_nloc, scvngi, cv_ele_type )
      end if

      ewrite(3,*)'cv_ele_type:', cv_ele_type
      ewrite(3,*)'cv_nloc, scvngi:', cv_nloc, scvngi

      cv_on_face = .false. ; cvfem_on_face = .false.
      if ( ( cv_ele_type == 1 ) .or. ( cv_ele_type == 2 ) ) then ! 1D
         do iloc = 1, cv_nloc
            cv_on_face( iloc, iloc ) = .true.
            cv_on_face( iloc, iloc + 1 ) = .true.
         end do
         cv_on_face = cvfem_on_face
      else
         do iloc = 1, cv_nloc
            do gi = 1, scvngi 
               ! ewrite(3,*)'cv_neiloc, cvfem_on_face:', iloc, gi, &
               !      cv_neiloc( iloc, gi ), cvfem_neiloc( iloc, gi ) 
               if ( cv_neiloc( iloc, gi ) == -1 ) &
                    cv_on_face( iloc, gi ) = .true.
               if ( cvfem_neiloc( iloc, gi ) == -1 ) &
                    cvfem_on_face( iloc, gi ) = .true. 
            end do
         end do

      end if

      do iloc = 1, u_nloc
         do gi = 1, scvngi 
            ufem_on_face( iloc, gi ) = ( ufem_neiloc( iloc, gi ) == -1 )
         end do
      end do

      do iloc = 1, cv_nloc
         ewrite(3,*)'iloc, cv_on_face:', iloc, ( cv_on_face( iloc, gi ), gi = 1, scvngi )
         ewrite(3,*)'iloc, cvfem_neiloc:', iloc, ( cvfem_neiloc( iloc, gi ), gi = 1, scvngi )
         ewrite(3,*)'iloc, cvfem_on_face:', iloc, ( cvfem_on_face( iloc, gi ), gi = 1, scvngi )
      end do

      deallocate( m )
      deallocate( mu )
      deallocate( cvn_dummy )
      deallocate( cvweigh_dummy )
      deallocate( cvfem_neiloc )

      return
    end subroutine shapesv_fem_plus


    SUBROUTINE FV_1D_QUAD( SCVNGI, CV_NLOC, SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
         SCVFENLX, SCVFENLY, SCVFENLZ )
      !     
      !     - this subroutine generates the FE basis functions, weights and the
      !     - derivatives of the shape functions for a variety of elements. 
      !     - The routine also generates the shape functions and derivatives 
      !     - associated with the CV surfaces and also the FV basis functions.
      !     -------------------------------
      !     - date last modified : 24/05/2003
      !     -------------------------------

      IMPLICIT NONE

      INTEGER, intent( in ) :: SCVNGI, CV_NLOC
      REAL, DIMENSION( CV_NLOC, SCVNGI ), intent( inout ) :: SCVFEN, SCVFENSLX, SCVFENSLY, &  
           SCVFENLX, SCVFENLY, SCVFENLZ
      REAL, DIMENSION( SCVNGI ), intent( inout ) :: SCVFEWEIGH
      ! Local variables
      INTEGER, PARAMETER :: TWO = 2, THREE = 3, FOUR = 4
      REAL, DIMENSION( : ), allocatable :: LX, CV_NODPOS, WEI
      INTEGER :: NCV_BOU, GPOI, ILOC
      LOGICAL :: DIFF, NDIFF, GETNDP

      !ewrite(3,*) 'cv_nloc:',cv_nloc

      ALLOCATE( LX( SCVNGI ))
      ALLOCATE( CV_NODPOS( CV_NLOC ))
      ALLOCATE( WEI( CV_NLOC ))

      SCVFEWEIGH = 1.0
      IF(SCVNGI==2) THEN
         NCV_BOU = 2
         LX(1) = -1.0
         LX(2) = +1.0
      ELSE IF(SCVNGI==3) THEN
         NCV_BOU = 3
         LX(1) = -1.0
         LX(2) =  0.0
         LX(3) = +1.0
      ELSE IF(SCVNGI==4) THEN
         NCV_BOU = 4
         LX(1) = -1.0
         LX(2) = -0.5
         LX(3) = +0.5
         LX(4) = +1.0
      ELSE IF(SCVNGI==5) THEN
         NCV_BOU = 5
         LX(1) = -1.0
         LX(2) = -1.0 + 1./3.
         LX(3) = +0.0
         LX(4) = +1.0 - 1./3.
         LX(5) = +1.0
      ELSE
         FLAbort(" Wrong number of computed surface quadrature points for CV ")
      ENDIF

      DIFF=.TRUE.
      NDIFF=.FALSE.

      GETNDP=.TRUE.
      !      Compute standard Gauss quadrature. weits and points
      CALL LAGROT(WEI,CV_NODPOS,CV_NLOC,GETNDP) 
      EWRITE(3,*)'CV_NODPOS:',CV_NODPOS

      Loop_P2: DO GPOI = 1, NCV_BOU

         Loop_ILX2: DO  ILOC = 1, CV_NLOC 
            SCVFEN( ILOC, GPOI )  = LAGRAN( NDIFF, LX( GPOI ), ILOC, CV_NLOC, CV_NODPOS )
            SCVFENSLX( ILOC, GPOI ) = 1.0
            SCVFENSLY( ILOC, GPOI ) = 0.0
            SCVFENLX( ILOC, GPOI ) = LAGRAN( DIFF, LX( GPOI ), ILOC, CV_NLOC, CV_NODPOS )
            SCVFENLY( ILOC, GPOI ) = 0.0
            SCVFENLZ( ILOC, GPOI ) = 0.0
         END DO Loop_ILX2

      end do Loop_P2

      DEALLOCATE( LX )
      DEALLOCATE( CV_NODPOS )
      DEALLOCATE( WEI )

    END SUBROUTINE FV_1D_QUAD




    SUBROUTINE SHAPESV( ELETYP, NEILOC, FEM_NEILOC,  &
         NGI, NLOC,   &
         SVNGI, &
         SVN, &
         SVNLX, SVNLY,  &
         SVWEIGH,    &
         M, D1 )
      !     
      !     - this subroutine generates the FE basis functions, weights and the
      !     - derivatives of the shape functions for a variety of elements. 
      !     - The routine also generates the shape functions and derivatives 
      !     - associated with the CV surfaces and also the FV basis functions.
      !     -------------------------------
      !     - date last modified : 24/05/2003
      !     -------------------------------

      IMPLICIT NONE

      INTEGER, intent( in ) :: ELETYP
      INTEGER, intent( in ) :: NGI, NLOC, SVNGI
      INTEGER, DIMENSION( NLOC, SVNGI ), intent( inout ) :: NEILOC, FEM_NEILOC
      REAL, DIMENSION( NLOC, SVNGI ), intent( inout ) :: SVN, SVNLX, SVNLY
      REAL, DIMENSION( SVNGI ), intent( inout ) :: SVWEIGH
      REAL, DIMENSION( NLOC, NGI ), intent( inout ) :: M
      LOGICAL, intent( in ) :: D1
      ! Local variables
      INTEGER, PARAMETER :: TWO = 2, THREE = 3
      REAL, DIMENSION( : ), allocatable :: LX, XN, DXN
      INTEGER :: P, NQUAD, ILX, NL, GPOI

      ALLOCATE( LX( TWO ))
      ALLOCATE( XN( THREE ))
      ALLOCATE( DXN( THREE ))
      !     
      !     
      Cond_ShapeType: IF(ELETYP == 1) THEN ! calculate the shape functions on the control volume boundaries

         Cond_Dimension: IF( D1 ) THEN

            SVWEIGH = 1.0
            NQUAD = 2
            LX(1) = 1./3.
            LX(2) = 2./3.

            Loop_P: DO P = 1, NQUAD

               GPOI = P
               XN( 1 ) = 0.5 * LX( P ) * ( LX( P ) - 1. )
               XN( 2 ) = 1. - LX( P ) * LX( P )
               XN( 3 ) = 0.5 * LX( P ) * ( LX( P ) + 1. )
               DXN( 1 ) = 0.5 * ( 2. * LX( P ) - 1. )
               DXN( 2 ) = - 2. * LX( P )
               DXN( 3 ) = 0.5 * ( 2. * LX( P ) + 1. )

               Loop_ILX: DO  ILX = 1, THREE 
                  NL = ILX
                  SVN( NL, GPOI )  = XN( ILX )  
                  SVNLX( NL, GPOI ) = DXN( ILX ) 
                  SVNLY( NL, GPOI ) = 0.0
               END DO Loop_ILX

            END DO Loop_P

         ELSE

            IF(ELETYP == 1) THEN

               CALL FVQUAD( NGI,    NLOC, SVNGI,&
                    M,      SVN,  SVNLX,&
                    SVWEIGH              )

            ELSE IF(ELETYP == 2) THEN

               CALL FVTRI( NGI,    NLOC, SVNGI,&
                    M,      SVN,  SVNLX,&
                    SVWEIGH              )

            ELSE IF(ELETYP == 3) THEN

               CALL FVHEX( NGI,    NLOC,    SVNGI,&
                    M,      SVN,     SVNLX, &
                    SVNLY,  SVWEIGH         )

            ELSE IF(ELETYP == 4) THEN

               CALL FVTET( NGI,    NLOC,    SVNGI,&
                    M,      SVN,     SVNLX, &
                    SVNLY,  SVWEIGH         )

            ELSE IF(ELETYP == 5) THEN

               CALL FVQQUAD(  NGI,    NLOC,    SVNGI,&
                    M,      SVN,     SVNLX, &
                    SVWEIGH                 )

            ELSE IF(ELETYP == 6) THEN

               CALL FVQHEX( NGI,   NLOC,   SVNGI,&
                    M,     SVN,    SVNLX,&
                    SVNLY, SVWEIGH        )

            ELSE IF(ELETYP == 7) THEN ! quadratic triangles

               CALL FVQTRI( NGI,    NLOC, SVNGI,&
                    M,      SVN,  SVNLX,  &
                    SVWEIGH              )

            ELSE IF(ELETYP == 8) THEN

               CALL FVQTET( NGI,   NLOC,   SVNGI,&
                    M,     SVN,    SVNLX,&
                    SVNLY, SVWEIGH        )

            END IF


         ENDIF Cond_Dimension

      ELSE

         EWRITE(3,*)'NOT GOT THIS YET'

      END IF Cond_ShapeType

      CALL VOLNEI( NEILOC, FEM_NEILOC, NLOC, SVNGI, ELETYP )

      DEALLOCATE( LX )
      DEALLOCATE( XN )
      DEALLOCATE( DXN )

      RETURN
    END SUBROUTINE SHAPESV



    SUBROUTINE FVQUAD( NGI,    NLOC,    SVNGI,&
         M,      SVN,     SVNLX, &
         SVWEIGH                 )
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
    END SUBROUTINE FVQUAD
    !     
    !     
    !     
    !     
    SUBROUTINE FVTRI( NGI,    NLOC, SVNGI,&
         M,      SVN,  SVNLX,  &
         SVWEIGH              )
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
                          + CORN(IFACE,ICOORD,JLOC)&
                          *0.5*(1.+XI(JLOC)*XIGP(GJ))
                     !     
                     DPDXI(ICOORD) = DPDXI(ICOORD)&
                          + CORN(IFACE,ICOORD,JLOC)&
                          *0.5*XI(JLOC)
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
    END SUBROUTINE FVTRI

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
    END SUBROUTINE FVHEX
    !     
    SUBROUTINE FVTET( NGI,   NLOC, SVNGI,&
         M,     SVN,  SVNLX,&
         SVNLY, SVWEIGH      )
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
      !      use fldebug
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
      !      ewrite(2,*) 'SUBROUTINE FVTET()'

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
    END SUBROUTINE FVTET
    !     
    !     
    SUBROUTINE FVQQUAD( NGI,    NLOC,    SVNGI,&
                                !     - REALS
         M,      SVN,     SVNLX, &
         SVWEIGH                 )
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
    END SUBROUTINE FVQQUAD

    SUBROUTINE FVQHEX( NGI,   NLOC,   SVNGI,&
                                !     - REALS            
         M,     SVN,    SVNLX,&
         SVNLY, SVWEIGH        )
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
    END SUBROUTINE FVQHEX

    SUBROUTINE FVQTRI( NGI,    NLOC, SVNGI,&
         M,      SVN,  SVNLX,  &
         SVWEIGH              )
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
         IF((ILOC == 1).OR.(ILOC == 3).OR.(ILOC == 5)) THEN
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
         IF((ILOC == 1).OR.(ILOC == 3).OR.(ILOC == 5)) THEN
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
    END SUBROUTINE FVQTRI

    SUBROUTINE FVQTET( NGI,   NLOC,   SVNGI,&
                                !     - REALS            
         M,     SVN,    SVNLX,&
         SVNLY, SVWEIGH        )
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
         IF((ILOC == 1).OR.(ILOC == 3).OR.&
              &        (ILOC == 5).OR.(ILOC == 10)) THEN
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
         IF((ILOC == 1).OR.(ILOC == 3).OR.&
              &        (ILOC == 5).OR.(ILOC == 10)) THEN
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
    END SUBROUTINE FVQTET



!!!==============================================================

    SUBROUTINE VOLNEI( NEILOC, FEM_NEILOC, NLOC, SVNGI, CV_ELE_TYPE )
      !--------------------------------------------------------   
      !- this subroutine calculates NEILOC which is the       
      !- array containing information given a local node
      !- and an integration point what is the other opposing 
      !- local node. It contains -1 if on the boundary of
      !- of the element.  
      !--------------------------------------------------------
      !- date last modified : 25/11/2011
      !--------------------------------------------------------

      IMPLICIT NONE     
      INTEGER, intent( in ) :: NLOC, SVNGI, CV_ELE_TYPE   
      INTEGER, DIMENSION( NLOC, SVNGI ), intent( inout ) :: NEILOC
      INTEGER, DIMENSION( NLOC, SVNGI ), intent( inout ) :: FEM_NEILOC
      ! Local variables    
      INTEGER :: ILOC, IEXT, IGP

      NEILOC = 0
      FEM_NEILOC = 0

      Conditional_Type: SELECT CASE( CV_ELE_TYPE )

      CASE( 1, 2 ) ! 1D
         DO ILOC = 1, NLOC
            NEILOC( ILOC, ILOC ) = ILOC - 1
            NEILOC( ILOC, ILOC + 1 ) = ILOC + 1
         END DO
         NEILOC( 1, 1 ) = -1
         NEILOC( NLOC, SVNGI ) = -1
         FEM_NEILOC = NEILOC

      CASE( 3 ) ! Linear Triangle
         FLAbort( " Defined elsewhere -- it needs to be updated " )

         ILOC = 1 ! CV 1
         NEILOC( ILOC, 1 ) = 2
         NEILOC( ILOC, 3 ) = 3
         ! External
         neiloc( iloc, 4 ) = -1
         neiloc( iloc, 9 ) = -1

         ILOC = 2 ! CV 2
         NEILOC( ILOC, 1 ) = 1
         NEILOC( ILOC, 2 ) = 3
         ! External
         neiloc( iloc, 5 ) = -1
         neiloc( iloc, 6 ) = -1

         ILOC = 3 ! CV 3
         NEILOC( ILOC, 2) = 2
         NEILOC( ILOC, 3) = 1
         ! External
         neiloc( iloc, 7 ) = -1
         neiloc( iloc, 8 ) = -1

         fem_neiloc = neiloc

      CASE( 4 ) ! Quadratic Triangle
         FLAbort( " Defined elsewhere -- it needs to be updated " )
         ILOC = 1
         ! Face 1
         NEILOC( ILOC, 1 ) = 2
         NEILOC( ILOC, 2 ) = 2
         ! Face 2
         NEILOC( ILOC,3 ) = 6
         NEILOC( ILOC,4 ) = 6
         ! External
         neiloc( iloc, 19 ) = -1
         neiloc( iloc, 20 ) = -1
         neiloc( iloc, 35 ) = -1
         neiloc( iloc, 36 ) = -1

         ILOC = 2
         ! Face 1
         NEILOC( ILOC,1 ) = 1
         NEILOC( ILOC,2 ) = 1
         ! Face 7
         NEILOC( ILOC,13 ) = 6
         NEILOC( ILOC,14 ) = 6
         ! Face 8
         NEILOC( ILOC,15 ) = 4
         NEILOC( ILOC,16 ) = 4
         ! Face 3
         NEILOC( ILOC,5 ) = 3
         NEILOC( ILOC,6 ) = 3
         ! External
         neiloc( iloc, 21 ) = -1
         neiloc( iloc, 22 ) = -1

         ILOC = 3
         ! Face 3
         NEILOC( ILOC,5 ) = 2
         NEILOC( ILOC,6 ) = 2
         ! Face 4
         NEILOC( ILOC,7 ) = 4
         NEILOC( ILOC,8 ) = 4
         ! External
         neiloc( iloc, 23 ) = -1
         neiloc( iloc, 24 ) = -1
         neiloc( iloc, 25 ) = -1
         neiloc( iloc, 26 ) = -1

         ILOC = 4
         ! Face 4
         NEILOC( ILOC,7 ) = 3
         NEILOC( ILOC,8 ) = 3
         ! Face 8
         NEILOC( ILOC,15 ) = 2
         NEILOC( ILOC,16 ) = 2
         ! Face 9
         NEILOC( ILOC,17 ) = 6
         NEILOC( ILOC,18 ) = 6
         ! Face 5
         NEILOC( ILOC,9 )  = 5
         NEILOC( ILOC,10 ) = 5
         ! External
         neiloc( iloc, 27 ) = -1
         neiloc( iloc, 28 ) = -1

         ILOC = 5
         ! Face 5
         NEILOC( ILOC,9 )  = 4  
         NEILOC( ILOC,10 ) = 4
         ! Face 6
         NEILOC( ILOC,11 ) = 6 
         NEILOC( ILOC,12 ) = 6
         ! External
         neiloc( iloc, 29 ) = -1
         neiloc( iloc, 30 ) = -1
         neiloc( iloc, 31 ) = -1
         neiloc( iloc, 32 ) = -1

         ILOC = 6
         ! Face 6
         NEILOC( ILOC,11 ) = 5
         NEILOC( ILOC,12 ) = 5
         ! Face 9
         NEILOC( ILOC,17 ) = 4
         NEILOC( ILOC,18 ) = 4
         ! Face 7
         NEILOC( ILOC,13 ) = 2 
         NEILOC( ILOC,14 ) = 2
         ! Face 2
         NEILOC( ILOC,3 ) = 1
         NEILOC( ILOC,4 ) = 1
         ! External
         neiloc( iloc, 35 ) = -1
         neiloc( iloc, 36 ) = -1

      CASE( 5  ) ! Bi-linear Quadrilateral
         ILOC = 1
         NEILOC( ILOC, 1 ) = 3
         NEILOC( ILOC, 3 ) = 2
         ! External
         neiloc( iloc, 5 ) = -1
         neiloc( iloc, 6 ) = -1

         ILOC = 2
         NEILOC( ILOC, 2 ) = 4
         NEILOC( ILOC, 3 ) = 1
         ! External
         neiloc( iloc, 7 ) = -1
         neiloc( iloc, 8 ) = -1

         ILOC = 3
         NEILOC( ILOC, 1 ) = 1
         NEILOC( ILOC, 4 ) = 4
         ! External
         neiloc( iloc, 11 ) = -1
         neiloc( iloc, 12 ) = -1

         ILOC = 4
         NEILOC( ILOC, 2 ) = 2
         NEILOC( ILOC, 4 ) = 3
         ! External
         neiloc( iloc, 9 ) = -1
         neiloc( iloc, 10 ) = -1

      CASE( 6 ) ! Bi-quadratic Quadrilateral
         ILOC = 1
         NEILOC( ILOC, 1 ) = 2
         NEILOC( ILOC, 2 ) = 2
         NEILOC( ILOC, 5 ) = 4
         NEILOC( ILOC, 6 ) = 4
         ! External
         neiloc( iloc, 25 ) = -1
         neiloc( iloc, 26 ) = -1
         neiloc( iloc, 27 ) = -1
         neiloc( iloc, 28 ) = -1

         ILOC = 2
         NEILOC( ILOC, 1 ) = 1
         NEILOC( ILOC, 2 ) = 1
         NEILOC( ILOC, 7 ) = 5
         NEILOC( ILOC, 8 ) = 5
         NEILOC( ILOC, 3 ) = 3
         NEILOC( ILOC, 4 ) = 3 
         ! External
         neiloc( iloc, 47 ) = -1
         neiloc( iloc, 48 ) = -1     

         ILOC = 3
         NEILOC( ILOC, 3) = 2
         NEILOC( ILOC, 4) = 2
         NEILOC( ILOC, 9) = 6
         NEILOC( ILOC, 10 ) = 6
         ! External
         neiloc( iloc, 43 ) = -1
         neiloc( iloc, 44 ) = -1
         neiloc( iloc, 45 ) = -1
         neiloc( iloc, 46 ) = -1

         ILOC = 4
         NEILOC( ILOC, 5) = 1
         NEILOC( ILOC, 6) = 1
         NEILOC( ILOC, 11 ) = 5
         NEILOC( ILOC, 12 ) = 5
         NEILOC( ILOC, 15 ) = 7
         NEILOC( ILOC, 16 ) = 7
         ! External
         neiloc( iloc, 29 ) = -1
         neiloc( iloc, 30 ) = -1   

         ILOC = 5
         NEILOC( ILOC, 7) = 2
         NEILOC( ILOC, 8) = 2
         NEILOC( ILOC, 13 ) = 6
         NEILOC( ILOC, 14 ) = 6
         NEILOC( ILOC, 18 ) = 8
         NEILOC( ILOC, 17 ) = 8
         NEILOC( ILOC, 12 ) = 4
         NEILOC( ILOC, 11 ) = 4

         ILOC = 6
         NEILOC( ILOC, 9) = 3
         NEILOC( ILOC, 10 ) = 3
         NEILOC( ILOC, 13 ) = 5
         NEILOC( ILOC, 14 ) = 5
         NEILOC( ILOC, 19 ) = 9
         NEILOC( ILOC, 20 ) = 9
         ! External
         neiloc( iloc, 41 ) = -1
         neiloc( iloc, 42 ) = -1 

         ILOC = 7
         NEILOC( ILOC, 15 ) = 4
         NEILOC( ILOC, 16 ) = 4
         NEILOC( ILOC, 21 ) = 8
         NEILOC( ILOC, 22 ) = 8
         ! External
         neiloc( iloc, 31 ) = -1
         neiloc( iloc, 32 ) = -1
         neiloc( iloc, 33 ) = -1
         neiloc( iloc, 34 ) = -1

         ILOC = 8
         NEILOC( ILOC, 21 ) = 7
         NEILOC( ILOC, 22 ) = 7
         NEILOC( ILOC, 17 ) = 5
         NEILOC( ILOC, 18 ) = 5
         NEILOC( ILOC, 23 ) = 9
         NEILOC( ILOC, 24 ) = 9
         ! External
         neiloc( iloc, 35 ) = -1
         neiloc( iloc, 36 ) = -1 

         ILOC = 9
         NEILOC( ILOC, 23 ) = 8
         NEILOC( ILOC, 24 ) = 8
         NEILOC( ILOC, 19 ) = 6
         NEILOC( ILOC, 20 ) = 6
         ! External
         neiloc( iloc, 37 ) = -1
         neiloc( iloc, 38 ) = -1
         neiloc( iloc, 39 ) = -1
         neiloc( iloc, 40 ) = -1

      CASE( 7 ) ! Linear Tetrahedron
         iext = 6
         igp = 1
         ILOC = 1
         NEILOC( ILOC, 1 ) = 2
         NEILOC( ILOC, 3 ) = 3
         NEILOC( ILOC, 4 ) = 4
         ! External
         neiloc( iloc, iext + 1 : iext + 3 * igp ) = -1

         ILOC = 2
         NEILOC( ILOC, 1 ) = 1
         NEILOC( ILOC, 2 ) = 3
         NEILOC( ILOC, 5 ) = 4
         ! External
         neiloc( iloc, iext + 3 * igp + 1 : iext + 6 * igp ) = -1

         ILOC = 3
         NEILOC( ILOC, 2 ) = 2
         NEILOC( ILOC, 3 ) = 1
         NEILOC( ILOC, 6 ) = 4
         ! External
         neiloc( iloc, iext + 6 * igp + 1 : iext + 9 * igp ) = -1

         ILOC = 4
         NEILOC( ILOC, 4 ) = 1
         NEILOC( ILOC, 5 ) = 2
         NEILOC( ILOC, 6 ) = 3
         ! External
         neiloc( iloc, iext + 9 * igp + 1 : iext + 12 * igp ) = -1

      CASE( 8 ) ! Quadratic Tetrahedron
         iext = 96
         igp = 4
         ILOC = 1
         ! Face 1
         NEILOC( ILOC, 1 ) = 6
         NEILOC( ILOC, 2 ) = 6
         NEILOC( ILOC, 3 ) = 6
         NEILOC( ILOC, 4 ) = 6
         ! Face 2
         NEILOC( ILOC, 5 ) = 2
         NEILOC( ILOC, 6 ) = 2
         NEILOC( ILOC, 7 ) = 2
         NEILOC( ILOC, 8 ) = 2
         ! Face 3
         NEILOC( ILOC, 9 )  = 7
         NEILOC( ILOC, 10 ) = 7
         NEILOC( ILOC, 11 ) = 7
         NEILOC( ILOC, 12 ) = 7
         ! External
         neiloc( iloc, iext + 1 : iext + 3 * igp ) = -1

         ILOC = 2
         ! Face 2
         NEILOC( ILOC, 5 ) = 1
         NEILOC( ILOC, 6 ) = 1
         NEILOC( ILOC, 7 ) = 1
         NEILOC( ILOC, 8 ) = 1
         ! Face 5
         NEILOC( ILOC, 17 ) = 3
         NEILOC( ILOC, 18 ) = 3
         NEILOC( ILOC, 19 ) = 3
         NEILOC( ILOC, 20 ) = 3
         ! Face 20
         NEILOC( ILOC, 77 ) = 8
         NEILOC( ILOC, 78 ) = 8
         NEILOC( ILOC, 79 ) = 8
         NEILOC( ILOC, 80 ) = 8
         ! Face 21
         NEILOC( ILOC, 81 ) = 7
         NEILOC( ILOC, 82 ) = 7
         NEILOC( ILOC, 83 ) = 7
         NEILOC( ILOC, 84 ) = 7
         ! Face 23
         NEILOC( ILOC, 89 ) = 6
         NEILOC( ILOC, 80 ) = 6
         NEILOC( ILOC, 91 ) = 6
         NEILOC( ILOC, 92 ) = 6
         ! Face 24
         NEILOC( ILOC, 93 ) = 4
         NEILOC( ILOC, 94 ) = 4
         NEILOC( ILOC, 95 ) = 4
         NEILOC( ILOC, 96 ) = 4
         ! External
         neiloc( iloc, iext + 3 * igp + 1 : iext + 5 * igp ) = -1

         ILOC = 3
         ! Face 4
         NEILOC( ILOC, 13 ) = 4
         NEILOC( ILOC, 14 ) = 4
         NEILOC( ILOC, 15 ) = 4
         NEILOC( ILOC, 16 ) = 4
         ! Face 5
         NEILOC( ILOC, 17 ) = 2
         NEILOC( ILOC, 18 ) = 2
         NEILOC( ILOC, 19 ) = 2
         NEILOC( ILOC, 20 ) = 2
         ! Face 6
         NEILOC( ILOC, 21 ) = 8
         NEILOC( ILOC, 22 ) = 8
         NEILOC( ILOC, 23 ) = 8
         NEILOC( ILOC, 24 ) = 8
         ! External
         neiloc( iloc, iext + 5 * igp + 1 : iext + 8 * igp ) = -1

         ILOC = 4
         ! Face 4
         NEILOC( ILOC, 13 ) = 3
         NEILOC( ILOC, 14 ) = 3
         NEILOC( ILOC, 15 ) = 3
         NEILOC( ILOC, 16 ) = 3
         ! Face 8
         NEILOC( ILOC, 29 ) = 5
         NEILOC( ILOC, 30 ) = 5
         NEILOC( ILOC, 31 ) = 5
         NEILOC( ILOC, 32 ) = 5
         ! Face 18
         NEILOC( ILOC, 69 ) = 9
         NEILOC( ILOC, 70 ) = 9
         NEILOC( ILOC, 71 ) = 9
         NEILOC( ILOC, 72 ) = 9
         ! Face 19
         NEILOC( ILOC, 73 ) = 8
         NEILOC( ILOC, 74 ) = 8
         NEILOC( ILOC, 75 ) = 8
         NEILOC( ILOC, 76 ) = 8
         ! Face 22
         NEILOC( ILOC, 85 ) = 6
         NEILOC( ILOC, 86 ) = 6
         NEILOC( ILOC, 87 ) = 6
         NEILOC( ILOC, 88 ) = 6
         ! Face 24
         NEILOC( ILOC, 93 ) = 2
         NEILOC( ILOC, 94 ) = 2
         NEILOC( ILOC, 95 ) = 2
         NEILOC( ILOC, 96 ) = 2
         ! External
         neiloc( iloc, iext + 8 * igp + 1 : iext + 10 * igp ) = -1

         ILOC = 5
         ! Face 7
         NEILOC( ILOC, 25 ) = 6
         NEILOC( ILOC, 26 ) = 6
         NEILOC( ILOC, 27 ) = 6
         NEILOC( ILOC, 28 ) = 6
         ! Face 8
         NEILOC( ILOC, 29 ) = 4
         NEILOC( ILOC, 30 ) = 4
         NEILOC( ILOC, 31 ) = 4
         NEILOC( ILOC, 32 ) = 4
         ! Face 9
         NEILOC( ILOC, 33 ) = 9
         NEILOC( ILOC, 34 ) = 9
         NEILOC( ILOC, 35 ) = 9
         NEILOC( ILOC, 36 ) = 9
         ! External
         neiloc( iloc, iext + 10 * igp + 1 : iext + 13 * igp ) = -1

         ILOC = 6
         ! Face 1
         NEILOC( ILOC, 1 ) = 1
         NEILOC( ILOC, 2 ) = 1
         NEILOC( ILOC, 3 ) = 1
         NEILOC( ILOC, 4 ) = 1
         ! Face 7
         NEILOC( ILOC, 25 ) = 5
         NEILOC( ILOC, 26 ) = 5
         NEILOC( ILOC, 27 ) = 5
         NEILOC( ILOC, 28 ) = 5
         ! Face 16
         NEILOC( ILOC, 61 ) = 7
         NEILOC( ILOC, 62 ) = 7
         NEILOC( ILOC, 63 ) = 7
         NEILOC( ILOC, 64 ) = 7
         ! Face 17
         NEILOC( ILOC, 65 ) = 9
         NEILOC( ILOC, 66 ) = 9
         NEILOC( ILOC, 67 ) = 9
         NEILOC( ILOC, 68 ) = 9
         ! Face 22
         NEILOC( ILOC, 85 ) = 4
         NEILOC( ILOC, 86 ) = 4
         NEILOC( ILOC, 87 ) = 4
         NEILOC( ILOC, 88 ) = 4
         ! Face 23
         NEILOC( ILOC, 89 ) = 2
         NEILOC( ILOC, 90 ) = 2
         NEILOC( ILOC, 91 ) = 2
         NEILOC( ILOC, 92 ) = 2
         ! External
         neiloc( iloc, iext + 13 * igp + 1 : iext + 15 * igp ) = -1
         !     
         ILOC = 7
         ! Face 3
         NEILOC( ILOC, 9 )  = 1
         NEILOC( ILOC, 10 ) = 1
         NEILOC( ILOC, 11 ) = 1
         NEILOC( ILOC, 12 ) = 1
         ! Face 11
         NEILOC( ILOC, 41 ) = 10
         NEILOC( ILOC, 42 ) = 10
         NEILOC( ILOC, 43 ) = 10
         NEILOC( ILOC, 44 ) = 10
         ! Face 13
         NEILOC( ILOC, 49 ) = 9
         NEILOC( ILOC, 50 ) = 9
         NEILOC( ILOC, 51 ) = 9
         NEILOC( ILOC, 52 ) = 9
         ! Face 15
         NEILOC( ILOC, 57 ) = 8
         NEILOC( ILOC, 58 ) = 8
         NEILOC( ILOC, 59 ) = 8
         NEILOC( ILOC, 60 ) = 8
         ! Face 16
         NEILOC( ILOC, 61 ) = 6
         NEILOC( ILOC, 62 ) = 6
         NEILOC( ILOC, 63 ) = 6
         NEILOC( ILOC, 64 ) = 6
         ! Face 21
         NEILOC( ILOC, 81 ) = 2
         NEILOC( ILOC, 82 ) = 2
         NEILOC( ILOC, 83 ) = 2
         NEILOC( ILOC, 84 ) = 2
         ! External
         neiloc( iloc, iext + 15 * igp + 1 : iext + 17 * igp ) = -1

         ILOC = 8
         ! Face 6
         NEILOC( ILOC, 21 ) = 3
         NEILOC( ILOC, 22 ) = 3
         NEILOC( ILOC, 23 ) = 3
         NEILOC( ILOC, 24 ) = 3
         ! Face 12
         NEILOC( ILOC, 45 ) = 10
         NEILOC( ILOC, 46 ) = 10
         NEILOC( ILOC, 47 ) = 10
         NEILOC( ILOC, 48 ) = 10
         ! Face 14
         NEILOC( ILOC, 53 ) = 9
         NEILOC( ILOC, 54 ) = 9
         NEILOC( ILOC, 55 ) = 9
         NEILOC( ILOC, 56 ) = 9
         ! Face 15
         NEILOC( ILOC, 57 ) = 7
         NEILOC( ILOC, 58 ) = 7
         NEILOC( ILOC, 59 ) = 7
         NEILOC( ILOC, 60 ) = 7
         ! Face 19
         NEILOC( ILOC, 73 ) = 4
         NEILOC( ILOC, 74 ) = 4
         NEILOC( ILOC, 75 ) = 4
         NEILOC( ILOC, 76 ) = 4
         ! Face 20
         NEILOC( ILOC, 77 ) = 2
         NEILOC( ILOC, 78 ) = 2
         NEILOC( ILOC, 79 ) = 2
         NEILOC( ILOC, 80 ) = 2
         ! External
         neiloc( iloc, iext + 17 * igp + 1 : iext + 20 * igp ) = -1

         ILOC = 9
         ! Face 9
         NEILOC( ILOC, 33 ) = 5
         NEILOC( ILOC, 34 ) = 5
         NEILOC( ILOC, 35 ) = 5
         NEILOC( ILOC, 36 ) = 5
         ! Face 10
         NEILOC( ILOC, 37 ) = 10
         NEILOC( ILOC, 38 ) = 10
         NEILOC( ILOC, 39 ) = 10
         NEILOC( ILOC, 40 ) = 10
         ! Face 13
         NEILOC( ILOC, 49 ) = 7
         NEILOC( ILOC, 50 ) = 7
         NEILOC( ILOC, 51 ) = 7
         NEILOC( ILOC, 52 ) = 7
         ! Face 14
         NEILOC( ILOC, 53 ) = 8
         NEILOC( ILOC, 54 ) = 8
         NEILOC( ILOC, 55 ) = 8
         NEILOC( ILOC, 56 ) = 8
         ! Face 17
         NEILOC( ILOC, 65 ) = 6
         NEILOC( ILOC, 66 ) = 6
         NEILOC( ILOC, 67 ) = 6
         NEILOC( ILOC, 68 ) = 6
         ! Face 18
         NEILOC( ILOC, 69 ) = 4
         NEILOC( ILOC, 70 ) = 4
         NEILOC( ILOC, 71 ) = 4
         NEILOC( ILOC, 72 ) = 4
         ! External
         neiloc( iloc, iext + 20 * igp + 1 : iext + 22 * igp ) = -1

         ILOC = 10
         ! Face 10
         NEILOC( ILOC, 37 ) = 9
         NEILOC( ILOC, 38 ) = 9
         NEILOC( ILOC, 39 ) = 9
         NEILOC( ILOC, 40 ) = 9
         ! Face 11
         NEILOC( ILOC, 41 ) = 7
         NEILOC( ILOC, 42 ) = 7
         NEILOC( ILOC, 43 ) = 7
         NEILOC( ILOC, 44 ) = 7
         ! Face 12
         NEILOC( ILOC, 45 ) = 8
         NEILOC( ILOC, 46 ) = 8
         NEILOC( ILOC, 47 ) = 8
         NEILOC( ILOC, 48 ) = 8
         ! External
         neiloc( iloc, iext + 22 * igp + 1 : iext + 24 * igp ) = -1

      CASE( 9 ) ! Tri-linear Hexahedron
         ILOC = 1
         NEILOC( ILOC, 1 ) = 3
         NEILOC( ILOC, 3 ) = 2
         NEILOC( ILOC, 5 ) = 5
         ! External
         neiloc( iloc, 13 ) = -1
         neiloc( iloc, 17 ) = -1
         neiloc( iloc, 21 ) = -1

         ILOC = 2
         NEILOC( ILOC, 2 ) = 4
         NEILOC( ILOC, 3 ) = 1
         NEILOC( ILOC, 6 ) = 6
         ! External
         neiloc( iloc, 14 ) = -1
         neiloc( iloc, 18 ) = -1
         neiloc( iloc, 25 ) = -1

         ILOC = 3
         NEILOC( ILOC, 1 ) = 1
         NEILOC( ILOC, 4 ) = 4
         NEILOC( ILOC, 7 ) = 7
         ! External
         neiloc( iloc, 16 ) = -1
         neiloc( iloc, 22 ) = -1
         neiloc( iloc, 31 ) = -1

         ILOC = 4
         NEILOC( ILOC, 2 ) = 2
         NEILOC( ILOC, 4 ) = 3
         NEILOC( ILOC, 8 ) = 8
         ! External
         neiloc( iloc, 15 ) = -1
         neiloc( iloc, 26 ) = -1
         neiloc( iloc, 30 ) = -1

         ILOC = 5
         NEILOC( ILOC, 5 ) = 1
         NEILOC( ILOC, 9 ) = 7
         NEILOC( ILOC, 11 ) = 6
         ! External
         neiloc( iloc, 19 ) = -1
         neiloc( iloc, 24 ) = -1
         neiloc( iloc, 33 ) = -1

         ILOC = 6
         NEILOC( ILOC, 6 ) = 2
         NEILOC( ILOC, 10 ) = 8
         NEILOC( ILOC, 11 ) = 5
         ! External
         neiloc( iloc, 20 ) = -1
         neiloc( iloc, 27 ) = -1
         neiloc( iloc, 34 ) = -1

         ILOC = 7
         NEILOC( ILOC, 7 ) = 3
         NEILOC( ILOC, 9 ) = 5
         NEILOC( ILOC, 12 ) = 8
         ! External
         neiloc( iloc, 23 ) = -1
         neiloc( iloc, 32 ) = -1
         neiloc( iloc, 36 ) = -1

         ILOC = 8
         NEILOC( ILOC, 8 ) = 4
         NEILOC( ILOC, 10 ) = 6
         NEILOC( ILOC, 12 ) = 7
         ! External
         neiloc( iloc, 28 ) = -1
         neiloc( iloc, 29 ) = -1
         neiloc( iloc, 35 ) = -1

      CASE( 10 ) ! Tri-quadratic Hexahedron
         iext = 216
         igp = 4
         ILOC = 1
         ! Face 1
         NEILOC( ILOC, 1 ) = 2
         NEILOC( ILOC, 2 ) = 2
         NEILOC( ILOC, 3 ) = 2
         NEILOC( ILOC, 4 ) = 2
         ! Face 3
         NEILOC( ILOC, 9 ) = 4
         NEILOC( ILOC, 10 ) = 4
         NEILOC( ILOC, 11 ) = 4
         NEILOC( ILOC, 12 ) = 4
         ! Face 13
         NEILOC( ILOC, 49 ) = 10
         NEILOC( ILOC, 50 ) = 10
         NEILOC( ILOC, 51 ) = 10
         NEILOC( ILOC, 52 ) = 10
         ! External
         neiloc( iloc, iext + 1 : iext + 3 * igp ) = -1

         ILOC = 2
         ! Face 1
         NEILOC( ILOC, 1 ) = 1
         NEILOC( ILOC, 2 ) = 1
         NEILOC( ILOC, 3 ) = 1
         NEILOC( ILOC, 4 ) = 1
         ! Face 4
         NEILOC( ILOC, 13 ) = 5
         NEILOC( ILOC, 14 ) = 5
         NEILOC( ILOC, 15 ) = 5
         NEILOC( ILOC, 16 ) = 5
         ! Face 2
         NEILOC( ILOC, 5 ) = 3
         NEILOC( ILOC, 6 ) = 3
         NEILOC( ILOC, 7 ) = 3
         NEILOC( ILOC, 8 ) = 3
         ! Face 14
         NEILOC( ILOC, 53 ) = 11
         NEILOC( ILOC, 54 ) = 11
         NEILOC( ILOC, 55 ) = 11
         NEILOC( ILOC, 56 ) = 11
         ! External
         neiloc( iloc, iext + 3 * igp + 1 : iext + 5 * igp ) = -1

         ILOC = 3
         ! Face 2
         NEILOC( ILOC, 5 ) = 2
         NEILOC( ILOC, 6 ) = 2
         NEILOC( ILOC, 7 ) = 2
         NEILOC( ILOC, 8 ) = 2
         ! Face 5
         NEILOC( ILOC, 17 ) = 6
         NEILOC( ILOC, 18 ) = 6
         NEILOC( ILOC, 19 ) = 6
         NEILOC( ILOC, 20 ) = 6
         ! Face 15
         NEILOC( ILOC, 57 ) = 12
         NEILOC( ILOC, 58 ) = 12
         NEILOC( ILOC, 59 ) = 12
         NEILOC( ILOC, 60 ) = 12
         ! External
         neiloc( iloc, iext + 5 * igp + 1 : iext + 8 * igp ) = -1

         ILOC = 4
         ! Face 3
         NEILOC( ILOC, 9 ) = 1
         NEILOC( ILOC, 10 ) = 1
         NEILOC( ILOC, 11 ) = 1
         NEILOC( ILOC, 12 ) = 1
         ! Face 6
         NEILOC( ILOC, 21 ) = 5
         NEILOC( ILOC, 22 ) = 5
         NEILOC( ILOC, 23 ) = 5
         NEILOC( ILOC, 24 ) = 5
         ! Face 8
         NEILOC( ILOC, 29 ) = 7
         NEILOC( ILOC, 30 ) = 7
         NEILOC( ILOC, 31 ) = 7
         NEILOC( ILOC, 32 ) = 7
         ! Face 16
         NEILOC( ILOC, 61 ) = 13
         NEILOC( ILOC, 62 ) = 13
         NEILOC( ILOC, 63 ) = 13
         NEILOC( ILOC, 64 ) = 13
         ! External
         neiloc( iloc, iext + 8 * igp + 1 : iext + 10 * igp ) = -1

         ILOC = 5
         ! Face 4
         NEILOC( ILOC, 13 ) = 2
         NEILOC( ILOC, 14 ) = 2
         NEILOC( ILOC, 15 ) = 2
         NEILOC( ILOC, 16 ) = 2
         ! Face 7
         NEILOC( ILOC, 25 ) = 6
         NEILOC( ILOC, 26 ) = 6
         NEILOC( ILOC, 27 ) = 6
         NEILOC( ILOC, 28 ) = 6
         ! Face 9
         NEILOC( ILOC, 33 ) = 8
         NEILOC( ILOC, 34 ) = 8
         NEILOC( ILOC, 35 ) = 8
         NEILOC( ILOC, 36 ) = 8
         ! Face 6
         NEILOC( ILOC, 21 ) = 4
         NEILOC( ILOC, 22 ) = 4
         NEILOC( ILOC, 23 ) = 4
         NEILOC( ILOC, 24 ) = 4
         ! Face 17
         NEILOC( ILOC, 65 ) = 14
         NEILOC( ILOC, 66 ) = 14
         NEILOC( ILOC, 67 ) = 14
         NEILOC( ILOC, 68 ) = 14
         ! External
         neiloc( iloc, iext + 10 * igp + 1 : iext + 11 * igp ) = -1

         ILOC = 6
         ! Face 5
         NEILOC( ILOC, 17 ) = 3
         NEILOC( ILOC, 18 ) = 3
         NEILOC( ILOC, 19 ) = 3
         NEILOC( ILOC, 20 ) = 3
         ! Face 7
         NEILOC( ILOC, 25 ) = 5
         NEILOC( ILOC, 26 ) = 5
         NEILOC( ILOC, 27 ) = 5
         NEILOC( ILOC, 28 ) = 5
         ! Face 10
         NEILOC( ILOC, 37 ) = 9
         NEILOC( ILOC, 38 ) = 9
         NEILOC( ILOC, 39 ) = 9
         NEILOC( ILOC, 40 ) = 9
         ! Face 18
         NEILOC( ILOC, 69 ) = 15
         NEILOC( ILOC, 70 ) = 15
         NEILOC( ILOC, 71 ) = 15
         NEILOC( ILOC, 72 ) = 15
         ! External
         neiloc( iloc, iext + 11 * igp + 1 : iext + 13 * igp ) = -1

         ILOC = 7
         ! Face 8
         NEILOC( ILOC, 29 ) = 4
         NEILOC( ILOC, 30 ) = 4
         NEILOC( ILOC, 31 ) = 4
         NEILOC( ILOC, 32 ) = 4
         ! Face 11
         NEILOC( ILOC, 41 ) = 8
         NEILOC( ILOC, 42 ) = 8
         NEILOC( ILOC, 43 ) = 8
         NEILOC( ILOC, 44 ) = 8
         ! Face 19
         NEILOC( ILOC, 73 ) = 16
         NEILOC( ILOC, 74 ) = 16
         NEILOC( ILOC, 75 ) = 16
         NEILOC( ILOC, 76 ) = 16
         ! External
         neiloc( iloc, iext + 13 * igp + 1 : iext + 16 * igp ) = -1

         ILOC = 8
         ! Face 11
         NEILOC( ILOC, 41 ) = 7
         NEILOC( ILOC, 42 ) = 7
         NEILOC( ILOC, 43 ) = 7
         NEILOC( ILOC, 44 ) = 7
         ! Face 9
         NEILOC( ILOC, 33 ) = 5
         NEILOC( ILOC, 34 ) = 5
         NEILOC( ILOC, 35 ) = 5
         NEILOC( ILOC, 36 ) = 5
         ! Face 12
         NEILOC( ILOC, 45 ) = 9
         NEILOC( ILOC, 46 ) = 9
         NEILOC( ILOC, 47 ) = 9
         NEILOC( ILOC, 48 ) = 9
         ! Face 20
         NEILOC( ILOC, 77 ) = 17
         NEILOC( ILOC, 78 ) = 17
         NEILOC( ILOC, 79 ) = 17
         NEILOC( ILOC, 80 ) = 17
         ! External
         neiloc( iloc, iext + 16 * igp + 1 : iext + 18 * igp ) = -1

         ILOC = 9
         ! Face 10
         NEILOC( ILOC, 37 ) = 6
         NEILOC( ILOC, 38 ) = 6
         NEILOC( ILOC, 39 ) = 6
         NEILOC( ILOC, 40 ) = 6
         ! Face 12
         NEILOC( ILOC, 45 ) = 8
         NEILOC( ILOC, 46 ) = 8
         NEILOC( ILOC, 47 ) = 8
         NEILOC( ILOC, 48 ) = 8
         ! Face 21
         NEILOC( ILOC, 81 ) = 18
         NEILOC( ILOC, 82 ) = 18
         NEILOC( ILOC, 83 ) = 18
         NEILOC( ILOC, 84 ) = 18
         ! External
         neiloc( iloc, iext + 18 * igp + 1 : iext + 21 * igp ) = -1

         ILOC = 10
         ! Face 13
         NEILOC( ILOC, 49 ) = 1
         NEILOC( ILOC, 50 ) = 1
         NEILOC( ILOC, 51 ) = 1
         NEILOC( ILOC, 52 ) = 1
         ! Face 22
         NEILOC( ILOC, 85 ) = 11
         NEILOC( ILOC, 86 ) = 11
         NEILOC( ILOC, 87 ) = 11
         NEILOC( ILOC, 88 ) = 11
         ! Face 24
         NEILOC( ILOC, 93 ) = 13
         NEILOC( ILOC, 94 ) = 13
         NEILOC( ILOC, 95 ) = 13
         NEILOC( ILOC, 96 ) = 13
         ! Face 34
         NEILOC( ILOC, 133 ) = 19
         NEILOC( ILOC, 134 ) = 19
         NEILOC( ILOC, 135 ) = 19
         NEILOC( ILOC, 136 ) = 19
         ! External
         neiloc( iloc, iext + 21 * igp + 1 : iext + 23 * igp ) = -1

         ILOC = 11
         ! Face 14
         NEILOC( ILOC, 53 ) = 2
         NEILOC( ILOC, 54 ) = 2
         NEILOC( ILOC, 55 ) = 2
         NEILOC( ILOC, 56 ) = 2
         ! Face 22
         NEILOC( ILOC, 85 ) = 10
         NEILOC( ILOC, 86 ) = 10
         NEILOC( ILOC, 87 ) = 10
         NEILOC( ILOC, 88 ) = 10
         ! Face 25
         NEILOC( ILOC, 97 ) = 14
         NEILOC( ILOC, 98 ) = 14
         NEILOC( ILOC, 99 ) = 14
         NEILOC( ILOC, 100 ) = 14
         ! Face 23
         NEILOC( ILOC, 89 ) = 12
         NEILOC( ILOC, 90 ) = 12
         NEILOC( ILOC, 91 ) = 12
         NEILOC( ILOC, 92 ) = 12
         ! Face 35
         NEILOC( ILOC, 137 ) = 20
         NEILOC( ILOC, 138 ) = 20
         NEILOC( ILOC, 139 ) = 20
         NEILOC( ILOC, 140 ) = 20
         ! External
         neiloc( iloc, iext + 23 * igp + 1 : iext + 24 * igp ) = -1

         ILOC = 12
         ! Face 15
         NEILOC( ILOC, 57 ) = 3
         NEILOC( ILOC, 58 ) = 3  
         NEILOC( ILOC, 59 ) = 3 
         NEILOC( ILOC, 60 ) = 3
         ! Face 23
         NEILOC( ILOC, 89 ) = 11
         NEILOC( ILOC, 90 ) = 11
         NEILOC( ILOC, 91 ) = 11
         NEILOC( ILOC, 92 ) = 11
         ! Face 26
         NEILOC( ILOC, 101 ) = 15 
         NEILOC( ILOC, 102 ) = 15
         NEILOC( ILOC, 103 ) = 15
         NEILOC( ILOC, 104 ) = 15
         ! Face 36
         NEILOC( ILOC, 141 ) = 21 
         NEILOC( ILOC, 142 ) = 21
         NEILOC( ILOC, 143 ) = 21
         NEILOC( ILOC, 144 ) = 21
         ! External
         neiloc( iloc, iext + 24 * igp + 1 : iext + 26 * igp ) = -1

         ILOC = 13
         ! Face 16
         NEILOC( ILOC, 61 ) = 4
         NEILOC( ILOC, 62 ) = 4
         NEILOC( ILOC, 63 ) = 4
         NEILOC( ILOC, 64 ) = 4
         ! Face 24
         NEILOC( ILOC, 93 ) = 10
         NEILOC( ILOC, 94 ) = 10
         NEILOC( ILOC, 95 ) = 10
         NEILOC( ILOC, 96 ) = 10
         ! Face 27
         NEILOC( ILOC, 105 ) = 14
         NEILOC( ILOC, 106 ) = 14
         NEILOC( ILOC, 107 ) = 14
         NEILOC( ILOC, 108 ) = 14
         ! Face 29
         NEILOC( ILOC, 113 ) = 16
         NEILOC( ILOC, 114 ) = 16
         NEILOC( ILOC, 115 ) = 16
         NEILOC( ILOC, 116 ) = 16
         ! Face 37
         NEILOC( ILOC, 145 ) = 22
         NEILOC( ILOC, 146 ) = 22
         NEILOC( ILOC, 147 ) = 22
         NEILOC( ILOC, 148 ) = 22
         ! External
         neiloc( iloc, iext + 26 * igp + 1 : iext + 27 * igp ) = -1

         ILOC = 14
         ! Face 17
         NEILOC( ILOC, 65 ) = 5
         NEILOC( ILOC, 66 ) = 5
         NEILOC( ILOC, 67 ) = 5
         NEILOC( ILOC, 68 ) = 5
         ! Face 25
         NEILOC( ILOC, 97 ) = 11
         NEILOC( ILOC, 98 ) = 11
         NEILOC( ILOC, 99 ) = 11
         NEILOC( ILOC, 100 ) = 11
         ! Face 28
         NEILOC( ILOC, 109 ) = 15
         NEILOC( ILOC, 110 ) = 15
         NEILOC( ILOC, 111 ) = 15
         NEILOC( ILOC, 112 ) = 15
         ! Face 30
         NEILOC( ILOC, 117 ) = 17
         NEILOC( ILOC, 118 ) = 17
         NEILOC( ILOC, 119 ) = 17
         NEILOC( ILOC, 120 ) = 17
         ! Face 27
         NEILOC( ILOC, 105 ) = 13
         NEILOC( ILOC, 106 ) = 13
         NEILOC( ILOC, 107 ) = 13
         NEILOC( ILOC, 108 ) = 13
         ! Face 38
         NEILOC( ILOC, 149 ) = 23
         NEILOC( ILOC, 150 ) = 23
         NEILOC( ILOC, 151 ) = 23
         NEILOC( ILOC, 152 ) = 23

         ILOC = 15
         ! Face 18
         NEILOC( ILOC, 69 ) = 6
         NEILOC( ILOC, 70 ) = 6
         NEILOC( ILOC, 71 ) = 6
         NEILOC( ILOC, 72 ) = 6
         ! Face 26
         NEILOC( ILOC, 101 ) = 12
         NEILOC( ILOC, 102 ) = 12
         NEILOC( ILOC, 103 ) = 12
         NEILOC( ILOC, 104 ) = 12
         ! Face 28
         NEILOC( ILOC, 109 ) = 14
         NEILOC( ILOC, 110 ) = 14
         NEILOC( ILOC, 111 ) = 14
         NEILOC( ILOC, 112 ) = 14
         ! Face 31
         NEILOC( ILOC, 121 ) = 18
         NEILOC( ILOC, 122 ) = 18
         NEILOC( ILOC, 123 ) = 18
         NEILOC( ILOC, 124 ) = 18
         ! Face 39
         NEILOC( ILOC, 153 ) = 24
         NEILOC( ILOC, 154 ) = 24
         NEILOC( ILOC, 155 ) = 24
         NEILOC( ILOC, 156 ) = 24
         ! External
         neiloc( iloc, iext + 27 * igp + 1 : iext + 28 * igp ) = -1

         ILOC = 16
         ! Face 19
         NEILOC( ILOC, 73 ) = 7
         NEILOC( ILOC, 74 ) = 7
         NEILOC( ILOC, 75 ) = 7
         NEILOC( ILOC, 76 ) = 7
         ! Face 29
         NEILOC( ILOC, 113 ) = 13
         NEILOC( ILOC, 114 ) = 13
         NEILOC( ILOC, 115 ) = 13
         NEILOC( ILOC, 116 ) = 13
         ! Face 32
         NEILOC( ILOC, 125 ) = 17
         NEILOC( ILOC, 126 ) = 17
         NEILOC( ILOC, 127 ) = 17
         NEILOC( ILOC, 128 ) = 17
         ! Face 40
         NEILOC( ILOC, 157 ) = 25
         NEILOC( ILOC, 158 ) = 25
         NEILOC( ILOC, 159 ) = 25
         NEILOC( ILOC, 160 ) = 25
         ! External
         neiloc( iloc, iext + 28 * igp + 1 : iext + 30 * igp ) = -1

         ILOC = 17
         ! Face 20
         NEILOC( ILOC, 77 ) = 8
         NEILOC( ILOC, 78 ) = 8
         NEILOC( ILOC, 79 ) = 8
         NEILOC( ILOC, 80 ) = 8
         ! Face 32
         NEILOC( ILOC, 125 ) = 16
         NEILOC( ILOC, 126 ) = 16
         NEILOC( ILOC, 127 ) = 16
         NEILOC( ILOC, 128 ) = 16
         ! Face 30
         NEILOC( ILOC, 117 ) = 14
         NEILOC( ILOC, 118 ) = 14
         NEILOC( ILOC, 119 ) = 14
         NEILOC( ILOC, 120 ) = 14
         ! Face 33
         NEILOC( ILOC, 129 ) = 18
         NEILOC( ILOC, 130 ) = 18
         NEILOC( ILOC, 131 ) = 18
         NEILOC( ILOC, 132 ) = 18
         ! Face 41
         NEILOC( ILOC, 161 ) = 26
         NEILOC( ILOC, 162 ) = 26
         NEILOC( ILOC, 163 ) = 26
         NEILOC( ILOC, 164 ) = 26
         ! External
         neiloc( iloc, iext + 30 * igp + 1 : iext + 31 * igp ) = -1

         ILOC = 18
         ! Face 21
         NEILOC( ILOC, 81 ) = 9
         NEILOC( ILOC, 82 ) = 9
         NEILOC( ILOC, 83 ) = 9
         NEILOC( ILOC, 84 ) = 9
         ! Face 31
         NEILOC( ILOC, 121 ) = 15
         NEILOC( ILOC, 122 ) = 15
         NEILOC( ILOC, 123 ) = 15
         NEILOC( ILOC, 124 ) = 15
         ! Face 33
         NEILOC( ILOC, 129 ) = 17
         NEILOC( ILOC, 130 ) = 17
         NEILOC( ILOC, 131 ) = 17
         NEILOC( ILOC, 132 ) = 17
         ! Face 42
         NEILOC( ILOC, 165 ) = 27
         NEILOC( ILOC, 166 ) = 27
         NEILOC( ILOC, 167 ) = 27
         NEILOC( ILOC, 168 ) = 27
         ! External
         neiloc( iloc, iext + 31 * igp + 1 : iext + 33 * igp ) = -1

         ILOC = 19
         ! Face 34
         NEILOC( ILOC, 133 ) = 10
         NEILOC( ILOC, 134 ) = 10
         NEILOC( ILOC, 135 ) = 10
         NEILOC( ILOC, 136 ) = 10
         ! Face 43
         NEILOC( ILOC, 169 ) = 20
         NEILOC( ILOC, 170 ) = 20
         NEILOC( ILOC, 171 ) = 20
         NEILOC( ILOC, 172 ) = 20
         ! Face 45
         NEILOC( ILOC, 177 ) = 22
         NEILOC( ILOC, 178 ) = 22
         NEILOC( ILOC, 179 ) = 22
         NEILOC( ILOC, 180 ) = 22
         ! External
         neiloc( iloc, iext + 33 * igp + 1 : iext + 36 * igp ) = -1

         ILOC = 20
         ! Face 35
         NEILOC( ILOC, 137 ) = 11
         NEILOC( ILOC, 138 ) = 11
         NEILOC( ILOC, 139 ) = 11
         NEILOC( ILOC, 140 ) = 11
         ! Face 43
         NEILOC( ILOC, 169 ) = 19
         NEILOC( ILOC, 170 ) = 19
         NEILOC( ILOC, 171 ) = 19
         NEILOC( ILOC, 172 ) = 19
         ! Face 46
         NEILOC( ILOC, 181 ) = 23
         NEILOC( ILOC, 182 ) = 23
         NEILOC( ILOC, 183 ) = 23
         NEILOC( ILOC, 184 ) = 23
         ! Face 44
         NEILOC( ILOC, 173 ) = 21
         NEILOC( ILOC, 174 ) = 21
         NEILOC( ILOC, 175 ) = 21
         NEILOC( ILOC, 176 ) = 21
         ! External
         neiloc( iloc, iext + 36 * igp + 1 : iext + 38 * igp ) = -1

         ILOC = 21
         ! Face 36
         NEILOC( ILOC, 141 ) = 12
         NEILOC( ILOC, 142 ) = 12
         NEILOC( ILOC, 143 ) = 12
         NEILOC( ILOC, 144 ) = 12
         ! Face 44
         NEILOC( ILOC, 173 ) = 20
         NEILOC( ILOC, 174 ) = 20
         NEILOC( ILOC, 175 ) = 20
         NEILOC( ILOC, 176 ) = 20
         ! Face 47
         NEILOC( ILOC, 185 ) = 24
         NEILOC( ILOC, 186 ) = 24
         NEILOC( ILOC, 187 ) = 24
         NEILOC( ILOC, 188 ) = 24
         ! External
         neiloc( iloc, iext + 38 * igp + 1 : iext + 41 * igp ) = -1

         ILOC = 22
         ! Face 37
         NEILOC( ILOC, 145 ) = 13
         NEILOC( ILOC, 146 ) = 13
         NEILOC( ILOC, 147 ) = 13
         NEILOC( ILOC, 148 ) = 13
         ! Face 45
         NEILOC( ILOC, 177 ) = 19
         NEILOC( ILOC, 178 ) = 19
         NEILOC( ILOC, 179 ) = 19
         NEILOC( ILOC, 180 ) = 19
         ! Face 48
         NEILOC( ILOC, 189 ) = 23
         NEILOC( ILOC, 190 ) = 23
         NEILOC( ILOC, 191 ) = 23
         NEILOC( ILOC, 192 ) = 23
         ! Face 50
         NEILOC( ILOC, 197 ) = 25
         NEILOC( ILOC, 198 ) = 25
         NEILOC( ILOC, 199 ) = 25
         NEILOC( ILOC, 200 ) = 25
         ! External
         neiloc( iloc, iext + 41 * igp + 1 : iext + 43 * igp ) = -1

         ILOC = 23
         ! Face 38
         NEILOC( ILOC, 149 ) = 14
         NEILOC( ILOC, 150 ) = 14
         NEILOC( ILOC, 151 ) = 14
         NEILOC( ILOC, 152 ) = 14
         ! Face 46
         NEILOC( ILOC, 181 ) = 20
         NEILOC( ILOC, 182 ) = 20
         NEILOC( ILOC, 183 ) = 20
         NEILOC( ILOC, 184 ) = 20
         ! Face 49
         NEILOC( ILOC, 193 ) = 24
         NEILOC( ILOC, 194 ) = 24
         NEILOC( ILOC, 195 ) = 24
         NEILOC( ILOC, 196 ) = 24
         ! Face 51
         NEILOC( ILOC, 201 ) = 26
         NEILOC( ILOC, 202 ) = 26
         NEILOC( ILOC, 203 ) = 26
         NEILOC( ILOC, 204 ) = 26
         ! Face 48
         NEILOC( ILOC, 189 ) = 22
         NEILOC( ILOC, 190 ) = 22
         NEILOC( ILOC, 191 ) = 22
         NEILOC( ILOC, 192 ) = 22
         ! External
         neiloc( iloc, iext + 43 * igp + 1 : iext + 44 * igp ) = -1

         ILOC = 24
         ! Face 39
         NEILOC( ILOC, 153 ) = 15
         NEILOC( ILOC, 154 ) = 15
         NEILOC( ILOC, 155 ) = 15
         NEILOC( ILOC, 156 ) = 15
         ! Face 47
         NEILOC( ILOC, 185 ) = 21
         NEILOC( ILOC, 186 ) = 21
         NEILOC( ILOC, 187 ) = 21
         NEILOC( ILOC, 188 ) = 21
         ! Face 49
         NEILOC( ILOC, 193 ) = 23
         NEILOC( ILOC, 194 ) = 23
         NEILOC( ILOC, 195 ) = 23
         NEILOC( ILOC, 196 ) = 23
         ! Face 52
         NEILOC( ILOC, 205 ) = 27
         NEILOC( ILOC, 206 ) = 27
         NEILOC( ILOC, 207 ) = 27
         NEILOC( ILOC, 208 ) = 27
         ! External
         neiloc( iloc, iext + 44 * igp + 1 : iext + 46 * igp ) = -1

         ILOC = 25
         ! Face 40
         NEILOC( ILOC, 157 ) = 16
         NEILOC( ILOC, 158 ) = 16
         NEILOC( ILOC, 159 ) = 16
         NEILOC( ILOC, 160 ) = 16
         ! Face 50
         NEILOC( ILOC, 197 ) = 22
         NEILOC( ILOC, 198 ) = 22
         NEILOC( ILOC, 199 ) = 22
         NEILOC( ILOC, 200 ) = 22
         ! Face 53
         NEILOC( ILOC, 209 ) = 26
         NEILOC( ILOC, 210 ) = 26
         NEILOC( ILOC, 211 ) = 26
         NEILOC( ILOC, 212 ) = 26
         ! External
         neiloc( iloc, iext + 46 * igp + 1 : iext + 49 * igp ) = -1

         ILOC = 26
         ! Face 41
         NEILOC( ILOC, 161 ) = 17
         NEILOC( ILOC, 162 ) = 17
         NEILOC( ILOC, 163 ) = 17
         NEILOC( ILOC, 164 ) = 17
         ! Face 53
         NEILOC( ILOC, 209 ) = 25
         NEILOC( ILOC, 210 ) = 25
         NEILOC( ILOC, 211 ) = 25
         NEILOC( ILOC, 212 ) = 25
         ! Face 51
         NEILOC( ILOC, 201 ) = 23
         NEILOC( ILOC, 202 ) = 23
         NEILOC( ILOC, 203 ) = 23
         NEILOC( ILOC, 204 ) = 23
         ! Face 54
         NEILOC( ILOC, 213 ) = 27
         NEILOC( ILOC, 214 ) = 27
         NEILOC( ILOC, 215 ) = 27
         NEILOC( ILOC, 216 ) = 27
         ! External
         neiloc( iloc, iext + 49 * igp + 1 : iext + 51 * igp ) = -1

         ILOC = 27
         ! Face 42
         NEILOC( ILOC, 165 ) = 18
         NEILOC( ILOC, 166 ) = 18
         NEILOC( ILOC, 167 ) = 18
         NEILOC( ILOC, 168 ) = 18
         ! Face 52
         NEILOC( ILOC, 205 ) = 24
         NEILOC( ILOC, 206 ) = 24
         NEILOC( ILOC, 207 ) = 24
         NEILOC( ILOC, 208 ) = 24
         ! Face 54
         NEILOC( ILOC, 213 ) = 26
         NEILOC( ILOC, 214 ) = 26
         NEILOC( ILOC, 215 ) = 26
         NEILOC( ILOC, 216 ) = 26
         ! External
         neiloc( iloc, iext + 51 * igp + 1 : iext + 54 * igp ) = -1

      CASE DEFAULT; FLExit( " Invalid integer for cv_ele_type " )

      END SELECT Conditional_Type

      RETURN     
    END SUBROUTINE VOLNEI

    subroutine U_Volnei( cv_ele_type, cv_nloc, u_nloc, scvngi, &
         cv_neiloc,   &
         u_on_face )
      !-----------------------------------------------------------------!
      !- This subroutine calculates U_ON_FACE, a logical that works    -!
      !- in a similar way of CV_ON_FACE.                               -!
      !-----------------------------------------------------------------!
      implicit none
      integer, intent( in ) :: cv_ele_type, cv_nloc, u_nloc, scvngi
      integer, dimension( cv_nloc, scvngi ), intent( in ) :: cv_neiloc
      logical, dimension( u_nloc, scvngi ), intent( inout ) :: u_on_face
      ! Local variables
      integer :: cv_iloc, u_iloc, gi, u_jloc, u_nloc2

      u_on_face = .false.
      u_nloc2 = u_nloc / cv_nloc
      ewrite(3,*)u_nloc, cv_nloc, u_nloc2

      do cv_iloc = 1, cv_nloc
         do gi = 1, scvngi
            ! if ( cv_neiloc( cv_iloc, gi ) == -1 ) then
            ewrite(3,*)'cv_neiloc:', u_nloc2, cv_iloc, gi, cv_neiloc( cv_iloc, gi )
            !  endif
         end do
      end do

      Conditional_ElementType: if( cv_ele_type <= 2 ) then
         u_jloc = 0
         Loop_CV: do cv_iloc = 1, cv_nloc
            Loop_U: do u_iloc = 1, u_nloc2
               Loop_GI: do gi = 1, scvngi
                  if (  cv_neiloc( cv_iloc, gi ) == -1 ) then
                     u_jloc = ( cv_iloc - 1 ) * u_nloc2 + u_iloc
                     u_on_face( u_jloc, gi ) = .true.
                     if( u_iloc == 1 ) cycle Loop_CV
                  end if
               end do Loop_GI
            end do Loop_U
         end do Loop_CV
      else
         u_jloc = 0
         do cv_iloc = 1, cv_nloc
            do u_iloc = 1, u_nloc2
               do gi = 1, scvngi 
                  if ( cv_neiloc( cv_iloc, gi ) == -1 ) then
                     u_jloc = ( cv_iloc - 1 ) * u_nloc2 + u_iloc
                     u_on_face( u_jloc, gi ) = .true.
                  end if
               end do
            end do
         end do
      end if Conditional_ElementType

      return
    end subroutine U_Volnei



    SUBROUTINE GAUSSILOC( FINDGPTS, COLGPTS, NCOLGPTS,&
         NEILOC, NLOC, SVNGI )
      !     ----------------------------------------------------    
      !     
      ! This subroutine calculates FINDGPTS,COLGPTS,NCOLGPTS 
      ! which contains given a local node ILOC the Gauss pts 
      ! that are used to integrate around this local node. 
      !     
      !     -------------------------------
      !     - date last modified : 12/02/2002
      !     -------------------------------

      implicit none

      INTEGER, intent( in ) :: NLOC, SVNGI    
      INTEGER, DIMENSION( NLOC + 1 ), intent( inout ) :: FINDGPTS     
      ! We have overestimated the size of COLGPTS.
      INTEGER, DIMENSION( NLOC * SVNGI ), intent( inout ) :: COLGPTS
      INTEGER, DIMENSION( NLOC,  SVNGI ), intent( in ) :: NEILOC
      INTEGER, intent( inout ) :: NCOLGPTS
      ! Local     
      INTEGER :: ILOC, COUNT, GI

      COUNT = 0

      DO ILOC = 1, NLOC

         FINDGPTS( ILOC ) = COUNT + 1

         DO GI = 1, SVNGI

            IF( NEILOC( ILOC, GI ) /= 0 ) THEN
               COUNT = COUNT + 1
               COLGPTS( COUNT ) = GI
            END IF
            !ewrite(3,*)'iloc,gi,NEILOC( ILOC, GI ):',iloc,gi,NEILOC( ILOC, GI )

         END DO

      END DO

      FINDGPTS( NLOC + 1 ) = COUNT + 1
      NCOLGPTS = COUNT
      !stop 2821

      RETURN 
    END SUBROUTINE GAUSSILOC



    SUBROUTINE DETNLXMAI( ELE, X,Y,Z, XONDGL, TOTELE, XNONODS, NLOC, MLOC, NGI, &
         NLX, NLY, NLZ, MLX, MLY, MLZ, WEIGHT, DETWEI, &
         NX, NY, NZ,  MX, MY, MZ,&
         A11, A12, A13, A21, A22, A23, A31, A32, A33 )
      implicit none
      ! Calculate DETWEI,NX,NY,NZ, MX,MY,MZ for element ELE. For coefficient in 
      ! the inverse mat of the Jacobian.
      INTEGER, intent( in ):: ELE, TOTELE, XNONODS, NLOC, MLOC, NGI
      REAL, DIMENSION( XNONODS ), intent( in ) :: X, Y, Z
      INTEGER, DIMENSION( TOTELE * NLOC), intent( in ) :: XONDGL
      ! Volume shape functions
      REAL, DIMENSION( NLOC, NGI ), intent( in ) :: NLX, NLY, NLZ
      REAL, DIMENSION( MLOC, NGI ), intent( in ) :: MLX, MLY, MLZ
      REAL, DIMENSION( NGI ), intent( in ) :: WEIGHT
      REAL, DIMENSION( NGI ), intent( inout ) :: DETWEI
      REAL, DIMENSION( NLOC, NGI ), intent( inout ) :: NX, NY, NZ
      REAL, DIMENSION( MLOC, NGI ), intent( inout ) :: MX, MY, MZ
      REAL, DIMENSION( NGI ), intent( inout ):: A11, A12, A13, A21, A22, A23, &
           A31, A32, A33
      ! Local variables
      REAL :: AGI, BGI, CGI, DGI, EGI, FGI, GGI, HGI, KGI, DETJ, VOLUME
      INTEGER :: GI, L, IGLX

      VOLUME = 0.
      Loop_NGI: DO GI = 1, NGI

         AGI = 0.
         BGI = 0.
         CGI = 0.
         DGI = 0.
         EGI = 0.
         FGI = 0.
         GGI = 0.
         HGI = 0.
         KGI = 0.

         Loop_NLOC: DO L = 1, NLOC

            IGLX = XONDGL(( ELE - 1 ) * NLOC + L ) 

            AGI = AGI + NLX( L, GI ) * X( IGLX ) 
            BGI = BGI + NLX( L, GI ) * Y( IGLX ) 
            CGI = CGI + NLX( L, GI ) * Z( IGLX ) 
            DGI = DGI + NLY( L, GI ) * X( IGLX ) 
            EGI = EGI + NLY( L, GI ) * Y( IGLX ) 
            FGI = FGI + NLY( L, GI ) * Z( IGLX ) 
            GGI = GGI + NLZ( L, GI ) * X( IGLX ) 
            HGI = HGI + NLZ( L, GI ) * Y( IGLX ) 
            KGI = KGI + NLZ( L, GI ) * Z( IGLX ) 

         END DO Loop_NLOC

         DETJ = AGI * ( EGI * KGI - FGI * HGI) &
              -BGI * ( DGI * KGI - FGI * GGI ) &
              +CGI * ( DGI * HGI - EGI *GGI)

         DETWEI( GI ) = ABS( DETJ ) * WEIGHT( GI )
         VOLUME = VOLUME + DETWEI( GI )

         A11( GI )= ( EGI * KGI-FGI * HGI ) / DETJ
         A21( GI )=-( DGI * KGI-FGI * GGI ) / DETJ
         A31( GI )= ( DGI * HGI-EGI * GGI ) / DETJ

         A12( GI )=-( BGI * KGI-CGI * HGI ) / DETJ
         A22( GI )= ( AGI * KGI-CGI * GGI ) / DETJ
         A32( GI )=-( AGI * HGI-BGI * GGI ) / DETJ

         A13( GI )= ( BGI * FGI-CGI * EGI ) / DETJ
         A23( GI )=-( AGI * FGI-CGI * DGI ) / DETJ
         A33( GI )= ( AGI * EGI-BGI * DGI ) / DETJ

         DO L=1,NLOC
            NX( L, GI )= A11( GI ) * NLX( L, GI ) + A12( GI ) * NLY( L, GI ) + A13( GI ) * NLZ( L, GI )
            NY( L, GI )= A21( GI ) * NLX( L, GI ) + A22( GI ) * NLY( L, GI ) + A23( GI ) * NLZ( L, GI )
            NZ( L, GI )= A31( GI ) * NLX( L, GI ) + A32( GI ) * NLY( L, GI ) + A33( GI ) * NLZ( L, GI )
         END DO

         DO L=1,MLOC
            MX( L, GI )= A11( GI ) * MLX( L, GI ) + A12( GI ) * MLY( L, GI ) + A13( GI ) * MLZ( L, GI )
            MY( L, GI )= A21( GI ) * MLX( L, GI ) + A22( GI ) * MLY( L, GI ) + A23( GI ) * MLZ( L, GI )
            MZ( L, GI )= A31( GI ) * MLX( L, GI ) + A32( GI ) * MLY( L, GI ) + A33( GI ) * MLZ( L, GI )
         END DO

      END DO Loop_NGI

      RETURN

    END SUBROUTINE DETNLXMAI



    SUBROUTINE DETNLXR_PLUS_U( ELE, X, Y, Z, XONDGL, TOTELE, NONODS, &
         X_NLOC, CV_NLOC, NGI, &
         N, NLX, NLY, NLZ, WEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
         NX, NY, NZ, &
         U_NLOC, UNLX, UNLY, UNLZ, UNX, UNY, UNZ ) 
      implicit none
      INTEGER, intent( in ) :: ELE, TOTELE, NONODS, X_NLOC, NGI, CV_NLOC, U_NLOC
      INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) :: XONDGL
      REAL, DIMENSION( NONODS ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( CV_NLOC, NGI ), intent( in ) :: N, NLX, NLY, NLZ 
      REAL, DIMENSION( NGI ), intent( in ) :: WEIGHT
      REAL, DIMENSION( NGI ), intent( inout ) :: DETWEI, RA
      REAL, intent( inout ) :: VOLUME
      LOGICAL, intent( in ) :: D1, D3, DCYL
      REAL, DIMENSION( X_NLOC, NGI ), intent( inout ) :: NX, NY, NZ
      REAL, DIMENSION( U_NLOC, NGI ), intent( inout ) :: UNLX, UNLY, UNLZ
      REAL, DIMENSION( U_NLOC, NGI ), intent( inout ) :: UNX, UNY, UNZ

      ! Local variables
      REAL, PARAMETER :: PIE = 3.141592654
      REAL :: AGI, BGI, CGI, DGI, EGI, FGI, GGI, HGI, KGI, A11, A12, A13, A21, &
           A22, A23, A31, A32, A33, DETJ, TWOPIE, RGI, rsum
      INTEGER :: GI, L, IGLX

      ewrite(3,*)' In Detnlxr_Plus_U'

      VOLUME = 0.

      Conditional_D3: IF( D3 ) THEN
         Loop_GI1: DO GI = 1, NGI

            AGI = 0.
            BGI = 0.
            CGI = 0.
            DGI = 0.
            EGI = 0.
            FGI = 0.
            GGI = 0.
            HGI = 0.
            KGI = 0.

            Loop_L1: DO L = 1, X_NLOC ! NB R0 does not appear here although the z-coord might be Z+R0.
               IGLX = XONDGL(( ELE - 1 ) * X_NLOC + L)

               AGI = AGI + NLX( L, GI) * X( IGLX ) 
               BGI = BGI + NLX( L, GI) * Y( IGLX ) 
               CGI = CGI + NLX( L, GI) * Z( IGLX ) 
               DGI = DGI + NLY( L, GI) * X( IGLX ) 
               EGI = EGI + NLY( L, GI) * Y( IGLX ) 
               FGI = FGI + NLY( L, GI) * Z( IGLX ) 
               GGI = GGI + NLZ( L, GI) * X( IGLX ) 
               HGI = HGI + NLZ( L, GI) * Y( IGLX ) 
               KGI = KGI + NLZ( L, GI) * Z( IGLX )
            END DO Loop_L1

            DETJ = AGI * ( EGI * KGI - FGI * HGI ) &
                 -BGI * ( DGI * KGI - FGI * GGI ) &
                 +CGI * ( DGI * HGI - EGI * GGI )
            DETWEI( GI ) = ABS( DETJ ) * WEIGHT( GI )
            RA( GI ) = 1.0
            VOLUME = VOLUME + DETWEI( GI )

            ! For coefficient in the inverse mat of the jacobian. 
            A11=   ( EGI * KGI - FGI * HGI ) / DETJ
            A21= - ( DGI * KGI - FGI * GGI ) / DETJ
            A31=   ( DGI * HGI - EGI * GGI ) / DETJ
            A12= - ( BGI * KGI - CGI * HGI ) / DETJ
            A22=   ( AGI * KGI - CGI * GGI ) / DETJ
            A32= - ( AGI * HGI - BGI * GGI ) / DETJ
            A13=   ( BGI * FGI - CGI * EGI ) / DETJ
            A23= - ( AGI * FGI - CGI * DGI ) / DETJ
            A33=   ( AGI * EGI - BGI * DGI ) / DETJ

            Loop_L2: DO L = 1, X_NLOC
               NX( L, GI ) = A11 * NLX( L, GI) + A12 * NLY( L, GI ) + A13 * NLZ( L, GI )
               NY( L, GI ) = A21 * NLX( L, GI) + A22 * NLY( L, GI ) + A23 * NLZ( L, GI )
               NZ( L, GI ) = A31 * NLX( L, GI) + A32 * NLY( L, GI ) + A33 * NLZ( L, GI )
            END DO Loop_L2

            Loop_L3: DO L = 1, U_NLOC
               UNX( L, GI ) = A11 * UNLX( L, GI) + A12 * UNLY( L, GI ) + A13 * UNLZ( L, GI )
               UNY( L, GI ) = A21 * UNLX( L, GI) + A22 * UNLY( L, GI ) + A23 * UNLZ( L, GI )
               UNZ( L, GI ) = A31 * UNLX( L, GI) + A32 * UNLY( L, GI ) + A33 * UNLZ( L, GI )
            END DO Loop_L3

         END DO Loop_GI1

      ELSE IF(.NOT.D1) THEN

         TWOPIE = 1.0 
         IF( DCYL ) TWOPIE = 2. * PIE
         !         rsum=0.0

         Loop_GI2: DO GI = 1, NGI

            RGI = 0.
            AGI = 0.
            BGI = 0.
            CGI = 0.
            DGI = 0.

            Loop_L4: DO L = 1, X_NLOC
               IGLX = XONDGL(( ELE - 1 ) * X_NLOC + L )
               AGI = AGI + NLX( L, GI ) * X( IGLX ) 
               BGI = BGI + NLX( L, GI ) * Y( IGLX ) 
               CGI = CGI + NLY( L, GI ) * X( IGLX ) 
               DGI = DGI + NLY( L, GI ) * Y( IGLX ) 
               RGI = RGI + N( L, GI ) * Y( IGLX )
            END DO Loop_L4

            IF( .NOT. DCYL ) RGI = 1.0

            DETJ = AGI * DGI - BGI * CGI 
            RA( GI ) = RGI
            DETWEI( GI ) = TWOPIE * RGI * DETJ * WEIGHT( GI )
            VOLUME = VOLUME + DETWEI( GI )

            Loop_L5: DO L = 1, X_NLOC
               NX( L, GI ) = (  DGI * NLX( L, GI ) - BGI * NLY( L, GI )) / DETJ
               NY( L, GI ) = ( -CGI * NLX( L, GI ) + AGI * NLY( L, GI )) / DETJ
               NZ( L, GI ) = 0.0
            END DO Loop_L5


            Loop_L6: DO L = 1, U_NLOC
               UNX( L, GI ) = (  DGI * UNLX( L, GI ) - BGI * UNLY( L, GI )) / DETJ
               UNY( L, GI ) = ( -CGI * UNLX( L, GI ) + AGI * UNLY( L, GI )) / DETJ
               UNZ( L, GI ) = 0.0
            END DO Loop_L6

         END DO Loop_GI2

         !ewrite(3,*)'ngi,sum(weight),rsum:',ngi,sum(weight),rsum

      ELSE ! FOR 1D...

         Loop_GI14: DO GI = 1, NGI

            AGI = 0.

            Loop_L15: DO L = 1, X_NLOC
               IGLX = XONDGL(( ELE - 1 ) * X_NLOC + L )
               AGI = AGI + NLX( L, GI ) * X( IGLX ) 
            END DO Loop_L15

            DETJ = AGI 
            RA( GI ) = 1.0
            DETWEI( GI ) = DETJ * WEIGHT( GI )
            VOLUME = VOLUME + DETWEI( GI )

            Loop_L16: DO L = 1, X_NLOC
               NX( L, GI ) = NLX( L, GI ) / DETJ
               NY( L, GI ) = 0.0
               NZ( L, GI ) = 0.0
            END DO Loop_L16

            Loop_L17: DO L = 1, U_NLOC
               UNX( L, GI ) = UNLX( L, GI ) / DETJ
               UNY( L, GI ) = 0.0
               UNZ( L, GI ) = 0.0
            END DO Loop_L17

         END DO Loop_GI14
      ENDIF Conditional_D3

      RETURN
    END SUBROUTINE DETNLXR_PLUS_U



    SUBROUTINE DETNNN( ELE, X, Y, Z, XONDGL, TOTELE, XNONODS, NLOC, NGI, &
         N, NLX, NLY, NLZ, WEIGHT, DETWEI, VOLUME, D3, DCYL) 
      implicit none
      INTEGER, intent( in ) :: ELE, TOTELE, XNONODS, NLOC, NGI
      INTEGER, DIMENSION( TOTELE * NLOC), intent( in ) :: XONDGL
      REAL, DIMENSION( XNONODS ), intent( in ) :: X, Y, Z
      ! Volume shape functions
      REAL, DIMENSION( NLOC, NGI ), intent( in ) :: N, NLX, NLY, NLZ
      REAL, DIMENSION( NGI ), intent( in ) :: WEIGHT
      REAL, DIMENSION( NGI ), intent( inout ) :: DETWEI
      REAL, intent( inout ) :: VOLUME
      LOGICAL, intent( in ) :: D3, DCYL
      ! Local variables
      REAL, PARAMETER :: PIE = 3.141592654
      REAL :: AGI, BGI, CGI, DGI, EGI, FGI, GGI, HGI, KGI, DETJ, TWOPIE, RGI
      INTEGER :: GI, L, IGLX

      VOLUME=0. 

      Conditional_Dimension: IF( D3 ) THEN

         Loop_NGI: DO GI = 1, NGI

            AGI = 0.
            BGI = 0.
            CGI = 0.
            DGI = 0.
            EGI = 0.
            FGI = 0.
            GGI = 0.
            HGI = 0.
            KGI = 0.

            DO L = 1, NLOC
               IGLX = XONDGL(( ELE - 1) * NLOC + L )
               ! NB R0 does not appear here although the z-coord might be Z +R0. 
               AGI = AGI + NLX( L, GI ) * X( IGLX ) 
               BGI = BGI + NLX( L, GI ) * Y( IGLX ) 
               CGI = CGI + NLX( L, GI ) * Z( IGLX ) 

               DGI = DGI + NLY( L, GI ) * X( IGLX ) 
               EGI = EGI + NLY( L, GI ) * Y( IGLX ) 
               FGI = FGI + NLY( L, GI ) * Z( IGLX ) 

               GGI = GGI + NLZ( L, GI ) * X( IGLX ) 
               HGI = HGI + NLZ( L, GI ) * Y( IGLX ) 
               KGI = KGI + NLZ( L, GI ) * Z( IGLX )
            END DO

            DETJ = AGI * ( EGI * KGI - FGI * HGI )  &
                 - BGI * ( DGI * KGI - FGI * GGI ) &
                 + CGI * ( DGI * HGI - EGI * GGI )
            DETWEI( GI ) = ABS( DETJ ) * WEIGHT( GI )
            VOLUME = VOLUME + DETWEI( GI )

         END DO Loop_NGI

      ELSE

         TWOPIE = 1.0 
         IF( DCYL ) TWOPIE = 2. * PIE

         DO GI = 1, NGI

            RGI = 0.
            AGI = 0.
            BGI = 0.
            CGI = 0.
            DGI = 0.

            DO L = 1, NLOC
               IGLX = XONDGL(( ELE - 1) * NLOC + L )
               AGI = AGI  +  NLX( L, GI ) * X( IGLX ) 
               BGI = BGI  +  NLX( L, GI ) * Y( IGLX ) 
               CGI = CGI  +  NLY( L, GI ) * X( IGLX ) 
               DGI = DGI  +  NLY( L, GI ) * Y( IGLX ) 
               RGI = RGI + N( L, GI ) * Y( IGLX )
            END DO

            IF(.NOT.DCYL) RGI = 1.0

            DETJ =  AGI * DGI - BGI * CGI 
            DETWEI( GI ) = TWOPIE * RGI * DETJ * WEIGHT( GI )
            VOLUME = VOLUME + DETWEI( GI )
         END DO

      END IF Conditional_Dimension

      RETURN
    END subroutine DETNNN


    subroutine detnlxmsup( nloc, mloc, ngi, &
         nlx, nly, nlz, mlx, mly, mlz, weight, detwei, &
         nx, ny, nz, mx, my, mz, &
         ncloc, nclx, ncly, nclz, xl, yl, zl, &
         a11, a12, a13, a21, a22, a23, a31, a32, a33 )
      implicit none
      ! Calculate detwei,nx,ny,nz, mx,my,mz for element ele for coefficient in 
      ! the inverse mat of the Jacobian.  This works for a superparametric FEM 
      ! nc, nclx, ncly, nclz are basis functions associated with coordinates.
      ! (xl,yl,zl) are the local coordinates of the nodes.
      integer, intent( in ) :: nloc, mloc, ngi, ncloc
      real, dimension( nloc, ngi ), intent( in ) ::  nlx, nly, nlz, mlx, mly, mlz
      real, dimension( ngi ), intent( in ) :: weight
      real, dimension( ngi ), intent( inout ) :: detwei
      real, dimension( nloc, ngi ), intent( inout ) :: nx, ny, nz, mx, my, mz 
      real, dimension( nloc, ngi ), intent( in ) :: nclx, ncly, nclz
      real, dimension( ncloc ), intent( in ) :: xl, yl, zl
      real, dimension( ngi ), intent( inout ) :: a11, a12, a13, a21, a22, &
           a23, a31, a32, a33
      ! Local variables
      real :: agi, bgi, cgi, dgi, egi, fgi, ggi, hgi, kgi, detj, volume
      integer :: gi, l

      volume = 0.
      Loop_GI: do gi = 1,ngi

         agi = 0.
         bgi = 0.
         cgi = 0.
         dgi = 0.
         egi = 0.
         fgi = 0.
         ggi = 0.
         hgi = 0.
         kgi = 0.

         Loop_NCLOC: do l = 1,ncloc
            agi = agi + nclx( l, gi ) * xl( l ) 
            bgi = bgi + nclx( l, gi ) * yl( l ) 
            cgi = cgi + nclx( l, gi ) * zl( l ) 

            dgi = dgi + ncly( l, gi ) * xl( l ) 
            egi = egi + ncly( l, gi ) * yl( l ) 
            fgi = fgi + ncly( l, gi ) * zl( l ) 

            ggi = ggi + nclz( l, gi ) * xl( l ) 
            hgi = hgi + nclz( l, gi ) * yl( l ) 
            kgi = kgi + nclz( l, gi ) * zl( l ) 
         end do Loop_NCLOC

         detj = agi * (egi * kgi - fgi * hgi) &
              - bgi * (dgi * kgi - fgi * ggi) &
              + cgi * (dgi * hgi - egi * ggi)

         detwei( gi ) = abs( detj ) * weight( gi )
         volume = volume + detwei( gi )

         a11( gi ) =    (egi * kgi - fgi * hgi) / detj
         a21( gi ) =  - (dgi * kgi - fgi * ggi) / detj
         a31( gi ) =    (dgi * hgi - egi * ggi) / detj

         a12( gi ) =  - (bgi * kgi - cgi * hgi) / detj
         a22( gi ) =    (agi * kgi - cgi * ggi) / detj
         a32( gi ) =  - (agi * hgi - bgi * ggi) / detj

         a13( gi ) =    (bgi * fgi - cgi * egi) / detj
         a23( gi ) =  - (agi * fgi - cgi * dgi) / detj
         a33( gi ) =    (agi * egi - bgi * dgi) / detj

         do l = 1,nloc
            nx( l, gi ) =  a11( gi ) * nlx( l, gi ) + a12( gi ) * nly( l, gi ) + a13( gi ) * nlz( l, gi )
            ny( l, gi ) =  a21( gi ) * nlx( l, gi ) + a22( gi ) * nly( l, gi ) + a23( gi ) * nlz( l, gi )
            nz( l, gi ) =  a31( gi ) * nlx( l, gi ) + a32( gi ) * nly( l, gi ) + a33( gi ) * nlz( l, gi )
         end do

         do l = 1,mloc
            mx( l, gi ) =  a11( gi ) * mlx( l, gi ) + a12( gi ) * mly( l, gi ) + a13( gi ) * mlz( l, gi )
            my( l, gi ) =  a21( gi ) * mlx( l, gi ) + a22( gi ) * mly( l, gi ) + a23( gi ) * mlz( l, gi )
            mz( l, gi ) =  a31( gi ) * mlx( l, gi ) + a32( gi ) * mly( l, gi ) + a33( gi ) * mlz( l, gi )
         enddo

      end do Loop_GI

      return

    end subroutine detnlxmsup

  end module shape_functions

