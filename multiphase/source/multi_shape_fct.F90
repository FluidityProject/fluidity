   

!!!!========================================!!!!
!!!!     SHAPE FUNCTIONS SUBRTS             !!!!
!!!!========================================!!!!


module shape_functions

contains

!!!
!!!  N_SHAPE_FUN AND RELATED SUBRTS & FUNCTIONS
!!!

  SUBROUTINE N_SHAPE_FUN( N, NLX, WEIGHT, U_NLOC, NGI ) 

    IMPLICIT NONE
    INTEGER, intent( in ) :: U_NLOC, NGI
    REAL, DIMENSION( U_NLOC, NGI ), intent( inout ) :: N, NLX
    REAL, DIMENSION( NGI ), intent( inout ) :: WEIGHT
    ! Local variables
    REAL, DIMENSION( : ), ALLOCATABLE :: WEIT, LX, LXP
    LOGICAL :: GETNDP
    INTEGER :: GI, CORN, GPOI
    
    write(357,*) 'In N_SHAPE_FUN'

    ALLOCATE( WEIT( NGI ))
    ALLOCATE( LX( NGI )) 
    ALLOCATE( LX( U_NLOC )) 


    LXP( 1 ) = -1.
    LXP( 2 ) =  0.
    LXP( 3 ) =  1.
    GETNDP = .FALSE.

    CALL LAGROT( WEIT, LX, NGI, GETNDP )

    Loop_Gauss: DO GI = 1, NGI

       DO CORN = 1, U_NLOC
          GPOI = GI
          N( CORN, GPOI ) = 0.5 * ( 1. + LXP( CORN ) * LX( GI ))
          NLX( CORN, GPOI ) = 0.5 * LXP( CORN )
          WEIGHT( GPOI ) = WEIT( GI )
       END DO

    END DO Loop_Gauss

    DEALLOCATE( WEIT )
    DEALLOCATE( LX )
    DEALLOCATE( LXP )
    
    write(357,*) 'Leaving N_SHAPE_FUN'

    RETURN

  END SUBROUTINE N_SHAPE_FUN



  SUBROUTINE LAGROT( WEIT, QUAPOS, NDGI, GETNDP )

    IMPLICIT NONE
    ! This computes the weight and points for standard Gaussian quadrature.
    ! IF(GETNDP) then get the POSITION OF THE NODES 
    ! AND DONT BOTHER WITH THE WEITS.

    INTEGER, intent( in ) ::  NDGI
    REAL, DIMENSION( NDGI ), intent( inout ) ::  WEIT, QUAPOS
    LOGICAL, intent( in ) ::  GETNDP
    ! Local variables
    LOGICAL :: WEIGHT
    INTEGER :: IG
    
    write(357,*) 'In LAGROT'

    IF( .NOT. GETNDP ) THEN

       WEIGHT = .TRUE.
       DO IG = 1, NDGI
          WEIT( IG ) = RGPTWE( IG, NDGI, WEIGHT )
       END DO

       WEIGHT = .FALSE.
       DO IG = 1, NDGI
          QUAPOS( IG ) = RGPTWE( IG, NDGI, WEIGHT )
       END DO

    ELSE

       IF( NDGI == 1 ) THEN
          QUAPOS( 1 ) = 0.
       ELSE
          DO IG = 1, NDGI
             QUAPOS( IG ) = -1. + 2. * REAL( IG - 1 ) / REAL( NDGI - 1 )
          END DO
       ENDIF

    END IF
    
    write(357,*) 'Leaving LAGROT'

    RETURN

  END SUBROUTINE LAGROT



  REAL FUNCTION RGPTWE( IG, ND, WEIGHT )

    IMPLICIT NONE

    ! If WEIGHT is TRUE in function RGPTWE then return the Gauss-pt weight else return the Gauss-pt. 
    ! There are ND Gauss points we are looking for either the weight or the x-coord of the IG'th Gauss point. 

    INTEGER :: IG, ND
    LOGICAL :: WEIGHT

    Loop_Weight: IF( WEIGHT ) THEN ! Gauss points weights

       SELECT CASE( ND )

       CASE( 1 ) ; RGPTWE = 2.0

       CASE( 2 ) ; RGPTWE = 1.0

       CASE( 3 ) 
          SELECT CASE( IG )
          CASE( 1, 3 ) ; RGPTWE = 0.555555555555556
          CASE( 2 )    ; RGPTWE = 0.888888888888889
          END SELECT

       CASE( 4 )
          SELECT CASE( IG )
          CASE( 1, 4 ) ; RGPTWE = 0.347854845137454
          CASE( 2, 3 ) ; RGPTWE = 0.652145154862546
          END SELECT

       CASE( 5 )
          SELECT CASE( IG )
          CASE( 1, 5 ) ; RGPTWE = 0.236926885056189
          CASE( 2, 4 ) ; RGPTWE = 0.478628670499366
          CASE( 3 )    ; RGPTWE = 0.568888888888889
          END SELECT

       CASE( 6 )
          SELECT CASE( IG )
          CASE( 1, 6 ) ; RGPTWE = 0.171324492379170
          CASE( 2, 5 ) ; RGPTWE = 0.360761573048139
          CASE( 3, 4 ) ; RGPTWE = 0.467913934572691
          END SELECT

       CASE( 7 )
          SELECT CASE( IG )
          CASE( 1, 7 ) ; RGPTWE = 0.129484966168870
          CASE( 2, 6 ) ; RGPTWE = 0.279705391489277
          CASE( 3, 5 ) ; RGPTWE = 0.381830050505119
          CASE( 4 )    ; RGPTWE = 0.417959183673469
          END SELECT

       CASE( 8 )
          SELECT CASE( IG )
          CASE( 1, 8 ) ; RGPTWE = 0.101228536290376
          CASE( 2, 7 ) ; RGPTWE = 0.222381034453374
          CASE( 3, 6 ) ; RGPTWE = 0.313706645877877
          CASE( 4, 5 ) ; RGPTWE = 0.362683783378362
          END SELECT

       CASE( 9 )
          SELECT CASE( IG )
          CASE( 1, 9 ) ; RGPTWE = 0.081274388361574
          CASE( 2, 8 ) ; RGPTWE = 0.180648160694857
          CASE( 3, 7 ) ; RGPTWE = 0.260610696402935
          CASE( 4, 6 ) ; RGPTWE = 0.312347077040003
          CASE( 5 )    ; RGPTWE = 0.330239355001260
          END SELECT

       CASE( 10 )
          SELECT CASE( IG )
          CASE( 1, 10 ) ; RGPTWE = 0.066671344308688
          CASE( 2, 9 )  ; RGPTWE = 0.149451349150581
          CASE( 3, 8 )  ; RGPTWE = 0.219086362515982
          CASE( 4, 7 )  ; RGPTWE = 0.269266719309996
          CASE( 5, 6 )  ; RGPTWE = 0.295524224714753
          END SELECT

       END SELECT

    ELSE

       SELECT CASE( ND ) ! Gauss points

       CASE( 1 ) ; RGPTWE = 0.0

       CASE( 2 ) ; RGPTWE = 0.577350269189626

       CASE( 3 )
          SELECT CASE( IG )
          CASE( 1, 3 ) ; RGPTWE = 0.774596669241483
          CASE( 2 )    ; RGPTWE = 0.0
          END SELECT

       CASE( 4 )
          SELECT CASE( IG )
          CASE( 1, 4 ) ; RGPTWE = 0.861136311594953
          CASE( 2, 3 ) ; RGPTWE = 0.339981043584856
          END SELECT

       CASE( 5 )
          SELECT CASE( IG )
          CASE( 1, 5 ) ; RGPTWE = 0.906179845938664
          CASE( 2, 4 ) ; RGPTWE = 0.538469310105683
          CASE( 3 )    ; RGPTWE = 0.0
          END SELECT

       CASE( 6 )
          SELECT CASE( IG )
          CASE( 1, 6 ) ; RGPTWE = 0.932469514203152
          CASE( 2, 5 ) ; RGPTWE = 0.661209386466265
          CASE( 3, 4 ) ; RGPTWE = 0.238619186083197
          END SELECT

       CASE( 7 )
          SELECT CASE( IG )
          CASE( 1, 7 ) ; RGPTWE = 0.949107912342759
          CASE( 2, 6 ) ; RGPTWE = 0.741531185599394
          CASE( 3, 5 ) ; RGPTWE = 0.405845151377397
          CASE( 4 )    ; RGPTWE = 0.0
          END SELECT

       CASE( 8 )
          SELECT CASE( IG )
          CASE( 1, 8 ) ; RGPTWE = 0.960289856497536
          CASE( 2, 7 ) ; RGPTWE = 0.796666477413627
          CASE( 3, 6 ) ; RGPTWE = 0.525532409916329
          CASE( 4, 5 ) ; RGPTWE = 0.183434642495650
          END SELECT

       CASE( 9 )
          SELECT CASE( IG )
          CASE( 1, 9 ) ; RGPTWE = 0.968160239507626
          CASE( 2, 8 ) ; RGPTWE = 0.836031107326636
          CASE( 3, 7 ) ; RGPTWE = 0.613371432700590
          CASE( 4, 6 ) ; RGPTWE = 0.324253423403809
          CASE( 5)     ; RGPTWE = 0.0
          END SELECT

       CASE( 10 )
          SELECT CASE( IG )
          CASE( 1, 10 ) ; RGPTWE = 0.973906528517172
          CASE( 2, 9 )  ; RGPTWE = 0.865063366688985
          CASE( 3, 8 )  ; RGPTWE = 0.679409568299024
          CASE( 4, 7 )  ; RGPTWE = 0.433395394129247
          CASE( 5, 6 )  ; RGPTWE = 0.148874338981631
          END SELECT

       END SELECT

       IF( IG <= INT( ( ND / 2 ) + 0.1 )) RGPTWE = -RGPTWE

    END IF Loop_Weight

    RETURN

  END FUNCTION RGPTWE

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
    
    write(357,*) 'In SHAPESE'

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
    
    write(357,*) 'Leaving SHAPESE'


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
    
    write(357,*) 'In SURNEI'


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
       ! - This is a little note for the future if we are: using 16 point quadrature for the surface elements
       !   then change the loop DO SGI = COUNT, COUNT to DO SGI = COUNT, COUNT + 3. Also change the loop
       !   increment to COUNT = COUNT + 4. This is similar for triangles.

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
       WRITE(357,*) 'SELETYP .GE. 7' 
       STOP 1234

    END SELECT
    
    write(357,*) 'Leaving SURNEI'

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
    
    write(357,*) 'In SFELIN'

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
    
    write(357,*) 'Leaving SFELIN'

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
    
    write(357,*) 'In SFEQUAD'

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
    
    write(357,*) 'Leaving SFEQUAD'


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
    
    write(357,*) 'In SFETRI'

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
    
    write(357,*) 'Leaving SFETRI'

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
    
    write(357,*) 'In SFEQLIN'
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
    
    write(357,*) 'Leaving SFEQLIN'

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
    
    write(357,*) 'In SFEQQUAD'

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
    
    write(357,*) 'Leaving SFEQQUAD'

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
    
    write(357,*) 'In SFEQTRI'

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
    
    write(357,*) 'Leaving SFEQTRI'

    RETURN
    !
  END SUBROUTINE SFEQTRI
  !

!!!
!!!   SHAPESV AND RELATED SUBRTS & FUNCTIONS
!!!


  SUBROUTINE SHAPE_CV_N( NDIM, CV_ELE_TYPE,  & 
       CV_NGI, CV_NLOC, U_NLOC, CVN, CVWEIGH, N, NLX, NLY, NLZ, UN, UNLX, UNLY, UNLZ )
    ! Shape functions associated with volume integration using both CV basis 
    ! functions CVN as well as FEM basis functions N (and its derivatives NLX, NLY, NLZ)
    ! also for velocity basis functions UN, UNLX, UNLY, UNLZ
    IMPLICIT NONE
    INTEGER, intent( in ) :: NDIM, CV_ELE_TYPE, CV_NGI, CV_NLOC, U_NLOC
    REAL, DIMENSION( CV_NGI ), intent( inout ) :: CVWEIGH
    REAL, DIMENSION( CV_NLOC, CV_NGI ), intent( inout ) :: CVN, N, NLX, NLY, NLZ
    REAL, DIMENSION( U_NLOC, CV_NGI ), intent( inout ) :: UN, UNLX, UNLY, UNLZ

    ! A Spectal element using Legendra, Lagrange or Chebichef polynomials. 
    !        CALL SPECTR(NGI,NLOC,MLOC, &
    !             M,WEIGHT,N,NLX,NLY,NLZ,D3,.NOT.D3, IPOLY,IQADRA)
    
    write(357,*) 'In SHAPE_CV_N'

    IF( CV_ELE_TYPE == 1 ) THEN

       IF( NDIM == 1 ) THEN
          CALL QUAD_1D_SHAP( CV_NGI, CV_NLOC, U_NLOC, CVN, CVWEIGH, N, NLX, UN, UNLX )
          NLY = 0.0
          NLZ = 0.0
          UNLY = 0.0
          UNLZ = 0.0
       ENDIF

    ENDIF
    
    write(357,*) 'Leaving SHAPE_CV_N'

  END SUBROUTINE SHAPE_CV_N




  SUBROUTINE QUAD_1D_SHAP( CV_NGI, CV_NLOC, U_NLOC, CVN, CVWEIGH, N, NLX, UN, UNLX )
    ! for quadatic elements
    ! Shape functions associated with volume integration using both CV basis 
    ! functions CVN as well as FEM basis functions N (and its derivatives NLX, NLY, NLZ)
    IMPLICIT NONE
    INTEGER, intent( in ) :: CV_NGI, CV_NLOC, U_NLOC
    REAL, DIMENSION( CV_NGI ), intent( inout ) :: CVWEIGH
    REAL, DIMENSION( CV_NLOC, CV_NGI ), intent( inout ) :: CVN, N, NLX
    REAL, DIMENSION( U_NLOC, CV_NGI ), intent( inout ) :: UN, UNLX
    ! local variables
    INTEGER, PARAMETER :: THREE = 3
    INTEGER :: GPOI, ICV, P, NQUAD, ILOC
    REAL, DIMENSION( : ), allocatable :: LX, WEI, RDUMMY, XI_MIN, XI_MAX, &
         CV_NODPOS, U_NODPOS
    REAL :: VOL_CV, LX_TRAN, LXGP
    LOGICAL :: GETNDP,DIFF,NDIFF
    
    write(357,*) 'In QUAD_1D_SHAP'

    NQUAD = CV_NGI/CV_NLOC

    ! Allocating memory
    ALLOCATE( LX( NQUAD ))
    ALLOCATE( WEI( NQUAD ))
    ALLOCATE( RDUMMY( MAX(CV_NGI,U_NLOC,CV_NLOC) ))
    ALLOCATE( XI_MIN( NQUAD ))
    ALLOCATE( XI_MAX( NQUAD ))
    ALLOCATE( CV_NODPOS( CV_NLOC ))
    ALLOCATE( U_NODPOS( U_NLOC ))

    IF( (CV_NGI /= 18) .and. (CV_NGI /= 15) .and. (CV_NGI /= 12) &
         .and. (CV_NGI /= 9) .and. (CV_NGI /= 6) .and. (CV_NGI /= 1 )) THEN
       WRITE( 357, * )'NOT THE RIGHT NO OF QUADATRUE PTS CV_NGI=', CV_NGI
       STOP 292
    ENDIF

    Case_CV_NLOC: SELECT CASE( CV_NLOC )
    CASE ( 4 );
       XI_MIN(1) = -1.0
       XI_MAX(1) = -1.0 + 1./3.

       XI_MIN(2) = -1.0 + 1./3.
       XI_MAX(2) = -0.0

       XI_MIN(3) = -0.0
       XI_MAX(3) = +1.0 - 1./3.

       XI_MIN(4) = +1.0 - 1./3.
       XI_MAX(4) = +1.0

    CASE( 3 );
       IF( .TRUE. ) then
          XI_MIN(1) = -1.0
          XI_MAX(1) = -0.5

          XI_MIN(2) = -0.5
          XI_MAX(2) = +0.5

          XI_MIN(3) = +0.5
          XI_MAX(3) = +1.0
       ENDIF
       IF( .FALSE. ) then
          XI_MIN(1) = -1.0
          XI_MAX(1) = -1.0+2.0/3.0

          XI_MIN(2) = -1.0+2.0/3.0
          XI_MAX(2) = +1.0-2.0/3.0

          XI_MIN(3) = +1.0-2.0/3.0
          XI_MAX(3) = +1.0
       ENDIF

    CASE( 2 );
       XI_MIN(1) = -1.0
       XI_MAX(1) = -0.0

       XI_MIN(2) = -0.0
       XI_MAX(2) = +1.0

    CASE( 1 );
       XI_MIN(1) = -1.0
       XI_MAX(1) =  1.0

    END SELECT Case_CV_NLOC

    DIFF=.TRUE.
    NDIFF=.FALSE.
    ! Find the roots of the quadrature points and nodes
    ! also get the weights.
    GETNDP=.TRUE.
    !     Compute standard Gauss quadrature. weits and points
    CALL LAGROT(RDUMMY,CV_NODPOS,CV_NLOC,GETNDP) 
    CALL LAGROT(RDUMMY,U_NODPOS,U_NLOC,GETNDP) 
    GETNDP=.FALSE.
    !     Compute standard Gauss quadrature. weits and points
    CALL LAGROT(WEI,LX,NQUAD,GETNDP)
    write(357,*)'cv_nloc,u_nloc, NQUAD, CV_NGI:',cv_nloc,u_nloc, NQUAD, CV_NGI
    write(357,*)'wei=',wei
    write(357,*)'lx=',lx
    write(357,*)'CV_nodpos=',CV_NODPOS
    write(357,*)'U_nodpos=',U_NODPOS

    GPOI = 0

    CVN = 0.0 

    DO ICV = 1, CV_NLOC

       VOL_CV = XI_MAX(ICV) - XI_MIN(ICV)
       !
       DO P = 1, NQUAD! Was loop 23
          GPOI = GPOI + 1
          CVN(ICV, GPOI ) = 1.0 
          ! Map to a new local coord system...   
!!!!      lx(p)=   1.0
          LX_TRAN = 0.5*( XI_MAX( ICV ) + XI_MIN( ICV )) &
               + 0.5*( XI_MAX( ICV ) - XI_MIN( ICV ))* LX(P)
  LXGP=LX_TRAN
   !      write(357,*)'LX_TRAN =',LX_TRAN
   !      stop 67
   ! Map to a new local coord system...      
   !          LX_TRAN = ( 2.*LX(P) - ( XI_MAX( ICV ) + XI_MIN( ICV ))) / ( XI_MAX( ICV ) - XI_MIN( ICV ))
   !           WRITE(357,*)'ICV,P,GPOI,LX_TRAN:',ICV,P,LX_TRAN
   !
          CVWEIGH( GPOI ) = WEI( P ) * VOL_CV / 2.0

          DO ILOC=1,CV_NLOC
             N(ILOC,GPOI)  =LAGRAN(NDIFF,LXGP,ILOC,CV_NLOC,CV_NODPOS)
             NLX(ILOC,GPOI)=LAGRAN(DIFF,LXGP, ILOC,CV_NLOC,CV_NODPOS)
  END DO

          DO ILOC=1,U_NLOC
             UN(ILOC,GPOI)  =LAGRAN(NDIFF,LXGP,ILOC,U_NLOC,U_NODPOS)
             UNLX(ILOC,GPOI)=LAGRAN(DIFF,LXGP, ILOC,U_NLOC,U_NODPOS)
  END DO

   ! endof DO P=1,NQUAD...
       END DO
       ! endof DO ICV=1,3...
    END DO

    !    IF(U_NLOC == CV_NLOC) THEN
    !       UN  =N
    !       UNLX=NLX
    !    ENDIF

    !    STOP 8764


    DEALLOCATE( LX )
    DEALLOCATE( WEI )
    DEALLOCATE( RDUMMY )
    DEALLOCATE( XI_MIN )
    DEALLOCATE( XI_MAX )
    DEALLOCATE( CV_NODPOS )
    DEALLOCATE( U_NODPOS )
    
    write(357,*) 'Leaving QUAD_1D_SHAP'


    RETURN
    !
  END SUBROUTINE QUAD_1D_SHAP


  !     
  !     
  REAL FUNCTION LAGRAN(DIFF,LX,INOD,NDNOD,NODPOS)
    IMPLICIT NONE
    !     This return the Lagrange poly assocaited with node INOD at point LX
    !     If DIFF then send back the value of this poly differentiated. 
    LOGICAL DIFF
    INTEGER INOD,NDNOD
    REAL LX,NODPOS(0:NDNOD-1)
    REAL DENOMI,OVER,OVER1
    INTEGER N,K,I,JJ
    !     ewrite(3,*) 'inside lagran'
    !     ewrite(3,*) 'DIFF,LX,INOD,NDNOD,NODPOS:',
    !     &            DIFF,LX,INOD,NDNOD,NODPOS
    !
    
    N=NDNOD-1
    K=INOD-1
    !     
    DENOMI=1.
    do I=0,K-1
       DENOMI=DENOMI*(NODPOS(K)-NODPOS(I))
    END DO
    do I=K+1,N
       DENOMI=DENOMI*(NODPOS(K)-NODPOS(I))
    END DO
    !     
    IF(.NOT.DIFF) THEN
       OVER=1.
       do I=0,K-1
          OVER=OVER*(LX-NODPOS(I))
       END DO
       do I=K+1,N
          OVER=OVER*(LX-NODPOS(I))
       END DO
       LAGRAN=OVER/DENOMI
    ELSE
       OVER=0.
       do JJ=0,N
          IF(JJ.NE.K) THEN
             OVER1=1.
             do I=0,K-1
                IF(JJ.NE.I) OVER1=OVER1*(LX-NODPOS(I))
             END DO
             do I=K+1,N
                IF(JJ.NE.I) OVER1=OVER1*(LX-NODPOS(I))
             END DO
             OVER=OVER+OVER1
          ENDIF
       END DO
       LAGRAN=OVER/DENOMI
    ENDIF
    !     
    !     ewrite(3,*) 'FINISHED LAGRAN'
  END FUNCTION LAGRAN




  SUBROUTINE RETRIEVE_NGI( CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, NFACE, &
       NDIM, CV_ELE_TYPE, CV_NLOC, U_NLOC)
    IMPLICIT NONE 
    ! determine number of quadatutre points. 
    INTEGER, intent( inout ) :: CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, NFACE
    INTEGER, intent( in ) :: NDIM, CV_ELE_TYPE, CV_NLOC, U_NLOC
    
    write(357,*) 'In RETRIEVE_NGI'

    IF( NDIM == 1) THEN
       IF(CV_NLOC == 1) THEN
          !         CV_NGI = 6
          !         CV_NGI = 8
          CV_NGI = 1
          SCVNGI = 2
       ENDIF
       IF(CV_NLOC == 2) THEN
          !         CV_NGI = 6
          !         CV_NGI = 8
          CV_NGI = 12
          SCVNGI = 3
       ENDIF
       IF(CV_NLOC == 3) THEN
          !         CV_NGI = 9
          CV_NGI = 12
          SCVNGI = 4
          !     IF(U_NLOC==4) CV_NGI = 12
          IF(U_NLOC == 4) CV_NGI = 18
          !     IF(U_NLOC==5) CV_NGI = 15
          IF(U_NLOC == 5) CV_NGI = 18
       ENDIF
       SBCVNGI = 1
       NFACE = 2
       CV_NGI_SHORT = CV_NGI
    ENDIF

    IF(CV_ELE_TYPE==2) THEN
       CV_NGI = CV_NGI * CV_NLOC
    ENDIF
    
    write(357,*) 'Leaving RETRIEVE_NGI'

    RETURN
  END SUBROUTINE RETRIEVE_NGI




  SUBROUTINE CV_FEM_SHAPE_FUNS( &
                                ! Volume shape functions...
       NDIM, CV_ELE_TYPE, & 
       CV_NGI, CV_NGI_SHORT, CV_NLOC, U_NLOC, CVN, CVN_SHORT, &
       CVWEIGHT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
       CVWEIGHT_SHORT, CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
       UFEN, UFENLX, UFENLY, UFENLZ, &
                                ! Surface of each CV shape functions...
       SCVNGI, CV_NEILOC, CV_ON_FACE, &  
       SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
       SCVFENLX, SCVFENLY, SCVFENLZ,  &
       SUFEN, SUFENSLX, SUFENSLY,  &
       SUFENLX, SUFENLY, SUFENLZ,  &
                                ! Surface element shape funcs...
       U_ON_FACE, NFACE, & 
       SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
       SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, SBUFENLY, SBUFENLZ, &
       CV_SLOCLIST, U_SLOCLIST, CV_SNLOC, U_SNLOC, &
                                ! Define the gauss points that lie on the surface of the CV...
       FINDGPTS, COLGPTS, NCOLGPTS, &
                                ! for over lapping elements else=1.0
       SELE_OVERLAP_SCALE ) 
    !     ======= DEFINE THE SUB-CONTROL VOLUME & FEM SHAPE FUNCTIONS ========
    IMPLICIT NONE
    INTEGER, intent( in ) :: NDIM, CV_ELE_TYPE, CV_NGI, CV_NGI_SHORT, CV_NLOC, U_NLOC
    INTEGER, intent( in ) :: SCVNGI
    REAL, DIMENSION( CV_NGI ), intent( inout ) :: CVWEIGHT
    REAL, DIMENSION( CV_NGI_SHORT ), intent( inout ) :: CVWEIGHT_SHORT
    REAL, DIMENSION( CV_NLOC, CV_NGI ), intent( inout ) :: CVN,  &
         CVFEN, CVFENLX, CVFENLY, CVFENLZ
    REAL, DIMENSION( CV_NLOC, CV_NGI_SHORT ), intent( inout ) :: CVN_SHORT,  &
         CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT
    REAL, DIMENSION( U_NLOC, CV_NGI ), intent( inout ) :: UFEN, UFENLX, UFENLY, UFENLZ
    INTEGER, DIMENSION( CV_NLOC, SCVNGI ), intent( inout ) :: CV_NEILOC
    ! CV_ON_FACE(CV_KLOC,GI)=.TRUE. if CV_KLOC is on the face that GI is centred on.
    LOGICAL, DIMENSION( CV_NLOC, SCVNGI ), intent( inout ) :: CV_ON_FACE
    LOGICAL, DIMENSION( U_NLOC, SCVNGI ), intent( inout ) :: U_ON_FACE
    REAL, DIMENSION( CV_NLOC, SCVNGI ), intent( inout ) :: SCVFEN, SCVFENSLX, SCVFENSLY, &       
         SCVFENLX, SCVFENLY, SCVFENLZ
    REAL, DIMENSION( SCVNGI ), intent( inout ) :: SCVFEWEIGH
    REAL, DIMENSION( U_NLOC, SCVNGI ), intent( inout ) :: SUFEN, SUFENSLX, SUFENSLY, &       
         SUFENLX, SUFENLY, SUFENLZ
    INTEGER, intent( in ) :: NFACE, SBCVNGI, CV_SNLOC, U_SNLOC
    REAL, DIMENSION( CV_SNLOC, SBCVNGI ), intent( inout ) :: SBCVFEN, SBCVFENSLX, SBCVFENSLY, &       
         SBCVFENLX, SBCVFENLY, SBCVFENLZ
    REAL, DIMENSION( SBCVNGI ), intent( inout ) :: SBCVFEWEIGH
    REAL, DIMENSION( U_SNLOC, SBCVNGI ), intent( inout ) :: SBUFEN, SBUFENSLX, SBUFENSLY, &       
         SBUFENLX, SBUFENLY, SBUFENLZ
    INTEGER, DIMENSION( NFACE, CV_SNLOC ), intent( inout ) ::  CV_SLOCLIST
    INTEGER, DIMENSION( NFACE, U_SNLOC ), intent( inout ) ::  U_SLOCLIST
    INTEGER, DIMENSION( CV_NLOC + 1 ), intent( inout ) :: FINDGPTS     
    ! We have overestimated the size of COLGPTS.
    INTEGER, DIMENSION( CV_NLOC * SCVNGI ), intent( inout ) :: COLGPTS
    INTEGER, intent( inout ) :: NCOLGPTS
    REAL, DIMENSION( CV_NLOC ), intent( inout ) :: SELE_OVERLAP_SCALE
    ! local variables
    INTEGER :: CV_ELE_TYPE2,U_NLOC2,ILEV,U_SNLOC2,U_ELE_TYPE2
    REAL, DIMENSION( :, : ), ALLOCATABLE :: UFEN2, UFENLX2, UFENLY2, UFENLZ2, &
         SUFEN2, SUFENSLX2, SUFENSLY2, SUFENLX2, SUFENLY2, SUFENLZ2,  &
         SBUFEN2, SBUFENSLX2, SBUFENSLY2, SBUFENLX2, SBUFENLY2, SBUFENLZ2
    INTEGER, DIMENSION( :, : ), ALLOCATABLE :: U_SLOCLIST2
    LOGICAL, DIMENSION( :, : ), ALLOCATABLE :: U_ON_FACE2

    ! Shape functions associated with volume integration using both CV basis 
    ! functions CVN as well as FEM basis functions CVFEN (and its derivatives CVFENLX, CVFENLY, CVFENLZ)

    write(357,*) 'in  cv_fem_shape_funs subrt'
    !write(357,*)'-CV_ELE_TYPE=',CV_ELE_TYPE

    ! SELE_OVERLAP_SCALE(P_JNOD) is the scaling needed to convert to overlapping element surfaces.     
    SELE_OVERLAP_SCALE = 1.0

    IF( CV_ELE_TYPE == 2 ) THEN
       IF( NDIM == 1 ) THEN
          IF( CV_NLOC >= 3 ) THEN
             SELE_OVERLAP_SCALE = 4.0
          ELSEIF( CV_NLOC == 2 ) THEN
             SELE_OVERLAP_SCALE = 2.0
          ENDIF
       ENDIF
    ENDIF

    IF( CV_ELE_TYPE == 2 ) THEN
       ! Define basis functions (overlapping) for porious media flow
       CV_ELE_TYPE2 = 1
       U_NLOC2 = U_NLOC / CV_NLOC

       ALLOCATE( UFEN2( U_NLOC2, CV_NGI_SHORT ))
       ALLOCATE( UFENLX2( U_NLOC2, CV_NGI_SHORT ))
       ALLOCATE( UFENLY2( U_NLOC2, CV_NGI_SHORT ))
       ALLOCATE( UFENLZ2( U_NLOC2, CV_NGI_SHORT ))

       CALL SHAPE_CV_N( NDIM, CV_ELE_TYPE2,  & 
            CV_NGI_SHORT, CV_NLOC, U_NLOC2, CVN_SHORT, CVWEIGHT_SHORT, &
            CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
            UFEN2, UFENLX2, UFENLY2, UFENLZ2 )

       ! Define the basis functions for each level.
       UFEN   = 0.0
       UFENLX = 0.0
       UFENLY = 0.0
       UFENLZ = 0.0
       CVN=0.0
       DO ILEV = 1, CV_NLOC
          UFEN(1+(ILEV-1)*U_NLOC2 : ILEV*U_NLOC2 , 1+(ILEV-1)*CV_NGI_SHORT : ILEV*CV_NGI_SHORT) &
               = UFEN2(1:U_NLOC2,1:CV_NGI_SHORT)
          UFENLX(1+(ILEV-1)*U_NLOC2 : ILEV*U_NLOC2 , 1+(ILEV-1)*CV_NGI_SHORT : ILEV*CV_NGI_SHORT) &
               = UFENLX2(1:U_NLOC2,1:CV_NGI_SHORT)
          UFENLY(1+(ILEV-1)*U_NLOC2 : ILEV*U_NLOC2 , 1+(ILEV-1)*CV_NGI_SHORT : ILEV*CV_NGI_SHORT) &
               = UFENLY2(1:U_NLOC2,1:CV_NGI_SHORT)
          UFENLZ(1+(ILEV-1)*U_NLOC2 : ILEV*U_NLOC2 , 1+(ILEV-1)*CV_NGI_SHORT : ILEV*CV_NGI_SHORT) &
               = UFENLZ2(1:U_NLOC2,1:CV_NGI_SHORT)

          CVFEN(1:CV_NLOC, 1+(ILEV-1)*CV_NGI_SHORT : ILEV*CV_NGI_SHORT) &
               = CVFEN_SHORT(1:CV_NLOC, 1:CV_NGI_SHORT)
          CVFENLX(1:CV_NLOC, 1+(ILEV-1)*CV_NGI_SHORT : ILEV*CV_NGI_SHORT) &
               = CVFENLX_SHORT(1:CV_NLOC, 1:CV_NGI_SHORT)
          CVFENLY(1:CV_NLOC, 1+(ILEV-1)*CV_NGI_SHORT : ILEV*CV_NGI_SHORT) &
               = CVFENLY_SHORT(1:CV_NLOC, 1:CV_NGI_SHORT)
          CVFENLZ(1:CV_NLOC, 1+(ILEV-1)*CV_NGI_SHORT : ILEV*CV_NGI_SHORT) &
               = CVFENLZ_SHORT(1:CV_NLOC, 1:CV_NGI_SHORT)

          CVWEIGHT(1+(ILEV-1)*CV_NGI_SHORT : ILEV*CV_NGI_SHORT) = CVWEIGHT_SHORT(1:CV_NGI_SHORT)
          CVN(ILEV, 1+(ILEV-1)*CV_NGI_SHORT : ILEV*CV_NGI_SHORT) = 1.0 
       END DO
    ELSE
       CALL SHAPE_CV_N( NDIM, CV_ELE_TYPE,  & 
            CV_NGI, CV_NLOC, U_NLOC, CVN, CVWEIGHT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
            UFEN, UFENLX, UFENLY, UFENLZ )
       CVN_SHORT=CVN
       CVFEN_SHORT=CVFEN
       CVFENLX_SHORT=CVFENLX
       CVFENLY_SHORT=CVFENLY
       CVFENLZ_SHORT=CVFENLZ
       CVWEIGHT_SHORT=CVWEIGHT
    ENDIF
    write(357,*) 'out of shape_cv_n - CVWEIGHT', CVWEIGHT
    ! SCVFEN(CV_NLOC,SCVNGI)         - the shape function evaluated 
    !                                  for each node at each surface gauss point
    ! SCVFENSLX[X/Y](CV_NLOC,SCVNGI) - the surface derivatives of the shape 
    !                                  function for each node at those same points, and
    ! SCVFENLX[X/Y](CV_NLOC,SCVNGI)  - the derivatives of the shape 
    !                                  function for each node at those same points, and
    ! SCVFEWEIGH(SCVNGI)             - the Gauss weights to use when integrating around 
    !                                  the control volume surface
    ! CV_NEILOC(CV_NLOC,SCVNGI)      - neighbour node for a given node/gauss-point pair 
    !
    ! This also include quadature points around the element. 

    IF( CV_ELE_TYPE == 2 ) THEN
       ! Define basis functions (overlapping) for porious media flow
       CV_ELE_TYPE2 = 1
       U_ELE_TYPE2 = 1
       U_NLOC2 = U_NLOC / CV_NLOC

       ALLOCATE( U_ON_FACE2( U_NLOC2, SCVNGI ))

       ALLOCATE( SUFEN2( U_NLOC2, SCVNGI ))
       ALLOCATE( SUFENSLX2( U_NLOC2, SCVNGI ))
       ALLOCATE( SUFENSLY2( U_NLOC2, SCVNGI ))
       ALLOCATE( SUFENLX2( U_NLOC2, SCVNGI ))
       ALLOCATE( SUFENLY2( U_NLOC2, SCVNGI ))
       ALLOCATE( SUFENLZ2( U_NLOC2, SCVNGI ))

       CALL SHAPESV_FEM_PLUS( SCVNGI, CV_NEILOC, CV_ON_FACE, &  
            CV_ELE_TYPE2, CV_NLOC, SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
            SCVFENLX, SCVFENLY, SCVFENLZ,  &
            U_NLOC2, SUFEN2, SUFENSLX2, SUFENSLY2,  &
            SUFENLX2, SUFENLY2, SUFENLZ2,  &
            NDIM )
       SUFEN=0.0
       SUFENSLX=0.0
       SUFENSLY=0.0
       SUFENLX=0.0
       SUFENLY=0.0
       SUFENLZ=0.0
       DO ILEV=1,CV_NLOC
          U_ON_FACE( 1+(ILEV-1)*U_NLOC2 : ILEV*U_NLOC2, 1:SCVNGI )  &
               =U_ON_FACE2( 1:U_NLOC2, 1:SCVNGI )

          SUFEN(1+(ILEV-1)*U_NLOC2 : ILEV*U_NLOC2 , 1 : SCVNGI ) &
               = SUFEN2(1:U_NLOC2, 1 : SCVNGI)
          SUFENSLX(1+(ILEV-1)*U_NLOC2 : ILEV*U_NLOC2 , 1 : SCVNGI ) &
               = SUFENSLX2(1:U_NLOC2, 1 : SCVNGI)
          SUFENSLY(1+(ILEV-1)*U_NLOC2 : ILEV*U_NLOC2 , 1 : SCVNGI ) &
               = SUFENSLY2(1:U_NLOC2, 1 : SCVNGI)
          SUFENLX(1+(ILEV-1)*U_NLOC2 : ILEV*U_NLOC2 , 1 : SCVNGI ) &
               = SUFENLX2(1:U_NLOC2, 1 : SCVNGI)
          SUFENLY(1+(ILEV-1)*U_NLOC2 : ILEV*U_NLOC2 , 1 : SCVNGI ) &
               = SUFENLY2(1:U_NLOC2, 1 : SCVNGI)
          SUFENLZ(1+(ILEV-1)*U_NLOC2 : ILEV*U_NLOC2 , 1 : SCVNGI ) &
               = SUFENLZ2(1:U_NLOC2, 1 : SCVNGI)
       END DO
    ELSE   
       CALL SHAPESV_FEM_PLUS( SCVNGI, CV_NEILOC, CV_ON_FACE, &  
            CV_ELE_TYPE, CV_NLOC, SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
            SCVFENLX, SCVFENLY, SCVFENLZ,  &
            U_NLOC, SUFEN, SUFENSLX, SUFENSLY,  &
            SUFENLX, SUFENLY, SUFENLZ,  &
            NDIM )
    ENDIF
    write(357,*) 'out of shapesv_fem_+ - SCVFEWEIGH', SCVFEWEIGH
    write(357,*) 'CV_ELE_TYPE', CV_ELE_TYPE

    ! Determine the surface element shape functions from those 
    ! calculated in SHAPESV_FEM_PLUS and also CV_SLOCLIST( NFACE,CV_SNLOC )
    IF(CV_ELE_TYPE==2) THEN
       ! Define basis functions (overlapping) for porious media flow
       CV_ELE_TYPE2=1
       U_ELE_TYPE2=1
       U_NLOC2=U_NLOC/CV_NLOC
       U_SNLOC2=U_SNLOC/CV_NLOC
       write(357,*)'U_NLOC,CV_NLOC:',U_NLOC,CV_NLOC
       write(357,*)'U_NLOC2,U_SNLOC2:',U_NLOC2,U_SNLOC2
       !!       stop 3921

       ALLOCATE( SBUFEN2( U_SNLOC2, SBCVNGI ))
       ALLOCATE( SBUFENSLX2( U_SNLOC2, SBCVNGI ))
       ALLOCATE( SBUFENSLY2( U_SNLOC2, SBCVNGI ))
       ALLOCATE( SBUFENLX2( U_SNLOC2, SBCVNGI ))
       ALLOCATE( SBUFENLY2( U_SNLOC2, SBCVNGI ))
       ALLOCATE( SBUFENLZ2( U_SNLOC2, SBCVNGI ))

       ALLOCATE( U_SLOCLIST2( 1:NFACE, U_SNLOC2 ))

       CALL DET_SUF_ELE_SHAPE( SCVNGI, CV_ON_FACE, U_ON_FACE2, NFACE, & 
            CV_NLOC, SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
            SCVFENLX, SCVFENLY, SCVFENLZ,  &
            U_NLOC2,  SUFEN2, SUFENSLX2, SUFENSLY2,  &
            SUFENLX2, SUFENLY2, SUFENLZ2,  &
            SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
            SBUFEN2, SBUFENSLX2, SBUFENSLY2, SBUFENLX2, SBUFENLY2, SBUFENLZ2, &
            CV_SLOCLIST, U_SLOCLIST2, CV_SNLOC, U_SNLOC2 )

       DO ILEV=1,CV_NLOC
          SBUFEN(1+(ILEV-1)*U_SNLOC2 : ILEV*U_SNLOC2 , 1 : SBCVNGI ) &
               = SBUFEN2(1:U_SNLOC2, 1 : SBCVNGI)
          SBUFENSLX(1+(ILEV-1)*U_SNLOC2 : ILEV*U_SNLOC2 , 1 : SBCVNGI ) &
               = SBUFENSLX2(1:U_SNLOC2, 1 : SBCVNGI)
          SBUFENSLY(1+(ILEV-1)*U_SNLOC2 : ILEV*U_SNLOC2 , 1 : SBCVNGI ) &
               = SBUFENSLY2(1:U_SNLOC2, 1 : SBCVNGI)
          SBUFENLX(1+(ILEV-1)*U_SNLOC2 : ILEV*U_SNLOC2 , 1 : SBCVNGI ) &
               = SBUFENLX2(1:U_SNLOC2, 1 : SBCVNGI)
          SBUFENLY(1+(ILEV-1)*U_SNLOC2 : ILEV*U_SNLOC2 , 1 : SBCVNGI ) &
               = SBUFENLY2(1:U_SNLOC2, 1 : SBCVNGI)
          SBUFENLZ(1+(ILEV-1)*U_SNLOC2 : ILEV*U_SNLOC2 , 1 : SBCVNGI ) &
               = SBUFENLZ2(1:U_SNLOC2, 1 : SBCVNGI)

          U_SLOCLIST(1:NFACE, 1+(ILEV-1)*U_SNLOC2 : ILEV*U_SNLOC2 )  &
               = U_SLOCLIST2(1:NFACE, 1:U_SNLOC2) + (ILEV-1)*U_NLOC2
       END DO

    ELSE
       CALL DET_SUF_ELE_SHAPE( SCVNGI, CV_ON_FACE, U_ON_FACE, NFACE, & 
            CV_NLOC, SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
            SCVFENLX, SCVFENLY, SCVFENLZ,  &
            U_NLOC,  SUFEN, SUFENSLX, SUFENSLY,  &
            SUFENLX, SUFENLY, SUFENLZ,  &
            SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
            SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, SBUFENLY, SBUFENLZ, &
            CV_SLOCLIST, U_SLOCLIST, CV_SNLOC, U_SNLOC )
    ENDIF

    ! Define the gauss points that lie on the surface of the
    ! control volume surrounding a given local node (iloc)
    ! that is FINDGPTS, COLGPTS, NCOLGPTS
    CALL GAUSSILOC( FINDGPTS, COLGPTS, NCOLGPTS, &
         CV_NEILOC, CV_NLOC, SCVNGI )

    write(357,*) 'leaving cv_fem_shape_funs subrt', NCOLGPTS
    WRITE(357,*)'CV_ON_FACE:',CV_ON_FACE
    WRITE(357,*)'U_ON_FACE:',U_ON_FACE
    WRITE(357,*)'CV_SLOCLIST=',CV_SLOCLIST
    WRITE(357,*)'U_SLOCLIST=',U_SLOCLIST
    !    STOP 34

    RETURN
  END SUBROUTINE CV_FEM_SHAPE_FUNS




  SUBROUTINE DET_SUF_ELE_SHAPE( SCVNGI, CV_ON_FACE, U_ON_FACE, NFACE, & 
       CV_NLOC, SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
       SCVFENLX, SCVFENLY, SCVFENLZ,  &
       U_NLOC,  SUFEN, SUFENSLX, SUFENSLY,  &
       SUFENLX, SUFENLY, SUFENLZ,  &
       SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
       SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, SBUFENLY, SBUFENLZ, &
       CV_SLOCLIST, U_SLOCLIST, CV_SNLOC, U_SNLOC )
    !     
    !     - this subroutine generates the FE basis functions, weights and the
    !     - derivatives of the shape functions for a variety of elements on the 
    !     - control volume boundaries.  
    !     - The routine also generates the shape functions and derivatives 
    !     - associated with the CV surfaces and also the FV basis functions.
    !     -------------------------------
    !     - date last modified : 24/05/2003
    !     -------------------------------

    IMPLICIT NONE

    INTEGER, intent( in ) :: SCVNGI, CV_NLOC, U_NLOC, NFACE, &
         SBCVNGI, CV_SNLOC, U_SNLOC
    ! CV_ON_FACE(CV_KLOC,GI)=.TRUE. if CV_KLOC is on the face that GI is centred on.
    LOGICAL, DIMENSION( CV_NLOC, SCVNGI ), intent( in ) :: CV_ON_FACE
    LOGICAL, DIMENSION( U_NLOC, SCVNGI ), intent( in ) :: U_ON_FACE
    REAL, DIMENSION( CV_NLOC, SCVNGI ), intent( in ) :: SCVFEN, SCVFENSLX, SCVFENSLY, &       
         SCVFENLX, SCVFENLY, SCVFENLZ
    REAL, DIMENSION( SCVNGI ), intent( inout ) :: SCVFEWEIGH
    REAL, DIMENSION( U_NLOC, SCVNGI ), intent( in ) :: SUFEN, SUFENSLX, SUFENSLY, &       
         SUFENLX, SUFENLY, SUFENLZ

    REAL, DIMENSION( CV_SNLOC, SBCVNGI ), intent( inout ) :: SBCVFEN, SBCVFENSLX, SBCVFENSLY, &       
         SBCVFENLX, SBCVFENLY, SBCVFENLZ
    REAL, DIMENSION( SBCVNGI ), intent( inout ) :: SBCVFEWEIGH
    REAL, DIMENSION( U_SNLOC, SBCVNGI ), intent( inout ) :: SBUFEN, SBUFENSLX, SBUFENSLY, &       
         SBUFENLX, SBUFENLY, SBUFENLZ
    INTEGER, DIMENSION( NFACE, CV_SNLOC ), intent( inout ) ::  CV_SLOCLIST
    INTEGER, DIMENSION( NFACE, U_SNLOC ), intent( inout ) ::  U_SLOCLIST
    ! Local variables
    INTEGER :: GI, GI2, GIS, CV_KLOC, CV_SKLOC, U_KLOC, U_SKLOC, ii
    LOGICAL :: FOUND

    FOUND = .FALSE. 
    DO CV_KLOC = 1, CV_NLOC
       DO GI2 = 1, SCVNGI
          IF( .NOT. FOUND ) THEN
             IF( CV_ON_FACE( CV_KLOC, GI2 )) THEN
                GI = GI2
                FOUND= .TRUE.
             END IF
          END IF
       END DO
    END DO


    ! Use all quadrature pts on face with GI on it. Determine SBCVFEN: 
    CV_SKLOC = 0
    Loop_CVKLOC: DO CV_KLOC = 1, CV_NLOC
       IF( CV_ON_FACE( CV_KLOC, GI )) THEN ! We are on the correct face
          CV_SKLOC = CV_SKLOC + 1
          GIS = 0

          DO GI2 = 1, SCVNGI
             write(357,*) 'SCVFEN( CV_KLOC, GI2 ):', ( SCVFEN( CV_KLOC, ii ), ii=1,SCVNGI)
             IF( ABS( SCVFEN( CV_KLOC, GI2 )) > 1.E-6 ) THEN
                IF( GIS < SBCVNGI ) THEN
                   GIS = GIS + 1
                   !WRITE(357,*) 'CV_SNLOC, SBCVNGI, SCVNGI:',CV_SNLOC, SBCVNGI, SCVNGI
                   SBCVFEN( CV_SKLOC, GIS ) = SCVFEN( CV_KLOC, GI2 )
                   SBCVFENSLX( CV_SKLOC, GIS ) = SCVFENSLX( CV_KLOC, GI2 )
                   SBCVFENSLY( CV_SKLOC, GIS ) = SCVFENSLY( CV_KLOC, GI2 )     
                   SBCVFENLX( CV_SKLOC, GIS ) = SCVFENLX( CV_KLOC, GI2 )
                   SBCVFENLY( CV_SKLOC, GIS ) = SCVFENLY( CV_KLOC, GI2 )
                   SBCVFENLZ( CV_SKLOC, GIS ) = SCVFENLZ( CV_KLOC, GI2 )
                   SBCVFEWEIGH( GIS ) = SCVFEWEIGH( GI2 )
                END IF
             END IF
          END DO

          IF(GIS /= SBCVNGI ) THEN
             WRITE(357,*)'SOMETHING GONE WRONG WITH NO OF QUADRATURE PTS'
             WRITE(357,*)'GIS, SBCVNGI:',GIS, SBCVNGI
             STOP 2382
          ENDIF

       ENDIF
    END DO Loop_CVKLOC

    ! Use all quadrature pts on face with GI on it. Determine SBUFEN: 
    write(357,*)'SBCVNGI=',SBCVNGI
    IF(CV_SNLOC==1) THEN
       SBUFEN= SBCVFEN
       SBUFENSLX=SBCVFENSLX
       SBUFENSLY=SBCVFENSLY 
       SBUFENLX=SBCVFENLX
       SBUFENLY=SBCVFENLY
       SBUFENLZ=SBCVFENLZ
    ELSE
       WRITE(357,*)'NOT READY' 
       STOP 289
       U_SKLOC = 0
       Loop_UKLOC: DO U_KLOC = 1, U_NLOC
          !       IF( U_ON_FACE( U_KLOC, GI )) THEN ! We are on the correct face
          U_SKLOC = U_SKLOC + 1
          GIS = 0
          DO GI2 = 1, SCVNGI

             IF(ABS( SUFEN( U_KLOC, GI2 )) > 1.E-6 ) THEN
                IF( GIS < SBCVNGI ) THEN
                   GIS = GIS + 1
          ! WRITE(357,*)'U_SKLOC, GIS ,U_KLOC, GI2, SUFEN( U_KLOC, GI2 ):', &
                  !      U_SKLOC, GIS ,U_KLOC, GI2, SUFEN( U_KLOC, GI2 )
                   SBUFEN( U_SKLOC, GIS ) = SUFEN( U_KLOC, GI2 )
                   SBUFENSLX( U_SKLOC, GIS ) = SUFENSLX( U_KLOC, GI2 )
                   SBUFENSLY( U_SKLOC, GIS ) = SUFENSLY( U_KLOC,GI2 )      
                   SBUFENLX( U_SKLOC, GIS ) = SUFENLX( U_KLOC, GI2 )
                   SBUFENLY( U_SKLOC, GIS ) = SUFENLY( U_KLOC, GI2 )
                   SBUFENLZ( U_SKLOC, GIS ) = SUFENLZ( U_KLOC, GI2 )
                END IF
             END IF
          END DO

          !       ENDIF
       END DO Loop_UKLOC
    ENDIF

    ! Determine U_ON_FACE
  !  DO U_KLOC = 1, U_NLOC
  !     write(357,*)'SUFEN( U_KLOC, GI):',(SUFEN( U_KLOC, GI),gi=1,SCVNGI)
  !  end do

    ! determine CV_SLOCLIST & U_SLOCLIST
    CALL DETERMIN_SLOCLIST( CV_SLOCLIST, CV_NLOC, CV_SNLOC, SCVNGI, CV_ON_FACE, NFACE )
    IF(U_SNLOC==1) THEN
       U_SLOCLIST(1,1)=1
       U_SLOCLIST(2,1)=U_NLOC
    ELSE
       WRITE(357,*)'NOT YET READY'
       STOP 3831
       CALL DETERMIN_SLOCLIST( U_SLOCLIST, U_NLOC, U_SNLOC, SCVNGI, U_ON_FACE, NFACE )
    ENDIF
    write(357,*)'CV_SNLOC,U_SNLOC,SCVNGI:',CV_SNLOC,U_SNLOC,SCVNGI
    write(357,*)'CV_SLOCLIST:',CV_SLOCLIST
    write(357,*)'U_SLOCLIST:',U_SLOCLIST
    !    STOP 3893

    RETURN
  END SUBROUTINE DET_SUF_ELE_SHAPE




  SUBROUTINE DETERMIN_SLOCLIST( CV_SLOCLIST, CV_NLOC, CV_SNLOC, SCVNGI, CV_ON_FACE, NFACE )
    ! determine CV_SLOCLIST
    IMPLICIT NONE
    INTEGER, intent( in ) :: CV_NLOC, CV_SNLOC, SCVNGI, NFACE
    INTEGER, DIMENSION( NFACE, CV_SNLOC ), intent( inout ) :: CV_SLOCLIST
    LOGICAL, DIMENSION( CV_NLOC, SCVNGI ), intent( in ) :: CV_ON_FACE
    ! Local variables
    INTEGER, DIMENSION( : ), allocatable :: CV_SLOC2LOC
    LOGICAL, DIMENSION( : ), allocatable :: GI_ON_ELE_BOUND
    INTEGER :: NUM_FACE, GI, CV_SKLOC, CV_KLOC, CV_SILOC, CV_SJLOC, CV_ILOC, CV_JLOC, IFACE
    LOGICAL :: FOUND, FOUND_ALL, NEW_FACE

    ALLOCATE( CV_SLOC2LOC( CV_SNLOC ))
    !write(357,*)'**********************************CV_SNLOC=',CV_SNLOC
    !write(357,*)'**********************************CV_SNLOC=',CV_SNLOC
    !write(357,*)'**********************************CV_SNLOC=',CV_SNLOC

    ALLOCATE( GI_ON_ELE_BOUND( SCVNGI ))
    DO GI=1,SCVNGI
       GI_ON_ELE_BOUND(GI)=.FALSE.
       DO CV_KLOC = 1, CV_NLOC
  IF(CV_ON_FACE(CV_KLOC,GI)) GI_ON_ELE_BOUND(GI)=.NOT.GI_ON_ELE_BOUND(GI)
       END DO
    END DO

    write(357,*)'GI_ON_ELE_BOUND:',GI_ON_ELE_BOUND

    !  write(357,*) 'cvsnloc:', cv_snloc
    !do cv_kloc=1, cv_nloc
    !   WRITE(357,*) 'CV_ON_FACE:',(CV_ON_FACE(cv_kloc,gi),gi=1,scvngi)
       !     write(357,*) 'CV_SLOCLIST:',(CV_SLOCLIST(iface,cv_kloc), iface=1,nface)
    !end do

    NUM_FACE = 0
    Loop_GI: DO GI = 1, SCVNGI
       IF(GI_ON_ELE_BOUND(GI)) THEN
      !WRITE(357,*)'*************GI=',GI

          CV_SKLOC = 0
          DO CV_KLOC = 1, CV_NLOC
             WRITE(357,*)'GI,CV_KLOC,CV_ON_FACE( CV_KLOC, GI ):',GI,CV_KLOC,CV_ON_FACE( CV_KLOC, GI )
             IF( CV_ON_FACE( CV_KLOC, GI )) THEN
                CV_SKLOC = CV_SKLOC + 1
                CV_SLOC2LOC( CV_SKLOC ) = CV_KLOC
             END IF
          END DO
          !write(357,*)'num_face,gi=',num_face,gi
          !write(357,*)'CV_SLOC2LOC:',CV_SLOC2LOC
          !WRITE(357,*) 'CV_SLOC2LOC --',( CV_SKLOC /= 0 ), (CV_SLOC2LOC( CV_KLOC ), CV_KLOC=1,CV_SKLOC)
          !do iface=1,num_face
          !   WRITE(357,*) 'gi, iface, CV_SLOCList', gi, iface, (CV_SLOCLIST(iface, cv_kloc),cv_kloc=1,cv_snloc)
          !end do

          Conditional_CV_SKLOC: IF( CV_SKLOC /= 0 ) THEN ! is this a new face
             NEW_FACE = .TRUE.

             Loop_IFACE: DO IFACE = 1, NUM_FACE
                FOUND_ALL = .TRUE.
                !write(357,*) '===> iface:',iface

                Loop_CVLOC1: DO CV_SILOC = 1, CV_SNLOC
                   CV_ILOC = CV_SLOC2LOC( CV_SILOC )
                   FOUND = .FALSE.

                   Loop_CVLOC2: DO CV_SJLOC = 1, CV_SNLOC
                      CV_JLOC = CV_SLOCLIST( IFACE, CV_SJLOC )
                      !if( (gi == scvngi) .and. (iface == nface))write(357,*) 'here',cv_jloc
                      IF( CV_ILOC == CV_JLOC ) FOUND = .TRUE.
                      !write(357,*) 'from found:',found, iface, cv_sjloc, cv_sloclist(iface,cv_sjloc)
                   END DO Loop_CVLOC2

                   IF( .NOT. FOUND ) FOUND_ALL = .FALSE.

                END DO Loop_CVLOC1

                IF( FOUND_ALL ) NEW_FACE = .FALSE.
             END DO Loop_IFACE

     IF(CV_NLOC==1) NEW_FACE = .TRUE.

             !WRITE(357,*) '******new_face:', new_face
             !WRITE(357,*) '******CV_SLOC2LOC:',CV_SLOC2LOC
             IF( NEW_FACE ) THEN
                NUM_FACE = NUM_FACE + 1
                CV_SLOCLIST( NUM_FACE, : ) = CV_SLOC2LOC( : )
                !WRITE(357,*) 'NUM FACE:', NUM_FACE
             END IF
          END IF Conditional_CV_SKLOC

       ENDIF
    END DO Loop_GI

    IF(NUM_FACE /= NFACE) THEN
       WRITE(357,*)' Not got the number of faces correct ' 
       WRITE(357,*)' NUM_FACE, NFACE :', NUM_FACE, NFACE
       STOP 962
    END IF

    DEALLOCATE( CV_SLOC2LOC )
    DEALLOCATE( GI_ON_ELE_BOUND )

    RETURN
  END SUBROUTINE DETERMIN_SLOCLIST




  SUBROUTINE SHAPESV_FEM_PLUS( SCVNGI, CV_NEILOC, CV_ON_FACE, &  
       CV_ELE_TYPE, CV_NLOC, SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
       SCVFENLX, SCVFENLY, SCVFENLZ,  &
       U_NLOC,  SUFEN, SUFENSLX, SUFENSLY,  &
       SUFENLX, SUFENLY, SUFENLZ,  &
       NDIM )
    !     
    !     - this subroutine generates the FE basis functions, weights and the
    !     - derivatives of the shape functions for a variety of elements on the 
    !     - control volume boundaries.  
    !     - The routine also generates the shape functions and derivatives 
    !     - associated with the CV surfaces and also the FV basis functions.
    !     -------------------------------
    !     - date last modified : 24/05/2003
    !     -------------------------------

    IMPLICIT NONE

    INTEGER, intent( in ) :: SCVNGI,CV_ELE_TYPE,CV_NLOC,U_NLOC,NDIM
    INTEGER, DIMENSION( CV_NLOC, SCVNGI ), intent( inout ) :: CV_NEILOC
    ! CV_ON_FACE(CV_KLOC,GI)=.TRUE. if CV_KLOC is on the face that GI is centred on.
    LOGICAL, DIMENSION( CV_NLOC, SCVNGI ), intent( inout ) :: CV_ON_FACE
    REAL, DIMENSION( CV_NLOC, SCVNGI ), intent( inout ) :: SCVFEN, SCVFENSLX, SCVFENSLY, &       
         SCVFENLX, SCVFENLY, SCVFENLZ
    REAL, DIMENSION( SCVNGI ), intent( inout ) :: SCVFEWEIGH
    REAL, DIMENSION( U_NLOC, SCVNGI ), intent( inout ) :: SUFEN, SUFENSLX, SUFENSLY, &       
         SUFENLX, SUFENLY, SUFENLZ
    ! Local variables
    INTEGER, PARAMETER :: TWO = 2, THREE = 3
    REAL, DIMENSION( :,: ), allocatable :: M
    REAL, DIMENSION( : ), allocatable :: DXN!, LX, XN, 
    INTEGER :: ILOC

    !write(357,*) 'b4 hhhh', two, three, CV_NLOC,SCVNGI
    !ALLOCATE( LX( TWO ))
    !ALLOCATE( XN( THREE ))
    ALLOCATE( DXN( THREE ))
    ALLOCATE( M( CV_NLOC, SCVNGI ))

    Cond_ShapeType: IF(CV_ELE_TYPE == 1) THEN ! calculate the shape functions on the control volume boundaries

       Cond_Dimension: IF( NDIM == 1 ) THEN

          CALL FV_1D_QUAD(SCVNGI, CV_NLOC, SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
               SCVFENLX, SCVFENLY, SCVFENLZ) ! For T

          CALL FV_1D_QUAD(SCVNGI, U_NLOC, SUFEN, SUFENSLX, SUFENSLY, SCVFEWEIGH, &
               SUFENLX, SUFENLY, SUFENLZ) ! For U

          !           STOP 37

          CV_ON_FACE=.FALSE.
          DO ILOC = 1, CV_NLOC
             CV_ON_FACE( ILOC, ILOC ) = .TRUE.
             CV_ON_FACE( ILOC, ILOC + 1 ) =.TRUE.
          END DO

          !          U_ON_FACE=.FALSE.
          !          DO ILOC = 1, U_NLOC
          !             U_ON_FACE( ILOC, ILOC ) = .TRUE.
          !             U_ON_FACE( ILOC, ILOC + 1 ) =.TRUE.
          !          END DO

       ELSE


          !     
          IF(CV_ELE_TYPE == 1) THEN
             !     
             CALL FVQUAD( SCVNGI,    CV_NLOC, SCVNGI,&
                  M,      SCVFEN,  SCVFENSLX,&
                  SCVFEWEIGH              )
             !     
          ELSE IF(CV_ELE_TYPE == 2) THEN
             !     
             CALL FVTRI( SCVNGI,    CV_NLOC, SCVNGI,&
                  M,      SCVFEN,  SCVFENSLX,&
                  SCVFEWEIGH              )
             !     
          ELSE IF(CV_ELE_TYPE == 3) THEN
             !     
             CALL FVHEX( SCVNGI,    CV_NLOC,    SCVNGI,&
                  M,      SCVFEN,     SCVFENSLX, &
                  SCVFENSLY,  SCVFEWEIGH         )
             !     
          ELSE IF(CV_ELE_TYPE == 4) THEN
             !     
             CALL FVTET( SCVNGI,    CV_NLOC,    SCVNGI,&
                  M,      SCVFEN,     SCVFENSLX, &
                  SCVFENSLY,  SCVFEWEIGH         )
             !     
          ELSE IF(CV_ELE_TYPE == 5) THEN
             !     
             CALL FVQQUAD(  SCVNGI,    CV_NLOC,    SCVNGI,&
                  M,      SCVFEN,     SCVFENSLX, &
                  SCVFEWEIGH                 )
             !     
          ELSE IF(CV_ELE_TYPE == 6) THEN
             !     
             CALL FVQHEX( SCVNGI,   CV_NLOC,   SCVNGI,&
                  M,     SCVFEN,    SCVFENSLX,&
                  SCVFENSLY, SCVFEWEIGH        )
             !     
          ELSE IF(CV_ELE_TYPE == 7) THEN
             !     
             !     - quadratic triangles
             !     
             CALL FVQTRI( SCVNGI,    CV_NLOC, SCVNGI,&
                  M,      SCVFEN,  SCVFENSLX,  &
                  SCVFEWEIGH              )
             !     
          ELSE IF(CV_ELE_TYPE == 8) THEN
             !     
             CALL FVQTET( SCVNGI,   CV_NLOC,   SCVNGI,&
                  M,     SCVFEN,    SCVFENSLX,&
                  SCVFENSLY, SCVFEWEIGH        )
             !     
          END IF
          !     

       ENDIF Cond_Dimension

    ELSE

       WRITE(357,*)'NOT GOT THIS YET - 0'

    END IF Cond_ShapeType
    !     
    !     
    CALL VOLNEI( CV_NEILOC, CV_NLOC, SCVNGI, CV_ELE_TYPE )
    !     
    !DEALLOCATE( LX )
    !DEALLOCATE( XN )
    DEALLOCATE( DXN )
    DEALLOCATE( M )

    RETURN
  END SUBROUTINE SHAPESV_FEM_PLUS




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

    !write(357,*) 'cv_nloc:',cv_nloc

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
       WRITE(357,*)'NOYT AN OPT'
       STOP 38
    ENDIF

    DIFF=.TRUE.
    NDIFF=.FALSE.

    GETNDP=.TRUE.
    !      Compute standard Gauss quadrature. weits and points
    CALL LAGROT(WEI,CV_NODPOS,CV_NLOC,GETNDP) 
    WRITE(357,*)'CV_NODPOS:',CV_NODPOS
    !       STOP 29

    Loop_P2: DO GPOI = 1, NCV_BOU

       Loop_ILX2: DO  ILOC = 1, CV_NLOC 
          SCVFEN( ILOC, GPOI )  = LAGRAN(NDIFF,LX(GPOI),ILOC,CV_NLOC,CV_NODPOS)
          SCVFENSLX(ILOC, GPOI )=1.0
          SCVFENSLY(ILOC, GPOI )=0.0
          SCVFENLX( ILOC, GPOI ) = LAGRAN(DIFF,LX(GPOI),ILOC,CV_NLOC,CV_NODPOS)
          SCVFENLY( ILOC, GPOI ) = 0.0
          SCVFENLZ( ILOC, GPOI ) = 0.0
       END DO Loop_ILX2

    end do Loop_P2

    DEALLOCATE( LX )
    DEALLOCATE( CV_NODPOS )
    DEALLOCATE( WEI )

  END SUBROUTINE FV_1D_QUAD




  SUBROUTINE SHAPESV( ELETYP,  NEILOC,    &
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
    INTEGER, DIMENSION( NLOC, SVNGI ), intent( inout ) :: NEILOC
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

       WRITE(357,*)'NOT GOT THIS YET'

    END IF Cond_ShapeType
    !     
    !     
    CALL VOLNEI( NEILOC, NLOC, SVNGI, ELETYP )
    !     
    DEALLOCATE( LX )
    DEALLOCATE( XN )
    DEALLOCATE( DXN )
    !     
    RETURN
  END SUBROUTINE SHAPESV



  SUBROUTINE FVQUAD( NGI,    NLOC,    SVNGI,&
                                !     - REALS
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




  SUBROUTINE VOLNEI( NEILOC, NLOC, SVNGI, ELETYP )
    !--------------------------------------------------------   
    !- this subroutine calculates NEILOC which is the       
    !- array containing information given a local node
    !- and an integration point what is the other opposing 
    !- local node. 
    !- It contains -1 if on the boundary of the element     
    !--------------------------------------------------------
    !- date last modified : 11/10/2002
    !--------------------------------------------------------

    IMPLICIT NONE     
    INTEGER, intent( in ) :: ELETYP, NLOC, SVNGI    
    INTEGER, DIMENSION( NLOC, SVNGI ), intent( inout ) :: NEILOC
    ! Local variables    
    INTEGER :: ILOC

    NEILOC = 0

    IF(ELETYP == 1) THEN ! quadrilaterals

       NEILOC=0
       DO ILOC=1,NLOC
          NEILOC(ILOC,ILOC)  =ILOC-1
          NEILOC(ILOC,ILOC+1)=ILOC+1
       END DO
       NEILOC(1,1)=-1
       NEILOC(NLOC,SVNGI)=-1

    ELSE 

       write(357,*)'not an option'
       stop 292

    ENDIF


    RETURN     
  END SUBROUTINE VOLNEI




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
          !      write(357,*)'iloc,gi,NEILOC( ILOC, GI ):',iloc,gi,NEILOC( ILOC, GI )

       END DO

    END DO

    FINDGPTS( NLOC + 1 ) = COUNT + 1
    NCOLGPTS = COUNT
    !    stop 2821

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
            +CGI * ( DGI * HGI - EGI *GGI )

       DETWEI( GI ) = abs( DETJ ) * WEIGHT( GI )
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





  SUBROUTINE DETNLXR( ELE, X,Y,Z, XONDGL, TOTELE, NONODS, NLOC, NGI, &
       N, NLX, NLY, NLZ, WEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
       NX, NY, NZ) 
    IMPLICIT NONE
    INTEGER, intent( in ) :: ELE, TOTELE, NONODS, NLOC, NGI
    INTEGER, DIMENSION( TOTELE * NLOC ) :: XONDGL
    REAL, DIMENSION( NONODS ), intent( in ) :: X, Y, Z
    REAL, DIMENSION( NLOC, NGI ), intent( in ) :: N, NLX, NLY, NLZ 
    REAL, DIMENSION( NGI ), intent( in ) :: WEIGHT
    REAL, DIMENSION( NGI ), intent( inout ) :: DETWEI, RA
    REAL, intent( inout ) :: VOLUME
    LOGICAL, intent( in ) :: D1, D3, DCYL
    REAL, DIMENSION( NLOC, NGI ), intent( inout ) :: NX, NY, NZ
    ! Local variables
    REAL, PARAMETER :: PIE = 3.141592654
    REAL :: AGI, BGI, CGI, DGI, EGI, FGI, GGI, HGI, KGI, A11, A12, A13, A21, &
         A22, A23, A31, A32, A33, DETJ, TWOPIE, RGI
    INTEGER :: GI, L, IGLX
    !
    VOLUME = 0.
    !
    IF(D3) THEN
       do  GI=1,NGI! Was loop 331
          !
          AGI=0.
          BGI=0.
          CGI=0.
          !
          DGI=0.
          EGI=0.
          FGI=0.
          !
          GGI=0.
          HGI=0.
          KGI=0.
          !
          do  L=1,NLOC! Was loop 79
             IGLX=XONDGL((ELE-1)*NLOC+L)
             ! NB R0 does not appear here although the z-coord might be Z+R0. 
             AGI=AGI+NLX(L,GI)*X(IGLX) 
             BGI=BGI+NLX(L,GI)*Y(IGLX) 
             CGI=CGI+NLX(L,GI)*Z(IGLX) 
             !
             DGI=DGI+NLY(L,GI)*X(IGLX) 
             EGI=EGI+NLY(L,GI)*Y(IGLX) 
             FGI=FGI+NLY(L,GI)*Z(IGLX) 
             !
             GGI=GGI+NLZ(L,GI)*X(IGLX) 
             HGI=HGI+NLZ(L,GI)*Y(IGLX) 
             KGI=KGI+NLZ(L,GI)*Z(IGLX)
          end do ! Was loop 79
          !
          DETJ=AGI*(EGI*KGI-FGI*HGI)&
               -BGI*(DGI*KGI-FGI*GGI)&
               +CGI*(DGI*HGI-EGI*GGI)
          DETWEI(GI)=DETJ*WEIGHT(GI)
          RA(GI)=1.0
          VOLUME=VOLUME+DETWEI(GI)
          ! For coefficient in the inverse mat of the jacobian. 
          A11= (EGI*KGI-FGI*HGI) /DETJ
          A21=-(DGI*KGI-FGI*GGI) /DETJ
          A31= (DGI*HGI-EGI*GGI) /DETJ
          !
          A12=-(BGI*KGI-CGI*HGI) /DETJ
          A22= (AGI*KGI-CGI*GGI) /DETJ
          A32=-(AGI*HGI-BGI*GGI) /DETJ
          !
          A13= (BGI*FGI-CGI*EGI) /DETJ
          A23=-(AGI*FGI-CGI*DGI) /DETJ
          A33= (AGI*EGI-BGI*DGI) /DETJ
          do  L=1,NLOC! Was loop 373
             NX(L,GI)= A11*NLX(L,GI)+A12*NLY(L,GI)+A13*NLZ(L,GI)
             NY(L,GI)= A21*NLX(L,GI)+A22*NLY(L,GI)+A23*NLZ(L,GI)
             NZ(L,GI)= A31*NLX(L,GI)+A32*NLY(L,GI)+A33*NLZ(L,GI)
          end do ! Was loop 373
          !
       end do ! Was loop 331
       ! IF(D3) THEN...
    ELSE IF(.NOT.D1) THEN
       TWOPIE=1.0 
       IF(DCYL) TWOPIE=2.*PIE
       do  GI=1,NGI! Was loop 1331
          !
          RGI=0.
          !
          AGI=0.
          BGI=0.
          CGI=0.
          DGI=0.
          !
          do  L=1,NLOC! Was loop 179
             IGLX=XONDGL((ELE-1)*NLOC+L)

             AGI=AGI + NLX(L,GI)*X(IGLX) 
             BGI=BGI + NLX(L,GI)*Y(IGLX) 
             CGI=CGI + NLY(L,GI)*X(IGLX) 
             DGI=DGI + NLY(L,GI)*Y(IGLX) 
             !
             RGI=RGI+N(L,GI)*Y(IGLX)
          end do ! Was loop 179
          !
          IF(.NOT.DCYL) RGI=1.0
          !
          DETJ= AGI*DGI-BGI*CGI 
          RA(GI)=RGI
          DETWEI(GI)=TWOPIE*RGI*DETJ*WEIGHT(GI)
          VOLUME=VOLUME+DETWEI(GI)
          !
          do L=1,NLOC
             NX(L,GI)=(DGI*NLX(L,GI)-BGI*NLY(L,GI))/DETJ
             NY(L,GI)=(-CGI*NLX(L,GI)+AGI*NLY(L,GI))/DETJ
             NZ(L,GI)=0.0
          END DO
          !
       end do ! Was loop 1331
       ! ENDOF IF(D3) THEN ELSE...
    ELSE 
       ! For 1D...
       do  GI = 1, NGI
          !
          AGI = 0.
          !
          do  L = 1, NLOC
             IGLX = XONDGL(( ELE - 1 ) * NLOC + L )
             AGI = AGI + NLX( L, GI ) * X( IGLX ) 
          end do
          !
          DETJ = AGI 
          DETWEI( GI ) = DETJ * WEIGHT( GI )
          VOLUME = VOLUME + DETWEI( GI )
          !
          do L = 1, NLOC
             NX( L, GI ) = NLX( L, GI ) / DETJ
             NY( L, GI ) = 0.0
             NZ( L, GI ) = 0.0
          END DO
          !
       end do 
       ! ENDOF IF(D3) THEN ELSE...
    ENDIF
    !
    RETURN
  END SUBROUTINE DETNLXR




  SUBROUTINE DETNLXR_PLUS_U( ELE, X, Y, Z, XONDGL, TOTELE, NONODS, X_NLOC, NGI, &
       N, NLX, NLY, NLZ, WEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
       NX, NY, NZ, &
       U_NLOC, UNLX, UNLY, UNLZ, UNX, UNY, UNZ ) 
    implicit none
    INTEGER, intent( in ) :: ELE, TOTELE, NONODS, X_NLOC, NGI, U_NLOC
    INTEGER, DIMENSION( TOTELE * X_NLOC ) :: XONDGL
    REAL, DIMENSION( NONODS ), intent( in ) :: X, Y, Z
    REAL, DIMENSION( X_NLOC, NGI ), intent( in ) :: N, NLX, NLY, NLZ 
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
         A22, A23, A31, A32, A33, DETJ, TWOPIE, RGI
    INTEGER :: GI, L, IGLX

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
          DETWEI( GI ) = DETJ * WEIGHT( GI )
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
          write(357,*) 'gi:', gi,detj, detwei(gi)

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
    !
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

       detj =  agi * (egi * kgi - fgi * hgi) &
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

