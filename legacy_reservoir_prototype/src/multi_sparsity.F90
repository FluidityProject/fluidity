  
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
  
  module spact

    use fldebug
    integer, parameter :: nface_p1 = 3, internal_faces = 1

  contains


    SUBROUTINE GET_SPARS_PATS(   &
         NDIM, U_PHA_NONODS, CV_PHA_NONODS,   &
         U_NONODS, CV_NONODS, X_NONODS, &
         U_NLOC, CV_NLOC, X_NLOC, U_SNLOC, CV_SNLOC, X_SNLOC, NPHASE, TOTELE, &
         U_NDGLN, CV_NDGLN, X_NDGLN, &
                                ! CV multi-phase eqns (e.g. vol frac, temp)..***************
         MX_NCOLACV, NCOLACV, FINACV, COLACV, MIDACV,  &
                                ! Force balance plus cty multi-phase eqns....***************
         NLENMCY, MX_NCOLMCY, NCOLMCY,FINMCY, COLMCY, MIDMCY, &
                                ! Element connectivity...
         MXNELE, NCOLELE, MIDELE, FINELE, COLELE, &
                                ! Force balance sparsity...
         MX_NCOLDGM_PHA, NCOLDGM_PHA, COLDGM_PHA, FINDGM_PHA, MIDDGM_PHA, &
                                ! CT sparsity - global cty eqn...
         MX_NCT, NCT, FINDCT, COLCT, &
                                ! C sparsity operating on pressure in force balance...
         MX_NC, NC, FINDC, COLC, &
                                ! pressure matrix for projection method...
         MX_NCOLCMC, NCOLCMC, FINDCMC, COLCMC, MIDCMC, &
                                ! CV-FEM matrix...
         MX_NCOLM, NCOLM, FINDM, COLM, MIDM, U_ELE_TYPE )
      ! Obtain the sparsity patterns of the two types of 
      ! matricies for (momentum + cty) and for energy

      IMPLICIT NONE

      INTEGER, intent( in ) :: NDIM, U_PHA_NONODS, CV_PHA_NONODS, U_NONODS, CV_NONODS, X_NONODS, U_NLOC, CV_NLOC, &
           X_NLOC, U_SNLOC, CV_SNLOC, X_SNLOC, NPHASE, TOTELE, U_ELE_TYPE
      INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in )     :: U_NDGLN
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in )    :: CV_NDGLN
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in )    :: X_NDGLN

      INTEGER, intent( in )                                   :: MX_NCOLACV 
      INTEGER, intent( inout )                                :: NCOLACV

      INTEGER, DIMENSION( CV_PHA_NONODS + 1), intent( inout ) :: FINACV ! CV multi-phase eqns (e.g. vol frac, temp)
      INTEGER, DIMENSION( MX_NCOLACV ), intent( inout )       :: COLACV
      INTEGER, DIMENSION( CV_PHA_NONODS ), intent( inout )    :: MIDACV

      INTEGER, intent( in )                                   :: NLENMCY, MX_NCOLMCY
      INTEGER, intent( inout )                                :: NCOLMCY

      INTEGER, DIMENSION( NLENMCY + 1 ), intent( inout )      :: FINMCY ! Force balance plus cty multi-phase eqns
      INTEGER, DIMENSION( MX_NCOLMCY ), intent( inout )       :: COLMCY
      INTEGER, DIMENSION( NLENMCY ), intent( inout )          :: MIDMCY

      INTEGER, intent( in )                                   :: MXNELE
      INTEGER, intent( inout )                                :: NCOLELE

      INTEGER, DIMENSION( TOTELE + 1), intent( inout )        :: FINELE ! Element connectivity...
      INTEGER, DIMENSION( TOTELE ), intent( inout )           :: MIDELE
      INTEGER, DIMENSION( MXNELE ), intent( inout )           :: COLELE

      INTEGER, intent( in )                                   :: MX_NCOLDGM_PHA
      INTEGER, intent( inout )                                :: NCOLDGM_PHA

      INTEGER, DIMENSION( MX_NCOLDGM_PHA ), intent( inout )   :: COLDGM_PHA ! Force balance sparsity
      INTEGER, DIMENSION( U_PHA_NONODS + 1 ), intent( inout ) :: FINDGM_PHA
      INTEGER, DIMENSION( U_PHA_NONODS ), intent( inout )     :: MIDDGM_PHA

      INTEGER, intent( in )                                    :: MX_NC, MX_NCT
      INTEGER, intent( inout )                                 :: NC, NCT

      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( inout )     :: FINDCT
      INTEGER, DIMENSION( MX_NCT ), intent( inout )            :: COLCT !operating on pressure in force balance
      INTEGER, DIMENSION( U_NONODS + 1 ), intent( inout )      :: FINDC
      INTEGER, DIMENSION( MX_NC ), intent( inout )             :: COLC

      INTEGER, intent( in )                                   :: MX_NCOLCMC
      INTEGER, intent( inout )                                :: NCOLCMC
      INTEGER, DIMENSION( CV_NONODS+1 ), intent( inout )      :: FINDCMC
      INTEGER, DIMENSION( MX_NCOLCMC ), intent( inout )       :: COLCMC
      INTEGER, DIMENSION( CV_NONODS ), intent( inout )        :: MIDCMC

      INTEGER, intent( in )                                   :: MX_NCOLM
      INTEGER, intent( inout )                                :: NCOLM
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( inout )    :: FINDM
      INTEGER, DIMENSION( NCOLM ), intent( inout )            :: COLM
      INTEGER, DIMENSION( CV_NONODS ), intent( inout )        :: MIDM

      ! Local variables
      INTEGER :: NPHA_TOTELE, MXNACV_LOC, NACV_LOC, MX_NCOLELE_PHA, II

      INTEGER, ALLOCATABLE, DIMENSION( : ) :: MIDACV_LOC, FINACV_LOC, COLACV_LOC, COLELE_PHA, FINELE_PHA, &
           MIDELE_PHA, CENTCT
      LOGICAL, parameter :: OldSetUp = .false.
      logical :: presym
      integer, dimension( : ), allocatable :: tempvec1, tempvec2, tempvec3, dummyvec
      integer :: tempnct

      ewrite(3,*) 'In GET_SPARS_PATS'

      ! Extend the element sparsity to an element multiphase sparsity

      MX_NCOLELE_PHA = NPHASE * MXNELE + ( NPHASE - 1 ) * NPHASE * TOTELE
      NPHA_TOTELE = TOTELE * NPHASE

      ALLOCATE( COLELE_PHA( MX_NCOLELE_PHA ))
      ALLOCATE( FINELE_PHA( NPHA_TOTELE + 1 ))
      ALLOCATE( MIDELE_PHA( NPHA_TOTELE ))
      ALLOCATE( CENTCT( CV_NONODS ))

      if ( OldSetUp ) then
         CALL DEF_SPAR( 1, TOTELE, MXNELE, NCOLELE, &
              MIDELE, FINELE, COLELE )
         !  ewrite(3,*) colele( 1 : ncolele )
         !  ewrite(3,*) midele( 1 : totele )
         !  ewrite(3,*) ' '
         !  allocate( tempvec1( totele + 1 ))
         !  allocate( tempvec2( totele ))
         !  allocate( tempvec3( ncolele ))
         !  tempvec1 = finele
         !  tempvec2 = midele
         !  tempvec3( 1: ncolele ) = colele( 1: ncolele )
         !  finele = 0  
         !  midele = 0
         !  colele = 0
      else
         call getfinele( totele, cv_nloc, x_snloc, x_nonods, x_ndgln, nface_p1, &
              finele, colele, midele )
         ncolele = nface_p1 * totele
         ! ewrite(3,*) colele( 1 : ncolele )
         ! ewrite(3,*) midele( 1 : totele )
         ! ewrite(3,*) ncolele, totele * nface_p1, mxnele
         ! ewrite(3,*)'  -----  finele  -----  '
         ! call comparematrices( totele + 1,        tempvec1, finele )
         ! ewrite(3,*)' -----  midele -----  '
         ! call comparematrices( totele,            tempvec2, midele )
         ! ewrite(3,*)' -----  colele  -----  '
         ! call comparematrices( ncolele, tempvec3, colele )
      endif
      ! stop 11

      CALL EXTEN_SPARSE_MULTI_PHASE( TOTELE, MXNELE, FINELE, COLELE, &
           NPHASE, NPHA_TOTELE, MX_NCOLELE_PHA, FINELE_PHA, COLELE_PHA, MIDELE_PHA ) 


      CALL FORM_DGM_PHA_SPARSITY( U_PHA_NONODS, MX_NCOLDGM_PHA, NCOLDGM_PHA, COLDGM_PHA, FINDGM_PHA, MIDDGM_PHA, &
           TOTELE, NPHASE, U_NLOC, MX_NCOLELE_PHA, FINELE_PHA, COLELE_PHA ) 


      ! Now form the global matrix: FINMCY, COLMCY and MIDMCY for the 
      ! momentum and continuity eqns: 

      if ( OldSetUp .or. ( cv_nonods /= totele * cv_nloc )) then
         CALL DEF_SPAR_CT_DG( CV_NONODS, MX_NCT, NCT, FINDCT, COLCT, TOTELE, CV_NLOC, U_NLOC, U_NDGLN, U_ELE_TYPE)
         NC = NCT
         !   call swapprint( .true., cv_nonods, mx_nct, nct, findct, colct, centct )
      else
         call pousinmc2( totele, u_nonods, u_nloc, cv_nonods, cv_nloc, mx_nct, u_ndgln, cv_ndgln, &
              nct, findct, colct, centct )
         nc = nct
      endif

      ! Convert CT sparsity to C sparsity.
      CALL CONV_CT2C( CV_NONODS, NCT, FINDCT, COLCT, U_NONODS, MX_NC, FINDC, COLC )

      ewrite(3,*) 'NCOLDGM_PHA,CV_NONODS, TOTELE*CV_NLOC:',NCOLDGM_PHA,CV_NONODS, TOTELE*CV_NLOC

      IF(CV_NONODS /= TOTELE*CV_NLOC) THEN ! Continuous pressure

         if ( .not. OldSetUp ) then
            CALL DEF_SPAR( CV_NLOC - 1 , CV_NONODS, MX_NCOLCMC, NCOLCMC, &
                                !    CALL DEF_SPAR( CV_NLOC+2, CV_NONODS, MX_NCOLCMC, NCOLCMC, &
                 MIDCMC, FINDCMC, COLCMC )
         else
            allocate( dummyvec( u_nonods ))
            dummyvec = 0
            presym = .false.
            call poscmc( cv_nonods, u_nonods, mx_ncolcmc, nct, &
                 findct, colct, &
                 ncolcmc, findcmc, colcmc, midcmc, dummyvec, presym )
            deallocate( dummyvec )
         end if

      ELSE
         if ( OldSetUp ) then
            CALL DEF_SPAR( CV_NLOC+2, CV_NONODS, MX_NCOLCMC, NCOLCMC, &
                 MIDCMC, FINDCMC, COLCMC )
            ewrite(3,*)'findcmc: ', findcmc(1:cv_nonods+1)
            ewrite(3,*)'colcmc: ', colcmc( 1:ncolcmc )
            ewrite(3,*)'midcmc: ', midcmc( 1:cv_nonods)
         else
            allocate(tempvec1(cv_nonods+1))
            allocate(tempvec2(mx_ncolcmc))
            allocate(tempvec3(cv_nonods))
            tempvec1 = findcmc
            tempvec2 = colcmc
            tempvec3 = midcmc
            findcmc = 0 ; colcmc = 0 ; midcmc = 0
            allocate( dummyvec( u_nonods ))
            presym = .false.
            call poscmc( cv_nonods, u_nonods, mx_ncolcmc, nct, &
                 findct, colct, &
                 ncolcmc, findcmc, colcmc, midcmc, dummyvec, presym )
            ewrite(3,*)'findcmc: ', findcmc(1:cv_nonods+1)
            ewrite(3,*)'colcmc: ', colcmc(1:mx_ncolcmc )
            ewrite(3,*)'midcmc: ', midcmc(1:cv_nonods)
            deallocate( dummyvec )
            deallocate( tempvec1 )
            deallocate( tempvec2 )
            deallocate( tempvec3 )

         end if
         !stop 98
      ENDIF
      ewrite( 3, *)'CV_NONODS /= TOTELE*CV_NLOC: ', CV_NONODS /= TOTELE*CV_NLOC
      !stop 98
      if(MX_NCOLCMC < NCOLCMC) FLAbort(" Incorrect number of computed dimension of CMC sparcity matrix")

      ewrite(3,*)'defining sparsity of matrix'

      MIDMCY = 0
      ! Extend momentum sparsity to include Continuity & Pressure: 
      CALL EXTEN_SPARSE_MOM_CTY( NDIM, FINDGM_PHA, COLDGM_PHA, U_PHA_NONODS, NCOLDGM_PHA, NCT, &
           CV_NONODS, FINDCT, COLCT, &
           U_NONODS, NC, MX_NCOLMCY, &
           FINDC, COLC, FINMCY, COLMCY, MIDMCY, NLENMCY, &
           NCOLMCY, NPHASE, NCOLCMC, FINDCMC, COLCMC  )  
      !
      !        MX_NMCY
      ! Sparsity for energy equation...   
      ! NACV_LOC:                            
      MXNACV_LOC = 3 * CV_NONODS
      ALLOCATE( MIDACV_LOC( CV_NONODS ))
      ALLOCATE( FINACV_LOC( CV_NONODS + 1 ))
      ALLOCATE( COLACV_LOC( MXNACV_LOC ))

   !!   if( .not. OldSetUp ) then
         CALL DEF_SPAR( 1 , CV_NONODS, MXNACV_LOC, NACV_LOC, &
              MIDACV_LOC, FINACV_LOC, COLACV_LOC )

         NCOLACV = NPHASE * NACV_LOC + ( NPHASE - 1 ) * NPHASE * CV_NONODS
         NACV_LOC = NCOLACV
         ewrite(3,*)'totele,cv_nonods,cv_nloc:', totele,cv_nonods,cv_nloc
         ewrite(3,*)'finacv_loc:', size( finacv_loc ), finacv_loc( 1 : cv_nonods + 1 )
         ewrite(3,*)'colacv_loc:', size( colacv_loc ), colacv_loc( 1 : mxnacv_loc )
         ewrite(3,*)'midacv_loc:', size( midacv_loc ), midacv_loc( 1 : cv_nonods )

         !stop 98874
    !!     allocate( tempvec1( cv_nonods + 1 ))
     !!    allocate( tempvec2( mxnacv_loc ))
     !!    allocate( tempvec3( cv_nonods ))
     !!    tempvec1( 1 : cv_nonods + 1 ) = finacv_loc( 1 : cv_nonods + 1 )
     !!    tempvec2( 1 : mxnacv_loc ) = colacv_loc( 1 : mxnacv_loc )
     !!    tempvec3( 1 : cv_nonods ) = midacv_loc( 1 : cv_nonods )
        !! call swapprint( .true., cv_nonods, mxnacv_loc, nacv_loc, finacv_loc, colacv_loc, midacv_loc )

         ! else
         !!nacv_loc = 0
         !!nct = 0
         !!call posinm_colour( totele, cv_nonods, cv_nloc, colacv_loc, nct, mxnacv_loc, finacv_loc, midacv_loc, &
        !!                 u_nonods, u_nloc, u_ndgln, cv_ndgln )
       !!  call pousinmc2( totele, cv_nonods, cv_nloc, u_nonods, u_nloc, mxnacv_loc, cv_ndgln, u_ndgln, &
       !!       nct, finacv_loc, colacv_loc, midacv_loc )
         !!nacv_loc = nct
         !!call swapprint( .false., cv_nonods, mxnacv_loc, nacv_loc, finacv_loc, colacv_loc, midacv_loc )
         !!call comparematrices( cv_nonods + 1, tempvec1, finacv_loc )
         !!call comparematrices( nacv_loc,      tempvec2, colacv_loc )
         !!call comparematrices( cv_nonods,     tempvec3, midacv_loc )
  !!    end if
  !!    stop 9887

      ewrite(3,*) 'going into EXTEN_SPARSE_MULTI_PHASE sbrt 2nd time'
      CALL EXTEN_SPARSE_MULTI_PHASE( CV_NONODS, MXNACV_LOC, FINACV_LOC, COLACV_LOC, &
           NPHASE, CV_PHA_NONODS, MX_NCOLACV, FINACV, COLACV, MIDACV )  

      ! Sparsity of CV-FEM
      ewrite(3,*)'going to find sparsity of colm MX_NCOLM:',MX_NCOLM
      ewrite(3,*)'CV_NONODS,CV_NLOC:',CV_NONODS,CV_NLOC
      CALL DEF_SPAR( CV_NLOC - 1 , CV_NONODS, MX_NCOLM, NCOLM, &
           MIDM, FINDM, COLM )

      ewrite(3,*) 'NCOLMCY, NCOLACV, NCOLCMC:', NCOLMCY, NCOLACV, NCOLCMC

      ewrite(3,*) 'Leaving GET_SPARS_PATS'

      RETURN

    END SUBROUTINE GET_SPARS_PATS

!!!
!!! Def_Spar Subroutine
!!!

    SUBROUTINE DEF_SPAR( SEMI_BAND_WID, NONODS, MX_NCOLC, NCOLC, &
         CENTC, FINDC, COLC)
      ! define sparsity...
      ! SEMI_BAND_WID is the semi band width.
      IMPLICIT NONE
      INTEGER, intent( in ) :: SEMI_BAND_WID, NONODS, MX_NCOLC
      INTEGER, intent( inout) :: NCOLC
      INTEGER, DIMENSION( NONODS ), intent( inout) ::  CENTC
      INTEGER, DIMENSION( NONODS + 1 ), intent( inout) :: FINDC
      INTEGER, DIMENSION( MX_NCOLC ), intent( inout) :: COLC

      ! Local variables
      INTEGER :: NOD, II, COL, COUNT

      ewrite(3,*) 'In DEF_SPAR'

      COUNT = 0
      COL = 0
      loop_nods: DO NOD = 1, NONODS

         FINDC( NOD ) = COUNT + 1

         DO II = -SEMI_BAND_WID, SEMI_BAND_WID, 1

            COL = NOD + II

            IF( ( COL >= 1 ) .AND. ( COL <= NONODS ) ) THEN
               COUNT = COUNT + 1
               COLC( COUNT ) = COL
               IF( COL == NOD) CENTC( NOD ) = COUNT
            END IF

         END DO

      END DO loop_nods

      FINDC( NONODS + 1 ) = COUNT + 1

      NCOLC = COUNT

      ewrite(3,*) 'Leaving DEF_SPAR'

      RETURN
    END SUBROUTINE DEF_SPAR


    SUBROUTINE EXTEN_SPARSE_MULTI_PHASE( NONODS, MXNELE, FINM, COLM, &
         NPHASE, NPHA_NONODS, NCOLM_PHA, FINM_PHA, COLM_PHA, MIDM_PHA ) 

      ! Extend the sparsity to a multiphase sparsity

      IMPLICIT NONE 
      INTEGER, intent( in ) :: NONODS, MXNELE, NPHASE
      INTEGER, DIMENSION( NONODS + 1), intent (in ) :: FINM
      INTEGER, DIMENSION( MXNELE ), intent ( in ) :: COLM
      INTEGER, intent( in ) :: NPHA_NONODS, NCOLM_PHA
      INTEGER, DIMENSION( NPHA_NONODS + 1 ), intent( inout ) :: FINM_PHA
      INTEGER, DIMENSION( NPHA_NONODS ), intent( inout ) :: MIDM_PHA
      INTEGER, DIMENSION( NCOLM_PHA ), intent( inout ) :: COLM_PHA

      INTEGER :: COUNT, COUNT2, IPHASE, JPHASE, NOD 

      ewrite(3,*) 'In EXTEN_SPARSE_MULTI_PHASE'

      COUNT2 = 0

      ewrite(3,*) ' ########  In  EXTEN_SPARSE_MULTI_PHASE Subrt.  ####### '

      loop_phase1: DO IPHASE = 1, NPHASE

         loop_CVNODS: DO NOD = 1, NONODS

            loop_phase2: DO JPHASE = 1, NPHASE

               IF( JPHASE == 1) FINM_PHA(( IPHASE - 1 ) * NONODS + NOD ) = COUNT2 + 1

               IF( IPHASE == JPHASE ) THEN
                  DO COUNT = FINM( NOD ), FINM( NOD + 1 ) - 1, 1
                     COUNT2 = COUNT2 + 1
                     COLM_PHA( COUNT2 ) = COLM( COUNT ) + ( JPHASE - 1 ) * NONODS
                     IF( COLM( COUNT ) ==  NOD ) MIDM_PHA( NOD + ( JPHASE - 1 ) * NONODS ) = COUNT2
                  END DO
               ELSE
                  COUNT2 = COUNT2 + 1
                  COLM_PHA( COUNT2 ) = NOD + ( JPHASE - 1 ) * NONODS
               ENDIF

            END DO loop_phase2

         END DO loop_CVNODS

      END DO loop_phase1

      FINM_PHA( NPHASE * NONODS + 1 ) = COUNT2 + 1

      ewrite(3,*) 'Leaving EXTEN_SPARSE_MULTI_PHASE'

      RETURN
    END SUBROUTINE EXTEN_SPARSE_MULTI_PHASE



    SUBROUTINE FORM_DGM_PHA_SPARSITY( U_PHA_NONODS, MX_NCOLDGM_PHA, NCOLDGM_PHA, COLDGM_PHA, FINDGM_PHA, MIDDGM_PHA, &
         TOTELE, NPHASE, U_NLOC, NELE_PHA, FINELE_PHA, COLELE_PHA ) 
      ! Form the sparsity of the phase coupled DG discretised matrix 
      ! from the element-wise multi-phase sparsity matrix. 
      IMPLICIT NONE

      INTEGER, intent( in ) :: U_PHA_NONODS, MX_NCOLDGM_PHA
      INTEGER, intent( inout ) :: NCOLDGM_PHA
      INTEGER, DIMENSION( MX_NCOLDGM_PHA ), intent ( inout ) :: COLDGM_PHA
      INTEGER, DIMENSION( U_PHA_NONODS + 1 ), intent( inout ) :: FINDGM_PHA
      INTEGER, DIMENSION( U_PHA_NONODS ), intent( inout ) :: MIDDGM_PHA
      INTEGER, intent( in ) :: TOTELE, NPHASE, U_NLOC, NELE_PHA
      INTEGER, DIMENSION( NELE_PHA ), intent( inout ) :: COLELE_PHA 
      INTEGER, DIMENSION( TOTELE * NPHASE + 1), intent( in ) :: FINELE_PHA

      ! Local variables...
      INTEGER :: COUNT, COUNT2, ELE_PHA, ELE_PHA2, ILOC, JLOC, IROW, JROW

      ewrite(3,*) 'In FORM_DGM_PHA_SPARSITY'

      COUNT2 = 0

      LOOP_ELEMENT_PHASE: DO ELE_PHA = 1, TOTELE * NPHASE

         LOOP_ILOC: DO ILOC = 1, U_NLOC
            IROW = ( ELE_PHA - 1 ) * U_NLOC + ILOC
            FINDGM_PHA( IROW ) = COUNT2 + 1

            DO COUNT = FINELE_PHA( ELE_PHA ) , FINELE_PHA( ELE_PHA + 1 ) - 1, 1
               ELE_PHA2 = COLELE_PHA( COUNT ) 

               LOOP_JLOC: DO JLOC = 1, U_NLOC
                  JROW = ( ELE_PHA2 - 1 ) * U_NLOC + JLOC
                  COUNT2 = COUNT2 + 1
                  COLDGM_PHA( COUNT2 ) = JROW
                  IF( IROW == JROW ) MIDDGM_PHA( IROW ) = COUNT2

               END DO LOOP_JLOC

            END DO

         END DO LOOP_ILOC

      END DO LOOP_ELEMENT_PHASE

      FINDGM_PHA( U_PHA_NONODS + 1 ) = COUNT2 + 1
      NCOLDGM_PHA = COUNT2

      IF( NCOLDGM_PHA .GT. MX_NCOLDGM_PHA ) THEN
         EWRITE(3,*) ' -------------   IN FORM_DGM_PHA_SPARSITY SUBRT.   ------------------'
         EWRITE(3,*) 'these should be NCOLDGM_PHA .GT. MX_NCOLDGM_PHA',NCOLDGM_PHA,MX_NCOLDGM_PHA  
         FLAbort(" Incorrect number of computed dimension of sparcity matrix - coupling terms for multiphase flow")
      ENDIF

      ewrite(3,*) 'Leaving FORM_DGM_PHA_SPARSITY'

      RETURN
    END SUBROUTINE FORM_DGM_PHA_SPARSITY




    SUBROUTINE DEF_SPAR_CT_DG( CV_NONODS, MX_NCT, NCT, FINDCT, COLCT, TOTELE, CV_NLOC, U_NLOC, U_NDGLN, U_ELE_TYPE)
      ! define sparsity...
      ! SEMI_BAND_WID is the semi band width.
      IMPLICIT NONE
      INTEGER, intent( in ) :: CV_NONODS, MX_NCT
      INTEGER, intent( inout ) :: NCT
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent (inout ) :: FINDCT
      INTEGER, DIMENSION( MX_NCT ), intent( inout ) :: COLCT
      INTEGER, intent( in ) :: TOTELE, CV_NLOC, U_NLOC, U_ELE_TYPE
      INTEGER, DIMENSION ( U_NLOC * TOTELE ), intent( in ) :: U_NDGLN
      ! Local variables...
      INTEGER :: CV_NOD, U_NOD, JLOC, COUNT, ELE, ELE1, ELE2, CV_NODI, CV_ILOC, ILEV

      ewrite(3,*) 'In DEF_SPAR_CT_DG'

      COUNT = 0

      IF(CV_NONODS /= CV_NLOC*TOTELE ) THEN
         ! Have a cty CV_NOD
         loop_cvnod: DO CV_NOD = 1, CV_NONODS

            FINDCT( CV_NOD ) = COUNT + 1
            ELE1 = 1 + ( CV_NOD - 2 ) / ( CV_NLOC - 1 )
            ELE2 = 1 + ( CV_NOD - 1 ) / ( CV_NLOC - 1 )

            loop_elements: DO ELE = MAX( 1 , ELE1 ), MIN( TOTELE , ELE2 ), 1

               DO JLOC = 1 ,U_NLOC
                  U_NOD = U_NDGLN(( ELE - 1 ) * U_NLOC + JLOC )
                  COUNT = COUNT + 1
                  ewrite(3,*) u_nod
                  COLCT( COUNT ) = U_NOD
               END DO

            END DO loop_elements

         END DO loop_cvnod
      ELSE
         ! Have a discontinuous CV_NOD
         loop_elements2: DO ELE = 1,TOTELE 
            loop_cviloc: DO CV_ILOC = 1, CV_NLOC

               CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC

               FINDCT( CV_NOD ) = COUNT + 1

               IF(U_ELE_TYPE==2) THEN
                  IF((CV_ILOC==1).AND.(ELE /= 1)) THEN
                     ELE2=ELE-1
                     JLOC=U_NLOC/CV_NLOC
                     DO ILEV=1,CV_NLOC
                        U_NOD = U_NDGLN(( ELE2 - 1 ) * U_NLOC + JLOC + (ILEV-1)*U_NLOC/CV_NLOC)
                        COUNT = COUNT + 1
                        COLCT( COUNT ) = U_NOD
                     END DO
                  ENDIF

                  DO JLOC = 1 ,U_NLOC
                     U_NOD = U_NDGLN(( ELE - 1 ) * U_NLOC + JLOC )
                     COUNT = COUNT + 1
                     COLCT( COUNT ) = U_NOD
                  END DO

                  IF((CV_ILOC==CV_NLOC).AND.(ELE /= TOTELE)) THEN
                     ELE2=ELE+1
                     JLOC=1
                     DO ILEV=1,CV_NLOC
                        U_NOD = U_NDGLN(( ELE2 - 1 ) * U_NLOC + JLOC + (ILEV-1)*U_NLOC/CV_NLOC)
                        COUNT = COUNT + 1
                        COLCT( COUNT ) = U_NOD
                     END DO
                  ENDIF
               ELSE
                  IF((CV_ILOC==1).AND.(ELE /= 1)) THEN
                     ELE2=ELE-1
                     JLOC=U_NLOC
                     U_NOD = U_NDGLN(( ELE2 - 1 ) * U_NLOC + JLOC )
                     COUNT = COUNT + 1
                     COLCT( COUNT ) = U_NOD
                  ENDIF

                  DO JLOC = 1 ,U_NLOC
                     U_NOD = U_NDGLN(( ELE - 1 ) * U_NLOC + JLOC )
                     COUNT = COUNT + 1
                     COLCT( COUNT ) = U_NOD
                  END DO

                  IF((CV_ILOC==CV_NLOC).AND.(ELE /= TOTELE)) THEN
                     ELE2=ELE+1
                     JLOC=1
                     U_NOD = U_NDGLN(( ELE2 - 1 ) * U_NLOC + JLOC )
                     COUNT = COUNT + 1
                     COLCT( COUNT ) = U_NOD
                  ENDIF
               ENDIF

            END DO loop_cviloc

         END DO loop_elements2
      ENDIF

      FINDCT( CV_NONODS + 1) = COUNT + 1
      NCT = COUNT

      IF(NCT > MX_NCT ) THEN
         EWRITE(3,*)'MX_NCT is not long enough NCT,MX_NCT:',NCT,MX_NCT
      ENDIF

      !DO CV_NODI=1,CV_NONODS
      !   EWRITE(3,*)'CV_NODI=',CV_NODI
      !   EWRITE(3,*)(COLCT( COUNT ),COUNT = FINDCT( CV_NODI ), FINDCT( CV_NODI + 1 ) - 1 )
      !END DO

      !ewrite(3,*)'findct', (findct(cv_nod), cv_nod=1,cv_nonods)
      !ewrite(3,*)'colct', (colct(cv_nod), cv_nod=1,nct)

      ewrite(3,*) 'Leaving DEF_SPAR_CT_DG'

      RETURN
    END SUBROUTINE DEF_SPAR_CT_DG




    SUBROUTINE CONV_CT2C( CV_NONODS, NCT, FINDCT, COLCT, U_NONODS, MX_NC, FINDC, COLC )
      IMPLICIT NONE
      ! Convert CT sparsity to C sparsity.
      INTEGER, intent( in ) :: CV_NONODS, NCT
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( NCT ), intent( in ) ::  COLCT
      INTEGER, intent( in ) :: U_NONODS, MX_NC
      INTEGER, DIMENSION( U_NONODS + 1), intent( inout ) :: FINDC
      INTEGER, DIMENSION( MX_NC ), intent( inout ) :: COLC

      ! local variables...
      INTEGER :: COL, U_NOD, CV_NOD, COUNT, COUNT2
      INTEGER, DIMENSION( : ), ALLOCATABLE ::  IN_ROW_C

      ewrite(3,*) 'In CONV_CT2C'

      ! No. of non-zero's in row of C matrix. 
      ALLOCATE( IN_ROW_C( U_NONODS ))
      IN_ROW_C = 0

      DO COUNT = 1, NCT
         COL = COLCT( COUNT )
         IN_ROW_C( COL ) = IN_ROW_C( COL ) + 1
      END DO

      COUNT = 0
      DO U_NOD = 1, U_NONODS
         FINDC( U_NOD ) = COUNT + 1
         COUNT = COUNT + IN_ROW_C( U_NOD )
      END DO
      FINDC( U_NONODS + 1 ) = COUNT + 1

      IN_ROW_C( 1 : U_NONODS )   = 0

      DO CV_NOD = 1, CV_NONODS
         DO COUNT = FINDCT( CV_NOD ), FINDCT( CV_NOD + 1 ) - 1, 1
            COL = COLCT( COUNT )
            IN_ROW_C( COL ) = IN_ROW_C( COL ) + 1
            COUNT2 = FINDC( COL ) + IN_ROW_C( COL ) - 1
            COLC( COUNT2 ) = CV_NOD

         END DO
      END DO


      DEALLOCATE( IN_ROW_C )
      EWRITE(3,*) 'HERE'

      ewrite(3,*) 'Leaving CONV_CT2C'

      RETURN

    END SUBROUTINE CONV_CT2C




    SUBROUTINE EXTEN_SPARSE_MOM_CTY( NDIM, FINMCY2, COLMCY2, NLENMCY2, NMCY2, MXNCT, &
         CV_NONODS, FINDCT, COLCT, &
         U_NONODS, NCOLC, MX_NCOLMCY, &
         FINDC, COLC, FINMCY, COLMCY, MIDMCY, NLENMCY, &
         NCOLMCY, NPHASE, NCOLCMC, FINDCMC, COLCMC  )
      IMPLICIT NONE 
      ! Extend momentum sparsity to include CTY/pressure.
      INTEGER, intent( in ) :: NDIM, NLENMCY2, NMCY2, MXNCT, CV_NONODS, U_NONODS, NCOLC, &
           MX_NCOLMCY, NLENMCY, NCOLCMC
      INTEGER, DIMENSION( NLENMCY2 + 1 ), intent( in ) :: FINMCY2
      INTEGER, DIMENSION( NMCY2 ), intent( in  ) :: COLMCY2
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( MXNCT ), intent( in ) :: COLCT
      INTEGER, DIMENSION( U_NONODS + 1 ), intent( in ) :: FINDC
      INTEGER, DIMENSION( NCOLC ), intent( in ) :: COLC
      INTEGER, DIMENSION( NLENMCY + 1), intent( inout ) :: FINMCY
      INTEGER, DIMENSION( MX_NCOLMCY ), intent( inout ) :: COLMCY
      INTEGER, DIMENSION( NLENMCY ), intent( inout ) :: MIDMCY
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( NCOLCMC ), intent( in ) :: COLCMC
      INTEGER, intent( inout ) :: NCOLMCY
      INTEGER, intent( in ) :: NPHASE
      ! Local variables...
      INTEGER :: COUNT, COUNT2, IPHASE, JPHASE, U_NOD, U_PHA_NOD, CV_NOD, ICOL, JDIM

      ewrite(3,*) 'In EXTEN_SPARSE_MOM_CTY'

      COUNT2 = 0
      COUNT  = 0
      ! put the pressure part in...
      LOOP_PHASE: DO IPHASE = 1, NPHASE
         LOOP_NODES: DO U_NOD = 1, U_NONODS

            U_PHA_NOD = U_NOD + ( IPHASE - 1 ) * U_NONODS 

            FINMCY( U_PHA_NOD ) = COUNT2 + 1

            DO COUNT = FINMCY2( U_PHA_NOD ), FINMCY2( U_PHA_NOD + 1 ) - 1, 1
               COUNT2 = COUNT2 + 1
               COLMCY( COUNT2 ) = COLMCY2( COUNT )

            END DO

            DO COUNT = FINDC( U_NOD ), FINDC( U_NOD + 1 ) - 1, 1
               COUNT2 = COUNT2 + 1
               COLMCY( COUNT2 ) = COLC( COUNT ) + U_NONODS * NPHASE

            END DO

         END DO LOOP_NODES
      END DO LOOP_PHASE

      ! put the continuity part in... 
      Loop_CVNODS: DO CV_NOD = 1, CV_NONODS  
         U_PHA_NOD = CV_NOD + U_NONODS * NPHASE
         FINMCY( U_PHA_NOD ) = COUNT2 + 1    
         Loop_Phases: DO JPHASE = 1, NPHASE 
            Loop_Dimensions: DO JDIM = 1, NDIM
               DO COUNT = FINDCT( CV_NOD ), FINDCT( CV_NOD + 1 ) - 1, 1
                  COUNT2 = COUNT2 + 1
                  ICOL = COLCT(COUNT) + (JDIM-1)*U_NONODS + (JPHASE-1)*U_NONODS*NDIM
                  COLMCY( COUNT2 ) = ICOL
               END DO
            END DO Loop_Dimensions
         END DO Loop_Phases


         ! add a diagonal which will have zero values for incompressible, but non-zero for compressible
         DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1, 1
            COUNT2 = COUNT2 + 1
            ICOL = COLCMC(COUNT) + U_NONODS*NDIM*NPHASE
            COLMCY( COUNT2 ) = ICOL
            IF(U_PHA_NOD.EQ.ICOL) MIDMCY( U_PHA_NOD ) = COUNT2

         END DO
      END DO Loop_CVNODS
      NCOLMCY = COUNT2
      FINMCY( U_NONODS *NPHASE + CV_NONODS + 1 ) = COUNT2 + 1

      ewrite(3,*) 'Leaving EXTEN_SPARSE_MOM_CTY'

      RETURN
    END SUBROUTINE EXTEN_SPARSE_MOM_CTY


    subroutine getfinele( totele, nloc, snloc, nonods, ndglno, nface_p1, &
         finele, colele, midele )
      ! This sub caluculates COLELE the element
      ! connectivitiy list in order of faces.
      implicit none
      integer, intent( in ) :: totele, nloc, snloc, nonods
      integer, dimension( totele * nloc ), intent( in ) :: ndglno
      integer, intent( in ) :: nface_p1
      integer, dimension( totele * nface_p1 ), intent( inout ) :: colele
      integer, dimension( totele + 1 ), intent( inout ) :: finele
      integer, dimension( totele ), intent( inout ) :: midele
      ! Local variables
      integer :: nface, ele, iloc, jloc, iloc2, nod, inod, jnod, count, ele2, i, hit, &
           iface, iface2, kface, ncolel, itemp, count2
      logical :: found

      integer, allocatable, dimension( : ) :: fintran, coltran, icount


      allocate( fintran( nonods + 1 ))
      allocate( coltran( max( totele, nonods ) * nface_p1 ))
      allocate( icount( max( nonods, totele )))

      nface = nface_p1 - 1

      icount = 0
      do ele = 1, totele
         do iloc = 1, nloc
            nod = ndglno( ( ele - 1 ) * nloc + iloc )
            icount( nod ) = icount( nod ) + 1
         end do
      end do

      ! Not sure if this is correct -- double check later on.
      midele( 1 ) = 1
      do ele = 2, totele
         count = 0 
         do iloc = 1, nloc
            count = count + 1
         end do
         midele( ele ) = midele( ele - 1 ) + count 
      end do

      fintran = 0
      fintran( 1 ) = 1
      do nod = 1, nonods
         fintran( nod + 1 ) = fintran( nod ) + icount( nod )
      end do

      icount = 0
      coltran = 0
      do ele = 1, totele
         do iloc = 1, nloc
            nod = ndglno( ( ele - 1 ) * nloc + iloc )
            !   ewrite(3,*)'nod, filtran, icount, ele:', nod, fintran( nod ), ele, totele, nonods
            coltran( fintran( nod ) + icount( nod )) = ele
            icount( nod ) = icount( nod ) + 1
         end do
      end do
      ! ewrite(3,*)'coltran:', coltran( 1: max(totele, nonods ) * nface_p1 )
      ! ewrite(3,*)'fintran:', fintran( 1: nonods + 1 )
      !    ewrite(3,*)'X_NDGLN:', ndglno( 1: totele*nloc )

      icount = 0
      colele = 0
      Loop_Elements1: do ele = 1, totele

         Loop_Iloc: do iloc = 1, nloc

            nod = ndglno( ( ele - 1 ) * nloc + iloc )
            Loop_Count1: do count = fintran( nod ), fintran( nod + 1 ) - 1, 1

               ele2 = coltran( count )
               found = .false. ! Add ELE2 into list FINELE and COLELE
               do i = 1, icount( ele )
                  if( colele( ( ele - 1 ) * nface_p1 + i ) == ele2 ) found = .true.
               end do

               Conditional_Found: if ( .not. found ) then ! Do elements ELE and ELE2 share at least 3 nodes?

                  hit = 0
                  do iloc2 = 1, nloc
                     inod = ndglno( ( ele - 1 ) * nloc + iloc2 )
                     do jloc = 1, nloc
                        jnod = ndglno( ( ele2 - 1 ) * nloc + jloc )
                        if ( inod == jnod ) hit = hit + 1
                     end do
                  end do
                  if ( hit >= snloc ) then
                     icount( ele ) = icount( ele ) + 1
                     colele( ( ele - 1 ) * nface_p1 + icount( ele )) = ele2 
                  end if

               end if Conditional_Found

            end do Loop_Count1

         end do Loop_Iloc

      end do Loop_Elements1

      finele( 1 ) = 1
      do ele = 1, totele
         finele( ele + 1 ) = finele( ele ) + icount( ele )
      end do

      count = 0
      Loop_Elements2: do ele = 1, totele
         ! Shorten COLELE then perform a bubble sort to get the ordering right for.
         do iface = 1, nface_p1
            if ( colele( ( ele - 1 ) * nface_p1 + iface ) /= 0 ) then
               count = count + 1
               colele( count ) = colele( ( ele - 1 ) * nface_p1 + iface )
            end if
         end do
      end do Loop_Elements2

      Loop_BubbleSort: do ele = 1, totele
         do count = finele( ele ) , finele( ele + 1 ) - 2, 1
            do count2 = finele( ele ) , finele( ele + 1 ) - 1, 1
               if ( colele( count ) > colele( count + 1 )) then ! swop over
                  itemp = colele( count + 1 )
                  colele( count + 1 ) = colele( count )
                  colele( count ) = itemp
               end if
            end do
         end do
      end do Loop_BubbleSort

      deallocate( fintran )
      deallocate( coltran )
      deallocate( icount )

      return
    end subroutine getfinele


    subroutine pousinmc2( totele, nonods1, nloc1, nonods2, nloc2, nimem, ndglno1, ndglno2, &
         lencolm, findrm, colm, centrm )
      implicit none
      integer, intent( in ) :: totele, nonods1, nloc1, nonods2, nloc2, nimem
      integer, dimension( totele * nloc1 ), intent( in ) :: ndglno1
      integer, dimension( totele * nloc2 ), intent( in ) :: ndglno2
      integer, intent( inout ) :: lencolm
      integer, dimension( nonods2 + 1 ), intent( inout ) :: findrm
      integer, dimension( nimem ), intent( inout ) :: colm
      integer, dimension( nonods2 ), intent( inout ) :: centrm

      ! Local Variables
      integer :: ele, globi, globj, loci, locj, i, irow, ptr

      ! Defining derived data types for the linked list
      type node
         integer :: id                 ! id number of node
         type( node ), pointer :: next ! next node
      end type node

      type row
         type( node ), pointer :: row ! recursive data type
      end type row

      type( row ), dimension( : ), allocatable :: matrix
      type( node ), pointer :: list, current, next

      allocate( matrix( nonods2 ))
      do i = 1, nonods2
         allocate( list )
         list % id = -1
         nullify( list % next )
         matrix( i ) % row => list
         nullify( list )
      end do

      Loop_Elements1: do ele = 1, totele

         Loop_LocI: do loci = 1, nloc2

            globi = ndglno2( ( ele - 1 ) * nloc2 + loci )
            list => matrix( globi ) % row

            Loop_LocJ: do locj = 1, nloc1

               globj = ndglno1( ( ele - 1 ) * nloc1 + locj )

               if ( list % id == -1 ) then ! Check if the list is initalised
                  list % id = globj
                  cycle
               end if

         !      if ( list % id == -1 ) then ! Check if the list is initalised
         !         list % id = globj
         !         cycle
         !      end if

               Conditional1: if ( globj < list % id ) then ! Insert at start of list
                  allocate( current )
                  current % id = globj
                  current % next => list
                  matrix( globi ) % row => current
                  list => matrix( globi ) % row

               else ! Conditional1
                  current => list
                  Loop_While1: do while ( associated( current ))

                     if ( globj == current % id ) then ! Already have this node
                        exit

                     elseif ( .not. associated( current % next )) then  ! End of list - insert this node
                        allocate( current % next )
                        nullify( current % next % next )
                        current % next % id = globj
                        exit

                     elseif ( globj < current % next % id ) then ! Insert new node here
                        allocate( next )
                        next % id = globj
                        next % next => current % next
                        current % next => next
                        exit                

                     end if
                     current => current % next

                  end do Loop_While1

               end if Conditional1

            end do Loop_LocJ

         end do Loop_LocI

      end do Loop_Elements1

      ! From matrix write COLM, FINDRM and CENTRM
      ! linked list as we go
      ptr = 1
      Loop_Irow: do irow = 1, nonods2
         findrm( irow ) = ptr
         centrm( irow ) = -1

         current => matrix( irow ) % row
         Loop_While2: do while ( associated( current ))
            assert( ptr <= nimem )
            colm( ptr ) = current % id
            if ( current % id == -1 ) then
               ewrite(0,*) "ERROR: POSINM() seriously unhappy with node", IROW
               FLAbort( "ERROR: Mesh contains nodes that are not associated with any elements." )
            end if
            if ( current % id == irow ) then
               centrm( irow ) = ptr
            endif
            next => current % next
            deallocate( current )
            current => next
            ptr = ptr + 1
         end do Loop_While2

      end do Loop_Irow

      lencolm = ptr - 1
      findrm( nonods2 + 1 ) = lencolm + 1

      deallocate( matrix )

      return
    end subroutine pousinmc2

    subroutine poscmc( totele, nonods, nimem, nct, &
         findct, colct, &
         ncmc, fincmc, colcmc, midcmc, noinod, presym )
      implicit none
      integer, intent ( in ) :: totele, nonods, nimem, nct
      integer, dimension( totele + 1 ), intent( in ) :: findct
      integer, dimension( nct ), intent( in ) :: colct
      integer, intent( inout ) :: ncmc
      integer, dimension( totele + 1 ), intent( inout ) :: fincmc
      integer, dimension( nimem ), intent( inout ) :: colcmc
      integer, dimension( totele ), intent( inout ) :: midcmc
      integer, dimension( nonods ), intent( inout ) ::  noinod
      logical, intent( inout ) :: presym

      ! Local variables
      integer :: count, count2, globi, globj, i, irow, inod, jnod, jrow, ptr 

      ! Defining derived data types for the linked list
      type node
         integer :: id                 ! id number of node
         type( node ), pointer :: next ! next node, recursive data type
      end type node

      type row
         type( node ), pointer :: row
      end type row

      type( row ), dimension( : ), allocatable :: matrix, matrix2
      type( node ), pointer :: list, current, current2, next

      allocate( matrix( nonods ))
      do i = 1, nonods
         allocate( list )
         list % id = -1
         nullify( list % next )
         matrix( i ) % row => list
         nullify( list )
      end do

      noinod = 0

      Loop_Row1: do irow = 1, totele ! Given a vel node find the pressure nodes surrounding it 
         globj = irow
         Loop_Count1: do count = findct( irow ), findct( irow + 1 ) - 1, 1
            globi = colct( count )
            list => matrix( globi ) % row
! ewrite(3,*)'list%id:',globj,globi,list % id

            if ( list % id == -1 ) then ! Check if the list is initalised
               list % id = globj
               cycle
            end if
! ewrite(3,*)'true?:', list % id == -1, globj, list % id 

            Conditional1: if ( globj < list % id ) then ! Insert at start of list
               allocate( current )
               current % id = globj
               current % next => list
               matrix( globj ) % row => current
               list => matrix( globj ) % row
            else
               current => list
 !ewrite(3,*)'-->', globj, current % id, current % next % id
               Loop_While1: do while( associated( current ))
                  if ( globj == current % id ) then ! Already have this node
                     exit
                  elseif ( .not. associated( current % next )) then ! End of list - insert this node
                     allocate( current % next )
                     nullify( current % next % next )
                     current % next % id = globj
                     exit
                  elseif( globj < current % next % id ) then ! Insert new node here
                     allocate( next )
                     next % id = globj
                     next % next => current % next
                     current % next => next
                     exit
                  end if
                  current => current % next
               end do Loop_While1
            end if Conditional1

            noinod( globi ) = noinod( globi ) + 1

         end do Loop_Count1

      end do Loop_Row1

      allocate( matrix2( totele )) ! Initalise the linked lists
      do i = 1, totele
         allocate( list )
         list % id = -1
         nullify( list % next )
         matrix2( i ) % row => list
         nullify( list )
      end do

      Loop_Row2: do irow = 1, totele
         ! Find the pressure nodes surrounding pressure node IROW
         Loop_Count2: do count = findct( irow ), findct( irow + 1 ) - 1, 1 
            inod = colct( count )

            ! Find the pressure nodes surrounding node INOD 
            ! these will be connected to pressure node IROW.
            current => matrix( inod ) % row
            Loop_While2: do while( associated( current ))
               jrow = current % id

               Conditional2: if (( .not. presym ) .or. ( jrow >= irow )) then
                  list => matrix2( irow ) % row

                  if ( list % id == -1 ) then
                     list % id = jrow
                     cycle
                  end if

                  Conditional3: if ( jrow < list % id ) then ! Insert at start of list
                     allocate( current2 )
                     current2 % id = jrow
                     current2 % next => list
                     matrix2( irow ) % row => current2
                     list => matrix2( irow ) % row

                  else
                     current2 => list
                     Loop_While3: do while( associated( current2 ))

                        Conditional4: if ( jrow == current2 % id ) then ! Already have this node
                           exit
                        elseif( .not. associated( current2 % next )) then ! End of list - insert this node
                           allocate( current2 % next )
                           nullify(  current2 % next % next )
                           current2 % next % id = jrow
                           exit 
                        elseif( jrow < current2 % next % id ) then ! Insert new node here
                           allocate( next )
                           next % id = jrow
                           next % next => current2 % next
                           current2 % next => next
                           exit
                        end if Conditional4
                        current2 => current2 % next

                     end do Loop_While3

                  end if Conditional3

               end if Conditional2
               current => current % next          

            end do Loop_While2

         end do Loop_Count2

      end do Loop_Row2


      do irow = 1, nonods ! Delete Matrix
         current => matrix( irow ) % row
         do while ( associated( current ))
            next => current % next
            deallocate( current )
            current => next
         end do
      end do
      deallocate( matrix )

      ptr = 1
      Loop_Row3: do irow = 1, totele ! From matrix write COLCMC, FINCMC and MIDCMC
         midcmc( irow ) = -1
         fincmc( irow ) = ptr
         current => matrix2( irow ) % row ! Warning: overwriten 

         Loop_While4: do while ( associated( current ))

            if ( ptr > nimem ) then
               ewrite( -1, * ) 'nimem, ptr: ',nimem, ptr
               ewrite( -1, * ) 'totele, irow: ',totele, irow
               FLAbort( "Integer memory too small" )
            end if

            colcmc( ptr ) = current % id
            if( current % id == irow ) then
               midcmc( irow ) = ptr
            end if

            next => current % next
            deallocate( current )
            current => next
            ptr = ptr + 1
!ewrite(3,*)'ptr:', ptr

         end do Loop_While4
!ewrite(3,*)'irow:', irow
      end do Loop_Row3

      ncmc =  ptr - 1
      fincmc( totele + 1 ) = ncmc + 1

      return
    end subroutine poscmc

    subroutine comparematrices( n, vec1, vec2 )
      implicit none
      integer :: n
      integer, dimension( n ) :: vec1, vec2
      integer :: i 
      real :: res1, res2, diff

      res1 = 0.
      res2 = 0.
      diff = 0.
      do i = 1, n
         res1 = res1 + real( vec1( i ))**2 
         res2 = res2 + real( vec2( i ))**2 
         diff = diff + real( abs(vec1( i )) - abs(vec2( i )))**2
         !    ewrite(3,*) i, vec1( i ), vec2( i ), diff
      end do
      res1 = sqrt( res1 )
      res2 = sqrt( res2 )
      diff = sqrt( diff )

      ewrite(3,*) 'vec1, vec2, diff: ', res1, res2, diff

      return
    end subroutine comparematrices


    SUBROUTINE POSINMAT( POSMAT, GLOBI, GLOBJ,& 
         NONODS, FINDRM, COLM, NCOLM)
      ! Find position in matrix POSMAT which has column GLOBJ
      INTEGER, intent( in ) :: NONODS, NCOLM, GLOBI, GLOBJ
      INTEGER, intent( inout ) :: POSMAT
      INTEGER, DIMENSION( NONODS + 1), intent( in ) :: FINDRM
      INTEGER, DIMENSION( NCOLM ), intent( in ) :: COLM
      ! Local
      INTEGER :: INUM, LOWER, UPPER, COUNT
      INTEGER, PARAMETER :: NMAX = 100000

      !ewrite(3,*) 'In POSINMAT'

      LOWER = FINDRM( GLOBI ) 
      UPPER = FINDRM( GLOBI + 1 ) - 1


      COUNT = 1
      Loop_While: DO WHILE  ( COUNT <= NMAX )

         INUM = LOWER + ( UPPER - LOWER + 1 ) / 2

         IF(GLOBJ >= COLM( INUM ) ) THEN 
            LOWER = INUM
         ELSE
            UPPER = INUM
         ENDIF

         IF( ( UPPER - LOWER ) <= 1 ) THEN
            IF( GLOBJ == COLM( LOWER )) THEN
               POSMAT = LOWER
            ELSE
               POSMAT = UPPER
            ENDIF

            RETURN

         ENDIF

         COUNT = COUNT + 1

      END DO Loop_While

      !ewrite(3,*) 'Leaving POSINMAT'

      RETURN 
    END SUBROUTINE POSINMAT


    subroutine checksparsity( option_mid, unit, ncol2, nonods, ncol, find, mid, col )
      implicit none
      logical, intent( in ) :: option_mid
      integer, intent( in ) :: unit, ncol2, nonods, ncol
      integer, dimension( nonods + 1 ), intent( in ) :: find
      integer, dimension( nonods ), intent( in ) :: mid
      integer, dimension( ncol ), intent( in ) :: col
      ! Local variables
      integer :: inod, icol

      ewrite(3,*) 'In checksparsity'

      write( unit, * )'find:', ( find( inod ), inod = 1, nonods + 1 )
      if( option_mid )write( unit, * )'mid:', ( mid( inod ), inod = 1, nonods )
      write( unit, * )'col:', ( col( icol ), icol = 1, ncol2  )

      ewrite(3,*) 'Leaving checksparsity'

      return
    end subroutine checksparsity

    subroutine swapprint( first, n, nc, nc2, fin, col, mid )
      implicit none
      logical, intent( in ) :: first
      integer, intent( in ) :: n, nc
      integer, intent( inout ) :: nc2
      integer, dimension( n + 1 ), intent( inout ) :: fin
      integer, dimension( nc ), intent( inout ) :: col
      integer, dimension( n ), intent( inout ) :: mid

      ewrite(3,*)'fin:', size( fin ), fin( 1 : n+1 )
      ewrite(3,*)'col:', nc2, col( 1 : nc )
      ewrite(3,*)'mid:', size( mid ), mid( 1 : n )

      if( first ) then
         fin = 0
         col = 0
         mid = 0
      end if

      return
    end subroutine swapprint

  end module spact


  module spact2
    use fldebug
    use shape_functions

  type cv_faces_type
    ! loc = number of vertices, faces = number of faces
    ! degree = degree of polynomial, dim = dimensions of parent element
    integer :: loc, faces, coords, degree, dim
    integer :: sloc, sfaces, scoords
    ! corners = volume coordinates of corners of faces
    ! faces x loc x face vertices
    real, dimension(:,:,:), pointer :: corners, scorners
    ! neiloc = relates faces to faces and vice versa
    ! loc x faces
    integer, dimension(:,:), pointer :: neiloc, sneiloc
    ! shape = shape function used in quadrature of faces
    ! 1 dimension lower than parent element
  !  type(element_type) :: shape
  end type cv_faces_type

  contains

    subroutine sparse_cvdomain( ndim, totele, cv_nloc, u_nloc,  &
         cv_ele_type, cv_ndgln, u_ndgln, &
         ncolgpts, findgpts, colgpts )
      implicit none
      integer, intent( in ) :: ndim, totele, cv_nloc, u_nloc, cv_ele_type
      integer, dimension( totele * cv_nloc ), intent( in ) :: cv_ndgln
      integer, dimension( totele * u_nloc ), intent( in ) :: u_ndgln
      integer, intent( inout ) :: ncolgpts
      integer, dimension( cv_nloc + 1 ), intent( inout ) :: findgpts
      integer, dimension( cv_nloc * totele ), intent( inout ) :: colgpts

      ! Local variables: 
      integer, dimension( : , : ), allocatable :: cv_neiloc
      ! integer, dimension( : ), allocatable :: 
      integer :: cv_ele_type2, u_nloc2, cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, &
           ele, cv_iloc, cv_jloc, cv_nodi, gcount, gi

      Conditional_CVELETYPE: if( cv_ele_type == 2 ) then
         cv_ele_type2 = 1
         u_nloc2 = u_nloc / cv_nloc

         call retrieve_ngi2( ndim, cv_ele_type2, cv_nloc, u_nloc2, &
              cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface ) ! obtain cv_ngi and scvngi

         allocate( cv_neiloc( cv_nloc, scvngi ))
         call volnei( cv_neiloc, cv_nloc, scvngi, cv_ele_type2 ) ! obtain cv_neiloc

         call gaussiloc( findgpts, colgpts( 1 : cv_nloc * scvngi ), ncolgpts, &
              cv_neiloc, cv_nloc, scvngi ) 
         !  findgpts and colgpts has the address of Gauss point around iloc

         Loop_Elements1: do ele = 1, totele

            Loop_cviloc: do cv_iloc = 1, cv_nloc
               cv_nodi = cv_ndgln( ( ele - 1 ) * cv_nloc + cv_iloc )

               ! Loop over quadrature (gauss) points in ele neighbouring iloc
               Loop_GaussCount: do gcount = findgpts( cv_iloc ), findgpts( cv_iloc + 1 ), 1 
                  gi = colgpts( gcount ) ! colgpts stores the local Gauss-point number in ele
                  cv_jloc = cv_neiloc( cv_iloc, gi ) ! get the neighbouring node for node iloc and Gauss point gi

               end do Loop_GaussCount

            end do Loop_cviloc

         end do Loop_Elements1

      end if Conditional_CVELETYPE

      return
    end subroutine sparse_cvdomain



    subroutine retrieve_ngi2( ndim, cv_ele_type, cv_nloc, u_nloc, &
         cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface )
      implicit none
      integer, intent( in ) :: ndim, cv_ele_type, cv_nloc, u_nloc
      integer, intent( inout ) :: cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface
      !! It is necessary to extend it 2/3-D

      Conditional_NDIM: Select Case( ndim )
      case( 1 )

         Conditional_CV_NLOC: Select Case( cv_nloc )

         case( 1 )
            cv_ngi = 1
            scvngi = 2
         case( 2 )
            cv_ngi = 12
            scvngi = 3
         case( 3 )
            cv_ngi = 12
            scvngi = 4
            if( ( u_nloc == 4 ) .or. ( u_nloc == 5 )) cv_ngi = 18

         case default; FLExit("Invalid integer for u_nloc ")

         end Select Conditional_CV_NLOC
         sbcvngi = 1
         nface = 2
         cv_ngi_short = cv_ngi
         if( cv_ele_type == 2 ) cv_ngi = cv_ngi * cv_nloc

      case default ; FLExit("Invalid integer for problem dimension")
      end Select Conditional_NDIM

      return
    end subroutine retrieve_ngi2



  end module spact2


