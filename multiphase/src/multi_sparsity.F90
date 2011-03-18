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

module spact

contains


  SUBROUTINE GET_SPARS_PATS(   &
       NDIM, U_PHA_NONODS, CV_PHA_NONODS,   &
       U_NONODS, CV_NONODS, &
       U_NLOC, CV_NLOC, NPHASE, TOTELE, U_NDGLN,  &
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

    INTEGER, intent( in ) :: NDIM, U_PHA_NONODS, CV_PHA_NONODS, U_NONODS, CV_NONODS, U_NLOC, CV_NLOC, &
         NPHASE, TOTELE, U_ELE_TYPE
    INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in )     :: U_NDGLN

    INTEGER, intent( in )                                   :: MX_NCOLACV 
    INTEGER, intent( inout )                                :: NCOLACV

    INTEGER, DIMENSION( CV_PHA_NONODS + 1), intent( inout ) :: FINACV ! CV multi-phase eqns (e.g. vol frac, temp)
    INTEGER, DIMENSION( MX_NCOLACV ), intent( inout )          :: COLACV
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
         
    write(357,*) 'In GET_SPARS_PATS'

    ! Extend the element sparsity to an element multiphase sparsity

    MX_NCOLELE_PHA = NPHASE * MXNELE + ( NPHASE - 1 ) * NPHASE * TOTELE
    NPHA_TOTELE = TOTELE * NPHASE

    ALLOCATE( COLELE_PHA( MX_NCOLELE_PHA ))
    ALLOCATE( FINELE_PHA( NPHA_TOTELE + 1 ))
    ALLOCATE( MIDELE_PHA( NPHA_TOTELE ))
    ALLOCATE( CENTCT( CV_NONODS ))

    CALL DEF_SPAR( 1, TOTELE, MXNELE, NCOLELE, &
         MIDELE, FINELE, COLELE )


    CALL EXTEN_SPARSE_MULTI_PHASE( TOTELE, MXNELE, FINELE, COLELE, &
         NPHASE, NPHA_TOTELE, MX_NCOLELE_PHA, FINELE_PHA, COLELE_PHA, MIDELE_PHA ) 


    CALL FORM_DGM_PHA_SPARSITY( U_PHA_NONODS, MX_NCOLDGM_PHA, NCOLDGM_PHA, COLDGM_PHA, FINDGM_PHA, MIDDGM_PHA, &
         TOTELE, NPHASE, U_NLOC, MX_NCOLELE_PHA, FINELE_PHA, COLELE_PHA ) 


    ! Now form the global matrix: FINMCY, COLMCY and MIDMCY for the 
    ! momentum and continuity eqns: 

    write(357,*) 'u_ndgl:',(U_NDGLN(II), II=1,TOTELE * U_NLOC)
    CALL DEF_SPAR_CT_DG( CV_NONODS, MX_NCT, NCT, FINDCT, COLCT, TOTELE, CV_NLOC, U_NLOC, U_NDGLN, U_ELE_TYPE)
    NC = NCT
    write(357,*) 'nc, MX_NC, nct, MX_NCT:',nc, MX_NC, nct, MX_NCT

    ! Convert CT sparsity to C sparsity.
    CALL CONV_CT2C( CV_NONODS, NCT, FINDCT, COLCT, U_NONODS, MX_NC, FINDC, COLC )

      write(357,*) 'NCOLDGM_PHA,CV_NONODS, TOTELE*CV_NLOC:',NCOLDGM_PHA,CV_NONODS, TOTELE*CV_NLOC
      
    IF(CV_NONODS /= TOTELE*CV_NLOC) THEN
    CALL DEF_SPAR( CV_NLOC - 1 , CV_NONODS, MX_NCOLCMC, NCOLCMC, &
!    CALL DEF_SPAR( CV_NLOC+2, CV_NONODS, MX_NCOLCMC, NCOLCMC, &
         MIDCMC, FINDCMC, COLCMC )
    ELSE
     write(357,*)'CV_NLOC,MX_NCOLCMC=',CV_NLOC,MX_NCOLCMC
!    CALL DEF_SPAR( CV_NLOC - 1 , CV_NONODS, MX_NCOLCMC, NCOLCMC, &
    CALL DEF_SPAR( CV_NLOC+2, CV_NONODS, MX_NCOLCMC, NCOLCMC, &
         MIDCMC, FINDCMC, COLCMC )
    ENDIF
    if(MX_NCOLCMC < NCOLCMC) stop 272
    write(357,*)'defining sparsity of matrix'

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

    CALL DEF_SPAR( 1 , CV_NONODS, MXNACV_LOC, NACV_LOC, &
         MIDACV_LOC, FINACV_LOC, COLACV_LOC )

    NCOLACV = NPHASE * NACV_LOC + ( NPHASE - 1 ) * NPHASE * CV_NONODS
    NACV_LOC = NCOLACV


    write(357,*) 'going into EXTEN_SPARSE_MULTI_PHASE sbrt 2nd time'
    CALL EXTEN_SPARSE_MULTI_PHASE( CV_NONODS, MXNACV_LOC, FINACV_LOC, COLACV_LOC, &
         NPHASE, CV_PHA_NONODS, MX_NCOLACV, FINACV, COLACV, MIDACV )  

    ! Sparsity of CV-FEM
     write(357,*)'going to find sparsity of colm MX_NCOLM:',MX_NCOLM
     write(357,*)'CV_NONODS,CV_NLOC:',CV_NONODS,CV_NLOC
    CALL DEF_SPAR( CV_NLOC - 1 , CV_NONODS, MX_NCOLM, NCOLM, &
         MIDM, FINDM, COLM )

    write(357,*) 'NCOLMCY, NCOLACV, NCOLCMC:', NCOLMCY, NCOLACV, NCOLCMC

    write(357,*) 'Leaving GET_SPARS_PATS'

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
    
    write(357,*) 'In DEF_SPAR'

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
    
    write(357,*) 'Leaving DEF_SPAR'

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

    write(357,*) 'In EXTEN_SPARSE_MULTI_PHASE'

    COUNT2 = 0
  
    write(357,*) ' ########  In  EXTEN_SPARSE_MULTI_PHASE Subrt.  ####### '

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
    
    write(357,*) 'Leaving EXTEN_SPARSE_MULTI_PHASE'

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
    
    write(357,*) 'In FORM_DGM_PHA_SPARSITY'

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
       WRITE(357,*) '***************** SOMETHING WRONG -- REVIEW IT !! *******************'
       WRITE(357,*) ' -------------   IN FORM_DGM_PHA_SPARSITY SUBRT.   ------------------'
       WRITE(357,*) 'these should be NCOLDGM_PHA .GT. MX_NCOLDGM_PHA',NCOLDGM_PHA,MX_NCOLDGM_PHA  !!double check it 
       STOP 383
    ENDIF
    
    write(357,*) 'Leaving FORM_DGM_PHA_SPARSITY'

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
    
    write(357,*) 'In DEF_SPAR_CT_DG'

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
            write(357,*) u_nod
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
       WRITE(357,*)'MX_NCT is not long enough NCT,MX_NCT:',NCT,MX_NCT
    ENDIF

    DO CV_NODI=1,CV_NONODS
      WRITE(357,*)'CV_NODI=',CV_NODI
      WRITE(357,*)(COLCT( COUNT ),COUNT = FINDCT( CV_NODI ), FINDCT( CV_NODI + 1 ) - 1 )
    END DO

    write(357,*)'findct', (findct(cv_nod), cv_nod=1,cv_nonods)
    write(357,*)'colct', (colct(cv_nod), cv_nod=1,nct)
    
    write(357,*) 'Leaving DEF_SPAR_CT_DG'

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
    
    write(357,*) 'In CONV_CT2C'

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
    WRITE(357,*) 'HERE'
    
    write(357,*) 'Leaving CONV_CT2C'

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
    
    write(357,*) 'In EXTEN_SPARSE_MOM_CTY'

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
    
    write(357,*) 'Leaving EXTEN_SPARSE_MOM_CTY'

    RETURN
  END SUBROUTINE EXTEN_SPARSE_MOM_CTY



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
    
    !write(357,*) 'In POSINMAT'

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
    
    !write(357,*) 'Leaving POSINMAT'

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
    
    write(357,*) 'In checksparsity'

    write( unit, * )'find:', ( find( inod ), inod = 1, nonods + 1 )
    if( option_mid )write( unit, * )'mid:', ( mid( inod ), inod = 1, nonods )
    write( unit, * )'col:', ( col( icol ), icol = 1, ncol2  )
    
    write(357,*) 'Leaving checksparsity'

    return
  end subroutine checksparsity


end module spact


