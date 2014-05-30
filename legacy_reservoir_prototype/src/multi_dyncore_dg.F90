
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

module multiphase_1D_engine

    use state_module 
    use fields
    use field_options
    use spud
    use global_parameters, only: option_path_len, is_overlapping, is_compact_overlapping
    use futils, only: int2str

    use Fields_Allocates, only : allocate

    use solvers_module
    use mapping_for_ocvfem
    use cv_advection  
    use matrix_operations
    use shape_functions
    use spact
    use Copy_Outof_State
    use multiphase_EOS
    use Copy_Outof_State, only: as_vector
    use fldebug
    use solvers

    implicit none

    private :: UVW_2_ULONG, &
    CV_ASSEMB_FORCE_CTY, &
    PUT_MOM_C_IN_GLOB_MAT, &
    PUT_CT_IN_GLOB_MAT, &
    ASSEMB_FORCE_CTY, &
    DG_DIFFUSION, &
    ASSEM_CS, &
    AVESOU, &
    AVESIG, &
    LUMP_ENERGY_EQNS

    public  :: INTENERGE_ASSEM_SOLVE, &
    VolumeFraction_Assemble_Solve, &
    FORCE_BAL_CTY_ASSEM_SOLVE

contains

    SUBROUTINE INTENERGE_ASSEM_SOLVE( state, packed_state, &
         tracer, velocity, density, &
    NCOLACV, FINACV, COLACV, MIDACV, &
    SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
    block_to_global_acv, global_dense_block_acv, &
    NCOLCT, FINDCT, COLCT, &
    CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
    U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE, &
    NPHASE,  &
    CV_NLOC, U_NLOC, X_NLOC,  &
    CV_NDGLN, X_NDGLN, U_NDGLN, &
    CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
    T, TOLD, &
    MAT_NLOC,MAT_NDGLN,MAT_NONODS, TDIFFUSION, IGOT_THERM_VIS, THERM_U_DIFFUSION, &
    T_DISOPT, T_DG_VEL_INT_OPT, DT, T_THETA, T_BETA, &
    SUF_SIG_DIAGTEN_BC, &
    DERIV, &
    T_SOURCE, T_ABSORB, VOLFRA_PORE, &
    NDIM, &
    NCOLM, FINDM, COLM, MIDM, &
    XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
    OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
    T_FEMT, DEN_FEMT, &
    IGOT_T2, T2, T2OLD, igot_theta_flux,SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
    THETA_GDIFF, &
    IN_ELE_UPWIND, DG_ELE_UPWIND, &
    NOIT_DIM, &
    MEAN_PORE_CV, &
    option_path, &
    mass_ele_transp, &
    thermal, THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, &
    StorageIndexes, icomp, saturation )

        ! Solve for internal energy using a control volume method.

        implicit none
        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type(tensor_field), intent(inout) :: tracer
        type(tensor_field), intent(in) :: velocity, density

        INTEGER, intent( in ) :: NCOLACV, NCOLCT, CV_NONODS, U_NONODS, X_NONODS, MAT_NONODS, TOTELE, &
        U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE, NPHASE, CV_NLOC, U_NLOC, X_NLOC,  MAT_NLOC, &
        CV_SNLOC, U_SNLOC, STOTEL, XU_NLOC, NDIM, NCOLM, NCOLELE, &
        NOPT_VEL_UPWIND_COEFS, &
        IGOT_T2, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND, igot_theta_flux
        LOGICAL, intent( in ) :: GET_THETA_FLUX, USE_THETA_FLUX
        LOGICAL, intent( in ), optional ::THERMAL
        INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: U_SNDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: FINACV
        INTEGER, DIMENSION( : ), intent( in ) :: COLACV
        INTEGER, DIMENSION( : ), intent( in ) :: MIDACV
        INTEGER, DIMENSION( : ), intent( in ) :: SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV
        integer, dimension(:)    :: block_to_global_acv
        integer, dimension (:,:) :: global_dense_block_acv
        INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
        INTEGER, DIMENSION( : ), intent( in ) :: COLCT
        REAL, DIMENSION( : ), intent( inout ) :: T, T_FEMT, DEN_FEMT
        REAL, DIMENSION( : ), intent( in ) :: TOLD
        REAL, DIMENSION( : ), intent( in ) :: T2, T2OLD
        REAL, DIMENSION( :, : ), intent( inout ) :: THETA_GDIFF
        REAL, DIMENSION( :,: ), intent( inout ), optional :: THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J
        REAL, DIMENSION( :,:,:, : ), intent( in ) :: TDIFFUSION
        INTEGER, intent( in ) :: IGOT_THERM_VIS
        REAL, DIMENSION(NDIM,NDIM,NPHASE,MAT_NONODS*IGOT_THERM_VIS), intent( in ) :: THERM_U_DIFFUSION
        INTEGER, intent( in ) :: T_DISOPT, T_DG_VEL_INT_OPT
        REAL, intent( in ) :: DT, T_THETA
        REAL, intent( in ) :: T_BETA
        REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
        REAL, DIMENSION( NPHASE, CV_NONODS ), intent( in ) :: DERIV
        REAL, DIMENSION( : ), intent( in ) :: T_SOURCE
        REAL, DIMENSION( : , : , : ), intent( in ) :: T_ABSORB
        REAL, DIMENSION( : ), intent( in ) :: VOLFRA_PORE
        INTEGER, DIMENSION( : ), intent( in ) :: FINDM
        INTEGER, DIMENSION( : ), intent( in ) :: COLM
        INTEGER, DIMENSION( : ), intent( in ) :: MIDM
        INTEGER, DIMENSION( : ), intent( in ) :: FINELE
        INTEGER, DIMENSION( : ), intent( in ) :: COLELE
        REAL, DIMENSION( : ), intent( in ) :: OPT_VEL_UPWIND_COEFS
        INTEGER, INTENT(IN) :: NOIT_DIM
        REAL, DIMENSION( : ), intent( inout ) :: MEAN_PORE_CV
        character( len = * ), intent( in ), optional :: option_path
        real, dimension( : ), intent( inout ), optional :: mass_ele_transp
        integer, dimension(:), intent(inout) :: StorageIndexes
        type(tensor_field), intent(in), optional :: saturation
        
        integer, optional :: icomp
        ! Local variables
        LOGICAL, PARAMETER :: GETCV_DISC = .TRUE., GETCT= .FALSE., RETRIEVE_SOLID_CTY=.FALSE.
        integer :: nits_flux_lim, its_flux_lim
        logical :: lump_eqns
        REAL, DIMENSION( : ), allocatable :: ACV, CV_RHS, DIAG_SCALE_PRES, CT_RHS
        REAL, DIMENSION( : ), allocatable :: block_acv, mass_mn_pres
        REAL, DIMENSION( : , : , : ), allocatable :: dense_block_matrix, CT
        REAL, DIMENSION( : , : ), allocatable :: den_all, denold_all
        REAL, DIMENSION( : ), allocatable :: CV_RHS_SUB, ACV_SUB
        type( scalar_field ), pointer :: P
        INTEGER, DIMENSION( : ), allocatable :: COLACV_SUB, FINACV_SUB, MIDACV_SUB
        INTEGER :: NCOLACV_SUB, IPHASE, I, J
        REAL :: SECOND_THETA
        INTEGER :: STAT
        character( len = option_path_len ) :: path
        type(vector_field) :: rhs_field
        type( tensor_field ), pointer :: den_all2, denold_all2
        integer :: lcomp

        type(petsc_csr_matrix) :: petsc_acv
        type(vector_field)  :: vtracer

        if (present(icomp)) then
           lcomp=icomp
        else
           lcomp=0
        end if

        call allocate(rhs_field,nphase,tracer%mesh,"RHS")

        ALLOCATE( ACV( NCOLACV ) )
        ALLOCATE( mass_mn_pres( size(small_COLACV ) ))
        allocate( block_acv(size(block_to_global_acv) ) )
        allocate( dense_block_matrix (nphase,nphase,cv_nonods) ); dense_block_matrix=0;
        ALLOCATE( CV_RHS( CV_NONODS * NPHASE ) )


        allocate( den_all( nphase, cv_nonods ), denold_all( nphase, cv_nonods ) )


        if ( thermal ) then
           p => extract_scalar_field( packed_state, "CVPressure" )
           den_all2 => extract_tensor_field( packed_state, "PackedDensityHeatCapacity" )
           denold_all2 => extract_tensor_field( packed_state, "PackedOldDensityHeatCapacity" )
           den_all    = den_all2    % val ( 1, :, : )
           denold_all = denold_all2 % val ( 1, :, : )
        else if ( lcomp > 0 ) then
           p => extract_scalar_field( packed_state, "FEPressure" )
           den_all2 => extract_tensor_field( packed_state, "PackedComponentDensity" )
           denold_all2 => extract_tensor_field( packed_state, "PackedOldComponentDensity" )
           den_all = den_all2 % val ( icomp, :, : )
           denold_all = denold_all2 % val ( icomp, :, : )
        else
           p => extract_scalar_field( packed_state, "Pressure" )
              den_all=1.0
              denold_all=1.0
        end if


        if( present( option_path ) ) then

            if( trim( option_path ) == '/material_phase[0]/scalar_field::Temperature' ) then
                call get_option( '/material_phase[0]/scalar_field::Temperature/prognostic/temporal_discretisation/' // &
                'control_volumes/number_advection_iterations', nits_flux_lim, default = 3 )
            end if

            path='/material_phase[0]/scalar_field::Temperature/prognostic/temporal_discretisation' // &
            '/control_volumes/second_theta'
            call get_option( path, second_theta, default=1. )

        else

            call get_option( '/material_phase[' // int2str( nphase ) // ']/scalar_field::ComponentMassFractionPhase1/' // &
            'prognostic/temporal_discretisation/control_volumes/number_advection_iterations', nits_flux_lim, default = 1 )

            path= '/material_phase[' // int2str( nphase ) // ']/scalar_field::ComponentMassFractionPhase1/' // &
            'prognostic/temporal_discretisation/control_volumes/second_theta'

            call get_option( path, second_theta, default=1. )

        end if

        lump_eqns = have_option( '/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
        'spatial_discretisation/continuous_galerkin/mass_terms/lump_mass_matrix' )


        Loop_NonLinearFlux: DO ITS_FLUX_LIM = 1, NITS_FLUX_LIM


            call CV_ASSEMB( state, packed_state, &
                 tracer, velocity, density, &
            CV_RHS, &
            NCOLACV, block_acv, DENSE_BLOCK_MATRIX, FINACV, COLACV, MIDACV, &
            SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV,&
            NCOLCT, CT, DIAG_SCALE_PRES, CT_RHS, FINDCT, COLCT, &
            CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
            CV_ELE_TYPE, &
            NPHASE, &
            CV_NLOC, U_NLOC, X_NLOC, &
            CV_NDGLN, X_NDGLN, U_NDGLN, &
            CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
            DEN_ALL, DENOLD_ALL, &
            MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, IGOT_THERM_VIS, THERM_U_DIFFUSION, &
            T_DISOPT, T_DG_VEL_INT_OPT, DT, T_THETA, SECOND_THETA, T_BETA, &
            SUF_SIG_DIAGTEN_BC, &
            DERIV, P%val, &
            T_SOURCE, T_ABSORB, VOLFRA_PORE, &
            NDIM, GETCV_DISC, GETCT, &
            NCOLM, FINDM, COLM, MIDM, &
            XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
            OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
            DEN_FEMT, &
            IGOT_T2, T2, T2OLD,IGOT_THETA_FLUX ,SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
            THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, THETA_GDIFF, &
            IN_ELE_UPWIND, DG_ELE_UPWIND, &
            NOIT_DIM, &
            MEAN_PORE_CV, &
            SMALL_FINACV, SMALL_COLACV, size(small_colacv), mass_Mn_pres, THERMAL, RETRIEVE_SOLID_CTY, &
            mass_ele_transp,&
            StorageIndexes, -1, T_input = T, TOLD_input=TOLD, FEMT_input =  T_FEMT ,&
            saturation=saturation)
            t=0.

            Conditional_Lumping: IF ( LUMP_EQNS ) THEN
                ! Lump the multi-phase flow eqns together
                ALLOCATE( CV_RHS_SUB( CV_NONODS ) )

                CV_RHS_SUB = 0.0
                DO IPHASE = 1, NPHASE
                    CV_RHS_SUB( : ) = CV_RHS_SUB( : ) + CV_RHS( 1 +( IPHASE - 1) * CV_NONODS : &
                    IPHASE * CV_NONODS )
                END DO

                NCOLACV_SUB = FINACV( CV_NONODS + 1) - 1 - CV_NONODS *( NPHASE - 1 )

                ALLOCATE( ACV_SUB( NCOLACV_SUB ))
                ALLOCATE( COLACV_SUB( NCOLACV_SUB ))
                ALLOCATE( FINACV_SUB( CV_NONODS + 1 ))
                ALLOCATE( MIDACV_SUB( CV_NONODS ))

                CALL LUMP_ENERGY_EQNS( CV_NONODS, NPHASE, &
                NCOLACV, NCOLACV_SUB, &
                FINACV, COLACV, COLACV_SUB, FINACV_SUB, ACV_SUB )
                CALL SOLVER( ACV_SUB, T, CV_RHS_SUB, &
                FINACV_SUB, COLACV_SUB, &
                trim(option_path))

                DO IPHASE = 2, NPHASE
                   T( 1 + ( IPHASE - 1 ) * CV_NONODS : IPHASE * CV_NONODS ) = T ( 1 : CV_NONODS )
                END DO

             ELSE

                call assemble_global_multiphase_csr(acv,&
                     block_acv,dense_block_matrix,&
                     block_to_global_acv,global_dense_block_acv)

                call assemble_global_multiphase_petsc_csr(petsc_acv,&
                     block_acv,dense_block_matrix,&
                     finacv,colacv,p%mesh%halos)
            
                T([([(i+(j-1)*nphase,j=1,cv_nonods)],i=1,nphase)]) = T

                IF ( IGOT_T2 == 1) THEN
                   vtracer=as_vector(tracer,dim=2)
                   
                   call zero (vtracer)
                   rhs_field%val(:,:)=reshape(cv_rhs,[nphase,cv_nonods])
                   call petsc_solve(vtracer,petsc_acv,rhs_field,'/material_phase::Component1/scalar_field::ComponentMassFractionPhase1/prognostic')
!                   CALL SOLVER( ACV, T, CV_RHS, &
!                        FINACV, COLACV, &
!                        trim('/material_phase::Component1/scalar_field::ComponentMassFractionPhase1/prognostic') )

                   T([([(i+(j-1)*cv_nonods,j=1,nphase)],i=1,cv_nonods)]) = [vtracer%val]

                ELSE
                   vtracer=as_vector(tracer,dim=2)
                   call zero (vtracer)
                   rhs_field%val(:,:)=reshape(cv_rhs,[nphase,cv_nonods])
                   call petsc_solve(vtracer,petsc_acv,rhs_field,trim(option_path))
!                   CALL SOLVER( ACV, T, CV_RHS, &
!                        FINACV, COLACV, &
!                        trim(option_path) )

                   T([([(i+(j-1)*cv_nonods,j=1,nphase)],i=1,cv_nonods)]) = [tracer%val]

                END IF

!!                call dump_petsc_csr_matrix(petsc_acv)


                !ewrite(3,*)'cv_rhs:', cv_rhs
                !ewrite(3,*)'SUF_T_BC:',SUF_T_BC
                !ewrite(3,*)'ACV:',  (acv(i),i= FINACV(1), FINACV(2)-1)
                !ewrite(3,*)'T_ABSORB:',((T_ABSORB(1,i,j), i=1,nphase),j=1,nphase)
                !ewrite(3,*)

             END IF Conditional_Lumping

        END DO Loop_NonLinearFlux

        DEALLOCATE( ACV )
        deALLOCATE( mass_mn_pres )
        deallocate( block_acv, dense_block_matrix )
        DEALLOCATE( CV_RHS )
        call deallocate(RHS_FIELD)

        ewrite(3,*)'t:', t
        !ewrite(3,*)'told:', told

        ewrite(3,*) 'Leaving INTENERGE_ASSEM_SOLVE'

    END SUBROUTINE INTENERGE_ASSEM_SOLVE








    SUBROUTINE CV_ASSEMB_CV_DG( state, packed_state, &
         tracer, velocity, density, pressure, &
    CV_RHS, &
    NCOLACV, ACV, DENSE_BLOCK_MATRIX, FINACV, COLACV, MIDACV, &
    SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
    NCOLCT, CT, DIAG_SCALE_PRES, CT_RHS, FINDCT, COLCT, &
    CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
    CV_ELE_TYPE,  &
    NPHASE, &
    CV_NLOC, U_NLOC, X_NLOC,  &
    CV_NDGLN, X_NDGLN, U_NDGLN, &
    CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
    X, Y, Z,  &
    NU, NV, NW, NUOLD, NVOLD, NWOLD, UG, VG, WG, &
    T, TOLD, DEN, DENOLD, IDIVID_BY_VOL_FRAC, FEM_VOL_FRAC, &
    MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, &
    T_DISOPT, T_DG_VEL_INT_OPT, DT, T_THETA, SECOND_THETA, T_BETA, &
    SUF_SIG_DIAGTEN_BC, &
    DERIV, P,  &
    T_SOURCE, T_ABSORB, VOLFRA_PORE, &
    NDIM, &
    NCOLM, FINDM, COLM, MIDM, &
    XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
    OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
    DEN_FEMT, &
    IGOT_T2, T2, T2OLD, IGOT_THETA_FLUX, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
    THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, THETA_GDIFF, &
    IN_ELE_UPWIND, DG_ELE_UPWIND, &
    NOIT_DIM, &
    MEAN_PORE_CV, &
    THERMAL, &
    mass_ele_transp, &
    option_path, StorageIndexes, Field_selector )

        ! Solve for internal energy using a control volume method.

        implicit none
        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type(tensor_field), intent(inout) :: tracer
        type(tensor_field), intent(in) :: velocity, density
        type(scalar_field), intent(in) :: pressure

        INTEGER, intent( in ) :: NCOLACV, NCOLCT, CV_NONODS, U_NONODS, X_NONODS, MAT_NONODS, TOTELE, &
        CV_ELE_TYPE, NPHASE, CV_NLOC, U_NLOC, X_NLOC,  MAT_NLOC, &
        CV_SNLOC, U_SNLOC, STOTEL, XU_NLOC, NDIM, NCOLM, NCOLELE, &
        NOPT_VEL_UPWIND_COEFS, &
        IGOT_T2, IGOT_THETA_FLUX, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND, IDIVID_BY_VOL_FRAC, Field_selector

        LOGICAL, intent( in ) :: GET_THETA_FLUX, USE_THETA_FLUX, THERMAL
        INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: U_SNDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: FINACV
        INTEGER, DIMENSION( : ), intent( in ) :: COLACV
        INTEGER, DIMENSION( : ), intent( in ) :: MIDACV
        INTEGER, DIMENSION( : ), intent( in ) :: SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV
        INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
        INTEGER, DIMENSION( : ), intent( in ) :: COLCT
        REAL, DIMENSION( : ), allocatable, intent( inout ) :: ACV
        REAL, DIMENSION( :, :, : ), intent( inout ) :: DENSE_BLOCK_MATRIX
        REAL, DIMENSION( :, :, : ), intent( inout ) :: CV_RHS
        REAL, DIMENSION( : ), intent( inout ) :: DIAG_SCALE_PRES
        REAL, DIMENSION( : ), intent( inout ) :: CT_RHS
        REAL, DIMENSION( :, :, : ), intent( inout ) :: CT
        REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
        REAL, DIMENSION( : ), intent( in ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD, UG, VG, WG
        REAL, DIMENSION( : ), intent( inout ) :: T, DEN_FEMT
        REAL, DIMENSION( :), intent( in ) :: TOLD
        REAL, DIMENSION( :, : ), intent( in ) :: DEN, DENOLD
        REAL, DIMENSION( :, : ), intent( in ) :: FEM_VOL_FRAC
        REAL, DIMENSION( : ), intent( in ) :: T2, T2OLD
        REAL, DIMENSION( :, : ), intent( inout ) :: THETA_GDIFF
        REAL, DIMENSION(:, : ),  intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J
        REAL, DIMENSION( :, :, :, : ), intent( in ) :: TDIFFUSION
        INTEGER, intent( in ) :: T_DISOPT, T_DG_VEL_INT_OPT
        REAL, intent( in ) :: DT, T_THETA
        REAL, intent( in ) :: T_BETA
        REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
        REAL, DIMENSION( NPHASE, CV_NONODS ), intent( in ) :: DERIV
        REAL, DIMENSION( : ), intent( in ) :: P
        REAL, DIMENSION( : ), intent( in ) :: T_SOURCE
        REAL, DIMENSION( :, :, : ), intent( in ) :: T_ABSORB
        REAL, DIMENSION( : ), intent( in ) :: VOLFRA_PORE
        INTEGER, DIMENSION( : ), intent( in ) :: FINDM
        INTEGER, DIMENSION( : ), intent( in ) :: COLM
        INTEGER, DIMENSION( : ), intent( in ) :: MIDM
        INTEGER, DIMENSION( : ), intent( in ) :: FINELE
        INTEGER, DIMENSION( : ), intent( in ) :: COLELE
        REAL, DIMENSION( : ), intent( in ) :: OPT_VEL_UPWIND_COEFS
        INTEGER, INTENT( IN ) :: NOIT_DIM
        REAL, DIMENSION( : ), intent( inout ) :: MEAN_PORE_CV
        real, dimension( : ), intent( inout ) :: mass_ele_transp
        character( len = * ), intent( in ), optional :: option_path
        integer, dimension(:), intent(inout) :: StorageIndexes
        ! Local variables
        LOGICAL, PARAMETER :: GETCV_DISC = .TRUE., GETCT= .FALSE., RETRIEVE_SOLID_CTY= .FALSE.
        INTEGER :: ITS_FLUX_LIM, IGOT_THERM_VIS
        INTEGER :: NCOLACV_SUB, IPHASE, I, J
        REAL :: SECOND_THETA
        INTEGER :: STAT,U_ELE_TYPE
        LOGICAL :: CV_METHOD
        character( len = option_path_len ) :: path
      
        REAL, DIMENSION( : ), allocatable :: CV_RHS1
        REAL, DIMENSION( :,:,:,: ), allocatable :: THERM_U_DIFFUSION

        allocate(  cv_rhs1( cv_nonods * nphase ) ) ; cv_rhs1=0.0

        SECOND_THETA = 1.0
        U_ELE_TYPE = CV_ELE_TYPE
        path='/material_phase[0]/scalar_field::Temperature/prognostic/temporal_discretisation/control_volumes/second_theta'
        call get_option( path, second_theta, stat )
        CV_METHOD = .FALSE.

        IF(CV_METHOD) THEN ! cv method...

            IGOT_THERM_VIS=0
            ALLOCATE( THERM_U_DIFFUSION(NDIM,NDIM,NPHASE,MAT_NONODS*IGOT_THERM_VIS ) )

            CALL CV_ASSEMB( state, packed_state, &
                 tracer, velocity, density, &
            CV_RHS1, &
            NCOLACV, ACV, DENSE_BLOCK_MATRIX, FINACV, COLACV, MIDACV, &
            SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
            NCOLCT, CT, DIAG_SCALE_PRES, CT_RHS, FINDCT, COLCT, &
            CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
            CV_ELE_TYPE,  &
            NPHASE, &
            CV_NLOC, U_NLOC, X_NLOC,  &
            CV_NDGLN, X_NDGLN, U_NDGLN, &
            CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
            DEN, DENOLD, &
            MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, IGOT_THERM_VIS, THERM_U_DIFFUSION, &
            T_DISOPT, T_DG_VEL_INT_OPT, DT, T_THETA, SECOND_THETA, T_BETA, &
            SUF_SIG_DIAGTEN_BC, &
            DERIV, P, &
            T_SOURCE, T_ABSORB, VOLFRA_PORE, &
            NDIM, GETCV_DISC, GETCT, &
            NCOLM, FINDM, COLM, MIDM, &
            XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
            OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
            DEN_FEMT, &
            IGOT_T2, T2, T2OLD, IGOT_THETA_FLUX, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
            THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, THETA_GDIFF, &
            IN_ELE_UPWIND, DG_ELE_UPWIND, &
            NOIT_DIM, &
            MEAN_PORE_CV, &
            FINACv, COLACV, NCOLACV, ACV, THERMAL, RETRIEVE_SOLID_CTY, &
            mass_ele_transp , &
            StorageIndexes, Field_selector)!-1,  T_input = T, TOLD_input=TOLD )

        ELSE ! this is for DG...

            !TEMPORAL
            allocate(ACV(NCOLACV))


            CALL WRAPPER_ASSEMB_FORCE_CTY( state, packed_state, &
                 velocity,pressure, &
            NDIM, NPHASE, U_NLOC, X_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
            U_ELE_TYPE, CV_ELE_TYPE, &
            U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
            U_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
            STOTEL, U_SNDGLN, CV_SNDGLN, U_SNLOC, CV_SNLOC, &
            X, Y, Z, T_ABSORB, T_SOURCE, TDIFFUSION, &
            T, TOLD, &
            NU, NV, NW, NUOLD, NVOLD, NWOLD, &
            DEN, DENOLD, IDIVID_BY_VOL_FRAC, FEM_VOL_FRAC, &
            DT, &
            CV_RHS, &
            ACV, NCOLACV, FINACV, COLACV, & ! Force balance sparsity
            NCOLELE, FINELE, COLELE, & ! Element connectivity.
            XU_NLOC, XU_NDGLN, &
            option_path,&
            StorageIndexes=StorageIndexes )

        ENDIF

    END SUBROUTINE CV_ASSEMB_CV_DG




    SUBROUTINE WRAPPER_ASSEMB_FORCE_CTY( state, packed_state,&
         velocity,pressure, &
    NDIM, NPHASE, U_NLOC, X_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
    U_ELE_TYPE, CV_ELE_TYPE, &
    U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
    U_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
    STOTEL, U_SNDGLN, CV_SNDGLN, U_SNLOC, CV_SNLOC, &
    X, Y, Z, T_ABSORB, T_SOURCE, TDIFFUSION, &
    T, TOLD, &
    U, V, W, UOLD, VOLD, WOLD, &
    DEN_ALL, DENOLD_ALL, IDIVID_BY_VOL_FRAC, FEM_VOL_FRAC, &
    DT, &
    CV_RHS, &
    ACV, NCOLACV, FINACV, COLACV, & ! Force balance sparsity
    NCOLELE, FINELE, COLELE, & ! Element connectivity.
    XU_NLOC, XU_NDGLN, &
    option_path,&
    StorageIndexes )
        use shape_functions_NDim
        implicit none

        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type(tensor_field), intent(in) :: velocity
        type(scalar_field), intent(in) :: pressure
        INTEGER, intent( in ) :: NDIM, NPHASE, U_NLOC, X_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
        U_ELE_TYPE, CV_ELE_TYPE, U_NONODS, CV_NONODS, X_NONODS, &
        MAT_NONODS, STOTEL, U_SNLOC, CV_SNLOC, &
        NCOLACV, NCOLELE, XU_NLOC, IDIVID_BY_VOL_FRAC
        INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
        INTEGER, DIMENSION( : ), intent( in )  :: CV_NDGLN
        INTEGER, DIMENSION( : ), intent( in )  :: X_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: U_SNDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
        REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
        REAL, DIMENSION( :, :, : ), intent( in ) :: T_ABSORB
        REAL, DIMENSION( : ), intent( in ) :: T_SOURCE
        REAL, DIMENSION( : ), intent( in ) :: U, V, W, UOLD, VOLD, WOLD
        REAL, DIMENSION( : ), intent( in ) :: T, TOLD
        REAL, DIMENSION( :, : ), intent( in ) :: DEN_ALL, DENOLD_ALL
        REAL, intent( in ) :: DT
        REAL, DIMENSION( :, :, : ), intent( inout ) :: CV_RHS
        REAL, DIMENSION( : ),  intent( inout ) :: ACV
        INTEGER, DIMENSION( : ), intent( in ) :: FINACV
        INTEGER, DIMENSION( : ), intent( in ) :: COLACV
        INTEGER, DIMENSION( : ), intent( in ) :: FINELE
        INTEGER, DIMENSION( : ), intent( in ) :: COLELE
        REAL, DIMENSION( :, :, :, : ), intent( in ) :: TDIFFUSION
        REAL, DIMENSION( :, : ), intent( in ) :: FEM_VOL_FRAC
        character( len = * ), intent( in ), optional :: option_path
        integer, dimension(:), intent(inout) :: StorageIndexes
        ! Local  variables... none
        REAL, DIMENSION ( :, :, : ), allocatable :: RZERO, T_IN, TOLD_IN
        REAL, DIMENSION ( : ), allocatable :: RDUM
        REAL, DIMENSION ( :, : ), allocatable :: RDUM2
        REAL, DIMENSION ( :, :, : ), allocatable :: RDUM3
        INTEGER, DIMENSION ( : ), allocatable :: IDUM,IZERO

        INTEGER :: IPLIKE_GRAD_SOU, NDIM_IN, NPHASE_IN, X_ILOC, X_INOD, MAT_INOD, S, E
        LOGICAL :: JUST_BL_DIAG_MAT

        INTEGER :: U_NLOC2, ILEV, NLEV, ELE, U_ILOC, U_INOD, IPHASE, IDIM, I, SELE, U_SILOC
        REAL, DIMENSION( :, :, : ), allocatable :: U_ALL, UOLD_ALL, T_SOURCE_ALL, T_ABSORB_ALL
        REAL, DIMENSION( :, : ), allocatable :: X_ALL

        INTEGER, DIMENSION( :, : ), allocatable :: IZERO2


        ndim_in = 1 ; nphase_in = 1

        ALLOCATE( U_ALL( NDIM, NPHASE, U_NONODS ), UOLD_ALL( NDIM, NPHASE, U_NONODS ), &
        X_ALL( NDIM, X_NONODS ) )
        U_ALL = 0. ; UOLD_ALL = 0. ; X_ALL = 0.

        IF(U_NLOC.NE.CV_NLOC) THEN
            ewrite(3,*) 'u_nloc, cv_nloc:', u_nloc, cv_nloc
            FLAbort( 'Only working for u_nloc == cv_nloc ' )
        END IF

        ALLOCATE(RZERO(TOTELE , U_NLOC * NPHASE * NDIM , U_NLOC * NPHASE * NDIM))
        RZERO=0.0
        ALLOCATE(IZERO(TOTELE * U_NLOC * NPHASE * NDIM * U_NLOC * NPHASE * NDIM))
        IZERO=0
        ALLOCATE(RDUM(TOTELE *  U_NLOC * NPHASE * NDIM  *  U_NLOC * NPHASE * NDIM))
        RDUM=0.0
        ALLOCATE(IDUM(TOTELE * U_NLOC * NPHASE * NDIM * U_NLOC * NPHASE * NDIM))
        IDUM=0
        ALLOCATE(RDUM2(NPHASE,CV_NONODS))
        RDUM2=0.0

        IPLIKE_GRAD_SOU=0

        IF ( IS_OVERLAPPING ) THEN
            NLEV = CV_NLOC
            U_NLOC2 = MAX( 1, U_NLOC / CV_NLOC )
        ELSE
            NLEV = 1
            U_NLOC2 = U_NLOC
        END IF
        DO ELE = 1, TOTELE
            DO ILEV = 1, NLEV
                DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
                    U_INOD = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_ILOC )
                    DO IPHASE = 1, NPHASE
                        DO IDIM = 1, NDIM
                            IF ( IDIM==1 ) THEN
                                U_ALL( IDIM, IPHASE, U_INOD ) = U( U_INOD + (IPHASE-1)*U_NONODS )
                                UOLD_ALL( IDIM, IPHASE, U_INOD ) = UOLD( U_INOD + (IPHASE-1)*U_NONODS )
                            ELSE IF ( IDIM==2 ) THEN
                                U_ALL( IDIM, IPHASE, U_INOD ) = V( U_INOD + (IPHASE-1)*U_NONODS )
                                UOLD_ALL( IDIM, IPHASE, U_INOD ) = VOLD( U_INOD + (IPHASE-1)*U_NONODS )
                            ELSE
                                U_ALL( IDIM, IPHASE, U_INOD ) = W( U_INOD + (IPHASE-1)*U_NONODS )
                                UOLD_ALL( IDIM, IPHASE, U_INOD ) = WOLD( U_INOD + (IPHASE-1)*U_NONODS )
                            END IF
                        END DO
                    END DO
                END DO
            END DO
        END DO
        DO IDIM = 1, NDIM
            IF ( IDIM==1 ) THEN
                X_ALL( IDIM, : ) = X
            ELSE IF ( IDIM==2 ) THEN
                X_ALL( IDIM, : ) = Y
            ELSE
                X_ALL( IDIM, : ) = Z
            END IF
        END DO

        ALLOCATE( T_SOURCE_ALL( NDIM_IN, NPHASE_IN, U_NONODS ) )
        DO IPHASE = 1, NPHASE
            DO IDIM = 1, NDIM_IN
                S = 1 + (IDIM-1)*U_NONODS + (IPHASE-1)*NDIM_IN*U_NONODS
                E = IDIM*U_NONODS + (IPHASE-1)*NDIM_IN*U_NONODS
                T_SOURCE_ALL( IDIM, IPHASE, : ) = T_SOURCE( S:E )
            END DO
        END DO

        ALLOCATE( T_ABSORB_ALL( NDIM_IN * NPHASE_IN, NDIM_IN * NPHASE_IN, MAT_NONODS ) )
        DO MAT_INOD = 1, MAT_NONODS
            T_ABSORB_ALL( :, :, MAT_INOD ) = T_ABSORB( MAT_INOD, :, : )
        END DO




        allocate( t_in( ndim_in, nphase_in, u_nonods ) ) ; t_in(1,1,:) = t
        allocate( told_in( ndim_in, nphase_in, u_nonods ) ) ; told_in(1,1,:) = t

        ALLOCATE( IZERO2( NPHASE_IN,STOTEL ) ) ; IZERO2 = 0
        ALLOCATE( RDUM3( NPHASE_IN,CV_SNLOC,STOTEL ) ) ; RDUM3 = 0.

        CALL ASSEMB_FORCE_CTY( state, packed_state, &
             velocity,pressure, &
        NDIM, NPHASE_IN, U_NLOC, X_NLOC, CV_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
        U_ELE_TYPE, CV_ELE_TYPE, &
        U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
        U_NDGLN, CV_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
        STOTEL, U_SNDGLN, CV_SNDGLN, CV_SNDGLN, U_SNLOC, CV_SNLOC, CV_SNLOC, &
        X_ALL, RZERO, T_ABSORB_ALL, T_SOURCE_ALL, RDUM3, &
        T_IN, TOLD_IN, &
        U_ALL, UOLD_ALL, &
        DEN_ALL, DENOLD_ALL, IDIVID_BY_VOL_FRAC, FEM_VOL_FRAC, &
        DT, &
        CV_RHS, &
        RDUM3, 0, IDUM, IDUM, &
        ACV, NCOLACV, FINACV, COLACV, &! Force balance sparsity
        NCOLELE, FINELE, COLELE, & ! Element connectivity.
        XU_NLOC, XU_NDGLN, &
        RZERO, JUST_BL_DIAG_MAT,  &
        TDIFFUSION, DEN_ALL, DENOLD_ALL, .FALSE., & ! TDiffusion need to be obtained down in the tree according to the option_path
        IPLIKE_GRAD_SOU, RDUM2, RDUM2, &
        RDUM, NDIM_IN, &
        StorageIndexes )



        DEALLOCATE( U_ALL, UOLD_ALL, X_ALL, RZERO, IZERO, RDUM, IDUM, T_IN, TOLD_IN )

    END SUBROUTINE WRAPPER_ASSEMB_FORCE_CTY




    SUBROUTINE SIMPLE_SOLVER( CMC, P, RHS,  &
    NCMC, NONODS, FINCMC, COLCMC, MIDCMC,  &
    ERROR, RELAX, RELAX_DIAABS, RELAX_DIA, N_LIN_ITS )
        !
        ! Solve CMC * P = RHS for RHS.
        ! RELAX: overall relaxation coeff; =1 for no relaxation.
        ! RELAX_DIAABS: relaxation of the absolute values of the sum of the row of the matrix;
        !               - recommend >=2 for hard problems, =0 for easy
        ! RELAX_DIA: relaxation of diagonal; =1 no relaxation (normally applied).
        ! N_LIN_ITS = no of linear iterations
        ! ERROR= solver tolerence between 2 consecutive iterations
        implicit none
        REAL, intent( in ) :: ERROR, RELAX, RELAX_DIAABS, RELAX_DIA
        INTEGER, intent( in ) ::  N_LIN_ITS, NCMC, NONODS
        REAL, DIMENSION( : ), intent( in ) ::  CMC
        REAL, DIMENSION( : ), intent( inout ) ::  P
        REAL, DIMENSION( : ), intent( in ) :: RHS
        INTEGER, DIMENSION( : ), intent( in ) :: FINCMC
        INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
        INTEGER, DIMENSION( : ), intent( in ) :: MIDCMC
        ! Local variables
        INTEGER :: ITS, ILOOP, ISTART, IFINI, ISTEP, NOD, COUNT
        REAL :: R, SABS_DIAG, RTOP, RBOT, POLD, MAX_ERR

        ewrite(3,*) 'In Solver'

        Loop_Non_Linear_Iter: DO ITS = 1, N_LIN_ITS

            MAX_ERR = 0.0
            Loop_Internal: DO ILOOP = 1, 2
                IF( ILOOP == 1 ) THEN
                    ISTART = 1
                    IFINI = NONODS
                    ISTEP = 1
                ELSE
                    ISTART = NONODS
                    IFINI = 1
                    ISTEP = -1
                ENDIF

                Loop_Nods: DO NOD = ISTART, IFINI, ISTEP
                    R = RELAX_DIA * CMC( MIDCMC( NOD )) * P( NOD ) + RHS( NOD )
                    SABS_DIAG = 0.0
                    DO COUNT = FINCMC( NOD ), FINCMC( NOD + 1 ) - 1
                        R = R - CMC( COUNT ) * P( COLCMC( COUNT ))
                        SABS_DIAG = SABS_DIAG + ABS( CMC( COUNT ))
                    END DO
                    RTOP = R + RELAX_DIAABS * SABS_DIAG * P( NOD )
                    RBOT = RELAX_DIAABS * SABS_DIAG + RELAX_DIA * CMC( MIDCMC( NOD ))
                    POLD = P( NOD )
                    P( NOD ) = RELAX * ( RTOP / RBOT ) + ( 1.0 - RELAX ) * P( NOD )
                    MAX_ERR = MAX( MAX_ERR, ABS( POLD - P( NOD )))
                END DO Loop_Nods
            END DO Loop_Internal

            IF( MAX_ERR < ERROR ) CYCLE

        END DO Loop_Non_Linear_Iter

        ewrite(3,*) 'Leaving Solver'

        RETURN
    END SUBROUTINE SIMPLE_SOLVER



    subroutine VolumeFraction_Assemble_Solve( state,packed_state, &
    NCOLACV, FINACV, COLACV, MIDACV, &
    SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
    block_to_global_acv, global_dense_block_acv, &
    NCOLCT, FINDCT, COLCT, &
    CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
    CV_ELE_TYPE,  &
    NPHASE, &
    CV_NLOC, U_NLOC, X_NLOC, &
    CV_NDGLN, X_NDGLN, U_NDGLN, &
    CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
    MAT_NLOC,MAT_NDGLN,MAT_NONODS, &
    V_DISOPT, V_DG_VEL_INT_OPT, DT, V_THETA, V_BETA, &
    SUF_SIG_DIAGTEN_BC, &
    DERIV, &
    V_SOURCE, V_ABSORB, VOLFRA_PORE, &
    NDIM, &
    NCOLM, FINDM, COLM, MIDM, &
    XU_NLOC, XU_NDGLN ,FINELE, COLELE, NCOLELE, &
    OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
    DEN_FEMT, &
    igot_theta_flux, SCVNGI_THETA, USE_THETA_FLUX, &
    IN_ELE_UPWIND, DG_ELE_UPWIND, &
    NOIT_DIM, &
    option_path, &
    mass_ele_transp,&
    THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, &
    StorageIndexes)

        implicit none
        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ) :: packed_state
        INTEGER, intent( in ) :: NCOLACV, NCOLCT, &
        CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
        CV_ELE_TYPE, &
        NPHASE, CV_NLOC, U_NLOC, X_NLOC, &
        CV_SNLOC, U_SNLOC, STOTEL, XU_NLOC, NDIM, &
        NCOLM, NCOLELE, NOPT_VEL_UPWIND_COEFS, &
        MAT_NLOC, MAT_NONODS, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND,igot_theta_flux
        LOGICAL, intent( in ) :: USE_THETA_FLUX
        INTEGER, DIMENSION(: ), intent( in ) :: CV_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: U_SNDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: FINACV
        INTEGER, DIMENSION( : ), intent( in ) :: COLACV
        INTEGER, DIMENSION( : ), intent( in ) :: MIDACV
        integer, dimension(:), intent(in)  :: small_finacv,small_colacv,small_midacv
        integer, dimension(:), intent(in)  :: block_to_global_acv
        integer, dimension(:,:), intent(in) :: global_dense_block_acv
        INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
        INTEGER, DIMENSION( : ), intent( in ) :: COLCT
        REAL, DIMENSION( : ), intent( inout ) :: DEN_FEMT
        REAL, DIMENSION( :, :), intent( inout ), optional :: THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J
        INTEGER, intent( in ) :: V_DISOPT, V_DG_VEL_INT_OPT
        REAL, intent( in ) :: DT, V_THETA
        REAL, intent( inout ) :: V_BETA
        REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
        REAL, DIMENSION( NPHASE, CV_NONODS ), intent( in ) :: DERIV
        REAL, DIMENSION( : ), intent( in ) :: V_SOURCE
        REAL, DIMENSION( :, :, : ), intent( in ) :: V_ABSORB
        REAL, DIMENSION( : ), intent( in ) :: VOLFRA_PORE
        INTEGER, DIMENSION( : ), intent( in ) :: FINDM
        INTEGER, DIMENSION( : ), intent( in ) :: COLM
        INTEGER, DIMENSION( : ), intent( in ) :: MIDM
        INTEGER, DIMENSION( : ), intent( in ) :: FINELE
        INTEGER, DIMENSION( : ), intent( in ) :: COLELE
        REAL, DIMENSION( : ), intent( in ) :: OPT_VEL_UPWIND_COEFS
        INTEGER, INTENT( IN ) :: NOIT_DIM
        character(len= * ), intent(in), optional :: option_path
        real, dimension( : ), intent( inout ) :: mass_ele_transp
        integer, dimension(:), intent(inout) :: StorageIndexes
        ! Local Variables
        LOGICAL, PARAMETER :: THERMAL= .false.
        integer :: nits_flux_lim, its_flux_lim, igot_t2
        REAL, DIMENSION( : ), allocatable :: ACV, mass_mn_pres, block_ACV, CV_RHS, DIAG_SCALE_PRES, CT_RHS
        REAL, DIMENSION( :,:,: ), allocatable :: dense_block_matrix, CT
        REAL, DIMENSION( :,:,:,: ), allocatable :: TDIFFUSION
        REAL, DIMENSION( :, : ), allocatable :: THETA_GDIFF, DEN_ALL, DENOLD_ALL
        REAL, DIMENSION( : ), allocatable :: T2, T2OLD, MEAN_PORE_CV
        REAL, DIMENSION( : ), allocatable :: DENSITY_OR_ONE, DENSITYOLD_OR_ONE
        REAL, DIMENSION( :, :, :, : ), allocatable :: THERM_U_DIFFUSION
        LOGICAL :: GET_THETA_FLUX
        REAL :: SECOND_THETA
        INTEGER :: STAT, i,j, IGOT_THERM_VIS
        character( len = option_path_len ) :: path
        LOGICAL, PARAMETER :: GETCV_DISC = .TRUE., GETCT= .FALSE., RETRIEVE_SOLID_CTY= .FALSE.
        real, dimension(:), allocatable :: X
        type( tensor_field ), pointer :: den_all2, denold_all2
        !type( scalar_field ), pointer :: p
        !Working pointers
        real, dimension(:), pointer :: p
        real, dimension(:,:), pointer :: satura,saturaold
        type(tensor_field), pointer :: tracer, velocity, density

        call get_var_from_packed_state(packed_state,FEPressure = P,&
        PhaseVolumeFraction = satura,OldPhaseVolumeFraction = saturaold)
        GET_THETA_FLUX = .FALSE.
        IGOT_T2 = 0

        ALLOCATE( T2( CV_NONODS * NPHASE * IGOT_T2 ))
        ALLOCATE( T2OLD( CV_NONODS * NPHASE * IGOT_T2 ))
        ALLOCATE( THETA_GDIFF( NPHASE * IGOT_T2, CV_NONODS * IGOT_T2 ))

        ewrite(3,*) 'In VOLFRA_ASSEM_SOLVE'

        ALLOCATE( ACV( NCOLACV ) ) ; ACV = 0.
        ALLOCATE( block_ACV( size(block_to_global_acv) ) ) ; block_ACV = 0.
        ALLOCATE( mass_mn_pres(size(small_colacv)) ) ; mass_mn_pres = 0.
        ALLOCATE( dense_block_matrix( nphase , nphase , cv_nonods) ); dense_block_matrix=0;
        ALLOCATE( CV_RHS( CV_NONODS * NPHASE ) ) ; CV_RHS = 0.
        ALLOCATE( CT( NDIM,NPHASE,NCOLCT ) )
        ALLOCATE( DIAG_SCALE_PRES( CV_NONODS ) )
        ALLOCATE( CT_RHS( CV_NONODS ) )
        ALLOCATE( TDIFFUSION( MAT_NONODS, NDIM, NDIM, NPHASE ) )
        ALLOCATE( MEAN_PORE_CV( CV_NONODS ) )



        ALLOCATE( DEN_ALL( NPHASE, CV_NONODS  ), DENOLD_ALL( NPHASE, CV_NONODS ) )
        IF ( IGOT_THETA_FLUX == 1 ) THEN
            ! use DEN=1 because the density is already in the theta variables
            DEN_ALL=1.0 ; DENOLD_ALL=1.0
        ELSE
            DEN_ALL2 => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedDensity" )
            DENOLD_ALL2 => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedOldDensity" )
            DEN_ALL = DEN_ALL2%VAL( 1, :, : ) ; DENOLD_ALL = DENOLD_ALL2%VAL( 1, :, : )
        END IF


        TDIFFUSION = 0.0
        V_BETA = 1.0

        path = '/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/temporal_discretisation/' // &
        'control_volumes/'
        call get_option( trim( path ) // 'second_theta', second_theta, stat , default = 1.0)
        call get_option( trim( path ) // 'number_advection_iterations', nits_flux_lim, default = 1 )

        ! THIS DOES NOT WORK FOR NITS_FLUX_LIM>1 (NOBODY KNOWS WHY)
        IGOT_THERM_VIS=0
        ALLOCATE( THERM_U_DIFFUSION(NDIM,NDIM,NPHASE,MAT_NONODS*IGOT_THERM_VIS ) )

!         p => extract_scalar_field( packed_state, "FEPressure" )
        allocate(X(size(CV_RHS,1)))

        tracer=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
        velocity=>extract_tensor_field(packed_state,"PackedVelocity")
        density=>extract_tensor_field(packed_state,"PackedDensity")

        Loop_NonLinearFlux: DO ITS_FLUX_LIM = 1, 1 !nits_flux_lim

            call CV_ASSEMB( state, packed_state, &
            tracer, velocity, density, &
            CV_RHS, &
            NCOLACV, block_acv, DENSE_BLOCK_MATRIX, FINACV, COLACV, MIDACV, &
            SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV,&
            NCOLCT, CT, DIAG_SCALE_PRES, CT_RHS, FINDCT, COLCT, &
            CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
            CV_ELE_TYPE,  &
            NPHASE, &
            CV_NLOC, U_NLOC, X_NLOC, &
            CV_NDGLN, X_NDGLN, U_NDGLN, &
            CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
            DEN_ALL, DENOLD_ALL, &
            MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, IGOT_THERM_VIS, THERM_U_DIFFUSION, &
            V_DISOPT, V_DG_VEL_INT_OPT, DT, V_THETA, SECOND_THETA, V_BETA, &
            SUF_SIG_DIAGTEN_BC, &
            DERIV, P, &
            V_SOURCE, V_ABSORB, VOLFRA_PORE, &
            NDIM, GETCV_DISC, GETCT, &
            NCOLM, FINDM, COLM, MIDM, &
            XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
            OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
            DEN_FEMT, &
            IGOT_T2, T2, T2OLD, igot_theta_flux, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
            THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, THETA_GDIFF, &
            IN_ELE_UPWIND, DG_ELE_UPWIND, &
            NOIT_DIM, &
            MEAN_PORE_CV, &
            SMALL_FINACV, SMALL_COLACV, size(small_colacv), mass_Mn_pres, THERMAL, RETRIEVE_SOLID_CTY, &
            mass_ele_transp,&
            StorageIndexes, 3 )

!            satura=0.0 !saturaold([([(i+(j-1)*cv_nonods,j=1,nphase)],i=1,cv_nonods)])

            X = 0.
            call assemble_global_multiphase_csr(acv,&
            block_acv,dense_block_matrix,&
            block_to_global_acv,global_dense_block_acv)

            CALL SOLVER( ACV, X, CV_RHS, &
            FINACV, COLACV, &
            trim(option_path) )

            !Copy to phaseVolumeFraction in packed_state
            do j = 1, cv_nonods
                satura(:,j) = x(1+(j-1)*NPHASE : j*NPHASE)
            end do

!            satura([([(i+(j-1)*cv_nonods,j=1,nphase)],i=1,cv_nonods)])=satura
        END DO Loop_NonLinearFlux
        deallocate(X)
        !Set saturation to be between bounds
        satura = min(max(satura,0.0), 1.0)

        DEALLOCATE( ACV )
        DEALLOCATE( mass_mn_pres )
        deallocate( block_acv )
        deallocate( dense_block_matrix )
        DEALLOCATE( CV_RHS )
        DEALLOCATE( CT )
        DEALLOCATE( DIAG_SCALE_PRES )
        DEALLOCATE( CT_RHS )
        DEALLOCATE( TDIFFUSION )
        DEALLOCATE( T2 )
        DEALLOCATE( T2OLD )
        DEALLOCATE( THETA_GDIFF )

        ewrite(3,*) 'Leaving VOLFRA_ASSEM_SOLVE'

        RETURN
    end subroutine VolumeFraction_Assemble_Solve


    SUBROUTINE FORCE_BAL_CTY_ASSEM_SOLVE( state, packed_state, &
         velocity,pressure, &
    NDIM, NPHASE, NCOMP, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
    U_ELE_TYPE, P_ELE_TYPE, &
    U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
    U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
    STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
    U_SNLOC, P_SNLOC, CV_SNLOC, &
    U_ABS_STAB, MAT_ABSORB, U_ABSORBIN, U_SOURCE, U_SOURCE_CV, &
    DERIV, IDIVID_BY_VOL_FRAC, FEM_VOL_FRAC, &
    DT, &
    NCOLC, FINDC, COLC, & ! C sparcity - global cty eqn
    NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, MIDDGM_PHA, &! Force balance sparcity
    NCOLELE, FINELE, COLELE, & ! Element connectivity.
    NCOLCMC, FINDCMC, COLCMC, MIDCMC, & ! pressure matrix for projection method
    NCOLACV, FINACV, COLACV, MIDACV, & ! For CV discretisation method
    NCOLSMALL,SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
    NLENMCY, NCOLMCY, FINMCY, COLMCY, MIDMCY, & ! Force balance plus cty multi-phase eqns
    NCOLCT, FINDCT, COLCT, & ! CT sparcity - global cty eqn.
    CV_ELE_TYPE, &
    V_DISOPT, V_DG_VEL_INT_OPT, V_THETA, &
    SUF_SIG_DIAGTEN_BC, &
    V_SOURCE, V_ABSORB, VOLFRA_PORE, &
    NCOLM, FINDM, COLM, MIDM, & ! Sparsity for the CV-FEM
    XU_NLOC, XU_NDGLN, &
    UDIFFUSION, &
    OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
    IGOT_THETA_FLUX, SCVNGI_THETA, USE_THETA_FLUX, &
    THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, &
    IN_ELE_UPWIND, DG_ELE_UPWIND, &
    NOIT_DIM, &
    IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, &
    scale_momentum_by_volume_fraction, &
    StorageIndexes )

        IMPLICIT NONE
        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type( tensor_field ), intent(in) :: velocity
        type( scalar_field ), intent(in) :: pressure
        INTEGER, intent( in ) :: NDIM, NPHASE, NCOMP, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, &
        TOTELE, U_ELE_TYPE, P_ELE_TYPE, &
        U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
        STOTEL, U_SNLOC, P_SNLOC, &
        CV_SNLOC, &
        NCOLC, NCOLDGM_PHA, NCOLELE, NCOLCMC, NCOLACV, ncolsmall, NLENMCY, NCOLMCY, NCOLCT, &
        CV_ELE_TYPE, V_DISOPT, V_DG_VEL_INT_OPT, NCOLM, XU_NLOC, &
        NOPT_VEL_UPWIND_COEFS, IGOT_THETA_FLUX, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND, &
        IPLIKE_GRAD_SOU, IDIVID_BY_VOL_FRAC
        LOGICAL, intent( in ) :: USE_THETA_FLUX, scale_momentum_by_volume_fraction
        INTEGER, DIMENSION(  :  ), intent( in ) :: U_NDGLN
        INTEGER, DIMENSION(  :  ), intent( in ) :: P_NDGLN
        INTEGER, DIMENSION(  :  ), intent( in ) :: CV_NDGLN
        INTEGER, DIMENSION(  :  ), intent( in ) :: X_NDGLN
        INTEGER, DIMENSION(  :  ), intent( in ) :: MAT_NDGLN
        INTEGER, DIMENSION(  :  ), intent( in ) :: U_SNDGLN
        INTEGER, DIMENSION(  :  ), intent( in ) :: P_SNDGLN

        INTEGER, DIMENSION(  : ), intent( in ) :: CV_SNDGLN
        INTEGER, DIMENSION(  : ), intent( in ) :: XU_NDGLN
        REAL, DIMENSION(  :, :, :  ), intent( inout ) :: U_ABS_STAB, U_ABSORBIN, MAT_ABSORB
        REAL, DIMENSION(  :  ), intent( in ) :: U_SOURCE
        REAL, DIMENSION(  :  ), intent( inout ) :: U_SOURCE_CV

!        REAL, DIMENSION(  :  ), intent( in ) :: SATURAOLD
!        REAL, DIMENSION(  :  ), intent( inout ) :: SATURA
        REAL, DIMENSION(  NPHASE, CV_NONODS ), intent( in ) :: DERIV
        REAL, DIMENSION(  NPHASE*IDIVID_BY_VOL_FRAC, CV_NONODS *IDIVID_BY_VOL_FRAC ), intent( in ) :: FEM_VOL_FRAC
        REAL, DIMENSION(  : , :  ), intent( in ) :: SUF_SIG_DIAGTEN_BC
        REAL, intent( in ) :: DT
        INTEGER, DIMENSION(  :  ), intent( in ) :: FINDC
        INTEGER, DIMENSION(  :  ), intent( in ) :: COLC
        INTEGER, DIMENSION(  :  ), intent( in ) :: FINDGM_PHA
        INTEGER, DIMENSION(  :  ), intent( in ) :: COLDGM_PHA
        INTEGER, DIMENSION(  :  ), intent( in ) :: MIDDGM_PHA

        INTEGER, DIMENSION(  :  ), intent( in ) :: FINELE
        INTEGER, DIMENSION(  :  ), intent( in ) :: COLELE
        INTEGER, DIMENSION(  :  ), intent( in ) :: FINDCMC
        INTEGER, DIMENSION(  :  ), intent( in ) :: COLCMC
        INTEGER, DIMENSION(  :  ), intent( in ) :: MIDCMC
        INTEGER, DIMENSION(  :  ), intent( in ) :: FINACV
        INTEGER, DIMENSION(  :  ), intent( in ) :: COLACV
        INTEGER, DIMENSION(  :  ), intent( in ) :: MIDACV
        integer, dimension(  :  ), intent( in ) :: small_finacv
        integer, dimension(  :  ), intent( in ) :: small_colacv
        integer, dimension(  :  ), intent( in ) :: small_midacv
        INTEGER, DIMENSION(  :  ), intent( in ) :: FINMCY
        INTEGER, DIMENSION(  :  ), intent( in ) :: COLMCY
        INTEGER, DIMENSION(  :  ), intent( in ) :: MIDMCY
        INTEGER, DIMENSION(  :  ), intent( in ) :: FINDCT
        INTEGER, DIMENSION(  :  ), intent( in ) :: COLCT
        REAL, intent( in ) :: V_THETA
        REAL, DIMENSION(  :  ), intent( in ) :: V_SOURCE
        REAL, DIMENSION(  : ,  : ,: ), intent( in ) :: V_ABSORB
        REAL, DIMENSION(  :  ), intent( in ) :: VOLFRA_PORE
        INTEGER, DIMENSION(  :  ), intent( in ) :: FINDM
        INTEGER, DIMENSION(  :  ), intent( in ) :: COLM
        INTEGER, DIMENSION(  :  ), intent( in ) :: MIDM
        REAL, DIMENSION(  : ,  : ,  : ,  :  ), intent( inout ) :: UDIFFUSION
        REAL, DIMENSION(  :  ), intent( in ) :: OPT_VEL_UPWIND_COEFS
        REAL, DIMENSION( : ,  :  ), intent( inout ) :: &
        THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J
        INTEGER, INTENT( IN ) :: NOIT_DIM
        REAL, DIMENSION( :  ), intent( in ) :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD
        integer, dimension(:), intent(inout) :: StorageIndexes
        ! Local Variables
        LOGICAL, PARAMETER :: GLOBAL_SOLVE = .FALSE.
        ! If IGOT_CMC_PRECON=1 use a sym matrix as pressure preconditioner,=0 else CMC as preconditioner as well.
        INTEGER, PARAMETER :: IGOT_CMC_PRECON = 0
! Gidaspow model B - can use conservative from of
        LOGICAL :: SOLID_FLUID_MODEL_B = .TRUE.
! switch on solid fluid coupling (THE ONLY SWITCH THAT NEEDS TO BE SWITCHED ON FOR SOLID-FLUID COUPLING)...
        LOGICAL :: RETRIEVE_SOLID_CTY = .FALSE.
        character( len = option_path_len ) :: opt

        REAL, DIMENSION( : ), allocatable :: CT_RHS, DIAG_SCALE_PRES, &
        MCY_RHS, MCY, &
        CMC, CMC_PRECON, MASS_MN_PRES, MASS_CV, P_RHS, UP, U_RHS_CDP, DP, &
        UP_VEL, DGM_PHA, DIAG_P_SQRT, ACV
        REAL, DIMENSION( :, :, : ), allocatable :: PIVIT_MAT, C, CDP, CT, U_RHS, DU_VEL, U_RHS_CDP2
        INTEGER :: CV_NOD, COUNT, CV_JNOD, IPHASE, ele, x_nod1, x_nod2, x_nod3, cv_iloc, &
        cv_nod1, cv_nod2, cv_nod3, mat_nod1, u_iloc, u_nod, u_nod_pha, ndpset
        REAL :: der1, der2, der3, uabs, rsum, xc, yc
        LOGICAL :: JUST_BL_DIAG_MAT, NO_MATRIX_STORE, SCALE_P_MATRIX, LINEARISE_DENSITY

        INTEGER :: I, J, IDIM, U_INOD

        !TEMPORARY VARIABLES, ADAPT FROM OLD VARIABLES TO NEW
        INTEGER :: U_NLOC2, ILEV, NLEV, X_ILOC, X_INOD, MAT_INOD, S, E, sele, p_sjloc, u_siloc
        REAL, DIMENSION( :, :, : ), allocatable :: U_ALL, UOLD_ALL, U_SOURCE_ALL, U_SOURCE_CV_ALL, U_ABSORB_ALL, U_ABS_STAB_ALL, U_ABSORB
        REAL, DIMENSION( :, : ), allocatable :: X_ALL, UDEN_ALL, UDENOLD_ALL, DEN_ALL, DENOLD_ALL, PLIKE_GRAD_SOU_COEF_ALL, PLIKE_GRAD_SOU_GRAD_ALL
        REAL, DIMENSION( :, :, :, : ), allocatable :: UDIFFUSION_ALL

        type( tensor_field ), pointer :: u_all2, uold_all2, den_all2, denold_all2
        type( vector_field ), pointer :: x_all2
        type( scalar_field ), pointer :: p_all, cvp_all, Pressure_State, sf, soldf

        real, dimension(:,:), pointer :: den_fem

        ALLOCATE( U_ALL( NDIM, NPHASE, U_NONODS ), UOLD_ALL( NDIM, NPHASE, U_NONODS ), &
        X_ALL( NDIM, X_NONODS ), UDEN_ALL( NPHASE, CV_NONODS ), UDENOLD_ALL( NPHASE, CV_NONODS ), &
        DEN_ALL( NPHASE, CV_NONODS ), DENOLD_ALL( NPHASE, CV_NONODS ) )
        U_ALL = 0. ; UOLD_ALL = 0. ; X_ALL = 0. ; UDEN_ALL = 0. ; UDENOLD_ALL = 0.
        DEN_ALL = 0. ; DENOLD_ALL = 0.

        ewrite(3,*) 'In FORCE_BAL_CTY_ASSEM_SOLVE'

        ALLOCATE( CT( NDIM, NPHASE, NCOLCT )) ; CT=0.
        ALLOCATE( CT_RHS( CV_NONODS )) ; CT_RHS=0.
        ALLOCATE( DIAG_SCALE_PRES( CV_NONODS )) ; DIAG_SCALE_PRES=0.
        ALLOCATE( U_RHS( NDIM, NPHASE, U_NONODS )) ; U_RHS=0.
        ALLOCATE( MCY_RHS( NDIM * NPHASE * U_NONODS + CV_NONODS )) ; MCY_RHS=0.
        ALLOCATE( C( NDIM, NPHASE, NCOLC )) ; C=0.
        ALLOCATE( MCY( NCOLMCY )) ; MCY=0.
        ALLOCATE( CMC( NCOLCMC )) ; CMC=0.
        ALLOCATE( CMC_PRECON( NCOLCMC*IGOT_CMC_PRECON)) ; IF(IGOT_CMC_PRECON.NE.0) CMC_PRECON=0.
        ALLOCATE( MASS_MN_PRES( NCOLCMC )) ;MASS_MN_PRES=0.
        ALLOCATE( MASS_CV( CV_NONODS )) ; MASS_CV=0.
        ALLOCATE( P_RHS( CV_NONODS )) ; P_RHS=0.
        ALLOCATE( UP( NLENMCY )) ; UP=0.
        ALLOCATE( U_RHS_CDP( NDIM * NPHASE * U_NONODS )) ; U_RHS_CDP=0.
        ALLOCATE( U_RHS_CDP2( NDIM, NPHASE, U_NONODS )) ; U_RHS_CDP2=0.

        ALLOCATE( DP( CV_NONODS )) ; DP = 0.
        ALLOCATE( CDP( NDIM, NPHASE, U_NONODS )) ; CDP = 0.
        ALLOCATE( DU_VEL( NDIM,  NPHASE, U_NONODS )) ; DU_VEL = 0.
        ALLOCATE( UP_VEL( NDIM * NPHASE * U_NONODS )) ; UP_VEL = 0.




        ALLOCATE( PIVIT_MAT( NDIM * NPHASE * U_NLOC, NDIM * NPHASE * U_NLOC, TOTELE )) ; PIVIT_MAT=0.0
        ALLOCATE( DGM_PHA( NCOLDGM_PHA )) ; DGM_PHA=0.
        ALLOCATE( ACV( NCOLACV )) ; ACV = 0.

        !################TEMPORARY ADAPT FROM OLD VARIABLES TO NEW###############
    


        U_ALL2 => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedVelocity" )
        U_ALL = U_ALL2%VAL

        UOLD_ALL2 => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedOldVelocity" )
        UOLD_ALL = UOLD_ALL2%VAL

        X_ALL2 => EXTRACT_VECTOR_FIELD( PACKED_STATE, "PressureCoordinate" )
        X_ALL = X_ALL2%VAL

        P_ALL => EXTRACT_SCALAR_FIELD( PACKED_STATE, "FEPressure" )
        CVP_ALL => EXTRACT_SCALAR_FIELD( PACKED_STATE, "CVPressure" )

        linearise_density = have_option( '/material_phase[0]/linearise_density' )

        DEN_ALL2 => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedDensity" )
        DENOLD_ALL2 => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedOldDensity" )

        DEN_ALL = DEN_ALL2%VAL( 1, :, : )
        DENOLD_ALL = DENOLD_ALL2%VAL( 1, :, : )

        !Calculate gravity source terms
        if ( have_option ( '/material_phase[0]/multiphase_properties/relperm_type' ) )then
           UDEN_ALL=0.0; UDENOLD_ALL=0.0
           !Calculate gravity effects for porous media
           !Compact overlapping calculation is after the conversion from old to new
           if (is_overlapping) &
                call calculate_u_source_cv( state, cv_nonods, ndim, nphase, DEN_ALL, U_Source_CV )
        else
           if ( linearise_density ) then
              call linearise_field( DEN_ALL2, UDEN_ALL )
              call linearise_field( DENOLD_ALL2, UDENOLD_ALL )
           else
              UDEN_ALL = DEN_ALL2%VAL( 1, :, : )
              UDENOLD_ALL = DENOLD_ALL2%VAL( 1, :, : )
           end if
            call calculate_u_source_cv( state, cv_nonods, ndim, nphase, uden_all, U_Source_CV )
        end if


         if ( have_option( '/blasting' ) ) then
            RETRIEVE_SOLID_CTY = .true.
            call get_option( '/blasting/Gidaspow_model', opt )
            if ( trim( opt ) == "A" ) SOLID_FLUID_MODEL_B = .false.
         end if


        IF(RETRIEVE_SOLID_CTY) THEN
!        IF(.TRUE.) THEN
! if model B and solid-fluid coupling: 
           sf => EXTRACT_SCALAR_FIELD( PACKED_STATE, "SolidConcentration" )
           soldf => EXTRACT_SCALAR_FIELD( PACKED_STATE, "OldSolidConcentration" )

           IF(SOLID_FLUID_MODEL_B) THEN ! Gidaspow model B - can use conservative from of momentum
              DO IPHASE=1,NPHASE
                 UDEN_ALL(IPHASE,:) = UDEN_ALL(IPHASE,:) * ( 1. - sf%val)
                 UDENOLD_ALL(IPHASE,:) = UDENOLD_ALL(IPHASE,:) * ( 1. - soldf%val)
              END DO
           ENDIF
        ENDIF


        ! calculate the viscosity for the momentum equation...
        !if ( its == 1 ) 
        call calculate_viscosity( state, ncomp, nphase, ndim, mat_nonods, mat_ndgln, uDiffusion )

        ! stabilisation for high aspect ratio problems - switched off
        call calculate_u_abs_stab( U_ABS_STAB, MAT_ABSORB, &
           opt_vel_upwind_coefs, nphase, ndim, totele, cv_nloc, mat_nloc, mat_nonods, mat_ndgln )

        allocate( U_ABSORB( mat_nonods, ndim * nphase, ndim * nphase ) )

        U_ABSORB = U_ABSORBIN + MAT_ABSORB


        ALLOCATE( U_SOURCE_ALL( NDIM, NPHASE, U_NONODS ) )
        ALLOCATE( U_SOURCE_CV_ALL( NDIM, NPHASE, CV_NONODS ) )
        DO IPHASE = 1, NPHASE
            DO IDIM = 1, NDIM
                S = 1 + (IDIM-1)*U_NONODS + (IPHASE-1)*NDIM*U_NONODS
                E = IDIM*U_NONODS + (IPHASE-1)*NDIM*U_NONODS
                U_SOURCE_ALL( IDIM, IPHASE, : ) = U_SOURCE( S:E )

                S = 1 + (IDIM-1)*CV_NONODS + (IPHASE-1)*NDIM*CV_NONODS
                E = IDIM*CV_NONODS + (IPHASE-1)*NDIM*CV_NONODS
                U_SOURCE_CV_ALL( IDIM, IPHASE, : ) = U_SOURCE_CV( S:E )
            END DO
        END DO

        ALLOCATE( U_ABSORB_ALL( NDIM * NPHASE, NDIM * NPHASE, MAT_NONODS ) )
        ALLOCATE( U_ABS_STAB_ALL( NDIM * NPHASE, NDIM * NPHASE, MAT_NONODS ) )
        ALLOCATE( UDIFFUSION_ALL( NDIM, NDIM, NPHASE, MAT_NONODS ) )

        DO MAT_INOD = 1, MAT_NONODS
            U_ABSORB_ALL( :, :, MAT_INOD ) = U_ABSORB( MAT_INOD, :, : )
            U_ABS_STAB_ALL( :, :, MAT_INOD ) = U_ABS_STAB( MAT_INOD, :, : )
            UDIFFUSION_ALL( :, :, :, MAT_INOD ) = UDIFFUSION( MAT_INOD, :, :, : )
        END DO


        ALLOCATE( PLIKE_GRAD_SOU_COEF_ALL( NPHASE, CV_NONODS ) )
        ALLOCATE( PLIKE_GRAD_SOU_GRAD_ALL( NPHASE, CV_NONODS ) )
        DO IPHASE = 1, NPHASE
            PLIKE_GRAD_SOU_COEF_ALL( IPHASE, : ) = PLIKE_GRAD_SOU_COEF( 1 + (IPHASE-1)*CV_NONODS : IPHASE*CV_NONODS )
            PLIKE_GRAD_SOU_GRAD_ALL( IPHASE, : ) = PLIKE_GRAD_SOU_GRAD( 1 + (IPHASE-1)*CV_NONODS : IPHASE*CV_NONODS )
        END DO
        !##########TEMPORARY ADAPT FROM OLD VARIABLES TO NEW############

        !Calculate the RHS for compact_overlapping
        if (is_compact_overlapping) then
           !FEM representation
           call get_var_from_packed_state(packed_state, FEDensity = den_fem)
           call calculate_u_source( state, den_fem, U_SOURCE_ALL )
        end if


        CALL CV_ASSEMB_FORCE_CTY( state, packed_state, &
             velocity,pressure, &
        NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
        U_ELE_TYPE, P_ELE_TYPE, &
        U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
        U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
        STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
        U_SNLOC, P_SNLOC, CV_SNLOC, &
        X_ALL, U_ABS_STAB_ALL, U_ABSORB_ALL, U_SOURCE_ALL, U_SOURCE_CV_ALL, &
        U_ALL, UOLD_ALL, &
        P_ALL%VAL, CVP_ALL%VAL, DEN_ALL, DENOLD_ALL, DERIV, IDIVID_BY_VOL_FRAC, FEM_VOL_FRAC, &
        DT, &
        NCOLC, FINDC, COLC, & ! C sparcity - global cty eqn
        DGM_PHA, NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, &! Force balance sparcity
        NCOLELE, FINELE, COLELE, & ! Element connectivity.
        NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, & ! pressure matrix for projection method
        NCOLACV, FINACV, COLACV, MIDACV, & ! For CV discretisation method
        SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
        NCOLCT, FINDCT, COLCT, &
        CV_ELE_TYPE, &
        V_DISOPT, V_DG_VEL_INT_OPT, V_THETA, &
        SUF_SIG_DIAGTEN_BC, &
        V_SOURCE, V_ABSORB, VOLFRA_PORE, &
        NCOLM, FINDM, COLM, MIDM, &
        XU_NLOC, XU_NDGLN, &
        U_RHS, MCY_RHS, C, CT, CT_RHS, DIAG_SCALE_PRES, GLOBAL_SOLVE, &
        NLENMCY, NCOLMCY, MCY, FINMCY, PIVIT_MAT, JUST_BL_DIAG_MAT, &
        UDEN_ALL, UDENOLD_ALL, UDIFFUSION_ALL, &
        OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
        IGOT_THETA_FLUX, SCVNGI_THETA, USE_THETA_FLUX, &
        THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, &
        IN_ELE_UPWIND, DG_ELE_UPWIND, &
        NOIT_DIM, RETRIEVE_SOLID_CTY, &
        IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF_ALL, PLIKE_GRAD_SOU_GRAD_ALL,scale_momentum_by_volume_fraction ,&
        StorageIndexes)

        IF ( .NOT.GLOBAL_SOLVE ) THEN
            ! form pres eqn.

            CALL PHA_BLOCK_INV( PIVIT_MAT, TOTELE, U_NLOC * NPHASE * NDIM )

            CALL COLOR_GET_CMC_PHA( CV_NONODS, U_NONODS, NDIM, NPHASE, &
            NCOLC, FINDC, COLC, &
            PIVIT_MAT, &
            TOTELE, U_NLOC, U_NDGLN, &
            NCOLCT, FINDCT, COLCT, DIAG_SCALE_PRES, &
            CMC, CMC_PRECON, IGOT_CMC_PRECON, NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, &
            C, CT, state, StorageIndexes(31) )

        END IF

        DEALLOCATE( ACV )
        NO_MATRIX_STORE = ( NCOLDGM_PHA <= 1 )

        IF ( GLOBAL_SOLVE ) THEN
            ! Global solve
            IF ( JUST_BL_DIAG_MAT ) THEN
                EWRITE(-1,*) 'OPTION NOT READY YET WITH A GLOBAL SOLVE'
                STOP 8331
            END IF
         
            UP = 0.0
            CALL SOLVER( MCY, UP, MCY_RHS, &
            FINMCY, COLMCY, &
            option_path = '/material_phase[0]/vector_field::Velocity')

            U_ALL2 % val = reshape( UP( 1 : U_NONODS * NDIM * NPHASE ), (/ ndim, nphase, u_nonods /) )

            P_ALL % val = UP( U_NONODS * NDIM * NPHASE + 1 : U_NONODS * NDIM * NPHASE + CV_NONODS )

        ELSE ! solve using a projection method

            ! Put pressure in rhs of force balance eqn: CDP = C * P
            CALL C_MULT2( CDP, P_ALL%val , CV_NONODS, U_NONODS, NDIM, NPHASE, C, NCOLC, FINDC, COLC) 

            IF ( JUST_BL_DIAG_MAT .OR. NO_MATRIX_STORE ) THEN

                U_RHS_CDP2 = U_RHS + CDP

                ! DU = BLOCK_MAT * CDP
                CALL PHA_BLOCK_MAT_VEC_old( UP_VEL, PIVIT_MAT, U_RHS_CDP2, U_NONODS, NDIM, NPHASE, &
                TOTELE, U_NLOC, U_NDGLN )

            ELSE

                U_RHS_CDP = RESHAPE( U_RHS + CDP, (/ NDIM * NPHASE * U_NONODS /) )

                UP_VEL = 0.0
                CALL SOLVER( DGM_PHA, UP_VEL, U_RHS_CDP, &
                FINDGM_PHA, COLDGM_PHA, &
                option_path = '/material_phase[0]/vector_field::Velocity', &
                block_size = NDIM*NPHASE*U_NLOC )

            END IF

            U_ALL2 % VAL = RESHAPE( UP_VEL, (/ NDIM, NPHASE, U_NONODS /) )


            !ewrite(3,*) 'u::', u
            !ewrite(3,*) 'v::', v
            !ewrite(3,*) 'w::', w
            !ewrite(3,*) 'ct::', ct
            !ewrite(3,*) 'c::', c
            !ewrite(3,*) 'ct_rhs::', ct_rhs

            ! put on rhs the cty eqn; put most recent pressure in RHS of momentum eqn
            ! NB. P_RHS = -CT * U + CT_RHS
            CALL CT_MULT2( P_RHS, UP_VEL, CV_NONODS, U_NONODS, NDIM, NPHASE, &
            CT, NCOLCT, FINDCT, COLCT )

            P_RHS = -P_RHS + CT_RHS

            ! Matrix vector involving the mass diagonal term
            DO CV_NOD = 1, CV_NONODS
                DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                    CV_JNOD = COLCMC( COUNT )
                    P_RHS( CV_NOD ) = P_RHS( CV_NOD ) &
                    -DIAG_SCALE_PRES( CV_NOD ) * MASS_MN_PRES( COUNT ) * P_ALL%VAL( CV_JNOD )
                END DO
            END DO


           ! Sf => extract_scalar_field( packed_state, "SolidConcentration" )
           !p_rhs = p_rhs * ( 1. - sf%val )
           !p_rhs = p_rhs * ( 0.5 )


            call get_option( '/material_phase[0]/scalar_field::Pressure/' // &
            'prognostic/reference_node', ndpset, default = 0 )
            if ( ndpset /= 0 ) p_rhs( ndpset ) = 0.0

            !ewrite(3,*) 'P_RHS2::', p_rhs
            !ewrite(3,*) 'CT_RHS::', ct_rhs

            ! solve for pressure correction DP that is solve CMC*DP=P_RHS...
            ewrite(3,*)'about to solve for pressure'

            ! Print cmc
            if( .false. ) then
                DO CV_NOD = 1, CV_NONODS
                    ewrite(3,*) 'cv_nod=',cv_nod, &
                    'findcmc=', FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                    rsum = 0.0
                    DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                        CV_JNOD = COLCMC( COUNT )
                        ewrite(3,*) 'count,CV_JNOD,cmc(count):', count, CV_JNOD, cmc( count )
                        if ( cv_nod /= cv_jnod ) rsum = rsum + abs( cmc( count ) )
                    END DO
                    ewrite(3,*) 'off_diag, diag=',rsum,cmc(midcmc(cv_nod))
                END DO
               !stop 1244
            end if

            ewrite(3,*)'b4 pressure solve P_RHS:' !, P_RHS
            DP = 0.

            ! Add diffusion to DG version of CMC to try and encourage a continuous formulation...
            ! the idea is to stabilize pressure without effecting the soln i.e. the rhs of the eqns as
            ! pressure may have some singularities associated with it.
            if ( cv_nonods/=x_nonods .and. .false. ) then !DG only...
                CALL ADD_DIFF_CMC(CMC, &
                NCOLCMC, cv_NONODS, FINDCMC, COLCMC, MIDCMC, &
                totele, cv_nloc, x_nonods, cv_ndgln, x_ndgln, p_all%val )
            end if

            if( cv_nonods == x_nonods .or. .true. ) then ! a continuous pressure

               CALL SOLVER( CMC, DP, P_RHS, &
                    FINDCMC, COLCMC, &
                    option_path = '/material_phase[0]/scalar_field::Pressure' )
            else ! a discontinuous pressure multi-grid solver
               CALL PRES_DG_MULTIGRID(CMC, CMC_PRECON, IGOT_CMC_PRECON, DP, P_RHS, &
                    NCOLCMC, cv_NONODS, FINDCMC, COLCMC, MIDCMC, &
                    totele, cv_nloc, x_nonods, cv_ndgln, x_ndgln )
            end if

            ewrite(3,*) 'after pressure solve DP:', DP

            P_all % val = P_all % val + DP

            ! Use a projection method
            ! CDP = C * DP
            CALL C_MULT2( CDP, DP, CV_NONODS, U_NONODS, NDIM, NPHASE, C, NCOLC, FINDC, COLC )

            ! Correct velocity...
            ! DU = BLOCK_MAT * CDP
            CALL PHA_BLOCK_MAT_VEC2( DU_VEL, PIVIT_MAT, CDP, U_NONODS, NDIM, NPHASE, &
            TOTELE, U_NLOC, U_NDGLN )
            U_ALL2 % VAL = U_ALL2 % VAL + DU_VEL

        END IF

        ! Calculate control volume averaged pressure CV_P from fem pressure P
        CVP_ALL %VAL = 0.0
        MASS_CV = 0.0
        DO CV_NOD = 1, CV_NONODS
            DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                CVP_all%val( CV_NOD ) = CVP_all%val( CV_NOD ) + MASS_MN_PRES( COUNT ) * P_all%val( COLCMC( COUNT ) )
                MASS_CV( CV_NOD ) = MASS_CV( CV_NOD ) + MASS_MN_PRES( COUNT )
            END DO
        END DO
        CVP_all%val = CVP_all%val / MASS_CV

               Pressure_State => extract_scalar_field( state( 1 ), 'Pressure' )
               Pressure_State % val = CVP_all%val



        DEALLOCATE( CT )
        DEALLOCATE( CT_RHS )
        DEALLOCATE( DIAG_SCALE_PRES )
        DEALLOCATE( U_RHS )
        DEALLOCATE( MCY_RHS )
        DEALLOCATE( C )
        DEALLOCATE( MCY )
        DEALLOCATE( CMC )
        DEALLOCATE( MASS_MN_PRES )
        DEALLOCATE( P_RHS )
        DEALLOCATE( UP )
        DEALLOCATE( U_RHS_CDP )
        DEALLOCATE( DP )
        DEALLOCATE( CDP )
        DEALLOCATE( DU_VEL )
        DEALLOCATE( UP_VEL )
        DEALLOCATE( PIVIT_MAT )

        ewrite(3,*) 'Leaving FORCE_BAL_CTY_ASSEM_SOLVE'

    END SUBROUTINE FORCE_BAL_CTY_ASSEM_SOLVE






    ! Add diffusion to CMC to try and encourage a continuous formulation...
    SUBROUTINE ADD_DIFF_CMC(CMC, &
    NCOLCMC, cv_NONODS, FINDCMC, COLCMC, MIDCMC, &
    totele, cv_nloc, x_nonods, cv_ndgln, x_ndgln, p )
        ! Add diffusion to CMC to try and encourage a continuous formulation...
        !
        implicit none
        INTEGER, intent( in ) ::  NCOLCMC, CV_NONODS, totele, cv_nloc, x_nonods
        REAL, DIMENSION( : ), intent( inout ) ::  CMC
        REAL, DIMENSION( : ), intent( inout ) ::  p
        INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
        INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
        INTEGER, DIMENSION( : ), intent( in ) :: MIDCMC
        INTEGER, DIMENSION( : ), intent( in ) :: cv_ndgln, x_ndgln

        ! local variables...

        integer, dimension( : ), allocatable :: dg_nods, MAP_DG2CTY
        real, dimension( : ), allocatable :: diag_lum, P_TEMP
        integer :: ele, cv_iloc, dg_nod, cty_nod, CV_NOD, CV_JNOD
        integer :: count
        real :: alpha

        ! works...
        alpha=1.e-3
        ! can also be used...
        !    alpha=1.e-1
        !    alpha=1.e-2


        allocate( MAP_DG2CTY(cv_nonods) )
        allocate( p_TEMP(X_nonods) )
        allocate( diag_lum(x_nonods) )
        allocate( dg_nods(x_nonods) )

        ! lump the pressure nodes to take away the discontinuity...
        DO ELE = 1, TOTELE
            DO CV_ILOC = 1, CV_NLOC
                !          dg_nod = (ele-1) * cv_nloc + cv_iloc
                dg_nod = cv_ndgln( (ele-1) * cv_nloc + cv_iloc )
                cty_nod = x_ndgln( (ele-1) * cv_nloc + cv_iloc)
                MAP_DG2CTY(dg_nod) = cty_nod
            END DO
        END DO

        diag_lum=0.0
        dg_nods=0
        P_TEMP=0.0
        DO ELE = 1, TOTELE
            DO CV_ILOC = 1, CV_NLOC
                !          dg_nod = (ele-1) * cv_nloc + cv_iloc
                dg_nod = cv_ndgln( (ele-1) * cv_nloc + cv_iloc )
                cty_nod = x_ndgln( (ele-1) * cv_nloc + cv_iloc )
                diag_lum(cty_nod)=diag_lum(cty_nod) + abs( cmc(midcmc(dg_nod)) )
                dg_nods(cty_nod)=dg_nods(cty_nod)+1
                P_TEMP(cty_nod)=P_TEMP(cty_nod)+P(DG_NOD)
            END DO
        END DO
        P_TEMP=p_TEMP/DG_NODS


        DO ELE = 1, TOTELE
            DO CV_ILOC = 1, CV_NLOC
                !          dg_nod = (ele-1) * cv_nloc + cv_iloc
                dg_nod = cv_ndgln( (ele-1) * cv_nloc + cv_iloc )
                cty_nod = x_ndgln( (ele-1) * cv_nloc + cv_iloc )
            ! uncomment to get a cty pressure...
               !P(DG_NOD)=P_TEMP(cty_nod)
            END DO
        END DO

        DO ELE = 1, TOTELE
            DO CV_ILOC = 1, CV_NLOC
                !          dg_nod = (ele-1) * cv_nloc + cv_iloc
                dg_nod = cv_ndgln( (ele-1) * cv_nloc + cv_iloc )
                cty_nod = x_ndgln( (ele-1) * cv_nloc + cv_iloc )
                CV_NOD=DG_NOD
                DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                    CV_JNOD = COLCMC( COUNT )
                    IF(CV_JNOD==CV_NOD) THEN ! on the diagonal...
                        CMC( COUNT ) = CMC( COUNT ) + alpha*diag_lum(cty_nod)
                    ELSE
                        IF(MAP_DG2CTY(CV_JNOD)==cty_nod) THEN ! off diagonal... 
                            CMC( COUNT ) = CMC( COUNT ) - alpha*diag_lum(cty_nod)/real(dg_nods(cty_nod)-1)
                        ENDIF 
                    ENDIF
                END DO

            END DO
        END DO
        RETURN
    END SUBROUTINE ADD_DIFF_CMC






    SUBROUTINE UVW_2_ULONG( U, V, W, UP, U_NONODS, NDIM, NPHASE )
        implicit none
        INTEGER, intent( in ) :: U_NONODS, NDIM, NPHASE
        REAL, DIMENSION( : ), intent( in ) :: U, V, W
        REAL, DIMENSION( : ), intent( inout ) :: UP
        ! Local variables
        INTEGER :: IPHASE

        DO IPHASE = 1, NPHASE
            UP( 1 + ( IPHASE - 1 ) * NDIM * U_NONODS : U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS ) = &
            U( 1 + ( IPHASE - 1 ) * U_NONODS : U_NONODS + ( IPHASE - 1 ) * U_NONODS )
            IF( NDIM >= 2 ) &
            UP( 1 + U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS : 2 * U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS ) = &
            V( 1 + ( IPHASE - 1 ) * U_NONODS : U_NONODS + ( IPHASE - 1 ) * U_NONODS )
            IF( NDIM >= 3 ) &
            UP( 1 + 2 * U_NONODS + ( IPHASE - 1) * NDIM * U_NONODS : 3 * U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS ) = &
            W( 1 + ( IPHASE - 1 ) * U_NONODS : U_NONODS + ( IPHASE - 1 ) * U_NONODS )
        END DO

    END SUBROUTINE UVW_2_ULONG







    SUBROUTINE CV_ASSEMB_FORCE_CTY( state, packed_state, &
         velocity,pressure, &
    NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
    U_ELE_TYPE, P_ELE_TYPE, &
    U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
    U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
    STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
    U_SNLOC, P_SNLOC, CV_SNLOC, &
    X_ALL, U_ABS_STAB_ALL, U_ABSORB_ALL, U_SOURCE_ALL, U_SOURCE_CV_ALL, &
    U_ALL, UOLD_ALL, &
    P, CV_P, DEN_ALL, DENOLD_ALL, DERIV, IDIVID_BY_VOL_FRAC, FEM_VOL_FRAC, &
    DT, &
    NCOLC, FINDC, COLC, & ! C sparcity - global cty eqn
    DGM_PHA, NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, &! Force balance sparcity
    NCOLELE, FINELE, COLELE, & ! Element connectivity.
    NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, & ! pressure matrix for projection method
    NCOLACV, FINACV, COLACV, MIDACV, & ! For CV discretisation method
    SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
    NCOLCT, FINDCT, COLCT, &
    CV_ELE_TYPE, &
    V_DISOPT, V_DG_VEL_INT_OPT, V_THETA, &
    SUF_SIG_DIAGTEN_BC, &
    V_SOURCE, V_ABSORB, VOLFRA_PORE, &
    NCOLM, FINDM, COLM, MIDM, &
    XU_NLOC, XU_NDGLN, &
    U_RHS, MCY_RHS, C, CT, CT_RHS, DIAG_SCALE_PRES, GLOBAL_SOLVE, &
    NLENMCY, NCOLMCY, MCY, FINMCY, PIVIT_MAT, JUST_BL_DIAG_MAT, &
    UDEN_ALL, UDENOLD_ALL, UDIFFUSION_ALL, &
    OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
    IGOT_THETA_FLUX, SCVNGI_THETA, USE_THETA_FLUX, &
    THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, &
    IN_ELE_UPWIND, DG_ELE_UPWIND, &
    NOIT_DIM, RETRIEVE_SOLID_CTY, &
    IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF_ALL, PLIKE_GRAD_SOU_GRAD_ALL ,scale_momentum_by_volume_fraction,&
    StorageIndexes)
        use printout
        implicit none

        ! Form the global CTY and momentum eqns and combine to form one large matrix eqn.

        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type( tensor_field ), intent(in) :: velocity
        type( scalar_field ), intent(in) :: pressure

        INTEGER, intent( in ) :: NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, &
        TOTELE, U_ELE_TYPE, P_ELE_TYPE, &
        U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
        STOTEL, U_SNLOC, P_SNLOC, &
        CV_SNLOC, &
        NCOLC, NCOLDGM_PHA, NCOLELE, NCOLCMC, NCOLACV, NCOLCT, &
        CV_ELE_TYPE, V_DISOPT, V_DG_VEL_INT_OPT, NCOLM, XU_NLOC, &
        NLENMCY, NCOLMCY, NOPT_VEL_UPWIND_COEFS, IGOT_THETA_FLUX, SCVNGI_THETA, &
        IN_ELE_UPWIND, DG_ELE_UPWIND, IPLIKE_GRAD_SOU,  IDIVID_BY_VOL_FRAC
        LOGICAL, intent( in ) :: USE_THETA_FLUX,scale_momentum_by_volume_fraction, RETRIEVE_SOLID_CTY
        INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
        INTEGER, DIMENSION( :  ), intent( in ) :: P_NDGLN
        INTEGER, DIMENSION(  :  ), intent( in ) :: CV_NDGLN
        INTEGER, DIMENSION(  :  ), intent( in ) ::  X_NDGLN
        INTEGER, DIMENSION(  :  ), intent( in ) ::  MAT_NDGLN
        INTEGER, DIMENSION(  :  ), intent( in ) :: CV_SNDGLN
        INTEGER, DIMENSION(  :  ), intent( in ) :: U_SNDGLN
        INTEGER, DIMENSION(  :  ), intent( in ) :: P_SNDGLN
        INTEGER, DIMENSION(  :  ), intent( in ) :: XU_NDGLN
        real, dimension(:,:), intent(in) :: X_ALL
        REAL, DIMENSION(  : ,  : ,  :  ), intent( in ) :: U_ABS_STAB_ALL
        REAL, DIMENSION(  : ,  : ,  :  ), intent( in ) :: U_ABSORB_ALL
        REAL, DIMENSION(  :, :, :  ), intent( in ) :: U_SOURCE_ALL
        REAL, DIMENSION(  :, :, :  ), intent( in ) :: U_SOURCE_CV_ALL
        REAL, DIMENSION(  : ,:,: ), intent( in ) :: U_ALL, UOLD_ALL
        REAL, DIMENSION(  :  ), intent( in ) :: CV_P, P
!        REAL, DIMENSION(  :  ), intent( in ) :: SATURA, SATURAOLD
        REAL, DIMENSION(  :, :  ), intent( in ) :: FEM_VOL_FRAC, DEN_ALL, DENOLD_ALL
        REAL, DIMENSION(  NPHASE, CV_NONODS  ), intent( in ) :: DERIV
        REAL, DIMENSION(  : ,  :   ), intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J
        REAL, intent( in ) :: DT
        INTEGER, DIMENSION(  :  ), intent( in ) :: FINDC
        INTEGER, DIMENSION(  :  ), intent( in ) :: COLC
        REAL, DIMENSION(  :  ), intent( inout ) :: DGM_PHA
        INTEGER, DIMENSION(  :  ), intent( in ) :: FINDGM_PHA
        INTEGER, DIMENSION(  :  ), intent( in ) :: COLDGM_PHA
        INTEGER, DIMENSION(  :  ), intent( in ) :: FINELE
        INTEGER, DIMENSION(  :  ), intent( in ) :: COLELE
        INTEGER, DIMENSION(  :  ), intent( in ) :: FINDCMC
        INTEGER, DIMENSION(  :  ), intent( in ) :: COLCMC
        INTEGER, DIMENSION(  :  ), intent( in ) :: FINACV
        INTEGER, DIMENSION(  :  ), intent( in ) :: COLACV
        INTEGER, DIMENSION(  :  ), intent( in ) :: MIDACV
        integer, dimension(:), intent(in) :: SMALL_FINACV, SMALL_COLACV, small_midacv
        INTEGER, DIMENSION(  :  ), intent( in ) :: FINDCT
        INTEGER, DIMENSION(  :  ), intent( in ) :: COLCT
        REAL, intent( in ) :: V_THETA
        REAL, DIMENSION(  : , : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
        REAL, DIMENSION(  :  ), intent( in ) :: V_SOURCE
        REAL, DIMENSION( :, :, : ), intent( in ) :: V_ABSORB
        REAL, DIMENSION( : ), intent( in ) :: VOLFRA_PORE
        INTEGER, DIMENSION( : ), intent( in ) :: FINDM
        INTEGER, DIMENSION( : ), intent( in ) :: COLM
        INTEGER, DIMENSION( : ), intent( in ) :: MIDM
        REAL, DIMENSION( :, :, : ), intent( inout ) :: U_RHS
        REAL, DIMENSION( : ), intent( inout ) :: MCY_RHS
        REAL, DIMENSION( :, :, : ), intent( inout ) :: C
        REAL, DIMENSION( :, :, : ), intent( inout ) :: CT
        REAL, DIMENSION( : ), intent( inout ) :: MASS_MN_PRES
        REAL, DIMENSION( : ), intent( inout ) :: CT_RHS
        REAL, DIMENSION( : ), intent( inout ) :: DIAG_SCALE_PRES
        LOGICAL, intent( in ) :: GLOBAL_SOLVE
        INTEGER, DIMENSION( : ), intent( in ) :: FINMCY
        REAL, DIMENSION( : ), intent( inout ) :: MCY
        REAL, DIMENSION( :, :,: ), intent( out ) :: PIVIT_MAT
        REAL, DIMENSION( :, : ), intent( in ) :: UDEN_ALL, UDENOLD_ALL
        REAL, DIMENSION( :, :, :, : ), intent( in ) :: UDIFFUSION_ALL
        LOGICAL, intent( inout ) :: JUST_BL_DIAG_MAT
        REAL, DIMENSION( : ), intent( in ) :: OPT_VEL_UPWIND_COEFS
        INTEGER, INTENT( IN ) :: NOIT_DIM
        REAL, DIMENSION( :, :), intent( in ) :: PLIKE_GRAD_SOU_COEF_ALL, PLIKE_GRAD_SOU_GRAD_ALL
        integer, dimension(:), intent(inout) :: StorageIndexes
        ! Local variables
        REAL, PARAMETER :: V_BETA = 1.0
! NEED TO CHANGE RETRIEVE_SOLID_CTY TO MAKE AN OPTION
        REAL :: SECOND_THETA
        LOGICAL, PARAMETER :: GETCV_DISC = .FALSE., GETCT= .TRUE., THERMAL= .FALSE.
        REAL, DIMENSION( : ), allocatable :: ACV, Block_acv, CV_RHS, SUF_VOL_BC_ROB1, SUF_VOL_BC_ROB2, &
        SAT_FEMT, DEN_FEMT, dummy_transp
        REAL, DIMENSION( :,:,:), allocatable :: DENSE_BLOCK_MATRIX
        REAL, DIMENSION( :,:,:,: ), allocatable :: TDIFFUSION
        REAL, DIMENSION( : ), allocatable :: SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, SUF_T2_BC
        INTEGER, DIMENSION( : ), allocatable :: WIC_T2_BC
        REAL, DIMENSION( :, : ), allocatable :: THETA_GDIFF, DEN_OR_ONE, DENOLD_OR_ONE
        REAL, DIMENSION( : ), allocatable :: T2, T2OLD, MEAN_PORE_CV !, DEN_OR_ONE, DENOLD_OR_ONE
        REAL, DIMENSION( :,:,:,: ), allocatable :: THERM_U_DIFFUSION
        LOGICAL :: GET_THETA_FLUX
        INTEGER :: IGOT_T2, I, P_SJLOC, SELE, U_SILOC, IGOT_THERM_VIS

        INTEGER :: U_NLOC2, ILEV, NLEV, ELE, U_ILOC, U_INOD, IPHASE, IDIM, X_ILOC, X_INOD, MAT_INOD, S, E


        type(tensor_field), pointer :: tracer, density

        ewrite(3,*)'In CV_ASSEMB_FORCE_CTY'

        GET_THETA_FLUX = .FALSE.
        IGOT_T2 = 0

!        ALLOCATE( DEN_OR_ONE( NPHASE * CV_NONODS )) ; DEN_OR_ONE = 0.
!        ALLOCATE( DENOLD_OR_ONE( NPHASE * CV_NONODS )) ; DENOLD_OR_ONE = 0.
 
        ALLOCATE( DEN_OR_ONE( NPHASE, CV_NONODS )) ; DEN_OR_ONE = 0.
        ALLOCATE( DENOLD_OR_ONE( NPHASE, CV_NONODS )) ; DENOLD_OR_ONE = 0.
 
 
        ALLOCATE( T2( CV_NONODS * NPHASE * IGOT_T2 )) ; T2 = 0.
        ALLOCATE( T2OLD( CV_NONODS * NPHASE * IGOT_T2 )) ; T2OLD =0.
        ALLOCATE( THETA_GDIFF( NPHASE * IGOT_T2, CV_NONODS * IGOT_T2 )) ; THETA_GDIFF = 0.
        ALLOCATE( ACV( NCOLACV )) ; ACV = 0.
        ALLOCATE( BLOCK_ACV( NPHASE*size(SMALL_COLACV )))  ; BLOCK_ACV = 0.
        ALLOCATE( DENSE_BLOCK_MATRIX( NPHASE,nphase,cv_nonods))  ; DENSE_BLOCK_MATRIX = 0.
        ALLOCATE( CV_RHS( CV_NONODS * NPHASE )) ; CV_RHS = 0.
        ALLOCATE( TDIFFUSION( MAT_NONODS, NDIM, NDIM, NPHASE )) ; TDIFFUSION = 0.
        ALLOCATE( SUF_VOL_BC_ROB1( STOTEL * CV_SNLOC * NPHASE )) ; SUF_VOL_BC_ROB1 = 0.
        ALLOCATE( SUF_VOL_BC_ROB2( STOTEL * CV_SNLOC * NPHASE )) ; SUF_VOL_BC_ROB2 = 0.
        ALLOCATE( MEAN_PORE_CV( CV_NONODS )) ; MEAN_PORE_CV = 0.
        ALLOCATE( SAT_FEMT( NPHASE * CV_NONODS ) ) ; SAT_FEMT = 0.
        ALLOCATE( DEN_FEMT( NPHASE * CV_NONODS ) ) ; DEN_FEMT = 0.
        allocate( dummy_transp( totele ) ) ; dummy_transp = 0.


        TDIFFUSION = 0.0

        IF( GLOBAL_SOLVE ) MCY = 0.0

        ! Obtain the momentum and C matricies
        CALL ASSEMB_FORCE_CTY( state, packed_state, &
             velocity,pressure, &
        NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
        U_ELE_TYPE, P_ELE_TYPE, &
        U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
        U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
        STOTEL, U_SNDGLN, P_SNDGLN, CV_SNDGLN, U_SNLOC, P_SNLOC, CV_SNLOC, &
        X_ALL, U_ABS_STAB_ALL, U_ABSORB_ALL, U_SOURCE_ALL, U_SOURCE_CV_ALL, &
        U_ALL, UOLD_ALL, &
        U_ALL, UOLD_ALL, &    ! This is nu...
        UDEN_ALL, UDENOLD_ALL, IDIVID_BY_VOL_FRAC, FEM_VOL_FRAC, &
        DT, &
        U_RHS, &
        C, NCOLC, FINDC, COLC, & ! C sparsity - global cty eqn
        DGM_PHA, NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, &! Force balance sparsity
        NCOLELE, FINELE, COLELE, & ! Element connectivity.
        XU_NLOC, XU_NDGLN, &
        PIVIT_MAT, JUST_BL_DIAG_MAT, &
        UDIFFUSION_ALL,  DEN_ALL, DENOLD_ALL, RETRIEVE_SOLID_CTY, &
        IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF_ALL, PLIKE_GRAD_SOU_GRAD_ALL, &
        P, NDIM, StorageIndexes=StorageIndexes )
        ! scale the momentum equations by the volume fraction / saturation for the matrix and rhs

        IF ( GLOBAL_SOLVE ) THEN
            ! put momentum and C matrices into global matrix MCY...

            MCY_RHS = 0.0
            DO ELE = 1, TOTELE
                DO U_ILOC = 1, U_NLOC
                    U_INOD = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_ILOC )
                    DO IPHASE = 1, NPHASE
                        DO IDIM = 1, NDIM
                            I = U_INOD + (IDIM-1)*U_NONODS + (IPHASE-1)*NDIM*U_NONODS
                            MCY_RHS( I ) = U_RHS( IDIM, IPHASE, U_INOD )
                        END DO
                    END DO
                END DO
            END DO




            CALL PUT_MOM_C_IN_GLOB_MAT( NPHASE,NDIM, &
            NCOLDGM_PHA, DGM_PHA, FINDGM_PHA, &
            NLENMCY, NCOLMCY, MCY, FINMCY, &
            U_NONODS, NCOLC, C, FINDC )
        END IF

!        IF ( USE_THETA_FLUX ) THEN ! We have already put density in theta...
!            DEN_OR_ONE = 1.0
!            DENOLD_OR_ONE = 1.0
!        ELSE
!            DEN_OR_ONE = DEN
!            DENOLD_OR_ONE = DENOLD
!        END IF



        IF ( USE_THETA_FLUX ) THEN ! We have already put density in theta...
            DEN_OR_ONE = 1.0
            DENOLD_OR_ONE = 1.0
        ELSE
            DEN_OR_ONE = DEN_ALL
            DENOLD_OR_ONE = DENOLD_ALL
        END IF




        ! unused at this stage
        second_theta = 0.0

        IGOT_THERM_VIS=0
        ALLOCATE( THERM_U_DIFFUSION(NDIM,NDIM,NPHASE,MAT_NONODS*IGOT_THERM_VIS ) )


        tracer=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction") 
        density=>extract_tensor_field(packed_state,"PackedDensity")

        call CV_ASSEMB( state, packed_state, &
             tracer, velocity, density, &
        CV_RHS, &
        NCOLACV,  ACV, DENSE_BLOCK_MATRIX, FINACV, COLACV, MIDACV, &
        SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV,&
        NCOLCT, CT, DIAG_SCALE_PRES, CT_RHS, FINDCT, COLCT, &
        CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
        CV_ELE_TYPE,  &
        NPHASE,  &
        CV_NLOC, U_NLOC, X_NLOC, &
        CV_NDGLN, X_NDGLN, U_NDGLN, &
        CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
        DEN_OR_ONE, DENOLD_OR_ONE, &
        MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, IGOT_THERM_VIS, THERM_U_DIFFUSION, &
        V_DISOPT, V_DG_VEL_INT_OPT, DT, V_THETA, SECOND_THETA, V_BETA, &
        SUF_SIG_DIAGTEN_BC, &
        DERIV, CV_P, &
        V_SOURCE, V_ABSORB, VOLFRA_PORE, &
        NDIM, GETCV_DISC, GETCT, &
        NCOLM, FINDM, COLM, MIDM, &
        XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
        OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
        DEN_FEMT, &
        IGOT_T2, T2, T2OLD, IGOT_THETA_FLUX, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
        THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, THETA_GDIFF, &
        IN_ELE_UPWIND, DG_ELE_UPWIND, &
        NOIT_DIM, &
        MEAN_PORE_CV, &
        FINDCMC, COLCMC, NCOLCMC, MASS_MN_PRES, THERMAL,  RETRIEVE_SOLID_CTY,&
        dummy_transp, &
        StorageIndexes, 3 )

        ewrite(3,*)'Back from cv_assemb'

        IF ( GLOBAL_SOLVE ) THEN
            ! Put CT into global matrix MCY...
            MCY_RHS( U_NONODS * NDIM * NPHASE + 1 : U_NONODS * NDIM * NPHASE + CV_NONODS ) = &
            CT_RHS( 1 : CV_NONODS )

            CALL PUT_CT_IN_GLOB_MAT( NPHASE, NDIM, U_NONODS, &
            NLENMCY, NCOLMCY, MCY, FINMCY, &
            CV_NONODS, NCOLCT, CT, DIAG_SCALE_PRES, FINDCT, &
            FINDCMC, NCOLCMC, MASS_MN_PRES )
        END IF

        DEALLOCATE( T2 )
        DEALLOCATE( T2OLD )
        DEALLOCATE( THETA_GDIFF )
        DEALLOCATE( ACV )
        DEALLOCATE( BLOCK_ACV )
        DEALLOCATE( DENSE_BLOCK_MATRIX )
        DEALLOCATE( CV_RHS )
        DEALLOCATE( TDIFFUSION )
        DEALLOCATE( SUF_VOL_BC_ROB1 )
        DEALLOCATE( SUF_VOL_BC_ROB2 )
        DEALLOCATE( MEAN_PORE_CV )
        DEALLOCATE( SAT_FEMT )
        DEALLOCATE( DEN_FEMT )

        ewrite(3,*) 'Leaving CV_ASSEMB_FORCE_CTY'

    END SUBROUTINE CV_ASSEMB_FORCE_CTY




































    SUBROUTINE PUT_MOM_C_IN_GLOB_MAT( NPHASE, NDIM, &
    NCOLDGM_PHA, DGM_PHA, FINDGM_PHA, &
    NLENMCY, NCOLMCY, MCY, FINMCY, &
    U_NONODS, NCOLC, C, FINDC )
        implicit none
        ! put momentum and C matrices into global matrix MCY

        INTEGER, intent( in ) :: NPHASE, NDIM, U_NONODS, NCOLDGM_PHA, &
        NCOLC, NLENMCY, NCOLMCY
        INTEGER, DIMENSION( : ), intent( in ) ::  FINDGM_PHA
        REAL, DIMENSION( : ), intent( in ) ::  DGM_PHA
        INTEGER, DIMENSION( : ), intent( in ) :: FINMCY
        INTEGER, DIMENSION( : ), intent( in ) :: FINDC
        REAL, DIMENSION( : ), intent( inout ) :: MCY
        REAL, DIMENSION( :, :, : ), intent( in ) :: C
        ! Local variables...
        INTEGER :: U_NOD_PHA, IWID, I, U_NOD, IPHASE, IDIM, U_NOD_PHA_I, COUNT, COUNT2

        ewrite(3,*) 'In PUT_MOM_C_IN_GLOB_MAT'

        MCY = 0.0
        ! Put moment matrix DGM_PHA into global matrix MCY
        DO U_NOD_PHA = 1, U_NONODS  * NDIM * NPHASE
            IWID = FINDGM_PHA( U_NOD_PHA + 1 ) - FINDGM_PHA( U_NOD_PHA )

            DO I = 1, IWID
                MCY( FINMCY( U_NOD_PHA ) - 1 + I ) = DGM_PHA( FINDGM_PHA( U_NOD_PHA ) - 1 + I )
            END DO

        END DO

        ! Put C matrix into global matrix MCY

        Loop_IPHASE: DO IPHASE = 1, NPHASE

            Loop_IDIM: DO IDIM = 1, NDIM
                Loop_UNOD: DO U_NOD = 1, U_NONODS

                    U_NOD_PHA_I = U_NOD + ( IDIM - 1 ) * U_NONODS + ( IPHASE - 1 ) * U_NONODS * NDIM
                    IWID = FINDC( U_NOD + 1 ) - FINDC( U_NOD )

                    DO I = 1, IWID
                        COUNT2 = FINMCY( U_NOD_PHA_I + 1 ) - I
                        COUNT = FINDC( U_NOD + 1 ) - I + ( IDIM - 1 ) * NCOLC + ( IPHASE - 1 ) * NCOLC * NDIM
                        MCY( COUNT2 ) = C( IDIM, IPHASE, COUNT )
                    END DO

                END DO Loop_UNOD
            END DO Loop_IDIM
        END DO Loop_IPHASE

        ewrite(3,*) 'Leaving PUT_MOM_C_IN_GLOB_MAT'

    END SUBROUTINE PUT_MOM_C_IN_GLOB_MAT




    SUBROUTINE PUT_CT_IN_GLOB_MAT( NPHASE, NDIM, U_NONODS, &
    NLENMCY, NCOLMCY, MCY, FINMCY, &
    CV_NONODS, NCOLCT, CT, DIAG_SCALE_PRES, FINDCT, &
    FINDCMC, NCOLCMC, MASS_MN_PRES )
        implicit none
        ! Put CT into global matrix MCY

        INTEGER, intent( in ) ::  NPHASE, NDIM, U_NONODS, NLENMCY, NCOLMCY, CV_NONODS, NCOLCT, &
        NCOLCMC
        REAL, DIMENSION( : ), intent( inout ) :: MCY
        INTEGER, DIMENSION( : ), intent( in ) ::  FINMCY
        REAL, DIMENSION( :, :, : ), intent( in ) :: CT
        REAL, DIMENSION( : ), intent( in ) :: DIAG_SCALE_PRES
        INTEGER, DIMENSION( : ), intent( in ) :: FINDCT, FINDCMC
        REAL, DIMENSION( : ), intent( in ) :: MASS_MN_PRES
        ! Local variables...
        INTEGER CV_NOD, IWID, COUNT, IPHASE, COUNT_MCY1, &
        COUNT_MCY, COUNT_CMC, COUNT_TAKE, IDIM, I

        ewrite(3,*) 'In PUT_CT_IN_GLOB_MAT'

        Loop_CVNOD: DO CV_NOD = 1, CV_NONODS
            IWID = FINDCT( CV_NOD + 1 ) - FINDCT( CV_NOD )

            Loop_COUNT: DO COUNT = FINDCT( CV_NOD ), FINDCT( CV_NOD + 1 ) - 1

                Loop_PHASE: DO IPHASE = 1, NPHASE
                    Loop_DIM: DO IDIM = 1, NDIM
                        COUNT_MCY1 = FINMCY( U_NONODS * NPHASE * NDIM + CV_NOD ) - 1 + (COUNT - FINDCT( CV_NOD ) +1) &
                        + ( IPHASE - 1 ) * IWID * NDIM &
                        + IWID*(IDIM-1)
                        MCY( COUNT_MCY1 ) = CT( IDIM, IPHASE, COUNT )

                    END DO Loop_DIM
                END DO Loop_PHASE

            END DO Loop_COUNT

        END DO Loop_CVNOD

        DO CV_NOD = 1, CV_NONODS
            IWID = FINDCMC( CV_NOD + 1 )- FINDCMC( CV_NOD )
            DO I = 1, IWID
                COUNT_CMC = FINDCMC( CV_NOD + 1) - I
                COUNT_MCY = FINMCY( NDIM * NPHASE * U_NONODS + CV_NOD + 1 ) - I
                MCY( COUNT_MCY ) = DIAG_SCALE_PRES( CV_NOD ) * MASS_MN_PRES( COUNT_CMC )
            END DO
        END DO

        ewrite(3,*) 'Leaving PUT_CT_IN_GLOB_MAT'

        RETURN

    END SUBROUTINE PUT_CT_IN_GLOB_MAT





    SUBROUTINE ASSEMB_FORCE_CTY( state, packed_state,&
         velocity,pressure, &
    NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
    U_ELE_TYPE, P_ELE_TYPE, &
    U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
    U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN, &
    STOTEL, U_SNDGLN, P_SNDGLN, CV_SNDGLN, U_SNLOC, P_SNLOC, CV_SNLOC, &
     
    X_ALL, U_ABS_STAB, U_ABSORB, U_SOURCE, U_SOURCE_CV, &
    U_ALL, UOLD_ALL, &
    NU_ALL, NUOLD_ALL, &
    UDEN, UDENOLD,  IDIVID_BY_VOL_FRAC, FEM_VOL_FRAC, &
    DT, &      
    U_RHS, &
    C, NCOLC, FINDC, COLC, & ! C sparsity - global cty eqn
    DGM_PHA, NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, &! Force balance sparsity
    NCOLELE, FINELE, COLELE, & ! Element connectivity.
    XU_NLOC, XU_NDGLN, &
    PIVIT_MAT, JUST_BL_DIAG_MAT,  &
    UDIFFUSION, DEN_ALL, DENOLD_ALL, RETRIEVE_SOLID_CTY, &
    IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, &
         
    P, NDIM_VEL,&
    StorageIndexes )

        implicit none

        type( state_type ), dimension( : ), intent( inout ) :: state
        type( state_type ), intent( inout ) :: packed_state
        type(tensor_field), intent(in) :: velocity
        type(scalar_field), intent(in) :: pressure
! If IGOT_VOL_X_PRESSURE=1 then have a voln fraction in the pressure term and multiply density by volume fraction...
        INTEGER, PARAMETER :: IGOT_VOL_X_PRESSURE = 0
        INTEGER, intent( in ) :: NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
        U_ELE_TYPE, P_ELE_TYPE, U_NONODS, CV_NONODS, X_NONODS, &
        MAT_NONODS, STOTEL, U_SNLOC, P_SNLOC, CV_SNLOC, &
        NCOLC, NCOLDGM_PHA, NCOLELE, XU_NLOC, IPLIKE_GRAD_SOU, NDIM_VEL, IDIVID_BY_VOL_FRAC
! If IDIVID_BY_VOL_FRAC==1 then modify the stress term to take into account dividing through by volume fraction. 
        ! NDIM_VEL
        INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
        INTEGER, DIMENSION( : ), intent( in )  :: P_NDGLN
        INTEGER, DIMENSION(: ), intent( in )  :: CV_NDGLN
        INTEGER, DIMENSION( :), intent( in )  :: X_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: U_SNDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: P_SNDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
        REAL, DIMENSION( :, : ), intent( in ) :: X_ALL

        REAL, DIMENSION( :, :, : ), intent( in ) :: U_ABS_STAB
        REAL, DIMENSION( :, :, : ), intent( in ) :: U_ABSORB
        REAL, DIMENSION( :, :, : ), intent( in ) :: U_SOURCE
        REAL, DIMENSION( :, :, : ), intent( in ) :: U_SOURCE_CV

        REAL, DIMENSION ( :, :, : ), intent( in ) :: U_ALL, UOLD_ALL, NU_ALL, NUOLD_ALL

        REAL, DIMENSION( :, : ), intent( in ) :: UDEN, UDENOLD
        REAL, DIMENSION( NPHASE, CV_NONODS*max(1,IDIVID_BY_VOL_FRAC+IGOT_VOL_X_PRESSURE) ), intent( in ) :: FEM_VOL_FRAC
        REAL, intent( in ) :: DT
        REAL, DIMENSION( :, :, : ), intent( inout ) :: U_RHS
        REAL, DIMENSION( :, :, : ), intent( inout ) :: C
        INTEGER, DIMENSION( : ), intent( in ) :: FINDC
        INTEGER, DIMENSION( : ), intent( in ) :: COLC
        REAL, DIMENSION( : ), intent( inout ) :: DGM_PHA
        INTEGER, DIMENSION( :), intent( in ) :: FINDGM_PHA
        INTEGER, DIMENSION( :), intent( in ) :: COLDGM_PHA
        INTEGER, DIMENSION(: ), intent( in ) :: FINELE
        INTEGER, DIMENSION( : ), intent( in ) :: COLELE
        REAL, DIMENSION( : , : , : ), intent( out ) :: PIVIT_MAT
        REAL, DIMENSION( :, :, :, : ), intent( in ) :: UDIFFUSION
        LOGICAL, intent( inout ) :: JUST_BL_DIAG_MAT
        REAL, DIMENSION( :, : ), intent( in ) :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD
        REAL, DIMENSION( : ), intent( in ) :: P
        REAL, DIMENSION(  :, :  ), intent( in ) :: DEN_ALL, DENOLD_ALL
        LOGICAL, intent( in ) :: RETRIEVE_SOLID_CTY
        integer, dimension(:), intent(inout) :: StorageIndexes
        ! Local Variables
        ! This is for decifering WIC_U_BC & WIC_P_BC
        type( tensor_field ), pointer :: tensorfield
        character( len = option_path_len ) :: option_path
        LOGICAL, PARAMETER :: VOL_ELE_INT_PRES = .TRUE., STRESS_FORM=.true., STAB_VISC_WITH_ABS=.FALSE.
!        LOGICAL, PARAMETER :: POROUS_VEL = .false. ! For reduced variable porous media treatment.
        ! if STAB_VISC_WITH_ABS then stabilize (in the projection mehtod) the viscosity using absorption.
        !      REAL, PARAMETER :: WITH_NONLIN = 1.0, TOLER = 1.E-10, ZERO_OR_TWO_THIRDS=2.0/3.0
        REAL, PARAMETER :: WITH_NONLIN = 1.0, TOLER = 1.E-10, ZERO_OR_TWO_THIRDS=0.0
        !  perform Roe averaging
        LOGICAL, PARAMETER :: ROE_AVE = .false.
        ! NON_LIN_DGFLUX = .TRUE. non-linear DG flux for momentum - if we have an oscillation use upwinding else use central scheme.
        ! UPWIND_DGFLUX=.TRUE. Upwind DG flux.. Else use central scheme. if NON_LIN_DGFLUX = .TRUE. then this option is ignored.
        LOGICAL :: NON_LIN_DGFLUX, UPWIND_DGFLUX
        ! Storage for pointers to the other side of the element.
        ! Switched off for now until this is hooked up.
        LOGICAL, PARAMETER :: STORED_OTHER_SIDE = .FALSE.
        INTEGER, PARAMETER :: ISTORED_OTHER_SIDE = 0
        ! This is for rapid access to the C matrix...
        LOGICAL, PARAMETER :: STORED_AC_SPAR_PT=.FALSE.
        INTEGER, PARAMETER :: IDO_STORE_AC_SPAR_PT=0
        ! re-calculate C matrix...
        LOGICAL :: got_c_matrix

        INTEGER, DIMENSION( :, : ), allocatable ::  FACE_ELE
        INTEGER, DIMENSION( : ), allocatable :: CV_SLOC2LOC, U_SLOC2LOC, &
        U_ILOC_OTHER_SIDE, U_OTHER_LOC, MAT_OTHER_LOC
        REAL, DIMENSION( : ),    ALLOCATABLE ::  &
        SNORMXN, SNORMYN, SNORMZN, SDETWE, NXUDN, VLN,VLN_OLD, &
        XSL,YSL,ZSL, MASS_ELE, VLK
        REAL, DIMENSION( :, : ),    ALLOCATABLE :: XL_ALL, XL2_ALL, XSL_ALL, SNORMXN_ALL, GRAD_SOU_GI_NMX
        REAL, DIMENSION( : ),    ALLOCATABLE :: NORMX_ALL
!        REAL, DIMENSION( :, : ), ALLOCATABLE :: CVN, CVN_SHORT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
!        CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
!        CVFENX_SHORT, CVFENY_SHORT, CVFENZ_SHORT, &
!        UFEN, UFENLX, UFENLY, UFENLZ, SCVFEN, SCVFENSLX, SCVFENSLY, &
!        SCVFENLX, SCVFENLY, SCVFENLZ, &
!        SUFEN, SUFENSLX, SUFENSLY, SUFENLX, SUFENLY, SUFENLZ, &
!        SBCVN, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
!        SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, SBUFENLY, SBUFENLZ
        REAL, DIMENSION( : ), allocatable :: X, Y, Z
        REAL, DIMENSION ( : , :,  : ), allocatable :: SIGMAGI, SIGMAGI_STAB, SIGMAGI_STAB_SOLID_RHS, &
        DUX_ELE, DUY_ELE, DUZ_ELE, DUOLDX_ELE, DUOLDY_ELE, DUOLDZ_ELE, &
        DVX_ELE, DVY_ELE, DVZ_ELE, DVOLDX_ELE, DVOLDY_ELE, DVOLDZ_ELE, &
        DWX_ELE, DWY_ELE, DWZ_ELE, DWOLDX_ELE, DWOLDY_ELE, DWOLDZ_ELE, &
        WORK_ELE_ALL, &
        DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX, FTHETA, SNDOTQ_IN, SNDOTQ_OUT, &
        SNDOTQOLD_IN, SNDOTQOLD_OUT, UD, UDOLD, UD_ND, UDOLD_ND
        REAL, DIMENSION ( : , : ), allocatable :: MAT_M,  &
        DENGI, DENGIOLD,GRAD_SOU_GI, &
        SNDOTQ, SNDOTQOLD, SNDOTQ_ROE, SNDOTQOLD_ROE, SINCOME, SINCOMEOLD, SDEN, SDENOLD, &
        SDEN_KEEP, SDENOLD_KEEP, SDEN2_KEEP, SDENOLD2_KEEP, &
        SNDOTQ_KEEP, SNDOTQ2_KEEP, SNDOTQOLD_KEEP, SNDOTQOLD2_KEEP, &
        N_DOT_DU, N_DOT_DU2, N_DOT_DUOLD, N_DOT_DUOLD2, RHS_U_CV, RHS_U_CV_OLD, UDEN_VFILT, UDENOLD_VFILT
        REAL, DIMENSION ( : , :, : ), allocatable :: SUD_ALL, SUDOLD_ALL, SUD2_ALL, SUDOLD2_ALL, SUD_ALL_KEEP, &
        SUDOLD_ALL_KEEP, SUD2_ALL_KEEP, SUDOLD2_ALL_KEEP
        REAL, DIMENSION ( : ), allocatable :: vel_dot, vel_dot2, velold_dot, velold_dot2, grad_fact

        ! Nonlinear Petrov-Galerkin stuff...
        REAL, DIMENSION ( : , : ), allocatable ::LOC_MASS_INV, LOC_MASS, &
        U_DX, U_DY, U_DZ, V_DX, V_DY, V_DZ, W_DX, W_DY, W_DZ, &
        UOLD_DX, UOLD_DY, UOLD_DZ, VOLD_DX, VOLD_DY, VOLD_DZ, &
        WOLD_DX, WOLD_DY, WOLD_DZ, &
        P_DX

        REAL, DIMENSION ( : ), allocatable :: VLK_UVW, U_R2_COEF, U_GRAD_N_MAX2
        REAL, DIMENSION ( :, :, : ), allocatable :: RESID, &
        MAT_ELE, DIFFGI_U, RHS_DIFF_U, DIFF_VEC_U, SOUGI_X, RESID_U, U_DT, &
        DIF_STAB_U, U_GRAD_NORM2, U_GRAD_NORM, A_DOT_U, STAR_U_COEF, P_STAR_U
        REAL, DIMENSION ( :, :, :, :, : ), allocatable :: UDIFF_SUF_STAB

        LOGICAL, DIMENSION (:,:), ALLOCATABLE :: CV_ON_FACE, CVFEM_ON_FACE,&
                U_ON_FACE,  UFEM_ON_FACE
    !***Pointers for Shape function calculation***
      integer, pointer :: NCOLGPTS
      integer, dimension(:), pointer ::  FINDGPTS, COLGPTS
      integer, dimension(:,:), pointer ::  CV_NEILOC, CV_SLOCLIST, U_SLOCLIST

      real, dimension(:), pointer :: CVWEIGHT, CVWEIGHT_SHORT, SCVFEWEIGH, SBCVFEWEIGH,&
        SELE_OVERLAP_SCALE
      REAL, DIMENSION( : , :, : ), pointer :: CVFENLX_ALL, CVFENLX_SHORT_ALL, UFENLX_ALL,&
      SCVFENLX_ALL, SUFENLX_ALL, SBCVFENLX_ALL, SBUFENLX_ALL
      REAL, DIMENSION( : , : ), pointer :: CVN, CVN_SHORT, CVFEN, CVFEN_SHORT, &
           UFEN, SCVFEN, SCVFENSLX, SCVFENSLY, &
           SUFEN, SUFENSLX, SUFENSLY,  &
           SBCVN,SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
           SBUFEN, SBUFENSLX, SBUFENSLY, &
           DUMMY_ZERO_NDIM_NDIM
    !###Pointers for Shape function calculation###
        real, pointer, dimension(:,:,:) :: CVFENX_ALL, UFENX_ALL
        real, pointer, dimension(:) :: RA, DETWEI
        real, pointer :: VOLUME

! Local variables...
            INTEGER, PARAMETER :: LES_DISOPT=0
! LES_DISOPT is LES option e.g. =0 No LES
!                               =1 Anisotropic element length scale
!                               =2 Take the average length scale h
!                               =3 Take the min length scale h
!                               =4 Take the max length scale h
            REAL, PARAMETER :: LES_THETA=1.0
! LES_THETA =1 is backward Euler for the LES viscocity.
! COEFF_SOLID_FLUID is the coeffficient that determins the magnitude of the relaxation to the solid vel...
            REAL, PARAMETER :: COEFF_SOLID_FLUID = 1.0

        !
        ! Variables used to reduce indirect addressing...
        !INTEGER, DIMENSION ( :, :, : ), allocatable :: WIC_U_BC_ALL
        REAL, DIMENSION ( :, :, : ), allocatable :: LOC_U_RHS
        REAL, DIMENSION ( :, :, :, : ), allocatable :: UFENX_JLOC_U
        REAL, DIMENSION ( :, :, : ), allocatable :: LOC_U, LOC_UOLD, LOC_US, LOC_U_ABS_STAB_SOLID_RHS
        REAL, DIMENSION ( :, :, : ), allocatable :: LOC_NU, LOC_NUOLD
        REAL, DIMENSION ( :, :, : ), allocatable :: LOC_U_ABSORB, LOC_U_ABS_STAB
        REAL, DIMENSION ( :, :, :, : ), allocatable :: LOC_UDIFFUSION, U_DX_ALL, UOLD_DX_ALL, DIFF_FOR_BETWEEN_U
        !REAL, DIMENSION ( :, :, :, : ), allocatable :: SUF_U_BC_ALL, SUF_MOM_BC_ALL, SUF_NU_BC_ALL, SUF_ROB1_UBC_ALL, SUF_ROB2_UBC_ALL, TEN_XX
        REAL, DIMENSION ( :, :, :, : ), allocatable :: TEN_XX

        !REAL, DIMENSION ( :, :, : ), allocatable :: SUF_P_BC_ALL
        REAL, DIMENSION ( :, :, :, :, :, : ), allocatable :: LOC_DGM_PHA
        REAL, DIMENSION ( :, : ), allocatable :: LOC_UDEN,  LOC_UDENOLD
        REAL, DIMENSION ( : ), allocatable :: LOC_P
        REAL, DIMENSION ( :, : ), allocatable :: LOC_PLIKE_GRAD_SOU_COEF, LOC_PLIKE_GRAD_SOU_GRAD
        REAL, DIMENSION ( :, :, : ), allocatable :: LOC_U_SOURCE, LOC_U_SOURCE_CV


        REAL, DIMENSION ( :, :, :,   :, :, :,   : ), allocatable :: DIAG_BIGM_CON, BIGM_CON

        ! memory for fast retreval of surface info...
        INTEGER, DIMENSION ( :, :, : ), allocatable :: STORED_U_ILOC_OTHER_SIDE, STORED_U_OTHER_LOC, STORED_MAT_OTHER_LOC
        INTEGER, DIMENSION ( :, :, : ), allocatable :: POSINMAT_C_STORE
        INTEGER, DIMENSION ( :, :, :, : ), allocatable :: POSINMAT_C_STORE_SUF_DG
        ! To memory access very local...
        REAL, DIMENSION ( :, :, : ), allocatable :: SLOC_U, SLOC_UOLD, SLOC2_U, SLOC2_UOLD
        REAL, DIMENSION ( :, :, : ), allocatable :: SLOC_NU, SLOC_NUOLD, SLOC2_NU, SLOC2_NUOLD
        REAL, DIMENSION ( :, :, :, : ), allocatable :: SLOC_DUX_ELE_ALL, SLOC2_DUX_ELE_ALL, SLOC_DUOLDX_ELE_ALL, SLOC2_DUOLDX_ELE_ALL
        REAL, DIMENSION ( :, : ), allocatable :: SLOC_UDEN, SLOC2_UDEN, SLOC_UDENOLD, SLOC2_UDENOLD
        REAL, DIMENSION ( :, :, :, : ), allocatable :: SLOC_UDIFFUSION, SLOC2_UDIFFUSION
        REAL, DIMENSION ( :, :, : ), allocatable :: SLOC_DIFF_FOR_BETWEEN_U, SLOC2_DIFF_FOR_BETWEEN_U

        REAL, DIMENSION ( :, :, : ), allocatable :: U_NODI_SGI_IPHASE_ALL, U_NODJ_SGI_IPHASE_ALL, UOLD_NODI_SGI_IPHASE_ALL, UOLD_NODJ_SGI_IPHASE_ALL
        ! For derivatives...
        REAL, DIMENSION ( : ), allocatable :: NMX_ALL, VNMX_ALL,  RNMX_ALL

        LOGICAL :: D1, D3, DCYL, GOT_DIFFUS, GOT_UDEN, DISC_PRES, QUAD_OVER_WHOLE_ELE, &
        have_oscillation, have_oscillation_old
        INTEGER :: CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, NFACE
        INTEGER :: IPHASE, ELE, GI, ILOC, GLOBI, GLOBJ, U_NOD, IU_NOD, JCV_NOD, &
        COUNT, COUNT2, IPHA_IDIM, JPHA_JDIM, COUNT_PHA, IU_PHA_NOD, MAT_NOD, SGI, SELE, &
        U_INOD_IDIM_IPHA, U_JNOD_JDIM_IPHA, U_JNOD_JDIM_JPHA, U_SILOC, P_SJLOC, SUF_P_SJ_IPHA, &
        ICV_NOD, IFACE, U_ILOC, U_JLOC, I, J, MAT_ILOC, MAT_NODI, &
        IDIM, P_ILOC, P_JLOC, CV_KLOC, CV_NODK, CV_NODK_PHA, CV_SKLOC, ELE2, ELE3, SELE2, &
        JU_NOD, JU_NOD_PHA, JU_NOD_DIM_PHA, JU_NOD2, JU_NOD2_PHA, JU_NOD2_DIM_PHA, &
        SUF_U_SJ2, SUF_U_SJ2_IPHA, U_ILOC2, U_INOD, U_INOD2, U_JLOC2, U_KLOC, U_NOD_PHA, &
        IU_NOD_PHA, IU_NOD_DIM_PHA, U_NODI_IPHA, U_NODK, U_NODK_PHA, U_SKLOC, X_INOD, X_INOD2, &
        U_NODJ, U_NODJ2, U_NODJ_IPHA, U_SJLOC, X_ILOC, MAT_ILOC2, MAT_INOD, MAT_INOD2, MAT_SILOC, &
        CV_ILOC, CV_JLOC, CV_NOD, CV_NOD_PHA, U_JNOD_IDIM_IPHA, COUNT_PHA2, P_JLOC2, P_JNOD, P_JNOD2, &
        CV_SILOC, JDIM, JPHASE, ILEV, U_NLOC2, CV_KLOC2, CV_NODK2, CV_NODK2_PHA, GI_SHORT, NLEV, STAT, &
        GLOBI_CV, U_INOD_jDIM_jPHA, u_nod2, u_nod2_pha, cv_inod, COUNT_ELE, CV_ILOC2, CV_INOD2
        REAL    :: NN, NXN, NNX, NXNX, NMX, NMY, NMZ, SAREA, &
        VNMX, VNMY, VNMZ, NM, R
        REAL    :: MN, XC, YC, ZC, XC2, YC2, ZC2, HDC, VLM, VLM_NEW,VLM_OLD, NN_SNDOTQ_IN,NN_SNDOTQ_OUT, &
        NN_SNDOTQOLD_IN,NN_SNDOTQOLD_OUT, NORMX, NORMY, NORMZ, RNN, RN, RNMX(3), c1(NDIM), c2(NDIM)
        REAL    :: MASSE, MASSE2, rsum
        ! Nonlinear Petrov-Galerkin stuff...
        INTEGER :: RESID_BASED_STAB_DIF
        REAL :: U_NONLIN_SHOCK_COEF,RNO_P_IN_A_DOT
        REAL :: JTT_INV
        REAL :: VLKNN, U_N

        REAL :: CENT_RELAX,CENT_RELAX_OLD
        INTEGER :: P_INOD, U_INOD_IPHA, U_JNOD, U_KLOC2, U_NODK2, U_NODK2_PHA, GLOBJ_IPHA, IDIM_VEL
        logical firstst,NO_MATRIX_STORE
        character( len = 100 ) :: name

        character( len = option_path_len ) :: overlapping_path
        logical :: mom_conserv, lump_mass, GOT_OTHER_ELE, BETWEEN_ELE_STAB
        real :: beta

        INTEGER :: FILT_DEN, J2, JU2_NOD_DIM_PHA
        LOGICAL :: GOTDEC
        REAL :: NCVM, UFENX_JLOC, UFENY_JLOC, UFENZ_JLOC
        REAL :: FEN_TEN_XX, FEN_TEN_XY,FEN_TEN_XZ
        REAL :: FEN_TEN_YX, FEN_TEN_YY,FEN_TEN_YZ
        REAL :: FEN_TEN_ZX, FEN_TEN_ZY,FEN_TEN_ZZ
        REAL :: MASS_U(U_NLOC,U_NLOC),STORE_MASS_U(U_NLOC,U_NLOC),MASS_U_CV(U_NLOC,CV_NLOC)
        integer :: IPIV(U_NLOC)

        !Variables to improve PIVIT_MAT creation speed
        REAL, DIMENSION ( :, :, :, : ), allocatable :: NN_SIGMAGI_ELE, NN_SIGMAGI_STAB_ELE, NN_SIGMAGI_STAB_SOLID_RHS_ELE, NN_MASS_ELE, NN_MASSOLD_ELE
        REAL, DIMENSION ( :, :, :, :, : ), allocatable :: STRESS_IJ_ELE, DUX_ELE_ALL, DUOLDX_ELE_ALL
        REAL, DIMENSION ( :, :, : ), allocatable :: VLK_ELE
        REAL, DIMENSION ( :, :, :, : ), allocatable :: UDIFFUSION_ALL, LES_UDIFFUSION

! for the option where we divid by voln fraction...
        REAL, DIMENSION ( :, : ), allocatable :: VOL_FRA_GI
        REAL, DIMENSION ( :, :, : ), allocatable :: VOL_FRA_GI_DX_ALL
        REAL, DIMENSION ( :, : ), allocatable :: SLOC_VOL_FRA, SLOC2_VOL_FRA,  SVOL_FRA, SVOL_FRA2, VOL_FRA_NMX_ALL

        logical :: capillary_pressure_activated
        !Variables to store things in state
        type(mesh_type), pointer :: fl_mesh
        type(mesh_type) :: Auxmesh
        type(scalar_field), target :: Targ_C_Mat
        real, dimension(:,:,:), pointer :: Point_C_Mat

!! femdem
        type( vector_field ), pointer :: delta_u_all, us_all
        type(scalar_field), pointer :: sf
        integer :: cv_nodi

        !! Boundary_conditions

        INTEGER, DIMENSION ( ndim , nphase , surface_element_count(velocity) )  :: WIC_U_BC_ALL, WIC_MOMU_BC_ALL, WIC_NU_BC_ALL
        INTEGER, DIMENSION ( surface_element_count(pressure) ) :: WIC_P_BC_ALL
        REAL, DIMENSION ( :, :, : ), pointer  :: SUF_U_BC_ALL, SUF_MOMU_BC_ALL, SUF_NU_BC_ALL
        REAL, DIMENSION ( :, :, : ), pointer :: SUF_U_ROB2_BC_ALL
        REAL, DIMENSION ( : ), pointer :: SUF_P_BC_ALL

        type(scalar_field) :: pressure_BCs
        type(tensor_field) :: velocity_BCs, velocity_BCs_robin2
        type(tensor_field) :: momentum_BCs

        INTEGER, DIMENSION( 4 ), PARAMETER :: ELEMENT_CORNERS=(/1,3,6,10/)


        capillary_pressure_activated = have_option( '/material_phase[0]/multiphase_properties/capillary_pressure' )

        !If we do not have an index where we have stored C, then we need to calculate it
        got_c_matrix  = StorageIndexes(32)/=0
        if (got_c_matrix) then
            !Get from state
            Point_C_Mat(1:size(C,1),1:size(C,2),1:size(C,3)) =>&
            state(1)%scalar_fields(StorageIndexes(32))%ptr%val
            C = Point_C_Mat
        else
            !Prepare stuff to store C in state
            if (has_scalar_field(state(1), "C_MAT")) then
                !If we are recalculating due to a mesh modification then
                !we return to the original situation
                call remove_scalar_field(state(1), "C_MAT")
            end if
            !Get mesh file just to be able to allocate the fields we want to store
            fl_mesh => extract_mesh( state(1), "CoordinateMesh" )
            Auxmesh = fl_mesh
            !The number of nodes I want does not coincide
            Auxmesh%nodes = size(C,1) * size(C,2) * size(C,3)
            call allocate (Targ_C_Mat, Auxmesh)

            !Now we insert them in state and store the index
            call insert(state(1), Targ_C_Mat, "C_MAT")
            StorageIndexes(32) = size(state(1)%scalar_fields)

            !Get from state
            Point_C_Mat(1:size(C,1),1:size(C,2),1:size(C,3)) =>&
            state(1)%scalar_fields(StorageIndexes(32))%ptr%val
            Point_C_Mat = 0.
        end if



        ewrite(3,*) 'In ASSEMB_FORCE_CTY'

        !! get boundary information
        
        call get_entire_boundary_condition(pressure,&
           ['dirichlet'],&
           pressure_BCs,WIC_P_BC_ALL)
        call get_entire_boundary_condition(velocity,&
           ['dirichlet','robin    '],&
           velocity_BCs,WIC_U_BC_ALL,boundary_second_value=velocity_BCs_robin2)
        call get_entire_boundary_condition(velocity,&
           ['momentum     ','momentuminout'],&
           momentum_BCs,WIC_MOMU_BC_ALL)

        WIC_NU_BC_ALL=WIC_U_BC_ALL
        SUF_P_BC_ALL=>pressure_BCs%val
        SUF_U_BC_ALL=>velocity_BCs%val
        SUF_NU_BC_ALL=>velocity_BCs%val
        SUF_MOMU_BC_ALL=>momentum_BCs%val
        SUF_U_ROB2_BC_ALL=>velocity_BCs_robin2%val
        

        !ewrite(3,*) 'Just double-checking sparsity patterns memory allocation:'
        !ewrite(3,*) 'FINDC with size,', size( FINDC ), ':', FINDC( 1 :  size( FINDC ) )
        !ewrite(3,*) 'COLC with size,', size( COLC ), ':', COLC( 1 :  size( COLC ) )
        !ewrite(3,*) 'FINDGM_PHA with size,', size( FINDGM_PHA ), ':', FINDGM_PHA( 1 :  size( FINDGM_PHA ) )
        !ewrite(3,*) 'COLDGM_PHA with size,', size( COLDGM_PHA ), ':', COLDGM_PHA( 1 :  size( COLDGM_PHA ) )
        !ewrite(3,*) 'FINELE with size,', size( FINELE ), ':', FINELE( 1 :  size( FINELE ) )
        !ewrite(3,*) 'COLELE with size,', size( COLELE ), ':', COLELE( 1 :  size( COLELE ) )

        !ewrite(3,*)'UDEN=',uden
        !ewrite(3,*)'UDENOLD=',udenold
        !ewrite(3,*)'u_absorb=',u_absorb
        !ewrite(3,*)'u_abs_stab=',u_abs_stab
        !stop 2921

        mom_conserv=.false.
        call get_option( &
        '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/conservative_advection', &
        beta )
        if (beta>=.999) mom_conserv=.true.
        ewrite(3,*) 'mom_conserv:', mom_conserv

        lump_mass = .false.
        if ( have_option( &
        '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/mass_terms/lump_mass_matrix') &
        ) lump_mass = .true.

        ! This applies a non-linear shock capturing scheme which
        ! may be used to reduce oscillations in velocity or
        ! perform implicit LES modelling of turbulence.
        ! In all residual approaches do not apply Petrov-Galerkin
        ! dissipation on the 1st non-linear iteration within a
        ! time step as there is no good guess of the (U^{n+1}-U^n)/DT.
        ! RESID_BASED_STAB_DIF decides what type of Petrov-Galerkin
        ! method to use.
        ! =1 is the residual squared approach.
        ! =2 is max(0, A . grad U * residual ).
        ! =3 is the max of 1 and 2 (the most dissipative).
        ! U_NONLIN_SHOCK_COEF \in [0,1] is the magnitude of the non-linear
        ! dissipation
        ! =0.25 is small
        ! =1.0 is large
        ! RNO_P_IN_A_DOT \in [0,1] decides if we include the pressure term in
        ! A . grad soln if
        ! =0.0 dont include pressure term.
        ! =1.0 include the pressure term.

        call get_option('/material_phase[0]/vector_field::Velocity/prognostic/' // &
        'spatial_discretisation/discontinuous_galerkin/stabilisation/method', &
        RESID_BASED_STAB_DIF, default=0 )
        BETWEEN_ELE_STAB = RESID_BASED_STAB_DIF/=0 ! Always switch on between element diffusion if using non-linear

        call get_option('/material_phase[0]/vector_field::Velocity/prognostic/' // &
        'spatial_discretisation/discontinuous_galerkin/stabilisation/nonlinear_velocity_coefficient', &
        U_NONLIN_SHOCK_COEF, default=1.)
        call get_option('/material_phase[0]/vector_field::Velocity/prognostic/' // &
        'spatial_discretisation/discontinuous_galerkin/stabilisation/include_pressure', &
        RNO_P_IN_A_DOT, default=1.)

        ewrite(3,*) 'RESID_BASED_STAB_DIF, U_NONLIN_SHOCK_COEF, RNO_P_IN_A_DOT:', &
        RESID_BASED_STAB_DIF, U_NONLIN_SHOCK_COEF, RNO_P_IN_A_DOT

        QUAD_OVER_WHOLE_ELE = is_overlapping.OR.is_compact_overlapping ! Do NOT divide element into CV's to form quadrature.
        call retrieve_ngi( ndim, u_ele_type, cv_nloc, u_nloc, &
        cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE )
        if ( is_overlapping ) then
            nlev = cv_nloc
            U_NLOC2 = max( 1, U_NLOC/CV_NLOC )
        else
            nlev = 1
            U_NLOC2 = U_NLOC
        end if

        GOT_DIFFUS = .FALSE.

        ! is this the 1st iteration of the time step.
        FIRSTST = ( SUM( (U_ALL(1,:,:) - UOLD_ALL(1,:,:) ) **2) < 1.e-10 )
        IF(NDIM_VEL>=2) FIRSTST = FIRSTST .OR. ( SUM( ( U_ALL(2,:,:) - UOLD_ALL(2,:,:) )**2 ) < 1.e-10 )
        IF(NDIM_VEL>=3) FIRSTST = FIRSTST .OR. ( SUM( ( U_ALL(3,:,:) - UOLD_ALL(3,:,:) )**2 ) < 1.e-10 )

        UPWIND_DGFLUX = .TRUE.
        if ( have_option( &
        '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/advection_scheme/central_differencing') &
        ) UPWIND_DGFLUX = .FALSE.
        NON_LIN_DGFLUX = .FALSE.
        if ( have_option( &
        '/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/advection_scheme/nonlinear_flux') &
        ) NON_LIN_DGFLUX = .TRUE.

        ALLOCATE( UD( NDIM_VEL, NPHASE, CV_NGI ))
        ALLOCATE( UDOLD( NDIM_VEL, NPHASE, CV_NGI ))

        ALLOCATE( UD_ND( NDIM, NPHASE, CV_NGI ))
        ALLOCATE( UDOLD_ND( NDIM, NPHASE, CV_NGI ))

        ALLOCATE( DENGI( NPHASE, CV_NGI ))
        ALLOCATE( DENGIOLD( NPHASE, CV_NGI ))
        ALLOCATE( GRAD_SOU_GI( NPHASE, CV_NGI ))

        ALLOCATE( RHS_U_CV( NPHASE, U_NLOC ))
        ALLOCATE( RHS_U_CV_OLD( NPHASE, U_NLOC ))
        ALLOCATE( UDEN_VFILT( NPHASE, U_NLOC ))
        ALLOCATE( UDENOLD_VFILT( NPHASE, U_NLOC ))

        ALLOCATE( SIGMAGI( NDIM_VEL * NPHASE, NDIM_VEL * NPHASE, CV_NGI ))
        ALLOCATE( SIGMAGI_STAB( NDIM_VEL * NPHASE, NDIM_VEL * NPHASE, CV_NGI ))
        ALLOCATE( MAT_M( MAT_NLOC, CV_NGI ))
        ALLOCATE( SNORMXN( SBCVNGI ))
        ALLOCATE( SNORMYN( SBCVNGI ))
        ALLOCATE( SNORMZN( SBCVNGI ))

        ALLOCATE( XL_ALL(NDIM,CV_NLOC), XL2_ALL(NDIM,CV_NLOC), XSL_ALL(NDIM,CV_SNLOC) )
        ALLOCATE( NORMX_ALL(NDIM), SNORMXN_ALL(NDIM,SBCVNGI) )

        !Variables to improve PIVIT_MAT creation speed
        ALLOCATE(NN_SIGMAGI_ELE( NDIM_VEL * NPHASE, U_NLOC, NDIM_VEL * NPHASE,U_NLOC ))
        ALLOCATE(NN_SIGMAGI_STAB_ELE( NDIM_VEL * NPHASE, U_NLOC, NDIM_VEL * NPHASE,U_NLOC ))
        ALLOCATE(NN_MASS_ELE( NDIM_VEL * NPHASE, U_NLOC, NDIM_VEL * NPHASE,U_NLOC ))
        ALLOCATE(NN_MASSOLD_ELE( NDIM_VEL * NPHASE, U_NLOC, NDIM_VEL * NPHASE,U_NLOC ))
        ALLOCATE( STRESS_IJ_ELE( NDIM, NDIM, NPHASE, U_NLOC,U_NLOC ))
        ALLOCATE( VLK_ELE( NPHASE, U_NLOC, U_NLOC ))


        IF(RETRIEVE_SOLID_CTY) THEN
           ALLOCATE( SIGMAGI_STAB_SOLID_RHS( NDIM_VEL * NPHASE, NDIM_VEL * NPHASE, CV_NGI ))
           ALLOCATE(NN_SIGMAGI_STAB_SOLID_RHS_ELE( NDIM_VEL * NPHASE, U_NLOC, NDIM_VEL * NPHASE,U_NLOC ))
        ENDIF

        ALLOCATE( NXUDN( SCVNGI ))

        ALLOCATE( SDETWE( SBCVNGI ))

        ALLOCATE( CV_SLOC2LOC( CV_SNLOC ))
        ALLOCATE( U_SLOC2LOC( U_SNLOC ))
        ALLOCATE( U_ILOC_OTHER_SIDE(U_SNLOC))

        ALLOCATE( CV_ON_FACE( CV_NLOC, SCVNGI ))
        ALLOCATE( CVFEM_ON_FACE( CV_NLOC, SCVNGI ))
        ALLOCATE( U_ON_FACE( U_NLOC, SCVNGI ))
        ALLOCATE( UFEM_ON_FACE( U_NLOC, SCVNGI ))
        ALLOCATE( U_OTHER_LOC( U_NLOC ))
        ALLOCATE( MAT_OTHER_LOC( MAT_NLOC ))

        ALLOCATE( TEN_XX( NDIM, NDIM, NPHASE, CV_NGI ))

        ALLOCATE( VLN( NPHASE ))
        ALLOCATE( VLN_OLD( NPHASE ))
        ALLOCATE( VLK( NPHASE ))

        ALLOCATE( SUD_ALL(NDIM,NPHASE,SBCVNGI) )
        ALLOCATE( SUDOLD_ALL(NDIM,NPHASE,SBCVNGI) )
        ALLOCATE( SUD2_ALL(NDIM,NPHASE,SBCVNGI) )
        ALLOCATE( SUDOLD2_ALL(NDIM,NPHASE,SBCVNGI) )

        ALLOCATE( SNDOTQ(NPHASE,SBCVNGI) )
        ALLOCATE( SNDOTQOLD(NPHASE,SBCVNGI) )
        ALLOCATE( SNDOTQ_ROE(NPHASE,SBCVNGI) )
        ALLOCATE( SNDOTQOLD_ROE(NPHASE,SBCVNGI) )
        ALLOCATE( SINCOME(NPHASE,SBCVNGI) )
        ALLOCATE( SINCOMEOLD(NPHASE,SBCVNGI) )
        ALLOCATE( SDEN(NPHASE,SBCVNGI) )
        ALLOCATE( SDENOLD(NPHASE,SBCVNGI) )

        ALLOCATE( SDEN_KEEP(NPHASE,SBCVNGI) )
        ALLOCATE( SDENOLD_KEEP(NPHASE,SBCVNGI) )
        ALLOCATE( SDEN2_KEEP(NPHASE,SBCVNGI) )
        ALLOCATE( SDENOLD2_KEEP(NPHASE,SBCVNGI) )

        ALLOCATE( SUD_ALL_KEEP(NDIM,NPHASE,SBCVNGI) )
        ALLOCATE( SUDOLD_ALL_KEEP(NDIM,NPHASE,SBCVNGI) )
        ALLOCATE( SUD2_ALL_KEEP(NDIM,NPHASE,SBCVNGI) )
        ALLOCATE( SUDOLD2_ALL_KEEP(NDIM,NPHASE,SBCVNGI) )

        ALLOCATE( SNDOTQ_KEEP(NPHASE,SBCVNGI) )
        ALLOCATE( SNDOTQ2_KEEP(NPHASE,SBCVNGI) )
        ALLOCATE( SNDOTQOLD_KEEP(NPHASE,SBCVNGI) )
        ALLOCATE( SNDOTQOLD2_KEEP(NPHASE,SBCVNGI) )

        ALLOCATE( N_DOT_DU(NPHASE,SBCVNGI) )
        ALLOCATE( N_DOT_DU2(NPHASE,SBCVNGI) )
        ALLOCATE( N_DOT_DUOLD(NPHASE,SBCVNGI) )
        ALLOCATE( N_DOT_DUOLD2(NPHASE,SBCVNGI) )

        ALLOCATE( vel_dot(SBCVNGI), vel_dot2(SBCVNGI), velold_dot(SBCVNGI), velold_dot2(SBCVNGI), grad_fact(SBCVNGI) )

        ALLOCATE( DIFF_COEF_DIVDX(NDIM_VEL,NPHASE,SBCVNGI) )
        ALLOCATE( DIFF_COEFOLD_DIVDX(NDIM_VEL,NPHASE,SBCVNGI) )
        ALLOCATE( FTHETA(NDIM_VEL,NPHASE,SBCVNGI) )
        ALLOCATE( SNDOTQ_IN(NDIM_VEL,NPHASE,SBCVNGI) )
        ALLOCATE( SNDOTQ_OUT(NDIM_VEL,NPHASE,SBCVNGI) )
        ALLOCATE( SNDOTQOLD_IN(NDIM_VEL,NPHASE,SBCVNGI) )
        ALLOCATE( SNDOTQOLD_OUT(NDIM_VEL,NPHASE,SBCVNGI) )

        ALLOCATE( XSL(CV_SNLOC) )
        ALLOCATE( YSL(CV_SNLOC) )
        ALLOCATE( ZSL(CV_SNLOC) )

!        ALLOCATE( SELE_OVERLAP_SCALE( CV_NLOC ) )


        ALLOCATE( GRAD_SOU_GI_NMX( NDIM_VEL, NPHASE ))


        ALLOCATE( MASS_ELE( TOTELE ))
        MASS_ELE=0.0

        ! Allocating for non-linear Petrov-Galerkin diffusion stabilization...
        ALLOCATE( LOC_MASS_INV(U_NLOC, U_NLOC) )
        ALLOCATE( LOC_MASS(U_NLOC, U_NLOC) )
        ALLOCATE( RHS_DIFF_U(NDIM_VEL,NPHASE,U_NLOC) )

        ALLOCATE( DIFF_VEC_U(NDIM_VEL,NPHASE,U_NLOC) )

        ALLOCATE( DIFFGI_U(NDIM_VEL,NPHASE,CV_NGI) )

        ALLOCATE( U_DT(NDIM_VEL, NPHASE,CV_NGI) )


        ALLOCATE( U_DX_ALL( NDIM, NDIM_VEL, NPHASE, CV_NGI ) )
        ALLOCATE( UOLD_DX_ALL( NDIM, NDIM_VEL, NPHASE, CV_NGI ) )

        ALLOCATE( UOLD_DX(CV_NGI,NPHASE), UOLD_DY(CV_NGI,NPHASE), UOLD_DZ(CV_NGI,NPHASE) )
        ALLOCATE( VOLD_DX(CV_NGI,NPHASE), VOLD_DY(CV_NGI,NPHASE), VOLD_DZ(CV_NGI,NPHASE) )
        ALLOCATE( WOLD_DX(CV_NGI,NPHASE), WOLD_DY(CV_NGI,NPHASE), WOLD_DZ(CV_NGI,NPHASE) )

        ALLOCATE( SOUGI_X(NDIM_VEL,NPHASE,CV_NGI) )


        ALLOCATE( RESID_U(NDIM_VEL,NPHASE,CV_NGI) )
        ALLOCATE( P_DX(NDIM, CV_NGI)  )

        ALLOCATE( U_GRAD_NORM2(NDIM_VEL,NPHASE,CV_NGI), U_GRAD_NORM(NDIM_VEL,NPHASE,CV_NGI) )


        ALLOCATE( A_DOT_U(NDIM_VEL,NPHASE,CV_NGI) )
        ALLOCATE( STAR_U_COEF(NDIM_VEL,NPHASE,CV_NGI) )
        ALLOCATE( P_STAR_U(NDIM_VEL,NPHASE,CV_NGI) )
        ALLOCATE( DIF_STAB_U(NDIM_VEL,NPHASE, CV_NGI) )


        ALLOCATE( U_R2_COEF( NDIM_VEL ) )
        ALLOCATE( U_GRAD_N_MAX2( NDIM_VEL ) )
        ALLOCATE( VLK_UVW(NDIM_VEL) )

        ! Variables used to reduce indirect addressing...
        ALLOCATE( LOC_U(NDIM_VEL, NPHASE, U_NLOC),  LOC_UOLD(NDIM_VEL, NPHASE, U_NLOC) )
        IF(RETRIEVE_SOLID_CTY) THEN
           ALLOCATE( LOC_US(NDIM_VEL, NPHASE, U_NLOC))
           ALLOCATE( LOC_U_ABS_STAB_SOLID_RHS(NDIM_VEL* NPHASE, NDIM_VEL* NPHASE, MAT_NLOC))
        ENDIF
        ALLOCATE( LOC_NU(NDIM, NPHASE, U_NLOC),  LOC_NUOLD(NDIM, NPHASE, U_NLOC) )
        ALLOCATE( LOC_UDEN(NPHASE, CV_NLOC),  LOC_UDENOLD(NPHASE, CV_NLOC) )
        ALLOCATE( LOC_P(P_NLOC) )
        ALLOCATE( LOC_PLIKE_GRAD_SOU_COEF(NPHASE, CV_NLOC) )
        ALLOCATE( LOC_PLIKE_GRAD_SOU_GRAD(NPHASE, CV_NLOC) )
        ALLOCATE( LOC_U_SOURCE(NDIM_VEL, NPHASE, U_NLOC) )
        ALLOCATE( LOC_U_SOURCE_CV(NDIM_VEL, NPHASE, CV_NLOC) )
        ALLOCATE( LOC_U_ABSORB  (NDIM_VEL* NPHASE, NDIM_VEL* NPHASE, MAT_NLOC) )
        ALLOCATE( LOC_U_ABS_STAB(NDIM_VEL* NPHASE, NDIM_VEL* NPHASE, MAT_NLOC) )
        ALLOCATE( LOC_UDIFFUSION(NDIM, NDIM, NPHASE, MAT_NLOC) )
        ALLOCATE( LOC_U_RHS( NDIM_VEL, NPHASE, U_NLOC ) )
        ALLOCATE( UFENX_JLOC_U(NDIM,NDIM,CV_NGI,U_NLOC) )

        ! To memory access very local...
        ALLOCATE( SLOC_U(NDIM_VEL,NPHASE,U_SNLOC) )
        ALLOCATE( SLOC_UOLD(NDIM_VEL,NPHASE,U_SNLOC) )
        ALLOCATE( SLOC2_U(NDIM_VEL,NPHASE,U_SNLOC) )
        ALLOCATE( SLOC2_UOLD(NDIM_VEL,NPHASE,U_SNLOC) )

        ALLOCATE( SLOC_NU(NDIM,NPHASE,U_SNLOC) )
        ALLOCATE( SLOC_NUOLD(NDIM,NPHASE,U_SNLOC) )
        ALLOCATE( SLOC2_NU(NDIM,NPHASE,U_SNLOC) )
        ALLOCATE( SLOC2_NUOLD(NDIM,NPHASE,U_SNLOC) )

        ALLOCATE( SLOC_DUX_ELE_ALL( NDIM_VEL, NDIM , NPHASE, U_SNLOC ) )
        ALLOCATE( SLOC2_DUX_ELE_ALL( NDIM_VEL, NDIM , NPHASE, U_SNLOC ) )
        ALLOCATE( SLOC_DUOLDX_ELE_ALL( NDIM_VEL, NDIM , NPHASE, U_SNLOC ) )
        ALLOCATE( SLOC2_DUOLDX_ELE_ALL( NDIM_VEL, NDIM , NPHASE, U_SNLOC ) )

        ALLOCATE( SLOC_DIFF_FOR_BETWEEN_U(NDIM_VEL,NPHASE,U_SNLOC) )
        ALLOCATE( SLOC2_DIFF_FOR_BETWEEN_U(NDIM_VEL,NPHASE,U_SNLOC) )

        ALLOCATE( SLOC_UDEN(NPHASE, CV_SNLOC)  )
        ALLOCATE( SLOC2_UDEN(NPHASE, CV_SNLOC)  )
        ALLOCATE( SLOC_UDENOLD(NPHASE, CV_SNLOC)  )
        ALLOCATE( SLOC2_UDENOLD(NPHASE, CV_SNLOC) )

        ALLOCATE( SLOC_UDIFFUSION(NDIM, NDIM, NPHASE, CV_SNLOC) )
        ALLOCATE( SLOC2_UDIFFUSION(NDIM, NDIM, NPHASE, CV_SNLOC) )

        ! Derivatives...
        ALLOCATE( NMX_ALL(NDIM) )
        ALLOCATE( VNMX_ALL(NDIM) )
        ALLOCATE( RNMX_ALL(NDIM) )

        ALLOCATE( U_NODI_SGI_IPHASE_ALL(NDIM_VEL,NPHASE,SBCVNGI) )
        ALLOCATE( U_NODJ_SGI_IPHASE_ALL(NDIM_VEL,NPHASE,SBCVNGI) )
        ALLOCATE( UOLD_NODI_SGI_IPHASE_ALL(NDIM_VEL,NPHASE,SBCVNGI) )
        ALLOCATE( UOLD_NODJ_SGI_IPHASE_ALL(NDIM_VEL,NPHASE,SBCVNGI) )

        ALLOCATE( X(X_NONODS), Y(X_NONODS), Z(X_NONODS) ) !; X=0. ; Y=0. ; Z=0.

        IF(IDIVID_BY_VOL_FRAC+IGOT_VOL_X_PRESSURE.GE.1) THEN
            ALLOCATE( VOL_FRA_GI(NPHASE, CV_NGI_SHORT), VOL_FRA_GI_DX_ALL( NDIM, NPHASE, CV_NGI_SHORT) )
            ALLOCATE( SLOC_VOL_FRA(NPHASE, CV_SNLOC), SLOC2_VOL_FRA(NPHASE, CV_SNLOC))
            ALLOCATE( SVOL_FRA(NPHASE, SBCVNGI), SVOL_FRA2(NPHASE, SBCVNGI) )
            ALLOCATE( VOL_FRA_NMX_ALL(NDIM,NPHASE) )
        ENDIF

        DO IDIM = 1, NDIM
            IF ( IDIM == 1 ) THEN
                X = X_ALL( IDIM, : )
            ELSE IF ( IDIM == 2 ) THEN
                Y = X_ALL( IDIM, : )
            ELSE
                Z = X_ALL( IDIM, : )
            END IF
        END DO

        GOT_DIFFUS = ( R2NORM( UDIFFUSION, MAT_NONODS * NDIM * NDIM * NPHASE ) /= 0.0 )  &
        .OR. BETWEEN_ELE_STAB
        IF(LES_DISOPT.NE.0) GOT_DIFFUS=.TRUE.

        GOT_UDEN = .FALSE.
        DO IPHASE = 1, NPHASE
            GOT_UDEN = GOT_UDEN .OR. ( R2NORM( UDEN( IPHASE, : ), CV_NONODS ) /= 0.0 )
        END DO

        JUST_BL_DIAG_MAT=( ( .NOT. GOT_DIFFUS ) .AND. ( .NOT. GOT_UDEN ) )

        ALLOCATE( UDIFF_SUF_STAB( NDIM_VEL, NDIM, NDIM, NPHASE, SBCVNGI ) )
        UDIFF_SUF_STAB = 0.0

        IF ( BETWEEN_ELE_STAB ) THEN
            ! Calculate stabilization diffusion coefficient between elements...
            ALLOCATE( DIFF_FOR_BETWEEN_U( NDIM_VEL, NPHASE, U_NLOC, TOTELE ) ) ; DIFF_FOR_BETWEEN_U = 0.0
            ALLOCATE( MAT_ELE( U_NLOC, U_NLOC, TOTELE ) ) ; MAT_ELE = 0.0
        END IF

        IF ( GOT_DIFFUS ) THEN
            ALLOCATE( DUX_ELE_ALL( NDIM_VEL, NDIM, NPHASE, U_NLOC, TOTELE ) )
            ALLOCATE( DUOLDX_ELE_ALL( NDIM_VEL, NDIM, NPHASE, U_NLOC, TOTELE ) )
            ALLOCATE( WORK_ELE_ALL( U_NLOC, NPHASE, TOTELE ) )
        ENDIF

        D1   = ( NDIM == 1 )
        DCYL = ( NDIM ==-2 )
        D3   = ( NDIM == 3 )

        NO_MATRIX_STORE = NCOLDGM_PHA<=1
        IF( (.NOT.JUST_BL_DIAG_MAT) .AND. (.NOT.NO_MATRIX_STORE) ) DGM_PHA = 0.0
        if (.not.got_c_matrix) C = 0.0
        U_RHS = 0.0
        IF (.NOT.NO_MATRIX_STORE ) THEN
            ALLOCATE( DIAG_BIGM_CON( NDIM_VEL, NDIM_VEL, NPHASE, NPHASE, U_NLOC, U_NLOC, TOTELE ) )
            ALLOCATE( BIGM_CON( NDIM_VEL, NDIM_VEL, NPHASE, NPHASE, U_NLOC, U_NLOC, NCOLELE ) )
            DIAG_BIGM_CON = 0.0
            BIGM_CON = 0.0
        END IF

        !======= DEFINE THE SUB-CONTROL VOLUME SHAPE FUNCTIONS, ETC ========

        ! Shape functions associated with volume integration using both CV basis
        ! functions CVN as well as FEM basis functions CVFEN (and its derivatives CVFENLX, CVFENLY, CVFENLZ)

        !======= DEFINE THE SUB-CONTROL VOLUME & FEM SHAPE FUNCTIONS ========
        CALL cv_fem_shape_funs_plus_storage( &
                             ! Volume shape functions...
           NDIM, P_ELE_TYPE,  &
           CV_NGI, CV_NGI_SHORT, CV_NLOC, U_NLOC, CVN, CVN_SHORT, &
           CVWEIGHT, CVFEN, CVFENLX_ALL, &
           CVWEIGHT_SHORT, CVFEN_SHORT, CVFENLX_SHORT_ALL, &
           UFEN, UFENLX_ALL, &
                                ! Surface of each CV shape functions...
           SCVNGI, CV_NEILOC, CV_ON_FACE, CVFEM_ON_FACE, &
           SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
           SCVFENLX_ALL,  &
           SUFEN, SUFENSLX, SUFENSLY,  &
           SUFENLX_ALL,  &
                                ! Surface element shape funcs...
           U_ON_FACE, UFEM_ON_FACE, NFACE, &
           SBCVNGI, SBCVN, SBCVFEN,SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SBCVFENLX_ALL, &
           SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX_ALL, &
           CV_SLOCLIST, U_SLOCLIST, CV_SNLOC, U_SNLOC, &
                             ! Define the gauss points that lie on the surface of the CV...
        FINDGPTS, COLGPTS, NCOLGPTS, &
        SELE_OVERLAP_SCALE, QUAD_OVER_WHOLE_ELE,&
        state, 'ASEMCTY', StorageIndexes(41))


        ! Memory for rapid retreval...
        ! Storage for pointers to the other side of the element.
        ALLOCATE( STORED_U_ILOC_OTHER_SIDE( U_SNLOC, NFACE, TOTELE*ISTORED_OTHER_SIDE ) )
        ALLOCATE( STORED_U_OTHER_LOC( U_NLOC, NFACE, TOTELE*ISTORED_OTHER_SIDE ) )
        ALLOCATE( STORED_MAT_OTHER_LOC( MAT_NLOC, NFACE, TOTELE*ISTORED_OTHER_SIDE ) )

        ALLOCATE( POSINMAT_C_STORE( U_NLOC,P_NLOC, TOTELE*IDO_STORE_AC_SPAR_PT) )
        ALLOCATE( POSINMAT_C_STORE_SUF_DG( U_SNLOC,P_SNLOC,NFACE,TOTELE*IDO_STORE_AC_SPAR_PT ) )

        ALLOCATE( FACE_ELE( NFACE, TOTELE ) )
        ! Calculate FACE_ELE
        CALL CALC_FACE_ELE( FACE_ELE, TOTELE, STOTEL, NFACE, &
        NCOLELE, FINELE, COLELE, CV_NLOC, CV_SNLOC, CV_NONODS, CV_NDGLN, CV_SNDGLN, &
        CV_SLOCLIST, X_NLOC, X_NDGLN )

        IF( GOT_DIFFUS ) THEN
            CALL DG_DERIVS_ALL( U_ALL, UOLD_ALL, &
            DUX_ELE_ALL, DUOLDX_ELE_ALL, &
            NDIM, NPHASE, NDIM_VEL, U_NONODS, TOTELE, U_NDGLN, &
            XU_NDGLN, X_NLOC, X_NDGLN, &
            CV_NGI, U_NLOC, CVWEIGHT, &
            UFEN, UFENLX_ALL(1,:,:), UFENLX_ALL(2,:,:), UFENLX_ALL(3,:,:), &
            CVFEN, CVFENLX_ALL(1,:,:), CVFENLX_ALL(2,:,:), CVFENLX_ALL(3,:,:), &
            X_NONODS, X, Y, Z, &
            NFACE, FACE_ELE, U_SLOCLIST, CV_SLOCLIST, STOTEL, U_SNLOC, CV_SNLOC, WIC_U_BC_ALL, SUF_U_BC_ALL, &
            SBCVNGI, SBUFEN, SBUFENSLX, SBUFENSLY, SBCVFEWEIGH, &
            SBCVFEN, SBCVFENSLX, SBCVFENSLY ,&
            state, "CTY", StorageIndexes(23))
        ENDIF


! LES VISCOCITY CALC.
        IF(GOT_DIFFUS ) THEN
           ALLOCATE(UDIFFUSION_ALL(NDIM,NDIM,NPHASE,MAT_NONODS))
           IF(LES_DISOPT.NE.0) THEN
              ALLOCATE(LES_UDIFFUSION(NDIM,NDIM,NPHASE,MAT_NONODS))
              CALL VISCOCITY_TENSOR_LES_CALC(LES_UDIFFUSION, LES_THETA*DUX_ELE_ALL + (1.-LES_THETA)*DUOLDX_ELE_ALL, &
                                             NDIM,NPHASE, U_NLOC,X_NLOC,TOTELE, X_NONODS, &
                                             X_ALL, X_NDGLN,  MAT_NONODS, MAT_NLOC, MAT_NDGLN, LES_DISOPT)
              UDIFFUSION_ALL=UDIFFUSION + LES_UDIFFUSION
           ELSE
              UDIFFUSION_ALL=UDIFFUSION
           ENDIF
        ENDIF

        if( RETRIEVE_SOLID_CTY ) THEN
           sf=> extract_scalar_field( packed_state, "SolidConcentration" )
           delta_u_all => extract_vector_field( packed_state, "delta_U" ) ! this is delta_u
           us_all => extract_vector_field( packed_state, "solid_U" )
        endif


        !This term is obtained from the surface tension and curvature
        !For capillary pressure we are using the entry pressure method instead of
        !calculating the entry pressure from the surface tension and curvature
        IF ( capillary_pressure_activated ) GRAD_SOU_GI = 1.0

        Loop_Elements: DO ELE = 1, TOTELE ! Volume integral

            ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
            CALL DETNLXR_PLUS_U( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, &
            X_NLOC, CV_NLOC, CV_NGI, &
            CVFEN, CVFENLX_ALL(1,:,:), CVFENLX_ALL(2,:,:), CVFENLX_ALL(3,:,:), CVWEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
            CVFENX_ALL, &
            U_NLOC, UFENLX_ALL(1,:,:), UFENLX_ALL(2,:,:), UFENLX_ALL(3,:,:), UFENX_ALL , &
            state ,"C_1", StorageIndexes(25))


            ! Adjust the volume according to the number of levels.
            VOLUME = VOLUME / REAL( NLEV )
            MASS_ELE( ELE ) = VOLUME


            ! *********subroutine Determine local vectors...

            LOC_U_RHS = 0.0

            DO ILEV = 1, NLEV
                DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
                    U_INOD = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_ILOC )
                    DO IPHASE = 1, NPHASE
                        DO IDIM = 1, NDIM_VEL
                            LOC_U( IDIM, IPHASE, U_ILOC ) = U_ALL( IDIM, IPHASE, U_INOD )
                            LOC_UOLD( IDIM, IPHASE, U_ILOC ) = UOLD_ALL( IDIM, IPHASE, U_INOD )
                            LOC_U_SOURCE( IDIM, IPHASE, U_ILOC ) = U_SOURCE( IDIM, IPHASE, U_INOD )
     IF(RETRIEVE_SOLID_CTY) LOC_US( IDIM, IPHASE, U_ILOC ) = us_all%val( IDIM, U_INOD )
                        END DO
                    END DO
                END DO
            END DO

            DO ILEV = 1, NLEV
                DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
                    U_INOD = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_ILOC )
                    DO IPHASE = 1, NPHASE
                        DO IDIM = 1, NDIM
                            LOC_NU( IDIM, IPHASE, U_ILOC ) = NU_ALL( IDIM, IPHASE, U_INOD )
                            LOC_NUOLD( IDIM, IPHASE, U_ILOC ) = NUOLD_ALL( IDIM, IPHASE, U_INOD )
                        END DO
                    END DO
                END DO
            END DO

            DO CV_ILOC = 1, CV_NLOC
                CV_INOD = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_ILOC )
                IF(IGOT_VOL_X_PRESSURE==1) THEN
                   LOC_UDEN( :, CV_ILOC ) = UDEN( :, CV_INOD ) * FEM_VOL_FRAC( :, CV_INOD )
                   LOC_UDENOLD( :, CV_ILOC) = UDENOLD( :, CV_INOD ) * FEM_VOL_FRAC( :, CV_INOD )
                ELSE
                   LOC_UDEN( :, CV_ILOC ) = UDEN( :, CV_INOD )
                   LOC_UDENOLD( :, CV_ILOC) = UDENOLD( :, CV_INOD )
                ENDIF

                DO IPHASE = 1, NPHASE
                    IF ( IPLIKE_GRAD_SOU /= 0 ) THEN
                        LOC_PLIKE_GRAD_SOU_COEF( IPHASE, CV_ILOC ) = PLIKE_GRAD_SOU_COEF( IPHASE, CV_INOD )
                        LOC_PLIKE_GRAD_SOU_GRAD( IPHASE, CV_ILOC ) = PLIKE_GRAD_SOU_GRAD( IPHASE, CV_INOD )
                    END IF
                    DO IDIM = 1, NDIM_VEL
                        LOC_U_SOURCE_CV( IDIM, IPHASE, CV_ILOC ) = U_SOURCE_CV( IDIM, IPHASE, CV_INOD )
                    END DO
                END DO
            END DO

            DO P_ILOC = 1, P_NLOC
                P_INOD = P_NDGLN( ( ELE - 1 ) * P_NLOC + P_ILOC )
                LOC_P( P_ILOC ) = P( P_INOD )
            END DO

            LOC_U_ABSORB = 0.0
            IF(RETRIEVE_SOLID_CTY) LOC_U_ABS_STAB_SOLID_RHS=0.0

            DO MAT_ILOC = 1, MAT_NLOC
                MAT_INOD = MAT_NDGLN( ( ELE - 1 ) * MAT_NLOC + MAT_ILOC )

                IF(is_compact_overlapping) THEN ! Set to the identity - NOT EFFICIENT BUT GOOD ENOUGH FOR NOW AS ITS SIMPLE...
                   DO I=1,NDIM_VEL* NPHASE
                      LOC_U_ABSORB( I, I, MAT_ILOC ) = 1.0 
                   END DO
                ELSE

                   LOC_U_ABSORB( :, :, MAT_ILOC ) = U_ABSORB( :, :, MAT_INOD )

! Switch on for solid fluid-coupling...
                   IF(RETRIEVE_SOLID_CTY) THEN
                      CV_INOD = CV_NDGLN( ( ELE - 1 ) * MAT_NLOC + MAT_ILOC )
                      DO IDIM=1,NDIM
                         DO IPHASE=1,NPHASE
                            I=IDIM + (IPHASE-1)*NDIM
                            LOC_U_ABSORB( I, I, MAT_ILOC ) = LOC_U_ABSORB( I, I, MAT_ILOC ) + &
                                   COEFF_SOLID_FLUID * ( DEN_ALL( IPHASE, cv_inod ) / dt ) * sf%val( cv_inod )

                            LOC_U_ABS_STAB_SOLID_RHS( I, I, MAT_ILOC ) = LOC_U_ABS_STAB_SOLID_RHS( I, I, MAT_ILOC )  &
                                  + COEFF_SOLID_FLUID * ( DEN_ALL( IPHASE, cv_inod ) / dt ) 
                         END DO
                      END DO
                   ENDIF

                END IF
                LOC_U_ABS_STAB( :, :, MAT_ILOC ) = U_ABS_STAB( :, :, MAT_INOD )

! Switch on for solid fluid-coupling apply stabilization term...
!                   IF(RETRIEVE_SOLID_CTY) THEN
!                      CV_INOD = CV_NDGLN( ( ELE - 1 ) * MAT_NLOC + MAT_ILOC )
!                      DO IDIM=1,NDIM
!                         DO IPHASE=1,NPHASE
!                            I=IDIM + (IPHASE-1)*NDIM
!                            LOC_U_ABS_STAB( I, I, MAT_ILOC ) = LOC_U_ABS_STAB( I, I, MAT_ILOC ) + &
!                                   COEFF_SOLID_FLUID * ( DEN_ALL( IPHASE, cv_inod ) / dt ) * sf%val( cv_inod )
!   
!                            LOC_U_ABS_STAB_SOLID_RHS( I, I, MAT_ILOC ) = LOC_U_ABS_STAB_SOLID_RHS( I, I, MAT_ILOC )  &
!                                  + COEFF_SOLID_FLUID * ( DEN_ALL( IPHASE, cv_inod ) / dt ) 
!                         END DO
!                      END DO
!                        
!                   ENDIF


                IF ( GOT_DIFFUS ) THEN
                   LOC_UDIFFUSION( :, :, :, MAT_ILOC ) = UDIFFUSION_ALL( :, :, :, MAT_INOD )
                ELSE
                   LOC_UDIFFUSION( :, :, :, MAT_ILOC ) = 0.0
                ENDIF
            END DO

            ! *********subroutine Determine local vectors...

            UD = 0.0 ; UDOLD = 0.0
            UD_ND = 0.0 ; UDOLD_ND = 0.0

            DO ILEV = 1, NLEV
                DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
                    DO GI = 1 + (ILEV-1)*CV_NGI_SHORT, ILEV*CV_NGI_SHORT
                        UD( :, :, GI ) = UD( :, :, GI ) + UFEN( U_ILOC, GI ) * LOC_NU( :, :, U_ILOC )
                        UDOLD( :, :, GI ) = UDOLD( :, :, GI ) + UFEN( U_ILOC, GI ) * LOC_NUOLD( :, :, U_ILOC )
                    END DO
                END DO
            END DO
            UD_ND( 1:NDIM_VEL, :, : ) = UD
            UDOLD_ND( 1:NDIM_VEL, :, : ) = UDOLD


            IF(IDIVID_BY_VOL_FRAC+IGOT_VOL_X_PRESSURE.GE.1) THEN

               VOL_FRA_GI_DX_ALL=0.0
               VOL_FRA_GI=0.0
               DO CV_ILOC = 1, CV_NLOC
                  CV_INOD = CV_NDGLN( (ELE-1)*CV_NLOC + CV_ILOC ) 
                  DO GI = 1, CV_NGI_SHORT
                     DO IPHASE=1,NPHASE
                        VOL_FRA_GI( IPHASE, GI )           = VOL_FRA_GI( IPHASE,GI )            + CVFEN_SHORT( CV_ILOC, GI )       * FEM_VOL_FRAC( IPHASE, CV_INOD )
                        VOL_FRA_GI_DX_ALL( :, IPHASE, GI ) = VOL_FRA_GI_DX_ALL( :, IPHASE, GI ) + CVFENX_ALL( 1:NDIM, CV_ILOC, GI )* FEM_VOL_FRAC( IPHASE, CV_INOD )
                     END DO
                  END DO
               END DO
               VOL_FRA_GI=MAX(VOL_FRA_GI, 0.0)

            ENDIF


            DENGI = 0.0 ; DENGIOLD = 0.0
            GRAD_SOU_GI = 0.0
            DO CV_ILOC = 1, CV_NLOC
                DO GI = 1, CV_NGI_SHORT
                    IF ( .FALSE. ) then ! FEM DEN...
                        DENGI( :, GI ) = DENGI( :, GI ) + CVFEN_SHORT( CV_ILOC, GI ) * LOC_UDEN( :, CV_ILOC )
                        DENGIOLD( :, GI ) = DENGIOLD( :, GI ) &
                        + CVFEN_SHORT( CV_ILOC, GI ) * LOC_UDENOLD( :, CV_ILOC )
                    ELSE ! CV DEN...
                        DENGI( :, GI ) = DENGI( :, GI ) + CVN_SHORT( CV_ILOC, GI ) * LOC_UDEN( :, CV_ILOC )
                        DENGIOLD( :, GI ) = DENGIOLD( :, GI ) &
                        + CVN_SHORT( CV_ILOC, GI ) * LOC_UDENOLD( :, CV_ILOC )
                    END IF
                    IF ( IPLIKE_GRAD_SOU == 1 ) THEN
                        GRAD_SOU_GI( :, GI ) = GRAD_SOU_GI( :, GI ) &
                        + CVFEN_SHORT( CV_ILOC, GI ) * LOC_PLIKE_GRAD_SOU_COEF( :, CV_ILOC )
                    END IF
                END DO
            END DO

            ! Start filtering density
            !FILT_DEN = 1
            !FILT_DEN = 2 ! best option to use
            FILT_DEN = 0
            IF ( FILT_DEN /= 0 ) THEN ! Filter the density...
                DENGI = 0.0 ; DENGIOLD = 0.0
                MASS_U = 0.0 ; MASS_U_CV = 0.0
                DO U_ILOC = 1, U_NLOC
                    DO U_JLOC = 1, U_NLOC
                        NN = SUM( UFEN( U_ILOC, : ) * UFEN( U_JLOC, : ) * DETWEI(:) )
                        IF ( FILT_DEN==2 ) THEN ! Lump the mass matrix for the filter - positive density...
                            MASS_U( U_ILOC, U_ILOC ) = MASS_U( U_ILOC, U_ILOC ) + NN
                        ELSE
                            MASS_U( U_ILOC, U_JLOC ) = MASS_U( U_ILOC, U_JLOC ) + NN
                        END IF
                    END DO
                END DO
                DO U_ILOC = 1, U_NLOC
                    DO CV_JLOC = 1, CV_NLOC
                        NCVM = SUM( UFEN( U_ILOC, : ) * CVN_SHORT( CV_JLOC, : ) * DETWEI(:) )
                        MASS_U_CV( U_ILOC, CV_JLOC ) = MASS_U_CV( U_ILOC, CV_JLOC ) + NCVM
                    END DO
                END DO

                STORE_MASS_U=MASS_U
                ! Store the LU decomposition...
                GOTDEC = .FALSE.

                RHS_U_CV = 0.0 ; RHS_U_CV_OLD = 0.0
                DO CV_JLOC = 1, CV_NLOC
                    DO U_ILOC = 1, U_NLOC
                        RHS_U_CV( :, U_ILOC ) = RHS_U_CV( :, U_ILOC ) + MASS_U_CV( U_ILOC, CV_JLOC ) * LOC_UDEN( :, CV_JLOC )
                        RHS_U_CV_OLD( :, U_ILOC ) = RHS_U_CV_OLD( :, U_ILOC ) + MASS_U_CV( U_ILOC, CV_JLOC ) * LOC_UDENOLD( :, CV_JLOC )
                    END DO
                END DO

                DO IPHASE = 1, NPHASE
                    CALL SMLINNGOT( STORE_MASS_U, UDEN_VFILT( IPHASE, : ), RHS_U_CV( IPHASE, : ), U_NLOC, U_NLOC, IPIV, GOTDEC )
                    GOTDEC = .TRUE.
                    CALL SMLINNGOT( STORE_MASS_U, UDENOLD_VFILT( IPHASE, : ), RHS_U_CV_OLD( IPHASE, : ), U_NLOC, U_NLOC, IPIV, GOTDEC )
                END DO

                DO U_ILOC = 1, U_NLOC
                    DO GI = 1, CV_NGI_SHORT
                        DENGI( :, GI ) = DENGI( :, GI ) + UFEN( U_ILOC, GI ) * UDEN_VFILT( :, U_ILOC )
                        DENGIOLD( :, GI ) = DENGIOLD( :, GI ) + UFEN( U_ILOC, GI ) * UDENOLD_VFILT( :, U_ILOC )
                    END DO
                END DO
            END IF

            ! not good to have -ve density at quadature pt...
            DENGI = MAX( 0.0, DENGI )
            DENGIOLD = MAX( 0.0, DENGIOLD )

            SIGMAGI = 0.0 ; SIGMAGI_STAB = 0.0
            IF(RETRIEVE_SOLID_CTY) SIGMAGI_STAB_SOLID_RHS=0.0
            TEN_XX  = 0.0
            if (is_compact_overlapping) then
               DO IPHA_IDIM = 1, NDIM_VEL * NPHASE
                    SIGMAGI( IPHA_IDIM, IPHA_IDIM, : ) = 1.0
               end do
            else
                DO MAT_ILOC = 1, MAT_NLOC
                    DO GI = 1, CV_NGI
                        DO IPHA_IDIM = 1, NDIM_VEL * NPHASE
                            DO JPHA_JDIM = 1, NDIM_VEL * NPHASE

                   IF(RETRIEVE_SOLID_CTY) THEN
                                SIGMAGI( IPHA_IDIM, JPHA_JDIM, GI ) = SIGMAGI( IPHA_IDIM, JPHA_JDIM, GI ) &
                                      !+ CVFEN( MAT_ILOC, GI ) * LOC_U_ABSORB( IPHA_IDIM, JPHA_JDIM, MAT_ILOC )
!                                + CVN( MAT_ILOC, GI ) * LOC_U_ABSORB( IPHA_IDIM, JPHA_JDIM, MAT_ILOC )
                                + CVFEN( MAT_ILOC, GI ) * LOC_U_ABSORB( IPHA_IDIM, JPHA_JDIM, MAT_ILOC )

                                SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM, GI ) = SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM, GI ) &
                                      !+ CVFEN( MAT_ILOC, GI ) * LOC_U_ABS_STAB( IPHA_IDIM, JPHA_JDIM, MAT_ILOC )
                                + CVFEN( MAT_ILOC, GI ) * LOC_U_ABS_STAB( IPHA_IDIM, JPHA_JDIM, MAT_ILOC )

                                SIGMAGI_STAB_SOLID_RHS( IPHA_IDIM, JPHA_JDIM, GI ) = SIGMAGI_STAB_SOLID_RHS( IPHA_IDIM, JPHA_JDIM, GI ) &
                           !     + CVN( MAT_ILOC, GI ) * LOC_U_ABS_STAB_SOLID_RHS( IPHA_IDIM, JPHA_JDIM, MAT_ILOC ) 
                                + CVFEN( MAT_ILOC, GI ) * LOC_U_ABS_STAB_SOLID_RHS( IPHA_IDIM, JPHA_JDIM, MAT_ILOC ) 

                                SIGMAGI( IPHA_IDIM, JPHA_JDIM, GI ) = max(0.0, SIGMAGI( IPHA_IDIM, JPHA_JDIM, GI ) ) 
                                SIGMAGI_STAB_SOLID_RHS( IPHA_IDIM, JPHA_JDIM, GI ) = max(0.0, SIGMAGI_STAB_SOLID_RHS( IPHA_IDIM, JPHA_JDIM, GI ) )
                   ELSE
                                SIGMAGI( IPHA_IDIM, JPHA_JDIM, GI ) = SIGMAGI( IPHA_IDIM, JPHA_JDIM, GI ) &
                                      !+ CVFEN( MAT_ILOC, GI ) * LOC_U_ABSORB( IPHA_IDIM, JPHA_JDIM, MAT_ILOC )
                                + CVN( MAT_ILOC, GI ) * LOC_U_ABSORB( IPHA_IDIM, JPHA_JDIM, MAT_ILOC )

                                SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM, GI ) = SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM, GI ) &
                                      !+ CVFEN( MAT_ILOC, GI ) * LOC_U_ABS_STAB( IPHA_IDIM, JPHA_JDIM, MAT_ILOC )
                                + CVN( MAT_ILOC, GI ) * LOC_U_ABS_STAB( IPHA_IDIM, JPHA_JDIM, MAT_ILOC )
                   ENDIF

                            END DO
                        END DO
                        TEN_XX( :, :, :, GI ) = TEN_XX( :, :, :, GI ) + CVFEN( MAT_ILOC, GI ) * LOC_UDIFFUSION( :, :, :, MAT_ILOC )
                    END DO
                END DO
            end if
            RHS_DIFF_U=0.0

            Loop_ilev_DGNods1: DO ILEV = 1, NLEV
                NN_SIGMAGI_ELE = 0.0
                NN_SIGMAGI_STAB_ELE = 0.0
                IF(RETRIEVE_SOLID_CTY) NN_SIGMAGI_STAB_SOLID_RHS_ELE = 0.0
                NN_MASS_ELE = 0.0
                NN_MASSOLD_ELE = 0.0
                VLK_ELE = 0.0
                STRESS_IJ_ELE = 0.0
                !Prepare data
                DO U_JLOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2

!                    IF ( STAB_VISC_WITH_ABS ) THEN
                    IF ( GOT_DIFFUS ) THEN

                        DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
                            DO GI = 1 + (ILEV-1)*CV_NGI_SHORT, ILEV*CV_NGI_SHORT
                                DO IPHASE = 1, NPHASE
                                    IF ( STRESS_FORM ) THEN ! stress form of viscosity...
                                        IF(IDIVID_BY_VOL_FRAC==1) THEN
                                           CALL CALC_STRESS_TEN( STRESS_IJ_ELE( :, :, IPHASE, U_ILOC, U_JLOC ), ZERO_OR_TWO_THIRDS, NDIM, &
        ( -UFEN( U_ILOC, GI )*VOL_FRA_GI_DX_ALL(1:NDIM,IPHASE,GI) + UFENX_ALL( 1:NDIM, U_ILOC, GI )*VOL_FRA_GI(IPHASE,GI) ),  UFENX_ALL( 1:NDIM, U_JLOC, GI )* DETWEI( GI ), TEN_XX( :, :, IPHASE, GI ) )
                                        ELSE
                                           CALL CALC_STRESS_TEN( STRESS_IJ_ELE( :, :, IPHASE, U_ILOC, U_JLOC ), ZERO_OR_TWO_THIRDS, NDIM, &
                                           UFENX_ALL( 1:NDIM, U_ILOC, GI ), UFENX_ALL( 1:NDIM, U_JLOC, GI )* DETWEI( GI ), TEN_XX( :, :, IPHASE, GI ) )
                                        ENDIF
                                    ELSE
                                        DO IDIM = 1, NDIM
                                            VLK_ELE( IPHASE, U_ILOC, U_JLOC ) = VLK_ELE( IPHASE, U_ILOC, U_JLOC ) + &
                                            UFENX_ALL( IDIM, U_ILOC, GI ) * SUM( UFENX_ALL( 1:NDIM, U_JLOC, GI ) * TEN_XX( IDIM, :, IPHASE, GI ) ) * DETWEI( GI )
                                        END DO
                                    END IF
                                END DO
                            END DO
                        END DO
                    END IF


                    DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
                        DO GI = 1 + (ILEV-1)*CV_NGI_SHORT, ILEV*CV_NGI_SHORT

                            GI_SHORT = GI - (ILEV-1)*CV_NGI_SHORT

                            RNN = UFEN( U_ILOC, GI ) * UFEN( U_JLOC,  GI ) * DETWEI( GI )

                            DO JPHASE = 1, NPHASE
                                DO JDIM = 1, NDIM_VEL
                                    JPHA_JDIM = JDIM + (JPHASE-1)*NDIM

                                    NN_MASS_ELE( JPHA_JDIM, U_ILOC, JPHA_JDIM, U_JLOC ) = NN_MASS_ELE( JPHA_JDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                    + DENGI(JPHASE,GI_SHORT) * RNN
                                    NN_MASSOLD_ELE( JPHA_JDIM, U_ILOC, JPHA_JDIM, U_JLOC ) = NN_MASSOLD_ELE( JPHA_JDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                    + DENGIOLD(JPHASE, GI_SHORT) * RNN

                                    ! Stabilization for viscosity...
                                    IF ( STAB_VISC_WITH_ABS ) THEN
                                        IF ( STRESS_FORM ) THEN
                                            NN_SIGMAGI_STAB_ELE( JPHA_JDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                            = NN_SIGMAGI_STAB_ELE( JPHA_JDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                            + MAX( 0.0, STRESS_IJ_ELE( JDIM, JDIM, JPHASE, U_ILOC, U_JLOC ) )
                                        ELSE
                                            NN_SIGMAGI_STAB_ELE( JPHA_JDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                            = NN_SIGMAGI_STAB_ELE( JPHA_JDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                            + MAX( 0.0, VLK_ELE( JPHASE, U_ILOC, U_JLOC ) )
                                        END IF
                                    END IF


                                    DO IPHASE = 1, NPHASE
                                        DO IDIM = 1, NDIM_VEL
                                            IPHA_IDIM = IDIM + (IPHASE-1)*NDIM

                                            NN_SIGMAGI_ELE( IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                            = NN_SIGMAGI_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) + RNN *  &
                                            SIGMAGI( IPHA_IDIM, JPHA_JDIM, GI )

                                            NN_SIGMAGI_STAB_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                            = NN_SIGMAGI_STAB_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) + RNN *  &
                                            SIGMAGI_STAB( IPHA_IDIM, JPHA_JDIM, GI )
                   IF(RETRIEVE_SOLID_CTY) THEN
                                            NN_SIGMAGI_STAB_SOLID_RHS_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                            = NN_SIGMAGI_STAB_SOLID_RHS_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) + RNN *  &
                                            SIGMAGI_STAB_SOLID_RHS( IPHA_IDIM, JPHA_JDIM, GI )
                   ENDIF
                                        END DO
                                    END DO
                                END DO

                            END DO
                        END DO
                    END DO

                END DO

                DO U_JLOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
                    DO JPHASE = 1, NPHASE
                        DO JDIM = 1, NDIM_VEL
                            JPHA_JDIM = JDIM + (JPHASE-1)*NDIM
                            J = JDIM+(JPHASE-1)*NDIM_VEL+(U_JLOC-1)*NDIM_VEL*NPHASE
                            DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
                                DO IPHASE = 1, NPHASE
                                    DO IDIM = 1, NDIM_VEL
                                        IPHA_IDIM = IDIM + (IPHASE-1)*NDIM
                                        I = IDIM+(IPHASE-1)*NDIM_VEL+(U_ILOC-1)*NDIM_VEL*NPHASE
                                        !Assemble
                                        IF ( LUMP_MASS ) THEN
                                            PIVIT_MAT( I, I, ELE ) =   &
                                            NN_SIGMAGI_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                            + NN_SIGMAGI_STAB_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                            + NN_MASS_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC )/DT

                                        ELSE
                                            PIVIT_MAT( I, J, ELE ) =  &
                                            NN_SIGMAGI_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                            + NN_SIGMAGI_STAB_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                            + NN_MASS_ELE(IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC )/DT
                                        END IF

                                        IF ( .NOT.NO_MATRIX_STORE ) THEN
                                            IF ( .NOT.JUST_BL_DIAG_MAT ) THEN
                                                IF ( LUMP_MASS ) THEN
                                                    DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_ILOC, ELE ) =  &
                                                    DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_ILOC, ELE )  &
                                                    + NN_SIGMAGI_ELE( IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                                    + NN_SIGMAGI_STAB_ELE( IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                                    + NN_MASS_ELE( IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) / DT
                                                ELSE
                                                    DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) = &
                                                    DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE )  &
                                                    + NN_SIGMAGI_ELE( IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                                    + NN_SIGMAGI_STAB_ELE( IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) &
                                                    + NN_MASS_ELE( IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) / DT
                                                END IF
                                            END IF
                                        END IF

                                    END DO
                                END DO
                            END DO
                        END DO
                    END DO
                END DO


                Loop_DGNods1: DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
                    !GLOBI = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_ILOC )
                    !IF ( NLEV==1 .AND. LUMP_MASS ) GLOBI_CV = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + U_ILOC )
                    ! put CV source in...
                    Loop_CVNods2: DO CV_JLOC = 1, CV_NLOC
                        NM = SUM( UFEN( U_ILOC, : ) * CVN( CV_JLOC,  : ) * DETWEI( : ) )
                        IF ( LUMP_MASS ) THEN
                            ! This is bugged...fix me!
                            STOP 6378
                            IF ( CV_NLOC==6 .OR. (CV_NLOC==10 .AND. NDIM==3) ) THEN
                                IF ( CV_JLOC==1 .OR. CV_JLOC==3 .OR. CV_JLOC==6 .OR. CV_JLOC==10 ) THEN
                                    CV_ILOC = ELEMENT_CORNERS( U_ILOC )
                                    LOC_U_RHS( :, :, U_ILOC ) = LOC_U_RHS( :, :, U_ILOC ) + NM * LOC_U_SOURCE_CV( :, :, CV_ILOC )
                                END IF
                            ELSE
                                CV_ILOC = U_ILOC
                                LOC_U_RHS( :, :, U_ILOC ) = LOC_U_RHS( :, :, U_ILOC ) + NM * LOC_U_SOURCE_CV( :, :, CV_ILOC )
                            END IF
                        ELSE
                            LOC_U_RHS( :, :, U_ILOC ) = LOC_U_RHS( :, :, U_ILOC ) + NM * LOC_U_SOURCE_CV( :, :, CV_JLOC )
                        END IF
                    END DO LOOP_CVNODS2


            !        IF ( LUMP_MASS ) THEN
            !           IF ( CV_NLOC==6 .OR. (CV_NLOC==10 .AND. NDIM==3) ) THEN
            !              CV_ILOC = ELEMENT_CORNERS( U_ILOC )
            !           ELSE
            !              CV_ILOC = U_ILOC
            !           END IF
            !           DO CV_JLOC = 1, CV_NLOC
            !              NM = SUM( UFEN( U_ILOC, : ) * CVN( CV_JLOC,  : ) * DETWEI( : ) )
            !              LOC_U_RHS( :, :, U_ILOC ) = LOC_U_RHS( :, :, U_ILOC ) + NM * LOC_U_SOURCE_CV( :, :, CV_ILOC )
            !           END DO
            !        ELSE
            !           DO CV_JLOC = 1, CV_NLOC
            !              NM = SUM( UFEN( U_ILOC, : ) * CVN( CV_JLOC,  : ) * DETWEI( : ) )
            !              LOC_U_RHS( :, :, U_ILOC ) = LOC_U_RHS( :, :, U_ILOC ) + NM * LOC_U_SOURCE_CV( :, :, CV_JLOC )
            !           END DO
            !        END IF


                    Loop_DGNods2: DO U_JLOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2

                        NN = 0.0
                        VLN = 0.0
                        VLN_OLD = 0.0
!                        VLK = 0.0

                        Loop_Gauss2: DO GI = 1 + (ILEV-1)*CV_NGI_SHORT, ILEV*CV_NGI_SHORT

                            RNN = UFEN( U_ILOC, GI ) * UFEN( U_JLOC,  GI ) * DETWEI( GI )
                            NN = NN + RNN

                            Loop_IPHASE: DO IPHASE = 1, NPHASE ! Diffusion tensor

                                IF ( MOM_CONSERV ) THEN
                                    VLN( IPHASE ) = VLN( IPHASE ) - &
                                    DENGI( IPHASE, GI ) * SUM( UD( :, IPHASE, GI ) * UFENX_ALL( 1:NDIM, U_ILOC, GI ) )  &
                                    * UFEN( U_JLOC, GI ) * DETWEI( GI ) * WITH_NONLIN

                                    VLN_OLD( IPHASE ) = VLN_OLD( IPHASE ) - &
                                    DENGI( IPHASE, GI ) * SUM( UDOLD( :, IPHASE, GI ) * UFENX_ALL( 1:NDIM, U_ILOC, GI ) )  &
                                    * UFEN( U_JLOC, GI ) * DETWEI( GI ) * WITH_NONLIN
                                ELSE
                                    VLN( IPHASE ) = VLN( IPHASE ) + &
                                    UFEN( U_ILOC, GI ) * DENGI( IPHASE, GI ) * SUM( UD( :, IPHASE, GI ) * UFENX_ALL(1:NDIM, U_JLOC, GI ) ) &
                                    * DETWEI( GI ) * WITH_NONLIN

                                    VLN_OLD( IPHASE ) = VLN_OLD( IPHASE ) + &
                                    UFEN( U_ILOC, GI ) * DENGI( IPHASE, GI ) * SUM( UDOLD( :, IPHASE, GI ) * UFENX_ALL( 1:NDIM, U_JLOC, GI ) ) &
                                    * DETWEI( GI ) * WITH_NONLIN
                                END IF

!                                DO IDIM = 1, NDIM
!                                   VLK( IPHASE ) = VLK( IPHASE ) + &
!                                            UFENX_ALL( IDIM, U_ILOC, GI ) * SUM( UFENX_ALL( 1:NDIM, U_JLOC, GI ) * TEN_XX( IDIM, :, IPHASE, GI ) ) * DETWEI( GI )
!                                END DO

                            END DO Loop_IPHASE

                        END DO Loop_Gauss2

                        LOC_U_RHS( :, :, U_ILOC ) =  LOC_U_RHS( :, :, U_ILOC ) + NN * LOC_U_SOURCE( :, :, U_JLOC  )

                        DO JPHASE = 1, NPHASE
                            DO JDIM = 1, NDIM_VEL

                                JPHA_JDIM = (JPHASE-1)*NDIM_VEL + JDIM

                                DO IPHASE = 1, NPHASE
                                    DO IDIM = 1, NDIM_VEL

                                        IPHA_IDIM = (IPHASE-1)*NDIM_VEL + IDIM

                                        IF ( MOM_CONSERV ) THEN
                                            IF ( LUMP_MASS ) THEN
                                                LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                + NN_SIGMAGI_STAB_ELE( IPHA_IDIM, U_ILOC,JPHA_JDIM, U_JLOC ) * LOC_U( JDIM, JPHASE, U_JLOC )     &
                                                + ( NN_MASSOLD_ELE( IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) / DT ) * LOC_UOLD( JDIM, JPHASE, U_ILOC )
                                                IF(RETRIEVE_SOLID_CTY) LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                + NN_SIGMAGI_STAB_SOLID_RHS_ELE( IPHA_IDIM, U_ILOC,JPHA_JDIM, U_JLOC ) * LOC_US( JDIM, JPHASE, U_JLOC )    
                                            ELSE
                                                LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                + NN_SIGMAGI_STAB_ELE( IPHA_IDIM, U_ILOC,JPHA_JDIM, U_JLOC ) * LOC_U( JDIM, JPHASE, U_JLOC )  &
                                                + ( NN_MASSOLD_ELE( IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) / DT ) * LOC_UOLD( JDIM, JPHASE, U_JLOC )
                                                IF(RETRIEVE_SOLID_CTY) LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                + NN_SIGMAGI_STAB_SOLID_RHS_ELE( IPHA_IDIM, U_ILOC,JPHA_JDIM, U_JLOC ) * LOC_US( JDIM, JPHASE, U_JLOC )  
                                            END IF
                                        ELSE
                                            IF ( LUMP_MASS ) THEN
                                                LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                + NN_SIGMAGI_STAB_ELE( IPHA_IDIM, U_ILOC,JPHA_JDIM, U_JLOC ) * LOC_U( JDIM, JPHASE, U_JLOC )  &
                                                + ( NN_MASS_ELE( IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) / DT ) * LOC_UOLD( JDIM, JPHASE, U_ILOC )
                                                IF(RETRIEVE_SOLID_CTY) LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                + NN_SIGMAGI_STAB_SOLID_RHS_ELE( IPHA_IDIM, U_ILOC,JPHA_JDIM, U_JLOC ) * LOC_US( JDIM, JPHASE, U_JLOC ) 
                                            ELSE
                                                LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                + NN_SIGMAGI_STAB_ELE( IPHA_IDIM, U_ILOC,JPHA_JDIM, U_JLOC ) * LOC_U( JDIM, JPHASE, U_JLOC )  &
                                                + ( NN_MASS_ELE( IPHA_IDIM, U_ILOC, JPHA_JDIM, U_JLOC ) / DT ) * LOC_UOLD( JDIM, JPHASE, U_JLOC )
                                                IF(RETRIEVE_SOLID_CTY) LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                + NN_SIGMAGI_STAB_SOLID_RHS_ELE( IPHA_IDIM, U_ILOC,JPHA_JDIM, U_JLOC ) * LOC_US( JDIM, JPHASE, U_JLOC )  
                                            END IF
                                        END IF

                                    END DO
                                END DO
                            END DO
                        END DO

                        IF ( .NOT.JUST_BL_DIAG_MAT ) THEN
                            IF ( STRESS_FORM ) THEN
                                DO IPHASE = 1, NPHASE
                                    JPHASE = IPHASE
                                    DO IDIM = 1, NDIM_VEL
                                        DO JDIM = 1, NDIM_VEL

                                            IF ( NO_MATRIX_STORE ) THEN
                                                LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                                - STRESS_IJ_ELE( IDIM, JDIM,  IPHASE, U_ILOC, U_JLOC ) * LOC_U( JDIM, IPHASE, U_JLOC )
                                            ELSE
                                                DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE )  &
                                                = DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) &
                                                + STRESS_IJ_ELE( IDIM, JDIM, IPHASE, U_ILOC, U_JLOC )
                                            END IF

                                            RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) = RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) + &
                                            STRESS_IJ_ELE( IDIM, JDIM, IPHASE, U_ILOC, U_JLOC ) * LOC_U( JDIM, IPHASE, U_JLOC )
                                        END DO
                                    END DO

                                END DO
                            END IF

                            DO IDIM = 1, NDIM_VEL
                                DO IPHASE = 1, NPHASE
                                    JDIM = IDIM
                                    JPHASE = IPHASE

                                    IF ( NO_MATRIX_STORE ) THEN
                                        LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC )  &
                                        - VLN( IPHASE ) * LOC_U( IDIM, IPHASE, U_JLOC )
                                    ELSE
                                        DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) &
                                        = DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) + VLN( IPHASE )
                                    END IF

                                    IF ( .NOT.STRESS_FORM ) THEN
                                        IF ( NO_MATRIX_STORE ) THEN
                                            LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                            - VLK_ELE( IPHASE, U_ILOC, U_JLOC ) * LOC_U( IDIM, IPHASE, U_JLOC )
!                                            LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
!                                            - VLK( IPHASE ) * LOC_U( IDIM, IPHASE, U_JLOC )
                                        ELSE
                                            DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) &
                                            = DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) + VLK_ELE( IPHASE, U_ILOC, U_JLOC )
!                                            DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) &
!                                            = DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) + VLK( IPHASE )

                                        END IF

                                        RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) = RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) + &
                                        VLK_ELE( IPHASE, U_ILOC, U_JLOC ) * LOC_U( IDIM, IPHASE, U_JLOC )
!                                        RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) = RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) + &
!                                        VLK( IPHASE ) * LOC_U( IDIM, IPHASE, U_JLOC )
                                    END IF

                                END DO
                            END DO

                        END IF ! .NOT.JUST_BL_DIAG_MAT

                    END DO Loop_DGNods2

                END DO Loop_DGNods1
            END DO Loop_ilev_DGNods1
            ! **********REVIEWER 1-END**********************


            ! **********REVIEWER 2-START**********************
            !ewrite(3,*)'just after Loop_DGNods1'

            ! Add-in surface contributions.

            ! Find diffusion contributions at the surface
            !CALL DG_DIFFUSION( ELE, U_NLOC, U_NONODS, TOTELE, LMMAT1, LMMAT, LNXNMAT1, LNNXMAT, LINVMMAT1, &
            !LINVMNXNMAT1, AMAT )

            ! Add in C matrix contribution: (DG velocities)
            Loop_ILEV1: DO ILEV = 1, NLEV
                Loop_U_ILOC1: DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
                    if(.not.got_c_matrix) IU_NOD = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_ILOC )

                    Loop_P_JLOC1: DO P_JLOC = 1, P_NLOC
                        if(.not.got_c_matrix) JCV_NOD = P_NDGLN( ( ELE - 1 ) * P_NLOC + P_JLOC )

                        NMX_ALL = 0.0
                        GRAD_SOU_GI_NMX = 0.0
                        IF(IGOT_VOL_X_PRESSURE==1) VOL_FRA_NMX_ALL(:,:) = 0.0
                        Loop_GaussPoints1: DO GI = 1 + (ILEV-1)*CV_NGI_SHORT, ILEV*CV_NGI_SHORT
                            RN = UFEN( U_ILOC, GI ) * DETWEI( GI )

                            RNMX_ALL( : ) = RN * CVFENX_ALL( 1:NDIM, P_JLOC, GI )

                            NMX_ALL( : ) = NMX_ALL( : ) + RNMX_ALL( : )

                            IF ( IPLIKE_GRAD_SOU == 1 .OR. CAPILLARY_PRESSURE_ACTIVATED ) THEN
                                DO IDIM = 1, NDIM_VEL
                                    GRAD_SOU_GI_NMX( IDIM, : ) = GRAD_SOU_GI_NMX( IDIM, : ) &
                                    + GRAD_SOU_GI( :, GI ) * RNMX_ALL( IDIM )
                                END DO
                            END IF
                            IF(IGOT_VOL_X_PRESSURE==1) THEN
                               DO IPHASE = 1, NPHASE
                                  VOL_FRA_NMX_ALL( :, IPHASE ) = VOL_FRA_NMX_ALL( :, IPHASE ) + VOL_FRA_GI( IPHASE,GI ) * RNMX_ALL( : )
                               END DO
                            ENDIF

                        END DO Loop_GaussPoints1

                        ! Put into matrix
                        IF ( .NOT.GOT_C_MATRIX ) THEN
                            CALL USE_POSINMAT_C_STORE( COUNT, IU_NOD, JCV_NOD,  &
                            U_NONODS, FINDC, COLC, NCOLC, &
                            IDO_STORE_AC_SPAR_PT, STORED_AC_SPAR_PT, POSINMAT_C_STORE, ELE, U_ILOC, P_JLOC, &
                            TOTELE, U_NLOC, P_NLOC )
                        END IF

                        Loop_Phase1: DO IPHASE = 1, NPHASE

                            ! Put into matrix
                            IF ( .NOT.GOT_C_MATRIX ) THEN
                                DO IDIM = 1, NDIM_VEL
                                    IF(IGOT_VOL_X_PRESSURE==1) THEN
                                       C( IDIM, IPHASE, COUNT ) = C( IDIM, IPHASE, COUNT ) - VOL_FRA_NMX_ALL( IDIM, IPHASE )
                                    ELSE
                                       C( IDIM, IPHASE, COUNT ) = C( IDIM, IPHASE, COUNT ) - NMX_ALL( IDIM )
                                    ENDIF
                                END DO
                            END IF

                            IF ( IPLIKE_GRAD_SOU == 1 .OR. CAPILLARY_PRESSURE_ACTIVATED ) THEN ! Capillary pressure for example terms...
                                DO IDIM = 1, NDIM_VEL
                                    LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                    - GRAD_SOU_GI_NMX( IDIM, IPHASE ) * LOC_PLIKE_GRAD_SOU_GRAD( IPHASE, P_JLOC )
                                END DO

                            END IF
                        END DO Loop_Phase1

                    END DO Loop_P_JLOC1

                END DO Loop_U_ILOC1
            END DO Loop_ILEV1

            !ewrite(3,*)'just after Loop_U_ILOC1'

            ! **********REVIEWER 2-END**********************


            ! **********REVIEWER 2-START**********************

            IF( (.NOT.FIRSTST) .AND. (RESID_BASED_STAB_DIF/=0) ) THEN
                !! *************************INNER ELEMENT STABILIZATION****************************************
                !! *************************INNER ELEMENT STABILIZATION****************************************

                DO U_ILOC = 1, U_NLOC
                    DO U_JLOC = 1, U_NLOC
                        ! Sum over quadrature pts...
                        LOC_MASS( U_ILOC, U_JLOC ) = SUM( UFEN( U_ILOC, : ) * UFEN( U_JLOC,  : ) * DETWEI( : ) )
                    END DO
                END DO

                LOC_MASS_INV = LOC_MASS
                !CALL INVERT(LOC_MASS_INV)
                CALL MATDMATINV( LOC_MASS, LOC_MASS_INV, U_NLOC )

                DO U_ILOC = 1, U_NLOC
                    DO IPHASE = 1, NPHASE
                        DO IDIM = 1, NDIM_VEL
                            ! sum cols of matrix * rows of vector...
                            DIFF_VEC_U( IDIM, IPHASE, U_ILOC ) = SUM( LOC_MASS_INV( U_ILOC, : ) * RHS_DIFF_U( IDIM, IPHASE, : ) )
                        END DO
                    END DO
                END DO

                DIFFGI_U = 0.0
                U_DX_ALL = 0.0 ; UOLD_DX_ALL = 0.0
                SOUGI_X = 0.0

                DO U_ILOC = 1, U_NLOC
                    DO GI = 1, CV_NGI
                        DO IPHASE = 1, NPHASE

                            DIFFGI_U( :, IPHASE, GI ) = DIFFGI_U( :, IPHASE, GI ) + &
                            UFEN( U_ILOC, GI ) * DIFF_VEC_U( :, IPHASE, U_ILOC )

                            SOUGI_X( :, IPHASE, GI ) = SOUGI_X( :, IPHASE, GI ) + &
                            UFEN( U_ILOC, GI ) * LOC_U_SOURCE( :, IPHASE, U_ILOC )

                            DO JDIM = 1, NDIM_VEL
                                DO IDIM = 1, NDIM
                                    U_DX_ALL( IDIM, JDIM, IPHASE, GI ) = U_DX_ALL( IDIM, JDIM, IPHASE, GI ) + &
                                    LOC_U( JDIM, IPHASE, U_ILOC ) * UFENX_ALL( IDIM, U_ILOC, GI )
                                    UOLD_DX_ALL( IDIM, JDIM, IPHASE, GI ) = UOLD_DX_ALL( IDIM, JDIM, IPHASE, GI ) + &
                                    LOC_UOLD( JDIM, IPHASE, U_ILOC ) * UFENX_ALL( IDIM, U_ILOC, GI )
                                END DO
                            END DO

                        END DO
                    END DO
                END DO

                U_DT = ( UD - UDOLD ) / DT

                RESID_U = 0.0
                DO GI = 1, CV_NGI
                    DO IPHASE = 1, NPHASE
                        DO IDIM = 1, NDIM_VEL
                            IPHA_IDIM = (IPHASE-1)*NDIM_VEL + IDIM
                            DO JPHASE = 1, NPHASE
                                DO JDIM = 1, NDIM_VEL
                                    JPHA_JDIM = (JPHASE-1)*NDIM_VEL + JDIM
                                    RESID_U( IDIM, IPHASE, GI ) = RESID_U( IDIM, IPHASE, GI ) + &
                                    SIGMAGI( IPHA_IDIM, JPHA_JDIM, GI ) * UD( JDIM, IPHASE, GI )
                                END DO
                            END DO
                        END DO
                    END DO
                END DO


                P_DX = 0.0

                DO P_ILOC = 1, P_NLOC
                    DO GI = 1, CV_NGI

                        P_DX( :, GI ) = P_DX( :, GI ) + CVFENX_ALL(1:NDIM, P_ILOC, GI ) * LOC_P( P_ILOC )

                        IF ( IPLIKE_GRAD_SOU == 1 .OR. CAPILLARY_PRESSURE_ACTIVATED ) THEN ! Capillary pressure for example terms...
                            DO IPHASE = 1, NPHASE

                                R = GRAD_SOU_GI( IPHASE, GI ) * LOC_PLIKE_GRAD_SOU_GRAD( IPHASE, P_ILOC )
                                DO IDIM = 1, NDIM_VEL
                                    RESID_U( IDIM, IPHASE, GI ) = RESID_U( IDIM, IPHASE, GI ) + R * CVFENX_ALL( IDIM, P_ILOC, GI )
                                END DO

                            END DO
                        END IF
                    END DO
                END DO



                DO GI = 1, CV_NGI
                    DO IPHASE = 1, NPHASE

                        DO IDIM = 1, NDIM_VEL
                            RESID_U( IDIM, IPHASE, GI ) = RESID_U( IDIM, IPHASE, GI ) + &
                            DENGI( IPHASE, GI ) * SUM( UD_ND( :, IPHASE, GI ) * U_DX_ALL( :, IDIM, IPHASE, GI ) ) &
                            * WITH_NONLIN &
                            + DENGI( IPHASE, GI ) * U_DT( IDIM, IPHASE, GI )   &
                            - SOUGI_X( IDIM, IPHASE, GI ) - DIFFGI_U( IDIM, IPHASE, GI ) + P_DX( IDIM, GI )

                            U_GRAD_NORM2( IDIM, IPHASE, GI ) = U_DT( IDIM, IPHASE, GI )**2 + SUM( U_DX_ALL( :, IDIM, IPHASE, GI )**2 )
                            U_GRAD_NORM( IDIM, IPHASE, GI ) = MAX( TOLER, SQRT( U_GRAD_NORM2( IDIM, IPHASE, GI ) ) )
                            U_GRAD_NORM2( IDIM, IPHASE, GI ) = MAX( TOLER, U_GRAD_NORM2( IDIM, IPHASE, GI ) )

                            A_DOT_U( IDIM, IPHASE, GI ) = DENGI( IPHASE, GI ) * ( SUM( UD_ND( :, IPHASE, GI ) * U_DX_ALL( :, IDIM, IPHASE, GI ) ) &
                            * WITH_NONLIN + U_DT( IDIM, IPHASE, GI ) ) + P_DX( IDIM, GI ) * RNO_P_IN_A_DOT

                            STAR_U_COEF( IDIM, IPHASE, GI ) = A_DOT_U( IDIM, IPHASE, GI ) / U_GRAD_NORM2( IDIM, IPHASE, GI )

                        END DO

                        JTT_INV = 2. / DT

                        U_GRAD_N_MAX2=0.0
                        DO U_ILOC = 1, U_NLOC
                            DO IDIM = 1, NDIM_VEL
                                U_GRAD_N_MAX2( IDIM ) = MAX( U_GRAD_N_MAX2( IDIM ), &
                                ( JTT_INV * U_DT( IDIM, IPHASE, GI ) )**2 &
                                + 4. * SUM( ( UFENX_ALL( 1:NDIM, U_ILOC, GI ) * U_DX_ALL( 1:NDIM, IDIM, IPHASE, GI ) )**2 ) )
                            END DO
                        END DO

                        DO IDIM = 1, NDIM_VEL
                            P_STAR_U( IDIM, IPHASE, GI ) = U_NONLIN_SHOCK_COEF / MAX( TOLER, SQRT( STAR_U_COEF( IDIM, IPHASE, GI )**2 * U_GRAD_N_MAX2( IDIM ) ) )
                        END DO

                        IF ( RESID_BASED_STAB_DIF==1 ) THEN

                            U_R2_COEF( : ) = RESID_U( :, IPHASE, GI )**2

                        ELSE IF ( RESID_BASED_STAB_DIF==2 ) THEN

                            U_R2_COEF( : ) = MAX( 0.0, A_DOT_U( :, IPHASE, GI ) * RESID_U( :, IPHASE, GI ) )

                        ELSE IF ( RESID_BASED_STAB_DIF==3 ) THEN ! Max of two previous methods.

                            U_R2_COEF( : ) = MAX( RESID_U( :, IPHASE, GI )**2, A_DOT_U( :, IPHASE, GI ) * RESID_U( :, IPHASE, GI ) )

                        END IF

                        DIF_STAB_U( :, IPHASE, GI ) = U_R2_COEF( : ) * P_STAR_U( :, IPHASE, GI ) / U_GRAD_NORM2( :, IPHASE, GI )

                    END DO
                END DO


                ! Place the diffusion term into matrix...
                DO U_ILOC = 1, U_NLOC
                    DO U_JLOC = 1, U_NLOC
                        DO IPHASE = 1, NPHASE
                            JPHASE = IPHASE

                            VLK_UVW = 0.0
                            DO GI = 1, CV_NGI
                                VLKNN = SUM( UFENX_ALL( 1:NDIM, U_ILOC, GI ) * UFENX_ALL( 1:NDIM, U_JLOC, GI ) ) * DETWEI( GI )
                                VLK_UVW( : ) = VLK_UVW( : ) + DIF_STAB_U( :, IPHASE, GI ) * VLKNN
                            END DO

                            DO IDIM = 1, NDIM_VEL
                                JDIM = IDIM
                                IF ( .NOT.JUST_BL_DIAG_MAT ) THEN
                                    IF ( NO_MATRIX_STORE ) THEN
                                        LOC_U_RHS( IDIM, IPHASE, U_ILOC ) = LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                        - VLK_UVW( IDIM ) * LOC_U( IDIM, IPHASE, U_JLOC )
                                    ELSE
                                        DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE )  &
                                        = DIAG_BIGM_CON( IDIM, JDIM, IPHASE, JPHASE, U_ILOC, U_JLOC, ELE ) + VLK_UVW( IDIM )
                                    END IF
                                END IF
                            END DO

                        END DO
                    END DO
                END DO

                ! Place the diffusion term into matrix for between element diffusion stabilization...
                IF ( BETWEEN_ELE_STAB ) THEN
                    ! we store these vectors in order to try and work out the between element
                    ! diffusion/viscocity.
                    DO U_ILOC = 1, U_NLOC
                        DO U_JLOC = 1, U_NLOC
                            MAT_ELE( U_ILOC, U_JLOC, ELE ) = MAT_ELE( U_ILOC, U_JLOC, ELE ) + &
                            SUM( UFEN( U_ILOC, : ) * UFEN( U_JLOC,  : ) * DETWEI( : ) )
                        END DO
                    END DO

                    DO U_ILOC = 1, U_NLOC
                        DO IPHASE = 1, NPHASE
                            DO GI = 1, CV_NGI
                                ! we store these vectors in order to try and work out the between element
                                ! diffusion/viscocity.
                                DIFF_FOR_BETWEEN_U( :, IPHASE, U_ILOC, ELE ) = DIFF_FOR_BETWEEN_U( :, IPHASE, U_ILOC, ELE ) &
                                + UFEN( U_ILOC, GI ) * DETWEI( GI ) * DIF_STAB_U( :, IPHASE, GI )
                            END DO
                        END DO
                    END DO
                   ! End of IF(BETWEEN_ELE_STAB) THEN...
                END IF

               !! *************************INNER ELEMENT STABILIZATION****************************************
               !! *************************INNER ELEMENT STABILIZATION****************************************
               ! endof IF(RESID_BASED_STAB_DIF.NE.0) THEN
            END IF
            ! **********REVIEWER 2-END**********************

            ! copy local memory
            DO U_ILOC = 1, U_NLOC
                U_INOD = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_ILOC )
                DO IPHASE = 1, NPHASE
                    DO IDIM = 1, NDIM_VEL
                        I = U_INOD + (IDIM-1)*U_NONODS + (IPHASE-1)*NDIM_VEL*U_NONODS
                        U_RHS( IDIM, IPHASE, U_INOD ) = U_RHS( IDIM, IPHASE, U_INOD ) + LOC_U_RHS( IDIM, IPHASE, U_ILOC )
                    END DO
                END DO
            END DO

        END DO Loop_Elements


        !!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!!
        !!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!!


        !! *************************loop over surfaces*********************************************
        ! at some pt we need to merge these 2 loops but there is a bug when doing that!!!!!

        ! **********REVIEWER 3-START**********************
        DISC_PRES = ( CV_NONODS == TOTELE * CV_NLOC )

        Loop_Elements2: DO ELE = 1, TOTELE

            ! for copy local memory copying...
            LOC_U_RHS = 0.0

            Between_Elements_And_Boundary: DO IFACE = 1, NFACE
                ELE2  = FACE_ELE( IFACE, ELE )
                SELE2 = MAX( 0, - ELE2 )
                SELE  = SELE2
                ELE2  = MAX( 0, + ELE2 )

                ! Find COUNT_ELE
                IF(.NOT.NO_MATRIX_STORE) THEN
                    DO COUNT=FINELE(ELE), FINELE(ELE+1)-1
                        IF(ELE2==COLELE(COUNT)) COUNT_ELE=COUNT
                    END DO
                ENDIF

                ! The surface nodes on element face IFACE.
                U_SLOC2LOC( : ) = U_SLOCLIST( IFACE, : )
                CV_SLOC2LOC( : ) = CV_SLOCLIST( IFACE, : )

                ! Recalculate the normal...
                DO CV_ILOC = 1, CV_NLOC
                    X_INOD = X_NDGLN( (ELE-1)*X_NLOC + CV_ILOC )
                    XL_ALL(:,CV_ILOC) = X_ALL( :, X_INOD )
                END DO

                ! Recalculate the normal...
                DO CV_SILOC = 1, CV_SNLOC
                    CV_ILOC = CV_SLOC2LOC( CV_SILOC )
                    X_INOD = X_NDGLN( (ELE-1)*X_NLOC + CV_ILOC )
                    XSL_ALL( :, CV_SILOC ) = X_ALL( :, X_INOD )
                END DO

                CALL DGSIMPLNORM_ALL( CV_NLOC, CV_SNLOC, NDIM, &
                XL_ALL, XSL_ALL, NORMX_ALL )

                CALL DGSDETNXLOC2_ALL( CV_SNLOC, SBCVNGI, NDIM, &
                XSL_ALL, &
                SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SDETWE, SAREA, &
                SNORMXN_ALL, &
                NORMX_ALL )

                If_ele2_notzero_1: IF(ELE2 /= 0) THEN
                    ! ***********SUBROUTINE DETERMINE_OTHER_SIDE_FACE - START************

                    If_stored: IF(STORED_OTHER_SIDE) THEN

                        U_ILOC_OTHER_SIDE( : ) = STORED_U_ILOC_OTHER_SIDE( :, IFACE, ELE )
                        U_OTHER_LOC( : )       = STORED_U_OTHER_LOC( :, IFACE, ELE )
                        MAT_OTHER_LOC( : )     = STORED_MAT_OTHER_LOC( :, IFACE, ELE )

                    ELSE If_stored

                        IF ( IS_OVERLAPPING ) THEN
                            U_OTHER_LOC=0
                            U_ILOC_OTHER_SIDE=0
                            IF( XU_NLOC == 1 ) THEN ! For constant vel basis functions...
                                DO ILEV = 1, CV_NLOC
                                    U_ILOC_OTHER_SIDE( 1 +(ILEV-1)*U_SNLOC/CV_NLOC) &
                                    = 1 + (ILEV-1)*U_NLOC/CV_NLOC
                                    U_OTHER_LOC( 1 + (ILEV-1)*U_NLOC/CV_NLOC) &
                                    = 1 + (ILEV-1)*U_NLOC/CV_NLOC
                                END DO
                            ELSE
                                DO U_SILOC = 1, U_SNLOC/CV_NLOC
                                    U_ILOC = U_SLOC2LOC( U_SILOC )
                                    U_INOD = XU_NDGLN(( ELE - 1 ) * XU_NLOC + U_ILOC )
                                    DO U_ILOC2 = 1, U_NLOC/CV_NLOC
                                        U_INOD2 = XU_NDGLN(( ELE2 - 1 ) * XU_NLOC + U_ILOC2 )
                                        IF ( U_INOD2 == U_INOD ) THEN
                                            DO ILEV = 1, CV_NLOC
                                                U_ILOC_OTHER_SIDE( U_SILOC +(ILEV-1)*U_SNLOC/CV_NLOC) &
                                                = U_ILOC2 + (ILEV-1)*U_NLOC/CV_NLOC
                                                U_OTHER_LOC( U_ILOC + (ILEV-1)*U_NLOC/CV_NLOC) &
                                                = U_ILOC2 + (ILEV-1)*U_NLOC/CV_NLOC
                                            END DO
                                        END IF
                                    END DO
                                END DO
                            END IF
                        ELSE ! Not overlapping...
                            U_OTHER_LOC=0
                            U_ILOC_OTHER_SIDE=0
                            IF( XU_NLOC == 1 ) THEN ! For constant vel basis functions...
                                U_ILOC_OTHER_SIDE( 1 ) = 1
                                U_OTHER_LOC( 1 )= 1
                            ELSE
                                DO U_SILOC = 1, U_SNLOC
                                    U_ILOC = U_SLOC2LOC( U_SILOC )
                                    U_INOD = XU_NDGLN( ( ELE - 1 ) * U_NLOC + U_ILOC )
                                    DO U_ILOC2 = 1, U_NLOC
                                        U_INOD2 = XU_NDGLN(( ELE2 - 1 ) * U_NLOC + U_ILOC2 )
                                        IF ( U_INOD2 == U_INOD ) THEN
                                            U_ILOC_OTHER_SIDE( U_SILOC ) = U_ILOC2
                                            U_OTHER_LOC( U_ILOC )=U_ILOC2
                                        END IF
                                    END DO
                                END DO
                            ENDIF
                        ENDIF

                        MAT_OTHER_LOC=0
                        DO MAT_SILOC = 1, CV_SNLOC
                            MAT_ILOC = CV_SLOC2LOC( MAT_SILOC )
                            MAT_INOD = X_NDGLN(( ELE - 1 ) * MAT_NLOC + MAT_ILOC )
                            DO MAT_ILOC2 = 1, MAT_NLOC
                                MAT_INOD2 = X_NDGLN(( ELE2 - 1 ) * MAT_NLOC + MAT_ILOC2 )
                                IF ( MAT_INOD2 == MAT_INOD ) THEN
                                    MAT_OTHER_LOC( MAT_ILOC )=MAT_ILOC2
                                END IF
                            END DO
                        END DO

                        IF ( ISTORED_OTHER_SIDE.NE.0 ) THEN

                            STORED_U_ILOC_OTHER_SIDE( :, IFACE, ELE ) = U_ILOC_OTHER_SIDE( : )
                            STORED_U_OTHER_LOC( :, IFACE, ELE )       = U_OTHER_LOC( : )
                            STORED_MAT_OTHER_LOC( :, IFACE, ELE )     = MAT_OTHER_LOC( : )

                        END IF


                    END IF If_stored

                   ! ***********SUBROUTINE DETERMINE_OTHER_SIDE_FACE - END************
                END IF If_ele2_notzero_1



                ! ********Mapping to local variables****************
                ! CV variables...
                DO CV_SILOC = 1, CV_SNLOC
                    CV_ILOC = CV_SLOC2LOC( CV_SILOC )
                    CV_INOD = CV_NDGLN( (ELE-1)*CV_NLOC + CV_ILOC )
                    IF ( ELE2 /= 0) THEN
                        CV_ILOC2 = MAT_OTHER_LOC( CV_ILOC )
                        CV_INOD2 = CV_NDGLN( (ELE2-1)*CV_NLOC + CV_ILOC2 )
                    ELSE
                        CV_ILOC2 = CV_ILOC
                        CV_INOD2 = CV_INOD
                    END IF

                    IF(IGOT_VOL_X_PRESSURE==1) THEN
                       SLOC_UDEN( :, CV_SILOC )  = UDEN( :, CV_INOD ) * FEM_VOL_FRAC( :, CV_INOD )
                       SLOC2_UDEN( :, CV_SILOC ) = UDEN( :, CV_INOD2 ) * FEM_VOL_FRAC( :, CV_INOD2 )
                       SLOC_UDENOLD( :, CV_SILOC ) = UDENOLD( :, CV_INOD ) * FEM_VOL_FRAC( :, CV_INOD )
                       SLOC2_UDENOLD( :, CV_SILOC ) = UDENOLD( :, CV_INOD2 ) * FEM_VOL_FRAC( :, CV_INOD2 )
                    ELSE
                       SLOC_UDEN( :, CV_SILOC )  = UDEN( :, CV_INOD )
                       SLOC2_UDEN( :, CV_SILOC ) = UDEN( :, CV_INOD2 )
                       SLOC_UDENOLD( :, CV_SILOC ) = UDENOLD( :, CV_INOD )
                       SLOC2_UDENOLD( :, CV_SILOC ) = UDENOLD( :, CV_INOD2 )
                    ENDIF

                    IF(IDIVID_BY_VOL_FRAC+IGOT_VOL_X_PRESSURE.GE.1) THEN
                       SLOC_VOL_FRA( :, CV_SILOC )  = FEM_VOL_FRAC( :, CV_INOD )
                       SLOC2_VOL_FRA( :, CV_SILOC ) = FEM_VOL_FRAC( :, CV_INOD2 )
                    ENDIF

                    IF ( GOT_DIFFUS ) THEN
                        DO IPHASE = 1, NPHASE
                           SLOC_UDIFFUSION( 1:NDIM, 1:NDIM, IPHASE, CV_SILOC ) = UDIFFUSION_ALL( 1:NDIM, 1:NDIM, IPHASE, CV_INOD )
                           SLOC2_UDIFFUSION( 1:NDIM, 1:NDIM, IPHASE, CV_SILOC ) = UDIFFUSION_ALL( 1:NDIM, 1:NDIM, IPHASE, CV_INOD2 )
                        END DO
                    END IF
                END DO

                DO U_SILOC = 1, U_SNLOC
                    U_ILOC = U_SLOC2LOC( U_SILOC )
                    U_INOD = U_NDGLN( (ELE-1)*U_NLOC + U_ILOC )
                    IF ( ELE2 /= 0 ) THEN
                        U_ILOC2 = U_ILOC_OTHER_SIDE( U_SILOC )
                        U_INOD2 = U_NDGLN( (ELE2-1)*U_NLOC + U_ILOC2 )
                    ELSE
                        U_ILOC2 = U_ILOC
                        U_INOD2 = U_INOD
                    END IF
                    DO IPHASE = 1, NPHASE
                        IF ( GOT_DIFFUS ) THEN
                            SLOC_DUX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_SILOC ) = DUX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_ILOC, ELE )
                            SLOC_DUOLDX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_SILOC ) = DUOLDX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_ILOC, ELE )
                            IF(ELE2 /= 0) THEN
                                SLOC2_DUX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_SILOC ) = DUX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_ILOC2, ELE2 )
                                SLOC2_DUOLDX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_SILOC ) = DUOLDX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_ILOC2, ELE2 )
                            ELSE
                                SLOC2_DUX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_SILOC ) = DUX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_ILOC, ELE )
                                SLOC2_DUOLDX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_SILOC ) = DUOLDX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_ILOC, ELE )
                            END IF
                        END IF
                    END DO
                END DO

                ! velocity variables...
                DO U_SILOC = 1, U_SNLOC
                    U_ILOC = U_SLOC2LOC( U_SILOC )
                    U_INOD = U_NDGLN( (ELE-1)*U_NLOC + U_ILOC )
                    IF ( ELE2 /= 0 ) THEN
                        U_ILOC2 = U_ILOC_OTHER_SIDE( U_SILOC )
                        U_INOD2 = U_NDGLN( (ELE2-1)*U_NLOC + U_ILOC2 )
                    ELSE
                        U_ILOC2 = U_ILOC
                        U_INOD2 = U_INOD
                    END IF
                    ! for normal calc...
                    DO IPHASE = 1, NPHASE

                        IF ( GOT_DIFFUS ) THEN
                            SLOC_DUX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_SILOC ) = DUX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_ILOC, ELE )
                            SLOC_DUOLDX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_SILOC ) = DUOLDX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_ILOC, ELE )
                            IF(ELE2 /= 0) THEN
                                SLOC2_DUX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_SILOC ) = DUX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_ILOC2, ELE2 )
                                SLOC2_DUOLDX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_SILOC ) = DUOLDX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_ILOC2, ELE2 )
                            ELSE
                                SLOC2_DUX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_SILOC ) = DUX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_ILOC, ELE )
                                SLOC2_DUOLDX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_SILOC ) = DUOLDX_ELE_ALL( 1:NDIM_VEL, 1:NDIM, IPHASE, U_ILOC, ELE )
                            END IF
                        END IF

                        IF ( BETWEEN_ELE_STAB ) THEN
                            ! Calculate stabilization diffusion coefficient...
                            DO IDIM_VEL = 1, NDIM_VEL
                                SLOC_DIFF_FOR_BETWEEN_U( IDIM_VEL, IPHASE, U_SILOC ) = DIFF_FOR_BETWEEN_U( IDIM_VEL, IPHASE, U_ILOC, ELE )
                                IF ( ELE2 /= 0 ) THEN
                                    SLOC2_DIFF_FOR_BETWEEN_U( IDIM_VEL, IPHASE, U_SILOC ) = DIFF_FOR_BETWEEN_U( IDIM_VEL, IPHASE, U_ILOC2, ELE2 )
                                ELSE
                                    SLOC2_DIFF_FOR_BETWEEN_U( IDIM_VEL, IPHASE, U_SILOC ) = DIFF_FOR_BETWEEN_U( IDIM_VEL, IPHASE, U_ILOC, ELE )
                                END IF
                            END DO
                        END IF

                        ! U:
                        DO IDIM = 1, NDIM_VEL
                            SLOC_U( IDIM, IPHASE, U_SILOC ) = U_ALL( IDIM, IPHASE, U_INOD )
                            SLOC_UOLD( IDIM, IPHASE, U_SILOC ) = UOLD_ALL( IDIM, IPHASE, U_INOD )
                            SLOC2_U( IDIM, IPHASE, U_SILOC ) = U_ALL( IDIM, IPHASE, U_INOD2 )
                            SLOC2_UOLD( IDIM, IPHASE, U_SILOC ) = UOLD_ALL( IDIM, IPHASE, U_INOD2 )
                        END DO

                        DO IDIM = 1, NDIM
                            SLOC_NU( IDIM, IPHASE, U_SILOC ) = NU_ALL( IDIM, IPHASE, U_INOD )
                            SLOC_NUOLD( IDIM, IPHASE, U_SILOC ) = NUOLD_ALL( IDIM, IPHASE, U_INOD )
                            SLOC2_NU( IDIM, IPHASE, U_SILOC ) = NU_ALL( IDIM, IPHASE, U_INOD2 )
                            SLOC2_NUOLD( IDIM, IPHASE, U_SILOC ) = NUOLD_ALL( IDIM, IPHASE, U_INOD2 )
                        END DO

                    END DO
                END DO



                If_diffusion_or_momentum1: IF(GOT_DIFFUS .OR. GOT_UDEN) THEN
                    SDEN=0.0
                    SDENOLD=0.0
                    SDEN_KEEP=0.0 ; SDEN2_KEEP=0.0
                    SDENOLD_KEEP=0.0 ; SDENOLD2_KEEP=0.0
                    IF(IDIVID_BY_VOL_FRAC+IGOT_VOL_X_PRESSURE.GE.1) THEN
                       SVOL_FRA =0.0
                       SVOL_FRA2=0.0
                    ENDIF
                    DO CV_SILOC=1,CV_SNLOC
                        DO SGI=1,SBCVNGI
                            DO IPHASE=1, NPHASE
                                SDEN(IPHASE,SGI)=SDEN(IPHASE,SGI) + SBCVFEN(CV_SILOC,SGI) &
                                *0.5*(SLOC_UDEN(IPHASE,CV_SILOC)+SLOC2_UDEN(IPHASE,CV_SILOC)) *WITH_NONLIN
                                SDENOLD(IPHASE,SGI)=SDENOLD(IPHASE,SGI) + SBCVFEN(CV_SILOC,SGI) &
                                *0.5*(SLOC_UDENOLD(IPHASE,CV_SILOC)+SLOC2_UDENOLD(IPHASE,CV_SILOC)) *WITH_NONLIN

                                SDEN_KEEP(IPHASE,SGI)=SDEN_KEEP(IPHASE,SGI) + SBCVFEN(CV_SILOC,SGI) &
                                *SLOC_UDEN(IPHASE,CV_SILOC)*WITH_NONLIN
                                SDEN2_KEEP(IPHASE,SGI)=SDEN2_KEEP(IPHASE,SGI) + SBCVFEN(CV_SILOC,SGI) &
                                *SLOC2_UDEN(IPHASE,CV_SILOC)*WITH_NONLIN

                                SDENOLD_KEEP(IPHASE,SGI)=SDENOLD_KEEP(IPHASE,SGI) + SBCVFEN(CV_SILOC,SGI) &
                                *SLOC_UDENOLD(IPHASE,CV_SILOC)*WITH_NONLIN
                                SDENOLD2_KEEP(IPHASE,SGI)=SDENOLD2_KEEP(IPHASE,SGI) + SBCVFEN(CV_SILOC,SGI) &
                                *SLOC2_UDENOLD(IPHASE,CV_SILOC)*WITH_NONLIN
                                IF(IDIVID_BY_VOL_FRAC+IGOT_VOL_X_PRESSURE.GE.1) THEN
                                   SVOL_FRA(IPHASE,SGI) =SVOL_FRA(IPHASE,SGI) + SBCVFEN(CV_SILOC,SGI) *SLOC_VOL_FRA(IPHASE,CV_SILOC)
                                   SVOL_FRA2(IPHASE,SGI)=SVOL_FRA2(IPHASE,SGI)+ SBCVFEN(CV_SILOC,SGI) *SLOC2_VOL_FRA(IPHASE,CV_SILOC)
                                ENDIF
                            END DO
                        END DO
                    END DO

                    SUD_ALL=0.0
                    SUDOLD_ALL=0.0
                    DO U_SILOC=1,U_SNLOC
                        DO SGI=1,SBCVNGI
                            DO IPHASE=1, NPHASE
                                SUD_ALL(:,IPHASE,SGI)   =SUD_ALL(:,IPHASE,SGI)    + SBUFEN(U_SILOC,SGI)*SLOC_NU(:,IPHASE,U_SILOC)
                                SUDOLD_ALL(:,IPHASE,SGI)=SUDOLD_ALL(:,IPHASE,SGI) + SBUFEN(U_SILOC,SGI)*SLOC_NUOLD(:,IPHASE,U_SILOC)
                            END DO
                        END DO
                    END DO

                    SUD_ALL_KEEP=SUD_ALL

                    SUDOLD_ALL_KEEP=SUDOLD_ALL


                ENDIF If_diffusion_or_momentum1



                If_on_boundary_domain: IF(SELE /= 0) THEN
                    ! ***********SUBROUTINE DETERMINE_SUF_PRES - START************
                    ! Put the surface integrals in for pressure b.c.'s
                    ! that is add into C matrix and U_RHS. (DG velocities)

                    U_NLOC2 = MAX( 1, U_NLOC/CV_NLOC )
                    Loop_ILOC2: DO U_SILOC = 1, U_SNLOC
                        U_ILOC = U_SLOC2LOC( U_SILOC )
                        ILEV = (U_ILOC-1)/U_NLOC2 + 1

                        if( .not. is_overlapping ) ilev = 1

                        if(.not.got_c_matrix) IU_NOD = U_SNDGLN(( SELE - 1 ) * U_SNLOC + U_SILOC )

                        Loop_JLOC2: DO P_SJLOC = 1, P_SNLOC
                            P_JLOC = CV_SLOC2LOC( P_SJLOC )


                            if( ( .not. is_overlapping ) .or. ( p_jloc == ilev ) ) then
                                if(.not.got_c_matrix) JCV_NOD = P_SNDGLN(( SELE - 1 ) * P_SNLOC + P_SJLOC )

                                NMX_ALL = 0.0
                                IF(IGOT_VOL_X_PRESSURE==1) VOL_FRA_NMX_ALL(:,:) = 0.0
                                Loop_GaussPoints2: DO SGI = 1, SBCVNGI
                                    NMX_ALL(:) = NMX_ALL(:) + SNORMXN_ALL( :, SGI ) *SBUFEN( U_SILOC, SGI ) * SBCVFEN( P_SJLOC, SGI ) * SDETWE( SGI )
                                    IF(IGOT_VOL_X_PRESSURE==1) THEN
                                       DO IPHASE = 1, NPHASE
                                          VOL_FRA_NMX_ALL( :, IPHASE ) = VOL_FRA_NMX_ALL( :, IPHASE ) + SVOL_FRA( IPHASE,SGI ) * RNMX_ALL( : )
                                       END DO
                                    ENDIF
                                END DO Loop_GaussPoints2


                                ! Put into matrix

                                ! Find COUNT - position in matrix : FINMCY, COLMCY
                                IF ( .NOT.GOT_C_MATRIX ) THEN
                                    CALL USE_POSINMAT_C_STORE( COUNT, IU_NOD, JCV_NOD, &
                                    U_NONODS, FINDC, COLC, NCOLC, &
                                    IDO_STORE_AC_SPAR_PT, STORED_AC_SPAR_PT, POSINMAT_C_STORE, ELE, U_ILOC, P_JLOC, &
                                    TOTELE, U_NLOC, P_NLOC )
                                END IF

                                Loop_Phase2: DO IPHASE = 1, NPHASE
                                    IF( WIC_P_BC_ALL( SELE ) == WIC_P_BC_DIRICHLET ) THEN

                                        DO IDIM = 1, NDIM_VEL
                                            IF(IGOT_VOL_X_PRESSURE==1) THEN
                                               IF ( .NOT.GOT_C_MATRIX ) THEN
                                                   C( IDIM, IPHASE, COUNT ) = C( IDIM, IPHASE, COUNT ) &
                                                   + VOL_FRA_NMX_ALL( IDIM, IPHASE ) * SELE_OVERLAP_SCALE( P_JLOC )
                                               END IF
                                               LOC_U_RHS( IDIM, IPHASE, U_ILOC) =  LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                               - VOL_FRA_NMX_ALL( IDIM, IPHASE ) * SUF_P_BC_ALL( P_SJLOC + P_SNLOC * ( SELE - 1) ) * SELE_OVERLAP_SCALE( P_JLOC )
                                            ELSE
                                               IF ( .NOT.GOT_C_MATRIX ) THEN
                                                   C( IDIM, IPHASE, COUNT ) = C( IDIM, IPHASE, COUNT ) &
                                                   + NMX_ALL( IDIM ) * SELE_OVERLAP_SCALE( P_JLOC )
                                               END IF
                                               LOC_U_RHS( IDIM, IPHASE, U_ILOC) =  LOC_U_RHS( IDIM, IPHASE, U_ILOC ) &
                                               - NMX_ALL( IDIM ) * SUF_P_BC_ALL( P_SJLOC + P_SNLOC* ( SELE - 1 ) ) * SELE_OVERLAP_SCALE( P_JLOC )
                                            ENDIF
                                        END DO

                                    END IF

                                END DO Loop_Phase2
                            ENDIF



                        END DO Loop_JLOC2

                    END DO Loop_ILOC2
                   ! ***********SUBROUTINE DETERMINE_SUF_PRES - END************
                ENDIF If_on_boundary_domain





                If_ele2_notzero: IF(ELE2 /= 0) THEN

                    got_c_matrix1: if(.not.got_c_matrix) then
                        discontinuous_pres: IF(DISC_PRES) THEN

                            DO P_SJLOC = 1, CV_SNLOC
                                P_JLOC = CV_SLOC2LOC( P_SJLOC )
                                P_JNOD = P_NDGLN(( ELE - 1 ) * P_NLOC + P_JLOC )
                                P_JLOC2 = MAT_OTHER_LOC(P_JLOC)
                                P_JNOD2 = P_NDGLN(( ELE2 - 1 ) * P_NLOC + P_JLOC2 )
                                DO U_SILOC = 1, U_SNLOC
                                    U_ILOC = U_SLOC2LOC( U_SILOC )
                                    U_NLOC2 = MAX(1,U_NLOC/CV_NLOC)
                                    ILEV = (U_ILOC-1)/U_NLOC2 + 1
                                    IF( .NOT. IS_OVERLAPPING ) ILEV = 1

                                    IF( ( .NOT. IS_OVERLAPPING ) .OR. &
                                    (( MAT_OTHER_LOC( ILEV ) /= 0 )) ) THEN
                                        U_INOD = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_ILOC )
                                        VNMX_ALL = 0.0
                                        DO SGI = 1, SBCVNGI
                                            RNN = SDETWE( SGI ) * SBUFEN( U_SILOC, SGI ) * SBCVFEN( P_SJLOC, SGI )
                                            VNMX_ALL = VNMX_ALL + SNORMXN_ALL( :, SGI ) * RNN
                                        END DO

                                        CALL USE_POSINMAT_C_STORE( COUNT, U_INOD, P_JNOD,  &
                                        U_NONODS, FINDC, COLC, NCOLC, &
                                        IDO_STORE_AC_SPAR_PT, STORED_AC_SPAR_PT, POSINMAT_C_STORE, ELE, U_ILOC, P_JLOC, &
                                        TOTELE, U_NLOC, P_NLOC )

                                        CALL USE_POSINMAT_C_STORE_SUF_DG( COUNT2, U_INOD, P_JNOD2,  &
                                        U_NONODS, FINDC, COLC, NCOLC, &
                                        IDO_STORE_AC_SPAR_PT, STORED_AC_SPAR_PT, POSINMAT_C_STORE_SUF_DG, ELE, IFACE, U_SILOC, P_SJLOC,  &
                                        TOTELE, NFACE, U_SNLOC, P_SNLOC )

                                        Loop_Phase5: DO IPHASE = 1, NPHASE
                                            COUNT_PHA  = COUNT  + ( IPHASE - 1 ) * NDIM_VEL * NCOLC
                                            COUNT_PHA2 = COUNT2 + ( IPHASE - 1 ) * NDIM_VEL * NCOLC
                                            ! weight integral according to non-uniform mesh spacing otherwise it will go unstable.
                                            IF( VOL_ELE_INT_PRES ) THEN
                                                ! bias the weighting towards bigger eles - works with 0.25 and 0.1 and not 0.01.
                                                MASSE = MASS_ELE( ELE ) + 0.25 * MASS_ELE( ELE2 )
                                                MASSE2 = MASS_ELE( ELE2 ) + 0.25 * MASS_ELE( ELE )
                                            ELSE ! Simple average (works well with IN_ELE_UPWIND=DG_ELE_UPWIND=2)...
                                                MASSE = 1.0
                                                MASSE2 = 1.0
                                            END IF

                                            ! SELE_OVERLAP_SCALE(P_JNOD) is the scaling needed to convert to overlapping element surfaces.
                                            IF ( .NOT.GOT_C_MATRIX ) THEN
                                                DO IDIM = 1, NDIM_VEL
                                                    C( IDIM, IPHASE, COUNT ) = C( IDIM, IPHASE, COUNT ) &
                                                    + VNMX_ALL( IDIM ) * SELE_OVERLAP_SCALE( P_JLOC ) * MASSE / ( MASSE + MASSE2 )

                                                    C( IDIM, IPHASE, COUNT2 ) = C( IDIM, IPHASE, COUNT2 ) &
                                                    - VNMX_ALL( IDIM ) * SELE_OVERLAP_SCALE( P_JLOC ) * MASSE / ( MASSE + MASSE2 )
                                                END DO
                                            END IF
                                        END DO Loop_Phase5
                                    ENDIF
                                END DO
                            END DO
                           !STOP 383
                        ENDIF discontinuous_pres
                    ENDIF got_c_matrix1

                    If_diffusion_or_momentum2: IF(GOT_DIFFUS .OR. GOT_UDEN) THEN
                        ! Calculate distance between centres of elements HDC
                        DO CV_ILOC = 1, CV_NLOC
                            X_INOD = X_NDGLN( (ELE2-1)*X_NLOC + CV_ILOC )
                            XL2_ALL(:,CV_ILOC) = X_ALL( :, X_INOD )
                        END DO

                        DO IDIM = 1, NDIM
                            C1 ( IDIM ) = SUM( XL_ALL( IDIM, : ) ) / REAL( X_NLOC )
                            C2 ( IDIM ) = SUM( XL2_ALL( IDIM, : ) ) / REAL( X_NLOC )
                        END DO
                        HDC = SQRT( SUM( ( C1 - C2 )**2 ) )

                        SUD2_ALL=0.0
                        SUDOLD2_ALL=0.0
                        DO U_SILOC=1,U_SNLOC
                            DO IPHASE=1, NPHASE
                                DO SGI=1,SBCVNGI
                                    SUD2_ALL(:,IPHASE,SGI)=SUD2_ALL(:,IPHASE,SGI) + SBUFEN(U_SILOC,SGI)*SLOC_NU(:,IPHASE,U_SILOC)
                                    SUDOLD2_ALL(:,IPHASE,SGI)=SUDOLD2_ALL(:,IPHASE,SGI) + SBUFEN(U_SILOC,SGI)*SLOC_NUOLD(:,IPHASE,U_SILOC)
                                END DO
                            END DO
                        END DO

                        SUD2_ALL_KEEP=SUD2_ALL

                        SUDOLD2_ALL_KEEP=SUDOLD2_ALL

                        IF(MOM_CONSERV) THEN
                            SUD_ALL=0.5*(SUD_ALL+SUD2_ALL)
                            SUDOLD_ALL=0.5*(SUDOLD_ALL+SUDOLD2_ALL)
                        ENDIF

                    ENDIF If_diffusion_or_momentum2

                ELSE

                    DO IDIM = 1, NDIM
                       C1 ( IDIM ) = SUM( XL_ALL( IDIM, : ) ) / REAL( X_NLOC )
                       C2 ( IDIM ) = SUM( XSL_ALL( IDIM, : ) ) / REAL( CV_SNLOC )
                    END DO
                    HDC = SQRT( SUM( ( C1 - C2 )**2 ) )

                END IF If_ele2_notzero


                IF(GOT_UDEN) THEN
                    IF(MOM_CONSERV) THEN
                        IF(SELE2 /= 0) THEN
                            SUD2_ALL=0.0
                            SUDOLD2_ALL=0.0
                            DO IPHASE=1, NPHASE
                                IF( WIC_U_BC_ALL( 1, IPHASE, SELE2 ) == WIC_U_BC_DIRICHLET) THEN
                                    DO U_SILOC=1,U_SNLOC
                                        DO SGI=1,SBCVNGI
                                            SUD2_ALL(:,IPHASE,SGI)=SUD2_ALL(:,IPHASE,SGI) + SBUFEN(U_SILOC,SGI) * suf_nu_bc_all( :,iphase,u_siloc + u_SNLOC * ( sele2 - 1 ) )

                                            SUDOLD2_ALL(:,IPHASE,SGI)=SUDOLD2_ALL(:,IPHASE,SGI) + SBUFEN(U_SILOC,SGI) * suf_nu_bc_all( :,iphase,u_siloc + u_SNLOC * ( sele2 - 1 ) )
                                        END DO
                                    END DO

                                    DO SGI=1,SBCVNGI
                                        IF( SUM(SUD_ALL(:,IPHASE,SGI)*SNORMXN_ALL(:,SGI)) < 0.0) THEN
                                            SUD_ALL(:,IPHASE,SGI)=0.5*(SUD_ALL(:,IPHASE,SGI)+SUD2_ALL(:,IPHASE,SGI))
                                        ENDIF
                                        IF( SUM(SUDOLD_ALL(:,IPHASE,SGI)*SNORMXN_ALL(:,SGI)) < 0.0) THEN
                                            SUDOLD_ALL(:,IPHASE,SGI)=0.5*(SUDOLD_ALL(:,IPHASE,SGI)+SUDOLD2_ALL(:,IPHASE,SGI))
                                        ENDIF
                                    END DO
                                ENDIF
                            END DO

                        ENDIF
                    ENDIF
                ENDIF

                If_diffusion_or_momentum3: IF(GOT_DIFFUS .OR. GOT_UDEN) THEN

                    IF(BETWEEN_ELE_STAB) THEN
                        ! Calculate stabilization diffusion coefficient...

                        UDIFF_SUF_STAB=0.0
                        DO U_SILOC = 1, U_SNLOC
                            DO IPHASE=1,NPHASE
                                DO IDIM_VEL=1,NDIM_VEL
                                    DO IDIM=1,NDIM
                                        UDIFF_SUF_STAB(IDIM_VEL,IDIM,IDIM,IPHASE,: ) = UDIFF_SUF_STAB(IDIM_VEL,IDIM,IDIM,IPHASE,: )  &
                                        +SBUFEN(U_SILOC,:)*0.5*(  SLOC_DIFF_FOR_BETWEEN_U(IDIM_VEL, IPHASE, U_SILOC) &
                                        + SLOC2_DIFF_FOR_BETWEEN_U(IDIM_VEL, IPHASE, U_SILOC)  )

                                    END DO
                                END DO
                            END DO
                        END DO
                    END IF

                    DO IPHASE = 1, NPHASE
                        DO SGI = 1, SBCVNGI
                            SNDOTQ(IPHASE,SGI)    = SUM( SUD_ALL(:,IPHASE,SGI)*SNORMXN_ALL(:,SGI) )
                            SNDOTQOLD(IPHASE,SGI) = SUM( SUDOLD_ALL(:,IPHASE,SGI)*SNORMXN_ALL(:,SGI) )
                        END DO
                    END DO

                    SINCOME = 0.5 + 0.5 * SIGN( 1.0, -SNDOTQ )
                    SINCOMEOLD = 0.5 + 0.5 * SIGN( 1.0, -SNDOTQOLD )

                    SNDOTQ_IN  = 0.0
                    SNDOTQ_OUT = 0.0
                    SNDOTQOLD_IN  = 0.0
                    SNDOTQOLD_OUT = 0.0


                    IF( NON_LIN_DGFLUX ) THEN
                        DO IPHASE=1, NPHASE
                            DO SGI=1,SBCVNGI
                                SNDOTQ_KEEP(IPHASE,SGI)   = SUM( SUD_ALL_KEEP(:,IPHASE,SGI)*SNORMXN_ALL(:,SGI) )
                                SNDOTQ2_KEEP(IPHASE,SGI)   =SUM( SUD2_ALL_KEEP(:,IPHASE,SGI)*SNORMXN_ALL(:,SGI)  )

                                SNDOTQOLD_KEEP(IPHASE,SGI)   = SUM( SUDOLD_ALL_KEEP(:,IPHASE,SGI)*SNORMXN_ALL(:,SGI) )
                                SNDOTQOLD2_KEEP(IPHASE,SGI)   =SUM( SUDOLD2_ALL_KEEP(:,IPHASE,SGI)*SNORMXN_ALL(:,SGI)  )
                            END DO
                        END DO



                        IF ( ROE_AVE ) THEN ! perform Roe averaging....
                            do iphase = 1, nphase
                                do sgi = 1, SBCVNGI
                                    !  consider momentum normal to the element only...
                                    ! that is the ( (\rho u_n u_n)_left - (\rho u_n u_n)_right ) / ( (u_n)_left - (u_n)_right )
                                    SNDOTQ_ROE(IPHASE,SGI) =( SDEN_KEEP(IPHASE,SGI) * SNDOTQ_KEEP(IPHASE,SGI)**2 - &
                                    SDEN2_KEEP(IPHASE,SGI) * SNDOTQ2_KEEP(IPHASE,SGI)**2 ) &
                                    / tolfun(  SNDOTQ_KEEP(IPHASE,SGI) -  SNDOTQ2_KEEP(IPHASE,SGI) )

                                    SNDOTQOLD_ROE(IPHASE,SGI) =( SDENOLD_KEEP(IPHASE,SGI) * SNDOTQOLD_KEEP(IPHASE,SGI)**2 - &
                                    SDENOLD2_KEEP(IPHASE,SGI) * SNDOTQOLD2_KEEP(IPHASE,SGI)**2 ) &
                                    / tolfun(  SNDOTQOLD_KEEP(IPHASE,SGI) -  SNDOTQOLD2_KEEP(IPHASE,SGI) )
                                end do
                            end do
                            SINCOME = 0.5 + 0.5 * SIGN( 1.0, -SNDOTQ_ROE )
                            SINCOMEOLD = 0.5 + 0.5 * SIGN( 1.0, -SNDOTQOLD_ROE )
                        END IF


                        ELE3 = ELE2
                        IF ( ELE2==0 ) ELE3 = ELE

                        N_DOT_DU=0.0
                        N_DOT_DU2=0.0
                        N_DOT_DUOLD=0.0
                        N_DOT_DUOLD2=0.0
                        DO U_SILOC = 1, U_SNLOC

                            DO IPHASE = 1, NPHASE

                                do sgi = 1, SBCVNGI

                                    vel_dot(sgi)  =  sum( SUD_ALL(:,IPHASE,SGI) *snormxn_all(:,sgi) )
                                    vel_dot2(sgi) =  sum( SUD2_ALL(:,IPHASE,SGI)*snormxn_all(:,sgi) )

                                    velold_dot(sgi)  = sum( SUDOLD_ALL(:,IPHASE,SGI) *snormxn_all(:,sgi) )
                                    velold_dot2(sgi) = sum( SUDOLD2_ALL(:,IPHASE,SGI) *snormxn_all(:,sgi) )

                                    grad_fact(sgi) = sum( UFENX_ALL(1:NDIM,U_ILOC,1)*snormxn_ALL(:,SGI) )
                                end do

                                N_DOT_DU(iphase,:)  = N_DOT_DU(iphase,:)  + grad_fact(:)*vel_dot(:)
                                N_DOT_DU2(iphase,:) = N_DOT_DU2(iphase,:) + grad_fact(:)*vel_dot2(:)

                                N_DOT_DUOLD(iphase,:) = N_DOT_DUOLD(iphase,:)  + grad_fact(:)*velold_dot(:)
                                N_DOT_DUOLD2(iphase,:) = N_DOT_DUOLD2(iphase,:) + grad_fact(:)*velold_dot2(:)
                            END DO

                        END DO
                    END IF

                    ! Have a surface integral on element boundary...
                    ! Calculate the velocities either side of the element...
                    U_NODJ_SGI_IPHASE_ALL=0.0 ; U_NODI_SGI_IPHASE_ALL=0.0
                    UOLD_NODJ_SGI_IPHASE_ALL=0.0 ; UOLD_NODI_SGI_IPHASE_ALL=0.0

                    DO U_SILOC = 1, U_SNLOC
                        DO SGI=1,SBCVNGI
                            DO IPHASE=1, NPHASE
                                U_NODI_SGI_IPHASE_ALL(:,IPHASE,SGI) = U_NODI_SGI_IPHASE_ALL(:,IPHASE,SGI) + SBUFEN(U_SILOC,SGI) * SLOC_U(:,IPHASE,U_SILOC)
                                U_NODJ_SGI_IPHASE_ALL(:,IPHASE,SGI) = U_NODJ_SGI_IPHASE_ALL(:,IPHASE,SGI) + SBUFEN(U_SILOC,SGI) * SLOC2_U(:,IPHASE,U_SILOC)
                                UOLD_NODI_SGI_IPHASE_ALL(:,IPHASE,SGI) = UOLD_NODI_SGI_IPHASE_ALL(:,IPHASE,SGI) + SBUFEN(U_SILOC,SGI) * SLOC_UOLD(:,IPHASE,U_SILOC)
                                UOLD_NODJ_SGI_IPHASE_ALL(:,IPHASE,SGI) = UOLD_NODJ_SGI_IPHASE_ALL(:,IPHASE,SGI) + SBUFEN(U_SILOC,SGI) * SLOC2_UOLD(:,IPHASE,U_SILOC)
                            END DO
                        END DO
                    END DO





                    ! This sub should be used for stress and tensor viscocity replacing the rest...
                    If_GOT_DIFFUS2: IF(GOT_DIFFUS) THEN
                        CALL DIFFUS_CAL_COEFF_STRESS_OR_TENSOR( DIFF_COEF_DIVDX, &
                        DIFF_COEFOLD_DIVDX, STRESS_FORM, ZERO_OR_TWO_THIRDS, &
                        U_SNLOC, U_NLOC, CV_SNLOC, CV_NLOC, MAT_NLOC, NPHASE, &
                        SBCVFEN,SBCVNGI, NDIM_VEL, NDIM, SLOC_UDIFFUSION, SLOC2_UDIFFUSION, UDIFF_SUF_STAB, &
                        HDC, &
                        U_NODJ_SGI_IPHASE_ALL,    U_NODI_SGI_IPHASE_ALL, &
                        UOLD_NODJ_SGI_IPHASE_ALL, UOLD_NODI_SGI_IPHASE_ALL, &
                        ELE, ELE2, SNORMXN_ALL,  &
                        SLOC_DUX_ELE_ALL, SLOC2_DUX_ELE_ALL,   SLOC_DUOLDX_ELE_ALL, SLOC2_DUOLDX_ELE_ALL,  &
                        SELE, STOTEL, WIC_U_BC_ALL(1,:,: ), WIC_U_BC_DIRICHLET  )
                    ELSE If_GOT_DIFFUS2
                        DIFF_COEF_DIVDX   =0.0
                        DIFF_COEFOLD_DIVDX=0.0
                    END IF If_GOT_DIFFUS2
                    ! *************REVIEWER 3-END*************


                    ! *************REVIEWER 4-START*************

                    SNDOTQ_IN  = 0.0
                    SNDOTQ_OUT = 0.0
                    SNDOTQOLD_IN  = 0.0
                    SNDOTQOLD_OUT = 0.0

                    DO SGI=1,SBCVNGI
                        DO IPHASE=1, NPHASE
                            DO IDIM=1, NDIM_VEL

                                !FTHETA( SGI,IDIM,IPHASE )=0.5 !1.0  - should be 1. as there is no theta set for the internal part of an element.
                                FTHETA( IDIM,IPHASE,SGI )=1.0 ! 0.5

                                ! CENT_RELAX=1.0 (central scheme) =0.0 (upwind scheme).
                                IF( NON_LIN_DGFLUX ) THEN
                                    ! non-linear DG flux - if we have an oscillation use upwinding else use central scheme.
                                    CENT_RELAX = dg_oscilat_detect( SNDOTQ_KEEP(IPHASE,SGI), SNDOTQ2_KEEP(IPHASE,SGI), &
                                    N_DOT_DU(IPHASE,SGI), N_DOT_DU2(IPHASE,SGI), SINCOME(IPHASE,SGI), MASS_ELE(ELE), MASS_ELE(ELE2) )
                                    CENT_RELAX_OLD = dg_oscilat_detect( SNDOTQOLD_KEEP(IPHASE,SGI), SNDOTQOLD2_KEEP(IPHASE,SGI), &
                                    N_DOT_DUOLD(IPHASE,SGI), N_DOT_DUOLD2(IPHASE,SGI), SINCOMEOLD(IPHASE,SGI), MASS_ELE(ELE), MASS_ELE(ELE2) )
                                ELSE
                                    IF( UPWIND_DGFLUX ) THEN
                                        ! Upwind DG flux...
                                        CENT_RELAX    =0.0
                                        CENT_RELAX_OLD=0.0
                                    ELSE
                                        ! Central diff DG flux...
                                        CENT_RELAX    =1.0
                                        CENT_RELAX_OLD=1.0
                                    ENDIF
                                ENDIF
                                ! CENT_RELAX=1.0 (central scheme) =0.0 (upwind scheme).

                                SNDOTQ_IN(IDIM,IPHASE,SGI)    =SNDOTQ_IN(IDIM,IPHASE,SGI)  &
                                +FTHETA(IDIM,IPHASE,SGI)*SDEN(IPHASE,SGI)*SNDOTQ(IPHASE,SGI)  &
                                * (0.5 * CENT_RELAX + SINCOME(IPHASE,SGI)*(1.-CENT_RELAX))
                                SNDOTQ_OUT(IDIM,IPHASE,SGI)   =SNDOTQ_OUT(IDIM,IPHASE,SGI)  &
                                +FTHETA(IDIM,IPHASE,SGI)*SDEN(IPHASE,SGI)*SNDOTQ(IPHASE,SGI) &
                                * (0.5* CENT_RELAX + (1.-SINCOME(IPHASE,SGI))*(1.-CENT_RELAX))

                                SNDOTQOLD_IN(IDIM,IPHASE,SGI) =SNDOTQOLD_IN(IDIM,IPHASE,SGI)  &
                                +(1.-FTHETA(IDIM,IPHASE,SGI))*SDEN(IPHASE,SGI)*SNDOTQOLD(IPHASE,SGI)  &
                                * (0.5* CENT_RELAX_OLD + SINCOMEOLD(IPHASE,SGI)*(1.-CENT_RELAX_OLD))
                                SNDOTQOLD_OUT(IDIM,IPHASE,SGI)=SNDOTQOLD_OUT(IDIM,IPHASE,SGI)  &
                                +(1.-FTHETA(IDIM,IPHASE,SGI))*SDEN(IPHASE,SGI)*SNDOTQOLD(IPHASE,SGI) &
                                * (0.5* CENT_RELAX_OLD + (1.-SINCOMEOLD(IPHASE,SGI))*(1.-CENT_RELAX_OLD))


                            END DO
                        END DO

                    END DO




                    DO U_SILOC=1,U_SNLOC
                        U_ILOC   =U_SLOC2LOC(U_SILOC)
                        DO U_SJLOC=1,U_SNLOC
                            U_JLOC =U_SLOC2LOC(U_SJLOC)
                            IF(SELE2 /= 0) THEN
                                U_JLOC2=U_JLOC
                            ELSE
                                U_JLOC2=U_ILOC_OTHER_SIDE(U_SJLOC)
                            ENDIF



                            ! add diffusion term...
                            DO IPHASE = 1, NPHASE
                                JPHASE = IPHASE
                                DO IDIM = 1, NDIM_VEL
                                    JDIM = IDIM



                              I=IDIM + (IPHASE-1)*NDIM_VEL + (U_ILOC-1)*NDIM_VEL*NPHASE
                              J=JDIM + (JPHASE-1)*NDIM_VEL + (U_JLOC-1)*NDIM_VEL*NPHASE
                              IU_NOD_DIM_PHA = I + (ELE-1)*NDIM_VEL*NPHASE*U_NLOC
                              JU_NOD_DIM_PHA = J + (ELE-1)*NDIM_VEL*NPHASE*U_NLOC

                              J2=JDIM + (JPHASE-1)*NDIM_VEL + (U_JLOC2-1)*NDIM_VEL*NPHASE
                              JU2_NOD_DIM_PHA = J2 + (ELE2-1)*NDIM_VEL*NPHASE*U_NLOC

                              IF(.NOT.NO_MATRIX_STORE) THEN
                                 CALL POSINMAT( COUNT, IU_NOD_DIM_PHA, JU_NOD_DIM_PHA, &
                                   U_NONODS * NPHASE * NDIM_VEL, FINDGM_PHA, COLDGM_PHA, NCOLDGM_PHA )
                                 CALL POSINMAT( COUNT2, IU_NOD_DIM_PHA, JU2_NOD_DIM_PHA, &
                                   U_NONODS * NPHASE * NDIM_VEL, FINDGM_PHA, COLDGM_PHA, NCOLDGM_PHA )
                              ENDIF



                                    VLM=0.0
                                    VLM_NEW=0.0
                                    VLM_OLD=0.0
                                    NN_SNDOTQ_IN    = 0.0
                                    NN_SNDOTQ_OUT   = 0.0
                                    NN_SNDOTQOLD_IN = 0.0
                                    NN_SNDOTQOLD_OUT= 0.0
                                    ! Have a surface integral on element boundary...
                                    DO SGI=1,SBCVNGI

                                        RNN=SDETWE(SGI)*SBUFEN(U_SILOC,SGI)*SBUFEN(U_SJLOC,SGI)

                                        VLM=VLM+RNN

                                        IF(IDIVID_BY_VOL_FRAC==1) THEN ! We are dividing by vol fract.
                                           VLM_NEW = VLM_NEW + FTHETA( IDIM,IPHASE,SGI ) * RNN &
                                           * DIFF_COEF_DIVDX( IDIM,IPHASE,SGI ) * SVOL_FRA2(IPHASE,SGI)
                                           VLM_OLD = VLM_OLD + (1.-FTHETA( IDIM,IPHASE,SGI )) * RNN &
                                           * DIFF_COEFOLD_DIVDX( IDIM,IPHASE,SGI ) * SVOL_FRA2(IPHASE,SGI)
                                        ELSE
                                           VLM_NEW = VLM_NEW + FTHETA( IDIM,IPHASE,SGI ) * RNN &
                                           * DIFF_COEF_DIVDX( IDIM,IPHASE,SGI )
                                           VLM_OLD = VLM_OLD + (1.-FTHETA( IDIM,IPHASE,SGI )) * RNN &
                                           * DIFF_COEFOLD_DIVDX( IDIM,IPHASE,SGI )
                                        ENDIF

                                        NN_SNDOTQ_IN    = NN_SNDOTQ_IN     + SNDOTQ_IN(IDIM,IPHASE,SGI)    *RNN
                                        NN_SNDOTQ_OUT   = NN_SNDOTQ_OUT    + SNDOTQ_OUT(IDIM,IPHASE,SGI)   *RNN
                                        NN_SNDOTQOLD_IN = NN_SNDOTQOLD_IN  + SNDOTQOLD_IN(IDIM,IPHASE,SGI) *RNN
                                        NN_SNDOTQOLD_OUT= NN_SNDOTQOLD_OUT + SNDOTQOLD_OUT(IDIM,IPHASE,SGI)*RNN

                                    END DO


                                    IF(SELE2 == 0) THEN

                                        IF(NO_MATRIX_STORE) THEN
                                            IF(MOM_CONSERV) THEN
                                                LOC_U_RHS( IDIM,IPHASE,U_ILOC ) = LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                - NN_SNDOTQ_OUT*SLOC_U( IDIM,IPHASE,U_SJLOC ) -  NN_SNDOTQ_IN*SLOC2_U(IDIM,IPHASE,U_SJLOC)
                                            ELSE
                                                LOC_U_RHS( IDIM,IPHASE,U_ILOC ) = LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                + NN_SNDOTQ_IN*SLOC_U( IDIM,IPHASE,U_SJLOC ) -  NN_SNDOTQ_IN*SLOC2_U(IDIM,IPHASE,U_SJLOC)
                                            ENDIF
                                            ! viscosity...
                                            LOC_U_RHS( IDIM,IPHASE,U_ILOC ) = LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                            - VLM_NEW*SLOC_U( IDIM,IPHASE,U_SJLOC ) +  VLM_NEW*SLOC2_U(IDIM,IPHASE,U_SJLOC)
                                        ELSE

                                            IF(MOM_CONSERV) THEN
                                         !       DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)  &
                                         !       =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)       + NN_SNDOTQ_OUT
                                         !       BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)  &
                                         !       =BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)     + NN_SNDOTQ_IN
                               DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  + NN_SNDOTQ_OUT
                               DGM_PHA( COUNT2 )  =  DGM_PHA( COUNT2 )  +NN_SNDOTQ_IN
                                            ELSE
                                          !      DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)  &
                                          !      =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)     - NN_SNDOTQ_IN
                                          !      BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)  &
                                          !      =BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)   + NN_SNDOTQ_IN
                               DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  - NN_SNDOTQ_IN
                               DGM_PHA( COUNT2 )  =  DGM_PHA( COUNT2 )  +NN_SNDOTQ_IN
                                            ENDIF
                                            ! viscosity...
                                    !        DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)  &
                                    !        =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)       + VLM_NEW
                                    !        BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)  &
                                    !        =BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC2,COUNT_ELE)     - VLM_NEW
                               DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  + VLM_NEW
                               DGM_PHA( COUNT2 )  =  DGM_PHA( COUNT2 )  - VLM_NEW
                                        ENDIF


                                        RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) = RHS_DIFF_U( IDIM, IPHASE, U_ILOC ) &
                                        - VLM_OLD * SLOC_UOLD( IDIM, IPHASE, U_SJLOC ) + VLM_OLD * SLOC2_UOLD( IDIM, IPHASE, U_SJLOC ) &
                                        - VLM_NEW * SLOC_U( IDIM, IPHASE, U_SJLOC )    + VLM_NEW * SLOC2_U( IDIM, IPHASE, U_SJLOC )

                                        IF(MOM_CONSERV) THEN
                                            LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                            -(+NN_SNDOTQOLD_OUT) * SLOC_UOLD( IDIM,IPHASE,U_SJLOC )   -(+NN_SNDOTQOLD_IN) * SLOC2_UOLD( IDIM,IPHASE,U_SJLOC )
                                        ELSE
                                            LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                            -(-NN_SNDOTQOLD_IN) * SLOC_UOLD( IDIM,IPHASE,U_SJLOC )  -(+NN_SNDOTQOLD_IN) * SLOC2_UOLD( IDIM,IPHASE,U_SJLOC )
                                        ENDIF
                                        ! Viscosity...
                                        LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) -VLM_OLD * SLOC_UOLD( IDIM,IPHASE,U_SJLOC )
                                        LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) +VLM_OLD * SLOC2_UOLD( IDIM,IPHASE,U_SJLOC )

                                    ELSE


                                        IF( WIC_U_BC_ALL( IDIM, IPHASE, SELE2 ) == WIC_U_BC_DIRICHLET ) THEN

                                    !       DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)  &
                                    !       =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) + VLM_NEW

                               DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  + VLM_NEW

                                           LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) -VLM_OLD * SLOC_UOLD( IDIM,IPHASE,U_SJLOC )

                                           LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                           =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) + VLM_NEW * SUF_U_BC_ALL( IDIM,IPHASE,U_SJLOC + U_SNLOC* ( SELE2 - 1 ) )

                                        ELSE IF( (WIC_U_BC_ALL( IDIM, IPHASE, SELE2 ) == WIC_U_BC_ROBIN) .OR. &
                                        (WIC_U_BC_ALL( IDIM, IPHASE, SELE2 ) == WIC_U_BC_DIRI_ADV_AND_ROBIN )) THEN

                                            IF(NO_MATRIX_STORE) THEN
                                                LOC_U_RHS( IDIM,IPHASE,U_ILOC ) = LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                - VLM * SUF_U_BC_ALL( IDIM,IPHASE,U_SJLOC + U_SNLOC* ( SELE2 - 1 ) )*SLOC_U( IDIM,IPHASE,U_SJLOC )

                                            ELSE
                                          !      DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) &
                                          !      =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)+ VLM * SUF_U_BC_ALL( IDIM,IPHASE,U_SJLOC + U_SNLOC * ( SELE2 - 1 ) )
                               DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  + VLM * SUF_U_BC_ALL( IDIM,IPHASE,U_SJLOC + U_SNLOC * ( SELE2 - 1 ) )
                                            ENDIF
                                            LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) - VLM * SUF_U_ROB2_BC_ALL( IDIM,IPHASE,U_SJLOC + U_SNLOC * ( SELE2 - 1 ) )

                                            RHS_DIFF_U( IDIM,IPHASE,U_ILOC ) = RHS_DIFF_U( IDIM,IPHASE,U_ILOC ) &
                                            - VLM * SUF_U_BC_ALL( IDIM, IPHASE, U_SJLOC + U_SNLOC * (  SELE2 - 1 ) ) * SLOC_U( IDIM, IPHASE, U_SJLOC ) &
                                            - VLM * SUF_U_ROB2_BC_ALL( IDIM, IPHASE, U_SJLOC + U_SNLOC * (  SELE2 - 1 ) )

                                        ENDIF
                                        ! BC for incoming momentum...
                                        IF( WIC_MOMU_BC_ALL( IDIM, IPHASE, SELE2 ) == WIC_U_BC_DIRICHLET ) THEN

                                            IF(MOM_CONSERV) THEN

                                                IF(.NOT.NO_MATRIX_STORE) THEN
                                                    DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) &
                                                    =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE)+ NN_SNDOTQ_OUT
                                                ELSE
                                                    LOC_U_RHS( IDIM,IPHASE,U_ILOC ) = LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                    - NN_SNDOTQ_OUT * SLOC_U( IDIM,IPHASE,U_SJLOC )
                                                ENDIF
                                                LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN )*SUF_MOMU_BC_ALL( IDIM,IPHASE,U_SJLOC + U_SNLOC * ( SELE2 - 1 ) ) &
                                                - NN_SNDOTQOLD_OUT * SLOC_UOLD(IDIM,IPHASE,U_SJLOC)

                                               ! ENDOF IF(MOM_CONSERV) THEN...
                                            ELSE

                                                IF(.NOT.NO_MATRIX_STORE) THEN
                                       !             DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) &
                                       !             =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) - NN_SNDOTQ_IN
                               DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  - NN_SNDOTQ_IN
                                                ELSE
                                                    LOC_U_RHS( IDIM,IPHASE,U_ILOC ) = LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                    + NN_SNDOTQ_IN * SLOC_U( IDIM,IPHASE,U_SJLOC )
                                                ENDIF
                                                LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN )*SUF_MOMU_BC_ALL( IDIM,IPHASE,U_SJLOC + U_SNLOC * ( SELE2 - 1 ) ) &
                                                + NN_SNDOTQOLD_IN * SLOC_UOLD(IDIM,IPHASE,U_SJLOC)

                                               ! END OF IF(MOM_CONSERV) THEN ELSE...
                                            ENDIF

                                           ! BC for incoming and outgoing momentum (NO leaking of momentum into or out of domain for example)...
                                        ELSE IF( WIC_MOMU_BC_ALL( IDIM, IPHASE, SELE2 ) == WIC_U_BC_DIRICHLET_INOUT ) THEN

                                            IF(MOM_CONSERV) THEN

                                                LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN + NN_SNDOTQ_OUT + NN_SNDOTQOLD_OUT)*SUF_MOMU_BC_ALL( IDIM,IPHASE,U_SJLOC + U_SNLOC * ( SELE2 - 1 ) )

                                               ! ENDOF IF(MOM_CONSERV) THEN...
                                            ELSE

                                                IF(.NOT.NO_MATRIX_STORE) THEN
                                      !              DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) &
                                      !              =DIAG_BIGM_CON(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC,ELE) - (NN_SNDOTQ_IN + NN_SNDOTQ_OUT)
                               DGM_PHA( COUNT )  =  DGM_PHA( COUNT )  - (NN_SNDOTQ_IN + NN_SNDOTQ_OUT)
                                                ELSE
                                                    LOC_U_RHS( IDIM,IPHASE,U_ILOC ) = LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                    + (NN_SNDOTQ_IN + NN_SNDOTQ_OUT) * SLOC_U( IDIM,IPHASE,U_SJLOC )
                                                ENDIF

                                                LOC_U_RHS( IDIM,IPHASE,U_ILOC ) =  LOC_U_RHS( IDIM,IPHASE,U_ILOC ) &
                                                - ( NN_SNDOTQ_IN + NN_SNDOTQOLD_IN + NN_SNDOTQ_OUT + NN_SNDOTQOLD_OUT)*SUF_MOMU_BC_ALL( IDIM,IPHASE,U_SJLOC + U_SNLOC * ( SELE2 - 1 ) ) &
                                                + (NN_SNDOTQOLD_IN + NN_SNDOTQOLD_OUT) * SLOC_UOLD(IDIM,IPHASE,U_SJLOC)

                                               ! END OF IF(MOM_CONSERV) THEN ELSE...
                                            ENDIF
                                           ! END OF IF( WIC_MOMU_BC(SELE2+(IPHASE-1)*STOTEL) == WIC_U_BC_DIRICHLET) THEN ELSE...
                                        ENDIF

                                    ENDIF

                                END DO
                            END DO
                        END DO
                    END DO
                ENDIF If_diffusion_or_momentum3

            END DO Between_Elements_And_Boundary

            !      END DO Loop_Elements2
            !! *************************end loop over surfaces*********************************************




            ! ideally insert inner element stabilization here...

            ! copy local memory
            DO U_ILOC = 1, U_NLOC
                U_INOD = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_ILOC )
                DO IPHASE = 1, NPHASE
                    DO IDIM = 1, NDIM_VEL
                        I = U_INOD + (IDIM-1)*U_NONODS + (IPHASE-1)*NDIM_VEL*U_NONODS
                        U_RHS( IDIM, IPHASE, U_INOD ) = U_RHS( IDIM, IPHASE, U_INOD ) + LOC_U_RHS( IDIM, IPHASE, U_ILOC )
                    END DO
                END DO
            END DO



           !      END DO Loop_Elements
        END DO Loop_Elements2
        ! **********REVIEWER 4-END**********************


        ! This subroutine combines the distributed and block diagonal for an element
        ! into the matrix DGM_PHA.
        IF(.NOT.NO_MATRIX_STORE) THEN
            CALL COMB_VEL_MATRIX_DIAG_DIST(DIAG_BIGM_CON, BIGM_CON, &
            DGM_PHA, NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, & ! Force balance sparsity
            NCOLELE, FINELE, COLELE,  NDIM_VEL, NPHASE, U_NLOC, U_NONODS, TOTELE )  ! Element connectivity.
            DEALLOCATE( DIAG_BIGM_CON )
            DEALLOCATE( BIGM_CON)
        ENDIF

        !If C was not stored in state, after its calculation we store it.
        if (.not.got_c_matrix) then
            Point_C_Mat = C
        end if


        !ewrite(3,*)'p=',p
        !ewrite(3,*)'U_RHS:',U_RHS
        !stop 222
        !do i=1, ndim*nphase*u_nonods
        !   ewrite(3,*) i, sum(DGM_PHA(FINDGM_PHA(i):FINDGM_PHA(i+1)-1))
        !end do

        !EWRITE(3,*)'-STOTEL, U_SNLOC, P_SNLOC:', STOTEL, U_SNLOC, P_SNLOC
        !EWRITE(3,*)'-WIC_P_BC:', WIC_P_BC( 1 : STOTEL * NPHASE )
        !EWRITE(3,*)'-SUF_P_BC:', SUF_P_BC( 1 : STOTEL * P_SNLOC * NPHASE )
        !ewrite(3,*)'pqp'
        !stop 242

        !do i=1,ncolc
        !  ewrite(3,*)'i,c:',i,c(i)
        !end do
        !ewrite(3,*)'U_RHS:',u_rhs
        !ewrite(3,*)'PIVIT_MAT:', PIVIT_MAT
        !ewrite(3,*)'JUST_BL_DIAG_MAT:',JUST_BL_DIAG_MAT
        !stop 27


        DEALLOCATE( UD, UD_ND )
        DEALLOCATE( UDOLD, UDOLD_ND )
        DEALLOCATE( DENGI )
        DEALLOCATE( DENGIOLD )
        DEALLOCATE( GRAD_SOU_GI )

        DEALLOCATE( SIGMAGI )
        DEALLOCATE( SIGMAGI_STAB )
        DEALLOCATE( MAT_M )
        DEALLOCATE( SNORMXN )
        DEALLOCATE( SNORMYN )
        DEALLOCATE( SNORMZN )

        DEALLOCATE( NN_SIGMAGI_ELE )
        DEALLOCATE( NN_SIGMAGI_STAB_ELE )
        DEALLOCATE( NN_MASS_ELE )
        DEALLOCATE( NN_MASSOLD_ELE )

        DEALLOCATE( NXUDN )

        DEALLOCATE( CV_SLOC2LOC )
        DEALLOCATE( U_SLOC2LOC )

        DEALLOCATE(SDETWE)

        DEALLOCATE( U_ILOC_OTHER_SIDE )

        DEALLOCATE( CV_ON_FACE )
        DEALLOCATE( U_ON_FACE )
        DEALLOCATE(CVFEM_ON_FACE)
        DEALLOCATE(UFEM_ON_FACE)

        DEALLOCATE( U_OTHER_LOC )
        DEALLOCATE( MAT_OTHER_LOC )

        DEALLOCATE( TEN_XX )

        DEALLOCATE( VLN )
        DEALLOCATE( VLN_OLD )
        DEALLOCATE( VLK )

        DEALLOCATE( STRESS_IJ_ELE )
        DEALLOCATE( VLK_ELE )

        DEALLOCATE( SUD_ALL )
        DEALLOCATE( SUDOLD_ALL )
        DEALLOCATE( SUD2_ALL )
        DEALLOCATE( SUDOLD2_ALL )
        DEALLOCATE( SNDOTQ )
        DEALLOCATE( SNDOTQOLD )
        DEALLOCATE( SNDOTQ_ROE )
        DEALLOCATE( SNDOTQOLD_ROE )

        DEALLOCATE( SINCOME )
        DEALLOCATE( SINCOMEOLD )
        DEALLOCATE( SDEN )
        DEALLOCATE( SDENOLD )

        DEALLOCATE( SDEN_KEEP )
        DEALLOCATE( SDENOLD_KEEP )
        DEALLOCATE( SDEN2_KEEP )
        DEALLOCATE( SDENOLD2_KEEP )

        DEALLOCATE( DIFF_COEF_DIVDX )
        DEALLOCATE( DIFF_COEFOLD_DIVDX )
        DEALLOCATE( FTHETA )
        DEALLOCATE( SNDOTQ_IN )
        DEALLOCATE( SNDOTQ_OUT )
        DEALLOCATE( SNDOTQOLD_IN )
        DEALLOCATE( SNDOTQOLD_OUT )

        DEALLOCATE( XSL )
        DEALLOCATE( YSL )
        DEALLOCATE( ZSL )

        DEALLOCATE( GRAD_SOU_GI_NMX )

        DEALLOCATE( MASS_ELE )
        DEALLOCATE( FACE_ELE )

        ! Deallocating for non-linear Petrov-Galerkin diffusion stabilization...
        DEALLOCATE( LOC_MASS_INV )
        DEALLOCATE( LOC_MASS )
        DEALLOCATE( RHS_DIFF_U )


        DEALLOCATE( DIFF_VEC_U )


        DEALLOCATE( DIFFGI_U )

        DEALLOCATE( U_DT )

        DEALLOCATE( SOUGI_X )
        DEALLOCATE( RESID_U)
        DEALLOCATE( P_DX )

        DEALLOCATE( U_GRAD_NORM2, U_GRAD_NORM )

        DEALLOCATE( A_DOT_U )
        DEALLOCATE( STAR_U_COEF )
        DEALLOCATE( P_STAR_U )
        DEALLOCATE( DIF_STAB_U )

        DEALLOCATE( UDIFF_SUF_STAB )

        DEALLOCATE( VLK_UVW )


        call deallocate(velocity_BCs)
        call deallocate(velocity_BCs_robin2)
        call deallocate(momentum_BCs)
        call deallocate(pressure_BCs)


        ewrite(3,*)'Leaving assemb_force_cty'

        RETURN

    END SUBROUTINE ASSEMB_FORCE_CTY





            SUBROUTINE VISCOCITY_TENSOR_LES_CALC(LES_UDIFFUSION, DUX_ELE_ALL, &
                                                 NDIM,NPHASE, U_NLOC,X_NLOC,TOTELE, X_NONODS, &
                                                 X_ALL, X_NDGLN,  MAT_NONODS, MAT_NLOC, MAT_NDGLN, LES_DISOPT)
! This subroutine calculates a tensor of viscocity LES_UDIFFUSION. 
            IMPLICIT NONE
            INTEGER, intent( in ) :: NDIM, NPHASE, U_NLOC, X_NLOC, TOTELE, X_NONODS, MAT_NONODS, MAT_NLOC, LES_DISOPT
            INTEGER, DIMENSION( X_NLOC * TOTELE  ), intent( in ) :: X_NDGLN
            INTEGER, DIMENSION( MAT_NLOC * TOTELE  ), intent( in ) :: MAT_NDGLN
            REAL, DIMENSION( NDIM, X_NONODS  ), intent( in ) :: X_ALL
            REAL, DIMENSION( NDIM, NDIM, NPHASE, MAT_NONODS  ), intent( inout ) :: LES_UDIFFUSION
            REAL, DIMENSION( NDIM, NDIM, NPHASE, U_NLOC, TOTELE  ), intent( in ) :: DUX_ELE_ALL
! Local variables...
!            INTEGER, PARAMETER :: LES_DISOPT=1
! LES_DISOPT is LES option e.g. =1 Anisotropic element length scale
!                               =2 Take the average length scale h
!                               =3 Take the min length scale h
!                               =4 Take the max length scale h
            integer :: ele, MAT_iloc, MAT_INOD
            real, dimension( :, :, :, :, : ), allocatable :: LES_U_UDIFFUSION, LES_MAT_UDIFFUSION
            integer, dimension( : ), allocatable :: NOD_COUNT

            ALLOCATE(LES_U_UDIFFUSION(NDIM,NDIM,NPHASE,U_NLOC,TOTELE))
            ALLOCATE(LES_MAT_UDIFFUSION(NDIM,NDIM,NPHASE,MAT_NLOC,TOTELE))
            ALLOCATE(NOD_COUNT(MAT_NONODS))

            CALL VISCOCITY_TENSOR_LES_CALC_U(LES_U_UDIFFUSION, DUX_ELE_ALL, NDIM,NPHASE, U_NLOC,X_NLOC,TOTELE, X_NONODS, &
                                                 X_ALL, X_NDGLN, LES_DISOPT)

            IF(MAT_NLOC==U_NLOC) THEN
               LES_MAT_UDIFFUSION=LES_U_UDIFFUSION
            ELSE IF( (U_NLOC==3.AND.MAT_NLOC==6) .OR. (U_NLOC==4.AND.MAT_NLOC==10) ) THEN
! 
               LES_MAT_UDIFFUSION(:,:,:,1,:) = LES_U_UDIFFUSION(:,:,:,1,:)
               LES_MAT_UDIFFUSION(:,:,:,3,:) = LES_U_UDIFFUSION(:,:,:,2,:)
               LES_MAT_UDIFFUSION(:,:,:,6,:) = LES_U_UDIFFUSION(:,:,:,3,:)

               LES_MAT_UDIFFUSION(:,:,:,2,:) = 0.5 * ( LES_U_UDIFFUSION(:,:,:,1,:) + LES_U_UDIFFUSION(:,:,:,2,:) )
               LES_MAT_UDIFFUSION(:,:,:,4,:) = 0.5 * ( LES_U_UDIFFUSION(:,:,:,1,:) + LES_U_UDIFFUSION(:,:,:,3,:) )
               LES_MAT_UDIFFUSION(:,:,:,5,:) = 0.5 * ( LES_U_UDIFFUSION(:,:,:,2,:) + LES_U_UDIFFUSION(:,:,:,3,:) )

               if( MAT_NLOC == 10 ) then
                  LES_MAT_UDIFFUSION(:,:,:,7,:) = 0.5 * ( LES_U_UDIFFUSION(:,:,:,1,:) + LES_U_UDIFFUSION(:,:,:,10,:) )
                  LES_MAT_UDIFFUSION(:,:,:,8,:) = 0.5 * ( LES_U_UDIFFUSION(:,:,:,3,:) + LES_U_UDIFFUSION(:,:,:,10,:) )
                  LES_MAT_UDIFFUSION(:,:,:,9,:) = 0.5 * ( LES_U_UDIFFUSION(:,:,:,6,:) + LES_U_UDIFFUSION(:,:,:,10,:) )

                  LES_MAT_UDIFFUSION(:,:,:,10,:) = LES_U_UDIFFUSION(:,:,:,4,:)
               end if

            ELSE IF( (U_NLOC==6.AND.MAT_NLOC==3) .OR. (U_NLOC==10.AND.MAT_NLOC==4) ) THEN
               LES_MAT_UDIFFUSION(:,:,:,1,:) = LES_U_UDIFFUSION(:,:,:,1,:)
               LES_MAT_UDIFFUSION(:,:,:,2,:) = LES_U_UDIFFUSION(:,:,:,3,:)
               LES_MAT_UDIFFUSION(:,:,:,3,:) = LES_U_UDIFFUSION(:,:,:,6,:)
               if( MAT_nloc == 10 ) then
                  LES_MAT_UDIFFUSION(:,:,:,4,:) = LES_U_UDIFFUSION(:,:,:,10,:)
               end if
            ELSE
              PRINT *,'not ready to onvert between these elements'
              STOP 2211
            ENDIF


! Now map to nodal variables from element variables...
         NOD_COUNT=0
         LES_UDIFFUSION=0.0
         do ele = 1, totele
            do MAT_iloc = 1, MAT_nloc
               MAT_Inod = MAT_ndgln( ( ele - 1 ) * MAT_nloc + MAT_iloc )
               LES_UDIFFUSION(:,:,:,MAT_Inod) = LES_UDIFFUSION(:,:,:,MAT_Inod) + LES_MAT_UDIFFUSION(:,:,:,MAT_ILOC,ELE)
               NOD_COUNT(MAT_INOD) = NOD_COUNT(MAT_INOD) + 1
            END DO
         END DO

         DO MAT_INOD=1,MAT_NONODS
               LES_UDIFFUSION(:,:,:,MAT_INOD) = LES_UDIFFUSION(:,:,:,MAT_INOD)/REAL( NOD_COUNT(MAT_INOD) )
         END DO

         RETURN
         END SUBROUTINE VISCOCITY_TENSOR_LES_CALC




            SUBROUTINE VISCOCITY_TENSOR_LES_CALC_U(LES_U_UDIFFUSION, DUX_ELE_ALL, NDIM,NPHASE, U_NLOC,X_NLOC,TOTELE, X_NONODS, &
                                                 X_ALL, X_NDGLN, LES_DISOPT)
! This subroutine calculates a tensor of viscocity LES_UDIFFUSION. 
            IMPLICIT NONE
            INTEGER, intent( in ) :: NDIM, NPHASE, U_NLOC, X_NLOC, TOTELE, X_NONODS, LES_DISOPT
! LES_DISOPT is LES option e.g. =1 Anisotropic element length scale
            INTEGER, DIMENSION( X_NLOC * TOTELE  ), intent( in ) :: X_NDGLN
            REAL, DIMENSION( NDIM, X_NONODS  ), intent( in ) :: X_ALL
            REAL, DIMENSION( NDIM, NDIM, NPHASE, U_NLOC, TOTELE  ), intent( inout ) :: LES_U_UDIFFUSION
            REAL, DIMENSION( NDIM, NDIM, NPHASE, U_NLOC, TOTELE  ), intent( in ) :: DUX_ELE_ALL
! Local variables...
            LOGICAL, PARAMETER :: ONE_OVER_H2=.TRUE.
            !     SET to metric which has 1/h^2 in it
            REAL, PARAMETER :: CS=0.1
            REAL :: LOC_X_ALL(NDIM, X_NLOC), TENSXX_ALL(NDIM, NDIM), RSUM, FOURCS, CS2, VIS 
            INTEGER :: ELE, X_ILOC, U_ILOC, IPHASE, X_NODI, IDIM, JDIM, KDIM

            CS2=CS**2
            FOURCS=4.*CS2

            DO ELE=1,TOTELE

               DO X_ILOC=1,X_NLOC
                  X_NODI = X_NDGLN((ELE-1)*X_NLOC+X_ILOC) 
                  LOC_X_ALL(:,X_ILOC) = X_ALL(:,X_NODI)
               END DO

               CALL ONEELETENS_ALL( LOC_X_ALL, LES_DISOPT, ONE_OVER_H2, TENSXX_ALL, X_NLOC, NDIM )

               DO U_ILOC=1,U_NLOC
                  DO IPHASE=1,NPHASE 

                     RSUM=0.0
                     DO IDIM=1,NDIM
                        DO JDIM=1,NDIM
                           RSUM=RSUM + (0.5*( DUX_ELE_ALL(IDIM,JDIM,IPHASE,U_ILOC,ELE) + DUX_ELE_ALL(JDIM,IDIM,IPHASE,U_ILOC,ELE) ))**2
                        END DO
                     END DO
                     RSUM=SQRT(RSUM) 
                     VIS=RSUM 

! THEN FIND TURBULENT 'VISCOSITIES'

               ! Put a bit in here which multiplies E by FOURCS*VIS 
                     LES_U_UDIFFUSION(:,:,IPHASE,U_ILOC,ELE)= FOURCS*VIS*TENSXX_ALL(:,:)

                  END DO ! DO IPHASE=1,NPHASE
               END DO ! DO U_ILOC=1,U_NLOC

            END DO
            RETURN
            END SUBROUTINE VISCOCITY_TENSOR_LES_CALC_U





       SUBROUTINE ONEELETENS_ALL( LOC_X_ALL, LES_DISOPT, ONE_OVER_H2, TENSXX_ALL, X_NLOC, NDIM )
         !     This sub calculates the ELEMENT-WISE TENSOR TENS
         !     REPRESENTS THE SIZE AND SHAPE OF THE SURROUNDING ELEMENTS.
         !     LES_DISOPT=LES option.
         IMPLICIT NONE
         INTEGER, intent( in ) ::  X_NLOC, NDIM
         LOGICAL, intent( in ) ::  ONE_OVER_H2
         INTEGER, intent( in ) ::  LES_DISOPT

         REAL, intent( inout ) ::  TENSXX_ALL(NDIM,NDIM) 
         REAL, intent( in ) ::  LOC_X_ALL(NDIM,X_NLOC)

         !     HX,HY-characteristic length scales in x,y directions.
         !     Local variables...
         ! IF ONE_OVER_H2=.TRUE. then SET to metric which has 1/h^2 in it
         REAL RN
         REAL AA(NDIM,NDIM),V(NDIM,NDIM),D(NDIM),A(NDIM,NDIM)

         REAL UDL_ALL(NDIM, X_NLOC*X_NLOC)
         REAL GAMMA(X_NLOC*X_NLOC)

         INTEGER ELE,ILOC,L,L1,L2,IGLX1,IGLX2,IGLX,ID,NID,IDIM,JDIM,KDIM

         REAL HOVERQ
         REAL RWIND, D_SCALAR
         REAL RFACT,RT1,RT2,RT3,D1,D2,D3,VOLUME

         RWIND =1./REAL(6)
         NID=X_NLOC*X_NLOC

         TENSXX_ALL=0.0

         !     This subroutine forms a contabution to the Right Hand Side
         !     of Poissons pressure equation, as well as  F1 & F2.


            !     C The first is the old filter term, the second the new one MDP getting
            !     c different results and stabiltiy for tidal applications ????
            RWIND =1./REAL(6)
            NID=X_NLOC*X_NLOC
            !     **********calculate normalised velocitys across element...  
            ID=0
            do L1=1,X_NLOC
               do L2=1,X_NLOC
                  ID=ID+1
                  if(l1.eq.l2) then
                     UDL_ALL(:,ID)=0.0
                     GAMMA(ID)=0.0
                  else
                     UDL_ALL(:,ID)=LOC_X_ALL(:,L1)-LOC_X_ALL(:,L2)

                     !     Normalise 
                     RN=SQRT( SUM(UDL_ALL(:,ID)**2) )
                     UDL_ALL(:,ID)=UDL_ALL(:,ID)/RN
                     !     HX,HY are the characteristic length scales in x,y directions. 
                     HOVERQ=RN
                     GAMMA(ID)=RWIND*HOVERQ 
                  endif
               END DO
            END DO
            !     **********calculate normalised velocitys across element... 


         do  ID=1,NID

               RFACT=GAMMA(ID)/REAL(X_NLOC) 

               DO IDIM=1,NDIM
                  DO JDIM=1,NDIM
                     TENSXX_ALL(IDIM,JDIM)=TENSXX_ALL(IDIM,JDIM) + RFACT*UDL_ALL(IDIM,ID)*UDL_ALL(JDIM,ID)
                  END DO
               END DO

               !     USE THE COMPONENT OF DIFLIN THE X,Y & Z-DIRECTIONS 
               !     RESPECTIVELY FOR C1T,C2T,C3T.
         end do

         !     nb we want 1/L^2 - at the moment we have L on the diagonal.
         !     Make sure the eigen-values are positive...
         AA=TENSXX_ALL

         CALL JACDIA(AA,V,D,NDIM,A,.FALSE.)

         IF(LES_DISOPT==1) THEN ! Take the anisotropic length scales
            D(:)=D(:)
         ELSE IF(LES_DISOPT==2) THEN ! Take the average length scale h
            D_SCALAR=SUM(D(:))/REAL(NDIM)
            D(:)=D_SCALAR
         ELSE IF(LES_DISOPT==3) THEN ! Take the min length scale h
            D_SCALAR=MINVAL(D(:))
            D(:)=D_SCALAR
         ELSE IF(LES_DISOPT==4) THEN ! Take the max length scale h
            D_SCALAR=MAXVAL(D(:))
            D(:)=D_SCALAR
         ELSE
            !            ERROR("NOT A VALID OPTION FOR LES ASSEMBLED EQNS")
            STOP 9331
         ENDIF

         IF(ONE_OVER_H2) THEN
            !     SET to metric which has 1/h^2 in it...
            D(:)=1./MAX(1.E-16,D(:)**2)
         ELSE 
            !     set to inverse of metric which is a multiple of the tensor
            D(:)=MAX(1.E-16,D(:)**2)
         ENDIF

         TENSXX_ALL=0.0
         DO IDIM=1,NDIM
            DO JDIM=1,NDIM

               DO KDIM=1,NDIM
! TENSOR=V^T D V
                     TENSXX_ALL(IDIM,JDIM)=TENSXX_ALL(IDIM,JDIM) + V(KDIM,IDIM) * D(KDIM) * V(KDIM,JDIM)
               END DO

            END DO
         END DO

       RETURN
       END SUBROUTINE ONEELETENS_ALL


!
!
          SUBROUTINE JACDIA(AA,V,D,N, &
! Working arrays...
     &       A,PRISCR) 
! This sub performs Jacobi rotations of a symmetric matrix in order to 
! find the eigen-vectors V and the eigen values A so 
! that AA=V^T D V & D is diagonal. 
! It uses the algorithm of Matrix Computations 2nd edition, p196. 
          IMPLICIT NONE
          REAL TOLER,CONVEG
          PARAMETER(TOLER=1.E-14,CONVEG=1.E-7) 
          INTEGER N
          REAL AA(N,N),V(N,N),D(N), A(N,N)
          LOGICAL PRISCR
! Local variables...
          REAL R,ABSA,MAXA,COSAL2,COSALF,SINAL2,SINALF,MAXEIG
          INTEGER ITS,NITS,Q,P,QQ,PP 
! 
          NITS=9*(N*N-N)/2
!
!          CALL RCLEAR(V,N*N)
!          CALL TOCOPY(A,AA,N*N) 
      do P=1,N
      do Q=1,N
            V(P,Q)=0.
            A(P,Q)=AA(P,Q)
!             ewrite(2,*) 'P,Q,AA:',P,Q,AA(P,Q) 
          END DO
          END DO
!
!
!
!     Check first whether matrix is diagonal
            IF(Q.EQ.0) THEN

      do PP=1,N
                  D(PP) = A(PP,PP)
      do QQ=1,N
                     IF(PP.EQ.QQ) THEN
                        V(PP,QQ) = 1.0
                     ELSE
                        V(PP,QQ) = 0.0
                     END IF
                  END DO
               END DO
               RETURN
            END IF
!
!
!            
          MAXEIG=0.
      do P=1,N
            V(P,P)=1.0
            MAXEIG=MAX(MAXEIG,ABS(A(P,P)))
          END DO
          IF(MAXEIG.LT.TOLER) THEN
            D(1:N) = 0.0
            GOTO 2000
          ENDIF
!           ewrite(2,*) 'maxeig=',maxeig
!
      do  ITS=1,NITS! Was loop 10
! Find maximum on upper diagonal of matrix. 
! QQ is the coln; PP is the row. 
            Q=0
            P=0
            MAXA=0.
      do PP=1,N-1
      do QQ=PP+1,N
                  ABSA=ABS(A(PP,QQ)) 
                  IF(ABSA.GT.MAXA) THEN
                     MAXA=ABSA
                     Q=QQ
                     P=PP
                  ENDIF
               END DO
            END DO

!            IF(PRISCR) ewrite(2,*) 'MAXA,MAXEIG,its=',MAXA,MAXEIG,its
            IF(MAXA/MAXEIG.LT.CONVEG) GOTO 2000
! Rotate with (Q,P) postions.
            R=MAX(TOLER,SQRT( (A(P,P)-A(Q,Q))**2 + 4.*A(P,Q)**2 ) ) 
            IF(A(P,P).GT.A(Q,Q)) THEN
              COSAL2=0.5+0.5*(A(P,P)-A(Q,Q))/R
              COSALF=SQRT(COSAL2)
              IF(ABS(COSALF).LT.TOLER) COSALF=TOLER
              SINALF=A(Q,P)/(R*COSALF)
            ELSE
              SINAL2=0.5-0.5*(A(P,P)-A(Q,Q))/R
              SINALF=SQRT(SINAL2)
              IF(ABS(SINALF).LT.TOLER) SINALF=TOLER
              COSALF=A(Q,P)/(R*SINALF)
            ENDIF
! Pre and Post multiply of A=R^T A R  by rotation matrix. 
            CALL JACPRE(-SINALF,COSALF,P,Q,A,N)
            CALL JACPOS( SINALF,COSALF,P,Q,A,N)
! Accumulate rotations V=R^T V
            CALL JACPRE(-SINALF,COSALF,P,Q,V,N)
      end do ! Was loop 10
!
2000      CONTINUE 
! Put e-values in a vector...
      do Q=1,N
            D(Q)=A(Q,Q)
          END DO
!
          RETURN
          END
!
!
!
!
          SUBROUTINE JACPRE(SINALF,COSALF,P,Q,A,N)
! This sub performs matrix-matrix multiplication A=R*A. 
! PRE-MULTIPLY matrix A by transpose of Rotation matrix 
! is realised by passing -SINALF down into SINALF. 
          IMPLICIT NONE
          INTEGER N
          REAL SINALF,COSALF,A(N,N)
          LOGICAL TRANSP
          INTEGER P,Q
! Local variables...
          INTEGER I
          REAL P1I
!
! Premultiply by rotation matrix...
      do I=1,N
! Row P 1st...
              P1I   =COSALF*A(P,I)-SINALF*A(Q,I)
! Row 2nd put strait in A...
              A(Q,I)=SINALF*A(P,I)+COSALF*A(Q,I)
              A(P,I)=P1I 
            END DO
          RETURN
          END    
!
!
!
!
          SUBROUTINE JACPOS(SINALF,COSALF,P,Q,A,N)
! This sub performs matrix-matrix multiplication A=A*R. 
! POST-MULTIPLY matrix A by transpose of Rotation matrix 
! is realised by passing -SINALF down into SINALF. 
          IMPLICIT NONE
          INTEGER N
          REAL SINALF,COSALF,A(N,N)
          INTEGER P,Q
! Local variables...
          INTEGER I
          REAL IP1
!
! Post multiply by rotation matrix...
      do I=1,N
! Column P 1st...
              IP1   = COSALF*A(I,P)+SINALF*A(I,Q)
! column 2nd put strait in A...
              A(I,Q)=-SINALF*A(I,P)+COSALF*A(I,Q)
              A(I,P)=IP1
            END DO
!
          RETURN
          END    
!
! 




    SUBROUTINE COMB_VEL_MATRIX_DIAG_DIST(DIAG_BIGM_CON, BIGM_CON, &
    DGM_PHA, NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, & ! Force balance sparsity
    NCOLELE, FINELE, COLELE,  NDIM_VEL, NPHASE, U_NLOC, U_NONODS, TOTELE )  ! Element connectivity.
        ! This subroutine combines the distributed and block diagonal for an element
        ! into the matrix DGM_PHA.
        IMPLICIT NONE
        INTEGER, intent( in ) :: NDIM_VEL, NPHASE, U_NLOC, U_NONODS, TOTELE, NCOLDGM_PHA, NCOLELE
        !
        REAL, DIMENSION( :,:,:, :,:,:, : ), intent( in ) :: DIAG_BIGM_CON
        REAL, DIMENSION( :,:,:, :,:,:, : ), intent( in ) :: BIGM_CON
        REAL, DIMENSION( : ), intent( inout ) :: DGM_PHA
        INTEGER, DIMENSION( :), intent( in ) :: FINDGM_PHA
        INTEGER, DIMENSION( :), intent( in ) :: COLDGM_PHA
        INTEGER, DIMENSION(: ), intent( in ) :: FINELE
        INTEGER, DIMENSION( : ), intent( in ) :: COLELE
        ! NEW_ORDERING then order the matrix: IDIM,IPHASE,UILOC,ELE
        ! else use the original ordering...
        LOGICAL, PARAMETER :: NEW_ORDERING = .false.
        LOGICAL, PARAMETER :: tempory_order=.true.

        INTEGER :: ELE,ELE_ROW_START,ELE_ROW_START_NEXT,ELE_IN_ROW
        INTEGER :: U_ILOC,U_JLOC, IPHASE,JPHASE, IDIM,JDIM, I,J, GLOBI, GLOBJ, U_INOD_IDIM_IPHA, U_JNOD_JDIM_JPHA
        INTEGER :: COUNT,COUNT_ELE,JCOLELE
        real, dimension(:,:,:, :,:,:), allocatable :: LOC_DGM_PHA


        ALLOCATE(LOC_DGM_PHA(NDIM_VEL,NDIM_VEL,NPHASE,NPHASE,U_NLOC,U_NLOC))

        Loop_Elements20: DO ELE = 1, TOTELE

            ELE_ROW_START=FINELE(ELE)
            ELE_ROW_START_NEXT=FINELE(ELE+1)
            ELE_IN_ROW = ELE_ROW_START_NEXT - ELE_ROW_START

            ! Block diagonal and off diagonal terms...
            Between_Elements_And_Boundary20: DO COUNT_ELE=ELE_ROW_START, ELE_ROW_START_NEXT-1

                JCOLELE=COLELE(COUNT_ELE)

                IF(JCOLELE==ELE) THEN
                    ! Block diagonal terms (Assume full coupling between the phases and dimensions)...
                    LOC_DGM_PHA(:,:,:, :,:,:) = DIAG_BIGM_CON(:,:,:, :,:,:, ELE) + BIGM_CON(:,:,:, :,:,:, COUNT_ELE)
                ELSE
                    LOC_DGM_PHA(:,:,:, :,:,:) = BIGM_CON(:,:,:, :,:,:, COUNT_ELE)
                ENDIF

                DO U_ILOC=1,U_NLOC
                    DO U_JLOC=1,U_NLOC
                        DO IPHASE=1,NPHASE
                            DO JPHASE=1,NPHASE
                                DO IDIM=1,NDIM_VEL
                                    DO JDIM=1,NDIM_VEL


                                        IF ( NEW_ORDERING ) THEN

                                            ! New for rapid code ordering of variables...
                                            I=IDIM + (IPHASE-1)*NDIM_VEL + (U_ILOC-1)*NDIM_VEL*NPHASE
                                            J=JDIM + (JPHASE-1)*NDIM_VEL + (U_JLOC-1)*NDIM_VEL*NPHASE
                                            COUNT = (COUNT_ELE-1)*(NDIM_VEL*NPHASE)**2 + (I-1)*NDIM_VEL*NPHASE*U_NLOC + J
                                            DGM_PHA(COUNT) = LOC_DGM_PHA(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC)

                                        ELSE

                                            ! Old ordering of the variables BIGM...
                                            GLOBI=(ELE-1)*U_NLOC + U_ILOC
                                            GLOBJ=(JCOLELE-1)*U_NLOC + U_JLOC

                                            if ( tempory_order ) then
                                                I=IDIM + (IPHASE-1)*NDIM_VEL + (U_ILOC-1)*NDIM_VEL*NPHASE
                                                J=JDIM + (JPHASE-1)*NDIM_VEL + (U_JLOC-1)*NDIM_VEL*NPHASE
                                                U_INOD_IDIM_IPHA = I + (ELE-1)*NDIM_VEL*NPHASE*U_NLOC
                                                U_JNOD_JDIM_JPHA = J + (JCOLELE-1)*NDIM_VEL*NPHASE*U_NLOC
                                            else
                                                U_INOD_IDIM_IPHA = GLOBI + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM_VEL*U_NONODS
                                                U_JNOD_JDIM_JPHA = GLOBJ + (JDIM-1)*U_NONODS + ( JPHASE - 1 ) * NDIM_VEL*U_NONODS
                                            end if

                                            COUNT=0
                                            CALL POSINMAT( COUNT, U_INOD_IDIM_IPHA, U_JNOD_JDIM_JPHA, &
                                            U_NONODS * NPHASE * NDIM_VEL, FINDGM_PHA, COLDGM_PHA, NCOLDGM_PHA )
                                            IF(COUNT.NE.0) THEN

                                                !DGM_PHA(COUNT) = LOC_DGM_PHA(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC)
                                                DGM_PHA(COUNT) = DGM_PHA(COUNT)+LOC_DGM_PHA(IDIM,JDIM,IPHASE,JPHASE,U_ILOC,U_JLOC)

                                            END IF

                                        END IF

                                    END DO
                                END DO
                            END DO
                        END DO
                    END DO
                END DO


            END DO Between_Elements_And_Boundary20

        END DO Loop_Elements20


        RETURN
    END SUBROUTINE COMB_VEL_MATRIX_DIAG_DIST





    REAL FUNCTION dg_oscilat_detect(SNDOTQ_KEEP, SNDOTQ2_KEEP, &
    N_DOT_DU, N_DOT_DU2, SINCOME, MASS_ELE, MASS_ELE2 )
        ! Determine if we have an oscillation in the normal direction...
        ! dg_oscilat_detect=1.0- CENTRAL SCHEME.
        ! dg_oscilat_detect=0.0- UPWIND SCHEME.
        real SNDOTQ_KEEP, SNDOTQ2_KEEP, N_DOT_DU, N_DOT_DU2, SINCOME
        REAL MASS_ELE, MASS_ELE2
        ! If cons_oscillation then apply upwinding as often as possible...
        LOGICAL, PARAMETER :: cons_oscillation = .false.
        REAL H1,H2, U1,U2,U3
        !              REAL TOLFUN

        if(cons_oscillation) then

            dg_oscilat_detect = 1.0

            if( SINCOME> 0.5 ) then
                ! velcity comming into element ELE...
                if( (SNDOTQ_KEEP - SNDOTQ2_KEEP)*N_DOT_DU2 > 0.0 ) dg_oscilat_detect = 0.0
            !                   if( (SNDOTQ_KEEP - SNDOTQ2_KEEP)*N_DOT_DU2 > 0.0 ) dg_oscilat_detect = 0.333
            !                   if( (SNDOTQ_KEEP - SNDOTQ2_KEEP)*N_DOT_DU2 > 0.0 ) dg_oscilat_detect = 0.5
            else
                ! velcity pointing out of the element ELE...
                if( (SNDOTQ2_KEEP - SNDOTQ_KEEP)*N_DOT_DU < 0.0 ) dg_oscilat_detect = 0.0
            !                   if( (SNDOTQ2_KEEP - SNDOTQ_KEEP)*N_DOT_DU < 0.0 ) dg_oscilat_detect = 0.333
            !                   if( (SNDOTQ2_KEEP - SNDOTQ_KEEP)*N_DOT_DU < 0.0 ) dg_oscilat_detect = 0.5
            end if
        else
            ! tvd in the means...

            dg_oscilat_detect = 1.0

            if( SINCOME> 0.5 ) then
                ! velcity comming into element ELE...
                if( (SNDOTQ_KEEP - SNDOTQ2_KEEP)*N_DOT_DU2 > 0.0 ) then
                    H1=MASS_ELE
                    H2=MASS_ELE2
                    U1=SNDOTQ_KEEP - N_DOT_DU* H1
                    U2=SNDOTQ2_KEEP+ N_DOT_DU2* H2
                    U3=SNDOTQ2_KEEP+ N_DOT_DU2* 3*H2
                    ! Have oscillations...
                    IF( (U1-U2)/TOLFUN(U2-U3) .LE. 0.0) dg_oscilat_detect = 0.0
                endif
            else
                ! velcity pointing out of the element ELE...
                if( (SNDOTQ2_KEEP - SNDOTQ_KEEP)*N_DOT_DU < 0.0 ) then
                    H1=MASS_ELE
                    H2=MASS_ELE2
                    U1=SNDOTQ_KEEP - N_DOT_DU* 3.*H1
                    U2=SNDOTQ_KEEP - N_DOT_DU* H1
                    U3=SNDOTQ2_KEEP+ N_DOT_DU2* H2
                    ! Have oscillations...
                    IF( (U1-U2)/TOLFUN(U2-U3) .LE. 0.0) dg_oscilat_detect = 0.0
                endif
            end if

        endif

        return
    end function dg_oscilat_detect





    SUBROUTINE USE_POSINMAT_C_STORE(COUNT, U_INOD, P_JNOD,  &
    U_NONODS, FINDC, COLC, NCOLC, &
    IDO_STORE_AC_SPAR_PT,STORED_AC_SPAR_PT, POSINMAT_C_STORE,ELE,U_ILOC,P_JLOC, &
    TOTELE,U_NLOC,P_NLOC)
        INTEGER, intent( inout ) :: COUNT
        INTEGER, intent( in ) :: U_INOD, P_JNOD, U_NONODS,  NCOLC
        INTEGER, intent( in ) :: ELE,U_ILOC,P_JLOC,  TOTELE,U_NLOC,P_NLOC
        INTEGER, intent( in ) :: IDO_STORE_AC_SPAR_PT
        LOGICAL, intent( in ) :: STORED_AC_SPAR_PT
        INTEGER, DIMENSION( U_NLOC,P_NLOC, TOTELE*IDO_STORE_AC_SPAR_PT), intent( inout ) :: POSINMAT_C_STORE
        INTEGER, DIMENSION( U_NONODS + 1), intent( in ) :: FINDC
        INTEGER, DIMENSION( NCOLC), intent( in ) :: COLC

        ! Find COUNT - position in matrix : FINMCY, COLMCY
        IF(STORED_AC_SPAR_PT) THEN
            COUNT=POSINMAT_C_STORE(U_ILOC,P_JLOC,ELE)
        ELSE
            CALL POSINMAT( COUNT, U_INOD, P_JNOD,  &
            U_NONODS, FINDC, COLC, NCOLC )
            IF(IDO_STORE_AC_SPAR_PT.NE.0) POSINMAT_C_STORE(U_ILOC,P_JLOC,ELE) = COUNT
        ENDIF
        RETURN
    END SUBROUTINE USE_POSINMAT_C_STORE




    SUBROUTINE USE_POSINMAT_C_STORE_SUF_DG(COUNT, U_INOD, P_JNOD,  &
    U_NONODS, FINDC, COLC, NCOLC, &
    IDO_STORE_AC_SPAR_PT,STORED_AC_SPAR_PT, POSINMAT_C_STORE_SUF_DG, ELE,IFACE,U_SILOC,P_SJLOC,  &
    TOTELE,NFACE,U_SNLOC,P_SNLOC)
        INTEGER, intent( inout ) :: COUNT
        INTEGER, intent( in ) :: U_INOD, P_JNOD, U_NONODS,  NCOLC
        INTEGER, intent( in ) :: ELE,IFACE,U_SILOC,P_SJLOC,  TOTELE,NFACE,U_SNLOC,P_SNLOC
        INTEGER, intent( in ) :: IDO_STORE_AC_SPAR_PT
        LOGICAL, intent( in ) :: STORED_AC_SPAR_PT
        INTEGER, DIMENSION( U_SNLOC,P_SNLOC,NFACE,TOTELE*IDO_STORE_AC_SPAR_PT ), intent( inout ) :: POSINMAT_C_STORE_SUF_DG
        INTEGER, DIMENSION( U_NONODS + 1), intent( in ) :: FINDC
        INTEGER, DIMENSION( NCOLC), intent( in ) :: COLC
        ! Find COUNT2 - position in matrix : FINMCY, COLMCY
        IF(STORED_AC_SPAR_PT) THEN
            COUNT=POSINMAT_C_STORE_SUF_DG(U_SILOC,P_SJLOC,IFACE,ELE)
        ELSE
            CALL POSINMAT( COUNT, U_INOD, P_JNOD,  &
            U_NONODS, FINDC, COLC, NCOLC )
            IF(IDO_STORE_AC_SPAR_PT.NE.0) POSINMAT_C_STORE_SUF_DG(U_SILOC,P_SJLOC,IFACE,ELE)=COUNT
        ENDIF
        RETURN
    END SUBROUTINE USE_POSINMAT_C_STORE_SUF_DG




    SUBROUTINE DG_DIFFUSION( ELE, U_NLOC, NONODS, LMMAT1, LINVMMAT1, LMMAT, LNNXMAT, LNXNMAT1, LINVMNXNMAT1, AMAT )
        ! Find diffusion contributions at the surface
        implicit none

        INTEGER, intent( in ) :: ELE, U_NLOC, NONODS
        REAL, DIMENSION( :, : ), intent( inout ) :: LMMAT1, LINVMMAT1
        REAL, DIMENSION( :, : ), intent( inout ) :: LMMAT, LNNXMAT
        REAL, DIMENSION( :, : ), intent( inout ) :: LNXNMAT1, LINVMNXNMAT1
        REAL, DIMENSION( :, : ), intent( inout ) :: AMAT
        ! Local
        INTEGER :: ILOC, GLOBI

        ewrite(3,*) 'In DG_DIFFUSION'

        ! LMMAT1
        LMMAT1( 1 : U_NLOC + 1, 1 : U_NLOC + 1 ) = 0.0
        LMMAT1( 1 : U_NLOC    , 1 : U_NLOC ) = LMMAT( 1 : U_NLOC, 1 : U_NLOC )
        LMMAT1( 2 : U_NLOC + 1, 2 : U_NLOC + 1 ) = LMMAT1( 2 : U_NLOC + 1 , 2 : U_NLOC + 1 ) + &
        LMMAT( 1 : U_NLOC     , 1 : U_NLOC )

        ! LNXNMAT1 - surface integral
        LNXNMAT1( 1 : U_NLOC    , 1 : U_NLOC )    =  LNNXMAT( 1 : U_NLOC , 1 : U_NLOC )
        LNXNMAT1( 2 : U_NLOC + 1, 3 : U_NLOC + 2 )=  LNNXMAT( 1 : U_NLOC , 1 : U_NLOC )

        LNXNMAT1( 2, 2 ) = LNXNMAT1( 2, 2 ) - 1.0
        LNXNMAT1( 2, 3 ) = LNXNMAT1( 2, 3 ) + 1.0

        ! Find inverse:
        CALL MATDMATINV( LMMAT1, LINVMMAT1, 2 * U_NLOC - 1 ) ! is the size of LMMAT1 right? DOUBLE CHECK THIS LATER

        ! Matrix X Matrix:
        CALL ABMATRIXMUL( LINVMNXNMAT1, LINVMMAT1, 2 * U_NLOC - 1, 2 * U_NLOC - 1, &
        LNXNMAT1, 2 * U_NLOC - 1, 2 * U_NLOC )

        ! RHS OF ELEMENT:
        ILOC = U_NLOC
        GLOBI = ( ELE - 1 ) * U_NLOC + ILOC
        AMAT( GLOBI, GLOBI - 1)  = AMAT( GLOBI, GLOBI - 1 ) - LINVMNXNMAT1( 2, 1 )
        AMAT( GLOBI, GLOBI )     = AMAT( GLOBI, GLOBI )     - LINVMNXNMAT1( 2, 2 )
        AMAT( GLOBI, GLOBI + 1 ) = AMAT( GLOBI, GLOBI + 1 ) - LINVMNXNMAT1( 2, 3 )
        AMAT( GLOBI, GLOBI + 2 ) = AMAT( GLOBI, GLOBI + 2 ) - LINVMNXNMAT1( 2, 4 )

        ! LHS OF ELEMENT:
        ILOC = 1
        GLOBI = ( ELE - 1 ) * U_NLOC + ILOC
        AMAT( GLOBI, GLOBI - 2 )= AMAT( GLOBI, GLOBI - 2 ) + LINVMNXNMAT1( 2, 1 )
        AMAT( GLOBI, GLOBI - 1 )= AMAT( GLOBI, GLOBI - 1 ) + LINVMNXNMAT1( 2, 2 )
        AMAT( GLOBI, GLOBI )    = AMAT( GLOBI, GLOBI )     + LINVMNXNMAT1( 2, 3 )
        AMAT( GLOBI, GLOBI + 1 )= AMAT( GLOBI, GLOBI + 1 ) + LINVMNXNMAT1( 2, 4 )

        ewrite(3,*) 'Leaving DG_DIFFUSION'

    END SUBROUTINE DG_DIFFUSION




    SUBROUTINE ASSEM_CS( CTP, CT, CTYRHS, FREDOP, NONODS, NCOLCT, FINDCT, COLCT, U, DEN, &
    UBOT, UTOP, DEN_IN_TOP, DEN_IN_BOT,  &
    BOT_BC_TYPE, TOP_BC_TYPE )
        implicit none
        ! assemble CTP (eqn 3.22 without time term) & CT operating on P in eqn 3.21
        ! and also CTYRHS which is the rhs of the cty eqn.
        ! Local variables...
        ! 2 types of B.C's:
        ! BOT_BC_TYPE or TOP_BC_TYPE =3 is a specified inlet velocity & density b.c.
        ! BOT_BC_TYPE or TOP_BC_TYPE =2 is a specified inlet velocity & No density b.c.
        ! BOT_BC_TYPE or TOP_BC_TYPE =1 is No velocity b.c (ZERO PRESSURE BC)& but have density b.c.
        ! BOT_BC_TYPE or TOP_BC_TYPE =0 is an open zero pressure b.c.

        INTEGER, intent( in ) ::  FREDOP, NONODS, NCOLCT
        REAL, DIMENSION( : ), intent( inout ) :: CTP, CT
        INTEGER, DIMENSION( : ), intent( inout ) :: FINDCT
        INTEGER, DIMENSION( : ), intent( inout ) :: COLCT
        REAL, DIMENSION( : ), intent( inout ) :: U
        REAL, DIMENSION( : ), intent( inout ) :: DEN, CTYRHS
        REAL, intent( in ) :: UBOT, UTOP, DEN_IN_TOP, DEN_IN_BOT
        INTEGER, intent( in ) :: BOT_BC_TYPE, TOP_BC_TYPE

        ! Local
        REAL :: NORMX, DENSITY, VEL
        INTEGER PNOD, II, COL, COUNT, COUNT2
        LOGICAL BOT_BC_VEL, BOT_BC_DEN, TOP_BC_VEL, TOP_BC_DEN

        ewrite(3,*) 'In ASSEM_CS'

        BOT_BC_VEL = .FALSE.
        BOT_BC_DEN = .FALSE.
        TOP_BC_VEL = .FALSE.
        TOP_BC_DEN = .FALSE.

        Case_TOP_BC_TYPE: SELECT CASE( TOP_BC_TYPE )
            CASE( 1 ) ; TOP_BC_DEN = .TRUE.
            CASE( 2 ) ; TOP_BC_VEL = .TRUE.
            CASE( 3 )
                TOP_BC_DEN = .TRUE.
                TOP_BC_VEL = .TRUE.
        END SELECT Case_TOP_BC_TYPE

        Case_BOT_BC_TYPE: SELECT CASE( BOT_BC_TYPE )
            CASE( 1 ) ; BOT_BC_DEN = .TRUE.
            CASE( 2 ) ; BOT_BC_VEL = .TRUE.
            CASE( 3 )
                BOT_BC_DEN = .TRUE.
                BOT_BC_VEL = .TRUE.
        END SELECT CASE_BOT_BC_TYPE

        CTP = 0.0
        CT = 0.0
        CTYRHS = 0.0

        ! internal node discretisation
        Loop_Disc: DO PNOD = 2, FREDOP - 1
            Loop_II: DO II = 0, 1
                NORMX = REAL( II * 2 - 1 )
                VEL = U( PNOD + II )

                IF( VEL * NORMX >= 0.0 ) THEN
                    DENSITY = DEN( PNOD )
                ELSE
                    DENSITY = DEN( PNOD + II * 2 - 1 )
                ENDIF

                COL = PNOD + II
                COUNT = 0
                DO COUNT2 = FINDCT( PNOD ) , FINDCT( PNOD + 1 ) - 1
                    IF( COLCT( COUNT2 ) == COL ) COUNT = COUNT2
                END DO
                CT(  COUNT ) = CT(  COUNT ) + NORMX
                CTP( COUNT ) = CTP( COUNT ) + NORMX * DENSITY
            END DO Loop_II
        END DO Loop_Disc

        ! Part of 1st row
        NORMX = 1.0
        VEL = U( 2 )
        IF( VEL * NORMX >= 0.0 ) THEN
            DENSITY = DEN( 1 )
        ELSE
            DENSITY = DEN( 2 )
        ENDIF
        CT( 2 ) = CT( 2 ) + NORMX
        CTP( 2 ) = CTP( 2 ) + NORMX * DENSITY

        ! Part of last st row
        NORMX = -1.0
        VEL = U( NONODS - 1 )
        IF( VEL *NORMX >= 0.0 ) THEN
            DENSITY = DEN( FREDOP )
        ELSE
            DENSITY = DEN( FREDOP - 1 )
        ENDIF
        CT(  NCOLCT - 1 ) = CT(  NCOLCT - 1 ) + NORMX
        CTP( NCOLCT - 1 ) = CTP( NCOLCT - 1 ) + NORMX * DENSITY

        ! Left boundary
        NORMX = -1.0
        VEL = U( 1 )
        IF( BOT_BC_VEL ) VEL = UBOT
        DENSITY = DEN( 1 )
        IF(( VEL * NORMX < 0.0 ) .AND. BOT_BC_DEN ) DENSITY = DEN_IN_BOT

        IF( BOT_BC_VEL ) THEN
            CTYRHS( 1 ) = CTYRHS( 1 ) - DENSITY * NORMX * UBOT
        ELSE
            CT(  1 ) = CT(  1 ) + NORMX
            CTP( 1 ) = CTP( 1 ) + NORMX * DENSITY
        ENDIF

        ! Right boundary
        NORMX = 1.0
        VEL = U( NONODS )
        IF( TOP_BC_VEL ) VEL = UTOP
        DENSITY = DEN( FREDOP )
        IF(( VEL * NORMX <  0.0 ) .AND. TOP_BC_DEN ) DENSITY = DEN_IN_TOP

        IF( TOP_BC_VEL ) THEN
            CTYRHS( FREDOP ) = CTYRHS( FREDOP ) - DENSITY * NORMX *UTOP
        ELSE
            CT(  NCOLCT ) = CT(  NCOLCT ) + NORMX
            CTP( NCOLCT ) = CTP( NCOLCT ) + NORMX * DENSITY
        ENDIF

        ewrite(3,*) 'Leaving ASSEM_CS'

    END SUBROUTINE ASSEM_CS




    SUBROUTINE AVESOU( S2AVE, S2, FREDOP )
        implicit none

        INTEGER, intent( in ) :: FREDOP
        REAL, DIMENSION( : ),     intent( inout ) :: S2AVE, S2
        ! Local
        REAL, DIMENSION( : ), allocatable :: SOURCE
        ! Local variables
        INTEGER :: ELE

        ALLOCATE( SOURCE( FREDOP + 1 ))

        SOURCE( 1 ) = S2( 1 )
        DO ELE= 2, FREDOP
            SOURCE( ELE ) = 0.5 * ( S2( ELE - 1 ) + S2( ELE ))
        END DO
        SOURCE( FREDOP + 1 ) = S2( FREDOP )

        DO ELE= 1, FREDOP
            S2AVE( ELE ) = 0.5 * ( SOURCE( ELE ) + SOURCE( ELE + 1 ))
        END DO

        DEALLOCATE( SOURCE )

    END SUBROUTINE AVESOU




    SUBROUTINE AVESIG( SIGMA2AVE, SIGMA2, FREDOP)
        implicit none

        INTEGER, intent( in ) :: FREDOP
        REAL, DIMENSION( : ), intent( inout ) :: SIGMA2AVE
        REAL, DIMENSION( : ), intent( in ) :: SIGMA2
        ! Local variables
        REAL, DIMENSION( : ), allocatable :: SIGMA
        INTEGER :: ELE

        ALLOCATE( SIGMA( FREDOP + 1 ))

        SIGMA( 1 ) = SIGMA2( 1 )
        DO ELE = 2, FREDOP
            SIGMA( ELE ) = 0.5 * ( SIGMA2( ELE - 1 ) + SIGMA2( ELE ))
        END DO
        SIGMA( FREDOP + 1 ) = SIGMA2( FREDOP )

        DO ELE = 1, FREDOP
            SIGMA2AVE( ELE ) = 0.5 * ( SIGMA( ELE ) + SIGMA( ELE + 1 ))
        END DO

        DEALLOCATE( SIGMA )

    END SUBROUTINE AVESIG




    SUBROUTINE LUMP_ENERGY_EQNS( CV_NONODS, NPHASE, &
    NCOLACV, NCOLACV_SUB, &
    FINACV, COLACV, COLACV_SUB, FINACV_SUB, ACV_SUB )
        implicit none


        INTEGER, intent( in ) :: CV_NONODS, NPHASE, NCOLACV, NCOLACV_SUB
        INTEGER, DIMENSION( : ), intent( in ) :: FINACV
        INTEGER, DIMENSION( : ), intent( in ) :: COLACV
        INTEGER, DIMENSION( : ), intent( inout ) :: COLACV_SUB
        INTEGER, DIMENSION( : ), intent( inout ) :: FINACV_SUB
        REAL, DIMENSION( :), intent( inout ) :: ACV_SUB
        ! Local Variables
        INTEGER :: COUNT, COUNT2, CV_NOD, ICOL, ICOL_PHA, CV_NOD_PHA

        ewrite(3,*) 'In LUMP_ENERGY_EQNS'

        COUNT2 = 0

        DO CV_NOD = 1, CV_NONODS
            FINACV_SUB( CV_NOD ) = COUNT2 + 1
            DO COUNT = FINACV( CV_NOD ), FINACV( CV_NOD + 1 ) - 1, 1
                ICOL = COLACV( COUNT )
                IF(ICOL <= CV_NONODS) THEN
                    COUNT2 = COUNT2 + 1
                    COLACV_SUB( COUNT2 ) = ICOL
                END IF
            END DO
        END DO
        FINACV_SUB( CV_NONODS + 1 ) = COUNT2 + 1


        ACV_SUB = 0.
        DO CV_NOD_PHA = 1, CV_NONODS * NPHASE
            DO COUNT = 1, FINACV( CV_NOD_PHA + 1 ) - 1
                CV_NOD = MOD( CV_NOD_PHA, CV_NONODS )
                ICOL_PHA = COLACV( COUNT )
                ICOL = MOD ( ICOL_PHA, CV_NONODS )

                CALL POSINMAT( COUNT2, CV_NOD, ICOL, &
                CV_NONODS, FINACV_SUB, COLACV_SUB, NCOLACV_SUB )

            END DO
        END DO

        ewrite(3,*) 'Leaving LUMP_ENERGY_EQNS'

    END SUBROUTINE LUMP_ENERGY_EQNS


    SUBROUTINE CALCULATE_SURFACE_TENSION( state, packed_state, nphase, ncomp, &
    PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, IPLIKE_GRAD_SOU, &
    U_SOURCE_CV, U_SOURCE, &
    COMP, &
    NCOLACV, FINACV, COLACV, MIDACV, &
    SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
    block_to_global_acv, global_dense_block_acv, &
    NCOLCT, FINDCT, COLCT, &
    CV_NONODS, U_NONODS, X_NONODS, TOTELE, STOTEL, &
    CV_ELE_TYPE, CV_SELE_TYPE, U_ELE_TYPE, &
    CV_NLOC, U_NLOC, X_NLOC, CV_SNLOC, U_SNLOC, &
    CV_NDGLN, CV_SNDGLN, X_NDGLN, U_NDGLN, U_SNDGLN, &
    MAT_NLOC, MAT_NDGLN, MAT_NONODS,  &
    NDIM,  &
    NCOLM, FINDM, COLM, MIDM, &
    XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
    StorageIndexes )

        IMPLICIT NONE

        real, dimension( cv_nonods * nphase ), intent( inout ) :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD
        integer, intent( inout ) :: IPLIKE_GRAD_SOU
        real, dimension( cv_nonods * nphase * ndim ), intent( inout ) :: U_SOURCE_CV
        real, dimension( u_nonods * nphase * ndim ), intent( inout ) :: U_SOURCE

        type(state_type), dimension( : ), intent( inout ) :: state
        type(state_type), intent( inout ) :: packed_state
        integer, intent( in ) :: nphase, ncomp, cv_nonods, U_NONODS, X_NONODS, MAT_NONODS, &
        NCOLACV, NCOLCT, TOTELE, CV_ELE_TYPE, CV_SELE_TYPE, U_ELE_TYPE, &
        CV_NLOC, U_NLOC, X_NLOC, MAT_NLOC, CV_SNLOC, U_SNLOC, NDIM, &
        NCOLM, XU_NLOC, NCOLELE, STOTEL
        integer, dimension( : ), intent( in ) :: CV_NDGLN
        integer, dimension( :), intent( in )  :: CV_SNDGLN
        integer, dimension( : ), intent( in ) ::  X_NDGLN
        integer, dimension( : ), intent( in ) :: U_NDGLN
        integer, dimension( : ), intent( in ) :: U_SNDGLN
        integer, dimension( : ), intent( in ) :: XU_NDGLN
        integer, dimension( : ), intent( in ) :: MAT_NDGLN
        integer, dimension( : ), intent( in ) :: FINACV
        integer, dimension( : ), intent( in ) :: COLACV
        integer, dimension( : ), intent( in ) :: MIDACV
        integer, dimension(:), intent(in) :: small_finacv,small_colacv,small_midacv
        integer, dimension(:), intent(in) :: block_to_global_acv
        integer, dimension(:,:), intent(in) :: global_dense_block_acv
        integer, dimension( : ), intent( in ) :: FINDCT
        integer, dimension( : ), intent( in ) :: COLCT

        real, dimension( : ), intent( in ) :: COMP

        integer, dimension( : ), intent( in ) :: FINDM
        integer, dimension( : ), intent( in ) :: COLM
        integer, dimension( : ), intent( in ) :: MIDM
        integer, dimension( : ), intent( in ) :: FINELE
        integer, dimension( : ), intent( in ) :: COLELE
        integer, dimension(:), intent(inout) ::  StorageIndexes
        !Local variables
        real, dimension( : ), allocatable :: U_FORCE_X_SUF_TEN, U_FORCE_Y_SUF_TEN, U_FORCE_Z_SUF_TEN, &
        CV_U_FORCE_X_SUF_TEN, CV_U_FORCE_Y_SUF_TEN, CV_U_FORCE_Z_SUF_TEN, X, Y, Z
        real, dimension( STOTEL * CV_SNLOC ) :: DUMMY_SUF_COMP_BC
        integer, dimension( STOTEL ) :: DUMMY_WIC_COMP_BC

        integer :: iphase, icomp
        real :: coefficient
        logical :: surface_tension, use_pressure_force, use_smoothing

        type( vector_field ), pointer :: x_all


        ewrite(3,*) 'Entering CALCULATE_SURFACE_TENSION'

        allocate( X(  X_NONODS ) ) ; X = 0.0
        allocate( Y(  X_NONODS ) ) ; Y = 0.0
        allocate( Z(  X_NONODS ) ) ; Z = 0.0

        x_all => extract_vector_field( packed_state, "PressureCoordinate" )
        x = x_all % val( 1, : )
        if (ndim >=2 ) y = x_all % val( 2, : )
        if (ndim >=3 ) z = x_all % val( 3, : )


        ! Initialise...
        IPLIKE_GRAD_SOU = 0
        PLIKE_GRAD_SOU_COEF = 0.0
        !For capillary pressure these terms already have a value, so overwritting is a problem
        if( .not. have_option( '/material_phase[0]/multiphase_properties/capillary_pressure' ) )  then
            PLIKE_GRAD_SOU_GRAD = 0.0
        end if
        U_SOURCE_CV = 0.0

        DUMMY_SUF_COMP_BC = 0.0
        DUMMY_WIC_COMP_BC = 0

        do icomp = 1, ncomp

            surface_tension = have_option( '/material_phase[' // int2str( nphase - 1 + icomp ) // &
            ']/is_multiphase_component/surface_tension' )

            if ( surface_tension ) then

                ewrite(3,*) 'Calculating surface tension for component ', icomp

                call get_option( '/material_phase[' // int2str( nphase - 1 + icomp ) // &
                ']/is_multiphase_component/surface_tension/coefficient', coefficient )

                use_smoothing = have_option( '/material_phase[' // int2str( nphase - 1 + icomp ) // &
                ']/is_multiphase_component/surface_tension/smooth' )

                allocate( U_FORCE_X_SUF_TEN( U_NONODS) ) ; U_FORCE_X_SUF_TEN = 0.0
                allocate( U_FORCE_Y_SUF_TEN( U_NONODS) ) ; U_FORCE_Y_SUF_TEN = 0.0
                allocate( U_FORCE_Z_SUF_TEN( U_NONODS) ) ; U_FORCE_Z_SUF_TEN = 0.0

                allocate( CV_U_FORCE_X_SUF_TEN( CV_NONODS) ) ; CV_U_FORCE_X_SUF_TEN = 0.0
                allocate( CV_U_FORCE_Y_SUF_TEN( CV_NONODS) ) ; CV_U_FORCE_Y_SUF_TEN = 0.0
                allocate( CV_U_FORCE_Z_SUF_TEN( CV_NONODS) ) ; CV_U_FORCE_Z_SUF_TEN = 0.0

                USE_PRESSURE_FORCE = .TRUE.

                if ( USE_PRESSURE_FORCE ) then
                    IPLIKE_GRAD_SOU = 1
                else
                    IPLIKE_GRAD_SOU = 0
                end if

                do iphase = 1, nphase

                    CALL SURFACE_TENSION_WRAPPER( state, packed_state, &
                    U_FORCE_X_SUF_TEN, U_FORCE_Y_SUF_TEN, U_FORCE_Z_SUF_TEN, &
                    CV_U_FORCE_X_SUF_TEN, CV_U_FORCE_Y_SUF_TEN, CV_U_FORCE_Z_SUF_TEN, &
                    PLIKE_GRAD_SOU_COEF( 1+CV_NONODS*(IPHASE-1) : CV_NONODS*IPHASE ), & 
                    PLIKE_GRAD_SOU_GRAD( 1+CV_NONODS*(IPHASE-1) : CV_NONODS*IPHASE ), &
                    COEFFICIENT, &
                    COMP( 1 + (IPHASE-1)*CV_NONODS + (ICOMP-1)*NPHASE*CV_NONODS : &
                    IPHASE*CV_NONODS + (ICOMP-1)*NPHASE*CV_NONODS ), &
                    NCOLACV, FINACV, COLACV, MIDACV, &
                    SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
                    block_to_global_acv, global_dense_block_acv, &
                    NCOLCT, FINDCT, COLCT, &
                    CV_NONODS, U_NONODS, X_NONODS, TOTELE, STOTEL, &
                    CV_ELE_TYPE, CV_SELE_TYPE, U_ELE_TYPE, &
                    CV_NLOC, U_NLOC, X_NLOC, CV_SNLOC, U_SNLOC, &
                    CV_NDGLN, CV_SNDGLN, X_NDGLN, U_NDGLN, U_SNDGLN, &
                    X, Y, Z, &
                    MAT_NLOC, MAT_NDGLN, MAT_NONODS,  &
                    NDIM, USE_PRESSURE_FORCE, &
                    NCOLM, FINDM, COLM, MIDM, &
                    XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
                    USE_SMOOTHING,&
                    StorageIndexes )

                end do

                if ( .not.USE_PRESSURE_FORCE ) then

                    !U_SOURCE_CV(1:cv_nonods) = CV_U_FORCE_X_SUF_TEN
                    !U_SOURCE_CV(1+cv_nonods:2*cv_nonods) = CV_U_FORCE_Y_SUF_TEN

                    U_SOURCE(1:U_nonods) = U_FORCE_X_SUF_TEN
                    U_SOURCE(1+U_nonods:2*U_nonods) = U_FORCE_Y_SUF_TEN

                end if

                deallocate( U_FORCE_X_SUF_TEN, U_FORCE_Y_SUF_TEN, U_FORCE_Z_SUF_TEN )
                deallocate( CV_U_FORCE_X_SUF_TEN, CV_U_FORCE_Y_SUF_TEN, CV_U_FORCE_Z_SUF_TEN )

            end if

        end do

        ewrite(3,*) 'Leaving CALCULATE_SURFACE_TENSION'

        RETURN
    END SUBROUTINE CALCULATE_SURFACE_TENSION





    SUBROUTINE SURFACE_TENSION_WRAPPER( state, packed_state, &
    U_FORCE_X_SUF_TEN, U_FORCE_Y_SUF_TEN, U_FORCE_Z_SUF_TEN, &
    CV_U_FORCE_X_SUF_TEN, CV_U_FORCE_Y_SUF_TEN, CV_U_FORCE_Z_SUF_TEN, &
    PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, &
    SUF_TENSION_COEF, VOLUME_FRAC, &
    NCOLACV, FINACV, COLACV, MIDACV, &
    SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
    block_to_global_acv, global_dense_block_acv, &
    NCOLCT, FINDCT, COLCT, &
    CV_NONODS, U_NONODS, X_NONODS, TOTELE, STOTEL, &
    CV_ELE_TYPE, CV_SELE_TYPE, U_ELE_TYPE, &
    CV_NLOC, U_NLOC, X_NLOC, CV_SNLOC, U_SNLOC, &
    CV_NDGLN, CV_SNDGLN, X_NDGLN, U_NDGLN, U_SNDGLN, &
    X, Y, Z, &
    MAT_NLOC, MAT_NDGLN, MAT_NONODS,  &
    NDIM, USE_PRESSURE_FORCE, &
    NCOLM, FINDM, COLM, MIDM, &
    XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
    USE_SMOOTHING,&
    StorageIndexes )

        ! Calculate the surface tension force: U_FORCE_X_SUF_TEN,U_FORCE_X_SUF_TEN,U_FORCE_X_SUF_TEN
        ! or PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD,
        ! for a given volume fraction field VOLUME_FRAC
        ! SUF_TENSION_COEF is the surface tension coefficient.

        !  =====================================================================
        !     In this subroutine the advection terms in the advection-diffusion
        !     equation (in the matrix and RHS) are calculated as ACV and CV_RHS.
        !
        !     This routine uses a Control Volume (CV) formulation to compute
        !     the advection terms. The general procedure is as follows:
        !
        !        1. For each node-pair, define which node is the donor, which is
        !           the receptor and define an "upwind" value of the field being
        !           advected (and the accompanying "density"; see note below)
        !        2. Calculate the volume flux across the CV face that separates
        !           these two nodes
        !        3. Estimate the value of the advected variable at the control
        !           volume face.
        !        4. Using information from the donor, receptor and upwind nodes,
        !           limit the field face value (removes oscillations from the
        !           solution)
        !        5. Assemble the fluxes to form the matrix and rhs of the
        !           advection equation
        !
        !     This procedure is implemented by considering the CV to be made up
        !     of a number of sub-control-volumes, which represent the part of
        !     the control volume within a given element.  The assembly of terms
        !     considers each of these sub-CVs in turn, calculating (and limiting)
        !     the flux across sub-CV faces that are external to the CV...
        !
        !     NOTE: Add in note about what density is in this sub!!!
        !
        !     To define the "upwind" value of the field variable, which is
        !     necessary for the limiting scheme, either:
        !
        !        A. The upwind value of the field variable to be advected is
        !           found by interpolation and stored in a matrix (TUPWIND)
        !        B. The neighbouring nodes are searched for the local maximum
        !           and minimum
        !
        !     The subroutine has several options...
        !
        !     Discretisation option
        !     ---------------------
        !      - The estimate of the face value may be determined in one of
        !        several ways.
        !      - The face value may be centered in time by either a specified
        !        CV_THETA value, or a non-linear CV_THETA value that is determined
        !        automatically.
        !      - The face value may be limited using a univeral-limiter-type
        !        scheme, or a limited-downwind scheme that is ideal for INTERFACE
        !        TRACKING.  Alternatively no limiting can be applied.
        !
        !     These options are defined by the value of CV_DISOPT, which corresponds
        !     to the clast digit of the GEM option NDISOT for the field in question.
        !
        !     CV_DISOPT=discretisation option in space and time
        !     ------------------------------------------------------------------
        !     CV_DISOPT   Method for face-value est.    Time-stepping     Limiting
        !     ------------------------------------------------------------------
        !       =0      1st order in space          Theta=specified    UNIVERSAL
        !       =1      1st order in space          Theta=non-linear   UNIVERSAL
        !       =2      Trapezoidal rule in space   Theta=specified    UNIVERSAL
        !       =2 if isotropic limiter then FEM-quadratic & stratification adjust. Theta=non-linear
        !       =3      Trapezoidal rule in space   Theta=non-linear   UNIVERSAL
        !       =4      Finite elements in space    Theta=specified    UNIVERSAL
        !       =5      Finite elements in space    Theta=non-linear   UNIVERSAL
        !       =6      Finite elements in space    Theta=specified    NONE
        !       =7      Finite elements in space    Theta=non-linear   NONE
        !       =8      Finite elements in space    Theta=specified    DOWNWIND+
        !       =9      Finite elements in space    Theta=non-linear   DOWNWIND+
        !
        !     CV_DG_VEL_INT_OPT=interface velocity calculation option between elements
        !
        !     Limiting scheme
        !     ---------------
        !     The limiting scheme is defined in the subroutine NVDFUNNEW;
        !     the limited values are computed in subroutine ANONVDLIM/ONVDLIM.
        !
        !     ONVDLIM is the original limiting algorithm
        !
        !     ANONVDLIM is a new anisoptropic limiting algorithm, which is
        !     called if either ALOLIM=1 (where ALOLIM is an option flag set
        !     in this subroutine), or if the interface tracking limiting option
        !     is selected (CV_DISOPT=8/9).  ***In general ALOLIM appears to be set to 1 (GSC)
        !
        !     NOTE: ANONVDLIM only works for TETS; for all other element types
        !     ONVDLIM is used.
        !
        !
        !     IMPORTANT INPUTS:
        !     ----------------
        !
        !     ACV   - Matrix for assembling the advection terms (empty on input)
        !     CV_RHS      - Right-hand side vector for advection-diffusion terms
        !     X,Y,Z    - Node co-ordinates
        !     NU       - Nodal velocity component
        !     T,TOLD   - New and old advected field values at nodes
        !     DEN,  - New and old "density" at nodes, which is actually a constant
        !     DENOLD     multiplying the advection diffusion equation for the field
        !     CV_DISOPT   - The discretisation/limiting option (see above)
        !     DT       - The time step
        !     CV_THETA    - The time-stepping discretisation parameter
        !     CV_BETA     - Conservative(1.)/non-conservative(0.) flag
        !     ELE_TYP   - Integer flag definining element type
        !
        !     IMPORTANT OUTPUTS:
        !     -----------------
        !
        !     ACV   - Matrix updated to include the advection terms
        !     CV_RHS      - Right-hand side vector updated to include advection terms
        !
        !
        !     IMPORTANT LOCAL PARAMETERS:
        !     --------------------------
        !
        !     TIMOPT    - Temporal discretisation option, derived from CV_DISOPT.
        !                (1 for non-linear theta; 0 for theta specified (THETA))
        !
        !
        !***********************************************************************
        use shape_functions
        use matrix_operations
        use printout
        ! Inputs/Outputs
        IMPLICIT NONE
        type(state_type), dimension( : ), intent( inout ) :: state
        type(state_type), intent( inout ) :: packed_state

        INTEGER, PARAMETER :: NPHASE = 1
        INTEGER, PARAMETER :: SMOOTH_NITS = 0 ! smoothing iterations, 10 seems good.
        INTEGER, intent( in ) :: NCOLACV, NCOLCT, CV_NONODS, U_NONODS, X_NONODS, MAT_NONODS, &
        TOTELE, STOTEL, &
        CV_ELE_TYPE, CV_SELE_TYPE, U_ELE_TYPE, &
        CV_NLOC, U_NLOC, X_NLOC, MAT_NLOC, &
        CV_SNLOC, U_SNLOC, NDIM, &
        NCOLM, XU_NLOC, NCOLELE
        INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
        INTEGER, DIMENSION( : ), intent( in )  :: CV_SNDGLN
        INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: U_SNDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
        INTEGER, DIMENSION( : ), intent( in ) :: FINACV
        INTEGER, DIMENSION( : ), intent( in ) :: COLACV
        INTEGER, DIMENSION( : ), intent( in ) :: MIDACV
        integer, dimension(:), intent(in) :: small_finacv,small_colacv,small_midacv
        integer, dimension(:), intent(in) :: block_to_global_acv
        integer, dimension(:, :), intent(in) :: global_dense_block_acv

        INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
        INTEGER, DIMENSION( : ), intent( in ) :: COLCT

        REAL, intent( in ) ::  SUF_TENSION_COEF

        REAL, DIMENSION( : ), intent( inout ) :: U_FORCE_X_SUF_TEN,U_FORCE_Y_SUF_TEN,U_FORCE_Z_SUF_TEN
        REAL, DIMENSION( : ), intent( inout ) :: CV_U_FORCE_X_SUF_TEN,CV_U_FORCE_Y_SUF_TEN,CV_U_FORCE_Z_SUF_TEN
        REAL, DIMENSION( : ), intent( inout ) :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD

        REAL, DIMENSION( : ), intent( in ) :: VOLUME_FRAC

        REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
        INTEGER, DIMENSION( : ), intent( in ) :: FINDM
        INTEGER, DIMENSION( : ), intent( in ) :: COLM
        INTEGER, DIMENSION( : ), intent( in ) :: MIDM
        INTEGER, DIMENSION( : ), intent( in ) :: FINELE
        INTEGER, DIMENSION( : ), intent( in ) :: COLELE
        LOGICAL, intent( in ) :: USE_PRESSURE_FORCE, USE_SMOOTHING
        integer, dimension(:), intent(inout) ::  StorageIndexes
        ! Local variables
        LOGICAL, DIMENSION( : ), allocatable :: X_SHARE,LOG_ON_BOUND
        LOGICAL, DIMENSION( :, : ), allocatable :: CV_ON_FACE, U_ON_FACE, &
        CVFEM_ON_FACE, UFEM_ON_FACE
        INTEGER, DIMENSION( : ), allocatable :: FINDGPTS, &
        CV_OTHER_LOC, U_OTHER_LOC, MAT_OTHER_LOC, &
        JCOUNT_KLOC, JCOUNT_KLOC2, COLGPTS, CV_SLOC2LOC, U_SLOC2LOC, &
        TMAX_NOD, TMIN_NOD, TOLDMAX_NOD, &
        TOLDMIN_NOD, DENMAX_NOD, DENMIN_NOD, DENOLDMAX_NOD, DENOLDMIN_NOD, &
        T2MAX_NOD, T2MIN_NOD, T2OLDMAX_NOD, T2OLDMIN_NOD, IDUM, IZERO, DG_CV_NDGLN
        INTEGER, DIMENSION( : , : ), allocatable :: CV_SLOCLIST, U_SLOCLIST, &
        FACE_ELE, CV_NEILOC
        REAL, DIMENSION( : ), allocatable :: CVWEIGHT, CVWEIGHT_SHORT, SCVFEWEIGH, SBCVFEWEIGH, &
        CVNORMX, &
        CVNORMY, CVNORMZ, MASS_CV, MASS_ELE, SNDOTQ, SNDOTQOLD,  &
        FEMT, SHARP_FEMT,FEMTOLD, FEMTOLD2,FEMT2, FEMT2OLD, FEMDEN, FEMDENOLD, XC_CV, YC_CV, ZC_CV, &
        SCVDETWEI, SRA, UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
        UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2,  &
        SUM_CV, ONE_PORE, SELE_OVERLAP_SCALE, &
        T2MAX, T2MIN, T2OLDMAX, &
        T2OLDMIN, &
        T2MAX_2ND_MC, T2MIN_2ND_MC, T2OLDMAX_2ND_MC, &
        T2OLDMIN_2ND_MC, &
        UP_WIND_NOD, DU, DV, DW, RDUM, RZERO, CURVATURE, CV_ONE
        REAL, DIMENSION( : ), allocatable :: CV_FORCE_X_SUF_TEN, CV_FORCE_Y_SUF_TEN, CV_FORCE_Z_SUF_TEN
        REAL, DIMENSION( : , : ), allocatable :: CVN, CVN_SHORT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
        CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT,  &
        UFEN, UFENLX, UFENLY, UFENLZ, SCVFEN, SCVFENSLX, SCVFENSLY, &
        SCVFENLX, SCVFENLY, SCVFENLZ,  &
        SUFEN, SUFENSLX, SUFENSLY, SUFENLX, SUFENLY, SUFENLZ, &
        SBCVN, SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
        SBCVFENLX, SBCVFENLY, SBCVFENLZ, SBUFEN, SBUFENSLX, SBUFENSLY, &
        SBUFENLX, SBUFENLY, SBUFENLZ, &
        DUMMY_ZERO_NDIM_NDIM, RZERO_DIAGTEN
        REAL, DIMENSION( : , : ), allocatable :: MASS, STORE_MASS
        integer, dimension(:), allocatable :: IPIV
        REAL, DIMENSION( : , :, : ), allocatable :: DTX_ELE,DTY_ELE,DTZ_ELE, &
        SHARP_DTX_ELE,SHARP_DTY_ELE,SHARP_DTZ_ELE, &
        DTOLDX_ELE,DTOLDY_ELE,DTOLDZ_ELE
        REAL, DIMENSION( : ), allocatable :: B_CV_X,B_CV_Y,B_CV_Z, &
        RHS_U_SHORT_X,RHS_U_SHORT_Y,RHS_U_SHORT_Z, &
        U_SOL_X,U_SOL_Y,U_SOL_Z, &
        DIF_TX, DIF_TY, DIF_TZ, &
        MASS_NORMALISE, &
        TAU_XX, TAU_XY, TAU_XZ, &
        TAU_YX, TAU_YY, TAU_YZ, &
        TAU_ZX, TAU_ZY, TAU_ZZ
        REAL, DIMENSION( :, :, : ), allocatable ::&
        DX_TAU_XX, DY_TAU_XY, DZ_TAU_XZ, &
        DX_TAU_YX, DY_TAU_YY, DZ_TAU_YZ, &
        DX_TAU_ZX, DY_TAU_ZY, DZ_TAU_ZZ, &
        DX_DIFF_X, DY_DIFF_X, DZ_DIFF_X, &
        DX_DIFF_Y, DY_DIFF_Y, DZ_DIFF_Y, &
        DX_DIFF_Z, DY_DIFF_Z, DZ_DIFF_Z, rzero3
        REAL, DIMENSION( :,:,:,: ), allocatable :: THERM_U_DIFFUSION

        !        ===> INTEGERS <===
        INTEGER :: CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, COUNT, JCOUNT, &
        ELE, ELE2, GI, GCOUNT, SELE, &
        NCOLGPTS, &
        CV_SILOC, U_ILOC, U_JLOC, U_KLOC, &
        CV_ILOC, CV_JLOC, IPHASE, JPHASE, &
        CV_NODJ, CV_NODJ_IPHA, &
        CV_NODI, CV_NODI_IPHA, CV_NODI_JPHA, U_NODK, TIMOPT, &
        JCOUNT_IPHA, IMID_IPHA, &
        NFACE, X_NODI,  U_INOD, U_NOD, &
        CV_INOD, CV_JNOD, MAT_NODI, FACE_ITS, NFACE_ITS, &
        CVNOD, XNOD, CV_NOD, DG_CV_NOD, IDIM, IGOT_T2, &
        nopt_vel_upwind_coefs, DG_CV_NONODS, IGOT_THERM_VIS
        !        ===>  REALS  <===
        REAL :: NDOTQ, NDOTQOLD,  &
        INCOME, INCOMEOLD, HDC, FVT, FVTOLD, FVT2, FVT2OLD, &
        FVD, FVDOLD, LIMT, LIMTOLD, LIMT2, LIMT2OLD,&
        LIMD, LIMDOLD, FTHETA, VTHETA, &
        LIMDT, LIMDTOLD, LIMDTT2, LIMDTT2OLD, &
        FEMDGI, FEMTGI,FEMT2GI, FEMDOLDGI, FEMTOLDGI, FEMT2OLDGI, &
        TMID, TOLDMID, &
        DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX, BCZERO, ROBIN1, ROBIN2, &
        SUM, &
        SUM_LIMT, SUM_LIMTOLD, FTHETA_T2, ONE_M_FTHETA_T2OLD, THERM_FTHETA, &
        W_SUM_ONE1, W_SUM_ONE2, NDOTQNEW, NN, NM, DT, T_THETA, T_BETA, RDIF, RR, &
        RSUM, RRSUM, rr2, grad_c_x,grad_c_y,grad_c_z

        REAL, PARAMETER :: W_SUM_ONE = 1.0, TOLER=1.0E-10

        integer :: cv_inod_ipha, IGETCT, U_NODK_IPHA, NOIT_DIM, &
        CV_DG_VEL_INT_OPT, IN_ELE_UPWIND, DG_ELE_UPWIND, &
        CV_DISOPT, IGOT_THETA_FLUX, scvngi_theta,SMOOTH_ITS
        ! Functions...
        !REAL :: R2NORM, FACE_THETA
        !        ===>  LOGICALS  <===
        LOGICAL :: GETMAT, LIMIT_USE_2ND, &
        D1, D3, DCYL, GOT_DIFFUS, INTEGRAT_AT_GI, &
        NORMALISE, SUM2ONE, GET_GTHETA, QUAD_OVER_WHOLE_ELE, GETCT
        LOGICAL :: GET_THETA_FLUX, USE_THETA_FLUX, THERMAL, LUMP_EQNS, &
        SIMPLE_LINEAR_SCHEME, GOTDEC, STRESS_FORM

        REAL FEMT_CV_NOD(CV_NLOC)

        CHARACTER(LEN=OPTION_PATH_LEN) :: OPTION_PATH
        REAL, allocatable, DIMENSION(:) :: DUMMY_ELE

        real, dimension(0,0,0,0):: tflux
        real, allocatable, dimension(:,:,:) :: T_ABSORB
        real, allocatable, dimension(:,:,:,:) :: tdiffusion
        real, dimension(0,0) :: ALIMTOLD,ALIMT2OLD,ALIMDOLD,ALIMDTOLD,ALIMDTT2OLD,ANDOTQOLD
        !Pointer
        real, pointer, dimension(:,:,:) :: CVFENX_ALL, UFENX_ALL
        real, pointer, dimension(:) :: DETWEI, RA
        real, pointer :: VOLUME

        REAL, DIMENSION(1,1) :: DUMMY_THETA_GDIFF
        type(tensor_field) :: tfield
      

        ALLOCATE(DUMMY_ELE(TOTELE))
        DUMMY_ELE = 0

        IGOT_T2=0
        CV_DISOPT=0
        CV_DG_VEL_INT_OPT=0
        IN_ELE_UPWIND=0
        DG_ELE_UPWIND=0
        GETCT=.FALSE.
        GET_THETA_FLUX=.FALSE.
        USE_THETA_FLUX=.FALSE.
        THERMAL=.FALSE.
        LIMIT_USE_2ND=.FALSE.

        ALLOCATE(RDUM(MAX(U_NLOC,CV_NLOC)*TOTELE)) ; RDUM = 0.0
        ALLOCATE(IDUM(MAX(U_NLOC,CV_NLOC)*TOTELE)) ; IDUM = 0
        ALLOCATE(RZERO(MAX(U_NLOC,CV_NLOC)*TOTELE)) ; RZERO=0.0
        ALLOCATE(RZERO3(MAX(U_NLOC,CV_NLOC),1,TOTELE)); RZERO3=0.0
        ALLOCATE(RZERO_DIAGTEN(CV_SNLOC*STOTEL*NPHASE, NDIM)) ; RZERO_DIAGTEN=0.0
        ALLOCATE(IZERO(MAX(U_NLOC,CV_NLOC)*TOTELE))  ; IZERO=0
        ALLOCATE(CV_ONE(CV_NONODS)) ; CV_ONE=1.0
        ALLOCATE(CURVATURE(CV_NONODS))
        NOPT_VEL_UPWIND_COEFS=0

        ndotq = 0. ; ndotqold = 0.

        QUAD_OVER_WHOLE_ELE=.FALSE.
        ! If QUAD_OVER_WHOLE_ELE=.true. then dont divide element into CV's to form quadrature.
        call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
        cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE )

        GOT_DIFFUS = .true.
        ALLOCATE(CV_FORCE_X_SUF_TEN(CV_NONODS))
        ALLOCATE(CV_FORCE_Y_SUF_TEN(CV_NONODS))
        ALLOCATE(CV_FORCE_Z_SUF_TEN(CV_NONODS))

        ! Allocate memory for the control volume surface shape functions, etc.
        ALLOCATE( JCOUNT_KLOC(  U_NLOC )) ; jcount_kloc = 0
        ALLOCATE( JCOUNT_KLOC2(  U_NLOC )) ; jcount_kloc2 = 0

        ALLOCATE( CVNORMX( SCVNGI ))
        ALLOCATE( CVNORMY( SCVNGI ))
        ALLOCATE( CVNORMZ( SCVNGI ))
        ALLOCATE( COLGPTS( CV_NLOC * SCVNGI )) !The size of this vector is over-estimated
        ALLOCATE( FINDGPTS( CV_NLOC + 1 ))
        ALLOCATE( SNDOTQ( SCVNGI ))
        ALLOCATE( SNDOTQOLD( SCVNGI ))
        ALLOCATE( CV_ON_FACE( CV_NLOC, SCVNGI ))
        ALLOCATE( CVFEM_ON_FACE( CV_NLOC, SCVNGI ))
        ALLOCATE( U_ON_FACE( U_NLOC, SCVNGI ))
        ALLOCATE( UFEM_ON_FACE( U_NLOC, SCVNGI ))
        ALLOCATE( CV_OTHER_LOC( CV_NLOC ))
        ALLOCATE( U_OTHER_LOC( U_NLOC ))
        ALLOCATE( MAT_OTHER_LOC( MAT_NLOC ))
        ALLOCATE( X_SHARE( X_NONODS ))
        ALLOCATE( CVWEIGHT( CV_NGI ))
        ALLOCATE( CVN( CV_NLOC, CV_NGI ))
        ALLOCATE( CVFEN( CV_NLOC, CV_NGI ))
        ALLOCATE( CVFENLX( CV_NLOC, CV_NGI ))
        ALLOCATE( CVFENLY( CV_NLOC, CV_NGI ))
        ALLOCATE( CVFENLZ( CV_NLOC, CV_NGI ))


        ALLOCATE( CVWEIGHT_SHORT( CV_NGI_SHORT ))
        ALLOCATE( CVN_SHORT( CV_NLOC, CV_NGI_SHORT ))
        ALLOCATE( CVFEN_SHORT( CV_NLOC, CV_NGI_SHORT))
        ALLOCATE( CVFENLX_SHORT( CV_NLOC, CV_NGI_SHORT ))
        ALLOCATE( CVFENLY_SHORT( CV_NLOC, CV_NGI_SHORT ))
        ALLOCATE( CVFENLZ_SHORT( CV_NLOC, CV_NGI_SHORT ))

        ALLOCATE( UFEN( U_NLOC, CV_NGI))
        ALLOCATE( UFENLX( U_NLOC, CV_NGI ))
        ALLOCATE( UFENLY( U_NLOC, CV_NGI ))
        ALLOCATE( UFENLZ( U_NLOC, CV_NGI ))


        ALLOCATE( SCVFEN( CV_NLOC, SCVNGI ))
        ALLOCATE( SCVFENSLX( CV_NLOC, SCVNGI ))
        ALLOCATE( SCVFENSLY( CV_NLOC, SCVNGI ))
        ALLOCATE( SCVFENLX( CV_NLOC, SCVNGI ))
        ALLOCATE( SCVFENLY( CV_NLOC, SCVNGI ))
        ALLOCATE( SCVFENLZ( CV_NLOC, SCVNGI ))
        ALLOCATE( SCVFEWEIGH( SCVNGI ))

        ALLOCATE( SUFEN( U_NLOC, SCVNGI ))
        ALLOCATE( SUFENSLX( U_NLOC, SCVNGI ))
        ALLOCATE( SUFENSLY( U_NLOC, SCVNGI ))
        ALLOCATE( SUFENLX( U_NLOC, SCVNGI ))
        ALLOCATE( SUFENLY( U_NLOC, SCVNGI ))
        ALLOCATE( SUFENLZ( U_NLOC, SCVNGI ))

        ALLOCATE( SCVDETWEI( SCVNGI )) ; SCVDETWEI = 0.
        ALLOCATE( SRA( SCVNGI ))
        ALLOCATE( LOG_ON_BOUND(CV_NONODS))

        ALLOCATE( SBCVN( CV_SNLOC, SBCVNGI ))
        ALLOCATE( SBCVFEN( CV_SNLOC, SBCVNGI ))
        ALLOCATE( SBCVFENSLX( CV_SNLOC, SBCVNGI ))
        ALLOCATE( SBCVFENSLY( CV_SNLOC, SBCVNGI ))
        ALLOCATE( SBCVFEWEIGH( SBCVNGI ))
        ALLOCATE( SBCVFENLX( CV_SNLOC, SBCVNGI ))
        ALLOCATE( SBCVFENLY( CV_SNLOC, SBCVNGI ))
        ALLOCATE( SBCVFENLZ( CV_SNLOC, SBCVNGI ))
        ALLOCATE( SBUFEN( U_SNLOC, SBCVNGI ))
        ALLOCATE( SBUFENSLX( U_SNLOC, SBCVNGI ))
        ALLOCATE( SBUFENSLY( U_SNLOC, SBCVNGI ))
        ALLOCATE( SBUFENLX( U_SNLOC, SBCVNGI ))
        ALLOCATE( SBUFENLY( U_SNLOC, SBCVNGI ))
        ALLOCATE( SBUFENLZ( U_SNLOC, SBCVNGI ))
        ALLOCATE( DUMMY_ZERO_NDIM_NDIM(NDIM,NDIM))
        DUMMY_ZERO_NDIM_NDIM=0.0

        ALLOCATE( CV_SLOC2LOC( CV_SNLOC ))
        ALLOCATE( U_SLOC2LOC( U_SNLOC ))
        ALLOCATE( CV_SLOCLIST( NFACE, CV_SNLOC ))
        ALLOCATE( U_SLOCLIST( NFACE, U_SNLOC ))
        ALLOCATE( CV_NEILOC( CV_NLOC, SCVNGI ))

        ALLOCATE( SELE_OVERLAP_SCALE(CV_NLOC) )

        ALLOCATE( UGI_COEF_ELE(U_NLOC),  VGI_COEF_ELE(U_NLOC),  WGI_COEF_ELE(U_NLOC) )
        ALLOCATE( UGI_COEF_ELE2(U_NLOC), VGI_COEF_ELE2(U_NLOC), WGI_COEF_ELE2(U_NLOC) )
        ! The procity mapped to the CV nodes
        ALLOCATE( SUM_CV( CV_NONODS ))
        ALLOCATE( UP_WIND_NOD( CV_NONODS * NPHASE )) ; UP_WIND_NOD = 0.0

        D1 = ( NDIM == 1 )
        D3 = ( NDIM == 3 )
        DCYL= ( NDIM == -2 )

        GETMAT = .TRUE.

        X_SHARE = .FALSE.

        ! If using the original limiting scheme, the first step is to estimate
        ! the upwind field value from the surrounding nodes

        ! Allocate memory for terms needed by GETGXYZ OR ONVDLIM

        !     ======= DEFINE THE SUB-CONTROL VOLUME & FEM SHAPE FUNCTIONS ========

        CALL CV_FEM_SHAPE_FUNS( &
                             ! Volume shape functions...
        NDIM, CV_ELE_TYPE,  &
        CV_NGI, CV_NGI_SHORT, CV_NLOC, U_NLOC, CVN, CVN_SHORT, &
        CVWEIGHT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
        CVWEIGHT_SHORT, CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
        UFEN, UFENLX, UFENLY, UFENLZ, &
                             ! Surface of each CV shape functions...
        SCVNGI, CV_NEILOC, CV_ON_FACE, CVFEM_ON_FACE, &
        SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
        SCVFENLX, SCVFENLY, SCVFENLZ,  &
        SUFEN, SUFENSLX, SUFENSLY,  &
        SUFENLX, SUFENLY, SUFENLZ,  &
                             ! Surface element shape funcs...
        U_ON_FACE, UFEM_ON_FACE, NFACE, &
        SBCVNGI, SBCVN, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
        SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, SBUFENLY, SBUFENLZ, &
        CV_SLOCLIST, U_SLOCLIST, CV_SNLOC, U_SNLOC, &
                             ! Define the gauss points that lie on the surface of the CV...
        FINDGPTS, COLGPTS, NCOLGPTS, &
        SELE_OVERLAP_SCALE, QUAD_OVER_WHOLE_ELE )

        ! Determine FEMT (finite element wise) etc from T (control volume wise)
        ! Also determine the CV mass matrix MASS_CV and centre of the CV's XC_CV,YC_CV,ZC_CV.
        ! This is for projecting to finite element basis functions...
        ALLOCATE( FEMT( CV_NONODS * NPHASE ))
        ALLOCATE( SHARP_FEMT( CV_NONODS * NPHASE ))
        ALLOCATE( FEMTOLD( CV_NONODS * NPHASE ))
        ALLOCATE( FEMTOLD2( CV_NONODS * NPHASE ))
        ALLOCATE( MASS_CV( CV_NONODS ))
        ALLOCATE( MASS_ELE( TOTELE ))
        ALLOCATE( XC_CV( CV_NONODS ))
        ALLOCATE( YC_CV( CV_NONODS ))
        ALLOCATE( ZC_CV( CV_NONODS ))
        ALLOCATE( DTX_ELE( CV_NLOC, NPHASE, TOTELE ))
        ALLOCATE( DTY_ELE( CV_NLOC, NPHASE, TOTELE ))
        ALLOCATE( DTZ_ELE( CV_NLOC, NPHASE, TOTELE ))
        ALLOCATE( SHARP_DTX_ELE( CV_NLOC, NPHASE, TOTELE ))
        ALLOCATE( SHARP_DTY_ELE( CV_NLOC, NPHASE, TOTELE ))
        ALLOCATE( SHARP_DTZ_ELE( CV_NLOC, NPHASE, TOTELE ))
        ALLOCATE( DTOLDX_ELE( CV_NLOC, NPHASE, TOTELE ))
        ALLOCATE( DTOLDY_ELE( CV_NLOC, NPHASE, TOTELE ))
        ALLOCATE( DTOLDZ_ELE( CV_NLOC, NPHASE, TOTELE ))

        IGETCT=0
        IF(GETCT) IGETCT=1

        option_path='/material_phase[0]/scalar_field::Pressure'

        CALL PROJ_CV_TO_FEM( FEMT, VOLUME_FRAC, 1, NDIM, &
        RDUM,0, RDUM,0, MASS_ELE, &
        CV_NONODS, TOTELE, CV_NDGLN, X_NLOC, X_NDGLN, &
        CV_NGI_SHORT, CV_NLOC, CVN_SHORT, CVWEIGHT_SHORT, &
        CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
        X_NONODS, X, Y, Z, NCOLM, FINDM, COLM, MIDM, &
        IGETCT, RDUM, IDUM, IDUM, 0, OPTION_PATH )

        FEMT=1.0-VOLUME_FRAC
        FEMTOLD=0.0

        SHARP_FEMT=FEMT


        if(.false.) then ! mide side node average...
            DO ELE=1,TOTELE
                DO CV_ILOC=1,CV_NLOC
                    FEMT_CV_NOD(CV_ILOC)=SHARP_FEMT(CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC))
                END DO
                FEMT_CV_NOD(2)=0.5*(FEMT_CV_NOD(1)+FEMT_CV_NOD(3))
                FEMT_CV_NOD(4)=0.5*(FEMT_CV_NOD(1)+FEMT_CV_NOD(6))
                FEMT_CV_NOD(5)=0.5*(FEMT_CV_NOD(3)+FEMT_CV_NOD(6))
                DO CV_ILOC=1,CV_NLOC
                    SHARP_FEMT(CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC))=FEMT_CV_NOD(CV_ILOC)
                END DO
            END DO
            FEMT=SHARP_FEMT
        endif

        !       Smooth FEMT...
        if(.false.) then
            !         DO SMOOTH_ITS=1,SMOOTH_NITS
            DO SMOOTH_ITS=1,3
                !     DO SMOOTH_ITS=1,20
                DO CV_NOD=1,CV_NONODS
                    RSUM=0.0
                    RRSUM=0.0
                    DO COUNT=FINACV(CV_NOD),FINACV(CV_NOD+1)-1
                        IF(COLACV(COUNT).LE.CV_NONODS) THEN
                            !                     RSUM=RSUM+FEMT(COLACV(COUNT))
                            RSUM=RSUM+SHARP_FEMT(COLACV(COUNT))
                            RRSUM=RRSUM+1.0
                        ENDIF
                    END DO
                    !               FEMTOLD(CV_NOD)=0.5*FEMT(CV_NOD)+0.5*RSUM/RRSUM
                    FEMTOLD(CV_NOD)=0.5*SHARP_FEMT(CV_NOD)+0.5*RSUM/RRSUM
                   !         FEMTOLD(CV_NOD)=0.75*FEMT(CV_NOD)+0.25*RSUM/RRSUM
                   !         FEMTOLD(CV_NOD)=0.9*FEMT(CV_NOD)+0.1*RSUM/RRSUM
                END DO
                !            FEMT=FEMTOLD
                SHARP_FEMT=FEMTOLD
                FEMTOLD=0.0
            END DO
        endif


        ALLOCATE( FACE_ELE( NFACE, TOTELE ) ) ; FACE_ELE = 0
        ! Calculate FACE_ELE
        CALL CALC_FACE_ELE( FACE_ELE, TOTELE, STOTEL, NFACE, &
        NCOLELE, FINELE, COLELE, CV_NLOC, CV_SNLOC, CV_NONODS, CV_NDGLN, CV_SNDGLN, &
        CV_SLOCLIST, X_NLOC, X_NDGLN )

        CALL DG_DERIVS( FEMT, FEMTOLD, &
        DTX_ELE, DTY_ELE, DTZ_ELE, DTOLDX_ELE, DTOLDY_ELE, DTOLDZ_ELE, &
        NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, &
        X_NDGLN, X_NLOC, X_NDGLN, &
        CV_NGI_SHORT, CV_NLOC, CVWEIGHT_SHORT, &
        CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
        CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
        X_NONODS, X, Y, Z,  &
        NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, IZERO, &
        RZERO, &
        1, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, &
        SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
        state, "wrap1", StorageIndexes(19) )


        CALL DG_DERIVS( SHARP_FEMT, FEMTOLD, &
        SHARP_DTX_ELE, SHARP_DTY_ELE, SHARP_DTZ_ELE, DTOLDX_ELE, DTOLDY_ELE, DTOLDZ_ELE, &
        NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, &
        X_NDGLN, X_NLOC, X_NDGLN, &
        CV_NGI_SHORT, CV_NLOC, CVWEIGHT_SHORT, &
        CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
        CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
        X_NONODS, X, Y, Z,  &
        NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, IZERO, &
        RZERO, &
        1, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, &
        SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
        state, "Surf_ten_wrap2", StorageIndexes(20))

        ! determine the curvature by solving a simple eqn...

        ALLOCATE( TDIFFUSION( NDIM, NDIM, CV_NONODS,nphase ) ) ; TDIFFUSION=0.0
        ALLOCATE( MASS_NORMALISE( CV_NONODS ) ) ; MASS_NORMALISE=0.0
        DO ELE=1,TOTELE
            DO CV_ILOC=1,CV_NLOC
                CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                MASS_NORMALISE(CV_NOD) = MASS_NORMALISE(CV_NOD) + MASS_ELE(ELE)
            END DO
        END DO

        ! smooth...
        if ( USE_SMOOTHING ) then
            femtold=0.0
            DO ELE=1,TOTELE
                DO CV_ILOC=1,CV_NLOC
                    CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                    femtold(cv_nod)=femtold(cv_nod)+DTX_ELE(CV_ILOC, 1, ELE) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
                END DO
            END DO
            DO ELE=1,TOTELE
                DO CV_ILOC=1,CV_NLOC
                    CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                    DTX_ELE(CV_ILOC, 1, ELE) = femtold(cv_nod)
                END DO
            END DO

            femtold=0.0
            DO ELE=1,TOTELE
                DO CV_ILOC=1,CV_NLOC
                    CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                    femtold(cv_nod)=femtold(cv_nod)+DTY_ELE(CV_ILOC, 1, ELE) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
                END DO
            END DO
            DO ELE=1,TOTELE
                DO CV_ILOC=1,CV_NLOC
                    CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                    DTY_ELE(CV_ILOC, 1, ELE) = femtold(cv_nod)
                END DO
            END DO
        endif

        ! smooth sharp...
        if(.false.) then
            femtold=0.0
            DO ELE=1,TOTELE
                DO CV_ILOC=1,CV_NLOC
                    CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                    femtold(cv_nod)=femtold(cv_nod)+SHARP_DTX_ELE(CV_ILOC, 1, ELE) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
                END DO
            END DO
            DO ELE=1,TOTELE
                DO CV_ILOC=1,CV_NLOC
                    CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                    SHARP_DTX_ELE(CV_ILOC, 1, ELE) = femtold(cv_nod)
                END DO
            END DO

            femtold=0.0
            DO ELE=1,TOTELE
                DO CV_ILOC=1,CV_NLOC
                    CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                    femtold(cv_nod)=femtold(cv_nod)+SHARP_DTY_ELE(CV_ILOC, 1, ELE) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
                END DO
            END DO
            DO ELE=1,TOTELE
                DO CV_ILOC=1,CV_NLOC
                    CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                    SHARP_DTY_ELE(CV_ILOC, 1, ELE) = femtold(cv_nod)
                END DO
            END DO
        endif

        DO ELE=1,TOTELE
            DO CV_ILOC=1,CV_NLOC
                CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                RR = DTX_ELE(CV_ILOC, 1, ELE)**2
                IF(NDIM.GE.2) RR = RR+ DTY_ELE(CV_ILOC, 1, ELE)**2
                IF(NDIM.GE.3) RR = RR+ DTZ_ELE(CV_ILOC, 1, ELE)**2
                RDIF = 1.0 / MAX( TOLER, SQRT(RR) )
                DO IDIM=1,NDIM
                    TDIFFUSION(IDIM,IDIM,CV_NOD,1) = TDIFFUSION(IDIM,IDIM,CV_NOD,1) + &
                    RDIF * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
                END DO
            END DO
        END DO

        SIMPLE_LINEAR_SCHEME=.TRUE.
        STRESS_FORM=.false.

        IF ( SIMPLE_LINEAR_SCHEME ) THEN

            ! Direct linear scheme

            DG_CV_NONODS=CV_NLOC*TOTELE

            ALLOCATE(DIF_TX(DG_CV_NONODS)) ; DIF_TX=0.0
            ALLOCATE(DIF_TY(DG_CV_NONODS)) ; DIF_TY=0.0
            ALLOCATE(DIF_TZ(DG_CV_NONODS)) ; DIF_TZ=0.0

            if ( stress_form ) then
                ALLOCATE(TAU_XX(DG_CV_NONODS), TAU_XY(DG_CV_NONODS), TAU_XZ(DG_CV_NONODS))
                ALLOCATE(TAU_YX(DG_CV_NONODS), TAU_YY(DG_CV_NONODS), TAU_YZ(DG_CV_NONODS))
                ALLOCATE(TAU_ZX(DG_CV_NONODS), TAU_ZY(DG_CV_NONODS), TAU_ZZ(DG_CV_NONODS))

                TAU_XX=0.0 ; TAU_XY=0.0 ; TAU_XZ=0.0
                TAU_YX=0.0 ; TAU_YY=0.0 ; TAU_YZ=0.0
                TAU_ZX=0.0 ; TAU_ZY=0.0 ; TAU_ZZ=0.0
            end if

            !print *,'SUF_TENSION_COEF:',SUF_TENSION_COEF
            !stop 822

            DO ELE=1,TOTELE
                DO CV_ILOC=1,CV_NLOC
                    CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                    DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC

                    RR = DTX_ELE(CV_ILOC, 1, ELE)**2
                    IF(NDIM.GE.2) RR = RR+ DTY_ELE(CV_ILOC, 1, ELE)**2
                    IF(NDIM.GE.3) RR = RR+ DTZ_ELE(CV_ILOC, 1, ELE)**2
                    RDIF = 1.0 / MAX( TOLER, SQRT(RR) )
                    !               RDIF = 1.0 / MAX( 1.e-5, SQRT(RR) )

                    DIF_TX(DG_CV_NOD)=RDIF * DTX_ELE(CV_ILOC, 1, ELE)
                    IF(NDIM.GE.2) DIF_TY(DG_CV_NOD)=RDIF * DTY_ELE(CV_ILOC, 1, ELE)
                    IF(NDIM.GE.3) DIF_TZ(DG_CV_NOD)=RDIF * DTZ_ELE(CV_ILOC, 1, ELE)

                    ! for stress form...
                    if ( stress_form ) then

                        TAU_XX(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                        * DTX_ELE(CV_ILOC, 1, ELE ) * DTX_ELE(CV_ILOC, 1, ELE) - SQRT(RR) )
                        IF(NDIM.GE.2) TAU_XY(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                        * DTX_ELE(CV_ILOC, 1, ELE ) * DTY_ELE(CV_ILOC, 1, ELE ) )
                        IF(NDIM.GE.3) TAU_XZ(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                        * DTX_ELE(CV_ILOC, 1, ELE ) * DTZ_ELE(CV_ILOC, 1, ELE) )

                        IF(NDIM.GE.2) THEN
                            TAU_YX(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                            * DTY_ELE(CV_ILOC, 1, ELE ) * DTX_ELE(CV_ILOC, 1, ELE) )
                            TAU_YY(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                            * DTY_ELE(CV_ILOC, 1, ELE ) * DTY_ELE(CV_ILOC, 1, ELE ) - SQRT(RR) )
                            IF(NDIM.GE.3) TAU_YZ(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                            * DTY_ELE(CV_ILOC, 1, ELE ) * DTZ_ELE( CV_ILOC, 1, ELE) )
                        ENDIF
                        IF(NDIM.GE.3) THEN
                            TAU_ZX(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                            * DTZ_ELE(CV_ILOC, 1, ELE ) * DTX_ELE(CV_ILOC, 1, ELE) )
                            TAU_ZY(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                            * DTZ_ELE(CV_ILOC, 1, ELE ) * DTY_ELE(CV_ILOC, 1, ELE) )
                            TAU_ZZ(DG_CV_NOD)=-SUF_TENSION_COEF*(RDIF &
                            * DTZ_ELE(CV_ILOC, 1, ELE ) * DTZ_ELE( CV_ILOC, 1, ELE) - SQRT(RR) )
                        ENDIF

                    end if

                END DO
            END DO

            ALLOCATE( DG_CV_NDGLN( DG_CV_NONODS ) )
            DO ELE=1,TOTELE
                DO CV_ILOC=1,CV_NLOC
                    DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC
                    DG_CV_NDGLN(DG_CV_NOD)=DG_CV_NOD
                END DO
            END DO

            if ( stress_form ) then

                ALLOCATE(DX_TAU_XX(CV_NLOC,nphase,TOTELE), DY_TAU_XY(CV_NLOC,nphase,TOTELE), DZ_TAU_XZ(CV_NLOC,nphase,TOTELE))
                ALLOCATE(DX_TAU_YX(CV_NLOC,nphase,TOTELE), DY_TAU_YY(CV_NLOC,nphase,TOTELE), DZ_TAU_YZ(CV_NLOC,nphase,TOTELE))
                ALLOCATE(DX_TAU_ZX(CV_NLOC,nphase,TOTELE), DY_TAU_ZY(CV_NLOC,nphase,TOTELE), DZ_TAU_ZZ(CV_NLOC,nphase,TOTELE))

                DX_TAU_XX=0.0 ; DY_TAU_XY=0.0 ; DZ_TAU_XZ=0.0
                DX_TAU_YX=0.0 ; DY_TAU_YY=0.0 ; DZ_TAU_YZ=0.0
                DX_TAU_ZX=0.0 ; DY_TAU_ZY=0.0 ; DZ_TAU_ZZ=0.0

                CALL DG_DERIVS_UVW( TAU_XX, TAU_XX, TAU_XY, TAU_XY, TAU_XZ, TAU_XZ, &
                DX_TAU_XX, RZERO3, RZERo3, RZERo3, Rzero3, Rzero3, &
                Rzero3, DY_TAU_XY, Rzero3, Rzero3, Rzero3, Rzero3, &
                rzero3, Rzero3, DZ_TAU_XZ, rzero3, rzero3, rzero3, &
                NDIM, NDIM, NPHASE, DG_CV_NONODS, TOTELE, DG_CV_NDGLN, &
                X_NDGLN, X_NLOC, X_NDGLN, &
                CV_NGI, CV_NLOC, CVWEIGHT, &
                CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
                CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
                X_NONODS, X, Y, Z, &
                NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, IZERO,  &
                RZERO,RZERO,RZERO, &
                1, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, &
                SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
                state, "wrapp1", StorageIndexes(4:6))

                U_FORCE_X_SUF_TEN = pack(DX_TAU_XX(:,1,:) + DY_TAU_XY(:,1,:) + DZ_TAU_XZ(:,1,:),.true.)
                !!$            femtold=0.0
                !!$            DO ELE=1,TOTELE
                !!$               DO CV_ILOC=1,CV_NLOC
                !!$                  CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                !!$                  dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
                !!$                  femtold(cv_nod)=femtold(cv_nod)+U_FORCE_X_SUF_TEn(dg_cv_nod) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
                !!$               END DO
                !!$            END DO
                !!$            DO ELE=1,TOTELE
                !!$               DO CV_ILOC=1,CV_NLOC
                !!$                  CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                !!$                  dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
                !!$                  U_FORCE_X_SUF_TEn(dg_cv_nod) = femtold(cv_nod)
                !!$               END DO
                !!$            END DO

                IF(NDIM.GE.2) THEN
                    CALL DG_DERIVS_UVW( TAU_YX, TAU_YX, TAU_YY, TAU_YY, TAU_YZ, TAU_YZ, &
                    DX_TAU_YX, RZERO3, RZERo3, RZERo3, Rzero3, Rzero3, &
                    Rzero3, DY_TAU_YY, Rzero3, Rzero3, Rzero3, Rzero3, &
                    rzero3, Rzero3, DZ_TAU_YZ, rzero3, rzero3, rzero3,&
                    NDIM, NDIM, NPHASE, DG_CV_NONODS, TOTELE, DG_CV_NDGLN, &
                    X_NDGLN, X_NLOC, X_NDGLN, &
                    CV_NGI, CV_NLOC, CVWEIGHT, &
                    CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
                    CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
                    X_NONODS, X, Y, Z, &
                    NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, IZERO,  &
                    RZERO,RZERO,RZERO, &
                    1, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, & 
                    SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
                    state, "wrapp2", StorageIndexes(7:9))

                    U_FORCE_Y_SUF_TEN = pack(DX_TAU_YX(:,1,:)+ DY_TAU_YY(:,1,:) + DZ_TAU_YZ(:,1,:),.true.)


                !!$               femtold=0.0
                !!$               DO ELE=1,TOTELE
                !!$                  DO CV_ILOC=1,CV_NLOC
                !!$                     CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                !!$                     dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
                !!$                     femtold(cv_nod)=femtold(cv_nod)+U_FORCE_Y_SUF_TEn(dg_cv_nod) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
                !!$               END DO
                !!$            END DO
                !!$            DO ELE=1,TOTELE
                !!$               DO CV_ILOC=1,CV_NLOC
                !!$                  CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                !!$                  dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
                !!$                  U_FORCE_Y_SUF_TEn(dg_cv_nod) = femtold(cv_nod)
                !!$               END DO
                !!$            END DO



                ENDIF

                IF(NDIM.GE.3) THEN
                    CALL DG_DERIVS_UVW( TAU_ZX, TAU_ZX, TAU_ZY, TAU_ZY, TAU_ZZ, TAU_ZZ, &
                    DX_TAU_ZX,  RZERO3, RZERo3, RZERo3, Rzero3, Rzero3, &
                    Rzero3, DY_TAU_ZY, Rzero3, Rzero3, Rzero3, Rzero3, &
                    rzero3, Rzero3, DZ_TAU_ZZ, rzero3, rzero3, rzero3, &
                    NDIM, NDIM, NPHASE, DG_CV_NONODS, TOTELE, DG_CV_NDGLN, &
                    X_NDGLN, X_NLOC, X_NDGLN, &
                    CV_NGI, CV_NLOC, CVWEIGHT, &
                    CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
                    CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
                    X_NONODS, X, Y, Z, &
                    NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, IZERO,  &
                    RZERO,RZERO,RZERO, &
                    1, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, & 
                    SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
                    state, "wrapp3", StorageIndexes(10:12))

                    U_FORCE_Z_SUF_TEN = pack(DX_TAU_ZX(:,1,:) + DY_TAU_ZY(:,1,:) + DZ_TAU_ZZ(:,1,:),.true.)

                !!$               femtold=0.0
                !!$               DO ELE=1,TOTELE
                !!$                  DO CV_ILOC=1,CV_NLOC
                !!$                     CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                !!$                     dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
                !!$                     femtold(cv_nod)=femtold(cv_nod)+U_FORCE_Z_SUF_TEn(dg_cv_nod) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
                !!$               END DO
                !!$            END DO
                !!$            DO ELE=1,TOTELE
                !!$               DO CV_ILOC=1,CV_NLOC
                !!$                  CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                !!$                  dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
                !!$                  U_FORCE_Z_SUF_TEn(dg_cv_nod) = femtold(cv_nod)
                !!$               END DO
                !!$            END DO


                ENDIF


                DEALLOCATE(DX_TAU_XX, DY_TAU_XY, DZ_TAU_XZ, &
                DX_TAU_YX, DY_TAU_YY, DZ_TAU_YZ, &
                DX_TAU_ZX, DY_TAU_ZY, DZ_TAU_ZZ)

                DEALLOCATE(TAU_XX, TAU_XY, TAU_XZ, &
                TAU_YX, TAU_YY, TAU_YZ, &
                TAU_ZX, TAU_ZY, TAU_ZZ)

            else ! non stress form

                ALLOCATE(DX_DIFF_X(CV_NLOC,nphase,TOTELE), DY_DIFF_X(CV_NLOC,nphase,TOTELE), DZ_DIFF_X(CV_NLOC,nphase,TOTELE))
                ALLOCATE(DX_DIFF_Y(CV_NLOC,nphase,TOTELE), DY_DIFF_Y(CV_NLOC,nphase,TOTELE), DZ_DIFF_Y(CV_NLOC,nphase,TOTELE))
                ALLOCATE(DX_DIFF_Z(CV_NLOC,nphase,TOTELE), DY_DIFF_Z(CV_NLOC,nphase,TOTELE), DZ_DIFF_Z(CV_NLOC,nphase,TOTELE))


                DX_DIFF_X=0. ;  DY_DIFF_X=0. ; DZ_DIFF_X=0.
                DX_DIFF_Y=0. ;  DY_DIFF_Y=0. ; DZ_DIFF_Y=0.
                DX_DIFF_Z=0. ;  DY_DIFF_Z=0. ; DZ_DIFF_Z=0.

                if(.true.) then

                    CALL DG_DERIVS_UVW( DIF_TX, DIF_TX, DIF_TY, DIF_TY, DIF_TZ, DIF_TZ, &
                    DX_DIFF_X, DY_DIFF_X, DZ_DIFF_X, rzero3, rzero3, rzero3, &
                    DX_DIFF_Y, DY_DIFF_Y, DZ_DIFF_Y, rzero3, rzero3, rzero3, &
                    DX_DIFF_Z, DY_DIFF_Z, DZ_DIFF_Z, rzero3, rzero3, rzero3, &
                    NDIM, NDIM, NPHASE, DG_CV_NONODS, TOTELE, DG_CV_NDGLN, &
                    X_NDGLN, X_NLOC, X_NDGLN, &
                    CV_NGI, CV_NLOC, CVWEIGHT, &
                    CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
                    CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
                    X_NONODS, X, Y, Z, &
                    NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, IZERO,  &
                    RZERO,RZERO,RZERO, &
                    1, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, & 
                    SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
                    state, "wrapp4", StorageIndexes(13:15))

                else

                    femtold=0.0
                    DO ELE=1,TOTELE
                        DO CV_ILOC=1,CV_NLOC
                            CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                            dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
                            femtold(cv_nod)=femtold(cv_nod)+DIF_TX(dg_cv_nod) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
                        END DO
                    END DO

                    CALL DG_DERIVS( FEMTOLD, rzero, &
                    DToldX_ELE, rzero3, rzero3,   rzero3, rzero3, rzero3, &
                    NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, &
                    X_NDGLN, X_NLOC, X_NDGLN, &
                    CV_NGI_SHORT, CV_NLOC, CVWEIGHT_SHORT, &
                    CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
                    CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
                    X_NONODS, X, Y, Z,  &
                    NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, IZERO, &
                    RZERO, &
                    1, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, &
                    SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
                    state, "wrap4", StorageIndexes(21))

                    DO ELE=1,TOTELE
                        DO CV_ILOC=1,CV_NLOC
                            dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
                            DX_DIFF_X(cv_ILOC,1,ele)=DToldX_ELE(CV_ILOC, 1, ELE)
                        END DO
                    END DO


                    femtold=0.0
                    DO ELE=1,TOTELE
                        DO CV_ILOC=1,CV_NLOC
                            CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                            dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
                            femtold(cv_nod)=femtold(cv_nod)+DIF_TY(dg_cv_nod) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
                        END DO
                    END DO

                    CALL DG_DERIVS( FEMTOLD, rzero, &
                    rzero3, DToldY_ELE, rzero3,   rzero3, rzero3, rzero3, & 
                    NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, &
                    X_NDGLN, X_NLOC, X_NDGLN, &
                    CV_NGI_SHORT, CV_NLOC, CVWEIGHT_SHORT, &
                    CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
                    CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
                    X_NONODS, X, Y, Z,  &
                    NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, IZERO, &
                    RZERO, &
                    1, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, &
                    SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
                    state, "wrap5", StorageIndexes(22))

                    DO ELE=1,TOTELE
                        DO CV_ILOC=1,CV_NLOC
                            dg_cv_nod=(ELE-1)*CV_NLOC+CV_ILOC
                            DY_DIFF_Y(CV_ILOC, 1, ELE)=DToldY_ELE(CV_ILOC, 1, ELE)
                        END DO
                    END DO

                endif

                CURVATURE=0.0
                DO ELE=1,TOTELE
                    DO CV_ILOC=1,CV_NLOC
                        CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                        DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC
                        RR=DX_DIFF_X(CV_ILOC,1,ele)
                        IF(NDIM.GE.2) RR=RR + DY_DIFF_Y(CV_ILOC,1,ele)
                        IF(NDIM.GE.3) RR=RR + DZ_DIFF_Z(CV_ILOC,1,ele)
                        CURVATURE(CV_NOD) = CURVATURE(CV_NOD) + RR * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
                    END DO
                END DO

               !DEALLOCATE(DIF_TX, DIF_TY, DIF_TZ)
               !DEALLOCATE(DG_CV_NDGLN)
               !DEALLOCATE(DX_DIFF_X, DY_DIFF_X, DZ_DIFF_X)
               !DEALLOCATE(DX_DIFF_Y, DY_DIFF_Y, DZ_DIFF_Y)
               !DEALLOCATE(DX_DIFF_Z, DY_DIFF_Z, DZ_DIFF_Z)

            end if

        ELSE
            ALLOCATE(T_ABSORB(CV_NONODS,nphase,nphase)) ; T_ABSORB=1.0
            DT=1.0
            T_THETA=0.0
            T_BETA=0.0
            NOIT_DIM=1
            LUMP_EQNS=.FALSE.

            IGOT_THERM_VIS=0
            ALLOCATE( THERM_U_DIFFUSION(NDIM,NDIM,NPHASE,MAT_NONODS*IGOT_THERM_VIS ) )


            CALL INTENERGE_ASSEM_SOLVE( state, packed_state, &
                 tfield, tfield,tfield,&
            NCOLACV, FINACV, COLACV, MIDACV, &
            SMALL_FINACV, SMALL_COLACV, SMALL_MIDACV, &
            block_to_global_acv, global_dense_block_acv, &
            NCOLCT, FINDCT, COLCT, &
            CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
            U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE,  &
            NPHASE,  &
            CV_NLOC, U_NLOC, X_NLOC,  &
            CV_NDGLN, X_NDGLN, U_NDGLN, &
            CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
            CURVATURE, VOLUME_FRAC, &
            MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, IGOT_THERM_VIS, THERM_U_DIFFUSION, &
            CV_DISOPT, CV_DG_VEL_INT_OPT, DT, T_THETA, T_BETA, &
            RZERO_DIAGTEN, &
            RZERO, &
            RZERO, T_ABSORB, RZERO, &
            NDIM,  &
            NCOLM, FINDM, COLM, MIDM, &
            XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
            RDUM, NOPT_VEL_UPWIND_COEFS, &
            RDUM, CV_ONE, &
            IGOT_T2, CURVATURE, VOLUME_FRAC,IGOT_THETA_FLUX, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
            DUMMY_THETA_GDIFF, &
            IN_ELE_UPWIND, DG_ELE_UPWIND, &
            NOIT_DIM, &
            ! nits_flux_lim_t
            RZERO, &
            option_path = '/material_phase[0]/scalar_field::Pressure', &
            mass_ele_transp = dummy_ele, &
            thermal = .FALSE.,&
            StorageIndexes=StorageIndexes, icomp=-1)

            DEALLOCATE(T_ABSORB)

        END IF

        IF_USE_PRESSURE_FORCE: IF ( USE_PRESSURE_FORCE ) THEN

            ! should be minus because is discretised as a pressure term

            !PLIKE_GRAD_SOU_COEF = PLIKE_GRAD_SOU_COEF - SUF_TENSION_COEF * ABS( CURVATURE )
            !         PLIKE_GRAD_SOU_COEF = PLIKE_GRAD_SOU_COEF + SUF_TENSION_COEF * max(0.0,CURVATURE)
            PLIKE_GRAD_SOU_COEF = PLIKE_GRAD_SOU_COEF + SUF_TENSION_COEF * CURVATURE

            !PLIKE_GRAD_SOU_GRAD = PLIKE_GRAD_SOU_GRAD + VOLUME_FRAC
            !PLIKE_GRAD_SOU_GRAD = PLIKE_GRAD_SOU_GRAD + FEMT
            PLIKE_GRAD_SOU_GRAD = PLIKE_GRAD_SOU_GRAD + sharp_FEMT

           !ewrite(3,*) 'MASS_ELE:', MASS_ELE
           !ewrite(3,*) 'MASS_NORMALISE:', MASS_NORMALISE

           !ewrite(3,*) 'CURVATURE:', CURVATURE
           !ewrite(3,*) 'PLIKE_GRAD_SOU_COEF:', PLIKE_GRAD_SOU_COEF
           !ewrite(3,*) 'PLIKE_GRAD_SOU_GRAD:', PLIKE_GRAD_SOU_GRAD
           !stop 2481

        ELSE

            if ( .not.stress_form ) then

                ! determine the curvature by solving a simple eqn...
                CV_FORCE_X_SUF_TEN=0.0
                CV_FORCE_Y_SUF_TEN=0.0
                CV_FORCE_Z_SUF_TEN=0.0

                U_FORCE_X_SUF_TEN=0.0
                U_FORCE_Y_SUF_TEN=0.0
                U_FORCE_Z_SUF_TEN=0.0
                ! smooth...
                if(.true.) then
                    femtold=0.0
                    DO ELE=1,TOTELE
                        DO CV_ILOC=1,CV_NLOC
                            CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                            DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC
                            femtold(cv_nod)=femtold(cv_nod)+Dx_DIFF_x(DG_CV_NOD,1, ele) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
                        END DO
                    END DO
                    DO SMOOTH_ITS=1,SMOOTH_NITS ! this produces better results but a complex scheme
                        DO CV_NOD=1,CV_NONODS
                            RSUM=0.0
                            RRSUM=0.0
                            DO COUNT=FINACV(CV_NOD),FINACV(CV_NOD+1)-1
                                IF(COLACV(COUNT).LE.CV_NONODS) THEN
                                    RSUM=RSUM+FEMTold(COLACV(COUNT))
                                    RRSUM=RRSUM+1.0
                                ENDIF
                            END DO
                            FEMTOLD2(CV_NOD)=0.5*FEMTold(CV_NOD)+0.5*RSUM/RRSUM
                           !         FEMTOLD(CV_NOD)=0.75*FEMT(CV_NOD)+0.25*RSUM/RRSUM
                           !         FEMTOLD(CV_NOD)=0.9*FEMT(CV_NOD)+0.1*RSUM/RRSUM
                        END DO
                        FEMTOLD=FEMTOLD2
                    END DO
                    DO ELE=1,TOTELE
                        DO CV_ILOC=1,CV_NLOC
                            CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                            DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC
                            Dx_DIFF_x(DG_CV_NOD,1 ,ele) = femtold(cv_nod)
                        END DO
                    END DO

                    femtold=0.0
                    DO ELE=1,TOTELE
                        DO CV_ILOC=1,CV_NLOC
                            CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                            DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC
                            femtold(cv_nod)=femtold(cv_nod)+DY_DIFF_Y(DG_CV_NOD,1,ele) * MASS_ELE(ELE) / MASS_NORMALISE(CV_NOD)
                        END DO
                    END DO
                    DO SMOOTH_ITS=1,SMOOTH_NITS! this produces better results but a complex scheme
                        DO CV_NOD=1,CV_NONODS
                            RSUM=0.0
                            RRSUM=0.0
                            DO COUNT=FINACV(CV_NOD),FINACV(CV_NOD+1)-1
                                IF(COLACV(COUNT).LE.CV_NONODS) THEN
                                    RSUM=RSUM+FEMTold(COLACV(COUNT))
                                    RRSUM=RRSUM+1.0
                                ENDIF
                            END DO
                            FEMTOLD2(CV_NOD)=0.5*FEMTold(CV_NOD)+0.5*RSUM/RRSUM
                           !         FEMTOLD(CV_NOD)=0.75*FEMT(CV_NOD)+0.25*RSUM/RRSUM
                           !         FEMTOLD(CV_NOD)=0.9*FEMT(CV_NOD)+0.1*RSUM/RRSUM
                        END DO
                        FEMTOLD=FEMTOLD2
                    END DO
                    DO ELE=1,TOTELE
                        DO CV_ILOC=1,CV_NLOC
                            CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                            DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC
                            DY_DIFF_Y(DG_CV_NOD,1,ele) = femtold(cv_nod)
                        END DO
                    END DO
                endif


                DO ELE=1,TOTELE
                    DO CV_ILOC=1,CV_NLOC

                        CV_NOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                        DG_CV_NOD=(ELE-1)*CV_NLOC+CV_ILOC

                        !RR =  - SUF_TENSION_COEF * CURVATURE(CV_NOD)

                        if(.true.) then ! make the direction of the force pt towards the smooth gradient
                            ! but keep the magnitude the same.
                            rr =sqrt(SHARP_DTX_ELE(CV_ILOC, 1, ELE)**2+SHARP_DTY_ELE(CV_ILOC, 1, ELE)**2)
                            rr2=sqrt(DTX_ELE(CV_ILOC, 1, ELE)**2+DTY_ELE(CV_ILOC, 1, ELE)**2)
                            grad_c_x=DTX_ELE(CV_ILOC, 1, ELE)*rr/max(1.e-10,rr2)
                            grad_c_y=DTY_ELE(CV_ILOC, 1, ELE)*rr/max(1.e-10,rr2)
                            grad_c_z=DTZ_ELE(CV_ILOC, 1, ELE)*rr/max(1.e-10,rr2)
                        else
                            grad_c_x=SHARP_DTX_ELE(CV_ILOC, 1, ELE)
                            grad_c_y=SHARP_DTY_ELE(CV_ILOC, 1, ELE)
                            grad_c_z=SHARP_DTZ_ELE(CV_ILOC, 1, ELE)
                           !                    grad_c_x=DTX_ELE(CV_ILOC, 1, ELE)
                           !                    grad_c_y=DTY_ELE(CV_ILOC, 1, ELE)
                           !                    grad_c_z=DTZ_ELE(CV_ILOC, 1, ELE)
                        endif

                        !CV_FORCE_X_SUF_TEN(CV_NOD)=CV_FORCE_X_SUF_TEN(CV_NOD)+RR*DTX_ELE(CV_ILOC, 1, ELE)
                        !IF(NDIM.GE.2) CV_FORCE_Y_SUF_TEN(CV_NOD)=CV_FORCE_Y_SUF_TEN(CV_NOD)+RR*DTY_ELE(CV_ILOC, 1, ELE)
                        !IF(NDIM.GE.3) CV_FORCE_Z_SUF_TEN(CV_NOD)=CV_FORCE_Z_SUF_TEN(CV_NOD)+RR*DTZ_ELE(CV_ILOC, 1, ELE)

                        RR=DX_DIFF_X(DG_CV_NOD,1 ,ele)
                        IF(NDIM.GE.2) RR=RR + DY_DIFF_Y(DG_CV_NOD,1 ,ele)
                        IF(NDIM.GE.3) RR=RR + DZ_DIFF_Z(DG_CV_NOD,1, ele)
                        !                 if(rr.ne.0.0) print *,'ele,cv_iloc,rr=',ele,cv_iloc,rr, &
                        !                     sqrt(DTX_ELE(CV_ILOC, 1, ELE)**2+DTX_ELE(CV_ILOC, 1, ELE)**2)
                        !                  RR = - SUF_TENSION_COEF * max(RR,0.0)
                        !                  rr=0.5*16.6666
                        RR = - SUF_TENSION_COEF * RR

                        U_FORCE_X_SUF_TEN(DG_CV_NOD) = U_FORCE_X_SUF_TEN(DG_CV_NOD) + RR * grad_c_x
                        !                  U_FORCE_X_SUF_TEN(DG_CV_NOD) = U_FORCE_X_SUF_TEN(DG_CV_NOD) + RR * DTX_ELE(CV_ILOC, 1, ELE)
                        IF(NDIM.GE.2) U_FORCE_Y_SUF_TEN(DG_CV_NOD) = U_FORCE_Y_SUF_TEN(DG_CV_NOD) + RR *grad_c_y
                        IF(NDIM.GE.3) U_FORCE_Z_SUF_TEN(DG_CV_NOD) = U_FORCE_Z_SUF_TEN(DG_CV_NOD) + RR *grad_c_z

                    END DO
                END DO
                !             stop 121

                DEALLOCATE(DIF_TX, DIF_TY, DIF_TZ)
                DEALLOCATE(DG_CV_NDGLN)
                DEALLOCATE(DX_DIFF_X, DY_DIFF_X, DZ_DIFF_X)
                DEALLOCATE(DX_DIFF_Y, DY_DIFF_Y, DZ_DIFF_Y)
                DEALLOCATE(DX_DIFF_Z, DY_DIFF_Z, DZ_DIFF_Z)

            end if

            !CV_U_FORCE_X_SUF_TEN = CV_FORCE_X_SUF_TEN
            !IF(NDIM.GE.2) CV_U_FORCE_Y_SUF_TEN = CV_FORCE_Y_SUF_TEN
            !IF(NDIM.GE.3) CV_U_FORCE_Z_SUF_TEN = CV_FORCE_Z_SUF_TEN

            if (.false.) then

                ! Convert force to velocity space...
                ALLOCATE(MASS(U_NLOC,U_NLOC))
                ALLOCATE(STORE_MASS(U_NLOC,U_NLOC))
                ALLOCATE(IPIV(U_NLOC))
                ALLOCATE(B_CV_X(CV_NLOC), B_CV_Y(CV_NLOC), B_CV_Z(CV_NLOC))
                ALLOCATE(RHS_U_SHORT_X(U_NLOC), RHS_U_SHORT_Y(U_NLOC), RHS_U_SHORT_Z(U_NLOC))
                ALLOCATE(U_SOL_X(U_NLOC), U_SOL_Y(U_NLOC), U_SOL_Z(U_NLOC))
                DO ELE=1,TOTELE
                    ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
                    CALL DETNLXR_PLUS_U( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, &
                    X_NLOC, CV_NLOC, CV_NGI, &
                    CVFEN, CVFENLX, CVFENLY, CVFENLZ, CVWEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
                    CVFENX_ALL, &
                    U_NLOC, UFENLX, UFENLY, UFENLZ, UFENX_ALL ,&
                    state, "wrapper", StorageIndexes(26))

                    MASS=0.0
                    DO U_ILOC=1,U_NLOC
                        DO U_JLOC=1,U_NLOC
                            NN=0.0
                            DO GI=1,CV_NGI
                                NN = NN + UFEN( U_ILOC, GI ) * UFEN( U_JLOC, GI ) * DETWEI(GI)
                            END DO
                            MASS(U_ILOC,U_JLOC)=MASS(U_ILOC,U_JLOC)+NN
                        END DO
                    END DO

                    DO CV_JLOC=1,CV_NLOC
                        CV_JNOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_JLOC)
                        B_CV_X(CV_JLOC)=CV_FORCE_X_SUF_TEN(CV_JNOD)
                        IF(NDIM.GE.2) B_CV_Y(CV_JLOC)=CV_FORCE_Y_SUF_TEN(CV_JNOD)
                        IF(NDIM.GE.3) B_CV_Z(CV_JLOC)=CV_FORCE_Z_SUF_TEN(CV_JNOD)
                    END DO

                    RHS_U_SHORT_X=0.0
                    RHS_U_SHORT_Y=0.0
                    RHS_U_SHORT_Z=0.0
                    DO U_ILOC=1,U_NLOC
                        DO CV_JLOC=1,CV_NLOC
                            NM=0.0
                            DO GI=1,CV_NGI
                                NM=NM+UFEN( U_ILOC, GI ) * CVFEN( CV_JLOC, GI ) *DETWEI(GI)
                            END DO
                            RHS_U_SHORT_X(U_ILOC)=RHS_U_SHORT_X(U_ILOC)+NM*B_CV_X(CV_JLOC)
                            IF(NDIM.GE.2) RHS_U_SHORT_Y(U_ILOC)=RHS_U_SHORT_Y(U_ILOC)+NM*B_CV_Y(CV_JLOC)
                            IF(NDIM.GE.3) RHS_U_SHORT_Z(U_ILOC)=RHS_U_SHORT_Z(U_ILOC)+NM*B_CV_Z(CV_JLOC)
                        END DO
                    END DO
                    ! Invert mass matrix...
                    ! Solve STORE_MASS *U_SOL_X = RHS_U_SHORT_X
                    ! STORE_MASS is overwritten by lu decomposition which used after the 1st solve.
                    STORE_MASS=MASS
                    GOTDEC = .FALSE.
                    CALL SMLINNGOT( STORE_MASS, U_SOL_X, RHS_U_SHORT_X, U_NLOC, U_NLOC, IPIV,GOTDEC)
                    GOTDEC =.TRUE.
                    IF(NDIM.GE.2) CALL SMLINNGOT( STORE_MASS, U_SOL_Y, RHS_U_SHORT_Y, U_NLOC, U_NLOC, IPIV,GOTDEC)
                    IF(NDIM.GE.3) CALL SMLINNGOT( STORE_MASS, U_SOL_Z, RHS_U_SHORT_Z, U_NLOC, U_NLOC, IPIV,GOTDEC)

                    ! Solve mass matrix systems...
                    DO U_ILOC=1,U_NLOC
                        U_NOD=U_NDGLN((ELE-1)*U_NLOC+U_ILOC)
                        U_FORCE_X_SUF_TEN(U_INOD)=U_SOL_X(U_ILOC)
                        IF(NDIM.GE.2) U_FORCE_Y_SUF_TEN(U_INOD)=U_SOL_Y(U_ILOC)
                        IF(NDIM.GE.3) U_FORCE_Z_SUF_TEN(U_INOD)=U_SOL_Z(U_ILOC)
                    END DO
                END DO

                DEALLOCATE( MASS, STORE_MASS, B_CV_X, B_CV_Y, B_CV_Z, &
                RHS_U_SHORT_X, RHS_U_SHORT_Y, RHS_U_SHORT_Z, &
                U_SOL_X, U_SOL_Y, U_SOL_Z, IPIV)

            end if

        END IF IF_USE_PRESSURE_FORCE


        DEALLOCATE( TDIFFUSION, MASS_NORMALISE, FACE_ELE )
        DEALLOCATE( FEMT, FEMTOLD, MASS_CV, MASS_ELE, &
        XC_CV, YC_CV, ZC_CV, DTX_ELE, DTY_ELE, &
        DTZ_ELE, DTOLDX_ELE, DTOLDY_ELE, DTOLDZ_ELE )
        DEALLOCATE( JCOUNT_KLOC, JCOUNT_KLOC2 )
        DEALLOCATE( CVNORMX, CVNORMY, CVNORMZ )
        DEALLOCATE( COLGPTS, FINDGPTS )
        DEALLOCATE( SNDOTQ, SNDOTQOLD )
        DEALLOCATE( CV_ON_FACE, CVFEM_ON_FACE, &
        U_ON_FACE, UFEM_ON_FACE )
        DEALLOCATE( CV_OTHER_LOC,  U_OTHER_LOC, MAT_OTHER_LOC )
        DEALLOCATE( X_SHARE )
        DEALLOCATE( CVWEIGHT, CVN, CVFEN, &
        CVFENLX, CVFENLY, CVFENLZ )
        DEALLOCATE( CVWEIGHT_SHORT, CVN_SHORT, CVFEN_SHORT, &
        CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT )
        DEALLOCATE( UFEN, UFENLX, UFENLY, UFENLZ )
        DEALLOCATE( SCVFEN, SCVFENSLX, SCVFENSLY, &
        SCVFENLX, SCVFENLY, SCVFENLZ, SCVFEWEIGH )
        DEALLOCATE( SUFEN, SUFENSLX, SUFENSLY, &
        SUFENLX, SUFENLY, SUFENLZ )
        DEALLOCATE( SCVDETWEI, SRA, LOG_ON_BOUND )
        DEALLOCATE( SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
        SBCVFEWEIGH, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
        SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, &
        SBUFENLY, SBUFENLZ, DUMMY_ZERO_NDIM_NDIM )
        DEALLOCATE( CV_SLOC2LOC, U_SLOC2LOC , &
        CV_SLOCLIST, U_SLOCLIST, CV_NEILOC )
        DEALLOCATE( SELE_OVERLAP_SCALE )
        DEALLOCATE( UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE )
        DEALLOCATE( UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2 )
        DEALLOCATE( SUM_CV, UP_WIND_NOD )
        DEALLOCATE( CV_FORCE_X_SUF_TEN, &
        CV_FORCE_Y_SUF_TEN, CV_FORCE_Z_SUF_TEN )
        DEALLOCATE( RDUM, IDUM, RZERO, &
        IZERO, CV_ONE, CURVATURE )

    END SUBROUTINE SURFACE_TENSION_WRAPPER


    subroutine linearise_field( field_in, field_out )
        implicit none
        type( tensor_field ), intent( in ) :: field_in
        real, dimension( :, : ), intent( inout ) :: field_out

        integer, dimension( : ), pointer :: ndglno, cv_nods
        integer :: n, totele, cv_nloc, ncomp, nphase, cv_nonods, ele, cv_iloc, cv_nod
        real, dimension( :, :, : ), allocatable :: field_tmp, field_cv_nod

        ! This sub will linearise a p2 field

        field_out = 0.0

        n = field_in%mesh%shape%degree

        if ( n==2 ) then

            ndglno => get_ndglno( field_in%mesh )

            totele = field_in%mesh%elements
            cv_nloc = field_in%mesh%shape%loc

            ncomp = size( field_in%val, 1 )
            nphase = size( field_in%val, 2 )

            allocate( field_tmp( ncomp, nphase, cv_nonods ) ) ; field_tmp = field_in % val
            allocate( field_cv_nod( ncomp, nphase, cv_nloc ) ) ; field_cv_nod = 0.0

            do ele = 1, totele

                cv_nods => ndglno( ( ele - 1 ) * cv_nloc + 1 : ele * cv_nloc )
                field_cv_nod =  field_tmp( :, :, cv_nods )

                field_cv_nod( :, :, 2 ) = 0.5 * ( field_cv_nod( :, :, 1 ) + field_cv_nod( :, :, 3 ) )
                field_cv_nod( :, :, 4 ) = 0.5 * ( field_cv_nod( :, :, 1 ) + field_cv_nod( :, :, 6 ) )
                field_cv_nod( :, :, 5 ) = 0.5 * ( field_cv_nod( :, :, 3 ) + field_cv_nod( :, :, 6 ) )

                if ( cv_nloc == 10 ) then
                    field_cv_nod( :, :, 7 ) = 0.5 * ( field_cv_nod( :, :, 1 ) + field_cv_nod( :, :, 10 ) )
                    field_cv_nod( :, :, 8 ) = 0.5 * ( field_cv_nod( :, :, 3 ) + field_cv_nod( :, :, 10 ) )
                    field_cv_nod( :, :, 9 ) = 0.5 * ( field_cv_nod( :, :, 6 ) + field_cv_nod( :, :, 10 ) )
                end if

                do cv_iloc = 1, cv_nloc
                    cv_nod = ndglno( ( ele - 1 ) * cv_nloc + cv_iloc )
                    field_out( :, cv_nod ) = field_cv_nod( 1, :, cv_iloc )
                end do

            end do

            deallocate( field_tmp, field_cv_nod )

        end if

        return
    end subroutine linearise_field


end module multiphase_1D_engine
