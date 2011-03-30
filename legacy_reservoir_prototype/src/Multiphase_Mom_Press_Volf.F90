
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

module multiphase_mom_press_volf

  use spact
  use multiphase_1D_engine
  use multiphase_EOS
  use solvers_module
  use matrix_operations
  use shape_functions
  use printout
  use Compositional_Terms
  use fldebug

  IMPLICIT NONE

  private 

  public  :: solve_multiphase_mom_press_volf

contains

  subroutine  solve_multiphase_mom_press_volf( nphase, ncomp, totele, ndim, &
                                ! Nodes et misc
       u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, mat_nloc, &
       cv_snloc, u_snloc, p_snloc, stotel, &
       ncoef, nuabs_coefs, &
                                ! Element types
       u_ele_type, p_ele_type, cv_ele_type, &
       cv_sele_type, u_sele_type, &
                                ! Total time loop and initialisation parameters
       ntime, ntime_dump, nits, nits_internal, &
       nits_flux_lim_volfra, nits_flux_lim_comp, & 
       ndpset, &
       dt, &
                                ! Discretisation parameters
       v_beta, v_theta, &
       v_disopt, &
       v_dg_vel_int_opt, &
       t_beta, t_theta, t_disopt, t_dg_vel_int_opt, lump_eqns, &
       volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, comp_get_theta_flux, &
       opt_vel_upwind_coefs, nopt_vel_upwind_coefs, &
       noit_dim, &
       sat_error_relax2_noit, t_error_relax2_noit, gl_error_relax2_noit, &
       u_error_relax2_noit, p_error_relax2_noit, mass_error_relax2_noit, &
       in_ele_upwind, dg_ele_upwind, &
                                ! Total nodes for different meshes
       domain_length, &
       cv_nonods, u_nonods, cv_pha_nonods, u_pha_nonods, mat_nonods, &
       x_nonods, xu_nonods, &
       u_ndgln, xu_ndgln, cv_ndgln, x_ndgln, p_ndgln, mat_ndgln, &
       u_sndgln, cv_sndgln, p_sndgln, &
                                ! Boundary conditions and surface elements
       wic_vol_bc, wic_d_bc, wic_u_bc, wic_p_bc, wic_t_bc, wic_comp_bc, &
       suf_vol_bc, suf_d_bc, suf_p_bc, suf_t_bc, &
       suf_u_bc, suf_v_bc, suf_w_bc, suf_comp_bc,  &
       suf_u_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob1, suf_v_bc_rob2, &
       suf_w_bc_rob1, suf_w_bc_rob2, suf_comp_bc_rob1, suf_comp_bc_rob2, &
       suf_t_bc_rob1, suf_t_bc_rob2, suf_vol_bc_rob1, suf_vol_bc_rob2, &
                                ! Positions and grid velocities
       x, y, z, nu, nv, nw, ug, vg, wg, &
                                ! Absorption and source terms and coefficients
       uabs_option, uabs_coefs, u_abs_stab, Mobility, &
       u_absorb, v_absorb, comp_absorb, &
       u_source, v_source, comp_source, &
       t_absorb, t_source, Comp_Sum2One, &
                                ! Diffusion parameters
       udiffusion, &
       tdiffusion, &
       comp_diffusion_opt, ncomp_diff_coef, comp_diffusion, comp_diff_coef, &
                                ! Velocities and scalar fields
       u, v, w, &
       den, satura, comp, t, p, cv_p, cv_one,  volfra_pore, Viscosity, &
       perm, &
       uold, vold, wold, denold, saturaold, compold, uden, udenold, deriv, &
       told, pold, cv_pold, nuold, nvold, nwold, &
                                ! EOS terms
       eos_option, &
       eos_coefs, &
       KComp_Sigmoid, K_Comp, alpha_beta, &
                                ! Capillarity pressure terms
       capil_pres_opt, ncapil_pres_coef, capil_pres_coef, &
                                ! Matrices sparsity
       mx_ncolacv, ncolacv, finacv, colacv, midacv, & ! CV multi-phase eqns (e.g. vol frac, temp)
       nlenmcy, mx_ncolmcy, ncolmcy, finmcy, colmcy, midmcy, & ! Force balance plus cty multi-phase eqns
       mxnele, ncolele, finele, colele, & ! Element connectivity 
       mx_ncoldgm_pha, ncoldgm_pha, coldgm_pha, findgm_pha, middgm_pha, & ! Force balance sparsity  
       mx_nct, ncolct, findct, colct, & ! CT sparsity - global cty eqn
       mx_nc, ncolc, findc, colc, & ! C sparsity operating on pressure in force balance
       mx_ncolcmc, ncolcmc, findcmc, colcmc, midcmc, & ! pressure matrix for projection method
       mx_ncolm, ncolm, findm, colm, midm ) ! CV-FEM matrix

    implicit none

    integer, intent( in ) :: nphase, ncomp, totele, ndim, &
         u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, mat_nloc, &
         cv_snloc, u_snloc, p_snloc, stotel, &
         ncoef, nuabs_coefs, &
         u_ele_type, p_ele_type, cv_ele_type, &
         cv_sele_type, u_sele_type, &
         ntime, ntime_dump, nits, nits_internal, &
         nits_flux_lim_volfra, nits_flux_lim_comp, & 
         ndpset 
    real, intent( in ) :: dt
    ! The following need to be changed later in the other subrts as it should be controlled by the 
    ! input files
    real, intent( inout ) :: v_beta, v_theta
    integer, intent( in )  :: v_disopt, v_dg_vel_int_opt
    real, intent( inout ) :: t_beta, t_theta 
    integer, intent( in ) :: t_disopt, t_dg_vel_int_opt
    logical, intent( in ) :: lump_eqns
    logical, intent( inout ) :: volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, comp_get_theta_flux
    integer, intent( in ) :: nopt_vel_upwind_coefs
    real, dimension( nopt_vel_upwind_coefs ), intent( inout ) :: opt_vel_upwind_coefs
    integer, intent( in ) :: noit_dim
    real, dimension( noit_dim ), intent( in ) :: sat_error_relax2_noit, t_error_relax2_noit, gl_error_relax2_noit, &
         u_error_relax2_noit, p_error_relax2_noit, mass_error_relax2_noit
    integer, intent( in ) :: in_ele_upwind, dg_ele_upwind

    real, intent( in )  :: domain_length 
    integer, intent( in ) :: cv_nonods, u_nonods, cv_pha_nonods, u_pha_nonods, mat_nonods, x_nonods, xu_nonods
    integer, dimension( totele * u_nloc ), intent( in )  :: u_ndgln
    integer, dimension( totele * xu_nloc ), intent( in )  :: xu_ndgln
    integer, dimension( totele * cv_nloc ), intent( in )  :: cv_ndgln
    integer, dimension( totele * cv_nloc ), intent( in )  :: x_ndgln
    integer, dimension( totele * p_nloc ), intent( in )  :: p_ndgln
    integer, dimension( totele * mat_nloc ), intent( in )  :: mat_ndgln
    integer, dimension( stotel * u_snloc ), intent( in )  :: u_sndgln
    integer, dimension( stotel * cv_snloc ), intent( in )  :: cv_sndgln

    integer, dimension( stotel * p_snloc ), intent( in )  :: p_sndgln
    integer, dimension( stotel * nphase ), intent( in ) :: wic_vol_bc, wic_d_bc, wic_u_bc, wic_p_bc, wic_t_bc, wic_comp_bc
    real, dimension( stotel * cv_snloc * nphase ), intent( in ) :: suf_vol_bc, suf_d_bc
    real, dimension( stotel * p_snloc * nphase ), intent( in ) :: suf_p_bc
    real, dimension( stotel * cv_snloc * nphase ), intent( in ) :: suf_t_bc
    real, dimension( stotel * u_snloc * nphase ), intent( in ) :: suf_u_bc, suf_v_bc, suf_w_bc
    real, dimension( stotel * cv_snloc * nphase * ncomp ), intent( in ) :: suf_comp_bc
    real, dimension( stotel * u_snloc * nphase ), intent( in ) :: suf_u_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob1, suf_v_bc_rob2, &
         suf_w_bc_rob1, suf_w_bc_rob2
    real, dimension( stotel * cv_snloc * nphase ), intent( in ) :: suf_comp_bc_rob1, suf_comp_bc_rob2, suf_t_bc_rob1, suf_t_bc_rob2, &
         suf_vol_bc_rob1, suf_vol_bc_rob2

    ! Variables bellow may change on this subroutine, however this should be changed later
    real, dimension( x_nonods ), intent( inout ) :: x, y, z
    real, dimension( u_nonods * nphase ), intent( inout ) :: nu, nv, nw, ug, vg, wg
    integer, dimension( nphase ), intent( in ) :: uabs_option
    real, dimension( nphase, nuabs_coefs ), intent( in ) :: uabs_coefs
    real, dimension( mat_nonods, ndim * nphase, ndim * nphase ), intent( inout ) :: u_abs_stab
    real, intent( in ) :: Mobility
    real, dimension( mat_nonods, ndim * nphase, ndim * nphase ), intent( inout ) :: u_absorb
    real, dimension( cv_nonods, nphase, nphase ), intent( inout ) :: v_absorb
    real, dimension( cv_nonods, nphase, nphase ), intent( inout ) :: comp_absorb
    real, dimension( u_pha_nonods ), intent( inout ) :: u_source
    real, dimension( cv_pha_nonods ), intent( inout ) :: v_source, comp_source
    real, dimension( cv_nonods, nphase, nphase ), intent( inout ) :: t_absorb
    real, dimension( cv_pha_nonods ), intent( inout ) :: t_source
    logical, intent( in ) :: Comp_Sum2One

    real, dimension( mat_nonods, ndim, ndim, nphase ), intent( inout ) :: udiffusion, tdiffusion
    integer, intent( in ) :: comp_diffusion_opt, ncomp_diff_coef
    real, dimension( mat_nonods, ndim, ndim, nphase ), intent( inout ) :: comp_diffusion
    real, dimension( ncomp, ncomp_diff_coef, nphase ), intent( in ) :: comp_diff_coef

    ! Variables initialised before but will change during computational time
    real, dimension( u_pha_nonods ), intent( inout ) :: u, v, w
    real, dimension( cv_pha_nonods ), intent( inout ) :: den, satura
    real, dimension( cv_pha_nonods * ncomp ), intent( inout ) :: comp
    real, dimension( cv_pha_nonods ), intent( inout ) :: t
    real, dimension( cv_nonods ), intent( inout ) :: p, cv_p
    real, dimension( cv_pha_nonods ), intent( inout ) :: cv_one
    real, dimension( cv_nonods * nphase ), intent( in ) :: Viscosity
    real, dimension( totele ), intent( inout ) :: volfra_pore
    real, dimension( totele, ndim, ndim ), intent( inout ) :: perm
    real, dimension( u_pha_nonods ), intent( inout ) :: uold, vold, wold
    real, dimension( cv_pha_nonods ), intent( inout ) :: denold, saturaold
    real, dimension( cv_pha_nonods * ncomp ), intent( inout ) :: compold
    real, dimension( cv_pha_nonods ), intent( inout ) :: uden, udenold, deriv, told
    real, dimension( cv_nonods ), intent( inout ) :: pold, cv_pold
    real, dimension( u_pha_nonods ), intent( inout ) :: nuold, nvold, nwold
    integer, dimension( nphase ), intent( in ) :: eos_option
    real, dimension( nphase, ncoef ), intent( in ) :: eos_coefs
    real, dimension( ncomp, nphase, nphase ), intent( inout ) :: K_Comp
    logical, intent( in ) :: KComp_Sigmoid
    real, intent( inout ) :: alpha_beta

    integer, intent( in ) :: capil_pres_opt, ncapil_pres_coef
    real, dimension( ncapil_pres_coef, nphase, nphase ), intent( in ) :: &
         capil_pres_coef

    integer, intent ( in ) :: mx_ncolacv, ncolacv
    integer, dimension( cv_pha_nonods + 1 ), intent (in ) :: finacv
    integer, dimension( mx_ncolacv ), intent (in ) :: colacv
    integer, dimension( cv_pha_nonods ), intent (in ) :: midacv
    integer, intent ( in ) :: nlenmcy, mx_ncolmcy, ncolmcy
    integer, dimension( nlenmcy + 1 ), intent (in ) :: finmcy
    integer, dimension( mx_ncolmcy ), intent (in ) :: colmcy
    integer, dimension( nlenmcy ), intent (in ) :: midmcy
    integer, intent ( in ) :: mxnele, ncolele

    integer, dimension( totele + 1 ), intent (in ) :: finele
    integer, dimension( mxnele ), intent (in ) :: colele
    integer, intent ( in ) :: mx_ncoldgm_pha, ncoldgm_pha
    integer, dimension( mx_ncoldgm_pha ), intent (in ) :: coldgm_pha
    integer, dimension( u_pha_nonods + 1 ), intent (in ) :: findgm_pha
    integer, dimension( u_pha_nonods ), intent (in ) :: middgm_pha
    integer, intent ( in ) :: mx_nct, ncolct
    integer, dimension( cv_nonods + 1 ), intent (in ) :: findct
    integer, dimension( mx_nct ), intent (in ) :: colct
    integer, intent ( in ) :: mx_nc, ncolc
    integer, dimension( u_nonods + 1 ), intent (in ) :: findc
    integer, dimension( mx_nc ), intent (in ) :: colc
    integer, intent ( in ) :: mx_ncolcmc, ncolcmc
    integer, dimension( cv_nonods + 1 ), intent (in ) :: findcmc
    integer, dimension( mx_ncolcmc ), intent (in ) :: colcmc
    integer, dimension( cv_nonods ), intent (in ) :: midcmc
    integer, intent ( in ) :: mx_ncolm, ncolm
    integer, dimension( cv_nonods + 1 ), intent (in ) :: findm
    integer, dimension( mx_ncolm ), intent (in ) :: colm
    ! integer, dimension( ncolm ), intent (in ) :: colm
    integer, dimension( cv_nonods ), intent (in ) :: midm

    ! Local variables
    real :: acctim ! Accumulated time
    integer :: itime, iphase, jphase, its, its2, icomp, icomp2, ncomp2
    real :: dx
    integer, parameter :: izero = 0, rzero = 0.

    ! Absorption terms
    real, dimension( :, :, : ), allocatable :: sigma

    ! For assembling and solving
    real, dimension( : ), allocatable :: rhs, rhs_cv, diag_pres, femt_print, femt_dummy

    ! For output:
    real, dimension( : ), allocatable :: Sat_FEMT, Den_FEMT, Comp_FEMT, SumConc_FEMT, MEAN_PORE_CV
    ! For capillary pressure:
    real, dimension( : ), allocatable :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD

    ! For compositional
    real, dimension( : ), allocatable :: theta_gdiff, V_SOURCE_STORE, V_SOURCE_COMP
    real, dimension( : , : , : , : ), allocatable :: theta_flux, one_m_theta_flux, &
         sum_theta_flux, sum_one_m_theta_flux
    integer :: igot_t2, igot_theta_flux
    integer :: scvngi_theta, cv_ngi, cv_ngi_short, sbcvngi, nface, cv_nodi, IPLIKE_GRAD_SOU


    character( len = 500 ) :: dummy_string_phase, dummy_string_dump, dummy_string_comp, dump_name, file_format
    integer :: output_channel

    ewrite(3,*) 'In solve_multiphase_mom_press_volf'

    allocate( sigma( mat_nonods, ndim * nphase, ndim * nphase ))

    allocate( rhs( cv_nonods ))
    allocate( rhs_cv( cv_nonods ))
    allocate( diag_pres( cv_nonods ))
    allocate( femt_print( cv_nonods * nphase ))
    allocate( femt_dummy( cv_nonods * nphase ))
    allocate( Sat_FEMT( cv_nonods * nphase ))
    allocate( Den_FEMT( cv_nonods * nphase ))
    allocate( Comp_FEMT( cv_nonods * nphase * ncomp ))
    allocate( SumConc_FEMT( cv_nonods * ncomp ))
    allocate( MEAN_PORE_CV( cv_nonods ))

    allocate( V_SOURCE_STORE( cv_nonods * nphase ))
    allocate( V_SOURCE_COMP( cv_nonods * nphase ))
    V_SOURCE_COMP = 0.0

    IPLIKE_GRAD_SOU = 0
    IF( CAPIL_PRES_OPT /= 0 ) IPLIKE_GRAD_SOU = 1
    allocate( PLIKE_GRAD_SOU_COEF( IPLIKE_GRAD_SOU * cv_nonods * nphase ))
    allocate( PLIKE_GRAD_SOU_GRAD( IPLIKE_GRAD_SOU * cv_nonods * nphase ))

    ! Determine scvngi_theta:
    igot_t2 = 0
    igot_theta_flux = 0
    if( ncomp /= 0 )then
       igot_t2 = 1 
       igot_theta_flux = 1
    end if

    call retrieve_ngi(cv_ngi, cv_ngi_short, scvngi_theta, sbcvngi, nface, &
         ndim, cv_ele_type, cv_nloc, u_nloc )
    allocate( theta_flux( totele * igot_theta_flux, cv_nloc, scvngi_theta, nphase ))
    allocate( one_m_theta_flux( totele * igot_theta_flux, cv_nloc, scvngi_theta, nphase ))
    allocate( sum_theta_flux( totele * igot_theta_flux, cv_nloc, scvngi_theta, nphase ))
    allocate( sum_one_m_theta_flux( totele * igot_theta_flux, cv_nloc, scvngi_theta, nphase ))
    allocate( theta_gdiff( cv_nonods * nphase ))

    sum_theta_flux = 1.
    sum_one_m_theta_flux = 0.

    compold = comp

    !             ewrite(3,*)'satura:',satura
    !      stop 383

    DX = DOMAIN_LENGTH / REAL(TOTELE)

    !       print *,'NTIME, ITIME: ',NTIME, ITIME

    ACCTIM = 0.

    Loop_Time: DO ITIME = 1, NTIME

       !          print *,'NITS, ITS: ',NITS, ITS

       ACCTIM = ACCTIM + DT
       UOLD = U
       NU = U
       NUOLD = U
       DENOLD = DEN
       POLD = P
       CV_POLD = CV_P
       SATURAOLD = SATURA
       TOLD = T
       COMPOLD = COMP

       ! Non linear its:
       Loop_ITS: DO ITS = 1, NITS

          CALL CAL_BULK_DENSITY( NPHASE, NCOMP, CV_NONODS, CV_PHA_NONODS, NCOEF, DEN, DERIV, &
               T, CV_P, COMP, EOS_COEFS, EOS_OPTION )

          ! Calculate absorption for momentum eqns    
          CALL CAL_U_ABSORB( MAT_NONODS, CV_NONODS, NPHASE, NDIM, SATURA, TOTELE, CV_NLOC, MAT_NLOC, &
               CV_NDGLN, MAT_NDGLN, &
               NUABS_COEFS, UABS_COEFS, UABS_OPTION, U_ABSORB, VOLFRA_PORE, PERM, MOBILITY, VISCOSITY, &
               X, X_NLOC, X_NONODS, X_NDGLN, &
               OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS )

          IF( CAPIL_PRES_OPT /= 0 ) THEN 
             ! CAPIL_PRES_OPT, NCAPIL_PRES_COEF, CAPIL_PRES_COEF need to be sent down into this sub. 
             CALL CALC_CAPIL_PRES( CAPIL_PRES_OPT, CV_NONODS, NPHASE, NCAPIL_PRES_COEF, &
                  CAPIL_PRES_COEF, IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, SATURA )

          ENDIF

          V_SOURCE_STORE = V_SOURCE + V_SOURCE_COMP 

          ewrite(3,*)'ITIME,NTIME,ITS,NITS:',ITIME,NTIME,ITS,NITS

          ewrite(3,*)'satura:',satura
          ewrite(3,*)'saturaold:',saturaold

          IF( NCOMP <= 1 ) THEN
             VOLFRA_USE_THETA_FLUX = .false.
          else 
             VOLFRA_USE_THETA_FLUX = .true.
          END IF


          CALL FORCE_BAL_CTY_ASSEM_SOLVE( &
               NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
               U_ELE_TYPE, P_ELE_TYPE, &
               U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
               U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN,&
               STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
               U_SNLOC, P_SNLOC, CV_SNLOC, &
               X, Y, Z, U_ABS_STAB, U_ABSORB, U_SOURCE, &
               U, V, W, UOLD, VOLD, WOLD, &
               P, CV_P, DEN, DENOLD, SATURA, SATURAOLD, DERIV, &
               DT, &
               NCOLC, FINDC, COLC, & ! C sparsity - global cty eqn 
               NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, MIDDGM_PHA, &! Force balance sparsity
               NCOLELE, FINELE, COLELE, & ! Element connectivity.
               NCOLCMC, FINDCMC, COLCMC, MIDCMC, & ! pressure matrix for projection method
               NCOLACV, FINACV, COLACV, MIDACV, & ! For CV discretisation method
               NLENMCY, NCOLMCY, FINMCY, COLMCY, MIDMCY, & ! Force balance plus cty multi-phase eqns
               NCOLCT, FINDCT, COLCT, & ! CT sparsity - global cty eqn.
               CV_ELE_TYPE, &
               NU, NV, NW, NUOLD, NVOLD, NWOLD, &
               V_DISOPT, V_DG_VEL_INT_OPT, V_THETA,  &
               SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_P_BC, &
               SUF_U_BC_ROB1, SUF_U_BC_ROB2, SUF_V_BC_ROB1, SUF_V_BC_ROB2,  &
               SUF_W_BC_ROB1, SUF_W_BC_ROB2, &       
               WIC_VOL_BC, WIC_D_BC, WIC_U_BC, WIC_P_BC,  &
               V_SOURCE_STORE, V_ABSORB, VOLFRA_PORE, &
               NCOLM, FINDM, COLM, MIDM, & ! Sparsity for the CV-FEM
               XU_NLOC, XU_NDGLN, &
               UDEN, UDENOLD, UDIFFUSION, NDPSET, &
               OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
               IGOT_THETA_FLUX, SCVNGI_THETA, VOLFRA_USE_THETA_FLUX, &
               SUM_THETA_FLUX, SUM_ONE_M_THETA_FLUX, &
               IN_ELE_UPWIND, DG_ELE_UPWIND, &
               NOIT_DIM, &
               GL_ERROR_RELAX2_NOIT, U_ERROR_RELAX2_NOIT, P_ERROR_RELAX2_NOIT, &
               MASS_ERROR_RELAX2_NOIT, IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD ) 


          CALL CAL_BULK_DENSITY( NPHASE, NCOMP, CV_NONODS, CV_PHA_NONODS, NCOEF, DEN, DERIV, &
               T, CV_P, COMP, EOS_COEFS, EOS_OPTION )

          NU = U
          NV = V
          NW = W

          NUOLD = UOLD
          NVOLD = VOLD
          NWOLD = WOLD

          CALL VOLFRA_ASSEM_SOLVE( &
               NCOLACV, FINACV, COLACV, MIDACV, &
               NCOLCT, FINDCT, COLCT, &
               CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
               CV_ELE_TYPE,  &
               NPHASE,  &
               CV_NLOC, U_NLOC, X_NLOC,  &
               CV_NDGLN, X_NDGLN, U_NDGLN, &
               CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
               X, Y, Z, &
               NU, NV, NW, NUOLD, NVOLD, NWOLD, SATURA, SATURAOLD, DEN, DENOLD, &
               MAT_NLOC,MAT_NDGLN,MAT_NONODS, &
               V_DISOPT, V_DG_VEL_INT_OPT, DT, V_THETA, V_BETA, &
               SUF_VOL_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, &
               WIC_VOL_BC, WIC_D_BC, WIC_U_BC, &
               DERIV, P, &
               V_SOURCE_STORE, V_ABSORB, VOLFRA_PORE, &
               NDIM, &
               NCOLM, FINDM, COLM, MIDM, &
               XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
               OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, & 
               Sat_FEMT, Den_FEMT, &
               IGOT_THETA_FLUX, SCVNGI_THETA, VOLFRA_USE_THETA_FLUX, &
               SUM_THETA_FLUX, SUM_ONE_M_THETA_FLUX, &
               IN_ELE_UPWIND, DG_ELE_UPWIND, &
               NOIT_DIM, &
               SAT_ERROR_RELAX2_NOIT, MASS_ERROR_RELAX2_NOIT, NITS_FLUX_LIM_VOLFRA )


          SUM_THETA_FLUX = 0.0
          SUM_ONE_M_THETA_FLUX = 0.0
          V_SOURCE_COMP = 0.0
          IF( NCOMP <= 1 ) THEN
             NCOMP2 = 0
          ELSE 
             NCOMP2 = NCOMP
          END IF
          ewrite(3,*)'COMP: ',COMP

          Loop_COMPONENTS: DO ICOMP = 1, NCOMP2

             !             COMP_SOURCE = 0.0

             ! Use values from the previous time step DENOLD,SATURAOLD so its easier to converge
             ! ALPHA_BETA is an order 1 acaling coefficient set up in the input file

             CALL CALC_COMP_ABSORB( ICOMP, NCOMP, DT, ALPHA_BETA, &
                                ! NPHASE, CV_NONODS, KCOMP_SIGMOID, K_COMP, DENOLD, SATURAOLD, &
                  NPHASE, CV_NONODS, &
                                !KCOMP_SIGMOID, K_COMP, &
                  .false., K_COMP, &
                  DENOLD, SATURAOLD, &
                  VOLFRA_PORE, COMP_ABSORB, &
                  TOTELE, CV_NLOC, CV_NDGLN )

             IF( KCOMP_SIGMOID ) THEN
                DO CV_NODI = 1, CV_NONODS
                   !                   IF( SATURAOLD( CV_NODI ) > 0.8 ) THEN
                   !                      COMP_ABSORB( CV_NODI, 1, 2 ) = COMP_ABSORB( CV_NODI, 1, 2 ) * &
                   !                           max( 0.01, 5. * ( 1.0 - SATURAOLD( CV_NODI )))
                   !                      COMP_ABSORB( CV_NODI, 2, 2 ) = COMP_ABSORB( CV_NODI, 2, 2 ) * &
                   !                           max( 0.01, 5. * ( 1.0 - SATURAOLD( CV_NODI )))
                   !                   ENDIF
                   !                   IF( SATURAOLD( CV_NODI ) > 0.9 ) THEN
                   !                      COMP_ABSORB( CV_NODI, 1, 2 ) = COMP_ABSORB( CV_NODI, 1, 2 ) * &
                   !                           max( 0.01, 10. * ( 1.0 - SATURAOLD( CV_NODI )))
                   !                      COMP_ABSORB( CV_NODI, 2, 2 ) = COMP_ABSORB( CV_NODI, 2, 2 ) * &
                   !                           max( 0.01, 10. * ( 1.0 - SATURAOLD( CV_NODI )))
                   !                   ENDIF
                   IF( SATURAOLD( CV_NODI ) > 0.95 ) THEN
                      COMP_ABSORB( CV_NODI, 1, 2 ) = COMP_ABSORB( CV_NODI, 1, 2 ) * &
                           max( 0.01, 20. * ( 1.0 - SATURAOLD( CV_NODI )))
                      COMP_ABSORB( CV_NODI, 2, 2 ) = COMP_ABSORB( CV_NODI, 2, 2 ) * &
                           max( 0.01, 20. * ( 1.0 - SATURAOLD( CV_NODI )))
                   ENDIF
                END DO
             END IF

             ! Calculate the diffusion COMP_DIFFUSION...        
             CALL CALC_COMP_DIF( NDIM, NPHASE, COMP_DIFFUSION_OPT, MAT_NONODS, &
                  TOTELE, MAT_NLOC, CV_NLOC, U_NLOC, CV_SNLOC, U_SNLOC, &
                  COMP_DIFFUSION, NCOMP_DIFF_COEF, COMP_DIFF_COEF( ICOMP, :, : ), &
                  X_NONODS, X, Y, Z, NU, NV, NW, U_NONODS, CV_NONODS, &
                  MAT_NDGLN, U_NDGLN, X_NDGLN, &
                  U_ELE_TYPE, P_ELE_TYPE ) 


             ! bear in mind there are 3 internal iterations
             Loop_ITS2: DO ITS2 = 1, nits_internal

                COMP_GET_THETA_FLUX = .TRUE. ! This will be set up in the input file
                COMP_USE_THETA_FLUX = .FALSE.

                CALL INTENERGE_ASSEM_SOLVE(  &
                     NCOLACV, FINACV, COLACV, MIDACV, & ! CV sparcity pattern matrix
                     NCOLCT, FINDCT, COLCT, &
                     CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
                     U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE,  &
                     NPHASE,  &
                     CV_NLOC, U_NLOC, X_NLOC,  &
                     CV_NDGLN, X_NDGLN, U_NDGLN, &
                     CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
                     X, Y, Z, &
                     NU, NV, NW, NUOLD, NVOLD, NWOLD, &
                     UG, VG, WG, &
                     COMP(( ICOMP - 1 ) * NPHASE * CV_NONODS + 1 ), COMPOLD(( ICOMP - 1 ) * NPHASE * CV_NONODS + 1 ), &
                     DEN, DENOLD,  &
                     MAT_NLOC, MAT_NDGLN, MAT_NONODS, COMP_DIFFUSION, &
                     V_DISOPT, V_DG_VEL_INT_OPT, DT, V_THETA, V_BETA, &
                     SUF_COMP_BC( 1 + STOTEL * CV_SNLOC * NPHASE *( ICOMP - 1 )), SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, &
                     SUF_COMP_BC_ROB1, SUF_COMP_BC_ROB2,  &
                     WIC_COMP_BC, WIC_D_BC, WIC_U_BC, &
                     DERIV, P, &
                     COMP_SOURCE, COMP_ABSORB, VOLFRA_PORE, &
                     NDIM,  &
                     NCOLM, FINDM, COLM, MIDM, &
                     XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, LUMP_EQNS, &
                                !  OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, Comp_FEMT(( ICOMP - 1 ) * NPHASE * CV_NONODS + 1 ), Den_FEMT, & 
                     OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
                     Comp_FEMT(( ICOMP - 1 ) * NPHASE * CV_NONODS + 1 : ICOMP * NPHASE * CV_NONODS ), Den_FEMT, & 
                     IGOT_T2, SATURA, SATURAOLD, IGOT_THETA_FLUX, SCVNGI_THETA, COMP_GET_THETA_FLUX, COMP_USE_THETA_FLUX, &
                     THETA_FLUX, ONE_M_THETA_FLUX, THETA_GDIFF, &
                     SUF_VOL_BC, SUF_VOL_BC_ROB1, SUF_VOL_BC_ROB2, WIC_VOL_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
                     NOIT_DIM, &
                     T_ERROR_RELAX2_NOIT, MASS_ERROR_RELAX2_NOIT, NITS_FLUX_LIM_COMP, &
                     MEAN_PORE_CV )

                do iphase = 1, nphase
                   ewrite(3,*) 'comp:', icomp, iphase, comp(( ICOMP - 1 ) * NPHASE * CV_NONODS + (iphase - 1) * cv_nonods + 1 )
                end do

             END DO Loop_ITS2

             SUM_THETA_FLUX = SUM_THETA_FLUX + THETA_FLUX
             SUM_ONE_M_THETA_FLUX = SUM_ONE_M_THETA_FLUX + ONE_M_THETA_FLUX

             V_SOURCE_COMP = V_SOURCE_COMP + THETA_GDIFF
             DO IPHASE = 1, NPHASE
                DO JPHASE = 1, NPHASE
                   DO CV_NODI = 1, CV_NONODS
                      V_SOURCE_COMP( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) = V_SOURCE_COMP( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) &
                           -COMP_ABSORB( CV_NODI, IPHASE, JPHASE ) &
                           * COMP( CV_NODI + ( JPHASE - 1 ) * CV_NONODS + ( ICOMP - 1 ) * NPHASE * CV_NONODS )
                   END DO
                END DO
             END DO


          END DO Loop_COMPONENTS

          IF( COMP_SUM2ONE .AND. ( NCOMP2 > 1 ) ) THEN 
             ! make sure the composition sums to 1.0 by putting constraint into V_SOURCE_COMP. 
             CALL CAL_COMP_SUM2ONE_SOU( V_SOURCE_COMP, CV_NONODS, NPHASE, NCOMP2, DT, ITS, NITS,  &  
                  MEAN_PORE_CV, SATURA, SATURAOLD, DEN, DENOLD, COMP, COMPOLD ) 
          ENDIF

          ewrite(3,*)'Finished VOLFRA_ASSEM_SOLVE ITS,nits,ITIME,NTIME:',ITS,nits,ITIME,NTIME

       END DO Loop_ITS

       Conditional_TIMDUMP: if( ( mod( itime, ntime_dump ) == 0 ) .or. ( itime == 1 ) .or. &
            ( itime == ntime )) then

          ! output the saturation, density and velocity for each phase
          Phase_output_loop: do iphase = 1,nphase

             ! form string for phase number
             call write_integer_to_string(iphase,dummy_string_phase,len(dummy_string_phase))  

             ! form string for dump number
             call write_integer_to_string(itime,dummy_string_dump,len(dummy_string_dump))  

             ! form the output file name for saturation --------------------------
             dump_name = 'saturation_FEM_phase_'//trim(dummy_string_phase)//'.d.'//trim(dummy_string_dump)

             ! open the output file for saturation   
             output_channel = 1
             file_format = 'formatted'
             call open_output_file(output_channel,dump_name,len(dump_name),file_format)

             ! output the fem interpolation of the CV solution for saturation of this phase
             call output_fem_sol_of_cv(output_channel,totele,cv_nonods,x_nonods,nphase,cv_nloc,x_nloc,cv_ndgln, &
                  x_ndgln,x,Sat_FEMT( ((iphase-1)*cv_nonods+1):((iphase-1)*cv_nonods+cv_nonods) ))

             ! close the file for saturation   
             close(output_channel)  


             ! form the output file name for density --------------------------
             dump_name = 'density_FEM_phase_'//trim(dummy_string_phase)//'.d.'//trim(dummy_string_dump)

             ! open the output file for density   
             output_channel = 1
             file_format = 'formatted'
             call open_output_file(output_channel,dump_name,len(dump_name),file_format)

             ! output the fem interpolation of the CV solution for density of this phase
             call output_fem_sol_of_cv(output_channel,totele,cv_nonods,x_nonods,nphase,cv_nloc,x_nloc,cv_ndgln, &
                  x_ndgln,x,DEN_FEMT( ((iphase-1)*cv_nonods+1):((iphase-1)*cv_nonods+cv_nonods) ))

             ! close the file for density   
             close(output_channel)  

             ! form the output file name for velocity --------------------------
             dump_name = 'velocity_phase_'//trim(dummy_string_phase)//'.d.'//trim(dummy_string_dump)

             ! open the output file for velocity   
             output_channel = 1
             file_format = 'formatted'
             call open_output_file(output_channel,dump_name,len(dump_name),file_format)

             ! output the fem velocity
             call printing_veloc_field( output_channel, totele, xu_nonods, xu_nloc, xu_ndgln, u_nloc, u_ndgln, &
                  x, u_nonods, u_nonods * nphase, u, iphase )

             ! close the file for velocity   
             close(output_channel)

          end do Phase_output_loop


          Loop_Comp_Print: do icomp = 1, ncomp

             Loop_Phase_Print: do iphase = 1, nphase

                ! form string for phase number
                call write_integer_to_string( iphase, dummy_string_phase, len( dummy_string_phase )) 

                ! form string for component number
                call write_integer_to_string( icomp, dummy_string_comp, len( dummy_string_comp ))  

                ! form string for dump number
                call write_integer_to_string( itime, dummy_string_dump, len( dummy_string_dump ))  

                ! form the output file name for saturation --------------------------
                dump_name = 'MassFraction_FEM_Phase_'//trim( dummy_string_phase )// '_Comp_' // &
                     trim( dummy_string_comp ) // '.d.' // trim( dummy_string_dump )

                ! open the output file for Mass/Molar Fraction   
                output_channel = 1
                file_format = 'formatted'
                call open_output_file( output_channel, dump_name, len( dump_name ), file_format )

                ! output the fem interpolation of the CV solution for Mass Fraction of this phase and this component
                call output_fem_sol_of_cv( output_channel, totele, cv_nonods, x_nonods, nphase, cv_nloc, x_nloc, cv_ndgln, &
                     x_ndgln, x, Comp_FEMT(( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1 : &
                     ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods ))

                ! close the file   
                close( output_channel )  

             end do Loop_Phase_Print

          end do Loop_Comp_Print



          Loop_Comp_Print2: do icomp = 1, ncomp

             ! form string for component number
             call write_integer_to_string( icomp, dummy_string_comp, len( dummy_string_comp )) 

             ! form string for dump number
             call write_integer_to_string( itime, dummy_string_dump, len( dummy_string_dump ))  

             ! form the output file name for saturation --------------------------
             dump_name = 'Conc_FEM_Comp_' // trim( dummy_string_comp ) // '.d.' // trim( dummy_string_dump )

             ! open the output file for Mass/Molar Fraction   
             output_channel = 1
             file_format = 'formatted'
             call open_output_file( output_channel, dump_name, len( dump_name ), file_format )

             SumConc_FEMT = 0.

             Loop_Phase_Print2: do iphase = 1, nphase

                SumConc_FEMT( ( icomp - 1 ) * cv_nonods + 1 : ( icomp - 1 ) * cv_nonods + cv_nonods ) = &
                     SumConc_FEMT( ( icomp - 1 ) * cv_nonods + 1 : ( icomp - 1 ) * cv_nonods + cv_nonods ) + &
                     Comp_FEMT(( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1 : &
                     ( icomp - 1 ) * nphase * cv_nonods + iphase * cv_nonods )  * &
                     Sat_FEMT(( ( iphase - 1 ) * cv_nonods + 1 ) : (( iphase - 1 ) * cv_nonods + cv_nonods ))

             end do Loop_Phase_Print2

             ! output the fem interpolation of the CV solution for Mass Fraction of this phase and this component
             call output_fem_sol_of_cv( output_channel, totele, cv_nonods, x_nonods, ncomp, cv_nloc, x_nloc, cv_ndgln, &
                  x_ndgln, x, SumConc_FEMT(( ( icomp - 1 ) * cv_nonods + 1 ) : (( icomp - 1 ) * cv_nonods + cv_nonods )))

             ! close the file   
             close( output_channel )  

          end do Loop_Comp_Print2


          ! output the global pressure --------------------------
          dump_name = 'pressure'//'.d.'//trim(dummy_string_dump)

          ! open the output file for pressure   
          output_channel = 1
          file_format = 'formatted'
          call open_output_file(output_channel,dump_name,len(dump_name),file_format) 

          ! output the pressure      
          call printing_field_array( output_channel, totele, cv_nonods, x_nonods, x_nloc, x_ndgln, cv_nloc, cv_ndgln, &
               x, cv_nonods, cv_p, 1 )

          ! close the file for pressure   
          close(output_channel)   

       end if Conditional_TIMDUMP

    END DO Loop_Time

    ewrite(3,*) 'Leaving solve_multiphase_mom_press_volf'



    deallocate( sigma )
    deallocate( rhs )
    deallocate( rhs_cv )
    deallocate( diag_pres )
    deallocate( femt_print )
    deallocate( femt_dummy )
    deallocate( Sat_FEMT )
    deallocate( Den_FEMT )
    deallocate( Comp_FEMT )

  end subroutine solve_multiphase_mom_press_volf


end module multiphase_mom_press_volf
