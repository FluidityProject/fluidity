
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
  use matrix_operations
  use shape_functions
  use printout
  use Compositional_Terms
  use copy_outof_into_state
  use write_state_module
  use diagnostic_variables
  use diagnostic_fields_wrapper
  use diagnostic_fields_new, only : &
         & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
         & check_diagnostic_dependencies
  
  use fldebug
  use state_module
  use spud
  use signal_vars
  use populate_state_module
  use mapping_for_ocvfem

  IMPLICIT NONE

  private 

  public  :: solve_multiphase_mom_press_volf

contains

  subroutine  solve_multiphase_mom_press_volf( state, nphase, ncomp, totele, ndim, &
                                ! Nodes et misc
       u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, mat_nloc, &
       cv_snloc, u_snloc, p_snloc, stotel, &
                                ! Element types
       u_ele_type, p_ele_type, cv_ele_type, &
       cv_sele_type, u_sele_type, &
                                ! Total time loop and initialisation parameters
       ntime_dump, nits, nits_internal, dump_no, &
       nits_flux_lim_volfra, nits_flux_lim_comp, & 
       ndpset, &
                                ! Discretisation parameters
       v_beta, v_theta, &
       v_disopt, &
       v_dg_vel_int_opt, &
       t_beta, t_theta, t_disopt, t_dg_vel_int_opt, lump_eqns, &
       volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, comp_get_theta_flux, &
       opt_vel_upwind_coefs, nopt_vel_upwind_coefs, &
       noit_dim, &
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
       x, y, z, xu, yu, zu, nu, nv, nw, ug, vg, wg, &
                                ! Absorption and source terms and coefficients
       u_abs_stab, Mobility, &
       u_absorb, v_absorb, comp_absorb, &
       u_source, v_source, comp_source, &
       t_absorb, t_source, &
                                ! Diffusion parameters
       udiffusion, &
       tdiffusion, &
       comp_diffusion_opt, ncomp_diff_coef, comp_diffusion, comp_diff_coef, &
                                ! Velocities and scalar fields
       u, v, w, &
       den, satura, comp, t, p, cv_p, volfra_pore, &
       perm, &
       uold, vold, wold, denold, saturaold, compold, uden, udenold, deriv, &
       told, pold, cv_pold, nuold, nvold, nwold, &
                                ! EOS terms
       K_Comp, alpha_beta, &
                                ! Matrices sparsity
       mx_ncolacv, ncolacv, finacv, colacv, midacv, & ! CV multi-phase eqns (e.g. vol frac, temp)
       nlenmcy, mx_ncolmcy, ncolmcy, finmcy, colmcy, midmcy, & ! Force balance plus cty multi-phase eqns
       mxnele, ncolele, finele, colele, & ! Element connectivity 
       mx_ncoldgm_pha, ncoldgm_pha, coldgm_pha, findgm_pha, middgm_pha, & ! Force balance sparsity  
       mx_nct, ncolct, findct, colct, & ! CT sparsity - global cty eqn
       mx_nc, ncolc, findc, colc, & ! C sparsity operating on pressure in force balance
       mx_ncolcmc, ncolcmc, findcmc, colcmc, midcmc, & ! pressure matrix for projection method
       mx_ncolm, ncolm, findm, colm, midm, & ! CV-FEM matrix
       have_temperature_fields, scale_momentum_by_volume_fraction ,cv_one, nits_flux_lim_t, t_use_theta_flux, &
       t_get_theta_flux )
    implicit none

    type( state_type ), dimension( : ), intent( inout ) :: state

    integer, intent( in ) :: nphase, ncomp, totele, ndim, &
         u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, mat_nloc, &
         cv_snloc, u_snloc, p_snloc, stotel, &
         u_ele_type, p_ele_type, cv_ele_type, &
         cv_sele_type, u_sele_type, &
         ntime_dump, nits, nits_internal, &
         nits_flux_lim_volfra, nits_flux_lim_comp, & 
         ndpset 
    real :: dt, current_time
    integer :: dump_no

    ! The following need to be changed later in the other subrts as it should be controlled by the 
    ! input files
    real, intent( inout ) :: v_beta, v_theta
    integer, intent( in )  :: v_disopt, v_dg_vel_int_opt
    real, intent( inout ) :: t_beta, t_theta 
    integer, intent( in ) :: t_disopt, t_dg_vel_int_opt
    logical, intent( in ) :: lump_eqns
    logical, intent( inout ) :: volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, comp_get_theta_flux
    logical, intent( inout ) :: t_use_theta_flux, t_get_theta_flux    
    integer, intent( in ) :: nopt_vel_upwind_coefs
    real, dimension( nopt_vel_upwind_coefs ), intent( inout ) :: opt_vel_upwind_coefs
    integer, intent( in ) :: noit_dim
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
    real, dimension( stotel * u_snloc * nphase ), intent( in ) :: suf_u_bc_rob1, suf_u_bc_rob2,&
         suf_v_bc_rob1, suf_v_bc_rob2, suf_w_bc_rob1, suf_w_bc_rob2
    real, dimension( stotel * cv_snloc * nphase ), intent( in ) :: suf_comp_bc_rob1, &
         suf_comp_bc_rob2, suf_t_bc_rob1, suf_t_bc_rob2, suf_vol_bc_rob1, suf_vol_bc_rob2

    ! Variables below may change on this subroutine, however this should be changed later
    real, dimension( x_nonods ), intent( inout ) :: x, y, z
    real, dimension( xu_nonods ), intent( inout ) :: xu, yu, zu
    real, dimension( u_nonods * nphase ), intent( inout ) :: nu, nv, nw, ug, vg, wg
    real, dimension( mat_nonods, ndim * nphase, ndim * nphase ), intent( inout ) :: u_abs_stab
    real, intent( in ) :: Mobility
    real, dimension( mat_nonods, ndim * nphase, ndim * nphase ), intent( inout ) :: u_absorb
    real, dimension( cv_nonods, nphase, nphase ), intent( inout ) :: v_absorb
    real, dimension( cv_nonods, nphase, nphase ), intent( inout ) :: comp_absorb
    real, dimension( u_nonods * nphase * ndim ), intent( inout ) :: u_source
    real, dimension( cv_pha_nonods ), intent( inout ) :: v_source, comp_source
    real, dimension( cv_nonods, nphase, nphase ), intent( inout ) :: t_absorb
    real, dimension( cv_pha_nonods ), intent( inout ) :: t_source

    real, dimension( mat_nonods, ndim, ndim, nphase ), intent( inout ) :: udiffusion, tdiffusion
    integer, intent( in ) :: comp_diffusion_opt, ncomp_diff_coef
    real, dimension( mat_nonods, ndim, ndim, nphase ), intent( inout ) :: comp_diffusion
    real, dimension( ncomp, ncomp_diff_coef, nphase ), intent( in ) :: comp_diff_coef

    ! Variables initialised before but will change during computational time
    real, dimension( u_nonods * nphase ), intent( inout ) :: u, v, w
    real, dimension( cv_pha_nonods ), intent( inout ) :: den, satura
    real, dimension( cv_pha_nonods * ncomp ), intent( inout ) :: comp
    real, dimension( cv_pha_nonods ), intent( inout ) :: t
    real, dimension( cv_nonods ), intent( inout ) :: p, cv_p
    real, dimension( totele ), intent( inout ) :: volfra_pore
    real, dimension( totele, ndim, ndim ), intent( inout ) :: perm
    real, dimension( u_nonods * nphase ), intent( inout ) :: uold, vold, wold
    real, dimension( cv_pha_nonods ), intent( inout ) :: denold, saturaold
    real, dimension( cv_pha_nonods * ncomp ), intent( inout ) :: compold
    real, dimension( cv_pha_nonods ), intent( inout ) :: uden, udenold, deriv, told
    real, dimension( cv_nonods ), intent( inout ) :: pold, cv_pold
    real, dimension( u_nonods * nphase ), intent( inout ) :: nuold, nvold, nwold

    real, dimension( ncomp, nphase, nphase ), intent( inout ) :: K_Comp
    real, intent( inout ) :: alpha_beta

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
    integer, dimension( cv_nonods ), intent (in ) :: midm

    logical, intent(in) :: have_temperature_fields, scale_momentum_by_volume_fraction 
    real, dimension( cv_nonods * nphase ), intent( inout ) :: cv_one
    integer, intent(in) :: nits_flux_lim_t

    ! Local variables
    real :: acctim ! Accumulated time
    integer :: itime, iphase, jphase, its, its2, icomp, icomp2, ncomp2
    integer :: nstates
    real :: dx
    integer, parameter :: izero = 0
    real, parameter :: rzero = 0.

    ! Absorption terms and dg velocity
    real, dimension( :, :, : ), allocatable :: sigma, velocity_dg

    ! For assembling and solving
    real, dimension( : ), allocatable :: rhs, rhs_cv, diag_pres, femt_print, femt_dummy

    ! For output:
    real, dimension( : ), allocatable :: Sat_FEMT, Den_FEMT, Comp_FEMT, &
         SumConc_FEMT, MEAN_PORE_CV
    ! For capillary pressure:
    real, dimension( : ), allocatable :: PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD

    ! For compositional
    real, dimension( : ), allocatable :: theta_gdiff, V_SOURCE_STORE, V_SOURCE_COMP
    real, dimension( : , : , : , : ), allocatable :: theta_flux, one_m_theta_flux, &
         sum_theta_flux, sum_one_m_theta_flux
    integer :: igot_t2, igot_theta_flux
    integer :: scvngi_theta, cv_ngi, cv_ngi_short, sbcvngi, nface, cv_nodi, IPLIKE_GRAD_SOU

    real, dimension( : ), allocatable :: T_FEMT, mass_ele, dummy_ele


    real, dimension(cv_nonods * nphase * ndim) :: u_source_cv

    real :: finish_time
    integer :: dump_period_in_timesteps
    integer :: final_timestep
    integer :: stat, velocity_max_iterations, saturation_max_iterations

    character( len = 500 ) :: dummy_string_phase, dummy_string_dump, &
         dummy_string_comp, dump_name, file_format
    integer :: output_channel
    real :: norm_satura1, norm_satura2
    logical :: solve_force_balance, solve_saturation

    type(scalar_field), pointer :: pressure, temperature

    allocate( sigma( mat_nonods, ndim * nphase, ndim * nphase ))

    allocate( rhs( cv_nonods ))
    allocate( velocity_dg( cv_nloc*totele,nphase,ndim ))
    allocate( rhs_cv( cv_nonods ))
    allocate( diag_pres( cv_nonods ))
    allocate( femt_print( cv_nonods * nphase ))
    allocate( femt_dummy( cv_nonods * nphase ))
    allocate( Sat_FEMT( cv_nonods * nphase ))
    allocate( T_FEMT( cv_nonods * nphase ))
    allocate( Den_FEMT( cv_nonods * nphase ))
    allocate( Comp_FEMT( cv_nonods * nphase * ncomp ))
    allocate( SumConc_FEMT( cv_nonods * ncomp ))
    allocate( MEAN_PORE_CV( cv_nonods ))

    allocate( V_SOURCE_STORE( cv_nonods * nphase )) ; V_SOURCE_STORE = 0.0
    allocate( V_SOURCE_COMP( cv_nonods * nphase )) ; V_SOURCE_COMP = 0.0

    allocate( PLIKE_GRAD_SOU_GRAD( cv_nonods * nphase )) ; PLIKE_GRAD_SOU_GRAD = 0.0
    allocate( PLIKE_GRAD_SOU_COEF( cv_nonods * nphase )) ; PLIKE_GRAD_SOU_COEF = 0.0

    allocate( mass_ele( totele ) ) ; mass_ele = 0.
    allocate( dummy_ele( totele ) ) ; dummy_ele = 0.

    ! Determine scvngi_theta:
    igot_t2 = 0
    igot_theta_flux = 0
    if( ncomp /= 0 )then
       igot_t2 = 1 
       igot_theta_flux = 1
    end if

    call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
         cv_ngi, cv_ngi_short, scvngi_theta, sbcvngi, nface, .false. )

    allocate( theta_flux( totele * igot_theta_flux, cv_nloc, scvngi_theta, nphase ))
    allocate( one_m_theta_flux( totele * igot_theta_flux, cv_nloc, scvngi_theta, nphase ))
    allocate( sum_theta_flux( totele * igot_theta_flux, cv_nloc, scvngi_theta, nphase ))
    allocate( sum_one_m_theta_flux( totele * igot_theta_flux, cv_nloc, scvngi_theta, nphase ))
    allocate( theta_gdiff( cv_nonods * nphase ))

    sum_theta_flux = 1.
    sum_one_m_theta_flux = 0.

    compold = comp

    DX = DOMAIN_LENGTH / REAL(TOTELE)

!!!
!!! Check later for a more well-defined way to set up an advection problem
!!!

    call get_option( "/material_phase[0]/vector_field::Velocity/prognostic/solver/max_iterations", &
         velocity_max_iterations,  default =  500 )

    call get_option( "/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/solver/max_iterations", &
         saturation_max_iterations,  default =  500 )

    solve_force_balance = .false.
    if( velocity_max_iterations /= 0 ) solve_force_balance = .true.
    solve_saturation = .false.
    if( saturation_max_iterations /= 0 ) solve_saturation = .true. 

    !ewrite(3,*) 'suf_u_bc:'
    !do iphase = 1, nphase
    !   ewrite(3,*) suf_u_bc( ( iphase - 1 ) * stotel * u_snloc  + 1 : &
    !        iphase * stotel * u_snloc  )
    !end do

    !ewrite(3,*)'suf_comp_bc:', stotel, cv_snloc
    !do icomp = 1, ncomp
    !   do iphase = 1, nphase
    !      ewrite(3,*) icomp, iphase, suf_comp_bc( ( icomp - 1 ) * nphase * stotel * cv_snloc + &
    !           ( iphase - 1 ) * stotel * cv_snloc + 1 : &
    !           ( icomp - 1 ) * nphase * stotel * cv_snloc + &
    !           ( iphase - 1 ) * stotel * cv_snloc + stotel * cv_snloc )
    !      ewrite(3,*)' '
    !   end do
    !end do

    call get_option("/timestepping/current_time", acctim)
    call get_option("/timestepping/timestep", dt)
    call get_option("/timestepping/finish_time", finish_time)
    call get_option("/io/dump_period_in_timesteps/constant", dump_period_in_timesteps, default=1)

    nstates = option_count("/material_phase")

    itime = 0

    Loop_Time: DO 

       itime = itime + 1

       ACCTIM = ACCTIM + DT
       call set_option("/timestepping/current_time", acctim)

       if ( ACCTIM > finish_time ) then 

          ewrite(1, *) "Passed final time"

          exit Loop_Time

       end if

       call get_option("/timestepping/final_timestep", final_timestep, stat)

       if(stat == SPUD_NO_ERROR) then

          if(itime > final_timestep) then

             ewrite(1, *) "Passed final timestep"

             exit Loop_Time

          end if

       end if

       UOLD = U
       NU = U
       NUOLD = U

       VOLD = V
       NV = V
       NVOLD = V

       WOLD = W
       NW = W
       NWOLD = W

       DENOLD = DEN
       POLD = P
       CV_POLD = CV_P
       SATURAOLD = SATURA
       TOLD = T
       COMPOLD = COMP

       ! update state memory
       if (have_temperature_fields) then
          do iphase = 1, nphase
             temperature => extract_scalar_field( state( iphase ), 'Temperature' )
             temperature % val = t( 1+(iphase-1)*cv_nonods : iphase*cv_nonods )
          end do
       end if
       pressure => extract_scalar_field( state( 1 ), 'Pressure' )
       pressure % val = p

       ewrite(1,*)' New Time Step:', itime

       ! evaluate prescribed fields at time = current_time+dt
       call set_prescribed_field_values(state, exclude_interpolated=.true., &
            exclude_nonreprescribed=.true., time=ACCTIM)

       ! Non linear its:
       Loop_ITS: DO ITS = 1, NITS
          ewrite(1,*)' New Non-Linear Iteration:', its

          CALL Calculate_Phase_Component_Densities( state, DEN, DERIV )

          if ( its == 1 ) DENOLD = DEN

          ! solve temperature fields if found in input and prognostic
          solve_temp: if (have_temperature_fields .and. &
               have_option('/material_phase[0]/scalar_field::Temperature/prognostic')) then

             NU = U
             NV = V
             NW = W

             call INTENERGE_ASSEM_SOLVE( state, &
                  NCOLACV, FINACV, COLACV, MIDACV, & 
                  NCOLCT, FINDCT, COLCT, &
                  CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
                  U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE,  &
                  NPHASE,  &
                  CV_NLOC, U_NLOC, X_NLOC,  &
                  CV_NDGLN, X_NDGLN, U_NDGLN, &
                  CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
                  X, Y, Z, &
                  NU, NV, NW, NUOLD, NVOLD, NWOLD, UG, VG, WG, &
                  T, TOLD, &
                  DEN, DENOLD, &
                  MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, &
                  T_DISOPT, T_DG_VEL_INT_OPT, DT, T_THETA, T_BETA, &
                  SUF_T_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, &
                  SUF_T_BC_ROB1, SUF_T_BC_ROB2,  &
                  WIC_T_BC, WIC_D_BC, WIC_U_BC, &
                  DERIV, CV_P, &
                  T_SOURCE, T_ABSORB, VOLFRA_PORE, &
                  NDIM,  &
                  NCOLM, FINDM, COLM, MIDM, &
                  XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, LUMP_EQNS, &
                  OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
                  T_FEMT, CV_ONE, &
                  IGOT_T2,T,TOLD, IGOT_THETA_FLUX, SCVNGI_THETA, T_GET_THETA_FLUX, T_USE_THETA_FLUX, &
                  T, T, T, &
                  SUF_T_BC, SUF_T_BC_ROB1, SUF_T_BC_ROB2, WIC_T_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
                  NOIT_DIM, &
                  nits_flux_lim_t, &
                  MEAN_PORE_CV, &
                  option_path = '/material_phase[0]/scalar_field::Temperature', &
                  mass_ele_transp = dummy_ele, &
                  thermal = .true. )

             ! update state memory
             do iphase = 1, nphase
                temperature => extract_scalar_field( state( iphase ), 'Temperature' )
                temperature % val = t( 1+(iphase-1)*cv_nonods : iphase*cv_nonods )
             end do

             CALL Calculate_Phase_Component_Densities( state, DEN, DERIV )

             if (SIG_INT) exit Loop_ITS

          end if solve_temp

          if (SIG_INT) exit Loop_ITS

          if (solve_force_balance) then
             ! Calculate absorption for momentum eqns    
             CALL calculate_absorption( MAT_NONODS, CV_NONODS, NPHASE, NDIM, &
                  SATURA, TOTELE, CV_NLOC, MAT_NLOC, &
                  CV_NDGLN, MAT_NDGLN, &
                  U_ABSORB, PERM, MOBILITY, &
                  OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS )
          end if

          if (SIG_INT) exit Loop_ITS

          IF( have_option("/material_phase[0]/multiphase_properties/capillary_pressure") ) THEN
             CALL calculate_capillary_pressure( state, CV_NONODS, NPHASE, PLIKE_GRAD_SOU_GRAD, SATURA )
          END IF

          V_SOURCE_STORE = V_SOURCE + V_SOURCE_COMP

          IF( NCOMP <= 1 ) THEN
             VOLFRA_USE_THETA_FLUX = .false.
          ELSE 
             VOLFRA_USE_THETA_FLUX = .true.
          END IF

          if(.not. have_option("/material_phase[0]/multiphase_properties/relperm_type")) then
             uden = den
             udenold = denold
          end if

          if (solve_force_balance) then
             NU = U
             NV = V
             NW = W

             NUOLD = UOLD
             NVOLD = VOLD
             NWOLD = WOLD

             ! this calculates u_source_cv - it is the buoyancy term
             ! as the name suggests it's a cv source term for u
             call calculate_u_source_cv(state, cv_nonods, ndim, nphase, den, u_source_cv)

             CALL FORCE_BAL_CTY_ASSEM_SOLVE( state, &
                  NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
                  U_ELE_TYPE, P_ELE_TYPE, &
                  U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
                  U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN,&
                  STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
                  U_SNLOC, P_SNLOC, CV_SNLOC, &
                  X, Y, Z, U_ABS_STAB, U_ABSORB, U_SOURCE, U_SOURCE_CV, &
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
                  IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, &
                  scale_momentum_by_volume_fraction )
          end if

          if (SIG_INT) exit Loop_ITS

          if (solve_saturation) then
             CALL VOLFRA_ASSEM_SOLVE( state, &
                  NCOLACV, FINACV, COLACV, MIDACV, &
                  NCOLCT, FINDCT, COLCT, &
                  CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
                  CV_ELE_TYPE,  &
                  NPHASE,  &
                  CV_NLOC, U_NLOC, X_NLOC,  &
                  CV_NDGLN, X_NDGLN, U_NDGLN, &
                  CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
                  X, Y, Z, U, V, W, &
                  NU, NV, NW, NUOLD, NVOLD, NWOLD, &
                  SATURA, SATURAOLD, &
                  DEN, DENOLD, &
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
                  NITS_FLUX_LIM_VOLFRA, &
                  option_path = '/material_phase[0]/scalar_field::PhaseVolumeFraction', &
                  mass_ele_transp = mass_ele )

          end if

          if (SIG_INT) exit Loop_ITS

          SUM_THETA_FLUX = 1.
          SUM_ONE_M_THETA_FLUX = 0.
          V_SOURCE_COMP = 0.
          IF( NCOMP <= 1 ) THEN
             NCOMP2 = 0
          ELSE 
             NCOMP2 = NCOMP
          END IF

          Loop_COMPONENTS: DO ICOMP = 1, NCOMP2

             !COMP_SOURCE = 0.0

             ! Use values from the previous time step DENOLD,SATURAOLD so its easier to converge
             ! ALPHA_BETA is an order 1 scaling coefficient set up in the input file

             CALL CALC_COMP_ABSORB( ICOMP, NCOMP, DT, ALPHA_BETA, &
                  NPHASE, CV_NONODS, &
                  !.true., K_COMP, &
                  .false., K_COMP, &
                  DENOLD, SATURAOLD, &
                  VOLFRA_PORE, COMP_ABSORB, &
                  TOTELE, CV_NLOC, CV_NDGLN, mass_ele )

             IF( have_option("/material_phase[" // int2str(nstates-ncomp) // &
                  "]/is_multiphase_component/KComp_Sigmoid" ) .and. nphase>1) THEN

                ewrite(3,*) "+++++++++KComp_Sigmoid"

                DO CV_NODI = 1, CV_NONODS
                   !IF( SATURAOLD( CV_NODI ) > 0.8 ) THEN
                   !   COMP_ABSORB( CV_NODI, 1, 2 ) = COMP_ABSORB( CV_NODI, 1, 2 ) * &
                   !        max( 0.01, 5. * ( 1.0 - SATURAOLD( CV_NODI )))
                   !   COMP_ABSORB( CV_NODI, 2, 2 ) = COMP_ABSORB( CV_NODI, 2, 2 ) * &
                   !        max( 0.01, 5. * ( 1.0 - SATURAOLD( CV_NODI )))
                   !ENDIF
                   !IF( SATURAOLD( CV_NODI ) > 0.9 ) THEN
                   !   COMP_ABSORB( CV_NODI, 1, 2 ) = COMP_ABSORB( CV_NODI, 1, 2 ) * &
                   !        max( 0.01, 10. * ( 1.0 - SATURAOLD( CV_NODI )))
                   !   COMP_ABSORB( CV_NODI, 2, 2 ) = COMP_ABSORB( CV_NODI, 2, 2 ) * &
                   !        max( 0.01, 10. * ( 1.0 - SATURAOLD( CV_NODI )))
                   !ENDIF
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
                  TOTELE, MAT_NLOC, CV_NLOC, U_NLOC, X_NLOC, CV_SNLOC, U_SNLOC, &
                  COMP_DIFFUSION, NCOMP_DIFF_COEF, COMP_DIFF_COEF( ICOMP, :, : ), &
                  X_NONODS, X, Y, Z, NU, NV, NW, U_NONODS, CV_NONODS, &
                  MAT_NDGLN, U_NDGLN, X_NDGLN, &
                  U_ELE_TYPE, P_ELE_TYPE )

             ewrite(3,*)'COMP', COMP(( ICOMP - 1 ) * NPHASE * CV_NONODS + 1 : ICOMP * NPHASE * CV_NONODS )
             ewrite(3,*)'COMPold', COMPold(( ICOMP - 1 ) * NPHASE * CV_NONODS + 1 : ICOMP * NPHASE * CV_NONODS )

             ! nits_internal=1 check Copy_Outof_Into_State.F90 for more
             Loop_ITS2: DO ITS2 = 1, 3 !nits_internal

                COMP_GET_THETA_FLUX = .TRUE. ! This will be set up in the input file
                COMP_USE_THETA_FLUX = .FALSE.       

                !ewrite(3,*)'suf_comp_bc:',  SUF_COMP_BC( 1 + STOTEL * CV_SNLOC * NPHASE *( ICOMP - 1 ) : &
                !     &                                                            STOTEL * CV_SNLOC * NPHASE * ICOMP )
                !ewrite(3,*)'SUF_D_BC:', SUF_D_BC
                !ewrite(3,*)'SUF_U_BC:', SUF_U_BC
                !ewrite(3,*)'SUF_V_BC:', SUF_V_BC
                !ewrite(3,*)'SUF_COMP_BC_ROB1:',SUF_COMP_BC_ROB1
                !ewrite(3,*)'SUF_COMP_BC_ROB2:',SUF_COMP_BC_ROB2
                !ewrite(3,*)'COMP_SOURCE:', COMP_SOURCE
                !ewrite(3,*)'COMP_ABSORB:', COMP_ABSORB
                !ewrite(3,*)'COMP_DIFFUSION:', COMP_DIFFUSION

                CALL INTENERGE_ASSEM_SOLVE( state, &
                     NCOLACV, FINACV, COLACV, MIDACV, & ! CV sparsity pattern matrix
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
                     COMP(( ICOMP - 1 ) * NPHASE * CV_NONODS + 1 : ICOMP * NPHASE * CV_NONODS ), &
                     COMPOLD(( ICOMP - 1 ) * NPHASE * CV_NONODS + 1 : ICOMP * NPHASE * CV_NONODS ), &
                     DEN, DENOLD,  &
                     MAT_NLOC, MAT_NDGLN, MAT_NONODS, COMP_DIFFUSION, &
                     V_DISOPT, V_DG_VEL_INT_OPT, DT, V_THETA, V_BETA, &
                     SUF_COMP_BC( 1 + STOTEL * CV_SNLOC * NPHASE * ( ICOMP - 1 ) : STOTEL * CV_SNLOC * NPHASE * ICOMP ), &
                     SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, &
                     SUF_COMP_BC_ROB1, SUF_COMP_BC_ROB2,  &
                     WIC_COMP_BC, WIC_D_BC, WIC_U_BC, &
                     DERIV, P, &
                     COMP_SOURCE, COMP_ABSORB, VOLFRA_PORE, &
                     NDIM,  &
                     NCOLM, FINDM, COLM, MIDM, &
                     XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, LUMP_EQNS, &
                     OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
                     Comp_FEMT(( ICOMP - 1 ) * NPHASE * CV_NONODS + 1 : &
                     ICOMP * NPHASE * CV_NONODS ), Den_FEMT, & 
                     IGOT_T2, SATURA, SATURAOLD, IGOT_THETA_FLUX, SCVNGI_THETA, COMP_GET_THETA_FLUX, COMP_USE_THETA_FLUX, &
                     THETA_FLUX, ONE_M_THETA_FLUX, THETA_GDIFF, &
                     SUF_VOL_BC, SUF_VOL_BC_ROB1, SUF_VOL_BC_ROB2, WIC_VOL_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
                     NOIT_DIM, &
                     NITS_FLUX_LIM_COMP, &
                     MEAN_PORE_CV, &
                     option_path = '', &
                     mass_ele_transp = dummy_ele, &
                     thermal=.false. ) ! the false means that we don't add an extra source term

             END DO Loop_ITS2

             SUM_THETA_FLUX = SUM_THETA_FLUX + THETA_FLUX
             SUM_ONE_M_THETA_FLUX = SUM_ONE_M_THETA_FLUX + ONE_M_THETA_FLUX

!!!!
!!!! AND HERE ZEROED THE COMP_ABSORB TERM WITHIN THE V_SOURCE_COMP
!!!!
             V_SOURCE_COMP = V_SOURCE_COMP + THETA_GDIFF
             DO IPHASE = 1, NPHASE
                DO JPHASE = 1, NPHASE
                   DO CV_NODI = 1, CV_NONODS
                      V_SOURCE_COMP( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) = &
                           V_SOURCE_COMP( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) - &
                           COMP_ABSORB( CV_NODI, IPHASE, JPHASE ) * &
                           COMP( CV_NODI + ( JPHASE - 1 ) * CV_NONODS + &
                           ( ICOMP - 1 ) * NPHASE * CV_NONODS )
                   END DO
                END DO
             END DO

             if (SIG_INT) exit Loop_COMPONENTS

             !ewrite(3,*) 'absorb, source b4-2:', r2norm( comp_absorb, cv_nonods * nphase * &
             !     nphase ), r2norm(V_SOURCE_COMP, cv_nonods * nphase)

          END DO Loop_COMPONENTS

          !ewrite(3,*)'V_SOURCE_COMP b4:', r2norm(V_SOURCE_COMP, cv_nonods * nphase)

          IF( have_option("/material_phase[" // int2str(nstates-ncomp) // &
               "]/is_multiphase_component/Comp_Sum2One") .AND. ( NCOMP2 > 1 ) ) THEN
             ! make sure the composition sums to 1.0 by putting constraint into V_SOURCE_COMP.
             CALL CAL_COMP_SUM2ONE_SOU( V_SOURCE_COMP, CV_NONODS, NPHASE, NCOMP2, DT, ITS, NITS,  &
                  MEAN_PORE_CV, SATURA, SATURAOLD, DEN, DENOLD, COMP, COMPOLD ) 
          ENDIF


          !ewrite(3,*)'V_SOURCE_COMP after:', r2norm(V_SOURCE_COMP, cv_nonods * nphase)

          !ewrite(3,*)'Finished VOLFRA_ASSEM_SOLVE ITS,nits,ITIME:',ITS,nits,ITIME

          !ewrite(3,*)'after all loops -- comp with its=', its
          !ewrite(3,*)'saturation:'
          !do iphase = 1, nphase
          !   ewrite(3,*) iphase, ( satura( ( iphase - 1 ) * cv_nonods + cv_nodi ), &
          !        cv_nodi = 1, cv_nonods )
          !end do

          !ewrite(3,*)''
          !ewrite(3,*)'composition:'
          !do iphase = 1, nphase
          !   do icomp2 = 1, ncomp2
          !      ewrite(3,*) iphase, icomp2, ( comp( ( icomp2 - 1 ) * nphase * cv_nonods + &
          !           ( iphase - 1 ) * cv_nonods + cv_nodi ), cv_nodi = 1, cv_nonods )
          !   end do
          !end do

          !ewrite(3,*)''
          !ewrite(3,*)'composition FEMT:'
          !do iphase = 1, nphase
          !   do icomp2 = 1, ncomp2
          !      ewrite(3,*) iphase, icomp2, ( comp_femt( ( icomp2 - 1 ) * nphase * cv_nonods + &
          !           ( iphase - 1 ) * cv_nonods + cv_nodi ), cv_nodi = 1, cv_nonods )
          !   end do
          !end do

          ! norm_satura2 = r2norm( satura, nphase * cv_nonods )  
          ! ewrite(3,*)'norm_satura b4 new ITS:', norm_satura2

          !if( abs( norm_satura1 - norm_satura2 ) >= 1.e-7 ) stop 9786

          !if (SIG_INT) exit Loop_ITS

       END DO Loop_ITS


       ewrite(3,*)'just after satura subrt:', satura( 1 : nphase * cv_nonods )

       call set_option("/timestepping/current_time", ACCTIM)
       call set_option("/timestepping/timestep", dt)

       Conditional_TIMDUMP: if( ( mod( itime, dump_period_in_timesteps ) == 0 ) .or. ( itime == 1 ) ) then

          if (.false.) then

             ! calculate the quadratic DG representation of the overlapping CVFEM velocity fields
             call overlapping_to_quadratic_dg( &
                  cv_nonods, x_nonods,u_nonods,  totele, &
                  cv_ele_type,   nphase,  cv_nloc, u_nloc, x_nloc, &
                  cv_ndgln,  u_ndgln, x_ndgln, cv_snloc, u_snloc, stotel, cv_sndgln, u_sndgln, &
                  x, y, z,  u, v, w, uold, vold, wold,velocity_dg,ndim,p_ele_type )

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
                     xu, u_nonods, u_nonods * nphase, u, iphase )

                ! close the file for velocity   
                close(output_channel)

                dump_name = 'velocity_phase_proj'//trim(dummy_string_phase)//'.d.'//trim(dummy_string_dump)

                ! open the output file for velocity   
                output_channel = 1
                file_format = 'formatted'
                call open_output_file(output_channel,dump_name,len(dump_name),file_format)

                ! output the fem velocity
                call printing_veloc_field( output_channel, totele, xu_nonods, xu_nloc, xu_ndgln, u_nloc, u_ndgln, &
                     xu, u_nonods, u_nonods * nphase, u, iphase )


                call printing_field_array_veloc(output_channel, totele,   cv_nonods, x_nonods, x_nloc, x_ndgln, cv_nloc, cv_ndgln, &
                     x, velocity_dg(:,iphase,1))
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

          end if

          !! Output vtus from state
          ! Start by copying the interesting files back into state:
          call copy_into_state( state, satura, t, p, u, v, w, den, comp, & 
               ncomp, nphase, cv_ndgln, p_ndgln, u_ndgln, ndim )

          ! find the current time - that reached by the prototype
          call get_option("/timestepping/current_time", current_time)

          ! calc diagnostic fields 
          call calculate_diagnostic_variables(state, exclude_nonrecalculated = .true.)
          call calculate_diagnostic_variables_new(state, exclude_nonrecalculated = .true.)

          ! Call the modern and significantly less satanic version of study
          call write_diagnostics(state, current_time, dt, itime)

          ! Make the vtu dump number dump_no the same as the prototype output number
          dump_no = itime

          ! Write state out to vtu
          call write_state(dump_no, state)

       end if Conditional_TIMDUMP

       if (SIG_INT) exit Loop_Time

    END DO Loop_Time

    ewrite(3,*) 'Leaving solve_multiphase_mom_press_volf'



    deallocate( sigma )
    deallocate( rhs )
    deallocate( rhs_cv )
    deallocate( diag_pres )
    deallocate( femt_print )
    deallocate( femt_dummy )
    deallocate( Sat_FEMT )
    deallocate( T_FEMT )
    deallocate( Den_FEMT )
    deallocate( Comp_FEMT )

  end subroutine solve_multiphase_mom_press_volf


end module multiphase_mom_press_volf
