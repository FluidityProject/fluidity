
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

  module multiphase_time_loop

    use write_state_module
    use diagnostic_variables
    use diagnostic_fields_wrapper
    use diagnostic_fields_new, only : &
         calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
         check_diagnostic_dependencies
    use global_parameters, only: timestep, simulation_start_time, simulation_start_cpu_time, &
                               simulation_start_wall_time, &
                               topology_mesh_name, current_time, is_overlapping, is_compact_overlapping
    use fldebug
    use state_module
    use fields
    use field_options
    use fields_allocates
    use spud
    use signal_vars
    use populate_state_module
    use vector_tools
    use global_parameters

!!$ Modules required by adaptivity
    use qmesh_module
    use adapt_state_module
    use adapt_state_prescribed_module!, only: do_adapt_state_prescribed, adapt_state_prescribed
    use populate_sub_state_module
    use fluids_module!, only: pre_adapt_tasks, update_state_post_adapt
    use parallel_tools

!!$ Modules indigenous to the prototype Code
    use cv_advection, only : cv_count_faces
    use multiphase_1D_engine
    use spact
    use multiphase_EOS
    use multiphase_caching, only: set_caching_level, cache_level
    use shape_functions_Linear_Quadratic
    use Compositional_Terms
    use Copy_Outof_State
    use Copy_BackTo_State
    use checkpoint
    use boundary_conditions

    use multiphase_fractures
    use boundary_conditions_from_options


#ifdef HAVE_ZOLTAN
  use zoltan
#endif
    !use mapping_for_ocvfem
    !use matrix_operations
    !use shape_functions
    !use printout

    implicit none
    private
    public :: MultiFluids_SolveTimeLoop

  contains

    subroutine MultiFluids_SolveTimeLoop( state, &
         dt, nonlinear_iterations, dump_no )
      implicit none
      type( state_type ), dimension( : ), intent( inout ) :: state
      integer, intent( inout ) :: dump_no, nonlinear_iterations
      real, intent( inout ) :: dt

!!$ additional state variables for multiphase & multicomponent

      type(state_type) :: packed_state
      type(state_type), dimension(:), pointer :: multiphase_state, multicomponent_state


!!$ Primary scalars
      integer :: nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods
      real :: dx

!!$ Node global numbers
      integer, dimension( : ), pointer :: x_ndgln_p1, x_ndgln, cv_ndgln, p_ndgln, &
           mat_ndgln, u_ndgln, xu_ndgln, cv_sndgln, p_sndgln, u_sndgln

!!$ Sparsity patterns
      integer :: nlenmcy, mx_nface_p1, mx_ncolacv, mxnele, mx_ncoldgm_pha, &
           mx_ncolmcy, mx_nct, mx_nc, mx_ncolcmc, mx_ncolm, &
           ncolacv, ncolmcy, ncolele, ncoldgm_pha, ncolct, ncolc, ncolcmc, ncolm
      integer, dimension( : ), allocatable :: finacv, midacv, finmcy,  midmcy, &
           finele, midele, findgm_pha, middgm_pha, findct, &
           findc, findcmc, midcmc, findm, &
           midm
      integer, dimension(:), pointer :: colacv, colmcy, colele, colct,colm,colc,colcmc,coldgm_pha
      integer, dimension(:), pointer :: small_finacv, small_colacv, small_midacv
      integer, dimension(:), pointer :: block_to_global_acv
      integer, dimension(:,:), allocatable :: global_dense_block_acv

!!$ Defining element-pair type and discretisation options and coefficients
      integer :: cv_ele_type, p_ele_type, u_ele_type, mat_ele_type, u_sele_type, cv_sele_type, &
           t_disopt, v_disopt, t_dg_vel_int_opt, u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
           comp_diffusion_opt, ncomp_diff_coef, in_ele_upwind, dg_ele_upwind, nopt_vel_upwind_coefs, &
           nits_flux_lim_t, nits_flux_lim_volfra, nits_flux_lim_comp,  IDIVID_BY_VOL_FRAC
      logical :: volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, comp_get_theta_flux, &
           t_use_theta_flux, t_get_theta_flux, scale_momentum_by_volume_fraction, q_scheme
      real :: t_beta, v_beta, t_theta, v_theta, u_theta
      real, dimension( : ), allocatable :: opt_vel_upwind_coefs

!!$ Defining time- and nonlinear interations-loops variables
      integer :: itime, dump_period_in_timesteps, final_timestep, &
           NonLinearIteration, NonLinearIteration_Components, dtime
      real :: acctim, finish_time

!!$ Defining problem that will be solved
      logical :: have_temperature_field, have_component_field, have_extra_DiffusionLikeTerm, &
           solve_force_balance, solve_PhaseVolumeFraction

!!$ Defining solver options
      integer :: velocity_max_iterations, PhaseVolumeFraction_max_iterations

!!$ Shape function related fields:
      integer :: cv_ngi, cv_ngi_short, scvngi_theta, sbcvngi, nface, igot_t2, igot_theta_flux, IGOT_THERM_VIS

!!$ For output:
      real, dimension( : ), allocatable :: PhaseVolumeFraction_FEMT, Temperature_FEMT, Density_FEMT, &
           Component_FEMT, Mean_Pore_CV, SumConc_FEMT, Dummy_PhaseVolumeFraction_FEMT
      type( scalar_field ), pointer :: Component_State

!!$ Variables that can be effectively deleted as they are not used anymore:
      integer :: noit_dim

!!$ Variables used in the diffusion-like term: capilarity and surface tension:
      integer :: iplike_grad_sou
      real, dimension( : ), allocatable :: plike_grad_sou_grad, plike_grad_sou_coef

!!$ Adaptivity related fields and options:
      type( tensor_field ) :: metric_tensor
      type( state_type ), dimension( : ), pointer :: sub_state => null()
      integer :: nonlinear_iterations_adapt
      logical :: do_reallocate_fields, not_to_move_det_yet = .false., initialised

!!$ Working arrays:
      real, dimension( : ), pointer :: &
           Temperature, PhaseVolumeFraction, &
           Component, Temperature_Old, &
           PhaseVolumeFraction_Old, Component_Old, &
           Velocity_U_Source, Velocity_U_Source_CV, Temperature_Source, PhaseVolumeFraction_Source, &
           ScalarField_Source, Component_Source, ScalarAdvectionField_Source, &
           ScalarField_Source_Store, ScalarField_Source_Component, &
           mass_ele, dummy_ele

      real, dimension( :, :, :, : ), allocatable :: THERM_U_DIFFUSION
      real, dimension( :, : ), allocatable :: THERM_U_DIFFUSION_VOL

      real, dimension( :, : ), pointer :: THETA_GDIFF

      real, dimension( :, : ), pointer ::  DRhoDPressure, FEM_VOL_FRAC
!!$
      real, dimension( :, :, : ), allocatable :: Permeability, Material_Absorption, Material_Absorption_Stab, &
           Velocity_Absorption, ScalarField_Absorption, Component_Absorption, Temperature_Absorption, &
!!$
           Component_Diffusion_Operator_Coefficient
      real, dimension( :, :, :, : ), allocatable :: Momentum_Diffusion, ScalarAdvectionField_Diffusion, &
           Component_Diffusion
      real, dimension( :, : ), allocatable :: Momentum_Diffusion_Vol
           
      real, dimension( :, : ), allocatable ::theta_flux, one_m_theta_flux, theta_flux_j, one_m_theta_flux_j, &
                               sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j

!!$ Material_Absorption_Stab = u_abs_stab; Material_Absorption = u_absorb; ScalarField_Absorption = v_absorb
!!$ Component_Absorption = comp_absorb; ScalarAdvectionField_Absorption = t_absorb
!!$ Velocity_U_Source = u_source, ScalarField_Source = v_source, Component_Source = comp_source,
!!$ ScalarAdvectionField_Source = Temperature_Source = t_source; Component_Diffusion_Operator_Coefficient = comp_diff_coef
!!$ Momentum_Diffusion = udiffusion; ScalarAdvectionField_Diffusion = tdiffusion, 
!!$ Component_Diffusion = comp_diffusion

      !character( len = option_path_len ) :: eos_option_path( 1 )

      integer :: stat, istate, iphase, jphase, icomp, its, its2, cv_nodi, adapt_time_steps, cv_inod
      real, dimension( : ), allocatable :: rsum

      real, dimension(:, :), allocatable :: SUF_SIG_DIAGTEN_BC

      type( scalar_field ), pointer :: cfl, rc_field
      real :: c, rc, minc, maxc, ic
      !Variables for adaptive time stepping based on non-linear iterations
      logical :: nonLinearAdaptTs, Repeat_time_step, ExitNonLinearLoop
      real, dimension(:,:,:), allocatable  :: reference_field

      type( tensor_field ), pointer :: NU_s, NUOLD_s, U_s, UOLD_s, D_s, DOLD_s, DC_s, DCOLD_s
      type( tensor_field ), pointer :: MFC_s, MFCOLD_s!, MFC_FEMT_s, MFCOLD_FEMT_s

      !! face value storage
      integer :: ncv_faces
      real::  second_theta

      integer :: checkpoint_number

      !Variable to store where we store things. Do not oversize this array, the size has to be the last index in use
      integer, dimension (31) :: StorageIndexes
      !Distribution of the indexes of StorageIndexes:
      !cv_fem_shape_funs_plus_storage: 1 (ASSEMB_FORCE_CTY), 13 (CV_ASSEMB)
      !CALC_ANISOTROP_LIM            : 2 (DETNLXR_PLUS_U_WITH_STORAGE in the inside, maybe 14 as well?)
      !DG_DERIVS_ALL2                : 3 (DETNLXR_PLUS_U_WITH_STORAGE in the inside, maybe 14 as well?)
      !DETNLXR_INVJAC                : 4
      !UNPACK_LOC                    : 5,6,7,8,9,10 (disabled)
      !COLOR_GET_CMC_PHA             : 11 (can be optimised, now it is not using only pointers)
      !Matrix C                      : 12
      !DG_DERIVS_ALL                 : 14 (DETNLXR_PLUS_U_WITH_STORAGE in the inside)
      !DETNLXR_PLUS_U_WITH_STORAGE   : 14
      !Indexes used in SURFACE_TENSION_WRAPPER (deprecated and will be removed):[15,30]
      !PROJ_CV_TO_FEM_state          : 31 (disabled)

      !Working pointers
      real, dimension(:,:), pointer :: SAT_s, OldSAT_s, FESAT_s

      type( tensor_field ), pointer :: tracer_field, velocity_field, density_field, saturation_field, old_saturation_field
      type(scalar_field), pointer :: pressure_field, porosity_field, tracer_field2

      !Dummy to print FEM saturation
      real, dimension(:,:), allocatable :: dummy_to_print_FEM

      logical :: write_all_stats=.true.

#ifdef HAVE_ZOLTAN
      real(zoltan_float) :: ver
      integer(zoltan_int) :: ierr
      
      ierr = Zoltan_Initialize(ver)  
      assert(ierr == ZOLTAN_OK)
#endif

      !Initially we set to use Stored data and that we have a new mesh
      StorageIndexes = 0!Initialize them as zero !


      !Read info for adaptive timestep based on non_linear_iterations

    !! JRP changes to make a multiphasic state
      call pack_multistate(state,packed_state,multiphase_state,&
           multicomponent_state)

      call set_caching_level()

    !Get from packed_state
    call get_var_from_packed_state(packed_state,PhaseVolumeFraction = SAT_s,&
    OldPhaseVolumeFraction=OldSAT_s,FEPhaseVolumeFraction = FESAT_s )
      IDIVID_BY_VOL_FRAC=0
      !call print_state( packed_state )
      !stop 78

      !  Access boundary conditions via a call like
      !  call get_entire_boundary_condition(extract_tensor_field(packed_state,"Packed"//name),["weakdirichlet"],tfield,bc_type_list)
      !  where tfield is type(tensor_field) and bc_type_list is integer, dimension(tfield%dim(1),tfield%dim(2),nonods)
      !  Then values are in tfield%val(1/ndim/ncomp,nphase,nonods) 
      !  Type ids are in bc_type_list(1/ndim/ncomp,nphase,stotel) 
      !
      !A deallocate tfield when finished!!

      Repeat_time_step = .false.!Initially has to be false
      nonLinearAdaptTs = have_option(  '/timestepping/nonlinear_iterations/nonlinear_iterations_automatic/adaptive_timestep_nonlinear')

!     !If adaptive time_stepping then we need to create backup_state
!    if (nonLinearAdaptTs)  call pack_multistate(state,backup_state,multiphase_state,&
!           multicomponent_state)

!!$ Compute primary scalars used in most of the code
      call Get_Primary_Scalars( state, &         
           nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, dx )

!!$ Calculating Global Node Numbers
      allocate( cv_sndgln( stotel * cv_snloc ), p_sndgln( stotel * p_snloc ), &
           u_sndgln( stotel * u_snloc ) )


!      x_ndgln_p1 = 0 ; x_ndgln = 0 ; cv_ndgln = 0 ; p_ndgln = 0 ; mat_ndgln = 0 ; u_ndgln = 0 ; xu_ndgln = 0 ; &
           cv_sndgln = 0 ; p_sndgln = 0 ; u_sndgln = 0

      call Compute_Node_Global_Numbers( state, &
           totele, stotel, x_nloc, x_nloc_p1, cv_nloc, p_nloc, u_nloc, xu_nloc, &
           cv_snloc, p_snloc, u_snloc, &
           cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, xu_ndgln, mat_ndgln, &
           cv_sndgln, p_sndgln, u_sndgln )

!!$
!!$ Computing Sparsity Patterns Matrices
!!$ 
!!$ Defining lengths and allocating space for the matrices
      call Defining_MaxLengths_for_Sparsity_Matrices( ndim, nphase, totele, u_nloc, cv_nloc, cv_nonods, &
           mx_nface_p1, mxnele, mx_nct, mx_nc, mx_ncolcmc, mx_ncoldgm_pha, mx_ncolmcy, &
           mx_ncolacv, mx_ncolm )
      nlenmcy = u_nonods * nphase * ndim + cv_nonods
      allocate( finacv( cv_nonods * nphase + 1 ), colacv( mx_ncolacv ), midacv( cv_nonods * nphase ), &
           finmcy( nlenmcy + 1 ), colmcy( mx_ncolmcy ), midmcy( nlenmcy ), &
           finele( totele + 1 ), colele( mxnele ), midele( totele ), &
           findgm_pha( u_nonods * nphase * ndim + 1 ), coldgm_pha( mx_ncoldgm_pha ), &
           middgm_pha( u_nonods * nphase * ndim ), &
           findct( cv_nonods + 1 ), colct( mx_nct ), &
           findc( u_nonods + 1 ), colc( mx_nc ), &
           findcmc( cv_nonods + 1 ), colcmc( 0 ), midcmc( cv_nonods ), &
           findm( cv_nonods + 1 ), colm( mx_ncolm ), midm( cv_nonods ) )


      allocate( global_dense_block_acv( nphase , cv_nonods ) )

      finacv = 0 ; colacv = 0 ; midacv = 0 ; finmcy = 0 ; colmcy = 0 ; midmcy = 0 ; finele = 0
      colele = 0 ; midele = 0 ; findgm_pha = 0 ; coldgm_pha = 0 ; middgm_pha = 0 ; findct = 0
      colct = 0 ; findc = 0 ; colc = 0 ; findcmc = 0 ; colcmc = 0 ; midcmc = 0 ; findm = 0
      colm = 0 ; midm = 0

!!$ Defining element-pair type 
      call Get_Ele_Type( x_nloc, cv_ele_type, p_ele_type, u_ele_type, &
           mat_ele_type, u_sele_type, cv_sele_type )

!!$ Sparsity Patterns Matrices 
      call Get_Sparsity_Patterns( state, &
!!$ CV multi-phase eqns (e.g. vol frac, temp)
           mx_ncolacv, ncolacv, finacv, colacv, midacv, &
           small_finacv, small_colacv, small_midacv, &
           block_to_global_acv, global_dense_block_acv, &
!!$ Force balance plus cty multi-phase eqns
           nlenmcy, mx_ncolmcy, ncolmcy, finmcy, colmcy, midmcy, &
!!$ Element connectivity
           mxnele, ncolele, midele, finele, colele, &
!!$ Force balance sparsity
           mx_ncoldgm_pha, ncoldgm_pha, coldgm_pha, findgm_pha, middgm_pha, &
!!$ CT sparsity - global continuity eqn
           mx_nct, ncolct, findct, colct, &
!!$ C sparsity operating on pressure in force balance
           mx_nc, ncolc, findc, colc, &
!!$ pressure matrix for projection method
           mx_ncolcmc, ncolcmc, findcmc, colcmc, midcmc, &
!!$ CV-FEM matrix
           mx_ncolm, ncolm, findm, colm, midm, mx_nface_p1 )

    call temp_mem_hacks()

      Q_SCHEME = have_option( '/material_phase[0]/scalar_field::Temperature/prognostic/spatial_discretisation/control_volumes/q_scheme' )
      IGOT_THERM_VIS = 0
      IF ( Q_SCHEME ) IGOT_THERM_VIS = 1

!!$ Allocating space for various arrays:
      allocate( &
!!$
           Temperature( nphase * cv_nonods ), &
           PhaseVolumeFraction( nphase * cv_nonods ), Component( nphase * cv_nonods * ncomp ), &
           DRhoDPressure( nphase, cv_nonods ), FEM_VOL_FRAC( nphase, cv_nonods ),&
!!$
           Temperature_Old( nphase * cv_nonods ), &
           PhaseVolumeFraction_Old( nphase * cv_nonods ), Component_Old( nphase * cv_nonods * ncomp ), &
!!$
           suf_sig_diagten_bc( stotel * cv_snloc * nphase, ndim ), &
           PhaseVolumeFraction_FEMT( cv_nonods * nphase ), Temperature_FEMT( cv_nonods * nphase ), &
           Density_FEMT( cv_nonods * nphase ), Component_FEMT( cv_nonods * nphase * ncomp ), &
           Mean_Pore_CV( cv_nonods ),  SumConc_FEMT( cv_nonods * ncomp ), &
           Dummy_PhaseVolumeFraction_FEMT( cv_nonods * nphase ), dummy_ele( totele ), mass_ele( totele ), &
!!$
           Temperature_Source( nphase * cv_nonods ), &
           PhaseVolumeFraction_Source( cv_nonods * nphase ), Velocity_U_Source( u_nonods * nphase * ndim ), &
           Velocity_U_Source_CV( cv_nonods * nphase * ndim ), Component_Source( cv_nonods * nphase ), &
           ScalarField_Source( cv_nonods * nphase ), ScalarAdvectionField_Source( cv_nonods * nphase ), &
!!$
           Permeability( totele, ndim, ndim ), &
!!$
           Material_Absorption( mat_nonods, ndim * nphase, ndim * nphase ), &
           Velocity_Absorption( mat_nonods, ndim * nphase, ndim * nphase ), &
           Material_Absorption_Stab( mat_nonods, ndim * nphase, ndim * nphase ), & 
           ScalarField_Absorption( cv_nonods, nphase, nphase ), Component_Absorption( cv_nonods, nphase, nphase ), &
           Temperature_Absorption( cv_nonods, nphase, nphase ), &
           Momentum_Diffusion( mat_nonods, ndim, ndim, nphase ), &
           Momentum_Diffusion_Vol( mat_nonods, nphase ), &
           ScalarAdvectionField_Diffusion( mat_nonods, ndim, ndim, nphase ), &
           Component_Diffusion( mat_nonods, ndim, ndim, nphase ), &
!!$ Variables used in the diffusion-like term: capilarity and surface tension:
           plike_grad_sou_grad( cv_nonods * nphase ), &
           plike_grad_sou_coef( cv_nonods * nphase ), &
           THERM_U_DIFFUSION(NDIM,NDIM,NPHASE,MAT_NONODS*IGOT_THERM_VIS ), THERM_U_DIFFUSION_VOL(NPHASE,MAT_NONODS*IGOT_THERM_VIS ) )
 
      ncv_faces=CV_count_faces( packed_state, CV_ELE_TYPE, stotel, cv_sndgln, u_sndgln)

!!$
      Temperature=0.
      PhaseVolumeFraction=0. ; Component=0.
      DRhoDPressure=0.
!!$
      Temperature_Old=0.
      PhaseVolumeFraction_Old=0. ; Component_Old=0.
!!$
      Temperature_Source=0.
      suf_sig_diagten_bc=0.
!!$
      PhaseVolumeFraction_FEMT=0. ; Temperature_FEMT=0.
      Density_FEMT=1. ; Component_FEMT=0.
      Mean_Pore_CV=0. ; SumConc_FEMT=0.
      Dummy_PhaseVolumeFraction_FEMT=0. ; dummy_ele=0. ; mass_ele=0.
!!$
      PhaseVolumeFraction_Source=0. ; Velocity_U_Source=0.
      Velocity_U_Source_CV=0. ; Component_Source=0.
      ScalarField_Source=0. ; ScalarAdvectionField_Source=0.
!!$
      Permeability=0.
!!$
      Material_Absorption=0.
      Velocity_Absorption=0.
      Material_Absorption_Stab=0.
      ScalarField_Absorption=0. ; Component_Absorption=0.
      Temperature_Absorption=0.
      Momentum_Diffusion=0.
      Momentum_Diffusion_Vol=0.
      ScalarAdvectionField_Diffusion=0.
      Component_Diffusion=0.
      THERM_U_DIFFUSION=0.
      THERM_U_DIFFUSION_VOL=0.
!!$
      plike_grad_sou_grad=0.
      plike_grad_sou_coef=0.
      iplike_grad_sou=0 


      tracer_field=>extract_tensor_field(packed_state,"PackedTemperature",stat)
      if(stat==0)then
         do iphase = 1, nphase
            tracer_field2=>extract_scalar_field(state(iphase),"DummyT",stat)
            if(stat==0)then
               tracer_field2%val = tracer_field%val(1,iphase,:)
            end if
         end do
      end if



!!$ Extracting Mesh Dependent Fields
      initialised = .false.
      call Extracting_MeshDependentFields_From_State( state, packed_state, initialised, &
           SAT_s, PhaseVolumeFraction_Source,&
           Component, Component_Source, &
           Velocity_U_Source, Velocity_Absorption, &
           Temperature, Temperature_Source, &
           Permeability )
      FESAT_s = 0; OldSAT_s = 0.
!!$ Calculate diagnostic fields
      call calculate_diagnostic_variables( state, exclude_nonrecalculated = .true. )
      call calculate_diagnostic_variables_new( state, exclude_nonrecalculated = .true. )

!!$ Dummy field used in the scalar advection option:
      Dummy_PhaseVolumeFraction_FEMT = 1.

!!$
!!$ Initialising Robin boundary conditions --  this still need to be defined in the schema:
!!$

!!$
!!$ Initialising Absorption terms that do not appear in the schema
!!$
      ScalarField_Absorption = 0. ; Component_Absorption = 0. ; Temperature_Absorption = 0.

!!$ Variables that can be effectively deleted as they are not used anymore:
      noit_dim = 0

!!$ Computing shape function scalars
      igot_t2 = 0 ; igot_theta_flux = 0
      if( ncomp /= 0 )then
         igot_t2 = 1 ; igot_theta_flux = 1
      end if

      call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
           cv_ngi, cv_ngi_short, scvngi_theta, sbcvngi, nface, .false. )

      allocate( theta_flux( nphase, ncv_faces * igot_theta_flux ), &
           one_m_theta_flux( nphase, ncv_faces * igot_theta_flux ), &
           theta_flux_j( nphase, ncv_faces * igot_theta_flux ), &
           one_m_theta_flux_j( nphase, ncv_faces * igot_theta_flux ), &
           sum_theta_flux( nphase, ncv_faces * igot_theta_flux ), &
           sum_one_m_theta_flux( nphase, ncv_faces * igot_theta_flux ), &
           sum_theta_flux_j( nphase, ncv_faces * igot_theta_flux ), &
           sum_one_m_theta_flux_j( nphase, ncv_faces * igot_theta_flux ), &
           theta_gdiff( nphase, cv_nonods ), ScalarField_Source_Store( cv_nonods * nphase ), &
           ScalarField_Source_Component( cv_nonods * nphase ) )

      sum_theta_flux = 1. ; sum_one_m_theta_flux = 0.
      sum_theta_flux_j = 1. ; sum_one_m_theta_flux_j = 0.
      ScalarField_Source_Store=0. ; ScalarField_Source_Component=0.

!!$ Defining discretisation options
      call Get_Discretisation_Options( state, &
           t_disopt, v_disopt, t_beta, v_beta, t_theta, v_theta, u_theta, &
           t_dg_vel_int_opt, u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
           comp_diffusion_opt, ncomp_diff_coef, in_ele_upwind, dg_ele_upwind, &
           nits_flux_lim_t, nits_flux_lim_volfra, nits_flux_lim_comp, &
           volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, comp_get_theta_flux, &
           t_use_theta_flux, t_get_theta_flux, scale_momentum_by_volume_fraction )

      allocate( Component_Diffusion_Operator_Coefficient( ncomp, ncomp_diff_coef, nphase ) )
      Component_Diffusion_Operator_Coefficient = 0.

!!$ Option not currently set up in the schema and zeroed from the begining. It is used to control 
!!$ the upwinding rate (in the absorption term) during advection/assembling.
      nopt_vel_upwind_coefs = mat_nonods * nphase * ndim * ndim * 2
      allocate( opt_vel_upwind_coefs( nopt_vel_upwind_coefs ) ) ; opt_vel_upwind_coefs = 0.

!!$ Defining problem to be solved:
      call get_option( '/material_phase[0]/vector_field::Velocity/prognostic/solver/max_iterations', &
           velocity_max_iterations,  default =  500 )
      call get_option( '/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/solver/max_iterations', &
           PhaseVolumeFraction_max_iterations,  default =  500 )

      solve_force_balance = .false. ; solve_PhaseVolumeFraction = .false.
      if( velocity_max_iterations /= 0 ) solve_force_balance = .true.
      if( PhaseVolumeFraction_max_iterations /= 0 ) solve_PhaseVolumeFraction = .true. 

!!$ Setting up variables for the Time- and NonLinear Iterations-Loops:
      call get_option( '/timestepping/current_time', acctim )
      call get_option( '/timestepping/timestep', dt )
      call get_option( '/timestepping/finish_time', finish_time )
      call get_option( '/io/dump_period_in_timesteps/constant', dump_period_in_timesteps, default = 1 )
      call get_option( '/timestepping/nonlinear_iterations', NonLinearIteration, default = 3 )
!      call get_option( '/timestepping/nonlinear_iterations/nonlinear_iterations_automatic', tolerance_between_non_linear, default = -1. )
!!$
      have_temperature_field = .false. ; have_component_field = .false. ; have_extra_DiffusionLikeTerm = .false.
      do istate = 1, nstate
         if( have_option( '/material_phase[' // int2str( istate - 1 ) // ']/scalar_field::Temperature' ) ) &
              have_temperature_field = .true.
         if( have_option( '/material_phase[' // int2str( istate - 1 ) // ']/is_multiphase_component' ) ) &
              have_component_field = .true.
!!$
 if( have_temperature_field ) then
            call Calculate_All_Rhos( state, packed_state, ncomp, nphase, ndim, cv_nonods, cv_nloc, totele, &
                 cv_ndgln, DRhoDPressure )
      end if

         if( have_component_field ) then
            call get_option( '/material_phase[' // int2str( istate - 1 ) // 'scalar_field::' // &
                 'ComponentMassFractionPhase1/prognostic/temporal_discretisation/control_volumes' // &
                 '/number_advection_iterations', NonLinearIteration_Components, default = 3 )
         end if
      end do

      if( have_option( '/material_phase[0]/multiphase_properties/capillary_pressure' ) ) &
           have_extra_DiffusionLikeTerm = .true.

      if ( have_option( '/mesh_adaptivity/hr_adaptivity' ) ) then
         call allocate( metric_tensor, extract_mesh(state(1), topology_mesh_name), 'ErrorMetric' )
      end if


   !   print *,'u_nonods, cv_nonods, totele:',u_nonods, cv_nonods, totele
   !   print *,'mx_ncolacv, ncolacv:',mx_ncolacv, ncolacv
   !   print *,'nlenmcy, mx_ncolmcy, ncolmcy,mxnele, ncolele:',nlenmcy, mx_ncolmcy, ncolmcy,mxnele, ncolele
   !   print *,'mx_ncoldgm_pha, ncoldgm_pha:',mx_ncoldgm_pha, ncoldgm_pha
   !   print *,'mx_nct, ncolct,mx_nc, ncolc, mx_ncolcmc, ncolcmc:',mx_nct, ncolct,mx_nc, ncolc, mx_ncolcmc, ncolcmc
   !   print *,'mx_ncolm, ncolm:',mx_ncolm, ncolm
   !   stop 282

    if (have_component_field) then
        !######TEMPORARY CONVERSION FROM OLD PhaseVolumeFraction TO PACKED######
        do cv_inod = 1, size(SAT_s,2)
            do iphase = 1, size(SAT_s,1)
                phaseVolumeFraction(cv_inod +(iphase-1)*size(SAT_s,2)) = SAT_s(iphase,cv_inod)
                PhaseVolumeFraction_Old(cv_inod +(iphase-1)*size(SAT_s,2)) = OldSAT_s(iphase,cv_inod)
            end do
        end do
        !#############################################################
    end if
!!$ Starting Time Loop
      itime = 0
      dtime = 0
      checkpoint_number=1
      Loop_Time: do
!!$

!print *, '    NEW DT', itime+1

         itime = itime + 1
         timestep = itime
         call get_option( '/timestepping/timestep', dt )


         acctim = acctim + dt
         call set_option( '/timestepping/current_time', acctim )
         new_lim = .true.

         if ( acctim > finish_time ) then 
            ewrite(1,*) "Passed final time"
            exit Loop_Time
         end if

         call get_option( '/timestepping/final_timestep', final_timestep, stat )
         if( stat == spud_no_error ) then
            if( itime > final_timestep ) then
               ewrite(1,*) "Passed final timestep"
               exit Loop_Time
            end if
         end if

        ExitNonLinearLoop = .false.
        !Store backup to be able to repeat a timestep
         if (nonLinearAdaptTs) call Adaptive_NonLinear(packed_state, reference_field, its, &
        Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs,1)


!!$ Update all fields from time-step 'N - 1'
         if (have_component_field) PhaseVolumeFraction_Old = PhaseVolumeFraction
         Temperature_Old = Temperature ; Component_Old = Component

         U_s  => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedVelocity" )
         UOLD_s  => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedOldVelocity" )

         NU_s => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedNonlinearVelocity" )
         NUOLD_s => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedOldNonlinearVelocity" )

         D_s  => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedDensity" )
         DOLD_s  => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedOldDensity" )
         if( have_component_field ) then
            DC_s  => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedComponentDensity" )
            DCOLD_s  => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedOldComponentDensity" )
            MFC_s  => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedComponentMassFraction" )
            MFCOLD_s  => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedOldComponentMassFraction" )
         end if

         porosity_field=>extract_scalar_field(packed_state,"Porosity")

!! Temporal working array:
         DO CV_INOD = 1, CV_NONODS
            DO IPHASE = 1, NPHASE
               DO ICOMP = 1, NCOMP
                  MFC_s%val(ICOMP, IPHASE, CV_INOD) = &
                      COMPONENT( CV_INOD + ( IPHASE - 1 ) * CV_NONODS + ( ICOMP - 1 ) * NPHASE * CV_NONODS )
                  MFCOLD_s%val(ICOMP, IPHASE, CV_INOD) = &         
                      COMPONENT_OLD( CV_INOD + ( IPHASE - 1 ) * CV_NONODS + ( ICOMP - 1 ) * NPHASE * CV_NONODS )
               END DO
            END DO
         END DO

         ! evaluate prescribed fields at time = current_time+dt
         call set_prescribed_field_values( state, exclude_interpolated = .true., &
              exclude_nonreprescribed = .true., time = acctim )
        !! Update all fields from time-step 'N - 1'
         call copy_packed_new_to_old( packed_state )

         ! update velocity absorption
         call update_velocity_absorption( state, ndim, nphase, mat_nonods, velocity_absorption )

!!$ FEMDEM...
#ifdef USING_FEMDEM
         if ( have_option( '/blasting' ) ) then
            call blasting( packed_state, nphase )
            call update_blasting_memory( packed_state, state, timestep )
         end if
#endif

!!$ Start non-linear loop
         Loop_NonLinearIteration: do  its = 1, NonLinearIteration

!print *, '  NEW ITS', its

        !To force the recalculation of all the stored variables uncomment the following line:
!        call Clean_Storage(state, StorageIndexes)


         if( have_temperature_field .and. &
              have_option( '/material_phase[0]/scalar_field::Temperature/prognostic' ) ) then
            call get_option( '/material_phase[0]/scalar_field::Temperature/prognostic/temporal_discretisation' // &
              '/control_volumes/second_theta', second_theta, default=1. )
         end if

         if( have_component_field ) then
            call get_option( '/material_phase[' // int2str( nphase ) // ']/scalar_field::ComponentMassFractionPhase1/' // &
              'prognostic/temporal_discretisation/control_volumes/second_theta', second_theta, default=1. )
         end if

            !Store the field we want to compare with to check how are the computations going
            call Adaptive_NonLinear(packed_state, reference_field, its, &
            Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs,2)

            call Calculate_All_Rhos( state, packed_state, ncomp, nphase, ndim, cv_nonods, cv_nloc, totele, &
                 cv_ndgln, DRhoDPressure )


            if( its == 1 ) then
               DOLD_s%val = D_s%val
               if( have_component_field ) then
                  DCOLD_s%val = DC_s%val
                  MFCOLD_s%val = MFC_s%val
               end if

            end if

            if( solve_force_balance ) then
               call Calculate_AbsorptionTerm( state, packed_state,&
                    cv_ndgln, mat_ndgln, &
                    nopt_vel_upwind_coefs, opt_vel_upwind_coefs, Material_Absorption )

               ! calculate SUF_SIG_DIAGTEN_BC this is \sigma_in^{-1} \sigma_out
               ! \sigma_in and \sigma_out have the same anisotropy so SUF_SIG_DIAGTEN_BC
               ! is diagonal
               if( is_overlapping .or. is_compact_overlapping ) then
                  call calculate_SUF_SIG_DIAGTEN_BC( packed_state, suf_sig_diagten_bc, totele, stotel, cv_nloc, &
                       cv_snloc, nphase, ndim, nface, mat_nonods, cv_nonods, x_nloc, ncolele, cv_ele_type, &
                       finele, colele, cv_ndgln, cv_sndgln, x_ndgln, mat_ndgln, material_absorption, &
                       state, x_nonods )
               end if
            end if
!!$ Solve advection of the scalar 'Temperature':
            Conditional_ScalarAdvectionField: if( have_temperature_field .and. &
                 have_option( '/material_phase[0]/scalar_field::Temperature/prognostic' ) ) then
               ewrite(3,*)'Now advecting Temperature Field'

               NU_s % val = U_s % val
               NUOLD_s % val = UOLD_s % val

               call calculate_diffusivity( state, ncomp, nphase, ndim, cv_nonods, mat_nonods, &
                                           mat_nloc, totele, mat_ndgln, ScalarAdvectionField_Diffusion )


               tracer_field=>extract_tensor_field(packed_state,"PackedTemperature")
               velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
               density_field=>extract_tensor_field(packed_state,"PackedDensity",stat)

               call INTENERGE_ASSEM_SOLVE( state, packed_state, &
                    tracer_field,velocity_field,density_field,&
                    NCOLACV, FINACV, COLACV, MIDACV, &
                    small_FINACV, small_COLACV, small_MIDACV, &
                    block_to_global_acv, global_dense_block_acv, &
                    NCOLCT, FINDCT, COLCT, &
                    CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
                    U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE,  &
                    NPHASE, &
                    CV_NLOC, U_NLOC, X_NLOC, &
                    CV_NDGLN, X_NDGLN, U_NDGLN, &
                    CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
!!$
                    Temperature, Temperature_Old, &
!!$
                    MAT_NLOC, MAT_NDGLN, MAT_NONODS, ScalarAdvectionField_Diffusion, IGOT_THERM_VIS, THERM_U_DIFFUSION, THERM_U_DIFFUSION_VOL, &
                    t_disopt, t_dg_vel_int_opt, dt, t_theta, t_beta, &
                    suf_sig_diagten_bc,&
                    DRhoDPressure, &
                    Temperature_Source, Temperature_Absorption, Porosity_field%val, &
                    ndim, &
!!$
                    NCOLM, FINDM, COLM, MIDM, &
!!$
                    XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
!!$
                    opt_vel_upwind_coefs, nopt_vel_upwind_coefs, &
                    Temperature_FEMT, Dummy_PhaseVolumeFraction_FEMT, &
                    0,Temperature, Temperature_Old,igot_theta_flux, scvngi_theta, &
                    t_get_theta_flux, t_use_theta_flux, &
                    THETA_GDIFF, &
                    in_ele_upwind, dg_ele_upwind, &
!!$                    
                    NOIT_DIM, & ! This need to be removed as it is already deprecated
!!$
                    Mean_Pore_CV, &
                    option_path = '/material_phase[0]/scalar_field::Temperature', &
                    mass_ele_transp = dummy_ele, &
                    thermal = have_option( '/material_phase[0]/scalar_field::Temperature/prognostic/equation::InternalEnergy'),&
                    StorageIndexes=StorageIndexes )

!!$  Update state memory
!!$               do iphase = 1, nphase
!!$                  Temperature_State => extract_scalar_field( state( iphase ), 'Temperature' )
!!$                  Temperature_State % val = Temperature( 1 + ( iphase - 1 ) * cv_nonods : iphase * cv_nonods )
!!$               end do

               do iphase = 1, nphase
                  tracer_field2=>extract_scalar_field(state(iphase),"DummyT")
                  tracer_field2%val = tracer_field%val(1,iphase,:)
               end do

               call Calculate_All_Rhos( state, packed_state, ncomp, nphase, ndim, cv_nonods, cv_nloc, totele, &
                    cv_ndgln, DRhoDPressure )

            end if Conditional_ScalarAdvectionField

            ScalarField_Source_Store = ScalarField_Source + ScalarField_Source_Component

            volfra_use_theta_flux = .true.
            if( ncomp <= 1 ) volfra_use_theta_flux = .false.

!!$ Now solving the Momentum Equation ( = Force Balance Equation )
            Conditional_ForceBalanceEquation: if ( solve_force_balance ) then

!!$ Updating velocities:
               NU_s % val = U_s % val
               NUOLD_s % val = UOLD_s % val

!!$ Diffusion-like term -- here used as part of the capillary pressure for porous media. It can also be
!!$ extended to surface tension -like term.
               iplike_grad_sou = 0
               plike_grad_sou_grad = 0


               CALL CALCULATE_SURFACE_TENSION( state, packed_state, nphase, ncomp, &
                    PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, IPLIKE_GRAD_SOU, &
                    Velocity_U_Source_CV, Velocity_U_Source, Component, &
                    NCOLACV, FINACV, COLACV, MIDACV, &
                    small_FINACV, small_COLACV, small_MIDACV, &
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
                    StorageIndexes=StorageIndexes )

               if( have_option( '/material_phase[0]/multiphase_properties/capillary_pressure' ) )then
                  call calculate_capillary_pressure( state, packed_state, .true.)
               end if


               velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
               pressure_field=>extract_scalar_field(state(1),"Pressure")



               CALL FORCE_BAL_CTY_ASSEM_SOLVE( state, packed_state, &
                    velocity_field, pressure_field, &
                    NDIM, NPHASE, NCOMP, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
                    U_ELE_TYPE, P_ELE_TYPE, &
                    U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
                    U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN,&
                    STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
                    U_SNLOC, P_SNLOC, CV_SNLOC, &
!!$
                    Material_Absorption_Stab, Material_Absorption, Velocity_Absorption, Velocity_U_Source, Velocity_U_Source_CV, &
                    DRhoDPressure, IDIVID_BY_VOL_FRAC, FEM_VOL_FRAC, &
                    dt, &
!!$
                    NCOLC, FINDC, COLC, & ! C sparsity - global cty eqn 
                    NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, MIDDGM_PHA, &! Force balance sparsity
                    NCOLELE, FINELE, COLELE, & ! Element connectivity.
                    NCOLCMC, FINDCMC, COLCMC, MIDCMC, & ! pressure matrix for projection method
                    NCOLACV, FINACV, COLACV, MIDACV, & ! For CV discretisation method
                    size(small_colacv),small_FINACV, small_COLACV, small_MIDACV, &
                    NLENMCY, NCOLMCY, FINMCY, COLMCY, MIDMCY, & ! Force balance plus cty multi-phase eqns
                    NCOLCT, FINDCT, COLCT, & ! CT sparsity - global cty eqn.
                    CV_ELE_TYPE, &
!!$
                    v_disopt, v_dg_vel_int_opt, v_theta, &
                    SUF_SIG_DIAGTEN_BC, &
                    ScalarField_Source_Store, ScalarField_Absorption, Porosity_field%val, &
!!$
                    NCOLM, FINDM, COLM, MIDM, & ! Sparsity for the CV-FEM
                    XU_NLOC, XU_NDGLN, &
!!$
                    Momentum_Diffusion, Momentum_Diffusion_Vol, THERM_U_DIFFUSION, THERM_U_DIFFUSION_VOL, &
                    opt_vel_upwind_coefs, nopt_vel_upwind_coefs, &
                    igot_theta_flux, scvngi_theta, volfra_use_theta_flux, &
                    sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j, &
                    in_ele_upwind, dg_ele_upwind, &
!!$
                    NOIT_DIM, & ! This need to be removed as it is already deprecated
!!$
                    iplike_grad_sou, plike_grad_sou_coef, plike_grad_sou_grad, &
                    scale_momentum_by_volume_fraction,&
                    StorageIndexes=StorageIndexes )
!!$ Calculate Density_Component for compositional
               if( have_component_field ) &
                    call Calculate_Component_Rho( state, packed_state, &
                    ncomp, nphase, cv_nonods )

            end if Conditional_ForceBalanceEquation

            Conditional_PhaseVolumeFraction: if ( solve_PhaseVolumeFraction ) then
               call VolumeFraction_Assemble_Solve( state, packed_state, &
                    NCOLACV, FINACV, COLACV, MIDACV, &
                    small_FINACV, small_COLACV, small_MIDACV, &
                    block_to_global_acv, global_dense_block_acv, &
                    NCOLCT, FINDCT, COLCT, &
                    CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
                    CV_ELE_TYPE, &
                    NPHASE, &
                    CV_NLOC, U_NLOC, X_NLOC,  &
                    CV_NDGLN, X_NDGLN, U_NDGLN, &
                    CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
!!$
                    MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
!!$
                    v_disopt, v_dg_vel_int_opt, dt, v_theta, v_beta, &
                    SUF_SIG_DIAGTEN_BC, &
                    DRhoDPressure, &
                    ScalarField_Source_Store, ScalarField_Absorption, Porosity_field%val, &
!!$
                    NDIM, &
                    NCOLM, FINDM, COLM, MIDM, &
                    XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
!!$
                    opt_vel_upwind_coefs, nopt_vel_upwind_coefs, &
                    Density_FEMT, &
                    igot_theta_flux,scvngi_theta, volfra_use_theta_flux, &
                    in_ele_upwind, dg_ele_upwind, &
!!$                    
                    NOIT_DIM, & ! This need to be removed as it is already deprecated
!!$
                    option_path = '/material_phase[0]/scalar_field::PhaseVolumeFraction', &
                    mass_ele_transp = mass_ele,&
                    theta_flux=sum_theta_flux, one_m_theta_flux=sum_one_m_theta_flux, &
                    theta_flux_j=sum_theta_flux_j, one_m_theta_flux_j=sum_one_m_theta_flux_j,&
                    StorageIndexes=StorageIndexes )

            end if Conditional_PhaseVolumeFraction


!!$ Starting loop over components
            sum_theta_flux = 0. ; sum_one_m_theta_flux = 0. ; sum_theta_flux_j = 0. ; sum_one_m_theta_flux_j = 0. ; ScalarField_Source_Component = 0.

            velocity_field=>extract_tensor_field(packed_state,"PackedVelocity")
            saturation_field=>extract_tensor_field(packed_state,"PackedPhaseVolumeFraction")
            old_saturation_field=>extract_tensor_field(packed_state,"PackedOldPhaseVolumeFraction")
               

            Conditional_Components:if( have_component_field ) then

               Loop_Components: do icomp = 1, ncomp

                  tracer_field=>extract_tensor_field(multicomponent_state(icomp),"PackedComponentMassFraction")
                  density_field=>extract_tensor_field(multicomponent_state(icomp),"PackedComponentDensity",stat)

!!$ Computing the absorption term for the multi-components equation
                  call Calculate_ComponentAbsorptionTerm( state, packed_state, &
                       icomp, cv_ndgln, & 
                       D_s%val, Porosity_field%val, mass_ele, &
                       Component_Absorption )

                  Conditional_SmoothAbsorption: if( have_option( '/material_phase[' // int2str( nstate - ncomp ) // &
                       ']/is_multiphase_component/KComp_Sigmoid' ) .and. nphase > 1 ) then
                     do cv_nodi = 1, cv_nonods
                        if( PhaseVolumeFraction( cv_nodi ) > 0.95 ) then
                           do iphase = 1, nphase
                              do jphase = min( iphase + 1, nphase ), nphase
                                 Component_Absorption( cv_nodi, iphase, jphase ) = &
                                      Component_Absorption( cv_nodi, iphase, jphase ) * max( 0.01, &
                                      20. * ( 1. - PhaseVolumeFraction( cv_nodi ) ) )
                              end do
                           end do
                        end if

                    !    if( PhaseVolumeFraction( cv_nodi ) > 0.90 ) then
                    !       do iphase = 1, nphase
                    !          do jphase = min( iphase + 1, nphase ), nphase
                    !             Component_Absorption( cv_nodi, iphase, jphase ) = &
                    !                  Component_Absorption( cv_nodi, iphase, jphase ) * max( 0.00001, &
                    !                  20. * ( 1. - (PhaseVolumeFraction( cv_nodi )-0.05)  ) )
                    !          end do
                    !       end do
                    !    end if

                     end do
                  end if Conditional_SmoothAbsorption

!!$ Computing diffusion term for the component conservative equation:
                  call Calculate_ComponentDiffusionTerm( state, packed_state, &
                       mat_ndgln, u_ndgln, x_ndgln, &
                       u_ele_type, p_ele_type, ncomp_diff_coef, comp_diffusion_opt, &
                       Component_Diffusion_Operator_Coefficient( icomp, :, : ), &
                       Component_Diffusion ,&
                       StorageIndexes=StorageIndexes )

!!$ NonLinear iteration for the components advection:
                  Loop_NonLinearIteration_Components: do its2 = 1, NonLinearIteration_Components
                     comp_use_theta_flux = .false. ; comp_get_theta_flux = .true.

                     call INTENERGE_ASSEM_SOLVE( state, multicomponent_state(icomp), &
                          tracer_field,velocity_field,density_field,&
                          NCOLACV, FINACV, COLACV, MIDACV, & ! CV sparsity pattern matrix
                          SMALL_FINACV, SMALL_COLACV, small_MIDACV,&
                          block_to_global_acv, global_dense_block_acv, &
                          NCOLCT, FINDCT, COLCT, &
                          CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
                          U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE,  &
                          NPHASE,  &
                          CV_NLOC, U_NLOC, X_NLOC,  &
                          CV_NDGLN, X_NDGLN, U_NDGLN, &
                          CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
!!$
                          Component( ( icomp - 1 ) * nphase * cv_nonods + 1 : icomp * nphase * cv_nonods ), &
                          Component_Old( ( icomp - 1 ) * nphase * cv_nonods + 1 : icomp * nphase * cv_nonods ), &
!!$
                          MAT_NLOC, MAT_NDGLN, MAT_NONODS, Component_Diffusion, 0, THERM_U_DIFFUSION, THERM_U_DIFFUSION_VOL,&
                          v_disopt, v_dg_vel_int_opt, dt, v_theta, v_beta, &
                          SUF_SIG_DIAGTEN_BC,&
                          DRhoDPressure, &
                          Component_Source, Component_Absorption, Porosity_field%val, &
!!$
                          NDIM,  &
                          NCOLM, FINDM, COLM, MIDM, &
                          XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
!!$
                          opt_vel_upwind_coefs, nopt_vel_upwind_coefs, &
                          Component_FEMT( ( icomp - 1 ) * nphase * cv_nonods + 1 : icomp * nphase * cv_nonods ), &
                          Density_FEMT, &
                          igot_t2, PhaseVolumeFraction, PhaseVolumeFraction_Old, igot_theta_flux, scvngi_theta, &
                          comp_get_theta_flux, comp_use_theta_flux, &
                          theta_gdiff, &
                          in_ele_upwind, dg_ele_upwind, &
!!$
                          NOIT_DIM, & ! This need to be removed as it is already deprecated
!!$
                          Mean_Pore_CV, &
                          mass_ele_transp = dummy_ele, &
                          thermal = .false.,& ! the false means that we don't add an extra source term
                          theta_flux=theta_flux, one_m_theta_flux=one_m_theta_flux, theta_flux_j=theta_flux_j, one_m_theta_flux_j=one_m_theta_flux_j,&
                          StorageIndexes=StorageIndexes, icomp=icomp, saturation=saturation_field )

                   Component( ( icomp - 1 ) * nphase * cv_nonods + 1 : icomp * nphase * cv_nonods )  &
!                       =min(  max(Component( ( icomp - 1 ) * nphase * cv_nonods + 1 : icomp * nphase * cv_nonods ),0.0), 0.95) 
                       =min(  max(Component( ( icomp - 1 ) * nphase * cv_nonods + 1 : icomp * nphase * cv_nonods ),0.0), 1.0) 

                  end do Loop_NonLinearIteration_Components

                  sum_theta_flux = sum_theta_flux + theta_flux
                  sum_one_m_theta_flux = sum_one_m_theta_flux + one_m_theta_flux

                  sum_theta_flux_j = sum_theta_flux_j + theta_flux_j
                  sum_one_m_theta_flux_j = sum_one_m_theta_flux_j + one_m_theta_flux_j


                  ! We have divided through by density 
                  do cv_inod=1,cv_nonods
                     do iphase=1,nphase
                  ScalarField_Source_Component((iphase-1)*cv_nonods+cv_inod) = ScalarField_Source_Component((iphase-1)*cv_nonods+cv_inod) + THETA_GDIFF(iphase,cv_inod)
                     end do
                  end do

               end do Loop_Components

!! Temporal working array:
         DO CV_INOD = 1, CV_NONODS
            DO IPHASE = 1, NPHASE
               DO ICOMP = 1, NCOMP
                  MFC_s%val(ICOMP, IPHASE, CV_INOD) = &
                      COMPONENT( CV_INOD + ( IPHASE - 1 ) * CV_NONODS + ( ICOMP - 1 ) * NPHASE * CV_NONODS )
               END DO
            END DO
         END DO


               if( have_option( '/material_phase[' // int2str( nstate - ncomp ) // & 
                    ']/is_multiphase_component/Comp_Sum2One/Enforce_Comp_Sum2One' ) ) then
                  ! Initially clip and then ensure the components sum to unity so we don't get surprising results...
!                  DO I = 1, CV_NONODS * NPHASE * NCOMP
!                     COMPONENT( I ) = MIN( MAX( COMPONENT( I ), 0. ), 1. )
!                  END DO
                   MFC_s % val = min ( max ( MFC_s % val, 0.0), 1.0)

                  ALLOCATE( RSUM( NPHASE ) )
                  DO CV_INOD = 1, CV_NONODS
                     RSUM = 0.0
                     DO IPHASE = 1, NPHASE
                        RSUM( IPHASE ) = SUM (MFC_s % val (:, IPHASE, CV_INOD) )
!                        DO ICOMP = 1, NCOMP
!                           RSUM( IPHASE ) = RSUM( IPHASE ) + &
!                                COMPONENT( CV_INOD + ( IPHASE - 1 ) * CV_NONODS + ( ICOMP - 1 ) * NPHASE * CV_NONODS )
!                        END DO
                     END DO
                     DO IPHASE = 1, NPHASE
                        MFC_s % val (:, IPHASE, CV_INOD) = MFC_s % val (:, IPHASE, CV_INOD) / RSUM( IPHASE )
!                        DO ICOMP = 1, NCOMP
!                           COMPONENT( CV_INOD + ( IPHASE - 1 ) * CV_NONODS + ( ICOMP - 1 ) * NPHASE * CV_NONODS ) = &
!                                COMPONENT( CV_INOD + ( IPHASE - 1 ) * CV_NONODS + ( ICOMP - 1 ) * NPHASE * CV_NONODS ) / RSUM( IPHASE )
!                        END DO
                     END DO
                  END DO
                  DEALLOCATE( RSUM )
               end if

               DO ICOMP = 1, NCOMP

                  call Calculate_ComponentAbsorptionTerm( state, packed_state,&
                       icomp, cv_ndgln, & 
                       D_s%val, Porosity_field%val, mass_ele, &
                       Component_Absorption )

                  if( have_option( '/material_phase[' // int2str( nstate - ncomp ) // &
                       ']/is_multiphase_component/KComp_Sigmoid' ) .and. nphase > 1 ) then
                     do cv_nodi = 1, cv_nonods
                        if( PhaseVolumeFraction( cv_nodi ) > 0.95 ) then
                           do iphase = 1, nphase
                              do jphase = min( iphase + 1, nphase ), nphase
                                 Component_Absorption( cv_nodi, iphase, jphase ) = &
                                      Component_Absorption( cv_nodi, iphase, jphase ) * max( 0.01, &
                                      20. * ( 1. - SAT_s(1, cv_nodi ) ) )
                              end do
                           end do
                        end if
                     end do
                  end if

                  Loop_Phase_SourceTerm1: do iphase = 1, nphase
                     Loop_Phase_SourceTerm2: do jphase = 1, nphase
                        DO CV_NODI = 1, CV_NONODS
                           ScalarField_Source_Component( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) = &
                                ScalarField_Source_Component( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) - &
                                Component_Absorption( CV_NODI, IPHASE, JPHASE ) * &
!                                Component( CV_NODI + ( JPHASE - 1 ) * CV_NONODS + ( ICOMP - 1 ) * NPHASE * CV_NONODS ) / &
                                MFC_s%val(ICOMP, JPHASE, CV_NODI) / &
                                DC_s%val( icomp, iphase, cv_nodi  )
                        END DO
                     end do Loop_Phase_SourceTerm2
                  end do Loop_Phase_SourceTerm1

                  ! For compressibility
                  DO IPHASE = 1, NPHASE
                     DO CV_NODI = 1, CV_NONODS
                        ScalarField_Source_Component( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) = &
                             ScalarField_Source_Component( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) &
!                             + Mean_Pore_CV( CV_NODI ) * Component_Old( CV_NODI + ( IPHASE - 1 ) * CV_NONODS + ( ICOMP - 1 ) * NPHASE * CV_NONODS ) &
                             + Mean_Pore_CV( CV_NODI ) * MFCOLD_s%val(ICOMP, IPHASE, CV_NODI) &
                             * ( DCOLD_s%val( ICOMP, IPHASE, CV_NODI ) - DC_s%val( ICOMP, IPHASE, CV_NODI) ) &
                             * OldSAT_s( IPHASE, CV_NONODS ) &
                             / ( DC_s%val( ICOMP, IPHASE, CV_NODI ) * DT )
                     END DO
                  END DO

               END DO ! ICOMP

               if( have_option( '/material_phase[' // int2str( nstate - ncomp ) // & 
                    ']/is_multiphase_component/Comp_Sum2One' ) .and. ( ncomp > 1 ) ) then
                  call Cal_Comp_Sum2One_Sou( packed_state, ScalarField_Source_Component, cv_nonods, nphase, ncomp, dt, its, &
                       NonLinearIteration, &
                       Mean_Pore_CV )
               end if

!! Temporal working array:
               DO CV_INOD = 1, CV_NONODS
                  DO IPHASE = 1, NPHASE
                     DO ICOMP = 1, NCOMP
                      COMPONENT( CV_INOD + ( IPHASE - 1 ) * CV_NONODS + ( ICOMP - 1 ) * NPHASE * CV_NONODS ) = &
                           MFC_s%val(ICOMP, IPHASE, CV_INOD)
                     END DO
                  END DO
               END DO

!!$               ! Update state memory
!!$               do icomp = 1, ncomp
!!$                  do iphase = 1, nphase
!!$                     Component_State => extract_scalar_field( state( icomp + nphase ), & 
!!$                          'ComponentMassFractionPhase' // int2str( iphase ) )
!!$                     Component_State % val = component( 1 + ( iphase - 1 ) * cv_nonods + ( icomp - 1 ) * &
!!$                          nphase * cv_nonods : iphase * cv_nonods + ( icomp - 1 ) * nphase * cv_nonods )
!!$                  end do
!!$               end do

            end if Conditional_Components

            !Check if the results are good so far and act in consequence, only does something if requested by the user
            call Adaptive_NonLinear(packed_state, reference_field, its,&
            Repeat_time_step, ExitNonLinearLoop,nonLinearAdaptTs,3)
            if (ExitNonLinearLoop) exit Loop_NonLinearIteration

         end do Loop_NonLinearIteration

        if (nonLinearAdaptTs) then
            !As the value of dt and acctim may have changed we retrieve their values
            !to make sure that everything is coherent
            call get_option( '/timestepping/current_time', acctim )
            call get_option( '/timestepping/timestep', dt)
            !If repeat timestep we don't want to adapt mesh or dump results
            if ( Repeat_time_step ) then
                itime = itime - 1
                cycle Loop_Time
            end if
        end if

         call set_option( '/timestepping/current_time', acctim )
         call set_option( '/timestepping/timestep', dt)

         current_time = acctim

!!$ Copying fields back to state:
!         call copy_into_state( state, & ! Copying main fields into state
!              PhaseVolumeFraction, Temperature, &
!              Component, ncomp, nphase, cv_ndgln )

!!$ Calculate diagnostic fields
         call calculate_diagnostic_variables( state, exclude_nonrecalculated = .true. )
         call calculate_diagnostic_variables_new( state, exclude_nonrecalculated = .true. )


         if (write_all_stats) call write_diagnostics( state, current_time, dt, itime )  ! Write stat file

         Conditional_TimeDump: if( ( mod( itime, dump_period_in_timesteps ) == 0 ) ) then

            dtime=dtime+1
            if (do_checkpoint_simulation(dtime)) then
               call checkpoint_simulation(state,cp_no=checkpoint_number,&
                    protect_simulation_name=.true.,file_type='.mpml')
               checkpoint_number=checkpoint_number+1
            end if

            if ( have_option( "/io/output_scalars_fem" ) ) then
                 !As SAT_s is pointing into state, we need to use a backup to print the FE saturation.
                 allocate(dummy_to_print_FEM(size(SAT_s,1),size(SAT_s,2)))
                 dummy_to_print_FEM = SAT_s!<=create backup
                 call copy_into_state( state, & ! Copying main fields into state
                 FESAT_s, Temperature_FEMT, &
                 Component_FEMT, ncomp, nphase, cv_ndgln )
            end if

            call get_option( '/timestepping/current_time', current_time ) ! Find the current time 

            if (.not. write_all_stats)call write_diagnostics( state, current_time, dt, itime/dump_period_in_timesteps )  ! Write stat file
            not_to_move_det_yet = .false. ; dump_no = itime/dump_period_in_timesteps ! Sync dump_no with itime
            call write_state( dump_no, state ) ! Now writing into the vtu files
           
            if ( have_option( "/io/output_scalars_fem" ) ) then
                 SAT_s = dummy_to_print_FEM!<= retrieve backup
                 deallocate(dummy_to_print_FEM)

                 call copy_into_state( state, & ! Copying main fields into state
                 SAT_s, Temperature, &
                 Component, ncomp, nphase, cv_ndgln )
            end if
         end if Conditional_TimeDump

!!$! ******************
!!$! *** Mesh adapt ***
!!$! ******************

         do_reallocate_fields = .false.
         Conditional_Adaptivity_ReallocatingFields: if( have_option( '/mesh_adaptivity/hr_adaptivity') ) then
            if( have_option( '/mesh_adaptivity/hr_adaptivity/period_in_timesteps') ) then
               call get_option( '/mesh_adaptivity/hr_adaptivity/period_in_timesteps', &
                    adapt_time_steps, default=5 )
            end if
            if( mod( itime, adapt_time_steps ) == 0 ) do_reallocate_fields = .true.
         elseif( have_option( '/mesh_adaptivity/prescribed_adaptivity' ) ) then
            if( do_adapt_state_prescribed( current_time ) ) do_reallocate_fields = .true.
         end if Conditional_Adaptivity_ReallocatingFields

         new_mesh = do_reallocate_fields

         Conditional_ReallocatingFields: if( do_reallocate_fields ) then

            Conditional_Adaptivity: if( have_option( '/mesh_adaptivity/hr_adaptivity ') ) then

               Conditional_Adapt_by_TimeStep: if( mod( itime, adapt_time_steps ) == 0 ) then

                  ! linearise compositional fields:
                  Conditional_Components_Linearisation2: if ( ncomp > 1 ) then

                     do icomp = 1, ncomp
                        do iphase = 1, nphase

                           Component_State => extract_scalar_field( state( icomp + nphase ), & 
                                'ComponentMassFractionPhase' // int2str( iphase ) )
                           if (.not. have_option(trim(Component_State%option_path)//"/prognostic/consistent_interpolation")) then

                              call Updating_Linearised_Components( totele, ndim, cv_nloc, cv_nonods, cv_ndgln,& 
                           component ( 1 + ( iphase - 1 ) * cv_nonods + ( icomp - 1 ) * &
                                nphase * cv_nonods : iphase * cv_nonods + ( icomp - 1 ) * nphase * cv_nonods ) )
                           Component_State % val = component( 1 + ( iphase - 1 ) * cv_nonods + ( icomp - 1 ) * &
                                nphase * cv_nonods : iphase * cv_nonods + ( icomp - 1 ) * nphase * cv_nonods )
                        end if
                        end do
                     end do
                  end if Conditional_Components_Linearisation2


                  call pre_adapt_tasks( sub_state )

                  call qmesh( state, metric_tensor )

                  if( have_option( '/io/stat/output_before_adapts' ) ) call write_diagnostics( state, current_time, dt, &
                       itime, not_to_move_det_yet = .true. )
!!$                  if( have_option( '/io/stat/output_before_adapts' ) ) call write_diagnostics( state, current_time, dt, &
!!$                       timestep, not_to_move_det_yet = .true. )

                  call run_diagnostics( state )

                  call adapt_state( state, metric_tensor, suppress_reference_warnings=.true. )

                  call update_state_post_adapt( state, metric_tensor, dt, sub_state, nonlinear_iterations, &
                       nonlinear_iterations_adapt )

                  if( have_option( '/io/stat/output_after_adapts' ) ) call write_diagnostics( state, current_time, dt, &
                       itime, not_to_move_det_yet = .true. )
!!$                  if( have_option( '/io/stat/output_after_adapts' ) ) call write_diagnostics( state, current_time, dt, &
!!$                       timestep, not_to_move_det_yet = .true. )

                  call run_diagnostics( state )
               end if Conditional_Adapt_by_TimeStep

            elseif( have_option( '/mesh_adaptivity/prescribed_adaptivity' ) ) then !!$ Conditional_Adaptivity:

               Conditional_Adapt_by_Time: if( do_adapt_state_prescribed( current_time ) ) then

                  call pre_adapt_tasks( sub_state )

                  if( have_option( '/io/stat/output_before_adapts' ) ) call write_diagnostics( state, current_time, dt, &
                       timestep, not_to_move_det_yet = .true. )

                  call run_diagnostics( state )

                  call adapt_state_prescribed( state, current_time )

                  call update_state_post_adapt( state, metric_tensor, dt, sub_state, nonlinear_iterations, &
                       nonlinear_iterations_adapt)

                  if(have_option( '/io/stat/output_after_adapts' ) ) call write_diagnostics( state, current_time, dt, &
                       timestep, not_to_move_det_yet = .true. )

                  call run_diagnostics( state )

               end if Conditional_Adapt_by_Time

               not_to_move_det_yet = .false.

            end if Conditional_Adaptivity

            call nullify( packed_state )
            call deallocate(packed_state)
            call nullify( multicomponent_state )
            call deallocate(multicomponent_state)
            call pack_multistate(state,packed_state,&
                 multiphase_state,multicomponent_state)

!        !If we are using adaptive time stepping, backup_state needs also to be redone
!        if (nonLinearAdaptTs) then
!            call deallocate(backup_state)
!            call pack_multistate(state,backup_state,&
!                 multiphase_state,multicomponent_state)
!        end if

            !The storaged variables must be recalculated
            call Clean_Storage(state, StorageIndexes)

!!$ Deallocating array variables:
            deallocate( &
!!$ Node glabal numbers
                 cv_sndgln, p_sndgln, u_sndgln, &
!!$ Sparsity patterns
                 finacv, colacv, midacv,&
                 small_finacv, small_colacv, small_midacv, &
                 finmcy, colmcy, midmcy, &
                 block_to_global_acv, global_dense_block_acv, &
                 finele, colele, midele, findgm_pha, coldgm_pha, middgm_pha, findct, &
                 colct, findc, colc, findcmc, colcmc, midcmc, findm, &
                 colm, midm, &
!!$ Defining element-pair type and discretisation options and coefficients
                 opt_vel_upwind_coefs, &
!!$ For output:
                 PhaseVolumeFraction_FEMT, Temperature_FEMT, Density_FEMT, &
                 Component_FEMT, Mean_Pore_CV, SumConc_FEMT, Dummy_PhaseVolumeFraction_FEMT, &
!!$ Variables used in the diffusion-like term: capilarity and surface tension:
                 plike_grad_sou_grad, plike_grad_sou_coef, &
!!$ Working arrays
                 Temperature, PhaseVolumeFraction, SAT_s,oldsat_s, &
                 Component, &
                 Temperature_Old, &
                 PhaseVolumeFraction_Old, Component_Old, &
                 DRhoDPressure, &
                 Velocity_U_Source, Velocity_U_Source_CV, Temperature_Source, PhaseVolumeFraction_Source, &
                 ScalarField_Source, Component_Source, ScalarAdvectionField_Source, &
                 suf_sig_diagten_bc, &
                 theta_gdiff,  ScalarField_Source_Store, ScalarField_Source_Component, &
                 mass_ele, dummy_ele, &
                 Permeability, Material_Absorption, Material_Absorption_Stab, &
                 Velocity_Absorption, ScalarField_Absorption, Component_Absorption, Temperature_Absorption, &
                 Component_Diffusion_Operator_Coefficient, &
                 Momentum_Diffusion, Momentum_Diffusion_Vol, ScalarAdvectionField_Diffusion, &
                 Component_Diffusion, &
                 theta_flux, one_m_theta_flux, theta_flux_j, one_m_theta_flux_j, sum_theta_flux, &
                 sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j )


!!$  Compute primary scalars used in most of the code
            call Get_Primary_Scalars( state, &         
                 nphase, nstate, ncomp, totele, ndim, stotel, &
                 u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
                 x_snloc, cv_snloc, u_snloc, p_snloc, &
                 cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, dx )
!!$ Calculating Global Node Numbers
            allocate( cv_sndgln( stotel * cv_snloc ), p_sndgln( stotel * p_snloc ), &
                 u_sndgln( stotel * u_snloc ) )

  !          x_ndgln_p1 = 0 ; x_ndgln = 0 ; cv_ndgln = 0 ; p_ndgln = 0 ; mat_ndgln = 0 ; u_ndgln = 0 ; xu_ndgln = 0 ; &
                 cv_sndgln = 0 ; p_sndgln = 0 ; u_sndgln = 0

            call Compute_Node_Global_Numbers( state, &
                 totele, stotel, x_nloc, x_nloc_p1, cv_nloc, p_nloc, u_nloc, xu_nloc, &
                 cv_snloc, p_snloc, u_snloc, &
                 cv_ndgln, u_ndgln, p_ndgln, x_ndgln, x_ndgln_p1, xu_ndgln, mat_ndgln, &
                 cv_sndgln, p_sndgln, u_sndgln )
!!$
!!$ Computing Sparsity Patterns Matrices
!!$

!!$ Defining lengths and allocating space for the matrices
            call Defining_MaxLengths_for_Sparsity_Matrices( ndim, nphase, totele, u_nloc, cv_nloc, cv_nonods, &
                 mx_nface_p1, mxnele, mx_nct, mx_nc, mx_ncolcmc, mx_ncoldgm_pha, mx_ncolmcy, &
                 mx_ncolacv, mx_ncolm )
            nlenmcy = u_nonods * nphase * ndim + cv_nonods
            allocate( finacv( cv_nonods * nphase + 1 ), colacv( mx_ncolacv ), midacv( cv_nonods * nphase ), &
                 finmcy( nlenmcy + 1 ), colmcy( mx_ncolmcy ), midmcy( nlenmcy ), &
                 finele( totele + 1 ), colele( mxnele ), midele( totele ), &
                 findgm_pha( u_nonods * nphase * ndim + 1 ), coldgm_pha( mx_ncoldgm_pha ), &
                 middgm_pha( u_nonods * nphase * ndim ), &
                 findct( cv_nonods + 1 ), colct( mx_nct ), &
                 findc( u_nonods + 1 ), colc( mx_nc ), &
                 findcmc( cv_nonods + 1 ), colcmc( mx_ncolcmc ), midcmc( cv_nonods ), &
                 findm( cv_nonods + 1 ), colm( mx_ncolm ), midm( cv_nonods ) )

            allocate( global_dense_block_acv (nphase,cv_nonods) )
                 finacv = 0 ; colacv = 0 ; midacv = 0 ; finmcy = 0 ; colmcy = 0 ; midmcy = 0 ; finele = 0 ; &
                 colele = 0 ; midele = 0 ; findgm_pha = 0 ; coldgm_pha = 0 ; middgm_pha = 0 ; findct = 0 ; &
                 colct = 0 ; findc = 0 ; colc = 0 ; findcmc = 0 ; colcmc = 0 ; midcmc = 0 ; findm = 0 ; &
                 colm = 0 ; midm = 0

!!$ Defining element-pair type 
            call Get_Ele_Type( x_nloc, cv_ele_type, p_ele_type, u_ele_type, &
                 mat_ele_type, u_sele_type, cv_sele_type )

!!$ Sparsity Patterns Matrices 
            call Get_Sparsity_Patterns( state, &
!!$ CV multi-phase eqns (e.g. vol frac, temp)
                 mx_ncolacv, ncolacv, finacv, colacv, midacv, &
                 small_finacv, small_colacv, small_midacv, &
                 block_to_global_acv, global_dense_block_acv, &
!!$ Force balance plus cty multi-phase eqns
                 nlenmcy, mx_ncolmcy, ncolmcy, finmcy, colmcy, midmcy, &
!!$ Element connectivity
                 mxnele, ncolele, midele, finele, colele, &
!!$ Force balance sparsity
                 mx_ncoldgm_pha, ncoldgm_pha, coldgm_pha, findgm_pha, middgm_pha, &
!!$ CT sparsity - global continuity eqn
                 mx_nct, ncolct, findct, colct, &
!!$ C sparsity operating on pressure in force balance
                 mx_nc, ncolc, findc, colc, &
!!$ pressure matrix for projection method
                 mx_ncolcmc, ncolcmc, findcmc, colcmc, midcmc, &
!!$ CV-FEM matrix
                 mx_ncolm, ncolm, findm, colm, midm, mx_nface_p1 )

            call temp_mem_hacks()

!!$ Allocating space for various arrays:
            allocate( &
!!$
                 Temperature( nphase * cv_nonods ), &
                 PhaseVolumeFraction( nphase * cv_nonods ), SAT_s( nphase , cv_nonods ), Component( nphase * cv_nonods * ncomp ), &
                 oldSAT_s( nphase , cv_nonods ), &
                 DRhoDPressure( nphase, cv_nonods ), &
!!$
                 Temperature_Old( nphase * cv_nonods ), &
                 PhaseVolumeFraction_Old( nphase * cv_nonods ), Component_Old( nphase * cv_nonods * ncomp ), &
!!$             
                 suf_sig_diagten_bc( stotel * cv_snloc * nphase, ndim ), &
                 PhaseVolumeFraction_FEMT( cv_nonods * nphase ), Temperature_FEMT( cv_nonods * nphase ), &
                 Density_FEMT( cv_nonods * nphase ), Component_FEMT( cv_nonods * nphase * ncomp ), &
                 Mean_Pore_CV( cv_nonods ), SumConc_FEMT( cv_nonods * ncomp ), &
                 Dummy_PhaseVolumeFraction_FEMT( cv_nonods * nphase ), dummy_ele( totele ), mass_ele( totele ), &
!!$
                 Temperature_Source( cv_nonods * nphase ), &
                 PhaseVolumeFraction_Source( cv_nonods * nphase ), Velocity_U_Source( u_nonods * nphase * ndim ), &
                 Velocity_U_Source_CV( cv_nonods * nphase * ndim ), Component_Source( cv_nonods * nphase ), &
                 ScalarField_Source( cv_nonods * nphase ), ScalarAdvectionField_Source( cv_nonods * nphase ), &
!!$
                 Permeability( totele, ndim, ndim ), &
!!$
                 Material_Absorption( mat_nonods, ndim * nphase, ndim * nphase ), &
                 Velocity_Absorption( mat_nonods, ndim * nphase, ndim * nphase ), &
                 Material_Absorption_Stab( mat_nonods, ndim * nphase, ndim * nphase ), & 
                 ScalarField_Absorption( cv_nonods, nphase, nphase ), Component_Absorption( cv_nonods, nphase, nphase ), &
                 Temperature_Absorption( cv_nonods, nphase, nphase ), &
                 Momentum_Diffusion( mat_nonods, ndim, ndim, nphase ), &
                 Momentum_Diffusion_Vol( mat_nonods, nphase ), &
                 ScalarAdvectionField_Diffusion( mat_nonods, ndim, ndim, nphase ), & 
                 Component_Diffusion( mat_nonods, ndim, ndim, nphase ), &
!!$ Variables used in the diffusion-like term: capilarity and surface tension:
                 plike_grad_sou_grad( cv_nonods * nphase ), &
                 plike_grad_sou_coef( cv_nonods * nphase ) )    
!!$
            Velocity_U_Source = 0. ; Velocity_Absorption = 0. ; Velocity_U_Source_CV = 0. 
            Momentum_Diffusion=0.
            Momentum_Diffusion_Vol=0.
!!$
            Temperature=0. ; Temperature_Source=0. ; 
            Temperature_FEMT=0. ; Temperature_Absorption=0.
!!$
            Component=0. ; Component_Source=0.
            Component_Diffusion=0. ; Component_Absorption=0.
!!$
            Permeability=0.
!!$
!!$
            PhaseVolumeFraction=0. ; PhaseVolumeFraction_Old=0. ; PhaseVolumeFraction_Source=0.
            PhaseVolumeFraction_FEMT=0. ; Dummy_PhaseVolumeFraction_FEMT=0.
!!$
            ScalarAdvectionField_Diffusion=0. ; ScalarField_Absorption=0.
            ScalarField_Source=0. ; ScalarAdvectionField_Source=0.
!!$
            Material_Absorption=0. ; Material_Absorption_Stab=0.
!!$
            plike_grad_sou_grad=0. ; plike_grad_sou_coef=0.
!!$
            suf_sig_diagten_bc=0.
!!$


!!$ Extracting Mesh Dependent Fields
            initialised = .true.
            call Extracting_MeshDependentFields_From_State( state, packed_state, initialised, &
                 SAT_s, PhaseVolumeFraction_Source, &
                 Component, Component_Source, &
                 Velocity_U_Source, Velocity_Absorption, &
                 Temperature,  Temperature_Source, &
                 Permeability )

            call get_var_from_packed_state(packed_state,PhaseVolumeFraction = SAT_s,&
                 OldPhaseVolumeFraction=OldSAT_s,FEPhaseVolumeFraction = FESAT_s )

if (have_component_field) then
        !######TEMPORARY CONVERSION FROM OLD PhaseVolumeFraction TO PACKED######
        do cv_inod = 1, size(SAT_s,2)
            do iphase = 1, size(SAT_s,1)
                phaseVolumeFraction(cv_inod +(iphase-1)*size(SAT_s,2)) = SAT_s(iphase,cv_inod)
                PhaseVolumeFraction_Old(cv_inod +(iphase-1)*size(SAT_s,2)) = OldSAT_s(iphase,cv_inod)
            end do
        end do
        !#############################################################
    end if


!!$ Dummy field used in the scalar advection option:
            Dummy_PhaseVolumeFraction_FEMT = 1.

            ncv_faces=CV_count_faces( packed_state, CV_ELE_TYPE, stotel, cv_sndgln, u_sndgln )


!!$
!!$ Initialising Absorption terms that do not appear in the schema
!!$
            ScalarField_Absorption = 0. ; Component_Absorption = 0. ; Temperature_Absorption = 0.


!!$ Variables that can be effectively deleted as they are not used anymore:
            noit_dim = 0

!!$ Computing shape function scalars
            igot_t2 = 0 ; igot_theta_flux = 0
            if( ncomp /= 0 )then
               igot_t2 = 1 ; igot_theta_flux = 1
            end if

            call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
                 cv_ngi, cv_ngi_short, scvngi_theta, sbcvngi, nface, .false. )

            allocate( theta_flux( nphase, scvngi_theta*cv_nloc*totele * igot_theta_flux ), &
                 one_m_theta_flux( nphase, scvngi_theta*cv_nloc*totele * igot_theta_flux ), & 
                 theta_flux_j( nphase, scvngi_theta*cv_nloc*totele * igot_theta_flux ), &
                 one_m_theta_flux_j( nphase, scvngi_theta*cv_nloc*totele * igot_theta_flux ), & 
                 sum_theta_flux( nphase, scvngi_theta*cv_nloc*totele * igot_theta_flux ), &
                 sum_one_m_theta_flux( nphase, scvngi_theta*cv_nloc*totele * igot_theta_flux ), &
                 sum_theta_flux_j( nphase, scvngi_theta*cv_nloc*totele * igot_theta_flux ), &
                 sum_one_m_theta_flux_j( nphase, scvngi_theta*cv_nloc*totele * igot_theta_flux ), &
                 theta_gdiff( nphase, cv_nonods ), ScalarField_Source_Store( cv_nonods * nphase ), &
                 ScalarField_Source_Component( cv_nonods * nphase ) )

            sum_theta_flux = 1. ; sum_one_m_theta_flux = 0.  
            sum_theta_flux_j = 1. ; sum_one_m_theta_flux_j = 0.  
            ScalarField_Source_Store=0. ; ScalarField_Source_Component=0.

            allocate( Component_Diffusion_Operator_Coefficient( ncomp, ncomp_diff_coef, nphase ) )  
            nopt_vel_upwind_coefs = mat_nonods * nphase * ndim * ndim * 2
            allocate( opt_vel_upwind_coefs( nopt_vel_upwind_coefs ) ) ; opt_vel_upwind_coefs = 0.

         end if Conditional_ReallocatingFields

!!$ Simple adaptive time stepping algorithm
         if ( have_option( '/timestepping/adaptive_timestep' ) ) then
            c = -66.6 ; minc = 0. ; maxc = 66.e6 ; ic = 66.e6
            call get_option( '/timestepping/adaptive_timestep/requested_cfl', rc )
            call get_option( '/timestepping/adaptive_timestep/minimum_timestep', minc, stat )
            call get_option( '/timestepping/adaptive_timestep/maximum_timestep', maxc, stat )
            call get_option( '/timestepping/adaptive_timestep/increase_tolerance', ic, stat )

            do iphase = 1, nphase
               ! requested cfl
               rc_field => extract_scalar_field( state( iphase ), 'RequestedCFL', stat )
               if ( stat == 0 ) rc = min( rc, minval( rc_field % val ) )
               ! max cfl
               cfl => extract_scalar_field( state( iphase ), 'CFLNumber' )
               c = max ( c, maxval( cfl % val ) )
            end do

            call get_option( '/timestepping/timestep', dt )
            dt = max( min( min( dt * rc / c, ic * dt ), maxc ), minc )
            call allmin(dt)
            call set_option( '/timestepping/timestep', dt )
         end if

         call set_boundary_conditions_values(state, shift_time=.true.)

      end do Loop_Time

!!$ Now deallocating arrays:
      deallocate( &
!!$ Node glabal numbers
           cv_sndgln, p_sndgln, u_sndgln, &
!!$ Sparsity patterns
           finacv, colacv, midacv,&
           small_finacv, small_colacv, small_midacv, &
           finmcy, colmcy, midmcy, &
           block_to_global_acv, global_dense_block_acv, &
           finele, colele, midele, findgm_pha, coldgm_pha, middgm_pha, findct, &
           colct, findc, colc, findcmc, colcmc, midcmc, findm, &
           colm, midm, &
!!$ Defining element-pair type and discretisation options and coefficients
           opt_vel_upwind_coefs, &
!!$ For output:
           PhaseVolumeFraction_FEMT, Temperature_FEMT, Density_FEMT, &
           Component_FEMT, Mean_Pore_CV, SumConc_FEMT, Dummy_PhaseVolumeFraction_FEMT, &
!!$ Variables used in the diffusion-like term: capilarity and surface tension:
           plike_grad_sou_grad, plike_grad_sou_coef, &
!!$ Working arrays
           Temperature, PhaseVolumeFraction, &
           Component, &
           Temperature_Old, &
           PhaseVolumeFraction_Old, Component_Old, &
           DRhoDPressure, FEM_VOL_FRAC, &
           Velocity_U_Source, Velocity_U_Source_CV, Temperature_Source, PhaseVolumeFraction_Source, &
           ScalarField_Source, Component_Source, ScalarAdvectionField_Source, &
           theta_gdiff,  ScalarField_Source_Store, ScalarField_Source_Component, &
           mass_ele, dummy_ele, &
           Permeability, Material_Absorption, Material_Absorption_Stab, &
           Velocity_Absorption, ScalarField_Absorption, Component_Absorption, Temperature_Absorption, &
           Component_Diffusion_Operator_Coefficient, &
           Momentum_Diffusion, Momentum_Diffusion_Vol, ScalarAdvectionField_Diffusion, &
           Component_Diffusion, &
           theta_flux, one_m_theta_flux, theta_flux_j, one_m_theta_flux_j, &
           sum_theta_flux, sum_one_m_theta_flux, sum_theta_flux_j, sum_one_m_theta_flux_j )

      return

      contains 

        
        subroutine temp_mem_hacks()
          
!!! routine puts various CSR sparsities into packed_state

          use sparse_tools
          
          type(csr_sparsity) :: sparsity
          type(scalar_field), pointer :: sfield

          integer ic

          sparsity=wrap(finele,midele,colm=colele,name='ElementConnectivity')
          call insert(packed_state,sparsity,'ElementConnectivity')

          sparsity=wrap(small_finacv,small_midacv,colm=small_colacv,name='SinglePhaseAdvectionSparsity')
          call insert(packed_state,sparsity,'SinglePhaseAdvectionSparsity')
          sparsity=wrap(finacv,midacv,colm=colacv,name='PackedAdvectionSparsity')
          call insert(packed_state,sparsity,'PackedAdvectionSparsity')
          sparsity=wrap(findc,colm=colc,name='CMatrixSparsity')
          call insert(packed_state,sparsity,'CMatrixSparsity')
          sparsity=wrap(findct,colm=colct,name='CTMatrixSparsity')
          call insert(packed_state,sparsity,'CTMatrixSparsity')
          sparsity=wrap(findcmc,colm=colcmc,name='CMCMatrixSparsity')
          call insert(packed_state,sparsity,'CMCMatrixSparsity')
          sparsity=wrap(findm,midm,colm=colm,name='CVFEMSparsity')
          call insert(packed_state,sparsity,'CVFEMSparsity')
          
          sfield=>extract_scalar_field(packed_state,"Pressure")
          sparsity=make_sparsity(sfield%mesh,sfield%mesh,&
               "PressureMassMatrixSparsity")
          call insert(packed_state,sparsity,"PressureMassMatrixSparsity")
          do ic=1,size(multicomponent_state)
             call insert(multicomponent_state(ic),sparsity,"PressureMassMatrixSparsity")
          end do 
          call deallocate(sparsity)

        end subroutine temp_mem_hacks


    end subroutine MultiFluids_SolveTimeLoop

    subroutine Updating_Linearised_Components( totele, ndim, cv_nloc, cv_nonods, cv_ndgln, &
         component )
      implicit none
      integer, intent( in ) :: totele, ndim, cv_nloc, cv_nonods
      integer, dimension( : ), intent( in ) :: cv_ndgln
      real, dimension ( : ), intent( inout ) :: component
!!$Local variables
      integer :: ele, cv_iloc, cv_nod
      real, dimension( : ), allocatable :: density_tmp, den_cv_nod

      allocate( density_tmp( cv_nonods ), den_cv_nod( cv_nloc ) )
      density_tmp = 0. ; den_cv_nod = 0.

      Conditional_CV_Number: if( cv_nloc == 6 .or. (cv_nloc == 10 .and. ndim==3) ) then ! P2 triangle or tet
         density_tmp = component
         Loop_Elements: do ele = 1, totele
            Loop_CV: do cv_iloc = 1, cv_nloc
               cv_nod = cv_ndgln( ( ele - 1 ) * cv_nloc + cv_iloc )
               den_cv_nod( cv_iloc ) = density_tmp( cv_nod )
            end do Loop_CV

            den_cv_nod( 2 ) = 0.5 * ( den_cv_nod( 1 ) + den_cv_nod( 3 ) )
            den_cv_nod( 4 ) = 0.5 * ( den_cv_nod( 1 ) + den_cv_nod( 6 ) )
            den_cv_nod( 5 ) = 0.5 * ( den_cv_nod( 3 ) + den_cv_nod( 6 ) )
            if( cv_nloc == 10 ) then
               den_cv_nod( 7 ) = 0.5 * ( den_cv_nod( 1 ) + den_cv_nod( 10 ) )
               den_cv_nod( 8 ) = 0.5 * ( den_cv_nod( 3 ) + den_cv_nod( 10 ) )
               den_cv_nod( 9 ) = 0.5 * ( den_cv_nod( 6 ) + den_cv_nod( 10 ) )
            end if

            Loop_CV2: do cv_iloc = 1, cv_nloc
               cv_nod = cv_ndgln( ( ele - 1 ) * cv_nloc + cv_iloc )
               component( cv_nod ) = den_cv_nod( cv_iloc )
            end do Loop_CV2

         end do Loop_Elements

      end if Conditional_CV_Number

      deallocate( density_tmp, den_cv_nod )

      return
    end subroutine Updating_Linearised_Components

   subroutine copy_packed_new_to_old(packed_state)
     type(state_type), intent(inout) :: packed_state
     
     type(scalar_field), pointer :: sfield, nsfield
     type(vector_field), pointer :: vfield, nvfield
     type(tensor_field), pointer :: tfield, ntfield

     integer :: i

     do i=1,size(packed_state%scalar_fields)
        sfield=>packed_state%scalar_fields(i)%ptr
        if (sfield%name(1:9)=="PackedOld") then
           nsfield=>extract_scalar_field(packed_state,"Packed"//sfield%name(10:))
           sfield%val=nsfield%val
        end if
     end do

     do i=1,size(packed_state%vector_fields)
        vfield=>packed_state%vector_fields(i)%ptr
        if (vfield%name(1:9)=="PackedOld") then
           nvfield=>extract_vector_field(packed_state,"Packed"//vfield%name(10:))
           vfield%val=nvfield%val
        end if
     end do
     
     do i=1,size(packed_state%tensor_fields)
        tfield=>packed_state%tensor_fields(i)%ptr
        if (tfield%name(1:9)=="PackedOld") then
           ntfield=>extract_tensor_field(packed_state,"Packed"//tfield%name(10:))
           tfield%val=ntfield%val
        end if
     end do

     sfield=>extract_scalar_field(packed_state,"OldFEPressure")
     nsfield=>extract_scalar_field(packed_state,"FEPressure")
     sfield%val=nsfield%val

     sfield=>extract_scalar_field(packed_state,"OldCVPressure")
     nsfield=>extract_scalar_field(packed_state,"CVPressure")
     sfield%val=nsfield%val

   end subroutine copy_packed_new_to_old



  end module multiphase_time_loop
