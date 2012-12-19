
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
                               topology_mesh_name
    use fldebug
    use state_module
    use fields
    use field_options
    use fields_allocates
    use spud
    use signal_vars
    use populate_state_module
    use vector_tools

!!$ Modules required by adaptivity
    use qmesh_module
    use adapt_state_module
    use adapt_state_prescribed_module!, only: do_adapt_state_prescribed, adapt_state_prescribed
    use populate_sub_state_module
    use fluids_module!, only: pre_adapt_tasks, update_state_post_adapt

!!$ Modules indigenous to the prototype Code
    use multiphase_1D_engine
    use spact
    use multiphase_EOS
    use shape_functions_Linear_Quadratic
    use Compositional_Terms
    use Copy_Outof_State
    use Copy_BackTo_State


    !use mapping_for_ocvfem
    !use matrix_operations
    !use shape_functions
    !use printout

    implicit none
    private
    public :: MultiFluids_SolveTimeLoop

  contains

    subroutine MultiFluids_SolveTimeLoop( state, &
         dt, nonlinear_iterations, nonlinear_iteration_tolerance, &
         dump_no )
      implicit none
      type( state_type ), dimension( : ), intent( inout ) :: state
      integer, intent( inout ) :: dump_no, nonlinear_iterations 
      real, intent( inout ) :: dt 
      real, intent( inout ) :: nonlinear_iteration_tolerance

!!$ Primary scalars
      integer :: nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods
      real :: dx
      logical :: is_overlapping

!!$ Node global numbers
      integer, dimension( : ), allocatable :: x_ndgln_p1, x_ndgln, cv_ndgln, p_ndgln, &
           mat_ndgln, u_ndgln, xu_ndgln, cv_sndgln, p_sndgln, u_sndgln

!!$ Sparsity patterns
      integer :: nlenmcy, mx_nface_p1, mx_ncolacv, mxnele, mx_ncoldgm_pha, &
           mx_ncolmcy, mx_nct, mx_nc, mx_ncolcmc, mx_ncolm, &
           ncolacv, ncolmcy, ncolele, ncoldgm_pha, ncolct, ncolc, ncolcmc, ncolm
      integer, dimension( : ), allocatable :: finacv, colacv, midacv, finmcy, colmcy, midmcy, &
           finele, colele, midele, findgm_pha, coldgm_pha, middgm_pha, findct, &
           colct, findc, colc, findcmc, colcmc, midcmc, findm, &
           colm, midm

!!$ Defining element-pair type and discretisation options and coefficients
      integer :: cv_ele_type, p_ele_type, u_ele_type, mat_ele_type, u_sele_type, cv_sele_type, &
           t_disopt, v_disopt, t_dg_vel_int_opt, u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
           comp_diffusion_opt, ncomp_diff_coef, in_ele_upwind, dg_ele_upwind, nopt_vel_upwind_coefs, &
           nits_flux_lim_t, nits_flux_lim_volfra, nits_flux_lim_comp
      logical :: volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, comp_get_theta_flux, &
           t_use_theta_flux, t_get_theta_flux, scale_momentum_by_volume_fraction
      real :: t_beta, v_beta, t_theta, v_theta, u_theta
      real, dimension( : ), allocatable :: opt_vel_upwind_coefs

!!$ Defining time- and nonlinear interations-loops variables
      integer :: itime, dump_period_in_timesteps, final_timestep, &
           NonLinearIteration, NonLinearIteration_Components
      real :: acctim, finish_time, current_time

!!$ Defining problem that will be solved
      logical :: have_temperature_field, have_component_field, have_extra_DiffusionLikeTerm, &
           solve_force_balance, solve_PhaseVolumeFraction

!!$ Defining solver options
      integer :: velocity_max_iterations, PhaseVolumeFraction_max_iterations

!!$ Shape function related fields:
      integer :: cv_ngi, cv_ngi_short, scvngi_theta, sbcvngi, nface, igot_t2, igot_theta_flux

!!$ For output:
      real, dimension( : ), allocatable :: PhaseVolumeFraction_FEMT, Temperature_FEMT, Density_FEMT, &
           Component_FEMT, Mean_Pore_CV, SumConc_FEMT, Dummy_PhaseVolumeFraction_FEMT
      type( scalar_field ), pointer :: Pressure_State, Temperature_State, Component_State

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
      integer, dimension( : ), allocatable :: PhaseVolumeFraction_BC_Spatial, Pressure_FEM_BC_Spatial, &
           Density_BC_Spatial, Component_BC_Spatial, Velocity_U_BC_Spatial, Temperature_BC_Spatial
      real, dimension( : ), allocatable :: xu, yu, zu, x, y, z, ug, vg, wg, &
           Velocity_U, Velocity_V, Velocity_W, Velocity_U_Old, Velocity_V_Old, Velocity_W_Old, &
           Velocity_NU, Velocity_NV, Velocity_NW, Velocity_NU_Old, Velocity_NV_Old, Velocity_NW_Old, &
           Pressure_FEM, Pressure_CV, Temperature, Density, Density_Component, PhaseVolumeFraction, &
           Component, U_Density, Pressure_FEM_Old, Pressure_CV_Old, Temperature_Old, Density_Old, &
           Density_Component_Old, PhaseVolumeFraction_Old, Component_Old, U_Density_Old, DRhoDPressure, &
           Porosity, &
           Velocity_U_Source, Velocity_U_Source_CV, Temperature_Source, PhaseVolumeFraction_Source, &
           ScalarField_Source, Component_Source, ScalarAdvectionField_Source, &
           PhaseVolumeFraction_BC, Pressure_FEM_BC, &
           Density_BC, Component_BC, Velocity_U_BC, Velocity_V_BC, Velocity_W_BC, Temperature_BC, &
           suf_u_bc_rob1, suf_v_bc_rob1, suf_w_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob2, suf_w_bc_rob2, &
           suf_t_bc_rob1, suf_t_bc_rob2, suf_vol_bc_rob1, suf_vol_bc_rob2, suf_comp_bc_rob1, suf_comp_bc_rob2, &
           theta_gdiff,  ScalarField_Source_Store, ScalarField_Source_Component, &
           mass_ele, dummy_ele, density_tmp, density_old_tmp

!!$
      real, dimension( :, :, : ), allocatable :: Permeability, Material_Absorption, Material_Absorption_Stab, &
           Velocity_Absorption, ScalarField_Absorption, Component_Absorption, Temperature_Absorption, &
!!$
           Component_Diffusion_Operator_Coefficient, Momentum_Diffusion_tmp
      real, dimension( :, :, :, : ), allocatable :: Momentum_Diffusion, ScalarAdvectionField_Diffusion, &
           Component_Diffusion, &
           theta_flux, one_m_theta_flux, sum_theta_flux, sum_one_m_theta_flux

!!$ Material_Absorption_Stab = u_abs_stab; Material_Absorption = u_absorb; ScalarField_Absorption = v_absorb
!!$ Component_Absorption = comp_absorb; ScalarAdvectionField_Absorption = t_absorb
!!$ Velocity_U_Source = u_source, ScalarField_Source = v_source, Component_Source = comp_source, #
!!$ ScalarAdvectionField_Source = Temperature_Source = t_source; Component_Diffusion_Operator_Coefficient = comp_diff_coef
!!$ Momentum_Diffusion = udiffusion; ScalarAdvectionField_Diffusion = tdiffusion, 
!!$ Component_Diffusion = comp_diffusion

      type( tensor_field ), pointer :: t_field
      character( len = option_path_len ) :: option_path
      integer :: stat, istate, iphase, jphase, icomp, its, its2, cv_nodi, adapt_time_steps, cv_inod
      real :: rsum

      real, dimension(:, :), allocatable :: DEN_CV_NOD
      integer :: CV_NOD, CV_NOD_PHA, CV_ILOC, ELE

!!$ Compute primary scalars used in most of the code
      call Get_Primary_Scalars( state, &         
           nphase, nstate, ncomp, totele, ndim, stotel, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
           x_snloc, cv_snloc, u_snloc, p_snloc, &
           cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, dx, &
           is_overlapping )

!!$ Calculating Global Node Numbers
      allocate( x_ndgln_p1( totele * x_nloc_p1 ), x_ndgln( totele * x_nloc ), cv_ndgln( totele * cv_nloc ), &
           p_ndgln( totele * p_nloc ), mat_ndgln( totele * mat_nloc ), u_ndgln( totele * u_nloc ), &
           xu_ndgln( totele * xu_nloc ), cv_sndgln( stotel * cv_snloc ), p_sndgln( stotel * p_snloc ), &
           u_sndgln( stotel * u_snloc ) )

      x_ndgln_p1 = 0 ; x_ndgln = 0 ; cv_ndgln = 0 ; p_ndgln = 0 ; mat_ndgln = 0 ; u_ndgln = 0 ; xu_ndgln = 0 ; &
           cv_sndgln = 0 ; p_sndgln = 0 ; u_sndgln = 0

      call Compute_Node_Global_Numbers( state, &
           is_overlapping, totele, stotel, x_nloc, x_nloc_p1, cv_nloc, p_nloc, u_nloc, xu_nloc, &
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

!!$ Allocating space for various arrays:
      allocate( xu( xu_nonods ), yu( xu_nonods ), zu( xu_nonods ), &
           x( x_nonods ), y( x_nonods ), z( x_nonods ), &
           ug( u_nonods * nphase ), vg( u_nonods * nphase ), wg( u_nonods * nphase ), &
!!$
           Velocity_U( u_nonods * nphase ), Velocity_V( u_nonods * nphase ), Velocity_W( u_nonods * nphase ), &
           Velocity_U_Old( u_nonods * nphase ), Velocity_V_Old( u_nonods * nphase ), Velocity_W_Old( u_nonods * nphase ), &
           Velocity_NU( u_nonods * nphase ), Velocity_NV( u_nonods * nphase ), Velocity_NW( u_nonods * nphase ), &
           Velocity_NU_Old( u_nonods * nphase ), Velocity_NV_Old( u_nonods * nphase ), Velocity_NW_Old( u_nonods * nphase ), &
!!$
           Pressure_FEM( cv_nonods ), Pressure_CV( cv_nonods ), &
           Temperature( nphase * cv_nonods ), Density( nphase * cv_nonods ), &
           Density_Component(nphase * cv_nonods * ncomp), &
           PhaseVolumeFraction( nphase * cv_nonods ), Component( nphase * cv_nonods * ncomp ), &
           U_Density( nphase * cv_nonods ), DRhoDPressure( nphase * cv_nonods ), &
!!$
           Pressure_FEM_Old( cv_nonods ), Pressure_CV_Old( cv_nonods ), &
           Temperature_Old( nphase * cv_nonods ), Density_Old( nphase * cv_nonods ), &
           Density_Component_Old(nphase * cv_nonods * ncomp), &
           PhaseVolumeFraction_Old( nphase * cv_nonods ), Component_Old( nphase * cv_nonods * ncomp ), &
           U_Density_Old( nphase * cv_nonods ), &
!!$
           PhaseVolumeFraction_BC_Spatial( stotel * nphase ), Pressure_FEM_BC_Spatial( stotel * nphase ), &
           Density_BC_Spatial( stotel * nphase ), Component_BC_Spatial( stotel * nphase ), &
           Velocity_U_BC_Spatial( stotel * nphase ), Temperature_BC_Spatial( stotel * nphase ), &
           PhaseVolumeFraction_BC( stotel * cv_snloc * nphase ), Pressure_FEM_BC( stotel * p_snloc * nphase ), &
           Density_BC( stotel * cv_snloc * nphase ), Temperature_BC( stotel * cv_snloc * nphase ), &
           Component_BC( stotel * cv_snloc * nphase * ncomp ), &
           Velocity_U_BC( stotel * u_snloc * nphase ), Velocity_V_BC( stotel * u_snloc * nphase ), &
           Velocity_W_BC( stotel * u_snloc * nphase ), Temperature_Source( cv_nonods * nphase ), &
           suf_u_bc_rob1( stotel * u_snloc * nphase ), suf_v_bc_rob1( stotel * u_snloc * nphase ), &
           suf_w_bc_rob1( stotel * u_snloc * nphase ), suf_u_bc_rob2( stotel * u_snloc * nphase ), &
           suf_v_bc_rob2( stotel * u_snloc * nphase ), suf_w_bc_rob2( stotel * u_snloc * nphase ), &
           suf_t_bc_rob1( stotel * cv_snloc * nphase ), suf_t_bc_rob2( stotel * cv_snloc * nphase ), &
           suf_vol_bc_rob1( stotel * cv_snloc * nphase ), suf_vol_bc_rob2( stotel * cv_snloc * nphase ), &
           suf_comp_bc_rob1( stotel * cv_snloc * nphase ), suf_comp_bc_rob2( stotel * cv_snloc * nphase ), &
!!$
           Porosity( totele ), &
           PhaseVolumeFraction_FEMT( cv_nonods * nphase ), Temperature_FEMT( cv_nonods * nphase ), &
           Density_FEMT( cv_nonods * nphase ), Component_FEMT( cv_nonods * nphase * ncomp ), &
           Mean_Pore_CV( cv_nonods ),  SumConc_FEMT( cv_nonods * ncomp ), &
           Dummy_PhaseVolumeFraction_FEMT( cv_nonods * nphase ), dummy_ele( totele ), mass_ele( totele ), &
!!$
           PhaseVolumeFraction_Source( cv_nonods * nphase ), Velocity_U_Source( u_nonods * nphase * ndim ), &
           Velocity_U_Source_CV( cv_nonods * nphase * ndim ), Component_Source( cv_nonods * nphase ), &
           ScalarField_Source( cv_nonods * nphase ), ScalarAdvectionField_Source( cv_nonods * nphase ), &
!!$
           Permeability( totele, ndim, ndim ), &
!!$
           Material_Absorption( mat_nonods, ndim * nphase, ndim * nphase ), &
           Velocity_Absorption( u_nloc * totele, ndim * nphase, ndim * nphase ), &
           Material_Absorption_Stab( mat_nonods, ndim * nphase, ndim * nphase ), & 
           ScalarField_Absorption( cv_nonods, nphase, nphase ), Component_Absorption( cv_nonods, nphase, nphase ), &
           Temperature_Absorption( cv_nonods, nphase, nphase ), &
           Momentum_Diffusion( mat_nonods, ndim, ndim, nphase ), &
           ScalarAdvectionField_Diffusion( mat_nonods, ndim, ndim, nphase ), & 
           Component_Diffusion( mat_nonods, ndim, ndim, nphase ), &
!!$ Variables used in the diffusion-like term: capilarity and surface tension:
           plike_grad_sou_grad( cv_nonods * nphase ), &
           plike_grad_sou_coef( cv_nonods * nphase ) )    
!!$
      xu=0. ; yu=0. ; zu=0.
      x=0. ; y=0. ; z=0.
      ug=0. ; vg=0. ; wg=0.
!!$
      Velocity_U=0. ; Velocity_V=0 ; Velocity_W=0.
      Velocity_U_Old=0. ; Velocity_V_Old=0. ; Velocity_W_Old=0.
      Velocity_NU=0. ; Velocity_NV=0. ; Velocity_NW=0.
      Velocity_NU_Old=0. ; Velocity_NV_Old=0. ; Velocity_NW_Old=0.
!!$
      Pressure_FEM=0. ; Pressure_CV=0.
      Temperature=0. ; Density=0.
      Density_Component=0.
      PhaseVolumeFraction=0. ; Component=0.
      U_Density=0. ; DRhoDPressure=0.
!!$
      Pressure_FEM_Old=0. ; Pressure_CV_Old=0.
      Temperature_Old=0. ; Density_Old=0.
      Density_Component_Old=0.
      PhaseVolumeFraction_Old=0. ; Component_Old=0.
      U_Density_Old=0.
!!$
      PhaseVolumeFraction_BC_Spatial=0 ; Pressure_FEM_BC_Spatial=0
      Density_BC_Spatial=0 ; Component_BC_Spatial=0
      Velocity_U_BC_Spatial=0 ; Temperature_BC_Spatial=0
      PhaseVolumeFraction_BC=0. ; Pressure_FEM_BC=0.
      Density_BC=0. ; Temperature_BC=0.
      Component_BC=0.
      Velocity_U_BC=0. ; Velocity_V_BC=0.
      Velocity_W_BC=0. ; Temperature_Source=0.
      suf_u_bc_rob1=0. ; suf_v_bc_rob1=0.
      suf_w_bc_rob1=0. ; suf_u_bc_rob2=0.
      suf_v_bc_rob2=0. ; suf_w_bc_rob2=0.
      suf_t_bc_rob1=0. ; suf_t_bc_rob2=0.
      suf_vol_bc_rob1=0. ; suf_vol_bc_rob2=0.
      suf_comp_bc_rob1=0. ; suf_comp_bc_rob2=0.
!!$
      Porosity=0.
      PhaseVolumeFraction_FEMT=0. ; Temperature_FEMT=0.
      Density_FEMT=0. ; Component_FEMT=0.
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
      ScalarAdvectionField_Diffusion=0.
      Component_Diffusion=0.
!!$
      plike_grad_sou_grad=0.
      plike_grad_sou_coef=0.
      iplike_grad_sou=0 

! dummy densities for testing
 allocate( density_tmp(cv_nonods) , density_old_tmp(cv_nonods) )
 density_tmp=0. ; density_old_tmp=0.

!!$ Extracting Mesh Dependent Fields
      initialised = .false.
      call Extracting_MeshDependentFields_From_State( state, initialised, &
           xu, yu, zu, x, y, z, &
           PhaseVolumeFraction, PhaseVolumeFraction_BC_Spatial, PhaseVolumeFraction_BC, PhaseVolumeFraction_Source, &
           Pressure_CV, Pressure_FEM, Pressure_FEM_BC_Spatial, Pressure_FEM_BC, &
           Density, Density_BC_Spatial, Density_BC, &
           Component, Component_BC_Spatial, Component_BC, Component_Source, &
           Velocity_U, Velocity_V, Velocity_W, Velocity_NU, Velocity_NV, Velocity_NW, &
           Velocity_U_BC_Spatial, Velocity_U_BC, Velocity_V_BC, Velocity_W_BC, Velocity_U_Source, Velocity_Absorption, &
           Temperature, Temperature_BC_Spatial, Temperature_BC, Temperature_Source, &
           Porosity, Permeability )


!print *, Component_BC_Spatial
!print *, '@@@@@'
!print *, Component_BC
!stop 777

!!$ Dummy field used in the scalar advection option:
      Dummy_PhaseVolumeFraction_FEMT = 1.

!!$
!!$ Initialising Robin boundary conditions --  this still need to be defined in the schema:
!!$
      suf_u_bc_rob1 = 0. ; suf_v_bc_rob1 = 0. ; suf_w_bc_rob1 = 0. ; suf_u_bc_rob2 = 0. ; suf_v_bc_rob2 = 0.
      suf_w_bc_rob2 = 0. ; suf_t_bc_rob1 = 0. ; suf_t_bc_rob2 = 0. ; suf_vol_bc_rob1 = 0. ; suf_vol_bc_rob2 = 0.
      suf_comp_bc_rob1 = 0. ; suf_comp_bc_rob2 = 0.

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

      allocate( theta_flux( totele * igot_theta_flux, cv_nloc, scvngi_theta, nphase ), &
           one_m_theta_flux( totele * igot_theta_flux, cv_nloc, scvngi_theta, nphase ), &
           sum_theta_flux( totele * igot_theta_flux, cv_nloc, scvngi_theta, nphase ), &
           sum_one_m_theta_flux( totele * igot_theta_flux, cv_nloc, scvngi_theta, nphase ), &
           theta_gdiff( cv_nonods * nphase ), ScalarField_Source_Store( cv_nonods * nphase ), &
           ScalarField_Source_Component( cv_nonods * nphase ) )

      sum_theta_flux = 1. ; sum_one_m_theta_flux = 0.
      ScalarField_Source_Store=0. ; ScalarField_Source_Component=0.


!!$ Defining discretisation options
      call Get_Discretisation_Options( state, is_overlapping, &
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
!!$
      have_temperature_field = .false. ; have_component_field = .false. ; have_extra_DiffusionLikeTerm = .false.
      do istate = 1, nstate
         if( have_option( '/material_phase[' // int2str( istate - 1 ) // ']/scalar_field::Temperature' ) ) &
              have_temperature_field = .true.
         if( have_option( '/material_phase[' // int2str( istate - 1 ) // ']/is_multiphase_component' ) ) &
              have_component_field = .true.
!!$
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

!!$ Starting Time Loop 
      itime = 0
      Loop_Time: do
!!$
         itime = itime + 1
         acctim = acctim + dt
         call set_option( '/timestepping/current_time', acctim )

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

!!$ Update all fields from time-step 'N - 1'
         Velocity_U_Old = Velocity_U ; Velocity_V_Old = Velocity_V ; Velocity_W_Old = Velocity_W
         Velocity_NU = Velocity_U ; Velocity_NV = Velocity_V ; Velocity_NW = Velocity_W
         Velocity_NU_Old = Velocity_U ; Velocity_NV_Old = Velocity_V ; Velocity_NW_Old = Velocity_W
         Density_Old = Density ; Pressure_FEM_Old = Pressure_FEM ; Pressure_CV_Old = Pressure_CV
         PhaseVolumeFraction_Old = PhaseVolumeFraction ; Temperature_Old = Temperature ; Component_Old = Component
         Density_Old_tmp = Density_tmp

!!$ Update state memory
!!$         if ( have_temperature_field ) then
!!$            do iphase = 1, nphase
!!$               Temperature_State => extract_scalar_field( state( iphase ), 'Temperature' )
!!$               Temperature_State % val = Temperature( 1 + ( iphase - 1 ) * cv_nonods : iphase * cv_nonods )
!!$            end do
!!$         end if
!!$         Pressure_State => extract_scalar_field( state( 1 ), 'Pressure' )
!!$         Pressure_State % val = Pressure_CV
!!$         do icomp = 1, ncomp
!!$            do iphase = 1, nphase
!!$               Component_State => extract_scalar_field( state( icomp + nphase ), & 
!!$                    'ComponentMassFractionPhase' // int2str( iphase ) )
!!$               Component_State % val = component( 1 + ( iphase - 1 ) * cv_nonods + ( icomp - 1 ) * &
!!$                    nphase * cv_nonods : iphase * cv_nonods + ( icomp - 1 ) * nphase * cv_nonods )
!!$            end do
!!$         end do

         ! evaluate prescribed fields at time = current_time+dt
         call set_prescribed_field_values( state, exclude_interpolated = .true., &
              exclude_nonreprescribed = .true., time = acctim )

         IF( .false. ) THEN
            ! Make sure the FEM representation sums to unity so we don't get surprising results...
            DO CV_INOD = 1, CV_NONODS
               RSUM = 0.0
               DO IPHASE = 1, NCOMP
                  RSUM = RSUM + COMPONENT( CV_INOD + ( IPHASE - 1 ) * CV_NONODS )
               END DO
               DO IPHASE = 1, NPHASE
                  COMPONENT( CV_INOD + ( IPHASE - 1 ) * CV_NONODS ) = COMPONENT( CV_INOD + ( IPHASE - 1 ) &
                       * CV_NONODS ) / RSUM
               END DO
            END DO



if ( .false. .and. itime==1 ) then
               allocate( DEN_CV_NOD( CV_NLOC, NPHASE) ) 
               if ( cv_nloc==6 .or. cv_nloc==10 ) then ! P2 triangle or tet

                  density_tmp = component

                  DO ELE = 1, TOTELE
                     DO CV_ILOC = 1, CV_NLOC
                        CV_NOD = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_ILOC )
                        DO IPHASE = 1,NPHASE
                           CV_NOD_PHA = CV_NOD +( IPHASE - 1) * CV_NONODS
                           DEN_CV_NOD( CV_ILOC, IPHASE ) = density_tmp( CV_NOD_PHA )
                        END DO
                     END DO

                     DEN_CV_NOD(2, :) = 0.5 * ( DEN_CV_NOD(1, :) + DEN_CV_NOD(3, :) )
                     DEN_CV_NOD(4, :) = 0.5 * ( DEN_CV_NOD(1, :) + DEN_CV_NOD(6, :) )
                     DEN_CV_NOD(5, :) = 0.5 * ( DEN_CV_NOD(3, :) + DEN_CV_NOD(6, :) )

                     if ( cv_nloc==10 ) then
                        DEN_CV_NOD(7, :) = 0.5 * ( DEN_CV_NOD(1, :) + DEN_CV_NOD(10, :) )
                        DEN_CV_NOD(8, :) = 0.5 * ( DEN_CV_NOD(3, :) + DEN_CV_NOD(10, :) )
                        DEN_CV_NOD(9, :) = 0.5 * ( DEN_CV_NOD(6, :) + DEN_CV_NOD(10, :) )
                     end if

                     DO CV_ILOC = 1, CV_NLOC
                        CV_NOD = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_ILOC )
                        DO IPHASE = 1, NPHASE
                           CV_NOD_PHA = CV_NOD +( IPHASE - 1) * CV_NONODS
                           component( CV_NOD_PHA ) = DEN_CV_NOD( CV_ILOC, IPHASE )
                        END DO
                     END DO
                  END DO
          
               end if
               deallocate( DEN_CV_NOD ) 

               component_old = component
end if




            do icomp = 1, ncomp
               do iphase = 1, nphase
                  Component_State => extract_scalar_field( state( icomp + nphase ), & 
                       'ComponentMassFractionPhase' // int2str( iphase ) )
                  Component_State % val = component( 1 + ( iphase - 1 ) * cv_nonods + ( icomp - 1 ) * &
                       nphase * cv_nonods : iphase * cv_nonods + ( icomp - 1 ) * nphase * cv_nonods )


                  Component_State => extract_scalar_field( state( icomp + nphase ), & 
                       'ComponentMassFractionPhase' // int2str( iphase ) // 'Old' )
                  Component_State % val = component_old( 1 + ( iphase - 1 ) * cv_nonods + ( icomp - 1 ) * &
                       nphase * cv_nonods : iphase * cv_nonods + ( icomp - 1 ) * nphase * cv_nonods )

               end do
            end do
         END IF





         Loop_NonLinearIteration: do its = 1, NonLinearIteration

            call Calculate_Phase_Component_Densities( state, &
                 Density, DRhoDPressure )

            if( its == 1 ) Density_Old = Density


!!$ Solve advection of the scalar 'Temperature':
            Conditional_ScalarAdvectionField: if( have_temperature_field .and. &
                 have_option( '/material_phase[0]/scalar_field::Temperature/prognostic' ) ) then

               ewrite(3,*)'Now advecting Temperature Field'

               Velocity_NU = Velocity_U ; Velocity_NV = Velocity_V ; Velocity_NW = Velocity_W


               !call calculate_diffusivity( state, ncomp, nphase, ndim, cv_nonods, mat_nonods, ScalarAdvectionField_Diffusion  )

               call INTENERGE_ASSEM_SOLVE( state, &
                    NCOLACV, FINACV, COLACV, MIDACV, & 
                    NCOLCT, FINDCT, COLCT, &
                    CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
                    U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE,  &
                    NPHASE,  &
                    CV_NLOC, U_NLOC, X_NLOC, &
                    CV_NDGLN, X_NDGLN, U_NDGLN, &
                    CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
                    X, Y, Z, &
!!$
                    Velocity_NU, Velocity_NV, Velocity_NW, Velocity_NU_Old, Velocity_NV_Old, Velocity_NW_Old, &
                    ug, vg, wg, &
                    Temperature, Temperature_Old, &
                    Density, Density_Old, &
!!$
                    MAT_NLOC, MAT_NDGLN, MAT_NONODS, ScalarAdvectionField_Diffusion, &
                    t_disopt, t_dg_vel_int_opt, dt, t_theta, t_beta, &
                    Temperature_BC, Density_BC, Velocity_U_BC, Velocity_V_BC, Velocity_W_BC, &
                    suf_t_bc_rob1, suf_t_bc_rob2, &
                    Temperature_BC_Spatial, Density_BC_Spatial, Velocity_U_BC_Spatial, &
                    DRhoDPressure, Pressure_CV, &
                    Temperature_Source, Temperature_Absorption, Porosity, &
                    ndim, &
!!$
                    NCOLM, FINDM, COLM, MIDM, &
!!$
                    XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
!!$
                    opt_vel_upwind_coefs, nopt_vel_upwind_coefs, &
                    Temperature_FEMT, Dummy_PhaseVolumeFraction_FEMT, &
                    igot_t2, Temperature, Temperature_Old, igot_theta_flux, scvngi_theta, &
                    t_get_theta_flux, t_use_theta_flux, &
                    Temperature, Temperature, Temperature, &
                    Temperature_BC, suf_t_bc_rob1, suf_t_bc_rob2, Temperature_BC_Spatial, &
                    in_ele_upwind, dg_ele_upwind, &
!!$                    
                    NOIT_DIM, & ! This need to be removed as it is already deprecated
!!$
!!$                 nits_flux_lim_t, &
                    Mean_Pore_CV, &
                    option_path = '/material_phase[0]/scalar_field::Temperature', &
                    mass_ele_transp = dummy_ele, &
                    thermal = .true. )

!!$ Update state memory
               do iphase = 1, nphase
                  Temperature_State => extract_scalar_field( state( iphase ), 'Temperature' )
                  Temperature_State % val = Temperature( 1 + ( iphase - 1 ) * cv_nonods : iphase * cv_nonods )
               end do

               call Calculate_Phase_Component_Densities( state, &
                    Density, DRhoDPressure )

            end if Conditional_ScalarAdvectionField

            if( solve_force_balance ) then
               call Calculate_AbsorptionTerm( state, &
                    cv_ndgln, mat_ndgln, &
                    PhaseVolumeFraction, Permeability, &
                    nopt_vel_upwind_coefs, opt_vel_upwind_coefs, Material_Absorption )
            end if

!!$ Diffusion-like term -- here used as part of the capillary pressure for porous media. It can also be 
!!$ extended to surface tension -like term.
            iplike_grad_sou = 0
            if( have_option( '/material_phase[0]/multiphase_properties/capillary_pressure' ) )then
               iplike_grad_sou = 1
               call calculate_capillary_pressure( state, cv_nonods, nphase, plike_grad_sou_grad, &
                    PhaseVolumeFraction )
            end if

            CALL CALCULATE_SURFACE_TENSION( state, nphase, ncomp, &
                 PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, IPLIKE_GRAD_SOU, &
                 Velocity_U_Source_CV, Velocity_U_Source, Component, &
                 NCOLACV, FINACV, COLACV, MIDACV, &
                 NCOLCT, FINDCT, COLCT, &
                 CV_NONODS, U_NONODS, X_NONODS, TOTELE, STOTEL, &
                 CV_ELE_TYPE, CV_SELE_TYPE, U_ELE_TYPE, &
                 CV_NLOC, U_NLOC, X_NLOC, CV_SNLOC, U_SNLOC, &
                 CV_NDGLN, CV_SNDGLN, X_NDGLN, U_NDGLN, U_SNDGLN, &
                 X, Y, Z, &
                 MAT_NLOC, MAT_NDGLN, MAT_NONODS,  &
                 NDIM,  &
                 NCOLM, FINDM, COLM, MIDM, &
                 XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
                 Component_BC_Spatial, Component_BC )

            ScalarField_Source_Store = ScalarField_Source + ScalarField_Source_Component
            volfra_use_theta_flux = .true.
            if( ncomp <= 1 ) volfra_use_theta_flux = .false.

            U_Density = 0. ; U_Density_Old = 0.
            if( .not. have_option( '/material_phase[0]/multiphase_properties/relperm_type' ) )then
               U_Density = Density ; U_Density_Old = Density_Old
            end if

!!$ Now solving the Momentum Equation ( = Force Balance Equation )
            Conditional_ForceBalanceEquation: if ( solve_force_balance ) then

!!$ Updating velocities:
               Velocity_NU = Velocity_U ; Velocity_NV = Velocity_V ; Velocity_NW = Velocity_W 
               Velocity_NU_Old = Velocity_U_Old ; Velocity_NV_Old = Velocity_V_Old ; Velocity_NW_Old = Velocity_W_Old

!!$ This calculates u_source_cv = ScalarField_Source_CV -- ie, the buoyancy term and as the name
!!$ suggests it's a CV source term for the velocity field

               !call calculate_u_source_cv( state, cv_nonods, ndim, nphase, Density, Velocity_U_Source_CV )

               density_tmp = density

               ! make mid side nodes the average of the 2 corner nodes...

               allocate( DEN_CV_NOD( CV_NLOC, NPHASE) ) 
               if ( cv_nloc==6 .or. cv_nloc==10 ) then ! P2 triangle or tet

                  DO ELE = 1, TOTELE
                     DO CV_ILOC = 1, CV_NLOC
                        CV_NOD = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_ILOC )
                        DO IPHASE = 1,NPHASE
                           CV_NOD_PHA = CV_NOD +( IPHASE - 1) * CV_NONODS
                           DEN_CV_NOD( CV_ILOC, IPHASE ) = density_tmp( CV_NOD_PHA )
                        END DO
                     END DO

                     DEN_CV_NOD(2, :) = 0.5 * ( DEN_CV_NOD(1, :) + DEN_CV_NOD(3, :) )
                     DEN_CV_NOD(4, :) = 0.5 * ( DEN_CV_NOD(1, :) + DEN_CV_NOD(6, :) )
                     DEN_CV_NOD(5, :) = 0.5 * ( DEN_CV_NOD(3, :) + DEN_CV_NOD(6, :) )

                     if ( cv_nloc==10 ) then
                        DEN_CV_NOD(7, :) = 0.5 * ( DEN_CV_NOD(1, :) + DEN_CV_NOD(10, :) )
                        DEN_CV_NOD(8, :) = 0.5 * ( DEN_CV_NOD(3, :) + DEN_CV_NOD(10, :) )
                        DEN_CV_NOD(9, :) = 0.5 * ( DEN_CV_NOD(6, :) + DEN_CV_NOD(10, :) )
                     end if

                     DO CV_ILOC = 1, CV_NLOC
                        CV_NOD = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_ILOC )
                        DO IPHASE = 1, NPHASE
                           CV_NOD_PHA = CV_NOD +( IPHASE - 1) * CV_NONODS
                           Density_tmp( CV_NOD_PHA ) = DEN_CV_NOD( CV_ILOC, IPHASE )
                        END DO
                     END DO
                  END DO
          
               end if
               deallocate( DEN_CV_NOD ) 

               !if ( (cv_nloc==6 .or. cv_nloc==10) .and. &
               !     .not. have_option( '/material_phase[0]/multiphase_properties/relperm_type' ) &
               !     ) then
               !   U_Density = Density_tmp
               !   U_Density_Old = Density_Old_tmp
               !   if ( its == 1 ) U_Density_Old = Density_tmp
               !end if

               call calculate_u_source_cv( state, cv_nonods, ndim, nphase, density_tmp, Velocity_U_Source_CV )

!!$ Calculate diffusion
               if ( have_option( '/physical_parameters/mobility' ) ) then

                  ! if solving for porous media and mobility is calculated
                  ! through the viscosity ratio this code will fail
                  momentum_diffusion=0.

               else

                  momentum_diffusion=0.

                  t_field => extract_tensor_field( state( 1 ), 'Viscosity', stat )
                  if (stat == 0) then

                     if ( ncomp>1 ) then

!!$                        allocate( momentum_diffusion_tmp( mat_nonods, ndim, ndim ) ) ; momentum_diffusion_tmp = 0.
!!$                        allocate( constant( 1, 1 ) ) ; constant = 0.
!!$
!!$                           do icomp = 1, ncomp
!!$
!!$                                    do iphase = 1, nphase
!!$
!!$
!!$                              option_path = '/material_phase[' // int2str( nphase + icomp - 1 ) // &
!!$                                   ']/vector_field::Velocity/prognostic/tensor_field::ComponentViscosityPhase'// int2str( iphase )
!!$                              option_path=option_path // '/prescribed/value[0]/isotropic/constant'
!!$
!!$         
!!$               call get_option( trim( option_path ), constant( 1, 1 ) )
!!$               do idim = 1, ndim
!!$                 momentum_diffusion_tmp( : , idim, idim ) = constant( 1, 1 )
!!$               end do
!!$               
!!$
!!$
!!$                              momentum_diffusion(:, :, :, iphase) = momentum_diffusion(:, :, :, iphase) + &
!!$                                   momentum_diffusion_tmp * component( 1 + ( iphase - 1 ) * cv_nonods + ( icomp - 1 ) * &
!!$                                   nphase * cv_nonods : iphase * cv_nonods + ( icomp - 1 ) * nphase * cv_nonods
!!$
!!$                           end do
!!$                        end do
!!$
!!$                        deallocate( constant )
!!$                        deallocate( momentum_diffusion_tmp )

                        momentum_diffusion=0.
                        momentum_diffusion(:, 1,1,1) = 0. !1.e-5
                        momentum_diffusion(:, 2,2,1) = 0. !1.e-5

                     else

                        do iphase = 1, nphase

                           t_field => extract_tensor_field( state( iphase ), 'Viscosity', stat )

                           option_path = '/material_phase[' // int2str( iphase - 1 ) // &
                                ']/vector_field::Velocity/prognostic/tensor_field::Viscosity'
                           call Extract_TensorFields_Outof_State( state, iphase, &
                                t_field, option_path, &
                                momentum_diffusion( :, :, :, iphase ), &
                                mat_ndgln )
                        end do

                     end if

                  end if

               end if

               CALL FORCE_BAL_CTY_ASSEM_SOLVE( state, &
                    NDIM, NPHASE, U_NLOC, X_NLOC, P_NLOC, CV_NLOC, MAT_NLOC, TOTELE, &
                    U_ELE_TYPE, P_ELE_TYPE, &
                    U_NONODS, CV_NONODS, X_NONODS, MAT_NONODS, &
                    U_NDGLN, P_NDGLN, CV_NDGLN, X_NDGLN, MAT_NDGLN,&
                    STOTEL, CV_SNDGLN, U_SNDGLN, P_SNDGLN, &
                    U_SNLOC, P_SNLOC, CV_SNLOC, &
!!$
                    x, y, z, Material_Absorption_Stab, Material_Absorption, Velocity_U_Source, Velocity_U_Source_CV, &
                    Velocity_U, Velocity_V, Velocity_W, Velocity_U_Old, Velocity_V_Old, Velocity_W_Old, &
                    Pressure_FEM, Pressure_CV, Density, Density_Old, PhaseVolumeFraction, PhaseVolumeFraction_Old, & 
                    DRhoDPressure, &
                    dt, &
!!$
                    NCOLC, FINDC, COLC, & ! C sparsity - global cty eqn 
                    NCOLDGM_PHA, FINDGM_PHA, COLDGM_PHA, MIDDGM_PHA, &! Force balance sparsity
                    NCOLELE, FINELE, COLELE, & ! Element connectivity.
                    NCOLCMC, FINDCMC, COLCMC, MIDCMC, & ! pressure matrix for projection method
                    NCOLACV, FINACV, COLACV, MIDACV, & ! For CV discretisation method
                    NLENMCY, NCOLMCY, FINMCY, COLMCY, MIDMCY, & ! Force balance plus cty multi-phase eqns
                    NCOLCT, FINDCT, COLCT, & ! CT sparsity - global cty eqn.
                    CV_ELE_TYPE, &
!!$
                    Velocity_NU, Velocity_NV, Velocity_NW, Velocity_NU_Old, Velocity_NV_Old, Velocity_NW_Old, &
                    v_disopt, v_dg_vel_int_opt, v_theta, &
                    PhaseVolumeFraction_BC, Density_BC, Velocity_U_BC, Velocity_V_BC, Velocity_W_BC, Pressure_FEM_BC, &
                    suf_u_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob1, suf_v_bc_rob2, &
                    suf_w_bc_rob1, suf_w_bc_rob2, &
                    PhaseVolumeFraction_BC_Spatial, Density_BC_Spatial, Velocity_U_BC_Spatial, Pressure_FEM_BC_Spatial, &
                    ScalarField_Source_Store, ScalarField_Absorption, Porosity, &
!!$
                    NCOLM, FINDM, COLM, MIDM, & ! Sparsity for the CV-FEM
                    XU_NLOC, XU_NDGLN, &
!!$
                    U_Density, U_Density_Old, Momentum_Diffusion, &
                    opt_vel_upwind_coefs, nopt_vel_upwind_coefs, &
                    igot_theta_flux, scvngi_theta, volfra_use_theta_flux, &
                    sum_theta_flux, sum_one_m_theta_flux, &
                    in_ele_upwind, dg_ele_upwind, &
!!$                    
                    NOIT_DIM, & ! This need to be removed as it is already deprecated
!!$
                    iplike_grad_sou, plike_grad_sou_coef, plike_grad_sou_grad, &
                    scale_momentum_by_volume_fraction )

               Pressure_State => extract_scalar_field( state( 1 ), 'Pressure' )
               Pressure_State % val = Pressure_CV

            end if Conditional_ForceBalanceEquation

            Conditional_PhaseVolumeFraction: if ( solve_PhaseVolumeFraction ) then

               call VolumeFraction_Assemble_Solve( state, &
                    NCOLACV, FINACV, COLACV, MIDACV, &
                    NCOLCT, FINDCT, COLCT, &
                    CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
                    CV_ELE_TYPE,  &
                    NPHASE,  &
                    CV_NLOC, U_NLOC, X_NLOC,  &
                    CV_NDGLN, X_NDGLN, U_NDGLN, &
                    CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
!!$
                    x, y, z, Velocity_U, Velocity_V, Velocity_W, &
                    Velocity_NU, Velocity_NV, Velocity_NW, Velocity_NU_Old, Velocity_NV_Old, Velocity_NW_Old, &
                    PhaseVolumeFraction, PhaseVolumeFraction_Old, &
                    Density, Density_Old, &
!!$
                    MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
!!$
                    v_disopt, v_dg_vel_int_opt, dt, v_theta, v_beta, &
                    PhaseVolumeFraction_BC, Density_BC, Velocity_U_BC, Velocity_V_BC, Velocity_W_BC, &
                    PhaseVolumeFraction_BC_Spatial, Density_BC_Spatial, Velocity_U_BC_Spatial, &
                    DRhoDPressure, Pressure_FEM, &
                    ScalarField_Source_Store, ScalarField_Absorption, Porosity, &
!!$
                    NDIM, &
                    NCOLM, FINDM, COLM, MIDM, &
                    XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
!!$
                    opt_vel_upwind_coefs, nopt_vel_upwind_coefs, &
                    PhaseVolumeFraction_FEMT, Density_FEMT, &
                    igot_theta_flux, scvngi_theta, volfra_use_theta_flux, &
                    sum_theta_flux, sum_one_m_theta_flux, &
                    in_ele_upwind, dg_ele_upwind, &
!!$                    
                    NOIT_DIM, & ! This need to be removed as it is already deprecated
!!$
!!$                 nits_flux_lim_volfra, &
                    option_path = '/material_phase[0]/scalar_field::PhaseVolumeFraction', &
                    mass_ele_transp = mass_ele )

            end if Conditional_PhaseVolumeFraction

!!$ Starting loop over components
            sum_theta_flux = 0. ; sum_one_m_theta_flux = 0. ; ScalarField_Source_Component = 0.

            Conditional_Components: if( have_component_field ) then
               Loop_Components: do icomp = 1, ncomp

                  ! NEEDS GENERALIZING *********************************
                  IF(ICOMP==1) THEN
                     DENSITY_COMPONENT(( ICOMP - 1 ) * NPHASE * CV_NONODS + 1 : ICOMP * NPHASE * CV_NONODS )=1000.0
                     DENSITY_COMPONENT_OLD(( ICOMP - 1 ) * NPHASE * CV_NONODS + 1 : ICOMP * NPHASE * CV_NONODS )=1000.0
                  ELSE
                     DENSITY_COMPONENT(( ICOMP - 1 ) * NPHASE * CV_NONODS + 1 : ICOMP * NPHASE * CV_NONODS )=1.0
                     DENSITY_COMPONENT_OLD(( ICOMP - 1 ) * NPHASE * CV_NONODS + 1 : ICOMP * NPHASE * CV_NONODS )=1.0
                  ENDIF
                  ! NEEDS GENERALIZING *********************************



!!$ generalisation for above... i.e. calc DENSITY_COMPONENT
!!$
!!$           allocate( Density_tmp( cv_nonods ), DRhoDPressure_tmp( cv_nonods ) )
!!$           Density_tmp=0. ; DRhoDPressure_tmp=0.
!!$
!!$
!!$      do iphase = 1, nphase
!!$
!!$         eos_option_path = trim( '/material_phase[' // int2str( icomp + nphase - 1 ) // ']/equation_of_state' )
!!$
!!$              call Assign_Equation_of_State( eos_option_path )
!!$
!!$
!!$                  call Computing_Perturbation_Density( state, &
!!$                       iphase, icomp, icomp, eos_option_path, &
!!$                       Density_tmp, &
!!$                       DRhoDPressure_tmp, .true. )
!!$
!!$                Density_Component( ( icomp - 1 ) * nphase * cv_nonods + 1 : icomp * nphase * cv_nonods ) = Density_tmp
!!$                Density_Component_Old( ( icomp - 1 ) * nphase * cv_nonods + 1 : icomp * nphase * cv_nonods ) = Density_tmp
!!$         
!!$         ( icomp - 1 ) * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods + 1 :    icomp * nphase * cv_nonods + ( iphase - 1 ) * cv_nonods
!!$
!!$
!!$     end do
!!$    deallocate( Density_tmp, DRhoDPressure_tmp )
!!$


!!$ Computing the absorption term for the multi-components equation
                  call Calculate_ComponentAbsorptionTerm( state, &
                       icomp, cv_ndgln, & 
                       Density_Old, PhaseVolumeFraction_Old, Porosity, mass_ele, &
                       Component_Absorption )

                  Conditional_SmoothAbsorption: if( have_option( '/material_phase[' // int2str( nstate - ncomp ) // &
                       ']/is_multiphase_component/KComp_Sigmoid' ) .and. nphase > 1 ) then
                     do cv_nodi = 1, cv_nonods
                        if( PhaseVolumeFraction_Old( cv_nodi ) > 0.95 ) then
                           do iphase = 1, nphase
                              do jphase = min( iphase + 1, nphase ), nphase
                                 Component_Absorption( cv_nodi, iphase, jphase ) = &
                                      Component_Absorption( cv_nodi, iphase, jphase ) * max( 0.01, &
                                      20. * PhaseVolumeFraction_Old( cv_nodi ) )
                              end do
                           end do
                        end if
                     end do
                  end if Conditional_SmoothAbsorption

!!$ Computing diffusion term for the component conservative equation:
                  call Calculate_ComponentDiffusionTerm( state, &
                       mat_ndgln, u_ndgln, x_ndgln, &
                       x, y, z, Velocity_NU, Velocity_NV, Velocity_NW, &
                       u_ele_type, p_ele_type, ncomp_diff_coef, comp_diffusion_opt, &
                       Component_Diffusion_Operator_Coefficient( icomp, :, : ), &
                       Component_Diffusion )

!!$ NonLinear iteration for the components advection:
                  Loop_NonLinearIteration_Components: do its2 = 1, NonLinearIteration_Components
                     comp_use_theta_flux = .false. ; comp_get_theta_flux = .true.

                     call INTENERGE_ASSEM_SOLVE( state, &
                          NCOLACV, FINACV, COLACV, MIDACV, & ! CV sparsity pattern matrix
                          NCOLCT, FINDCT, COLCT, &
                          CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
                          U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE,  &
                          NPHASE,  &
                          CV_NLOC, U_NLOC, X_NLOC,  &
                          CV_NDGLN, X_NDGLN, U_NDGLN, &
                          CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
                          X, Y, Z, &
!!$
                          Velocity_NU, Velocity_NV, Velocity_NW, Velocity_NU_Old, Velocity_NV_Old, Velocity_NW_Old, &
                          ug, vg, wg, &
                          Component( ( icomp - 1 ) * nphase * cv_nonods + 1 : icomp * nphase * cv_nonods ), &
                          Component_Old( ( icomp - 1 ) * nphase * cv_nonods + 1 : icomp * nphase * cv_nonods ), &
                          DENSITY_COMPONENT(( ICOMP - 1 ) * NPHASE * CV_NONODS + 1 : ICOMP * NPHASE * CV_NONODS ), &
                          DENSITY_COMPONENT_OLD(( ICOMP - 1 ) * NPHASE * CV_NONODS + 1 : ICOMP * NPHASE * CV_NONODS ), &
!!$
                          MAT_NLOC, MAT_NDGLN, MAT_NONODS, Component_Diffusion, &
                          v_disopt, v_dg_vel_int_opt, dt, v_theta, v_beta, &
                          Component_BC( 1 + stotel * cv_snloc * nphase * ( icomp - 1 ) : stotel * cv_snloc * nphase * icomp ), &
                          Density_BC, Velocity_U_BC, Velocity_V_BC, Velocity_W_BC, &
                          suf_comp_bc_rob1, suf_comp_bc_rob2, &
                          Component_BC_Spatial, Density_BC_Spatial, Velocity_U_BC_Spatial, &
                          DRhoDPressure, Pressure_FEM, &
                          Component_Source, Component_Absorption, Porosity, &
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
                          theta_flux, one_m_theta_flux, theta_gdiff, &
                          PhaseVolumeFraction_BC, suf_vol_bc_rob1, suf_vol_bc_rob2, PhaseVolumeFraction_BC_Spatial, &
                          in_ele_upwind, dg_ele_upwind, &
!!$
                          NOIT_DIM, & ! This need to be removed as it is already deprecated
!!$
!!$                          nits_flux_lim_comp, &
                          Mean_Pore_CV, &
                                !option_path = '', &
                          mass_ele_transp = dummy_ele, &
                          thermal = .false. ) ! the false means that we don't add an extra source term

                  end do Loop_NonLinearIteration_Components

                  sum_theta_flux = sum_theta_flux + theta_flux
                  sum_one_m_theta_flux = sum_one_m_theta_flux + one_m_theta_flux

!!$
!!$ And Here, zeroed Component_Absorption within ScalarField_Source_Component
!!$
!!!!! We have divided through by density 
                  ScalarField_Source_Component = ScalarField_Source_Component + THETA_GDIFF

                  Loop_Phase_SourceTerm1: do iphase = 1, nphase
                     Loop_Phase_SourceTerm2: do jphase = 1, nphase
                        DO CV_NODI = 1, CV_NONODS
                           ScalarField_Source_Component( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) = &
                                ScalarField_Source_Component( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) - &
                                Component_Absorption( CV_NODI, IPHASE, JPHASE ) * &
                                Component( CV_NODI + ( JPHASE - 1 ) * CV_NONODS + &
                                ( ICOMP - 1 ) * NPHASE * CV_NONODS )  &
                                / DENSITY_COMPONENT( CV_NODI + ( IPHASE - 1 ) * CV_NONODS + &
                                ( ICOMP - 1 ) * NPHASE * CV_NONODS )
                        END DO
                     end do Loop_Phase_SourceTerm2
                  end do Loop_Phase_SourceTerm1

                  ! For compressability
                  DO IPHASE = 1, NPHASE
                     DO CV_NODI = 1, CV_NONODS
                        ScalarField_Source_Component( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) = &
                             ScalarField_Source_Component( CV_NODI + ( IPHASE - 1 ) * CV_NONODS )  &
                             + Mean_Pore_CV(CV_NODI)* Component( CV_NODI + ( IPHASE - 1 ) * CV_NONODS + &
                             ( ICOMP - 1 ) * NPHASE * CV_NONODS ) &
                             * ( DENSITY_COMPONENT_OLD( CV_NODI + ( IPHASE - 1 ) * CV_NONODS + &
                             ( ICOMP - 1 ) * NPHASE * CV_NONODS ) &
                             -DENSITY_COMPONENT( CV_NODI + ( IPHASE - 1 ) * CV_NONODS + &
                             ( ICOMP - 1 ) * NPHASE * CV_NONODS )   ) &
                             *PhaseVolumeFraction_Old(CV_NODI + ( IPHASE - 1 ) * CV_NONODS) &
                             / DENSITY_COMPONENT( CV_NODI + ( IPHASE - 1 ) * CV_NONODS + &
                             ( ICOMP - 1 ) * NPHASE * CV_NONODS )
                     END DO
                  END DO

               end do Loop_Components

               if( have_option( '/material_phase[' // int2str( nstate - ncomp ) // & 
                    ']/is_multiphase_component/Comp_Sum2One' ) .and. ( ncomp > 1 ) ) then
                  call Cal_Comp_Sum2One_Sou( ScalarField_Source_Component, cv_nonods, nphase, ncomp, dt, its, &
                       NonLinearIteration, &
                       Porosity, PhaseVolumeFraction, PhaseVolumeFraction_Old, Density_Component, Density_Component_Old, &
                       Component, Component_Old )
               end if

               ! Update state memory
               do icomp = 1, ncomp
                  do iphase = 1, nphase
                     Component_State => extract_scalar_field( state( icomp + nphase ), & 
                          'ComponentMassFractionPhase' // int2str( iphase ) )
                     Component_State % val = component( 1 + ( iphase - 1 ) * cv_nonods + ( icomp - 1 ) * &
                          nphase * cv_nonods : iphase * cv_nonods + ( icomp - 1 ) * nphase * cv_nonods )
                  end do
               end do

            end if Conditional_Components

         end do Loop_NonLinearIteration

         call set_option( '/timestepping/current_time', acctim )
         call set_option( '/timestepping/timestep', dt)

!!$ Copying fields back to state:
         call copy_into_state( state, & ! Copying main fields into state
              PhaseVolumeFraction, Temperature, Pressure_CV, Velocity_U, Velocity_V, Velocity_W, &
!!$              PhaseVolumeFraction, Temperature, Pressure_FEM, Velocity_U, Velocity_V, Velocity_W, &
              Density, Component, &
              ncomp, nphase, cv_ndgln, p_ndgln, u_ndgln, ndim )

         Conditional_TimeDump: if( ( mod( itime, dump_period_in_timesteps ) == 0 ) .or. ( itime == 1 ) ) then
            call get_option( '/timestepping/current_time', current_time ) ! Find the current time 
!!$ Calculate diagnostic fields and write into the stat file
            call calculate_diagnostic_variables( state, exclude_nonrecalculated = .true. )
            call calculate_diagnostic_variables_new( state, exclude_nonrecalculated = .true. )
            call write_diagnostics( state, current_time, dt, itime )
            not_to_move_det_yet = .false. ; dump_no = itime ! Sync dump_no with itime
            call write_state( dump_no, state ) ! Now writing into the vtu files
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

         Conditional_ReallocatingFields: if( do_reallocate_fields ) then

            Conditional_Adaptivity: if( have_option( '/mesh_adaptivity/hr_adaptivity ') ) then

               Conditional_Adapt_by_TimeStep: if( mod( itime, adapt_time_steps ) == 0 ) then
!!$               Conditional_Adapt_by_TimeStep: if( do_adapt_mesh( current_time, itime ) ) then
!!$               Conditional_Adapt_by_TimeStep: if( do_adapt_mesh( current_time, timestep ) ) then

                  call pre_adapt_tasks( sub_state )

                  call qmesh( state, metric_tensor )

                  if( have_option( '/io/stat/output_before_adapts' ) ) call write_diagnostics( state, current_time, dt, &
                       itime, not_to_move_det_yet = .true. )
!!$                  if( have_option( '/io/stat/output_before_adapts' ) ) call write_diagnostics( state, current_time, dt, &
!!$                       timestep, not_to_move_det_yet = .true. )

                  call run_diagnostics( state )

                  call adapt_state( state, metric_tensor )

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

!!$ Deallocating array variables:
            deallocate( &
!!$ Node glabal numbers
                 x_ndgln_p1, x_ndgln, cv_ndgln, p_ndgln, &
                 mat_ndgln, u_ndgln, xu_ndgln, cv_sndgln, p_sndgln, u_sndgln, &
!!$ Sparsity patterns
                 finacv, colacv, midacv, finmcy, colmcy, midmcy, &
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
                 PhaseVolumeFraction_BC_Spatial, Pressure_FEM_BC_Spatial, &
                 Density_BC_Spatial, Component_BC_Spatial, Velocity_U_BC_Spatial, Temperature_BC_Spatial, &
                 xu, yu, zu, x, y, z, ug, vg, wg, &
                 Velocity_U, Velocity_V, Velocity_W, Velocity_U_Old, Velocity_V_Old, Velocity_W_Old, &
                 Velocity_NU, Velocity_NV, Velocity_NW, Velocity_NU_Old, Velocity_NV_Old, Velocity_NW_Old, &
                 Pressure_FEM, Pressure_CV, Temperature, Density, Density_Component, PhaseVolumeFraction, &
                 Component, U_Density, &
                 Pressure_FEM_Old, Pressure_CV_Old, Temperature_Old, Density_Old, Density_Component_Old, &
                 PhaseVolumeFraction_Old, Component_Old, &
                 U_Density_Old, DRhoDPressure, &
                 Porosity, &
                 Velocity_U_Source, Velocity_U_Source_CV, Temperature_Source, PhaseVolumeFraction_Source, &
                 ScalarField_Source, Component_Source, ScalarAdvectionField_Source, &
                 PhaseVolumeFraction_BC, Pressure_FEM_BC, &
                 Density_BC, Component_BC, Velocity_U_BC, Velocity_V_BC, Velocity_W_BC, Temperature_BC, &
                 suf_u_bc_rob1, suf_v_bc_rob1, suf_w_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob2, suf_w_bc_rob2, &
                 suf_t_bc_rob1, suf_t_bc_rob2, suf_vol_bc_rob1, suf_vol_bc_rob2, suf_comp_bc_rob1, suf_comp_bc_rob2, &
                 theta_gdiff,  ScalarField_Source_Store, ScalarField_Source_Component, &
                 mass_ele, dummy_ele, &
                 Permeability, Material_Absorption, Material_Absorption_Stab, &
                 Velocity_Absorption, ScalarField_Absorption, Component_Absorption, Temperature_Absorption, &
                 Component_Diffusion_Operator_Coefficient, &
                 Momentum_Diffusion, ScalarAdvectionField_Diffusion, &
                 Component_Diffusion, &
                 theta_flux, one_m_theta_flux, sum_theta_flux, sum_one_m_theta_flux )

!!$  Compute primary scalars used in most of the code
            call Get_Primary_Scalars( state, &         
                 nphase, nstate, ncomp, totele, ndim, stotel, &
                 u_nloc, xu_nloc, cv_nloc, x_nloc, x_nloc_p1, p_nloc, mat_nloc, &
                 x_snloc, cv_snloc, u_snloc, p_snloc, &
                 cv_nonods, mat_nonods, u_nonods, xu_nonods, x_nonods, x_nonods_p1, p_nonods, dx, &
                 is_overlapping )
!!$ Calculating Global Node Numbers
            allocate( x_ndgln_p1( totele * x_nloc_p1 ), x_ndgln( totele * x_nloc ), cv_ndgln( totele * cv_nloc ), &
                 p_ndgln( totele * p_nloc ), mat_ndgln( totele * mat_nloc ), u_ndgln( totele * u_nloc ), &
                 xu_ndgln( totele * xu_nloc ), cv_sndgln( stotel * cv_snloc ), p_sndgln( stotel * p_snloc ), &
                 u_sndgln( stotel * u_snloc ) )

            x_ndgln_p1 = 0 ; x_ndgln = 0 ; cv_ndgln = 0 ; p_ndgln = 0 ; mat_ndgln = 0 ; u_ndgln = 0 ; xu_ndgln = 0 ; &
                 cv_sndgln = 0 ; p_sndgln = 0 ; u_sndgln = 0

            call Compute_Node_Global_Numbers( state, &
                 is_overlapping, totele, stotel, x_nloc, x_nloc_p1, cv_nloc, p_nloc, u_nloc, xu_nloc, &
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

!!$ Allocating space for various arrays:
            allocate( xu( xu_nonods ), yu( xu_nonods ), zu( xu_nonods ), &
                 x( x_nonods ), y( x_nonods ), z( x_nonods ), &
                 ug( u_nonods * nphase ), vg( u_nonods * nphase ), wg( u_nonods * nphase ), &
!!$
                 Velocity_U( u_nonods * nphase ), Velocity_V( u_nonods * nphase ), Velocity_W( u_nonods * nphase ), &
                 Velocity_U_Old( u_nonods * nphase ), Velocity_V_Old( u_nonods * nphase ), Velocity_W_Old( u_nonods * nphase ), &
                 Velocity_NU( u_nonods * nphase ), Velocity_NV( u_nonods * nphase ), Velocity_NW( u_nonods * nphase ), &
                 Velocity_NU_Old( u_nonods * nphase ), Velocity_NV_Old( u_nonods * nphase ), Velocity_NW_Old( u_nonods * nphase ), &
!!$
                 Pressure_FEM( cv_nonods ), Pressure_CV( cv_nonods ), &
                 Temperature( nphase * cv_nonods ), Density( nphase * cv_nonods ), &
                 Density_Component(nphase * cv_nonods * ncomp), &
                 PhaseVolumeFraction( nphase * cv_nonods ), Component( nphase * cv_nonods * ncomp ), &
                 U_Density( nphase * cv_nonods ), DRhoDPressure( nphase * cv_nonods ), &
!!$
                 Pressure_FEM_Old( cv_nonods ), Pressure_CV_Old( cv_nonods ), &
                 Temperature_Old( nphase * cv_nonods ), Density_Old( nphase * cv_nonods ), &
                 Density_Component_Old(nphase * cv_nonods * ncomp), &
                 PhaseVolumeFraction_Old( nphase * cv_nonods ), Component_Old( nphase * cv_nonods * ncomp ), &
                 U_Density_Old( nphase * cv_nonods ), &
!!$
                 PhaseVolumeFraction_BC_Spatial( stotel * nphase ), Pressure_FEM_BC_Spatial( stotel * nphase ), &
                 Density_BC_Spatial( stotel * nphase ), Component_BC_Spatial( stotel * nphase ), &
                 Velocity_U_BC_Spatial( stotel * nphase ), Temperature_BC_Spatial( stotel * nphase ), &
                 PhaseVolumeFraction_BC( stotel * cv_snloc * nphase ), Pressure_FEM_BC( stotel * p_snloc * nphase ), &
                 Density_BC( stotel * cv_snloc * nphase ), Temperature_BC( stotel * cv_snloc * nphase ), &
                 Component_BC( stotel * cv_snloc * nphase * ncomp ), &
                 Velocity_U_BC( stotel * u_snloc * nphase ), Velocity_V_BC( stotel * u_snloc * nphase ), &
                 Velocity_W_BC( stotel * u_snloc * nphase ), Temperature_Source( cv_nonods * nphase ), &
                 suf_u_bc_rob1( stotel * u_snloc * nphase ), suf_v_bc_rob1( stotel * u_snloc * nphase ), &
                 suf_w_bc_rob1( stotel * u_snloc * nphase ), suf_u_bc_rob2( stotel * u_snloc * nphase ), &
                 suf_v_bc_rob2( stotel * u_snloc * nphase ), suf_w_bc_rob2( stotel * u_snloc * nphase ), &
                 suf_t_bc_rob1( stotel * cv_snloc * nphase ), suf_t_bc_rob2( stotel * cv_snloc * nphase ), &
                 suf_vol_bc_rob1( stotel * cv_snloc * nphase ), suf_vol_bc_rob2( stotel * cv_snloc * nphase ), &
                 suf_comp_bc_rob1( stotel * cv_snloc * nphase ), suf_comp_bc_rob2( stotel * cv_snloc * nphase ), &
!!$
                 Porosity( totele ), &
                 PhaseVolumeFraction_FEMT( cv_nonods * nphase ), Temperature_FEMT( cv_nonods * nphase ), &
                 Density_FEMT( cv_nonods * nphase ), Component_FEMT( cv_nonods * nphase * ncomp ), &
                 Mean_Pore_CV( cv_nonods ),  SumConc_FEMT( cv_nonods * ncomp ), &
                 Dummy_PhaseVolumeFraction_FEMT( cv_nonods * nphase ), dummy_ele( totele ), mass_ele( totele ), &
!!$
                 PhaseVolumeFraction_Source( cv_nonods * nphase ), Velocity_U_Source( u_nonods * nphase * ndim ), &
                 Velocity_U_Source_CV( cv_nonods * nphase * ndim ), Component_Source( cv_nonods * nphase ), &
                 ScalarField_Source( cv_nonods * nphase ), ScalarAdvectionField_Source( cv_nonods * nphase ), &
!!$
                 Permeability( totele, ndim, ndim ), &
!!$
                 Material_Absorption( mat_nonods, ndim * nphase, ndim * nphase ), &
                 Velocity_Absorption( u_nloc * totele, ndim * nphase, ndim * nphase ), &
                 Material_Absorption_Stab( mat_nonods, ndim * nphase, ndim * nphase ), & 
                 ScalarField_Absorption( cv_nonods, nphase, nphase ), Component_Absorption( cv_nonods, nphase, nphase ), &
                 Temperature_Absorption( cv_nonods, nphase, nphase ), &
                 Momentum_Diffusion( mat_nonods, ndim, ndim, nphase ), &
                 ScalarAdvectionField_Diffusion( mat_nonods, ndim, ndim, nphase ), & 
                 Component_Diffusion( mat_nonods, ndim, ndim, nphase ), &
!!$ Variables used in the diffusion-like term: capilarity and surface tension:
                 plike_grad_sou_grad( cv_nonods * nphase ), &
                 plike_grad_sou_coef( cv_nonods * nphase ) )    
!!$
            Velocity_U=0. ; Velocity_V=0. ; Velocity_W=0.
            Velocity_U_Old=0. ; Velocity_V_Old=0. ; Velocity_W_Old=0.
            Velocity_NU=0. ; Velocity_NV=0. ; Velocity_NW=0.
            Velocity_NU_Old=0. ; Velocity_NV_Old=0. ; Velocity_NW_Old=0.
            UG=0. ; VG=0. ; WG=0.
            Velocity_U_Source = 0. ; Velocity_Absorption = 0. ; Velocity_U_Source_CV = 0. 
            Velocity_U_BC_Spatial=0 ; Velocity_U_BC=0. ; Velocity_V_BC=0. ; Velocity_W_BC=0.
            Momentum_Diffusion=0.
!!$
            Temperature=0. ; Temperature_Source=0. ; Temperature_BC_Spatial=0 ; Temperature_BC=0.
            Temperature_FEMT=0. ; Temperature_Absorption=0.
!!$
            Component=0. ; Component_BC_Spatial=0 ; Component_BC=0. ; Component_Source=0.
            Component_Diffusion=0. ; Component_Absorption=0.
!!$
            Porosity=0. ; Permeability=0.
!!$
            Pressure_CV=0. ; Pressure_FEM=0. ; 
            Pressure_FEM_Old=0. ; Pressure_CV_Old=0.
            Pressure_FEM_BC_Spatial=0 ; Pressure_FEM_BC=0.
!!$
            PhaseVolumeFraction=0. ; PhaseVolumeFraction_Old=0. ; PhaseVolumeFraction_BC_Spatial=0
            PhaseVolumeFraction_BC=0. ; PhaseVolumeFraction_Source=0.
            PhaseVolumeFraction_FEMT=0. ; Dummy_PhaseVolumeFraction_FEMT=0.
!!$
            Density=0. ; Density_Old=0. ; Density_BC_Spatial=0 ; Density_BC=0.
            U_Density=0. ; U_Density_Old=0.
!!$
            ScalarAdvectionField_Diffusion=0. ; ScalarField_Absorption=0.
            ScalarField_Source=0. ; ScalarAdvectionField_Source=0.
!!$
            Material_Absorption=0. ; Material_Absorption_Stab=0.
!!$
            plike_grad_sou_grad=0. ; plike_grad_sou_coef=0.
!!$
            suf_u_bc_rob1=0. ; suf_v_bc_rob1=0. 
            suf_w_bc_rob1=0. ; suf_u_bc_rob2=0. 
            suf_v_bc_rob2=0. ; suf_w_bc_rob2=0. 
            suf_t_bc_rob1=0. ; suf_t_bc_rob2=0. 
            suf_vol_bc_rob1=0. ; suf_comp_bc_rob1=0. ; suf_comp_bc_rob2=0.
!!$

!!$ Extracting Mesh Dependent Fields
            initialised = .true.
            call Extracting_MeshDependentFields_From_State( state, initialised, &
                 xu, yu, zu, x, y, z, &
                 PhaseVolumeFraction, PhaseVolumeFraction_BC_Spatial, PhaseVolumeFraction_BC, PhaseVolumeFraction_Source, &
                 Pressure_CV, Pressure_FEM, Pressure_FEM_BC_Spatial, Pressure_FEM_BC, &
                 Density, Density_BC_Spatial, Density_BC, &
                 Component, Component_BC_Spatial, Component_BC, Component_Source, &
                 Velocity_U, Velocity_V, Velocity_W, Velocity_NU, Velocity_NV, Velocity_NW, &
                 Velocity_U_BC_Spatial, Velocity_U_BC, Velocity_V_BC, Velocity_W_BC, Velocity_U_Source, Velocity_Absorption, &
                 Temperature, Temperature_BC_Spatial, Temperature_BC, Temperature_Source, &
                 Porosity, Permeability )

!!$ Dummy field used in the scalar advection option:
            Dummy_PhaseVolumeFraction_FEMT = 1.

!!$
!!$ Initialising Robin boundary conditions --  this still need to be defined in the schema:
!!$
            suf_u_bc_rob1 = 0. ; suf_v_bc_rob1 = 0. ; suf_w_bc_rob1 = 0. ; suf_u_bc_rob2 = 0. ; suf_v_bc_rob2 = 0. ; &
                 suf_w_bc_rob2 = 0. ; suf_t_bc_rob1 = 0. ; suf_t_bc_rob2 = 0. ; suf_vol_bc_rob1 = 0. ; suf_vol_bc_rob2 = 0. ; &
                 suf_comp_bc_rob1 = 0. ; suf_comp_bc_rob2 = 0.

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

            allocate( theta_flux( totele * igot_theta_flux, cv_nloc, scvngi_theta, nphase ), &
                 one_m_theta_flux( totele * igot_theta_flux, cv_nloc, scvngi_theta, nphase ), &
                 sum_theta_flux( totele * igot_theta_flux, cv_nloc, scvngi_theta, nphase ), &
                 sum_one_m_theta_flux( totele * igot_theta_flux, cv_nloc, scvngi_theta, nphase ), &
                 theta_gdiff( cv_nonods * nphase ), ScalarField_Source_Store( cv_nonods * nphase ), &
                 ScalarField_Source_Component( cv_nonods * nphase ) )

            sum_theta_flux = 1. ; sum_one_m_theta_flux = 0.  
            ScalarField_Source_Store=0. ; ScalarField_Source_Component=0.

            allocate( Component_Diffusion_Operator_Coefficient( ncomp, ncomp_diff_coef, nphase ) )  
            nopt_vel_upwind_coefs = mat_nonods * nphase * ndim * ndim * 2
            allocate( opt_vel_upwind_coefs( nopt_vel_upwind_coefs ) ) ; opt_vel_upwind_coefs = 0.

!!$            !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!$            call copy_into_state( state, &
!!$                 PhaseVolumeFraction, Temperature, Pressure_CV, Velocity_U, Velocity_V, Velocity_W, &
!!$                 Density, Component, &
!!$                 ncomp, nphase, cv_ndgln, p_ndgln, u_ndgln, ndim )
!!$            dump_no=666
!!$            call write_state( dump_no, state ) ! Now writing into the vtu files
!!$            stop 777
!!$            !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

         end if Conditional_ReallocatingFields

      end do Loop_Time

!!$ Now deallocating arrays:
      deallocate( &
!!$ Node glabal numbers
           x_ndgln_p1, x_ndgln, cv_ndgln, p_ndgln, &
           mat_ndgln, u_ndgln, xu_ndgln, cv_sndgln, p_sndgln, u_sndgln, &
!!$ Sparsity patterns
           finacv, colacv, midacv, finmcy, colmcy, midmcy, &
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
           PhaseVolumeFraction_BC_Spatial, Pressure_FEM_BC_Spatial, &
           Density_BC_Spatial, Component_BC_Spatial, Velocity_U_BC_Spatial, Temperature_BC_Spatial, &
           xu, yu, zu, x, y, z, ug, vg, wg, &
           Velocity_U, Velocity_V, Velocity_W, Velocity_U_Old, Velocity_V_Old, Velocity_W_Old, &
           Velocity_NU, Velocity_NV, Velocity_NW, Velocity_NU_Old, Velocity_NV_Old, Velocity_NW_Old, &
           Pressure_FEM, Pressure_CV, Temperature, Density, Density_Component, PhaseVolumeFraction, &
           Component, U_Density, &
           Pressure_FEM_Old, Pressure_CV_Old, Temperature_Old, Density_Old, Density_Component_Old, &
           PhaseVolumeFraction_Old, Component_Old, &
           U_Density_Old, DRhoDPressure, &
           Porosity, &
           Velocity_U_Source, Velocity_U_Source_CV, Temperature_Source, PhaseVolumeFraction_Source, &
           ScalarField_Source, Component_Source, ScalarAdvectionField_Source, &
           PhaseVolumeFraction_BC, Pressure_FEM_BC, &
           Density_BC, Component_BC, Velocity_U_BC, Velocity_V_BC, Velocity_W_BC, Temperature_BC, &
           suf_u_bc_rob1, suf_v_bc_rob1, suf_w_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob2, suf_w_bc_rob2, &
           suf_t_bc_rob1, suf_t_bc_rob2, suf_vol_bc_rob1, suf_vol_bc_rob2, suf_comp_bc_rob1, suf_comp_bc_rob2, &
           theta_gdiff,  ScalarField_Source_Store, ScalarField_Source_Component, &
           mass_ele, dummy_ele, &
           Permeability, Material_Absorption, Material_Absorption_Stab, &
           Velocity_Absorption, ScalarField_Absorption, Component_Absorption, Temperature_Absorption, &
           Component_Diffusion_Operator_Coefficient, &
           Momentum_Diffusion, ScalarAdvectionField_Diffusion, &
           Component_Diffusion, &
           theta_flux, one_m_theta_flux, sum_theta_flux, sum_one_m_theta_flux )

      return
    end subroutine MultiFluids_SolveTimeLoop

  end module multiphase_time_loop
