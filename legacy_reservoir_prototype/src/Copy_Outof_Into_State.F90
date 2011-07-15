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

module copy_outof_into_state
  !! This module enables the multiphase prototype code to interact with state
  !! First copying everything required from state to the old variable space
  !! Second copying what is needed back to state for output 

  use fldebug
  use state_module
  use fields
  use spud
  use populate_state_module
  use diagnostic_variables
  use diagnostic_fields
  use diagnostic_fields_wrapper
  use global_parameters, only: option_path_len
  use diagnostic_fields_wrapper_new
!  use diagnostic_fields_new, only : &
!    & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
!    & check_diagnostic_dependencies

  use boundary_conditions

  implicit none
  
  private
  
  public :: copy_outof_state, copy_into_state
  
  contains
  
    subroutine copy_outof_state(state, dt, &
         nonlinear_iterations, nonlinear_iteration_tolerance, &
                                ! Begin here all the variables from read_scalar
         problem, nphases, ncomps, totele, ndim, nlev, &
         u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, &
         cv_snloc, u_snloc, p_snloc, stotel, &
         ncoef, nuabs_coefs, &
         u_ele_type, p_ele_type, mat_ele_type, cv_ele_type, &
         cv_sele_type, u_sele_type, &
         ntime, ntime_dump, nits, nits_internal, ndpset, noit_dim, &
         nits_flux_lim_volfra, nits_flux_lim_comp, &
         t_disopt, u_disopt, v_disopt, t_dg_vel_int_opt, &
         u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
         capil_pres_opt, ncapil_pres_coef, comp_diffusion_opt, ncomp_diff_coef, &
         patmos, p_ini, t_ini, t_beta, v_beta, t_theta, v_theta, &
         u_theta, domain_length, &
         lump_eqns, volfra_use_theta_flux, volfra_get_theta_flux, &
         comp_use_theta_flux, comp_get_theta_flux, &
                                ! Now the variables from read_all
         volfra_relax_number_iterations, scalar_relax_number_iterations, &
         global_relax_number_iterations,  velocity_relax_number_iterations, &
         pressure_relax_number_iterations, mass_matrix_relax_number_iterations, &
         in_ele_upwind, dg_ele_upwind, &
         volfra_error, volfra_relax, volfra_relax_diag, volfra_relax_row,  & 
         scalar_error, scalar_relax, scalar_relax_diag, scalar_relax_row, & 
         global_error, global_relax, global_relax_diag, global_relax_row, & 
         velocity_error, velocity_relax, velocity_relax_diag, velocity_relax_row, & 
         pressure_error, pressure_relax, pressure_relax_diag, pressure_relax_row, &
         mass_matrix_error, mass_matrix_relax, mass_matrix_relax_diag, mass_matrix_relax_row, &
         Mobility, alpha_beta, &
         KComp_Sigmoid, Comp_Sum2One, &
         wic_vol_bc, wic_d_bc, wic_u_bc, wic_p_bc, wic_t_bc, &
         wic_comp_bc, uabs_option, eos_option, cp_option, &
         suf_vol_bc, suf_d_bc, suf_cpd_bc, suf_t_bc, suf_p_bc, &
         suf_u_bc, suf_v_bc, suf_w_bc, suf_one_bc, suf_comp_bc, &
         suf_u_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob1, suf_v_bc_rob2, &
         suf_w_bc_rob1, suf_w_bc_rob2, suf_t_bc_rob1, suf_t_bc_rob2, &
         suf_vol_bc_rob1, suf_vol_bc_rob2, &
         suf_comp_bc_rob1, suf_comp_bc_rob2, &
         x, y, z, xu, yu, zu, nu, nv, nw, ug, vg, wg, &
         u_source, t_source, v_source, comp_source, &
         u, v, w, &
         den, satura, volfra, comp, t, p, cv_p, volfra_pore, &
         cv_one, Viscosity, & 
                                ! cv_one corresponds to density in single-phase advection problem
         uabs_coefs, &
         eos_coefs, cp_coefs, &
         comp_diff_coef, capil_pres_coef, &
         u_abs_stab, u_absorb, comp_absorb, &
         t_absorb, v_absorb, &
         perm, K_Comp, &
         comp_diffusion, &
                                ! Now adding other things which we have taken inside this routine to define
         cv_nonods, p_nonods, u_nonods, x_nonods, xu_nonods)

      !! New variables

      type(state_type), dimension(:), pointer :: state
      type(mesh_type) :: cmesh, vmesh, pmesh !! coordinate, velocity and pressure meshes

      type(scalar_field), pointer :: porosity, permeability, density, pressure, &
           phasevolumefraction, pvf_source, &
           scalarfield, scalarfield_source, &
           componentmassfraction, &
           phasevolumefraction_bc, density_bc, pressure_bc, scalarfield_bc

      type(vector_field), pointer :: velocity, velocity_source
      type(vector_field), pointer :: velocity_bc

      type(tensor_field), pointer :: viscosity_ph1, viscosity_ph2

      integer :: nonlinear_iterations, &  !! equal to nits in prototype code
           stat, nstates

      real :: dt
      real :: nonlinear_iteration_tolerance

      character(len=OPTION_PATH_LEN) :: field_name, material_phase_name, option_path

      !! temporary variables only needed for interfacing purposes

      type(vector_field), pointer :: positions

      integer :: i, j, k, l, nscalar_fields, cv_nonods, p_nonods, &
           x_nonods, xu_nonods, u_nonods

      real :: coord_min, coord_max, &
           eos_value!, viscosity_ph1, viscosity_ph2

      !! Variables needed by the prototype code
      !! and therefore needing to be pulled out of state or
      !! derived from information in state:

      ! Scalars (from read_scalar())
      integer :: problem, nphases, ncomps, totele, ndim, nlev, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, &
           cv_snloc, u_snloc, p_snloc, stotel, &
           ncoef, nuabs_coefs, &
           u_ele_type, p_ele_type, mat_ele_type, cv_ele_type, &
           cv_sele_type, u_sele_type

      integer :: ntime, ntime_dump, nits, nits_internal, ndpset, noit_dim, &
           nits_flux_lim_volfra, nits_flux_lim_comp

      integer :: t_disopt, u_disopt, v_disopt, t_dg_vel_int_opt, &
           u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt

      integer :: capil_pres_opt, ncapil_pres_coef, comp_diffusion_opt, ncomp_diff_coef

      real :: patmos, p_ini, t_ini, t_beta, v_beta, t_theta, v_theta, &
           u_theta, domain_length

      logical :: lump_eqns, volfra_use_theta_flux, volfra_get_theta_flux, &
           comp_use_theta_flux, comp_get_theta_flux


      ! Others (from read_all())

      integer :: volfra_relax_number_iterations, scalar_relax_number_iterations, &
           global_relax_number_iterations,  velocity_relax_number_iterations, &
           pressure_relax_number_iterations, mass_matrix_relax_number_iterations, &
           in_ele_upwind, dg_ele_upwind

      real :: volfra_error, volfra_relax, volfra_relax_diag, volfra_relax_row,  & 
           scalar_error, scalar_relax, scalar_relax_diag, scalar_relax_row, & 
           global_error, global_relax, global_relax_diag, global_relax_row, & 
           velocity_error, velocity_relax, velocity_relax_diag, velocity_relax_row, & 
           pressure_error, pressure_relax, pressure_relax_diag, pressure_relax_row, &
           mass_matrix_error, mass_matrix_relax, mass_matrix_relax_diag, mass_matrix_relax_row, &
           Mobility, alpha_beta

      logical :: KComp_Sigmoid, Comp_Sum2One

      integer, dimension( : ), allocatable :: wic_vol_bc, wic_d_bc, wic_u_bc, wic_p_bc, wic_t_bc, &
           wic_comp_bc, uabs_option, eos_option, cp_option

      real, dimension( : ), allocatable :: suf_vol_bc, suf_d_bc, suf_cpd_bc, suf_t_bc, suf_p_bc, &
           suf_u_bc, suf_v_bc, suf_w_bc, suf_one_bc, suf_comp_bc, &
           suf_u_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob1, suf_v_bc_rob2, &
           suf_w_bc_rob1, suf_w_bc_rob2, suf_t_bc_rob1, suf_t_bc_rob2, &
           suf_vol_bc_rob1, suf_vol_bc_rob2, &
           suf_comp_bc_rob1, suf_comp_bc_rob2, &
           x, y, z, xu, yu, zu, nu, nv, nw, ug, vg, wg, &
           u_source, t_source, v_source, comp_source, &
           u, v, w, &
           den, satura, volfra, comp, t, p, cv_p, volfra_pore, &
           cv_one, Viscosity

      real, dimension( :, : ), allocatable :: uabs_coefs
      real, dimension( :, : ), allocatable :: eos_coefs, cp_coefs

      real, dimension( : , : , : ), allocatable :: comp_diff_coef, capil_pres_coef, &
           u_abs_stab, u_absorb, comp_absorb, &
           t_absorb, v_absorb, &
           perm, K_Comp

      real, dimension( : , : , : , : ), allocatable :: comp_diffusion

      integer :: Velocity_BC_Type, Pressure_BC_Type, density_bc_type, &
           component_bc_type, pvf_bc_type, Temperature_bc_type, shape_option(2)

      real :: component_suf_bc

      integer, dimension(:), allocatable :: Velocity_SufID_BC, Pressure_SufID_BC, density_sufid_bc, &
           component_sufid_bc, pvf_sufid_bc, Temperature_sufid_bc

      ! Gravity terms to be linked with u_source
      logical :: have_gravity
      !!  real :: gravity_magnitude, gravity_direction, delta_den
      real :: gravity_magnitude, delta_den, grm
      type(vector_field) :: gravity_direction
      !      type( vector_field ), pointer :: gravity
      !      type(vector_field), pointer :: dummyvector

      integer :: nobcs    

      !! Finish declaration of variables needed from user input

      ewrite(3,*) 'In copy_outof_state'

      !! Here are all of the needed items
      !! I suggest we work through them one by one
      !! according to Alex and Brendan's markup of the input file
      !! Shuffle the order as necessary, but it would be worth keeping
      !! related items together as much as possible. JHS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! these does not matter -- most of them can be deleted quite soon
      problem = 1
      nlev = 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      nstates = option_count("/material_phase")
      ewrite(3,*) ' nstates:', nstates
      !! Assume there are the same number of components in each phase (will need to check this eventually)
      nscalar_fields = option_count("/material_phase[0]/scalar_field")
      ncomps=0
      do i=1,nstates
         if (have_option("/material_phase[" // int2str(i-1) // "]/is_multiphase_component")) then
            ncomps=ncomps+1
         end if
      end do
      nphases=nstates-ncomps
      !! Going to assume for now that phases have been inserted into state before components
      ewrite(3,*) 'nphases, ncomps:',nphases, ncomps

      !! Let's get the meshes out here, as we're going to need them for a whole
      !! load of things:
      cmesh = extract_mesh(state, "CoordinateMesh")
      if (have_option("/geometry/mesh::VelocityMesh")) then
         vmesh = extract_mesh(state, "VelocityMesh")
      else
         vmesh = cmesh
      endif
      if (have_option("/geometry/mesh::PressureMesh")) then
         pmesh = extract_mesh(state, "PressureMesh")
      else
         pmesh = cmesh
      endif
      positions => extract_vector_field(state, "Coordinate")

      totele = ele_count(cmesh)

      call get_option("/geometry/dimension",ndim)

      ! nlev = 
      xu_nloc = cmesh%shape%loc
      cv_nloc = pmesh%shape%loc
      x_nloc = pmesh%shape%loc
      p_nloc = pmesh%shape%loc ! 3
      u_nloc = vmesh%shape%loc
      !      if (have_option("/material_phase::phase1/vector_field::Velocity/prognostic/" // &
      !           "spatial_discretisation/discontinuous_galerkin/overlapping")) then
      !         u_nloc=u_nloc*p_nloc
      !      endif
      if (pmesh%continuity>=0) then
         ! Continuous pressure mesh
         cv_nonods = ( cv_nloc - 1 ) * totele + 1
      else
         ! Discontinuous pressure mesh
         cv_nonods = cv_nloc * totele
      endif
      p_nonods = cv_nonods
      cv_snloc = 1
      u_snloc = 3
      p_snloc = 1
      ! x_snloc = 1
      stotel = surface_element_count(cmesh)

      !! EoS things are going to be done very differently
      !! Currently the default value of ncoef is 10 so I'm going
      !! to put that here until we have a more concrete idea of what
      !! we're dealing with
      ncoef = 10

      !! I don't understand why absorption needs a coefficient, or what it is
      !! but it's equal to 1 in all the test cases
      nuabs_coefs = 1

      call get_option('/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/polynomial_degree', u_ele_type, default=1)
      call get_option('/geometry/mesh::PressureMesh/from_mesh/mesh_shape/polynomial_degree', p_ele_type, default=1)
!      ewrite(3,*) 'u_ele_type', u_ele_type
      u_ele_type = 2 ! this will need to be changed later -- switcher for the 
      !! Sufficient to set this to 1 for now? It is 1 in all test cases
      mat_ele_type = 1
      !! Sufficient to set this to 2 for now? It is 2 in all test cases
      !! which presumably means discontinuous
      !! Will need to update once schema is changed
      cv_ele_type = p_ele_type !2
      !! These aren't used 
      cv_sele_type = 1
      u_sele_type = 1

      if( nphases == 1 ) then
         problem = 0 ! Single-phase advection (continuous)
         if( pmesh%continuity < 0 ) problem = -1  
         cv_ele_type = 1
         u_ele_type = 1
      end if


      !! Time options
      ewrite(3,*) ' Getting time options'
      call get_option('/io/max_dump_file_count', ntime, default=160)
      if (have_option('/io/dump_period_in_timesteps/constant')) then
         call get_option('/io/dump_period_in_timesteps/constant',ntime_dump)
      else
         ! This way might be inaccurate due to rounding errors:
         call get_option('/io/dump_period/constant',ntime_dump)
         ntime_dump = int(ntime_dump/dt)
      end if

      ewrite(3,*) ' Getting iteration info'
      call get_option( '/timestepping/nonlinear_iterations', nonlinear_iterations, &
           default = 3 )
      nits = nonlinear_iterations

      if( have_option( '/timestepping/nonlinear_iterations/tolerance' )) then
         call get_option( '/timestepping/nonlinear_iterations/tolerance', &
              nonlinear_iteration_tolerance )
      else
         nonlinear_iteration_tolerance = 1.e-6
      end if

      ! This one is only for compositional problems
      call get_option('/material_phase[0]/scalar_field::component1/&
           &prognostic/temporal_discretisation/control_volumes/number_advection_iterations',&
           &nits_internal, default=1)

      ndpset = 0
      noit_dim = 5
      !! These two are always equal to 3 except for the 2phase, 3comp, non-equilibrium test, where they are 1
      !! Change in schema ??
      !! still need to do this!
      call get_option("/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/" // &
           "max_iterations_flux_limiter", nits_flux_lim_volfra, default=3)
      call get_option("/material_phase[" // int2str(nphases) // "]/scalar_field::Phase1ComponentMassFraction/" // &
           "prognostic/max_iterations_flux_limiter", nits_flux_lim_comp, default=3)


      !! disopt options: going to need to change the schema I think
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

      if (have_option('/material_phase[0]/scalar_field::Temperature/prognostic/' // &
           'spatial_discretisation/control_volumes/face_value::FiniteElement/limit_face_value')) then
         t_disopt = 8 !! Unless all the other options, but need to be able to get 8 here
      else
         t_disopt = 1
      endif

      u_disopt = 1 ! Ditto, except that this probably IS used, being velocity. Hmm.

      if (have_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'spatial_discretisation/control_volumes/face_value::FirstOrderUpwind')) then
         v_disopt = 0 ! Unless theta=non_linear ???
      elseif (have_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'spatial_discretisation/control_volumes/face_value::Trapezoidal')) then
         v_disopt = 2 ! Unless theta=non_linear ???
      elseif (have_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'spatial_discretisation/control_volumes/face_value::FiniteElement/do_not_limit_face_value')) then
         v_disopt = 6 ! Unless theta=non_linear ???
      elseif (have_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'spatial_discretisation/control_volumes/face_value::FiniteElement/limit_face_value')) then
         v_disopt = 8 !! Unless all the other options, but need to be able to get 8 here
      endif

      t_dg_vel_int_opt = 0
      u_dg_vel_int_opt = 4 ! Not used -- it can be deleted
      v_dg_vel_int_opt = 4
      w_dg_vel_int_opt = 0 ! Not used -- it can be deleted

      ewrite(3,*) ' Getting capillary pressure options'
      if (have_option('/porous_media/multiphase_parameters/cp_A')) then
         ! 1 is the only option available at the moment...
         capil_pres_opt = 1
         ! and even this one doesn't depend on any coefficients!
         ncapil_pres_coef = 0
      else
         capil_pres_opt = 0
         ncapil_pres_coef = 0
      end if

      !! These are not currently used in any of the code
      comp_diffusion_opt = 0
      ncomp_diff_coef = 0

      call get_option('/material_phase[0]/scalar_field::Pressure/prognostic/&
           &atmospheric_pressure', patmos, default=0.0)
      call get_option('/material_phase[0]/scalar_field::Pressure/prognostic/&
           &initial_condition::WholeMesh/constant',p_ini, default=0.0)

      t_ini = 0.0
      if (nscalar_fields>2) then ! we might have an extra scalar_field so might need t
         call get_option('/material_phase[0]/scalar_field::Temperature/prognostic/' // &
              'spatial_discretisation/conservative_advection', t_beta, default=0.)
      end if

      call get_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'spatial_discretisation/conservative_advection', v_beta)

      if (nscalar_fields>2) then ! we might have an extra scalar_field so might need t
         call get_option('/material_phase[0]/scalar_field::Temperature/prognostic/' // &
              'temporal_discretisation/theta', t_theta, default=0.)
      end if
      call get_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'temporal_discretisation/theta', v_theta)
      call get_option('/material_phase[0]/vector_field::Velocity/prognostic/' // &
           'temporal_discretisation/theta', u_theta)

      ! This is a strictly 1d property, so will have a suitably 1d method for finding it in state!
      coord_min=1.0e9
      coord_max=-1.0e-9
      do i=1,node_count(positions)
         coord_min=min(coord_min,positions%val(X_,i))
         coord_max=max(coord_max,positions%val(X_,i))
      end do
      domain_length = coord_max - coord_min

      !! I'm going to get this one from the PhaseVolumeFraction scalar_field
      lump_eqns = have_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'spatial_discretisation/continuous_galerkin/mass_terms/lump_mass_matrix')

      volfra_use_theta_flux = .FALSE.
      volfra_get_theta_flux = .TRUE.
      comp_use_theta_flux = .FALSE.
      comp_get_theta_flux = .TRUE.

      ewrite(3,*) 'Finished stuff from read_scalar'


      ! Others (from read_all())

      !!
      !! Most of the options below should be replaced by PETSc options soon, i.e., still at the second part of
      !! of the first stage of the integration.  So let's keep as it is for the moment and remove them
      !! as soon as we have all PETSc data structure enabled.

      Conditional_VolumeFraction_Solver: if( have_option( '/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic' )) then
         call get_option( '/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
              'solver/max_iterations', volfra_relax_number_iterations, default = 100 )
         call get_option( '/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
              'solver/relative_error', volfra_error, default = 1.e-5 )
         volfra_relax = 1.
         volfra_relax_diag = 0.
         volfra_relax_row = 1.
      endif Conditional_VolumeFraction_Solver

      Conditional_Pressure_Solver: if( have_option( '/material_phase[0]/scalar_field::Pressure' )) then
         call get_option( '/material_phase[0]/scalar_field::Pressure/prognostic/solver/max_iterations', &
              pressure_relax_number_iterations, default = 4000 )
         call get_option( "/material_phase[0]/scalar_field::Pressure/prognostic/solver/relative_error", &
              pressure_error, default = 1.e-3 )
         pressure_relax = 1.
         pressure_relax_diag = 0.
         pressure_relax_row = 1.
      endif Conditional_Pressure_Solver

      Conditional_Velocity_Solver: if( have_option( "/material_phase[0]/vector_field::Velocity" )) then
         call get_option( "/material_phase[0]/vector_field::Velocity/prognostic/solver/max_iterations", &
              velocity_relax_number_iterations, default = 100 )
         call get_option( "/material_phase[0]/vector_field::Velocity/prognostic/solver/relative_error", &
              velocity_error, default = 1.e-5 ) 
         velocity_relax = 1.
         velocity_relax_diag = 0.
         velocity_relax_row = 1.
      end if Conditional_Velocity_Solver

      ewrite(3,*) "Got first lot of solver options"

      scalar_relax_number_iterations = 100
      global_relax_number_iterations = 100
      mass_matrix_relax_number_iterations = 200 

      scalar_error = 1.e-10
      global_error = 1.e-10
      mass_matrix_error = 1.e-10

      scalar_relax = 1.
      global_relax = 1.
      mass_matrix_relax = 1. 

      scalar_relax_diag = 0.
      global_relax_diag = 0.
      mass_matrix_relax_diag = 0. 

      scalar_relax_row = 1.
      global_relax_row = 1.
      mass_matrix_relax_row = 1. 

      ! IN/DG_ELE_UPWIND are options for optimisation of upwinding across faces in the overlapping
      ! formulation. The data structure and options for this formulation need to be added later. 
      in_ele_upwind = 3
      dg_ele_upwind = 3

      ! Calculating Mobility
      ewrite(3,*) "going to get viscosities"

      ! Also have to allocate and initialise Viscosity
      allocate(Viscosity( cv_nonods * nphases ))
      Viscosity = 0.

      viscosity_ph1 => extract_tensor_field(state(1), "Viscosity")

!!!
!!! This will be changed later, as we just need the viscosity of one phase 
!!! and the mobility. All test cases will need to be updated, when we move to 
!!! this. Also for the single phase advection, both mobility and viscosity need
!!! to be initialised although will not be used.     
      if( have_option( "/physical_parameters/mobility" ))then
         call get_option( "/physical_parameters/mobility", Mobility )
         Viscosity( 1 : cv_nonods ) = Viscosity_Ph1%val( 1, 1, 1 )
      elseif ( have_option( "/material_phase[1]/vector_field::Velocity/prognostic/" // &
           "tensor_field::Viscosity/prescribed/value::WholeMesh/" // &
           "isotropic")) then
         viscosity_ph2 => extract_tensor_field(state(2), "Viscosity")
         Mobility =  Viscosity_Ph2%val( 1, 1, 1 ) / Viscosity_Ph1%val( 1, 1, 1 )         
      elseif( nphases == 1 ) then
         Mobility = 0. 
      end if

      ewrite(3,*) "Got viscosities & mobility"

!!!
!!! Options below are for the multi-component flow model
!!! 
      call get_option("/material_phase[" // int2str(nstates-ncomps) // &
           "]/is_multiphase_component/alpha_beta", alpha_beta, default=1.0)
      KComp_Sigmoid = have_option("/material_phase[" // int2str(nstates-ncomps) // &
           "]/is_multiphase_component/KComp_Sigmoid")
      Comp_Sum2One = have_option("/material_phase[" // int2str(nstates-ncomps) // &
           "]/is_multiphase_component/Comp_Sum2One")

!!!
!!! Porosity and Permeability: it WILL be necessary to change the permeability as it
!!! is defined in the PC as a tensor with dimension ( totele, ndim, ndim )
!!!
      porosity => extract_scalar_field(state, "Porosity")
      allocate(volfra_pore(totele))
      ! The porosity will be in element order in 1d
      do i=1,totele
         volfra_pore(i)=porosity%val(i)
      enddo
      ewrite(3,*) "Got porosity"

      if (have_option("/porous_media/scalar_field::Permeability")) then
         permeability => extract_scalar_field(state(1), "Permeability")
         allocate(perm(totele, ndim, ndim))
         perm = 0.
         do i=1,totele
            do j=1,ndim
               do k=1,ndim
                  perm(i,j,k)=permeability%val(i)
               end do
            end do
         end do
      elseif (have_option("/porous_media/tensor_field::Permeability")) then
         FLAbort("Have not coded up tensor permeability yet! Try scalar instead")
      endif
      ewrite(3,*) "Got permeability"

!!!
!!! WIC_X_BC (in which X = D, U, V, W, P, T, COMP and VOL) controls the boundary conditions
!!! type applied. == 1 (Dirichlet), = 2 (Robin), = 3 (Newman) 
!!!
!!!
!!! Components Boundary Conditions

      allocate( wic_comp_bc( stotel * nphases ))
      allocate( suf_comp_bc( stotel * 1 * nphases * ncomps ))
      wic_comp_bc = 0
      suf_comp_bc = 0.

      Loop_Component_BC: do i=nphases, nphases+ncomps-1
         do j=1,nphases
            if( have_option("/material_phase[" // int2str(i) // "]/scalar_field::Phase" // int2str(j) // &
                 "ComponentMassFraction/prognostic/" // &
                 "boundary_conditions[0]/type::dirichlet" )) then

               shape_option=option_shape("/material_phase[" // int2str(i) // &
                    "]/scalar_field::Phase" // int2str(j) // &
                    "ComponentMassFraction/prognostic/boundary_conditions[0]/surface_ids")
               allocate(component_sufid_bc(1:shape_option(1)))
               Component_BC_Type = 1
               call get_option( "/material_phase[" // int2str(i) // "]/scalar_field::Phase" // int2str(j) // &
                    "ComponentMassFraction/prognostic/" // &
                    "boundary_conditions[0]/surface_ids", component_sufid_bc )
               call get_option( "/material_phase[" // int2str(i) // "]/scalar_field::Phase" // int2str(j) // &
                    "ComponentMassFraction/prognostic/" // &
                    "boundary_conditions[0]/type::dirichlet/constant", Component_Suf_BC )
                    
               ewrite(3,*) 'csufid', component_sufid_bc
               ewrite(3,*) 'nphases, j, ncomps, i', nphases, j, ncomps, i

               do k=1,shape_option(1)
                  wic_comp_bc( component_sufid_bc(1) + nphases*(j-1) ) = Component_BC_Type
                  suf_comp_bc( component_sufid_bc(1) + nphases*(j-1) + nphases*ncomps*(i-nphases) ) = Component_Suf_BC
               enddo
               deallocate(Component_sufid_bc)
            endif

         enddo

      enddo Loop_Component_BC

      allocate( suf_cpd_bc( stotel * 1 * nphases ))
      suf_cpd_bc = 0.      

      ! Extra allocation to stop the code dying when also running everything through the old io
      allocate( suf_one_bc( stotel * cv_snloc * nphases ))
      suf_one_bc=0.

      ewrite(3,*) "Done with boundary conditions"

!!!
!!! Robin Boundary conditions need to be added at a later stage
!!!
      allocate( suf_u_bc_rob1( stotel * 3 * nphases ))
      allocate( suf_v_bc_rob1( stotel * 3 * nphases ))
      allocate( suf_w_bc_rob1( stotel * 3 * nphases ))
      allocate( suf_u_bc_rob2( stotel * 3 * nphases ))
      allocate( suf_v_bc_rob2( stotel * 3 * nphases ))
      allocate( suf_w_bc_rob2( stotel * 3 * nphases ))
      allocate( suf_t_bc_rob1( stotel * nphases ))
      allocate( suf_t_bc_rob2( stotel * nphases ))
      allocate( suf_vol_bc_rob1( stotel * nphases ))
      allocate( suf_vol_bc_rob2( stotel * nphases ))
      allocate( suf_comp_bc_rob1( stotel * nphases ))
      allocate( suf_comp_bc_rob2( stotel * nphases ))
      suf_u_bc_rob1 = 0.
      suf_u_bc_rob2 = 0.
      suf_v_bc_rob1 = 0.
      suf_v_bc_rob2 = 0.
      suf_w_bc_rob1 = 0.
      suf_w_bc_rob2 = 0.
      suf_t_bc_rob1 = 0.
      suf_t_bc_rob2 = 0.
      suf_vol_bc_rob1 = 0.
      suf_vol_bc_rob2 = 0.
      suf_comp_bc_rob1 = 0.
      suf_comp_bc_rob2 = 0.

!!!
!!!  End of Boundary Conditions Section
!!!


!!!
!!! Initial conditions for all fields:
!!! Density and pressure might need re-ordering because of the
!!! peculiarity of the quadratic element setup
!!! see copy_into_state() below
!!!

!!!
!!!  Density
!!!
      ewrite(3,*) "going to get density..."

      Loop_Density: do i = 1, nphases
         density => extract_scalar_field( state( i ), "Density")
         if ( .not. allocated( den )) allocate( den( nphases * node_count( density )))

         do j = node_count( density ), 1, -1
            den(( i - 1 ) * node_count( density ) + j ) = density%val( j )
            ! This will make sure that the fields in the PC *does not* contain any boundary conditions
            ! elements. This need to be changed along with the future data structure to take into 
            ! account the overlapping formulation.
            if( j == 1 ) den(( i - 1 ) * node_count( density ) + j ) = &
                 den(( i - 1 ) * node_count( density ) + j + 1)     
         end do

         if( .not. ( allocated( wic_d_bc ) .and. allocated( suf_d_bc ))) then 
            allocate( wic_d_bc( stotel * nphases ))
            allocate( suf_d_bc( stotel * 1 * nphases ))
            wic_d_bc = 0
            suf_d_bc = 0.
         end if

         Conditional_Density_BC: if( have_option( '/material_phase[' // int2str(i-1) // &
              ']/scalar_field::Density/prognostic/' // &
              'boundary_conditions[0]/type::dirichlet' )) then

            shape_option=option_shape('/material_phase[' // int2str(i-1) // &
                 ']/scalar_field::Density/prognostic/boundary_conditions[0]/' // &
                 'surface_ids' )
            allocate( density_sufid_bc( 1 : shape_option( 1 )))

            Density_BC_Type = 1

            call get_option( '/material_phase[' // int2str(i-1) // &
                 ']/scalar_field::Density/prognostic/' // &
                 'boundary_conditions[0]/surface_ids', Density_SufID_BC )

            do j=1,shape_option(1)
               wic_d_bc( density_sufid_bc( 1 ) + ( i - 1 ) * nphases ) = density_bc_type
            enddo

            nobcs = get_boundary_condition_count( density )
            do j = 1, nobcs
               density_bc => extract_surface_field( density, j, "value" )
            end do
            do j = 1, node_count( density_bc )
               suf_d_bc( ( i - 1 ) * stotel + j ) = density_bc%val( j )
            end do

            deallocate(density_sufid_bc)

         endif Conditional_Density_BC

      enddo Loop_Density

!!!
!!! Pressure
!!!
      ewrite(3,*) "pressure..."
      pressure => extract_scalar_field( state( 1 ), "Pressure" )
      allocate( p( node_count( pressure )))

      ! This will make sure that the fields in the PC *does not* contain any boundary conditions
      ! elements. This need to be changed along with the future data structure to take into 
      ! account the overlapping formulation.
      do i= node_count( pressure ), 1, -1
         p( i ) = pressure%val( i )
         if( i == 1 ) p ( i ) = p( i + 1 )
      enddo

      !! Control-volume pressure used for shock tube initialisation
      allocate( cv_p( node_count( pressure )))
      cv_p = p

      Loop_Pressure: do i = 1, nphases

         if( .not. ( allocated( wic_p_bc ) .and. allocated( suf_p_bc ))) then
            allocate( wic_p_bc( stotel * nphases ))
            allocate( suf_p_bc( stotel * 1 * nphases ))
            wic_p_bc = 0
            suf_p_bc = 0.
         end if

         Conditional_Pressure_BC: if( have_option( '/material_phase[' // int2str(i-1) // &
              ']/scalar_field::Pressure/prognostic/' // &
              'boundary_conditions[0]/type::dirichlet' )) then

            shape_option=option_shape('/material_phase[' // int2str(i-1) // &
                 ']/scalar_field::Pressure/' // &
                 'prognostic/boundary_conditions[0]/surface_ids')

            allocate( pressure_sufid_bc( 1 : shape_option( 1 )))

            Pressure_BC_Type = 1
            call get_option( '/material_phase[' // int2str(i-1) // &
                 ']/scalar_field::Pressure/prognostic/' // &
                 'boundary_conditions[0]/surface_ids', Pressure_SufID_BC )

            do j = 1, shape_option( 1 )
               wic_p_bc( pressure_sufid_bc( 1 ) + ( i - 1 ) * nphases ) = pressure_bc_type
               ! The bellow is done as pressure for phase 2 is aliased therefore the same info for nodes
               ! for the phase 1 should be copied for phase 2. This need to be changed later.
               if( nphases > 1 ) &
                    wic_p_bc( pressure_sufid_bc( 1 ) + i * nphases ) = pressure_bc_type
            enddo

            nobcs = get_boundary_condition_count( pressure )
            do j = 1, nobcs
               pressure_bc => extract_surface_field( pressure, j, "value" )
            end do
            do j = 1, node_count( pressure_bc )
               suf_p_bc( ( i - 1 ) * stotel + j ) = pressure_bc%val( j )
            end do

         end if Conditional_Pressure_BC

      end do Loop_Pressure

!!!
!!! Volume Fraction (or Saturation) and associated boundary conditions:
!!!
      ewrite(3,*) "phasevolumefraction..."

      Loop_VolumeFraction: do i = 1, nphases

         phasevolumefraction => extract_scalar_field(state(i), "PhaseVolumeFraction")

         if ( .not. allocated( satura )) allocate( satura( nphases * node_count( phasevolumefraction )))

         ! This will make sure that the fields in the PC *does not* contain any boundary conditions
         ! elements. This need to be changed along with the future data structure to take into 
         ! account the overlapping formulation.
         do j = node_count( phasevolumefraction ), 1, -1
            satura( ( i - 1 ) * node_count( phasevolumefraction ) + j ) = phasevolumefraction%val( j )
            if( j == 1 ) satura( ( i - 1 ) * node_count( phasevolumefraction ) + j ) = &
                 satura( ( i - 1 ) * node_count( phasevolumefraction ) + j + 1)
         enddo

         if (.not. allocated(wic_vol_bc)) then
            allocate( wic_vol_bc( stotel * nphases ))
            wic_vol_bc = 0.
         endif

         if (.not. allocated(suf_vol_bc)) then
            allocate( suf_vol_bc( stotel * 1 * nphases ))
            suf_vol_bc = 0.
         endif

         Conditional_VolumeFraction_BC: if( have_option( "/material_phase[" // int2str(i-1) // &
              "]/scalar_field::PhaseVolumeFraction/" // &
              "prognostic/boundary_conditions[0]/type::dirichlet" )) then

            shape_option = option_shape( "/material_phase[" // int2str(i-1) // "]/scalar_field::" // &
                 "PhaseVolumeFraction/prognostic/boundary_conditions[0]/surface_ids" )
            allocate( pvf_sufid_bc( 1 : shape_option( 1 )))

            pvf_bc_type = 1

            call get_option( "/material_phase[" // int2str(i-1) //"]/scalar_field::" // &
                 "PhaseVolumeFraction/prognostic/boundary_conditions[0]/surface_ids", pvf_sufid_bc )
            do j = 1, shape_option(1)
               wic_vol_bc( pvf_sufid_bc(1) + ( i - 1 ) * nphases ) = pvf_bc_type
            enddo

            nobcs = get_boundary_condition_count( phasevolumefraction )
            do j = 1, nobcs
               phasevolumefraction_bc => extract_surface_field( phasevolumefraction, j, "value" )
            end do
            do j = 1, node_count( phasevolumefraction_bc )
               suf_vol_bc( ( i - 1 ) * stotel + j ) = phasevolumefraction_bc%val( j )
            end do

            deallocate(pvf_sufid_bc)

         endif Conditional_VolumeFraction_BC

      enddo Loop_VolumeFraction

      ! Also need to allocate and initialise volfra
      allocate(volfra( cv_nonods * nphases ))
      volfra=0.

!!!
!!! Velocity and associated boundary conditions:
!!!
      ewrite(3,*) "velocity..."

      Loop_Velocity: do i = 1, nphases

         velocity => extract_vector_field(state(i), "Velocity")
         if (.not. allocated(u)) then
            allocate(u(nphases*node_count(velocity)))
            allocate(v(nphases*node_count(velocity)))
            allocate(w(nphases*node_count(velocity)))
         endif
         u=0.
         v=0.
         w=0.
         ewrite(3,*)'nodes veloc:', node_count(velocity)

         ! This will make sure that the fields in the PC *does not* contain any boundary conditions
         ! elements. This need to be changed along with the future data structure to take into 
         ! account the overlapping formulation.
         do j = node_count(velocity), 1, -1
            u((i-1)*node_count(velocity)+j)=velocity%val(X_, j)
            if( j == 1 ) u( ( i - 1 ) * node_count( velocity ) + j ) = &
                 u( ( i - 1 ) * node_count( velocity ) + j + 1)
            if (ndim>1) then
               v((i-1)*node_count(velocity)+j)=velocity%val(Y_, j)
               if( j == 1 ) v( ( i - 1 ) * node_count( velocity ) + j ) = &
                    v( ( i - 1 ) * node_count( velocity ) + j + 1)
            endif
            if (ndim>2) then
               w((i-1)*node_count(velocity)+j)=velocity%val(Z_, j)
               if( j == 1 ) w( ( i - 1 ) * node_count( velocity ) + j ) = &
                    w( ( i - 1 ) * node_count( velocity ) + j + 1)
            endif
         enddo

         if (.not.allocated(wic_u_bc)) then
            allocate( wic_u_bc( stotel * nphases ))
            wic_u_bc = 0
         endif
         if (.not.allocated(suf_u_bc)) then
            allocate( suf_u_bc( stotel * 3 * nphases ))
            suf_u_bc = 0.
         endif
         if (.not.allocated(suf_v_bc)) then
            allocate( suf_v_bc( stotel * 3 * nphases ))
            suf_v_bc = 0.
         endif
         if (.not.allocated(suf_w_bc)) then
            allocate( suf_w_bc( stotel * 3 * nphases ))
            suf_w_bc = 0.
         endif

         Conditional_Velocity_BC: if( have_option( '/material_phase[' // int2str(i-1) // &
              ']/vector_field::Velocity/prognostic/' // &
              'boundary_conditions[0]/type::dirichlet' )) then

            shape_option=option_shape("/material_phase[" // int2str(i-1) // "]/vector_field::Velocity/&
                 &prognostic/boundary_conditions[0]/surface_ids")
            if( .not. allocated( velocity_sufid_bc ))allocate(velocity_sufid_bc(1:shape_option(1)))

            Velocity_BC_Type = 1
            call get_option( "/material_phase[" // int2str(i-1) // "]/vector_field::Velocity/" // &
                 "prognostic/boundary_conditions[0]/surface_ids", Velocity_SufID_BC )

            do j = 1, shape_option(1)
               wic_u_bc( velocity_sufid_bc(1) + ( i - 1 ) * nphases ) = Velocity_BC_Type
            enddo

            nobcs = get_boundary_condition_count( velocity )
            do j = 1, nobcs
               velocity_bc => extract_surface_field( velocity, j, "value" )
            end do
            do j = 1, node_count( velocity_bc )
               suf_u_bc( ( i - 1 ) * stotel + j ) = velocity_bc%val( 1, j )
               if (velocity%dim>1) suf_v_bc( ( i - 1 ) * stotel + j ) = velocity_bc%val( 2, j )
               if (velocity%dim>2) suf_w_bc( ( i - 1 ) * stotel + j ) = velocity_bc%val( 3, j )
            end do

         endif Conditional_Velocity_BC

      enddo Loop_Velocity


      allocate(uabs_option(nphases))
      allocate(eos_option(nphases))
      allocate(cp_option(nphases))
      do i=1,nphases
         option_path = "/material_phase[" // int2str(i-1) // "]/multiphase_properties/relperm_type"
         if (have_option(trim(option_path)//"/Corey")) uabs_option(i)=3
         if (have_option(trim(option_path)//"/Corey/boost_at_zero_saturation")) uabs_option(i)=5

         option_path = "/material_phase[" // int2str(i-1) // "]/equation_of_state"
         if (have_option(trim(option_path)//"/incompressible/linear")) then
            eos_option(i) = 2
         elseif (have_option(trim(option_path)//"/compressible/stiffened_gas")) then
            eos_option(i) = 1
         elseif (have_option(trim(option_path)//"/compressible/exponential_oil_gas")) then
            eos_option(i) = 3
         else
            FLAbort("Unknown EoS option for phase "// int2str(i))
         endif
         cp_option(i) = 1
      enddo

      !! uabs_coefs is currently only used in rel perm options which aren't
      !! selected in any of the test cases, ie the 'standard polynomial
      !! representation of relative permeability'
      allocate(uabs_coefs(nphases, nuabs_coefs))
      allocate(eos_coefs(nphases, ncoef))
      eos_coefs=0.
      !! Capillary pressure isn't used at all at the moment
      allocate(cp_coefs(nphases, nphases))
      do i=1,nphases
         uabs_coefs(i,1) = 1.
         if (eos_option(i)==2) then
            if (have_option("/material_phase[" // int2str(i-1) // "]/equation_of_state/incompressible/linear/all_equal")) then
               call get_option("/material_phase[" // int2str(i-1) // "]/equation_of_state/incompressible/linear/all_equal", eos_value)
               eos_coefs(i,1) = eos_value
            else
               call get_option("/material_phase[" // int2str(i-1) // "]/equation_of_state/incompressible/linear/specify_all", eos_coefs(i, :))
            endif
         elseif (eos_option(i)==1) then
            call get_option("/material_phase[" // int2str(i-1) // "]/equation_of_state/compressible/stiffened_gas/eos_option1", eos_coefs(i, 1))
            call get_option("/material_phase[" // int2str(i-1) // "]/equation_of_state/compressible/stiffened_gas/eos_option2", eos_coefs(i, 2))
            eos_coefs(i, 3:ncoef) = 0.
         endif
         ewrite(3,*) 'i, eos_coefs', i, eos_coefs(i, :)
      enddo
      cp_coefs = 1.

      ! These are (nearly) all initialised to zero in the input files
      x_nonods = max(( x_nloc - 1 ) * totele + 1, totele )
      xu_nonods = max(( xu_nloc - 1 ) * totele + 1, totele )
      allocate(x(x_nonods))
      allocate(y(x_nonods))
      allocate(z(x_nonods))
      x=0.
      y=0.
      z=0.
      allocate(xu(xu_nonods))
      allocate(yu(xu_nonods))
      allocate(zu(xu_nonods))
      xu=0.
      yu=0.
      zu=0.
      u_nonods = u_nloc * totele

      allocate( nu( u_nonods * nphases ))
      allocate( nv( u_nonods * nphases ))
      allocate( nw( u_nonods * nphases ))
      allocate( ug( u_nonods * nphases ))
      allocate( vg( u_nonods * nphases ))
      allocate( wg( u_nonods * nphases ))
      nu=1.
      nv=0.
      nw=0.
      if (ndim>1) nv=1.
      if (ndim>2) nw=1.
      ug=0.
      vg=0.
      wg=0.

      ewrite(3,*) 'Getting source terms -- gravity '
      ! Gravity is associated with the u_source term
      call get_option( "/physical_parameters/gravity/magnitude", gravity_magnitude, stat )
      have_gravity = ( stat == 0 )

      if( have_gravity ) then
         gravity_direction = extract_vector_field(state(1), 'GravityDirection', stat )
         ! Normalise direction vector
         grm=0
         do i=1,ndim
            grm=grm + gravity_direction%val(i,1)**2
         end do
         do i=1,ndim
            gravity_direction%val(i,:) = gravity_direction%val(i,:)/sqrt(grm)
         end do
      end if

      !!if( have_gravity ) then
      !!   if( have_option( '/physical_parameters/gravity/vector_field::' // &
      !!        'GravityDirection/prescribed/value::WholeMesh' ))then
      !!      call get_option( '/physical_parameters/gravity/vector_field::' // &
      !!           'GravityDirection/prescribed/value::WholeMesh/constant', &
      !!           gravity_direction, stat )
      !!   end if
      !!end if

      ewrite(3, *)"Getting source terms -- velocity "
      allocate( u_source( ndim * u_nonods * nphases ))
      u_source = 0.
      Conditional_VelocitySource: if( have_option( '/material_phase[0]/vector_field::Velocity/' // &
           'prognostic/vector_field::Source' )) then 
         ! This is still not working as the length of node_count(velocity_source) =
         ! node_count(velocity) /= u_nonods
         do i=1,nphases
            velocity_source => extract_vector_field(state(i), "VelocitySource", stat)
            if (stat==0) then
               do j=1,node_count(velocity_source)
                  u_source((i-1)*node_count(velocity_source)+j)=velocity_source%val(X_, j)
                  if (ndim>1) u_source(u_nonods*nphases + (i-1)*node_count(velocity_source)+j) = velocity_source%val(Y_, j)
                  if (ndim>2) u_source(2*u_nonods*nphases + (i-1)*node_count(velocity_source)+j) = velocity_source%val(Z_, j)                  
               enddo
            else
               u_source = 0.
            endif
         enddo
      end if Conditional_VelocitySource

      if (have_gravity) then
         do i = 1, nphases - 1, 1
            delta_den = 0.
            do j = 2, node_count( density )
               delta_den = delta_den + ( den((i)*node_count(density)+j) - &
                    den((i-1)*node_count(density)+j) ) / &
                    real( ( node_count( density ) - 1 ) * max( 1, ( nphases - 1 )))
            end do
            do j = 1, u_nonods
               do k = 1, ndim
                  u_source( ( k - 1 ) * u_nonods * nphases + ( i - 1 ) * u_nonods + j  ) = &
                       u_source( ( k - 1 ) * u_nonods * nphases + ( i - 1 ) * u_nonods + j  ) + &
                       delta_den * gravity_magnitude * gravity_direction%val(k,1) * &
                       domain_length / ( 5.5 * real( totele ))
               end do
            end do
         end do
      end if

      ewrite(3,*) 'Getting PVF Source'
      do i=1,nphases
         pvf_source => extract_scalar_field(state(i), "PhaseVolumeFractionSource", stat)
         if (.not.allocated(v_source)) allocate(v_source(cv_nonods*nphases))
         if (stat==0) then
            do j=1,node_count(pvf_source)
               v_source((i-1)*node_count(pvf_source)+j)=pvf_source%val(j)
            enddo
         else
            v_source = 0.
         endif
      enddo

      ewrite(3,*) 'Getting dummy temperature field'
      if ( nscalar_fields > 2 ) then ! we might have an extra scalar_field so might need t
         ! Assuming that it's a temperature field
         Conditional_ExtraScalarField: if( have_option( "/material_phase[0]/" // &
              "scalar_field::Temperature" ))then

            Loop_Temperature: do i = 1, nphases
               scalarfield => extract_scalar_field(state(i), "Temperature")
               if (.not.allocated(t)) allocate( t( cv_nonods * nphases ))

               do j = node_count( scalarfield ), 1, -1
                  t(( i - 1 ) * node_count( scalarfield ) + j ) = scalarfield%val( j )
                  if( j == 1 ) t(( i - 1 ) * node_count( scalarfield ) + j ) = &
                       t(( i - 1 ) * node_count( scalarfield ) + j + 1 )
               enddo

               scalarfield_source => extract_scalar_field( state( i ), trim( field_name ) // &
                    "Source", stat)

               if ( .not. allocated( t_source )) allocate( t_source( cv_nonods * nphases ))
               if ( stat == 0 ) then
                  do j = 1, node_count( scalarfield_source )
                     t_source(( i - 1 ) * node_count( scalarfield_source) + j ) = &
                          scalarfield_source%val( j )
                  enddo
               else
                  t_source = 0.
               endif

               if( ( .not. allocated( wic_t_bc )) .and. ( .not. allocated( suf_t_bc ))) then
                  allocate( wic_t_bc( stotel * nphases ))
                  allocate( suf_t_bc( stotel * 1 * nphases ))
                  wic_t_bc = 0.
                  suf_t_bc = 0.
               endif

               Conditional_Temperature_BC: if( have_option( "/material_phase[" // int2str(i-1) // &
                    "]/scalar_field::Temperature/prognostic/boundary_conditions[0]/" // &
                    "type::dirichlet" )) then

                  shape_option=option_shape("/material_phase[" // int2str(i-1) // &
                       "]/scalar_field::Temperature/prognostic/boundary_conditions[0]/" // &
                       "surface_ids")

                  if( .not. allocated( Temperature_sufid_bc )) &
                       allocate( Temperature_sufid_bc( 1 : shape_option( 1 )))
                  Temperature_bc_type = 1

                  call get_option( "/material_phase[" // int2str(i-1) // "]/scalar_field::" // &
                       "Temperature/prognostic/" // &
                       "boundary_conditions[0]/surface_ids", Temperature_sufid_bc )

                  do j=1,shape_option(1)
                     wic_t_bc( Temperature_sufid_bc(1) + (i-1)*nphases ) = Temperature_bc_type
                  enddo

                  nobcs = get_boundary_condition_count( scalarfield )
                  do j = 1, nobcs
                     scalarfield_bc => extract_surface_field( scalarfield, j, "value" )
                  end do
                  do j = 1, node_count( scalarfield_bc )
                     suf_t_bc( ( i - 1 ) * stotel + j ) = scalarfield_bc%val( j )
                  end do

                  deallocate(Temperature_sufid_bc)

               endif Conditional_Temperature_BC

            enddo Loop_Temperature

         end if Conditional_ExtraScalarField

      end if

      ewrite(3,*)'temperature field:',t

      ewrite(3,*) 'Getting component source'
      allocate(comp_source(cv_nonods*nphases))
      comp_source=0.
      ! comp is stored in the order
      !   comp1 phase1
      !   comp1 phase2
      !   comp2 phase1
      !   comp2 phase2
      !   etc
      ! Component order doesn't really matter but the material_phase names
      ! should match up
      if (ncomps>0) then
         allocate(comp(cv_nonods*nphases*ncomps))
         !! Assume for now that components have been inserted in state
         !! AFTER all the phases
         do i=nstates-ncomps+1,nstates
            do j=1,option_count("/material_phase[" // int2str(i-1) //"]/scalar_field")
               componentmassfraction => extract_scalar_field(state(i), j)
               if (componentmassfraction%name(7:27) == "ComponentMassFraction") then
                  call get_option("/material_phase[" // int2str(i-1) //"]/scalar_field[" // int2str(j-1) //&
                       &"]/material_phase_name", material_phase_name)
                  do k=1,nphases
                     ewrite(3,*) 'mp_name, state_name', trim(material_phase_name), state(k)%name
                     if (trim(material_phase_name) == state(k)%name) then
                        ewrite(3,*) 'here'
                        do l=1,node_count(componentmassfraction)
                           comp( ((i-(1+nphases))*nphases+(k-1))*node_count(componentmassfraction)+l ) = componentmassfraction%val(l)
                        enddo
                     endif
                  enddo
               endif
            enddo
         enddo
      endif

      ewrite(3,*) 'Getting capillary pressure options and absorptions'
      if (ncapil_pres_coef>0) then
         allocate(capil_pres_coef(ncapil_pres_coef,nphases,nphases))
         capil_pres_coef=0.
      endif
      ! Not sure what this one does yet
      allocate(u_abs_stab(cv_nloc*totele, ndim*nphases, ndim*nphases))
      u_abs_stab=0.
      ! or this one
      allocate(u_absorb(cv_nloc*totele, ndim*nphases, ndim*nphases))
      u_absorb=0.
      allocate(comp_absorb(cv_nonods, nphases, nphases))
      comp_absorb=0.
      allocate( t_absorb( cv_nonods, nphases, nphases ))
      t_absorb=0.
      allocate( v_absorb( cv_nonods, nphases, nphases ))
      v_absorb=0.
      allocate( k_comp( ncomps, nphases, nphases ))
      k_comp=0.
      if (KComp_Sigmoid) then
         do i=1, ncomps
            call get_option('material_phase['// int2str(i+nphases-1) //']/is_multiphase_component/' // &
                 'KComp_Sigmoid/k_comp', k_comp(i, 1, 1))
            ewrite(3,*) 'i, kcomp', i, k_comp(i, 1, 1)
            k_comp(i, 1:nphases, 1:nphases) = k_comp(i, 1, 1)
         end do
      end if

      allocate( comp_diffusion( cv_nloc*totele, ndim, ndim, nphases ))
      comp_diffusion=0.
      allocate( comp_diff_coef( ncomps, ncomp_diff_coef, nphases ))
      comp_diff_coef=0.
      allocate( cv_one( nphases * cv_nonods ))
      cv_one = 0.

      ewrite(3,*) "Leaving copy_outof_state"

    end subroutine copy_outof_state



    subroutine copy_into_state(state, saturations, proto_pressure, nphase, cv_ndgln)
      
      !!< Copy prototype saturations and pressure into fluidity state array for output
      
      type(state_type), dimension(:), intent(inout) :: state
      real, dimension(:), intent(in) :: saturations
      real, dimension(:), intent(in) :: proto_pressure
      integer, intent(in) :: nphase
      integer, dimension(:), intent(in) :: cv_ndgln
      
      ! local variables        
      integer :: stat
      integer :: i,j,p
      integer, dimension(:), pointer :: element_nodes
      type(scalar_field), pointer :: phasevolumefraction
      type(scalar_field), pointer :: pressure

      ewrite(3,*) "In copy_into_state"
      
      ! set volume fraction fields for each phase into state
      
      assert(size(state) >= nphase)
      
      phase_loop: do p = 1,nphase
         
         phasevolumefraction => extract_scalar_field(state(p), "PhaseVolumeFraction", stat=stat)
         
         if (stat /= 0) then 
            
            ewrite(1,*) 'Issue in prototype interface for phase ',p
            
            FLAbort('Failed to extract phase volume fraction from state in copy_into_state')
         
         end if
         
         volf_ele_loop: do i = 1,element_count(phasevolumefraction)
            
            element_nodes => ele_nodes(phasevolumefraction,i)
            
            volf_node_loop: do j = 1,size(element_nodes)
               
               ! this is hard wired for quadratic elements 1d with regard to prototype array
               call set(phasevolumefraction, &
                        element_nodes(j), &
!!!                        saturations((2*(i-1)+j) + (p-1)*node_count(phasevolumefraction)))
                        saturations((cv_ndgln((i-1)*size(element_nodes)+j)) + (p-1)*node_count(phasevolumefraction)))
            
            end do volf_node_loop
            
         end do volf_ele_loop
      
      end do phase_loop
         
      pressure => extract_scalar_field(state(1), "Pressure")    
      
      press_ele_loop: do i = 1,ele_count(pressure)
        
        element_nodes => ele_nodes(pressure,i)
        
        press_node_loop: do j = 1,size(element_nodes)
          
          ! this is hard wired for quadratic elements 1d with regard to prototype array
          call set(pressure, &
                   element_nodes(j), &
                   proto_pressure(2*(i-1)+j))
        
        end do press_node_loop
      
      end do press_ele_loop
    
      ewrite(3,*) "Leaving copy_into_state"
  
    end subroutine copy_into_state

end module copy_outof_into_state
