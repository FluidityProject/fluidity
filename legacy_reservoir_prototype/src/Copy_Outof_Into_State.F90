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
  use diagnostic_fields_wrapper
  use global_parameters, only: option_path_len
  use diagnostic_fields_new, only : &
    & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
    & check_diagnostic_dependencies
  implicit none
  
  private
  
  public :: copy_outof_state, copy_into_state
  
  contains
  
    subroutine copy_outof_state(state, dt, current_time, finish_time, &
         nonlinear_iterations, nonlinear_iteration_tolerance)

      !! New variables

      type(state_type), dimension(:), pointer :: state
      type(mesh_type) :: cmesh, vmesh, pmesh !! coordinate, velocity and pressure meshes

      integer :: nonlinear_iterations  !! equal to nits in prototype code

      real :: dt, current_time, finish_time
      real :: nonlinear_iteration_tolerance

      !! temporary variables only needed for interfacing purposes

      type(vector_field), pointer :: positions

      integer :: i, nscalar_fields

      real :: coord_min, coord_max

      !! Variables needed by the prototype code
      !! and therefore needing to be pulled out of state or
      !! derived from information in state:

      ! Scalars (from read_scalar())
      integer :: problem, nphase, ncomp, totele, ndim, nlev, &
           u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, &
           cv_snloc, u_snloc, p_snloc, x_snloc, stotel, &
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

      integer :: U_BC_Type, P_BC_Type, SufID_BC_U, SufID_BC_P
      real :: Suf_BC_U



      !! Finish declaration of variables needed from user input

      ewrite(3,*) 'In copy_outof_state'

      !! Here are all of the needed items
      !! I suggest we work through them one by one
      !! according to Alex and Brendan's markup of the input file
      !! Shuffle the order as necessary, but it would be worth keeping
      !! related items together as much as possible. JHS

      ! problem = 
      nphase = option_count("/material_phase")
      ewrite(3,*) ' nphase:', nphase
      !! Assume there are the same number of components in each phase (will need to check this eventually)
      nscalar_fields = option_count("/material_phase[0]/scalar_field")
      ncomp=0
      do i=1,nscalar_fields
         if (have_option("/material_phase[0]/scalar_field["//int2str(i-1)//"]/prognostic/is_multiphase_component")) then
            ncomp=ncomp+1
         end if
      end do
      ewrite(3,*) ' ncomp:',ncomp

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
      ! u_nloc = 
      ! xu_nloc = 
      ! cv_nloc = 
      ! x_nloc = 
      ! p_nloc = 
      ! cv_snloc = 
      ! u_snloc = 
      ! p_snloc = 
      ! x_snloc = 
      ! stotel = 

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
      !! Sufficient to set this to 1 for now? It is 1 in all test cases
      mat_ele_type = 1
      !! Sufficient to set this to 2 for now? It is 2 in all test cases
      cv_ele_type = 2
      !! These aren't used 
      cv_sele_type = 0
      u_sele_type = 0

      !! Time options
      ewrite(3,*) ' Getting time options'
      call get_option('/io/max_dump_file_count', ntime, default=500)
      if (have_option('/io/dump_period_in_timesteps/constant')) then
         call get_option('/io/dump_period_in_timesteps/constant',ntime_dump)
      else
         ! This way might be inaccurate due to rounding errors:
         call get_option('/io/dump_period',ntime_dump)
         ntime_dump = ntime_dump/dt
      end if

      ewrite(3,*) ' Getting iteration info'
      nits = nonlinear_iterations
      ! This one is only for compositional problems
      call get_option('/material_phase::phase1/scalar_field::component1/' // &
           'prognostic/temporal_discretisation/control_volumes/number_advection_iterations',&
           nits_internal, default=1)
      ! ndpset = 
      ! noit_dim = 
      ! nits_flux_lim_volfra = 
      ! nits_flux_lim_comp = 
      ! t_disopt = 
      ! u_disopt = 
      ! v_disopt = 
      ! t_dg_vel_int_opt = 
      ! u_dg_vel_int_opt = 
      ! v_dg_vel_int_opt = 
      ! w_dg_vel_int_opt = 

      ewrite(3,*) ' Getting capillary pressure options'
      if (have_option('/porous_media/multiphase_parameters/cp_A')) then
         ! 1 is the only option available at the moment...
         capil_pres_opt = 1
         ! and even this one doesn't depend on any coefficients!
         ncapil_pres_coef = 0
      else
         capil_pres_opt = 0
      end if

      !! These are not currently used in any of the code
      comp_diffusion_opt = 0
      ncomp_diff_coef = 0

      patmos = 0.0
      p_ini = 0.0

      t_ini = 0.0
      t_beta = 0

      call get_option('/material_phase::phase1/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'spatial_discretisation/conservative_advection', v_beta)
      ! t_theta = 
      call get_option('/material_phase::phase1/scalar_field::PhaseVolumeFraction/prognostic/' // &
           'temporal_discretisation/theta', v_theta)
      call get_option('/material_phase::phase1/vector_field::Velocity/prognostic/' // &
           'temporal_discretisation/theta', u_theta)

      ! This is a strictly 1d property, so will have a suitably 1d method for finding it in state!
      coord_min=1.0e9
      coord_max=-1.0e-9
      do i=1,node_count(positions)
         ewrite(3,*) 'positions:', positions%val(X_,i)
         coord_min=min(coord_min,positions%val(X_,i))
         coord_max=max(coord_max,positions%val(X_,i))
      end do
      domain_length = coord_max - coord_min
      ewrite(3,*) 'domain_length:',domain_length

      ! lump_eqns = 
      ! volfra_use_theta_flux = 
      ! volfra_get_theta_flux = 
      ! comp_use_theta_flux = 
      ! comp_get_theta_flux


      ! Others (from read_all())

      !!
      !! All options bellow should be replaced by PETSc options soon, i.e., still at the second part of
      !! of the first stage of the integration.  So let's keep as it is for the moment and remove them
      !! as soon as we have all PETSc data structure enabled.

      volfra_relax_number_iterations = 100
      scalar_relax_number_iterations = 100
      global_relax_number_iterations = 100
      velocity_relax_number_iterations = 100
      pressure_relax_number_iterations = 4000
      mass_matrix_relax_number_iterations = 100 

      volfra_error = 1.e-5
      scalar_error = 1.e-5
      global_error = 1.e-5
      velocity_error = 1.e-5
      pressure_error = 1.e-3
      mass_matrix_error = 1.-5

      volfra_relax = 1.
      scalar_relax = 1.
      global_relax = 1.
      velocity_relax = 1.
      pressure_relax = 1.
      mass_matrix_relax = 1. 

      volfra_relax_diag = 0.
      scalar_relax_diag = 0.
      global_relax_diag = 0.
      velocity_relax_diag = 0.
      pressure_relax_diag = 0.
      mass_matrix_relax_diag = 0. 

      volfra_relax_row = 1.
      scalar_relax_row = 1.
      global_relax_row = 1.
      velocity_relax_row = 1.
      pressure_relax_row = 1.
      mass_matrix_relax_row = 1. 

      ! IN/DG_ELE_UPWIND are options for optimisisation of upwinding across faces in the overlapping
      ! formulation. The data structure and options for this formulation need to be added later. 
      in_ele_upwind = 3
      dg_ele_upwind = 3

      ! Calculating Mobility
      Viscosity = 0.
      Viscosity_Ph1 = 0.
      Viscosity_Ph2 = 0.
      call get_option( '/material_phase::phase1/vector_field::Velocity/prognostic/' // &
           'tensor_field::Viscosity/prescribed/value::WholeMesh/' // &
           'isotropic/constant', Viscosity_Ph1 )
      call get_option( '/material_phase::phase2/vector_field::Velocity/prognostic/' // &
           'tensor_field::Viscosity/prescribed/value::WholeMesh/' // &
           'isotropic/constant', Viscosity_Ph2 )

      Viscosity( 1 : cv_nonods ) = Viscosity_Ph1
      Viscosity( cv_nonods + 1 : cv_nonods * nphase ) = Viscosity_Ph2

      ! Maybe should be worthy to add Mobility to schema and Viscosity_Ph2 become a prognostic field
      Mobility = Viscosity_Ph1 / Viscosity_Ph2

!!!
!!! Options bellow are for the multi-component flow model, still needed to be added into the schema
!!! 
      alpha_beta = 1.
      KComp_Sigmoid = .true. 
      Comp_Sum2One = .false.

!!!
!!! Porosity and Permeability: it may be neecessary to change the permeability as it
!!! is defined in the PC as a tensor with dimension ( totele, ndim, ndim )in
!!!
      call get_option( '/porous_media/scalar_field::Porosity/prescribed/' // &
           'value::WholeMesh/constant', volfra_pore )

      call get_option( '/porous_media/scalar_field::Permeability/prescribed/' // &
           'value::WholeMesh/constant', perm )

!!!
!!! WIC_X_BC (in which X = D, U, V, W, P, T, COMP and VOL) controls the boundary conditions
!!! type applied. == 1 (Dirichlet), = 2 (Robin), = 3 (Newman) 

      Conditional_Velocity_BC: if( have_option( '/material_phase[0]/vector_field::Velocity/prognostic/' // &
           'boundary_conditions[0]/type::dirichlet' )) then
         U_BC_Type = 1
         call get_option( '/material_phase[0]/vector_field::Velocity/' // &
              'prognostic/boundary_conditions::/surface_ids', SufID_BC_U )
         call get_option( '/material_phase[0]/vector_field::Velocity/' // &
              'prognostic/boundary_conditions::/type::dirichlet/' // &
              'align_bc_with_cartesian/x_component/constant', Suf_BC_U )
      elseif( have_option( '/material_phase[0]/vector_field::Velocity/prognostic/' // &
           'boundary_conditions[0]/type::neumann' )) then
         U_BC_Type = 3
         call get_option( '/material_phase[0]/vector_field::Velocity/' // &
              'prognostic/boundary_conditions::/surface_ids', SufID_BC_U )
         call get_option( '/material_phase[0]/vector_field::Velocity/' // &
              'prognostic/boundary_conditions::/type::neumann/' // &
              'align_bc_with_cartesian/x_component/constant', Suf_BC_U )
      endif Conditional_Velocity_BC

      Conditional_Pressure_BC:if( have_option( '/material_phase[0]/scalar_field::Pressure/prognostic/' // &
           'boundary_conditions[0]/type::dirichlet' )) then
         P_BC_Type = 1
         call get_option( '/material_phase[0]/scalar_field::Pressure/prognostic/' // &
              'boundary_conditions[0]/surface_ids', SufID_BC_P )
         call get_option( '/material_phase[0]/scalar_field::Pressure/prognostic/' // &
              'boundary_conditions[0]/type::dirichlet/constant', Suf_BC_P )
      endif Conditional_Pressure_BC


      ! wic_vol_bc
      ! wic_d_bc
      ! wic_u_bc
      ! wic_p_bc
      ! wic_t_bc
      ! wic_comp_bc
      ! uabs_option
      ! eos_option
      ! cp_option

      ! suf_vol_bc
      ! suf_d_bc
      ! suf_cpd_bc
      ! suf_t_bc
      ! suf_p_bc
      ! suf_u_bc
      ! suf_v_bc
      ! suf_w_bc
      ! suf_one_bc
      ! suf_comp_bc
      ! suf_u_bc_rob1
      ! suf_u_bc_rob2
      ! suf_v_bc_rob1
      ! suf_v_bc_rob2
      ! suf_w_bc_rob1
      ! suf_w_bc_rob2
      ! suf_t_bc_rob1
      ! suf_t_bc_rob2
      ! suf_comp_bc_rob1
      ! suf_comp_bc_rob2
      ! x
      ! y
      ! z
      ! xu
      ! yu
      ! zu
      ! nu
      ! nv
      ! nw
      ! ug
      ! vg
      ! wg
      ! u_source
      ! t_source
      ! v_source
      ! comp_source
      ! u
      ! v
      ! w
      ! den
      ! satura
      ! volfra
      ! comp
      ! t
      ! p
      ! cv_p
      ! cv_one

      ! uabs_coefs
      ! eos_coefs
      ! cp_coefs

      ! comp_diff_coef
      ! capil_pres_coef
      ! u_abs_stab
      ! u_absorb
      ! comp_absorb
      ! t_absorb
      ! v_absorb
      ! perm
      ! K_Comp

      ! comp_diffusion

      ewrite(3,*) 'Leaving copy_outof_state'

    end subroutine copy_outof_state
  
  subroutine copy_into_state()
  
    ewrite(3,*) 'In copy_into_state'
    
    
    ewrite(3,*) 'Leaving copy_into_state'
  
  end subroutine copy_into_state

end module copy_outof_into_state
