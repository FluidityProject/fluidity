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
  use diagnostic_fields_wrapper_new
!  use diagnostic_fields_new, only : &
!    & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
!    & check_diagnostic_dependencies
  implicit none
  
  private
  
  public :: copy_outof_state, copy_into_state
  
  contains
  
    subroutine copy_outof_state(state, dt, current_time, finish_time, &
                              nonlinear_iterations, nonlinear_iteration_tolerance)
  
  !! New variables
  
      type(state_type), dimension(:), pointer :: state
      type(mesh_type) :: cmesh, vmesh, pmesh !! coordinate, velocity and pressure meshes
    
      type(scalar_field), pointer :: porosity, density, pressure, &
                                     phasevolumefraction

      type(vector_field), pointer :: velocity

      type(tensor_field), pointer :: viscosity_ph1, viscosity_ph2
      
      integer :: nonlinear_iterations  !! equal to nits in prototype code
      
      real :: dt, current_time, finish_time
      real :: nonlinear_iteration_tolerance
  
  !! temporary variables only needed for interfacing purposes
  
      type(vector_field), pointer :: positions
  
      integer :: i, j, k, nscalar_fields, cv_nonods, &
                 U_BC_Type, P_BC_Type, SufID_BC_U, SufID_BC_P
      real :: Suf_BC_U, suf_bc_p, &
              coord_min, coord_max, &
              permeability!, viscosity_ph1, viscosity_ph2
              
      integer, dimension(:), allocatable :: elenodes
              
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

      integer :: Velocity_BC_Type, Pressure_BC_Type, density_bc_type, shape_option(2)
      real :: Velocity_Suf_BC_U, Velocity_Suf_BC_V, Velocity_Suf_BC_W, &
              Pressure_Suf_BC, density_suf_bc
              
      integer, dimension(:), allocatable :: Velocity_SufID_BC, Pressure_SufID_BC, density_sufid_bc



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
      u_nloc = vmesh%shape%loc ! 6, but only when we've got the right options in the schema
      xu_nloc = cmesh%shape%loc
      cv_nloc = pmesh%shape%loc
      x_nloc = 3*ndim
      p_nloc = pmesh%shape%loc ! 3
    ! cv_snloc = 1
    ! u_snloc = -1
    ! p_snloc = 1
    ! x_snloc = 1
      stotel = surface_element_count(cmesh)
    
      ewrite(3,*) 'u_nloc, x_nloc, p_nloc', u_nloc, x_nloc, p_nloc
      ewrite(3,*) 'stotel:', stotel
    
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
        ntime_dump = int(ntime_dump/dt)
      end if
    
      ewrite(3,*) ' Getting iteration info'
      nits = nonlinear_iterations
    ! This one is only for compositional problems
      call get_option('/material_phase[0]/scalar_field::component1/&
                      &prognostic/temporal_discretisation/control_volumes/number_advection_iterations',&
                      &nits_internal, default=1)
                    
      ndpset = 0
      noit_dim = 5
    
    !! These two are always equal to 3 except for the 2phase, 3comp, non-equilibrium test, where they are 1
    !! Change in schema ??
    ! nits_flux_lim_volfra = 
    ! nits_flux_lim_comp = 


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

      t_disopt = 1 ! I don't know what the theta=non-linear means and I think t field isn't used at the moment anyway
    
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
    ! u_dg_vel_int_opt =  ! Not used
      v_dg_vel_int_opt = 4
    ! w_dg_vel_int_opt =  ! Not used
    
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

      call get_option('/material_phase[0]/scalar_field::Pressure/prognostic/&
                      &atmospheric_pressure', patmos, default=0.0)
      call get_option('/material_phase[0]/scalar_field::Pressure/prognostic/&
                      &initial_condition::WholeMesh/constant',p_ini, default=0.0)
    
      t_ini = 0.0
      if (nscalar_fields>2) then ! we might have an extra scalar_field so might need t
        call get_option('/material_phase[0]/scalar_field[2]/prognostic/' // &
                        'spatial_discretisation/conservative_advection', t_beta)
      end if

      call get_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
                      'spatial_discretisation/conservative_advection', v_beta)
    
      if (nscalar_fields>2) then ! we might have an extra scalar_field so might need t
        call get_option('/material_phase[0]/scalar_field[2]/prognostic/' // &
                        'temporal_discretisation/theta', t_theta)
      end if
      call get_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/' // &
                      'temporal_discretisation/theta', v_theta)
      call get_option('/material_phase[0]/vector_field::Velocity/prognostic/' // &
                      'temporal_discretisation/theta', u_theta)
    
    ! This is a strictly 1d property, so will have a suitably 1d method for finding it in state!
      coord_min=1.0e9
      coord_max=-1.0e-9
      do i=1,node_count(positions)
      !ewrite(3,*) 'positions:', positions%val(X_,i)
        coord_min=min(coord_min,positions%val(X_,i))
        coord_max=max(coord_max,positions%val(X_,i))
      end do
      domain_length = coord_max - coord_min
      ewrite(3,*) 'domain_length:',domain_length

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
         call get_option( '/material_phase[0]/scalar_field::Pressure/prognostic/solver/relative_error', &
              pressure_error, default = 1.e-3 )
         pressure_relax = 1.
         pressure_relax_diag = 0.
         pressure_relax_row = 1.
      endif Conditional_Pressure_Solver

      Conditional_Velocity_Solver: if( have_option( '/material_phase[0]/vector_field::Velocity' )) then
         call get_option( '/material_phase[0]/vector_field::Velocity/prognostic/solver/max_iterations', &
              velocity_relax_number_iterations, default = 100 )
         call get_option( '/material_phase[0]/vector_field::Velocity/prognostic/solver/relative_error', &
              velocity_error, default = 1.e-5 ) 
         velocity_relax = 1.
         velocity_relax_diag = 0.
         velocity_relax_row = 1.
      end if Conditional_Velocity_Solver
      
      ewrite(3,*) 'Got first lot of solver options'

      scalar_relax_number_iterations = 100
      global_relax_number_iterations = 100
      mass_matrix_relax_number_iterations = 100 

      scalar_error = 1.e-5
      global_error = 1.e-5
      mass_matrix_error = 1.-5

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
      ewrite(3,*) 'going to get viscosities'

      viscosity_ph1 => extract_tensor_field(state(1), "Viscosity")
      viscosity_ph2 => extract_tensor_field(state(2), "Viscosity")
      
      ewrite(3,*) 'Got viscosities'

      ! Maybe should be worth adding Mobility to schema and Viscosity_Ph2 become a diagnostic field
      if (have_option('/material_phase[1]/vector_field::Velocity/prognostic/&
                      &tensor_field::Viscosity/prescribed/value::WholeMesh/&
                      &isotropic')) then
        Mobility = Viscosity_Ph1%val(1,1,1) / Viscosity_Ph2%val(1,1,1)
      endif
      
      ewrite(3,*) 'sorted mobility'

!!!
!!! Options below are for the multi-component flow model, still needed to be added into the schema
!!! 
      alpha_beta = 1.
      KComp_Sigmoid = .true. 
      Comp_Sum2One = .false.

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
      ewrite(3,*) 'Got porosity'

      if (have_option('/porous_media/scalar_field::Permeability')) then
        call get_option( '/porous_media/scalar_field::Permeability/prescribed/&
                         &value::WholeMesh/constant', permeability )
        allocate(perm(totele, ndim, ndim))
        do i=1,totele
          do j=1,ndim
            do k=1,ndim
               perm(i,j,k)=permeability
            end do
          end do
        end do
      elseif (have_option('/porous_media/tensor_field::Permeability')) then
        FLAbort('Have not coded up tensor permeability yet! Try scalar instead')
      endif
      
      ewrite(3,*) 'Got permeability'

!!!
!!! WIC_X_BC (in which X = D, U, V, W, P, T, COMP and VOL) controls the boundary conditions
!!! type applied. == 1 (Dirichlet), = 2 (Robin), = 3 (Newman) 

      Conditional_Velocity_BC_U: if( have_option( '/material_phase[0]/vector_field::Velocity/prognostic/' // &
           'boundary_conditions[0]/type::dirichlet' )) then

         shape_option=option_shape('/material_phase[0]/vector_field::Velocity/&
                                   &prognostic/boundary_conditions[0]/surface_ids')
         allocate(velocity_sufid_bc(1:shape_option(1)))
         ewrite(3,*) 'allocated vel suf id'

         Velocity_BC_Type = 1
         call get_option( '/material_phase[0]/vector_field::Velocity/' // &
              'prognostic/boundary_conditions[0]/surface_ids', Velocity_SufID_BC )
         ewrite(3,*) 'velocity_sufid_bc', velocity_sufid_bc
         if( have_option( '/material_phase[0]/vector_field::Velocity/' // &
              'prognostic/boundary_conditions[0]/type::dirichlet/' // &
              'align_bc_with_cartesian/x_component' )) then
            call get_option( '/material_phase[0]/vector_field::Velocity/' // &
                 'prognostic/boundary_conditions[0]/type::dirichlet/' // &
                 'align_bc_with_cartesian/x_component/constant', Velocity_Suf_BC_U )
         endif
         if( have_option( '/material_phase[0]/vector_field::Velocity/' // &
              'prognostic/boundary_conditions[0]/type::dirichlet/' // &
              'align_bc_with_cartesian/y_component' )) then
            call get_option( '/material_phase[0]/vector_field::Velocity/' // &
                 'prognostic/boundary_conditions[0]/type::dirichlet/' // &
                 'align_bc_with_cartesian/y_component/constant', Velocity_Suf_BC_V )
         endif
         if( have_option( '/material_phase[0]/vector_field::Velocity/' // &
              'prognostic/boundary_conditions[0]/type::dirichlet/' // &
              'align_bc_with_cartesian/z_component' )) then
            call get_option( '/material_phase[0]/vector_field::Velocity/' // &
                 'prognostic/boundary_conditions[0]/type::dirichlet/' // &
                 'align_bc_with_cartesian/z_component/constant', Velocity_Suf_BC_W )
         endif
         ewrite(3,*) 'vel suf bc u', velocity_suf_bc_u

      elseif( have_option( '/material_phase[0]/vector_field::Velocity/prognostic/' // &
           'boundary_conditions[0]/type::neumann' )) then
         Velocity_BC_Type = 3
         call get_option( '/material_phase[0]/vector_field::Velocity/' // &
              'prognostic/boundary_conditions[0]/surface_ids', Velocity_SufID_BC )
         if( have_option( '/material_phase[0]/vector_field::Velocity/' // &
              'prognostic/boundary_conditions[0]/type::neumann/' // &
              'align_bc_with_cartesian/x_component' )) then
            call get_option( '/material_phase[0]/vector_field::Velocity/' // &
                 'prognostic/boundary_conditions[0]/type::neumann/' // &
                 'align_bc_with_cartesian/x_component/constant', Velocity_Suf_BC_U )
         endif
         if( have_option( '/material_phase[0]/vector_field::Velocity/' // &
              'prognostic/boundary_conditions[0]/type::neumann/' // &
              'align_bc_with_cartesian/y_component' )) then
            call get_option( '/material_phase[0]/vector_field::Velocity/' // &
                 'prognostic/boundary_conditions[0]/type::neumann/' // &
                 'align_bc_with_cartesian/y_component/constant', Velocity_Suf_BC_V )
         endif
         if( have_option( '/material_phase[0]/vector_field::Velocity/' // &
              'prognostic/boundary_conditions[0]/type::neumann/' // &
              'align_bc_with_cartesian/y_component' )) then
            call get_option( '/material_phase[0]/vector_field::Velocity/' // &
                 'prognostic/boundary_conditions[0]/type::neumann/' // &
                 'align_bc_with_cartesian/z_component/constant', Velocity_Suf_BC_W )
         endif

      endif Conditional_Velocity_BC_U
      
      allocate( wic_u_bc( stotel * nphase ))
      allocate( suf_u_bc( stotel * 3 * nphase ))
      allocate( suf_v_bc( stotel * 3 * nphase ))
      allocate( suf_w_bc( stotel * 3 * nphase ))
      wic_u_bc = 0
      suf_u_bc = 0.
      suf_v_bc = 0.
      suf_w_bc = 0.
      do i=1,shape_option(1)
        wic_u_bc( Velocity_SufID_BC(i) ) = Velocity_BC_Type
        suf_u_bc( Velocity_SufID_BC(i) ) = Velocity_Suf_BC_U
        suf_v_bc( Velocity_SufID_BC(i) ) = Velocity_Suf_BC_V
        suf_w_bc( Velocity_SufID_BC(i) ) = Velocity_Suf_BC_W
      enddo

      deallocate(velocity_sufid_bc)
      ewrite(3,*) 'Done with velocity boundary conditions'

      Conditional_Pressure_BC: if( have_option( '/material_phase[0]/scalar_field::Pressure/prognostic/' // &
                                                'boundary_conditions[0]/type::dirichlet' )) then

        shape_option=option_shape('/material_phase[0]/scalar_field::Pressure/&
                                  &prognostic/boundary_conditions[0]/surface_ids')
        allocate(pressure_sufid_bc(1:shape_option(1)))
        Pressure_BC_Type = 1
        call get_option( '/material_phase[0]/scalar_field::Pressure/prognostic/' // &
                'boundary_conditions[0]/surface_ids', Pressure_SufID_BC )
        call get_option( '/material_phase[0]/scalar_field::Pressure/prognostic/' // &
                'boundary_conditions[0]/type::dirichlet/constant', Pressure_Suf_BC )

        allocate( wic_p_bc( stotel * nphase ))
        allocate( suf_p_bc( stotel * 1 * nphase ))
        wic_p_bc = 0
        suf_p_bc = 0.
        do i=1,shape_option(1)
          wic_p_bc( Pressure_SufID_BC(1) ) = Pressure_BC_Type
          suf_p_bc( Pressure_SufID_BC(1) ) = Pressure_Suf_BC
        enddo

        deallocate(pressure_sufid_bc)

      endif Conditional_Pressure_BC
      
      Conditional_Density_BC: if( have_option( '/material_phase[0]/scalar_field::Density/prognostic/' // &
           'boundary_conditions[0]/type::dirichlet' )) then

        shape_option=option_shape('/material_phase[0]/scalar_field::Density/&
                                  &prognostic/boundary_conditions[0]/surface_ids')
        allocate(density_sufid_bc(1:shape_option(1)))
        Density_BC_Type = 1
        call get_option( '/material_phase[0]/scalar_field::Density/prognostic/' // &
             'boundary_conditions[0]/surface_ids', Density_SufID_BC )
        call get_option( '/material_phase[0]/scalar_field::Density/prognostic/' // &
             'boundary_conditions[0]/type::dirichlet/constant', Density_Suf_BC )

        allocate( wic_d_bc( stotel * nphase ))
        allocate( suf_d_bc( stotel * 1 * nphase ))
        wic_d_bc = 0
        suf_d_bc = 0.
        do i=1,shape_option(1)
          wic_d_bc( Density_SufID_BC(1) ) = Density_BC_Type
          suf_d_bc( Density_SufID_BC(1) ) = Density_Suf_BC
        enddo
        
        deallocate(density_sufid_bc)
  
      endif Conditional_Density_BC
      
      ewrite(3,*) 'Done with boundary conditions'

!!!
!!! Robin Boundary conditions need to be added at a later stage
!!!
      allocate( suf_u_bc_rob1( stotel * 3 * nphase ))
      allocate( suf_v_bc_rob1( stotel * 3 * nphase ))
      allocate( suf_w_bc_rob1( stotel * 3 * nphase ))
      allocate( suf_u_bc_rob2( stotel * 3 * nphase ))
      allocate( suf_v_bc_rob2( stotel * 3 * nphase ))
      allocate( suf_w_bc_rob2( stotel * 3 * nphase ))
      allocate( suf_t_bc_rob1( stotel * nphase ))
      allocate( suf_t_bc_rob2( stotel * nphase ))
      allocate( suf_vol_bc_rob1( stotel * nphase ))
      allocate( suf_vol_bc_rob2( stotel * nphase ))
      allocate( suf_comp_bc_rob1( stotel * nphase ))
      allocate( suf_comp_bc_rob2( stotel * nphase ))
      suf_u_bc_rob1 = 0.
      suf_u_bc_rob2 = 0.
      suf_v_bc_rob1 = 0.
      suf_v_bc_rob2 = 0.
      suf_w_bc_rob1 = 0.
      suf_w_bc_rob2 = 0.
      suf_t_bc_rob1 = 0.
      suf_t_bc_rob2 = 0.
      suf_comp_bc_rob1 = 0.
      suf_comp_bc_rob2 = 0.

!!!==============================================!!!
!!! wic_comp_bc ==> These still need to be done  !!!
!!! suf_comp_bc                                  !!!
!!!==============================================!!!

!!!
!!! Initial conditions for all fields
!!!
    ! Density and pressure might need re-ordering because of the
    ! peculiarity of the quadratic element setup
    ! see copy_into_state() below
      ewrite(3,*) 'going to get density...'
      do i=1,nphase
        density => extract_scalar_field(state(i), "Density")
!       ewrite(3,*) 'size of density', node_count(density)
        if (.not.allocated(den)) allocate(den(nphase*node_count(density)))
        do j=1,node_count(density)
           den((i-1)*node_count(density)+j)=density%val(j)
        enddo
      enddo
      
      ewrite(3,*) 'pressure...'
      pressure => extract_scalar_field(state(1), "Pressure")
      allocate(p(node_count(pressure)))
      do i=1,node_count(pressure)
         p(i)=pressure%val(i)
      enddo
      
!!! Need to add cv_p (the cv representation of pressure field used in the interpolation
!!! (overlapping) formulation
      ewrite(3,*) 'phasevolumefraction...'
      do i=1,nphase
        phasevolumefraction => extract_scalar_field(state(i), "PhaseVolumeFraction")
        if (.not.allocated(satura)) allocate(satura(nphase*node_count(phasevolumefraction)))
        do j=1,node_count(phasevolumefraction)
          satura((i-1)*node_count(phasevolumefraction)+j)=phasevolumefraction%val(j)
        enddo
      enddo
      
      ewrite(3,*) 'velocity...'
      do i=1,nphase
        velocity => extract_vector_field(state(i), "Velocity")
        if (.not. allocated(u)) then
           allocate(u(nphase*node_count(velocity)))
           allocate(v(nphase*node_count(velocity)))
           allocate(w(nphase*node_count(velocity)))
        endif
        u=0.
        v=0.
        w=0.
        do j=1,node_count(velocity)
          u((i-1)*node_count(velocity)+j)=velocity%val(X_, j)
          if (ndim>1) v((i-1)*node_count(velocity)+j)=velocity%val(Y_, j)
          if (ndim>2) w((i-1)*node_count(velocity)+j)=velocity%val(Z_, j)
        enddo
      enddo

!!===================================================!!
!! How are the components of the velocity defined ?? !!
!! u, v, w (t=0) = Velocity ?                        !!
!!===================================================!!

      allocate(uabs_option(nphase))
      allocate(eos_option(nphase))
      allocate(cp_option(nphase))
      do i=1,nphase
        call get_option('/material_phase[' // int2str(i-1) // ']/multiphase_options/relperm_option',uabs_option(i))
        if (have_option('/material_phase[' // int2str(i-1) // ']/equation_of_state/incompressible/linear')) then
          eos_option(i) = 2
        elseif (have_option('/material_phase[' // int2str(i-1) // ']/equation_of_state/compressible/stiffened_gas')) then
          eos_option(i) = 1
        else
          FLAbort('Unknown EoS option for phase '// int2str(i))
        endif
        cp_option(i) = 1
      enddo
      
      !! uabs_coefs is currently only used in rel perm options which aren't
      !! selected in any of the test cases, ie the 'standard polynomial
      !! representation of relative permeability'
      allocate(uabs_coefs(nphase, nuabs_coefs))
      allocate(eos_coefs(nphase, ncoef))
      !! Capillary pressure isn't used at all at the moment
      allocate(cp_coefs(nphase, nphase))
      do i=1,nphase
        uabs_coefs(i,1) = 1.
        if (eos_option(i)==2) then
          if (have_option('/material_phase[' // int2str(i-1) // ']/equation_of_state/incompressible/linear/all_equal')) then
             call get_option('/material_phase[' // int2str(i-1) // ']/equation_of_state/incompressible/linear/all_equal', eos_coefs(i,1:ncoef))
          else
             call get_option('/material_phase[' // int2str(i-1) // ']/equation_of_state/incompressible/linear/specify_all', eos_coefs(i, :))
          endif
        else
          call get_option('/material_phase[' // int2str(i-1) // ']/equation_of_state/compressible/stiffened_gas/eos_option1', eos_coefs(i, 1))
          call get_option('/material_phase[' // int2str(i-1) // ']/equation_of_state/compressible/stiffened_gas/eos_option2', eos_coefs(i, 2))
          eos_coefs(i, 3:ncoef) = 0.
        endif
      enddo
      cp_coefs = 1.


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
      ! comp
      ! t

      ! comp_diff_coef
      ! capil_pres_coef
      ! u_abs_stab
      ! u_absorb
      ! comp_absorb
      ! t_absorb
      ! v_absorb
      ! K_Comp

      ! comp_diffusion

      ewrite(3,*) 'Leaving copy_outof_state'

    end subroutine copy_outof_state
  
    subroutine copy_into_state(state, saturations, proto_pressure)
  
      type(state_type), dimension(:), pointer :: state
      type(scalar_field), pointer :: galerkinprojection, phasevolumefraction, pressure
      type(vector_field), pointer :: positions

      integer :: i,j
      integer, dimension(:), allocatable :: elenodes
      real, dimension(:), intent(in) :: saturations, proto_pressure

      ewrite(3,*) 'In copy_into_state'

      ewrite(3,*) 'size of satura', size(saturations)
      ewrite(3,*) 'saturations:', saturations ! This is cv_nonods * nphase long, and in order

      positions => extract_vector_field(state, "Coordinate")
    
      ! The plan is to copy the first half of saturations into PhaseVolumeFraction:
      phasevolumefraction => extract_scalar_field(state(1), "PhaseVolumeFraction")
      ewrite(3,*) 'size of pvf:', node_count(phasevolumefraction)
      do i=1,node_count(phasevolumefraction)
        call set(phasevolumefraction, i, saturations(i))
      enddo
    
      ! Then to get the projection of PhaseVolumeFraction onto the CoordinateMesh
      ! so that it can be compared to the analytical solution
      galerkinprojection => extract_scalar_field(state(1), "GalerkinProjection")
      call calculate_diagnostic_variable(state, 1, galerkinprojection)


      ! Let's practice with getting the pressure out
      ! This is designed for the quadratic elements
      pressure => extract_scalar_field(state(1), "Pressure")    
      allocate(elenodes(3*positions%dim))
      do i=1,ele_count(pressure)
        elenodes = ele_nodes(pressure,i)
        do j=1,size(elenodes)
          call set(pressure, elenodes(j), proto_pressure(2*(i-1)+j))
        end do
      end do  
      deallocate(elenodes)
    
      ewrite(3,*) 'Leaving copy_into_state'
  
    end subroutine copy_into_state

end module copy_outof_into_state
