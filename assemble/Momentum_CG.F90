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
!    C.Pain@Imperial.ac.uk
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

  module momentum_cg

    use fields
    use state_module
    use spud
    use fldebug
    use sparse_tools
    use boundary_conditions
    use boundary_conditions_from_options
    use solvers, only: petsc_solve
    use sparse_tools_petsc
    use sparse_matrices_fields
    use field_options
    use halos
    use global_parameters, only: FIELD_NAME_LEN
    use elements
    use transform_elements, only: transform_to_physical
    use coriolis_module
    use vector_tools
    use fetools
    use upwind_stabilisation
    use les_viscosity_module
    use metric_tools
    use field_derivatives
    use state_fields_module
    use state_matrices_module
    use sparsity_patterns_meshes
    use fefields
    use rotated_boundary_conditions
    implicit none

    private
    public :: construct_momentum_cg, correct_masslumped_velocity, &
              correct_velocity_cg, assemble_masslumped_poisson_rhs, &
              add_kmk_matrix, add_kmk_rhs, assemble_kmk_matrix, &
              deallocate_cg_mass

    ! are we lumping the mass, absorption or source
    logical :: lump_mass, lump_absorption, lump_source
    ! is the pressure correction included in the absorption term?
    ! if so, lump_absorption gets set equal to lump_mass
    logical :: pressure_corrected_absorption
    ! do we have isotropic viscosity?
    logical :: isotropic_viscosity
    ! do we have diagonal viscosity?
    logical :: diagonal_viscosity
    ! are we using the stress form of the viscosity terms?
    logical :: stress_form
    ! do we want to integrate the continuity matrix by parts?
    logical :: integrate_continuity_by_parts
    ! assemble the viscosity as it used to be done in diff3d
    logical :: legacy_stress
    ! exclude the advection or mass terms from the equation
    logical :: exclude_advection, exclude_mass
    ! integrate the advection term by parts
    logical :: integrate_advection_by_parts
    ! do we need the inverse lumped mass to assemble a lumped cmc preconditioner
    logical :: cmc_lump_mass
    ! use the sub mesh to lump the mass
    logical :: vel_lump_on_submesh, cmc_lump_on_submesh, abs_lump_on_submesh
    ! integrate the surface tension by parts
    logical :: integrate_surfacetension_by_parts
    
    ! which terms do we have?
    logical :: have_source
    logical :: have_gravity
    logical :: have_absorption
    logical :: have_viscosity
    logical :: have_surfacetension
    logical :: have_coriolis
    logical :: have_geostrophic_pressure
    logical :: have_les
    logical :: les_fourth_order
    
    logical :: move_mesh
    
    ! assemble mass or inverse lumped mass?
    logical :: assemble_mass_matrix
    logical :: assemble_inverse_masslump

    ! implicitness parameter, timestep, conservation parameter
    real :: theta, dt, beta, gravity_magnitude

    ! Stabilisation schemes.
    integer :: stabilisation_scheme
    integer, parameter :: STABILISATION_NONE=0
    integer, parameter :: STABILISATION_STREAMLINE_UPWIND=1, &
      & STABILISATION_SUPG=2
    integer :: nu_bar_scheme
    real :: nu_bar_scale = 1.0
    
    ! LES coefficients
    real :: smagorinsky_coefficient

  contains

    subroutine construct_momentum_cg(u, p, density, x, &
                                     big_m, rhs, ct_m, ct_rhs, mass, inverse_masslump, &
                                     state, assemble_ct_matrix, cg_pressure)
      !!< Assembles the momentum matrix and rhs for the LinearMomentum,
      !!< Boussinesq and Drainage equation types such that
      !!< big_m*u = rhs + ct_m*p
      !!<
      !!< This subroutine is intended to replace assnav and all new code added to it
      !!< should be in new format and be compatible with both 2 and 3 dimensions.
      !!<
      !!< For clarity big_m is assumed to always be a dim x dim block_csr_matrix even 
      !!< when velocities aren't coupled

      ! velocity and coordinate
      type(vector_field), intent(inout) :: u, x
      ! pressure and density
      type(scalar_field), intent(inout) :: p, density
      ! the lhs matrix
      type(petsc_csr_matrix), intent(inout) :: big_m
      
      ! the mass matrix
      ! NOTE: see the logical assemble_mass below to see when this is actually assembled
      type(petsc_csr_matrix), intent(inout) :: mass
      ! the lumped mass matrix (may vary per component as absorption could be included)
      ! NOTE: see the logical assemble_inverse_masslump below to see when this is actually assembled
      type(vector_field), intent(inout) :: inverse_masslump
      ! NOTE: you have to call deallocate_cg_mass after you're done
      ! with mass and inverse_masslump
      
      ! the pressure gradient matrix (might be null if assemble_ct_matrix=.false.)
      type(block_csr_matrix), pointer :: ct_m
      ! the pressure gradient rhs
      type(scalar_field), intent(inout), optional :: ct_rhs
      ! the rhs
      type(vector_field), intent(inout) :: rhs
      ! bucket full of fields
      type(state_type), intent(inout) :: state
      ! do we want to get a cg pressure gradient matrix?
      logical, intent(in) :: assemble_ct_matrix, cg_pressure

      type(scalar_field), pointer :: buoyancy
      type(scalar_field), pointer :: gp
      type(vector_field), pointer :: gravity
      type(vector_field), pointer :: oldu, nu, ug, source, absorption
      type(tensor_field), pointer :: viscosity
      type(tensor_field), pointer :: surfacetension
      type(vector_field), pointer :: x_old, x_new

      ! dummy fields in case state doesn't contain the above fields
      type(scalar_field), pointer :: dummyscalar
      type(vector_field), pointer :: dummyvector
      type(tensor_field), pointer :: dummytensor

      ! single component of lumped mass
      type(scalar_field) :: masslump_component        
      ! sparsity for mass matrices
      type(csr_sparsity), pointer :: u_sparsity

      ! bc arrays
      type(vector_field) :: velocity_bc
      type(scalar_field) :: pressure_bc
      integer, dimension(:,:), allocatable :: velocity_bc_type
      integer, dimension(:), allocatable :: pressure_bc_type

      ! fields for the assembly of absorption when
      ! lumping on the submesh
      type(vector_field) :: abslump
      type(scalar_field) :: absdensity, abslump_component, abs_component

      ! for 4th order les:
      type(tensor_field):: grad_u

      integer :: stat, dim, ele, sele, dim2
      
      ewrite(1,*) 'entering construct_momentum_cg'
    
      assert(continuity(u)>=0)

      nu=>extract_vector_field(state, "NonlinearVelocity")
      oldu=>extract_vector_field(state, "OldVelocity")

      allocate(dummyscalar)
      call allocate(dummyscalar, u%mesh, "DummyScalar", field_type=FIELD_TYPE_CONSTANT)
      call zero(dummyscalar)
      dummyscalar%option_path=""

      allocate(dummyvector)
      call allocate(dummyvector, u%dim, u%mesh, "DummyVector", field_type=FIELD_TYPE_CONSTANT)
      call zero(dummyvector)
      dummyvector%option_path=""

      allocate(dummytensor)
      call allocate(dummytensor, u%mesh, "DummyTensor", field_type=FIELD_TYPE_CONSTANT)
      call zero(dummytensor)
      dummytensor%option_path=""

      source=>extract_vector_field(state, "VelocitySource", stat)
      have_source = stat == 0
      if(.not. have_source) source=>dummyvector
      do dim = 1, source%dim
        ewrite_minmax(source%val(dim)%ptr(:))
      end do

      absorption=>extract_vector_field(state, "VelocityAbsorption", stat)
      have_absorption = stat == 0
      if(.not. have_absorption) absorption=>dummyvector
      do dim = 1, absorption%dim
        ewrite_minmax(absorption%val(dim)%ptr(:))
      end do

      call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude, &
          stat=stat)
      have_gravity = stat == 0
      if(have_gravity) then
        buoyancy=>extract_scalar_field(state, "VelocityBuoyancyDensity")
        gravity=>extract_vector_field(state, "GravityDirection", stat)
      else
        buoyancy=>dummyscalar
        gravity=>dummyvector
        gravity_magnitude = 0.0
      end if
      ewrite_minmax(buoyancy%val)

      viscosity=>extract_tensor_field(state, "Viscosity", stat)
      have_viscosity = stat == 0
      if(.not. have_viscosity) then
         viscosity=>dummytensor
      else
         do dim = 1, viscosity%dim
            do dim2 = 1, viscosity%dim
              if(dim2<dim) cycle
              ewrite_minmax(viscosity%val(dim,dim2,:))
            end do
         end do
      end if
      
      surfacetension=>extract_tensor_field(state, "VelocitySurfaceTension", stat)
      have_surfacetension = stat == 0
      if(.not. have_surfacetension) then
         surfacetension=>dummytensor
      else
         do dim = 1, surfacetension%dim
            do dim2 = 1, surfacetension%dim
              if(dim2<dim) cycle
              ewrite_minmax(surfacetension%val(dim,dim2,:))
            end do
         end do
      end if
      
      have_coriolis = have_option("/physical_parameters/coriolis")
      have_les = have_option(trim(u%option_path)//"/prognostic/spatial_discretisation/&
         &/continuous_galerkin/les_model")
      if (have_les) then
         call get_option(trim(u%option_path)//"/prognostic/spatial_discretisation/&
            &/continuous_galerkin/les_model/smagorinsky_coefficient", &
            smagorinsky_coefficient)
         les_fourth_order=have_option(trim(u%option_path)//"/prognostic/spatial_discretisation/&
            &/continuous_galerkin/les_model/order/fourth_order")
         if (les_fourth_order) then
           call allocate( grad_u, u%mesh, "VelocityGradient")
           call differentiate_field_lumped( nu, x, grad_u)
         end if
      end if
      
      have_geostrophic_pressure = has_scalar_field(state, "GeostrophicPressure")
      if(have_geostrophic_pressure) then
        gp => extract_scalar_field(state, "GeostrophicPressure")
        
        ewrite_minmax(gp%val)
      else
        gp => dummyscalar
      end if

      call get_option("/timestepping/timestep", dt)
      call get_option(trim(u%option_path)//"/prognostic/temporal_discretisation/theta", &
                      theta)
      call get_option(trim(u%option_path)//"/prognostic/spatial_discretisation/&
           &conservative_advection", beta)

      lump_mass=have_option(trim(u%option_path)//&
          &"/prognostic/spatial_discretisation&
          &/continuous_galerkin/mass_terms/lump_mass_matrix")
      lump_absorption=have_option(trim(u%option_path)//&
          &"/prognostic/vector_field::Absorption&
          &/lump_absorption")
      abs_lump_on_submesh = have_option(trim(u%option_path)//&
          &"/prognostic/vector_field::Absorption&
          &/lump_absorption/use_submesh")
      pressure_corrected_absorption=have_option(trim(u%option_path)//&
          &"/prognostic/vector_field::Absorption&
          &/include_pressure_correction")
      if (pressure_corrected_absorption) then
         ! as we add the absorption into the mass matrix
         ! lump_absorption needs to match lump_mass
         lump_absorption = lump_mass
      end if
      lump_source=have_option(trim(u%option_path)//&
          &"/prognostic/vector_field::Source&
          &/lump_source")
      if(have_viscosity) then
         isotropic_viscosity = have_viscosity .and. &
           & isotropic_field(viscosity)
         diagonal_viscosity = have_viscosity .and. &
           & diagonal_field(viscosity)
         stress_form=have_option(trim(u%option_path)//&
             &"/prognostic/spatial_discretisation/continuous_galerkin&
             &/stress_terms/stress_form")
      else
         isotropic_viscosity = .false.
         diagonal_viscosity = .false.
         stress_form = .false.
      end if
      integrate_continuity_by_parts=have_option(trim(p%option_path)//&
          &"/prognostic/spatial_discretisation/continuous_galerkin&
          &/integrate_continuity_by_parts")
      legacy_stress = have_option(trim(u%option_path)//&
          &"/prognostic/spatial_discretisation&
          &/continuous_galerkin/stress_terms/stress_form/legacy_stress_form")
      integrate_advection_by_parts = have_option(trim(u%option_path)//&
          &"/prognostic/spatial_discretisation&
          &/continuous_galerkin/advection_terms/integrate_advection_by_parts")
      exclude_advection = have_option(trim(u%option_path)//&
          &"/prognostic/spatial_discretisation&
          &/continuous_galerkin/advection_terms/exclude_advection_terms")
      exclude_mass = have_option(trim(u%option_path)//&
          &"/prognostic/spatial_discretisation&
          &/continuous_galerkin/mass_terms/exclude_mass_terms")
      vel_lump_on_submesh = have_option(trim(u%option_path)//&
          &"/prognostic/spatial_discretisation&
          &/continuous_galerkin/mass_terms&
          &/lump_mass_matrix/use_submesh")
      if (pressure_corrected_absorption) then
         ! as we add the absorption into the mass matrix
         ! the meshes need to be the same
         abs_lump_on_submesh = vel_lump_on_submesh
      end if
      cmc_lump_mass = have_option(trim(p%option_path)//&
          &"/prognostic/scheme&
          &/use_projection_method/full_schur_complement&
          &/preconditioner_matrix::LumpedSchurComplement")
      cmc_lump_on_submesh = have_option(trim(p%option_path)//&
          &"/prognostic/scheme&
          &/use_projection_method/full_schur_complement&
          &/preconditioner_matrix[0]/lump_on_submesh")
      assemble_inverse_masslump = lump_mass .or. cmc_lump_mass
      assemble_mass_matrix = have_option(trim(p%option_path)//&
          "/prognostic/scheme/use_projection_method&
          &/full_schur_complement/inner_matrix::FullMassMatrix")
      if(have_option(trim(u%option_path)//"/prognostic/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind")) then
        stabilisation_scheme = STABILISATION_STREAMLINE_UPWIND
        call get_upwind_options(trim(u%option_path) // "/prognostic/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind", &
          & nu_bar_scheme, nu_bar_scale)
      else if(have_option(trim(u%option_path)//"/prognostic/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind_petrov_galerkin")) then
        stabilisation_scheme = STABILISATION_SUPG
        call get_upwind_options(trim(u%option_path) // "/prognostic/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind_petrov_galerkin", &
          & nu_bar_scheme, nu_bar_scale)
      else
        stabilisation_scheme = STABILISATION_NONE
      end if
      integrate_surfacetension_by_parts = have_option(trim(u%option_path)//&
          &"/prognostic/tensor_field::SurfaceTension&
          &/diagnostic/integrate_by_parts")
          
      if (assemble_inverse_masslump) then
        ! construct the inverse of the lumped mass matrix
        call allocate( inverse_masslump, u%dim, u%mesh, "InverseLumpedMass")
        call zero(inverse_masslump)
      end if
      if (assemble_mass_matrix) then
        ! construct mass matrix instead
        u_sparsity => get_csr_sparsity_firstorder(state, u%mesh, u%mesh)
        
        call allocate( mass, u_sparsity, (/ u%dim, u%dim /), &
            diagonal=.true., name="MassMatrix")
            
        call zero( mass )
      end if
      
      move_mesh = (have_option("/mesh_adaptivity/mesh_movement").and.(.not.exclude_mass))
      if(move_mesh) then
        ewrite(2,*) 'Moving mesh'
        x_old => extract_vector_field(state, "OldCoordinate")
        x_new => extract_vector_field(state, "IteratedCoordinate")
        ug=>extract_vector_field(state, "GridVelocity")
      else
        ewrite(2,*) 'Not moving mesh'
      end if
      
      ! ----- Volume integrals over elements -------------
      
      element_loop: do ele=1, element_count(u)
         call construct_momentum_element_cg(ele, big_m, rhs, ct_m, mass, inverse_masslump, &
              x, x_old, x_new, u, oldu, nu, ug, &
              density, p, &
              source, absorption, buoyancy, gravity, &
              viscosity, grad_u, &
              gp, surfacetension, &
              assemble_ct_matrix, cg_pressure)
         
      end do element_loop
        
      ! ----- Surface integrals over boundaries -----------
      
      if((integrate_advection_by_parts.and.(.not.exclude_advection)).or.&
           (integrate_continuity_by_parts)) then
         allocate(velocity_bc_type(u%dim, surface_element_count(u)))
         call get_entire_boundary_condition(u, &
           & (/ &
             "weakdirichlet      ", &
             "no_normal_flow     ", &
             "periodic           ", &
             "free_surface       " &
           & /), velocity_bc, velocity_bc_type)
           
         allocate(pressure_bc_type(surface_element_count(p)))
         call get_entire_boundary_condition(p, &
           & (/ &
              "weakdirichlet", &
              "dirichlet    " /), &
              pressure_bc, pressure_bc_type)

         surface_element_loop: do sele=1, surface_element_count(u)
            
            ! if no_normal flow and no other condition in the tangential directions, or if periodic
            ! but not if there's a pressure bc
            if(((velocity_bc_type(1,sele)==2 .and. sum(velocity_bc_type(:,sele))==2) &
                 .or. any(velocity_bc_type(:,sele)==3)) &
                 .and. pressure_bc_type(sele)==0) cycle
            
            ele = face_ele(x, sele)
            
            call construct_momentum_surface_element_cg(sele, ele, big_m, rhs, ct_m, ct_rhs, &
                 x, u, nu, ug, density, p, &
                 velocity_bc, velocity_bc_type, &
                 pressure_bc, pressure_bc_type, &
                 assemble_ct_matrix, cg_pressure, viscosity, oldu)
            
         end do surface_element_loop
         
         call deallocate(velocity_bc)
         deallocate(velocity_bc_type)
         call deallocate(pressure_bc)
         deallocate(pressure_bc_type)
      end if
      
      if(abs_lump_on_submesh) then
        
        call allocate(abslump, inverse_masslump%dim, inverse_masslump%mesh, "LumpedAbsorption")
        call allocate(absdensity, absorption%mesh, "AbsorptionComponentTimesDensity")
        
        do dim = 1, inverse_masslump%dim
          call remap_field(density, absdensity)
          abs_component = extract_scalar_field(absorption, dim)
          call scale(absdensity, abs_component)
      
          abslump_component = extract_scalar_field(abslump, dim)              
          call compute_lumped_mass_on_submesh(state, abslump_component, density=absdensity)
        end do
        
        call deallocate(absdensity)
        
        if(assemble_inverse_masslump.and.pressure_corrected_absorption) then
          call addto(inverse_masslump, abslump, theta)
        end if

        call addto_diag(big_m, abslump, dt*theta)
        
        call scale(abslump, oldu)
        call addto(rhs, abslump, -1.0)
          
        call deallocate(abslump)
      end if

      if (assemble_inverse_masslump) then
        
        if(vel_lump_on_submesh .or. cmc_lump_on_submesh) then
          if(move_mesh) then
            FLAbort("Can't move the mesh and lump on the submesh yet.")
          end if
          ! we still have to make the lumped mass if this is true
          masslump_component=extract_scalar_field(inverse_masslump, 1)

          call compute_lumped_mass_on_submesh(state, masslump_component, density=density)

          ! copy over to other components
          do dim = 2, inverse_masslump%dim
            call set(inverse_masslump, dim, masslump_component)
          end do

          if(vel_lump_on_submesh) then
            call addto_diag(big_m, masslump_component)
          end if
        end if
        
        ! thus far we have just assembled the lumped mass in inverse_masslump
        ! now invert it:
        call invert(inverse_masslump)
        ! apply boundary conditions (zeroing out strong dirichl. rows)
        call apply_dirichlet_conditions_inverse_mass(inverse_masslump, u)
        
        do dim = 1, rhs%dim
          ewrite_minmax(inverse_masslump%val(dim)%ptr)
        end do
      end if
      
      if (assemble_mass_matrix) then
        call apply_dirichlet_conditions(matrix=mass, field=u)
      end if
            
      do dim = 1, rhs%dim
        ewrite_minmax(rhs%val(dim)%ptr)
      end do

      if (les_fourth_order) then
        call deallocate(grad_u)
      end if

      call deallocate(dummytensor)
      deallocate(dummytensor)
      call deallocate(dummyvector)
      deallocate(dummyvector)
      call deallocate(dummyscalar)
      deallocate(dummyscalar)

    end subroutine construct_momentum_cg

    subroutine construct_momentum_surface_element_cg(sele, ele, big_m, rhs, ct_m, ct_rhs, &
                                                     x, u, nu, ug, density, p, &
                                                     velocity_bc, velocity_bc_type, &
                                                     pressure_bc, pressure_bc_type, &
                                                     assemble_ct_matrix, cg_pressure, viscosity, &
                                                     oldu)

      integer, intent(in) :: sele, ele

      type(petsc_csr_matrix), intent(inout) :: big_m
      type(vector_field), intent(inout) :: rhs

      type(block_csr_matrix), pointer :: ct_m
      type(scalar_field), intent(inout) :: ct_rhs

      type(vector_field), intent(in) :: x, oldu
      type(vector_field), intent(in) :: u, nu
      type(vector_field), pointer :: ug
      type(scalar_field), intent(in) :: density, p
      type(tensor_field), intent(in) :: viscosity

      type(vector_field), intent(in) :: velocity_bc
      integer, dimension(:,:), intent(in) :: velocity_bc_type

      type(scalar_field), intent(in) :: pressure_bc
      integer, dimension(:), intent(in) :: pressure_bc_type
      
      logical, intent(in) :: assemble_ct_matrix, cg_pressure

      ! local
      integer :: dim

      integer, dimension(face_loc(u, sele)) :: u_nodes_bdy
      integer, dimension(face_loc(p, sele)) :: p_nodes_bdy
      type(element_type), pointer :: u_shape, p_shape

      real, dimension(face_ngi(u, sele)) :: detwei_bdy
      real, dimension(u%dim, face_ngi(u, sele)) :: normal_bdy
      real, dimension(u%dim, face_loc(p, sele), face_loc(u, sele)) :: ct_mat_bdy
      real, dimension(face_loc(u, sele), face_loc(u, sele)) :: adv_mat_bdy

      real, dimension(u%dim, face_ngi(u, sele)) :: relu_gi

      u_shape=> face_shape(u, sele)
      p_shape=> face_shape(p, sele)

      u_nodes_bdy = face_global_nodes(u, sele)
      p_nodes_bdy = face_global_nodes(p, sele)

      call transform_facet_to_physical(X, sele, &
           detwei_f=detwei_bdy, normal=normal_bdy)
                                     
      ! Note that with SUPG the surface element test function is not modified
            
      ! first the advection (dirichlet) bcs:
      
      ! if no no_normal_flow or free_surface
      if (velocity_bc_type(1,sele)/=2) then
         if(integrate_advection_by_parts.and.(.not.exclude_advection)) then
            
            relu_gi = face_val_at_quad(nu, sele)
            if(move_mesh) then
              relu_gi = relu_gi - face_val_at_quad(ug, sele)
            end if
            
            adv_mat_bdy = shape_shape(u_shape, u_shape, &
                 detwei_bdy*sum(relu_gi*normal_bdy,1)*&
                 face_val_at_quad(density, sele))
            do dim = 1, u%dim
               
               if(velocity_bc_type(dim, sele)==1) then

                  call addto(rhs, dim, u_nodes_bdy, -matmul(adv_mat_bdy, &
                       ele_val(velocity_bc, dim, sele)))
               else

                  call addto(big_m, dim, dim, u_nodes_bdy, u_nodes_bdy, &
                       dt*theta*adv_mat_bdy)

                  call addto(rhs, dim, u_nodes_bdy, -matmul(adv_mat_bdy, face_val(oldu, dim, sele)))

               end if
            end do
         end if
      end if
      
      ! now do surface integrals for divergence/pressure gradient matrix
      if(integrate_continuity_by_parts.and.cg_pressure) then
         
        if (velocity_bc_type(1,sele)/=2 .and. velocity_bc_type(1,sele)/=4) then

          ct_mat_bdy = shape_shape_vector(p_shape, u_shape, detwei_bdy, normal_bdy)
          do dim = 1, u%dim
             if(velocity_bc_type(dim, sele)==1 )then
                call addto(ct_rhs, p_nodes_bdy, &
                     -matmul(ct_mat_bdy(dim,:,:), ele_val(velocity_bc, dim, sele)))
             else if (assemble_ct_matrix) then
                call addto(ct_m, 1, dim, p_nodes_bdy, u_nodes_bdy, ct_mat_bdy(dim,:,:))
             end if
             if(pressure_bc_type(sele)>0) then
                ! for both weak and strong pressure dirichlet bcs:
                !      /
                ! add -|  N_i M_j \vec n p_j, where p_j are the prescribed bc values
                !      /
                call addto(rhs, dim, u_nodes_bdy, -matmul( ele_val(pressure_bc, sele), &
                                                            ct_mat_bdy(dim,:,:) ))
             end if
          end do
        end if
        
      end if

    end subroutine construct_momentum_surface_element_cg


    subroutine construct_momentum_element_cg(ele, big_m, rhs, ct_m, &
                                            mass, masslump, &
                                            x, x_old, x_new, u, oldu, nu, ug, &
                                            density, p, &
                                            source, absorption, buoyancy, gravity, &
                                            viscosity, grad_u, &
                                            gp, surfacetension, &
                                            assemble_ct_matrix, cg_pressure)
    !!< Assembles the local element matrix contributions and places them in big_m
    !!< and rhs for the continuous galerkin momentum equations

      ! current element
      integer, intent(in) :: ele
      type(petsc_csr_matrix), intent(inout) :: big_m
      type(vector_field), intent(inout) :: rhs
      type(block_csr_matrix), pointer :: ct_m
      type(petsc_csr_matrix), intent(inout) :: mass
      ! above we supply inverse_masslump, but we start assembling the non-inverted
      ! lumped mass matrix in it:
      type(vector_field), intent(inout) :: masslump

      type(vector_field), intent(in) :: x, u, oldu, nu 
      type(vector_field), pointer :: x_old, x_new, ug
      type(scalar_field), intent(in) :: density, p, buoyancy
      type(vector_field), intent(in) :: source, absorption, gravity
      type(tensor_field), intent(in) :: viscosity
      type(scalar_field), intent(in) :: gp
      type(tensor_field), intent(in) :: surfacetension
      type(tensor_field), intent(in) :: grad_u

      logical, intent(in) :: assemble_ct_matrix, cg_pressure

      integer, dimension(:), pointer :: u_ele, p_ele
      real, dimension(u%dim, ele_loc(u, ele)) :: oldu_val
      type(element_type), pointer :: u_shape, p_shape
      real, dimension(ele_ngi(u, ele)) :: detwei, detwei_old, detwei_new
      real, dimension(u%dim, u%dim, ele_ngi(u,ele)) :: J_mat
      real, dimension(ele_loc(u, ele), ele_ngi(u, ele), u%dim) :: du_t
      real, dimension(ele_loc(u, ele), ele_ngi(u, ele), u%dim) :: dug_t
      real, dimension(ele_loc(p, ele), ele_ngi(p, ele), u%dim) :: dp_t

      real, dimension(u%dim, ele_ngi(u, ele)) :: relu_gi
      real, dimension(u%dim, ele_loc(p, ele), ele_loc(u, ele)) :: grad_p_u_mat
      
      ! What we will be adding to the matrix and RHS - assemble these as we
      ! go, so that we only do the calculations we really need
      real, dimension(u%dim, ele_loc(u, ele)) :: big_m_diag_addto, rhs_addto
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)) :: big_m_tensor_addto
      logical, dimension(u%dim, u%dim) :: block_mask ! control whether the off diagonal entries are used
      integer :: dim
      type(element_type) :: test_function

      if(move_mesh) then
        ! we've assumed the following in the declarations
        ! above so we better make sure they're true!
        assert(ele_loc(ug, ele)==ele_loc(u,ele))
        assert(ele_ngi(ug, ele)==ele_ngi(u,ele))
        assert(ug%dim==u%dim)
      end if
      
      big_m_diag_addto = 0.0
      big_m_tensor_addto = 0.0
      rhs_addto = 0.0
      ! we always want things added to the diagonal blocks
      ! but we must check if we have_coriolis to add things to the others
      if(have_coriolis.or.(have_viscosity.and.stress_form)) then
        block_mask = .true.
      else
        block_mask = .false.
        do dim = 1, u%dim
          block_mask(dim, dim) = .true.
        end do
      end if

      u_ele=>ele_nodes(u, ele)
      u_shape=>ele_shape(u, ele)

      p_ele=>ele_nodes(p, ele)
      p_shape=>ele_shape(p, ele)

      oldu_val = ele_val(oldu, ele)
      ! Step 1: Transform

      ! transform the velocity derivatives into physical space
      ! (and get detwei)
      if(stabilisation_scheme==STABILISATION_NONE) then
        call transform_to_physical(X, ele, &
                                  u_shape, dshape=du_t, detwei=detwei)
      !  J_mat = 0.0
      else
        call transform_to_physical(x, ele, &
                                  u_shape, dshape=du_t, detwei=detwei, J=J_mat)
      end if

      if(assemble_ct_matrix.and.cg_pressure.and.integrate_continuity_by_parts) then
        ! transform the pressure derivatives into physical space
        call transform_to_physical(x, ele, &
                                  p_shape, dshape=dp_t)
      !else
      !  dp_t = 0.0
      end if
      
      if(move_mesh) then
        call transform_to_physical(x_old, ele, detwei=detwei_old)
        call transform_to_physical(x_new, ele, detwei=detwei_new)
        if(.not.exclude_advection.and..not.integrate_advection_by_parts) then
          call transform_to_physical(x, ele, &
                                    ele_shape(ug, ele), dshape=dug_t)
        end if
      end if
      
      ! Step 2: Set up test function
    
      select case(stabilisation_scheme)
        case(STABILISATION_SUPG)
          relu_gi = ele_val_at_quad(nu, ele)
          if(move_mesh) then
            relu_gi = relu_gi - ele_val_at_quad(ug, ele)
          end if
          if(have_viscosity) then
            test_function = make_supg_shape(u_shape, du_t, relu_gi, j_mat, diff_q = ele_val_at_quad(viscosity, ele), &
              & nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
          else
            test_function = make_supg_shape(u_shape, du_t, relu_gi, j_mat, &
              & nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
          end if
        case default
          test_function = u_shape
          call incref(test_function)
      end select
      ! Important note: the test function derivatives have not been modified -
      ! i.e. du_t is currently used everywhere. This is fine for P1, but is not
      ! consistent for P>1.

      if(assemble_ct_matrix.and.cg_pressure) then
        if(integrate_continuity_by_parts) then
          grad_p_u_mat = -dshape_shape(dp_t, u_shape, detwei)
        else
          grad_p_u_mat = shape_dshape(p_shape, du_t, detwei)
        end if
      !else
      !  grad_p_u_mat = 0.0
      end if

      ! Step 3: Assemble contributions

      ! Mass terms
      if(assemble_inverse_masslump .or. assemble_mass_matrix .or. &
        (.not. exclude_mass)) then
        call add_mass_element_cg(ele, test_function, u, oldu_val, density, detwei, detwei_old, detwei_new, big_m_diag_addto, big_m_tensor_addto, rhs_addto, mass, masslump)
      end if

      ! Advection terms
      if(.not. exclude_advection) then
        call add_advection_element_cg(ele, test_function, u, oldu_val, nu, ug, density, viscosity, du_t, dug_t, detwei, J_mat, big_m_tensor_addto, rhs_addto)
      end if

      ! Source terms
      if(have_source) then
        call add_sources_element_cg(ele, test_function, u, density, source, detwei, rhs_addto)
      end if
      
      ! Buoyancy terms
      if(have_gravity) then
        call add_buoyancy_element_cg(ele, test_function, u, buoyancy, gravity, detwei, rhs_addto)
      end if
      
      ! Surface tension
      if(have_surfacetension) then
        call add_surfacetension_element_cg(ele, test_function, u, surfacetension, du_t, detwei, rhs_addto)
      end if

      ! Absorption terms (sponges)
      if(have_absorption) then
       call add_absorption_element_cg(ele, test_function, u, oldu_val, density, &
         absorption, detwei, big_m_diag_addto, big_m_tensor_addto, rhs_addto, &
         masslump, mass)
      end if
      
      ! Viscous terms
      if(have_viscosity .or. have_les) then
        call add_viscosity_element_cg(ele, u, oldu_val, nu, x, viscosity, grad_u, &
           du_t, detwei, big_m_tensor_addto, rhs_addto)
      end if
      
      ! Coriolis terms
      if(have_coriolis) then
        call add_coriolis_element_cg(ele, test_function, x, u, oldu_val, density, detwei, big_m_tensor_addto, rhs_addto)
      end if
      
      ! Geostrophic pressure
      if(have_geostrophic_pressure) then
        call add_geostrophic_pressure_element_cg(ele, test_function, x, u, gp, detwei, rhs_addto)
      end if

      ! Step 4: Insertion

      ! add lumped terms to the diagonal of the matrix
      call add_diagonal_to_tensor(big_m_diag_addto, big_m_tensor_addto)
      ! add to the matrix
      call addto(big_m, u_ele, u_ele, big_m_tensor_addto, block_mask=block_mask)
      ! add to the rhs
      call addto(rhs, u_ele, rhs_addto)
      
      if(assemble_ct_matrix.and.cg_pressure) then
        call addto(ct_m, p_ele, u_ele, spread(grad_p_u_mat, 1, 1))
      end if
      
      call deallocate(test_function)
      
    contains
    
      subroutine add_diagonal_to_tensor(big_m_diag_addto, big_m_tensor_addto)
        real, dimension(u%dim, ele_loc(u, ele)), intent(in) :: big_m_diag_addto
        real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)), intent(inout) :: big_m_tensor_addto
        
        integer :: dim, loc
        
        forall(dim = 1:size(big_m_diag_addto, 1), loc = 1:size(big_m_diag_addto, 2))
          big_m_tensor_addto(dim, dim, loc, loc) = big_m_tensor_addto(dim, dim, loc, loc) + big_m_diag_addto(dim, loc)
        end forall
        
      end subroutine add_diagonal_to_tensor
             
    end subroutine construct_momentum_element_cg
    
    subroutine add_mass_element_cg(ele, test_function, u, oldu_val, density, detwei, detwei_old, detwei_new, big_m_diag_addto, big_m_tensor_addto, rhs_addto, mass, masslump)
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: u
      real, dimension(:,:), intent(in) :: oldu_val
      type(scalar_field), intent(in) :: density
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei, detwei_old, detwei_new
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: big_m_diag_addto
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)), intent(inout) :: big_m_tensor_addto
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
      type(petsc_csr_matrix), intent(inout), optional :: mass
      type(vector_field), intent(inout) :: masslump
      
      integer :: dim
      integer, dimension(:), pointer :: u_ele
      logical:: compute_lumped_mass_here
      real, dimension(ele_loc(u, ele)) :: mass_lump
      real, dimension(ele_ngi(u, ele)) :: density_gi
      real, dimension(ele_loc(u, ele), ele_loc(u, ele)) :: mass_mat
      type(element_type), pointer :: u_shape
      
      u_shape => ele_shape(u, ele)
      u_ele=>ele_nodes(u, ele)
      
      density_gi=ele_val_at_quad(density, ele)
            
      ! element mass matrix
      !  /
      !  | N_A N_B rho dV
      !  /
      if(move_mesh) then
        mass_mat = shape_shape(test_function, u_shape, detwei_new*density_gi)
      else
        mass_mat = shape_shape(test_function, u_shape, detwei*density_gi)
      end if
      mass_lump = sum(mass_mat, 2)
      
      ! if we're lumping on the submesh, this is done later:
      compute_lumped_mass_here=.not. (vel_lump_on_submesh .or. cmc_lump_on_submesh)
      
      if(.not.exclude_mass) then
        if(lump_mass) then
          if (compute_lumped_mass_here) then
            do dim = 1, u%dim
              big_m_diag_addto(dim, :) = big_m_diag_addto(dim, :) + mass_lump
            end do
          end if
        else
          do dim = 1, u%dim
            big_m_tensor_addto(dim, dim, :, :) = big_m_tensor_addto(dim, dim, :, :) + mass_mat
          end do
        end if
      end if
            
      if(assemble_inverse_masslump .and. compute_lumped_mass_here) then
        ! store the lumped mass as field, the same for each component
        do dim = 1, u%dim
           call addto(masslump, dim, u_ele, mass_lump)
        end do
      end if
      
      if(assemble_mass_matrix) then
         do dim=1, u%dim
            call addto(mass, dim, dim, u_ele, u_ele, mass_mat)
         end do
      end if
      
      if(move_mesh) then
        ! In the unaccelerated form we solve:
        !  /
        !  |  N^{n+1} u^{n+1}/dt - N^{n} u^n/dt + ... = f
        !  /
        ! so in accelerated form:
        !  /
        !  |  N^{n+1} du + (N^{n+1}- N^{n}) u^n/dt + ... = f
        !  /
        ! where du=(u^{n+1}-u^{n})/dt is the acceleration.
        ! Put the (N^{n+1}-N^{n}) u^n term on the rhs
        mass_mat = shape_shape(test_function, u_shape, (detwei_new-detwei_old)*density_gi)
        if(lump_mass) then
          if(compute_lumped_mass_here) then
            mass_lump = sum(mass_mat, 2)
            do dim = 1, u%dim
              rhs_addto(dim,:) = rhs_addto(dim,:) - mass_lump*oldu_val(dim,:)/dt
            end do
          end if
        else
          do dim = 1, u%dim
            rhs_addto(dim,:) = rhs_addto(dim,:) - matmul(mass_mat, oldu_val(dim,:))/dt
          end do
        end if
      end if
      
    end subroutine add_mass_element_cg
    
    subroutine add_advection_element_cg(ele, test_function, u, oldu_val, nu, ug,  density, viscosity, du_t, dug_t, detwei, J_mat, big_m_tensor_addto, rhs_addto)
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: u
      real, dimension(:,:), intent(in) :: oldu_val
      type(vector_field), intent(in) :: nu
      type(vector_field), pointer :: ug
      type(scalar_field), intent(in) :: density
      type(tensor_field), intent(in) :: viscosity
      real, dimension(ele_loc(u, ele), ele_ngi(u, ele), u%dim), intent(in) :: du_t
      real, dimension(ele_loc(u, ele), ele_ngi(u, ele), u%dim), intent(in) :: dug_t
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei
      real, dimension(u%dim, u%dim, ele_ngi(u,ele)) :: J_mat
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)), intent(inout) :: big_m_tensor_addto
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
    
      integer :: dim
      real, dimension(ele_ngi(u, ele)) :: density_gi, div_relu_gi
      real, dimension(ele_loc(u, ele), ele_loc(u, ele)) :: advection_mat
      real, dimension(u%dim, ele_ngi(u, ele)) :: relu_gi
      type(element_type), pointer :: u_shape
      
      u_shape=>ele_shape(u, ele)
      
            
      density_gi=ele_val_at_quad(density, ele)
      relu_gi = ele_val_at_quad(nu, ele)
      if(move_mesh) then
        relu_gi = relu_gi - ele_val_at_quad(ug, ele)
      end if
      div_relu_gi = ele_div_at_quad(nu, ele, du_t)
            
      if(integrate_advection_by_parts) then
        ! element advection matrix
        !    /                                            /
        !  - | (grad N_A dot nu) N_B rho dV - (1. - beta) | N_A ( div nu ) N_B rho dV
        !    /                                            /
        advection_mat = -dshape_dot_vector_shape(du_t, relu_gi, u_shape, detwei*density_gi)  &
                      -(1.-beta)*shape_shape(test_function, u_shape, div_relu_gi*detwei*density_gi)
      else
        ! element advection matrix
        !  /                                     /
        !  | N_A (nu dot grad N_B) rho dV + beta | N_A ( div nu ) N_B rho dV
        !  /                                     /
        advection_mat = shape_vector_dot_dshape(test_function, relu_gi, du_t, detwei*density_gi)  &
                      +beta*shape_shape(test_function, u_shape, div_relu_gi*detwei*density_gi)
        if(move_mesh) then
          advection_mat = advection_mat - shape_shape(test_function, u_shape, ele_div_at_quad(ug, ele, dug_t)*detwei*density_gi)
        end if
      end if
      
      select case(stabilisation_scheme)
      case(STABILISATION_STREAMLINE_UPWIND)
        if(have_viscosity) then
          advection_mat = advection_mat + &
            & element_upwind_stabilisation(u_shape, du_t, relu_gi, J_mat, detwei, &
            & diff_q = ele_val_at_quad(viscosity, ele), nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
        else
           advection_mat = advection_mat + &
            & element_upwind_stabilisation(u_shape, du_t, relu_gi, J_mat, detwei, &
            & nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
        end if
      end select
      
      do dim = 1, u%dim
        big_m_tensor_addto(dim, dim, :, :) = big_m_tensor_addto(dim, dim, :, :) + dt*theta*advection_mat
        rhs_addto(dim, :) = rhs_addto(dim, :) - matmul(advection_mat, oldu_val(dim,:))
      end do
      
    end subroutine add_advection_element_cg
    
    subroutine add_sources_element_cg(ele, test_function, u, density, source, detwei, rhs_addto)
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: u
      type(scalar_field), intent(in) :: density
      type(vector_field), intent(in) :: source
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
      
      integer :: dim
      real, dimension(ele_ngi(u, ele)) :: density_gi
      real, dimension(ele_loc(u, ele)) :: source_lump
      real, dimension(ele_loc(u, ele), ele_loc(source, ele)) :: source_mat
      
      density_gi=ele_val_at_quad(density, ele)

      ! element source matrix
      !  /
      !  | N_A N_B rho dV
      !  /
      source_mat = shape_shape(test_function, ele_shape(source, ele), detwei*density_gi)
      if(lump_source) then
        assert(ele_loc(source, ele)==ele_loc(u, ele))
        source_lump = sum(source_mat, 2)
        do dim = 1, u%dim
          ! lumped source
          rhs_addto(dim, :) = rhs_addto(dim, :) + source_lump*ele_val(source, dim, ele)
        end do
      else
        do dim = 1, u%dim
          rhs_addto(dim, :) = rhs_addto(dim, :) + matmul(source_mat, ele_val(source, dim, ele))
        end do
      end if
      
    end subroutine add_sources_element_cg
    
    subroutine add_buoyancy_element_cg(ele, test_function, u, buoyancy, gravity, detwei, rhs_addto)
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: u
      type(scalar_field), intent(in) :: buoyancy
      type(vector_field), intent(in) :: gravity
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
      
      rhs_addto = rhs_addto + &
                  shape_vector_rhs(test_function, &
                                   ele_val_at_quad(gravity, ele), &
                                   detwei*gravity_magnitude*ele_val_at_quad(buoyancy, ele))
      
    end subroutine add_buoyancy_element_cg
    
    subroutine add_surfacetension_element_cg(ele, test_function, u, surfacetension, du_t, detwei, rhs_addto)
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: u
      type(tensor_field), intent(in) :: surfacetension
      real, dimension(ele_loc(u, ele), ele_ngi(u, ele), u%dim), intent(in) :: du_t
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
      
      real, dimension(u%dim, ele_ngi(u, ele)) :: dtensiondj
      real, dimension(u%dim, u%dim, ele_ngi(u, ele)) :: tension
            
      if(integrate_surfacetension_by_parts) then
        tension = ele_val_at_quad(surfacetension, ele)
        
        rhs_addto = rhs_addto - dshape_dot_tensor_rhs(du_t, tension, detwei)
      else
        dtensiondj = ele_div_at_quad_tensor(surfacetension, ele, du_t)
        
        rhs_addto = rhs_addto + shape_vector_rhs(test_function,dtensiondj,detwei)
      end if
      
    end subroutine add_surfacetension_element_cg
    
    subroutine add_absorption_element_cg(ele, test_function, u, oldu_val, density, absorption, detwei, big_m_diag_addto, big_m_tensor_addto, rhs_addto, masslump, mass)
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: u
      real, dimension(:,:), intent(in) :: oldu_val
      type(scalar_field), intent(in) :: density
      type(vector_field), intent(in) :: absorption
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: big_m_diag_addto
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)), intent(inout) :: big_m_tensor_addto
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
      type(vector_field), intent(inout) :: masslump
      type(petsc_csr_matrix), optional, intent(inout) :: mass
    
      integer :: dim
      real, dimension(ele_ngi(u, ele)) :: density_gi
      real, dimension(u%dim , ele_loc(u, ele)) :: absorption_lump
      real, dimension(u%dim, ele_ngi(u, ele)) :: absorption_gi
      real, dimension(u%dim, ele_loc(u, ele), ele_loc(u, ele)) :: absorption_mat
      
      density_gi=ele_val_at_quad(density, ele)
      absorption_gi = ele_val_at_quad(absorption, ele)
      
      ! element absorption matrix
      !  /
      !  | N_A N_B abs rho dV
      !  /
      absorption_mat = shape_shape_vector(test_function, ele_shape(u, ele), detwei*density_gi, absorption_gi)
      if(lump_absorption) then
        if(.not.abs_lump_on_submesh) then
          absorption_lump = sum(absorption_mat, 3)
          do dim = 1, u%dim
            big_m_diag_addto(dim, :) = big_m_diag_addto(dim, :) + dt*theta*absorption_lump(dim,:)
            rhs_addto(dim, :) = rhs_addto(dim, :) - absorption_lump(dim,:)*oldu_val(dim,:)
          end do
        end if
      else
        do dim = 1, u%dim
          big_m_tensor_addto(dim, dim, :, :) = big_m_tensor_addto(dim, dim, :, :) + &
            & dt*theta*absorption_mat(dim,:,:)
          rhs_addto(dim, :) = rhs_addto(dim, :) - matmul(absorption_mat(dim,:,:), oldu_val(dim,:))
        end do
        absorption_lump = 0.0
      end if
      if (pressure_corrected_absorption) then
        if (assemble_inverse_masslump.and.(.not.(abs_lump_on_submesh))) then
          call addto(masslump, ele_nodes(u, ele), theta*absorption_lump)
        end if
        if (assemble_mass_matrix) then
          do dim = 1, u%dim
            call addto(mass, dim, dim, ele_nodes(u, ele), ele_nodes(u,ele), &
               theta*absorption_mat(dim,:,:))
          end do
        end if
      end if
      
    end subroutine add_absorption_element_cg
      
    subroutine add_viscosity_element_cg(ele, u, oldu_val, nu, x, viscosity, grad_u, &
         du_t, detwei, big_m_tensor_addto, rhs_addto)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: u, nu
      real, dimension(:,:), intent(in) :: oldu_val
      type(vector_field), intent(in) :: x
      type(tensor_field), intent(in) :: viscosity
      type(tensor_field), intent(in) :: grad_u
      real, dimension(ele_loc(u, ele), ele_ngi(u, ele), u%dim), intent(in) :: du_t
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)), intent(inout) :: big_m_tensor_addto
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
    
      integer :: dim, dimj, gi, iloc
      real, dimension(u%dim, u%dim, ele_ngi(u, ele)) :: viscosity_gi
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)) :: viscosity_mat
      real, dimension(x%dim, x%dim, ele_ngi(u,ele)) :: les_tensor_gi
      real, dimension(ele_ngi(u, ele)) :: les_coef_gi
      real, dimension(x%dim, ele_loc(u,ele), ele_loc(u,ele)) :: div_les_viscosity
      real, dimension(x%dim, x%dim, ele_loc(u,ele)) :: grad_u_nodes
      
      if (have_viscosity) then
         viscosity_gi = ele_val_at_quad(viscosity, ele)
      else
         ! if we don't have viscosity but maybe LES
         viscosity_gi = 0.0
      end if
      
      if (have_les) then
         ! add in LES viscosity
         les_tensor_gi=les_length_scale_tensor(du_t, ele_shape(u, ele))
         les_coef_gi=les_viscosity_strength(du_t, ele_val(nu, ele))
         do gi=1, size(les_coef_gi)
            les_tensor_gi(:,:,gi)=les_coef_gi(gi)*les_tensor_gi(:,:,gi)* &
                 smagorinsky_coefficient**2
         end do
         if (les_fourth_order) then
           div_les_viscosity=dshape_dot_tensor_shape(du_t, les_tensor_gi, ele_shape(u, ele), detwei)
           grad_u_nodes=ele_val(grad_u, ele)
           do dim=1, u%dim
             do iloc=1, ele_loc(u, ele)
               rhs_addto(dim,iloc)=rhs_addto(dim,iloc)+ &
                 sum(div_les_viscosity(:,:,iloc)*grad_u_nodes(:,dim,:))
             end do
           end do
         end if
         viscosity_gi=viscosity_gi+les_tensor_gi
      end if
      
      ! element viscosity matrix - tensor form
      !  /
      !  | gradN_A^T viscosity gradN_B dV
      !  /
      ! only valid when incompressible and viscosity tensor is isotropic
      viscosity_mat = 0.0
      if(stress_form) then
        ! add in the stress form entries of the element viscosity matrix
        !  /
        !  | B_A^T C B_B dV
        !  /
        viscosity_mat = stiffness_matrix(du_t, viscosity_gi, du_t, detwei, legacy_stress)
      else
        if(isotropic_viscosity .and. .not. have_les) then
          assert(u%dim > 0)
          viscosity_mat(1, 1, :, :) = dshape_dot_dshape(du_t, du_t, detwei * viscosity_gi(1, 1, :))
          do dim = 2, u%dim
            viscosity_mat(dim, dim, :, :) = viscosity_mat(1, 1, :, :)
          end do
        else if(diagonal_viscosity .and. .not. have_les) then
          assert(u%dim > 0)
          viscosity_mat(1, 1, :, :) = dshape_diagtensor_dshape(du_t, viscosity_gi, du_t, detwei)
          do dim = 2, u%dim
            viscosity_mat(dim, dim, :, :) = viscosity_mat(1, 1, :, :)
          end do
        else
          do dim = 1, u%dim
            viscosity_mat(dim, dim, :, :) = &
                        dshape_tensor_dshape(du_t, viscosity_gi, du_t, detwei)
          end do
        end if
      end if
      
      big_m_tensor_addto = big_m_tensor_addto + dt*theta*viscosity_mat
      
      do dim = 1, u%dim
        rhs_addto(dim, :) = rhs_addto(dim, :) - matmul(viscosity_mat(dim,dim,:,:), oldu_val(dim,:))
      
        ! off block diagonal viscosity terms
        if(stress_form) then
          do dimj = 1, u%dim

            if (dim==dimj) cycle ! already done this

            rhs_addto(dim, :) = rhs_addto(dim, :) - matmul(viscosity_mat(dim,dimj,:,:), oldu_val(dimj,:))
          end do
        end if
      end do
      
    end subroutine add_viscosity_element_cg
    
    subroutine add_coriolis_element_cg(ele, test_function, x, u, oldu_val, density, detwei, big_m_tensor_addto, rhs_addto)
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: x
      type(vector_field), intent(in) :: u
      real, dimension(:,:), intent(in) :: oldu_val
      type(scalar_field), intent(in) :: density
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei
      real, dimension(u%dim, u%dim, ele_loc(u, ele), ele_loc(u, ele)), intent(inout) :: big_m_tensor_addto
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
      
      real, dimension(ele_ngi(u, ele)) :: coriolis_gi, density_gi
      real, dimension(ele_loc(u, ele), ele_loc(u, ele)) :: coriolis_mat
      
      density_gi = ele_val_at_quad(density, ele)
      
      ! element coriolis matrix
      !  /
      !  | N_A N_B rho omega dV
      !  /
      !
      ! scaling factor (omega, or f_0+\beta y, etc. depending on options):
      coriolis_gi=coriolis(ele_val_at_quad(x,ele))
      coriolis_mat = shape_shape(test_function, ele_shape(u, ele), density_gi*coriolis_gi*detwei)
      
      ! cross terms in U_ and V_ for coriolis
      big_m_tensor_addto(U_, V_, :, :) = big_m_tensor_addto(U_, V_, :, :) - dt*theta*coriolis_mat
      big_m_tensor_addto(V_, U_, :, :) = big_m_tensor_addto(V_, U_, :, :) + dt*theta*coriolis_mat
      
      rhs_addto(U_, :) = rhs_addto(U_, :) + matmul(coriolis_mat, oldu_val(V_,:))
      rhs_addto(V_, :) = rhs_addto(V_, :) - matmul(coriolis_mat, oldu_val(U_,:))
      
    end subroutine add_coriolis_element_cg
    
    subroutine add_geostrophic_pressure_element_cg(ele, test_function, x, u, gp,  detwei, rhs_addto)
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function
      type(vector_field), intent(in) :: x
      type(vector_field), intent(in) :: u
      type(scalar_field), intent(in) :: gp
      real, dimension(ele_ngi(u, ele)), intent(in) :: detwei
      real, dimension(u%dim, ele_loc(u, ele)), intent(inout) :: rhs_addto
            
      real, dimension(ele_loc(gp, ele), ele_ngi(gp, ele), mesh_dim(gp)) :: dgp_t
      
      ! We assume here that gp is usually on a different mesh to u or p
      call transform_to_physical(x, ele, ele_shape(gp, ele), &
        & dshape = dgp_t)
        
      rhs_addto = rhs_addto - shape_vector_rhs(test_function, transpose(ele_grad_at_quad(gp, ele, dgp_t)), detwei)
      
    end subroutine add_geostrophic_pressure_element_cg
    
    function stiffness_matrix(dshape1, tensor, dshape2, detwei, legacy) result (matrix)
      !!< Calculates the stiffness matrix.
      !!< 
      !!<          /
      !!< matrix = | b_a^T c b_b dV
      !!<          /
      !!<
      !!< where
      !!< b_a = / N_a,x   0     0   \   c = /  4/3*mu_xx  -2/3*mu_xy -2/3*mu_xz  0    0    0   \
      !!<       |  0    N_a,y   0   |       | -2/3*mu_yx   4/3*mu_yy -2/3*mu_yz  0    0    0   |
      !!<       |   0     0   N_a,z |       | -2/3*mu_zx  -2/3*mu_zy  4/3*mu_zz  0    0    0   |
      !!<       | N_a,y N_a,x   0   |       |     0           0          0     mu_xy  0    0   |
      !!<       | N_a,z   0   N_a,x |       |     0           0          0       0  mu_xz  0   |
      !!<       \   0   N_a,z N_a,y /       \     0           0          0       0    0  mu_yz /
      !!< which results in:
      !!< b_a^T c b_b - I gradN_a^T diag(mu) gradN_b =
      !!<               /  N_a,x*N_b,x*mu_xx - 2/3*N_a,x*N_b,x*mu_xx + (N_a,x*N_b,x*mu_xx + N_a,y*N_b,y*mu_xy + N_a,z*N_b,z*mu_xz)
      !!<              |   N_a,x*N_b,y*mu_xy - 2/3*N_a,y*N_b,x*mu_yx   ...
      !!<              \   N_a,x*N_b,z*mu_xz - 2/3*N_a,z*N_b,x*mu_zx
      !!<
      !!<                  N_a,y*N_b,x*mu_xy - 2/3*N_a,x*N_b,y*mu_xy
      !!<              ... N_a,y*N_b,y*mu_yy - 2/3*N_a,y*N_b,y*mu_yy + (N_a,x*N_b,x*mu_xy + N_a,y*N_b,y*mu_yy + N_a,z*N_b,z*mu_yz) ...
      !!<                  N_a,y*N_b,z*mu_yz - 2/3*N_a,z*N_b,y*mu_zy
      !!<
      !!<                  N_a,z*N_b,x*mu_xz - 2/3*N_a,x*N_b,z*mu_xz                                                               \
      !!<              ... N_a,z*N_b,y*mu_yz - 2/3*N_a,y*N_b,z*mu_yz                                                               |
      !!<                  N_a,z*N_b,z*mu_zz - 2/3*N_a,z*N_b,z*mu_zz + (N_a,x*N_b,x*mu_xz + N_a,y*N_b,y*mu_yz + N_a,z*N_b,z*mu_zz) /
      !!< where the terms in brackets correspond to the tensor form entries I gradN_a^T row(symm(mu)) gradN_b (see below).
      !!<
      !!< The optional switch legacy controls whether the assembly follows the above multiplication
      !!< form or uses the legacy version from subroutine diff3d, which cannot be assembled from
      !!< any variation on the multiplication above.
      !!<
      !!< legacy implementation:
      !!< b_a^T c b_b - I gradN_a^T diag(mu) gradN_b =
      !!<               /  N_a,x*N_b,x*mu_xx - 2/3*N_a,x*N_b,x*mu_xx + (N_a,x*N_b,x*mu_xx + N_a,y*N_b,y*mu_yy + N_a,z*N_b,z*mu_zz)
      !!<              |   N_a,x*N_b,y*mu_yx - 2/3*N_a,y*N_b,x*mu_yy   ...
      !!<              \   N_a,x*N_b,z*mu_zx - 2/3*N_a,z*N_b,x*mu_zz
      !!<
      !!<                  N_a,y*N_b,x*mu_xy - 2/3*N_a,x*N_b,y*mu_xx
      !!<              ... N_a,y*N_b,y*mu_yy - 2/3*N_a,y*N_b,y*mu_yy + (N_a,x*N_b,x*mu_xx + N_a,y*N_b,y*mu_yy + N_a,z*N_b,z*mu_zz) ...
      !!<                  N_a,y*N_b,z*mu_zy - 2/3*N_a,z*N_b,y*mu_zz
      !!<
      !!<                  N_a,z*N_b,x*mu_xz - 2/3*N_a,x*N_b,z*mu_xx                                                               \
      !!<              ... N_a,z*N_b,y*mu_yz - 2/3*N_a,y*N_b,z*mu_yy                                                               |
      !!<                  N_a,z*N_b,z*mu_zz - 2/3*N_a,z*N_b,z*mu_zz + (N_a,x*N_b,x*mu_xx + N_a,y*N_b,y*mu_yy + N_a,z*N_b,z*mu_zz) /
      !!< where the terms in brackets correspond to the tensor form entries I gradN_a^T diag(mu) gradN_b (see below).

      real, dimension(:,:,:), intent(in) :: dshape1, dshape2
      real, dimension(size(dshape1,3),size(dshape1,3),size(dshape1,2)), intent(in) :: tensor
      real, dimension(size(dshape1,2)), intent(in) :: detwei
      logical, intent(in), optional :: legacy

      real, dimension(size(dshape1,3),size(dshape1,3),size(dshape1,1),size(dshape2,1)) :: matrix

      real, dimension(size(dshape1,3),size(dshape1,2)) :: tensor_diag, tensor_entries

      integer :: iloc,jloc, gi, i, j
      integer :: loc1, loc2, ngi, dim

      logical :: l_legacy

      if(present(legacy)) then
        l_legacy = legacy
      else
        l_legacy = .false.
      end if

      loc1=size(dshape1,1)
      loc2=size(dshape2,1)
      ngi=size(dshape1,2)
      dim=size(dshape1,3)

      assert(loc1==loc2)

      tensor_diag = 0.0
      tensor_entries = 0.0

      matrix=0.0
      if(l_legacy) then

        !            /
        ! matrix = I| gradN_a^T offdiag(mu) gradN_b dV
        !           /
        do i=1,dim
          matrix(i,i,:,:) = dshape_diagtensor_dshape(dshape1, tensor, dshape2, detwei)
        end do
        
        forall(i=1:dim)
          tensor_diag(i,:) = tensor(i,i,:)
        end forall

        ! matrix = matrix +  b_a^T c b_b - I gradN_a^T diag(mu) gradN_b =
        !          matrix +  /  N_a,x*N_b,x*mu_xx - 2/3*N_a,x*N_b,x*mu_xx
        !                    |   N_a,x*N_b,y*mu_yx - 2/3*N_a,y*N_b,x*mu_yy   ...
        !                    \   N_a,x*N_b,z*mu_zx - 2/3*N_a,z*N_b,x*mu_zz
        !
        !                        N_a,y*N_b,x*mu_xy - 2/3*N_a,x*N_b,y*mu_xx
        !                    ... N_a,y*N_b,y*mu_yy - 2/3*N_a,y*N_b,y*mu_yy   ...
        !                        N_a,y*N_b,z*mu_zy - 2/3*N_a,z*N_b,y*mu_zz
        !
        !                        N_a,z*N_b,x*mu_xz - 2/3*N_a,x*N_b,z*mu_xx   \
        !                    ... N_a,z*N_b,y*mu_yz - 2/3*N_a,y*N_b,z*mu_yy   |
        !                        N_a,z*N_b,z*mu_zz - 2/3*N_a,z*N_b,z*mu_zz  /
        do gi=1,ngi
          forall(iloc=1:loc1,jloc=1:loc2)
              matrix(:,:,iloc,jloc) = matrix(:,:,iloc,jloc) &
                                      +(spread(dshape1(iloc,gi,:), 1, dim) &
                                       *spread(dshape2(jloc,gi,:), 2, dim) &
                                       *tensor(:,:,gi) &
                                       -spread(dshape1(iloc,gi,:), 2, dim) &
                                       *spread(dshape2(jloc,gi,:), 1, dim) &
                                       *(2./3.)*spread(tensor_diag(:,gi), 2, dim)) &
                                      *detwei(gi)
          end forall
        end do

      else

        !            /
        ! matrix = I| gradN_a^T row(symm(mu)) gradN_b dV
        !           /
        do i=1,dim
          ! extract the relevent tensor entries into a vector
          do j = 1, i-1
            tensor_entries(j,:) = tensor(j,i,:)
          end do
          do j = i, dim
            tensor_entries(j,:) = tensor(i,j,:)
          end do
          matrix(i,i,:,:) = dshape_vector_dshape(dshape1, tensor_entries, dshape2, detwei)
        end do

        ! matrix = matrix + b_a^T c b_b - I gradN_a^T row(symm(mu)) gradN_b =
        !          matrix +  /  N_a,x*N_b,x*mu_xx - 2/3*N_a,x*N_b,x*mu_xx
        !                    |   N_a,x*N_b,y*mu_xy - 2/3*N_a,y*N_b,x*mu_yx   ...
        !                    \   N_a,x*N_b,z*mu_xz - 2/3*N_a,z*N_b,x*mu_zx
        !
        !                        N_a,y*N_b,x*mu_xy - 2/3*N_a,x*N_b,y*mu_xy
        !                    ... N_a,y*N_b,y*mu_yy - 2/3*N_a,y*N_b,y*mu_yy   ...
        !                        N_a,y*N_b,z*mu_yz - 2/3*N_a,z*N_b,y*mu_zy
        !
        !                        N_a,z*N_b,x*mu_xz - 2/3*N_a,x*N_b,z*mu_xz   \
        !                    ... N_a,z*N_b,y*mu_yz - 2/3*N_a,y*N_b,z*mu_yz   |
        !                        N_a,z*N_b,z*mu_zz - 2/3*N_a,z*N_b,z*mu_zz  /
        do gi=1,ngi
          forall(iloc=1:loc1,jloc=1:loc2)
              matrix(:,:,iloc,jloc) = matrix(:,:,iloc,jloc) &
                                      +(spread(dshape1(iloc,gi,:), 1, dim) &
                                       *spread(dshape2(jloc,gi,:), 2, dim) &
                                       *tensor(:,:,gi) &
                                       -spread(dshape1(iloc,gi,:), 2, dim) &
                                       *spread(dshape2(jloc,gi,:), 1, dim) &
                                       *(2./3.)*tensor(:,:,gi)) &
                                      *detwei(gi)
          end forall
        end do

      end if

    end function stiffness_matrix
    
    subroutine deallocate_cg_mass(mass, inverse_masslump)
      !!< Deallocates mass and/or inverse_masslump
      !!< if they are assembled in construct_momentum_cg()
      type(petsc_csr_matrix), intent(inout):: mass
      type(vector_field), intent(inout):: inverse_masslump
      
      if (assemble_mass_matrix) then
        call deallocate(mass)
      end if
      if (assemble_inverse_masslump) then
        call deallocate(inverse_masslump)
      end if
      
    end subroutine deallocate_cg_mass

    subroutine correct_masslumped_velocity(u, inverse_masslump, ct_m, delta_p)
      !!< Given the pressure correction delta_p, correct the velocity.
      !!<
      !!< U_new = U_old + M_l^{-1} * C * delta_P
      type(vector_field), intent(inout) :: u
      type(vector_field), intent(inout) :: inverse_masslump
      type(block_csr_matrix), intent(in) :: ct_m
      type(scalar_field), intent(in) :: delta_p

      ! Correction to u one dimension at a time.
      type(scalar_field) :: delta_u, inverse_masslump_component

      integer :: dim, i

      ewrite(1,*) 'correct_masslumped_velocity'

      call allocate(delta_u, u%mesh, "Delta_U")

      do dim=1,u%dim
        call mult_t(delta_u, block(ct_m,1,dim), delta_p)
        inverse_masslump_component = extract_scalar_field(inverse_masslump, dim)

        call scale(delta_u, inverse_masslump_component)
        call addto(u, dim, delta_u)
      end do

      call halo_update(u)
      do i = 1, u%dim
        ewrite_minmax(u%val(i)%ptr(:))
      end do

      call deallocate(delta_u)

    end subroutine correct_masslumped_velocity

    subroutine correct_velocity_cg(u, mass, ct_m, delta_p)
      !!< Given the pressure correction delta_p, correct the velocity.
      !!<
      !!< U_new = U_old + M_l^{-1} * C * delta_P
      type(vector_field), intent(inout) :: u
      type(petsc_csr_matrix), intent(inout) :: mass
      type(block_csr_matrix), intent(in) :: ct_m
      type(scalar_field), intent(in) :: delta_p

      ! Correction to u one dimension at a time.
      type(vector_field) :: delta_u1, delta_u2

      integer :: i

      ewrite(1,*) 'correct_velocity_cg'

      call allocate(delta_u1, u%dim, u%mesh, "Delta_U1")
      call allocate(delta_u2, u%dim, u%mesh, "Delta_U2")
      delta_u2%option_path = trim(delta_p%option_path)//&
                                  "/prognostic/scheme/use_projection_method&
                                  &/full_schur_complement/inner_matrix[0]"
      
      ! compute delta_u1=grad delta_p
      call mult_t(delta_u1, ct_m, delta_p)
      
      ! compute M^{-1} delta_u1
      call zero(delta_u2)
      call petsc_solve(delta_u2, mass, delta_u1)
      
      call addto(u, delta_u2)
      
      call halo_update(u)
      do i = 1, u%dim
        ewrite_minmax(u%val(i)%ptr(:))
      end do

      call deallocate(delta_U1)
      call deallocate(delta_U2)

    end subroutine correct_velocity_cg

    subroutine assemble_masslumped_poisson_rhs(poisson_rhs, &
      ctp_m, mom_rhs, ct_rhs, inverse_masslump, velocity, dt, theta_pg)

      type(scalar_field), intent(inout) :: poisson_rhs
      type(block_csr_matrix), intent(in) :: ctp_m
      type(vector_field), intent(in) :: mom_rhs
      type(scalar_field), intent(in) :: ct_rhs
      type(vector_field), intent(in) :: inverse_masslump
      type(vector_field), intent(in) :: velocity
      real, intent(in) :: dt, theta_pg

      type(vector_field) :: l_mom_rhs

      ewrite(1,*) 'Entering assemble_masslumped_poisson_rhs'

      call allocate(l_mom_rhs, mom_rhs%dim, mom_rhs%mesh, name="AssemblePoissonMomRHS")
      
      ! poisson_rhs = ct_rhs/dt - C^T ( M_L^-1 mom_rhs + velocity/dt )
      
      ! compute M_L^-1 mom_rhs + velocity/dt
      call set(l_mom_rhs, mom_rhs)
      call scale(l_mom_rhs, inverse_masslump)
      call addto(l_mom_rhs, velocity, scale=1.0/dt/theta_pg)
      
      ! need to update before the mult, as halo of mom_rhs may not be valid
      ! (although it probably is in halo 1 - let's be safe anyway)
      call halo_update(l_mom_rhs)

      call mult(poisson_rhs, ctp_m, l_mom_rhs)
      call scale(poisson_rhs, -1.0)

      call addto(poisson_rhs, ct_rhs, scale=1.0/dt/theta_pg)

      call deallocate(l_mom_rhs)

    end subroutine assemble_masslumped_poisson_rhs

    subroutine assemble_kmk_matrix(state, pressure_mesh, coordinates, &
      theta_pg)
    ! Assemble P1-P1 stabilisation term in the pressure matrix.
      type(state_type), intent(inout) :: state
      type(mesh_type), intent(inout) :: pressure_mesh
      type(vector_field), intent(in) :: coordinates
      ! the required term is K^T M^-1 K (theta dt dp), the variable we're
      ! solving for in the pressure equation however is theta**2 dt dp
      ! thus we have to divide kmk by theta
      real, intent(in) :: theta_pg

      type(csr_matrix), pointer :: kmk  
      type(csr_sparsity), pointer :: p_sparsity

      integer :: ele
      type(csr_matrix) :: kt
      real, dimension(mesh_dim(pressure_mesh), mesh_dim(pressure_mesh), ele_ngi(pressure_mesh, 1)) :: h_bar
      real, dimension(mesh_dim(pressure_mesh), mesh_dim(pressure_mesh)) :: ele_tensor
      type(element_type), pointer :: p_shape
      real, dimension(ele_ngi(pressure_mesh, 1)) :: detwei
      real, dimension(ele_loc(pressure_mesh, 1), ele_ngi(pressure_mesh, 1), coordinates%dim) :: dp_t
      real, dimension(ele_loc(pressure_mesh, 1), ele_loc(pressure_mesh, 1)) :: little_stiff_matrix
      type(scalar_field) :: scaled_p_masslump
      type(scalar_field), pointer :: p_masslump

      p_shape => ele_shape(pressure_mesh, 1)

      kmk => get_pressure_stabilisation_matrix(state)
      
      p_sparsity => get_csr_sparsity_firstorder(state, pressure_mesh, pressure_mesh)
      call allocate(kt, p_sparsity, name="PressureDiffusionMatrix")
      call zero(kt)
      p_masslump => get_lumped_mass(state, pressure_mesh)

      ! Assemble the pressure diffusion matrix k. The diffusion parameter is
      ! given by a tensor describing the element length scales in physical space
      ! (h_bar). Simplex_tensor gives the metric that would make that element
      ! the ideal element.
      do ele=1,ele_count(pressure_mesh)
        call transform_to_physical(coordinates, ele, p_shape, dshape=dp_t, detwei=detwei)
        ele_tensor = simplex_tensor(coordinates, ele)
        h_bar = spread(edge_length_from_eigenvalue(ele_tensor), 3, size(h_bar, 3))
        little_stiff_matrix = dshape_tensor_dshape(dp_t, h_bar, dp_t, detwei)
        call addto(kt, ele_nodes(pressure_mesh, ele), ele_nodes(pressure_mesh, ele), 0.5 * little_stiff_matrix)
      end do
        
      ! by scaling masslump with theta, we divide kmk by theta
      if(abs(theta_pg - 1.0) < epsilon(0.0)) then
        call mult_div_invscalar_div_T(kmk, kt, p_masslump, kt)
      else
        call allocate(scaled_p_masslump, p_masslump%mesh, trim(p_masslump%name) // "Scaled")
        call set(scaled_p_masslump, p_masslump)
        call scale(scaled_p_masslump, theta_pg)
      
        ! Compute kmk, the stabilisation term.
        call mult_div_invscalar_div_T(kmk, kt, scaled_p_masslump, kt)
        
        call deallocate(scaled_p_masslump)
      end if
      call deallocate(kt)
      
    end subroutine assemble_kmk_matrix

    subroutine add_kmk_matrix(state, cmc_m)
    ! Add kmk (P1-P1 stabilisation term in the pressure matrix) to cmc_m.
      type(state_type), intent(inout) :: state
      type(csr_matrix), intent(inout) :: cmc_m
      type(csr_matrix), pointer :: kmk

      kmk => get_pressure_stabilisation_matrix(state)
      call addto(cmc_m, kmk)

    end subroutine add_kmk_matrix

    subroutine add_kmk_rhs(state, rhs, pressure, dt)
      type(state_type), intent(inout) :: state
      type(scalar_field), intent(inout) :: rhs
      type(scalar_field), intent(in) :: pressure
      real, intent(in) :: dt

      type(csr_matrix), pointer :: kmk

      kmk => get_pressure_stabilisation_matrix(state)
      call mult(rhs, kmk, pressure)
      call scale(rhs, dt)
    end subroutine add_kmk_rhs

  end module momentum_cg
