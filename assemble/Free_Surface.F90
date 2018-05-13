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

module free_surface_module
use fldebug
use integer_set_module
use data_structures
use spud
use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
use futils, only: int2str, present_and_true
use parallel_tools
use sparse_tools
use parallel_fields
use eventcounter
use cv_faces
use transform_elements
use fetools
use fields
use sparse_tools_petsc
use state_module
use sparse_matrices_fields
use boundary_conditions
use vertical_extrapolation_module
use halos
use field_options
use physics_from_options
use tidal_module, only: calculate_diagnostic_equilibrium_pressure
use sparsity_patterns
use sparsity_patterns_meshes
use solvers
use cv_shape_functions
implicit none

private

public move_mesh_free_surface, add_free_surface_to_cmc_projection, &
  vertical_prolongator_from_free_surface, &
  free_surface_nodes, calculate_diagnostic_free_surface, &
  add_free_surface_to_poisson_rhs, copy_poisson_solution_to_interior, &
  calculate_diagnostic_wettingdrying_alpha, insert_original_distance_to_bottom, &
  calculate_volume_by_surface_integral
public get_extended_pressure_mesh_for_viscous_free_surface, copy_to_extended_p, &
  get_extended_velocity_divergence_matrix, get_extended_pressure_poisson_matrix, &
  get_extended_schur_auxillary_sparsity, &
  update_pressure_and_viscous_free_surface,  &
  add_implicit_viscous_free_surface_integrals, &
  add_implicit_viscous_free_surface_integrals_cv, &
  add_implicit_viscous_free_surface_scaled_mass_integrals, update_prognostic_free_surface, &
  update_implicit_scaled_free_surface, has_implicit_viscous_free_surface_bc, &
  has_explicit_viscous_free_surface_bc, has_standard_free_surface_bc, &
  add_explicit_viscous_free_surface_integrals, &
  add_explicit_viscous_free_surface_integrals_cv

public free_surface_module_check_options

contains

  function has_standard_free_surface_bc(u) result(standard_free_surface)
  !!< whether velocity has a 'standard' free surface boundary condition,
  !!< i.e. a fs bc without the no_normal_stress option
  type(vector_field), intent(in):: u
  logical:: standard_free_surface

   character(len=FIELD_NAME_LEN):: bctype
   character(len=OPTION_PATH_LEN) :: bc_option_path
   integer:: i    

   standard_free_surface = .false.
   do i=1, get_boundary_condition_count(u)
     call get_boundary_condition(u, i, type=bctype, &
          option_path=bc_option_path)
     if (bctype=="free_surface") then

       if(.not.have_option(trim(bc_option_path)//"/type[0]/no_normal_stress")) then
         standard_free_surface = .true.
         return
       end if
     end if
   end do

  end function has_standard_free_surface_bc

  function has_implicit_viscous_free_surface_bc(u) result(implicit_free_surface)
  !!< whether velocity has a free surface boundary condition with
  !!< the no_normal_stress option but without the explicit option
  type(vector_field), intent(in):: u
  logical:: implicit_free_surface

   character(len=FIELD_NAME_LEN):: bctype
   character(len=OPTION_PATH_LEN) :: bc_option_path
   integer:: i    

   implicit_free_surface = .false.
   do i=1, get_boundary_condition_count(u)
     call get_boundary_condition(u, i, type=bctype, &
          option_path=bc_option_path)
     if (bctype=="free_surface") then

       if(have_option(trim(bc_option_path)//"/type[0]/no_normal_stress") .and. &
          (.not.have_option(trim(bc_option_path)//"/type[0]/no_normal_stress/explicit"))) then
         implicit_free_surface = .true.
         return
       end if
     end if
   end do

  end function has_implicit_viscous_free_surface_bc

  function has_explicit_viscous_free_surface_bc(u) result(explicit_free_surface)
  !!< whether velocity has a free surface boundary condition with
  !!< the no_normal_stress option and the explicit option
  type(vector_field), intent(in):: u
  logical:: explicit_free_surface

   character(len=FIELD_NAME_LEN):: bctype
   character(len=OPTION_PATH_LEN) :: bc_option_path
   integer:: i    

   explicit_free_surface = .false.
   do i=1, get_boundary_condition_count(u)
     call get_boundary_condition(u, i, type=bctype, &
          option_path=bc_option_path)
     if (bctype=="free_surface") then

       if(have_option(trim(bc_option_path)//"/type[0]/no_normal_stress") .and. &
          have_option(trim(bc_option_path)//"/type[0]/no_normal_stress/explicit")) then
         explicit_free_surface = .true.
         return
       end if
     end if
   end do

  end function has_explicit_viscous_free_surface_bc

  subroutine update_implicit_scaled_free_surface(states)
  !!< Set OldScaledFreeSurface to ScaledFreeSurface, these
  !!< ScaledFreeSurface is the surface fields \Delta\rho g\eta that we solve for 
  !!< with the implicit viscous fs method
  type(state_type), dimension(:), intent(in) :: states

    type(scalar_field), pointer:: free_surface, scaled_fs, old_scaled_fs
    integer:: i, fs_stat

    do i=1,size(states)
      free_surface => extract_scalar_field(states(i), "FreeSurface", stat=fs_stat)
      if(fs_stat==0) then
        if(have_option(trim(free_surface%option_path)//"/prognostic")) then
          if (has_boundary_condition_name(free_surface, "_implicit_free_surface")) then
            ewrite(2,*) 'Updating OldScaledFreeSurface surface field for the FreeSurface field in state '//trim(states(i)%name)//'.' 

            scaled_fs => extract_surface_field(free_surface, "_implicit_free_surface", "ScaledFreeSurface")
            old_scaled_fs => extract_surface_field(free_surface, "_implicit_free_surface", "OldScaledFreeSurface")
            
            call set(old_scaled_fs, scaled_fs)

          end if
        end if
      end if
    end do

  end subroutine update_implicit_scaled_free_surface

  subroutine insert_original_distance_to_bottom(state)
    !!< Adds the OriginalDistanceToBottom field into the state. 
    !!< Note: In order to to get the correct values, this subroutine 
    !!< has to be called before the first timestep.
    type(state_type), intent(inout) :: state
    type(mesh_type), pointer :: p_mesh
    type(scalar_field), pointer :: bottomdist
    type(scalar_field) :: original_bottomdist, original_bottomdist_remap

    if (.not. has_scalar_field(state, "OriginalDistanceToBottom")) then
       ewrite(2, *) "Inserting OriginalDistanceToBottom field into state."   
       bottomdist => extract_scalar_field(state, "DistanceToBottom")
       call allocate(original_bottomdist, bottomdist%mesh, "OriginalDistanceToBottom")
       call zero(original_bottomdist)
       call addto(original_bottomdist, bottomdist)
       call insert(state, original_bottomdist, name="OriginalDistanceToBottom")
       call deallocate(original_bottomdist)

       ! We also cache  the OriginalDistanceToBottom on the pressure mesh
       ewrite(2, *) "Inserting OriginalDistanceToBottomPressureMesh field into state."   
       p_mesh => extract_pressure_mesh(state)
       call allocate(original_bottomdist_remap, p_mesh, "OriginalDistanceToBottomPressureMesh")
       call remap_field(original_bottomdist, original_bottomdist_remap)
       call insert(state, original_bottomdist_remap, name="OriginalDistanceToBottomPressureMesh")
       call deallocate(original_bottomdist_remap)
    end if
         
  end subroutine insert_original_distance_to_bottom
  
  subroutine add_free_surface_to_cmc_projection(state, cmc, dt, &
                                                theta_pressure_gradient, theta_divergence, &
                                                assemble_cmc, rhs)
  !!< Adds a boundary integral to the continuity equation
  !!< that weakly enforces the kinematic boundary condition.
  !!<
  !!< The free surface is combined in with the pressure field such that
  !!< *at* the free surface p=g rho0 \eta. This has the advantage of combining
  !!< the pressure gradient and free surface gradient terms in the momentum
  !!< equation. It solves the continuity equation directly coupled with
  !!< the free surface.
  !!< With this approach all pressures are considered at the integer time 
  !!< levels, i.e. we apply a theta weighting for the pressure gradient term
  !!< in the momentum equation:
  !!<   M (u^n+1-u^n) - dt C p^{n+theta_pressure_gradient} + ... = 0
  !!< We're solving the continuity equation:
  !!<   C^T u^{n+theta_divergence}+ M_fs p^{n+1}-p^n = 0
  !!< which leads to a projection equation of:
  !!<   ( C^T M^-1 C dt dp + coef M_fs ) phi = 
  !!<     theta_divergence C^T u* + (1-theta_divergence) C^T u^n - 
  !!<     alpha M_fs (p*-p^n)
  !!< where M_fs is the free surface integral of M_i M_j, 
  !!< alpha=1/(g dt), coef=alpha/(theta_divergence theta_pressure_gradient dt)
  !!< and phi=dp theta_divergence theta_pressure_gradient dt. Note however,
  !!< that dp in the routine stands for phi, see correct_pressure in 
  !!< Momentum_Equation.F90
  
    type(state_type), intent(inout) :: state
    type(csr_matrix), intent(inout) :: cmc
    real, intent(in) :: dt
    real, intent(in) :: theta_pressure_gradient
    real, intent(in) :: theta_divergence
    !! only add in to the matrix if assemble_cmc==.true.
    !! if .not. assemble_cmc we still need to add in a correction 
    !! if the timestep has changed since last call
    logical, intent(in):: assemble_cmc
    type(scalar_field), optional, intent(inout) :: rhs
    
      type(integer_hash_table):: sele_to_fs_ele
      type(vector_field), pointer:: positions, u, gravity_normal, old_positions
      type(scalar_field), pointer:: p, prevp
      type(scalar_field), pointer:: density, old_surface_density
      type(scalar_field), pointer:: free_surface, scaled_fs, old_scaled_fs
      type(scalar_field), pointer :: original_bottomdist_remap
      type(scalar_field), pointer :: external_density
      type(mesh_type), pointer:: surface_mesh, fs_mesh, embedded_fs_mesh
      type(mesh_type):: dens_surface_mesh
      integer, dimension(:), pointer:: dens_surface_node_list
      character(len=FIELD_NAME_LEN):: bctype
      character(len=OPTION_PATH_LEN) :: fs_option_path
      real:: g, rho0, delta_rho, alpha, alpha_old, coef, coef_old, d0
      integer, dimension(:), pointer:: surface_element_list, fs_surface_element_list
      integer:: i, j, grav_stat, dens_stat, external_density_stat
      logical:: include_normals, move_mesh
      logical:: addto_cmc, variable_density, any_variable_density, have_density
      logical:: have_wd, have_wd_node_int, have_external_density
      logical:: implicit_prognostic_fs, use_fs_mesh
      
      real, save :: dt_old = 0.0
      
      ewrite(1,*) 'Entering add_free_surface_to_cmc_projection routine'
      ewrite(2,*) "Are we adding free-surface contribution to RHS:",present(rhs)

      ! gravity acceleration
      call get_option('/physical_parameters/gravity/magnitude', g, stat=grav_stat)
      if (grav_stat/=0) then
        FLExit("For a free surface you need gravity")
      end if
      
      ! get the pressure, and the pressure at the beginning of the time step
      p => extract_scalar_field(state, "Pressure")
      prevp => extract_scalar_field(state, "OldPressure")
      u => extract_vector_field(state, "Velocity")

      ! the prognostic free surface can only be used with the no_normal_stress option
      implicit_prognostic_fs=has_implicit_viscous_free_surface_bc(u)
      if(implicit_prognostic_fs) then
        free_surface => extract_scalar_field(state, "FreeSurface")
        assert(have_option(trim(free_surface%option_path)//"/prognostic"))
        if (.not. has_boundary_condition_name(free_surface, "_implicit_free_surface")) then
          call initialise_implicit_prognostic_free_surface(state, free_surface, u)
        end if
        ! obtain the f.s. surface mesh that has been stored under the
        ! "_implicit_free_surface" boundary condition
        call get_boundary_condition(free_surface, "_implicit_free_surface", &
           surface_mesh=fs_mesh, surface_element_list=fs_surface_element_list)
        ! create a map from face numbers to element numbers in fs_mesh
        call invert_set(fs_surface_element_list, sele_to_fs_ele)
        scaled_fs => extract_surface_field(free_surface, "_implicit_free_surface", "ScaledFreeSurface")
        old_scaled_fs => extract_surface_field(free_surface, "_implicit_free_surface", "OldScaledFreeSurface")
        embedded_fs_mesh => extract_mesh(state, "_embedded_free_surface_mesh")
      end if
      
      ! reference density
      call get_fs_reference_density_from_options(rho0, state%option_path)
      density => extract_scalar_field(state, "Density", stat=dens_stat)
      have_density = (dens_stat==0)

      have_wd=have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying")
      have_wd_node_int=have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/conserve_geometric_volume")
      if (have_wd) then
        if (.not. assemble_cmc) then
             ewrite(-1,*) "Wetting and drying needs to be reassembled at each timestep at the moment. Switch it on "//&
                   &"in diamond under .../Pressure/prognostic/scheme/update_discretised_equation"
             FLExit("Error in user options")
        end if
        call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/d0", d0)
        original_bottomdist_remap=>extract_scalar_field(state, "OriginalDistanceToBottomPressureMesh")
      end if

      move_mesh = have_option("/mesh_adaptivity/mesh_movement/free_surface")
      ! only include the inner product of gravity and surface normal
      ! if the free surface nodes are actually moved (not necessary
      ! for large scale ocean simulations) - otherwise the free surface is assumed flat
      ! (and on top unless we use a prognostic fs, which may be on the bottom)
      include_normals = move_mesh .or. implicit_prognostic_fs
      if (include_normals) then
        ewrite(2,*) 'Including inner product of normals in kinematic bc'
        gravity_normal => extract_vector_field(state, "GravityDirection")
      end if
      
      ewrite_minmax(p)
      ewrite_minmax(prevp)
      if(present(rhs)) then
         ewrite_minmax(rhs)
      end if
      
      if (move_mesh) then
        positions => extract_vector_field(state, "IteratedCoordinate")
        old_positions => extract_vector_field(state, "OldCoordinate")
      else
        positions => extract_vector_field(state, "Coordinate")
      end if
      
      ! only add to cmc if we're reassembling it (assemble_cmc = .true.)
      ! or if the timestep has changed
      addto_cmc = assemble_cmc.or.&
        (have_option("/timestepping/adaptive_timestep").and.(dt/=dt_old))
      any_variable_density = .false.
      do i=1, get_boundary_condition_count(u)
        call get_boundary_condition(u, i, type=bctype, &
           surface_mesh=surface_mesh, &
           surface_element_list=surface_element_list, &
           option_path=fs_option_path)
        if (bctype=="free_surface" .and. size(surface_element_list)>0) then

          if (have_option(trim(fs_option_path)//"/type[0]/no_normal_stress").and. &
              (.not.have_option(trim(fs_option_path)//"/type[0]/no_normal_stress/explicit"))) then
            assert(implicit_prognostic_fs)
            use_fs_mesh=.true.
          else
            use_fs_mesh=.false.
          end if

          variable_density = have_option(trim(fs_option_path)//"/type[0]/variable_density") &
                            .and. (.not.move_mesh)
          if (variable_density .and. (.not. have_density)) then
            FLExit("Variable density free surface requires a Density field.")
          end if
          if (variable_density) then
            if(.not.has_scalar_surface_field(u, i, "OldSurfaceDensity")) then

              if(.not.assemble_cmc) then
                FLAbort("Adding free surface to cmc with variable density but no old surface density field available.")
              end if

              call create_surface_mesh(dens_surface_mesh, dens_surface_node_list, &
                                       density%mesh, surface_element_list, "FreeSurfaceDensityMesh")

              allocate(old_surface_density)
              call allocate(old_surface_density, dens_surface_mesh, "OldSurfaceDensity")
              call deallocate(dens_surface_mesh)

              ! shouldn't be used this time around anyway, thanks to coef_old=0.0
              call remap_field_to_surface(density, old_surface_density, surface_element_list)

              call insert_surface_field(u, i, old_surface_density)

              call deallocate(old_surface_density)
              deallocate(old_surface_density)
            end if
            old_surface_density => extract_scalar_surface_field(u, i, "OldSurfaceDensity")

          end if
          any_variable_density = any_variable_density .or. variable_density

          external_density => extract_scalar_surface_field(u, i, "ExternalDensity", stat=external_density_stat)
          have_external_density = external_density_stat==0
          ! delta_rho is only used when .not. variable_density:
          if (have_external_density) then
            delta_rho = rho0 - node_val(external_density, 1)
          else
            delta_rho = rho0
          end if

          alpha=1.0/g/dt                        ! delta_rho included in alpha and coeff within element loop
          coef = alpha/(theta_pressure_gradient*theta_divergence*dt)

          if (assemble_cmc) then
            coef_old=0.0
          else
            alpha_old=1.0/g/dt_old
            coef_old = alpha_old/(theta_pressure_gradient*theta_divergence*dt_old)
          end if

          do j=1, size(surface_element_list)
            call add_free_surface_element(j, surface_element_list(j))
          end do

          if(variable_density) then
            ! save the current density at the surface for next time
            call remap_field_to_surface(density, old_surface_density, surface_element_list)
          end if

        end if
      end do
      if(addto_cmc.or.any_variable_density) then
        ! cmc has been modified (most likely by changing the timestep)
        ! therefore we need to invalidate the solver context
        call destroy_solver_cache(cmc)
      end if
      if (implicit_prognostic_fs) then
        call deallocate(sele_to_fs_ele)
      end if      
      
      ! save the current timestep, so we know how much to add when it changes
      dt_old = dt
    
      if(present(rhs)) then
         ewrite_minmax(rhs)
      end if
      
    contains 
    
    subroutine add_free_surface_element(surface_mesh_ele, sele)
      integer, intent(in):: surface_mesh_ele, sele ! element number in surface mesh and facet number in full mesh
      
      integer, dimension(face_loc(p, sele)):: nodes
      integer :: i

      real, dimension(face_loc(p, sele)):: top_pressures, old_top_pressures
      real, dimension(positions%dim, face_ngi(positions, sele)):: normals
      real, dimension(face_loc(p, sele), face_loc(p, sele)):: mass_ele, mass_ele_wd, mass_ele_old, mass_ele_old_wd
      real, dimension(face_ngi(p, sele)):: detwei, alpha_wetdry_quad, alpha_wetdry_quad_prevp
      real, dimension(face_loc(p, sele)):: alpha_wetdry, alpha_wetdry_prevp

      real, dimension(face_ngi(p, sele)):: inv_delta_rho_quad, old_inv_delta_rho_quad

      if(include_normals) then
        call transform_facet_to_physical(positions, sele, detwei_f=detwei,&
            & normal=normals)
        ! at each gauss point multiply with inner product of gravity and surface normal
        detwei=detwei*(-1.0)*sum(face_val_at_quad(gravity_normal,sele)*normals, dim=1)
      else
        call transform_facet_to_physical(positions, sele, detwei_f=detwei)
      end if
      
      if(variable_density) then
        if (have_external_density) then
          inv_delta_rho_quad = 1.0/(face_val_at_quad(density, sele) - ele_val_at_quad(external_density, surface_mesh_ele))
          old_inv_delta_rho_quad = 1.0/(ele_val_at_quad(old_surface_density, surface_mesh_ele) - &
            ele_val_at_quad(external_density, surface_mesh_ele))
        else
          inv_delta_rho_quad = 1.0/face_val_at_quad(density, sele)
          old_inv_delta_rho_quad = 1.0/ele_val_at_quad(old_surface_density, surface_mesh_ele)
        end if
      end if
      
      if (have_wd) then
        if (have_wd_node_int) then
           call compute_alpha_wetdry(p, sele, alpha_wetdry)
           call compute_alpha_wetdry(prevp, sele, alpha_wetdry_prevp)
        else
           call compute_alpha_wetdry_quad(p, sele, alpha_wetdry_quad)
           call compute_alpha_wetdry_quad(prevp, sele, alpha_wetdry_quad_prevp)
        end if
      end if
      
      if (have_wd .and. .not. have_wd_node_int) then
        mass_ele=shape_shape(face_shape(p, sele), face_shape(p, sele), detwei*(1.0-alpha_wetdry_quad))
        mass_ele_wd=shape_shape(face_shape(p, sele), face_shape(p, sele), detwei*alpha_wetdry_quad)
      else
        if (variable_density) then   ! this excludes the case of a moving mesh so we use mass_ele_old here too
          mass_ele=shape_shape(face_shape(p, sele), face_shape(p, sele), detwei*inv_delta_rho_quad)
          mass_ele_old=shape_shape(face_shape(p, sele), face_shape(p, sele), detwei*old_inv_delta_rho_quad)
        else
          mass_ele=shape_shape(face_shape(p, sele), face_shape(p, sele), detwei)
        end if
        if (have_wd .and. have_wd_node_int) then
          mass_ele_wd=mass_ele
          do i=1,size(mass_ele,1)
            mass_ele(i,:)=mass_ele(i,:)*(1.0-alpha_wetdry)
            mass_ele_wd(i,:)=mass_ele_wd(i,:)*alpha_wetdry
          end do
        end if
      end if
      
      if (use_fs_mesh) then
        nodes = ele_nodes(embedded_fs_mesh, fetch(sele_to_fs_ele, sele))
        top_pressures = ele_val(scaled_fs, fetch(sele_to_fs_ele, sele))
        old_top_pressures = ele_val(old_scaled_fs, fetch(sele_to_fs_ele, sele))
      else
        nodes = face_global_nodes(p, sele)
        top_pressures = face_val(p, sele)
        old_top_pressures = face_val(prevp, sele)
      end if

      if (addto_cmc.or.variable_density) then
        ! we consider the projection equation to solve for 
        ! phi=theta_pressure_gradient theta_divergence dt dp, so that the f.s. integral
        ! alpha M_fs dp=alpha M_fs phi/(theta_pressure_gradient theta_divergence g dt**2)
        !              =coef M_fs phi
        if(variable_density) then
          call addto(cmc, nodes, nodes, &
            (coef*mass_ele - coef_old*mass_ele_old))   ! here, delta_rho is already incorporated into mass_ele and mass_ele_old
        else
          call addto(cmc, nodes, nodes, &
            ((coef-coef_old)/delta_rho)*mass_ele)
        end if
      end if
      if (move_mesh) then
        ! detwei and normals at the begin of the time step
        call transform_facet_to_physical(old_positions, sele, detwei_f=detwei,&
           & normal=normals)
        ! at each gauss point multiply with inner product of gravity and surface normal
        detwei=detwei*(-1.0)*sum(face_val_at_quad(gravity_normal,sele)*normals, dim=1)
        if (have_wd .and. .not. have_wd_node_int) then
             mass_ele_old=shape_shape(face_shape(prevp, sele), face_shape(prevp, sele), detwei*(1.0-alpha_wetdry_quad_prevp))
             mass_ele_old_wd=shape_shape(face_shape(prevp, sele), face_shape(prevp, sele), detwei*alpha_wetdry_quad_prevp)
        else
             mass_ele_old=shape_shape(face_shape(p, sele), face_shape(p, sele), detwei)
             if (have_wd .and. have_wd_node_int) then
                mass_ele_old_wd=mass_ele_old
                do i=1,size(mass_ele_old,1)
                  mass_ele_old(i,:)=mass_ele_old(i,:)*(1.0-alpha_wetdry_prevp)
                  mass_ele_old_wd(i,:)=mass_ele_old_wd(i,:)*alpha_wetdry_prevp
                end do
            end if
        end if

        if(present(rhs)) then
           call addto(rhs, nodes, &
                -(matmul(mass_ele, top_pressures) &
                -matmul(mass_ele_old, old_top_pressures))*alpha)
        end if
        if (have_wd .and. present(rhs)) then
           call addto(rhs, nodes, &
                +(matmul(mass_ele_wd, face_val(original_bottomdist_remap, sele)-d0) &
                -matmul(mass_ele_old_wd, face_val(original_bottomdist_remap,sele)-d0))*alpha*g)
        end if

      else
        ! no mesh movement - just use the same mass matrix as above
         if(present(rhs)) then
            if (variable_density) then
              call addto(rhs, nodes, &
                   -1.0*matmul(mass_ele, top_pressures-old_top_pressures)*alpha)  ! here, delta_rho is in mass_ele
            else
              call addto(rhs, nodes, &
                   -1.0*matmul(mass_ele, top_pressures-old_top_pressures)*alpha/delta_rho)
            end if
         end if
      end if
      
    end subroutine add_free_surface_element

    ! Computes alpha_wetdry. The resulting array is 0 if the node point is wet (p > -g d_0) and 1 if the node point is dry (p <= -g d_0)
    subroutine compute_alpha_wetdry(p, sele, alpha_wetdry)
      type(scalar_field), pointer, intent(in) :: p
      integer, intent(in) :: sele
      real, dimension(:), intent(inout) :: alpha_wetdry
      integer :: i

      alpha_wetdry = -face_val(p, sele)-face_val(original_bottomdist_remap, sele)*g+d0*g
       do i=1, size(alpha_wetdry)
           if (alpha_wetdry(i)>0.0)  then
               alpha_wetdry(i)=1.0
           else
               alpha_wetdry(i)=0.0
           end if
      end do
    end subroutine compute_alpha_wetdry
 
    ! Computes alpha_wetdry at each quadrature point. The resulting array is 0 if the quad point is wet (p > -g d_0) and 1 if the quad point is dry (p <= -g d_0)
    subroutine compute_alpha_wetdry_quad(p, sele, alpha_wetdry_quad)
      type(scalar_field), pointer, intent(in) :: p
      integer, intent(in) :: sele
      real, dimension(:), intent(inout) :: alpha_wetdry_quad
      integer :: i

      alpha_wetdry_quad = -face_val_at_quad(p, sele)-face_val_at_quad(original_bottomdist_remap, sele)*g+d0*g
       do i=1, size(alpha_wetdry_quad)
           if (alpha_wetdry_quad(i)>0.0)  then
               alpha_wetdry_quad(i)=1.0
           else
               alpha_wetdry_quad(i)=0.0
           end if
      end do
    end subroutine compute_alpha_wetdry_quad
 
  end subroutine add_free_surface_to_cmc_projection
    
  subroutine update_prognostic_free_surface(state, fs, implicit_prognostic_fs, explicit_prognostic_fs)
    !!< For the viscous free surface method, update the prognostic surface field
    !!< from the scaled free surface (\Delta\rho g\eta) that has just been solved with
    !!< the pressure projection (implicit) or time-integrated explicitly.
    !!< This is done via a small, surface Galerkin projection equation that is
    !!< assembled and solved for. If /geometry/ocean_boundaries are specified the updated
    !!< values are also extrapolated from the top surface.
    type(state_type), intent(inout):: state
    type(scalar_field), intent(inout):: fs
    logical, intent(in):: implicit_prognostic_fs, explicit_prognostic_fs

    type(vector_field), pointer:: u, x, gravity_normal
    type(scalar_field), pointer:: topdis
    character(len=FIELD_NAME_LEN):: bctype
    character(len=OPTION_PATH_LEN):: fs_option_path
    real:: g, rho0, dt
    integer:: i, j, stat

    type(integer_hash_table):: sele_to_fs_ele, sele_to_implicit_fs_ele
    integer, dimension(:), pointer:: surface_node_list, surface_element_list
    type(mesh_type), pointer:: surface_mesh
    type(csr_sparsity), pointer :: surface_sparsity
    type(csr_matrix):: fs_matrix
    type(scalar_field):: fs_rhs, surface_fs
    type(scalar_field), pointer :: density, external_density
    type(scalar_field), pointer :: scaled_fs, old_fs
    integer external_density_stat
    logical :: have_density, variable_density, move_mesh, have_external_density

    ewrite(1,*) 'Entering update_prognostic_free_surface'

    assert(implicit_prognostic_fs.or.explicit_prognostic_fs)

    u => extract_vector_field(state, "Velocity")
    assert(have_option(trim(u%option_path)//"/prognostic"))

    if(implicit_prognostic_fs) then
      if(.not. has_boundary_condition_name(fs, "_implicit_free_surface")) then
        call initialise_implicit_prognostic_free_surface(state, fs, u)
      end if

      call get_boundary_condition(fs, "_implicit_free_surface", &
                                  surface_element_list=surface_element_list)
      call invert_set(surface_element_list, sele_to_implicit_fs_ele)
      ! obtain the scaled f.s. that has been stored under the
      ! "_implicit_free_surface" boundary condition
      scaled_fs => extract_surface_field(fs, "_implicit_free_surface", "ScaledFreeSurface")

    end if

    if(explicit_prognostic_fs) then

      old_fs => extract_scalar_field(state, "OldFreeSurface")

    end if

    if (.not. has_boundary_condition_name(fs, "_free_surface")) then
      call initialise_prognostic_free_surface(fs, u)
    end if

    call get_boundary_condition(fs, "_free_surface", surface_mesh=surface_mesh, &
                                surface_element_list=surface_element_list, &
                                surface_node_list=surface_node_list)
    call invert_set(surface_element_list, sele_to_fs_ele)

    move_mesh = have_option("/mesh_adaptivity/mesh_movement/free_surface")
    x => extract_vector_field(state, "Coordinate")
    ! gravity acceleration
    call get_option('/physical_parameters/gravity/magnitude', g, stat=stat)
    if (stat/=0) then
      FLExit("For a free surface you need gravity")
    end if
    gravity_normal => extract_vector_field(state, "GravityDirection")
    density => extract_scalar_field(state, "Density", stat=stat)
    have_density = (stat==0)
    call get_fs_reference_density_from_options(rho0, state%option_path)
    call get_option('/timestepping/timestep', dt)

    surface_sparsity => get_csr_sparsity_firstorder(state, surface_mesh, surface_mesh)
    call allocate(fs_matrix, surface_sparsity, name="FSMatrix")
    call zero(fs_matrix)
    call allocate(fs_rhs, surface_mesh, name="FSRHS")
    call zero(fs_rhs)
    call allocate(surface_fs, surface_mesh, name="SurfaceFS")
    call zero(surface_fs)

    do i=1, get_boundary_condition_count(u)
      call get_boundary_condition(u, i, type=bctype, &
           surface_element_list=surface_element_list, &
           option_path=fs_option_path)
      if (bctype=="free_surface") then
        if(have_option(trim(fs_option_path)//"/type[0]/no_normal_stress")) then

          external_density => extract_scalar_surface_field(u, i, "ExternalDensity", stat=external_density_stat)
          have_external_density = external_density_stat==0

          variable_density = have_option(trim(fs_option_path)//"/type[0]/variable_density") &
                            .and. (.not.move_mesh)
          if (variable_density.and.(.not.have_density)) then
            FLExit("Variable density free surface requires a Density field.")
          end if

          if(have_option(trim(fs_option_path)//"/type[0]/no_normal_stress/explicit")) then

            do j=1, size(surface_element_list)
              call add_explicit_boundary_integral_sele(surface_element_list(j))
            end do

          else

            do j=1, size(surface_element_list)
              call add_implicit_boundary_integral_sele(j, surface_element_list(j))
            end do
          end if

        end if
      end if
    end do

    if(implicit_prognostic_fs) call deallocate(sele_to_implicit_fs_ele)
    call deallocate(sele_to_fs_ele)

    call petsc_solve(surface_fs, fs_matrix, fs_rhs, option_path=trim(fs%option_path))
    ewrite_minmax(surface_fs)

    call deallocate(fs_matrix)
    call deallocate(fs_rhs)

    call set(fs, surface_node_list, surface_fs%val)
    ewrite_minmax(fs)
    call deallocate(surface_fs)

    ! if /geometry/ocean_boundaries are specified take the new fs values
    ! at the top and extrapolate them downwards
    topdis => extract_scalar_field(state, "DistanceToTop", stat=stat)
    if (stat==0) then
       ! note we're not using the actual free_surface bc here, as 
       ! that may be specified in parts, or not cover the whole area
       call get_boundary_condition(topdis, 1, &
         surface_element_list=surface_element_list)
         
       x => extract_vector_field(state, "Coordinate")
 
       ! vertically extrapolate pressure values at the free surface downwards
       ! (reuse projected horizontal top surface mesh cached under DistanceToTop)
       call VerticalExtrapolation(fs, fs, x, &
         gravity_normal, surface_element_list=surface_element_list, &
         surface_name="DistanceToTop")
    end if

  contains

    subroutine add_implicit_boundary_integral_sele(surface_mesh_ele, sele)
      integer, intent(in):: surface_mesh_ele, sele ! element number in surface mesh, and facet number in u%mesh

      real, dimension(face_ngi(fs, sele)) :: detwei_bdy
      real, dimension(face_ngi(fs, sele)):: inv_delta_rho_g_quad

      if(variable_density .and. have_external_density) then
        inv_delta_rho_g_quad = 1.0/g/(face_val_at_quad(density, sele)-ele_val_at_quad(external_density, surface_mesh_ele))
      elseif (variable_density) then
        inv_delta_rho_g_quad = 1.0/g/face_val_at_quad(density, sele)
      elseif (have_external_density) then
        inv_delta_rho_g_quad = 1.0/g/(rho0-node_val(external_density, 1))
      else 
        inv_delta_rho_g_quad = 1.0/g/rho0
      end if
      
      
      call transform_facet_to_physical(x, sele, &
           detwei_f=detwei_bdy)

      call addto(fs_matrix, ele_nodes(surface_mesh, fetch(sele_to_fs_ele, sele)), &
                            ele_nodes(surface_mesh, fetch(sele_to_fs_ele, sele)), &
                            shape_shape(face_shape(fs, sele), face_shape(fs, sele), detwei_bdy))

      call addto(fs_rhs, ele_nodes(surface_mesh, fetch(sele_to_fs_ele, sele)),  &
                                               shape_rhs(face_shape(fs, sele),  &
                                                         detwei_bdy* &
                                                         ele_val_at_quad(scaled_fs, fetch(sele_to_implicit_fs_ele, sele))* &
                                                         inv_delta_rho_g_quad))

    end subroutine add_implicit_boundary_integral_sele

    subroutine add_explicit_boundary_integral_sele(sele)
      integer, intent(in):: sele

      real, dimension(face_ngi(u, sele)) :: detwei_bdy
      real, dimension(x%dim, face_ngi(u, sele)) :: normal_bdy

      call transform_facet_to_physical(x, sele, &
           detwei_f=detwei_bdy, normal=normal_bdy)

      call addto(fs_matrix, ele_nodes(surface_mesh, fetch(sele_to_fs_ele, sele)), &
                            ele_nodes(surface_mesh, fetch(sele_to_fs_ele, sele)), &
                            shape_shape(face_shape(fs, sele), face_shape(fs, sele), detwei_bdy))

      call addto(fs_rhs, ele_nodes(surface_mesh, fetch(sele_to_fs_ele, sele)),  &
                                               shape_rhs(face_shape(fs, sele),  &
                                                         detwei_bdy* &
                                                         (face_val_at_quad(old_fs, sele) - &
                                                         dt*sum(face_val_at_quad(u,sele)*normal_bdy, dim=1)/ &
                                                         sum(face_val_at_quad(gravity_normal,sele)*normal_bdy, dim=1))))

    end subroutine add_explicit_boundary_integral_sele

  end subroutine update_prognostic_free_surface
    
  subroutine add_free_surface_to_poisson_rhs(poisson_rhs, state, dt, theta_pg)
    !!< Add the rhs contributions of the fs terms in the continuity equation
    !!< to the initial Poisson equation
    type(scalar_field), intent(inout) :: poisson_rhs
    type(state_type), intent(inout) :: state
    real, intent(in) :: dt, theta_pg

    type(integer_hash_table):: sele_to_fs_ele
    type(vector_field), pointer:: positions, u, gravity_normal
    type(scalar_field), pointer:: p, free_surface, scaled_fs, density, external_density
    type(mesh_type), pointer:: embedded_fs_mesh
    character(len=FIELD_NAME_LEN):: bctype
    character(len=OPTION_PATH_LEN):: fs_option_path
    real g, coef, rho0
    integer, dimension(:), pointer:: surface_element_list, fs_surface_element_list
    integer i, j, grav_stat, dens_stat, external_density_stat
    logical:: include_normals
    logical:: implicit_prognostic_fs, use_fs_mesh
    logical:: move_mesh, variable_density, have_density, have_external_density

    ewrite(1,*) 'Entering assemble_masslumped_poisson_rhs_free_surface'

    ! gravity acceleration
    call get_option('/physical_parameters/gravity/magnitude', g, stat=grav_stat)
    if (grav_stat/=0) then
      FLExit("For a free surface you need gravity")
    end if

    ! with a free surface the initial condition prescribed for pressure
    ! is used at the free surface nodes only
    p => extract_scalar_field(state, "Pressure")
    u => extract_vector_field(state, "Velocity")
    
    ! the prognostic free surface can only be used with the no_normal_stress option
    implicit_prognostic_fs=has_implicit_viscous_free_surface_bc(u)
    if(implicit_prognostic_fs) then
      free_surface => extract_scalar_field(state, "FreeSurface")
      assert(have_option(trim(free_surface%option_path)//"/prognostic"))
      if (.not. has_boundary_condition_name(free_surface, "_implicit_free_surface")) then
        call initialise_implicit_prognostic_free_surface(state, free_surface, u)
      end if
      ! obtain the f.s. surface mesh that has been stored under the
      ! "_implicit_free_surface" boundary condition
      call get_boundary_condition(free_surface, "_implicit_free_surface", &
        surface_element_list=fs_surface_element_list)
      ! create a map from face numbers to element numbers in the fs mesh
      call invert_set(fs_surface_element_list, sele_to_fs_ele)
      scaled_fs => extract_surface_field(free_surface, "_implicit_free_surface", "ScaledFreeSurface")
      embedded_fs_mesh => extract_mesh(state, "_embedded_free_surface_mesh")
    end if
    
    ! reference density
    call get_fs_reference_density_from_options(rho0, state%option_path)
    density => extract_scalar_field(state, "Density", stat=dens_stat)
    have_density = (dens_stat==0)

    ! only include the inner product of gravity and surface normal
    ! if the free surface nodes are actually moved (not necessary
    ! for large scale ocean simulations)
    include_normals = have_option("/mesh_adaptivity/mesh_movement/free_surface")
    if (include_normals) then
      ewrite(2,*) 'Including inner product of normals in kinematic bc'
      gravity_normal => extract_vector_field(state, "GravityDirection")
    end if
      
    move_mesh = have_option("/mesh_adaptivity/mesh_movement/free_surface")
    ! adding in the free surface integral using the free surface
    ! elevation (p/rho0/g) specified by the inital pressure at the surface nodes
    ! or, with prognostic fs, use the free surface node values directly
    positions => extract_vector_field(state, "Coordinate")
      
    do i=1, get_boundary_condition_count(u)
      call get_boundary_condition(u, i, type=bctype, &
          surface_element_list=surface_element_list, &
          option_path=fs_option_path)
      if (bctype=="free_surface") then
        if (have_option(trim(fs_option_path)//"/type[0]/no_normal_stress") .and. &
            (.not.have_option(trim(fs_option_path)//"/type[0]/no_normal_stress/explicit"))) then
          assert(implicit_prognostic_fs)
          use_fs_mesh=.true.
        else
          use_fs_mesh=.false.
        end if
        coef=g*theta_pg**2*dt**2

        variable_density = have_option(trim(fs_option_path)//"/type[0]/variable_density") &
                          .and. (.not.move_mesh)
        if(variable_density.and.(.not.have_density)) then
          FLExit("Variable density free surface requires a Density field.")
        end if

        external_density => extract_scalar_surface_field(u, i, "ExternalDensity", stat=external_density_stat)
        have_external_density = external_density_stat==0

        do j=1, size(surface_element_list)
          call add_free_surface_element(j, surface_element_list(j))
        end do
      end if
    end do
      
    if (implicit_prognostic_fs) then
      call deallocate(sele_to_fs_ele)
    end if
      
    contains
    
    subroutine add_free_surface_element(surface_mesh_ele, sele)
    integer, intent(in):: surface_mesh_ele, sele
      
      real, dimension(face_loc(p, sele)):: values
      real, dimension(positions%dim, face_ngi(positions, sele)):: normals
      real, dimension(face_loc(p, sele), face_loc(p, sele)):: mass_ele
      real, dimension(face_ngi(p, sele)):: detwei
      real, dimension(face_ngi(p, sele)):: inv_delta_rho_quad
      integer, dimension(face_loc(p, sele)):: nodes
      integer:: ele
      
      ele=face_ele(positions, sele)
      call transform_facet_to_physical(positions, sele, detwei_f=detwei,&
           & normal=normals)
      if (include_normals) then
         ! at each gauss point multiply with inner product of gravity and surface normal
         detwei=detwei*(-1.0)*sum(face_val_at_quad(gravity_normal,sele)*normals, dim=1)
      end if

      if(variable_density .and. have_external_density) then
        inv_delta_rho_quad = 1.0/(face_val_at_quad(density, sele) - ele_val_at_quad(external_density, surface_mesh_ele))
        mass_ele=shape_shape(face_shape(p, sele), face_shape(p, sele), detwei*inv_delta_rho_quad)
      else if (variable_density) then
        mass_ele=shape_shape(face_shape(p, sele), face_shape(p, sele), detwei/face_val_at_quad(density, sele))
      else if (have_external_density) then 
        mass_ele = shape_shape(face_shape(p, sele), face_shape(p, sele), detwei) / (rho0 - node_val(external_density, 1))
      else
        mass_ele = shape_shape(face_shape(p, sele), face_shape(p, sele), detwei) / rho0
      end if
      
      if (use_fs_mesh) then
        nodes = ele_nodes(embedded_fs_mesh, fetch(sele_to_fs_ele, sele))
        values = ele_val(scaled_fs, fetch(sele_to_fs_ele, sele))
      else
        nodes = face_global_nodes(p, sele)
        values = face_val(p, sele)
      end if
      
      call addto(poisson_rhs, nodes, matmul(mass_ele, values)/coef)
      
    end subroutine add_free_surface_element
    
  end subroutine add_free_surface_to_poisson_rhs
    
  subroutine copy_poisson_solution_to_interior(state, p_theta, p, old_p, u)
    !!< Copy the solved for initial Poisson solution p_theta into p and old_p
    !!< but maintain initial condition for free surface nodes.
  type(state_type), intent(in):: state
  type(scalar_field), intent(inout), target:: p_theta, p, old_p
  type(vector_field), intent(in):: u
    
    type(scalar_field), pointer:: free_surface
    character(len=FIELD_NAME_LEN):: bctype
    character(len=OPTION_PATH_LEN):: fs_option_path
    logical:: prognostic_fs
    integer, dimension(:), pointer:: surface_element_list 
    integer:: fs_stat
    integer:: i, j, sele 
    
    ! the prognostic free surface can only be used with the no_normal_stress option
    free_surface => extract_scalar_field(state, "FreeSurface", stat=fs_stat)
    prognostic_fs=.false.
    if (fs_stat==0) then
      prognostic_fs=have_option(trim(free_surface%option_path)//"/prognostic")
    end if

    if (prognostic_fs) then
      ! only copy over solved for pressure values
      ! we keep the initial free surface
      call set_all(p, p_theta%val(1:node_count(p)))
    else
    
      ! first copy initial free surface elevations (p/g) at free surface nodes
      ! to p_theta
      do i=1, get_boundary_condition_count(u)
        call get_boundary_condition(u, i, type=bctype, &
            surface_element_list=surface_element_list, &
            option_path=fs_option_path)
        if (bctype=="free_surface") then

          if (have_option(trim(fs_option_path)//"/type[0]/no_normal_stress")) then
            ! this should have been options checked
            FLAbort("No normal stress free surface without prognostic free surface")
          end if

          do j=1, size(surface_element_list)
            sele=surface_element_list(j)
            call set(p_theta, face_global_nodes(p_theta,sele), face_val(p, sele))
          end do

        end if
      end do

      ! then copy everything (including interior) back from p_theta to p
      call set(p, p_theta)
    end if

    ! p and old_p should be the same (as we're in the first non-linear iteration)
    ! but they might be different fields (if #nonlinear iterations>1)
    call set(old_p, p)
    ewrite_minmax(p)

  end subroutine copy_poisson_solution_to_interior

  subroutine initialise_implicit_prognostic_free_surface(state, fs, u)
    !!< Setup ScaledFreeSurface and OldScaledFreeSurface surface fields (to contain
    !!< \Delta\rho g\eta values). Initialise them with initial condition (requires
    !!< little mass matrix solve).
    type(state_type), intent(inout):: state
    type(scalar_field), intent(inout):: fs
    type(vector_field), intent(in):: u

    character(len=OPTION_PATH_LEN):: fs_option_path
    character(len=FIELD_NAME_LEN):: bctype
    type(integer_set):: surface_elements
    type(mesh_type), pointer:: surface_mesh
    type(scalar_field):: scaled_fs, old_scaled_fs
    type(integer_hash_table):: sele_to_fs_ele
    integer, dimension(:), pointer:: surface_element_list, surface_node_list
    integer:: i


    ewrite(1,*) 'Entering initialise_implicit_prognostic_free_surface'

    call allocate(surface_elements)
    do i=1, get_boundary_condition_count(u)
      call get_boundary_condition(u, i, type=bctype, &
          surface_element_list=surface_element_list, &
          option_path=fs_option_path)
      if (bctype=="free_surface") then
         if(have_option(trim(fs_option_path)//"/type[0]/no_normal_stress") .and. &
            (.not.have_option(trim(fs_option_path)//"/type[0]/no_normal_stress/explicit"))) then
           ! include only implicit free surfaces
           call insert(surface_elements, surface_element_list)
         end if
      end if
    end do

   allocate(surface_element_list(1:key_count(surface_elements)))
   surface_element_list=set2vector(surface_elements)
   call add_boundary_condition_surface_elements(fs, "_implicit_free_surface", "free_surface", &
      surface_element_list)
   call invert_set(surface_element_list, sele_to_fs_ele)
   deallocate(surface_element_list)
   call deallocate(surface_elements)

   call get_boundary_condition(fs, "_implicit_free_surface", surface_mesh=surface_mesh, &
     surface_node_list=surface_node_list)
   if (IsParallel()) then
     call generate_surface_mesh_halos(fs%mesh, surface_mesh, surface_node_list)
   end if

   call allocate(scaled_fs, surface_mesh, "ScaledFreeSurface")
   call zero(scaled_fs)

   call solve_initial_scaled_free_surface(state, scaled_fs, u, sele_to_fs_ele, fs)

   call insert_surface_field(fs, "_implicit_free_surface", scaled_fs)

   call allocate(old_scaled_fs, surface_mesh, "OldScaledFreeSurface")
   call set(old_scaled_fs, scaled_fs)
   call deallocate(scaled_fs)

   call insert_surface_field(fs, "_implicit_free_surface", old_scaled_fs)
   call deallocate(old_scaled_fs)

   call deallocate(sele_to_fs_ele)


  end subroutine initialise_implicit_prognostic_free_surface

  subroutine solve_initial_scaled_free_surface(state, scaled_fs, u, sele_to_fs_ele, fs)
    !!< Solve for the initial value of the "ScaledFreeSurface" surface field, this is
    !!< obtained by solving scaled_fs = g*delta_rho*fs as a mass matrix equation on the surface
    !!< where fs is the initial condition (or interpolated after adaptivity) value of the prognostic fs field
    type(state_type), intent(inout):: state
    type(scalar_field), intent(inout):: scaled_fs ! the surface field solved for
    type(vector_field), intent(in):: u
    ! a map between surface element nos of the fs%mesh, and elements in scaled_fs%mesh:
    type(integer_hash_table), intent(in):: sele_to_fs_ele 
    type(scalar_field), intent(in):: fs ! the full, prognostic fs 

    character(len=FIELD_NAME_LEN):: bctype
    character(len=OPTION_PATH_LEN):: fs_option_path
    real:: rho0, g
    integer, dimension(:), pointer:: surface_element_list
    integer:: i, j, stat, external_density_stat

    type(csr_sparsity), pointer :: surface_sparsity
    type(csr_matrix):: fs_matrix
    type(scalar_field):: fs_rhs
    type(vector_field), pointer :: x
    type(scalar_field), pointer :: density, external_density
    logical :: have_density, variable_density, move_mesh, have_external_density

    ewrite(1,*) "Inside solve_initial_scaled_free_surface"

    move_mesh = have_option("/mesh_adaptivity/mesh_movement/free_surface")
    x => extract_vector_field(state, "Coordinate")
    call get_option('/physical_parameters/gravity/magnitude', g, stat=stat)
    if (stat/=0) then
      FLExit("For a free surface you need gravity")
    end if
    density => extract_scalar_field(state, "Density", stat=stat)
    have_density = (stat==0)
    call get_fs_reference_density_from_options(rho0, state%option_path)

    surface_sparsity => get_csr_sparsity_firstorder(state, scaled_fs%mesh, scaled_fs%mesh)
    call allocate(fs_matrix, surface_sparsity, name="FSMatrix")
    call zero(fs_matrix)
    call allocate(fs_rhs, scaled_fs%mesh, name="FSRHS")
    call zero(fs_rhs)

    do i=1, get_boundary_condition_count(u)
      call get_boundary_condition(u, i, type=bctype, &
        surface_element_list=surface_element_list, &
        option_path=fs_option_path)
      if (bctype=="free_surface") then
        if(have_option(trim(fs_option_path)//"/type[0]/no_normal_stress") .and. &
          (.not.have_option(trim(fs_option_path)//"/type[0]/no_normal_stress/explicit"))) then

          variable_density = have_option(trim(fs_option_path)//"/type[0]/variable_density") &
            .and. (.not.move_mesh)
          if (variable_density.and.(.not.have_density)) then
            FLExit("Variable density free surface requires a Density field.")
          end if

          external_density => extract_scalar_surface_field(u, i, "ExternalDensity", stat=external_density_stat)
          have_external_density = external_density_stat==0

          do j=1, size(surface_element_list)
            call add_boundary_integral_sele(j, surface_element_list(j))
          end do
        end if
      end if
    end do

    call petsc_solve(scaled_fs, fs_matrix, fs_rhs, option_path=trim(fs%option_path))
    ewrite_minmax(fs)
    ewrite_minmax(scaled_fs)

    call deallocate(fs_matrix)
    call deallocate(fs_rhs)

  contains

    subroutine add_boundary_integral_sele(surface_mesh_ele, sele)
      integer, intent(in):: surface_mesh_ele, sele

      real, dimension(face_ngi(fs, sele)) :: detwei_bdy
      real, dimension(face_ngi(fs, sele)):: delta_rho_quad
      integer :: fs_ele

      if(variable_density .and. have_external_density) then
        delta_rho_quad = face_val_at_quad(density, sele)-ele_val_at_quad(external_density, surface_mesh_ele)
      else if (variable_density) then
        delta_rho_quad = face_val_at_quad(density, sele)
      else if (have_external_density) then
        delta_rho_quad = rho0-node_val(external_density, 1)
      else
        delta_rho_quad = rho0
      end if

      call transform_facet_to_physical(x, sele, &
          detwei_f=detwei_bdy)

      fs_ele = fetch(sele_to_fs_ele, sele)
      call addto(fs_matrix, ele_nodes(scaled_fs, fs_ele), &
          ele_nodes(scaled_fs, fs_ele), &
          shape_shape(face_shape(fs, sele), face_shape(fs, sele), detwei_bdy))

      call addto(fs_rhs, ele_nodes(scaled_fs, fs_ele), &
          shape_rhs(face_shape(fs,sele),  &
              detwei_bdy * face_val_at_quad(fs,sele) * g * delta_rho_quad))

    end subroutine add_boundary_integral_sele

  end subroutine solve_initial_scaled_free_surface

  subroutine initialise_prognostic_free_surface(fs, u)
    type(scalar_field), intent(inout):: fs
    type(vector_field), intent(in):: u

    type(mesh_type), pointer:: surface_mesh
    type(integer_set):: surface_elements
    character(len=FIELD_NAME_LEN):: bctype
    character(len=OPTION_PATH_LEN):: fs_option_path
    integer, dimension(:), pointer:: surface_element_list, surface_node_list
    integer:: i

    ewrite(1,*) 'Entering initialise_prognostic_free_surface'

    call allocate(surface_elements)
    do i=1, get_boundary_condition_count(u)
      call get_boundary_condition(u, i, type=bctype, &
          surface_element_list=surface_element_list, &
          option_path=fs_option_path)
      if (bctype=="free_surface") then
         if(have_option(trim(fs_option_path)//"/type[0]/no_normal_stress")) then
           ! include all explicit and implicit free surfaces
           call insert(surface_elements, surface_element_list)
         end if
      end if
    end do

   allocate(surface_element_list(1:key_count(surface_elements)))
   surface_element_list=set2vector(surface_elements)
   call add_boundary_condition_surface_elements(fs, "_free_surface", "free_surface", &
      surface_element_list)
   deallocate(surface_element_list)
   call deallocate(surface_elements)
   if (IsParallel()) then
     call get_boundary_condition(fs, "_free_surface", surface_mesh=surface_mesh, &
       surface_node_list=surface_node_list)
     call generate_surface_mesh_halos(fs%mesh, surface_mesh, surface_node_list)
   end if

  end subroutine initialise_prognostic_free_surface

  function get_extended_pressure_mesh_for_viscous_free_surface(state, pressure_mesh, fs) result (extended_pressure_mesh)
    ! extend the pressure mesh to contain the extra free surface dofs
    ! (doubling the pressure nodes at the free surface). The returned mesh uses
    ! a new node numbering that includes the pressure and separate fs nodes, but has
    ! the same elements as the pressure_mesh. This routine also creates, seperately, a surface
    ! mesh (the "embedded" fs mesh) that uses the same combined pressure and fs node numbering, but only contains the surface
    ! elements with fs nodes as its elements. Both are stored in state.

    ! this must be the state containing the prognostic pressure and free surface
    type(state_type), intent(inout):: state
    ! the origional pressure mesh (without fs nodes)
    type(mesh_type), intent(in):: pressure_mesh
    ! the prognostic free surface field
    type(scalar_field), intent(inout):: fs
    type(mesh_type), pointer:: extended_pressure_mesh

    type(integer_vector), dimension(2):: new_node_map
    type(vector_field), pointer:: u
    type(mesh_type), pointer:: fs_mesh
    type(mesh_type):: extended_mesh, embedded_fs_mesh

    integer, dimension(ele_loc(pressure_mesh, 1)):: new_nodes
    integer, dimension(face_loc(fs%mesh, 1)):: new_fs_nodes
    integer, dimension(:), pointer:: nodes
    integer :: ihalo, ele, j, pnodes
    logical :: parallel

    if (.not. has_mesh(state, "_extended_pressure_mesh")) then

      if (.not. has_boundary_condition_name(fs, "_implicit_free_surface")) then
        u => extract_vector_field(state, "Velocity")
        assert(have_option(trim(u%option_path)//"/prognostic"))
        call initialise_implicit_prognostic_free_surface(state, fs, u)
      end if
      ! obtain the f.s. surface mesh that has been stored under the
      ! "_implicit_free_surface" boundary condition (this is a surface
      ! mesh with a separate fs node numbering)
      call get_boundary_condition(fs, "_implicit_free_surface", &
          surface_mesh=fs_mesh)

      if (associated(pressure_mesh%halos) .and. associated(fs_mesh%halos)) then
        assert(size(pressure_mesh%halos)==2 .and. size(fs_mesh%halos)==2)

        ! create a node numbering for the combined pressure and fs system
        ! that is trailing receives consistent
        new_node_map = create_combined_numbering_trailing_receives( &
           (/ pressure_mesh%halos(1), fs_mesh%halos(1) /), &
           (/ pressure_mesh%halos(2), fs_mesh%halos(2) /))
        parallel = .true.
      else
        assert(.not. associated(pressure_mesh%halos))
        assert(.not. associated(fs_mesh%halos))
        parallel = .false.
      end if

      call allocate(extended_mesh, node_count(pressure_mesh)+node_count(fs_mesh), &
          element_count(pressure_mesh), pressure_mesh%shape, &
          "Extended"//trim(pressure_mesh%name))
      extended_mesh%periodic = pressure_mesh%periodic

      do ele=1, element_count(pressure_mesh)
        nodes => ele_nodes(pressure_mesh,ele)
        if (parallel) then
          do j=1, size(nodes)
            new_nodes(j) = new_node_map(1)%ptr(nodes(j))
          end do
          call set_ele_nodes(extended_mesh, ele, new_nodes)
        else 
          call set_ele_nodes(extended_mesh, ele, nodes)
        end if
      end do

      call add_faces(extended_mesh, model=pressure_mesh)

      if (parallel) then
        ! Allocate extended mesh halos:
        allocate(extended_mesh%halos(2))

        ! Derive extended_mesh nodal halos:
        do ihalo = 1, 2

           extended_mesh%halos(ihalo) = combine_halos( &
             (/ pressure_mesh%halos(ihalo), fs_mesh%halos(ihalo) /), &
             new_node_map, &
             name="Extended"//trim(pressure_mesh%halos(ihalo)%name))

           assert(trailing_receives_consistent(extended_mesh%halos(ihalo)))
           assert(halo_valid_for_communication(extended_mesh%halos(ihalo)))
           call create_global_to_universal_numbering(extended_mesh%halos(ihalo))
           call create_ownership(extended_mesh%halos(ihalo))
           
        end do
      end if

      call insert(state, extended_mesh, "_extended_pressure_mesh")
      call deallocate(extended_mesh)

      ! now create an auxillary surface mesh, that has the same topology as surface_mesh
      ! but uses the node numbering for free surface nodes of the extended mesh
      call allocate(embedded_fs_mesh, node_count(extended_mesh), &
        element_count(fs_mesh), fs_mesh%shape, &
        name="Emmbedded"//trim(fs_mesh%name))

      pnodes = node_count(pressure_mesh)
      do ele=1, element_count(fs_mesh)
        nodes => ele_nodes(fs_mesh,ele)
        if (parallel) then
          do j=1, size(nodes)
            new_fs_nodes(j) = new_node_map(2)%ptr(nodes(j))
          end do
        else
          new_fs_nodes = nodes + pnodes ! fs nodes come after the pressure nodes in the extended mesh
        end if
        call set_ele_nodes(embedded_fs_mesh, ele, new_fs_nodes)
      end do

      call insert(state, embedded_fs_mesh, "_embedded_free_surface_mesh")
      call deallocate(embedded_fs_mesh)

      if (parallel) then
        ! the node maps are no longer needed: the correspondence between
        ! pressure mesh/fs surface mesh and extended_mesh/embedded_fs_mesh resp.
        ! is maintained by looping over the elements of both simultaneously
        deallocate(new_node_map(1)%ptr)
        deallocate(new_node_map(2)%ptr)
      end if
      
    end if

    extended_pressure_mesh => extract_mesh(state, "_extended_pressure_mesh")

  end function get_extended_pressure_mesh_for_viscous_free_surface
  
  subroutine copy_to_extended_p(p, fs, theta_pg, p_theta)
    type(scalar_field), intent(in):: p, fs
    real, intent(in):: theta_pg
    type(scalar_field), intent(inout):: p_theta
    ! copy p and fs into p_theta that is allocated on the extended pressure mesh

    type(scalar_field), pointer:: scaled_fs, old_scaled_fs
    integer:: powned_nodes, fsowned_nodes

    ! obtain the scaled f.s. that has been stored under the
    ! "_implicit_free_surface" boundary condition
    scaled_fs => extract_surface_field(fs, "_implicit_free_surface", "ScaledFreeSurface")
    old_scaled_fs => extract_surface_field(fs, "_implicit_free_surface", "OldScaledFreeSurface")

    powned_nodes = nowned_nodes(p)
    fsowned_nodes = nowned_nodes(scaled_fs)
    
    ! p is not theta weighted (as usual for incompressible)
    p_theta%val(1:powned_nodes) = p%val(1:powned_nodes)

    ! setting this to the old scaled fs values means we're calculating the change in the free surface across the timestep
    ! not just since the last solve - unlike pressure where delta_p is since the last solve!
    p_theta%val(powned_nodes+1:powned_nodes+fsowned_nodes) = theta_pg*scaled_fs%val(1:fsowned_nodes) + (1.-theta_pg)*old_scaled_fs%val(1:fsowned_nodes)

    ! we only copy the owned nodes into to combined pressure/fs field p_theta
    ! because the renumbering of the receive nodes is a little more complicated
    ! therefore, do a halo update to give the receive nodes its correct values as well
    call halo_update(p_theta)

  end subroutine copy_to_extended_p

  subroutine update_pressure_and_viscous_free_surface(state, p, fs, delta_p, theta_pg)
    type(state_type), intent(inout):: state
    ! after the pressure+fs projection add in the solved for delta_p
    ! into the pressure and prognostic fs
    type(scalar_field), intent(inout):: p, fs
    type(scalar_field), intent(in):: delta_p
    real, intent(in):: theta_pg

    type(vector_field), pointer:: u
    type(scalar_field), pointer :: scaled_fs
    integer:: powned_nodes, fsowned_nodes
    real:: dt

    ewrite(1,*) 'Entering update_pressure_and_viscous_free_surface'

    u => extract_vector_field(state, "Velocity")
    if (.not. has_boundary_condition_name(fs, "_implicit_free_surface")) then
      assert(have_option(trim(u%option_path)//"/prognostic"))
      call initialise_implicit_prognostic_free_surface(state, fs, u)
    end if
    ! obtain the scaled f.s. that has been stored under the
    ! "_implicit_free_surface" boundary condition
    scaled_fs => extract_surface_field(fs, "_implicit_free_surface", "ScaledFreeSurface")

    call get_option('/timestepping/timestep', dt)

    powned_nodes = nowned_nodes(p)
    fsowned_nodes = nowned_nodes(scaled_fs)

    p%val(1:powned_nodes) = p%val(1:powned_nodes) + delta_p%val(1:powned_nodes)/dt
    call halo_update(p)

    scaled_fs%val(1:fsowned_nodes) = scaled_fs%val(1:fsowned_nodes) + delta_p%val(powned_nodes+1:powned_nodes+fsowned_nodes)/(dt*theta_pg)
    call halo_update(scaled_fs)
    ewrite_minmax(scaled_fs)

  end subroutine update_pressure_and_viscous_free_surface

  function get_extended_velocity_divergence_matrix(state, u, fs, &
      extended_mesh, get_ct, ct_m_name) result (ct_m_ptr)
    ! returns the velocity divergence matrix from state (or creates a new one if none present)
    ! with extra rows associated with prognostic fs nodes
    ! this routine mimicks the behaviour of get_velocity_divergence_matrix()
    type(state_type), intent(inout):: state
    type(vector_field), intent(in):: u
    type(scalar_field), intent(in):: fs
    type(mesh_type), intent(in):: extended_mesh
    ! returns .true. if the matrix needs to be reassembled (because it has just been allocated or because of mesh movement)
    logical, optional, intent(out):: get_ct
    character(len=*), intent(in), optional :: ct_m_name
    type(block_csr_matrix), pointer:: ct_m_ptr

    type(block_csr_matrix):: new_ct_m
    type(csr_sparsity), pointer:: sparsity
    integer:: stat
    integer, save:: last_mesh_movement = -1
    logical:: mesh_moved
    character(len=FIELD_NAME_LEN) :: l_ct_m_name

    mesh_moved = eventcount(EVENT_MESH_MOVEMENT)/=last_mesh_movement
    last_mesh_movement = eventcount(EVENT_MESH_MOVEMENT)

    ! Form the ct_m_name dependent on interface argument
    if (present(ct_m_name)) then
       l_ct_m_name = trim(ct_m_name)
    else
       l_ct_m_name = "ExtendedVelocityDivergenceMatrix"
    end if

    ct_m_ptr => extract_block_csr_matrix(state, l_ct_m_name, stat=stat)
    if (stat==0) then
      if (present(get_ct)) then
        get_ct = mesh_moved
      end if
      return
    end if

    sparsity => get_extended_velocity_divergence_sparsity(state, u, fs, extended_mesh)

    call allocate(new_ct_m, sparsity, blocks=(/ 1, u%dim /), name=l_ct_m_name)
    call insert(state, new_ct_m, new_ct_m%name)
    call deallocate(new_ct_m)

    ct_m_ptr => extract_block_csr_matrix(state, l_ct_m_name)

    if (present(get_ct)) then
      get_ct = .true.
    end if

  end function get_extended_velocity_divergence_matrix

  function get_extended_velocity_divergence_sparsity(state, u, fs, extended_mesh) result (sparsity)
    ! returns the sparsity of the velocity divergence matrix from state (or creates a new one if none present)
    ! with extra rows associated with prognostic fs nodes
    ! this routine mimicks the behaviour of get_velocity_divergence_matrix()
    type(state_type), intent(inout):: state
    type(vector_field), intent(in):: u
    type(scalar_field), intent(in):: fs
    type(mesh_type), intent(in):: extended_mesh
    type(csr_sparsity), pointer:: sparsity

    type(mesh_type), pointer:: embedded_fs_mesh
    type(mesh_type):: u_surface_mesh
    type(csr_sparsity):: new_sparsity, fs_sparsity, p_sparsity
    integer, dimension(:), pointer :: fs_surface_element_list
    integer, dimension(:), pointer :: p_row, fs_row, new_row
    integer:: entries, stat, i, sele

    sparsity => extract_csr_sparsity(state, "ExtendedVelocityDivergenceSparsity", stat=stat)
    if (stat==0) return

    ! this will create a sparsity with only the rows associated with the pressure 
    ! nodes filled in (in the extended mesh numbering), fs nodes have zero row length
    p_sparsity = make_sparsity(extended_mesh, u%mesh, name="TempPartialSparsityPressure")

    ! now do the same for rows associated with fs nodes in the combined mesh
    ! to make this sparsity we need a velocity surface mesh that uses the
    ! velocity dof numbering of the entire mesh
    call get_boundary_condition(fs, "_implicit_free_surface", &
        surface_element_list=fs_surface_element_list)
    call allocate(u_surface_mesh, node_count(u), size(fs_surface_element_list), &
      face_shape(u,1), &
      name="TempVelocityFreeSurfaceMesh")
    do i=1, size(fs_surface_element_list)
      sele = fs_surface_element_list(i)
      call set_ele_nodes(u_surface_mesh, i, face_global_nodes(u, sele))
    end do
    embedded_fs_mesh => extract_mesh(state, "_embedded_free_surface_mesh")
    fs_sparsity = make_sparsity(embedded_fs_mesh, u_surface_mesh, name="TempPartialSparsityFreeSurface")
    call deallocate(u_surface_mesh)

    ! now we simply combine the two sparsities
    entries = size(p_sparsity%colm) + size(fs_sparsity%colm)
    call allocate(new_sparsity, node_count(extended_mesh), node_count(u), &
      entries, name="ExtendedVelocityDivergenceSparsity")
    new_sparsity%findrm(1) = 1
    do i=1, node_count(extended_mesh)
      new_sparsity%findrm(i+1) = new_sparsity%findrm(i) + &
        row_length(p_sparsity, i) + row_length(fs_sparsity, i)
    end do
    do i=1, node_count(extended_mesh)
      p_row => row_m_ptr(p_sparsity, i)
      fs_row => row_m_ptr(fs_sparsity, i)
      new_row => row_m_ptr(new_sparsity, i)
      if (size(new_row)==size(p_row)) then
        new_row = p_row
      else if (size(new_row)==size(fs_row)) then
        new_row = fs_row
      else
        ! we assume that each row is either associated with a fs node, in which case size(p_row)==0 
        ! and size(new_row)=size(new_row), or associated with a p node, in which case size(fs_row)==0
        ! and size(new_row)=size(p_row)
        FLAbort("Node in combined mesh that is neither fs or pressure node (or both).")
      end if
    end do
    new_sparsity%sorted_rows = p_sparsity%sorted_rows .and. fs_sparsity%sorted_rows
    call deallocate(p_sparsity)
    call deallocate(fs_sparsity)

    if (associated(extended_mesh%halos)) then
      allocate(new_sparsity%row_halo)
      new_sparsity%row_halo = extended_mesh%halos(1)
      call incref(new_sparsity%row_halo)
      allocate(new_sparsity%column_halo)
      new_sparsity%column_halo = u%mesh%halos(1)
      call incref(new_sparsity%column_halo)
    end if
    call insert(state, new_sparsity, new_sparsity%name)
    call deallocate(new_sparsity)

    sparsity => extract_csr_sparsity(state, "ExtendedVelocityDivergenceSparsity")

  end function get_extended_velocity_divergence_sparsity

  function get_extended_pressure_poisson_matrix(state, ct_m, extended_mesh, get_cmc) result (cmc_m_ptr)
    ! returns the pressure poisson matrix from state (or creates a new one if none present)
    ! with extra rows and columns associated with prognostic fs nodes
    ! this routine mimicks the behaviour of get_pressure_poisson_matrix()
    type(state_type), intent(inout):: state
    type(block_csr_matrix), intent(in):: ct_m
    type(mesh_type), intent(in):: extended_mesh
    ! returns .true. if the matrix needs to be reassembled (because it has just been allocated or because of mesh movement)
    logical, optional, intent(out):: get_cmc
    type(csr_matrix), pointer:: cmc_m_ptr

    type(csr_matrix):: cmc_m
    type(csr_sparsity):: grad_sparsity, cmc_sparsity
    integer, save:: last_mesh_movement=-1
    integer:: stat
    logical:: mesh_moved

    mesh_moved = eventcount(EVENT_MESH_MOVEMENT)/=last_mesh_movement
    last_mesh_movement = eventcount(EVENT_MESH_MOVEMENT)

    cmc_m_ptr => extract_csr_matrix(state, "ExtendedPressurePoissonMatrix", stat=stat)
    if (stat==0) then
      if (present(get_cmc)) then
        get_cmc = mesh_moved
      end if
      return
    end if

    grad_sparsity=transpose(ct_m%sparsity)
    cmc_sparsity = matmul(ct_m%sparsity, grad_sparsity)

    if (associated(extended_mesh%halos)) then
      allocate(cmc_sparsity%row_halo)
      cmc_sparsity%row_halo = extended_mesh%halos(2)
      call incref(cmc_sparsity%row_halo)
      allocate(cmc_sparsity%column_halo)
      cmc_sparsity%column_halo = extended_mesh%halos(2)
      call incref(cmc_sparsity%column_halo)
    end if

    call allocate(cmc_m, cmc_sparsity, name="ExtendedPressurePoissonMatrix")
    call insert(state, cmc_m, cmc_m%name)
    call deallocate(cmc_m)
    call deallocate(cmc_sparsity)
    call deallocate(grad_sparsity)

    cmc_m_ptr => extract_csr_matrix(state, "ExtendedPressurePoissonMatrix")
    if (present(get_cmc)) then
      get_cmc = mesh_moved
    end if

  end function get_extended_pressure_poisson_matrix

  function get_extended_schur_auxillary_sparsity(state, ct_m, extended_mesh) result (aux_sparsity_ptr)
    ! returns the sparsity of the schur auxillary matrix from state (or creates a new one if none present)
    ! with extra rows and columns associated with prognostic fs nodes
    ! here it is assumed that the matrix has the same sparsity as the pressure poisson matrix has, obtained
    ! by sparsity multiplication of C^T * C
    type(state_type), intent(inout):: state
    type(block_csr_matrix), intent(in):: ct_m
    type(mesh_type), intent(in):: extended_mesh
    type(csr_sparsity), pointer:: aux_sparsity_ptr

    type(csr_sparsity):: grad_sparsity, aux_sparsity

    if (.not. has_csr_sparsity(state, "ExtendedSchurAuxillarySparsity")) then
      grad_sparsity=transpose(ct_m%sparsity)
      aux_sparsity = matmul(ct_m%sparsity, grad_sparsity)
      aux_sparsity%name="ExtendedSchurAuxillarySparsity"

      if (associated(extended_mesh%halos)) then
        allocate(aux_sparsity%row_halo)
        aux_sparsity%row_halo = extended_mesh%halos(2)
        call incref(aux_sparsity%row_halo)
        allocate(aux_sparsity%column_halo)
        aux_sparsity%column_halo = extended_mesh%halos(2)
        call incref(aux_sparsity%column_halo)
      end if

      call insert(state, aux_sparsity, aux_sparsity%name)
      call deallocate(aux_sparsity)
      call deallocate(grad_sparsity)
    end if

    aux_sparsity_ptr => extract_csr_sparsity(state, "ExtendedSchurAuxillarySparsity")

  end function get_extended_schur_auxillary_sparsity

  subroutine add_implicit_viscous_free_surface_integrals(state, ct_m, u, &
      p_mesh, fs)
    ! This routine adds in the boundary conditions for the viscous free surface
    ! (that is the free_surface bc with the no_normal_stress option)
    ! ct_m has been extended with some extra rows (corresponding to free surface
    ! nodes) that are used to enforce the kinematic bc, the transpose of that
    ! will produce the \rho_0 g\eta term in the no_normal_stress boundary condition:
    !   n\cdot\tau\cdot n + p - \rho_0 g\eta = 0
    ! if pressure is also extended to contain \rho g\eta in the extra nodes.
    
    type(state_type), intent(inout):: state
    type(block_csr_matrix), intent(inout):: ct_m
    type(vector_field), intent(in):: u
    type(mesh_type), intent(in):: p_mesh ! the extended pressure mesh
    type(scalar_field), intent(inout):: fs

    type(vector_field), pointer:: x
    type(mesh_type), pointer:: fs_mesh, embedded_fs_mesh
    type(integer_hash_table):: sele_to_fs_ele
    character(len=FIELD_NAME_LEN):: bc_type
    character(len=OPTION_PATH_LEN):: bc_option_path
    integer, dimension(:), pointer:: surface_element_list
    integer:: i, j

    assert(have_option(trim(fs%option_path)//"/prognostic"))

    if (.not. has_boundary_condition_name(fs, "_implicit_free_surface")) then
      call initialise_implicit_prognostic_free_surface(state, fs, u)
    end if
    ! obtain the f.s. surface mesh that has been stored under the
    ! "_implicit_free_surface" boundary condition
    call get_boundary_condition(fs, "_implicit_free_surface", &
        surface_mesh=fs_mesh, surface_element_list=surface_element_list)
    ! create a map from face numbers to element numbers in fs_mesh
    call invert_set(surface_element_list, sele_to_fs_ele)
    embedded_fs_mesh => extract_mesh(state, "_embedded_free_surface_mesh")

    x => extract_vector_field(state, "Coordinate")

    do i=1, get_boundary_condition_count(u)
      call get_boundary_condition(u, i, type=bc_type, option_path=bc_option_path, &
          surface_element_list=surface_element_list)
      if (bc_type=="free_surface") then
        if(have_option(trim(bc_option_path)//"/type[0]/no_normal_stress") .and. &
           (.not.have_option(trim(bc_option_path)//"/type[0]/no_normal_stress/explicit"))) then
          do j=1, size(surface_element_list)
            call add_boundary_integral_sele(surface_element_list(j))
          end do
        end if
      end if
    end do

  contains

    subroutine add_boundary_integral_sele(sele)
      integer, intent(in):: sele

      real, dimension(u%dim, face_loc(p_mesh, sele), face_loc(u, sele)) :: ct_mat_bdy
      real, dimension(u%dim, face_loc(fs, sele), face_loc(u, sele)) :: ht_mat_bdy
      real, dimension(face_ngi(u, sele)) :: detwei_bdy
      real, dimension(u%dim, face_ngi(u, sele)) :: normal_bdy
      integer:: dim
      
      call transform_facet_to_physical(x, sele, &
           detwei_f=detwei_bdy, normal=normal_bdy)
      ct_mat_bdy = shape_shape_vector(face_shape(p_mesh, sele), face_shape(u, sele), &
           detwei_bdy, normal_bdy)
      ht_mat_bdy = shape_shape_vector(face_shape(fs, sele), face_shape(u, sele), &
           detwei_bdy, normal_bdy)
      do dim=1, u%dim
        ! we've integrated continuity by parts, but not yet added in the resulting
        ! surface integral - for the non-viscous free surface this is left
        ! out to enforce the kinematic bc
        call addto(ct_m, 1, dim, face_global_nodes(p_mesh,sele), &
             face_global_nodes(u,sele), ct_mat_bdy(dim,:,:))
        ! for the viscous bc however we add this bc in the extra rows at the bottom of ct_m
        ! this integral will also enforce the \rho_0 g\eta term in the no_normal_stress bc:
        !   n\cdot\tau\cdot n + p - (\rho_0-\rho_external) g\eta = 0
        call addto(ct_m, 1, dim, ele_nodes(embedded_fs_mesh, fetch(sele_to_fs_ele, sele)), &
             face_global_nodes(u,sele), -ht_mat_bdy(dim,:,:))
      end do

    end subroutine add_boundary_integral_sele

  end subroutine add_implicit_viscous_free_surface_integrals

  subroutine add_implicit_viscous_free_surface_integrals_cv(state, ct_m, u, &
      p_mesh, fs)
    ! This routine adds in the boundary conditions for the viscous free surface
    ! (that is the free_surface bc with the no_normal_stress option)
    ! ct_m has been extended with some extra rows (corresponding to free surface
    ! nodes) that are used to enforce the kinematic bc, the transpose of that
    ! will produce the \rho_0 g\eta term in the no_normal_stress boundary condition:
    !   n\cdot\tau\cdot n + p - \rho_0 g\eta = 0
    ! if pressure is also extended to contain \rho g\eta in the extra nodes.
    
    type(state_type), intent(inout):: state
    type(block_csr_matrix), intent(inout):: ct_m
    type(vector_field), intent(in):: u
    type(mesh_type), intent(in):: p_mesh ! the extended pressure mesh
    type(scalar_field), intent(inout):: fs

    type(vector_field), pointer:: x
    type(mesh_type), pointer:: fs_mesh, embedded_fs_mesh
    type(integer_hash_table):: sele_to_fs_ele
    character(len=FIELD_NAME_LEN):: bc_type
    character(len=OPTION_PATH_LEN):: bc_option_path
    integer, dimension(:), pointer:: surface_element_list
    integer:: i, j

    ! information about cv faces
    type(cv_faces_type) :: cvfaces
    ! shape functions for region and surface
    type(element_type) :: x_cvbdyshape
    type(element_type) :: p_cvbdyshape
    type(element_type) :: u_cvbdyshape
    integer:: quaddegree, ele, sele

    assert(have_option(trim(fs%option_path)//"/prognostic"))

    if (.not. has_boundary_condition_name(fs, "_implicit_free_surface")) then
      call initialise_implicit_prognostic_free_surface(state, fs, u)
    end if
    ! obtain the f.s. surface mesh that has been stored under the
    ! "_implicit_free_surface" boundary condition
    call get_boundary_condition(fs, "_implicit_free_surface", &
        surface_mesh=fs_mesh, surface_element_list=surface_element_list)
    ! create a map from face numbers to element numbers in fs_mesh
    call invert_set(surface_element_list, sele_to_fs_ele)
    embedded_fs_mesh => extract_mesh(state, "_embedded_free_surface_mesh")

    x => extract_vector_field(state, "Coordinate")

    call get_option("/geometry/quadrature/controlvolume_surface_degree", &
                   quaddegree, default=1)

    cvfaces=find_cv_faces(vertices=ele_vertices(p_mesh, 1), &
                          dimension=mesh_dim(p_mesh), &
                          polydegree=p_mesh%shape%degree, &
                          quaddegree=quaddegree)

    x_cvbdyshape=make_cvbdy_element_shape(cvfaces, x%mesh%faces%shape)
    p_cvbdyshape=make_cvbdy_element_shape(cvfaces, p_mesh%faces%shape)
    u_cvbdyshape=make_cvbdy_element_shape(cvfaces, u%mesh%faces%shape)

    do i=1, get_boundary_condition_count(u)
      call get_boundary_condition(u, i, type=bc_type, option_path=bc_option_path, &
          surface_element_list=surface_element_list)
      if (bc_type=="free_surface") then
        if(have_option(trim(bc_option_path)//"/type[0]/no_normal_stress") .and. &
           (.not.have_option(trim(bc_option_path)//"/type[0]/no_normal_stress/explicit"))) then
          do j=1, size(surface_element_list)
            sele = surface_element_list(j)
            ele = face_ele(x, sele)
            call add_boundary_integral_sele(sele, ele)
          end do
        end if
      end if
    end do

    call deallocate(x_cvbdyshape)
    call deallocate(p_cvbdyshape)
    call deallocate(u_cvbdyshape)
    call deallocate(cvfaces)

  contains

    subroutine add_boundary_integral_sele(sele, ele)
      integer, intent(in):: sele, ele

      real, dimension(u%dim, face_loc(p_mesh, sele), face_loc(u, sele)) :: ct_mat_bdy
      real, dimension(u%dim, face_loc(fs, sele), face_loc(u, sele)) :: ht_mat_bdy
      real, dimension(face_ngi(u, sele)) :: detwei_bdy
      real, dimension(x_cvbdyshape%ngi) :: detwei_bdy_cv
      real, dimension(x%dim, face_ngi(x, sele)) :: normal_bdy
      real, dimension(x%dim, x_cvbdyshape%ngi) :: normal_bdy_cv
      real, dimension(x%dim, ele_loc(x, ele)) :: x_ele
      real, dimension(x%dim, face_loc(x, sele)) :: x_ele_bdy
      integer:: dim, gi, ggi, iloc, jloc, face
      
      x_ele = ele_val(x, ele)
      x_ele_bdy = face_val(x, sele)

      call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, &
                            x_cvbdyshape, normal_bdy_cv, detwei_bdy_cv)
  
      ct_mat_bdy = 0.0
      surface_nodal_loop_i: do iloc = 1, face_loc(p_mesh, sele)

        surface_face_loop: do face = 1, cvfaces%sfaces
          if(cvfaces%sneiloc(iloc,face)/=0) then

            surface_quadrature_loop: do gi = 1, cvfaces%shape%ngi

              ggi = (face-1)*cvfaces%shape%ngi + gi

              surface_nodal_loop_j: do jloc = 1, face_loc(u, sele)

                surface_inner_dimension_loop: do dim = 1, size(normal_bdy_cv,1)

                   ct_mat_bdy(dim, iloc, jloc) =  ct_mat_bdy(dim, iloc, jloc) + &
                         u_cvbdyshape%n(jloc,ggi)*detwei_bdy_cv(ggi)*normal_bdy_cv(dim, ggi)

                end do surface_inner_dimension_loop

              end do surface_nodal_loop_j

            end do surface_quadrature_loop

          end if

        end do surface_face_loop

      end do surface_nodal_loop_i

      call transform_facet_to_physical(x, sele, &
           detwei_f=detwei_bdy, normal=normal_bdy)
      ht_mat_bdy = shape_shape_vector(face_shape(fs, sele), face_shape(u, sele), &
           detwei_bdy, normal_bdy)


      do dim=1, u%dim
        ! we've integrated continuity by parts, but not yet added in the resulting
        ! surface integral - for the non-viscous free surface this is left
        ! out to enforce the kinematic bc
        call addto(ct_m, 1, dim, face_global_nodes(p_mesh,sele), &
             face_global_nodes(u,sele), ct_mat_bdy(dim,:,:))
        ! for the viscous bc however we add this bc in the extra rows at the bottom of ct_m
        ! this integral will also enforce the \rho_0 g\eta term in the no_normal_stress bc:
        !   n\cdot\tau\cdot n + p - (\rho_0-\rho_external) g\eta = 0
        call addto(ct_m, 1, dim, ele_nodes(embedded_fs_mesh, fetch(sele_to_fs_ele, sele)), &
             face_global_nodes(u,sele), -ht_mat_bdy(dim,:,:))
      end do

    end subroutine add_boundary_integral_sele

  end subroutine add_implicit_viscous_free_surface_integrals_cv

  subroutine add_explicit_viscous_free_surface_integrals(state, mom_rhs, ct_m, &
        reassemble_ct_m, u, p_mesh, fs)
   type(state_type), intent(inout):: state
   type(vector_field), intent(inout):: mom_rhs 
   type(block_csr_matrix), intent(inout):: ct_m
   logical, intent(in):: reassemble_ct_m
   type(vector_field), intent(in):: u
   type(mesh_type), intent(in):: p_mesh
   type(scalar_field), intent(in):: fs

   type(scalar_field), pointer:: it_fs, old_fs, density, external_density
   type(vector_field), pointer:: x
   character(len=FIELD_NAME_LEN):: bctype
   character(len=OPTION_PATH_LEN):: fs_option_path
   real:: itheta, rho0, g
   integer, dimension(:), pointer:: surface_element_list
   integer:: i, j, dens_stat, stat, external_density_stat
   logical:: variable_density, have_density, move_mesh, have_external_density

   assert(have_option(trim(fs%option_path)//"/prognostic"))

   it_fs => extract_scalar_field(state, "IteratedFreeSurface")
   old_fs => extract_scalar_field(state, "OldFreeSurface")

   move_mesh = have_option("/mesh_adaptivity/mesh_movement/free_surface")
   x => extract_vector_field(state, "Coordinate")
   ! reference density
   call get_fs_reference_density_from_options(rho0, state%option_path)
   density => extract_scalar_field(state, "Density", stat=dens_stat)
   have_density = (dens_stat==0)
   call get_option('/physical_parameters/gravity/magnitude', g, stat=stat)
   if (stat/=0) then
     FLExit("For a free surface you need gravity")
   end if

   call get_option(trim(u%option_path)//"/prognostic/temporal_discretisation/relaxation", &
                  itheta)
    
   do i=1, get_boundary_condition_count(u)
     call get_boundary_condition(u, i, type=bctype, &
         surface_element_list=surface_element_list, &
         option_path=fs_option_path)
     if (bctype=="free_surface") then
        if(have_option(trim(fs_option_path)//"/type[0]/no_normal_stress") .and. &
           have_option(trim(fs_option_path)//"/type[0]/no_normal_stress/explicit")) then

          external_density => extract_scalar_surface_field(u, i, "ExternalDensity", stat=external_density_stat)
          have_external_density = external_density_stat==0

          variable_density = have_option(trim(fs_option_path)//"/type[0]/variable_density") &
                            .and. (.not.move_mesh)
          if (variable_density.and.(.not.have_density)) then
            FLExit("Variable density free surface requires a Density field.")
          end if

          do j=1, size(surface_element_list)
            call add_boundary_integral_sele(j, surface_element_list(j))
          end do
        end if
     end if
   end do

  contains

    subroutine add_boundary_integral_sele(surface_mesh_ele, sele)
      integer, intent(in):: surface_mesh_ele, sele

      real, dimension(u%dim, face_loc(p_mesh, sele), face_loc(u, sele)) :: ct_mat_bdy
      real, dimension(face_ngi(u, sele)) :: detwei_bdy
      real, dimension(u%dim, face_ngi(u, sele)) :: normal_bdy
      real, dimension(face_ngi(fs, sele)):: delta_rho_g_quad
      integer, dimension(face_loc(u, sele)) :: u_nodes_bdy
      integer, dimension(face_loc(p_mesh, sele)) :: p_nodes_bdy
      integer:: dim

      u_nodes_bdy = face_global_nodes(u, sele)

      if(variable_density .and. have_external_density) then
        delta_rho_g_quad = g*(face_val_at_quad(density, sele)-ele_val_at_quad(external_density, surface_mesh_ele))
      else if (variable_density) then
        delta_rho_g_quad = g*face_val_at_quad(density, sele)
      else if (have_external_density) then
        delta_rho_g_quad = g*(rho0-node_val(external_density,1))
      else
        delta_rho_g_quad = g*rho0
      end if

      call transform_facet_to_physical(x, sele, &
           detwei_f=detwei_bdy, normal=normal_bdy)
                                     
      call addto(mom_rhs, u_nodes_bdy, shape_vector_rhs(face_shape(u, sele), normal_bdy, &
                                    -detwei_bdy*delta_rho_g_quad* &
                                    (itheta*face_val_at_quad(it_fs, sele) + &
                                     (1.0-itheta)*face_val_at_quad(old_fs, sele))))

      if(reassemble_ct_m) then
        p_nodes_bdy = face_global_nodes(p_mesh, sele)

        ct_mat_bdy = shape_shape_vector(face_shape(p_mesh, sele), face_shape(u, sele), &
             detwei_bdy, normal_bdy)

        do dim=1, u%dim
          ! we've integrated continuity by parts, but not yet added in the resulting
          ! surface integral - for the non-viscous free surface this is left
          ! out to enforce the kinematic bc
          call addto(ct_m, 1, dim, p_nodes_bdy, u_nodes_bdy, ct_mat_bdy(dim,:,:))
        end do
      end if

    end subroutine add_boundary_integral_sele

  end subroutine add_explicit_viscous_free_surface_integrals

  subroutine add_explicit_viscous_free_surface_integrals_cv(state, ct_m, &
        reassemble_ct_m, u, p_mesh, fs, mom_rhs)
    type(state_type), intent(inout):: state
    type(block_csr_matrix), intent(inout):: ct_m
    logical, intent(in):: reassemble_ct_m
    type(vector_field), intent(in):: u
    type(mesh_type), intent(in):: p_mesh
    type(scalar_field), intent(inout):: fs
    type(vector_field), intent(inout), optional:: mom_rhs 

    type(scalar_field), pointer:: it_fs, old_fs, density, external_density
    type(vector_field), pointer:: x
    character(len=FIELD_NAME_LEN):: bctype
    character(len=OPTION_PATH_LEN):: fs_option_path
    real:: itheta, rho0, g
    integer, dimension(:), pointer:: surface_element_list
    integer:: i, j, dens_stat, stat, external_density_stat
    logical:: variable_density, have_density, move_mesh, have_external_density

    ! information about cv faces
    type(cv_faces_type) :: cvfaces
    ! shape functions for region and surface
    type(element_type) :: x_cvbdyshape
    type(element_type) :: p_cvbdyshape
    type(element_type) :: u_cvbdyshape
    integer:: quaddegree, ele, sele

    if (.not.present(mom_rhs) .and. .not. reassemble_ct_m) return

    assert(have_option(trim(fs%option_path)//"/prognostic"))

    it_fs => extract_scalar_field(state, "IteratedFreeSurface")
    old_fs => extract_scalar_field(state, "OldFreeSurface")

    move_mesh = have_option("/mesh_adaptivity/mesh_movement/free_surface")
    x => extract_vector_field(state, "Coordinate")
    ! reference density
    call get_fs_reference_density_from_options(rho0, state%option_path)
    density => extract_scalar_field(state, "Density", stat=dens_stat)
    have_density = (dens_stat==0)
    call get_option('/physical_parameters/gravity/magnitude', g, stat=stat)
    if (stat/=0) then
      FLExit("For a free surface you need gravity")
    end if

    call get_option(trim(u%option_path)//"/prognostic/temporal_discretisation/relaxation", &
                   itheta)

    call get_option("/geometry/quadrature/controlvolume_surface_degree", &
                   quaddegree, default=1)

    cvfaces=find_cv_faces(vertices=ele_vertices(p_mesh, 1), &
                          dimension=mesh_dim(p_mesh), &
                          polydegree=p_mesh%shape%degree, &
                          quaddegree=quaddegree)

    x_cvbdyshape=make_cvbdy_element_shape(cvfaces, x%mesh%faces%shape)
    p_cvbdyshape=make_cvbdy_element_shape(cvfaces, p_mesh%faces%shape)
    u_cvbdyshape=make_cvbdy_element_shape(cvfaces, u%mesh%faces%shape)

    do i=1, get_boundary_condition_count(u)
      call get_boundary_condition(u, i, type=bctype, option_path=fs_option_path, &
          surface_element_list=surface_element_list)
      if (bctype=="free_surface") then
        if(have_option(trim(fs_option_path)//"/type[0]/no_normal_stress") .and. &
           have_option(trim(fs_option_path)//"/type[0]/no_normal_stress/explicit")) then
          external_density => extract_scalar_surface_field(u, i, "ExternalDensity", stat=external_density_stat)
          have_external_density = external_density_stat==0

          variable_density = have_option(trim(fs_option_path)//"/type[0]/variable_density") &
                            .and. (.not.move_mesh)
          if (variable_density.and.(.not.have_density)) then
            FLExit("Variable density free surface requires a Density field.")
          end if

          do j=1, size(surface_element_list)
            sele = surface_element_list(j)
            ele = face_ele(x, sele)
            call add_boundary_integral_sele(j, sele, ele)
          end do
        end if
      end if
    end do

    call deallocate(x_cvbdyshape)
    call deallocate(p_cvbdyshape)
    call deallocate(u_cvbdyshape)
    call deallocate(cvfaces)

  contains

    subroutine add_boundary_integral_sele(surface_mesh_ele, sele, ele)
      integer, intent(in):: surface_mesh_ele, sele, ele

      real, dimension(u%dim, face_loc(p_mesh, sele), face_loc(u, sele)) :: ct_mat_bdy
      real, dimension(face_ngi(u, sele)) :: detwei_bdy
      real, dimension(x_cvbdyshape%ngi) :: detwei_bdy_cv
      real, dimension(x%dim, face_ngi(x, sele)) :: normal_bdy
      real, dimension(x%dim, x_cvbdyshape%ngi) :: normal_bdy_cv
      real, dimension(face_ngi(fs, sele)):: delta_rho_g_quad
      integer, dimension(face_loc(u, sele)) :: u_nodes_bdy
      integer, dimension(face_loc(p_mesh, sele)) :: p_nodes_bdy
      real, dimension(x%dim, ele_loc(x, ele)) :: x_ele
      real, dimension(x%dim, face_loc(x, sele)) :: x_ele_bdy
      integer:: dim, gi, ggi, iloc, jloc, face
      
      u_nodes_bdy = face_global_nodes(u, sele)

      if (present(mom_rhs)) then
        call transform_facet_to_physical(x, sele, &
             detwei_f=detwei_bdy, normal=normal_bdy)
  
        if(variable_density .and. have_external_density) then
          delta_rho_g_quad = g*(face_val_at_quad(density, sele)-ele_val_at_quad(external_density, surface_mesh_ele))
        else if (variable_density) then
          delta_rho_g_quad = g*face_val_at_quad(density, sele)
        else if (have_external_density) then
          delta_rho_g_quad = g*(rho0-node_val(external_density, 1))
        else
          delta_rho_g_quad = g*rho0
        end if

        call addto(mom_rhs, u_nodes_bdy, shape_vector_rhs(face_shape(u, sele), normal_bdy, &
                                      -detwei_bdy*delta_rho_g_quad* &
                                      (itheta*face_val_at_quad(it_fs, sele) + &
                                       (1.0-itheta)*face_val_at_quad(old_fs, sele))))
      end if

      if(reassemble_ct_m) then

        p_nodes_bdy = face_global_nodes(p_mesh, sele)
        x_ele = ele_val(x, ele)
        x_ele_bdy = face_val(x, sele)

        call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, &
                              x_cvbdyshape, normal_bdy_cv, detwei_bdy_cv)
        ct_mat_bdy = 0.0
        surface_nodal_loop_i: do iloc = 1, face_loc(p_mesh, sele)

          surface_face_loop: do face = 1, cvfaces%sfaces
            if(cvfaces%sneiloc(iloc,face)/=0) then

              surface_quadrature_loop: do gi = 1, cvfaces%shape%ngi

                ggi = (face-1)*cvfaces%shape%ngi + gi

                surface_nodal_loop_j: do jloc = 1, face_loc(u, sele)

                  surface_inner_dimension_loop: do dim = 1, size(normal_bdy_cv,1)

                     ct_mat_bdy(dim, iloc, jloc) =  ct_mat_bdy(dim, iloc, jloc) + &
                           u_cvbdyshape%n(jloc,ggi)*detwei_bdy_cv(ggi)*normal_bdy_cv(dim, ggi)

                  end do surface_inner_dimension_loop

                end do surface_nodal_loop_j

              end do surface_quadrature_loop

            end if

          end do surface_face_loop

        end do surface_nodal_loop_i

        do dim=1, u%dim
          ! we've integrated continuity by parts, but not yet added in the resulting
          ! surface integral - for the non-viscous free surface this is left
          ! out to enforce the kinematic bc
          call addto(ct_m, 1, dim, face_global_nodes(p_mesh,sele), &
               face_global_nodes(u,sele), ct_mat_bdy(dim,:,:))
        end do
      end if

    end subroutine add_boundary_integral_sele

  end subroutine add_explicit_viscous_free_surface_integrals_cv

  subroutine add_implicit_viscous_free_surface_scaled_mass_integrals(state, mass, u, p_mesh, fs, dt)
    ! This routine adds in the boundary conditions for the viscous free surface
    ! (that is the free_surface bc with the no_normal_stress option)
    ! to the "scaled mass matrix" (pressure mass scaled with inverse of viscosity used as stokes preconditioner)
    ! mass has been extended with some extra rows and columns (corresponding to free surface
    ! nodes) that are used to enforce the kinematic bc.
    ! Here we fill in those terms with the scaled mass on the free surface.
    
    type(state_type), intent(inout):: state
    type(csr_matrix), intent(inout):: mass
    type(vector_field), intent(in):: u
    type(scalar_field), intent(in):: p_mesh
    type(scalar_field), intent(inout):: fs
    real, intent(in) :: dt

    type(tensor_field), pointer :: viscosity
    type(vector_field), pointer:: x
    type(scalar_field) :: viscosity_component
    type(mesh_type), pointer:: fs_mesh, embedded_fs_mesh
    type(integer_hash_table):: sele_to_fs_ele
    character(len=FIELD_NAME_LEN):: bc_type
    character(len=OPTION_PATH_LEN):: bc_option_path
    integer, dimension(:), pointer:: surface_element_list, fs_surface_element_list
    integer:: i, j

    assert(have_option(trim(fs%option_path)//"/prognostic"))

    if (.not. has_boundary_condition_name(fs, "_implicit_free_surface")) then
      call initialise_implicit_prognostic_free_surface(state, fs, u)
    end if
    ! obtain the f.s. surface mesh that has been stored under the
    ! "_implicit_free_surface" boundary condition
    call get_boundary_condition(fs, "_implicit_free_surface", &
        surface_mesh=fs_mesh, surface_element_list=fs_surface_element_list)
    ! create a map from face numbers to element numbers in fs_mesh
    call invert_set(fs_surface_element_list, sele_to_fs_ele)
    embedded_fs_mesh => extract_mesh(state, "_embedded_free_surface_mesh")

    x => extract_vector_field(state, "Coordinate")

    ! Extract viscosity tensor from state:
    viscosity => extract_tensor_field(state,'Viscosity')

    ! Extract first component of viscosity tensor from full tensor:
    viscosity_component = extract_scalar_field(viscosity,1,1)

    do i=1, get_boundary_condition_count(u)
      call get_boundary_condition(u, i, type=bc_type, &
           option_path=bc_option_path, &
           surface_element_list=surface_element_list)
      if (bc_type=="free_surface") then
        if (have_option(trim(bc_option_path)//"/type[0]/no_normal_stress") .and. &
            (.not.have_option(trim(bc_option_path)//"/type[0]/no_normal_stress/explicit"))) then
          do j=1, size(surface_element_list)
            call add_boundary_integral_sele(surface_element_list(j))
          end do
        end if
      end if
    end do

  contains

    subroutine add_boundary_integral_sele(sele)
      integer, intent(in):: sele

      real, dimension(face_loc(p_mesh, sele), face_loc(p_mesh, sele)) :: mat_bdy
      real, dimension(face_ngi(p_mesh, sele)) :: detwei_bdy
      
      call transform_facet_to_physical(x, sele, &
           detwei_f=detwei_bdy)
      mat_bdy = shape_shape(face_shape(fs, sele), face_shape(fs, sele), &
           detwei_bdy/(face_val_at_quad(viscosity_component, sele)*dt))
      call addto(mass, ele_nodes(embedded_fs_mesh, fetch(sele_to_fs_ele, sele)), &
                       ele_nodes(embedded_fs_mesh, fetch(sele_to_fs_ele, sele)), &
                       mat_bdy)

    end subroutine add_boundary_integral_sele

  end subroutine add_implicit_viscous_free_surface_scaled_mass_integrals

  
  subroutine move_mesh_free_surface(states, initialise, nonlinear_iteration)
    type(state_type), dimension(:), intent(inout) :: states
    ! if present_and_true: zero gridvelocity and compute OldCoordinate=Coordinate=IteratedCoordinate
    logical, intent(in), optional :: initialise
    ! only supply if total number nonlinear_iterations>1, in which case we do something else for the first nonlinear_iteration
    integer, intent(in), optional:: nonlinear_iteration

    type(vector_field), pointer :: velocity
    real :: itheta
    integer :: i, its

    logical :: complete
    
    if (present(nonlinear_iteration)) then
      ! we do something different in the first nonlinear iteration, see below
      its=nonlinear_iteration
    else
      ! if we don't have a non-linear loop, we do the same as
      ! in the 2nd nonlinear iteration if we would have non-linear iterations, i.e.:
      ! OldCoordinate gets set to IteratedCoordinate at the end of last timestep
      ! and we compute a new IteratedCoordinate and therefore Coordinate and GridVelocity
      its=2
    end if

    complete = .false.

    do i=1, size(states)
      velocity => extract_vector_field(states(i), "Velocity")
      
      if (aliased(velocity)) cycle
      
      if (has_boundary_condition(velocity, "free_surface") .and. &
            & have_option('/mesh_adaptivity/mesh_movement/free_surface')) then
          
        if(complete) then
          FLExit("Two velocity fields with free_surface boundary conditions are not permitted.")
        end if
        
        call get_option( trim(velocity%option_path)//'/prognostic/temporal_discretisation/relaxation', &
            itheta, default=0.5)
            
        if (its==1) then
          
          ! The first nonlinear iteration we'll keep using the OldCoordinate and IteratedCoordinate
          ! and GridVelocityfrom last timestep. Only Coordinate has been set to IteratedCoordinate 
          ! at the end of last timestep, so we need to weight it again
          call interpolate_coordinate_with_theta(states(i), itheta)
          
        else
        
          ewrite(1,*) "Going into move_free_surface_nodes to compute new node coordinates"
          
          if (its==2) then
            ! the first nonlinear iteration we've used OldCoordinate,IteratedCoordinate,GridVelocity
            ! from previous timestep. Now we recompute IteratedCoordinate and GridVelocity, based on
            ! the new free surface approx. calculated in the first nonlinear iteration. OldCoordinate
            ! should be set to IteratedCoordinate of last timestep
            call set_vector_field_in_state(states(i), "OldCoordinate", "IteratedCoordinate")
          end if
          
          call move_free_surface_nodes(states(i), itheta, initialise = initialise)
          
          ! need to update ocean boundaries again if you've just moved the mesh
          if (has_scalar_field(states(i), "DistanceToTop")) then
            if (.not. have_option('/geometry/ocean_boundaries')) then
                FLExit("ocean_boundaries required under geometry for mesh movement with a free_surface")
            end if
            call CalculateTopBottomDistance(states(i))
          end if
          call update_wettingdrying_alpha(states(i))

          
          
        end if
        
        complete = .true.
        
      end if
    end do
  
  end subroutine move_mesh_free_surface

  subroutine interpolate_coordinate_with_theta(state, theta)
  type(state_type), intent(inout) :: state
  real, intent(in):: theta
  
    type(vector_field), pointer:: positions, old_positions, iterated_positions
    
    call IncrementEventCounter(EVENT_MESH_MOVEMENT)
    
    positions => extract_vector_field(state, "Coordinate")
    old_positions => extract_vector_field(state, "OldCoordinate")
    iterated_positions => extract_vector_field(state, "IteratedCoordinate")
    call set(positions, old_positions, iterated_positions, theta)
    
  end subroutine interpolate_coordinate_with_theta
      
  subroutine move_free_surface_nodes(state, theta, initialise)
  type(state_type), intent(inout) :: state
  real, intent(in):: theta
  logical, intent(in), optional :: initialise
    
    type(vector_field), pointer:: positions, u, original_positions
    type(vector_field), pointer:: gravity_normal, old_positions, grid_u
    type(vector_field), pointer:: iterated_positions
    type(scalar_field), pointer:: p
    type(scalar_field), pointer:: topdis, bottomdis
    type(vector_field), target :: local_grid_u
    type(scalar_field), target:: fs_mapped_to_coordinate_space, local_fs
    type(scalar_field):: local_fracdis, scaled_fs
    character(len=FIELD_NAME_LEN):: bctype
    real dt 
    integer, dimension(:), allocatable:: face_nodes
    integer, dimension(:), pointer:: surface_element_list
    integer i, j, k, node, sele, stat

    ! some fields for when moving the entire mesh
    type(scalar_field), pointer :: fracdis
    type(scalar_field), pointer :: fs
    
    logical :: l_initialise, have_prognostic_fs
    
    ewrite(1,*) 'Entering move_free_surface_nodes'
    
    ! increase event counter, so position caching know the mesh has moved
    call IncrementEventCounter(EVENT_MESH_MOVEMENT)
    
    l_initialise = present_and_true(initialise)
    
    ! gravity acceleration
    call get_option('/timestepping/timestep', dt)
 
    positions => extract_vector_field(state, "Coordinate")
    original_positions => extract_vector_field(state, "OriginalCoordinate")
    iterated_positions => extract_vector_field(state, "IteratedCoordinate")
    old_positions => extract_vector_field(state, "OldCoordinate")

    gravity_normal => extract_vector_field(state, "GravityDirection")
    ! it's alright for gravity to be on a DG version of the CoordinateMesh:
    assert( face_loc(gravity_normal,1)==face_loc(positions,1) )

    u => extract_vector_field(state, "Velocity")
    p => extract_scalar_field(state, "Pressure")
    
    fs => extract_scalar_field(state, "FreeSurface", stat=stat)
    have_prognostic_fs=.false.
    if (stat==0) then
      have_prognostic_fs = have_option(trim(fs%option_path)//"/prognostic")
    else
      call allocate(local_fs, p%mesh, "LocalFreeSurface")
      fs => local_fs
    end if
    if (.not. have_prognostic_fs) then
      ! make sure the fs is up-to-date with the latest pressure values
      call calculate_diagnostic_free_surface(state, fs)
    end if

    if (.not. fs%mesh==positions%mesh) then
      call allocate(fs_mapped_to_coordinate_space, positions%mesh)
      call remap_field(fs, fs_mapped_to_coordinate_space, stat=stat)
      if(stat==REMAP_ERR_DISCONTINUOUS_CONTINUOUS) then
        ewrite(-1,*) "Just remapped from a discontinuous to a continuous field when using free_surface mesh movement."
        ewrite(-1,*) "This suggests the FreeSurface is discontinuous, which isn't supported."
        FLExit("Discontinuous pressure not permitted.")
      else if(stat/=0 .and. stat/=REMAP_ERR_UNPERIODIC_PERIODIC .and. stat/=REMAP_ERR_HIGHER_LOWER_CONTINUOUS) then
        FLAbort("Something went wrong mapping FreeSurface to the CoordinateMesh")
      end if
      ! we've allowed it to remap from periodic to unperiodic and from higher order to lower order
      
      if (associated(fs, local_fs)) then
        call deallocate(local_fs)
      end if
      fs => fs_mapped_to_coordinate_space
    end if
   
    if(.not.l_initialise) then
    ! if we're initialising then the grid velocity stays as zero
      grid_u => extract_vector_field(state, "GridVelocity")
      
      if(.not. grid_u%mesh==positions%mesh) then
        ! allocate this on the positions mesh to calculate the values
        call allocate(local_grid_u, grid_u%dim, positions%mesh, "LocalGridVelocity")
        call zero(local_grid_u)
        grid_u => local_grid_u
      end if
    end if
    
    if (have_option("/mesh_adaptivity/mesh_movement/free_surface/move_whole_mesh")) then
    
      if (.not. have_option("/geometry/ocean_boundaries")) then
        ! ensure we have ocean boundaries so the fs has been extrapolated downwards
        ewrite(-1,*) "With /mesh_adaptivity/mesh_movement/free_surface/move_whole_mesh" // &
           "you need /geometry/ocean_boundaries"
        FLExit("Missing option")
      end if
      
      topdis => extract_scalar_field(state, "DistanceToTop", stat=stat)
      bottomdis => extract_scalar_field(state, "DistanceToBottom", stat=stat)

      ! allocate a scaled free surface, that is equal to fs at the top
      ! and linearly decreases to 0 at the bottom
      call allocate(scaled_fs, fs%mesh, "ScaledFreeSurface")
      ! start by setting it equal to fs everywhere
      call set(scaled_fs, fs)

      ! Now we need to scale it by its fractional distance from the bottom
      ! Since the fractional distance is constant in time, we compute it once and save it.
      if (.not. has_scalar_field(state, "FractionalDistance")) then
        call allocate(local_fracdis, fs%mesh, "FractionalDistance")     
        call set(local_fracdis, topdis)
        call addto(local_fracdis, bottomdis)
        call invert(local_fracdis)
        call scale(local_fracdis, bottomdis)
        call insert(state, local_fracdis, name="FractionalDistance")
        call deallocate(local_fracdis)
      end if
      fracdis => extract_scalar_field(state, "FractionalDistance")

      call scale(scaled_fs, fracdis)
        
      do node=1, node_count(positions)
        call set(iterated_positions, node, &
                  node_val(original_positions, node)- &
                  node_val(scaled_fs, node)*node_val(gravity_normal, node))
                  
        if(.not.l_initialise) call set(grid_u, node, &
                  (node_val(iterated_positions, node)-node_val(old_positions,node))/dt)
      end do
      
      call deallocate(scaled_fs)
    
    else
      
      ! we assume no p-refinement on the coordinate mesh
      allocate( face_nodes(1:face_loc(positions,1)) )
      
      do i=1, get_boundary_condition_count(u)
        call get_boundary_condition(u, i, type=bctype, &
          surface_element_list=surface_element_list)
        if (bctype=="free_surface") then
          
          face_loop: do j=1, size(surface_element_list)

            sele=surface_element_list(j)
            face_nodes=face_global_nodes(positions, sele)
            
            node_loop: do k=1, size(face_nodes)
                node=face_nodes(k)
                ! compute new surface node position:
                call set(iterated_positions, node, &
                    node_val(original_positions, node)- &
                    node_val(fs, node)*node_val(gravity_normal, node))

                ! compute new surface node grid velocity:
                if(.not.l_initialise) call set(grid_u, node, &
                  (node_val(iterated_positions, node)-node_val(old_positions,node))/dt)
            end do node_loop
            
          end do face_loop
            
        end if
      end do
      
    end if
   

    if (l_initialise) then
      call set(positions, iterated_positions)
      call set(old_positions, iterated_positions)
    else
      call set(positions, iterated_positions, old_positions, theta)
      if(associated(grid_u, local_grid_u)) then
        grid_u => extract_vector_field(state, "GridVelocity")
        call remap_field(local_grid_u, grid_u)
        call deallocate(local_grid_u)
      end if
      ewrite_minmax(grid_u)
    end if
    
    if (associated(fs, fs_mapped_to_coordinate_space)) then
       call deallocate(fs_mapped_to_coordinate_space)
    else if (associated(fs, local_fs)) then
       call deallocate(local_fs)
    end if

  end subroutine move_free_surface_nodes

  function vertical_prolongator_from_free_surface(state, mesh) result (vertical_prolongator)
  !! Creates a prolongation operator from the free surface mesh to
  !! the full mesh which when provided to petsc_solve is used
  !! to implement vertical lumping if the "mg" preconditioner is selected.
  type(state_type), intent(in):: state
  type(mesh_type), intent(in):: mesh
  type(petsc_csr_matrix):: vertical_prolongator
    
    type(csr_matrix):: csr_vertical_prolongator
    type(scalar_field), pointer:: topdis
    type(vector_field), pointer:: positions, vertical_normal
    type(integer_set):: owned_surface_nodes
    integer, dimension(:), pointer:: surface_element_list, surface_node_list
    integer:: stat, i, face, ncols
    
    ewrite(1, *) "Constructing vertical_prolongator_from_free_surface to be used in mg"
    topdis => extract_scalar_field(state, "DistanceToTop", stat=stat)
    if (stat/=0) then
       FLExit("For vertical lumping you need to specify the ocean_boundaries under /geometry")
    end if

    positions => extract_vector_field(state, "Coordinate")
    vertical_normal => extract_vector_field(state, "GravityDirection")

    call get_boundary_condition(topdis, 1, &
      surface_element_list=surface_element_list, &
      surface_node_list=surface_node_list)
    

    csr_vertical_prolongator=VerticalProlongationOperator( &
         mesh, positions, vertical_normal, surface_element_list)
#ifdef DDEBUG
    ! note that in surface_positions the non-owned free surface nodes may be inbetween
    ! the reduce_columns option should have removed those however
    ! with debugging perform test to check if this is the case:
    if (IsParallel()) then
      assert( associated(mesh%halos) )      
      ! count n/o owned surface nodes
      call allocate(owned_surface_nodes)
      do i=1, size(surface_element_list)
        face=surface_element_list(i)
        call insert(owned_surface_nodes, face_global_nodes(mesh, face))
      end do
      ! this should be size(csr_vertical_prolongator,2) but intel with debugging
      ! seems to choke on it:
      ncols = csr_vertical_prolongator%sparsity%columns
      ewrite(2,*) "Number of owned surface nodes:", key_count(owned_surface_nodes)
      ewrite(2,*) "Number of columns in vertical prolongator:", ncols
      if (ncols>key_count(owned_surface_nodes)) then
        ewrite(-1,*) "Vertical prolongator seems to be using more surface nodes than the number"
        ewrite(-1,*) "of surface nodes within completely owned surface elements. This indicates"
        ewrite(-1,*) "the parallel decomposition is not done along columns. You shouldn't be using"
        ewrite(-1,*) "mg with vertical_lumping in that case."
        FLExit("Vertical lumping requires 2d decomposition along columns")
      end if
      call deallocate(owned_surface_nodes)
    end if
#endif

    vertical_prolongator=csr2petsc_csr(csr_vertical_prolongator)
    call deallocate(csr_vertical_prolongator)

  end function vertical_prolongator_from_free_surface
  
  function free_surface_nodes(state, mesh)
  !! Returns the list of the nodes on the free surface of the given mesh.
  !! Returns a pointer to an allocated array to be deallocated by the caller.
  integer, dimension(:), pointer:: free_surface_nodes
  type(state_type), intent(in):: state
  type(mesh_type), intent(in):: mesh
    
    type(scalar_field), pointer:: topdis
    type(mesh_type) surface_mesh
    integer, dimension(:), pointer:: surface_element_list, surface_node_list
    integer stat
    
    ewrite(1, *) "Extracting list of nodes on the free surface"
    topdis => extract_scalar_field(state, "DistanceToTop", stat=stat)
    if (stat/=0) then
       FLExit("Need to specify the ocean_boundaries under /geometry")
    end if
    
    if (mesh==topdis%mesh) then
      ! we can just copy this info, from the coordinate mesh
      call get_boundary_condition(topdis, 1, &
        surface_node_list=surface_node_list)
      allocate( free_surface_nodes(1: size(surface_node_list)) )
      free_surface_nodes=surface_node_list
    else
      ! by creating a temporary surface mesh we get exactly this information
      call get_boundary_condition(topdis, 1, &
         surface_element_list=surface_element_list)
      call create_surface_mesh(surface_mesh, surface_node_list, &
         mesh, surface_element_list, name=trim(mesh%name)//'FreeSurface')
      free_surface_nodes => surface_node_list
      call deallocate(surface_mesh)
    end if
  
  end function free_surface_nodes
  
  subroutine calculate_diagnostic_free_surface(state, free_surface)
  !!< calculates a 3D field (constant over the vertical) of the free surface elevation
  !!< This can be added as a diagnostic field in the flml.
  type(state_type), intent(inout):: state
  type(scalar_field), target, intent(inout):: free_surface
    
     integer, dimension(:), pointer:: surface_element_list
     type(vector_field), pointer:: x, u, vertical_normal
     type(scalar_field), pointer:: p, topdis, external_density
     type(scalar_field), pointer :: original_bottomdist_remap
     character(len=OPTION_PATH_LEN):: fs_option_path
     type(scalar_field) :: p_min
     type(scalar_field), target :: p_capped
     character(len=FIELD_NAME_LEN):: bctype
     real:: g, rho0, delta_rho, d0, p_atm 
     integer:: i, j, sele, stat, external_density_stat
     logical :: have_wd, have_external_density

     ! the prognostic free surface is calculated elsewhere (this is used
     ! in combination with the viscous free surface)
     if (have_option(trim(free_surface%option_path)//'/prognostic')) return
     
     x => extract_vector_field(state, "Coordinate")
     p => extract_scalar_field(state, "Pressure")
     assert(free_surface%mesh==p%mesh)

     call get_option('/physical_parameters/gravity/magnitude', g)
     call get_option(trim(p%option_path)//'/prognostic/atmospheric_pressure', &
       p_atm, default=0.0)

     have_wd=have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying")

     u => extract_vector_field(state, "Velocity")
     call get_fs_reference_density_from_options(rho0, state%option_path)

     if (have_option(trim(u%option_path)//'/prognostic/equation::ShallowWater')) then
       ! For the shallow water equations we only have a 2D, horizontal mesh
       ! and f.s. is simply p/g everywhere:
       call set(free_surface, p)
       call scale(free_surface, 1./g)
       return
     end if

     !
     ! first we compute the right free surface values at the free surface
     ! nodes only
     !

     ! Do the wetting and drying corrections: 
     ! In dry regions, the free surface is not coupled to the pressure but is fixed to -OriginalCoordinate+d0.
     ! Hence, we create temporary pressure that is capped to d0-bottom_depth on the surface before extruding it downwards.
     if (have_wd) then
        call allocate(p_min, p%mesh, "MinimumSurfacePressure")
        call allocate(p_capped, p%mesh, "CappedPressure")
        call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/d0", d0)
        original_bottomdist_remap=>extract_scalar_field(state, "OriginalDistanceToBottomPressureMesh")

        ! We are looking for p_capped = min(p, g*rho0*(d0-bottom_depth)) on the surface
        call set(p_min, original_bottomdist_remap)
        call addto(p_min, -d0)
        call scale(p_min, -g*rho0)
        call set(p_capped, p)
        call bound(p_capped, lower_bound=p_min)
        p=>p_capped

        call deallocate(p_min)
     end if      

     ! make sure other nodes are zeroed
     call zero(free_surface)    

     do i=1, get_boundary_condition_count(u)
        call get_boundary_condition(u, i, type=bctype, &
           surface_element_list=surface_element_list, &
           option_path=fs_option_path)
        if (bctype=="free_surface") then

          if (have_option(trim(fs_option_path)//"/type[0]/variable_density")) then
            ! options checked below with an flexit
            FLAbort("Cannot use a diagnostic free surface field with a variable density free surface bc")
          end if
          external_density => extract_scalar_surface_field(u, i, "ExternalDensity", stat=external_density_stat)
          have_external_density = external_density_stat==0

          if (size(surface_element_list)==0) cycle

          if (have_external_density) then
            delta_rho=rho0-node_val(external_density, 1)
          else
            delta_rho=rho0
          end if
      
          face_loop: do j=1, size(surface_element_list)

            sele=surface_element_list(j)
             
            call set(free_surface, &
               face_global_nodes(free_surface, sele), &
               (face_val(p, sele)-p_atm)/delta_rho/g)
            if (have_wd) then
              ! bound free surface from below by -orig_bottomdist+d0
              call set(free_surface, face_global_nodes(free_surface, sele), &
                max(face_val(free_surface, sele), -face_val(original_bottomdist_remap, sele)))
            end if
             
          end do face_loop
             
        end if
        
     end do

     topdis => extract_scalar_field(state, "DistanceToTop", stat=stat)
     ! if /geometry/ocean_boundaries are not specified we leave at this
     if (stat/=0) return

     !
     ! otherwise we continue to extrapolate vertically from the free
     ! surface values calculated above
     !


     ! note we're not using the actual free_surface bc here, as 
     ! that may be specified in parts, or not cover the whole area
     call get_boundary_condition(topdis, 1, &
         surface_element_list=surface_element_list)
     vertical_normal => extract_vector_field(state, "GravityDirection")
 
     ! vertically extrapolate pressure values at the free surface downwards
     ! (reuse projected horizontal top surface mesh cached under DistanceToTop)
     call VerticalExtrapolation(free_surface, free_surface, x, &
       vertical_normal, surface_element_list=surface_element_list, &
       surface_name="DistanceToTop")
       
     if (have_wd) then
       call deallocate(p_capped)
     end if

  end subroutine calculate_diagnostic_free_surface

  subroutine update_wettingdrying_alpha(state)
  !!< calculates and updates the alpha coefficients for wetting and drying.
  type(state_type), intent(in):: state
  type(scalar_field), pointer:: scalar_surface_field

  integer, dimension(:), pointer :: surface_element_list
  type(vector_field), pointer:: u
  type(scalar_field), pointer:: p, original_bottomdist_remap
  character(len=FIELD_NAME_LEN):: bctype
  character(len=OPTION_PATH_LEN) fs_option_path
  real:: rho0, g, d0
  integer:: i, j, sele
  real, dimension(:), allocatable :: alpha


  u => extract_vector_field(state, "Velocity")
  p => extract_scalar_field(state, "Pressure")
  allocate(alpha(face_loc(p, 1)))

  do i=1, get_boundary_condition_count(u)
       call get_boundary_condition(u, i, type=bctype, &
          surface_element_list=surface_element_list, option_path=fs_option_path)
       if (bctype=="free_surface" .and. has_scalar_surface_field(u, i, "WettingDryingAlpha")) then
             scalar_surface_field => extract_scalar_surface_field(u, i, "WettingDryingAlpha")
             ! Update WettingDryingAlpha
             call get_fs_reference_density_from_options(rho0, state%option_path)
             original_bottomdist_remap => extract_scalar_field(state, "OriginalDistanceToBottomPressureMesh") 
             call get_option('/physical_parameters/gravity/magnitude', g)
             call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/d0", d0)

             ! Calculate alpha for each surface element
             face_loop: do j=1, size(surface_element_list)
                sele=surface_element_list(j)
                call calculate_alpha(sele, alpha)
                call set(scalar_surface_field, &
                     ele_nodes(scalar_surface_field, j), &
                     alpha)
             end do face_loop
       end if
   end do
   deallocate(alpha)

  contains
   subroutine calculate_alpha(sele, alpha)
        integer, intent(in) :: sele
        real, dimension(face_loc(p, sele)), intent(inout) :: alpha
        integer :: i
        
        alpha = g * face_val(original_bottomdist_remap, sele) -g * d0 + face_val(p, sele)
        alpha = alpha/g * (-d0)
        do i=1, size(alpha)
          if (alpha(i)<=0.0) then
              alpha(i)=0.0
          else
              alpha(i)=1.0
          end if
        end do
    end subroutine calculate_alpha
  end subroutine update_wettingdrying_alpha


  subroutine calculate_diagnostic_wettingdrying_alpha(state, wettingdrying_alpha)
    !!< calculates the alpha coefficient of the wetting and drying algorithm.
    !!< This can be added as a diagnostic field in the flml.
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: wettingdrying_alpha
    type(scalar_field), pointer :: scalar_surface_field

    integer, dimension(:), pointer :: surface_element_list
    type(vector_field), pointer :: x, u, gravity_normal
    type(scalar_field), pointer :: topdis
    character(len=FIELD_NAME_LEN):: bctype
    integer :: i, j, sele, stat

    u => extract_vector_field(state, "Velocity")
    if (.not. wettingdrying_alpha%mesh==extract_pressure_mesh(state)) then
        FLExit("The WettingDryingAlpha diagnostic field must live on the PressureMesh.")
    end if
  
    call zero(wettingdrying_alpha)
    do i = 1, get_boundary_condition_count(u)
        call get_boundary_condition(u, i, type=bctype, &
                surface_element_list = surface_element_list)
        if (bctype=="free_surface") then
            scalar_surface_field => extract_scalar_surface_field(u, i, "WettingDryingAlpha")
            do j = 1, size(surface_element_list)
                sele = surface_element_list(j)
                call set(wettingdrying_alpha, &
                face_global_nodes(wettingdrying_alpha, sele), &
                ele_val(scalar_surface_field, j))
            end do
        end if
    end do
    ! Extrapolate values down the horizontal if possible.
    topdis => extract_scalar_field(state, "DistanceToTop", stat = stat)
    if (stat == 0) then
        ! note we're not using the actual wettingdrying_alpha bc here, as
        ! that may be specified in parts, or not cover the whole area
        call get_boundary_condition(topdis, 1, &
        surface_element_list = surface_element_list)
        gravity_normal => extract_vector_field(state, "GravityDirection")
        x => extract_vector_field(state, "Coordinate")

        ! vertically extrapolate pressure values at the free surface downwards
        ! (reuse projected horizontal top surface mesh cached under DistanceToTop)
        call VerticalExtrapolation(wettingdrying_alpha, wettingdrying_alpha, x, &
        gravity_normal, surface_element_list = surface_element_list)
    end if

  end subroutine calculate_diagnostic_wettingdrying_alpha


  function calculate_volume_by_surface_integral(state) result(volume)
      type(state_type), intent(in) :: state
      real :: dt
      
      type(vector_field), pointer:: positions, u, gravity_normal
      type(scalar_field), pointer:: p, original_bottomdist_remap
      character(len=FIELD_NAME_LEN):: bctype
      character(len=OPTION_PATH_LEN) :: fs_option_path
      real:: g, rho0, alpha, volume, d0, delta_rho, external_density
      integer, dimension(:), pointer:: surface_element_list
      integer:: i, j, grav_stat
      logical:: include_normals, move_mesh
      logical:: have_wd
      
      ! gravity acceleration
      call get_option('/physical_parameters/gravity/magnitude', g, stat=grav_stat)
        
      ! get the pressure, and the pressure at the beginning of the time step
      p => extract_scalar_field(state, "Pressure")
      u => extract_vector_field(state, "Velocity")
      have_wd=have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying")
      if (have_wd) then
        call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/d0", d0)
        ! original_bottomdist is needed on the pressure mesh
        original_bottomdist_remap=>extract_scalar_field(state, "OriginalDistanceToBottomPressureMesh")
      end if
      
      ! reference density
      call get_fs_reference_density_from_options(rho0, state%option_path)
      ! Timestep
      call get_option('/timestepping/timestep', dt)
      move_mesh = have_option("/mesh_adaptivity/mesh_movement/free_surface")
      ! only include the inner product of gravity and surface normal
      ! if the free surface nodes are actually moved (not necessary
      ! for large scale ocean simulations) - otherwise the free surface is assumed flat
      include_normals = move_mesh
      if (include_normals) then
        gravity_normal => extract_vector_field(state, "GravityDirection")
      end if
      positions => extract_vector_field(state, "Coordinate")
      
      volume=0.0
      do i=1, get_boundary_condition_count(u)
        call get_boundary_condition(u, i, type=bctype, &
           surface_element_list=surface_element_list,option_path=fs_option_path)
        if (bctype=="free_surface") then
          call get_option(trim(fs_option_path)//"/type[0]/external_density", &
             external_density, default=0.0)
          delta_rho=rho0-external_density
          alpha=1.0/g/delta_rho/dt
          do j=1, size(surface_element_list)
            volume=volume+calculate_volume_by_surface_integral_element(surface_element_list(j))
          end do
        end if
      end do

    contains 

    function calculate_volume_by_surface_integral_element(sele) result(volume)
        integer, intent(in) :: sele
        integer :: i

        real, dimension(positions%dim, face_ngi(positions, sele)):: normals
        real, dimension(face_loc(p, sele), face_loc(p, sele)):: mass_ele, mass_ele_wd
        real, dimension(face_ngi(p, sele)):: detwei, alpha_wetdry_quad 
        real, dimension(face_loc(p, sele)) :: one
        real :: volume

        one = 1.0

        if(include_normals) then
          call transform_facet_to_physical(positions, sele, detwei_f=detwei,&
              & normal=normals)
          ! at each gauss point multiply with inner product of gravity and surface normal
          detwei=detwei*(-1.0)*sum(face_val_at_quad(gravity_normal,sele)*normals, dim=1)
        else
          call transform_facet_to_physical(positions, sele, detwei_f=detwei)
        end if
        

        if (have_wd) then
            ! Calculate alpha_wetdry_quad. The resulting array is 0 if the quad point is wet (p > -g d_0)  
            ! and 1 if the quad point is dry (p <= -g d_0)
            alpha_wetdry_quad = -face_val_at_quad(p, sele)-face_val_at_quad(original_bottomdist_remap, sele)*g + d0 * g
            do i=1, size(alpha_wetdry_quad)
                if (alpha_wetdry_quad(i)>0.0)  then
                    alpha_wetdry_quad(i)=1.0
                else
                    alpha_wetdry_quad(i)=0.0
                end if
            end do
        end if
        
        if (have_wd) then
          mass_ele=shape_shape(face_shape(p, sele), face_shape(p, sele), detwei*(1.0-alpha_wetdry_quad))
          mass_ele_wd=shape_shape(face_shape(p, sele), face_shape(p, sele), detwei*alpha_wetdry_quad)
        else
          mass_ele=shape_shape(face_shape(p, sele), face_shape(p, sele), detwei)
          mass_ele_wd=0.0
        end if
        
        volume=dot_product(one,matmul(mass_ele, face_val(original_bottomdist_remap, sele)) &
          +matmul(mass_ele_wd, face_val(original_bottomdist_remap, sele)))
   
        volume=volume-dt*dot_product(one,-matmul(mass_ele, face_val(p, sele))*alpha &
          +matmul(mass_ele_wd, face_val(original_bottomdist_remap, sele)-d0)*alpha*g)
       
    end function calculate_volume_by_surface_integral_element
  end function calculate_volume_by_surface_integral

  subroutine free_surface_module_check_options
    
    character(len=OPTION_PATH_LEN):: option_path, phase_path, pressure_path, pade_path, bc_option_path
    character(len=FIELD_NAME_LEN):: fs_meshname, p_meshname, bctype
    logical:: have_free_surface, have_explicit_free_surface, have_viscous_free_surface
    logical::  have_standard_free_surface, have_variable_density
    logical:: local_have_explicit_free_surface, local_have_viscous_free_surface
    logical:: have_wd, have_swe
    integer i, p

    do p=1, option_count('/material_phase')
      phase_path='/material_phase['//int2str(p-1)//']'
      pressure_path=trim(phase_path)//'/scalar_field::Pressure/prognostic'      
      ! check if we have a free_surface bc
      option_path=trim(phase_path)//'/vector_field::Velocity/prognostic'
      if (have_option(trim(option_path))) then
        have_free_surface=.false.
        have_viscous_free_surface=.false.
        have_explicit_free_surface=.false.
        have_standard_free_surface=.false.
        have_variable_density=.false.
        do i=1, option_count(trim(option_path)//'/boundary_conditions')
          call get_option(trim(option_path)//'/boundary_conditions['// &
             int2str(i-1)//']/type[0]/name', bctype)
          if (bctype=='free_surface') then
            have_free_surface=.true.
            bc_option_path = trim(option_path)//'/boundary_conditions['// &
                int2str(i-1)//']/type[0]'
            local_have_viscous_free_surface = have_option(trim(bc_option_path)//'/no_normal_stress').and. &
                .not.have_option(trim(bc_option_path)//'/no_normal_stress/explicit')
            local_have_explicit_free_surface = have_option(trim(bc_option_path)//'/no_normal_stress').and. &
                have_option(trim(bc_option_path)//'/no_normal_stress/explicit')
            have_standard_free_surface=have_standard_free_surface.or. &
                ((.not.local_have_viscous_free_surface).and.(.not.local_have_explicit_free_surface))
            have_viscous_free_surface=have_viscous_free_surface.or.local_have_viscous_free_surface
            have_explicit_free_surface=have_explicit_free_surface.or.local_have_explicit_free_surface
            have_variable_density = have_variable_density .or. have_option(trim(bc_option_path)//'/variable_density')

            if (have_option(trim(bc_option_path)//'/external_density') .and. &
              .not. have_option(trim(bc_option_path)//'/external_density/constant') .and. &
              .not. have_option(trim(bc_option_path)//'/variable_density')) then
              ewrite(-1,*) "Under the free surface boundary condition at "//trim(bc_option_path)
              FLExit("With a non-constant external density you also need the variable_density option")
            end if
          end if
        end do
      else
        ! no prognostic velocity, no free_surface bc
        have_free_surface=.false.
        have_standard_free_surface=.false.
        have_explicit_free_surface=.false.
        have_viscous_free_surface = .false.
      end if
      have_wd=have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying")

      have_swe=have_option(trim(option_path)//'/equation::ShallowWater')
      
      if (have_standard_free_surface) then
         ewrite(2,*) "You have a standard free surface boundary condition, checking its options"
      end if
      
      if (have_explicit_free_surface) then
         ewrite(2,*) "You have an explicit viscous free surface boundary condition, checking its options"
      end if
      
      if (have_viscous_free_surface) then
         ewrite(2,*) "You have an implicit viscous free surface boundary condition, checking its options"
      end if
      
      ! first check we're using the new code path (cg_test or dg)
      if (have_free_surface .and. .not. have_option(trim(option_path)// &
           '/spatial_discretisation/continuous_galerkin') .and. &
           .not. have_option(trim(option_path)// &
           '/spatial_discretisation/discontinuous_galerkin')) then
         ewrite(-1,*) "With the free_surface boundary condition"
         FLExit("you have to use continuous_galerkin or discontinuous_galerkin Velocity")
      end if

      ! check pressure options
      ! first check we have a progn. pressure at all
      if (have_free_surface .and. .not. have_option(pressure_path)) then
         ewrite(-1,*) "With the free_surface boundary condition"
         FLExit("You need a prognostic pressure")
      end if

      if (have_standard_free_surface .and. .not. have_option(trim(pressure_path)// &
        '/spatial_discretisation/continuous_galerkin')) then
        ewrite(-1,*) "With standard free_surface boundary condition"
        FLExit("only a continuous_galerkin spatial_discretisation works for Pressure.")
      end if
      
      if (have_standard_free_surface .and. have_option(trim(pressure_path)// &
        '/spatial_discretisation/continuous_galerkin/test_continuity_with_cv_dual')) then
        ewrite(-1,*) "With standard free_surface boundary condition"
        FLExit("you cannot use test_continuity_with_cv_dual under Pressure.")
      end if
      
      if (have_free_surface .and. .not. have_option(trim(pressure_path)// &
        '/spatial_discretisation/continuous_galerkin/integrate_continuity_by_parts')) then
         ewrite(-1,*) "With the free_surface boundary condition"
         FLExit("you have to use the integrate_continuity_by_parts option under Pressure")
      end if
      
      ! check diagnostic FreeSurface options:
      option_path=trim(phase_path)//'/scalar_field::FreeSurface/diagnostic'
      if (have_option(trim(option_path))) then
        call get_option(trim(option_path)//'/mesh[0]/name', fs_meshname)
        call get_option(trim(pressure_path)//'/mesh[0]/name', p_meshname)
        if (.not. (have_free_surface .or. have_swe)) then
          ewrite(-1,*) "The diagnostic FreeSurface field has to be used in combination " // &
            "with the free_surface boundary condition under Velocity, or with " // &
            "equation type ShallowWater for Velocity."
          FLExit("Exit")
        end if
        if (.not. fs_meshname==p_meshname) then
          FLExit("The diagnostic FreeSurface field and the Pressure field have to be on the same mesh")
        end if
        if (.not. have_option('/geometry/ocean_boundaries') .and. .not. have_swe) then
          ewrite(0,*) "Warning: your diagnostic free surface will only be " // &
            "defined at the free surface nodes and not extrapolated downwards, " // &
            "because you didn't specify geometry/ocean_boundaries."
        end if
        if (have_variable_density) then
          FLExit("The diagnostic FreeSurface field cannot be used in combination with a variable density")
        end if
      end if
      
      ! check prognostic FreeSurface options:
      option_path=trim(phase_path)//'/scalar_field::FreeSurface/prognostic'
      if ((have_viscous_free_surface.or.have_explicit_free_surface) .and. .not. have_option(option_path)) then
        FLExit("For a free surface with no_normal_stress you need a prognostic FreeSurface")
      end if
      if (have_option(trim(option_path))) then
        call get_option(trim(option_path)//'/mesh[0]/name', fs_meshname)
        call get_option(trim(pressure_path)//'/mesh[0]/name', p_meshname)
        if (.not. have_viscous_free_surface .and. .not. have_explicit_free_surface) then
          ewrite(-1,*) "The prognostic FreeSurface field has to be used in combination " // &
            "with the free_surface boundary condition under Velocity with the." //&
            "no_normal_stress option underneath it."
          FLExit("Exit")
        end if
        if (.not. fs_meshname==p_meshname) then
          FLExit("The prognostic FreeSurface field and the Pressure field have to be on the same mesh")
        end if
        if (.not. have_option('/geometry/ocean_boundaries')) then
          ewrite(0,*) "Warning: your prognostic free surface will only be " // &
            "defined at the free surface nodes and not extrapolated downwards, " // &
            "because you didn't specify geometry/ocean_boundaries."
        end if
      end if
      
      option_path=trim(phase_path)//'/equation_of_state/fluids/linear/subtract_out_hydrostatic_level'
      pade_path=trim(phase_path)//'/equation_of_state/fluids/ocean_pade_approximation'
      if (have_free_surface .and. .not.(have_option(option_path)) .and. .not.(have_option(pade_path))) then
        ewrite(-1,*) "Missing option: ", trim(option_path)
        FLExit("With the free surface you need to subtract out the hydrostatic level.")
      end if
      
      option_path=trim(phase_path)//'/scalar_field::Pressure/prognostic/reference_node'
      if (have_free_surface .and. have_option(option_path)) then
        FLExit("With the free surface you shouldn't set a reference node for Pressure")
      end if
      
      option_path=trim(phase_path)//'/scalar_field::Pressure/prognostic/solver/remove_null_space'
      if (have_free_surface .and. have_option(option_path)) then
        FLExit("With the free surface you shouldn't set remove the null space ")
      end if

      if ((have_viscous_free_surface .or. have_explicit_free_surface) .and. have_wd) then
        ! feel free to try and add a test case
        ! if you find it indeed doesn't work please change to FLExit
        ewrite(0,*) "WARNING: the combination no_normal_stress under and wetting and drying is completely untested."
      end if
    end do
    
  end subroutine free_surface_module_check_options

end module free_surface_module

