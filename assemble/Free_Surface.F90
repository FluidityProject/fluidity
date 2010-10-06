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

module free_surface_module
use fields
use state_module
use sparse_tools
use sparse_matrices_fields
use boundary_conditions
use spud
use vertical_extrapolation_module
use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
use parallel_tools
use halos
use eventcounter
use integer_set_module
use field_options
implicit none

private

public move_mesh_free_surface, add_free_surface_to_cmc_projection, &
  vertical_prolongator_from_free_surface, &
  free_surface_nodes, calculate_diagnostic_free_surface, &
  add_free_surface_to_poisson_rhs, copy_poisson_solution_to_interior, &
  calculate_diagnostic_wettingdrying_alpha, initialize_wetting_and_drying
  

public free_surface_module_check_options

contains

  subroutine initialize_wetting_and_drying(state)
     ! This needs to be called before the first timestep in order to get the correct d0 fields,
     ! since their calculation is based on the initial(!) DistanceToBottom field.
     type(state_type), intent(in) :: state
     type(scalar_field), pointer :: d0
     type(vector_field), pointer :: u
     logical :: have_wd
     
     u => extract_vector_field(state, "Velocity")
     have_wd=have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying")
     if (have_wd) then
       d0=>get_wettingdrying_d0(state)
     end if
  end subroutine initialize_wetting_and_drying
  
  subroutine add_free_surface_to_cmc_projection(state, cmc, dt, theta, get_cmc, rhs)
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
  !!<   M (u^n+1-u^n) - dt C p^{n+theta} + ... = 0
  !!< We're solving the continuity equation:
  !!<   C^T u^{n+theta}+ M_fs p^{n+1}-p^n = 0
  !!< which leads to a projection equation of:
  !!<   theta^2 C^T M^-1 C dt dp + alpha M_fs dp = 
  !!<    -theta C^T u* - (1-theta) C^T u^n -alpha M_fs (p*-p^n)
  !!< where M_fs is the free surface integral of M_i M_j
  !!< and alpha=1/(g dt)
  
    type(state_type), intent(inout) :: state
    type(csr_matrix), intent(inout) :: cmc
    real, intent(in) :: dt, theta
    !! only add in to the matrix if get_cmc==.true.
    logical, intent(in):: get_cmc
    type(scalar_field), optional, intent(inout) :: rhs      
    
      type(vector_field), pointer:: positions, u, gravity_normal, old_positions
      type(scalar_field), pointer:: p, prevp, d0, d0_pfield
      character(len=FIELD_NAME_LEN):: bctype
      character(len=OPTION_PATH_LEN) :: fs_option_path
      real:: g, rho0, alpha, coef
      integer, dimension(:), pointer:: surface_element_list
      integer:: i, j, grav_stat
      logical:: include_normals, move_mesh
      logical:: addto_cmc
      logical:: have_wd
      
      real, save :: coef_old = 0.0
      
      ewrite(1,*) 'Entering add_free_surface_to_cmc_projection routine'
      ewrite(2,*) "Are we adding free-surface contribution to RHS:",present(rhs)

      ! gravity acceleration
      call get_option('/physical_parameters/gravity/magnitude', g, stat=grav_stat)
        
      ! if we've assembled cmc from scratch then this should be zeroed
      ! so that the addition to the matrix isn't incremental
      if(get_cmc) coef_old = 0.0
      ! if we haven't assembled cmc from scratch then the addition involves
      ! the change in the timestep (i.e. it increments from the previous timestep)
      
      ! get the pressure, and the pressure at the beginning of the time step
      p => extract_scalar_field(state, "Pressure")
      prevp => extract_scalar_field(state, "OldPressure")
      u => extract_vector_field(state, "Velocity")
      
      ! reference density
      call get_reference_density_from_options(state, u, rho0)

      move_mesh = have_option("/mesh_adaptivity/mesh_movement/free_surface")
      ! only include the inner product of gravity and surface normal
      ! if the free surface nodes are actually moved (not necessary
      ! for large scale ocean simulations) - otherwise the free surface is assumed flat
      include_normals = move_mesh
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
      
      alpha=1.0/g/rho0/dt
      coef = alpha/(theta**2*dt)
      addto_cmc = .false. ! assume we don't need to add to cmc
      ! but it will be set to true if we have adaptive timestepping turned
      ! on and that results in a different coefficient
      do i=1, get_boundary_condition_count(u)
        call get_boundary_condition(u, i, type=bctype, &
           surface_element_list=surface_element_list,option_path=fs_option_path)
        if (bctype=="free_surface") then
          if (grav_stat/=0) then
             FLExit("For a free surface you need gravity")
          end if
          have_wd=have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying")
          if (have_wd) then
            if (.not. get_cmc) then
                 FLExit("Wetting and drying needs to be reassembled at each timestep at the moment!")
            end if
            ! d0 is needed on the pressure mesh
            d0=>get_wettingdrying_d0(state)
            allocate(d0_pfield)
            call allocate(d0_pfield, p%mesh, "d0_pfield")
            call remap_field(d0, d0_pfield)
          end if
          ! only add to cmc if we're reassembling it (get_cmc = .true.)
          ! or if the timestep has changed (using adaptive timestepping and coef/=coef_old)
          addto_cmc = get_cmc.or.&
                      (have_option("/timestepping/adaptive_timestep").and.(coef/=coef_old))
          do j=1, size(surface_element_list)
            call add_free_surface_element(surface_element_list(j), have_wd, d0_pfield)
          end do
        end if
      end do
      if(addto_cmc) then
        ! cmc has been modified (most likely by changing the timestep)
        ! therefore we need to invalidate the solver context
        call destroy_solver_cache(cmc)
      end if
      
      ! save the coefficient with the current timestep for the next time round
      coef_old = coef
    
      if(present(rhs)) then
         ewrite_minmax(rhs)
      end if
      if (have_wd) then
            call deallocate(d0_pfield)
            deallocate(d0_pfield)
      end if
      
    contains 
    
    subroutine add_free_surface_element(sele,  have_wd, d0_pfield)
      integer, intent(in):: sele
      logical :: have_wd
      type(scalar_field), pointer:: d0_pfield
      integer :: i

      real, dimension(positions%dim, face_ngi(positions, sele)):: normals
      real, dimension(face_loc(p, sele), face_loc(p, sele)):: mass_ele, mass_ele_wd, mass_ele_old, mass_ele_old_wd
      real, dimension(face_ngi(p, sele)):: detwei, alpha_wetdry_quad, alpha_wetdry_quad_prevp
      real :: smoothmax
      
      if(include_normals) then
        call transform_facet_to_physical(positions, sele, detwei_f=detwei,&
            & normal=normals)
        ! at each gauss point multiply with inner product of gravity and surface normal
        detwei=detwei*(-1.0)*sum(face_val_at_quad(gravity_normal,sele)*normals, dim=1)
      else
        call transform_facet_to_physical(positions, sele, detwei_f=detwei)
      end if
      

      if (have_wd) then
          assert(coef_old==0.0)  ! wetting and drying needs to reassemble constantly at the moment.
          if (have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/smooth_max_discretisation")) then
              ! Calculate alpha_wetdry_quad. The resulting array is 0 if the quad point is wet (p > -g (d_0-smootherof)) and increases linearly to 1 if the quad point reaches p <= -g d_0
              ! Note: smooth_max_discretisation is deprectated and will be removed soon.
             call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/smooth_max_discretisation", smoothmax)
             alpha_wetdry_quad = -face_val_at_quad(p, sele)-(face_val_at_quad(d0_pfield, sele)-smoothmax)*g
             alpha_wetdry_quad = alpha_wetdry_quad/(smoothmax*g)
             alpha_wetdry_quad = max(min(alpha_wetdry_quad, 1.0), 0.0)
             alpha_wetdry_quad_prevp = alpha_wetdry_quad
          else
              ! Calculat alpha_wetdry_quad. The resulting array is 0 if the quad point is wet (p > -g d_0) and 1 if the quad point is dry (p <= -g d_0)
             alpha_wetdry_quad = -face_val_at_quad(p, sele)-face_val_at_quad(d0_pfield, sele)*g
             do i=1, size(alpha_wetdry_quad)
                 if (alpha_wetdry_quad(i)>0.0)  then
                     alpha_wetdry_quad(i)=1.0
                 else
                     alpha_wetdry_quad(i)=0.0
                 end if
             end do
             ! Calculat alpha_wetdry_quad_prev. Same as alpha_wetdry_quad but for the previous pressure.
             alpha_wetdry_quad_prevp = -face_val_at_quad(prevp, sele)-face_val_at_quad(d0_pfield, sele)*g
             do i=1, size(alpha_wetdry_quad_prevp)
                 if (alpha_wetdry_quad_prevp(i)>0.0)  then
                     alpha_wetdry_quad_prevp(i)=1.0
                 else
                     alpha_wetdry_quad_prevp(i)=0.0
                 end if
             end do
          end if

      end if
      
      if (have_wd) then
        mass_ele=shape_shape(face_shape(p, sele), face_shape(p, sele), detwei*(1.0-alpha_wetdry_quad))
        mass_ele_wd=shape_shape(face_shape(p, sele), face_shape(p, sele), detwei*alpha_wetdry_quad)
      else
        mass_ele=shape_shape(face_shape(p, sele), face_shape(p, sele), detwei)
      end if

      if (addto_cmc) then
        ! we consider the projection equation to solve for 
        ! phi=theta^2 dt dp, so that the f.s. integral
        ! alpha M_fs dp=alpha M_fs phi/(theta^2 dt)
        call addto(cmc, &
          face_global_nodes(p, sele), face_global_nodes(p,sele), &
          (coef-coef_old)*mass_ele)
      end if
      if (move_mesh) then
        ! detwei and normals at the begin of the time step
        call transform_facet_to_physical(old_positions, sele, detwei_f=detwei,&
           & normal=normals)
        ! at each gauss point multiply with inner product of gravity and surface normal
        detwei=detwei*(-1.0)*sum(face_val_at_quad(gravity_normal,sele)*normals, dim=1)
        if (have_wd) then
             mass_ele_old=shape_shape(face_shape(prevp, sele), face_shape(prevp, sele), detwei*(1.0-alpha_wetdry_quad_prevp))
             mass_ele_old_wd=shape_shape(face_shape(prevp, sele), face_shape(prevp, sele), detwei*alpha_wetdry_quad_prevp)
        else
             mass_ele_old=shape_shape(face_shape(p, sele), face_shape(p, sele), detwei)
        end if
        
        if(present(rhs)) then
           call addto(rhs, face_global_nodes(p, sele), &
                -(matmul(mass_ele, face_val(p, sele)) &
                -matmul(mass_ele_old, face_val(prevp,sele)))*alpha)
        end if
        if (have_wd .and. present(rhs)) then
           call addto(rhs, face_global_nodes(p, sele), &
                +(matmul(mass_ele_wd, face_val(d0_pfield, sele)) &
                -matmul(mass_ele_old_wd, face_val(d0_pfield,sele)))*alpha*g)
        end if

      else
        ! no mesh movement - just use the same mass matrix as above
         if(present(rhs)) then
            call addto(rhs, face_global_nodes(p, sele), &
                 -1.0*matmul(mass_ele, face_val(p, sele)-face_val(prevp,sele))*alpha)
         end if
      end if
      
    end subroutine add_free_surface_element
    
  end subroutine add_free_surface_to_cmc_projection
    
  subroutine add_free_surface_to_poisson_rhs(poisson_rhs, state, dt, theta_pg)

    type(scalar_field), intent(inout) :: poisson_rhs
    type(state_type), intent(in) :: state
    real, intent(in) :: dt, theta_pg

    type(vector_field), pointer:: positions, u, gravity_normal
    type(scalar_field), pointer:: p
    character(len=FIELD_NAME_LEN):: bctype
    real g, coef, rho0
    integer, dimension(:), pointer:: surface_element_list
    integer i, j, grav_stat
    logical:: include_normals

    ewrite(1,*) 'Entering assemble_masslumped_poisson_rhs_free_surface'

    ! gravity acceleration
    call get_option('/physical_parameters/gravity/magnitude', g, stat=grav_stat)
      
    ! with a free surface the initial condition prescribed for pressure
    ! is used at the free surface nodes only
    p => extract_scalar_field(state, "Pressure")
    u => extract_vector_field(state, "Velocity")

    ! reference density
    call get_reference_density_from_options(state, u, rho0)

    ! only include the inner product of gravity and surface normal
    ! if the free surface nodes are actually moved (not necessary
    ! for large scale ocean simulations)
    include_normals = have_option("/mesh_adaptivity/mesh_movement/free_surface")
    if (include_normals) then
      ewrite(2,*) 'Including inner product of normals in kinematic bc'
      gravity_normal => extract_vector_field(state, "GravityDirection")
    end if
    
      
    ! adding in the free surface integral using the free surface
    ! elevation (p/g) specified by the inital pressure at the surface nodes
    positions => extract_vector_field(state, "Coordinate")
    coef=g*rho0*theta_pg**2*dt**2
      
    do i=1, get_boundary_condition_count(u)
      call get_boundary_condition(u, i, type=bctype, &
          surface_element_list=surface_element_list)
      if (bctype=="free_surface") then
        if (grav_stat/=0) then
           FLExit("For a free surface you need gravity")
        end if
        do j=1, size(surface_element_list)
          call add_free_surface_element(surface_element_list(j))
        end do
      end if
    end do
      
    contains
    
    subroutine add_free_surface_element(sele)
    integer, intent(in):: sele
      
      real, dimension(positions%dim, face_ngi(positions, sele)):: normals
      real, dimension(face_loc(p, sele), face_loc(p, sele)):: mass_ele
      real, dimension(face_ngi(p, sele)):: detwei
      integer:: ele
      
      ele=face_ele(positions, sele)
      call transform_facet_to_physical(positions, sele, detwei_f=detwei,&
           & normal=normals)
      if (include_normals) then
         ! at each gauss point multiply with inner product of gravity and surface normal
         detwei=detwei*(-1.0)*sum(face_val_at_quad(gravity_normal,sele)*normals, dim=1)
      end if
      mass_ele=shape_shape(face_shape(p, sele), face_shape(p, sele), detwei)
      call addto(poisson_rhs, face_global_nodes(p, sele), &
        matmul(mass_ele, face_val(p, sele))/coef)
      
    end subroutine add_free_surface_element
    
  end subroutine add_free_surface_to_poisson_rhs
    
  subroutine copy_poisson_solution_to_interior(p_theta, p, old_p, u)
  type(scalar_field), intent(inout):: p_theta, p, old_p
  type(vector_field), intent(in):: u
    
    character(len=FIELD_NAME_LEN):: bctype
    integer, dimension(:), pointer:: surface_element_list
    integer:: i, j, sele
    
    ! first copy initial free surface elevations (p/g) at free surface nodes
    ! to p_theta
    do i=1, get_boundary_condition_count(u)
      call get_boundary_condition(u, i, type=bctype, &
          surface_element_list=surface_element_list)
      if (bctype=="free_surface") then
        do j=1, size(surface_element_list)
          sele=surface_element_list(j)
          call set(p_theta, face_global_nodes(p,sele), face_val(p, sele))
        end do
      end if
    end do
      
    ! then copy everything (including interior) back from p_theta to p
    call set(p, p_theta)
    
    ! p and old_p are the same (as we're in the first non-linear iteration)
    ! but they might be different fields (if #nonlinear iterations>1)
    call set(old_p, p_theta)
      
  end subroutine copy_poisson_solution_to_interior
  
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
    type(vector_field), target :: local_grid_u
    type(scalar_field), target:: linear_p
    character(len=FIELD_NAME_LEN):: bctype
    real g, dt, rho0, gravity_magnitude
    integer, dimension(:), allocatable:: face_nodes
    integer, dimension(:), pointer:: surface_element_list
    integer i, j, k, node, sele, stat

    ! some fields for when moving the entire mesh
    type(scalar_field), pointer :: topdis, bottomdis
    !type(scalar_field), save :: fracdis
    type(scalar_field) :: fracdis
    type(scalar_field) :: extrapolated_p
    
    logical :: l_initialise, have_wd
    
    ewrite(1,*) 'Entering move_free_surface_nodes'
    
    ! increase event counter, so position caching know the mesh has moved
    call IncrementEventCounter(EVENT_MESH_MOVEMENT)
    
    l_initialise = present_and_true(initialise)
    
    ! gravity acceleration
    call get_option('/physical_parameters/gravity/magnitude', g)
    call get_option('/timestepping/timestep', dt)
    
    positions => extract_vector_field(state, "Coordinate")
    original_positions => extract_vector_field(state, "OriginalCoordinate")
    iterated_positions => extract_vector_field(state, "IteratedCoordinate")
    old_positions => extract_vector_field(state, "OldCoordinate")

    gravity_normal => extract_vector_field(state, "GravityDirection")
    ! it's alright for gravity to be on a DG version of the CoordinateMesh:
    assert( face_loc(gravity_normal,1)==face_loc(positions,1) )

    u => extract_vector_field(state, "Velocity")
    
    ! eos/fluids/linear/subtract_out_hydr.level is options checked below
    ! so the ref. density should be present
    call get_reference_density_from_options(state, u, rho0)
    
    p => extract_scalar_field(state, "Pressure")
    if (.not. p%mesh==positions%mesh) then
      call allocate(linear_p, positions%mesh)
      call remap_field(p, linear_p, stat=stat)
      if(stat==REMAP_ERR_DISCONTINUOUS_CONTINUOUS) then
        ewrite(-1,*) "Just remapped from a discontinuous to a continuous field when using free_surface mesh movement."
        ewrite(-1,*) "This suggests the pressure is discontinuous, which isn't supported."
        FLExit("Discontinuous pressure not permitted.")
      else if(stat/=0 .and. stat/=REMAP_ERR_UNPERIODIC_PERIODIC .and. stat/=REMAP_ERR_HIGHER_LOWER_CONTINUOUS) then
        FLAbort("Something went wrong mapping pressure to the CoordinateMesh")
      end if
      ! we've allowed it to remap from periodic to unperiodic and from higher order to lower order
      p => linear_p
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
    
      topdis => extract_scalar_field(state, "DistanceToTop")
      bottomdis => extract_scalar_field(state, "DistanceToBottom")
      
      ! first we need to extrapolate the pressure down from the surface
      call allocate(extrapolated_p, p%mesh, "ExtrapolatedPressure")
    
      call get_boundary_condition(topdis, 1, &
        surface_element_list=surface_element_list)
        
      ! Vertically extrapolate pressure values at the free surface downwards
      ! (reuse projected horizontal top surface mesh cached under DistanceToTop)
      ! The use of Coordinate here (as opposed to IteratedCoordinate or OldCoordinate)
      ! doesn't affect anything as nodes only move in the vertical.
      call VerticalExtrapolation(p, extrapolated_p, positions, &
        gravity_normal, surface_element_list=surface_element_list, &
        surface_name="DistanceToTop")

      have_wd=have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying")
      if (have_wd) then
         call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
         do node=1, node_count(extrapolated_p)
            call set(extrapolated_p, node, &
                  max(node_val(extrapolated_p, node),-gravity_magnitude*node_val(get_wettingdrying_d0(state), node)))
          end do
      end if

      ! then we need to scale it by its fractional distance from the bottom
      call allocate(fracdis, topdis%mesh, "FractionalDistance")      
      call set(fracdis, topdis)
      call addto(fracdis, bottomdis)
      call invert(fracdis)
      call scale(fracdis, bottomdis)
      call scale(extrapolated_p, fracdis)
      call deallocate(fracdis)
    
      do node=1, node_count(positions)
        call set(iterated_positions, node, &
                  node_val(original_positions, node)- &
                  node_val(extrapolated_p, node)*node_val(gravity_normal, node)/g/rho0)
                  
        if(.not.l_initialise) call set(grid_u, node, &
                  (node_val(iterated_positions, node)-node_val(old_positions,node))/dt)
      end do
      
      call deallocate(extrapolated_p)
    
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
                  node_val(p, node)*node_val(gravity_normal, node)/g/rho0)
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
      do i = 1, grid_u%dim
        ewrite_minmax(grid_u%val(i)%ptr)
      end do      
    end if
    
    if (associated(p, linear_p)) then
       call deallocate(p)
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
    integer:: stat, i, face
    
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
      ewrite(2,*) "Number of owned surface nodes:", key_count(owned_surface_nodes)
      ewrite(2,*) "Number of columns in vertical prolongator:", size(vertical_prolongator,2)
      if (size(vertical_prolongator,2)>key_count(owned_surface_nodes)) then
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
  type(state_type), intent(in):: state
  type(scalar_field), intent(inout):: free_surface
    
     integer, dimension(:), pointer:: surface_element_list
     type(vector_field), pointer:: x, u, vertical_normal
     type(scalar_field), pointer:: p, topdis, d0, d0_pfield
     character(len=FIELD_NAME_LEN):: bctype
     real:: g, rho0
     integer:: i, j, sele, stat
     logical :: have_wd

     ! the prognostic free surface is calculated elsewhere (this is the
     ! separate free surface equation approach in the old code path)
     if (have_option(trim(free_surface%option_path)//'/prognostic')) return
     
     x => extract_vector_field(state, "Coordinate")
     p => extract_scalar_field(state, "Pressure")
     assert(free_surface%mesh==p%mesh)

     call get_option('/physical_parameters/gravity/magnitude', g)

     u => extract_vector_field(state, "Velocity")
     call get_reference_density_from_options(state, u, rho0)
          
     topdis => extract_scalar_field(state, "DistanceToTop", stat=stat)
     if (stat==0) then
       ! note we're not using the actual free_surface bc here, as 
       ! that may be specified in parts, or not cover the whole area
       call get_boundary_condition(topdis, 1, &
         surface_element_list=surface_element_list)
         
       vertical_normal => extract_vector_field(state, "GravityDirection")
     
       ! vertically extrapolate pressure values at the free surface downwards
       ! (reuse projected horizontal top surface mesh cached under DistanceToTop)
       call VerticalExtrapolation(p, free_surface, x, &
         vertical_normal, surface_element_list=surface_element_list, &
         surface_name="DistanceToTop")
       ! divide by rho0 g
       call scale(free_surface, 1/g/rho0)
       
     else
     
       ! if no vertical extrapolation is available, only copy
       ! the values at the free surface nodes and divide by rho0 g
       
       ! make sure other nodes are zeroed
       call zero(free_surface)    

       do i=1, get_boundary_condition_count(u)
          call get_boundary_condition(u, i, type=bctype, &
             surface_element_list=surface_element_list)
          if (bctype=="free_surface") then
        
             face_loop: do j=1, size(surface_element_list)

               sele=surface_element_list(j)
               
               call set(free_surface, &
                 face_global_nodes(free_surface, sele), &
                 face_val(p, sele)/rho0/g)
               
             end do face_loop
               
          end if
          
       end do      
     end if

     ! Do the wetting and drying corrections
     have_wd=have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying")
     if (have_wd) then
        d0=>get_wettingdrying_d0(state)
        allocate(d0_pfield)
        call allocate(d0_pfield, p%mesh, "d0_pfield")
        call remap_field(d0, d0_pfield)
        call scale(free_surface, -1.0)
        call bound(free_surface, upper_bound=d0_pfield)
        call scale(free_surface, -1.0)
        call deallocate(d0_pfield)
        deallocate(d0_pfield)
     end if


    
  end subroutine calculate_diagnostic_free_surface

  subroutine get_reference_density_from_options(state, u, rho0)

    type(state_type), intent(in) :: state
    type(vector_field), intent(in) :: u
    real, intent(out) :: rho0

    if (have_option('/material_phase::'//trim(state%name)//'/equation_of_state/fluids/linear/')) then
      if (have_option(trim(u%option_path)//"/prognostic/equation::Boussinesq")) then
        ! for Boussinesq the pressure is actually p/rho0 already so we don't want to divide by
        ! the reference density again
        rho0 = 1.0
      else
        ! for other equation types the pressure is really the pressure so we do need to
        ! divide by the real reference density
        call get_option('/material_phase::'//trim(state%name)// &
          '/equation_of_state/fluids/linear/reference_density', rho0)
      end if
    else if (have_option('/material_phase::'//trim(state%name)//'/equation_of_state/fluids/ocean_pade_approximation/')) then
      ! The reference density is hard-coded to be 1.0 for the Ocean Pade Approximation
      rho0=1.0
    else 
      ewrite(-1,*) "Unless using Boussinesq Velocity, you must specify a"
      ewrite(-1,*) "linear or pade equation of state for the free surface."
      FLExit("Error retrieving reference density from options.")
    endif

  end subroutine get_reference_density_from_options

 subroutine update_wettingdrying_alpha(state)
  !!< calculates and updates the alpha coefficients for wetting and drying.
  type(state_type), intent(in):: state
  type(scalar_field), pointer:: scalar_surface_field

  integer, dimension(:), pointer :: surface_element_list
  type(vector_field), pointer:: u
  type(scalar_field), pointer:: p, d0_field
  type(scalar_field), pointer:: d0_pfield
  character(len=FIELD_NAME_LEN):: bctype
  character(len=OPTION_PATH_LEN) fs_option_path
  real:: rho0, gravity_magnitude, d0offset
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
             call get_reference_density_from_options(state, u, rho0)
             d0_field => get_wettingdrying_d0(state)
             call get_option('/physical_parameters/gravity/magnitude', gravity_magnitude)
             call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/d0", d0offset)

             allocate(d0_pfield)
             call allocate(d0_pfield, p%mesh, "d0_pmesh")
             call remap_field(d0_field, d0_pfield)
             ! Calculate alpha for each surface element
             face_loop: do j=1, size(surface_element_list)
                sele=surface_element_list(j)
                call calculate_alpha(sele, gravity_magnitude, p, d0_pfield, d0offset, alpha )
                call set(scalar_surface_field, &
                     ele_nodes(scalar_surface_field, j), &
                     alpha)
             end do face_loop
             call deallocate(d0_pfield)
            deallocate(d0_pfield)
       end if
   end do
   deallocate(alpha)

  contains
   subroutine calculate_alpha(sele, gravity_magnitude, p, d0_pfield, d0offset, alpha)
        type(scalar_field), pointer, intent(in) :: p, d0_pfield
        integer, intent(in) :: sele
        real, intent(in) :: gravity_magnitude, d0offset
        real, dimension(face_loc(p, sele)), intent(inout) :: alpha

        !!! The calculation of alpha is based on pressure, thus pressure and free surface don't have to be coupled.
        alpha = gravity_magnitude * face_val(d0_pfield, sele) + face_val(p, sele)
        alpha = alpha/gravity_magnitude * (-d0offset)
        alpha = max(0.0, alpha)
    end subroutine calculate_alpha
 end subroutine update_wettingdrying_alpha

 ! Returns a scalar field on the CoordinateMesh with the auxiliary wetting and drying d0 variable, saving the distance between
 ! steady state free surface height and bathymetry + d0 offset.
 ! In non-wetting and drying regions the value of d0_pfield is INFINITY.
 ! This field multiplied by -1.0*gravity_magnitude gives the free surface pressure on the d0 level.
 function get_wettingdrying_d0(state) result(d0)
    type(state_type), intent(in):: state

  logical,save  :: initialized=.false.
  real :: d0offset
  real, dimension(:), allocatable :: d0offset_loc
  type(scalar_field), pointer :: bottomdist, topdis
  type(vector_field), pointer :: x, vertical_normal
  type(scalar_field), save, target :: d0_field
  type(scalar_field), pointer :: d0
  type(vector_field), pointer:: u
  character(len=FIELD_NAME_LEN):: bctype
  integer:: i, j, sele, stat
  integer, dimension(:), pointer:: surface_element_list

  if (.not. initialized) then
     bottomdist => extract_scalar_field(state, "DistanceToBottom")
     u => extract_vector_field(state, "Velocity")
     allocate(d0offset_loc(face_loc(bottomdist, 1)))
     call allocate(d0_field, bottomdist%mesh, "d0")
     call set(d0_field, INFINITY)
     do i = 1, get_boundary_condition_count(u)
            call get_boundary_condition(u, i, type=bctype, &
            surface_element_list = surface_element_list)
            if (bctype=="free_surface") then
                call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/d0", d0offset)
                d0offset_loc=d0offset
                do j = 1, size(surface_element_list)
                    sele = surface_element_list(j)
                    call set(d0_field, &
                    face_global_nodes(d0_field, sele), &
                    -d0offset_loc)
                end do
            end if
     end do
     call addto(d0_field, bottomdist)

     ! Vertical Extrapolation since the values are only valid on the free surface for now.
     topdis => extract_scalar_field(state, "DistanceToTop", stat=stat)
     assert(stat==0)
     if (stat==0) then
       call get_boundary_condition(topdis, 1, &
         surface_element_list=surface_element_list)

       vertical_normal => extract_vector_field(state, "GravityDirection")
       x => extract_vector_field(state, "Coordinate")

       ! vertically extrapolate pressure values at the free surface downwards
       ! (reuse projected horizontal top surface mesh cached under DistanceToTop)
       call VerticalExtrapolation(d0_field, d0_field, x, &
         vertical_normal, surface_element_list=surface_element_list)
     end if

     initialized=.true.
     deallocate(d0offset_loc)
  end if
  d0=>d0_field
  
end function get_wettingdrying_d0


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


  subroutine free_surface_module_check_options
    
    character(len=OPTION_PATH_LEN):: option_path, phase_path, pressure_path, pade_path
    character(len=FIELD_NAME_LEN):: fs_meshname, p_meshname, bctype
    logical:: have_free_surface
    integer i, p
    
    do p=1, option_count('/material_phase')
      phase_path='/material_phase['//int2str(p-1)//']'
      pressure_path=trim(phase_path)//'/scalar_field::Pressure/prognostic'
      
      ! check if we have a free_surface bc
      option_path=trim(phase_path)//'/vector_field::Velocity/prognostic'
      if (have_option(trim(option_path))) then
        have_free_surface=.false.
        do i=1, option_count(trim(option_path)//'/boundary_conditions')
          call get_option(trim(option_path)//'/boundary_conditions['// &
             int2str(i-1)//']/type[0]/name', bctype)
          have_free_surface=have_free_surface .or. (bctype=='free_surface')
        end do
      else
        ! no prognostic velocity, no free_surface bc
        have_free_surface=.false.
      end if
      
      if (have_free_surface) then
         ewrite(2,*) "You have a free surface, checking its options"
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
        if (.not. have_free_surface) then
          ewrite(-1,*) "The diagnostic FreeSurface field has to be used in combination " // &
            "with the free_surface boundary condition under Velocity."
          FLExit("Exit")
        end if
        if (.not. fs_meshname==p_meshname) then
          FLExit("The diagnostic FreeSurface field and the Pressure field have to be on the same mesh")
        end if
        if (.not. have_option('/geometry/ocean_boundaries')) then
          ewrite(0,*) "Warning: your diagnostic free surface will only be " // &
            "defined at the free surface nodes and not extrapolated downwards, " // &
            "because you didn't specify geometry/ocean_boundaries."
        end if
      end if
      
      ! check we're not combining old and new free surface method
      option_path=trim(phase_path)//'/scalar_field::FreeSurface/prognostic'
      if (have_free_surface .and. have_option(trim(option_path))) then
        ewrite(-1,*) "Trying to combine free_surface boundary condition (new method) " // &
          "with prognostic FreeSurface (old method)."
        FLExit("Cannot use both old and new free surface method")
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
    end do
    
  end subroutine free_surface_module_check_options

end module free_surface_module

