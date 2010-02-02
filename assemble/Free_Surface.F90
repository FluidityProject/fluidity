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
use sparse_matrices_fields
use boundary_conditions
use spud
use vertical_extrapolation_module
use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
use parallel_tools
use halos
use eventcounter
implicit none

private

public add_free_surface_to_cmc_projection, &
  move_free_surface_nodes, vertical_prolongator_from_free_surface, &
  free_surface_nodes, calculate_diagnostic_free_surface, &
  add_free_surface_to_poisson_rhs, copy_poisson_solution_to_interior
  

public free_surface_module_check_options

contains

  subroutine add_free_surface_to_cmc_projection(state, cmc, rhs, dt, theta, get_cmc)
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
    type(scalar_field), intent(inout) :: rhs
    real, intent(in) :: dt, theta
    !! only add in to the matrix if get_cmc==.true.
    logical, intent(in):: get_cmc
    
      type(vector_field), pointer:: positions, u, gravity_normal
      type(scalar_field), pointer:: p, prevp
      character(len=FIELD_NAME_LEN):: bctype
      real:: g, rho0
      integer, dimension(:), pointer:: surface_element_list
      integer:: i, j, grav_stat
      logical:: include_normals
      
      ewrite(1,*) 'Entering add_free_surface_to_cmc_projection'
      
      ! gravity acceleration
      call get_option('/physical_parameters/gravity/magnitude', g, stat=grav_stat)
      ! reference density
      call get_option('/material_phase::'//trim(state%name)// &
        '/equation_of_state/fluids/linear/reference_density', rho0)
      
      ! get the pressure, and the pressure at the beginning of the time step
      p => extract_scalar_field(state, "Pressure")
      prevp => extract_scalar_field(state, "OldPressure")
      u => extract_vector_field(state, "Velocity")
      
      ! only include the inner product of gravity and surface normal
      ! if the free surface nodes are actually moved (not necessary
      ! for large scale ocean simulations)
      include_normals = have_option("/mesh_adaptivity/mesh_movement/free_surface")
      if (include_normals) then
        ewrite(2,*) 'Including inner product of normals in kinematic bc'
        gravity_normal => extract_vector_field(state, "GravityDirection")
      end if
      
      ewrite_minmax(p)
      ewrite_minmax(prevp)
      ewrite_minmax(rhs)
      
      positions => extract_vector_field(state, "Coordinate")
      
      do i=1, get_boundary_condition_count(u)
        call get_boundary_condition(u, i, type=bctype, &
           surface_element_list=surface_element_list)
        if (bctype=="free_surface") then
          if (grav_stat/=0) then
             FLAbort("For a free surface you need gravity")
          end if
          do j=1, size(surface_element_list)
            call add_free_surface_element(surface_element_list(j))
          end do
        end if
      end do
    
      ewrite_minmax(rhs)
      
    contains 
    
    subroutine add_free_surface_element(sele)
    integer, intent(in):: sele
      
      real, dimension(positions%dim, face_ngi(positions, sele)):: normals
      real, dimension(face_loc(p, sele), face_loc(p, sele)):: mass_ele
      real, dimension(face_ngi(p, sele)):: detwei
      real:: alpha
      integer:: ele
      
      ele=face_ele(positions, sele)
      call transform_facet_to_physical(positions, sele, detwei_f=detwei,&
           & normal=normals)
      if (include_normals) then
        ! at each gauss point multiply with inner product of gravity and surface normal
        detwei=detwei*(-1.0)*sum(face_val_at_quad(gravity_normal,sele)*normals, dim=1)
      end if
      
      mass_ele=shape_shape(face_shape(p, sele), face_shape(p, sele), detwei)
      alpha=1.0/g/rho0/dt
      if (get_cmc) then
        ! we consider the projection equation to solve for 
        ! phi=theta^2 dt dp, so that the f.s. integral
        ! alpha M_fs dp=alpha M_fs phi/(theta^2 dt)
        call addto(cmc, &
         face_global_nodes(p, sele), face_global_nodes(p,sele), &
         alpha*mass_ele/(theta**2*dt) )
      end if
      call addto(rhs, face_global_nodes(p, sele), &
        -1.0*matmul(mass_ele, face_val(p, sele)-face_val(prevp,sele))*alpha)
      
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
    ! reference density
    call get_option('/material_phase::'//trim(state%name)// &
        '/equation_of_state/fluids/linear/reference_density', rho0)
      
    ! with a free surface the initial condition prescribed for pressure
    ! is used at the free surface nodes only
    p => extract_scalar_field(state, "Pressure")
    u => extract_vector_field(state, "Velocity")
    
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
           FLAbort("For a free surface you need gravity")
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
      
  subroutine move_free_surface_nodes(state, theta)
  type(state_type), intent(inout) :: state
  real, intent(in), optional :: theta
    
    type(vector_field), pointer:: positions, u, original_positions
    type(vector_field), pointer:: gravity_normal, old_positions, grid_u
    type(vector_field), pointer:: iterated_positions
    type(scalar_field), pointer:: p
    type(vector_field) :: local_grid_u, old_grid_u, iterated_grid_u
    type(scalar_field), target:: linear_p
    character(len=FIELD_NAME_LEN):: bctype
    real g, dt, rho0
    integer, dimension(:), allocatable:: face_nodes
    integer, dimension(:), pointer:: surface_element_list
    integer i, j, k, node, sele

    ! some fields for when moving the entire mesh
    type(scalar_field), pointer :: topdis, bottomdis
    type(scalar_field) :: fracdis
    type(scalar_field) :: extrapolated_p
    
    ewrite(1,*) 'Entering move_free_surface_nodes'
    
    ! increase event counter, so position caching know the mesh has moved
    call IncrementEventCounter(EVENT_MESH_MOVEMENT)
    
    ! gravity acceleration
    call get_option('/physical_parameters/gravity/magnitude', g)
    call get_option('/timestepping/timestep', dt)
    ! eos/fluids/linear/subtract_out_hydr.level is options checked below
    ! so the ref. density should be present
    call get_option('/material_phase::'//trim(state%name)// &
      '/equation_of_state/fluids/linear/reference_density', rho0)
    
    positions => extract_vector_field(state, "Coordinate")
    original_positions => extract_vector_field(state, "OriginalCoordinate")
    iterated_positions => extract_vector_field(state, "IteratedCoordinate")
    old_positions => extract_vector_field(state, "OldCoordinate")
    gravity_normal => extract_vector_field(state, "GravityDirection")
    
      
    p => extract_scalar_field(state, "Pressure")
    if (.not. p%mesh==positions%mesh) then
      call allocate(linear_p, positions%mesh)
      call remap_field(p, linear_p)
      p => linear_p
    end if
    ! it's alright for gravity to be on a DG version of the CoordinateMesh:
    assert( face_loc(gravity_normal,1)==face_loc(positions,1) )
    u => extract_vector_field(state, "Velocity")
    grid_u => extract_vector_field(state, "GridVelocity")
    ! allocate this on the grid_u mesh to save its values
    call allocate(old_grid_u, grid_u%dim, grid_u%mesh, "TempOldGridVelocity")
    call set(old_grid_u, grid_u)
    ! allocate this on the grid_u mesh to remap to
    call allocate(iterated_grid_u, grid_u%dim, grid_u%mesh, "TempIteratedGridVelocity")
    
    ! allocate this on the positions mesh to calculate the values
    call allocate(local_grid_u, grid_u%dim, positions%mesh, "LocalGridVelocity")
    call zero(local_grid_u)
    
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
                  
        call set(local_grid_u, node, &
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
                call set(local_grid_u, node, &
                  (node_val(iterated_positions, node)-node_val(old_positions,node))/dt)
            end do node_loop
            
          end do face_loop
            
        end if
      end do
      
    end if
    
    call remap_field(local_grid_u, iterated_grid_u)

    if(present(theta)) then
      call set(positions, iterated_positions, old_positions, theta)
      call set(grid_u, iterated_grid_u, old_grid_u, theta)
    else
      call set(positions, iterated_positions)
      call set(grid_u, iterated_grid_u)
    end if
    
    do i = 1, grid_u%dim
      ewrite_minmax(grid_u%val(i)%ptr)
    end do
    
    if (associated(p, linear_p)) then
       call deallocate(p)
    end if
    call deallocate(local_grid_u)
    call deallocate(iterated_grid_u)
    call deallocate(old_grid_u)
          
  end subroutine move_free_surface_nodes
  
  function vertical_prolongator_from_free_surface(state, mesh) result (vertical_prolongator)
  !! Creates a prolongation operator from the free surface mesh to
  !! the full mesh which when provided to petsc_solve is used
  !! to implement vertical lumping if the "mg" preconditioner is selected.
  type(state_type), intent(in):: state
  type(mesh_type), intent(in):: mesh
  type(csr_matrix):: vertical_prolongator
    
    type(mesh_type), pointer:: coordinate_surface_mesh
    type(mesh_type):: surface_mesh
    type(scalar_field), pointer:: topdis
    type(vector_field), pointer:: positions
    type(vector_field) mesh_positions, surface_positions
    integer, dimension(:), pointer:: surface_element_list, surface_node_list
    integer stat, i, nowned_surface_nodes, node
    
    ewrite(1, *) "Constructing vertical_prolongator_from_free_surface to be used in mg"
    topdis => extract_scalar_field(state, "DistanceToTop", stat=stat)
    if (stat/=0) then
       FLExit("For vertical lumping you need to specify the ocean_boundaries under /geometry")
    end if
       
    positions => extract_vector_field(state, "Coordinate")
    if (.not. (positions%mesh==mesh)) then
      call allocate(mesh_positions, positions%dim, mesh)
      call remap_field(positions, mesh_positions)
    else
      mesh_positions=positions
    end if

    call get_boundary_condition(topdis, 1, &
      surface_mesh=coordinate_surface_mesh, &
      surface_element_list=surface_element_list, &
      surface_node_list=surface_node_list)
    if (coordinate_surface_mesh%shape%degree/=mesh%shape%degree) then
      ! note we replace surface_node_list here:
      call create_surface_mesh(surface_mesh, surface_node_list, &
        mesh, surface_element_list, name="SurfaceMesh_VerticalLumping")
    else
      surface_mesh=coordinate_surface_mesh
      ! note we keep surface_node_list associated with the coordinate mesh
    end if
    call allocate(surface_positions, positions%dim, surface_mesh, &
       name="FreeSurfacePositions")
    ! this does a remap if coordinate shape/=mesh shape degree
    call remap_field_to_surface(positions, surface_positions, surface_element_list)
    
    if (IsParallel()) then
      assert( associated(mesh%halos) )
      vertical_prolongator=VerticalProlongationOperator( &
         mesh_positions, surface_positions, reduce_columns=.true., &
         owned_nodes=halo_nowned_nodes(mesh%halos(1)) )
      ! note that in surface_positions the non-owned free surface nodes may be inbetween
      ! the reduce_columns option should have removed those however
      ! with debugging perform test to check if this is the case:
#ifdef DDEBUG
      ! count n/o owned surface nodes
      nowned_surface_nodes=0
      do i=1, node_count(surface_mesh)
        node=surface_node_list(i)
        if (node_owned(mesh, node)) then
          nowned_surface_nodes=nowned_surface_nodes+1
        end if
      end do
      ewrite(2,*) "Number of owned surface nodes:", nowned_surface_nodes
      ewrite(2,*) "Number of columns in vertical prolongator:", size(vertical_prolongator,2)
      if (size(vertical_prolongator,2)>nowned_surface_nodes) then
        ewrite(-1,*) "Vertical prolongator seems to be using more surface nodes than the number"
        ewrite(-1,*) "of surface nodes within completely owned surface elements. This indicates"
        ewrite(-1,*) "the parallel decomposition is not done along columns. You shouldn't be using"
        ewrite(-1,*) "mg with vertical_lumping in that case."
        FLAbort("Vertical lumping requires 2d decomposition along columns")
      end if
#endif
    else
      vertical_prolongator=VerticalProlongationOperator( &
         mesh_positions, surface_positions, reduce_columns=.true.)
    end if

    if (.not. (positions%mesh==mesh)) then
      call deallocate(mesh_positions)
    end if
    if (coordinate_surface_mesh%shape%degree/=mesh%shape%degree) then
      call deallocate(surface_mesh)
    end if
    call deallocate(surface_positions)

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
       FLAbort("Need to specify the ocean_boundaries under /geometry")
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
     type(scalar_field), pointer:: p, topdis
     character(len=FIELD_NAME_LEN):: bctype
     real:: g, rho0
     integer:: i, j, sele, stat
     
     ! the prognostic free surface is calculated elsewhere (this is the
     ! separate free surface equation approach in the old code path)
     if (have_option(trim(free_surface%option_path)//'/prognostic')) return
     
     x => extract_vector_field(state, "Coordinate")
     p => extract_scalar_field(state, "Pressure")
     assert(free_surface%mesh==p%mesh)

     call get_option('/physical_parameters/gravity/magnitude', g)
     call get_option('/material_phase::'//trim(state%name)// &
      '/equation_of_state/fluids/linear/reference_density', rho0)
          
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
       
       u => extract_vector_field(state, "Velocity")
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
    
  end subroutine calculate_diagnostic_free_surface
    
  subroutine free_surface_module_check_options
    
    character(len=OPTION_PATH_LEN):: option_path, phase_path, pressure_path
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
         FLExit("you have to use continuous_galerkin or discontinuous_galerkin")
      end if

      ! check pressure options
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
          ewrite(0,*) "Warning: your diagnostic free surface will only be" // &
            "defined at the free surface nodes and not extrapolated downwards," // &
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
      if (have_free_surface .and. .not.(have_option(option_path))) then
        ewrite(-1,*) "Missing option: ", trim(option_path)
        FLExit("With the free surface you need to subtract out the hydrostatic level.")
      end if
    end do
    
  end subroutine free_surface_module_check_options

end module free_surface_module

