!     Copyright (C) 2006 Imperial College London and others.
!     
!     Please see the AUTHORS file in the main source directory for a full list
!     of copyright holders.
!     
!     Prof. C Pain
!     Applied Modelling and Computation Group
!     Department of Earth Science and Engineering
!     Imperial College London
!     
!     C.Pain@Imperial.ac.uk
!     
!     This library is free software; you can redistribute it and/or
!     modify it under the terms of the GNU Lesser General Public
!     License as published by the Free Software Foundation,
!     version 2.1 of the License.
!     
!     This library is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!     Lesser General Public License for more details.
!     
!     You should have received a copy of the GNU Lesser General Public
!     License along with this library; if not, write to the Free Software
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!     USA

#include "fdebug.h"

module advection_diffusion_cg
  
  ! keep in this order, please:
  use quadrature
  use elements
  use sparse_tools
  use fields
  !
  use boundary_conditions
  use boundary_conditions_from_options
  use field_options
  use fldebug
  use global_parameters, only : FIELD_NAME_LEN, OPTION_PATH_LEN
  use spud
  use solvers
  use state_module
  use upwind_stabilisation
  use sparsity_patterns_meshes
  
  implicit none
  
  private
  
  public :: solve_field_equation_cg, advection_diffusion_cg_check_options
  
  character(len = *), parameter, public :: advdif_cg_m_name = "AdvectionDiffusionCGMatrix"
  character(len = *), parameter, public :: advdif_cg_rhs_name = "AdvectionDiffusionCGRHS"
  character(len = *), parameter, public :: advdif_cg_delta_t_name = "AdvectionDiffusionCGChange"
  character(len = *), parameter, public :: advdif_cg_velocity_name = "AdvectionDiffusionCGVelocity"
  
  ! Stabilisation schemes
  integer, parameter :: STABILISATION_NONE = 0, &
    & STABILISATION_STREAMLINE_UPWIND = 1, STABILISATION_SUPG = 2
  
  ! Boundary condition types
  integer, parameter :: BC_TYPE_NEUMANN = 1, BC_TYPE_WEAKDIRICHLET = 2
  
  ! Global variables, set by assemble_advection_diffusion_cg for use by
  ! assemble_advection_diffusion_element_cg and
  ! assemble_advection_diffusion_face_cg
  
  ! Local timestep
  real :: local_dt
  ! Implicitness/explicitness factor * timestep
  real :: dt_theta
  ! Implicitness/explicitness factor
  real :: theta
  ! Conservative/non-conservative discretisation factor
  real :: beta
  ! Stabilisation scheme
  integer :: stabilisation_scheme
  integer :: nu_bar_scheme
  real :: nu_bar_scale
  
  ! equation type
  integer :: equation_type
  ! Implicitness/explicitness factor for density
  real :: density_theta
  ! Which terms do we have?
  
  ! Mass term?
  logical :: have_mass
  ! Lump mass?
  logical :: lump_mass
  ! Advection?
  logical :: have_advection
  ! Integrate advection by parts?
  logical :: integrate_advection_by_parts
  ! Source?
  logical :: have_source
  ! Absorption?
  logical :: have_absorption
  ! Diffusivity?
  logical :: have_diffusivity
  ! Isotropic diffusivity?
  logical :: isotropic_diffusivity
  ! Is the mesh moving?
  logical :: move_mesh

contains

  subroutine solve_field_equation_cg(field_name, state, dt, velocity_name)
    !!< Construct and solve the advection-diffusion equation for the given
    !!< field using a continuous Galerkin discretisation. Based on
    !!< Advection_Diffusion_DG and Momentum_CG.
    
    character(len = *), intent(in) :: field_name
    type(state_type), intent(inout) :: state
    real, intent(in) :: dt
    character(len = *), optional, intent(in) :: velocity_name
    
    type(csr_matrix) :: matrix
    type(scalar_field) :: delta_t, rhs
    type(scalar_field), pointer :: t
    
    ewrite(1, *) "In solve_field_equation_cg"
    
    ewrite(2, *) "Solving advection-diffusion equation for field " // &
      & trim(field_name) // " in state " // trim(state%name)

    call initialise_advection_diffusion_cg(field_name, t, delta_t, matrix, rhs, state)

    call assemble_advection_diffusion_cg(t, matrix, rhs, state, dt, velocity_name = velocity_name)

    call solve_advection_diffusion_cg(t, delta_t, matrix, rhs)

    call apply_advection_diffusion_cg_change(t, delta_t, dt)
    
    call finalise_advection_diffusion_cg(delta_t, matrix, rhs)
    
    ewrite(1, *) "Exiting solve_field_equation_cg"
    
  end subroutine solve_field_equation_cg
  
  subroutine initialise_advection_diffusion_cg(field_name, t, delta_t, matrix, rhs, state)
    character(len = *), intent(in) :: field_name
    type(scalar_field), pointer :: t
    type(scalar_field), intent(out) :: delta_t
    type(csr_matrix), intent(out) :: matrix
    type(scalar_field), intent(out) :: rhs
    type(state_type), intent(inout) :: state
    
    integer :: stat
    type(csr_sparsity), pointer :: sparsity
    type(scalar_field), pointer :: t_old
        
    t => extract_scalar_field(state, field_name)
    if(t%mesh%continuity /= 0) then
      FLAbort("CG advection-diffusion requires a continuous mesh")
    end if
        
    t_old => extract_scalar_field(state, "Old" // field_name, stat = stat)
    if(stat == 0) then
      assert(t_old%mesh == t%mesh)
      ! Reset t to value at the beginning of the timestep
      call set(t, t_old)
    end if
    
    sparsity => get_csr_sparsity_firstorder(state, t%mesh, t%mesh)
    
    call allocate(matrix, sparsity, name = advdif_cg_m_name)
    call allocate(rhs, t%mesh, name = advdif_cg_rhs_name)
    call allocate(delta_t, t%mesh, name = advdif_cg_delta_t_name)
    
    call set_advection_diffusion_cg_initial_guess(delta_t)
    
  end subroutine initialise_advection_diffusion_cg
  
  subroutine finalise_advection_diffusion_cg(delta_t, matrix, rhs)
    type(scalar_field), intent(inout) :: delta_t
    type(csr_matrix), intent(inout) :: matrix
    type(scalar_field), intent(inout) :: rhs
    
    call deallocate(matrix)
    call deallocate(rhs)
    call deallocate(delta_t)
    
  end subroutine finalise_advection_diffusion_cg
  
  subroutine set_advection_diffusion_cg_initial_guess(delta_t)
    type(scalar_field), intent(inout) :: delta_t
    
    call zero(delta_t)
    
  end subroutine set_advection_diffusion_cg_initial_guess
  
  subroutine assemble_advection_diffusion_cg(t, matrix, rhs, state, dt, velocity_name)
    type(scalar_field), intent(inout) :: t
    type(csr_matrix), intent(inout) :: matrix
    type(scalar_field), intent(inout) :: rhs
    type(state_type), intent(in) :: state
    real, intent(in) :: dt
    character(len = *), optional, intent(in) :: velocity_name
    
    character(len = FIELD_NAME_LEN) :: lvelocity_name
    integer :: i, j, stat
    integer, dimension(:), allocatable :: t_bc_types
    type(scalar_field) :: t_bc
    type(scalar_field), pointer :: absorption, sinking_velocity, source
    type(tensor_field), pointer :: diffusivity
    type(vector_field) :: velocity
    type(vector_field), pointer :: gravity_direction, velocity_ptr, grid_velocity
    type(vector_field), pointer :: positions, old_positions, new_positions
    type(scalar_field), target :: dummydensity
    type(scalar_field), pointer :: density, olddensity
    character(len = FIELD_NAME_LEN) :: density_name
    type(scalar_field), pointer :: pressure
        
    ewrite(1, *) "In assemble_advection_diffusion_cg"
    
    assert(mesh_dim(rhs) == mesh_dim(t))
    assert(ele_count(rhs) == ele_count(t))
    
    if(present(velocity_name)) then
      lvelocity_name = velocity_name
    else
      lvelocity_name = "NonlinearVelocity"
    end if
    
    ! Step 1: Pull fields out of state
    
    ! Coordinate
    positions => extract_vector_field(state, "Coordinate")
    do i = 1, positions%dim
      ewrite_minmax(positions%val(i)%ptr)
    end do
    assert(positions%dim == mesh_dim(t))
    assert(ele_count(positions) == ele_count(t))
    
    ! Velocity    
    velocity_ptr => extract_vector_field(state, lvelocity_name, stat = stat)
    if(stat == 0) then
      assert(velocity_ptr%dim == mesh_dim(t))
      assert(ele_count(velocity_ptr) == ele_count(t))
      
      ewrite(2, *) "Velocity:"
      do i = 1, velocity_ptr%dim
        ewrite_minmax(velocity_ptr%val(i)%ptr)
      end do
      
      call allocate(velocity, velocity_ptr%dim, velocity_ptr%mesh, name = advdif_cg_velocity_name)
      call set(velocity, velocity_ptr)
    else
      ewrite(2, *) "No velocity"
      call allocate(velocity, mesh_dim(t), t%mesh, name = advdif_cg_velocity_name)
      call zero(velocity)
    end if
        
    ! Source
    source => extract_scalar_field(state, trim(t%name) // "Source", stat = stat)
    have_source = stat == 0
    if(have_source) then
      assert(mesh_dim(source) == mesh_dim(t))
      assert(ele_count(source) == ele_count(t))
    
      ewrite_minmax(source%val)
    else
      ewrite(2, *) "No source"
    end if
    
    ! Absorption
    absorption => extract_scalar_field(state, trim(t%name) // "Absorption", stat = stat)
    have_absorption = stat == 0
    if(have_absorption) then
      assert(mesh_dim(absorption) == mesh_dim(t))
      assert(ele_count(absorption) == ele_count(t))
    
      ewrite_minmax(absorption%val)
    else
      ewrite(2, *) "No absorption"
    end if
    
    ! Sinking velocity
    sinking_velocity => extract_scalar_field(state, trim(t%name) // "SinkingVelocity", stat = stat)
    if(stat == 0) then
      ewrite_minmax(sinking_velocity%val)
      
      gravity_direction => extract_vector_field(state, "GravityDirection")
      ! this may perform a "remap" internally from CoordinateMesh to VelocitMesh
      call addto(velocity, gravity_direction, scale = sinking_velocity)
    else
      ewrite(2, *) "No sinking velocity"
    end if
    
    ! Diffusivity
    diffusivity => extract_tensor_field(state, trim(t%name) // "Diffusivity", stat = stat)
    have_diffusivity = stat == 0
    if(have_diffusivity) then
      assert(diffusivity%dim == mesh_dim(t))
      assert(ele_count(diffusivity) == ele_count(t))
      
      isotropic_diffusivity = option_count(complete_field_path(diffusivity%option_path)) &
        & == option_count(trim(complete_field_path(diffusivity%option_path)) // "/value/isotropic")
        
      if(isotropic_diffusivity) then
        ewrite(2, *) "Isotropic diffusivity"
        assert(diffusivity%dim > 0)
        ewrite_minmax(diffusivity%val(1, 1, :))
      else
        do i = 1, diffusivity%dim
          do j = 1, diffusivity%dim
            ewrite_minmax(diffusivity%val(i, j, :))
          end do
        end do
      end if
    else
      isotropic_diffusivity = .false.
      ewrite(2, *) "No diffusivity"
    end if
    
    ! Step 2: Pull options out of the options tree
    
    call get_option(trim(t%option_path) // "/prognostic/temporal_discretisation/theta", theta)
    assert(theta >= 0.0 .and. theta <= 1.0)
    ewrite(2, *) "Theta = ", theta
    dt_theta = dt * theta
    local_dt = dt
    
    call get_option(trim(t%option_path) // "/prognostic/spatial_discretisation/conservative_advection", beta)
    assert(beta >= 0.0 .and. beta <= 1.0)
    ewrite(2, *) "Beta = ", beta
    
    have_advection = .not. have_option(trim(t%option_path) // "/prognostic/spatial_discretisation/continuous_galerkin/advection_terms/exclude_advection_terms")
    if(have_advection) then
      ewrite(2, *) "Including advection"
      
      integrate_advection_by_parts = have_option(trim(t%option_path) // "/prognostic/spatial_discretisation/continuous_galerkin/advection_terms/integrate_advection_by_parts")
      if(integrate_advection_by_parts) then
        ewrite(2, *) "Integrating advection terms by parts"
      end if
    else
      integrate_advection_by_parts = .false.
      ewrite(2, *) "Excluding advection"
    end if
    
    have_mass = .not. have_option(trim(t%option_path) // "/prognostic/spatial_discretisation/continuous_galerkin/mass_terms/exclude_mass_terms")
    if(have_mass) then
      ewrite(2, *) "Including mass"
      
      lump_mass = have_option(trim(t%option_path) // "/prognostic/spatial_discretisation/continuous_galerkin/mass_terms/lump_mass_matrix")
      if(lump_mass) then
        ewrite(2, *) "Lumping mass"
      end if
    else
      lump_mass = .false.
      ewrite(2, *) "Excluding mass"
    end if
    
    if(have_option(trim(t%option_path) // "/prognostic/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind")) then
      ewrite(2, *) "Streamline upwind stabilisation"
      stabilisation_scheme = STABILISATION_STREAMLINE_UPWIND
      call get_upwind_options(trim(t%option_path) // "/prognostic/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind", &
          & nu_bar_scheme, nu_bar_scale)
    else if(have_option(trim(t%option_path) // "/prognostic/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind_petrov_galerkin")) then
      ewrite(2, *) "SUPG stabilisation"
      stabilisation_scheme = STABILISATION_SUPG
      call get_upwind_options(trim(t%option_path) // "/prognostic/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind_petrov_galerkin", &
          & nu_bar_scheme, nu_bar_scale)
    else
      ewrite(2, *) "No stabilisation"
      stabilisation_scheme = STABILISATION_NONE
    end if
    
    ! are we moving the mesh?
    move_mesh = (have_option("/mesh_adaptivity/mesh_movement") .and. have_mass)
    if(move_mesh) then
      ewrite(2,*) "Moving the mesh"
      old_positions => extract_vector_field(state, "OldCoordinate")
      do i = 1, old_positions%dim
        ewrite_minmax(old_positions%val(i)%ptr)
      end do
      new_positions => extract_vector_field(state, "IteratedCoordinate")
      do i = 1, new_positions%dim
        ewrite_minmax(new_positions%val(i)%ptr)
      end do
      
      ! Grid velocity
      grid_velocity => extract_vector_field(state, "GridVelocity")
      assert(grid_velocity%dim == mesh_dim(t))
      assert(ele_count(grid_velocity) == ele_count(t))
      
      ewrite(2, *) "Grid velocity:"    
      do i = 1, grid_velocity%dim
        ewrite_minmax(grid_velocity%val(i)%ptr)
      end do
    else
      ewrite(2,*) "Not moving the mesh"
    end if

    
    call allocate(dummydensity, t%mesh, "DummyDensity", field_type=FIELD_TYPE_CONSTANT)
    call set(dummydensity, 1.0)
    ! find out equation type and hence if density is needed or not
    equation_type=equation_type_index(trim(t%option_path))
    select case(equation_type)
    case(FIELD_EQUATION_ADVECTIONDIFFUSION)
      ewrite(2,*) "Solving advection-diffusion equation"
      ! density not needed so use a constant field for assembly
      density => dummydensity
      olddensity => dummydensity
      density_theta = 1.0
      pressure => dummydensity
    case(FIELD_EQUATION_INTERNALENERGY)
      ewrite(2,*) "Solving internal energy equation"
      call get_option(trim(t%option_path)//'/prognostic/equation[0]/density[0]/name', &
                      density_name)
      density=>extract_scalar_field(state, trim(density_name))
      ewrite_minmax(density%val)
      
      olddensity=>extract_scalar_field(state, "Old"//trim(density_name))
      ewrite_minmax(olddensity%val)
      
      call get_option(trim(density%option_path)//"/prognostic/temporal_discretisation/theta", &
                      density_theta)
                      
      pressure=>extract_scalar_field(state, "Pressure")
      ewrite_minmax(pressure%val)
    case default
      FLExit("Unknown field equation type for cg advection diffusion.")
    end select
    
    ! Step 3: Assembly
    
    call zero(matrix)
    call zero(rhs)
    
    do i = 1, ele_count(t)
      call assemble_advection_diffusion_element_cg(i, t, matrix, rhs, &
                                        positions, old_positions, new_positions, &
                                        velocity, grid_velocity, &
                                        source, absorption, diffusivity, &
                                        density, olddensity, pressure)
    end do
    
    ! Step 4: Boundary conditions
    
    if( &
      & (integrate_advection_by_parts .and. have_advection) &
      & .or. have_diffusivity &
      & ) then
    
      allocate(t_bc_types(surface_element_count(t)))
      call get_entire_boundary_condition(t, (/ &
        "neumann      ", &
        "weakdirichlet"/), t_bc, t_bc_types)

      if(any(t_bc_types /= 0)) then
        call ewrite_bc_counts(2, t_bc_types)
      end if
    
      do i = 1, surface_element_count(t)
        call assemble_advection_diffusion_face_cg(i, t_bc_types(i), t, t_bc,  &
                                                  matrix, rhs, &
                                                  positions, velocity, grid_velocity, &
                                                  density, olddensity)
      end do
    
      call deallocate(t_bc)
      deallocate(t_bc_types)
    
    end if
    
    ewrite(2, *) "Applying strong Dirichlet boundary conditions"
    call apply_dirichlet_conditions(matrix, rhs, t, dt)
    
    ewrite_minmax(rhs%val)
    
    call deallocate(velocity)
    call deallocate(dummydensity)
    
    ewrite(1, *) "Exiting assemble_advection_diffusion_cg"
    
  end subroutine assemble_advection_diffusion_cg
  
  subroutine ewrite_bc_counts(debug_level, bc_types)
    !!< A simple subroutine to count and output the number of elements with
    !!< each boundary conditions (combines counts into a single surface
    !!< element loop).
    
    integer, intent(in) :: debug_level
    integer, dimension(:), intent(in) :: bc_types
    
    integer :: i, nneumann, nweak_dirichlet
    
    if(debug_level > current_debug_level) return
    
    nneumann = 0
    nweak_dirichlet = 0
    do i = 1, size(bc_types)
      select case(bc_types(i))
        case(BC_TYPE_NEUMANN)
          nneumann = nneumann + 1
        case(BC_TYPE_WEAKDIRICHLET)
          nweak_dirichlet = nweak_dirichlet + 1
        case(0)
        case default
          ewrite(-1, *) "For boundary condition type: ", bc_types(i)
          FLAbort("Unrecognised boundary condition type")
      end select
    end do
    
    ewrite(debug_level, *) "Surface elements with Neumann boundary condition: ", nneumann
    ewrite(debug_level, *) "Surface elements with weak Dirichlet boundary condition: ", nweak_dirichlet
    
  end subroutine ewrite_bc_counts
  
  subroutine assemble_advection_diffusion_element_cg(ele, t, matrix, rhs, &
                                      positions, old_positions, new_positions, &
                                      velocity, grid_velocity, &
                                      source, absorption, diffusivity, &
                                      density, olddensity, pressure)
    integer, intent(in) :: ele
    type(scalar_field), intent(in) :: t
    type(csr_matrix), intent(inout) :: matrix
    type(scalar_field), intent(inout) :: rhs
    type(vector_field), intent(in) :: positions
    type(vector_field), pointer :: old_positions, new_positions
    type(vector_field), intent(in) :: velocity
    type(vector_field), pointer :: grid_velocity
    type(scalar_field), intent(in) :: source
    type(scalar_field), intent(in) :: absorption
    type(tensor_field), intent(in) :: diffusivity
    type(scalar_field), intent(in) :: density
    type(scalar_field), intent(in) :: olddensity
    type(scalar_field), intent(in) :: pressure
    
    integer, dimension(:), pointer :: element_nodes
    real, dimension(ele_ngi(t, ele)) :: detwei, detwei_old, detwei_new
    real, dimension(ele_loc(t, ele), ele_ngi(t, ele), mesh_dim(t)) :: dt_t
    real, dimension(ele_loc(density, ele), ele_ngi(density, ele), mesh_dim(density)) :: drho_t
    real, dimension(ele_loc(velocity, ele), ele_ngi(velocity, ele), mesh_dim(t)) :: du_t 
    real, dimension(ele_loc(velocity, ele), ele_ngi(velocity, ele), mesh_dim(t)) :: dug_t 
    real, dimension(mesh_dim(t), mesh_dim(t), ele_ngi(t, ele)) :: j_mat 
    type(element_type) :: test_function
    type(element_type), pointer :: t_shape
    
    ! What we will be adding to the matrix and RHS - assemble these as we
    ! go, so that we only do the calculations we really need
    real, dimension(ele_loc(t, ele)) :: rhs_addto
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)) :: matrix_addto
    
#ifdef DDEBUG
    assert(ele_ngi(positions, ele) == ele_ngi(t, ele))
    assert(ele_ngi(velocity, ele) == ele_ngi(t, ele))
    if(have_diffusivity) then
      assert(ele_ngi(diffusivity, ele) == ele_ngi(t, ele))
    end if
    if(have_source) then
      assert(ele_ngi(source, ele) == ele_ngi(t, ele))
    end if
    if(have_absorption) then
      assert(ele_ngi(absorption, ele) == ele_ngi(t, ele))
    end if
    if(move_mesh) then
      ! the following has been assumed in the declarations above
      assert(ele_loc(grid_velocity, ele) == ele_loc(velocity, ele))
      assert(ele_ngi(grid_velocity, ele) == ele_ngi(velocity, ele))
    end if
#endif

    matrix_addto = 0.0
    rhs_addto = 0.0
    
    t_shape => ele_shape(t, ele)
      
    ! Step 1: Transform
    
    if(.not. have_advection .and. .not. have_diffusivity) then
      call transform_to_physical(positions, ele, detwei = detwei)
    else if(any(stabilisation_scheme == (/STABILISATION_STREAMLINE_UPWIND, STABILISATION_SUPG/))) then
      call transform_to_physical(positions, ele, t_shape, &
        & dshape = dt_t, detwei = detwei, j = j_mat)
    else
      call transform_to_physical(positions, ele, t_shape, &
        & dshape = dt_t, detwei = detwei)
    end if
    
    if(have_advection.or.(equation_type==FIELD_EQUATION_INTERNALENERGY)) then
      call transform_to_physical(positions, ele, &
           & ele_shape(velocity, ele), dshape = du_t)
      if(move_mesh) then
        call transform_to_physical(positions, ele, &
            & ele_shape(grid_velocity, ele), dshape = dug_t)
      end if
    end if
    
    if(move_mesh) then
      call transform_to_physical(old_positions, ele, detwei=detwei_old)
      call transform_to_physical(new_positions, ele, detwei=detwei_new)
    end if
    
    if(have_advection.and.(equation_type==FIELD_EQUATION_INTERNALENERGY)) then
      if(ele_shape(density, ele)==t_shape) then
        drho_t = dt_t
      else
        call transform_to_physical(positions, ele, &
          & ele_shape(density, ele), dshape = drho_t)
      end if
    end if
    
    ! Step 2: Set up test function
    
    select case(stabilisation_scheme)
      case(STABILISATION_SUPG)
        if(have_diffusivity) then
          test_function = make_supg_shape(t_shape, dt_t, ele_val_at_quad(velocity, ele), j_mat, diff_q = ele_val_at_quad(diffusivity, ele), &
            & nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
        else
          test_function = make_supg_shape(t_shape, dt_t, ele_val_at_quad(velocity, ele), j_mat, &
            & nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
        end if
      case default
        test_function = t_shape
        call incref(test_function)
    end select
    ! Important note: with SUPG the test function derivatives have not been
    ! modified - i.e. dt_t is currently used everywhere. This is fine for P1,
    ! but is not consistent for P>1.

    ! Step 3: Assemble contributions
    
    ! Mass
    if(have_mass) call add_mass_element_cg(ele, test_function, t, density, olddensity, detwei, detwei_old, detwei_new, matrix_addto, rhs_addto)
    
    ! Advection
    if(have_advection) call add_advection_element_cg(ele, test_function, t, &
                                        velocity, grid_velocity, diffusivity, &
                                        density, olddensity, &
                                        dt_t, du_t, dug_t, drho_t, detwei, j_mat, matrix_addto, rhs_addto)
        
    ! Absorption
    if(have_absorption) call add_absorption_element_cg(ele, test_function, t, absorption, detwei, matrix_addto, rhs_addto)
    
    ! Diffusivity
    if(have_diffusivity) call add_diffusivity_element_cg(ele, t, diffusivity, dt_t, detwei, matrix_addto, rhs_addto)
    
    ! Source
    if(have_source) call add_source_element_cg(ele, test_function, t, source, detwei, rhs_addto)
    
    ! Pressure
    if(equation_type==FIELD_EQUATION_INTERNALENERGY) call add_pressurediv_element_cg(ele, test_function, t, &
                                                                                  velocity, pressure, &
                                                                                  du_t, detwei, rhs_addto)
    
    ! Step 4: Insertion
            
    element_nodes => ele_nodes(t, ele)
    call addto(matrix, element_nodes, element_nodes, matrix_addto)
    call addto(rhs, element_nodes, rhs_addto)

    call deallocate(test_function)
      
  end subroutine assemble_advection_diffusion_element_cg
  
  subroutine add_mass_element_cg(ele, test_function, t, density, olddensity, detwei, detwei_old, detwei_new, matrix_addto, rhs_addto)
    integer, intent(in) :: ele
    type(element_type), intent(in) :: test_function
    type(scalar_field), intent(in) :: t
    type(scalar_field), intent(in) :: density, olddensity
    real, dimension(ele_ngi(t, ele)), intent(in) :: detwei, detwei_old, detwei_new
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)), intent(inout) :: matrix_addto
    real, dimension(ele_loc(t, ele)), intent(inout) :: rhs_addto
    
    integer :: i
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)) :: mass_matrix
    
    real, dimension(ele_ngi(density,ele)) :: density_at_quad
    
    assert(have_mass)
    
    select case(equation_type)
    case(FIELD_EQUATION_INTERNALENERGY)
      assert(ele_ngi(density, ele)==ele_ngi(olddensity, ele))
      
      density_at_quad = ele_val_at_quad(olddensity, ele)

      if(move_mesh) then
        ! needs to be evaluated at t+dt
        mass_matrix = shape_shape(test_function, ele_shape(t, ele), detwei_new*density_at_quad)
      else
        mass_matrix = shape_shape(test_function, ele_shape(t, ele), detwei*density_at_quad)
      end if
    case default
    
      if(move_mesh) then
        ! needs to be evaluated at t+dt
        mass_matrix = shape_shape(test_function, ele_shape(t, ele), detwei_new)
      else
        mass_matrix = shape_shape(test_function, ele_shape(t, ele), detwei)
      end if
      
    end select
    
    if(lump_mass) then
      do i = 1, size(matrix_addto, 1)
        matrix_addto(i, i) = matrix_addto(i, i) + sum(mass_matrix(i, :))
      end do
    else
      matrix_addto = matrix_addto + mass_matrix
    end if
  
    if(move_mesh) then
      ! In the unaccelerated form we solve:
      !  /
      !  |  N^{n+1} T^{n+1}/dt - N^{n} T^n/dt + ... = f
      !  /
      ! so in accelerated form:
      !  /
      !  |  N^{n+1} dT + (N^{n+1}- N^{n}) T^n/dt + ... = f
      !  /
      ! where dT=(T^{n+1}-T^{n})/dt is the acceleration.
      ! Put the (N^{n+1}-N^{n}) T^n term on the rhs
      mass_matrix = shape_shape(test_function, ele_shape(t, ele), (detwei_new-detwei_old))
      if(lump_mass) then
        rhs_addto = rhs_addto - sum(mass_matrix, 2)*ele_val(t, ele)/local_dt
      else
        rhs_addto = rhs_addto - matmul(mass_matrix, ele_val(t, ele))/local_dt
      end if
    end if

  end subroutine add_mass_element_cg
  
  subroutine add_advection_element_cg(ele, test_function, t, &
                                velocity, grid_velocity, diffusivity, &
                                density, olddensity, &
                                dt_t, du_t, dug_t, drho_t, detwei, j_mat, matrix_addto, rhs_addto)
    integer, intent(in) :: ele
    type(element_type), intent(in) :: test_function
    type(scalar_field), intent(in) :: t
    type(vector_field), intent(in) :: velocity
    type(vector_field), pointer :: grid_velocity
    type(tensor_field), intent(in) :: diffusivity
    type(scalar_field), intent(in) :: density, olddensity
    real, dimension(ele_loc(t, ele), ele_ngi(t, ele), mesh_dim(t)), intent(in) :: dt_t
    real, dimension(ele_loc(velocity, ele), ele_ngi(velocity, ele), mesh_dim(t)) :: du_t
    real, dimension(ele_loc(velocity, ele), ele_ngi(velocity, ele), mesh_dim(t)) :: dug_t
    real, dimension(ele_loc(density, ele), ele_ngi(density, ele), mesh_dim(density)), intent(in) :: drho_t
    real, dimension(ele_ngi(t, ele)), intent(in) :: detwei
    real, dimension(mesh_dim(t), mesh_dim(t), ele_ngi(t, ele)), intent(in) :: j_mat 
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)), intent(inout) :: matrix_addto
    real, dimension(ele_loc(t, ele)), intent(inout) :: rhs_addto
    
    real, dimension(ele_loc(t, ele), ele_loc(t,ele)) :: advection_mat
    real, dimension(velocity%dim, ele_ngi(velocity, ele)) :: velocity_at_quad
    real, dimension(ele_ngi(velocity, ele)) :: velocity_div_at_quad
    type(element_type), pointer :: t_shape
    
    real, dimension(ele_ngi(density, ele)) :: density_at_quad
    real, dimension(ele_ngi(density, ele), velocity%dim) :: densitygrad_at_quad
    real, dimension(ele_ngi(density, ele)) :: udotgradrho_at_quad
    
    assert(have_advection)
    
    t_shape => ele_shape(t, ele)
    
    velocity_at_quad = ele_val_at_quad(velocity, ele)
    if(move_mesh) then
      velocity_at_quad = velocity_at_quad - ele_val_at_quad(grid_velocity, ele)
    end if
    
    select case(equation_type)
    case(FIELD_EQUATION_INTERNALENERGY)
      assert(ele_ngi(density, ele)==ele_ngi(olddensity, ele))
      
      density_at_quad = density_theta*ele_val_at_quad(density, ele)&
                       +(1.-density_theta)*ele_val_at_quad(olddensity, ele)
      densitygrad_at_quad = density_theta*ele_grad_at_quad(density, ele, drho_t) &
                           +(1.-density_theta)*ele_grad_at_quad(olddensity, ele, drho_t)
      udotgradrho_at_quad = sum(transpose(densitygrad_at_quad)*velocity_at_quad, 1)
    end select
                
    if(integrate_advection_by_parts) then
      ! element advection matrix
      !    /                                        /
      !  - | (grad N_A dot nu) N_B dV - (1. - beta) | N_A ( div nu ) N_B dV
      !    /                                        /
      select case(equation_type)
      case(FIELD_EQUATION_INTERNALENERGY)
        advection_mat = -dshape_dot_vector_shape(dt_t, velocity_at_quad, t_shape, detwei*density_at_quad)
        if(abs(1.0 - beta) > epsilon(0.0)) then
          velocity_div_at_quad = ele_div_at_quad(velocity, ele, du_t)
          if(move_mesh) then
            velocity_div_at_quad = velocity_div_at_quad - ele_div_at_quad(grid_velocity, ele, dug_t)
          end if
          advection_mat = advection_mat &
                    - (1.0-beta) * shape_shape(test_function, t_shape, (velocity_div_at_quad*density_at_quad &
                                                                      +udotgradrho_at_quad)* detwei)
        end if
      case default
        advection_mat = -dshape_dot_vector_shape(dt_t, velocity_at_quad, t_shape, detwei)
        if(abs(1.0 - beta) > epsilon(0.0)) then
          velocity_div_at_quad = ele_div_at_quad(velocity, ele, du_t)
          if(move_mesh) then
            velocity_div_at_quad = velocity_div_at_quad - ele_div_at_quad(grid_velocity, ele, dug_t)
          end if
          advection_mat = advection_mat &
                    - (1.0-beta)*shape_shape(test_function, t_shape, velocity_div_at_quad*detwei)
        end if
      end select
    else
      ! element advection matrix
      !  /                                 /
      !  | N_A (nu dot grad N_B) dV + beta | N_A ( div nu ) N_B dV
      !  /                                 /
      select case(equation_type)
      case(FIELD_EQUATION_INTERNALENERGY)
        advection_mat = shape_vector_dot_dshape(test_function, velocity_at_quad, dt_t, detwei*density_at_quad)
        if(abs(beta) > epsilon(0.0)) then
          velocity_div_at_quad = ele_div_at_quad(velocity, ele, du_t)
          if(move_mesh) then
            velocity_div_at_quad = velocity_div_at_quad - ele_div_at_quad(grid_velocity, ele, dug_t)
          end if
          advection_mat = advection_mat &
                    + beta*shape_shape(test_function, t_shape, (velocity_div_at_quad*density_at_quad &
                                                                +udotgradrho_at_quad)*detwei)
        end if
      case default
        advection_mat = shape_vector_dot_dshape(test_function, velocity_at_quad, dt_t, detwei)
        if(abs(beta) > epsilon(0.0)) then
          velocity_div_at_quad = ele_div_at_quad(velocity, ele, du_t)
          if(move_mesh) then
            velocity_div_at_quad = velocity_div_at_quad - ele_div_at_quad(grid_velocity, ele, dug_t)
          end if
          advection_mat = advection_mat &
                    + beta*shape_shape(test_function, t_shape, velocity_div_at_quad*detwei)
        end if
      end select
    end if
    
    ! Stabilisation
    select case(stabilisation_scheme)
      case(STABILISATION_STREAMLINE_UPWIND)
        if(have_diffusivity) then
          advection_mat = advection_mat + &
            & element_upwind_stabilisation(t_shape, dt_t, velocity_at_quad, j_mat, detwei, &
            & diff_q = ele_val_at_quad(diffusivity, ele), nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
        else
          advection_mat = advection_mat + &
            & element_upwind_stabilisation(t_shape, dt_t, velocity_at_quad, j_mat, detwei, &
            & nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
        end if
      case default
    end select
      
    if(abs(dt_theta) > epsilon(0.0)) then
      matrix_addto = matrix_addto + dt_theta * advection_mat
    end if
    
    rhs_addto = rhs_addto - matmul(advection_mat, ele_val(t, ele))
    
  end subroutine add_advection_element_cg
  
  subroutine add_source_element_cg(ele, test_function, t, source, detwei, rhs_addto)
    integer, intent(in) :: ele
    type(element_type), intent(in) :: test_function
    type(scalar_field), intent(in) :: t
    type(scalar_field), intent(in) :: source
    real, dimension(ele_ngi(t, ele)), intent(in) :: detwei
    real, dimension(ele_loc(t, ele)), intent(inout) :: rhs_addto
   
    assert(have_source)
   
    rhs_addto = rhs_addto + shape_rhs(test_function, detwei * ele_val_at_quad(source, ele))
    
  end subroutine add_source_element_cg
  
  subroutine add_absorption_element_cg(ele, test_function, t, absorption, detwei, matrix_addto, rhs_addto)
    integer, intent(in) :: ele
    type(element_type), intent(in) :: test_function
    type(scalar_field), intent(in) :: t
    type(scalar_field), intent(in) :: absorption
    real, dimension(ele_ngi(t, ele)), intent(in) :: detwei
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)), intent(inout) :: matrix_addto
    real, dimension(ele_loc(t, ele)), intent(inout) :: rhs_addto
    
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)) ::  absorption_mat
    
    assert(have_absorption)
    
    absorption_mat = shape_shape(test_function, ele_shape(t, ele), detwei * ele_val_at_quad(absorption, ele))
    
    if(abs(dt_theta) > epsilon(0.0)) matrix_addto = matrix_addto + dt_theta * absorption_mat
    
    rhs_addto = rhs_addto - matmul(absorption_mat, ele_val(t, ele))
    
  end subroutine add_absorption_element_cg
  
  subroutine add_diffusivity_element_cg(ele, t, diffusivity, dt_t, detwei, matrix_addto, rhs_addto)
    integer, intent(in) :: ele
    type(scalar_field), intent(in) :: t
    type(tensor_field), intent(in) :: diffusivity
    real, dimension(ele_loc(t, ele), ele_ngi(t, ele), mesh_dim(t)), intent(in) :: dt_t
    real, dimension(ele_ngi(t, ele)), intent(in) :: detwei
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)), intent(inout) :: matrix_addto
    real, dimension(ele_loc(t, ele)), intent(inout) :: rhs_addto
    
    real, dimension(diffusivity%dim, diffusivity%dim, ele_ngi(diffusivity, ele)) :: diffusivity_gi
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)) :: diffusivity_mat
    
    assert(have_diffusivity)
    
    diffusivity_gi = ele_val_at_quad(diffusivity, ele)
    if(isotropic_diffusivity) then
      assert(size(diffusivity_gi, 1) > 0)
      diffusivity_mat = dshape_dot_dshape(dt_t, dt_t, detwei * diffusivity_gi(1, 1, :))
    else
      diffusivity_mat = dshape_tensor_dshape(dt_t, diffusivity_gi, dt_t, detwei)
    end if
    
    if(abs(dt_theta) > epsilon(0.0)) matrix_addto = matrix_addto + dt_theta * diffusivity_mat
    
    rhs_addto = rhs_addto - matmul(diffusivity_mat, ele_val(t, ele))
    
  end subroutine add_diffusivity_element_cg
  
  subroutine add_pressurediv_element_cg(ele, test_function, t, velocity, pressure, du_t, detwei, rhs_addto)
  
    integer, intent(in) :: ele
    type(element_type), intent(in) :: test_function
    type(scalar_field), intent(in) :: t
    type(vector_field), intent(in) :: velocity
    type(scalar_field), intent(in) :: pressure
    real, dimension(ele_loc(velocity, ele), ele_ngi(velocity, ele), mesh_dim(t)) :: du_t
    real, dimension(ele_ngi(t, ele)), intent(in) :: detwei
    real, dimension(ele_loc(t, ele)), intent(inout) :: rhs_addto
    
    assert(equation_type==FIELD_EQUATION_INTERNALENERGY)
    assert(ele_ngi(pressure, ele)==ele_ngi(t, ele))
    
    rhs_addto = rhs_addto - &
                shape_rhs(test_function, ele_div_at_quad(velocity, ele, du_t) * ele_val_at_quad(pressure, ele) * detwei)
    
  end subroutine add_pressurediv_element_cg
  
  subroutine assemble_advection_diffusion_face_cg(face, bc_type, t, t_bc, matrix, rhs, positions, velocity, grid_velocity, density, olddensity)
    integer, intent(in) :: face
    integer, intent(in) :: bc_type
    type(scalar_field), intent(in) :: t
    type(scalar_field), intent(in) :: t_bc
    type(csr_matrix), intent(inout) :: matrix
    type(scalar_field), intent(inout) :: rhs
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: velocity
    type(vector_field), pointer :: grid_velocity
    type(scalar_field), intent(in) :: density
    type(scalar_field), intent(in) :: olddensity
    
    integer :: ele
    integer, dimension(face_loc(t, face)) :: face_nodes
    real, dimension(face_ngi(t, face)) :: detwei
    real, dimension(mesh_dim(t), face_ngi(t, face)) :: normal
    
    ! What we will be adding to the matrix and RHS - assemble these as we
    ! go, so that we only do the calculations we really need
    real, dimension(face_loc(t, face)) :: rhs_addto
    real, dimension(face_loc(t, face), face_loc(t, face)) :: matrix_addto
    
    assert(any(bc_type == (/0, BC_TYPE_NEUMANN, BC_TYPE_WEAKDIRICHLET/)))
    assert(face_ngi(positions, face) == face_ngi(t, face))
    assert(face_ngi(velocity, face) == face_ngi(t, face))

    matrix_addto = 0.0
    rhs_addto = 0.0

    ele = face_ele(t, face)
    
    ! Step 1: Transform
    
    if(have_advection .and. integrate_advection_by_parts) then
      call transform_facet_to_physical(positions, face, &
        & detwei_f = detwei, normal = normal)
    else if(have_diffusivity.and.(bc_type == BC_TYPE_NEUMANN)) then
      call transform_facet_to_physical(positions, face, &
        & detwei_f = detwei)
    end if
    
    ! Note that with SUPG the surface element test function is not modified
          
    ! Step 2: Assemble contributions
    
    ! Advection
    if(have_advection .and. integrate_advection_by_parts) &
      call add_advection_face_cg(face, bc_type, t, t_bc, velocity, grid_velocity, density, olddensity, detwei, normal, matrix_addto, rhs_addto)
    
    ! Diffusivity
    if(have_diffusivity) call add_diffusivity_face_cg(face, bc_type, t, t_bc, detwei, rhs_addto)
    
    ! Step 3: Insertion
    
    face_nodes = face_global_nodes(t, face)
    call addto(matrix, face_nodes, face_nodes, matrix_addto)
    call addto(rhs, face_nodes, rhs_addto)
    
  end subroutine assemble_advection_diffusion_face_cg
  
  subroutine add_advection_face_cg(face, bc_type, t, t_bc, velocity, grid_velocity, density, olddensity, detwei, normal, matrix_addto, rhs_addto)
    integer, intent(in) :: face
    integer, intent(in) :: bc_type
    type(scalar_field), intent(in) :: t
    type(scalar_field), intent(in) :: t_bc
    type(vector_field), intent(in) :: velocity
    type(vector_field), pointer :: grid_velocity
    type(scalar_field), intent(in) :: density
    type(scalar_field), intent(in) :: olddensity
    real, dimension(face_ngi(t, face)), intent(in) :: detwei
    real, dimension(mesh_dim(t), face_ngi(t, face)), intent(in) :: normal
    real, dimension(face_loc(t, face), face_loc(t, face)), intent(inout) :: matrix_addto
    real, dimension(face_loc(t, face)), intent(inout) :: rhs_addto
    
    real, dimension(velocity%dim, face_ngi(velocity, face)) :: velocity_at_quad
    real, dimension(face_loc(t, face), face_loc(t,face)) :: advection_mat
    type(element_type), pointer :: t_shape
    
    real, dimension(face_ngi(density, face)) :: density_at_quad
    
    assert(have_advection)
    assert(integrate_advection_by_parts)
    
    t_shape => face_shape(t, face)
    
    velocity_at_quad = face_val_at_quad(velocity, face)
    if(move_mesh) then
      velocity_at_quad = velocity_at_quad - face_val_at_quad(grid_velocity, face)
    end if
    select case(equation_type)
    case(FIELD_EQUATION_INTERNALENERGY)
      density_at_quad = density_theta*face_val_at_quad(density, face) &
                       +(1.0-density_theta)*face_val_at_quad(olddensity, face)

      advection_mat = shape_shape(t_shape, t_shape, detwei * sum(velocity_at_quad * normal, 1) * density_at_quad)
    case default
      
      advection_mat = shape_shape(t_shape, t_shape, detwei * sum(velocity_at_quad * normal, 1))
      
    end select
    
    if(abs(dt_theta) > epsilon(0.0)) then
      if(bc_type == BC_TYPE_WEAKDIRICHLET) then
        rhs_addto = rhs_addto - theta * matmul(advection_mat, ele_val(t_bc, face) - face_val(t, face))
      else
        matrix_addto = matrix_addto + dt_theta * advection_mat
      end if
    end if
    
    rhs_addto = rhs_addto - matmul(advection_mat, face_val(t, face))

  end subroutine add_advection_face_cg
  
  subroutine add_diffusivity_face_cg(face, bc_type, t, t_bc, detwei, rhs_addto)
    integer, intent(in) :: face
    integer, intent(in) :: bc_type
    type(scalar_field), intent(in) :: t
    type(scalar_field), intent(in) :: t_bc
    real, dimension(face_ngi(t, face)), intent(in) :: detwei
    real, dimension(face_loc(t, face)), intent(inout) :: rhs_addto
    
    assert(have_diffusivity)

    if(bc_type == BC_TYPE_NEUMANN) then
      rhs_addto = rhs_addto + shape_rhs(face_shape(t, face), detwei * ele_val_at_quad(t_bc, face))
    else if(bc_type == BC_TYPE_WEAKDIRICHLET) then
      ! Need to add stuff here once transform_to_physical can supply gradients
      ! on faces to ensure that weak bcs work
      FLAbort("Weak Dirichlet boundary conditions with diffusivity are not supported by CG advection-diffusion")
    end if

  end subroutine add_diffusivity_face_cg
  
  subroutine solve_advection_diffusion_cg(t, delta_t, matrix, rhs)
    type(scalar_field), intent(in) :: t
    type(scalar_field), intent(inout) :: delta_t
    type(csr_matrix), intent(in) :: matrix
    type(scalar_field), intent(in) :: rhs
    
    call petsc_solve(delta_t, matrix, rhs, option_path = t%option_path)
    
    ewrite_minmax(delta_t%val)
    
  end subroutine solve_advection_diffusion_cg
  
  subroutine apply_advection_diffusion_cg_change(t, delta_t, dt)
    type(scalar_field), intent(inout) :: t
    type(scalar_field), intent(in) :: delta_t
    real, intent(in) :: dt
    
    ewrite_minmax(t%val)
    
    call addto(t, delta_t, dt)
    
    ewrite_minmax(t%val)
    
  end subroutine apply_advection_diffusion_cg_change
    
  subroutine advection_diffusion_cg_check_options
    !!< Check CG advection-diffusion specific options
    
    character(len = FIELD_NAME_LEN) :: field_name, state_name
    character(len = OPTION_PATH_LEN) :: path
    integer :: i, j, stat
    real :: beta, l_theta
    
    if(option_count("/material_phase/scalar_field/prognostic/spatial_discretisation/continuous_galerkin") == 0) then
      ! Nothing to check
      return
    end if
    
    ewrite(2, *) "Checking CG advection-diffusion options"
            
    if(option_count("/material_phase/scalar_field::" // advdif_cg_rhs_name) > 0) then
      FLExit("The scalar field name " // advdif_cg_rhs_name // " is reserved")
    end if
    
    if(option_count("/material_phase/scalar_field::" // advdif_cg_delta_t_name) > 0) then
      FLExit("The scalar field name " // advdif_cg_delta_t_name // " is reserved")
    end if
    
    do i = 0, option_count("/material_phase") - 1
      path = "/material_phase[" // int2str(i) // "]"
      call get_option(trim(path) // "/name", state_name)
      
      do j = 0, option_count(trim(path) // "/scalar_field") - 1
        path = "/material_phase[" // int2str(i) // "]/scalar_field[" // int2str(j) // "]"
        call get_option(trim(path) // "/name", field_name)
        
        if(field_name /= "Pressure") then
        
          path = trim(path) // "/prognostic"
          
          if(have_option(trim(path) // "/spatial_discretisation/continuous_galerkin").and.&
             have_option(trim(path) // "/equation[0]")) then       
            call get_option(trim(path) // "/spatial_discretisation/conservative_advection", beta, stat)
            if(stat == SPUD_NO_ERROR) then
              if(beta < 0.0 .or. beta > 1.0) then
              
                call field_error(state_name, field_name, &
                  & "Conservative advection factor (beta) must be >= 0.0 and <= 1.0")
              end if
            else
              call field_error(state_name, field_name, &
                & "Conservative advection factor (beta) required")
            end if
            
            call get_option(trim(path) // "/temporal_discretisation/theta", l_theta, stat)
            if(stat == SPUD_NO_ERROR) then
              if(l_theta < 0. .or. l_theta > 1.0) then
                call field_error(state_name, field_name, &
                  &"Implicitness factor (theta) must be >= 0.0 and <= 1.0")
              end if
            else
              call field_error(state_name, field_name, &
                & "Implicitness factor (theta) required")
            end if
            if(have_option(trim(path) // "/spatial_discretisation/continuous_galerkin/mass_terms/exclude_mass_terms") .and. &
              & abs(l_theta - 1.0) > epsilon(0.0)) then
              call field_warning(state_name, field_name, &
                & "Implicitness factor (theta) should = 1.0 when excluding mass")
            end if
  
            if(have_option(trim(path) // "/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind_petrov_galerkin") .and. &
              & have_option(trim(path) // "/spatial_discretisation/continuous_galerkin/advection_terms/integrate_advection_by_parts")) then
              call field_warning(state_name, field_name, &
                & "SUPG stabilisation should only be used with advection not integrated by parts")
            end if
  
            if(option_count(trim(path) // "/boundary_conditions/type::dirichlet/apply_weakly") > 0 &
              & .and. have_option(trim(path) // "/tensor_field::Diffusivity")) then
              call field_error(state_name, field_name, &
                & "Weak Dirichlet boundary conditions with diffusivity are not supported by CG advection-diffusion")
            end if
            
            if(have_option(trim(path) // "/spatial_discretisation/continuous_galerkin/advection_terms/exclude_advection_terms")) then
              if(have_option(trim(path) // "/scalar_field::SinkingVelocity")) then
                call field_warning(state_name, field_name, &
                  & "SinkingVelocity set, but advection terms have been excluded - SinkingVelocity will have no effect")
              end if
            end if
  
            if(option_count(trim(path) // "/boundary_conditions/type::neumann") > 0 &
              & .and. .not. (have_option(trim(path) // "/tensor_field::Diffusivity") &
              & .or. have_option(trim(path) // "/subgridscale_parameterisation::GLS"))) then
                call field_warning(state_name, field_name, &
                & "Neumann boundary condition set, but have no diffusivity - boundary condition will not be applied")
            end if
          end if
        end if
      end do
    end do
    
    ewrite(2, *) "Finished checking CG advection-diffusion options"
    
  contains
  
    subroutine field_warning(state_name, field_name, msg)
      character(len = *), intent(in) :: state_name
      character(len = *), intent(in) :: field_name
      character(len = *), intent(in) :: msg
      
      ewrite(0, *) "Warning: For field " // trim(field_name) // " in state " // trim(state_name)
      ewrite(0, *) trim(msg)
    
    end subroutine field_warning
  
    subroutine field_error(state_name, field_name, msg)
      character(len = *), intent(in) :: state_name
      character(len = *), intent(in) :: field_name
      character(len = *), intent(in) :: msg
      
      ewrite(-1, *) "For field " // trim(field_name) // " in state " // trim(state_name)
      FLExit(trim(msg))
    
    end subroutine field_error
  
  end subroutine advection_diffusion_cg_check_options

end module advection_diffusion_cg
