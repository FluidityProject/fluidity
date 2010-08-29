! James Maddison

#include "fdebug.h"

!! Quick and dirty shallow water solver. Based on Shallow_Water.F90 and
!! Linear_Shallow_Water.F90. This program solves the shallow water equations by
!! assembling and solving the full three dimensional system, without solving an
!! intermediate wave equation. This is slow, but works with more femtools
!! supported element types.
program simple_shallow_water

#ifndef HAVE_PETSC
#error PETSc required
#endif
#ifdef HAVE_PETSC_MODULES
  use petsc 
  use petscvec 
  use petscmat 
  use petscksp 
  use petscpc
#else
#error PETSc modules required
#endif
#ifdef HAVE_ZOLTAN
  use zoltan
#endif

  use adapt_state_prescribed_module
  use checkpoint
  use coriolis_module
  use diagnostic_fields_new, only : calculate_diagnostic_variables_new => &
    & calculate_diagnostic_variables, check_diagnostic_dependencies
  use diagnostic_fields_wrapper
  use diagnostic_variables
  use fields
  use fldebug
  use global_parameters, only : current_time, dt, OPTION_PATH_LEN, &
    & simulation_start_cpu_time, simulation_start_wall_time, timestep
  use memory_diagnostics
  use parallel_tools
  use populate_state_module
  use reference_counting
  use reserve_state_module
  use signals
  use solvers
  use sparse_matrices_fields
  use sparse_tools
  use sparsity_patterns_meshes
  use spud
  use state_fields_module
  use state_module
  use tictoc
  use timeloop_utilities
  use timers
  use transform_elements
  use vtk_interfaces
  use write_state_module

  implicit none
  
  interface check_options
    subroutine check_options()
    end subroutine check_options
  end interface check_options
  
  interface petscfinalize
    subroutine petscfinalize(ierr)
      implicit none
      integer, intent(out) :: ierr
    end subroutine petscfinalize
  end interface petscfinalize
  
  interface petscinitialize
    subroutine petscinitialize(c, ierr)
      implicit none
      character(len = *), intent(in) :: c
      integer, intent(out) :: ierr
    end subroutine petscinitialize
  end interface petscinitialize
  
  interface python_end
    subroutine python_end()
    end subroutine python_end
  end interface python_end
  
  interface python_init
    subroutine python_init()
    end subroutine python_init
  end interface python_init
  
  interface set_global_debug_level
    subroutine set_global_debug_level(level)
      implicit none
      integer, intent(in) :: level
    end subroutine set_global_debug_level
  end interface set_global_debug_level
  
  character(len = *), parameter :: temp_solver_path = "/temporary/solver/path"
  
#ifdef HAVE_ZOLTAN
  real(zoltan_float) :: ver
  integer(zoltan_int) :: zierr
#endif
  integer :: ierr
 
  integer :: adapt_no, dump_no, nonlinear_iteration, nonlinear_iterations
  real :: last_dump_timestep
  
  character(len = OPTION_PATH_LEN) :: simulation_name
  logical :: lump_m1, slump_m1, lump_m2, slump_m2, lump_l, include_advection
  real :: beta, g, h, theta
  
  type(block_csr_matrix) :: c_m, ct_m
  type(csr_matrix) :: linear_big_m, big_m, l_m, n_m
  type(csr_matrix), pointer :: m1, m2
  type(scalar_field) :: big_solution, big_rhs, lumped_l_m
  type(scalar_field), pointer :: lumped_m1, lumped_m2
  type(vector_field) :: n_rhs
  
  type(scalar_field), pointer :: p
  type(state_type), pointer :: state
  type(state_type), dimension(:), pointer :: states
  type(vector_field), pointer :: u, u_nl, x
  
  ! Initialise external libraries
#ifdef HAVE_MPI
  call mpi_init(ierr)
  assert(ierr == MPI_SUCCESS)
#endif
#ifdef HAVE_PETSC
  call petscinitialize(PETSC_NULL_CHARACTER, ierr)
#endif
  call python_init()
#ifdef HAVE_ZOLTAN
  zierr = zoltan_initialize(ver)  
  assert(zierr == ZOLTAN_OK)
#endif

  ! Initialise the signal handlers
  call initialise_signals()

  ! Initialise the model timer
  call tictoc_reset()
  call tic(TICTOC_ID_SIMULATION)

  ! Read and check input options
  call read_command_line()
  call simple_shallow_water_check_options()
  
  ! Populate the system state and read options
  call populate_state(states)
  call auxilliary_state_allocation()
  call read_options()
  
  ! Perform system state checks
  call initialise_pointers()
  call simple_shallow_water_check_state()
    
  ! Remaining initialisation
  timestep = 0;  adapt_no = 0;  dump_no = 0
#ifdef HAVE_MEMORY_STATS
  call reset_memory_logs()
#endif
  call initialise_diagnostics(simulation_name, states)
  call initialise_write_state()
  
  ! Assemble the system matrices
  call initialise_matrices()
  call initialise_timestep()
  
  if(have_option("/material_phase::Fluid/vector_field::Velocity/prognostic/initial_condition::WholeMesh/balanced")) then
    ! Balanced Velocity initialisation
    call initialise_balanced_velocity()
    if(nonlinear_iterations > 1) call copy_to_stored_values(states, "Old")
  end if

  call calculate_diagnostic_variables(states, exclude_nonrecalculated = .false.)
  call calculate_diagnostic_variables_new(states, exclude_nonrecalculated = .false.)
  call write_diagnostics(states, current_time, dt, timestep)
  call write_state(dump_no, states);  last_dump_timestep = timestep

  ! Execute the timestep loop
  if(.not. simulation_completed(current_time, timestep)) then
    timestep_loop: do
      timestep = timestep + 1
      
      ewrite(1, *) "Executing timestep: ", timestep
      ewrite(1, *) "Model time: ", current_time
      
      nonlinear_loop: do nonlinear_iteration = 1, nonlinear_iterations
        ewrite(1, *) "Executing non-linear iteration ", nonlinear_iteration, " of ", nonlinear_iterations
      
        if(nonlinear_iterations > 1) then
          call copy_to_stored_values(states, "Iterated")
          call relax_to_nonlinear(states)
          if(nonlinear_iteration > 1) call copy_from_stored_values(states, "Old")
        end if
        
        call assemble_big_m()
        call assemble_big_rhs()
        call petsc_solve(big_solution, big_m, big_rhs, option_path = trim(p%option_path) // "/prognostic")
        call increment_fields()
      
      end do nonlinear_loop
      
      if(nonlinear_iterations > 1) call copy_to_stored_values(states, "Old")
        
      current_time = current_time + dt
      
      call calculate_diagnostic_variables(states, exclude_nonrecalculated = .true.)
      call calculate_diagnostic_variables_new(states, exclude_nonrecalculated = .true.)
      call write_diagnostics(states, current_time, dt, timestep)
      
      if(do_write_state(current_time, timestep)) then
         if(do_checkpoint_simulation(dump_no)) then
          call checkpoint_simulation(states, cp_no = dump_no)
        end if
        call write_state(dump_no, states);  last_dump_timestep = timestep
      end if
      
      if(simulation_completed(current_time, timestep + 1)) exit timestep_loop
      
      if(have_option("/mesh_adaptivity/prescribed_adaptivity")) then
        if(do_adapt_state_prescribed(current_time)) then   
          ! Execute a prescribed adapt
          call prescribed_adapt()
        end if
      end if
    end do timestep_loop
  end if
 
  if(last_dump_timestep /= timestep) then
    if(have_option("/io/checkpointing/checkpoint_at_end")) then
      call checkpoint_simulation(states, cp_no = dump_no)
    end if
    call write_state(dump_no, states)
  end if
  
  ! Cleanup
  call deallocate_matrices()
  call finalise_timestep()
  call deallocate(states)
  deallocate(states)
  call close_diagnostic_files()
  call deallocate_reserve_state()
  call deallocate_transform_cache()
  
  ! Memory leak reporting
  call print_references(0)
#ifdef HAVE_MEMORY_STATS
  call print_current_memory_stats(0)
#endif

  ! Stop and report the model timer
  call toc(TICTOC_ID_SIMULATION)
  call tictoc_report(1, TICTOC_ID_SIMULATION)

  ! Model interrupt
  if(SIG_INT) then
    FLExit("Interrupt signal received")
  end if

  ! Finalise external libraries
  call python_end()
#ifdef HAVE_PETSC
  call petscfinalize(ierr)
#endif
#ifdef HAVE_MPI
  call mpi_finalize(ierr)
  assert(ierr == MPI_SUCCESS)
#endif

contains

! --- From Shallow_Water.F90 (with edition) ---

  subroutine read_command_line()
  
    character(len = OPTION_PATH_LEN) :: argument
    integer :: status, argi, argc, level

    level = 0

    argc = command_argument_count()
    if(argc == 0) then
      goto 42
    end if
    
    argi = 1
    do while(argi <= argc - 1)
       call get_command_argument(argi, value = argument, status = status)
       argi = argi + 1
       if(status /= 0) then
         ewrite(-1, *) "Failed to read argument: ", argi
         goto 666
       end if

       select case(argument)
        case("-h", "--help")
          goto 42
        case("-v")
          if(argi > argc - 1) then
            ewrite(-1, *) "Failed to find verbosity value"
            goto 666
          end if
          call get_command_argument(argi, value = argument, status = status)
          argi = argi + 1
          if(status /= 0) then
            ewrite(-1, *) "Failed to read argument: ", argi
            goto 666
          end if

          read(argument, "(i1)", err = 666) level
          level = max(level, 0)
          level = min(level, 3)
        case("-v0")
          level = 0
        case("-v1")
          level = 1
        case("-v2")
          level = 2
        case("-v3")
          level = 3
        case default
          ewrite(-1, *) "Unrecognised argument: " // trim(argument)
          goto 666
       end select
    end do
    
    call set_global_debug_level(level)
    
    call get_command_argument(argc, value = argument, status = status)
    if(status /= 0) then
      ewrite(-1, *) "Failed to read argument: ", argc
      goto 666
    else if(trim(argument) == "-h" .or. trim(argument) == "--help") then
      goto 42
    end if

    call load_options(argument)
    if(.not. have_option("/simulation_name")) then
      ewrite(-1, *) "Failed to find simulation_name"
      goto 666
    end if

    return

42  call usage()
    stop

666 call usage()
    FLExit("Failed to read command line arguments")

  end subroutine read_command_line

  subroutine usage

    print *, "Usage: linear_shallow_water [-h] [-v LEVEL] SSWML"
    print *, ""
    print *, "-h        Display this help"
    print *, "-v LEVEL  Verbosity level"
    
  end subroutine usage
  
! ---End of from Shallow_Water.F90 ---

  subroutine read_options()
    !!< Read model parameters
  
    call get_option("/simulation_name", simulation_name)
    
    ! Timestep options
    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/timestep", dt)
    assert(dt > 0.0)
    call cpu_time(simulation_start_cpu_time)
    simulation_start_wall_time = wall_time()
    call get_option("/timestepping/nonlinear_iterations", nonlinear_iterations)
    assert(nonlinear_iterations >= 0)
    
    ! Physical parameters
    call get_option("/physical_parameters/gravity/magnitude", g)
    call get_option("/material_phase::Fluid/scalar_field::LayerThickness/prognostic/mean_layer_thickness", h)
    call get_option("/material_phase::Fluid/scalar_field::LayerThickness/prognostic/temporal_discretisation/theta", theta)
    assert(theta >= 0.0 .and. theta <= 1.0)
    
    ! Discretisation options
    lump_m1 = have_option("/material_phase::Fluid/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/mass_terms/lump_mass_matrix") &
      & .or. have_option("/material_phase::Fluid/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/mass_terms/lump_mass_matrix")
    slump_m1 = have_option("/material_phase::Fluid/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/mass_terms/lump_mass_matrix/use_submesh") &
      & .or. have_option("/material_phase::Fluid/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/mass_terms/lump_mass_matrix/use_submesh")
       
    lump_l = have_option("/material_phase::Fluid/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/coriolis_terms/lump_mass_matrix") &
      & .or. have_option("/material_phase::Fluid/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/coriolis_terms/lump_mass_matrix")

    lump_m2 = have_option("/material_phase::Fluid/scalar_field::LayerThickness/prognostic/spatial_discretisation/continuous_galerkin/mass_terms/lump_mass_matrix")
    slump_m2 = have_option("/material_phase::Fluid/scalar_field::LayerThickness/prognostic/spatial_discretisation/continuous_galerkin/mass_terms/lump_mass_matrix/use_submesh")
      
    include_advection = .not. have_option("/material_phase::Fluid/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/advection_terms/exclude_advection_terms") &
      & .and. .not. have_option("/material_phase::Fluid/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/advection_scheme/none")
    call get_option("/material_phase::Fluid/vector_field::Velocity/prognostic/spatial_discretisation/conservative_advection", beta, default = 1.0)
    assert(beta >= 0.0 .and. beta <= 1.0)
          
  end subroutine read_options
  
  subroutine auxilliary_state_allocation()
    !!< Initialise all auxilliary fields in state
    
    call allocate_and_insert_auxilliary_fields(states)
    if(nonlinear_iterations > 1) then
      call copy_to_stored_values(states, "Old")
      call copy_to_stored_values(states, "Iterated")
      call relax_to_nonlinear(states)
    end if
  
  end subroutine auxilliary_state_allocation
  
  subroutine initialise_pointers()
    !!< Initialise state and field pointers
  
    state => states(1)
    x => extract_vector_field(state, "Coordinate")
    p => extract_scalar_field(state, "LayerThickness")  
    u => extract_vector_field(state, "Velocity")
    if(nonlinear_iterations > 1) then
      u_nl => extract_vector_field(state, "NonlinearVelocity")
    else
      u_nl => u
    end if
  
  end subroutine initialise_pointers
  
  subroutine initialise_matrices()
    !!< Allocate the individual term matrices. Assemble the linear term
    !!< matrices.
  
    integer :: i
    type(csr_sparsity), pointer :: m1_sparsity, c_sparsity, ct_sparsity
  
    ! Velocity mass matrix
    if(lump_m1) then
      ewrite(2, *) "Using lumped Velocity mass matrix"
      if(slump_m1) then
        ewrite(2, *) "Lumping Velocity mass matrix on sub-mesh"
        lumped_m1 => get_lumped_mass_on_submesh(state, u%mesh)
      else
        lumped_m1 => get_lumped_mass(state, u%mesh)
      end if
    else
      ewrite(2, *) "Using consistent Velocity mass matrix"
      m1 => get_mass_matrix(state, u%mesh)
    end if
    m1_sparsity => get_csr_sparsity_firstorder(state, u%mesh, u%mesh)
    
    ! Advection matrix
    if(include_advection) then 
      call allocate(n_m, m1_sparsity, name = "N")
      call allocate(n_rhs, 2, u%mesh, name = "NRHS")
    end if
    
    ! Coriolis matrix
    if(lump_l) then
      ewrite(2, *) "Using lumped Coriolis"
      call allocate(lumped_l_m, u%mesh, name = "LumpedL")
      call zero(lumped_l_m)
    else
      call allocate(l_m, m1_sparsity, name = "L")
      call zero(l_m)
    end if
    
    ! Layer thickness mass matrix
    if(lump_m2) then
      ewrite(2, *) "Using lumped LayerThickness mass matrix"
      if(slump_m2) then
        ewrite(2, *) "Lumping LayerThickness mass matrix on sub-mesh"
        lumped_m2 => get_lumped_mass_on_submesh(state, p%mesh)
      else
        lumped_m2 => get_lumped_mass(state, p%mesh)
      end if       
    else
      ewrite(2, *) "Using consistent LayerThickness mass matrix"
      m2 => get_mass_matrix(state, p%mesh)
    end if
    
    ! Divergence matrix
    ct_sparsity => get_csr_sparsity_firstorder(state, p%mesh, u%mesh)
    call allocate(ct_m, ct_sparsity, blocks = (/1, 2/), name = "CT")    
    call zero(ct_m)
    
    ! Gradient matrix
    c_sparsity => get_csr_sparsity_firstorder(state, u%mesh, p%mesh)
    call allocate(c_m, c_sparsity, blocks = (/2, 1/), name = "C")
    call zero(c_m)
    
    do i = 1, ele_count(u)
      call assemble_matrices_ele(i, x, u%mesh, p%mesh, l_m, lumped_l_m, ct_m, c_m)
    end do
    
    if(lump_l) then
      ewrite_minmax(lumped_l_m%val)
    else
      ewrite_minmax(l_m%val)
    end if
    ewrite_minmax(ct_m%val(1, 1)%ptr)
    ewrite_minmax(ct_m%val(1, 2)%ptr)
    ewrite_minmax(c_m%val(1, 1)%ptr)
    ewrite_minmax(c_m%val(2, 1)%ptr)
    
  end subroutine initialise_matrices
  
  subroutine assemble_matrices_ele(ele, positions, u_mesh, p_mesh, l_m, lumped_l_m, ct_m, c_m)
    !!< Assemble the linear term matrices
  
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(mesh_type), intent(in) :: u_mesh
    type(mesh_type), intent(in) :: p_mesh
    type(csr_matrix), intent(inout) :: l_m
    type(scalar_field), intent(inout) :: lumped_l_m
    type(block_csr_matrix), intent(inout) :: ct_m
    type(block_csr_matrix), intent(inout) :: c_m
    
    integer, dimension(:), pointer :: u_nodes, p_nodes
    real, dimension(ele_ngi(positions, ele)) :: detwei, f_gi
    real, dimension(ele_loc(u_mesh, ele), ele_loc(u_mesh, ele)) :: little_l
    real, dimension(ele_loc(p_mesh, ele), ele_ngi(positions, ele), positions%dim) :: dp_t
    type(element_type), pointer :: u_shape
    
    u_shape => ele_shape(u_mesh, ele)
    u_nodes => ele_nodes(u_mesh, ele)
    p_nodes => ele_nodes(p_mesh, ele)
    
    call transform_to_physical(positions, ele, ele_shape(p_mesh, ele), &
      & dshape = dp_t, detwei = detwei)
    
    f_gi = coriolis(ele_val_at_quad(positions, ele))
    little_l = shape_shape(u_shape, u_shape, detwei * f_gi)
    if(lump_l) then
      call addto(lumped_l_m, u_nodes, sum(little_l, 2))
    else
      call addto(l_m, u_nodes, u_nodes, little_l)
    end if
    
    call addto(ct_m, p_nodes, u_nodes, &
      & spread(-dshape_shape(dp_t, u_shape, detwei), 1, 1))
    ! This is inefficient and really not required. Present in the absence of a
    ! convenient block_csr transpose.
    call addto(c_m, u_nodes, p_nodes, &
      & spread(-shape_dshape(u_shape, dp_t, detwei), 2, 1))
    
  end subroutine assemble_matrices_ele
  
  subroutine initialise_timestep()
    !!< Assemble the linear timestep matrix and allocate the timestep matrix
    !!< and solution vector
  
    integer :: i
    type(csr_sparsity), pointer :: uu_sparsity, up_sparsity, pu_sparsity, pp_sparsity
    type(csr_matrix), dimension(3, 3) :: e_big_m
    type(mesh_type) :: big_mesh
    
    ! Step 1: Assemble the linear system matrix as an array of CSR matrices
    
    uu_sparsity => get_csr_sparsity_firstorder(state, u%mesh, u%mesh)
    up_sparsity => get_csr_sparsity_firstorder(state, u%mesh, p%mesh)
    pu_sparsity => get_csr_sparsity_firstorder(state, p%mesh, u%mesh)
    pp_sparsity => get_csr_sparsity_firstorder(state, p%mesh, p%mesh)
    
    ! Optimisation possible here when using lumped matrices
    call allocate(e_big_m(1, 1), uu_sparsity, name = "SystemMatrix11")
    call allocate(e_big_m(1, 2), uu_sparsity, name = "SystemMatrix12")
    call allocate(e_big_m(1, 3), up_sparsity, name = "SystemMatrix13")
    call allocate(e_big_m(2, 1), uu_sparsity, name = "SystemMatrix21")
    call allocate(e_big_m(2, 3), up_sparsity, name = "SystemMatrix23")
    call allocate(e_big_m(3, 1), pu_sparsity, name = "SystemMatrix31")
    call allocate(e_big_m(3, 2), pu_sparsity, name = "SystemMatrix32")
    call allocate(e_big_m(3, 3), pp_sparsity, name = "SystemMatrix33")
        
    if(lump_m1) then                                         ! u   u   component (u mass)
      call zero(e_big_m(1, 1))
      do i = 1, node_count(u)
         call addto_diag(e_big_m(1, 1), i, node_val(lumped_m1, i))
      end do
    else
      e_big_m(1, 1)%val = m1%val
    end if
    if(lump_l) then                                          ! u   v   component (x Coriolis)
      call zero(e_big_m(1, 2))
      do i = 1, node_count(u)
        call addto_diag(e_big_m(1, 2), i, dt * theta * (-node_val(lumped_l_m, 1)))
      end do
    else
      e_big_m(1, 2)%val = dt * theta * (-l_m%val)
    end if
    e_big_m(1, 3)%val = -g * dt * theta * c_m%val(1, 1)%ptr  ! u   eta component (x gradient)
    if(lump_l) then                                          ! v   u   component (y Coriolis)
      call zero(e_big_m(2, 1))
      do i = 1, node_count(u)
        call addto_diag(e_big_m(2, 1), i, dt * theta * node_val(lumped_l_m, 1))
      end do
    else
      e_big_m(2, 1)%val = dt * theta * l_m%val
    end if
    e_big_m(2, 2) = e_big_m(1, 1)                            ! v   v   component (u mass)
    e_big_m(2, 3)%val = -g * dt * theta * c_m%val(2, 1)%ptr  ! v   eta component (y gradient)
    e_big_m(3, 1)%val = h * dt * theta * ct_m%val(1, 1)%ptr  ! eta u   component (x divergence)
    e_big_m(3, 2)%val = h * dt * theta * ct_m%val(1, 2)%ptr  ! eta v   component (y divergence)
    if(lump_m2) then                                         ! eta eta component (eta mass)     
      call zero(e_big_m(3, 3))
      do i = 1, node_count(p)
         call addto_diag(e_big_m(3, 3), i, node_val(lumped_m2, i))
      end do   
    else
      e_big_m(3, 3)%val = m2%val
    end if
    
    ! Step 2: Flatten the array of CSR matrices to form a single, big, CSR
    ! matrix
    
    call flatten_big_m(node_count(u), node_count(p), e_big_m, linear_big_m)
    
    ! Step 3: Cleanup
    
    call deallocate(e_big_m(1, 1))
    call deallocate(e_big_m(1, 2))
    call deallocate(e_big_m(1, 3))
    call deallocate(e_big_m(2, 1))
    call deallocate(e_big_m(2, 3))
    call deallocate(e_big_m(3, 1))
    call deallocate(e_big_m(3, 2))
    call deallocate(e_big_m(3, 3))
    
    ! Step 4: Allocate the system matrix, RHS and solution
    
    if(include_advection) then
      call allocate(big_m, linear_big_m%sparsity, name = "SystemMatrix")
    else
      big_m = linear_big_m
      call incref(big_m)
    end if
    
    ! This is a dummy mesh, for the purposes of the PETSc interfaces
    call allocate(big_mesh, nodes = 2 * node_count(u) + node_count(p), elements = 0, shape = ele_shape(u, 1), name = "SystemMesh")
    call allocate(big_rhs, big_mesh, name = "SystemRHS")
    call allocate(big_solution, big_mesh, name = "SystemSolution")
    call zero(big_solution)
    call deallocate(big_mesh)
  
  end subroutine initialise_timestep
  
  subroutine assemble_big_m()
    !!< Assemble the system matrix
        
    integer :: i, un
        
    if(include_advection) then
      select case(continuity(u))
        case(-1)
          call assemble_n_dg()
        case(0)
          FLExit("Continuous Galerkin advection not supported")
        case default
          ewrite(-1, *) "For mesh continuity: ", continuity(u)
          FLAbort("Unrecognised mesh continuity")
      end select
        
      ! Linear terms
      call set(big_m, linear_big_m)
      
      ! Non-linear advection
      un = node_count(u)
      do i = 1, un
        call addto(big_m, (/i/), row_m_ptr(n_m, i), spread(row_val(n_m, i), 1, 1))
        call addto(big_m, (/i + un/), row_m_ptr(n_m, i) + un, spread(row_val(n_m, i), 1, 1))
      end do
    !else
    !  ! Only linear terms
    end if
  
  end subroutine assemble_big_m
  
  subroutine assemble_n_dg()
    !!< Assemble the advection terms using a discontinuous Galerkin formulation
    
    integer :: i
    
    ewrite(1, *) "Assembling advection terms (DG)"
    
    call zero(n_m)
    call zero(n_rhs)
    
    do i = 1, ele_count(u)
      call assemble_n_dg_ele(i, x, u_nl, n_m, n_rhs)
    end do
    
    ewrite_minmax(n_m%val)
    ewrite(2, *) "sum(n_m%val) = ", sum(n_m%val)
    ewrite_minmax(n_rhs%val(1)%ptr)
    ewrite(2, *) "sum(n_rhs%val(1)%ptr) = ", sum(n_rhs%val(1)%ptr)
    ewrite_minmax(n_rhs%val(2)%ptr)
    ewrite(2, *) "sum(n_rhs%val(2)%ptr) = ", sum(n_rhs%val(2)%ptr)
    
    ewrite(1, *) "Finished assembling advection terms"
  
  end subroutine assemble_n_dg
  
  subroutine assemble_n_dg_ele(ele, positions, u, n_m, n_rhs)
    !!< Assemble the advection terms using a discontinuous Galerkin formulation
    
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: u
    type(csr_matrix), intent(inout) :: n_m
    type(vector_field), intent(inout) :: n_rhs
        
! --- From Advection_Diffusion_DG.F90 (with edition) ---
    integer :: ele_2, face, face_2, i
    integer, dimension(:), pointer :: neigh, nodes
    real, dimension(ele_ngi(positions, ele)) :: detwei
    real, dimension(ele_loc(u, ele), ele_ngi(positions, ele), positions%dim) :: du_t
    real, dimension(ele_loc(u, ele), ele_loc(u, ele)) :: little_n
    real, dimension(ele_ngi(u, ele)) :: div_u_gi
    real, dimension(u%dim, ele_ngi(u, ele)) :: u_gi
    type(element_type), pointer :: u_shape
    
    u_shape => ele_shape(u, ele)
    nodes => ele_nodes(u, ele)
    
    call transform_to_physical(positions, ele, u_shape, &
      & dshape = du_t, detwei = detwei)
    
    u_gi = ele_val_at_quad(u, ele)    
    div_u_gi = ele_div_at_quad(u, ele, du_t)
    
    little_n = -dshape_dot_vector_shape(du_t, u_gi, u_shape, detwei)
    if(abs(1.0 - beta) > epsilon(0.0)) then
      little_n = little_n - (1.0 - beta) * shape_shape(u_shape, u_shape, detwei * div_u_gi)
    end if
    
    call addto(n_m, nodes, nodes, dt * theta * little_n)
    call addto(n_rhs, 1, nodes, -matmul(little_n, ele_val(u, 1, ele)))
    call addto(n_rhs, 2, nodes, -matmul(little_n, ele_val(u, 2, ele)))
    
    neigh => ele_neigh(u, ele)
    do i = 1, size(neigh)
      ele_2 = neigh(i)
      face = ele_face(u, ele, ele_2)
      if(ele_2 > 0) then
        ! Internal
        face_2 = ele_face(u, ele_2, ele)
      else
        ! External
        face_2 = face
      end if
      
      call assemble_n_dg_face(face, face_2, positions, u, n_m)
    end do
! ---End of from Advection_Diffusion_DG.F90 ---
    
  end subroutine assemble_n_dg_ele
  
  subroutine assemble_n_dg_face(face, face_2, positions, u, n_m)
    !!< Assemble the advection terms using a discontinuous Galerkin formulation
  
    integer, intent(in) :: face
    integer, intent(in) :: face_2
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: u
    type(csr_matrix), intent(inout) :: n_m
  
! --- From Advection_Diffusion_DG.F90 (with edition) ---
    integer, dimension(face_loc(u, face)) :: nodes
    integer, dimension(face_loc(u, face_2)) :: nodes_2
    real, dimension(face_ngi(positions, face)) :: detwei
    real, dimension(face_ngi(u, face)) :: inflow, u_gi_n
    real, dimension(u%dim, face_ngi(u, face)) :: normal, u_gi, u_gi_2, u_gi_bar
    real, dimension(face_loc(u, face), face_loc(u, face)) :: little_n_out
    real, dimension(face_loc(u, face), face_loc(u, face_2)) :: little_n_in
    type(element_type), pointer :: u_shape, u_shape_2
    
    u_shape => face_shape(u, face)
    u_shape_2 => face_shape(u, face_2)
    nodes = face_global_nodes(u, face)
    nodes_2 = face_global_nodes(u, face_2)
    
    call transform_facet_to_physical(positions, max(face, face_2), &
      & detwei_f = detwei, normal = normal)
    if(face_2 > face) normal = -normal

    u_gi = face_val_at_quad(u, face)
    u_gi_2 = face_val_at_quad(u, face_2)
    u_gi_bar = 0.5 * (u_gi + u_gi_2)
    
    u_gi_n = sum(u_gi_bar * normal, 1)
    inflow = merge(1.0, 0.0, u_gi_n < 0.0)
    
    little_n_out = shape_shape(u_shape, u_shape,  detwei * (1.0 - inflow) * u_gi_n) 
    little_n_in = shape_shape(u_shape, u_shape_2, detwei * inflow * u_gi_n) 

    ! Outflow boundary integral
    call addto(n_m, nodes, nodes, dt * theta * little_n_out)
    ! Inflow boundary integral
    call addto(n_m, nodes, nodes_2, dt * theta * little_n_in)
    call addto(n_rhs, 1, nodes, &
      ! Outflow boundary integral
      & -matmul(little_n_out,face_val(u, 1, face)) &
      ! Inflow boundary integral
      & -matmul(little_n_in, face_val(u, 1, face_2)))
    call addto(n_rhs, 2, nodes, &
      ! Outflow boundary integral
      & -matmul(little_n_out,face_val(u, 2, face)) &
      ! Inflow boundary integral
      & -matmul(little_n_in, face_val(u, 2, face_2)))
! ---End of from Advection_Diffusion_DG.F90 ---
  
  end subroutine assemble_n_dg_face
    
  subroutine assemble_big_rhs()
    !!< Assemble the system RHS
  
    type(scalar_field) :: u_x, u_y, big_rhs_u_addto
    type(scalar_field), dimension(3) :: e_big_rhs
    
    ! Split up the system RHS into an array of scalar fields    
    e_big_rhs(1) = wrap_scalar_field(u%mesh, big_rhs%val(:node_count(u)), name = "SystemRHS1")
    e_big_rhs(2) = wrap_scalar_field(u%mesh, big_rhs%val(node_count(u) + 1:2 * node_count(u)), name = "SystemRHS2")
    e_big_rhs(3) = wrap_scalar_field(p%mesh, big_rhs%val(2 * node_count(u) + 1:), name = "SystemRHS3")
    
    call allocate(big_rhs_u_addto, u%mesh, "VelocityAddto")
    u_x = extract_scalar_field(u, 1)
    u_y = extract_scalar_field(u, 2)
    
    ! Advection
    if(include_advection) then
      call addto(e_big_rhs(1), extract_scalar_field(n_rhs, 1), scale = dt)
      call addto(e_big_rhs(2), extract_scalar_field(n_rhs, 2), scale = dt)
    end if
    
    ! Coriolis
    if(lump_l) then
      e_big_rhs(1)%val = dt * lumped_l_m%val * u_y%val
      e_big_rhs(2)%val = dt * lumped_l_m%val * (-u_x%val)
    else
      call mult(big_rhs_u_addto, l_m, u_y);  call scale(big_rhs_u_addto, dt);   call set(e_big_rhs(1), big_rhs_u_addto)
      call mult(big_rhs_u_addto, l_m, u_x);  call scale(big_rhs_u_addto, -dt);  call set(e_big_rhs(2), big_rhs_u_addto)
    end if
    
    ! Gradient
    call mult(big_rhs_u_addto, block(c_m, 1, 1), p);  call scale(big_rhs_u_addto, dt * g);  call addto(e_big_rhs(1), big_rhs_u_addto)
    call mult(big_rhs_u_addto, block(c_m, 2, 1), p);  call scale(big_rhs_u_addto, dt * g);  call addto(e_big_rhs(2), big_rhs_u_addto)

    call deallocate(big_rhs_u_addto)    
    
    ! Divergence
    call mult(e_big_rhs(3), ct_m, u);  call scale(e_big_rhs(3), -dt * h)
    
    call deallocate(e_big_rhs(1))
    call deallocate(e_big_rhs(2))
    call deallocate(e_big_rhs(3))
      
  end subroutine assemble_big_rhs
  
  subroutine flatten_big_m(un, pn, e_big_m, big_m)
    !!< Flatten the array of system CSR matrix blocks to form a single system
    !!< CSR matrix
  
    integer, intent(in) :: un
    integer, intent(in) :: pn
    type(csr_matrix), dimension(3, 3), intent(in) :: e_big_m
    type(csr_matrix), intent(out) :: big_m
  
    integer :: i, n, col_offset, row_offset
    type(dynamic_csr_matrix) :: d_big_m
    
    ewrite(1, *) "Flattening system matrix"
    
    ! Computing the sparsity pattern of the full system matrix is hard, so we
    ! make use of a dynamic CSR matrix here. This is going to be slow...
    
    ! Allocate a dynamic CSR matrix object
    n = 2 * un + pn
    call allocate(d_big_m, rows = n, columns = n, name = "SystemMatrix")
    
    ! Flatten each block in turn
    
    ewrite(2, *) "u   u   component"
    do i = 1, un
      call set(d_big_m, i, row_m_ptr(e_big_m(1, 1), i), row_val(e_big_m(1, 1), i))
    end do
  
    ewrite(2, *) "u   v   component"
    col_offset = un
    do i = 1, un
      call set(d_big_m, i, row_m_ptr(e_big_m(1, 2), i) + col_offset, row_val(e_big_m(1, 2), i))
    end do
    
    ewrite(2, *) "u   eta component"
    col_offset = 2 * un
    do i = 1, un
      call set(d_big_m, i, row_m_ptr(e_big_m(1, 3), i) + col_offset, row_val(e_big_m(1, 3), i))
    end do
    
    ewrite(2, *) "v   u   component"
    row_offset = un
    do i = 1, un
      call set(d_big_m, i + row_offset, row_m_ptr(e_big_m(2, 1), i), row_val(e_big_m(2, 1), i))
    end do 
    
    ewrite(2, *) "v   v   component"
    col_offset = un    
    do i = 1, un
      call set(d_big_m, i + row_offset, row_m_ptr(e_big_m(2, 2), i) + col_offset, row_val(e_big_m(2, 2), i))
    end do
    
    ewrite(2, *) "v   eta component"
    col_offset = 2 * un    
    do i = 1, un
      call set(d_big_m, i + row_offset, row_m_ptr(e_big_m(2, 3), i) + col_offset, row_val(e_big_m(2, 3), i))
    end do
    
    ewrite(2, *) "eta u   component"
    row_offset = 2 * un
    do i = 1, pn
      call set(d_big_m, i + row_offset, row_m_ptr(e_big_m(3, 1), i), row_val(e_big_m(3, 1), i))
    end do  
    
    ewrite(2, *) "eta v   component"
    col_offset = un    
    do i = 1, pn
      call set(d_big_m, i + row_offset, row_m_ptr(e_big_m(3, 2), i) + col_offset, row_val(e_big_m(3, 2), i))
    end do
    
    ewrite(2, *) "eta eta   component"
    col_offset = 2 * un    
    do i = 1, pn
      call set(d_big_m, i + row_offset, row_m_ptr(e_big_m(3, 3), i) + col_offset, row_val(e_big_m(3, 3), i))
    end do
    
    ! Statify the dynamic CSR matrix to form the system matrix
    ewrite(2, *) "Statifying"
    big_m = dcsr2csr(d_big_m)
    call deallocate(d_big_m)
    
    ewrite_minmax(big_m%val)
    ewrite(2, *) "sum(big_m%val) = ", sum(big_m%val)
    
    ewrite(1, *) "Flattening complete"
  
  end subroutine flatten_big_m
  
  subroutine increment_fields()
    !!< Increment the system fields using the solution vector
  
    u%val(1)%ptr = u%val(1)%ptr + big_solution%val(:node_count(u))
    u%val(2)%ptr = u%val(2)%ptr + big_solution%val(node_count(u) + 1:2 * node_count(u))
    p%val = p%val + big_solution%val(2 * node_count(u) + 1:)
    
  end subroutine increment_fields
  
  subroutine deallocate_matrices()
    !!< Deallocate the term matrices
  
    if(include_advection) then
      call deallocate(n_m)
      call deallocate(n_rhs)
    end if
    if(lump_l) then
      call deallocate(lumped_l_m)
    else
      call deallocate(l_m)
    end if    
    call deallocate(ct_m)
    call deallocate(c_m)
  
  end subroutine deallocate_matrices
  
  subroutine finalise_timestep()
    !!< Deallocate the system matrix and RHS and solution vectors
  
    call deallocate(linear_big_m)
    call deallocate(big_m)
    call deallocate(big_rhs)
    call deallocate(big_solution)
  
  end subroutine finalise_timestep
  
  subroutine initialise_balanced_velocity()
    !!< Initialise the velocity field to be in discrete geostrophic balance with
    !!< the layer thickness field
    
    type(scalar_field) :: rhs, u_x, u_y
        
    u_x = extract_scalar_field(u, 2)
    u_y = extract_scalar_field(u, 1)
    call allocate(rhs, u%mesh, "RHS")
        
    call mult(rhs, block(c_m, 1, 1), p)
    call scale(rhs, -g)
    if(lump_l) then
      u_x%val = rhs%val / lumped_l_m%val
    else
      call set_solver_options(temp_solver_path, ksptype = "cg", pctype = "sor", rtol = 0.0, atol = epsilon(0.0), max_its = 10000)
      call petsc_solve(u_x, l_m, rhs, option_path = temp_solver_path)
    end if
    
    call mult(rhs, block(c_m, 2, 1), p)
    call scale(rhs, g)
    if(lump_l) then
      u_y%val = rhs%val / lumped_l_m%val
    else
      call petsc_solve(u_y, l_m, rhs, option_path = temp_solver_path)
      call delete_option(temp_solver_path)
    end if
    
    call deallocate(rhs)
    
  end subroutine initialise_balanced_velocity
  
  subroutine prescribed_adapt()
    !!< Run a prescribed adapt
  
    integer :: stat
  
    ! Pre adapt debug output
    call vtk_write_fields("pre_adapt_u", adapt_no, position = x, model = u%mesh, vfields = (/u/), stat = stat)
    if(stat /= 0) then
      ewrite(0, *) "WARNING: Error returned by vtk_write_fields: ", stat
    end if
    call vtk_write_fields("pre_adapt_p", adapt_no, position = x, model = p%mesh, sfields = (/p/), stat = stat) 
    if(stat /= 0) then
      ewrite(0, *) "WARNING: Error returned by vtk_write_fields: ", stat
    end if
     
    call deallocate_matrices()
    call finalise_timestep()
                      
    call adapt_state_prescribed(states, current_time)
     
    call auxilliary_state_allocation()
    call initialise_pointers()
    call initialise_matrices()
    call initialise_timestep()
    
    ! Post adapt diagnostics
    call calculate_diagnostic_variables(states, exclude_nonrecalculated = .true.)
    call calculate_diagnostic_variables_new(states, exclude_nonrecalculated = .true.)
    call write_diagnostics(states, current_time, dt, timestep)
    
    ! Post adapt debug output
    call vtk_write_fields("post_adapt_u", adapt_no, position = x, model = u%mesh, vfields = (/u/), stat = stat)
    if(stat /= 0) then
      ewrite(0, *) "WARNING: Error returned by vtk_write_fields: ", stat
    end if
    call vtk_write_fields("post_adapt_p", adapt_no, position = x, model = p%mesh, sfields = (/p/), stat = stat)
    if(stat /= 0) then
      ewrite(0, *) "WARNING: Error returned by vtk_write_fields: ", stat
    end if
    
    adapt_no = adapt_no + 1 
  
  end subroutine prescribed_adapt

  subroutine simple_shallow_water_check_state()
  
    if(continuity(p) /= 0) then
      FLExit("LayerThickness must be on a continuous mesh")
    end if
    
    select case(continuity(u))
      case(0)
        if(.not. have_option("/material_phase::Fluid/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin")) then
          FLExit("Continuous Velocity mesh selected, but continuous_galerkin spatial_discretisation not selected")
        end if
      case(-1)
        if(.not. have_option("/material_phase::Fluid/vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin")) then
          FLExit("Discontinuous Velocity mesh selected, but discontinuous_galerkin spatial_discretisation not selected")
        end if
      case default
        ewrite(-1, *) "For mesh continuity: ", continuity(u)
        FLAbort("Unrecognised mesh continuity")
    end select
    
    call check_diagnostic_dependencies(states)
    
  end subroutine simple_shallow_water_check_state

  subroutine simple_shallow_water_check_options()
    !!< Check simple shallow water options
      
    character(len = OPTION_PATH_LEN) :: model
      
    if(isparallel()) then
      FLExit("Serial job required")
    end if
    
    if(option_count("/material_phase::Fluid/vector_field::Velocity/prognostic/initial_condition/balanced") > 0 .and. &
      & option_count("/material_phase::Fluid/vector_field::Velocity/prognostic/initial_condition::WholeMesh/balanced") == 0) then
      FLExit("Balanced Velocity initialisation must be for WholeMesh")
    end if
    
    call get_option("/model", model, default = "unknown")
    if(model == "simple_shallow_water") then
      ! sswml input
      call check_options()
    else
      ! swml input
      
      if(.not. have_option("/material_phase::Fluid/vector_field::Velocity/prognostic/")) then
        FLExit("Prognostic Velocity required")
      end if   
      if(.not. have_option("/material_phase::Fluid/scalar_field::LayerThickness/prognostic/spatial_discretisation/continuous_galerkin/advection_terms/exclude_advection_terms")) then
        ewrite(-1, *) "Disable advection for LayerThickness"
        FLExit("Advection not supported for Layer Thickness")
      end if
      
      if(have_option("/debug/check_inverse_coriolis_matrix")) then
        ewrite(0, *) "Warning: No inverse Coriolis matrix check will be performed"
      end if
      if(have_option("/debug/check_wave_matrix")) then
        ewrite(0, *) "Warning: No wave matrix check will be performed (no such matrix will be assembled!)"
      end if
      if(have_option("/debug/check_solution")) then
        ewrite(0, *) "Warning: No solution check will be performed"
      end if
    end if
  
  end subroutine simple_shallow_water_check_options

end program simple_shallow_water
