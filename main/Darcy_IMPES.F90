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
#include "confdefs.h"

program Darcy_IMPES
   
   ! Intrinsic Fortran modules
   use ISO_FORTRAN_ENV, only: OUTPUT_UNIT, ERROR_UNIT
   
   ! Fluidity modules
   use AuxilaryOptions
   use MeshDiagnostics
   use signal_vars
   use spud
   use equation_of_state
   use timers
   use adapt_state_module
   use adapt_state_prescribed_module
   use FLDebug
   use sparse_tools
   use elements
   use fields
   use field_options
   use reserve_state_module
   use vtk_interfaces
   use Diagnostic_variables
   use diagnostic_fields_new, only : &
     & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
     & check_diagnostic_dependencies
   use diagnostic_fields_wrapper
   use diagnostic_children
   use field_equations_cv, only: solve_field_eqn_cv, initialise_advection_convergence, coupled_cv_field_eqn
   use vertical_extrapolation_module
   use qmesh_module
   use synthetic_bc
   use goals
   use adaptive_timestepping
   use conformity_measurement
   use adjacency_lists
   use parallel_tools
   use write_triangle
   use timeloop_utilities
   use free_surface_module
   use field_priority_lists
   use boundary_conditions
   use discrete_properties_module
   use halos
   use memory_diagnostics
   use global_parameters, only: current_time, &
                                dt, &
                                timestep, &
                                OPTION_PATH_LEN, &
                                simulation_start_time, &
                                simulation_start_cpu_time, &
                                simulation_start_wall_time, &
                                topology_mesh_name
   use eventcounter
   use reduced_model_runtime
   use detector_parallel, only: sync_detector_coordinates, deallocate_detector_list_array
#ifdef HAVE_ZOLTAN
   use zoltan
#endif
   use tictoc
   use signals
   use write_state_module
   use populate_state_module
   use boundary_conditions_from_options
   use checkpoint
   use sparsity_patterns_meshes
   use cv_shape_functions
   use cv_faces
   use cv_options
   use FEFields
   
   ! *** Use Darcy IMPES module ***
   use darcy_impes_assemble_module
      
   implicit none

#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif

   interface
      subroutine set_global_debug_level(n)
         integer, intent(in) :: n
      end subroutine set_global_debug_level

      subroutine mpi_init(ierr)
         integer, intent(out) :: ierr
      end subroutine mpi_init

      subroutine mpi_finalize(ierr)
         integer, intent(out) :: ierr
      end subroutine mpi_finalize

      subroutine python_init
      end subroutine python_init

      subroutine petscinitialize(s, i)
         character(len=*), intent(in) :: s
         integer, intent(out) :: i
      end subroutine petscinitialize
   end interface

   interface
      subroutine check_options
      end subroutine check_options
   end interface

   character(len = OPTION_PATH_LEN) :: filename

   real :: finish_time

   type(state_type), dimension(:), pointer :: state => null()    
   
   type(tensor_field) :: metric_tensor

   integer :: dump_no = 0
   
   character(len=OPTION_PATH_LEN) :: option_buffer

   integer :: i, ierr

   logical :: not_to_move_det_yet = .false.
   
   logical :: output_log_file = .false.
   
   ! Non linear iteration variables
   integer :: its, nonlinear_iterations, nonlinear_iterations_adapt
   real :: change, chaold, nonlinear_iteration_tolerance
   
   ! *** Data associated with the Darcy IMPES solver ***  
   type(darcy_impes_type) :: di
   integer :: darcy_debug_log_unit, darcy_debug_err_unit   
   
#ifdef HAVE_ZOLTAN
   real(zoltan_float) :: ver
   integer(zoltan_int) :: ierrz
#endif

   call tictoc_reset()     

   call tic(TICTOC_ID_SIMULATION)

#ifdef HAVE_MPI
   call mpi_init(ierr)
   assert(ierr == MPI_SUCCESS)
#endif

#ifdef HAVE_PETSC
   call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

   call python_init()
    
   call read_command_line_load_options_set_simulation_name(darcy_debug_log_unit, darcy_debug_err_unit)
   
   ewrite(1,*) "***************"    
   ewrite(1,*) "* Darcy IMPES *"
   ewrite(1,*) "***************"

#ifdef HAVE_ZOLTAN
   ierrz = Zoltan_Initialize(ver)  
   assert(ierrz == ZOLTAN_OK)
#endif

   call initialise_signals()

   call get_option("/simulation_name",filename)

   call set_simulation_start_times()
    
   call initialise_walltime
    
   timestep = 0

   call initialise_qmesh
    
   call initialise_write_state
   
   call populate_state(state)
   
   call check_diagnostic_dependencies(state)

   call allocate_and_insert_auxilliary_fields(state)   

   default_stat%zoltan_drive_call=.false.

   ! set the nonlinear timestepping options, needs to be before the adapt at first timestep
   call get_option("/timestepping/nonlinear_iterations", &
                  & nonlinear_iterations, &
                  & default = 1)
   
   call get_option("/timestepping/nonlinear_iterations/tolerance", &
                  & nonlinear_iteration_tolerance, &
                  & default = 0.0)

   call get_option("/timestepping/timestep", dt)   
   
   if(have_option("/timestepping/adaptive_timestep/at_first_timestep")) then
      call calc_cflnumber_field_based_dt(state, dt, force_calculation = .true.)
      call set_option("/timestepping/timestep", dt)
   end if

   call get_option("/timestepping/current_time", current_time)
   call get_option("/timestepping/finish_time", finish_time)
 
   call get_option("/timestepping/nonlinear_iterations", &
                  & nonlinear_iterations, &
                  & default = 1)
   
   call get_option("/timestepping/nonlinear_iterations/tolerance", &
                  & nonlinear_iteration_tolerance, &
                  & default = 0.0)

   ! *** Initialise data used in IMPES solver *** 
   call darcy_impes_initialise(di, state, dt)
          
   if(have_option("/mesh_adaptivity/hr_adaptivity/adapt_at_first_timestep")) then

      if(have_option("/timestepping/nonlinear_iterations/nonlinear_iterations_at_adapt")) then
         call get_option('/timestepping/nonlinear_iterations/nonlinear_iterations_at_adapt',nonlinear_iterations_adapt)
         nonlinear_iterations = nonlinear_iterations_adapt
      end if

      call adapt_state_first_timestep(state)
      
      ! *** Update Darcy IMPES post spatial adapt ***
      call darcy_impes_update_post_spatial_adapt(di)
            
      ! Ensure that checkpoints do not adapt at first timestep.
      call delete_option("/mesh_adaptivity/hr_adaptivity/adapt_at_first_timestep")
   end if
   
   ! Initial diagnostics output via ewrite
   call run_diagnostics(state)

   if (have_option("/mesh_adaptivity/hr_adaptivity")) then
      call allocate(metric_tensor, extract_mesh(state(1), topology_mesh_name), "ErrorMetric")
   end if

   ! Determine the output format.
   call get_option('/io/dump_format', option_buffer)
   if(trim(option_buffer) /= "vtk") then
      ewrite(-1,*) "You must specify a dump format and it must be vtk."
      FLExit("Rejig your: /io/dump_format")
   end if
   
   call initialise_diagnostics(filename, state)
    
   ! Checkpoint at start
   if(do_checkpoint_simulation(dump_no)) call checkpoint_simulation(state, cp_no = dump_no)
   ! Dump at start
   if( &
       ! if this is not a zero timestep simulation (otherwise, there would
       ! be two identical dump files)
       & current_time < finish_time &
       ! unless explicitly disabled
       & .and. .not. have_option("/io/disable_dump_at_start") &
       & ) then
       call write_state(dump_no, state)
   end if

   call initialise_convergence(filename, state)
   call initialise_steady_state(filename, state)

   if(have_option("/io/stat/output_at_start")) call write_diagnostics(state, current_time, dt, timestep, not_to_move_det_yet=.true.)

   not_to_move_det_yet=.false.
         
   ! ******************************
   ! *** Start of timestep loop ***
   ! ******************************

   timestep_loop: do
      timestep = timestep + 1

      ewrite(1, *) "********************"
      ewrite(1, *) "*** NEW TIMESTEP ***"
      ewrite(1, *) "********************"
      ewrite(1, *) "Current simulation time: ", current_time
      ewrite(1, *) "Timestep number: ", timestep
      ewrite(1, *) "Timestep size (dt): ", dt
      if(.not. allfequals(dt)) then
         ewrite(-1, *) "Timestep size (dt): ", dt
         FLAbort("The timestep is not global across all processes!")
      end if

      if(simulation_completed(current_time, timestep)) exit timestep_loop

      if( &
                               ! Do not dump at the start of the simulation (this is handled by write_state call earlier)
           & current_time > simulation_start_time &
                               ! Do not dump at the end of the simulation (this is handled by later write_state call)
           & .and. current_time < finish_time &
                               ! Test write_state conditions
           & .and. do_write_state(current_time, timestep) &
           & ) then

         ! Intermediate dumps
         if(do_checkpoint_simulation(dump_no)) then
            call checkpoint_simulation(state, cp_no = dump_no)
         end if
         call write_state(dump_no, state)
      end if

      call copy_to_stored_values(state,"Old")
      
      ! this may already have been done in populate_state, but now
      ! we evaluate at the correct "shifted" time level:
      call set_boundary_conditions_values(state, shift_time=.true.)
      
      ! evaluate prescribed fields at time = current_time+dt
      call set_prescribed_field_values(state, exclude_interpolated=.true., &
           exclude_nonreprescribed=.true., time=current_time+dt)

      nonlinear_iteration_loop: do its = 1,nonlinear_iterations

         ewrite(1,*)'###################'
         ewrite(1,*)'Start of another nonlinear iteration; its,nonlinear_iterations=',its,nonlinear_iterations
         ewrite(1,*)'###################'         
         
         call copy_to_stored_values(state, "Iterated")
      
         ! *** Solve the Darcy equations using IMPES ***
         call darcy_impes_assemble_and_solve(di)

         if(nonlinear_iterations > 1) then
            
            ! Check for convergence between non linear iteration loops
            call test_and_write_convergence(state, current_time + dt, dt, its, change)
            
            if(its == 1) chaold = change

            if (have_option("/timestepping/nonlinear_iterations/tolerance")) then
              
              ewrite(2, *) "Nonlinear iteration change = ", change
              ewrite(2, *) "Nonlinear iteration tolerance = ", nonlinear_iteration_tolerance

              if(change < abs(nonlinear_iteration_tolerance)) then
                 ewrite(1, *) "Nonlinear iteration tolerance has been reached"
                 ewrite(1, "(a,i0,a)") "Exiting nonlinear iteration loop after ", its, " iterations"
                 exit nonlinear_iteration_loop
              endif
            
            end if
         
         end if

      end do nonlinear_iteration_loop

      ! Reset the number of nonlinear iterations in case it was overwritten by nonlinear_iterations_adapt
      call get_option('/timestepping/nonlinear_iterations', nonlinear_iterations, default=1)

      if(have_option("/timestepping/nonlinear_iterations/terminate_if_not_converged")) then
         if(its >= nonlinear_iterations .and. change >= abs(nonlinear_iteration_tolerance)) then
            ewrite(0, *) "Nonlinear iteration tolerance not reached - termininating"
            exit timestep_loop
         end if
      end if
            
      current_time = current_time + DT

      ! if strong bc or weak that overwrite then enforce the bc on the fields
      call set_dirichlet_consistent(state)
      
      ! calculate and write diagnostics before the timestep gets changed
      call calculate_diagnostic_variables(state, exclude_nonrecalculated=.true.)
      call calculate_diagnostic_variables_new(state, exclude_nonrecalculated = .true.)
       
      ! Call the modern and significantly less satanic version of study
      call write_diagnostics(state, current_time, dt, timestep)
      
      if(have_option("/timestepping/adaptive_timestep")) call calc_cflnumber_field_based_dt(state, dt)
      
      call set_option("/timestepping/timestep", dt)
      call set_option("/timestepping/current_time", current_time)            
         
      ! ******************
      ! *** Mesh adapt ***
      ! ******************

      adapt_if: if(have_option("/mesh_adaptivity/hr_adaptivity")) then

         do_adapt_if: if(do_adapt_mesh(current_time, timestep)) then

            call qmesh(state, metric_tensor)
            if(have_option("/io/stat/output_before_adapts")) call write_diagnostics(state, current_time, dt, timestep, not_to_move_det_yet=.true.)
            call run_diagnostics(state)

            call adapt_state(state, metric_tensor)
            
            ! Overwrite the number of nonlinear iterations if the option is switched on
            if(have_option("/timestepping/nonlinear_iterations/nonlinear_iterations_at_adapt")) then
               call get_option('/timestepping/nonlinear_iterations/nonlinear_iterations_at_adapt',nonlinear_iterations_adapt)
               nonlinear_iterations = nonlinear_iterations_adapt
            end if

            ! re-allocate the adaptivity metric
            if(have_option("/mesh_adaptivity/hr_adaptivity")) then
               call allocate(metric_tensor, extract_mesh(state(1), topology_mesh_name), "ErrorMetric")
            end if
            
            ! *** Update Darcy IMPES post spatial adapt ***
            call darcy_impes_update_post_spatial_adapt(di)

            if (have_option("/mesh_adaptivity/hr_adaptivity/adaptive_timestep_at_adapt")) then
               if (have_option("/timestepping/adaptive_timestep/minimum_timestep")) then
                  call get_option("/timestepping/adaptive_timestep/minimum_timestep", dt)
                  call set_option("/timestepping/timestep", dt)
               else
                  ewrite(-1,*) "Warning: you have adaptive timestep adjustment after &&
                                && adapt, but have not set a minimum timestep"
               end if
            else
               ! Timestep adapt
               if(have_option("/timestepping/adaptive_timestep")) then
                  call calc_cflnumber_field_based_dt(state, dt, force_calculation = .true.)
                  call set_option("/timestepping/timestep", dt)
               end if
            end if
                        
            if(have_option("/io/stat/output_after_adapts")) call write_diagnostics(state, current_time, dt, timestep, not_to_move_det_yet=.true.)
            
            call run_diagnostics(state)
                        
         end if do_adapt_if

      end if adapt_if

      ! *** Set the darcy impes time step ***
      di%dt = dt

   end do timestep_loop
   
   ! ****************************
   ! *** END OF TIMESTEP LOOP ***
   ! ****************************
      
   ! Checkpoint at end, if enabled
   if(have_option("/io/checkpointing/checkpoint_at_end")) then
      call checkpoint_simulation(state, cp_no = dump_no)
   end if
   ! Dump at end, unless explicitly disabled
   if(.not. have_option("/io/disable_dump_at_end")) then
      call write_state(dump_no, state)
   end if
   
   ! closing .stat, .convergence and .detector files
   call close_diagnostic_files()

   ! deallocate the array of all detector lists
   call deallocate_detector_list_array()

   ewrite(1, *) "Printing references before final deallocation"
   call print_references(1)

   ! Deallocate the metric tensor
   if(have_option("/mesh_adaptivity/hr_adaptivity")) call deallocate(metric_tensor)
   
   ! *** Finalise darcy impes variables ***
   call darcy_impes_finalise(di)
    
   ! Deallocate state
   do i = 1, size(state)
      call deallocate(state(i))
   end do
    
   ! Deallocate the reserve state
   call deallocate_reserve_state()

   deallocate(state)

   ! Clean up registered diagnostics
   call destroy_registered_diagnostics 

   ! Delete the transform_elements cache.
   call deallocate_transform_cache

   ewrite(1, *) "Printing references after final deallocation"
   call print_references(0)

#ifdef HAVE_MEMORY_STATS
    call print_current_memory_stats(0)
#endif
  
   call toc(TICTOC_ID_SIMULATION)
   call tictoc_report(2, TICTOC_ID_SIMULATION)

   ewrite(1,*) "************************"    
   ewrite(1,*) "* Finished Darcy IMPES *"
   ewrite(1,*) "************************"

   if(output_log_file) close(darcy_debug_log_unit)
   if(output_log_file) close(darcy_debug_err_unit)

#ifdef HAVE_MPI
    call mpi_finalize(ierr)
    assert(ierr == MPI_SUCCESS)
#endif

   if(SIG_INT) then
     FLExit("Interrupt signal received")
   end if

contains

! ----------------------------------------------------------------------------

   subroutine darcy_impes_initialise(di, &
                                     state, &
                                     dt)
      
      !!< Initialise the Darcy IMPES type from options and state
      
      type(darcy_impes_type),                       intent(inout) :: di
      type(state_type),       dimension(:), target, intent(in)    :: state
      real,                                         intent(in)    :: dt
      
      ! Local variables
      integer :: p
      
      di%dt = dt
      
      di%state => state
      
      di%pressure                  => extract_scalar_field(di%state(1), "Pressure")
      di%old_pressure              => extract_scalar_field(di%state(1), "OldPressure")
      di%gradient_pressure         => extract_vector_field(di%state(1), "GradientPressure")
      di%old_gradient_pressure     => extract_vector_field(di%state(1), "OldGradientPressure")
      di%porosity                  => extract_scalar_field(di%state(1), "Porosity")
      di%old_porosity              => extract_scalar_field(di%state(1), "OldPorosity")
      di%absolute_permeability     => extract_scalar_field(di%state(1), "AbsolutePermeability")
      di%positions                 => extract_vector_field(di%state(1), "Coordinate")
      di%total_darcy_velocity      => extract_vector_field(di%state(1), "TotalDarcyVelocity")
      di%sum_saturation            => extract_scalar_field(di%state(1), "SumSaturation")
      di%old_sum_saturation        => extract_scalar_field(di%state(1), "OldSumSaturation")
      di%div_total_darcy_velocity  => extract_scalar_field(di%state(1), "DivergenceTotalDarcyVelocity")

      ! Form the positions on the pressure mesh
      di%positions_pressure_mesh = get_coordinate_field(di%state(1), di%pressure%mesh)

      ! Determine the pressure matrix sparsity
      di%sparsity => get_csr_sparsity_firstorder(di%state(1), di%pressure%mesh, di%pressure%mesh)

      ! Allocate the pressure matrix and lhs and rhs to use for saturations
      call allocate(di%pressure_matrix, di%sparsity)

      call allocate(di%lhs, di%pressure%mesh)
      call allocate(di%rhs, di%pressure%mesh)
      call allocate(di%rhs_adv, di%pressure%mesh)
      call allocate(di%rhs_time, di%pressure%mesh)
      call allocate(di%old_saturation_subcycle, di%pressure%mesh)
      call allocate(di%cv_mass_pressure_mesh_with_porosity, di%pressure%mesh)
      call allocate(di%cv_mass_pressure_mesh_with_old_porosity, di%pressure%mesh)

      ! Calculate the latest CV mass on the pressure mesh with porosity
      call compute_cv_mass(di%positions, di%cv_mass_pressure_mesh_with_porosity, di%porosity)      

      ! deduce the number of phase
      di%number_phase = size(di%state)

      ! Pull phase dependent fields from state

      allocate(di%saturation(di%number_phase))
      allocate(di%old_saturation(di%number_phase))
      allocate(di%relative_permeability(di%number_phase))
      allocate(di%viscosity(di%number_phase))      
      allocate(di%inter_velocity_porosity(di%number_phase))
      allocate(di%old_inter_velocity_porosity(di%number_phase))
      allocate(di%cfl(di%number_phase))
      allocate(di%old_cfl(di%number_phase))
      allocate(di%darcy_velocity(di%number_phase))
      allocate(di%fractional_flow(di%number_phase))

      do p = 1,di%number_phase

         di%saturation(p)%ptr                  => extract_scalar_field(di%state(p), "Saturation")
         di%old_saturation(p)%ptr              => extract_scalar_field(di%state(p), "OldSaturation")
         di%relative_permeability(p)%ptr       => extract_scalar_field(di%state(p), "RelativePermeability")
         di%viscosity(p)%ptr                   => extract_scalar_field(di%state(p), "Viscosity")
         di%inter_velocity_porosity(p)%ptr     => extract_vector_field(di%state(p), "InterstitialVelocityPorosity")
         di%old_inter_velocity_porosity(p)%ptr => extract_vector_field(di%state(p), "OldInterstitialVelocityPorosity")
         di%cfl(p)%ptr                         => extract_scalar_field(di%state(p), "InterstitialVelocityCFL")
         di%old_cfl(p)%ptr                     => extract_scalar_field(di%state(p), "OldInterstitialVelocityCFL")
         di%darcy_velocity(p)%ptr              => extract_vector_field(di%state(p), "DarcyVelocity")
         di%fractional_flow(p)%ptr             => extract_vector_field(di%state(p), "FractionalFlow")

      end do

      ! Allocate field used for the subcycle time step cfl
      call allocate(di%cfl_subcycle, di%cfl(1)%ptr%mesh)

      di%phase_one_saturation_diagnostic = have_option(trim(di%saturation(1)%ptr%option_path)//'/diagnostic')

      ! Determine the inverse cv mass matrix of cfl mesh
      call allocate(di%inverse_cv_mass_cfl_mesh, &
                    di%cfl(1)%ptr%mesh)

      call compute_cv_mass(di%positions, di%inverse_cv_mass_cfl_mesh)

      call invert(di%inverse_cv_mass_cfl_mesh)   

      ! Determine the inverse cv mass matrix of pressure mesh
      call allocate(di%inverse_cv_mass_pressure_mesh, &
                    di%pressure%mesh)

      call compute_cv_mass(di%positions, di%inverse_cv_mass_pressure_mesh)

      call invert(di%inverse_cv_mass_pressure_mesh)   

      ! Get the saturation cv options from the first phase
      di%saturation_cv_options = get_cv_options(di%saturation(1)%ptr%option_path, &
                                               &di%saturation(1)%ptr%mesh%shape%numbering%family, &
                                               &mesh_dim(di%saturation(1)%ptr))

      ! Determine the max courant number per saturation subcycle
      call get_option(trim(complete_field_path(di%saturation(1)%ptr%option_path))//'/temporal_discretisation/control_volumes/maximum_courant_number_per_subcycle', &
                      di%saturation_max_courant_per_subcycle, &
                      default = 0.25)

      ! Determine the CV surface degree to use when integrating functions across them
      call get_option("/geometry/quadrature/controlvolume_surface_degree", &
                      di%quaddegree, default = 1)

      ! Determine the CV faces information for the pressure mesh
      ! assuming all elements are the same type
      di%cvfaces = find_cv_faces(vertices   = ele_vertices(di%pressure, 1), &
                                 dimension  = mesh_dim(di%pressure), &
                                 polydegree = di%pressure%mesh%shape%degree, &
                                 quaddegree = di%quaddegree)

      ! Generate the CV shape function that contains the derivatives with respect
      ! to the parent elements canonical coordinates evaluated at the 
      ! control volume faces for the positions mesh, assuming all elements 
      ! are the same type.      
      di%x_cvshape_full = make_cv_element_shape(di%cvfaces, di%positions%mesh%shape, &
                                                type = ELEMENT_CONTROLVOLUME_SURFACE_BODYDERIVATIVES)

      ! Generate the CV shape function that contains the derivatives with respect
      ! to the parent elements canonical coordinates evaluated at the 
      ! control volume faces for the pressure mesh, assuming all elements 
      ! are the same type.      
      di%p_cvshape_full = make_cv_element_shape(di%cvfaces, di%pressure%mesh%shape, &
                                                type = ELEMENT_CONTROLVOLUME_SURFACE_BODYDERIVATIVES)

      ! Generate the CV shape function with reduced number of derivatives 
      ! evaluated at the control volume faces for the positions mesh, 
      ! assuming all elements are the same type. Reduced refers to 
      ! the CV suface being 1 dimension lower than the problem dimension. 
      ! For example in 3d the CV surface is 2d, so this shape function 
      ! contains the 2d derivatives across the CV face.
      di%x_cvshape = make_cv_element_shape(di%cvfaces, di%positions%mesh%shape)

      ! Generate the CV shape function with reduced number of derivatives 
      ! evaluated at the control volume faces for the pressure mesh, 
      ! assuming all elements are the same type.     
      di%p_cvshape = make_cv_element_shape(di%cvfaces, di%pressure%mesh%shape) 

      ! Generate the CV shape function with reduced number of derivatives
      ! evaluated at the control volume faces located on the domain boundary
      ! for the positions mesh assuming all boundary elements are the same type.
      di%x_cvbdyshape = make_cvbdy_element_shape(di%cvfaces, di%positions%mesh%faces%shape)

      ! If the first phase saturation is diagnostic then calculate it
      if (di%phase_one_saturation_diagnostic) call darcy_impes_calculate_phase_one_saturation_diagnostic(di)

      ! Calculate the relative permeabilities - and other generic diagnostic fields
      ! (NOTE this is required before the darcy impes specific diagnostic fields below are calculated)
      call calculate_diagnostic_variables(di%state)
      call calculate_diagnostic_variables_new(di%state)

      ! Copy field values to Old - required before below
      call copy_to_stored_values(di%state,"Old")

      ! Initialise the Darcy IMPES specific diagnostic fields (gradient pressure, inter_velocity_porosity ... etc)
      call darcy_impes_calculate_gradient_pressure_etc(di)

      ! Calculate the sum of the saturations
      call darcy_impes_calculate_sum_saturation(di)

      ! Copy ALL latest fields to Old and Iterated
      call copy_to_stored_values(di%state,"Old")
      call copy_to_stored_values(di%state,"Iterated")
      call relax_to_nonlinear(di%state)
   
   end subroutine darcy_impes_initialise

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_finalise(di)
      
      !!< Finalise (ie deallocate) the Darcy IMPES data
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! Deallocate data types that are NOT pointers to data in di%state
      
      call deallocate(di%positions_pressure_mesh)
      call deallocate(di%pressure_matrix)
      call deallocate(di%lhs)
      call deallocate(di%rhs)
      call deallocate(di%rhs_adv)
      call deallocate(di%rhs_time)
      call deallocate(di%old_saturation_subcycle)      
      call deallocate(di%inverse_cv_mass_cfl_mesh)
      call deallocate(di%inverse_cv_mass_pressure_mesh)
      call deallocate(di%cv_mass_pressure_mesh_with_porosity)
      call deallocate(di%cv_mass_pressure_mesh_with_old_porosity)
      call deallocate(di%cfl_subcycle)
      deallocate(di%saturation)
      deallocate(di%old_saturation)
      deallocate(di%relative_permeability)
      deallocate(di%viscosity)
      deallocate(di%inter_velocity_porosity)
      deallocate(di%old_inter_velocity_porosity)      
      deallocate(di%cfl)
      deallocate(di%old_cfl)
      deallocate(di%darcy_velocity) 
      deallocate(di%fractional_flow) 
      call deallocate(di%cvfaces)
      call deallocate(di%x_cvshape_full)
      call deallocate(di%p_cvshape_full)
      call deallocate(di%x_cvshape)
      call deallocate(di%p_cvshape)
      call deallocate(di%x_cvbdyshape)
            
   end subroutine darcy_impes_finalise

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_update_post_spatial_adapt(di)
      
      !!< Update the Darcy IMPES data post spatial adapt
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! Auxilliary fields.
      call allocate_and_insert_auxilliary_fields(di%state)

      ! If the first phase saturation is diagnostic then calculate it
      if (di%phase_one_saturation_diagnostic) call darcy_impes_calculate_phase_one_saturation_diagnostic(di)

      ! Calculate the relative permeabilities - and other generic diagnostic fields
      ! (NOTE this is required before the darcy impes specific diagnostic fields below are calculated)
      call calculate_diagnostic_variables(di%state)
      call calculate_diagnostic_variables_new(di%state)

      ! Copy field values to Old - required before below
      call copy_to_stored_values(di%state,"Old")

      ! Initialise the Darcy IMPES specific diagnostic fields (gradient pressure, inter_velocity_porosity ... etc)
      call darcy_impes_calculate_gradient_pressure_etc(di)

      ! Calculate the sum of the saturations
      call darcy_impes_calculate_sum_saturation(di)

      ! Copy ALL latest fields to Old and Iterated
      call copy_to_stored_values(di%state,"Old")
      call copy_to_stored_values(di%state,"Iterated")
      call relax_to_nonlinear(di%state)

      ! Re-allocate IMPES data
      call deallocate(di%positions_pressure_mesh)
      call deallocate(di%pressure_matrix)
      call deallocate(di%lhs)
      call deallocate(di%rhs)
      call deallocate(di%rhs_adv)
      call deallocate(di%rhs_time)
      call deallocate(di%old_saturation_subcycle)      
      call deallocate(di%inverse_cv_mass_cfl_mesh)
      call deallocate(di%inverse_cv_mass_pressure_mesh)
      call deallocate(di%cv_mass_pressure_mesh_with_porosity)
      call deallocate(di%cv_mass_pressure_mesh_with_old_porosity)
      call deallocate(di%cfl_subcycle)

      call compute_cv_mass(di%positions, di%cv_mass_pressure_mesh_with_porosity, di%porosity)      

      di%positions_pressure_mesh = get_coordinate_field(di%state(1), di%pressure%mesh)

      di%sparsity => get_csr_sparsity_firstorder(di%state(1), di%pressure%mesh, di%pressure%mesh)

      call allocate(di%pressure_matrix, di%sparsity)

      call allocate(di%lhs, di%pressure%mesh)   
      call allocate(di%rhs, di%pressure%mesh)
      call allocate(di%rhs_adv, di%pressure%mesh)
      call allocate(di%rhs_time, di%pressure%mesh)
      call allocate(di%old_saturation_subcycle, di%pressure%mesh)
      call allocate(di%cv_mass_pressure_mesh_with_porosity, di%pressure%mesh)
      call allocate(di%cv_mass_pressure_mesh_with_old_porosity, di%pressure%mesh)

      call allocate(di%cfl_subcycle, di%cfl(1)%ptr%mesh)

      call allocate(di%inverse_cv_mass_cfl_mesh, &
                    di%cfl(1)%ptr%mesh)

      call compute_cv_mass(di%positions, di%inverse_cv_mass_cfl_mesh)

      call invert(di%inverse_cv_mass_cfl_mesh)   

      call allocate(di%inverse_cv_mass_pressure_mesh, &
                    di%pressure%mesh)

      call compute_cv_mass(di%positions, di%inverse_cv_mass_pressure_mesh)

      call invert(di%inverse_cv_mass_pressure_mesh)      
      
   end subroutine darcy_impes_update_post_spatial_adapt

! --------------------------------------------------------------------------------

   subroutine set_simulation_start_times()
      !!< Set the simulation start times

      call get_option("/timestepping/current_time", simulation_start_time)

      call cpu_time(simulation_start_cpu_time)
      call allmax(simulation_start_cpu_time)

      simulation_start_wall_time = wall_time()
      call allmax(simulation_start_wall_time)

   end subroutine set_simulation_start_times

! --------------------------------------------------------------------------------

   subroutine read_command_line_load_options_set_simulation_name(darcy_debug_log_unit, darcy_debug_err_unit)

      !!< Read the input filename and debug level from the command line
      !!< then read the options files. Check for a simulation name in the
      !!< options file and if not found add one using the filename.diml
      !!< There is also a check on the command line for the debug verbosity
      !!< and if this should be written to a log file.
      
      integer, intent(inout) :: darcy_debug_log_unit, darcy_debug_err_unit
      
      ! local variables
      character(len=1024) :: argument, output_filename
      integer :: status, argn, debug_level, ierror, rank, nprocs
      
      ! Initialise the global debug level
      call set_global_debug_level(0)
      
      ! loop each argument found on the command line
      argn=1
      read_loop: do 

         call get_command_argument(argn, value=argument, status=status)
         argn=argn+1

         if (status/=0) call usage()

         have_verbose: if (argument(1:2)=='-v') then

            read(argument(3:3), '(i1)', err=666) debug_level
            
            call set_global_debug_level(debug_level)

            ! Go back to pick up the command line.
            cycle read_loop
         
         end if have_verbose
         
         have_log: if (argument(1:2) == '-l') then 
            
            output_log_file = .true.
           
            ! Go back to pick up the command line.
            cycle read_loop
         
         end if have_log  

         ! Found everything we need with last argument being the simulation input filename
         exit read_loop
         
      end do read_loop
      
      ! Open log file if required
      open_log: if (output_log_file) then
           
         darcy_debug_log_unit = OUTPUT_UNIT
#ifdef HAVE_MPI         
         rank = getrank()
         
         nprocs = getnprocs()
         
         if (nprocs > 1) then
            output_filename = argument(1:len_trim(argument)-4)//"log_"//int2str(rank)
         else
            output_filename = argument(1:len_trim(argument)-4)//"log"
         end if
#else
         output_filename = argument(1:len_trim(argument)-4)//"log"
#endif
         open(unit   = darcy_debug_log_unit, &
              file   = trim(output_filename), &
              status = 'replace', &
              form   = 'FORMATTED', &
              action = 'write', &
              iostat = ierror)
         
         if (ierror /= 0) then
            ewrite(0,*) 'iostat from opening output log file :',ierror
            FLExit('Error opening output log')
         end if
           
         darcy_debug_err_unit = ERROR_UNIT
#ifdef HAVE_MPI         
         rank = getrank()
         
         nprocs = getnprocs()
         
         if (nprocs > 1) then
            output_filename = argument(1:len_trim(argument)-4)//"err_"//int2str(rank)
         else
            output_filename = argument(1:len_trim(argument)-4)//"err"
         end if
#else
         output_filename = argument(1:len_trim(argument)-4)//"err"
#endif         
         open(unit   = darcy_debug_err_unit, &
              file   = trim(output_filename), &
              status = 'replace', &
              form   = 'FORMATTED', &
              action = 'write', &
              iostat = ierror)
         
         if (ierror /= 0) then
            ewrite(0,*) 'iostat from opening output err file :',ierror
            FLExit('Error opening output err')
         end if
      
      end if open_log
        
      ! Load the SPUD options
      call load_options(argument)
      
      ! If there is no simulation name in the options then add one
      ! using the filename which is currently what argument is.
      add_sim_name: if (.not. have_option('/simulation_name')) then
         
         ewrite(3,*) 'Options dictionary has no simulation_name so adding one using filename'
         
         ! check filename is long enough to have appended .diml
         if (len_trim(argument) < 6) call usage()
         
         ! check the last five characters of filename are .diml
         if (argument(len_trim(argument)-4:) /= '.diml' .and. argument(len_trim(argument)-4:) /= '.diml') call usage()
         
         call add_option('/simulation_name', status) 

         if (status /= SPUD_NEW_KEY_WARNING) then
            ewrite(-1,*) 'SPUD status', status
            FLAbort('Issue adding option for simulation name')
         end if
         
         call set_option('/simulation_name', argument(:len_trim(argument)-5), status)
         
         if (status /= SPUD_NO_ERROR) then
            ewrite(-1,*) 'SPUD status', status
            FLAbort('Issue setting option for simulation name')
         end if
         
      end if add_sim_name
      
      return

666   call usage()

   end subroutine read_command_line_load_options_set_simulation_name

! --------------------------------------------------------------------------------

   subroutine usage()
      write (0,*) 'Darcy IMPES usage: '
      write (0,*) ''
      write (0,*) '   darcy_impes [-vINT] [-l] <options_file>.diml'
      write (0,*) ''
      write (0,*) 'where: -vINT sets the verbosity of debugging'
      write(0,*) '         -l outputs to a log file <options_file>.log'
      FLExit('Incorrect command line usage')
   end subroutine usage

! --------------------------------------------------------------------------------

end program Darcy_IMPES
