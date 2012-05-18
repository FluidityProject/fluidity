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
   use timers
   use adapt_state_module
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
   use vertical_extrapolation_module
   use qmesh_module
   use synthetic_bc
   use adaptive_timestepping
   use conformity_measurement
   use adjacency_lists
   use parallel_tools
   use write_triangle
   use timeloop_utilities
   use boundary_conditions
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
   use state_module
   
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

   type(scalar_field), pointer :: time_field

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
   
   ewrite(1,*) "***********************"    
   ewrite(1,*) "* Fluidity-DarcyIMPES *"
   ewrite(1,*) "***********************"

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
   
   default_stat%zoltan_drive_call=.false.

   call get_option("/timestepping/timestep", dt)   
   
   ! Set non linear iteration data - required before first adapt
   call get_option("/timestepping/nonlinear_iterations", &
                  & nonlinear_iterations, &
                  & default = 1)
   
   call get_option("/timestepping/nonlinear_iterations/tolerance", &
                  & nonlinear_iteration_tolerance, &
                  & default = 0.0)
   
   call get_option("/timestepping/current_time", current_time)
   call get_option("/timestepping/finish_time", finish_time)

   ! *** Initialise data used in IMPES solver *** 
   call darcy_impes_initialise(di, state, dt, current_time)

   ! Set Time field
   time_field => extract_scalar_field(state(1), 'Time')
   call set(time_field, current_time)
   
   ! Adapt time at first time step - required before first adapt
   if(di%adaptive_dt_options%have .and. di%adaptive_dt_options%at_first_dt) then
      call darcy_impes_calculate_cflnumber_field_based_dt(di)
      call set_option("/timestepping/timestep", di%dt)
      dt = di%dt
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

      if(simulation_completed(current_time, timestep)) then
         
         ewrite(1,*) 'Simulation completed. Exiting timestep loop'
         
         exit timestep_loop
      
      end if

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
      
      ! *** Darcy IMPES copy non state data to old ***
      call darcy_impes_copy_to_old(di)
      
      ! this may already have been done in populate_state, but now
      ! we evaluate at the correct "shifted" time level:
      call set_boundary_conditions_values(state, shift_time=.true.)
      
      ! evaluate prescribed fields at time = current_time+dt
      call set_prescribed_field_values(state, exclude_interpolated=.true., &
           exclude_nonreprescribed=.true., time=current_time+dt)
      
      ! *** Darcy IMPES calculate latest gravity field - direction field could be defined by python function ***
      if (di%have_gravity) then
         call set(di%gravity, di%gravity_direction)
         call scale(di%gravity, di%gravity_magnitude)
         ewrite_minmax(di%gravity)
      end if
      
      ! *** Darcy IMPES set max number non linear iterations for this time step ***
      di%max_nonlinear_iter_this_timestep = nonlinear_iterations

      ! *** Darcy IMPES Calculate the relperm and density first face values ***
      call darcy_impes_calculate_relperm_den_first_face_values(di)
      
      nonlinear_iteration_loop: do its = 1,nonlinear_iterations

         ewrite(1,*)'###################'
         ewrite(1,*)'Start of another nonlinear iteration; its,nonlinear_iterations=',its,nonlinear_iterations
         ewrite(1,*)'###################'         
         
         call copy_to_stored_values(state, "Iterated")
         
         ! *** Set the Darcy IMPES nonlinear_iter
         di%nonlinear_iter = its
         
         ! *** Solve the Darcy equations using IMPES ***
         call darcy_impes_assemble_and_solve(di)

         ! calculate generic diagnostics - DO NOT ADD 'calculate_diagnostic_variables'
         call calculate_diagnostic_variables_new(di%state, exclude_nonrecalculated = .true.)
         
         ! *** Calculate the Darcy IMPES Velocity, Mobilities, Fractional flow and CFL fields
         !     needs to be after the python fields ***
         call darcy_impes_calculate_vel_mob_ff_and_cfl_fields(di)

         if(nonlinear_iterations > 1) then
            
            ! Check for convergence between non linear iteration loops
            call test_and_write_convergence(state, current_time + dt, dt, its, change)
            
            ewrite(1,*) 'Nonlinear iteration, change: ',its, change
            
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
         
         if (SIG_INT) exit nonlinear_iteration_loop

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

      ! *** Update DarcyIMPES current_time ***
      di%current_time = current_time
      
      ! Set Time field
      time_field => extract_scalar_field(state(1), 'Time')
      call set(time_field, current_time)
       
      ! Call the modern and significantly less satanic version of study
      call write_diagnostics(state, current_time, dt, timestep)
      
      ! *** Darcy IMPES adaptive time stepping choice ***
      if(di%adaptive_dt_options%have) call darcy_impes_calculate_cflnumber_field_based_dt(di)
      
      dt = di%dt
      call set_option("/timestepping/timestep", di%dt)
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
            call darcy_impes_update_post_spatial_adapt(di, &
                                                       state, &
                                                       dt, &
                                                       current_time)

            if (have_option("/mesh_adaptivity/hr_adaptivity/adaptive_timestep_at_adapt")) then
               if (have_option("/timestepping/adaptive_timestep/minimum_timestep")) then
                  call get_option("/timestepping/adaptive_timestep/minimum_timestep", dt)
                  call set_option("/timestepping/timestep", dt)
                  ! *** Set Darcy IMPES dt
                  di%dt = dt
               else
                  ewrite(-1,*) "Warning: you have adaptive timestep adjustment after &&
                                && adapt, but have not set a minimum timestep"
               end if
            else
               ! Timestep adapt
               if(di%adaptive_dt_options%have) then
                  call darcy_impes_calculate_cflnumber_field_based_dt(di)
                  call set_option("/timestepping/timestep", di%dt)
                  dt = di%dt
               end if
            end if
                        
            if(have_option("/io/stat/output_after_adapts")) call write_diagnostics(state, current_time, dt, timestep, not_to_move_det_yet=.true.)
            
            call run_diagnostics(state)
                        
         end if do_adapt_if

      end if adapt_if
      
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

   ewrite(1,*) "********************************"    
   ewrite(1,*) "* Finished Fluidity-DarcyIMPES *"
   ewrite(1,*) "********************************"

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
                                     dt, &
                                     current_time)
      
      !!< Initialise the Darcy IMPES type from options and state
      
      type(darcy_impes_type),                       intent(inout) :: di
      type(state_type),       dimension(:), target, intent(inout) :: state
      real,                                         intent(in)    :: dt
      real,                                         intent(in)    :: current_time
      
      ! Local variables
      integer                                     :: p, stat
      real,                          dimension(2) :: tmp_option_shape
      character(len=OPTION_PATH_LEN)              :: tmp_char_option
      
      ewrite(1,*) 'Initialise Darcy IMPES data'
      
      ! Auxilliary fields.
      call allocate_and_insert_auxilliary_fields(state, &
                                                 force_prescribed_diagnositc_allocate_old_iterated = .true.)
      
      di%dt = dt

      di%current_time = current_time
      
      di%nonlinear_iter = 0
      
      call get_option('/geometry/dimension', di%ndim)
      
      di%state => state
      
      ! deduce the number of phase
      di%number_phase = size(di%state)
      
      di%pressure_mesh             => extract_mesh(di%state(1), "PressureMesh")
      di%gradient_pressure_mesh    => extract_mesh(di%state(1), "GradientPressureMesh")
      di%velocity_mesh             => extract_mesh(di%state(1), "VelocityMesh")
      di%elementwise_mesh          => extract_mesh(di%state(1), "ElementWiseMesh")
      
      di%number_vele       = element_count(di%pressure_mesh)
      di%number_sele       = surface_element_count(di%pressure_mesh)
      di%number_pmesh_node = node_count(di%pressure_mesh)
      
      di%average_pressure         => extract_scalar_field(di%state(1), "AveragePressure")
      di%porosity                 => extract_scalar_field(di%state(1), "Porosity")
      di%absolute_permeability    => extract_scalar_field(di%state(1), "AbsolutePermeability")
      di%positions                => extract_vector_field(di%state(1), "Coordinate")
      di%total_darcy_velocity     => extract_vector_field(di%state(1), "TotalDarcyVelocity")
      di%total_mobility           => extract_scalar_field(di%state(1), "TotalMobility")
      di%sum_saturation           => extract_scalar_field(di%state(1), "SumSaturation")
      di%div_total_darcy_velocity => extract_scalar_field(di%state(1), "DivergenceTotalDarcyVelocity")
      
      if (have_option('/physical_parameters/gravity')) then
         
         di%have_gravity = .true.
         
         di%gravity_direction => extract_vector_field(di%state(1), "GravityDirection")
     
         call get_option('/physical_parameters/gravity/magnitude', di%gravity_magnitude)

         call allocate(di%gravity, di%positions%dim, di%gravity_direction%mesh)
            
         call set(di%gravity, di%gravity_direction)
         call scale(di%gravity, di%gravity_magnitude)
      
      else

         di%have_gravity = .false.
         
         nullify(di%gravity_direction)
         
         di%gravity_magnitude = 0.0
         
         call allocate(di%gravity, di%positions%dim, di%elementwise_mesh, field_type = FIELD_TYPE_CONSTANT)
         
         call zero(di%gravity)
      
      end if
            
      ! Form the positions on the pressure mesh
      di%positions_pressure_mesh = get_coordinate_field(di%state(1), di%pressure_mesh)

      ! Determine the pressure matrix sparsity
      di%sparsity_pmesh_pmesh => get_csr_sparsity_firstorder(di%state(1), di%pressure_mesh, di%pressure_mesh)
      
      ! Allocate the pressure matrix and lhs and rhs to use for saturations
      call allocate(di%pressure_matrix, di%sparsity_pmesh_pmesh)
      call allocate(di%lhs, di%pressure_mesh)
      call allocate(di%rhs, di%pressure_mesh)
      call allocate(di%rhs_adv, di%pressure_mesh)
      call allocate(di%rhs_time, di%pressure_mesh)
      call allocate(di%cv_mass_pressure_mesh_with_porosity, di%pressure_mesh)
      call allocate(di%cv_mass_pressure_mesh_with_old_porosity, di%pressure_mesh)
      call allocate(di%cv_mass_pressure_mesh, di%pressure_mesh)
      
      ! Allocate a field that is always zero on the pressure mesh to be pointed 
      ! at by capilliary pressure and saturation source if not required.
      ! NOTE using a FIELD_TYPE_CONSTANT causes issues.
      allocate(di%constant_zero_sfield_pmesh)
      call allocate(di%constant_zero_sfield_pmesh, di%pressure_mesh)
      
      call zero(di%constant_zero_sfield_pmesh)
      
      ! Calculate the latest CV mass on the pressure mesh with porosity
      call compute_cv_mass(di%positions, di%cv_mass_pressure_mesh_with_porosity, di%porosity)      

      ! Calculate the latest CV mass on the pressure mesh
      call compute_cv_mass(di%positions, di%cv_mass_pressure_mesh)      
      
      ! Allocate phase flags for special fields
      allocate(di%have_capilliary_pressure(di%number_phase))
      allocate(di%have_saturation_source(di%number_phase))
      
      ! Pull phase dependent fields from state
      allocate(di%pressure(di%number_phase))
      allocate(di%capilliary_pressure(di%number_phase))
      allocate(di%gradient_pressure(di%number_phase))
      allocate(di%iterated_gradient_pressure(di%number_phase))
      allocate(di%saturation(di%number_phase))
      allocate(di%old_saturation(di%number_phase))
      allocate(di%saturation_source(di%number_phase))
      allocate(di%relative_permeability(di%number_phase))
      allocate(di%viscosity(di%number_phase))      
      allocate(di%darcy_velocity(di%number_phase))
      allocate(di%cfl(di%number_phase))
      allocate(di%mobility(di%number_phase))
      allocate(di%fractional_flow(di%number_phase))
      allocate(di%density(di%number_phase))
      allocate(di%old_density(di%number_phase))

      do p = 1,di%number_phase

         di%pressure(p)%ptr                   => extract_scalar_field(di%state(p), "Pressure")
         
         if (p > 1) then
            di%capilliary_pressure(p)%ptr     => extract_scalar_field(di%state(p), "CapilliaryPressure", stat = stat)            
            if (stat == 0) then
               di%have_capilliary_pressure    = .true.
            else
               di%have_capilliary_pressure    = .false.
               di%capilliary_pressure(p)%ptr  => di%constant_zero_sfield_pmesh
            end if
         else
            ! Cannot have capilliary pressure for phase 1 as it is special
            di%have_capilliary_pressure       = .false.
            di%capilliary_pressure(p)%ptr     => di%constant_zero_sfield_pmesh
         end if
         
         di%gradient_pressure(p)%ptr          => extract_vector_field(di%state(p), "GradientPressure")
         di%iterated_gradient_pressure(p)%ptr => extract_vector_field(di%state(p), "IteratedGradientPressure")
         di%saturation(p)%ptr                 => extract_scalar_field(di%state(p), "Saturation")
         di%old_saturation(p)%ptr             => extract_scalar_field(di%state(p), "OldSaturation")
         
         di%saturation_source(p)%ptr          => extract_scalar_field(di%state(p), "SaturationSource", stat = stat)
         if (stat == 0) then
            di%have_saturation_source         = .true.
         else
            di%have_saturation_source         = .false.
            di%saturation_source(p)%ptr       => di%constant_zero_sfield_pmesh
         end if
         
         di%relative_permeability(p)%ptr      => extract_scalar_field(di%state(p), "RelativePermeability")
         di%viscosity(p)%ptr                  => extract_scalar_field(di%state(p), "Viscosity")
         di%darcy_velocity(p)%ptr             => extract_vector_field(di%state(p), "DarcyVelocity")
         di%cfl(p)%ptr                        => extract_scalar_field(di%state(p), "DarcyVelocityCFL")
         di%mobility(p)%ptr                   => extract_scalar_field(di%state(p), "Mobility")
         di%fractional_flow(p)%ptr            => extract_scalar_field(di%state(p), "FractionalFlow")
         di%density(p)%ptr                    => extract_scalar_field(di%state(p), "Density")
         di%old_density(p)%ptr                => extract_scalar_field(di%state(p), "OldDensity")
         
      end do
      
      ! Determine if the first phase pressure is prognostic, else it is prescribed
      di%first_phase_pressure_prognostic = have_option(trim(di%pressure(1)%ptr%option_path)//'/prognostic')
                  
      di%phase_one_saturation_diagnostic = have_option(trim(di%saturation(1)%ptr%option_path)//'/diagnostic')
      
      ! Determine the inverse cv mass matrix of velocity mesh
      call allocate(di%inverse_cv_mass_velocity_mesh, &
                    di%velocity_mesh)

      call compute_cv_mass(di%positions, di%inverse_cv_mass_velocity_mesh)

      call invert(di%inverse_cv_mass_velocity_mesh)   

      ! Determine the inverse cv mass matrix of pressure mesh
      call allocate(di%inverse_cv_mass_pressure_mesh, &
                    di%pressure_mesh)

      call compute_cv_mass(di%positions, di%inverse_cv_mass_pressure_mesh)

      call invert(di%inverse_cv_mass_pressure_mesh)   
      
      ! Allocate the v, pressure and saturation BC value mesh, field and flag
      
      di%bc_surface_mesh = make_mesh(di%pressure_mesh%faces%surface_mesh, &
                                     di%pressure_mesh%faces%shape, &
                                     continuity=-1, &
                                     name='DGSurfaceMesh')
               
      call allocate(di%v_bc_value, &
                   &di%bc_surface_mesh, &
                   &name='DarcyVelocityBCValue')         
      
      allocate(di%v_bc_flag(surface_element_count(di%pressure_mesh)))

      call allocate(di%pressure_bc_value, &
                   &di%bc_surface_mesh, &
                   &name='PressureBCValue')         
      
      allocate(di%pressure_bc_flag(surface_element_count(di%pressure_mesh)))
      
      ! Find the weak_pressure_bc_coefficient option
      call get_option('/weak_pressure_bc_coefficient', &
                      di%weak_pressure_bc_coeff, &
                      default = 1.0)
            
      ! Initialise a surface field used to cache the inverse characteristic 
      ! length required for weak dirichlet BCs for pressure
      call allocate(di%inverse_characteristic_length, &
                   &di%bc_surface_mesh, &
                   &name='InverseCharacteristicLength')
      
      call darcy_impes_calculate_inverse_characteristic_length(di)
      
      ! Get the relperm correlation options from the first phase field
      call get_option(trim(di%relative_permeability(1)%ptr%option_path)//&
                     &'/diagnostic/correlation/name', &
                      tmp_char_option)
      
      ! Allocate and initialise the relperm exponents for PowerLaw correlation
      allocate(di%relperm_corr_options%exponents(di%number_phase))
      di%relperm_corr_options%exponents = 0

      ! Allocate and initialise the relperm residual saturations
      allocate(di%relperm_corr_options%residual_saturations(di%number_phase))
      
      if (have_option(trim(di%relative_permeability(1)%ptr%option_path)//&
                     &'/diagnostic/correlation/residual_saturations')) then
      
         ! get the number of option residual saturations to check it is of length number_phase
         tmp_option_shape = option_shape(trim(di%relative_permeability(1)%ptr%option_path)//&
                                        &'/diagnostic/correlation/residual_saturations')

         if (tmp_option_shape(1) /= di%number_phase) then
            FLExit('To specify the residual saturations a value for each phase must be given')
         end if

         call get_option(trim(di%relative_permeability(1)%ptr%option_path)//&
                        &'/diagnostic/correlation/residual_saturations', &
                         di%relperm_corr_options%residual_saturations)         
      
      else
      
         di%relperm_corr_options%residual_saturations = 0.0
      
      end if
      
      if (trim(tmp_char_option) == 'PowerLaw') then
         
         di%relperm_corr_options%type = RELPERM_CORRELATION_POWER
                  
         ! get the number of option exponents to check it is of length number_phase
         tmp_option_shape = option_shape(trim(di%relative_permeability(1)%ptr%option_path)//&
                                        &'/diagnostic/correlation/exponents')
         
         if (tmp_option_shape(1) /= di%number_phase) then
            FLExit('To use the relative permeability correlation PowerLaw, a exponent for each phase must be specified')
         end if
         
         call get_option(trim(di%relative_permeability(1)%ptr%option_path)//&
                        &'/diagnostic/correlation/exponents', &
                         di%relperm_corr_options%exponents)         
         
      else if (trim(tmp_char_option) == 'Corey2Phase') then
         
         di%relperm_corr_options%type = RELPERM_CORRELATION_COREY2PHASE
         
         ! Check there are 2 phases
         if (di%number_phase /= 2) then
            FLExit('Cannot use the relative permeability correlation Corey2Phase if not a 2 phase simulation')
         end if                  

      else if (trim(tmp_char_option) == 'Corey2PhaseOpposite') then
         
         di%relperm_corr_options%type = RELPERM_CORRELATION_COREY2PHASEOPPOSITE
         
         ! Check there are 2 phases
         if (di%number_phase /= 2) then
            FLExit('Cannot use the relative permeability correlation Corey2PhaseOpposite if not a 2 phase simulation')
         end if
                  
      end if
      
      ! Get the relative permeability and density darcy impes cv options
      di%relperm_cv_options = darcy_impes_get_cv_options(di%relative_permeability(1)%ptr%option_path, &
                                                        &di%relative_permeability(1)%ptr%mesh%shape%numbering%family, &
                                                        &mesh_dim(di%relative_permeability(1)%ptr))
      
      ! If the relperm face value scheme is requires it get the saturation cv options
      ! Also get the minimum saturation value to use for the denominator of modrelperm
      ! and set a flag that saturation face values are require to be determined.
      ! If the relperm is to be limited allocate a field to find the modrelperm.
      if ((di%relperm_cv_options%facevalue == DARCY_IMPES_CV_FACEVALUE_RELPERMOVERSATUPWIND) .or. &
          (di%relperm_cv_options%facevalue == DARCY_IMPES_CV_FACEVALUE_RELPERMOVERSATFINITEELEMENT) .or. &
          (di%relperm_cv_options%facevalue == DARCY_IMPES_CV_FACEVALUE_RELPERMOVERSATCORRELATION) ) then
         
         ! Check there are saturation face value options
         if (.not. have_option(trim(di%saturation(1)%ptr%option_path)//'/diagnostic/face_value')) then
            FLExit('You are using RelativePermeability face_value options that require Saturation face_value options')
         end if
         
         di%saturation_cv_options = darcy_impes_get_cv_options(di%saturation(1)%ptr%option_path, &
                                                              &di%saturation(1)%ptr%mesh%shape%numbering%family, &
                                                              &mesh_dim(di%saturation(1)%ptr))

         call get_option(trim(di%relative_permeability(1)%ptr%option_path)//'/diagnostic/face_value[0]/minimum_denominator_saturation_value', &
                         di%minimum_denominator_saturation_value, &
                         default = 1.0e-06)
         
         di%determine_saturation_face_values = .true.
         
         if(di%relperm_cv_options%limit_facevalue) then
         
            call allocate(di%modified_relative_permeability, di%pressure_mesh)
         
         end if
         
      else
                  
         di%determine_saturation_face_values = .false. 
         
      end if
      
      ! Get the density darcy impes cv options from the first phase field
      di%density_cv_options = darcy_impes_get_cv_options(di%density(1)%ptr%option_path, &
                                                        &di%density(1)%ptr%mesh%shape%numbering%family, &
                                                        &mesh_dim(di%density(1)%ptr))
     
      ! Allocate crs matrices used to store the upwind scalar field values in CV assemble
      if(di%relperm_cv_options%limit_facevalue) then
         
         call allocate(di%relperm_upwind, di%sparsity_pmesh_pmesh) 

      end if
         
      if (di%determine_saturation_face_values .and. di%saturation_cv_options%limit_facevalue) then
               
         call allocate(di%saturation_upwind, di%sparsity_pmesh_pmesh) 
            
      end if
            
      if(di%density_cv_options%limit_facevalue) then
      
         call allocate(di%density_upwind, di%sparsity_pmesh_pmesh) 
      
      end if
      
      ! Determine the saturation advection subcycle options
      di%subcy_opt_sat%have = have_option(trim(complete_field_path(di%saturation(1)%ptr%option_path))//&
                             &'/number_advection_subcycle')
      
      if (di%subcy_opt_sat%have) then

         call get_option(trim(complete_field_path(di%saturation(1)%ptr%option_path))//&
                        &'/number_advection_subcycle', &
                        &di%subcy_opt_sat%number)
         
         di%subcy_opt_sat%consistent = have_option(trim(complete_field_path(di%saturation(1)%ptr%option_path))//&
                                                  &'/number_advection_subcycle/consistent_global_continuity') 

      else
         
         di%subcy_opt_sat%number     = 1
         di%subcy_opt_sat%consistent = .false.
      
      end if
      
      ! allocate tmp subcycle fields         
      allocate(di%old_saturation_subcycle(di%number_phase))         
      do p = 1, di%number_phase            
         allocate(di%old_saturation_subcycle(p)%ptr)
         call allocate(di%old_saturation_subcycle(p)%ptr, di%pressure_mesh)
      end do
            
      ! Get the maximum number of non linear iter
      call get_option('/timestepping/nonlinear_iterations', &
                      di%max_nonlinear_iter, &
                      default = 1)
      
      ! Get the maximum number of non linear iter after adapt
      call get_option('/timestepping/nonlinear_iterations/nonlinear_iterations_at_adapt', &
                      di%max_nonlinear_iter_first_timestep_after_adapt, &
                      default = 1)
            
      ! Determine the adaptive time stepping options
      di%adaptive_dt_options%have = have_option('/timestepping/adaptive_timestep')
      
      if (di%adaptive_dt_options%have) then 
         
         call get_option('/timestepping/adaptive_timestep/requested_cfl', &
                        &di%adaptive_dt_options%requested_cfl, &
                        &default = 0.5)
         
         call get_option('/timestepping/adaptive_timestep/minimum_timestep', &
                        &di%adaptive_dt_options%min_dt, &
                        &default = tiny(0.0))
         
         call get_option('/timestepping/adaptive_timestep/maximum_timestep', &
                         &di%adaptive_dt_options%max_dt, &
                         &default = huge(0.0))
         
         call get_option('/timestepping/adaptive_timestep/increase_tolerance', &
                        &di%adaptive_dt_options%increase_tolerance, &
                        &default = huge(0.0) * epsilon(0.0))
            
         di%adaptive_dt_options%min_dt_terminate_if_reached = &
        &have_option('/timestepping/adaptive_timestep/minimum_timestep/terminate_if_reached')

         di%adaptive_dt_options%at_first_dt = &
        &have_option('/timestepping/adaptive_timestep/at_first_dt')
      
      end if
      
      ! Determine the EoS options for each phase
      allocate(di%eos_options(di%number_phase))
      do p = 1, di%number_phase
         di%eos_options(p)%have_fluids_linear = &
        &have_option(trim(di%density(p)%ptr%option_path)//'/diagnostic/equation_of_state::IncompressibleLinear')
         
         if (di%eos_options(p)%have_fluids_linear) then
            call get_option(trim(di%density(p)%ptr%option_path)//'/diagnostic/equation_of_state::IncompressibleLinear/reference_density', &
                            di%eos_options(p)%fluids_linear_reference_density)
         else
            di%eos_options(p)%fluids_linear_reference_density = 0.0
         end if
      end do
            
      ! Determine the CV surface degree to use when integrating functions across them
      call get_option('/geometry/quadrature/controlvolume_surface_degree', &
                      di%cv_surface_quaddegree)
    
      ! Determine the CV faces information for the pressure mesh
      ! assuming all elements are the same type
      di%cvfaces = find_cv_faces(vertices   = ele_vertices(di%pressure_mesh, 1), &
                                 dimension  = mesh_dim(di%pressure_mesh), &
                                 polydegree = di%pressure_mesh%shape%degree, &
                                 quaddegree = di%cv_surface_quaddegree)

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
      di%p_cvshape_full = make_cv_element_shape(di%cvfaces, di%pressure_mesh%shape, &
                                                type = ELEMENT_CONTROLVOLUME_SURFACE_BODYDERIVATIVES)

      ! Generate the CV shape function that contains the derivatives with respect
      ! to the parent elements canonical coordinates evaluated at the 
      ! control volume faces for the gradient pressure mesh, assuming all elements 
      ! are the same type.      
      di%gradp_cvshape_full = make_cv_element_shape(di%cvfaces, di%gradient_pressure(1)%ptr%mesh%shape, &
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
      di%p_cvshape = make_cv_element_shape(di%cvfaces, di%pressure_mesh%shape) 

      ! Generate the CV shape function with reduced number of derivatives 
      ! evaluated at the control volume faces for the gradient pressure mesh, 
      ! assuming all elements are the same type.     
      di%gradp_cvshape = make_cv_element_shape(di%cvfaces, di%gradient_pressure(1)%ptr%mesh%shape) 

      ! Generate the CV shape function with reduced number of derivatives
      ! evaluated at the control volume faces located on the domain boundary
      ! for the positions mesh assuming all boundary elements are the same type.
      di%x_cvbdyshape = make_cvbdy_element_shape(di%cvfaces, di%positions%mesh%faces%shape)

      ! Generate the CV shape function with reduced number of derivatives
      ! evaluated at the control volume faces located on the domain boundary
      ! for the pressure mesh assuming all boundary elements are the same type.
      di%p_cvbdyshape = make_cvbdy_element_shape(di%cvfaces, di%pressure_mesh%faces%shape)

      ! Generate the CV shape function with reduced number of derivatives
      ! evaluated at the control volume faces located on the domain boundary
      ! for the gradient pressure mesh assuming all boundary elements are the same type.
      di%gradp_cvbdyshape = make_cvbdy_element_shape(di%cvfaces, di%gradient_pressure(1)%ptr%mesh%faces%shape)

      ! Get the cache detwei_normal and p_dshape options
      di%cached_face_value%cached_detwei_normal =  have_option('/cache_detwei_normal')
      di%cached_face_value%cached_p_dshape      =  have_option('/cache_p_dshape')

      ! Initialise the arrays used to cache the face values
      call darcy_impes_initialise_cached_face_value(di)
      
      ! If the first phase saturation is diagnostic then calculate it
      if (di%phase_one_saturation_diagnostic) call darcy_impes_calculate_phase_one_saturation_diagnostic(di)
      
      ! Calculate the relative permeabilities of each phase
      call darcy_impes_calculate_relperm_fields(di)
            
      ! Calculate the gradient pressure for each phase
      call darcy_impes_calculate_gradient_pressures(di)

      ! calculate generic diagnostics - DO NOT ADD 'calculate_diagnostic_variables'
      call calculate_diagnostic_variables_new(di%state)
      
      ! Calculate the non first phase pressure's, needs to be after the python fields
      call darcy_impes_calculate_non_first_phase_pressures(di)       

      ! Calculate the density field of each phase
      call darcy_impes_calculate_densities(di)
      
      ! Calculate the relperm and density first face values
      call darcy_impes_calculate_relperm_den_first_face_values(di)
      
      ! calculate the Darcy IMPES Velocity, Mobilities, Fractional flow and CFL fields, needs to be after the python fields
      call darcy_impes_calculate_vel_mob_ff_and_cfl_fields(di)

      ! Calculate the sum of the saturations
      call darcy_impes_calculate_sum_saturation(di)

      ! Copy ALL latest fields to Old and Iterated
      call copy_to_stored_values(di%state,"Old")
      call copy_to_stored_values(di%state,"Iterated")
      call relax_to_nonlinear(di%state)
      
      ewrite(1,*) 'Finished initialising Darcy IMPES data'
      
   end subroutine darcy_impes_initialise

! ----------------------------------------------------------------------------

   function darcy_impes_get_cv_options(option_path, element_family, dim) result(darcy_impes_cv_options)
      
      !!< Find the Darcy IMPES CV necessary discretisation options for relperm and density
      
      character(len=*),                  intent(in)    :: option_path
      integer,                           intent(in)    :: element_family
      integer,                           intent(in)    :: dim
      
      type(darcy_impes_cv_options_type) :: darcy_impes_cv_options
      
      ! local variables
      character(len=FIELD_NAME_LEN) :: tmpstring

      call get_option(trim(option_path)//'/diagnostic/face_value[0]/name', tmpstring)

      darcy_impes_cv_options%facevalue = darcy_impes_cv_facevalue_integer(tmpstring)

      darcy_impes_cv_options%limit_facevalue = &
      have_option(trim(option_path)//'/diagnostic/face_value[0]/limit_face_value')

      call get_option(trim(option_path)//&
                      '/diagnostic/face_value[0]/limit_face_value/limiter[0]/name', &
                      tmpstring, &
                      default = 'None')

      darcy_impes_cv_options%limiter = cv_limiter_integer(tmpstring)

      call get_option(trim(option_path)//&
                      '/diagnostic/face_value[0]/limit_face_value/limiter[0]/slopes/lower', &
                      darcy_impes_cv_options%limiter_slopes(1), &
                      default = 1.0)

      call get_option(trim(option_path)//&
                      '/diagnostic/face_value[0]/limit_face_value/limiter[0]/slopes/upper', &
                      darcy_impes_cv_options%limiter_slopes(2), &
                      default = 2.0)

      darcy_impes_cv_options%upwind_scheme = cv_upwind_scheme(option_path, element_family, dim)

  end function darcy_impes_get_cv_options

! ----------------------------------------------------------------------------

   integer function darcy_impes_cv_facevalue_integer(face_discretisation)

      character(len=*) :: face_discretisation

      select case(trim(face_discretisation))
      case ("FirstOrderUpwind")
        darcy_impes_cv_facevalue_integer = DARCY_IMPES_CV_FACEVALUE_FIRSTORDERUPWIND
      case ("FiniteElement")
        darcy_impes_cv_facevalue_integer = DARCY_IMPES_CV_FACEVALUE_FINITEELEMENT
      case ("RelPermOverSatUpwind")
        darcy_impes_cv_facevalue_integer = DARCY_IMPES_CV_FACEVALUE_RELPERMOVERSATUPWIND
      case ("RelPermOverSatFiniteElement")
        darcy_impes_cv_facevalue_integer = DARCY_IMPES_CV_FACEVALUE_RELPERMOVERSATFINITEELEMENT
      case ("RelPermOverSatCorrelation")
        darcy_impes_cv_facevalue_integer = DARCY_IMPES_CV_FACEVALUE_RELPERMOVERSATCORRELATION
      case default
        FLAbort("Unknown control volume face value scheme.")
      end select

   end function darcy_impes_cv_facevalue_integer

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_finalise(di)
      
      !!< Finalise (ie deallocate) the Darcy IMPES data
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! Deallocate, nullify or zero Darcy IMPES data
      !  - pointers to data in state are nullified
      !  - objects with memory only in darcy_impes_type are deallocated
      !  - variables in darcy_impes_type are zeroed
      
      ! The darcy_impes_type components are finalised in the same order they are initialised
      ! (almost, apart from _upwind matrices)
      
      ! local variables
      integer :: p
      
      ewrite(1,*) 'Finalise Darcy IMPES data'
      
      di%dt = 0.0

      di%current_time = 0.0
      
      di%nonlinear_iter = 0
      
      di%ndim = 0
      
      nullify(di%state)

      ! number_phase is zeroed at end as it is used for looping
      
      nullify(di%pressure_mesh)
      nullify(di%gradient_pressure_mesh)
      nullify(di%velocity_mesh)
      nullify(di%elementwise_mesh)

      di%number_vele       = 0
      di%number_sele       = 0
      di%number_pmesh_node = 0
      
      nullify(di%average_pressure)      
      nullify(di%porosity)
      nullify(di%absolute_permeability)
      nullify(di%positions)
      nullify(di%total_darcy_velocity)
      nullify(di%total_mobility)
      nullify(di%sum_saturation)
      nullify(di%div_total_darcy_velocity)
      nullify(di%gravity_direction)
      
      di%gravity_magnitude = 0.0
      
      call deallocate(di%gravity)
            
      call deallocate(di%positions_pressure_mesh)
      nullify(di%sparsity_pmesh_pmesh)
      call deallocate(di%pressure_matrix)
      call deallocate(di%lhs)
      call deallocate(di%rhs)
      call deallocate(di%rhs_adv)
      call deallocate(di%rhs_time)
      call deallocate(di%cv_mass_pressure_mesh_with_porosity)
      call deallocate(di%cv_mass_pressure_mesh_with_old_porosity)
      call deallocate(di%cv_mass_pressure_mesh)

      call deallocate(di%constant_zero_sfield_pmesh)
      deallocate(di%constant_zero_sfield_pmesh)
      
      deallocate(di%have_capilliary_pressure)
      deallocate(di%have_saturation_source)
      
      deallocate(di%pressure)
      deallocate(di%capilliary_pressure)
      deallocate(di%gradient_pressure)
      deallocate(di%iterated_gradient_pressure)
      deallocate(di%saturation)
      deallocate(di%old_saturation)
      deallocate(di%saturation_source)
      deallocate(di%relative_permeability)
      deallocate(di%viscosity)
      deallocate(di%darcy_velocity)
      deallocate(di%cfl)
      deallocate(di%mobility) 
      deallocate(di%fractional_flow) 
      deallocate(di%density) 
      deallocate(di%old_density) 
      
      di%first_phase_pressure_prognostic = .false.
                  
      di%phase_one_saturation_diagnostic = .false.
            
      call deallocate(di%inverse_cv_mass_velocity_mesh)
      call deallocate(di%inverse_cv_mass_pressure_mesh)
      
      call deallocate(di%bc_surface_mesh)      
      call deallocate(di%v_bc_value)
      deallocate(di%v_bc_flag)      
      call deallocate(di%pressure_bc_value)
      deallocate(di%pressure_bc_flag)      
      di%weak_pressure_bc_coeff = 0.0
      call deallocate(di%inverse_characteristic_length)
      
      di%relperm_corr_options%type     = 0
      deallocate(di%relperm_corr_options%exponents)
      deallocate(di%relperm_corr_options%residual_saturations)

      if(di%relperm_cv_options%limit_facevalue) then
         call deallocate(di%relperm_upwind)
         if (di%determine_saturation_face_values) then
            call deallocate(di%modified_relative_permeability)
         end if
      end if
      if (di%determine_saturation_face_values .and. di%saturation_cv_options%limit_facevalue) then
         call deallocate(di%saturation_upwind)
      end if
      if(di%density_cv_options%limit_facevalue) then
         call deallocate(di%density_upwind)
      end if
      
      di%minimum_denominator_saturation_value = 0.0
      
      di%determine_saturation_face_values = .false.
      
      di%relperm_cv_options%facevalue       = 0
      di%relperm_cv_options%limit_facevalue = .false.
      di%relperm_cv_options%limiter         = 0
      di%relperm_cv_options%limiter_slopes  = 0.0
      di%relperm_cv_options%upwind_scheme   = 0.0

      di%saturation_cv_options%facevalue       = 0
      di%saturation_cv_options%limit_facevalue = .false.
      di%saturation_cv_options%limiter         = 0
      di%saturation_cv_options%limiter_slopes  = 0.0
      di%saturation_cv_options%upwind_scheme   = 0.0
            
      di%density_cv_options%facevalue       = 0
      di%density_cv_options%limit_facevalue = .false.
      di%density_cv_options%limiter         = 0
      di%density_cv_options%limiter_slopes  = 0.0
      di%density_cv_options%upwind_scheme   = 0.0
            
      di%subcy_opt_sat%have       = .false.
      di%subcy_opt_sat%number     = 0      
      di%subcy_opt_sat%consistent = .false.      
              
      do p = 1, di%number_phase            
         call deallocate(di%old_saturation_subcycle(p)%ptr)
      end do         
      deallocate(di%old_saturation_subcycle)
            
      di%max_nonlinear_iter                            = 0
      di%max_nonlinear_iter_first_timestep_after_adapt = 0
      
      di%adaptive_dt_options%have                        = .false.
      di%adaptive_dt_options%requested_cfl               = 0.0
      di%adaptive_dt_options%min_dt                      = 0.0
      di%adaptive_dt_options%max_dt                      = 0.0
      di%adaptive_dt_options%increase_tolerance          = 0.0
      di%adaptive_dt_options%min_dt_terminate_if_reached = .false.
      di%adaptive_dt_options%at_first_dt                 = .false.
      
      deallocate(di%eos_options)
            
      di%cv_surface_quaddegree = 0
      
      call deallocate(di%cvfaces)
      call deallocate(di%x_cvshape_full)
      call deallocate(di%p_cvshape_full)
      call deallocate(di%gradp_cvshape_full)
      call deallocate(di%x_cvshape)
      call deallocate(di%p_cvshape)
      call deallocate(di%gradp_cvshape)
      call deallocate(di%x_cvbdyshape)
      call deallocate(di%p_cvbdyshape)
      call deallocate(di%gradp_cvbdyshape)

      di%cached_face_value%cached_detwei_normal = .false.
      di%cached_face_value%cached_p_dshape      = .false.
      
      deallocate(di%cached_face_value%relperm)
      deallocate(di%cached_face_value%relperm_bdy)
      deallocate(di%cached_face_value%den)
      deallocate(di%cached_face_value%den_bdy)
      if (associated(di%cached_face_value%detwei))       deallocate(di%cached_face_value%detwei)
      if (associated(di%cached_face_value%detwei_bdy))   deallocate(di%cached_face_value%detwei_bdy)
      if (associated(di%cached_face_value%normal))       deallocate(di%cached_face_value%normal)
      if (associated(di%cached_face_value%normal_bdy))   deallocate(di%cached_face_value%normal_bdy)
      if (associated(di%cached_face_value%p_dshape))     deallocate(di%cached_face_value%p_dshape)
!!!!! *** THIS IS NOT POSSIBLE YET ***
!!!!!      if (associated(di%cached_face_value%p_dshape_bdy)) deallocate(di%cached_face_value%p_dshape_bdy)
          
      ! This must be last as it is used in loops above
      di%number_phase = 0
      
      ewrite(1,*) 'Finished finalising Darcy IMPES data'
         
   end subroutine darcy_impes_finalise

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_update_post_spatial_adapt(di, &
                                                    state, &
                                                    dt, &
                                                    current_time)
      
      !!< Update the Darcy IMPES data post spatial adapt
      
      type(darcy_impes_type),                       intent(inout) :: di
      type(state_type),       dimension(:), target, intent(inout) :: state
      real,                                         intent(in)    :: dt
      real,                                         intent(in)    :: current_time
      
      ewrite(1,*) 'Update Darcy IMPES data post spatial adapt'
      
      ! ALL THE IMPORTANT SOLUTION DATA IS IN STATE
      ! WHICH IS NOT DEALLOCATED IN THE FOLLOWING PROCEDURE
      
      call darcy_impes_finalise(di)
            
      call darcy_impes_initialise(di, &
                                  state, &
                                  dt, &
                                  current_time)
      
      ewrite(1,*) 'Finished updating Darcy IMPES data post spatial adapt'
      
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
