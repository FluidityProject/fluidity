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


  ! Ordering for Intel compiler
   use quadrature
   use elements 
   use sparse_tools
   use shape_functions
   use fields
   use state_module
  

   use sparse_tools_petsc
   use AuxilaryOptions
   use MeshDiagnostics
   use signal_vars
   use spud
   use timers
   use adapt_state_module
   use FLDebug
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
   use write_gmsh
   use merge_tensors, only: merge_tensor_fields
   use form_metric_field, only: bound_metric
     
   ! *** Use Darcy IMPES module ***
   use darcy_impes_assemble_module
   
   !use Leaching chemical model***Lcai***04 July 2014****
   use darcy_impes_leaching_chemical_model
   use darcy_impes_assemble_type


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

   logical :: timing
   real :: time1, time2

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
   integer :: nonlinear_iterations_at_first_timestep, nonlinear_iterations_this_timestep
   real    :: change, chaold, nonlinear_iteration_tolerance
   integer :: number_adapts_first_timestep
   logical :: adapt_first_timestep
   
   ! *** Data associated with the Darcy IMPES solver ***  
   type(darcy_impes_type) :: di
   type(darcy_impes_type) :: di_dual
   type(state_type) , dimension(:), pointer :: state_prime, state_dual
   logical :: have_dual, have_dual_perm, solve_dual_pressure
   integer :: darcy_debug_log_unit, darcy_debug_err_unit
   integer :: ele_prt ! LCai *****Count the porosity element number 
   integer :: number_phase_prime, number_phase_dual, number_of_first_adapts, p 
   character(len=OPTION_PATH_LEN) :: phase_name
   type(tensor_field) :: metric_tensor_dual
   
   type(vector_field), pointer :: output_positions
   
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

   ! Determine the output format.
   call get_option('/io/dump_format', option_buffer)
   if(trim(option_buffer) /= "vtk") then
      ewrite(-1,*) "You must specify a dump format and it must be vtk."
      FLExit("Rejig your: /io/dump_format")
   end if

   call get_option("/timestepping/timestep", dt)   
   
   ! Set non linear iteration data 
   call get_option("/timestepping/nonlinear_iterations", &
                  & nonlinear_iterations, &
                  & default = 1)


   ! Get the number of non linear iter for first timestep
   call get_option('/timestepping/nonlinear_iterations/nonlinear_iterations_at_first_timestep', &
                   nonlinear_iterations_at_first_timestep, &
                   default = nonlinear_iterations)
   
   ! Get the number of non linear iter after adapt 
   call get_option('/timestepping/nonlinear_iterations/nonlinear_iterations_at_adapt', &
                   nonlinear_iterations_adapt, &
                   default = nonlinear_iterations)
      
   call get_option("/timestepping/nonlinear_iterations/tolerance", &
                  & nonlinear_iteration_tolerance, &
                  & default = 0.0)
   
   call get_option("/timestepping/current_time", current_time)
   call get_option("/timestepping/finish_time", finish_time)
      
   ! Auxilliary fields.
   call allocate_and_insert_auxilliary_fields(state, &
                                              force_prescribed_diagnositc_allocate_old_iterated = .true.)
   
   ! Set Time field manually as we dont call old general diagnostics routines
   time_field => extract_scalar_field(state(1), 'Time')
   call set(time_field, current_time)

   ! ***** setting up dual *****
   have_dual =  have_option("/porous_media_dual")
   have_dual_perm =  have_option("/porous_media_dual/scalar_field::AbsolutePermeabilityDual")   
   ! Decide if the dual pressure first phase should be solved for
   if (have_dual_perm) then
      solve_dual_pressure = .not. have_option("/porous_media_dual/scalar_field::AbsolutePermeabilityDual/do_not_solve_for_dual_first_phase_pressure")
   else
      solve_dual_pressure = .false.
   end if
   if (have_dual) then
      ! check there are the same number of prime and dual phases
      ! - NOTE: due to the options schemas all prime phases come first
      number_phase_prime = 0
      number_phase_dual  = 0
      do i = 1, option_count("/material_phase")
         call get_option("/material_phase["//int2str(i-1)//"]/name", phase_name)
         if (phase_name(len_trim(phase_name)-3:len_trim(phase_name)) == "dual") then
            number_phase_dual  = number_phase_dual  + 1
         else
            number_phase_prime = number_phase_prime + 1
         end if
      end do
      
      if (number_phase_prime /= number_phase_dual) then
         ewrite(-1,*) "Number of prime phases: ",number_phase_prime
         ewrite(-1,*) "Number of dual phases: ",number_phase_dual
         FLExit("The dual model must have the same number of prime and dual phases")
      end if
         
      state_prime => state(:number_phase_prime)   
      state_dual  => state(number_phase_prime+1:)
            
      call darcy_impes_initialise(di, &
                                  state_prime, &
                                  state, &
                                  dt, &
                                  current_time, &
                                  have_dual, &
                                  have_dual_perm, &
                                  solve_dual_pressure, &
                                  this_is_dual = .false.)

      call darcy_impes_initialise(di_dual, &
                                  state_dual, &
                                  state, &
                                  dt, &
                                  current_time, &
                                  have_dual, &
                                  have_dual_perm, &
                                  solve_dual_pressure, &
                                  this_is_dual = .true.)
   else
      ! *** Initialise data used in IMPES solver *** 
      call darcy_impes_initialise(di, &
                                  state, &
                                  state, &
                                  dt, &
                                  current_time, &
                                  have_dual, &
                                  have_dual_perm, &
                                  solve_dual_pressure, &
                                  this_is_dual = .false.)
   end if
   
   ! *** Darcy impes Adapt time at first time step ***
   if(di%adaptive_dt_options%have .and. di%adaptive_dt_options%at_first_dt) then
      call darcy_impes_adaptive_timestep()
   end if
   
   ! initialise the stat files
   call initialise_diagnostics(filename, state)
   call initialise_convergence(filename, state)
   call initialise_steady_state(filename, state)
   
   ! Allocate the metrics for prime and dual
   if (have_option("/mesh_adaptivity/hr_adaptivity")) then
      call allocate(metric_tensor, extract_mesh(state(1), topology_mesh_name), "ErrorMetric")
      if (have_dual) then
         call allocate(metric_tensor_dual, extract_mesh(state_dual(1), topology_mesh_name), "ErrorMetric")      
      end if
   end if
   
   ! *** Darcy impes adapt at first time ***
   first_adapt_if: if (have_option("/mesh_adaptivity/hr_adaptivity/adapt_at_first_timestep")) then
      
      call get_option("/mesh_adaptivity/hr_adaptivity/adapt_at_first_timestep/number_of_adapts", number_of_first_adapts)
      
      first_adapt_iter: do i = 1, number_of_first_adapts
                  
         call darcy_impes_adapt_and_update(initialise_fields = .true.)
         
      end do first_adapt_iter

      if(have_option("/mesh_adaptivity/hr_adaptivity/adapt_at_first_timestep/output_adapted_mesh")) then
         output_positions => extract_vector_field(state(1), "Coordinate")
         if(isparallel()) then
           call write_gmsh_file(parallel_filename("first_timestep_adapted_mesh"), output_positions)
           call write_halos("first_timestep_adapted_mesh", output_positions%mesh)
         else
           call write_gmsh_file("first_timestep_adapted_mesh", output_positions)
         end if
      end if
      
   else
   
      ! write to screen useful problem diagnostics
      call run_diagnostics(state)
   
   end if first_adapt_if
   
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
   
   ! write initial stat file diagnostics
   if(have_option("/io/stat/output_at_start")) call write_diagnostics(state, current_time, dt, timestep, not_to_move_det_yet=.true.)

   not_to_move_det_yet=.false.
         
   ! ******************************
   ! *** Start of timestep loop ***
   ! ******************************

   ! Initialise the non linear iter actual
   nonlinear_iterations_this_timestep = nonlinear_iterations_at_first_timestep
   
   timing=(debug_level()>=2)

   timestep_loop: do
      timestep = timestep + 1
      if (timing) then
         call cpu_time(time1)
      end if

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
      if (have_dual) call darcy_impes_copy_to_old(di_dual)
      
      ! *** Darcy IMPES evaluate pre solve system fields - BCs, prescribed fields and gravity ***
      call darcy_impes_evaluate_pre_solve_fields()
      
      
      ! ****** 22 Aug 2013 ******** LCai********************************
      !calculate the porosity used to solve the Mobile-immobile model
      if (size(di%MIM_options%immobile_prog_sfield) > 0) then

           if (di%prt_is_constant) then
             di%porosity_cnt = ele_val(di%porosity, 1)

           else
             call zero(di%porosity_pmesh)
             ! then do galerkin projection to project the porosity on elementwise mesh to pressure mesh
             !note that this now only support the porosity which is originally based on dg elemernwise mesh
             if (continuity(di%porosity) < 0) then !check wether the porosity is dg
               do ele_prt=1,ele_count(di%porosity_pmesh)
               call darcy_trans_assemble_galerkin_projection_elemesh_to_pmesh(di%porosity_pmesh, di%porosity, di%positions, ele_prt)
               end do
             else
               FLExit("the mesh of porosity should be elementwise dg")
             end if
           end if
       end if
       ! ************Finish *** LCai **********************************
       
      ! Solve the system of equations with a non linear iteration
      call darcy_impes_solve_system_of_equations_with_non_linear_iteration()
      
      ! exit timestep loop if non linear iteration has not converged if required
      if(have_option("/timestepping/nonlinear_iterations/terminate_if_not_converged")) then
         if(its >= nonlinear_iterations_this_timestep .and. change >= abs(nonlinear_iteration_tolerance)) then
            ewrite(0, *) "Nonlinear iteration tolerance not reached - termininating"
            exit timestep_loop
         end if
      end if
      
      ! Update the current_time variable
      current_time = current_time + DT

      ! *** Update DarcyIMPES current_time ***
      di%current_time      = current_time
      di_dual%current_time = current_time     
      call set_option("/timestepping/current_time", current_time)            
      
      ! Set Time field manually as we dont call old general diagnostics routines
      time_field => extract_scalar_field(state(1), 'Time')
      call set(time_field, current_time)
       
      ! Call the modern and significantly less satanic version of study
      call write_diagnostics(state, current_time, dt, timestep)
               
      ! hr Mesh adapt if required
      adapt_if: if(have_option("/mesh_adaptivity/hr_adaptivity")) then

         do_adapt_if: if(do_adapt_mesh(current_time, timestep)) then
            
            ! *** Darcy impes adapt mesh and update
            call darcy_impes_adapt_and_update(initialise_fields = .false.)
                     
         else
         
            ! Set the number of non linear iterations to use next time step
            nonlinear_iterations_this_timestep = nonlinear_iterations            

            ! *** Darcy IMPES adaptive time stepping choice ***
            call darcy_impes_adaptive_timestep()
                          
         end if do_adapt_if
      
      else
      
         ! Set the number of non linear iterations to use next time step
         nonlinear_iterations_this_timestep = nonlinear_iterations       

         ! *** Darcy IMPES adaptive time stepping choice ***
         call darcy_impes_adaptive_timestep()
      
      end if adapt_if
      

      if (timing) then
         call cpu_time(time2)
         ewrite(2,*) "Time spent in One loop of time step: dBp", time2-time1
      end if
      
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
   if(have_option("/mesh_adaptivity/hr_adaptivity")) then
      call deallocate(metric_tensor)
      if (have_dual) then
         call deallocate(metric_tensor_dual)
      end if
   end if
   
   !****07 July 2014 Lcai****Deallocate Leaching Chemical model******!
   if (di%lc%have_leach_chem_model) then
     call finalize_leaching_chemical_model(di)
   end if

   ! ***** Finalise dual permeability model *****
   if (have_dual) then
      call darcy_impes_finalise(di_dual, &
                                have_dual, &
                                have_dual_perm, &
                                solve_dual_pressure, &
                                this_is_dual = .true.)
   end if
   
   ! *** Finalise darcy impes variables ***
   call darcy_impes_finalise(di, &
                             have_dual, &
                             have_dual_perm, &
                             solve_dual_pressure, &
                             this_is_dual = .false.)
    
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
                                     all_state, &
                                     dt, &
                                     current_time, &
                                     have_dual, &
                                     have_dual_perm, &
                                     solve_dual_pressure, &
                                     this_is_dual)
      
      !!< Initialise the Darcy IMPES type from options and state
      
      type(darcy_impes_type),                       intent(inout) :: di
      type(state_type),       dimension(:), target, intent(inout) :: state
      type(state_type),       dimension(:), target, intent(inout) :: all_state    
      real,                                         intent(in)    :: dt
      real,                                         intent(in)    :: current_time
      logical ,                                     intent(in)    :: have_dual
      logical,                                      intent(in)    :: have_dual_perm
      logical,                                      intent(in)    :: solve_dual_pressure
      logical ,                                     intent(in)    :: this_is_dual
      
      ! Local variables
      integer                                     :: p, stat, f_count, im_count, f, f_p, im_p, imsrc_p, ele 
                                                    ! modified on 16 Aug 2013 ** LCai
      real,                          dimension(2) :: tmp_option_shape
      character(len=OPTION_PATH_LEN)              :: tmp_char_option
      
      ewrite(1,*) 'Initialise Darcy IMPES data'
      
      di%dt = dt

      di%current_time = current_time
      
      di%nonlinear_iter = 0
      
      call get_option('/geometry/dimension', di%ndim)
      
      di%state => state
      
      ! deduce the number of phase
      di%number_phase = size(di%state)
      
      di%pressure_mesh    => extract_mesh(di%state(1), "PressureMesh")
      di%elementwise_mesh => extract_mesh(di%state(1), "ElementWiseMesh")
            
      di%number_vele       = element_count(di%pressure_mesh)
      di%number_sele       = surface_element_count(di%pressure_mesh)
      di%number_pmesh_node = node_count(di%pressure_mesh)
      
      ! Allocate a field that is always zero on the element wise mesh to be pointed 
      ! at by capilliary pressure and saturation source if not required.
      ! NOTE using a FIELD_TYPE_CONSTANT causes issues.
      allocate(di%constant_zero_sfield_elementwisemesh)
      call allocate(di%constant_zero_sfield_elementwisemesh, di%elementwise_mesh)
      
      call zero(di%constant_zero_sfield_elementwisemesh)
      
      di%average_pressure         => extract_scalar_field(di%state(1), "AveragePressure")
      if (have_dual .and. this_is_dual) then 
         di%porosity              => extract_scalar_field(di%state(1), "PorosityDual")
         if (have_dual_perm) then
            di%absolute_permeability => extract_scalar_field(di%state(1), "AbsolutePermeabilityDual")
         else
            di%absolute_permeability => di%constant_zero_sfield_elementwisemesh
         end if
      else
         di%porosity              => extract_scalar_field(di%state(1), "Porosity")
         di%old_porosity          => extract_scalar_field(di%state(1), "OldPorosity")
         di%absolute_permeability => extract_scalar_field(di%state(1), "AbsolutePermeability")
      end if            
      di%positions                => extract_vector_field(di%state(1), "Coordinate")
      di%total_darcy_velocity     => extract_vector_field(di%state(1), "TotalDarcyVelocity")
      di%total_mobility           => extract_scalar_field(di%state(1), "TotalMobility")
      di%sum_saturation           => extract_scalar_field(di%state(1), "SumSaturation")
      di%div_total_darcy_velocity => extract_scalar_field(di%state(1), "DivergenceTotalDarcyVelocity")
      di%bulk_darcy_velocity      => extract_vector_field(di%state(1), "BulkDarcyVelocity")
      
      ! Allocate the gradient pressure data
      
      ! make a shape which is one degree less than the pressure mesh
      di%gradient_pressure_shape = make_element_shape(vertices = ele_loc(di%positions,1), &
                                                      dim      = di%ndim, &
                                                      degree   = di%pressure_mesh%shape%degree - 1, &
                                                      quad     = di%pressure_mesh%shape%quadrature)
      
      ! make a mesh using the new shape that is discontinuous
      di%gradient_pressure_mesh = make_mesh(di%pressure_mesh, &
                                            di%gradient_pressure_shape, &
                                            continuity = -1, &
                                            name       = 'GradientPressureMesh')
      
      allocate(di%gradient_pressure(di%number_phase))
      allocate(di%iterated_gradient_pressure(di%number_phase))
      do p = 1,di%number_phase
         allocate(di%gradient_pressure(p)%ptr)
         allocate(di%iterated_gradient_pressure(p)%ptr)
         call allocate(di%gradient_pressure(p)%ptr, di%ndim, di%gradient_pressure_mesh)
         call allocate(di%iterated_gradient_pressure(p)%ptr, di%ndim, di%gradient_pressure_mesh)          
      end do
      
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
      
      ! Allocate the matrix and lhs and rhs to use for pressure, saturations and generic coupled fields
      call allocate(di%matrix, di%sparsity_pmesh_pmesh)
      ! Only allocate the pressure matrix for the prime di call here
      if (.not. this_is_dual) then
         if (have_dual .and. have_dual_perm .and. solve_dual_pressure) then
            call allocate(di%dual_block_pressure_matrix, di%sparsity_pmesh_pmesh, blocks=(/2,2/), name ='DualBlockPressureMatrix')
         else
            call allocate(di%pressure_matrix, di%sparsity_pmesh_pmesh, name ='PressureMatrix')      
         end if
      end if
      allocate(di%pressure_rhs)
      call allocate(di%pressure_rhs, di%pressure_mesh)      
      call allocate(di%lhs, di%pressure_mesh)
      call allocate(di%rhs, di%pressure_mesh)
      call allocate(di%rhs_full, di%pressure_mesh)
      call allocate(di%rhs_high_resolution, di%pressure_mesh)
      call allocate(di%cv_mass_pressure_mesh_with_source, di%pressure_mesh)
      call allocate(di%cv_mass_pressure_mesh_with_porosity, di%pressure_mesh)
      call allocate(di%cv_mass_pressure_mesh_with_old_porosity, di%pressure_mesh)
      call allocate(di%inverse_cv_sa_pressure_mesh, di%pressure_mesh)
      call allocate(di%work_array_of_size_pressure_mesh, di%pressure_mesh)
      
      if (have_dual) then
         call allocate(di%cv_mass_pressure_mesh_with_lambda_dual, di%pressure_mesh)
         call allocate(di%rhs_dual, di%pressure_mesh)
      end if
      
      ! Allocate a field that is always zero on the pressure mesh to be pointed 
      ! at by capilliary pressure and saturation source if not required.
      ! NOTE using a FIELD_TYPE_CONSTANT causes issues.
      allocate(di%constant_zero_sfield_pmesh)
      call allocate(di%constant_zero_sfield_pmesh, di%pressure_mesh)
      
      call zero(di%constant_zero_sfield_pmesh)
      
      ! Calculate the latest CV mass on the pressure mesh with porosity
      call compute_cv_mass(di%positions, di%cv_mass_pressure_mesh_with_porosity, di%porosity)
      
      ! Allocate phase flags for special fields
      allocate(di%have_capilliary_pressure(di%number_phase))
      allocate(di%have_saturation_source(di%number_phase))      
 
      ! Pull phase dependent fields from state
      allocate(di%pressure(di%number_phase))
      allocate(di%capilliary_pressure(di%number_phase))
      allocate(di%saturation(di%number_phase))
      allocate(di%old_saturation(di%number_phase))
      allocate(di%saturation_source(di%number_phase))
      allocate(di%relative_permeability(di%number_phase))
      allocate(di%old_relative_permeability(di%number_phase))
      allocate(di%viscosity(di%number_phase))      
      allocate(di%darcy_velocity(di%number_phase))
      allocate(di%cfl(di%number_phase))
      allocate(di%mobility(di%number_phase))
      allocate(di%fractional_flow(di%number_phase))
      allocate(di%density(di%number_phase))
      allocate(di%old_density(di%number_phase))
      
      !*****************************LCai 24 July 2013**************************!
      ! Allocate the MIM model
      allocate(di%MIM_options%immobile_saturation(di%number_phase))
      allocate(di%MIM_options%old_immobile_saturation(di%number_phase))
      allocate(di%MIM_options%mobile_saturation(di%number_phase))
      allocate(di%MIM_options%old_mobile_saturation(di%number_phase))
      allocate(di%MIM_options%mass_trans_coef(di%number_phase))
      allocate(di%MIM_options%old_mass_trans_coef(di%number_phase)) 
      !flag of MIM model and mass transfer coefficient
      allocate(di%MIM_options%have_MIM(di%number_phase))
      allocate(di%MIM_options%have_mass_trans_coef(di%number_phase))
      !flag of having immobile prognostic field
      !***NOT**USED**allocate(di%MIM_options%have_immobile_prog_sfield(di%number_phase))
      !**********Finish**************LCai***************************************! 
      

      ! *****************22 Aug 2013 ***** LCai *****************************************
      !check wether the porosity is set to be constant or not
      !Mind that since the bug inside 'allocate_and_insert_scalar_field'
      !which is using 'is_constant=allocate_tensor_field_as_constant' inside the subrouting for scalar field
      !it will never return the scalar field as Field_tyoe_constant, so we need to check that explicitly
      di%prt_is_constant = .false. !this is to check the porosity 
      if (option_count(trim(di%porosity%option_path) // "/prescribed/value") == 1) then
         di%prt_is_constant = have_option(trim(di%porosity%option_path) // "/prescribed/value[0]/constant")
      end if
      ewrite(1,*) 'Is porosity a constant?', di%prt_is_constant
      !****************************Finish*** LCai ***************************************
      
      if (have_dual) then
         allocate(di%transmissibility_lambda_dual(di%number_phase))
      end if
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
            
           !*****************************LCai 24 July 2013******************!
          di%MIM_options%immobile_saturation(p)%ptr  => extract_scalar_field(di%state(p), "ImmobileSaturation", stat = stat)
          if ((stat == 0) .and. (.not. this_is_dual)) then
               di%MIM_options%have_MIM(p) = .true.
               di%MIM_options%old_immobile_saturation(p)%ptr  => extract_scalar_field(di%state(p), "OldImmobileSaturation")
               di%MIM_options%mobile_saturation(p)%ptr        => extract_scalar_field(di%state(p), "MobileSaturation")
               di%MIM_options%old_mobile_saturation(p)%ptr    => extract_scalar_field(di%state(p), "OldMobileSaturation")
            else
               di%MIM_options%have_MIM(p) = .false.
               di%MIM_options%immobile_saturation(p)%ptr  => di%constant_zero_sfield_pmesh
               di%MIM_options%mobile_saturation(p)%ptr    => di%constant_zero_sfield_pmesh
               di%MIM_options%old_immobile_saturation(p)%ptr  => di%constant_zero_sfield_pmesh
               di%MIM_options%old_mobile_saturation(p)%ptr    => di%constant_zero_sfield_pmesh
          end if
           !**********Fisnish**************LCai ****************************!
 
         else
            ! Cannot have capilliary pressure for phase 1 as it is special
            di%have_capilliary_pressure       = .false.
            di%capilliary_pressure(p)%ptr     => di%constant_zero_sfield_pmesh
            
            !*******************LCai 24 July & 08 Aug 2013******************!
             ! Cannot have MIM for phase 1 and dual phase
             di%MIM_options%have_MIM(p) = .false.
             di%MIM_options%immobile_saturation(p)%ptr  => di%constant_zero_sfield_pmesh
             di%MIM_options%mobile_saturation(p)%ptr    => di%constant_zero_sfield_pmesh
             di%MIM_options%old_immobile_saturation(p)%ptr  => di%constant_zero_sfield_pmesh
             di%MIM_options%old_mobile_saturation(p)%ptr    => di%constant_zero_sfield_pmesh
            !*********Finish***************LCai******************************!
         end if
         
         !**********************LCai 25 July 2013****************************!
         !If the MIM exists, check whether there is mass transfer coefficient
         if (di%MIM_options%have_MIM(p) ) then
            di%MIM_options%mass_trans_coef(p)%ptr => extract_scalar_field(di%state(p), "MassTransferCoefficient", stat = stat)
   
            di%MIM_options%old_mass_trans_coef(p)%ptr   => extract_scalar_field(di%state(p), "OldMassTransferCoefficient", stat = stat)
         if (stat == 0) then
                di%MIM_options%have_mass_trans_coef(p) = .true.
            else
                di%MIM_options%have_mass_trans_coef(p) = .false.
                di%MIM_options%mass_trans_coef(p)%ptr  => di%constant_zero_sfield_pmesh
                di%MIM_options%old_mass_trans_coef(p)%ptr => di%constant_zero_sfield_pmesh
            end if
         else
            !Cannot have mass transfer coefficient without MIM model
            di%MIM_options%have_mass_trans_coef(p) = .false.
            di%MIM_options%mass_trans_coef(p)%ptr  => di%constant_zero_sfield_pmesh
            di%MIM_options%old_mass_trans_coef(p)%ptr  => di%constant_zero_sfield_pmesh  
         end if
         !********Finish*********LCai 25 July 2013****************************!
           
         
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
         di%old_relative_permeability(p)%ptr  => extract_scalar_field(di%state(p), "OldRelativePermeability")
         di%viscosity(p)%ptr                  => extract_scalar_field(di%state(p), "Viscosity")
         di%darcy_velocity(p)%ptr             => extract_vector_field(di%state(p), "DarcyVelocity")
         di%cfl(p)%ptr                        => extract_scalar_field(di%state(p), "DarcyVelocityOverPorosityCFL")
         di%mobility(p)%ptr                   => extract_scalar_field(di%state(p), "Mobility")
         di%fractional_flow(p)%ptr            => extract_scalar_field(di%state(p), "FractionalFlow")
         di%density(p)%ptr                    => extract_scalar_field(di%state(p), "Density")
         di%old_density(p)%ptr                => extract_scalar_field(di%state(p), "OldDensity")
         
         if (have_dual .and. this_is_dual) then
            di%transmissibility_lambda_dual(p)%ptr => extract_scalar_field(di%state(p), "TransmissibilityLambdaDual")        
         end if 
         
      end do
      
      !**************************26 July 2013 LCai ***************************************!
      !the flag to check wether the MIM exist in at least one phase
      di%MIM_options%have_MIM_phase = .false.
      do p = 1, di%number_phase
        if (di%MIM_options%have_MIM(p)) then
           di%MIM_options%have_MIM_phase = .true.
           exit
        end if
      end do 
      ewrite(1,*) 'if have MIM ', di%MIM_options%have_MIM_phase
      !*****Finish****************26 July 2013 LCai ***************************************!
           
      
      ! Determine if the first phase pressure is prognostic, else it is prescribed
      di%first_phase_pressure_prognostic = have_option(trim(di%pressure(1)%ptr%option_path)//'/prognostic')
                  
      di%phase_one_saturation_diagnostic = have_option(trim(di%saturation(1)%ptr%option_path)//'/diagnostic')
      
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

      ! Get the data associated with generic prognostic scalar fields
      
      f_count = 0
      im_count = 0 ! *** 08 Aug 2013**LCai **initicalize count the prognositic immobile field

      do p = 1, di%number_phase
         
         im_p = 0 !  *** 16 Aug 2013 ** LCai ** count the field within each phase
         f_p = 0 ! *** 16 Aug 2013 ** LCai ** count the field within each phase
         imsrc_p = 0 ! *** 22 Aug 2013 ** LCai *** count the generic prognostic field with immobile source option within each phase

         do f = 1, option_count('/material_phase['//int2str(p-1)//']/scalar_field')
         
            if (have_option('/material_phase['//int2str(p-1)//']/scalar_field['//int2str(f-1)//']/prognostic')) then
            
               call get_option('/material_phase['//int2str(p-1)//']/scalar_field['//int2str(f-1)//']/name', &
                               tmp_char_option)
               
               if ((trim(tmp_char_option) /= 'Pressure') .and. &
                   (trim(tmp_char_option) /= 'Saturation')) then
               
                  f_count = f_count + 1
                  f_p = f_p + 1 ! *** 16 Aug 2013 ** LCai ** count the field within each phase

                  !count the number of generic prog sfield with Immobile_source_option
                  if (have_option('/material_phase['//int2str(p-1)//']/scalar_field['//int2str(f-1)//']/prognostic/ImmobileSource'))  then
                    imsrc_p = imsrc_p + 1
                  end if

               end if 
               
            end if  
           
         end do 
         
        ! *************09 Aug 2013 LCai ***************************************  
        ! count the  prognostic immobile field          
        if (di%MIM_options%have_MIM(p)) then

          do f=1, option_count('/material_phase['//int2str(p-1)//']/MobileImmobileModel/scalar_field')

            if (have_option('/material_phase['//int2str(p-1)//']/MobileImmobileModel/scalar_field['//int2str(f-1)//']/prognostic'))  then
 
              im_count = im_count + 1
              im_p = im_p + 1 ! *** 16 Aug 2013 ** LCai ** count the field within each phase
                      
            end if

          end do

        end if

        !If there exist the Immobile prog field in this phase, check wether the number of those fields equal to the generic prog sfield
        if (im_p > 0) then
          if ((im_p /= f_p).or.(im_p /= imsrc_p)) then
            print *, "This is phase", p
            FLExit('The number of the immobile prognostic sfield should either be zero or equal to the generic prog sfield within this phase')
          end if
        end if
       ! *****************finish ***LCai *************************************
      end do 
      

      ! ************ 09 & 16 Aug 2013 LCai ****************************************
      ! get the data associate with the prognostic immobile field
      allocate(di%MIM_options%immobile_prog_sfield(im_count))
      
      if (size(di%MIM_options%immobile_prog_sfield) > 0) then

        im_count = 0
       
        ! allocate the MIM sorce terms for the matrix to solve the prognostic sfield
        call allocate(di%MIM_options%MIM_src, di%pressure_mesh)


        call allocate(di%MIM_options%MIM_src_s, di%pressure_mesh)
  

        ! cannot have immobile_prog_field in phase 1
        !***NOT**USED**di%MIM_options%have_immobile_prog_sfield(1)= .false.
        
        do p = 2, di%number_phase

           if (di%MIM_options%have_MIM(p)) then

             do f=1, option_count('/material_phase['//int2str(p-1)//']/MobileImmobileModel/scalar_field')

               if (have_option('/material_phase['//int2str(p-1)//']/MobileImmobileModel/scalar_field['//int2str(f-1)//']/prognostic')) then
   
                 im_count = im_count + 1 
                 
                 call get_option('/material_phase['//int2str(p-1)//']/MobileImmobileModel/scalar_field['//int2str(f-1)//']/name', &
                                  tmp_char_option)


                 di%MIM_options%immobile_prog_sfield(im_count)%phase = p
                   
                 tmp_char_option = 'Immobile'//trim(tmp_char_option)
                 
                 di%MIM_options%immobile_prog_sfield(im_count)%sfield => extract_scalar_field(di%state(p), &
                                                                                          trim(tmp_char_option)) 

                 di%MIM_options%immobile_prog_sfield(im_count)%old_sfield => extract_scalar_field(di%state(p), &
                                                                                   'Old'//trim(tmp_char_option))
                 ! get the source_field names
                 call  get_option(trim(di%MIM_options%immobile_prog_sfield(im_count)%sfield%option_path)//'/prognostic/source_field_name', &
                                   di%MIM_options%immobile_prog_sfield(im_count)%source_name )

               end if         
             
             end do
 
           end if
                   
        end do     

      end if
      ! ***************Finish*******LCai *************************************** 

      allocate(di%generic_prog_sfield(f_count))
      
      if (size(di%generic_prog_sfield) > 0) then
      
         f_count = 0

         do p = 1, di%number_phase

            do f = 1, option_count('/material_phase['//int2str(p-1)//']/scalar_field')

               if (have_option('/material_phase['//int2str(p-1)//']/scalar_field['//int2str(f-1)//']/prognostic')) then

                  call get_option('/material_phase['//int2str(p-1)//']/scalar_field['//int2str(f-1)//']/name', &
                                  tmp_char_option)

                  if ((trim(tmp_char_option) /= 'Pressure') .and. &
                      (trim(tmp_char_option) /= 'Saturation')) then

                     f_count = f_count + 1

                     di%generic_prog_sfield(f_count)%phase = p

                     di%generic_prog_sfield(f_count)%sfield => extract_scalar_field(di%state(p), &
                                                                                    trim(tmp_char_option))

                     di%generic_prog_sfield(f_count)%old_sfield => extract_scalar_field(di%state(p), &
                                                                                        'Old'//trim(tmp_char_option))

                     di%generic_prog_sfield(f_count)%sfield_diff => extract_tensor_field(di%state(p), &
                                                                                         trim(tmp_char_option)//'Diffusivity', &
                                                                                         stat = stat)
                     if (stat == 0) then

                        di%generic_prog_sfield(f_count)%have_diff = .true.

                     else

                        nullify(di%generic_prog_sfield(f_count)%sfield_diff)

                        di%generic_prog_sfield(f_count)%have_diff = .false.

                     end if

                     di%generic_prog_sfield(f_count)%sfield_abs => extract_scalar_field(di%state(p), &
                                                                                        trim(tmp_char_option)//'Absorption', &
                                                                                        stat = stat)
                     if (stat == 0) then

                        di%generic_prog_sfield(f_count)%have_abs = .true.

                     else

                        nullify(di%generic_prog_sfield(f_count)%sfield_abs)

                        di%generic_prog_sfield(f_count)%have_abs = .false.

                     end if

                     di%generic_prog_sfield(f_count)%sfield_src => extract_scalar_field(di%state(p), &
                                                                                        trim(tmp_char_option)//'Source', &
                                                                                        stat = stat)
                     if (stat == 0) then

                        di%generic_prog_sfield(f_count)%have_src = .true.

                     else

                        nullify(di%generic_prog_sfield(f_count)%sfield_src)

                        di%generic_prog_sfield(f_count)%have_src = .false.

                     end if

                     di%generic_prog_sfield(f_count)%sfield_cv_options = &
                     darcy_impes_get_cv_options(di%generic_prog_sfield(f_count)%sfield%option_path, &
                                                di%generic_prog_sfield(f_count)%sfield%mesh%shape%numbering%family, &
                                                mesh_dim(di%generic_prog_sfield(f_count)%sfield))
                     
                     if (di%generic_prog_sfield(f_count)%sfield_cv_options%facevalue == DARCY_IMPES_CV_FACEVALUE_NONE) then
                        di%generic_prog_sfield(f_count)%have_adv = .false.
                     else
                        di%generic_prog_sfield(f_count)%have_adv = .true.
                     end if
                    
                     ! *** 09 & 13 Aug 2013 LCai***************************************
                     tmp_char_option = trim(di%generic_prog_sfield(f_count)%sfield%option_path)//'/prognostic/ImmobileSource'
                     if (have_option(tmp_char_option)) then
                       di%generic_prog_sfield(f_count)%have_MIM_source = .true.
                       ! get the source_field names
                       call get_option(trim(tmp_char_option)//"/algorithm/source_field_name", di%generic_prog_sfield(f_count)%source_name)
                       di%generic_prog_sfield(f_count)%source_name = 'Immobile'//trim(di%generic_prog_sfield(f_count)%source_name)
                     else
                       di%generic_prog_sfield(f_count)%have_MIM_source = .false.
                     end if
                     ewrite(1,*) 'Does immobile source term exist in the generic scalar field', di%generic_prog_sfield(f_count)%sfield%name, &
                                                                                               di%generic_prog_sfield(f_count)%have_MIM_source
                 ! *** Finish *** LCai **************************************
                  end if 

               end if 

            end do 

         end do

         ! ******* 22 Aug 2013 ****** LCai ***********************************
         ! check that is the immobile prog sfield exist
         ! if exist, then allocate the porosity used to calcalate MIM
         !constant or projected to pmesh
         if (size(di%MIM_options%immobile_prog_sfield) > 0) then

           if (di%prt_is_constant) then
             di%porosity_cnt = ele_val(di%porosity, 1)
             
           else
             call allocate(di%porosity_pmesh, di%pressure_mesh)
             call allocate(di%old_porosity_pmesh, di%pressure_mesh)
             call zero(di%old_porosity_pmesh)
             call zero(di%porosity_pmesh)
             ! then do galerkin projection to project the porosity on elementwise mesh to pressure mesh
             !note that this now only support the porosity which is originally based on dg elemernwise mesh
             if (continuity(di%porosity) < 0) then !check wether the porosity is dg
               do ele_prt=1,ele_count(di%porosity_pmesh)
               call darcy_trans_assemble_galerkin_projection_elemesh_to_pmesh(di%porosity_pmesh, di%porosity, di%positions, ele_prt)
               end do
               print *, di%porosity_pmesh%val
             else
               FLExit("the mesh of porosity should be elementwise dg")
             end if
           end if
         end if

         ! *********** Finish **********LCai *********************************

         call allocate(di%sfield_bc_value, &
                       &di%bc_surface_mesh, &
                       &name='GenericPrognosticScalarField'//int2str(f-1)//'BCValue')         

         allocate(di%sfield_bc_flag(surface_element_count(di%pressure_mesh)))
      
      end if
      
      ! Allocate an upwind value matrix that is used for generic scalar 
      ! fields, saturations and density if required
      call allocate(di%sfield_upwind, di%sparsity_pmesh_pmesh) 
      
      ! Get the relperm correlation options from the first phase field
      call get_option(trim(di%relative_permeability(1)%ptr%option_path)//&
                     &'/diagnostic/correlation/name', &
                      tmp_char_option)
      
      ! Allocate and initialise the relperm exponents for PowerLaw correlation
      allocate(di%relperm_corr_options%exponents(di%number_phase))
      di%relperm_corr_options%exponents = 0

      ! Allocate and initialise the relperm residual saturations
      allocate(di%relperm_corr_options%residual_saturations(di%number_phase))
      di%relperm_corr_options%residual_saturations = 0.0
      
      ! Allocate and initialise the relperm cut off saturations
      allocate(di%relperm_corr_options%cutoff_saturations(di%number_phase))
      di%relperm_corr_options%cutoff_saturations = 0.0

      ! Allocate and initialise the relperm scaling_coefficients
      allocate(di%relperm_corr_options%scaling_coefficients(di%number_phase))
      di%relperm_corr_options%scaling_coefficients = 0.0

      if (have_option(trim(di%relative_permeability(1)%ptr%option_path)//&
                     &'/diagnostic/correlation/exponents')) then

         ! get the number of option exponents to check it is of length number_phase
         tmp_option_shape = option_shape(trim(di%relative_permeability(1)%ptr%option_path)//&
                                        &'/diagnostic/correlation/exponents')

         if (tmp_option_shape(1) /= di%number_phase) then
            FLExit('To specify the relperm exponents a value for each phase must be specified')
         end if

         call get_option(trim(di%relative_permeability(1)%ptr%option_path)//&
                        &'/diagnostic/correlation/exponents', &
                         di%relperm_corr_options%exponents)

      else

         di%relperm_corr_options%exponents = 0.0

      end if
      
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

      if (have_option(trim(di%relative_permeability(1)%ptr%option_path)//&
                     &'/diagnostic/correlation/cutoff_saturations')) then

         ! get the number of option cut off saturations to check it is of length number_phase
         tmp_option_shape = option_shape(trim(di%relative_permeability(1)%ptr%option_path)//&
                                        &'/diagnostic/correlation/cutoff_saturations')

         if (tmp_option_shape(1) /= di%number_phase) then
            FLExit('To specify the cut off saturation a value for each phase must be specified')
         end if

         call get_option(trim(di%relative_permeability(1)%ptr%option_path)//&
                        &'/diagnostic/correlation/cutoff_saturations', &
                         di%relperm_corr_options%cutoff_saturations)

      else

         di%relperm_corr_options%cutoff_saturations = di%relperm_corr_options%residual_saturations

      end if

      if (have_option(trim(di%relative_permeability(1)%ptr%option_path)//&
                     &'/diagnostic/correlation/scaling_coefficients')) then

         ! get the number of option scaling_coefficients to check it is of length number_phase
         tmp_option_shape = option_shape(trim(di%relative_permeability(1)%ptr%option_path)//&
                                        &'/diagnostic/correlation/scaling_coefficients')

         if (tmp_option_shape(1) /= di%number_phase) then
            FLExit('To specify the scaling_coefficients a value for each phase must be specified')
         end if

         call get_option(trim(di%relative_permeability(1)%ptr%option_path)//&
                        &'/diagnostic/correlation/scaling_coefficients', &
                         di%relperm_corr_options%scaling_coefficients)

      else

         di%relperm_corr_options%scaling_coefficients = 1.0

      end if
      
      if (trim(tmp_char_option) == 'PowerLaw') then
         
         di%relperm_corr_options%type = RELPERM_CORRELATION_POWER
                  
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
      
      else if (trim(tmp_char_option) == 'Mineral') then
 
         di%relperm_corr_options%type = RELPERM_CORRELATION_MINERAL

      else if (trim(tmp_char_option) == 'VanGenuchten') then
 
         di%relperm_corr_options%type = RELPERM_CORRELATION_VANGENUCHTEN
         
         ! Check there are 2 phases
         if (di%number_phase /= 2) then
            FLExit('Cannot use the relative permeability correlation VanGenuchten if not a 2 phase simulation')
         end if

      else if (trim(tmp_char_option) == 'Jackson2Phase') then
         
         di%relperm_corr_options%type = RELPERM_CORRELATION_JACKSON2PHASE
         
         ! Check there are 2 phases
         if (di%number_phase /= 2) then
            FLExit('Cannot use the relative permeability correlation Jackson2Phase if not a 2 phase simulation')
         end if                  

      else if (trim(tmp_char_option) == 'Jackson2PhaseOpposite') then
         
         di%relperm_corr_options%type = RELPERM_CORRELATION_JACKSON2PHASEOPPOSITE
         
         ! Check there are 2 phases
         if (di%number_phase /= 2) then
            FLExit('Cannot use the relative permeability correlation Jackson2PhaseOpposite if not a 2 phase simulation')
         end if
                       
      end if
      
      ! Get the relative permeability and density darcy impes cv options
      di%relperm_cv_options = darcy_impes_get_cv_options(di%relative_permeability(1)%ptr%option_path, &
                                                        &di%relative_permeability(1)%ptr%mesh%shape%numbering%family, &
                                                        &mesh_dim(di%relative_permeability(1)%ptr))
      
      ! If the relperm face value scheme requires it get the saturation cv options
      ! Also get the minimum saturation value to use for the denominator of modrelperm
      ! and set a flag that saturation face values require to be determined.
      ! If the relperm is to be limited allocate a field to find the modrelperm.
      if ((di%relperm_cv_options%facevalue == DARCY_IMPES_CV_FACEVALUE_RELPERMOVERSATUPWIND) .or. &
          (di%relperm_cv_options%facevalue == DARCY_IMPES_CV_FACEVALUE_RELPERMOVERSATFINITEELEMENT) .or. &
          (di%relperm_cv_options%facevalue == DARCY_IMPES_CV_FACEVALUE_RELPERMOVERSATCORRELATION) ) then
         
         ! Check there are saturation face value options
         if (.not. have_option(trim(complete_field_path(trim(di%saturation(1)%ptr%option_path)))//'/face_value')) then
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
           
      ! Allocate crs matrice used to store the upwind relperm field values in CV assemble
      if(di%relperm_cv_options%limit_facevalue) then
         
         call allocate(di%relperm_upwind, di%sparsity_pmesh_pmesh) 

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
                        &default = 1.1)
            
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
      

      !*******04 July 2014 Lcai*Leaching chemical model*************!
      !Initialise the Leaching chemical model
      if (have_option('/Leaching_Chemical_Model')) then

        di%lc%have_leach_chem_model= .true.
        call initialize_leaching_chemical_model(di)

      else

        di%lc%have_leach_chem_model= .false.

      end if
      !********finish*******leaching chemistry model**************!

      ! If the first phase saturation is diagnostic then calculate it
      if (di%phase_one_saturation_diagnostic) call darcy_impes_calculate_phase_one_saturation_diagnostic(di)
      
      ! Calculate the relative permeabilities of each phase
      call darcy_impes_calculate_relperm_fields(di)
            
      ! Calculate the gradient pressure for each phase
      call darcy_impes_calculate_gradient_pressures(di)
      
      ! Copy gradient pressure to iterated field
      call darcy_impes_copy_to_iterated(di)

      ! calculate generic diagnostics - DO NOT ADD 'calculate_diagnostic_variables'
      call calculate_diagnostic_variables_new(all_state)
      
      !*********************26 July 2013 LCai*****************************!
      ! calculate the Mobile saturations if MIM is used
      if (di%MIM_options%have_MIM_phase) call darcy_trans_MIM_assemble_and_solve_mobile_saturation(di)
      
      !*******Finish********26 July 2013 LCai*****************************!
      
      
      ! Calculate the non first phase pressure's, needs to be after the python fields
      call darcy_impes_calculate_non_first_phase_pressures(di)       

      ! Calculate the density field of each phase
      call darcy_impes_calculate_densities(di)
      
      ! Copy field values to Old - required for saturation and relperm to calc start face values
      call copy_to_stored_values(di%state,"Old")
      
      ! Calculate the relperm and density first face values
      call darcy_impes_calculate_relperm_den_first_face_values(di)
      
      ! Calculate the latest inverse of the CV surface area on the pressure mesh
      call darcy_impes_calculate_inverse_cv_sa(di)
      
      ! calculate the Darcy IMPES Velocity, Mobilities, Fractional flow and CFL fields 
      ! - needs to be after the python fields and calculate the cv sa field
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

      call get_option(trim(complete_field_path(trim(option_path)))//'/face_value[0]/name', tmpstring)

      darcy_impes_cv_options%facevalue = darcy_impes_cv_facevalue_integer(tmpstring)

      call get_option(trim(complete_field_path(trim(option_path)))//&
                      '/face_value[0]/number_face_value_iteration', &
                      darcy_impes_cv_options%number_face_value_iteration, &
                      default = 1)

      darcy_impes_cv_options%limit_facevalue = &
      have_option(trim(complete_field_path(trim(option_path)))//'/face_value[0]/limit_face_value')

      call get_option(trim(complete_field_path(trim(option_path)))//&
                      '/face_value[0]/limit_face_value/limiter[0]/name', &
                      tmpstring, &
                      default = 'None')

      darcy_impes_cv_options%limiter = cv_limiter_integer(tmpstring)

      call get_option(trim(complete_field_path(trim(option_path)))//&
                      '/face_value[0]/limit_face_value/limiter[0]/slopes/lower', &
                      darcy_impes_cv_options%limiter_slopes(1), &
                      default = 1.0)

      call get_option(trim(complete_field_path(trim(option_path)))//&
                      '/face_value[0]/limit_face_value/limiter[0]/slopes/upper', &
                      darcy_impes_cv_options%limiter_slopes(2), &
                      default = 2.0)

      darcy_impes_cv_options%upwind_scheme = cv_upwind_scheme(option_path, element_family, dim)

  end function darcy_impes_get_cv_options

! ----------------------------------------------------------------------------

   integer function darcy_impes_cv_facevalue_integer(face_discretisation)

      character(len=*) :: face_discretisation

      select case(trim(face_discretisation))
      case ("None")
         darcy_impes_cv_facevalue_integer = DARCY_IMPES_CV_FACEVALUE_NONE
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
   
   subroutine darcy_impes_finalise(di, &
                                   have_dual, &
                                   have_dual_perm, &
                                   solve_dual_pressure, &
                                   this_is_dual)
      
      !!< Finalise (ie deallocate) the Darcy IMPES data
      
      type(darcy_impes_type), intent(inout) :: di
      logical,                intent(in)    :: have_dual
      logical,                intent(in)    :: have_dual_perm
      logical,                intent(in)    :: solve_dual_pressure      
      logical,                intent(in)    :: this_is_dual
          
      ! Deallocate, nullify or zero Darcy IMPES data
      !  - pointers to data in state are nullified
      !  - objects with memory only in darcy_impes_type are deallocated
      !  - variables in darcy_impes_type are zeroed
      
      ! The darcy_impes_type components are finalised in the same order they are initialised
      ! (almost, apart from _upwind matrices)
      
      ! local variables
      integer :: p, f
      
      ewrite(1,*) 'Finalise Darcy IMPES data'
      
      di%dt = 0.0

      di%current_time = 0.0
      
      di%nonlinear_iter = 0
      
      di%ndim = 0
      
      nullify(di%state)

      ! number_phase is zeroed at end as it is used for looping
      
      nullify(di%pressure_mesh)
      nullify(di%elementwise_mesh)
            
      di%number_vele       = 0
      di%number_sele       = 0
      di%number_pmesh_node = 0
      
      call deallocate(di%constant_zero_sfield_elementwisemesh)
      deallocate(di%constant_zero_sfield_elementwisemesh)
      
      nullify(di%average_pressure)      
      nullify(di%porosity)
      nullify(di%old_porosity)
      nullify(di%absolute_permeability)
      nullify(di%positions)
      nullify(di%total_darcy_velocity)
      nullify(di%total_mobility)
      nullify(di%sum_saturation)
      nullify(di%div_total_darcy_velocity)
      nullify(di%bulk_darcy_velocity)
      nullify(di%gravity_direction)

      call deallocate(di%gradient_pressure_shape)
      call deallocate(di%gradient_pressure_mesh)
      
      do p = 1,di%number_phase
         call deallocate(di%gradient_pressure(p)%ptr)
         call deallocate(di%iterated_gradient_pressure(p)%ptr)          
         deallocate(di%gradient_pressure(p)%ptr)
         deallocate(di%iterated_gradient_pressure(p)%ptr)
         ! *** 09 Aug 2013 LCai *************************
         nullify(di%capilliary_pressure(p)%ptr)
         nullify(di%saturation(p)%ptr)
         nullify(di%old_saturation(p)%ptr)
         nullify(di%saturation_source(p)%ptr)
         nullify(di%relative_permeability(p)%ptr)
         nullify(di%old_relative_permeability(p)%ptr)
         nullify(di%viscosity(p)%ptr)
         nullify(di%darcy_velocity(p)%ptr)
         nullify(di%cfl(p)%ptr)
         nullify(di%mobility(p)%ptr)
         nullify(di%fractional_flow(p)%ptr)
         nullify(di%density(p)%ptr)
         nullify(di%old_density(p)%ptr)
         nullify(di%MIM_options%old_immobile_saturation(p)%ptr)
         nullify(di%MIM_options%immobile_saturation(p)%ptr)
         nullify(di%MIM_options%mobile_saturation(p)%ptr)
         nullify(di%MIM_options%old_mobile_saturation(p)%ptr)
         nullify(di%MIM_options%mass_trans_coef(p)%ptr)
         nullify(di%MIM_options%old_mass_trans_coef(p)%ptr)
         ! *** Finish ***LCai ***************************
      end do
      deallocate(di%gradient_pressure)
      deallocate(di%iterated_gradient_pressure)
      
      di%gravity_magnitude = 0.0
      
      call deallocate(di%gravity)
            
      call deallocate(di%positions_pressure_mesh)
      nullify(di%sparsity_pmesh_pmesh)
      call deallocate(di%matrix)
      if (.not. this_is_dual) then 
         if (have_dual .and. have_dual_perm .and. solve_dual_pressure) then
            call deallocate(di%dual_block_pressure_matrix)
         else
            call deallocate(di%pressure_matrix)     
         end if
      end if
      call deallocate(di%pressure_rhs)
      deallocate(di%pressure_rhs)
      call deallocate(di%lhs)
      call deallocate(di%rhs)
      call deallocate(di%rhs_full)
      call deallocate(di%rhs_high_resolution)
      call deallocate(di%cv_mass_pressure_mesh_with_source)     
      call deallocate(di%cv_mass_pressure_mesh_with_porosity)
      call deallocate(di%cv_mass_pressure_mesh_with_old_porosity)
      call deallocate(di%inverse_cv_sa_pressure_mesh)
      call deallocate(di%work_array_of_size_pressure_mesh)
      
      if (have_dual) then
         call deallocate(di%cv_mass_pressure_mesh_with_lambda_dual)
         call deallocate(di%rhs_dual)
      end if

      call deallocate(di%constant_zero_sfield_pmesh)
      deallocate(di%constant_zero_sfield_pmesh)
      
      deallocate(di%have_capilliary_pressure)
      deallocate(di%have_saturation_source)
     
       ! **** 09 Aug 2013 *** LCai *************************
       if (size(di%MIM_options%immobile_prog_sfield) > 0) then
 
         do f = 1, size(di%MIM_options%immobile_prog_sfield)
           nullify(di%MIM_options%immobile_prog_sfield(f)%sfield)
           nullify(di%MIM_options%immobile_prog_sfield(f)%old_sfield)
         end do
         call deallocate(di%MIM_options%MIM_src)
         call deallocate(di%MIM_options%MIM_src_s)
         
         if (.not.(di%prt_is_constant)) then
           call deallocate(di%porosity_pmesh)
           call deallocate(di%old_porosity_pmesh)
         end if

       end if
       deallocate(di%MIM_options%immobile_prog_sfield) 
      ! ***** Finish **** LCai **************************** 

 
      !*****************************LCai 25 July & 08 & 22 Aug 2013******************!
      di%prt_is_constant= .false.
      ! Deallocate the MIM model
      di%MIM_options%have_mass_trans_coef = .false.
      !***NOT**USED**di%MIM_options%have_immobile_prog_sfield = .false.
      di%MIM_options%have_MIM = .false.
      di%MIM_options%have_MIM_phase = .false.
      deallocate(di%MIM_options%immobile_saturation)
      deallocate(di%MIM_options%old_immobile_saturation)
      deallocate(di%MIM_options%old_mobile_saturation)
      deallocate(di%MIM_options%mobile_saturation)
      deallocate(di%MIM_options%mass_trans_coef)
      deallocate(di%MIM_options%have_MIM)
      deallocate(di%MIM_options%have_mass_trans_coef)
      deallocate(di%MIM_options%old_mass_trans_coef)
      !***NOT**USED**deallocate(di%MIM_options%have_immobile_prog_sfield) 
      !**********Finish**************LCai **************************************!

      deallocate(di%pressure)
      deallocate(di%capilliary_pressure)
      deallocate(di%saturation)
      deallocate(di%old_saturation)
      deallocate(di%saturation_source)
      deallocate(di%relative_permeability)
      deallocate(di%old_relative_permeability)
      deallocate(di%viscosity)
      deallocate(di%darcy_velocity)
      deallocate(di%cfl)
      deallocate(di%mobility) 
      deallocate(di%fractional_flow) 
      deallocate(di%density) 
      deallocate(di%old_density) 
      
      if (have_dual .and. this_is_dual) then
         deallocate(di%transmissibility_lambda_dual)
      end if
      
      di%first_phase_pressure_prognostic = .false.
                  
      di%phase_one_saturation_diagnostic = .false.
            
      call deallocate(di%inverse_cv_mass_pressure_mesh)
      
      call deallocate(di%bc_surface_mesh)      
      call deallocate(di%v_bc_value)
      deallocate(di%v_bc_flag)      
      call deallocate(di%pressure_bc_value)
      deallocate(di%pressure_bc_flag)      
      di%weak_pressure_bc_coeff = 0.0
      call deallocate(di%inverse_characteristic_length)

      if (size(di%generic_prog_sfield) > 0) then
         call deallocate(di%sfield_bc_value)
         deallocate(di%sfield_bc_flag)      
      
         do f = 1, size(di%generic_prog_sfield)
            nullify(di%generic_prog_sfield(f)%sfield)
            nullify(di%generic_prog_sfield(f)%old_sfield)
            nullify(di%generic_prog_sfield(f)%sfield_diff)
            nullify(di%generic_prog_sfield(f)%sfield_abs)
            nullify(di%generic_prog_sfield(f)%sfield_src)
            
            di%generic_prog_sfield(f)%have_diff = .false.
            di%generic_prog_sfield(f)%have_abs  = .false.
            di%generic_prog_sfield(f)%have_src  = .false.
            di%generic_prog_sfield(f)%have_adv  = .false.

            ! *** 09 Aug 2013 LCai *************************
            di%generic_prog_sfield(f)%have_MIM_source  = .false.
            ! *** Finish ******LCai
            
            di%generic_prog_sfield(f)%sfield_cv_options%facevalue                   = 0
            di%generic_prog_sfield(f)%sfield_cv_options%number_face_value_iteration = 0
            di%generic_prog_sfield(f)%sfield_cv_options%limit_facevalue             = .false.
            di%generic_prog_sfield(f)%sfield_cv_options%limiter                     = 0
            di%generic_prog_sfield(f)%sfield_cv_options%limiter_slopes              = 0.0
            di%generic_prog_sfield(f)%sfield_cv_options%upwind_scheme               = 0
         end do

      end if
      deallocate(di%generic_prog_sfield)
      call deallocate(di%sfield_upwind)
      di%relperm_corr_options%type     = 0

      
      deallocate(di%relperm_corr_options%exponents)
      deallocate(di%relperm_corr_options%residual_saturations)
      deallocate(di%relperm_corr_options%cutoff_saturations)
      deallocate(di%relperm_corr_options%scaling_coefficients)

      if(di%relperm_cv_options%limit_facevalue) then
         call deallocate(di%relperm_upwind)
         if (di%determine_saturation_face_values) then
            call deallocate(di%modified_relative_permeability)
         end if
      end if
      
      di%minimum_denominator_saturation_value = 0.0
      
      di%determine_saturation_face_values = .false.
      
      di%relperm_cv_options%facevalue                   = 0
      di%relperm_cv_options%number_face_value_iteration = 0
      di%relperm_cv_options%limit_facevalue             = .false.
      di%relperm_cv_options%limiter                     = 0
      di%relperm_cv_options%limiter_slopes              = 0.0
      di%relperm_cv_options%upwind_scheme               = 0

      di%saturation_cv_options%facevalue                   = 0
      di%saturation_cv_options%number_face_value_iteration = 0
      di%saturation_cv_options%limit_facevalue             = .false.
      di%saturation_cv_options%limiter                     = 0
      di%saturation_cv_options%limiter_slopes              = 0.0
      di%saturation_cv_options%upwind_scheme               = 0
            
      di%density_cv_options%facevalue                   = 0
      di%density_cv_options%number_face_value_iteration = 0
      di%density_cv_options%limit_facevalue             = .false.
      di%density_cv_options%limiter                     = 0
      di%density_cv_options%limiter_slopes              = 0.0
      di%density_cv_options%upwind_scheme               = 0
            
      di%subcy_opt_sat%have       = .false.
      di%subcy_opt_sat%number     = 0      
      di%subcy_opt_sat%consistent = .false.      
              
      do p = 1, di%number_phase            
         call deallocate(di%old_saturation_subcycle(p)%ptr)
      end do         
      deallocate(di%old_saturation_subcycle)
            
      di%nonlinear_iter                   = 0
      di%max_nonlinear_iter_this_timestep = 0
      
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
                                                    all_state, &
                                                    dt, &
                                                    current_time, &
                                                    have_dual, &
                                                    have_dual_perm, &
                                                    solve_dual_pressure, &
                                                    this_is_dual)
      
      !!< Update the Darcy IMPES data post spatial adapt
      
      type(darcy_impes_type),                       intent(inout) :: di
      type(state_type),       dimension(:), target, intent(inout) :: state
      type(state_type),       dimension(:), target, intent(inout) :: all_state
      real,                                         intent(in)    :: dt
      real,                                         intent(in)    :: current_time
      logical,                                      intent(in)    :: have_dual  
      logical,                                      intent(in)    :: have_dual_perm
      logical,                                      intent(in)    :: solve_dual_pressure
      logical,                                      intent(in)    :: this_is_dual  
        
      ewrite(1,*) 'Update Darcy IMPES data post spatial adapt'
      
      ! ALL THE IMPORTANT SOLUTION DATA IS IN STATE
      ! WHICH IS NOT DEALLOCATED IN THE FOLLOWING PROCEDURE
      
      call darcy_impes_finalise(di, &
                                have_dual, &
                                have_dual_perm, &
                                solve_dual_pressure, &
                                this_is_dual)
            
      call darcy_impes_initialise(di, &
                                  state, &
                                  all_state, &
                                  dt, &
                                  current_time, &
                                  have_dual, &
                                  have_dual_perm, &
                                  solve_dual_pressure, &
                                  this_is_dual)
      
      ewrite(1,*) 'Finished updating Darcy IMPES data post spatial adapt'
      
   end subroutine darcy_impes_update_post_spatial_adapt

! --------------------------------------------------------------------------------

   subroutine darcy_impes_adapt_and_update(initialise_fields)
      
      !!< Adapt the dary impes model and update as required
      
      logical, intent(in) :: initialise_fields

      ! Form metric to adapt mesh to for prime and dual      
      if (have_dual) then
         call qmesh(state_prime, metric_tensor)
         call qmesh(state_dual, metric_tensor_dual)
         call merge_tensor_fields(metric_tensor, metric_tensor_dual)
         call bound_metric(metric_tensor, state(1))
      else
         call qmesh(state, metric_tensor)
      end if
      
      ! write to screen useful problem diagnostics
      call run_diagnostics(state)

      ! adapt state including mesh to mesh interpolation
      call adapt_state(state, metric_tensor, initialise_fields = initialise_fields)

      ! Set the number of non linear iterations to use next time step
      nonlinear_iterations_this_timestep = nonlinear_iterations_adapt

      ! re-allocate the adaptivity metric
      if(have_option("/mesh_adaptivity/hr_adaptivity")) then
         call allocate(metric_tensor, extract_mesh(state(1), topology_mesh_name), "ErrorMetric")
      end if

      ! re allocate and insert Old and Iterated fields (plus a few specials)
      call allocate_and_insert_auxilliary_fields(state, &
                                                 force_prescribed_diagnositc_allocate_old_iterated = .true.)

      if (have_dual) then
         ! ***** Update DUAL Darcy IMPES post spatial adapt *****
         state_prime => state(:number_phase_prime)   
         state_dual  => state(number_phase_prime+1:)
         call darcy_impes_update_post_spatial_adapt(di, &
                                                    state_prime, &
                                                    state, &
                                                    dt, &
                                                    current_time, &
                                                    have_dual, &
                                                    have_dual_perm, &
                                                    solve_dual_pressure, &
                                                    this_is_dual = .false.)

         call darcy_impes_update_post_spatial_adapt(di_dual, &
                                                    state_dual, &
                                                    state, &
                                                    dt, &
                                                    current_time, &
                                                    have_dual, &
                                                    have_dual_perm, &
                                                    solve_dual_pressure, &
                                                    this_is_dual = .true.)
      else                                                                 
         ! *** Update Darcy IMPES post spatial adapt ***
         call darcy_impes_update_post_spatial_adapt(di, &
                                                    state, &
                                                    state, &
                                                    dt, &
                                                    current_time, &
                                                    have_dual, &
                                                    have_dual_perm, &
                                                    solve_dual_pressure, &
                                                    this_is_dual = .false.)
      end if  

      ! *** Darcy IMPES adaptive time stepping choice ***
      call darcy_impes_adaptive_timestep()

      ! write to screen useful problem diagnostics          
      call run_diagnostics(state)      
      
   end subroutine darcy_impes_adapt_and_update

! --------------------------------------------------------------------------------
   
   subroutine darcy_impes_adaptive_timestep()
      
      !!< Adapt the timestep and make sure all variables are synced

      if(di%adaptive_dt_options%have) then
         call darcy_impes_calculate_cflnumber_field_based_dt(di)
         if (have_dual) then
            call darcy_impes_calculate_cflnumber_field_based_dt(di_dual)
            dt = min(di%dt, di_dual%dt)  
         else 
            dt = di%dt  
         end if
         call set_option("/timestepping/timestep", dt)
      end if
      
      ! Constrain the timestep such that the current time will not overshoot the
      ! finish time.  But spread this operation over two timesteps so as not to
      ! end up with a very small time step.  Neglecting to do this may cause
      ! subtle divide-by-huge-number bugs.
      if ( (current_time + 2*dt > finish_time) .and. (current_time + dt < finish_time) ) then
         ! penultimate time step
         ewrite(1,*) 'Constrain timestep size such as to not overshoot the finish time nor end up too small'
         dt = 0.5*(finish_time - current_time) + epsilon(0.0)
         ewrite(1,*) 'Constrained timestep size: ',dt
      else if ( current_time + dt > finish_time ) then
         ! last time step
         ewrite(1,*) 'Constrain timestep size such as to not overshoot the finish time'
         dt = finish_time - current_time + epsilon(0.0)
         ewrite(1,*) 'Constrained timestep size: ',dt
      end if      
      
      ! *** setting di time step ***
      if (have_dual) then     
         di_dual%dt = dt
      end if
      di%dt = dt
      
      call set_option("/timestepping/timestep", dt)
      
   end subroutine darcy_impes_adaptive_timestep

! --------------------------------------------------------------------------------
   
   subroutine darcy_impes_evaluate_pre_solve_fields()
      
      !!< Evaluate pre system solve fields - BCs, prescribed fields and gravity 

      ! evaluate BC fields at shifted time level
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
      if (di_dual%have_gravity .and. have_dual) then
         call set(di_dual%gravity, di_dual%gravity_direction)
         call scale(di_dual%gravity, di_dual%gravity_magnitude)
         ewrite_minmax(di_dual%gravity)
      end if  
            
   end subroutine darcy_impes_evaluate_pre_solve_fields

! --------------------------------------------------------------------------------

   subroutine darcy_impes_solve_system_of_equations_with_non_linear_iteration()
      !!< Solve the system of equations with a non linear iteration

      ! *** Darcy IMPES set max number non linear iterations for this time step ***
      di%max_nonlinear_iter_this_timestep = nonlinear_iterations_this_timestep
      di_dual%max_nonlinear_iter_this_timestep = nonlinear_iterations_this_timestep     
      
      nonlinear_iteration_loop: do its = 1,nonlinear_iterations_this_timestep

         ewrite(1,*)'###################'
         ewrite(1,*)'Start of another nonlinear iteration; its,nonlinear_iterations_this_timestep=',its,nonlinear_iterations_this_timestep
         ewrite(1,*)'###################'         
         
         call copy_to_stored_values(state, "Iterated")
         
         ! *** Darcy IMPES copy non state data to iterated ***
         call darcy_impes_copy_to_iterated(di)                  
         if (have_dual) call darcy_impes_copy_to_iterated(di_dual)         
         
         ! *** Set the Darcy IMPES nonlinear_iter
         di%nonlinear_iter      = its
         di_dual%nonlinear_iter = its
         
         ! ***** Point other porous media pressures and transmissibility_lambda as required, only used in assemble *****
         if (have_dual) then
            di%pressure_other_porous_media      => di_dual%pressure
            di_dual%pressure_other_porous_media => di%pressure
            
            do p = 1,di%number_phase
               di%transmissibility_lambda_dual(p)%ptr => di_dual%transmissibility_lambda_dual(p)%ptr
            end do
         end if
         
         ! *** Darcy IMPES Calculate the relperm and density first face values (depend on upwind direction) ***
         call darcy_impes_calculate_relperm_den_first_face_values(di)
         if (have_dual) call darcy_impes_calculate_relperm_den_first_face_values(di_dual)
         
         !*********11 Aug 2014 *****Lcai ******leaching chemical model************!
         !calculate the chemical reactions explicitly based on concentrations of the previous time step
         !before calculating the prognostic scalar fields
         call calculate_leaching_chemical_model(di)
         !*********finish****Lcai********Leaching chemical model

         ! *** Solve the Darcy equations using IMPES in three parts ***
         call darcy_impes_assemble_and_solve_part_one(di, have_dual)
         if (have_dual) call darcy_impes_assemble_and_solve_part_one(di_dual, have_dual)
         
         ! This one solves for the pressure phase 1
         call darcy_impes_assemble_and_solve_part_two(di, &
                                                      di_dual, &
                                                      have_dual, &
                                                      solve_dual_pressure)

         call darcy_impes_assemble_and_solve_part_three(di, have_dual)
         if (have_dual) call darcy_impes_assemble_and_solve_part_three(di_dual, have_dual)
          
         ! calculate generic diagnostics - DO NOT ADD 'calculate_diagnostic_variables'
         call calculate_diagnostic_variables_new(state, exclude_nonrecalculated = .true.)
         
         ! *** Calculate the Darcy IMPES Velocity, Mobilities, Fractional flow and CFL fields
         !     needs to be after the python fields ***
         call darcy_impes_calculate_vel_mob_ff_and_cfl_fields(di)
         if (have_dual) call darcy_impes_calculate_vel_mob_ff_and_cfl_fields(di_dual)

         if(nonlinear_iterations_this_timestep > 1) then
            
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
      
   end subroutine darcy_impes_solve_system_of_equations_with_non_linear_iteration

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
!    Copyright (C) 2006 Imperial College London and others.
!
!    Please see the AUTHORS file in the main source directory for a full list
