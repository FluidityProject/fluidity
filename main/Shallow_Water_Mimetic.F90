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
  program shallow_water
    use spud
    use signals
    use bubble_tools
    use fields
    use state_module
    use FLDebug
    use populate_state_module
    use write_state_module
    use vtk_interfaces
    use timeloop_utilities
    use sparsity_patterns_meshes
    use sparse_matrices_fields
    use solvers
    use diagnostic_variables
    use diagnostic_fields_wrapper
    use hybridized_helmholtz
    use advection_local_dg
    use global_parameters, only: option_path_len, python_func_len, current_time, dt
    use memory_diagnostics
    use reserve_state_module
    use boundary_conditions_from_options
      use diagnostic_fields_new, only : &
    & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
    & check_diagnostic_dependencies
    use iso_c_binding
    use mangle_options_tree
    use manifold_tools
    use FEFields
    use field_copies_diagnostics
    use slope_limiters_dg
    implicit none
#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif

    ! Interface blocks for the initialisation routines we need to call
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

    type(state_type), dimension(:), pointer :: states
    character(len = OPTION_PATH_LEN) :: simulation_name

    integer :: timestep
    integer :: ierr
    integer, save :: dump_no=0
    real :: energy

#ifdef HAVE_MPI
    call mpi_init(ierr)
    assert(ierr == MPI_SUCCESS)
#endif

#ifdef HAVE_PETSC
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

    call python_init
    call read_command_line
    call mangle_options_tree_forward
    ! Establish signal handlers
    call initialise_signals()

    timestep=0

    call populate_state(states)

!    call test_mesh(states(1))

    ! No support for multiphase or multimaterial at this stage.
    if (size(states)/=1) then
       FLExit("Multiple material_phases are not supported")
    end if

    call insert_time_in_state(states)

    ! Check the diagnostic field dependencies for circular dependencies
    call check_diagnostic_dependencies(states)

    call get_option('/simulation_name',simulation_name)
    call initialise_diagnostics(trim(simulation_name),states)
    
    call setup_fields(states(1))
    call setup_pv(states(1))

    call calculate_diagnostic_variables(states)
    call calculate_diagnostic_variables_new(states)

    call get_option("/timestepping/timestep", dt)
    call write_diagnostics(states, current_time, dt, timestep)

    ! Always output the initial conditions.
    
    call output_state(states)
    call compute_energy_hybridized(states(1),energy)

    timestep_loop: do
       timestep=timestep+1
       if (simulation_completed(current_time, timestep)) exit timestep_loop
       ewrite (1,*) "SW: start of timestep ", timestep, current_time
       call execute_timestep(states(1))

       call project_local_to_cartesian(states(1))
       call calculate_diagnostic_variables(states,&
            & exclude_nonrecalculated = .true.)
       call calculate_diagnostic_variables_new(states,&
            & exclude_nonrecalculated = .true.)

       call set_prescribed_field_values(states, &
            exclude_interpolated=.true., &
            exclude_nonreprescribed=.true., time=current_time+dt)

       call advance_current_time(current_time, dt)
       if (simulation_completed(current_time, timestep)) exit timestep_loop

       if (do_write_state(current_time, timestep)) then
          call output_state(states)
       end if

       !Update the variables
       call compute_energy_hybridized(states(1),energy)
       call write_diagnostics(states,current_time, dt, timestep)

       ewrite(1,*) 'END OF TIMESTEP   TIME=', current_time
    end do timestep_loop

    call compute_energy_hybridized(states(1),energy)
    ewrite(2,*) 'Energy = ',energy

    ! One last dump
    call project_local_to_cartesian(states(1))
    call output_state(states)
    call write_diagnostics(states,current_time, dt, timestep)

    call deallocate(states)
    call deallocate_transform_cache
    call deallocate_reserve_state
    call close_diagnostic_files
    call uninitialise_diagnostics
    ! Clean up registered diagnostics
    call destroy_registered_diagnostics 

    call print_references(0)

#ifdef HAVE_MEMORY_STATS
    call print_current_memory_stats(0)
#endif

#ifdef HAVE_MPI
    call mpi_finalize(ierr)
    assert(ierr == MPI_SUCCESS)
#endif
  contains

    subroutine test_mesh(state)

      type(state_type), intent(inout) :: state

      type(scalar_field), pointer :: test_field
      integer:: ele,ele2,face,face2,ni
      integer, dimension(:), pointer:: neigh
      
      test_field => extract_scalar_field(state,"P1DGtest")

      do ele = 1, element_count(test_field)
         neigh => ele_neigh(test_field,ele)

         do ni=1, size(neigh)
            ele2 = neigh(ni)
            face = ele_face(test_field,ele,ele2)
            face2 = ele_face(test_field,ele2,ele)
            call test_mesh_face(test_field,ele,ele2,face,face2)
         end do

      end do
    end subroutine test_mesh

    subroutine test_mesh_face(test_field,ele,ele2,face,face2)

      type(scalar_field), intent(in) :: test_field
      integer, intent(in) :: ele, ele2, face, face2

      real, dimension(face_ngi(test_field,face)) :: v1,v2
    
      v1=face_val_at_quad(test_field,face)
      v2=face_val_at_quad(test_field,face2)
  
      if (maxval(abs(v1-v2))>1.0e-10) then
         ewrite(0,*) "Hello dham, why am I back-to-front?"
         ewrite(0,*) "Mesh info: ele, ele2, face, face2", ele, ele2, face, face2
         ewrite(0,*) "face_val_at_quad(test_field,face) = ", v1
         ewrite(0,*) "face_val_at_quad(test_field,face2) = ", v2

      end if
 
    end subroutine test_mesh_face

    subroutine execute_timestep(state)
      type(state_type), intent(inout) :: state
      !
      type(vector_field), pointer :: U, advecting_u, U_old
      type(scalar_field), pointer :: D_old, D, &
           & Coriolis, PV, &
           & PV_old
      type(vector_field) :: newU, MassFlux, PVFlux, UResidual
      type(scalar_field) :: newD, DResidual
      type(scalar_field), pointer :: PVtracer, pvtracer_old

      integer :: nonlinear_iterations, nits, stat
      logical :: have_pv_tracer, &
           & have_PV=.false., have_PV_old=.false.
      real :: theta
      real :: total_divergence

      call get_option("/timestepping/nonlinear_iterations"&
           &,nonlinear_iterations)
      call get_option("/timestepping/theta",theta)

      !Set up iterative solutions and provide initial guess
      !From previous timestep
      
      U => extract_vector_field(state, "LocalVelocity")
      U_old => extract_vector_field(state, "OldLocalVelocity")
      D => extract_scalar_field(state, "LayerThickness")
      D_old => extract_scalar_field(state, "OldLayerThickness")
      Coriolis => extract_scalar_field(state, "Coriolis")
      PV => extract_scalar_field(state, "PotentialVorticity",stat)
      have_pv = (stat==0)
      if(have_pv) then
         PV_old => extract_scalar_field(state, "OldPotentialVorticity",stat)
         have_pv_old = (stat==0)
      end if
      PVtracer => extract_scalar_field(state, "PotentialVorticityTracer",stat)
      if(stat==0) then
         have_pv_tracer = .true.
         pvtracer_old => extract_scalar_field(state, &
              "OldPotentialVorticityTracer")
      else
         have_pv_tracer = .false.
      end if
      advecting_u=>extract_vector_field(state, "NonlinearVelocity")

      total_divergence = 0.
      call get_total_divergence(total_divergence,state)
      ewrite(2,*) 'TOTAL DIVERGENCE = ', total_divergence

      !Set up old values
      call set(U_old,U)
      call set(D_old,D)
      if(.not.have_option('/material_phase::Fluid/vector_field::Velocity/pro&
           &gnostic/wave_equation/')) then
         FLExit('Need wave equation option')
      end if
      if(have_option('/material_phase::Fluid/vector_field::Velocity/prognostic/wave_equation/no_wave_equation_step')) then
         !   !Just advance the D and PV fields.
         call set(advecting_u, u)
         call allocate(MassFlux,mesh_dim(U),u%mesh,'MassFlux')
         call zero(massFlux)
         ! no flux for now - just testing advection (jemma 7/3/13)
!         call solve_advection_dg_subcycle("LayerThickness", state, &
!              "NonlinearVelocity",continuity=.true.)
         call solve_advection_dg_subcycle("LayerThickness", state, &
              "NonlinearVelocity",continuity=.true.,Flux=MassFlux)
         if(have_pv_tracer) then
            call set(pvtracer_old,pvtracer)
            call allocate(PVFlux,mesh_dim(U),u%mesh,'PVFlux')
            call zero(PVFlux)
            call solve_advection_cg_tracer(PVtracer,D,&
                 d_old,MassFlux,advecting_u,PVFlux,state)
            call deallocate(PVFlux)
         end if
         call deallocate(MassFlux)
      else
         call allocate(newU,U%dim,U%mesh,"NewLocalVelocity")
         call allocate(newD,D%mesh,"NewLayerThickness")
         call set(newD,D)
         call set(newU,U)
         
         if(have_option('/material_phase::Fluid/vector_field::Velocity/prognostic/wave_equation/just_wave_equation_step')) then
            !Just solve the linear equations
            ewrite(2,*) 'CJC newU newD', maxval(abs(newU%val)), maxval(abs(newD%val))
            call solve_hybridised_timestep_residual(state,newU,newD)
         else
            !Solve the fully coupled SWE

            !Allocate some temporary fields
            call allocate(MassFlux,mesh_dim(U),u%mesh,'MassFlux')
            call allocate(PVFlux,mesh_dim(U),u%mesh,'PVFlux')
            call zero(MassFlux)
            call zero(PVFlux)
            call allocate(UResidual, mesh_dim(U), u%mesh,&
                 & 'VelocityResidual')
            call allocate(DResidual, D%mesh, 'LayerThicknessResidual')
            call zero(UResidual)
            call zero(DResidual)

            !PV calculation
            if(.not.have_pv) then
               FLExit('Need PV field')
            end if
            call get_PV(state,PV,U_old,D_old)
            call set(PV_old,PV)

            !Start of Newton iteration loop
            newton_iteration: do nits = 1, nonlinear_iterations
               ewrite(1,*) 'NONLINEAR ITERATION', nits
               !Set up advecting velocity 
               call set(advecting_u, u)
               call scale(advecting_u, 1-theta)
               call addto(advecting_u, newU, theta)

               !Compute D residual
               call solve_advection_dg_subcycle("LayerThickness", state, &
                    "NonlinearVelocity",continuity=.true.,Flux=MassFlux)
               !Should be newD-D
               call set(DResidual,newD)
               call addto(DResidual,D,-1.0)
               call apply_dg_mass(DResidual,state)

               !Update PV
               call set(PV,PV_old)
               call solve_advection_cg_tracer(PV,D,D_old,&
                    MassFlux,advecting_u,PVFlux,state)

               ewrite(1,*) 'PVFlux', maxval(abs(PVFlux%val))

               ! Compute U residual
               call compute_U_residual(UResidual,U_old,D_old,&
                    newU,newD,PVFlux,state)

               !Perform (quasi)Newton iteration
               ewrite(2,*) 'CJC newU newD', maxval(abs(newU%val)), maxval(abs(newD%val))
               call solve_hybridised_timestep_residual(state,newU,newD,&
                    & UResidual,DResidual)
            end do newton_iteration

            !Final PV calculation for output
            call get_PV(state,PV,newU,newD)         

            call deallocate(MassFlux)
            call deallocate(PVFlux)
            call deallocate(UResidual)
            call deallocate(DResidual)
            
         end if

         ewrite(1,*) 'jump in D', maxval(abs(d%val-newd%val))
         ewrite(1,*) 'jump in U', maxval(abs(U%val-newU%val))
         
         call set(D,newD)
         call set(U,newU)

         call project_to_constrained_space(state,U)

         call deallocate(newU)
         call deallocate(newD)

      end if

    end subroutine execute_timestep

    subroutine setup_fields(state)
      type(state_type), intent(inout) :: state
      type(vector_field), pointer :: v_field,U,X,Z
      type(scalar_field), pointer :: s_field,D,f_ptr
      type(mesh_type), pointer :: v_mesh
      type(vector_field) :: U_local, advecting_u, old_U, linear_orography_term
      character(len=PYTHON_FUNC_LEN) :: coriolis
      type(scalar_field) :: f, new_s_field
      integer :: stat

      X=>extract_vector_field(state, "Coordinate")
      U=>extract_vector_field(state, "Velocity")
      Z=>extract_vector_field(state, "Vorticity",stat)
      D=>extract_scalar_field(state, "LayerThickness")

      if(have_option('/geometry/mesh::CoordinateMesh/recompute_coordinate_f&
           &ield/python')) then
         call recompute_coordinate_field(state)
         call vtk_write_fields('RecomputedCoordinates', position=X, &
              model=D%mesh, &
              vfields=(/X/))
      end if

      !SET UP LOCAL VELOCITY
      !This needs an option to switch on as we don't always want to do it.
      !   !project velocity into div-conforming space
      call allocate(U_local, mesh_dim(U), U%mesh, "LocalVelocity")
      call zero(U_local)
      call insert(state, U_local, "LocalVelocity")
      call deallocate(U_local)
      v_field => extract_vector_field(state, "Velocity")
      call project_cartesian_to_local(state, v_field)

      !Advecting velocity (just used for courant number)
      call allocate(advecting_u, mesh_dim(U), U%mesh, "NonlinearVelocity")
      call zero(advecting_u)
      call insert(state, advecting_u, "NonlinearVelocity")
      call deallocate(advecting_u)

      !Old velocity (used for residual calculations)
      V_field => extract_vector_field(state, "LocalVelocity")
      call allocate(old_U, v_field%dim, v_field%mesh, "OldLocalVelocity")
      call zero(old_U)
      call insert(state, old_U, "OldLocalVelocity")
      call deallocate(old_U)

      !Old layer thickness (used for advection dg)
      D => extract_scalar_field(state, "LayerThickness")
      call allocate(new_s_field, D%mesh, "OldLayerThickness")
      call zero(new_s_field)
      call insert(state, new_s_field, "OldLayerThickness")
      call deallocate(new_s_field)

      s_field => extract_scalar_field(&
           state, "PotentialVorticityTracer",stat)
      if(stat==0) then
         call allocate(new_s_field, s_field%mesh,&
              "OldPotentialVorticityTracer")
         call zero(new_s_field)
         call insert(state, new_s_field, "OldPotentialVorticityTracer")
         call deallocate(new_s_field)
      end if

      s_field => extract_scalar_field(&
           state, "PotentialVorticity",stat)
      if(stat==0) then
         call allocate(new_s_field, s_field%mesh,&
              "OldPotentialVorticity")
         call zero(new_s_field)
         call insert(state, new_s_field, "OldPotentialVorticity")
         call deallocate(new_s_field)
      end if
      
      !SET UP CORIOLIS FORCE
      f_ptr => extract_scalar_field(state,"Coriolis",stat=stat)
      if(stat.ne.0) then
         v_mesh => extract_mesh(state,"CoordinateMesh")
         call allocate(f, v_mesh, "Coriolis")
         call get_option("/physical_parameters/coriolis", coriolis, stat)
         if(stat==0) then
            call set_from_python_function(f, coriolis, X, time=0.0)
         else
            call zero(f)
         end if
         call insert(state, f, "Coriolis")
         if(f%mesh%shape%numbering%type==ELEMENT_BUBBLE) then
            call fix_bubble_component(f)
         end if
         call deallocate(f)
      end if

      !sort out streamfunction in bubble space
      s_field => extract_scalar_field(state,"Streamfunction",stat)
      if(stat==0) then
         if(have_option('/material_phase::Fluid/scalar_field::Streamfunction&
              &/prescribed/initialise_from_quadrature_points')) then
            call initialise_from_quadrature_points(state,s_field)
         else
            if(s_field%mesh%shape%numbering%type==ELEMENT_BUBBLE) then
               call fix_bubble_component(s_field)
            end if
         end if
      end if

      !VARIOUS BALANCED INITIAL OPTIONS
      ! Geostrophic balanced initial condition, if required
      if(have_option("/material_phase::Fluid/vector_field::Velocity/prognostic&
           &/initial_condition::WholeMesh/balanced")) then
         call set_velocity_from_geostrophic_balance_hybridized(state)
      end if
      
      !Set velocity from streamfunction
      if(have_option("/material_phase::Fluid/vector_field::Velocity/prognost&
           &ic/initial_condition::WholeMesh/from_streamfunction")) then
         call set_velocity_from_streamfunction(state)
      end if

      !Set velocity from Galerkin projection
      if(have_option("/material_phase::Fluid/vector_field::Velocity/prognost&
           &ic/initial_condition::WholeMesh/galerkin_projection")) then
         call set_velocity_galerkin_projection(state)
      end if

      v_field => extract_vector_field(state, "LocalVelocity")
      call project_to_constrained_space(state,v_field)

      if(have_option("/material_phase::Fluid/scalar_field::LayerThickness/pr&
           &ognostic/initial_condition::ProjectionFromPython")) then
         call set_layerthickness_projection(state)
      end if
      if(have_option("/material_phase::Fluid/scalar_field::PrescribedLayerDe&
           &pthFromProjection")) then
         call set_layerthickness_projection(state,&
              &"PrescribedLayerDepthFromProjection")
      end if
      if(have_option("/material_phase::Fluid/scalar_field::Orography/prescribe&
           &d/projected")) then
         call set_layerthickness_projection(state,"Orography")
      end if
      if(have_option("/material_phase::Fluid/scalar_field::Orography/prescribe&
           &d/limit_values")) then
         s_field => extract_scalar_field(state, "Orography")
         call limit_vb_manifold(state,s_field)
      end if
      call fix_layerdepth_mean(state)
      if(have_option("/material_phase::Fluid/scalar_field::Orography/prescribe&
           &d/subtract_from_layer_thickness")) then
         s_field => extract_scalar_field(state, "Orography")
         D => extract_scalar_field(state, "LayerThickness")
         D%val = D%val - s_field%val
      end if
      s_field => extract_scalar_field(state, "Orography",stat)
      print*, "JEMMA:", stat, have_option("/material_phase::Fluid/scalar_field::Orography")
      print*, have_option('/material_phase::Fluid/vector_field::Velocity/prognostic/wave_equation/just_wave_equation_step')
      if(stat==0 .and. have_option('/material_phase::Fluid/vector_field::Velocity/prognostic/wave_equation/just_wave_equation_step')) then
         print*, 'calculating linear orography term'
         ! calculate LinearOrographyTerm field = g<div w,h_orog>
         call allocate(linear_orography_term, mesh_dim(U), U%mesh, "LinearOrographyTerm")
         call calculate_linear_orography(state, s_field, linear_orography_term)
         ewrite_minmax(linear_orography_term)
         call insert(state, linear_orography_term, "LinearOrographyTerm")
         call deallocate(linear_orography_term)
      end if

      s_field => extract_scalar_field(state,"PotentialVorticityTracer",stat)
      if(stat==0) then
         if(s_field%mesh%shape%numbering%type==ELEMENT_BUBBLE) then
            call fix_bubble_component(s_field)
         end if
      end if
      
      !Initial layer thickness
      s_field=>extract_scalar_field(state, "InitialLayerThickness"&
           &,stat)    
      if(stat==0) then
         s_field%val = D%val
      end if

    end subroutine setup_fields

    subroutine recompute_coordinate_field(state)
      type(state_type), intent(inout) :: state
      !
      type(vector_field), pointer :: X
      integer :: ele
      character(len=PYTHON_FUNC_LEN) :: Python_Function
      
      X => extract_vector_field(state,"Coordinate")
      call get_option('/geometry/mesh::CoordinateMesh/recompute_coordinate_f&
           &ield/python',Python_Function)

      do ele = 1, ele_count(X)
         call recompute_coordinate_field_ele(X,Python_Function,ele)
      end do
    end subroutine recompute_coordinate_field

    subroutine recompute_coordinate_field_ele(X,Python_Function,ele)
      type(vector_field), intent(inout) :: X
      character(len=PYTHON_FUNC_LEN), intent(in) :: Python_Function
      integer, intent(in) :: ele
      !
      real, dimension(X%dim,ele_loc(X,ele)) :: X_ele_val,X_ele_val_2
      integer :: stat

      X_ele_val = ele_val(X,ele)
      if(X%dim.ne.3) then
         FLExit('Option only available for 3 dimensional coordinates.')
      end if

      call set_vector_field_from_python(python_function, len(python_function),&
           & dim=3,nodes=ele_loc(X,ele),x=X_ele_val(1,:),y=X_ele_val(2,:)&
           &,z=x_ele_val(3,:),t=0.0,result_dim=3,&
           & result_x=X_ele_val_2(1,:),&
           & result_y=X_ele_val_2(2,:),&
           & result_z=X_ele_val_2(3,:),&
           & stat=stat)
    if(stat /= 0) then
       FLAbort('Failed to set new coordinate values from Python.')
    end if

    call set(X,ele_nodes(X,ele),X_ele_val_2)

    end subroutine recompute_coordinate_field_ele

    subroutine set_velocity_galerkin_projection(state)
      type(state_type), intent(in) :: state
      !
      type(vector_field), pointer  :: X, U
      integer :: ele
      character(len=PYTHON_FUNC_LEN) :: Python_Function

      ewrite(1,*) 'set_velocity_galerkin_projection'

      X=>extract_vector_field(state,"Coordinate")
      U=>extract_vector_field(state,"LocalVelocity")
      call get_option('/material_phase::Fluid/vector_field::Velocity/prognos&
           &tic/initial_condition::WholeMesh/galerkin_projection/python'&
           &,Python_Function)

      do ele = 1, ele_count(X)
         call set_velocity_galerkin_projection_ele(U,Python_Function,X,ele)
      end do
    end subroutine set_velocity_galerkin_projection

    subroutine set_velocity_galerkin_projection_ele(U,Python_Function,X,ele)
      !! Subroutine to 
      type(vector_field), intent(in) :: X
      character(len=PYTHON_FUNC_LEN), intent(in) :: Python_Function
      type(vector_field), intent(inout) :: U
      integer, intent(in) :: ele
      !
      real, dimension(ele_ngi(X,ele)) :: detwei, detJ
      real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
      real, dimension(mesh_dim(U),ele_loc(U,ele)) :: U_val
      real, dimension(X%dim,ele_ngi(X,ele)) :: X_quad, U_quad
      real, dimension(mesh_dim(U),ele_ngi(X,ele)) :: Ul_quad
      real, dimension(mesh_dim(U)*ele_loc(U,ele)) :: U_rhs
      real, dimension(mesh_dim(U)*ele_loc(U,ele),&
           mesh_dim(U)*ele_loc(U,ele)) :: U_mat
      real, dimension(mesh_dim(U), mesh_dim(U), ele_ngi(U,ele)) :: Metric
      integer :: dim1, dim2, stat, gi
      real, dimension(mesh_dim(U),mesh_dim(U),&
           & ele_loc(U,ele),ele_loc(U,ele)) :: l_u_mat
      type(element_type), pointer :: u_shape
      integer, dimension(mesh_dim(U)) :: U_start, U_end

      if(X%dim.ne.3) then
         FLExit('Option only available for 3 dimensional coordinates.')
      end if
      u_shape => ele_shape(U,ele)
      X_quad = ele_val_at_quad(X,ele)

      call set_vector_field_from_python(python_function, len(python_function),&
           & dim=3,nodes=ele_ngi(X,ele),x=X_quad(1,:),y=X_quad(2,:)&
           &,z=x_quad(3,:),t=0.0,result_dim=3,&
           & result_x=U_quad(1,:),&
           & result_y=U_quad(2,:),&
           & result_z=U_quad(3,:),&
           & stat=stat)

      call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), detwei=detwei&
           &,J=J, detJ=detJ)

      do gi=1,ele_ngi(U,ele)
         Metric(:,:,gi)=matmul(J(:,:,gi), transpose(J(:,:,gi)))/detJ(gi)
         ul_quad(:,gi) = matmul(J(:,:,gi),U_quad(:,gi))
      end do

      !velocity mass matrix (done in local coordinates)
      l_u_mat = shape_shape_tensor(u_shape, u_shape, &
           u_shape%quadrature%weight, Metric)
      !test function times ul_quad (stored temporarily in u_val)
      u_val = shape_vector_rhs(u_shape,ul_quad,u_shape%quadrature%weight)

      !indices
      do dim1 = 1, mesh_dim(U)
         u_start(dim1) = ele_loc(U,ele)*(dim1-1)+1
         u_end(dim1) = ele_loc(U,ele)*dim1
      end do
      
      do dim1 = 1, mesh_dim(U)
         u_rhs(u_start(dim1):u_end(dim1)) = u_val(dim1,:)
         do dim2 = 1, mesh_dim(U)
            u_mat(u_start(dim1):u_end(dim1),&
                 u_start(dim2):u_end(dim2))=&
                 & l_u_mat(dim1,dim2,:,:)
         end do
      end do

      call solve(u_mat,u_rhs)

      do dim1 = 1, mesh_dim(U)
         u_val(dim1,:) = u_rhs(u_start(dim1):u_end(dim1))
      end do

      call set(u,ele_nodes(u,ele),u_val)
      
    end subroutine set_velocity_galerkin_projection_ele

    subroutine fix_layerdepth_mean(state)
      type(state_type), intent(inout) :: state
      real :: D0, h_mean, area
      integer :: ele
      type(vector_field), pointer  :: X
      type(scalar_field), pointer :: D

      ewrite(1,*) 'fix_layerdepth_mean(state)'

      X=>extract_vector_field(state,"Coordinate")
      D=>extract_scalar_field(state,"LayerThickness")
      
      h_mean = 0.
      area = 0.
      do ele = 1, element_count(D)
         call assemble_mean_ele(D,X,h_mean,area,ele)
      end do
      h_mean = h_mean/area

      if(have_option('/material_phase::Fluid/scalar_field::LayerThickness/pr&
           &ognostic/mean_layer_thickness/reset_this_value')) then
         ewrite(2,*) 'Setting D0 to mean value computed from field',h_mean
         call set_option('/material_phase::Fluid/scalar_field::LayerThicknes&
              &s/prognostic/mean_layer_thickness',h_mean)
         D0 = h_mean
      else
         D%val = D%val - h_mean

         h_mean = 0.
         area = 0.
         do ele = 1, element_count(D)
            call assemble_mean_ele(D,X,h_mean,area,ele)
         end do
         h_mean = h_mean/area
         ewrite(1,*) 'area cjc', area, h_mean
         assert(abs(h_mean)<1.0e-8)

         !Add back on the correct mean depth
         call get_option("/material_phase::Fluid/scalar_field::LayerThickness/&
              &prognostic/mean_layer_thickness",D0)
         ewrite(2,*) 'Setting field mean value to D0',D0
         
         D%val = D%val + D0
      end if
      h_mean = 0.
      area = 0.
      do ele = 1, element_count(D)
         call assemble_mean_ele(D,X,h_mean,area,ele)
      end do
      ewrite(2,*) 'fix_layerdepth_mean H mean',h_mean/area,D0
      assert(abs(h_mean/area-D0)/(max(1.0,D0))<1.0e-8)
      ewrite(1,*) 'END fix_layerdepth_mean(state)'
  end subroutine fix_layerdepth_mean

  subroutine setup_pv(state)
      type(state_type), intent(inout) :: state
      !
      type(vector_field), pointer :: U
      type(scalar_field), pointer :: D,Q
      integer :: stat
      !
      Q => extract_scalar_field(state,"PotentialVorticity",stat)
      if(stat==0) then
         U => extract_vector_field(state,"LocalVelocity")
         D => extract_scalar_field(state,"LayerThickness")
         call get_PV(state,Q,U,D)         
      end if
    end subroutine setup_pv
    
    subroutine set_velocity_from_streamfunction(state)
      type(state_type), intent(inout) :: state
      !
      type(vector_field), pointer :: U, down, X
      type(scalar_field), pointer :: psi
      integer :: ele

      X => extract_vector_field(state,"Coordinate")
      psi => extract_scalar_field(state,"Streamfunction")
      u => extract_vector_field(state,"LocalVelocity")
      down=>extract_vector_field(state, "GravityDirection")

      do ele = 1, ele_count(u)
         call set_velocity_from_streamfunction_ele(U,psi,down,X,ele)
      end do

    end subroutine set_velocity_from_streamfunction

    subroutine set_velocity_from_streamfunction_ele(U,psi,down,X,ele)
      type(vector_field), intent(inout) :: U
      type(vector_field), intent(in) :: down,X
      type(scalar_field), intent(in) :: psi
      integer, intent(in) :: ele
      !
      real, dimension(ele_loc(U,ele),ele_loc(U,ele)) :: &
           & l_mass_mat
      type(element_type) :: u_shape, psi_shape
      real, dimension(down%dim, ele_ngi(down,ele)) :: up_gi
      integer :: orientation, dim1
      real, dimension(ele_loc(psi,ele)) :: psi_loc
      real, dimension(mesh_dim(psi),ele_ngi(psi,ele)) :: dpsi_gi, grad_psi_gi
      real, dimension(mesh_dim(U),ele_loc(U,ele)) :: U_loc

      !Set U = grad perp psi
      !Can do in local coordinates
      !Need to evaluate at U DOFS, easiest to do by
      !doing L2 projection which is equivalent to identity

      u_shape = ele_shape(U,ele)
      l_mass_mat = shape_shape(u_shape,u_shape,U_shape%quadrature%weight)

      psi_shape = ele_shape(psi,ele)
      up_gi = -ele_val_at_quad(down,ele)
      call get_up_gi(X,ele,up_gi,orientation)

      !Streamfunction at node values
      psi_loc = ele_val(psi,ele)
      !Skew gradient of streamfunction at quadrature points
      select case(mesh_dim(psi))
      case (2)
         grad_psi_gi = ele_grad_at_quad(psi,ele,psi_shape%dn)
         dpsi_gi(1,:) = -grad_psi_gi(2,:)
         dpsi_gi(2,:) = grad_psi_gi(1,:)
      case default
         FLAbort('Exterior derivative not implemented for given mesh dimension')
      end select
      dpsi_gi = orientation*dpsi_gi
      U_loc = shape_vector_rhs(u_shape,dpsi_gi,U_shape%quadrature%weight)
      
      do dim1 = 1, U%dim
         call solve(l_mass_mat,U_loc(dim1,:))
      end do

      call set(U,ele_nodes(U,ele),U_loc)

    end subroutine set_velocity_from_streamfunction_ele

    subroutine get_total_divergence(total_divergence,state)
      real, intent(inout) :: total_divergence
      type(state_type), intent(in) :: state
      !
      type(vector_field), pointer :: U
      type(scalar_field), pointer :: D
      integer :: ele

      U => extract_vector_field(state, "LocalVelocity")
      D => extract_scalar_field(state, "LayerThickness")

      do ele = 1, ele_count(U)
         call get_total_divergence_ele(total_divergence,U,D,ele)
      end do
    end subroutine get_total_divergence
    subroutine get_total_divergence_ele(total_divergence,U,D,ele)
      type(vector_field), intent(in) :: U
      type(scalar_field), intent(in) :: D
      integer, intent(in) :: ele
      real, intent(inout) :: total_divergence
      !
      real, dimension(mesh_dim(U),ele_loc(U,ele)) ::&
           & U_vals
      real, dimension(mesh_dim(U),ele_loc(U,ele),ele_loc(D,ele)) ::&
           & l_div_mat
      
      type(element_type), pointer :: u_shape,d_shape
      integer :: dim1, i, j
      
      U_vals = ele_val(U,ele)
      u_shape => ele_shape(U,ele)
      d_shape => ele_shape(D,ele)

      l_div_mat = dshape_shape(u_shape%dn,d_shape,D_shape%quadrature%weight)      
      !real, dimension(mesh_dim(U),ele_loc(U,ele),ele_loc(D,ele)) ::&
      !     & l_div_mat
      
      do dim1 = 1, mesh_dim(U)
         do i = 1, ele_loc(U,ele)
            do j = 1, ele_loc(D,ele)
               total_divergence = total_divergence + &
                    & l_div_mat(dim1,i,j)*U_vals(dim1,i)
            end do
         end do
      end do
    end subroutine get_total_divergence_ele
       
       subroutine apply_dg_mass(s_field,state)
      type(state_type), intent(inout) :: state
      type(scalar_field), intent(inout) :: s_field
      !
      integer :: ele      
      type(vector_field), pointer :: X

      X=>extract_vector_field(state, "Coordinate")
      do ele = 1, ele_count(s_field)
         call apply_dg_mass_ele(s_field,X,ele)
      end do
    end subroutine apply_dg_mass

    subroutine apply_dg_mass_ele(s_field,X,ele)
      type(scalar_field), intent(inout) :: s_field
      type(vector_field), intent(in) :: X
      integer, intent(in) :: ele
      !
      real, dimension(ele_ngi(s_field,ele)) :: s_gi, detwei
      real, dimension(ele_loc(s_field,ele)) :: s_rhs
      real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J

      call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), detwei=detwei&
           &,J=J)
      s_gi = ele_val_at_quad(s_field,ele)
      s_rhs = shape_rhs(ele_shape(s_field,ele),s_gi*detwei)
      call set(s_field,ele_nodes(s_field,ele),s_rhs)
    end subroutine apply_dg_mass_ele
      
    subroutine fix_bubble_component(s_field)
      type(scalar_field), intent(inout) :: s_field
      !
      integer :: ele
      real, dimension(7) :: N_vals
      !!Debugging
      real, dimension(7,7) :: lmat
      real, dimension(2,2,7,7) :: ldmat
      type(element_type), pointer :: sshape

      sshape => ele_shape(s_field,1)

      lmat = shape_shape(sshape,sshape,sshape%quadrature%weight)
      ldmat = dshape_outer_dshape(sshape%dn,sshape%dn,sshape%quadrature%weight)
      !Basis functions evaluated at bubble node.
      N_vals = eval_shape(ele_shape(s_field,1), (/1.0/3.0,1.0/3.0,1.0/3.0/))

      assert(s_field%mesh%shape%numbering%type==ELEMENT_BUBBLE)
      ewrite(2,*) 'CJC FIXING', s_field%name
      do ele = 1, ele_count(s_field)
         call fix_bubble_component_ele(s_field,N_vals,ele)
      end do
    end subroutine fix_bubble_component

    subroutine fix_bubble_component_ele(s_field,N_vals,ele)
      type(scalar_field), intent(inout) :: s_field
      real, dimension(:), intent(in) :: N_vals
      integer, intent(in) :: ele
      !
      integer, pointer, dimension(:) :: f_nodes
      real, dimension(ele_loc(s_field,ele)) :: f_vals

      f_nodes => ele_nodes(s_field,ele)
      f_vals = ele_val(s_field,ele)
      f_vals(7) = f_vals(7)-sum(f_vals(1:6)*N_vals(1:6))
      f_vals(7) = f_vals(7)/N_vals(7)
      
      call set(s_field,f_nodes(7),f_vals(7))

    end subroutine fix_bubble_component_ele
    
    subroutine initialise_from_quadrature_points(state,T)
      type(scalar_field), intent(inout) :: T
      type(state_type), intent(inout) :: state
      !
      integer :: ele
      type(csr_sparsity), pointer :: sparsity
      type(csr_matrix) :: init_mat
      type(scalar_field) :: T_rhs
      type(vector_field), pointer :: X      
      character(len=PYTHON_FUNC_LEN) :: Python_Function

      sparsity => get_csr_sparsity_firstorder(state, T%mesh, T%mesh)
      call allocate(init_mat,sparsity)
      call zero(init_mat)
      call allocate(T_rhs,T%mesh,'T_rhs')
      call zero(T_rhs)
      X=>extract_vector_field(state, "Coordinate")

       call get_option(&
            trim(T%option_path)//"/prescribed/value/python",Python_Function)
      
      do ele = 1, ele_count(T)
         call initialise_from_quadrature_points_ele(&
              T_rhs,init_mat,X,Python_Function,ele)
      end do
      ewrite(1,*) maxval(abs(T_rhs%val)), 'T_rhs'
      call petsc_solve(T,init_mat,T_rhs,option_path=trim(T%option_path)//"/prescribed/")
      call deallocate(T_rhs)
      call deallocate(init_mat)
    end subroutine initialise_from_quadrature_points
    
    subroutine initialise_from_quadrature_points_ele(T_rhs,init_mat,X,&
         Python_Function,ele)
      type(scalar_field), intent(inout) :: T_rhs
      type(csr_matrix), intent(inout) :: init_mat
      integer, intent(in) :: ele
      character(len=PYTHON_FUNC_LEN), intent(in) :: Python_Function
      type(vector_field), intent(in) :: X
      !
      real, dimension(X%dim, ele_ngi(X, ele)) :: X_quad
      real, dimension(ele_ngi(X, ele)) :: T_quad, detwei
      real, dimension(ele_loc(T_rhs, ele)) :: l_t_rhs
      real, dimension(ele_loc(T_rhs, ele), ele_loc(T_rhs,ele)) :: l_mat
      real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
      integer :: stat
      X_quad = ele_val_at_quad(X,ele)

      call set_scalar_field_from_python(python_function, len(python_function),&
           & dim=3,nodes=ele_ngi(X,ele),x=X_quad(1,:),y=X_quad(2,:)&
           &,z=x_quad(3,:),t=0.0,&
           & result=T_quad,&
           & stat=stat)
      if(stat /= 0) then
         FLAbort('Failed to set face values from Python.')
      end if

      call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J=J, &
           detwei=detwei)
      
      l_mat = shape_shape(Ele_shape(T_rhs,Ele),Ele_shape(T_rhs,Ele),detwei)
      l_t_rhs = shape_rhs(Ele_shape(T_rhs,Ele),detwei*T_quad)

      call addto(init_mat,ele_nodes(T_rhs,ele),ele_nodes(T_rhs,ele),l_mat)
      call addto(T_rhs,ele_nodes(T_rhs,ele),l_t_rhs)
      
    end subroutine initialise_from_quadrature_points_ele

    subroutine advance_current_time(current_time, dt)
      implicit none
      real, intent(inout) :: current_time
      real, intent(in) :: dt

      ! Adaptive timestepping could go here.

      current_time=current_time + dt
      call set_option("/timestepping/current_time", current_time)

    end subroutine advance_current_time

    subroutine insert_time_in_state(state)
      type(state_type), dimension(:), intent(inout) :: state

      type(scalar_field) :: aux_sfield
      type(mesh_type), pointer :: x_mesh
      real :: current_time

      ! Disgusting and vomitous hack to ensure that time is output in
      ! vtu files.
      x_mesh => extract_mesh(state, "CoordinateMesh")
      call allocate(aux_sfield, x_mesh, "Time", field_type=FIELD_TYPE_CONSTANT)
      call get_option("/timestepping/current_time", current_time)
      call set(aux_sfield, current_time)
      aux_sfield%option_path = ""
      call insert(state, aux_sfield, trim(aux_sfield%name))
      call deallocate(aux_sfield)

    end subroutine insert_time_in_state

    subroutine output_state(state)
      implicit none
      type(state_type), dimension(:), intent(inout) :: state
      !
      type(mesh_type), pointer :: zeta_mesh
      type(scalar_field), pointer :: sfield
      type(scalar_field), dimension(:), allocatable :: BubbleFields
      integer :: stat, n_bubble_fields, n_scalar_fields,i
      ! project the local velocity to cartesian coordinates
      call project_local_to_cartesian(state(1))
      ! Now we're ready to call write_state
      zeta_mesh => extract_mesh(state(1),"VorticityMesh",stat)
      if(stat==0) then
         if(zeta_mesh%shape%numbering%type==ELEMENT_BUBBLE) then
            n_scalar_fields = scalar_field_count(state(1))
            n_bubble_fields = 0
            do i = 1, n_scalar_fields
               sfield => extract_scalar_field(state(1),i)
               if(sfield%mesh==zeta_mesh) then
                  n_bubble_fields = n_bubble_fields+1
               end if
            end do
            allocate(BubbleFields(n_bubble_fields))
            n_bubble_fields = 0
            do i = 1, n_scalar_fields
               sfield => extract_scalar_field(state(1),i)
               if(sfield%mesh==zeta_mesh) then
                  n_bubble_fields = n_bubble_fields+1
                  BubbleFields(n_bubble_fields) = sfield
               end if
            end do
            if(stat==0) then
               call bubble_field_to_vtk(state(1),BubbleFields,&
                    &trim(simulation_name)//"Bubble",dump_no)
            end if
            deallocate(BubbleFields)
         end if
      end if
      call write_state(dump_no, state)
    end subroutine output_state

    subroutine read_command_line()
      implicit none
      ! Read the input filename.

      character(len=1024) :: argument
      integer :: status, argn, level

      call set_global_debug_level(0)

      argn=1
      do

         call get_command_argument(argn, value=argument, status=status)
         argn=argn+1

         if (status/=0) then
            call usage
            stop
         end if

         if (argument=="-v") then
            call get_command_argument(argn, value=argument, status=status)
            argn=argn+1

            if (status/=0) then
               call usage
               stop
            end if

            read(argument, "(i1)", err=666) level
            call set_global_debug_level(level)

            ! Go back to pick up the command line.
            cycle
         end if

         exit
      end do

      call load_options(argument)
      if(.not. have_option("/simulation_name")) goto 666

      return

666   call usage
      stop

    end subroutine read_command_line

    subroutine usage
      implicit none

      write (0,*) "usage: shallow_water [-v n] <options_file>"
      write (0,*) ""
      write (0,*) "-v n sets the verbosity of debugging"
    end subroutine usage


  end program shallow_water

! TODO LIST FOR BDFM1 SWE

! Matrix-ify Helmholtz solver - CODED and TESTED
! Newton iteration for linear equations - CODED and TESTED
! Call to Newton iteration from main code - CODED and TESTED
! Check that timestepping produces some output - DONE and TESTED
! Extract fluxes from DG -- DONE and TESTED
! Vorticity calculation -- DONE and TESTED
! Mass lumping for P2b -- DONE and TESTED
! Visualisation of P2b by mapping back to P2 -- DONE and TESTED
! Mass mapping from P1dg to P2b -- DONE and TESTED
! PV calculation -- DONE
! Timestepping for PV -- DONE, TESTED (advection of 1)
! Computation of PV flux -- DONE, TESTED
! Nonlinear residual calculation from Velocity and DG advection, DONE, TESTED
! Check on spherical mesh, DONE, TESTED
! stabilisation for PV, DONE, TESTED
! discontinuity capturing for PV
! stabilisation for divergence
! Improve DG timestepping
! Improve slope limiter
! Visualisation of vorticity with bubbles, DONE, TESTED
