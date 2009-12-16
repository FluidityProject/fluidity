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
  program shallow_water
    use advection_diffusion_dg
    use advection_diffusion_cg
    use linear_shallow_water
    use spud
    use fields
    use state_module
    use FLDebug
    use populate_state_module
    use write_state_module
    use populate_state_module
    use timeloop_utilities
    use sparsity_patterns_meshes
    use sparse_matrices_fields
    use solvers
    use diagnostic_variables
    use diagnostic_fields_wrapper
    use assemble_cmc
    use global_parameters, only: option_path_len
    implicit none
#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif

    type(state_type), dimension(:), pointer :: state
    real :: current_time, dt, f0, D0, g, theta, itheta
    logical :: exclude_velocity_advection, exclude_pressure_advection
    real, dimension(:), allocatable :: beta
    integer :: timestep, nonlinear_iterations
    integer :: ierr

    !! Sparsity for matrices.
    type(csr_sparsity) :: ct_sparsity,u_sparsity,wave_sparsity

    !! Mass matrices
    type(csr_matrix) :: h_mass_mat, u_mass_mat
    !! Coriolis matrix
    type(csr_matrix) :: coriolis_mat
    !! div matrix
    type(block_csr_matrix) :: div_mat
    !! Wave matrix
    type(csr_matrix) :: wave_mat
    !! U momentum matrix 
    type(block_csr_matrix) :: big_mat
    character(len = OPTION_PATH_LEN) :: simulation_name
#ifdef HAVE_MPI
    call mpi_init(ierr)
#endif

#ifdef HAVE_PETSC
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

    call python_init()
    call read_command_line()

    call populate_state(state)

    call get_option('/simulation_name',simulation_name)
    Call initialise_diagnostics(trim(simulation_name),state)

    call get_parameters()

    ! No support for multiphase or multimaterial at this stage.
    if (size(state)/=1) then
       FLAbort("Multiple material_phases are not supported")
    end if

    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/timestep", dt)
    call setup_wave_matrices(state(1),u_sparsity,wave_sparsity,ct_sparsity, &
         h_mass_mat,u_mass_mat,coriolis_mat,div_mat,wave_mat,big_mat, &
         dt,theta,D0,g,f0,beta)
    call get_option("/timestepping/nonlinear_iterations"&
         &,nonlinear_iterations)
    exclude_pressure_advection = &
         have_option("/material_phase::Fluid/scalar_field::LayerThickness/pro&
         &gnostic/spatial_discretisation/continuous_galerkin/advection_terms&
         &/exclude_advection_terms")
    exclude_velocity_advection = &
         have_option("/material_phase::Fluid/vector_field::Velocity/prognost&
         &ic/spatial_discretisation/discontinuous_galerkin/advection_scheme/&
         &none")
    call get_option("/material_phase::Fluid/scalar_field::LayerThickness/pro&
         &gnostic/relaxation",itheta)
    ! Geostrophic balanced initial condition, if required    
    if(have_option("/material_phase::Fluid/vector_field::Velocity/prognostic&
         &/initial_condition::WholeMesh/balanced")) then       
       call set_velocity_from_geostrophic_balance(state(1), &
            div_mat,coriolis_mat)
    end if

    ! Always output the initial conditions.
    call output_state(state)

    timestep=0
    timestep_loop: do 
       timestep=timestep+1
       ewrite (1,*) "SW: start of timestep ",timestep, current_time

       call execute_timestep(state(1), dt)

       call calculate_diagnostic_variables(state)

       if (simulation_completed(current_time, timestep)) exit timestep_loop     

       call advance_current_time(current_time, dt)

       if (do_write_state(current_time, timestep)) then
          call output_state(state)
       end if

       call write_diagnostics(state,timestep*dt,dt)

    end do timestep_loop

    ! One last dump
    call output_state(state)

#ifdef HAVE_MPI
    call mpi_finalize(ierr)
#endif

  contains

    subroutine get_parameters()
      implicit none
      !Get some parameters
      !Coriolis
      allocate(beta(2))
      if(have_option("/physical_parameters/coriolis")) then
         if(have_option("/physical_parameters/coriolis/f_plane")) then
            call get_option("/physical_parameters/coriolis/f_plane/f",f0)
            beta = 0.0
         else if(have_option("/physical_parameters/coriolis/beta_plane")) then
            call get_option("/physical_parameters/coriolis/f_plane/f_0",f0)
            call get_option("/physical_parameters/coriolis/f_plane/beta",beta)
         else
            FLAbort('Your chosen Coriolis option is not supported')
         end if
      else
         f0 = 0.
         beta = 0.
      end if
      !gravity
      call get_option("/physical_parameters/gravity/magnitude",g)
      !theta
      call get_option("/material_phase::Fluid/scalar_field::LayerThickness/prognostic/temporal_discretisation/theta",theta)
      !itheta
      
      !D0
      call get_option("/material_phase::Fluid/scalar_field::LayerThickness/p&
           &rognostic/mean_layer_thickness",D0)
    end subroutine get_parameters

    subroutine execute_timestep(state, dt)
      implicit none
      type(state_type), intent(inout) :: state
      real, intent(in) :: dt

      !! Layer thickness
      type(scalar_field), pointer :: D
      !! velocity.
      type(vector_field), pointer :: U
      !! Coordinate field.
      type(vector_field), pointer :: X

      !!Intermediate fields
      type(scalar_field) :: d_rhs, delta_d, old_d
      type(vector_field) :: u_rhs, delta_u, advecting_u, old_u
      type(scalar_field) :: velocity_cpt, old_velocity_cpt
      integer :: i,  dim, nit, d1

      !Pull the fields out of state
      D=>extract_scalar_field(state, "LayerThickness")
      X=>extract_vector_field(state, "Coordinate")
      U=>extract_vector_field(state, "Velocity")

      dim = mesh_dim(u)

      call execute_timestep_setup(D,X,U,d_rhs,u_rhs,advecting_u, &
           old_u,old_d,delta_d,delta_u,state)

      call insert(state,advecting_u,"NonlinearVelocity")

      call set(advecting_u,u)
      call set(old_u,u)
      call set(old_d,d)

      do nit = 1, nonlinear_iterations

         call set(u,old_u)
         call set(d,old_d)

         !velocity advection step
         if(.not.exclude_velocity_advection) then
            do d1 = 1, dim
               velocity_cpt = extract_scalar_field(u,d1)
               old_velocity_cpt = extract_scalar_field(old_u,d1)
               call insert(state,velocity_cpt,"VelocityComponent")
               call insert(state,old_velocity_cpt,"VelocityComponentOld")
               call solve_advection_diffusion_dg("VelocityComponent", state)
            end do
         end if
         !pressure advection
         if(.not.exclude_pressure_advection) then
            call solve_field_equation_cg("LayerThickness", state, dt, &
                 "NonlinearVelocity")
         end if
         
         !Wave equation step
         ! M\Delta u + \Delta t F(\theta\Delta u + u^n) + \Delta t C(\theta
         ! \Delta\eta + \eta^n) = 0
         ! M\Delta h - \Delta t HC^T(\theta\Delta u + u^n) = 0
         ! SO
         ! (M+\theta\Delta t F)\Delta u = 
         !   -\theta\Delta t g C\Delta\eta + \Delta t(-Fu^n - gC\eta^n)
         ! SET r = \Delta t(-Fu^n - gC\eta^n)
         ! THEN SUBSTITUTION GIVES
         ! (M + \theta^2\Delta t^2 gH C^T(M+\theta\Delta t F)^{-1}C)\Delta\eta
         !  = \Delta t HC^T(u^n + \theta(M+\theta\Delta t F)^{-1}r)
         
         !Construct explicit parts of u rhs
         call get_u_rhs(u_rhs,U,D,dt,g, &
              coriolis_mat,div_mat)
         
         !Construct explicit parts of h rhs in wave equation
         call get_d_rhs(d_rhs,u_rhs,D,U,h_mass_mat,div_mat,big_mat,D0,dt,theta)
         !Solve wave equation for D update
         delta_d%option_path = d%option_path
         call petsc_solve(delta_d, wave_mat, d_rhs)
         
         !Add the new D contributions into the RHS for u
         call update_u_rhs(u_rhs,U,delta_D,div_mat,theta,dt,g)
         
         !Solve momentum equation for U update
         call zero(delta_u)
         call mult(delta_u, big_mat, u_rhs)
         
         !Check the equation was solved correctly
         !call check_solution(delta_u,delta_d,d,u,dt,theta,g,D0,u_mass_mat&
         !     &,h_mass_mat, coriolis_mat,div_mat)

         call addto(u,delta_u)
         call addto(d,delta_d)

         call set(advecting_u,old_u)
         call scale(advecting_u,(1-itheta))
         call addto(advecting_u,u,scale=itheta)

      end do

      !Update the variables


      call get_energy(u,d,d0,g,u_mass_mat,h_mass_mat)

      call deallocate(d_rhs)
      call deallocate(u_rhs)
      call deallocate(delta_d)
      call deallocate(delta_u)
      call deallocate(old_u)
      call deallocate(old_d)
      call deallocate(advecting_u)

    end subroutine execute_timestep

    subroutine get_energy(u,d,d0,g,u_mass_mat,h_mass_mat)
      type(scalar_field), intent(inout) :: d
      type(vector_field), intent(inout) :: u
      real, intent(in) :: d0,g
      type(csr_matrix), intent(in) :: u_mass_mat,h_mass_mat
      !
      type(scalar_field) :: Md
      type(vector_field) :: Mu
      real :: D_l2,u_l2, energy
      integer :: d1, dim

      dim = mesh_dim(U)
      call allocate(Md,D%mesh,'Md')
      call allocate(Mu,mesh_dim(U),u%mesh,'Mu')
      !
      call mult(Md,h_mass_mat,D)
      call mult(Mu,u_mass_mat,U)
      D_l2 = sum(Md%val*d%val) 
      U_l2 = 0.
      do d1 = 1, dim
         U_l2 = U_l2 + sum(Mu%val(d1)%ptr*u%val(d1)%ptr)
      end do
      energy = 0.5*g*D_l2 + 0.5*D0*U_l2
      ewrite(2,*) 'SW: energy = ', energy
      !
      call deallocate(Md)
      call deallocate(Mu)
    end subroutine get_energy

    subroutine execute_timestep_setup(D,X,U,d_rhs,u_rhs,advecting_u, &
         old_u,old_d,delta_d,delta_u,state)
      implicit none
      type(state_type), intent(inout) :: state
      type(scalar_field), intent(inout), pointer :: D
      type(scalar_field), intent(inout) :: D_rhs, delta_d, old_d
      type(vector_field), intent(inout), pointer :: U,X
      type(vector_field), intent(inout) :: U_rhs, delta_u, advecting_u, old_u
      integer :: dim

      dim = mesh_dim(D)

      !allocate RHS variables
      call allocate(d_rhs, D%mesh, "LayerThickness_RHS")
      call zero(d_rhs)
      call allocate(u_rhs, dim, U%mesh, "Velocity_RHS")
      call zero(u_rhs)
      !allocate update variables
      call allocate(delta_d, D%mesh, "LayerThickness_update")
      call zero(delta_d)
      call allocate(delta_u, dim, U%mesh, "Velocity_update")
      call zero(delta_u)
      !allocate advecting velocity
      call allocate(advecting_u,dim,U%mesh, "Advecting_velocity")
      call zero(advecting_u)
      !allocate previous timestep variables
      call allocate(old_d, D%mesh, "LayerThickness_old")
      call zero(old_d)
      call allocate(old_u, dim, U%mesh, "Velocity_old")
      call zero(old_u)

    end subroutine execute_timestep_setup

    subroutine advance_current_time(current_time, dt)
      implicit none
      real, intent(inout) :: current_time, dt

      ! Adaptive timestepping could go here.

      current_time=current_time + dt

    end subroutine advance_current_time

    subroutine output_state(state)
      implicit none
      type(state_type), dimension(:), intent(inout) :: state

      integer, save :: dump_no=0

      call write_state(dump_no, state)

    end subroutine output_state

    subroutine set_velocity_from_geostrophic_balance(state,div_mat,coriolis_mat)
      implicit none
      type(state_type), intent(inout) :: state
      type(block_csr_matrix), intent(in) :: div_mat
      type(csr_matrix), intent(in) :: coriolis_mat
      !
      type(scalar_field), pointer :: D
      type(vector_field), pointer :: U,X
      integer :: dim, ele
      type(vector_field) :: u_tmp
      type(scalar_field) :: u1,u2

      ewrite(1,*) '    subroutine set_velocity_from_geostrophic_balance()'

      !Pull the fields out of state
      D=>extract_scalar_field(state, "LayerThickness")
      U=>extract_vector_field(state, "Velocity")
      X=>extract_vector_field(state, "Coordinate")
      
      ewrite(2,*) 'inside', sum(D%val)

      u1=extract_scalar_field(u,1)
      u2=extract_scalar_field(u,2)

      call allocate(u_tmp,mesh_dim(U),U%mesh,name="tempmem")

      call mult_T(u_tmp,div_mat,D)

      call set(u1,u_tmp,2)
      call set(u2,u_tmp,1)
      call scale(u1,-g)
      call scale(u2,g)

      do ele = 1, ele_count(U)
         call invert_coriolis(U,X,ele)
      end do

      call get_u_rhs(u_tmp,U,D,dt,g, &
         coriolis_mat,div_mat)

      ewrite(3,*) 'SW: TESTING BALANCED INITIAL CONDITION'
      ewrite(3,*) 'SW: infty norm of u_tmp(1)=', maxval(abs(u_tmp%val(1)%ptr))
      ewrite(3,*) 'SW: infty norm of u_tmp(2)=', maxval(abs(u_tmp%val(2)%ptr))
      ewrite(3,*) 'SW: TESTING BALANCED INITIAL CONDITION'

      call deallocate(u_tmp)

      ewrite(1,*) 'END subroutine set_velocity_from_geostrophic_balance()'

    end subroutine set_velocity_from_geostrophic_balance

    subroutine invert_coriolis(U,X,ele)
      type(vector_field), intent(inout) :: U
      type(vector_field), intent(in) :: X
      integer, intent(in) :: ele
      !
      real, dimension(ele_ngi(U, ele)) :: detwei
      integer, dimension(:), pointer :: U_ele
      type(element_type) :: u_shape
      real, dimension(ele_loc(U,ele),ele_loc(U,ele)) :: c_mat
      integer :: dim
      real, dimension(ele_ngi(U, ele)) :: f_gi
      real, dimension(mesh_dim(U), ele_ngi(x,ele)) :: x_gi
      real, dimension(mesh_dim(U), ele_loc(u,ele)) :: U_val
      integer :: d1

      dim = mesh_dim(U)

      x_gi = ele_val_at_quad(X,ele)
      f_gi = f0 + matmul(beta,x_gi)

      u_shape=ele_shape(u, ele)
      U_ele => ele_nodes(U, ele)
      U_val = ele_val(U, ele)

      call transform_to_physical(X,ele,detwei=detwei)
      
      c_mat = shape_shape(u_shape,u_shape,detwei*f_gi)
      call invert(c_mat)
      
      do d1 = 1, dim
         call set(U,d1,u_ele,matmul(c_mat,U_val(d1,:)))
      end do
      
    end subroutine invert_coriolis

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
