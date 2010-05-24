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
program test_advection_diffusion
  use spud
  use fields
  use state_module
  use FLDebug
  use populate_state_module
  use write_state_module
  use populate_state_module
  use timeloop_utilities
  use sparsity_patterns_meshes
  use solvers
  use diagnostic_fields_wrapper
  implicit none
#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif

  type(state_type), dimension(:), pointer :: state
  real :: dt
  integer :: timestep
  integer :: ierr

#ifdef HAVE_MPI
  call mpi_init(ierr)
#endif

#ifdef HAVE_PETSC
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

  call python_init()
  call read_command_line()

  call populate_state(state)

  ! No support for multiphase or multimaterial at this stage.
  if (size(state)/=1) then
     FLAbort("Multiple material_phases are not supported")
  end if

  ! Always output the initial conditions.
  call output_state(state)

  call get_option("/timestepping/current_time", current_time)
  call get_option("/timestepping/timestep", dt)

  timestep=0
  timestep_loop: do 
     timestep=timestep+1
     ewrite (1,'(a,i0)') "Start of timestep ",timestep
     
     call execute_timestep(state(1), dt)

     call calculate_diagnostic_variables(state)

     if (simulation_completed(current_time, timestep)) exit timestep_loop     

     call advance_current_time(current_time, dt)

     if (do_write_state(current_time, timestep)) then
        call output_state(state)
     end if

  end do timestep_loop
  
  ! One last dump
  call output_state(state)

#ifdef HAVE_MPI
  call mpi_finalize(ierr)
#endif

contains

  subroutine execute_timestep(state, dt)
    type(state_type), intent(inout) :: state
    logical :: have_advection, have_diffusion
    integer :: stat
    real, intent(in) :: dt

    !! The tracer we're solving for.
    type(scalar_field), pointer :: T
    !! Advecting velocity.
    type(vector_field), pointer :: U
    !! Coordinate field.
    type(vector_field), pointer :: X
    !! Diffusivity of the tracer.
    type(tensor_field), pointer :: mu

    !! Intermediate results.
    type(scalar_field), dimension(4) :: T_RK
    type(scalar_field) :: rhs, T_update
    real, dimension(4) :: rk_coefs=(/1./6.,1./3.,1./3.,1./6./)

    !! Sparsity for matrices.
    type(csr_sparsity) :: sparsity
    !! Mass matrix
    type(csr_matrix) :: Mass
    !! Diffusion matrix
    type(csr_matrix) :: Diffusion_mat

    integer :: i, ele

    T=>extract_scalar_field(state, "Tracer")

    X=>extract_vector_field(state, "Coordinate")

    U=>extract_vector_field(state, "Velocity", stat=stat)
    have_advection=(stat==0)

    mu=>extract_tensor_field(state, "TracerDiffusivity", stat=stat)
    have_diffusion=(stat==0)

    if ((.not.have_advection).and.(.not.have_diffusion)) then
       FLAbort("Running with no advection or diffusion is pointless")
    end if

    sparsity=get_csr_sparsity_firstorder(state, T%mesh, T%mesh)

    call allocate(rhs, T%mesh, "RHS")
    call zero(rhs)

    !! Advection step.
    if (have_advection) then
       call allocate(Mass, sparsity, name="Mass")

       do i=1,4
          call allocate(T_RK(i), T%mesh, "Tracer_RK")
          ! Ensure T_RK inherits solver options from T.
          T_RK(i)%option_path=T%option_path
       end do
       
       call rk_advection_step(T_RK(1), Mass, dt, X, U, T)
       
       call rk_advection_step(T_RK(2), Mass, dt, X, U, T, 0.5, T_RK(1))
       
       call rk_advection_step(T_RK(3), Mass, dt, X, U, T, 0.5, T_RK(2))
       
       call rk_advection_step(T_RK(4), Mass, dt, X, U, T, 1.0, T_RK(3))
       
       call allocate(T_update, T%mesh, "T_update")
       call zero(T_update)

       do i=1,4
          call addto(T_update, T_RK(i), scale=rk_coefs(i))

          ewrite_minmax(T_RK(i)%val)
          
          call deallocate(T_RK(i))
       end do
    
       call addto(T, T_update)
       
       call deallocate(T_update)
       call deallocate(Mass)
    end if

    !! Diffusion step.
    if (have_diffusion) then
       call zero(rhs)
       
       call allocate(Diffusion_mat, sparsity, name="DiffusionMatrix")
       call zero(Diffusion_mat)
       
       do ele=1,element_count(T)
          call construct_diffusion_ele(ele, Diffusion_mat, rhs, dt, X,&
               & mu, T)
       end do
       
       call petsc_solve(T, Diffusion_mat, rhs)
       
       call deallocate(Diffusion_mat)

    end if
       
    call deallocate(rhs)

  end subroutine execute_timestep

  subroutine construct_diffusion_ele(ele, Diffusion_mat, rhs, dt, X, mu, T)
    integer, intent(in) :: ele
    type(csr_matrix), intent(inout) :: Diffusion_mat
    type(scalar_field), intent(inout) :: rhs
    real, intent(in) :: dt
    type(vector_field), intent(in) :: X
    type(tensor_field), intent(in) :: mu
    type(scalar_field), intent(in) :: T

    real, dimension(ele_ngi(T, ele)) :: detwei
    real, dimension(ele_loc(T, ele), ele_ngi(T, ele), mesh_dim(T)) :: dt_t
    ! Diffusivity at quadrature points.
    real, dimension(mu%dim, mu%dim, ele_ngi(mu, ele)) :: mu_q
    integer, dimension(:), pointer :: T_ele
    real, dimension(ele_loc(T,ele)) :: T_val
    real, dimension(ele_loc(T, ele), ele_loc(T, ele)) :: mass, diffusion

    type(element_type) :: t_shape

    t_shape=ele_shape(T,ele)
    
    call transform_to_physical(X,ele,&
         & t_shape , dshape=dt_t, detwei=detwei)
    
    T_val=ele_val(T, ele)
    T_ele=>ele_nodes(T, ele)
    mu_q=ele_val_at_quad(mu, ele)

    mass=shape_shape(t_shape, t_shape, detwei)
    diffusion=-dt*dshape_tensor_dshape(dt_t, mu_q, dt_t, detwei)

    call addto(Diffusion_mat, T_ele, T_ele, mass-0.5*diffusion)

    call addto(rhs, T_ele, matmul(mass+0.5*diffusion, T_val))

  end subroutine construct_diffusion_ele

  subroutine rk_advection_step(T_RK, Mass, dt, X, U, T, T_scale,&
       & T_change)
    type(scalar_field), intent(inout) :: T_RK
    type(csr_matrix), intent(inout) :: Mass
    real, intent(in) :: dt
    type(vector_field), intent(in) :: X, U
    type(scalar_field), intent(in) :: T
    real, intent(in), optional :: T_scale
    type(scalar_field), intent(in), optional :: T_change

    type(scalar_field) :: rhs
    integer :: ele

    call allocate(rhs, T_RK%mesh, "RHS_RK")
    call zero(rhs)
    call zero(T_RK)
    call zero(Mass)
    
    do ele=1, element_count(T_RK)
       call rk_advection_step_ele(ele, Mass, rhs, dt, X, U, T, T_scale,&
       & T_change)
    end do

    ewrite_minmax(rhs)
    call petsc_solve(T_RK, Mass, rhs)

    ewrite_minmax(T_RK)
    call deallocate(rhs)
    
  end subroutine rk_advection_step

  subroutine rk_advection_step_ele(ele, Mass, rhs, dt, X, U, T, T_scale,&
       & T_change)
    integer, intent(in) :: ele
    type(csr_matrix), intent(inout) :: Mass
    type(scalar_field), intent(inout) :: rhs
    real, intent(in) :: dt
    type(vector_field), intent(in) :: X, U
    type(scalar_field), intent(in) :: T
    real, intent(in), optional :: T_scale
    type(scalar_field), intent(in), optional :: T_change

    real, dimension(ele_ngi(T, ele)) :: detwei
    real, dimension(ele_loc(T, ele), ele_ngi(T, ele), mesh_dim(T)) :: dt_t
    ! Velocity at quadrature points.
    real, dimension(U%dim, ele_ngi(U, ele)) :: u_q
    real, dimension(ele_ngi(U, ele)) :: u_div_q
    real, dimension(ele_loc(T,ele)) :: T_val
    integer, dimension(:), pointer :: T_ele

    type(element_type) :: t_shape

    t_shape=ele_shape(T, ele)
    
    call transform_to_physical(X,ele,&
         & t_shape , dshape=dt_t, detwei=detwei)
    
    T_val=ele_val(T, ele)
    if (present(T_change)) then
       T_val=T_val + T_scale * ele_val(T_change, ele)
    end if

    ! Advecting velocity at quadrature points.
    U_q=ele_val_at_quad(U,ele)
    U_div_q=ele_div_at_quad(U,ele,dt_t)

    T_ele=>ele_nodes(T, ele)

    call addto(rhs, T_ele, &
         dt* matmul((dshape_dot_vector_shape(dt_t, U_q, t_shape, detwei)&
         &-shape_shape(t_shape, t_shape, U_div_q*detwei)&
         ),          T_val) &
               )
    
    call addto(Mass, T_ele, T_ele, shape_shape(t_shape, t_shape, detwei))
!    call addto_diag(Mass, T_ele, sum(shape_shape(t_shape, t_shape, detwei),1))

  end subroutine rk_advection_step_ele

  subroutine advance_current_time(current_time, dt)
    real, intent(inout) :: current_time, dt
    
    ! Adaptive timestepping could go here.

    current_time=current_time + dt

  end subroutine advance_current_time

  subroutine output_state(state)
    type(state_type), dimension(:), intent(inout) :: state

    integer, save :: dump_no=0

    call write_state(dump_no, state)
    
  end subroutine output_state

  subroutine read_command_line()
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

          ! Go back to picj up the command line.
          cycle
       end if

       exit
    end do

    call load_options(argument)

    return

666 call usage
    stop

  end subroutine read_command_line

  subroutine usage
    
    write (0,*) "usage: test_advection_diffusion [-v n] <options_file>"
    write (0,*) ""
    write (0,*) "-v n sets the verbosity of debugging"
  end subroutine usage

end program test_advection_diffusion
