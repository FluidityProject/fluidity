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

module shallow_water
  !!< This module solves the linear shallow water equations. 
  use spud
  use global_parameters, only: FIELD_NAME_LEN
  use state_module
  use boundary_conditions_from_options
  use populate_state_module
  use sparse_matrices_fields
  use diagnostic_variables
  use diagnostic_fields_wrapper
  use sparsity_patterns
  use momentum_dg
  use fetools
  use solvers
  implicit none

  private
  
  public :: shallow_water_icom

  !! Implicitness and timestep
  real, save :: theta, dt

  !! Current and finish time
  real, save :: current_time, finish_time

  !! Index of the current timestep
  integer, save :: timestep=0

  !! Index of the current output dump
  integer, save :: dump_no=0

  !! Time at which next dump should occur
  real, save :: dump_time
  !! Period between dumps in time units
  real, save :: dump_period

  !! Timestep at which next dump should occur
  integer, save :: dump_step
  !! Period between dumps in timesteps
  integer, save :: dump_period_steps

  !! Acceleration due to gravity.
  real, save :: gravity=9.8

contains

  subroutine shallow_water_icom
    type(state_type), dimension(:), pointer :: state
    character(len=FIELD_NAME_LEN) :: simulation_name

    call get_option("/simulation_name", simulation_name)

    ! Read in and populate all the fields.
    call populate_state(state)
    call allocate_shallow_water_rhs(state(1))

    ! Set up the matrices we'll need.

    ! Generic fluidity sparsity routines.
    call form_matrix_sparsity(state)
    ! Extra sparsity patterns needed by shallow water.
    call form_shallow_water_sparsity(state(1))
    ! Allocate corresponding matrices.
    call allocate_shallow_water_matrices(state(1))

    ! Set up the .stat file for output.
    call initialise_diagnostics(simulation_name, state)

    ! Extract some essential values
    call set_global_variables

    call assemble_shallow_water_equations(state(1))

    ! Time zero output
    call calculate_diagnostic_variables(state)
    call write_diagnostics(simulation_name, state(1))
    ! Don't dump now for zero length simulations: this will be done during
    ! cleanup. 
    if (current_time<finish_time) then
       call write_state(dump_no, state)
       dump_no=dump_no+1
    end if

    call shallow_water_timestep_loop(state)
    
    ! Final output
    call write_state(dump_no, state)

    ! Call 
    !call cleanup_shallow_water(state)

  end subroutine shallow_water_icom
 

  subroutine set_global_variables
    
    call get_option("/material_phase[0]/vector_field::Velocity"//&
         &"/prognostic/temporal_discretisation/theta", theta)
    call get_option("/timestepping/timestep", dt)
    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/finish_time", finish_time)
    if (have_option("/io/dump_period")) then
       call get_option("/io/dump_period", dump_period)
       dump_time=dump_period
       dump_step=huge(0)
    else if (have_option("/io/dump_period_in_timesteps") then
       call get_option("/io/dump_period_in_timesteps", dump_period_steps)
       dump_step=dump_period_steps
       dump_time=huge(0.0)
    end if

  end subroutine set_global_variables

  
  subroutine shallow_water_timestep_loop(state)
    type(state_type), dimension(:), intent(inout) :: state
    real, dimension(:), pointer :: rhs=>null()
    
    timeloop: do
       if (current_time>=finish_time) exit timeloop
       timestep=timestep+1
       
       call assemble_rhs_shallow_water(state(1), rhs)
       
       call solve_shallow_water_equation(state(1), rhs)

       call shallow_water_advance_time(state)
    end do timeloop
  
  end subroutine shallow_water_timestep_loop


  subroutine assemble_shallow_water_equations(state)
    type(state_type), intent(inout) :: state
    type(block_csr_matrix) :: U_mat, CT
    type(csr_matrix) ::  Hmass, Umass
    type(vector_field) :: U_rhs
    type(scalar_field) :: h_rhs
    type(vector_field) :: X
    type(scalar_field) :: H
    type(scalar_field), pointer :: p, rho
    type(vector_field), pointer :: u, x_coord

    ! get the velocity
    u=>extract_vector_field(state, "Velocity")
    x_coord=>extract_vector_field(state, "Coordinate")
    p=>extract_scalar_field(state, "Pressure")
    rho=>extract_scalar_field(state, "Density")

    U_mat=extract_block_csr_matrix(state, "CoupledMomentumMatrix")
    CT=extract_block_csr_matrix(state, "PressureGradientSparsity")
    Hmass=extract_csr_matrix(state, "FreeSurfaceMatrix")
    Umass=extract_csr_matrix(state, "VelocityMassMatrix")

    U_rhs=extract_vector_field(state, "VelocityRHS")
    H_rhs=extract_scalar_field(state, "FreeSurfaceRHS")
    X=extract_vector_field(state, "Position")
    H=extract_scalar_field(state, "FreeSurface")

    call zero(U_mat)
    call zero(CT)
    call zero(Hmass)
    call zero(Umass)

    call zero(U_rhs)
    call zero(H_rhs)
        
    call construct_momentum_dg(u, p, rho, x_coord, U_mat, CT, U_rhs, H_rhs, &
         state=state, mass=Umass, acceleration_form=.false.)

    call assemble_mass(Hmass, H, X)
    
    call assemble_full_matrix(state, U_mat, CT, Umass, Hmass)

  end subroutine assemble_shallow_water_equations

  subroutine assemble_rhs_shallow_water(state, rhs)
    !!< Generate the RHS of the equations. Note that RHS is now a vector
    !!< not a field as femtools does not yet support mixed fields.
    type(state_type), intent(in) :: state
    real, dimension(:), pointer :: rhs

    real, dimension(:), pointer :: u_rhs, h_rhs
    integer :: u_nodes, h_nodes

    type(csr_matrix) :: Hmass, Umass
    type(block_csr_matrix) :: U_mat, CT
    type(scalar_field) :: H_boundary, h_tmp, H
    type(vector_field) :: U_boundary, u_tmp, U

    U_mat=extract_block_csr_matrix(state, "CoupledMomentumMatrix")
    CT=extract_block_csr_matrix(state, "PressureGradientSparsity")
    Hmass=extract_csr_matrix(state, "FreeSurfaceMatrix")
    Umass=extract_csr_matrix(state, "VelocityMassMatrix")

    U=extract_vector_field(state, "Velocity")
    H=extract_scalar_field(state, "FreeSurface")
    U_boundary=extract_vector_field(state, "VelocityRHS")
    H_boundary=extract_scalar_field(state, "FreeSurfaceRHS")

    u_nodes=node_count(U_boundary)
    h_nodes=node_count(H_boundary)

    if (.not.associated(rhs)) then
       allocate(rhs(u_nodes*2+h_nodes))
    end if
    rhs=0.0

    u_rhs=>rhs(:u_nodes)
    h_rhs=>rhs(u_nodes+1:)

    call allocate(u_tmp, 2, U_boundary%mesh, "VelocityTmp")
    call allocate(h_tmp, H_boundary%mesh, "FreeSurfaceTmp")

    ! Boundary values
    u_rhs(:u_nodes/2)=U_boundary%val(1)%ptr
    u_rhs(u_nodes/2+1:)=U_boundary%val(2)%ptr
    h_rhs=H_boundary%val

    ! Coriolis and advection matrix.
    call mult(u_tmp, U_mat, U)
    u_rhs(:u_nodes/2)=u_rhs(:u_nodes/2)+dt*(1-theta)*u_tmp%val(1)%ptr
    u_rhs(u_nodes/2+1:)=u_rhs(u_nodes/2+1:)+dt*(1-theta)*u_tmp%val(2)%ptr

    ! Velocity mass.
    call mult(u_tmp%val(1)%ptr, Umass, U%val(1)%ptr)
    call mult(u_tmp%val(2)%ptr, Umass, U%val(2)%ptr)
    u_rhs(:u_nodes/2)=u_rhs(:u_nodes/2)+u_tmp%val(1)%ptr
    u_rhs(u_nodes/2+1:)=u_rhs(u_nodes/2+1:)+u_tmp%val(2)%ptr
    
    ! Pressure gradient.
    call mult_T(u_tmp, CT, H)
    u_rhs(:u_nodes/2)=u_rhs(:u_nodes/2)+dt*(1-theta)*u_tmp%val(1)%ptr&
         &*gravity
    u_rhs(u_nodes/2+1:)=u_rhs(u_nodes/2+1:)+dt*(1-theta)*u_tmp%val(2)%ptr&
         &*gravity

    ! Free surface mass
    call mult(h_tmp, Hmass, H)
    h_rhs=h_rhs+h_tmp%val

    ! Velocity divergence.
    call mult(h_tmp, CT, u)
    h_rhs=h_rhs+dt*(1-theta)*h_tmp%val
    
    call deallocate(u_tmp)
    call deallocate(h_tmp)
    
  end subroutine assemble_rhs_shallow_water


  subroutine solve_shallow_water_equation(state, rhs)
    !!< Solve the equations previously assembled.
    type(state_type), intent(in) :: state
    real, dimension(:), pointer :: rhs

    type(csr_matrix) :: Big_m
    real, dimension(:), allocatable :: x

    type(vector_field) :: U
    type(scalar_field) :: H
    integer :: H_nodes, U_nodes

    U=extract_vector_field(state, "Velocity")
    H=extract_scalar_field(state, "FreeSurface")
    Big_m=extract_csr_matrix(state, "FullSystemMatrix")

    u_nodes=node_count(U)
    h_nodes=node_count(H)
    
    allocate(x(size(rhs)))

    x(1:u_nodes)=U%val(1)%ptr
    x(u_nodes+1:2*u_nodes)=U%val(2)%ptr
    x(2*u_nodes+1:)=H%val

    call petsc_solve(x, Big_m, rhs, U%option_path)

    U%val(1)%ptr=x(1:u_nodes)
    U%val(2)%ptr=x(u_nodes+1:2*u_nodes)
    H%val=x(2*u_nodes+1:)

  end subroutine solve_shallow_water_equation

  subroutine shallow_water_output(state)
    !!< Output diagnostics and if it's time for a dump, do that too.
    
    call calculate_diagnostic_variables(state)
    call write_diagnostics(simulation_name, state(1))
        
    if (current_time>=dump_time)

  end subroutine shallow_water_output

  subroutine shallow_water_advance_time(state)
    !!< Move forward the timestep and update prescribed fields and boundary
    !!< values.
    type(state_type), dimension(:), intent(inout) :: state
    
    current_time=current_time+dt
    call set_option("/timestepping/current_time", current_time)

    call set_prescribed_field_values(state)
    call set_boundary_conditions_values(state)
    ! if strong bc or weak that overwrite then enforce the bc on the fields
    call set_dirichlet_consistent(state)

  end subroutine shallow_water_advance_time
  

  subroutine assemble_full_matrix(state, U_mat, CT, Umass, Hmass)
    !!< Assemble the full left hand side matrix for the shallow water
    !!< equations.
    type(state_type), intent(inout) :: state
    type(block_csr_matrix), intent(in) :: U_mat, CT
    type(csr_matrix), intent(in) :: Umass, Hmass

    type(dynamic_csr_matrix) :: dbig_m
    type(csr_matrix) :: big_m
    
    integer :: total_rows, u_rows, h_rows, row, u_nodes
    integer :: offset1, offset2, dim1, dim2
    integer, dimension(:), pointer :: row_m

    u_nodes=block_size(U_mat,1)
    u_rows=size(U_mat,1)
    h_rows=size(Hmass,1)
    total_rows=u_rows+h_rows

    call allocate(dbig_m, total_rows, total_rows)
    
    do row = 1, u_nodes
       do dim1 = 1, 2
          offset1=(dim1-1)*u_nodes
          
          row_m => row_m_ptr(U_mat,row)          

          ! Coriolis(, advection) and diffusion part.
          do dim2= 1, 2
             offset2=(dim2-1)*u_nodes

             call addto(dbig_m, row+offset1, row_m+offset2, &
                  -1*row_val_ptr(U_mat,dim1,dim2,row)*theta*dt)
          end do

          ! Velocity mass.
          row_m => row_m_ptr(Umass,row)          
          
          call addto(dbig_m, row+offset1, row_m+offset1, &
                  row_val_ptr(Umass,row))
       end do
    end do
    
    ! Pressure gradient and div u
    do row = 1, h_rows

       row_m => row_m_ptr(CT, row)

       do dim2 = 1, 2
          offset2=(dim2-1)*u_nodes
          
          ! Pressure gradient
          call addto(dbig_m, row+u_rows, row_m+offset2, &
               -1*row_val_ptr(CT, 1, dim2, row)*gravity*theta*dt)

          ! Div u
          call addto(dbig_m, row_m+offset2, row+u_rows,  &
               -1*row_val_ptr(CT, 1, dim2, row)*theta*dt)

       end do

       ! dH/dt term.
       row_m => row_m_ptr(Hmass, row)
       
       call addto(dbig_m, row+u_rows, row_m+u_rows, &
            row_val(Hmass, row))

    end do

    big_m=dcsr2csr(dbig_m)

    call deallocate(dbig_m)

    call insert(state, big_m, "FullSystemMatrix") 
    ! Drop the extra reference.
    call deallocate(big_m)

  end subroutine assemble_full_matrix


  subroutine assemble_mass(mass, field, positions)
    !!< Contstruct the mass matrix implied by field.
    !!< This routine should be somewhere more general.
    type(csr_matrix), intent(inout) :: mass
    type(scalar_field), intent(in) :: field
    type(vector_field), intent(in) :: positions

    integer :: ele

    do ele = 1, element_count(field)
       
       call assemble_mass_element(mass, field, positions, ele)

    end do

  contains

    subroutine assemble_mass_element(mass, field, positions, ele)
      type(csr_matrix), intent(inout) :: mass
      type(scalar_field), intent(in) :: field
      type(vector_field), intent(in) :: positions
      integer, intent(in) :: ele

      real, dimension(ele_ngi(field,ele)) :: detwei
      integer, dimension(:), pointer :: field_ele
      type(element_type), pointer :: field_shape

      call transform_to_physical(ele_val(positions, ele),&
           & ele_shape(positions, ele), detwei=detwei)

      field_shape=>ele_shape(field, ele)
      field_ele=>ele_nodes(field, ele)

      call addto(mass, field_ele, field_ele, &
           shape_shape(field_shape, field_shape, detwei))

    end subroutine assemble_mass_element

  end subroutine assemble_mass
  

  subroutine form_shallow_water_sparsity(state)
    !!< Construct those sparsities which are not general to fluidity as a
    !!< whole.
    type(state_type), intent(inout) :: state

    type(csr_sparsity) :: tmp_sparsity
    type(mesh_type) :: tmp_mesh

    tmp_mesh=extract_mesh(state, "PressureMesh")
    
    tmp_sparsity=make_sparsity(tmp_mesh, tmp_mesh, "PressureSparsity")

  end subroutine form_shallow_water_sparsity


  subroutine allocate_shallow_water_matrices(state)
    type(state_type), intent(inout) :: state

    type(csr_sparsity) :: tmp_sparsity
    type(csr_matrix) :: tmp_matrix
    type(block_csr_matrix) :: block_tmp_matrix

    tmp_sparsity=extract_csr_sparsity(state, "VelocitySparsity")
    call allocate(tmp_matrix, tmp_sparsity, name="VelocityMassMatrix")
    call insert(state, tmp_matrix, "VelocityMassMatrix")
    ! Drop the extra reference.
    call deallocate(tmp_matrix)
    
    call allocate(block_tmp_matrix, tmp_sparsity, blocks=(/2,2/), &
         name="CoupledMomentumMatrix1")
    call insert(state, block_tmp_matrix, "CoupledMomentumMatrix")
    ! Drop the extra reference.
    call deallocate(block_tmp_matrix)

    tmp_sparsity=extract_csr_sparsity(state, "PressureGradientSparsity")    
    call allocate(block_tmp_matrix, tmp_sparsity, blocks=(/1,2/), &
         name="PressureGradientMatrix")
    call insert(state, block_tmp_matrix, "PressureGradientSparsity")
    ! Drop the extra reference.
    call deallocate(tmp_matrix)
    
    tmp_sparsity=extract_csr_sparsity(state, "PressureSparsity")
    call allocate(tmp_matrix, tmp_sparsity, name="FreeSurfaceMatrix")
    call insert(state, tmp_matrix, "FreeSurfaceMatrix")
    ! Drop the extra reference.
    call deallocate(tmp_matrix)
    
  end subroutine allocate_shallow_water_matrices


  subroutine allocate_shallow_water(state)
    type(state_type), intent(inout) :: state

    type(scalar_field) :: H_rhs, H
    type(vector_field) :: U_rhs
    type(mesh_type) :: tmp_mesh
    
    tmp_mesh=extract_mesh(state, "VelocityMesh")

    call allocate(U_rhs, 2, tmp_mesh, name="VelocityRHS")
    call insert(State, U_rhs, "VelocityRHS")
    ! Drop excess reference.
    call deallocate(U_rhs)

    tmp_mesh=extract_mesh(state, "PressureMesh")

    call allocate(H_rhs, tmp_mesh, name="FreeSurfaceRHS")
    call insert(State, U_rhs, "FreeSurfaceRHS")
    ! Drop excess reference.
    call deallocate(H_rhs)

    ! Duplicate pressure as copy of FreeSurface. This is needed to keep
    ! Momentum_DG happy.
    H=extract_scalar_field(state,"FreeSurface")
    call insert(state, H, "Pressure")

  end subroutine allocate_shallow_water


  subroutine shallow_water_check_options_disabled
    character(len=FIELD_NAME_LEN) :: buffer
    integer :: itmp

    call get_option('/problem_type', buffer)
    if (trim(buffer)/="shallow water") return

    if (option_count("/material_phase")/=1) then
       FLExit("Shallow water requires exactly 1 material phase.")
    end if

    call get_option("/geometry/dimension", itmp)
    if (itmp/=2) then
       FLExit("Shallow water only works in 2D.")
    end if

    if (.not.have_option("/material_phase[0]/scalar_field::FreeSurface"))&
         & then 
       FLExit("Shallow water requires a free surface field")
    end if

    if (mesh_name("/material_phase[0]/scalar_field::FreeSurface") &
         /= "PressureMesh") then
       FLExit("Shallow water requires that free surface be on the pressure mesh")
    end if

    if (have_option("/material_phase[0]/scalar_field::Pressure"))&
         & then 
       FLExit("Do not specify a pressure field for shallow water.")
    end if
    
    if (.not.have_option("/material_phase[0]/vector_field::Velocity"))&
         & then 
       FLExit("Shallow water requires a velocity field")
    end if

  end subroutine shallow_water_check_options_disabled

end module shallow_water
