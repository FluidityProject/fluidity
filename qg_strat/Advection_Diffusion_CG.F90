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

module advection_diffusion_CG
  !!< This module contains the Continuous Galerkin form of the advection
  !!< -diffusion equation for scalars.
  use elements
  use sparse_tools
  use fetools
  use fields
  use state_module
  use shape_functions
  use global_numbering
  use transform_elements
  use vector_tools
  use fldebug
  use vtk_interfaces
  use solvers
  use boundary_conditions
  use spud

  implicit none

  ! Buffer for output messages.
  character(len=255), private :: message

  private
  public solve_advection_diffusion_cg, construct_advection_diffusion_cg

  ! Local private control parameters. These are module-global parameters
  ! because it would be expensive and/or inconvenient to re-evaluate them
  ! on a per-element or per-face basis
  real :: dt, theta, b_l

contains

  subroutine solve_advection_diffusion_cg(field_name, state, sparsity, b)
    !!< Construct and solve the advection-diffusion equation for the given
    !!< field using continuous elements.
    
    !! Name of the field to be solved for.
    character(len=*), intent(in) :: field_name
    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state
    !! Sparsity of advection_diffusion matrix.    
    type(csr_sparsity), intent(inout) :: sparsity
    real, intent(in), optional :: b
    
    !! Tracer to be solved for.
    type(scalar_field), pointer :: T, T_old

    !! Change in T over one timestep.
    type(scalar_field) :: delta_T

    !! System matrix.
    type(csr_matrix) :: matrix

    !! Right hand side vector.
    type(scalar_field) :: rhs

    if(present(b)) then
       b_l = b
    else
       b_l = 1.
    end if

    T=>extract_scalar_field(state, field_name)
    T_old=>extract_scalar_field(state, "Old"//field_name)

    ! Reset T to value at the beginning of the timestep.
    call set(T, T_old)

    call allocate(matrix, sparsity) ! Add data space to the sparsity
    ! pattern.

    ! Ensure delta_T inherits options from T.
    call allocate(delta_T, T%mesh, "delta_T")
    delta_T%option_path = T%option_path
    call allocate(rhs, T%mesh, trim(field_name)//" RHS")
    
    call construct_advection_diffusion_cg(matrix, rhs, field_name, state)

    ! Apply boundary conditions.
    ! This is for big spring boundary conditions.
    !call apply_dirichlet_conditions(matrix, rhs, T, dt)
 
    call zero(delta_T) ! Impose zero initial guess.
    ! Solve for the change in T.
    call petsc_solve(delta_T, matrix, rhs)

    ! Add the change in T to T.
    call addto(T, delta_T, b_l*dt)

    call deallocate(delta_T)
    call deallocate(matrix)
    call deallocate(rhs)

  end subroutine solve_advection_diffusion_cg

  subroutine construct_advection_diffusion_cg(big_m, rhs, field_name,&
       & state, mass)
    !!< Construct the advection_diffusion equation for continuous elements in
    !!< acceleration form.
    !!< 
    !!< If mass is provided then the mass matrix is not added into big_m or
    !!< rhs. It is instead returned as mass. This may be useful for testing
    !!< or for solving equations otherwise than in acceleration form.
    
    !! Main advection_diffusion matrix.    
    type(csr_matrix), intent(inout) :: big_m
    !! Right hand side vector.
    type(scalar_field), intent(inout) :: rhs
    
    !! Name of the field to be advected.
    character(len=*), intent(in) :: field_name
    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state
    !! Optional separate mass matrix.
    type(csr_matrix), intent(inout), optional :: mass

    !! Position, and velocity fields.
    type(vector_field) :: X, U_nl, U_mesh
    !! Tracer to be solved for and density.
    type(scalar_field) :: T, Rho
    !! Diffusivity
    type(tensor_field) :: Diffusivity

    !! Fake source field
    type(scalar_field), pointer, save :: Source=>null(), Abs=>null()

    !! Element index
    integer :: ele

    !! Status variable for field extraction.
    integer :: stat

    !! Shape function for T
    type(element_type), pointer :: T_shape=>null()

    ewrite(1,*) "Writing advection-diffusion equation for "&
         &//trim(field_name)

    ! These names are based on the CGNS SIDS.
    T=extract_scalar_field(state, field_name)
    X=extract_vector_field(state, "Coordinate")
    U_nl=extract_vector_field(state, "NonlinearVelocity")
!    U_mesh=extract_vector_field(state, "GridVelocity")
!    Abs=>extract_scalar_field(state, "Absorption")
    Rho=extract_scalar_field(state, "Density", stat=stat)
    if (stat/=0) then
       call allocate(Rho, T%mesh, "Density", FIELD_TYPE_CONSTANT)
       call set(Rho, 1.0)
    else
       ! Grab an extra reference to cause the deallocate below to be safe.
       call incref(Rho)
    end if

    Diffusivity=extract_tensor_field(state, trim(field_name)//"Diffusivity"&
         &, stat=stat)
    if (stat/=0) then
       call allocate(Diffusivity, T%mesh, trim(field_name)//"Diffusivity",&
            FIELD_TYPE_CONSTANT)
       call zero(Diffusivity)
    else
       ! Grab an extra reference to cause the deallocate below to be safe.
       call incref(Diffusivity)
    end if

    ! Retrieve scalar options from the options dictionary.
    call get_option(trim(T%option_path)//&
         &"/prognostic/temporal_discretisation/theta", theta)
    call get_option("/timestepping/timestep", dt)

    call zero(big_m)
    call zero(RHS)
    if (present(mass)) call zero(mass)

    element_loop: do ELE=1,element_count(T)
       
       call construct_adv_diff_element_cg(ele, big_m, rhs, &
            & state, X, T, U_nl, U_mesh, Source, Abs,&
            & Diffusivity, Rho, mass) 
       
    end do element_loop

    ! Drop any extra field references.
    call deallocate(Diffusivity)
    call deallocate(Rho)

  end subroutine construct_advection_diffusion_cg

  subroutine construct_adv_diff_element_cg(ele, big_m,rhs,&
       & state, X, T, U_nl, U_mesh, Source, Abs, Diffusivity,&
       & Rho, mass)  
    !!< Construct the advection_diffusion equation for continuous elements in
    !!< acceleration form on a single element.
    implicit none
    !! Index of current element
    integer :: ele
    !! Main advection_diffusion matrix.
    type(csr_matrix), intent(inout) :: big_m
    !! Right hand side vector.
    type(scalar_field), intent(inout) :: rhs
    !! Optional separate mass matrix.
    type(csr_matrix), intent(inout), optional :: mass

    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state   
    
    !! Position and velocity.
    type(vector_field) :: X, U_nl, U_mesh

    type(scalar_field) :: T, Rho, Source, Abs
    !! Diffusivity
    type(tensor_field) :: Diffusivity

    
    ! Bilinear forms.
    real, dimension(ele_loc(T,ele), ele_loc(T,ele)) :: &
         Rho_mat
    real, dimension(ele_loc(T,ele), ele_loc(T,ele)) :: &
         Advection_mat
    real, dimension(ele_and_faces_loc(T,ele),ele_and_faces_loc(T,ele)) ::&
         & Diffusivity_mat
    real, dimension(Diffusivity%dim, Diffusivity%dim, &
         & ele_loc(Diffusivity,ele)) :: Diffusivity_ele
    
    ! Local assembly matrices.
    real, dimension(ele_loc(T,ele), ele_loc(T,ele)) :: l_T_mat

    ! Local variables.
    
    ! Count variable for loops over dimension.
    integer :: dim1, dim2
    ! Loops over faces.
    integer :: ni
    ! Array bounds for faces of the 2nd order element.
    integer :: start, finish
    
    ! Variable transform times quadrature weights.
    real, dimension(ele_ngi(T,ele)) :: detwei
    ! Transformed gradient function for tracer.
    real, dimension(ele_loc(T, ele), ele_ngi(T, ele), mesh_dim(T)) :: dt_t
    ! Transformed gradient function for velocity.
    real, dimension(ele_loc(U_nl, ele), ele_ngi(U_nl, ele), mesh_dim(T)) ::&
         & du_t 
    ! Density at quadrature points.
    real, dimension(ele_ngi(T, ele)) :: Rho_q
    ! Different velocities at quad points.
    real, dimension(U_nl%dim, ele_ngi(U_nl, ele)) :: u_nl_q, u_mesh_q
    real, dimension(ele_ngi(U_nl, ele)) :: u_div_q, u_nl_div_q

    ! Node and shape pointers.
    integer, dimension(:), pointer :: t_ele
    type(element_type), pointer :: t_shape, u_shape
    ! Neighbours of this element.
    integer, dimension(:), pointer :: neigh

    !----------------------------------------------------------------------
    ! Establish local node lists
    !----------------------------------------------------------------------

    T_ele=>ele_nodes(T,ele)  ! Velocity

    !----------------------------------------------------------------------
    ! Establish local shape functions
    !----------------------------------------------------------------------

    t_shape=>ele_shape(T, ele)
    u_shape=>ele_shape(U_nl, ele)

    ! Transform Tracer derivatives and weights into physical space.
    call transform_to_physical(X, ele, &
         & t_shape , dshape=dt_t, detwei=detwei)

    ! Transform U_nl derivatives and weights into physical space.
    call transform_to_physical(X, ele,&
         & u_shape , dshape=du_t, detwei=detwei)

    !----------------------------------------------------------------------
    ! Construct element-wise quantities.
    !----------------------------------------------------------------------

    ! Introduce grid velocities in non-linear terms. 
    ! Advecting velocity at quadrature points.
    U_nl_q=ele_val_at_quad(U_nl,ele)
    ! Divergence of advecting velocity.
    U_nl_div_q=ele_div_at_quad(U_nl, ele, du_t)

    ! Mesh velocity at quadrature points.
    !U_mesh_q=ele_val_at_quad(U_mesh,ele)
    ! Divergence of mesh movement.
    !U_div_q=ele_div_at_quad(U_mesh, ele, du_t)
    
    Rho_q=ele_val_at_quad(Rho, ele)

    Diffusivity_ele = ele_val(Diffusivity,ele)

    !----------------------------------------------------------------------
    ! Construct bilinear forms.
    !----------------------------------------------------------------------

    ! Element density matrix.
    !  /
    !  | T T Rho dV
    !  / 
    Rho_mat = shape_shape(T_shape, T_shape, detwei*Rho_q)

    ! Element advection matrix
    !  /                            /
    ! -| div(U_nl Q)T Rho dV + beta | Q ( div U_nl ) T Rho dV
    !  /                            /
    !  /                               /
    !=-| U_nl. grad Q T Rho dV+(beta-1)| Q ( div U_nl ) T Rho dV
    !  /                               /

    Advection_mat = -dshape_dot_vector_shape(dt_t, U_nl_q, t_shape, detwei&
         &*Rho_q)

    ! Mesh movement terms still to be added.

    ! Source matrix.
!    Source_mat = shape_shape(T_shape, Source_shape, detwei)

    ! Absorption matrix.
!    Abs_mat = shape_shape(T_shape, T_shape, detwei*ele_val_at_quad(Abs,ele))

    ! Diffusion.
    Diffusivity_mat=0

    Diffusivity_mat(:size(T_ele),:size(T_ele))= &
         dshape_tensor_dshape(dt_t, ele_val_at_quad(Diffusivity,ele), &
         &                    dt_t, detwei)

    !----------------------------------------------------------------------
    ! Perform global assembly.
    !----------------------------------------------------------------------

    ! Right hand side field.
    call addto(RHS, t_ele, &
         ! Advection and absorbtion
         - matmul(Advection_mat &
!         &        + Abs_mat(:,:)&
                                , ele_val(T,ele)) &
!         + matmul(Source_mat, ele_val(Source, ele))&
                                )
    
    ! Assemble matrix.
    
    ! Advection.
    l_T_mat= Advection_mat*theta*b_l*dt    !&
         ! Absorption.
 !        + Abs_mat(:,:)*theta*dt

    if (present(mass)) then
       ! Return mass separately.
       call addto(mass, t_ele, t_ele, Rho_mat)
    else
       ! Put mass in the matrix.
       l_T_mat=l_T_mat+Rho_mat
    end if

    call addto(big_m, t_ele, t_ele, l_T_mat)

    call addto(Big_m, t_ele, t_ele, &
         & Diffusivity_mat(:size(T_ele),:size(T_ele))*theta*b_l*dt)

    call addto(RHS, t_ele, &
         & -matmul(Diffusivity_mat(:size(T_ele),:size(T_ele)),ele_val(T,ele)))

  end subroutine construct_adv_diff_element_cg

end module advection_diffusion_CG
