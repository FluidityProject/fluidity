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

module foam_flow_module
  ! This module contains the options used to solve Laplace's equation for arbitrary geometries and BC's in Fluidity and to calculate the foam velocity
  use fldebug
  use spud
  use global_parameters, only : OPTION_PATH_LEN, FIELD_NAME_LEN
  use vector_tools
  use quadrature
  use elements
  use sparse_tools
  use transform_elements
  use fetools
  use fields
  use state_module
  use field_options
  use vtk_interfaces
  use sparsity_patterns
  use sparse_matrices_fields
  use solvers
  use boundary_conditions
  use sparsity_patterns_meshes
  use petsc_solve_state_module
  use field_derivatives
  use gradient_matrix_cg, only: assemble_gradient_matrix_cg
  use state_matrices_module

  implicit none

  private
  public :: calculate_potential_flow, calculate_foam_velocity

  character(len = *), parameter, public :: phi_name = "FoamVelocityPotential"
  character(len = *), parameter, public :: phi_sparsity_name = "FoamVelocityPotentialSparsity"
  character(len = *), parameter, public :: phi_rhs_name = "FoamVelocityPotentialRhs"
  character(len = *), parameter, public :: phi_m_name = "FoamVelocityPotentialMatrix"
  character(len = *), parameter, public :: foamvel_name = "FoamVelocity"

contains

  subroutine calculate_potential_flow(state, phi)

    type(state_type), intent(inout) :: state
    type(scalar_field), optional, intent(out) :: phi
    type(csr_matrix) :: phi_m
    type(scalar_field) :: phi_rhs


    ! Step 1: Allocate and insert objects into state / extract objects from
    ! state and take references
    call initialise_potential_flow(phi, phi_rhs, state, phi_m)


    ! Step 2: Assemble
    call assemble_potential_flow_cg(phi_rhs, state, phi_m, phi)

    ! Step 3: Solve
    call solve_potential_flow(phi_m, phi_rhs, phi, state)

    ! Step 4: Drop references
    call deallocate(phi_m)
    call deallocate(phi_rhs)


  end subroutine calculate_potential_flow


  subroutine initialise_potential_flow(phi, phi_rhs, state, phi_m)
    ! Allocate / extract FoamVelocityPotential variables. phi, phi_m and
    ! phi_rhs all take references in this routine and, if new objects are
    ! constructed, are inserted into state.

    type(scalar_field), target, intent(out) :: phi
    type(csr_matrix), intent(out) :: phi_m
    type(scalar_field), intent(out) :: phi_rhs
    type(state_type), intent(inout) :: state

    type(csr_sparsity), pointer :: phi_sparsity => null()

    phi = extract_scalar_field(state, phi_name)

    ! Matrix sparsity
    phi_sparsity => get_csr_sparsity_firstorder(state, phi%mesh, phi%mesh)
    call allocate(phi_m, phi_sparsity, name = phi_m_name)

    ! RHS
    call allocate(phi_rhs, phi%mesh, phi_rhs_name)

  end subroutine initialise_potential_flow


  subroutine assemble_potential_flow_cg(phi_rhs, state, phi_m, phi)
    ! Assemble Laplace's equation for FoamVelocityPotential
    type(scalar_field), target, intent(inout) :: phi_rhs
    type(state_type), intent(inout) :: state
    type(csr_matrix), optional, intent(inout) :: phi_m
    type(scalar_field), intent(inout) :: phi
    type(scalar_field) :: foamvelocitypotential_bc
    integer, dimension(:), allocatable :: foamvelocitypotential_bc_type
    integer :: ele, sele
    type(vector_field), pointer :: positions => null()

    ewrite(1, *) "In assemble_potential_flow_cg"

    positions => extract_vector_field(state, "Coordinate")

    call zero(phi_m)

    call zero(phi_rhs)


    ! element loop
    element_loop: do ele=1, element_count(phi)
      call assemble_potential_flow_element_cg(ele, phi_m, phi_rhs, positions, phi)
    end do element_loop
    ! end element loop

      allocate(foamvelocitypotential_bc_type(surface_element_count(phi)))
      call get_entire_boundary_condition(phi, (/ "neumann       "/), foamvelocitypotential_bc, foamvelocitypotential_bc_type)


    ! surface element loop
    surface_element_loop: do sele=1, surface_element_count(phi)
      ele = face_ele(positions, sele)
      call assemble_potential_flow_surface_element_cg(sele, ele, positions, &
        & phi_rhs, phi, foamvelocitypotential_bc)
    end do surface_element_loop
    ! end surface element loop

      call deallocate(foamvelocitypotential_bc)

    ! Set a reference node, zero at the first node of the first process
    ! (should be called by all processes though)
    ! This needs to be done every time, to zero the rhs
    call set_reference_node(phi_m, 1, phi_rhs, 0.0)

    ewrite_minmax(phi_rhs)

    ewrite(1, *) "Exiting assemble_potential_flow_cg"

  end subroutine assemble_potential_flow_cg


  subroutine assemble_potential_flow_element_cg(ele, phi_m, phi_rhs, positions, phi)
    !!< Assemble the element matrix contributions
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: phi
    type(csr_matrix), optional, intent(inout) :: phi_m
    type(scalar_field), intent(inout) :: phi_rhs

    ! Node numbers of phi element:
    integer, dimension(:), pointer :: phi_ele
    ! Locations of nodes:
    real, dimension(positions%dim, ele_loc(positions, ele)) :: x_val
    ! Shape functions:
    type(element_type), pointer :: phi_shape, x_shape
    ! Coordinate transform*quadrature weights :
    real, dimension(ele_ngi(positions, ele)) :: detwei

    ! Locations of quadrature points
    real, dimension(positions%dim,ele_ngi(positions,ele)) :: x_quad

    !Derivatives of shape function:
    real, dimension(ele_loc(phi, ele), ele_ngi(phi, ele), positions%dim) :: dphi_t

    ! Local Laplacian matrix
    real, dimension(ele_loc(phi, ele), ele_loc(phi, ele)) :: phi_mat

    phi_ele => ele_nodes(phi, ele)
    phi_shape=>ele_shape(phi, ele)

    x_shape=> ele_shape(positions, ele)
    ! Location of local vertices
    x_val = ele_val(positions,ele)

    ! Locations of quadrature points
    x_quad=ele_val_at_quad(positions, ele)


    ! Transform derivatives and weights into physical space
    call transform_to_physical(positions, ele, phi_shape, &
      & dshape = dphi_t, detwei = detwei)

    ! Local assembly:
    ! /
    ! | grad N_A dot grad N_B dV
    ! /
    phi_mat=dshape_dot_dshape(dphi_t, dphi_t, detwei)

    ! Global assembly:
    call addto(phi_m, phi_ele, phi_ele, phi_mat)

    ! The rhs is always zero for Laplace's equation
    call zero(phi_rhs)

  end subroutine assemble_potential_flow_element_cg


  subroutine assemble_potential_flow_surface_element_cg(sele, ele, positions, &
        & phi_rhs, phi, foamvelocitypotential_bc)
    integer, intent(in) :: sele, ele
    type(scalar_field), intent(inout) :: phi_rhs

    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: phi

    type(scalar_field), intent(in) :: foamvelocitypotential_bc

    integer, dimension(face_loc(phi, sele)) :: phi_nodes_bdy
    real, dimension(positions%dim, ele_loc(positions, ele)) :: x_ele
    real, dimension(positions%dim, face_loc(positions, sele)) :: x_ele_bdy
    type(element_type), pointer :: x_shape, x_shape_bdy, phi_shape

    real, dimension(face_ngi(positions, sele)) :: detwei_bdy
    real, dimension(face_loc(phi, sele), face_loc(phi, sele)) :: face_mat

    x_ele = ele_val(positions, ele)
    x_ele_bdy = face_val(positions, sele)

    x_shape=> ele_shape(positions, ele)
    x_shape_bdy=>face_shape(positions, sele)


    ! Transform derivatives and weights into physical space (calculate quadrature weights)
    call transform_facet_to_physical(positions, sele, detwei_f=detwei_bdy)


    ! integral over the face of the form \int N_i N_j
    ! where N_i and N_j are shape functions of phi
    phi_shape=> face_shape(phi, sele)
    face_mat=shape_shape(phi_shape, phi_shape, detwei_bdy)

    ! global node numbers of nodes of this face
    phi_nodes_bdy = face_global_nodes(phi, sele)

    ! global node numbers of nodes of this face in phi%mesh
    ! this implements a Neumann bc with values given by foamvelocitypotential_bc
    call addto(phi_rhs, phi_nodes_bdy, matmul(face_mat, ele_val(foamvelocitypotential_bc, sele) ))

  end subroutine assemble_potential_flow_surface_element_cg



  subroutine solve_potential_flow(phi_m, phi_rhs, phi, state)

    type(csr_matrix), intent(in) :: phi_m
    type(scalar_field), intent(in) :: phi_rhs
    type(scalar_field), intent(inout) :: phi
    type(state_type), intent(inout) :: state

    call petsc_solve(phi, phi_m, phi_rhs, state)

    ewrite_minmax(phi)

  end subroutine solve_potential_flow



  subroutine calculate_foam_velocity(state, foamvel)

    type(state_type), intent(inout) :: state
    type(vector_field), pointer, intent(out) :: foamvel
    type(vector_field), pointer :: liquidvelocity
    type(scalar_field), pointer :: field
    type(csr_sparsity) :: gradient_sparsity
    type(block_csr_matrix) :: C_m
    type(csr_sparsity) :: mass_sparsity
    type(csr_matrix) :: mass
    type(vector_field) :: cfield

    integer :: i, stat
    
    foamvel => extract_vector_field(state, foamvel_name)
    if (have_option(trim(foamvel%option_path)//'/diagnostic')) then
      field=>extract_scalar_field(state, "FoamVelocityPotential")

      !The following does the same as the subroutine calculate_grad_fe in Diagnostic_Fields_Matrices.F90
      call allocate(cfield, foamvel%dim, foamvel%mesh, name="CField")

      ! Sparsity of C^T - the transpose of the gradient operator.
      gradient_sparsity=make_sparsity(foamvel%mesh, field%mesh, "GradientSparsity")
      call allocate(C_m, gradient_sparsity, (/foamvel%dim, 1/), name="GradientMatrix" )

      mass_sparsity=make_sparsity(foamvel%mesh, foamvel%mesh, "MassSparsity")
      call allocate(mass, mass_sparsity, name="MassMatrix")

      call assemble_gradient_matrix_cg(C_m, state, &
                                         test_mesh=foamvel%mesh, field=field, &
                                         option_path=trim(foamvel%option_path), &
                                         grad_mass=mass)

      call mult(cfield, C_m, field)

      call zero(foamvel)
      call petsc_solve(foamvel, mass, cfield)

      call deallocate(gradient_sparsity)
      call deallocate(C_m)
      call deallocate(mass_sparsity)
      call deallocate(mass)
      call deallocate(cfield)

      !Changing the sign since foamvel= -grad phi
      do i=1, node_count(foamvel)
        call set(foamvel, i, ( -node_val(foamvel, i) ) )
      enddo

      ewrite_minmax(foamvel)

    endif



  end subroutine calculate_foam_velocity



end module foam_flow_module


