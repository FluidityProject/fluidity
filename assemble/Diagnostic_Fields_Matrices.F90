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

module diagnostic_fields_matrices
  !!< A module to link to diagnostic variable calculations.

  use global_parameters, only:FIELD_NAME_LEN,ACCTIM, new_options
  use fields
  use sparse_matrices_fields
  use field_derivatives
  use state_module
  use futils
  use fetools
  use spud
  use parallel_fields, only: zero_non_owned
  use divergence_matrix_cv, only: assemble_divergence_matrix_cv
  use divergence_matrix_cg, only: assemble_divergence_matrix_cg
  use gradient_matrix_cg, only: assemble_gradient_matrix_cg
  use parallel_tools
  use sparsity_patterns, only: make_sparsity
  use state_fields_module
  use solvers
  use transform_elements
  use sparse_tools
  use sparsity_patterns_meshes
  use state_matrices_module
  use halos

  implicit none

  private
  
  public :: calculate_divergence_cv, calculate_divergence_fe, &
            calculate_div_t_cv, calculate_div_t_fe, &
            calculate_grad_fe

contains

  subroutine calculate_divergence_cv(state, div)

      type(state_type), intent(inout) :: state
      type(scalar_field), intent(inout) :: div

      type(vector_field), pointer :: field, x

      type(block_csr_matrix), pointer :: CT_m

      character(len=FIELD_NAME_LEN) :: field_name

      type(scalar_field), pointer :: lumpedmass
      type(scalar_field) :: inverse_lumpedmass
      type(scalar_field) :: ctfield, ct_rhs

      call get_option(trim(div%option_path)//"/diagnostic/field_name", field_name)

      field=>extract_vector_field(state, trim(field_name))
      x=>extract_vector_field(state, "Coordinate")

      call allocate(ctfield, div%mesh, name="CTField")
      call allocate(ct_rhs, div%mesh, name="CTRHS")

      if(element_degree(div, 1)>1) then
        lumpedmass => get_lumped_mass_on_submesh(state, div%mesh)
      else
        lumpedmass => get_lumped_mass(state, div%mesh)
      end if

      CT_m => get_divergence_matrix_cv(state, test_mesh=div%mesh, field=field, div_rhs=ct_rhs)

      call mult(ctfield, CT_m, field)
      call addto(ctfield, ct_rhs, -1.0)

      call allocate(inverse_lumpedmass, lumpedmass%mesh, "InverseLumpedMass")
      call invert(lumpedmass, inverse_lumpedmass)
      call set(div, ctfield)
      call scale(div, inverse_lumpedmass)
      
      call deallocate(inverse_lumpedmass)
      call deallocate(ctfield)
      call deallocate(ct_rhs)

      call halo_update(div)

  end subroutine calculate_divergence_cv

  subroutine calculate_div_t_cv(state, grad)

      type(state_type), intent(inout) :: state
      type(vector_field), intent(inout) :: grad

      type(scalar_field), pointer :: field

      type(block_csr_matrix), pointer :: CT_m

      character(len=FIELD_NAME_LEN) :: field_name

      type(scalar_field), pointer :: lumped_mass
      type(scalar_field) :: inverse_lumped_mass
      type(csr_matrix), pointer :: mass, inverse_mass
      type(vector_field) :: cfield
      type(scalar_field) :: mag, inverse_mag

      logical :: lump_mass, dg, normalise

      call get_option(trim(grad%option_path)//"/diagnostic/field_name", field_name)
      
      dg = (continuity(grad)<0)
      lump_mass = have_option(trim(grad%option_path)//&
                  &"/diagnostic/lump_mass_matrix")
      normalise = have_option(trim(grad%option_path)//&
                  &"/diagnostic/normalise")
      
      field=>extract_scalar_field(state, trim(field_name))

      call allocate(cfield, grad%dim, grad%mesh, name="CField")

      CT_m => get_divergence_matrix_cv(state, test_mesh=field%mesh, field=grad, exclude_boundaries=.true.)

      call mult_T(cfield, CT_m, field)
      call scale(cfield, -1.0)
      
      call zero(grad)
      
      if(lump_mass) then
        
        lumped_mass => get_lumped_mass(state, grad%mesh)
        call allocate(inverse_lumped_mass, lumped_mass%mesh, "InverseLumpedMass")
        
        call invert(lumped_mass, inverse_lumped_mass)
        call set(grad, cfield)
        call scale(grad, inverse_lumped_mass)
       
        call deallocate(inverse_lumped_mass)
      else if(dg) then
        inverse_mass => get_dg_inverse_mass(state, grad%mesh)
        call mult(grad, inverse_mass, cfield)
      else
        mass => get_mass_matrix(state, grad%mesh)
        call petsc_solve(grad, mass, cfield)
      end if

      if(normalise) then
        mag = magnitude(grad)
        call allocate(inverse_mag, mag%mesh, "InverseMagnitude")
        
        call invert(mag, inverse_mag, tolerance=epsilon(0.0))
        call scale(grad, inverse_mag)
        
        call deallocate(inverse_mag)
        call deallocate(mag)
      end if
      
      call deallocate(cfield)

  end subroutine calculate_div_t_cv

  subroutine calculate_divergence_fe(state, div)

      type(state_type), intent(inout) :: state
      type(scalar_field), intent(inout) :: div

      type(vector_field), pointer :: field, x

      type(csr_sparsity) :: divergence_sparsity
      type(block_csr_matrix) :: CT_m

      character(len=FIELD_NAME_LEN) :: field_name

      type(csr_sparsity) :: mass_sparsity
      type(csr_matrix) :: mass
      type(scalar_field) :: ctfield, ct_rhs

      call get_option(trim(div%option_path)//"/diagnostic/field_name", field_name)

      field=>extract_vector_field(state, trim(field_name))
      x=>extract_vector_field(state, "Coordinate")

      call allocate(ctfield, div%mesh, name="CTField")

      divergence_sparsity=make_sparsity(div%mesh, field%mesh, "DivergenceSparsity")
      call allocate(CT_m, divergence_sparsity, (/1, field%dim/), name="DivergenceMatrix" )
      call allocate(ct_rhs, div%mesh, name="CTRHS")


      mass_sparsity=make_sparsity(div%mesh, div%mesh, "MassSparsity")
      call allocate(mass, mass_sparsity, name="MassMatrix")
      call zero(mass)

      call assemble_divergence_matrix_cg(CT_m, state, ct_rhs=ct_rhs, &
                                         test_mesh=div%mesh, field=field, &
                                         option_path=div%option_path, div_mass=mass)

      call mult(ctfield, CT_m, field)
      call addto(ctfield, ct_rhs, -1.0)

      call zero(div)
      call petsc_solve(div, mass, ctfield)

      call deallocate(CT_m)
      call deallocate(ct_rhs)
      call deallocate(ctfield)
      call deallocate(divergence_sparsity)
      call deallocate(mass_sparsity)
      call deallocate(mass)

  end subroutine calculate_divergence_fe

  subroutine calculate_div_t_fe(state, grad)

      type(state_type), intent(inout) :: state
      type(vector_field), intent(inout) :: grad

      type(scalar_field), pointer :: field

      type(csr_sparsity) :: divergence_sparsity
      type(block_csr_matrix) :: CT_m

      character(len=FIELD_NAME_LEN) :: field_name

      type(csr_sparsity) :: mass_sparsity
      type(csr_matrix) :: mass
      type(vector_field) :: cfield


      call get_option(trim(grad%option_path)//"/diagnostic/field_name", field_name)

      field=>extract_scalar_field(state, trim(field_name))

      call allocate(cfield, grad%dim, grad%mesh, name="CField")

      ! Sparsity of C^T - the transpose of the gradient operator.
      divergence_sparsity=make_sparsity(field%mesh, grad%mesh, "DivergenceSparsity")
      call allocate(CT_m, divergence_sparsity, (/1, grad%dim/), name="DivergenceMatrix" )

      mass_sparsity=make_sparsity(grad%mesh, grad%mesh, "MassSparsity")
      call allocate(mass, mass_sparsity, name="MassMatrix")

      call assemble_divergence_matrix_cg(CT_m, state, &
                                         test_mesh=field%mesh, field=grad, &
                                         grad_mass=mass)

      call mult_T(cfield, CT_m, field)
      call scale(cfield, -1.0)

      call zero(grad)
      call petsc_solve(grad, mass, cfield)

      call deallocate(divergence_sparsity)
      call deallocate(CT_m)
      call deallocate(mass_sparsity)
      call deallocate(mass)
      call deallocate(cfield)

  end subroutine calculate_div_t_fe

  subroutine calculate_grad_fe(state, grad)

      type(state_type), intent(inout) :: state
      type(vector_field), intent(inout) :: grad

      type(scalar_field), pointer :: field

      type(csr_sparsity) :: gradient_sparsity
      type(block_csr_matrix) :: C_m

      character(len=FIELD_NAME_LEN) :: field_name

      type(csr_sparsity) :: mass_sparsity
      type(csr_matrix) :: mass
      type(vector_field) :: cfield


      call get_option(trim(grad%option_path)//"/diagnostic/field_name", field_name)

      field=>extract_scalar_field(state, trim(field_name))

      call allocate(cfield, grad%dim, grad%mesh, name="CField")

      ! Sparsity of C^T - the transpose of the gradient operator.
      gradient_sparsity=make_sparsity(grad%mesh, field%mesh, "GradientSparsity")
      call allocate(C_m, gradient_sparsity, (/grad%dim, 1/), name="GradientMatrix" )

      mass_sparsity=make_sparsity(grad%mesh, grad%mesh, "MassSparsity")
      call allocate(mass, mass_sparsity, name="MassMatrix")

      call assemble_gradient_matrix_cg(C_m, state, &
                                         test_mesh=grad%mesh, field=field, &
                                         option_path=trim(grad%option_path), &
                                         grad_mass=mass)

      call mult(cfield, C_m, field)

      call zero(grad)
      call petsc_solve(grad, mass, cfield)

      call deallocate(gradient_sparsity)
      call deallocate(C_m)
      call deallocate(mass_sparsity)
      call deallocate(mass)
      call deallocate(cfield)

  end subroutine calculate_grad_fe

end module diagnostic_fields_matrices
