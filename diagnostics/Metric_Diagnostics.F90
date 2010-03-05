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
!    Found

#include "fdebug.h"

module metric_diagnostics

  use conformity_measurement
  use diagnostic_source_fields
  use edge_length_module
  use field_derivatives
  use field_options
  use fields
  use form_metric_field
  use fldebug
  use metric_tools
  use quicksort
  use state_module
  use spud
  use vector_tools

  implicit none
  
  private
  
  public :: calculate_scalar_edge_lengths, calculate_field_tolerance, &
    & calculate_eigenvalues_symmetric

contains

  subroutine calculate_field_tolerance(state, t_field)
    type(state_type), intent(in) :: state
    type(tensor_field), intent(inout) :: t_field
    
    type(scalar_field), pointer :: source_field
    type(vector_field), pointer :: positions
    type(tensor_field) :: mesh_metric, hessian
    logical :: allocated
    integer :: node, p, stat
    
    ewrite(1, *) "In calculate_field_tolerance"

    positions => extract_vector_field(state, "Coordinate")

    call compute_mesh_metric(positions, mesh_metric)
    
    source_field => scalar_source_field(state, t_field, allocated = allocated)

    call allocate(hessian, source_field%mesh, "Hessian")
    call compute_hessian(source_field, positions, hessian)
    call get_option(trim(complete_field_path(t_field%option_path)) // "/algorithm/p_norm", p, stat = stat)
    if(stat == SPUD_NO_ERROR) call p_norm_scale_metric(hessian, p)

    assert(hessian%mesh == mesh_metric%mesh)
    assert(hessian%mesh == t_field%mesh)
    do node = 1, node_count(mesh_metric)
      call set(t_field, node, matmul(inverse(node_val(mesh_metric, node)), absolutify_tensor(node_val(hessian, node))))
    end do

    if(allocated) deallocate(source_field)
    call deallocate(hessian)
    call deallocate(mesh_metric)
    
    ewrite(1, *) "Exiting calculate_field_tolerance"
  
  end subroutine calculate_field_tolerance

  subroutine calculate_scalar_edge_lengths(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    
    type(tensor_field) :: metric
    type(vector_field), pointer :: positions => null()
    
    positions => extract_vector_field(state, "Coordinate")
    assert(positions%mesh == s_field%mesh)
    
    call compute_mesh_metric(positions, metric)
    call get_edge_lengths(metric, s_field)
    
    call deallocate(metric)
    
  end subroutine calculate_scalar_edge_lengths
  
  subroutine calculate_eigenvalues_symmetric(state, v_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: v_field
    
    integer :: i
    integer, dimension(v_field%dim) :: permutation
    real, dimension(v_field%dim) :: evals
    real, dimension(v_field%dim, v_field%dim) :: evecs
    type(tensor_field), pointer :: source_field
    
    source_field => tensor_source_field(state, v_field)
    if(.not. source_field%mesh == v_field%mesh) then
      ewrite(-1, *) trim(v_field%name) // " mesh: " // trim(v_field%mesh%name)
      ewrite(-1, *) trim(source_field%name) // " mesh: " // trim(source_field%mesh%name)
      FLExit("Eigendecomposition mesh must match source tensor field mesh")
    end if
    
    do i = 1, node_count(source_field)
      call eigendecomposition_symmetric(node_val(source_field, i), evecs, evals)
      call qsort(evals, permutation)
      call set(v_field, i, evals(permutation))
    end do
    
  end subroutine calculate_eigenvalues_symmetric

end module metric_diagnostics
