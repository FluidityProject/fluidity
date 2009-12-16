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
  use edge_length_module
  use fields
  use fldebug
  use state_module
  use field_derivatives
  use global_parameters, only: FIELD_NAME_LEN
  use spud
  use vector_tools
  use metric_tools

  implicit none
  
  private
  
  public :: calculate_edge_lengths, calculate_field_tolerance
  
  interface calculate_edge_lengths
    module procedure calculate_edge_lengths_scalar
  end interface calculate_edge_lengths

contains

  subroutine calculate_field_tolerance(state, field_tolerance)
    type(state_type), intent(in) :: state
    type(tensor_field), intent(inout) :: field_tolerance
    type(tensor_field) :: mesh_metric, hessian
    type(vector_field), pointer :: positions
    type(scalar_field), pointer :: field

    logical :: allocated
    character(len=FIELD_NAME_LEN) :: field_name

    integer :: node
    
    ewrite(1, *) "In calculate_field_tolerance"

    positions => extract_vector_field(state, "Coordinate")

    call compute_mesh_metric(positions, mesh_metric)
    call get_option(trim(field_tolerance%option_path) // "/diagnostic/source_field_name", field_name)
    ewrite(2, *) "Source field: " // trim(field_name)
    field => extract_scalar_field(state, trim(field_name), allocated=allocated)

    call allocate(hessian, field%mesh, trim(field_name) // "Hessian")
    call compute_hessian(field, positions, hessian)

    assert(hessian%mesh == mesh_metric%mesh)
    assert(hessian%mesh == field_tolerance%mesh)
    do node=1,node_count(mesh_metric)
      call set(field_tolerance, node, matmul(inverse(node_val(mesh_metric, node)), absolutify_tensor(node_val(hessian, node))))
    end do

    if (allocated) then
      deallocate(field)
    end if

    call deallocate(hessian)
    call deallocate(mesh_metric)
    
    ewrite(1, *) "Exiting calculate_field_tolerance"
  
  end subroutine calculate_field_tolerance

  subroutine calculate_edge_lengths_scalar(state, edge_lengths)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: edge_lengths
    
    type(tensor_field) :: metric
    type(vector_field), pointer :: positions => null()
    
    positions => extract_vector_field(state, "Coordinate")
    assert(positions%mesh == edge_lengths%mesh)
    
    call compute_mesh_metric(positions, metric)
    call get_edge_lengths(metric, edge_lengths)
    
    call deallocate(metric)
    
  end subroutine calculate_edge_lengths_scalar

end module metric_diagnostics
