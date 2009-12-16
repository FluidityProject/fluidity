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
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
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

subroutine pseudo_mesh(filename, filename_len, target_elements)

  use adapt_state_module
  use conformity_measurement
  use fields
  use fldebug
  use global_parameters, only : current_debug_level
  use interpolation_module
  use limit_metric_module
  use read_triangle
  use parallel_tools
  use reference_counting
  use write_triangle

  implicit none
  
  integer, intent(in) :: filename_len
  
  character(len = filename_len), intent(in) :: filename
  integer, intent(in) :: target_elements
  
  integer :: i
  integer, parameter :: adapt_iterations = 5
  real, dimension(:, :), allocatable :: unit_matrix
  type(tensor_field) :: metric, target_metric
  type(vector_field) :: input_positions, output_positions, stage_2_input_positions
  
  ewrite(1, *) "In pseudo_mesh"
  
  if(isparallel()) then
    FLAbort("pseudo_mesh does not work in parallel")
  end if
  
  ewrite(2, *) "Target input elements? ", target_elements /= 0
  
  ! Step 1: Input
  input_positions = read_triangle_files(filename, quad_degree = 4)
  
  ! Step 2: Adapt. We do this in two stages - the first to pull us a long way
  ! from the target mesh, and the second to pull us back towards it. This
  ! ensures that adaptivity actually does something, rather than just giving us
  ! a copy of the input mesh.
  
  ! Adapt stage 1: Coarsen as much as possible everywhere
  ewrite(2, *) "Adapt stage 1: Coarsening"
  
  call allocate(metric, input_positions%mesh, "Metric", field_type = FIELD_TYPE_CONSTANT)
  allocate(unit_matrix(metric%dim, metric%dim))
  unit_matrix = 0.0
  do i = 1, metric%dim
    unit_matrix(i, i) = 1.0
  end do
  call set(metric, unit_matrix)
  deallocate(unit_matrix)
  call limit_metric(input_positions, metric, target_nodes = 1)
  
  call adapt_mesh(input_positions, metric, output_positions)
  call deallocate(metric)
  
  call write_triangle_files("pseudo_mesh_stage_1", output_positions)
  ewrite(2, *) "Node count = ", node_count(output_positions)
  ewrite(2, *) "Element count = ", ele_count(output_positions)
  
  ! Adapt stage 2: Now adapt to the actual target metric
  ewrite(2, *) "Adapt stage 2: Targetting input mesh metric"
  
  call compute_mesh_metric(input_positions, target_metric)
  
  ! Subcycle, as each time we're targetting an interpolated metric
  do i = 1, adapt_iterations
    ewrite(2, "(a,i0,a,i0)") "Iteration ", i, " of ", adapt_iterations
  
    call allocate(metric, output_positions%mesh, "InterpolatedTargetMetric")
    call linear_interpolation(target_metric, input_positions, metric, output_positions)
    if(target_elements /= 0) then
      ! Target elements, to ensure convergence of repeated pseudo-meshing
      call limit_metric_elements(output_positions, metric, target_eles = ele_count(input_positions))
    else
      ! If we're not targetting, limit the metric anyway (but loosly) as it's
      ! possible that there are very small or large edge lengths at the
      ! boundaries of the input mesh
      call limit_metric(output_positions, metric, min_nodes = max(int(0.5 * node_count(input_positions)), 1), max_nodes = 2 * node_count(input_positions))
    end if
    
    ! Copy the output from the last adapt, for use as input for the next adapt
    call allocate(stage_2_input_positions, output_positions%dim, output_positions%mesh, output_positions%name)
    call set(stage_2_input_positions, output_positions)
    call deallocate(output_positions)
    
    call adapt_mesh(stage_2_input_positions, metric, output_positions)
    call deallocate(metric)
    call deallocate(stage_2_input_positions)
    
    call write_triangle_files("pseudo_mesh_stage_2_" // int2str(i), output_positions)
    ewrite(2, *) "Node count = ", node_count(output_positions)
    ewrite(2, *) "Element count = ", ele_count(output_positions)
  end do
  
  call deallocate(target_metric)
  call deallocate(input_positions)
  
  ! Step 3: Output
  
  call write_triangle_files("pseudo_mesh", output_positions)
  call deallocate(output_positions)
  
  call print_references(0)
  
  ewrite(1, *) "Exiting pseudo_mesh"
  
end subroutine pseudo_mesh
