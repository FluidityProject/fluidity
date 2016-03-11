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

module femtools

  use integer_hash_table_module
  use integer_set_module
  use data_structures
  use global_parameters
  use futils
  use element_numbering
  use grundmann_moeller_quadrature
  use polynomials
  use reference_counting
  use vector_tools
  use wandzura_quadrature
  use quadrature
  use elements
  use embed_python
  use mpi_interfaces
  use parallel_tools
  use memory_diagnostics
  use superconvergence
  use ieee_arithmetic
  use shape_functions
  use quicksort
  use halo_data_types
  use halos_base
  use halos_debug
  use halos_allocates
  use sparse_tools
  use tensors
  use picker_data_types
  use fields_data_types
  use fields_base
  use adjacency_lists
  use cv_faces
  use element_set
  use eventcounter
  use transform_elements
  use fetools
  use linked_lists
  use halos_communications
  use halos_numbering
  use halos_ownership
  use global_numbering
  use unittest_tools
  use halos_repair
  use pickers_base
  use pickers_deallocates
  use fields_allocates
  use fields_manipulation
  use metric_tools
  use parallel_fields
  use tetrahedron_intersection_module
  use unify_meshes_module
  use supermesh_construction
  use intersection_finder_module
  use fields_calculations
  use profiler
  use petsc_tools
  use sparse_tools_petsc
  use state_module
  use boundary_conditions
  use field_options
  use cvtools
  use cv_options
  use cv_shape_functions
  use vtk_interfaces
  use halos_diagnostics
  use halos_derivation
  use halos_registration
  use merge_tensors
  use node_boundary
  use vector_set
  use field_derivatives
  use cv_face_values
  use cv_upwind_values
  use c_interfaces
  use detector_data_types
  use dgtools
  use sparse_matrices_fields
  use fefields
  use meshdiagnostics
  use node_owner_finder
  use pickers_allocates
  use pickers_inquire
  use node_ownership
  use interpolation_module
  use multigrid
  use signal_vars
  use solvers
  use sparsity_patterns
  use sparsity_patterns_meshes
  use state_fields_module
  use streamfunction
  use diagnostic_fields
  use adaptive_timestepping
  use adaptive_interpolation_module
  use auxilaryoptions
  use bound_field_module
  use cgal_tools
  use cv_fields
  use detector_tools
  use detector_parallel
  use detector_move_lagrangian
  use mixing_statistics
  use surface_integrals
  use timers
  use write_state_module
  use diagnostic_variables
  use mesh_files
  use checkpoint
  use colouring
  use conservative_interpolation_module
  use coordinates
  use dg_interpolation_module
  use dynamic_bin_sort_module
  use exodusii_common
  use gmsh_common
  use lagrangian_remap
  use matrix_norms
  use pseudo_consistent_interpolation
  use read_exodusii
  use read_gmsh
  use read_triangle
  use rotated_boundary_conditions
  use smoothing_module
  use supermesh_assembly
  use vertical_extrapolation_module
  use write_gmsh
  use python_state

  implicit none

  public

end module femtools
