!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineeringp
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

subroutine project_to_continuous_parallel(vtuname, vtuname_len, trianglename,&
     & trianglename_len) 
  !!< Given a vtu file containing fields on a discontinuous mesh, and the
  !!< triangle files for the corresponding continuous mesh, produce a vtu
  !!< with its fields projected onto the continuous mesh.
  
  use fefields
  use fields
  use flcomms_io
  use halos
  use global_parameters, only : current_debug_level, halo_tag, halo_tag_p
  use parallel_tools
  use read_triangle
  use solvers
  use sparsity_patterns
  use sparse_tools
  use state_module
  use vtk_interfaces
  
  implicit none
  
  integer, intent(in):: vtuname_len, trianglename_len
  character(len=vtuname_len), intent(in):: vtuname
  character(len=trianglename_len), intent(in):: trianglename

  character(len = *), parameter :: solver_option_path = "/dummy"
  integer :: i, nfields
  ! VTK can store at most P2 fields. Hence degree 4 allows complete quadrature
  ! in the Galerkin projection.
  integer, parameter :: quad_degree = 4
  type(csr_matrix) :: mass
  type(csr_sparsity) :: sparsity
  type(mesh_type), pointer :: cg_mesh
  type(scalar_field) :: cg_s_field
  type(scalar_field), dimension(:), allocatable :: rhs
  type(scalar_field), dimension(:), pointer :: cg_state_fields, dg_state_fields
  type(scalar_field), pointer :: dg_s_field
  type(state_type) :: cg_state, dg_state
  type(tensor_field) :: cg_t_field
  type(tensor_field), pointer :: dg_t_field
  type(vector_field) :: cg_v_field
  type(vector_field), pointer :: dg_v_field
  type(vector_field), target :: cg_coordinate
  
  !current_debug_level = 2
  
  if(.not. isparallel()) then
    FLExit("project_to_continuous_parallel should only be used in parallel")
  end if

  call vtk_read_state(parallel_filename(vtuname, ".vtu"), dg_state, quad_degree = quad_degree)
  
  cg_coordinate = read_triangle_files(trianglename, quad_degree = quad_degree)
  call read_halos(trianglename)
  cg_mesh => cg_coordinate%mesh
  allocate(cg_mesh%halos(2))
  call import_halo(halo_tag, cg_mesh%halos(1))
  call import_halo(halo_tag_p, cg_mesh%halos(2))
  allocate(cg_mesh%element_halos(2))
  call derive_element_halo_from_node_halo(cg_mesh, &
    & ordering_scheme = HALO_ORDER_TRAILING_RECEIVES)

  do i = 1, scalar_field_count(dg_state)
    dg_s_field => extract_scalar_field(dg_state, i)
    call allocate(cg_s_field, cg_mesh, name = dg_s_field%name)
    call insert(cg_state, cg_s_field, cg_s_field%name)
    call deallocate(cg_s_field)
  end do
  
  do i = 1, vector_field_count(dg_state)
    dg_v_field => extract_vector_field(dg_state, i)
    call allocate(cg_v_field, dg_v_field%dim, cg_mesh, name = dg_v_field%name)
    call insert(cg_state, cg_v_field, cg_v_field%name)
    call deallocate(cg_v_field)
  end do
  
  do i = 1, tensor_field_count(dg_state)
    dg_t_field => extract_tensor_field(dg_state, i)
    call allocate(cg_t_field, cg_mesh, name = dg_t_field%name)
    call insert(cg_state, cg_t_field, cg_t_field%name)
    call deallocate(cg_t_field)
  end do

  sparsity = make_sparsity(cg_mesh, cg_mesh, name = "MassSparsity")
  call allocate(mass, sparsity, name = "MassMatrix")
  call deallocate(sparsity)
  call compute_mass(cg_coordinate, cg_mesh, mass)

  call collapse_state(dg_state, dg_state_fields)
  nfields = size(dg_state_fields)
  
  allocate(rhs(nfields))
  do i = 1, nfields
    call allocate(rhs(i), cg_mesh, name = trim(dg_state_fields(i)%name) // "Rhs")
    call zero(rhs(i))
  end do

  do i = 1, ele_count(cg_mesh)
    call assemble_gp_rhs_ele(i, cg_mesh, cg_coordinate, dg_state_fields, rhs)
  end do
  deallocate(dg_state_fields)
  call deallocate(dg_state)
  
  call collapse_state(cg_state, cg_state_fields)
  call set_solver_options(solver_option_path, &
    & ksptype = "cg", pctype = "sor", atol = epsilon(0.0), rtol = 0.0, max_its = 2000, start_from_zero = .true.)
  call petsc_solve(cg_state_fields, mass, rhs, option_path = solver_option_path)
  
  deallocate(cg_state_fields)
  call deallocate(mass)
  do i = 1, nfields
    call deallocate(rhs(i))
  end do
  deallocate(rhs)
  
  call insert(cg_state, cg_coordinate, "Coordinate")
  call insert(cg_state, cg_mesh, "CoordinateMesh")
  call deallocate(cg_coordinate)
  
  call vtk_write_state(trim(vtuname) //"_continuous.pvtu", state = (/cg_state/))
  call deallocate(cg_state)
  
  call print_references(0)
       
contains

  subroutine assemble_gp_rhs_ele(ele, mesh, positions, fields, rhs)
    integer, intent(in) :: ele
    type(mesh_type), intent(in) :: mesh
    type(vector_field), intent(in) :: positions
    type(scalar_field), dimension(:), intent(in) :: fields
    type(scalar_field), dimension(size(fields)), intent(inout) :: rhs
       
    integer :: i
    integer, dimension(:), pointer :: element_nodes
    real, dimension(ele_ngi(mesh, ele)) :: detwei
    type(element_type), pointer :: shape
    
    call transform_to_physical(positions, ele, detwei = detwei)
    
    shape => ele_shape(mesh, ele)
    element_nodes => ele_nodes(mesh, ele)

    do i = 1, size(fields)
      call addto(rhs(i), element_nodes, shape_rhs(shape, detwei * ele_val_at_quad(fields(i), ele)))
    end do
    
  end subroutine assemble_gp_rhs_ele
       
end subroutine project_to_continuous_parallel
