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

subroutine streamfunction_2d(input_basename, input_basename_len, &
  & output_basename, output_basename_len)
  
  use fields
  use fldebug
  use global_parameters, only : FIELD_NAME_LEN
  use reference_counting
  use solvers
  use sparse_tools
  use sparsity_patterns_meshes
  use state_module
  use vtk_interfaces
  
  implicit none
  
  integer :: input_basename_len
  integer :: output_basename_len
  
  character(len = input_basename_len), intent(in) :: input_basename
  character(len = output_basename_len), intent(in) :: output_basename
  
  character(len = FIELD_NAME_LEN) :: model
  integer :: i, stat
  type(csr_matrix) :: matrix
  type(csr_sparsity), pointer :: sparsity
  type(scalar_field) :: psi, rhs
  type(state_type), pointer :: state
  type(state_type), dimension(1), target :: states
  type(vector_field), pointer :: positions, velocity
  
  ewrite(1, *) "In streamfunction_2d"
  
  state => states(1)
  
  call vtk_read_state(input_basename, state)
  
  positions => extract_vector_field(state, "Coordinate")
  if(positions%dim /= 2) then
    ewrite(-1,*) "Your problem is of dimension ", positions%dim
    FLExit("streamfunction_2d requires a 2D input vtu")
  end if
  if(positions%mesh%continuity /= 0) then
    ewrite(-1,*) "Your Coordinates mesh is not continuous"
    FLExit("streamfunction_2d requires a continuous input vtu")
  end if
  
  velocity => extract_vector_field(state, "Velocity")
  assert(velocity%dim == positions%dim)
  assert(ele_count(velocity) == ele_count(positions))
  do i = 1, velocity%dim
    ewrite_minmax(velocity%val(i,:))
  end do
  
  psi = extract_scalar_field(state, "StreamFunction", stat)
  ! Horrible hack - actually need to allocate a new mesh here and add the
  ! faces (but this is ok as we cheat below) - mesh%faces really needs to be a
  ! pointer to a pointer to make this work nicely
  if(stat == 0) then
    if(.not. associated(psi%mesh%faces)) call add_faces(psi%mesh)
    call incref(psi)
  else
    call allocate(psi, positions%mesh, "StreamFunction")
    call zero(psi)
    if(.not. associated(psi%mesh%faces)) call add_faces(psi%mesh)
    call insert(state, psi, psi%name)
  end if
  assert(mesh_dim(psi) == positions%dim)
  assert(ele_count(psi) == ele_count(positions))
  ewrite_minmax(psi%val)
  
  sparsity => get_csr_sparsity_firstorder(state, psi%mesh, psi%mesh)
  call allocate(matrix, sparsity, name = trim(psi%name) // "Matrix")
  call allocate(rhs, psi%mesh, name = trim(psi%name) // "Rhs")
  
  call zero(matrix)
  call zero(rhs)
  do i = 1, ele_count(psi)
    call assemble_streamfunction_2d_element(i, matrix, rhs, positions, velocity)
  end do
  ewrite_minmax(rhs%val)
  
  do i = 1, surface_element_count(psi)
    call set_inactive(matrix, face_global_nodes(rhs, i))
    call set(rhs, face_global_nodes(rhs, i), spread(0.0, 1, face_loc(rhs, i)))
  end do
  
  call set_solver_options(psi, ksptype = "cg", pctype = "sor", rtol = 1.0e-10, max_its = 3000)
  call petsc_solve(psi, matrix, rhs)
  ewrite_minmax(psi%val)
  
  if(has_mesh(state, "CoordinateMesh")) then
    model = "CoordinateMesh"
  else
    model = "Mesh"
  end if
  call vtk_write_state(output_basename, model = model, state = states)

  ! A hack so that we can manually deallocate the mesh, and hence deallocate the
  ! faces - James was too busy to go digging in vtk_read_state
  call incref(psi%mesh)
  
  call deallocate(psi)
  call deallocate(matrix)  
  call deallocate(rhs)
  call deallocate(state)
  call deallocate(psi%mesh)
  
  call print_references(0)
  
  ewrite(1, *) "Exiting streamfunction_2d"
  
contains

  subroutine assemble_streamfunction_2d_element(ele, matrix, rhs, positions, velocity)
    integer, intent(in) :: ele
    type(csr_matrix), intent(inout) :: matrix
    type(scalar_field), intent(inout) :: rhs
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: velocity
    
    integer, dimension(:), pointer :: element_nodes
    real, dimension(ele_ngi(rhs, ele)) :: detwei, vorticity_gi
    real, dimension(ele_loc(rhs, ele), ele_ngi(rhs, ele), mesh_dim(rhs)) :: dn_t
    real, dimension(ele_loc(velocity, ele), ele_ngi(rhs, ele), mesh_dim(rhs)) :: du_t
    type(element_type), pointer :: positions_shape, psi_shape, velocity_shape
    
    assert(ele_ngi(positions, ele) == ele_ngi(rhs, ele))
    assert(ele_ngi(velocity, ele) == ele_ngi(rhs, ele))
    
    positions_shape => ele_shape(positions, ele)
    psi_shape => ele_shape(rhs, ele)
    velocity_shape => ele_shape(velocity, ele)
    
    call transform_to_physical(positions, ele, psi_shape, &
      & dshape = dn_t, detwei = detwei)
     
    assert(sum(abs(detwei)) > epsilon(0.0))
      
    if(psi_shape == velocity_shape) then
      du_t = dn_t
    else
      call transform_to_physical(positions, ele, velocity_shape, &
        & dshape = du_t)
    end if
    
    vorticity_gi = ele_2d_curl_at_quad(velocity, ele, du_t)
    
    element_nodes => ele_nodes(rhs, ele)
    
    call addto(matrix, element_nodes, element_nodes, dshape_dot_dshape(dn_t, dn_t, detwei))
    call addto(rhs, element_nodes, shape_rhs(psi_shape, vorticity_gi * detwei))
    
  end subroutine assemble_streamfunction_2d_element

end subroutine streamfunction_2d
