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

module distance_diagnostics

  use iso_c_binding
  use diagnostic_source_fields
  use field_options
  use fields_manipulation
  use initialise_fields_module
  use fields
  use fldebug
  use global_parameters, only : timestep, OPTION_PATH_LEN, current_time
  use spud
  use boundary_conditions
  use state_fields_module
  use state_module
  use parallel_tools

implicit none

  private

  public  :: calculate_scalar_distance_function
  public  :: calculate_vector_distance_function



  contains

subroutine calculate_scalar_distance_function(states,state_index,s_field)

  type(state_type), dimension(:) :: states
  integer, intent(in) :: state_index
  type(scalar_field), intent(inout):: s_field

  integer :: i, stat, diagnostic_count, n 
  type(vector_field) :: X
  type(scalar_field) :: p1_field
  logical :: diagnostic
  real :: zero_level
  character( len = OPTION_PATH_LEN ) :: algorithm 

  integer :: sele, snodes, seles
  integer, dimension(2) :: shape
  integer, dimension(:), allocatable :: surface_ids
  real :: max_distance

  type(mesh_type), pointer :: surface_mesh
  integer, dimension(:), pointer :: surface_node_list

  real, dimension(:,:), allocatable :: surface_X
  integer, dimension(:), allocatable :: sndglno

  interface

     subroutine vtk_distance(dim, nodes, snodes, eles, seles, &
          pts, ndglno,  surface_X, sndglno, distance) bind(c)
       use iso_c_binding
       integer(c_int), value :: dim, nodes, snodes, eles, seles
       real(c_double) :: pts(dim*nodes), distance(nodes), surface_X(dim*snodes)
       integer(c_int) :: ndglno((dim+1)*eles), sndglno(dim*seles)
     end subroutine vtk_distance
  end interface

  if (get_boundary_condition_count(s_field) == 0) then
     shape = option_shape(trim(s_field%option_path) // "/diagnostic/algorithm/surface_ids")
     assert(shape(1) >= 0)
     allocate(surface_ids(shape(1)))
     call get_option(trim(s_field%option_path) // "/diagnostic/algorithm/surface_ids", surface_ids)
     call add_boundary_condition(s_field, 'wall', 'distance', surface_ids)
  end if

  call get_boundary_condition(s_field, 'wall', &
       surface_node_list=surface_node_list,&
       surface_mesh=surface_mesh)

  X = get_coordinate_field(states(state_index), s_field%mesh)

  call serialise_boundary(X, surface_mesh, surface_node_list, surface_X, sndglno, snodes, seles)

  call vtk_distance(X%dim, node_count(X),&
       snodes, element_count(X),&
       seles, X%val,&
       X%mesh%ndglno, surface_X, sndglno, &
       s_field%val)

  call deallocate(X)
  
end subroutine calculate_scalar_distance_function

subroutine calculate_vector_distance_function(states,state_index,v_field)

  type(state_type), dimension(:) :: states
  integer, intent(in) :: state_index
  type(vector_field), intent(inout):: v_field

  integer :: i, stat, diagnostic_count, n 
  type(vector_field) :: X
  type(scalar_field) :: p1_field
  logical :: diagnostic
  real :: zero_level
  character( len = OPTION_PATH_LEN ) :: algorithm 

  integer :: sele, snodes, seles
  integer, dimension(2) :: shape
  integer, dimension(:), allocatable :: surface_ids
  real :: max_distance

  type(mesh_type), pointer :: surface_mesh
  integer, dimension(:), pointer :: surface_node_list

  real, dimension(:,:), allocatable :: surface_X
  integer, dimension(:), allocatable :: sndglno
  

  interface

     subroutine vtk_direction(dim, nodes, snodes, eles, seles, &
          pts, ndglno,  spts, sndglno, direction) bind(c)
       use iso_c_binding
       integer(c_int), value :: dim, nodes, snodes, eles, seles
       real(c_double) :: pts(dim*nodes), direction(dim*nodes), spts(dim*snodes)
       integer(c_int) :: ndglno((dim+1)*eles), sndglno(dim*seles)
     end subroutine vtk_direction
  end interface

  if (get_boundary_condition_count(v_field) == 0) then
     shape = option_shape(trim(v_field%option_path) // "/diagnostic/algorithm/surface_ids")
     assert(shape(1) >= 0)
     allocate(surface_ids(shape(1)))
     call get_option(trim(v_field%option_path) // "/diagnostic/algorithm/surface_ids", surface_ids)
     call add_boundary_condition(v_field, 'wall', 'distance', surface_ids)
  end if

  call get_boundary_condition(v_field, 'wall', &
       surface_node_list=surface_node_list,&
       surface_mesh=surface_mesh)

  X = get_coordinate_field(states(state_index), v_field%mesh)

  call serialise_boundary(X, surface_mesh, surface_node_list, surface_X, sndglno, snodes, seles)

  call vtk_direction(X%dim, node_count(X),&
       snodes, element_count(X),&
       seles, X%val,&
       X%mesh%ndglno,surface_X, sndglno, &
       v_field%val)

  call deallocate(X) 
  
end subroutine calculate_vector_distance_function

subroutine serialise_boundary(X,mesh,node_list, surface_X, sndglno, nodes, eles)
  type(vector_field) :: X
  type(mesh_type) :: mesh
  integer, dimension(:) :: node_list
  real, dimension(:,:), allocatable :: surface_X
  integer, dimension(:), allocatable :: sndglno
  integer :: nodes, eles

  integer, dimension(:), allocatable :: offset, recvsize
  integer :: i, s, ierr
  real, dimension(:,:), allocatable :: X_loc
  integer, dimension(:), allocatable :: sndglno_loc
 
  allocate(offset(getnprocs()), recvsize(getnprocs()), &
       X_loc(X%dim, size(node_list)), &
       sndglno_loc(element_count(mesh)*ele_loc(mesh,1)))

  nodes = node_count(mesh)

  s = mesh_dim(X)*node_count(mesh)
  call mpi_allgather(s,1,getpinteger(),recvsize,&
       1,getpinteger(), MPI_COMM_FEMTOOLS, ierr)
  offset(1)=0
  do i=1,getnprocs()-1
     offset(i+1) = offset(i) + recvsize(i)
  end do

  eles = element_count(mesh)  
  call allsum(eles)
  allocate(sndglno(eles*ele_loc(mesh,1)))
  call allsum(nodes)
  allocate(surface_X(mesh_dim(X),nodes))

  x_loc = X%val(:,node_list)
  sndglno_loc = mesh%ndglno + offset(getprocno())/mesh_dim(X)

  call mpi_allgatherv(X_loc, s , getpreal(),&
       surface_X,recvsize, offset, getpreal(), MPI_COMM_FEMTOOLS, ierr)



  s  = element_count(mesh) * ele_loc(mesh,1)
  call mpi_allgather(s,1,getpinteger(),recvsize,&
       1,getpinteger(), MPI_COMM_FEMTOOLS, ierr)
  offset(1)=0
  do i=1,getnprocs()-1
     offset(i+1) = offset(i) + recvsize(i)
  end do

  call mpi_allgatherv(sndglno_loc, s , getpinteger(),&
       sndglno,recvsize, offset, getpinteger(), MPI_COMM_FEMTOOLS, ierr)
  

end subroutine serialise_boundary

end module distance_diagnostics
    
