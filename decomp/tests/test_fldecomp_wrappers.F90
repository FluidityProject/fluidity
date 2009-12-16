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

subroutine test_fldecomp_wrappers
  !!< Tests wrapped routines used by fldecomp
  
  use fields
  use fields_data_types
  use read_triangle
  use unittest_tools
  use vtk_interfaces
  
  implicit none
  
  interface
    subroutine query_mesh(filename, filename_len, meshtype, meshtype_len, dim, nnodes, nelements, nloc, snelements, snloc)
      implicit none
      integer, intent(in) :: filename_len
      integer, intent(in) :: meshtype_len
      character(len = filename_len), intent(in) :: filename
      character(len = meshtype_len), intent(in) :: meshtype
      integer, intent(out) :: dim
      integer, intent(out) :: nnodes
      integer, intent(out) :: nelements
      integer, intent(out) :: nloc
      integer, intent(out) :: snelements
      integer, intent(out) :: snloc
    end subroutine query_mesh
    
    function read_mesh(filename, filename_len, meshtype, meshtype_len, x, dim, nnodes, enlist, region_ids, nelements, nloc, senlist, boundary_ids, snelements, snloc)
      implicit none
      integer, intent(in) :: filename_len
      integer, intent(in) :: meshtype_len
      integer, intent(in) :: dim
      integer, intent(in) :: nnodes
      integer, intent(in) :: nelements
      integer, intent(in) :: nloc
      integer, intent(in) :: snelements
      integer, intent(in) :: snloc
      character(len = filename_len), intent(in) :: filename
      character(len = meshtype_len), intent(in) :: meshtype
      real, dimension(nnodes * dim), intent(out) :: x
      integer, dimension(nelements * nloc), intent(out) :: enlist
      integer, dimension(nelements), intent(out) :: region_ids
      integer, dimension(snelements * snloc), intent(out) :: senlist
      integer, dimension(snelements), intent(out) :: boundary_ids
      integer :: read_mesh
    end function read_mesh
    
    function write_mesh(filename, filename_len, meshtype, meshtype_len, x, dim, nnodes, enlist, region_ids, nelements, nloc, senlist, boundary_ids, snelements, snloc)
      implicit none
      integer, intent(in) :: filename_len
      integer, intent(in) :: meshtype_len
      integer, intent(in) :: dim
      integer, intent(in) :: nnodes
      integer, intent(in) :: nelements
      integer, intent(in) :: nloc
      integer, intent(in) :: snelements
      integer, intent(in) :: snloc
      character(len = filename_len), intent(in) :: filename
      character(len = meshtype_len), intent(in) :: meshtype
      real, dimension(nnodes * dim), intent(in) :: x
      integer, dimension(nelements * nloc), intent(in) :: enlist
      integer, dimension(nelements), intent(in) :: region_ids
      integer, dimension(snelements * snloc), target, intent(in) :: senlist
      integer, dimension(snelements), intent(in) :: boundary_ids
      integer :: write_mesh
    end function write_mesh
  end interface
  
  integer :: i
  type(vector_field) :: mesh_field
  
  ! read/write_mesh variables
  character(len = 255) :: filename, meshtype
  integer :: dim, dim_read, nelements, nelements_read, nloc, nloc_read, nnodes, nnodes_read, snelements, snelements_read, snloc, snloc_read
  integer, dimension(:), allocatable :: boundary_ids, boundary_ids_read, enlist, enlist_read, region_ids, region_ids_read, senlist, senlist_read
  real, dimension(:), allocatable :: x, x_read
  
  ! Read a mesh
  mesh_field = read_triangle_files("data/square-cavity", quad_degree = 1)
  
  filename = "data/test_fldecomp_wrappers_out"
  meshtype = "triangle"
    
  ! Extract the mesh data
  dim = mesh_field%dim
  nnodes = node_count(mesh_field)
  nelements = ele_count(mesh_field)
  nloc = ele_loc(mesh_field, 1)
  snelements = surface_element_count(mesh_field)
  snloc = face_loc(mesh_field, 1)
  allocate(enlist(nelements * nloc))
  allocate(region_ids(nelements))
  allocate(senlist(snelements * snloc))
  allocate(boundary_ids(snelements))
  allocate(x(nnodes * dim))
  do i = 1, nnodes 
    x((i - 1) * dim + 1:i * dim) = node_val(mesh_field, i)
  end do
  enlist = mesh_field%mesh%ndglno
  if(associated(mesh_field%mesh%region_ids)) then
    region_ids = mesh_field%mesh%region_ids
  else
    region_ids = 0
  end if
  call getsndgln(mesh_field%mesh, senlist)
  boundary_ids = mesh_field%mesh%faces%boundary_ids

  ! Write the mesh
  call report_test("[write_mesh]", &
    & write_mesh(trim(filename), len_trim(filename), trim(meshtype), len_trim(meshtype), &
    & x, dim, nnodes, enlist, region_ids, nelements, nloc, senlist, boundary_ids, snelements, snloc) /= 0, &
    & .false., "Failed to write mesh")
  
  ! Query the written mesh
  call query_mesh(trim(filename), len_trim(filename), trim(meshtype), len_trim(meshtype), &
    & dim_read, nnodes_read, nelements_read, nloc_read, snelements_read, snloc_read)
  call report_test("[Query mesh returned correct dim]", dim_read /= dim, .false., "Returned dim = " // int2str(dim_read))
  call report_test("[Query mesh returned correct nnodes]", nnodes_read /= nnodes, .false., "Returned nnodes = " // int2str(nnodes_read))
  call report_test("[Query mesh returned correct nloc]", nloc_read /= nloc, .false., "Returned nloc = " // int2str(nloc_read))
  call report_test("[Query mesh returned correct nelements]", nelements_read /= nelements, .false., "Returned nelements = " // int2str(nelements_read))
  call report_test("[Query mesh returned correct snloc]", snloc_read /= snloc, .false., "Returned snloc = " // int2str(snloc_read))
  call report_test("[Query mesh returned correct snelements]", snelements_read /= snelements, .false., "Returned snelements = " // int2str(snelements_read))
  
  ! Read the written mesh
  allocate(x_read(nnodes_read * dim_read))
  allocate(enlist_read(nelements_read * nloc_read))
  allocate(region_ids_read(nelements_read))
  allocate(senlist_read(snelements_read * snloc_read))
  allocate(boundary_ids_read(snelements_read))
  call report_test("[read_mesh]", &
    & read_mesh(trim(filename), len_trim(filename), trim(meshtype), len_trim(meshtype), &
    & x_read, dim_read, nnodes_read, enlist_read, region_ids_read, nelements_read, nloc_read, senlist_read, boundary_ids_read, snelements_read, snloc_read) /= 0, &
    & .false., "Failed to read mesh")
  call report_test("[Correct x read]", x_read .fne. x, .false., "Incorrect x read")
  call report_test("[Correct enlist read]", count(enlist_read /= enlist) /= 0, .false., "Incorrect enlist read")
  call report_test("[Correct region_ids read]", count(region_ids_read /= region_ids) /= 0, .false., "Incorrect region_ids read")
  call report_test("[Correct senlist read]", count(senlist_read /= senlist) /= 0, .false., "Incorrect senlist read")
  call report_test("[Correct boundary_ids read]", count(boundary_ids_read /= boundary_ids) /= 0, .false., "Incorrect boundary_ids read")
  
  deallocate(x)
  deallocate(enlist)
  deallocate(region_ids)
  deallocate(senlist)
  deallocate(boundary_ids)
  deallocate(x_read)
  deallocate(enlist_read)
  deallocate(region_ids_read)
  deallocate(senlist_read)
  deallocate(boundary_ids_read)
  call deallocate(mesh_field)
  
end subroutine test_fldecomp_wrappers
