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

subroutine triangle2vtu(filename_, filename_len) bind(c)
  !!< Read in a triangle mesh and output a vtu mesh.
  
  use fields
  use read_triangle
  use vtk_interfaces
  use iso_c_binding
  implicit none
  
  integer(kind=c_size_t), value :: filename_len
  character(kind=c_char, len=1) :: filename_(*)
  
  character(len=filename_len) :: filename
  
  integer :: stat, i
  type(vector_field), target :: positions
  type(scalar_field) :: mapA, mapB, regions

  do i=1, filename_len
    filename(i:i)=filename_(i)
  end do

  positions=read_triangle_files(filename, quad_degree=3, no_faces=.true.)

  ! For supermesh stuff.
  ! It tests for the existence of the mapping files
  ! mapping from elements in the supermesh to elements
  ! in the original mesh.
  ! If they're not there, nothing changes in the output.
  mapA = read_elemental_mappings(positions, filename, "mapCA", stat = stat)
  if(stat == 0) mapB = read_elemental_mappings(positions, filename, "mapCB", stat = stat)
  
  if (stat == 0) then
    call vtk_write_fields(filename, position=positions, &
         model=positions%mesh, sfields=(/mapA, mapB/), vfields=(/positions/))
  else if (associated(positions%mesh%region_ids)) then
    regions=piecewise_constant_field(positions%mesh, name="Regions")
    regions%val=float(positions%mesh%region_ids)
    call vtk_write_fields(filename, position=positions, &
         model=positions%mesh, vfields=(/positions/), sfields=(/ regions /))
    call deallocate(regions)
  else
    call vtk_write_fields(filename, position=positions, model=positions%mesh)
  end if
  
  call deallocate(positions)
  if (associated(mapA%val)) then
    call deallocate(mapA)
  end if
  if (associated(mapB%val)) then
    call deallocate(mapB)
  end if

end subroutine triangle2vtu
