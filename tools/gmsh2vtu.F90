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

subroutine gmsh2vtu(filename_, filename_len) bind(c)
  !!< Read in a gmsh mesh and output a vtu mesh.
  
  use fields
  use read_gmsh
  use vtk_interfaces
  use iso_c_binding
  implicit none
  
  integer(kind=c_size_t), value :: filename_len
  character(kind=c_char, len=1) :: filename_(*)
  
  character(len=filename_len) :: filename
  
  integer :: stat, i
  type(vector_field), target :: positions
  type(scalar_field) :: regions

  do i=1, filename_len
    filename(i:i)=filename_(i)
  end do

  positions=read_gmsh_file(filename, quad_degree=3)

  if (associated(positions%mesh%region_ids)) then
     regions=piecewise_constant_field(positions%mesh, name="Regions")
     regions%val=float(positions%mesh%region_ids)
     call vtk_write_fields(filename, position=positions, &
          model=positions%mesh, vfields=(/positions/), sfields=(/ regions /))
     call deallocate(regions)
  else
    call vtk_write_fields(filename, position=positions, model=positions%mesh)
  end if
  
  call deallocate(positions)

end subroutine gmsh2vtu
