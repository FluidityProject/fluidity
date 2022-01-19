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

subroutine vtu2gmsh(filename_, filename_len) bind(c)
  !!< Read in a vtu and output a gmsh mesh.

  use fields
  use state_module
  use mesh_files
  use vtk_interfaces
  use iso_c_binding
  implicit none

  character(kind=c_char, len=1) :: filename_(*)
  integer(kind=c_size_t), value :: filename_len

  character(len=filename_len) :: filename
  type(vector_field), pointer :: positions
  type(state_type) :: state
  integer :: i

  do i=1, filename_len
    filename(i:i)=filename_(i)
  end do

  call vtk_read_state(filename // ".vtu", state)
  positions => extract_vector_field(state, "Coordinate")
  call add_faces(positions%mesh)
  call write_mesh_files(filename, format="gmsh", positions=positions)

  call deallocate(state)

end subroutine vtu2gmsh
