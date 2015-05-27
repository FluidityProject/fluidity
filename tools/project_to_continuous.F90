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
#include "confdefs.h"
#include "fdebug.h"

subroutine project_to_continuous(vtuname_, vtuname_len, meshname_,&
     & meshname_len) bind(c)
  !!< Given a vtu file containing fields on a discontinuous mesh, and the
  !!< mesh files for the corresponding continuous mesh, produce a vtu
  !!< with its fields projected onto the continuous mesh.
  use state_module
  use elements
  use fields
  use mesh_files
  use vtk_interfaces
  use sparse_tools
  use fefields
  use sparse_matrices_fields
  use iso_c_binding
  implicit none

  character(kind=c_char, len=1) :: vtuname_(*)
  integer(kind=c_size_t), value :: vtuname_len
  character(kind=c_char, len=1) :: meshname_(*)
  integer(kind=c_size_t), value :: meshname_len

  character(len=vtuname_len):: vtuname
  character(len=meshname_len):: meshname

  type(state_type) :: dg_state, cg_state
  type(vector_field) :: cg_coordinate
  type(csr_matrix) :: P
  type(scalar_field) :: lumped_mass
  type(mesh_type) :: dg_mesh, cg_mesh
  type(scalar_field) :: dg_scalar, cg_scalar
  type(vector_field) :: dg_vector, cg_vector
  type(tensor_field) :: dg_tensor, cg_tensor

  integer :: i, j, k

! now turn into proper fortran strings (is there an easier way to do this?)
  do i=1, vtuname_len
    vtuname(i:i)=vtuname_(i)
  end do
  do i=1, meshname_len
    meshname(i:i)=meshname_(i)
  end do

  call vtk_read_state(vtuname, dg_state, quad_degree=6)
  
  cg_coordinate= read_mesh_files(meshname, quad_degree=6, format="gmsh")
  cg_mesh=cg_coordinate%mesh

  call allocate(lumped_mass, cg_mesh, "LumpedMass")
  
  call compute_lumped_mass(cg_coordinate, lumped_mass)
  ! Invert lumped mass.
  lumped_mass%val=1./lumped_mass%val

  dg_mesh=extract_mesh(dg_state, "Mesh")
  
  P=compute_projection_matrix(cg_mesh, dg_mesh, cg_coordinate)
  
  do i=1,size(dg_state%scalar_fields)
     dg_scalar=dg_state%scalar_fields(i)%ptr
     call allocate(cg_scalar, cg_mesh, name=dg_scalar%name)
     
     ! Perform projection.
     call mult(cg_scalar, P, dg_scalar)
     ! Apply inverted lumped mass to projected quantity.
     call scale(cg_scalar, lumped_mass)

     call insert(cg_state, cg_scalar, cg_scalar%name)
     
     ! Drop the additional reference.
     call deallocate(cg_scalar)
  end do

  do i=1,size(dg_state%vector_fields)
     dg_vector=dg_state%vector_fields(i)%ptr
     call allocate(cg_vector, dg_vector%dim, cg_mesh, name=dg_vector%name)
     
     ! Perform projection.
     do j=1,cg_vector%dim
        cg_scalar=extract_scalar_field_from_vector_field(cg_vector, j)
        dg_scalar=extract_scalar_field_from_vector_field(dg_vector, j)
        call mult(cg_scalar, P, dg_scalar)
     end do
     ! Apply inverted lumped mass to projected quantity.
     call scale(cg_vector, lumped_mass)

     call insert(cg_state, cg_vector, cg_vector%name)
     
     ! Drop the additional reference.
     call deallocate(cg_vector)
  end do

  do i=1,size(dg_state%tensor_fields)
     dg_tensor=dg_state%tensor_fields(i)%ptr
     call allocate(cg_tensor, cg_mesh, name=dg_tensor%name)
     
     ! Perform projection.
     do j=1,cg_tensor%dim(1)
        do k=1,cg_tensor%dim(2)
           cg_scalar=extract_scalar_field_from_tensor_field(cg_tensor, j, k)
           dg_scalar=extract_scalar_field_from_tensor_field(dg_tensor, j, k)
           call mult(cg_scalar, P, dg_scalar)
        end do
     end do
     ! Apply inverted lumped mass to projected quantity.
     call scale(cg_tensor, lumped_mass)

     call insert(cg_state, cg_tensor, cg_tensor%name)
     
     ! Drop the additional reference.
     call deallocate(cg_tensor)
  end do

  ! We do this insertion last because otherwise we end up with a projection
  ! of coordinate and that gives us a wonky mesh.
  call insert(cg_state, cg_coordinate, "Coordinate")
  call insert(cg_state, cg_mesh, "CoordinateMesh")

  call vtk_write_state(vtuname(1:len_trim(vtuname)-4)//"_continuous", &
       state=(/cg_state/))

end subroutine project_to_continuous
