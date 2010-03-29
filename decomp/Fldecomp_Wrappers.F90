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

subroutine query_mesh(filename, filename_len, meshtype, meshtype_len, dim, nnodes, nelements, nloc, snelements, snloc)
  !!< Determine properties of the supplied mesh file

  use elements
  use fields
  use fields_data_types
  use read_triangle

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

  type(element_type) :: shape
  type(quadrature_type) :: quad

  ewrite(1, *) "In query_mesh"

  select case(meshtype)
    case("triangle")
      call identify_triangle_file(filename, dim = dim, loc = nloc, nodes = nnodes, elements = nelements, selements = snelements)
      quad = make_quadrature(nloc, dim, degree = 1)
      shape = make_element_shape(nloc, dim, 1, quad)
      call deallocate(quad)
      select case(ele_numbering_family(shape))
        case(FAMILY_SIMPLEX)
          if(nloc /= dim + 1) then
            ewrite(0, *) "Warning: triangle boundary markers not supported for quadratic space elements"
            assert(snelements == 0)
            snloc = 0
          else
            snloc = nloc - 1
          end if
        case(FAMILY_CUBE)
          snloc = nloc / 2
        case default
          ewrite(-1, "(a,i0)") "For element family ", shape%numbering%family
          FLAbort("Unrecognised element family")
      end select
      call deallocate(shape)
    case default
      FLAbort("Invalid mesh type supplied")
  end select

  ewrite(1, *) "Exiting query_mesh"

end subroutine query_mesh

function read_mesh(filename, filename_len, meshtype, meshtype_len, x, dim, nnodes, enlist, region_ids, nelements, nloc, senlist, boundary_ids, snelements, snloc)
  !!< A wrapper to enable fldecomp to read meshes without needing to port code to C++

  use fields
  use fields_data_types
  use read_triangle
    
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
  real, dimension(nnodes * dim) :: x
  integer, dimension(nelements * nloc), intent(out) :: enlist
  integer, dimension(nelements), intent(out) :: region_ids
  integer, dimension(snelements * snloc), intent(out) :: senlist
  integer, dimension(snelements), intent(out) :: boundary_ids

  integer :: read_mesh

  integer :: i
  type(vector_field) :: mesh_field
  
  ewrite(1, *) "In read_mesh"
  
  ! Read the mesh.
  select case(meshtype)
    case("triangle")
      ! Use a quadrature degree of one as we won't be integrating
      mesh_field = read_triangle_files(filename, quad_degree = 1)
    case default
      FLAbort("Invalid mesh type supplied")
  end select
    
  ! Input check (doesn't protect against bad calls)
  assert(mesh_field%dim == dim)
  assert(node_count(mesh_field) == nnodes)
  assert(ele_count(mesh_field) == nelements)
  assert(ele_loc(mesh_field, 1) == nloc)
  if(surface_element_count(mesh_field) /= snelements) then
    ewrite(-1, "(a,i0)") "Expected surface elements: ", surface_element_count(mesh_field)
    ewrite(-1, "(a,i0)") "Surface elements found: ", snelements
    FLExit("All surface elements must be marked")
  end if
  assert(face_loc(mesh_field, 1) == snloc)

  ! Copy data into format expected by fldecomp
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
  
  read_mesh = 0
  
  ! Deallocate memory used to read the mesh
  call deallocate(mesh_field)
  
  ewrite(1, *) "Exiting read_mesh"

end function read_mesh

function write_mesh(filename, filename_len, meshtype, meshtype_len, x, dim, nnodes, enlist, region_ids, nelements, nloc, senlist, boundary_ids, snelements, snloc)
  !!< A wrapper to enable fldecomp to write meshes without needing to port code to C++
  
  use fields
  use fields_data_types
  use write_triangle
  
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
  integer, dimension(snelements * snloc), intent(in) :: senlist
  integer, dimension(snelements), intent(in) :: boundary_ids
  
  integer :: write_mesh
  
  integer :: i
  type(element_type) :: shape
  type(mesh_type) :: mesh
  type(quadrature_type) :: quad
  type(vector_field) :: mesh_field
  
  ewrite(1, *) "In write_mesh"
        
  ! Allocate the Coordinate field
  ! Use a quadrature degree of one as we won't be integrating.
  quad = make_quadrature(nloc, dim, degree = 1)
  shape = make_element_shape(nloc, dim, 1, quad)
  call allocate(mesh, nnodes, nelements, shape, name = filename)
  call allocate(mesh_field, dim, mesh, "Coordinate")

  ! Copy data into mesh
  do i = 1, nnodes
    call set(mesh_field, i, x((i - 1) * dim + 1:i * dim))
  end do
  mesh_field%mesh%ndglno = enlist
  allocate(mesh_field%mesh%region_ids(size(region_ids)))
  mesh_field%mesh%region_ids = region_ids
  ! incomplete_surface_mesh=.true. to prevent legacy behaviour of
  ! "completing" the surface mesh (this would add faces at the end of the halo
  ! to the surface mesh)
  call add_faces(mesh_field%mesh, sndgln = senlist, boundary_ids = boundary_ids, &
    incomplete_surface_mesh=.true.)
  
  ! Write the mesh
  select case(meshtype)
    case("triangle")
      call write_triangle_files(filename, mesh_field)
    case default
      FLAbort("Invalid mesh type supplied")
  end select
  
  ! Deallocate memory used to write the mesh
  call deallocate(quad)
  call deallocate(shape)
  call deallocate(mesh)
  call deallocate(mesh_field)
  
  write_mesh = 0
  
  ewrite(1, *) "Exiting write_mesh"
  
end function write_mesh
