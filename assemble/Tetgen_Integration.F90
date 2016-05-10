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
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
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

module tetgen_integration

  use iso_c_binding
  use fldebug
  use futils, only: present_and_true
  use data_structures
  use quadrature
  use elements
  use spud
  use parallel_tools
  use fields
  use halos
  use meshdiagnostics
  use tictoc

  implicit none
  
  private

  public tetgenerate_mesh

  contains

    subroutine tetgenerate_mesh(input_positions, output_positions)

      type(vector_field), intent(in) :: input_positions
      type(vector_field), target, intent(out) :: output_positions

      integer :: n_surface_elements
      integer, parameter :: dim = 3, nloc = 4, snloc = 3
      integer, dimension(:), allocatable :: snlist

      type(mesh_type), pointer :: output_mesh
      integer(kind=c_int) :: neles
      type(c_ptr) :: context

      interface 
         subroutine tetgen_cleanup(context) bind(c)
           use iso_c_binding
           implicit none
           type(c_ptr), value :: context
         end subroutine tetgen_cleanup

         subroutine set_from_tetgenio(neles,element_list,context) bind(c)
           use iso_c_binding
           implicit none
           integer(c_int), value, intent(in) :: neles
           integer(c_int) :: element_list(4*neles)
           type(c_ptr), value :: context
         end subroutine set_from_tetgenio
      end interface

      n_surface_elements = surface_element_count(input_positions)

      allocate(snlist(n_surface_elements * snloc))
      if(n_surface_elements > 0) then
         call getsndgln(input_positions%mesh, snlist)
      end if

#ifdef HAVE_LIBTET
      call pass_to_tetgen(input_positions%val, snlist, neles, context)
#endif

      allocate(output_mesh)
      call allocate(output_mesh, node_count(input_positions), &
           neles, input_positions%mesh%shape, &
           name = input_positions%mesh%name)
      output_mesh%shape%refcount%tagged = .false.
      output_mesh%shape%quadrature%refcount%tagged = .false.
    
      call set_from_tetgenio(neles, output_mesh%ndglno, context) 
      output_mesh%option_path = input_positions%mesh%option_path 

      ! Construct the new positions
      call allocate(output_positions, dim, output_mesh, name = input_positions%name)
      call deallocate(output_mesh)
      deallocate(output_mesh)
      output_mesh => output_positions%mesh
      output_positions%option_path = input_positions%option_path
      output_positions%val(:,:) = input_positions%val

      call add_faces(output_mesh, sndgln = snlist, boundary_ids = input_positions%mesh%faces%boundary_ids)
      if(associated(input_positions%mesh%faces%coplanar_ids)) then
         allocate(output_mesh%faces%coplanar_ids(n_surface_elements))
         output_mesh%faces%coplanar_ids = input_positions%mesh%faces%coplanar_ids
      end if

      deallocate(snlist)
#ifdef HAVE_LIBTET
      call tetgen_cleanup(context)
#endif

    end subroutine tetgenerate_mesh

    subroutine pass_to_tetgen(nodes, facets, neles, context)

      real(kind=c_double), intent(in)  :: nodes(:,:)
      integer(kind=c_int), intent(in) :: facets(:)
      integer(kind=c_int), intent(out) :: neles
      type(c_ptr) :: context


      interface 
         subroutine tetgen_binding(nnodes,nodes,nfacets,&
                                         facets,nelements, context) bind(c)
           use iso_c_binding
           implicit none
           integer(kind=c_int), value, intent(in) :: nnodes, nfacets
           real(kind=c_double), intent(in) :: nodes(3,nnodes)
           integer(kind=c_int), intent(in) :: facets(3*nfacets)
           integer(kind=c_int) :: nelements
           type(c_ptr), intent(out) :: context
         end subroutine tetgen_binding
      end interface 

#ifdef HAVE_LIBTET
      call tetgen_binding(size(nodes,2), nodes,&
                          size(facets)/3, facets,&
                          neles, context)
#else
      FLExit("Fluidity compiled without tetgen support")
#endif

    end subroutine pass_to_tetgen

  end module tetgen_integration


                          
