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
  use futils, only: present_and_true, int2str
  use data_structures
  use quadrature
  use elements
  use spud
  use fefields
  use parallel_tools
  use fields
  use halos
  use meshdiagnostics
  use tictoc

  implicit none
  
  private

  type, bind(c) :: mesh_data
     integer(c_int) :: nnodes, nelements, nfacets, nholes
     integer(c_int) :: lnode_ids, lregion_ids, lface_ids
     type(c_ptr)    :: nodes 
     type(c_ptr)    :: node_ids
     type(c_ptr)    :: ndglno
     type(c_ptr)    :: region_ids
     type(c_ptr)    :: facets
     type(c_ptr)    :: face_ids
     type(c_ptr)    :: holes
  end type mesh_data

  public tetgenerate_mesh

  contains

    subroutine tetgenerate_mesh(input_positions, output_positions)

      type(vector_field), intent(in) :: input_positions
      type(vector_field), target, intent(out) :: output_positions

      integer :: nhalos
      type(halo_type), pointer :: old_halo, new_halo
      integer :: i, proc
      integer, dimension(:), allocatable :: renumber_permutation

      integer :: stat
      integer, dimension(2) :: shape_option
      integer, dimension(:), allocatable :: region_id_nos

      ewrite(1,*) "In tetgenerate_mesh"

      stat = 1
      !! Work out which regions are to be remeshed
      if( have_option("/mesh_adaptivity/delaunay_adaptivity/region_ids")) then
         shape_option=option_shape("/mesh_adaptivity/delaunay_adaptivity/region_ids")
         allocate(region_id_nos(1:shape_option(1)))
         call get_option("/mesh_adaptivity/delaunay_adaptivity/region_ids",&
              region_id_nos, stat)
      elseif(.not. have_option("/mesh_adaptivity/delaunay_adaptivity/ignore_region_ids")) then
         call count_regions(input_positions%mesh,region_id_nos)
      else
         allocate(region_id_nos(1))
         region_id_nos=-1
      end if

      nhalos = halo_count(input_positions)
      assert(any(nhalos == (/0, 1, 2/)))
      if(nhalos > 0) then
         old_halo => input_positions%mesh%halos(nhalos)
      end if

      call external_remesh(input_positions, output_positions,region_id_nos)

      deallocate(region_id_nos)

      if(nhalos > 0) then

         allocate(output_positions%mesh%halos(nhalos))
         new_halo => output_positions%mesh%halos(nhalos)
      
         ! halo is the same in terms of n/o sends and receives, name, ordering_type, numbering, etc.
         call allocate(new_halo, old_halo)
         ! except (maybe) for n/o owned nodes
         call set_halo_nowned_nodes(new_halo, node_count(output_positions)-halo_all_receives_count(new_halo))
         do proc=1, halo_proc_count(new_halo)
            call set_halo_sends(new_halo, proc, halo_sends(old_halo, proc))
         end do
         do proc=1, halo_proc_count(new_halo)
            call set_halo_receives(new_halo, proc, halo_receives(old_halo, proc))
         end do
         
         allocate(renumber_permutation(element_count(output_positions)))
      
         if(nhalos == 2) then
            ! Derive remaining halos
            call derive_l1_from_l2_halo(output_positions%mesh, ordering_scheme = HALO_ORDER_GENERAL, create_caches = .false.)
            ! Reorder the nodes for trailing receives consistency
            call renumber_positions_trailing_receives(output_positions)
            
            allocate(output_positions%mesh%element_halos(2))
            ! Reorder the elements for trailing receives consistency
            call derive_element_halo_from_node_halo(output_positions%mesh, &
                 & ordering_scheme = HALO_ORDER_GENERAL, create_caches = .false.)
            call renumber_positions_elements_trailing_receives(output_positions, permutation=renumber_permutation)
         else
            ! Reorder the nodes for trailing receives consistency
            call renumber_positions_trailing_receives(output_positions)
        
            allocate(output_positions%mesh%element_halos(1))
            ! Reorder the elements for trailing receives consistency
            call derive_element_halo_from_node_halo(output_positions%mesh, &
                 & ordering_scheme = HALO_ORDER_GENERAL, create_caches = .false.)
            call renumber_positions_elements_trailing_receives(output_positions, permutation=renumber_permutation)
         end if

         deallocate(renumber_permutation)

         ! Adaptivity is not guaranteed to return halo elements in the same
         ! order in which they went in. We therefore need to fix this order.
         call reorder_element_numbering(output_positions)

#ifdef DDEBUG
         do i = 1, nhalos
            assert(trailing_receives_consistent(output_positions%mesh%halos(i)))
            assert(halo_valid_for_communication(output_positions%mesh%halos(i)))
            assert(halo_verifies(output_positions%mesh%halos(i), output_positions))
         end do
#endif
      end if

      ewrite(1,*) "Leaving tetgenerate_mesh"

    end subroutine tetgenerate_mesh

    subroutine external_remesh(input_positions, output_positions, regions)
      !! This routine does the hard work in extracting the regions to be
      !! remeshed, then reconstructing the final output after tetgen has
      !! been called.
      !! The input positional vector field
      type(vector_field), intent(in) :: input_positions
      !! The output positional vector field
      type(vector_field), intent(out), target :: output_positions
      !! An array of region ids to be considered, or -1 to skip,
      integer, intent(in), dimension(:) :: regions

      integer :: dim, nloc, snloc, reg
      type(mesh_data) :: input, output(size(regions))
      integer(c_int), pointer, dimension(:,:) :: idata_2d
      real(c_double), pointer, dimension(:,:) :: rdata_2d, holes

      integer :: i, j, n_surface_elements, nholes, fixed_eles, nelements
      integer, dimension(size(regions)) :: extra_nodes
      integer, dimension(:), pointer :: snlist
      type int_ptr
         integer, dimension(:), pointer :: ptr
      end type int_ptr
      type(int_ptr), dimension(size(regions)) :: node_list

      real, dimension(:,:), pointer :: pnodes

      type(mesh_type), pointer :: output_mesh
      type(mesh_type) :: submesh(size(regions))
      logical, dimension(element_count(input_positions)) :: validity
      logical :: use_submesh

      interface
         
         function tetgen( mesh, command) bind(c)
           use iso_c_binding
           implicit none
           type, bind(c) :: mesh_data
              integer(c_int) :: nnodes, nelements, nfacets, nholes
              integer(c_int) :: lnode_ids, lregion_ids, lface_ids
              type(c_ptr)    :: nodes 
              type(c_ptr)    :: node_ids
              type(c_ptr)    :: ndglno
              type(c_ptr)    :: region_ids
              type(c_ptr)    :: facets
              type(c_ptr)    :: face_ids
              type(c_ptr)    :: holes
           end type mesh_data
           type(mesh_data) :: mesh
           character(kind=c_char) :: command(*)
           type(mesh_data) :: tetgen
         end function tetgen

          function triangle( mesh, command) bind(c)
           use iso_c_binding

           implicit none
           type, bind(c) :: mesh_data
              integer(c_int) :: nnodes, nelements, nfacets, nholes
              integer(c_int) :: lnode_ids, lregion_ids, lface_ids
              type(c_ptr)    :: nodes 
              type(c_ptr)    :: node_ids
              type(c_ptr)    :: ndglno
              type(c_ptr)    :: region_ids
              type(c_ptr)    :: facets
              type(c_ptr)    :: face_ids
              type(c_ptr)    :: holes
           end type mesh_data
           type(mesh_data) :: mesh
           character(kind=c_char) :: command(*)
           type(mesh_data) :: triangle
         end function triangle

         subroutine mesh_data_cleanup(mesh, dim) bind(c)
           use iso_c_binding
           implicit none
           type, bind(c) :: mesh_data
              integer(c_int) :: nnodes, nelements, nfacets, nholes
              integer(c_int) :: lnode_ids, lregion_ids, lface_ids
              type(c_ptr)    :: nodes 
              type(c_ptr)    :: node_ids
              type(c_ptr)    :: ndglno
              type(c_ptr)    :: region_ids
              type(c_ptr)    :: facets
              type(c_ptr)    :: face_ids
              type(c_ptr)    :: holes
           end type mesh_data
           type(mesh_data) :: mesh
           integer(c_int) :: dim
         end subroutine mesh_data_cleanup

      end interface

      dim = mesh_dim(input_positions)
      nloc = ele_loc(input_positions,1)
      snloc = face_loc(input_positions,1)

      ewrite(1,*) 'Remeshing regions:', regions

      if (isparallel() .and. .false.) then
         nholes=0
         allocate(holes(dim,0))
      else
         nholes = option_count("/mesh_adaptivity/delaunay_adaptivity/mesh_hole")
         allocate(holes(dim, nholes))
         do i=0,nholes-1
            call get_option("/mesh_adaptivity/delaunay_adaptivity/mesh_hole["&
                 //int2str(i)//"]",holes(:,i+1))
         end do
      end if

      n_surface_elements = surface_element_count(input_positions)

      use_submesh = isparallel() .or. (regions(1) >= 0)

      if (use_submesh) then
         call get_validity(input_positions%mesh, regions, validity)
      end if

      fixed_eles = element_count(input_positions)
      do reg = 1, size(regions) 
         if (use_submesh) then
            call get_free_subdomain(input_positions%mesh,submesh(reg),&
                 validity, node_list(reg)%ptr, snlist,regions(reg))
            allocate(pnodes(dim,node_count(submesh(reg))))
            pnodes = input_positions%val(:,node_list(reg)%ptr)
            input = pack_mesh_data(pnodes,submesh(reg), snlist, holes)
            fixed_eles = fixed_eles&
                 - element_count(submesh(reg))
         else
            allocate(snlist(n_surface_elements * snloc))
            if(n_surface_elements > 0) then
               call getsndgln(input_positions%mesh, snlist)
            end if
            input = pack_mesh_data(input_positions%val, input_positions%mesh,&
                 snlist, holes)
            fixed_eles = 0
         end if

         select case(dim)
         case(2)
#ifdef HAVE_LIBTRIANGLE
            output(reg) = triangle(input, "YYYCS0ps5ss20"//C_NULL_CHAR)
#else
            FLExit("Fluidity compiled without triangle support")
#endif
         case(3)
#ifdef HAVE_LIBTET
            output(reg) = tetgen(input, "YYYpABMFO2/1"//C_NULL_CHAR)
#else
            FLExit("Fluidity compiled without tetgen support")
#endif
         end select

         extra_nodes(reg) = output(reg)%nnodes - input%nnodes
         deallocate(snlist)
         if (use_submesh) deallocate(pnodes)
      end do

      ewrite(1,*) "nodes:", node_count(input_positions)+sum(extra_nodes)
      nelements = count_interior_elements(output)+fixed_eles
      ewrite(1,*) "elements:", nelements
      
      allocate(output_mesh)
      call allocate(output_mesh,&
           node_count(input_positions)+sum(extra_nodes), &
           nelements, input_positions%mesh%shape, &
           name = input_positions%mesh%name)
      output_mesh%shape%refcount%tagged = .false.
      output_mesh%shape%quadrature%refcount%tagged = .false.
      if (regions(1)>=0) then
         allocate(output_mesh%region_ids(element_count(output_mesh)))
      end if
      if (use_submesh) then
         do reg=1, size(regions)
            call c_f_pointer(output(reg)%ndglno, idata_2d,&
                 [nloc,output(reg)%nelements])
            call add_free_elements(input_positions%mesh,&
                 output(reg), output_mesh,&
                 count_interior_elements(output(1:reg-1)),&
                 node_count(submesh(reg)),&
                 node_count(input_positions)+sum(extra_nodes(1:reg-1)),&
                 idata_2d, node_list(reg)%ptr, regions(reg))
            call deallocate(submesh(reg))
            deallocate(node_list(reg)%ptr)
         end do
         call add_fixed_elements(input_positions%mesh,&
                 output_mesh, count_interior_elements(output), validity)
      else
         j=0
         call c_f_pointer(output(1)%ndglno, idata_2d,&
              [nloc,output(1)%nelements])
         call add_free_elements(input_positions%mesh,&
              output(1), output_mesh,&
              0, node_count(input_positions),&
              node_count(input_positions),&
              idata_2d, node_list(1)%ptr, -1)
      end if
      output_mesh%option_path = input_positions%mesh%option_path 

      ! Construct the new positions
      call allocate(output_positions, dim, output_mesh, &
           name = input_positions%name)
      call deallocate(output_mesh)
      deallocate(output_mesh)
      output_mesh => output_positions%mesh
      output_positions%option_path = input_positions%option_path

      !!! The preexising nodes go back where they came from
      output_positions%val(:,1:node_count(input_positions))&
           = input_positions%val
      j=node_count(input_positions)
      do reg=1,size(regions)
         call c_f_pointer(output(reg)%nodes, rdata_2d, [dim, output(reg)%nnodes])
         !! extra nodes are put on the end for now.
         output_positions%val(:,j+1:j+extra_nodes(reg))&
              = rdata_2d(:,output(reg)%nnodes-extra_nodes(reg)+1:output(reg)%nnodes)
         j = j+extra_nodes(reg)
      end do

      allocate(snlist(n_surface_elements * snloc))
      call getsndgln(input_positions%mesh,&
           snlist(1:n_surface_elements * snloc))
      call add_faces(output_mesh,&
           sndgln = snlist,&
           boundary_ids = input_positions%mesh%faces%boundary_ids)
      if(associated(input_positions%mesh%faces%coplanar_ids)) then
         allocate(output_mesh%faces%coplanar_ids(n_surface_elements))
         output_mesh%faces%coplanar_ids = input_positions%mesh%faces%coplanar_ids
      end if

      assert(n_surface_elements == surface_element_count(output_positions))

      deallocate(snlist)

      do reg=1,size(regions)
         call mesh_data_cleanup(output(reg), dim)
      end do

    end subroutine external_remesh

    function pack_mesh_data(nodes, inmesh, snlist, holes) result(mesh)
      !! put the data for the mesh into a C compatible type
      real, dimension(:,:), pointer :: nodes
      type(mesh_type), intent(in) :: inmesh
      integer(c_int), dimension(:), pointer :: snlist
      real(c_double), dimension(:,:), pointer :: holes
      type(mesh_data) :: mesh

      mesh%nnodes    = node_count(inmesh)
      mesh%nelements = element_count(inmesh)
      mesh%nfacets   = size(snlist)/face_loc(inmesh,1)
      mesh%nholes = size(holes,2)

      mesh%lnode_ids = 0
      mesh%lregion_ids = 0
      mesh%lface_ids = 0

      mesh%nodes  = c_loc(nodes)
      mesh%node_ids = c_null_ptr
      mesh%ndglno = c_loc(inmesh%ndglno)
      mesh%region_ids = c_null_ptr
      mesh%facets = c_loc(snlist)
      mesh%face_ids = c_null_ptr
      if (mesh%nholes>0) then      
         mesh%holes = c_loc(holes)
      else
         mesh%holes = c_null_ptr
      end if
      

    end function pack_mesh_data

    subroutine get_validity(xmesh, regions, validity)
      ! full mesh to take free elements from
      type(mesh_type), intent(in) :: xmesh
      integer, dimension(:), intent(in) :: regions
      ! elements that will make up the free submesh
      logical, dimension(:) :: validity

      integer :: ele, nhalos, k
      integer, dimension(:), pointer :: nodes

      validity = .true.

      if (regions(1)>=0) then
         do ele=1,element_count(xmesh)
            if (.not. any(regions==xmesh%region_ids(ele))) then
               validity(ele) = .false.
            end if
         end do
      end if

      nhalos = halo_count(xmesh)
      if (nhalos>0) then
         do ele=1,element_count(xmesh)
            nodes => ele_nodes(xmesh,ele)
            do k=1, size(nodes)
               if (.not. node_owned(xmesh%halos(nhalos), nodes(k))) then
                  validity(ele) = .false.
                  exit
               end if
            end do
         end do
      end if
    end subroutine get_validity
    
    subroutine get_free_subdomain(xmesh, submesh, validity, node_list, snlist, region_id)
      ! full mesh to take submesh from
      type(mesh_type), intent(in) :: xmesh
      ! submesh created
      type(mesh_type), intent(out) :: submesh
      ! elements that will make up the submesh
      logical, dimension(:) :: validity
      ! list of nodes in submesh (also functions as node map from submesh to full mesh)
      integer, dimension(:), pointer :: node_list
      ! surface face list
      integer, dimension(:), pointer :: snlist
      ! region_id to pick out
      integer, intent(in) :: region_id
     
      integer :: k
      
      if (region_id>=0) then
         call create_subdomain_mesh(xmesh,&
              pack([(k,k=1,element_count(xmesh))],&
                     validity .and. xmesh%region_ids==region_id),&
              "Free region", submesh, node_list)
      else
         call create_subdomain_mesh(xmesh,&
              pack([(k,k=1,element_count(xmesh))],validity),&
              "Free subdomain", submesh, node_list)
      end if
         
      call get_surfaces(submesh,snlist)
      
    end subroutine get_free_subdomain

    subroutine get_surfaces(xmesh, snlist)
      type(mesh_type), intent(in) :: xmesh
      integer, dimension(:), pointer :: snlist

      integer :: ele, face, k, opface, sloc
      integer :: face_list(face_count(xmesh))
      integer, dimension(:), pointer :: faces

      sloc = face_loc(xmesh,1)

      face_list = 0

      do ele=1,element_count(xmesh)
         faces => ele_faces(xmesh, ele)
         do k=1, size(faces)
               opface = face_opposite(xmesh, faces(k))
               if (opface<0) face_list(faces(k)) = 1
         end do
      end do


     allocate(snlist(3*(count(face_list==1))))

     k=0

      do face=1,size(face_list)
         if (face_list(face)==1) then
            snlist(sloc*k+1:sloc*(k+1)) = face_global_nodes(xmesh, face)
            k=k+1
         end if
      end do

    end subroutine get_surfaces

    subroutine add_free_elements(xmesh, imesh, output_mesh, offset,&
         free_nodes, original_nodes, data, node_list, region)
      type(mesh_type), intent(in) :: xmesh
      type(mesh_data), intent(in) :: imesh
      type(mesh_type), intent(inout) :: output_mesh
      integer, intent(in) :: offset, free_nodes, original_nodes, region
      integer, intent(in), dimension(:,:) :: data
      integer, intent(in), dimension(:) :: node_list

      logical, dimension(imesh%nelements) :: interior
      
      integer :: i, j, k, nloc

      print*, 'offset:', offset

      nloc =ele_loc(xmesh,1)

      interior = get_interior_elements(imesh)

      k=0
      do j = 1,imesh%nelements
         if (.not. interior(j)) cycle
         do i=1,nloc
            if (data(i,j)<=size(node_list)) then
               output_mesh%ndglno(i+nloc*(offset+k))=node_list(data(i,j))
            else
               output_mesh%ndglno(i+nloc*(offset+k))=data(i,j)&
                    -free_nodes+original_nodes
            end if
         end do
         k=k+1 
      end do

      if (associated(output_mesh%region_ids)) then
         output_mesh%region_ids(offset+1:offset+count(interior)) = region
      end if
      
    end subroutine add_free_elements

    subroutine add_fixed_elements(xmesh, output_mesh, offset, validitylist)

      type(mesh_type), intent(in) :: xmesh
      type(mesh_type), intent(inout) :: output_mesh
      integer, intent(in) :: offset
      logical, intent(in), dimension(:) :: validitylist

      integer :: j, ele, nloc

      nloc =ele_loc(xmesh,1)
      j = offset
      if (j == element_count(output_mesh)) return

      do ele=1,element_count(xmesh)
         if (validitylist(ele)) cycle
         output_mesh%ndglno(nloc*j+1:nloc*(j+1))=ele_nodes(xmesh,ele)
         if (associated(output_mesh%region_ids)) then
            output_mesh%region_ids(j+1) = xmesh%region_ids(ele)
         end if
         j = j + 1
      end do

    end subroutine add_fixed_elements

    subroutine count_regions(mesh,region_id_nos)
      type(mesh_type), intent(in) :: mesh
      integer, dimension(:), allocatable, intent(out) :: region_id_nos

      integer :: ele, j
      integer, dimension(element_count(mesh)) :: tmp_region_id_nos

      j=0

      do ele=1, element_count(mesh)
         if (.not. any(tmp_region_id_nos(1:j)==mesh%region_ids(ele))) then
            j = j+1
            tmp_region_id_nos(j) = mesh%region_ids(ele)
         end if
      end do

      allocate(region_id_nos(j))
      region_id_nos = tmp_region_id_nos(1:j)
    end subroutine count_regions

    function count_interior_elements(mesh) result(out)
      type(mesh_data), intent(in), dimension(:) :: mesh
      integer :: i, out

      out = 0
      do i=1, size(mesh)
         out = out + count(get_interior_elements(mesh(i)))
      end do

    end function count_interior_elements
    
    real function get_max_region(mesh)

      type(mesh_data), intent(in) :: mesh

      real(c_double), pointer, dimension(:,:) :: region_ids

      call c_f_pointer(mesh%region_ids,region_ids,&
           [mesh%lregion_ids,mesh%nelements])

      get_max_region = maxval(region_ids(1,:))

    end function get_max_region

    function get_interior_elements(mesh) result(out)
      type(mesh_data), intent(in) :: mesh
      logical, dimension(mesh%nelements) :: out

      real(c_double), pointer, dimension(:,:) :: region_ids

      call c_f_pointer(mesh%region_ids,region_ids,&
           [mesh%lregion_ids,mesh%nelements])
      out = region_ids(1,:)>0.0

     end function get_interior_elements
      

  end module tetgen_integration


                          
