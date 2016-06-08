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
  use tetgen_data_types
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
  use pickers

  implicit none
  
  private

  public tetgenerate_mesh

  contains

    subroutine tetgenerate_mesh(input_positions, output_positions)

      type(vector_field), intent(inout) :: input_positions
      type(vector_field), target, intent(out) :: output_positions

      integer :: nhalos
      type(halo_type), pointer :: old_halo, new_halo
      integer :: i, proc
      integer, dimension(:), allocatable :: renumber_permutation

      integer :: stat
      integer, dimension(2) :: shape_option
      integer, dimension(:), allocatable :: region_id_nos
      logical :: keep_regions, renumber_nodes

      ewrite(1,*) "In tetgenerate_mesh"

      stat = 1
      !! Work out which regions are to be remeshed
!!      if( have_option("/mesh_adaptivity/delaunay_adaptivity/region_ids")) then
!!         shape_option=option_shape("/mesh_adaptivity/delaunay_adaptivity/region_ids")
!!         allocate(region_id_nos(1:shape_option(1)))
!!         call get_option("/mesh_adaptivity/delaunay_adaptivity/region_ids",&
!!              region_id_nos, stat)
      if(.not. have_option("/mesh_adaptivity/delaunay_adaptivity/ignore_region_ids")) then
!!         call count_regions(input_positions%mesh,region_id_nos)
         keep_regions = .true.
      else
!         allocate(region_id_nos(1))
         keep_regions = .false.
      end if

      nhalos = halo_count(input_positions)
      assert(any(nhalos == (/0, 1, 2/)))
      if(nhalos > 0) then
         old_halo => input_positions%mesh%halos(nhalos)
      end if

      call external_remesh(input_positions, output_positions, keep_regions)

      renumber_nodes = node_count(input_positions) /= node_count(output_positions)

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


            if (renumber_nodes) then
               ! Reorder the nodes for trailing receives consistency
               call renumber_positions_trailing_receives(output_positions)
            else
               do i=1, nhalos
                  call create_ownership(output_positions%mesh%halos(i))
                  call create_global_to_universal_numbering(output_positions%mesh%halos(i))
               end do
            end if
            
            allocate(output_positions%mesh%element_halos(2))
            ! Reorder the elements for trailing receives consistency
            call derive_element_halo_from_node_halo(output_positions%mesh, &
                 & ordering_scheme = HALO_ORDER_GENERAL)
 !           call renumber_positions_elements_trailing_receives(output_positions, permutation=renumber_permutation)
         else
            if (renumber_nodes) then
               ! Reorder the nodes for trailing receives consistency
               call renumber_positions_trailing_receives(output_positions)
            else
               do i=1, nhalos
                  call create_ownership(output_positions%mesh%halos(i))
                  call create_global_to_universal_numbering(output_positions%mesh%halos(i))
               end do
            end if

        
            allocate(output_positions%mesh%element_halos(1))
            ! Reorder the elements for trailing receives consistency
            call derive_element_halo_from_node_halo(output_positions%mesh, &
                 & ordering_scheme = HALO_ORDER_GENERAL)
!            call renumber_positions_elements_trailing_receives(output_positions, permutation=renumber_permutation)
         end if

         deallocate(renumber_permutation)

         ! Adaptivity is not guaranteed to return halo elements in the same
         ! order in which they went in. We therefore need to fix this order.
!         call reorder_element_numbering(output_positions)

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

    subroutine external_remesh(input_positions, output_positions, keep_regions)
      !! This routine does the hard work in extracting the regions to be
      !! remeshed, then reconstructing the final output after tetgen has
      !! been called.
      !! The input positional vector field
      type(vector_field), intent(inout) :: input_positions
      !! The output positional vector field
      type(vector_field), intent(out), target :: output_positions
      !! An array of region ids to be considered, or -1 to skip,
      logical, intent(in) :: keep_regions

      integer :: dim, nloc, snloc, reg
      type(mesh_data) :: input, output, temp
      integer(c_int), pointer, dimension(:,:) :: idata_2d
      real(c_double), pointer, dimension(:,:) :: rdata_2d, holes

      integer :: i, j, n_surface_elements, nholes, fixed_eles, nelements
      integer, dimension(:), pointer :: snlist
      integer(c_int), dimension(:), pointer :: eelist, faceelist, faces
      real(c_double), pointer, dimension(:) :: region_attributes=> null()

      real, dimension(:,:), pointer :: pnodes

      type(mesh_type), pointer :: output_mesh
      integer, dimension(element_count(input_positions)) :: validity
      logical :: use_submesh

      interface
         
         function tetgen( mesh, command) bind(c)
           use iso_c_binding
           use tetgen_data_types

           implicit none

           type(mesh_data) :: mesh
           character(kind=c_char) :: command(*)
           type(mesh_data) :: tetgen
         end function tetgen

          function triangle( mesh, command) bind(c)
           use iso_c_binding
           use tetgen_data_types

           implicit none

           type(mesh_data) :: mesh
           character(kind=c_char) :: command(*)
           type(mesh_data) :: triangle
         end function triangle

         subroutine mesh_data_cleanup(mesh, dim) bind(c)
           use iso_c_binding
           use tetgen_data_types

           implicit none

           type(mesh_data) :: mesh
           integer(c_int) :: dim
         end subroutine mesh_data_cleanup

      end interface

      dim = mesh_dim(input_positions)
      nloc = ele_loc(input_positions,1)
      snloc = face_loc(input_positions,1)

      nholes = option_count("/mesh_adaptivity/delaunay_adaptivity/mesh_hole")
      allocate(holes(dim, nholes))
      do i=0,nholes-1
         call get_option("/mesh_adaptivity/delaunay_adaptivity/mesh_hole["&
              //int2str(i)//"]",holes(:,i+1))
      end do

      n_surface_elements = surface_element_count(input_positions)

      use_submesh = isparallel() .or. keep_regions 

      if (use_submesh) then
         call get_validity(input_positions%mesh, keep_regions, validity)
         fixed_eles = count(validity == 0)
         call get_free_subdomain(input_positions%mesh, validity, snlist)
         call paint_regions(input_positions, validity, region_attributes)
         input = pack_mesh_data(input_positions%val,input_positions%mesh,&
              snlist, holes, region_attributes)
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
         output = triangle(input, "YYYCS0ps5ss20"//C_NULL_CHAR)
#else
         FLExit("Fluidity compiled without triangle support")
#endif
      case(3)
#ifdef HAVE_LIBTET
         output = tetgen(input, "YYYS0pABMfnnO2/1"//C_NULL_CHAR)
#else
         FLExit("Fluidity compiled without tetgen support")
#endif
      end select

      deallocate(snlist)

      ewrite(1,*) "nodes:", output%nnodes
      nelements = count_interior_elements([output])+fixed_eles
      ewrite(1,*) "elements:", nelements
      ewrite(1,*) "free elements:", count_interior_elements([output])
      
      allocate(output_mesh)
      call allocate(output_mesh,&
           node_count(input_positions), &
           nelements, input_positions%mesh%shape, &
           name = input_positions%mesh%name)
      output_mesh%shape%refcount%tagged = .false.
      output_mesh%shape%quadrature%refcount%tagged = .false.
      if (keep_regions) then
         allocate(output_mesh%region_ids(element_count(output_mesh)))
      end if
      call c_f_pointer(output%ndglno, idata_2d,&
           [nloc,output%nelements])
      call add_free_elements(input_positions%mesh,&
           output, output_mesh,&
           node_count(input_positions),&
           node_count(input_positions),&
           idata_2d)
      if (use_submesh) then
         call add_fixed_elements(input_positions%mesh,&
              output_mesh, count_interior_elements([output]), validity)
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
      call c_f_pointer(output%nodes, rdata_2d, [dim, output%nnodes])
      !! extra nodes are put on the end for now.
      output_positions%val(:,j+1:node_count(output_positions))&
           = rdata_2d(:,j+1:output%nnodes)

      allocate(snlist(n_surface_elements * snloc))
      call getsndgln(input_positions%mesh,&
           snlist(1:n_surface_elements * snloc))
      if (isparallel()) then
         call mesh_data_cleanup(output, dim)
         input = pack_mesh_data(output_positions%val, output_positions%mesh)
#ifdef HAVE_LIBTET
         output = tetgen(input, "YYYS0rBNMfnn"//C_NULL_CHAR)
#else
         FLExit("Fluidity compiled without tetgen support")
#endif
 !                 call add_faces(output_mesh,&
 !             sndgln = snlist,&
 !             boundary_ids = input_positions%mesh%faces%boundary_ids)
      end if

         call c_f_pointer(output%element_adjacency, eelist,&
              [nloc*element_count(output_positions)])
         call c_f_pointer(output%faces_adjacency, faceelist,&
              [2*output%nfaces])
         call c_f_pointer(output%faces, faces,&
              [3*output%nfaces])
         call add_faces(output_mesh,&
              sndgln = snlist,&
              boundary_ids = input_positions%mesh%faces%boundary_ids,&
              known_faces = faces,&
              face_adjacency = faceelist,&
              known_eelist = eelist)
!      end if
      if(associated(input_positions%mesh%faces%coplanar_ids)) then
         allocate(output_mesh%faces%coplanar_ids(n_surface_elements))
         output_mesh%faces%coplanar_ids = input_positions%mesh%faces%coplanar_ids
      end if

      assert(n_surface_elements == surface_element_count(output_positions))

      deallocate(snlist, holes)

      call mesh_data_cleanup(output, dim)

    end subroutine external_remesh

    subroutine match_regions(positions, validity, outmesh)

      type(vector_field), intent(inout) :: positions
      integer, dimension(:), intent(in) :: validity
      type(mesh_data), intent(inout) :: outmesh

      real(c_double), pointer, dimension(:,:) :: region_ids
      integer(c_int), pointer, dimension(:,:) :: ndglno

      integer :: max_id, ele, ele_b, nloc, i
      real, dimension(positions%dim) :: coord
      

      integer, dimension(:), allocatable:: map 

      nloc = positions%dim + 1

      call c_f_pointer(outmesh%region_ids,region_ids,&
           [outmesh%lregion_ids,outmesh%nelements])
      call c_f_pointer(outmesh%ndglno, ndglno,&
           [nloc,outmesh%nelements])

      max_id = int(maxval(region_ids(1,:)))

      !! only one region of interest
      
      if (max_id==1) return

      allocate(map(max_id))

      do i=1, max_id
         do ele = 1, outmesh%nelements
            if (int(region_ids(1,ele)) == i) then
               coord = get_centre(ele)
               call picker_inquire(positions, coord, ele_b, global=.false.)
               map(i) = validity(ele_b)
               exit
            end if
         end do
      end do

      do ele = 1, outmesh%nelements
         region_ids(1,ele)= map(nint(region_ids(1,ele)))
      end do

      contains

        function get_centre(e)
          integer, intent(in) :: e

          real, dimension(positions%dim) :: get_centre
          integer :: j

          get_centre = 0

          do j=1, nloc
             get_centre = get_centre + node_val(positions, ndglno(j,e))
          end do

          get_centre = get_centre / nloc

        end function get_centre

    end subroutine match_regions

    subroutine paint_regions(positions, validity, region_attributes)
      type(vector_field), intent(inout) :: positions
      integer, dimension(:), intent(in) :: validity
      real(c_double), dimension(:), pointer :: region_attributes

      integer :: ncoloured, ncols, nfront1,nfront2, ele
      integer, dimension(element_count(positions)) :: ele_colour, front_list1,&
           front_list2, ele_num
      integer, dimension(:), pointer :: nodes
      logical :: do_colouring

      ele_colour = 0
      ncoloured = 0
      ncols = 0
      
      do_colouring = .true.

      if (.not. isParallel() .and. associated(positions%mesh%region_ids)) then
         if (size(positions%mesh%region_ids)==1) then
            ncols =0
            do_colouring = .false.
         end if
      end if   

      if (do_colouring) then
         do while (ncoloured<element_count(positions))
            nfront1=0
            front_list1=0
            do ele=1, element_count(positions)
               if (ele_colour(ele) == 0) then
                  ncols = ncols + 1
                  ele_colour(ele) = ncols
                  ncoloured = ncoloured + 1
                  ele_num(ncols) = ele               
                  call paint_neighs(ele)
                  exit
               end if
            end do
            do while(nfront1>0) 
               nfront2 = nfront1
               front_list2 = front_list1
               nfront1=0
               front_list1=0
               do ele=1,nfront2
                  call paint_neighs(front_list2(ele))
               end do
            end do
         end do
      end if

      allocate(region_attributes(5*ncols))

      do ele=1, ncols
         region_attributes(5*(ele-1)+1:5*(ele-1)+3) = &
              sum(ele_val(positions,ele_num(ele)), dim=2)/4.0
         region_attributes(5*(ele-1)+4) = 1.0*validity(ele_num(ele))
         region_attributes(5*(ele-1)+5) = 0.0
      end do
        
      contains

        subroutine paint_neighs(e)

          integer, intent(in) :: e

          integer :: i
          integer, dimension(:), pointer :: neigh

          neigh => ele_neigh(positions%mesh, e)
          do i=1,size(neigh)
             if (neigh(i)<0) cycle
             if (ele_colour(neigh(i))>0) cycle
             if (validity(e) /= validity(neigh(i))) cycle
             
             ele_colour(neigh(i)) = ncols
             ncoloured = ncoloured + 1
             nfront1 = nfront1+1
             front_list1(nfront1) = neigh(i)
          end do

        end subroutine paint_neighs
        
    end subroutine paint_regions

    function pack_mesh_data(nodes, inmesh, snlist, holes, region_attributes) result(mesh)
      !! put the data for the mesh into a C compatible type
      real, dimension(:,:), pointer :: nodes
      type(mesh_type), intent(in) :: inmesh
      integer(c_int), dimension(:), pointer, optional :: snlist
      real(c_double), dimension(:,:), pointer, optional :: holes
      real(c_double), dimension(:), pointer, optional :: region_attributes
      type(mesh_data) :: mesh
      integer(c_int), pointer :: iptr
      real(c_double), pointer :: dptr

      mesh%nnodes    = node_count(inmesh)
      mesh%nelements = element_count(inmesh)
      if (present(snlist)) then
         mesh%nfacets = size(snlist)/face_loc(inmesh,1)
      else
         mesh%nfacets = 0
      end if
      if (present(holes)) then
         mesh%nholes = size(holes,2)
      else
         mesh%nholes = 0
      end if

      mesh%lnode_ids = 0
      if (present(region_attributes)) then
         mesh%nattributes = size(region_attributes)/5
         if (size(region_attributes)>0) then
            dptr => region_attributes(1)
         else
            nullify(dptr) 
         end if
         mesh%region_attributes = c_loc(dptr)
      else
         mesh%nattributes = 0
         mesh%region_attributes = c_null_ptr
      end if
      mesh%lregion_ids = 0
      mesh%lface_ids = 0
      
      dptr => nodes(1,1)
      mesh%nodes  = c_loc(dptr)
      mesh%node_ids = c_null_ptr
      iptr => inmesh%ndglno(1)
      mesh%ndglno = c_loc(iptr)
      mesh%region_ids = c_null_ptr
      if (mesh%nfacets>0) then
         iptr => snlist(1)
         mesh%facets = c_loc(iptr)
      else
         mesh%facets = c_null_ptr
      end if
      mesh%face_ids = c_null_ptr
      if (mesh%nholes>0) then 
         dptr => holes(1,1)
         mesh%holes = c_loc(dptr)
      else
         mesh%holes = c_null_ptr
      end if
      mesh%element_adjacency = c_null_ptr

      mesh%nfaces = 0
      mesh%faces = c_null_ptr
      mesh%faces_adjacency = c_null_ptr

    end function pack_mesh_data

    subroutine get_validity(xmesh, keep_regions, validity)
      ! full mesh to take free elements from
      type(mesh_type), intent(in) :: xmesh
      logical, intent(in) :: keep_regions
      ! elements that will make up the free submesh
      integer, dimension(:) :: validity

      integer :: ele, nhalos, k
      integer, dimension(:), pointer :: nodes

      if (keep_regions) then
         validity = xmesh%region_ids
      else
         validity = 1
      end if

      nhalos = halo_count(xmesh)
      if (nhalos>0) then
         do ele=1,element_count(xmesh)
            nodes => ele_nodes(xmesh,ele)
            do k=1, size(nodes)
               if (.not. node_owned(xmesh%halos(nhalos), nodes(k))) then
                  validity(ele) = 0
                  exit
               end if
            end do
         end do
      end if

    end subroutine get_validity
    
    subroutine get_free_subdomain(xmesh, validity, snlist)
      ! full mesh to take submesh from
      type(mesh_type), intent(in) :: xmesh
      ! elements that will make up the regions of the submesh
      integer, dimension(:) :: validity
      integer, dimension(:), pointer :: snlist

      integer :: ele, face, k, opface, sloc
      integer :: face_list(face_count(xmesh))
      integer, dimension(:), pointer :: faces, neigh

      sloc = face_loc(xmesh,1)

      face_list = 0

      do ele=1,element_count(xmesh)
         faces => ele_faces(xmesh, ele)
         neigh => ele_neigh(xmesh, ele)
         do k=1, size(faces)
               if (ele<neigh(k)) cycle
               if (neigh(k)<0) then
                  face_list(faces(k)) = 1
                  cycle
               end if
               if (validity(ele) /= validity(neigh(k))) then
                  face_list(faces(k)) = 1
                  cycle
               end if
         end do
      end do


     allocate(snlist(sloc*(count(face_list==1))))

     k=0

      do face=1,size(face_list)
         if (face_list(face)==1) then
            snlist(sloc*k+1:sloc*(k+1)) = face_global_nodes(xmesh, face)
            k=k+1
         end if
      end do

    end subroutine get_free_subdomain

    subroutine add_free_elements(xmesh, imesh, output_mesh,&
         free_nodes, original_nodes, data)
      type(mesh_type), intent(in) :: xmesh
      type(mesh_data), intent(in) :: imesh
      type(mesh_type), intent(inout) :: output_mesh
      integer, intent(in) :: free_nodes, original_nodes
      integer, intent(in), dimension(:,:) :: data

      logical, dimension(imesh%nelements) :: interior
      
      integer :: i, j, k, nloc

      nloc =ele_loc(xmesh,1)

      interior = get_interior_elements(imesh)

      k=0
      do j = 1,imesh%nelements
         if (.not. interior(j)) cycle
         do i=1,nloc
            output_mesh%ndglno(i+nloc*k)=data(i,j)
         end do
         k=k+1 
      end do

      assert(k==count(interior))

      if (associated(output_mesh%region_ids)) then
         output_mesh%region_ids(1:count(interior)) = 1
      end if
      
    end subroutine add_free_elements

    subroutine add_fixed_elements(xmesh, output_mesh, offset, validitylist)

      type(mesh_type), intent(in) :: xmesh
      type(mesh_type), intent(inout) :: output_mesh
      integer, intent(in) :: offset
      integer, intent(in), dimension(:) :: validitylist

      integer :: j, ele, nloc

      nloc =ele_loc(xmesh,1)
      j = offset
      if (j == element_count(output_mesh)) return

      do ele=1,element_count(xmesh)
         if (validitylist(ele) > 0) cycle
         output_mesh%ndglno(nloc*j+1:nloc*(j+1)) = ele_nodes(xmesh,ele)
         if (associated(output_mesh%region_ids)) then
            output_mesh%region_ids(j+1) = xmesh%region_ids(ele)
         end if
         j = j + 1
      end do

      assert( j == element_count(output_mesh))

    end subroutine add_fixed_elements

    subroutine count_regions(mesh,region_id_nos)
      type(mesh_type), intent(in) :: mesh
      integer, dimension(:), allocatable, intent(out) :: region_id_nos

      integer :: ele, j
      integer, dimension(element_count(mesh)) :: tmp_region_id_nos

      if (.not. associated(mesh%region_ids)) then
         allocate(region_id_nos(1))
         region_id_nos(1) = -1
         return
      end if

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


                          
