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
module global_numbering
  ! **********************************************************************
  ! Module to construct the global node numbering map for elements of a
  ! given degree.
  use adjacency_lists
  use elements
  use sparse_tools
  use fldebug
  use halo_data_types
  use halos_allocates
  use halos_base
  use halos_debug
  use halos_numbering
  use halos_ownership
  use parallel_tools
  use integer_set_module
  use linked_lists
  use mpi_interfaces
  use fields_base
  
  implicit none

  private
  
  public :: make_global_numbering_DG, make_boundary_numbering,&
       & make_global_numbering, element_halo_communicate_visibility

contains
  
  subroutine make_global_numbering_DG(new_nonods, new_ndglno, Totele,&
       & element, element_halos, new_halos)
    ! Construct a global node numbering for the solution variables in a
    ! Discontinuous Galerkin simulation. This is trivial.
    !
    ! Note that this code is broken for mixed element meshes.
    integer, intent(in) :: totele
    type(element_type), intent(in) :: element

    integer, dimension(:), intent(out) :: new_ndglno
    integer, intent(out) :: new_nonods
    type(halo_type), dimension(:), intent(in), optional :: element_halos
    type(halo_type), dimension(:), intent(out), optional :: new_halos

    integer :: i

    new_nonods=totele*element%loc

    forall (i=1:new_nonods)
       new_ndglno(i)=i
    end forall

    if (.not.present(element_halos)) return
    assert(present(new_halos))
    assert(size(element_halos)==size(new_halos))

    do i=1,size(new_halos)
       call make_halo_dg(element, element_halos(i), new_halos(i))
    end do

  contains
    
    subroutine make_halo_dg(element, element_halo, new_halo)
      !!< This routine constructs a node halo given an element halo.
      type(element_type), intent(in) :: element
      type(halo_type), intent(in), optional :: element_halo
      type(halo_type), intent(out), optional :: new_halo
      
      integer, dimension(size(element_halo%sends)) :: nsends
      integer, dimension(size(element_halo%receives)) :: nreceives
      
      integer :: i,j,k, nloc

      nloc=element%loc

      do i=1, size(nsends)
         nsends(i)=nloc*size(element_halo%sends(i)%ptr)
      end do
      do i=1, size(nreceives)
         nreceives(i)=nloc*size(element_halo%receives(i)%ptr)
      end do

      call allocate(new_halo, &
              nsends, &
              nreceives, &
!! Query what is the naming convention for halos.
              name=trim(halo_name(element_halo)) // "DG", &
              nprocs=element_halo%nprocs, &
              nowned_nodes=element_halo%nowned_nodes*element%loc, &
              data_type=HALO_TYPE_DG_NODE, &
              ordering_scheme=halo_ordering_scheme(element_halo))

      do i=1, size(nsends)
         do j=1,size(element_halo%sends(i)%ptr)
            
            new_halo%sends(i)%ptr((j-1)*nloc+1:j*nloc)&
                 =(element_halo%sends(i)%ptr(j)-1)*nloc + (/(k,k=1,nloc)/)

         end do
      end do

      do i=1, size(nreceives)
         do j=1,size(element_halo%receives(i)%ptr)
            
            new_halo%receives(i)%ptr((j-1)*nloc+1:j*nloc)&
                 =(element_halo%receives(i)%ptr(j)-1)*nloc + (/(k,k=1,nloc)/)

         end do
      end do

      call create_global_to_universal_numbering(new_halo)
      call create_ownership(new_halo)
      
    end subroutine make_halo_dg


  end subroutine make_global_numbering_DG
  

  subroutine make_global_numbering &
       (new_nonods, new_ndglno, Nonods, Totele, NDGLNO, element, halos,&
       & element_halo, new_halos) 
    ! Construct the global node numbering based on the node numbering for
    ! linear tets given in NDGLNO. 
    integer, intent(in) :: nonods, totele
    integer, intent(in), dimension(:), target :: ndglno
    type(element_type), intent(in) :: element
    !! The level 1 and 2 halos associated with the incoming mesh.
    type(halo_type), dimension(:), intent(in), optional :: halos
    !! The full element halo associated with these meshes.
    type(halo_type), intent(in), optional :: element_halo
    !! The level 1 and 2 halos associated with the new mesh.
    type(halo_type), dimension(:), intent(out), optional :: new_halos

    integer, dimension(:), intent(out) :: new_ndglno
    integer, intent(out) :: new_nonods

    ! Adjacency lists.
    type(csr_sparsity) :: NEList, NNList, EEList

    logical :: D3
    integer :: dim, vertex_count, facet_vertex_count, interior_facets, exterior_facets, owned_nodes

    ! Number of nodes associated with each object.
    integer :: face_len, edge_len

    ! Total nodes associated with an element
    integer :: element_tot_len

    integer :: ele, ele2, node, node2, new_node, new_node2, i, j, k, halo

    integer, dimension(:), allocatable :: n

    ! Process number of this processor
    integer :: rank

    ! Owner and halo_level of all the nodes
    integer, dimension(:), allocatable :: node_owner, receive_halo_level, &
         & new_receive_halo_level
    integer, dimension(:,:), allocatable :: new_node_owner
    integer :: this_node_owner, this_receive_halo_level

    ! Cache for positions in ndglno which are currently being worked on and
    ! new values to go in them
    integer, dimension(:), allocatable :: ndglno_pos, ndglno_val, face_nodes

    ! Map from old node numbers to new ones.
    integer, dimension(:), allocatable :: node_map
    
    ! In each case, halos is the last dimension of halos.
    type(integer_set), dimension(:), allocatable :: this_send_targets, node_targets,&
         & node2_targets
    type(integer_set), dimension(:,:), allocatable :: old_send_targets, new_send_targets

    ! Nodes in current element
    integer, pointer, dimension(:) :: ele_node

    ! Nodes and elements adjacent to current node
    integer, pointer, dimension(:) :: adjacent_nodes, adjacent_elements

    ! Flag for whether halos are being calculated:
    logical :: have_halos

    ! List of vertex numbers in an element.
    integer, dimension(:), pointer :: ele1_vertices, ele2_vertices

!DEBUG
    character(len=10) :: format

    ! Ascertain whether we are calculating halos or not.
    if (present(halos)) then
       assert(present(new_halos))
       assert(present(element_halo))
       allocate(this_send_targets(size(halos)), node_targets(size(halos)),&
         & node2_targets(size(halos)))

       have_halos=.true.
    else
       allocate(this_send_targets(0))
       have_halos=.false.
    end if

    ! Dimensionality flag.
    dim = element%numbering%dimension
    D3=dim==3
    
    ! Vertices per element.
    vertex_count=element%numbering%vertices
    ! Vertices per surface element
    facet_vertex_count = vertex_count - 1

    call MakeLists_Dynamic(Nonods, Totele, Vertex_count, ndglno, D3, NEList,&
       & NNList, EEList)

    new_ndglno=0

    rank=getprocno()
    call halo_lists(halos, node_owner, receive_halo_level,&
         & old_send_targets) 

    !----------------------------------------------------------------------
    ! Calculate the total number of nodes by adding the vertices, edge
    ! elements, face elements and interior elements.
    !----------------------------------------------------------------------

    
    ! Add any interior elements
    new_nonods=totele*element%numbering%nodes_per(dim)

    if (dim>0) then
       ! For dim 0 there are no facets.

       ! Interior faces
       interior_facets=0.5*count(EEList%colm/=0)

       exterior_facets=size(EEList%colm)-2*interior_facets

       new_nonods=new_nonods+element%numbering%nodes_per(dim-1)&
            &*(interior_facets+exterior_facets)

       if (dim>1) then
          ! Vertices.
          new_nonods=new_nonods+element%numbering%nodes_per(0)&
               *nonods
          
          if (dim>2) then
             ! For dim==3 only, still need edges.

             ! Edges. 
             ! 0.5*size(NNList%colm) is the number of edges in the mesh.
             new_nonods=new_nonods+0.5*size(NNList%colm)&
                  &*element%numbering%nodes_per(1)
             
          end if
       end if

    end if


    ! Total nodes per element
    element_tot_len=element%numbering%nodes

    if (have_halos) then
       allocate(new_node_owner(new_nonods, 0:size(halos)), &
            &   new_receive_halo_level(size(new_ndglno)), &
            &   new_send_targets(new_nonods, size(halos)))
       do i = 1, size(new_send_targets,1)
          do j = 1,  size(new_send_targets,2)
             call allocate(new_send_targets(i,j))
          end do
       end do
    end if

    ! n is the current maximum node number.

    ! We need one n plus one for each halo so as to separately number the
    ! nodes in each halo. We then transplant them into one long list
    ! afterwards.
    if(have_halos) then
       allocate(n(0:size(halos)))
    else
       allocate(n(0:0))
    end if
    n=0

    !----------------------------------------------------------------------
    ! Vertex numbers
    !----------------------------------------------------------------------
    
    ! ndglno_pos stores the postion of the vertices of the current element
    ! in the new_ndglno list. ndglno_val stores the corresponding values.
    allocate(ndglno_pos(vertex_count), ndglno_val(vertex_count))
    allocate(node_map(nonods))
    node_map=0
    
    do i=1,totele
       ele_node=>NDGLNO(vertex_count*(i-1)+1:vertex_count*i)
       ndglno_pos=(i-1)*element_tot_len &
            + vertex_num(ele_node, &
                         ele_node, &
                         element%numbering)
       

       ! Pick up those values which have been done already.
       ndglno_val=node_map(ele_node)

       do j=1,vertex_count
          if (ndglno_val(j)==0) then
             n(receive_halo_level(ele_node(j)))=&
                  n(receive_halo_level(ele_node(j)))+1 
             ndglno_val(j)=n(receive_halo_level(ele_node(j)))

             node_map(ele_node(j))=ndglno_val(j)
             
             if(receive_halo_level(ele_node(j)) == 0 .and. have_halos) then
               call copy(new_send_targets(n(0),:), &
                    &    old_send_targets(ele_node(j),:))
             end if
          end if
       end do
       
       new_ndglno(ndglno_pos) = ndglno_val
       if (have_halos) then
          do j=1,vertex_count
             new_node_owner(ndglno_val(j),receive_halo_level(ele_node(j)))&
                  =node_owner(ele_node(j))
          end do
          new_receive_halo_level(ndglno_pos)=receive_halo_level(ele_node)
       end if
    end do

    deallocate(ndglno_pos, ndglno_val)

    !----------------------------------------------------------------------
    ! Edge numbers.
    !----------------------------------------------------------------------
    ! edge_len is the number of non-vertex nodes on each edge.
    edge_len=element%numbering%nodes_per(1)

    if (edge_len>0) then
       allocate(ndglno_pos(edge_len))

       do node=1,size(NNList,1)
          adjacent_nodes=>row_m_ptr(NNlist, node)

          do j=1, size(adjacent_nodes)
             node2=adjacent_nodes(j)
             
             new_node=node_map(node)
             new_node2=node_map(node2)

             ! Listen very carefully, I shall do each edge only once!
             if (node2<=node) cycle
             
             call conduct_halo_voting((/node,node2/), this_node_owner,&
                  & this_receive_halo_level, this_send_targets)
             
             ! Now double loop over all adjacent elements and number in the
             ! elements which border this edge.
             adjacent_elements=>row_m_ptr(NEList, node2)
             
             do k=1, size(adjacent_elements)
                ele=adjacent_elements(k)

                if (any(ele==row_m(NEList, node))) then
                   ! This element contains this edge.

                   ! Nodes in this element

                   ! This horrible mess finds the appropriate nodes in this
                   ! element and assigns the next edge_len indices to them.
                   ndglno_pos=(ele-1)*element_tot_len &
                        +edge_num((/node,node2/), &
                        &  NDGLNO(vertex_count*(ele-1)+1:vertex_count*ele), &
                        &  element%numbering, interior=.true.)

                   new_ndglno(ndglno_pos) &
                        = sequence(n(this_receive_halo_level)+1, edge_len)

                   if (have_halos) then
                      new_receive_halo_level(ndglno_pos)= &
                           this_receive_halo_level
                      new_node_owner(new_ndglno(ndglno_pos),&
                           & this_receive_halo_level)=this_node_owner
                   end if
                      
                end if
             end do

             ! Because we by definition own send nodes, we can set these
             ! directly rather than via the node element list.
             if (any(key_count(this_send_targets)/=0)) then
                do i=n(this_receive_halo_level)+1,&
                     n(this_receive_halo_level)+edge_len

                   call copy(new_send_targets(i,:), this_send_targets)
                end do
             end if
             
             ! Move on the node count.
             n(this_receive_halo_level)=n(this_receive_halo_level)+edge_len
             
             if (have_halos) then
                ! Clean up send targets
                call deallocate(this_send_targets)
             end if
          end do

       end do

       deallocate(ndglno_pos)

    end if

    !----------------------------------------------------------------------
    ! Interior face numbers - only for 3D.
    ! 
    ! 2D interior face number are just interior numbers and are dealt with below.
    !----------------------------------------------------------------------
    face_len=element%numbering%nodes_per(2)

    if (D3.and.face_len>0) then
       allocate(face_nodes(facet_vertex_count), ndglno_pos(face_len))

       do ele=1,size(EEList,1)
          adjacent_elements=>row_m_ptr(EEList, ele)
          

          do j=1, size(adjacent_elements)
             ele2=adjacent_elements(j)

             ! Skip exterior faces.
             if (ele2==0) cycle

             ! Listen very carefully, I shall do each face only once!
             if (ele2<=ele) cycle

             ele1_vertices=>NDGLNO(vertex_count*(ele-1)+1:vertex_count*ele)
             ele2_vertices=>NDGLNO(vertex_count*(ele2-1)+1:vertex_count*ele2)

             ! Find the vertices on the common face.
             face_nodes=common(ele1_vertices,ele2_vertices)

             ! Work out who should own these nodes.
             call conduct_halo_voting(face_nodes, this_node_owner,&
                  & this_receive_halo_level, this_send_targets)

             ! Set new_ndglno for this face in ele1.
             ndglno_pos=(ele-1)*element_tot_len&
                  +face_num(face_nodes, ele1_vertices, &
                  element%numbering, interior=.true.)
             
             new_ndglno(ndglno_pos) &
                  = sequence(n(this_receive_halo_level)+1, face_len)

             if(have_halos) then
                new_receive_halo_level(ndglno_pos)=this_receive_halo_level
                new_node_owner(ndglno_pos, this_receive_halo_level)&
                     =this_node_owner
             end if

             ! Set new_ndglno for this face in ele2.
             ndglno_pos=(ele2-1)*element_tot_len&
                  +face_num(face_nodes, ele2_vertices, &
                  element%numbering, interior=.true.)

             new_ndglno(ndglno_pos) &
                  = sequence(n(this_receive_halo_level)+1, face_len)

             if (have_halos) then
                new_receive_halo_level(ndglno_pos)=this_receive_halo_level
                new_node_owner(ndglno_pos, this_receive_halo_level)&
                     =this_node_owner
             end if

             ! Because we by definition own send nodes, we can set these
             ! directly rather than via the node element list.
             if (any(key_count(this_send_targets)/=0)) then
                do i=n(this_receive_halo_level)+1,&
                     n(this_receive_halo_level)+face_len

                   call copy(new_send_targets(i,:), this_send_targets)
                end do
             end if
             
             ! Move on the node count.
             n(this_receive_halo_level)=n(this_receive_halo_level)+face_len

             if (have_halos) then
                ! Clean up send targets
                call deallocate(this_send_targets)
             end if

          end do
       end do

       deallocate(face_nodes, ndglno_pos)
    end if

    !----------------------------------------------------------------------
    ! Remaining numbers.
    !----------------------------------------------------------------------

    allocate (ndglno_pos(element%loc))
    this_receive_halo_level = 0

    ! This is the internal nodes of the elements plus the external faces.
    do ele=1,size(EEList,1)
              
       ndglno_pos=sequence((ele-1)*element%loc+1, element%loc)

       ! Work out who should own these nodes.
       ele1_vertices=>NDGLNO(vertex_count*(ele-1)+1:vertex_count*ele)
       call conduct_halo_voting(ele1_vertices, this_node_owner,&
            & this_receive_halo_level, this_send_targets)
              
       do i=1, element%loc
          if (new_ndglno(ndglno_pos(i))==0) then
             n(this_receive_halo_level)=n(this_receive_halo_level)+1
             new_ndglno(ndglno_pos(i))=n(this_receive_halo_level)

             if(have_halos) then
                new_receive_halo_level(ndglno_pos(i))=this_receive_halo_level
                new_node_owner(ndglno_pos(i), this_receive_halo_level)&
                     =this_node_owner
             end if

             if (any(key_count(this_send_targets)/=0)) then
                call copy(new_send_targets(n(this_receive_halo_level),:), &
                     &    this_send_targets)
                
             end if
          end if
       end do

       if (have_halos) then
          ! Clean up send targets
          call deallocate(this_send_targets)
       end if
       
    end do

    deallocate(ndglno_pos)
    call deallocate(NEList)
    call deallocate(NNList)
    call deallocate(EEList)

    write(format,'(a,i0,a)') '(',element%loc,'i4)'

    ASSERT(sum(n)==new_nonods)

    owned_nodes=n(0)
    do i=size(n)-1,0,-1
       ! Work out the offset for halo nodes.
       n(i)=sum(n(0:i-1))
    end do
    n(0)=0

    if (have_halos) then
       new_ndglno=new_ndglno+n(new_receive_halo_level)

       ! Repack new_node_owner into a single list.
       do halo=1,size(halos)-1
          new_node_owner(n(halo)+1:n(halo+1),0)&
               =new_node_owner(:n(halo+1)-n(halo),halo)
       end do
       halo=size(halos)
       new_node_owner(n(halo)+1:,0)&
            =new_node_owner(:new_nonods-n(halo),halo)

!!$       call remove_spurious_sends(new_send_targets, new_ndglno,&
!!$            & element_tot_len, element_halo)

       call generate_new_halos(new_halos, new_ndglno, new_node_owner(:,0)&
            &, new_receive_halo_level, totele,&
            & element%loc, owned_nodes, new_send_targets) 

       do i = 1, size(old_send_targets,1)
          do j = 1, size(old_send_targets,2)
             call deallocate(old_send_targets(i, j))
          end do
       end do
       do i = 1, size(new_send_targets,1)
          do j = 1, size(old_send_targets,2)
             call deallocate(new_send_targets(i, j))
          end do
       end do
       deallocate(old_send_targets)
       deallocate(new_send_targets)
    end if

  contains

    subroutine conduct_halo_voting(vertices, this_node_owner,&
         & this_receive_halo_level, this_send_targets) 
      !!< Given a list of vertices of a topological object (element, face,
      !!< edge etc.), apply the voting algorithm to determine which
      !!< processor owns that object and therefore what the halo properties
      !!< of any new nodes placed on that object should be.
      integer, dimension(:), intent(in) :: vertices
      integer, intent(out) :: this_node_owner
      integer, intent(out) :: this_receive_halo_level
      type(integer_set), dimension(:), intent(out) :: this_send_targets

      integer :: halo, i

      ! Defaults which may be overwritten in the following if block
      this_receive_halo_level=0                
      this_node_owner=maxval(node_owner(vertices))

      ! Ensure the send targets pointer is null
      call allocate(this_send_targets)
             
      if (.not.have_halos) return
      
      if (any(node_owner(vertices)/=this_node_owner)) then
         ! Contested object. Nodes have the same halo properties as
         ! the winning side.
       
         if (this_node_owner/=rank) then
            ! We don't own the new nodes so they go on the receive list.
            this_receive_halo_level=maxval(receive_halo_level(vertices))

         else
            ! We do own the new nodes so we form send lists for them.
            ! 
            ! The send sets are the intersection of the send
            ! sets. This is because if one vertex is not a sender at the
            ! given halo level, the obect does not belong to that halo.
            do halo = 1, size(halos) 
               call set_intersection(this_send_targets(halo),&
                    old_send_targets(vertices, halo), &
                    mask=(node_owner(vertices)==rank))
            end do
         end if

         !! If we get here, the object is uncontested. It's either all ours
         !! or all someone else's.
      else if(all(receive_halo_level(vertices)>0)) then
         ! It's someone else's.
         
         this_receive_halo_level=maxval(receive_halo_level(vertices))

      else if (all(key_count(&
           pack(old_send_targets(vertices,:),mask=.true.))>0)) then
         ! It's ours and there's something to send.
         
         do halo = 1, size(halos) 
            call set_intersection(this_send_targets(halo),&
                 old_send_targets(vertices, halo), &
                 mask=(node_owner(vertices)==rank))
         end do
         
      end if

    end subroutine conduct_halo_voting

!!$    subroutine remove_spurious_sends(send_targets, ndglno, nloc, element_halo)
!!$      !!< Given the node ownership and a complete element halo, remove any
!!$      !!< sends which apply only to elements about which the receiving
!!$      !!< processor is unaware.
!!$
!!$      !! Send_targets provides the list of processors to which each
!!$      !! node is broadcast at each halo level. It is nonods x halos
!!$      type(ilist), dimension(:,:), intent(inout) :: send_targets
!!$      !! The element node list and number of nodes per element for the NEW
!!$      !! mesh.
!!$      integer, dimension(:), intent(in), target :: ndglno
!!$      integer, intent(in) :: nloc
!!$      type(halo_type), intent(in) :: element_halo
!!$      
!!$      !! For each node, a list of processors which can see that node.
!!$      type(ilist), dimension(:), allocatable :: visible_to
!!$      !! For each processor, a list of halo elements that processor can see.
!!$      type(integer_vector), dimension(:), allocatable :: known_element_lists
!!$
!!$      integer :: e, ele, element_count
!!$      integer :: n, node_count, proc, h, halo_count
!!$      integer, dimension(:), pointer :: this_ele
!!$      type(ilist) :: tmplist
!!$
!!$      element_count=size(ndglno)/nloc
!!$      node_count=size(send_targets,1)
!!$      halo_count=size(send_targets,2)
!!$
!!$      allocate(visible_to(node_count))
!!$      allocate(known_element_lists(halo_proc_count(element_halo)))
!!$
!!$      ! Check that nloc does actually divide size(ndglno)
!!$      assert(size(ndglno)==element_count*nloc)
!!$
!!$      ! Retrieve the list of elements which each processor can see.
!!$      call element_halo_communicate_visibility(element_halo,&
!!$           & known_element_lists)
!!$
!!$      ! Mark all nodes in elements visible to a processor as themselves
!!$      ! visible to that processor
!!$      do proc=1,halo_proc_count(element_halo)
!!$         do e = 1, size(known_element_lists(proc)%ptr)
!!$            ele = known_element_lists(proc)%ptr(e)
!!$            this_ele=>ndglno((ele-1)*nloc+1:ele*nloc)
!!$
!!$            do n=1,size(this_ele)
!!$               
!!$               call insert_ascending(visible_to(this_ele(n)), proc)
!!$            
!!$            end do
!!$         end do
!!$      end do
!!$
!!$      do proc=1,halo_proc_count(element_halo)
!!$         deallocate(known_element_lists(proc)%ptr)
!!$      end do
!!$
!!$      do n=1,node_count
!!$         
!!$         do h=1,halo_count
!!$            ! Do nothing for non-send nodes.
!!$            if (send_targets(n,h)%length==0) cycle
!!$         
!!$            
!!$            tmplist=intersect_ascending(send_targets(n,h),visible_to(n))
!!$         
!!$            call flush_list(send_targets(n,h))
!!$            send_targets(n,h)=tmplist
!!$         end do
!!$
!!$         call flush_list(visible_to(n))
!!$      end do
!!$      deallocate(visible_to)
!!$
!!$
!!$    end subroutine remove_spurious_sends
    
    subroutine generate_new_halos(new_halos, new_ndglno, new_node_owner&
            &, new_receive_halo_level, elements,&
            & nloc, owned_nodes, send_targets)
      !!< Given the node ownership and halo level information on the new
      !!< mesh, construct halos for the new mesh.
      !!< Note that these halos are unsorted as sorting the halos requires
      !!< the coordinate field which is not available here.
      type(halo_type), dimension(:), intent(out) :: new_halos
      integer, dimension(:), intent(in), target :: new_ndglno
      integer, dimension(:), intent(in), target :: new_node_owner
      integer, dimension(:), intent(in), target :: new_receive_halo_level
      integer, intent(in) :: elements, nloc, owned_nodes
      ! Send_targets provides the list of processors to which each
      ! node is broadcast at each halo level. It is nonods x halos
      type(integer_set), dimension(:,:), intent(in) :: send_targets

      integer :: processors, this_proc, halo, n, proc
      
      type(integer_set), dimension(:), allocatable :: sends, receives
      integer, dimension(:), pointer :: ele_nodes

      processors=getnprocs()
      this_proc=getprocno()

      allocate(sends(processors))
      allocate(receives(processors))
      do n=1,size(sends)
         call allocate(sends(n))
         call allocate(receives(n))
      end do
      
      halo_loop: do halo = 1, size(new_halos)
         
         do n=1, size(new_receive_halo_level)
            ! Receive node
            if (new_receive_halo_level(n)>0.and.new_receive_halo_level(n)<=halo) then
               call insert(receives(new_node_owner(new_ndglno(n))), new_ndglno(n))
            end if
         end do            

         do n = 1, new_nonods
            ! Send node
            do i = 1, key_count(send_targets(n,halo))
               call insert(sends(fetch(send_targets(n,halo),i)),n)
            end do
         end do
         
         call allocate(new_halos(halo), nsends=key_count(sends), &
              nreceives=key_count(receives), nprocs=processors, &
              nowned_nodes=owned_nodes)
         
         do proc = 1, processors
            
            new_halos(halo)%sends(proc)%ptr=set2vector(sends(proc))
            call deallocate(sends(proc))
            
            new_halos(halo)%receives(proc)%ptr=set2vector(receives(proc))
            call deallocate(receives(proc))
            
         end do
         
         call print_halo(new_halos(halo), 0)
    
         assert(trailing_receives_consistent(new_halos(halo)))
         assert(halo_valid_for_communication(new_halos(halo)))

      end do halo_loop
      
    end subroutine generate_new_halos

    function sequence(start, len)
      ! Return len consecutive integers starting at start.
      integer, intent(in) :: start, len
      integer, dimension(len) :: sequence
      
      integer :: i

      forall (i=1:len)
         sequence(i)=start+i-1
      end forall

    end function sequence

    function common(list1, list2) result (list)
      ! Return the common members of list1 and list2.
      integer, dimension(:), intent(in) :: list1, list2
      integer, dimension(count(&
           (spread(list1,2,size(list2))-spread(list2,1,size(list1)))==0)) &
           :: list

      integer :: i, j

      j=0

      do i=1,size(list2)
         if (any(list1(i)==list2)) then
            j=j+1
            list(j)=list1(i)
         end if
      end do

      ASSERT(j==size(list))

    end function common
    
    subroutine halo_lists(halos, node_owner, receive_halo_level, old_send_targets)
      type(halo_type), dimension(:), intent(in), optional :: halos
      integer, dimension(:), allocatable, intent(out) :: node_owner
      integer, dimension(:), allocatable, intent(out) :: receive_halo_level
      !! Targets to broadcast each node to. Nonods x halos
      type(integer_set), dimension(:,:), allocatable, intent(out) :: old_send_targets

      integer :: h, n, p

      allocate(node_owner(nonods))
      allocate(receive_halo_level(nonods))

      node_owner = rank
      receive_halo_level = 0
      
      if(.not. present(halos)) return
      allocate(old_send_targets(nonods, size(halos)))
      do h=1,size(old_send_targets,1)
         do p=1,size(old_send_targets,2)
            call allocate(old_send_targets(h,p))
         end do
      end do

      ! Count down through the halos as halo 1 is a subset of halo 2 and
      ! we therefore need halo 1 to overwrite.
      halo_loop: do h = size(halos), 1, -1
         if(.not. associated(halos(h)%receives)) cycle halo_loop
         assert(associated(halos(h)%sends))
         
         proc_loop: do p = 1, halo_proc_count(halos(h))

#ifdef DDEBUG
            if(h < size(halos)) then
              ! Check that this halo really is a subset of all higher level
              ! halos
              assert(all(receive_halo_level(halo_receives(halos(h), p)) /= 0))
            end if
#endif
            receive_halo_level(halo_receives(halos(h), p)) = h

            do n = 1, halo_send_count(halos(h), p)
              call insert(old_send_targets(halo_send(halos(h), p, n),h), p)
            end do
            
            if(h == size(halos)) then
              ! If we're on the highest level halo, set the node owners
              node_owner(halo_receives(halos(h), p)) = p
#ifdef DDEBUG
            else
              ! Otherwise, these should already have been set by the higher
              ! level halo
              assert(all(node_owner(halo_receives(halos(h), p)) == p))
#endif
            end if
            
         end do proc_loop

      end do halo_loop

    end subroutine halo_lists

  end subroutine make_global_numbering

  subroutine make_boundary_numbering(boundary_list, boundary_n_lno, &
       boundary_m_lno, EEList, xnonod&
       &, xndglno, ele_n, boundary_n, ele_m, boundary_m)
    ! Generate boundary numberings to facilitate the evaluation of boundary
    ! integrals. As is usual in fluidity, the suffix n refers to velocity
    ! elements while the suffix m refers to pressure elements.
    !
    ! If f=boundary_list(i,j) then boundary_n_lno((f-1)*boundary%loc+1:f*boundary%loc) is the
    ! vector of local indices of the nodes on the boundary between i and j. 
    ! 
    ! Let f2=boundary_list(j,i). Then:
    ! xndglno(boundary_list((f-1)*boundary%loc+1:f*boundary%loc)) = 
    !          xndglno(boundary_list((f2-1)*boundary%loc+1:f2*boundary%loc))
    type(csr_sparsity), intent(in) :: EEList
    type(element_type), intent(in) :: ele_n, boundary_n
    type(element_type), optional, intent(in) :: ele_m, boundary_m
    integer, intent(in) :: xnonod
    integer, dimension(size(EEList,1)*ele_n%loc), intent(in), target ::&
         & xndglno 

    type(csr_matrix), intent(out) :: boundary_list
    integer, dimension(boundary_n%loc*entries(EEList)), intent(out), optional :: & 
         boundary_n_lno
    integer, dimension(:), optional, intent(out) :: & 
         boundary_m_lno

    integer, dimension(ele_n%numbering%boundaries) :: neigh
    integer, dimension(ele_n%numbering%vertices) :: vertices
    integer, dimension(:), pointer :: ele_i, ele_j
    
    integer :: boundary_cnt, i, j, m, n, p, rlen
    logical :: logtest
    integer, dimension(boundary_n%numbering%vertices) :: boundary_i, boundary_j

    ewrite(2,*) "subroutine make_boundary_numbering"

    !CHECK(size(neigh))

    call allocate(boundary_list, EEList, type=CSR_INTEGER)
    if (present(boundary_m_lno)) boundary_m_lno=0

    logtest = associated(boundary_list%val)
    !CHECK(logtest)
    logtest = associated(boundary_list%ival)
    !CHECK(logtest)

    if(present(boundary_m_lno)) then
!       MSG("Checking size of boundary_m_lno")
!       ASSERT(size(boundary_m_lno)==(boundary_m%loc*entries(EEList)))
    end if

    ewrite(2,*) "zeroing boundary_list"
    call zero(boundary_list)

    if (present(boundary_n_lno)) boundary_n_lno=0
    if (present(boundary_m_lno)) boundary_m_lno=0

    boundary_cnt=0

    vertices=local_vertices(ele_n%numbering)
       
    do i=1,size(boundary_list,1)
       ele_i=>xndglno((i-1)*ele_n%loc+1:i*ele_n%loc)
       
       ! Sanity check EEList form.
       rlen = row_length(boundary_List,i)
       ASSERT(rlen==size(neigh))

       neigh=row_m(EEList,i)

       neighbourloop: do j=1, size(neigh)
          ! Exclude boundary boundaries.
          if (neigh(j)==0) then
             cycle neighbourloop 
          end if
          
          ! EEList includes the current ele_n, which is obviously wrong.
          if (neigh(j)==i) cycle neighbourloop

          ! Check to see if the boundary has been done.
          if (val(boundary_list,i,neigh(j))/=0) cycle neighbourloop 

          ele_j=>xndglno((neigh(j)-1)*ele_n%loc+1:neigh(j)*ele_n%loc)

          p=0

          ! Look for common boundaries.
          do m=1,size(vertices)
             do n=1,size(vertices)
                if (ele_i(vertices(m))==ele_j(vertices(n))) then
                   p=p+1
                   boundary_i(p)=m
                   boundary_j(p)=n
                end if
             end do
          end do

          ! Check that we really have found two boundaries.
          !ASSERT(p==boundary_n%numbering%vertices)

          ! Put the boundaries we have found in the next two spots.
          boundary_cnt=boundary_cnt+1

          ! Velocity element boundaries.
          if(present(boundary_n_lno)) then
             boundary_n_lno((boundary_cnt-1)*boundary_n%loc+1: &
                  boundary_cnt*boundary_n%loc)= &
                  & boundary_local_num(boundary_i, ele_n%numbering)
          end if

          ! Pressure element boundaries.
          if(present(boundary_m_lno)) then
             ASSERT(present(boundary_m))
             ASSERT(present(boundary_m))
             ASSERT(present(ele_m))
             boundary_m_lno((boundary_cnt-1)*boundary_m%loc+1: &
                  boundary_cnt*boundary_m%loc)= &
                  boundary_local_num(boundary_i, ele_m%numbering)
          end if

          call set(boundary_list,i,neigh(j),boundary_cnt)

          boundary_cnt=boundary_cnt+1
          
          ! Velocity element boundaries.
          if(present(boundary_n_lno)) then
             boundary_n_lno((boundary_cnt-1)*boundary_n%loc+1: &
                  boundary_cnt*boundary_n%loc)= &
                  boundary_local_num(boundary_j, ele_n%numbering)
          end if

          ! Pressure element boundaries.
          if(present(boundary_m_lno)) then
             boundary_m_lno((boundary_cnt-1)*boundary_m%loc+1: &
                  boundary_cnt*boundary_m%loc)= &
                  boundary_local_num(boundary_j, ele_m%numbering)
          end if

          call set(boundary_list,neigh(j),i,boundary_cnt)
       end do neighbourloop

    end do

    logtest = (boundary_cnt==entries(EEList))
    ASSERT(logtest)

    ewrite(2,*) "END subroutine make_boundary_numbering"

  end subroutine make_boundary_numbering


  subroutine element_halo_communicate_visibility(element_halo,&
       & known_element_lists)
    !!< Given the element halo, for each processor we know about, send
    !!< the list of halo elements we know about.
    type(halo_type), intent(in) :: element_halo
    type(integer_vector), dimension(halo_proc_count(element_halo)),&
         & intent(out) :: known_element_lists

#ifdef HAVE_MPI
    ! For each element, the list of other processors it is sent to.
    type(ilist), dimension(:), allocatable :: visible_list
    ! For each processor, the list of extra elements it can see.
    type(ilist), dimension(:), allocatable :: visible_elements
    type(inode), pointer :: this_item
    
    ! For each processor, the visible_list which is to be sent to it.
    type(integer_vector), dimension(halo_proc_count(element_halo)) ::&
         & send_lists
    integer, dimension(:), allocatable :: requests, receive_list
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer :: proc, e, ierr, rank, communicator, count, nprocs, sends,&
         & receives, pos, p, sendproc, tag
    

    nprocs=halo_proc_count(element_halo)

    allocate(visible_list(node_count(element_halo)))

    do proc=1, nprocs
       do e=1, halo_send_count(element_halo, proc)

          call insert_ascending(&
               visible_list(halo_send(element_halo, proc, e)), &
               proc)
       end do
    end do

    do proc=1, nprocs
       allocate(send_lists(proc)%ptr(&
            sum(visible_list(halo_sends(element_halo,proc))%length)))
    
       ! Use the first halo_send_count places to indicate how many
       ! processors (other than the receiver) know about each element.
       send_lists(proc)%ptr(:halo_send_count(element_halo, proc)) &
            = visible_list(halo_sends(element_halo,proc))%length - 1 
       
       sends = halo_send_count(element_halo, proc)       
       pos=sends
       do e=1, sends
          this_item=>visible_list(halo_send(element_halo, proc, e))&
               &%firstnode 
          do while(associated(this_item))
             ! Eliminate self-cites.
             if (this_item%value/=proc) then
                pos=pos+1
                send_lists(proc)%ptr(pos)=this_item%value
             end if
             
             this_item=>this_item%next
          end do
          
       end do
       assert(pos==size(send_lists(proc)%ptr))
    end do
    
    call flush_lists(visible_list)
    deallocate(visible_list)

    ! Set up non-blocking communications
    communicator = halo_communicator(element_halo)
    allocate(requests(nprocs))
    rank = getrank(communicator)  
    tag = next_mpi_tag()

    do proc=1, nprocs

       call mpi_isend(send_lists(proc)%ptr, &
            size(send_lists(proc)%ptr), MPI_INTEGER,&
            proc-1, tag, communicator, requests(proc), ierr)
       assert(ierr == MPI_SUCCESS)

    end do

    allocate(visible_elements(nprocs))

    ! Wait for incoming data.
    do proc=1, nprocs

       call mpi_probe(MPI_ANY_SOURCE, tag, communicator, status,&
            & ierr)
       assert(ierr == MPI_SUCCESS)

       call mpi_get_count(status, MPI_INTEGER, count, ierr)
       assert(ierr == MPI_SUCCESS)

       allocate(receive_list(count))

       call mpi_recv(receive_list, count, MPI_INTEGER, status(MPI_SOURCE),&
            tag, communicator, MPI_STATUS_IGNORE, ierr)
       assert(ierr == MPI_SUCCESS)

       sendproc=status(MPI_SOURCE)+1
       ! For each element on the recieves list for this processor, record
       ! which other elements know about it.
       receives = halo_receive_count(element_halo, sendproc)
       pos = receives
       do e=1, receives
          ! Recall that the first halo_receive_count entries in the list
          ! simply say how many extra processors know about this element.
          do p=1,receive_list(e)
             pos=pos+1
             call insert_ascending(visible_elements(receive_list(pos)), &
                  halo_receive(element_halo, sendproc, e))
          end do
       end do
       
       assert(pos==size(receive_list))
       deallocate(receive_list)
    end do
    
    ! Wait for sends to complete
    call mpi_waitall(size(requests), requests, MPI_STATUSES_IGNORE, ierr)

    ! Now actually sort out the lists. The elements which proc and this
    ! processor can both see comprise the send and receive elements between
    ! these two processors, plus the elements which other processors send
    ! to proc and which we were told about above.
    do proc=1, nprocs
       sends=halo_send_count(element_halo,proc)
       receives=halo_receive_count(element_halo,proc)
       allocate(known_element_lists(proc)%ptr(&
            sends+receives+visible_elements(proc)%length))

       known_element_lists(proc)%ptr(:sends) &
            = halo_sends(element_halo, proc)

       known_element_lists(proc)%ptr(sends+1:sends+receives) &
            = halo_receives(element_halo, proc)
       
       known_element_lists(proc)%ptr(sends+receives+1:) &
            = list2vector(visible_elements(proc))
    end do

    call flush_lists(visible_elements)
#else
    FLAbort("Communicating halo visibility makes no sense without MPI.")
#endif

  end subroutine element_halo_communicate_visibility

end module global_numbering
