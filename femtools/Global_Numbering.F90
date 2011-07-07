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
  use memory_diagnostics
  
  implicit none

  private
  
  public :: make_global_numbering_DG, make_global_numbering, &
       & make_global_numbering_trace, make_global_numbering_new


contains

  subroutine make_global_numbering_new(mesh)
    ! Construct a new global numbering for mesh using the ordering from
    !  topology.
    type(mesh_type), intent(inout), target :: mesh

    type(csr_sparsity) :: facet_list, edge_list
    type(csr_matrix) :: facet_numbers, edge_numbers
    integer, dimension(0:mesh_dim(mesh)) :: entity_counts, dofs_per
    type(cell_type), pointer :: cell
    type(element_type), pointer :: element, topo_element
    integer :: d, e, ele, entity
    integer, dimension(:), pointer :: ele_dofs,topo_dofs,facets
    integer, dimension(2) :: edge

    element=>mesh%shape
    cell=>element%cell    
    topo_element=>mesh%topology%shape
    
    call makelists(mesh%topology, NNlist=edge_list)
    facet_list=mesh%topology%faces%face_list%sparsity

    ! Number the topological entities.
    call number_topology
    call calculate_dofs_per_entity

    do ele=1,element_count(mesh)
       ele_dofs=>ele_nodes(mesh, ele)
       topo_dofs=>ele_nodes(mesh%topology, ele)
       facets=>row_m_ptr(facet_list,ele)
#ifdef DDEBUG
       ele_dofs=0
#endif

       do d=0,mesh_dim(mesh)
          do e=1,cell%entity_counts(d)
             if (d==0) then
                entity=topo_dofs(e)
             else if (d==cell%dimension-1) then
                entity=facets(e)
             else if (d==cell%dimension) then
                entity=ele
             else if (d==1) then
                ! This case only gets hit for the edges of 3d elements.

                ! The edge consists of the global dofs in the topology.
                edge=topo_dofs(entity_vertices(cell,[1,e]))
                ! Now look up the corresponding edge number.
                entity=val(edge_numbers,edge(1),edge(2))
             end if
  
             ele_dofs(element%entity2dofs(d,e)%dofs)=&
                  entity_dofs(d,entity, dofs_per)
          end do
       end do
       
       assert(all(ele_dofs>0))

    end do

    mesh%nodes=sum(entity_counts*dofs_per)

    assert(mesh%nodes==maxval(mesh%ndglno))

    if (mesh_dim(mesh)>1) call deallocate(facet_numbers)
    if (mesh_dim(mesh)>2) call deallocate(edge_numbers)
    call deallocate(edge_list)    
    
  contains
   
    function entity_dofs(dim, entity, dofs_per) result (dofs)
      ! Global dofs associated with local entity
      integer, intent(in) :: dim, entity
      integer, dimension(0:) :: dofs_per
      integer, dimension(dofs_per(dim)) :: dofs
      
      integer :: i,d,e

      dofs=0
      do d=0,dim-1
         dofs=dofs+dofs_per(d)*entity_counts(d)
      end do
      dofs=dofs+(entity-1)*dofs_per(dim)+[(i, i=1,dofs_per(dim))]
      
    end function entity_dofs
    
    subroutine calculate_dofs_per_entity
      integer :: i
      
      do i=0,mesh_dim(mesh)
         if (allocated(element%entity2dofs(i,1)%dofs)) then
            dofs_per(i)=size(element%entity2dofs(i,1)%dofs)
         else
            dofs_per(i)=0
         end if
      end do

    end subroutine calculate_dofs_per_entity

    subroutine number_topology
      integer :: ele1, ele2, node1, node2      
      integer :: facets, edges,n
      integer, dimension(:), pointer :: neigh

      ! No need to number vertices: it comes from the topology.
      entity_counts(0)=node_count(mesh%topology)

      if (mesh_dim(mesh)>0) then
         entity_counts(mesh_dim(mesh))=element_count(mesh)
      end if

      ! Number the facets
      if (mesh_dim(mesh)>1) then
         call allocate(facet_numbers, facet_list, type=CSR_INTEGER)
         call zero(facet_numbers)

         facets=0
         do ele1=1, element_count(mesh)
            neigh=>ele_neigh(mesh%topology,ele1)
            do n=1, size(neigh)
               ele2=neigh(n)
               if (ele2<0) then
                  ! Exterior facet
                  facets=facets+1
                  call set(facet_numbers,ele1,ele2,facets)
               else if (ele2>ele1) then
                  facets=facets+1
                  call set(facet_numbers,ele1,ele2,facets)
                  call set(facet_numbers,ele2,ele1,facets)
               end if
            end do
         end do
         
         entity_counts(mesh_dim(mesh)-1)=facets
      end if

      ! Number the edges
      if (mesh_dim(mesh)>2) then
         ! If mesh_dim(mesh)==1 then this is subsumed in the facets.
         call allocate(edge_numbers, edge_list, type=CSR_INTEGER)
         call zero(edge_numbers)

         edges=0
         do node1=1, node_count(mesh%topology)
            neigh=>row_m_ptr(edge_list, node1)
            do n=1, size(neigh)
               node2=neigh(n)
               if (node2>node1) then
                  edges=edges+1
                  call set(edge_numbers,node1,node2,edges)
                  call set(edge_numbers,node2,node1,edges)
               end if
            end do
         end do
         entity_counts(1)=edges
      end if

    end subroutine number_topology
    

  end subroutine make_global_numbering_new
  
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

    new_nonods=totele*element%ndof

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
      type(halo_type), intent(in) :: element_halo
      type(halo_type), intent(out) :: new_halo
      
      integer, dimension(size(element_halo%sends)) :: nsends
      integer, dimension(size(element_halo%receives)) :: nreceives
      
      integer :: i,j,k, nloc

      nloc=element%ndof

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
              nowned_nodes=element_halo%nowned_nodes*element%ndof, &
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

  subroutine make_global_numbering_trace(mesh)
    ! Construct a global node numbering for a trace mesh
    !
    ! Note that this code is broken for mixed element meshes.
    type(mesh_type), intent(inout) :: mesh
    !
    integer :: ele, totele, ni, ele_2, current_global_index
    integer, pointer, dimension(:) :: neigh
    integer :: face_1,face_2,nfaces,i,face_loc, nloc

    totele = mesh%elements
    face_loc = mesh%faces%shape%ndof
    nloc = mesh%shape%ndof

    !count up how many faces there are
    nfaces = 0
    do ele = 1, totele
       neigh => ele_neigh(mesh,ele)
       do ni = 1, size(neigh)
          ele_2 = neigh(ni)
          if(ele_2<ele) then
             nfaces=nfaces+1
          end if
       end do
    end do
    mesh%nodes = nfaces*face_loc

    !construct mesh%ndglno
    mesh%ndglno = 0
    current_global_index = 0
    do ele = 1, totele
       neigh => ele_neigh(mesh,ele)
       do ni = 1, size(neigh)
          ele_2 = neigh(ni)
          if(ele_2<ele) then
             face_1=ele_face(mesh, ele, ele_2)
             mesh%ndglno((ele-1)*nloc+face_local_nodes(mesh,face_1))&
                  &=current_global_index+(/(i, i=1,face_loc)/)
             if(ele_2>0) then
                !it's not a domain boundary
                !not quite sure how this works in parallel
                face_2=ele_face(mesh, ele_2, ele)
                mesh%ndglno((ele_2-1)*nloc+face_local_nodes(mesh,face_2))&
                     &=current_global_index+(/(i, i=1,face_loc)/)
             end if
             current_global_index = current_global_index + &
                  & mesh%faces%shape%ndof
          end if
       end do
    end do
    if(current_global_index /= mesh%nodes) then
       FLAbort('bad global index count in make_global_numbering_trace')
    end if
    if(any(mesh%ndglno==0)) then
       FLAbort('Failed to fully populate trace mesh ndglno')
    end if

  end subroutine make_global_numbering_trace
  

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

    allocate (ndglno_pos(element%ndof))
    this_receive_halo_level = 0

    ! This is the internal nodes of the elements plus the external faces.
    do ele=1,size(EEList,1)
              
       ndglno_pos=sequence((ele-1)*element%ndof+1, element%ndof)

       ! Work out who should own these nodes.
       ele1_vertices=>NDGLNO(vertex_count*(ele-1)+1:vertex_count*ele)
       call conduct_halo_voting(ele1_vertices, this_node_owner,&
            & this_receive_halo_level, this_send_targets)
              
       do i=1, element%ndof
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

       call generate_new_halos(new_halos, new_ndglno, new_node_owner(:,0)&
            &, new_receive_halo_level, element%ndof, owned_nodes,&
            & new_send_targets) 

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

      integer :: halo

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

    
    subroutine generate_new_halos(new_halos, new_ndglno, new_node_owner&
            &, new_receive_halo_level, nloc, owned_nodes, send_targets)
      !!< Given the node ownership and halo level information on the new
      !!< mesh, construct halos for the new mesh.
      !!< Note that these halos are unsorted as sorting the halos requires
      !!< the coordinate field which is not available here.
      type(halo_type), dimension(:), intent(out) :: new_halos
      integer, dimension(:), intent(in), target :: new_ndglno
      integer, dimension(:), intent(in), target :: new_node_owner
      integer, dimension(:), intent(in), target :: new_receive_halo_level
      integer, intent(in) :: nloc, owned_nodes
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

end module global_numbering
