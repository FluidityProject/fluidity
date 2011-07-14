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
  use quicksort
  
  implicit none

  private
  
  public :: make_global_numbering_new, make_global_numbering_dg


contains

  subroutine make_global_numbering_new(mesh, uid)
    ! Construct a new global numbering for mesh using the ordering from
    !  topology.
    type(mesh_type), intent(inout), target :: mesh
    ! The universal identifier. Required for parallel.
    type(scalar_field), intent(in), optional :: uid
    
    type(csr_sparsity) :: facet_list, edge_list
    type(csr_matrix) :: facet_numbers, edge_numbers
    integer, dimension(0:mesh_dim(mesh)) :: entity_counts, dofs_per
    type(cell_type), pointer :: cell
    type(element_type), pointer :: element, topo_element
    type(mesh_type), pointer :: topology
    integer :: d, e, ele, entity, max_vertices
    integer, dimension(:), pointer :: ele_dofs,topo_dofs, facets
    integer, dimension(2) :: edge
    logical :: have_facets, have_halos

    type(integer_set), dimension(:,:), allocatable :: entity_send_targets
    integer, dimension(:), allocatable :: entity_owner
    integer, dimension(:), allocatable :: entity_receive_level
    
    ! Array to enable all the entities to be sorted together according to
    !  halo levels and universal numbers.
    integer, dimension(:,:), allocatable :: entity_sort_list
    ! Order in which entities should be visited
    integer, dimension(:), allocatable :: visit_order
    ! Dofs associated with each entity.
    integer, dimension(:), allocatable :: entity_dof_starts

    element=>mesh%shape
    cell=>element%cell    
    topo_element=>mesh%topology%shape
    topology=>mesh%topology

    have_facets=associated(topology%faces)
    have_halos=associated(topology%halos)
    ! Halos require uid.
    assert((.not.have_halos).or.(.not.present(uid)))

    call makelists(topology, NNlist=edge_list)
    if (have_facets) then
       facet_list=topology%faces%face_list%sparsity
    end if

    entity_counts=count_topological_entities(topology, edge_list, facet_list)
    max_vertices=size(entity_vertices(cell,[cell%dimension,1]))
    if (have_halos) then
       allocate(entity_sort_list(sum(entity_counts),0:max_vertices))
       entity_sort_list=-1
    end if
    allocate(visit_order(sum(entity_counts)))
    allocate(entity_dof_starts(sum(entity_counts)))

    ! Number the topological entities.
    call number_topology
    call calculate_dofs_per_entity
    if (have_halos) call create_topology_halos(topology%halos)
    call entity_order
    call topology_dofs

    do ele=1,element_count(mesh)
       ele_dofs=>ele_nodes(mesh, ele)
       topo_dofs=>ele_nodes(topology, ele)
       if (have_facets.and.mesh_dim(mesh)>1) facets=>row_ival_ptr(facet_numbers,ele)
#ifdef DDEBUG
       ele_dofs=0
#endif

       do d=0,mesh_dim(mesh)
          do e=1,cell%entity_counts(d)
             if (d==0) then
                entity=topo_dofs(e)
             else if (d==cell%dimension) then
                entity=sum(entity_counts(:d-1))+ele
             else if (d==cell%dimension-1) then
                if (have_facets) entity=facets(e)
             else if (d==1) then
                ! This case only gets hit for the edges of 3d elements.

                ! The edge consists of the global dofs in the topology.
                edge=topo_dofs(entity_vertices(cell,[1,e]))
                ! Now look up the corresponding edge number.
                entity=ival(edge_numbers,edge(1),edge(2))
             end if
  
             ele_dofs(element%entity2dofs(d,e)%dofs)=&
                  entity_dofs(d,entity, dofs_per)
             
          end do
       end do
       
       assert(all(ele_dofs>0))

    end do

    mesh%nodes=sum(entity_counts*dofs_per)

    assert(mesh%nodes==maxval(mesh%ndglno))

    if(have_facets) then
       if (mesh_dim(mesh)>1) call deallocate(facet_numbers)
    end if
    if (mesh_dim(mesh)>2) call deallocate(edge_numbers)
    call deallocate(edge_list)    
    call deallocate_topology_halos

  contains
   
    function entity_dofs(dim, entity, dofs_per) result (dofs)
      ! Global dofs associated with local entity
      integer, intent(in) :: dim, entity
      integer, dimension(0:) :: dofs_per
      integer, dimension(dofs_per(dim)) :: dofs
      
      integer :: i,d,e

      dofs=entity_dof_starts(entity)&
           +[(i, i=1,dofs_per(dim))]
      
    end function entity_dofs
    
    subroutine calculate_dofs_per_entity
      ! Note that this will fail for meshes where not every topological
      !  entity of a given type has the same number of dofs. EG. wedge
      !  elements.
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
      integer :: ele1, ele2, vertex1, vertex2      
      integer :: entity,n
      integer, dimension(:), pointer :: neigh

      ! No need to number vertices: it comes from the topology.
      entity=entity_counts(0)

      ! Number the edges
      if (mesh_dim(mesh)>2) then
         ! If mesh_dim(mesh)==2 then this is subsumed in the facets.
         call allocate(edge_numbers, edge_list, type=CSR_INTEGER)
         call zero(edge_numbers)

         do vertex1=1, node_count(topology)
            neigh=>row_m_ptr(edge_list, vertex1)
            do n=1, size(neigh)
               vertex2=neigh(n)
               if (vertex2>vertex1) then
                  entity=entity+1
                  call set(edge_numbers,vertex1,vertex2,entity)
                  call set(edge_numbers,vertex2,vertex1,entity)
               end if
            end do
         end do
      end if

      ! Number the facets
      if (mesh_dim(mesh)>1.and.have_facets) then
         call allocate(facet_numbers, facet_list, type=CSR_INTEGER)
         call zero(facet_numbers)

         do ele1=1, element_count(mesh)
            neigh=>ele_neigh(topology,ele1)
            do n=1, size(neigh)
               ele2=neigh(n)
               if (ele2<0) then
                  ! Exterior facet
                  entity=entity+1
                  call set(facet_numbers,ele1,ele2,entity)
               else if (ele2>ele1) then
                  entity=entity+1
                  call set(facet_numbers,ele1,ele2,entity)
                  call set(facet_numbers,ele2,ele1,entity)
               end if
            end do
         end do
      end if
      
      ! There's nothing to do for cells.
      if (mesh_dim(mesh)>0) then
         entity=entity+entity_counts(mesh_dim(mesh))
      end if

      assert(sum(entity_counts)==entity)

    end subroutine number_topology

    subroutine create_topology_halos(halos)
      ! Armed with vertex, edge and facet numberings, we can work out who
      !  owns each of these objects and to which halos they may belong.
      type(halo_type), dimension(:), intent(in) :: halos
      
      integer :: face, ele1, ele2
      integer :: vertex1, vertex2
      integer :: vertices, entity, n
      integer, dimension(:), pointer :: neigh

      allocate(entity_send_targets(node_count(topology), size(halos)))
      allocate(entity_owner(node_count(topology)))
      allocate(entity_receive_level(node_count(topology)))
      vertices=entity_counts(0)
      call invert_halos(halos, entity_owner(:vertices), &
           entity_receive_level(:vertices), &
           entity_send_targets(:vertices,:))

      entity=0
      do vertex1=1,node_count(uid)
         entity=entity+1
         entity_sort_list(entity,0)=entity_owner(vertex1)
         entity_sort_list(entity,1)=node_val(uid, vertex1)
      end do
      
      ! Loop over the edges.
      if (mesh_dim(mesh)>2) then
         ! If mesh_dim(mesh)==2 then this is subsumed in the facets.

         do vertex1=1, node_count(topology)
            neigh=>row_m_ptr(edge_list, vertex1)
            do n=1, size(neigh)
               vertex2=neigh(n)
               if (vertex2>vertex1) then

                  entity=entity+1

                  ! Set the owner and any sends/receives for this edge.
                  call conduct_halo_voting(entity_owner(:vertices), &
                       entity_receive_level(:vertices), &
                       entity_send_targets(:vertices,:), &
                       [vertex1,vertex2], &
                       entity_owner(entity), &
                       entity_receive_level(entity), &
                       entity_send_targets(entity,:))
                  
                  entity_sort_list(entity,0)=entity_owner(entity)
                  entity_sort_list(entity,1:2)=&
                       node_val(uid, [vertex1,vertex2])

               end if
            end do
         end do
      end if

      ! Loop over the facets
      if (mesh_dim(mesh)>1) then
         facets=entity_counts(mesh_dim(mesh)-1)

         do ele1=1, element_count(mesh)
            neigh=>ele_neigh(topology,ele1)
            do n=1, size(neigh)
               ele2=neigh(n)
               if (ele2>ele1) then
                  ! Only need to hit each facet once. 
                  cycle
               else
                  face=ele_face(topology,ele1,ele2)

                  entity=entity+1

                  ! Set the owner and any sends/receives for this facet.
                  call conduct_halo_voting(entity_owner(:vertices), &
                       entity_receive_level(:vertices), &
                       entity_send_targets(:vertices,:), &
                       face_global_nodes(topology,face), &
                       entity_owner(entity), &
                       entity_receive_level(entity), &
                       entity_send_targets(entity,:))
                  
                  entity_sort_list(entity,0)=entity_owner(entity)
                  entity_sort_list(entity,1:face_loc(topology,face))=&
                       sorted(int(face_val(uid, face)))

               end if
            end do
         end do
      end if

      ! Loop over the cells
      do ele1=1, element_count(mesh)
         entity=entity+1

         ! Set the owner and any sends/receives for this cell.
         call conduct_halo_voting(entity_owner(:vertices), &
              entity_receive_level(:vertices), &
              entity_send_targets(:vertices,:), &
              ele_nodes(topology,ele1), &
              entity_owner(entity), &
              entity_receive_level(entity), &
              entity_send_targets(entity,:))

         entity_sort_list(entity,0)=entity_owner(entity)
         entity_sort_list(entity,1:ele_loc(topology,ele1))=&
                       sorted(int(ele_val(uid, ele1)))
      end do

      assert(entity==sum(entity_counts))

    end subroutine create_topology_halos

    subroutine deallocate_topology_halos
      
      if (.not.have_halos) return
      call deallocate(entity_send_targets)

      deallocate(entity_owner, entity_receive_level, entity_send_targets)

    end subroutine deallocate_topology_halos

    subroutine entity_order
      ! Form an orderly queue of topological entities. With halos, this
      !  ensures that all non-owned entities follow all owned entities.
      integer :: rank, i
      
      if (have_halos) then
         rank=getprocno()
         do i=1,size(entity_sort_list,1)
            ! ensure local entities come first.
            if (entity_sort_list(i,0)==rank) then
               entity_sort_list(i,0)=-1
            end if
         end do
         call sort(entity_sort_list, visit_order)
      else
         ! In the serial case, we just run through the list in order.
         visit_order=[(i,i=1,size(visit_order))]
      end if

    end subroutine entity_order

    subroutine topology_dofs
      ! For each topological entity, calculate the dofs which will lie on
      !  it.
      integer :: i, dof

      dof=0
      do i=1,size(visit_order)
         entity_dof_starts(visit_order(i))=dof
         dof=dof+dofs_per(entity_dim(visit_order(i)))
      end do
      
    end subroutine topology_dofs
    
    function entity_dim(entity)
      integer, intent(in) :: entity
      integer :: entity_dim

      integer d

      do d=0,ubound(entity_counts,1)
         if (entity<=sum(entity_counts(0:d))) then
            entity_dim=d
            return
         end if
      end do

      FLAbort("illegal entity")
      
    end function entity_dim

  end subroutine make_global_numbering_new

  function count_topological_entities(topology, edge_list, facet_list) result&
       & (entities) 
    ! Calculate the number of entities of each dimension in topology.
    type(mesh_type), intent(in) :: topology
    type(csr_sparsity), intent(in) :: edge_list, facet_list
    
    integer, dimension(0:mesh_dim(topology)) :: entities

    integer :: dim

    dim=mesh_dim(topology)

    entities(0)=node_count(topology)
    entities(mesh_dim(topology))=element_count(topology)

    if (dim>1) then
       entities(1)=entries(edge_list)/2
    end if
    if (dim>2) then
       ! (all_facets + exterior_facets)/2
       entities(2)=(entries(facet_list)+count(facet_list%colm<0))/2
    end if

  end function count_topological_entities

  subroutine invert_halos(halos, vertex_owner, receive_level, send_targets)
    ! Topology halos consist of lists of vertices to be sent or received at
    !  various halo levels from different processors. 
    !
    ! This routine inverts this: for each vertex it specifies to or from
    !  where it is sent at each halo level. It also specifies who owns the
    !   vertex.
    type(halo_type), dimension(:), intent(in), optional :: halos
    integer, dimension(:), intent(out) :: vertex_owner
    integer, dimension(:), intent(out) :: receive_level
    !! Targets to broadcast each node to. Nonods x halos
    type(integer_set), dimension(:,:), intent(out) :: send_targets

    integer :: h, n, p, rank

    rank=getprocno()

    vertex_owner = rank
    receive_level = 0

    if(.not. present(halos)) return
    do h=1,size(send_targets,1)
       do p=1,size(send_targets,2)
          call allocate(send_targets(h,p))
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
             assert(all(receive_level(halo_receives(halos(h), p)) /= 0))
          end if
#endif
          receive_level(halo_receives(halos(h), p)) = h

          do n = 1, halo_send_count(halos(h), p)
             call insert(send_targets(halo_send(halos(h), p, n),h), p)
          end do

          if(h == size(halos)) then
             ! If we're on the highest level halo, set the node owners
             vertex_owner(halo_receives(halos(h), p)) = p
#ifdef DDEBUG
          else
             ! Otherwise, these should already have been set by the higher
             ! level halo
             assert(all(vertex_owner(halo_receives(halos(h), p)) == p))
#endif
          end if

       end do proc_loop

    end do halo_loop

  end subroutine invert_halos

  
  subroutine conduct_halo_voting(vertex_owner, vertex_receive_level, vertex_send_targets,&
       & vertices, entity_owner, entity_receive_level, entity_send_targets) 
    !!< Given a list of vertices of a topological entity (element, face,
    !!< edge etc.), apply the voting algorithm to determine which
    !!< processor owns that object and therefore what the halo properties
    !!< of any new nodes placed on that object should be.
    ! List owner of each vertex.
    integer, dimension(:), intent(in) :: vertex_owner
    ! Halo receive level of each vertex.
    integer, dimension(:), intent(in) :: vertex_receive_level
    ! Send targets of each vertex at each halo level.
    type(integer_set), dimension(:,:), intent(in) :: vertex_send_targets
    integer, dimension(:), intent(in) :: vertices
    integer, intent(out) :: entity_owner
    integer, intent(out) :: entity_receive_level
    type(integer_set), dimension(:), intent(out) :: entity_send_targets

    integer :: halo, rank

    ! Defaults which may be overwritten in the following if block
    entity_receive_level=0                
    entity_owner=maxval(vertex_owner(vertices))
    rank=getprocno()

    ! Ensure the send target pointer is null
    call allocate(entity_send_targets)

    if (any(vertex_owner(vertices)/=entity_owner)) then
       ! Contested object. Nodes have the same halo properties as
       ! the winning side.

       if (entity_owner/=rank) then
          ! We don't own the new nodes so they go on the receive list.
          entity_receive_level=maxval(vertex_receive_level(vertices))

       else
          ! We do own the new nodes so we form send lists for them.
          ! 
          ! The send sets are the intersection of the send
          ! sets. This is because if one vertex is not a sender at the
          ! given halo level, the obect does not belong to that halo.
          do halo = 1, size(entity_send_targets) 
             call set_intersection(entity_send_targets(halo),&
                  vertex_send_targets(vertices, halo), &
                  mask=(vertex_owner(vertices)==rank))
          end do
       end if

       !! If we get here, the object is uncontested. It's either all ours
       !! or all someone else's.
    else if(all(vertex_receive_level(vertices)>0)) then
       ! It's someone else's.

       entity_receive_level=maxval(vertex_receive_level(vertices))

    else if (all(key_count(&
         pack(vertex_send_targets(vertices,:),mask=.true.))>0)) then
       ! It's ours and there's something to send.

       do halo = 1, size(entity_send_targets) 
          call set_intersection(entity_send_targets(halo),&
               vertex_send_targets(vertices, halo), &
               mask=(vertex_owner(vertices)==rank))
       end do

    end if

  end subroutine conduct_halo_voting

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
!!$
!!$  subroutine make_global_numbering_trace(mesh)
!!$    ! Construct a global node numbering for a trace mesh
!!$    !
!!$    ! Note that this code is broken for mixed element meshes.
!!$    type(mesh_type), intent(inout) :: mesh
!!$    !
!!$    integer :: ele, totele, ni, ele_2, current_global_index
!!$    integer, pointer, dimension(:) :: neigh
!!$    integer :: face_1,face_2,nfaces,i,face_loc, nloc
!!$
!!$    totele = mesh%elements
!!$    face_loc = mesh%faces%shape%ndof
!!$    nloc = mesh%shape%ndof
!!$
!!$    !count up how many faces there are
!!$    nfaces = 0
!!$    do ele = 1, totele
!!$       neigh => ele_neigh(mesh,ele)
!!$       do ni = 1, size(neigh)
!!$          ele_2 = neigh(ni)
!!$          if(ele_2<ele) then
!!$             nfaces=nfaces+1
!!$          end if
!!$       end do
!!$    end do
!!$    mesh%nodes = nfaces*face_loc
!!$
!!$    !construct mesh%ndglno
!!$    mesh%ndglno = 0
!!$    current_global_index = 0
!!$    do ele = 1, totele
!!$       neigh => ele_neigh(mesh,ele)
!!$       do ni = 1, size(neigh)
!!$          ele_2 = neigh(ni)
!!$          if(ele_2<ele) then
!!$             face_1=ele_face(mesh, ele, ele_2)
!!$             mesh%ndglno((ele-1)*nloc+face_local_nodes(mesh,face_1))&
!!$                  &=current_global_index+(/(i, i=1,face_loc)/)
!!$             if(ele_2>0) then
!!$                !it's not a domain boundary
!!$                !not quite sure how this works in parallel
!!$                face_2=ele_face(mesh, ele_2, ele)
!!$                mesh%ndglno((ele_2-1)*nloc+face_local_nodes(mesh,face_2))&
!!$                     &=current_global_index+(/(i, i=1,face_loc)/)
!!$             end if
!!$             current_global_index = current_global_index + &
!!$                  & mesh%faces%shape%ndof
!!$          end if
!!$       end do
!!$    end do
!!$    if(current_global_index /= mesh%nodes) then
!!$       FLAbort('bad global index count in make_global_numbering_trace')
!!$    end if
!!$    if(any(mesh%ndglno==0)) then
!!$       FLAbort('Failed to fully populate trace mesh ndglno')
!!$    end if
!!$
!!$  end subroutine make_global_numbering_trace
!!$  
!!$
!!$  subroutine make_global_numbering &
!!$       (new_nonods, new_ndglno, Nonods, Totele, NDGLNO, element, halos,&
!!$       & element_halo, new_halos) 
!!$    ! Construct the global node numbering based on the node numbering for
!!$    ! linear tets given in NDGLNO. 
!!$    integer, intent(in) :: nonods, totele
!!$    integer, intent(in), dimension(:), target :: ndglno
!!$    type(element_type), intent(in) :: element
!!$    !! The level 1 and 2 halos associated with the incoming mesh.
!!$    type(halo_type), dimension(:), intent(in), optional :: halos
!!$    !! The full element halo associated with these meshes.
!!$    type(halo_type), intent(in), optional :: element_halo
!!$    !! The level 1 and 2 halos associated with the new mesh.
!!$    type(halo_type), dimension(:), intent(out), optional :: new_halos
!!$
!!$    integer, dimension(:), intent(out) :: new_ndglno
!!$    integer, intent(out) :: new_nonods
!!$
!!$    ! Adjacency lists.
!!$    type(csr_sparsity) :: NEList, NNList, EEList
!!$
!!$    logical :: D3
!!$    integer :: dim, vertex_count, facet_vertex_count, interior_facets, exterior_facets, owned_nodes
!!$
!!$    ! Number of nodes associated with each object.
!!$    integer :: face_len, edge_len
!!$
!!$    ! Total nodes associated with an element
!!$    integer :: element_tot_len
!!$
!!$    integer :: ele, ele2, node, node2, new_node, new_node2, i, j, k, halo
!!$
!!$    integer, dimension(:), allocatable :: n
!!$
!!$    ! Process number of this processor
!!$    integer :: rank
!!$
!!$    ! Owner and halo_level of all the nodes
!!$    integer, dimension(:), allocatable :: node_owner, receive_halo_level, &
!!$         & new_receive_halo_level
!!$    integer, dimension(:,:), allocatable :: new_node_owner
!!$    integer :: this_node_owner, this_receive_halo_level
!!$
!!$    ! Cache for positions in ndglno which are currently being worked on and
!!$    ! new values to go in them
!!$    integer, dimension(:), allocatable :: ndglno_pos, ndglno_val, face_nodes
!!$
!!$    ! Map from old node numbers to new ones.
!!$    integer, dimension(:), allocatable :: node_map
!!$    
!!$    ! In each case, halos is the last dimension of halos.
!!$    type(integer_set), dimension(:), allocatable :: this_send_targets, node_targets,&
!!$         & node2_targets
!!$    type(integer_set), dimension(:,:), allocatable :: old_send_targets, new_send_targets
!!$
!!$    ! Nodes in current element
!!$    integer, pointer, dimension(:) :: ele_node
!!$
!!$    ! Nodes and elements adjacent to current node
!!$    integer, pointer, dimension(:) :: adjacent_nodes, adjacent_elements
!!$
!!$    ! Flag for whether halos are being calculated:
!!$    logical :: have_halos
!!$
!!$    ! List of vertex numbers in an element.
!!$    integer, dimension(:), pointer :: ele1_vertices, ele2_vertices
!!$
!!$    ! Ascertain whether we are calculating halos or not.
!!$    if (present(halos)) then
!!$       assert(present(new_halos))
!!$       assert(present(element_halo))
!!$       allocate(this_send_targets(size(halos)), node_targets(size(halos)),&
!!$         & node2_targets(size(halos)))
!!$
!!$       have_halos=.true.
!!$    else
!!$       allocate(this_send_targets(0))
!!$       have_halos=.false.
!!$    end if
!!$
!!$    ! Dimensionality flag.
!!$    dim = element%numbering%dimension
!!$    D3=dim==3
!!$    
!!$    ! Vertices per element.
!!$    vertex_count=element%numbering%vertices
!!$    ! Vertices per surface element
!!$    facet_vertex_count = vertex_count - 1
!!$
!!$    call MakeLists_Dynamic(Nonods, Totele, Vertex_count, ndglno, D3, NEList,&
!!$       & NNList, EEList)
!!$
!!$    new_ndglno=0
!!$
!!$    rank=getprocno()
!!$    call halo_lists(halos, node_owner, receive_halo_level,&
!!$         & old_send_targets) 
!!$
!!$    !----------------------------------------------------------------------
!!$    ! Calculate the total number of nodes by adding the vertices, edge
!!$    ! elements, face elements and interior elements.
!!$    !----------------------------------------------------------------------
!!$
!!$    
!!$    ! Add any interior elements
!!$    new_nonods=totele*element%numbering%nodes_per(dim)
!!$
!!$    if (dim>0) then
!!$       ! For dim 0 there are no facets.
!!$
!!$       ! Interior faces
!!$       interior_facets=0.5*count(EEList%colm/=0)
!!$
!!$       exterior_facets=size(EEList%colm)-2*interior_facets
!!$
!!$       new_nonods=new_nonods+element%numbering%nodes_per(dim-1)&
!!$            &*(interior_facets+exterior_facets)
!!$
!!$       if (dim>1) then
!!$          ! Vertices.
!!$          new_nonods=new_nonods+element%numbering%nodes_per(0)&
!!$               *nonods
!!$          
!!$          if (dim>2) then
!!$             ! For dim==3 only, still need edges.
!!$
!!$             ! Edges. 
!!$             ! 0.5*size(NNList%colm) is the number of edges in the mesh.
!!$             new_nonods=new_nonods+0.5*size(NNList%colm)&
!!$                  &*element%numbering%nodes_per(1)
!!$             
!!$          end if
!!$       end if
!!$
!!$    end if
!!$
!!$
!!$    ! Total nodes per element
!!$    element_tot_len=element%numbering%nodes
!!$
!!$    if (have_halos) then
!!$       allocate(new_node_owner(new_nonods, 0:size(halos)), &
!!$            &   new_receive_halo_level(size(new_ndglno)), &
!!$            &   new_send_targets(new_nonods, size(halos)))
!!$       do i = 1, size(new_send_targets,1)
!!$          do j = 1,  size(new_send_targets,2)
!!$             call allocate(new_send_targets(i,j))
!!$          end do
!!$       end do
!!$    end if
!!$
!!$    ! n is the current maximum node number.
!!$
!!$    ! We need one n plus one for each halo so as to separately number the
!!$    ! nodes in each halo. We then transplant them into one long list
!!$    ! afterwards.
!!$    if(have_halos) then
!!$       allocate(n(0:size(halos)))
!!$    else
!!$       allocate(n(0:0))
!!$    end if
!!$    n=0
!!$
!!$    !----------------------------------------------------------------------
!!$    ! Vertex numbers
!!$    !----------------------------------------------------------------------
!!$    
!!$    ! ndglno_pos stores the postion of the vertices of the current element
!!$    ! in the new_ndglno list. ndglno_val stores the corresponding values.
!!$    allocate(ndglno_pos(vertex_count), ndglno_val(vertex_count))
!!$    allocate(node_map(nonods))
!!$    node_map=0
!!$    
!!$    do i=1,totele
!!$       ele_node=>NDGLNO(vertex_count*(i-1)+1:vertex_count*i)
!!$       ndglno_pos=(i-1)*element_tot_len &
!!$            + vertex_num(ele_node, &
!!$                         ele_node, &
!!$                         element%numbering)
!!$       
!!$
!!$       ! Pick up those values which have been done already.
!!$       ndglno_val=node_map(ele_node)
!!$
!!$       do j=1,vertex_count
!!$          if (ndglno_val(j)==0) then
!!$             n(receive_halo_level(ele_node(j)))=&
!!$                  n(receive_halo_level(ele_node(j)))+1 
!!$             ndglno_val(j)=n(receive_halo_level(ele_node(j)))
!!$
!!$             node_map(ele_node(j))=ndglno_val(j)
!!$             
!!$             if(receive_halo_level(ele_node(j)) == 0 .and. have_halos) then
!!$               call copy(new_send_targets(n(0),:), &
!!$                    &    old_send_targets(ele_node(j),:))
!!$             end if
!!$          end if
!!$       end do
!!$       
!!$       new_ndglno(ndglno_pos) = ndglno_val
!!$       if (have_halos) then
!!$          do j=1,vertex_count
!!$             new_node_owner(ndglno_val(j),receive_halo_level(ele_node(j)))&
!!$                  =node_owner(ele_node(j))
!!$          end do
!!$          new_receive_halo_level(ndglno_pos)=receive_halo_level(ele_node)
!!$       end if
!!$    end do
!!$
!!$    deallocate(ndglno_pos, ndglno_val)
!!$
!!$    !----------------------------------------------------------------------
!!$    ! Edge numbers.
!!$    !----------------------------------------------------------------------
!!$    ! edge_len is the number of non-vertex nodes on each edge.
!!$    edge_len=element%numbering%nodes_per(1)
!!$
!!$    if (edge_len>0) then
!!$       allocate(ndglno_pos(edge_len))
!!$
!!$       do node=1,size(NNList,1)
!!$          adjacent_nodes=>row_m_ptr(NNlist, node)
!!$
!!$          do j=1, size(adjacent_nodes)
!!$             node2=adjacent_nodes(j)
!!$             
!!$             new_node=node_map(node)
!!$             new_node2=node_map(node2)
!!$
!!$             ! Listen very carefully, I shall do each edge only once!
!!$             if (node2<=node) cycle
!!$             
!!$             call conduct_halo_voting((/node,node2/), this_node_owner,&
!!$                  & this_receive_halo_level, this_send_targets)
!!$             
!!$             ! Now double loop over all adjacent elements and number in the
!!$             ! elements which border this edge.
!!$             adjacent_elements=>row_m_ptr(NEList, node2)
!!$             
!!$             do k=1, size(adjacent_elements)
!!$                ele=adjacent_elements(k)
!!$
!!$                if (any(ele==row_m(NEList, node))) then
!!$                   ! This element contains this edge.
!!$
!!$                   ! Nodes in this element
!!$
!!$                   ! This horrible mess finds the appropriate nodes in this
!!$                   ! element and assigns the next edge_len indices to them.
!!$                   ndglno_pos=(ele-1)*element_tot_len &
!!$                        +edge_num((/node,node2/), &
!!$                        &  NDGLNO(vertex_count*(ele-1)+1:vertex_count*ele), &
!!$                        &  element%numbering, interior=.true.)
!!$
!!$                   new_ndglno(ndglno_pos) &
!!$                        = sequence(n(this_receive_halo_level)+1, edge_len)
!!$
!!$                   if (have_halos) then
!!$                      new_receive_halo_level(ndglno_pos)= &
!!$                           this_receive_halo_level
!!$                      new_node_owner(new_ndglno(ndglno_pos),&
!!$                           & this_receive_halo_level)=this_node_owner
!!$                   end if
!!$                      
!!$                end if
!!$             end do
!!$
!!$             ! Because we by definition own send nodes, we can set these
!!$             ! directly rather than via the node element list.
!!$             if (any(key_count(this_send_targets)/=0)) then
!!$                do i=n(this_receive_halo_level)+1,&
!!$                     n(this_receive_halo_level)+edge_len
!!$
!!$                   call copy(new_send_targets(i,:), this_send_targets)
!!$                end do
!!$             end if
!!$             
!!$             ! Move on the node count.
!!$             n(this_receive_halo_level)=n(this_receive_halo_level)+edge_len
!!$             
!!$             if (have_halos) then
!!$                ! Clean up send targets
!!$                call deallocate(this_send_targets)
!!$             end if
!!$          end do
!!$
!!$       end do
!!$
!!$       deallocate(ndglno_pos)
!!$
!!$    end if
!!$
!!$    !----------------------------------------------------------------------
!!$    ! Interior face numbers - only for 3D.
!!$    ! 
!!$    ! 2D interior face number are just interior numbers and are dealt with below.
!!$    !----------------------------------------------------------------------
!!$    face_len=element%numbering%nodes_per(2)
!!$
!!$    if (D3.and.face_len>0) then
!!$       allocate(face_nodes(facet_vertex_count), ndglno_pos(face_len))
!!$
!!$       do ele=1,size(EEList,1)
!!$          adjacent_elements=>row_m_ptr(EEList, ele)
!!$          
!!$
!!$          do j=1, size(adjacent_elements)
!!$             ele2=adjacent_elements(j)
!!$
!!$             ! Skip exterior faces.
!!$             if (ele2==0) cycle
!!$
!!$             ! Listen very carefully, I shall do each face only once!
!!$             if (ele2<=ele) cycle
!!$
!!$             ele1_vertices=>NDGLNO(vertex_count*(ele-1)+1:vertex_count*ele)
!!$             ele2_vertices=>NDGLNO(vertex_count*(ele2-1)+1:vertex_count*ele2)
!!$
!!$             ! Find the vertices on the common face.
!!$             face_nodes=common(ele1_vertices,ele2_vertices)
!!$
!!$             ! Work out who should own these nodes.
!!$             call conduct_halo_voting(face_nodes, this_node_owner,&
!!$                  & this_receive_halo_level, this_send_targets)
!!$
!!$             ! Set new_ndglno for this face in ele1.
!!$             ndglno_pos=(ele-1)*element_tot_len&
!!$                  +face_num(face_nodes, ele1_vertices, &
!!$                  element%numbering, interior=.true.)
!!$             
!!$             new_ndglno(ndglno_pos) &
!!$                  = sequence(n(this_receive_halo_level)+1, face_len)
!!$
!!$             if(have_halos) then
!!$                new_receive_halo_level(ndglno_pos)=this_receive_halo_level
!!$                new_node_owner(ndglno_pos, this_receive_halo_level)&
!!$                     =this_node_owner
!!$             end if
!!$
!!$             ! Set new_ndglno for this face in ele2.
!!$             ndglno_pos=(ele2-1)*element_tot_len&
!!$                  +face_num(face_nodes, ele2_vertices, &
!!$                  element%numbering, interior=.true.)
!!$
!!$             new_ndglno(ndglno_pos) &
!!$                  = sequence(n(this_receive_halo_level)+1, face_len)
!!$
!!$             if (have_halos) then
!!$                new_receive_halo_level(ndglno_pos)=this_receive_halo_level
!!$                new_node_owner(ndglno_pos, this_receive_halo_level)&
!!$                     =this_node_owner
!!$             end if
!!$
!!$             ! Because we by definition own send nodes, we can set these
!!$             ! directly rather than via the node element list.
!!$             if (any(key_count(this_send_targets)/=0)) then
!!$                do i=n(this_receive_halo_level)+1,&
!!$                     n(this_receive_halo_level)+face_len
!!$
!!$                   call copy(new_send_targets(i,:), this_send_targets)
!!$                end do
!!$             end if
!!$             
!!$             ! Move on the node count.
!!$             n(this_receive_halo_level)=n(this_receive_halo_level)+face_len
!!$
!!$             if (have_halos) then
!!$                ! Clean up send targets
!!$                call deallocate(this_send_targets)
!!$             end if
!!$
!!$          end do
!!$       end do
!!$
!!$       deallocate(face_nodes, ndglno_pos)
!!$    end if
!!$
!!$    !----------------------------------------------------------------------
!!$    ! Remaining numbers.
!!$    !----------------------------------------------------------------------
!!$
!!$    allocate (ndglno_pos(element%ndof))
!!$    this_receive_halo_level = 0
!!$
!!$    ! This is the internal nodes of the elements plus the external faces.
!!$    do ele=1,size(EEList,1)
!!$              
!!$       ndglno_pos=sequence((ele-1)*element%ndof+1, element%ndof)
!!$
!!$       ! Work out who should own these nodes.
!!$       ele1_vertices=>NDGLNO(vertex_count*(ele-1)+1:vertex_count*ele)
!!$       call conduct_halo_voting(ele1_vertices, this_node_owner,&
!!$            & this_receive_halo_level, this_send_targets)
!!$              
!!$       do i=1, element%ndof
!!$          if (new_ndglno(ndglno_pos(i))==0) then
!!$             n(this_receive_halo_level)=n(this_receive_halo_level)+1
!!$             new_ndglno(ndglno_pos(i))=n(this_receive_halo_level)
!!$
!!$             if(have_halos) then
!!$                new_receive_halo_level(ndglno_pos(i))=this_receive_halo_level
!!$                new_node_owner(ndglno_pos(i), this_receive_halo_level)&
!!$                     =this_node_owner
!!$             end if
!!$
!!$             if (any(key_count(this_send_targets)/=0)) then
!!$                call copy(new_send_targets(n(this_receive_halo_level),:), &
!!$                     &    this_send_targets)
!!$                
!!$             end if
!!$          end if
!!$       end do
!!$
!!$       if (have_halos) then
!!$          ! Clean up send targets
!!$          call deallocate(this_send_targets)
!!$       end if
!!$       
!!$    end do
!!$
!!$    deallocate(ndglno_pos)
!!$    call deallocate(NEList)
!!$    call deallocate(NNList)
!!$    call deallocate(EEList)
!!$
!!$    ASSERT(sum(n)==new_nonods)
!!$
!!$    owned_nodes=n(0)
!!$    do i=size(n)-1,0,-1
!!$       ! Work out the offset for halo nodes.
!!$       n(i)=sum(n(0:i-1))
!!$    end do
!!$    n(0)=0
!!$
!!$    if (have_halos) then
!!$       new_ndglno=new_ndglno+n(new_receive_halo_level)
!!$
!!$       ! Repack new_node_owner into a single list.
!!$       do halo=1,size(halos)-1
!!$          new_node_owner(n(halo)+1:n(halo+1),0)&
!!$               =new_node_owner(:n(halo+1)-n(halo),halo)
!!$       end do
!!$       halo=size(halos)
!!$       new_node_owner(n(halo)+1:,0)&
!!$            =new_node_owner(:new_nonods-n(halo),halo)
!!$
!!$       call generate_new_halos(new_halos, new_ndglno, new_node_owner(:,0)&
!!$            &, new_receive_halo_level, element%ndof, owned_nodes,&
!!$            & new_send_targets) 
!!$
!!$       do i = 1, size(old_send_targets,1)
!!$          do j = 1, size(old_send_targets,2)
!!$             call deallocate(old_send_targets(i, j))
!!$          end do
!!$       end do
!!$       do i = 1, size(new_send_targets,1)
!!$          do j = 1, size(old_send_targets,2)
!!$             call deallocate(new_send_targets(i, j))
!!$          end do
!!$       end do
!!$       deallocate(old_send_targets)
!!$       deallocate(new_send_targets)
!!$    end if
!!$
!!$  contains
!!$
!!$    subroutine conduct_halo_voting(vertices, this_node_owner,&
!!$         & this_receive_halo_level, this_send_targets) 
!!$      !!< Given a list of vertices of a topological object (element, face,
!!$      !!< edge etc.), apply the voting algorithm to determine which
!!$      !!< processor owns that object and therefore what the halo properties
!!$      !!< of any new nodes placed on that object should be.
!!$      integer, dimension(:), intent(in) :: vertices
!!$      integer, intent(out) :: this_node_owner
!!$      integer, intent(out) :: this_receive_halo_level
!!$      type(integer_set), dimension(:), intent(out) :: this_send_targets
!!$
!!$      integer :: halo
!!$
!!$      ! Defaults which may be overwritten in the following if block
!!$      this_receive_halo_level=0                
!!$      this_node_owner=maxval(node_owner(vertices))
!!$
!!$      ! Ensure the send targets pointer is null
!!$      call allocate(this_send_targets)
!!$             
!!$      if (.not.have_halos) return
!!$      
!!$      if (any(node_owner(vertices)/=this_node_owner)) then
!!$         ! Contested object. Nodes have the same halo properties as
!!$         ! the winning side.
!!$       
!!$         if (this_node_owner/=rank) then
!!$            ! We don't own the new nodes so they go on the receive list.
!!$            this_receive_halo_level=maxval(receive_halo_level(vertices))
!!$
!!$         else
!!$            ! We do own the new nodes so we form send lists for them.
!!$            ! 
!!$            ! The send sets are the intersection of the send
!!$            ! sets. This is because if one vertex is not a sender at the
!!$            ! given halo level, the obect does not belong to that halo.
!!$            do halo = 1, size(halos) 
!!$               call set_intersection(this_send_targets(halo),&
!!$                    old_send_targets(vertices, halo), &
!!$                    mask=(node_owner(vertices)==rank))
!!$            end do
!!$         end if
!!$
!!$         !! If we get here, the object is uncontested. It's either all ours
!!$         !! or all someone else's.
!!$      else if(all(receive_halo_level(vertices)>0)) then
!!$         ! It's someone else's.
!!$         
!!$         this_receive_halo_level=maxval(receive_halo_level(vertices))
!!$
!!$      else if (all(key_count(&
!!$           pack(old_send_targets(vertices,:),mask=.true.))>0)) then
!!$         ! It's ours and there's something to send.
!!$         
!!$         do halo = 1, size(halos) 
!!$            call set_intersection(this_send_targets(halo),&
!!$                 old_send_targets(vertices, halo), &
!!$                 mask=(node_owner(vertices)==rank))
!!$         end do
!!$         
!!$      end if
!!$
!!$    end subroutine conduct_halo_voting
!!$
!!$    
!!$    subroutine generate_new_halos(new_halos, new_ndglno, new_node_owner&
!!$            &, new_receive_halo_level, nloc, owned_nodes, send_targets)
!!$      !!< Given the node ownership and halo level information on the new
!!$      !!< mesh, construct halos for the new mesh.
!!$      !!< Note that these halos are unsorted as sorting the halos requires
!!$      !!< the coordinate field which is not available here.
!!$      type(halo_type), dimension(:), intent(out) :: new_halos
!!$      integer, dimension(:), intent(in), target :: new_ndglno
!!$      integer, dimension(:), intent(in), target :: new_node_owner
!!$      integer, dimension(:), intent(in), target :: new_receive_halo_level
!!$      integer, intent(in) :: nloc, owned_nodes
!!$      ! Send_targets provides the list of processors to which each
!!$      ! node is broadcast at each halo level. It is nonods x halos
!!$      type(integer_set), dimension(:,:), intent(in) :: send_targets
!!$
!!$      integer :: processors, this_proc, halo, n, proc
!!$      
!!$      type(integer_set), dimension(:), allocatable :: sends, receives
!!$      integer, dimension(:), pointer :: ele_nodes
!!$
!!$      processors=getnprocs()
!!$      this_proc=getprocno()
!!$
!!$      allocate(sends(processors))
!!$      allocate(receives(processors))
!!$      do n=1,size(sends)
!!$         call allocate(sends(n))
!!$         call allocate(receives(n))
!!$      end do
!!$      
!!$      halo_loop: do halo = 1, size(new_halos)
!!$         
!!$         do n=1, size(new_receive_halo_level)
!!$            ! Receive node
!!$            if (new_receive_halo_level(n)>0.and.new_receive_halo_level(n)<=halo) then
!!$               call insert(receives(new_node_owner(new_ndglno(n))), new_ndglno(n))
!!$            end if
!!$         end do            
!!$
!!$         do n = 1, new_nonods
!!$            ! Send node
!!$            do i = 1, key_count(send_targets(n,halo))
!!$               call insert(sends(fetch(send_targets(n,halo),i)),n)
!!$            end do
!!$         end do
!!$         
!!$         call allocate(new_halos(halo), nsends=key_count(sends), &
!!$              nreceives=key_count(receives), nprocs=processors, &
!!$              nowned_nodes=owned_nodes)
!!$         
!!$         do proc = 1, processors
!!$            
!!$            new_halos(halo)%sends(proc)%ptr=set2vector(sends(proc))
!!$            call deallocate(sends(proc))
!!$            
!!$            new_halos(halo)%receives(proc)%ptr=set2vector(receives(proc))
!!$            call deallocate(receives(proc))
!!$            
!!$         end do
!!$         
!!$         call print_halo(new_halos(halo), 0)
!!$    
!!$         assert(trailing_receives_consistent(new_halos(halo)))
!!$         assert(halo_valid_for_communication(new_halos(halo)))
!!$
!!$      end do halo_loop
!!$      
!!$    end subroutine generate_new_halos
!!$
!!$    function sequence(start, len)
!!$      ! Return len consecutive integers starting at start.
!!$      integer, intent(in) :: start, len
!!$      integer, dimension(len) :: sequence
!!$      
!!$      integer :: i
!!$
!!$      forall (i=1:len)
!!$         sequence(i)=start+i-1
!!$      end forall
!!$
!!$    end function sequence
!!$
!!$    function common(list1, list2) result (list)
!!$      ! Return the common members of list1 and list2.
!!$      integer, dimension(:), intent(in) :: list1, list2
!!$      integer, dimension(count(&
!!$           (spread(list1,2,size(list2))-spread(list2,1,size(list1)))==0)) &
!!$           :: list
!!$
!!$      integer :: i, j
!!$
!!$      j=0
!!$
!!$      do i=1,size(list2)
!!$         if (any(list1(i)==list2)) then
!!$            j=j+1
!!$            list(j)=list1(i)
!!$         end if
!!$      end do
!!$
!!$      ASSERT(j==size(list))
!!$
!!$    end function common
!!$    
!!$
!!$  end subroutine make_global_numbering

end module global_numbering
