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

module SurfaceLabels
  !!< These IDs are used to indicate how surface
  !!< elements can be coarsened or refined.

  use vector_tools
  use fldebug
  use linked_lists
  use merge_tensors
  use parallel_tools
  use sparse_tools
  use adjacency_lists
  use fields
  use elements
  use vtk_interfaces
  use adjacency_lists
  use global_parameters, only : halo_tag_p
  use data_structures
  use halos
  use mpi_interfaces
  
  implicit none
  
  private
  public :: FindGeometryConstraints, &
       get_coplanar_ids, reset_coplanar_ids, &
       minimum_distance_to_line_segment
       
  ! The magic numbers corresponds to what's used in libadapt
  real, parameter:: COPLANAR_MAGIC_NUMBER=0.999999
  

contains

  subroutine FindGeometryConstraints(positions, gconstraint)
    type(vector_field), target, intent(in):: positions
    real, dimension(:), intent(out):: gconstraint(:)
    
    integer :: NNodes, NElements, SNLOC
    integer, dimension(:), allocatable:: ENList
    integer, dimension(:), pointer:: SurfaceIds
    real, dimension(:), pointer:: X, Y, Z

    integer i,j,nid,npatches,p
    type(ilist), allocatable, dimension(:)::neigh
    type(ilist) :: border, corner
    type(elist) :: geom_edges
    type(edgenode), pointer :: edge
    type(inode), pointer::node, node_a, node_b
    real bbox(6), V(9), A(3)
    real length, u1, u2, u3, min_eigenvalue, dist, max_dist, min_dist
    integer start_ele, n1, n2, ele
    real, pointer, dimension(:, :)::tensor1, gtensor

    if (positions%dim/=3) then
       FLAbort("FindGeometryConstraints currently only works in 3D.")
    end if
    
    NNodes=node_count(positions)
    NElements=surface_element_count(positions)
    SNLOC=face_loc(positions, 1)
    X => positions%val(X_)%ptr
    Y => positions%val(Y_)%ptr
    Z => positions%val(Z_)%ptr
    SurfaceIds => positions%mesh%faces%coplanar_ids
    ! get surface element (global) node list
    allocate(ENList(1:NElements*SNLOC))
    call getsndgln(positions%mesh, ENLIST)
    
    ! What surfaces meet at each node.
    allocate(neigh(NNodes))
    do i=1, NElements
       do j=1, snloc
          nid = ENList((i-1)*snloc+j)
          call insert_ascending(neigh(nid), SurfaceIds(i))
       end do
    end do

    ! Find bounding box
    bbox(1) = X(1) ; bbox(2) = X(1)
    bbox(3) = Y(1) ; bbox(4) = Y(1)
    bbox(5) = Z(1) ; bbox(6) = Z(1)
    do i=2, NNodes
       bbox(1) = min(X(i), bbox(1)) ; bbox(2) = max(X(i), bbox(2))
       bbox(3) = min(Y(i), bbox(3)) ; bbox(4) = max(Y(i), bbox(4))
       bbox(5) = min(Z(i), bbox(5)) ; bbox(6) = max(Z(i), bbox(6))
    end do
    
    ! Find the minimum eigenvalue that will be used for constraining
    ! the metric.
    max_dist = max(max(bbox(2)-bbox(1), bbox(4)-bbox(3)), bbox(6)-bbox(5))
    min_eigenvalue = 1.0/(max_dist**2)
    
    ! Initialise geometry constraints
    gconstraint = 0.0
    do i=0, NNodes-1
       gconstraint(i*9+1) = 1.0/(bbox(2)-bbox(1))**2
       gconstraint(i*9+5) = 1.0/(bbox(4)-bbox(3))**2
       gconstraint(i*9+9) = 1.0/(bbox(6)-bbox(5))**2
    end do

    allocate(tensor1(3, 3), gtensor(3, 3))
    npatches = maxval(SurfaceIds)

    ! Loop over each patch, merging in the constraints from each patch
    ! into the nodes that lie along geometry edges.
    do p=1, npatches

      !ewrite(-1,*) "p == ", p
      start_ele = -1
      ! identify border and corner nodes of this patch, na ja?
      do ele=1,NElements
        if(SurfaceIds(ele).ne.p) cycle
        if (start_ele == -1) then
          start_ele = ele
          !ewrite(-1,*) "start_ele == ", start_ele
        end if
        do j=1,snloc
          n1 = ENList((ele-1)*snloc+j)
          if (neigh(n1)%length > 1) then
            call insert(border, n1)
            if (neigh(n1)%length > 2) then
              call insert(corner, n1)
            end if
          end if
        end do
      end do
        
      ! in parallel we have global coplanar ids, that are non-consecutive on the local process
      if (start_ele == -1) cycle

      ! from the corner nodes, identify the geometry edges
      node_A => corner%firstnode
      do while(associated(node_A))
        node_B => corner%firstnode
        do while(associated(node_B))
          if (node_A%value /= node_B%value) then
            if (size_intersection(neigh(node_A%value), neigh(node_B%value)) > 1) then
              if (.not. has_value(geom_edges, node_A%value, node_B%value) .and. &
                  .not. has_value(geom_edges, node_B%value, node_A%value)) then
                !ewrite(-1,'(a, i0, a, i0, a)') "edge: (", node_A%value, ", ", node_B%value, ")"
                call insert(geom_edges, node_A%value, node_B%value)
              end if
            end if
          end if
          node_B => node_B%next
        end do
        node_A => node_A%next
      end do

      ! Let the first eigenvector be the first edge of the surface
      V(1) = x(enlist((start_ele-1)*snloc+1))-x(enlist((start_ele-1)*snloc+2))
      V(2) = y(enlist((start_ele-1)*snloc+1))-y(enlist((start_ele-1)*snloc+2))
      V(3) = z(enlist((start_ele-1)*snloc+1))-z(enlist((start_ele-1)*snloc+2))
      length = sqrt(V(1)**2+V(2)**2+V(3)**2)
      V(1) = V(1)/length
      V(2) = V(2)/length
      V(3) = V(3)/length

      !ewrite(-1,*) "1st eigenvector: ", (/v(1), v(2), v(3)/)

      ! The second eigenvector is the surface normal       
      u1 = x(enlist((start_ele-1)*snloc+1))-x(enlist((start_ele-1)*snloc+3))
      u2 = y(enlist((start_ele-1)*snloc+1))-y(enlist((start_ele-1)*snloc+3))
      u3 = z(enlist((start_ele-1)*snloc+1))-z(enlist((start_ele-1)*snloc+3))
      V(4) =  u2*V(3) - u3*V(2)
      V(5) = -u1*V(3) + u3*V(1)
      V(6) =  u1*V(2) - u2*V(1)
      length = sqrt(V(4)**2+V(5)**2+V(6)**2)
      V(4) = V(4)/length
      V(5) = V(5)/length
      V(6) = V(6)/length

      !ewrite(-1,*) "2nd eigenvector: ", (/v(4), v(5), v(6)/)

      ! The third eigenvector is the cross product of the other two
      V(7:9) = cross_product(V(1:3), V(4:6))

      !ewrite(-1,*) "3rd eigenvector: ", (/v(7), v(8), v(9)/)

      ! Loop through all boundary nodes and get distance to the
      ! geometry edges which it's not a element of.
      node => border%firstnode
      do while(associated(node))
        !ewrite(-1,*) "---------------------------------------------------"
        !ewrite(-1,*) "node == ", node%value
        !ewrite(-1,*) "position == ", x(node%value), y(node%value), z(node%value)
        min_dist = max_dist
        edge => geom_edges%firstnode
        do while(associated(edge))
          n1 = edge%i
          n2 = edge%j
          if (node%value /= n1 .and. node%value /= n2) then
            dist = minimum_distance_to_line_segment( &
                 (/X(node%value),  Y(node%value),  Z(node%value)/), &
                 (/X(n1), Y(n1), Z(n1)/), &
                 (/X(n2), Y(n2), Z(n2)/))

            !ewrite(-1,*) " -- n1: == ", n1, ";", (/X(n1), Y(n1), Z(n1)/)
            !ewrite(-1,*) " -- n2: == ", n2, ";", (/X(n2), Y(n2), Z(n2)/)
            !ewrite(-1,*) " -- dist == ", dist

            if (dist > epsilon(0.0_4)) then
              min_dist = min(min_dist, dist)
            end if
            !ewrite(-1,*) " -- min_dist == ", min_dist
          end if
          edge => edge%next
        end do

        !ewrite(-1,*) " -- min_dist == ", min_dist
        min_dist = min_dist
        ! Form eigenvalues
        A(1) = 1.0/min_dist**2
        A(2) = 1.0/min_dist**2
        !A(2) = min_eigenvalue;
        A(3) = 1.0/min_dist**2

        ! Form g-constraint tensor
        call eigenrecomposition(gtensor, reshape(V, (/3, 3/)), A)
        
        ! Merge constraints
        tensor1 = reshape(gconstraint((node%value-1)*9+1:node%value*9), (/3, 3/))
        call merge_tensor(tensor1, gtensor)

        gconstraint((node%value-1)*9+1:node%value*9) = reshape(tensor1, (/9/))          
        node => node%next
      end do
    
      call flush_list(border)
      call flush_list(corner)
      call flush_list(geom_edges)
   end do
    
    ! Delete neighbour list
    do i=1, NNodes
       call flush_list(neigh(i))
    end do
    deallocate(neigh)

    ! Delete various arrays
    deallocate(tensor1, gtensor)
    deallocate(ENList)
    
  end subroutine FindGeometryConstraints

  subroutine get_coplanar_ids(mesh, positions, coplanar_ids)
  !!< Returns a pointer to an array of coplanar ids that assigns
  !!< a unique id to each coplanar patch of the surface mesh of "mesh".
  !!< The calculated coplanar_ids are caches inside the mesh, so that if
  !!< this routine is called again they are exactly the same.
  type(mesh_type), intent(inout):: mesh
  type(vector_field), intent(in):: positions
  integer, dimension(:), pointer:: coplanar_ids
    type(mesh_type), pointer:: surface_mesh
    type(csr_sparsity), pointer :: eelist
    type(ilist) front
    real, allocatable, dimension(:,:):: normals, normalgi
    real, allocatable, dimension(:):: detwei_f
    integer, allocatable, dimension(:):: face_nodes
    real coplanar
    integer, dimension(:), pointer:: neigh
    integer current_id, sngi, stotel
    integer j, k, sele, pos, ele
    
    if (.not. has_faces(mesh)) then
       call add_faces(mesh)
    end if
    
    stotel = surface_element_count(mesh)
    if(stotel == 0) then
      sngi = 0
    else
      sngi = face_ngi(mesh, 1)
    end if
    surface_mesh => mesh%faces%surface_mesh
    
    if (.not. associated(mesh%faces%coplanar_ids)) then
       allocate(mesh%faces%coplanar_ids(1:stotel))
       coplanar_ids => mesh%faces%coplanar_ids
    else
       ! we assume they have been calculated already:
       coplanar_ids => mesh%faces%coplanar_ids
       return
    end if
    coplanar_ids => mesh%faces%coplanar_ids
    allocate( normalgi(positions%dim,sngi), &
       detwei_f(sngi), face_nodes(1:face_loc(mesh,1)) )
       
    ! Calculate element normals for all surface elements
    allocate(normals(positions%dim, stotel))
    do sele=1, stotel
       ele=face_ele(mesh, sele)
       ! we don't want to use face_val on positions 
       ! as it may not have a faces object (as opposed to the mesh)
       face_nodes=face_global_nodes(mesh, sele)
       call transform_facet_to_physical( &
          positions, sele, detwei_f, normalgi)
       ! average over gauss points:
       normals(:,sele)=matmul(normalgi, detwei_f)/sum(detwei_f)
    end do
    deallocate(normalgi, detwei_f, face_nodes)
    
    ! create element-element list for surface mesh
    eelist => extract_eelist(surface_mesh)
    
    coplanar_ids = 0
    current_id = 1
    pos = 1
    do while (.true.)
       ! Create a new starting point
       do sele=pos, stotel
          if(coplanar_ids(sele)==0) then
             ! This is the first element in the new patch
             pos = sele
             coplanar_ids(pos) = current_id
             exit
          end if
       end do
          
       ! Jump out of this while loop if we are finished
       if (sele>stotel) exit
          
       ! Initialise the front
       call insert_ascending(front, pos)
          
       ! Advance this front
       do while (front%length.ne.0)
          sele = pop(front)
          
          ! surrounding surface elements:
          neigh => row_m_ptr(eelist, sele)
          do j=1, size(neigh)
             k = neigh(j)
             if (k>0) then
                if(coplanar_ids(k)==0) then
                   coplanar = dot_product(normals(:,pos), normals(:,k))
                   if(coplanar>=COPLANAR_MAGIC_NUMBER) then
                   
                      call insert_ascending(front, k)
                      coplanar_ids(k) = current_id
                   end if
                end if
            end if
          end do
       end do
         
       current_id = current_id + 1
       pos = pos + 1
    end do
    deallocate(normals)
    
    call merge_surface_ids_flcomms(mesh, coplanar_ids, max_id = current_id - 1)
    !call vtk_write_coplanar_ids("coplanar_ids", positions, coplanar_ids)

  end subroutine get_coplanar_ids
  
  subroutine vtk_write_coplanar_ids(filename, positions, coplanar_ids)
    character(len = *), intent(in) :: filename
    type(vector_field), intent(inout) :: positions
    integer, dimension(surface_element_count(positions)), intent(in) :: coplanar_ids
    
    integer, dimension(:), allocatable :: old_coplanar_ids
    
    assert(has_faces(positions%mesh))
    if(.not. associated(positions%mesh%faces%coplanar_ids)) then
      allocate(positions%mesh%faces%coplanar_ids(size(coplanar_ids)))
      positions%mesh%faces%coplanar_ids = coplanar_ids
    
      call vtk_write_surface_mesh(filename, position = positions)
      
      deallocate(positions%mesh%faces%coplanar_ids)
      nullify(positions%mesh%faces%coplanar_ids)
    else    
      allocate(old_coplanar_ids(size(coplanar_ids)))
      old_coplanar_ids = positions%mesh%faces%coplanar_ids
      positions%mesh%faces%coplanar_ids = coplanar_ids
      
      call vtk_write_surface_mesh(filename, position = positions)
      
      positions%mesh%faces%coplanar_ids = old_coplanar_ids
      deallocate(old_coplanar_ids)
    end if
    
  end subroutine vtk_write_coplanar_ids
  
  subroutine merge_surface_ids_flcomms(mesh, surface_ids, max_id)
    !!< Given a local set of surface IDs on a mesh, merge the surface IDs
    !!< across all processes (using FLComms)
  
    type(mesh_type), intent(in) :: mesh
    integer, dimension(surface_element_count(mesh)), intent(inout) :: surface_ids
    integer, optional, intent(in) :: max_id
  
#ifdef HAVE_MPI
    integer :: ierr, lmax_id, offset, snloc, stotel
    integer, dimension(:), allocatable :: sndgln, surface_ids_copy
    
    interface
      subroutine flcomms_merge_surface_patches(tag, nselements, nlocal, senlist, ids)
        implicit none
        integer, intent(in) :: tag
        integer, intent(in) :: nselements
        integer, intent(in) :: nlocal
        integer, dimension(nselements * nlocal), intent(in) :: senlist
        integer, dimension(nselements), intent(inout) :: ids
      end subroutine flcomms_merge_surface_patches
    end interface
    
    if(IsParallel()) then
       if(present(max_id)) then
         lmax_id = max_id
       else
         lmax_id = maxval(surface_ids)
       end if
    
       stotel = surface_element_count(mesh)
       if(stotel == 0) then
         snloc = 0
       else
         snloc = face_loc(mesh, 1)
       end if
    
       ! First, ensure that each process is using a different set of
       ! numbers for numbering.
       offset = 0
       call MPI_Scan(lmax_id, offset, 1, getpinteger(), MPI_SUM, MPI_COMM_WORLD, ierr)
       assert(ierr == MPI_SUCCESS)
       
       offset = offset - lmax_id
       
       surface_ids=surface_ids + offset

       ! Ensure that co-planar patches between domains have the same surface id
       allocate(sndgln(1:stotel*snloc))
       call getsndgln(mesh, sndgln)
       ! make a copy as we can't pass a fortran pointer directly to C
       allocate(surface_ids_copy(1:stotel))
       surface_ids_copy=surface_ids
       call flcomms_merge_surface_patches(halo_tag_p, stotel, snloc, sndgln, surface_ids_copy)
       surface_ids=surface_ids_copy
       deallocate(sndgln, surface_ids_copy)
    end if
#endif

  end subroutine merge_surface_ids_flcomms
  
  subroutine merge_surface_ids(mesh, surface_ids, max_id)
    !!< Given a local set of surface IDs on a mesh, merge the surface IDs
    !!< across all processes
    
    type(mesh_type), intent(in) :: mesh
    integer, dimension(surface_element_count(mesh)), intent(inout) :: surface_ids
    integer, optional, intent(in) :: max_id
    
#ifdef HAVE_MPI
    integer :: comm, communicated_nsele, communicator, ele, face, i, id_base, &
      & ierr, j, lmax_id, new_id, nhalos, nprocs, nsele, old_id, proc, procno
    integer, dimension(:), allocatable :: nreceives, requests, statuses
    integer, dimension(:), pointer :: faces
    integer, parameter :: max_comm_count = 100, vals_per_comm = 3
    logical :: complete
    type(integer_hash_table) :: gens, id_map
    type(integer_set), dimension(:), allocatable :: send_paint, sends, sent_ids
    type(integer_vector), dimension(:), allocatable :: receive_buffer, &
      & send_buffer
    type(halo_type), pointer :: ele_halo, node_halo
    
    ewrite(1, *) "In merge_surface_ids"
    
    nhalos = halo_count(mesh)
    if(nhalos == 0) return
    node_halo => mesh%halos(nhalos)
    if(serial_storage_halo(node_halo)) return
    
    assert(element_halo_count(mesh) > 0)
    ele_halo => mesh%element_halos(element_halo_count(mesh))
    
    communicator = halo_communicator(node_halo)
    nprocs = halo_proc_count(node_halo)
    procno = getprocno(communicator = communicator)
    nsele = surface_element_count(mesh)
    
    ! First things first: Make sure all IDs are unique across all processes
    
    if(present(max_id)) then
      lmax_id = max_id
    else
      lmax_id = maxval(surface_ids)
    end if
    call mpi_scan(lmax_id, id_base, 1, getpinteger(), MPI_SUM, communicator, ierr)
    assert(ierr == MPI_SUCCESS)
    id_base = id_base - lmax_id
    surface_ids = surface_ids + id_base
    
    ! Paint the sends
    
    allocate(send_paint(ele_count(mesh)))
    call allocate(send_paint)
    do i = 1, nprocs    
      do j = 1, halo_send_count(ele_halo, i)
        call insert(send_paint(halo_send(ele_halo, i, j)), i)
      end do
    end do
    
    ! Construct lists of owned faces to send to other processes ("sends") and
    ! the number of non-owned faces to receive from other processes
    ! ("nreceives"). This is maximal surface element halo data.
    
    allocate(sends(nprocs))
    call allocate(sends)
    allocate(nreceives(nprocs))
    nreceives = 0
    do i = 1, nsele
      ele = face_ele(mesh, i)
      proc = ele_owner(ele, mesh, node_halo)
      if(proc == procno) then
        do j = 1, key_count(send_paint(ele))
          proc = fetch(send_paint(ele), j)
          assert(proc /= procno)
          call insert(sends(proc), i)
        end do
      else
        assert(proc /= procno)
        nreceives(proc) = nreceives(proc) + 1
      end if
    end do
    call deallocate(send_paint)
    deallocate(send_paint)
    
    allocate(send_buffer(nprocs))
    allocate(receive_buffer(nprocs))
    do i = 1, nprocs
      communicated_nsele = key_count(sends(i))
      allocate(send_buffer(i)%ptr(communicated_nsele * vals_per_comm))
      allocate(receive_buffer(i)%ptr(nreceives(i) * vals_per_comm))
    end do
    deallocate(nreceives)
    call get_universal_numbering_inverse(ele_halo, gens)
    allocate(sent_ids(nprocs))
    comm = 0
    comm_loop: do
      ! We loop until all new surface IDs match the incoming old surface IDs.
      ! This can take multiple communications, as we may need to merge areas of
      ! the surface on non-adjacent processes.
    
      comm = comm + 1
      if(comm > max_comm_count) then
        ! Congratulations, you have two processes more than max_comm_count
        ! partitions apart on a single surface. Increase max_comm_count or write
        ! a divide and conquer algorithm for the indirect merges.
        FLAbort("Maximum communication count encountered in merge_surface_ids")
      end if
      ewrite(2, *) "Performing surface merge ", comm
    
      ! Pack send data. For each face this consists of three values:
      !   1. The old surface ID of the face
      !   2. The universal element number of the owning element
      !   3. The local face number of the face in its owning element
      ! These are packed with the surface_ids at the head of the buffer, as
      ! only these need communicating for indirect merge checking
            
      call allocate(sent_ids)
      do i = 1, nprocs
        communicated_nsele = key_count(sends(i))
        if(comm == 1) then
          do j = 1, communicated_nsele
            face = fetch(sends(i), j)
            ! Pack surface element send buffer
            old_id = surface_ids(face)
            send_buffer(i)%ptr(j) = old_id
            call insert(sent_ids(i), old_id)
            send_buffer(i)%ptr(communicated_nsele + (j - 1) * (vals_per_comm - 1) + 1) = halo_universal_number(ele_halo, face_ele(mesh, face))
            send_buffer(i)%ptr(communicated_nsele + (j - 1) * (vals_per_comm - 1) + 2) = local_face_number(mesh, face)
          end do
        else
          do j = 1, communicated_nsele
            face = fetch(sends(i), j)
            ! Update surface element send buffer
            old_id = surface_ids(face)
            send_buffer(i)%ptr(j) = old_id
            call insert(sent_ids(i), old_id)
          end do
        end if
      end do
          
      ! Communicate the face data
                    
      allocate(requests(nprocs * 2))
      requests = MPI_REQUEST_NULL
      call mpi_barrier(communicator, ierr)
      assert(ierr == MPI_SUCCESS)
      if(comm == 1) then
        ! This is the first merge. Send everything in the buffer.
        
        do i = 1, nprocs          
          ! Non-blocking sends
          if(size(send_buffer(i)%ptr) > 0) then
            call mpi_isend(send_buffer(i)%ptr, size(send_buffer(i)%ptr), getpinteger(), i - 1, procno - 1, communicator, requests(i), ierr)
            assert(ierr == MPI_SUCCESS)
          end if
          
          ! Non-blocking receives
          if(size(receive_buffer(i)%ptr) > 0) then
            call mpi_irecv(receive_buffer(i)%ptr, size(receive_buffer(i)%ptr), getpinteger(), i - 1, i - 1, communicator, requests(i + nprocs), ierr)
            assert(ierr == MPI_SUCCESS)
          end if
        end do   
      else 
        ! This is an indirect merge check. Send just the updated surface IDs.
      
        do i = 1, nprocs          
          ! Non-blocking sends
          if(size(send_buffer(i)%ptr) > 0) then
            call mpi_isend(send_buffer(i)%ptr, size(send_buffer(i)%ptr) / vals_per_comm, getpinteger(), i - 1, procno - 1, communicator, requests(i), ierr)
            assert(ierr == MPI_SUCCESS)
          end if
          
          ! Non-blocking receives
          if(size(receive_buffer(i)%ptr) > 0) then
            call mpi_irecv(receive_buffer(i)%ptr, size(receive_buffer(i)%ptr) / vals_per_comm, getpinteger(), i - 1, i - 1, communicator, requests(i + nprocs), ierr)
            assert(ierr == MPI_SUCCESS)
          end if
        end do   
      end if
          
      ! Wait for all non-blocking communications to complete
      allocate(statuses(MPI_STATUS_SIZE * size(requests)))
      call mpi_waitall(size(requests), requests, statuses, ierr)
      assert(ierr == MPI_SUCCESS)
      deallocate(statuses)
      deallocate(requests)
      
      ! Now we invert the received data. This generates a map "id_map", mapping
      ! surface IDs to their new (merged) values. The new surface ID chosen is
      ! the lowest ID received overlaying the old surface ID.
      
      call allocate(id_map)
      do i = 1, nprocs
        communicated_nsele = size(receive_buffer(i)%ptr) / vals_per_comm
        do j = 1, communicated_nsele
          ! Unpack surface element receive buffer
          new_id = receive_buffer(i)%ptr(j)
          if(comm == 1) then
            ! Invert the received uen in the buffer, as this won't change for
            ! subsequent merges
            receive_buffer(i)%ptr(communicated_nsele + (j - 1) * (vals_per_comm - 1) + 1) = &
              & fetch(gens, receive_buffer(i)%ptr(communicated_nsele + (j - 1) * (vals_per_comm - 1) + 1))
          end if
          ele = receive_buffer(i)%ptr(communicated_nsele + (j - 1) * (vals_per_comm - 1) + 1)
          face = receive_buffer(i)%ptr(communicated_nsele + (j - 1) * (vals_per_comm - 1) + 2)
          
          faces => ele_faces(mesh, ele)
          old_id = surface_ids(faces(face))
          
          if(new_id == old_id) then
            ! This has already been merged
            cycle
          else if(new_id > old_id &
            ! A nasty corner case - it's possible for the incoming ID to be
            ! larger than the existing ID, but we didn't send the existing ID to
            ! the incoming ID sender. This can happen when partition boundaries
            ! align with surface ID boundaries.
            & .and. has_value(sent_ids(i), old_id)) then
            ! The incoming ID is larger than the existing ID, and we sent the
            ! existing ID to the incoming ID sender. The sender should be
            ! swapping out its corresponding ID for the one on this process.
            cycle
          else if(has_key(id_map, old_id)) then
            if(fetch(id_map, old_id) <= new_id) then
              ! This is already being mapped to a lower or equal ID
              cycle
            end if
          end if
          
          call insert(id_map, old_id, new_id)
        end do
      end do
      call deallocate(sent_ids)
      
      ! This is where a divide and conquer algorithm would live. We need to
      ! communicate (for comms > 1) information from the map id_map:
      !   old_id -> new_id
      ! between each process sharing a given common old_id.
      
      ewrite(2, *) "Number of merged surface IDs: ", key_count(id_map)
      complete = (key_count(id_map) == 0)
      call alland(complete, communicator = communicator)
      if(complete) then
        ! We have no more indirect merges - we're done
        call deallocate(id_map)
        exit comm_loop
      end if
      ! We're changing IDs. Hence we have to check for indirect merges.
      
      ! Remap the surface IDs
      
      do i = 1, nsele
        if(has_key(id_map, surface_ids(i))) then
          surface_ids(i) = fetch(id_map, surface_ids(i))
        end if
      end do
      call deallocate(id_map)
            
      ! We have to check for indirect merges (merges with processes that are not
      ! adjacent to this one). Let's go around again ...
    end do comm_loop
    call deallocate(sends)
    deallocate(sends)
    do i = 1, nprocs
      deallocate(send_buffer(i)%ptr)
      deallocate(receive_buffer(i)%ptr)
    end do
    deallocate(send_buffer) 
    deallocate(receive_buffer) 
    call deallocate(gens)
    deallocate(sent_ids)  
    
    ewrite(1, *) "Exiting merge_surface_ids"
#endif
    
  end subroutine merge_surface_ids
    
  subroutine reset_coplanar_ids(mesh)
  type(mesh_type), intent(inout):: mesh
    
    if (.not. has_faces(mesh)) then
       FLAbort("Need to have faces to reset_coplanar_ids")
    end if
    nullify(mesh%faces%coplanar_ids)
    
  end subroutine reset_coplanar_ids

  function distance_to_line(point, line_start, line_end) result(dist)
    real, dimension(3), intent(in) :: point, line_start, line_end
    real :: u, dist
    real, dimension(3) :: intersection, line

    line = line_end - line_start

    u = (((point(1) - line_start(1))*(line_end(1) - line_start(1))) + &
       & ((point(2) - line_start(2))*(line_end(2) - line_start(2))) + &
       & ((point(3) - line_start(3))*(line_end(3) - line_start(3)))) /  &
       dot_product(line, line)

    if (u < 0.0 .or. u > 1.0) then
      FLAbort("the perpendicular projection of the point is not on the line segment")
    end if

    intersection = line_start + u*(line_end - line_start)
    dist = norm2(point - intersection)
  end function distance_to_line

  function minimum_distance_to_line_segment(point, line_start, line_end) result(dist)
    real, dimension(3), intent(in) :: point, line_start, line_end
    real, dimension(3) :: line, vec1, vec2
    real :: epsilon, dist

    epsilon = 1e-3

    vec1 = point - line_start
    vec2 = point - line_end
    line = line_end - line_start

    if (dot_product(vec1, line) < epsilon) then
      dist = norm2(line_start - point)
    else if (-1 * dot_product(vec2, line) < epsilon) then
      dist = norm2(line_end - point)
    else
      dist = distance_to_line(point, line_start, line_end)
    end if
  end function minimum_distance_to_line_segment

end module SurfaceLabels
