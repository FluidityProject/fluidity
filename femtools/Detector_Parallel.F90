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

module detector_parallel
  use state_module
  use fields
  use spud
  use integer_hash_table_module
  use detector_data_types
  use detector_tools
  use halo_data_types
  use halos_numbering
  use mpi_interfaces

  implicit none
  
  private

  public :: distribute_detectors, exchange_detectors, register_detector_list, &
            get_num_detector_lists, get_registered_detector_lists, &
            deallocate_detector_list_array

  type(detector_list_ptr), dimension(:), allocatable, target, save :: detector_list_array
  integer :: num_detector_lists = 0

contains

  subroutine register_detector_list(detector_list)
    ! Register a detector list, so that adaptivity and Zoltan will distribute detectors on adapts
    ! This adds a pointer to the given detector list to detector_list_array and assign detector list ID
    type(detector_linked_list), target, intent(inout) :: detector_list

    type(detector_list_ptr), dimension(:), allocatable :: tmp_list_array
    integer :: i, old_size

    ! Allocate a new detector list
    if (allocated(detector_list_array)) then
       old_size = size(detector_list_array)
       allocate(tmp_list_array(old_size+1))
       do i=1, old_size
          tmp_list_array(i)%ptr=>detector_list_array(i)%ptr
       end do
       deallocate(detector_list_array)
       allocate(detector_list_array(old_size+1))
       do i=1, old_size
          detector_list_array(i)%ptr=>tmp_list_array(i)%ptr
       end do
       detector_list_array(old_size+1)%ptr=>detector_list
    else
       ! Allocate and return first detector list
       allocate(detector_list_array(1))
       detector_list_array(1)%ptr=>detector_list
    end if

    ! Advance counter and assign list ID
    num_detector_lists=num_detector_lists+1
    detector_list%id=num_detector_lists

  end subroutine register_detector_list

  subroutine get_registered_detector_lists(all_registered_lists)
    ! Return a pointer to a lists of pointers to all registered detector lists
    type(detector_list_ptr), dimension(:), pointer, intent(out) :: all_registered_lists

    all_registered_lists=>detector_list_array
  end subroutine get_registered_detector_lists

  function get_num_detector_lists()
    ! Return the number of registered detector lists
    integer :: get_num_detector_lists

    get_num_detector_lists=num_detector_lists
  end function get_num_detector_lists

  subroutine deallocate_detector_list_array()

    if (allocated(detector_list_array)) then 
       deallocate(detector_list_array)
       num_detector_lists = 0
    end if

  end subroutine deallocate_detector_list_array

  subroutine distribute_detectors(state, detector_list)
    ! Loop over all the detectors in the list and check that I own the element they are in. 
    ! If not, they need to be sent to the processor owner before adaptivity happens
    type(state_type), dimension(:), intent(in) :: state
    type(detector_linked_list), intent(inout) :: detector_list

    type(detector_linked_list), dimension(:), allocatable :: send_list_array
    type(detector_type), pointer :: detector, node_to_send
    type(vector_field), pointer :: xfield
    integer :: i, k, nprocs, all_send_lists_empty, processor_owner

    ewrite(2,*) "In distribute_detectors"  

    xfield => extract_vector_field(state(1),"Coordinate")

    ! We allocate a sendlist for every processor
    nprocs=getnprocs()
    allocate(send_list_array(nprocs))

    detector => detector_list%firstnode
    do i = 1, detector_list%length
       assert(detector%element>0)
       processor_owner=element_owner(xfield%mesh,detector%element)

       if (processor_owner/= getprocno()) then
          node_to_send => detector
          detector => detector%next

          call move(detector_list,node_to_send,send_list_array(processor_owner))
       else
          detector => detector%next
       end if
    end do

    ! Exchange detectors if there are any detectors to exchange globally
    all_send_lists_empty=0
    do k=1, nprocs
       if (send_list_array(k)%length/=0) then
          all_send_lists_empty=1
       end if
    end do
    call allmax(all_send_lists_empty)
    if (all_send_lists_empty/=0) then
       call exchange_detectors(state(1),detector_list,send_list_array)
    end if

    ! Make sure send lists are empty and deallocate them
    do k=1, nprocs
       assert(send_list_array(k)%length==0)
    end do
    deallocate(send_list_array)

    detector => detector_list%firstnode
    do i = 1, detector_list%length
       assert(element_owner(xfield%mesh,detector%element)==getprocno())
       detector=>detector%next
    end do

  end subroutine distribute_detectors

  subroutine exchange_detectors(state, detector_list, send_list_array)
    ! This subroutine serialises send_list_array, sends it, 
    ! receives serialised detectors from all procs and unpacks them.
    type(state_type), intent(in) :: state
    type(detector_linked_list), intent(inout) :: detector_list
    type(detector_linked_list), dimension(:), intent(inout) :: send_list_array

    real, dimension(:,:), allocatable :: detector_buffer
    type(detector_type), pointer :: detector, detector_received
    type(vector_field), pointer :: xfield
    type(halo_type), pointer :: ele_halo
    type(integer_hash_table) :: ele_numbering_inverse
    integer :: j, dim, count, n_stages, target_proc, receive_proc, &
               det_size, ndet_to_send, ndet_received, &
               halo_level, nprocs, IERROR
    integer, parameter ::TAG=12
    integer, dimension(:), allocatable :: sendRequest, status
    logical :: have_update_vector

    ewrite(2,*) "In exchange_detectors"  

    ! We want a sendlist for every processor
    nprocs=getnprocs()
    assert(size(send_list_array)==nprocs)

    xfield => extract_vector_field(state,"Coordinate")
    dim=xfield%dim
    allocate( sendRequest(nprocs) )

    ! Get the element halo 
    halo_level = element_halo_count(xfield%mesh)
    if (halo_level /= 0) then
       ele_halo => xfield%mesh%element_halos(halo_level)
    else
       ewrite(-1,*) "In exchange_detectors: No element halo found to translate detector%element"
       FLAbort("Exchanging detectors requires halo_level > 0")
    end if

    ! Get buffer size, depending on whether RK-GS parameters are still allocated
    have_update_vector=associated(detector_list%move_parameters)
    if(have_update_vector) then
       n_stages=detector_list%move_parameters%n_stages
       det_size=detector_buffer_size(dim,have_update_vector,n_stages)
    else
       det_size=detector_buffer_size(dim,have_update_vector)
    end if
    
    ! Send to all procs
    do target_proc=1, nprocs
       ndet_to_send=send_list_array(target_proc)%length
       allocate(detector_buffer(ndet_to_send,det_size))

       if (ndet_to_send>0) then
          ewrite(2,*) " Sending", ndet_to_send, "detectors to process", target_proc
       end if

       detector => send_list_array(target_proc)%firstnode
       if (ndet_to_send>0) then
          do j=1, send_list_array(target_proc)%length

             ! translate detector element to universal element
             assert(detector%element>0)
             detector%element = halo_universal_number(ele_halo, detector%element)

             if (have_update_vector) then
                call pack_detector(detector, detector_buffer(j,1:det_size), dim, nstages=n_stages)
             else
                call pack_detector(detector, detector_buffer(j,1:det_size), dim)
             end if

             ! delete also advances detector
             call delete(send_list_array(target_proc), detector)
          end do
       end if

       ! getprocno() returns the rank of the processor + 1, hence target_proc-1
       call MPI_ISEND(detector_buffer,size(detector_buffer), &
            & getpreal(), target_proc-1, TAG, MPI_COMM_FEMTOOLS, sendRequest(target_proc), IERROR)
       assert(ierror == MPI_SUCCESS)

       ! deallocate buffer after sending
       deallocate(detector_buffer)
    end do

    allocate( status(MPI_STATUS_SIZE) )
    call get_universal_numbering_inverse(ele_halo, ele_numbering_inverse)

    ! Receive from all procs
    do receive_proc=1, nprocs
       call MPI_PROBE(receive_proc-1, TAG, MPI_COMM_FEMTOOLS, status(:), IERROR) 
       assert(ierror == MPI_SUCCESS)

       call MPI_GET_COUNT(status(:), getpreal(), count, IERROR) 
       assert(ierror == MPI_SUCCESS)

       ndet_received=count/det_size
       allocate(detector_buffer(ndet_received,det_size))

       if (ndet_received>0) then
          ewrite(2,*) " Receiving", ndet_received, "detectors from process", receive_proc
       end if

       call MPI_Recv(detector_buffer,count, getpreal(), status(MPI_SOURCE), TAG, MPI_COMM_FEMTOOLS, MPI_STATUS_IGNORE, IERROR)
       assert(ierror == MPI_SUCCESS)

       do j=1, ndet_received
          allocate(detector_received)

          ! Unpack routine uses ele_numbering_inverse to translate universal element 
          ! back to local detector element
          if (have_update_vector) then
             call unpack_detector(detector_received,detector_buffer(j,1:det_size),dim,&
                    global_to_local=ele_numbering_inverse,coordinates=xfield,nstages=n_stages)
          else
             call unpack_detector(detector_received,detector_buffer(j,1:det_size),dim,&
                    global_to_local=ele_numbering_inverse,coordinates=xfield)
          end if

          ! If there is a list of detector names, use it, otherwise set ID as name
          if (allocated(detector_list%detector_names)) then
             detector_received%name=detector_list%detector_names(detector_received%id_number)
          else
             detector_received%name=int2str(detector_received%id_number)
          end if

          call insert(detector_list, detector_received)           
       end do
       deallocate(detector_buffer)
    end do    

    call MPI_WAITALL(nprocs, sendRequest, MPI_STATUSES_IGNORE, IERROR)
    assert(ierror == MPI_SUCCESS)

    call deallocate(ele_numbering_inverse)

    ewrite(2,*) "Exiting exchange_detectors"  

  end subroutine exchange_detectors

end module detector_parallel
