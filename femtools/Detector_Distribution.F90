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

module detector_distribution
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

  public :: distribute_detectors, serialise_lists_exchange_receive, register_detector_list

  type(detector_linked_list), dimension(:), allocatable, save, target :: detector_list_array

contains

  subroutine register_detector_list(detector_list_ptr)
    type(detector_linked_list), pointer, intent(out) :: detector_list_ptr

    type(detector_linked_list), dimension(:), allocatable :: temp_list_array
    integer :: i, old_size

    ! Allocate a new detector list
    if (allocated(detector_list_array)) then
       old_size = size(detector_list_array)
       allocate(temp_list_array(old_size+1))
       do i=1, old_size
          temp_list_array(i)=detector_list_array(i)
       end do
       detector_list_array=temp_list_array
       detector_list_ptr=>detector_list_array(old_size+1)
    else
       ! Allocate and return first detector list
       allocate(detector_list_array(1))
       detector_list_ptr=>detector_list_array(1)
    end if

  end subroutine register_detector_list

  subroutine distribute_detectors(state, detector_list, ihash, detector_names)
    ! Loop over all the detectors in the list and check that I own the element they are in. 
    ! If not, they need to be sent to the processor owner before adaptivity happens
    type(state_type), intent(in) :: state
    type(detector_linked_list), intent(inout) :: detector_list
    type(integer_hash_table), intent(in) :: ihash
    character(len = FIELD_NAME_LEN), dimension(:), intent(in), optional :: detector_names

    type(detector_linked_list), dimension(:), allocatable :: send_list_array, receive_list_array
    type(detector_type), pointer :: detector, node_to_send
    type(vector_field), pointer :: vfield
    integer :: i, k, all_send_lists_empty, list_neigh_processor, number_neigh_processors, processor_owner

    vfield => extract_vector_field(state,"Velocity")
    number_neigh_processors=key_count(ihash)
    allocate(send_list_array(number_neigh_processors))
    allocate(receive_list_array(number_neigh_processors))

    detector => detector_list%firstnode
    do i = 1, detector_list%length
       if (detector%element>0) then
          processor_owner=element_owner(vfield%mesh,detector%element)

          if (processor_owner/= getprocno()) then

             list_neigh_processor=fetch(ihash,processor_owner)
             node_to_send => detector
             detector => detector%next

             call move(detector_list,node_to_send,send_list_array(list_neigh_processor))
          else
             detector => detector%next
          end if
       else
          detector => detector%next
       end if
    end do

    all_send_lists_empty=number_neigh_processors
    do k=1, number_neigh_processors
       if (send_list_array(k)%length==0) then
          all_send_lists_empty=all_send_lists_empty-1
       end if
    end do

    call allmax(all_send_lists_empty)

    if (all_send_lists_empty/=0) then
       if (present(detector_names)) then
          call serialise_lists_exchange_receive(state,send_list_array,receive_list_array,number_neigh_processors,ihash,detector_names)
       else
          call serialise_lists_exchange_receive(state,send_list_array,receive_list_array,number_neigh_processors,ihash)
       end if
    end if

    do i=1, number_neigh_processors
       if  (receive_list_array(i)%length/=0) then      
          call move_all(receive_list_array(i),detector_list)
       end if
    end do

!!! BEFORE DEALLOCATING THE LISTS WE SHOULD MAKE SURE THEY ARE EMPTY
    do k=1, number_neigh_processors
       if (send_list_array(k)%length/=0) then 
          call delete_all(send_list_array(k))
       end if
    end do

    do k=1, number_neigh_processors
       if (receive_list_array(k)%length/=0) then  
          call delete_all(receive_list_array(k))
       end if
    end do

    deallocate(send_list_array)
    deallocate(receive_list_array)

  end subroutine distribute_detectors

  subroutine serialise_lists_exchange_receive(&
       state,send_list_array,receive_list_array,number_neigh_processors,ihash,detector_names)
    !This subroutine serialises send_list_array,
    !sends it, receives serialised receive_list_array,
    !unserialises that.
    type(state_type), intent(in) :: state
    type(detector_linked_list), dimension(:), &
         &intent(inout) :: send_list_array, receive_list_array
    integer, intent(inout) :: number_neigh_processors
    type(integer_hash_table), intent(in) :: ihash
    character(len = FIELD_NAME_LEN), dimension(:), intent(in), optional :: detector_names

    type array_ptr
       real, dimension(:,:), pointer :: ptr
    end type array_ptr

    type(array_ptr), dimension(:), allocatable :: &
         &send_list_array_serialise, receive_list_array_serialise
    type(detector_type), pointer :: detector, detector_received
    type(vector_field), pointer :: vfield, xfield
    type(halo_type), pointer :: ele_halo
    type(integer_hash_table) :: gens
    integer :: global_ele, univ_ele, &
         &number_detectors_to_send, number_of_columns, &
         &number_detectors_received, target_proc, count, IERROR, dim, i, j
    integer, PARAMETER ::TAG=12

    integer, ALLOCATABLE, DIMENSION(:) :: sendRequest
    integer, ALLOCATABLE, DIMENSION(:) :: status
 
    type(integer_hash_table) :: ihash_inverse
    integer :: halo_level, gensaa, gensaaa, target_proc_a, mapped_val_a
    integer :: ki 

    type(element_type), pointer :: shape

    type(mesh_type) :: pwc_mesh
    type(vector_field) :: pwc_positions

    !stuff for sending RK info
    logical :: have_update_vector
    integer :: update_start, k_start,n_stages

    xfield => extract_vector_field(state,"Coordinate")
    shape=>ele_shape(xfield,1)
    vfield => extract_vector_field(state,"Velocity")
    dim=vfield%dim
 
    !set up sendrequest tags because we are going to do an MPI_Isend
    !these are used by wait_all
    allocate( sendRequest(0:number_neigh_processors-1) )

    !Allocate an array of pointers for the serialised send list
    allocate(send_list_array_serialise(number_neigh_processors))

    !Get the element halo 
    halo_level = element_halo_count(vfield%mesh)
    if (halo_level /= 0) then
       ele_halo => vfield%mesh%element_halos(halo_level)
    end if

    !Get the inverse of the hash table mapping between 
    !processor numbers and numbering in the send list array
    call allocate(ihash_inverse) 
    do i=1, key_count(ihash)
       call fetch_pair(ihash, i, target_proc_a, mapped_val_a)
       call insert(ihash_inverse, mapped_val_a, target_proc_a)
    end do

    !==============================================
    !This is some kind of wierd debugging check.
    !Should we remove it? CJC
    if (halo_level /= 0) then
       pwc_mesh = piecewise_constant_mesh(&
            xfield%mesh, "PiecewiseConstantMesh")
       call allocate(pwc_positions, xfield%dim, pwc_mesh, "Coordinate")
       call deallocate(pwc_mesh)
       call remap_field(xfield, pwc_positions)
       assert(halo_verifies(ele_halo, pwc_positions))
       call deallocate(pwc_positions)
    end if
    !==============================================

    have_update_vector = &
         &have_option("/io/detectors/lagrangian_timestepping/explicit_runge_k&
         &utta_guided_search")
    if(have_update_vector) then
       call get_option("/io/detectors/lagrangian_timestepping/explicit_runge&
            &_kutta_guided_search/n_stages",n_stages)
    else
       n_stages = 1
    end if

    number_of_columns=dim+4
    if(have_update_vector) then
       !we need to include the RK update vector in serial array
       update_start = number_of_columns
       number_of_columns=number_of_columns+dim
       !we need to include the RK k values in the serial array
       k_start = number_of_columns
       number_of_columns=number_of_columns+dim*n_stages
    end if
    
    do i=1, number_neigh_processors

       detector => send_list_array(i)%firstnode
       number_detectors_to_send=send_list_array(i)%length

       allocate(send_list_array_serialise(i)%ptr(number_detectors_to_send,number_of_columns))

       if(number_detectors_to_send>0) then
          do j=1, send_list_array(i)%length

             global_ele=detector%element

             if (detector%element>0) then
                univ_ele = halo_universal_number(ele_halo, global_ele)
             else

!!! added this line below since I had to add extra code to cope with issues after adapt+zoltan. 
!!! In particular, after adapt + zoltan, the element that owns a detector can be negative, i.e., 
!!! this current proc does not see it/own it. This can happen due to floating errors if det in element boundary
                univ_ele =-1

             end if
             send_list_array_serialise(i)%ptr(j,1:dim)=detector%position
             send_list_array_serialise(i)%ptr(j,dim+1)=univ_ele
             send_list_array_serialise(i)%ptr(j,dim+2)=detector%dt
             send_list_array_serialise(i)%ptr(j,dim+3)=detector%id_number
             send_list_array_serialise(i)%ptr(j,dim+4)=detector%type
             if(have_update_vector) then
               send_list_array_serialise(i)%ptr(&
                     &j,update_start+1:update_start+dim)=&
                     &detector%update_vector
                do ki = 1, n_stages
                   send_list_array_serialise(i)%ptr(j,&
                        k_start+(ki-1)*dim+1:k_start+ki*dim) =&
                        &detector%k(ki,:)
                end do
             end if
             detector => detector%next

          end do
       end if
       target_proc=fetch(ihash_inverse, i)

       call MPI_ISEND(send_list_array_serialise(i)%ptr,size(send_list_array_serialise(i)%ptr), &
            & getpreal(), target_proc-1, TAG, MPI_COMM_FEMTOOLS, sendRequest(i-1), IERROR)
       assert(ierror == MPI_SUCCESS)
       !!!getprocno() returns the rank of the processor + 1, hence, for 4 proc, we have 1,2,3,4 whereas 
       !!!the ranks are 0,1,2,3. That is why I am using target_proc-1, so that for proc 4, it sends to 
       !!!proc with rank 3.

    end do
  
    allocate(receive_list_array_serialise(number_neigh_processors))

    allocate( status(MPI_STATUS_SIZE) )

    call get_universal_numbering_inverse(ele_halo, gens)

    ewrite(1,*) "gens length is:", key_count(gens)    

    !I THINK THIS DOES NOTHING - CJC
    do i=1, key_count(gens)
      call fetch_pair(gens, i, gensaa, gensaaa)
    end do

    do i=1, number_neigh_processors
              
       call MPI_PROBE(MPI_ANY_SOURCE, TAG, MPI_COMM_FEMTOOLS, status(:), IERROR) 
       assert(ierror == MPI_SUCCESS)

       call MPI_GET_COUNT(status(:), getpreal(), count, IERROR) 
       assert(ierror == MPI_SUCCESS)

       number_detectors_received=count/number_of_columns

       allocate(receive_list_array_serialise(i)%ptr(number_detectors_received,number_of_columns))

       call MPI_Recv(receive_list_array_serialise(i)%ptr,count, getpreal(), status(MPI_SOURCE), TAG, MPI_COMM_FEMTOOLS, MPI_STATUS_IGNORE, IERROR)
       assert(ierror == MPI_SUCCESS)

       do j=1, number_detectors_received

          univ_ele=receive_list_array_serialise(i)%ptr(j,dim+1); 

          !!! added this line below since I had to add extra code to cope
          !!! with issues after adapt+zoltan. In particular, after adapt +
          !!! zoltan, the element that owns a detector can be negative,
          !!! i.e., this current proc does not see it/own it. This can
          !!! happen due to floating errors if det in element boundary
          !!! Ana SG wrote the above, I think it can be ditched now -- cjc
          if (univ_ele/=-1) then
             global_ele=fetch(gens,univ_ele)
          else        
             global_ele=-1
          end if

          allocate(detector_received)

          allocate(detector_received%position(vfield%dim))

          detector_received%position=&
               receive_list_array_serialise(i)%ptr(j,1:dim)
          detector_received%element=global_ele
          detector_received%dt=receive_list_array_serialise(i)%ptr(j,dim+2)
          detector_received%type = &
               receive_list_array_serialise(i)%ptr(j,dim+4)
          detector_received%local = .true. 
          detector_received%id_number=&
               receive_list_array_serialise(i)%ptr(j,dim+3)
          if(have_update_vector) then
             allocate(detector_received%update_vector(vfield%dim))
             allocate(detector_received%k(n_stages,vfield%dim))
             detector_received%update_vector = &
                  receive_list_array_serialise(i)%ptr(j,update_start+1:&
                  &update_start+dim)
             do ki = 1, n_stages
                detector_received%k(ki,:) = &
                  receive_list_array_serialise(i)%ptr(&
                  j,k_start+dim*(ki-1)+1:&
                  &k_start+dim*ki)
             end do
             detector_received%search_complete=.false.
          end if
          detector_received%initial_owner=getprocno()
 
          allocate(detector_received%local_coords(local_coord_count(shape)))    
  
        !!! added this line below since I had to add extra code to cope with issues after adapt+zoltan. 
        !!! In particular, after adapt + zoltan, the element that owns a detector can be negative, i.e., this current proc does not see it/own it. 
        !!! This can happen due to floating errors if det in element boundary   

          if (detector_received%element/=-1) then 
              detector_received%local_coords=local_coords(xfield,detector_received%element,detector_received%position)
          end if

          if (present(detector_names)) then
             detector_received%name=detector_names(detector_received%id_number)
          else
             detector_received%name=int2str(detector_received%id_number)
          end if

          call insert(receive_list_array(i),detector_received) 
          
       end do

    end do    

    call MPI_WAITALL(number_neigh_processors, sendRequest, MPI_STATUSES_IGNORE, IERROR)
    assert(ierror == MPI_SUCCESS)

    call deallocate(gens)

    do i=1, number_neigh_processors

       deallocate(send_list_array_serialise(i)%ptr)
       deallocate(receive_list_array_serialise(i)%ptr)

    end do    

    deallocate(send_list_array_serialise)
    deallocate(receive_list_array_serialise)
    call deallocate(ihash_inverse) 

  end subroutine serialise_lists_exchange_receive

end module detector_distribution
