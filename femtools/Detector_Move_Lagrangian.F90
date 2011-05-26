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

module detector_move_lagrangian
  use state_module
  use spud
  use fields
  use integer_hash_table_module
  use halo_data_types
  use halos_base
  use detector_data_types
  use detector_tools
  use detector_distribution
  use detector_move_rk_guided_search

  implicit none
  
  private

  public :: move_lagrangian_detectors, read_detector_move_options, check_any_lagrangian

contains

  ! Subroutine to allocate the detector parameters, 
  ! including RK stages and update vector
  subroutine read_detector_move_options(detector_list, detector_path)
    type(detector_linked_list), intent(inout) :: detector_list
    character(len=*), intent(in) :: detector_path

    type(detector_parameters), pointer :: parameters
    integer :: i,j,k
    real, allocatable, dimension(:) :: stage_weights
    integer, dimension(2) :: option_rank

    if (associated(detector_list%move_parameters)) then
       deallocate(detector_list%move_parameters)
    end if
    allocate(detector_list%move_parameters)
    parameters => detector_list%move_parameters

    if(have_option(trim(detector_path)//"/lagrangian_timestepping/explicit_runge_kutta_guided_search")) then

       call get_option(trim(detector_path)//"/lagrangian_timestepping/explicit_runge_kutta_guided_searc&
            h/n_stages",parameters%n_stages)
       call get_option(trim(detector_path)//"/lagrangian_timestepping/explicit_runge_kutta_guided_searc&
            h/subcycles",parameters%n_subcycles)
       call get_option(trim(detector_path)//"/lagrangian_timestepping/explicit_runge_kutta_guided_searc&
            &h/search_tolerance",parameters%search_tolerance)

       ! Allocate and read stage_matrix from options
       allocate(stage_weights(parameters%n_stages*(parameters%n_stages-1)/2))
       option_rank = option_shape(trim(detector_path)//"/lagrangian_timestepping/explicit_runge_kutta_g&
            &uided_search/stage_weights")
       if(option_rank(2).ne.-1) then
          FLExit('Stage Array wrong rank')
       end if
       if(option_rank(1).ne.size(stage_weights)) then
          ewrite(-1,*) 'size expected was', size(stage_weights)
          ewrite(-1,*) 'size actually was', option_rank(1)
          FLExit('Stage Array wrong size')
       end if
       call get_option(trim(detector_path)//"/lagrangian_timestepping/explicit_runge_kutta_guid&
            &ed_search/stage_weights",stage_weights)
       allocate(parameters%stage_matrix(parameters%n_stages,parameters%n_stages))
       parameters%stage_matrix = 0.
       k = 0
       do i = 1, parameters%n_stages
          do j = 1, parameters%n_stages
             if(i>j) then
                k = k + 1
                parameters%stage_matrix(i,j) = stage_weights(k)
             end if
          end do
       end do

       ! Allocate and read timestep_weights from options
       allocate(parameters%timestep_weights(parameters%n_stages))
       option_rank = option_shape(trim(detector_path)//"/lagrangian_timestepping/explicit_runge_kutta_g&
            &uided_search/timestep_weights")
       if(option_rank(2).ne.-1) then
          FLExit('Timestep Array wrong rank')
       end if
       if(option_rank(1).ne.size(parameters%timestep_weights)) then
          FLExit('Timestep Array wrong size')
       end if
       call get_option(trim(detector_path)//"/lagrangian_timestepping/explicit_runge_kutta_guid&
            &ed_search/timestep_weights",parameters%timestep_weights)

    else
       if (check_any_lagrangian(detector_list)) then
          ewrite(-1,*) "Found lagrangian detectors, but no timstepping options"
          FLExit('No lagrangian timestepping specified')
       end if
    end if

  end subroutine read_detector_move_options

  subroutine move_lagrangian_detectors(state, detector_list, dt, timestep, detector_names)
    type(state_type), intent(in) :: state
    type(detector_linked_list), intent(inout) :: detector_list
    real, intent(in) :: dt
    integer, intent(in) :: timestep
    character(len = FIELD_NAME_LEN), dimension(:), intent(in), optional :: detector_names

    type(detector_parameters), pointer :: parameters
    type(vector_field), pointer :: vfield, xfield
    type(detector_type), pointer :: detector
    type(detector_linked_list), dimension(:), allocatable :: send_list_array, receive_list_array
    type(integer_hash_table) :: ihash
    type(halo_type), pointer :: ele_halo
    integer :: i, j, k, num_proc, dim, number_neigh_processors, all_send_lists_empty, nprocs, &
               halo_level, stage, cycle
    logical :: any_lagrangian
    real :: rk_dt

    ewrite(2,*) "In move_lagrangian_detectors,", detector_list%length, "local detectors before moving"

    parameters => detector_list%move_parameters

    !Pull some information from state
    xfield=>extract_vector_field(state, "Coordinate")
    vfield => extract_vector_field(state,"Velocity")
    halo_level = element_halo_count(vfield%mesh)

    !making a hash table of {processor_number,count}
    !where processor_number is the number of a processor
    !which overlaps the domain of this processor
    !and count goes from 1 up to total number of overlapping processors
    !This is used to compute number_neigh_processors
    number_neigh_processors=0
    call allocate(ihash) 
    if (halo_level /= 0.) then
       ele_halo => vfield%mesh%element_halos(halo_level)
       nprocs = halo_proc_count(ele_halo)
       num_proc=1
       do i = 1, nprocs 
          if ((halo_send_count(ele_halo, i) + &
               halo_receive_count(ele_halo, i) > 0)&
               .and.(.not.has_key(ihash, i))) then
             call insert(ihash, i, num_proc)
             num_proc=num_proc+1
          end if
       end do
       number_neigh_processors=key_count(ihash)
    end if

    allocate(send_list_array(number_neigh_processors))
    allocate(receive_list_array(number_neigh_processors))

    call allocate_rk_guided_search(detector_list, xfield%dim, parameters%n_stages)
    rk_dt = dt/parameters%n_subcycles

    subcycling_loop: do cycle = 1, parameters%n_subcycles
       RKstages_loop: do stage = 1, parameters%n_stages

          call set_stage(detector_list,vfield,xfield,rk_dt,stage,parameters%n_stages, &
                       parameters%stage_matrix,parameters%timestep_weights)

          !this loop continues until all detectors have completed their
          ! timestep this is measured by checking if the send and receive
          ! lists are empty in all processors
          detector_timestepping_loop: do  

             !check if detectors any are Lagrangian
             any_lagrangian=check_any_lagrangian(detector_list)

             if (any_lagrangian) then
                !This is the actual call to move the detectors
                !The hash table is used to work out which processor
                !corresponds to which entry in the send_list_array

                !Detectors leaving the domain from non-owned elements
                !are entering a domain on another processor rather 
                !than leaving the physical domain. In this subroutine
                !such detectors are removed from the detector list
                !and added to the send_list_array

                call move_detectors_guided_search(detector_list,&
                        vfield,xfield,ihash,send_list_array,parameters%search_tolerance)

                !Work out whether all send lists are empty,
                !in which case exit.
                !This is slightly Byzantine, I think it would
                !also work if it was initially zero, and then
                !set to 1 if any of the lists were not empty. CJC
                all_send_lists_empty=number_neigh_processors
                do k=1, number_neigh_processors
                   if (send_list_array(k)%length==0) then
                      all_send_lists_empty=all_send_lists_empty-1
                   end if
                end do
                call allmax(all_send_lists_empty)
                if (all_send_lists_empty==0) exit

                !This call serialises send_list_array,
                !sends it, receives serialised receive_list_array,
                !unserialises that.
                if (present(detector_names)) then
                   call serialise_lists_exchange_receive(&
                     state,send_list_array,receive_list_array,&
                     number_neigh_processors,ihash,detector_names)
                else
                   call serialise_lists_exchange_receive(&
                     state,send_list_array,receive_list_array,&
                     number_neigh_processors,ihash)
                end if

                !This call moves detectors into the detector_list
                !I'm still unsure about how detectors are removed
                !from the detector list if they are sent.
                do i=1, number_neigh_processors
                   if  (receive_list_array(i)%length/=0) then      
                      call move_all(receive_list_array(i),detector_list)
                   end if
                end do

                !Flush the detector lists
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
             end if

          end do detector_timestepping_loop
       end do RKstages_loop
    end do subcycling_loop

    deallocate(send_list_array)
    deallocate(receive_list_array)

    !!! at the end of write_detectors subroutine I need to loop over all the detectors 
    !!! in the list and check that I own them (the element where they are). 
    !!! If not, they need to be sent to the processor owner before adaptivity happens
    if (timestep/=0) then
       call distribute_detectors(state, detector_list, ihash)
    end if

    ! This needs to be called after distribute_detectors because, since the exchange  
    ! routine serialises det%k and det%update_vector if it finds the RK-GS option
    call deallocate_rk_guided_search(detector_list)

    call deallocate(ihash) 

  end subroutine move_lagrangian_detectors

  !function to check if there are any lagrangian particles in the given list
  function check_any_lagrangian(detector_list0)
    logical :: check_any_lagrangian
    type(detector_linked_list), intent(inout) :: detector_list0
    type(detector_type), pointer :: det0
    integer :: i
    integer :: checkint 
      
    checkint = 0
    det0 => detector_list0%firstnode
    do i = 1, detector_list0%length         
       if (det0%type==LAGRANGIAN_DETECTOR) then
          checkint = 1
          exit
       end if
       det0 => det0%next
    end do
    call allmax(checkint)
    check_any_lagrangian = .false.
    if(checkint>0) check_any_lagrangian = .true.
  end function check_any_lagrangian

end module detector_move_lagrangian
