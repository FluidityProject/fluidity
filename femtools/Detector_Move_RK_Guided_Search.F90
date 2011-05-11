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

module detector_move_rk_guided_search
  use fldebug
  use fields
  use spud
  use detector_data_types
  use detector_tools
  use integer_hash_table_module

  implicit none
  
  private

  public :: allocate_rk_guided_search, deallocate_rk_guided_search, &
            set_stage, move_detectors_guided_search

contains

  ! Subroutine to allocate the RK stages and update vector
  subroutine allocate_rk_guided_search(detector_list0, dim, n_stages)
    type(detector_linked_list), intent(inout) :: detector_list0
    integer, intent(in) :: n_stages, dim
    type(detector_type), pointer :: det0

    det0 => detector_list0%firstnode
    do while (associated(det0))
       if(det0%type==LAGRANGIAN_DETECTOR) then
          if(allocated(det0%k)) then
             deallocate(det0%k)
          end if
          if(allocated(det0%update_vector)) then
             deallocate(det0%update_vector)
          end if
          allocate(det0%k(n_stages,dim))
          det0%k = 0.
          allocate(det0%update_vector(dim))
          det0%update_vector=0.
       end if
       det0 => det0%next
    end do
  end subroutine allocate_rk_guided_search

  ! Subroutine to deallocate the RK stages and update vector
  subroutine deallocate_rk_guided_search(detector_list0)
    type(detector_linked_list), intent(inout) :: detector_list0
      
    type(detector_type), pointer :: det0
    integer :: j0
      
    det0 => detector_list0%firstnode
    do j0=1, detector_list0%length
       if(det0%type==LAGRANGIAN_DETECTOR) then
          if(allocated(det0%k)) then
             deallocate(det0%k)
          end if
          if(allocated(det0%update_vector)) then
             deallocate(det0%update_vector)
          end if
       end if
       det0 => det0%next
    end do
  end subroutine deallocate_rk_guided_search

  !Subroutine to compute the vector to search for the next RK stage
  subroutine set_stage(detector_list0,vfield,xfield,dt0,stage0,n_stages, &
                       stage_matrix,timestep_weights)
    type(detector_linked_list), intent(inout) :: detector_list0
    type(vector_field), pointer, intent(in) :: vfield, xfield
    real, allocatable, dimension(:), intent(in) :: timestep_weights
    real, allocatable, dimension(:,:), intent(in) :: stage_matrix
    real, intent(in) :: dt0
    integer, intent(in) :: stage0, n_stages
    
    type(detector_type), pointer :: det0
    integer :: det_count,j0
    real, dimension(mesh_dim(xfield)+1) :: stage_local_coords
    
    det0 => detector_list0%firstnode
    do det_count=1, detector_list0%length 
       if(det0%type==LAGRANGIAN_DETECTOR) then
          det0%search_complete = .false.
          if(stage0.eq.1) then
             det0%update_vector = det0%position
          end if
          !stage vector is computed by evaluating velocity at current
          !position
          stage_local_coords=&
               local_coords(xfield,det0%element,det0%update_vector)
          det0%k(stage0,:) = &
               & eval_field(det0%element, vfield, stage_local_coords)
          if(stage0<n_stages) then
             !update vector maps from current position to place required
             !for computing next stage vector
             det0%update_vector = det0%position
             do j0 = 1, stage0
                det0%update_vector = det0%update_vector + &
                     dt0*stage_matrix(stage0+1,j0)*det0%k(j0,:)
             end do
          else
             !update vector maps from current position to final position
             det0%update_vector = det0%position
             do j0 = 1, n_stages
                det0%update_vector = det0%update_vector + &
                     dt0*timestep_weights(j0)*det0%k(j0,:)
             end do
             det0%position = det0%update_vector
          end if
       end if
       det0 => det0%next
    end do
  end subroutine set_stage

    !Subroutine to find the element containing 
    !detector%update_vector -- CJC
    !before leaving the processor or computational domain
    !detectors leaving the computational domain are set to STATIC
    !detectors leaving the processor domain are added to the list 
    !of detectors to communicate to the other processor
    !This works by searching for the element containing the next point 
    !in the RK through element faces
    !This is done by computing the local coordinates of the target point,
    !finding the local coordinate closest to -infinity
    !and moving to the element through that face
    subroutine move_detectors_guided_search(detector_list0,vfield,xfield,ihash,send_list_array,search_tolerance)
      type(detector_linked_list), intent(inout) :: detector_list0
      type(detector_linked_list), dimension(:), intent(inout) :: send_list_array
      type(vector_field), pointer, intent(in) :: vfield,xfield
      type(integer_hash_table), intent(in) :: ihash
      real, intent(in) :: search_tolerance

      type(detector_type), pointer :: det0, det_send
      integer :: det_count
      logical :: owned
      real, dimension(mesh_dim(vfield)+1) :: arrival_local_coords
      integer, dimension(:), pointer :: neigh_list
      integer :: neigh, proc_local_number
      logical :: make_static
      !
      !Loop over all the detectors
      det0 => detector_list0%firstnode
      do det_count=1, detector_list0%length
         !Only move Lagrangian detectors
         if(det0%type==LAGRANGIAN_DETECTOR.and..not.det0%search_complete) then
            search_loop: do
               !compute the local coordinates of the arrival point
               !with respect to this element
               arrival_local_coords=&
                    local_coords(xfield,det0%element,det0%update_vector)
               if(minval(arrival_local_coords)>-search_tolerance) then
                  !the arrival point is in this element
                  det0%search_complete = .true.
                  !move on to the next detector
                  det0 => det0%next
                  exit search_loop
               end if
               !the arrival point is not in this element, try to get
               ! closer to
               !it by searching in the coordinate direction in which it is
               !furthest away
               neigh = minval(minloc(arrival_local_coords))
               neigh_list=>ele_neigh(xfield,det0%element)
               if(neigh_list(neigh)>0) then
                  !the neighbouring element is also on this domain
                  !so update the element and try again
                  det0%element = neigh_list(neigh)
               else
                  !check if this element is owned (to decide where
                  !to send particles leaving the processor domain)
                  if(element_owned(vfield,det0%element)) then
                     !this face goes outside of the computational domain
                     !try all of the faces with negative local coordinate
                     !just in case we went through a corner
                     make_static=.true.
                     face_search: do neigh = 1, size(arrival_local_coords)
                        if(arrival_local_coords(neigh)<-search_tolerance&
                             & .and. neigh_list(neigh)>0) then
                           make_static = .false.
                           det0%element = neigh_list(neigh)
                           exit face_search
                        end if
                     end do face_search
                     if (make_static) then
                        det0%type=STATIC_DETECTOR
                        !move on to the next detector
                        det0%position = det0%update_vector
                        !move on to the next detector
                        det0 => det0%next
                        exit search_loop
                     end if
                  else
                     det_send => det0
                     det0 => det0%next
                     !this face goes into another computational domain
                     proc_local_number=fetch(ihash,&
                          &element_owner(vfield%mesh,det_send%element))

                     call move(detector_list0,det_send&
                          &,send_list_array(proc_local_number))
                     !move on to the next detector
                     exit search_loop
                  end if
               end if
            end do search_loop
         else
            !move on to the next detector
            det0 => det0%next
         end if
      end do
    end subroutine move_detectors_guided_search

end module detector_move_rk_guided_search
