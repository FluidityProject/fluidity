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

  implicit none
  
  private

  public :: initialise_rk_guided_search, deallocate_rk_guided_search, &
            set_stage

  real, allocatable, dimension(:) :: timestep_weights
  real, allocatable, dimension(:,:) :: stage_matrix

contains

  ! Subroutine to allocate the RK stages, and update vector
  subroutine initialise_rk_guided_search(detector_list0, n_stages, dim)
    type(detector_linked_list), intent(inout) :: detector_list0
    integer, intent(in) :: n_stages, dim
      
    type(detector_type), pointer :: det0
    integer :: i,j,k,j0

    real, allocatable, dimension(:) :: stage_weights
    integer, dimension(2) :: option_rank

    ! First allocate and read stage_matrix and timestep_weights from options
    allocate(stage_weights(n_stages*(n_stages-1)/2))
    option_rank = option_shape("/io/detectors/lagrangian_timestepping/explicit_runge_kutta_g&
         &uided_search/stage_weights")
    if(option_rank(2).ne.-1) then
       FLExit('Stage Array wrong rank')
    end if
    if(option_rank(1).ne.size(stage_weights)) then
       ewrite(-1,*) 'size expected was', size(stage_weights)
       ewrite(-1,*) 'size actually was', option_rank(1)
       FLExit('Stage Array wrong size')
    end if
    call get_option("/io/detectors/lagrangian_timestepping/explicit_runge_kutta_guid&
        &ed_search/stage_weights",stage_weights)
    allocate(stage_matrix(n_stages,n_stages))
    stage_matrix = 0.
    k = 0
    do i = 1, n_stages
       do j = 1, n_stages
          if(i>j) then
             k = k + 1
             stage_matrix(i,j) = stage_weights(k)
          end if
       end do
    end do
    allocate(timestep_weights(n_stages))
    option_rank = option_shape("/io/detectors/lagrangian_timestepping/explicit_runge_kutta_g&
         &uided_search/timestep_weights")
    if(option_rank(2).ne.-1) then
       FLExit('Timestep Array wrong rank')
    end if
    if(option_rank(1).ne.size(timestep_weights)) then
       FLExit('Timestep Array wrong size')
    end if
    call get_option("/io/detectors/lagrangian_timestepping/explicit_runge_kutta_guid&
         &ed_search/timestep_weights",timestep_weights)
      
    det0 => detector_list0%firstnode
    do j0=1, detector_list0%length
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
  end subroutine initialise_rk_guided_search

  ! Subroutine to deallocate the RK stages and update vector - CJC
  subroutine deallocate_rk_guided_search(detector_list0)
    type(detector_linked_list), intent(inout) :: detector_list0
      
    type(detector_type), pointer :: det0
    integer :: j0

    if (allocated(stage_matrix)) then
       deallocate(stage_matrix)
    end if

    if (allocated(timestep_weights)) then
       deallocate(timestep_weights)
    end if
      
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
  subroutine set_stage(detector_list0,vfield,xfield,dt0,stage0,n_stages)
    type(detector_linked_list), intent(inout) :: detector_list0
    type(vector_field), pointer, intent(in) :: vfield, xfield
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

end module detector_move_rk_guided_search
