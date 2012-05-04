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
  use global_parameters, only: OPTION_PATH_LEN
  use integer_hash_table_module
  use halo_data_types
  use halos_base
  use detector_data_types
  use detector_tools
  use detector_parallel
  use detector_python
  use iso_c_binding
  use transform_elements
  use Profiler
  use linked_lists
  use integer_hash_table_module
  use pickers
  use parallel_fields, only: element_owner
  use ieee_arithmetic, only: ieee_is_nan

  implicit none
  
  private

  public :: move_lagrangian_detectors, read_detector_move_options, check_any_lagrangian, &
            read_random_walk_options

  character(len=OPTION_PATH_LEN), parameter :: rk_gs_path="/explicit_runge_kutta_guided_search"

contains

  subroutine read_detector_move_options(detector_list, options_path)
    ! Subroutine to allocate the detector parameters, 
    ! including RK stages and update vector
    type(detector_linked_list), intent(inout) :: detector_list
    character(len=*), intent(in) :: options_path

    integer :: i, j, k
    real, allocatable, dimension(:) :: stage_weights
    integer, dimension(2) :: option_rank

    if (have_option(trim(options_path))) then

       call get_option(trim(options_path)//"/subcycles",detector_list%n_subcycles)
       call get_option(trim(options_path)//"/search_tolerance",detector_list%search_tolerance)

       ! Forward Euler options
       if (have_option(trim(options_path)//"/forward_euler_guided_search")) then
          detector_list%velocity_advection = .true.
          detector_list%n_stages = 1
          allocate(detector_list%timestep_weights(detector_list%n_stages))
          detector_list%timestep_weights = 1.0
       end if

       ! Parameters for classical Runge-Kutta
       if (have_option(trim(options_path)//"/rk4_guided_search")) then
          detector_list%velocity_advection = .true.
          detector_list%n_stages = 4
          allocate(stage_weights(detector_list%n_stages*(detector_list%n_stages-1)/2))
          stage_weights = (/0.5, 0., 0.5, 0., 0., 1./)
          allocate(detector_list%stage_matrix(detector_list%n_stages,detector_list%n_stages))
          detector_list%stage_matrix = 0.
          k = 0
          do i = 1, detector_list%n_stages
             do j = 1, detector_list%n_stages
                if(i>j) then
                   k = k + 1
                   detector_list%stage_matrix(i,j) = stage_weights(k)
                end if
             end do
          end do
          allocate(detector_list%timestep_weights(detector_list%n_stages))
          detector_list%timestep_weights = (/ 1./6., 1./3., 1./3., 1./6. /)
       end if

       ! Generic Runge-Kutta options
       if (have_option(trim(options_path)//trim(rk_gs_path))) then
          detector_list%velocity_advection = .true.
          call get_option(trim(options_path)//trim(rk_gs_path)//"/n_stages",detector_list%n_stages)

          ! Allocate and read stage_matrix from options
          allocate(stage_weights(detector_list%n_stages*(detector_list%n_stages-1)/2))
          option_rank = option_shape(trim(options_path)//trim(rk_gs_path)//"/stage_weights")
          if (option_rank(2).ne.-1) then
             FLExit('Stage Array wrong rank')
          end if
          if (option_rank(1).ne.size(stage_weights)) then
             ewrite(-1,*) 'size expected was', size(stage_weights)
             ewrite(-1,*) 'size actually was', option_rank(1)
             FLExit('Stage Array wrong size')
          end if
          call get_option(trim(options_path)//trim(rk_gs_path)//"/stage_weights",stage_weights)
          allocate(detector_list%stage_matrix(detector_list%n_stages,detector_list%n_stages))
          detector_list%stage_matrix = 0.
          k = 0
          do i = 1, detector_list%n_stages
             do j = 1, detector_list%n_stages
                if(i>j) then
                   k = k + 1
                   detector_list%stage_matrix(i,j) = stage_weights(k)
                end if
             end do
          end do

          ! Allocate and read timestep_weights from options
          allocate(detector_list%timestep_weights(detector_list%n_stages))
          option_rank = option_shape(trim(options_path)//trim(rk_gs_path)//"/timestep_weights")
          if (option_rank(2).ne.-1) then
             FLExit('Timestep Array wrong rank')
          end if
          if (option_rank(1).ne.size(detector_list%timestep_weights)) then
             FLExit('Timestep Array wrong size')
          end if
          call get_option(trim(options_path)//trim(rk_gs_path)//"/timestep_weights",detector_list%timestep_weights)
       end if

       ! Boundary reflection
       if (have_option(trim(options_path)//trim("/reflect_on_boundary"))) then
          detector_list%reflect_on_boundary=.true.
       end if

       if (have_option(trim(options_path)//trim("/parametric_guided_search"))) then
          detector_list%tracking_method = GUIDED_SEARCH_TRACKING
       elseif (have_option(trim(options_path)//trim("/geometric_tracking"))) then
          detector_list%tracking_method = GEOMETRIC_TRACKING
       else
          if (check_any_lagrangian(detector_list)) then
             ewrite(-1,*) "Found lagrangian detectors, but no tracking options"
             FLExit('No lagrangian particle tracking method specified')
          end if
       end if

    else
       if (check_any_lagrangian(detector_list)) then
          ewrite(-1,*) "Found lagrangian detectors, but no timstepping options"
          FLExit('No lagrangian timestepping specified')
       end if
    end if

  end subroutine read_detector_move_options

  subroutine read_random_walk_options(detector_list, options_path)
    ! Subroutine to set the meta-parameters concerning random walks
    type(detector_linked_list), intent(inout) :: detector_list
    character(len=*), intent(in) :: options_path

    character(len=OPTION_PATH_LEN) :: options_buffer
    integer :: i, n_random_walks

    n_random_walks = option_count(trim(options_path)//"/random_walk")
    allocate(detector_list%random_walks(n_random_walks))

    do i=1, n_random_walks
       write(options_buffer, "(a,i0,a)") trim(options_path)//"/random_walk[",i-1,"]"
       call get_option(trim(options_buffer)//"/name", detector_list%random_walks(i)%name)

       ! Python Random Walk
       if (have_option(trim(options_buffer)//"/python")) then 
          detector_list%random_walks(i)%python_random_walk = .true.
          call get_option(trim(options_buffer)//"/python", detector_list%random_walks(i)%python_code)
       end if

       ! Internal Naive Random Walk
       if (trim(detector_list%random_walks(i)%name) == "Naive") then 
          detector_list%random_walks(i)%naive_random_walk = .true.
          call get_option(trim(options_buffer)//"/diffusivity_field", detector_list%random_walks(i)%diffusivity_field)
       end if

       ! Internal Diffusive Random Walk
       if (trim(detector_list%random_walks(i)%name) == "Diffusive") then 
          detector_list%random_walks(i)%diffusive_random_walk = .true.
          call get_option(trim(options_buffer)//"/diffusivity_field", detector_list%random_walks(i)%diffusivity_field)
          call get_option(trim(options_buffer)//"/diffusivity_gradient", detector_list%random_walks(i)%diffusivity_grad)

          ! Flag if we want automatic subcycling
          if (have_option(trim(options_buffer)//"/auto_subcycle")) then 
             detector_list%random_walks(i)%auto_subcycle = .true.
             call get_option(trim(options_buffer)//"/auto_subcycle/diffusivity_2nd_gradient", &
                      detector_list%random_walks(i)%diffusivity_2nd_grad)
             call get_option(trim(options_buffer)//"/auto_subcycle/scale_factor", &
                      detector_list%random_walks(i)%subcycle_scale_factor)
          end if
       end if
    end do

  end subroutine read_random_walk_options

  subroutine move_lagrangian_detectors(state, detector_list, dt, timestep)
    type(state_type), dimension(:), intent(in) :: state
    type(detector_linked_list), intent(inout) :: detector_list
    real, intent(in) :: dt
    integer, intent(in) :: timestep

    type(vector_field), pointer :: vfield, xfield
    type(detector_type), pointer :: detector, move_detector
    type(halo_type), pointer :: ele_halo
    integer :: dim, stage, cycle, rw, i
    real :: sub_dt

    ! Random Walk velocity source
    real, dimension(:), allocatable :: rw_displacement
    type(scalar_field), pointer :: diffusivity_field
    type(vector_field), pointer :: diffusivity_grad, diffusivity_2nd_grad

    call profiler_tic(trim(detector_list%name)//"::movement")

    ewrite(1,*) "In move_lagrangian_detectors for detectors list: ", detector_list%name
    ewrite(2,*) "Detector list", detector_list%id, "has", detector_list%length, &
         "local and", detector_list%total_num_det, "global detectors"

    ! Pull some information from state
    xfield=>extract_vector_field(state(1), "Coordinate")
    vfield=>extract_vector_field(state(1), "Velocity")
    allocate(rw_displacement(xfield%dim))

    ! Setup for Random Walk schemes
    if (allocated(detector_list%random_walks)) then
       do rw=1, size(detector_list%random_walks)

          ! For Python Random Walk run the user code that pulls fields from state
          if (detector_list%random_walks(rw)%python_random_walk) then
             call python_run_detector_string(trim(detector_list%random_walks(rw)%python_code), &
                      trim(detector_list%name), trim(detector_list%random_walks(rw)%name))
          end if

          ! For hardcoded Random Walks pull the relevant fields from state
          if (detector_list%random_walks(rw)%naive_random_walk) then
             diffusivity_field=>extract_scalar_field(state(1), trim(detector_list%random_walks(rw)%diffusivity_field))
          end if

          if (detector_list%random_walks(rw)%diffusive_random_walk) then
             diffusivity_field=>extract_scalar_field(state(1), trim(detector_list%random_walks(rw)%diffusivity_field))
             diffusivity_grad=>extract_vector_field(state(1), trim(detector_list%random_walks(rw)%diffusivity_grad))

             if (detector_list%random_walks(rw)%auto_subcycle) then
                diffusivity_2nd_grad=>extract_vector_field(state(1), trim(detector_list%random_walks(rw)%diffusivity_2nd_grad))
             end if
          end if
       end do
    end if

    ! Allocate the update_vector for each detector
    ! and det%k if RK velocity advection is switched on
    detector => detector_list%first
    do while (associated(detector))
       if(detector%type==LAGRANGIAN_DETECTOR) then

          if(allocated(detector%update_vector)) then
             deallocate(detector%update_vector)
          end if
          allocate(detector%update_vector(xfield%dim))
          detector%update_vector = 0.

          if (detector_list%velocity_advection) then
             if (allocated(detector%k)) then
                deallocate(detector%k)
             end if
             allocate(detector%k(detector_list%n_stages,xfield%dim))
             detector%k = 0.
          end if

          if (detector_list%tracking_method == GEOMETRIC_TRACKING) then
             allocate(detector%ray_o(xfield%dim))
             detector%ray_o = detector%position
             allocate(detector%ray_d(xfield%dim))
          end if

       end if
       detector => detector%next
    end do

    ! This is the outer, user-defined subcycling loop
    sub_dt = dt / detector_list%n_subcycles
    subcycling_loop: do cycle = 1, detector_list%n_subcycles

       ! Reset update_vector to position
       detector => detector_list%first
       do while (associated(detector))
          if(detector%type==LAGRANGIAN_DETECTOR) then
             detector%update_vector = detector%position

             do i=1, size(detector%update_vector)
                if (ieee_is_nan(detector%update_vector(i))) then
                   FLAbort("Encountered NaN detector position")
                end if
             end do
          end if
          detector => detector%next
       end do

       ! Explicit Runge-Kutta iterations
       if (detector_list%velocity_advection) then
          call profiler_tic(trim(detector_list%name)//"::movement::Runge-Kutta-"//trim(int2str(detector_list%n_stages)))
          RKstages_loop: do stage = 1, detector_list%n_stages

             ! Compute the update vector for the current stage
             call set_stage(detector_list,vfield,xfield,sub_dt,stage)

             ! Move update parametric detector coordinates according to update_vector
             ! If this takes a detector across parallel domain boundaries
             ! the routine will also send the detector
             call track_detectors(detector_list,xfield)

          end do RKstages_loop
          call profiler_toc(trim(detector_list%name)//"::movement::Runge-Kutta-"//trim(int2str(detector_list%n_stages)))
       end if

       if (allocated(detector_list%random_walks)) then
          do rw=1, size(detector_list%random_walks)
             call profiler_tic(trim(detector_list%name)//"::movement::"//trim(detector_list%random_walks(rw)%name))

             ! Apply single-cycle Random Walks
             if (.not. detector_list%random_walks(rw)%auto_subcycle) then

                detector => detector_list%first
                do while (associated(detector))
                   if (detector%type==LAGRANGIAN_DETECTOR) then

                      ! Internal Naive
                      if (detector_list%random_walks(rw)%naive_random_walk) then
                         call naive_random_walk(detector, sub_dt, xfield, diffusivity_field, rw_displacement)

                      ! Internal Diffusive
                      elseif (detector_list%random_walks(rw)%diffusive_random_walk) then
                         call diffusive_random_walk(detector, sub_dt, xfield, diffusivity_field, &
                                  diffusivity_grad, rw_displacement)

                      ! Python
                      elseif (detector_list%random_walks(rw)%python_random_walk) then
                         call python_run_random_walk(detector,detector_list%fgroup,xfield,sub_dt, trim(detector_list%name), &
                                  trim(detector_list%random_walks(rw)%name),rw_displacement)
                      end if

                      detector%update_vector=detector%update_vector + rw_displacement
                   end if
                   detector => detector%next
                end do

                call track_detectors(detector_list,xfield)

             ! Internal Diffusive Random Walk with automated sub-cycling
             elseif (detector_list%random_walks(rw)%diffusive_random_walk &
                         .and. detector_list%random_walks(rw)%auto_subcycle) then

                call auto_subcycle_random_walk(detector_list,sub_dt,detector_list%random_walks(rw)%subcycle_scale_factor,&
                          xfield,diffusivity_field,diffusivity_grad,diffusivity_2nd_grad)

             else
                FLExit("Auto-subcycling can only be used with internal diffusive Random Walk.")
             end if

             call profiler_toc(trim(detector_list%name)//"::movement::"//trim(detector_list%random_walks(rw)%name))
          end do
       end if 

       ! After everything is done, we update detector%position
       detector => detector_list%first
       do while (associated(detector))
          if (detector%type==LAGRANGIAN_DETECTOR) then
             detector%position = detector%update_vector
          end if
          detector => detector%next
       end do

    end do subcycling_loop

    deallocate(rw_displacement)

    ! Make sure all local detectors are owned and distribute the ones that 
    ! stoppped moving in a halo element
    call distribute_detectors(state(1), detector_list)

    ! Deallocate update_vector and det%k after distribute_detectors 
    ! because the exchange routine serialises them if it finds the RK-GS option
    detector => detector_list%first
    do while (associated(detector))
       if(detector%type==LAGRANGIAN_DETECTOR) then
          if(allocated(detector%k)) then
             deallocate(detector%k)
          end if
          if(allocated(detector%update_vector)) then
             deallocate(detector%update_vector)
          end if
          if(allocated(detector%ray_o)) then
             deallocate(detector%ray_o)
          end if
          if(allocated(detector%ray_d)) then
             deallocate(detector%ray_d)
          end if
       end if
       detector => detector%next
    end do

    ewrite(2,*) "After moving and distributing we have", detector_list%length, &
         "local and", detector_list%total_num_det, "global detectors"
    ewrite(1,*) "Exiting move_lagrangian_detectors"

    call profiler_toc(trim(detector_list%name)//"::movement")

  end subroutine move_lagrangian_detectors

  subroutine set_stage(detector_list,vfield,xfield,dt,stage)
    ! Compute the vector to search for in the next RK stage
    ! If this is the last stage, update detector position
    type(detector_linked_list), intent(inout) :: detector_list
    type(vector_field), pointer, intent(in) :: vfield, xfield
    real, intent(in) :: dt
    integer, intent(in) :: stage
    
    type(detector_type), pointer :: detector
    integer :: j0
    real, dimension(mesh_dim(xfield)+1) :: stage_local_coords
    
    detector => detector_list%first
    do while (associated(detector))

       if(detector%type==LAGRANGIAN_DETECTOR) then

          ! Evaluate velocity at update_vector and set k
          if (detector_list%velocity_advection) then
             ! stage vector is computed by evaluating velocity at current position
             stage_local_coords=local_coords(xfield,detector%element,detector%update_vector)
             detector%k(stage,:)=eval_field(detector%element, vfield, stage_local_coords)
          else
             ! do not advect detector with the velocity field
             detector%k(stage,:)=0.0
          end if

          if(stage<detector_list%n_stages) then
             ! Update vector maps from current position to place required
             ! for computing next stage vector
             detector%update_vector = detector%position
             do j0 = 1, stage
                detector%update_vector = detector%update_vector + dt*detector_list%stage_matrix(stage+1,j0)*detector%k(j0,:)
             end do
          else
             ! Update vector maps from current position to final position
             detector%update_vector = detector%position
             do j0 = 1, detector_list%n_stages
                detector%update_vector = detector%update_vector + dt*detector_list%timestep_weights(j0)*detector%k(j0,:)
             end do
          end if

       end if
       detector => detector%next
    end do
  end subroutine set_stage

  subroutine track_detectors(detector_list,xfield)
    !Subroutine to find the element containing the update vector:
    ! - Detectors leaving the computational domain are set to STATIC
    ! - Detectors leaving the processor domain are added to the list 
    !   of detectors to communicate to the other processor.
    !   This works by searching for the element containing the next point 
    !   in the RK through element faces.
    !   This is done by computing the local coordinates of the target point,
    !   finding the local coordinate closest to -infinity
    !   and moving to the element through that face.
    type(detector_linked_list), intent(inout) :: detector_list
    type(vector_field), pointer, intent(in) :: xfield

    type(detector_type), pointer :: detector, send_detector
    type(detector_linked_list), dimension(:), allocatable :: send_list_array
    real :: search_tolerance
    integer :: k, nprocs, new_owner, all_send_lists_empty
    logical :: outside_domain, any_lagrangian, have_ray

    call profiler_tic(trim(detector_list%name)//"::movement::tracking")

    ! We allocate a sendlist for every processor
    nprocs=getnprocs()
    allocate(send_list_array(nprocs))

    search_tolerance=detector_list%search_tolerance

    if (detector_list%tracking_method == GEOMETRIC_TRACKING) then
       have_ray = .true.

       detector => detector_list%first
       do while (associated(detector))

          ! Calcualte and normalise ray direction
          detector%ray_d = detector%update_vector - detector%ray_o
          detector%target_distance = sum(detector%ray_d**2)
          if (detector%target_distance > 0.0) then
             detector%target_distance = sqrt(detector%target_distance)
             detector%ray_d = detector%ray_d / detector%target_distance
          else
             detector%target_distance = 0.0
             detector%ray_d = 0.0
          end if
          detector%current_t = 0.0

          call flush_list(detector%ele_path_list)
          call flush_list(detector%ele_dist_list)

          detector => detector%next
       end do
    else
       have_ray = .false.
    end if

    ! This loop continues until all detectors have completed their
    ! timestep this is measured by checking if the send and receive
    ! lists are empty in all processors
    detector_timestepping_loop: do  

       ! Make sure we still have lagrangian detectors
       any_lagrangian=check_any_lagrangian(detector_list)
       if (any_lagrangian) then

          !Loop over all the detectors
          detector => detector_list%first
          do while (associated(detector))

             !Only move Lagrangian detectors
             if (.not. detector%type==LAGRANGIAN_DETECTOR) then
                detector => detector%next
                cycle
             end if

             if (detector_list%tracking_method == GUIDED_SEARCH_TRACKING) then
                call local_guided_search(xfield,detector%update_vector,detector%element,&
                        search_tolerance,new_owner,detector%local_coords)

             elseif (detector_list%tracking_method == GEOMETRIC_TRACKING) then
                call geometric_ray_tracing(xfield,detector%ray_o,detector%ray_d,detector%target_distance,&
                        detector%element,new_owner,detector%current_t,search_tolerance, &
                        detector%ele_path_list, detector%ele_dist_list)
                detector%local_coords = local_coords(xfield,detector%element,detector%update_vector)

                if (allocated(detector%ele_path)) then
                   deallocate(detector%ele_path)
                end if
                allocate(detector%ele_path(detector%ele_path_list%length))
                detector%ele_path = list2vector(detector%ele_path_list)

                if (allocated(detector%ele_dist)) then
                   deallocate(detector%ele_dist)
                end if
                allocate(detector%ele_dist(detector%ele_dist_list%length))
                detector%ele_dist = list2vector(detector%ele_dist_list)

             else
                FLExit('No lagrangian particle tracking method specified')
             end if

             if (new_owner==-1) then
                if (detector_list%reflect_on_boundary) then
                   ! We reflect the detector path at the face we just went through
                   call reflect_on_boundary(xfield,detector%update_vector,detector%element)

                else
                   ! Turn detector static inside the domain
                   ewrite(1,*) "WARNING: detector attempted to leave computational domain;"
                   ewrite(1,*) "Turning detector static, ID:", detector%id_number, "element:", detector%element
                   detector%type=STATIC_DETECTOR
                   ! move on to the next detector, without updating det0%position, 
                   ! because det0%update_vector is by now outside of the computational domain
                   detector => detector%next
                end if
             end if

             if (new_owner/=getprocno().and.new_owner>0) then
                send_detector => detector
                detector => detector%next
                call move(send_detector, detector_list, send_list_array(new_owner))
             end if

             if (new_owner==getprocno()) then
                !move on to the next detector
                detector => detector%next
             end if
          end do

          ! Work out whether all send lists are empty, in which case exit.
          all_send_lists_empty=0
          do k=1, nprocs
             if (send_list_array(k)%length/=0) then
                all_send_lists_empty=1
             end if
          end do
          call allmax(all_send_lists_empty)
          if (all_send_lists_empty==0) exit

          !This call serialises send_list_array, sends it, 
          !receives serialised receive_list_array, and unserialises that.
          call exchange_detectors(detector_list,xfield,send_list_array,have_rk=.true.,have_ray=have_ray)
       else
          ! If we run out of lagrangian detectors for some reason, exit the loop
          exit
       end if
    end do detector_timestepping_loop

    deallocate(send_list_array)

    if (detector_list%tracking_method == GEOMETRIC_TRACKING) then
       detector => detector_list%first
       do while (associated(detector))
          detector%ray_o = detector%update_vector
          detector => detector%next
       end do
    end if

    call profiler_toc(trim(detector_list%name)//"::movement::tracking")

  end subroutine track_detectors

  subroutine local_guided_search(xfield,coordinate,element,search_tolerance,new_owner,l_coords,ele_path)
    ! Do a local guided search until we either hit the new element 
    ! or have to leave the local domain. The new_owner argument indicates
    ! whether we have found the new element, where to continue searching
    ! or flags that we have gone outside the domain with -1.
    type(vector_field), pointer, intent(in) :: xfield
    real, dimension(mesh_dim(xfield)), intent(inout) :: coordinate
    integer, intent(inout) :: element
    real, intent(in) :: search_tolerance
    real, dimension(mesh_dim(xfield)+1), intent(out) :: l_coords
    integer, intent(out) :: new_owner
    type(ilist), intent(inout), optional :: ele_path

    integer, dimension(:), pointer :: neigh_list
    integer :: neigh, face
    logical :: outside_domain

    search_loop: do
       ! Compute the local coordinates of the arrival point with respect to this element
       l_coords=local_coords(xfield,element,coordinate)

       if (minval(l_coords)>-search_tolerance) then
          !The arrival point is in this element, we're done
          new_owner=getprocno()
          exit search_loop
       end if

       ! The arrival point is not in this element, try to get closer to it by 
       ! searching in the coordinate direction in which it is furthest away
       neigh = minval(minloc(l_coords))
       neigh_list=>ele_neigh(xfield,element)
       if (neigh_list(neigh)>0) then
          ! The neighbouring element is also on this domain
          ! so update the element and try again
          element = neigh_list(neigh)

          ! Record the elements along the path travelled
          if (present(ele_path)) then
             call insert(ele_path, element)
          end if
       else
          ! Next element in coordinate direction is not on this domain.
          ! Two cases arise: 
          !    * If we don't own the current element we're on a Halo.
          !      We need to continue searching on another processor.
          !    * If we own the current element the detector will move outside the domain.
          !      But we need to check if we've gone through a corner.

          ! So check if this element is owned
          if (element_owned(xfield,element)) then
             ! This face goes outside of the computational domain.
             ! Try all of the faces with negative local coordinate
             ! just in case we went through a corner.
             outside_domain=.true.
             face_search: do face = 1, size(l_coords)
                if (l_coords(face)<-search_tolerance.and.neigh_list(face)>0) then
                   outside_domain = .false.
                   element = neigh_list(face)
                   exit face_search
                end if
             end do face_search

             if (outside_domain) then
                new_owner=-1
                exit search_loop
             end if
          else
             ! The current element is on a Halo, we need to send it to the owner.
             new_owner=element_owner(xfield%mesh,element)
             exit search_loop
          end if
       end if
    end do search_loop

  end subroutine local_guided_search

  subroutine geometric_ray_tracing(xfield,r_o,r_d,target_distance,new_element,new_owner,current_t,search_tolerance,ele_path,ele_dist)
    ! This tracking method is based on a standard Ray-tracing algorithm using planes and half-spaces.
    ! Reference: 
    type(vector_field), pointer, intent(in) :: xfield
    real, dimension(mesh_dim(xfield)), intent(in) :: r_d, r_o
    integer, intent(inout) :: new_element
    integer, intent(out) :: new_owner
    real, intent(in) :: target_distance, search_tolerance
    real, intent(inout) :: current_t
    type(ilist), intent(inout), optional :: ele_path
    type(rlist), intent(inout), optional :: ele_dist

    real :: face_t, ele_t
    integer :: i, neigh_face, next_face
    integer, dimension(:), pointer :: face_list

    ! Exit if r_d is zero
    if (target_distance < search_tolerance) then
       new_owner=getprocno()
       return
    end if

    call insert(ele_path, new_element)

    search_loop: do 

       ! Go through all faces of the element and look for the smallest t in the ray's direction
       next_face = -1
       ele_t = huge(1.0)
       face_list=>ele_faces(xfield,new_element)
       do i=1, size(face_list)
          call ray_intersetion_distance(face_list(i), face_t)

          if (face_t < ele_t .and. current_t < face_t) then
             ele_t = face_t
             next_face = face_list(i)
          end if
       end do

       if (ele_t < target_distance) then
          neigh_face = face_neigh(xfield, next_face)
          if (neigh_face /= next_face) then

             if (present(ele_dist)) then
                call insert(ele_dist, ele_t - current_t)
             end if

             ! Recurse on the next element
             current_t = ele_t
             new_element = face_ele(xfield, neigh_face)

             ! Record the elements along the path travelled
             ! and the distance travelled within them
             if (present(ele_path)) then
                call insert(ele_path, new_element)
             end if
          else
             if (element_owned(xfield,new_element)) then
                ! Detector is going outside domain
                new_owner=-1
                exit search_loop
             else
                ! The current element is on a Halo, we need to send it to the owner.
                new_owner=element_owner(xfield%mesh,new_element)
                exit search_loop
             end if
          end if
       else
          ! The arrival point is in this element, we're done
          if (present(ele_dist)) then
             call insert(ele_dist, target_distance - current_t)
          end if

          new_owner=getprocno()
          exit search_loop
       end if
       
    end do search_loop
    
  contains

    subroutine ray_intersetion_distance(face, t)
      ! Calculate the distance t of the intersection 
      ! between the face and our ray (r_o, r_d)
      integer, intent(in) :: face
      real, intent(out) :: t

      real :: d, v_n, v_d
      integer, dimension(face_loc(xfield, face)) :: face_nodes
      real, dimension(xfield%dim,xfield%mesh%faces%shape%ngi) :: facet_normals
      real, dimension(xfield%mesh%faces%shape%ngi) :: detwei_f
      real, dimension(xfield%dim) :: face_normal, face_node_val

      ! Get face normal
      call transform_facet_to_physical(xfield, face, detwei_f, facet_normals)
      face_normal = facet_normals(:,1)

      ! Establish d using a point on the plane
      face_nodes = face_global_nodes(xfield, face)
      face_node_val = node_val(xfield,face_nodes(1))
      d = - sum(face_normal * face_node_val)

      ! Calculate t, the distance along the ray to the face intersection
      v_n = dot_product(face_normal, r_o) + d
      v_d = dot_product(face_normal, r_d)
      t = - v_n / v_d

    end subroutine ray_intersetion_distance

  end subroutine geometric_ray_tracing

  subroutine reflect_on_boundary(xfield, coordinate, element)
    ! Reflect the coordinate in the according boundary face of the element.
    ! This uses simple reflection equations and assumes the face to be flat.
    type(vector_field), pointer, intent(in) :: xfield
    real, dimension(mesh_dim(xfield)), intent(inout) :: coordinate
    integer, intent(in) :: element

    real, dimension(mesh_dim(xfield)+1) :: arrival_local_coords
    integer :: i, neigh, neigh_face
    real, dimension(xfield%dim,xfield%mesh%faces%shape%ngi) :: facet_normals
    real, dimension(xfield%mesh%faces%shape%ngi) :: detwei_f
    integer, dimension(:), pointer :: neigh_list
    integer, dimension(:), allocatable :: face_nodes 
    real, dimension(xfield%dim) :: face_node_val, face_normal
    real :: offset, p, D

    arrival_local_coords=local_coords(xfield,element,coordinate)
    neigh = minval(minloc(arrival_local_coords))
    neigh_list=>ele_neigh(xfield,element)

    ! First get the face we went through and the coordinates of one node on it
    neigh_face = ele_face(xfield, element, neigh_list(neigh))
    allocate(face_nodes(face_loc(xfield, neigh_face)))
    face_nodes = face_global_nodes(xfield, neigh_face)    
    face_node_val=node_val(xfield,face_nodes(1))

    ! Now we get the face normal from the transform (detwei_f is a dummy)
    call transform_facet_to_physical(xfield, neigh_face, detwei_f, facet_normals)
    face_normal = facet_normals(:,1)

    ! p = - n . x_f / sqrt(a**2 + b**2 + c**2), where x_f is a point of the face
    ! http://mathworld.wolfram.com/Plane.html, eqs. 3 and 11
    p = 0.0
    do i=1, xfield%dim
        p = p + face_normal(i)**2
    end do
    p = - dot_product(face_node_val,face_normal) / sqrt(p)

    ! D = n . x + p, where x is the point we want to reflect
    ! http://mathworld.wolfram.com/Plane.html, eq. 13
    D = dot_product(face_normal,coordinate) + p

    ! x' = x - 2Dn, where x' is the reflected point
    ! http://mathworld.wolfram.com/Reflection.html, eq. 7
    do i=1, xfield%dim
       coordinate(i) = coordinate(i) - (2 * D * face_normal(i))
    end do

    deallocate(face_nodes)

  end subroutine reflect_on_boundary

  subroutine auto_subcycle_random_walk(detector_list,sub_dt,scale_factor,xfield,diff_field,grad_field,grad2_field)
    type(detector_linked_list), intent(inout) :: detector_list
    real, intent(in) :: sub_dt, scale_factor
    type(scalar_field), pointer, intent(in) :: diff_field
    type(vector_field), pointer, intent(in) :: xfield, grad_field, grad2_field

    type(detector_linked_list) :: subcycle_detector_list
    type(detector_type), pointer :: detector, move_detector
    real, dimension(:), allocatable :: rw_displacement
    integer, dimension(:), allocatable :: auto_subcycle_per_ele
    real :: total_subsubcycle
    integer :: ele, subsubcycle, max_subsubcycle

    max_subsubcycle = 1
    total_subsubcycle = 0
    allocate(rw_displacement(xfield%dim))

    ! Setup our detector work list
    subcycle_detector_list%search_tolerance = detector_list%search_tolerance
    subcycle_detector_list%reflect_on_boundary = detector_list%reflect_on_boundary
    subcycle_detector_list%name = detector_list%name

    ! Create a list of auto-subcycle numbers per element
    allocate(auto_subcycle_per_ele(element_count(xfield)))
    do ele=1, element_count(xfield)
       call element_rw_subcycling(ele, sub_dt, xfield, diff_field, grad_field, &
               grad2_field, scale_factor, auto_subcycle_per_ele(ele))
       end do

    ! First pass establishes how many sub-sub-cycles are required
    detector => detector_list%first
    do while (associated(detector))
       if (detector%type==LAGRANGIAN_DETECTOR) then

          ! Establish how many subcycles are needed
          detector%rw_subsubcycles = auto_subcycle_per_ele(detector%element)
          total_subsubcycle = total_subsubcycle + detector%rw_subsubcycles
          if (detector%rw_subsubcycles > max_subsubcycle) then
             max_subsubcycle = detector%rw_subsubcycles
          end if

          call diffusive_random_walk(detector, sub_dt/detector%rw_subsubcycles, &
                   xfield, diff_field, grad_field, rw_displacement)
          detector%update_vector=detector%update_vector + rw_displacement

          ! Separate the ones that have more to do...
          if (detector%rw_subsubcycles > 1) then
             move_detector=>detector
             detector => detector%next
             call move(move_detector, detector_list, subcycle_detector_list)
          else
             detector => detector%next
          end if
       else
          detector => detector%next
       end if
    end do

    deallocate(auto_subcycle_per_ele)

    ! Update elements in the subcycle list
    if (subcycle_detector_list%length > 0) then
       call track_detectors(subcycle_detector_list,xfield)
    end if

    ! Remaining subcycle passes until subcycle list is empty
    subsubcycle = 2
    do while (subcycle_detector_list%length > 0)
       detector => subcycle_detector_list%first
       do while (associated(detector))
          call diffusive_random_walk(detector, sub_dt/detector%rw_subsubcycles, &
                   xfield, diff_field, grad_field, rw_displacement)
          detector%update_vector=detector%update_vector + rw_displacement

          ! Separate the ones that have more to do...
          if (detector%rw_subsubcycles <= subsubcycle) then
             move_detector=>detector
             detector => detector%next
             call move(move_detector, subcycle_detector_list, detector_list)
          else
             detector => detector%next
          end if
       end do
       call track_detectors(subcycle_detector_list,xfield)
       subsubcycle = subsubcycle + 1
    end do

    ! Update detectors after the final move
    call track_detectors(detector_list,xfield)

    ewrite(2,*) "Auto-subcycling: max. subcycles:", max_subsubcycle
    ewrite(2,*) "Auto-subcycling: avg. subcycles:", total_subsubcycle / detector_list%length

    deallocate(rw_displacement)

  end subroutine auto_subcycle_random_walk

  subroutine element_rw_subcycling(element,dt,xfield,diff_field,grad_field,grad2_field,scale_factor,subcycles)
    integer, intent(in) :: element
    type(scalar_field), pointer, intent(in) :: diff_field
    type(vector_field), pointer, intent(in) :: xfield, grad_field, grad2_field
    real, intent(in) :: dt, scale_factor
    integer, intent(out) :: subcycles

    type(ilist) :: neigh_ele_list, ele_in_range
    type(integer_hash_table) :: visited_eles
    type(inode), pointer :: ele
    integer :: i, current_ele, local_subcycling, key, value
    real :: k0, d_z, min_z, max_z, ele_min_z, ele_max_z, neigh_min_z, neigh_max_z, min_dt
    integer, dimension(:), pointer :: neighbours
    real, dimension(xfield%mesh%shape%loc) :: coords
    real, dimension(grad_field%mesh%shape%loc) :: k_grad
    real, dimension(grad2_field%mesh%shape%loc) :: k_grad2

    ! Calculate the vertical range for the criteria
    k0 = maxval(abs(ele_val(diff_field, element)))
    k_grad = ele_val(grad_field, grad_field%dim, element)
    d_z = sqrt(6 * K0 * dt) + maxval(abs(k_grad)) * dt

    coords = ele_val(xfield, xfield%dim, element)
    min_z = minval(coords) - d_z
    max_z = maxval(coords) + d_z

    call insert(neigh_ele_list, element)
    call allocate(visited_eles)

    ! Now unfold neighbours into neigh_ele_list until we are out of range
    do while (neigh_ele_list%length >= 1)
       current_ele = ipop(neigh_ele_list)

       if (has_key(visited_eles, current_ele)) cycle

       call insert(visited_eles, current_ele, 0)

       coords = ele_val(xfield, xfield%dim, current_ele)
       ele_min_z = minval(coords)
       ele_max_z = maxval(coords)

       ! if one of the two z coordinates is in the range add the neighbours to our search space
       if ((ele_min_z>min_z .and. ele_min_z<max_z).or.(ele_max_z<max_z .and. ele_max_z>min_z)) then
          call insert(ele_in_range, current_ele)

          ! add neighbours to search space
          neighbours => ele_neigh(xfield, current_ele)
          do i=1, size(neighbours)
             if (neighbours(i)>0) then
                call insert(neigh_ele_list, neighbours(i))
             end if
          end do
       end if       
    end do

    ! Now we have all elements in our range and we can find the maximum sub-cycle number
    subcycles = 1
    ele => ele_in_range%firstnode
    do while(associated(ele))
       k_grad2 = ele_val(grad2_field, grad2_field%dim, ele%value)
       k_grad2 = 1.0 / abs(k_grad2)
       min_dt = minval(k_grad2) / scale_factor
       local_subcycling = ceiling(dt/min_dt)
       if (local_subcycling > subcycles) then
          subcycles = local_subcycling
       end if
       ele => ele%next
    end do
    call flush_list(ele_in_range)
    call deallocate(visited_eles)

  end subroutine element_rw_subcycling

  subroutine diffusive_random_walk(detector, dt, xfield, diff_field, grad_field, displacement)
    type(detector_type), pointer, intent(in) :: detector
    type(scalar_field), pointer, intent(in) :: diff_field
    type(vector_field), pointer, intent(in) :: xfield, grad_field
    real, intent(in) :: dt
    real, dimension(xfield%dim), intent(out) :: displacement

    real, dimension(1) :: rnd
    real :: K
    real, dimension(xfield%dim) :: position, K_grad
    real, dimension(xfield%dim+1) :: lcoord
    integer :: new_owner, offset_element

    call random_number(rnd)
    rnd = (rnd * 2.0) - 1.0
    K_grad=detector_value(grad_field, detector)

    position = detector_value(xfield, detector)
    position(xfield%dim)=position(xfield%dim) + 0.5*dt*K_grad(xfield%dim)
    offset_element=detector%element
    call local_guided_search(xfield,position,offset_element,1e-10,new_owner,lcoord)

    ! In case the offset sampling point is outside the domain we extrapolate
    if (new_owner < 0) then
       K=eval_field(detector%element, diff_field, detector%local_coords)
    end if

    if (new_owner/=getprocno() .and. new_owner > 0) then
       ! Offset sampling point is on another parallel domain
       ewrite(-1,*) "Detected non-local element in internal Diffusive Random Walk;"
       ewrite(-1,*) "offset_element", offset_element, "detector%element", detector%element
       FLAbort("Guided search in internal Random Walk function detected non-local element")
    else
       ! Evaluate K at the offset sampling point
       K=eval_field(offset_element, diff_field, lcoord)
    end if

    ! Make sure we don't use negative K values to prevent NaN in the sqrt()
    K = abs(K)

    displacement(:)=0.0
    displacement(xfield%dim)=K_grad(xfield%dim)*dt + rnd(1)*sqrt(6*K*dt)

  end subroutine diffusive_random_walk

  subroutine naive_random_walk(detector, dt, xfield, diff_field, displacement)
    type(detector_type), pointer, intent(in) :: detector
    type(scalar_field), pointer, intent(in) :: diff_field
    type(vector_field), pointer, intent(in) :: xfield
    real, intent(in) :: dt
    real, dimension(xfield%dim), intent(out) :: displacement

    real, dimension(1) :: rnd
    real :: K

    call random_number(rnd)
    rnd(1) = (rnd(1) * 2.0) - 1.0
    K=abs(detector_value(diff_field, detector))

    displacement(:)=0.0
    displacement(xfield%dim)=rnd(1)*sqrt(6*K*dt)

  end subroutine naive_random_walk


  function check_any_lagrangian(detector_list)
    ! Check if there are any lagrangian detectors in the given list
    ! across all processors
    logical :: check_any_lagrangian
    type(detector_linked_list), intent(inout) :: detector_list

    type(detector_type), pointer :: detector
    integer :: checkint 
      
    checkint = 0
    detector => detector_list%first
    do while (associated(detector))       
       if (detector%type==LAGRANGIAN_DETECTOR) then
          checkint = 1
          exit
       end if
       detector => detector%next
    end do
    call allmax(checkint)
    check_any_lagrangian = .false.
    if(checkint>0) check_any_lagrangian = .true.

  end function check_any_lagrangian

end module detector_move_lagrangian
