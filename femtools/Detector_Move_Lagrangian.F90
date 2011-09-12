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
  use python_state
  use iso_c_binding
  use transform_elements

  implicit none
  
  private

  public :: move_lagrangian_detectors, read_detector_move_options, check_any_lagrangian

  character(len=OPTION_PATH_LEN), parameter :: rk_gs_path="/lagrangian_timestepping/explicit_runge_kutta_guided_search"

contains

  subroutine read_detector_move_options(detector_list, detector_path)
    ! Subroutine to allocate the detector parameters, 
    ! including RK stages and update vector
    type(detector_linked_list), intent(inout) :: detector_list
    character(len=*), intent(in) :: detector_path

    type(rk_gs_parameters), pointer :: parameters
    integer :: i,j,k
    real, allocatable, dimension(:) :: stage_weights
    integer, dimension(2) :: option_rank

    if (associated(detector_list%move_parameters)) then
       deallocate(detector_list%move_parameters)
    end if
    allocate(detector_list%move_parameters)
    parameters => detector_list%move_parameters

    if (have_option(trim(detector_path)//trim("/lagrangian_timestepping/reflect_on_boundary"))) then
       parameters%reflect_on_boundary=.true.
    end if

    if(have_option(trim(detector_path)//"/lagrangian_timestepping")) then

       call get_option(trim(detector_path)//"/lagrangian_timestepping/subcycles",parameters%n_subcycles)
       call get_option(trim(detector_path)//"/lagrangian_timestepping/search_tolerance",parameters%search_tolerance)

       ! Forward Euler options
       if (have_option(trim(detector_path)//"/lagrangian_timestepping/forward_euler_guided_search")) then
          parameters%n_stages = 1
          allocate(parameters%timestep_weights(parameters%n_stages))
          parameters%timestep_weights = 1.0
       end if

       ! Parameters for classical Runge-Kutta
       if (have_option(trim(detector_path)//"/lagrangian_timestepping/rk4_guided_search")) then
          parameters%n_stages = 4
          allocate(stage_weights(parameters%n_stages*(parameters%n_stages-1)/2))
          stage_weights = (/0.5, 0., 0.5, 0., 0., 1./)
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
          allocate(parameters%timestep_weights(parameters%n_stages))
          parameters%timestep_weights = (/ 1./6., 1./3., 1./3., 1./6. /)
       end if

       ! Generic Runge-Kutta options
       if (have_option(trim(detector_path)//trim(rk_gs_path))) then
          call get_option(trim(detector_path)//trim(rk_gs_path)//"/n_stages",parameters%n_stages)

          ! Allocate and read stage_matrix from options
          allocate(stage_weights(parameters%n_stages*(parameters%n_stages-1)/2))
          option_rank = option_shape(trim(detector_path)//trim(rk_gs_path)//"/stage_weights")
          if (option_rank(2).ne.-1) then
             FLExit('Stage Array wrong rank')
          end if
          if (option_rank(1).ne.size(stage_weights)) then
             ewrite(-1,*) 'size expected was', size(stage_weights)
             ewrite(-1,*) 'size actually was', option_rank(1)
             FLExit('Stage Array wrong size')
          end if
          call get_option(trim(detector_path)//trim(rk_gs_path)//"/stage_weights",stage_weights)
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
          option_rank = option_shape(trim(detector_path)//trim(rk_gs_path)//"/timestep_weights")
          if (option_rank(2).ne.-1) then
             FLExit('Timestep Array wrong rank')
          end if
          if (option_rank(1).ne.size(parameters%timestep_weights)) then
             FLExit('Timestep Array wrong size')
          end if
          call get_option(trim(detector_path)//trim(rk_gs_path)//"/timestep_weights",parameters%timestep_weights)
       end if

    else
       if (check_any_lagrangian(detector_list)) then
          ewrite(-1,*) "Found lagrangian detectors, but no timstepping options"
          FLExit('No lagrangian timestepping specified')
       end if
    end if

  end subroutine read_detector_move_options

  subroutine move_lagrangian_detectors(state, detector_list, dt, timestep)
    type(state_type), dimension(:), intent(in) :: state
    type(detector_linked_list), intent(inout) :: detector_list
    real, intent(in) :: dt
    integer, intent(in) :: timestep

    type(rk_gs_parameters), pointer :: parameters
    type(vector_field), pointer :: vfield, xfield
    type(detector_type), pointer :: detector
    type(detector_linked_list), dimension(:), allocatable :: send_list_array
    type(halo_type), pointer :: ele_halo
    integer :: i, j, k, num_proc, dim, all_send_lists_empty, nprocs, stage, cycle, det
    logical :: any_lagrangian
    real :: rk_dt

    ! Random Walk velocity source
    real, dimension(:), allocatable :: rw_velocity_source
    real, dimension(:,:), allocatable :: rw_displacement
    type(scalar_field), pointer :: diffusivity_field
    type(vector_field), pointer :: diffusivity_grad

    ewrite(1,*) "In move_lagrangian_detectors for detectors list: ", detector_list%name
    ewrite(2,*) "Detector list", detector_list%id, "has", detector_list%length, &
         "local and", detector_list%total_num_det, "global detectors"

    parameters => detector_list%move_parameters

    ! For Random Walk first run the user code, so we can pull fields from state
    if (parameters%do_random_walk.and. .not.parameters%use_internal_rw) then
       ! Run the user's code and store val object in "random_walk" dict
       call python_run_detector_string(trim(parameters%rw_pycode),&
              trim(detector_list%name),trim("random_walk"))
    end if

    ! Pull some information from state
    xfield=>extract_vector_field(state(1), "Coordinate")
    vfield=>extract_vector_field(state(1),"Velocity")
    allocate(rw_velocity_source(xfield%dim))

    if (parameters%use_internal_rw) then
       diffusivity_field=>extract_scalar_field(state(1), trim(parameters%diffusivity_field))
       diffusivity_grad=>extract_vector_field(state(1), trim(parameters%diffusivity_grad))
    end if

    ! We allocate a sendlist for every processor
    nprocs=getnprocs()
    allocate(send_list_array(nprocs))

    ! Allocate det%k and det%update_vector
    call allocate_rk_guided_search(detector_list, xfield%dim, parameters%n_stages)
    rk_dt = dt/parameters%n_subcycles

    subcycling_loop: do cycle = 1, parameters%n_subcycles

       ! Reset update_vector to position
       detector => detector_list%first
       do det = 1, detector_list%length
          if(detector%type==LAGRANGIAN_DETECTOR) then
             detector%update_vector = detector%position
          end if
          detector => detector%next
       end do

       RKstages_loop: do stage = 1, parameters%n_stages

          ! Compute the update vector
          call set_stage(detector_list,vfield,xfield,rk_dt,stage)

          ! This loop continues until all detectors have completed their
          ! timestep this is measured by checking if the send and receive
          ! lists are empty in all processors
          detector_timestepping_loop: do  

             ! Make sure we still have lagrangian detectors
             any_lagrangian=check_any_lagrangian(detector_list)
             if (any_lagrangian) then

                !Detectors leaving the domain from non-owned elements
                !are entering a domain on another processor rather 
                !than leaving the physical domain. In this subroutine
                !such detectors are removed from the detector list
                !and added to the send_list_array
                call move_detectors_guided_search(detector_list,&
                        vfield,xfield,send_list_array,parameters%search_tolerance)

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
                call exchange_detectors(state(1),detector_list, send_list_array)
             else
                ! If we run out of lagrangian detectors for some reason, exit the loop
                exit
             end if
          end do detector_timestepping_loop
       end do RKstages_loop

       ! Add the Random Walk displacement 
       if (parameters%do_random_walk) then
          allocate(rw_displacement(xfield%dim,detector_list%length))

          detector => detector_list%first
          do det = 1, detector_list%length
             if (detector%type==LAGRANGIAN_DETECTOR) then
                ! Evaluate the RW python function and add to update_vector
                if (parameters%use_internal_rw) then
                   call calc_diffusive_rw(detector, rk_dt, xfield, diffusivity_field, &
                          diffusivity_grad, rw_velocity_source)
                else
                   call python_run_detector_val_function(detector,xfield,rk_dt, &
                          trim(detector_list%name),trim("random_walk"),rw_velocity_source)
                end if
                detector%update_vector=detector%update_vector + (rw_velocity_source)  
                detector%search_complete=.false.
             end if
             detector => detector%next
          end do

          !call python_evaluate_random_walk(detector_list, xfield, rk_dt, rw_displacement)

          !detector => detector_list%first
          !do det = 1, detector_list%length
          !   detector%update_vector=detector%update_vector + rw_displacement(:,det)
          !   detector%search_complete=.false.

          !   detector => detector%next
          !end do

          call move_detectors_guided_search(detector_list,vfield,xfield, &
                 send_list_array,parameters%search_tolerance)

          ! If needed exchange_detectors
          all_send_lists_empty=0
          do k=1, nprocs
             if (send_list_array(k)%length/=0) then
                all_send_lists_empty=1
             end if
          end do
          call allmax(all_send_lists_empty)
          if (all_send_lists_empty>0) then
             call exchange_detectors(state(1),detector_list, send_list_array)
          end if

          deallocate(rw_displacement)
       end if

       ! After everything is done, we update detector%position
       detector => detector_list%first
       do det = 1, detector_list%length
          if (detector%type==LAGRANGIAN_DETECTOR) then
             detector%position = detector%update_vector
          end if
          detector => detector%next
       end do

    end do subcycling_loop

    deallocate(send_list_array)
    deallocate(rw_velocity_source)

    ! Make sure all local detectors are owned and distribute the ones that 
    ! stoppped moving in a halo element
    call distribute_detectors(state, detector_list)

    ! This needs to be called after distribute_detectors because the exchange  
    ! routine serialises det%k and det%update_vector if it finds the RK-GS option
    call deallocate_rk_guided_search(detector_list)

    ewrite(2,*) "After moving and distributing we have", detector_list%length, &
         "local and", detector_list%total_num_det, "global detectors"
    ewrite(1,*) "Exiting move_lagrangian_detectors"

  end subroutine move_lagrangian_detectors

  function check_any_lagrangian(detector_list0)
    ! Check if there are any lagrangian detectors in the given list
    ! across all processors
    logical :: check_any_lagrangian
    type(detector_linked_list), intent(inout) :: detector_list0
    type(detector_type), pointer :: det0
    integer :: i
    integer :: checkint 
      
    checkint = 0
    det0 => detector_list0%first
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

  subroutine allocate_rk_guided_search(detector_list, dim, n_stages)
    ! Allocate the RK stages and update vector
    type(detector_linked_list), intent(inout) :: detector_list
    integer, intent(in) :: n_stages, dim
    type(detector_type), pointer :: det0

    det0 => detector_list%first
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

  subroutine deallocate_rk_guided_search(detector_list)
    ! Deallocate the RK stages and update vector
    type(detector_linked_list), intent(inout) :: detector_list
      
    type(detector_type), pointer :: det0
    integer :: j0
      
    det0 => detector_list%first
    do j0=1, detector_list%length
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

  subroutine set_stage(detector_list,vfield,xfield,dt0,stage0)
    ! Compute the vector to search for in the next RK stage
    ! If this is the last stage, update detector position
    type(detector_linked_list), intent(inout) :: detector_list
    type(vector_field), pointer, intent(in) :: vfield, xfield
    real, intent(in) :: dt0
    integer, intent(in) :: stage0
    
    type(rk_gs_parameters), pointer :: parameters
    type(detector_type), pointer :: det0
    integer :: det_count,j0
    real, dimension(mesh_dim(xfield)+1) :: stage_local_coords

    parameters => detector_list%move_parameters
    
    det0 => detector_list%first
    do det_count=1, detector_list%length 

       if(det0%type==LAGRANGIAN_DETECTOR) then

          !! Set det%update_vector for lagrangian advection
          det0%search_complete = .false.

          ! Evaluate velocity at update_vector and set k
          if (parameters%do_velocity_advect) then
             ! stage vector is computed by evaluating velocity at current position
             stage_local_coords=local_coords(xfield,det0%element,det0%update_vector)
             det0%k(stage0,:)=eval_field(det0%element, vfield, stage_local_coords)
          else
             ! do not advect detector with the velocity field
             det0%k(stage0,:)=0.0
          end if

          if(stage0<parameters%n_stages) then
             ! Update vector maps from current position to place required
             ! for computing next stage vector
             det0%update_vector = det0%position
             do j0 = 1, stage0
                det0%update_vector = det0%update_vector + dt0*parameters%stage_matrix(stage0+1,j0)*det0%k(j0,:)
             end do
          else
             ! Update vector maps from current position to final position
             det0%update_vector = det0%position
             do j0 = 1, parameters%n_stages
                det0%update_vector = det0%update_vector + dt0*parameters%timestep_weights(j0)*det0%k(j0,:)
             end do
          end if

       end if
       det0 => det0%next
    end do
  end subroutine set_stage


  subroutine move_detectors_guided_search(detector_list,vfield,xfield,send_list_array,search_tolerance)
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
    type(detector_linked_list), dimension(:), intent(inout) :: send_list_array
    type(vector_field), pointer, intent(in) :: vfield,xfield
    real, intent(in) :: search_tolerance

    type(detector_type), pointer :: det0, det_send
    integer :: det_count
    logical :: owned
    real, dimension(mesh_dim(vfield)+1) :: arrival_local_coords
    integer, dimension(:), pointer :: neigh_list
    integer :: neigh, proc_local_number, face, i, j
    logical :: outside_domain

    ! boundary reflection
    integer :: neigh_face
    real, dimension(xfield%dim,xfield%mesh%faces%shape%ngi) :: facet_normals
    real, dimension(xfield%mesh%faces%shape%ngi) :: detwei_f
    integer, dimension(:), allocatable :: face_nodes 
    real, dimension(xfield%dim) :: face_node_val, face_normal
    real :: offset, p, D

    ewrite(2,*) "In move_detectors_guided_search"

    !Loop over all the detectors
    det0 => detector_list%first
    do det_count=1, detector_list%length

       !Only move Lagrangian detectors
       if (det0%type==LAGRANGIAN_DETECTOR.and..not.det0%search_complete) then
          search_loop: do

             !Compute the local coordinates of the arrival point with respect to this element
             arrival_local_coords=local_coords(xfield,det0%element,det0%update_vector)
             if (minval(arrival_local_coords)>-search_tolerance) then
                !the arrival point is in this element
                det0%search_complete = .true.
                !move on to the next detector
                det0 => det0%next
                exit search_loop
             end if

             !The arrival point is not in this element, try to get closer to it by 
             !searching in the coordinate direction in which it is furthest away
             neigh = minval(minloc(arrival_local_coords))
             neigh_list=>ele_neigh(xfield,det0%element)
             if (neigh_list(neigh)>0) then
                !the neighbouring element is also on this domain
                !so update the element and try again
                det0%element = neigh_list(neigh)
             else
                !check if this element is owned (to decide where
                !to send particles leaving the processor domain)
                if (element_owned(vfield,det0%element)) then
                   !this face goes outside of the computational domain
                   !try all of the faces with negative local coordinate
                   !just in case we went through a corner
                   outside_domain=.true.
                   face_search: do face = 1, size(arrival_local_coords)
                      if (arrival_local_coords(face)<-search_tolerance.and.neigh_list(face)>0) then
                         outside_domain = .false.
                         det0%element = neigh_list(face)
                         exit face_search
                      end if
                   end do face_search

                   if (outside_domain) then
                      if (detector_list%move_parameters%reflect_on_boundary) then
                         ! We reflect the detector path at the face we just went through

                         ! First get the face we went through and the coordinates of one node on it
                         neigh_face = ele_face(xfield, det0%element, neigh_list(neigh))
                         allocate(face_nodes(face_loc(xfield, neigh_face)))
                         face_nodes = face_global_nodes(xfield, neigh_face)    
                         face_node_val=node_val(xfield,face_nodes(1))

                         ! Now we get the face normal from the transform (detwei_f is a dummy)
                         call transform_facet_to_physical(xfield, neigh_face, detwei_f, facet_normals)
                         face_normal = facet_normals(:,1)

                         ! p = - n . x_f / sqrt(a**2 + b**2 + c**2), where x_f is a point of the face
                         ! http://mathworld.wolfram.com/Plane.html, eqs. 3 and 11
                         p = 0.0
                         do j=1, xfield%dim
                            p = p + face_normal(j)**2
                         end do
                         p = - dot_product(face_node_val,face_normal) / sqrt(p)

                         ! D = n . x + p, where x is the point we want to reflect
                         ! http://mathworld.wolfram.com/Plane.html, eq. 13
                         D = dot_product(face_normal,det0%update_vector) + p

                         ! x' = x - 2Dn, where x' is the reflected point
                         ! http://mathworld.wolfram.com/Reflection.html, eq. 7
                         do j=1, xfield%dim
                            det0%update_vector(j) = det0%update_vector(j) - (2 * D * face_normal(j))
                         end do

                         deallocate(face_nodes)
                      else
                         ! Turn detector static inside the domain
                         ewrite(1,*) "WARNING: detector attempted to leave computational &
                           domain; making it static, detector ID:", det0%id_number, "detector element:", det0%element
                         det0%type=STATIC_DETECTOR
                         ! move on to the next detector, without updating det0%position, 
                         ! because det0%update_vector is by now outside of the computational domain
                         det0 => det0%next
                         exit search_loop
                      end if
                   end if
                else
                   det_send => det0
                   det0 => det0%next
                   !this face goes into another computational domain
                   proc_local_number=element_owner(vfield%mesh,det_send%element)

                   call move(det_send, detector_list, send_list_array(proc_local_number))
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

  subroutine calc_diffusive_rw(detector, dt, xfield, diff_field, grad_field, displacement)
    type(detector_type), pointer, intent(in) :: detector
    type(scalar_field), pointer, intent(in) :: diff_field
    type(vector_field), pointer, intent(in) :: xfield, grad_field
    real, intent(in) :: dt
    real, dimension(xfield%dim), intent(out) :: displacement

    real, dimension(1) :: rnd
    real :: K
    real, dimension(xfield%dim) :: position, K_grad
    real, dimension(xfield%dim+1) :: lcoord

    call random_number(rnd)
    rnd = (rnd * 2.0) - 1.0
    K_grad=eval_field(detector%element, grad_field, detector%local_coords)

    position(:)=eval_field(detector%element, xfield, detector%local_coords)
    position(xfield%dim)=position(xfield%dim) + 0.5*dt*K_grad(xfield%dim)
    lcoord=local_coords(xfield, detector%element, position)
    K=eval_field(detector%element, diff_field, lcoord)

    displacement(:)=0.0
    displacement(xfield%dim)=K_grad(xfield%dim)*dt + rnd(1)*sqrt(6*K*dt)

  end subroutine calc_diffusive_rw

end module detector_move_lagrangian
