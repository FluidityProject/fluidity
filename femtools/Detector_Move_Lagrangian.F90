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
    type(detector_type), pointer :: detector, move_detector
    type(detector_linked_list) :: subcycle_detector_list
    type(halo_type), pointer :: ele_halo
    integer :: i, j, num_proc, dim, stage, cycle, subsubcycle
    logical :: any_lagrangian
    real :: sub_dt

    ! Random Walk velocity source
    real, dimension(:), allocatable :: rw_velocity_source
    type(scalar_field), pointer :: diffusivity_field
    type(vector_field), pointer :: diffusivity_grad, diffusivity_2nd_grad
    integer :: auto_subcycles = 0

    ewrite(1,*) "In move_lagrangian_detectors for detectors list: ", detector_list%name
    ewrite(2,*) "Detector list", detector_list%id, "has", detector_list%length, &
         "local and", detector_list%total_num_det, "global detectors"

    parameters => detector_list%move_parameters
    subcycle_detector_list%move_parameters => parameters
    subcycle_detector_list%move_parameters => parameters

    ! For Random Walk first run the user code, so we can pull fields from state
    if (parameters%do_random_walk.and. .not.parameters%use_internal_rw) then
       ! Run the user's code and store val object in "random_walk" dict
       call python_run_detector_string(trim(parameters%rw_pycode),&
              trim(detector_list%name),trim("random_walk"))
    end if

    ! Pull some information from state
    xfield=>extract_vector_field(state(1), "Coordinate")
    vfield=>extract_vector_field(state(1), "Velocity")
    allocate(rw_velocity_source(xfield%dim))

    if (parameters%use_internal_rw) then
       diffusivity_field=>extract_scalar_field(state(1), trim(parameters%diffusivity_field))
       diffusivity_grad=>extract_vector_field(state(1), trim(parameters%diffusivity_grad))

       if (parameters%auto_subcycle) then
          diffusivity_2nd_grad=>extract_vector_field(state(1), trim(parameters%diffusivity_2nd_grad))
       end if
    end if

    ! Allocate det%k and det%update_vector
    call allocate_rk_guided_search(detector_list, xfield%dim, parameters%n_stages)
    sub_dt = dt/parameters%n_subcycles

    ! This is the outer, user-defined subcycling loop
    subcycling_loop: do cycle = 1, parameters%n_subcycles

       ! Reset update_vector to position
       detector => detector_list%first
       do while (associated(detector))
          if(detector%type==LAGRANGIAN_DETECTOR) then
             detector%update_vector = detector%position
          end if
          detector => detector%next
       end do

       ! Explicit Runge-Kutta iterations
       RKstages_loop: do stage = 1, parameters%n_stages

          ! Compute the update vector for the current stage
          call set_stage(detector_list,vfield,xfield,sub_dt,stage)

          ! Move update parametric detector coordinates according to update_vector
          ! If this takes a detector across parallel domain boundaries
          ! the routine will also send the detector
          call move_detectors_guided_search(detector_list,xfield)

       end do RKstages_loop

       ! Add the Random Walk displacement 
       if (parameters%do_random_walk) then

          if (parameters%use_internal_rw .and. parameters%auto_subcycle) then
             ! Sub-cycle the RW 

             ! First pass establishes how many sub-sub-cycles are required
             detector => detector_list%first
             do while (associated(detector))
                if (detector%type==LAGRANGIAN_DETECTOR) then
                   call calc_auto_subcycling_rw(detector, sub_dt, xfield, diffusivity_field, &
                          diffusivity_2nd_grad, detector%rw_subsubcycles)         

                   call calc_diffusive_rw(detector, sub_dt/detector%rw_subsubcycles, xfield, &
                          diffusivity_field, diffusivity_grad, rw_velocity_source)
                   detector%update_vector=detector%update_vector + (rw_velocity_source)

                   ! Separate the ones that have more to do...
                   if (detector%rw_subsubcycles > 1) then
                      move_detector=>detector
                      detector => detector%next
                      call move(move_detector, detector_list, subcycle_detector_list)
                   else
                      detector => detector%next
                   end if
                end if
             end do

             ! Update elements in the subcycle list
             if (subcycle_detector_list%length > 0) then
                call move_detectors_guided_search(subcycle_detector_list,xfield)
             end if

             subsubcycle = 2
             do while (subcycle_detector_list%length > 0)
                ewrite(2,*) "ml805 subsubcycle", subsubcycle
                detector => subcycle_detector_list%first
                do while (associated(detector))
                   call calc_diffusive_rw(detector, sub_dt/detector%rw_subsubcycles, xfield, diffusivity_field, &
                          diffusivity_grad, rw_velocity_source)
                   detector%update_vector=detector%update_vector + (rw_velocity_source)

                   ! Separate the ones that have more to do...
                   if (detector%rw_subsubcycles <= subsubcycle) then
                      move_detector=>detector
                      detector => detector%next
                      call move(move_detector, subcycle_detector_list, detector_list)
                   else
                      detector => detector%next
                   end if
                end do
                call move_detectors_guided_search(subcycle_detector_list,xfield)
                subsubcycle = subsubcycle + 1
             end do

             call move_detectors_guided_search(detector_list,xfield)

          else
             ! Evaluate RW for one cycle
             detector => detector_list%first
             do while (associated(detector))
                if (detector%type==LAGRANGIAN_DETECTOR) then
                   ! Evaluate the RW python function and add to update_vector
                   if (parameters%use_internal_rw) then
                      call calc_diffusive_rw(detector, sub_dt, xfield, diffusivity_field, &
                          diffusivity_grad, rw_velocity_source)
                   else
                      call python_run_detector_val_function(detector,xfield,sub_dt, &
                          trim(detector_list%name),trim("random_walk"),rw_velocity_source)
                   end if
                   detector%update_vector=detector%update_vector + (rw_velocity_source)
                end if
                detector => detector%next
             end do

             call move_detectors_guided_search(detector_list,xfield)
          end if
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

    deallocate(rw_velocity_source)

    ! Make sure all local detectors are owned and distribute the ones that 
    ! stoppped moving in a halo element
    call distribute_detectors(state(1), detector_list)

    ! This needs to be called after distribute_detectors because the exchange  
    ! routine serialises det%k and det%update_vector if it finds the RK-GS option
    call deallocate_rk_guided_search(detector_list)

    ewrite(2,*) "After moving and distributing we have", detector_list%length, &
         "local and", detector_list%total_num_det, "global detectors"
    ewrite(1,*) "Exiting move_lagrangian_detectors"

  end subroutine move_lagrangian_detectors

  subroutine set_stage(detector_list,vfield,xfield,dt,stage)
    ! Compute the vector to search for in the next RK stage
    ! If this is the last stage, update detector position
    type(detector_linked_list), intent(inout) :: detector_list
    type(vector_field), pointer, intent(in) :: vfield, xfield
    real, intent(in) :: dt
    integer, intent(in) :: stage
    
    type(rk_gs_parameters), pointer :: parameters
    type(detector_type), pointer :: detector
    integer :: j0
    real, dimension(mesh_dim(xfield)+1) :: stage_local_coords

    parameters => detector_list%move_parameters
    
    detector => detector_list%first
    do while (associated(detector))

       if(detector%type==LAGRANGIAN_DETECTOR) then

          ! Evaluate velocity at update_vector and set k
          if (parameters%do_velocity_advect) then
             ! stage vector is computed by evaluating velocity at current position
             stage_local_coords=local_coords(xfield,detector%element,detector%update_vector)
             detector%k(stage,:)=eval_field(detector%element, vfield, stage_local_coords)
          else
             ! do not advect detector with the velocity field
             detector%k(stage,:)=0.0
          end if

          if(stage<parameters%n_stages) then
             ! Update vector maps from current position to place required
             ! for computing next stage vector
             detector%update_vector = detector%position
             do j0 = 1, stage
                detector%update_vector = detector%update_vector + dt*parameters%stage_matrix(stage+1,j0)*detector%k(j0,:)
             end do
          else
             ! Update vector maps from current position to final position
             detector%update_vector = detector%position
             do j0 = 1, parameters%n_stages
                detector%update_vector = detector%update_vector + dt*parameters%timestep_weights(j0)*detector%k(j0,:)
             end do
          end if

       end if
       detector => detector%next
    end do
  end subroutine set_stage

  subroutine move_detectors_guided_search(detector_list,xfield)
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
    logical :: outside_domain, any_lagrangian

    ! We allocate a sendlist for every processor
    nprocs=getnprocs()
    allocate(send_list_array(nprocs))

    search_tolerance=detector_list%move_parameters%search_tolerance

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

             call local_guided_search(xfield,detector%update_vector,detector%element, &
                        search_tolerance,new_owner,detector%local_coords)

             if (new_owner==-1) then
                if (detector_list%move_parameters%reflect_on_boundary) then
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
          call exchange_detectors(detector_list,xfield,send_list_array)
       else
          ! If we run out of lagrangian detectors for some reason, exit the loop
          exit
       end if
    end do detector_timestepping_loop

    deallocate(send_list_array)

  end subroutine move_detectors_guided_search

  subroutine local_guided_search(xfield,coordinate,element,search_tolerance,new_owner)
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

    integer, dimension(:), pointer :: neigh_list
    integer :: neigh, face
    logical :: outside_domain

    search_loop: do
       ! Compute the local coordinates of the arrival point with respect to this element
       l_coords=local_coords(xfield,element,coordinate)
       if (minval(l_coords)>-search_tolerance) then
          !The arrival point is in this element, we're done
          new_owner=getprocno()
          return
       end if

       ! The arrival point is not in this element, try to get closer to it by 
       ! searching in the coordinate direction in which it is furthest away
       neigh = minval(minloc(l_coords))
       neigh_list=>ele_neigh(xfield,element)
       if (neigh_list(neigh)>0) then
          ! The neighbouring element is also on this domain
          ! so update the element and try again
          element = neigh_list(neigh)
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
                return
             end if
          else
             ! The current element is on a Halo, we need to send it to the owner.
             new_owner=element_owner(xfield%mesh,element)
             return
          end if
       end if
    end do search_loop

  end subroutine local_guided_search

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

  subroutine calc_auto_subcycling_rw(detector, dt, xfield, diff_field, grad2_field, subcycles)
    type(detector_type), pointer, intent(in) :: detector
    type(scalar_field), pointer, intent(in) :: diff_field
    type(vector_field), pointer, intent(in) :: xfield, grad2_field
    real, intent(in) :: dt
    integer, intent(out) :: subcycles

    integer, dimension(:), pointer :: current_ele_nodes
    real, dimension(:), allocatable :: node_vals

    ! Find maximum K'' node value in current element
    ! dt << MIN(1/|K''|), for all nodes
    current_ele_nodes=>ele_nodes(grad2_field, detector%element)
    allocate(node_vals(size(current_ele_nodes)))
    node_vals = abs(node_val(grad2_field, grad2_field%dim, current_ele_nodes))
    node_vals = 1.0 / node_vals
    subcycles = ceiling(dt/minval(node_vals))
    deallocate(node_vals)

  end subroutine calc_auto_subcycling_rw

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
    integer :: new_owner, offset_element

    call random_number(rnd)
    rnd = (rnd * 2.0) - 1.0
    K_grad=detector_value(grad_field, detector)

    position = detector_value(xfield, detector)
    position(xfield%dim)=position(xfield%dim) + 0.5*dt*K_grad(xfield%dim)
    offset_element=detector%element
    call local_guided_search(xfield,position,offset_element,1e-10,new_owner,lcoord)
    if (new_owner/=getprocno()) then
       ewrite(-1,*) "Detected non-local element in internal Diffusive Random Walk;"
       ewrite(-1,*) "offset_element", offset_element, "detector%element", detector%element
       FLAbort("Guided search in internal Random Walk function detected non-local element")
    end if
    K=eval_field(offset_element, diff_field, lcoord)

    ! Make sure we don't use negative K values to prevent NaN in the sqrt()
    K = abs(K)

    displacement(:)=0.0
    displacement(xfield%dim)=K_grad(xfield%dim)*dt + rnd(1)*sqrt(6*K*dt)

  end subroutine calc_diffusive_rw

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

  subroutine allocate_rk_guided_search(detector_list, dim, n_stages)
    ! Allocate the RK stages and update vector
    type(detector_linked_list), intent(inout) :: detector_list
    integer, intent(in) :: n_stages, dim

    type(detector_type), pointer :: detector

    detector => detector_list%first
    do while (associated(detector))
       if(detector%type==LAGRANGIAN_DETECTOR) then
          if(allocated(detector%k)) then
             deallocate(detector%k)
          end if
          if(allocated(detector%update_vector)) then
             deallocate(detector%update_vector)
          end if
          allocate(detector%k(n_stages,dim))
          detector%k = 0.
          allocate(detector%update_vector(dim))
          detector%update_vector=0.
       end if
       detector => detector%next
    end do

  end subroutine allocate_rk_guided_search

  subroutine deallocate_rk_guided_search(detector_list)
    ! Deallocate the RK stages and update vector
    type(detector_linked_list), intent(inout) :: detector_list
      
    type(detector_type), pointer :: detector
      
    detector => detector_list%first
    do while (associated(detector))
       if(detector%type==LAGRANGIAN_DETECTOR) then
          if(allocated(detector%k)) then
             deallocate(detector%k)
          end if
          if(allocated(detector%update_vector)) then
             deallocate(detector%update_vector)
          end if
       end if
       detector => detector%next
    end do
  end subroutine deallocate_rk_guided_search

end module detector_move_lagrangian
