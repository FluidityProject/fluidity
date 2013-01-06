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
  use iso_c_binding
  use transform_elements
  use Profiler
  use linked_lists
  use integer_hash_table_module
  use pickers
  use parallel_fields, only: element_owner
  use ieee_arithmetic, only: ieee_is_nan
  use lebiology_python
  use path_element_module
  use element_path_list, only: elepath_list_insert => list_insert_head, &
                               elepath_list_create => list_create, &
                               elepath_list_destroy => list_destroy

  implicit none
  
  private

  public :: move_lagrangian_detectors, read_random_walk_options, initialise_rw_subcycling

contains

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
    type(scalar_field), pointer :: diffusivity_field, subcycle_field
    type(vector_field), pointer :: diffusivity_grad

    ! Periodic tracking
    type(mesh_type), pointer :: periodic_mesh
    integer, dimension(:), allocatable :: boundary_ids
    integer, dimension(2) :: shape_option
    integer :: bc, bid, n_periodic_bcs, boundary_id, mapping_count

    if (.not. check_any_lagrangian(detector_list)) return

    call profiler_tic(trim(detector_list%name)//"::movement")

    ewrite(1,*) "In move_lagrangian_detectors for detectors list: ", detector_list%name
    ewrite(2,*) "Detector list", detector_list%id, "has", detector_list%length, &
         "local and", detector_list%total_num_det, "global detectors"

    ! Pull some information from state
    xfield=>extract_vector_field(state(1), "Coordinate")
    vfield=>extract_vector_field(state(1), "Velocity")
    allocate(rw_displacement(xfield%dim))

    ! Prepare for tracking in periodic domains by storing 
    ! boundary_ids and corresponding Python mapping functions (only do this once...)
    ! and extracting the periodic mesh
    if (detector_list%track_periodic) then
       periodic_mesh => extract_mesh(state, trim(detector_list%periodic_tracking_mesh))

       ! Note: The following should be done during initialisation really...
       if (.not. allocated(detector_list%boundary_mappings)) then
          mapping_count = 1
          n_periodic_bcs = option_count(trim(periodic_mesh%option_path)//"/from_mesh/periodic_boundary_conditions")
          allocate( detector_list%boundary_mappings(n_periodic_bcs*2) )
          do bc=0, n_periodic_bcs-1
             ! Aliased BID uses coordinate_map
             shape_option = option_shape(trim(periodic_mesh%option_path)//"/from_mesh/periodic_boundary_conditions["//int2str(bc)//"]/aliased_boundary_ids")
             allocate( boundary_ids(shape_option(1)) )
             call get_option(trim(periodic_mesh%option_path)//"/from_mesh/periodic_boundary_conditions["//int2str(bc)//"]/aliased_boundary_ids",boundary_ids)
             do bid=1, size(boundary_ids)
                call insert(detector_list%bid_to_boundary_mapping, boundary_ids(bid), mapping_count)
             end do
             deallocate(boundary_ids)
             call get_option(trim(periodic_mesh%option_path)//"/from_mesh/periodic_boundary_conditions["//int2str(bc)//"]/coordinate_map", detector_list%boundary_mappings(mapping_count))
             mapping_count = mapping_count + 1

             ! Physical BIDs uses inverse_coordinate_map
             shape_option = option_shape(trim(periodic_mesh%option_path)//"/from_mesh/periodic_boundary_conditions["//int2str(bc)//"]/physical_boundary_ids")
             allocate( boundary_ids(shape_option(1)) )
             call get_option(trim(periodic_mesh%option_path)//"/from_mesh/periodic_boundary_conditions["//int2str(bc)//"]/physical_boundary_ids",boundary_ids)
             do bid=1, size(boundary_ids)
                call insert(detector_list%bid_to_boundary_mapping, boundary_ids(bid), mapping_count)
             end do
             deallocate(boundary_ids)
             call get_option(trim(periodic_mesh%option_path)//"/from_mesh/periodic_boundary_conditions["//int2str(bc)//"]/inverse_coordinate_map", detector_list%boundary_mappings(mapping_count))
             mapping_count = mapping_count + 1
          end do

       end if
    else
       periodic_mesh => null()
    end if

    ! Setup for Random Walk schemes
    if (allocated(detector_list%random_walks)) then
       do rw=1, size(detector_list%random_walks)

          ! For Python Random Walk run the user code that pulls fields from state
          if (detector_list%random_walks(rw)%python_random_walk) then
             call lebiology_prepare_pyfunc(detector_list%fgroup, trim(detector_list%random_walks(rw)%name), &
                       trim(detector_list%random_walks(rw)%python_code) )
          end if

          ! For hardcoded Random Walks pull the relevant fields from state
          if (detector_list%random_walks(rw)%naive_random_walk) then
             diffusivity_field=>extract_scalar_field(state(1), trim(detector_list%random_walks(rw)%diffusivity_field))
          end if

          if (detector_list%random_walks(rw)%diffusive_random_walk) then
             diffusivity_field=>extract_scalar_field(state(1), trim(detector_list%random_walks(rw)%diffusivity_field))
             diffusivity_grad=>extract_vector_field(state(1), trim(detector_list%random_walks(rw)%diffusivity_grad))
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

          if (detector_list%tracking_method == GEOMETRIC_TRACKING .or. &
             detector_list%tracking_method == PURE_GS) then
             if (allocated(detector%ray_o)) deallocate(detector%ray_o)
             allocate(detector%ray_o(xfield%dim))
             detector%ray_o = detector%position

             if (associated(detector%path_elements)) then
                call elepath_list_destroy(detector%path_elements)
             end if
          end if

          if (detector_list%tracking_method == GEOMETRIC_TRACKING) then
             if (allocated(detector%ray_d)) deallocate(detector%ray_d)
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
          RKstages_loop: do stage = 1, detector_list%n_stages

             ! Compute the update vector for the current stage
             call set_stage(detector_list,vfield,xfield,sub_dt,stage)

             ! Move update parametric detector coordinates according to update_vector
             ! If this takes a detector across parallel domain boundaries
             ! the routine will also send the detector
             call track_detectors(detector_list,xfield,periodic_mesh)

          end do RKstages_loop
       end if

       if (allocated(detector_list%random_walks)) then
          do rw=1, size(detector_list%random_walks)

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
                         call lebiology_move_agent(detector_list%fgroup, trim(detector_list%random_walks(rw)%name), &
                                   detector, sub_dt, rw_displacement)
                      end if

                      detector%update_vector=detector%update_vector + rw_displacement
                   end if
                   detector => detector%next
                end do

                call track_detectors(detector_list,xfield,periodic_mesh)

             ! Internal Diffusive Random Walk with automated sub-cycling
             elseif (detector_list%random_walks(rw)%diffusive_random_walk &
                         .and. detector_list%random_walks(rw)%auto_subcycle) then

                subcycle_field => extract_scalar_field(state(1), "RandomWalkSubcycling")
                call auto_subcycle_random_walk(detector_list,dt,&
                          xfield,diffusivity_field,diffusivity_grad,subcycle_field,periodic_mesh)

             else
                FLExit("Auto-subcycling can only be used with internal diffusive Random Walk.")
             end if

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

  subroutine track_detectors(detector_list,xfield,periodic_mesh)
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
    type(mesh_type), pointer, intent(in) :: periodic_mesh

    type(detector_type), pointer :: detector, send_detector
    type(detector_linked_list), dimension(:), allocatable :: send_list_array
    real :: search_tolerance
    integer :: k, nprocs, new_owner, all_send_lists_empty
    logical :: outside_domain, any_lagrangian

    ! We allocate a sendlist for every processor
    nprocs=getnprocs()
    allocate(send_list_array(nprocs))

    search_tolerance=detector_list%search_tolerance

    if (detector_list%tracking_method == GEOMETRIC_TRACKING) then
       detector => detector_list%first
       do while (associated(detector))

           !Only initialise Lagrangian detectors
           if (.not. detector%type==LAGRANGIAN_DETECTOR) then
              detector => detector%next
              cycle
           end if

          ! Calcualte and normalise ray direction
          call initialise_ray(detector, detector%update_vector)

          detector => detector%next
       end do
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
                if (associated(periodic_mesh)) then
                   call geometric_ray_tracing(xfield,detector,detector_list,new_owner,search_tolerance,periodic_mesh=periodic_mesh)
                else
                   call geometric_ray_tracing(xfield,detector,detector_list,new_owner,search_tolerance)
                end if
             elseif (detector_list%tracking_method == PURE_GS) then
                call parametric_guided_search(xfield,detector,new_owner,search_tolerance)
             else
                FLExit('No lagrangian particle tracking method specified')
             end if

             if (new_owner==-1) then
                if (detector_list%reflect_on_boundary) then
                   ! We reflect the detector path at the face we just went through
                   if (detector_list%tracking_method == GEOMETRIC_TRACKING) then
                      call reflect_on_boundary(xfield,detector%update_vector,detector%element, face=detector%current_face)

                      ! After reflection our geometry changes:
                      ! We move our origin to the last known point on our ray (ele_t), 
                      ! and calculate a new direction
                      detector%ray_o = detector%ray_o + (detector%current_t * detector%ray_d)
                      call initialise_ray(detector, detector%update_vector)
                   elseif (detector_list%tracking_method == PURE_GS) then
                      call reflect_on_boundary(xfield,detector%update_vector,detector%element, face=detector%current_face)
                   else
                      call reflect_on_boundary(xfield,detector%update_vector,detector%element)
                   end if
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
          call exchange_detectors(detector_list,xfield,send_list_array,tracking=.true.)
       else
          ! If we run out of lagrangian detectors for some reason, exit the loop
          exit
       end if
    end do detector_timestepping_loop

    deallocate(send_list_array)

    detector => detector_list%first
    do while (associated(detector))

       !Only update Lagrangian detectors
       if (.not. detector%type==LAGRANGIAN_DETECTOR) then
          detector => detector%next
          cycle
       end if

       if (detector_list%tracking_method == GEOMETRIC_TRACKING) then
          detector%ray_o = detector%update_vector

          detector%local_coords = local_coords(xfield,detector%element,detector%update_vector)
       end if

       if (minval(detector%local_coords) < -detector_list%search_tolerance) then
          ewrite(-1,*) "Detector", detector%id_number, ", in element", detector%element, " has local coordinates: ", detector%local_coords
          FLAbort("Negative local coordinate for lagrangian detector after tracking!")
       end if

       detector => detector%next
    end do

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

  subroutine initialise_ray(detector, target_coord, origin)
    type(detector_type), pointer, intent(inout) :: detector
    real, dimension(:), intent(in) :: target_coord
    real, dimension(:), intent(in), optional :: origin

    if (present(origin)) detector%ray_o = origin

    ! Calculate and normalise ray direction
    detector%ray_d = target_coord - detector%ray_o
    detector%target_distance = sum(detector%ray_d**2)
    if (detector%target_distance > 0.0) then
       detector%target_distance = sqrt(detector%target_distance)
       detector%ray_d = detector%ray_d / detector%target_distance
    else
       detector%target_distance = 0.0
       detector%ray_d = 0.0
    end if
    detector%current_t = 0.0

  end subroutine initialise_ray

  subroutine geometric_ray_tracing(xfield,detector,detector_list,new_owner,search_tolerance,periodic_mesh)
    ! This tracking method is based on a standard Ray-tracing algorithm using planes and half-spaces.
    ! Reference: 
    type(vector_field), pointer, intent(in) :: xfield
    type(detector_type), pointer, intent(inout) :: detector
    type(detector_linked_list), intent(inout) :: detector_list
    integer, intent(out) :: new_owner
    real, intent(in) :: search_tolerance
    type(mesh_type), pointer, intent(in), optional :: periodic_mesh

    real :: face_t, ele_t
    integer :: i, neigh_face, next_face, boundary_id, py_map_index
    integer, dimension(:), pointer :: face_list, element_nodes
    real, dimension(xfield%dim,1):: old_position, new_position
    real, dimension(xfield%dim) :: face_normal, face_node_val
    type(path_element) :: path_ele

    ! Exit if r_d is zero
    if (detector%target_distance < search_tolerance) then
       new_owner=getprocno()
       return
    end if

    search_loop: do 

       ! Go through all faces of the element and look for the smallest t in the ray's direction
       next_face = -1
       ele_t = huge(1.0)
       face_list => ele_faces(xfield,detector%element)
       element_nodes => ele_nodes(xfield,detector%element)
       do i=1, size(face_list)
          call ray_intersetion_distance(face_list(i), element_nodes, face_normal, face_node_val, face_t)

          ! The second condition is to ensure that we do not trace backwards
          if (face_t < ele_t .and. detector%current_t < face_t .and. search_tolerance < face_t) then
             ele_t = face_t
             next_face = face_list(i)
          end if
       end do

       ! This is an extreme corner case, where we ended the last timestep 
       ! exactly on a domain boundary, which is not periodic.
       ! In this case we go through the faces again, to find the one we're on, 
       ! so we can set ele_t and next_face to get the update_vector reflected
       if (ele_t == huge(1.0)) then

          ele_t = detector%current_t
          do i=1, size(face_list)
            call ray_intersetion_distance(face_list(i), element_nodes, face_normal, face_node_val, face_t)
            if (abs(face_t) < search_tolerance) then
               next_face = face_list(i)
            end if
         end do
       end if

       if (ele_t < detector%target_distance + search_tolerance) then
          neigh_face = face_neigh(xfield, next_face)

          ! Store current face, in case we need it for reflection
          detector%current_face = next_face

          ! Record the elements along the path travelled
          ! and the distance travelled within them
          path_ele%ele = detector%element
          path_ele%dist = ele_t - detector%current_t
          if (associated(detector%path_elements)) then
             call elepath_list_insert( detector%path_elements, path_ele )
          else
             call elepath_list_create( detector%path_elements, path_ele )
          end if
          detector%current_t = ele_t

          if (neigh_face /= next_face) then
             ! Recurse on the next element
             detector%element = face_ele(xfield, neigh_face)

          else
             ! If we're tracking on a periodic mesh...
             !  * we find out which boundary we have hit
             !  * translate our udpate_vector via the python function
             !  * and keep tracking on the periodic (connected mesh).
             !    Since we don't use the update vector once the ray parameters are established,
             !    this is legitimate and should give the right result
             !    Note: This assumes that the Python mapping function works universally!
             if (present(periodic_mesh)) then
                boundary_id = surface_element_id(xfield, next_face)
                py_map_index = fetch(detector_list%bid_to_boundary_mapping, boundary_id)

                ! Apply the python map
                old_position(:,1) = detector%update_vector
                call set_from_python_function(new_position, trim(detector_list%boundary_mappings(py_map_index)), old_position, time=0.0)
                detector%update_vector = new_position(:,1)

                ! New update_vector is the target, and the origin has changed
                old_position(:,1) = detector%ray_o + detector%ray_d * detector%current_t
                call set_from_python_function(new_position, trim(detector_list%boundary_mappings(py_map_index)), old_position, time=0.0)
                detector%ray_o = new_position(:,1)
                call initialise_ray(detector,detector%update_vector)

                neigh_face = face_neigh(periodic_mesh, next_face)
                if (neigh_face /= next_face) then
                   ! Recurse on the next element
                   detector%element = face_ele(xfield, neigh_face)
                   cycle search_loop
                end if
             end if

             if (element_owned(xfield,detector%element)) then
                ! Detector is going outside domain
                new_owner=-1
                exit search_loop
             else
                ! The current element is on a Halo, we need to send it to the owner.
                new_owner=element_owner(xfield%mesh,detector%element)
                exit search_loop
             end if
          end if
       else
          ! The arrival point is in this element, we're done
          path_ele%ele = detector%element
          path_ele%dist = detector%target_distance - detector%current_t
          if (associated(detector%path_elements)) then
             call elepath_list_insert( detector%path_elements, path_ele )
          else
             call elepath_list_create( detector%path_elements, path_ele )
          end if

          new_owner=getprocno()
          exit search_loop
       end if
       
    end do search_loop
    
  contains

    subroutine ray_intersetion_distance(face, ele_nodes, face_normal, face_node_val, t)
      ! Calculate the distance t of the intersection 
      ! between the face and our ray (r_o, r_d)
      integer, intent(in) :: face
      integer, dimension(:), pointer :: ele_nodes 
      ! Alocate these outside to increase performance
      real, dimension(:), intent(inout) :: face_normal, face_node_val
      real, intent(out) :: t

      real :: d, v_n, v_d, detJ
      integer, dimension(face_loc(xfield, face)) :: face_nodes
      logical :: face_cache_valid

      ! Get face normal
      face_cache_valid = .false.
      face_cache_valid = retrieve_cached_face_transform(xfield, face, face_normal, detJ)
      assert(face_cache_valid)

      ! Establish d using a point on the plane
      face_nodes = face_local_nodes(xfield, face)
      face_node_val = node_val(xfield,ele_nodes(face_nodes(1)))
      d = - sum(face_normal * face_node_val)

      ! Calculate t, the distance along the ray to the face intersection
      v_n = dot_product(face_normal, detector%ray_o) + d
      v_d = dot_product(face_normal, detector%ray_d)
      t = - v_n / v_d

    end subroutine ray_intersetion_distance

  end subroutine geometric_ray_tracing

  subroutine parametric_guided_search(xfield, detector, new_owner, search_tolerance)
    ! Implementation of Guided Search tracking algorithm purely in parametric space
    ! This works generaically for any detector displacement, 
    ! by transforming the displcement vector to parametric space via the inverse Jacobian,
    ! computing the intersection with the element boundary, 
    ! and updating the local coordinates accordingly
    type(vector_field), pointer, intent(in) :: xfield
    type(detector_type), pointer, intent(inout) :: detector
    integer, intent(out) :: new_owner
    real, intent(in) :: search_tolerance

    ! Inverse of the local coordinate change matrix.
    real, dimension(mesh_dim(xfield), mesh_dim(xfield), ele_ngi(xfield, 1)) :: invJ
    ! Displacement vector in physical space
    real, dimension(xfield%dim) :: vdisp, last_position, displacement
    real, dimension(xfield%dim+1) :: lc_update, lc_vdisp, dist_intersect
    real, dimension(ele_loc(xfield, detector%element), ele_ngi(xfield, detector%element), xfield%dim) :: dshape
    real, dimension(ele_ngi(xfield, detector%element)) :: old_detJ, new_detJ
    integer :: i, dim, neigh, local_face
    integer, dimension(:), pointer :: neigh_list, face_list
    type(path_element) :: path_ele
    
    dim = xfield%dim   
    search_loop: do 
       ! Displacement vector, where ray_o is the origin of the tracking step
       vdisp = detector%update_vector - detector%ray_o

       ! Get the inverse Jacobian and detJ for the old element
       call transform_to_physical(xfield, detector%element, xfield%mesh%shape, dshape, invJ=invJ, detJ=old_detJ)

       ! Compute equivalent of detector%update_vector in parametric space
       lc_update(1:dim) = detector%local_coords(1:dim) + matmul(vdisp, invJ(:,:,1))
       lc_update(dim+1) = 1. - sum(lc_update(1:dim))

       if (minval(lc_update)>-search_tolerance) then
          !The arrival point is in this element, we're done
          detector%local_coords = lc_update

          ! Record the elements along the path travelled
          ! and the distance travelled within them
          path_ele%ele = detector%element
          displacement = detector%update_vector - detector%ray_o
          path_ele%dist = sqrt( sum(displacement**2) )
          if (associated(detector%path_elements)) then
             call elepath_list_insert( detector%path_elements, path_ele )
          else
             call elepath_list_create( detector%path_elements, path_ele )
          end if

          detector%ray_o = detector%update_vector
          new_owner=getprocno()

          exit search_loop
       end if

       ! Find closest intersection with an element boundary in parametric space, 
       ! and derive the according element neighbour
       lc_vdisp = lc_update - detector%local_coords
       dist_intersect = - detector%local_coords / lc_vdisp
       neigh = minval(minloc(dist_intersect, dist_intersect > search_tolerance))
       detector%local_coords = detector%local_coords + lc_vdisp * dist_intersect(neigh)
       last_position = detector%ray_o
       detector%ray_o = detector%ray_o + vdisp * dist_intersect(neigh)
       
       ! Record the elements along the path travelled
       ! and the distance travelled within them
       path_ele%ele = detector%element
       displacement = detector%ray_o - last_position
       path_ele%dist = sqrt( sum(displacement**2) )
       if (associated(detector%path_elements)) then
          call elepath_list_insert( detector%path_elements, path_ele )
       else
          call elepath_list_create( detector%path_elements, path_ele )
       end if

       ! Record the face we went through
       neigh_list=>ele_neigh(xfield,detector%element)
       face_list=>ele_faces(xfield,detector%element)
       detector%current_face = face_list(neigh)
       
       if (neigh > 0 .and. neigh_list(neigh)>0) then
          ! The neighbouring element is also on this domain:
          ! Now construct the local_coords of the intersection point 
          ! on the element boundary wrt. to the new element
          detector%element = neigh_list(neigh) 

          ! The detJ seems to indicate whether we need to 
          ! inverse the internal element numbering or not ...
          ! Note, however, that the transform is virtually for free thanks to caching
          call transform_to_physical(xfield, detector%element, xfield%mesh%shape, dshape, invJ=invJ, detJ=new_detJ)
          local_face = local_face_number(xfield, face_neigh(xfield, detector%current_face))

          ! For some reason, if the detJ is negative the element connectivity is reversed?
          if (old_detJ(1) * new_detJ(1) > search_tolerance) then 
             detector%local_coords = detector%local_coords(ubound(detector%local_coords,1):lbound(detector%local_coords,1):-1)
             neigh = dim+2 - neigh
          end if

          ! Now we can shift the local coords according to the difference in local number
          detector%local_coords = cshift(detector%local_coords, shift=neigh-local_face)
       else
          ! Update vector is not on this domain, 
          ! so check if the last element is owned
          if (element_owned(xfield,detector%element)) then
             ! Detector leaving local domain
             new_owner=-1
          else
             ! The current element is on a Halo, we need to send it to the owner.
             new_owner=element_owner(xfield%mesh,detector%element)
          end if
          exit search_loop
       end if

    end do search_loop
  end subroutine parametric_guided_search

  subroutine reflect_on_boundary(xfield, coordinate, element, face)
    ! Reflect the coordinate in the according boundary face of the element.
    ! This uses simple reflection equations and assumes the face to be flat.
    type(vector_field), pointer, intent(in) :: xfield
    real, dimension(mesh_dim(xfield)), intent(inout) :: coordinate
    integer, intent(in) :: element
    integer, intent(in), optional :: face

    real, dimension(mesh_dim(xfield)+1) :: arrival_local_coords
    integer :: i, neigh, neigh_face
    real, dimension(xfield%dim,xfield%mesh%faces%shape%ngi) :: facet_normals
    integer, dimension(:), pointer :: neigh_list
    integer, dimension(:), allocatable :: face_nodes 
    real, dimension(xfield%dim) :: face_node_val, face_normal
    real :: offset, p, D, arrival_lcoord

    ! First, if its not given, get the face we went through
    if (present(face)) then
       neigh_face = face
    else
       neigh_list=>ele_neigh(xfield,element)
       arrival_local_coords=local_coords(xfield,element,coordinate)
       neigh = minval(minloc(arrival_local_coords))
       neigh_face = ele_face(xfield, element, neigh_list(neigh))
    end if

    ! Now get the coordinates of one node on the face we went through
    allocate(face_nodes(face_loc(xfield, neigh_face)))
    face_nodes = face_global_nodes(xfield, neigh_face)    
    face_node_val=node_val(xfield,face_nodes(1))

    ! Now we get the face normal from the transform (detwei_f is a dummy)
    call transform_facet_to_physical(xfield, neigh_face, normal=facet_normals)
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

  subroutine auto_subcycle_random_walk(detector_list,dt,xfield,diff_field,grad_field,subcycle_field,periodic_mesh)
    type(detector_linked_list), intent(inout) :: detector_list
    real, intent(in) :: dt
    type(scalar_field), pointer, intent(in) :: diff_field, subcycle_field
    type(vector_field), pointer, intent(in) :: xfield, grad_field
    type(mesh_type), pointer, intent(in) :: periodic_mesh

    type(detector_linked_list) :: subcycle_detector_list
    type(detector_type), pointer :: detector, move_detector
    real, dimension(xfield%dim) :: rw_displacement
    real :: total_subsubcycle
    integer :: ele, subsubcycle, max_subsubcycle , remaining_detectors
    integer, dimension(1) :: subc

    max_subsubcycle = 1
    total_subsubcycle = 0

    ! Setup our detector work list
    subcycle_detector_list%search_tolerance = detector_list%search_tolerance
    subcycle_detector_list%reflect_on_boundary = detector_list%reflect_on_boundary
    subcycle_detector_list%tracking_method = detector_list%tracking_method
    subcycle_detector_list%name = detector_list%name

    ! First pass establishes how many sub-sub-cycles are required
    detector => detector_list%first
    do while (associated(detector))
       if (detector%type==LAGRANGIAN_DETECTOR) then

          ! Establish how many subcycles are needed
          subc = ele_val(subcycle_field, detector%element)
          detector%rw_subsubcycles = subc(1)
          total_subsubcycle = total_subsubcycle + detector%rw_subsubcycles
          if (detector%rw_subsubcycles > max_subsubcycle) then
             max_subsubcycle = detector%rw_subsubcycles
          end if

          call diffusive_random_walk(detector, dt/detector%rw_subsubcycles, &
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

    ! Update elements in the subcycle list
    if (subcycle_detector_list%length > 0) then
       call track_detectors(subcycle_detector_list,xfield,periodic_mesh)
    end if


    ! Remaining subcycle passes until subcycle list is empty
    subsubcycle = 2
    remaining_detectors = subcycle_detector_list%length
    call allmax(remaining_detectors)
    do while (remaining_detectors > 0)
       detector => subcycle_detector_list%first
       do while (associated(detector))
          call diffusive_random_walk(detector, dt/detector%rw_subsubcycles, &
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
       call track_detectors(subcycle_detector_list,xfield,periodic_mesh)

       remaining_detectors = subcycle_detector_list%length
       call allmax(remaining_detectors)

       subsubcycle = subsubcycle + 1
    end do

    ! Update detectors after the final move
    call track_detectors(detector_list,xfield,periodic_mesh)

    ewrite(2,*) "Auto-subcycling: max. subcycles:", max_subsubcycle
    ewrite(2,*) "Auto-subcycling: avg. subcycles:", total_subsubcycle / detector_list%length

  end subroutine auto_subcycle_random_walk

  subroutine initialise_rw_subcycling(state, xfield, option_path, dt)
    type(state_type), intent(inout) :: state
    type(vector_field), pointer, intent(inout) :: xfield
    character(len=OPTION_PATH_LEN), intent(in) :: option_path
    real, intent(in) :: dt

    type(scalar_field), pointer :: diffusivity, subcycle_field
    type(tensor_field), pointer :: diffusivity_hessian
    integer :: n, ele
    integer, dimension(:), pointer :: nodes
    real :: scale_factor, subcycles

    ewrite(2,*) "In initialise_rw_subcycling"

    call get_option(trim(option_path)//"/scalar_field::RandomWalkSubcycling/scale_factor", scale_factor)
    subcycle_field => extract_scalar_field(state, "RandomWalkSubcycling")
    call zero(subcycle_field)

    diffusivity => extract_scalar_field(state, "ScalarDiffusivity")
    diffusivity_hessian => extract_tensor_field(state, "DiffusivityHessian")

    do ele=1, element_count(subcycle_field)
       call element_rw_subcycling(ele, dt, scale_factor, subcycles)
       nodes => ele_nodes(subcycle_field, ele)
       call addto(subcycle_field, ele, subcycles)
    end do

  contains

    subroutine element_rw_subcycling(element, dt, scale_factor, subcycles)
      integer, intent(in) :: element
      real, intent(in) :: dt, scale_factor
      real, intent(out) :: subcycles

      integer :: i, subc, local_subcycling
      real :: k0, d_z, min_z, max_z, min_dt
      real, dimension(xfield%mesh%shape%loc) :: coords
      real, dimension(diffusivity_hessian%mesh%shape%loc) :: k_grad2

      real, dimension(xfield%dim) :: low_coord, high_coord
      integer, dimension(:), pointer :: picked_elements
      integer :: dim

      ! Calculate the vertical range for the criteria
      k0 = maxval(abs(ele_val(diffusivity, element)))
      d_z = sqrt(2 * k0 * dt)

      coords = ele_val(xfield, xfield%dim, element)
      min_z = minval(coords) - d_z
      max_z = maxval(coords) + d_z

      dim = xfield%dim
      do i=1, dim-1
         low_coord(i) = minval(ele_val(xfield, i, element))
         high_coord(i) = maxval(ele_val(xfield, i, element))
      end do
      low_coord(dim) = min_z
      high_coord(dim) = max_z

      ! Query the picker for all elements in this region
      call picker_inquire_region(xfield, low_coord, high_coord, picked_elements)

      ! Now we have all elements in our range and we can find the maximum sub-cycle number
      subc = 1
      do i=1, size(picked_elements)
         k_grad2 = ele_val(diffusivity_hessian, dim, dim, picked_elements(i))
         k_grad2 = 1.0 / abs(k_grad2)
         min_dt = minval(k_grad2) / scale_factor
         local_subcycling = ceiling(dt/min_dt)
         if (local_subcycling > subc) then
            subc = local_subcycling
         end if
      end do

      deallocate(picked_elements)
      subcycles = real(subc)

    end subroutine element_rw_subcycling

  end subroutine initialise_rw_subcycling

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

end module detector_move_lagrangian
