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

module detector_move_bisection
  use state_module
  use fields
  use integer_hash_table_module
  use detector_data_types
  use detector_tools
  use halo_data_types
  use halos_numbering

  implicit none
  
  private

  public :: move_detectors_bisection_method

contains

  subroutine move_detectors_bisection_method(state, detector_list, dt, ihash, send_list_array)
    !!< Move Lagrangian detectors using a bisecting method, i.e., dividing the dt in smaller values that add to dt. 
    !Each smaller value is such that the detector is moved from its previous position to the boundary face between    
    !the element where it belonged to and the neighbouring element in the direction of the flow, through the minimum local 
    !coordinate. Sometimes going through the minimum coordinate is not right (does not follow the flow direction) and in 
    !that case the detector has been placed back into the previous position and previous element and care has been taken for the 
    !detector to follow the flow until the next boundary face
    type(state_type), dimension(:), intent(in) :: state
    type(detector_linked_list), intent(inout) :: detector_list
    real, intent(in) :: dt
    type(integer_hash_table), intent(in) :: ihash
    type(detector_linked_list), dimension(:), allocatable, intent(inout) :: send_list_array

    integer :: j, number_neigh_processors
    real, dimension(:), allocatable :: vel, old_vel, old_pos
    type(vector_field), pointer :: vfield, xfield, old_vfield
    type(detector_type), pointer :: this_det
    integer, dimension(:), pointer :: ele_number_ptr
    real :: dt_temp
    integer :: index_next_face, current_element, previous_element, & 
               & cont_repeated_situation, processor_number

    ewrite(1,*) "Inside move_detectors_bisection_method subroutine"

    vfield => extract_vector_field(state(1),"Velocity")
    old_vfield => extract_vector_field(state(1),"OldVelocity")
    xfield => extract_vector_field(state(1),"Coordinate")

    allocate(vel(vfield%dim))
    allocate(old_vel(old_vfield%dim))
    allocate(old_pos(old_vfield%dim))

    number_neigh_processors=key_count(ihash)

    this_det => detector_list%firstnode

    do j=1, detector_list%length
             
       if (this_det%type /= LAGRANGIAN_DETECTOR) then            
          this_det => this_det%next
          cycle
       end if

       dt_temp=0.0

       previous_element = this_det%element

       cont_repeated_situation=0.0

       !do_until_whole_dt_is_reached or detector leaves the domain towards another processor or through an outflow: 
       do  
         
          if (this_det%dt<1.0e-3*dt) exit
             
          vel =  detector_value(vfield, this_det)
          old_vel = detector_value(old_vfield, this_det)
          old_pos = this_det%position
          current_element = this_det%element

          call move_detectors_subtime_step(this_det, xfield, dt, dt_temp, old_pos, vel, old_vel, vfield, old_vfield,previous_element,index_next_face)

! From the previous subroutine, the detector is moved until the boundary between two elements is found or all smaller time steps have reached dt. 
! If the detector leaves the domain or we lose track of it, then it will be converted into a static one

          processor_number=getprocno()

          if (this_det%type==STATIC_DETECTOR) exit

          ele_number_ptr=>ele_neigh(xfield,this_det%element)
          previous_element=this_det%element
          this_det%element=ele_number_ptr(index_next_face) 

          processor_number=getprocno()

          call check_if_det_gone_through_domain_boundary(state, this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, vfield, old_vfield, index_next_face, current_element, &
                                         cont_repeated_situation,detector_list,send_list_array,ihash,processor_number)

          if (processor_number /= getprocno()) exit
     
          if (this_det%type==STATIC_DETECTOR) exit

          this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)
   
          dt_temp=dt_temp+this_det%dt
          this_det%dt=dt-dt_temp;

          if (this_det%dt<1.0e-3*dt) then
             this_det%element=previous_element
             this_det%local_coords=local_coords(xfield,this_det%element,this_det%position) 
          end if

       end do 
       !do_until_whole_dt_is_reached or detector leaves the domain towards another processor or through an outflow: 

       if (processor_number == getprocno()) then 
       !if proc_number /= getprocno() this_det has already been made to point to this_det%next inside 
       !in check_if_det_gone_through_domain_bound.

          this_det => this_det%next
       end if

    end do  !do_for_each_detector
   
    deallocate(vel)
    deallocate(old_vel)
    deallocate(old_pos)

  end subroutine move_detectors_bisection_method

  subroutine move_detectors_subtime_step(this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, vfield, old_vfield, previous_element,index_next_face)
  !!< Subroutine that makes sure the Lagrangian detector ends up on one of its boundary faces (within tolerance of
  !+/-10.0e-8) up to where it has reached using one of the smaller time steps (sutime step). 
  !Once the detector is on the boundary (within tolerance of +/-10.0e-8), in the subroutine that calls this one, the   
  !detector is asigned to the neighbouring element through that face and again calling this subroutine, the detector is 
  !moved using another smaller time step until the boundary of this new element, always following the direction of the flow. 
  !These steps are repeated until the sum of all the smaller time steps used to go from elemenet to element is equal to the 
  !total time step.
!    type(state_type), dimension(:), intent(in) :: state

    type(detector_type), pointer :: this_det
    type(vector_field), intent(inout) :: xfield
    real, intent(in) :: dt
    real, intent(inout) :: dt_temp
    real,  dimension(:), intent(inout) :: old_pos, vel, old_vel
    type(vector_field), pointer :: vfield, old_vfield
    integer, intent(in) :: previous_element
    integer, intent(inout) :: index_next_face
    real,  dimension(1:size(this_det%local_coords)) :: bbb

    real :: keep_value_this_det_dt, keep_value_this_det_dt_a
    integer :: cont, cont_a, cont_b, cont_c, cont_d, cont_p, cont_r, &
               cont_check, index_minloc_current, cont_check_a, &
               cont_static_det, bound_elem_iteration

    ewrite(1,*) "Inside move_detectors_subtime_step subroutine"
    ewrite(-1,*) 'WARNING, BISECTION METHOD NOT RECOMMENDED!'

    cont=0
    cont_a=0
    cont_b=0
    cont_c=0
    cont_d=0
    cont_p=0
    cont_r=0
    cont_check=0
    cont_check_a=0

    do
   
       this_det%dt=this_det%dt/2.0
       this_det%position=((vel+old_vel)/2.0)*this_det%dt+old_pos

       this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

       keep_value_this_det_dt=this_det%dt

       cont=0

       do 
   
          if (all(this_det%local_coords>-10.0e-8)) exit

          this_det%dt=this_det%dt/2.0
       
          call update_detector_position_bisect(this_det, xfield, old_pos, vel, old_vel)

          cont=cont+1

          bbb=this_det%local_coords

          !If no matter how much this_det%dt is reduced, still the detector is not inside the new element when start moving it 
          !in the direction of the flow, it is most likely because the detector left the previous element through the wrong 
          !face. The condition of leaving the element is normally through the face with respect to which the local coordinate 
          !has the minimum value. This can lead sometimes to the detector leaving the element through a face that does not intersect 
          !the direction of the flow. In the next lines, the detector is returned to the previous element and previous position and 
          !it is checked that moving it from there in the direction of the flow, the detector stays in the element at least for a specific 
          !value of this_det%dt (smaller time step).
      
          if ((cont>1.0e6).or.(this_det%dt<1.0e-5*keep_value_this_det_dt)) then

             this_det%element=previous_element
             this_det%position=old_pos
             this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

             vel =  detector_value(vfield, this_det)
             old_vel = detector_value(old_vfield, this_det)

             index_minloc_current=minloc(this_det%local_coords, dim=1)
             cont_check=cont_check+1                 
             this_det%dt=keep_value_this_det_dt

             this_det%position=(vel+old_vel)/2*this_det%dt+old_pos
             this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

             cont_r=0
             cont_check_a=0

             do 

                if (all(this_det%local_coords>-10.0e-8)) exit

                this_det%dt=this_det%dt/2
         
                call update_detector_position_bisect(this_det, xfield, old_pos, vel, old_vel)

                cont_r=cont_r+1

                if ((cont_r>=1.0e6).or.(this_det%dt<1.0e-5*keep_value_this_det_dt)) then
                     cont_check_a=cont_check_a+1
                end if

                if (cont_check_a/=0) exit

             end do
          end if

          if  (cont_check_a/=0) exit

       end do

       vel = detector_value(vfield, this_det)
       old_vel = detector_value(old_vfield, this_det)
       this_det%dt=this_det%dt*2.0
       this_det%position=((vel+old_vel)/2.0)*this_det%dt+old_pos

       keep_value_this_det_dt_a=this_det%dt
       this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

       if (all(this_det%local_coords>-10.0e-8)) exit

       cont_d=cont_d+1

       this_det%dt=this_det%dt/2.0
       this_det%position=old_pos
       this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)
       vel = detector_value(vfield, this_det)
       old_vel = detector_value(old_vfield, this_det)

       if  (cont_check_a/=0) exit

    end do

    !dt_var=this_det%dt

    !In the next if, since the detector seems to have gone through the wrong face and still when returning it to its   
    !previous element and moving it with the flow, the detector does not stay in the element for a specific value of dt_var 
    !(unless it is for dt_var practically zero, so the detector is in its initial position in the element), it could be that 
    !it is going through a vortex. In this case, in order to find the next element to which the detector belongs when moving 
    !in the direction of the flow, all the elements that share a node are checked. The node chosen is the one opposite to the 
    !face with respect to which the local coordinate is the maximum.

    call scenario_det_gone_through_element_vortex(this_det, xfield, dt, dt_temp, keep_value_this_det_dt, old_pos, &
                                         vel, old_vel, vfield, old_vfield, previous_element,index_next_face,cont_check,cont_check_a,bound_elem_iteration)

    !In the next if, since the detector seems to have gone through the wrong face and when returning it to its   
    !previous element and moving it with the flow, the detector stays in the element for a specific value of dt_var,
    !now the detector is moved towards the appropriate face of the element.

    call scenario_det_gone_through_wrong_face(this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, index_next_face,cont_check,cont_check_a,bound_elem_iteration,index_minloc_current)
 
    !In the next if, the most straight forward scenario is dealed with, i.e., the detector previously on the boundary of an element, made to move 
    !in the direction of the flow towards the next element across that boundary, belongs to the next element for a particular value of dt_var
    !(next smaller time step).

    call scenario_det_gone_through_right_face(this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, index_next_face,cont_check,bound_elem_iteration)
 
    if (bound_elem_iteration>=33554432) then 
       ewrite(-1,*) 'Warning: converting detector to static'
       this_det%type=STATIC_DETECTOR
              
    end if     

    cont_static_det=0
    
    if (this_det%type==STATIC_DETECTOR) then

       cont_static_det=cont_static_det+1
          
       this_det%position=old_pos
       this_det%element=previous_element
       this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

    end if

    if (this_det%type==STATIC_DETECTOR) return

  end subroutine move_detectors_subtime_step

  subroutine check_if_det_gone_through_domain_boundary(state, this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, vfield, old_vfield, index_next_face, current_element, &
                                         cont_repeated_situation,detector_list,send_list_array,ihash,processor_number)

    type(state_type), dimension(:), intent(in) :: state

    type(detector_type), pointer :: this_det, node_to_send
    type(vector_field), intent(inout) :: xfield
    real, intent(in) :: dt
    real, intent(inout) :: dt_temp
    real,  dimension(:), intent(inout) :: old_pos, vel, old_vel
    type(vector_field), pointer :: vfield, old_vfield
    integer, intent(inout) :: index_next_face, current_element, cont_repeated_situation
    type(detector_linked_list), intent(inout) :: detector_list
    type(detector_linked_list), dimension(:), intent(inout) :: send_list_array
    type(integer_hash_table), intent(in) :: ihash
    integer, intent(inout) :: processor_number

    integer :: k,h, univ_ele, halo_level, ele
    real :: dt_var_temp, dt_temp_temp, dt_var_value
    integer :: old_index_minloc, cont_aaa, cont_bbb, cont_ccc, cont_ddd, cont_eee, cont_fff, & 
               cont_ggg, cont_times_same_element, list_neigh_processor, processor_number_curr_ele
    integer, dimension(:), allocatable :: element_next_to_boundary

    type(halo_type), pointer :: ele_halo 

    cont_aaa=0.0
    cont_bbb=0.0
    cont_ccc=0.0
    cont_ddd=0.0
    cont_eee=0.0
    cont_fff=0.0
    cont_ggg=0.0

    vfield => extract_vector_field(state(1),"Velocity")
    halo_level = element_halo_count(vfield%mesh)

    if (halo_level /= 0) then

    ele_halo => vfield%mesh%element_halos(halo_level)
  
    univ_ele = halo_universal_number(ele_halo, current_element)

    end if
    
    processor_number_curr_ele=element_owner(vfield,current_element)

    cont_times_same_element=0.0

    ele = current_element

    if (this_det%element<0.0) then

       if (element_owned(vfield,current_element)) then

          write(1,*) "CURRENT PROC OWNS the previous element:", current_element 

          !If I own the previous element where the detector was before it left the domain, it means it was not an halo element
          !of another processor, and hence, it has left through a proper boundary (outflow) and it is converted into static.
          !Or it can also be that by some issue in the geometry of the element for example, the detector is not being moved in the right 
          !direction and is leaving through the wrong face.
          !If after some checks done below, the detector is still leaving the domain through that face, it is converted into static.

          !current_element is the previous element where the detector was before moving into a negative element (outside the domain)

          cont_aaa=cont_aaa+1

          cont_repeated_situation=cont_repeated_situation+1
          if (cont_repeated_situation>99999) then
             cont_repeated_situation=1
          end if
       
          allocate(element_next_to_boundary(cont_repeated_situation))

          element_next_to_boundary(cont_repeated_situation)=current_element 

          this_det%element=current_element
          dt_var_value=this_det%dt
          this_det%dt=this_det%dt*0.6

          call update_detector_position_bisect(this_det, xfield, old_pos, vel, old_vel)
        
          do
            
            if ((all(this_det%local_coords>=0.0)).and.(this_det%local_coords(index_next_face)>1.0e-5)) exit

            this_det%dt=this_det%dt/2
          
            call update_detector_position_bisect(this_det, xfield, old_pos, vel, old_vel)

            cont_bbb=cont_bbb+1

            !If a lagrangian detector is initially in a boundary it does not work because the velocity is zero 
            !and gets stuck in this loop
            !Better not to put any detector in a boundary but in case it happens a check has been included so that 
            !the code does not get stuck.
            !The Lagrangian detector would be converted into a static one.

            if ((cont_bbb>1.0e6).or.(this_det%dt<1.0e-5*dt_var_value))  exit

          end do        

          if ((cont_bbb>1.0e6).or.(this_det%dt<1.0e-5*dt_var_value)) then

            cont_fff=cont_fff+1

            this_det%position=old_pos
            this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)
            ewrite(-1,*) 'Warning: converting detector to static'
            this_det%type=STATIC_DETECTOR
        
          end if 

         !Below a check is included in case the detector keeps going back to the same element next to the boundary. 
         !This means the velocity in that area is not moving the detector inside the domain away from the boundary 
         !element, and the detector continuously ends up in the boundary. In this situation the code will keep 
         !entering here to place the detector inside the same element next to the boundary and the detector 
         !does not progress any more. 

         do h=1, size(element_next_to_boundary)

            if (element_next_to_boundary(h)==current_element) then

              cont_times_same_element=cont_times_same_element+1

            end if

         end do
     
         if ((cont_repeated_situation>50.0).and.(cont_times_same_element>50.0)) then

            cont_ggg=cont_ggg+1

            this_det%position=old_pos
            this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)
            ewrite(-1,*) 'Warning: converting detector to static'
            this_det%type=STATIC_DETECTOR

         end if

         if (this_det%type==STATIC_DETECTOR) return

         vel =  detector_value(vfield, this_det)
         old_vel = detector_value(old_vfield, this_det)
         old_pos = this_det%position

!Next we start moving a bit the detector after placing it back into the previous element and check that this time it does not go through the same boundary and hence, ends outside the domain. If that is the case, then we make the detector static already now.             

         dt_temp_temp=dt_temp+this_det%dt
         dt_var_temp=dt-dt_temp_temp;

!As part of the RK second order algorithm, we move first the particle or detector using half of the current time step (dt_var_temp)

         this_det%position=(vel+old_vel)/2*dt_var_temp/2+old_pos

         this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

         old_index_minloc=index_next_face

!Next we check if the detector ends up having again negative local coord. through the face it came out of the domain before.
!If so, we check then if another of the local coords. is also negative and in that case it has gone out through the element via another face or boundary, but if no other local coord is also negative, then the detector is ending in the same place as before, that was outside the domain and if this is the case we make it static.

         if (this_det%local_coords(index_next_face)<0.) then

            cont_ccc=cont_ccc+1
           
            do k=1, size(this_det%local_coords)

               if (k /= old_index_minloc) then

                 cont_ddd=cont_ddd+1

                 if (this_det%local_coords(k)<0.) then

                 cont_eee=cont_eee+1

                 index_next_face=k

                 end if

               end if

            end do

            if (index_next_face == old_index_minloc) then

               cont_fff=cont_fff+1

               this_det%position=old_pos
               this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)
               ewrite(-1,*) 'Warning: converting detector to static'
               this_det%type=STATIC_DETECTOR
        
            end if     
   
         end if 
    
         this_det%position=old_pos
         this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

         deallocate(element_next_to_boundary)

       else

         !If I don't own the previous element where the detector was before it left the domain, it means it was an halo element
         !of another processor, and hence, the detector is sent to the send list associated with that processor.

         !Node in linked_list that contains this detector is removed and incorporated into the send_list 
         !corresponding to the neighbouring processor

         this_det%element=current_element 

         if (halo_level /= 0) then
  
         univ_ele = halo_universal_number(ele_halo, current_element)

         end if

         processor_number=element_owner(vfield,this_det%element)

         this_det%position=old_pos
         this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

         dt_temp=dt_temp+this_det%dt
         this_det%dt=dt-dt_temp

         list_neigh_processor=fetch(ihash,processor_number)

         node_to_send => this_det

         this_det => this_det%next

         call move(detector_list,node_to_send,send_list_array(list_neigh_processor))

      end if

    end if

  end subroutine check_if_det_gone_through_domain_boundary 

  subroutine scenario_det_gone_through_element_vortex(this_det, xfield, dt, dt_temp, keep_value_this_det_dt, old_pos, &
                                          vel, old_vel, vfield, old_vfield, previous_element, index_next_face,cont_check,cont_check_a,bound_elem_iteration)
  !!< Subroutine that makes sure the Lagrangian detector ends up on one of its boundary faces (within tolerance of
  !+/-10.0e-8) up to where it has reached using one of the smaller time steps (sutime step). 
  !Once the detector is on the boundary (within tolerance of +/-10.0e-8), in the subroutine that calls this one, the   
  !detector is asigned to the neighbouring element through that face and again calling this subroutine, the detector is 
  !moved using another smaller time step until the boundary of this new element, always following the direction of the flow. 
  !These steps are repeated until the sum of all the smaller time steps used to go from elemenet to element is equal to the 
  !total time step.
!    type(state_type), dimension(:), intent(in) :: state

    type(detector_type), pointer :: this_det
    type(vector_field), intent(inout) :: xfield
    real, intent(in) :: dt, keep_value_this_det_dt
    real, intent(inout) :: dt_temp
    real,  dimension(:), intent(inout) :: old_pos, vel, old_vel
    type(vector_field), pointer :: vfield, old_vfield
    integer, intent(in) :: previous_element, cont_check, cont_check_a
    integer, intent(inout) :: index_next_face, bound_elem_iteration

    real :: maxvalue
    integer :: i, k, cont_d, index_temp_maxvalue, cont_loop, node, current_element_number, number_of_elem, &
               cont_elem_neg, cont_static_det, det_inside_an_ele

    integer, dimension(:), pointer :: nodes, elements
    type(csr_sparsity), pointer :: nelist

    ewrite(-1,*) 'WARNING, BISECTION METHOD NOT RECOMMENDED!'
    if ((cont_check /= 0).and.(cont_check_a /= 0)) then

       ewrite(1,*) "Inside scenario_det_gone_through_element_vortex subroutine"

       this_det%position=old_pos
       this_det%local_coords=local_coords(xfield,this_det%element,this_det%position) 
       vel = detector_value(vfield, this_det)
       old_vel = detector_value(old_vfield, this_det)
       
       nelist => extract_nelist(xfield)

       nodes => ele_nodes(xfield, this_det%element) !pointer to the nodes of a given element number this_det%element

       maxvalue=0.0

       do k=1, size(this_det%local_coords)

          if (this_det%local_coords(k)>maxvalue) then

             maxvalue=this_det%local_coords(k)

             index_temp_maxvalue=k

          end if

       end do  

       current_element_number=this_det%element

       node=nodes(index_temp_maxvalue) !node number of the node that is opposite to the face with respect to which the local 
                                       !coodinate has the highest value
       elements => row_m_ptr(nelist, node) 
 
       !pointer to the row number corresponding to the node number that returns the index of 
       !the columns that are different than zero. These indexes are the numbers of the elements 
       !that share/contain that node number
  
       number_of_elem=size(elements)

       
       !In the loop below, it is checked if the detector belongs to any of the elements that share the node
      
       do i=1, size(elements)
        
          !if (elements(i)==current_element_number) cycle 

          !If one of the elements is negative

          cont_elem_neg=0

          if ((elements(i)<0.0)) then

             cont_elem_neg=cont_elem_neg+1        

          end if

          if (cont_elem_neg/=0) cycle

          this_det%element=elements(i)

          this_det%dt=keep_value_this_det_dt

          cont_d=0

          do
  
             this_det%dt=this_det%dt/2

             this_det%position=(vel+old_vel)/2*this_det%dt+old_pos

             this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

             cont_loop=0

             do 

                if (all(this_det%local_coords>-10.0e-8)) exit

                this_det%dt=this_det%dt/2
       
                call update_detector_position_bisect(this_det, xfield, old_pos, vel, old_vel)

                cont_loop=cont_loop+1

                if ((cont_loop>1.0e6).or.(this_det%dt<1.0e-5*keep_value_this_det_dt)) exit

             end do

             if ((cont_loop>1.0e6).or.(this_det%dt<1.0e-5*keep_value_this_det_dt)) exit

             vel =  detector_value(vfield, this_det)
             old_vel = detector_value(old_vfield, this_det)
             this_det%dt=this_det%dt*2
             this_det%position=(vel+old_vel)/2*this_det%dt+old_pos

             this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

             if (all(this_det%local_coords>-10.0e-8)) exit

             cont_d=cont_d+1

             this_det%dt=this_det%dt/2
             this_det%position=old_pos
             this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)
             vel =  detector_value(vfield, this_det)
             old_vel = detector_value(old_vfield, this_det)

          end do

!          if (all(this_det%local_coords>-10.0e-8)) exit

          if ((all(this_det%local_coords>-10.0e-8)).and.(this_det%dt>1.0e-5*keep_value_this_det_dt)) exit

       end do

       det_inside_an_ele=0.0

!       if  (all(this_det%local_coords>-10.0e-8)) then

       if ((all(this_det%local_coords>-10.0e-8)).and.(this_det%dt>1.0e-5*keep_value_this_det_dt)) then

           det_inside_an_ele=det_inside_an_ele+1

       end if

       if ((i==size(elements)).and.(det_inside_an_ele==0.0)) then
          ewrite(-1,*) 'Warning: converting detector to static'
           this_det%type=STATIC_DETECTOR

       end if

       cont_static_det=0
    
       if (this_det%type==STATIC_DETECTOR) then

          cont_static_det=cont_static_det+1
          
          this_det%position=old_pos
          this_det%element=previous_element
          this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

       end if

       if (this_det%type==STATIC_DETECTOR) return

       call placing_det_in_boundary_between_elem(this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, index_next_face, bound_elem_iteration)
  
    end if

  end subroutine scenario_det_gone_through_element_vortex

  subroutine placing_det_in_boundary_between_elem(this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, index_next_face, bound_elem_iteration)

    type(detector_type), pointer :: this_det
    type(vector_field), intent(inout) :: xfield
    real, intent(in) :: dt
    real, intent(inout) :: dt_temp
    real,  dimension(:), intent(inout) :: old_pos, vel, old_vel
    integer, intent(inout) :: index_next_face, bound_elem_iteration
    
    real :: dt_var_check, dt_temp_check, dt_a, dt_b
    integer :: cont_p, cont_ready_out

    ewrite(1,*) "Inside placing_det_in_boundary_between_elem subroutine"
    ewrite(-1,*) 'WARNING, BISECTION METHOD NOT RECOMMENDED!'
    bound_elem_iteration=0

       do

       !dt_a is the value of the bisected dt with which the detector falls inside the element for first time

          dt_temp_check=dt_temp+this_det%dt
          dt_var_check=dt-dt_temp_check

          index_next_face=minloc(this_det%local_coords, dim=1)

          cont_ready_out=0

          cont_p=0

          !Next it is checked if the detector is on one of the element faces (within tolerance +/-10.0e-8). 

          !If not, it needs to be moved further in the direction of the flow until it hits the element face.  
  
          !Other conditions where placed before regarding if the ratio between the other local coordinates after and before moving 
          !the detector backwards was less than 1, i.e., the other local coords are generally increasing as we aproach the right 
          !face (with respect to which the local coord is minimum).

          !However there was some ambiguity when using those other conditions as well since for a few geometric cases it was not 
          !satisfied and some detectors were converted into static in the middle of the domain. Hence, they were removed. If the 
          !detector leaves through the wrong face, it is taken care of in other parts of this subroutine, making the detector
          !to go back to the previous element.

          if ((minval(this_det%local_coords)<10.0e-8).and.(minval(this_det%local_coords)>-10.0e-8)) then

                cont_ready_out=cont_ready_out+1
           
          end if 


          if (cont_ready_out /= 0)  index_next_face=minloc(this_det%local_coords, dim=1)
   
          if (cont_ready_out /= 0)  exit   

          !If dt_var_check is practically zero means that all the smaller dt_var used add up to dt so it is the final position
          !of the detector for that dt

          if (dt_var_check<1.0e-3*dt) exit 

          if (bound_elem_iteration>=33554432) exit

          dt_a=this_det%dt
          dt_b=2*dt_a

          !If not yet on one of the faces of the element (within tolerance), it needs to be moved further 
          !in the direction of the flow until it hits the element face.

          if (dt_b-dt_a<1.0e-8) exit

          do

             !This loop terminates when the detector is on one of the element faces (within tolerance +/-10.0e-8). 

             !When that happens, the detector is moved backwards a bit, to check that the ratio of the local coordinates 
             !associated to the min local coordinate before moving backwards the detector, is not equal to 1. Same explanation 
             !as before
         
             if ((minval(this_det%local_coords)<10.0e-8).and.(minval(this_det%local_coords)>-10.0e-8)) then

                 cont_ready_out=cont_ready_out+1 

             end if 

             if (cont_ready_out /= 0)  exit  

             if (dt_var_check<1.0e-3*dt) exit 

             if (dt_b-dt_a<1.0e-8) exit

             if (bound_elem_iteration>=33554432) exit

             call iterating_for_det_in_bound_elem(this_det, xfield, old_pos, vel, old_vel, dt_a, dt_b, bound_elem_iteration)

          end do

          cont_p=cont_p+1 

       end do

  end subroutine placing_det_in_boundary_between_elem

  subroutine iterating_for_det_in_bound_elem(this_det, xfield, old_pos, vel, old_vel, dt_a, dt_b, bound_elem_iteration)

    type(detector_type), pointer :: this_det
    type(vector_field), intent(inout) :: xfield
    real, intent(inout) :: dt_a, dt_b
    real,  dimension(:), intent(inout) :: old_pos, vel, old_vel
    integer, intent(inout) :: bound_elem_iteration
    
            
    integer :: cont_a, cont_b, cont_c

       ewrite(1,*) "Inside iterating_for_det_in_bound_elem subroutine"
       ewrite(-1,*) 'WARNING, BISECTION METHOD NOT RECOMMENDED!'       

       cont_a=0
       cont_b=0
       cont_c=0

       this_det%dt=dt_a+((dt_b-dt_a)/2)

       bound_elem_iteration=2

       cont_a=cont_a+1
  
       do 
              
          call update_detector_position_bisect(this_det, xfield, old_pos, vel, old_vel)

          if (all(this_det%local_coords>-10.0e-8)) exit 

          bound_elem_iteration=bound_elem_iteration*2

          this_det%dt=dt_a+((dt_b-dt_a)/bound_elem_iteration)

          cont_b=cont_b+1

          !if i becomes too big is because no matter how much dt_var is reduced the detector is not found inside the element

          if (bound_elem_iteration>=33554432) exit

       end do

       ewrite(1,*) "cont_b:", cont_b

       ewrite(1,*) "bound_elem_iteration:", bound_elem_iteration

       ewrite(1,*) "dt_a:", dt_a

       ewrite(1,*) "dt_b:", dt_b

       dt_b=dt_a+((dt_b-dt_a)/(bound_elem_iteration/2))
       dt_a=this_det%dt

       ewrite(1,*) "dt_a:", dt_a

       ewrite(1,*) "dt_b:", dt_b

       cont_c=cont_c+1  

  end subroutine iterating_for_det_in_bound_elem

  subroutine scenario_det_gone_through_wrong_face(this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, index_next_face,cont_check,cont_check_a,bound_elem_iteration, index_minloc_current)

    type(detector_type), pointer :: this_det
    type(vector_field), intent(inout) :: xfield
    real, intent(in) :: dt
    real, intent(inout) :: dt_temp
    real,  dimension(:), intent(inout) :: old_pos, vel, old_vel
    integer, intent(in) :: cont_check, cont_check_a
    integer, intent(inout) :: index_next_face, bound_elem_iteration, index_minloc_current

    real :: dt_var_check, dt_temp_check, dt_a, dt_b, minvalue
    integer :: k, cont_p, &
               cont_apq, index_temp, cont_ready

    ewrite(-1,*) 'WARNING, BISECTION METHOD NOT RECOMMENDED!'
    if ((cont_check /= 0).and.(cont_check_a==0)) then

       ewrite(1,*) "Inside scenario_det_gone_through_wrong_face subroutine"

       bound_elem_iteration=0

       cont_ready=0

       do

          dt_temp_check=dt_temp+this_det%dt
          dt_var_check=dt-dt_temp_check
    
          !This loop terminates when the detector is on one of the element faces (within tolerance +/-10.0e-8) and the min local
          !coordinate occurs with respect to a different face than before when the detector left through the wrong face. 

          !if ((minval(this_det%local_coords)<10.0e-8).and.(minval(this_det%local_coords)>-10.0e-8).and.(minloc(this_det%local_coords, dim=1)/=index_minloc_current)) exit 

          !In the next if, if the min local coords still occurs in the same index as before, then it is checked if 
          !there is another local coord that is also within tolerance and the face associated with it is the new face through 
          !which the detector will leave the element.


           if (cont_ready /= 0) exit

           if (dt_var_check<1.0e-3*dt) exit 

          !The check below means we have lost the detector at some point. We dont know what has happened so we exit the loop and at the end of 
          !the subroutine we convert the detector into a static one.

          if (bound_elem_iteration>=33554432) exit

          if ((minval(this_det%local_coords)<10.0e-8).and.(minval(this_det%local_coords)>-10.0e-8).and.(minloc(this_det%local_coords, dim=1)/=index_minloc_current)) exit 

          index_next_face=minloc(this_det%local_coords, dim=1)

          cont_apq=0

          dt_a=this_det%dt
          dt_b=2*dt_a

          !If the detector is not on one of the element faces (within tolerance +/-10.0e-8), then it needs to keep moving in the
          !direction of the flow until it hits a face.

!          cont_cont=0

          do

          !As before this loop terminates when the detector is on one of the element faces (within tolerance +/-10.0e-8) 
          !and the min local coordinate occurs with respect to a different face than before when the detector left through 
          !the wrong face. Same explanations as before apply
 
             cont_ready=0

          !This loop terminates when the detector is on one of the element faces (within tolerance +/-10.0e-8) and the min local
          !coordinate occurs with respect to a different face than before when the detector left through the wrong face. 

             if ((minval(this_det%local_coords)<10.0e-8).and.(minval(this_det%local_coords)>-10.0e-8).and.(minloc(this_det%local_coords, dim=1)/=index_minloc_current)) exit 


          !In the next if, if the min local coords still occurs in the same index as before, then it is checked if 
          !there is another local coord that is also within tolerance and the face associated with it is the new face through 
          !which the detector will leave the element.

             if ((minval(this_det%local_coords)<10.0e-8).and.(minval(this_det%local_coords)>-10.0e-8).and.(minloc(this_det%local_coords, dim=1)==index_minloc_current)) then

                 minvalue=1000.0

                 do k=1, size(this_det%local_coords)

                    if (k /= index_minloc_current) then

                        if (this_det%local_coords(k)<minvalue) then

                            minvalue=this_det%local_coords(k)

                            index_temp=k

                            ewrite(1,*) "index_temp", index_temp

                        end if

                    end if

                 end do  

                 cont_apq=cont_apq+1

                 cont_ready=0

                 if ((this_det%local_coords(index_temp)<10.0e-8).and.(this_det%local_coords(index_temp)>-10.0e-8))  then

                    cont_ready=cont_ready+1

                    index_next_face=index_temp
              
                 end if

                 ewrite(1,*) "index_temp", index_temp

                 ewrite(1,*) "index_next_face", index_next_face

             end if
 
             if (cont_ready /= 0) exit
          
             dt_temp_check=dt_temp+this_det%dt
             dt_var_check=dt-dt_temp_check
               
             if (dt_var_check<1.0e-3*dt) exit 

             if (dt_b-dt_a<1.0e-8) exit

             if (bound_elem_iteration>=33554432) exit

             call iterating_for_det_in_bound_elem(this_det, xfield, old_pos, vel, old_vel, dt_a, dt_b, bound_elem_iteration)

          end do

          cont_p=cont_p+1 

       end do

    end if
  
  end subroutine scenario_det_gone_through_wrong_face

  subroutine scenario_det_gone_through_right_face(this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, index_next_face, cont_check,bound_elem_iteration)

    type(detector_type), pointer :: this_det
    type(vector_field), intent(inout) :: xfield
    real, intent(in) :: dt
    real, intent(inout) :: dt_temp
    real,  dimension(:), intent(inout) :: old_pos, vel, old_vel
!    type(vector_field), pointer :: vfield, old_vfield
    integer, intent(in) :: cont_check
    integer, intent(inout) :: index_next_face, bound_elem_iteration

    ewrite(-1,*) 'WARNING, BISECTION METHOD NOT RECOMMENDED!'    
    if (cont_check == 0) then

       ewrite(1,*) "Inside scenario_det_gone_through_right_face subroutine"

       call placing_det_in_boundary_between_elem(this_det, xfield, dt, dt_temp, old_pos, &
                                         vel, old_vel, index_next_face, bound_elem_iteration)
      
    end if 

  end subroutine scenario_det_gone_through_right_face

  subroutine update_detector_position_bisect(this_det, xfield, old_pos, vel, old_vel)
    !!< Moves the detector in the direction of the flow from an initial position (old_pos) during a small time step equal to dt_var 
    !(a bisection of the total time step) and updates the local coordinates 
!    type(state_type), dimension(:), intent(in) :: state

    type(detector_type), pointer :: this_det
    type(vector_field), intent(inout) :: xfield
    real,  dimension(:), intent(in) :: old_pos, vel, old_vel

       this_det%position=((vel+old_vel)/2.0)*this_det%dt+old_pos

       this_det%local_coords=local_coords(xfield,this_det%element,this_det%position)

  end subroutine update_detector_position_bisect

end module detector_move_bisection
