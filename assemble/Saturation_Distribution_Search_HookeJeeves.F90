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

module saturation_distribution_search_hookejeeves
  !!< This module contains an algorithm for searching through a number of
  !!< permutations of saturation distribution (assuming multiphase flow) and
  !!< calculating the electrical potential resulting from it.
  !!< The immediate application for this is as a quasi-inversion method for
  !!< determining water distribution within a reservoir during production
  
  use diagnostic_fields_wrapper_new
  use fields
  use fields_base, only: ele_nodes
  use field_options
  use fldebug
  use global_parameters, only: OPTION_PATH_LEN
  use spontaneous_potentials
  use spud
  use state_module
  use write_state_module

  implicit none

  private

  public :: search_saturations_hookejeeves

  contains

  subroutine search_saturations_hookejeeves(state, state_no)

    ! declare interface variables
    type(state_type), dimension(:), intent(inout) :: state
    type(state_type), dimension(:), pointer :: working_state
    integer, intent(in) :: state_no

    ! declare working variables for this subrtn
    character(len=50) :: positions
    character(len=OPTION_PATH_LEN) :: option_buffer
    real, allocatable :: base_point(:), widths(:)
    real :: step_length, nought, search_min, search_max
    type(scalar_field) :: saturation
    type(scalar_field), pointer :: Cv
    type(vector_field), pointer :: coordinates

    integer :: sections, stat, j
    logical :: is_vwell, do_PS

    ewrite(3,*) 'In search_saturations'

    option_buffer = '/material_phase['//int2str(state_no-1)//']/electrical_properties/Saturation_Distribution_Search/'

    if (have_option(trim(option_buffer)//'search_criteria_vertical_well')) then
       is_vwell=.true.
       option_buffer = '/material_phase['//int2str(state_no-1)//']/electrical_properties/Saturation_Distribution_Search/search_criteria_vertical_well/'
    elseif (have_option(trim(option_buffer)//'search_criteria_horizontal_well')) then
       is_vwell=.false.
       option_buffer = '/material_phase['//int2str(state_no-1)//']/electrical_properties/Saturation_Distribution_Search/search_criteria_horizontal_well/'
    else
       FLExit('Cannot find well geometry')
    endif
    
    !! Optimisation based on Hooke and Jeeves algorithm (Basic Optimisation Methods, BD Bunday, 1984)
    !set step_length
    call get_option(trim(option_buffer)//'initial_step_length', step_length)
    !set dimension of coordinates (number of variables to be optimised)
    call get_option(trim(option_buffer)//'sections', sections)
    allocate(base_point(sections))
    allocate(widths(sections))
    
    if (is_vwell) then
       call get_option(trim(option_buffer)//'y_min',search_min)
       call get_option(trim(option_buffer)//'y_max',search_max)
    else
       call get_option(trim(option_buffer)//'z_min',search_min)
       call get_option(trim(option_buffer)//'z_max',search_max)
    endif
    call get_option(trim(option_buffer)//'width',widths(1))
    
    do j=1,sections
       widths(j)=widths(1)
    end do
    ewrite(3,*) 'Initial widths',widths(1)

    if (have_option(trim(option_buffer)//'initial_base_point')) then
       call get_option(trim(option_buffer)//'initial_base_point',base_point)
       ewrite(3,*) 'Initial base point', base_point
    else
       do j=1,sections
          base_point(j)=(search_max-search_min)/2
       end do
    endif
    do_PS=.false.
    nought=0.0
    call set_base_point(state, state_no, do_PS, base_point, widths, step_length, base_point, nought, &
                        search_min, search_max)
    
    deallocate(base_point)
    deallocate(widths)
   
  end subroutine search_saturations_hookejeeves
    
  recursive subroutine set_base_point(state, state_no, do_PS, base_point, widths, step_length, previous_base_point, &
                                      current_error, search_min, search_max)
  
     type(state_type), dimension(:), intent(inout) :: state
     real, dimension(:) :: base_point, widths
     real, dimension(size(base_point)) :: copy_base_point
     real :: step_length, minimum_step_length, new_error, search_min, search_max
     logical :: do_PS, is_improved, finished
     ! optional argument previous_base_point allows us to make pattern step
     real, dimension(:) :: previous_base_point
     ! optional argument current_error records error at this base_point
     real :: current_error, error
!     integer, save :: dump_no
     integer, intent(in) :: state_no
     
     ewrite(3,*) 'setting base point', base_point
     
     ! Set minimum_step_length for now
     minimum_step_length = 5.0
  
     if (do_PS) then
!        call write_state(dump_no, state)
        call make_pattern_step(state, state_no, base_point, previous_base_point, widths, step_length, minimum_step_length, &
                               current_error, search_min, search_max)
     else
        ! This is the first time through so need to calculate a first error
        ! WRONG WRONG WRONG! WE CAN GET HERE IF WE RESET TO AN EARLIER BASEPOINT
        ! IN WHICH CASE WE NEED TO BE PASSING BACK THE SMALLEST ERROR
        call run_model(state, state_no, base_point, widths, step_length, error, is_improved)
!        dump_no=1
!        call write_state(dump_no, state)
        ! Copy current base point for pattern search in case base_point is updated by exploration
        copy_base_point=base_point
        
        ! THE ERROR WE PASS IN HERE HAS TO BE THE SMALLEST ERROR
        call make_exploration(state, state_no, base_point, widths, step_length, is_improved, error, search_min, search_max)
        if(.not.is_improved) then
           call reduce_step_length(state, state_no, base_point, base_point, widths, step_length, minimum_step_length, &
                                   error, finished, search_min, search_max)
           if(.not.finished) then
              call set_base_point(state, state_no, do_PS, base_point, widths, step_length, base_point, error, &
                                  search_min, search_max)
           else
              ewrite(3,*) 'Finishing search'
              return ! FINISHED! PRINT RESULTS ETC ETC
           endif
        else
           call set_base_point(state, state_no, .TRUE., base_point, widths, step_length, copy_base_point, error, &
                               search_min, search_max)
        endif
     endif
  
  end subroutine set_base_point
  
  subroutine make_pattern_step(state, state_no, base_point, previous_base_point, widths, step_length, minimum_step_length, &
                               current_error, search_min, search_max)
  
     type(state_type), dimension(:), intent(inout) :: state
     real, dimension(:) :: base_point, previous_base_point, widths
     real, dimension(size(base_point)) :: pattern_point
     real :: step_length, current_error, minimum_step_length, search_min, search_max
     integer :: i
     integer, intent(in) :: state_no
     logical :: is_improved, finished
     
     ewrite(3,*) 'making pattern step'
  
     ! calculate pattern point based on previous success
     do i=1,size(base_point)
        pattern_point(i) = previous_base_point(i) + 2*(base_point(i)-previous_base_point(i))
        pattern_point(i) = max(pattern_point(i), search_min)
        pattern_point(i) = min(pattern_point(i), search_max)
     end do
     ewrite(3,*) 'base_point, pattern_point', base_point, pattern_point
     
     ! current error records error at base point we're going to try to leap over
     
     ! explore around this new point
     call make_exploration(state, state_no, pattern_point, widths, step_length, is_improved, current_error, &
                           search_min, search_max)
     finished=.false.
     if (.not.is_improved) then
        call reduce_step_length(state, state_no, base_point, pattern_point, widths, step_length, minimum_step_length, &
                                current_error, finished, search_min, search_max)
     else
        call set_base_point(state, state_no, .TRUE., pattern_point, widths, step_length, base_point, current_error, &
                            search_min, search_max)
     endif
  
  end subroutine make_pattern_step
  
  subroutine make_exploration(state, state_no, base_point, widths, step_length, global_improved, error, &
                              search_min, search_max)
     
     type(state_type), dimension(:), intent(inout) :: state
     real, dimension(:) :: base_point, widths
     real, allocatable, dimension(:) :: base_point_in
     real :: step_length, error, new_error, search_min, search_max
     logical :: is_improved, global_improved
     integer :: i
     integer, intent(in) :: state_no
     
     ewrite(3,*) 'making exploration, step ', step_length
     
     ! return new base point in base_point if it's been improved
     ! ie overwrite
  
     ! cycle through co-ordinates:
     is_improved=.false.
     global_improved=.false.
     
     allocate(base_point_in(size(base_point)))
     base_point_in=base_point
     do i=1,size(base_point)
        ! increase by step length
        base_point(i)=base_point(i)+step_length
        base_point(i)=max(base_point(i), search_min)
        base_point(i)=min(base_point(i), search_max)
        call run_model(state, state_no, base_point, widths, step_length, new_error, is_improved)
!        ewrite(3,*) 'here 1, new_error, error', new_error, error
        if (.not.is_improved) then
           !decrease by step length
           base_point(i)=base_point_in(i)-step_length
           base_point(i)=max(base_point(i), search_min)
           base_point(i)=min(base_point(i), search_max)
           call run_model(state, state_no, base_point, widths, step_length, new_error, is_improved)
!           ewrite(3,*) 'here 2, new_error, error', new_error, error
           if (.not.is_improved) then
              ! keep original
              base_point(i)=base_point_in(i)
           else
              ! set coordinate to new value
              ewrite(3,*) 'here 3'
              ewrite(3,*) 'Error at', base_point(i), new_error, 'less than previous error', error
              error=new_error
              global_improved=.true.
           endif
        else
           ! set coordinate to new value
           ewrite(3,*) 'here 4'
           ewrite(3,*) 'Error at', base_point(i), new_error, 'less than previous error', error
           error=new_error
           global_improved=.true.
        endif
     end do
     deallocate(base_point_in)
     
     ! Vary width of each layer to find a better fit
     call improve_width(state, state_no, base_point, widths, step_length, error, is_improved)
     global_improved=global_improved.or.is_improved
  
  end subroutine make_exploration
  
  subroutine run_model(state, state_no, point, widths, step_length, error, is_improved)
  
     type(state_type), dimension(:), intent(inout) :: state
     type(scalar_field) :: saturation, temp_saturation
     type(scalar_field), pointer :: Cv, electrical_potential
     type(vector_field), pointer :: coordinates
     real, dimension(:) :: point, widths
     real, allocatable, dimension(:,:), save :: target_potential, used_list
     real, allocatable, dimension(:) :: model_potential
     real :: step_length, error
     real :: x_min, x_max, y_min, y_max, z_min, z_max
     real :: bh_x, bh_y, bh_z, err
     real, save :: smallest_error
     logical :: is_vwell, used, is_improved
     integer :: stat, i, j
     integer, intent(in) :: state_no
     integer, save :: have_target, samples, dump_no
     integer, allocatable, dimension(:), save :: node_list
     character(len=OPTION_PATH_LEN) :: option_buffer, target_filename
     character(len=25) :: fstr
     
     ewrite(3,*) 'running model, point, widths:',point, widths
     
     ! Check to see if we've run this point width combo before
     used = .false.
     if (.not.allocated(used_list)) then
        allocate(used_list(1000,size(point)+size(widths)+1))
        used_list=0.0
     endif
     call check_already_used(point, widths, used_list, used, error)
     if (used) then
        ewrite(3,*) 'Used this point, width combination before'
        is_improved=.false.
        return
     endif
     
     ! Get all the limits etc relevant - do this fresh each time from flml
     option_buffer = '/material_phase['//int2str(state_no-1)//']/electrical_properties/Saturation_Distribution_Search/'

     if (have_option(trim(option_buffer)//'search_criteria_vertical_well')) then
        is_vwell=.true.
        option_buffer = '/material_phase['//int2str(state_no-1)//']/electrical_properties/Saturation_Distribution_Search/search_criteria_vertical_well/'
     elseif (have_option(trim(option_buffer)//'search_criteria_horizontal_well')) then
        is_vwell=.false.
        option_buffer = '/material_phase['//int2str(state_no-1)//']/electrical_properties/Saturation_Distribution_Search/search_criteria_horizontal_well/'
     else
        FLExit('Cannot find well geometry')
     endif
     ewrite(3,*) 'Is_vwell: ', is_vwell

     call get_option(trim(option_buffer)//'x_min', x_min)
     call get_option(trim(option_buffer)//'x_max', x_max)
     call get_option(trim(option_buffer)//'y_min', y_min)
     call get_option(trim(option_buffer)//'y_max', y_max)
     call get_option(trim(option_buffer)//'z_min', z_min)
     call get_option(trim(option_buffer)//'z_max', z_max)
     call get_option(trim(option_buffer)//'borehole_x',bh_x)
     
     ! Set saturations to point
     saturation = extract_scalar_field(state, "PhaseVolumeFraction", stat)
     call allocate(temp_saturation, saturation%mesh, name="OriginalSaturation")
     ewrite(3,*) 'copied to temporary saturation field'
     call set(temp_saturation,saturation)
     coordinates=>extract_vector_field(state, "Coordinate")
     
     call set_saturations(saturation, point, widths, bh_x, is_vwell, &
                          x_min, x_max, y_min, y_max, z_min, z_max, &
                          coordinates)

     ! get target potential curve if it hasn't been got before
     if (.not.allocated(target_potential)) then
        samples=100
        have_target=8
        call get_option(trim(option_buffer)//'target_filename', target_filename)
        allocate(target_potential(samples,2))
        open(unit = 910, file = trim(target_filename), action = "read")
        i=1
        err=1.0
        do while(i<=1000)
           read(910,*,end=69) target_potential(i,1:2)
!           ewrite(3,*) 'target_potential line',i,target_potential(i,1), target_potential(i,2)
           if (is_vwell) then
              if ((target_potential(i,1).ge.z_min-err).and.(target_potential(i,1).le.z_max+err)) then
                 i=i+1
              endif
           else
              if ((target_potential(i,1).ge.y_min-err).and.(target_potential(i,1).le.y_max+err)) then
                 i=i+1
              endif
           endif
        end do
69      ewrite(3,*) 'target file read in'
        samples=i-1
        close(910)
        allocate(node_list(samples))
     ! Find nodes corresponding to target potential col 1
        do i=1,samples
!           ewrite(3,*) 'i: ',i
           if (is_vwell) then
              call get_option(trim(option_buffer)//'borehole_y',bh_y)
              do j=1,node_count(coordinates)
                 if ((abs(coordinates%val(X_)%ptr(j)-bh_x)<err).and.(abs(coordinates%val(Y_)%ptr(j)-bh_y)<err).and.&
                    (abs(coordinates%val(Z_)%ptr(j)-target_potential(i,1))<err)) then
                    node_list(i)=j
                 end if
              end do
           else
              call get_option(trim(option_buffer)//'borehole_x',bh_x)
              call get_option(trim(option_buffer)//'borehole_z',bh_z)
              do j=1,node_count(coordinates)
                 if ((abs(coordinates%val(X_)%ptr(j)-bh_x)<err).and.(abs(coordinates%val(Z_)%ptr(j)-bh_z)<err).and.&
                    (abs(coordinates%val(Y_)%ptr(j)-target_potential(i,1))<err)) then
                    node_list(i)=j
                 endif
              end do
           endif
        end do
     end if
!     ewrite(3,*) 'node_list: ', node_list
     
     ! Update coupling coefficient
     Cv => extract_scalar_field(state(state_no), 'Electrokinetic[0]', stat=stat)
     if (stat/=0) then
        FLExit('Did not find an Electrokinetic coupling coefficient scalar field.')
     end if

     call calculate_diagnostic_variable(state, state_no, Cv) 
     ! Get streaming potential
     call calculate_electrical_potential(state(state_no), state_no)
     ewrite(3,*) 'out of calculate_electrical_potential'
     
     ! get potential curve from new model
     electrical_potential=>extract_scalar_field(state(state_no), "ElectricalPotential", stat)
     allocate(model_potential(samples))
     
     ewrite(3,*) 'getting potentials from latest model'
     do i=1,samples
        model_potential(i)=electrical_potential%val(node_list(i))
     end do

     ! Compare and calculate error
!     ewrite(3,*) 'error=',error
!     ewrite(3,*) 'target_potential, model_potential',target_potential(1:samples,2), model_potential
     call curve_error(target_potential(1:samples,2), model_potential, error)
     
     deallocate(model_potential)
     
     ! add point to used_list
     call add_to_used_list(point, widths, error, used_list)
     
     ewrite(3,*) 'Smallest error, error',smallest_error, error
     is_improved=.false.
     
     dump_no=max(1,dump_no)
     if (dump_no.eq.1) then
        ewrite(3,*) 'here 5'
        smallest_error=error
        call write_state(dump_no, state)
        ewrite(3,*) 'point, widths, error: ', point, widths, error
        open(911, file = 'output.data', action="write", position="append")
        write(911,10) point, widths, error
     10 FORMAT(16(F9.4,2x),16(F5.1,2x),F10.6)
        close(911)
        is_improved=.true.
     elseif (error<smallest_error) then
        ewrite(3,*) 'here 6'
        smallest_error=error
        call write_state(dump_no, state)
        ewrite(3,*) 'point, widths, error: ', point, widths, error
        open(913, file = 'output.data', action="write", position="append")
        write(913,12) point, widths, error
     12 FORMAT(16(F9.4,2x),16(F5.1,2x),F10.6)
        close(913)
        is_improved=.true.
     else
        error=smallest_error
        is_improved=.false.
     endif
     
     ! Set saturation back to the initial condition ready for next time
     call set(saturation,temp_saturation)
     call deallocate(temp_saturation)

  end subroutine run_model
  
  subroutine set_saturations(saturation, point, widths, borehole_x, is_vwell, &
                             x_min, x_max, y_min, y_max, z_min, z_max, &
                             coordinates)
    
    type(scalar_field) :: saturation
    type(vector_field), intent(in) :: coordinates
    real, dimension (:) :: point, widths
    real, intent(in) :: x_min, x_max, y_min, y_max, z_min, z_max
    real :: along_step, borehole_x
    real :: zone_along_min, err
    real, dimension(coordinates%dim, coordinates%mesh%shape%loc) :: coords
    integer :: ele, node, in_zone_counter
    integer :: i, n, coord_along, coord_perp
    logical :: cell_centred, in_zone, is_vwell
    
    !!< Set saturations in each of the zones to match the positions string
    ewrite(3,*) 'In set_saturations'
    
    ! determine scalar field continuity
    if (saturation%mesh%continuity==-1) then
      ! mesh is discontinuous and saturation is cell-centred
      n = ele_count(saturation)
      cell_centred=.true.
    else
      ! mesh is continuous and saturation is node-centred
      n = node_count(saturation)
      cell_centred=.false.
    end if
!    ewrite(3,*) 'Saturation cell-centred:',cell_centred
    
    if (is_vwell) then
       zone_along_min=z_min
       coord_along=3
       coord_perp=2
       along_step=(z_max-z_min)/size(point)
    else
       zone_along_min=y_min
       coord_along=2
       coord_perp=3
       along_step=(y_max-y_min)/size(point)
    endif
    err=0.1
    
    do i=1,size(point)
      if (cell_centred) then
        do ele=1,n
          ! coords is 3x8 array in 3d on a hex mesh
          ! organised so coords(1,4) is x-coordinate of 4th node
          ! and coords(3,2) is z coordinate of 2nd node etc
          coords = ele_val(coordinates, ele)
          in_zone=.false.
          in_zone_counter=0
          if ((coords(1,1).gt.borehole_x-widths(i)/2-err).and.(coords(1,1).lt.borehole_x+widths(i)/2+err).and.&
                (coords(coord_along,1).gt.zone_along_min-err).and.(coords(coord_along,1).lt.zone_along_min+along_step+err).and.&
                (coords(coord_perp,1).lt.point(i)+err)) then
            do node=2,coordinates%mesh%shape%loc
              if ((coords(1,node).gt.borehole_x-widths(i)/2-err).and.(coords(1,node).lt.borehole_x+widths(i)/2+err).and.&
                   (coords(coord_along,node).gt.zone_along_min-err).and.(coords(coord_along,node).lt.zone_along_min+along_step+err).and.&
                   (coords(coord_perp,node).lt.point(i)+err)) then
                in_zone_counter=in_zone_counter+1
              end if
            end do
            ! If they're ALL in the zone
            if (in_zone_counter==coordinates%mesh%shape%loc-1) then
               in_zone=.true.
            endif            
          end if
          if (in_zone) then
            call set(saturation, ele, 0.7)
!            ewrite(3,*) 'coords, sat, swc, sor: ', coords, saturation%val(ele), swc%val(ele), sor%val(ele)
          end if
        end do
!      else
 !       do node=1,n
 !       end do
      end if
      ! move limits to next zone:
      zone_along_min=zone_along_min+along_step
    end do
    
  end subroutine set_saturations
  
  recursive subroutine reduce_step_length(state, state_no, base_point, pattern_point, widths, step_length, minimum_step_length, &
                                          current_error, finished, search_min, search_max)
  
     type(state_type), dimension(:), intent(inout) :: state
     real, dimension(:) :: base_point, pattern_point, widths
     real :: step_length, minimum_step_length, current_error, search_min, search_max
     logical :: finished, is_improved, from_pattern_step
     integer :: i
     integer, intent(in) :: state_no
     character(len=12) :: fstr
     
     ewrite(3,*) 'reducing step length'
  
     ! reduce step length by factor of 2? or to minimum step length
     finished=.false.
     step_length=step_length/2
     if (step_length<minimum_step_length) then
        step_length=minimum_step_length
        finished = .true.
        from_pattern_step=.false.
        do i=1,size(base_point)
           if (base_point(i).ne.pattern_point(i)) then
              from_pattern_step=.true.
           endif
        end do
        if (from_pattern_step) then
           ! ignore pattern_point and return to last base_point
           ewrite(3,*) 'Reached minimum step length, returning to last base point'
           call set_base_point(state, state_no, .FALSE., base_point, widths, step_length, base_point, &
                               current_error, search_min, search_max)
        else
           ewrite(3,*) 'Reached minimum step length, one last try'
           ! Abusing file names here!:
           pattern_point=base_point
           call make_exploration(state, state_no, base_point, widths, step_length, is_improved, &
                                 current_error, search_min, search_max)
           if (is_improved) then
              call set_base_point(state, state_no, .TRUE., base_point, widths, step_length, pattern_point, &
                                  current_error, search_min, search_max)
           else
              ewrite(3,*) 'Failed to improve with minimum step length, finished'
              ewrite(3,*) 'Final base point and RESULT:',base_point
              open(913, file = 'output.data', action="write", position="append")
!              write(fstr,'("(",i4,"f7.2,2x,",i4,"f5.1)")') size(base_point)
!              ewrite(3,*) 'fstr: ',fstr
              write(913,11) 'Result: ',base_point, widths 
           11 FORMAT((A8),32(F9.4,2x))
              close(913)
              return
           endif
        endif
     else
        ewrite(3,*) 'New step length: ', step_length
        call make_exploration(state, state_no, pattern_point, widths, step_length, is_improved, &
                              current_error, search_min, search_max)
        if(is_improved) then
           call set_base_point(state, state_no, .TRUE., pattern_point, widths, step_length, base_point, &
                               current_error, search_min, search_max)
        else
           call reduce_step_length(state, state_no, base_point, pattern_point, widths, step_length, minimum_step_length, &
                                   current_error, finished, search_min, search_max)
        endif
     endif
  
  end subroutine reduce_step_length
  
  subroutine curve_error(curve1, curve2, error)
  
     real, dimension(:), intent(in) :: curve1, curve2
     real, allocatable, dimension(:) :: curve1q, curve2q, ncurve1, ncurve2
     real :: error, maxi1, maxi2, mini1, mini2
     integer :: k
     
     ewrite(3,*) 'calculating error'
     error=0.0
     if (size(curve1).ne.size(curve2)) then
        ewrite(3,*) 'Test curves are not the same length'
     endif
     
     allocate(curve1q(size(curve1)))
     allocate(curve2q(size(curve2)))
     allocate(ncurve1(size(curve1)))
     allocate(ncurve2(size(curve2)))
     
     ! Get rid of solver noise
     curve1q = float(int(curve1 * 1e8 + 0.5))/1e8
     curve2q = float(int(curve2 * 1e8 + 0.5))/1e8
     
     ! Normalise curves to run between 0 and 1
     mini1=minval(curve1q)
     mini2=minval(curve2q)
     curve1q=curve1q-mini1
     curve2q=curve2q-mini2
     
     maxi1=maxval(curve1q)
     maxi2=maxval(curve2q)
     !ewrite(3,*) 'maxi1, maxi2', maxi1, maxi2
     if ((maxi1>0.0).and.(maxi2>0.0)) then
        ncurve1=curve1q/maxi1
        ncurve2=curve2q/maxi2
        do k=1,min(size(ncurve1),size(ncurve2))
           error = error + abs(ncurve1(k)-ncurve2(k))
        end do
     else
        error=1000.0
     endif
!     ewrite(3,*) 'error is ', error
!     ewrite(3,*) 'target',curve1
!     ewrite(3,*) 'curve2',curve2
!     ewrite(3,*) 'targetq',curve1q
!     ewrite(3,*) 'curve2q',curve2q
!     ewrite(3,*) 'maxi1, maxi2', maxi1, maxi2
!     ewrite(3,*) 'targetn', ncurve1
!     ewrite(3,*) 'ncurve2', ncurve2
     
     deallocate(curve1q)
     deallocate(curve2q)
     deallocate(ncurve1)
     deallocate(ncurve2)
  
  end subroutine curve_error
  
  recursive subroutine improve_width(state, state_no, base_point, widths, step_length, error, improved)
  
     type(state_type), dimension(:), intent(inout) :: state
     real, dimension(:) :: base_point, widths
     real, allocatable :: widths_in(:)
     real :: width_step, error, new_error, x_min, x_max, step_length
     logical :: local_improved, improved
     integer :: i
     integer, intent(in) :: state_no
     character(len=OPTION_PATH_LEN) :: option_buffer
     
     ! Get all the limits etc relevant - do this fresh each time from flml
     option_buffer = '/material_phase['//int2str(state_no-1)//']/electrical_properties/Saturation_Distribution_Search/'
     if (have_option(trim(option_buffer)//'search_criteria_vertical_well')) then
        if (have_option(trim(option_buffer)//'search_criteria_vertical_well/search_for_width')) then
           call get_option(trim(option_buffer)//'search_criteria_vertical_well/x_min', x_min)
           call get_option(trim(option_buffer)//'search_criteria_vertical_well/x_max', x_max)
        else
           return
        endif
     elseif (have_option(trim(option_buffer)//'search_criteria_horizontal_well')) then
        if (have_option(trim(option_buffer)//'search_criteria_horizontal_well/search_for_width')) then
           call get_option(trim(option_buffer)//'search_criteria_horizontal_well/x_min', x_min)
           call get_option(trim(option_buffer)//'search_criteria_horizontal_well/x_max', x_max)
        else
           return
        endif
     endif
     
     ewrite(3,*) 'trying to improve the widths'
     width_step = 20.0
     improved=.false.
     
     ! cycle through co-ordinates:
     allocate(widths_in(size(widths)))
     widths_in=widths
     do i=1,size(widths)
        ! increase by step length
        widths(i)=widths(i)+width_step
        widths(i)=max(widths(i), x_min)
        widths(i)=min(widths(i), x_max)
        call run_model(state, state_no, base_point, widths, step_length, new_error, local_improved)
        if (.not.local_improved) then
           !decrease by step length
           widths(i)=widths_in(i)-width_step
           widths(i)=max(widths(i), x_min)
           widths(i)=min(widths(i), x_max)
           call run_model(state, state_no, base_point, widths, step_length, new_error, local_improved)
           if (.not.local_improved) then
              ! keep original
              widths(i)=widths_in(i)
           else
              ! set coordinate to new value
              ewrite(3,*) 'Error with width', widths(i), new_error, 'less than previous error', error
              error=new_error
              improved=.true.
           endif
        else
           ! set coordinate to new value
           ewrite(3,*) 'Error with width', widths(i), new_error, 'less than previous error', error
           error=new_error
           improved=.true.
        endif
     end do
     deallocate(widths_in)
     
  end subroutine improve_width
  
  subroutine check_already_used(point, widths, used_list, used, error)
  
     real, dimension(:) :: point, widths
     real, dimension(:,:) :: used_list
     real :: error
     logical :: used, points_match, widths_match
     integer :: i, j
     
!     ewrite(3,*) 'size(used_list), size(used_list,1), size(used_list,2)', size(used_list), size(used_list,1), size(used_list,2)
     ewrite(3,*) "Checking if we've already used this point width combo"
     
     do i=1,size(used_list,1)
        points_match = .true.
        widths_match = .true.
        if (used_list(i,1)==0.0) then
           return
        else
!           ewrite(3,*) 'checking points'
           do j=1,size(point)
              if (abs(used_list(i,j)-point(j)).gt.0.01) then
                 points_match=.false.
                 exit
!              else
!                 ewrite(3,*) 'points match'
              endif
           end do
!           ewrite(3,*) 'checking widths'
           do j=1,size(widths)
              if (abs(used_list(i,size(point)+j)-widths(j)).gt.0.01) then
                 widths_match=.false.
                 exit
!              else
!                 ewrite(3,*) 'widths match'
              endif
           end do
           if ((points_match).and.(widths_match)) then
              used=.true.
              error=used_list(i,size(used_list,2))
              exit
           endif
        endif
     end do
     
  end subroutine check_already_used
  
  subroutine add_to_used_list(point, widths, error, used_list)
  
     real, dimension(:) :: point, widths
     real, dimension(:,:) :: used_list
     real :: error
     integer :: i, j
     
     ewrite(3,*) 'Adding point, width, error to used_list'
     
     do i=1,size(used_list,1)
        if (used_list(i,1)==0.0) then
           ! insert new value
           used_list(i,1:size(point))=point
           used_list(i,size(point)+1:size(point)+size(widths))=widths
           used_list(i,size(point)+size(widths)+1)=error
!           ewrite(3,*) used_list(i,:)
           return
!        else
!           ewrite(3,*) used_list(i,:)
        endif
     end do

  end subroutine add_to_used_list

end module saturation_distribution_search_hookejeeves