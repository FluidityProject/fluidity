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
!    C.Pain@Imperial.ac.uk
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

  subroutine search_saturations_hookejeeves(state,i)

    ! declare interface variables
    type(state_type), dimension(:), intent(inout) :: state
    type(state_type), dimension(:), pointer :: working_state
    integer, intent(in) :: i

    ! declare working variables for this subrtn
    character(len=50) :: positions
    character(len=OPTION_PATH_LEN) :: option_buffer
    real, allocatable :: base_point(:)
    real :: step_length, nought, search_min, search_max
    type(scalar_field) :: saturation
    type(scalar_field), pointer :: Cv
    type(vector_field), pointer :: coordinates

    integer :: sections, stat, j
    logical :: is_vwell, do_PS

    ewrite(3,*) 'In search_saturations'

    option_buffer = '/material_phase['//int2str(i-1)//']/electrical_properties/Saturation_Distribution_Search/'

    if (have_option(trim(option_buffer)//'search_criteria_vertical_well')) then
       is_vwell=.true.
       option_buffer = '/material_phase['//int2str(i-1)//']/electrical_properties/Saturation_Distribution_Search/search_criteria_vertical_well/'
    elseif (have_option(trim(option_buffer)//'search_criteria_horizontal_well')) then
       is_vwell=.false.
       option_buffer = '/material_phase['//int2str(i-1)//']/electrical_properties/Saturation_Distribution_Search/search_criteria_horizontal_well/'
    else
       FLExit('Cannot find well geometry')
    endif
    
    !! Optimisation based on Hooke and Jeeves algorithm (Basic Optimisation Methods, BD Bunday, 1984)
    !set step_length
    call get_option(trim(option_buffer)//'initial_step_length', step_length)
    !set dimension of coordinates (number of variables to be optimised)
    call get_option(trim(option_buffer)//'sections', sections)
    allocate(base_point(sections))
    
    if (is_vwell) then
       call get_option(trim(option_buffer)//'y_min',search_min)
       call get_option(trim(option_buffer)//'y_max',search_max)
    else
       call get_option(trim(option_buffer)//'z_min',search_min)
       call get_option(trim(option_buffer)//'z_max',search_max)
    endif

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
    call set_base_point(state, do_PS, base_point, step_length, base_point, nought, search_min, search_max)
    
    deallocate(base_point)
   
  end subroutine search_saturations_hookejeeves
    
  recursive subroutine set_base_point(state, do_PS, base_point, step_length, previous_base_point, &
                                      current_error, search_min, search_max)
  
     type(state_type), dimension(:), intent(inout) :: state
     real, dimension(:) :: base_point
     real, dimension(size(base_point)) :: copy_base_point
     real :: step_length, minimum_step_length, new_error, search_min, search_max
     logical :: do_PS, is_improved, finished
     ! optional argument previous_base_point allows us to make pattern step
     real, dimension(:) :: previous_base_point
     ! optional argument current_error records error at this base_point
     real :: current_error, error
!     integer, save :: dump_no
     
     ewrite(3,*) 'setting base point', base_point
     
     ! Set minimum_step_length for now
     minimum_step_length = 30.0
  
     if (do_PS) then
!        call write_state(dump_no, state)
        call make_pattern_step(state, base_point, previous_base_point, step_length, minimum_step_length, &
                               current_error, search_min, search_max)
     else
        ! This is the first time through so need to calculate a first error
        call run_model(state, base_point, step_length, error)
!        dump_no=1
!        call write_state(dump_no, state)
        ! Copy current base point for pattern search in case base_point is updated by exploration
        copy_base_point=base_point
        call make_exploration(state, base_point, step_length, is_improved, error, search_min, search_max)
        if(.not.is_improved) then
           call reduce_step_length(state, base_point, base_point, step_length, minimum_step_length, &
                                   error, finished, search_min, search_max)
           if(.not.finished) then
              call set_base_point(state, do_PS, base_point, step_length, base_point, error, search_min, search_max)
           else
              ewrite(3,*) 'Finishing search'
              return ! FINISHED! PRINT RESULTS ETC ETC
           endif
        else
           call set_base_point(state, .TRUE., base_point, step_length, copy_base_point, error, search_min, search_max)
        endif
     endif
  
  end subroutine set_base_point
  
  subroutine make_pattern_step(state, base_point, previous_base_point, step_length, minimum_step_length, &
                               current_error, search_min, search_max)
  
     type(state_type), dimension(:), intent(inout) :: state
     real, dimension(:) :: base_point, previous_base_point
     real, dimension(size(base_point)) :: pattern_point
     real :: step_length, current_error, minimum_step_length, search_min, search_max
     integer :: i
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
     call make_exploration(state, pattern_point, step_length, is_improved, current_error, search_min, search_max)
     finished=.false.
     if (.not.is_improved) then
        call reduce_step_length(state, base_point, pattern_point, step_length, minimum_step_length, &
                                current_error, finished, search_min, search_max)
     else
        call set_base_point(state, .TRUE., pattern_point, step_length, base_point, current_error, search_min, search_max)
     endif
  
  end subroutine make_pattern_step
  
  subroutine make_exploration(state, base_point, step_length, is_improved, error, search_min, search_max)
     
     type(state_type), dimension(:), intent(inout) :: state
     real, dimension(:) :: base_point
     real, allocatable, dimension(:) :: base_point_in
     real :: step_length, error, new_error, search_min, search_max
     logical :: is_improved
     integer :: i
     
     ewrite(3,*) 'making exploration, step ', step_length
     
     ! return new base point in base_point if it's been improved
     ! ie overwrite
  
     ! cycle through co-ordinates:
     is_improved=.false.
     
     allocate(base_point_in(size(base_point)))
     do i=1,size(base_point)
        ! increase by step length
        base_point(i)=base_point(i)+step_length
        base_point(i)=max(base_point(i), search_min)
        base_point(i)=min(base_point(i), search_max)
        call run_model(state, base_point, step_length, new_error)
        if (new_error>=error) then
           !decrease by step length
           base_point(i)=base_point_in(i)-step_length
           base_point(i)=max(base_point(i), search_min)
           base_point(i)=min(base_point(i), search_max)
           call run_model(state, base_point, step_length, new_error)
           if (new_error>=error) then
              ! keep original
              base_point(i)=base_point_in(i)
           else
              ! set coordinate to new value
              ewrite(3,*) 'Error at', base_point(i), new_error, 'less than error at', base_point(i)+step_length, error
              error=new_error
              is_improved=.true.
           endif
        else
           ! set coordinate to new value
           ewrite(3,*) 'Error at', base_point(i), new_error, 'less than error at', base_point(i)-step_length, error
           error=new_error
           is_improved=.true.
        endif
     end do
     deallocate(base_point_in)
  
  end subroutine make_exploration
  
  subroutine run_model(state, point, step_length, error)
  
     type(state_type), dimension(:), intent(inout) :: state
     type(scalar_field) :: saturation, temp_saturation
     type(scalar_field), pointer :: Cv, electrical_potential
     type(vector_field), pointer :: coordinates
     real, dimension(:) :: point
     real, allocatable, dimension(:,:), save :: target_potential
     real, allocatable, dimension(:) :: model_potential
     real :: step_length, error
     real :: x_min, x_max, y_min, y_max, z_min, z_max
     real :: bh_x, bh_y, bh_z, err
     logical :: is_vwell
     integer :: stat, i, j
     integer, save :: have_target, samples, dump_no
     integer, allocatable, dimension(:), save :: node_list
     character(len=OPTION_PATH_LEN) :: option_buffer, target_filename
     
     ewrite(3,*) 'running model, base_point:',base_point
     
     ! Get all the limits etc relevant - do this fresh each time from flml
     option_buffer = '/material_phase['//int2str(i-1)//']/electrical_properties/Saturation_Distribution_Search/'

     if (have_option(trim(option_buffer)//'search_criteria_vertical_well')) then
        is_vwell=.true.
        option_buffer = '/material_phase['//int2str(i-1)//']/electrical_properties/Saturation_Distribution_Search/search_criteria_vertical_well/'
     elseif (have_option(trim(option_buffer)//'search_criteria_horizontal_well')) then
        is_vwell=.false.
        option_buffer = '/material_phase['//int2str(i-1)//']/electrical_properties/Saturation_Distribution_Search/search_criteria_horizontal_well/'
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
     
     ! Set saturations to point
     saturation = extract_scalar_field(state, "PhaseVolumeFraction", stat)
     call allocate(temp_saturation, saturation%mesh, name="OriginalSaturation")
     ewrite(3,*) 'copied to temporary saturation field'
     call set(temp_saturation,saturation)
     coordinates=>extract_vector_field(state, "Coordinate")
     
     call set_saturations(saturation, point, is_vwell, &
                          x_min, x_max, y_min, y_max, z_min, z_max, &
                          coordinates)

     ! get target potential curve if it hasn't been got before
     if (have_target.ne.8) then
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
              call get_option(trim(option_buffer)//'borehole_x',bh_x)
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
     Cv => extract_scalar_field(state, 'Electrokinetic[0]', stat=stat)
     if (stat/=0) then
        FLExit('Did not find an Electrokinetic coupling coefficient scalar field.')
     end if

     call calculate_diagnostic_variable(state, 1, Cv) 
     ! Get streaming potential
     call calculate_electrical_potential(state(1), 1)
     ewrite(3,*) 'out of calculate_electrical_potential'
     
     ! get potential curve from new model
     electrical_potential=>extract_scalar_field(state, "ElectricalPotential", stat)
     allocate(model_potential(samples))
     
     ewrite(3,*) 'getting potentials from latest model'
     do i=1,samples
        model_potential(i)=electrical_potential%val(node_list(i))
     end do

     ! Compare and calculate error
!     ewrite(3,*) 'error=',error
!     ewrite(3,*) 'target_potential, model_potential',target_potential(1:samples,2), model_potential
     call curve_error(target_potential(1:samples,2), model_potential, error)
     
     ewrite(3,*) 'point, error: ', point, error
     deallocate(model_potential)
     
     dump_no=max(1,dump_no)
     call write_state(dump_no, state)
     
     ! Set saturation back to the initial condition ready for next time
     call set(saturation,temp_saturation)
     call deallocate(temp_saturation)

  end subroutine run_model
  
  subroutine set_saturations(saturation, point, is_vwell, &
                             x_min, x_max, y_min, y_max, z_min, z_max, &
                             coordinates)
    
    type(scalar_field) :: saturation
    type(vector_field), intent(in) :: coordinates
    real, dimension (:) :: point
    real, intent(in) :: x_min, x_max, y_min, y_max, z_min, z_max
    real :: along_step
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
          if ((coords(1,1).gt.x_min-err).and.(coords(1,1).lt.x_max+err).and.&
                (coords(coord_along,1).gt.zone_along_min-err).and.(coords(coord_along,1).lt.zone_along_min+along_step+err).and.&
                (coords(coord_perp,1).lt.point(i)+err)) then
            do node=2,coordinates%mesh%shape%loc
              if ((coords(1,node).gt.x_min-err).and.(coords(1,node).lt.x_max+err).and.&
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
  
  recursive subroutine reduce_step_length(state, base_point, pattern_point, step_length, minimum_step_length, &
                                          current_error, finished, search_min, search_max)
  
     type(state_type), dimension(:), intent(inout) :: state
     real, dimension(:) :: base_point, pattern_point
     real :: step_length, minimum_step_length, current_error, search_min, search_max
     logical :: finished, is_improved, from_pattern_step
     integer :: i
     
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
           call set_base_point(state, .FALSE., base_point, step_length, base_point, &
                               current_error, search_min, search_max)
        else
           ewrite(3,*) 'Reached minimum step length, one last try'
           ! Abusing file names here!:
           pattern_point=base_point
           call make_exploration(state, base_point, step_length, is_improved, &
                                 current_error, search_min, search_max)
           if (is_improved) then
              call set_base_point(state, .TRUE., base_point, step_length, pattern_point, &
                                  current_error, search_min, search_max)
           else
              ewrite(3,*) 'Failed to improve with minimum step length, finished'
              return
           endif
        endif
     else
        ewrite(3,*) 'New step length: ', step_length
        call make_exploration(state, pattern_point, step_length, is_improved, &
                              current_error, search_min, search_max)
        if(is_improved) then
           call set_base_point(state, .TRUE., pattern_point, step_length, base_point, &
                               current_error, search_min, search_max)
        else
           call reduce_step_length(state, base_point, pattern_point, step_length, minimum_step_length, &
                                   current_error, finished, search_min, search_max)
        endif
     endif
  
  end subroutine reduce_step_length
  
  subroutine curve_error(curve1, curve2, error)
  
     real, dimension(:), intent(in) :: curve1, curve2
     real, allocatable, dimension(:) :: curve1q, curve2q, ncurve1, ncurve2
     real :: error, maxi1, maxi2
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
     
     maxi1=maxval(curve1q)
     maxi2=maxval(curve2q)
     if ((maxi1>0.0).and.(maxi2>0.0)) then
        ncurve1=curve1q/maxi1
        ncurve2=curve2q/maxi2
        do k=1,min(size(ncurve1),size(ncurve2))
           error = error + abs(ncurve1(k)-ncurve2(k))
        end do
     else
        error=1000.0
     endif
     ewrite(3,*) 'error is ', error
     ewrite(3,*) 'target',curve1
     ewrite(3,*) 'curve2',curve2
     ewrite(3,*) 'targetq',curve1q
     ewrite(3,*) 'curve2q',curve2q
     ewrite(3,*) 'maxi1, maxi2', maxi1, maxi2
     ewrite(3,*) 'targetn', ncurve1
     ewrite(3,*) 'ncurve2', ncurve2
     
     deallocate(curve1q)
     deallocate(curve2q)
     deallocate(ncurve1)
     deallocate(ncurve2)
  
  end subroutine curve_error
    
end module saturation_distribution_search_hookejeeves