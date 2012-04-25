!    Copyright (C) 2008 Imperial College London and others.
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

module lagrangian_biology_pm
  use fldebug
  use spud
  use state_module
  use fields
  use global_parameters, only: PYTHON_FUNC_LEN
  use detector_data_types
  use detector_tools
  use detector_parallel, only: get_next_detector_id
  use detector_python
  use integer_hash_table_module
  use Profiler

implicit none

  private

  public :: particle_management

  integer, parameter :: BIOVAR_STAGE=1, BIOVAR_SIZE=2

contains

  subroutine particle_management(state, agent_list)
    type(state_type), intent(inout) :: state
    type(detector_linked_list), intent(inout) :: agent_list

    type(detector_type), pointer :: agent, agent_to_move, merge_target, agent_to_merge
    type(scalar_field), pointer :: agent_density_field
    type(vector_field), pointer :: xfield
    character(len=PYTHON_FUNC_LEN) :: pm_pycode
    type(ilist) :: elements_split_list, elements_merge_list
    type(integer_hash_table) :: split_index_mapping, merge_index_mapping, element_minima, element_maxima
    type(detector_linked_list), dimension(:), allocatable :: split_lists, merge_lists
    integer, dimension(:), allocatable :: elements_split, elements_merge, xbiggest_ids, xsmallest_ids
    real, dimension(:), allocatable :: xbiggest_sizes, xsmallest_sizes
    integer :: i, ele, index, splits_wanted, merges_wanted, found_agent, element_minimum, element_maximum
    real :: agent_density, ele_volume, minimum_density, maximum_density
    logical :: do_min_python, do_max_python

    call profiler_tic(trim(agent_list%name)//"::particle_management")

    ewrite(2,*) "Performing Particle Management for list: ", trim(agent_list%name)

    xfield=>extract_vector_field(state, "Coordinate")
    agent_density_field=>extract_scalar_field(state, trim(agent_list%fgroup%name)//"Agents"//trim(agent_list%stage_name))

    ! Initialise the Python functions to determine PM limits
    if (have_option(trim(agent_list%stage_options)//"/particle_management/minimum")) then
       if (have_option(trim(agent_list%stage_options)//"/particle_management/minimum/python")) then
          call get_option(trim(agent_list%stage_options)//"/particle_management/minimum/python", pm_pycode)
          call python_run_detector_string(trim(pm_pycode), trim(agent_list%name), trim("pm_minimum"))
          do_min_python = .true.
       else
          call get_option(trim(agent_list%stage_options)//"/particle_management/minimum/constant", minimum_density)
          do_min_python = .false.
       end if
    else
       do_min_python = .false.
       minimum_density = 0.0
    end if

    if (have_option(trim(agent_list%stage_options)//"/particle_management/maximum")) then
       if (have_option(trim(agent_list%stage_options)//"/particle_management/maximum/python")) then
          call get_option(trim(agent_list%stage_options)//"/particle_management/maximum/python", pm_pycode)
          call python_run_detector_string(trim(pm_pycode), trim(agent_list%name), trim("pm_maximum"))
          do_max_python = .true.
       else
          do_max_python = .false.
          call get_option(trim(agent_list%stage_options)//"/particle_management/maximum/constant", maximum_density)
       end if
    else
       do_max_python = .false.
       maximum_density = huge(1.0)
    end if

    ! First establish in which elements we need to split/merge
    ! We skip elements with no agents...
    ewrite(2,*) "Particle Management: Establishing elements to split/merge in"
    call allocate(element_minima)
    call allocate(element_maxima)
    do ele=1, ele_count(agent_density_field)
       agent_density = node_val(agent_density_field, ele)
       ele_volume = element_volume(xfield, ele)

       ! Here we either use the global minima/maxima from the options,
       ! or run Python code to determine the minimum/maximum for each element
       if (do_min_python) then
          call python_get_element_limit(ele, xfield, trim(agent_list%name), trim("pm_minimum"), minimum_density)
       end if

       if (do_max_python) then
          call python_get_element_limit(ele, xfield, trim(agent_list%name), trim("pm_maximum"), maximum_density)
       end if

       ! Store element numbers and minima/maxima for split and merge
       if (agent_density < minimum_density .and. agent_density > 0) then
          call insert(elements_split_list, ele)
          call insert(element_minima, ele, nint(minimum_density * ele_volume))
       end if

       if (agent_density > maximum_density .and. agent_density > 0) then
          call insert(elements_merge_list, ele)
          call insert(element_maxima, ele, nint(maximum_density * ele_volume))
       end if
    end do

    ! Convert linked lists of elements to split/merge to vectors
    allocate(elements_split(elements_split_list%length))
    elements_split = list2vector(elements_split_list)
    allocate(elements_merge(elements_merge_list%length))
    elements_merge = list2vector(elements_merge_list)

    ! Create a hashtable from element -> index into split_lists
    call allocate(split_index_mapping)
    do i=1, size(elements_split)
       call insert(split_index_mapping, elements_split(i), i)
    end do

    ! Create a hashtable from element -> index into merge_lists
    call allocate(merge_index_mapping)
    do i=1, size(elements_merge)
       call insert(merge_index_mapping, elements_merge(i), i)
    end do

    ! Put agents into temporary agent lists for each element to split/merge
    ewrite(2,*) "Particle Management: Moving agents to split/merge into temporary lists"
    allocate(split_lists(elements_split_list%length))
    allocate(merge_lists(elements_merge_list%length))
    agent=>agent_list%first
    do while(associated(agent))
       ! 
       if (has_key(split_index_mapping, agent%element)) then
          index = fetch(split_index_mapping,agent%element)
          agent_to_move=>agent
          agent=>agent%next
          call move(agent_to_move, agent_list, split_lists(index))
       ! 
       elseif (has_key(merge_index_mapping, agent%element)) then
          index = fetch(merge_index_mapping,agent%element)
          agent_to_move=>agent
          agent=>agent%next
          call move(agent_to_move, agent_list, merge_lists(index))

       else
          agent=>agent%next
       end if
    end do 

    ! Perform splits over the temporary agent arrays
    ! Each array now represents all agents in one element
    ! the according element number is stored in elements_split(i)
    ewrite(2,*) "Particle Management: Splitting agents per element"
    do i=1, size(split_lists)
       element_minimum = fetch(element_minima, elements_split(i))

       !!!!!!!!!!!!!!!!!!!!!
       ! Original VEW split algorithm:
       do while(split_lists(i)%length < element_minimum)
          ! 1) Find x, the number of splits we want
          splits_wanted = element_minimum - split_lists(i)%length
          splits_wanted = min(splits_wanted, split_lists(i)%length)

          ! 2) Get the x biggest agents
          allocate(xbiggest_ids(splits_wanted))
          allocate(xbiggest_sizes(splits_wanted))
          xbiggest_sizes = 0.0

          agent=>split_lists(i)%first
          do while(associated(agent))
             if (minval(xbiggest_sizes) < agent%biology(BIOVAR_SIZE)) then
                index = minval(minloc(xbiggest_sizes))
                xbiggest_ids(index) = agent%id_number 
                xbiggest_sizes(index) = agent%biology(BIOVAR_SIZE)
             end if

             agent=>agent%next
          end do

          ! 3) Split each of the x agents...
          agent=>split_lists(i)%first
          do while(associated(agent))
             if (find(xbiggest_ids, agent%id_number) > 0) then
                call pm_split(agent, split_lists(i))
             end if
             agent=>agent%next
          end do

          ! ...and recurse
          deallocate(xbiggest_ids)
          deallocate(xbiggest_sizes)
       end do
       !!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!
       ! Simple algorithm: always split first agent in the list
       !do while(agent_density < element_minimum)
       !   call pm_split(split_lists(i)%first, split_lists(i))
       !   agent_density = split_lists(i)%length / ele_volume
       !end do
       !!!!!!!!!!!!!!!!!!!!!

    end do

    ! Copy all agents back to the original list
    do i=1, size(split_lists)
       call move_all(split_lists(i), agent_list)
    end do

    ! Deallocate everythin split related
    deallocate(split_lists)
    deallocate(elements_split)
    call flush_list(elements_split_list)
    call deallocate(split_index_mapping)
    call deallocate(element_minima)

    ! Perform merges over the temporary agent arrays
    ewrite(2,*) "Particle Management: Merging agents per element"
    do i=1, size(merge_lists)
       element_maximum = fetch(element_maxima, elements_merge(i))

       !!!!!!!!!!!!!!!!!!!!!
       ! Original VEW merge algorithm:
       do while(merge_lists(i)%length > element_maximum)
          ! 1) Find x, the number of merges we want
          merges_wanted = merge_lists(i)%length - element_maximum
          ! We need to make sure we don't have more merge requests than agent pairs
          merges_wanted = floor(min(real(merges_wanted), real(merge_lists(i)%length) / 2.0))

          ! 2) Get the 2*x smallest agents (x pairs)
          allocate(xsmallest_ids(2*merges_wanted))
          allocate(xsmallest_sizes(2*merges_wanted))
          xsmallest_sizes = huge(1.0)

          agent=>merge_lists(i)%first
          do while(associated(agent))
             if (maxval(xsmallest_sizes) > agent%biology(BIOVAR_SIZE)) then
                index = minval(maxloc(xsmallest_sizes))
                xsmallest_ids(index) = agent%id_number 
                xsmallest_sizes(index) = agent%biology(BIOVAR_SIZE)
             end if

             agent=>agent%next
          end do

          ! 3) Merge pairwise...
          ! Note: We don't sort, but select merge pairs as we find them
          merge_target=>null()
          agent=>merge_lists(i)%first
          do while(associated(agent))
             found_agent = find(xsmallest_ids, agent%id_number)
             ! First hit is the target...
             if (found_agent > 0 .and. .not.associated(merge_target)) then
                merge_target=>agent
                agent=>agent%next
             ! ... second hit; we merge and reset the target
             elseif (found_agent > 0 .and. associated(merge_target)) then
                agent_to_merge=>agent
                agent=>agent%next
                call pm_merge(xfield, merge_target, agent_to_merge, merge_lists(i))
                merge_target=>null()
             else
                agent=>agent%next
             end if
          end do

          ! ...and recurse
          deallocate(xsmallest_ids)
          deallocate(xsmallest_sizes)
       end do
       !!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!
       ! Simple algorithm: always merge the first two agents in the list
       !do while(agent_density > agent_list%pm_max)
       !   call pm_merge(merge_lists(i)%first, merge_lists(i)%first%next, merge_lists(i))
       !   agent_density = merge_lists(i)%length / ele_volume
       !end do
       !!!!!!!!!!!!!!!!!!!!!

    end do

    ! Copy all agents back to the original list
    do i=1, size(merge_lists)
       call move_all(merge_lists(i), agent_list)
    end do

    ! Deallocate eveything merge related
    deallocate(merge_lists)
    deallocate(elements_merge)
    call flush_list(elements_merge_list)
    call deallocate(merge_index_mapping)
    call deallocate(element_maxima)

    call profiler_toc(trim(agent_list%name)//"::particle_management")

  contains

    function find(array, val) result(loc)
      !!< Find the first instance of val in array.
      integer, intent(in), dimension(:) :: array
      integer, intent(in) :: val
      integer :: i, loc

      loc = -1
      do i=1,size(array)
        if (array(i) == val) then
          loc = i
          return
        end if
      end do
    end function find

  end subroutine particle_management

  subroutine pm_split(agent, agent_list)
    type(detector_type), pointer, intent(inout) :: agent
    type(detector_linked_list), intent(inout) :: agent_list

    type(detector_type), pointer :: new_agent
    integer :: i

    ! Allocate and insert agent
    new_agent=>null()
    call allocate(new_agent, agent)
    call insert(new_agent,agent_list)

    ! Populate new agent
    call get_next_detector_id(new_agent%id_number)
    new_agent%name=trim(int2str(agent_list%total_num_det))
    new_agent%position=agent%position
    new_agent%element=agent%element
    new_agent%local_coords=agent%local_coords
    new_agent%type=agent%type
    new_agent%list_id=agent%list_id

    ! Populate new agent's biology
    allocate(new_agent%biology(size(agent%biology)))
    do i=1, size(agent%biology)
       new_agent%biology(i)=agent%biology(i)
    end do

    ! Distribute biomass
    agent%biology(BIOVAR_SIZE)=agent%biology(BIOVAR_SIZE) / 2.0
    new_agent%biology(BIOVAR_SIZE)=new_agent%biology(BIOVAR_SIZE) / 2.0

    ! Update the list
    agent_list%total_num_det=agent_list%total_num_det + 1

  end subroutine pm_split

  subroutine pm_merge(xfield, agent1, agent2, agent_list)
    type(vector_field), pointer, intent(in) :: xfield
    type(detector_type), pointer, intent(inout) :: agent1, agent2
    type(detector_linked_list), intent(inout) :: agent_list

    integer :: i

    ! Both agents are in the same element, 
    ! so we weight-average each physical coordinate
    ! and re-set the local coordinates
    do i=1, size(agent1%position)
       agent1%position(i)=wtavg(agent1%position(i),agent1%biology(BIOVAR_SIZE),agent2%position(i),agent2%biology(BIOVAR_SIZE))
    end do
    agent1%local_coords=local_coords(xfield,agent1%element,agent1%position)

    ! Weight-average the biology variables
    if (size(agent1%biology)>2) then
       do i=3, size(agent1%biology)
          agent1%biology(i) = wtavg(agent1%biology(i),agent1%biology(BIOVAR_SIZE),agent2%biology(i),agent2%biology(BIOVAR_SIZE))
       end do
    end if

    ! Add the sizes
    agent1%biology(BIOVAR_SIZE) = agent1%biology(BIOVAR_SIZE) + agent2%biology(BIOVAR_SIZE)

    call delete(agent2, agent_list)
    agent_list%total_num_det=agent_list%total_num_det - 1

  contains

    function wtavg(v1, w1, v2, w2) result(val)
      real, intent(in) :: v1, w1, v2, w2
      real :: val

      val = ((v1*w1)+(v2*w2))/(w1+w2)
    end function wtavg

  end subroutine pm_merge

end module lagrangian_biology_pm
