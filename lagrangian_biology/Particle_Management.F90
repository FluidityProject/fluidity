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

module proxy_agent_module
  use detector_data_types
  type proxy_agent
     type(detector_type), pointer :: agent
     integer :: id
     real :: size
  end type proxy_agent
end module proxy_agent_module

module proxy_agent_list
  use proxy_agent_module, LIST_DATA => proxy_agent

  include "../flibs/linkedlist.f90"
end module proxy_agent_list

module lagrangian_biology_pm
  use fldebug
  use spud
  use state_module
  use fields
  use global_parameters, only: PYTHON_FUNC_LEN
  use detector_data_types
  use detector_tools
  use detector_parallel, only: get_next_detector_id
  use integer_hash_table_module
  use pickers
  use Profiler
  use proxy_agent_module
  use proxy_agent_list, only: &
       proxy_list => linked_list, &
       proxy_list_create => list_create, &
       proxy_list_insert => list_insert, &
       proxy_list_next => list_next, &
       proxy_list_count => list_count, &
       proxy_list_remove => list_delete_element, &
       proxy_list_destroy => list_destroy

implicit none

  private

  public :: particle_management, pm_strip_insignificant

  integer, parameter :: BIOVAR_STAGE=1, BIOVAR_SIZE=2

  type proxy_list_ptr
     type(proxy_list), pointer :: ptr => null()
  end type proxy_list_ptr

contains

  subroutine particle_management(state, agent_list)
    type(state_type), intent(inout) :: state
    type(detector_linked_list), intent(inout) :: agent_list

    type(scalar_field), pointer :: agent_field, agent_min_field, agent_max_field
    type(vector_field), pointer :: xfield
    integer :: local_maximum

    call profiler_tic(trim(agent_list%name)//"::particle_management")

    ewrite(2,*) "Performing Particle Management for list: ", trim(agent_list%name)

    xfield=>extract_vector_field(state, "Coordinate")
    agent_field=>extract_scalar_field(state, trim(agent_list%fgroup%name)//"Agents"//trim(agent_list%stage_name))

    ! Split agents per element until the prescribed minimum is achieved
    if (have_option(trim(agent_list%stage_options)//"/particle_management/scalar_field::AgentsMin")) then
       agent_min_field=>extract_scalar_field(state, trim(agent_list%fgroup%name)//"AgentsMin"//trim(agent_list%stage_name))

       call pm_split_by_element(xfield, agent_field, agent_min_field, agent_list)
    end if

    ! Merge agents per element until the prescribed maximum is achieved
    if (have_option(trim(agent_list%stage_options)//"/particle_management/scalar_field::AgentsMax")) then
       agent_max_field=>extract_scalar_field(state, trim(agent_list%fgroup%name)//"AgentsMax"//trim(agent_list%stage_name))

       call pm_merge_by_element(xfield, agent_field, agent_max_field, agent_list)
    end if

    ! After the element-wise split/merges are done,
    ! we check for global agent maxima and enforce them
    if (have_option(trim(agent_list%stage_options)//"/particle_management/local_maximum")) then
       FLExit("local_maximum merging is temporarily disabled")
       ewrite(2,*) "Particle Management: Enforcing local maximum"
       call get_option(trim(agent_list%stage_options)//"/particle_management/local_maximum", local_maximum)

       ! For now we are simply enforcing a local maximum here
       ! Ideally we would enforce one global maximum across all procs
       !call pm_merge_list(xfield, agent_list, local_maximum)
    end if

    call profiler_toc(trim(agent_list%name)//"::particle_management")

  end subroutine particle_management

  subroutine pm_split_by_element(xfield, agent_field, agent_min_field, agent_list)
    type(vector_field), pointer, intent(inout) :: xfield
    type(scalar_field), pointer, intent(inout) :: agent_field, agent_min_field
    type(detector_linked_list), intent(inout) :: agent_list

    type(detector_type), pointer :: agent
    type(proxy_list_ptr), dimension(:), allocatable :: proxy_lists
    type(proxy_agent) :: proxy
    type(integer_hash_table) :: split_index_mapping, element_minima
    type(ilist) :: elements_split_list
    integer, dimension(:), allocatable :: elements_split
    integer :: i, ele, index, element_minimum
    real :: agent_minimum, agent_total

    ewrite(2,*) "Particle Management: Splitting agents per element"

    ! First establish in which elements we need to split/merge
    ! We skip elements with no agents...
    call allocate(element_minima)
    do ele=1, ele_count(agent_field)
       agent_total = integral_element(agent_field, xfield, ele)
       agent_minimum = integral_element(agent_min_field, xfield, ele)

       if (agent_total < agent_minimum .and. agent_total > 0.0) then
          call insert(elements_split_list, ele)
          call insert(element_minima, ele, nint(agent_minimum))
       end if
    end do

    ! Convert linked lists of elements to split/merge to vectors
    allocate(elements_split(elements_split_list%length))
    elements_split = list2vector(elements_split_list)

    ! Create a hashtable from element -> index into split_lists
    call allocate(split_index_mapping)
    do i=1, size(elements_split)
       call insert(split_index_mapping, elements_split(i), i)
    end do

    ! Create proxy agents in temporary proxy lists for each element to split/merge
    allocate(proxy_lists(size(elements_split)))
    agent=>agent_list%first
    do while(associated(agent))
       if (has_key(split_index_mapping, agent%element)) then
          index = fetch(split_index_mapping,agent%element)

          proxy%agent => agent
          proxy%id = agent%id_number
          proxy%size = agent%biology(BIOVAR_SIZE)

          if (associated(proxy_lists(index)%ptr)) then
             call proxy_list_insert( proxy_lists(index)%ptr, proxy )
          else
             call proxy_list_create( proxy_lists(index)%ptr, proxy )
          end if
       end if
       agent=>agent%next
    end do 

    ! Perform splits over the temporary proxy agents.
    ! Each proxy array represents all agents in one element.
    do i=1, size(proxy_lists)
       element_minimum = fetch(element_minima, elements_split(i))
       call pm_split_list(proxy_lists(i)%ptr, element_minimum, agent_list)
       call proxy_list_destroy(proxy_lists(i)%ptr)
    end do

    call deallocate(split_index_mapping)
    deallocate(proxy_lists)
    deallocate(elements_split)
    call flush_list(elements_split_list)
    call deallocate(element_minima)
  end subroutine pm_split_by_element

  subroutine pm_split_list(split_list, minimum, agent_list)
    type(proxy_list), pointer, intent(inout) :: split_list
    type(detector_linked_list), intent(inout) :: agent_list
    integer, intent(in) :: minimum

    type(detector_type), pointer :: agent
    type(proxy_list), pointer :: proxy
    type(proxy_agent) :: new_proxy
    integer, dimension(:), allocatable :: xbiggest_ids
    real, dimension(:), allocatable :: xbiggest_sizes
    integer :: index, splits_wanted

    ! Original VEW split-x-biggest algorithm
    do while(proxy_list_count(split_list) < minimum)
       ! 1) Find x, the number of splits we want
       splits_wanted = minimum - proxy_list_count(split_list) !agent_list%length
       splits_wanted = min(splits_wanted, proxy_list_count(split_list)) !agent_list%length)

       ! 2) Get the x biggest agents
       allocate(xbiggest_ids(splits_wanted))
       allocate(xbiggest_sizes(splits_wanted))
       xbiggest_sizes = 0.0

       proxy => split_list
       do while(associated(proxy))
          if (minval(xbiggest_sizes) <= proxy%data%size) then
             index = minval(minloc(xbiggest_sizes))
             xbiggest_ids(index) = proxy%data%id 
             xbiggest_sizes(index) = proxy%data%size  !agent%biology(BIOVAR_SIZE)
          end if

          proxy => proxy_list_next(proxy)
       end do

       ! 3) Split each of the x agents...
       proxy => split_list
       do while(associated(proxy))
          if (find(xbiggest_ids, proxy%data%id) > 0) then
             call pm_split(proxy%data%agent, agent_list)
             
             ! Add new proxy to list
             new_proxy%agent => agent_list%last
             new_proxy%id = agent_list%last%id_number
             new_proxy%size = agent_list%last%biology(BIOVAR_SIZE)
             call proxy_list_insert( proxy, new_proxy )             
          end if
          proxy => proxy_list_next(proxy)
       end do

       ! ...and recurse
       deallocate(xbiggest_ids)
       deallocate(xbiggest_sizes)
    end do
  end subroutine pm_split_list

  subroutine pm_merge_by_element(xfield, agent_field, agent_max_field, agent_list)
    type(vector_field), pointer, intent(inout) :: xfield
    type(scalar_field), pointer, intent(inout) :: agent_field, agent_max_field
    type(detector_linked_list), intent(inout) :: agent_list

    type(detector_type), pointer :: agent, agent_to_move, merge_target, agent_to_merge
    type(proxy_list_ptr), dimension(:), allocatable :: proxy_lists
    type(proxy_agent) :: proxy
    type(integer_hash_table) :: element_maxima, merge_index_mapping
    type(ilist) :: elements_split_list, elements_merge_list
    integer, dimension(:), allocatable :: elements_merge
    integer :: i, ele, index, element_maximum
    real :: agent_maximum, agent_total

    ewrite(2,*) "Particle Management: Merging agents per element"

    ! First establish in which elements we need to split/merge
    ! We skip elements with no agents...
    call allocate(element_maxima)
    do ele=1, ele_count(agent_field)
       agent_total = integral_element(agent_field, xfield, ele)
       agent_maximum = integral_element(agent_max_field, xfield, ele)

       if (agent_total > agent_maximum .and. agent_total > 0.0) then
          call insert(elements_merge_list, ele)
          call insert(element_maxima, ele, nint(agent_maximum))
       end if
    end do

    allocate(elements_merge(elements_merge_list%length))
    elements_merge = list2vector(elements_merge_list)

    ! Create a hashtable from element -> index into merge_lists
    call allocate(merge_index_mapping)
    do i=1, size(elements_merge)
       call insert(merge_index_mapping, elements_merge(i), i)
    end do

    ! Create agent proxiess in temporary proxy lists for each element to split/merge
    allocate(proxy_lists(size(elements_merge)))
    agent=>agent_list%first
    do while(associated(agent))
       if (has_key(merge_index_mapping, agent%element)) then
          index = fetch(merge_index_mapping,agent%element)

          proxy%agent => agent
          proxy%id = agent%id_number
          proxy%size = agent%biology(BIOVAR_SIZE)

          if (associated(proxy_lists(index)%ptr)) then
             call proxy_list_insert( proxy_lists(index)%ptr, proxy )
          else
             call proxy_list_create( proxy_lists(index)%ptr, proxy )
          end if
       end if
       agent=>agent%next
    end do 

    ! Perform merges over the temporary agent arrays
    do i=1, size(proxy_lists)
       element_maximum = fetch(element_maxima, elements_merge(i))
       call pm_merge_list(xfield, proxy_lists(i)%ptr, element_maximum, agent_list)
       call proxy_list_destroy(proxy_lists(i)%ptr)
    end do

    call deallocate(merge_index_mapping)
    deallocate(proxy_lists)
    deallocate(elements_merge)
    call flush_list(elements_merge_list)
    call deallocate(element_maxima)

  end subroutine pm_merge_by_element

  subroutine pm_merge_list(xfield, merge_list, maximum, agent_list)
    type(vector_field), pointer, intent(inout) :: xfield
    type(proxy_list), pointer, intent(inout) :: merge_list
    type(detector_linked_list), intent(inout) :: agent_list
    integer, intent(in) :: maximum

    type(detector_type), pointer :: agent
    type(proxy_list), pointer :: proxy, merge_target
    integer, dimension(:), allocatable :: xsmallest_ids
    real, dimension(:), allocatable :: xsmallest_sizes
    integer :: index, merges_wanted, found_agent

    ! Original VEW merge-x-smallest algorithm
    do while(proxy_list_count(merge_list) > maximum .and. proxy_list_count(merge_list) > 1)
       ! 1) Find x, the number of merges we want
       merges_wanted = proxy_list_count(merge_list) - maximum
       ! We need to make sure we don't have more merge requests than agent pairs
       merges_wanted = floor(min(real(merges_wanted), real(proxy_list_count(merge_list)) / 2.0))

       ! 2) Get the 2*x smallest agents (x pairs)
       allocate(xsmallest_ids(2*merges_wanted))
       allocate(xsmallest_sizes(2*merges_wanted))
       xsmallest_sizes = huge(1.0)

       proxy => merge_list
       do while(associated(proxy))
          if (maxval(xsmallest_sizes) >= proxy%data%size) then
             index = minval(maxloc(xsmallest_sizes))
             xsmallest_ids(index) = proxy%data%id
             xsmallest_sizes(index) = proxy%data%size
          end if
          proxy => proxy_list_next(proxy)
       end do

       ! 3) Merge pairwise...
       ! Note: We don't sort, but select merge pairs as we find them
       merge_target=>null()
       proxy => merge_list
       do while(associated(proxy))
          found_agent = find(xsmallest_ids, proxy%data%id)
          ! First hit is the target...
          if (found_agent > 0 .and. .not.associated(merge_target)) then
             merge_target => proxy
             proxy => proxy_list_next(proxy)
          ! ... second hit; we merge and reset the target
          elseif (found_agent > 0 .and. associated(merge_target)) then
             call pm_merge(xfield, proxy%data%agent, merge_target%data%agent, agent_list)

             ! Merge proxies
             proxy%data%size = proxy%data%size + merge_target%data%size
             call proxy_list_remove(merge_list, merge_target)
             merge_target=>null()
             proxy => proxy_list_next(proxy)
          else
             proxy => proxy_list_next(proxy)
          end if
       end do

       ! ...and recurse
       deallocate(xsmallest_ids)
       deallocate(xsmallest_sizes)
    end do

  end subroutine pm_merge_list

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

    if (allocated(agent%env_samples)) then
       allocate(new_agent%env_samples(size(agent%env_samples)))
    end if

    ! Distribute biomass
    agent%biology(BIOVAR_SIZE)=agent%biology(BIOVAR_SIZE) / 2.0
    new_agent%biology(BIOVAR_SIZE)=new_agent%biology(BIOVAR_SIZE) / 2.0

    ! Update the list
    agent_list%total_num_det=agent_list%total_num_det + 1

  end subroutine pm_split

  subroutine pm_merge(xfield, agent1, agent2, agent_list)
    type(vector_field), pointer, intent(inout) :: xfield
    type(detector_type), pointer, intent(inout) :: agent1, agent2
    type(detector_linked_list), intent(inout) :: agent_list

    integer :: i

    ! Weight-average each physical coordinate,
    ! re-establish element with the picker if the agents do not share an element,
    ! and re-set the local coordinates
    do i=1, size(agent1%position)
       agent1%position(i)=wtavg(agent1%position(i),agent1%biology(BIOVAR_SIZE),agent2%position(i),agent2%biology(BIOVAR_SIZE))
    end do
    if (agent1%element /= agent2%element) then
       call picker_inquire(xfield, agent1%position, agent1%element, local_coord=agent1%local_coords, global=.false.)
    end if
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

  subroutine pm_strip_insignificant(agent_list)
    type(detector_linked_list), intent(inout) :: agent_list

    type(detector_type), pointer :: agent, agent_to_delete
    integer :: kill_count

    kill_count = 0
    agent=>agent_list%first
    do while(associated(agent))
       if (agent%biology(BIOVAR_SIZE) < 1.0e-10) then
          agent_to_delete => agent
          agent => agent%next
          call delete(agent_to_delete, agent_list)
          kill_count = kill_count + 1
       else
          agent => agent%next
       end if
    end do 

    if (kill_count > 0) then
       ewrite(2,*) "Lagrangian biology: Stripped", kill_count, " insignificant agents"
    end if
  end subroutine pm_strip_insignificant

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

end module lagrangian_biology_pm
