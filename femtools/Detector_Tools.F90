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

module detector_tools
  use spud
  use fldebug
  use detector_data_types
  use fields
  use integer_hash_table_module
  use element_path_list, only: elepath_list_destroy => list_destroy
  implicit none
  
  private

  public :: insert, allocate, deallocate, copy, move, move_all, remove, &
            delete, delete_all, pack_detector, unpack_detector, &
            detector_value, set_detector_coords_from_python, &
            detector_buffer_size

  interface insert
     module procedure insert_into_detector_list
  end interface

  ! Removes detector from a list without deallocating it
  interface remove
     module procedure remove_detector_from_list
  end interface

  ! Removes detector from a list and deallocates it
  interface delete
     module procedure delete_detector
  end interface

  ! Deletes all detectors from a given list
  interface delete_all
     module procedure delete_all_detectors
  end interface

  ! Move detector from one list to another
  interface move
     module procedure move_detector
  end interface

  ! Move all detectors in a list from one that list to another
  interface move_all
     module procedure move_all_detectors
  end interface

  interface copy
     module procedure detector_copy
  end interface

  interface allocate
     module procedure detector_allocate_from_params, detector_allocate_from_detector
  end interface

  interface deallocate
     module procedure detector_deallocate, detector_list_deallocate
  end interface

  ! Evaluate field at the location of the detector.
  interface detector_value
     module procedure detector_value_scalar, detector_value_vector
  end interface

contains 

  subroutine detector_allocate_from_params(new_detector, ndims, local_coord_count)
    type(detector_type), pointer, intent(out) :: new_detector
    integer, intent(in) :: ndims, local_coord_count
    
    assert(.not. associated(new_detector))
      
    ! allocate the memory for the new detector
    if (.not. associated(new_detector)) then
       allocate(new_detector)
    end if
    allocate(new_detector%position(ndims))
    allocate(new_detector%local_coords(local_coord_count))
      
    assert(associated(new_detector))

    new_detector%path_elements => null()
      
  end subroutine detector_allocate_from_params
    
  subroutine detector_allocate_from_detector(new_detector, old_detector)
    type(detector_type), pointer, intent(in) :: old_detector
    type(detector_type), pointer, intent(out) :: new_detector
      
    integer :: ndims, local_coord_count
      
    ndims = size(old_detector%position)
    local_coord_count = size(old_detector%local_coords)
      
    ! allocate the memory for the new detector
    call detector_allocate_from_params(new_detector, ndims, local_coord_count)
      
  end subroutine detector_allocate_from_detector
    
  subroutine detector_deallocate(detector)
    type(detector_type), pointer :: detector
      
    if(associated(detector)) then
       if(allocated(detector%local_coords)) then
          deallocate(detector%local_coords)
       end if
       if(allocated(detector%position)) then
          deallocate(detector%position)
       end if
       if(allocated(detector%k)) then
          deallocate(detector%k)
       end if
       if(allocated(detector%update_vector)) then
          deallocate(detector%update_vector)
       end if
       if(allocated(detector%biology)) then
          deallocate(detector%biology)
       end if
       if(allocated(detector%food_requests)) then
          deallocate(detector%food_requests)
       end if
       if(allocated(detector%food_ingests)) then
          deallocate(detector%food_ingests)
       end if
       if(allocated(detector%ray_o)) then
          deallocate(detector%ray_o)
       end if
       if(allocated(detector%ray_d)) then
          deallocate(detector%ray_d)
       end if
       if (associated(detector%path_elements)) then
          call elepath_list_destroy(detector%path_elements)
       end if

       deallocate(detector)
    end if
    detector => null()

  end subroutine detector_deallocate

  subroutine detector_list_deallocate(detector_list)
    type(detector_linked_list), pointer :: detector_list

    ! Delete detectors
    if (detector_list%length > 0) then
       call delete_all(detector_list)
    end if

    ! Deallocate list information
    if (allocated(detector_list%detector_names)) then
       deallocate(detector_list%detector_names)
    end if
    if (allocated(detector_list%sfield_list)) then
       deallocate(detector_list%sfield_list)
    end if
    if (allocated(detector_list%vfield_list)) then
       deallocate(detector_list%vfield_list)
    end if

    ! Deallocate move_parameters
    if (allocated(detector_list%timestep_weights)) then
       deallocate(detector_list%timestep_weights)
    end if
    if (allocated(detector_list%stage_matrix)) then
       deallocate(detector_list%stage_matrix)
    end if

  end subroutine detector_list_deallocate
    
  subroutine detector_copy(new_detector, old_detector)
    ! Copies all the information from the old detector to
    ! the new detector
    type(detector_type), pointer, intent(in) :: old_detector
    type(detector_type),  pointer :: new_detector
      
    new_detector%position = old_detector%position
    new_detector%element = old_detector%element
    new_detector%id_number = old_detector%id_number
    new_detector%type = old_detector%type
    new_detector%name = old_detector%name
    new_detector%local_coords=old_detector%local_coords
      
  end subroutine detector_copy

  subroutine insert_into_detector_list(detector, current_list)
    ! Inserts detector at the end of a list
    type(detector_linked_list), intent(inout) :: current_list
    type(detector_type), pointer :: detector

    if (current_list%length == 0) then
       current_list%first => detector 
       current_list%last => detector 
       current_list%first%previous => null()
       current_list%last%next => null()
       current_list%length = 1
    else
       detector%previous => current_list%last
       current_list%last%next => detector
       current_list%last => detector
       current_list%last%next => null()
       current_list%length = current_list%length+1
    end if

  end subroutine insert_into_detector_list

  subroutine remove_detector_from_list(detector, detector_list)
    !! Removes the detector from the list, 
    !! but does not deallocated it
    type(detector_linked_list), intent(inout) :: detector_list
    type(detector_type), pointer :: detector

    if (detector_list%length==1) then
       detector_list%first => null()
       detector_list%last => null()
    else
       if (associated(detector%previous)) then
          detector%previous%next => detector%next
       else
          detector_list%first => detector%next
       end if

       if (associated(detector%next)) then
          detector%next%previous => detector%previous
       else
          detector_list%last => detector%previous
       end if
    end if

    detector_list%length = detector_list%length-1

  end subroutine remove_detector_from_list

  subroutine delete_detector(detector, detector_list)
    ! Removes and deallocates the given detector 
    ! and outputs the next detector in the list as detector
    type(detector_type), pointer :: detector
    type(detector_linked_list), intent(inout), optional :: detector_list
    
    type(detector_type), pointer :: temp_detector
    
    if (present(detector_list)) then
       temp_detector => detector
       detector => detector%next
       call remove(temp_detector, detector_list)
       call deallocate(temp_detector)
    else
       call deallocate(detector)
    end if
      
  end subroutine delete_detector

  subroutine move_detector(detector, from_list, to_list)
    ! Move detector from one list to the other
    type(detector_linked_list), intent(inout) :: from_list
    type(detector_type), pointer :: detector
    type(detector_linked_list), intent(inout) :: to_list

    call remove(detector, from_list)
    call insert(detector, to_list)  

  end subroutine move_detector

  subroutine move_all_detectors(from_list,to_list)
    ! Move all detectors from one list to the other
    type(detector_linked_list), intent(inout) :: from_list
    type(detector_linked_list), intent(inout) :: to_list
    type(detector_type), pointer :: detector

    do while (associated(from_list%first))
       detector => from_list%first
       call move(detector, from_list, to_list)   
    end do

  end subroutine move_all_detectors

  subroutine delete_all_detectors(detector_list)
    ! Remove and deallocate all detectors in a list   
    type(detector_linked_list), intent(inout) :: detector_list
    type(detector_type), pointer :: detector

    detector => detector_list%first
    do while (associated(detector))
       call delete(detector,detector_list) 
    end do

  end subroutine delete_all_detectors

  function detector_buffer_size(ndims, detector, list, tracking) result(buffer_size)
    ! Returns the number of reals we need to pack a detector
    integer, intent(in) :: ndims
    type(detector_type), pointer, intent(inout), optional :: detector
    type(detector_linked_list), intent(in), optional :: list
    logical, intent(in), optional :: tracking
    integer :: buffer_size

    ! Basics: Position(ndims) + element + id + type 
    buffer_size = ndims+3

    if (present(list)) then

       if (present_and_true(tracking)) then
          ! update_vector
          buffer_size = buffer_size + ndims

          if (list%n_stages > 0) then
             ! RK advection: k(nstages)
             buffer_size = buffer_size + list%n_stages*ndims
          end if

          if (list%tracking_method == GEOMETRIC_TRACKING) then
             ! Ray tracking: ray_o(ndims) + ray_d(ndims) + target_distance + current_t
             buffer_size = buffer_size + 2*ndims + 2
          end if
          if (list%tracking_method == PURE_GS) then
             ! Ray tracking: ray_o(ndims)
             buffer_size = buffer_size + ndims
          end if
       end if

       ! LEBiology
       if (associated(list%fgroup)) then
          buffer_size = buffer_size + size(list%fgroup%variables)

          if (size(list%fgroup%food_sets) > 0) then
             buffer_size = buffer_size + 3 * size(list%fgroup%food_sets(1)%varieties)
          end if
       end if

    else
       ! Zoltan sends a list ID with the detector
       buffer_size = buffer_size + 1
    end if

  end function detector_buffer_size

  subroutine pack_detector(detector, buffer, ind, dim, list, tracking)
    ! Packs (serialises) detector into buff
    ! Basic fields are: element, position, id_number and type
    ! If nstages is given, the detector is still moving
    ! and we also pack update_vector and k
    type(detector_type), pointer, intent(inout) :: detector
    real, dimension(:), intent(inout) :: buffer
    integer, intent(inout) :: ind, dim
    type(detector_linked_list), intent(in), optional :: list
    logical, intent(in), optional :: tracking

    integer :: nvars

    assert(size(detector%position)==dim)
    assert(size(buffer)>=dim+3)

    ! Basic fields: dim+3
    buffer(ind:ind+dim-1) = detector%position
    ind = ind + dim
    buffer(ind)   = real( detector%element )
    buffer(ind+1) = real( detector%id_number )
    buffer(ind+2) = real( detector%type )
    ind = ind + 3

    ! Lagrangian advection fields: (nstages+1)*dim
    if (present(list)) then

       if (present_and_true(tracking)) then
          assert(allocated(detector%update_vector))
          buffer(ind:ind+dim-1) = detector%update_vector
          ind = ind + dim

          if (list%n_stages > 0) then
             assert(allocated(detector%k))
             buffer(ind:ind+list%n_stages*dim-1) = reshape(detector%k,(/list%n_stages*dim/))
             ind = ind + list%n_stages*dim
          end if

          if (list%tracking_method == GEOMETRIC_TRACKING) then
             ! Ray tracking: ray_o + ray_d
             buffer(ind) = detector%target_distance
             ind = ind+1
             buffer(ind) = detector%current_t
             ind = ind+1
             buffer(ind:ind+dim-1) = detector%ray_o
             ind = ind + dim
             buffer(ind:ind+dim-1) = detector%ray_d
             ind = ind + dim
          end if

          if (list%tracking_method == PURE_GS) then
             ! GS tracking: ray_o 
             buffer(ind:ind+dim-1) = detector%ray_o
             ind = ind + dim
          end if
       end if

       ! LEBiology
       if (associated(list%fgroup)) then
          nvars = size(list%fgroup%variables)
          buffer(ind:ind+nvars-1) = detector%biology
          ind = ind + nvars

          if (size(list%fgroup%food_sets) > 0) then
             nvars = size(list%fgroup%food_sets(1)%varieties)
             buffer(ind:ind+nvars-1) = detector%food_requests
             ind = ind + nvars
             buffer(ind:ind+nvars-1) = detector%food_ingests
             ind = ind + nvars
             buffer(ind:ind+nvars-1) = detector%food_thresholds
             ind = ind + nvars
          end if
       end if

    else
       buffer(ind) = detector%list_id
       ind = ind+1
    end if

  end subroutine pack_detector

  subroutine unpack_detector(detector, buffer, ind, dim, global_to_local, coordinates, list, tracking)
    ! Unpacks the detector from buffer and fills in the blanks
    type(detector_type), pointer, intent(inout) :: detector
    real, dimension(:), intent(inout) :: buffer
    integer, intent(inout) :: ind, dim
    type(integer_hash_table), intent(in), optional :: global_to_local
    type(vector_field), intent(in), optional :: coordinates
    type(detector_linked_list), intent(in), optional :: list
    logical, intent(in), optional :: tracking

    integer :: nvars

    assert(size(buffer)>=dim+3)

    ! Basic fields: dim+3
    if (.not. allocated(detector%position)) allocate(detector%position(dim))
    detector%position = buffer(ind:ind+dim-1)
    ind = ind + dim
    detector%element   = nint( buffer(ind) )
    detector%id_number = nint( buffer(ind+1) )
    detector%type      = nint( buffer(ind+2) )
    ind = ind + 3

    ! Reconstruct element number if global-to-local mapping is given
    if (present(global_to_local)) then
       assert(has_key(global_to_local, detector%element))
       detector%element=fetch(global_to_local,detector%element)

       ! Update local coordinates if coordinate field is given
       if (present(coordinates)) then
          if (.not. allocated(detector%local_coords)) then
             allocate(detector%local_coords(local_coord_count(ele_shape(coordinates,1))))
          end if
          detector%local_coords=local_coords(coordinates,detector%element,detector%position)
       end if
    end if

    ! Lagrangian advection fields: (nstages+1)*dim
    if (present(list)) then

       if (present_and_true(tracking)) then
          ! update_vector, dimension(ndim)
          if (.not.allocated(detector%update_vector)) allocate(detector%update_vector(dim))  
          detector%update_vector = reshape(buffer(ind:ind+dim-1),(/dim/))
          ind = ind + dim

          if (list%n_stages > 0) then
             ! k, dimension(nstages:ndim)
             if (.not. allocated(detector%k)) allocate(detector%k(list%n_stages,dim))
             detector%k = reshape(buffer(ind:ind+list%n_stages*dim-1),(/list%n_stages,dim/))
             ind = ind + list%n_stages*dim
          end if

          if (list%tracking_method == GEOMETRIC_TRACKING) then
             ! Ray tracking: ray_o + ray_d
             detector%target_distance = buffer(ind)
             ind = ind+1
             detector%current_t = buffer(ind)
             ind = ind+1

             if (.not.allocated(detector%ray_o)) allocate(detector%ray_o(dim))
             detector%ray_o = buffer(ind:ind+dim-1)
             ind = ind + dim

             if (.not. allocated(detector%ray_d)) allocate(detector%ray_d(dim))
             detector%ray_d = buffer(ind:ind+dim-1)
             ind = ind + dim
          end if

          if (list%tracking_method == PURE_GS) then
             ! GS tracking: ray_o
             if (.not.allocated(detector%ray_o)) allocate(detector%ray_o(dim))
             detector%ray_o = buffer(ind:ind+dim-1)
             ind = ind + dim

             if (present(coordinates)) then  
                detector%local_coords=local_coords(coordinates,detector%element,detector%ray_o)
             end if
          end if
       end if

       ! LEBiology
       if (associated(list%fgroup)) then
          nvars = size(list%fgroup%variables)
          if (.not.allocated(detector%biology)) allocate(detector%biology(nvars))
          detector%biology = buffer(ind:ind+nvars-1)
          ind = ind + nvars

          if (size(list%fgroup%food_sets) > 0) then
             nvars = size(list%fgroup%food_sets(1)%varieties)

             if (.not.allocated(detector%food_requests)) allocate(detector%food_requests(nvars))
             detector%food_requests = buffer(ind:ind+nvars-1)
             ind = ind + nvars

             if (.not.allocated(detector%food_ingests)) allocate(detector%food_ingests(nvars))
             detector%food_ingests = buffer(ind:ind+nvars-1)
             ind = ind + nvars

             if (.not.allocated(detector%food_thresholds)) allocate(detector%food_thresholds(nvars))
             detector%food_thresholds = buffer(ind:ind+nvars-1)
             ind = ind + nvars
          end if
       end if

    else
       detector%list_id = buffer(ind)
       ind = ind+1
    end if

  end subroutine unpack_detector

  function detector_value_scalar(sfield, detector) result(value)
    !!< Evaluate field at the location of the detector.
    real :: value
    type(scalar_field), intent(in) :: sfield
    type(detector_type), intent(in) :: detector
    
    assert(detector%element>0)
    value = eval_field(detector%element, sfield, detector%local_coords)

  end function detector_value_scalar

  function detector_value_vector(vfield, detector) result(value)
    !!< Evaluate field at the location of the detector.
    type(vector_field), intent(in) :: vfield
    type(detector_type), intent(in) :: detector
    real, dimension(vfield%dim) :: value
    
    assert(detector%element>0)
    value = eval_field(detector%element, vfield, detector%local_coords)

  end function detector_value_vector

  subroutine set_detector_coords_from_python(values, ndete, func, time)
    !!< Given a list of positions and a time, evaluate the python function
    !!< specified in the string func at those points. 
    real, dimension(:,:), target, intent(inout) :: values
    !! Func may contain any python at all but the following function must
    !! be defiled:
    !!  def val(t)
    !! where t is the time. The result must be a float. 
    character(len=*), intent(in) :: func
    real :: time
    
    real, dimension(:), pointer :: lvx,lvy,lvz
    real, dimension(0), target :: zero
    integer :: stat, dim, ndete

    call get_option("/geometry/dimension",dim)

    lvx=>values(1,:)
    lvy=>zero
    lvz=>zero
    if(dim>1) then
       lvy=>values(2,:)
       if(dim>2) then
          lvz => values(3,:)
       end if
    end if

    call set_detectors_from_python(func, len(func), dim, &
         ndete, time, dim,                               &
         lvx, lvy, lvz, stat)

    if (stat/=0) then
      ewrite(-1, *) "Python error, Python string was:"
      ewrite(-1 , *) trim(func)
      FLExit("Dying")
    end if

  end subroutine set_detector_coords_from_python

end module detector_tools
