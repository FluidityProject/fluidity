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
  use elements, only: local_coord_count
  use detector_data_types
  use embed_python, only: set_detectors_from_python
  use integer_hash_table_module
  use fields
  use state_module, only: state_type, extract_scalar_field
  use futils, only: int2str
  use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
  use pickers
  use parallel_tools
  use parallel_fields, only: element_owned
  
  implicit none
  
  private

  public :: insert, allocate, deallocate, copy, move, move_all, remove, &
            delete, delete_all, pack_detector, unpack_detector, &
            detector_value, set_detector_coords_from_python, &
            detector_buffer_size, set_particle_attribute_from_python, &
            set_particle_fields_from_python, set_particle_constant_from_options

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

  subroutine detector_allocate_from_params(new_detector, ndims, local_coord_count, attribute_dims)
    type(detector_type),  pointer, intent(out) :: new_detector
    integer, intent(in) :: ndims, local_coord_count
    integer, optional, intent(in) :: attribute_dims
      
    assert(.not. associated(new_detector))
      
    ! allocate the memory for the new detector
    if (.not. associated(new_detector)) then
       allocate(new_detector)
    end if
    allocate(new_detector%position(ndims))
    allocate(new_detector%local_coords(local_coord_count))
    if (present(attribute_dims)) then
       allocate(new_detector%attributes(attribute_dims))
    end if
      
    assert(associated(new_detector))
      
  end subroutine detector_allocate_from_params
    
  subroutine detector_allocate_from_detector(new_detector, old_detector)
    type(detector_type), pointer, intent(in) :: old_detector
    type(detector_type),  pointer, intent(out) :: new_detector
      
    integer :: ndims, local_coord_count, attribute_dims
      
    ndims = size(old_detector%position)
    local_coord_count = size(old_detector%local_coords)
    attribute_dims = size(old_detector%attributes)
      
    ! allocate the memory for the new detector
    call detector_allocate_from_params(new_detector, ndims, local_coord_count, attribute_dims)
      
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
       if(allocated(detector%attributes)) then
          deallocate(detector%attributes)
       end if
       deallocate(detector)
    end if
    detector => null()

  end subroutine detector_deallocate

  subroutine detector_list_deallocate(detector_list)
    type(detector_linked_list), pointer :: detector_list

    type(rk_gs_parameters), pointer :: parameters

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
    parameters => detector_list%move_parameters
    if (associated(parameters)) then 
       if (allocated(parameters%timestep_weights)) then
          deallocate(parameters%timestep_weights)
       end if
       if (allocated(parameters%stage_matrix)) then
          deallocate(parameters%stage_matrix)
       end if
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

  function detector_buffer_size(ndims, have_update_vector, nstages, attribute_dims)
    ! Returns the number of reals we need to pack a detector
    integer, intent(in) :: ndims
    logical, intent(in) :: have_update_vector
    integer, intent(in), optional :: nstages
    integer :: detector_buffer_size
    integer, optional, intent(in) :: attribute_dims
    if (present(attribute_dims)) then
       if (have_update_vector) then
          assert(present(nstages))
          detector_buffer_size=(nstages+2)*ndims+3+attribute_dims
       else
          detector_buffer_size=ndims+4+attribute_dims
       end if
    else
       if (have_update_vector) then
          assert(present(nstages))
          detector_buffer_size=(nstages+2)*ndims+3
       else
          detector_buffer_size=ndims+4
       end if 
    end if

  end function detector_buffer_size

  subroutine pack_detector(detector,buff,ndims,nstages, attribute_dims)
    ! Packs (serialises) detector into buff
    ! Basic fields are: element, position, id_number and type
    ! If nstages is given, the detector is still moving
    ! and we also pack update_vector and k
    type(detector_type), pointer, intent(in) :: detector
    real, dimension(:), intent(out) :: buff
    integer, intent(in) :: ndims
    integer, intent(in), optional :: nstages
    integer, optional, intent(in) :: attribute_dims
    
    assert(size(detector%position)==ndims)
    if (present(attribute_dims)) then
       assert(size(buff)>=ndims+3+attribute_dims)

       ! Basic fields: ndims+3
       buff(1:ndims) = detector%position
       buff(ndims+1) = detector%element
       buff(ndims+2) = detector%id_number
       buff(ndims+3) = detector%type
       if (attribute_dims.ne.0) then
          buff(ndims+4:ndims+3+attribute_dims) = detector%attributes
       end if
       ! Lagrangian advection fields: (nstages+1)*ndims
       if (present(nstages)) then
          assert(size(buff)==(nstages+2)*ndims+3+attribute_dims)
          assert(allocated(detector%update_vector))
          assert(allocated(detector%k))
          
          buff(ndims+4+attribute_dims:2*ndims+3+attribute_dims) = detector%update_vector
          buff(2*ndims+4+attribute_dims:(nstages+2)*ndims+3+attribute_dims) = reshape(detector%k,(/nstages*ndims/))
       else
          assert(size(buff)==ndims+4+attribute_dims)
          buff(ndims+4+attribute_dims) = detector%list_id
       end if
    else
       assert(size(buff)>=ndims+3)

       ! Basic fields: ndims+3
       buff(1:ndims) = detector%position
       buff(ndims+1) = detector%element
       buff(ndims+2) = detector%id_number
       buff(ndims+3) = detector%type

       ! Lagrangian advection fields: (nstages+1)*ndims
       if (present(nstages)) then
          assert(size(buff)==(nstages+2)*ndims+3)
          assert(allocated(detector%update_vector))
          assert(allocated(detector%k))
          
          buff(ndims+4:2*ndims+4) = detector%update_vector
          buff(2*ndims+4:(nstages+2)*ndims+4) = reshape(detector%k,(/nstages*ndims/))
       else
          assert(size(buff)==ndims+4)
          buff(ndims+4) = detector%list_id
       end if
    end if
    
  end subroutine pack_detector

  subroutine unpack_detector(detector,buff,ndims,global_to_local,coordinates,nstages, attribute_dims)
    ! Unpacks the detector from buff and fills in the blanks
    type(detector_type), pointer :: detector
    real, dimension(:), intent(in) :: buff
    integer, intent(in) :: ndims
    type(integer_hash_table), intent(in), optional :: global_to_local
    type(vector_field), intent(in), optional :: coordinates
    integer, intent(in), optional :: nstages
    integer, optional, intent(in) :: attribute_dims
    if (present(attribute_dims)) then
       assert(size(buff)>=ndims+3+attribute_dims)

       if (.not. allocated(detector%position)) then
          allocate(detector%position(ndims))
       end if
       if (attribute_dims.ne.0) then
          allocate(detector%attributes(attribute_dims))
       end if
       
       ! Basic fields: ndims+3
       detector%position = reshape(buff(1:ndims),(/ndims/))
       detector%element = buff(ndims+1)
       detector%id_number = buff(ndims+2)
       detector%type = buff(ndims+3)
       if (attribute_dims.ne.0) then  
          detector%attributes = reshape(buff(ndims+4:ndims+3+attribute_dims),(/attribute_dims/))
       end if
       
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
       
       ! Lagrangian advection fields: (nstages+1)*ndims
       if (present(nstages)) then
          assert(size(buff)==(nstages+2)*ndims+3+attribute_dims)
          
          ! update_vector, dimension(ndim)
          if (.not. allocated(detector%update_vector)) then
             allocate(detector%update_vector(ndims))
          end if
          detector%update_vector = reshape(buff(ndims+4+attribute_dims:2*ndims+3+attribute_dims),(/ndims/))
          
          ! k, dimension(nstages:ndim)
          if (.not. allocated(detector%k)) then
             allocate(detector%k(nstages,ndims))
          end if
          detector%k = reshape(buff(2*ndims+4+attribute_dims:(nstages+2)*ndims+3+attribute_dims),(/nstages,ndims/))
          
          ! If update_vector still exists, we're not done moving
          detector%search_complete=.false.
       else
          assert(size(buff)==ndims+4+attribute_dims)
          
          detector%list_id = buff(ndims+4+attribute_dims)
          detector%search_complete=.true.
       end if
    else
       assert(size(buff)>=ndims+3)

       if (.not. allocated(detector%position)) then
          allocate(detector%position(ndims))
       end if
       
       ! Basic fields: ndims+3
       detector%position = reshape(buff(1:ndims),(/ndims/))
       detector%element = buff(ndims+1)
       detector%id_number = buff(ndims+2)
       detector%type = buff(ndims+3)
       
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
       
       ! Lagrangian advection fields: (nstages+1)*ndims
       if (present(nstages)) then
          assert(size(buff)==(nstages+2)*ndims+3)
          
          ! update_vector, dimension(ndim)
          if (.not. allocated(detector%update_vector)) then
             allocate(detector%update_vector(ndims))
          end if
          detector%update_vector = reshape(buff(ndims+4:2*ndims+3),(/ndims/))
          
          ! k, dimension(nstages:ndim)
          if (.not. allocated(detector%k)) then
             allocate(detector%k(nstages,ndims))
          end if
          detector%k = reshape(buff(2*ndims+4:(nstages+2)*ndims+3),(/nstages,ndims/))
          
          ! If update_vector still exists, we're not done moving
          detector%search_complete=.false.
       else
          assert(size(buff)==ndims+4)
          
          detector%list_id = buff(ndims+4)
          detector%search_complete=.true.
       end if
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

  subroutine set_particle_constant_from_options(attributes, ndete, constant)
    !!< Given a constant value, set full array of attributes from this value
    real, dimension(:), intent(inout) :: attributes
    integer, intent(in) :: ndete
    real, intent(in) :: constant

    integer :: i

    do i = 1,ndete
       attributes(i) = constant
    end do

  end subroutine set_particle_constant_from_options

  subroutine set_particle_attribute_from_python(attributes, positions, ndete, func, time)
    !!< Given a particle position and time, evaluate the python function
    !!< specified in the string func at that location. 
    real, dimension(:), intent(inout) :: attributes
    !! Func may contain any python at all but the following function must
    !! be defined::
    !!  def val(X,t)
    !! where X is position and t is the time. The result must be a float. 
    character(len=*), intent(in) :: func
    integer, intent(in) :: ndete
    real, intent(in) :: time
    real, dimension(:,:), target, intent(in) :: positions
    real, dimension(:), pointer :: lvx,lvy,lvz
    real, dimension(0), target :: zero
    integer :: stat, dim
    call get_option("/geometry/dimension",dim)

    select case(dim)
    case(1)
      lvx=>positions(1,:)
      lvy=>zero
      lvz=>zero
    case(2)
      lvx=>positions(1,:)
      lvy=>positions(2,:)
      lvz=>zero
    case(3)
      lvx=>positions(1,:)
      lvy=>positions(2,:)
      lvz=>positions(3,:)
    end select
    call set_particles_from_python(func, len(func), dim, &
         ndete, lvx, lvy, lvz, time, attributes, stat)
    if (stat/=0) then
      ewrite(-1, *) "Python error, Python string was:"
      ewrite(-1 , *) trim(func)
      FLExit("Dying")
    end if
  end subroutine set_particle_attribute_from_python

  subroutine set_particle_fields_from_python(state, xfield, dim, positions, ndete, ele, lcoords, attributes, func, time, field_name)
    type(state_type), dimension(:), intent(in) :: state
    real, dimension(:), intent(inout) :: attributes
    character(len=*), intent(in) :: func
    character(len=*), dimension(:), intent(in) :: field_name
    integer, intent(in) :: ndete
    real, intent(in) :: time
    integer, intent(in) :: dim
    real, dimension(:,:), target, intent(in) :: positions
    real, dimension(:,:), intent(in) :: lcoords
    integer, dimension(:), intent(in) :: ele
    real, allocatable, dimension(:,:) :: fields
    type(vector_field), pointer :: xfield
    real :: value
    real, dimension(:), pointer :: lvx,lvy,lvz
    real, dimension(0), target :: zero
    character(len=FIELD_NAME_LEN) :: buffer !set len as number

    type(scalar_field), pointer :: sfield
    character(len=FIELD_NAME_LEN) :: name
    integer :: phase, i, j, nfields
    integer :: p, f, stat, num_fields, k
    logical :: particles_f

    !get positions of particles for function

    select case(dim)
    case(1)
      lvx=>positions(1,:)
      lvy=>zero
      lvz=>zero
    case(2)
      lvx=>positions(1,:)
      lvy=>positions(2,:)
      lvz=>zero
    case(3)
      lvx=>positions(1,:)
      lvy=>positions(2,:)
      lvz=>positions(3,:)
    end select
    num_fields=0
    allocate(fields(size(field_name),ndete))
    do phase=1,size(state)
       nfields = option_count('/material_phase[' &
            //int2str(phase-1)//']/scalar_field')
       do f = 1, nfields
          call get_option('material_phase['//int2str(phase-1)//']/scalar_field['//int2str(f-1)//']/name', name)
          sfield => extract_scalar_field(state(phase),name)
          if (have_option(trim(sfield%option_path)//"/prescribed/particles/include_in_particles").or. &
               have_option(trim(sfield%option_path)//"/diagnostic/particles/include_in_particles").or. &
               have_option(trim(sfield%option_path)//"/prognostic/particles/include_in_particles")) then
             do j=1,size(field_name)
                if (name==field_name(j)) then
                   do i = 1,ndete
                      value = eval_field(ele(i), sfield, lcoords(:,i))
                      fields(j,i) = value
                   end do
                   num_fields = num_fields+1
                end if
             end do
          end if
       end do
    end do
    if (size(field_name).ne.num_fields) then
       ewrite(2,*) "number of fields is not consistent"
       FLExit("Dying")
    end if
    call set_particles_fields_from_python(func, len(func), dim, ndete, &
         lvx, lvy, lvz, time, num_fields, fields, attributes, stat)
    if (stat/=0) then
       ewrite(-1, *) "Python error, Python string was:"
       ewrite(-1 , *) trim(func)
       FLExit("Dying")
    end if
  end subroutine set_particle_fields_from_python

end module detector_tools
