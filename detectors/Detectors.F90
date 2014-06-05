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

module detectors
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use fields
  use field_options
  use state_module
  use spud
  use pickers, only: picker_inquire
  use node_owner_finder, only: node_owner_finder_find
  use fldebug
  use iso_c_binding, only: c_ptr, C_NULL_PTR, C_NULL_CHAR, C_LOC

  implicit none
  
  private

  public :: initialise_detectors_from_options, detectors_write_state

  ! We only need one global detector list (for now)
  type(c_ptr), save :: detector_list = C_NULL_PTR

  ! Reference to the output writer object
  type(c_ptr), save :: detector_writer = C_NULL_PTR

  ! We need the coordinate field to establish the enclosing cell of a particle
  type(vector_field), pointer, save :: xfield => null()  

  interface

     function create_particle_list() result(list) bind (c)
       use iso_c_binding
       implicit none
       type(c_ptr) :: list
     end function create_particle_list

     subroutine particle_list_view(plist) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: plist
     end subroutine particle_list_view

     subroutine particle_list_new_particle(plist, coordinates, dim, name) bind (c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: plist
       real (c_double), dimension(dim), intent(in) :: coordinates
       integer(c_int), intent(in), value :: dim
       character(kind=c_char) :: name(*)
     end subroutine particle_list_new_particle

     function create_particle_writer_stat(plist, filename, binary_format) &
          result(writer) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: plist
       character(kind=c_char) :: filename(*)
       integer(c_int), value :: binary_format
       type(c_ptr) :: writer
     end function create_particle_writer_stat

     subroutine particle_writer_register_field(writer, field, field_ptr, &
          mphase, components) bind(c)
       use iso_c_binding
       implicit none
       character(kind=c_char) :: field(*), mphase(*)
       type(c_ptr), value :: writer, field_ptr
       integer(c_int), value :: components
     end subroutine particle_writer_register_field

     subroutine particle_writer_write_header(writer) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: writer
     end subroutine particle_writer_write_header

     subroutine particle_writer_write_state(writer, time, dt) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: writer
       real (c_double), value :: time, dt
     end subroutine particle_writer_write_state

  end interface

  interface detector_field
     module procedure detector_field_scalar, detector_field_vector
  end interface

contains

  subroutine find_enclosing_cell(coordinates, dim, element, lcoords) bind(c)
    !! Callback routine to establish enclosing element of a particle
    use iso_c_binding
    implicit none
    real (c_double), dimension(dim), intent(in) :: coordinates
    integer(c_int), intent(in), value :: dim
    integer(c_int), intent(out) :: element
    real (c_double), dimension(dim+1), intent(out) :: lcoords

    call picker_inquire(xfield, coordinates, element, local_coord=lcoords, global=.true.)
  end subroutine find_enclosing_cell

  subroutine evaluate_scalar_field(field_ptr, dim, element, local_coords, \
    scalar_value) bind(c)
    !! Callback routine to sample a scalar field at a particle position
    use iso_c_binding
    implicit none
    type(c_ptr), value :: field_ptr
    integer(c_int), intent(in), value :: dim, element
    real(c_double), dimension(dim+1), intent(in) :: local_coords
    real (c_double), intent(inout) :: scalar_value

    type(scalar_field), pointer :: sfield
    call C_F_POINTER(field_ptr, sfield)

    scalar_value = eval_field(element, sfield, local_coords)
  end subroutine evaluate_scalar_field

  subroutine evaluate_vector_field(field_ptr, dim, element, local_coords, \
    vector_value) bind(c)
    !! Callback routine to sample a vector field at a particle position
    use iso_c_binding
    implicit none
    type(c_ptr), value :: field_ptr
    integer(c_int), intent(in), value :: dim, element
    real(c_double), dimension(dim+1), intent(in) :: local_coords
    real(c_double), dimension(dim), intent(inout) :: vector_value

    type(vector_field), pointer :: vfield
    call C_F_POINTER(field_ptr, vfield)

    vector_value = eval_field(element, vfield, local_coords)
  end subroutine evaluate_vector_field

  subroutine initialise_detectors_from_options(state)
    type(state_type), dimension(:), intent(in) :: state
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    integer :: i, ndet, phase, binary_format
    integer, dimension(2) :: shape_option
    character(len=OPTION_PATH_LEN) :: filename, optpath
    character(len=FIELD_NAME_LEN) :: name, field_name, material_phase_name
    real, dimension(:), allocatable :: location

    ! Store a safe reference to the coordinate field, because
    ! we need it to establish a particle's position in the mesh
    xfield => extract_vector_field(state(1), "Coordinate")

    ewrite(1,*) "Initialising detectors from options..."
    detector_list = create_particle_list()

    allocate(location(xfield%dim))

    ! Initialise static detectors
    ndet = option_count("/io/detectors/static_detector")
    do i=1,ndet
       write(optpath, "(a,i0,a)") "/io/detectors/static_detector[",i-1,"]"
       shape_option=option_shape(trim(optpath)//"/location")
       assert(xfield%dim==shape_option(1))
       call get_option(trim(optpath)//"/location", location)
       call get_option(trim(optpath)//"/name", name)

       call particle_list_new_particle(detector_list, &
            location, xfield%dim, trim(name)//C_NULL_CHAR)
    end do

    call particle_list_view(detector_list)

    ! Create a detector writer object
    call get_option("/simulation_name", filename)
    binary_format = 1  ! False
    detector_writer = create_particle_writer_stat(detector_list, &
         trim(filename)//".detectors"//C_NULL_CHAR, binary_format)

    ! Register fields to be included in detector output
    phaseloop: do phase=1,size(state)
       material_phase_name=trim(state(phase)%name)

       ! Scalar fields
       do i = 1, size(state(phase)%scalar_names)
          field_name = trim(state(phase)%scalar_names(i))
          sfield => extract_scalar_field(state(phase), field_name)   
          if (detector_field(sfield)) then
             call particle_writer_register_field(detector_writer, &
                  trim(field_name)//C_NULL_CHAR, C_LOC(sfield), &
                  trim(material_phase_name)//C_NULL_CHAR, 1)
          end if
       end do

       ! Vector fields
       do i = 1, size(state(phase)%vector_names)
          field_name = trim(state(phase)%vector_names(i))
          vfield => extract_vector_field(state(phase), field_name)   
          if (detector_field(vfield)) then
             call particle_writer_register_field(detector_writer, &
                  trim(field_name)//C_NULL_CHAR, C_LOC(vfield), &
                  trim(material_phase_name)//C_NULL_CHAR, vfield%dim)
          end if
       end do
    end do phaseloop

    call particle_writer_write_header(detector_writer)

    deallocate(location)
    ewrite(1,*) "Done initialising detectors"
    
  end subroutine initialise_detectors_from_options

  subroutine detectors_write_state(time, dt)
    real, intent(in) :: time, dt

    call particle_writer_write_state(detector_writer, time, dt)
  end subroutine detectors_write_state

  function detector_field_scalar(sfield)
    !!< Return whether the supplied field should be included in the .detector file
    logical :: detector_field_scalar
    type(scalar_field), target, intent(in) :: sfield
    
    if (sfield%option_path=="".or.aliased(sfield)) then
       detector_field_scalar=.false.
    else
       detector_field_scalar = have_option(&
            trim(complete_field_path(sfield%option_path)) // &
            "/detectors/include_in_detectors")
    end if

  end function detector_field_scalar

  function detector_field_vector(vfield)
    !!< Return whether the supplied field should be included in the .detector file
    logical :: detector_field_vector
    type(vector_field), target, intent(in) :: vfield

    if (vfield%option_path=="".or.aliased(vfield)) then
       detector_field_vector=.false.
    else
       detector_field_vector = have_option(&
            trim(complete_field_path(vfield%option_path)) // &
            "/detectors/include_in_detectors")
    end if

  end function detector_field_vector

end module detectors
