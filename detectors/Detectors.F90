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
  use fields_data_types, only: vector_field
  use state_module, only: state_type, extract_vector_field
  use pickers, only: picker_inquire
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use spud, only: option_count, option_shape, get_option
  use node_owner_finder, only: node_owner_finder_find
  use fldebug
  use iso_c_binding, only: c_ptr, C_NULL_PTR

  implicit none
  
  private

  public :: initialise_detectors_from_options

  ! We only need one global detector list (for now)
  type(c_ptr), save :: detector_list = C_NULL_PTR

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

     subroutine particle_list_new_particle(plist, coordinates, dim) bind (c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: plist
       real (c_double), dimension(dim), intent(in) :: coordinates
       integer(c_int), intent(in), value :: dim
     end subroutine particle_list_new_particle

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

  subroutine initialise_detectors_from_options(state)
    type(state_type), dimension(:), intent(in) :: state
    integer :: i, ndet
    integer, dimension(2) :: shape_option
    character(len=OPTION_PATH_LEN) :: optpath
    character(len=FIELD_NAME_LEN) :: name
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

       call particle_list_new_particle(detector_list, location, xfield%dim)
    end do

    call particle_list_view(detector_list)

    deallocate(location)
    ewrite(1,*) "Done initialising detectors"
    
  end subroutine initialise_detectors_from_options

end module detectors
