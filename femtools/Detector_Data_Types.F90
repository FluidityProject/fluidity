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

module detector_data_types

  use fldebug
  use global_parameters, only : FIELD_NAME_LEN
  
  implicit none
  
  private
  
  public :: detector_type, rk_gs_parameters, detector_linked_list, &
            detector_list_ptr, stringlist, attr_names_type, attr_write_type, field_phase_type, &
            allocate, deallocate

  type stringlist
    !!< Container type for a list of strings.
    character(len=FIELD_NAME_LEN), dimension(:), pointer :: ptr
  end type stringlist

  type attr_names_type
    !< A bundling of names of scalar, vector, and tensor attributes.
    character(len=FIELD_NAME_LEN), dimension(:), allocatable :: s, v, t
    !! The array length of each individual attribute.
    !! A value of 0 indicates a scalar attribute, otherwise it is
    !! array-valued with the specified dimension.
    integer, dimension(:), allocatable :: sn, vn, tn
  end type attr_names_type

  type field_phase_type
    !< A bundling of material_phase number for each
    !! scalar, vector, and tensor field that is to be
    !! included on particles.
    integer, dimension(:), allocatable :: s, v, t
  end type field_phase_type

  type attr_write_type
    !< A bundling of whether to include attributes in
    !! the output file, grouped by scalar, vector, and tensor to match
    !! other attribute-related datatypes.
    logical, dimension(:), allocatable :: s, v, t
  end type attr_write_type

  interface allocate
    module procedure allocate_attr_names, allocate_field_phases
  end interface allocate

  interface deallocate
    module procedure deallocate_attr_names
  end interface deallocate

  !! Type for caching detector position and search information.
  type detector_type
     !! Physical location of the detector.
     real, dimension(:), allocatable :: position
     !! Name of the detector in input and output.
     character(len=FIELD_NAME_LEN) :: name
     !! Element number in which the detector lies.
     integer :: element
     !! Local coordinates of the detector in that element.
     real, dimension(:), allocatable :: local_coords
     !! Identification number indicating the order in which the detectors are read
     integer :: id_number
     !! Identification number indicating parent processor when detector was created
     integer :: proc_id
     !! ID of the parent list, needed for Zoltan to map the detector back
     integer :: list_id
     !! RK timestepping stages (first index is stage no., second index is dim)
     real, dimension(:,:), allocatable :: k
     !! RK update destination vector (size dim)
     real, dimension(:), allocatable :: update_vector
     !! Attributes carried by particles.
     real, dimension(:), allocatable :: attributes
     !! Attributes carried by particles at the previous timestep.
     real, dimension(:), allocatable :: old_attributes
     !! Interpolated field values at the particle position at the previous timestep.
     real, dimension(:), allocatable :: old_fields
     !! Have we completed the search?
     logical :: search_complete
     !! Pointers for detector linked lists
     TYPE (detector_type), POINTER :: next=> null()
     TYPE (detector_type), POINTER :: previous=> null()
     TYPE (detector_type), POINTER :: temp_next => null()
     TYPE (detector_type), POINTER :: temp_previous => null()
  end type detector_type

  ! Parameters for lagrangian detector movement
  type rk_gs_parameters
    !! Runge-Kutta Guided Search parameters
    integer :: n_stages, n_subcycles
    !! Timestep_weights give the weights of temporal positions
    real, allocatable, dimension(:) :: timestep_weights
    !! Timestep_nodes give the locations of temporal positions
    real, allocatable, dimension(:) :: timestep_nodes
    !! Stage_matrix gives the weights to RK function values
    real, allocatable, dimension(:,:) :: stage_matrix
    real :: search_tolerance
  end type rk_gs_parameters

  type detector_linked_list
     !! Doubly linked list implementation
     integer :: length=0
     TYPE (detector_type), pointer :: first => null()
     TYPE (detector_type), pointer :: last => null()

     !! Internal ID used for packing/unpacking detectors
     integer :: id  ! IDs are counted from 1
     integer :: proc_part_count = 0!Counter for the number of particles spawned on the current processor

     !! Parameters for lagrangian movement (n_stages, stage_matrix, etc)
     type(rk_gs_parameters), pointer :: move_parameters => null()
     logical :: move_with_mesh = .false.

     !! Optional array for detector names; names are held in read order
     character(len = FIELD_NAME_LEN), dimension(:), allocatable :: detector_names

     !! List of scalar/vector fields to include in detector output
     type(stringlist), dimension(:), allocatable :: sfield_list
     type(stringlist), dimension(:), allocatable :: vfield_list
     integer :: num_sfields = 0   ! Total number of scalar fields across all phases
     integer :: num_vfields = 0   ! Total number of vector fields across all phases

     !! Total number of arrays stored for attributes and fields on a particle subgroup
     integer, dimension(3) :: total_attributes
     !! Whether attributes should be written or not
     type(attr_write_type) :: attr_write
     !! Names of attributes and fields stored in a particle subgroup
     type(attr_names_type) :: attr_names, old_attr_names, field_names, old_field_names
     !! The phase of each field that is used in particle attribute calculations
     type(field_phase_type) :: field_phases, old_field_phases

     !! I/O parameters
     logical :: write_nan_outside = .false.
     integer(kind=8) :: h5_id = -1 ! H5hut output identifier
     integer :: total_num_det = 0        ! Global number of detectors in this list
  end type detector_linked_list

  type detector_list_ptr
     type(detector_linked_list), pointer :: ptr
  end type detector_list_ptr

contains

  !> Allocate the attribute name type, given an array of
  !! the number of scalar, vector, and tensor components.
  subroutine allocate_attr_names(attr_names, counts)
    type(attr_names_type), intent(out) :: attr_names
    integer, dimension(3), intent(in) :: counts

    allocate(attr_names%s(counts(1)))
    allocate(attr_names%v(counts(2)))
    allocate(attr_names%t(counts(3)))

    allocate(attr_names%sn(counts(1)))
    allocate(attr_names%vn(counts(2)))
    allocate(attr_names%tn(counts(3)))

    attr_names%sn(:) = 0
    attr_names%vn(:) = 0
    attr_names%tn(:) = 0
  end subroutine allocate_attr_names

  !> Allocate the field phase type, given an array of
  !! the number of scalar, vector, and tensor components.
  subroutine allocate_field_phases(field_phases, counts)
    type(field_phase_type), intent(out) :: field_phases
    integer, dimension(3), intent(in) :: counts

    allocate(field_phases%s(counts(1)))
    allocate(field_phases%v(counts(2)))
    allocate(field_phases%t(counts(3)))
  end subroutine allocate_field_phases

  subroutine deallocate_attr_names(attr_names)
    type(attr_names_type), intent(inout) :: attr_names

    deallocate(attr_names%s)
    deallocate(attr_names%v)
    deallocate(attr_names%t)

    deallocate(attr_names%sn)
    deallocate(attr_names%vn)
    deallocate(attr_names%tn)
  end subroutine deallocate_attr_names

end module detector_data_types
