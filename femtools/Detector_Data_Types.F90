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
  use global_parameters, only : FIELD_NAME_LEN, OPTION_PATH_LEN, PYTHON_FUNC_LEN
  
  implicit none
  
  private
  
  public :: detector_type, detector_linked_list, &
            detector_list_ptr, stringlist, &
            random_walk, biovar, functional_group, &
            STATIC_DETECTOR, LAGRANGIAN_DETECTOR, &
            GUIDED_SEARCH_TRACKING, RTREE_TRACKING, GEOMETRIC_TRACKING

  integer, parameter :: STATIC_DETECTOR=1, LAGRANGIAN_DETECTOR=2
  integer, parameter :: GUIDED_SEARCH_TRACKING=1, RTREE_TRACKING=2, GEOMETRIC_TRACKING=3

  type stringlist
     !!< Container type for a list of strings.
     character(len=FIELD_NAME_LEN), dimension(:), pointer :: ptr
  end type stringlist

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
     !! Whether the detector is static or Lagrangian.
     integer :: type = STATIC_DETECTOR
     !! Identification number indicating the order in which the detectors are read
     integer :: id_number
     !! ID of the parent list, needed for Zoltan to map the detector back
     integer :: list_id

     !! RK timestepping stages (first index is stage no., second index is dim)
     real, dimension(:,:), allocatable :: k
     !! RK update destination vector (size dim)
     real, dimension(:), allocatable :: update_vector

     ! Definition of the ray for geometric tracking
     real, dimension(:), allocatable :: ray_o, ray_d
     ! Distance to the target coordinate
     real :: target_distance, current_t

     !! Biology variables
     real, dimension(:), allocatable :: biology

     !! Counter for internal Random Walk sub-cycling
     integer :: rw_subsubcycles = 1

     !! Pointers for detector linked lists
     TYPE (detector_type), POINTER :: next=> null()
     TYPE (detector_type), POINTER :: previous=> null() 
  end type detector_type

  !! Type for holding the parameters of random walk schemes, 
  !! so that we may combine multiple random walks.
  type random_walk
    !! Name for logging
    character(len=FIELD_NAME_LEN) :: name

    !! Internal Random Walk schemes
    logical :: naive_random_walk = .false.
    logical :: diffusive_random_walk = .false.

    !! Field names for internal Diffusive Random Walk scheme
    character(len=FIELD_NAME_LEN) :: diffusivity_field, diffusivity_grad
    character(len=FIELD_NAME_LEN) :: diffusivity_2nd_grad

    !! Auto-subcycling (for internal Diffusive Random Walk)
    logical :: auto_subcycle = .false.
    real :: subcycle_scale_factor = 1.0

    !! Python Random Walk function
    logical :: python_random_walk = .false.
    character(len=PYTHON_FUNC_LEN) :: python_code
  end type random_walk

  ! Type holding meta-information about biology variables of LE agents
  type biovar
    character(len=FIELD_NAME_LEN) :: name
    ! Type of variable, ie. diagnostic, uptake, release
    integer :: field_type = 0

    ! Name and option path of the primary diagnostic field
    character(len=FIELD_NAME_LEN) :: field_name
    character(len=OPTION_PATH_LEN) :: field_path
    ! Option path for the depletion field (Uptake vars only)
    character(len=OPTION_PATH_LEN):: depletion_field_path
    ! Name of chemical field associated with uptake/release variables
    character(len=FIELD_NAME_LEN) :: chemfield

    ! Flag indicating whether this variables gets written to the output file
    logical :: write_to_file = .false.
    ! Flag indicating whether to aggregate by stage
    logical :: stage_aggregate = .false.
    ! Variable index of the corresponding pool variable for uptake/release variables
    integer :: pool_index
  end type biovar

  type functional_group
    character(len=FIELD_NAME_LEN) :: name
    ! List of variables that define each agent of this group
    type(biovar), dimension(:), allocatable :: variables
    ! List of all stages within this FG
    type(stringlist) :: stage_names
    ! Option path for the Agents diagnostic field
    character(len=OPTION_PATH_LEN) :: agents_field_path
  end type functional_group

  type detector_linked_list
     !! Doubly linked list implementation
     integer :: length=0
     TYPE (detector_type), pointer :: first => null()
     TYPE (detector_type), pointer :: last => null()

     !! Internal ID used for packing/unpacking detectors
     integer :: id  ! IDs are counted from 1
     !! Name as specified in options (Lagrangian agents only)
     character(len=FIELD_NAME_LEN) :: name

     !! Biology options
     type(functional_group), pointer :: fgroup => null()
     real :: stage_id
     character(len=FIELD_NAME_LEN) :: stage_name
     character(len=OPTION_PATH_LEN) :: stage_options
     character(len=PYTHON_FUNC_LEN) :: biovar_pycode

     character(len=FIELD_NAME_LEN), dimension(:), allocatable :: env_field_name

     !!!!!!!!!!!!!!!!!!!!!!!!!
     !!! Detector movement !!!
     !!!!!!!!!!!!!!!!!!!!!!!!!
     ! Number of user-defined sub-cycles per timestep
     integer :: n_subcycles = 1
     ! Flag indicating whether to apply lagrangian advection
     logical :: velocity_advection = .false.
     ! Flag indicating whether we reflect detectors at the domain boundary
     logical :: reflect_on_boundary = .false.
     ! Flag indicating whether we update physical coordinates when the mesh moves
     logical :: move_with_mesh = .false.

     ! Runk-Kutta Guided Search parameters
     integer :: n_stages = 0
     real, allocatable, dimension(:) :: timestep_weights
     real, allocatable, dimension(:,:) :: stage_matrix

     integer :: tracking_method
     real :: search_tolerance = 1.0e-10

     ! List of random walk schemes to apply
     type(random_walk), dimension(:), allocatable :: random_walks


     !!!!!!!!!!!!!!!!!!!!!!!!!
     !!! Detector I/O      !!!
     !!!!!!!!!!!!!!!!!!!!!!!!!
     !! Optional array for detector names; names are held in read order
     character(len=FIELD_NAME_LEN), dimension(:), allocatable :: detector_names

     !! List of scalar/vector fields to include in detector output
     type(stringlist), dimension(:), allocatable :: sfield_list
     type(stringlist), dimension(:), allocatable :: vfield_list
     integer :: num_sfields = 0   ! Total number of scalar fields across all phases
     integer :: num_vfields = 0   ! Total number of vector fields across all phases

     !! I/O parameters
     logical :: binary_output = .false.
     logical :: write_nan_outside = .false.
     integer :: output_unit = 0          ! Output filehandle, assumed non-opened if 0
     integer :: mpi_write_offset = 0     ! Offset in MPI file
     integer :: total_num_det = 0        ! Global number of detectors in this list
  end type detector_linked_list

  type detector_list_ptr
     type(detector_linked_list), pointer :: ptr
  end type detector_list_ptr

end module detector_data_types
