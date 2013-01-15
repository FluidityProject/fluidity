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

module path_element_module
  type path_element
     integer :: ele
     real :: dist
  end type path_element
end module path_element_module

module element_path_list
  use path_element_module, LIST_DATA => path_element

  include "../flibs/linkedlist.f90"
end module element_path_list

module detector_data_types
  use fldebug
  use global_parameters, only : FIELD_NAME_LEN, OPTION_PATH_LEN, PYTHON_FUNC_LEN
  use linked_lists
  use integer_hash_table_module
  use element_path_list, only: ele_path_list => linked_list
  
  implicit none
  
  private
  
  public :: detector_type, detector_linked_list, detector_list_ptr, stringlist, detector_buffer, &
            random_walk, le_variable, functional_group, food_set, food_variety, &
            STATIC_DETECTOR, LAGRANGIAN_DETECTOR, &
            GUIDED_SEARCH_TRACKING, GEOMETRIC_TRACKING, PURE_GS

  integer, parameter :: STATIC_DETECTOR=1, LAGRANGIAN_DETECTOR=2
  integer, parameter :: GUIDED_SEARCH_TRACKING=1, RTREE_TRACKING=2, GEOMETRIC_TRACKING=3, PURE_GS=4

  type stringlist
     !!< Container type for a list of strings.
     character(len=FIELD_NAME_LEN), dimension(:), pointer :: ptr
  end type stringlist

  type detector_list_ptr
     ! Container type for detector_linked_list
     type(detector_linked_list), pointer :: ptr
  end type detector_list_ptr

  type detector_buffer
    real, dimension(:), pointer :: ptr
  end type detector_buffer

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
     ! We record face numbers for use with reflection
     integer :: current_face

     ! Element list and distances for path integration
     type(ele_path_list), pointer :: path_elements => null()

     !! Biology variables
     real, dimension(:), allocatable :: biology
     real, dimension(:), allocatable :: food_requests, food_ingests, food_thresholds

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

    !! Auto-subcycling (for internal Diffusive Random Walk)
    logical :: auto_subcycle = .false.

    !! Python Random Walk function
    logical :: python_random_walk = .false.
    character(len=PYTHON_FUNC_LEN) :: python_code
  end type random_walk

  ! Type holding meta-information about biology variables of LE agents
  type le_variable
    character(len=FIELD_NAME_LEN) :: name
    ! Type of variable, ie. diagnostic, uptake, release
    integer :: field_type = 0

    ! Name and option path of the primary diagnostic field
    character(len=FIELD_NAME_LEN) :: field_name
    character(len=OPTION_PATH_LEN) :: field_path
    ! Option path for the depletion field (Uptake vars only)
    character(len=OPTION_PATH_LEN) :: depletion_field_path
    ! Name of chemical field associated with uptake/release variables
    character(len=FIELD_NAME_LEN) :: chemfield
    integer :: i_chemfield

    ! Flag indicating whether this variables gets written to the output file
    logical :: write_to_file = .false.
    ! Variable index of the corresponding pool variable for uptake/release variables
    integer :: pool_index
    ! Variable index of the corresponding 'Ingested' variable for pool variables
    integer :: ingest_index
    ! Flag signalling path integration
    logical :: path_integration = .false.

    logical :: stage_diagnostic = .false.
  end type le_variable

  type food_variety
    ! Name is equivalent to the target's stage name
    character(len=FIELD_NAME_LEN) :: name

    ! Target agent list
    type(detector_list_ptr) :: target_list
    ! Target's concentration field
    character(len=FIELD_NAME_LEN) :: conc_field

    ! Indices for the associated request/ingest variables
    type(le_variable) :: vrequest, vingest
  end type food_variety

  type food_set
    character(len=FIELD_NAME_LEN) :: name
    ! Functional Group of our target
    character(len=FIELD_NAME_LEN) :: target_fgroup

    ! List of target varieties
    type(food_variety), dimension(:), allocatable :: varieties

    ! Variable indices of the chemical pools to ingest
    integer, dimension(:), allocatable :: ingest_chem_inds

    logical :: path_integrate = .false.
  end type food_set

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
     character(len=FIELD_NAME_LEN) :: stage_name
     character(len=OPTION_PATH_LEN) :: stage_options
     character(len=PYTHON_FUNC_LEN) :: biovar_pycode

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

     ! Periodic tracking
     logical :: track_periodic = .false.
     character(len=FIELD_NAME_LEN) :: periodic_tracking_mesh
     ! Map of boundary ID to the Python mapping function
     ! Note: This should not be here, but somewhere on the mesh_type...
     type(integer_hash_table) :: bid_to_boundary_mapping
     character(len=PYTHON_FUNC_LEN), dimension(:), allocatable :: boundary_mappings

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
     integer :: output_unit = 0          ! Assumed non-opened as long this is 0
     real :: output_period               ! Output period (simulation time)
     integer :: write_index = 0          ! Current position in output file
     integer :: total_num_det = 0        ! Global number of detectors in this list
  end type detector_linked_list

  type functional_group
    character(len=FIELD_NAME_LEN) :: name
    character(len=OPTION_PATH_LEN) :: option_path

    ! List of all stages within this FG
    type(stringlist) :: stage_names
    ! Agent arrays for each stage
    type(detector_linked_list), dimension(:), allocatable :: agent_arrays

    ! List of variables that define each agent of this group
    type(le_variable), dimension(:), allocatable :: variables
    ! List of aggregated agent sets to feed on
    type(food_set), dimension(:), allocatable :: food_sets

    ! List of environment fields to sample before the agent update
    character(len=FIELD_NAME_LEN), dimension(:), allocatable :: envfield_names
    logical, dimension(:), allocatable :: envfield_integrate

    integer, dimension(:), allocatable :: ivars_uptake, ivars_release

    ! Indices of all variables that are exposed to the Python motion functions
    integer, dimension(:), allocatable :: motion_var_inds
    ! Indices of all variables that are history buffers
    integer, dimension(:), allocatable :: history_var_inds

    ! Initialisation options
    character(len=OPTION_PATH_LEN) :: init_options
    ! Option path for the Agents diagnostic field
    character(len=OPTION_PATH_LEN) :: agents_field_path

    ! External FGroups to represent Eulerian predators
    logical :: is_external = .false.
    character(len=PYTHON_FUNC_LEN) :: external_python
  end type functional_group

end module detector_data_types
