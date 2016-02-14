
#include "confdefs.h"

module zoltan_global_variables

#ifdef HAVE_ZOLTAN

  ! Needed for zoltan_cb_owned_node_count
  use halos, only: halo_type

  ! Needed for zoltan_cb_get_owned_nodes
  use sparse_tools, only: csr_sparsity

  ! Needed for zoltan_cb_get_edge_list
  use fields, only: scalar_field, vector_field, mesh_type

  ! Needed for zoltan_cb_pack_node_sizes
  use zoltan, only: zoltan_int, zoltan_float
  use data_structures, only: integer_set, integer_hash_table

  ! Needed for zoltan_cb_pack_field_size
  use state_module, only: state_type
  use detector_data_types, only: detector_linked_list

  use global_parameters, only: OPTION_PATH_LEN

  implicit none

  public

  ! Needed for zoltan_cb_owned_node_count
  type(halo_type), save, pointer :: zoltan_global_zz_halo

  ! Needed for zoltan_cb_get_owned_nodes
  type(csr_sparsity), save :: zoltan_global_columns_sparsity
  logical, save :: zoltan_global_migrate_extruded_mesh
  logical, save :: zoltan_global_field_weighted_partitions
  type(scalar_field), save :: zoltan_global_field_weighted_partition_values

  ! Needed for zoltan_cb_get_num_edges
  type(csr_sparsity), save, pointer :: zoltan_global_zz_sparsity_one

  ! Needed for zoltan_cb_get_edge_list
  logical, save :: zoltan_global_calculate_edge_weights
  ! elements with quality greater than this value are ok
  ! those with element quality below it need to be adapted
  real, save :: zoltan_global_quality_tolerance
  type(scalar_field), save :: zoltan_global_element_quality
  type(scalar_field), save, pointer :: zoltan_global_max_edge_weight_on_node
  logical, save :: zoltan_global_output_edge_weights = .false.
  type(csr_sparsity), save, pointer :: zoltan_global_zz_nelist

  ! Needed for zoltan_cb_pack_node_sizes
  ! - added vector_field to use fields
  type(vector_field), save :: zoltan_global_zz_positions
  integer, parameter :: integer_size = bit_size(0_zoltan_int)/8
  logical, save :: zoltan_global_preserve_columns=.false.
  logical, save :: zoltan_global_preserve_mesh_regions
  type(csr_sparsity), save, pointer :: zoltan_global_zz_sparsity_two
  type(integer_set), save, dimension(:), allocatable :: zoltan_global_old_snelist


  ! Needed for zoltan_cb_pack_nodes
  type(integer_hash_table), save :: zoltan_global_universal_element_number_to_region_id
  type(integer_hash_table), save :: zoltan_global_universal_surface_number_to_element_owner
  type(integer_hash_table), save :: zoltan_global_universal_surface_number_to_surface_id
  integer, dimension(:), allocatable, save :: zoltan_global_universal_columns
  type(halo_type), save, pointer :: zoltan_global_zz_ele_halo 


  ! Needed for zoltan_cb_unpack_nodes
  type(vector_field), save :: zoltan_global_new_positions
  integer, save :: zoltan_global_new_positions_mesh_nhalos
  type(mesh_type), save :: zoltan_global_zz_mesh
  type(integer_hash_table), save :: zoltan_global_nodes_we_are_sending ! in old local numbers
  type(integer_set), save :: zoltan_global_nodes_we_are_keeping ! in old local numbers
  type(integer_hash_table), save :: zoltan_global_universal_to_new_local_numbering
  type(integer_hash_table), save :: zoltan_global_universal_to_old_local_numbering
  type(integer_set), save :: zoltan_global_new_nodes
  type(integer_hash_table), save :: zoltan_global_universal_to_new_local_numbering_m1d
  type(integer_set), save, dimension(:), allocatable :: zoltan_global_new_snelist
  type(integer_set), save :: zoltan_global_new_surface_elements
  type(integer_set), save, dimension(:), allocatable :: zoltan_global_new_nelist
  type(integer_set), save :: zoltan_global_new_elements
  integer(zoltan_int), save, dimension(:), pointer :: zoltan_global_my_import_procs => null()
  integer(zoltan_int), save, dimension(:), pointer :: zoltan_global_my_import_global_ids => null()
  integer(zoltan_int), save :: zoltan_global_my_num_import
  type(integer_set), save, dimension(:), allocatable :: zoltan_global_receives


  ! Needed for prepare_detectors_for_packing
  type(integer_hash_table), save :: zoltan_global_uen_to_new_local_numbering
  type(integer_hash_table), save :: zoltan_global_uen_to_old_local_numbering
  type(integer_hash_table), save :: zoltan_global_old_local_numbering_to_uen


  ! Needed for zoltan_cb_pack_field_sizes
  type(state_type), save, dimension(:), allocatable :: zoltan_global_source_states, zoltan_global_target_states
  integer, save, dimension(:), allocatable :: zoltan_global_ndets_in_ele  
  integer, save :: zoltan_global_ndata_per_det
  type(detector_linked_list), dimension(:), allocatable, target, save :: zoltan_global_to_pack_detectors_list

  ! Needed for zoltan_cb_pack_fields
  integer, save :: zoltan_global_ndims

  ! Needed for zoltan_cb_unpack_fields
  type(detector_linked_list), target, save :: zoltan_global_unpacked_detectors_list

  ! Option path set based on whether being called from adaptivity or flredecomp
  character(len = OPTION_PATH_LEN), save :: zoltan_global_base_option_path

#endif

end module zoltan_global_variables
