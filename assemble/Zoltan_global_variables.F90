module zoltan_global_variables

  ! Needed for zoltan_cb_owned_node_count
  use halos, only: halo_type

  ! Needed for zoltan_cb_get_owned_nodes
  use sparse_tools, only: csr_sparsity

  ! Needed for zoltan_cb_get_edge_list
  use fields, only: scalar_field, vector_field

  ! Needed for zoltan_cb_pack_node_sizes
  use zoltan, only: zoltan_int
  use data_structures, only: integer_set, integer_hash_table

  implicit none

  public

  ! Needed for zoltan_cb_owned_node_count
  type(halo_type), save, pointer :: zoltan_global_zz_halo

  ! Needed for zoltan_cb_get_owned_nodes
  type(csr_sparsity), save :: zoltan_global_columns_sparsity
  logical, save :: zoltan_global_migrate_extruded_mesh

  ! Needed for zoltan_cb_get_num_edges
  type(csr_sparsity), save, pointer :: zoltan_global_zz_sparsity_one

  ! Needed for zoltan_cb_get_edge_list
  integer, save :: zoltan_global_zoltan_iteration, zoltan_global_zoltan_max_adapt_iteration
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

end module zoltan_global_variables
