module zoltan_global_variables

  ! Needed for zoltan_cb_owned_node_count
  use halos, only: halo_type

  ! Needed for zoltan_cb_get_owned_nodes
  use sparse_tools, only: csr_sparsity

  implicit none

  public

  ! Needed for zoltan_cb_owned_node_count
  type(halo_type), save, pointer :: zoltan_global_zz_halo

  ! Needed for zoltan_cb_get_owned_nodes
  type(csr_sparsity), save :: zoltan_global_columns_sparsity
  logical, save :: zoltan_global_migrate_extruded_mesh

  ! Needed for zoltan_cb_get_num_edges
  type(csr_sparsity), pointer :: zoltan_global_zz_sparsity_one

end module zoltan_global_variables
