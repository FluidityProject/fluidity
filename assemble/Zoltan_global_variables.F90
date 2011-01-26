module zoltan_global_variables

  use halos

  public

  type(halo_type), save, pointer :: zoltan_global_zz_halo

end module zoltan_global_variables
