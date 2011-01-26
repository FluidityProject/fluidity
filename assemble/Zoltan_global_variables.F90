module zoltan_global_variables

  use halos, only: halo_type

  implicit none

  public

  type(halo_type), save, pointer :: zoltan_global_zz_halo

end module zoltan_global_variables
