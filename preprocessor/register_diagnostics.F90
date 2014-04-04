
subroutine register_diagnostics

  use implicit_solids, only: implicit_solids_register_diagnostic
  use linear_shallow_water, only: linear_shallow_water_register_diagnostic
  use forward_main_loop, only: forward_main_loop_register_diagnostic

   continue
     call implicit_solids_register_diagnostic
  call linear_shallow_water_register_diagnostic
  call forward_main_loop_register_diagnostic

end subroutine register_diagnostics
