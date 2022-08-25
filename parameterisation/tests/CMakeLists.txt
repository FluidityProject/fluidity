set(unittest_sources
  test_pade_equation_of_state.F90
  test_tensor_inner_product.F90
)

cmake_language(CALL build_unittest_executable "${unittest_sources}")
