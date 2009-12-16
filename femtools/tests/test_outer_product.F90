subroutine test_outer_product

  use vector_tools
  use unittest_tools
  implicit none

  real, dimension(4) :: ones = 1.0
  real, dimension(4, 4) :: output, correct
  logical :: fail

  correct = 1.0

  output = 0.0
  call outer_product(ones, ones, output)

  fail = .false.
  if (any(output /= correct)) fail = .true.

  call report_test("[outer product]", fail, .false., "Outer product should give known good values.")

end subroutine test_outer_product
