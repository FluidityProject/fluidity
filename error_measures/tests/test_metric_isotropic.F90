subroutine test_metric_isotropic

  use metric_tools
  use unittest_tools
  use vector_tools
  implicit none

  real, dimension(3, 3) :: mat, vecs
  real, dimension(3) :: vals, isotropic_vals, anisotropic_vals
  logical :: fail

  mat = random_posdef_matrix(3)
  call eigendecomposition_symmetric(mat, vecs, vals)

  isotropic_vals = 1.0
  call eigenrecomposition(mat, vecs, isotropic_vals)

  fail = .false.
  if (.not. metric_isotropic(mat)) fail = .true.
  call report_test("[isotropic metric]", fail, .false., "isotropic_metric should &
    & report true for an isotropic metric.")

  anisotropic_vals(1) = 1.0
  anisotropic_vals(2) = 2.0
  anisotropic_vals(3) = 3.0
  call eigenrecomposition(mat, vecs, anisotropic_vals)

  fail = .false.
  if (metric_isotropic(mat)) fail = .true.
  call report_test("[anisotropic metric]", fail, .false., "isotropic_metric should &
    & report false for an anisotropic metric.")
end subroutine test_metric_isotropic
