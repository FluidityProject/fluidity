subroutine test_aspect_ratio

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
  call eigenrecomposition(mat, vecs, eigenvalue_from_edge_length(isotropic_vals))

  fail = .false.
  if (.not. fequals(aspect_ratio(mat), 1.0)) fail = .true.
  call report_test("[aspect ratio 1]", fail, .false., "Aspect ratio known value.")

  anisotropic_vals(1) = 1.0
  anisotropic_vals(2) = 2.0
  anisotropic_vals(3) = 3.0
  call eigenrecomposition(mat, vecs, eigenvalue_from_edge_length(anisotropic_vals))

  fail = .false.
  if (aspect_ratio(anisotropic_vals) .fne. 1.0/sqrt(3.0)) fail = .true.
  call report_test("[aspect ratio 2]", fail, .false., "Aspect ratio known value.")
end subroutine test_aspect_ratio
