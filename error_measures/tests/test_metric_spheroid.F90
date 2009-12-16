subroutine test_metric_spheroid

  use metric_tools
  use unittest_tools
  use vector_tools
  implicit none

  real, dimension(3, 3) :: mat, vecs
  real, dimension(3) :: vals, spheroid_vals, nonspheroid_vals
  integer :: idx
  logical :: fail

  mat = random_posdef_matrix(3)
  call eigendecomposition_symmetric(mat, vecs, vals)

  spheroid_vals = 1.0; spheroid_vals(1) = 0.5
  call eigenrecomposition(mat, vecs, spheroid_vals)

  fail = .false.
  if (metric_spheroid(mat)) fail = .true.
  call report_test("[spheroid metric]", fail, .false., "spheroid_metric should &
   & report false for a flattened spheroid metric.")

  fail = .false.
  idx = get_spheroid_index(mat, vecs, spheroid_vals)
  if (idx == 1) fail = .true.
  call report_test("[spheroid index]", fail, .false., "get_spheroid_index should &
   & not return 1.")

  fail = .false.
  idx = get_polar_index(spheroid_vals)
  if (idx /= 1) fail = .true.
  call report_test("[polar index]", fail, .false., "get_polar_index should &
   & not return 1.")

  ! ----------

  mat = random_posdef_matrix(3)
  call eigendecomposition_symmetric(mat, vecs, vals)

  spheroid_vals = 1.0; spheroid_vals(1) = 1.5
  call eigenrecomposition(mat, vecs, spheroid_vals)

  fail = .false.
  if (.not. metric_spheroid(spheroid_vals)) fail = .true.
  call report_test("[spheroid metric]", fail, .false., "spheroid_metric should &
   & report true for a spheroid metric.")

  fail = .false.
  idx = get_spheroid_index(mat, vecs, spheroid_vals)
  if (idx == 1) fail = .true.
  call report_test("[spheroid index]", fail, .false., "get_spheroid_index should &
   & not return 1.")

  fail = .false.
  idx = get_polar_index(spheroid_vals)
  if (idx /= 1) fail = .true.
  call report_test("[polar index]", fail, .false., "get_polar_index should &
   & not return 1.")

  ! -----------

  nonspheroid_vals(1) = 1.0
  nonspheroid_vals(2) = 2.0
  nonspheroid_vals(3) = 3.0
  call eigenrecomposition(mat, vecs, nonspheroid_vals)

  fail = .false.
  if (metric_spheroid(mat)) fail = .true.
  call report_test("[nonspheroid metric]", fail, .false., "spheroid_metric should &
   & report false for a nonspheroid metric.")
end subroutine test_metric_spheroid
