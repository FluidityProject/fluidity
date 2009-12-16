subroutine test_project_to_subspace

  use metric_tools
  use unittest_tools
  implicit none

  real, dimension(3) :: a, b, c, d, e
  real, dimension(3, 2) :: basis
  logical :: fail
  integer :: i
  character(len=20) :: buf

  a = (/1.0, 1.0, 1.0/) ! the vector to project
  b = (/1.0, 0.0, 0.0/) ! basis vectors
  c = (/0.0, 1.0, 0.0/)
  e = (/1.0, 1.0, 0.0/) ! the correct answer 

  basis(:, 1) = b; basis(:, 2) = c

  d = project_to_subspace(a, basis)

  fail = .false.
  if (d .fne. e) fail = .true.
  call report_test("[project to subspace known values]", fail, .false., "Projecting &
                   & should give known good values.")

  a = (/0.0, 0.0, 1.0/) ! a is orthogonal to basis
  e = (/0.0, 0.0, 0.0/)

  d = project_to_subspace(a, basis)

  fail = .false.
  if (d .fne. e) fail = .true.
  call report_test("[project to subspace known values]", fail, .false., "Projecting &
                   & should give known good values.")

  a = (/1.0, 1.0, 0.0/) ! this is now in the subspace
  e = (/1.0, 1.0, 0.0/)

  d = project_to_subspace(a, basis)

  fail = .false.
  if (d .fne. e) fail = .true.
  call report_test("[project to subspace known values]", fail, .false., "Projecting &
                   & should give known good values.")

  do i=1,5
    write(buf, '(i0)') i

    a = random_vector(3)
    b = random_vector(3) ; b = b / norm(b)
    c = random_vector(3)
    c = c - dot_product(c, b) * b ; c = c / norm(c)
    basis(:, 1) = b; basis(:, 2) = c

    d = project_to_subspace(a, basis)
    e = project_to_subspace(d, basis)

    fail = .false.
    if (d .fne. e) fail = .true.
    call report_test("[project to subspace idempotent " // trim(buf) // "]", fail, .false., &
                     "Projecting to a subspace is an idempotent linear operation.")
  end do
end subroutine test_project_to_subspace
