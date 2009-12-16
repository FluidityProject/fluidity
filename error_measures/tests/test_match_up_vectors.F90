subroutine test_match_up_vectors

  use metric_tools
  use gradation_metric
  use unittest_tools
  use vector_tools
  implicit none

  real, dimension(3, 3) :: P, Q, vec_P, vec_Q
  real, dimension(3) :: eigenvals
  integer, dimension(3) :: permutation_P, permutation_Q
  real, dimension(3, 1) :: vec_P_short
  integer, dimension(1) :: perm_P_short
  integer :: stat

  logical :: fail
  integer :: i, j
  character(len=20) :: buf

  vec_P = get_matrix_identity(3)

  vec_Q(:, 1) = (/1.0, 0.0, 0.0/)
  vec_Q(:, 2) = (/0.0, 1.0, -0.1/) ; vec_Q(:, 2) = vec_Q(:, 2) / norm(vec_Q(:, 2))
  vec_Q(:, 3) = (/0.0, 0.1, 1.0/) ; vec_Q(:, 3) = vec_Q(:, 3) / norm(vec_Q(:, 3))

  call eigenrecomposition(P, vec_P, (/1.0, 1.0, 1.0/))
  call eigenrecomposition(Q, vec_Q, (/1.0, 1.0, 1.0/))

  call match_up_vectors(vec_P, permutation_P, vec_Q, permutation_Q)

  fail = .false.
  do j=1,3
    if (permutation_P(j) /= permutation_Q(j)) fail = .true.
  end do
  call report_test("[match up known vectors]", fail, .false., "Any reasonable &
                   & algorithm for matching up vectors must yield the identity  &
                   & permutation in this case.")

  vec_Q = get_matrix_identity(3)
  vec_Q(:, 3) =  vec_Q(:, 2)
  vec_Q(:, 2) =  vec_Q(:, 1)
  vec_Q(:, 1) = (/0.0, 0.0, 1.0/)

  call match_up_vectors(vec_P, permutation_P, vec_Q, permutation_Q)

  fail = .false.
  do j=1,3
    if (permutation_P(j) == 1 .and. permutation_Q(j) /= 3) fail = .true.
    if (permutation_P(j) == 2 .and. permutation_Q(j) /= 1) fail = .true.
    if (permutation_P(j) == 3 .and. permutation_Q(j) /= 2) fail = .true.
  end do
  call report_test("[match up known vectors]", fail, .false., "Any reasonable &
                     & algorithm for matching up vectors must yield the (3, 1, 2)  &
                     & permutation in this case.")

  vec_P_short(:, 1) = (/0.0, 0.0, 1.0/)
  vec_Q = get_matrix_identity(3)
  call match_up_vectors(vec_P_short, perm_P_short, vec_Q, permutation_Q)
  fail = .false.
  if (perm_P_short(1) /= 3) fail = .true.
  if (permutation_Q(1) /= 1) fail = .true.
  call report_test("[match up known vectors]", fail, .false., "Match up vectors &
                   & should work for just one vector, too.")

  do i=1,5
    write(buf, '(i0)') i

    P = random_posdef_matrix(3)
    Q = random_posdef_matrix(3)
    call eigendecomposition(P, vec_P, eigenvals)
    call eigendecomposition(Q, vec_Q, eigenvals)

    call match_up_vectors(vec_P, permutation_P, vec_Q, permutation_Q)

    fail = .false.
    if (sum(permutation_P) /= 6) fail = .true. ! 1 + 2 + 3 = 6
    if (sum(permutation_Q) /= 6) fail = .true. ! 1 + 2 + 3 = 6
    if (product(permutation_P) /= 6) fail = .true. ! 1 * 2 * 3 = 6
    if (product(permutation_Q) /= 6) fail = .true. ! 1 * 2 * 3 = 6
    call report_test("[match up vectors gives a permutation " // trim(buf) // "]", fail, .false., "Matching up &
                     & vectors should give a permutation.")
  end do

  vec_P(1, :) = (/0.96864848337775710796E+00, -0.24842986328315525002E+00, -0.16488417400240069823E-02/)
  vec_P(2, :) = (/0.24843533494648231685E+00, 0.96862714938199356851E+00, 0.64288280535788229486E-02/)
  vec_P(3, :) = (/0.00000000000000000000E+00, -0.66369050939521093482E-02, 0.99997797550284772683E+00/)

  vec_Q(1, :) = (/0.96864848337775721898E+00, -0.24842986326122473706E+00, -0.16488450442796352084E-02/)
  vec_Q(2, :) = (/0.24843533494648234461E+00, 0.96862714929648652262E+00, 0.64288409368597095039E-02/)
  vec_Q(3, :) = (/0.00000000000000000000E+00, -0.66369183942164558174E-02, 0.99997797541457311699E+00/)

  fail = .false.
  call match_up_vectors(vec_P, permutation_P, vec_Q, permutation_Q)
  call check_perm(permutation_P, stat)
  if (stat /= 0) then
    fail = .true.
    call write_vector(permutation_P, "perm_P")
  end if

  call check_perm(permutation_Q, stat)
  if (stat /= 0) then
    fail = .true.
    call write_vector(permutation_Q, "perm_Q")
  end if

  call report_test("[match up vectors regression]", fail, .false., "Previous bugs &
  & shouldn't happen again.")

end subroutine test_match_up_vectors
