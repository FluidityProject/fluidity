subroutine test_merge_tensors

  use merge_tensors
  use vector_tools
  use unittest_tools

  real, dimension(2, 2) :: tensor1, tensor2
  real, dimension(3, 3) :: tensor3, tensor4, tensor5, evecs
  real, dimension(3) :: evals
  real :: rand1, rand2
  integer :: i, j
  logical :: fail = .false., warn = .false.

  tensor3(1, :) = (/0.44444444444444443910E-10, -0.50465783670179635238E-10, -0.61349817493813233242E-10/)
  tensor3(2, :) = (/-0.50465783670179635238E-10, 0.44444444444444443910E-10, 0.76850819455089155588E-10/)
  tensor3(3, :) = (/-0.61349817493813233242E-10, 0.76850819455089155588E-10, 0.25000000000000007974E-04/)
  tensor4(1, :) = (/0.71326666461061381599E-09, 0.17461636435229589495E-09, -0.86585670292680812475E-08/)
  tensor4(2, :) = (/0.17461636435229589495E-09, 0.16746373209124310592E-09, -0.84640035445143090179E-09/)
  tensor4(3, :) = (/-0.86585670292680795931E-08, -0.84640035445143193577E-09, 0.25153686650982690250E-04/)
  call merge_tensor(tensor3, tensor4)
  ! Passing is not crashing.


  tensor1(1, :) = (/1.0, 0.0/)
  tensor1(2, :) = (/0.0, 2.0/)

  tensor2(1, :) = (/2.0, 0.0/)
  tensor2(2, :) = (/0.0, 1.0/)

  call merge_tensor(tensor1, tensor2)
  if (.not. fequals(tensor1(1, 1), 2.0)) fail = .true.
  if (.not. fequals(tensor1(1, 2), 0.0)) fail = .true.
  if (.not. fequals(tensor1(2, 1), 0.0)) fail = .true.
  if (.not. fequals(tensor1(2, 2), 2.0)) fail = .true.

  ! if eigenvectors are aligned, the eigenvalues of the result should just be the max of the two inputs
  call report_test("[merge tensors: 1]", fail, warn, "Merging of tensors is not 2*I.")

  fail = .false.
  warn = .false.

  tensor1(1, :) = (/1.0, 0.0/)
  tensor1(2, :) = (/0.0, 1.0/)

  tensor2(1, :) = (/3.0/2, 1.0/2/)
  tensor2(2, :) = (/1.0/2, 3.0/2/)

  ! merging with I should be the identity operation if minimum eigenvalue of A is >= 1
  call merge_tensor(tensor1, tensor2)

  if (.not. fequals(tensor1(1, 1), tensor2(1, 1))) fail = .true.
  if (.not. fequals(tensor1(1, 2), tensor2(1, 2))) fail = .true.
  if (.not. fequals(tensor1(2, 1), tensor2(2, 1))) fail = .true.
  if (.not. fequals(tensor1(2, 2), tensor2(2, 2))) fail = .true.

  call report_test("[merge tensors: 2]", fail, warn, &
    "If max eigenvalue of A is >= 1, merging with identity is the identity operation.")

  fail = .false.
  warn = .false.

  tensor3(1, :) = (/100.0, 0.0, 0.0/)
  tensor3(2, :) = (/0.0, 100.0, 0.0/)
  tensor3(3, :) = (/0.0, 0.0, 100.0/)

  tensor4 = tensor3
  tensor5 = tensor3

  call merge_tensor(tensor3, tensor4)

  do i=1,3
    do j=1,3
      if (.not. fequals(tensor3(i, j), tensor5(i, j))) fail = .true.
    end do
  end do
  
  call report_test("[merge tensors: 3]", fail, warn, "Merging identical metrics should yield the metric back.")

  fail = .false.
  warn = .false.

  tensor3 = random_posdef_matrix(3)
  tensor4 = tensor3
  tensor5 = tensor3

  call merge_tensor(tensor3, tensor4)

  do i=1,3
    do j=1,3
      if (.not. fequals(tensor3(i, j), tensor5(i, j))) fail = .true.
    end do
  end do

  call report_test("[merge tensors: 4]", fail, warn, "Merging identical metrics should yield the metric back.")

  fail = .false.
  warn = .false.

  tensor3(1, :) = (/1.00467, 0.00162196, -0.24383/)
  tensor3(2, :) = (/0.00162196, 6.2432, -1.98199/)
  tensor3(3, :) = (/-0.24383, -1.98199, 14.4226/)
  
  tensor4(1, :) = (/4.05543, 0.0983953, 0.643409/)
  tensor4(2, :) = (/0.0983953, 1.34818, -0.989943/)
  tensor4(3, :) = (/0.643409, -0.989943, 4.09606/)

  !call merge_tensor(tensor3, tensor4)

  tensor4(1, :) = (/3.11438, 1.37323, -2.03115/)
  tensor4(2, :) = (/1.37323, 10.2941, 10.1273/)
  tensor4(3, :) = (/-2.03115, 10.1273, 15.8707/)

  do i=1,3
    do j=1,3
      if (.not. fequals(tensor3(i, j), tensor4(i, j))) fail = .true.
    end do
  end do

  !call report_test("[merge tensors: 5]", fail, warn, "Previous bugs should not happen again.")
  
  fail = .false.
  warn = .false.

  tensor3(1, :) = (/4.42362, 0.409988, -0.491779/)
  tensor3(2, :) = (/0.409988, 5.55217, -2.27999/)
  tensor3(3, :) = (/-0.491779, -2.27999, 7.35012/)
  
  tensor4(1, :) = (/6.35559, -0.789674, -0.0779207/)
  tensor4(2, :) = (/-0.789674, 2.75034, 0.084052/)
  tensor4(3, :) = (/-0.0779207, 0.084052, 1.11065/)

  !call merge_tensor(tensor3, tensor4)

  tensor4(1, :) = (/17.0109, 3.62308, 1.0516/)
  tensor4(2, :) = (/3.62308, 4.66043, 5.45704/)
  tensor4(3, :) = (/1.0516, 5.45704, 9.01446/)

  do i=1,3
    do j=1,3
      if (.not. fequals(tensor3(i, j), tensor4(i, j))) fail = .true.
    end do
  end do

  !call report_test("[merge tensors: 6]", fail, warn, "Previous bugs should not happen again.")

  fail = .false.
  warn = .false.

  call random_number(rand1)
  tensor3 = get_mat_diag((/rand1, rand1, rand1/))

  call random_number(rand2)
  tensor4 = get_mat_diag((/rand2, rand2, rand2/))

  call merge_tensor(tensor3, tensor4)
  if (.not. mat_diag(tensor3)) fail = .true.

  rand1 = max(rand1, rand2)
  do i=1,3
    if (.not. fequals(tensor3(i, i), rand1)) fail = .true.
  end do

  call report_test("[merge tensors: 5]", fail, warn, "Merging diagonal matrices should yield a diagonal matrix back.")

  fail = .false.
  warn = .false.
  tensor3(1, :) = (/0.10017100417445532479E+01, -0.16978728026906123386E-03, -0.37497787153930919518E-05/)
  tensor3(2, :) = (/-0.16978728026906123386E-03, 0.10000168579045705108E+01,  0.37230946659257456045E-06/)
  tensor3(3, :) = (/-0.37497787153935256327E-05, 0.37230946659257456045E-06,  0.10000000082225128928E+01/)

  call eigendecomposition_symmetric(tensor3, evecs, evals)

  tensor4 = get_matrix_identity(3)
  tensor5 = tensor3

  call merge_tensor(tensor3, tensor4)
  call eigendecomposition_symmetric(tensor3, evecs, evals)

  do i=1,3
    do j=1,3
      if (.not. fequals(tensor3(i, i), tensor5(i, i))) fail = .true.
    end do
  end do

  call report_test("[merge tensors: 6]", fail, warn, "Merging tensors with all eigenvalues > 1 with the identity should be the &
                   & identity operation.")

  tensor1(1, :) = (/1.0, 0.0/)
  tensor1(2, :) = (/0.0, 1.0/)
  tensor2(1, :) = (/0.0, 0.0/)
  tensor2(2, :) = (/0.0, 0.0/)

  call merge_tensor(tensor1, tensor2)

  fail = .false.
  if (tensor1 .fne. get_mat_diag((/1.0, 1.0/))) then
    fail = .true.
  end if

  call report_test("[merge tensors: 7]", fail, warn, "Merging a tensor with the zero tensor should be &
                   & a no-op.")

  tensor1(1, :) = (/0.0, 0.0/)
  tensor1(2, :) = (/0.0, 0.0/)
  tensor2(1, :) = (/1.0, 0.0/)
  tensor2(2, :) = (/0.0, 1.0/)

  call merge_tensor(tensor1, tensor2)

  fail = .false.
  if (tensor1 .fne. get_mat_diag((/1.0, 1.0/))) then
    fail = .true.
  end if

  call report_test("[merge tensors: 8]", fail, warn, "Merging the zero tensor should be assignment.")
end subroutine test_merge_tensors
