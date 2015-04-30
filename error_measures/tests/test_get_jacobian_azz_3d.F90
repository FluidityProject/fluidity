subroutine test_get_jacobian_azz_3d

  use fields
  use mesh_files
  use anisotropic_zz_module
  use vector_tools
  use metric_tools
  use unittest_tools
  use limit_metric_module
  use meshdiagnostics
  implicit none

  type(vector_field) :: positions
  real, dimension(3, 3) :: j_k, invj_k
  real, dimension(3, 4) :: new_pos, old_pos
  logical :: fail
  real :: edge_length, ideal_edge_length
  integer :: iloc, jloc
  real :: transformed_volume
  real :: ideal_volume
  real, dimension(3) :: t_k
  integer :: i

  positions = read_mesh_files("data/tet", quad_degree=4, format="gmsh")
  old_pos = ele_val(positions, 1)
  j_k = get_jacobian_azz(positions, 1)
  invj_k = inverse(j_k)
  do i=1,3
    t_k(i) = sum(old_pos(i, :))
  end do
  t_k = t_k/4.0
  do i=1,4
    old_pos(:, i) = old_pos(:, i) - t_k
  end do
  new_pos = apply_transform(ele_val(positions, 1), invj_k)

  ideal_edge_length = 2 * sqrt(2.0/3.0)
  fail = .false.
  do iloc=1,ele_loc(positions, 1)
    do jloc=iloc+1,ele_loc(positions, 1)
      edge_length = norm2(new_pos(:, iloc) - new_pos(:, jloc))
      fail = fail .or. (edge_length .fne. ideal_edge_length)
    end do
  end do

  call report_test("[get_jacobian_azz_3d edge lengths]", fail, .false., "Edge lengths should be 2*sqrt(2.0/3.0)")

  transformed_volume = simplex_volume(positions, 1) * determinant(invj_k)
  ideal_volume = 8.0 / (9.0 * sqrt(3.0))
  fail = (transformed_volume .fne. ideal_volume)

  call report_test("[get_jacobian_azz_3d volume]", fail, .false., "Volume should be 8.0/(9.0*sqrt(3.0))")

end subroutine test_get_jacobian_azz_3d
