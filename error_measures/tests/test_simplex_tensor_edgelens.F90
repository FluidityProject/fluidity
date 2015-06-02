subroutine test_simplex_tensor_edgelens

  use fields
  use mesh_files
  use transform_elements
  use meshdiagnostics
  use conformity_measurement
  use vector_tools
  use unittest_tools
  use metric_tools
  implicit none

  integer :: ele
  type(vector_field) :: positions

  positions = read_mesh_files("data/laplacian_grid.3", quad_degree=4, format="gmsh")

  do ele=1,ele_count(positions)
    call transform_ele(positions, ele)
  end do

  contains
    subroutine transform_ele(positions, ele)
      type(vector_field), intent(inout) :: positions
      integer, intent(in) :: ele

      real, dimension(positions%dim, positions%dim) :: tens
      real, dimension(positions%dim, ele_loc(positions, ele)) :: pos, t_pos

      integer :: loc, iloc, jloc

      logical :: fail
      real :: edge_length

      pos = ele_val(positions, ele)
      tens = simplex_tensor(positions, ele, power=0.5)

      do loc=1,ele_loc(positions, ele)
        t_pos(:, loc) = matmul(tens, pos(:, loc))
      end do

      fail = .false.
      do iloc=1,ele_loc(positions, ele)
        do jloc=iloc+1,ele_loc(positions, ele)
          edge_length = norm2(t_pos(:, iloc) - t_pos(:, jloc))
          fail = fail .or. fnequals(edge_length, 1.0, tol = 200.0 * epsilon(0.0))
        end do
      end do

      call report_test("[simplex_tensor]", fail, .false., "Edge lengths should be 1, na ja?")
    end subroutine

end subroutine test_simplex_tensor_edgelens
