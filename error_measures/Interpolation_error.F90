#include "fdebug.h"

module interpolation_error

  use transform_elements, only: transform_to_physical
  use fields
  use vtk_interfaces, only: vtk_write_fields
  implicit none

  contains

  function compute_interpolation_error_l2(solution, field, positions) result(l2)
    !!< Compute the l2-norm of the error between the function solution
    !!< and its interpolant field.
    interface
      function solution(pos)
        real, dimension(:), intent(in) :: pos
        real :: solution
      end function solution
    end interface
    type(scalar_field), intent(in) :: field
    type(vector_field), intent(in) :: positions
    type(element_type), pointer :: shape_field
    real :: l2, ele_int
    type(scalar_field) :: debug

    integer :: ele
    real, dimension(ele_ngi(positions,1)) :: detwei

    debug = piecewise_constant_field(field%mesh, "L2 error")

    l2 = 0.0
    shape_field => ele_shape(field, 1)

    do ele=1,element_count(field)
      call transform_to_physical(positions, ele, detwei=detwei)
      ele_int = dot_product((abs(function_val_at_quad_scalar(solution, positions, ele) - ele_val_at_quad(field, ele)))**2, detwei)
      debug%val(ele) = ele_int
      l2 = l2 + ele_int
    end do
    l2 = sqrt(l2)

    call vtk_write_fields("l2_error", 0, positions, field%mesh, sfields=(/debug/))
  
  end function compute_interpolation_error_l2

  function compute_interpolation_error_inf(solution, field, positions) result(maxn)
    !!< Compute the inf-norm of the error between the function solution
    !!< and its interpolant field.
    interface
      function solution(pos)
        real, dimension(:), intent(in) :: pos
        real :: solution
      end function solution
    end interface
    type(scalar_field), intent(in) :: field
    type(vector_field), intent(in) :: positions
    real :: maxn, current_max
    type(scalar_field) :: debug

    integer :: ele

    debug = piecewise_constant_field(field%mesh, "Inf error")

    maxn = 0.0

    do ele=1,element_count(field)
      current_max = maxval(abs(function_val_at_quad_scalar(solution, positions, ele) - ele_val_at_quad(field, ele)))
      debug%val(ele) = current_max
      maxn = max(maxn, current_max)
    end do

    call vtk_write_fields("inf_error", 0, positions, field%mesh, sfields=(/debug/))
  
  end function compute_interpolation_error_inf

  function compute_interpolation_error_h1(gradsoln, field, positions) result(h1)
    !!< Compute the h1-norm of the error between the function gradsoln
    !!< and its interpolant field.
    interface
      function gradsoln(pos)
        real, dimension(:), intent(in) :: pos
        real, dimension(size(pos)) :: gradsoln
      end function gradsoln
    end interface
    type(scalar_field), intent(in) :: field
    type(vector_field), intent(in) :: positions
    type(element_type), pointer :: shape_field
    real :: h1, ele_int
    type(scalar_field) :: debug
    real, dimension(ele_loc(field,1), ele_ngi(field, 1), positions%dim) :: dm_t
    real, dimension(positions%dim, ele_ngi(field, 1)) :: grad_at_quad
    real, dimension(ele_ngi(field, 1)) :: err_at_quad
    integer :: dim


    integer :: ele
    real, dimension(ele_ngi(positions,1)) :: detwei

    debug = piecewise_constant_field(field%mesh, "H1 error")

    h1 = 0.0
    shape_field => ele_shape(field, 1)

    do ele=1,element_count(field)
      grad_at_quad = function_val_at_quad(gradsoln, positions, ele)
      call transform_to_physical(positions, ele, shape_field, dshape=dm_t, detwei=detwei)
      ele_int = 0.0
      do dim=1,positions%dim
        err_at_quad = (grad_at_quad(dim, :) - matmul(ele_val(field, ele), dm_t(:, :, dim)))**2
        ele_int = ele_int + dot_product(err_at_quad, detwei)
      end do
      debug%val(ele) = ele_int
      h1 = h1 + ele_int
    end do
    h1 = sqrt(h1)

    call vtk_write_fields("h1_error", 0, positions, field%mesh, sfields=(/debug/))
    call deallocate(debug)
  
  end function compute_interpolation_error_h1

end module interpolation_error
