#include "fdebug.h"

module recovery_estimator
!!< This module implements the class of recovery-based
!!< error estimators (see Ainesworth & Oden, chapter 4).
!!< Here the error is estimated by comparing the original
!!< direct derivative of the gradient and the post-processed
!!< derivative computed with some recovery technique.
!!< The elementwise error is then estimated as
!!< int(|G(u_x) - grad(u_x)|^2) over the element.

  use transform_elements, only: transform_to_physical
  use fields
  use field_derivatives, only: grad
  implicit none

  contains

  subroutine form_recovery_estimator(infield, positions, estimator)
    type(scalar_field), intent(in) :: infield
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(inout) :: estimator ! piecewise constant basis

    type(vector_field) :: grad_infield
    real, dimension(ele_ngi(infield, 1)) :: detwei
    real, dimension(ele_loc(infield, 1), ele_ngi(infield, 1), mesh_dim(infield)) :: dt_t
    real, dimension(ele_ngi(infield, 1)) :: r
    type(element_type), pointer :: t_shape
    integer :: ele, dim

    call zero(estimator)

    call allocate(grad_infield, mesh_dim(infield), infield%mesh, "Gradient")
    call grad(infield, positions, grad_infield)

    t_shape => ele_shape(infield, 1)

    do ele=1,element_count(infield)
      call transform_to_physical(positions, ele, t_shape, dshape=dt_t, detwei=detwei)
      do dim=1,mesh_dim(infield)
        r = (matmul(ele_val(grad_infield, dim, ele), t_shape%n) - matmul(ele_val(infield, ele), dt_t(:, :, dim)))**2
        call addto(estimator, ele, dot_product(r, detwei))
      end do
    end do

    call deallocate(grad_infield)
  end subroutine form_recovery_estimator

end module recovery_estimator
