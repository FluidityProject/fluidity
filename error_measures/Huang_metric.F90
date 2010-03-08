#include "fdebug.h"

module huang_metric
! 10.1016/j.jcp.2004.10.024

  use fields
  use state_module
  use vtk_interfaces
  use form_metric_field
  use vector_tools
  use global_parameters, only: domain_volume
  use fldebug
  use spud
  use metric_tools
  implicit none

  contains

  subroutine form_huang_metric(hessian, field, weight, state)
    type(tensor_field), intent(inout) :: hessian
    type(scalar_field), intent(in) :: field, weight
    type(state_type), intent(in) :: state

    integer :: n, m, l, p, q
    integer, dimension(2) :: tmp
    real :: gamma, rho, sigma, alpha, eps, power_A, power_B, power_C
    integer :: node
    real, dimension(hessian%dim, hessian%dim) :: identity, abs_h, metric

    assert(field%mesh%shape%degree == 1) ! this can be generalised, but currently I don't have the time
    assert(weight%field_type == FIELD_TYPE_CONSTANT) ! this can't

    n = mesh_dim(field)
    l = field%mesh%shape%degree + 1
    p = 1

    call get_option(trim(field%option_path) // "/adaptivity_options/huang_metric/seminorm", tmp)
    m = tmp(1)
    q = tmp(2)

    alpha = 1.0 ! FIXME :-(
    eps = node_val(weight, 1)
    gamma = (float(n) / q) + (2 - m)
    power_A = 1 + float((n * (p - 1)))/(p * gamma) + max(0.0, float(n)/(p*gamma) - 1)
    sigma = (2 ** power_A) * domain_volume

    identity = get_matrix_identity(n)
    power_B = (1.0/n) * ((2.0/gamma) - 1)
    power_C = float(n) / (2 - m)

    do node=1,node_count(hessian)
      abs_h = absolutify(node_val(hessian, node))

      metric = ((1.0/sigma) * (alpha / eps)**power_C)**(2.0/n) * &
               (det(identity + (1.0/alpha) * abs_h) ** power_B) * &
               (identity + (1.0/alpha) * abs_h)

      call set(hessian, node, metric)
    end do

    call bound_metric(hessian, state)

    contains 

      function absolutify(hessian) result(abs_h)
        real, dimension(:, :), intent(in) :: hessian
        real, dimension(size(hessian, 1), size(hessian, 2)) :: abs_h, evecs
        real, dimension(size(hessian, 1)) :: evals
        
        call eigendecomposition_symmetric(hessian, evecs, evals)
        evals = abs(evals)
        call eigenrecomposition(abs_h, evecs, evals)
      end function absolutify
  end subroutine form_huang_metric

end module huang_metric
