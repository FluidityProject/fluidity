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

  type(tensor_field), save :: m_hessian
  type(vector_field), pointer :: m_coordinate
  real :: m_gamma, m_sigma

  private
  public :: form_huang_metric

  contains

  subroutine form_huang_metric(hessian, field, weight, state)
    type(tensor_field), intent(inout) :: hessian
    type(scalar_field), intent(in) :: field, weight
    type(state_type), intent(in) :: state

    integer :: n, m, l, p, q
    integer, dimension(2) :: tmp
    real :: gamma, sigma, alpha, eps, power_A, power_B, power_C
    integer :: node
    real, dimension(hessian%dim, hessian%dim) :: identity, abs_h, metric

    assert(field%mesh%shape%degree == 1) ! this can be generalised, but currently I don't have the time
    assert(weight%field_type == FIELD_TYPE_CONSTANT) ! this can't

    m_coordinate => extract_vector_field(state, "Coordinate")

    n = mesh_dim(field)
    l = field%mesh%shape%degree + 1
    p = 1

    call get_option(trim(field%option_path) // "/adaptivity_options/huang_metric/seminorm", tmp)
    m = tmp(1)
    q = tmp(2)

    eps = node_val(weight, 1)
    gamma = (float(n) / q) + (2 - m)
    power_A = 1 + float((n * (p - 1)))/(p * gamma) + max(0.0, float(n)/(p*gamma) - 1)
    sigma = (2 ** power_A) * domain_volume
    m_sigma = sigma

    identity = get_matrix_identity(n)
    power_B = (1.0/n) * ((2.0/gamma) - 1)
    power_C = float(n) / (2 - m)

    do node=1,node_count(hessian)
      abs_h = absolutify(node_val(hessian, node))
      call set(hessian, node, abs_h)
    end do

    m_hessian = hessian
    m_gamma = gamma

    alpha = secant_method(beta, 1.0e-6, 1.0e6, 1.0e-6)

    do node=1,node_count(hessian)
      abs_h = node_val(hessian, node)

      metric = ((1.0/sigma) * (alpha / eps)**power_C)**(2.0/n) * &
               (det(identity + (1.0/alpha) * abs_h) ** power_B) * &
               (identity + (1.0/alpha) * abs_h)

      call set(hessian, node, metric)
    end do

    call bound_metric(hessian, state)

    contains 

      function secant_method(func, x0, x1, eps) result(root)
        real, intent(in) :: x0, x1, eps
        real :: root
        interface
          function func(alpha) result(output)
            real, intent(in) :: alpha
            real :: output
          end function func
        end interface

        real :: fnminus, fn, fnplus
        real :: xnminus, xn, xnplus

        xnminus = x0
        xn = x1
        fnminus = func(x0)
        fn = func(x1)

        if (abs(fnminus) < eps) then
          root = xnminus
          return
        end if
        if (abs(fn) < eps) then
          root = xn
          return
        end if

        do while(.true.)
          xnplus = xn - ( (xn - xnminus) / (fn - fnminus) ) * fn
          fnplus = func(xnplus)
          if (abs(fnplus) < eps) then
            root = xnplus
            return
          end if

          xnminus = xn; fnminus = fn
          xn = xnplus; fn = fnplus
        end do
      end function secant_method

      function absolutify(hessian) result(abs_h)
        real, dimension(:, :), intent(in) :: hessian
        real, dimension(size(hessian, 1), size(hessian, 2)) :: abs_h, evecs
        real, dimension(size(hessian, 1)) :: evals
        
        call eigendecomposition_symmetric(hessian, evecs, evals)
        evals = abs(evals)
        call eigenrecomposition(abs_h, evecs, evals)
      end function absolutify
  end subroutine form_huang_metric

  function beta(alpha) result(p)
    real, intent(in) :: alpha
    real :: p
    real, dimension(ele_ngi(m_hessian, 1)) :: detwei, integrand_s
    real, dimension(m_hessian%dim, m_hessian%dim, ele_ngi(m_hessian, 1)) :: integrand_t, id
    integer :: j, ele

    do j=1,ele_ngi(m_hessian, 1)
      id(:, :, j) = get_matrix_identity(m_hessian%dim)
    end do

    p = 0
    do ele=1,ele_count(m_hessian)
      integrand_t = (ele_val_at_quad(m_hessian, ele) / alpha) + id
      do j=1,ele_ngi(m_hessian, 1)
        integrand_s(j) = det(integrand_t(:, :, j)) ** (1.0/m_gamma)
      end do

      call transform_to_physical(m_coordinate, ele, detwei=detwei)
      p = p + dot_product(integrand_s, detwei)
    end do

    p = p - m_sigma

    write(0,*) "alpha: ", alpha, "; beta: ", p
  end function beta

end module huang_metric
