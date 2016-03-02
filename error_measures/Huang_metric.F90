#include "fdebug.h"

module huang_metric_module
! 10.1016/j.jcp.2004.10.024

  use fldebug
  use vector_tools
  use unittest_tools, only: is_nan
  use global_parameters, only: domain_volume, OPTION_PATH_LEN
  use spud
  use metric_tools
  use transform_elements
  use fields
  use state_module
  use vtk_interfaces
  use field_options
  use form_metric_field

  implicit none

  type(tensor_field), save :: m_hessian
  type(vector_field), pointer :: m_coordinate
  real :: m_gamma, m_sigma

  private
  public :: form_huang_metric

  contains

  subroutine form_huang_metric(hessian, field, coordinate, weight)
    type(tensor_field), intent(inout) :: hessian
    type(scalar_field), intent(in) :: field
    type(vector_field), intent(in), target :: coordinate
    real, intent(in) :: weight

    integer :: n, m, l, p, q
    integer, dimension(2) :: tmp
    real :: gamma, sigma, alpha, eps, power_A, power_B, power_C
    integer :: node
    real, dimension(hessian%dim(1), hessian%dim(2)) :: identity, abs_h, metric
    character(len=OPTION_PATH_LEN) :: path

    assert(field%mesh%shape%degree == 1) ! this can be generalised, but currently I don't have the time

    m_coordinate => coordinate

    n = mesh_dim(field)
    l = field%mesh%shape%degree + 1
    p = 1

    path = trim(complete_field_path(trim(field%option_path))) // "/adaptivity_options/huang_metric"
    call get_option(trim(path) // "/seminorm", tmp)
    m = tmp(1)
    q = tmp(2)

    eps = weight
    gamma = (real(n) / q) + (2 - m)
    power_A = 1 + real((n * (p - 1)))/(p * gamma) + max(0.0, real(n)/(p*gamma) - 1)
    assert(domain_volume > 0)
    sigma = (2 ** power_A) * domain_volume
    m_sigma = sigma

    identity = get_matrix_identity(n)
    power_B = (1.0/n) * ((2.0/gamma) - 1)
    power_C = real(n) / (l - m)

    do node=1,node_count(hessian)
      abs_h = absolutify(node_val(hessian, node))
      call set(hessian, node, abs_h)
    end do

    m_hessian = hessian
    m_gamma = gamma

    alpha = find_root(beta, 0.05, 2.5, 1.0e-6)

    do node=1,node_count(hessian)
      abs_h = node_val(hessian, node)

      metric = ((1.0/sigma) * ((alpha / eps)**power_C))**(2.0/n) * &
               (det(identity + (1.0/alpha) * abs_h) ** power_B) * &
               (identity + (1.0/alpha) * abs_h) * 4

      call set(hessian, node, metric)
    end do

    contains 

      function find_root(func, x0, x1, eps) result(root)
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

        if (fnminus * fn < 0) then
          do while(.true.)
            xnplus = (fn*xnminus - fnminus*xn) / (fn - fnminus)
            fnplus = func(xnplus)
            if (abs(fnplus) < eps) then
              root = xnplus
              return
            end if

            if (sign(1.0, fnplus) == sign(1.0, fnminus)) then
              fnminus = fnplus
              xnminus = xnplus
            else
              fn = fnplus
              xn = xnplus
            end if
          end do
        else
          do while(.true.)
            xnplus = xn - ( (xn - xnminus) / (fn - fnminus) ) * fn
            assert(.not. is_nan(xnplus))
            fnplus = func(xnplus)
            if (abs(fnplus) < eps) then
              root = xnplus
              return
            end if

            xnminus = xn; fnminus = fn
            xn = xnplus; fn = fnplus
          end do
        end if
      end function find_root

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
    real, dimension(m_hessian%dim(1), m_hessian%dim(2), ele_ngi(m_hessian, 1)) :: integrand_t, id
    integer :: j, ele

    do j=1,ele_ngi(m_hessian, 1)
      id(:, :, j) = get_matrix_identity(m_hessian%dim(1))
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

    write(0,*) "alpha: ", alpha, "; \int rho: ", p, "; sigma: ", m_sigma, "; output: ", p - m_sigma
    p = p - m_sigma

  end function beta

end module huang_metric_module
