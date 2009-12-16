#include "fdebug.h"

module anisotropic_zz_module
! Formaggia, Perotto, Micheletti et al.
! You really need to read the paper to understand this.
! If you understand the algorithm, the code is trivial,
! but not vice versa!

  use fields
  use sparse_tools
  use adjacency_lists
  use vector_tools
  use transform_elements
  use fetools
  use tensors
  use conformity_measurement
  use metric_tools
  use quicksort
  use sparsity_patterns
  use diagnostic_fields
  use state_module
  use global_parameters
  use field_options
  use meshdiagnostics
  use vtk_interfaces
  use bounding_box_metric
  use state_module
  use merge_tensors
  use edge_length_module
  use form_metric_field
  use vtk_interfaces
  use halos
  use limit_metric_module
  
  implicit none

  public :: compute_anisotropic_zz_metric, compute_g_hat, form_anisotropic_zz_metric, get_jacobian_azz

  contains

  subroutine form_anisotropic_zz_metric(state, metric)
    type(state_type), dimension(:), intent(in) :: state
    type(tensor_field), intent(inout) :: metric

    type(state_type) :: azz_fields
    type(scalar_field) :: field_s
    type(tensor_field) :: tmp_metric
    type(vector_field), pointer :: positions
    integer :: i, j

    logical :: debug_metric
    integer, save :: adaptcnt = 0
    character(len=100) :: buf
    type(scalar_field) :: edgelen

    do i=1,size(state)
      do j=1,scalar_field_count(state(i))
        field_s = extract_scalar_field(state(i), j)
        if (aliased(field_s)) cycle
        if (have_adapt_opt(trim(field_s%option_path), "/adaptivity_options/anisotropic_zienkiewicz_zhu")) then
          call insert(azz_fields, field_s, trim(field_s%name))
        end if
      end do
    end do

    if (scalar_field_count(azz_fields) == 0) return

    call allocate(tmp_metric, metric%mesh, "TemporaryMetric")
    call zero(tmp_metric)

    positions => extract_vector_field(state(1), "Coordinate")
    debug_metric = have_option("/mesh_adaptivity/hr_adaptivity/debug/write_metric_stages")
    if (debug_metric) then
      call allocate(edgelen, metric%mesh, "EdgeLengths")
    end if

    do i=1,scalar_field_count(azz_fields)
      write(buf, '(i0)') i
      field_s = extract_scalar_field(azz_fields, i)
      call halo_update(field_s)
      call compute_anisotropic_zz_metric(field_s, positions, tmp_metric)
      if (debug_metric) then
        call get_edge_lengths(tmp_metric, edgelen)
        call vtk_write_fieldS(trim("azz_metric_unbounded") // trim(buf), adaptcnt, positions, positions%mesh, &
                              sfields=(/field_s, edgelen/), tfields=(/tmp_metric/))
      end if
      call bound_metric(tmp_metric, state(1))
      if (debug_metric) then
        call get_edge_lengths(tmp_metric, edgelen)
        call vtk_write_fieldS(trim("azz_metric_bounded") // trim(buf), adaptcnt, positions, positions%mesh, &
                              sfields=(/field_s, edgelen/), tfields=(/tmp_metric/))
      end if
      call merge_tensor_fields(metric, tmp_metric)
    end do

    call deallocate(azz_fields)
    call deallocate(tmp_metric)

    if (debug_metric) then
      call deallocate(edgelen)
    end if
    adaptcnt = adaptcnt + 1

    call halo_update(metric)

  end subroutine form_anisotropic_zz_metric

  subroutine compute_anisotropic_zz_metric(field, positions, metric, eta_estimate)
    type(scalar_field), intent(in) :: field
    type(vector_field), intent(inout) :: positions
    type(tensor_field), intent(inout) :: metric
    real, intent(out), optional :: eta_estimate

    type(mesh_type) :: pwc_mesh
    type(tensor_field) :: pwc_metric

    integer :: ele
    real :: tau
    character(len=OPTION_PATH_LEN) :: path
    real :: eta, ele_eta

    path = trim(complete_field_path(trim(field%option_path))) // "/adaptivity_options/anisotropic_zienkiewicz_zhu"
    call get_option(trim(path) // "/tau", tau)
    ewrite(2,*) "Using tau = ", tau
    assert(tau > 0.0)

    call add_nelist(positions%mesh)

    assert(field%mesh%shape%degree == 1)
    assert(associated(metric%val))
    
    call zero(metric)

    ! Build the P0 metric
    pwc_mesh = piecewise_constant_mesh(metric%mesh, "PiecewiseConstantMesh")
    call allocate(pwc_metric, pwc_mesh, "PiecewiseConstantMetric")

    eta = 0.0
    do ele=1,ele_count(positions)
      call anisotropic_zz_element_metric(field, positions, pwc_metric, tau, ele, ele_eta)
      eta = eta + ele_eta
    end do
    eta = sqrt(eta)
    if (present(eta_estimate)) then
      eta_estimate = eta
    end if

    ewrite(2,*) "Current mesh error estimate: ", eta

    ! Project to P1 here
    call project_p0_metric_p1(positions, pwc_metric, metric)

    call deallocate(pwc_metric)
    call deallocate(pwc_mesh)
  end subroutine compute_anisotropic_zz_metric

  subroutine anisotropic_zz_element_metric(field, positions, pwc_metric, tau, ele, eta)
    type(scalar_field), intent(in) :: field
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(inout) :: pwc_metric
    real, intent(in) :: tau
    integer, intent(in) :: ele
    real, intent(out) :: eta

    real, dimension(positions%dim) :: lambda_k, g_evals, sorted_g_evals, out_evals, out_edges 
    real, dimension(positions%dim, positions%dim) :: g_hat, m_k, out_k, rt_k, vt_k, g_evecs, &
                                                  & sorted_g_evecs, out_evecs, g, j_k
    real :: patch_volume, transformed_patch_volume
    integer, dimension(positions%dim) :: sorted_idx
    integer :: i, dim
    real, parameter :: omega = 1.0
    real :: gamma, ele_vol, weight
    real :: sum_eta
    real :: dim_error
    real, dimension(positions%dim) :: evec

    dim = positions%dim
    ! Compute volume of ideal element
    if (dim == 2) then
      gamma = 3.0*sqrt(3.0) / 4.0
      weight = 1.0/3.0 ! 1/h**2, where h is ideal element length
    else
      gamma = 8.0 / (9.0 * sqrt(3.0))
      weight = 1.0 / (2 * sqrt(2.0/3.0))**2
    end if

    ! (a) and (b): compute g_hat
    g_hat = compute_g_hat(field, positions, ele_shape(pwc_metric, ele), ele, patch_vol=patch_volume)
    j_k = get_jacobian_azz(positions, ele)
    call svd(j_k, rt_k, lambda_k, vt_k)
    m_k = matmul(matmul(rt_k, get_mat_diag(1.0/(lambda_k**2))), transpose(rt_k))

    ! (e): compute transformed patch volume
    transformed_patch_volume = patch_volume / product(lambda_k)
    assert(transformed_patch_volume > gamma)

    ! (f): sorted eigenvalue decomposition of g_hat, in descending order
    call eigendecomposition_symmetric(g_hat, g_evecs, g_evals)
    call qsort(g_evals, sorted_idx)

    dim = positions%dim
    forall(i=1:dim)
      sorted_g_evals(i) = g_evals(sorted_idx(dim - i + 1))
      sorted_g_evecs(:, i) = g_evecs(:, sorted_idx(dim - i + 1))
    end forall

    ! Compute error estimator
    sum_eta = 0.0
    g = g_hat * patch_volume
    ele_vol = simplex_volume(positions, ele)
    do i=1,dim
      evec = rt_k(:, i)
      dim_error = lambda_k(i)**2 * dot_product(matmul(evec, g), evec)
      sum_eta = sum_eta + dim_error
    end do
    eta = (sum_eta * (1.0/product(lambda_k))**(2.0/dim))
    assert(.not. is_nan(eta))

    ! step (g): compute r
    out_evecs = compute_out_evecs(sorted_g_evecs)

    ! and compute lambda
    out_edges = compute_edge_lengths(tau, ele_count(positions), transformed_patch_volume, sorted_g_evals)
    assert(all(out_edges > 0.0))
    out_evals = eigenvalue_from_edge_length(out_edges)

    call eigenrecomposition(out_k, out_evecs, out_evals)
   
    out_k = weight * (omega * out_k + (1.0-omega)*m_k)

    call set(pwc_metric, ele, out_k)
  end subroutine anisotropic_zz_element_metric

  function compute_edge_lengths(tau, eles, vol, g_evals) result(edges)
    real, intent(in) :: tau
    integer, intent(in) :: eles
    real, intent(in) :: vol
    real, dimension(:), intent(in) :: g_evals
    real, dimension(size(g_evals)) :: edges, modified_g_evals, min_g
    real :: factor_a, factor_b, factor_c
    real :: max_length

    integer :: dim, i, j

    dim = size(g_evals)

    ! FIXME: make this more general by computing
    ! max_length as the diameter of the domain.
    assert(bounding_box_initialised)
    max_length = maxval(domain_width)
    factor_a = ((tau**2) / (dim*eles*vol))**(1.0/dim)
    assert(factor_a > 0.0)
    min_g = (max_length / factor_a)**(-1.0 * dim)
    modified_g_evals = max(g_evals, min_g)

    do i=1,dim
      j = dim - i + 1
      ! take a deep breath:
      factor_b = 1.0/sqrt(modified_g_evals(j))
      factor_c = (product(modified_g_evals))**((dim - 2.0) / (2.0*dim*dim))
      if (dim == 2) then
        assert(factor_c == 1.0)
      end if

      assert(factor_b > 0.0)
      assert(factor_c > 0.0)
      edges(i) = factor_a * factor_b * factor_c
      if (is_nan(edges(i)) .or. edges(i) > max_length) then
        edges(i) = max_length
      end if
    end do
  end function compute_edge_lengths

  function compute_out_evecs(g_evecs) result(out_evecs)
    real, dimension(:, :), intent(in) :: g_evecs
    real, dimension(size(g_evecs, 1), size(g_evecs, 1)) :: out_evecs
    integer :: dim
    integer :: i, j

    dim = size(g_evecs, 1)
    do i=1,dim
      j = dim - i + 1
      out_evecs(:, i) = g_evecs(:, j)
    end do
  end function compute_out_evecs

  function compute_g_hat(field, positions, pwc_shape, ele, patch_vol) result(g_hat)
    type(scalar_field), intent(in) :: field
    type(vector_field), intent(in) :: positions
    type(element_type), intent(in) :: pwc_shape
    integer, intent(in) :: ele
    real, optional, intent(out) :: patch_vol

    real, dimension(positions%dim) :: recovered_gradient
    real, dimension(positions%dim, ele_ngi(field, ele)) :: error_at_quad
    real, dimension(positions%dim, positions%dim) :: g, g_hat, g_int
    real, dimension(positions%dim, ele_loc(field, ele)) :: r


    integer :: i, neighbour
    real :: patch_volume
    integer :: gi, idim, jdim
    integer :: dim, loc, ngi
    logical, dimension(ele_count(positions)) :: seen_elements
    real, dimension(:, :), allocatable :: ele_gradient, detwei
    real, dimension(:, :, :, :), allocatable :: dm_t
    integer :: patch_size
    integer, dimension(:), allocatable :: patch_elements
    integer, dimension(:), pointer :: neighbours, nodes
    integer :: j

    ! Oh! For the want of a fecking set union in Fortran.

    ! Let's compute the patch around.

    seen_elements = .false.
    nodes => ele_nodes(positions, ele)
    assert(associated(positions%mesh%adj_lists))
    assert(associated(positions%mesh%adj_lists%nelist))
    do i=1,ele_loc(positions, ele)
      neighbours => row_m_ptr(positions%mesh%adj_lists%nelist, nodes(i))
      do j=1,size(neighbours)
        seen_elements(neighbours(j)) = .true.
      end do
    end do

    patch_size = count(seen_elements)
    allocate(patch_elements(patch_size))
    j = 1
    do i=1,size(seen_elements)
      if (seen_elements(i)) then
        patch_elements(j) = i
        j = j + 1
        if (j == patch_size + 1) then
          exit
        end if
      end if
    end do

    dim = positions%dim
    loc = ele_loc(field, ele)
    ngi = ele_ngi(field, ele)
    
    allocate(detwei(patch_size, ngi))
    allocate(dm_t(patch_size, loc, ngi, dim))
    allocate(ele_gradient(patch_size, dim))

    ! and compute their shape function derivatives etc.
    ! Note: this is grossly inefficient, but I don't see how to do it any better
    ! yet.
    do i=1,size(patch_elements)
      neighbour = patch_elements(i)
      assert(neighbour <= ele_count(field))
      call transform_to_physical(positions, neighbour, &
                                 ele_shape(field, neighbour), &
                                 detwei=detwei(i, :), dshape=dm_t(i, :, :, :))
    end do

    ! Step (a): compute recovered gradient 
    recovered_gradient = 0.0
    patch_volume = 0.0
    do i=1,size(patch_elements)
      neighbour = patch_elements(i)
      patch_volume = patch_volume + sum(detwei(i, :))

      r = dshape_rhs(dm_t(i, :, :, :), detwei(i, :))
      ele_gradient(i, :) = matmul(r, ele_val(field, neighbour)) / sum(detwei(i, :))

      ! and now do a weighted average to recover the approximate gradient:
      recovered_gradient = recovered_gradient + sum(detwei(i, :)) * ele_gradient(i,:)
    end do
    recovered_gradient = recovered_gradient / patch_volume

    ! Step (b): compute e_k and g_k and g_hat_k
    g = 0.0
    do i=1,size(patch_elements)
      neighbour = patch_elements(i)
      error_at_quad = spread(recovered_gradient, 2, size(error_at_quad, 2)) - &
                      matmul(reshape(ele_gradient(i, :), (/dim, 1/)), pwc_shape%n)

      do gi=1,ngi
        forall(idim=1:dim,jdim=1:dim)
          g_int(idim, jdim) = error_at_quad(idim, gi) * error_at_quad(jdim, gi)
        end forall
        g = g + g_int * detwei(i, gi)
      end do
    end do
    g_hat = g / patch_volume

    if (present(patch_vol)) then
      patch_vol = patch_volume
    end if
  end function compute_g_hat

  function get_jacobian_azz(positions, ele) result(J)
  !! This returns the jacobian of the affine map mapping from
  !! the ideal element used in these calculations, which
  !! is the element inscribed on the unit circle/sphere.
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: ele
    real, dimension(mesh_dim(positions), mesh_dim(positions)) :: J, mapped_pos, ideal_pos, ideal_pos_inv
    integer :: dim
    real, dimension(mesh_dim(positions), ele_loc(positions, ele)) :: val
    real, dimension(mesh_dim(positions)) :: t_k
    integer :: i

    real :: a, b, c

    val = ele_val(positions, ele)
    dim = mesh_dim(positions)
    if (dim == 2) then
      J(1, 1) = (val(1, 2) - val(1, 1))/sqrt(3.0)
      J(1, 2) = (2*val(1, 3) - val(1, 1) - val(1, 2))/3.0
      J(2, 1) = (val(2, 2) - val(2, 1))/sqrt(3.0)
      J(2, 2) = (2*val(2, 3) - val(2, 1) - val(2, 2))/3.0
    else
      ! We are trying to compute the Jacobian mapping
      ! the ideal element to the element we have, see?
      ! So, something like:
      ! (physical coords) = J * (ideal coords) + translation
      ! here, we write down the locations of 3 points of the tet
      ! in ideal coordinates:
      a = sqrt(2.0/3.0); b = sqrt(2.0)/3.0; c = 1.0/3.0
      ideal_pos(:, 1) = (/-a, -b, -c/)
      ideal_pos(:, 2) = (/a, -b, -c/)
      ideal_pos(:, 3) = (/0.0, 2.0*b, -c/)

      ! And its inverse, symbolically calculated thanks to sympy:
      ideal_pos_inv(:, 1) = (/-0.5/a, 0.5/a, 0.0/)
      ideal_pos_inv(:, 2) = (/(-1.0/6.0)/b, (-1.0/6.0)/b, (1.0/3.0)/b/)
      ideal_pos_inv(:, 3) = (/(-1.0/3.0)/c, (-1.0/3.0)/c, (-1.0/3.0)/c/)

      ! We compute the translation taking the centre to the origin:
      do i=1,dim
        t_k(i) = sum(val(i, :)) / 4.0
      end do

      ! And then perform the inverse translation:
      do i=1,dim
        mapped_pos(:, i) = val(:, i) - t_k
      end do

      ! And now we have the eqn:
      ! J * (ideal_coords) = (physical coords) - translation
      ! so, to solve this, we multiply on the right by
      ! the inverse of the ideal coordinates matrix:
      J = matmul(mapped_pos, ideal_pos_inv)
    end if
  end function get_jacobian_azz

end module anisotropic_zz_module
