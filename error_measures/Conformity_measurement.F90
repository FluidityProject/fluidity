#include "fdebug.h"

module conformity_measurement
!! See 
!! A universal measure of the conformity of a mesh with respect to an anisotropic metric field
!! P. Labbe and J. Dompierre and M.-G. Vallet and F. Guibault and J.-Y. Trépanier
!! International Journal for Numerical Methods in Engineering
!! 2004, vol 61, issue 15

  use vector_tools
  use transform_elements
  use sparse_tools
  use unittest_tools
  use fetools
  use metric_tools
  use fields
  use state_module
  use fefields
  use meshdiagnostics
  use merge_tensors
  use limit_metric_module

  implicit none

  private
  public :: elemental_metric, compute_mesh_conformity, insert_mesh_conformity,  &
          & compute_mesh_metric, metric_ratio_bounds, project_p0_metric_p1

  contains

  function elemental_metric(metric, positions, ele) result(m)
    type(tensor_field), intent(in) :: metric
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    real, dimension(metric%dim(1), metric%dim(2)) :: m
    real, dimension(metric%dim(1), metric%dim(2), ele_ngi(metric, ele)) :: A
    real, dimension(ele_ngi(metric, ele)) :: detwei
    type(element_type), pointer :: t_shape
    integer :: i, j, dim

    t_shape => ele_shape(metric, ele)

    call transform_to_physical(positions, ele, detwei=detwei)

    A = ele_val_at_quad(metric, ele)

    dim = positions%dim
    do i=1,dim
      do j=1,dim
        m(i, j) = dot_product(A(i, j, :), detwei)
      end do
    end do

    ! sum(detwei) is the volume.

    m = m / sum(detwei)

  end function elemental_metric

  subroutine compute_mesh_conformity(metric, positions, conformity)
    type(tensor_field), intent(in) :: metric
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(inout) :: conformity

    integer :: ele
    real, dimension(positions%dim, positions%dim) :: A, B, R, invA, invB, I

    ! conformity is a piecewise constant scalar field:
    assert(ele_loc(conformity, 1) == 1)

    I = get_matrix_identity(positions%dim)

    do ele=1,ele_count(positions)
      A = simplex_tensor(positions, ele)
      B = elemental_metric(metric, positions, ele)
      invA = A; call invert(A)
      invB = B; call invert(B)

      R = matmul(invA, B) + matmul(invB, A) - 2 * I
      call set(conformity, ele, frob(R))
    end do
  end subroutine compute_mesh_conformity

!  subroutine compute_mesh_metric(positions, mesh_metric)
!    type(vector_field), intent(in) :: positions
!    type(tensor_field), intent(inout) :: mesh_metric
!
!    real, dimension(positions%dim, positions%dim) :: A
!    real, dimension(positions%dim, positions%dim, ele_ngi(positions, 1)) :: A_tmp
!    real, dimension(ele_ngi(positions, 1)) :: detwei
!    integer :: ele, node, gi
!    type(element_type), pointer :: shape
!    type(scalar_field) :: lumped_mass
!    real, dimension(ele_loc(positions, 1), ele_loc(positions, 1)) :: mass_matrix
!
!    call zero(mesh_metric)
!    call allocate(lumped_mass, mesh_metric%mesh, "LumpedMass")
!    call zero(lumped_mass)
!    shape => ele_shape(positions, 1)
!
!    do ele=1,ele_count(positions)
!      call transform_to_physical(ele_val(positions, ele), shape, detwei=detwei)
!      A = simplex_tensor(positions, ele)
!      do gi=1,size(A_tmp, 3)
!        A_tmp(:, :, gi) = A
!      end do
!
!      call addto(mesh_metric, ele_nodes(positions, ele), shape_tensor_rhs(shape, A_tmp, detwei))
!
!      mass_matrix = shape_shape(shape, shape, detwei)
!      call addto(lumped_mass, ele_nodes(positions, ele), sum(mass_matrix, 2))
!    end do
!
!    do node=1,node_count(positions)
!      call set(mesh_metric, node, node_val(mesh_metric, node) / node_val(lumped_mass, node))
!    end do
!
!    call deallocate(lumped_mass)
!  end subroutine compute_mesh_metric

  function frob(R)
    real, dimension(:, :), intent(in) :: R
    real :: frob

    real, dimension(size(R, 1), size(R, 2)) :: X
    integer :: i
    real :: trace

    X = matmul(transpose(R), R)
    trace = 0.0
    do i=1,size(R, 1)
      trace = trace + X(i, i)
    end do
    frob = sqrt(trace)
  end function frob

  subroutine insert_mesh_conformity(state, error_metric)
    type(state_type), dimension(:), intent(inout) :: state
    type(tensor_field), intent(in), target, optional :: error_metric

    type(scalar_field) :: conformity
    type(vector_field), pointer :: positions
    type(tensor_field), pointer :: metric

    if (present(error_metric)) then
      metric => error_metric
    else
      metric => extract_tensor_field(state(1), "ErrorMetric")
    end if

    positions => extract_vector_field(state(1), "Coordinate")

    conformity = piecewise_constant_field(positions%mesh, "MeshConformity")
    call compute_mesh_conformity(metric, positions, conformity)
    call insert(state(1), conformity, "MeshConformity")
    call deallocate(conformity)
  end subroutine insert_mesh_conformity

  function metric_ratio_bounds(positions, metric) result(minmax)
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(in) :: metric
    
    real, dimension(2) :: minmax

    integer :: ele
    real, dimension(2) :: ele_minmax
    real, dimension(metric%dim(1), metric%dim(2)) :: A, B

    assert(ele_count(positions) > 0)
    A = edge_length_from_eigenvalue(simplex_tensor(positions, 1))
    B = edge_length_from_eigenvalue(elemental_metric(metric, positions, 1))
    minmax = anisotropic_ratio(A, B)
        
    do ele = 2, ele_count(positions)
      A = edge_length_from_eigenvalue(simplex_tensor(positions, ele))
      B = edge_length_from_eigenvalue(elemental_metric(metric, positions, ele))
      ele_minmax = anisotropic_ratio(A, B)
      minmax(1) = min(minmax(1), ele_minmax(1))
      minmax(2) = max(minmax(2), ele_minmax(2))
    end do

  end function metric_ratio_bounds

  function anisotropic_ratio(tensor1, tensor2) result(minmax)
    real, dimension(:, :), intent(in) :: tensor1, tensor2
    real, dimension(2) :: minmax  ! (/min, max/) of ratios of edge lengths
    real, dimension(size(tensor1, 1), size(tensor1, 1)) :: F, T, Finv, evecs, tensor3
    real, dimension(size(tensor1, 1)) :: evals
    integer :: dim

    dim = size(tensor1, 1)

    call eigendecomposition_symmetric(tensor2, evecs, evals)
    F = get_deformation_matrix(tensor2, evecs, evals)
    Finv = inverse(F)

    T = transpose(Finv)
    tensor3 = matmul(matmul(T, tensor1), transpose(T))
    call eigendecomposition_symmetric(tensor3, evecs, evals)

    minmax = (/minval(evals), maxval(evals)/)

  end function anisotropic_ratio

  subroutine compute_mesh_metric(mesh, metric)
    type(vector_field), intent(in) :: mesh
    type(tensor_field), intent(out) :: metric

    type(tensor_field) :: pwc_metric
    type(mesh_type) :: pwc_mesh
    integer :: ele

    call allocate(metric, mesh%mesh, "MeshSizingMetric")
    pwc_mesh = piecewise_constant_mesh(mesh%mesh, "PiecewiseConstantMesh")
    call allocate(pwc_metric, pwc_mesh, "PWCMeshSizingMetric")

    do ele=1,ele_count(mesh)
      call set(pwc_metric, ele, simplex_tensor(mesh, ele))
    end do

    call project_p0_metric_p1(mesh, pwc_metric, metric, target=ele_count(mesh))
    call deallocate(pwc_metric)
    call deallocate(pwc_mesh)
  end subroutine compute_mesh_metric

  subroutine project_p0_metric_p1(positions, pwc_metric, metric, target)
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(in) :: pwc_metric
    type(tensor_field), intent(inout) :: metric
    integer, intent(in), optional :: target

    integer :: node
    integer :: ele
    integer, dimension(:), pointer :: eles
    integer :: i
    real, dimension(positions%dim, positions%dim) :: tmp_metric
    integer :: dim
    real :: sum_vol, ele_vol
    integer :: pwc_xpct, metric_xpct
    real :: beta
    type(csr_sparsity), pointer :: nelist

    dim = mesh_dim(positions)
    nelist => extract_nelist(metric)

    do node=1,node_count(metric)
      tmp_metric = 0.0
      eles => row_m_ptr(nelist, node)
      sum_vol = 0.0
      do i=1,size(eles)
        ele = eles(i)
        ele_vol = simplex_volume(positions, ele)
        sum_vol = sum_vol + ele_vol
        tmp_metric = tmp_metric + ele_vol * node_val(pwc_metric, ele)
      end do
      tmp_metric = tmp_metric / sum_vol
      call set(metric, node, tmp_metric)
    end do

    if (present(target)) then
      pwc_xpct = target
    else
      pwc_xpct = expected_elements(positions, pwc_metric) * 1.25
    end if
    metric_xpct = expected_elements(positions, metric)
    beta = ((1.0 / metric_xpct) * pwc_xpct) ** (2.0 / dim)
    call scale(metric, beta)

  end subroutine project_p0_metric_p1


end module conformity_measurement
