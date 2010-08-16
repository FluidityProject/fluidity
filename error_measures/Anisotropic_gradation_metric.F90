#include "fdebug.h"

module anisotropic_gradation

  use fields
  use spud
  use initialise_fields_module
  use adjacency_lists
  use sparse_tools
  use linked_lists
  use gradation_metric
  use merge_tensors
  use vector_tools
  use metric_tools
  use state_module
  use form_metric_field

  implicit none

  public :: initialise_anisotropic_gradation
  public :: form_anisotropic_gradation_metric

  logical :: use_anisotropic_gradation = .true.

  contains

  subroutine initialise_anisotropic_gradation
    use_anisotropic_gradation = have_option("/mesh_adaptivity/hr_adaptivity/anisotropic_gradation")
  end subroutine

  subroutine form_anisotropic_gradation_metric(metric, positions, state, noits, gamma_field)
    type(tensor_field), intent(inout), target :: metric
    type(vector_field), intent(in)            :: positions
    type(state_type), intent(in) :: state
    integer, optional, intent(out)            :: noits
    type(tensor_field), optional, intent(in)  :: gamma_field

    type(tensor_field) :: gamma
    type(csr_sparsity), pointer :: nnlist_sparsity
    type(csr_matrix)   :: nnlist
    type(elist)        :: edge_list

    integer :: count, global_its, end_marker, p, q
    logical :: changed_p, changed_q, is_constant
    real    :: dist

    type(mesh_type), pointer :: mesh
    character(len=*), parameter :: path = "/mesh_adaptivity/hr_adaptivity/anisotropic_gradation/tensor_field::Gamma"

    real, dimension(positions%dim, positions%dim) :: val_gamma, val_p, val_q, grad, const_gamma, prev_grad

    integer :: stat

    ewrite(1,*) "Using anisotropic gradation algorithm"
    mesh => metric%mesh

    if (present(gamma_field)) then
      gamma = gamma_field
      is_constant = (gamma%field_type == FIELD_TYPE_CONSTANT)
    else
      is_constant = (have_option(path // "/anisotropic_symmetric/constant"))
      if (is_constant) then
        call allocate(gamma, mesh, "Gamma", field_type=FIELD_TYPE_CONSTANT)
      else
        call allocate(gamma, mesh, "Gamma")
      end if
      call initialise_field(gamma, path, positions)
    end if

    const_gamma = node_val(gamma, 1)

    ! The same as the old gradation algorithm. We store whether a pair of nodes
    ! is in the edge_list by marking nnlist(i, j) < 0 if it is not in the list
    ! and > 0 if it is.

    nnlist_sparsity => extract_nnlist(mesh)
    call allocate(nnlist, sparsity=nnlist_sparsity, type=CSR_INTEGER, name="NodeNodeList")
    nnlist%ival = -1
    call construct_edge_list(mesh, nnlist, edge_list)

    end_marker = edge_list%length
    global_its = 0

    do while (edge_list%length /= 0)

      end_marker = end_marker - 1
      if (end_marker == 0) then
        global_its = global_its + 1
        end_marker = edge_list%length
      end if

      changed_p = .false.
      changed_q = .false.

      call wrap_pop(nnlist, edge_list, p, q, count)

      dist = distance(positions, p, q)

      val_p = edge_lengths_from_metric(node_val(metric, p))
      val_q = edge_lengths_from_metric(node_val(metric, q))

      if (is_constant) then
        val_gamma = const_gamma
      else
        val_gamma = (node_val(gamma, p) + node_val(gamma, q)) / 2.0
      end if
      prev_grad = (val_q - val_p) / dist

      grad = anisotropic_min(prev_grad, val_gamma)
      changed_q = (prev_grad .fne. grad)
      if (changed_q) then
        val_q = dist * grad + val_p
        call set(metric, q, metric_from_edge_lengths(val_q))
      end if

      if (is_constant) then
        val_gamma = const_gamma
      else
        val_gamma = (node_val(gamma, p) + node_val(gamma, q)) / 2.0
      end if

      prev_grad = (val_p - val_q) / dist
      grad = anisotropic_min(prev_grad, val_gamma)
      changed_p = (prev_grad .fne. grad)
      if (changed_p) then
        val_p = dist * grad + val_q
        call set(metric, p, metric_from_edge_lengths(val_p))
      end if

      if (changed_p) then
        call tag_edges(nnlist, edge_list, p, q, count)
      end if
      if (changed_q) then
        call tag_edges(nnlist, edge_list, q, p, count)
      end if
    end do

    if (present(noits)) then
      noits = global_its
    end if

    if (.not. present(gamma_field)) then
      call deallocate(gamma)
    end if
    call deallocate(nnlist)
    
    call bound_metric(metric, state, stat=stat)
  end subroutine

  function anisotropic_min(tensor1, tensor2) result(tensor3)
    real, dimension(:, :), intent(in) :: tensor1, tensor2
    real, dimension(size(tensor1, 1), size(tensor1, 1)) :: tensor3, F, T, Finv, evecs
    real, dimension(size(tensor1, 1)) :: evals, ones
    integer :: i, dim

    dim = size(tensor1, 1)
    ones = 1

    if (all(tensor2 == 0.0)) then
      call eigendecomposition_symmetric(tensor1, evecs, evals)
      do i=1,dim
        evals(i) = min(evals(i), 0.0)
      end do
      call eigenrecomposition(tensor3, evecs, evals)
      return
    end if

    ! So we are dealing with the non-degenerate case.

    call eigendecomposition_symmetric(tensor2, evecs, evals)
    F = get_deformation_matrix(tensor2, evecs, evals)
    Finv = inverse(F)

    T = transpose(Finv)
    tensor3 = matmul(matmul(T, tensor1), transpose(T))
    call eigendecomposition_symmetric(tensor3, evecs, evals)

    evals = min(evals, ones)
    call eigenrecomposition(tensor3, evecs, evals)
    T = F
    tensor3 = matmul(matmul(transpose(T), tensor3), T)
  end function


end module anisotropic_gradation
