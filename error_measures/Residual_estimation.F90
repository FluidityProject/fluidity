#include "fdebug.h"

module residual_estimation

  use fields
  use field_derivatives
  use state_module
  use global_parameters
  implicit none

  contains

  subroutine advection_diffusion_residual(field, old_field, positions, nl_velocity, grid_velocity, diffusivity, phase, residual)
    type(scalar_field), intent(in) :: field
    type(scalar_field), intent(in) :: old_field
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: nl_velocity
    type(vector_field), intent(in) :: grid_velocity
    type(tensor_field), intent(in) :: diffusivity
    integer, intent(in) :: phase
    type(scalar_field), intent(inout) :: residual

    type(mesh_type), pointer :: mesh

    type(vector_field) :: grad_field
    type(vector_field) :: grad_old_field
    type(scalar_field) :: divergence

    integer :: node, j
    real, dimension(mesh_dim(residual)) :: tmp

    mesh => residual%mesh
    call allocate(divergence, mesh, "Divergence")
    call allocate(grad_field, mesh_dim(mesh), mesh, "Gradient of field")
    call allocate(grad_old_field, mesh_dim(mesh), mesh, "Gradient of old field")
    call zero(residual)

    call grad(field, positions, grad_field)
    call grad(old_field, positions, grad_old_field)

    do node=1,node_count(mesh)
      ! Time discretisation
      call addto(residual, node, (node_val(field, node) - node_val(old_field, node)) / dt)
      call addto(residual, node, theta(phase) * dot_product(node_val(nl_velocity, node), node_val(grad_field, node)))
      call addto(residual, node, (1 - theta(phase)) * dot_product(node_val(nl_velocity, node), node_val(grad_old_field, node)))

      ! Advective term
      call addto(residual, node, dot_product(node_val(nl_velocity, node) - node_val(grid_velocity, node), node_val(grad_field, node)))
    end do
    call deallocate(grad_old_field)

    ! Now we multiply grad_field by the diffusivity tensor, in-place, to save on memory.
    do node=1,node_count(mesh)
      tmp = matmul(node_val(diffusivity, node), node_val(grad_field, node))
      do j=1,size(tmp)
        grad_field%val(j,node) = tmp(j)
      end do
    end do

    ! Now grad_field is diffusivity * grad(T).
    ! Since we want the divergence:

    call div(grad_field, positions, divergence)

    ! Now for the diffusive term.
    do node=1,node_count(mesh)
      ! Diffusive term
      call addto(residual, node, -1 * node_val(divergence, node))
    end do

    call deallocate(divergence)
    call deallocate(grad_field)
  end subroutine advection_diffusion_residual

  subroutine interpolation_residual(field, positions, hessian, residual)
    type(scalar_field), intent(in) :: field
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(in) :: hessian
    type(scalar_field), intent(inout) :: residual

    type(mesh_type), pointer :: mesh
    type(patch_type) :: patch
    integer :: node, nnode, i
    real :: err
    real, dimension(hessian%dim, hessian%dim) :: avg_hessian, evecs
    real, dimension(hessian%dim) :: edge, evals

    mesh => field%mesh
    call zero(residual)

    do node=1,node_count(mesh)
      patch = get_patch_node(mesh, node, level=1)
      err = 0.0

      do i=1,patch%count
        nnode = patch%elements(i)
        edge = node_val(positions, node) - node_val(positions, nnode)
        avg_hessian = (node_val(hessian, node) + node_val(hessian, nnode)) / 2.0
        call eigendecomposition_symmetric(avg_hessian, evecs, evals)
        evals = abs(evals)
        call eigenrecomposition(avg_hessian, evecs, evals)
        err = max(err, dot_product(edge, matmul(avg_hessian, edge)))
      end do
      
      call addto(residual, node, err)
      deallocate(patch%elements)
    end do
  end subroutine interpolation_residual

end module residual_estimation
