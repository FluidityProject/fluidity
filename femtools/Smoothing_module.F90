#include "fdebug.h"

module smoothing_module
  use fields
  use sparse_tools
  use sparsity_patterns
  use solvers
  use global_parameters
  implicit none

  contains

  subroutine smooth_scalar(field_in,positions,field_out,alpha, path)
    implicit none
    !smoothing length tensor
    real, dimension(:,:), intent(in) :: alpha
    !input field
    type(scalar_field), intent(inout) :: field_in
    !coordinates field
    type(vector_field), intent(in) :: positions
    !output field, should have same mesh as input field
    type(scalar_field), intent(inout) :: field_out
    character(len=OPTION_PATH_LEN), intent(in) :: path
    
    !local variables
    type(csr_matrix) :: M
    type(csr_sparsity) :: M_sparsity
    type(scalar_field) :: RHSFIELD
    integer :: ele

    !allocate smoothing matrix
    M_sparsity=make_sparsity(field_in%mesh, &
         & field_in%mesh, name='HelmholtzSparsity')
    call allocate(M, M_sparsity, name="HelmholtzSmoothingMatrix")    
    call zero(M)
    
    !allocate RHSFIELD
    call allocate(rhsfield, field_in%mesh, "HelmholtzSmoothingRHS")
    call zero(rhsfield)

    ! Assemble M element by element.
    do ele=1, element_count(field_in)
       call assemble_smooth_scalar(M, rhsfield, positions, field_in, alpha, ele)
    end do

    call zero(field_out)
    call petsc_solve(field_out, M, rhsfield, option_path=trim(path))

    call deallocate(rhsfield)
    call deallocate(M)

  end subroutine smooth_scalar

  subroutine assemble_smooth_scalar(M, rhsfield, positions, field_in, alpha, ele)
    type(csr_matrix), intent(inout) :: M
    type(scalar_field), intent(inout) :: RHSFIELD
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: field_in
    real, dimension(:,:), intent(in) :: alpha
    integer, intent(in) :: ele

    ! value of field_in at quad points
    real, dimension(ele_ngi(positions,ele)) :: field_in_quad
    ! Locations of quadrature points
    real, dimension(positions%dim,ele_ngi(positions,ele)) :: X_quad
    ! smoothing tensor at quadrature points real,
    real,dimension(positions%dim,positions%dim,ele_ngi(positions,ele)) &
         &:: alpha_quad
    ! Derivatives of shape function:
    real, dimension(ele_loc(field_in,ele), &
         ele_ngi(field_in,ele), positions%dim) :: dshape_field_in
    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,ele)) :: detwei    
    ! Node numbers of field_in element.
    integer, dimension(:), pointer :: ele_field_in
    ! Shape functions.
    type(element_type), pointer :: shape_field_in
    ! Local Helmholtz matrix 
    real, dimension(ele_loc(field_in, ele), ele_loc(field_in, ele)) &
         & :: field_in_mat
    ! Local right hand side.
    real, dimension(ele_loc(field_in, ele)) :: lrhsfield
    integer :: i

    ele_field_in=>ele_nodes(field_in, ele)
    shape_field_in=>ele_shape(field_in, ele)

    ! Locations of quadrature points.
    X_quad=ele_val_at_quad(positions, ele)

    !value of tensor at quads
    forall(i=1:ele_ngi(positions,ele))
       alpha_quad(:,:,i) = alpha
    end forall

    field_in_quad = ele_val_at_quad(field_in, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, shape_field_in, dshape&
         &=dshape_field_in, detwei=detwei)

    ! Local assembly:
    field_in_mat=dshape_tensor_dshape(dshape_field_in, alpha_quad, &
         dshape_field_in, detwei) + shape_shape(shape_field_in&
         &,shape_field_in, detwei)

    lrhsfield=shape_rhs(shape_field_in, field_in_quad*detwei)

    ! Global assembly:
    call addto(M, ele_field_in, ele_field_in, field_in_mat)

    call addto(rhsfield, ele_field_in, lrhsfield)

  end subroutine assemble_smooth_scalar
end module smoothing_module

