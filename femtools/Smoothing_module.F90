#include "fdebug.h"

module smoothing_module
  use state_module
  use fields
  use sparse_tools
  use sparsity_patterns
  use solvers
  use global_parameters, only : OPTION_PATH_LEN
  implicit none

  private
  
  public :: smooth_scalar, smooth_vector
  public :: anisotropic_smooth_scalar, anisotropic_smooth_vector
  public :: anisotropic_smooth_tensor, length_scale_tensor

contains

  subroutine smooth_scalar(field_in,positions,field_out,alpha, path)

    !smoothing length tensor
    real, dimension(:,:), intent(in) :: alpha
    !input field
    type(scalar_field), intent(inout) :: field_in
    !coordinates field
    type(vector_field), intent(in) :: positions
    !output field, should have same mesh as input field
    type(scalar_field), intent(inout) :: field_out
    character(len=*), intent(in) :: path
    
    !local variables
    type(csr_matrix) :: M
    type(csr_sparsity) :: M_sparsity
    type(scalar_field) :: RHSFIELD
    integer :: ele

    !allocate smoothing matrix
    M_sparsity=make_sparsity(field_in%mesh, &
         & field_in%mesh, name='HelmholtzScalarSparsity')
    call allocate(M, M_sparsity, name="HelmholtzScalarSmoothingMatrix")   
    call deallocate(M_sparsity) 
    call zero(M)
    
    !allocate RHSFIELD
    call allocate(rhsfield, field_in%mesh, "HelmholtzScalarSmoothingRHS")
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

  subroutine smooth_vector(field_in,positions,field_out,alpha, path)

    !smoothing length tensor
    real, dimension(:,:), intent(in) :: alpha
    !input field
    type(vector_field), intent(inout) :: field_in
    !coordinates field
    type(vector_field), intent(in) :: positions
    !output field, should have same mesh as input field
    type(vector_field), intent(inout) :: field_out
    character(len=*), intent(in) :: path
    
    !local variables
    type(csr_matrix) :: M
    type(csr_sparsity) :: M_sparsity
    type(vector_field) :: RHSFIELD
    integer :: ele

    !allocate smoothing matrix
    M_sparsity=make_sparsity(field_in%mesh, &
         & field_in%mesh, name='HelmholtzVectorSparsity')
    call allocate(M, M_sparsity, name="HelmholtzVectorSmoothingMatrix")   
    call deallocate(M_sparsity) 
    call zero(M)
    
    !allocate RHSFIELD
    call allocate(rhsfield, field_in%dim, field_in%mesh, "HelmholtzVectorSmoothingRHS")
    call zero(rhsfield)

    ! Assemble M element by element.
    do ele=1, element_count(field_in)
       call assemble_smooth_vector(M, rhsfield, positions, field_in, alpha, ele)
    end do

    call zero(field_out)
    call petsc_solve(field_out, M, rhsfield, option_path=trim(path))

    call deallocate(rhsfield)
    call deallocate(M)

  end subroutine smooth_vector

  subroutine anisotropic_smooth_scalar(field_in,positions,field_out,alpha,path)

    !smoothing length tensor
    real, intent(in)                           :: alpha
    !input field
    type(scalar_field), intent(inout)    :: field_in
    !coordinates field
    type(vector_field), intent(in)             :: positions
    !output field, should have same mesh as input field
    type(scalar_field), intent(inout) :: field_out
    character(len=*), intent(in)               :: path
    
    !local variables
    type(csr_matrix) :: M
    type(csr_sparsity) :: M_sparsity
    type(scalar_field) :: RHSFIELD
    integer :: ele

    !allocate smoothing matrix
    M_sparsity=make_sparsity(field_in%mesh, &
         & field_in%mesh, name='HelmholtzScalarSparsity')
    call allocate(M, M_sparsity, name="HelmholtzScalarSmoothingMatrix")   
    call deallocate(M_sparsity) 
    call zero(M)
    
    !allocate RHSFIELD
    call allocate(rhsfield, field_in%mesh, "HelmholtzScalarSmoothingRHS")
    call zero(rhsfield)

    ! Assemble M element by element.
    do ele=1, element_count(field_in)
       call assemble_anisotropic_smooth_scalar(M, rhsfield, positions, field_in, alpha, ele)
    end do

    call zero(field_out)
    call petsc_solve(field_out, M, rhsfield, option_path=trim(path))

    call deallocate(rhsfield)
    call deallocate(M)

  end subroutine anisotropic_smooth_scalar

  subroutine anisotropic_smooth_vector(field_in,positions,field_out,alpha,path)

    !smoothing length tensor
    real, intent(in) :: alpha
    !input field
    type(vector_field), intent(inout) :: field_in
    !coordinates field
    type(vector_field), intent(in) :: positions
    !output field, should have same mesh as input field
    type(vector_field), intent(inout) :: field_out
    character(len=*), intent(in) :: path
    
    !local variables
    type(csr_matrix) :: M
    type(csr_sparsity) :: M_sparsity
    type(vector_field) :: RHSFIELD
    integer :: ele

    !allocate smoothing matrix
    M_sparsity=make_sparsity(field_in%mesh, &
         & field_in%mesh, name='HelmholtzVectorSparsity')
    call allocate(M, M_sparsity, name="HelmholtzVectorSmoothingMatrix")   
    call deallocate(M_sparsity) 
    call zero(M)
    
    !allocate RHSFIELD
    call allocate(rhsfield, field_in%dim, field_in%mesh, "HelmholtzVectorSmoothingRHS")
    call zero(rhsfield)

    ! Assemble M element by element.
    do ele=1, element_count(field_in)
       call assemble_anisotropic_smooth_vector(M, rhsfield, positions, field_in, alpha, ele)
    end do

    call zero(field_out)
    call petsc_solve(field_out, M, rhsfield, option_path=trim(path))

    call deallocate(rhsfield)
    call deallocate(M)

  end subroutine anisotropic_smooth_vector

  subroutine anisotropic_smooth_tensor(field_in,positions,field_out,alpha, path)

    !smoothing length tensor
    real, intent(in) :: alpha
    !input field
    type(tensor_field), intent(inout) :: field_in
    !coordinates field
    type(vector_field), intent(in) :: positions
    !output field, should have same mesh as input field
    type(tensor_field), intent(inout) :: field_out
    character(len=*), intent(in) :: path
    
    !local variables
    type(csr_matrix) :: M
    type(csr_sparsity) :: M_sparsity
    type(tensor_field) :: RHSFIELD
    integer :: ele

    !allocate smoothing matrix
    M_sparsity=make_sparsity(field_in%mesh, &
         & field_in%mesh, name='HelmholtzTensorSparsity')
    call allocate(M, M_sparsity, name="HelmholtzTensorSmoothingMatrix")   
    call deallocate(M_sparsity) 
    call zero(M)
    
    !allocate RHSFIELD
    call allocate(rhsfield, field_in%mesh, "HelmholtzTensorSmoothingRHS")
    call zero(rhsfield)

    ! Assemble M element by element.
    do ele=1, element_count(field_in)
       call assemble_anisotropic_smooth_tensor(M, rhsfield, positions, field_in, alpha, ele)
    end do

    call zero(field_out)
    call petsc_solve(field_out, M, rhsfield, option_path=trim(path))

    call deallocate(rhsfield)
    call deallocate(M)

  end subroutine anisotropic_smooth_tensor

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

  subroutine assemble_smooth_vector(M, rhsfield, positions, field_in, alpha, ele)

    type(csr_matrix), intent(inout) :: M
    type(vector_field), intent(inout) :: RHSFIELD
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: field_in
    real, dimension(:,:), intent(in) :: alpha
    integer, intent(in) :: ele

    ! value of field_in at quad points
    real, dimension(positions%dim,ele_ngi(positions,ele)) :: field_in_quad
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
    real, dimension(positions%dim, ele_loc(field_in, ele)) :: lrhsfield
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

    ! Local assembly:    assert(field_in%mesh==field_out%mesh))
    field_in_mat=dshape_tensor_dshape(dshape_field_in, alpha_quad, &
         dshape_field_in, detwei) + shape_shape(shape_field_in&
         &,shape_field_in, detwei)

    lrhsfield=shape_vector_rhs(shape_field_in, field_in_quad, detwei)

    ! Global assembly:
    call addto(M, ele_field_in, ele_field_in, field_in_mat)

    call addto(rhsfield, ele_field_in, lrhsfield)

  end subroutine assemble_smooth_vector

  subroutine assemble_anisotropic_smooth_scalar(M, rhsfield, positions, field_in, alpha, ele)
    type(csr_matrix), intent(inout) :: M
    type(scalar_field), intent(inout) :: RHSFIELD
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: field_in
    real, intent(in) :: alpha
    integer, intent(in) :: ele
    real :: beta
    real, dimension(ele_ngi(positions,ele))                                      :: field_in_quad
    real, dimension(positions%dim,ele_ngi(positions,ele))                        :: X_quad
    real,dimension(positions%dim,positions%dim,ele_ngi(positions,ele))           :: mesh_tensor_quad
    real, dimension(ele_loc(field_in,ele), ele_ngi(field_in,ele), positions%dim) :: dshape_field_in
    real, dimension(ele_ngi(positions,ele))                                      :: detwei
    integer, dimension(:), pointer                                               :: ele_field_in
    type(element_type), pointer                                                  :: shape_field_in
    real, dimension(ele_loc(field_in, ele), ele_loc(field_in, ele))              :: field_in_mat
    real, dimension(ele_loc(field_in, ele))                                      :: lrhsfield

    ele_field_in=>ele_nodes(field_in, ele)
    shape_field_in=>ele_shape(field_in, ele)

    ! Locations of quadrature points. Fields evaluated at quads.
    X_quad=ele_val_at_quad(positions, ele)
    field_in_quad = ele_val_at_quad(field_in, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, shape_field_in, dshape&
         &=dshape_field_in, detwei=detwei)

    ! mesh size tensor = gamma * (mesh size)**2
    ! gamma is a function of mesh shape e.g. hexahedral, unstructured triangular etc.
    ! Helmholtz smoothing factor = alpha * 1/24 * mesh size tensor
    ! Factor of 1/24 comes from Geurts & Holm, 2002
    beta=1.0/24.0
    mesh_tensor_quad = alpha*beta*length_scale_tensor(dshape_field_in, shape_field_in)

    ! Local assembly:
    field_in_mat=dshape_tensor_dshape(dshape_field_in, mesh_tensor_quad, &
         dshape_field_in, detwei) + shape_shape(shape_field_in&
         &,shape_field_in, detwei)

    lrhsfield=shape_rhs(shape_field_in, field_in_quad*detwei)

    ! Global assembly:
    call addto(M, ele_field_in, ele_field_in, field_in_mat)

    call addto(rhsfield, ele_field_in, lrhsfield)

  end subroutine assemble_anisotropic_smooth_scalar

  subroutine assemble_anisotropic_smooth_vector(M, rhsfield, positions, field_in, alpha, ele)
    type(csr_matrix), intent(inout) :: M
    type(vector_field), intent(inout) :: RHSFIELD
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: field_in
    real, intent(in) :: alpha
    integer, intent(in) :: ele
    real :: beta

    real, dimension(positions%dim,ele_ngi(positions,ele))                        :: field_in_quad
    real, dimension(positions%dim,ele_ngi(positions,ele))                        :: X_quad
    real,dimension(positions%dim,positions%dim,ele_ngi(positions,ele))           :: mesh_tensor_quad
    real, dimension(ele_loc(field_in,ele), ele_ngi(field_in,ele), positions%dim) :: dshape_field_in
    real, dimension(ele_ngi(positions,ele))                                      :: detwei
    integer, dimension(:), pointer                                               :: ele_field_in
    type(element_type), pointer                                                  :: shape_field_in
    real, dimension(ele_loc(field_in, ele), ele_loc(field_in, ele))              :: field_in_mat
    real, dimension(positions%dim, ele_loc(field_in, ele))                       :: lrhsfield

    ele_field_in=>ele_nodes(field_in, ele)
    shape_field_in=>ele_shape(field_in, ele)

    ! Locations of quadrature points. Fields evaluated at quads.
    X_quad=ele_val_at_quad(positions, ele)
    field_in_quad = ele_val_at_quad(field_in, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, shape_field_in, dshape&
         &=dshape_field_in, detwei=detwei)

    ! mesh size tensor = gamma * (mesh size)**2
    ! gamma is a function of mesh shape e.g. hexahedral, unstructured triangular etc.
    ! Helmholtz smoothing factor = alpha * 1/24 * mesh size tensor
    ! Factor of 1/24 comes from Geurts & Holm, 2002
    beta=1.0/24.0
    mesh_tensor_quad = alpha*beta*length_scale_tensor(dshape_field_in, shape_field_in)

    ! Local assembly:
    field_in_mat=dshape_tensor_dshape(dshape_field_in, mesh_tensor_quad, &
         dshape_field_in, detwei) + shape_shape(shape_field_in&
         &,shape_field_in, detwei)

    lrhsfield=shape_vector_rhs(shape_field_in, field_in_quad, detwei)

    ! Global assembly:
    call addto(M, ele_field_in, ele_field_in, field_in_mat)

    call addto(rhsfield, ele_field_in, lrhsfield)

  end subroutine assemble_anisotropic_smooth_vector

  subroutine assemble_anisotropic_smooth_tensor(M, rhsfield, positions, field_in, alpha, ele)
    type(csr_matrix), intent(inout) :: M
    type(tensor_field), intent(inout) :: RHSFIELD
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(in) :: field_in
    real, intent(in) :: alpha
    integer, intent(in) :: ele
    real :: beta

    real, dimension(positions%dim,positions%dim,ele_ngi(positions,ele))          :: field_in_quad
    real, dimension(positions%dim,ele_ngi(positions,ele))                        :: X_quad
    real,dimension(positions%dim,positions%dim,ele_ngi(positions,ele))           :: mesh_tensor_quad
    real, dimension(ele_loc(field_in,ele), ele_ngi(field_in,ele), positions%dim) :: dshape_field_in
    real, dimension(ele_ngi(positions,ele))                                      :: detwei    
    integer, dimension(:), pointer                                               :: ele_field_in
    type(element_type), pointer                                                  :: shape_field_in
    real, dimension(ele_loc(field_in, ele), ele_loc(field_in, ele))              :: field_in_mat
    real, dimension(positions%dim, positions%dim, ele_loc(field_in, ele))        :: lrhsfield

    ele_field_in=>ele_nodes(field_in, ele)
    shape_field_in=>ele_shape(field_in, ele)

    ! Locations of quadrature points. Fields evaluated at quads.
    X_quad=ele_val_at_quad(positions, ele)
    field_in_quad = ele_val_at_quad(field_in, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, shape_field_in, dshape&
         &=dshape_field_in, detwei=detwei)

    ! mesh size tensor = gamma * (mesh size)**2
    ! gamma is a function of mesh geometry
    ! Helmholtz smoothing factor = alpha * 1/24 * mesh size tensor
    ! Factor of 1/24 comes from Geurts & Holm, 2002
    beta=1.0/24.0
    mesh_tensor_quad = alpha*beta*length_scale_tensor(dshape_field_in, shape_field_in)

    ! Local assembly: (1+alpha^2.M) or (1-alpha^2.M)?
    field_in_mat=dshape_tensor_dshape(dshape_field_in, mesh_tensor_quad, &
         dshape_field_in, detwei) + shape_shape(shape_field_in&
         &,shape_field_in, detwei)

    lrhsfield=shape_tensor_rhs(shape_field_in, field_in_quad, detwei)

    ! Global assembly:
    call addto(M, ele_field_in, ele_field_in, field_in_mat)

    call addto(rhsfield, ele_field_in, lrhsfield)

  end subroutine assemble_anisotropic_smooth_tensor

  function length_scale_tensor(du_t, shape) result(t)
    !! Computes a length scale tensor to be used in LES (units are in length^2)
    !! derivative of velocity shape function (nloc x ngi x dim)
    real, dimension(:,:,:), intent(in):: du_t
    !! the resulting tensor (dim x dim x ngi)
    real, dimension(size(du_t,3),size(du_t,3),size(du_t,2)) :: t
    !! for a simplex if degree==1 the tensor is the same for all gaussian points
    type(element_type), intent(in):: shape

    real, dimension(size(t,1), size(t,2)):: M
    real r
    integer gi, loc, i
    integer dim, ngi, nloc, compute_ngi

    t=0.0

    nloc=size(du_t,1)
    ngi=size(du_t,2)
    dim=size(du_t,3)

    if (shape%degree<=1 .and. shape%numbering%family==FAMILY_SIMPLEX) then
       compute_ngi=1
    else
       compute_ngi=ngi
    end if

    do gi=1, compute_ngi
       do loc=1, nloc
          M=outer_product( du_t(loc,gi,:), du_t(loc,gi,:) )
          r=sum( (/ ( M(i,i), i=1, dim) /) )
          t(:,:,gi)=t(:,:,gi)+M/(r**2)
       end do
    end do

    ! copy the rest
    do gi=compute_ngi+1, ngi
       t(:,:,gi)=t(:,:,1)
    end do

  end function length_scale_tensor

end module smoothing_module
