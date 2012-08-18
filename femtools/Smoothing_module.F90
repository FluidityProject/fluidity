#include "fdebug.h"

module smoothing_module
  use state_module
  use fields
  use sparse_tools
  use sparsity_patterns
  use solvers
  use metric_tools
  use boundary_conditions, only: apply_dirichlet_conditions, get_entire_boundary_condition
  use global_parameters, only : OPTION_PATH_LEN
  use elements
  use transform_elements
  use merge_tensors
  use vector_tools
  use unittest_tools
  implicit none

  private
  
  public :: smooth_scalar, smooth_vector, smooth_tensor
  public :: anisotropic_smooth_scalar, anisotropic_smooth_vector
  public :: anisotropic_smooth_tensor, length_scale, average_length_scale, length_scale_tensor
  public :: length_scale_tensor_isotropic, equivalent_radius_ellipsoid
contains

  subroutine smooth_scalar(field_in,positions,field_out,alpha, path)

    !smoothing length tensor
    real, intent(in) :: alpha
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

    ! Boundary conditions
    ewrite(2,*) "Applying strong Dirichlet boundary conditions to filtered field"
    call apply_dirichlet_conditions(M, rhsfield, field_in)

    call zero(field_out)
    call petsc_solve(field_out, M, rhsfield, option_path=trim(path))

    call deallocate(rhsfield)
    call deallocate(M)

  end subroutine smooth_scalar

  subroutine smooth_vector(field_in,positions,field_out,alpha,path,bc_field)

    !smoothing length tensor
    real, intent(in) :: alpha
    !input field
    type(vector_field), intent(inout) :: field_in
    !coordinates field
    type(vector_field), intent(in) :: positions
    !output field, should have same mesh as input field
    type(vector_field), intent(inout) :: field_out
    character(len=*), intent(in) :: path
    ! For dynamic LES, need to specify the velocity field from which BCs are taken
    type(vector_field), intent(inout), optional :: bc_field
    
    !local variables
    type(csr_matrix) :: M
    type(csr_sparsity) :: M_sparsity
    type(vector_field) :: RHSFIELD
    integer :: ele, dim

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

    ! Boundary conditions
    ewrite(2,*) "Applying strong Dirichlet boundary conditions to filtered field:"
    if(present(bc_field)) then
      ewrite(2,*) trim(bc_field%name)
      do dim=1, bc_field%dim
        call apply_dirichlet_conditions(matrix=M, rhs=rhsfield, field=bc_field, dim=dim)
      end do
    else
      ewrite(2,*) trim(field_in%name)
      do dim=1, field_in%dim
        call apply_dirichlet_conditions(matrix=M, rhs=rhsfield, field=field_in, dim=dim)
      end do
    end if

    call zero(field_out)
    call petsc_solve(field_out, M, rhsfield, option_path=trim(path))

    call deallocate(rhsfield)
    call deallocate(M)

  end subroutine smooth_vector

  subroutine smooth_tensor(field_in,positions,field_out,alpha, path)

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
    type(tensor_field) :: rhsfield
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
       call assemble_smooth_tensor(M, rhsfield, positions, field_in, alpha, ele)
    end do

    call zero(field_out)
    call petsc_solve(field_out, M, rhsfield, option_path=trim(path))

    call deallocate(rhsfield)
    call deallocate(M)

  end subroutine smooth_tensor

  subroutine anisotropic_smooth_scalar(field_in,positions,velocity,field_out,alpha,path)

    !smoothing length tensor
    real, intent(in)                  :: alpha
    !input field
    type(scalar_field), intent(inout) :: field_in
    !coordinates and velocity fields
    type(vector_field), pointer, intent(in) :: positions, velocity
    !output field, should have same mesh as input field
    type(scalar_field), intent(inout) :: field_out
    character(len=*), intent(in)      :: path
    
    !local variables
    type(csr_matrix) :: M
    type(csr_sparsity) :: M_sparsity
    type(scalar_field) :: rhsfield
    integer :: ele
    integer, dimension(:), pointer :: neighs

    !allocate smoothing matrix, RHS
    M_sparsity=make_sparsity(field_in%mesh, field_in%mesh, name='HelmholtzScalarSparsity')
    call allocate(M, M_sparsity, name="HelmholtzScalarSmoothingMatrix")
    call allocate(rhsfield, field_in%mesh, "HelmholtzScalarSmoothingRHS")
    call zero(M); call zero(rhsfield); call zero(field_out)

    ! Assemble M element by element.
    do ele=1, element_count(field_in)
      call assemble_anisotropic_smooth_scalar(M, rhsfield, positions, field_in, alpha, ele)
    end do

    ! Boundary conditions
    ewrite(2,*) "Applying strong Dirichlet boundary conditions to filtered field"
    call apply_dirichlet_conditions(M, rhsfield, field_in)

    call petsc_solve(field_out, M, rhsfield, option_path=trim(path))

    call deallocate(rhsfield); call deallocate(M); call deallocate(M_sparsity)

  end subroutine anisotropic_smooth_scalar

  subroutine anisotropic_smooth_vector(field_in,positions,field_out,alpha,path,bc_field)

    !smoothing length tensor
    real, intent(in) :: alpha
    !input field
    type(vector_field), intent(inout) :: field_in
    !coordinates field
    type(vector_field), intent(in) :: positions
    !output field, should have same mesh as input field
    type(vector_field), intent(inout) :: field_out
    character(len=*), intent(in) :: path
    ! For dynamic LES, need to specify the velocity field from which BCs are taken
    type(vector_field), intent(inout), optional :: bc_field
    ! For weak BCs, this is the surface field containing BC values
    type(vector_field) :: field_bc
    ! For weak BCs, arrays of size(field_bc)
    integer, dimension(:,:), allocatable :: field_bc_type, field_bc_number

    !local variables
    type(csr_matrix) :: M
    type(csr_sparsity) :: M_sparsity
    type(vector_field) :: rhsfield
    integer :: ele, sele, dim
    integer, dimension(:), pointer :: neighs

    !allocate smoothing matrix
    M_sparsity=make_sparsity(field_in%mesh, field_in%mesh, name='HelmholtzVectorSparsity')
    call allocate(M, M_sparsity, name="HelmholtzVectorSmoothingMatrix")
    call allocate(rhsfield, field_in%dim, field_in%mesh, "HelmholtzVectorSmoothingRHS")
    call zero(M); call zero(rhsfield); call zero(field_out)

    ! Assemble M element by element.
    do ele=1, element_count(field_in)
      call assemble_anisotropic_smooth_vector(M, rhsfield, positions, field_in, alpha, ele)
    end do

    ! Boundary conditions
    ewrite(2,*) "Applying strong Dirichlet boundary conditions to filtered field:"
    if(present(bc_field)) then
      ewrite(2,*) trim(bc_field%name)
      do dim=1, bc_field%dim
        call apply_dirichlet_conditions(matrix=M, rhs=rhsfield, field=bc_field, dim=dim)
      end do
    else
      ewrite(2,*) trim(field_in%name)
      do dim=1, field_in%dim
        call apply_dirichlet_conditions(matrix=M, rhs=rhsfield, field=field_in, dim=dim)
      end do
    end if

    ! Weak Dirichlet conditions
    !allocate(field_bc_type(field_out%dim, surface_element_count(field_out)))
    !allocate(field_bc_number(field_out%dim, surface_element_count(field_out)))
    !call get_entire_boundary_condition(field_out, (/ "weakdirichlet" /), field_bc, field_bc_type, field_bc_number)
    !ewrite(2,*) 'applying weak Dirichlet BCs to filtered field: ', size(field_bc_number)
    !do sele=1, surface_element_count(field_out)
       !ele = face_ele(field_out, sele)
       !call assemble_anisotropic_smooth_vector_surface(M, rhsfield, positions, field_in, alpha, ele, sele, field_bc)
    !end do

    call petsc_solve(field_out, M, rhsfield, option_path=trim(path))

    call deallocate(rhsfield); call deallocate(M); call deallocate(M_sparsity)
    !call deallocate(field_bc); deallocate(field_bc_type); deallocate(field_bc_number)

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
    type(tensor_field) :: rhsfield
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
    real, intent(in) :: alpha
    integer, intent(in) :: ele

    ! value of field_in at quad points
    real, dimension(ele_ngi(positions,ele)) :: field_in_quad
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
    real :: w
    integer :: i,j

    ele_field_in=>ele_nodes(field_in, ele)
    shape_field_in=>ele_shape(field_in, ele)
    field_in_quad = ele_val_at_quad(field_in, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, shape_field_in, dshape&
         &=dshape_field_in, detwei=detwei)

    ! Calculate filter width (using Deardorff's proposal?)
    ! width = (volume)^(1/3)
    w = average_length_scale(positions, ele)
    alpha_quad=0.
    !value of tensor at quads
    forall(i=1:ele_ngi(positions,ele))
       forall(j=1:positions%dim)
          alpha_quad(j,j,i) = alpha**2*w/24.
       end forall
    end forall
    ! Test of isotropic tensor filter width definition
    alpha_quad=0.
    alpha_quad = alpha**2 / 24. * length_scale_tensor_isotropic(dshape_field_in, shape_field_in)

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
    real, intent(in) :: alpha
    integer, intent(in) :: ele

    ! value of field_in at quad points
    real, dimension(positions%dim,ele_ngi(positions,ele)) :: field_in_quad
    ! Locations of quadrature points
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
    real :: w
    integer :: i,j

    ele_field_in=>ele_nodes(field_in, ele)
    shape_field_in=>ele_shape(field_in, ele)
    field_in_quad = ele_val_at_quad(field_in, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, shape_field_in, dshape&
         &=dshape_field_in, detwei=detwei)

    ! Calculate filter width (using Deardorff's proposal?)
    ! width = (volume)^(1/3)
    w = length_scale(positions, ele)
    !value of tensor at quads
    alpha_quad=0.
    forall(i=1:ele_ngi(positions,ele))
       forall(j=1:positions%dim)
          alpha_quad(j,j,i) = alpha**2*w/24.
       end forall
    end forall
    alpha_quad=0.
    alpha_quad = alpha**2 / 24. * length_scale_tensor_isotropic(dshape_field_in, shape_field_in)

    ! Local assembly:    assert(field_in%mesh==field_out%mesh))
    field_in_mat=dshape_tensor_dshape(dshape_field_in, alpha_quad, &
         dshape_field_in, detwei) + shape_shape(shape_field_in&
         &,shape_field_in, detwei)

    lrhsfield=shape_vector_rhs(shape_field_in, field_in_quad, detwei)

    ! Global assembly:
    call addto(M, ele_field_in, ele_field_in, field_in_mat)

    call addto(rhsfield, ele_field_in, lrhsfield)

  end subroutine assemble_smooth_vector

  subroutine assemble_smooth_tensor(M, rhsfield, positions, field_in, alpha, ele)

    type(csr_matrix), intent(inout) :: M
    type(tensor_field), intent(inout) :: RHSFIELD
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(in) :: field_in
    real, intent(in) :: alpha
    integer, intent(in) :: ele

    ! value of field_in at quad points
    real, dimension(positions%dim,positions%dim,ele_ngi(positions,ele)) :: field_in_quad
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
    real, dimension(positions%dim,positions%dim,ele_loc(field_in, ele)) :: lrhsfield
    real :: w
    integer :: i,j

    ele_field_in=>ele_nodes(field_in, ele)
    shape_field_in=>ele_shape(field_in, ele)
    field_in_quad = ele_val_at_quad(field_in, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, shape_field_in, dshape&
         &=dshape_field_in, detwei=detwei)

    ! Calculate filter width (using Deardorff's proposal?)
    ! width = (volume)^(1/3)
    w = length_scale(positions, ele)

    !value of tensor at quads
    alpha_quad=0.
    forall(i=1:ele_ngi(positions,ele))
       forall(j=1:positions%dim)
          alpha_quad(j,j,i) = alpha**2*w/24.
       end forall
    end forall

    ! Local assembly:    assert(field_in%mesh==field_out%mesh))
    field_in_mat=dshape_tensor_dshape(dshape_field_in, alpha_quad, &
         dshape_field_in, detwei) + shape_shape(shape_field_in&
         &,shape_field_in, detwei)

    lrhsfield=shape_tensor_rhs(shape_field_in, field_in_quad, detwei)

    ! Global assembly:
    call addto(M, ele_field_in, ele_field_in, field_in_mat)

    call addto(rhsfield, ele_field_in, lrhsfield)

  end subroutine assemble_smooth_tensor

  subroutine assemble_anisotropic_smooth_scalar(M, rhsfield, positions, field_in, alpha, ele)
    type(csr_matrix), intent(inout) :: M
    type(scalar_field), intent(inout) :: rhsfield
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: field_in
    real, intent(in) :: alpha
    integer, intent(in) :: ele
    real, dimension(ele_ngi(positions,ele))                                      :: field_in_quad
    real,dimension(positions%dim,positions%dim,ele_ngi(positions,ele))           :: mesh_tensor_quad
    real, dimension(ele_loc(field_in,ele), ele_ngi(field_in,ele), positions%dim) :: dshape_field_in
    real, dimension(ele_ngi(positions,ele))                                      :: detwei
    integer, dimension(:), pointer                                               :: ele_field_in
    type(element_type), pointer                                                  :: shape_field_in, f_shape
    real, dimension(ele_loc(field_in, ele), ele_loc(field_in, ele))              :: field_in_mat
    real, dimension(ele_loc(field_in, ele))                                      :: lrhsfield
    !real, dimension(positions%dim,face_ngi(field_in,sele))                      :: normal_bdy
    !real, dimension(face_ngi(field_in,sele))                                    :: detwei_bdy
    !real, dimension(positions%dim, positions%dim, ele_ngi(field_in, sele))        :: invJ
    !real, dimension(positions%dim, positions%dim, face_ngi(field_in, sele))       :: invJ_face
    !type(element_type)                                     :: augmented_shape
    !real, dimension(ele_loc(field_in, ele), face_ngi(field_in, sele),positions%dim):: vol_dshape_face
    !integer                       :: l_face_number
    !real, dimension(face_loc(field_in, sele), face_loc(field_in, sele)) :: face_mat

    ele_field_in=>ele_nodes(field_in, ele)
    shape_field_in=>ele_shape(field_in, ele)
    field_in_quad = ele_val_at_quad(field_in, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, shape_field_in, dshape&
         &=dshape_field_in, detwei=detwei)

    ! mesh size tensor=(edge lengths)**2
    ! Helmholtz smoothing lengthscale = alpha**2 * 1/24 * mesh size tensor
    ! factor 1/24 derives from 2nd moment of filter (see Pope 2000, Geurts&Holm 2002)
    mesh_tensor_quad = alpha**2 / 24. * length_scale_tensor(dshape_field_in, shape_field_in)

    ! face wizardry if on a boundary
    !if(sele>0) then
      !f_shape   =>face_shape(field_in, sele)
      !l_face_number = local_face_number(field_in, sele)
      !call transform_facet_to_physical(positions, sele, detwei_f=detwei_bdy, normal=normal_bdy)
      !augmented_shape = make_element_shape(shape_field_in%loc, shape_field_in%dim, &
      !                  shape_field_in%degree, shape_field_in%quadrature, quad_s=f_shape%quadrature )
      !call compute_inverse_jacobian( ele_val(positions, ele), shape_field_in, invJ )
      !invJ_face = spread(invJ(:, :, 1), 3, size(invJ_face, 3))
      ! nloc x sngi x dim
      !vol_dshape_face = eval_volume_dshape_at_face_quad(augmented_shape, l_face_number, invJ_face)
      ! Local assembly with Neumann term
      !field_in_mat=dshape_tensor_dshape(dshape_field_in, mesh_tensor_quad, dshape_field_in, detwei) &
      !             + shape_shape(shape_field_in,shape_field_in, detwei)
      !             !+ dshape_dot_vector_shape(dshape_field_in, normal_bdy, shape_field_in, detwei)
      !call deallocate(augmented_shape)
    !else
      ! Local assembly if not on a boundary
    field_in_mat=dshape_tensor_dshape(dshape_field_in, mesh_tensor_quad, dshape_field_in, detwei) &
                 + shape_shape(shape_field_in,shape_field_in, detwei)
    !end if

    lrhsfield=shape_rhs(shape_field_in, field_in_quad*detwei)

    ! Global assembly
    call addto(M, ele_field_in, ele_field_in, field_in_mat)

    call addto(rhsfield, ele_field_in, lrhsfield)

  end subroutine assemble_anisotropic_smooth_scalar

  subroutine assemble_anisotropic_smooth_vector(M, rhsfield, positions, field_in, alpha, ele)
    type(csr_matrix), intent(inout) :: M
    type(vector_field), intent(inout) :: rhsfield
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: field_in
    real, intent(in) :: alpha
    integer, intent(in) :: ele
    integer :: dim
    real, dimension(positions%dim,ele_ngi(positions,ele))                        :: field_in_quad
    real,dimension(positions%dim,positions%dim,ele_ngi(field_in,ele))           :: mesh_tensor_quad
    real, dimension(ele_loc(field_in,ele), ele_ngi(field_in,ele), positions%dim) :: dshape_field_in
    real, dimension(ele_ngi(field_in,ele))                                      :: detwei
    integer, dimension(:), pointer                                               :: ele_field_in
    type(element_type), pointer                                                  :: shape_field_in, f_shape
    real, dimension(ele_loc(field_in, ele), ele_loc(field_in, ele))              :: field_in_mat
    real, dimension(positions%dim, ele_loc(field_in, ele))                       :: lrhsfield
    !real, dimension(positions%dim,face_ngi(field_in,sele))                      :: normal_bdy
    !real, dimension(face_ngi(field_in,sele))                                    :: detwei_bdy
    !real, dimension(positions%dim, positions%dim, ele_ngi(field_in, sele))        :: invJ
    !real, dimension(positions%dim, positions%dim, face_ngi(field_in, sele))       :: invJ_face
    !type(element_type)                                     :: augmented_shape
    !real, dimension(ele_loc(field_in, ele), face_ngi(field_in, sele),positions%dim):: vol_dshape_face
    !integer                       :: l_face_number

    ele_field_in=>ele_nodes(field_in, ele)
    shape_field_in=>ele_shape(field_in, ele)
    field_in_quad = ele_val_at_quad(field_in, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, shape_field_in, dshape&
         &=dshape_field_in, detwei=detwei)

    ! mesh size tensor=(edge lengths)**2
    ! Helmholtz smoothing lengthscale = alpha**2 * 1/24 * mesh size tensor
    mesh_tensor_quad = alpha**2 / 24. * length_scale_tensor(dshape_field_in, shape_field_in)

    ! face wizardry in on a boundary
    !if(sele>0) then
      !f_shape   =>face_shape(field_in, sele)
      !l_face_number = local_face_number(field_in, sele)
      !call transform_facet_to_physical(positions, sele, detwei_f=detwei_bdy, normal=normal_bdy)
      !augmented_shape = make_element_shape(shape_field_in%loc, shape_field_in%dim, &
      !                   shape_field_in%degree, shape_field_in%quadrature, quad_s=f_shape%quadrature )
      !call compute_inverse_jacobian( ele_val(positions, ele), shape_field_in, invJ )
      !invJ_face = spread(invJ(:, :, 1), 3, size(invJ_face, 3))
      ! nloc x sngi x dim
      !vol_dshape_face = eval_volume_dshape_at_face_quad(augmented_shape, l_face_number, invJ_face)
      ! Local assembly with Neumann term
      !field_in_mat=dshape_tensor_dshape(dshape_field_in, mesh_tensor_quad, dshape_field_in, detwei) &
      !             + shape_shape(shape_field_in,shape_field_in, detwei)
      !             !- dshape_dot_vector_shape(vol_dshape_face, normal_bdy, f_shape, detwei_bdy)
      !call deallocate(augmented_shape)
    !else
      ! Local assembly if not on a boundary
    field_in_mat=dshape_tensor_dshape(dshape_field_in, mesh_tensor_quad, dshape_field_in, detwei) &
                 + shape_shape(shape_field_in,shape_field_in, detwei)
    !end if

    lrhsfield=shape_vector_rhs(shape_field_in, field_in_quad, detwei)

    ! Global assembly:
    call addto(M, ele_field_in, ele_field_in, field_in_mat)
    call addto(rhsfield, ele_field_in, lrhsfield)

  end subroutine assemble_anisotropic_smooth_vector

  subroutine assemble_anisotropic_smooth_vector_surface(M, rhsfield, positions, field_in, alpha, ele, sele, field_bc)
    type(csr_matrix), intent(inout) :: M
    type(vector_field), intent(inout) :: rhsfield
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: field_in, field_bc
    real, intent(in) :: alpha
    integer, intent(in) :: ele, sele

    integer :: i, j, dim, l_face_number, inode
    real, dimension(positions%dim,face_ngi(positions,sele))                      :: field_in_quad
    real,dimension(positions%dim,positions%dim,ele_ngi(field_in,ele))          :: mesh_tensor_quad
    type(element_type), pointer                                                  :: shape_field_in, shape_face
    type(element_type)                                                            :: augmented_shape
    real, dimension(positions%dim, ele_loc(field_in, ele))                       :: bc_vals
    real, dimension(ele_loc(field_in, ele), ele_loc(field_in, ele))          :: field_in_mat
    real, dimension(positions%dim,face_ngi(field_in,sele))                      :: normal_bdy
    real, dimension(face_ngi(field_in,sele))                                    :: detwei_bdy
    real, dimension(ele_ngi(field_in,ele))                                    :: detwei
    integer, dimension(face_loc(field_in, sele))                                 :: nodes_bdy
    integer, dimension(ele_loc(field_in, ele))                                 :: nodes_ele
    real, dimension(positions%dim, positions%dim, ele_ngi(field_in, sele))        :: invJ
    real, dimension(positions%dim, positions%dim, face_ngi(field_in, sele))       :: invJ_face
    real, dimension(ele_loc(field_in,ele), ele_ngi(field_in,ele), positions%dim) :: dshape_field_in
    real, dimension(ele_loc(field_in, ele), face_ngi(field_in, sele),positions%dim):: dshape_face
    integer, dimension(:), pointer       :: faces_ele
    shape_field_in => ele_shape(field_in, ele)
    shape_face => face_shape(field_in, sele)
    nodes_bdy = face_global_nodes(field_in, sele)
    nodes_ele = ele_nodes(field_in, ele)
    !ewrite(2,*) 'nodes_bdy, nodes_ele: ', nodes_bdy, nodes_ele
    !ewrite(2,*) 'nodes_sele: ', ele_nodes(field_bc, sele)
    !ewrite(2,*) 'local_nodes: ', face_local_nodes(field_in, sele)
    l_face_number = local_face_number(field_in, sele)
    !ewrite(2,*) 'l_face_number: ', l_face_number
    field_in_quad = face_val_at_quad(field_in, sele)
    !ewrite(2,*) 'face shape: ', shape_face%ngi, shape_face%dim, shape_face%loc
    !ewrite(2,*) 'ele shape: ', shape_field_in%ngi, shape_field_in%dim, shape_field_in%loc

    call transform_to_physical(positions, ele, shape=shape_field_in, dshape=dshape_field_in, detwei=detwei, invJ=invJ)
    call transform_facet_to_physical(positions, sele, detwei_f=detwei_bdy, normal=normal_bdy)
    augmented_shape = make_element_shape(shape_field_in%loc, shape_field_in%dim, &
                      & shape_field_in%degree, shape_field_in%quadrature, quad_s=shape_face%quadrature)

    ! From Surface_Integrals.F90:
    if (shape_field_in%degree == 1 .and. shape_field_in%numbering%family == FAMILY_SIMPLEX) then
       invJ_face = spread(invJ(:, :, 1), 3, size(invJ_face, 3))
    else
       ewrite(-1,*) "If positions are nonlinear, then you have to compute"
       ewrite(-1,*) "the inverse Jacobian of the volume element at the surface"
       ewrite(-1,*) "quadrature points. Sorry ..."
       FLExit("Calculating the body drag not supported for nonlinear coordinates.")
    end if

    ! Get dshape transformed onto surface element: nloc x sngi x dim
    dshape_face = eval_volume_dshape_at_face_quad(augmented_shape, l_face_number, invJ_face)
    !ewrite(2,*) 'face dshape: ', size(dshape_face,1), size(dshape_face,2), size(dshape_face,3)

    ! mesh size tensor=(edge lengths)**2
    ! Helmholtz smoothing lengthscale = alpha**2 * 1/24 * mesh size tensor
    ! 2D OR 3D?
    mesh_tensor_quad = alpha**2 / 24. * length_scale_tensor(dshape_field_in, shape_field_in)
    !ewrite(2,*) 'mesh tensor: ', size(mesh_tensor_quad,1), size(mesh_tensor_quad,2), size(mesh_tensor_quad,3)
    !do gi=1,sngi
    !   filter_gi(:,:,gi) = matmul(dshape_face(:,gi,:), matmul(mesh_tensor_quad, dshape_face(:,gi,:)))
    !end do
    !field_in_mat = matmul(filter_gi, detwei_bdy)

    field_in_mat = dshape_tensor_dshape(dshape_field_in, mesh_tensor_quad, dshape_field_in, detwei) &
                   & + shape_shape(shape_field_in, shape_field_in, detwei)
!    field_in_mat = shape_shape(shape_face, shape_face, detwei_bdy)

    ! Vector length (nloc) comprising bc values and off-wall value:
    bc_vals=0.
!    do i = 1, size(nodes_bdy)
!       do j = 1, size(nodes_ele)
!          ewrite(2,*) 'i,j,nodes_ele,nodes_bdy: ', i,j,nodes_bdy(i),nodes_ele(j)
!          if (nodes_ele(j)==nodes_bdy(i)) then
!             do dim=1,positions%dim
!                bc_vals(dim,j)=node_val(field_bc, dim, nodes_bdy(i))
!             end do
!          end if
!       end do
!    end do
    do dim=1,positions%dim
       bc_vals(dim,face_local_nodes(field_in, sele)) = ele_val(field_bc, dim, sele)
    end do
    ! find element node inside the domain
    faces_ele => ele_faces(field_in, ele)
    do i = 1, size(faces_ele)
       if (faces_ele(i)==sele) exit
    end do
    inode = nodes_ele(i)
    do dim=1,positions%dim
       !bc_vals(dim,i)=sum(ele_val(field_bc, dim, sele))/size(nodes_bdy)
       bc_vals(dim,i)=node_val(field_in,dim,inode)
    end do
    ewrite(2,*) 'field_bc: ', ele_val(field_bc, sele)
    ewrite(2,*) 'bc_vals: ', bc_vals
    ! Global assembly of weak Dirichlet BC:
    do dim = 1, field_in%dim
       call addto(rhsfield, dim, ele_nodes(field_in, ele), -matmul(field_in_mat, bc_vals(dim,:)))
    end do

  end subroutine assemble_anisotropic_smooth_vector_surface

  subroutine assemble_anisotropic_smooth_tensor(M, rhsfield, positions, field_in, alpha, ele)
    type(csr_matrix), intent(inout) :: M
    type(tensor_field), intent(inout) :: rhsfield
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(in) :: field_in
    real, intent(in) :: alpha
    integer, intent(in) :: ele
    real, dimension(positions%dim,positions%dim,ele_ngi(positions,ele))          :: field_in_quad
    real,dimension(positions%dim,positions%dim,ele_ngi(positions,ele))           :: mesh_tensor_quad
    real, dimension(ele_loc(field_in,ele), ele_ngi(field_in,ele), positions%dim) :: dshape_field_in
    real, dimension(ele_ngi(positions,ele))                                      :: detwei    
    integer, dimension(:), pointer                                               :: ele_field_in
    type(element_type), pointer                                                  :: shape_field_in
    real, dimension(ele_loc(field_in, ele), ele_loc(field_in, ele))              :: field_in_mat
    real, dimension(positions%dim, positions%dim, ele_loc(field_in, ele))        :: lrhsfield

    ele_field_in=>ele_nodes(field_in, ele)
    shape_field_in=>ele_shape(field_in, ele)
    field_in_quad = ele_val_at_quad(field_in, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, shape_field_in, dshape&
         &=dshape_field_in, detwei=detwei)

    ! mesh size tensor=(edge lengths)**2
    ! Helmholtz smoothing lengthscale = alpha**2 * 1/24 * mesh size tensor
    mesh_tensor_quad = alpha**2 / 24. * length_scale_tensor(dshape_field_in, shape_field_in)

    ! Local assembly:
    field_in_mat=dshape_tensor_dshape(dshape_field_in, mesh_tensor_quad, &
         dshape_field_in, detwei) + shape_shape(shape_field_in&
         &,shape_field_in, detwei)

    lrhsfield=shape_tensor_rhs(shape_field_in, field_in_quad, detwei)

    ! Global assembly:
    call addto(M, ele_field_in, ele_field_in, field_in_mat)

    call addto(rhsfield, ele_field_in, lrhsfield)

  end subroutine assemble_anisotropic_smooth_tensor

  function length_scale(positions, ele) result(s)
    ! Computes a scalar length scale for LES models
    ! Preserves element volume. (units are in length^2)
    type(vector_field), intent(in) :: positions
    real :: s
    integer, intent(in) :: ele
    integer :: dim
    dim=positions%dim

    ! filter is a square/cube (width=side length) a la Deardorff:
    s=element_volume(positions, ele)
    s=s**(2./dim)

  end function length_scale

  function average_length_scale(positions, ele) result(s)
    ! Computes a scalar length scale for LES models
    ! Preserves element volume. (units are in length^2)
    ! Averages over neighbouring elements to reduce anisotropy
    type(vector_field), intent(in) :: positions
    real :: s
    integer, intent(in) :: ele
    integer :: dim, i, ele2
    integer, dimension(:), pointer :: patch

    dim=positions%dim
    s=element_volume(positions, ele)**(2./dim)
    ewrite(2,*) 's', s
    patch => ele_neigh(positions, ele)
    ewrite(2,*) 'patch', patch
    do i=1, size(patch)
      ele2=patch(i)
      if (ele2>0) then
        ewrite(2,*) 'ele2', ele2
        s=s+element_volume(positions, ele2)**(2./dim)
        ewrite(2,*) 's2', s
      end if
    end do
    
    s=s/size(patch)
    ewrite(2,*) 's3', s

  end function average_length_scale

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
    integer gi, loc, i, dim, nloc, compute_ngi

    t=0.0
    nloc=size(du_t,1)
    dim=size(du_t,3)

    if (.not.(shape%degree==1 .and. shape%numbering%family==FAMILY_SIMPLEX)) then
      ! for non-linear compute on all gauss points
      compute_ngi=shape%ngi
    else
      ! for linear: compute only the first and copy the rest
      compute_ngi=1
    end if

    do gi=1, compute_ngi
       do loc=1, nloc
          ! eigenvalues of metric
          M=outer_product( du_t(loc,gi,:), du_t(loc,gi,:) )
          ! determinant of M
          r=sum( (/ ( M(i,i), i=1, dim) /) )
          ! M^-1 = 1/det(M)*adj(M) = 1/det(M)*M
          ! Why det(M)^2?
          if (.not. r==0.0) then
             t(:,:,gi)=t(:,:,gi)+M/(r**2)
          end if
          ! rotation of M^-1 by V^T V comes later:
          ! construction of viscosity matrix happens in local coords
       end do
    end do

    ! copy the rest
    do gi=compute_ngi+1, shape%ngi
       t(:,:,gi)=t(:,:,1)
    end do

  end function length_scale_tensor

  function length_scale_tensor_isotropic(du_t, shape) result(t2)
    !! Computes a length scale tensor to be used in LES (units are in length^2)
    !! derivative of velocity shape function (nloc x ngi x dim)
    real, dimension(:,:,:), intent(in):: du_t
    !! the resulting tensor (dim x dim x ngi)
    real, dimension(size(du_t,3),size(du_t,3),size(du_t,2)) :: t, t2
    !! for a simplex if degree==1 the tensor is the same for all gaussian points
    type(element_type), intent(in):: shape

    real, dimension(size(t,1), size(t,2)):: M
    real, dimension(size(t, 1), size(t, 1)), target :: v1, v2, Finv, F ! eigenvectors, rotations
    real, dimension(size(t, 1)), target :: a1, a2 ! eigenvalues
    real r
    integer gi, loc, i, dim, nloc, compute_ngi

    t2=0.0
    nloc=size(du_t,1)
    dim=size(du_t,3)

    if (.not.(shape%degree==1 .and. shape%numbering%family==FAMILY_SIMPLEX)) then
      ! for non-linear compute on all gauss points
      compute_ngi=shape%ngi
    else
      ! for linear: compute only the first and copy the rest
      compute_ngi=1
    end if

    ! Get anisotropic tensor
    t = length_scale_tensor(du_t, shape)

    ! Distort metric t to sphere with same volume as ellipsoid
    do gi=1, compute_ngi
       call eigendecomposition_symmetric(t(:,:,gi), v1, a1)
       call vec_clean(a1, 1e-12)
       a2 = equivalent_radius_ellipsoid(dim, a1)
       v2 = get_matrix_identity(dim)
       call eigenrecomposition(t2(:,:,gi), v2, a2)
    end do

    ! copy the rest
    do gi=compute_ngi+1, shape%ngi
       t2(:,:,gi)=t2(:,:,1)
    end do

  end function length_scale_tensor_isotropic

  function equivalent_radius_ellipsoid(dim, evals) result(r)
    ! given the eigenvalues of an ellipsoid, compute the eigenvalue
    ! of the sphere with equal volume
    real, dimension(dim), target :: evals
    real r
    integer i, dim

    r = 1.0
    do i=1, dim
       r = r*evals(i)
    end do
    r = r**(1./dim)

  end function equivalent_radius_ellipsoid

end module smoothing_module
