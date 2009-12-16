c#include "fdebug.h"

program make_cube_test_data
  ! A small program to construct a vtu file
  ! for testing mayavi diagnostics

  use read_triangle
  use fields
  use FEtools
  use DGtools
  use elements
  use sparse_tools
  use vtk_interfaces
  use transform_elements
  use solvers
  use petsc_tools
  use sparsity_patterns
  use signal_vars

  implicit none
  type(vector_field), target :: u, positions
  type(scalar_field), target :: h
  integer :: degree, quad_degree
  type(quadrature_type), target :: quad
  type(element_type), target :: shape
  type(dynamic_csr_matrix) :: u1_inverse_mass, u2_inverse_mass
  type(mesh_type) :: h_mesh,u_mesh
  type(dynamic_csr_matrix) :: C1T, C2T, CMC1, CMC, CM1, CM2

  ! Arguments for handling the command line
  character(len=256) :: filename, degree_buffer
  integer :: status, i, count, dumcount, ndump,ele,j
  integer, dimension(:), pointer :: row_j
  character(len=100) :: dumcount_str, fmt,buffer

  !bc bits
  integer, dimension(:), allocatable :: u1_dirichlet, u2_dirichlet
  integer, dimension(2) :: bc_size
  integer :: dirichlet_flag 

  ewrite(1,*) 'program inspect_CMC'

  call Initialize_Petsc()

  call get_command_argument(1, value=filename, status=status)
  filename=trim('bin/tests/data/'//filename)
  
  select case(status)
  case(1:)
     call usage
     stop
  case(:-1)
     write(0,*) "Warning: truncating filename"
  end select

  ewrite(1,*) 'getting quadrature'

  quad_degree = 4
  quad=make_quadrature(loc=3, dimension=2, degree=quad_degree)
  f_quad=make_quadrature(loc=2, dimension=1, degree=quad_degree-1)
  
  ewrite(1,*) 'Getting shape functions'

  ! Shape functions for positions (linear)
  X_shape=make_element_shape(loc=3, dimension=2, degree=1, quad=quad)

  ewrite(1,*) 'reading mesh'

  positions=read_triangle_files(filename, X_shape)

  ewrite(1,*) 'getting shapes'

  ! Shape functions for velocity and height
  u_shape=make_element_shape(loc=3, dimension=2, degree=1, quad=quad)
  u_face_shape=make_element_shape(loc=2, dimension=1, degree=1, quad=f_quad)
  h_shape=make_element_shape(loc=3, dimension=2, degree=2, quad=quad)

  !connectivity for velocity and height
  ewrite(1,*) 'Getting u connectivity'
  u_mesh = make_mesh(positions%mesh,u_shape,-1,'u_mesh')
  ewrite(1,*) 'Getting h connectivity'
  h_mesh = make_mesh(positions%mesh,h_shape,0,'h_mesh')

  !fields for velocity and height
  call allocate(u,2,u_mesh,'velocity')
  call allocate(h,h_mesh,'height')

  ewrite(1,*) 'getting u boundary conditions'
  call get_u_bcs_size(bc_size,u,positions)
  allocate( u1_dirichlet(bc_size(1)), u2_dirichlet(bc_size(2)) )
  call get_u_bcs(u1_dirichlet, u2_dirichlet,u,positions,bc_size)

  !get inverse mass matrix for u
  ewrite(1,*) 'getting inverse mass matrix for u'

  call allocate(u1_inverse_mass,node_count(u),node_count(u))
  call allocate(u2_inverse_mass,node_count(u),node_count(u))

  dirichlet_flag = DIRICHLET_NONE
  !dirichlet_flag = DIRICHLET_BIG_SPRING
  !dirichlet_flag = DIRICHLET_ONES_ON_DIAGONAL

  if(.false.) then
     call get_dg_inverse_mass_matrix(u1_inverse_mass,u_mesh,positions, &
          u1_dirichlet,dirichlet_flag)
     call get_dg_inverse_mass_matrix(u2_inverse_mass,u_mesh,positions, &
          u2_dirichlet,dirichlet_flag)
  else
     call get_lumped_mass(u1_inverse_mass,u_mesh,positions, &
          u1_dirichlet,dirichlet_flag)
     call get_lumped_mass(u2_inverse_mass,u_mesh,positions, &
          u2_dirichlet,dirichlet_flag)
     do i = 1, node_count(u_mesh)
        u1_inverse_mass%val(i)%ptr = u1_inverse_mass%val(i)%ptr
        u2_inverse_mass%val(i)%ptr = u2_inverse_mass%val(i)%ptr
     end do
  end if

  call allocate(C1T,node_count(h),node_count(u))
  call allocate(C2T,node_count(h),node_count(u))

  call assemble_CT(positions,u,h,C1T,C2T)

  if(dirichlet_flag == DIRICHLET_ONES_ON_DIAGONAL) then
     do i = 1, size(u1_dirichlet)
        call zero_column(C1T,u1_dirichlet(i))
     end do
     do i = 1, size(u2_dirichlet)
        call zero_column(C2T,u2_dirichlet(i))
     end do
  end if

  CM1 = matmul_T(C1T,u1_inverse_mass)
  CM2 = matmul_T(C2T,u2_inverse_mass)
  CMC1 = matmul_T(CM1,C1T)
  CMC = matmul_T(CM2,C2T)

  do i = 1, size(CMC1,1)
     row_j => CMC1%colm(i)%ptr
     do j = 1, size(row_j)
        call addto(CMC,i,row_j(j),CMC1%val(i)%ptr(j))
     end do
  end do

  ewrite(1,*) 'Writing file'
  
  open(file='cmc.dat',action='write',unit=69)
  write(buffer,'(a,i0,a)') '(',size(cmc,1),'g22.8)'
  write(69,buffer) transpose(dense(cmc))
  close(69)
  
  ewrite(1,*) 'END program inspect_CMC'

contains
  
  subroutine assemble_CT(positions,u,h,C1T,C2T)
    type(vector_field), intent(in) :: positions, u
    type(scalar_field), intent(in) :: h
    type(dynamic_csr_matrix), intent(inout) :: C1T,C2T

    !locals
    integer :: ele

    do ele = 1, element_count(u)
       call assemble_CT_elemental(ele,positions,u,h,C1T,C2T)
    end do

  end subroutine assemble_CT
  
  subroutine assemble_CT_elemental(ele,positions,u,h,C1T,C2T)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: u
    type(scalar_field), intent(in) :: h
    type(dynamic_csr_matrix), intent(inout) :: C1T,C2T

    ! Locations of nodes.
    real, dimension(positions%dim,ele_loc(positions,ele)) :: X_ele
    ! Locations of quadrature points
    real, dimension(positions%dim,ele_ngi(positions,ele)) :: X_quad
    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,ele)) :: detwei    
    ! Derivatives of shape function:
    real, dimension(ele_loc(h,ele), &
         ele_ngi(h,ele), positions%dim) :: dshape_h
    ! Node numbers of field element.
    integer, dimension(:), pointer :: ele_u, ele_h
    ! Shape functions.
    type(element_type), pointer :: shape_u, shape_h, shape_X
    ! gradient matrix
    real, dimension(2,ele_loc(h,ele),ele_loc(u,ele)) :: grad_mat

    ele_u=>ele_nodes(u, ele)
    shape_u=>ele_shape(u, ele)
    ele_h=>ele_nodes(h, ele)
    shape_h=>ele_shape(h, ele)
    shape_X=>ele_shape(positions, ele)

    ! Locations of local vertices.
    X_ele=ele_val(positions, ele)

    ! Locations of quadrature points.
    X_quad=ele_val_at_quad(positions, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(X_ele, shape_X, m=shape_h, &
         dm_t=dshape_h, detwei=detwei)

    grad_mat = -dshape_shape(dshape_h, shape_u, detwei)

    call addto(C1T,ele_h,ele_u,grad_mat(1,:,:))
    call addto(C2T,ele_h,ele_u,grad_mat(2,:,:))

  end subroutine assemble_CT_elemental

  subroutine get_u_bcs_size(bc_size,u,positions)
    type(vector_field), intent(in) :: u
    type(vector_field), intent(in) :: positions
    integer, dimension(2), intent(out) :: bc_size

    !locals
    integer :: ele, n_u1_bcs, n_u2_bcs, i, u1_count=0, u2_count=0
    real, dimension(2,3) :: X
    integer, dimension(:), pointer :: u_ele_nodes

    write(1,*) 'Getting size of boundary conditions for u. This is a big hack'
    write(1,*) 'A x-cpt dirichlet condition is set if x<0.0 or x>1.0'
    write(1,*) 'A y-cpt dirichlet condition is set if y<0.0 or y>1.0'
    write(1,*) 'Also only works for dg linear elements'
    
    !count number of u1 and u2 conditions
    bc_size = 0
    do ele = 1, element_count(u)
       u_ele_nodes => ele_nodes(u,ele)
       X = ele_val(positions,ele)

       do i = 1, size(u_ele_nodes)
          if(X(1,i)<0.0.or.X(1,i)>1.0) then
             bc_size(1) = bc_size(1) + 1
          end if
          if(X(2,i)<0.0.or.X(2,i)>1.0) then
             bc_size(2) = bc_size(2) + 1
          end if
       end do
    end do

  end subroutine get_u_bcs_size

  subroutine get_u_bcs(u1_dirichlet, u2_dirichlet,u,positions,bc_size)
    type(vector_field), intent(in) :: u
    type(vector_field), intent(in) :: positions
    integer, intent(in), dimension(2) :: bc_size
    integer, dimension(bc_size(1)), intent(out) :: u1_dirichlet
    integer, dimension(bc_size(2)), intent(out) :: u2_dirichlet

    !locals
    integer :: ele, n_u1_bcs, n_u2_bcs, i, u1_count=0, u2_count=0
    real, dimension(2,3) :: X
    integer, dimension(:), pointer :: u_ele_nodes

    write(1,*) 'Getting boundary conditions for u. This is a big hack'
    write(1,*) 'A x-cpt dirichlet condition is set if x<0.0 or x>1.0'
    write(1,*) 'A y-cpt dirichlet condition is set if y<0.0 or y>1.0'
    write(1,*) 'Also only works for dg linear elements'
    
    u1_count = 0
    u2_count = 0
    do ele = 1, element_count(u)
       u_ele_nodes => ele_nodes(u,ele)
       X = ele_val(positions,ele)
       
       do i = 1, size(u_ele_nodes)
          if(X(1,i)<0.0.or.X(1,i)>1.0) then
             u1_count = u1_count + 1
             u1_dirichlet(u1_count) = u_ele_nodes(i)
          end if
          if(X(2,i)<0.0.or.X(2,i)>1.0) then
             u2_count = u2_count + 1
             u2_dirichlet(u2_count) = u_ele_nodes(i)
          end if
       end do
    end do

    assert(u1_count==bc_size(1))
    assert(u2_count==bc_size(2))

  end subroutine get_u_bcs

  subroutine usage
    
    write (0,*) "usage: inspect_CMC <triangle_file_name>"
    
  end subroutine usage

end program inspect_CMC
