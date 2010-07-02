#include "matrixtypes.h"
#include "fdebug.h"

program solve_CMC
  ! A small program to construct the CMC matrix
  ! using P1dg triangular elements for velocity
  ! and P2 triangular elements for height
  ! We solve              <v,u> = -<v,grad p> + (v.n,p)
  !      -<grad g,u> + < gn, u> = <g,f>
  ! 
  use mesh_files
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
  use solvers
  use global_parameters, only : current_debug_level

  implicit none
  type(vector_field), target :: u, positions, positions_h
  type(scalar_field), target :: h, RHS, h0, error
  integer :: degree, quad_degree
  type(quadrature_type), target :: quad,f_quad
  type(element_type), target :: X_shape, u_shape, h_shape, &
       X_shape_f,u_shape_f,h_shape_f
  type(dynamic_csr_matrix) :: u1_inverse_mass
  type(csr_matrix) :: u_inverse_mass
  type(dynamic_csr_matrix) :: h_mass
  type(mesh_type) :: h_mesh,u_mesh
  type(dynamic_csr_matrix) :: C1T, C2T, C3T
  type(dynamic_csr_matrix) :: CMC, C2MC2T, C3MC3T, &
       MC1T,MC2T,MC3T
  type(csr_matrix) :: CMC_static
  type(csr_matrix) :: h_mass_static
  real, dimension(:), allocatable :: tmp

  ! Arguments for handling the command line
  character(len=256) :: filename, degree_buffer
  integer :: status, i, count, dumcount, ndump,ele,j
  integer, dimension(:), pointer :: row_j
  character(len=100) :: dumcount_str, fmt,buffer,buf

  !bc bits
  integer :: dirichlet_flag, output_flag 
  integer :: u_continuity, h_order
  integer :: dim, loc, nnodes, nelements, node_attributes  
  real,parameter :: pi = 3.141592654
  integer, parameter :: fix_pressure_value = 1

  !number of boundary elements
  integer :: n_boundary_elements = 0
  integer, dimension(:), allocatable :: bc_marker

  !debugging bits
  type(vector_field), target :: positions_u  

  current_debug_level = 3

  ewrite(1,*) 'program solve_CMC'

  call Initialize_Petsc()

  call get_command_argument(1, value=filename, status=status)
  select case(status)
  case(1:)
     call usage
     stop
  case(:-1)
     write(0,*) "Warning: truncating filename"
  end select
  filename=trim(filename)

  u_continuity = 0
  call get_command_argument(2, value=buf, status=status)
  if(status>0) then
     call usage
     stop
  end if
  if(buf=='-1') u_continuity = -1  
  
  h_order = 1
  call get_command_argument(3, value=buf, status=status)
  if(status>0) then
     call usage
     stop
  end if
  if(buf=='2') h_order = 2

  dirichlet_flag = DIRICHLET_NONE
  call get_command_argument(4, value=buf, status=status)
  if(status>0) then
     call usage
     stop
  end if
  select case(buf)
  case ('1')
     dirichlet_flag = DIRICHLET_ONES_ON_DIAGONAL
     FLAbort('not supported in code')
  case ('2')
     dirichlet_flag = DIRICHLET_BIG_SPRING
     FLAbort('not supported in code')
  case('3') 
     dirichlet_flag = DIRICHLET_WEAK
  end select
  
  output_flag = 0
  call get_command_argument(5, value=buf, status=status)
  if(status>0) then
     call usage
     stop
  end if
  select case(buf)
  case('1')
     output_flag = 1
  end select

  ewrite(1,*) 'filename = ', filename
  ewrite(1,*) 'u continuity = ', u_continuity
  ewrite(1,*) 'h_order = ', h_order
  ewrite(1,*) 'dirichlet_flag = ',dirichlet_flag
  ewrite(1,*) 'output flag = ',output_flag
  
  ewrite(1,*) 'getting quadrature'

  call identify_mesh_file(filename, dim, loc, nnodes, nelements, &
       node_attributes)

  ewrite(1,*) 'dim = ', dim

  quad_degree = 6

  quad=make_quadrature(loc=loc, dimension=dim, degree=quad_degree)
  f_quad=make_quadrature(loc=loc-1, dimension=dim-1, degree=quad_degree-1)
  
  ewrite(1,*) 'Getting shape functions'

  ! Shape functions for positions (linear)
  X_shape=make_element_shape(loc=loc, dimension=dim, &
       degree=1, quad=quad)
  X_shape_f=make_element_shape(loc=loc-1, dimension=dim-1, &
       degree=1, quad=f_quad)

  ewrite(1,*) 'reading mesh'
  ewrite(1,*) 'loc = ',loc,'dim = ',dim

  positions=read_triangle_files(filename, X_shape)

  ewrite(1,*) node_count(positions)

  call add_faces(positions%mesh)

  ewrite(1,*) 'getting shapes'

  ! Shape functions for velocity and height
  u_shape=make_element_shape(loc=loc, &
       dimension=dim, degree=1, quad=quad)
  u_shape_f=make_element_shape(loc=loc-1, &
       dimension=dim-1, degree=1, quad=f_quad)
  h_shape=make_element_shape(loc=loc, &
       dimension=dim, degree=h_order, quad=quad)
  h_shape_f=make_element_shape(loc=loc-1, &
       dimension=dim-1, degree=h_order, quad=f_quad)

  !connectivity for velocity and height
  ewrite(1,*) 'Getting u connectivity'
  !u_mesh = make_mesh(positions%mesh,u_shape,0,'u_mesh')
  u_mesh = make_mesh(positions%mesh,u_shape,u_continuity,'u_mesh')
  call add_faces(u_mesh, model=positions%mesh)
  ewrite(1,*) 'Getting h connectivity'
  h_mesh = make_mesh(positions%mesh,h_shape,0,'h_mesh')
  call add_faces(h_mesh, model=positions%mesh)

  !fields for velocity and height
  ewrite(1,*) 'Allocating fields'

  call allocate(h,h_mesh,name='height')
  call allocate(RHS,h_mesh,name='RHS')
  call allocate(positions_h,dim,h_mesh,name='positions_h')  
  call allocate(positions_u,dim,u_mesh,name='positions_u')
  call allocate(u,dim,u_mesh,name='velocity')
  call allocate(h0,h%mesh,name='height0')
  call allocate(error,h%mesh,name='error')

  call remap_vector_field(positions,positions_h)
  call remap_vector_field(positions,positions_u)

  ewrite(1,*) node_count(u)
  ewrite(1,*) node_count(h)

  !get inverse mass matrix for u
  ewrite(1,*) 'allocating inverse mass matrix for u'

  call allocate(u1_inverse_mass,node_count(u),node_count(u))

  ewrite(1,*) 'getting inverse mass matrix for u'

  call get_dg_inverse_mass_matrix(u1_inverse_mass,u_mesh,positions)

  ewrite(1,*) 'allocating mass matrix for h'
  call allocate(h_mass,node_count(h),node_count(h))
  ewrite(1,*) 'assembling mass for h'
  call assemble_mass(positions,h,h_mass)

  ewrite(1,*) 'allocating CT'

  call allocate(C1T,node_count(h),node_count(u))
  call allocate(C2T,node_count(h),node_count(u))
  call allocate(C3T,node_count(h),node_count(u))

  ewrite(1,*) 'assembling CT'

  call assemble_CT(positions,u,h,C1T,C2T,C3T)

  ewrite(1,*) size(c1t,1)
  !if(dim<3) then
  !   call addto(C3T,1,1,0.0)
  !end if

  ewrite(1,*) 'assembling CMC'

  ewrite(1,*) 'static-ising u_inverse_mass'
  u_inverse_mass = dcsr2csr(u1_inverse_mass)

  ewrite(1,*) 'calling matmul_T mc1t'
  MC1T = matmul_T(C1T,u1_inverse_mass,check=.true.)

  ewrite(1,*) 'calling matmul_T mc2t'
  MC2T = matmul_T(C2T,u1_inverse_mass,check=.true.)
  call matcheck(MC1t,MC2T)
  if(dim==3) then
     ewrite(1,*) 'calling matmul_T mc3t'
     MC3T = matmul_T(C3T,u1_inverse_mass,check=.true.)
  end if
  ewrite(1,*) 'calling matmul_T cmc'
  CMC = matmul_T(C1T,MC1T,check=.true.)
  ewrite(1,*) 'calling matmul_T c2mc2t'
  C2MC2T = matmul_T(C2T,MC2T,check=.true.)
  if(dim==3) then
     ewrite(1,*) 'calling matmul_T c3mc3t'
     C3MC3T = matmul_T(C3T,MC3T,check=.true.)
  end if
  ewrite(1,*) 'calling addto'
  call addto(CMC,C2MC2T)
  if(dim==3) then
     call addto(CMC,C3MC3T)
  end if

  ewrite(1,*) 'staticising h_mass and cmc'

  h_mass_static = dcsr2csr(h_mass)
  cmc_static = dcsr2csr(cmc)

  ewrite(1,*) 'setting RHS'

  call set_RHS(positions,RHS,dirichlet_flag)

  if(dirichlet_flag==DIRICHLET_WEAK) then
     ewrite(1,*) 'fixing zero value'
     call addto(CMC_static,fix_pressure_value,fix_pressure_value,1.0e50)
  end if

  if(dirichlet_flag==DIRICHLET_NONE) then
     ewrite(2,*) ' Writing bcs'
     allocate(bc_marker(node_count(h)) )
     bc_marker = 0
     call get_bc_list(positions,h,bc_marker)
     if(.true.) then
        call lift_bcs(CMC_static,rhs,bc_marker)
        call lift_bcs(h_mass_static,rhs,bc_marker)
     else
        do i = 1, node_count(h)
           if(bc_marker(i)==1) then
              call addto(CMC_static,i,i,1.0e50)
           end if
        end do
     end if
  end if

  select case(output_flag)
  case(0)
     ewrite(1,*) 'calling petsc'
     
     call petsc_solve(h%val,CMC_Static, rhs%val, &
          MATRIX_SYMMETRIC, 100000, 1.0e-12, &
          .true.,.true.)
     
     ewrite(1,*) 'checking error'
     
     call check_error(positions,h,h_mass_static,h0,error,dirichlet_flag)
     
     ewrite(1,*) 'dumping data'
     
     call vtk_write_fields('solve_CMC_out', index=0, position=positions, &
          model=positions%mesh, sfields=(/RHS,h,h0,error/))
  case(1)
     call matrix2file('hmass_solve.dat',h_mass_static)
     call matrix2file('cmc_solve.dat',cmc_static)
  case default
     FLAbort('no such output option')
  end select

  ewrite(1,*) 'U dof =', node_count(u)
  ewrite(1,*) 'H dof =', node_count(h)

  ewrite(1,*) 'END program solve_CMC'

contains

  subroutine matcheck(M1,M2)
    type(dynamic_csr_matrix), intent(in) :: M1,M2
    integer ::i,j

    do i = 1, size(M1,1)
       if(size(M1%colm(i)%ptr).ne.size(M2%colm(i)%ptr)) then
          print *, size(M1%colm(i)%ptr), size(M2%colm(i)%ptr)
          print *, M1%colm(i)%ptr
          print *, M2%colm(i)%ptr
          print *,i,j
          FLAbort('wrong size colm')
       end if
       if(any(M1%colm(i)%ptr.ne.M2%colm(i)%ptr)) then
          FLAbort('wrong colm values')
       end if

       !do j = 1, size(M1%colm(i)%ptr)
       !   if(abs(M1%val(i)%ptr(j)-M2%val(i)%ptr(j))>1.0e-10) then
       !      print *, M1%val(i)%ptr(j), M2%val(i)%ptr(j)
       !      FLAbort('wrong values')
       !   end if
       !end do
    end do
    if(size(M2,1)>1) then
       do i = 1, size(M2,1)
          if(size(M2%colm(i)%ptr)>1) then
             do j = 2, size(M2%colm(i)%ptr)
                if(M2%colm(i)%ptr(j).le.(M2%colm(i)%ptr(j-1))) then
                   FLAbort('bad ordering in M2')
                end if
             end do
          end if
       end do
    end if
    
  end subroutine matcheck
  
  subroutine get_bc_list(positions,h,bc_marker)
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: h
    integer, dimension(:), intent(out) :: bc_marker

    !locals
    integer, dimension(:), pointer :: neigh
    integer :: ni,ele,i,nod, ele_2, face

    bc_marker = 0

    do ele = 1, element_count(h)
    
       !local variables
       
       neigh=>ele_neigh(U, ele)
       
       neighbourloop: do ni=1,size(neigh)
          
          ele_2=neigh(ni)
          
          face=ele_face(h, ele, ele_2)
          
          if (ele_2>0) cycle
          
          call get_bc_list_face(positions,h,bc_marker,face)
       
       end do neighbourloop
    end do

  end subroutine get_bc_list
  
  subroutine get_bc_list_face(positions,h,bc_marker,face)
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: h
    integer, intent(in) :: face
    integer, dimension(:), intent(inout) :: bc_marker
    !    
    integer, dimension(face_loc(H,face)) :: h_face
    
    h_face=face_global_nodes(H, face)
    bc_marker(h_face) = 1

  end subroutine get_bc_list_face

  subroutine assemble_CT(positions,u,h,C1T,C2T,C3T)
    type(vector_field), intent(in) :: positions, u
    type(scalar_field), intent(in) :: h
    type(dynamic_csr_matrix), intent(inout) :: C1T,C2T,C3T

    !locals
    integer :: ele

    do ele = 1, element_count(u)
       call assemble_CT_elemental(ele,positions,u,h,C1T,C2T,C3T)
    end do

  end subroutine assemble_CT
  
  subroutine assemble_CT_elemental(ele,positions,u,h,C1T,C2T,C3T)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: u
    type(scalar_field), intent(in) :: h
    type(dynamic_csr_matrix), intent(inout) :: C1T,C2T,C3T

    ! Locations of nodes.
    real, dimension(positions%dim,ele_loc(positions,ele)) :: X_ele
    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,ele)) :: detwei    
    ! Derivatives of shape function:
    real, dimension(ele_loc(h,ele), &
         ele_ngi(h,ele), positions%dim) :: dshape_h
    real, dimension(ele_loc(u,ele), &
         ele_ngi(u,ele), positions%dim) :: dshape_u
    ! Node numbers of field element.
    integer, dimension(:), pointer :: ele_u, ele_h
    ! Shape functions.
    type(element_type), pointer :: shape_u, shape_h, shape_X
    ! gradient matrix
    real, dimension(positions%dim,ele_loc(h,ele),ele_loc(u,ele)) :: grad_mat
    integer, dimension(:), pointer :: neigh
    integer :: ele_2, ni, face

    ele_u=>ele_nodes(u, ele)
    shape_u=>ele_shape(u, ele)
    ele_h=>ele_nodes(h, ele)
    shape_h=>ele_shape(h, ele)
    shape_X=>ele_shape(positions, ele)

    ! Locations of local vertices.
    X_ele=ele_val(positions, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(X_ele, shape_X, m=shape_h, &
         dm_t=dshape_h, detwei=detwei)

    grad_mat = -dshape_shape(dshape_h,shape_u,detwei)

    call addto(C1T,ele_h,ele_u,grad_mat(1,:,:))
    call addto(C2T,ele_h,ele_u,grad_mat(2,:,:))
    if(positions%dim==3) then
       call addto(C3T,ele_h,ele_u,grad_mat(3,:,:))
    end if

    !face integrals

    if(.false.) then
       if(dirichlet_flag==DIRICHLET_NONE) then
          neigh=>ele_neigh(U, ele)

          neighbourloop: do ni=1,size(neigh)

             !------------------------------------------------------
             ! Find the relevant faces.
             !------------------------------------------------------

             ! These finding routines are outside the inner loop 
             ! so as to allow
             ! for local stack variables of the right size in
             ! construct_momentum_interface_dg.
             ele_2=neigh(ni)

             face=ele_face(U, ele, ele_2)

             if (ele_2>0) cycle

             n_boundary_elements = n_boundary_elements + 1

             call assemble_CT_interface(C1T,C2T,C3T, &
                  & ele,face,u,h,positions)

          end do neighbourloop
       end if
    end if

  end subroutine assemble_CT_elemental

  subroutine assemble_CT_interface(C1T,C2T,C3T,ele,face,u,h,positions)
    type(dynamic_csr_matrix), intent(inout) :: C1T, C2T,C3T
    integer, intent(in) :: ele,face
    type(vector_field), intent(in) :: positions, u
    type(scalar_field), intent(in) :: h

    !local stuff
    real, dimension(face_ngi(U,face)) :: detwei
    real, dimension(U%dim, face_ngi(U, face)) :: normal
    type(element_type), pointer :: u_shape, h_shape
    integer, dimension(face_loc(U,face)) :: u_face
    integer, dimension(face_loc(H,face)) :: h_face
    real, dimension(U%dim, &
         & face_loc(h,face),face_loc(U,face)) :: mnCT
    !stuff for debugging
    real, dimension(positions%dim, face_loc(positions,face)) :: pos_face_val
    real, dimension(U%dim, face_loc(U,face)) :: U_face_val
    real, dimension(U%dim, face_loc(h,face)) :: H_face_val
    integer :: I

    u_shape=>face_shape(U, face)
    u_face=face_global_nodes(U, face)
    h_face=face_global_nodes(H, face)
    h_shape=>face_shape(H, face)

    call transform_facet_to_physical(positions, face, &
         &                          detwei_f=detwei,&
         &                          normal=normal) 

    mnCT = shape_shape_vector(h_shape, U_shape, detwei, normal)

    call addto(C1T,h_face,u_face, mnCT(1,:,:))
    call addto(C2T,h_face,u_face, mnCT(2,:,:))
    if(positions%dim==3) then
       call addto(C3T,h_face,u_face, mnCT(3,:,:))
    end if

  end subroutine assemble_CT_interface

  subroutine assemble_mass(positions,h,mass)
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: h
    type(dynamic_csr_matrix), intent(inout) :: mass

    !locals
    integer :: ele

    do ele = 1, element_count(h)
       call assemble_mass_elemental(ele,positions,h,mass)
    end do

  end subroutine assemble_mass

  subroutine assemble_mass_elemental(ele,positions,h,mass)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: h
    type(dynamic_csr_matrix), intent(inout) :: mass

    ! Locations of nodes.
    real, dimension(positions%dim,ele_loc(positions,ele)) :: X_ele
    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,ele)) :: detwei    
    ! Node numbers of field element.
    integer, dimension(:), pointer :: ele_h
    ! Shape functions.
    type(element_type), pointer :: shape_h, shape_X
    ! local mass matrix
    real, dimension(ele_loc(h,ele),ele_loc(h,ele)) :: mass_mat

    ele_h=>ele_nodes(h, ele)
    shape_h=>ele_shape(h, ele)
    shape_X=>ele_shape(positions, ele)

    ! Locations of local vertices.
    X_ele=ele_val(positions, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(X_ele, shape_X, detwei=detwei)

    mass_mat = shape_shape(shape_h,shape_h,detwei)

    call addto(mass,ele_h,ele_h,mass_mat)

  end subroutine assemble_mass_elemental

  subroutine set_RHS(positions,RHS,dirichlet_flag)
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(inout) :: RHS
    integer, intent(in) :: dirichlet_flag
    !locals
    integer :: ele

    call zero(rhs)

    do ele = 1,element_count(RHS)
       call set_RHS_elemental(ele,positions,RHS,dirichlet_flag)
    end do

    ewrite(1,*) 'rhs integral', sum(RHS%val)

  end subroutine set_RHS

  subroutine set_RHS_elemental(ele,positions,RHS,dirichlet_flag)
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(inout) :: RHS
    integer, intent(in) :: ele, dirichlet_flag

    ! Locations of nodes.
    real, dimension(positions%dim,ele_loc(positions,ele)) :: X_ele
    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,ele)) :: detwei    
    ! Node numbers of field element.
    integer, dimension(:), pointer :: ele_h
    ! Shape functions.
    type(element_type), pointer :: shape_h, shape_X
    ! local mass matrix
    real, dimension(ele_loc(rhs,ele),ele_loc(rhs,ele)) :: mass_mat
    real, dimension(positions%dim,ele_ngi(positions,ele)) :: X_quad

    ele_h=>ele_nodes(rhs, ele)
    shape_h=>ele_shape(rhs, ele)
    shape_X=>ele_shape(positions, ele)

    ! Locations of local vertices.
    X_ele=ele_val(positions, ele)
    ! Locations of quadrature points.
    X_quad=ele_val_at_quad(positions, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(X_ele, shape_X, detwei=detwei)
    
    call addto(rhs,ele_h,shape_rhs(shape_h, &
         detwei*rhs_fun(X_quad,dirichlet_flag)))

  end subroutine set_RHS_elemental

  subroutine check_error(positions,h,h_mass,h0,error,dirichlet_flag)
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(inout) :: h,h0,error
    type(csr_matrix), intent(in) :: h_mass
    integer, intent(in) :: dirichlet_flag
    !locals
    integer :: ele
    type(vector_field) :: X_h !positions mesh 
    ewrite(1,*) 'THIS WILL ONLY WORK IF MESH IS IN A 1 by 1 SQUARE OR CUBE'

    call allocate(X_h,positions%dim,h%mesh,name='X_h')

    call remap_vector_field(positions,X_h)
    
    select case(X_h%dim)
    case(2)
       select case(dirichlet_flag)
       case (0)
          h0%val = sin(2*pi*X_h%val(1)%ptr)*sin(2*pi*X_h%val(2)%ptr)
       case (3)
          h0%val = cos(2*pi*X_h%val(1)%ptr)*cos(2*pi*X_h%val(2)%ptr)
          h%val = h%val + h0%val(fix_pressure_value)
       case default
          FLAbort('Boundary condition option not supported')
       end select
    case(3)
       select case(dirichlet_flag)
       case (0)
          h0%val = sin(2*pi*X_h%val(1)%ptr)*sin(2*pi*X_h%val(2)%ptr)* &
               sin(2*pi*X_h%val(3)%ptr)
       case (3)
          h0%val = cos(2*pi*X_h%val(1)%ptr)*cos(2*pi*X_h%val(2)%ptr)* &
               cos(2*pi*X_h%val(3)%ptr)
          h%val = h%val + h0%val(fix_pressure_value)
       case default
          FLAbort('Boundary condition option not supported')
       end select
    case default
       FLAbort('dimension not supported')
    end select

    error%val = h0%val-h%val

    ewrite(1,*) 'maxval error =',maxval(error%val)    
    ewrite(1,*) 'minval error =',minval(error%val)
    ewrite(1,*) 'maxval h =',maxval(h%val)    
    ewrite(1,*) 'minval h =',minval(h%val)
    ewrite(1,*) 'maxval h0 =',maxval(h0%val)    
    ewrite(1,*) 'minval h0 =',minval(h0%val)

    ewrite(1,*) 'error is', maxval(abs(error%val))

  end subroutine check_error

  !function to return right-hand side 
  function rhs_fun(X,dirichlet_flag)
    real, dimension(:,:), intent(in) :: X
    real, dimension(size(X,2)) :: rhs_fun
    integer, intent(in) :: dirichlet_flag

    select case(size(X,1))
    case(2)
       select case(dirichlet_flag)
       case (DIRICHLET_NONE)
          rhs_fun = 8.0*pi*pi*sin(2*pi*X(1,:))*sin(2*pi*X(2,:))
       case (DIRICHLET_WEAK)
          rhs_fun = 8.0*pi*pi*cos(2*pi*X(1,:))*cos(2*pi*X(2,:))
       case default
          FLAbort('Bad boundary condition option')
       end select
    case(3)
       select case(dirichlet_flag)
       case (DIRICHLET_NONE)
          rhs_fun = 12.0*pi*pi*sin(2*pi*X(1,:))*sin(2*pi*X(2,:))* &
               sin(2*pi*X(3,:))
       case (DIRICHLET_WEAK)
          rhs_fun = 12.0*pi*pi*cos(2*pi*X(1,:))*cos(2*pi*X(2,:))* &
               cos(2*pi*X(3,:))
       case default
          FLAbort('Bad boundary condition option')
       end select
    case default
       FLAbort('bad stuff')
    end select
    
  end function rhs_fun

  subroutine usage
    
    write (0,*) "usage: solve_CMC <triangle_file_name> <u continuity> <h order> <dirichlet flag> <output_flag>"
    
  end subroutine usage

  subroutine single_element(u_shape,u_shape_f,h_shape,h_shape_f)
    type(element_type), intent(in) :: u_shape,u_shape_f,h_shape,h_shape_f

    FLAbort('Ending after single_element')
  end subroutine single_element

  subroutine lift_bcs(CMC,rhs,bc_marker)
    type(csr_matrix), intent(inout) :: CMC
    type(scalar_field), intent(inout) :: rhs
    integer, dimension(:), intent(in) :: bc_marker
    !
    integer :: i,jrow,j
    integer, dimension(:), pointer :: row
    real, dimension(:), pointer :: row_val

    !set RHS to bc value (zero bcs)
    do i = 1, node_count(rhs)
       if(bc_marker(i)==1) rhs%val(i) = 0.
    end do

    !move columns to RHS
    !(We don't need to do this for zero bcs)

    !zero rows
    do i = 1, node_count(rhs)
       if(bc_marker(i)==1) then
          row_val => row_val_ptr(CMC,i)
          row_val = 0.
       end if
    end do

    !zero columns
    do i = 1, node_count(rhs)
       row => row_m_ptr(CMC,i)
       if(any(bc_marker(row)==1)) then
          do jrow = 1, size(row)
             j = row(jrow)
             if(bc_marker(j)==1) then
                call set(CMC,i,j,0.)
             end if
          end do
       end if
    end do

    !put ones on diagonals
    do i = 1, node_count(rhs)
       if(bc_marker(i)==1) then
          call set(CMC,i,i,1000000.0)
       end if
    end do

  end subroutine lift_bcs

end program solve_CMC
