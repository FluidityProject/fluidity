#include "fdebug.h"

program main

  !a wrapper for P1nc solver for helmholtz equation

  !coded by Colin Cotter and David Ham May-June 2006
  !for last modification see CVS

  !links against libdfluidity

  use elements
  use sparse_tools
  use quadrature
  use global_numbering
  use shape_functions
  use global_parameters_gallopede
  use adjacency_lists
  use transform_elements
  use dgtools
  use text_io
  use helmholtz

  implicit none

  !local variables

  real, allocatable, dimension(:)::X,Y !vertex coordinates
  real, allocatable::m1(:),m2(:) !momentum values
  real, allocatable::u1(:),u2(:) !velocity values
  integer, allocatable::EVList_h(:) !element-vertex list for height
  integer, allocatable::EVList_u(:) !element-vertex list for vels
  real, allocatable:: u1_rhs(:), u2_rhs(:) !rhs for mom-vel equation
  !type(csr_matrix):: Mass ! Mass Matrix
  !type(block_csr_matrix)::Helm !Helmholtz matrix
  type(csr_matrix)::Helm !Helmholtz matrix
  type(element_type) :: nh, nh_f   ! Reference height element
  type(element_type) :: nu  ! Reference velocity element
  type(quadrature_type) :: g, g_f   ! Gauss quadrature
  type(csr_matrix)::EEList, B_SEG_list ! Element-Element list, B_SEG list
  integer, dimension(:), allocatable :: b_seg_nh_lno
  integer :: i, globi, globj
  logical :: nc
  real :: alpha1

  alpha1 = 10.0/32.0

  !load mesh data

  !nc = .false.
  nc = .true.

  MSG("Loading mesh data");
  call Get_Mesh(X,Y,EVList_h,EVList_u)

  MSG("allocating memory");
  if(nc) then
     allocate( u1(n_vels) )
     allocate( u2(n_vels) )
     allocate( u1_rhs(n_vels) )
     allocate( u2_rhs(n_vels) )
  else
     allocate( u1(n_verts) )
     allocate( u2(n_verts) )
     allocate( u1_rhs(n_verts) )
     allocate( u2_rhs(n_verts) )
  end if

  !load momentum data
  if(nc) then
     MSG("Loading source data");
     allocate( m1(n_vels) )
     allocate( m2(n_vels) )
     call read_field(m1,N_Vels,'u1_initial.dat',14)
     call read_field(m2,N_Vels,'u2_initial.dat',14)
  else
     allocate( m1(n_verts) )
     allocate( m2(n_verts) )
  end if

  !load height data
  !call Get_Height(D)

  !get quadrature
  !the 2d elements
  MSG("Making 2D quadrature");
  g=make_quadrature(loc=3,dimension=2,degree=3)
  MSG("Making 1D quadrature");
  g_f=make_quadrature(loc=2,dimension=1,degree=2)

  !get_basis functions
  !2D

  MSG("Making 2D P1 elements");
  nh=make_element_shape(loc=3, dimension=2, degree=1, quad=g)
  MSG("Making 2D P1nc elements");
  nu=make_element_shape(loc=3, dimension=2, degree=1, quad=g, &
       type=ELEMENT_NONCONFORMING)

  !1D -- continuous basis functions
  MSG("Making 1D P1 elements");
  nh_f=make_element_shape(loc=2, dimension=1, degree=1, quad=g_f)

  !get Helmholtz matrix structure
  MSG("Making matrix structure");
  call MakeLists(N_Verts,N_Elements, Nloc, EVList_h, .false., EEList=EEList)
  CHECK(entries(EEList))

  MSG("Allocating b_seg_nh_lno");
  !stores local node numbers for h for face numbers
  allocate(b_seg_nh_lno(entries(EEList)*nh_f%loc))
  CHECK(size(b_seg_nh_lno))
  MSG("Making boundary_seg numbering");
  call make_boundary_seg_numbering_nc(b_seg_list, b_seg_nh_lno, & 
       EElist = EEList, nonod = N_Verts, ndglno = EVList_h, &
       Ele_n = nh, b_seg_n = nh_f)

  ! Set up matrix for eliptic equation.

  if(nc) then
     MSG("Calling posinm_dg_nc");
     call POSINM_DG_nc(Helm,N_Elements, N_Vels, nu%loc, EVList_u,&
          N_Vels, nu%loc, EVList_u, eelist)
  else
     MSG("calling posinm")
     call POSINM(Helm,N_Elements, N_Verts, nh%loc, EVList_h,&
          N_Verts, nh%loc, EVList_h)
  end if

  !Helm = block_clone(Mass, (/2,2/))
  
  !assemble momentum-vel equation
  MSG("Assembling Helmholtz equation");
  call Assembl_helmholtz_eqn(Helm,X,Y,EVList_h,EVList_u, &
       m1,m2,u1_rhs,u2_rhs, &
       nu,nh,nh_f, &
       b_seg_list,b_seg_nh_lno,nc,alpha1)

  CHECK(sum(u1_rhs*u1_rhs))
  CHECK(sum(Helm%val*Helm%val))

  !solve momentum-vel equation (using cg)
  MSG("Solving Helmholtz equation");
  call solve_helmholtz_eqn(u1, Helm, u1_rhs, 1.0E-12, 1000)

  CHECK(n_vels)
  MSG("Done")

  !output velocity solution
  call write_vel(u1)

  deallocate( u1 )
  deallocate( u2 )
  deallocate( u1_rhs )
  deallocate( u2_rhs )
  deallocate( b_seg_nh_lno )
  call deallocate( EEList )
  call deallocate( B_seg_List )
  call deallocate( Helm )

contains

end program main
