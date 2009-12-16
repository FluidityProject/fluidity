#include "fdebug.h"

program main

  !solver for the GN operator

  !coded by Colin Cotter and David Ham May-Sep 2006
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
  use data_structures
  use fetools
  use vector_tools
  use gallopede_solvers
  use vtk_io
  use mesh_tools
  use GN_operator
  use gallopede_solvers

  implicit none

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"

  !local variables

  real, allocatable, target::X(:),Y(:) !vertex coordinates
  real, allocatable, dimension(:), target ::u,m
  !velocity vector storing all velocities
  real, pointer, dimension(:) :: u1, u2, m1, m2
  !individual components
  real, allocatable::Din(:), D(:) !height field
  !integer, allocatable, target::EVList_h(:) !element-vertex list for height
  integer, allocatable, target::EVList_X(:) !element-vertex list for points
  integer, allocatable, target::EVList_u(:) !element-vertex list for vels
  integer, allocatable, target::EVList_h(:) !element-vertex list for height
  type(csr_matrix),target:: Mass_h,Mass_u !mass matrix
  type(block_csr_matrix) :: mommat !differentiation matrix
  type(csr_matrix) :: mom11, mom12, mom21, mom22 !components of diff mat
  type(element_type) :: nu,nu_f,nh,nh_f
  type(quadrature_type) :: g, g_f   ! Gauss quadrature
  type(csr_matrix),target ::EEList, bdy_list, bdy_ordering
  ! Element-Element list, B_SEG list
  type(dg_mesh) :: mesh
  type(rotations) :: rots
  integer, dimension(:), allocatable, target :: bdy_nu_lno, bdy_nh_lno
  integer :: i,globi,globj

  PetscErrorCode :: ierr
  KSPType :: ksp_type
  PCType :: pc_type

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr);
  !==============================================
  !initial conditions

  !load mesh data
  MSG("Loading mesh data")
  call Get_Mesh(X,Y,EVList_X,EVList_u)

  MSG("allocating memory");
  allocate( u(N_vels*2) )
  u1 => u(1:N_vels)
  u2 => u(N_vels+1:2*N_vels)
  allocate( m(N_vels*2) )
  m1 => m(1:N_vels)
  m2 => m(N_vels+1:2*N_vels)

  allocate( Din(n_verts) )

  MSG("Loading initial conditions");

  call read_field(dr=m1,filename='m1.dat')
  call read_field(dr=m2,filename='m2.dat')
  call read_field(dr=u1,filename='u1.dat')
  call read_field(dr=u2,filename='u2.dat')
  call read_field(dr=Din,filename='D_initial.dat')

  !==============================================
  !Initialisation of mesh data

  !get quadrature
  !the 2d elements
  MSG("Making 2D quadrature");
  g=make_quadrature(loc=3,dimension=2,degree=4)
  MSG("Making 1D quadrature");
  g_f=make_quadrature(loc=2,dimension=1,degree=3)

  !get_basis functions
  !2D

  MSG("Making 2D P1 elements");
  nu=make_element_shape(loc=3, dimension=2, degree=1, quad=g)
  MSG("Making elements for div u");
  nh=make_element_shape(loc=3, dimension=2, degree=2, quad=g)

  !1D -- continuous basis functions
  MSG("Making 1D P1 elements");
  nu_f=make_element_shape(loc=2, dimension=1, degree=1, quad=g_f)
  MSG("Making elements for div u");
  nh_f=make_element_shape(loc=2, dimension=1, degree=2, quad=g_f)

  !get Helmholtz matrix structure
  MSG("Making matrix structure");

  CHECK(N_verts)
  CHECK(N_elements)
  CHECK(nloc)
  call MakeLists(N_Verts,N_Elements, Nloc, EVList_X, &
       .false., EEList=EEList)
  CHECK(entries(EEList))

  allocate( D(n_dens) )
  allocate( EVList_h(6*N_elements) )
  call get_D_mesh(Din,D,EVList_h,EVList_X,EEList)
  deallocate( Din )

  allocate(bdy_nu_lno(entries(EEList)*nu_f%loc))
  allocate(bdy_nh_lno(entries(EEList)*nh_f%loc))

  MSG("Making bdy numbering");
  CHECK(nu%numbering%boundaries)
  call make_boundary_numbering(bdy_list, bdy_nu_lno, &
       boundary_m_lno = bdy_nh_lno, &
       EElist = EEList, xnonod = N_Verts, xndglno = EVlist_X, &
       ele_n = nu, boundary_n = nu_f, ele_m = nh, boundary_m = nh_f)

  ! Set up Mass matrices
  call POSINM_DG(Mass_u, N_elements, N_vels, nu%loc, evlist_u,&
       n_vels, nu%loc, evlist_u, bdy_list,nu_f%loc,bdy_nu_lno)

  CHECK(size(mass_u))
  CHECK(size(mass_u,1))
  call zero(mass_u)

  !call getmass(mass_h,nh,evlist_h,nh,evlist_h)
  !call getmass(mass_u,nu,evlist_u,nu,evlist_u)

  MSG("setting up nc mesh data")
  mesh%n_elements = n_elements
  mesh%n_verts = n_verts
  mesh%n_vels = n_vels
  mesh%X => X
  mesh%Y => Y
  mesh%Evlist_u => Evlist_u
  mesh%Evlist_X => Evlist_X
  mesh%Evlist_h => Evlist_h
  mesh%EElist => EElist
  mesh%bdy_ordering => bdy_ordering
  mesh%mass_u => mass_u
  !mesh%mass_h => mass_h
  mesh%bdy_list => bdy_list
  mesh%bdy_nu_lno => bdy_nu_lno
  mesh%bdy_nh_lno => bdy_nh_lno
  mesh%nu = nu
  mesh%nu_f = nu_f
  mesh%nh = nh
  mesh%nh_f = nh_f

  call zero(mesh%bdy_ordering)
  !call get_bdy_ordering(mesh)

  MSG("cloning matrix")
  mommat = block_clone(mass_u,(/2,2/) )
  MSG("Getting blocks")
  mom11 = block(mommat,1,1)
  mom12 = block(mommat,1,2)
  mom21 = block(mommat,2,1)
  mom22 = block(mommat,2,2)
  MSG("zeroing matrix")
  call zero(mom11)
  call zero(mom12)
  call zero(mom21)
  call zero(mom22)

  call Get_Parameters()

  call get_vels(Mesh, D, u, m)

!!$  rhs = 0.
!!$
!!$  MSG("call assemble_GN_operator")
!!$  !assemble GN operator and rhs
!!$  call assemble_GN_D_operator( D, mesh, &
!!$       mom11=mom11,mom12=mom12,mom21=mom21,mom22=mom22, &
!!$       rhs1=rhs1,rhs2=rhs2, &
!!$       m1=m1,m2=m2)
!!$
!!$  MSG("solving GN operator equation -- first pass")
!!$  !solve GN equation
!!$  ksp_type = KSPCG
!!$  pc_type = PCICC
!!$  u = 0.
!!$  call gallopede_block_solve(u, mommat, rhs, ksp_type, pc_type, 1.0e-5, 500)
!!$
!!$  call add_penalty_term(mom11,mom12,mom21,mom22,u1,u2,mesh)
!!$
!!$  MSG("solving GN operator equation -- second pass")
!!$  !solve GN equation
!!$  ksp_type = KSPCG
!!$  !pc_type = PCSOR
!!$  pc_type = PCICC
!!$  !u = 0.
!!$
!!$  call gallopede_block_solve(u, mommat, rhs, ksp_type, pc_type, 1.0e-8, 500)

  !data output
  call dump_field(u1,N_vels,'u1_',1,1,9)
  call dump_field(u2,N_vels,'u2_',1,1,9)
  !call dump_field_vtk(u1,N_vels,EVList_u,EVLIST_X,X,Y,N_elements,N_verts,3, &
  !     'u1_',1,1,9)
  !call dump_field_vtk(u2,N_vels,EVList_u,EVLIST_X,X,Y,N_elements,N_verts,3, &
  !     'u2_',1,1,9)

  deallocate( u )
  deallocate( m )
  deallocate( D )
  deallocate( bdy_nu_lno )
  deallocate( bdy_nh_lno )
  call deallocate( EEList )
  call deallocate( bdy_List )

end program main
