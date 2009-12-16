#include "fdebug.h"

subroutine GN_operator_main()

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

  type(quadrature_type) :: g, g_f
  real, allocatable, dimension(:), target ::u,m
  !velocity vector storing all velocities
  real, pointer, dimension(:) :: u1, u2, m1, m2
  !individual components
  real, allocatable::Din(:), D(:) !height field
  type(block_csr_matrix) :: mommat !differentiation matrix
  type(csr_matrix) :: mom11, mom12, mom21, mom22 !components of diff mat
  ! Element-Element list, B_SEG list
  type(dg_mesh) :: mesh
  type(bc_info) :: bcs,quadbcs

  PetscErrorCode :: ierr
  KSPType :: ksp_type
  PCType :: pc_type

  !load mesh data
  ewrite(2,*)("Loading mesh data")
  call Get_Mesh(mesh%X,mesh%Y,mesh%EVList_X,mesh%EVList_m,bcs)
  mesh%EVList_u => mesh%EVList_X
  N_vels = N_verts
  call get_bc_marker(bcs)

  ewrite(2,*)("allocating memory");
  allocate( u(N_vels*2) )
  u1 => u(1:N_vels)
  u2 => u(N_vels+1:2*N_vels)
  allocate( m(N_moms*2) )
  m1 => m(1:N_moms)
  m2 => m(N_moms+1:2*N_moms)

  allocate( Din(n_verts) )

  ewrite(2,*)("Loading initial conditions");

  call read_field(dr=m1,filename='m1.dat')
  call read_field(dr=m2,filename='m2.dat')
  call read_field(dr=Din,filename='D_initial.dat')

  !==============================================
  !Initialisation of mesh data
  !get quadrature
  !the 2d elements
  ewrite(2,*)("Making 2D quadrature");
  g=make_quadrature(loc=3,dimension=2,degree=4)
  ewrite(2,*)("Making 1D quadrature");
  g_f=make_quadrature(loc=2,dimension=1,degree=3)

  !get_basis functions
  !2D
  mesh%nm=make_element_shape(loc=3, dimension=2, degree=1, quad=g)
  mesh%nu=make_element_shape(loc=3, dimension=2, degree=1, quad=g)
  mesh%nh=make_element_shape(loc=3, dimension=2, degree=2, quad=g)

  call POSINM(mesh%mass_u,N_elements,N_vels,mesh%nu%loc,mesh%Evlist_u,&
       N_vels, mesh%nu%loc, mesh%EVList_u)
  call MakeLists(N_Verts,N_Elements,Nloc,mesh%EVList_X, &
       .false., EEList=mesh%EEList)
  allocate( mesh%EVList_h(6*N_elements) )
  call get_quad_mesh(mesh%EVList_h,mesh%EVList_X,mesh%EEList,bcs,quadbcs, &
       mesh%X,mesh%Y,N_dens)
  allocate( D(N_dens) )
  call project_to(D,Din,mesh%nh,mesh%nm,mesh%nm, &
       mesh%EVList_h,mesh%EVList_X,mesh%EVList_X, &
       mesh%Mass_h,.true.,.true.,.true.)
  deallocate( Din )

  !call Get_Parameters()

  call get_vels_GN(Mesh, u, m, D, bcs)

  !data output
  call dump_data_vtk(u1,u2,D,mesh%EVList_U,mesh%EVList_h,mesh%EVList_X, &
       mesh%X,mesh%Y,N_Elements,'gn_op_out',1,1,15)

  deallocate( u )
  deallocate( m )
  deallocate( D )
  !deallocate( bdy_nu_lno )
  !deallocate( bdy_nh_lno )
  call deallocate( mesh%EEList )
  deallocate( mesh%EVList_X )
  deallocate( mesh%EVList_h )
  deallocate( mesh%EVList_m )

end subroutine GN_operator_main
