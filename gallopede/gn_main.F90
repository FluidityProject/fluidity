#include "fdebug.h"

subroutine gn_main()

  !this stage: a P1nc-P1 solver for Green-Naghdi

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
  use vtk_io
  use density_equation
  use gn_momentum_equation
  use gn_operator
  use data_structures
  use gallopede_solvers
  use mesh_tools
  use fldebug

  implicit none

  !local variables

  real, allocatable, dimension(:), target :: u,m
  !velocity vector storing all velocities
  real, pointer, dimension(:) :: u1, u2, m1, m2 !individual components
  real, allocatable, dimension(:) ::D, Din !height field
  type(quadrature_type) :: g, g_f   ! Gauss quadrature
  type(dg_mesh) :: mesh
  type(bc_info), target :: bcs,quadbcs
  type(bc_info), pointer :: usebcs
  integer, parameter :: u_deg = 1
  integer :: i, pad

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"

  !==============================================
  !initial conditions

  !load mesh data
  ewrite(2,*)("Loading mesh data");
  call Get_Mesh(mesh%X,mesh%Y,mesh%EVList_X,mesh%EVList_m,bcs)
  call get_bc_marker(bcs,N_verts)

  ewrite(2,*)("allocating memory");
  select case (u_deg)
  case (1)
     N_vels = N_verts
     usebcs => bcs
  case (2)
     N_vels = N_dens
     usebcs => quadbcs
  end select

  allocate( u(N_vels*2) )
  u1 => u(1:N_vels)
  u2 => u(N_vels+1:2*N_vels)
  allocate( m(N_moms*2) )
  m1 => m(1:N_moms)
  m2 => m(N_moms+1:2*N_moms)
  
  allocate( Din(n_verts) )

  ewrite(2,*)("Loading initial conditions");
  call read_field(dr=m1,filename='m1_initial.dat')
  call read_field(dr=m2,filename='m2_initial.dat')
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

  ewrite(2,*)("Making 2D P2 elements");
  mesh%nh=make_element_shape(loc=3, dimension=2, degree=2, quad=g)
  ewrite(2,*)("Making 2D P1 elements");
  mesh%nu=make_element_shape(loc=3, dimension=2, degree=u_deg, quad=g)
  mesh%nm=make_element_shape(loc=3, dimension=2, degree=1, quad=g)
  mesh%nx=make_element_shape(loc=3, dimension=2, degree=1, quad=g)

  !1D -- continuous basis functions
  ewrite(2,*)("Making 1D P1 elements");
  mesh%nu_f=make_element_shape(loc=2, dimension=1, degree=u_deg, quad=g_f)
  mesh%nm_f=make_element_shape(loc=2, dimension=1, degree=1, quad=g_f)
  ewrite(2,*)("Making 1D P2 elements");
  mesh%nh_f=make_element_shape(loc=2, dimension=1, degree=2, quad=g_f)

  ewrite(2,*)("Making matrix structure");
  call MakeLists(N_Verts,N_Elements, Nloc, &
       mesh%EVList_X, .false., EEList=mesh%EEList)
  call MakeLists(N_Verts,N_Elements,Nloc,mesh%EVList_X, &
       .false., EEList=mesh%EEList)
  allocate( mesh%EVList_h(6*N_elements) )
  call get_quad_mesh(mesh%EVList_h,mesh%EVList_X,mesh%EEList,mesh%X,mesh%Y,N_dens,bcs,quadbcs)
  allocate( D(N_dens) )
  ewrite(2,*) maxval( abs(Din) )
  call project_to(D,Din,mesh%nh,mesh%nm,mesh%nm, &
       mesh%EVList_h,mesh%EVList_X,mesh%EVList_X, &
       mesh%Mass_h,.true.,.true.,.true.)
  deallocate( Din )
  select case (u_deg)
  case(1)
     mesh%EVList_u => mesh%EVList_X
     call POSINM(mesh%Mass_U,N_elements,N_vels, &
          mesh%nu%loc,mesh%Evlist_X,&
          N_vels,mesh%nu%loc,mesh%EVList_X)
  case(2)
     mesh%EVList_u => mesh%EVList_h
     mesh%Mass_u = clone(mesh%Mass_h)
  end select

  ewrite(2,*)("Allocating bdy_n_lno");
  !stores local node numbers for h for face numbers
  allocate(mesh%bdy_nm_lno(entries(mesh%EEList)*mesh%nm_f%loc))
  allocate(mesh%bdy_nh_lno(entries(mesh%EEList)*mesh%nh_f%loc))
  select case (u_deg)
  case (1)
     mesh%bdy_nu_lno => mesh%bdy_nm_lno
  case (2)
     mesh%bdy_nu_lno => mesh%bdy_nh_lno
  end select

  ewrite(2,*)("Making bdy numbering");
  call make_boundary_numbering(mesh%bdy_list,mesh%bdy_nm_lno, &
       boundary_m_lno = mesh%bdy_nh_lno, &
       EElist = mesh%EEList, xnonod = N_Verts, &
       xndglno = mesh%EVlist_X, ele_n = mesh%nm, &
       boundary_n = mesh%nm_f, ele_m = mesh%nh, boundary_m = mesh%nh_f)
  
  ! Set up Mass matrices
  ewrite(2,*)("Calling posinm_dg");
  call POSINM_DG(mesh%Mass_m, N_elements, N_moms, mesh%nm%loc, mesh%evlist_m, &
       n_moms, mesh%nm%loc, mesh%evlist_m,mesh%bdy_list, &
       mesh%nm_f%loc,mesh%bdy_nm_lno)

  ! Set up Mass matrices
  ewrite(2,*)("Calling posinm");
  call POSINM(mesh%Mass_h, N_elements, N_dens, mesh%nh%loc, mesh%evlist_h,&
       N_dens, mesh%nh%loc, mesh%evlist_h)

  !flags to switch on different parts of the momentum equation
  ADVECTION_FLAG = .true.
  NONLINEAR_FLAG = .true.
  PRESSURE_FLAG = .true.

  !timestepping settings
  t = 0
  kappa = 0.
  mu = 0.
  dumpcount = 0
  filecount = 0
  call Get_Parameters()

  call get_vels_GN(Mesh, u, m, D, usebcs)

  timestep_loop: do

     !time increment
     t = t + dt
     ewrite(2,*)(t)
     ewrite(2,*)(tmax)

     ewrite(2,*) shape(mesh%nu%dn)

     call solve_density_equation(D,u1,u2,mesh)

     ewrite(2,*) shape(mesh%nu%dn)
     ewrite(2,*) 'go=',g0

     call solve_gn_momentum_equation(m,u,D,mesh,usebcs)

     !call get_vels(Mesh, D, u, m)

     !data output
     dumpcount = dumpcount + 1
     ewrite(2,*)(dumpcount)
     if(dumpcount==ndump) then
        filecount = filecount + 1
        pad = floor(log10(1.0*filecount))
        ewrite(2,*)(filecount)
        ewrite(2,*)(pad)

        call dump_data_vtk(u1,u2,m1,m2,D,mesh, &
             'gallopede_out_',filecount,4-pad,23)
        dumpcount = 0
     end if

     if(t>tmax) exit

  end do timestep_loop
  ewrite(2,*)'Time loop finished!'

end subroutine gn_main
