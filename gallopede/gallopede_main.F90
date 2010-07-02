#include "fdebug.h"

subroutine gallopede_main()

  !this stage: a solver for MLCM

  !coded by Colin Cotter, David Ham and James Percival (c) 2006-8
  !for last modification see CVS

  !links against libdfluidity
  use elements
  use sparse_tools
  use quadrature
  use global_numbering
  use shape_functions
  use global_parameters_gallopede
  use transform_elements
  use dgtools
  use text_io
  use vtk_io
 ! use density_equation
!  use mlcm_momentum_equation
!  use gallopede_barotropic_equation
!  use mlcm_operator
  use data_structures
!  use gallopede_solvers
  use mesh_tools
  use multilayer_tools
  use fields
  use sparsity_patterns
  use MLCM_scalar_formulation

  use mesh_files

  implicit none

  !local variables
  
  type(vector_field), pointer , dimension(:) :: u,m
  type(scalar_field), pointer,  dimension(:) :: D !height field
  type(vector_field), target :: u_bar
  type(scalar_field), target :: H,bottom !bottom boundary
  type(element_type), target :: X_shape, m_shape,u_shape, h_shape
  type(block_csr_matrix)  :: MLCM_mat
  real, dimension(:), allocatable :: rho, n1,n2,p
  real, dimension(:), allocatable :: nn1,nn2,np
  type(quadrature_type) :: g, g_f   ! Gauss quadrature
  type(dg_mesh) :: mesh
  type(bc_info), target :: bcs
  integer :: i,j,k, pad,ios,vertices(3)
  type(csr_sparsity) :: D_mass_sparsity, u_mass_sparsity,m_mass_sparsity
  real :: dt_out,dt_count
  real, dimension (:,:), allocatable :: test_mass
  CHARACTER(len=32) :: arg

  integer :: dim,loc,nodes,nelements,n_attributes,u_deg,jx

  real :: start_time, now_time

  namelist/command_flags/ ADVECTION_FLAG, NONLINEAR_FLAG,&
       HYDROSTATIC_PRESSURE_FLAG,NONFLUX_FLAG,&
       NONHYDROSTATIC_PRESSURE_FLAG, DENSITY_FLAG,&
       MOMENTUM_FLAG, NOUPWIND_FLAG, ROTATION_FLAG,&
       GN_FLAG

  call cpu_time(start_time)

  call Initialize_Petsc()

  !==============================================
  !process command-line options

  !Default values
          
  i = 0 ; j=2
  ADVECTION_FLAG = .true.
  NONLINEAR_FLAG = .true.
  HYDROSTATIC_PRESSURE_FLAG = .true.
  NONHYDROSTATIC_PRESSURE_FLAG = .true.
  NONFLUX_FLAG = .true.
  DENSITY_FLAG = .true.
  MOMENTUM_FLAG = .true.
  NOUPWIND_FLAG = .false.
  ROTATION_FLAG=.false.
  GN_FLAG=.false.
  CONSTANT_U=.false.
  TIDAL_FORCING_FLAG=.false.
  RIGID_LID=.false.
  BAROTROPIC_SPLIT=.false.
  
  
  do
     CALL get_command_argument(i, arg)
     print*, arg
     IF (LEN_TRIM(arg) == 0) EXIT
     select case(arg(1:1))
        
     case('-')
        j=1
        do
           if(j>len_trim(arg)) exit
           select case(arg(j:j))
           case('a')
              ADVECTION_FLAG=.false.
           case('b')
              BAROTROPIC_SPLIT=.true.
           case('f')
              NONFLUX_FLAG=.true.
           case('n')
              NONHYDROSTATIC_PRESSURE_FLAG=.false.
           case('d')
              DENSITY_FLAG=.false.
           case('g')
              GN_FLAG=.true.
           case('u')
              NOUPWIND_FLAG=.true.
           case('p')
              HYDROSTATIC_PRESSURE_FLAG=.false.
           case('m')
              MOMENTUM_FLAG=.false.
           case('l')
              NONLINEAR_FLAG=.false.
           case('r')
              RIGID_LID=.true.
           case('R')
              ROTATION_FLAG=.true.
           case('c')
              CONSTANT_U=.true.
           case('t')
              TIDAL_FORCING_FLAG=.true.
           case('s')
              SHALLOW=.true.
           case default

           end select
           j=j+1
        end do
     end select
     WRITE (*,*) TRIM(arg)
     i = i+1
  END DO
      

  ewrite(1,*) '****Parameters passed*********'

  if(.not. ADVECTION_FLAG) then
     ewrite(1,*) 'Advection terms off'
  end if
  if(.not. NONFLUX_FLAG) then
     ewrite(1,*) 'Flux version'
  end if
  if(.not. HYDROSTATIC_PRESSURE_FLAG) then 
     ewrite(1,*) 'Hydrostatic pressure Terms off'
  end if
  if(u_deg==2) then
     ewrite(1,*) 'Quadratic Velocities'
  end if
  if(.not. DENSITY_FLAG) then 
     ewrite(1,*) 'Density stepping off'
  end if
  if(.not. MOMENTUM_FLAG) then 
     ewrite(1,*) 'Momentum stepping off'
  end if
  if(.not. NOUPWIND_FLAG) then 
     ewrite(1,*) 'Upwinding stepping on'
  end if
  if(.not. NONLINEAR_FLAG) then 
     ewrite(1,*) 'Nonlinear Term off'
  end if
  if(.not. NONHYDROSTATIC_PRESSURE_FLAG) then 
     ewrite(1,*) 'Nonhydrostatic Term off'
  end if
  if(constant_u) then 
     ewrite(1,*) 'Constant U'
  end if
  if(ROTATION_FLAG) then 
     ewrite(1,*) 'Coriolis on'
  end if
  if(GN_FLAG) then 
     ewrite(1,*) 'GN terms only'
  end if
  if(TIDAL_FORCING_FLAG) then
     ewrite(1,*) 'Tidal Forcing on'
  end if
  if(RIGID_LID) then
     ewrite(1,*) 'Rigid Lid Imposed'
  end if
  if(BAROTROPIC_SPLIT) then
     ewrite(1,*) 'Timestep splitting on'
  end if
  if(SHALLOW) then
     ewrite(1,*) 'Shallow Water'
  end if

  open(unit=25,file='Input_flags.dat',iostat=ios)
  if(ios.ne.0) then
     write(0,*) ' Cannot write to Input_flags.dat'
     stop
  end if
  write (25,nml=command_flags)
  close(25)


  call identify_mesh_file('Positions',dim,loc,nodes,nelements,n_attributes)
  print*, dim,loc,nodes,nelements,n_attributes
  call identify_mesh_file('Variables',dim,loc,nodes,n_elements,n_attributes)
  u_deg = loc/3
  print*, 'u_deg=', u_deg
  N_vels = 3*nelements
  mesh%n_vels = N_vels
  N_dens = nodes
  print*, 'N_vels=', nodes
  N_moms = 3*nelements
  N_layers = (n_attributes-1)/3
  print*, 'N_layers=', n_layers

  !get quadrature
  !the 2d elements
  ewrite(2,*)("Making 2D quadrature");
  g=make_quadrature(loc=3,dimension=2,degree=8)
  ewrite(2,*)("Making 1D quadrature");
  g_f=make_quadrature(loc=2,dimension=1,degree=5)

  !get_basis functions
  !2D
  ewrite(2,*)("Making 2D P2 elements");
  allocate(mesh%nu,mesh%nm,mesh%nh,mesh%nx)
  mesh%nh=make_element_shape(loc=3, dimension=2, degree=2, quad=g,&
       quad_s=g_f)
  mesh%nu=make_element_shape(loc=3, dimension=2, degree=1, quad=g,&
       quad_s=g_f)

  !1D -- continuous basis functions
  ewrite(2,*)("Making 1D P1 elements");
  allocate(mesh%nu_f,mesh%nm_f,mesh%nh_f)
  mesh%nu_f=make_element_shape(loc=2, dimension=1, degree=1, quad=g_f)
  ewrite(2,*)("Making 1D P2 elements");
  mesh%nh_f=make_element_shape(loc=2, dimension=1, degree=2, quad=g_f)

  !load mesh data
  call Get_Parameters()
  allocate(m(n_layers),D(n_layers),u(n_layers))
  ewrite(2,*)("Loading mesh data");
  call Get_Mesh(mesh,bcs,m,u,d,u_bar,h,bottom) 


  ! get mass matrix for h

  ewrite(1,*) 'h sparsity'
  D_mass_sparsity=make_sparsity(D(1)%mesh,D(1)%mesh,name='Height_sparsity')
  ewrite(1,*) 'P sparsity'
  m_mass_sparsity=make_sparsity_transpose(m(1)%mesh,m(1)%mesh,name='Height_sparsity')
  u_mass_sparsity=make_sparsity(u(1)%mesh,u(1)%mesh,name='Height_sparsity')
  ewrite(1,*) 'mass matrices'
  call allocate(mesh%mass_h,D_mass_sparsity)
  call allocate(mesh%mass_u,u_mass_sparsity)
  call allocate(mesh%u_mat,u_mass_sparsity)

  call get_mass(mesh%mass_h,mesh%positions,s_field=D(1))
  call get_mass(mesh%mass_u,mesh%positions,v_field=u(1))

  t = 0
  kappa = 0.
  mu = 0.
  tdump = 0.
  filecount = 0
  ewrite(2,*)("Loading layer densities")
  allocate( rho(N_Layers) )
  call read_field(dr=rho, filename='rho.dat')


  ewrite(2,*)("Making matrix structure");

  call get_bc_marker(bcs,N_vels)

  call allocate(mesh%M_inv,(/1,2*n_layers+2/),(/node_count(u(1))/),&
       (/(node_count(u(1)),jx=1,2*n_layers+2)/))
  call allocate(mesh%mass_u_inv,node_count(u(1)),&
       node_count(u(1)))
  call allocate(mesh%C,(/2*n_layers,2*n_layers/),&
       (/(node_count(m(1)),jx=1,2*n_layers)/),&
       (/(node_count(u(1)),jx=1,2*n_layers)/))

  call get_ops(mesh%C,mesh%positions,&
       rho,D,m(1),u(1),bottom)

  do k=1,n_layers
     call get_weighted_dg_inverse_mass_matrix(&
          mesh%m_inv%blocks(1,1+2*(k-1)),&
          u(1)%mesh,mesh%positions,D(k),rho(k))
  end do

  call get_dg_inverse_mass_matrix(&
          mesh%mass_u_inv,&
          u(1)%mesh,mesh%positions)

  mesh%CMT=matmul_T(mesh%C%blocks(1,1),mesh%mass_u_inv)
  mesh%CMC=matmul_T(mesh%C%blocks(1,1),mesh%CMT)
 
  call allocate(mesh%D_mat,m_mass_sparsity)
  call allocate(mesh%MLCM_mat_block,m_mass_sparsity)

  call make_D_matrix(mesh,D(1),minval(D(1)%val))

  if (RIGID_LID) then
     call allocate(MLCM_mat,m_mass_sparsity,&
          (/2*n_layers+1,2*n_layers+1/))
  else
     call allocate(MLCM_mat,m_mass_sparsity,&
          (/2*n_layers,2*n_layers/))
  end if

! Output the initial conditions

!  call dump_data_vtk_layer_quad( &
!       u,m,D,bottom,mesh,'gallopede_out_00000',0,1,23)
!  call  energy_find(D,u,bottom,mesh,rho)
!  call get_model_u(mesh,u,D,(/50.0,350.0/),1.6,rho)
!  call get_barotropic_mode(mesh,D,u,H,u_bar,rho)

  ! temporary barotropic fiddle

!  allocate(np(2*node_count(D(1))))
!  np(1:node_count(D(1)))=d(1)%val
!  np(1+node_count(d(1)):2*node_count(D(1)))=d(2)%val

!  D(1)%val=D(1)%val+(300.0 +sqrt(300.0**2+4.0*rho(1)/rho(2)*350.0*50.0))&
!       /(2.0*rho(1)/rho(2)*350.0)*D(2)%val
!  D(2)%val= np(1:node_count(D(1)))+(300.0-sqrt(300.0**2+4.0*rho(1)/rho(2)*350.0*50.0))&
!       /(2.0*rho(1)/rho(2)*350.0)*np(1+node_count(d(1)):2*node_count(D(1)))

   call py_plot_write_state('gallopede_',&
               mesh%positions,m,u,u_bar,D,H,bottom,0,0.0) 

   ewrite(1,*) 'c_fast=', maxval(sqrt(g0*bottom%val))
   ewrite(1,*) 'CFL number =', get_cfl(mesh,maxval(sqrt(g0*bottom%val)),dt)   

  timestep_loop: do

     ewrite(1,*) '*******NEW TIMESTEP*********'


     call cpu_time(now_time)
     if(now_time-start_time>cpu_lim) t = tmax

     !Adapt timestep

     if (t/dt .ge. 0) then
        dt_count=1.0
     else
        dt_count=0.1
     end if

     !time increment
     t = t + dt_count*dt
     ewrite(1,*) 'T =', t
     ewrite(3,*)(tmax)
     if(t>tmax) exit

     ewrite(3,*)(N_DENS)
     ewrite(3,*)(N_ELEMENTS)

     ewrite(1,*) 'max D', maxval(D(n_layers)%val), 'max bottom', maxval(bottom%val)


!    Call main solution routine for a single timestep

!     call solve_MLCM_momentum_equation(m,u,D,bottom,rho,mesh,bcs)   

     call MLCM_timestep(MLCM_mat,mesh,bcs,u,u_bar,D,H,m,bottom,rho,dt*dt_count)


     !data output
     tdump = tdump + dt_count*dt
     ewrite(3,*)(tdump)
     ewrite(1,*) 'Writing to file'
     if(tdump.ge.tdumpmax) then
        filecount = filecount + 1
        pad = 4-floor(log10(1.0*filecount))
        ewrite(3,*)(filecount)
        ewrite(3,*)(pad)

!        call dump_data_vtk_layer_quad( &
!             u,m,D,bottom, &
!             mesh,'gallopede_out_',filecount,pad,23)
!        call get_barotropic_mode(mesh,D,u,H,u_bar,rho)
        call py_plot_write_state('gallopede_',&
               mesh%positions,m,u,u_bar,D,H,bottom,filecount,t)
        do
           tdump = tdump - tdumpmax
           if(tdump<tdumpmax) exit
        end do
     end if
  
  end do timestep_loop

  ewrite(1,*) 'INTEGRATION COMPLETE.'
  
end subroutine gallopede_main
