module Global_parameters_gallopede

  !debugging level
  use global_parameters, only: current_debug_level

  !parameters set at beginning of program
  integer, parameter::nloc = 3
  integer, parameter::basis_degree = 1,quad_degree = 3, quad_degree_l = 2
  real, parameter::alpha1 = 8.0/64.0
  integer, parameter::X_=1,Y_=2

  !time limit
  real :: cpu_lim

  !tidal parameters
  real :: Tidal_T, Tidal_mag

  !global values
  !===================================================================

  integer :: N_Moms,N_Vels,N_Verts,N_Elements
  integer :: N_Dens,N_bc, N_layers, n_pres

  !timestepping stuff
  real :: tmax, t, dt, dtmax
  integer :: filecount=0
  real :: tdump, tdumpmax
  real :: theta , alph_up=1.0
  integer :: mom_maxnits, den_maxnits

  !physical parameters
  real:: kappa !diffusion in density equation
  real:: mu !viscosity in momentum equation
  real:: g0 !gravitational acceleration
  real :: f_0 =1e-4, beta =0e-11, b_0=1e-4

  !flags
  logical :: ADVECTION_FLAG
  logical :: NONLINEAR_FLAG, HYDROSTATIC_PRESSURE_FLAG
  logical :: NONHYDROSTATIC_PRESSURE_FLAG
  logical :: NONFLUX_FLAG, NOUPWIND_FLAG
  logical :: DENSITY_FLAG, MOMENTUM_FLAG
  logical :: constant_u, ROTATION_FLAG
  logical :: GN_FLAG, TIDAL_FORCING_FLAG
  logical :: RIGID_LID, barotropic_Split
  logical :: SHALLOW=.false., LINEAR=.false.



  !boundary condition locator
  integer, dimension(3,2), target :: offnods
  data offnods(1,:)/ 2, 3/
  data offnods(2,:)/ 1, 3/
  data offnods(3,:)/ 1, 2/
  
  !penalty term constant


  real :: eta

  !smoothing parameter
  real :: alpha 

  !Implicit_factor
  real :: IMPLICIT_FACTOR=0.0

end module Global_parameters_gallopede
