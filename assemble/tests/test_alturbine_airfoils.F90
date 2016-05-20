subroutine test_alturbine_airfoils
  !use spud
  !use populate_state_module
  !use global_parameters
  !use state_module
  !use unittest_tools
  !use vtk_interfaces
  !use advection_diffusion_dg
  !use sparsity_patterns
    
  use alturbine_airfoils
  implicit none

  type(AirfoilType) :: airfoils(3)
  real :: alpha75, alpha5, adotnorm, Re, umach, CL, CD, CN, CT, CLCirc, thtoc ,CM25  
  integer :: Iair,i
  write(*,*) "Hello from the Airfoil test"
 
  airfoils(1)%afname="NACA_0015.dat"
  airfoils(2)%afname="NACA_0018.dat"
  airfoils(3)%afname="NACA_0021.dat"

  ! Set the values of interest
  alpha75 = 22*conrad
  alpha5  = 22*conrad
  adotnorm = 0
  Re = 4e4
  umach =0.0
  thtoc= 0.16

  do Iair=1,3
  call airfoil_init(airfoils(Iair))
  end do
  

  write(*,*) airfoils(1:3)%tc
  do i=1,10
  call compute_aeroCoeffs_one_airfoil(airfoils(1),alpha75,alpha5,Re,adotnorm,umach,CL,CD,CN,CT,CLCirc,CM25)
  write(*,*) CL, CD 
  end do
end subroutine test_alturbine_airfoils
