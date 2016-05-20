subroutine test_alturbine_airfoils
    
  use alturbine_airfoils
  implicit none

  type(AirfoilType) :: airfoil
  real :: alpha75, alpha5, adotnorm, Re, umach, CL, CD, CN, CT, CLCirc, thtoc ,CM25  
  real :: CLa,CLCritP,CLCritN
  integer :: Iair, iStall
  write(*,*) "Hello from the Airfoil test"
 
  airfoil%afname="NACA_0015.dat"

  ! Set the values of interest
  alpha75 = 22*conrad
  alpha5  = 22*conrad
  adotnorm = 0
  Re = 2e4
  umach =0.0
  thtoc= 0.16
  
  iStall=1

  call airfoil_init_data(airfoil) 
  call compute_aeroCoeffs(airfoil,iStall,alpha75,alpha5,Re,adotnorm,umach,CL,CD,CN,CT,CLCirc,CM25)
  call CalcLBStallAOALim(airfoil,Re,CLa,CLCritP,CLCritN)

  write(*,*) CL, CD, CLa, CLCritP, CLCritN 

end subroutine test_alturbine_airfoils
