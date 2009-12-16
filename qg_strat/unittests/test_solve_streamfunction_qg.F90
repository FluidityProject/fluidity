subroutine test_solve_streamfunction_qg()

  use fields_calculations
  use fields_data_types
  use boundary_conditions
  use read_triangle
  use solvers
  use spud
  use state_module
  use unittest_tools
  use vtk_interfaces
  use pv_inversion

  character, parameter :: NEWLINE_CHAR=achar(10)
  character(len=*), parameter :: python_function1= &
       "def val(X,t):"//NEWLINE_CHAR// &
       "  from math import exp, cos, pi"//NEWLINE_CHAR// &
       "  cx = 100000.0"//NEWLINE_CHAR// &
       "  cy = 100000.0"//NEWLINE_CHAR// &
       "  H = 1000.0 "//NEWLINE_CHAR// &
       "  f = 1.0E-4"//NEWLINE_CHAR// &
       "  N = 0.1"//NEWLINE_CHAR// &
       "  cxsq=cx**2"//NEWLINE_CHAR// &
       "  cysq=cy**2"//NEWLINE_CHAR// &
       "  fsq = f**2"//NEWLINE_CHAR// &
       "  Nsq = N**2"//NEWLINE_CHAR// &
       "  Hsq = H**2"//NEWLINE_CHAR// &
       "  pisq = pi**2"//NEWLINE_CHAR// &
       "  return (-2.0/cxsq-2.0/cysq + 4.0*(X[0]/cxsq)**2 + 4.0*(X[1]/cysq)**2-(fsq*pisq)/(Hsq*Nsq))*exp(-(X[0]**2)/(cxsq)-(X[1]**2)/(cysq))*cos(pi*X[2]/H)"
  character(len=*), parameter :: solution1= &
       "def val(X,t):"//NEWLINE_CHAR// &
       "  from math import exp, cos, pi"//NEWLINE_CHAR// &
       "  cx = 100000.0"//NEWLINE_CHAR// &
       "  cy = 100000.0"//NEWLINE_CHAR// &
       "  H = 1000.0 "//NEWLINE_CHAR// &
       "  cxsq=cx**2"//NEWLINE_CHAR// &
       "  cysq=cy**2"//NEWLINE_CHAR// &
       "  return exp(-(X[0]**2)/(cxsq)-(X[1]**2)/(cysq))*cos(pi*X[2]/H)"
  character(len=*), parameter :: solution2= &
       "def val(X,t):"//NEWLINE_CHAR// &
       "   from math import sin, cos, pi"//NEWLINE_CHAR// &
       "   beta=-2.0e-11"//NEWLINE_CHAR// &
       "   H=0.1"//NEWLINE_CHAR// &
       "   v=(sin(pi*X[0]))**5*(sin(2.0*pi*X[1]))**5*cos(pi*X[2]/H)"//NEWLINE_CHAR// &
       "   return v+beta*X[1]"
  character(len=*), parameter :: python_function2= &
       "def val(X,t):"//NEWLINE_CHAR// &
       "   from math import sin, cos, pi"//NEWLINE_CHAR// &
       "   H=0.1"//NEWLINE_CHAR// &
       "   f=1.0E-4"//NEWLINE_CHAR// &
       "   N=1.0E-1"//NEWLINE_CHAR// &
       "   vxx=(20.0*pi**2*(sin(pi*X[0]))**3*(cos(pi*X[0]))**2-5.0*pi**2*(sin(pi*X[0]))**5)*sin(2.0*pi*X[1])**5*cos(pi*X[2]/H)"//NEWLINE_CHAR// &
       "   vyy=(80.0*pi**2*(sin(2.0*pi*X[1]))**3*(cos(2.0*pi*X[1]))**2-20.0*pi**2*(sin(2.0*pi*X[1]))**5)*sin(pi*X[0])**5*cos(pi*X[2]/H)"//NEWLINE_CHAR// &
       "   vzz=-(pi/H)**2*(sin(pi*X[0]))**5*(sin(2.0*pi*X[1]))**5*cos(pi*X[2]/H)"//NEWLINE_CHAR// &
       "   return vxx+vyy+((f/N)**2)*vzz"

  real :: l2error
  logical fail

  print*, "in test_solve_streamfunction_qg"

  call test_solve_streamfunction_qg_from_file(python_function2, &
       solution2, l2error)
  print *, l2error
  fail=l2error>1e-2
  call report_test("[test_solve_streamfunction_qg]", fail, .false., &
       "incorrect streamfunction calculated from PV")

contains

subroutine test_solve_streamfunction_qg_from_file(python_function, &
     solution, l2error)

  character(len=*), intent(in) :: python_function
  character(len=*), intent(in) :: solution
  real, intent(out) :: l2error

  type(state_type) :: state
  type(scalar_field) :: PV, streamfunction, solution_field, error_field
  type(scalar_field) :: streamfunction_surface_field
  type(vector_field) :: positions
  type(mesh_type), pointer :: surface_mesh
  integer :: stat

  ! Read positions and insert Coordinate field into state
  positions=read_triangle_files("data/box", quad_degree=4)
  call insert(state, positions, "Coordinate")

  ! Allocate and insert PV and streamfunction
  call allocate(PV, positions%mesh, "PotentialVorticity")
  call insert(state, PV, "PotentialVorticity")
  call allocate(streamfunction, positions%mesh, "Streamfunction")
  call set_solver_options(streamfunction, ksptype='gmres', pctype='sor', rtol=1.0e-7, max_its=3000)
  call add_boundary_condition(field=streamfunction, name="zero_streamfn_on_sides", type="dirichlet", boundary_ids=(/33/))
  call get_boundary_condition(field=streamfunction, name="zero_streamfn_on_sides", surface_mesh=surface_mesh)
  call allocate(streamfunction_surface_field, surface_mesh, name="value")
  streamfunction_surface_field%val(:)=0.0
  call insert_surface_field(streamfunction, 1, streamfunction_surface_field)
  call deallocate(streamfunction_surface_field)

  call insert(state, streamfunction, "Streamfunction")

  ! Initialise PV and zero streamfunction
  call set_from_python_function(PV, python_function, &
       positions, time=0.0)
  call zero(streamfunction)

  ! Set required options
  call set_option("/physical_parameters/coriolis/f_plane/f_0", 1.0E-4, stat)
  call set_option("/physical_parameters/buoyancy_frequency", 0.1, stat)

  call solve_streamfunction_qg(state, PV)

  ! Compare with solution
  call allocate(error_field, positions%mesh, &
       name="ErrorField")
  call set_from_python_function(error_field, solution, &
       positions, time=0.0)
  call allocate(solution_field, positions%mesh, &
       name="SolutionField")
  call set_from_python_function(solution_field, solution, &
       positions, time=0.0)
  call addto(error_field, streamfunction, scale=-1.0)

  l2error = norm2(error_field, positions)

  call vtk_write_fields("data/streamfn_out", 0, positions, positions%mesh, &
       sfields=(/PV, streamfunction, solution_field, error_field/))

end subroutine test_solve_streamfunction_qg_from_file

end subroutine test_solve_streamfunction_qg
