#include "fdebug.h"

subroutine test_streamfunction2velocity()

  use fields_base
  use fields_manipulation
  use fields_calculations
  use global_parameters, only: PYTHON_FUNC_LEN
  use read_triangle
  use state_module
  use unittest_tools
  use vtk_interfaces
  use PV_inversion

  real, dimension(2) :: l2error
  integer i, j
  logical fail

  print*, "in test_streamfunction2velocity"

  call test_streamfunction2velocity_from_file( &
       mesh_file="meshes/basin", &
       streamfn_file="x-periodic_streamfunction.py", &
       solution_file="x-periodic_velocity.py", l2error=l2error, j=1)
  do i=1,2
     fail=l2error(i)>0.4
     print*, "testing component ", i, " of error"
     print*, l2error(i)
     call report_test("[test_streamfunction2velocity]", fail, .false., &
          "incorrect velocity calculated from streamfunction")
  end do

  call test_streamfunction2velocity_from_file( &
       mesh_file="meshes/refined_basin", &
       streamfn_file="x-periodic_streamfunction.py", &
       solution_file="x-periodic_velocity.py", l2error=l2error, j=2)
  do i=1,2
     fail=l2error(i)>0.25
     print*, "testing component ", i, " of error"
     print*, l2error(i)
     call report_test("[test_streamfunction2velocity]", fail, .false., &
          "incorrect velocity calculated from streamfunction")
  end do

contains

subroutine test_streamfunction2velocity_from_file(mesh_file, & 
     streamfn_file, solution_file, l2error, j)

  character(len=*), intent(in) :: mesh_file
  character(len=*), intent(in) :: streamfn_file
  character(len=*), intent(in) :: solution_file
  real, dimension(2) , intent(out) :: l2error
  integer, intent(in) :: j

  type(state_type) :: state
  type(scalar_field) :: streamfunction
  type(scalar_field), dimension(2) :: e, s
  type(vector_field) :: positions, velocity, error_field, solution_field
  type(mesh_type) :: v_mesh
  integer :: i
  character(len=PYTHON_FUNC_LEN) :: streamfn, solution

  ! Read positions
  positions=read_triangle_files(trim(mesh_file), quad_degree=4)

  ! Insert coordinate field into state
  call insert(state, positions, "Coordinate")

  ! Allocate and insert streamfunction
  call allocate(streamfunction, positions%mesh, "Streamfunction")
  call insert(state, streamfunction, "Streamfunction")

  ! Make DG mesh for velocity
  v_mesh=make_mesh(positions%mesh, positions%mesh%shape, continuity=-1, &
       name="VelocityMesh")

  ! Allocate and insert velocity
  call allocate(velocity, mesh_dim(positions), v_mesh, "NonlinearVelocity")
  call insert(state, velocity, "NonlinearVelocity")

  ! Initialise streamfunction and velocity
  streamfn=python_function(streamfn_file)
  call set_from_python_function(streamfunction, streamfn, &
       positions, time=0.0)
  call zero(velocity)

  call streamfunction2velocity(state)

  ! Compare with solution
  call allocate(error_field,  mesh_dim(positions), v_mesh, name="ErrorField")
  call allocate(solution_field,  mesh_dim(positions), v_mesh, name="SolutionField")
  solution=python_function(solution_file)
  call set_from_python_function(error_field, solution, &
       positions, time=0.0)
  call set_from_python_function(solution_field, solution, &
       positions, time=0.0)
  call addto(error_field, velocity, scale=-1.0)

  do i=1,2
     e(i) = extract_scalar_field(error_field, i)
     s(i) = extract_scalar_field(solution_field, i)
     l2error(i) = norm2(e(i), positions)/norm2(s(i), positions)
  end do

  call vtk_write_fields("output/streamfunction2velocity_out", j, positions, positions%mesh, &
       sfields=(/streamfunction/), vfields=(/velocity, solution_field, error_field/))

end subroutine test_streamfunction2velocity_from_file

function python_function(filename)

  !!< This function returns a python function read in from a file (assumed 
  !!< to live in the "python_functions" directory.
  !!< Newline characters are added at the end of each line as required.
  character(len=*), intent(in) :: filename

  character(len=PYTHON_FUNC_LEN) :: python_function, buffer
  logical :: file_exists
  integer :: unit

  inquire(file="python_functions/"//trim(filename), exist=file_exists)  
  if (.not.file_exists) then
     ewrite(-1,*) "Couldn't find file: ", "python_functions/"//trim(filename) 
     FLAbort('Couldnt find file')
  end if
  unit=free_unit() 
  open(unit, file="python_functions/"//trim(filename), action="read",&
       & status="old")
  read(unit, '(a)', end=43) python_function
  ! Read all the lines of the file and put in newlines between them.
  do
     read(unit, '(a)', end=43) buffer
     python_function=trim(python_function)//achar(10)//trim(buffer)
  end do
43 python_function=trim(python_function)
  close(unit)

end function python_function

end subroutine test_streamfunction2velocity
