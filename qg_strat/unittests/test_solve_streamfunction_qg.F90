#include "fdebug.h"

subroutine test_solve_streamfunction_qg()

  use boundary_conditions
  use fields_data_types
  use fields_calculations
  use global_parameters, only: PYTHON_FUNC_LEN
  use populate_state_module
  use pv_inversion
  use read_triangle
  use solvers
  use spud
  use state_module
  use vtk_interfaces
  use unittest_tools

  implicit none

  type(state_type), pointer, dimension(:) :: state

  real :: l2error
  logical fail

  print*, "In test_solve_streamfunction_qg"

  allocate(state(1))

  ! box
  print*, "Box with psi=0 on all sides:"
  call setup_state(state, mesh_filename="meshes/box", &
       pv_filename="pv_basin.py")
  call setup_dirichlet_boundary_conditions(state(1), ids=(/3, 4, 5, 6/), &
       bc_name="zero_streamfunction", constant=0.0)
  call setup_physical_parameters(state(1), f=1.0e-4, N=0.1, &
       beta_vector=(/0.0, -2.0e-11, 0.0/))
  call calculate_and_compare_streamfunction_qg(state(1), &
       "streamfunction_basin.py", l2error)
  print *, l2error
  fail=l2error>0.15
  call report_test("[test_solve_streamfunction_qg]", fail, .false., &
       "incorrect streamfunction calculated from PV")
  call vtk_write_state("output/streamfn_basin", state=state)

  call deallocate(state(1))
  allocate(state(1))

  ! periodic f-plane channel
  print*, "Periodic channel on f-plane:"
  call setup_state(state, mesh_filename="meshes/box", &
       pv_filename="periodic_f-plane_PV.py", &
       first_periodic_mapping="x-periodic_mapping.py", &
       first_physical_boundary_ids=(/5/), first_aliased_boundary_ids=(/6/))
  call setup_dirichlet_boundary_conditions(state(1), ids=(/3, 4/), &
       bc_name="zero_streamfunction", constant=0.0)
  call setup_physical_parameters(state(1), f=1.0e-4, N=0.1, &
       beta_vector=(/0.0, 0.0, 0.0/))
  call calculate_and_compare_streamfunction_qg(state(1), &
       "x-periodic_streamfunction.py", l2error)
  print *, l2error
  fail=l2error>1e-2
  call report_test("[test_solve_streamfunction_qg]", fail, .false., &
       "incorrect streamfunction calculated from PV")
  call vtk_write_state("output/streamfn_periodic_f-plane_channel", state=state)

  call deallocate(state(1))
  allocate(state(1))

  ! doubly periodic f-plane
  print*, "Doubly periodic problem on f-plane:"
  call setup_state(state, mesh_filename="meshes/basin", &
       pv_filename="doubly-periodic_f-plane_PV.py", &
       first_periodic_mapping="x-periodic_mapping.py", &
       first_physical_boundary_ids=(/5/), first_aliased_boundary_ids=(/6/), &
       second_periodic_mapping="y-periodic_mapping.py", &
       second_physical_boundary_ids=(/3/), second_aliased_boundary_ids=(/4/))
  call setup_physical_parameters(state(1), f=1.0e-4, N=0.1, &
       beta_vector=(/0.0, 0.0, 0.0/))
  call calculate_and_compare_streamfunction_qg(state(1), &
       "doubly-periodic_streamfunction.py", l2error)
  print *, l2error
  fail=l2error>1e-2
  call report_test("[test_solve_streamfunction_qg]", fail, .false., &
       "incorrect streamfunction calculated from PV")
  call vtk_write_state("output/streamfn_doubly-periodic_f-plane", state=state)

  call deallocate(state(1))

contains

subroutine calculate_and_compare_streamfunction_qg(state, &
     filename, l2error)

  type(state_type), intent(inout) :: state
  character(len=*), intent(in) :: filename
  real, intent(out) :: l2error

  type(scalar_field) :: solution_field, error_field
  type(scalar_field) :: beta_term
  type(scalar_field), pointer :: streamfunction, PV
  type(vector_field), pointer :: positions, beta
  character(len=PYTHON_FUNC_LEN) :: solution

  streamfunction=>extract_scalar_field(state, "Streamfunction")
  positions=>extract_vector_field(state, "Coordinate")

  call solve_streamfunction_qg(state)

  ! Compare with solution
  solution=trim(python_function(filename))
  call allocate(error_field, streamfunction%mesh, &
       name="ErrorField")
  call set_from_python_function(error_field, solution, &
       positions, time=0.0)
  call allocate(solution_field, streamfunction%mesh, &
       name="SolutionField")
  call set_from_python_function(solution_field, solution, &
       positions, time=0.0)
  call insert(state, solution_field, "Solution")
  call deallocate(solution_field)
  call addto(error_field, streamfunction, scale=-1.0)
  call insert(state, error_field, "Error")
  call deallocate(error_field)

  if(norm2(solution_field, positions) >= 1.0e-12) then
     l2error = norm2(error_field, positions)/norm2(solution_field, positions)
  else
     l2error = norm2(error_field, positions)
  end if

end subroutine calculate_and_compare_streamfunction_qg


subroutine setup_state(state, mesh_filename, pv_filename, first_periodic_mapping, first_physical_boundary_ids, first_aliased_boundary_ids, second_periodic_mapping, second_physical_boundary_ids, second_aliased_boundary_ids)

  !!< This subroutine gets the coordinate mesh and makes any required
  !! periodic meshes. It also allocates, initialises and inserts the PV
  !! and streamfunction fields.

  type(state_type), dimension(:), intent(inout) :: state
  character(len=*), intent(in) :: mesh_filename, pv_filename

  character(len=*), intent(in), optional :: first_periodic_mapping
  character(len=*), intent(in), optional :: second_periodic_mapping
  integer, dimension(:), intent(in), optional :: first_physical_boundary_ids
  integer, dimension(:), intent(in), optional :: first_aliased_boundary_ids
  integer, dimension(:), intent(in), optional :: second_physical_boundary_ids
  integer, dimension(:), intent(in), optional :: second_aliased_boundary_ids

  type(mesh_type) :: periodic_mesh, doubly_periodic_mesh
  type(scalar_field) :: PV, streamfunction, PV_perturbation
  type(vector_field) :: positions, periodic_positions
  character(len=PYTHON_FUNC_LEN) :: pv_function, periodic_mapping

  ! nullify state and clear options tree
  call nullify(state)
  call clear_options()

  ! Read and insert positions and mesh
  positions=read_triangle_files(mesh_filename, quad_degree=4)
  call insert(state, positions, "Coordinate")
  call insert(state, positions%mesh, "CoordinateMesh")

  ! Make other meshes if required
  if(present(first_periodic_mapping)) then
     periodic_mapping=python_function(first_periodic_mapping)
     periodic_mesh=make_mesh_periodic(positions%mesh, positions, first_physical_boundary_ids, first_aliased_boundary_ids, periodic_mapping)

     if(present(second_periodic_mapping)) then
        call allocate(periodic_positions, positions%dim, periodic_mesh, "periodic_positions")
        periodic_mapping=python_function(second_periodic_mapping)
        call remap_field(positions, periodic_positions)
        doubly_periodic_mesh=make_mesh_periodic(periodic_mesh, periodic_positions, second_physical_boundary_ids, second_aliased_boundary_ids, periodic_mapping)
        PV%mesh=doubly_periodic_mesh
        streamfunction%mesh=doubly_periodic_mesh
        call deallocate(periodic_positions)
     else
        PV%mesh=periodic_mesh
        streamfunction%mesh=periodic_mesh
     end if
  else
     PV%mesh=positions%mesh
     streamfunction%mesh=positions%mesh
  end if

  ! Allocate PV and streamfunction
  call allocate(PV, PV%mesh, "PotentialVorticity")
  call allocate(streamfunction, streamfunction%mesh, "Streamfunction")
  call set_solver_options(streamfunction, ksptype='cg', &
       pctype='eisenstat', rtol=1.0e-7, max_its=3000)

  ! Initialise PV and zero streamfunction
  pv_function=python_function(pv_filename)
  call set_from_python_function(PV, pv_function, &
       positions, time=0.0)
  call zero(streamfunction)

  call insert(state, PV, "PotentialVorticity")
  call insert(state, streamfunction, "Streamfunction")

  call deallocate(PV)
  call deallocate(streamfunction)

  call deallocate(positions)

end subroutine setup_state

subroutine setup_dirichlet_boundary_conditions(state, ids, bc_name, constant, filename)  

  type(state_type), intent(inout) :: state
  integer, dimension(:), intent(in) :: ids
  character(len=*), intent(in) :: bc_name
  ! Optional input
  real, intent(in), optional :: constant
  ! Name of file where function setting boundary condition is defined
  character(len=*), intent(in), optional :: filename

  type(mesh_type), pointer :: surface_mesh
  type(scalar_field), pointer :: streamfunction
  type(scalar_field) :: streamfunction_surface_field
  type(vector_field), pointer :: positions
  type(vector_field) :: bc_positions
  integer, pointer, dimension(:) :: surface_element_list
  character(len=PYTHON_FUNC_LEN) :: bc_function

  positions => extract_vector_field(state, "Coordinate")
  streamfunction => extract_scalar_field(state, "Streamfunction")

  ! Constant Dirichlet conditions
  if(present(constant)) then
     call add_boundary_condition(field=streamfunction, name=bc_name, &
          type="dirichlet", boundary_ids=ids)
     call get_boundary_condition(field=streamfunction, name=bc_name, &
          surface_mesh=surface_mesh)
     call allocate(streamfunction_surface_field, surface_mesh, name="value")
     streamfunction_surface_field%val(:)=constant
     call insert_surface_field(streamfunction, bc_name, &
          streamfunction_surface_field)
     call deallocate(streamfunction_surface_field)
  end if

  ! Dirichlet conditions specified by python function
  if(present(filename)) then
     bc_function=python_function(filename)
     call add_boundary_condition(field=streamfunction, name=bc_name, &
          type="dirichlet", boundary_ids=ids)
     call get_boundary_condition(field=streamfunction, name=bc_name, &
          surface_mesh=surface_mesh, surface_element_list=surface_element_list)
     call allocate(bc_positions, positions%dim, surface_mesh)
     call remap_field_to_surface(positions, bc_positions, surface_element_list)
     call allocate(streamfunction_surface_field, surface_mesh, name="value")
     call set_from_python_function(streamfunction_surface_field, bc_function, &
          bc_positions, time=0.0)
     call insert_surface_field(streamfunction, bc_name, &
          streamfunction_surface_field)
     call deallocate(streamfunction_surface_field)
     call deallocate(bc_positions)
  end if

end subroutine setup_dirichlet_boundary_conditions

subroutine setup_physical_parameters(state, f, N, beta_vector, beta_filename)

  !!< The subroutine sets up the Coriolis and buoyancy frequency options.
  type(state_type), intent(inout) :: state
  real, intent(in) :: f, N
  ! Optional input:
  real, dimension(:), intent(in), optional :: beta_vector
  character(len=*), intent(in), optional :: beta_filename

  type(vector_field), pointer :: positions
  type(vector_field) :: beta
  integer :: stat
  character(len=PYTHON_FUNC_LEN) :: beta_function

  ! Set buoyancy frequency
  call set_option("/physical_parameters/buoyancy_frequency", N, stat)

  ! Setup beta vector field
  positions => extract_vector_field(state, "Coordinate")
  call allocate(beta, mesh_dim(positions), positions%mesh, 'Beta')
  call zero(beta)

  ! Set Coriolis options
  if(present(beta_vector)) then
     call set_option("/geometry/dimension", 3, stat)
     call set_option("/physical_parameters/coriolis/beta_plane/f_0", f, stat)
     call set_option("/physical_parameters/coriolis/beta_plane/beta/vector_field/prescribed/value/constant", beta_vector, stat)
     call initialise_field(beta, "/physical_parameters/coriolis/beta_plane/beta/vector_field/prescribed/value/", positions)
  else if(present(beta_filename)) then
     beta_function=python_function(beta_filename)
     call set_option("/geometry/dimension", 3, stat)
     call set_option("/timestepping/current_time", 0.0, stat)
     call set_option("/physical_parameters/coriolis/beta_plane/f_0", f, stat)
     call set_option("/physical_parameters/coriolis/beta_plane/beta/vector_field/prescribed/value/python", beta_function, stat)
     call initialise_field(beta, "/physical_parameters/coriolis/beta_plane/beta/vector_field/prescribed/value/", positions)
  else
     call set_option("/physical_parameters/coriolis/f_plane/f", f, stat)
  end if

  ! Insert and deallocate beta field
  call insert(state, beta, 'Beta')
  call deallocate(beta)

end subroutine setup_physical_parameters

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

end subroutine test_solve_streamfunction_qg
