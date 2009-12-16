#include "fdebug.h"

program mesh_readnadapt
  ! A small program to read in a triangle mesh, and adapt to an (isotropic)
  ! edgelength function provided by the user
  use read_triangle
  use fields
  use FEtools
  use DGtools
  use elements
  use sparse_tools
  use vtk_interfaces
  use transform_elements
  use solvers
  use petsc_tools
  use sparsity_patterns
  use signal_vars
  use state_module
  !use diagnostic_variables
  use solvers
  use python_state
  use global_parameters, only : current_debug_level, PYTHON_FUNC_LEN
  use spud 
  use adapt_integration
  use write_triangle

  implicit none
  type(vector_field), target :: positions, new_positions
  integer :: degree, quad_degree
  type(quadrature_type), target :: quad
  type(element_type), target :: X_shape
  type(scalar_field) :: h, new_h
  type(tensor_field) :: metric
  ! Arguments for handling the command line
  character(len=256) :: filename
  character(len=100) :: fmt,buffer
  character(len=PYTHON_FUNC_LEN), save :: func

  logical :: file_exists
  integer :: dim, loc, nnodes, nelements, node_attributes, unit, io1, node,&
       & status
  real, dimension(3,3) :: EYE

  current_debug_level = 2

  ewrite(1,*) 'program  '

  call Initialize_Petsc()

  call python_init

  call get_command_argument(1, value=filename, status=status)
  select case(status)
  case(1:)
     print *, 'ERK'
     stop
  case(:-1)
     write(0,*) "Warning: truncating filename"
  end select
  filename=trim(filename)

  ewrite(2,*) 'Getting triangle file information'

  call identify_triangle_file(trim(filename), dim, loc, nnodes, nelements, &
       node_attributes)

  quad_degree = 6

  ewrite(2,*) 'Getting quadrature'

  quad=make_quadrature(loc=loc, dimension=dim, degree=quad_degree)
  
  ewrite(2,*) 'Getting shape functions'

  ! Shape functions for positions (linear)
  X_shape=make_element_shape(loc=loc, dimension=dim, &
       degree=1, quad=quad)
  
  ewrite(2,*) 'reading mesh'
  ewrite(2,*) 'loc = ',loc,'dim = ',dim

  positions=read_triangle_files(trim(filename), X_shape)
  call allocate(h,positions%mesh,'h')

  ewrite(2,*) 'Setting h'

  ewrite(2,*) 'setting h from python'
  inquire(file=trim(filename)//".py",exist=file_exists)  
  if (.not.file_exists) FLAbort('Couldnt find ' // trim(filename) // '.py file')
  unit=free_unit() 
  open(unit, file=trim(filename)//".py", action="read",&
       & status="old")
  read(unit, '(a)', end=43) func
  ! Read all the lines of the file and put in newlines between them.
  do
     read(unit, '(a)', end=43) buffer
     func=trim(func)//achar(10)//trim(buffer)
  end do
43 func=trim(func)//achar(10)
  close(unit)

  ewrite(2,*) func

  call set_from_python_function(h,trim(func), positions, 0.0)
  
  call allocate(metric,positions%mesh,name='metric')
  call zero(metric)
  
  EYE = 0.0
  EYE(1,1) = 1.0
  EYE(2,2) = 1.0
  EYE(3,3) = 1.0

  do node = 1, node_count(h)
     call set(metric, node, 1/(h%val(node)**2)*EYE)
  end do

  call vtk_write_fields(trim(filename)//'_b4', &
       index=0, position=positions, &
       sfields=(/h/), vfields=(/positions/), & !tfields= (/metric/), &
       & model=positions%mesh)

  call adapt_mesh(positions, metric, new_positions)
  write(0,*) "node_count(positions) == ", node_count(positions)
  write(0,*) "node_count(new_positions) == ", node_count(new_positions)

  call deallocate(h)
  call allocate(h,new_positions%mesh,'h')
  
  call set_from_python_function(h,trim(func), new_positions, 0.0)

  call write_triangle_files(trim(filename)//'_out', & 
       new_positions)
  call vtk_write_fields(trim(filename)//'_out', &
       index=0, position=new_positions, &
       sfields=(/h/), vfields=(/new_positions/),model=new_positions%mesh)

end program mesh_readnadapt
