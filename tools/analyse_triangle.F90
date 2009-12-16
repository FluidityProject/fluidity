#include "fdebug.h"
program analyse_triangle
  !!< Read in a triangle mesh and output stats
  use fields
  use read_triangle
  use vtk_interfaces
  use fldebug
  implicit none
  type(vector_field), target :: positions
  character(len=256) :: filename
  integer :: status
  logical :: d3
  type(csr_matrix) :: NNList

  call get_command_argument(1, value=filename, status=status)
  
  select case(status)
  case(1:)
     call usage
     stop
  case(:-1)
     write(0,*) "Warning: truncating filename"
  end select

  positions=read_triangle_files(filename, quad_degree=3)
  d3 = .false.
  if(positions%dim==3) d3 = .true.

  call MakeLists(node_count(Positions),element_count(Positions), &
       positions%mesh%shape%loc,&
       positions%mesh%ndglno,&
       d3, NNList=NNList)

  ewrite(1,*) 'N vertices =', node_count(positions)
  ewrite(1,*) 'N elements =', element_count(positions)
  ewrite(1,*) 'N edges =', size(NNList%colm)/2
  
contains
  
  subroutine usage
    
    write (0,*) "usage: triangle2vtu <triangle_file_name>"
    write (0,*) "The triangle file name should be without the .node suffix."
    
  end subroutine usage

end program analyse_triangle
