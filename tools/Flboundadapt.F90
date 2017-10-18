!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

subroutine flboundadapt(input_flmlname_c, &
     input_geometryname_c, &
     output_meshname_c)  bind(c)
  !!< Peforms a mesh adapt based on the supplied input options file.
  !!< Outputs the resulting mesh.
 
  use iso_c_binding
  use futils
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN,&
       topology_mesh_name
  use fldebug
  use parallel_tools, only: isparallel
  use reference_counting
  use spud
  use fields
  use edge_length_module
  use state_module
  use field_options
  use vtk_interfaces
  use diagnostic_fields_wrapper, only: calculate_diagnostic_variables
  use limit_metric_module
  use mesh_files
  use Read_gmsh
  use populate_state_module
  use metric_assemble
  use adapt_state_module
  use diagnostic_fields_new, only : &
    & calculate_diagnostic_variables_new => calculate_diagnostic_variables

  implicit none
  
  interface
    subroutine check_options()
    end subroutine check_options

#ifdef HAVE_PYTHON
    subroutine python_init()
    end subroutine python_init
#endif
  end interface
  
  type(c_ptr), value :: input_flmlname_c, input_geometryname_c,&
       output_meshname_c

  character(len=FIELD_NAME_LEN):: input_flmlname
  character(len=FIELD_NAME_LEN):: input_geometryname
  character(len=FIELD_NAME_LEN):: output_meshname
  integer :: i, stat
  type(mesh_type), pointer :: old_mesh
  type(state_type), dimension(:), pointer :: states
  type(vector_field) :: new_mesh_field
  type(vector_field), pointer :: new_mesh_field_ptr, position
  type(tensor_field) :: metric, t_edge_lengths
  character(len=FIELD_NAME_LEN) :: mesh_format
  character(len=OPTION_PATH_LEN) :: mesh_file_name
  real, dimension(:,:,:), allocatable :: metric_data 

  ! now turn into proper fortran string
  call copy_c_string_to_fortran(input_flmlname_c, input_flmlname)
  call copy_c_string_to_fortran(input_geometryname_c, input_geometryname)
  call copy_c_string_to_fortran(output_meshname_c, output_meshname)
  
  ewrite(1, *) "In flboundadapt"
 
#ifdef HAVE_PYTHON
  call python_init()
#endif

  ewrite(2, *) "Input flml name: " // trim(input_flmlname)
  ewrite(2, *) "Input geometry name: " // trim(input_geometryname)
  ewrite(2, *) "Output mesh name: " // trim(output_meshname)

  ! Load the options tree
  call load_options(trim(input_flmlname) // ".flml")
  if(.not. have_option("/simulation_name")) then
    FLExit("Failed to find simulation name after loading options file")
  end if
  if(debug_level() >= 1) then
    ewrite(1, *) "Options tree:"
    call print_options()
  end if
#ifdef DDEBUG
  ewrite(1, *) "Performing options sanity check"
  call check_options()
  ewrite(1, *) "Options sanity check successful"
#endif

  ! Populate the system state
  call populate_state(states)

  ! Calculate diagnostic fields
  call calculate_diagnostic_variables(states)
  call calculate_diagnostic_variables_new(states)

  ! Find the external mesh field
  
  ! Note that Gmsh expects metric data on a mesh covering the geometry
  ! so this may not work for periodic data.

  position => extract_vector_field(states(1), "Coordinate")
  old_mesh => extract_mesh(states(1), topology_mesh_name)

  ! Assemble the error metric
  call allocate(metric, old_mesh, "ErrorMetric")
  call assemble_metric(states, metric)
  
  ewrite(0, *) "Expected nodes = ", expected_nodes(position, metric)
  
  call allocate(t_edge_lengths, metric%mesh, "TensorEdgeLengths")
  call get_edge_lengths(metric, t_edge_lengths)
  call vtk_write_fields(trim(output_meshname) // "EdgeLengths", position = position, model = position%mesh, &
    & tfields = (/t_edge_lengths, metric/))
  call deallocate(t_edge_lengths)
  
  ! Generate and the mesh

  new_mesh_field = read_gmsh_file(trim(input_geometryname),&
       format_string="geo", &
       quad_ngi=ele_ngi(old_mesh,1),&
       position = position, metric = metric)
  
  ! Write the output mesh
  call write_mesh_files(output_meshname, "gmsh", new_mesh_field)
  
  ! Deallocate
  do i = 1, size(states)
    call deallocate(states(i))
  end do
  call deallocate(metric)
  call deallocate(new_mesh_field)
  
  call print_references(0)
  
  ewrite(1, *) "Exiting flboundadapt"
   
end subroutine flboundadapt
