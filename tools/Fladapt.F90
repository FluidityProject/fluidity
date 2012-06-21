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

subroutine fladapt(input_basename_, input_basename_len, &
  & output_basename_, output_basename_len)  bind(c)
  !!< Peforms a mesh adapt based on the supplied input options file.
  !!< Outputs the resulting mesh.
 
  use iso_c_binding
  use adapt_state_module
  use diagnostic_fields_wrapper
  use diagnostic_fields_new, only : &
    & calculate_diagnostic_variables_new => calculate_diagnostic_variables
  use edge_length_module
  use field_options
  use fields
  use fldebug
  use limit_metric_module
  use metric_assemble
  use populate_state_module
  use reference_counting
  use spud
  use state_module
  use vtk_interfaces
  use mesh_files
  use populate_state_module

  implicit none
  
  interface
    subroutine check_options()
    end subroutine check_options

#ifdef HAVE_PYTHON
    subroutine python_init()
    end subroutine python_init
#endif
  end interface
  
  character(kind=c_char, len=1) :: input_basename_(*)
  integer(kind=c_size_t), value :: input_basename_len
  character(kind=c_char, len=1) :: output_basename_(*)
  integer(kind=c_size_t), value :: output_basename_len

  character(len=input_basename_len):: input_basename
  character(len=output_basename_len):: output_basename
  integer :: i
  type(mesh_type), pointer :: old_mesh
  type(state_type), dimension(:), pointer :: states
  type(vector_field) :: new_mesh_field
  type(vector_field), pointer :: new_mesh_field_ptr, old_mesh_field
  type(tensor_field) :: metric, t_edge_lengths
  character(len=FIELD_NAME_LEN) :: mesh_format

  ! now turn into proper fortran strings (is there an easier way to do this?)
  do i=1, input_basename_len
    input_basename(i:i)=input_basename_(i)
  end do
  do i=1, output_basename_len
    output_basename(i:i)=output_basename_(i)
  end do
  
  ewrite(1, *) "In fladapt"
 
#ifdef HAVE_PYTHON
  call python_init()
#endif

  ewrite(2, *) "Input base name: " // trim(input_basename)
  ewrite(2, *) "Output base name: " // trim(output_basename)

  ! Load the options tree
  call load_options(trim(input_basename) // ".flml")
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
  !call find_mesh_field_to_adapt(states(1), old_mesh_field)
  old_mesh_field => extract_vector_field(states(1), "Coordinate")
  old_mesh => old_mesh_field%mesh

  ! Assemble the error metric
  call allocate(metric, old_mesh, "ErrorMetric")
  call assemble_metric(states, metric)
  
  ewrite(0, *) "Expected nodes = ", expected_nodes(old_mesh_field, metric)
  
  call allocate(t_edge_lengths, metric%mesh, "TensorEdgeLengths")
  call get_edge_lengths(metric, t_edge_lengths)
  call vtk_write_fields(trim(output_basename) // "EdgeLengths", position = old_mesh_field, model = metric%mesh, &
    & tfields = (/t_edge_lengths, metric/))
  call deallocate(t_edge_lengths)
  
  ! Adapt the mesh
  call allocate_metric_limits(states(1))
  if(isparallel()) then
    call adapt_state(states, metric)
    call find_mesh_field_to_adapt(states(1), new_mesh_field_ptr)
    new_mesh_field = new_mesh_field_ptr
    call incref(new_mesh_field)
    new_mesh_field_ptr => null()
  else
    call adapt_mesh(old_mesh_field, metric, new_mesh_field)
  end if
  
  call get_option(trim(old_mesh_field%mesh%option_path)//"/from_file/format/name", mesh_format)
  ! Write the output mesh
  call write_mesh_files(output_basename, mesh_format, new_mesh_field)
  
  ! Deallocate
  do i = 1, size(states)
    call deallocate(states(i))
  end do
  call deallocate(metric)
  call deallocate(new_mesh_field)
  
  call print_references(0)
  
  ewrite(1, *) "Exiting fladapt"
  
contains

  subroutine find_mesh_field_to_adapt(state, mesh_field)
    !!< Find the external mesh field to be used by adaptivity

    type(state_type), intent(in) :: state
    type(vector_field), pointer :: mesh_field

    character(len = FIELD_NAME_LEN) :: mesh_field_name    
    type(mesh_type), pointer :: mesh
    
    call find_mesh_to_adapt(state, mesh)
    if(trim(mesh%name) == "CoordinateMesh") then
       mesh_field_name = "Coordinate"
    else
       mesh_field_name = trim(mesh%name) // "Coordinate"
    end if
    
    if(.not. has_vector_field(state, mesh_field_name)) then
      FLAbort("External mesh field " // trim(mesh_field_name) // " not found in the system state")
    end if
      
    mesh_field => extract_vector_field(state, mesh_field_name)
    
  end subroutine find_mesh_field_to_adapt
  
end subroutine fladapt
