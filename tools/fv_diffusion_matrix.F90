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
!    C.Pain@Imperial.ac.uk
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

program fv_diffusion_matrix
  use spud
  use populate_state_module
  use global_parameters
  use state_module
  use unittest_tools
  use vtk_interfaces
  use advection_diffusion_fv
  use sparsity_patterns
  use fldebug
  use mpi_interfaces
  use boundary_conditions
  use reference_counting
  use dgtools

  implicit none

#ifdef HAVE_MPI
  integer :: ierr
#endif
  type(state_type), dimension(:), pointer :: states => null()

  character(len = 512) :: filename, simname
  type(csr_matrix) :: big_m, mass, inverse_mass
  type(csr_sparsity) :: sparsity, sparsity_1
  type(scalar_field) :: rhs, tracer
  type(vector_field), pointer :: x
  
#ifdef HAVE_MPI
  call MPI_Init(ierr)
  assert(ierr == MPI_SUCCESS)
#endif
  call python_init()

  call set_global_debug_level(0)

  call read_command_line(filename)

  call load_options(filename)
  call populate_state(states)
  call allocate_and_insert_auxilliary_fields(states)

  call get_option("/simulation_name", simname)

  tracer = extract_scalar_field(states(1), "Tracer")
  x => extract_vector_field(states(1), "Coordinate")

  call allocate(rhs, tracer%mesh, "RHS")

  sparsity = make_sparsity_transpose(tracer%mesh, tracer%mesh, "Sparsity")
  sparsity_1 = make_sparsity(tracer%mesh, tracer%mesh, "FirstOrderSparsity")

  call allocate(big_m, sparsity, name="Big_m")
  call zero(big_m)

  call allocate(mass, sparsity_1, name="Mass")
  call zero(mass)

  call construct_advection_diffusion_fv(big_m, rhs, "Tracer",&
       & states(1), mass)

  call get_dg_inverse_mass_matrix(inverse_mass, tracer%mesh, x)

  call apply_dirichlet_conditions(mass, rhs, tracer, dt)

  call mmwrite(trim(simname)//"_diffusion.mm", big_m)
  call mmwrite(trim(simname)//"_mass.mm", mass)
  call mmwrite(trim(simname)//"_inversemass.mm", inverse_mass)

  call deallocate(states(1))
  call deallocate(mass)
  call deallocate(big_m)
  call deallocate(rhs)
  call deallocate(sparsity)
  call deallocate(sparsity_1)
  call deallocate(inverse_mass)

  call print_references(0)

#ifdef HAVE_MPI
  call MPI_Finalize(ierr)
  assert(ierr == MPI_SUCCESS)
#endif

contains

  subroutine read_command_line(filename)
    character(len=*), intent(out) :: filename
    integer :: status

    call get_command_argument(1, value=filename, status=status)

    select case(status)
    case(1:)
       call usage
       stop
    case(:-1)
       write(0,*) "Warning: truncating filename"
    end select

  end subroutine read_command_line

  subroutine usage
    write(0,*) "usage:"
    write(0,*) "fv_diffusion_matrix flml-file"
    write(0,*) "Dumps the fv diffusion matrix."
  end subroutine usage

end program fv_diffusion_matrix
