#include "fdebug.h"

module Interpolation_ensemble_state_on_supermesh(ensemble_state_new, ensemble_state_old,supermesh)
  use fields
  use state_module
  use vtk_interfaces
  use pseudo_supermesh
  use interpolation_manager
  use spud
  implicit none

  type(state_type), dimension(:), intent(inout) :: ensemble_state_new,ensemble_state_old
  type(vector_field) :: initial_positions, out_positions
  character(len=255) :: filename
  character(len=255), dimension(:), allocatable :: files
  integer :: argc
  integer :: i, status
  type(state_type) :: initial_state
  integer :: ierr
  integer :: stat, mxnods
  integer :: nrens
  integer :: field_count

  subroutine interpolation_ensembles_on_supermesh

    nrens = size(ensemble_state_old)
    
    mxnods = 100000
    
    call mpi_init(ierr)
    call set_option('/mesh_adaptivity/hr_adaptivity/maximum_number_of_nodes', mxnods, stat=stat)
    
    allocate(files(nrens)
    do i=1, nrens
       !! Note that this won't work in parallel. Have to look for the pvtu in that case.
       write(files(i), '(a, i0, a)') trim(simulation_name)//'_', i, ".vtu"
    end do
    
    call vtk_read_state(trim(files(1)), initial_state)
    initial_positions = extract_vector_field(initial_state, "Coordinate")
    call add_faces(initial_positions%mesh)
    
    call compute_pseudo_supermesh(files, initial_positions, out_positions, mxnods=mxnods)
    
    !! call vtk_write_fields("pseudo_supermesh", 0, out_positions, out_positions%mesh)

    call insert(ensemble_state_new, out_positions, "Coordinate")
    call insert(ensemble_state_new, out_positions%mesh, "Mesh")

    !! Double check whether ensemble_state_new is setup correctly?? miss fields???
    
    do i=1, nrens
       !!sub linear_interpolation_scalars(old_fields, old_position, new_fields, new_position, map)
       !!sub quadratic_interpolation_qf(old_fields, old_position, new_fields, new_position)
       !!sub subroutine cubic_interpolation_cf_scalar(old_fields, old_position, new_fields, new_position)
       !!out:new_fields
       call linear_interpolation(ensemble_state_old(i), ensemble_state_new)
    enddo
    
  end subroutine interpolation_ensembles_on_supermesh
  
end module Interpolation_ensemble_state_on_supermesh
