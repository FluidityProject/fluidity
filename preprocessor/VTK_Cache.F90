!    Copyright (C) 2007 Imperial College London and others.
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
!! module that implements a cache for reading fields from several vtu files
!! in atomic calls where the state objects of vtu files that have been
!! read before are cached.
module vtk_cache_module
use fields
use state_module
use vtk_interfaces
implicit none

! the vtu states that have previously been read:
! the name of the state refers to the filename
type(state_type), dimension(:), pointer:: vtk_states_read => null()

private

public :: vtk_cache_finalise, vtk_cache_read_positions_field, &
  & vtk_cache_read_scalar_field, vtk_cache_read_vector_field, &
  & vtk_cache_read_tensor_field

interface vtk_cache_read_positions_field
  module procedure vtk_cache_read_positions_field_unnamed, &
    & vtk_cache_read_positions_field_named
end interface vtk_cache_read_positions_field

contains

subroutine vtk_cache_finalise()

  integer:: i
  
  if (.not. associated(vtk_states_read)) return
  
  do i=1, size(vtk_states_read)
    call deallocate(vtk_states_read(i))
  end do
    
  deallocate( vtk_states_read )
  
end subroutine vtk_cache_finalise
  
function vtk_cache_read_file(filename) result(state)
!! If not yet present in cache, reads vtu file to state and adds it 
!! to the cache. In either case returns a pointer to this state.
type(state_type), pointer:: state
character(len=*), intent(in):: filename

  type(state_type), dimension(:), pointer:: vtk_states_old
  integer:: stat, n
    
  ! increase n/o states in cache by 1 and point state to last added
  if (associated(vtk_states_read)) then
    ! read already?
    state => extract_state(vtk_states_read, filename, stat=stat)
    if (stat==0) return
  
    ! if not, add a new state
    n=size(vtk_states_read)
    vtk_states_old => vtk_states_read
    allocate( vtk_states_read(1:n+1) )
    vtk_states_read(1:n)=vtk_states_old
    deallocate(vtk_states_old)
    
    state => vtk_states_read(n+1)

  else
    ! no read states yet, so create one
    allocate( vtk_states_read(1) )
    state => vtk_states_read(1)
  end if
  
  ewrite(1, *) "In vtk_cache_read_file, reading file: " // trim(filename)
  call vtk_read_state(filename, state)
  state%name=filename
  
end function vtk_cache_read_file

function vtk_cache_read_positions_field_unnamed(filename) result(positions)
  !!< Searches cache for state read from vtu file filename, if not present reads
  !!< vtu file and add to caches. Returns a pointer to the requested field.
  !!< Borrows a references from the cached state, i.e. don't deallocate!
  
  character(len = *), intent(in) :: filename
  
  type(vector_field), pointer :: positions
  
  positions => vtk_cache_read_positions_field(filename, positions_fieldname = "Coordinate")
  
end function vtk_cache_read_positions_field_unnamed

function vtk_cache_read_positions_field_named(filename, positions_fieldname) result(positions)
  !!< Searches cache for state read from vtu file filename, if not present reads
  !!< vtu file and add to caches. Returns a pointer to the requested field.
  !!< Borrows a references from the cached state, i.e. don't deallocate!
  
  character(len = *), intent(in) :: filename
  character(len = *), intent(in) :: positions_fieldname
  
  type(vector_field), pointer :: positions
  
  positions => vtk_cache_read_vector_field(filename, positions_fieldname)
  
end function vtk_cache_read_positions_field_named

function vtk_cache_read_scalar_field(filename, fieldname) result(field)
!! Searches cache for state read from vtu file filename, if not present
!! reads vtu file and adds to cache. Returns a pointer to the requested
!! field. Borrows a reference from the cached state, i.e. don't deallocate!
type(scalar_field), pointer:: field
character(len=*), intent(in):: filename
character(len=*), intent(in):: fieldname
  
  type(state_type), pointer:: state
  integer:: stat
  
  ewrite(2,*) "vtk_cache_read_scalar_field - filename, fieldname: ", &
    trim(filename), ", ", trim(fieldname)
    
  state => vtk_cache_read_file(filename)
  
  field => extract_scalar_field(state, fieldname, stat=stat)
  if (stat/=0) then
    ewrite(-1,*) "In vtk_cache_read_scalar_field"
    ewrite(-1,*) "filename: ", trim(filename)
    ewrite(-1,*) "fieldname: ", trim(fieldname)
    FLAbort("Requested field is not in the vtu")
  end if
  
end function vtk_cache_read_scalar_field

function vtk_cache_read_vector_field(filename, fieldname) result(field)
!! Searches cache for state read from vtu file filename, if not present
!! reads vtu file and adds to cache. Returns a pointer to the requested
!! field. Borrows a reference from the cached state, i.e. don't deallocate!
type(vector_field), pointer:: field
character(len=*), intent(in):: filename
character(len=*), intent(in):: fieldname
  
  type(state_type), pointer:: state
  integer:: stat
  
  ewrite(2,*) "vtk_cache_read_vector_field - filename, fieldname: ", &
    trim(filename), ", ", trim(fieldname)
  
  state => vtk_cache_read_file(filename)
  
  field => extract_vector_field(state, fieldname, stat=stat)
  if (stat/=0) then
    ewrite(-1,*) "In vtk_cache_read_vector_field"
    ewrite(-1,*) "filename: ", trim(filename)
    ewrite(-1,*) "fieldname: ", trim(fieldname)
    FLAbort("Requested field is not in the vtu")
  end if
  
end function vtk_cache_read_vector_field

function vtk_cache_read_tensor_field(filename, fieldname) result(field)
!! Searches cache for state read from vtu file filename, if not present
!! reads vtu file and adds to cache. Returns a pointer to the requested
!! field. Borrows a reference from the cached state, i.e. don't deallocate!
type(tensor_field), pointer:: field
character(len=*), intent(in):: filename
character(len=*), intent(in):: fieldname
  
  type(state_type), pointer:: state
  integer:: stat
  
  ewrite(2,*) "vtk_cache_read_tensor_field - filename, fieldname: ", &
    trim(filename), ", ", trim(fieldname)
    
  state => vtk_cache_read_file(filename)
  
  field => extract_tensor_field(state, fieldname, stat=stat)
  if (stat/=0) then
    ewrite(-1,*) "In vtk_cache_read_tensor_field"
    ewrite(-1,*) "filename: ", trim(filename)
    ewrite(-1,*) "fieldname: ", trim(fieldname)
    FLAbort("Requested field is not in the vtu")
  end if
  
end function vtk_cache_read_tensor_field

end module vtk_cache_module
