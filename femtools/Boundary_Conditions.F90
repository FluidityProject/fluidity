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
module boundary_conditions
use fields_data_types
use fields
use sparse_tools
use sparse_tools_petsc
use state_module
use spud
use global_parameters, only: OPTION_PATH_LEN
implicit none
    
  interface add_boundary_condition
     module procedure add_scalar_boundary_condition, &
       add_vector_boundary_condition
  end interface add_boundary_condition

  interface remove_boundary_condition
     module procedure remove_scalar_boundary_condition, &
       remove_vector_boundary_condition
  end interface remove_boundary_condition

  interface add_boundary_condition_surface_elements
     module procedure add_scalar_boundary_condition_surface_elements, &
       add_vector_boundary_condition_surface_elements
  end interface add_boundary_condition_surface_elements

  interface get_boundary_condition
     module procedure get_scalar_boundary_condition_by_number, &
       get_vector_boundary_condition_by_number, &
       get_scalar_boundary_condition_by_name, &
       get_vector_boundary_condition_by_name
  end interface
    
  interface insert_surface_field
     module procedure insert_scalar_surface_field, &
       insert_vector_surface_field, insert_scalar_surface_field_by_name, &
       insert_vector_surface_field_by_name, &
       insert_vector_scalar_surface_field, &
       insert_vector_scalar_surface_field_by_name
  end interface
    
  interface extract_surface_field
     module procedure extract_scalar_surface_field_by_number, &
       extract_vector_surface_field_by_number, &
       extract_scalar_surface_field_by_name, &
       extract_vector_surface_field_by_name       
  end interface
    
  interface extract_scalar_surface_field
     module procedure extract_vector_scalar_surface_field, &
       extract_vector_scalar_surface_field_by_name
  end interface extract_scalar_surface_field
    
  interface has_surface_field
     module procedure has_scalar_surface_field_specific, &
       has_vector_surface_field_specific
  end interface
    
  interface get_boundary_condition_count
     module procedure get_scalar_boundary_condition_count, &
       get_vector_boundary_condition_count
  end interface
    
  interface get_entire_boundary_condition
     module procedure get_entire_scalar_boundary_condition, &
       get_entire_vector_boundary_condition
  end interface

  interface get_boundary_condition_nodes
     module procedure get_scalar_boundary_condition_nodes, &
       get_vector_boundary_condition_nodes
  end interface

  interface set_reference_node
     module procedure set_reference_node_scalar
  end interface set_reference_node

  interface has_boundary_condition
     module procedure has_boundary_condition_scalar, &
       has_boundary_condition_vector
  end interface has_boundary_condition

  interface has_boundary_condition_name
     module procedure has_boundary_condition_name_scalar, &
       has_boundary_condition_name_vector
  end interface has_boundary_condition_name
  
  interface set_dirichlet_consistent
    module procedure set_dirichlet_consistent, set_dirichlet_consistent_scalar, &
                     set_dirichlet_consistent_vector
  end interface

  interface apply_dirichlet_conditions
     module procedure apply_dirichlet_conditions_scalar, &
                      apply_dirichlet_conditions_scalar_lumped, &
                      apply_dirichlet_conditions_vector, &
                      apply_dirichlet_conditions_vector_petsc_csr, &
                      apply_dirichlet_conditions_vector_component, &
                      apply_dirichlet_conditions_vector_lumped, &
                      apply_dirichlet_conditions_vector_component_lumped
  end interface apply_dirichlet_conditions    

  private
  public add_boundary_condition, add_boundary_condition_surface_elements, &
    get_boundary_condition, get_boundary_condition_count, &
    insert_surface_field, extract_surface_field, has_surface_field, &
    extract_scalar_surface_field, get_entire_boundary_condition, &
    has_scalar_surface_field, get_boundary_condition_nodes, &
    get_dg_surface_mesh, has_boundary_condition, has_boundary_condition_name, &
    set_reference_node, &
    get_periodic_boundary_condition, remove_boundary_condition, &
    set_dirichlet_consistent, apply_dirichlet_conditions, &
    derive_collapsed_bcs

contains

  subroutine add_scalar_boundary_condition(field, name, type, boundary_ids, &
    option_path)
  !!< Add boundary condition to scalar field
  type(scalar_field), intent(inout):: field
  !! all things should have a name
  character(len=*), intent(in):: name
  !! type can be any of: ...
  character(len=*), intent(in):: type
  !! boundary ids indicating the part of the surface to apply this b.c. to
  integer, dimension(:), intent(in):: boundary_ids
  !! path to options for this b.c. in the options tree
  character(len=*), optional, intent(in) :: option_path
  
    logical, dimension(1:size(boundary_ids)):: boundary_id_used
    integer, dimension(:), allocatable:: surface_element_list
    integer i, ele_count
    
    assert(associated(field%mesh%faces))
    
    allocate( surface_element_list(1:surface_element_count(field)) )
    
    ! generate list of surface elements where this b.c. is applied
    ele_count=0
    boundary_id_used=.false.
    do i=1, surface_element_count(field)
      if (any(boundary_ids==surface_element_id(field, i))) then
        ele_count=ele_count+1
        surface_element_list(ele_count)=i
        where (boundary_ids==surface_element_id(field, i))
          boundary_id_used=.true.
        end where
      end if
    end do
      
    if (.not. IsParallel() .and. .not. all(boundary_id_used)) then
      ewrite(0,*) "WARNING: for boundary condition: ", trim(name)
      ewrite(0,*) "added to field: ", trim(field%name)
      ewrite(0,*) "The following boundary ids were specified, but they don't appear in the surface mesh:"
      ewrite(0,*) pack(boundary_ids, mask=.not. boundary_id_used)
    end if
    
    call add_scalar_boundary_condition_surface_elements(field, name, type, &
      surface_element_list(1:ele_count), option_path=option_path)

    deallocate(surface_element_list)
    
  end subroutine add_scalar_boundary_condition
  
  subroutine add_scalar_boundary_condition_surface_elements(field, name, type, surface_element_list, &
    option_path)
  !!< Add boundary condition to scalar field
  type(scalar_field), intent(inout):: field
  !! all things should have a name
  character(len=*), intent(in):: name
  !! type can be any of: ...
  character(len=*), intent(in):: type
  !! list of surface elements where this b.c. is to be applied
  integer, dimension(:), intent(in):: surface_element_list
  !! path to options for this b.c. in the options tree
  character(len=*), optional, intent(in) :: option_path
    
    type(scalar_boundary_condition), pointer:: tmp_boundary_condition(:)
    integer nobcs
    
    assert(associated(field%mesh%faces))
    
    if (.not. associated(field%bc%boundary_condition)) then
      allocate(field%bc%boundary_condition(1))
      nobcs=1
    else
      nobcs=size(field%bc%boundary_condition)+1
      ! save existing b.c.'s
      tmp_boundary_condition => field%bc%boundary_condition
      ! allocate new array with 1 new entry
      allocate(field%bc%boundary_condition(nobcs))
      ! copy back existing entries
      field%bc%boundary_condition(1:nobcs-1)=tmp_boundary_condition
      ! deallocate old b.c. array
      deallocate(tmp_boundary_condition)
    end if
    
    call allocate(field%bc%boundary_condition(nobcs), field%mesh, &
      surface_element_list=surface_element_list, &
      name=name, type=type)
        
    if (present(option_path)) then
      field%bc%boundary_condition(nobcs)%option_path=option_path
    end if
    
  end subroutine add_scalar_boundary_condition_surface_elements
    
  subroutine add_vector_boundary_condition(field, name, type, boundary_ids, &
      applies, option_path)
  !!< Add boundary condition to vector field
  type(vector_field), intent(inout):: field
  !! all things should have a name
  character(len=*), intent(in):: name
  !! type can be any of: ...
  character(len=*), intent(in):: type
  !! boundary ids indicating the part of the surface to apply this b.c. to
  integer, dimension(:), intent(in):: boundary_ids
  !! boundary condition only applies to component with applies==.true.
  logical, dimension(:), intent(in), optional:: applies
  !! path to options for this b.c. in the options tree
  character(len=*), optional, intent(in) :: option_path
    
    logical, dimension(1:size(boundary_ids)):: boundary_id_used
    integer, dimension(:), allocatable:: surface_element_list
    integer i, ele_count
    
    assert(associated(field%mesh%faces))
    
    allocate( surface_element_list(1:surface_element_count(field)) )
    
    ! generate list of surface elements where this b.c. is applied
    ele_count=0
    boundary_id_used=.false.
    do i=1, surface_element_count(field)
      if (any(boundary_ids==surface_element_id(field, i))) then
        ele_count=ele_count+1
        surface_element_list(ele_count)=i
        where (boundary_ids==surface_element_id(field, i))
          boundary_id_used=.true.
        end where
      end if
    end do
        
    if (.not. IsParallel() .and. .not. all(boundary_id_used)) then
      ewrite(0,*) "WARNING: for boundary condition: ", trim(name)
      ewrite(0,*) "added to field: ", trim(field%name)
      ewrite(0,*) "The following boundary ids were specified, but they don't appear in the surface mesh:"
      ewrite(0,*) pack(boundary_ids, mask=.not. boundary_id_used)
    end if
    
    call add_vector_boundary_condition_surface_elements(field, name, type, &
      surface_element_list(1:ele_count), applies=applies, option_path=option_path)
        
    deallocate(surface_element_list)
    
  end subroutine add_vector_boundary_condition
    
  subroutine add_vector_boundary_condition_surface_elements(field, name, type, surface_element_list, &
      applies, option_path)
      
  !!< Add boundary condition to vector field
  type(vector_field), intent(inout):: field
  !! all things should have a name
  character(len=*), intent(in):: name
  !! type can be any of: ...
  character(len=*), intent(in):: type
  !! list of surface elements where this b.c. is to be applied
  integer, dimension(:), intent(in):: surface_element_list
  !! boundary condition only applies to component with applies==.true.
  logical, dimension(:), intent(in), optional:: applies
  !! path to options for this b.c. in the options tree
  character(len=*), optional, intent(in) :: option_path
  
    type(vector_boundary_condition), pointer:: tmp_boundary_condition(:)
    integer nobcs
    
    assert(associated(field%mesh%faces))
    
    if (.not. associated(field%bc%boundary_condition)) then
      allocate(field%bc%boundary_condition(1))
      nobcs=1
    else
      nobcs=size(field%bc%boundary_condition)+1
      ! save existing b.c.'s
      tmp_boundary_condition => field%bc%boundary_condition
      ! allocate new array with 1 new entry
      allocate(field%bc%boundary_condition(nobcs))
      ! copy back existing entries
      field%bc%boundary_condition(1:nobcs-1)=tmp_boundary_condition
      ! deallocate old b.c. array
      deallocate(tmp_boundary_condition)
    end if
    
    call allocate(field%bc%boundary_condition(nobcs), field%mesh, &
      surface_element_list=surface_element_list, &
      name=name, type=type, applies=applies)
    
    if (present(option_path)) then
      field%bc%boundary_condition(nobcs)%option_path=option_path
    end if

  end subroutine add_vector_boundary_condition_surface_elements
    
  subroutine remove_scalar_boundary_condition(field, name)
  !!< Removed boundary condition from scalar field
  type(scalar_field), intent(inout):: field
  character(len=*), intent(in):: name
  
    type(scalar_boundary_condition), pointer:: tmp_boundary_condition(:)
    integer:: i, nobcs
    
    if (associated(field%bc%boundary_condition)) then
      nobcs=size(field%bc%boundary_condition)
      do i=1, nobcs
        if (field%bc%boundary_condition(i)%name==name) then
          call deallocate(field%bc%boundary_condition(i))
          ! save existing b.c.'s
          tmp_boundary_condition => field%bc%boundary_condition
          ! allocate new array with 1 less entry
          allocate(field%bc%boundary_condition(nobcs-1))
          ! copy back existing ones, except the i-th
          field%bc%boundary_condition(1:i-1)=tmp_boundary_condition(1:i-1)
          field%bc%boundary_condition(i:)=tmp_boundary_condition(i+1:)
          ! deallocate the old bcs array
          deallocate( tmp_boundary_condition )
          return
        end if
      end do
    end if
    ewrite(-1,*) 'In remove_scalar_boundary_condition'
    ewrite(-1,*) 'Unknown boundary condition: ', name
    FLAbort("Sorry!")
  
  end subroutine remove_scalar_boundary_condition
  
  subroutine remove_vector_boundary_condition(field, name)
  !!< Removed boundary condition from vector field
  type(vector_field), intent(inout):: field
  character(len=*), intent(in):: name
  
    type(vector_boundary_condition), pointer:: tmp_boundary_condition(:)
    integer:: i, nobcs
    
    if (associated(field%bc%boundary_condition)) then
      nobcs=size(field%bc%boundary_condition)
      do i=1, nobcs
        if (field%bc%boundary_condition(i)%name==name) then
          call deallocate(field%bc%boundary_condition(i))
          ! save existing b.c.'s
          tmp_boundary_condition => field%bc%boundary_condition
          ! allocate new array with 1 less entry
          allocate(field%bc%boundary_condition(nobcs-1))
          ! copy back existing ones, except the i-th
          field%bc%boundary_condition(1:i-1)=tmp_boundary_condition(1:i-1)
          field%bc%boundary_condition(i:)=tmp_boundary_condition(i+1:)
          ! deallocate the old bcs array
          deallocate( tmp_boundary_condition )
          return
        end if
      end do
    end if
    ewrite(-1,*) 'In remove_vector_boundary_condition'
    ewrite(-1,*) 'Unknown boundary condition: ', name
    FLAbort("Sorry!")
  
  end subroutine remove_vector_boundary_condition
  
  subroutine insert_scalar_surface_field(field, n, surface_field)
  !!< Adds a surface_field to a boundary condition: a field over the 
  !!< part of the surface mesh that this b.c. applies to. This can be used
  !!< to store b.c. values
  !! field to add surface_field to
  type(scalar_field), intent(in):: field
  !! add to n-th b.c.
  integer, intent(in):: n
  !! field to insert, callers of this routine should deallocate their copy
  !! of this field afterwards.
  type(scalar_field), intent(in):: surface_field
    
    type(scalar_field), dimension(:), pointer:: tmp_surface_fields
    type(scalar_boundary_condition), pointer:: bc
    integer i

    assert(n>=1 .and. n<=size(field%bc%boundary_condition))
    bc => field%bc%boundary_condition(n)
    
    if (.not. associated(bc%surface_fields)) then
      allocate(bc%surface_fields(1))
      i=1
    else
      ! save existing surface fields
      tmp_surface_fields => bc%surface_fields
      ! allocate one extra
      i=size(bc%surface_fields)+1
      allocate(bc%surface_fields(i))
      ! copy back existing surface fields
      bc%surface_fields(1:i-1)=tmp_surface_fields
    end if
    
    bc%surface_fields(i)=surface_field
    
    ! To remain consistent with insert_field for state we incref here
    ! so users of this routine should deallocate their copy of the
    ! surface_field. 
    call incref(surface_field)
    
  end subroutine insert_scalar_surface_field

  subroutine insert_vector_surface_field(field, n, surface_field)
  !!< Adds a surface_field to a boundary condition: a field over the 
  !!< part of the surface mesh that this b.c. applies to. This can be used
  !!< to store b.c. values
  !! field to add surface_field to
  type(vector_field), intent(in):: field
  !! add to n-th b.c.
  integer, intent(in):: n
  !! field to insert, callers of this routine should deallocate their copy
  !! of this field afterwards.
  type(vector_field), intent(in):: surface_field
    
    type(vector_field), dimension(:), pointer:: tmp_surface_fields
    type(vector_boundary_condition), pointer:: bc
    integer i

    assert(n>=1 .and. n<=size(field%bc%boundary_condition))
    bc => field%bc%boundary_condition(n)
    
    if (.not. associated(bc%surface_fields)) then
      allocate(bc%surface_fields(1))
      i=1
    else
      ! save existing surface fields
      tmp_surface_fields => bc%surface_fields
      ! allocate one extra
      i=size(bc%surface_fields)+1
      allocate(bc%surface_fields(i))
      ! copy back existing surface fields
      bc%surface_fields(1:i-1)=tmp_surface_fields
    end if
    
    bc%surface_fields(i)=surface_field
    
    ! To remain consistent with insert_field for state we incref here
    ! so users of this routine should deallocate their copy of the
    ! surface_field. 
    call incref(surface_field)
    
  end subroutine insert_vector_surface_field
    
  subroutine insert_vector_scalar_surface_field(field, n, surface_field)
  !!< Adds a surface_field to a boundary condition: a field over the 
  !!< part of the surface mesh that this b.c. applies to. This can be used
  !!< to store b.c. values
  !! field to add surface_field to
  type(vector_field), intent(in):: field
  !! add to n-th b.c.
  integer, intent(in):: n
  !! field to insert, callers of this routine should deallocate their copy
  !! of this field afterwards.
  type(scalar_field), intent(in):: surface_field
    
    type(scalar_field), dimension(:), pointer:: tmp_surface_fields
    type(vector_boundary_condition), pointer:: bc
    integer i

    assert(n>=1 .and. n<=size(field%bc%boundary_condition))
    bc => field%bc%boundary_condition(n)
    
    if (.not. associated(bc%scalar_surface_fields)) then
      allocate(bc%scalar_surface_fields(1))
      i=1
    else
      ! save existing surface fields
      tmp_surface_fields => bc%scalar_surface_fields
      ! allocate one extra
      i=size(bc%scalar_surface_fields)+1
      allocate(bc%scalar_surface_fields(i))
      ! copy back existing surface fields
      bc%scalar_surface_fields(1:i-1)=tmp_surface_fields
    end if
    
    bc%scalar_surface_fields(i)=surface_field
    
    ! To remain consistent with insert_field for state we incref here
    ! so users of this routine should deallocate their copy of the
    ! surface_field. 
    call incref(surface_field)
    
  end subroutine insert_vector_scalar_surface_field
    
  subroutine insert_scalar_surface_field_by_name(field, name, surface_field)
  !!< Adds a surface_field to a boundary condition: a field over the 
  !!< part of the surface mesh that this b.c. applies to. This can be used
  !!< to store b.c. values
  type(scalar_field), intent(in):: field
  !! add to b.c. with name:
  character(len=*), intent(in):: name
  !! field to insert, callers of this routine should deallocate their copy
  !! of this field afterwards.
  type(scalar_field), intent(in):: surface_field
    
    integer i
    do i=1, size(field%bc%boundary_condition)
      if (field%bc%boundary_condition(i)%name==name) then
        call insert_scalar_surface_field(field, i, surface_field)
        return
      end if
    end do
    ewrite(-1,*) 'Unknown boundary condition: ', name
    FLAbort("Hasta la vista")

  end subroutine insert_scalar_surface_field_by_name
    
  subroutine insert_vector_surface_field_by_name(field, name, surface_field)
  !!< Adds a surface_field to a boundary condition: a field over the 
  !!< part of the surface mesh that this b.c. applies to. This can be used
  !!< to store b.c. values
  type(vector_field), intent(in):: field
  !! add to b.c. with name:
  character(len=*), intent(in):: name
  !! field to insert, callers of this routine should deallocate their copy
  !! of this field afterwards.
  type(vector_field), intent(in):: surface_field
    
    integer i
    do i=1, size(field%bc%boundary_condition)
      if (field%bc%boundary_condition(i)%name==name) then
        call insert_vector_surface_field(field, i, surface_field)
        return
      end if
    end do
    ewrite(-1,*) 'Unknown boundary condition: ', name
    FLAbort("Hasta la vista")

  end subroutine insert_vector_surface_field_by_name
  
  subroutine insert_vector_scalar_surface_field_by_name(field, name, surface_field)
  !!< Adds a surface_field to a boundary condition: a field over the 
  !!< part of the surface mesh that this b.c. applies to. This can be used
  !!< to store b.c. values
  type(vector_field), intent(in):: field
  !! add to b.c. with name:
  character(len=*), intent(in):: name
  !! field to insert, callers of this routine should deallocate their copy
  !! of this field afterwards.
  type(scalar_field), intent(in):: surface_field
    
    integer i
    do i=1, size(field%bc%boundary_condition)
      if (field%bc%boundary_condition(i)%name==name) then
        call insert_vector_scalar_surface_field(field, i, surface_field)
        return
      end if
    end do
    ewrite(-1,*) 'Unknown boundary condition: ', name
    FLAbort("Hasta la vista")

  end subroutine insert_vector_scalar_surface_field_by_name
    
  function extract_scalar_surface_field_by_number(field, n, name) result (surface_field)
  !!< Extracts one of the surface_fields by name of the n-th b.c. of field
  type(scalar_field), pointer:: surface_field
  type(scalar_field), intent(in):: field
  integer, intent(in):: n
  character(len=*), intent(in):: name
  
    type(scalar_boundary_condition), pointer:: bc
    integer i

    assert(n>=1 .and. n<=size(field%bc%boundary_condition))
    bc => field%bc%boundary_condition(n)
    
    if (associated(bc%surface_fields)) then
      do i=1, size(bc%surface_fields)
        if (bc%surface_fields(i)%name==name) then
          surface_field => bc%surface_fields(i)
          return
        end if
      end do
    end if
    
    ewrite(-1, '(a," is not a surface_field of boundary condition n=",i0," of field ",a)') trim(name), n, trim(field%name)
    FLAbort("Sorry!")
    
  end function extract_scalar_surface_field_by_number

  function extract_vector_surface_field_by_number(field, n, name) result (surface_field)
  !!< Extracts one of the surface_fields by name of the n-th b.c. of field
  type(vector_field), pointer:: surface_field
  type(vector_field), intent(in):: field
  integer, intent(in):: n
  character(len=*), intent(in):: name
  
    type(vector_boundary_condition), pointer:: bc
    integer i

    assert(n>=1 .and. n<=size(field%bc%boundary_condition))
    bc => field%bc%boundary_condition(n)
    
    if (associated(bc%surface_fields)) then
      do i=1, size(bc%surface_fields)
        if (bc%surface_fields(i)%name==name) then
          surface_field => bc%surface_fields(i)
          return
        end if
      end do
    end if
    
    ewrite(-1, '(a," is not a surface_field of boundary condition n=",i0," of field ",a)') trim(name), n, trim(field%name)
    FLAbort("Sorry!")
    
  end function extract_vector_surface_field_by_number

  function extract_vector_scalar_surface_field(field, n, name) result (surface_field)
  !!< Extracts one of the surface_fields by name of the n-th b.c. of field
  type(scalar_field), pointer:: surface_field
  type(vector_field), intent(in):: field
  integer, intent(in):: n
  character(len=*), intent(in):: name
  
  type(vector_boundary_condition), pointer:: bc
    integer i

    assert(n>=1 .and. n<=size(field%bc%boundary_condition))
    bc => field%bc%boundary_condition(n)
    
    if (associated(bc%scalar_surface_fields)) then
      do i=1, size(bc%scalar_surface_fields)
        if (bc%scalar_surface_fields(i)%name==name) then
          surface_field => bc%scalar_surface_fields(i)
          return
        end if
      end do
    end if
    
    ewrite(-1, '(a," is not a surface_field of boundary condition n=",i0," of field ",a)') trim(name), n, trim(field%name)
    FLAbort("Sorry!")
    
  end function extract_vector_scalar_surface_field
  
  function extract_scalar_surface_field_by_name(field, bc_name, name) result (surface_field)
  !!< Extracts one of the surface_fields of the b.c. 'bc_name' by name of field
  type(scalar_field), pointer:: surface_field
  type(scalar_field), intent(in):: field
  character(len=*), intent(in):: bc_name, name
  
    integer i
    do i=1, size(field%bc%boundary_condition)
      if (field%bc%boundary_condition(i)%name==bc_name) then
        surface_field => extract_scalar_surface_field_by_number(field, i, name)
        return
      end if
    end do
    ewrite(-1,*) 'Unknown boundary condition: ', name
    FLAbort("Sorry!")
    
  end function extract_scalar_surface_field_by_name
  
  function extract_vector_surface_field_by_name(field, bc_name, name) result (surface_field)
  !!< Extracts one of the surface_fields of the b.c. 'bc_name' by name of field
  type(vector_field), pointer:: surface_field
  type(vector_field), intent(in):: field
  character(len=*), intent(in):: bc_name, name
  
    integer i
    do i=1, size(field%bc%boundary_condition)
      if (field%bc%boundary_condition(i)%name==bc_name) then
        surface_field => extract_vector_surface_field_by_number(field, i, name)
        return
      end if
    end do
    ewrite(-1,*) 'Unknown boundary condition: ', name
    FLAbort("Sorry!")
    
  end function extract_vector_surface_field_by_name

  function extract_vector_scalar_surface_field_by_name(field, bc_name, name) result (surface_field)
  !!< Extracts one of the surface_fields of the b.c. 'bc_name' by name of field
  type(scalar_field), pointer:: surface_field
  type(vector_field), intent(in):: field
  character(len=*), intent(in):: bc_name, name
  
    integer i
    do i=1, size(field%bc%boundary_condition)
      if (field%bc%boundary_condition(i)%name==bc_name) then
        surface_field => extract_vector_scalar_surface_field(field, i, name)
        return
      end if
    end do
    ewrite(-1,*) 'Unknown boundary condition: ', name
    FLAbort("Sorry!")
    
  end function extract_vector_scalar_surface_field_by_name
  
  function has_scalar_surface_field_specific(field, n, name)
  !!< Tells whether a surface_field with the given is present
  logical :: has_scalar_surface_field_specific
  type(scalar_field), intent(in):: field
  integer, intent(in):: n
  character(len=*), intent(in):: name
  
    type(scalar_boundary_condition), pointer:: bc
    integer i

    assert(n>=1 .and. n<=size(field%bc%boundary_condition))
    bc => field%bc%boundary_condition(n)
    
    if (associated(bc%surface_fields)) then
      do i=1, size(bc%surface_fields)
        if (bc%surface_fields(i)%name==name) then
          has_scalar_surface_field_specific=.true.
          return
        end if
      end do
    end if
    
    has_scalar_surface_field_specific=.false.
    
  end function has_scalar_surface_field_specific
  
  function has_vector_surface_field_specific(field, n, name)
  !!< Tells whether a surface_field with the given is present
  logical :: has_vector_surface_field_specific
  type(vector_field), intent(in):: field
  integer, intent(in):: n
  character(len=*), intent(in):: name
  
    type(vector_boundary_condition), pointer:: bc
    integer i

    assert(n>=1 .and. n<=size(field%bc%boundary_condition))
    bc => field%bc%boundary_condition(n)
    
    if (associated(bc%surface_fields)) then
      do i=1, size(bc%surface_fields)
        if (bc%surface_fields(i)%name==name) then
          has_vector_surface_field_specific=.true.
          return
        end if
      end do
    end if
    
    has_vector_surface_field_specific=.false.
    
  end function has_vector_surface_field_specific

  function has_scalar_surface_field(field, n, name)
  !!< Tells whether a surface_field with the given is present
  logical :: has_scalar_surface_field
  type(vector_field), intent(in):: field
  integer, intent(in):: n
  character(len=*), intent(in):: name
  
    type(vector_boundary_condition), pointer:: bc
    integer i

    assert(n>=1 .and. n<=size(field%bc%boundary_condition))
    bc => field%bc%boundary_condition(n)
    
    if (associated(bc%scalar_surface_fields)) then
      do i=1, size(bc%scalar_surface_fields)
        if (bc%scalar_surface_fields(i)%name==name) then
          has_scalar_surface_field=.true.
          return
        end if
      end do
    end if
    
    has_scalar_surface_field=.false.
    
  end function has_scalar_surface_field
  
  integer function get_scalar_boundary_condition_count(field)
  !!< Get number of boundary conditions of a scalar field
  type(scalar_field), intent(in):: field
    
    if (associated(field%bc%boundary_condition)) then
      get_scalar_boundary_condition_count=size(field%bc%boundary_condition)
    else
      get_scalar_boundary_condition_count=0
    end if
  
  end function get_scalar_boundary_condition_count
    
  integer function get_vector_boundary_condition_count(field)
  !!< Get number of boundary conditions of a scalar field
  type(vector_field), intent(in):: field
    
    if (associated(field%bc%boundary_condition)) then
      get_vector_boundary_condition_count=size(field%bc%boundary_condition)
    else
      get_vector_boundary_condition_count=0
    end if
  
  end function get_vector_boundary_condition_count

  subroutine get_scalar_boundary_condition_by_number(field, n, name, type, &
    surface_element_list, surface_node_list, surface_mesh, &
    option_path)
  !!< Get boundary condition of a scalar field
  type(scalar_field), intent(in):: field
  !! which boundary condition
  integer, intent(in):: n
  !! name of the b.c.
  character(len=*), intent(out), optional:: name
  !! type of b.c., any of: ...
  character(len=*), intent(out), optional:: type
  !! pointer to list of surface elements where this b.c. is applied
  integer, dimension(:), pointer, optional:: surface_element_list
  !! pointer to list of surface nodes where this b.c. is applied
  integer, dimension(:), pointer, optional:: surface_node_list
  !! surface mesh on which surface fields can be allocated
  type(mesh_type), pointer, optional:: surface_mesh
  !! option_path for the bc
  character(len=*), intent(out), optional:: option_path
  
    type(scalar_boundary_condition), pointer:: bc

    assert(n>=1 .and. n<=size(field%bc%boundary_condition))
    bc => field%bc%boundary_condition(n)
    
    if (present(name)) then
      name=bc%name
    end if
    
    if (present(type)) then
      type=bc%type
    end if
      
    if (present(surface_element_list)) then
      surface_element_list => bc%surface_element_list
    end if
    
    if (present(surface_node_list)) then
      surface_node_list => bc%surface_node_list
    end if
    
    if (present(surface_mesh)) then
      surface_mesh => bc%surface_mesh
    end if
    
    if (present(option_path)) then
      option_path = bc%option_path
    end if
    
  end subroutine get_scalar_boundary_condition_by_number

  subroutine get_vector_boundary_condition_by_number(field, n, name, type, &
    surface_element_list, surface_node_list, applies, surface_mesh, &
    option_path)
  !!< Get boundary condition of a vector field
  type(vector_field), intent(in):: field
  !! which boundary condition
  integer, intent(in):: n
  !! name of the b.c.
  character(len=*), intent(out), optional:: name
  !! type of b.c., any of: ...
  character(len=*), intent(out), optional:: type
  !! pointer to list of surface elements where this b.c. is applied
  integer, dimension(:), pointer, optional:: surface_element_list
  !! pointer to list of surface nodes where this b.c. is applied
  integer, dimension(:), pointer, optional:: surface_node_list
  !! vector components to which this b.c. applies
  logical, dimension(:), intent(out), optional:: applies
  !! surface mesh on which surface fields can be allocated
  type(mesh_type), pointer, optional:: surface_mesh
  !! option_path for the bc
  character(len=*), intent(out), optional:: option_path
  
    type(vector_boundary_condition), pointer:: bc

    assert(n>=1 .and. n<=size(field%bc%boundary_condition))
    bc => field%bc%boundary_condition(n)
    
    if (present(name)) then
      name=bc%name
    end if
    
    if (present(type)) then
      type=bc%type
    end if
    
    if (present(surface_element_list)) then
      surface_element_list => bc%surface_element_list
    end if
    
    if (present(surface_node_list)) then
      surface_node_list => bc%surface_node_list
    end if
    
    if (present(applies)) then
      applies=bc%applies(1:size(applies))
    end if
        
    if (present(surface_mesh)) then
      surface_mesh => bc%surface_mesh
    end if
    
    if (present(option_path)) then
      option_path = bc%option_path
    end if
    
  end subroutine get_vector_boundary_condition_by_number

  subroutine get_scalar_boundary_condition_by_name(field, name, &
    type, surface_node_list, surface_element_list, surface_mesh, &
    option_path)
  !!< Get boundary condition of a scalar field
  type(scalar_field), intent(in):: field
  !! which boundary condition
  character(len=*), intent(in):: name
  !! type of b.c., any of: ...
  character(len=*), intent(out), optional:: type
  !! pointer to list of surface elements where this b.c. is applied
  integer, dimension(:), pointer, optional:: surface_element_list
  !! pointer to list of surface nodes where this b.c. is applied
  integer, dimension(:), pointer, optional:: surface_node_list
  !! surface mesh on which surface fields can be allocated
  type(mesh_type), pointer, optional:: surface_mesh
  !! option_path for the bc
  character(len=*), intent(out), optional:: option_path
  
    integer i
    do i=1, size(field%bc%boundary_condition)
      if (field%bc%boundary_condition(i)%name==name) then
        call get_scalar_boundary_condition_by_number(field, i, &
            type=type, surface_element_list=surface_element_list, &
              surface_node_list=surface_node_list, &
              surface_mesh=surface_mesh, &
              option_path=option_path)
        return
      end if
    end do
    ewrite(-1,*) 'Unknown boundary condition: ', name
    FLAbort("Hasta la vista")
  
  end subroutine get_scalar_boundary_condition_by_name

  subroutine get_vector_boundary_condition_by_name(field, name, &
    type, surface_element_list, surface_node_list, applies, surface_mesh, &
    option_path)
  !!< Get boundary condition of a vector field
  type(vector_field), intent(in):: field
  !! which boundary condition
  character(len=*), intent(in):: name
  !! type of b.c., any of: ...
  character(len=*), intent(out), optional:: type
  !! pointer to list of surface elements where this b.c. is applied
  integer, dimension(:), pointer, optional:: surface_element_list
  !! pointer to list of surface nodes where this b.c. is applied
  integer, dimension(:), pointer, optional:: surface_node_list
  !! vector components to which this b.c. applies
  logical, dimension(:), intent(out), optional:: applies
  !! surface mesh on which surface fields can be allocated
  type(mesh_type), pointer, optional:: surface_mesh
  !! option_path for the bc
  character(len=*), intent(out), optional:: option_path
  
    integer i
    do i=1, size(field%bc%boundary_condition)
      if (field%bc%boundary_condition(i)%name==name) then
        call get_vector_boundary_condition_by_number(field, i, &
            type, surface_element_list=surface_element_list, &
              surface_node_list=surface_node_list, &
              applies=applies, surface_mesh=surface_mesh, &
              option_path=option_path)
        return
      end if
    end do
    ewrite(-1,*) 'Unknown boundary condition: ', name
    FLAbort("Hasta la vista")
  
  end subroutine get_vector_boundary_condition_by_name

  subroutine get_periodic_boundary_condition(mesh, periodic_bc_list)
    !!< Gets a list of the surface elements which are periodic
    
    !! mesh on which periodicness is to be evaluated
    type(mesh_type), intent(inout):: mesh
    !! For each surface_element returns true for periodic boundaries
    !! or false, if not:
    logical, dimension(:), intent(out):: periodic_bc_list
    
    integer sele
    
    integer, dimension(:), pointer :: neigh
    integer :: l_face_number, ele

    periodic_bc_list = .false.

    do sele = 1, surface_element_count(mesh)
      l_face_number = local_face_number(mesh, sele)
      ele=face_ele(mesh, sele)
      neigh => ele_neigh(mesh, ele)
      
      if(neigh(l_face_number)>0) then
        
        periodic_bc_list(sele)=.true.
        
      end if
    end do

  end subroutine get_periodic_boundary_condition
    
  subroutine get_entire_scalar_boundary_condition(field, &
     types, boundary_value, bc_type_list)
    !!< Gets the boundary conditions on the entire surface mesh for all
    !!< bc types requested
    
    !! field of which boundary conditions are retrieved
    type(scalar_field), intent(inout):: field
    !! list of bc types you want (others are ignored)
    character(len=*), dimension(:), intent(in):: types
    !! A field over the entire surface containing the boundary values
    !! for the bcs of the type requested. This field is defined on a 
    !! dg surface mesh so that it can deal with discontinuities between
    !! differen boundary conditions.
    !! The ordering of the (surface) elements in this mesh is the same
    !! as the ordering of the surface elements (faces) of the given
    !! field.
    !! This field should be deallocated after use.
    type(scalar_field), intent(out):: boundary_value
    !! For each surface_element returns the position in the types argument list,
    !! thus identifying the applied boundary condition type,
    !! or zero, if no bc of the requested types are applied to this face:
    integer, dimension(:), intent(out):: bc_type_list
    
    type(scalar_field), pointer:: surface_field
    type(mesh_type), pointer:: surface_mesh
    character(len=FIELD_NAME_LEN) bctype
    character(len=1024) name
    integer, dimension(:), pointer:: surface_element_list
    integer i, j, k, sele
    
    integer, dimension(:), pointer :: neigh
    integer :: l_face_number, ele

    surface_mesh => get_dg_surface_mesh(field%mesh)

    call allocate(boundary_value, surface_mesh, name=trim(field%name)//"EntireBC")
    call zero(boundary_value)
    bc_type_list=0
    
    do i=1, get_boundary_condition_count(field)
       call get_boundary_condition(field, i, type=bctype, &
          surface_element_list=surface_element_list,name=name)
          
       ! see if we're interested in this one, if not skip it
       do j=1, size(types)
          if (trim(types(j))==trim(bctype)) exit
       end do
       if (j>size(types)) cycle
       
       if (associated(field%bc%boundary_condition(i)%surface_fields)) then
          ! extract 1st surface field
          surface_field => field%bc%boundary_condition(i)%surface_fields(1)
       else
          nullify(surface_field)
       end if       
       
       do k=1, size(surface_element_list)
          sele=surface_element_list(k)
          
          if (bc_type_list(sele)/=0) then
             ewrite(0,*) 'Requested types:', types
             ewrite(0,*) 'Of these boundary condition types only one may be applied'
             ewrite(0,*) 'to each surface element.'
             ewrite(0,*) 'Surface element nr.:', sele
             ewrite(0,*) 'has types', types(bc_type_list(sele)), bctype
             ewrite(0,*) 'on field: ', field%name
             ewrite(0,*) 'with name: ',name
             FLAbort("Can't have that.")
          end if
          bc_type_list(sele)=j
          
          if (associated(surface_field)) then
             call set(boundary_value, ele_nodes(surface_mesh, sele), &
                ele_val(surface_field, k))
          end if
          
       end do
    end do
    
    do j=1, size(types)
      if (trim(types(j))=="periodic") exit
    end do
    if (j<=size(types)) then
      do sele = 1, surface_element_count(field)
        ele = face_ele(field, sele)
        l_face_number = local_face_number(field, sele)
        neigh => ele_neigh(field, ele)
        
        if(neigh(l_face_number)>0) then
          
          if (bc_type_list(sele)/=0) then
             ewrite(0,*) 'Requested types:', types
             ewrite(0,*) 'Of these boundary condition types only one may be applied'
             ewrite(0,*) 'to each surface element.'
             ewrite(0,*) 'Surface element nr.:', sele
             ewrite(0,*) 'has types', types(bc_type_list(sele)), 'periodic'
             FLAbort("Can't have that.")
          end if
          
          bc_type_list(sele)=j
          
        end if
      end do
    end if

  end subroutine get_entire_scalar_boundary_condition
    
  subroutine get_entire_vector_boundary_condition(field, &
     types, boundary_value, bc_type_list)
    !!< Gets the boundary conditions on the entire surface mesh for all
    !!< bc types requested
    
    !! field of which boundary conditions are retrieved
    type(vector_field), intent(inout):: field
    !! list of bc types you want (others are ignored)
    character(len=*), dimension(:), intent(in):: types
    !! A field over the entire surface containing the boundary values
    !! for the bcs of the type requested. This field is defined on a 
    !! dg surface mesh so that it can deal with discontinuities between
    !! differen boundary conditions.
    !! The ordering of the (surface) elements in this mesh is the same
    !! as the ordering of the surface elements (faces) of the given
    !! field.
    !! This field should be deallocated after use.
    type(vector_field), intent(out):: boundary_value
    !! For each surface_element returns the position in the types argument list,
    !! thus identifying the applied boundary condition type,
    !! or zero, if no bc of the requested types are applied to this face.
    !! BC can be set for each component separately, so ndim x surface_element_count()
    integer, dimension(:,:), intent(out):: bc_type_list
    
    type(vector_field), pointer:: surface_field
    type(mesh_type), pointer:: surface_mesh
    character(len=FIELD_NAME_LEN) bctype
    integer, dimension(:), pointer:: surface_element_list
    logical, dimension(field%dim):: applies
    integer i, j, k, n, sele

    integer, dimension(:), pointer :: neigh
    integer :: l_face_number, ele

    surface_mesh => get_dg_surface_mesh(field%mesh)

    call allocate(boundary_value, field%dim, surface_mesh, &
       name=trim(field%name)//"EntireBC")
    call zero(boundary_value)
    bc_type_list=0
    
    bcloop: do i=1, get_boundary_condition_count(field)
       call get_boundary_condition(field, i, type=bctype, applies=applies, &
          surface_element_list=surface_element_list)
          
       ! see if we're interested in this one, if not skip it
       do j=1, size(types)
          if (trim(types(j))==trim(bctype)) exit
       end do
       if (j>size(types)) cycle
       
       if (associated(field%bc%boundary_condition(i)%surface_fields)) then
          ! extract 1st surface field
          surface_field => field%bc%boundary_condition(i)%surface_fields(1)
       else
          nullify(surface_field)
       end if       
       
       faceloop: do k=1, size(surface_element_list)
          sele=surface_element_list(k)
          
          do n=1, field%dim
             if (applies(n)) then
                if (bc_type_list(n, sele)/=0) then
                   ewrite(0,*) 'Requested types:', types
                   ewrite(0,*) 'Of these boundary condition types'
                   ewrite(0,*) 'only one may be applied'
                   ewrite(0,*) 'to each surface element.'
                   ewrite(0,*) 'Surface element nr.:', sele
                   ewrite(0,*) 'Component nr:', n
                   ewrite(0,*) 'has types', types(bc_type_list(j, sele)), bctype
                   FLAbort("Can't have that.")
                end if
                bc_type_list(n, sele)=j
          
                if (associated(surface_field)) then
                   call set(boundary_value, n, ele_nodes(surface_mesh, sele), &
                      ele_val(surface_field, n, k))
                end if
             end if
          end do
          
       end do faceloop
    end do bcloop

    do j=1, size(types)
      if (trim(types(j))=="periodic") exit
    end do
    if (j<=size(types)) then
      do sele = 1, surface_element_count(field)
        ele = face_ele(field, sele)
        l_face_number = local_face_number(field, sele)
        neigh => ele_neigh(field, ele)
        
        if(neigh(l_face_number)>0) then
          
          if (any(bc_type_list(:,sele)/=0)) then
             ewrite(0,*) 'Requested types:', types
             ewrite(0,*) 'Of these boundary condition types only one may be applied'
             ewrite(0,*) 'to each surface element.'
             ewrite(0,*) 'Surface element nr.:', sele
             ewrite(0,*) 'has types', types(bc_type_list(:,sele)), 'periodic'
             FLAbort("Can't have that.")
          end if
          
          bc_type_list(:,sele)=j
          
        end if
      end do
    end if

  end subroutine get_entire_vector_boundary_condition
    
  function get_dg_surface_mesh(mesh)
  !! Returns a pointer a DG version of the entire surface mesh
  !! This DG surface mesh can be used to store boundary condition values on
  !! so it can deal with discontinuities in boundary condition values.
  !! If this dg surface mesh does not yet exist, it's created here. It
  !! will be cached on the mesh for later calls to this routine and only
  !! be deallocated upon deallocation of the entire mesh.
  !!
  !! Even in the case that the mesh is DG itself, it may still require a new dg surface mesh
  !! as the original surface mesh may not be dg at corners (which would prevent different
  !! bcs from being applied on elements either side of the corner).  Therefore, to be
  !! safe a new one is created.
  type(mesh_type), intent(inout):: mesh
  type(mesh_type), pointer:: get_dg_surface_mesh
  
    if (.not. has_faces(mesh)) then
      FLAbort("A mesh%faces structure is needed for get_dg_surface_mesh")
    end if
    
    if (.not. associated(mesh%faces%dg_surface_mesh)) then
      allocate(mesh%faces%dg_surface_mesh)
      mesh%faces%dg_surface_mesh=make_mesh(mesh%faces%surface_mesh, mesh%faces%shape, &
        continuity=-1, name=trim(mesh%name)//"DGSurfaceMesh")
    end if
    
    get_dg_surface_mesh => mesh%faces%dg_surface_mesh
    
  end function get_dg_surface_mesh
  
  subroutine get_scalar_boundary_condition_nodes(field, types, bc_type_node_list)
  ! gets the boundary conditions on the entire surface mesh for all
  ! bc types requested
     ! field of which boundary conditions are retrieved
     type(scalar_field), intent(in):: field
     ! bc type you want (others are ignored)
     character(len=*), dimension(:), intent(in):: types
     ! list of nodes on which boundary condition is applied
     integer, dimension(:), intent(out):: bc_type_node_list
     
     integer, dimension(:), pointer:: surface_node_list
     integer i, j, l_face_number, sele, ele
     integer, dimension(:), pointer :: neigh
     
     bc_type_node_list=0

     do i=1, get_boundary_condition_count(field)
        do j=1, size(types)
           if (types(j)==field%bc%boundary_condition(i)%type) exit
        end do

        if (j>size(types)) cycle

        call get_boundary_condition(field, i, &
           surface_node_list=surface_node_list)

        if(any((bc_type_node_list(surface_node_list)/=0) .and. &
               (bc_type_node_list(surface_node_list)/=j))) then
            ewrite(0,*) 'Requested types:', types
            ewrite(0,*) 'Of these boundary condition types only one'
            ewrite(0,*) 'may be applied to each node.'
            FLAbort("Sorry!")
        end if

        bc_type_node_list(surface_node_list) = j

     end do
      
      do j=1, size(types)
        if (trim(types(j))=="periodic") exit
      end do
      if (j<=size(types)) then
        do sele = 1, surface_element_count(field)
          ele = face_ele(field, sele)
          l_face_number = local_face_number(field, sele)
          neigh => ele_neigh(field, ele)
          
          if(neigh(l_face_number)>0) then
            
            if (any((bc_type_node_list(face_global_nodes(field, sele))/=0).and.&
                    (bc_type_node_list(face_global_nodes(field, sele))/=j))) then
                ewrite(0,*) 'Requested types:', types
                ewrite(0,*) 'Of these boundary condition types only one'
                ewrite(0,*) 'may be applied to each node.'
                FLAbort("Sorry!")
            end if
            
            bc_type_node_list(face_global_nodes(field, sele))=j
            
          end if
        end do
      end if

  end subroutine get_scalar_boundary_condition_nodes

  subroutine get_vector_boundary_condition_nodes(field, types, bc_type_node_list)
  ! gets the boundary conditions on the entire surface mesh for all
  ! bc types requested
     ! field of which boundary conditions are retrieved
     type(vector_field), intent(in):: field
     ! bc type you want (others are ignored)
     character(len=*), dimension(:), intent(in):: types
     ! list of nodes on which boundary condition is applied
     integer, dimension(:,:), intent(out):: bc_type_node_list
     
     logical, dimension(field%dim):: applies
     integer, dimension(:), pointer:: surface_node_list
     integer i, j, k, ele, sele, l_face_number
     integer, dimension(:), pointer :: neigh
     
     bc_type_node_list=0

     do i=1, get_boundary_condition_count(field)
        do j=1, size(types)
           if (types(j)==field%bc%boundary_condition(i)%type) exit
        end do

        if (j>size(types)) cycle

        call get_boundary_condition(field, i, &
           surface_node_list=surface_node_list, &
           applies=applies)

        do k = 1, field%dim
          if(applies(k)) then
            if(any((bc_type_node_list(k,surface_node_list)/=0) .and. &
                   (bc_type_node_list(k,surface_node_list)/=j))) then
                ewrite(0,*) 'Requested types:', types
                ewrite(0,*) 'Of these boundary condition types only one'
                ewrite(0,*) 'may be applied to each node.'
                FLAbort("Sorry!")
            end if

            bc_type_node_list(k,surface_node_list) = j
          end if
        end do

     end do
      
      do j=1, size(types)
        if (trim(types(j))=="periodic") exit
      end do
      if (j<=size(types)) then
        do sele = 1, surface_element_count(field)
          ele = face_ele(field, sele)
          l_face_number = local_face_number(field, sele)
          neigh => ele_neigh(field, ele)
          
          if(neigh(l_face_number)>0) then
            
            if (any((bc_type_node_list(:,face_global_nodes(field, sele))/=0).and.&
                    (bc_type_node_list(:,face_global_nodes(field, sele))/=j))) then
                ewrite(0,*) 'Requested types:', types
                ewrite(0,*) 'Of these boundary condition types only one'
                ewrite(0,*) 'may be applied to each node.'
                FLAbort("Sorry!")
            end if
            
            bc_type_node_list(:,face_global_nodes(field, sele))=j
            
          end if
        end do
      end if
    
  end subroutine get_vector_boundary_condition_nodes

  subroutine set_reference_node_scalar(matrix, node, rhs, reference_value)
    !!< Sets a reference node for which the value is fixed in the equation
    !!< This is typically done for a Poisson equation with all Neumann
    !!< bcs to eliminate the spurious freedom of adding a constant value
    !!< to the solution.
    type(csr_matrix), intent(inout) :: matrix
    integer, intent(in):: node
    !! if rhs is not provided, you have to make sure the rhs at 
    !! the reference node has the right value, usually 0, yourself:
    type(scalar_field), optional, intent(inout) :: rhs
    !! by default the field gets set to 0 at the reference node
    real, optional, intent(in) :: reference_value
    
    ! but only first should set the reference node
    if (GetProcNo()/=1) then
       ! Other processors still need to have the inactive mask, even if
       ! it's empty.
       call initialise_inactive(matrix)
       return
    end if
    
    call set_inactive(matrix, node)
    
    if (present(rhs)) then
       if (present(reference_value)) then
          call set(rhs, node, reference_value)
       else
          call set(rhs, node, 0.0)
       end if
    end if
    
  end subroutine set_reference_node_scalar

  logical function has_boundary_condition_scalar(field, type)
  !!< logical function that tells whether any of the bcs of a field
  !!< is of type 'type'
  type(scalar_field), intent(in):: field
  character(len=*), intent(in):: type
  
    type(scalar_boundary_condition), pointer :: this_bc
    integer i
    
    if (.not.associated(field%bc%boundary_condition)) then
      has_boundary_condition_scalar=.false.
      return
    end if
    
    bcloop: do i=1, size(field%bc%boundary_condition)
       this_bc=>field%bc%boundary_condition(i)

       if (this_bc%type==type) then
         
         has_boundary_condition_scalar=.true.
         return
         
       end if

    end do bcloop
      
    has_boundary_condition_scalar=.false.
      
  end function has_boundary_condition_scalar

  logical function has_boundary_condition_vector(field, type)
  !!< logical function that tells whether any of the bcs of a field
  !!< is of type 'type'
  type(vector_field), intent(in):: field
  character(len=*), intent(in):: type
  
    type(vector_boundary_condition), pointer :: this_bc
    integer i
    
    if (.not.associated(field%bc%boundary_condition)) then
      has_boundary_condition_vector=.false.
      return
    end if
    
    bcloop: do i=1, size(field%bc%boundary_condition)
       this_bc=>field%bc%boundary_condition(i)

       if (this_bc%type==type) then
         
         has_boundary_condition_vector=.true.
         return
         
       end if

    end do bcloop
      
    has_boundary_condition_vector=.false.
      
  end function has_boundary_condition_vector

  logical function has_boundary_condition_name_scalar(field, name)
  !!< logical function that tells whether any of the bcs of a field
  !!< has the name 'name'
  type(scalar_field), intent(in):: field
  character(len=*), intent(in):: name
  
    type(scalar_boundary_condition), pointer :: this_bc
    integer i
    
    if (.not.associated(field%bc%boundary_condition)) then
      has_boundary_condition_name_scalar=.false.
      return
    end if
    
    bcloop: do i=1, size(field%bc%boundary_condition)
       this_bc=>field%bc%boundary_condition(i)

       if (this_bc%name==name) then
         
         has_boundary_condition_name_scalar=.true.
         return
         
       end if

    end do bcloop
      
    has_boundary_condition_name_scalar=.false.
      
  end function has_boundary_condition_name_scalar

  logical function has_boundary_condition_name_vector(field, name)
  !!< logical function that tells whether any of the bcs of a field
  !!< has the name 'name'
  type(vector_field), intent(in):: field
  character(len=*), intent(in):: name
  
    type(vector_boundary_condition), pointer :: this_bc
    integer i
    
    if (.not.associated(field%bc%boundary_condition)) then
      has_boundary_condition_name_vector=.false.
      return
    end if
    
    bcloop: do i=1, size(field%bc%boundary_condition)
       this_bc=>field%bc%boundary_condition(i)

       if (this_bc%name==name) then
         
         has_boundary_condition_name_vector=.true.
         return
         
       end if

    end do bcloop
      
    has_boundary_condition_name_vector=.false.
      
  end function has_boundary_condition_name_vector
  
  subroutine set_dirichlet_consistent(states)
    !!< Once the fields and boundary conditions have been set, force the
    !!< boundary values of the fields to be consistent with any Dirichlet
    !!< conditions specified.
    type(state_type), dimension(:), intent(in):: states

    type(scalar_field), pointer:: sfield
    type(vector_field), pointer:: vfield

    character(len=OPTION_PATH_LEN) field_path
    integer :: p, nphases, f, nfields

    ewrite(1,*) "In set_dirichlet_consistent"

    nphases = size(states)
    do p = 0, nphases-1

       ! Scalar fields:
       nfields = scalar_field_count(states(p+1))
       do f = 1, nfields
          sfield => extract_scalar_field(states(p+1),f)
          field_path=sfield%option_path
          if (.not. have_option(trim(field_path)//'/prognostic')) cycle

          ! only prognostic fields from here:
          call set_dirichlet_consistent_scalar(sfield)

       end do

       ! Vector fields:

       nfields = vector_field_count(states(p+1))
       do f = 1, nfields
          vfield => extract_vector_field(states(p+1), f)
          field_path=vfield%option_path
          if (.not. have_option(trim(field_path)//'/prognostic')) cycle

          ! only prognostic fields from here:
          call set_dirichlet_consistent_vector(vfield)

       end do

    end do

  end subroutine set_dirichlet_consistent

  subroutine set_dirichlet_consistent_scalar(field)
    !!< Force the values of the boundary nodes of a scalar field to the
    !!< dirichlet boundary condition values.
    type(scalar_field), intent(inout) :: field

    type(scalar_field), pointer:: surface_field
    integer, dimension(:), pointer:: surface_node_list
    character(len=OPTION_PATH_LEN):: bc_option_path
    character(len=FIELD_NAME_LEN):: bc_type
    integer :: i, b

    do b=1, get_boundary_condition_count(field)
       call get_boundary_condition(field, b, &
          type=bc_type, option_path=bc_option_path, &
          surface_node_list=surface_node_list)
       

       if ((bc_type/="dirichlet").and.(bc_type/="weakdirichlet")) cycle
       
       ! Weak dirichlet bcs do not have consistent bcs by default.
       if (bc_type=="weakdirichlet".and. &
            .not.have_option(trim(bc_option_path)//&
            "/type[0]/apply_weakly/boundary_overwrites_initial_condition"))&
            & then
          cycle
       end if
       
       surface_field => extract_surface_field(field, b, "value")

       do i=1,size(surface_node_list)
          call set(field, surface_node_list(i), &
            node_val(surface_field, i))
       end do

    end do

  end subroutine set_dirichlet_consistent_scalar

  subroutine set_dirichlet_consistent_vector(field)
    !!< Force the values of the boundary nodes of a vector field to the
    !!< dirichlet boundary condition values.
    type(vector_field), intent(inout) :: field

    type(vector_field), pointer:: surface_field, normal, tangent_1, tangent_2
    integer, dimension(:), pointer:: surface_node_list
    character(len=OPTION_PATH_LEN):: bc_option_path
    character(len=FIELD_NAME_LEN):: bc_type
    real, dimension(1:field%dim, 1:field%dim):: rotation_mat
    real, dimension(1:field%dim):: rotated_vector
    logical, dimension(1:field%dim):: applies
    logical:: rotated
    integer :: i, b, d
    
    do b=1, get_boundary_condition_count(field)

       call get_boundary_condition(field, b, &
          type=bc_type, option_path=bc_option_path, &
          surface_node_list=surface_node_list, &
          applies=applies)
       

       if ((bc_type/="dirichlet").and.(bc_type/="weakdirichlet")) cycle
       
       ! Weak dirichlet bcs do not have consistent bcs by default.
       if (bc_type=="weakdirichlet".and. &
            .not.have_option(trim(bc_option_path)//&
            "/type[0]/apply_weakly/boundary_overwrites_initial_condition"))&
            & then
          cycle
       end if
       
       surface_field => extract_surface_field(field, b, "value")

       rotated=have_option(trim(bc_option_path)//'/type[0]/align_bc_with_surface')
       if (rotated) then
         
         normal => extract_surface_field(field, b, "normal")
         tangent_1 => extract_surface_field(field, b, "tangent1")
         tangent_2 => extract_surface_field(field, b, "tangent2")
         
         do i=1, size(surface_node_list)
           rotation_mat(:,1)=node_val(normal, i)
           if (field%dim>1) then
             rotation_mat(:,2)=node_val(tangent_1, i)
           end if
           if (field%dim>2) then
             rotation_mat(:,3)=node_val(tangent_2, i)
           end if
           ! first we rotate the existing vector into (normal, tangent1/2) coordinates
           rotated_vector=matmul( node_val(field, surface_node_list(i)), &
             rotation_mat )
           ! overwrite those components where the bc is applied
           do d=1, field%dim
             if (applies(d)) then
               rotated_vector(d)=node_val(surface_field, d, i)
             end if
           end do
           ! and rotate back to x,y,z coordinates
           call set(field, surface_node_list(i), &
              matmul( rotation_mat, rotated_vector ))
         end do
           
       else
       
         ! non-rotated, cartesian aligned case
       
         ! Loop over possible dimensions.
         do d=1, field%dim
            if (.not. applies(d)) cycle

            do i=1, size(surface_node_list)
               call set(field, d, surface_node_list(i), &
                  node_val(surface_field, d, i))
            end do
         end do
           
       end if
       
    end do

  end subroutine set_dirichlet_consistent_vector

  subroutine apply_dirichlet_conditions_scalar(matrix, rhs, field, dt)
    !!< Apply dirichlet boundary conditions from field to the problem
    !!< defined by matrix and rhs.
    !!<
    !!< If dt is supplied, this assumes that boundary 
    !!< conditions are applied in rate of change form.
    type(csr_matrix), intent(inout) :: matrix
    type(scalar_field), intent(inout) :: rhs
    type(scalar_field), optional, intent(in) :: field
    real, optional, intent(in) :: dt
    
    type(scalar_field), pointer:: surface_field
    integer, dimension(:), pointer:: surface_node_list
    character(len=FIELD_NAME_LEN):: bctype
    integer :: i,j

    bcloop: do i=1, get_boundary_condition_count(field)
       call get_boundary_condition(field, i, type=bctype, &
          surface_node_list=surface_node_list)

       if (bctype/="dirichlet") cycle bcloop
       
       call set_inactive(matrix, surface_node_list)
       
       surface_field => extract_surface_field(field, i, "value")
       
       if (present(dt)) then
          do j=1, size(surface_node_list)
            call set(rhs, surface_node_list(j), &
                 (node_val(surface_field, j)- &
                  node_val(field, surface_node_list(j)) &
                 )/dt)
          end do
       else
          do j=1, size(surface_node_list)
            call set(rhs, surface_node_list(j), &
                 node_val(surface_field, j))
          end do
       end if
       
    end do bcloop
    
  end subroutine apply_dirichlet_conditions_scalar

  subroutine apply_dirichlet_conditions_scalar_lumped(lhs, rhs, field, dt)
    !!< Apply dirichlet boundary conditions from field to the problem
    !!< defined by lhs and rhs.  lhs is the diagonal of the normal matrix
    !!< so this is only useful for fully explicit problems.
    !!<
    !!< This assumes that boundary conditions are applied in rate of change
    !!< form.
    type(scalar_field), intent(inout) :: lhs
    type(scalar_field), intent(inout) :: rhs
    type(scalar_field), intent(in) :: field
    real, intent(in) :: dt
    
    type(scalar_field), pointer:: surface_field
    integer, dimension(:), pointer:: surface_node_list
    character(len=FIELD_NAME_LEN):: bctype
    integer :: i,j

    bcloop: do i=1, get_boundary_condition_count(field)
       call get_boundary_condition(field, i, type=bctype, &
          surface_node_list=surface_node_list)

       if (bctype/="dirichlet") cycle bcloop
       
       surface_field => extract_surface_field(field, i, "value")
       
       do j=1,size(surface_node_list)
          call addto(rhs, surface_node_list(j), &
               ((node_val(surface_field, j)- &
                node_val(field, surface_node_list(j)) &
                ) /dt)*INFINITY)

          call addto(lhs, surface_node_list(j), &
               INFINITY)

       end do
       
    end do bcloop
    
  end subroutine apply_dirichlet_conditions_scalar_lumped

  subroutine apply_dirichlet_conditions_vector(matrix, rhs, field, dt)
    !!< Apply dirichlet boundary conditions from field to the problem
    !!< defined by matrix and rhs.
    !!<
    !!< This assumes that boundary conditions are applied in rate of change
    !!< form and that the matrix has dim x dim blocks.
    type(block_csr_matrix), intent(inout) :: matrix
    type(vector_field), intent(inout), optional :: rhs
    type(vector_field), intent(in) :: field
    real, intent(in), optional :: dt

    type(scalar_field) :: rhscomponent, bccomponent
    type(csr_matrix) :: matrixcomponent

    logical, dimension(field%dim):: applies
    character(len=FIELD_NAME_LEN):: bctype
    type(vector_field), pointer:: surface_field
    integer, dimension(:), pointer:: surface_node_list
    integer :: i,j,k

    bcloop: do i=1, get_boundary_condition_count(field)
       call get_boundary_condition(field, i, type=bctype, &
          surface_node_list=surface_node_list, applies=applies)

       if (bctype/="dirichlet") cycle bcloop
       
       surface_field => extract_surface_field(field, i, "value")
       
       do j=1,size(surface_node_list)
          do k = 1, field%dim
            if(applies(k)) then

              if(present(rhs).and.present(dt)) then
                rhscomponent = extract_scalar_field_from_vector_field(rhs, k)
                bccomponent = extract_scalar_field_from_vector_field(surface_field, k)

                call addto(rhscomponent, &
                           surface_node_list(j), &
                           ((node_val(bccomponent,j)&
                             -node_val(field, k, surface_node_list(j)) &
                             ) /dt)*INFINITY)
              end if

              matrixcomponent = block(matrix, k, k)
              call addto_diag(matrixcomponent, surface_node_list(j), &
                              INFINITY)
            end if
          end do
       end do

    end do bcloop

  end subroutine apply_dirichlet_conditions_vector

  subroutine apply_dirichlet_conditions_vector_petsc_csr(matrix, rhs, field, dt)
    !!< Apply dirichlet boundary conditions from field to the problem
    !!< defined by matrix and rhs.
    !!<
    !!< This assumes that boundary conditions are applied in rate of change
    !!< form and that the matrix has dim x dim blocks.
    type(petsc_csr_matrix), intent(inout) :: matrix
    type(vector_field), intent(inout), optional :: rhs
    type(vector_field), intent(in) :: field
    real, intent(in), optional :: dt

    type(scalar_field) :: rhscomponent, bccomponent

    logical, dimension(field%dim):: applies
    character(len=FIELD_NAME_LEN):: bctype
    type(vector_field), pointer:: surface_field
    integer, dimension(:), pointer:: surface_node_list
    integer :: i,j,k

    bcloop: do i=1, get_boundary_condition_count(field)
       call get_boundary_condition(field, i, type=bctype, &
          surface_node_list=surface_node_list, applies=applies)

       if (bctype/="dirichlet") cycle bcloop
       
       surface_field => extract_surface_field(field, i, "value")
       
       do j=1,size(surface_node_list)
          do k = 1, field%dim
            if(applies(k)) then

              if(present(rhs).and.present(dt)) then
                rhscomponent = extract_scalar_field_from_vector_field(rhs, k)
                bccomponent = extract_scalar_field_from_vector_field(surface_field, k)

                call addto(rhscomponent, &
                           surface_node_list(j), &
                           ((node_val(bccomponent,j)&
                             -node_val(field, k, surface_node_list(j)) &
                             ) /dt)*INFINITY)
              end if

              call addto(matrix, k, k, surface_node_list(j), &
                surface_node_list(j), INFINITY)
            end if
          end do
       end do

    end do bcloop

  end subroutine apply_dirichlet_conditions_vector_petsc_csr
  
  subroutine apply_dirichlet_conditions_vector_component(matrix, rhs, field, dt, dim)
    !!< Apply dirichlet boundary conditions from field to the problem
    !!< defined by matrix and rhs.
    !!<
    !!< This assumes that boundary conditions are applied in rate of change
    !!< form and that the matrix has dim x dim blocks.
    type(csr_matrix), intent(inout) :: matrix
    type(vector_field), intent(inout), optional :: rhs
    type(vector_field), intent(in) :: field
    real, intent(in), optional :: dt
    integer, intent(in) :: dim

    type(scalar_field) :: rhscomponent, bccomponent

    logical, dimension(field%dim):: applies
    character(len=FIELD_NAME_LEN):: bctype
    type(vector_field), pointer:: surface_field
    integer, dimension(:), pointer:: surface_node_list
    integer :: i,j

    bcloop: do i=1, get_boundary_condition_count(field)
       call get_boundary_condition(field, i, type=bctype, &
          surface_node_list=surface_node_list, applies=applies)

       if (bctype/="dirichlet") cycle bcloop
       
       surface_field => extract_surface_field(field, i, "value")
       
       do j=1,size(surface_node_list)
         
          if(applies(dim)) then
            
            call addto_diag(matrix, surface_node_list(j),&
                            INFINITY)

            if(present(rhs).and.present(dt)) then
              rhscomponent = extract_scalar_field_from_vector_field(rhs, dim)
              bccomponent = extract_scalar_field_from_vector_field(surface_field, dim)

              call addto(rhscomponent, &
                          surface_node_list(j), &
                          ((node_val(bccomponent, j) &
                            -node_val(field,dim,surface_node_list(j)) &
                           ) /dt)*INFINITY)

            end if

          end if
       end do

    end do bcloop

  end subroutine apply_dirichlet_conditions_vector_component

  subroutine apply_dirichlet_conditions_vector_lumped(lhs, rhs, field, dt)
    !!< Apply dirichlet boundary conditions from field to the problem
    !!< defined by lhs and rhs.  lhs is the diagonal of the normal matrix.
    !!<
    !!< This assumes that boundary conditions are applied in rate of change
    !!< form.
    type(vector_field), intent(inout) :: lhs
    type(vector_field), intent(inout), optional :: rhs
    type(vector_field), intent(in) :: field
    real, intent(in), optional :: dt

    type(scalar_field) :: rhscomponent, bccomponent, lhscomponent

    logical, dimension(field%dim):: applies
    character(len=FIELD_NAME_LEN):: bctype
    type(vector_field), pointer:: surface_field
    integer, dimension(:), pointer:: surface_node_list
    integer :: i,j,k

    bcloop: do i=1, get_boundary_condition_count(field)
       call get_boundary_condition(field, i, type=bctype, &
          surface_node_list=surface_node_list, applies=applies)

       if (bctype/="dirichlet") cycle bcloop
       
       surface_field => extract_surface_field(field, i, "value")
       
       do j=1,size(surface_node_list)
          do k = 1, field%dim
            if(applies(k)) then

              lhscomponent = extract_scalar_field_from_vector_field(lhs, k)

              if(present(rhs).and.present(dt)) then
                rhscomponent = extract_scalar_field_from_vector_field(rhs, k)
                bccomponent = extract_scalar_field_from_vector_field(surface_field, k)

                call addto(rhscomponent, &
                          surface_node_list(j), &
                          ((node_val(bccomponent,j)&
                            -node_val(field,k,surface_node_list(j)) &
                            ) /dt)*INFINITY)

              end if


              call addto(lhscomponent, surface_node_list(j),&
                         INFINITY)
            end if
         end do
       end do

    end do bcloop

  end subroutine apply_dirichlet_conditions_vector_lumped

  subroutine apply_dirichlet_conditions_vector_component_lumped(lhs, rhs, field, dt, dim)
    !!< Apply dirichlet boundary conditions from field to the problem
    !!< defined by lhs and rhs.  lhs is the diagonal of the normal matrix
    !!< so this is only useful for fully explicit problems.
    !!<
    !!< This assumes that boundary conditions are applied in rate of change
    !!< form.
    type(scalar_field), intent(inout) :: lhs
    type(vector_field), intent(inout), optional :: rhs
    type(vector_field), intent(in) :: field
    real, intent(in), optional :: dt
    integer, intent(in) :: dim

    type(scalar_field) :: rhscomponent, bccomponent

    logical, dimension(field%dim):: applies
    character(len=FIELD_NAME_LEN):: bctype
    type(vector_field), pointer:: surface_field
    integer, dimension(:), pointer:: surface_node_list
    integer :: i,j

    bcloop: do i=1, get_boundary_condition_count(field)
       call get_boundary_condition(field, i, type=bctype, &
          surface_node_list=surface_node_list, applies=applies)

       if (bctype/="dirichlet") cycle bcloop
       
       surface_field => extract_surface_field(field, i, "value")
       
       do j=1,size(surface_node_list)
         
          if(applies(dim)) then

            if(present(rhs).and.present(dt)) then
              bccomponent = extract_scalar_field_from_vector_field(surface_field, dim)
              rhscomponent = extract_scalar_field_from_vector_field(rhs, dim)

              call addto(rhscomponent, &
                        surface_node_list(j), &
                        ((node_val(bccomponent,j) &
                          -node_val(field, dim, surface_node_list(j)) &
                          ) /dt)*INFINITY)

            end if

            call addto(lhs, surface_node_list(j),&
                        INFINITY)
          end if
       end do

    end do bcloop

  end subroutine apply_dirichlet_conditions_vector_component_lumped
  
  subroutine derive_collapsed_bcs(input_states, collapsed_states, bctype)
    !!< For the collapsed state collapsed_states, containing the collapsed
    !!< components of fields in input_states, copy across the component
    !!< boundary conditions.
  
    type(state_type), dimension(:), intent(in) :: input_states
    type(state_type), dimension(size(input_states)), intent(inout) :: collapsed_states
    character(len = *), optional, intent(in) :: bctype
  
    character(len = FIELD_NAME_LEN) :: bcname, lbctype
    character(len = OPTION_PATH_LEN) :: bcoption_path
    integer :: i, j, k, l
    integer, dimension(:), allocatable :: bccount
    integer, dimension(:), pointer :: bcsurface_element_list
    logical, dimension(:), allocatable :: bcapplies
    type(mesh_type), pointer :: bcsurface_mesh
    type(scalar_field) :: bcvalue_comp
    type(scalar_field), target :: v_field_comp_wrap
    type(scalar_field_pointer), dimension(:), allocatable :: v_field_comps
    type(vector_field), pointer :: bcvalue, v_field
  
    ewrite(1, *) "In derive_collapsed_bcs"
  
    do i = 1, size(input_states)
      v_field_loop: do j = 1, vector_field_count(input_states(i))
        v_field => extract_vector_field(input_states(i), j)
        if(trim(v_field%name) == "Coordinate") cycle
        ewrite(2, *) "For vector field " // trim(v_field%name) // ", bc count = ", get_boundary_condition_count(v_field)
        if(get_boundary_condition_count(v_field) == 0) cycle v_field_loop
        
        allocate(v_field_comps(v_field%dim))
        do k = 1, v_field%dim
          v_field_comps(k)%ptr => extract_scalar_field(collapsed_states(i), trim(input_states(i)%vector_names(j)) // "%" // int2str(k))
         
          if(.not. associated(v_field_comps(k)%ptr%bc)) then
            ! This scalar field has not bc field, (probably a borrowed
            ! reference). Wrap it so that it isn't a borrowed reference.
            ! Wrapping is cheaper than a copy here.
            v_field_comp_wrap = wrap_scalar_field(v_field_comps(k)%ptr%mesh, v_field_comps(k)%ptr%val, v_field_comps(k)%ptr%name)
            call insert(collapsed_states(i), v_field_comp_wrap, v_field_comp_wrap%name)
            v_field_comps(k)%ptr => v_field_comp_wrap
            call deallocate(v_field_comp_wrap)
          end if
          
          assert(associated(v_field_comps(k)%ptr%bc))
        end do
        
        allocate(bcapplies(v_field%dim))
        allocate(bccount(v_field%dim))
        bccount = 0
        bc_loop: do k = 1, get_boundary_condition_count(v_field)
     
          call get_boundary_condition(v_field, k, name = bcname, type = lbctype, &
            & surface_element_list = bcsurface_element_list, applies = bcapplies, &
            & surface_mesh = bcsurface_mesh, option_path = bcoption_path)
          if(present(bctype)) then
            if(lbctype /= bctype) cycle bc_loop
          end if
          bcvalue => extract_surface_field(v_field, k, "value")
          ewrite(2, *) "Collapsing vector bc " // trim(bcname) // " for field " // trim(v_field%name)
          
          bc_dim_loop: do l = 1, v_field%dim
            ewrite_minmax(bcvalue%val(l)%ptr)  
            if(.not. bcapplies(l)) cycle bc_dim_loop
            
            call add_boundary_condition_surface_elements(v_field_comps(l)%ptr, bcname, bctype, &
              & bcsurface_element_list, option_path = bcoption_path)
              
            call allocate(bcvalue_comp, bcsurface_mesh, name = "value")
            call set(bcvalue_comp, extract_scalar_field(bcvalue, l))
            bccount(l) = bccount(l) + 1
            call insert_surface_field(v_field_comps(l)%ptr, bccount(l), bcvalue_comp)
            call deallocate(bcvalue_comp)
          end do bc_dim_loop
        
        end do bc_loop
        deallocate(bcapplies)
        deallocate(bccount)
        
        deallocate(v_field_comps)
      end do v_field_loop
    end do
    
    ewrite(1, *) "Exiting derive_collapsed_bcs"
  
  end subroutine derive_collapsed_bcs
  
end module boundary_conditions
