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
module plugin_functions

  use fields
  use global_parameters, only: OPTION_PATH_LEN, PYTHON_FUNC_LEN,&
       is_active_process
  use iso_c_binding
  use dl_fortran

  implicit none

  private

  public :: set_from_plugin_function

  interface set_from_plugin_function
     module procedure set_from_plugin_function_scalar,&
          set_from_plugin_function_vector, set_from_plugin_function_tensor
  end interface set_from_plugin_function


contains

  subroutine dynamic_load_function(handle,funptr,filename,function_name,isfortran)

    implicit none

    character (len=*),  intent( in )      :: filename,function_name
    type(c_ptr), intent(inout) :: handle
    type(c_funptr), intent(inout) :: funptr
    logical, intent( in ) :: isfortran

    character (len=1024), pointer :: fstring
    integer :: slen

    handle=DLOpen(trim(filename)//c_null_char,ior(RTLD_NOW,RTLD_GLOBAL))

    if (.not. c_associated(handle)) then
       call c_f_pointer(DLError(),fstring)
       slen=index(fstring,c_null_char) - 1
       FLAbort( 'error in DLOpen:'//fstring(1:slen) )
    end if
    
    if (isfortran) then
       funptr= DLSym(handle,trim(function_name)//'_'//c_null_char)
    else
       funptr= DLSym(handle,trim(function_name)//c_null_char)
    end if

    if (.not. c_associated(funptr)) then
       call c_f_pointer(DLError(),fstring)
       slen=index(fstring,c_null_char) - 1
       FLAbort( 'error in DLSym:'//fstring(1:slen) )
    end if

  end subroutine dynamic_load_function

  subroutine position_on_mesh(position,mesh_position,mesh)
    type(vector_field), intent(in)  :: position
    type(vector_field), intent(out) :: mesh_position
    type(mesh_type), intent(inout)  :: mesh
 
    integer :: stat
    
    ! Remap position onto the mesh provided, allowing interpolation
    ! Note that this routine allocates memory
    call allocate(mesh_position, position%dim,mesh, "Local Position")
    call remap_field(position, mesh_position, stat=stat)
    if(stat==REMAP_ERR_DISCONTINUOUS_CONTINUOUS) then
       ewrite(-1,*) 'Remapping of the coordinates just threw an error because'
       ewrite(-1,*) 'the input coordinates are discontinuous and you are trying'
       ewrite(-1,*) 'to remap them to a continuous field.'
       FLAbort("Why are your coordinates discontinuous?")
    end if
    ! we've just allowed remapping from a higher order to a lower order continuous field as this should be valid for
    ! coordinates
    ! also allowed to remap from unperiodic to periodic... hopefully the python function used will also be periodic!
  end subroutine position_on_mesh

  subroutine set_from_plugin_function_scalar(field, plugin_name,&
       function_name, position, time, interface_type)
    
    implicit none

    type( scalar_field ), intent( inout ) :: field
    type( vector_field ), intent( in), target     :: position
    real, intent (in )                    :: time
    
    character (len=*),  intent( in )      :: plugin_name,function_name
    character (len=*),  intent( in ), optional :: interface_type

    type( vector_field ), pointer :: fposition

    if (field%mesh==position%mesh) then
       fposition => position
    else
       allocate(fposition)
       call position_on_mesh(position,fposition,field%mesh)
    end if

    if (.not. present(interface_type)) then
       !! Default to Fortran 77 style interface
       call set_from_plugin_function_fortran77_cbind_scalar(field, plugin_name,&
       function_name, fposition, time, isfortran=.true.)
       return
    end if

    select case(interface_type)
    case("Fortran90")
       call set_from_plugin_function_fortran90_scalar(field, plugin_name,&
            function_name, fposition, time)
       return
    case("Fortran77")
       call set_from_plugin_function_fortran77_cbind_scalar(field, plugin_name,&
       function_name, fposition, time, isfortran=.true.)
       return
    case("C")
       call set_from_plugin_function_fortran77_cbind_scalar(field, plugin_name,&
       function_name, fposition, time, isfortran=.false.)
       return
    end select

    if (.not. (field%mesh == position%mesh)) then
       call deallocate(fposition)
       deallocate(fposition)
    end if

  end subroutine set_from_plugin_function_scalar

  subroutine set_from_plugin_function_fortran77_cbind_scalar(field, plugin_name,&
       function_name, position, time, isfortran) 
   
    implicit none

    type( scalar_field ), intent( inout ) :: field
    type( vector_field ), intent( in)     :: position
    real, intent (in )                    :: time
    logical, intent( in ) :: isfortran
    character (len=*),  intent( in )      :: plugin_name,function_name

    type(c_ptr) :: handle=C_NULL_PTR
    type(c_funptr) :: funptr=C_NULL_FUNPTR
    integer(c_int) :: status

    abstract interface
       subroutine val(dim,n,field,X,t)
         use iso_c_binding
         integer :: dim,n
         real, dimension(n) :: field
         real, dimension(dim,n) :: X
         real               :: t
       end subroutine val
    end interface
    procedure(val), pointer ::dll_sub

    call dynamic_load_function(handle,funptr,plugin_name,&
         function_name,isfortran)
    
    ! convert C function pointer to Fortran procedure pointer
    call c_f_procpointer(cptr=funptr,fptr=dll_sub)

    call dll_sub(mesh_dim(field),node_count(field),field%val,position%val,time)

    status=DLClose(handle)
    
  end subroutine set_from_plugin_function_fortran77_cbind_scalar

  subroutine set_from_plugin_function_fortran90_scalar(field, plugin_name,&
       function_name, position, time)
    
    implicit none

    type( scalar_field ), intent( inout ) :: field
    type( vector_field ), intent( in), target     :: position
    real, intent (in )                    :: time
    
    character (len=*),  intent( in )      :: plugin_name,function_name

    type(c_ptr) :: handle=C_NULL_PTR
    type(c_funptr) :: funptr=C_NULL_FUNPTR
    integer(c_int) :: status

    abstract interface
       subroutine val(field,X,t)
         use iso_c_binding
         real, dimension(:) :: field
         real, dimension(:,:) :: X
         real               :: t
       end subroutine val
    end interface
    procedure(val), pointer ::dll_sub

    call dynamic_load_function(handle,funptr,plugin_name,&
         function_name,isfortran=.true.)
    
    ! convert C function pointer to Fortran procedure pointer
    call c_f_procpointer(cptr=funptr,fptr=dll_sub)

    call dll_sub(field%val,position%val,time)

    status=DLClose(handle)
    
  end subroutine set_from_plugin_function_fortran90_scalar

  subroutine set_from_plugin_function_vector(field, plugin_name,&
       function_name, position, time, interface_type)
    
    implicit none

    type( vector_field ), intent( inout ) :: field
    type( vector_field ), intent( in), target     :: position
    real, intent (in )                    :: time
    
    character (len=*),  intent( in )      :: plugin_name,function_name
    character (len=*),  intent( in ), optional :: interface_type

    type( vector_field ), pointer :: fposition

    if (field%mesh==position%mesh) then
       fposition => position
    else
       allocate(fposition)
       call position_on_mesh(position,fposition,field%mesh)
    end if

    if (.not. present(interface_type)) then
       !! Default to Fortran 77 style interface
       call set_from_plugin_function_fortran77_cbind_vector(field, plugin_name,&
       function_name, fposition, time, isfortran=.true.)
       return
    end if

    select case(interface_type)
    case("Fortran90")
       call set_from_plugin_function_fortran90_vector(field, plugin_name,&
            function_name, fposition, time)
       return
    case("Fortran77")
       call set_from_plugin_function_fortran77_cbind_vector(field, plugin_name,&
       function_name, fposition, time, isfortran=.true.)
       return
    case("C")
       call set_from_plugin_function_fortran77_cbind_vector(field, plugin_name,&
       function_name, fposition, time, isfortran=.false.)
       return
    end select

    if (.not. (field%mesh == position%mesh)) then
       call deallocate(fposition)
       deallocate(fposition)
    end if

  end subroutine set_from_plugin_function_vector

  subroutine set_from_plugin_function_fortran77_cbind_vector(field, plugin_name,&
       function_name, position, time, isfortran) 
   
    implicit none

    type( vector_field ), intent( inout ) :: field
    type( vector_field ), intent( in)     :: position
    real, intent (in )                    :: time
    logical, intent( in ) :: isfortran
    character (len=*),  intent( in )      :: plugin_name,function_name

    type(c_ptr) :: handle=C_NULL_PTR
    type(c_funptr) :: funptr=C_NULL_FUNPTR
    integer(c_int) :: status

    abstract interface
       subroutine val(xdim,n,fdim,field,X,t)
         use iso_c_binding
         integer :: xdim,n, fdim
         real, dimension(fdim,n) :: field
         real, dimension(xdim,n) :: X
         real               :: t
       end subroutine val
    end interface
    procedure(val), pointer ::dll_sub

    call dynamic_load_function(handle,funptr,plugin_name,&
         function_name,isfortran)
    
    ! convert C function pointer to Fortran procedure pointer
    call c_f_procpointer(cptr=funptr,fptr=dll_sub)

    call dll_sub(mesh_dim(field),node_count(field),field%dim,field%val,position%val,time)

    status=DLClose(handle)
    
  end subroutine set_from_plugin_function_fortran77_cbind_vector

  subroutine set_from_plugin_function_fortran90_vector(field, plugin_name,&
       function_name, position, time)
    
    implicit none

    type( vector_field ), intent( inout ) :: field
    type( vector_field ), intent( in), target     :: position
    real, intent (in )                    :: time
    
    character (len=*),  intent( in )      :: plugin_name,function_name

    type(c_ptr) :: handle=C_NULL_PTR
    type(c_funptr) :: funptr=C_NULL_FUNPTR
    integer(c_int) :: status

    abstract interface
       subroutine val(field,X,t)
         use iso_c_binding
         real, dimension(:,:) :: field
         real, dimension(:,:) :: X
         real               :: t
       end subroutine val
    end interface
    procedure(val), pointer ::dll_sub

    call dynamic_load_function(handle,funptr,plugin_name,&
         function_name,isfortran=.true.)
    
    ! convert C function pointer to Fortran procedure pointer
    call c_f_procpointer(cptr=funptr,fptr=dll_sub)

    call dll_sub(field%val,position%val,time)

    status=DLClose(handle)
    
  end subroutine set_from_plugin_function_fortran90_vector

  subroutine set_from_plugin_function_tensor(field, plugin_name,&
       function_name, position, time, interface_type)
    
    implicit none

    type( tensor_field ), intent( inout ) :: field
    type( vector_field ), intent( in), target     :: position
    real, intent (in )                    :: time
    
    character (len=*),  intent( in )      :: plugin_name,function_name
    character (len=*),  intent( in ), optional :: interface_type

    type( vector_field ), pointer :: fposition

    if (field%mesh==position%mesh) then
       fposition => position
    else
       allocate(fposition)
       call position_on_mesh(position,fposition,field%mesh)
    end if

    if (.not. present(interface_type)) then
       !! Default to Fortran 77 style interface
       call set_from_plugin_function_fortran77_cbind_tensor(field, plugin_name,&
       function_name, fposition, time, isfortran=.true.)
       return
    end if

    select case(interface_type)
    case("Fortran90")
       call set_from_plugin_function_fortran90_tensor(field, plugin_name,&
            function_name, fposition, time)
       return
    case("Fortran77")
       call set_from_plugin_function_fortran77_cbind_tensor(field, plugin_name,&
       function_name, fposition, time, isfortran=.true.)
       return
    case("C")
       call set_from_plugin_function_fortran77_cbind_tensor(field, plugin_name,&
       function_name, fposition, time, isfortran=.false.)
       return
    end select

    if (.not. (field%mesh == position%mesh)) then
       call deallocate(fposition)
       deallocate(fposition)
    end if

  end subroutine set_from_plugin_function_tensor

  subroutine set_from_plugin_function_fortran77_cbind_tensor(field, plugin_name,&
       function_name, position, time, isfortran) 
   
    implicit none

    type( tensor_field ), intent( inout ) :: field
    type( vector_field ), intent( in)     :: position
    real, intent (in )                    :: time
    logical, intent( in ) :: isfortran
    character (len=*),  intent( in )      :: plugin_name,function_name

    type(c_ptr) :: handle=C_NULL_PTR
    type(c_funptr) :: funptr=C_NULL_FUNPTR
    integer(c_int) :: status

    abstract interface
       subroutine val(xdim,n,fdim,field,X,t)
         use iso_c_binding
         integer :: xdim,n
         integer, dimension(2) ::fdim
         real, dimension(fdim(1),fdim(2),n) :: field
         real, dimension(xdim,n) :: X
         real               :: t
       end subroutine val
    end interface
    procedure(val), pointer ::dll_sub

    call dynamic_load_function(handle,funptr,plugin_name,&
         function_name,isfortran)
    
    ! convert C function pointer to Fortran procedure pointer
    call c_f_procpointer(cptr=funptr,fptr=dll_sub)

    call dll_sub(mesh_dim(field),node_count(field),field%dim,field%val,&
         position%val,time)

    status=DLClose(handle)
    
  end subroutine set_from_plugin_function_fortran77_cbind_tensor

  subroutine set_from_plugin_function_fortran90_tensor(field, plugin_name,&
       function_name, position, time)
    
    implicit none

    type( tensor_field ), intent( inout ) :: field
    type( vector_field ), intent( in), target     :: position
    real, intent (in )                    :: time
    
    character (len=*),  intent( in )      :: plugin_name,function_name

    type(c_ptr) :: handle=C_NULL_PTR
    type(c_funptr) :: funptr=C_NULL_FUNPTR
    integer(c_int) :: status

    abstract interface
       subroutine val(field,X,t)
         use iso_c_binding
         real, dimension(:,:,:) :: field
         real, dimension(:,:) :: X
         real               :: t
       end subroutine val
    end interface
    procedure(val), pointer ::dll_sub

    call dynamic_load_function(handle,funptr,plugin_name,&
         function_name,isfortran=.true.)
    
    ! convert C function pointer to Fortran procedure pointer
    call c_f_procpointer(cptr=funptr,fptr=dll_sub)

    call dll_sub(field%val,position%val,time)

    status=DLClose(handle)
    
  end subroutine set_from_plugin_function_fortran90_tensor


end module plugin_functions
