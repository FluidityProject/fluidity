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

  implicit none

private

public :: set_from_plugin_function

contains

  subroutine set_from_plugin_function(field, plugin_name,&
       function_name, position, time, interface_type)

    
    use iso_c_binding
    use dl_fortran
    
    implicit none

    type( scalar_field ), intent( inout ) :: field
    type( vector_field ), intent( in), target     :: position
    real, intent (in )                    :: time
    
    character (len=*),  intent( in )      :: plugin_name,function_name

    character (len=*),  intent( in ), optional :: interface_type


    if (.not. present(interface_type)) then
       !! Default to Fortran 77 style interface
       call set_from_plugin_function_fortran77_cbind(field, plugin_name,&
       function_name, position, time, isfortran=.true.)
       return
    end if

    select case(interface_type)
    case("Fortran90")
       call set_from_plugin_function_fortran90(field, plugin_name,&
            function_name, position, time)
       return
    case("Fortran77")
       call set_from_plugin_function_fortran77_cbind(field, plugin_name,&
       function_name, position, time, isfortran=.true.)
       return
    case("C")
       call set_from_plugin_function_fortran77_cbind(field, plugin_name,&
       function_name, position, time, isfortran=.false.)
       return
    end select
  end subroutine set_from_plugin_function

  subroutine set_from_plugin_function_fortran77_cbind(field, plugin_name,&
       function_name, position, time, isfortran)

    use iso_c_binding
    use dl_fortran
    
    implicit none

    type( scalar_field ), intent( inout ) :: field
    type( vector_field ), intent( in), target     :: position
    real, intent (in )                    :: time
    
    character (len=*),  intent( in )      :: plugin_name,function_name

    logical, intent( in ) :: isfortran

character (len=1024), pointer :: fstring
    integer :: slen
    type(c_ptr) :: handle=C_NULL_PTR
    type(c_funptr) :: funptr=C_NULL_FUNPTR
    
    integer(c_int) :: status
    type( vector_field ), target :: lposition
    type( vector_field ), pointer :: fposition

    integer :: stat

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


    
    if (field%mesh==position%mesh) then
       
       fposition => position
    else
       ! Remap position first.
       call allocate(lposition, position%dim, field%mesh, "Local Position")
       call remap_field(position, lposition, stat=stat)
       if(stat==REMAP_ERR_DISCONTINUOUS_CONTINUOUS) then
          ewrite(-1,*) 'Remapping of the coordinates just threw an error because'
          ewrite(-1,*) 'the input coordinates are discontinuous and you are trying'
          ewrite(-1,*) 'to remap them to a continuous field.'
          FLAbort("Why are your coordinates discontinuous?")
       end if
       ! we've just allowed remapping from a higher order to a lower order continuous field as this should be valid for
       ! coordinates
       ! also allowed to remap from unperiodic to periodic... hopefully the python function used will also be periodic!

       fposition => lposition
    end if


    handle=DLOpen(trim(plugin_name)//C_NULL_CHAR,IOR(RTLD_NOW,RTLD_GLOBAL))
    
    if (.not. c_associated(handle)) then
       call c_f_pointer(DLError(),fstring)
       slen=index(fstring,c_null_char) - 1
       FLAbort( 'error in dlopen:'//fstring(1:slen) )
    end if

    
    if (isfortran) then
       funptr= DLSym(handle,trim(function_name)//'_'//C_NULL_CHAR)
    else
       funptr= DLSym(handle,trim(function_name)//C_NULL_CHAR)
    end if

    if (.not. c_associated(funptr)) then
       call c_f_pointer(DLError(),fstring)
       slen=index(fstring,c_null_char) - 1
       FLAbort( 'error in dlsym:'//fstring(1:slen) )
    end if

    
    ! convert C function pointer to Fortran procedure pointer
    call c_f_procpointer(cptr=funptr,fptr=dll_sub)

    call dll_sub(mesh_dim(field),node_count(field),field%val,fposition%val,time)

    status=DLClose(handle)
    
  end subroutine set_from_plugin_function_fortran77_cbind

  subroutine set_from_plugin_function_fortran90(field, plugin_name,&
       function_name, position, time)

    use iso_c_binding
    use dl_fortran
    
    implicit none

    type( scalar_field ), intent( inout ) :: field
    type( vector_field ), intent( in), target     :: position
    real, intent (in )                    :: time
    
    character (len=*),  intent( in )      :: plugin_name,function_name

    character (len=1024), pointer :: fstring
    integer :: slen
    type(c_ptr) :: handle=C_NULL_PTR
    type(c_funptr) :: funptr=C_NULL_FUNPTR
    
    integer(c_int) :: status
    type( vector_field ), target :: lposition
    type( vector_field ), pointer :: fposition

    integer :: stat

    abstract interface
       subroutine val(field,X,t)
         use iso_c_binding
         real, dimension(:) :: field
         real, dimension(:,:) :: X
         real               :: t
       end subroutine val
    end interface
    procedure(val), pointer ::dll_sub


    
    if (field%mesh==position%mesh) then
       
       fposition => position
    else
       ! Remap position first.
       call allocate(lposition, position%dim, field%mesh, "Local Position")
       call remap_field(position, lposition, stat=stat)
       if(stat==REMAP_ERR_DISCONTINUOUS_CONTINUOUS) then
          ewrite(-1,*) 'Remapping of the coordinates just threw an error because'
          ewrite(-1,*) 'the input coordinates are discontinuous and you are trying'
          ewrite(-1,*) 'to remap them to a continuous field.'
          FLAbort("Why are your coordinates discontinuous?")
       end if
       ! we've just allowed remapping from a higher order to a lower order continuous field as this should be valid for
       ! coordinates
       ! also allowed to remap from unperiodic to periodic... hopefully the python function used will also be periodic!

       fposition => lposition
    end if


    handle=DLOpen(trim(plugin_name)//C_NULL_CHAR,IOR(RTLD_NOW,RTLD_GLOBAL))
    
    if (.not. c_associated(handle)) then
       call c_f_pointer(DLError(),fstring)
       slen=index(fstring,c_null_char) - 1
       FLAbort( 'error in dlopen:'//fstring(1:slen) )
    end if

    
    funptr= DLSym(handle,trim(function_name)//'_'//C_NULL_CHAR)

    if (.not. c_associated(funptr)) then
       call c_f_pointer(DLError(),fstring)
       slen=index(fstring,c_null_char) - 1
       FLAbort( 'error in dlsym:'//fstring(1:slen) )
    end if

    
    ! convert C function pointer to Fortran procedure pointer
    call c_f_procpointer(cptr=funptr,fptr=dll_sub)


    call dll_sub(field%val,fposition%val,time)

    status=DLClose(handle)
    
  end subroutine set_from_plugin_function_fortran90


end module plugin_functions
