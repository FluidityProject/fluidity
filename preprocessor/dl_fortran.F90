
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


!     wrapper for c dlopen functionality.
!     based on examples at http://cims.nyu.edu/~donev/Fortran/DLL/DLL.Forum.txt

module dl_fortran
   use, intrinsic::  iso_c_binding
!   use iso_c_utilities
   implicit none
   private

   public :: DLOpen, DLSym, DLClose, DLError ! DL API
   
   ! Valid modes for mode in DLOpen:
   integer, parameter, public :: RTLD_LAZY=1, RTLD_NOW=2, RTLD_GLOBAL=256, RTLD_LOCAL=0
      ! Obtained from the output of the previously listed C program 
         
   interface ! all we need is interfaces for the prototypes in <dlfcn.h>
      function DLOpen(file,mode) RESULT(handle) BIND(C,NAME='dlopen')
         ! void *dlopen(const char *file, int mode);
        use, intrinsic::  iso_c_binding
         character(c_char), dimension(*), intent(in) :: file
            ! C strings should be declared as character arrays
         integer(c_int), value :: mode
         type(c_ptr) :: handle
      end function
      function DLSym(handle,name) result(funptr) bind(c,name='dlsym')
         ! void *dlsym(void *handle, const char *name);
        use, intrinsic::  iso_c_binding
         type(c_ptr), value :: handle
         character(c_char), dimension(*), intent(in) :: name
         type(c_funptr) :: funptr ! a function pointer
      end function
      function DLClose(handle) result(status) bind(c,name='dlclose')
         ! int dlclose(void *handle);
        use, intrinsic::  iso_c_binding
         type(c_ptr), value :: handle
         integer(c_int) :: status
      end function
      function dlerror() result(error) bind(c,name='dlerror')
         ! char *dlerror(void);
        use, intrinsic::  iso_c_binding
         type(c_ptr) :: error
      end function         
   end interface
      
end module
