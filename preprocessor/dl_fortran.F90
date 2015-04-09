!! file dl_fortran.F90

!! wrapper for c dlopen functionality.
!! based on examples at http://cims.nyu.edu/~donev/Fortran/DLL/DLL.Fortran.txt

MODULE DL_fortran
   USE, intrinsic::  ISO_C_BINDING
!   USE ISO_C_UTILITIES
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: DLOpen, DLSym, DLClose, DLError ! DL API
   
   ! Valid modes for mode in DLOpen:
   INTEGER, PARAMETER, PUBLIC :: RTLD_LAZY=1, RTLD_NOW=2, RTLD_GLOBAL=256, RTLD_LOCAL=0
      ! Obtained from the output of the previously listed C program 
         
   INTERFACE ! All we need is interfaces for the prototypes in <dlfcn.h>
      FUNCTION DLOpen(file,mode) RESULT(handle) BIND(C,NAME='dlopen')
         ! void *dlopen(const char *file, int mode);
        USE, intrinsic::  ISO_C_BINDING
         CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: file
            ! C strings should be declared as character arrays
         INTEGER(C_INT), VALUE :: mode
         TYPE(C_PTR) :: handle
      END FUNCTION
      FUNCTION DLSym(handle,name) RESULT(funptr) BIND(C,NAME='dlsym')
         ! void *dlsym(void *handle, const char *name);
        USE, intrinsic::  ISO_C_BINDING
         TYPE(C_PTR), VALUE :: handle
         CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: name
         TYPE(C_FUNPTR) :: funptr ! A function pointer
      END FUNCTION
      FUNCTION DLClose(handle) RESULT(status) BIND(C,NAME='dlclose')
         ! int dlclose(void *handle);
        USE, intrinsic::  ISO_C_BINDING
         TYPE(C_PTR), VALUE :: handle
         INTEGER(C_INT) :: status
      END FUNCTION
      FUNCTION DLError() RESULT(error) BIND(C,NAME='dlerror')
         ! char *dlerror(void);
        USE, intrinsic::  ISO_C_BINDING
         TYPE(C_PTR) :: error
      END FUNCTION         
   END INTERFACE
      
END MODULE
