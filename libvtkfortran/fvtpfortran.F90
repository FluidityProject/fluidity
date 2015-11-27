!!$ Copyright (C) 2004- Imperial College London and others.
!!$   
!!$   Please see the AUTHORS file in the main source directory for a full
!!$   list of copyright holders.
!!$   
!!$   Adrian Umpleby
!!$   Applied Modelling and Computation Group
!!$   Department of Earth Science and Engineering
!!$   Imperial College London
!!$   
!!$   adrian@imperial.ac.uk
!!$   
!!$   This library is free software; you can redistribute it and/or
!!$   modify it under the terms of the GNU Lesser General Public
!!$   License as published by the Free Software Foundation; either
!!$   version 2.1 of the License.
!!$   
!!$   This library is distributed in the hope that it will be useful,
!!$   but WITHOUT ANY WARRANTY; without even the implied warranty of
!!$   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!!$   Lesser General Public License for more details.
!!$   
!!$   You should have received a copy of the GNU Lesser General Public
!!$   License along with this library; if not, write to the Free Software
!!$   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!!$   USA

module vtpfortran
  !!< This module merely contains explicit interfaces to allow the
  !!< convenient use of vtkfortran in fortran.
  use iso_c_binding
  
  private

  ! Element types from VTK
  integer, public, parameter :: VTK_VERTEX=1
  integer, public, parameter :: VTK_POLY_VERTEX=2
  integer, public, parameter :: VTK_LINE=3
  integer, public, parameter :: VTK_POLY_LINE=4
  integer, public, parameter :: VTK_TRIANGLE=5
  integer, public, parameter :: VTK_TRIANGLE_STRIP=6
  integer, public, parameter :: VTK_POLYGON=7
  integer, public, parameter :: VTK_PIXEL=8
  integer, public, parameter :: VTK_QUAD=9
  integer, public, parameter :: VTK_TETRA=10
  integer, public, parameter :: VTK_VOXEL=11
  integer, public, parameter :: VTK_HEXAHEDRON=12
  integer, public, parameter :: VTK_WEDGE=13
  integer, public, parameter :: VTK_PYRAMID=14
  
  integer, public, parameter :: VTK_QUADRATIC_EDGE=21
  integer, public, parameter :: VTK_QUADRATIC_TRIANGLE=22
  integer, public, parameter :: VTK_QUADRATIC_QUAD=23
  integer, public, parameter :: VTK_QUADRATIC_TETRA=24
  integer, public, parameter :: VTK_QUADRATIC_HEXAHEDRON=25

  public :: vtpopen, vtpclose, vtppclose, vtpwritepoints, vtpwritesn,&
       & vtpwritesc, vtpwritevn, vtpwritevc, vtpwritetn, vtpwritetc, &
       & vtpsetactivescalars, vtpsetactivevectors, &
       & vtpsetactivetensors, vtpgetsizes, vtpread, vtpopentoread

  interface vtpopen
     subroutine vtpopen_c(outName, len1, vtkTitle, len2) bind(c,name="vtpopen")
       use iso_c_binding
       implicit none
       character(kind=c_char,len=1), dimension(*) :: outName
       integer(kind=c_int) :: len1
       character(kind=c_char,len=1), dimension(*) :: vtkTitle
       integer(kind=c_int) :: len2
     end subroutine vtpopen_c
     module procedure vtpopen_f90
  end interface

  interface vtpopentoread
     subroutine vtpopentoread_c(outName, len1, vtkTitle, len2) bind(c,name="vtpopentoread")
       use iso_c_binding
       implicit none
       character(kind=c_char,len=1), dimension(*) :: outName
       integer(kind=c_int) :: len1
       character(kind=c_char,len=1), dimension(*) :: vtkTitle
       integer(kind=c_int) :: len2
     end subroutine vtpopentoread_c
     module procedure vtpopentoread_f90
  end interface
  
  interface vtpclose
     ! Close the current vtp file.
     subroutine vtpclose() bind(c)
     end subroutine vtpclose
  end interface

  interface vtppclose
     ! Close the current vtp file - creates a parallel file.
     subroutine vtppclose(rank, npartitions) bind(c)
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: rank, npartitions
     end subroutine vtppclose
  end interface

  interface vtpwritepoints
     ! Write mesh information to the current vtp file.
     SUBROUTINE VTPWRITEPOINTS(NPoints, x, y, z) bind(c)
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: NPoints
       REAL(c_float) :: x(*)
       REAL(c_float) :: y(*)
       REAL(c_float) :: z(*)
     end SUBROUTINE VTPWRITEPOINTS
     SUBROUTINE VTPWRITEPOINTSD(NPoints, x, y, z)  bind(c)
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: NPoints
       REAL(c_double) :: x(*)
       REAL(c_double) :: y(*)
       REAL(c_double) :: z(*)
     end SUBROUTINE VTPWRITEPOINTSD
  end interface

  interface vtpwritesn
     ! Write a scalar field to the current vtp file.
     SUBROUTINE VTPWRITEISN_C(vect, name, len) bind(c,name="vtpwriteisn")
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: vect(*)
       character(kind=c_char,len=1), dimension(*) :: name
       integer(kind=c_int) :: len
     end SUBROUTINE VTPWRITEISN_C
     SUBROUTINE VTPWRITEFSN_C(vect, name, len) bind(c,name="vtpwritefsn")
       use iso_c_binding
       implicit none
       REAL(c_float) :: vect(*)
       character(kind=c_char,len=1), dimension(*) :: name
       integer(kind=c_int) :: len
     end SUBROUTINE VTPWRITEFSN_C
     SUBROUTINE VTPWRITEDSN_C(vect, name, len) bind(c,name="vtpwritedsn")
       use iso_c_binding
       implicit none
       REAL(c_double) :: vect(*)
       character(kind=c_char,len=1), dimension(*) :: name
       integer(kind=c_int) :: len
     end SUBROUTINE VTPWRITEDSN_C
     module procedure vtpwriteisn_f90, vtpwritefsn_f90, vtpwritedsn_f90
  end interface

  interface vtpwritesc
     ! Write a scalar field (cell-based) to the current vtp file.
     SUBROUTINE VTPWRITEISC_C(vect, name, len) bind(c,name="vtpwriteisc")
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: vect(*)
       character(kind=c_char,len=1), dimension(*) :: name
       integer(kind=c_int) :: len
     end SUBROUTINE VTPWRITEISC_C
     SUBROUTINE VTPWRITEFSC_C(vect, name, len) bind(c,name="vtpwritefsc")
       use iso_c_binding
       implicit none
       REAL(c_float) :: vect(*)
       character(kind=c_char,len=1), dimension(*) :: name
       integer(kind=c_int) :: len
     end SUBROUTINE VTPWRITEFSC_C
     SUBROUTINE VTPWRITEDSC_C(vect, name, len) bind(c,name="vtpwritedsc")
       use iso_c_binding
       implicit none
       REAL(c_double) :: vect(*)
       character(kind=c_char,len=1), dimension(*) :: name
       integer(kind=c_int) :: len
     end SUBROUTINE VTPWRITEDSC_C
     module procedure vtpwriteisc_f90, vtpwritefsc_f90, vtpwritedsc_f90
  end interface

  interface vtpwritevn
     ! Write a vector field to the current vtp file.
     SUBROUTINE VTPWRITEFVN_C(vx, vy, vz, name, len) bind(c,name="vtpwritefvn")
       use iso_c_binding
       implicit none
       REAL(c_float) :: vx(*)
       REAL(c_float) :: vy(*)
       REAL(c_float) :: vz(*)
       character(kind=c_char,len=1), dimension(*) :: name
       integer(kind=c_int) :: len
     end SUBROUTINE VTPWRITEFVN_C
     SUBROUTINE VTPWRITEDVN_C(vx, vy, vz, name, len) bind(c,name="vtpwritedvn")
       use iso_c_binding
       implicit none
       REAL(c_double) :: vx(*)
       REAL(c_double) :: vy(*)
       REAL(c_double) :: vz(*)
       character(kind=c_char,len=1), dimension(*) :: name
       integer(kind=c_int) :: len
     end SUBROUTINE VTPWRITEDVN_C
     module procedure vtpwritefvn_f90, vtpwritedvn_f90
  end interface

  interface vtpwritevc
     ! Write a vector field (cell-based) to the current vtp file.
     SUBROUTINE VTPWRITEFVC_C(vx, vy, vz, name, len) bind(c,name="vtpwritefvc")
       use iso_c_binding
       implicit none
       REAL(c_float) :: vx(*)
       REAL(c_float) :: vy(*)
       REAL(c_float) :: vz(*)
       character(kind=c_char,len=1), dimension(*) :: name
       integer(kind=c_int) :: len
     end SUBROUTINE VTPWRITEFVC_C
     SUBROUTINE VTPWRITEDVC_C(vx, vy, vz, name, len) bind(c,name="vtpwritedvc")
       use iso_c_binding
       implicit none
       REAL(c_double) :: vx(*)
       REAL(c_double) :: vy(*)
       REAL(c_double) :: vz(*)
       character(kind=c_char,len=1), dimension(*) :: name
       integer(kind=c_int) :: len
     end SUBROUTINE VTPWRITEDVC_C
     module procedure vtpwritefvc_f90, vtpwritedvc_f90
  end interface

  interface vtpwritetn
     ! Write a tensor field to the current vtp file.
     SUBROUTINE VTPWRITEFTN_C(v1, v2, v3, v4, v5, v6, v7, v8, v9, name,&
          & len) bind(c,name="vtpwriteftn")
       use iso_c_binding
       implicit none
       REAL(c_float) :: v1(*)
       REAL(c_float) :: v2(*)
       REAL(c_float) :: v3(*)
       REAL(c_float) :: v4(*)
       REAL(c_float) :: v5(*)
       REAL(c_float) :: v6(*)
       REAL(c_float) :: v7(*)
       REAL(c_float) :: v8(*)
       REAL(c_float) :: v9(*)
       character(kind=c_char,len=1), dimension(*) :: name
       integer(kind=c_int) :: len
     end SUBROUTINE VTPWRITEFTN_C
     SUBROUTINE VTPWRITEDTN_C(v1, v2, v3, v4, v5, v6, v7, v8, v9, name, len)&
          & bind(c,name="vtpwritedtn") 
       use iso_c_binding
       implicit none
       REAL(c_double) :: v1(*)
       REAL(c_double) :: v2(*)
       REAL(c_double) :: v3(*)
       REAL(c_double) :: v4(*)
       REAL(c_double) :: v5(*)
       REAL(c_double) :: v6(*)
       REAL(c_double) :: v7(*)
       REAL(c_double) :: v8(*)
       REAL(c_double) :: v9(*)
       character(kind=c_char,len=1), dimension(*) :: name
       integer(kind=c_int) :: len
     end SUBROUTINE VTPWRITEDTN_C
     module procedure vtpwriteftn_f90, vtpwritedtn_f90
  end interface

  interface vtpwritetc
     ! Write a tensor field (cell-based) to the current vtp file.
     SUBROUTINE VTPWRITEFTC_C(v1, v2, v3, v4, v5, v6, v7, v8, v9, name,&
          & len) bind(c,name="vtpwriteftc")
       use iso_c_binding
       implicit none
       REAL(c_float) :: v1(*)
       REAL(c_float) :: v2(*)
       REAL(c_float) :: v3(*)
       REAL(c_float) :: v4(*)
       REAL(c_float) :: v5(*)
       REAL(c_float) :: v6(*)
       REAL(c_float) :: v7(*)
       REAL(c_float) :: v8(*)
       REAL(c_float) :: v9(*)
       character(kind=c_char,len=1), dimension(*) :: name
       integer(kind=c_int) :: len
     end SUBROUTINE VTPWRITEFTC_C
     SUBROUTINE VTPWRITEDTC_C(v1, v2, v3, v4, v5, v6, v7, v8, v9, name,&
          & len) bind(c,name="vtpwritedtc") 
       use iso_c_binding
       implicit none
       REAL(c_double) :: v1(*)
       REAL(c_double) :: v2(*)
       REAL(c_double) :: v3(*)
       REAL(c_double) :: v4(*)
       REAL(c_double) :: v5(*)
       REAL(c_double) :: v6(*)
       REAL(c_double) :: v7(*)
       REAL(c_double) :: v8(*)
       REAL(c_double) :: v9(*)
       character(kind=c_char,len=1), dimension(*) :: name
       integer(kind=c_int) :: len
     end SUBROUTINE VTPWRITEDTC_C
     module procedure vtpwriteftc_f90, vtpwritedtc_f90
  end interface
 
  interface vtpsetactivescalars
    subroutine vtpsetactivescalars_c(name, len) bind(c,name="vtpsetactivescalars")
       use iso_c_binding
       implicit none
       character(kind=c_char,len=1), dimension(*) :: name
       integer(kind=c_int) :: len
     end subroutine vtpsetactivescalars_c
     module procedure vtpsetactivescalars_f90
  end interface

  interface vtpsetactivevectors
     subroutine vtpsetactivevectors_c(name, len) bind(c,name="vtpsetactivevectors")
       use iso_c_binding
       implicit none
       character(kind=c_char,len=1), dimension(*) :: name
       integer(kind=c_int) :: len
     end subroutine vtpsetactivevectors_c
     module procedure vtpsetactivevectors_f90
  end interface

  interface vtpsetactivetensors
     subroutine vtpsetactivetensors_c(name, len) bind(c,name="vtpsetactivetensors")
       use iso_c_binding
       implicit none
       character(kind=c_char,len=1), dimension(*) :: name
       integer(kind=c_int) :: len
     end subroutine vtpsetactivetensors_c
     module procedure vtpsetactivetensors_f90
  end interface

  interface vtpgetsizes
     subroutine vtpgetsizes(NPoints, NScalars, NVectors, NTensors) bind(c)
       use iso_c_binding
       implicit none
       integer(kind=c_int) :: NPoints, NScalars, NVectors, NTensors
     end subroutine vtpgetsizes
  end interface

  interface vtpread
     subroutine vtpread(x, y, z, scalars, vectors, tensors) bind(c)
       use iso_c_binding
       implicit none
       real(c_double) :: x(*), y(*), z(*), scalars(*), vectors(*), tensors(*)
     end subroutine vtpread
  end interface vtpread


contains

  subroutine vtpopen_f90(outName, vtkTitle)
    ! Wrapper routine with nicer interface.
    character(len=*), intent(in) :: outName, vtkTitle

    call vtpopen_c(outName, len(outName), vtkTitle, len(vtkTitle))

  end subroutine vtpopen_f90

  subroutine vtpopentoread_f90(outName, vtkTitle)
    ! Wrapper routine with nicer interface.
    character(len=*), intent(in) :: outName, vtkTitle

    call vtpopentoread_c(outName, len(outName), vtkTitle, len(vtkTitle))

  end subroutine vtpopentoread_f90
  
  subroutine vtpwriteisn_f90(vect, name)
    ! Wrapper routine with nicer interface.
    integer, intent(in) :: vect(*)
    character(len=*), intent(in) :: name

    call vtpwriteisn_c(vect, name, len(name))

  end subroutine vtpwriteisn_f90

  subroutine vtpwritefsn_f90(vect, name)
    ! Wrapper routine with nicer interface.
    real(c_float), intent(in) :: vect(*)
    character(len=*), intent(in) :: name

    call vtpwritefsn_c(vect, name, len(name))

  end subroutine vtpwritefsn_f90

  subroutine vtpwritedsn_f90(vect, name)
    ! Wrapper routine with nicer interface.
    real(c_double), intent(in) :: vect(*)
    character(len=*), intent(in) :: name

    call vtpwritedsn_c(vect, name, len(name))
  end subroutine vtpwritedsn_f90

  subroutine vtpwriteisc_f90(vect, name)
    ! Wrapper routine with nicer interface.
    integer, intent(in) :: vect(*)
    character(len=*), intent(in) :: name

    call vtpwriteisc_c(vect, name, len(name))

  end subroutine vtpwriteisc_f90

  subroutine vtpwritefsc_f90(vect, name)
    ! Wrapper routine with nicer interface.
    real(c_float), intent(in) :: vect(*)
    character(len=*), intent(in) :: name

    call vtpwritefsc_c(vect, name, len(name))

  end subroutine vtpwritefsc_f90

  subroutine vtpwritedsc_f90(vect, name)
    ! Wrapper routine with nicer interface.
    real(c_double), intent(in) :: vect(*)
    character(len=*), intent(in) :: name

    call vtpwritedsc_c(vect, name, len(name))
  end subroutine vtpwritedsc_f90

  subroutine vtpwritefvn_f90(vx, vy, vz, name)
    REAL(c_float), intent(in) :: vx(*), vy(*), vz(*)
    character(len=*) name

    call vtpwritefvn_c(vx, vy, vz, name, len(name))

  end subroutine vtpwritefvn_f90

  subroutine vtpwritedvn_f90(vx, vy, vz, name)
    REAL(c_double), intent(in) :: vx(*), vy(*), vz(*)
    character(len=*) name

    call vtpwritedvn_c(vx, vy, vz, name, len(name))

  end subroutine vtpwritedvn_f90  

  subroutine vtpwritefvc_f90(vx, vy, vz, name)
    REAL(c_float), intent(in) :: vx(*), vy(*), vz(*)
    character(len=*) name

    call vtpwritefvc_c(vx, vy, vz, name, len(name))

  end subroutine vtpwritefvc_f90

  subroutine vtpwritedvc_f90(vx, vy, vz, name)
    REAL(c_double), intent(in) :: vx(*), vy(*), vz(*)
    character(len=*) name

    call vtpwritedvc_c(vx, vy, vz, name, len(name))

  end subroutine vtpwritedvc_f90

  subroutine vtpwriteftn_f90(v1, v2, v3, v4, v5, v6, v7, v8, v9, name)
    REAL(c_float), intent(in) :: v1(*), v2(*), v3(*), v4(*), v5(*), v6(*), v7(*), v8(*), v9(*)
    character(len=*) name

    call vtpwriteftn_c(v1, v2, v3, v4, v5, v6, v7, v8, v9, name, len(name))

  end subroutine vtpwriteftn_f90

  subroutine vtpwritedtn_f90(v1, v2, v3, v4, v5, v6, v7, v8, v9, name)
    REAL(c_double), intent(in) :: v1(*), v2(*), v3(*), v4(*), v5(*), v6(*), v7(*), v8(*), v9(*)
    character(len=*) name

    call vtpwritedtn_c(v1, v2, v3, v4, v5, v6, v7, v8, v9, name, len(name))

  end subroutine vtpwritedtn_f90  
  
  subroutine vtpwriteftc_f90(v1, v2, v3, v4, v5, v6, v7, v8, v9, name)
    REAL(c_float), intent(in) :: v1(*), v2(*), v3(*), v4(*), v5(*), v6(*), v7(*), v8(*), v9(*)
    character(len=*) name

    call vtpwriteftc_c(v1, v2, v3, v4, v5, v6, v7, v8, v9, name, len(name))

  end subroutine vtpwriteftc_f90

  subroutine vtpwritedtc_f90(v1, v2, v3, v4, v5, v6, v7, v8, v9, name)
    REAL(c_double), intent(in) :: v1(*), v2(*), v3(*), v4(*), v5(*), v6(*), v7(*), v8(*), v9(*)
    character(len=*) name

    call vtpwritedtc_c(v1, v2, v3, v4, v5, v6, v7, v8, v9, name, len(name))

  end subroutine vtpwritedtc_f90  
  
  subroutine vtpsetactivescalars_f90(name)
    character(len=*) name
    
    call vtpsetactivescalars_c(name, len_trim(name))
    
  end subroutine vtpsetactivescalars_f90

  subroutine vtpsetactivevectors_f90(name)
    character(len=*) name
    
    call vtpsetactivevectors_c(name, len_trim(name))
    
  end subroutine vtpsetactivevectors_f90

  subroutine vtpsetactivetensors_f90(name)
    character(len=*) name
    
    call vtpsetactivetensors_c(name, len_trim(name))
    
  end subroutine vtpsetactivetensors_f90

end module vtpfortran
