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

! The element types are the same as the VTK element types. Support
! types are VTK_TETRA, VTK_HEXAHEDRON, VTK_TRIANGLE and VTK_QUAD
module picker

  use global_parameters, only : real_4, real_8

  implicit none

  real, parameter, public :: default_ownership_tolerance = 1.0e-3

  interface picker_create
    module procedure picker_create_sp
    
    subroutine picker_create(x, y, z, enlist, nelm, element_type, dim, fortran_zero_index, flat_earth, picker_id)
      use global_parameters, only : real_8
      implicit none
      real(kind = real_8), intent(in) :: x(*)
      real(kind = real_8), intent(in) :: y(*)
      real(kind = real_8), intent(in) :: z(*)
      integer, intent(in) :: enlist(*)
      integer, intent(in) :: nelm
      integer, intent(in) :: element_type
      integer, intent(in) :: dim
      integer, intent(in) :: fortran_zero_index
      integer, intent(in) :: flat_earth
      integer, intent(out)::picker_id
    end subroutine picker_create
  end interface picker_create

  interface picker_destroy
     subroutine picker_destroy(picker_id)
       implicit none
       integer, intent(in) :: picker_id
     end subroutine picker_destroy
  end interface picker_destroy

  interface picker_find2d
    module procedure picker_find2d_sp
  
    subroutine picker_find2d(picker_id, x, y, eid, eshape, tol)
      use global_parameters, only : real_8
      implicit none
      integer, intent(in) :: picker_id
      real(kind = real_8), intent(in) :: x
      real(kind = real_8), intent(in) :: y
      integer, intent(out) :: eid
      real(kind = real_8), intent(out) :: eshape(*)
      real(kind = real_8), intent(in) :: tol
    end subroutine picker_find2d
  end interface picker_find2d
  
  interface picker_find3d
    module procedure picker_find3d_sp
  
    subroutine picker_find3d(picker_id, x, y, z, eid, eshape, tol)
      use global_parameters, only : real_8
      implicit none
      integer, intent(in) :: picker_id
      real(kind = real_8), intent(in) :: x
      real(kind = real_8), intent(in) :: y
      real(kind = real_8), intent(in) :: z
      integer, intent(out) :: eid
      real(kind = real_8), intent(out) :: eshape(*)
      real(kind = real_8), intent(in) :: tol
    end subroutine picker_find3d
  end interface picker_find3d
    
  interface picker_pfind2d
    module procedure picker_pfind2d_sp
  
    subroutine picker_pfind2d(picker_id, x, y, nnodes, eid, eshape, rank, tol, different_domains)
      use global_parameters, only : real_8
      implicit none
      integer, intent(in) :: picker_id
      real(kind = real_8), intent(in) :: x(*)
      real(kind = real_8), intent(in) :: y(*)
      integer, intent(in) :: nnodes
      integer, intent(out) :: eid(*)
      real(kind = real_8), intent(out) :: eshape(*)
      integer, intent(out) :: rank(*)
      real(kind = real_8), intent(in) :: tol
      integer, intent(in) :: different_domains
    end subroutine picker_pfind2d
  end interface picker_pfind2d

  interface picker_pfind3d
    module procedure picker_pfind3d_sp
  
    subroutine picker_pfind3d(picker_id, x, y, z, nnodes, eid, eshape, rank, tol, different_domains)
      use global_parameters, only : real_8
      implicit none
      integer, intent(in) :: picker_id
      real(kind = real_8), intent(in) :: x(*)
      real(kind = real_8), intent(in) :: y(*)
      real(kind = real_8), intent(in) :: z(*)
      integer, intent(in) :: nnodes
      integer, intent(out) :: eid(*)
      real(kind = real_8), intent(out) :: eshape(*)
      integer, intent(out) :: rank(*)
      real(kind = real_8), intent(in) :: tol
      integer, intent(in) :: different_domains
    end subroutine picker_pfind3d
  end interface picker_pfind3d

  private
  
  public :: picker_create, picker_destroy, picker_find2d, picker_find3d, &
    & picker_pfind2d, picker_pfind3d

contains

  subroutine picker_create_sp(x, y, z, enlist, nelm, element_type, dim, fortran_zero_index, flat_earth, picker_id)
    real(kind = real_4), dimension(:), intent(in) :: x
    real(kind = real_4), dimension(size(x)), intent(in) :: y
    real(kind = real_4), dimension(size(x)), intent(in) :: z
    integer, dimension(:), intent(in) :: enlist
    integer, intent(in) :: nelm
    integer, intent(in) :: element_type
    integer, intent(in) :: dim
    integer, intent(in) :: fortran_zero_index
    integer, intent(in) :: flat_earth
    integer, intent(out) :: picker_id
    
    call picker_create(real(x, kind = real_8), real(y, kind = real_8), real(z, kind = real_8), enlist, nelm, element_type, dim, fortran_zero_index, flat_earth, picker_id)
    
  end subroutine picker_create_sp
  
  subroutine picker_find2d_sp(picker_id, x, y, eid, eshape, tol)
    integer, intent(in) :: picker_id
    real(kind = real_4), intent(in) :: x
    real(kind = real_4), intent(in) :: y
    integer, intent(out) :: eid
    real(kind = real_4), dimension(:), intent(out) :: eshape
    real(kind = real_4), intent(in) :: tol
    
    real(kind = real_8), dimension(size(eshape)) :: leshape
    
    call picker_find2d(picker_id, real(x, kind = real_8), real(y, kind = real_8), eid, leshape, real(tol, kind = real_8))
    eshape = leshape
    
  end subroutine picker_find2d_sp
  
  subroutine picker_find3d_sp(picker_id, x, y, z, eid, eshape, tol)
    integer, intent(in) :: picker_id
    real(kind = real_4), intent(in) :: x
    real(kind = real_4), intent(in) :: y
    real(kind = real_4), intent(in) :: z
    integer, intent(out) :: eid
    real(kind = real_4), dimension(:), intent(out) :: eshape
    real(kind = real_4), intent(in) :: tol
    
    real(kind = real_8), dimension(size(eshape)) :: leshape
    
    call picker_find3d(picker_id, real(x, kind = real_8), real(y, kind = real_8), real(z, kind = real_8), eid, leshape, real(tol, kind = real_8))
    eshape = leshape
    
  end subroutine picker_find3d_sp
  
  subroutine picker_pfind2d_sp(picker_id, x, y, nnodes, eid, eshape, rank, tol, different_domains)
    integer, intent(in) :: picker_id
    real(kind = real_4), dimension(:), intent(in) :: x
    real(kind = real_4), dimension(:), intent(in) :: y
    integer, intent(in) :: nnodes
    integer, dimension(:), intent(out) :: eid
    real(kind = real_4), dimension(:), intent(out) :: eshape
    integer, dimension(:), intent(out) :: rank
    real(kind = real_4), intent(in) :: tol
    integer, intent(in) :: different_domains
    
    real(kind = real_8), dimension(size(eshape)) :: leshape
    
    call picker_pfind2d(picker_id, real(x, kind = real_8), real(y, kind = real_8), nnodes, eid, leshape, rank, real(tol, kind = real_8), different_domains)
    eshape = leshape
    
  end subroutine picker_pfind2d_sp
  
  subroutine picker_pfind3d_sp(picker_id, x, y, z, nnodes, eid, eshape, rank, tol, different_domains)
    integer, intent(in) :: picker_id
    real(kind = real_4), dimension(:), intent(in) :: x
    real(kind = real_4), dimension(:), intent(in) :: y
    real(kind = real_4), dimension(:), intent(in) :: z
    integer, intent(in) :: nnodes
    integer, dimension(:), intent(out) :: eid
    real(kind = real_4), dimension(:), intent(out) :: eshape
    integer, dimension(:), intent(out) :: rank
    real(kind = real_4), intent(in) :: tol
    integer, intent(in) :: different_domains
    
    real(kind = real_8), dimension(size(eshape)) :: leshape
    
    call picker_pfind3d(picker_id, real(x, kind = real_8), real(y, kind = real_8), real(z, kind = real_8), nnodes, eid, leshape, rank, real(tol, kind = real_8), different_domains)
    eshape = leshape
    
  end subroutine picker_pfind3d_sp

end module picker
