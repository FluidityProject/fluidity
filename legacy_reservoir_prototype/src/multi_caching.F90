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
#include "fdebug.h"



module multiphase_caching

  use spud


  implicit none

  integer :: cache_level=0

  interface reshape_vector2pointer
      module procedure reshape_vector2pointer_A
      module procedure reshape_vector2pointer_B
      module procedure reshape_vector2pointer_C
      module procedure reshape_vector2pointer_D
  end interface reshape_vector2pointer
contains

  subroutine set_caching_level()

    integer :: i
    cache_level=0
    if (have_option('/caching/no_cache_shape_functions')) then
       return
    else
       do i=1,bit_size(cache_level)
          cache_level=ibset(cache_level,i)
       end do
    end if

  end subroutine set_caching_level


#ifdef USING_GFORTRAN
    !These subroutines are used to avoid problems with the reshaping under intel compilers.
    subroutine reshape_vector2pointer_A(vector, pointr, dim1,dim2)
    implicit none
    real, dimension(:,:), pointer, intent(inout) :: pointr
    real, dimension(:), intent(in), target :: vector
    integer, intent(in) :: dim1, dim2
        pointr(1:dim1,1:dim2) => vector
    end subroutine reshape_vector2pointer_A

    subroutine reshape_vector2pointer_B(vector, pointr,dim1,dim2,dim3)
    implicit none
    real, dimension(:,:,:), pointer, intent(inout) :: pointr
    real, dimension(:), intent(in), target :: vector
    integer, intent(in) :: dim1, dim2, dim3
        pointr(1:dim1,1:dim2,1:dim3) => vector
    end subroutine reshape_vector2pointer_B

    subroutine reshape_vector2pointer_C(vector, pointr, dim1,dim2)
    implicit none
    integer, dimension(:,:), pointer, intent(inout) :: pointr
    integer, dimension(:), intent(in), target :: vector
    integer, intent(in) :: dim1, dim2
        pointr(1:dim1,1:dim2) => vector
    end subroutine reshape_vector2pointer_C

    subroutine reshape_vector2pointer_D(vector, pointr,dim1,dim2,dim3)
    implicit none
    integer, dimension(:,:,:), pointer, intent(inout) :: pointr
    integer, dimension(:), intent(in), target :: vector
    integer, intent(in) :: dim1, dim2, dim3
        pointr(1:dim1,1:dim2,1:dim3) => vector
    end subroutine reshape_vector2pointer_D

#else
    subroutine reshape_vector2pointer_A(vector, pointr, dim1,dim2)
    implicit none
    real, dimension(:,:), pointer, intent(inout) :: pointr
    real, dimension(:), intent(in) :: vector
    integer, intent(in) :: dim1, dim2
        if (.not. associated(pointr)) allocate(pointr(dim1,dim2))
        pointr(1:dim1,1:dim2) = reshape(vector,[dim1,dim2])
    end subroutine reshape_vector2pointer_A

    subroutine reshape_vector2pointer_B(vector, pointr, dim1,dim2,dim3)
    implicit none
    real, dimension(:,:,:), pointer, intent(inout) :: pointr
    real, dimension(:), intent(in) :: vector
    integer, intent(in) :: dim1, dim2, dim3
        if (.not. associated(pointr)) allocate(pointr(dim1,dim2,dim3))
        pointr(1:dim1,1:dim2,1:dim3) = reshape(vector,[dim1,dim2,dim3])
    end subroutine reshape_vector2pointer_B

    subroutine reshape_vector2pointer_C(vector, pointr, dim1,dim2)
    implicit none
    integer, dimension(:,:), pointer, intent(inout) :: pointr
    integer, dimension(:), intent(in) :: vector
    integer, intent(in) :: dim1, dim2
        if (.not. associated(pointr)) allocate(pointr(dim1,dim2))
        pointr(1:dim1,1:dim2) = reshape(vector,[dim1,dim2])
    end subroutine reshape_vector2pointer_C

    subroutine reshape_vector2pointer_D(vector, pointr, dim1,dim2,dim3)
    implicit none
    integer, dimension(:,:,:), pointer, intent(inout) :: pointr
    integer, dimension(:), intent(in) :: vector
    integer, intent(in) :: dim1, dim2, dim3
        if (.not. associated(pointr)) allocate(pointr(dim1,dim2,dim3))
        pointr(1:dim1,1:dim2,1:dim3) = reshape(vector,[dim1,dim2,dim3])
    end subroutine reshape_vector2pointer_D

#endif



end module multiphase_caching
