!    Copyright (C) 2008 Imperial College London and others.
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

module les_viscosity_module
  !!< This module computes a viscosity term to implement LES
  use fields
  implicit none

  private

  public les_viscosity_strength, les_length_scale_tensor

contains

  function les_viscosity_strength(du_t, relu)
    !! Computes the strain rate for the LES model 
    !! derivative of velocity shape function (nloc x ngi x dim)
    real, dimension(:,:,:), intent(in):: du_t
    !! relative velocity (nonl. vel.- grid vel.) (dim x nloc)
    real, dimension(:,:), intent(in):: relu

    real, dimension( size(du_t,2) ):: les_viscosity_strength

    real, dimension(size(du_t,3),size(du_t,3)):: s
    real vis
    integer dim, ngi, nloc
    integer gi, loc, i, j

    nloc=size(du_t,1)
    ngi=size(du_t,2)
    dim=size(du_t,3)

    do gi=1, ngi

       s=0.5*matmul( relu, du_t(:,gi,:) )
       s=s+transpose(s)

       vis=sqrt( 2*sum( s**2 ) )

       les_viscosity_strength(gi)=4.0 * vis

    end do

  end function les_viscosity_strength

  function les_length_scale_tensor(du_t, shape) result(t)
    !! Computes a length scale tensor to be used in LES (units are in length^2)
    !! derivative of velocity shape function (nloc x ngi x dim)
    real, dimension(:,:,:), intent(in):: du_t
    !! the resulting tensor (dim x dim x ngi)
    real, dimension(size(du_t,3),size(du_t,3),size(du_t,2)) :: t
    !! for a simplex if degree==1 the tensor is the same for all gaussian points
    type(element_type), intent(in):: shape

    real, dimension(size(t,1), size(t,2)):: M
    real r
    integer gi, loc, i
    integer dim, ngi, nloc, compute_ngi

    t=0.0

    nloc=size(du_t,1)
    ngi=size(du_t,2)
    dim=size(du_t,3)

    if (shape%degree<=1 .and. shape%numbering%family==FAMILY_SIMPLEX) then
       compute_ngi=1
    else
       compute_ngi=ngi
    end if

    do gi=1, compute_ngi
       do loc=1, nloc
          call outer_product( du_t(loc,gi,:), du_t(loc,gi,:), M)
          r=sum( (/ ( M(i,i), i=1, dim) /) )
          t(:,:,gi)=t(:,:,gi)+M/(r**2)
       end do
    end do

    ! copy the rest
    do gi=compute_ngi+1, ngi
       t(:,:,gi)=t(:,:,1)
    end do

  end function les_length_scale_tensor

end module les_viscosity_module
