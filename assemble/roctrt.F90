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
  !
  !
  !
  !
module roctrt_module
  contains
  subroutine roctrt(c1t,c2t,c3t,findct,colct,&
       &         nct,fredop,nonods,&
       &         nx,ny,nz, t1x,t1y,t1z, t2x,t2y,t2z,  &
       &         nodrot,nnodro,d3)
    use fldebug
    implicit none
    ! this sub works out the matrix ct*r^t=(c1t,c2t,c3t)*r^t
    ! where r is the rotation matrix and puts the result back into ct. 
    ! nodrot contains the nodes to rotate. 
    ! nb the normal is pointing out of the domain. 
    !  the rotation matrix in 3-d is r=  
    !   t1x   t1y   t1z
    !   t2x   t2y   t2z
    !   nx    ny    nz
    !  the rotation matrix in 2-d is r=  
    !   t1x   t1y   
    !   nx    ny    
    ! 
    integer fredop,nonods,nct,nnodro
    real nx(nnodro),ny(nnodro),nz(nnodro)
    real t1x(nnodro),t1y(nnodro),t1z(nnodro)
    real t2x(nnodro),t2y(nnodro),t2z(nnodro)
    integer nodrot(nnodro)
    !
    real c1t(nct),c2t(nct),c3t(nct)
    integer, allocatable, dimension(:) :: tesnod(:) 
    ! if tesnod(nod)=0 then node nod is not to be rotated 
    ! if tesnod(nod)>0 then node nod is to be rotated, with 
    ! rotation index tesnod(nod)
    ! tesnod is real as it uses spare real storage. 
    integer findct(fredop+1),colct(nct)
    logical d3
    integer count,nod,nodcol
    integer i,ii,irow
    real rc1t,rc2t,rc3t

    allocate( tesnod(nonods) )

    ewrite(1,*)'just inside roctrt'
    !     
    do  i=1,nonods           ! was loop 30
       tesnod(i)=0
    end do                   ! was loop 30
    do  ii=1,nnodro          ! was loop 40
       nod=nodrot(ii)
       tesnod(nod)=ii
    end do                   ! was loop 40
    !     
    do  irow=1,fredop        ! was loop 50
       do  count=findct(irow),findct(irow+1)-1 ! was loop 50
          nodcol=colct(count)
          if(tesnod(nodcol).gt.0.5) then
             if(.not.d3) then
                rc1t=c1t(count)
                rc2t=c2t(count)
                ii=tesnod(nodcol)
                c1t(count)=rc1t*t1x(ii) + rc2t*t1y(ii)
                c2t(count)=rc1t*nx(ii)  + rc2t*ny(ii)
             endif
             if(d3) then
                rc1t=c1t(count)
                rc2t=c2t(count)
                rc3t=c3t(count)
                ii=tesnod(nodcol)
                c1t(count)=rc1t*t1x(ii) + rc2t*t1y(ii) + rc3t*t1z(ii)
                c2t(count)=rc1t*t2x(ii) + rc2t*t2y(ii) + rc3t*t2z(ii)
                c3t(count)=rc1t*nx(ii)  + rc2t*ny(ii)  + rc3t*nz(ii)
                !     
                !     c1t(count)=rc1t*t1x(ii) + rc2t*t2x(ii) + rc3t*nx(ii)
                !     c2t(count)=rc1t*t1y(ii) + rc2t*t2y(ii) + rc3t*ny(ii)
                !     c3t(count)=rc1t*t1z(ii)  + rc2t*t2z(ii)  + rc3t*nz(ii)
             endif
          endif
       end do                ! was loop 50
    end do                   ! was loop 50

    deallocate( tesnod )

    ewrite(1,*)'just leaving roctrt'
  end subroutine roctrt
end module roctrt_module

