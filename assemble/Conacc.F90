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

module conacc_module

  use fldebug
  
  implicit none
  
  private
  
  public :: conacc

contains

  SUBROUTINE CONACC(BCT1W,BCT1,BCT2,&
       &              T,DT,&
       &              cleart)
    ! Convert the boundary conditions to acceleration from a specified value.
    ! BCT1W=is the working acceleration actually used
    ! in NAVSTO or ADVDIF.

    REAL, dimension(:), intent(inout) :: BCT1W, bct1
    integer, dimension(:), intent(inout) :: bct2
    REAL, intent(in) :: DT
    REAL, intent(inout), dimension(:) :: T
    logical cleart
    !local variables
    INTEGER :: NOBCT,NONODS,I
    
    ewrite(1, *) "In conacc"

    nobct = size(bct1w)
    assert(size(bct1w)==size(bct1))
    assert(size(bct2)==size(bct1))
    nonods = size(T)

    !IF((TTPERI/100)==-666) THEN
    ! Constant boundary condition case
    if(cleart) then
       do i=1,nobct
          T(BCT2(I))=BCT1(I)
       end do
       cleart=.false.
       BCT1W= 0.0
    else
       do i=1,nobct
          BCT1W(i)= (BCT1(I)-T(BCT2(I)))/DT
       end do
    end if
    !else
    !     FLAbort('not supported outside new options')
    ! end IF

  end subroutine conacc

  Subroutine tideBC(BCT1W,BCT2,T,NONOD,NOBCT,DT,ACCTIM,TTPER2,x,y,z)
    ! for implementing tide elevation OBC
    !                   H. Liu, 23 Jun, 2007
    Integer, intent(in) :: nonod
    Integer, intent(inout) :: nobct
    Integer, intent(inout) :: bct2(nobct)
    Real, intent(in) :: t(nonod)
    Real, intent(in) :: x(nonod), y(nonod), z(nonod)
    Real, intent(in) :: ttper2, acctim, dt
    Real, intent(out) :: bct1w(nobct)

    Real, Allocatable :: tamp(:,:,:), tphase(:,:,:)
    Integer, Allocatable :: obcmask(:,:)
    Real, save, Allocatable :: amp(:,:), phase(:, :), tomega(:)
    Integer, save, Allocatable :: tmask(:)

    Real, Allocatable :: long(:), lat(:)

    Real :: xmin, ymin, xmax, ymax, dx, dy
    Real :: xtmp, ytmp
    Real :: xtmp11, ytmp11, xtmp22, ytmp22
    Real :: tmp1, tmp2
    Real :: acc, rrsin, temp
    Integer :: i, j, k, nod, n, m
    Integer :: istart, iend, istride
    Integer :: jstart, jend, jstride
    Integer :: ntmp, mtmp, ntmpp1, mtmpp1
    Integer :: sphere

    Integer, save :: ntc
    logical, save::initialised=.false.
    Logical :: mask=.true.
    !         Logical :: readinobc=.true.
    Logical :: readinobc=.false.

    integer ierr
    
    ewrite(1, *) "In tidebc"

    ewrite(3,*)  'intitialised==', initialised
    ewrite(3,*)  'nobct==', nobct

    if(.NOT.initialised) then
       ewrite(3,*)  'max x,y,z--->',maxval(x), maxval(y), maxval(z)
       ewrite(3,*)  'min x,y,z--->',minval(x), minval(y), minval(z)

       !      test the input obc type       
       inquire(file='obcnode.dat',exist=readinobc)

       if(.not.readinobc) then
          ewrite(1,*) 'start allocate long and lat'
          allocate(long(nobct), lat(nobct))
          ewrite(1,*) 'allocate finished'
          ewrite(1,*) 'open tideBC_grid.dat now'
          open(1,file='tideBC_grid.dat',status='old')
          rewind(1)
          read(1,*) sphere
          read(1,*) n,m,ntc
          read(1,*) xmin, ymin, xmax, ymax, dx, dy

          do k=1,ntc! Was loop 
             read(1,*) 
          end do

          do j=1,2
          do i=1,n
             read(1,*,iostat=ierr) xtmp, ytmp
             if(i.eq.1.and.j.eq.1) then
               xtmp11=xtmp
               ytmp11=ytmp
             elseif(i.eq.2.and.j.eq.2) then
               xtmp22=xtmp
               ytmp22=ytmp
             end if
           end do
           end do
           close(1)

          if(xtmp22.gt.xtmp11) then
            istart=1
            iend=n
            jstride=1
          else
            istart=n
            iend=1
            istride=-1
          endif

          if(ytmp22.gt.ytmp11) then
            jstart=1
            jend=m
            jstride=1
          else
            jstart=m
            jend=1
            jstride=-1
          endif


          open(1,file='tideBC_grid.dat',status='old')
          rewind(1)
          read(1,*) sphere
          read(1,*) n,m,ntc
          read(1,*) xmin, ymin, xmax, ymax, dx, dy

          ewrite(3,*) 'xmin, ymin, xmax, ymax, dx, dy'
          ewrite(3,*) xmin, ymin, xmax, ymax, dx, dy

          if(abs(xmax-xmin-(n-1)*dx).gt.0.01*dx) then
             ewrite(0,*) 'xmin, xmax, dx, xmax-xmin-(n-1)*dx'
             ewrite(0,*) xmin, xmax, dx, xmax-xmin-(n-1)*dx
             FLAbort ("xmax or dx wrong in tideBC_grid.dat!!!")
          end if
          if(abs(ymax-ymin-(m-1)*dy).gt.0.01*dy) then
             ewrite(0,*) 'ymin, ymax, dy, ymax-ymin-(m-1)*dy'
             ewrite(0,*) ymin, ymax, dy, ymax-ymin-(m-1)*dy
             FLAbort ("ymax or dy wrong in tideBC_grid.dat!!!")
          end if

          ewrite(2,*) 'allocate1'
          allocate (amp(nobct,ntc), phase(nobct,ntc), tmask(nobct))
          ewrite(2,*) 'allocate2'
          allocate (tomega(ntc),tamp(n,m,ntc), tphase(n,m,ntc), obcmask(n,m))
          ewrite(2,*) 'finished allocate2'

          do k=1,ntc! Was loop 
             read(1,*) tomega(k)
          end do

          do j=jstart,jend,jstride
             do i=istart,iend,istride
!          do j=1,m! Was loop 
!             do i=1,n! Was loop 
                read(1,*,iostat=ierr) xtmp, ytmp, obcmask(i,j), (tamp(i,j,k),&
                     &                                tphase(i,j,k), k=1,ntc)
                if(ierr.ne.0) then
                   FLAbort("tideBC_grid.dat is foobar")
                end if
                !                  do k=1,ntc
                !                     obcmask(i,j,k)=1
                !                  end do
             end do
          end do
          close(1)

          ewrite(3,*)  'finished reading tideBC_grid.dat'
          ewrite(3,*)  'sphere===', sphere
          ewrite(3,*)  'xmin, ymin==', xmin, ymin

          if(sphere.eq.1) then
             do i=1,nobct! Was loop 
                nod=bct2(i)
                long(i)=atan2(y(nod),x(nod))
                long(i)=long(i)*180./3.14159265
                lat(i)=asin(z(nod)/sqrt(x(nod)**2+y(nod)**2+z(nod)**2))
                lat(i)=lat(i)*180./3.14159265
             end do
          else
             do i=1,nobct! Was loop 
                nod=bct2(i)
                long(i)=x(nod)
                lat(i)=y(nod)
             end do
          end if


          do i=1,nobct! Was loop 
             nod=bct2(i)
             ntmp=int((long(i)-xmin)/dx)
             mtmp=int((lat(i)-ymin)/dy)

             if(ntmp.ge.n.or.mtmp.ge.m) then
                ewrite(3,*)  'i,nod,ntmp, mtmp==',i,nod,ntmp,mtmp
                ewrite(3,*)  'long(i), lat(i)==',long(i), lat(i)
                ewrite(3,*)  'x,y,z==',x(nod),y(nod),z(nod)
                FLAbort("wrong data range in (tideBC_grid.dat)")
             end if

             if(ntmp.lt.1.or.mtmp.lt.1) then
                ewrite(3,*)  'i,nod,ntmp, mtmp==',i,nod,ntmp,mtmp
                ewrite(3,*)  'long(i), lat(i)==',long(i), lat(i)
                ewrite(3,*)  'x,y,z==',x(nod),y(nod),z(nod)
                FLAbort("wrong data range in (tideBC_grid.dat)")
             end if

             ntmpp1=ntmp+1
             mtmpp1=mtmp+1


             xtmp=long(i)-xmin-ntmp*dx
             ytmp=lat(i)-ymin-mtmp*dy
             do k=1,ntc! Was loop 
                tmp1=tamp(ntmp,mtmp,k)+(tamp(ntmpp1,mtmp,k)&
                     &                  -tamp(ntmp,mtmp,k))/dx*xtmp
                tmp2=tamp(ntmp,mtmpp1,k)+(tamp(ntmpp1,mtmpp1,k)&
                     &                  -tamp(ntmp,mtmpp1,k))/dx*xtmp
                amp(i,k)=tmp1+(tmp2-tmp1)*ytmp/dy

                tmp1=tphase(ntmp,mtmp,k)+(tphase(ntmpp1,mtmp,k)&
                     &                  -tphase(ntmp,mtmp,k))/dx*xtmp
                tmp2=tphase(ntmp,mtmpp1,k)+(tphase(ntmpp1,mtmpp1,k)&
                     &                  -tphase(ntmp,mtmpp1,k))/dx*xtmp
                phase(i,k)=tmp1+(tmp2-tmp1)*ytmp/dy

                !               ewrite(3,*)  'i,nod, k,ntmp ==', i,nod,k,ntmp
                !               ewrite(3,*)  'long, lat==', long(i), lat(i)
                !               ewrite(3,*)  ' amp(i,k),phase(i,k)==', amp(i,k),phase(i,k)
             end do

             tmask(i)=obcmask(ntmp,mtmp)  +obcmask(ntmpp1,mtmp)+&
                  &                    obcmask(ntmp,mtmpp1)+obcmask(ntmpp1,mtmpp1)

             if((tmask(i).ne.4).and.(tmask(i).ne.-4)) then
                ewrite(3,*) 'ntmp,mtmp,ntmp+1, mtmp+1==',ntmp,mtmp,ntmpp1,mtmpp1
                ewrite(3,*) 'obcmask(n,m)==',obcmask(ntmp,mtmp)
                ewrite(3,*) 'obcmask(n+1,m)==',obcmask(ntmpp1,mtmp)
                ewrite(3,*) 'obcmask(n,m+1)==',obcmask(ntmp,mtmpp1)
                ewrite(3,*) 'obcmask(n+1,m+1)==',obcmask(ntmpp1,mtmpp1)
                ewrite(3,*)  'tmask(i)==', tmask(i)
                FLAbort("tide mask field wrong in (tideBC_grid.dat)")
             end if
          end do

          if(mask) then
             j=0
             do i=1,nobct! Was loop 
                ewrite(3,*) i,bct2(i),x(bct2(i)), y(bct2(i)), tmask(i)
                if(tmask(i).eq.4) then
                   j=j+1
                   bct2(j)=bct2(i)  
                   do k=1,ntc! Was loop 
                      amp(j,k)=amp(i,k)
                      phase(j,k)=phase(i,k)
                   end do
                else if(tmask(i).eq.-4) then

                else
                   ewrite(3,*) 'wrong tmask value, i, nod, tmask==',&
                        &                     i,bct2(i), tmask(i) 
                end if
                nobct=j
             end do
          end if

          deallocate(long, lat)
          deallocate (tamp, tphase)
          deallocate (obcmask)
       else
          ewrite(1,*) 'open obcnode.dat now'
          open(1,file='obcnode.dat',status='old')
          rewind(1)

          read(1,*) nobct,ntc

          allocate (tomega(ntc))
          allocate (amp(nobct*2,ntc), phase(nobct*2,ntc))

          do k=1,ntc! Was loop 
             read(1,*) tomega(k)
          end do

          do i=1,nobct! Was loop 
             read(1,*) bct2(i),(amp(i,k),phase(i,k),k=1,ntc)
          end do
          close(1)

          do i=1, nobct! Was loop 
             bct2(i+nobct)=bct2(i)-nonod/2
             do k=1,ntc! Was loop 
                amp(i+nobct,k)=amp(i,k)
                phase(i+nobct,k)=phase(i,k)
             end do
          end do

          nobct=nobct+nobct
       end if
    end if

    !         if(.not.initialised) then
    !            open(8888,file='obc.dat',status='unknown')
    !            rewind(8888)
    !            write(8888,*) 'nobct=',nobct
    !            write(8888,*) 'i, nod, nodbot, x, y,z, zbot, temp, acc'
    !         end if

    do i=1,nobct! Was loop 
       nod=bct2(i)
       rrsin=0.
       do k=1,ntc! Was loop 
          rrsin=rrsin+amp(i,k)*cos(tomega(k)*(acctim+dt)-phase(i,k))
       end do
       temp=rrsin+ttper2

       !          correction for wetting/drying obc
       !           only for the Mersey estuary test
       !           nodbot=mod(nod, nonod/2)
       !           temp=max(temp,z(nodbot)+0.01)
       !          end correction

       acc=(temp-t(nod))/dt
       !           write(8888,*) i, nod, nodbot, x(nod), y(nod),z(nod), z(nodbot), temp, acc
       bct1w(i)=acc
       !           ewrite(3,*)  'i, nod, rrsin, ttper2, acc==',i,nod, rrsin, ttper2, acc
    end do

    initialised = .true.
    
    ewrite(1, *) "Exiting tidebc"
    
  end subroutine tidebc

end module conacc_module
