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

module synthetic_bc

use FLDebug
use generic_interface
use spud
use state_module
use global_parameters, only: dt, option_path_len
use fields
use elements
use parallel_tools
implicit none

private
public synthetic_eddy_method, add_sem_bc, initialise_sem_memory

  
type eddy
   !!< Store eddy info, for synthetic eddy method
   real, dimension(:),pointer:: xeddy
   real, dimension(:),pointer:: yeddy
   real, dimension(:),pointer:: zeddy

   integer, dimension(:),pointer:: eu
   integer, dimension(:),pointer:: ev
   integer, dimension(:),pointer:: ew
end type eddy


integer, save:: sem_bc_count=0
type(eddy),dimension(:),pointer,save:: eddies

contains

  subroutine synthetic_eddy_method(surface_field, surface_field1,  & 
       surface_field2, surface_field3, bc_position, bc_path_i, ns)

      ! declarations
      type(vector_field),intent(inout)     :: surface_field
      type(vector_field),intent(in)        :: surface_field1,surface_field2
      type(vector_field),intent(in)        :: surface_field3
      type(vector_field),intent(in)        :: bc_position
      type(element_type),pointer           :: x_shape,vel_shape
      real,allocatable,dimension(:)        :: detwei,vel_int
      real,allocatable,dimension(:,:),save :: uav
      real,allocatable,dimension(:),save   :: area
      real                                 :: xmin,xmax,ymin,ymax,zmin,zmax,tarea
      real                                 :: lxmin,lxmax,lymin,lymax,lzmin,lzmax
      real                                 :: xminb,xmaxb,yminb,ymaxb,zminb,zmaxb
      real                                 :: dxb,dyb,dzb
      real                                 :: rdmax,rts_x,rts_y,rts_z
      real                                 :: dx,dy,dz,vol
      real                                 :: dxp,dyp,dzp
      integer                              :: nsem,ns,nots,bcnod,i,j,ele,ratio
      real,allocatable,dimension(:,:)      :: cords
      real,dimension(2)                    :: tmpcord
      real,dimension(:),pointer            :: x,y,z
      real,dimension(:),pointer            :: lx,ly,lz
      real                                 :: uf,vf,wf,ff
      real,allocatable,dimension(:)        :: ufl,vfl,wfl
      real                                 :: a11,a21,a22,a31,a32,a33
      real                                 :: resuu,resvv,resww
      real,dimension(3)                    :: vl
      real                                 :: minum,minvm,minwm,maxum,maxvm,maxwm
      real                                 :: reuumn,reuumx,revvmn,revvmx,rewwmn,rewwmx
      real                                 :: uflmn,uflmx,vflmn,vflmx,wflmn,wflmx
      logical,save                         :: initz=.false.
      logical,allocatable,dimension(:),save:: initeddy
      character(len=OPTION_PATH_LEN)       :: bc_path_i

#ifdef HAVE_MPI
      INCLUDE 'mpif.h'
      INTEGER IERR
#endif

      ! time varying turbulent inlet conditions
      ! input:: mean flow, Reynolds stress tensor, turbulence length-scale
      ! synthetic eddy method (see Jarrin et al. 2006)
      ! can be used to generate turbulent inlet bcs for any given surface
      ! aligned with the cartesian coordinate system

      ! known issues::
      ! 1.there is a bug when calculating the inflow plane surface area in parallel-tarea
      ! 2.there is a bug when calculating the inflow bulk velocity in parallel-uav
      ! both because of the halo

      ! field names::
      ! surface_field  - inlet boundary condition (calculated here)
      ! surface_field1 - mean velocity
      ! surface_field2 - Re_ij (uu,vv,ww)
      ! surface_field3 - turbulence length scale
      ! bc_position    - inlet plane mesh

      ewrite(3,*) 'setting a turbulent inlet boundary condition',ns
      nsem = sem_bc_count
      ! get additional turbulence related options
      call get_option(trim(bc_path_i),nots)

      ! calculate min&max turbulence lengthscale
      lx => surface_field3%val(1)%ptr  
      ly => surface_field3%val(2)%ptr
      lz => surface_field3%val(3)%ptr

      lxmin = minval(lx); lxmax = maxval(lx)
      lymin = minval(ly); lymax = maxval(ly)
      lzmin = minval(lz); lzmax = maxval(lz)

      call allmin(lxmin); call allmax(lxmax)
      call allmin(lymin); call allmax(lymax)
      call allmin(lzmin); call allmax(lzmax)

      ewrite(3,*) 'number of eddies',nots
      ewrite(3,*) 'turbulence lengthscale min&max',lxmin,lxmax,':',lymin,lymax,':',lzmin,lzmax

      ! number of nodes on boundary surface
      bcnod=surface_field%mesh%nodes

      ! coordinates of nodes
      x => bc_position%val(1)%ptr  
      y => bc_position%val(2)%ptr
      z => bc_position%val(3)%ptr

      ! work out min & max of boundary surface and orientation
      xmin = minval(x); xmax = maxval(x)
      ymin = minval(y); ymax = maxval(y)
      zmin = minval(z); zmax = maxval(z)

      call allmin(xmin)
      call allmin(ymin)
      call allmin(zmin)

      call allmax(xmax)
      call allmax(ymax)
      call allmax(zmax)

      dxp = xmax-xmin
      dyp = ymax-ymin
      dzp = zmax-zmin

      ! generate bounding box
      if (dxp<dyp .and. dxp<dzp) then
         ewrite(3,*) 'using a yz plane'

         xminb=-lxmax+xmax
         xmaxb= lxmax+xmax
         yminb= ymin-lymax
         ymaxb= ymax+lymax
         zminb= zmin-lzmax
         zmaxb= zmax+lzmax

         rdmax=lxmax

      elseif (dyp<dxp .and. dyp<dzp) then
         ewrite(3,*) 'using a xz plane'

         xminb= xmin-lxmax
         xmaxb= xmax+lxmax
         yminb=-lymax+ymax
         ymaxb= lymax+ymax
         zminb= zmin-lzmax
         zmaxb= zmax+lzmax

         rdmax=lymax

      elseif (dzp<dxp .and. dzp<dyp) then
         ewrite(3,*) 'using a xy plane'

         xminb= xmin-lxmax
         xmaxb= xmax+lxmax
         yminb= ymin-lymax
         ymaxb= ymax+lymax
         zminb=-lzmax+zmax
         zmaxb= lzmax+zmax

         rdmax=lzmax
      end if

      dxb= xmaxb-xminb
      dyb= ymaxb-yminb
      dzb= zmaxb-zminb

      ! generate initial coordinades, signs of turbulent spots
      ! calculate convection velocity
      if (.not.initz) then
         allocate(initeddy(nsem))
         initeddy = .true.

         allocate(uav(nsem,3))
         allocate(area(nsem))
         uav=0.; area=0.

         initz=.true.
      endif

      if(initeddy(ns))then
         if(GetRank()==0) then
            allocate(cords(nots,3))
            call random_seed()
            call random_number(cords)
            
            do i=1,nots
               eddies(ns)%xeddy(i)=xminb+(xmaxb-xminb)*cords(i,1)
               eddies(ns)%yeddy(i)=yminb+(ymaxb-yminb)*cords(i,2)
               eddies(ns)%zeddy(i)=zminb+(zmaxb-zminb)*cords(i,3)

               eddies(ns)%eu(i)=esign()
               eddies(ns)%ev(i)=esign()
               eddies(ns)%ew(i)=esign()
            enddo
            deallocate(cords)
         end if

#ifdef HAVE_MPI
         if(IsParallel())then
            call mpi_bcast(eddies(ns)%xeddy(1:nots),nots,GetPREAL(),0,mpi_comm_world,ierr)
            call mpi_bcast(eddies(ns)%yeddy(1:nots),nots,GetPREAL(),0,mpi_comm_world,ierr)
            call mpi_bcast(eddies(ns)%zeddy(1:nots),nots,GetPREAL(),0,mpi_comm_world,ierr)

            call mpi_bcast(eddies(ns)%eu(1:nots),nots,mpi_integer,0,mpi_comm_world,ierr)
            call mpi_bcast(eddies(ns)%ev(1:nots),nots,mpi_integer,0,mpi_comm_world,ierr)
            call mpi_bcast(eddies(ns)%ew(1:nots),nots,mpi_integer,0,mpi_comm_world,ierr)
         endif
#endif

         x_shape   => ele_shape(bc_position,1)
         vel_shape => ele_shape(surface_field1,1)

         allocate(detwei(1:ele_ngi(bc_position,1)))
         allocate(vel_int(1:ele_loc(surface_field1,1)))

         tarea=0.
         do ele=1,element_count(surface_field1)
            call transform_to_physical(bc_position,ele,detwei=detwei)
            vel_int=shape_rhs(vel_shape, detwei)
            tarea=tarea+sum(vel_int)
            do i=1,surface_field1%dim
               uav(ns,i)=uav(ns,i)+sum(ele_val(surface_field1,i,ele)*vel_int)
            enddo
         enddo

         deallocate(detwei,vel_int)

         call allsum(tarea)
         call allsumv(uav(ns,:))

         do i=1,surface_field1%dim
            uav(ns,i)=uav(ns,i)/tarea
         enddo

         area(ns)=tarea
         initeddy(ns)=.false.
      endif
      
      ewrite(3,*) 'surface info'
      ewrite(3,*) 'x min&max',xmin,xmax,dxp
      ewrite(3,*) 'y min&max',ymin,ymax,dyp
      ewrite(3,*) 'z min&max',zmin,zmax,dzp
      ewrite(3,*) 'surface area',area(ns)
      ewrite(3,*) 'number of nodes',bcnod
      ewrite(3,*) 'convection velocity:',uav(ns,1),uav(ns,2),uav(ns,3)

      ! bounding box volume
      vol=area(ns)*2*rdmax

      ewrite(3,*) 'bounding box info'
      ewrite(3,*) 'x min&max',xminb,xmaxb
      ewrite(3,*) 'y min&max',yminb,ymaxb
      ewrite(3,*) 'z min&max',zminb,zmaxb
      ewrite(3,*) 'box volume',vol

      if(GetRank()==0) then
         ! convect turbulent spots
         do i=1,nots
            eddies(ns)%xeddy(i)=eddies(ns)%xeddy(i)+uav(ns,1)*dt
            eddies(ns)%yeddy(i)=eddies(ns)%yeddy(i)+uav(ns,2)*dt
            eddies(ns)%zeddy(i)=eddies(ns)%zeddy(i)+uav(ns,3)*dt
         enddo

         ! regenerate
         do i=1,nots

            if(eddies(ns)%xeddy(i)>xmaxb)then
               call random_number(tmpcord)

               ratio= int((abs(eddies(ns)%xeddy(i)-xminb))/dxb)

               eddies(ns)%xeddy(i)=eddies(ns)%xeddy(i)-ratio*dxb
               eddies(ns)%yeddy(i)=yminb+(ymaxb-yminb)*tmpcord(1)
               eddies(ns)%zeddy(i)=zminb+(zmaxb-zminb)*tmpcord(2)

               eddies(ns)%eu(i)=esign()
               eddies(ns)%ev(i)=esign()
               eddies(ns)%ew(i)=esign()

            elseif(eddies(ns)%xeddy(i)<xminb)then
               call random_number(tmpcord)

               ratio= int((abs(eddies(ns)%xeddy(i)-xmaxb))/dxb)

               eddies(ns)%xeddy(i)=eddies(ns)%xeddy(i)+ratio*dxb
               eddies(ns)%yeddy(i)=yminb+(ymaxb-yminb)*tmpcord(1)
               eddies(ns)%zeddy(i)=zminb+(zmaxb-zminb)*tmpcord(2)

               eddies(ns)%eu(i)=esign()
               eddies(ns)%ev(i)=esign()
               eddies(ns)%ew(i)=esign()

            elseif(eddies(ns)%yeddy(i)>ymaxb)then
               call random_number(tmpcord)

               ratio= int((abs(eddies(ns)%yeddy(i)-yminb))/dyb)

               eddies(ns)%xeddy(i)=xminb+(xmaxb-xminb)*tmpcord(1)
               eddies(ns)%yeddy(i)=eddies(ns)%yeddy(i)-ratio*dyb
               eddies(ns)%zeddy(i)=zminb+(xmaxb-xminb)*tmpcord(2)

               eddies(ns)%eu(i)=esign()
               eddies(ns)%ev(i)=esign()
               eddies(ns)%ew(i)=esign()

            elseif(eddies(ns)%yeddy(i)<yminb)then
               call random_number(tmpcord)

               ratio= int((abs(eddies(ns)%yeddy(i)-ymaxb))/dyb)

               eddies(ns)%xeddy(i)=xmaxb+(xmaxb-xminb)*tmpcord(1)
               eddies(ns)%yeddy(i)=eddies(ns)%yeddy(i)+ratio*dyb
               eddies(ns)%zeddy(i)=zminb+(zmaxb-zminb)*tmpcord(2)

               eddies(ns)%eu(i)=esign()
               eddies(ns)%ev(i)=esign()
               eddies(ns)%ew(i)=esign()

            elseif(eddies(ns)%zeddy(i)>zmaxb)then
               call random_number(tmpcord)

               ratio= int((abs(eddies(ns)%zeddy(i)-zminb))/dzb)

               eddies(ns)%xeddy(i)=xminb+(xmaxb-xminb)*tmpcord(1)
               eddies(ns)%yeddy(i)=yminb+(ymaxb-yminb)*tmpcord(2)
               eddies(ns)%zeddy(i)=eddies(ns)%zeddy(i)-ratio*dzb

               eddies(ns)%eu(i)=esign()
               eddies(ns)%ev(i)=esign()
               eddies(ns)%ew(i)=esign()

            elseif(eddies(ns)%zeddy(i)<zminb)then
               call random_number(tmpcord)

               ratio= int((abs(eddies(ns)%zeddy(i)-zmaxb))/dzb)

               eddies(ns)%xeddy(i)=xmaxb+(xmaxb-xminb)*tmpcord(1)
               eddies(ns)%yeddy(i)=yminb+(ymaxb-yminb)*tmpcord(2)
               eddies(ns)%zeddy(i)=eddies(ns)%zeddy(i)+ratio*dzb
               
               eddies(ns)%eu(i)=esign()
               eddies(ns)%ev(i)=esign()
               eddies(ns)%ew(i)=esign()
            endif
         enddo
      end if

#ifdef HAVE_MPI
      if(IsParallel())then
         call mpi_bcast(eddies(ns)%xeddy(1:nots),nots,GetPREAL(),0,mpi_comm_world,ierr)
         call mpi_bcast(eddies(ns)%yeddy(1:nots),nots,GetPREAL(),0,mpi_comm_world,ierr)
         call mpi_bcast(eddies(ns)%zeddy(1:nots),nots,GetPREAL(),0,mpi_comm_world,ierr)

         call mpi_bcast(eddies(ns)%eu(1:nots),nots,mpi_integer,0,mpi_comm_world,ierr)
         call mpi_bcast(eddies(ns)%ev(1:nots),nots,mpi_integer,0,mpi_comm_world,ierr)
         call mpi_bcast(eddies(ns)%ew(1:nots),nots,mpi_integer,0,mpi_comm_world,ierr)        
      endif
#endif 

      !ewrite(3,*) 'eddy coordinates & signs'
      !do i=1,nots
      !   ewrite(3,*) i,                                                    &
      !        eddies(ns)%xeddy(i),eddies(ns)%yeddy(i),eddies(ns)%zeddy(i), & 
      !        eddies(ns)%eu(i),eddies(ns)%ev(i),eddies(ns)%ew(i)
      !enddo
      
      ! calculate fluctuating component
      allocate(ufl(bcnod))
      allocate(vfl(bcnod))
      allocate(wfl(bcnod))

      ufl=0.; vfl=0.; wfl=0.

      do i=1,bcnod

         uf=0.; vf=0.; wf=0.
         rts_x=node_val(surface_field3,1,i)
         rts_y=node_val(surface_field3,2,i)
         rts_z=node_val(surface_field3,3,i)

         do j=1,nots

            dx=abs(x(i)-eddies(ns)%xeddy(j))
            dy=abs(y(i)-eddies(ns)%yeddy(j))
            dz=abs(z(i)-eddies(ns)%zeddy(j))

            if (dx<rts_x .and. dy<rts_y .and. dz<rts_z) then

               ff=(1.-dx/rts_x)*(1.-dy/rts_y)*(1.-dz/rts_z)
               ff=ff/((sqrt(3./2.*rts_x))*(sqrt(3./2.*rts_y))*(sqrt(3./2.*rts_z)))

               uf=uf+eddies(ns)%eu(j)*ff
               vf=vf+eddies(ns)%ev(j)*ff
               wf=wf+eddies(ns)%ew(j)*ff

            end if
         end do

         uf=uf*sqrt(vol/nots)
         vf=vf*sqrt(vol/nots)
         wf=wf*sqrt(vol/nots)

         ! Cholesky decomposition of the Re_ij
         resuu=node_val(surface_field2,1,i)
         resvv=node_val(surface_field2,2,i)
         resww=node_val(surface_field2,3,i)

         if (resuu<0.)then
            ewrite(3,*) 'WARNING: there is a negative value in Re_uu'
            resuu=1.e-5
         end if
         if (resvv<0.)then
            ewrite(3,*) 'WARNING: there is a negative value in Re_vv'
            resvv=1.e-5
         end if
         if (resww<0.)then
            ewrite(3,*) 'WARNING: there is a negative value in Re_ww'
            resww=1.e-5
         end if
         
         a11 = sqrt(resuu)
         a21 = 0.
         a22 = sqrt(resvv)
         a31 = 0.
         a32 = 0.
         a33 = sqrt(resww)

         ufl(i) = a11* uf
         vfl(i) = a21*uf + a22*vf
         wfl(i) = a31*uf + a32*vf + a33*wf
      end do

      ! get min&max of mean velocities
      minum=minval(surface_field1%val(1)%ptr); maxum=maxval(surface_field1%val(1)%ptr)
      minvm=minval(surface_field1%val(2)%ptr); maxvm=maxval(surface_field1%val(2)%ptr)
      minwm=minval(surface_field1%val(3)%ptr); maxwm=maxval(surface_field1%val(3)%ptr)

      call allmin(minum); call allmax(maxum)
      call allmin(minvm); call allmax(maxvm)
      call allmin(minwm); call allmax(maxwm)

      ! calculate final velocity signal
      do i = 1, bcnod
         vl(1) = node_val(surface_field1,1,i)+ufl(i)
         vl(2) = node_val(surface_field1,2,i)+vfl(i)
         vl(3) = node_val(surface_field1,3,i)+wfl(i)
         call set(surface_field,i,vl)
      end do

      ! calculate min&max of Re_ij
      reuumn = minval(surface_field2%val(1)%ptr); reuumx=maxval(surface_field2%val(1)%ptr)
      revvmn = minval(surface_field2%val(2)%ptr); revvmx=maxval(surface_field2%val(2)%ptr)
      rewwmn = minval(surface_field2%val(3)%ptr); rewwmx=maxval(surface_field2%val(3)%ptr)

      call allmin(reuumn); call allmax(reuumx)
      call allmin(revvmn); call allmax(revvmx)
      call allmin(rewwmn); call allmax(rewwmx)

      ! calculate min&max of fluctuating component
      uflmn=minval(ufl); uflmx=maxval(ufl)
      vflmn=minval(vfl); vflmx=maxval(vfl)
      wflmn=minval(wfl); wflmx=maxval(wfl)

      call allmin(uflmn); call allmax(uflmx)
      call allmin(vflmn); call allmax(vflmx)
      call allmin(wflmn); call allmax(wflmx)

      ewrite(3,*) 'inlet condition diagnostics'
      ewrite(3,*) 'Re stresses min&max'
      ewrite(3,*) 'Re_uu',reuumn,':',reuumx
      ewrite(3,*) 'Re_vv',revvmn,':',revvmx
      ewrite(3,*) 'Re_ww',rewwmn,':',rewwmx
      ewrite(3,*) 'mean:fluctuating min&max'
      ewrite(3,*) 'u',minum,maxum,':',uflmn,uflmx
      ewrite(3,*) 'v',minvm,maxvm,':',vflmn,vflmx
      ewrite(3,*) 'w',minwm,maxwm,':',wflmn,wflmx

      deallocate(ufl,vfl,wfl)

      ewrite(3,*) 'leaving turbulent inlet boundary conditions routine'
    end subroutine synthetic_eddy_method

    !----------------------------------------------------------

    function esign()

      integer :: esign
      real    :: i

      call random_number(i)

      if (i<0.5) then
        esign = -1        
      else
        esign = +1
      endif

      return 
     end function esign
     
     !----------------------------------------------------------
     
     subroutine initialise_sem_memory(ns,nots)

       logical,save::initialise_memory=.false.
       logical,allocatable,dimension(:),save::initeddymem
       integer:: ns, nots, nsem

       nsem=sem_bc_count

       if(.not.initialise_memory)then
          allocate(initeddymem(nsem))
          initeddymem=.true.
          allocate(eddies(nsem))
          initialise_memory=.true.
       endif
       
       if(initeddymem(ns))then
          allocate(eddies(ns)%xeddy(nots));allocate (eddies(ns)%yeddy(nots));allocate(eddies(ns)%zeddy(nots))
          allocate(eddies(ns)%eu(nots))   ;allocate(eddies(ns)%ev(nots))    ;allocate(eddies(ns)%ew(nots))
          initeddymem(ns)=.false.
       endif
       
     end subroutine initialise_sem_memory
     
     !----------------------------------------------------------
     
     subroutine add_sem_bc(have_sem_bc)

       logical:: have_sem_bc

       if (have_sem_bc) then
          sem_bc_count=sem_bc_count+1
       end if
       
     end subroutine add_sem_bc
     
     
   end module synthetic_bc
   
