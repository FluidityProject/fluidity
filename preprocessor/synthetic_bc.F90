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

module synthetic_bc

use fldebug
use spud
use global_parameters, only: dt, option_path_len
use elements
use parallel_tools
use transform_elements
use fetools
use fields
use state_module

implicit none

private
public synthetic_eddy_method, add_sem_bc, initialise_sem_memory

  
type eddy
!!< Store eddy info, for synthetic eddy method
   real, dimension(:,:),    pointer :: eddy_positions

   integer, dimension(:,:), pointer :: eddy_polarities
end type eddy


integer, save:: sem_bc_count=0
type(eddy),dimension(:),pointer,save:: eddies

contains

  subroutine synthetic_eddy_method(eddy_inlet_velocity, mean_velocity,  & 
       reynolds_stresses, turbulence_lengthscales, bc_position, bc_path_i, ns)
    
      ! declarations
      type(vector_field),intent(inout)     :: eddy_inlet_velocity
      type(vector_field),intent(in)        :: mean_velocity, reynolds_stresses
      type(vector_field),intent(in)        :: turbulence_lengthscales
      type(vector_field),intent(in)        :: bc_position
      character(len=*),  intent(in)        :: bc_path_i
      integer, intent(in)                  :: ns

      type(element_type),pointer           :: x_shape,vel_shape

      real, dimension(ele_ngi(bc_position,1)) :: detwei
      real, dimension(ele_loc(mean_velocity,1)) :: vel_int
      real,allocatable,dimension(:,:),save :: averaged_velocity
      real,allocatable,dimension(:),save   :: area       

      real, dimension(mean_velocity%dim) :: xmin, xmax, lxmin, lxmax, &
           xminb, xmaxb, dxb, rts, dx, dxp
      real                                   :: rdmax, volume, total_area
      
      integer :: number_of_eddies, i, j, ele, ratio, dim
      integer, dimension(1)                :: idxp
      real,allocatable,dimension(:,:)      :: coords
      real,dimension(mean_velocity%dim):: tmpcord

      type(vector_field)                   :: ufl

      real,dimension(mean_velocity%dim):: uf, res, re_ii_min, re_ii_max,&
           uflmin, uflmax, min_mean_velocity, max_mean_velocity, ff
      
      real,dimension(mean_velocity%dim*(mean_velocity%dim+1)/2):: a

      logical,save                         :: initz=.false.
      logical,allocatable,dimension(:),save:: initeddy

#ifdef HAVE_MPI
      include 'mpif.h'
      integer ierr
#endif

      ! time varying turbulent inlet conditions
      ! input:: mean flow, Reynolds stress tensor, turbulence length-scale
      ! synthetic eddy method (see Jarrin et al. 2006)
      ! can be used to generate turbulent inlet bcs for any given surface
      ! aligned with the cartesian coordinate system


      ! field names::
      ! eddy_inlet_velocity  - inlet boundary condition (calculated here)
      ! mean_velocity - mean velocity
      ! reynolds_stresses - Re_ij (uu,vv,ww)
      ! turbulence_lengthscales - turbulence length scale
      ! bc_position    - inlet plane mesh

      ewrite(3,*) 'setting a turbulent inlet boundary condition',ns
      dim=mean_velocity%dim
      ! get additional turbulence related options
      call get_option(trim(bc_path_i),number_of_eddies)

      lxmin = minval(turbulence_lengthscales%val,dim=2)
      lxmax = maxval(turbulence_lengthscales%val,dim=2)

      do i=1,dim
         call allmin(lxmin(i))
         call allmax(lxmax(i))
      end do

      ewrite(3,*) 'number of eddies',number_of_eddies
      ewrite(3,*) 'turbulence lengthscale min&max',lxmin,':',lxmax

      ! work out min & max of boundary surface and calculate orientation
      xmin = minval(bc_position%val,dim=2)
      xmax = maxval(bc_position%val,dim=2)

      do i=1,dim
         call allmin(xmin(i))
         call allmax(xmax(i))
      end do

      dxp= xmax-xmin

      ! orient box of eddies to feed into domain along shortext axis.
      ! This assumes domain is oriented with coordinate system.
      idxp = minloc(dxp)
      select case(idxp(1))
      case(1)
         ewrite(3,*) 'using a yz plane'

         xminb(1)=-lxmax(1)+xmax(1)
         xmaxb(1)= lxmax(1)+xmax(1)

         xminb(2:)=xmin(2:) - lxmax(2:)
         xmaxb(2:)=xmax(2:) + lxmax(2:)

         rdmax=lxmax(1)
       
      case(2)
         ewrite(3,*) 'using an xz plane'

         xminb(1)=xmin(1) - lxmax(1)
         xmaxb(1)=xmax(1) + lxmax(1)

         xminb(2)=-lxmax(2)+xmax(2)
         xmaxb(2)= lxmax(2)+xmax(2)

         xminb(3:)=xmin(3:) - lxmax(3:)
         xmaxb(3:)=xmax(3:) + lxmax(3:)

         rdmax=lxmax(2)

      case(3)
         ewrite(3,*) 'using an xz plane'

         xminb(1:2)=xmin(1:2) - lxmax(1:2)
         xmaxb(1:2)=xmax(1:2) + lxmax(1:2)

         xminb(3)=-lxmax(3)+xmax(3)
         xmaxb(3)= lxmax(3)+xmax(3)

         rdmax=lxmax(3)

      end select

      dxb = xmaxb - xminb

      ! generate initial coordinades, signs of turbulent spots
      ! calculate convection velocity
      if (.not.initz) then
         allocate(initeddy(sem_bc_count))
         initeddy = .true.
         
         allocate(averaged_velocity(dim, sem_bc_count))
         allocate(area(sem_bc_count))
         averaged_velocity=0.; area=0.
       
         initz=.true.
      endif

      if(initeddy(ns))then
         ! This is only done on process zero. In theory we could load balance this.
         if(GetRank()==0) then
            allocate(coords(dim,number_of_eddies))
            call random_seed()
            call random_number(coords)
            
            !place each eddy randomly inside its box
            
            do i=1,number_of_eddies
               eddies(ns)%eddy_positions(:,i)=xminb+(xmaxb-xminb)*coords(:,i)
            end do
            deallocate(coords)

            eddies(ns)%eddy_polarities=esign(dim,number_of_eddies)
         end if

#ifdef HAVE_MPI
         if(IsParallel())then
            call mpi_bcast(eddies(ns)%eddy_positions(:,:),number_of_eddies*dim,GetPREAL(),0,MPI_COMM_FEMTOOLS,ierr)

            call mpi_bcast(eddies(ns)%eddy_polarities(:,:),number_of_eddies*dim,mpi_integer,0,MPI_COMM_FEMTOOLS,ierr)
         endif
#endif

         x_shape   => ele_shape(bc_position,1)
         vel_shape => ele_shape(mean_velocity,1)
         
         total_area=0.
         do ele=1,element_count(mean_velocity)
            if (.not. element_owned(bc_position,ele) ) cycle
            call transform_to_physical(bc_position,ele,detwei=detwei)
            vel_int=shape_rhs(vel_shape, detwei)
            total_area=total_area+sum(vel_int)
            averaged_velocity(:,ns)=averaged_velocity(:,ns)+matmul(ele_val(mean_velocity,ele),vel_int)
         enddo

         call allsum(total_area)
         call allsum(averaged_velocity(:,ns))

         averaged_velocity(:,ns)=averaged_velocity(:,ns)/total_area

         area(ns)=total_area
         initeddy(ns)=.false.
      endif
      
      ewrite(3,*) 'surface info'
      ewrite(3,*) 'x min&max',xmin,xmax,dxp
      ewrite(3,*) 'surface area',area(ns)
      ewrite(3,*) 'convection velocity:',averaged_velocity(:,ns)

      ! bounding box volume
      volume=area(ns)*2*rdmax

      ewrite(3,*) 'bounding box info'
      do i=1,dim
         ewrite(3,*) 'i min&max',i, xminb(i),xmaxb(i)
      end do
      ewrite(3,*) 'box volume',volume
      
      if(GetRank()==0) then
         ! convect turbulent spots by the representative averaged velocity
         do i=1,number_of_eddies
            eddies(ns)%eddy_positions(:,i)=eddies(ns)%eddy_positions(:,i)+averaged_velocity(:,ns)*dt
         end do

         ! regenerate positions of eddies which have left the box
         do i=1,number_of_eddies
            do j =1, dim
               if(eddies(ns)%eddy_positions(j,i)>xmaxb(j))then
                  call random_number(tmpcord)
                  ratio= int((abs(eddies(ns)%eddy_positions(j,i)-xminb(j)))/dxb(j))

                  eddies(ns)%eddy_positions(:,i)=xminb+(xmaxb-xminb)*tmpcord
                  eddies(ns)%eddy_positions(j,i)=eddies(ns)%eddy_positions(j,i)-ratio*dxb(j)
                  eddies(ns)%eddy_polarities(:,i:i)=esign(dim,1)

               elseif(eddies(ns)%eddy_positions(j,i)<xminb(j))then
                  call random_number(tmpcord)

                  ratio= int((abs(eddies(ns)%eddy_positions(j,i)-xmaxb(j)))/dxb(j))

                  eddies(ns)%eddy_positions(:,i)=xminb+(xmaxb-xminb)*tmpcord
                  eddies(ns)%eddy_positions(j,i)=eddies(ns)%eddy_positions(j,i)+ratio*dxb(j)
                   
                  eddies(ns)%eddy_polarities(:,i:i)=esign(dim,1)

               end if
            end do
         end do
      end if
#ifdef HAVE_MPI
      if(IsParallel())then
         call mpi_bcast(eddies(ns)%eddy_positions(:,:),number_of_eddies*dim,GetPREAL(),0,MPI_COMM_FEMTOOLS,ierr)

         call mpi_bcast(eddies(ns)%eddy_polarities(:,:),number_of_eddies*dim,mpi_integer,0,MPI_COMM_FEMTOOLS,ierr)
      endif
#endif

      ! calculate fluctuating eddy components on inlet velocity
      
      call allocate(ufl,eddy_inlet_velocity%dim,eddy_inlet_velocity%mesh,"FluctuatingVelocities")
      call zero(ufl)

      do i=1,node_count(ufl)

         uf=0.0
         rts=node_val(turbulence_lengthscales,i)

         do j=1,number_of_eddies
            dx=abs(node_val(bc_position,i)-eddies(ns)%eddy_positions(:,j))

            if (all(dx<rts)) then

               ff=product(1.-dx/rts)/product(sqrt(3./2.*rts))

               uf=uf+eddies(ns)%eddy_polarities(:,j)*ff
             end if
          end do

          uf=uf*sqrt(volume/number_of_eddies)

          res=node_val(reynolds_stresses,i)

          if (any(res<0.)) then
             ewrite(3,*) 'WARNING: there is a negative value in Reynolds stress matrix'
          end if

          where(res<0.)
             res =1.e-5
          end where


          ! A Cholesky decomposition of a Reynolds stress matrix assumed to
          ! be diagonal.

          ! The full form is:

          !       [ \sqrt{R_11}          0                         0              ]
          ! A  =  [  R_21/A_11  \sqrt{R_22-A_21^2}                 0              ]
          !       [  R_31/A_11  (R_32-A_21A_31/A_22 \sqrt{R_33 - A_31^2 - A_32^2} ]

          a=0.0
          do j=1,dim
             a(j*(j+1)/2)=sqrt(res(j))
          end do

          do j=1,dim
             call set(ufl,j,i,dot_product(a((j-1)*j/2+1:j*(j+1)/2),uf(:j)))
          end do

       end do

       min_mean_velocity=minval(mean_velocity%val,dim=2)
       max_mean_velocity=maxval(mean_velocity%val,dim=2)
            
       do i=1,dim
          call allmin(min_mean_velocity(i))
          call allmax(max_mean_velocity(i))
       end do

       call set(eddy_inlet_velocity,mean_velocity)
         
       call addto(eddy_inlet_velocity,ufl)

       re_ii_min = minval(reynolds_stresses%val,dim=2)
       re_ii_max = maxval(reynolds_stresses%val,dim=2)         
         
         do i=1,dim
            call allmin(re_ii_min(i))
            call allmax(re_ii_max(i))
         end do

         uflmin = minval(ufl%val,dim=2)
         uflmax = maxval(ufl%val,dim=2)

         do i=1,dim
            call allmin(uflmin(i))
            call allmax(uflmax(i))
         end do

         ewrite(3,*) 'inlet condition diagnostics'
         ewrite(3,*) 'Re stresses min&max'
      do i=1,dim
         ewrite(3,*) 'Re_ii',i,':',re_ii_min(i),':',re_ii_max(i)
      end do
      ewrite(3,*) 'mean:fluctuating min&max'
      do i=1,dim
         ewrite(3,*) i, ':', min_mean_velocity(i), max_mean_velocity(i),':',uflmin(i),uflmax(i)
      end do

      call deallocate(ufl)

      ewrite(3,*) 'leaving turbulent inlet boundary conditions routine'
    end subroutine synthetic_eddy_method

    !----------------------------------------------------------

    function esign(dim, num_eddies)

      integer, intent(in) :: dim, num_eddies

      integer :: esign(dim, num_eddies)
      real    :: i(dim, num_eddies)

      ! Assign polarity to the synthetic eddies on a given boundary. 

      call random_number(i)

      where(i<0.5)
        esign = -1        
      elsewhere
        esign = +1
      end where

      return 
     end function esign

     !----------------------------------------------------------
     
     subroutine initialise_sem_memory(ns,dim,number_of_eddies)

       logical,save::initialise_memory=.false.
       logical,allocatable,dimension(:),save::initeddymem
       integer:: ns, dim, number_of_eddies

       if(.not.initialise_memory)then
          allocate(initeddymem(sem_bc_count))
          initeddymem=.true.
          allocate(eddies(sem_bc_count))
          initialise_memory=.true.
       endif
       
       if(initeddymem(ns))then
          allocate(eddies(ns)%eddy_positions(dim, number_of_eddies))
          allocate(eddies(ns)%eddy_polarities(dim, number_of_eddies))  
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
   
