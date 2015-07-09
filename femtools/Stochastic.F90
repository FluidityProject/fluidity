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
#include "version.h"

module stochastic

  use fldebug
  use fields
  use parallel_tools

  implicit none

  integer, parameter :: UNIFORM=1, NORMAL=2
  logical, save :: initialised=.false. 

  private

  public :: get_random_points_in_mesh, initialise_stochastic_module,&
       partition_points_in_parallel

contains

  subroutine init_random_seed()
    ! Procedure is taken from the examples in the GCC online document
    use iso_fortran_env, only: int64
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid
    integer(int64) :: t
    
    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(t)
       if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
       end if
       pid = getpid()
       t = ieor(t, int(pid, kind(t)))
       do i = 1, n
          seed(i) = lcg(t)
       end do
    end if
    call random_seed(put=seed)
  contains
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    function lcg(s)
      integer :: lcg
      integer(int64) :: s
      if (s == 0) then
         s = 104729
      else
         s = mod(s, 4294967296_int64)
      end if
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
    end function lcg
  end subroutine init_random_seed

  subroutine initialise_stochastic_module()
    if (initialised) return
    call init_random_seed()
    initialised=.true.
  end subroutine initialise_stochastic_module

  subroutine get_random_from_distribution(harvest,distribution)
    real, intent(out), dimension(:) :: harvest
    integer, intent(in) :: distribution
    select case(distribution)
    case(UNIFORM)
       call uniform_random_number(harvest)
    case(NORMAL)
       call normal_random_number(harvest)
    end select
  end subroutine get_random_from_distribution

  subroutine uniform_random_number(harvest)
    real, dimension(:), intent(out) :: harvest
    
    call random_number(harvest)
  end subroutine uniform_random_number

  subroutine normal_random_number(harvest)

! basic Box-Muller algorithm for an arbitrary array.
! algorithm produces two variables from the correct distribution each time
! and we choose not to waste them.

    real, dimension(:), intent(out) :: harvest
    
    real :: x1(1),x2(1),w
    integer :: i

    do i=1,size(harvest)/2
       w=1.0
       do while (w>=0)
          call uniform_random_number(x1)
          call uniform_random_number(x2)
          x1(1)=2.0*x1(1)-1.0
          x2(1)=2.0*x2(1)-1.0
          w=x1(1)*x1(1)+x2(1)*x2(1)
       end do
       w=sqrt(-2.0*log(w)/w)
       harvest(2*i-1)=x1(1)*w
       harvest(2*i)=x2(1)*w
    end do
    
    if (mod(size(harvest),2)==1) then
       w=1.0
       do while (w>=0)
          call uniform_random_number(x1)
          call uniform_random_number(x2)
          x1(1)=2.0*x1(1)-1.0
          x2(1)=2.0*x2(1)-1.0
          w=x1(1)*x1(1)+x2(1)*x2(1)
       end do
       w=sqrt(-2.0*log(w)/w)
       harvest(size(harvest))=x1(1)*w
    end if

  end subroutine normal_random_number



    


  function get_random_point_in_simplex(X,p) result(point)
    real, dimension(:,:) :: X
    real, dimension(:), optional :: p
    real, dimension(size(X,1)) :: point

    ! locals
    real, dimension(size(X,2)-1) :: q
    real, dimension(size(X,2)) :: q_fixed
    integer ::nloc, i

    nloc=size(q_fixed)
    call uniform_random_number(q)

    if (present(p)) then
       FLAbort("I haven't finished this code yet")
    else
       do i=1,nloc-1
          q_fixed(i)=(1.0-sum(q_fixed(1:i-1)))*(1-q(i)**(1.0/(nloc-i)))
       end do
       q_fixed(nloc)=1.0-sum(q_fixed(1:nloc-1))
    end if

    point=matmul(X,q_fixed)

  end function get_random_point_in_simplex

  function get_random_points_in_mesh(X,npoints,pdf) result(pts)
    type(vector_field), intent(in) :: X
    integer, intent(in)            :: npoints
    type(scalar_field), intent(in), optional :: pdf

    real, dimension(X%dim,npoints) :: pts
    ! locals

    integer :: ele, i
    real, dimension(element_count(X)) :: cumulative_volume
    logical, dimension(npoints) :: mask
    real, dimension(npoints)    :: q
    integer, dimension(npoints)    :: elist

    !! if there's nothing to do, don't do it.
    if (npoints==0) return

    if (element_count(X)==0) then
       FLAbort("Trying to allocate stochastic detectors on a domain of zero measure. Check any surface ids are actually present in the mesh")
    end if

    if (present(pdf)) then
       FLAbort('I still need to write this.')
    else
       if (element_owned(X,1)) then
          cumulative_volume(1)=element_volume(X,1)
       else
          cumulative_volume(1)=0.0
       end if
       do ele=2,element_count(X)
          if (element_owned(X,ele)) then
             cumulative_volume(ele)=element_volume(x,ele)&
                  +cumulative_volume(ele-1)
          else
             cumulative_volume(ele)=cumulative_volume(ele-1)
          end if
       end do

       cumulative_volume=cumulative_volume&
            /cumulative_volume(size(cumulative_volume))
    end if

    call uniform_random_number(q)
    mask=.true.

    do ele=1,element_count(X)
       if( .not. element_owned(X,ele)) cycle
       where(mask .and. q<=cumulative_volume(ele))
          mask=.false.
          elist=ele
       end where
    end do

    do i=1,npoints
       if (present(pdf)) then
          pts(:,i)=get_random_point_in_simplex(ele_val(X,elist(i)),&
               ele_val(pdf,elist(i)))
       else
          pts(:,i)=get_random_point_in_simplex(ele_val(X,elist(i)))
       end if
    end do        
  end function get_random_points_in_mesh

  subroutine partition_points_in_parallel(X,total_points,npoints,offset,pdf)
    type(vector_field), intent(in) :: X
    integer, intent(in) :: total_points
    type(scalar_field), intent(in), optional :: pdf
    
    integer, intent(out):: npoints, offset

    ! local
    integer :: nprocs, procno
    integer, dimension(:), allocatable :: p_points
    real, dimension(:), allocatable :: vols
    real, dimension(:), allocatable :: q

    integer :: ele, i


    if (isparallel()) then
       nprocs=getnprocs()
       procno=getprocno()
       if (present(pdf)) then
          FLAbort('To be written')
       else
          npoints=total_points/nprocs
          if(procno==nprocs)&
             npoints=total_points-(nprocs-1)*npoints

          allocate(vols(nprocs),p_points(nprocs),q(npoints))

          vols=0.0
          do ele=1,element_count(X)
             if (element_owned(X,ele)) &
                  vols(procno)=vols(procno)+element_volume(X,ele)
          end do
          
          vols(procno+1:nprocs)=vols(procno)
          call allsum(vols)
          vols=vols/vols(nprocs)

          call uniform_random_number(q)
          
          do i=1,nprocs
             p_points(i)=count(q<vols(i))
             where(q<vols(i))
                q=2.0
             end where
          end do

          call allsum(p_points)
          npoints=p_points(procno)
          offset=sum(p_points(1:procno-1))
       end if
    else
       npoints=total_points
       offset=0
    end if
  end subroutine partition_points_in_parallel
       
end module stochastic
