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
  use data_structures

  use iso_c_binding
  use global_parameters, only: OPTION_PATH_LEN
  use spud
  use parallel_tools
  use mpi_interfaces
  
  implicit none

  integer, parameter :: UNIFORM=1, NORMAL=2
  logical, save :: initialised=.false.
  integer, allocatable, dimension(:) :: generic_seed
  integer :: seed_size

  type(string_hash_table) :: named_seed

  interface get_random
     module procedure get_random_1d, get_random_2d
  end interface get_random

  interface get_random_on_root
     module procedure get_random_on_root_1d, get_random_on_root_2d
  end interface get_random_on_root
     

  private

  public :: get_random_points_in_mesh, initialise_stochastic_module,&
       partition_points_in_parallel, deallocate_seeds, checkpoint_seeds,&
       get_random, get_random_on_root, finalise_stochastic_module, UNIFORM, NORMAL

contains

  subroutine init_random_seed()
    ! Procedure is based heavily on the examples in the GCC online document
    use iso_fortran_env, only: int64
    implicit none
    integer, dimension(:), allocatable :: seed
    integer :: i, random_unit, istat, dt(8), pid
    integer(int64) :: t
    
    random_unit=free_unit()
    
    allocate(seed(seed_size))
    ! First check if the OS provides a random number generator
    open(newunit=random_unit, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(random_unit) seed
       close(random_unit)
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
       do i = 1, seed_size
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

  subroutine load_seed(key)
    character(len=*) :: key
    integer, dimension(:), pointer :: seed

    call random_seed(get=generic_seed)
    
    if (has_key(named_seed,key)) then
       call random_seed(put=fetch(named_seed,key,[seed_size]))
    else
       allocate(seed(seed_size))
       call init_random_seed()
       call random_seed(get=seed)
       call insert(named_seed,key,c_loc(seed))
    end if

  end subroutine load_seed

  subroutine store_seed(key)
    character(len=*) :: key
    integer, dimension(:), pointer :: seed

    !!! store the named series value
    seed=>fetch(named_seed,key,[seed_size])
    call random_seed(get=seed)
    
    !!! get the generic series
    call random_seed(put=generic_seed)
  end subroutine store_seed

  subroutine initialize_from_options_file()
    integer :: i
    character(len = OPTION_PATH_LEN) :: series_name
    integer, dimension(:), pointer :: seed

    !!! read seeds out of the options file


    allocate(seed(seed_size))
    do i=1, option_count('/embedded_models/stochastic_routines/seed')
       call get_option('/embedded_models/stochastic_routines/seed['//int2str(i-1)&
            //']/name', series_name)
       call get_option('/embedded_models/stochastic_routines/seed['//int2str(i-1)&
            //']', seed)
       if (trim(series_name) == 'Generic') then
          generic_seed=seed
       else
          print*, seed
          call insert(named_seed,trim(series_name),c_loc(seed))
          nullify(seed)
       end if
    end do

  end subroutine initialize_from_options_file

  subroutine checkpoint_seeds()

    integer :: stat 
    type(c_ptr) :: ptr
    character(len=254) :: key

    integer, dimension(:), pointer :: fptr

    call random_seed(get=generic_seed)
    call add_option('/embedded_models/stochastic_routines/seed::Generic',stat)
    call set_option('/embedded_models/stochastic_routines/seed::Generic',&
         generic_seed,stat)

    !!! walk the seeds
    call get_first(named_seed,key,ptr,stat)
    do while (stat>0)
       call c_f_pointer(ptr,fptr,[seed_size])
       call add_option('/embedded_models/stochastic_routines/seed::'&
            //trim(key),stat)
       call set_option('/embedded_models/stochastic_routines/seed::'&
            //trim(key),fptr,stat)
       call get_next(named_seed,key,ptr,stat)
    end do

  end subroutine checkpoint_seeds
       
  subroutine deallocate_seeds()
    integer :: stat
    type(c_ptr) :: ptr
    character(len=254) :: key

    integer, dimension(:), pointer :: fptr

    !!! walk the seeds

    call get_first(named_seed,key,ptr,stat)
    do while (stat>0)
       call c_f_pointer(ptr,fptr,[seed_size])
       deallocate(fptr)
       call get_next(named_seed,key,ptr,stat)
    end do

    deallocate(generic_seed)

  end subroutine deallocate_seeds
  
  subroutine initialise_stochastic_module()
    !! if already initialised, do nothing
    if (initialised) return
    !! else set up the seed(s)
    call random_seed(size = seed_size)
    allocate(generic_seed(seed_size))
    if (have_option('/embedded_models/stochastic_routines')) then
       call initialize_from_options_file()
    else
       call init_random_seed()
    end if
    initialised=.true.
  end subroutine initialise_stochastic_module
  
  subroutine finalise_stochastic_module()
    call deallocate_seeds()
    initialised=.false.
  end subroutine finalise_stochastic_module

  subroutine get_random_on_root_1d(harvest,distribution,series,root,comm)
    real, intent(out), dimension(:) :: harvest
    integer, intent(in) :: distribution
    character(len=*), optional :: series
    integer, optional :: root, comm

    integer :: lroot, lcomm
    integer :: ierr

    !!! this routine fills harvest with random numbers from the chosen distribution on the root node (default 1), then broadcasts it across the communicator comm (default MPI_COMM_FEMTOOLS, i.e. the entire system.)

    if (present(root)) then
       lroot = root
    else
       lroot = 1
    end if

    if (present(comm)) then
       lcomm=comm
    else
       lcomm = MPI_COMM_FEMTOOLS
    end if
      
    if (getprocno()==lroot) then
       call get_random(harvest,distribution,series)
    end if
 
    if (isparallel()) then
       call mpi_bcast(harvest,size(harvest),getpreal(),lroot-1,lcomm,ierr)
    end if

  end subroutine get_random_on_root_1d

  subroutine get_random_on_root_2d(harvest,distribution,series,root,comm)
    real, intent(out), dimension(:,:) :: harvest
    integer, intent(in) :: distribution
    character(len=*), optional :: series
    integer, optional :: root, comm

    integer :: lroot, lcomm
    integer :: ierr

    !!! this routine fills harvest with random numbers from the chosen distribution on the root node (default 1), then broadcasts it across the communicator comm (default MPI_COMM_FEMTOOLS, i.e. the entire system.)

    if (present(root)) then
       lroot = root
    else
       lroot = 1
    end if

    if (present(comm)) then
       lcomm=comm
    else
       lcomm = MPI_COMM_FEMTOOLS
    end if
      
    if (getprocno()==lroot) then
       call get_random(harvest,distribution,series)
    end if
 
    if (isparallel()) then
       call mpi_bcast(harvest,size(harvest),getpreal(),lroot-1,lcomm,ierr)
    end if

  end subroutine get_random_on_root_2d

  subroutine get_random_1d(harvest,distribution,series)
    real, intent(out), dimension(:) :: harvest
    integer, intent(in) :: distribution
    character(len=*), optional :: series

    if (present(series)) call load_seed(series)

    select case(distribution)
    case(UNIFORM)
       call uniform_random_number(harvest,size(harvest))
    case(NORMAL)
       call normal_random_number(harvest,size(harvest))
    end select

    if (present(series)) call store_seed(series)

  end subroutine get_random_1d

  subroutine get_random_2d(harvest,distribution,series)
    real, intent(out), dimension(:,:) :: harvest
    integer, intent(in) :: distribution
    character(len=*), optional :: series

    if (present(series)) call load_seed(series)

    select case(distribution)
    case(UNIFORM)
       call uniform_random_number(harvest,size(harvest))
    case(NORMAL)
       call normal_random_number(harvest,size(harvest))
    end select

    if (present(series)) call store_seed(series)

  end subroutine get_random_2d

  subroutine uniform_random_number(harvest,n)
    integer :: n
    real, dimension(n), intent(out) :: harvest
    call random_number(harvest)
  end subroutine uniform_random_number

  subroutine normal_random_number(harvest,n)

! basic Box-Muller algorithm for an arbitrary array.
! algorithm produces two variables from the correct distribution each time
! and we choose not to waste them.

    integer :: n
    real, dimension(n), intent(out) :: harvest
    
    real :: x1(1),x2(1),w
    integer :: i

    do i=1,size(harvest)/2
       w=1.0
       do while (w>=0)
          call uniform_random_number(x1,size(x1))
          call uniform_random_number(x2,size(x2))
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
          call uniform_random_number(x1,size(x1))
          call uniform_random_number(x2,size(x2))
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
    call get_random(q,UNIFORM)

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

  function get_random_points_in_mesh(X,npoints) result(pts)
    type(vector_field), intent(in) :: X
    integer, intent(in)            :: npoints

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

    call get_random(q,UNIFORM)
    mask=.true.

    do ele=1,element_count(X)
       if( .not. element_owned(X,ele)) cycle
       where(mask .and. q<=cumulative_volume(ele))
          mask=.false.
          elist=ele
       end where
    end do

    do i=1,npoints
       pts(:,i)=get_random_point_in_simplex(ele_val(X,elist(i)))
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

          call get_random(q,UNIFORM)
          
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
