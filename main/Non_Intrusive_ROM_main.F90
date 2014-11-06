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
  module Nonintrusive_rom_module

  use spud
  use fields
  use state_module
  use write_state_module
  use timeloop_utilities
  use global_parameters, only: option_path_len, current_time, dt, FIELD_NAME_LEN ,&
                               timestep, OPTION_PATH_LEN, &
                               simulation_start_time, &
                               simulation_start_cpu_time, &
                               simulation_start_wall_time, &
                               topology_mesh_name,nd_total 
  use FLDebug
  use snapsvd_module
  use vtk_interfaces
  use memory_diagnostics
  use populate_state_module
  use python_state
  use vtk_interfaces
  use quadrature
  use diagnostic_fields_wrapper
  use field_options
  use signals
  use spud
  use tictoc
  use fldebug
  use fluids_module
  !use reduced_fluids_module
  use signals
  use spud
  use tictoc
!#ifdef HAVE_ZOLTAN
  use zoltan
  use nonlinear_conjugate_gradient
  use momentum_equation
  use momentum_equation_reduced_adjoint
  use momentum_equation_reduced
  use reduced_model_runtime
  use smolyak 
  use rbf_interp
  
!#endif
  implicit none

 
#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif 


  private


  public :: Non_Intrusive_ROM_main,Non_Intrusive_smolyak_main, nd_total 
 ! integer  nd_total
contains 

  subroutine Non_Intrusive_smolyak_main
   
   implicit none
!    We have functions that can generate a Lagrange interpolant to data
!    in M dimensions, with specified order or level in each dimension.
!
!    We use the Lagrange function as the inner evaluator for a sparse
!    grid procedure. 
!
!    The procedure computes sparse interpolants of levels 0 to SPARSE_MAX
!    to a given function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Local, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = 8 ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer ( kind = 4 ) C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer ( kind = 4 ) IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer ( kind = 4 ) L, the current Lagrange grid level.
!
!    Local, integer ( kind = 4 ) L_MAX, the current sparse grid level.
!
!    Local, integer ( kind = 4 ) MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ( kind = 4 ) ND, the number of points used in a Lagrange grid.
!
!    Local, integer ( kind = 4 ) ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer ( kind = 4 ) NI, the number of interpolant evaluation points.
!
!    Local, integer ( kind = 4 ) SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = 8 ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = 8 ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = 8 ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = 8 ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = 8 ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = 8 ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
 
    type(state_type), dimension(:), pointer :: state=> null()
    ! type(state_type), dimension(:), allocatable :: state_adj
    type(state_type), dimension(:,:,:), allocatable :: POD_state
    
    integer :: ierr,stat,i,j,k
    character(len = OPTION_PATH_LEN) :: simulation_name
    integer :: nsvd
    integer :: total_timestep,total_dump_no
    type(vector_field), pointer :: snapmean_velocity
    type(scalar_field), pointer :: snapmean_pressure
    type(vector_field), pointer :: u ,POD_velocity
    type(scalar_field), pointer :: p,POD_pressure
    real, dimension(:), allocatable :: pod_coef_optimised
    logical :: if_optimal,if_timetwo
    !  type(state_type), dimension(:,:,:) :: POD_state
    ! Change in pressure
    type(scalar_field) :: delta_p
    ! Change in velocity
    type(vector_field) :: delta_u
    real, dimension(:,:), allocatable :: pod_sol_velocity
    real, dimension(:), allocatable :: pod_sol_pressure

    ! non_intrusive_ROM
    type(state_type), dimension(:,:,:), allocatable :: Non_Intrusive_state
    real, dimension(:), allocatable ::  pod_coef_orig,pod_coef_new,pod_coef_old
    real, dimension(:,:), allocatable ::  DM0,DM1
    real, dimension(:,:,:), allocatable ::  DDM
    integer istate,NPERT,ntime
    integer dump_no_tmp, int_dump_period
    real current_time,finish_time,dt
    !smolyak
    real:: cpu_start1,cpu_end1,cpu_start2,cpu_end2

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error,app_error1
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  real ( kind = 8 ), allocatable :: fi(:) 
  logical more
  integer ( kind = 4 ) nd,nd_lq
  !integer ( kind = 4 ) nd_total
  integer ( kind = 4 ) ni
  ! real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ) r0
  real ( kind = 8 ), allocatable :: xd(:,:),xd_lq(:,:),pod_coef_all_obv(:,:),smoylak_all(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:),fd(:)
  real ( kind = 8 ), allocatable :: ze(:),temp_coef(:)
  real ( kind = 8 ), allocatable :: zi(:),minin(:),maxin(:)
  real ( kind = 8 ), allocatable :: zpi(:)
  real ( kind = 8 ), allocatable :: fe(:,:) !exact values of function produced from fluidity
  integer :: m,ll,kk 
  !real ( kind = 8 ) w_lq(total_timestep)
  !real ( kind = 8 ) fd(nd_lq)
  integer :: n1d,podnum
  real :: mean, positive,negtive,num_pos,num_neg
  ! real minin(podnum),maxin(podnum) 
  call get_option(&
            '/reduced_model/pod_basis_formation/pod_basis_count', nsvd) 
  call delete_option("/reduced_model/execute_reduced_model") ! switch to full model 
  call get_option("/timestepping/current_time", current_time)
  call get_option("/timestepping/finish_time", finish_time)       
  call get_option("/timestepping/timestep", dt)
  ntime=int((finish_time-current_time)/dt)
  total_timestep=int((finish_time-current_time)/dt)
  m = 3*nsvd
  podnum = 3*nsvd
  nd= total_timestep  
  nd_lq=200
  n1d=3
   call get_option(&
            '/reduced_model/smolyak_order', sparse_max) 
  !sparse_max=1
  allocate(pod_coef_all_obv(total_timestep,podnum))
  allocate(fd(nd_lq))
  allocate(temp_coef(m))
  allocate(minin(podnum))
  allocate(maxin(podnum))
  !allocate(pod_coef_all(total_timestep,podnum))
  print *, 'dd'
  open(unit=61,file='coef_pod_all_obv')
  read(61,*)((pod_coef_all_obv(j,k),k=1,podnum),j=1,total_timestep)
  close(61)
    if(have_option("/io/dump_period_in_timesteps/constant")) then
       call get_option("/io/dump_period_in_timesteps/constant", int_dump_period)
    else
       int_dump_period = 1
    endif
   
   Do i=1,m   ! cannot use minval because some of the vector is null.          
          minin(i)=pod_coef_all_obv(1,i)
          maxin(i)=pod_coef_all_obv(1,i)      
     do j=1,total_timestep
           if(minin(i)>pod_coef_all_obv(j,i)) then
              minin(i)=pod_coef_all_obv(j,i)          
           endif
           if(maxin(i)<pod_coef_all_obv(j,i)) then
              maxin(i)=pod_coef_all_obv(j,i)
           endif
      enddo
   ENDDO 
   
  positive=0
  negtive=0
  mean=0
  num_pos=0
  num_neg=0
  do j=1,m
    do k=1,nd
      if(pod_coef_all_obv(k,j)>0) then 
         positive=positive+pod_coef_all_obv(k,j)
         num_pos=num_pos+1
      else 
         negtive=negtive+pod_coef_all_obv(k,j)
         num_neg=num_neg+1
      endif
      mean=mean+pod_coef_all_obv(k,j)
    enddo
  enddo
    mean=mean/(num_pos+num_neg)
    positive=positive/num_pos
    negtive=negtive/num_neg
 ! do j=1,m
  ! do k=1,nd
     ! xd_lq(j,k)=pod_coef_all_obv(k,j) 
  ! enddo
 ! enddo
  
  m=podnum
  allocate ( ind(1:m) )  ! level of each dimension, so has the size of m
  
  !  Define the region.
  !
  allocate ( a(1:m) )
  allocate ( b(1:m) )

  a(1:m) =  minin(1:m)!!0.0D+00negtive
  b(1:m) =  maxin(1:m)!!2.0D+00positive

  ! Define the interpolation evaluation information.
  ! open(40,file='total', ACTION='read') 
  ! read(40,*) total_nd
  ! close(40)
   allocate (fe(1:m,1:total_nd)) 
   ni = 1
   allocate ( xi(1:m,1:ni) )
   allocate ( fi(1:ni) ) 
   
   allocate ( ze(1:ni) )
   ! call f_sinr ( m, ni, xi, ze ) 
   ! Compute a sequence of sparse grid interpolants of increasing level.
  
   allocate ( zpi(1:ni) )
   allocate ( zi(1:ni) )
    !----------------------------------
    ! get State
    !----------------------------------
    call initialise_write_state
    ! Read state from .flml file
    call populate_state(state)
    u => extract_vector_field(state(1), "Velocity", stat)
    p => extract_scalar_field(state(1), "Pressure") 
    call read_pod_basis_differntmesh(POD_state, state)
    allocate(pod_coef_orig((u%dim+1)*size(POD_state,1)))
    allocate(pod_coef_new((u%dim+1)*size(POD_state,1)))
    allocate(pod_coef_old((u%dim+1)*size(POD_state,1))) 
     ! Allocate the change in pressure field
    call allocate(delta_p, p%mesh, "DeltaP")
    delta_p%option_path = trim(p%option_path)
    call zero(delta_p)
    
    ! Allocate the change in velocity pressure field
    call allocate(delta_u, u%dim, u%mesh, "DeltaU")
    delta_u%option_path = trim(u%option_path)
    call zero(delta_u)
    !----------------------------------
    ! get the initial ROM solution coef.
    !----------------------------------
    istate = 1
    call project_from_full_to_pod(istate, pod_state, state, pod_coef_orig) 
      xi(:,1)=pod_coef_orig(:) 
      !print *, 'ddd',  xi(:,1)
    do kk =1, ntime,int_dump_period
       xi(:,1)=pod_coef_all_obv(kk,:) 
     do k=1, podnum
         call cpu_time(cpu_start1)
       sparse_min = 0
       lagsum=0
       nd_total = 0
       if(.not.(have_option("/reduced_model/Produce_Smolyak_Data"))) then
           open(44,file='pod_coef_smolyak')
        endif
     do l_max = sparse_min, sparse_max
      !print *, 'dddd'
      allocate ( c(0:l_max) )
      allocate ( w(0:l_max) )
      call smolyak_coefficients ( l_max, m, c, w )  
     !    Output, integer ( kind = 4 ) C(0:L_MAX), the coefficients for objects 
     !    at sublevels 0 through L_MAX.
     !    Output, integer ( kind = 4 ) W(0:L_MAX), the number of objects at 
     !    sublevels 0 through L_MAX.

    zi(1:ni) = 0.0D+00
    !nd_total = 0

    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.

      do 
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
       
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )
        allocate ( zd(1:nd))
        do i=1, nd
        ! if(nd_total .eq. 0) then
            open(10,file='xd', ACTION='WRITE') 
         !  else
         !   open(10,file='xd', position='append',ACTION='WRITE')
         ! endif 
          write(10,*)  xd(:,i)
          close(10)
          !call f_sinr ( m, nd, xd, zd )
          !get the values of interpolation points 
          !zd(1:nd)=fe(k,(nd_total+1):(nd_total + nd))
           if( (have_option("/reduced_model/Produce_Smolyak_Data"))) then
             call fluids()                   !open when produce data
           else 
             read(44,*) temp_coef(1:m) 
            endif
          nd_total = nd_total + 1
          lagsum=lagsum+1  
        ! print *, 'before call fluids and k^th podcoef', nd_total, k   
          zd(i)=temp_coef(k)
        enddo !do i=1, nd
     
        !nd_total = nd_total + nd
        !  Use the grid to evaluate the interpolant. 
        call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi ) 
        !nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)

        deallocate ( xd )
        deallocate ( zd )

        if ( .not. more ) then
          exit
        end if
      end do
    end do
    pod_coef_new(k)=zi(1)  
    deallocate ( c )
    deallocate ( w )

  end do  ! do l_max = sparse_min, sparse_max
   if(have_option("/reduced_model/Produce_Smolyak_Data")) then 
      print *,'finished producing data' 
      stop 2222 
       else
       close(44)
    endif
   call cpu_time(cpu_end1)
   print *,'cpuint',cpu_end1-cpu_start1
  enddo  !do k=1, podnum 
       !print *, 'pod_coef_new,timestep', kk!pod_coef_new!,! !size(pod_coef_new) pod_coef_new=pod_coef_all_obv(kk,:)!size(pod_coef_new)
       call cpu_time(cpu_start2) 
       call project_full(delta_u, delta_p, pod_sol_velocity, pod_sol_pressure, POD_state(:,:,istate), pod_coef_new)
       call cpu_time(cpu_end2)
        print *,'cpuprojection',cpu_end1-cpu_start1
         !stop 555
       snapmean_velocity=>extract_vector_field(POD_state(1,1,istate),"SnapmeanVelocity")
       snapmean_pressure=>extract_Scalar_field(POD_state(1,2,istate),"SnapmeanPressure")
       !call addto(u, delta_u, dt)
       u%val=snapmean_velocity%val
       call addto(u, delta_u)
       p%val=snapmean_pressure%val
       call addto(p, delta_p) 
       dump_no_tmp= int(kk/int_dump_period)       
       call write_state(dump_no_tmp, state)
        !call calculate_diagnostic_variables(state)
      ! call calculate_diagnostic_variables_new(state)          
       ! Call the modern and significantly less satanic version of study
      !  call write_diagnostics(state, real(kk), dt, timestep)
   
    enddo !kk,timestep loop
    !close(20)
    !----------------------------------
    ! deallocate variables
    !----------------------------------
   if (allocated(pod_state)) then
       do i=1, size(pod_state,1)
          do j=1,size(pod_state,2)
             do k=1,size(pod_state,3)
                call deallocate(pod_state(i,j,k))
             enddo
          enddo
       end do
    end if

    deallocate(pod_coef_orig)
    deallocate(pod_coef_new)
    deallocate(pod_coef_old)
    deallocate(maxin)
    deallocate(minin)
  end subroutine Non_Intrusive_smolyak_main


 subroutine test03 ( m, sparse_max )

!*****************************************************************************80
!
!! TEST01: sequence of sparse interpolants to an M-dimensional function.
!
!  Discussion:
!
!    We have functions that can generate a Lagrange interpolant to data
!    in M dimensions, with specified order or level in each dimension.
!
!    We use the Lagrange function as the inner evaluator for a sparse
!    grid procedure. 
!
!    The procedure computes sparse interpolants of levels 0 to SPARSE_MAX
!    to a given function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Local, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = 8 ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer ( kind = 4 ) C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer ( kind = 4 ) IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer ( kind = 4 ) L, the current Lagrange grid level.
!
!    Local, integer ( kind = 4 ) L_MAX, the current sparse grid level.
!
!    Local, integer ( kind = 4 ) MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ( kind = 4 ) ND, the number of points used in a Lagrange grid.
!
!    Local, integer ( kind = 4 ) ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer ( kind = 4 ) NI, the number of interpolant evaluation points.
!
!    Local, integer ( kind = 4 ) SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = 8 ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = 8 ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = 8 ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = 8 ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = 8 ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = 8 ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ) app_error
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) h
  integer ( kind = 4 ), allocatable :: ind(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_max
  integer ( kind = 4 ) l_min
  integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nd_total
  integer ( kind = 4 ) ni
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) sparse_max
  integer ( kind = 4 ) sparse_min
  integer ( kind = 4 ) t
  integer ( kind = 4 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: xd(:,:)
  real ( kind = 8 ), allocatable :: xi(:,:)
  real ( kind = 8 ), allocatable :: zd(:)
  real ( kind = 8 ), allocatable :: ze(:)
  real ( kind = 8 ), allocatable :: zi(:)
  real ( kind = 8 ), allocatable :: zpi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Sparse interpolation for a function f(x) of M-dimensional argument.'
  write ( *, '(a)' ) '  Use a sequence of sparse grids of levels 0 through SPARSE_MAX.'
  write ( *, '(a)' ) '  Invoke a general Lagrange interpolant function to do this.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Compare the exact function and the interpolants at a grid of points.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The "order" is the sum of the orders of all the product grids'
  write ( *, '(a)' ) '  used to make a particular sparse grid.'
!
!  User input.
!
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Spatial dimension M = ', m
  write ( *, '(a,i4)' ) '  Maximum sparse grid level = ', sparse_max

  allocate ( ind(1:m) )
!
!  Define the region.
!
  allocate ( a(1:m) )
  allocate ( b(1:m) )

  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00
!
!  Define the interpolation evaluation information.
!
  ni = 1
  seed = 123456789
  allocate ( xi(m,ni) )
  call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )

  write ( *, '(a,i6)' ) '  Number of interpolation points is NI = ', ni

  allocate ( ze(1:ni) )
  call f_sinr ( m, ni, xi, ze )
  print *, 'actual function value is, will interpolate this value later on,then do comparison', xi,ze
!
!  Compute a sequence of sparse grid interpolants of increasing level.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   L     Order    ApproxError'
  write ( *, '(a)' ) ''

  allocate ( zpi(1:ni) )
  allocate ( zi(1:ni) )

  sparse_min = 0

  do l_max = sparse_min, sparse_max

    allocate ( c(0:l_max) )
    allocate ( w(0:l_max) )
    call smolyak_coefficients ( l_max, m, c, w )

    zi(1:ni) = 0.0D+00
    nd_total = 0

    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.

      do
!
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )
        allocate ( zd(1:nd) )
        call f_sinr ( m, nd, xd, zd )
        print *, 'gridxd(1,1)', zd
!
!  Use the grid to evaluate the interpolant.
!
        call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi )
!
!  Weighted the interpolant values and add to the sparse grid interpolant.
!
        nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)

        deallocate ( xd )
        deallocate ( zd )

        if ( .not. more ) then
          exit
        end if

      end do

    end do
!
!  Compare sparse interpolant and exact function at interpolation points.
!
   ! app_error = r8vec_norm_affine ( ni, zi, ze ) / real ( ni, kind = 8 )

    !write ( *, '(2x,i2,2x,i8,2x,e8.2)' ) l, nd_total, app_error

    deallocate ( c )
    deallocate ( w )

  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )
  deallocate ( zpi )

  return
end


  SUBROUTINE Nonintrusive(NPERT,Ntime,DM0,DDM,TOLDKEEP_ORIG,told,tnew)
    IMPLICIT NONE
    INTEGER NPERT,NTIME
    REAL TNEW(NPERT),TOLD(NPERT),TOLDold(NPERT), TOLDkeep(NPERT)
    REAL TOLDKEEP_ORIG(NPERT)
    real, dimension(:,:), allocatable ::  M,MTM,MTM_INV,MTM_INV_M,M2,M1,M0,m00
    REAL DF(NPERT),DdF(NPERT),MTM_INV_F(NPERT)
    INTEGER ITIME,I,IPERT,JPERT, ISTART, IEND,kpert,k,ivar
    REAL TOLER
    real, dimension(NPERT,NPERT)::  DM0
    real, dimension(:,:), allocatable ::  DM1,DM2
    real, dimension(NPERT,NPERT,NPERT) ::  DDM
    real toler1,toler2
    real j
    
    allocate(DM2(NPERT,NPERT))
    
    ! ****the first two time steps...
    
    toldold=told
    told=tnew
    open(20,file='rom_solution_alltime.dat')
    write(20,*)(TOLDOLD(i),i=1,NPERT)
    write(20,*)(TNEW(i),i=1,NPERT)

    ! ****Perform the time stepping next...
    DO Itime=2,ntime
       tnew=told            
       do i=1,NPERT
          tnew(i)=tnew(i) + sum( DM0(i,:)*(told(:)-toldold(:)) )
       enddo
       
       do kpert=1,NPERT
          do i=1,NPERT
             tnew(i)=tnew(i) +(toldold(kpert)-TOLDKEEP_ORIG(kpert)  +   0.5*(told(kpert)-toldold(kpert)))  &
                  *sum( DDM(i,:,kpert)*(told(:)-toldold(:)) )
          end do
       end do
       
       toldold=told
       TOLD=tnew  

       write(20,*)(tnew(i),i=1,NPERT)

    END DO
    
    close(20)
    deallocate(DM2)
    
  END SUBROUTINE NONINTRUSIVE

  SUBROUTINE First_Order_DM(state,pod_state,TOLDKEEP_ORIG,DM0,if_timetwo) 
    type(state_type), dimension(:), intent(in) :: state 
    type(state_type), dimension(:,:,:), intent(in) :: pod_state
    INTEGER NONODS,NTIME,NPERT
    REAL,dimension(:), intent(inout) ::TOLDKEEP_ORIG
    real toler
    PARAMETER(TOLER=0.00001)
    REAL,dimension(:,:), intent(inout) :: DM0
    logical, optional :: if_timetwo
    
    ! local variables...
    integer ipert,istate,i,j,k
    logical if_optimal
    real, DIMENSION( : ), allocatable :: told,tnew
    real, DIMENSION( :, : ), allocatable :: m0,m1

    NPERT=size(TOLDKEEP_ORIG)
    allocate(told(NPERT))
    allocate(tnew(NPERT))
    allocate(M0(NPERT,NPERT))
    allocate(M1(NPERT,NPERT))
    

    TOLD=TOLDKEEP_ORIG ! no perturbations - a bit wasteful of CPU
    open(10,file='nonintrusive_coef_pod0.dat')
    write(10,*) (TOLD(i),i=1,NPERT)
    close(10)
    !------run fluids------
    call fluids()
    !------end run fluids------
    
    open(10,file='nonintrusive_coef_pod1.dat')
    read(10,*) (TNEW(i),i=1,NPERT)
    close(10)

    if(present(if_timetwo)) then
       open(10,file='nonintrusive_coef_pod2.dat')
       write(10,*) (TNEW(i),i=1,NPERT)
       close(10)   
       return
    endif
         
    DO IPERT=1,NPERT
       M0(:,IPERT) = tnew(:)
    END DO
    
    DO IPERT=1,NPERT            
       TOLD=TOLDKEEP_ORIG
       TOLD(IPERT)=TOLD(IPERT)+toler         
       open(10,file='nonintrusive_coef_pod0.dat')
       write(10,*) (TOLD(i),i=1,NPERT)
       close(10)
       !------run fluids------
       call fluids()
       !------end run fluids------
       open(10,file='nonintrusive_coef_pod1.dat')
       read(10,*) (TNEW(i),i=1,NPERT)
       close(10)
       M1(:,IPERT) = tnew(:)
    END DO
    
    DM0=(M1-M0)/toler
    
    deallocate(told)
    deallocate(tnew)
    deallocate(M0)
    deallocate(M1)

  END SUBROUTINE FIRST_ORDER_DM
  
  SUBROUTINE Hessian_DDM(state,pod_state,TOLDKEEP_ORIG,DM1,DDM) 
    type(state_type), dimension(:), intent(inout) :: state 
    type(state_type), dimension(:,:,:), intent(inout) :: pod_state
    INTEGER NONODS,NTIME,NPERT
    REAL TOLER
    REAL,dimension(:), intent(in) ::TOLDKEEP_ORIG
    REAL,dimension(:,:), intent(inout) ::DM1
    REAL,dimension(:,:,:), intent(inout) ::DDM
    PARAMETER(TOLER=0.00001)
    ! local variables...
    integer ipert,istate
    logical if_optimal
    real, DIMENSION( : ), allocatable :: told,tnew,TOLDkeep
    real, DIMENSION( :, : ), allocatable :: m0,m1
    real, dimension(:,:), allocatable ::  DM0,DM2
    real toler1,toler2
    integer KPERT,i,j,k
    
    ! DDM (Hessian) calculation: 
    NPERT=size(TOLDKEEP_ORIG)
    allocate(DM2(NPERT,NPERT))
    allocate(TOLDkeep(NPERT))
    
    DO KPERT=1,NPERT
       print*,'*************************KPERT111111111111111',KPERT
       DM1=0.0      
       TOLDkeep=TOLDKEEP_ORIG
       TOLDkeep(KPERT)=TOLDkeep(KPERT)-toler
       TOLD=TOLDkeep
       
       call First_Order_DM(state,pod_state,TOLD,DM1) 

       DM2=0.0          
       TOLDkeep=TOLDKEEP_ORIG
       TOLDkeep(KPERT)=TOLDkeep(KPERT)+toler           
       TOLD=TOLDkeep
       
       print*,'*************************KPERT2222222222222222222222222',KPERT
       call First_Order_DM(state,pod_state,TOLD,DM2) 
       DDM(:,:,kpert)=(DM2(:,:)-DM1(:,:))/(2*TOLER)
       print*,'DDM(:,:,kpert)',DDM(:,:,kpert)
    END DO
    close(10)

    deallocate(DM2)
    deallocate(TOLDkeep)
    
  END SUBROUTINE Hessian_DDM

   subroutine Non_Intrusive_ROM_main
    type(state_type), dimension(:), pointer :: state=> null()
    ! type(state_type), dimension(:), allocatable :: state_adj
    type(state_type), dimension(:,:,:), allocatable :: POD_state
    
    integer :: ierr,stat,i,j,k
    character(len = OPTION_PATH_LEN) :: simulation_name
    integer :: nsvd
    real :: valfun,fszero
    integer :: total_timestep,total_dump_no
    type(vector_field), pointer :: snapmean_velocity
    type(scalar_field), pointer :: snapmean_pressure
    type(vector_field), pointer :: u ,POD_velocity
    type(scalar_field), pointer :: p,POD_pressure
    real, dimension(:), allocatable :: pod_coef_optimised
    logical :: if_optimal,if_timetwo
    !  type(state_type), dimension(:,:,:) :: POD_state
    ! Change in pressure
    type(scalar_field) :: delta_p
    ! Change in velocity
    type(vector_field) :: delta_u
    real, dimension(:,:), allocatable :: pod_sol_velocity
    real, dimension(:), allocatable :: pod_sol_pressure

    ! non_intrusive_ROM
    type(state_type), dimension(:,:,:), allocatable :: Non_Intrusive_state
    real, dimension(:), allocatable ::  pod_coef_orig,pod_coef_new,pod_coef_old
    real, dimension(:,:), allocatable ::  DM0,DM1
    real, dimension(:,:,:), allocatable ::  DDM
    integer istate,NPERT,ntime
    integer dump_no_tmp, int_dump_period
    real current_time,finish_time,dt

    call delete_option("/reduced_model/execute_reduced_model") ! switch to full model
    !----------------------------------
    ! get State
    !----------------------------------
    call initialise_write_state
    ! Read state from .flml file
    call populate_state(state)
    u => extract_vector_field(state(1), "Velocity", stat)
    p => extract_scalar_field(state(1), "Pressure")

    !----------------------------------
    !----------------------------------

    call read_pod_basis_differntmesh(POD_state, state)
    allocate(pod_coef_orig((u%dim+1)*size(POD_state,1)))
    allocate(pod_coef_new((u%dim+1)*size(POD_state,1)))
    allocate(pod_coef_old((u%dim+1)*size(POD_state,1)))
    NPERT = (u%dim+1)*size(POD_state,1)
    allocate(DM0(NPERT,NPERT))
    allocate(DM1(NPERT,NPERT))
    allocate(DDM(NPERT,NPERT,NPERT))  

    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/finish_time", finish_time)       
    call get_option("/timestepping/timestep", dt)
    ntime=int((finish_time-current_time)/dt)

    !----------------------------------
    !get the initial ROM solution coef.
    !----------------------------------
    istate = 1
    call project_from_full_to_pod(istate, pod_state, state, pod_coef_orig)
    !----------------------------------------------
    !get the first order and second order matrices.
    !----------------------------------------------   
if(.true.) then
    call First_Order_DM(state,pod_state,pod_coef_orig,DM0) 
    open(10,file='DM.dat')
    do k=1,NPERT
       write(10,*) (DM0(i,k),i=1,NPERT)
    enddo
    close(10)
    DM1 = DM0
    call Hessian_DDM(state,pod_state,pod_coef_orig,DM1,DDM) 
    open(10,file='DDM.dat')
    do k=1,NPERT
       do j=1,NPERT
          write(10,*) (DDM(i,j,k),i=1,NPERT)
       enddo
    enddo
    close(10)
else
    open(10,file='DM.dat')
    do k=1,NPERT
       read(10,*) (DM0(i,k),i=1,NPERT)
    enddo
    close(10)

    open(10,file='DDM.dat')
    do k=1,NPERT
       do j=1,NPERT
          read(10,*) (DDM(i,j,k),i=1,NPERT)
       enddo
    enddo
    close(10)
endif
    !--------------------------------------------------------------
    ! take 2 time steps as we have an issue with initialization...
    !--------------------------------------------------------------
    if_timetwo=.true.
    call First_Order_DM(state,pod_state,pod_coef_orig,DM0, if_timetwo=if_timetwo)
    pod_coef_old = pod_coef_orig
    open(10,file='nonintrusive_coef_pod2.dat')
    read(10,*) (pod_coef_new(i),i=1,(u%dim+1)*size(POD_state,1))
    close(10)   

    call Nonintrusive(NPERT,ntime,DM0,DDM,pod_coef_orig,pod_coef_old,pod_coef_new)

             
    ! Allocate the change in pressure field
    call allocate(delta_p, p%mesh, "DeltaP")
    delta_p%option_path = trim(p%option_path)
    call zero(delta_p)
    
    ! Allocate the change in velocity pressure field
    call allocate(delta_u, u%dim, u%mesh, "DeltaU")
    delta_u%option_path = trim(u%option_path)
    call zero(delta_u)


    open(20,file='rom_solution_alltime.dat')

    if(have_option("/io/dump_period_in_timesteps/constant")) then
       call get_option("/io/dump_period_in_timesteps/constant", int_dump_period)
    else
       int_dump_period = 1
    endif

    do k = 0, ntime,int_dump_period
       print*,'**********************************',ntime,k
        read(20,*)(pod_coef_new(i),i=1,NPERT)
       call project_full(delta_u, delta_p, pod_sol_velocity, pod_sol_pressure, POD_state(:,:,istate), pod_coef_new)
       snapmean_velocity=>extract_vector_field(POD_state(1,1,istate),"SnapmeanVelocity")
       snapmean_pressure=>extract_Scalar_field(POD_state(1,2,istate),"SnapmeanPressure")
       !call addto(u, delta_u, dt)
       u%val=snapmean_velocity%val
       call addto(u, delta_u)
       p%val=snapmean_pressure%val
       call addto(p, delta_p)
       dump_no_tmp= int(k/int_dump_period)
       call write_state(dump_no_tmp, state)
    enddo
    close(20)
    !----------------------------------
    ! deallocate variables
    !----------------------------------
   if (allocated(pod_state)) then
       do i=1, size(pod_state,1)
          do j=1,size(pod_state,2)
             do k=1,size(pod_state,3)
                call deallocate(pod_state(i,j,k))
             enddo
          enddo
       end do
    end if

    deallocate(pod_coef_orig)
    deallocate(pod_coef_new)
    deallocate(pod_coef_old)
  end subroutine Non_Intrusive_ROM_main

 

subroutine function_val ( m, n, x, z ,p,total_timestep,podnum)

!*****************************************************************************80
!
!  F_SINR is a scalar function of an M-dimensional argument, to be interpolated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(M,N), the points.
!
!    Output, real ( kind = 8 ) Z(N), the value of the function at each point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) r(n)
  real ( kind = 8 ) x(m,n)
  real ( kind = 8 ) z(n)
      
  integer :: j,k,p
  integer :: total_timestep, podnum
  real :: pod_coef_all_obv(total_timestep, podnum)
   
  open(unit=61,file='coef_pod_all_obv')
  read(61,*)((pod_coef_all_obv(j,k),k=1,podnum),j=1,total_timestep)
  close(61)
   z(1:n)=pod_coef_all_obv(1:n,p)
  return
end
subroutine lagrange_interp_nd_value_smolyak ( m, n_1d, a, b, nd, zd, ni, xi, zi ,p)

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_VALUE evaluates an ND Lagrange interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N_1D(M), the order of the 1D rule to be used
!    in each dimension.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!
!    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
!
!    Input, real ( kind = 8 ) ZD(ND), the function evaluated at the points XD.
!
!    Input, integer ( kind = 4 ) NI, the number of points at which the 
!    interpolant is to be evaluated.
!
!    Input, real ( kind = 8 ) XI(M,NI), the points at which the interpolant is 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) ZI(NI), the interpolant evaluated at the 
!    points XI.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j 
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_1d(m)
  real ( kind = 8 ), allocatable :: value(:)
  real ( kind = 8 ) w(nd)
  real ( kind = 8 ), allocatable :: x_1d(:), pod_coef_all_obv(:,:)
  real :: xi(m,ni),xd(m,nd)
  real :: zd(nd)
  real ( kind = 8 ) zi(ni)
  integer ::k,p
  integer:: total_timestep, num_pod,nsvd
  real :: current_time,finish_time,dt
  call get_option(&
            '/reduced_model/pod_basis_formation/pod_basis_count', nsvd)  
  call get_option("/timestepping/current_time", current_time)
  call get_option("/timestepping/finish_time", finish_time)       
  call get_option("/timestepping/timestep", dt) 
  total_timestep= int((finish_time-current_time)/dt)
  allocate(pod_coef_all_obv(total_timestep, 3*nsvd))
  !total_timestep=2001
  num_pod=3*nsvd
  
   open(unit=61,file='coef_pod_all_obv')
   read(61,*)((pod_coef_all_obv(j,k),k=1,num_pod),j=1,total_timestep)
   close(61)
   n_1d(:)=3
   do j=1,m
   do k=1,nd
       xd(j,k)=pod_coef_all_obv(k,j) 
      !pod_coef_new(:)=pod_coef_all_obv(kk,:)
   enddo
   enddo
   do k=1,nd
       zd(k)=pod_coef_all_obv(K+1,p) 
   enddo
  do j = 1, ni

    w(1:nd) = 1.0D+00

    do i = 1, m
      n =  n_1d(i)
      allocate ( x_1d(1:n) )
      allocate ( value(1:n) )
      !call cc_compute_points ( n, x_1d )
      !x_1d(1:n) = 0.5D+00 * ( ( 1.0D+00 - x_1d(1:n) ) * a(i) &
      !                      + ( 1.0D+00 + x_1d(1:n) ) * b(i) )
     !  call lagrange_basis_1d ( n, x_1d, 1, xi(i,j), value )
      call lagrange_basis_1d ( n, xd(i,:), 1, xi(i,j), value ) 
      ! print *, 'value', value 
      call r8vec_direct_product2 ( i, n, value, m, nd, w )
     !  print *, 'value', value
      deallocate ( value )
      deallocate ( x_1d )
    end do
    !print *, 'weight', w
    !print *, 'zd', zd
    zi(j) = dot_product ( w, zd )
    print *, 'zj', zi(j)
  end do
  deallocate(pod_coef_all_obv)
  return
end 
END module Nonintrusive_rom_module

