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

subroutine test_harmonic_analysis

  use navsto_module
  use Tidal_module
  use fldebug
  use unittest_tools
  
  implicit none

  integer :: M, N, MM, i, j
  real :: sigma(11), amplitude(11), phase(11), amp(11), phs(11), dt
  real :: pi = 3.141592654, tol = 1.0e-8
  real, dimension(:,:), allocatable :: harmonic_A  ! for solving Ax=b system
  real, dimension(:), allocatable :: harmonic_x, harmonic_b, harmonic_time_series_vals_at_node, harmonic_times_reordered  
  character(len=3), dimension(11) :: constituents=(/"M2 ","S2 ","N2 ","K2 ","K1 ","O1 ","P1 ","Q1 ","Mf ","Mm ","SSa"/)            
  logical :: fail, warn  
  
  
  amplitude(1:11)=0.0
  amp(1:11)=0.0
  phase(1:11)=0.0
  phs(1:11)=0.0
  
  M = 2 ! Number of constituents we're analysing for 
  sigma(1) = get_tidal_frequency(constituents(1))/(2.*pi)   !M2
  sigma(2) = get_tidal_frequency(constituents(5))/(2.*pi)   !K1

!Values for the synthetic data
  amplitude(1) = 2.5
  amplitude(2) = 0.75
  phase(1) = 0.33
  phase(2) = 0.11
  
    
  N = 50 ! length of the time series        
  allocate(harmonic_A(2*M+1,2*M+1))
  allocate(harmonic_x(2*M+1))
  allocate(harmonic_b(2*M+1))
  allocate(harmonic_time_series_vals_at_node(N))
  allocate(harmonic_times_reordered(N))  
    
  
! synthetic time series  
  dt = 3600.0
  do i = 1,N
    harmonic_times_reordered(i) = (i-1)*dt
    harmonic_time_series_vals_at_node(i) = 0.0
    do j = 1,M
      harmonic_time_series_vals_at_node(i) = harmonic_time_series_vals_at_node(i) + &
                         amplitude(j)*cos(2.0*pi*(sigma(j)*harmonic_times_reordered(i) - phase(j)))
    end do 
  end do   

! Form and invert the least squares system
  call harmonic_analysis_at_single_node(N,harmonic_times_reordered,harmonic_time_series_vals_at_node,M,sigma,&
                                                harmonic_A,harmonic_x,harmonic_b)

! Extract the amplitude and phase from the solution to the linear system
  do i = 1,M
    amp(i) = sqrt( harmonic_x(i+1)**2 + harmonic_x(i+1+M)**2 )
    phs(i) = atan2(harmonic_x(i+1+M),harmonic_x(i+1))/(2.*pi)
    print *,amp(i),phs(i),amplitude(i),phase(i)
  end do  
   
! Perform and report the test  
  fail = .false.; warn = .false.
  if (any(abs(amp - amplitude)>tol)) fail = .true.
  call report_test("[Correct amplitudes]", fail, warn, "An amplitude is in error")
  fail = .false.; warn = .false.
  if (any(abs(phs - phase)>tol)) fail = .true.
  call report_test("[Correct phases]", fail, warn, "A phase is in error")  

  deallocate(harmonic_A)
  deallocate(harmonic_x)
  deallocate(harmonic_b)    
  deallocate(harmonic_time_series_vals_at_node)
  deallocate(harmonic_times_reordered)
    
end subroutine test_harmonic_analysis
