! Copyright (C) 2009 Imperial College London and others.
! 
! Please see the AUTHORS file in the main source directory for a full list
! of copyright holders.
! 
! Gerard Gorman
! Applied Modelling and Computation Group
! Department of Earth Science and Engineering
! Imperial College London
! 
! g.gorman@imperial.ac.uk
! 
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License.
! 
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
! USA

#include "confdefs.h"

module AdaptProgress
  implicit none
  
  private
  
  public::initialise, finalize, should_exit
#ifdef HAVE_MPI
  include 'mpif.h'
#endif
  logical::initialised=.false.
  integer, dimension(:), allocatable::load, rrequest
  integer::myrank, nprocs
  real::imbalance_tol=0.5
  integer win

contains
  subroutine initialise(count, tol)
    integer, intent(in)::count
    real, intent(in)::tol
#ifdef HAVE_MPI
    integer have_mpi_init, i, ierr
    if(.not.initialised) then
       call MPI_Initialized(have_mpi_init, ierr)
       if(have_mpi_init.eq.0) then
          nprocs=1
       else
          call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
       endif

       if(nprocs.gt.1) then
          call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
          allocate(load(0:nprocs-1))
          call MPI_Allgather(count, 1, MPI_INTEGER, load, 1, &
               MPI_INTEGER, MPI_COMM_WORLD, ierr) 
          
          allocate(rrequest(0:nprocs-1))
          
          do i=0, nprocs-1
             if(i.ne.myrank) then
                call MPI_Irecv(load(i), 1, MPI_INTEGER, i, 1, &
                     MPI_COMM_WORLD, rrequest(i), ierr)
             else
                rrequest(i) = MPI_REQUEST_NULL
             end if
          end do

          imbalance_tol = tol
       end if

       initialised = .true.
    end if
#endif
  end subroutine initialise

  subroutine finalize(count)
    integer, intent(in)::count
#ifdef HAVE_MPI
    integer i, ierr
    integer, allocatable, dimension(:)::request
    integer, allocatable, dimension(:, :)::status

    if(nprocs.gt.1) then
       allocate(request(0:nprocs-1))
       allocate(status(MPI_STATUS_SIZE, 0:nprocs-1))
       
       do i=0, nprocs-1
          if(i.ne.myrank) then
             call MPI_Isend(count, 1, MPI_INTEGER, i, 1, &
                  MPI_COMM_WORLD, request(i), ierr)
          else
             request(i) = MPI_REQUEST_NULL
          end if
       end do
       
       call MPI_Waitall(nprocs, rrequest, status, ierr)
       call MPI_Waitall(nprocs, request, status, ierr)
       
       deallocate(request, rrequest, status, load)
    end if
    initialised = .false.
#endif
  end subroutine finalize
  
  logical function should_exit(count)
    integer, intent(in)::count
#ifdef HAVE_MPI
    real imbalance
    integer ierr
    
    if(nprocs.gt.1) then
       load(myrank) = count
       
       imbalance = 1.0 - real(sum(load))/(nprocs*maxval(load))
       
       should_exit = (imbalance>imbalance_tol)
    else
       should_exit = .false.
    end if
#else
    should_exit = .false.
#endif
  end function should_exit
  
end module AdaptProgress
