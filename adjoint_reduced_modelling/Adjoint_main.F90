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
program Adjoint_main
  use spud
  use fields
  use state_module
  use write_state_module
  use timeloop_utilities
  use global_parameters, only: option_path_len, current_time, dt, FIELD_NAME_LEN
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

  implicit none
#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif

  type(state_type), dimension(:), pointer :: state
  type(state_type), dimension(:), allocatable :: state_adj

  integer :: timestep
  integer :: ierr,stat


  character(len = OPTION_PATH_LEN) :: simulation_name
!==============================================================================
!Initialize 
!==============================================================================

!--------------------------------------------------------------------------
!Initialize the L-BFGS optimization and Gradient generation
!--------------------------------------------------------------------------

  call mainfl_forward()

!  call python_init()
!  call read_command_line()
  !!state includes the fields i.c., b.c.and mesh info in flml file
!  call populate_state(state)
!  call get_option('/simulation_name',simulation_name)

!--------------------------------------------------------------------------
!Initial conditions to purturbation to generate observations
!--------------------------------------------------------------------------
!	CALL initialCondition
!	CALL perturbInitialCondition

!--------------------------------------------------------
!  Setup observational data
!--------------------------------------------------------

  !run the forward model to get the pseudo-observational data
	!CALL generateFullObervations

!---------------------------------------------------------------------------
!IC of full model with perturbations as initial guess of optmization
!---------------------------------------------------------------------------

	!eps_init = 5.0D-3
	!CALL initialGuessOptimization

!--------------------------------------------------------------------------
!From model solutions and observations, we get forcing terms and old cost
!--------------------------------------------------------------------------

!----------------------------------------------------------------------------
!From model solutions, forcing term and adjoint model, we compute the gradient
!----------------------------------------------------------------------------
  ! call adjoint_model

!--------------------------------------------------------------------------
!DATA ASSIMILLATION STARTS
!--------------------------------------------------------------------------

	!PRINT *, "DATA ASSIMILLATION STARTS HERE"
	!WRITE(*,200)

 !itera = 1

 !222  CONTINUE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The call to the L-BFGS-B optimization subroutine SETULB
! We use the L-BFGS-B implementation by Ciyou Zhu,
! Richard Byrd and Jorge Nocedal.:
!
!   subroutine setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa,
!     +                 task, iprint, csave, lsave, isave, dsave)
!
!
!     NEOS, November 1994. (Latest revision April 1997.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!     References:
!
!       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
!       memory algorithm for bound constrained optimization'',
!       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!
!       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
!       limited memory FORTRAN code for solving bound constrained
!       optimization problems'', Tech. Report, NAM-11, EECS Department,
!       Northwestern University, 1994.
!
!                           *  *  *
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!--------------------------------------------------------------------------
!Initialize and update the gradient of control variables for LBFGS-B
!--------------------------------------------------------------------------

!	ctrgrad(1:nvert) = psi_grad(1:nvert)
!	ctrgrad(nvert+1:2*nvert) = u_grad(1:nvert)
!	ctrgrad(2*nvert+1:3*nvert) = v_grad(1:nvert)
!--------------------------------------------------------------------------
!setulb requests the value of f and the gradient g at the current point
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!Update the psi,u,and v by new point of control variables 
!--------------------------------------------------------------------------


!--------------------------------------------------------------------------
!Run the model again to get the new state variables and new model solutions
!--------------------------------------------------------------------------
!  call forward model

!--------------------------------------------------------------------------
!From new model solutions, we update forcing terms and cost function
!--------------------------------------------------------------------------
! call subroutine cost function


!--------------------------------------------------------------------------
! if using L-BFGS-B, here run the adjoint modelling again
!From new model solutions, forcing terms and adjoint model,
!we compute the new gradient
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!Terminate the run if the total number of function/gradient
!Evaluation exceeds the threshold value MaxFGeval
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!Go back to the minimization subroutine with a new iterate
!--------------------------------------------------------------------------

!                itera = itera+1

!		GOTO 222

!--------------------------------------------------------------------------
!DATA ASSIMILLATION END
!--------------------------------------------------------------------------

contains

subroutine mainfl_forward() bind(C)
  !!< This program solves the Navier-Stokes, radiation, and/or
  !!< advection-diffusion types of equations
  use fldebug
  use fluids_module
  !use reduced_fluids_module
  use signals
  use spud
  use tictoc
#ifdef HAVE_ZOLTAN
  use zoltan
#endif
  
  implicit none

  ! We need to do this here because the fortran Zoltan initialisation
  ! routine does extra things on top of the C one. That wasn't a fun
  ! hour's debugging ...
#ifdef HAVE_ZOLTAN
  real(zoltan_float) :: ver
  integer(zoltan_int) :: ierr

  ierr = Zoltan_Initialize(ver)  
  assert(ierr == ZOLTAN_OK)
#endif
  
  ! Establish signal handlers
  call initialise_signals()

  call tictoc_reset()
 
  if(have_option("/model/fluids/pod")) then
     !######################################################
     !       Reduced Fluidity Model
     !######################################################
  
     FLExit("POD is disabled")      
     !call reducedfluids(filename, filename_len)
  else
     !######################################################
     !      Normal Fluidity Model
     !######################################################
     
     call tic(TICTOC_ID_SIMULATION)
     ewrite(1, *) "Calling fluids from mainfl"
     call fluids()
     ewrite(1, *) "Exited fluids"
     call toc(TICTOC_ID_SIMULATION)
     call tictoc_report(2, TICTOC_ID_SIMULATION)
     
  end if

  if(SIG_INT) then
    FLExit("Interrupt signal received")
  end if

end subroutine mainfl_forward

subroutine mainfl_adj() bind(C)
  !!< This program solves the Navier-Stokes, radiation, and/or
  !!< advection-diffusion types of equations
  use fldebug
  use fluids_module
  !use reduced_fluids_module
  use signals
  use spud
  use tictoc
#ifdef HAVE_ZOLTAN
  use zoltan
#endif
  
  implicit none

  ! We need to do this here because the fortran Zoltan initialisation
  ! routine does extra things on top of the C one. That wasn't a fun
  ! hour's debugging ...
#ifdef HAVE_ZOLTAN
  real(zoltan_float) :: ver
  integer(zoltan_int) :: ierr
  logical :: adjoint

  ierr = Zoltan_Initialize(ver)  
  assert(ierr == ZOLTAN_OK)
#endif
  
  ! Establish signal handlers
  call initialise_signals()

  call tictoc_reset()
 
  if(have_option("/model/fluids/pod")) then
     !######################################################
     !       Reduced Fluidity Model
     !######################################################
  
     FLExit("POD is disabled")      
     !call reducedfluids(filename, filename_len)
  else
     !######################################################
     !      Normal Fluidity Model
     !######################################################
     
     call tic(TICTOC_ID_SIMULATION)
     ewrite(1, *) "Calling fluids from mainfl"
#ifdef HAVE_MPI
     call mpi_init(ierr)
#endif
     
#ifdef HAVE_PETSC
     call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif
     
     call python_init()
     call read_command_line()
!     call get_option(&
!          '/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
     if (have_option("/reduced_model/execute_reduced_model")) adjoint = .true.
     call fluids(adjoint)
     ewrite(1, *) "Exited fluids"
     call toc(TICTOC_ID_SIMULATION)
     call tictoc_report(2, TICTOC_ID_SIMULATION)
     
  end if

  if(SIG_INT) then
    FLExit("Interrupt signal received")
  end if

end subroutine mainfl_adj

  subroutine read_command_line()
    implicit none
    ! Read the input filename.
    character(len=1024) :: argument, filename
    integer :: status, argn, level

    call set_global_debug_level(0)

    argn=1
    do 

       call get_command_argument(argn, value=argument, status=status)
       argn=argn+1

       if (status/=0) then
          call usage
          stop
       end if

       if (argument=="-v") then
          call get_command_argument(argn, value=argument, status=status)
          argn=argn+1

          if (status/=0) then
             call usage
             stop
          end if

          read(argument, "(i1)", err=666) level
          call set_global_debug_level(level)

       else

          ! Must be the filename
          filename=argument

       end if

       if (argn>=command_argument_count()) exit
    end do

    call load_options(filename)

    return

666 call usage
    stop

  end subroutine read_command_line


  subroutine usage
    implicit none

    write (0,*) "usage: form_pod_basis [-v n] <options_file>"
    write (0,*) ""
    write (0,*) "-v n sets the verbosity of debugging"
  end subroutine usage


end program Adjoint_main
