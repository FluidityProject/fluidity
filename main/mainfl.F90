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

subroutine mainfl() bind(C)
  !!< This program solves the Navier-Stokes, radiation, and/or
  !!< advection-diffusion types of equations
  use fldebug
  use fluids_module
  !use reduced_fluids_module
  use signals
  use spud
  use tictoc
  use Adjoint_ROM_main_module
  use Nonintrusive_rom_module
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

  if(have_option("/reduced_model/adjoint").and. have_option("/reduced_model/execute_reduced_model")) then
     !######################################################
     !       Reduced Fluidity Model
     !######################################################
    ewrite(1, *) 'entry adjoint model'

    call Adjoint_ROM_main()
  else if(have_option("/reduced_model/nonintrusive_rom")) then
     !######################################################
     !       Nonintrusive Reduced Model
     !######################################################
    ewrite(1, *) 'entry nonintrusive reduced order model'

    call Non_Intrusive_ROM_main()
    
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

end subroutine mainfl
