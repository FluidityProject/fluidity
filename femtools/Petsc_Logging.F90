!    Copyright (C) 2006 Imperial College London and others.
!
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
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
module petsc_logging
#ifdef HAVE_PETSC_MODULES
  use petsc 
#endif
  implicit none

#include "petsc_legacy.h"

  private
  public petsc_stage_register, petsc_event_begin, petsc_event_end, &
    petsc_stage_begin, petsc_stage_end, petsc_events_initialise

  PetscClassId, public :: petsc_class_id_fluidity = -1
  PetscLogEvent, public :: petsc_event_element_assembly = -1
  PetscLogEvent, public :: petsc_event_surface_element_assembly = -1
  PetscLogEvent, public :: petsc_event_halo_update = -1
  PetscLogEvent, public :: petsc_event_mesh_colouring = -1

  PetscLogEvent, public :: petsc_event_csr_mult = -1
  PetscLogEvent, public :: petsc_event_csr_mult_addto = -1
  PetscLogEvent, public :: petsc_event_csr_mult_t = -1
  PetscLogEvent, public :: petsc_event_csr_mult_t_addto = -1

  PetscLogEvent, public :: petsc_event_make_sparsity = -1
  PetscLogEvent, public :: petsc_event_make_sparsity_t = -1
  PetscLogEvent, public :: petsc_event_make_sparsity_mult = -1
  PetscLogEvent, public :: petsc_event_make_sparsity_dg_mass = -1
  PetscLogEvent, public :: petsc_event_make_sparsity_compactdgdouble = -1

  PetscLogEvent, public :: petsc_event_csr_zero = -1
  PetscLogEvent, public :: petsc_event_block_csr_zero = -1


  contains

  subroutine petsc_stage_register(name, stage)
    character(len=*), intent(in):: name
    PetscLogStage, intent(inout) ::  stage

    PetscErrorCode:: ierr

      if (stage<0) then
        call PetscLogStageRegister(name, stage, ierr)
      end if

  end subroutine petsc_stage_register

  subroutine petsc_stage_begin(stage)
     PetscLogStage, intent(in) :: stage

     PetscErrorCode :: ierr
     call PetscLogStagePush(stage, ierr)

  end subroutine petsc_stage_begin

  subroutine petsc_stage_end()

     PetscErrorCode :: ierr
     call PetscLogStagePop(ierr)

  end subroutine petsc_stage_end

  subroutine petsc_events_initialise

    PetscErrorCode :: ierr
    PetscClassId :: class_id

    call PetscClassIDRegister("Fluidity", petsc_class_id_fluidity, ierr)
    class_id = petsc_class_id_fluidity
    call PetscLogEventRegister("Element assembly", class_id, petsc_event_element_assembly, ierr)
    call PetscLogEventRegister("Surf. elem. ass.", class_id, petsc_event_surface_element_assembly, ierr)
    call PetscLogEventRegister("Halo update", class_id, petsc_event_halo_update, ierr)
    call PetscLogEventRegister("Mesh colouring", class_id, petsc_event_mesh_colouring, ierr)

    call PetscLogEventRegister("CSR Mult", class_id, petsc_event_csr_mult, ierr)
    call PetscLogEventRegister("CSR Mult Add", class_id, petsc_event_csr_mult_addto, ierr)
    call PetscLogEventRegister("CSR Mult_t", class_id, petsc_event_csr_mult_t, ierr)
    call PetscLogEventRegister("CSR Mult_t Add", class_id, petsc_event_csr_mult_t_addto, ierr)

    call PetscLogEventRegister("Make sparsity", class_id, petsc_event_make_sparsity, ierr)
    call PetscLogEventRegister("Mk sparsity_t", class_id, petsc_event_make_sparsity_t, ierr)
    call PetscLogEventRegister("Mk sparsity_mult", class_id, petsc_event_make_sparsity_mult, ierr)
    call PetscLogEventRegister("Mk spars.dg_mass", class_id, petsc_event_make_sparsity_dg_mass, ierr)
    call PetscLogEventRegister("Mk spars.cdg_dbl", class_id, petsc_event_make_sparsity_compactdgdouble, ierr)

    call PetscLogEventRegister("CSR zero", class_id, petsc_event_csr_zero, ierr)
    call PetscLogEventRegister("CSR block zero", class_id, petsc_event_block_csr_zero, ierr)

  end subroutine petsc_events_initialise

  subroutine petsc_event_begin(event)
     PetscLogEvent, intent(in) :: event

     PetscErrorCode :: ierr
     call PetscLogEventBegin(event, ierr)

  end subroutine petsc_event_begin

  subroutine petsc_event_end(event)
     PetscLogEvent, intent(in) :: event

     PetscErrorCode :: ierr
     call PetscLogEventEnd(event, ierr)

  end subroutine petsc_event_end
  
end module petsc_logging
