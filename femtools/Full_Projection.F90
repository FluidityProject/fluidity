!  Copyright (C) 2006 Imperial College London and others.
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

  module Full_Projection
    use FLDebug
    use elements
    use Petsc_tools
    use Solvers
    use Signal_Vars
    use Sparse_Tools
    use sparse_tools_petsc
    use Sparse_matrices_fields
    use Fields_Base
    use flcomms_module
    use Global_Parameters
    use spud
    use halos
    use Multigrid

#ifdef HAVE_PETSC_MODULES
    use petsc
    use petscvec
    use petscmat
    use petscksp
    use petscpc
#endif

    implicit none
    ! Module to provide solvers, preconditioners etc... for full_projection Solver.
    ! Not this is currently tested for Full CMC solves and Stokes flow:
#ifdef HAVE_PETSC
#ifdef HAVE_PETSC_MODULES
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscpcdef.h"
#else
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#endif
#include "petscversion.h"
#else
#define PetscReal real
#define PetscInt integer
#endif
    
    private
    
    public petsc_solve_full_projection
    
  contains
    
!--------------------------------------------------------------------------------------------------------------------
    subroutine petsc_solve_full_projection(x,ctp_m,inner_m,ct_m,rhs,pmat,kmk_matrix)
!--------------------------------------------------------------------------------------------------------------------

      ! Solve Schur complement problem the nice way (using petsc) !!
      ! Arguments that come in:
      type(scalar_field), intent(inout) :: x ! Fluidity Solution vector 
      type(scalar_field), intent(inout)  :: rhs ! Fluidity RHS vector
      type(petsc_csr_matrix), intent(inout) :: inner_m ! Inner matrix of schur complement
      type(block_csr_matrix), intent(in) :: ctp_m, ct_m ! Fluidity Compressible / Incompressible Divergence matrices.
      ! Fluidity Preconditioner matrix. If preconditioner matrix is set as LumpedSchurComplement, this array comes in as
      ! CMC where M is the inverse_lumped mass (ideal for Full_CMC solve). However, if preconditioner matrix is specified 
      ! as DiagonalSchurComplement, this comes in as C(Big_m_id)C, where Big_M_id is the inverse diagonal of the full 
      ! momentum matrix:
      type(csr_matrix), intent(inout) :: pmat 
      ! p1-p1 stabilization matrix:
      type(csr_matrix), pointer :: kmk_matrix


#ifdef HAVE_PETSC
#if PETSC_VERSION_MAJOR==3
      KSP ksp ! Object type for outer solve (i.e. A * delta_p = rhs)
      Mat A ! PETSc Schur complement matrix (i.e. G^t*m^-1*G) 
      Vec y, b ! PETSc solution vector (y), PETSc RHS (b)
      
      character(len=OPTION_PATH_LEN) solver_option_path, inner_solver_option_path, name
      integer literations
      type(petsc_numbering_type) petsc_numbering
      logical lstartfromzero
      
      ewrite(2,*) 'Entering Full Projection Solve:'
      ewrite(2,*) 'Inner matrix of Schur complement: ',inner_m%name 
      ewrite(2,*) 'Are we applying stabilization: ', associated(kmk_matrix)

      ! Start from zero? For now, always true.
      ewrite(2,*) 'Note: startfromzero hard-coded to .true.'
      ewrite(2,*) 'Ignoring setting from solver option.'
      lstartfromzero=.true.

      ! Convert matrices to PETSc format, setup PETSc object and PETSc numbering, then
      ! Build Schur complement and set KSP.
      ewrite(2,*) 'Entering PETSc setup for Full Projection Solve'
      call petsc_solve_setup_full_projection(y,A,b,ksp,petsc_numbering,name,solver_option_path, &
           inner_solver_option_path,lstartfromzero,inner_m,ctp_m,ct_m,x%option_path,pmat,rhs, &
           kmk_matrix)

      ewrite(2,*) 'Create RHS and solution Vectors in PETSc Format'
      ! create PETSc vec for rhs using above numbering:
      call field2petsc(rhs, petsc_numbering, b)
      ewrite(1, *) 'RHS assembly completed.'
      if (.not. lstartfromzero) then
         ewrite(1, *) 'Assembling initial guess.'
         ! create PETSc vec for initial guess and result using above numbering:
         call field2petsc(x, petsc_numbering, y)
         ewrite(1, *) 'Initial guess assembly completed.'
      end if

      ewrite(2,*) 'Entering Core PETSc Solve'
      ! Solve Ay = b using KSP and PC. Also check convergence. We call this the inner solve.
      call petsc_solve_core(y, A, b, ksp, petsc_numbering, name, solver_option_path, lstartfromzero, &
           literations, x0=x%val)

      ewrite(2,*) 'Copying PETSc solution vector into designated Fluidity array'
      ! Copy back the result into the fluidity solution array (x) using the PETSc numbering:
      call petsc2field(y, petsc_numbering, x, rhs)
      
      ewrite(2,*) 'Destroying all PETSc objects'
      ! Destroy all PETSc objects and the petsc_numbering:
      call petsc_solve_destroy(y, A, b, ksp, petsc_numbering)
#else
      FLAbort("PETSC_solve_full_projection only works with PETSc 3.0")
#endif
#else
      FLAbort("PETSc_solve_Full_Projection called while not configured with PETSc")
#endif
      
      ewrite(2,*) 'Leaving PETSc_solve_setup_full_projection'
      
    end subroutine petsc_solve_full_projection

#ifdef HAVE_PETSC
#if PETSC_VERSION_MAJOR==3
!--------------------------------------------------------------------------------------------------------
    subroutine petsc_solve_setup_full_projection(y,A,b,ksp,petsc_numbering_p,name,solver_option_path, &
         inner_solver_option_path,lstartfromzero,inner_m,div_matrix_comp, &
         div_matrix_incomp,option_path,preconditioner_matrix,rhs,kmk_matrix)
!--------------------------------------------------------------------------------------------------------

      !  Sets up things needed to call petsc_solve_core. The following arrays/arguments go out:
      !
      !  PETSc solution vector:
      Vec, intent(out) :: y
      !  PETSc Schur Complement Matrix:
      Mat, intent(out) :: A
      !  PETSc rhs vector:
      Vec, intent(out) :: b
      !  Solver object:
      Mat, intent(out) :: ksp
      ! Numbering from local to PETSc:
      type(petsc_numbering_type), intent(out) :: petsc_numbering_p
      ! Name of solve, to be printed on log output:
      character(len=*), intent(out) :: name
      ! Solver option paths:
      character(len=*), intent(out) :: solver_option_path
      character(len=*), intent(out) :: inner_solver_option_path

      ! Stuff that comes in:
      !
      ! Full matrix for CMC / Stokes solve:
      type(petsc_csr_matrix), intent(inout):: inner_m
      ! Divergence matrices:
      type(block_csr_matrix), intent(in):: div_matrix_comp
      type(block_csr_matrix), intent(in):: div_matrix_incomp
      ! Whether to start from zero initial guess (passed in):
      logical, intent(inout):: lstartfromzero
      ! Fluidity RHS:
      type(scalar_field), intent(inout) :: rhs
      ! Preconditioner matrix:
      type(csr_matrix), intent(inout) :: preconditioner_matrix
      ! Stabilization matrix:
      type(csr_matrix), pointer :: kmk_matrix
      ! Option path:
      character(len=*), intent(in):: option_path

      ! Additional arrays used internally:
      ! Additional numbering types:
      type(petsc_numbering_type) :: petsc_numbering_u

      integer, dimension(:), allocatable:: ghost_nodes

      PetscErrorCode ierr
      KSP ksp_schur ! ksp object for the inner solve in the Schur Complement
      Mat G_t_comp ! PETSc compressible divergence matrix
      Mat G_t_incomp ! PETSc incompressible divergence matrix
      Mat G ! PETSc Gradient matrix (transpose of G_t_incomp)
      Mat S ! PETSc Stabilization matrix (kmk_matrix)
      Mat pmat ! PETSc preconditioning matrix
      
      integer reference_node, stat
      logical parallel

      ! Sort option paths etc...
      solver_option_path=complete_solver_option_path(option_path)
      inner_solver_option_path= trim(option_path)//&
                                  "/prognostic/scheme/use_projection_method&
                                  &/full_schur_complement/inner_matrix[0]/solver"
      
      if (have_option(trim(option_path)//'/name')) then
         call get_option(trim(option_path)//'/name', name)
         ewrite(1,*) 'Inside petsc_solve_(block_)csr, solving for: ', trim(name)
      else
         ewrite(1,*) 'Inside petsc_solve_(block_)csr, solving using option_path: ', trim(option_path)
         name=option_path
      end if

      ! Impose reference pressure node on processor 1 only:
      call get_option(trim(option_path)//&
           '/prognostic/reference_node', reference_node, stat=stat)
      if (stat==0.AND.GetProcNo()==1) then
         ewrite(2,*) 'Imposing_reference_pressure_node'        
         allocate(ghost_nodes(1:1))
         ghost_nodes(1) = reference_node
         call set(rhs,reference_node,0.0) ! Modify RHS accordingly
      else
         allocate(ghost_nodes(1:0))
      end if

      ! set up numbering used in PETSc objects:
      call allocate(petsc_numbering_u, &
           nnodes=block_size(div_matrix_comp,2), nfields=blocks(div_matrix_comp,2), &
           halo=div_matrix_comp%sparsity%column_halo)
      call allocate(petsc_numbering_p, &
           nnodes=block_size(div_matrix_comp,1), nfields=1, &
           halo=preconditioner_matrix%sparsity%row_halo, ghost_nodes=ghost_nodes)

      ! Convert Divergence matrix (currently stored as block_csr matrix) to petsc format:   
      ! Create PETSc Div Matrix (comp & incomp) using this numbering:
      G_t_comp=block_csr2petsc(div_matrix_comp, petsc_numbering_p, petsc_numbering_u)
      G_t_incomp=block_csr2petsc(div_matrix_incomp, petsc_numbering_p, petsc_numbering_u)

      ! Scale G_t_comp to fit PETSc sign convention:
      call MatScale(G_t_comp,-1.,ierr)

      ! Determine transpose of G_t_incomp to form Gradient Matrix (G):
      call MatCreateTranspose(G_t_incomp,G,ierr)

      ! Convert Stabilization matrix --> PETSc format if required:
      if(associated(kmk_matrix)) then
         S=csr2petsc(kmk_matrix, petsc_numbering_p, petsc_numbering_p)
      end if
      
      ! Need to assemble the petsc matrix before we use it:
      call assemble(inner_M)

      ! Build Schur complement:
      ewrite(2,*) 'Building Schur complement'                
      if(associated(kmk_matrix)) then
         call MatCreateSchurComplement(inner_M%M,G,G_t_comp,S,A,ierr)
      else
         call MatCreateSchurComplement(inner_M%M,G,G_t_comp,PETSC_NULL_OBJECT,A,ierr)
      end if
      
      ! Set ksp for M block solver inside the Schur Complement (the inner, inner solve!). 
      call MatSchurComplementGetKSP(A,ksp_schur,ierr)
      call setup_ksp_from_options(ksp_schur, inner_M%M, inner_M%M, &
        inner_solver_option_path, startfromzero_in=.true.)

      ! Assemble preconditioner matrix in petsc format:
      pmat=csr2petsc(preconditioner_matrix, petsc_numbering_p, petsc_numbering_p)

      ! Set up RHS and Solution vectors (note these are loaded later):
      b = PetscNumberingCreateVec(petsc_numbering_p)
      call VecDuplicate(b,y,ierr) ! Duplicate vector b and form vector y

      ! Schur complement objects now fully set up. The next stage is to setup
      ! KSP and PC for the solution of Ay = b, where A = G_t*M^-1*G.

      parallel=IsParallel()

      call SetupKSP(ksp,A,pmat,solver_option_path,parallel,lstartfromzero)
      
      ! Destroy the matrices setup for the schur complement computation. While
      ! these matrices are destroyed here, they are still required for the inner solve,
      ! since the schur complement depends on them. As a result, even though they're destroyed
      ! here, the PETSc reference count ensures that they are not wiped from memory.
      ! In other words, references to these arrays remain intact. 

      call MatDestroy(G_t_comp,ierr) ! Destroy Compressible Divergence Operator.
      call MatDestroy(G_t_incomp, ierr) ! Destroy Incompressible Divergence Operator.
      call MatDestroy(G, ierr) ! Destroy Gradient Operator (i.e. transpose of incompressible div).
      call MatDestroy(pmat,ierr) ! Destroy preconditioning matrix if allocated.
      if(associated(kmk_matrix)) call MatDestroy(S,ierr) ! Destroy stabilization matrix
      
      call deallocate( petsc_numbering_u )
      ! petsc_numbering_p is passed back and destroyed there      

    end subroutine petsc_solve_setup_full_projection

#endif
#endif

  end module Full_Projection
  
