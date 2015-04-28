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
    use Global_Parameters
    use spud
    use halos
    use Multigrid
    use state_module
    use petsc_solve_state_module
    use boundary_conditions_from_options

#ifdef HAVE_PETSC_MODULES
    use petsc
#endif

    implicit none
    ! Module to provide solvers, preconditioners etc... for full_projection Solver.
    ! Not this is currently tested for Full CMC solves and Stokes flow:
#include "petsc_legacy.h"
    
    private
    
    public petsc_solve_full_projection
    
  contains
    
!--------------------------------------------------------------------------------------------------------------------
    subroutine petsc_solve_full_projection(x,ctp_m,inner_m,ct_m,rhs,pmat,&
      state, inner_mesh, auxiliary_matrix)
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
      ! momentum matrix. If preconditioner is set to ScaledPressureMassMatrix, this comes in as the pressure mass matrix,
      ! scaled by the inverse of viscosity.
      type(csr_matrix), intent(inout) :: pmat
      ! state, and inner_mesh are used to setup mg preconditioner of inner solve
      type(state_type), intent(in):: state
      type(mesh_type), intent(in):: inner_mesh
      ! p1-p1 stabilization matrix or free surface terms:
      type(csr_matrix), optional, intent(in) :: auxiliary_matrix

      KSP ksp ! Object type for outer solve (i.e. A * delta_p = rhs)
      Mat A ! PETSc Schur complement matrix (i.e. G^t*m^-1*G) 
      Vec y, b ! PETSc solution vector (y), PETSc RHS (b)
      
      character(len=OPTION_PATH_LEN) solver_option_path, name
      integer literations
      type(petsc_numbering_type) petsc_numbering
      logical lstartfromzero
      
      ewrite(2,*) 'Entering Full Projection Solve:'
      ewrite(2,*) 'Inner matrix of Schur complement: ',inner_m%name 

      ! Start from zero? For now, always true.
      ewrite(2,*) 'Note: startfromzero hard-coded to .true.'
      ewrite(2,*) 'Ignoring setting from solver option.'
      lstartfromzero=.true.

      ! Convert matrices to PETSc format, setup PETSc object and PETSc numbering, then
      ! Build Schur complement and set KSP.
      ewrite(2,*) 'Entering PETSc setup for Full Projection Solve'
      call petsc_solve_setup_full_projection(y,A,b,ksp,petsc_numbering,name,solver_option_path, &
           lstartfromzero,inner_m,ctp_m,ct_m,x%option_path,pmat, &
           rhs, state, inner_mesh, auxiliary_matrix)

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
      call petsc_solve_core(y, A, b, ksp, petsc_numbering, solver_option_path, lstartfromzero, &
           literations, sfield=x, x0=x%val, nomatrixdump=.true.)

      ewrite(2,*) 'Copying PETSc solution vector into designated Fluidity array'
      ! Copy back the result into the fluidity solution array (x) using the PETSc numbering:
      call petsc2field(y, petsc_numbering, x, rhs)
      
      ewrite(2,*) 'Destroying all PETSc objects'
      ! Destroy all PETSc objects and the petsc_numbering:
      call petsc_solve_destroy(y, A, b, ksp, petsc_numbering, solver_option_path)
      
      ewrite(2,*) 'Leaving PETSc_solve_full_projection'
      
    end subroutine petsc_solve_full_projection

!--------------------------------------------------------------------------------------------------------
    subroutine petsc_solve_setup_full_projection(y,A,b,ksp,petsc_numbering_p,name,solver_option_path, &
         lstartfromzero,inner_m,div_matrix_comp, div_matrix_incomp,option_path,preconditioner_matrix,rhs, &
         state, inner_mesh, auxiliary_matrix)
         
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
      type(csr_matrix), optional, intent(in) :: auxiliary_matrix
      ! Option path:
      character(len=*), intent(in):: option_path
      ! state, and inner_mesh are used to setup mg preconditioner of inner solve
      type(state_type), intent(in):: state
      type(mesh_type), intent(in):: inner_mesh

      type(vector_field), pointer :: positions, mesh_positions
      type(petsc_csr_matrix), pointer:: rotation_matrix

      ! Additional arrays used internally:
      ! Additional numbering types:
      type(petsc_numbering_type) :: petsc_numbering_u
      type(petsc_numbering_type) :: petsc_numbering_aux

      integer, dimension(:), allocatable:: ghost_nodes, ghost_nodes_aux
      ! stuff for "mg" preconditioner
      type(petsc_csr_matrix), dimension(:), pointer:: prolongators
      integer, dimension(:), pointer:: surface_nodes

      PetscObject:: myPETSC_NULL_OBJECT
      PetscErrorCode ierr
      KSP ksp_schur ! ksp object for the inner solve in the Schur Complement
      Mat G_t_comp ! PETSc compressible divergence matrix
      Mat G_t_incomp ! PETSc incompressible divergence matrix
      Mat G ! PETSc Gradient matrix (transpose of G_t_incomp)
      Mat S ! PETSc Stabilization matrix (auxiliary_matrix)
      Mat pmat ! PETSc preconditioning matrix
      
      character(len=OPTION_PATH_LEN) :: inner_option_path, inner_solver_option_path
      
      integer reference_node, stat, i, rotation_stat
      logical parallel, have_auxiliary_matrix, have_preconditioner_matrix

      logical :: apply_reference_node, apply_reference_node_from_coordinates, reference_node_owned

      ! Sort option paths etc...
      solver_option_path=complete_solver_option_path(option_path)
      inner_option_path= trim(option_path)//&
              "/prognostic/scheme/use_projection_method&
              &/full_schur_complement/inner_matrix[0]"

      if (have_option(trim(option_path)//'/name')) then
         call get_option(trim(option_path)//'/name', name)
         ewrite(1,*) 'Inside petsc_solve_(block_)csr, solving for: ', trim(name)
      else
         ewrite(1,*) 'Inside petsc_solve_(block_)csr, solving using option_path: ', trim(option_path)
         name=option_path
      end if

      ! Are we applying a reference pressure node?
      apply_reference_node = have_option(trim(option_path)//&
                 '/prognostic/reference_node')
      apply_reference_node_from_coordinates = have_option(trim(option_path)//&
                 '/prognostic/reference_coordinates')

      ! If so, impose reference pressure node:
      if(apply_reference_node) then

         call get_option(trim(option_path)//&
              '/prognostic/reference_node', reference_node)
         if (GetProcNo()==1) then
            ewrite(2,*) 'Imposing_reference_pressure_node'        
            allocate(ghost_nodes(1:1))
            ghost_nodes(1) = reference_node
            call set(rhs,reference_node,0.0) ! Modify RHS accordingly
         else
            allocate(ghost_nodes(1:0))
         end if

      elseif(apply_reference_node_from_coordinates) then

         ewrite(1,*) 'Imposing_reference_pressure_node from user-specified coordinates'
         positions => extract_vector_field(state, "Coordinate")
         call find_reference_node_from_coordinates(positions,rhs%mesh,option_path,reference_node,reference_node_owned)
         if(IsParallel()) then
            if (reference_node_owned) then
               allocate(ghost_nodes(1:1))
               ghost_nodes(1) = reference_node
               call set(rhs,reference_node,0.0)
            else
               allocate(ghost_nodes(1:0))
            end if
         else
            allocate(ghost_nodes(1:1))
            ghost_nodes(1) = reference_node
            call set(rhs,reference_node,0.0)
         end if

      else

         allocate(ghost_nodes(1:0))

      end if

      ! Is auxiliary matrix present?
      have_auxiliary_matrix = present(auxiliary_matrix)
      ewrite(2,*) 'Are we using an auxiliary matrix: ', have_auxiliary_matrix

      ! set up numbering used in PETSc objects:
      call allocate(petsc_numbering_u, &
           nnodes=block_size(div_matrix_comp,2), nfields=blocks(div_matrix_comp,2), &
           halo=div_matrix_comp%sparsity%column_halo)
      call allocate(petsc_numbering_p, &
           nnodes=block_size(div_matrix_comp,1), nfields=1, &
           halo=preconditioner_matrix%sparsity%row_halo, ghost_nodes=ghost_nodes)

           ! - why is this using the row halo of the preconditioner matrix when there might be rows missing?
           ! - same question about the nnodes use of the rows of the block of the divergence matrix?
           ! - and how can ghost_nodes be appropriate for both this and the auxiliary_matrix?
           ! this definitely appears to be inappropriate for the auxiliary matrix (hence there's a new one added
           ! below) so two questions:
           ! 1. is it definitely appropriate for all its other used (the divergence matrix and the pressure vectors)?
           ! 2. can it be made appropriate for the auxiliary matrix at the same time as being appropriate for the current uses?

      ! Convert Divergence matrix (currently stored as block_csr matrix) to petsc format:   
      ! Create PETSc Div Matrix (comp & incomp) using this numbering:
      G_t_comp=block_csr2petsc(div_matrix_comp, petsc_numbering_p, petsc_numbering_u)
      G_t_incomp=block_csr2petsc(div_matrix_incomp, petsc_numbering_p, petsc_numbering_u)

      ! Scale G_t_comp to fit PETSc sign convention:
      call MatScale(G_t_comp,real(-1.0, kind = PetscScalar_kind),ierr)

      ! Determine transpose of G_t_incomp to form Gradient Matrix (G):
      call MatCreateTranspose(G_t_incomp,G,ierr)
      call MatSetOption(G, MAT_USE_INODES, PETSC_FALSE, ierr)

      ! Convert Stabilization matrix --> PETSc format if required:
      if(have_auxiliary_matrix) then
         ! set up numbering used in PETSc objects:
         ! NOTE: we use size(auxiliary_matrix,2) here as halo rows may be absent
         allocate(ghost_nodes_aux(1:0))
         call allocate(petsc_numbering_aux, &
              nnodes=size(auxiliary_matrix,2), nfields=1, &
              halo=auxiliary_matrix%sparsity%column_halo, &
              ghost_nodes=ghost_nodes_aux)
         
         S=csr2petsc(auxiliary_matrix, petsc_numbering_aux, petsc_numbering_aux)
         
         call deallocate(petsc_numbering_aux)
         deallocate(ghost_nodes_aux)
      end if
      
      ! Need to assemble the petsc matrix before we use it:
      call assemble(inner_M)

      ! Build Schur complement:
      ewrite(2,*) 'Building Schur complement'                
      if(have_auxiliary_matrix) then
         call MatCreateSchurComplement(inner_M%M,inner_M%M,G,G_t_comp,S,A,ierr)
      else
         myPETSC_NULL_OBJECT=PETSC_NULL_OBJECT
         call MatCreateSchurComplement(inner_M%M,inner_M%M,G,G_t_comp,PETSC_NULL_OBJECT,A,ierr)
         if (myPETSC_NULL_OBJECT/=PETSC_NULL_OBJECT) then
           FLAbort("PETSC_NULL_OBJECT has changed please report to skramer")
         end if
      end if
      
      ! Set ksp for M block solver inside the Schur Complement (the inner, inner solve!). 
      call MatSchurComplementGetKSP(A,ksp_schur,ierr)
      call petsc_solve_state_setup(inner_solver_option_path, prolongators, surface_nodes, &
        state, inner_mesh, blocks(div_matrix_comp,2), inner_option_path, matrix_has_solver_cache=.false., &
        mesh_positions=mesh_positions)
      rotation_matrix => extract_petsc_csr_matrix(state, "RotationMatrix", stat=rotation_stat)
      if (associated(prolongators)) then
        if (rotation_stat==0) then
          FLExit("Rotated boundary conditions do not work with mg prolongators")
        end if

        if (associated(surface_nodes)) then
          FLExit("Internal smoothing not available for inner solve")
        end if
        if (associated(mesh_positions)) then
          call setup_ksp_from_options(ksp_schur, inner_M%M, inner_M%M, &
            inner_solver_option_path, petsc_numbering=petsc_numbering_u, startfromzero_in=.true., &
            prolongators=prolongators, positions=mesh_positions)
        else
          call setup_ksp_from_options(ksp_schur, inner_M%M, inner_M%M, &
            inner_solver_option_path, petsc_numbering=petsc_numbering_u, startfromzero_in=.true., &
            prolongators=prolongators)
        end if
        do i=1, size(prolongators)
          call deallocate(prolongators(i))
        end do
        deallocate(prolongators)
      else if (associated(mesh_positions) .and. rotation_stat==0) then
        call setup_ksp_from_options(ksp_schur, inner_M%M, inner_M%M, &
          inner_solver_option_path, petsc_numbering=petsc_numbering_u, startfromzero_in=.true., &
          positions=mesh_positions, rotation_matrix=rotation_matrix%M)
      else if (associated(mesh_positions)) then
        call setup_ksp_from_options(ksp_schur, inner_M%M, inner_M%M, &
          inner_solver_option_path, petsc_numbering=petsc_numbering_u, startfromzero_in=.true., &
          positions=mesh_positions)
      else if (rotation_stat==0) then
        call setup_ksp_from_options(ksp_schur, inner_M%M, inner_M%M, &
          inner_solver_option_path, petsc_numbering=petsc_numbering_u, startfromzero_in=.true., &
          rotation_matrix=rotation_matrix%M)
      else
        call setup_ksp_from_options(ksp_schur, inner_M%M, inner_M%M, &
          inner_solver_option_path, petsc_numbering=petsc_numbering_u, startfromzero_in=.true.)
      end if
      
      if (associated(mesh_positions)) then
        call deallocate(mesh_positions)
        deallocate(mesh_positions)
      end if
       
      ! leaving out petsc_numbering and mesh, so "iteration_vtus" monitor won't work!

      ! Assemble preconditioner matrix in petsc format (if required):
      have_preconditioner_matrix=.not.(have_option(trim(option_path)//&
              "/prognostic/scheme/use_projection_method&
              &/full_schur_complement/preconditioner_matrix::NoPreconditionerMatrix"))

      if(have_preconditioner_matrix) then
         pmat=csr2petsc(preconditioner_matrix, petsc_numbering_p, petsc_numbering_p)
      end if

      ! Set up RHS and Solution vectors (note these are loaded later):
      b = PetscNumberingCreateVec(petsc_numbering_p)
      call VecDuplicate(b,y,ierr) ! Duplicate vector b and form vector y

      ! Schur complement objects now fully set up. The next stage is to setup
      ! KSP and PC for the solution of Ay = b, where A = G_t*M^-1*G.

      parallel=IsParallel()

      if(have_preconditioner_matrix) then
         call SetupKSP(ksp,A,pmat,solver_option_path,parallel,petsc_numbering_p, lstartfromzero)
      else
         ! If preconditioner matrix is not required, send in A instead:
         call SetupKSP(ksp,A,A,solver_option_path,parallel,petsc_numbering_p, lstartfromzero)
      end if
      
      ! Destroy the matrices setup for the schur complement computation. While
      ! these matrices are destroyed here, they are still required for the inner solve,
      ! since the schur complement depends on them. As a result, even though they're destroyed
      ! here, the PETSc reference count ensures that they are not wiped from memory.
      ! In other words, references to these arrays remain intact. 

      call MatDestroy(G_t_comp,ierr) ! Destroy Compressible Divergence Operator.
      call MatDestroy(G_t_incomp, ierr) ! Destroy Incompressible Divergence Operator.
      call MatDestroy(G, ierr) ! Destroy Gradient Operator (i.e. transpose of incompressible div).
      if(have_preconditioner_matrix) call MatDestroy(pmat,ierr) ! Destroy preconditioning matrix.
      if(have_auxiliary_matrix) call MatDestroy(S,ierr) ! Destroy stabilization matrix
      
      call deallocate( petsc_numbering_u )
      ! petsc_numbering_p is passed back and destroyed there      

    end subroutine petsc_solve_setup_full_projection

  end module Full_Projection
  
