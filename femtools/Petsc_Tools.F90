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
module Petsc_Tools
  use FLDebug
  use Sparse_Tools
  use parallel_tools
  use fields_base
  use fields_manipulation
  use halo_data_types
  use halos_base
  use halos_communications
  use halos_numbering
  use Reference_Counting
  use profiler
#ifdef HAVE_PETSC_MODULES
  use petsc 
#endif
  implicit none

#include "petsc_legacy.h"

  PetscReal, parameter, private :: dummy_petsc_real = 0.0
  integer, parameter, public :: PetscReal_kind = kind(dummy_petsc_real)
  PetscScalar, parameter, private :: dummy_petsc_scalar = 0.0
  integer, parameter, public :: PetscScalar_kind = kind(dummy_petsc_scalar)
  
  ! see at the bottom of this file
  interface
    subroutine myMatGetInfo(A, flag, info, ierr)
       Mat, intent(in):: A
       MatInfoType, intent(in):: flag
       double precision, dimension(:), intent(out):: info
       PetscErrorCode, intent(out):: ierr
    end subroutine myMatGetInfo
  end interface

  type petsc_numbering_type
     type(halo_type), pointer :: halo => null()
     integer nprivatenodes
     ! length of a vector 
     integer universal_length
     integer offset
     ! mapping between "global" (fludity numbering inside each local domain)
     ! and "universal" numbering (truly global numbering over all processes
     ! used by PETSc), second index is for multi-component fields
     integer, dimension(:,:), pointer:: gnn2unn
     ! list of ghost nodes, these are skipped in copying from and to 
     ! PETSc vectors, and will have a zero row and column in the matrix
     ! with something suitable on the diagonal
     integer, dimension(:), pointer:: ghost_nodes => null()
     ! the universal numbering of the ghost_nodes
     ! (in gnn2unn the ghost_nodes are masked out with -1)
     integer, dimension(:,:), pointer:: ghost2unn => null()
     !! Reference counting
     type(refcount_type), pointer :: refcount => null()
     character(len=FIELD_NAME_LEN):: name=""
  end type petsc_numbering_type
  
  interface allocate
    module procedure allocate_petsc_numbering
  end interface
  
  interface deallocate
    module procedure deallocate_petsc_numbering
  end interface
    
  interface field2petsc
    module procedure VectorFields2Petsc, ScalarFields2Petsc, VectorField2Petsc, ScalarField2Petsc, ScalarFieldPtrs2Petsc
  end interface
    
  interface petsc2field
    module procedure Petsc2VectorFields, Petsc2ScalarFields, Petsc2VectorField, Petsc2ScalarField, Petsc2ScalarFieldPtrs
  end interface

  interface petsc_numbering_create_is
    module procedure petsc_numbering_create_is_dim
  end interface
    
#include "Reference_count_interface_petsc_numbering_type.F90"

  private

  public reorder, DumpMatrixEquation, Initialize_Petsc
  public csr2petsc, petsc2csr, block_csr2petsc, petsc2array, array2petsc
  public field2petsc, petsc2field, petsc_numbering_create_is
  public petsc_numbering_type, PetscNumberingCreateVec, allocate, deallocate
  public csr2petsc_CreateSeqAIJ, csr2petsc_CreateMPIAIJ
  public addup_global_assembly
  ! for petsc_numbering:
  public incref, decref, addref
  ! for unit-testing:
  logical, public, save :: petsc_test_error_handler_called = .false.
  public petsc_test_error_handler
#if PETSC_VERSION_MINOR>=3
  public MatCreateSeqAIJ, MatCreateMPIAIJ, MatCreateSeqBAIJ, MatCreateMPIBAIJ
#endif
#if PETSC_VERSION_MINOR<5
  public mykspgetoperators
#endif
contains

  ! Note about definitions in this module:
  !
  ! In this module the real arrays are assumed to only store 1 
  ! value per node continuously in memory. So if multiple fields, or fields
  ! with more than one component (vector fields) are stored in the arrays
  ! they start with nonods values for the first component of the first field
  ! followed by blocks of size nonods for each of the other components and
  ! fields. These blocks are refered to either as block or field (so one 
  ! component of one vector field counts as one field). The 
  ! (block_)csr_matrices are layed out correspondingly (using the word
  ! block in a similar way).
  
  ! In the produced petsc_numbering we allow for (some of) the components of a
  ! field to be grouped together per node in memory. Such a group forms
  ! a local block of values per node in memory, but will always be referred to
  ! as group to avoid confusion with the above definition.
  
  subroutine allocate_petsc_numbering(petsc_numbering, &
       nnodes, nfields, halo, ghost_nodes)
    !!< Set ups the 'universal'(what most people call global)
    !!< numbering used in PETSc. In serial this is trivial
    !!< but could still be used for reordering schemes.
    !! the numbering object created:
    type(petsc_numbering_type), intent(out):: petsc_numbering
    !! number of nodes and fields:
    !! (here nfields counts each scalar component of vector fields, so
    !!  e.g. for nphases velocity fields in 3 dimensions nfields=3*nphases)
    integer, intent(in):: nnodes, nfields
    !! for parallel: halo information
    type(halo_type), pointer, optional :: halo
    !! If supplied number these as -1, so they'll be skipped by Petsc
    integer, dimension(:), optional, intent(in):: ghost_nodes 
    integer, dimension(:), allocatable:: ghost_marker
    integer i, g, f, start, offset
    integer nuniversalnodes, ngroups, lgroup_size, ierr

    allocate( petsc_numbering%gnn2unn(1:nnodes, 1:nfields) )

    if (present(halo)) then
       if (associated(halo)) then

          allocate(petsc_numbering%halo)
          petsc_numbering%halo=halo
          call incref(petsc_numbering%halo)
       
       end if
    end if

    ngroups=nfields

    ! first we set up the petsc numbering for the first entry of each group only:

    if (.not.associated(petsc_numbering%halo)) then

       ! *** Serial case *or* parallel without halo

       ! standard, trivial numbering, starting at 0:
       start=0 ! start of each field -1
       do g=1, nfields
          petsc_numbering%gnn2unn(:, g )= &
               (/ ( start+i, i=0, nnodes-1 ) /)
          start=start+nnodes
       end do

       if (isParallel()) then

          ! universal numbering can now be worked out trivially
          ! by calculating the offset (start of the universal number
          ! range for each process)
          call mpi_scan(nnodes, offset, 1, MPI_INTEGER, &
               MPI_SUM, MPI_COMM_FEMTOOLS, ierr)
          offset=offset-nnodes
          petsc_numbering%gnn2unn=petsc_numbering%gnn2unn+offset

       end if
       
       petsc_numbering%nprivatenodes=nnodes

       ! the offset is the first universal number assigned to this process
       ! in the standard petsc numbering the universal number is equal to
       ! offset+local number
       petsc_numbering%offset=petsc_numbering%gnn2unn(1,1)

    else

       ! *** Parallel case with halo:

       ! get 'universal' numbering
       call get_universal_numbering(halo, petsc_numbering%gnn2unn)
       ! petsc uses base 0
       petsc_numbering%gnn2unn = petsc_numbering%gnn2unn-1
         
       petsc_numbering%nprivatenodes=halo_nowned_nodes(halo)

       petsc_numbering%offset=halo%my_owned_nodes_unn_base*nfields             

    end if

    if (isParallel()) then
       ! work out the length of global(universal) vector
       call mpi_allreduce(petsc_numbering%nprivatenodes, nuniversalnodes, 1, MPI_INTEGER, &
           MPI_SUM, MPI_COMM_FEMTOOLS, ierr)       

       petsc_numbering%universal_length=nuniversalnodes*nfields
    else
       ! trivial in serial case:
       petsc_numbering%universal_length=nnodes*nfields
    end if

    if (present(ghost_nodes)) then
       if (associated(petsc_numbering%halo)) then
          ! check whether any of the halo nodes have been marked as
          ! ghost nodes by the owner of that node
          allocate( ghost_marker( 1:nnodes ) )
          ghost_marker = 0
          ghost_marker( ghost_nodes ) = 1
          call halo_update( petsc_numbering%halo, ghost_marker )

          ! fill in ghost_nodes list, now including halo nodes
          g=count(ghost_marker/=0)
          allocate( petsc_numbering%ghost_nodes(1:g), &
             petsc_numbering%ghost2unn(1:g, 1:nfields) )
          g=0
          do i=1, nnodes
            if (ghost_marker(i)/=0) then
              g=g+1
              petsc_numbering%ghost_nodes(g)=i
              ! store the original universal number seperately
              petsc_numbering%ghost2unn(g,:)=petsc_numbering%gnn2unn(i,:)
              ! mask it out with -1 in gnn2unn
              petsc_numbering%gnn2unn(i,:)=-1
            end if
          end do
          assert(g == size(petsc_numbering%ghost_nodes))
       else
          ! serial case, or no halo in parallel
          g=size(ghost_nodes)
          allocate( petsc_numbering%ghost_nodes(1:g), &
             petsc_numbering%ghost2unn(1:g, 1:nfields) )
          petsc_numbering%ghost_nodes=ghost_nodes
          do g=1, size(ghost_nodes)
             i=ghost_nodes(g)
             ! store the original universal number seperately
             petsc_numbering%ghost2unn(g,:)=petsc_numbering%gnn2unn(i,:)
             ! mask it out with -1 in gnn2unn
             petsc_numbering%gnn2unn(i,:)=-1
          end do
       end if
       
    else
       nullify( petsc_numbering%ghost_nodes )
       nullify( petsc_numbering%ghost2unn )
    end if
    
    nullify(petsc_numbering%refcount) ! Hack for gfortran component initialisation
    !                         bug.
    call addref(petsc_numbering)

  end subroutine allocate_petsc_numbering
  
  subroutine deallocate_petsc_numbering(petsc_numbering)
  !!< Deallocate the petsc_numbering object
  type(petsc_numbering_type), intent(inout):: petsc_numbering
  
    call decref(petsc_numbering)
    if (has_references(petsc_numbering)) return
    
    deallocate(petsc_numbering%gnn2unn)
    if (associated(petsc_numbering%ghost_nodes)) then
      deallocate(petsc_numbering%ghost_nodes)
      deallocate(petsc_numbering%ghost2unn)
    end if
    
    if (associated(petsc_numbering%halo)) then
      call deallocate(petsc_numbering%halo)
      deallocate(petsc_numbering%halo)
    end if
    
  end subroutine deallocate_petsc_numbering
  
  subroutine reorder(petsc_numbering, sparsity, ordering_type)
  type(petsc_numbering_type), intent(inout):: petsc_numbering
  type(csr_sparsity), intent(in):: sparsity
  MatOrderingType, intent(in):: ordering_type
    
    Mat M
    IS rperm, cperm  
    PetscErrorCode ierr
    PetscInt, dimension(:), allocatable:: iperm
    integer nnodes, nfields, nprivatenodes
    integer b, g, ghost
    
    nnodes=size(petsc_numbering%gnn2unn,1)
    nfields=size(petsc_numbering%gnn2unn,2)
    nprivatenodes=petsc_numbering%nprivatenodes
    
    ! sets up sequential matrix with only private nodes,
    ! ignoring all halo entries:
    M=CreatePrivateMatrixFromSparsity(sparsity)
    
    call MatGetOrdering(M, ordering_type, rperm, cperm, ierr)
    
    call MatDestroy(M, ierr)

    allocate(iperm(1:nprivatenodes))
    call ISCopyIndices(rperm, iperm, ierr)
    iperm=iperm+1
    
    if (associated(petsc_numbering%ghost_nodes)) then
      ! fill ghost nodes back in (overwriting the -1s)
      ! so they get reordered as well
      do g=1, size(petsc_numbering%ghost_nodes)
        ghost=petsc_numbering%ghost_nodes(g)
        petsc_numbering%gnn2unn(ghost,:)=petsc_numbering%ghost2unn(g,:)
      end do
    end if
    
    do b=1, nfields
      petsc_numbering%gnn2unn(iperm, b)= &
        petsc_numbering%gnn2unn(1:nprivatenodes, b)
    end do
      
    if (associated(petsc_numbering%ghost_nodes)) then
      ! put the -1s back in and write new ghost2unn
      do g=1, size(petsc_numbering%ghost_nodes)
        ghost=petsc_numbering%ghost_nodes(g)
        petsc_numbering%ghost2unn(g,:)=petsc_numbering%gnn2unn(ghost,:)
        petsc_numbering%gnn2unn(ghost,:)=-1
      end do
    end if
    deallocate(iperm)
    
    call ISDestroy(rperm, ierr)
    call ISDestroy(cperm, ierr)
      
  end subroutine reorder

  subroutine Array2Petsc(array, petsc_numbering, vec)
    !!< Assembles petsc array using the specified numbering.
    !!< Allocates a petsc Vec that should be destroyed with VecDestroy
    real, dimension(:), intent(in):: array
    type(petsc_numbering_type), intent(in):: petsc_numbering
    Vec, intent(inout) :: vec
  
    integer ierr, nnodp, start, b, nfields, nnodes

    ! number of nodes owned by this process:
    nnodp=petsc_numbering%nprivatenodes
    
    ! number of nodes on this process, including halo/ghost nodes
    nnodes=size(petsc_numbering%gnn2unn, 1)
    
    ! number of fields:
    nfields=size(petsc_numbering%gnn2unn, 2)
    
    ! start of each field -1
    start=0

    if (associated(petsc_numbering%halo)) then
       if (.not. ((petsc_numbering%halo%data_type .eq. HALO_TYPE_CG_NODE) &
            .or. (petsc_numbering%halo%data_type .eq. HALO_TYPE_DG_NODE))) then
          FLAbort("Matrices can only be assembled on the basis of node halos")
       end if
    end if

    do b=1, nfields
      
#ifdef DOUBLEP
      call VecSetValues(vec, nnodp, &
        petsc_numbering%gnn2unn( 1:nnodp, b ), &
        array( start+1:start+nnodp ), INSERT_VALUES, ierr)
#else
      call VecSetValues(vec, nnodp, &
        petsc_numbering%gnn2unn( 1:nnodp, b ), &
        real(array( start+1:start+nnodp ), kind = PetscScalar_kind), INSERT_VALUES, ierr)
#endif
        
      ! go to next field:
      start=start+nnodes
    
    end do
    
    call VecAssemblyBegin(vec, ierr)
    call VecAssemblyEnd(vec, ierr)
    
  end subroutine Array2Petsc
  
  subroutine VectorFields2Petsc(fields, petsc_numbering, vec)
    !!< Assembles contiguous petsc array using the specified numbering from the given fields.
    !!< Allocates a petsc Vec that should be destroyed with VecDestroy
    type(vector_field), dimension(:), intent(in):: fields
    type(petsc_numbering_type), intent(in):: petsc_numbering
    Vec :: vec
  
    integer ierr, nnodp, b, nfields, nnodes
    integer i, j

    ! number of nodes owned by this process:
    nnodp=petsc_numbering%nprivatenodes
    
    ! number of nodes on this process, including halo/ghost nodes
    nnodes=size(petsc_numbering%gnn2unn, 1)
    
    ! number of "component" fields: (i.e. n/o vector fields *times* n/o components per vector field)
    nfields=size(petsc_numbering%gnn2unn, 2)
    assert( nfields==sum(fields%dim) )
    
    if (associated(petsc_numbering%halo)) then
       if (.not. ((petsc_numbering%halo%data_type .eq. HALO_TYPE_CG_NODE) &
            .or. (petsc_numbering%halo%data_type .eq. HALO_TYPE_DG_NODE))) then
          FLAbort("Matrices can only be assembled on the basis of node halos")
       end if
    end if

    b=1
    do i=1, size(fields)
      
       do j=1, fields(i)%dim
#ifdef DOUBLEP
         call VecSetValues(vec, nnodp, &
            petsc_numbering%gnn2unn( 1:nnodp, b ), &
            fields(i)%val(j, 1:nnodp ), INSERT_VALUES, ierr)
#else
         call VecSetValues(vec, nnodp, &
            petsc_numbering%gnn2unn( 1:nnodp, b ), &
            real(fields(i)%val(j, 1:nnodp ), kind = PetscScalar_kind), INSERT_VALUES, ierr)
#endif
        b=b+1
        
      end do
    
    end do
    
    call VecAssemblyBegin(vec, ierr)
    call VecAssemblyEnd(vec, ierr)
    
  end subroutine VectorFields2Petsc
  
  subroutine VectorField2Petsc(field, petsc_numbering, vec)
    !!< Assembles contiguous petsc array using the specified numbering from the given field.
    type(vector_field), intent(in):: field
    type(petsc_numbering_type), intent(in):: petsc_numbering
    Vec, intent(out) :: vec
  
    call VectorFields2Petsc( (/ field /), petsc_numbering, vec)
    
  end subroutine VectorField2Petsc
  
  subroutine ScalarFields2Petsc(fields, petsc_numbering, vec)
    !!< Assembles contiguous petsc array using the specified numbering from the given fields.
    type(scalar_field), dimension(:), intent(in):: fields
    type(petsc_numbering_type), intent(in):: petsc_numbering
    Vec, intent(inout) :: vec
  
    integer ierr, nnodp, b, nfields, nnodes

    ! number of nodes owned by this process:
    nnodp=petsc_numbering%nprivatenodes
    
    ! number of nodes on this process, including halo/ghost nodes
    nnodes=size(petsc_numbering%gnn2unn, 1)
    
    ! number of fields:
    nfields=size(petsc_numbering%gnn2unn, 2)
    assert(nfields==size(fields))

    if (associated(petsc_numbering%halo)) then
       if (.not. ((petsc_numbering%halo%data_type .eq. HALO_TYPE_CG_NODE) &
            .or. (petsc_numbering%halo%data_type .eq. HALO_TYPE_DG_NODE))) then
          FLAbort("Matrices can only be assembled on the basis of node halos")
       end if
    end if

    do b=1, nfields
      
#ifdef DOUBLEP
       call VecSetValues(vec, nnodp, &
            petsc_numbering%gnn2unn( 1:nnodp, b ), &
            fields(b)%val( 1:nnodp ), INSERT_VALUES, ierr)
#else
       call VecSetValues(vec, nnodp, &
            petsc_numbering%gnn2unn( 1:nnodp, b ), &
            real(fields(b)%val( 1:nnodp ), kind = PetscScalar_kind), INSERT_VALUES, ierr)
#endif
            
    end do
    
    call VecAssemblyBegin(vec, ierr)
    call VecAssemblyEnd(vec, ierr)

  end subroutine ScalarFields2Petsc
  
  subroutine ScalarField2Petsc(field, petsc_numbering, vec)
    !!< Assembles petsc array using the specified numbering.
    !!< Allocates a petsc Vec that should be destroyed with VecDestroy
    type(scalar_field), intent(in):: field
    type(petsc_numbering_type), intent(in):: petsc_numbering
    Vec, intent(inout) :: vec
  
    call ScalarFields2Petsc( (/ field /), petsc_numbering, vec)
  
  end subroutine ScalarField2Petsc

  subroutine ScalarFieldPtrs2Petsc(fields, petsc_numbering, vec)
    !!< Assembles contiguous petsc array using the specified numbering from the given field pointers.
    type(scalar_field_pointer), dimension(:), intent(in):: fields
    type(petsc_numbering_type), intent(in):: petsc_numbering
    Vec, intent(inout) :: vec
  
    integer ierr, nnodp, b, nfields, nnodes

    ! number of nodes owned by this process:
    nnodp=petsc_numbering%nprivatenodes
    
    ! number of nodes on this process, including halo/ghost nodes
    nnodes=size(petsc_numbering%gnn2unn, 1)
    
    ! number of fields:
    nfields=size(petsc_numbering%gnn2unn, 2)
    assert(nfields==size(fields))

    if (associated(petsc_numbering%halo)) then
       if (.not. ((petsc_numbering%halo%data_type .eq. HALO_TYPE_CG_NODE) &
            .or. (petsc_numbering%halo%data_type .eq. HALO_TYPE_DG_NODE))) then
          FLAbort("Matrices can only be assembled on the basis of node halos")
       end if
    end if

    do b=1, nfields
      
#ifdef DOUBLEP
       call VecSetValues(vec, nnodp, &
            petsc_numbering%gnn2unn( 1:nnodp, b ), &
            fields(b)%ptr%val( 1:nnodp ), INSERT_VALUES, ierr)
#else
       call VecSetValues(vec, nnodp, &
            petsc_numbering%gnn2unn( 1:nnodp, b ), &
            real(fields(b)%ptr%val( 1:nnodp ), kind = PetscScalar_kind), INSERT_VALUES, ierr)
#endif
            
    end do
    
    call VecAssemblyBegin(vec, ierr)
    call VecAssemblyEnd(vec, ierr)

  end subroutine ScalarFieldPtrs2Petsc
  
  function PetscNumberingCreateVec(petsc_numbering) result (vec)
    !!< Creates a petsc array with size corresponding to petsc_numbering.
    !!< After use it should be destroyed with VecDestroy. No vector values
    !!< are set in this function.
    type(petsc_numbering_type), intent(in):: petsc_numbering
    Vec vec
  
    
    integer ierr, nnodp, plength, ulength, nfields
    logical parallel

    ! number of nodes owned by this process:
    nnodp=petsc_numbering%nprivatenodes
    
    ! number of fields:
    nfields=size(petsc_numbering%gnn2unn, 2)
    
    ! length of local (private) vector
    plength=nnodp*nfields
    
    ! length of global (universal) vector
    ulength=petsc_numbering%universal_length
    
    ! whether this is a parallel vector:
    parallel= (associated(petsc_numbering%halo))
    
    if (parallel) then
       call VecCreateMPI(MPI_COMM_FEMTOOLS, plength, ulength, vec, ierr)
    else
      call VecCreateSeq(MPI_COMM_SELF, ulength, vec, ierr)
    end if

    call VecSetOption(vec, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE, ierr)
    call VecSetOption(vec, VEC_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)

  end function PetscNumberingCreateVec

  function petsc_numbering_create_is_dim(petsc_numbering, dim) result (index_set)
    IS:: index_set
    type(petsc_numbering_type), intent(in):: petsc_numbering
    integer, intent(in):: dim

    PetscErrorCode:: ierr
    integer:: nnodp

    nnodp = petsc_numbering%nprivatenodes

#if PETSC_VERSION_MINOR>=2
    call ISCreateGeneral(MPI_COMM_FEMTOOLS, nnodp, petsc_numbering%gnn2unn(:,dim), &
         PETSC_COPY_VALUES, index_set, ierr)
#else
    call ISCreateGeneral(MPI_COMM_FEMTOOLS, nnodp, petsc_numbering%gnn2unn(:,dim), &
         index_set, ierr)
#endif
       
  end function petsc_numbering_create_is_dim
  
  subroutine Petsc2Array(vec, petsc_numbering, array)
    !!< Copies the values of a PETSc Vec into an array. The PETSc Vec
    !!< must have been assembled using the same petsc_numbering.
    Vec, intent(in):: vec
    type(petsc_numbering_type), intent(in):: petsc_numbering
    real, dimension(:), intent(out) :: array
    
    integer ierr, nnodp, start, b, nfields, nnodes
#ifndef DOUBLEP
    PetscScalar, dimension(:), allocatable :: vals
#endif
    
    ! number of nodes owned by this process:
    nnodp=petsc_numbering%nprivatenodes
    
    ! number of nodes on this process, including halo/ghost nodes
    nnodes=size(petsc_numbering%gnn2unn, 1)
    
    ! number of fields:
    nfields=size(petsc_numbering%gnn2unn, 2)
    
    ! start of each field -1
    start=0
    
#ifdef DOUBLEP
    do b=1, nfields
      
      ! this check should be unnecessary but is a work around for a bug in petsc, fixed in 18ae1927 (pops up with intel 15)
      if (nnodp>0) then
        call VecGetValues(vec, nnodp, &
          petsc_numbering%gnn2unn( 1:nnodp, b ), &
          array( start+1:start+nnodp ), ierr)
      end if
        
      ! go to next field:
      start=start+nnodes

      ! Update the halo:
      if (associated(petsc_numbering%halo)) then
         call halo_update(petsc_numbering%halo, &
              array( start+1:start+nnodp ))
      end if
    
    end do
#else
    allocate(vals(nnodp))
    do b=1, nfields
      
      if (nnodp>0) then
        call VecGetValues(vec, nnodp, &
          petsc_numbering%gnn2unn( 1:nnodp, b ), &
          vals, ierr)
        array( start+1:start+nnodp ) = vals
      end if
        
      ! go to next field:
      start=start+nnodes

      ! Update the halo:
      if (associated(petsc_numbering%halo)) then
         call halo_update(petsc_numbering%halo, &
              array( start+1:start+nnodp ))
      end if
    
    end do
    deallocate(vals)
#endif

  end subroutine Petsc2Array

  subroutine Petsc2ScalarFields(vec, petsc_numbering, fields, rhs)
    !!< Copies the values of a PETSc Vec into scalar fields. The PETSc Vec
    !!< must have been assembled using the same petsc_numbering.
    Vec, intent(in):: vec
    type(petsc_numbering_type), intent(in):: petsc_numbering
    type(scalar_field), dimension(:), intent(inout) :: fields
    !! for ghost_nodes the value of the rhs gets copied into fields
    type(scalar_field), dimension(:), intent(in), optional :: rhs
    
    integer ierr, nnodp, b, nfields, nnodes
#ifndef DOUBLEP
    PetscScalar, dimension(:), allocatable :: vals
#endif
    
    ! number of nodes owned by this process:
    nnodp=petsc_numbering%nprivatenodes
    
    ! number of nodes on this process, including halo/ghost nodes
    nnodes=size(petsc_numbering%gnn2unn, 1)
    
    ! number of fields:
    nfields=size(petsc_numbering%gnn2unn, 2)
    
#ifdef DOUBLEP
    do b=1, nfields
      
      call profiler_tic(fields(b), "petsc2field")
      ! this check should be unnecessary but is a work around for a bug in petsc, fixed in 18ae1927 (pops up with intel 15)
      if (nnodp>0) then
        call VecGetValues(vec, nnodp, &
          petsc_numbering%gnn2unn( 1:nnodp, b ), &
          fields(b)%val( 1:nnodp ), ierr)
      end if
      call profiler_toc(fields(b), "petsc2field")
        
    end do
#else
    allocate(vals(nnodp))
    do b=1, nfields
      
      call profiler_tic(fields(b), "petsc2field")
      if (nnodp>0) then
        call VecGetValues(vec, nnodp, &
          petsc_numbering%gnn2unn( 1:nnodp, b ), &
          vals, ierr)
      end if
      fields(b)%val( 1:nnodp ) = vals
      call profiler_toc(fields(b), "petsc2field")
        
    end do
    deallocate(vals)
#endif
      
    if (associated(petsc_numbering%ghost_nodes) .and. present(rhs)) then
       
       do b=1, nfields
         call profiler_tic(fields(b), "petsc2field")
         call set(fields(b), petsc_numbering%ghost_nodes, &
            node_val(rhs(b), petsc_numbering%ghost_nodes))
         call profiler_toc(fields(b), "petsc2field")
       end do
       
    end if        
      

    if (associated(petsc_numbering%halo)) then
       if (size(fields)>1) then
          ewrite(2, *) 'Updating of halo of multiple scalar fields needs to be improved.'
       end if
       do b=1, size(fields)
          call profiler_tic(fields(b), "petsc2field")
          ! Always update on the level 2 halo to ensure that the whole
          ! field is well defined.
          call halo_update(fields(b), 2)
          call profiler_toc(fields(b), "petsc2field")
       end do
    end if

  end subroutine Petsc2ScalarFields
    
  subroutine Petsc2ScalarField(vec, petsc_numbering, field, rhs)
  !!< Copies the values of a PETSc Vec into a scalar field. The PETSc Vec
  !!< must have been assembled using the same petsc_numbering.
  Vec, intent(in):: vec
  type(petsc_numbering_type), intent(in):: petsc_numbering
  type(scalar_field), intent(inout) :: field
  !! for ghost_nodes the value of the rhs gets copied into fields
  type(scalar_field), intent(in), optional :: rhs
    
      type(scalar_field) fields(1), rhss(1)
        
      fields(1)=field
      if (present(rhs)) then
        rhss(1)=rhs
        call Petsc2ScalarFields(vec, petsc_numbering, fields, rhs=rhss)
      else
        call Petsc2ScalarFields(vec, petsc_numbering, fields)
      end if
    
  end subroutine Petsc2ScalarField
    
  subroutine Petsc2VectorFields(vec, petsc_numbering, fields)
  !!< Copies the values of a PETSc Vec into vector fields. The PETSc Vec
  !!< must have been assembled using the same petsc_numbering.
  Vec, intent(in):: vec
  type(petsc_numbering_type), intent(in):: petsc_numbering
  type(vector_field), dimension(:), intent(inout) :: fields
    
    integer ierr, nnodp, b, nfields, nnodes
    integer i, j
#ifndef DOUBLEP
    PetscScalar, dimension(:), allocatable :: vals
#endif
    
    ! number of nodes owned by this process:
    nnodp=petsc_numbering%nprivatenodes
    
    ! number of nodes on this process, including halo/ghost nodes
    nnodes=size(petsc_numbering%gnn2unn, 1)
    
    ! number of "component" fields: (i.e. n/o vector fields *times* n/o components per vector field)
    nfields=size(petsc_numbering%gnn2unn, 2)
    assert( nfields==sum(fields%dim) )
    
#ifdef DOUBLEP
    b=1
    do i=1, size(fields)
      
      call profiler_tic(fields(i), "petsc2field")
      do j=1, fields(i)%dim
         ! this check should be unnecessary but is a work around for a bug in petsc, fixed in 18ae1927 (pops up with intel 15)
         if (nnodp>0) then
           call VecGetValues(vec, nnodp, &
             petsc_numbering%gnn2unn( 1:nnodp, b ), &
             fields(i)%val(j, 1:nnodp ), ierr)
         end if
         b=b+1
      end do
      call profiler_toc(fields(i), "petsc2field")
        
    end do
#else
    allocate(vals(nnodp))
    b=1
    do i=1, size(fields)
      
      call profiler_tic(fields(i), "petsc2field")
      do j=1, fields(i)%dim
         if (nnodp>0) then
           call VecGetValues(vec, nnodp, &
             petsc_numbering%gnn2unn( 1:nnodp, b ), &
             vals, ierr)
         end if
         fields(i)%val(j, 1:nnodp ) = vals
         b=b+1
      end do
      call profiler_toc(fields(i), "petsc2field")
        
    end do
    deallocate(vals)
#endif

    ! Update the halo:
    if (associated(petsc_numbering%halo)) then
       ewrite(2, *) '*** Updating of halo of vector_fields needs to be improved.'
       do i=1, size(fields)
          ! Always update on the level 2 halo to ensure that the whole
          ! field is well defined.
          call profiler_tic(fields(i), "petsc2field")
          call halo_update(fields(i), 2)
          call profiler_toc(fields(i), "petsc2field")
       end do
    end if

  end subroutine Petsc2VectorFields
  
  subroutine Petsc2VectorField(vec, petsc_numbering, field)
  !!< Copies the values of a PETSc Vec into a vector field. The PETSc Vec
  !!< must have been assembled using the same petsc_numbering.
  Vec, intent(in):: vec
  type(petsc_numbering_type), intent(in):: petsc_numbering
  type(vector_field), intent(inout) :: field
    
      type(vector_field) fields(1)
        
      fields(1)=field
      call Petsc2VectorFields(vec, petsc_numbering, fields)
    
  end subroutine Petsc2VectorField

  subroutine Petsc2ScalarFieldPtrs(vec, petsc_numbering, fields)
    !!< Copies the values of a PETSc Vec into scalar fields. The PETSc Vec
    !!< must have been assembled using the same petsc_numbering.
    Vec, intent(in):: vec
    type(petsc_numbering_type), intent(in):: petsc_numbering
    type(scalar_field_pointer), dimension(:), intent(inout) :: fields
    
    integer ierr, nnodp, b, nfields, nnodes
#ifndef DOUBLEP
    PetscScalar, dimension(:), allocatable :: vals
#endif
    
    ! number of nodes owned by this process:
    nnodp=petsc_numbering%nprivatenodes
    
    ! number of nodes on this process, including halo/ghost nodes
    nnodes=size(petsc_numbering%gnn2unn, 1)
    
    ! number of fields:
    nfields=size(petsc_numbering%gnn2unn, 2)
    
#ifdef DOUBLEP
    do b=1, nfields
      
      call profiler_tic(fields(b)%ptr, "petsc2field")
      call VecGetValues(vec, nnodp, &
        petsc_numbering%gnn2unn( 1:nnodp, b ), &
        fields(b)%ptr%val( 1:nnodp ), ierr)
      call profiler_toc(fields(b)%ptr, "petsc2field")
        
    end do
#else
    allocate(vals(nnodp))
    do b=1, nfields
      
      call profiler_tic(fields(b)%ptr, "petsc2field")
      call VecGetValues(vec, nnodp, &
        petsc_numbering%gnn2unn( 1:nnodp, b ), &
        vals, ierr)
      fields(b)%ptr%val( 1:nnodp ) = vals
      call profiler_toc(fields(b)%ptr, "petsc2field")
        
    end do
    deallocate(vals)
#endif            

    if (associated(petsc_numbering%halo)) then
       if (size(fields)>1) then
          ewrite(2, *) 'Updating of halo of multiple scalar field pointers needs to be improved.'
       end if
       do b=1, size(fields)
          call profiler_tic(fields(b)%ptr, "petsc2field")
          ! Always update on the level 2 halo to ensure that the whole
          ! field is well defined.
          call halo_update(fields(b)%ptr, 2)
          call profiler_toc(fields(b)%ptr, "petsc2field")
       end do
    end if

  end subroutine Petsc2ScalarFieldPtrs
      
  function csr2petsc(A, petsc_numbering, column_petsc_numbering, &
       use_inodes) result(M)
  !!< Converts a csr_matrix from Sparse_Tools into a PETSc matrix.
  !!< Note: this function creates a PETSc matrix, it has to be deallocated
  !!< with MatDestroy by the user.
  type(csr_matrix), intent(in):: A
  !! If present use the following numbering, otherwise the standard numbering is
  !! set up and deallocated again.
  type(petsc_numbering_type), optional, intent(in):: petsc_numbering
  !! If present this numbering is used for the column numbering, and the previous
  !! for row numbering, otherwise the previous defines row and column numbering.
  !! In parallel they must be the same
  type(petsc_numbering_type), optional, intent(in):: column_petsc_numbering
  !! petsc's inodes don't work with certain preconditioners ("mg" and "eisenstat")
  !! that's why we default to not use them
  logical, intent(in), optional:: use_inodes
  Mat M
    
    type(block_csr_matrix) block_matrix
  
    block_matrix=wrap(A%sparsity, (/ 1, 1 /), A%val, name="TemporaryMatrix_csr2petsc")
    
    M=block_csr2petsc(block_matrix, petsc_numbering=petsc_numbering, &
        column_petsc_numbering=column_petsc_numbering, &
        use_inodes=use_inodes)
        
    call deallocate(block_matrix)
    
  end function csr2petsc
  
  function block_csr2petsc(A, petsc_numbering, column_petsc_numbering, &
       use_inodes) result(M)
  !!< Converts a block_csr_matrix from Sparse_Tools into a PETSc matrix.
  !!< Note: this function creates a PETSc matrix, it has to be deallocated
  !!< with MatDestroy by the user. 
  type(block_csr_matrix), intent(in):: A
  !! If present use the following numbering, otherwise the standard numbering is
  !! set up and deallocated again.
  type(petsc_numbering_type), optional, intent(in):: petsc_numbering
  !! If present this numbering is used for the column numbering, and the previous
  !! for row numbering, otherwise the previous defines row and column numbering.
  !! In parallel they must be the same.
  type(petsc_numbering_type), optional, intent(in):: column_petsc_numbering
  !! petsc's inodes don't work with certain preconditioners ("mg" and "eisenstat")
  !! that's why we default to not use them
  logical, intent(in), optional:: use_inodes
  Mat M
    
    type(petsc_numbering_type) row_numbering, col_numbering
    real, dimension(:), pointer:: vals
    real ghost_pivot, mindiag, maxdiag, diag
    integer, dimension(:), pointer:: cols 
    integer, dimension(:), allocatable:: colidx
    integer, dimension(:), allocatable:: row2ghost
    integer nbrows, nbcols, nblocksv, nblocksh
    integer nbrowsp, nbcolsp
    integer rows(1)
    integer len, bh, bv, i, g, row, ierr
   
    if (present(petsc_numbering)) then
      row_numbering=petsc_numbering
    else
       ! set up standard numbering for the rows:
       if (associated(A%sparsity%row_halo)) then
          call allocate(row_numbering, &
               nnodes=size(A, 1), nfields=blocks(A, 1), &
               halo=A%sparsity%row_halo)
       else
          call allocate(row_numbering, &
               nnodes=size(A, 1), nfields=blocks(A, 1))
       end if
    end if
    
    if (present(column_petsc_numbering)) then
      col_numbering=column_petsc_numbering
    else
      ! set up standard numbering for the columns:
       if (associated(A%sparsity%column_halo)) then
          call allocate(col_numbering, &
               nnodes=size(A, 2), nfields=blocks(A, 2), &
               halo=A%sparsity%column_halo)
       else
          call allocate(col_numbering, &
               nnodes=size(A, 2), nfields=blocks(A, 2))
       end if
    end if
    
    ! rows and cols per block:
    nbrows=size(row_numbering%gnn2unn, 1)
    nbcols=size(col_numbering%gnn2unn, 1)
    ! number of vertical and horizontal blocks:
    nblocksv=size(row_numbering%gnn2unn, 2)
    nblocksh=size(col_numbering%gnn2unn, 2)
    
    ! number of private rows and cols in each block
    nbrowsp=row_numbering%nprivatenodes
    nbcolsp=col_numbering%nprivatenodes
    
    ! setup reverse mapping from row no to ghost no
    allocate( row2ghost(1:nbrows) )
    row2ghost=0
    if (associated(row_numbering%ghost_nodes)) then
      ! only do something on the diagonal if the row numbering and
      ! column numbering have the same ghost nodes, otherwise the
      ! ghost rows and/or columns are just zeroed:
      if (size(row_numbering%ghost_nodes) > 0 .and. &
        associated(row_numbering%ghost_nodes, col_numbering%ghost_nodes) &
        ) then
       
         row2ghost( row_numbering%ghost_nodes )= &
            (/ ( i, i=1, size(row_numbering%ghost_nodes)) /)
            
         ! now find a suitable value to put on the diagonal
         mindiag=huge(0.0)
         maxdiag=-mindiag
         do i=1, nbrowsp
           do bv=1, nblocksv
              diag=abs(val(A, bv, bv, i, i))
              if (diag<mindiag) mindiag=diag
              if (diag>maxdiag) maxdiag=diag
           end do
         end do
         ghost_pivot=sqrt((maxdiag+mindiag)*maxdiag/2.0)
      end if
    end if
    
    ! collect the lengths of all rows 
    ! (the total horizontal row length across all blocks)
    if (.not. IsParallel()) then

      ! Create serial matrix:
      M=csr2petsc_CreateSeqAIJ(A%sparsity, row_numbering, col_numbering, A%diagonal, use_inodes=use_inodes)
      
      ! these should be the same, just to make sure:
      nbrowsp=nbrows
      
    else
    
      ! Create parallel matrix:
      M=csr2petsc_CreateMPIAIJ(A%sparsity, row_numbering, col_numbering, A%diagonal, use_inodes=use_inodes)

      call MatSetOption(M, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)
    endif 

    allocate(colidx(1:nbcols))
    
    if (associated(row_numbering%halo)) then
       if (.not. ((row_numbering%halo%data_type .eq. HALO_TYPE_CG_NODE) &
            .or. (row_numbering%halo%data_type .eq. HALO_TYPE_DG_NODE))) then
          FLAbort("Matrices can only be assembled on the basis of node halos")
       end if
    end if

    ! loop over rows within a block:
    do i=1, nbrowsp
      
      if (row2ghost(i)==0) then
         cols => row_m_ptr(A, i)
         len=size(cols)
        
         ! loop through all rows i in the different blocks:
         ! outer loop: from left to right
         do bh=1, nblocksh
        
            ! translate column indices to petsc
            colidx(1:len)=col_numbering%gnn2unn(cols, bh)
          
            ! we go down first as all the column indices stay the same
            do bv=1, nblocksv
               if (A%diagonal .and. bh/=bv) cycle
               ! row number in PETSc land:
               rows(1)=row_numbering%gnn2unn(i, bv)
               vals => row_val_ptr(A, bv, bh, i)
#ifdef DOUBLEP
               call MatSetValues(M, 1, rows, len, colidx(1:len), vals, &
                   INSERT_VALUES, ierr)
#else
               call MatSetValues(M, 1, rows, len, colidx(1:len), real(vals, kind = PetscScalar_kind), &
                   INSERT_VALUES, ierr)
#endif
            end do

         end do
           
        ! only set ghost pivot for owned rows:
      else if (i<=nbrowsp) then
      
         g=row2ghost(i)
         do bv=1, nblocksv
            row=row_numbering%ghost2unn(g, bv)
#ifdef DOUBLEP
            call MatSetValue(M, row, row, ghost_pivot, INSERT_VALUES, ierr)
#else
            call MatSetValue(M, row, row, real(ghost_pivot, kind = PetscScalar_kind), INSERT_VALUES, ierr)
#endif
         end do
          
      end if
      
    end do
    
    deallocate(colidx)
    
    call MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY, ierr)
    
    if (.not. present(petsc_numbering)) then
      call deallocate(row_numbering)
    endif
    if (.not. present(column_petsc_numbering)) then
      call deallocate(col_numbering)
    endif

  end function block_csr2petsc
  
  function CreatePrivateMatrixFromSparsity(sparsity) result (M)
  ! creates Petsc matrix containing only entries corresponding to private nodes
  type(csr_sparsity), intent(in):: sparsity
  
    Mat M
    PetscErrorCode ierr
    real, dimension(:), allocatable:: vals
    integer, dimension(:), pointer:: cols
    integer rows(1)
    integer, dimension(:), allocatable:: colidx, nnz
    integer ncols, nprows, npcols
    integer i, l
    
    ncols=size(sparsity,2)
    if (associated(sparsity%row_halo)) then
       nprows=halo_nowned_nodes(sparsity%row_halo)
    else
       nprows=size(sparsity,1)
    end if
    if (associated(sparsity%column_halo)) then
       npcols=halo_nowned_nodes(sparsity%column_halo)
    else
       npcols=size(sparsity,2)
    end if

    if (.not. IsParallel()) then
      nprows=size(sparsity,1)
      npcols=ncols
    end if
    
    allocate(nnz(1:nprows))
    ! calcute n/o nonzero entries of private rows:
    do i=1, nprows
      cols => row_m_ptr(sparsity, i)
      ! ignore halo column indices
      nnz(i)=count(cols<=npcols)

    end do
    
    call MatCreateSeqAIJ(MPI_COMM_SELF, nprows, npcols, &
      PETSC_NULL_INTEGER, nnz, M, ierr)
      
    call MatSetOption(M, MAT_USE_INODES, PETSC_FALSE, ierr)


    deallocate(nnz)
    
    allocate(colidx(1:ncols), vals(1:ncols))
    ! some random number for the matrix values:
    vals=1.0
      
    do i=1, nprows
      cols => row_m_ptr(sparsity,i)
      l=size(cols)
      where (cols<=npcols)
        colidx(1:l)=cols-1
      elsewhere
        ! -1 is skipped over by PETSc
        colidx(1:l)=-1
      end where
      rows(1)=i-1
#ifdef DOUBLEP
      call MatSetValues(M, 1, rows, l, colidx(1:l), vals(1:l), &
            INSERT_VALUES, ierr)
#else
      call MatSetValues(M, 1, rows, l, colidx(1:l), real(vals(1:l), kind = PetscScalar_kind), &
            INSERT_VALUES, ierr)
#endif
    end do    
    
    call MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY, ierr)
    
    deallocate(colidx, vals)

  end function CreatePrivateMatrixFromSparsity

  function csr2petsc_CreateSeqAIJ(sparsity, row_numbering, col_numbering, only_diagonal_blocks, use_inodes) result(M)
  !!< Creates a sequential PETSc Mat of size corresponding with
  !!< row_numbering and col_numbering.
  type(csr_sparsity), intent(in):: sparsity
  type(petsc_numbering_type), intent(in):: row_numbering, col_numbering
  logical, intent(in):: only_diagonal_blocks
  !! petsc's inodes don't work with certain preconditioners ("mg" and "eisenstat")
  !! that's why we default to not use them
  logical, intent(in), optional:: use_inodes
  Mat M

    integer, dimension(:), allocatable:: nnz
    integer nrows, ncols, nbrows, nbcols, nblocksv, nblocksh
    integer row, len, ierr
    integer bv, i

    ! total number of rows and cols:
    nrows=row_numbering%universal_length
    ncols=col_numbering%universal_length
    ! rows and cols per block:
    nbrows=size(row_numbering%gnn2unn, 1)
    nbcols=size(col_numbering%gnn2unn, 1)
    ! number of vertical and horizontal blocks:
    nblocksv=size(row_numbering%gnn2unn, 2)
    nblocksh=size(col_numbering%gnn2unn, 2)

    allocate(nnz(0:nrows-1))
    ! loop over complete horizontal rows within a block of rows
    nnz=1 ! ghost rows are skipped below and only have a diagonal
    do i=1, nbrows
      if (only_diagonal_blocks) then
        len=row_length(sparsity,i)
      else
        len=row_length(sparsity,i)*nblocksh
      end if
      ! loop over the blocks of rows
      do bv=1, nblocksv
        ! row in petsc numbering:
        row=row_numbering%gnn2unn(i,bv)
        if (row/=-1) then
          nnz(row)=len
        end if
      end do
    end do
      
    call MatCreateSeqAIJ(MPI_COMM_SELF, nrows, ncols, PETSC_NULL_INTEGER, &
      nnz, M, ierr)
      
    if (.not. present_and_true(use_inodes)) then
      call MatSetOption(M, MAT_USE_INODES, PETSC_FALSE, ierr)
    end if

    deallocate(nnz)
      
  end function csr2petsc_CreateSeqAIJ
  
  function csr2petsc_CreateMPIAIJ(sparsity, row_numbering, col_numbering, only_diagonal_blocks, use_inodes) result(M)
  !!< Creates a parallel PETSc Mat of size corresponding with
  !!< row_numbering and col_numbering.
  type(csr_sparsity), intent(in):: sparsity
  type(petsc_numbering_type), intent(in):: row_numbering, col_numbering
  logical, intent(in):: only_diagonal_blocks
  !! petsc's inodes don't work with certain preconditioners ("mg" and "eisenstat")
  !! that's why we default to not use them
  logical, intent(in), optional:: use_inodes
  Mat M
  
    integer, dimension(:), pointer:: cols 
    integer, dimension(:), allocatable:: d_nnz, o_nnz
    integer nrows, ncols, nbrows, nbcols, nblocksv, nblocksh
    integer nrowsp, ncolsp, nbrowsp, nbcolsp, row_offset
    integer row, len, private_len, ghost_len
    integer bv, i, ierr

    ! total number of rows and cols:
    nrows=row_numbering%universal_length
    ncols=col_numbering%universal_length
    ! rows and cols per block:
    nbrows=size(row_numbering%gnn2unn, 1)
    nbcols=size(col_numbering%gnn2unn, 1)
    ! number of vertical and horizontal blocks:
    nblocksv=size(row_numbering%gnn2unn, 2)
    nblocksh=size(col_numbering%gnn2unn, 2)
    
    ! number of private rows and cols in each block
    nbrowsp=row_numbering%nprivatenodes
    nbcolsp=col_numbering%nprivatenodes
    
    ! number of private rows and cols in total
    ! (this will be the number of local rows and cols for petsc)
    nrowsp=nbrowsp*nblocksv
    ncolsp=nbcolsp*nblocksh

    ! the universal numbers used by petsc for private nodes are in the range
    ! row_offset:row_offset+nrowsp-1
    row_offset=row_numbering%offset
    
    ! for each private row we have to count the number of column indices 
    ! refering to private nodes and refering to ghost/halos nodes
    allocate(d_nnz(row_offset:row_offset+nrowsp-1), &
       o_nnz(row_offset:row_offset+nrowsp-1))
    ! ghost rows are skipped below, and only have a diagonal
    d_nnz=1
    o_nnz=0
    ! loop over complete horizontal rows within a block of rows
    do i=1, nbrowsp
      ! this is just the row within a block
      cols => row_m_ptr(sparsity, i)
      ! the row length over all blocks from left to right
      len=size(cols)
      ! number of entries refering to private nodes:
      private_len=count(cols<=nbcolsp)
      if (.not. only_diagonal_blocks) then
        len=len*nblocksh
        private_len=private_len*nblocksh
      end if
      ! the rest refers to ghost/halo nodes:
      ghost_len=len-private_len
      
      ! loop over the blocks of rows
      do bv=1, nblocksv
        ! row in petsc numbering:
        row=row_numbering%gnn2unn(i,bv)
        if (row/=-1) then
           ASSERT(row>=row_offset .and. row<row_offset+nrowsp)
           d_nnz(row)=private_len
           o_nnz(row)=ghost_len
        end if
      end do
    end do

    call MatCreateMPIAIJ(MPI_COMM_FEMTOOLS, nrowsp, ncolsp, nrows, ncols, &
      PETSC_NULL_INTEGER, d_nnz, PETSC_NULL_INTEGER, o_nnz, M, ierr)
      
    if (.not. present_and_true(use_inodes)) then
      call MatSetOption(M, MAT_USE_INODES, PETSC_FALSE, ierr)
    end if

    deallocate(d_nnz, o_nnz)

  end function csr2petsc_CreateMPIAIJ

  function petsc2csr(matrix, column_numbering) result(A)
  !!< Converts a PETSc matrix into a csr_matrix from Sparse_Tools.
  !!< Note: this function allocates a csr_matrix, it has to be deallocated
  !!< by the user.
  !!< Note2: for parallel matrices this assumes standard numbering, i.e. the
  !!< numbering where each processor owns a consecutive range of 'nrows' local rows
  !!< from offset+0 to offset+nrows-1. If column_numbering is not provided
  !!< the column numbering is exactly the same as for rows.
  !!< Note3: for parallel matrices where column_numbering is not provided, 
  !!< only the entries for which both the row number *and* the column number 
  !!< are local are copied over, i.e. "remote
  !!< entries" in a local row are left out.
  !!< Note4: for parallel matrices where column_numbering *is* provided
  !!< only local rows are copied over (i.e. no halo *rows*). The translation
  !!< back from global to local indices for the halo columns is done via
  !!< a naive inverse map, i.e. we allocate an integer array of global length!
  !!< This should therefore not be used in production code, unless some more
  !!< (memory) efficient mapping is implemented - currently used for 
  !!< petsc_readnsolve only.
  type(csr_matrix) :: A
  Mat, intent(in):: matrix
  type(petsc_numbering_type), optional, intent(in):: column_numbering

    PetscErrorCode ierr
    type(csr_sparsity) :: sparsity
    double precision, dimension(MAT_INFO_SIZE):: matrixinfo
    PetscScalar, dimension(:), allocatable:: row_vals
    integer, dimension(:), allocatable:: row_cols, unn2gnn
    integer private_columns
    integer i, j, k, ui, rows, columns, entries, ncols, offset
    logical parallel
    
    ! get the necessary info about the matrix:
    call myMatGetInfo(matrix, MAT_LOCAL, matrixinfo, ierr)
    entries=matrixinfo(MAT_INFO_NZ_USED)
    ! note we're no longer using MAT_INFO for getting local n/o rows and cols
    ! as it's bugged in Petsc < 3.0 and obsoloted thereafter:
    call MatGetLocalSize(matrix, rows, columns, ierr)
    call MatGetOwnershipRange(matrix, offset, PETSC_NULL_INTEGER, ierr)
    parallel=IsParallel()

    if (present(column_numbering)) then
      ! only halo columns are copied over:
      private_columns=columns
      columns=size(column_numbering%gnn2unn,1)
      assert( private_columns==column_numbering%nprivatenodes )
      assert( private_columns<=columns )
      assert( size(column_numbering%gnn2unn,2)==1 )
      call allocate(sparsity, rows, columns, entries, &
        diag=.false., name="petsc2csrSparsity")
      if (associated(column_numbering%halo)) then
         allocate(sparsity%column_halo)
         sparsity%column_halo=column_numbering%halo
         call incref(sparsity%column_halo)
      end if
      ! allocate inverse mapping
      ! WARNING: this is a 'universal' length array
      allocate( unn2gnn(0:column_numbering%universal_length-1) )
      unn2gnn=0
      do i=1, columns
        ui=column_numbering%gnn2unn(i,1)
        if (ui>=0) then
          ! unused halo nodes are marked -1 (see allocate_petsc_numbering() )
          unn2gnn( ui ) = i
        end if
      end do
    else
      ! in parallel no halo rows or columns  are included, so
      ! in both the serial and the parallel case we allocate a matrix with
      ! nprivate_nodes==nnodes:
      call allocate(sparsity, rows, columns, entries, diag=.false., &
        name="petsc2csrSparsity")
    end if
    call allocate(A, sparsity)
    
    if (parallel) then
      ! in the parallel case we first copy in a temp. buffer
      ! and only insert the entries A_ij with both i *and* j local
      allocate(row_cols(1:entries), row_vals(1:entries))
      j=1
      do i=0, rows-1
        sparsity%findrm(i+1)=j
        call MatGetRow(matrix, offset+i, ncols, row_cols, row_vals, ierr)
        do k=1, ncols
          if (row_cols(k)>=offset .and. row_cols(k)<offset+rows) then
            ! only local->local entries get inserted
            sparsity%colm(j)=row_cols(k)-offset+1
            A%val(j)=row_vals(k)
            j=j+1
          else if (present(column_numbering)) then
            ! halo columns
            sparsity%colm(j)=unn2gnn(row_cols(k))
            A%val(j)=row_vals(k)
            assert( sparsity%colm(j)>0 )
            j=j+1
          end if
        end do
        ! This is stupid, we were given copies in MatGetRow so it could
        ! have restored its internal tmp arrays straight away, anyway:
        call MatRestoreRow(matrix, offset+i, ncols, row_cols, row_vals, ierr)
      end do
      A%sparsity%findrm(i+1)=j

      deallocate(row_cols, row_vals)
    else
      ! Serial case:
      j=1
      do i=0, rows-1
        sparsity%findrm(i+1)=j
#ifdef DOUBLEP
        call MatGetRow(matrix, offset+i, ncols, sparsity%colm(j:), A%val(j:), ierr)
        j=j+ncols
        ! This is stupid, we were given copies in MatGetRow so it could
        ! have restored its internal tmp arrays straight away, anyway:
        call MatRestoreRow(matrix, offset+i, ncols, sparsity%colm(j:), A%val(j:), ierr)
#else        
        allocate(row_vals(size(A%val) - j + 1))
        call MatGetRow(matrix, offset+i, ncols, sparsity%colm(j:), row_vals, ierr)
        A%val(j:) = row_vals
        j=j+ncols
        ! This is stupid, we were given copies in MatGetRow so it could
        ! have restored its internal tmp arrays straight away, anyway:
        call MatRestoreRow(matrix, offset+i, ncols, sparsity%colm(j:), row_vals, ierr)
        deallocate(row_vals)
#endif
      end do
      A%sparsity%findrm(i+1)=j
      
      ! matrix is indexed from offset+0, colm should be indexed from 1:
      sparsity%colm=sparsity%colm-offset+1

    end if
    call deallocate(sparsity)
    
    if (present(column_numbering)) then
      deallocate(unn2gnn)
    end if

  end function petsc2csr
  
  subroutine addup_global_assembly(vfield, halo)
  !!< adds up the local contributions in a globally assembled vfield
  !!< i.e. the non-owned contributions in halo nodes get added into
  !!< the non-halo nodes on the owning processes. This is followed by
  !!< a halo_update so that the added up values are communicated back to
  !!< the halo nodes
  type(vector_field), intent(inout):: vfield
  type(halo_type), pointer:: halo
  
    type(petsc_numbering_type):: petsc_numbering
    Vec:: vec
    PetscErrorCode:: ierr
    
    if (.not. IsParallel()) return
    
    call allocate(petsc_numbering, node_count(vfield), vfield%dim, &
      halo)
    vec=PetscNumberingCreateVec(petsc_numbering)
    ! assemble vfield into petsc Vec, this lets petsc do the adding up
    call field2petsc(vfield, petsc_numbering, vec)
    ! copy back (this includes the promised halo update):
    call petsc2field(vec, petsc_numbering, vfield)
    call VecDestroy(vec, ierr)
    call deallocate(petsc_numbering)
    
  end subroutine addup_global_assembly
    
  
  function FindrmFromRowSizes(sizes)
  !!< Auxilary routine to work out findrm from the row sizes
  integer, dimension(:), intent(in):: sizes
  integer, dimension(1:size(sizes)+1):: FindrmFromRowSizes

    integer i, j
    
    j=1
    do i=1, size(sizes)
      FindrmFromRowSizes(i)=j
      j=j+sizes(i)
    end do
    FindrmFromRowSizes(i)=j

  end function FindrmFromRowSizes

  subroutine DumpMatrixEquation(filename, x0, A, b)
  character(len=*), intent(in):: filename
  Mat, intent(in):: A
  Vec, intent(in):: x0, b
  
     PetscViewer :: viewer
     PetscErrorCode :: ierr
     
     ewrite(0, *) 'Dumping matrix equation in file called '//filename
     call PetscViewerBinaryOpen(MPI_COMM_FEMTOOLS, &
          filename, FILE_MODE_WRITE, &
          viewer, ierr)
     call MatView(A, viewer, ierr)
     call VecView(b, viewer, ierr)
     call VecView(x0, viewer, ierr)
     call PetscViewerDestroy(viewer, ierr)

  end subroutine DumpMatrixEquation

  subroutine Initialize_Petsc()
    PetscErrorCode :: ierr
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr);
  end subroutine Initialize_Petsc

! Simple dummy error handler that just tracks whether it's been called or not
! Useful for unittesting to see that petsc gives error messages at the right moment
#if PETSC_VERSION_MINOR>=2
subroutine petsc_test_error_handler(comm,line, func, file, dir, n, p, mess, ctx, ierr)
  MPI_Comm:: comm
#else
subroutine petsc_test_error_handler(line, func, file, dir, n, p, mess, ctx, ierr)
#endif
  PetscInt:: line
  character(len=*):: func, file, dir
  PetscErrorCode:: n
  PetscInt:: p
  character(len=*):: mess
  PetscInt:: ctx
  PetscErrorCode:: ierr
  

  petsc_test_error_handler_called = .true.
  
end subroutine petsc_test_error_handler

! In petsc-3.3 the MatCreate[B]{Seq|MPI}() routines have changed to MatCreate[B]Aij
! and MatSetup always needs to be called
#if PETSC_VERSION_MINOR>=3
  subroutine MatCreateSeqAIJ(MPI_Comm, nrows, ncols, &
      nz, nnz, M, ierr)
    integer, intent(in):: MPI_Comm
    PetscInt, intent(in):: nrows, ncols, nz
    PetscInt, dimension(:), intent(in):: nnz
    Mat, intent(out):: M
    PetscErrorCode, intent(out):: ierr

    call MatCreateAij(MPI_Comm, nrows, ncols, nrows, ncols, &
      nz, nnz, 0, PETSC_NULL_INTEGER, M, ierr)
    call MatSetup(M, ierr)

  end subroutine MatCreateSeqAIJ

  subroutine MatCreateMPIAIJ(MPI_Comm, nprows, npcols, &
      nrows, ncols, &
      dnz, dnnz, onz, onnz, M, ierr)
    integer, intent(in):: MPI_Comm
    PetscInt, intent(in):: nprows, npcols,nrows, ncols, dnz, onz
    PetscInt, dimension(:), intent(in):: dnnz, onnz
    Mat, intent(out):: M
    PetscErrorCode, intent(out):: ierr

    call MatCreateAij(MPI_Comm, nprows, npcols, nrows, ncols, &
      dnz, dnnz, onz, onnz, M, ierr)
    call MatSetup(M, ierr)

  end subroutine MatCreateMPIAIJ

  subroutine MatCreateSeqBAIJ(MPI_Comm, bs, nrows, ncols, &
      nz, nnz, M, ierr)
    integer, intent(in):: MPI_Comm
    PetscInt, intent(in):: bs, nrows, ncols, nz
    PetscInt, dimension(:), intent(in):: nnz
    Mat, intent(out):: M
    PetscErrorCode, intent(out):: ierr

    call MatCreateBAij(MPI_Comm, bs, nrows, ncols, nrows, ncols, &
      nz, nnz, 0, PETSC_NULL_INTEGER, M, ierr)
    call MatSetup(M, ierr)

  end subroutine MatCreateSeqBAIJ

  subroutine MatCreateMPIBAIJ(MPI_Comm, bs, nprows, npcols, &
      nrows, ncols, &
      dnz, dnnz, onz, onnz, M, ierr)
    integer, intent(in):: MPI_Comm
    PetscInt, intent(in):: bs, nprows, npcols,nrows, ncols, dnz, onz
    PetscInt, dimension(:), intent(in):: dnnz, onnz
    Mat, intent(out):: M
    PetscErrorCode, intent(out):: ierr

    call MatCreateBAij(MPI_Comm, bs, nprows, npcols, nrows, ncols, &
      dnz, dnnz, onz, onnz, M, ierr)
    call MatSetup(M, ierr)

  end subroutine MatCreateMPIBAIJ
#endif

! this is a wrapper around KSPGetOperators, that in petsc <3.5
! has an extra mat_structure flag. We need to wrap this because
! we need a local variable.
! in include/petsc_legacy.h we #define KSPGetOperators -> mykspgetoperators
#if PETSC_VERSION_MINOR<5
subroutine mykspgetoperators(ksp, amat, pmat, ierr)
  KSP, intent(in):: ksp
  Mat, intent(in):: amat, pmat
  PetscErrorCode, intent(out):: ierr

  MatStructure:: mat_structure
  
  ! need small caps, to avoid #define from include/petsc_legacy.h
  call  kspgetoperators(ksp, amat, pmat, mat_structure, ierr)

end subroutine mykspgetoperators
#endif

#include "Reference_count_petsc_numbering_type.F90"
end module Petsc_Tools

! this is a wrapper routine around MatGetInfo as it appears PETSc
! provides the wrong explicit interface for it. By putting it outside
! the module (and only including petsc headers and not use petsc modules)
! this routine calls MatGetInfo with an implicit interface.
subroutine myMatGetInfo(A, flag, info, ierr)
Mat, intent(in):: A
MatInfoType, intent(in):: flag
double precision, dimension(:), intent(out):: info
PetscErrorCode, intent(out):: ierr

  call MatGetInfo(A, flag, info, ierr)
  
end subroutine myMatGetInfo

