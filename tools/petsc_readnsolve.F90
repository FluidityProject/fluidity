#include "confdefs.h"
#include "fdebug.h"
!! Little program that reads in a matrix equation from a file called 
!! 'matrixdump' in the current directory containing a matrix, rhs vector and 
!! initial guess, written in PETSc binary format. It then subsequently solves 
!! the corresponding system of equations with PETSc using PETSc options set
!! in the PETSC_OPTIONS environment variable.
!!
!! Extra functionality is activated using the following PETSc options:
!! -prns_filename <file>   reads from specified file instead of 'matrixdump'
!! -prns_zero_init_guess   no init.guess is read from matrixdump
!!                          instead the init.guess is zeroed
!! -prns_scipy             writes out several files that can be read in scipy.
!! -prns_random_rhs        randomize rhs of the equation
!! -prns_verbosity <n>   sets the verbosity level -1..3 (similar to fluidity)
!!
!! (*) for not too big matrices an 'exact' answer can be obtained with a direct
!!     solver: put PETSc options -ksp_type preonly and -pc_type lu. Although
!!     for ill-conditioned matrices its accuracy may be limited.
subroutine petsc_readnsolve
use quadrature
use elements
use fields
use halos
use state_module
use Multigrid
use FLDebug
use Sparse_Tools
use Petsc_Tools
use solvers
use sparse_tools_petsc
use petsc_solve_state_module
use Global_Parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
use spud
use populate_state_module
use field_options
use halos_registration
use parallel_tools
#ifdef HAVE_PETSC_MODULES
  use petsc 
#endif
implicit none
#include "petsc_legacy.h"
  ! options read from command-line (-prns_... options)
  character(len=4096) filename, flml
  character(len=FIELD_NAME_LEN):: field
  logical zero_init_guess, scipy, random_rhs

  PetscViewer viewer
  PetscRandom pr
  Mat matrix
  Vec rhs, x
  PetscErrorCode ierr

  call petsc_readnsolve_options(filename, flml, field, &
     zero_init_guess, &
     scipy, random_rhs)
    
  ewrite(1,*) 'Opening: ', trim(filename)
  
  ! read in the matrix equation and init. guess:
  call PetscViewerBinaryOpen(MPI_COMM_FEMTOOLS, trim(filename), &
     FILE_MODE_READ, viewer, ierr)
#if PETSC_VERSION_MINOR>=2
  ! in petsc 3.2 MatLoad and VecLoad do no longer create a new vector
  call MatCreate(MPI_COMM_FEMTOOLS, matrix, ierr)
  call VecCreate(MPI_COMM_FEMTOOLS, rhs, ierr)
  if (IsParallel()) then
    call MatSetType(matrix, MATMPIAIJ, ierr)
    call VecSetType(rhs, VECMPI, ierr)
  else
    call MatSetType(matrix, MATSEQAIJ, ierr)
    call VecSetType(rhs, VECSEQ, ierr)
  end if
  call MatLoad(matrix, viewer, ierr)
  call VecLoad(rhs, viewer, ierr)
  if (zero_init_guess) then
    call VecDuplicate(rhs, x, ierr)
  else 
    call VecCreate(MPI_COMM_FEMTOOLS, x, ierr)
    if (IsParallel()) then
      call VecSetType(x, VECMPI, ierr)
    else
      call VecSetType(x, VECSEQ, ierr)
    end if
    call VecLoad(x, viewer, ierr)
  end if
#else
  if (IsParallel()) then
    call MatLoad(viewer, MATMPIAIJ, matrix, ierr)
    call VecLoad(viewer, VECMPI, rhs, ierr)
  else
    call MatLoad(viewer, MATSEQAIJ, matrix, ierr)
    call VecLoad(viewer, VECSEQ, rhs, ierr)
  end if
  if (zero_init_guess) then
    call VecDuplicate(rhs, x, ierr)
  else if (IsParallel()) then
    call VecLoad(viewer, VECMPI, x, ierr)
  else
    call VecLoad(viewer, VECSEQ, x, ierr)
  end if
#endif
  call PetscViewerDestroy(viewer, ierr)
  if (random_rhs) then
    call PetscRandomCreate(PETSC_COMM_WORLD, pr, ierr)
    call VecSetRandom(rhs, PETSC_NULL_OBJECT, ierr)
    call PetscRandomDestroy(pr, ierr)
  end if

  
  if (flml=='') then
    call petsc_readnsolve_old_style(filename, &
      zero_init_guess, scipy, &
      matrix, x, rhs)
  else
    call petsc_readnsolve_flml(flml, field, &
       matrix, x, rhs)
  end if
    
contains

  subroutine petsc_readnsolve_old_style(filename, &
      zero_init_guess, scipy, &
      matrix, x, rhs)
  !! this is the old way of running petsc_readnsolve, everything
  !! is handled with PETSc calls locally, i.e. nothing from petsc_solve()
  !! in fluidity is used
  ! options read from command line:
  character(len=*), intent(in):: filename
  logical, intent(in):: zero_init_guess, scipy
  ! PETSc matrix, rhs vector and initial guess vector read from matrixdump:
  Mat, intent(inout):: matrix
  Vec, intent(inout):: x, rhs
    
    real, allocatable:: xv(:), rv(:), rhsv(:), dv(:)
    KSPType krylov_method
    PCType pc_method
    KSP krylov
    PC  prec
    Vec y
    PetscViewer viewer
    PetscBool flag
    PetscErrorCode ierr
    KSPConvergedReason reason
    type(petsc_numbering_type) petsc_numbering
    type(csr_matrix) A
    real value
    real time1, time2, time3
    integer i, n, n1, n2, iterations
    
    ewrite(0,*) "Called petsc_readnsolve without specifying flml."
    ewrite(0,*) "The recommend way of calling petsc_readnsolve is now"//&
                 &" by specifying both the .flml and field to solve for"//&
                 &" on the command line, e.g.:"
    ewrite(0,*) "   petsc_readnsolve mycase.flml Pressure"
    ewrite(0,*)
        
    ! output initial basic statistics:
    call VecNorm(rhs, NORM_2, real(value, kind = PetscScalar_kind), ierr)
    ewrite(2,*) 'Right-hand side 2-norm:', value
    call VecNorm(rhs, NORM_INFINITY, real(value, kind = PetscScalar_kind), ierr)
    ewrite(2,*) 'Right-hand side inf-norm:', value
    call VecNorm(x, NORM_2, real(value, kind = PetscScalar_kind), ierr)
    ewrite(2,*) 'init. guess 2-norm:', value
    call VecNorm(x, NORM_INFINITY, real(value, kind = PetscScalar_kind), ierr)
    ewrite(2,*) 'init. guess inf-norm:', value
    
    ! including inital residual:
    call VecDuplicate(x,y, ierr)
    call MatMult(matrix, x, y, ierr)
    call VecAXPY(y, real(-1.0, kind = PetscScalar_kind), rhs, ierr)
    call VecNorm(y, NORM_2, real(value, kind = PetscScalar_kind), ierr)
    ewrite(2,*) 'init. residual 2-norm:', value
    call VecNorm(y, NORM_INFINITY, real(value, kind = PetscScalar_kind), ierr)
    ewrite(2,*) 'init. residual inf-norm:', value  

    ! get values to locate 'large sping' boundary conditions:
    call VecGetOwnershipRange(rhs, n1, n2, ierr)
    n=n2-n1
    allocate(rv(1:n), xv(1:n), rhsv(1:n), dv(1:n))
    call VecGetValues(y, n, (/ (i, i=0, n-1) /)+n1, &
           real(rv, kind = PetscScalar_kind),ierr)
    call VecGetValues(x, n, (/ (i, i=0, n-1) /)+n1, &
           real(xv, kind = PetscScalar_kind),ierr)
    call VecGetValues(rhs, n, (/ (i, i=0, n-1) /)+n1, &
           real(rhsv, kind = PetscScalar_kind),ierr)
    call MatGetDiagonal(matrix, y, ierr)
    call VecGetValues(y, n, (/ (i, i=0, n-1) /)+n1, &
           real(dv, kind = PetscScalar_kind),ierr)
    
    ! This can get massive so lower verbosity if you don't want it:
    ewrite(3,*) 'Large spring boundary conditions found at:'
    ewrite(3,*) '    rownumber, init. guess, rhs and init. residual'
    do i=1, n
      if (abs(dv(i))>1d15) then
        ewrite(3, *) i, xv(i), rhsv(i), rv(i)
      end if
    end do
   
    ! set up PETSc solver:
    ! default values:
    krylov_method=KSPGMRES
    pc_method=PCSOR
    
    call PetscOptionsGetString('', "-ksp_type", krylov_method, flag, ierr)
    call PetscOptionsGetString('', "-pc_type", pc_method, flag, ierr)
    call KSPCreate(MPI_COMM_FEMTOOLS, krylov, ierr)
    call KSPSetType(krylov, krylov_method, ierr)
    call KSPSetOperators(krylov, matrix, matrix, ierr)
    call KSPSetTolerances(krylov, 1.0d-100, 1d-12, PETSC_DEFAULT_REAL, &
      3000, ierr)
    if (zero_init_guess) then
      call KSPSetInitialGuessNonzero(krylov, PETSC_FALSE, ierr)
      ! also explicitly zero x vector for output purposes
      call VecZeroEntries(x, ierr)
    else
      call KSPSetInitialGuessNonzero(krylov, PETSC_TRUE, ierr)
    end if

    ! first timer includes set up time:
    call cpu_time(time1)
    call KSPGetPC(krylov, prec, ierr)
    call PCSetType(prec, pc_method, ierr)
    if (pc_method==PCMG) then
      call SetupSmoothedAggregation(prec, matrix, ierr)
      if (ierr/=0) then
         ewrite(0,*) 'WARNING: setup of mg preconditioner failed'
         stop
      end if
    end if  
    call PCSetFromOptions(prec, ierr)
    call PCSetUp(prec, ierr)
    call KSPSetFromOptions(krylov, ierr)
    call KSPSetUp(krylov, ierr)

    ! pure solving:
    call cpu_time(time2)
    call KSPSolve(krylov, rhs, x, ierr)
    call cpu_time(time3)
    call KSPGetConvergedReason(krylov, reason, ierr)
    call KSPGetIterationNumber(krylov, iterations, ierr)
    
    ! this also kills the pc:
    call KSPDestroy(krylov, ierr)
    
    ! most basic solver statistics:
    ewrite(2,*) 'Convergence reason:', reason
    ewrite(2,*) 'after iterations:  ', iterations
    ewrite(2,*) 'Total time including set up:   ', time3-time1
    ewrite(2,*) 'Time in solving:   ', time3-time2
    
    ! output final residual:
    call MatMult(matrix, x, y, ierr)
    call VecAXPY(y, real(-1.0, kind = PetscScalar_kind), rhs, ierr)
    call VecNorm(y, NORM_2, real(value, kind = PetscScalar_kind), ierr)
    ewrite(2,*) 'Final residual 2-norm:', value
    call VecNorm(y, NORM_INFINITY, real(value, kind = PetscScalar_kind), ierr)
    ewrite(2,*) 'Final residual inf-norm:', value  

    ! write matrix and rhs in scipy readable format
    if (scipy) then
      
      A=petsc2csr(matrix)
      call mmwrite(trim(filename)//'.mm', A)
      call deallocate(A)
      
      ! RHS
      call PetscViewerASCIIOpen(MPI_COMM_FEMTOOLS, &
          trim(filename)//'.rhs.vec', &
          viewer, ierr)
      call VecView(rhs, viewer, ierr)
      call PetscViewerDestroy(viewer, ierr)

      ! solution
      call PetscViewerASCIIOpen(MPI_COMM_FEMTOOLS, &
          trim(filename)//'.sol.vec', &
          viewer, ierr)
      call VecView(x, viewer, ierr)
      call PetscViewerDestroy(viewer, ierr)
      
    end if
    
    call VecDestroy(y, ierr)
    call MatDestroy(matrix, ierr)
    call VecDestroy(x, ierr)
    call VecDestroy(rhs, ierr)

  end subroutine petsc_readnsolve_old_style
    
  subroutine petsc_readnsolve_flml(flml, field, &
      matrix, x, rhs)
  character(len=*), intent(in):: flml, field
  Mat, intent(inout):: matrix
  Vec, intent(inout):: x, rhs
  
    type(petsc_numbering_type):: petsc_numbering
    type(element_type):: shape
    type(quadrature_type):: quad
    type(mesh_type), pointer:: linear_mesh
    type(mesh_type):: mesh
    type(state_type), pointer:: states(:)
    type(petsc_csr_matrix):: A
    type(scalar_field):: x_field, rhs_field
    character(len=OPTION_PATH_LEN):: option_path
    character(len=FIELD_NAME_LEN):: field_name, mesh_name
    logical read_state, fail
    integer i, n, istate, stat
    type(halo_type), pointer :: my_halo
    integer nstates, universal_nodes, components, dim

    ewrite(1,*) "Opening flml file ", trim(flml)
    call load_options(flml)
    if (.not. have_option('/simulation_name')) then
      ! backtrace not very useful, here:
      FLExit("Failed to load options tree from flml.")
    end if
    
    option_path=workout_option_path(field, fail)
    if (fail) then
      
      ewrite(1,*) "Since I can't work out the option path for the specified field"
      ewrite(1,*) "I'll try populating the states first."
      read_state=.true.
      
    else
      
      read_state=petsc_solve_needs_state(option_path)
      if (read_state) then
        ewrite(1,*) "The specified solver options require geometry information"
        ewrite(1,*) "Will read in all state information first"
      end if
    
    end if
    
    if (read_state) then
      
      call populate_state(states)
      
      ! work out option_path (possibly again)
      call workout_option_path_from_state(states, field, option_path, istate)
      
    else if (IsParallel()) then
      
      ! we have an option_path, so we can
      ! short track by only reading the meshes:
      
      ! Find out how many states there are
      nstates=option_count("/material_phase")
      allocate(states(1:nstates))
      do i = 1, nstates
         call nullify(states(i))
      end do

      call insert_external_mesh(states)

      call insert_derived_meshes(states)
      
      ! for meshes all states are the same:
      istate=1
    
    end if
    
    ! now we have an option_path
    ! let's find or construct a mesh:
    call VecGetSize(x, n, ierr)
    
    
    if (read_state .or. IsParallel()) then
      
      call get_option(trim(complete_field_path(option_path))// &
        '/mesh[0]/name', mesh_name)
      mesh=extract_mesh(states(istate), mesh_name)
          
      ! now work out the number of nodes according to the mesh
      if (IsParallel()) then
      
        assert(halo_count(mesh) > 0)
        my_halo => mesh%halos(halo_count(mesh))
        call allocate(petsc_numbering, node_count(mesh), 1, my_halo)
      else
      
        call allocate(petsc_numbering, node_count(mesh), 1)

      end if
      
      universal_nodes=petsc_numbering%universal_length
      
      ! and compare it with the size of the PETSc vector
      if (universal_nodes==n) then
        
        ewrite(1,*) "Node count of mesh agrees with that of the matrixdump: ", n
        components=1
        
      else if (mod(n, universal_nodes)==0) then
        
        components=universal_nodes/n
        call petsc_readnsolve_vector(mesh, n, universal_nodes, option_path, x, matrix, rhs, &
            states(istate), read_state)
              
      else
      
        ewrite(1,*) "Number of nodes in specified mesh: ", universal_nodes
        ewrite(1,*) "Vector length in matrixdump:", n
        FLExit("Mesh and matrixdump size don't agree")
        
      end if
      
    else
    
      ! allocate a dummy quadrature, shape and mesh
      ! all we need is a mesh with node_count(mesh)=n
      call get_option('/geometry/dimension/', dim)
      quad=make_quadrature(vertices=dim+1, dim=dim, degree=1)
      shape=make_element_shape(vertices=dim+1, dim=dim, degree=1, quad=quad)
      ! allocate the mesh with n nodes and 1 element
      call allocate(mesh, n, 1, shape, "Mesh")
      call deallocate(shape)
      call deallocate(quad)
      
      ! setup trivial petsc numbering
      call allocate(petsc_numbering, n, 1)
      universal_nodes=n
      components=1
      
    end if
    
    if (components==1) then
      
      if (IsParallel()) then
        
        call redistribute_matrix(matrix, x, rhs, petsc_numbering)

      end if    
    
      call allocate(A, matrix, petsc_numbering, petsc_numbering, "PetscReadNSolveMatrix")

      ! this might not be the full name, but it's only for log output:
      call get_option(trim(option_path)//'/name', field_name)
      call allocate(x_field, mesh, field_name)
      x_field%option_path=option_path
      call allocate(rhs_field, mesh, "RHS") 
      
      call petsc2field(x, petsc_numbering, x_field, rhs_field)
      call petsc2field(rhs, petsc_numbering, rhs_field, rhs_field)
      
      call VecDestroy(rhs, ierr)
      call VecDestroy(x, ierr)

      ! prevent rewriting the matrixdump on failure
      call add_option(trim(complete_solver_option_path(x_field%option_path))//'/no_matrixdump', stat=stat)
      
      ewrite(1,*) 'Going into petsc_solve'
      ewrite(1,*) '-------------------------------------------------------------'
      
      if (read_state) then
         call petsc_solve(x_field, A, rhs_field, states(istate))
      else
         call petsc_solve(x_field, A, rhs_field)
      end if
      
      ewrite_minmax(x_field)
      
      ewrite(1,*) '-------------------------------------------------------------'
      ewrite(1,*) 'Finished petsc_solve'
      
      call deallocate(A)
      call deallocate(x_field)
      call deallocate(rhs_field)
    
    end if

    if (associated(states)) then
      do istate=1, size(states)
        call deallocate(states(istate))
      end do
    else
      call deallocate(mesh)
    end if
    call deallocate(petsc_numbering)
    
  end subroutine petsc_readnsolve_flml
    
  subroutine petsc_readnsolve_vector(mesh, n, universal_nodes, option_path, x, matrix, rhs, &
      state, read_state)
    
    type(mesh_type), intent(inout):: mesh
    integer, intent(in):: n ! matrixdump size
    integer, intent(in):: universal_nodes ! mesh size
    character(len=*), intent(in):: option_path
    Mat, intent(inout):: matrix
    Vec, intent(inout):: x, rhs
    type(state_type), intent(in):: state
    logical, intent(in):: read_state ! have we actually called populate state fully?
    
    type(petsc_csr_matrix):: A
    type(halo_type), pointer:: halo
    type(vector_field):: x_field, rhs_field
    type(petsc_numbering_type):: petsc_numbering
    character(len=FIELD_NAME_LEN):: field_name
    PetscErrorCode:: ierr
    integer:: components, stat
    
    components=n/universal_nodes
    ewrite(1,*) "Number of nodes in the mesh is an integer multiple of the matrixdump size"
    ewrite(1,*) "Assuming it's a vector field with"
    ewrite(1,*) components, " components and ", universal_nodes, " nodes."
    
    ! redo the petsc_numbering, this time with the right n/o nodes and components
    if (isparallel()) then
       halo => mesh%halos(halo_count(mesh)) 
       call allocate(petsc_numbering, node_count(mesh), components, &
            & halo)
       call redistribute_matrix(matrix, x, rhs, petsc_numbering)
    else
       call allocate(petsc_numbering, node_count(mesh), components)
    end if
    
    call allocate(A, matrix, petsc_numbering, petsc_numbering, "PetscReadNSolveMatrix")

    ! this might not be the full name, but it's only for log output:
    call get_option(trim(option_path)//'/name', field_name)
    call allocate(x_field, components, mesh, name=field_name)
    x_field%option_path=option_path
    call allocate(rhs_field, components, mesh, "RHS") 
    
    call petsc2field(x, petsc_numbering, x_field)
    call petsc2field(rhs, petsc_numbering, rhs_field)
    
    call VecDestroy(rhs, ierr)
    call VecDestroy(x, ierr)

    ! prevent rewriting the matrixdump on failure
    call add_option(trim(complete_solver_option_path(x_field%option_path))//'/no_matrixdump', stat=stat)
    
    ewrite(1,*) 'Going into petsc_solve'
    ewrite(1,*) '-------------------------------------------------------------'
    
    if (read_state) then
      call petsc_solve(x_field, A, rhs_field, state=state)
    else
      call petsc_solve(x_field, A, rhs_field)
    end if
    
    ewrite_minmax(x_field)
    
    ewrite(1,*) '-------------------------------------------------------------'
    ewrite(1,*) 'Finished petsc_solve'
    
    call deallocate(A)
    
    call deallocate(x_field)
    call deallocate(rhs_field)
    call deallocate(petsc_numbering)
    
  end subroutine petsc_readnsolve_vector
    
  function workout_option_path(field, fail)
  character(len=OPTION_PATH_LEN):: workout_option_path
  character(len=*), intent(in):: field
  logical, optional, intent(out):: fail
  
    character(len=FIELD_NAME_LEN):: phase_name, field_name
    integer i, stat
    
    if (field(1:1)=='/') then
      ! complete option path is provided
      call get_option(trim(field)//'/name', field_name, stat=stat)
      if (stat/=0) then
        ewrite(-1,*) "Option path:", field
        FLExit("Provided option path is not valid")
      end if
      workout_option_path=field
    else
      ! assume it's the field name
      
      ! search for double colon
      do i=1, len_trim(field)-1
        if (field(i:i+1)=='::') exit
      end do
      if (i<len_trim(field)-1) then
        phase_name=field(1:i-1)
        field_name=field(i+2:)
      else
        if (option_count('/material_phase')>1) then
          ewrite(-1,*) "For multi-material/phase you need to provide&
            &the material_phase and the field name, e.g. SolidPhase::Pressure"
          FLExit("Missing material_phase name")
        end if
        call get_option('/material_phase[0]/name', phase_name)
        field_name=field
      end if
      ewrite(2,*) "Phase name: ", trim(phase_name)
      ! try scalar/vector/tensor field directly under material_phase
      workout_option_path='/material_phase::'//trim(phase_name)// &
        '/scalar_field::'//trim(field_name)
      if (.not. have_option(workout_option_path)) then
        workout_option_path='/material_phase::'//trim(phase_name)// &
          '/vector_field::'//trim(field_name)
        if (.not. have_option(workout_option_path)) then
          if (present(fail)) then
            fail=.true.
            return
          end if
          ewrite(-1,*) "Field name: ", trim(field_name)
          ewrite(-1,*) "Cannot find specified field directly under /material_phase"
          ewrite(-1,*) "If it is not, you have to specify the full option_path"
          FLExit("Unable to workout option path")
        end if
      end if
    end if
    ewrite(2,*) "Field name: ", trim(field_name)
    ewrite(2,*) "Full option_path: ", trim(workout_option_path)
    if (present(fail)) then
      fail=.false.
    end if
  
  end function workout_option_path
  
  subroutine workout_option_path_from_state(states, field, &
     option_path, istate)
  type(state_type), dimension(:), intent(in):: states
  ! the field name or option path specified on the command line
  character(len=*), intent(in):: field
  
  character(len=*), intent(out):: option_path
  integer, intent(out):: istate

    type(vector_field), pointer:: vfield
    type(scalar_field), pointer:: sfield
    character(len=FIELD_NAME_LEN) phase_name, field_name
    integer i, stat
    
    istate=0
    
    ! try to work out name of material_phase
    if (field(1:1)=='/') then
      ! we've been given an option_path
      option_path=field
      if (field(1:17)=='/material_phase::') then
        do i=18, len(field)
          if (field(i:i)=='/') exit
        end do
        phase_name=field(18:i-1)
      else if (size(states)==1) then
        istate=1
        return
      else
        ewrite(-1,*) "Try PhaseName::FieldName instead"
        FLAbort("Can't work out material_phase name from option path")
      end if
    else
      ! search for double colon
      do i=1, len_trim(field)-1
        if (field(i:i+1)=='::') exit
      end do
      if (i<len_trim(field)-1) then
        field_name=field(i+2:)
        phase_name=field(1:i-1)
      else if (size(states)==1) then
        field_name=field
        phase_name=states(1)%name
        istate=1
      else
        ewrite(-1,*) "For multi-material/phase you need to provide&
          &the material_phase and the field name, e.g. SolidPhase::Pressure"
        FLExit("Missing material_phase name")
      end if
    end if
    
    if (istate==0) then
      ! now find the right state
      do i=1, size(states)
        if (states(i)%name==phase_name) exit
      end do
      if (i>size(states)) then
        ewrite(-1,*) "Material_phase name:", trim(phase_name)
        FLExit("This phase name is not known in the flml")
      end if
      istate=i
    end if    
    
    if (field(i:i)=='/') then
      ! we have the option_path already
      return
    end if
    
    ! now pull out the field from state to find its option_path
    sfield => extract_scalar_field( states(istate), field_name, stat=stat)
    if (stat==0) then
      option_path=sfield%option_path
    else
      ! maybe a vector field then?
      vfield => extract_vector_field( states(istate), field_name, stat=stat)
      if (stat/=0) then
        ewrite(-1,*) "Field name:", trim(field_name)
        ewrite(-1,*) "Material_phase name:", trim(phase_name)
        FLExit("Not a field in this material_phase")
      end if
      option_path=vfield%option_path
    end if      

  end subroutine workout_option_path_from_state
    
  subroutine redistribute_matrix(matrix, x, rhs, petsc_numbering)
  Mat, intent(inout):: matrix
  Vec, intent(inout):: x, rhs
  type(petsc_numbering_type), intent(in):: petsc_numbering
    
    VecScatter scatter
    IS row_indexset, col_indexset
    Mat new_matrix
    Vec new_x, new_rhs
    PetscErrorCode ierr
    integer, dimension(:), allocatable:: allcols, unns
    integer i, n, m, ncomponents
    
    integer mm,nn
    
    n=petsc_numbering%nprivatenodes ! local length
    ncomponents=size(petsc_numbering%gnn2unn, 2)
    allocate(unns(1:n*ncomponents))
    unns=reshape( petsc_numbering%gnn2unn(1:n,:), (/ n*ncomponents /))
#if PETSC_VERSION_MINOR>=2
    ! for petsc 3.2 we have an extra PetscCopyMode argument
    call ISCreateGeneral(MPI_COMM_FEMTOOLS, &
       size(unns), unns, PETSC_COPY_VALUES, row_indexset, ierr)
#else
    call ISCreateGeneral(MPI_COMM_FEMTOOLS, &
       size(unns), unns, row_indexset, ierr)
#endif
       
    m=petsc_numbering%universal_length ! global length
       
    ! we only ask for owned columns (although presumably
    ! still all columns of owned rows are stored locally)
    ! we only deal with square matrices (same d.o.f. for rows and columns)
    ! in fluidity, so we can simply reuse row_indexset as the col_indexset
    call MatGetSubMatrix(matrix, row_indexset, row_indexset, &
       MAT_INITIAL_MATRIX, new_matrix, ierr)

    ! destroy the old read-in matrix and replace by the new one
    call MatDestroy(matrix, ierr)
    matrix=new_matrix
    
    call matgetsize(matrix, mm, nn, ierr)
    ewrite(2,*) "Matrix global size", mm, nn
    call matgetlocalsize(matrix, mm, nn, ierr)
    ewrite(2,*) "Matrix local size", mm, nn
    
    ! create a Vec according to the proper partioning:
    call VecCreateMPI(MPI_COMM_FEMTOOLS, n*ncomponents, m, new_x, ierr)
    ! fill it with values from the read x by asking for its row numbers
    call VecScatterCreate(x, row_indexset, new_x, PETSC_NULL_OBJECT, &
       scatter, ierr)
    call VecScatterBegin(scatter, x, new_x, INSERT_VALUES, &
       SCATTER_FORWARD, ierr)
    call VecScatterEnd(scatter, x, new_x, INSERT_VALUES, &
       SCATTER_FORWARD, ierr)
    ! destroy the read x and replace by new_x
    call VecDestroy(x, ierr)
    x=new_x
    
    ! do the same for the rhs Vec:
    call VecDuplicate(new_x, new_rhs, ierr)
    call VecScatterBegin(scatter, rhs, new_rhs, INSERT_VALUES, &
       SCATTER_FORWARD, ierr)
    call VecScatterEnd(scatter, rhs, new_rhs, INSERT_VALUES, &
       SCATTER_FORWARD, ierr)
    call VecDestroy(rhs, ierr)
    rhs=new_rhs
    
    call VecScatterDestroy(scatter, ierr)
    call ISDestroy(row_indexset, ierr)    
    
  end subroutine redistribute_matrix

  subroutine petsc_readnsolve_options(filename, flml, field, &
    zero_init_guess, &
    scipy, random_rhs)
  character(len=*), intent(out):: filename, flml, field
  logical, intent(out):: zero_init_guess, scipy, random_rhs

    PetscBool flag
    PetscErrorCode ierr
    
    call PetscOptionsGetString('prns_', '-filename', filename, flag, ierr)
    if (.not. flag) then
      filename='matrixdump'
    end if

    call PetscOptionsGetString('prns_', '-flml', flml, flag, ierr)
    if (.not. flag) then
      flml=''
    end if
    
    call PetscOptionsGetString('prns_', '-field', field, flag, ierr)
    if (.not. flag) then
      field=''
    end if
    
    call PetscOptionsHasName('prns_', '-zero_init_guess', zero_init_guess, ierr)

    call PetscOptionsHasName('prns_', '-scipy', scipy, ierr)

    call PetscOptionsHasName('prns_', '-random_rhs', random_rhs, ierr)
    
    call PetscOptionsGetInt('prns_', '-verbosity', current_debug_level, flag, ierr)
    if (.not.flag) then
      current_debug_level=3
    end if
      
  end subroutine petsc_readnsolve_options
  
end subroutine petsc_readnsolve

