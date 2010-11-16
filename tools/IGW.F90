#include "fdebug.h"

subroutine IGW
  ! A small program to solve for inertia-gravity waves
  ! using P1dg triangular elements for velocity
  ! and P2 triangular elements for height
  ! weak boundary conditions u.n = 0
  use mesh_files
  use fields
  use FEtools
  use DGtools
  use elements
  use sparse_tools
  use vtk_interfaces
  use transform_elements
  use solvers
  use petsc_tools
  use sparsity_patterns
  use signal_vars
  use state_module
  use solvers
  use global_parameters, only : current_debug_level, PYTHON_FUNC_LEN
  use spud 
  use conservative_interpolation_module
  use dg_interpolation_module
  use balanced_interpolation_module

  implicit none
#include "finclude/petsc.h"
#if PETSC_VERSION_MINOR==0
#include "finclude/petscmat.h"
#include "finclude/petscvec.h"
#include "finclude/petscviewer.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscmg.h"
#include "finclude/petscis.h"
#endif

  type(vector_field), target :: u,positions, positions_h,error_u, initial_u, &
       X_surface, exact_u
  type(scalar_field), target :: h, CTu,streamfunction, exact_h, error_h, &
       h_surface
  integer :: degree, quad_degree
  type(quadrature_type), target :: quad,f_quad
  type(element_type), target :: X_shape, u_shape, h_shape, &
       X_shape_f,u_shape_f,h_shape_f, vtk_shape
  type(dynamic_csr_matrix) :: h_mass, u_mass
  type(mesh_type) :: h_mesh,u_mesh,vtk_mesh
  type(dynamic_csr_matrix) :: C1T, C2T, C3T, big_m
  type(dynamic_csr_matrix) :: u_inverse_mass
  type(csr_matrix) :: u_inverse_mass_static
  type(dynamic_csr_matrix) :: CMC, C2MC2T, C3MC3T, &
       MC1T,MC2T,MC3T
  type(csr_matrix) :: CMC_static
  type(csr_matrix) :: u_mass_static, C1T_static, C2T_static, &
       C3T_static, h_mass_static
  type(csr_matrix) :: big_m_static
  type(csr_sparsity) :: psi_sparsity
  type(csr_matrix) :: psi_mat
  real, dimension(:), pointer :: Xu,Yu,Xh,Yh
  real, dimension(:), allocatable :: rhs, big_vec, tmpu,tmph
  character(len=100) :: tmpbuf
  real :: t, dt, tmax=1, tdump, Ro, Fr
  real,parameter :: pi = 3.141592654
  integer :: number_of_iterations
  integer :: dump = 0
  ! Arguments for handling the command line
  character(len=256) :: filename, u_input,h_input,mesh,exact_h_input, &
       exact_u_input
  integer :: status, i, count, dumcount, ndump,ele,j
  integer, dimension(:), pointer :: row_j
  character(len=100) :: dumcount_str, fmt,buffer,buf
  character(len=PYTHON_FUNC_LEN), save :: func
  integer :: u_degree=1, h_degree=2, u_cty=-1
  type(state_type) :: state
  logical :: patrick = .false.
  logical :: steady_state=.false., &
       balanced_u=.false., get_Streamfunction=.false.
  logical :: projection_test=.false.,compute_error=.false.
  namelist/IGW_data/tmax, tdump, dt, fr, ro, u_degree,h_degree, &
       u_cty, balanced_u, steady_state,mesh,u_input,h_input,exact_h_input, &
       get_streamfunction, projection_test, compute_error, exact_u_input, &
       matrices
  logical :: file_exists, matrices=.false.
  integer :: dim, loc, nnodes, nelements, node_attributes, unit, io1
  integer, allocatable, dimension(:) :: list
  real :: norm_error
  !debugging bits
  type(vector_field), target :: positions_u
  type(vector_field) :: source_positions, source_velocity, gravity
  type(scalar_field) :: source_height, vbd
  type(mesh_type) :: source_u_mesh, source_h_mesh
  type(state_type) :: source_state, target_state, &
       balance_state_old(1), balance_state_new(1)
  integer :: node
  integer :: stat

#ifdef HAVE_PETSC
  PetscTruth :: flag
  PetscErrorCode :: ierr
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

  current_debug_level = 2

  ewrite(1,*) 'subroutine IGW'

  call PetscOptionsGetString('igw_', '-filename', filename, flag, ierr)
  if (.not. flag) then
     call usage
     stop
  end if

  filename=trim(filename)

  ewrite(2,*) 'Reading parameters'

  unit = free_unit()
  open(unit=unit, file=trim(trim(filename)//".dat"), status='old', &
       iostat=io1)
          
  if(io1.ne.0) then
     ewrite(-1,*) 'Looked for ', trim(trim(filename)//".dat")
     FLExit('Could not read from .dat file')
  end if
  
  read(unit, IGW_data)
  close(unit) 

  ewrite(2,*) 'tmax, tdump, dt, fr, ro, steady_state', tmax, tdump, dt, fr, &
       ro, steady_state
  ewrite(2,*) 'balanced_u',balanced_u
  ewrite(2,*) 'u_degree, h_degree, u_cty', u_degree,h_degree, u_cty

  ewrite(2,*) 'Getting mesh file information'

  call identify_mesh_file(trim(mesh), dim, loc, nnodes, nelements, &
       node_attributes)

  ewrite(2,*) 'dim = ', dim

  quad_degree = 6

  ewrite(2,*) 'Getting quadrature'

  quad=make_quadrature(loc=loc, dimension=dim, degree=quad_degree)
  f_quad=make_quadrature(loc=loc-1, dimension=dim-1, degree=quad_degree-1)
  
  ewrite(2,*) 'Getting shape functions'

  ! Shape functions for positions (linear)
  X_shape=make_element_shape(loc=loc, dimension=dim, &
       degree=1, quad=quad)
  X_shape_f=make_element_shape(loc=loc-1, dimension=dim-1, &
       degree=1, quad=f_quad)

  ewrite(2,*) 'reading mesh'
  ewrite(2,*) 'loc = ',loc,'dim = ',dim

  positions=read_mesh_files(trim(mesh), X_shape)

  ewrite(2,*) 'adding faces'

  !call add_faces(positions%mesh)

  ewrite(2,*) 'getting shapes'

  ! Shape functions for velocity and height
  u_shape=make_element_shape(loc=loc, &
       dimension=dim, degree=u_degree, quad=quad)
  u_shape_f=make_element_shape(loc=loc-1, &
       dimension=dim-1, degree=u_degree, quad=f_quad)
  h_shape=make_element_shape(loc=loc, &
       dimension=dim, degree=h_degree, quad=quad)
  if(h_degree>1.or.u_degree>1) then
     vtk_shape=make_element_shape(loc=loc, &
          dimension=dim, degree=2, quad=quad)
  else
     vtk_shape=make_element_shape(loc=loc, &
          dimension=dim, degree=1, quad=quad)
  end if
  h_shape_f=make_element_shape(loc=loc-1, &
       dimension=dim-1, degree=h_degree, quad=f_quad)

  !connectivity for velocity and height
  ewrite(2,*) 'Getting u connectivity'
  if(u_cty==0 .and. u_degree==1) then
     u_mesh = positions%mesh
     u_mesh%name = "u_mesh"
     call incref(u_mesh)
  else
     u_mesh = make_mesh(positions%mesh,u_shape,u_cty,'u_mesh')
     !call add_faces(u_mesh, model=positions%mesh)
  end if
  ewrite(2,*) 'Getting h connectivity'
  if(h_degree==1) then
     h_mesh = positions%mesh
     h_mesh%name = "h_mesh"
     call incref(h_mesh)
  else
     h_mesh = make_mesh(positions%mesh,h_shape,0,'h_mesh')
     !call add_faces(h_mesh, model=positions%mesh)
  end if
  ewrite(2,*) 'Getting vtk connectivity'
  if(vtk_shape%degree==1) then
     vtk_mesh = positions%mesh
     vtk_mesh%name = "vtk_mesh"
     call incref(vtk_mesh)
  else
     vtk_mesh = make_mesh(positions%mesh,vtk_shape,0,'vtk_mesh')
     !call add_faces(vtk_mesh, model=positions%mesh)
  end if
  !fields for velocity and height
  ewrite(2,*) 'Allocating fields'

  allocate( tmpu(node_count(u_mesh)), tmph(node_count(h_mesh)) )
  allocate( big_vec(2*node_count(u_mesh)+node_count(h_mesh)), &
       rhs(2*node_count(u_mesh)+node_count(h_mesh)) )
  h = wrap_scalar_field(h_mesh, &
       big_vec(2*node_count(u_mesh)+1: &
       2*node_count(u_mesh)+node_count(h_mesh)), &
       'height')
  u = wrap_vector_field(u_mesh, &
       big_vec(1:node_count(u_mesh)), &
       big_vec(node_count(u_mesh)+1:2*node_count(u_mesh)), &
       name='velocity')
  call allocate(positions_h,dim,h_mesh,name='positions_h')  
  call allocate(positions_u,dim,u_mesh,name='positions_u')
  call remap_vector_field(positions,positions_h)
  call remap_vector_field(positions,positions_u)
  if(steady_state) then
     call allocate_vector_field(error_u,2,u%mesh,'u_error')
     call allocate_vector_field(initial_u,2,u%mesh,'u_initial')
  end if

  call allocate_vector_field(exact_u,2,u%mesh,'exact_u')
  call allocate_scalar_field(exact_h,h%mesh,'exact_h')

  call allocate_scalar_field(error_h,h%mesh,'error_h')
  exact_h%val = 0.
  error_h%val = 0.
  call allocate_scalar_field(ctu,h%mesh,'div u')

  ewrite(2,*) 'number of nodes for u,h', node_count(u), node_count(h)

  !get mass matrix for u
  ewrite(2,*) 'getting mass matrix for u'
  call allocate(u_mass,node_count(u),node_count(u))
  call assemble_mass_vec(positions,u,u_mass)

  ewrite(2,*) 'allocating mass matrix for h'
  call allocate(h_mass,node_count(h),node_count(h))
  ewrite(2,*) 'assembling mass for h'
  call assemble_mass(positions,h,h_mass)

  ewrite(2,*) 'allocating CT'

  call allocate(C1T,node_count(h),node_count(u))
  call allocate(C2T,node_count(h),node_count(u))
  call allocate(C3T,node_count(h),node_count(u))
  call allocate(big_m,node_count(u)*2+node_count(h), &
       node_count(u)*2+node_count(h))

  ewrite(2,*) 'assembling CT'

  call assemble_CT(positions,u,h,C1T,C2T,C3T)

  if(matrices) then
     u_mass_static = dcsr2csr(u_mass)
     h_mass_static = dcsr2csr(h_mass)
     C1T_static = dcsr2csr(C1T)
     C2T_static = dcsr2csr(C2T)
     call dump_matrices(u_mass_static,h_mass_static,&
          &C1T_static,C2T_static,h,u,positions)
     ewrite(0,*) 'stopping'
     stop
  end if

  ewrite(2,*) 'forming big matrix'
  
  call make_big_m(big_m,u,h,u_mass,h_mass,C1T,C2T,Ro,Fr,dt)

  ewrite(2,*) 'staticising big matrix'

  big_m_static = dcsr2csr(big_m)

  ewrite(2,*) 'staticising other matrices'
  u_mass_static = dcsr2csr(u_mass)
  h_mass_static = dcsr2csr(h_mass)
  C1T_static = dcsr2csr(C1T)
  C2T_static = dcsr2csr(C2T)
  
  ewrite(2,*) 'Initialising u and h'

  ewrite(2,*) 'setting h from python'
  inquire(file=trim(h_input),exist=file_exists)  
  if (.not.file_exists) FLExit('Couldnt find ' // trim(h_input) // ' file')
  unit=free_unit() 
  open(unit, file=trim(h_input), action="read",&
       & status="old")
  read(unit, '(a)', end=43) func
  ! Read all the lines of the file and put in newlines between them.
  do
     read(unit, '(a)', end=43) buffer
     func=trim(func)//achar(10)//trim(buffer)
  end do
43 func=trim(func)//achar(10)
  close(unit)

  ewrite(2,*) func

  if (patrick) then
    write(0,*) "patrick interpolating h"
    h%option_path = "/fields/height/prognostic/galerkin_projection/continuous"
    call set_solver_options(h, ksptype='cg', pctype='eisenstat', &
         rtol=1.0e-10, max_its=10000)
    h%option_path = "/fields/height"
    source_positions = read_mesh_files("source", X_shape)
    source_h_mesh = make_mesh(source_positions%mesh, h_shape, 0, 'h_mesh')
    source_u_mesh = make_mesh(source_positions%mesh, u_shape, u_cty, 'u_mesh')
    call allocate(source_height, source_h_mesh, 'height')
    call allocate(source_velocity, source_positions%dim, &
         source_u_mesh, 'velocity')

    call set_from_python_function(source_height, trim(func), &
         source_positions, 0.0)

    call insert(source_state, source_height, 'Pressure')
    call insert(source_state, source_positions, 'Coordinate')

    call insert(target_state, h, 'Pressure')
    call insert(target_state, positions, 'Coordinate')

    call interpolation_galerkin(source_state, target_state)
    call deallocate(source_state)
    call deallocate(target_state)

    do node=1,node_count(h)
      if (node_val(positions_h, node, 1) == 0.0) then
        call set(h, node, 0.0)
      end if
      if (node_val(positions_h, node, 2) == 0.0) then
        call set(h, node, 0.0)
      end if
    end do
  else
    call set_from_python_function(h,trim(func), positions, 0.0)
  end if

  if(compute_error) then
     
     ewrite(2,*) 'setting exact h from python'
     inquire(file=trim(exact_h_input),exist=file_exists)  
     if (.not.file_exists) FLExit('Couldnt find _h.py file')
     unit=free_unit() 
     open(unit, file=trim(exact_h_input), action="read",&
          & status="old")
     read(unit, '(a)', end=44) func
     ! Read all the lines of the file and put in newlines between them.
     do
        read(unit, '(a)', end=44) buffer
        func=trim(func)//achar(10)//trim(buffer)
     end do
44   func=trim(func)//achar(10)
     close(unit)
     
     ewrite(2,*) func
     call set_from_python_function(exact_h,trim(func), positions, 0.0)
  end if

  if(balanced_u) then
     ewrite(2,*) 'getting balanced u'
     call mult_T(tmpu,C2T_static,h%val)
     tmpbuf = "/temporary/fix"
     call set_solver_options(field_option_path=tmpbuf, &
          ksptype="preonly|", PCTYPE="lu", &
          max_its=1000, atol=1.0e-30, rtol=1.0e-14, &
          start_from_zero=.true.)
     tmpu = Ro*tmpu
     call petsc_solve(u%val(1,:),u_mass_static,tmpu, &
          option_path="/temporary/fix")
     ewrite(2,*) 'getting balanced v'
     call mult_T(tmpu,C1T_static,h%val)
     tmpu = -Ro*tmpu
     call petsc_solve(u%val(2,:),u_mass_static,tmpu, &
          option_path="/temporary/fix")
  else
     inquire(file=trim(u_input),exist=file_exists)  
     if (.not.file_exists) FLExit('Couldnt find _u.py file')
     unit=free_unit() 
     open(unit, file=trim(u_input), action="read",&
          & status="old")
     read(unit, '(a)', end=42) func
     ! Read all the lines of the file and put in newlines between them.
     do
        read(unit, '(a)', end=42) buffer
        func=trim(func)//achar(10)//trim(buffer)
     end do
42   func=trim(func)//achar(10)
     close(unit)
     call set_from_python_function(u,trim(func), positions,0.0)
   end if

   if (patrick) then

     write(0,*) "patrick interpolating u"

     ! First we interpolate from IGW's fields to source_*,
     ! then we interpolate back from source_* to IGW's fields

     u%option_path = "/fields/u"
     source_velocity%option_path = "/fields/u"
     call set_option("/fields/u/prognostic/galerkin_projection/balanced_interpolation", 1.0, stat=stat)
     tmpbuf = "/fields/u/prognostic/galerkin_projection/balanced_interpolation"
     call set_solver_options(tmpbuf, ksptype='cg', &
                             pctype='eisenstat', rtol=1.0e-10, max_its=10000)

     call set_option("/physical_parameters/coriolis/f_plane/f", ro/fr**2)

     call insert(source_state, u, 'Velocity')
     call insert(source_state, positions, 'Coordinate')
     call allocate(vbd, u%mesh, "VelocityBuoyancyDensity")
     call zero(vbd)
     call insert(source_state, vbd, "VelocityBuoyancyDensity")
     call deallocate(vbd)
     call allocate(gravity, u%dim, u%mesh, "GravityDirection")
     call zero(gravity)
     call insert(source_state, gravity, "GravityDirection")
     call deallocate(gravity)
     call insert(source_state, h, "Pressure")

     call insert(balance_state_old(1), u, 'Velocity')
     call insert(balance_state_old(1), positions, 'Coordinate')

     call insert(target_state, source_velocity, 'Velocity')
     call insert(target_state, source_velocity, 'Coordinate')

     call insert(balance_state_new(1), source_velocity, 'Velocity')
     call insert(balance_state_new(1), source_positions, 'Coordinate')
     call insert(balance_state_new(1), source_velocity%mesh, "VelocityMesh")

     call set_option("/mesh_adaptivity/hr_adaptivity/debug/write_metric_stages", 1.0, stat=stat)
     write(0,*) "balance"
     call balanced_interpolation_prepare_old((/source_state/), balance_state_old)
     call balanced_interpolation_prepare_new(balance_state_new)
     call print_state(balance_state_new(1), 0)

     call interpolation_galerkin(balance_state_old(1), balance_state_new(1))

     call print_state(balance_state_new(1), 0)
     call insert(balance_state_new(1), source_height, "Pressure")
     call balanced_interpolation_conclude(balance_state_new)

     ! -------------------- and now for the second go

     call insert(source_state, source_velocity, 'Velocity')
     call insert(source_state, source_positions, 'Coordinate')
     call allocate(vbd, source_velocity%mesh, "VelocityBuoyancyDensity")
     call zero(vbd)
     call insert(source_state, vbd, "VelocityBuoyancyDensity")
     call deallocate(vbd)
     call allocate(gravity, u%dim, u%mesh, "GravityDirection")
     call zero(gravity)
     call insert(source_state, gravity, "GravityDirection")
     call deallocate(gravity)
     call insert(source_state, source_height, "Pressure")

     call insert(balance_state_old(1), source_velocity, 'Velocity')
     call insert(balance_state_old(1), source_positions, 'Coordinate')

     call insert(target_state, u, 'Velocity')
     call insert(target_state, positions, 'Coordinate')

     call insert(balance_state_new(1), u, 'Velocity')
     call insert(balance_state_new(1), positions, 'Coordinate')
     call insert(balance_state_new(1), u%mesh, "VelocityMesh")

     call balanced_interpolation_prepare_old((/source_state/), balance_state_old)
     call balanced_interpolation_prepare_new(balance_state_new)

     call interpolation_galerkin(balance_state_old(1), balance_state_new(1))

     call insert(balance_state_new(1), h, "Pressure")
     call balanced_interpolation_conclude(balance_state_new)

     call deallocate(source_state)
     call deallocate(target_state)
     call deallocate(balance_state_old(1))
     call deallocate(balance_state_new(1))

     call deallocate(source_u_mesh)
     call deallocate(source_h_mesh)
     call deallocate(source_velocity)
     call deallocate(source_height)
     call deallocate(source_positions)
  end if

  if(get_streamfunction) then
     ewrite(2,*) 'getting stream function'
     call allocate_scalar_field(streamfunction,h%mesh,'StreamFunction')
     call zero(streamfunction)
     call calculate_stream_function_2d(streamfunction,positions,u,status)

     if(status.ne.0) then
        FLAbort('failed to make stream function')
     end if

     call vtk_write_fields('streamfunction', &
          index=0, position=positions, &
          model=vtk_mesh, sfields=(/h,streamfunction/),vfields=(/u/))

     ewrite(0,*) 'Finished after dumping streamfunction'
     stop
  end if

  if(projection_test) then
     call allocate(u_inverse_mass,node_count(u),node_count(u))
     call get_dg_inverse_mass_matrix(u_inverse_mass,u_mesh,positions)
     u_inverse_mass_static = dcsr2csr(u_inverse_mass)

     ewrite(2,*) 'calling matmul_T mc1t'
     MC1T = matmul_T(C1T,u_inverse_mass,check=.true.)
     
     ewrite(2,*) 'calling matmul_T mc2t'
     MC2T = matmul_T(C2T,u_inverse_mass,check=.true.)
     !call matcheck(MC1t,MC2T)
     if(dim==3) then
        ewrite(2,*) 'calling matmul_T mc3t'
        MC3T = matmul_T(C3T,u_inverse_mass,check=.true.)
     end if
     ewrite(2,*) 'calling matmul_T cmc'
     CMC = matmul_T(C1T,MC1T,check=.true.)
     ewrite(2,*) 'calling matmul_T c2mc2t'
     C2MC2T = matmul_T(C2T,MC2T,check=.true.)
     if(dim==3) then
        ewrite(2,*) 'calling matmul_T c3mc3t'
        C3MC3T = matmul_T(C3T,MC3T,check=.true.)
     end if
     ewrite(2,*) 'calling addto'
     call addto(CMC,C2MC2T)
     
     ewrite(2,*) 'staticising cmc'
     cmc_static = dcsr2csr(cmc)
     
     ewrite(2,*) 'fixing zero value'
     call addto(CMC_static,1,1,INFINITY)

     ewrite(2,*) 'computing div u'
     
     ctu%val = 0.0
     call mult(tmph,C1T_static,u%val(1,:))
     ctu%val = ctu%val + tmph
     call mult(tmph,C2T_static,u%val(2,:))
     ctu%val = ctu%val + tmph

     ewrite(2,*) 'max(abs(ctu))=',maxval(ctu%val), minval(ctu%val)

     ewrite(2,*) 'Projecting u to div u = 0'
     ewrite(2,*) 'solving'

     call petsc_solve(tmph,cmc_static,ctu%val, &
          startfromzero=.true., checkconvergence=.true., &
          iterations=number_of_iterations)

     ewrite(2,*) 'max(abs(p))=',maxval(tmph), minval(tmph)

     !want C^Tu=0, write u = u^* + M^-1Cp,
     !C^Tu = C^Tu^* + C^TM^-1Cp

     ewrite(2,*) 'adding correction'
     call mult_T(tmpu,C1T_static,tmph)
     call petsc_solve(tmpu,u_mass_static,tmpu, &
          startfromzero=.true., checkconvergence=.true., &
          iterations=number_of_iterations)
     !call mult(tmpu,u_inverse_mass_static,tmpu)
     u%val(1,:) = u%val(1,:) - tmpu
     ewrite(2,*) 'max u1 correction', maxval(abs(tmpu))
     call mult_T(tmpu,C2T_static,tmph)
     call petsc_solve(tmpu,u_mass_static,tmpu, &
          startfromzero=.true., checkconvergence=.true., &
          iterations=number_of_iterations)
     !call mult(tmpu,u_inverse_mass_static,tmpu)
     ewrite(2,*) 'max u2 correction', maxval(abs(tmpu))
     u%val(2,:) = u%val(2,:) - tmpu

     ewrite(2,*) 'computing div u'
     
     ctu%val = 0.0
     call mult(tmph,C1T_static,u%val(1,:))
     ctu%val = ctu%val + tmph
     call mult(tmph,C2T_static,u%val(2,:))
     ctu%val = ctu%val + tmph

     ewrite(2,*) 'max(abs(ctu))=',maxval(ctu%val), minval(ctu%val)

     ewrite(2,*) 'projecting to skew gradient of p'
     !if u = M^-1(-C_2,C_1)p
     !then -C_2u_1 + C_1u_2 = CMCp

     ewrite(2,*) 'computing curl of u'

     ctu%val = 0.0
     call mult(tmph,C2T_static,u%val(1,:))
     ctu%val = ctu%val - tmph
     call mult(tmph,C1T_static,u%val(2,:))
     ctu%val = ctu%val + tmph

     ewrite(2,*) 'max(abs(ctu))=',maxval(ctu%val), minval(ctu%val)
     
     ewrite(2,*) 'solving for p'
     call petsc_solve(tmph,cmc_static,ctu%val, &
          startfromzero=.true., checkconvergence=.true., &
          iterations=number_of_iterations)

     ewrite(2,*) 'max(abs(p))=',maxval(tmph), minval(tmph)

     ewrite(2,*) 'evaluating projection'
     call allocate_vector_field(error_u,2,u%mesh,'projected u')
     call mult_T(tmpu,C1T_static,tmph)
     call petsc_solve(tmpu,u_mass_static,tmpu, &
          startfromzero=.true., checkconvergence=.true., &
          iterations=number_of_iterations)
     ewrite(2,*) 'difference for u2', maxval(abs(tmpu-u%val(2,:)))
     error_u%val(2,:) = tmpu
     call mult_T(tmpu,C2T_static,tmph)
     call petsc_solve(tmpu,u_mass_static,tmpu, &
          startfromzero=.true., checkconvergence=.true., &
          iterations=number_of_iterations)
     ewrite(2,*) 'difference for u1', maxval(abs(tmpu+u%val(1,:)))
     error_u%val(1,:) = -tmpu

     call vtk_write_fields(trim('projection'), &
          index=0, position=positions, &
          model=vtk_mesh, sfields=(/h,ctu/), &
          vfields=(/u,error_u/))

     ewrite(2,*) 'Stopping after projection test'
     stop
  end if

  ctu%val = 0.0
  call mult(tmph,C1T_static,u%val(1,:))
  ctu%val = ctu%val + tmph
  call mult(tmph,C2T_static,u%val(2,:))
  ctu%val = ctu%val + tmph
  
  if(steady_state) then
     initial_u%val(1,:) = u%val(1,:)
     initial_u%val(2,:) = u%val(2,:)
     error_u%val(1,:) = 0.
     error_u%val(2,:) = 0.
     call vtk_write_fields(trim(filename), &
          index=0, position=positions, &
          model=vtk_mesh, sfields=(/h,ctu/), &
          vfields=(/u,error_u/))
  else 
     call vtk_write_fields(trim(filename), &
          index=0, position=positions, &
          model=vtk_mesh, sfields=(/h,ctu/), &
          vfields=(/u/))
  end if

  ewrite(2,*) 'entering timestepping loop'

  t = 0.

  do
     t = t + dt
     if(t>tmax) exit
     ewrite(2,*) 'calling petsc'

     rhs = 0.

     !u mass parts
     call mult(tmpu,u_mass_static,u%val(1,:))
     !d/dt
     rhs(1:node_count(u))=rhs(1:node_count(u)) + tmpu
     !coriolis
     rhs(1+node_count(u):2*node_count(u))= &
          rhs(1+node_count(u):2*node_count(u))&
          + 0.5*dt*tmpu*Fr*Fr/Ro

     call mult(tmpu,u_mass_static,u%val(2,:))
     !d/dt
     rhs(1+node_count(u):2*node_count(u))= &
          rhs(1+node_count(u):2*node_count(u))&
          + tmpu
     !coriolis
     rhs(1:node_count(u))=rhs(1:node_count(u)) &
          -0.5*dt*tmpu*Fr*Fr/Ro
     !C1,C2 parts
     if(.true.) then
        call mult_T(tmpu,C1T_static,h%val)
        rhs(1:node_count(u)) = rhs(1:node_count(u)) - &
             0.5*dt*tmpu
        call mult_T(tmpu,C2T_static,h%val)
        rhs(1+node_count(u):2*node_count(u)) = &
             rhs(1+node_count(u):2*node_count(u)) - &
             0.5*dt*tmpu
        
        call mult(tmph,C1T_static,u%val(1,:))
        rhs(2*node_count(u)+1:2*node_count(u)+node_count(h)) = &
             rhs(2*node_count(u)+1:2*node_count(u)+node_count(h)) &
             +0.5*dt*tmph
        call mult(tmph,C2T_static,u%val(2,:))
        rhs(2*node_count(u)+1:2*node_count(u)+node_count(h)) = &
             rhs(2*node_count(u)+1:2*node_count(u)+node_count(h)) &
             +0.5*dt*tmph
     end if

     ! h mass part
     call mult(tmph,h_mass_static,h%val)
     rhs(2*node_count(u)+1:2*node_count(u)+node_count(h)) = &
          rhs(2*node_count(u)+1:2*node_count(u)+node_count(h)) &
          + tmph

     ewrite(2,*) 'max(abs(rhs))=',maxval(abs(rhs))

     tmpbuf = "/temporary/fix"
     call set_solver_options(tmpbuf, start_from_zero=.false., rtol=1.0e-11, max_its=10000)
     call petsc_solve(big_vec, big_m_static, rhs, &
          option_path="/temporary/fix")

     if(t.ge.(dump+1)*tdump) then
        addingdumps: do
           dump = dump + 1
           if(t<(dump+1)*tdump) exit
        end do addingdumps

        ewrite(1,*) 'dumping t = ',t

        ctu%val = 0.0
        call mult(tmph,C1T_static,u%val(1,:))
        ctu%val = ctu%val + tmph
        call mult(tmph,C2T_static,u%val(2,:))
        ctu%val = ctu%val + tmph

        !call petsc_solve(ctu%val,h_mass_static, ctu%val, &
        !     startfromzero=.false., checkconvergence=.true., &
        !     iterations=number_of_iterations)
        
        if (steady_state) then
           error_u%val(1,:) = u%val(1,:) - initial_u%val(1,:)
           error_u%val(2,:) = u%val(2,:) - initial_u%val(2,:)
           call vtk_write_fields(trim(filename), &
                index=dump, position=positions, &
                model=vtk_mesh, sfields=(/h,ctu/), &
                vfields=(/u,error_u/))
        else
           call vtk_write_fields(trim(filename), &
                index=dump, position=positions, &
                model=vtk_mesh, sfields=(/h,ctu/), &
                vfields=(/u/))
        end if
     end if
  end do
  
  ewrite(2,*) 'computing statistics'
  ewrite(2,*) 'max h',maxval(h%val)
  call get_wave_position(h,positions)

  if(.false.) then
     
     !setting exact values ------------------------
     
     ewrite(2,*) 'setting exact_h from python'
     inquire(file=trim(exact_h_input),exist=file_exists)  
     if (.not.file_exists) FLExit('Couldnt find ' // trim(exact_h_input) // ' file')
     unit=free_unit() 
     open(unit, file=trim(exact_h_input), action="read",&
          & status="old")
     read(unit, '(a)', end=99) func
     ! Read all the lines of the file and put in newlines between them.
     do
        read(unit, '(a)', end=99) buffer
        func=trim(func)//achar(10)//trim(buffer)
     end do
99   func=trim(func)//achar(10)
     close(unit)
     
     ewrite(2,*) func

     call set_from_python_function(exact_h,trim(func), positions, t)
     
     ewrite(2,*) 'setting exact_u from python'
     inquire(file=trim(exact_u_input),exist=file_exists)  
     if (.not.file_exists) FLExit('Couldnt find _h.py file')
     unit=free_unit() 
     open(unit, file=trim(exact_u_input), action="read",&
          & status="old")
     read(unit, '(a)', end=999) func
     ! Read all the lines of the file and put in newlines between them.
     do
        read(unit, '(a)', end=999) buffer
        func=trim(func)//achar(10)//trim(buffer)
     end do
999  func=trim(func)//achar(10)
     close(unit)
     
     ewrite(2,*) func
     
     call set_from_python_function(exact_u,trim(func), positions, t)
  else
     exact_h%val = 0.
  end if

  !---------------------------------------------
     
  call get_L2_estimates(u,h,exact_u,exact_h,positions)

  error_h%val = h%val - exact_h%val

  call vtk_write_fields(trim('exact'), &
       index=0, position=positions, &
       model=vtk_mesh, sfields=(/h,ctu,exact_h,error_h/), &
       vfields=(/u,exact_u/))
  
  ewrite(1,*) 'END program IGW '

contains

  subroutine make_big_m(big_m,u,h,u_mass,h_mass,C1T,C2T,Ro,Fr,dt)
    type(dynamic_csr_matrix), intent(inout) :: big_m
    type(dynamic_csr_matrix), intent(in) :: u_mass,h_mass,C1T,C2T
    type(scalar_field), intent(in) :: h
    type(vector_field), intent(in) :: u
    real, intent(in) :: Ro,dt,Fr
    !locals
    integer :: i,j,k
    integer, dimension(:), pointer :: row_m
    real, dimension(:), pointer :: row_val

    ewrite(1,*) 'subroutine make_big_m'

    ! U mass matrix + Coriolis
    ewrite(2,*) 'u mass parts'

    do i = 1, node_count(u)
       row_m => row_m_ptr(u_mass,i)
       row_val => row_val_ptr(u_mass,i)

       !d/dt part
       call set(big_m,i,row_m,row_val*Fr*Fr)
       call set(big_m,node_count(u)+i,node_count(u)+row_m,row_val*Fr*Fr)

       !Coriolis part
       call set(big_m,node_count(u)+i,row_m,-0.5*dt*row_val/Ro*Fr*Fr)
       call set(big_m,i,node_count(u)+row_m, 0.5*dt*row_val/Ro*Fr*Fr)
    end do
    

    if(.true.) then
       ewrite(2,*) 'C1 parts'
       
       !C1 and C1T
       do i = 1, node_count(h)
          row_m => row_m_ptr(C1T,i)       
          row_val => row_val_ptr(C1T,i)
          
          call set(big_m,row_m,2*node_count(u)+i,0.5*dt*row_val)
          call set(big_m,2*node_count(u)+i,row_m,-0.5*dt*row_val)
          
       end do
    
       ewrite(2,*) 'C2 parts'
       
       !C2 and C2T
       do i = 1, node_count(h)
          row_m => row_m_ptr(C2T,i)
          row_val => row_val_ptr(C2T,i)
          
          call set(big_m,node_count(u)+row_m,2*node_count(u)+i,0.5*dt*row_val)
          call set(big_m,2*node_count(u)+i,node_count(u)+row_m,-0.5*dt*row_val)
          
       end do
    end if

    ewrite(2,*) 'h mass parts'

    !h mass part
    do i = 1, node_count(h)
       row_m => row_m_ptr(h_mass,i)       
       row_val => row_val_ptr(h_mass,i)

       call set(big_m,2*node_count(u)+i,2*node_count(u)+row_m, &
            row_val)

    end do

    ewrite(1,*) 'subroutine make_big_m'

  end subroutine make_big_m

  subroutine assemble_CT(positions,u,h,C1T,C2T,C3T)
    type(vector_field), intent(in) :: positions, u
    type(scalar_field), intent(in) :: h
    type(dynamic_csr_matrix), intent(inout) :: C1T,C2T,C3T

    !locals
    integer :: ele

    do ele = 1, element_count(u)
       call assemble_CT_elemental(ele,positions,u,h,C1T,C2T,C3T)
    end do

  end subroutine assemble_CT
  
  subroutine assemble_CT_elemental(ele,positions,u,h,C1T,C2T,C3T)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: u
    type(scalar_field), intent(in) :: h
    type(dynamic_csr_matrix), intent(inout) :: C1T,C2T,C3T

    ! Locations of nodes.
    real, dimension(positions%dim,ele_loc(positions,ele)) :: X_ele
    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,ele)) :: detwei    
    ! Derivatives of shape function:
    real, dimension(ele_loc(h,ele), &
         ele_ngi(h,ele), positions%dim) :: dshape_h
    real, dimension(ele_loc(u,ele), &
         ele_ngi(u,ele), positions%dim) :: dshape_u
    ! Node numbers of field element.
    integer, dimension(:), pointer :: ele_u, ele_h
    ! Shape functions.
    type(element_type), pointer :: shape_u, shape_h, shape_X
    ! gradient matrix
    real, dimension(positions%dim,ele_loc(h,ele),ele_loc(u,ele)) :: grad_mat
    integer, dimension(:), pointer :: neigh
    integer :: ele_2, ni, face

    ele_u=>ele_nodes(u, ele)
    shape_u=>ele_shape(u, ele)
    ele_h=>ele_nodes(h, ele)
    shape_h=>ele_shape(h, ele)
    shape_X=>ele_shape(positions, ele)

    ! Locations of local vertices.
    X_ele=ele_val(positions, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, shape_h, &
         dshape=dshape_h, detwei=detwei)

    grad_mat = -dshape_shape(dshape_h,shape_u,detwei)

    call addto(C1T,ele_h,ele_u,grad_mat(1,:,:))
    call addto(C2T,ele_h,ele_u,grad_mat(2,:,:))
    if(positions%dim==3) then
       call addto(C3T,ele_h,ele_u,grad_mat(3,:,:))
    end if

  end subroutine assemble_CT_elemental

  subroutine assemble_mass(positions,h,mass)
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: h
    type(dynamic_csr_matrix), intent(inout) :: mass

    !locals
    integer :: ele

    do ele = 1, element_count(h)
       call assemble_mass_elemental(ele,positions,h,mass)
    end do

  end subroutine assemble_mass

  subroutine assemble_mass_vec(positions,h,mass)
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: h
    type(dynamic_csr_matrix), intent(inout) :: mass

    !locals
    integer :: ele

    do ele = 1, element_count(h)
       call assemble_mass_elemental_vec(ele,positions,h,mass)
    end do

  end subroutine assemble_mass_vec

  subroutine assemble_mass_elemental(ele,positions,h,mass)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: h
    type(dynamic_csr_matrix), intent(inout) :: mass

    ! Locations of nodes.
    real, dimension(positions%dim,ele_loc(positions,ele)) :: X_ele
    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,ele)) :: detwei    
    ! Node numbers of field element.
    integer, dimension(:), pointer :: ele_h
    ! Shape functions.
    type(element_type), pointer :: shape_h, shape_X
    ! local mass matrix
    real, dimension(ele_loc(h,ele),ele_loc(h,ele)) :: mass_mat

    ele_h=>ele_nodes(h, ele)
    shape_h=>ele_shape(h, ele)
    shape_X=>ele_shape(positions, ele)

    ! Locations of local vertices.
    X_ele=ele_val(positions, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, detwei=detwei)

    mass_mat = shape_shape(shape_h,shape_h,detwei)

    call addto(mass,ele_h,ele_h,mass_mat)

  end subroutine assemble_mass_elemental

  subroutine assemble_mass_elemental_vec(ele,positions,h,mass)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: h
    type(dynamic_csr_matrix), intent(inout) :: mass

    ! Locations of nodes.
    real, dimension(positions%dim,ele_loc(positions,ele)) :: X_ele
    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,ele)) :: detwei    
    ! Node numbers of field element.
    integer, dimension(:), pointer :: ele_h
    ! Shape functions.
    type(element_type), pointer :: shape_h, shape_X
    ! local mass matrix
    real, dimension(ele_loc(h,ele),ele_loc(h,ele)) :: mass_mat

    ele_h=>ele_nodes(h, ele)
    shape_h=>ele_shape(h, ele)
    shape_X=>ele_shape(positions, ele)

    ! Locations of local vertices.
    X_ele=ele_val(positions, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, detwei=detwei)

    mass_mat = shape_shape(shape_h,shape_h,detwei)

    call addto(mass,ele_h,ele_h,mass_mat)

  end subroutine assemble_mass_elemental_vec

  subroutine usage
    
    write (0,*) "usage: IGW -igw_filename <data_file_name>"
    
  end subroutine usage

  subroutine calculate_stream_function_2d(streamfunc,X,U, stat)
    !!< Calculate the stream function for a 
    type(scalar_field), intent(inout) :: streamfunc
    type(vector_field), intent(in) :: X, U
    integer, intent(out), optional :: stat
    
    integer :: i, lstat, ele
    type(csr_sparsity) :: psi_sparsity
    type(csr_matrix) :: psi_mat
    type(scalar_field) :: rhs

    assert(X%dim==2)
    ! No discontinuous stream functions.
    assert(continuity(streamfunc)>=0)
    
    psi_sparsity = make_sparsity(streamfunc%mesh, streamfunc%mesh, &
         "StreamFunctionSparsity")
    
    call allocate(psi_mat, psi_sparsity, name="StreamFunctionMatrix")

    call zero(psi_mat)
    call allocate(rhs, streamfunc%mesh, "StreamFunctionRHS")
    call zero(rhs)

    do ele=1, element_count(streamfunc)
       
       call calculate_streamfunc_ele(psi_mat, rhs, ele, X, U)

    end do

    ewrite(2,*) 'max(abs(rhs))', maxval(abs(rhs%val))
    streamfunc%val = 0.
    !call set(psi_mat,1,1,INFINITY)
    call petsc_solve(streamfunc%val,psi_mat,rhs%val, &
         startfromzero=.true., checkconvergence=.true., &
         iterations=number_of_iterations)
    !call petsc_solve(streamfunc, psi_mat, rhs)

    call deallocate(rhs)
    call deallocate(psi_mat)
    call deallocate(psi_sparsity)
  end subroutine calculate_stream_function_2d
 
  subroutine calculate_streamfunc_ele(psi_mat, rhs, ele, X, U)
    type(csr_matrix), intent(inout) :: psi_mat
    type(scalar_field), intent(inout) :: rhs
    type(vector_field), intent(in) :: X,U
    integer, intent(in) :: ele
    
    ! Transformed gradient function for velocity.
    real, dimension(ele_loc(U, ele), ele_ngi(U, ele), mesh_dim(U)) :: du_t
    ! Ditto for the stream function, psi
    real, dimension(ele_loc(rhs, ele), ele_ngi(rhs, ele), mesh_dim(rhs))&
         & :: dpsi_t 
    
    ! Local vorticity_matrix
    real, dimension(2, ele_loc(rhs, ele), ele_loc(U, ele)) ::&
         & lvorticity_mat
    ! Local vorticity
    real, dimension(ele_loc(rhs, ele)) :: lvorticity
    
    ! Variable transform times quadrature weights.
    real, dimension(ele_ngi(U,ele)) :: detwei
    
    type(element_type), pointer :: U_shape, psi_shape
    integer, dimension(:), pointer :: psi_ele, neigh
    integer :: i, ni, face
    
    U_shape=> ele_shape(U, ele)
    psi_shape=> ele_shape(rhs, ele)
    psi_ele=>ele_nodes(rhs, ele)
    
    ! Transform U derivatives and weights into physical space.
    call transform_to_physical(positions, ele, &
         & U_shape , dshape=du_t, detwei=detwei)
    ! Ditto psi.
    call transform_to_physical(positions,ele, &
         & psi_shape, dshape=dpsi_t)
    
    
    call addto(psi_mat, psi_ele, psi_ele, &
         -dshape_dot_dshape(dpsi_t, dpsi_t, detwei))
    
    lvorticity_mat=shape_curl_shape_2d(psi_shape, du_t, detwei)
    
    lvorticity=0.0
    do i=1,2
       lvorticity=lvorticity &
            +matmul(lvorticity_mat(i,:,:), ele_val(U, ele, i))
    end do
    
    call addto(rhs, psi_ele, lvorticity)
    
    neigh=>ele_neigh(U, ele)

    if(.true.) then
       neighbourloop: do ni=1,size(neigh)
          !Find boundaries.
          if (neigh(ni)<=0) then
             
             face=ele_face(rhs, ele, neigh(ni))
             
             !Strong dirichlet condition (currently the only thing supported)
             call addto_diag(psi_mat, &
                  face_global_nodes(rhs, face), &
                  spread(INFINITY, 1, face_loc(rhs,face)))
          end if
       end do neighbourloop
    end if

  end subroutine calculate_streamfunc_ele

  subroutine get_max(h_surface,x_surface)
    implicit none
    type(scalar_field), intent(in) :: h_surface
    type(vector_field), intent(in) :: X_surface
    !
    integer :: ele
    real, dimension(2,ele_loc(x_surface,1)) :: X_ele
    real, dimension(ele_loc(h_surface,1)) :: h_ele
    real :: maxh, maxh_x(2), maxh_x_left(2), maxh_x_right(2), x0, h0
    real :: maxh_left,maxh_right, maxh_mid

    assert(size(X_ele,2).eq.2)
    assert(size(h_ele).eq.3)

    ewrite(1,*) 'Computing maximum on the surface'

    maxh = 0.

    do ele = 1, element_count(h_surface)
       
       X_ele=ele_val(x_surface, ele)
       h_ele=ele_val(h_surface, ele)

       x0 = (3*h_ele(1) - 4*h_ele(2) + h_ele(3))/ &
            (4*h_ele(1) - 8*h_ele(2) + 4*h_ele(3))

       if(x0.ge.0 .and. x0.le.1) then
          h0 = 0
          h0 = h_ele(1)*(x0**2 - 3*x0+1) + &
               h_ele(2)*(-4*x0**2 + 4*x0) + &
               h_ele(3)*(2*x0**2 - x0)
          if(h0>maxh) then
             maxh = h0
             maxh_x = x_ele(:,1)*(1-x0) + x_ele(:,2)*x0
             maxh_x_left = x_ele(:,1)
             maxh_x_right = x_ele(:,2)
             maxh_left = h_ele(1)
             maxh_mid = h_ele(2)
             maxh_right = h_ele(3)
          end if
          if(h_ele(1)>maxh) then
             maxh = h_ele(1)
             maxh_x = x_ele(:,1)
             maxh_x_left = x_ele(:,1)
             maxh_x_right = x_ele(:,1) 
             maxh_left = h_ele(1)
             maxh_mid = h_ele(1)
             maxh_right = h_ele(1)
         end if
          if(h_ele(3)>maxh) then
             maxh = h_ele(3)
             maxh_x = x_ele(:,2)
             maxh_x_left = x_ele(:,2)
             maxh_x_right = x_ele(:,2)
             maxh_left = h_ele(3)
             maxh_mid = h_ele(3)
             maxh_right = h_ele(3)
          end if
       end if
    end do

    !ewrite(2,*) 'maximum value is', maxh, '>', maxh_left, maxh_mid, maxh_right
    !ewrite(2,*) 'at coordinate', maxh_x_left, '<', maxh_x,'<',maxh_x_right

    ewrite(2,*) 'maximum value on surface is', maxh
    ewrite(2,*) 'at coordinate', maxh_x

  end subroutine get_max

  subroutine get_wave_position(h,positions)
    type(scalar_field), intent(in) :: h
    type(vector_field), intent(in) :: positions
    !locals
    real, dimension(2) :: X0
    real :: h_int
    integer :: ele

    h_int = 0.
    X0 = 0.
    do ele = 1, element_count(h)
       call get_wave_position_elemental(X0,h_int,h,positions,ele)
    end do

    ewrite(2,*) 'wave position integral = ', X0/h_int
    ewrite(2,*) 'n elements = ', element_count(h)

  end subroutine get_wave_position

  subroutine get_wave_position_elemental(X0,h_int,h,positions,ele)
    real, intent(inout), dimension(2) :: X0
    real, intent(inout) :: h_int
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: h
    
    ! Locations of nodes.
    real, dimension(positions%dim,ele_loc(positions,ele)) :: X_ele
    ! Locations of quadrature points
    real, dimension(positions%dim,ele_ngi(positions,ele)) :: X_quad
    real, dimension(ele_ngi(h,ele)) :: h_quad
    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,ele)) :: detwei    
    ! Node numbers of field element.
    integer, dimension(:), pointer :: ele_h
    ! Shape functions.
    type(element_type), pointer :: shape_h, shape_X
    ! local mass matrix
    real, dimension(ele_loc(h,ele),ele_loc(h,ele)) :: mass_mat

    ele_h=>ele_nodes(h, ele)
    shape_h=>ele_shape(h, ele)
    shape_X=>ele_shape(positions, ele)

    ! Locations of local vertices.
    X_ele=ele_val(positions, ele)

    ! Locations of quadrature points.
    X_quad=ele_val_at_quad(positions, ele)
    h_quad=ele_val_at_quad(h, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, detwei=detwei)

    X0 = X0 + sum(shape_vector_rhs(shape_h,X_quad,detwei*h_quad),2)
    h_int = h_int + sum(shape_rhs(shape_h,detwei*h_quad))

  end subroutine get_wave_position_elemental

  subroutine get_L2_estimates(u,h,exact_u,exact_h,positions)
    type(scalar_field), intent(in) :: h,exact_h
    type(vector_field), intent(in) :: positions, u, exact_u
    !locals
    real :: h_int, u_int
    integer :: ele

    h_int = 0.
    u_int = 0.
    do ele = 1, element_count(h)
       call get_L2_estimates_elemental(h_int,u_int,u,h,exact_u,exact_h, &
            positions,ele)
    end do

    ewrite(2,*) 'u integral = ', u_int
    ewrite(2,*) 'h integral = ', h_int

  end subroutine get_L2_estimates

    subroutine get_L2_estimates_elemental(h_int,u_int,u,h,exact_u,exact_h, &
         positions,ele)
      real, intent(inout) :: h_int,u_int
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: positions,u,exact_u
      type(scalar_field), intent(in) :: h,exact_h

      ! Locations of nodes.
      real, dimension(positions%dim,ele_loc(positions,ele)) :: X_ele
      ! Locations of quadrature points and values of u at quads
      real, dimension(positions%dim,ele_ngi(positions,ele)) :: X_quad, u_quad
      ! values of h at quad points
      real, dimension(ele_ngi(h,ele)) :: h_quad
      ! Coordinate transform * quadrature weights.
      real, dimension(ele_ngi(positions,ele)) :: detwei    
      ! Node numbers of field element.
      integer, dimension(:), pointer :: ele_h
      ! Shape functions.
      type(element_type), pointer :: shape_h, shape_X
      ! local mass matrix
      real, dimension(ele_loc(h,ele),ele_loc(h,ele)) :: mass_mat

    ele_h=>ele_nodes(h, ele)
    shape_h=>ele_shape(h, ele)
    shape_X=>ele_shape(positions, ele)

    ! Locations of local vertices.
    X_ele=ele_val(positions, ele)

    ! Locations of quadrature points.
    X_quad=ele_val_at_quad(positions, ele)
    U_quad=ele_val_at_quad(u, ele)-ele_val_at_quad(exact_u, ele)
    h_quad=ele_val_at_quad(h, ele)-ele_val_at_quad(exact_h, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions,ele, detwei=detwei)

    h_int = h_int + sum(detwei*h_quad**2)
    u_int = u_int + sum(detwei*(u_quad(1,:)**2+u_quad(2,:))**2)

  end subroutine get_L2_estimates_elemental

  subroutine dump_matrices(u_mass,h_mass,C1T,C2T,h,u,positions)
    type(csr_matrix), intent(in) :: U_mass, h_mass, C1T, C2T
    type(scalar_field) , intent(in) :: h
    type(vector_field), intent(in) :: u, positions
    !!
    integer :: ele
    integer, dimension(:), pointer :: local_nodes
    integer :: nod

    ewrite(0,*) 'DUMPING OUT ELEMENT-NODE LIST FOR POSITIONS'
    do ele = 1, element_count(positions)
       ewrite(0,*) ele
       ewrite(0,*) ele_val(positions, ele)
       local_nodes => ele_nodes(positions, ele)
       ewrite(0,*) local_nodes
    end do

    ewrite(0,*) 'DUMPING OUT ELEMENT-NODE LIST FOR U'
    do ele = 1, element_count(positions)
       ewrite(0,*) ele
       local_nodes => ele_nodes(u, ele)
       ewrite(0,*) local_nodes
    end do

    ewrite(0,*) 'DUMPING OUT ELEMENT-NODE LIST FOR H'
    do ele = 1, element_count(positions)
       ewrite(0,*) ele
       local_nodes => ele_nodes(h, ele)
       ewrite(0,*) local_nodes
    end do

    ewrite(0,*) 'DUMPING OUT U_MASS'
    ewrite(0,*) dense(u_mass)

    ewrite(0,*) 'DUMPING OUT H_MASS'
    ewrite(0,*) dense(h_mass)

    ewrite(0,*) 'DUMPING OUT C1T'
    ewrite(0,*) dense(C1T)

    ewrite(0,*) 'DUMPING OUT C2T'
    ewrite(0,*) dense(C2T)

  end subroutine dump_matrices

end subroutine IGW
