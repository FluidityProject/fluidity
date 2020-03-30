#include "fdebug.h"

module meshmovement

  use iso_c_binding
  use global_parameters
  use fldebug
  use vector_tools, only: solve, mat_diag_mat, eigendecomposition_symmetric
  use element_numbering
  use elements
  use shape_functions
  use spud
  use sparse_tools
  use fields_base
  use global_numbering
  use eventcounter
  use fetools
  use unittest_tools
  use parallel_fields
  use fields
  use profiler
  use state_module
  use boundary_conditions, only : get_boundary_condition, get_boundary_condition_count,&
       extract_surface_field
  use vtk_interfaces
  use sparse_matrices_fields
  use solvers
  use fefields
  use field_derivatives
  use sparsity_patterns
  use sparse_tools_petsc
  use sparsity_patterns_meshes

  implicit none
  integer,save :: MeshCount=0

  interface

     subroutine set_debug_level(level)
       implicit none
       integer, intent(in) :: level
     end subroutine set_debug_level

     subroutine reset_debug_level
     end subroutine reset_debug_level
  end interface

  interface
     subroutine lap_smoother(dim, num_nodes,&
          num_elements, num_surf_elements, num_owned_nodes,& 
          mapping, connectivity, phys_mesh, smooth_mesh,&
          comp_mesh, surf_connectivity, matrix,&
          options, debug_level) bind(c)
       use iso_c_binding  
#ifdef HAVE_PETSC_MODULES
       use petsc
#endif
       use solvers, only: solver_options_type
       implicit none
#include "petsc_legacy.h"
       
       integer (c_int), value :: dim, num_nodes, num_elements, &
            num_surf_elements, num_owned_nodes, debug_level
       integer (c_int), dimension(num_nodes) :: mapping
       integer (c_int), dimension((dim+1)*num_elements) :: connectivity
       real (c_double), dimension(num_nodes) :: phys_mesh, smooth_mesh, comp_mesh
       integer (c_int), dimension(dim*num_elements) :: surf_connectivity
       type(c_ptr), value :: matrix
       type(solver_options_type) :: options
       
     end subroutine lap_smoother
  end interface

  interface
     subroutine lin_smoother(dim, num_nodes,&
          num_elements, num_surf_elements, num_owned_nodes,& 
          mapping, connectivity, phys_mesh, smooth_mesh,&
          comp_mesh, surf_connectivity,findrm, colm, matrix,&
          options, debug_level) bind(c)
       use iso_c_binding
#ifdef HAVE_PETSC_MODULES
       use petsc
#endif
       use solvers, only: solver_options_type
       implicit none
#include "petsc_legacy.h"
       
       integer (c_int), value :: dim, num_nodes, num_elements, &
            num_surf_elements, num_owned_nodes, debug_level
       integer (c_int), dimension(num_nodes,dim) :: mapping
       integer (c_int), dimension((dim+1)*num_elements) :: connectivity
       real (c_double), dimension(num_nodes) :: phys_mesh, smooth_mesh, comp_mesh
       integer (c_int), dimension(dim*num_elements) :: surf_connectivity
       integer (c_int), dimension(num_nodes+1) :: findrm
       integer (c_int), dimension(findrm(num_nodes+1)-1) :: colm
       type(c_ptr), value :: matrix
       type(solver_options_type) :: options
       
     end subroutine lin_smoother
  end interface

   interface
    subroutine lin_tor_smoother(dim, num_nodes,&
          num_elements, num_surf_elements, num_owned_nodes,& 
          mapping, connectivity, phys_mesh, smooth_mesh,&
          comp_mesh, surf_connectivity,findrm, colm, matrix,&
          options, debug_level) bind(c)
       use iso_c_binding
#ifdef HAVE_PETSC_MODULES
       use petsc
#endif
       use solvers, only: solver_options_type
       implicit none
#include "petsc_legacy.h"
       
       integer (c_int), value :: dim, num_nodes, num_elements, &
            num_surf_elements, num_owned_nodes, debug_level
       integer (c_int), dimension(num_nodes,dim) :: mapping
       integer (c_int), dimension((dim+1)*num_elements) :: connectivity
       real (c_double), dimension(num_nodes) :: phys_mesh, smooth_mesh, comp_mesh
       integer (c_int), dimension(dim*num_elements) :: surf_connectivity
       integer (c_int), dimension(num_nodes+1) :: findrm
       integer (c_int), dimension(findrm(num_nodes+1)-1) :: colm  
       type(c_ptr), value :: matrix
       type(solver_options_type) :: options
     end subroutine lin_tor_smoother
  end interface

    interface
     subroutine cent_relaxer(dim, num_nodes,&
          num_elements, num_surf_elements, num_owned_nodes,& 
          mapping, connectivity, phys_mesh, smooth_mesh,&
          comp_mesh, surf_connectivity,findrm, colm) bind(c)
       use iso_c_binding

       implicit none
       
       integer (c_int), value :: dim, num_nodes, num_elements, &
            num_surf_elements, num_owned_nodes
       integer (c_int), dimension(num_nodes) :: mapping
       integer (c_int), dimension((dim+1)*num_elements) :: connectivity
       real (c_double), dimension(num_nodes) :: phys_mesh, smooth_mesh, comp_mesh
       integer (c_int), dimension(dim*num_elements) :: surf_connectivity
       integer (c_int), dimension(num_nodes+1) :: findrm
       integer (c_int), dimension(findrm(num_nodes+1)-1) :: colm
       
       
     end subroutine cent_relaxer
  end interface
  private
  
  public :: move_mesh_imposed_velocity, move_mesh_pseudo_lagrangian, &
       move_mesh_laplacian_smoothing, move_mesh_initialise_laplacian_smoothing,&
       move_mesh_lineal_smoothing, move_mesh_initialise_lineal_smoothing,&
       move_mesh_lineal_torsional_smoothing, move_mesh_initialise_lineal_torsional_smoothing,&
       move_mesh_centroid_relaxer,move_mesh_initialise_centroid_relaxer,&
       reset_mesh_positions

contains

  subroutine move_mesh_imposed_velocity(states)
  type(state_type), dimension(:), intent(inout) :: states
  
    type(vector_field), pointer :: coordinate, old_coordinate, new_coordinate
    type(vector_field), pointer :: velocity
    type(vector_field), pointer :: grid_velocity
    
    integer :: i, stat
    real :: itheta, dt
    logical :: found_velocity
    
    if(.not.have_option("/mesh_adaptivity/mesh_movement/imposed_grid_velocity")) return
    call IncrementEventCounter(EVENT_MESH_MOVEMENT)
    
    ewrite(1,*) 'Entering move_mesh_imposed_velocity'
    
    grid_velocity => extract_vector_field(states(1), "GridVelocity")
    
    coordinate => extract_vector_field(states(1), "Coordinate")
    old_coordinate => extract_vector_field(states(1), "OldCoordinate")
    new_coordinate => extract_vector_field(states(1), "IteratedCoordinate")
    
    call get_option("/timestepping/timestep", dt)
    
    found_velocity = .false.
    do i = 1, size(states)
      velocity => extract_vector_field(states(i), "Velocity", stat)
      if(stat==0 .and. .not. velocity%aliased) then
        call get_option(trim(velocity%option_path)//"/prognostic/temporal_discretisation/relaxation", itheta, stat)
        if(found_velocity.and.(stat==0)) then
          FLExit("Only one prognostic velocity allowed with imposed mesh movement.")
        else
          found_velocity = (stat==0)
        end if
      end if
    end do
    if(.not.found_velocity) then
      itheta = 0.5
    end if
    
    call set(new_coordinate, old_coordinate)

    if (have_option("/mesh_adaptivity/mesh_movement/transform_coordinates")) then
       call update(new_coordinate, grid_velocity, scale=dt)
    else
       call addto(new_coordinate, grid_velocity, scale=dt)
    end if
    
    call set(coordinate, new_coordinate, old_coordinate, itheta)

  end subroutine move_mesh_imposed_velocity

  subroutine update(X,u,scale)
    type(vector_field), intent(inout) :: X
    type(vector_field), intent(in)    :: u
    real, intent(in)                  :: scale

    integer :: node
    real    :: r, theta, pos(mesh_dim(X)), v(mesh_dim(X)), pos_new(mesh_dim(X))
    real, parameter :: epsilon=1.0e-16

    real :: axis(3), ihat(mesh_dim(X)), jhat(mesh_dim(X))
    real :: C(mesh_dim(X)),p1, p2, v1, v2, delta_r, delta_theta, delta_z
    integer:: dim

    dim=mesh_dim(x)

    ewrite(1,*) 'Transforming to cylindrical coordinates'
    
    C=0.0
    if (have_option("/mesh_adaptivity/mesh_movement/transform_coordinates/point_on_axis"))&
         call get_option("/mesh_adaptivity/mesh_movement/transform_coordinates/point_on_axis",&
         C)

    if (dim==1) then
       FLAbort("Trying to use cylindrical coordinates for 1D mesh movement. This isn't going to work.")
    else if (dim==2) then
       axis=[0.,0.,1.]
       ihat=[1.,0.]
       jhat=[0.,1.]
    else
       call get_option("/mesh_adaptivity/mesh_movement/transform_coordinates/axis_of_rotation",axis,default=[0.0,0.0,1.0])
       assert(sum(axis*axis)>0.0)
       axis=axis/sqrt(sum(axis*axis))
       if (abs(axis(3))>1.0e-16) then
          ihat=[axis(3),0.0,-axis(1)]
          ihat=ihat/sqrt(sum(ihat*ihat))
          jhat=[-axis(1)*axis(2),axis(3)**2+axis(1)**2,-axis(2)*axis(3)]
       else
          ihat=[-axis(2),axis(1),0.0]
          jhat=[-axis(1)*axis(3),-axis(2)*axis(3),axis(1)**2+axis(2)**2]
       end if
       ihat=ihat/sqrt(sum(ihat*ihat))
       jhat=jhat/sqrt(sum(jhat*jhat))
    end if       

    do node=1, node_count(X)
       pos=node_val(X,node)-C
       v=node_val(u,node)

       p1=dot_product(pos,ihat)
       p2=dot_product(pos,jhat)
       v1=dot_product(v,ihat)
       v2=dot_product(v,jhat)

       r=sqrt(sum((pos-dot_product(pos,axis(1:dim))*axis(1:dim))**2))
       theta=atan2(p1,p2)

       delta_r=r+dot_product(pos,v)*scale/(r+epsilon)
       delta_theta=theta+(v1*p2-v2*p1)*scale/(r+epsilon)**2
       delta_z=dot_product(pos+v*scale,axis(1:dim))

       pos_new=C+ihat*delta_r*sin(delta_theta)&
                +jhat*delta_r*cos(delta_theta)&
                +delta_z*axis(1:dim)

       call set(X,node,pos_new)
    end do


  end subroutine update


  subroutine move_mesh_pseudo_lagrangian(states)
  type(state_type), dimension(:), intent(inout) :: states
  
    type(vector_field), pointer :: coordinate, old_coordinate, new_coordinate
    type(vector_field), pointer :: velocity
    type(vector_field), pointer :: grid_velocity
    
    integer :: i, stat
    real :: itheta, dt
    logical :: found_velocity
    
    character(len=FIELD_NAME_LEN) :: state_name
    
    if(.not.have_option("/mesh_adaptivity/mesh_movement/pseudo_lagrangian")) return
    call IncrementEventCounter(EVENT_MESH_MOVEMENT)
    
    ewrite(1,*) 'Entering move_mesh_pseudo_lagrangian'
    
    grid_velocity => extract_vector_field(states(1), "GridVelocity")
    
    call get_option("/mesh_adaptivity/mesh_movement/pseudo_lagrangian/velocity_material_phase/material_phase_name", &
                    state_name, stat=stat)
    if(stat==0) then
      i = get_state_index(states, trim(state_name))
      velocity => extract_vector_field(states(i), "Velocity")
    else
      velocity => extract_vector_field(states(1), "Velocity")
    end if
    
    call set(grid_velocity, velocity)
    
    coordinate => extract_vector_field(states(1), "Coordinate")
    old_coordinate => extract_vector_field(states(1), "OldCoordinate")
    new_coordinate => extract_vector_field(states(1), "IteratedCoordinate")
    
    call get_option("/timestepping/timestep", dt)
    
    found_velocity = .false.
    do i = 1, size(states)
      velocity => extract_vector_field(states(i), "Velocity", stat)
      if(stat==0) then
        call get_option(trim(velocity%option_path)//"/prognostic/temporal_discretisation/relaxation", itheta, stat)
        if(found_velocity.and.(stat==0)) then
          FLExit("Only one prognostic velocity allowed with pseudo lagrangian mesh movement.")
        else
          found_velocity = (stat==0)
        end if
      end if
    end do
    if(.not.found_velocity) then
      itheta = 0.5
    end if
    
    call set(new_coordinate, old_coordinate)
    call addto(new_coordinate, grid_velocity, scale=dt)
    
    call set(coordinate, new_coordinate, old_coordinate, itheta)
  
  end subroutine move_mesh_pseudo_lagrangian


  subroutine move_mesh_laplacian_smoothing(states, diagnostic_only)
    type(state_type), dimension(:), intent(inout) :: states
    logical, optional, intent(in) :: diagnostic_only
  
    type(vector_field), pointer :: coordinate, old_coordinate, new_coordinate,&
         initial_coordinate
    type(vector_field), pointer :: velocity
    type(vector_field), pointer :: grid_velocity
    type(mesh_type), pointer    :: surface_mesh
    
    integer :: i, stat
    real :: itheta, dt
    logical :: found_velocity
    type(solver_options_type) :: solver_options
    type(petsc_csr_matrix), target :: A

    !!! This routine drives a call to C code which does the actual assembly and
    !!! solution for a Laplacian smoothed grid velocity.

    !!! The boundary conditions from the grid velocity field specified in diamond
    !!! should be satisfied by this solution.

    if (.not. have_option("/mesh_adaptivity/mesh_movement/laplacian_smoothing")) return

    ewrite(1,*) 'Entering move_mesh_laplacian_smoothing'
    
    grid_velocity => extract_vector_field(states(1), "GridVelocity")
    if (have_option('/mesh_adaptivity/mesh_movement/solver')) then
       call populate_solver_options_struct_from_path(solver_options,&
            '/mesh_adaptivity/mesh_movement/solver')
       solver_options%ksptype=trim(solver_options%ksptype)//c_null_char
       solver_options%pctype=trim(solver_options%pctype)//c_null_char
       solver_options%hypretype=trim(solver_options%hypretype)//c_null_char
    else 
       solver_options%ksptype=KSPGMRES//c_null_char
       solver_options%pctype=PCSOR//c_null_char
       solver_options%rtol=1.0e-7
       solver_options%atol=1.0e-7
       solver_options%max_its=1000
    end if
    
    call set_boundary_values(grid_velocity)
    
    coordinate => extract_vector_field(states(1), "Coordinate")
    old_coordinate => extract_vector_field(states(1), "OldCoordinate")
    new_coordinate => extract_vector_field(states(1), "IteratedCoordinate")
    initial_coordinate => extract_vector_field(states(1), "InitialCoordinate")
    surface_mesh => extract_mesh(states(1),"FullCoordinateSurfaceMesh")
    
    call profiler_tic("laplacian_smoother")

    call get_option("/timestepping/timestep", dt)

    
    found_velocity = .false.
    do i = 1, size(states)
       velocity => extract_vector_field(states(i), "Velocity", stat)
       if(stat==0 .and. .not. velocity%aliased) then
          call get_option(trim(velocity%option_path)//"/prognostic/temporal_discretisation/relaxation", itheta, stat)
          if(found_velocity.and.(stat==0)) then
             FLExit("Only one prognostic velocity allowed with imposed mesh movement.")
          else
             found_velocity = (stat==0)
          end if
       end if
    end do
    if(.not.found_velocity) then
       itheta = 0.5
    end if

    call set(coordinate, old_coordinate) 
    call addto(coordinate, grid_velocity, scale = dt)
    call allocate(A, &
         get_csr_sparsity_firstorder(states, coordinate%mesh, coordinate%mesh),&
         [1,1],"StiffnessMatrix")
    call zero(A)
  
    !!! actual call out to the C code.

!Figure out indexing between C and Fortran
    call lap_smoother(mesh_dim(coordinate), node_count(coordinate),&
          element_count(coordinate), surface_element_count(coordinate),&
          nowned_nodes(coordinate), universal_numbering(coordinate)-1,&
          coordinate%mesh%ndglno, coordinate%val, new_coordinate%val,&
          initial_coordinate%val, surface_mesh%ndglno,&
          c_loc(A), solver_options, debug_level())

    call deallocate(A)



    call halo_update(new_coordinate)

    !!! Convert the new coordinates returned into a grid velocity
    !!! and calculate the coordinate field at the theta time level.

    call set(grid_velocity,new_coordinate)
    call addto(grid_velocity, old_coordinate, scale = -1.0)
    call scale(grid_velocity, 1.0/dt)

    if (present_and_true(diagnostic_only)) then
       call set(coordinate, old_coordinate)
       call set(new_coordinate, old_coordinate)
    else
       call set(coordinate, new_coordinate, old_coordinate, itheta)
    end if

    call profiler_toc("laplacian_smoother")

  end subroutine move_mesh_laplacian_smoothing

  subroutine move_mesh_initialise_laplacian_smoothing(states)
    type(state_type), dimension(:), intent(inout) :: states

    type(vector_field), pointer :: coordinate
    type(vector_field) :: initial_coordinate
    type(mesh_type) :: surface_mesh

    !!! Store the initial mesh for use as the computational geometry for the
    !!  Laplacian mesh smoothing library.

    coordinate => extract_vector_field(states(1), "Coordinate")

    call allocate(initial_coordinate, mesh=coordinate%mesh, dim=coordinate%dim,&
         name="InitialCoordinate")
    call set(initial_coordinate, coordinate)
    call insert(states(1), initial_coordinate, "InitialCoordinate")
    call deallocate(initial_coordinate)

    !!! Add the surface mesh. This needs revisiting for parallel

    call extract_surface_mesh(surface_mesh, coordinate%mesh, "FullCoordinateSurfaceMesh")
    call insert(states(1), surface_mesh, "FullCoordinateSurfaceMesh")
    call deallocate(surface_mesh)

  end subroutine move_mesh_initialise_laplacian_smoothing

  subroutine move_mesh_lineal_smoothing(states, diagnostic_only)
    type(state_type), dimension(:), intent(inout) :: states
    logical, optional, intent(in) :: diagnostic_only
  
    type(vector_field), pointer :: coordinate, old_coordinate, new_coordinate,&
         initial_coordinate
    type(vector_field), pointer :: velocity
    type(vector_field), pointer :: grid_velocity
    type(mesh_type), pointer    :: surface_mesh
    
    integer :: i, stat
    real :: itheta, dt
    logical :: found_velocity
    type(solver_options_type) :: solver_options
    type(petsc_csr_matrix), target :: A

    !!! This routine drives a call to C code which does the actual assembly and
    !!! solution for a Lineal Spring smoothed grid velocity.

    !!! The boundary conditions from the grid velocity field specified in diamond
    !!! should be satisfied by this solution.

    if (.not. have_option("/mesh_adaptivity/mesh_movement/lineal_smoothing")) return

    ewrite(1,*) 'Entering move_mesh_lineal_smoothing'
    
    grid_velocity => extract_vector_field(states(1), "GridVelocity")
    if (have_option('/mesh_adaptivity/mesh_movement/solver')) then
       call populate_solver_options_struct_from_path(solver_options,&
            '/mesh_adaptivity/mesh_movement/solver')
       solver_options%ksptype=trim(solver_options%ksptype)//c_null_char
       solver_options%pctype=trim(solver_options%pctype)//c_null_char
       solver_options%hypretype=trim(solver_options%hypretype)//c_null_char
    else 
       solver_options%ksptype=KSPGMRES//c_null_char
       solver_options%pctype=PCSOR//c_null_char
       solver_options%rtol=1.0e-7
       solver_options%atol=1.0e-7
       solver_options%max_its=1000
    end if

    call set_boundary_values(grid_velocity)
    
    coordinate => extract_vector_field(states(1), "Coordinate")
    old_coordinate => extract_vector_field(states(1), "OldCoordinate")
    new_coordinate => extract_vector_field(states(1), "IteratedCoordinate")
    initial_coordinate => extract_vector_field(states(1), "InitialCoordinate")
    surface_mesh => extract_mesh(states(1),"FullCoordinateSurfaceMesh")
    
    call profiler_tic(coordinate, "lineal_smoother")

    call get_option("/timestepping/timestep", dt)

    
    found_velocity = .false.
    do i = 1, size(states)
       velocity => extract_vector_field(states(i), "Velocity", stat)
       if(stat==0 .and. .not. velocity%aliased) then
          call get_option(trim(velocity%option_path)//"/prognostic/temporal_discretisation/relaxation", itheta, stat)
          if(found_velocity.and.(stat==0)) then
             FLExit("Only one prognostic velocity allowed with imposed mesh movement.")
          else
             found_velocity = (stat==0)
          end if
       end if
    end do
    if(.not.found_velocity) then
       itheta = 0.5
    end if
       
    call set(coordinate, old_coordinate) 
    call addto(coordinate, grid_velocity, scale = dt)
    call allocate(A, &
         get_csr_sparsity_firstorder(states, coordinate%mesh, coordinate%mesh),&
         [mesh_dim(coordinate),mesh_dim(coordinate)],"StiffnessMatrix")
    call zero(A)

    !!! actual call out to the C code.


    if (.not. associated(coordinate%mesh%adj_lists%nnlist)) &
         call add_nnlist(coordinate%mesh)
     
!Figure out indexing between C and Fortran
    call lin_smoother(mesh_dim(coordinate), node_count(coordinate),&
          element_count(coordinate), surface_element_count(coordinate),&
          nowned_nodes(coordinate), A%column_numbering%gnn2unn,&
          coordinate%mesh%ndglno, coordinate%val, new_coordinate%val,&
          initial_coordinate%val, surface_mesh%ndglno,&
          coordinate%mesh%adj_lists%nnlist%findrm,&
          coordinate%mesh%adj_lists%nnlist%colm, c_loc(A),&
          solver_options, debug_level())
    

    call halo_update(new_coordinate)

    !!! Convert the new coordinates returned into a grid velocity
    !!! and calculate the coordinate field at the theta time level.

    call set(grid_velocity,new_coordinate)
    call addto(grid_velocity, old_coordinate, scale = -1.0)
    call scale(grid_velocity, 1.0/dt)

    if (present_and_true(diagnostic_only)) then
       call set(coordinate, old_coordinate)
       call set(new_coordinate, old_coordinate)
    else
       call set(coordinate, new_coordinate, old_coordinate, itheta)
    end if
    call deallocate(A)

    call profiler_toc(coordinate, "lineal_smoother")

  end subroutine move_mesh_lineal_smoothing

  subroutine move_mesh_initialise_lineal_smoothing(states)
    type(state_type), dimension(:), intent(inout) :: states

    type(vector_field), pointer :: coordinate
    type(vector_field) :: initial_coordinate
    type(mesh_type) :: surface_mesh

    !!! Store the initial mesh for use as the computational geometry for the
    !!  Lineal mesh smoothing library.

    coordinate => extract_vector_field(states(1), "Coordinate")

    call allocate(initial_coordinate, mesh=coordinate%mesh, dim=coordinate%dim,&
         name="InitialCoordinate")
    call set(initial_coordinate, coordinate)
    call insert(states(1), initial_coordinate, "InitialCoordinate")
    call deallocate(initial_coordinate)

    !!! Add the surface mesh. This needs revisiting for parallel

    call extract_surface_mesh(surface_mesh, coordinate%mesh, "FullCoordinateSurfaceMesh")
    call insert(states(1), surface_mesh, "FullCoordinateSurfaceMesh")
    call deallocate(surface_mesh)

  end subroutine move_mesh_initialise_lineal_smoothing

  subroutine move_mesh_lineal_torsional_smoothing(states, diagnostic_only)
    type(state_type), dimension(:), intent(inout) :: states
    logical, optional, intent(in) :: diagnostic_only
  
    type(vector_field), pointer :: coordinate, old_coordinate, new_coordinate,&
         initial_coordinate
    type(vector_field), pointer :: velocity
    type(vector_field), pointer :: grid_velocity
    type(mesh_type), pointer    :: surface_mesh
    
    integer :: i, stat
    real :: itheta, dt
    logical :: found_velocity
	type(solver_options_type) :: solver_options
    type(petsc_csr_matrix), target :: A

    !!! This routine drives a call to C code which does the actual assembly and
    !!! solution for a Lineal Torsional Spring smoothed grid velocity.

    !!! The boundary conditions from the grid velocity field specified in diamond
    !!! should be satisfied by this solution.

    if (.not. have_option("/mesh_adaptivity/mesh_movement/lineal_torsional_smoothing")) return

    ewrite(1,*) 'Entering move_mesh_lineal_torsional_smoothing'
    
    grid_velocity => extract_vector_field(states(1), "GridVelocity")
    if (have_option('/mesh_adaptivity/mesh_movement/solver')) then
       call populate_solver_options_struct_from_path(solver_options,&
            '/mesh_adaptivity/mesh_movement/solver')
       solver_options%ksptype=trim(solver_options%ksptype)//c_null_char
       solver_options%pctype=trim(solver_options%pctype)//c_null_char
       solver_options%hypretype=trim(solver_options%hypretype)//c_null_char
    else 
       solver_options%ksptype=KSPGMRES//c_null_char
       solver_options%pctype=PCSOR//c_null_char
       solver_options%rtol=1.0e-7
       solver_options%atol=1.0e-7
       solver_options%max_its=1000
    end if

    call set_boundary_values(grid_velocity)
    
    coordinate => extract_vector_field(states(1), "Coordinate")
    old_coordinate => extract_vector_field(states(1), "OldCoordinate")
    new_coordinate => extract_vector_field(states(1), "IteratedCoordinate")
    initial_coordinate => extract_vector_field(states(1), "InitialCoordinate")
    surface_mesh => extract_mesh(states(1),"FullCoordinateSurfaceMesh")
    
    call profiler_tic(coordinate, "lineal_torsional_smoother")

    call get_option("/timestepping/timestep", dt)

    
    found_velocity = .false.
    do i = 1, size(states)
       velocity => extract_vector_field(states(i), "Velocity", stat)
       if(stat==0 .and. .not. velocity%aliased) then
          call get_option(trim(velocity%option_path)//"/prognostic/temporal_discretisation/relaxation", itheta, stat)
          if(found_velocity.and.(stat==0)) then
             FLExit("Only one prognostic velocity allowed with imposed mesh movement.")
          else
             found_velocity = (stat==0)
          end if
       end if
    end do
    if(.not.found_velocity) then
       itheta = 0.5
    end if
       
    call set(coordinate, old_coordinate) 
    call addto(coordinate, grid_velocity, scale = dt)   
    call allocate(A, &
         get_csr_sparsity_firstorder(states, coordinate%mesh, coordinate%mesh),&
         [mesh_dim(coordinate),mesh_dim(coordinate)],"StiffnessMatrix")
    call zero(A)    
  
    !!! actual call out to the C code.
    if (.not. associated(coordinate%mesh%adj_lists%nnlist)) &
         call add_nnlist(coordinate%mesh)
     
!Figure out indexing between C and Fortran
    call lin_tor_smoother(mesh_dim(coordinate), node_count(coordinate),&
          element_count(coordinate), surface_element_count(coordinate),&
          nowned_nodes(coordinate), A%column_numbering%gnn2unn,&
          coordinate%mesh%ndglno, coordinate%val, new_coordinate%val,&
          initial_coordinate%val, surface_mesh%ndglno,&
          coordinate%mesh%adj_lists%nnlist%findrm,&
          coordinate%mesh%adj_lists%nnlist%colm, c_loc(A),&
          solver_options, debug_level())

    call halo_update(new_coordinate)

    !!! Convert the new coordinates returned into a grid velocity
    !!! and calculate the coordinate field at the theta time level.

    call set(grid_velocity,new_coordinate)
    call addto(grid_velocity, old_coordinate, scale = -1.0)
    call scale(grid_velocity, 1.0/dt)

    if (present_and_true(diagnostic_only)) then
       call set(coordinate, old_coordinate)
       call set(new_coordinate, old_coordinate)
    else
       call set(coordinate, new_coordinate, old_coordinate, itheta)
    end if
    call deallocate(A)

    call profiler_toc(coordinate, "lineal_torsional_smoother")

  end subroutine move_mesh_lineal_torsional_smoothing

  
  subroutine move_mesh_initialise_lineal_torsional_smoothing(states)
    type(state_type), dimension(:), intent(inout) :: states

    type(vector_field), pointer :: coordinate
    type(vector_field) :: initial_coordinate
    type(mesh_type) :: surface_mesh

    !!! Store the initial mesh for use as the computational geometry for the
    !!  centroid relaxation library.

    coordinate => extract_vector_field(states(1), "Coordinate")

    call allocate(initial_coordinate, mesh=coordinate%mesh, dim=coordinate%dim,&
         name="InitialCoordinate")
    call set(initial_coordinate, coordinate)
    call insert(states(1), initial_coordinate, "InitialCoordinate")
    call deallocate(initial_coordinate)

    !!! Add the surface mesh. This needs revisiting for parallel

    call extract_surface_mesh(surface_mesh, coordinate%mesh, "FullCoordinateSurfaceMesh")
    call insert(states(1), surface_mesh, "FullCoordinateSurfaceMesh")
    call deallocate(surface_mesh)

  end subroutine move_mesh_initialise_lineal_torsional_smoothing

  subroutine move_mesh_centroid_relaxer(states, diagnostic_only)
    type(state_type), dimension(:), intent(inout) :: states
    logical, optional, intent(in) :: diagnostic_only
  
    type(vector_field), pointer :: coordinate, old_coordinate, new_coordinate,&
         initial_coordinate
    type(vector_field), pointer :: velocity
    type(vector_field), pointer :: grid_velocity
    type(mesh_type), pointer    :: surface_mesh
    
    integer :: i, stat
    real :: itheta, dt
    logical :: found_velocity

    !!! This routine drives a call to C code which does the actual assembly and
    !!! solution for a Lineal Spring smoothed grid velocity.

    !!! The boundary conditions from the grid velocity field specified in diamond
    !!! should be satisfied by this solution.

    if (.not. have_option("/mesh_adaptivity/mesh_movement/centroid_relaxer")) return

    ewrite(1,*) 'Entering move_mesh_centroid_relaxer'
    
    grid_velocity => extract_vector_field(states(1), "GridVelocity")

    call set_boundary_values(grid_velocity)
    
    coordinate => extract_vector_field(states(1), "Coordinate")
    old_coordinate => extract_vector_field(states(1), "OldCoordinate")
    new_coordinate => extract_vector_field(states(1), "IteratedCoordinate")
    initial_coordinate => extract_vector_field(states(1), "InitialCoordinate")
    surface_mesh => extract_mesh(states(1),"FullCoordinateSurfaceMesh")
    
    call get_option("/timestepping/timestep", dt)

    
    found_velocity = .false.
    do i = 1, size(states)
       velocity => extract_vector_field(states(i), "Velocity", stat)
       if(stat==0 .and. .not. velocity%aliased) then
          call get_option(trim(velocity%option_path)//"/prognostic/temporal_discretisation/relaxation", itheta, stat)
          if(found_velocity.and.(stat==0)) then
             FLExit("Only one prognostic velocity allowed with imposed mesh movement.")
          else
             found_velocity = (stat==0)
          end if
       end if
    end do
    if(.not.found_velocity) then
       itheta = 0.5
    end if
       
    call set(coordinate, old_coordinate) 
    call addto(coordinate, grid_velocity, scale = dt)       
  
    !!! actual call out to the C code.


    if (.not. associated(coordinate%mesh%adj_lists%nnlist)) &
         call add_nnlist(coordinate%mesh)
     
!Figure out indexing between C and Fortran
    call cent_relaxer(mesh_dim(coordinate), node_count(coordinate),&
          element_count(coordinate), surface_element_count(coordinate),&
          nowned_nodes(coordinate), universal_numbering(coordinate)-1,&
          coordinate%mesh%ndglno, coordinate%val, new_coordinate%val,&
          initial_coordinate%val, surface_mesh%ndglno,&
          coordinate%mesh%adj_lists%nnlist%findrm,&
          coordinate%mesh%adj_lists%nnlist%colm)
    

    call halo_update(new_coordinate)

    !!! Convert the new coordinates returned into a grid velocity
    !!! and calculate the coordinate field at the theta time level.

    call set(grid_velocity,new_coordinate)
    call addto(grid_velocity, old_coordinate, scale = -1.0)
    call scale(grid_velocity, 1.0/dt)

    if (present_and_true(diagnostic_only)) then
       call set(coordinate, old_coordinate)
       call set(new_coordinate, old_coordinate)
    else
       call set(coordinate, new_coordinate, old_coordinate, itheta)
    end if


  end subroutine move_mesh_centroid_relaxer

  subroutine move_mesh_initialise_centroid_relaxer(states)
    type(state_type), dimension(:), intent(inout) :: states

    type(vector_field), pointer :: coordinate
    type(vector_field) :: initial_coordinate
    type(mesh_type) :: surface_mesh

    !!! Store the initial mesh for use as the computational geometry for the
    !!  Lineal mesh smoothing library.

    coordinate => extract_vector_field(states(1), "Coordinate")

    call allocate(initial_coordinate, mesh=coordinate%mesh, dim=coordinate%dim,&
         name="InitialCoordinate")
    call set(initial_coordinate, coordinate)
    call insert(states(1), initial_coordinate, "InitialCoordinate")
    call deallocate(initial_coordinate)

    !!! Add the surface mesh. This needs revisiting for parallel

    call extract_surface_mesh(surface_mesh, coordinate%mesh, "FullCoordinateSurfaceMesh")
    call insert(states(1), surface_mesh, "FullCoordinateSurfaceMesh")
    call deallocate(surface_mesh)

  end subroutine move_mesh_initialise_centroid_relaxer

  subroutine set_boundary_values(vfield)
    type(vector_field), intent(inout) :: vfield

    type(vector_field), pointer:: surface_field
    integer, dimension(:), pointer :: surface_node_list

    integer :: i, node

    do i=1,get_boundary_condition_count(vfield)

       call get_boundary_condition(vfield, i, surface_node_list= surface_node_list)
       surface_field => extract_surface_field(vfield,i,"value")

       do node=1, node_count(surface_field)
          call set(vfield, surface_node_list(node), node_val(surface_field,node))
       end do

    end do

  end subroutine set_boundary_values

  subroutine reset_mesh_positions(states)
    type(state_type), dimension(:), intent(inout) :: states
    type(vector_field), pointer :: coordinate, &
         old_coordinate, new_coordinate

    coordinate => extract_vector_field(states(1), "Coordinate")
    old_coordinate => extract_vector_field(states(1), "OldCoordinate")
    new_coordinate => extract_vector_field(states(1), "IteratedCoordinate")

    call set(new_coordinate, old_coordinate)
    call set(coordinate, old_coordinate)

    !! Should probably technically replace the velocity with
    !! velocity - new_grid_velocity + old_grid_velocity
    !! but that mostly falls out in the wash.

  end subroutine reset_mesh_positions


end module meshmovement


