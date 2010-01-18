!    Copyright (C) 2007 Imperial College London and others.
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
module populate_state_module
  use elements
  use state_module
  use FLDebug
  use spud
  use read_triangle
  use vtk_cache_module
  use global_parameters, only: OPTION_PATH_LEN, is_active_process, pi
  use field_options
  use reserve_state_module
  use fields_manipulation
  use diagnostic_variables, only: convergence_field, steady_state_field
  use field_options
  use surfacelabels
  use climatology
  use metric_tools
  use coordinates
  use halos
  use tictoc
  use hadapt_extrude
  use initialise_fields_module
  use transform_elements
  use parallel_tools
  use boundary_conditions_from_options
  use nemo_states_module

  implicit none

  private

  public populate_state
  public populate_state_module_check_options
  public insert_external_mesh, insert_derived_meshes, &
       allocate_field_as_constant, allocate_and_insert_fields, &
       initialise_prognostic_fields, set_prescribed_field_values, &
       alias_fields, mesh_name, &
       allocate_and_insert_auxilliary_fields, &
       initialise_field, allocate_metric_limits

  interface allocate_field_as_constant
    
    module procedure allocate_field_as_constant_path, &
          allocate_field_as_constant_scalar, allocate_field_as_constant_vector, &
          allocate_field_as_constant_tensor

  end interface allocate_field_as_constant
    
  !! A list of locations in which additional scalar fields
  !! are to be found. It is assumed that all additional scalar fields are
  !! in state 1.
  character(len=OPTION_PATH_LEN), dimension(5) :: field_locations=&
       (/ &
       "/ocean_biology/pznd                                                     ", &
       "/material_phase[0]/subgridscale_parameterisations/Mellor_Yamada         ", &
       "/material_phase[0]/subgridscale_parameterisations/prescribed_diffusivity", &
       "/material_phase[0]/subgridscale_parameterisations/GLS                   ", &
       "/porous_media                                                           " &
       /)

  !! Relative paths under a field that are searched for grandchildren
  !! (moved here because of extremely obscure intel ICE -Stephan)
  character(len=OPTION_PATH_LEN), dimension(1):: &
         grandchild_paths = (/&
         &    "/spatial_discretisation/inner_element" &
         /)

contains


  subroutine populate_state(states)
    type(state_type), pointer, dimension(:) :: states

    integer :: nstates ! number of states
    integer :: i

    ewrite(1,*) "In populate_state"
    call tictoc_clear(TICTOC_ID_IO_READ)

    ! Find out how many states there are
    nstates=option_count("/material_phase")
    allocate(states(1:nstates))
    do i = 1, nstates
       call nullify(states(i))
    end do

    call insert_external_mesh(states, save_vtk_cache = .true.)

    call insert_derived_meshes(states)

    call allocate_and_insert_fields(states)

    call initialise_prognostic_fields(states, save_vtk_cache=.true., &
      initial_mesh=.true.)

    call set_prescribed_field_values(states, initial_mesh=.true.)

    call populate_boundary_conditions(states)

    call set_boundary_conditions_values(states)

    call set_dirichlet_consistent(states)

    call alias_fields(states)

    call create_reserve_state(states)

    call tictoc_report(2, TICTOC_ID_IO_READ)
    ewrite(1, *) "Exiting populate_state"

  end subroutine populate_state

  subroutine insert_external_mesh(states, save_vtk_cache)
    !!< Read in external meshes from file as specified in options tree and
    !!< insert in  state
    type(state_type), intent(inout), dimension(:) :: states
    !! By default the vtk_cache, build up by the vtu mesh reads in this
    !! subroutine, is flushed at the end of this subroutine. This cache can be
    !! reused however in subsequent calls reading from vtu files.
    logical, intent(in), optional:: save_vtk_cache

    type(mesh_type) :: mesh
    type(vector_field) :: position
    type(vector_field), pointer :: position_ptr
    character(len=OPTION_PATH_LEN) :: mesh_path, mesh_file_name,&
         & mesh_file_format
    integer, dimension(:), pointer :: coplanar_ids
    integer :: i, j, nmeshes, nstates, quad_degree, stat
    type(element_type), pointer :: shape
    type(quadrature_type), pointer :: quad
    integer :: dim, loc
    integer :: quad_family
    
    call tic(TICTOC_ID_IO_READ)

    ! Find out how many states there are
    nstates=option_count("/material_phase")
    ! Get number of meshes
    nmeshes=option_count("/geometry/mesh")
    ewrite(2,*) "There are", nmeshes, "meshes."

    external_mesh_loop: do i=0, nmeshes-1

       ! Save mesh path
       mesh_path="/geometry/mesh["//int2str(i)//"]"

       if(have_option(trim(mesh_path)//"/from_file")) then

          ! Get file format
          ! Can remove stat test when mesh format data backwards compatibility is removed
          call get_option(trim(mesh_path)//"/from_file/format/name", mesh_file_format, stat)
          ! Can remove following when mesh format data backwards compatibility is removed
          if(stat /= 0) then
             ewrite(0, *) "Warning: Mesh format name attribute missing for mesh " // trim(mesh_path)
             call get_option(trim(mesh_path)//"/from_file/format", mesh_file_format)
          end if

          ! Get filename for mesh, and other options
          call get_option(trim(mesh_path)//"/from_file/file_name", mesh_file_name)
          call get_option("/geometry/quadrature/degree", quad_degree)
          quad_family = get_quad_family()

          if (.not. is_active_process) then
            ! is_active_process records whether we have data on disk or not
            ! see the comment in Global_Parameters. In this block, 
            ! we want to allocate an empty mesh and positions.
            assert(trim(mesh_file_format)=="triangle")

            allocate(quad)
            allocate(shape)
            call identify_triangle_file(trim(mesh_file_name) // "_0", dim, loc)
            quad = make_quadrature(loc, dim, degree=quad_degree, family=quad_family)
            shape=make_element_shape(loc, dim, 1, quad)
            call allocate(mesh, nodes=0, elements=0, shape=shape, name="EmptyMesh")
            call allocate(position, dim, mesh, "EmptyCoordinate")          
            call add_faces(mesh)

            ! Reference counting cleanups.
            call deallocate(mesh)
            call deallocate(quad)
            call deallocate(shape)
            
            deallocate(quad)
            deallocate(shape)

          else if(trim(mesh_file_format)=="triangle") then
             ! Read mesh from triangle file
             position=read_triangle_files(trim(mesh_file_name), quad_degree=quad_degree, quad_family=quad_family)
             mesh=position%mesh
          else if(trim(mesh_file_format) == "vtu") then
             position_ptr => vtk_cache_read_positions_field(mesh_file_name)
             ! No hybrid mesh support here
             assert(ele_count(position_ptr) > 0)
             dim = position_ptr%dim
             loc = ele_loc(position_ptr, 1)

             ! Generate a copy, and swap the quadrature degree
             ! Note: Even if positions_ptr has the correct quadrature degree, it
             ! won't have any faces and hence a copy is still required (as
             ! add_faces is a construction routine only)
             allocate(quad)
             allocate(shape)
             quad = make_quadrature(loc, dim, degree = quad_degree, family=quad_family)
             shape = make_element_shape(loc, dim, 1, quad)
             call allocate(mesh, nodes = node_count(position_ptr), elements = ele_count(position_ptr), shape = shape, name = position_ptr%mesh%name)
             do j = 1, ele_count(mesh)
               call set_ele_nodes(mesh, j, ele_nodes(position_ptr%mesh, j))
             end do
             call add_faces(mesh)
             call allocate(position, dim, mesh, position_ptr%name)
             call set(position, position_ptr)
             call deallocate(mesh)
             call deallocate(shape)
             call deallocate(quad)
             deallocate(quad)
             deallocate(shape)

             mesh = position%mesh
          else
             ewrite(-1,*) trim(mesh_file_format), " is not a valid format for a mesh file"
             FLAbort("Invalid format for mesh file")
          end if

          ! Get mesh name. This must be done after the triangle file has
          ! been read otherwise the filename is automatically inserted
          ! as the mesh name.
          call get_option(trim(mesh_path)//"/name", mesh%name)
          
          ! Set mesh option path.
          mesh%option_path = mesh_path
          
          ! Copy those changes back to the descriptor under position%mesh
          position%mesh=mesh

          if (mesh%name/="CoordinateMesh") then
             position%name=trim(mesh%name)//"Coordinate"
          else 
             position%name="Coordinate"
          end if
                       
          ! If running in parallel, additionally read in halo information and register the elements halo
          if(isparallel()) then
             call read_halos(mesh_file_name, position%mesh)
             call reorder_element_numbering(position)
             mesh = position%mesh
          end if
          
          ! coplanar ids are create here already and stored on the mesh, 
          ! so its derived meshes get the same coplanar ids
          ! (must be done after halo registration)
          call get_coplanar_ids(mesh, position, coplanar_ids)
          
          ! Insert mesh and position field into states(1) and
          ! alias it to all the others
          call insert(states, mesh, mesh%name)
          call insert(states, position, position%name)
         
          mesh%option_path=mesh_path
          
          call surface_id_stats(mesh, position)

          call deallocate(position)
       end if

    end do external_mesh_loop
    
    if(.not. present_and_true(save_vtk_cache)) then
       ! Flush the cache
       call vtk_cache_finalise()
    end if

    call toc(TICTOC_ID_IO_READ)

  end subroutine insert_external_mesh

  subroutine insert_derived_meshes(states)
    ! Insert derived meshes in state
    type(state_type), intent(inout), dimension(:) :: states

    type(element_type):: shape
    type(quadrature_type) :: quad
    type(mesh_type) :: mesh, model_mesh, from_mesh, periodic_mesh
    type(vector_field), pointer :: position, modelposition
    type(vector_field) :: coordinateposition, from_position, extrudedposition
    
    character(len=FIELD_NAME_LEN) :: mesh_name
    character(len=OPTION_PATH_LEN) :: mesh_path, model_mesh_name, from_mesh_name, bc_name
    character(len=OPTION_PATH_LEN) :: continuity_option
    character(len=OPTION_PATH_LEN) :: periodic_mapping_python
    ! logicals to find out if we have certain options
    logical :: new_shape, new_cont, extrusion
    integer, dimension(:), allocatable :: physical_boundary_ids, aliased_boundary_ids
    integer, dimension(2) :: shape_option
    integer loc, dim, quad_degree, continuity, poly_degree
    integer i, j, nmeshes, stat, n_periodic_bcs
    logical :: incomplete, updated

    ! Get number of meshes
    nmeshes=option_count("/geometry/mesh")

    outer_loop: do
       ! Updated becomes true if we manage to set up at least one mesh on
       ! this pass.
       updated=.false.
       ! Incomplete becomes true if we have to skip over at least one mesh
       ! on this pass.
       incomplete=.false.

       derived_mesh_loop: do i=0, nmeshes-1
          
          ! Save mesh path
          mesh_path="/geometry/mesh["//int2str(i)//"]"

          ! Get mesh name.
          call get_option(trim(mesh_path)//"/name", mesh_name)

          if (has_mesh(states(1), mesh_name)) then
             ! We already did this one.
             cycle derived_mesh_loop
          end if
          
          if(have_option(trim(mesh_path)//"/from_mesh")) then
             
             ! Get model mesh name
             call get_option(trim(mesh_path)//"/from_mesh/mesh[0]/name", model_mesh_name)
             
             ! Find out if the new mesh is different from the old mesh and if
             ! so, find out how it differs
             new_shape=have_option(trim(mesh_path)//"/from_mesh/mesh_shape/polynomial_degree")
             new_cont=have_option(trim(mesh_path)//"/from_mesh/mesh_continuity")
             extrusion=have_option(trim(mesh_path)//"/from_mesh/extrude")
             
             ! Extract model mesh
             model_mesh=extract_mesh(states(1), trim(model_mesh_name), stat=stat)
             if (stat/=0) then
                ! The mesh from which this mesh is derived is not yet
                ! present.
                incomplete=.true.
                cycle derived_mesh_loop
             end if
             ! We added at least one mesh on this pass.
             updated=.true.
                
             if(.not. new_shape .and. .not. new_cont .and. .not. extrusion) then
                
                ! copy mesh unchanged
                mesh=model_mesh
                
                call incref(mesh)
                
                mesh%name = mesh_name
                
                ! Set mesh option path.
                mesh%option_path = trim(mesh_path)
                
             else if (extrusion) then
                
                ! see if adaptivity has left us something:
                extrudedposition=extract_vector_field(states(1), &
                     "AdaptedExtrudedPositions", stat=stat)
                     
                if (stat==0) then

                   ! extrusion has already done by adaptivity
                  
                   ! we remove them here, as we want to insert them under different names
                   call incref(extrudedposition)
                   do j=1, size(states)
                      call remove_vector_field(states(j), "AdaptedExtrudedPositions")
                   end do
                   
                else
                   
                   ! extrusion by user specifed layer depths
                   
                   ! get the model positions from the external mesh... this assumes that we are only one
                   ! mesh removed from the external mesh... dangerous!
                   modelposition => extract_vector_field(states(1), trim(model_mesh_name)//"Coordinate")
                
                   call extrude(modelposition, mesh_path, extrudedposition)
                   
                end if
                
                mesh = extrudedposition%mesh
                mesh%name = mesh_name
                mesh%option_path = trim(mesh_path)
                
                ! the positions of this mesh have to be stored now 
                ! as it cannot be interpolated later.
                if (mesh_name=="CoordinateMesh") then
                   extrudedposition%name = "Coordinate"
                else
                   extrudedposition%name = trim(mesh_name)//"Coordinate"
                end if
                call insert(states, extrudedposition, extrudedposition%name)
                call deallocate(extrudedposition)                
                
                call incref(mesh)
                
             else
                
                ! Get new mesh shape information
                ! degree is the degree of the Lagrange polynomials
                call get_option(trim(mesh_path)//"/from_mesh/mesh_shape/polynomial_degree", poly_degree, stat)
                if(stat==0) then
                   ! loc is the number of vertices of the element
                   loc=model_mesh%shape%loc
                   ! dim is the dimension
                   dim=model_mesh%shape%dim
                   ! Make quadrature
                   call get_option("/geometry/quadrature/degree",&
                        & quad_degree)
                   quad=make_quadrature(loc, dim, degree=quad_degree, family=get_quad_family())
                   ! Make new mesh shape
                   shape=make_element_shape(loc, dim, poly_degree, quad)
                   call deallocate(quad) ! Really just drop a reference.
                else
                   shape=model_mesh%shape
                   call incref(shape)
                end if
                
                ! Get new mesh continuity
                call get_option(trim(mesh_path)//"/from_mesh/mesh_continuity", continuity_option, stat)
                if(stat==0) then
                   if(trim(continuity_option)=="discontinuous") then
                      continuity=-1
                   else if(trim(continuity_option)=="continuous") then
                      continuity=0
                   end if
                else
                   continuity=model_mesh%continuity
                end if
                
                ! Get mesh name.
                call get_option(trim(mesh_path)//"/name", mesh_name)
                
                ! Make new mesh
                mesh=make_mesh(model_mesh, shape, continuity, mesh_name)
                
                ! Set mesh option path
                mesh%option_path = trim(mesh_path)
                
                ! Drop one reference to shape
                call deallocate(shape)
             end if
             
             ! if this is the coordinate mesh then we should insert the coordinate field
             if (trim(mesh_name)=="CoordinateMesh" .and. .not. extrusion) then
                ! get the model positions from the external mesh... this assumes that we are only one
                ! mesh removed from the external mesh... dangerous!
                modelposition => extract_vector_field(states(1), trim(model_mesh_name)//"Coordinate")
                
                call allocate(coordinateposition, modelposition%dim, mesh, "Coordinate")
                ! remap the external mesh positions onto the CoordinateMesh... this requires that the space
                ! of the coordinates spans that of the external mesh
                call remap_field(from_field=modelposition, to_field=coordinateposition)
                ! insert into states(1) and alias to all others
                call insert(states, coordinateposition, coordinateposition%name)
                ! drop reference to the local copy of the Coordinate field
                call deallocate(coordinateposition)
             end if

             if (have_option(trim(mesh_path)//"/from_mesh/periodic_boundary_conditions")) then
                
                call get_option(trim(mesh_path)//"/from_mesh/mesh[0]/name", from_mesh_name)
                from_mesh = extract_mesh(states(1), trim(from_mesh_name), stat=stat)
                if (stat/=0) then
                   ! The mesh from which this mesh is derived is not yet
                   ! present.
                   incomplete=.true.
                   cycle derived_mesh_loop
                end if
                ! We added at least one mesh on this pass.
                updated=.true.

                position => extract_vector_field(states(1), "Coordinate")
                
                n_periodic_bcs=option_count(trim(mesh_path)//"/from_mesh/periodic_boundary_conditions")
                ewrite(2,*) "n_periodic_bcs=", n_periodic_bcs
                do j=0, n_periodic_bcs-1
                   
                   ! get some options
                   call get_option(trim(mesh_path)//"/from_mesh/periodic_boundary_conditions["//int2str(j)//"]/name", bc_name)
                   ewrite(1,*) "applying boundary condition: ", trim(bc_name)
                   shape_option = option_shape(trim(mesh_path)//"/from_mesh/periodic_boundary_conditions["//int2str(j)//"]/physical_boundary_ids")
                   allocate( physical_boundary_ids(shape_option(1)) )
                   call get_option(trim(mesh_path)//"/from_mesh/periodic_boundary_conditions["//int2str(j)//"]/physical_boundary_ids",physical_boundary_ids)
                   shape_option = option_shape(trim(mesh_path)//"/from_mesh/periodic_boundary_conditions["//int2str(j)//"]/aliased_boundary_ids")
                   allocate( aliased_boundary_ids(shape_option(1)) )
                   call get_option(trim(mesh_path)//"/from_mesh/periodic_boundary_conditions["//int2str(j)//"]/aliased_boundary_ids",aliased_boundary_ids)
                   call get_option(trim(mesh_path)//"/from_mesh/periodic_boundary_conditions["//int2str(j)//"]/coordinate_map",periodic_mapping_python)
                   
                   ewrite(2,*) 'Making periodic mesh'
                   
                   if (from_mesh==position%mesh) then
                     from_position=position
                     call incref(from_position)
                   else
                     call allocate(from_position, position%dim, from_mesh, "ModelPositions")
                     call remap_field(position, from_position)
                   end if
                   periodic_mesh=make_mesh_periodic(from_mesh,from_position,&
                      physical_boundary_ids,aliased_boundary_ids, &
                      periodic_mapping_python)
                   call deallocate(from_position)
                   from_mesh=periodic_mesh
                   
                   deallocate( physical_boundary_ids, aliased_boundary_ids )
                end do

                periodic_mesh%name = mesh%name
                periodic_mesh%option_path= mesh_path
                call deallocate(mesh)
                mesh=periodic_mesh
                
             end if
             
             ! Insert mesh into all states
             call insert(states, mesh, mesh%name)
             
             call deallocate(mesh)
             
          end if
          
       end do derived_mesh_loop

       ! If we didn't skip any fields then we are done.
       if (.not.incomplete) exit outer_loop

       ! If we did skip fields and didn't update any fields this pass, then
       ! we have unresolvable dependencies.
       if (.not.updated) then
          FLExit("Unresolvable mesh dependencies")
       end if
       
    end do outer_loop

  end subroutine insert_derived_meshes

  subroutine allocate_and_insert_fields(states, dont_allocate_prognostic_value_spaces)
    !!< allocates and inserts all fields present in the options tree
    !!< zeros field, but does not yet set initial conditions
    type(state_type), dimension(:), intent(inout):: states
    !! If provided and true will not allocate a full value space
    !! for those fields for which defer_allocation(option_path, mesh) is .true.
    !! but instead allocate them as constant fields. This is used
    !! for fields that are passed down to SAM in which case we want to be 
    !! able to one by one allocate them as we get them back from SAM.
    logical, optional, intent(in):: dont_allocate_prognostic_value_spaces

    character(len=OPTION_PATH_LEN) :: field_name
    integer :: i, ii ! counters
    integer :: nstates ! number of states
    integer :: ncars   ! number of vehicles
    character(len=255) :: tmp ! temporary string to make life a little easier

    nstates=option_count("/material_phase")

    ! Loop over states for the first time to get prognostic, prescribed and diagnostic fields.
    state_loop: do i=0, nstates-1

       ! Assign the material_phase name to state(i+1)%name
       call get_option('/material_phase['//int2str(i)//']/name', states(i+1)%name)

       call allocate_and_insert_one_phase(&
            '/material_phase['//int2str(i)//']', states(i+1), &
            dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
       
    end do state_loop

    ! special case fields outside material_phases:
    ! distance to top and bottom
    if (have_option('/geometry/ocean_boundaries')) then
       ! set up DistanceToTop field and insert in first state
       ! it is only allowed to be diagnostic by the schema, so not much to do
       call allocate_and_insert_scalar_field('/geometry/ocean_boundaries/scalar_field::DistanceToTop', &
            states(1))

       ! set up DistanceToBottom field and insert in first state
       call allocate_and_insert_scalar_field('/geometry/ocean_boundaries/scalar_field::DistanceToBottom', &
            states(1), &
            dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)

    end if

    ! direction of gravity
    if (have_option('/physical_parameters/gravity/vector_field::GravityDirection')) then
       call allocate_and_insert_vector_field('/physical_parameters/gravity/vector_field::GravityDirection', &
            states(1), parent_mesh='CoordinateMesh', backward_compatibility=.false., &
            dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    end if


    ! insert traffic tracer fields
    if (have_option('/traffic_model/scalar_field::TrafficTracerTemplate')) then
       call get_option('/traffic_model/number_of_vehicles', ncars)
       do i=1, ncars
          field_name='TrafficTracer'//int2str(i)
          call allocate_and_insert_scalar_field('/traffic_model/scalar_field::TrafficTracerTemplate', &
            states(1), parent_mesh='VelocityMesh', backward_compatibility=.true., &
            field_name=field_name, &
            dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
       end do
    end if

    ! insert porous media fields
    if (have_option('/porous_media')) then
       do i=1, nstates
          call allocate_and_insert_scalar_field('/porous_media/scalar_field::Porosity', &
             states(i), parent_mesh='VelocityMesh', field_name='Porosity')
          if (have_option("/porous_media/scalar_field::Permeability")) then
             call allocate_and_insert_scalar_field('/porous_media/scalar_field::Permeability', &
               states(i), parent_mesh='VelocityMesh', field_name='Permeability')
          elseif (have_option("/porous_media/vector_field::Permeability")) then
             call allocate_and_insert_vector_field('/porous_media/vector_field::Permeability', &
               states(i), parent_mesh='VelocityMesh')
          elseif (have_option("/porous_media/tensor_field::Permeability")) then
             call allocate_and_insert_tensor_field('/porous_media/tensor_field::Permeability', &
               states(i), parent_mesh='VelocityMesh')
          end if
       end do
    end if

    ! insert electrical property fields
    do i=1,nstates
      tmp = '/material_phase['//int2str(i-1)//']/electrical_properties/coupling_coefficients/'
      ! Electrokinetic coupling coefficient scalar field
      if (have_option(trim(tmp)//'scalar_field::Electrokinetic')) then
        call allocate_and_insert_scalar_field(trim(tmp)//'scalar_field::Electrokinetic', &
                                              states(i), &
                                              field_name='Electrokinetic['//int2str(i-1)//']')
      end if
      ! Thermoelectric coupling coefficient scalar field
      if (have_option(trim(tmp)//'scalar_field::Thermoelectric')) then
        call allocate_and_insert_scalar_field(trim(tmp)//'scalar_field::Thermoelectric', &
                                              states(i), &
                                              field_name='Thermoelectric['//int2str(i-1)//']')
      end if
      ! Electrochemical coupling coefficient scalar field
      if (have_option(trim(tmp)//'scalar_field::Electrochemical')) then
        call allocate_and_insert_scalar_field(trim(tmp)//'scalar_field::Electrochemical', &
                                              states(i), &
                                              field_name='Electrochemical['//int2str(i-1)//']')
      end if
    end do

    if ( has_scalar_field(states(1),'HarmonicAmplitudeM2') .or. has_scalar_field(states(1),'HarmonicPhaseM2') ) then
       do ii=1,50 ! this number will be an option from diamond, hard code to a reasonable value for now.

! if you give an option path these fields will be output to vtu                 
!          call allocate_and_insert_scalar_field('/material_phase[Fields]/scalar_field::harmonic', &
!             states(1), parent_mesh='VelocityMesh', field_name='harmonic'//int2str(ii))

          call allocate_and_insert_scalar_field('', &
             states(1), parent_mesh='VelocityMesh', field_name='harmonic'//int2str(ii))

       end do
    end if




    ! insert miscellaneous fields
    do i=1, size(field_locations)
       if (have_option(trim(field_locations(i)))) then

          call allocate_and_insert_one_phase(field_locations(i), states(1), &
             dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
          
       end if
    end do

    call allocate_metric_limits(states(1))
    
  contains
    
    subroutine allocate_and_insert_one_phase(state_path, state, dont_allocate_prognostic_value_spaces)
      !! Perform the allocation and insertion of the fields found under
      !! state_path into state.
      character(len=*), intent(in) :: state_path
      type(state_type), intent(inout) :: state
      logical, optional, intent(in):: dont_allocate_prognostic_value_spaces

      character(len=OPTION_PATH_LEN) :: path      
      integer :: nfields ! number of fields
      logical :: is_aliased
      
      integer :: j

      ! Get number of scalar fields that are children of this state
      nfields=option_count(trim(state_path)//"/scalar_field")

      ! Loop over scalar fields
      scalar_field_loop: do j=0, nfields-1

         ! Save path to field
         path=trim(state_path)//"/scalar_field["//int2str(j)//"]"

         ! Get field name
         call get_option(trim(path)//"/name", field_name)
         ! Reset path to have field name rather than index
         path=trim(state_path)//"/scalar_field::"//trim(field_name)

         ! If field is not aliased call allocate_and_insert_scalar_field
         is_aliased=have_option(trim(path)//"/aliased")
         if(.not.is_aliased) then
            call allocate_and_insert_scalar_field(path, state, &
              dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
         end if

      end do scalar_field_loop

      ! Get number of vector fields that are children of this state
      nfields=option_count(trim(state_path)//"/vector_field")

       ! Loop over vector fields
       vector_field_loop: do j=0, nfields-1

          ! Save path to field
          path=trim(state_path)//"/vector_field["//int2str(j)//"]"
          ! Get field name
          call get_option(trim(path)//"/name", field_name)
          ! Reset path to have field name rather than index
          path=trim(state_path)//"/vector_field::"//trim(field_name)

          ! If field is not aliased call allocate_and_insert_vector_field
          is_aliased=have_option(trim(path)//"/aliased")
          if(.not.is_aliased) then
             call allocate_and_insert_vector_field(path, state, &
               dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
          end if

       end do vector_field_loop

       ! Get number of tensor fields that are children of this state
       nfields=option_count(trim(state_path)//"/tensor_field")

       tensor_field_loop: do j=0, nfields-1

          ! Save path to field
          path=trim(state_path)//"/tensor_field["//int2str(j)//"]"
          ! Get field name
          call get_option(trim(path)//"/name", field_name)
          ! Reset path to have field name rather than index
          path=trim(state_path)//"/tensor_field::"//trim(field_name)

          ! If field is not aliased call allocate_and_insert_tensor_field
          is_aliased=have_option(trim(path)//"/aliased")
          if(.not.is_aliased) then
             call allocate_and_insert_tensor_field(path, state, &
                dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
          end if

       end do tensor_field_loop

       ! Sediment submodel
       if (have_option(trim(state_path)//"/sediment")) then 
          call allocate_and_insert_sediment(state_path, state)
       end if

    end subroutine allocate_and_insert_one_phase

    subroutine allocate_and_insert_sediment(state_path, state)
      !! Allocate all the sediment submodel fields.
      character(len=*), intent(in) :: state_path
      type(state_type), intent(inout) :: state
      
      integer :: nfields, j
      character(len=OPTION_PATH_LEN) :: class_path, field_name,class_name
      
      type(scalar_field), pointer :: sediment, sedimentflux

      nfields=option_count(trim(state_path)//"/sediment/sediment_class")

      
          sediment_class_loop: do j=0,nfields-1
             ! Note that this currently duplicates shared subfields such as
             ! diffusivity. This should be changed.

             class_path=trim(state_path)&
                  //"/sediment/sediment_class["//int2str(j)//"]"

             call get_option(trim(class_path)//"/name",class_name)
             field_name="SedimentConcentration"//trim(class_name)
             call allocate_and_insert_scalar_field(&
                  trim(state_path)&
                  //"/sediment/scalar_field::SedimentTemplate", &
                  state, field_name=field_name, &
                  dont_allocate_prognostic_value_spaces&
                  =dont_allocate_prognostic_value_spaces)
             
             sediment=>extract_scalar_field(state, field_name)
             
             ! Now replace the default children with any special ones.
             call allocate_and_insert_children(class_path, state, &
                  sediment%mesh%name, field_name, &
                  dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
             call allocate_and_insert_grandchildren(class_path, state, &
                  sediment%mesh%name, field_name, &
                  & dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)

             ! Now set up the diagnostic flux field.
             field_name="SedimentFlux"//trim(class_name)

             call allocate_and_insert_scalar_field(&
                  trim(state_path)&
                  //"/sediment/scalar_field::SedimentFluxTemplate", &
                  state, field_name=field_name, &
                  dont_allocate_prognostic_value_spaces&
                  =dont_allocate_prognostic_value_spaces)
             
             sedimentflux=>extract_scalar_field(state, field_name)

             call zero(sedimentflux)

          end do sediment_class_loop

    end subroutine allocate_and_insert_sediment

  end subroutine allocate_and_insert_fields

  subroutine alias_fields(states)
    type(state_type), dimension(:), intent(inout) :: states

    character(len=OPTION_PATH_LEN) :: path
    character(len=OPTION_PATH_LEN) :: state_name, aliased_field_name, field_name
    integer :: i, j, k ! counters
    integer :: nstates ! number of states
    integer :: nfields ! number of fields
    ! logicals to find out if we have certain options
    logical :: is_aliased
    type(scalar_field) :: sfield
    type(vector_field) :: vfield
    type(tensor_field) :: tfield

    nstates=option_count("/material_phase")

    state_loop: do i=0, nstates-1

       ! Get number of scalar fields that are children of this state
       nfields=option_count("/material_phase["//int2str(i)//"]/scalar_field")

       ! Loop over scalar fields
       scalar_field_loop: do j=0, nfields-1

          ! Save path to field
          path="/material_phase["//int2str(i)//"]/scalar_field["&
               &//int2str(j)//"]"

          ! If field is aliased, find which field it is aliased to, extract that field from the correct state and insert into current state
          is_aliased=have_option(trim(path)//"/aliased")
          if(is_aliased) then
             call get_option(trim(path)//"/name", field_name)
             call get_option(trim(path)//"/aliased/field_name", aliased_field_name)
             call get_option(trim(path)//"/aliased/material_phase_name", state_name)

             k=get_state_index(states, trim(state_name))
             sfield=extract_scalar_field(states(k), trim(aliased_field_name))
             sfield%name = trim(field_name)  ! this seems to be necessary
             ! to preserve the aliased field's original name
             sfield%aliased = .true.
             call insert(states(i+1), sfield, trim(field_name))
          end if

       end do scalar_field_loop

       ! Get number of vector fields that are children of this state
       nfields=option_count("/material_phase["//int2str(i)//"]/vecto&
            &r_field")

       ! Loop over vector fields
       vector_field_loop: do j=0, nfields-1

          ! Save path to field
          path="/material_phase["//int2str(i)//"]/vector_field["&
               &//int2str(j)//"]"

          ! If field is aliased, find which field it is aliased to, extract that field from the correct state and insert into current state
          is_aliased=have_option(trim(path)//"/aliased")
          if(is_aliased) then
             call get_option(trim(path)//"/name", field_name)
             call get_option(trim(path)//"/aliased/material_phase_name", state_name)
             call get_option(trim(path)//"/aliased/field_name", aliased_field_name)

             k=get_state_index(states, trim(state_name))
             vfield=extract_vector_field(states(k), trim(aliased_field_name))
             vfield%name = trim(field_name)  ! this seems to be necessary to preserve the aliased field's original name
             vfield%aliased = .true.
             call insert(states(i+1), vfield, trim(field_name))

          end if

       end do vector_field_loop

       ! Get number of tensor fields that are children of this state
       nfields=option_count("/material_phase["//int2str(i)//"]/tensor_field")

       tensor_field_loop: do j=0, nfields-1

          ! Save path to field
          path="/material_phase["//int2str(i)//"]/tensor_field["&
               &//int2str(j)//"]"

          ! If field is aliased, find which field it is aliased to, extract that field from the correct state and insert into current state
          is_aliased=have_option(trim(path)//"/aliased")
          if(is_aliased) then
             call get_option(trim(path)//"/name", field_name)
             call get_option(trim(path)//"/aliased/material_phase_name", state_name)
             call get_option(trim(path)//"/aliased/field_name", aliased_field_name)

             k=get_state_index(states, trim(state_name))
             tfield=extract_tensor_field(states(k), trim(aliased_field_name))
             tfield%name = trim(field_name)  ! this seems to be necessary to preserve the aliased field's original name
             tfield%aliased = .true.
             call insert(states(i+1), tfield, trim(field_name))

          end if

       end do tensor_field_loop

    end do state_loop

    ! special case fields outside material_phases:
    ! direction of gravity
    if (have_option('/physical_parameters/gravity/vector_field::GravityDirection')) then
       vfield=extract_vector_field(states(1), 'GravityDirection')
       vfield%aliased = .true.
       do i = 1,nstates-1

          call insert(states(i+1), vfield, 'GravityDirection')

       end do
    end if

    ! Deal with subgridscale parameterisations.
    call alias_diffusivity(states)

  end subroutine alias_fields

  subroutine alias_diffusivity(states)
    !!< Where fields get their diffusivity from a subgridscale
    !!< parameterisation, it is necessary to alias their diffusivity to the
    !!< diffusivity provided by the parameterisation.
    !!<
    !!< At this stage only prescribed diffusivity is handled via this
    !!< route. Gent-McWilliams diffusivity is calculated on the fly and
    !!< Mellor-Yamada is pending a rewrite.
    type(state_type), dimension(:), intent(inout) :: states
    type(scalar_field), pointer :: sfield
    type(tensor_field) :: tfield

    integer :: i, s, stat

    do i = 1, size(states)
       
       tfield=extract_tensor_field(states(i), "PrescribedDiffusivity", stat)

       if (stat/=0) cycle

       tfield%aliased=.True.

       do s = 1, scalar_field_count(states(i))

          sfield => extract_scalar_field(states(i), s)
          
          if (have_option(trim(sfield%option_path)//&
               "/prognostic/subgridscale_parameterisation&
               &::prescribed_diffusivity")) then
             
             tfield%name=trim(sfield%name)//"Diffusivity"
             call insert(states(i), tfield, tfield%name)

          end if

       end do
       
    end do


    do i = 1, size(states)
       
       tfield=extract_tensor_field(states(i), "GLSEddyDiffusivityKH", stat)

       if (stat/=0) cycle

       tfield%aliased=.True.

       do s = 1, scalar_field_count(states(i))

          sfield => extract_scalar_field(states(i), s)
          
          if (have_option(trim(sfield%option_path)//&
               "/prognostic/subgridscale_parameterisation&
               &::GLS")) then
             
             tfield%name=trim(sfield%name)//"Diffusivity"
             call insert(states(i), tfield, tfield%name)

          end if

       end do
       
    end do



  end subroutine alias_diffusivity

  function allocate_field_as_constant_path(option_path) result(is_constant)
    !!< Return whether the supplied option path signals a constant
    !!< field

    character(len = *), intent(in) :: option_path

    logical :: is_constant

    is_constant = .false.
    if(option_count(trim(option_path) // "/prescribed/value") == 1) then
       is_constant = have_option(trim(option_path) // "/prescribed/value[0]/constant")
    end if

  end function allocate_field_as_constant_path

  function allocate_field_as_constant_scalar(s_field) result(is_constant)
    !!< Return whether the options tree defines the supplied scalar field to
    !!< be constant

    type(scalar_field), intent(in) :: s_field

    logical :: is_constant

    is_constant = allocate_field_as_constant(s_field%option_path)

  end function allocate_field_as_constant_scalar

  function allocate_field_as_constant_vector(v_field) result(is_constant)
    !!< Return whether the options tree defines the supplied vector field to
    !!< be constant

    type(vector_field), intent(in) :: v_field

    logical :: is_constant

    is_constant = allocate_field_as_constant(v_field%option_path)

  end function allocate_field_as_constant_vector

  function allocate_field_as_constant_tensor(t_field) result(is_constant)
    !!< Return whether the options tree defines the supplied tensor field to
    !!< be constant

    type(tensor_field), intent(in) :: t_field

    logical :: is_constant

    if(trim(t_field%name) == "MinMetricEigenbound") then
      is_constant = have_option("/mesh_adaptivity/hr_adaptivity/tensor_field::MinimumEdgeLengths/anisotropic_symmetric/constant")
    else if(trim(t_field%name) == "MaxMetricEigenbound") then
      is_constant = have_option("/mesh_adaptivity/hr_adaptivity/tensor_field::MaximumEdgeLengths/anisotropic_symmetric/constant")
    else
      is_constant = .false.
    end if

  end function allocate_field_as_constant_tensor

  recursive subroutine allocate_and_insert_scalar_field(option_path, state, &
    parent_mesh, parent_name, backward_compatibility, field_name, &
    dont_allocate_prognostic_value_spaces)

    character(len=*), intent(in) :: option_path
    type(state_type), intent(inout) :: state
    character(len=*), intent(in), optional :: parent_mesh
    character(len=*), intent(in), optional :: parent_name
    logical, intent(in), optional :: backward_compatibility
    character(len=*), optional, intent(in):: field_name
    logical, optional, intent(in):: dont_allocate_prognostic_value_spaces

    logical :: is_prognostic, is_prescribed, is_diagnostic
    ! paths for options and child fields
    character(len=OPTION_PATH_LEN) :: path, adapt_path
    ! Strings for names
    character(len=OPTION_PATH_LEN) :: lfield_name, mesh_name
    type(scalar_field) :: field
    type(mesh_type), pointer :: mesh
    logical :: lbackward_compatibility, is_constant

    ! Save option_path
    path=trim(option_path)

    if (present(field_name)) then
       lfield_name=field_name
    else
       call get_option(trim(path)//"/name", lfield_name)
    end if
    
    if(present(parent_name)) then
       lfield_name=trim(parent_name)//trim(lfield_name)
    end if
    ewrite(1,*) "In allocate_and_insert_scalar_field, field is: ", trim(lfield_name)

    ! Do we need backward compatibility?
    ! If we need backward compatibility, then no matter how the field
    ! is described in XML, a value space will be allocated, for old-style
    ! code to use.
    ! If we do not need backward compatibility, we can make big savings
    ! on constant fields.
    if (present(backward_compatibility)) then
       lbackward_compatibility = backward_compatibility
    else
       lbackward_compatibility = .true.
    end if

    ! Find out what kind of field we have
    is_prognostic=have_option(trim(path)//"/prognostic")
    is_prescribed=have_option(trim(path)//"/prescribed")
    is_diagnostic=have_option(trim(path)//"/diagnostic")

    is_constant=allocate_field_as_constant(path)

    ewrite(1,*) "Is field prognostic? ", is_prognostic
    ewrite(1,*) "Is field prescribed? ", is_prescribed
    ewrite(1,*) "Is field constant? ", is_constant
    ewrite(1,*) "Is field diagnostic? ", is_diagnostic

    if (is_prognostic) then

       path=trim(path)//"/prognostic"

    else if(is_prescribed) then

       path=trim(path)//"/prescribed"

    else if(is_diagnostic) then

       path=trim(path)//"/diagnostic"

    end if

    ! Get mesh
    if(present(parent_mesh).and.&
         .not.have_option(trim(path)//"/mesh[0]/name")) then
       mesh => extract_mesh(state, trim(parent_mesh))
       mesh_name=parent_mesh
    else
       call get_option(trim(path)//"/mesh[0]/name", mesh_name)
       mesh => extract_mesh(state, trim(mesh_name))
    end if

    if (defer_allocation(option_path, mesh, dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)) then
       ! If we want to defer allocation (for sam), don't allocate the value space yet
       call allocate(field, mesh, name=trim(lfield_name), &
          field_type=FIELD_TYPE_DEFERRED)
    else if(is_constant .and. .not. lbackward_compatibility) then
         
       ! Allocate as constant field if possible (and we don't need backward compatibility)
       call allocate(field, mesh, name=trim(lfield_name), &
          field_type=FIELD_TYPE_CONSTANT)
       
       call zero(field)
    else
       ! If we have to keep backward compatibility, then
       ! just allocate the value space as normal,
       ! and don't try any funny tricks to save memory.

       ! Allocate field
       call allocate(field, mesh, name=trim(lfield_name))
       call zero(field)
    end if


    ewrite(2,*) trim(lfield_name), " is on mesh ", trim(mesh%name)

    ! Set field%option_path
    field%option_path=trim(option_path)

    ! Finally! Insert field into state!
    call insert(state, field, field%name)
    call deallocate(field)

    ! Check for fields that are children of this field:
    call allocate_and_insert_children(path, state, mesh_name, lfield_name, &
        dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    call allocate_and_insert_grandchildren(path, state, mesh_name,&
         & lfield_name, &
         & dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)

    ! Check for adaptivity weights associated with this field:
    adapt_path=trim(path)//"/adaptivity_options"
    if(have_option(trim(adapt_path)//"/absolute_measure")) then
       adapt_path=trim(adapt_path)//"/absolute_measure/scalar_field::InterpolationErrorBound"
       call allocate_and_insert_scalar_field(adapt_path, state, parent_mesh=mesh_name, &
          parent_name=lfield_name, backward_compatibility=.false., &
          dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    else if(have_option(trim(adapt_path)//"/relative_measure")) then
       adapt_path=trim(adapt_path)//"/relative_measure/scalar_field::InterpolationErrorBound"
       call allocate_and_insert_scalar_field(adapt_path, state, parent_mesh=mesh_name, &
          parent_name=lfield_name, backward_compatibility=.false., &
          dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    end if

  end subroutine allocate_and_insert_scalar_field

  recursive subroutine allocate_and_insert_vector_field(option_path, state, parent_mesh, parent_name, backward_compatibility, dont_allocate_prognostic_value_spaces)

    character(len=*), intent(in) :: option_path
    type(state_type), intent(inout) :: state
    character(len=*), intent(in), optional :: parent_mesh
    character(len=*), intent(in), optional :: parent_name
    logical, intent(in), optional :: backward_compatibility
    logical, intent(in), optional :: dont_allocate_prognostic_value_spaces
    
    integer :: dim
    logical :: is_prognostic, is_prescribed, is_diagnostic
    ! paths for options and child fields
    character(len=OPTION_PATH_LEN) :: path, adapt_path
    ! strings for names
    character(len=OPTION_PATH_LEN) :: field_name, mesh_name
    type(mesh_type), pointer :: mesh
    type(vector_field) :: field
    logical :: lbackward_compatibility, is_constant

    ! Save option_path
    path=trim(option_path)

    call get_option(trim(path)//"/name", field_name)
    if(present(parent_name)) then
       field_name=trim(parent_name)//trim(field_name)
    end if
    ewrite(1,*) "In allocate_and_insert_vector_field, field is: ", trim(field_name)

    ! Do we need backward compatibility?
    ! If we need backward compatibility, then no matter how the field
    ! is described in XML, a value space will be allocated, for old-style
    ! code to use.
    ! If we do not need backward compatibility, we can make big savings
    ! on constant fields.
    if (present(backward_compatibility)) then
       lbackward_compatibility = backward_compatibility
    else
       lbackward_compatibility = .true.
    end if

    ! Find out what kind of field we have
    is_prognostic=have_option(trim(path)//"/prognostic")
    is_prescribed=have_option(trim(path)//"/prescribed")
    is_diagnostic=have_option(trim(path)//"/diagnostic")

    is_constant=allocate_field_as_constant(path)

    ewrite(1,*) "Is field prognostic? ", is_prognostic
    ewrite(1,*) "Is field prescribed? ", is_prescribed
    ewrite(1,*) "Is field constant? ", is_constant
    ewrite(1,*) "Is field diagnostic? ", is_diagnostic

    ! Get dimension of vector - currently the dimension of the problem
    call get_option("/geometry/dimension", dim)

    if(is_prognostic) then
       path=trim(path)//"/prognostic"
    else if(is_prescribed) then
       path=trim(path)//"/prescribed"
    else if(is_diagnostic) then
       path=trim(path)//"/diagnostic"
    end if

    ! Get mesh
    if(present(parent_mesh).and.&
         .not.have_option(trim(path)//"/mesh[0]/name")) then
       mesh => extract_mesh(state, trim(parent_mesh))
       mesh_name=parent_mesh
    else
       call get_option(trim(path)//"/mesh[0]/name", mesh_name)
       mesh => extract_mesh(state, trim(mesh_name))
    end if

    if (defer_allocation(option_path, mesh, dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)) then
       ! If we want to defer allocation (for sam), don't allocate the value space yet
       call allocate(field, dim, mesh, name=trim(field_name), &
          field_type=FIELD_TYPE_DEFERRED)
    else if(is_constant .and. .not. lbackward_compatibility) then
         
       ! Allocate as constant field if possible (and we don't need backward compatibility)
       call allocate(field, dim, mesh, name=trim(field_name), &
          field_type=FIELD_TYPE_CONSTANT)
       call zero(field)
       
    else
       ! If we have to keep backward compatibility, then
       ! just allocate the value space as normal,
       ! and don't try any funny tricks to save memory.

       ! Allocate field
       call allocate(field, dim, mesh, trim(field_name))
       call zero(field)
    end if
    
    ewrite(2,*) trim(field_name), " is on mesh ", trim(mesh%name)

    ! Set field%option_path
    field%option_path=trim(option_path)

    ! Finally! Insert field into state!
    call insert(state, field, field%name)
    call deallocate(field)

    ! Check for fields that are children of this field:
    call allocate_and_insert_children(path, state, mesh_name, field_name, &
        dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    call allocate_and_insert_grandchildren(path, state, mesh_name, field_name, &
        dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)

    ! Check for adaptivity weights associated with this field:
    adapt_path=trim(path)//"/adaptivity_options"
    if(have_option(trim(adapt_path)//"/absolute_measure")) then
       adapt_path=trim(adapt_path)//"/absolute_measure/vector_field::InterpolationErrorBound"
       call allocate_and_insert_vector_field(adapt_path, state, mesh_name, field_name, backward_compatibility=.false., &
          dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    else if(have_option(trim(adapt_path)//"/relative_measure")) then
       adapt_path=trim(adapt_path)//"/relative_measure/vector_field::InterpolationErrorBound"
       call allocate_and_insert_vector_field(adapt_path, state, mesh_name, field_name, backward_compatibility=.false., &
          dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    end if

  end subroutine allocate_and_insert_vector_field

  recursive subroutine allocate_and_insert_tensor_field(option_path, state, parent_mesh, parent_name, &
      dont_allocate_prognostic_value_spaces)
    !!< This subroutine sets up the initial condition of a tensor field.
    !!< Note that the tensor dimensions are set to be the dimension of the 
    !!< problem.

    character(len=*), intent(in) :: option_path
    type(state_type), intent(inout) :: state
    character(len=*), intent(in), optional :: parent_mesh
    character(len=*), intent(in), optional :: parent_name
    logical, intent(in), optional :: dont_allocate_prognostic_value_spaces

    logical :: is_prescribed, is_diagnostic, is_constant
    ! paths for options and child fields
    character(len=OPTION_PATH_LEN) :: path, adapt_path
    character(len=OPTION_PATH_LEN) :: field_name, mesh_name
    type(tensor_field) :: field
    type(mesh_type), pointer:: mesh

    ! Save option_path
    path=trim(option_path)

    call get_option(trim(path)//"/name", field_name)
    if(present(parent_name)) then
       if(trim(field_name)/="Viscosity".and.trim(field_name)/="Elasticity") then
          field_name=trim(parent_name)//trim(field_name)
       end if
    end if
    ewrite(1,*) "In allocate_and_insert_tensor_field, field is: ", trim(field_name)

    ! Find out what kind of field we have
    is_prescribed=have_option(trim(path)//"/prescribed")
    is_diagnostic=have_option(trim(path)//"/diagnostic")
    is_constant=allocate_field_as_constant(path)

    ewrite(1,*) "Is field prescribed? ", is_prescribed
    ewrite(1,*) "Is field diagnostic? ", is_diagnostic
    ewrite(1,*) "Is field constant? ", is_constant

    if(is_prescribed) then

       path=trim(path)//"/prescribed"

    else if(is_diagnostic) then

       path=trim(path)//"/diagnostic"

    end if

    ! Get mesh
    if(present(parent_mesh).and.&
         .not.have_option(trim(path)//"/mesh[0]/name")) then
       mesh => extract_mesh(state, trim(parent_mesh))
       mesh_name=parent_mesh
    else
       call get_option(trim(path)//"/mesh[0]/name", mesh_name)
       mesh => extract_mesh(state, trim(mesh_name))
    end if

    if (defer_allocation(option_path, mesh, dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)) then
         
       ! If we want to defer allocation (for sam), don't allocate the value space yet
       call allocate(field, mesh, name=trim(field_name), &
          field_type=FIELD_TYPE_DEFERRED)
       
    else if(is_constant) then
         
       ! Allocate as constant field if possible (and we don't need backward compatibility)
       call allocate(field, mesh, name=trim(field_name), &
          field_type=FIELD_TYPE_CONSTANT)
       call zero(field)
    else

       ! Allocate field
       call allocate(field, mesh, trim(field_name))
       call zero(field)
    end if


    ! Set field%option_path
    field%option_path=trim(option_path)

    ! Finally! Insert field into state!
    call insert(state, field, field%name)
    call deallocate(field)

    ! Check for fields that are children of this field:
    call allocate_and_insert_children(path, state, mesh_name, field_name, &
         & dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    call allocate_and_insert_grandchildren(path, state, mesh_name, &
         & field_name, dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)

    ! Check for adaptivity weights associated with this field:
    adapt_path=trim(path)//"/adaptivity_options"
    if(have_option(trim(adapt_path)//"/absolute_measure")) then
       adapt_path=trim(adapt_path)//"/absolute_measure/tensor_field::InterpolationErrorBound"
       call allocate_and_insert_tensor_field(adapt_path, state, mesh_name, field_name, &
            dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    else if(have_option(trim(adapt_path)//"/relative_measure")) then
       adapt_path=trim(adapt_path)//"/relative_measure/tensor_field::InterpolationErrorBound"
       call allocate_and_insert_tensor_field(adapt_path, state, mesh_name, field_name, &
            dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    end if

  end subroutine allocate_and_insert_tensor_field

  subroutine allocate_and_insert_children(path, state, parent_mesh, parent_name, &
       dont_allocate_prognostic_value_spaces)
    character(len=*), intent(in) :: path !! option_path including prescribed/prognostic
    type(state_type), intent(inout) :: state
    character(len=*), intent(in) :: parent_mesh
    character(len=*), intent(in) :: parent_name
    logical, optional, intent(in) :: dont_allocate_prognostic_value_spaces

    character(len=OPTION_PATH_LEN) child_path, child_name
    character(len=FIELD_NAME_LEN) :: mesh_name
    integer i

    do i=0, option_count(trim(path)//"/scalar_field")-1
       child_path=trim(path)//"/scalar_field["//int2str(i)//"]"
       ! Reset path to have name instead of index
       call get_option(trim(child_path)//"/name", child_name)
       child_path=trim(path)//"/scalar_field::"//trim(child_name)
       call get_option(trim(complete_field_path(trim(child_path)))//"/mesh/name", &
                       mesh_name, default=trim(parent_mesh))
       call allocate_and_insert_scalar_field(child_path, state, &
            parent_mesh=mesh_name, parent_name=parent_name, &
            dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    end do
    do i=0, option_count(trim(path)//"/vector_field")-1
       child_path=trim(path)//"/vector_field["//int2str(i)//"]"
       ! Reset path to have name instead of index
       call get_option(trim(child_path)//"/name", child_name)
       child_path=trim(path)//"/vector_field::"//trim(child_name)
       call get_option(trim(complete_field_path(trim(child_path)))//"/mesh/name", &
                       mesh_name, default=trim(parent_mesh))
       call allocate_and_insert_vector_field(child_path, state, &
            parent_mesh=mesh_name, parent_name=parent_name, &
            dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    end do
    do i=0, option_count(trim(path)//"/tensor_field")-1
       child_path=trim(path)//"/tensor_field["//int2str(i)//"]"
       ! Reset path to have name instead of index
       call get_option(trim(child_path)//"/name", child_name)
       child_path=trim(path)//"/tensor_field::"//trim(child_name)
       call get_option(trim(complete_field_path(trim(child_path)))//"/mesh/name", &
                       mesh_name, default=trim(parent_mesh))
       call allocate_and_insert_tensor_field(child_path, state, &
            parent_mesh=mesh_name, parent_name=parent_name, &
            dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    end do

  end subroutine allocate_and_insert_children

  subroutine allocate_and_insert_grandchildren(path, state,&
       & parent_mesh, parent_name, dont_allocate_prognostic_value_spaces)
    !!< Allocate those fields contained in a field which are not direct
    !!< children. 
    character(len=*), intent(in) :: path !! option_path including prescribed/prognostic
    type(state_type), intent(inout) :: state
    character(len=*), intent(in) :: parent_mesh
    character(len=*), intent(in) :: parent_name
    logical, optional, intent(in) :: dont_allocate_prognostic_value_spaces

    
    integer :: i

    ! This is necessarily somewhat more sui generis than
    ! allocate_and_insert_children.

    do i = 1, size(grandchild_paths)
       call allocate_and_insert_children(trim(path)&
            &//trim(grandchild_paths(i)), state, parent_mesh, parent_name, &
            &dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    end do

  end subroutine allocate_and_insert_grandchildren
    
  logical function defer_allocation(option_path, mesh, &
     dont_allocate_prognostic_value_spaces)
    !!< Determines whether allocation of %val is deferred. 
    !!< This is used for fields that have been passed to SAM and
    !!< only are allocate one by one when we get them back from SAM.
    !!< Currently this is for all prognostic fields that are on a mesh
    !!< that is not excluded from mesh adaptivity.
    character(len=*), intent(in):: option_path
    type(mesh_type), intent(in):: mesh
    logical, optional, intent(in):: dont_allocate_prognostic_value_spaces
    
       defer_allocation=present_and_true(dont_allocate_prognostic_value_spaces) &
          .and. have_option(trim(option_path)//'/prognostic') &
          .and. .not. have_option(trim(mesh%option_path)//'/exclude_from_mesh_adaptivity')
  
  end function defer_allocation
  
  subroutine set_prescribed_field_values(states, &
    exclude_interpolated, exclude_nonreprescribed, initial_mesh)

    type(state_type), dimension(:), intent(in):: states
    !! don't prescribe the fields with interpolation options
    logical, intent(in), optional :: exclude_interpolated
    !! do not prescribe the fields that have requested not to be represcribed
    logical, intent(in), optional :: exclude_nonreprescribed
    !! indicates whether we're prescribing on the initial mesh, if not (default)
    !! the fields with needs_initial_mesh(field) are left untouched, they have to
    !! be interpolated (somewhere else)
    logical, intent(in), optional:: initial_mesh

    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield
    type(vector_field), pointer :: position
    character(len=OPTION_PATH_LEN):: phase_path
    logical :: mesh_changed
    integer :: p, f, nphases, nsfields, nvfields, ntfields

    ewrite(1,*) "In set_prescribed_field_values"

    mesh_changed = .not. present_and_true(initial_mesh)
    
    nphases = option_count('/material_phase')
    do p = 0, nphases-1

       phase_path = '/material_phase['//int2str(p)//']'

       ! Scalar fields:
       nsfields = scalar_field_count(states(p+1))
       do f = 1, nsfields
          sfield => extract_scalar_field(states(p+1),f)
          if (have_option(trim(sfield%option_path)//'/prescribed') .and. &
              .not. aliased(sfield) .and. &              
              .not. (present_and_true(exclude_interpolated) .and. &
                     interpolate_field(sfield)) .and. &
              .not. (present_and_true(exclude_nonreprescribed) .and. &
                     do_not_recalculate(sfield%option_path)) .and. &
              .not. (mesh_changed .and. needs_initial_mesh(sfield)) &
                     ) then
             
             position => get_external_coordinate_field(states(p+1), sfield%mesh)
             
             call zero(sfield)
             call initialise_field_over_regions(sfield, &
                trim(sfield%option_path)//'/prescribed/value', &
                position)
          end if
       end do

       nvfields = vector_field_count(states(p+1))
       do f = 1, nvfields
          vfield => extract_vector_field(states(p+1), f)
          if (have_option(trim(vfield%option_path)//'/prescribed') .and. &
              .not. aliased(vfield) .and. &              
              .not. (present_and_true(exclude_interpolated) .and. &
                     interpolate_field(vfield)) .and. &
              .not. (present_and_true(exclude_nonreprescribed) .and. &
                     do_not_recalculate(vfield%option_path)) .and. &
              .not. (mesh_changed .and. needs_initial_mesh(vfield)) &
                     ) then
                     
             position => get_external_coordinate_field(states(p+1), vfield%mesh)
             
             call zero(vfield)
             call initialise_field_over_regions(vfield, &
                trim(vfield%option_path)//'/prescribed/value', &
                position)
          end if
       end do
         
       ntfields = tensor_field_count(states(p+1))
       do f = 1, ntfields
          tfield => extract_tensor_field(states(p+1), f)
          if (have_option(trim(tfield%option_path)//'/prescribed') .and. &
              .not. aliased(tfield) .and. &              
              .not. (present_and_true(exclude_interpolated) .and. &
                     interpolate_field(tfield)) .and. &
              .not. (present_and_true(exclude_nonreprescribed) .and. &
                     do_not_recalculate(tfield%option_path)) .and. &
              .not. (mesh_changed .and. needs_initial_mesh(tfield)) &
                     ) then
                     
             position => get_external_coordinate_field(states(p+1), tfield%mesh)

             call zero(tfield)
             call initialise_field_over_regions(tfield, &
                trim(tfield%option_path)//'/prescribed/value', &
                position)
          end if
       end do
       
    end do

    if(have_option('/ocean_forcing/external_data_boundary_conditions')) then

      call set_nemo_fields(states(1))

    endif
      
    ! flush the cache
    call vtk_cache_finalise()

  end subroutine set_prescribed_field_values
  
  subroutine initialise_prognostic_fields(states, save_vtk_cache, &
    initial_mesh)
    !!< Set the values of prognostic fields with their initial conditions
    type(state_type), dimension(:), intent(in):: states
    !! By default the vtk_cache, build up by the from_file initialisations
    !! in this subroutine, is flushed at the end of this subroutine.  This
    !! cache can be reused however in subsequent calls reading from vtu files
    logical, intent(in), optional:: save_vtk_cache
    !! indicates whether we're initalising on the initial mesh, if not (default)
    !! the fields with needs_initial_mesh(field) are left untouched, they have to
    !! be interpolated (somewhere else)
    logical, intent(in), optional:: initial_mesh

    ! these must be pointers as bc's should be added to the original field
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield

    type(vector_field), pointer:: position
    character(len=OPTION_PATH_LEN):: phase_path
    integer p, f, nphases, nsfields, nvfields
    logical:: mesh_changed

    ewrite(1,*) "In initialise_prognostic_fields"
    
    mesh_changed = .not. present_and_true(initial_mesh)

    nphases = option_count('/material_phase')
    do p = 0, nphases-1

       phase_path = '/material_phase['//int2str(p)//']'

       position => extract_vector_field(states(p+1), "Coordinate")
       
       ! Scalar fields:
       nsfields = scalar_field_count(states(p+1))
       do f = 1, nsfields
          sfield => extract_scalar_field(states(p+1),f)
          if (mesh_changed .and. needs_initial_mesh(sfield)) cycle
          if (.not. aliased(sfield) .and. &
                have_option(trim(sfield%option_path)//'/prognostic')) then
             call zero(sfield)
             call initialise_field_over_regions(sfield, &
                trim(sfield%option_path)//'/prognostic/initial_condition', &
                position)
          end if
       end do

       nvfields = vector_field_count(states(p+1))
       do f = 1, nvfields
          vfield => extract_vector_field(states(p+1), f)
          if (mesh_changed .and. needs_initial_mesh(vfield)) cycle
          if (.not. aliased(vfield) .and. &
                have_option(trim(vfield%option_path)//'/prognostic')) then
             call zero(vfield)
             call initialise_field_over_regions(vfield, &
                trim(vfield%option_path)//'/prognostic/initial_condition', &
                position)
          end if
       end do

       if (have_option(trim(phase_path)//"/sediment")) then
          call initialise_prognostic_sediment
       end if

    end do
      
    if (.not. present_and_true(save_vtk_cache)) then
       ! flush the cache
       call vtk_cache_finalise()
    end if

  contains
    
    subroutine initialise_prognostic_sediment
      character(len=OPTION_PATH_LEN):: class_path, class_name, field_name

      nsfields=option_count(trim(phase_path)//"/sediment/sediment_class")
      do f=0,nsfields-1
         class_path=trim(phase_path)//"/sediment/sediment_class["&
              //int2str(f)//"]"

         ! Don't bother unless an additional initial condition is provided
         ! for this field (otherwise it defaults to the general case).
         if (.not.have_option(trim(class_path)//"/initial_condition")) then 
            cycle
         end if
         
         call get_option(trim(class_path)//"/name", class_name)
         field_name="SedimentConcentration"//trim(class_name)
         
         sfield => extract_scalar_field(states(p+1),field_name)

         if (mesh_changed .and. needs_initial_mesh(sfield)) cycle
         if (.not. aliased(sfield) .and. &
              have_option(trim(sfield%option_path)//'/prognostic')) then
            call zero(sfield)
            call initialise_field_over_regions(sfield, &
                 trim(class_path)//'/initial_condition', &
                 position)
         end if
         
      end do

    end subroutine initialise_prognostic_sediment

  end subroutine initialise_prognostic_fields

  subroutine allocate_and_insert_auxilliary_fields(states)
    ! Set up some auxilliary fields to the prognostic fields.
    ! i.e. old and iterated fields depending on options set
    type(state_type), dimension(:), intent(inout):: states

    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(vector_field), pointer :: positions

    type(scalar_field) :: aux_sfield
    type(vector_field) :: aux_vfield

    integer :: iterations
    logical :: steady_state_global, prognostic, prescribed, diagnostic, gravity

    character(len=FIELD_NAME_LEN) :: field_name
    character(len=OPTION_PATH_LEN) :: state_path, field_path

    integer :: nsfields, nvfields, p, f, p2, stat
    real :: current_time

    type(mesh_type), pointer :: u_mesh, x_mesh

    ewrite(1,*) "In allocate_and_insert_auxilliary_fields"

    call get_option("/timestepping/nonlinear_iterations", iterations, default=1)
    steady_state_global = have_option("/timestepping/steady_state")

    ! old and iterated fields
    do p = 1, size(states)

      ! Get number of scalar fields that are children of this state
      nsfields=scalar_field_count(states(p))

      ! Loop over scalar fields
      do f=1, nsfields

        sfield => extract_scalar_field(states(p), f)

        ! Save path to field
        field_path=trim(sfield%option_path)

        ! Get field name - this checks if the field has an option_path
        call get_option(trim(field_path)//"/name", field_name, stat)

        if((stat==0).and.(.not.aliased(sfield))) then

          prognostic=have_option(trim(sfield%option_path)//"/prognostic")
          prescribed=have_option(trim(sfield%option_path)//"/prescribed")
          diagnostic=have_option(trim(sfield%option_path)//"/diagnostic")

          ! if (prognostic or diagnostic) and (doing a steady state check on this field or doing more than 1 global iteration)
          if((prognostic.or.diagnostic)&
              .and.((steady_state_global.and.steady_state_field(sfield)).or.(iterations>1))) then

            call allocate(aux_sfield, sfield%mesh, "Old"//trim(sfield%name))
            call zero(aux_sfield)
            call insert(states(p), aux_sfield, trim(aux_sfield%name))
            call deallocate(aux_sfield)

          else

            aux_sfield = extract_scalar_field(states(p), trim(sfield%name))
            aux_sfield%name = "Old"//trim(sfield%name)
            aux_sfield%option_path=""  ! blank the option path so that it 
                                       ! doesn't get picked up in the next 
                                       ! aliased field loop
            aux_sfield%aliased=.true.
            call insert(states(p), aux_sfield, trim(aux_sfield%name))

          end if

          if((prognostic.or.diagnostic)&
              .and.(convergence_field(sfield).and.(iterations>1))) then

            call allocate(aux_sfield, sfield%mesh, "Iterated"//trim(sfield%name))
            call zero(aux_sfield)
            call insert(states(p), aux_sfield, trim(aux_sfield%name))
            call deallocate(aux_sfield)

          else

            aux_sfield = extract_scalar_field(states(p), trim(sfield%name))
            aux_sfield%name = "Iterated"//trim(sfield%name)
            aux_sfield%option_path=""  ! blank the option path so that it 
                                       ! doesn't get picked up in the next 
                                       ! aliased field loop
            aux_sfield%aliased=.true.
            call insert(states(p), aux_sfield, trim(aux_sfield%name))

          end if

        end if

      end do

      ! Get number of vector fields that are children of this state
      nvfields=vector_field_count(states(p))

      ! Loop over vector fields
      do f=1, nvfields

        vfield => extract_vector_field(states(p), f)

        ! Save path to field
        field_path=trim(vfield%option_path)

        ! Get field name - this checks if the field has an option_path
        call get_option(trim(field_path)//"/name", field_name, stat)

        if((stat==0).and.(.not.aliased(vfield))) then

          prognostic=have_option(trim(vfield%option_path)//"/prognostic")
          prescribed=have_option(trim(vfield%option_path)//"/prescribed")
          diagnostic=have_option(trim(vfield%option_path)//"/diagnostic")

          if((prognostic.or.diagnostic)&
             .and.((steady_state_global.and.steady_state_field(vfield)).or.(iterations>1))) then

            call allocate(aux_vfield, vfield%dim, vfield%mesh, "Old"//trim(vfield%name))
            call zero(aux_vfield)
            call insert(states(p), aux_vfield, trim(aux_vfield%name))
            call deallocate(aux_vfield)

          else

            aux_vfield = extract_vector_field(states(p), trim(vfield%name))
            aux_vfield%name = "Old"//trim(vfield%name)
            aux_vfield%option_path=""  ! blank the option path so that it 
                                       ! doesn't get picked up in the next 
                                       ! aliased field loop
            aux_vfield%aliased=.true.
            call insert(states(p), aux_vfield, trim(aux_vfield%name))

          end if

          if((prognostic.or.diagnostic)&
              .and.(convergence_field(vfield).and.(iterations>1))) then

            call allocate(aux_vfield, vfield%dim, vfield%mesh, "Iterated"//trim(vfield%name))
            call zero(aux_vfield)
            call insert(states(p), aux_vfield, trim(aux_vfield%name))
            call deallocate(aux_vfield)

          else

            aux_vfield = extract_vector_field(states(p), trim(vfield%name))
            aux_vfield%name = "Iterated"//trim(vfield%name)
            aux_vfield%option_path=""  ! blank the option path so that it 
                                       ! doesn't get picked up in the next 
                                       ! aliased field loop
            aux_vfield%aliased=.true.
            call insert(states(p), aux_vfield, trim(aux_vfield%name))

          end if

          if(trim(vfield%name)=="Velocity") then

            if(iterations>1) then

              call allocate(aux_vfield, vfield%dim, vfield%mesh, "Nonlinear"//trim(vfield%name))
              call zero(aux_vfield)
              call insert(states(p), aux_vfield, trim(aux_vfield%name))
              call deallocate(aux_vfield)

            else

              aux_vfield = extract_vector_field(states(p), trim(vfield%name))
              aux_vfield%name = "Nonlinear"//trim(vfield%name)
              aux_vfield%option_path=""
              aux_vfield%aliased = .true.
              call insert(states(p), aux_vfield, trim(aux_vfield%name))

            end if

            if(prognostic) then
              gravity = have_option("/physical_parameters/gravity")
              if(gravity) then
                sfield => extract_scalar_field(states(p), "Density", stat)
                if(stat==0) then
                  call allocate(aux_sfield, sfield%mesh, "VelocityBuoyancyDensity")
                else
                  call allocate(aux_sfield, vfield%mesh, "VelocityBuoyancyDensity")
                end if
                call zero(aux_sfield)
                aux_sfield%option_path=""
                call insert(states(p), aux_sfield, trim(aux_sfield%name))
                call deallocate(aux_sfield)
              end if
            end if
            
          end if

          if(trim(vfield%name)=="VelocityInnerElement") then

            if(iterations>1) then

              call allocate(aux_vfield, vfield%dim, vfield%mesh, "Nonlinear"//trim(vfield%name))
              call zero(aux_vfield)
              call insert(states(p), aux_vfield, trim(aux_vfield%name))
              call deallocate(aux_vfield)

            else

              aux_vfield = extract_vector_field(states(p), trim(vfield%name))
              aux_vfield%name = "Nonlinear"//trim(vfield%name)
              aux_vfield%option_path=""
              aux_vfield%aliased = .true.
              call insert(states(p), aux_vfield, trim(aux_vfield%name))

            end if

          end if

        end if

      end do

    end do

    ! old and iterated fields - aliased
    do p = 1, size(states)  ! now the aliased fields

      ! Get number of scalar fields that are children of this state
      nsfields=scalar_field_count(states(p))

      ! Loop over scalar fields
      do f=1, nsfields

        sfield => extract_scalar_field(states(p), f)

        ! Save path to field
        field_path=trim(sfield%option_path)

        ! Get field name - this checks if the field has an option_path
        ! but if it's aliased the name that it gets from the option path will be of the field it's aliased to!
        call get_option(trim(field_path)//"/name", field_name, stat)

        if((stat==0).and.aliased(sfield).and.(sfield%option_path(:15)=="/material_phase")) then

          prognostic=have_option(trim(sfield%option_path)//"/prognostic")
          prescribed=have_option(trim(sfield%option_path)//"/prescribed")
          diagnostic=have_option(trim(sfield%option_path)//"/diagnostic")

          if(prognostic&
             .or.(prescribed.and.(trim(field_name)==trim(sfield%name)))&
             .or.(diagnostic.and.(trim(field_name)==trim(sfield%name)))) then

            do p2 = 1, size(states)
              write(state_path, '(a,i0,a)') "/material_phase[",p2-1,"]"
              if(starts_with(trim(field_path), trim(state_path))) exit
            end do

            if(p2==size(states)+1) then
              FLAbort("scalar_field aliased but could not find to which material_phase")
            end if

            aux_sfield=extract_scalar_field(states(p2), "Old"//trim(field_name))
            aux_sfield%name = "Old"//trim(sfield%name)
            aux_sfield%aliased = .true.
            aux_sfield%option_path = ""  ! blank the option path for consistency
            call insert(states(p), aux_sfield, trim(aux_sfield%name))

            aux_sfield=extract_scalar_field(states(p2), "Iterated"//trim(field_name))
            aux_sfield%name = "Iterated"//trim(sfield%name)
            aux_sfield%aliased = .true.
            aux_sfield%option_path = ""  ! blank the option path for consistency
            call insert(states(p), aux_sfield, trim(aux_sfield%name))

          end if

        end if

      end do

      ! Get number of vector fields that are children of this state
      nvfields=vector_field_count(states(p))

      ! Loop over vector fields
      do f=1, nvfields

        vfield => extract_vector_field(states(p), f)

        ! Save path to field
        field_path=trim(vfield%option_path)

        ! Get field name - this checks if the field has an option_path
        ! but if it's aliased the name that it gets from the option path will be of the field it's aliased to!
        call get_option(trim(field_path)//"/name", field_name, stat)

        if((stat==0).and.aliased(vfield).and.(vfield%option_path(:15)=="/material_phase")) then

          prognostic=have_option(trim(vfield%option_path)//"/prognostic")
          prescribed=have_option(trim(vfield%option_path)//"/prescribed")
          diagnostic=have_option(trim(vfield%option_path)//"/diagnostic")

          do p2 = 1, size(states)
            write(state_path, '(a,i0,a)') "/material_phase[",p2-1,"]"
            if(starts_with(trim(field_path), trim(state_path))) exit
          end do

          if(p2==size(states)+1) then
            FLAbort("vector_field aliased but could not find to which material_phase")
          end if

          if(prognostic&
             .or.(prescribed.and.(trim(field_name)==trim(vfield%name)))&
             .or.(diagnostic.and.(trim(field_name)==trim(vfield%name)))) then

            aux_vfield=extract_vector_field(states(p2), "Old"//trim(field_name))
            aux_vfield%name = "Old"//trim(vfield%name)
            aux_vfield%aliased = .true.
            aux_vfield%option_path = ""  ! blank the option path for consistency
            call insert(states(p), aux_vfield, trim(aux_vfield%name))

            aux_vfield=extract_vector_field(states(p2), "Iterated"//trim(field_name))
            aux_vfield%name = "Iterated"//trim(vfield%name)
            aux_vfield%aliased = .true.
            aux_vfield%option_path = ""  ! blank the option path for consistency
            call insert(states(p), aux_vfield, trim(aux_vfield%name))

          end if

          if(trim(vfield%name)=="Velocity") then

            aux_vfield=extract_vector_field(states(p2), "Nonlinear"//trim(field_name))
            aux_vfield%name = "Nonlinear"//trim(vfield%name)
            aux_vfield%aliased = .true.
            aux_vfield%option_path = ""
            call insert(states(p), aux_vfield, trim(aux_vfield%name))

          end if

          if(trim(vfield%name)=="VelocityInnerElement") then

            aux_vfield=extract_vector_field(states(p2), "Nonlinear"//trim(field_name))
            aux_vfield%name = "Nonlinear"//trim(vfield%name)
            aux_vfield%aliased = .true.
            aux_vfield%option_path = ""
            call insert(states(p), aux_vfield, trim(aux_vfield%name))

          end if

        end if

      end do

    end do
      
    
    ! for mesh movement we need a "OriginalCoordinate" and "OldCoordinate" field
    ! inserted in each state similar to "Coordinate"
    if (have_option('/mesh_adaptivity/mesh_movement')) then 
       positions => extract_vector_field(states(1), name="Coordinate")
       
       ! first original coordinate field:
       call allocate(aux_vfield, positions%dim, positions%mesh, &
          name="OriginalCoordinate")
       ! at the moment original_positions is repointed to (XORIG, YORIG, ZORIG)
       ! after we've fixed that we need to copy it here:
       call set(aux_vfield, positions)
       aux_vfield%option_path=""
       ! insert in all states
       call insert(states, aux_vfield, name="OriginalCoordinate")
       call deallocate(aux_vfield)
       
       ! exactly the same for old coordinate field:
       call allocate(aux_vfield, positions%dim, positions%mesh, &
          name="OldCoordinate")
       ! at the moment OldPositions is repointed to (XOLD, YOLD, ZOLD)
       ! after we've fixed that we need to copy it here:
       call set(aux_vfield, positions)
       aux_vfield%option_path=""
       ! insert in all states
       call insert(states, aux_vfield, name="OldCoordinate")
       call deallocate(aux_vfield)
        
    end if

    u_mesh => extract_mesh(states(1), "VelocityMesh")
    x_mesh => extract_mesh(states(1), "CoordinateMesh")
    ! note that GridVelocity is on the velocity mesh
    ! rather than the coordinate mesh you might expect
    ! this is because its primary use is to subtract it from Velocity
    call allocate(aux_vfield, mesh_dim(u_mesh), u_mesh, "GridVelocity")
    call zero(aux_vfield)
    aux_vfield%option_path = ""
    call insert(states, aux_vfield, trim(aux_vfield%name))
    call deallocate(aux_vfield)

    ! Disgusting and vomitous hack to ensure that time is output in
    ! vtu files.
    call allocate(aux_sfield, x_mesh, "Time", field_type=FIELD_TYPE_CONSTANT)
    call get_option("/timestepping/current_time", current_time)
    call set(aux_sfield, current_time)
    aux_sfield%option_path = ""
    call insert(states, aux_sfield, trim(aux_sfield%name))
    call deallocate(aux_sfield)

  end subroutine allocate_and_insert_auxilliary_fields

  function mesh_name(field_path)
    !!< given a field path, establish the mesh that the field is on.
    use global_parameters, only: FIELD_NAME_LEN
    character(len=FIELD_NAME_LEN) :: mesh_name
    character(len=*), intent(in) :: field_path

    integer :: stat

    call get_option(trim(field_path)//'/prognostic/mesh[0]/name', &
         mesh_name, stat=stat)
    if (stat/=0) then
       call get_option(trim(field_path)//'/diagnostic/mesh[0]/name', &
            mesh_name, stat=stat)
       if (stat/=0) then
          call get_option(trim(field_path)//'/prescribed/mesh[0]/name', &
               mesh_name, stat=stat)
          if (stat/=0) then
             FLAbort("No mesh for field "//trim(field_path))
          end if
       end if
    end if

  end function mesh_name
    
  subroutine surface_id_stats(mesh, positions)
  type(mesh_type), target, intent(in):: mesh
  type(vector_field), intent(in):: positions
    
    real, dimension(1:face_ngi(mesh,1)):: detwei
    integer, dimension(:), pointer:: surface_ids
    real, dimension(:), allocatable:: area
    integer, dimension(:), allocatable:: no_elements
    integer i, sid, sidmin, sidmax
    
    if (current_debug_level<=1) return
    
    ewrite(2,*) "Surface id stats for mesh ", trim(mesh%name)
    
    surface_ids => mesh%faces%boundary_ids
    sidmin=minval(surface_ids)
    sidmax=maxval(surface_ids)
    
    allocate( no_elements(sidmin:sidmax), &
      area(sidmin:sidmax))
    no_elements=0
    area=0.0
    
    do i=1, surface_element_count(mesh)
       sid=surface_ids(i)
       no_elements(sid)=no_elements(sid)+1
       call transform_facet_to_physical(positions, i, detwei_f=detwei)
       area(sid)=area(sid)+sum(detwei)
    end do
      
    ewrite(2, *) 'Surface id,  n/o surface elements,       surface area'
    do i=sidmin, sidmax
       ewrite(2, "(i10,i23,es20.9)") i, no_elements(i), area(i)
    end do
      
    ewrite(2,*) 'Total number of surface elements:', surface_element_count(mesh)
    ewrite(2,'(a,es20.9)') 'Total surface area:', sum(area)
    
  end subroutine surface_id_stats
    
  subroutine populate_state_module_check_options

    character(len=OPTION_PATH_LEN) :: problem_type

    ! Check mesh options
    call check_mesh_options

    ! check problem specific options:
    call get_option("/problem_type", problem_type)
    select case (problem_type)
    case ("fluids")
    case ("oceans")
       call check_ocean_options
    case ("solids")
       call check_solid_options
    case ("multimaterial")
       call check_multimaterial_options
    case ("porous_media")
       call check_porous_media_options
    case ("stokes")
       call check_stokes_options
    case default
       ewrite(0,*) "Problem type:", trim(problem_type)
       FLExit("Error unknown problem_type")
    end select
    ewrite(2,*) 'Done with problem type choice'

  end subroutine populate_state_module_check_options

  character(len=OPTION_PATH_LEN) function AddFieldTypeOption(path)
    implicit none
    ! Auxillary function to add  prognostic/diagnostic/prescribed/aliased to field option path
    ! to be used with checking options as we die friendly with FLExit
    ! character(len=OPTION_PATH_LEN) :: AddFieldTypeOption
    character(len=*), intent(in) :: path

    if (have_option(trim(path) // "/prognostic")) then

       AddFieldTypeOption=trim(path) // "/prognostic"

    else if (have_option(trim(path) // "/diagnostic")) then

       AddFieldTypeOption=trim(path) // "/diagnostic"

    else if (have_option(trim(path) // "/prescribed")) then

       AddFieldTypeOption=trim(path) // "/prescribed"

    else if (have_option(trim(path) // "/aliased")) then

       AddFieldTypeOption=trim(path)// "/aliased"

    else

       ewrite(0,*) "For field:", trim(path)
       FLExit("Error: unknown field type")

    end if

  end function AddFieldTypeOption

  subroutine check_mesh_options

    character(len=OPTION_PATH_LEN) :: path
    character(len=OPTION_PATH_LEN) :: field_name, mesh_name, phase_name
    integer :: i, j ! counters
    integer :: nstates ! number of states
    integer :: nfields ! number of fields
    integer :: nmeshes ! number of meshes
    integer :: n_external_meshes ! number of meshes from file
    ! logicals to find out if we have certain options
    logical :: is_aliased

    ! Get number of meshes
    nmeshes=option_count("/geometry/mesh")

    ewrite(2,*) "Checking mesh options."
    ewrite(2,*) "There are", nmeshes, "meshes."

    n_external_meshes=0

    mesh_loop1: do i=0, nmeshes-1

       ! Save mesh path
       path="/geometry/mesh["//int2str(i)//"]"

       if(have_option(trim(path)//"/from_file")) then

          n_external_meshes=n_external_meshes+1

       else if (.not. have_option(trim(path)//"/from_mesh")) then

          call get_option(trim(path)//"/name", mesh_name)
          ewrite(0,*) "In options for /geometry/mesh ("//trim(mesh_name)//"):"
          FLExit("Error: unknown way of specifying mesh source.")

       end if

    end do mesh_loop1

    ! Check that at least one mesh is read in from a file.
    if(n_external_meshes==0) then
       FLExit("At least one mesh must come from a file.")
    end if
    if(isparallel() .and. n_external_meshes > 1) then
      FLExit("Only one mesh may be from_file in parallel.")
    end if

    ! Check that dimension of mesh is the same as the dimension defined in the options file
    ! ...that's not so easy: let's just do this in insert_external_mesh() with a nice FLEXit

    ! Check that the meshes required to make other meshes are present.
    mesh_loop2: do i=0, nmeshes-1

       ! Save mesh path
       path="/geometry/mesh["//int2str(i)//"]"

       if (have_option(trim(path)//"/from_mesh")) then

          call get_option(trim(path)//"/from_mesh/mesh[0]/name", mesh_name)
          if (.not. have_option("/geometry/mesh::"//trim(mesh_name))) then

             ewrite(0,*) "Unknown mesh: ", trim(mesh_name)
             call get_option(trim(path)//"/name", mesh_name)
             ewrite(0,*) "Specified as source (from_mesh) for ", trim(mesh_name)
             FLExit("Error in /geometry/mesh: unknown mesh.")

          end if
          
          if (have_option(trim(path)//"/from_mesh/extrude") .and. ( &
             
             have_option(trim(path)//"/from_mesh/mesh_shape") .or. &
             have_option(trim(path)//"/from_mesh/mesh_continuity") .or. &
             have_option(trim(path)//"/from_mesh/periodic_boundary_conditions") &
             ) ) then
             
             ewrite(0,*) "When extruding a mesh, you cannot at the same time"
             ewrite(0,*) "change its shape, continuity or add periodic bcs."
             ewrite(0,*) "Need to do this in seperate step (derivation)."
             FLExit("Error in /geometry/mesh with extrude option")
             
          end if

       end if

    end do mesh_loop2

    ! Check that mesh associated with each field exists

    nstates=option_count("/material_phase")

    state_loop: do i=0, nstates-1

       call get_option("/material_phase["//int2str(i)//"]/name", phase_name)

       ! Get number of scalar fields that are children of this state
       nfields=option_count("/material_phase["//int2str(i)//"]/scalar_field")

       ! Loop over scalar fields
       scalar_field_loop: do j=0, nfields-1

          ! Save path to field
          path="/material_phase["//int2str(i)//"]/scalar_field["&
               &//int2str(j)//"]"
          ! Get field name
          call get_option(trim(path)//"/name", field_name)
          ! Reset path to have field name rather than index
          path="/material_phase["//int2str(i)//"]/scalar_field::"//trim(field_name)

          ! If field is not aliased check mesh name
          is_aliased=have_option(trim(path)//"/aliased")
          if(.not.is_aliased) then
             call get_option(trim(AddFieldTypeOption(path))//"/mesh[0]/name", mesh_name)

             if (.not. have_option("/geometry/mesh::"//trim(mesh_name))) then

                ewrite(0,*) "Unknown mesh: ", trim(mesh_name)
                ewrite(0,*) "Specified as mesh for scalar_field ", trim(field_name)
                ewrite(0,*) "In material_phase ", trim(phase_name)
                FLExit("Error: unknown mesh.")

             end if

          end if

       end do scalar_field_loop

       ! Get number of vector fields that are children of this state
       nfields=option_count("/material_phase["//int2str(i)//"]/vecto&
            &r_field")

       ! Loop over vector fields
       vector_field_loop: do j=0, nfields-1

          ! Save path to field
          path="/material_phase["//int2str(i)//"]/vector_field["&
               &//int2str(j)//"]"
          ! Get field name
          call get_option(trim(path)//"/name", field_name)
          ! Reset path to have field name rather than index
          path="/material_phase["//int2str(i)//"]/vector_field::"//trim(field_name)

          ! If field is not aliased check mesh name
          is_aliased=have_option(trim(path)//"/aliased")
          if(.not.is_aliased) then
             call get_option(trim(AddFieldTypeOption(path))//"/mesh[0]/name", mesh_name)

             if (.not. have_option("/geometry/mesh::"//trim(mesh_name))) then

                ewrite(0,*) "Unknown mesh: ", trim(mesh_name)
                ewrite(0,*) "Specified as mesh for vector_field ", trim(field_name)
                ewrite(0,*) "In material_phase ", trim(phase_name)
                FLExit("Error: unknown mesh.")

             end if

          end if

       end do vector_field_loop

       ! Get number of tensor fields that are children of this state
       nfields=option_count("/material_phase["//int2str(i)//"]/tensor_field")

       tensor_field_loop: do j=0, nfields-1

          ! Save path to field
          path="/material_phase["//int2str(i)//"]/tensor_field["&
               &//int2str(j)//"]"
          ! Get field name
          call get_option(trim(path)//"/name", field_name)
          ! Reset path to have field name rather than index
          path="/material_phase["//int2str(i)//"]/tensor_field::"//trim(field_name)

          ! If field is not aliased check mesh name
          is_aliased=have_option(trim(path)//"/aliased")
          if(.not.is_aliased) then
             call get_option(trim(AddFieldTypeOption(path))//"/mesh[0]/name", mesh_name)

             if (.not. have_option("/geometry/mesh::"//trim(mesh_name))) then

                ewrite(0,*) "Unknown mesh: ", trim(mesh_name)
                ewrite(0,*) "Specified as mesh for tensor_field ", trim(field_name)
                ewrite(0,*) "In material_phase ", trim(phase_name)
                FLExit("Error: unknown mesh.")

             end if

          end if

       end do tensor_field_loop

    end do state_loop

  end subroutine check_mesh_options
  
  subroutine check_ocean_options

    character(len=OPTION_PATH_LEN) str, velocity_path, pressure_path, tmpstring
    character(len=OPTION_PATH_LEN) free_surface_path
    logical on_sphere, constant_gravity, new_navsto
    
    if (option_count('/material_phase')/=1) then
       FLExit("The checks for problem_type oceans only work for single phase.")
    endif
       
    ! from now on we may assume single material/phase    
    velocity_path="/material_phase[0]/vector_field::Velocity/prognostic"
    if (have_option(trim(velocity_path))) then
       new_navsto=have_option(trim(velocity_path)//'/spatial_discretisation/continuous_galerkin') .or. &
          have_option(trim(velocity_path)//'/spatial_discretisation/discontinuous_galerkin')
       ! Check that for ocean problems with prognostic velocity the mass is lumped
       ! in case of Continuous Galerkin:
       str=trim(velocity_path)//'/spatial_discretisation/legacy_continuous_galerkin'
       if (have_option(trim(str)) .and. .not. &
            have_option(trim(str)//"/lump_mass_matrix")) then
          ewrite(0,*) "Missing option spatial_discretisation/legacy_continuous_galerkin/lump_mass_matrix"
          ewrite(0,*) "under the prognostic velocity field."
          FLExit("For ocean problems you need to lump the mass matrix.")
       end if
       ! in case of legacy discretisation options:
       str=trim(velocity_path)//'/spatial_discretisation/legacy_discretisation'
       if (have_option(trim(str)).and. .not. &
            have_option(trim(str)//"/legacy_mlump")) then
          ewrite(0,*) "Missing option spatial_discretisation/legacy_discretisation/legacy_mlump"
          ewrite(0,*) "under the prognostic velocity field."
          FLExit("For ocean problems you need to lump the mass matrix.")
       end if
       
       ! check we have the right equation type for velocity
       if (.not. have_option(trim(velocity_path)//'/equation::Boussinesq')) then
          ewrite(0,*) "For ocean problems you need to set the equation type"
          ewrite(0,*) "for velocity to Boussinesq."
          FLExit("Wrong Velocity equation type")
       end if          

    end if

    pressure_path="/material_phase[0]/scalar_field::Pressure/prognostic"
    if(have_option("/material_phase[0]/scalar_field::Pressure/prognostic")) then
       if (.not.have_option(trim(pressure_path)//"/scheme/use_projection_method")) then
          FLExit("For ocean problems you should use the projection method under scheme for pressure")
       end if
       call get_option(trim(pressure_path)//"/scheme/poisson_pressure_solution", tmpstring)
       select case (tmpstring)
       case ("never", "every timestep")
          ewrite(0,*) ("WARNING: For ocean problems you should use the Poisson pressure solution at the first timestep only.")
       end select
    end if

    ! Warning about salinity options
    ! Density is only affected when you have salinity and either linear EoS with
    ! salinity dependency or Pade Approximation turned on
    if(have_option("/material_phase[0]/scalar_field::Salinity")&
         .and.(.not.(have_option("/material_phase[0]/equation_of_state/fluids/linear/salinity_dependency") .or.&
                     have_option("/material_phase[0]/equation_of_state/fluids/ocean_pade_approximation")))) then
       ewrite(0,*) "WARNING: You have a salinity field but it will not affect the density of the fluid."
    end if

    ! Check that vertical mesh movement is enabled if a free surface is present
    free_surface_path="/material_phase[0]/scalar_field::FreeSurface/prognostic"
    if(have_option(free_surface_path)) then
       if (.not. have_option("/mesh_adaptivity/mesh_movement")) then
          ewrite(0,*) "With a free surface you need to allow the mesh to move vertically"
          ewrite(0,*) "by adding the option /mesh_adaptivity/mesh_movement."
          FLExit("Missing /mesh_adaptivity/mesh_movement element")
       end if
    end if

    ! Check that the gravity field is not constant for spherical problems
    on_sphere=have_option("/geometry/spherical_earth")
    constant_gravity=have_option("/physical_parameters/gravity/vector_field::GravityDirection/prescribed/value[0]/constant")
    if(on_sphere .and. constant_gravity) then
       ewrite(0,*) "If you are using spherical geometry you cannot have"
       ewrite(0,*) "a constant gravity direction."
       ewrite(0,*) "See the waterworld test case for an example of how"
       ewrite(0,*) "to set this properly"
       FLExit("GravityDirection set incorrectly for spherical geometry.")
    end if

  end subroutine check_ocean_options

  subroutine check_multimaterial_options

    integer :: neos, nmat, i
    logical :: have_vfrac, have_dens
    
    integer :: cohesion_count, diagnosticvolumefraction_count, density_count, &
               frictionangle_count, viscosity_count, surfacetension_count, &
               elasticity_count

    neos = option_count("/material_phase/equation_of_state/multimaterial")
    nmat = option_count("/material_phase")

    if(neos>0) then
       if(nmat/=neos) then
          FLExit("Not all the material_phases have compressible equations of state.")
       end if
    end if

    do i = 0, nmat-1
       have_vfrac = have_option("/material_phase["//int2str(i)//&
            "]/scalar_field::MaterialVolumeFraction")
       have_dens = have_option("/material_phase["//int2str(i)//&
            "]/scalar_field::MaterialDensity").or.&
            have_option("/material_phase["//int2str(i)//&
            "]/equation_of_state/fluids/linear/reference_density")
       if((.not.have_vfrac).or.(.not.have_dens)) then
          FLExit("All material_phases need a MaterialVolumeFraction and either a MaterialDensity or an eos.")
       end if
    end do
    
    diagnosticvolumefraction_count = option_count('/material_phase/&
                &scalar_field::MaterialVolumeFraction/&
                &diagnostic')
    if(diagnosticvolumefraction_count>1) then
      ewrite(-1,*) diagnosticvolumefraction_count, 'diagnostic MaterialVolumeFractions.'
      FLExit("Only 1 diagnostic MaterialVolumeFraction is allowed")
    end if

    density_count = option_count('/material_phase/&
                &scalar_field::Density/diagnostic')
    if(density_count>1) then
      ewrite(-1,*) density_count, 'diagnostic bulk Densities.'
      FLExit("Only 1 diagnostic bulk Density is allowed")
    end if

    frictionangle_count = option_count('/material_phase/&
                &vector_field::Velocity/prognostic&
                &/scalar_field::FrictionAngle/diagnostic')
    if(frictionangle_count>1) then
      ewrite(-1,*) frictionangle_count, 'diagnostic bulk FrictionAngles.'
      FLExit("Only 1 diagnostic bulk FrictionAngle is allowed")
    end if

    cohesion_count = option_count('/material_phase/&
                &vector_field::Velocity/prognostic&
                &/scalar_field::Cohesion/diagnostic')
    if(cohesion_count>1) then
      ewrite(-1,*) cohesion_count, 'diagnostic bulk Cohesions.'
      FLExit("Only 1 diagnostic bulk Cohesion is allowed")
    end if

    viscosity_count = option_count('/material_phase/&
                &vector_field::Velocity/prognostic/&
                &tensor_field::Viscosity/diagnostic')
    if(viscosity_count>1) then
      ewrite(-1,*) viscosity_count, 'diagnostic bulk Viscosities.'
      FLExit("Only 1 diagnostic bulk Viscosity is allowed")
    end if

    elasticity_count = option_count('/material_phase/&
                &vector_field::Velocity/prognostic&
                &/tensor_field::Elasticity/diagnostic')
    if(elasticity_count>1) then
      ewrite(-1,*) elasticity_count, 'diagnostic bulk Elasticities.'
      FLExit("Only 1 diagnostic bulk Elasticity is allowed")
    end if
    
    surfacetension_count = option_count('/material_phase/&
                &vector_field::Velocity/prognostic&
                &/tensor_field::SurfaceTension/diagnostic')
    if(surfacetension_count>1) then
      ewrite(-1,*) surfacetension_count, 'diagnostic surface tensions.'
      FLExit("Only 1 diagnostic surface tension is allows")
    end if

  end subroutine check_multimaterial_options

  subroutine check_solid_options()

    integer :: nmat, i, counter

    nmat = option_count("/material_phase")

    do i = 0, nmat-1
       if(have_option("/material_phase["//int2str(i)//&
            "]/vector_field::Velocity/prognostic")) then
          if(.not.have_option("/material_phase["//int2str(i)//&
            "]/vector_field::Displacement")) then
            FLExit("Displacement field required for solid problems.")
          end if
          if(.not.have_option("/material_phase["//int2str(i)//&
            "]/vector_field::Velocity/prognostic/tensor_field::Elasticity")) then
            FLExit("Elasticity field required for solid problems.")
          end if
          if(have_option("/material_phase["//int2str(i)//&
            "]/vector_field::Velocity/prognostic/tensor_field::Elasticity")) then
            counter = option_count("/material_phase["//int2str(i)//&
                       "/vector_field::Velocity/prognostic/tensor_field::Elasticity/prescribed/value/isotropic")
            if(counter/=option_count("/material_phase["//int2str(i)//&
                       "/vector_field::Velocity/prognostic/tensor_field::Elasticity/prescribed/value")) then
              FLExit("Elasticity must be isotropic for solid problems.")
            end if
          end if
          if(have_option("/material_phase["//int2str(i)//&
            "]/vector_field::Velocity/prognostic/tensor_field::Viscosity")) then
            counter = option_count("/material_phase["//int2str(i)//&
                       "/vector_field::Velocity/prognostic/tensor_field::Viscosity/prescribed/value/isotropic")
            if(counter/=option_count("/material_phase["//int2str(i)//&
                       "/vector_field::Velocity/prognostic/tensor_field::Viscosity/prescribed/value")) then
              FLExit("Viscosity must be isotropic for solid problems.")
            end if
          end if
       end if
    end do

  end subroutine check_solid_options

  subroutine check_porous_media_options

    integer :: nmat, i
    logical :: have_vfrac, have_viscosity, have_porosity, have_permeability

    nmat = option_count("/material_phase")
    ewrite(2,*) 'nmat:',nmat

    have_porosity = have_option("/porous_media/scalar_field::Porosity")
    have_permeability = have_option("/porous_media/scalar_field::Permeability").or.&
                       &have_option("/porous_media/vector_field::Permeability").or.&
                       &have_option("/porous_media/tensor_field::Permeability")
    if((.not.have_porosity).or.(.not.have_permeability)) then
       FLExit("For porous media problems we need porosity and permeability.")
    end if
! Need to sort this out for multiphase!!!
    do i = 0, nmat-1
       if(have_option("/porous_media/multiphase_parameters")) then
          have_vfrac = have_option("/material_phase["//int2str(i)//&
               "]/scalar_field::PhaseVolumeFraction")
          have_viscosity = have_option("/materical_phase["//int2str(i)//&
               "]/tensor_field::MaterialViscosity")
          if((.not.have_vfrac).or.(.not.have_viscosity)) then
             FLExit("Need volume fractions and viscosities for each material phase.")
          endif
       endif
    end do

  end subroutine check_porous_media_options

  subroutine check_stokes_options

    ! Check options for Stokes flow simulations.

    character(len=OPTION_PATH_LEN) velocity_path, pressure_path

    velocity_path="/material_phase[0]/vector_field::Velocity/prognostic"
    if (have_option(trim(velocity_path))) then

       ! Check to see that mass terms are switched off:
       if(.not.have_option(trim(velocity_path)//'/spatial_discretisation/continuous_galerkin/mass_terms/exclude_mass_terms')) then
          ewrite(0,*) "Mass terms are invalid for a Stokes flow problem"
          FLExit("For Stokes problems you need to exclude the mass matrix.")
       end if

       ! Check to see that advection terms are switched off:
       if(.not.have_option(trim(velocity_path)//'/spatial_discretisation/continuous_galerkin/advection_terms/exclude_advection_terms')) then
          ewrite(0,*) "Advection terms are invalid for a Stokes flow problem"
          FLExit("For Stokes problems you need to exclude the advection terms.")    
       end if

    end if

    pressure_path="/material_phase[0]/scalar_field::Pressure/prognostic"

  end subroutine check_stokes_options

  subroutine allocate_metric_limits(state)
    type(state_type), intent(inout) :: state
    type(tensor_field) :: min_edge, max_edge
    type(tensor_field) :: min_eigen, max_eigen
    character(len=*), parameter :: path = &
    & "/mesh_adaptivity/hr_adaptivity/"
    logical :: is_constant
    type(mesh_type), pointer :: mesh
    type(vector_field), pointer :: X
    integer :: node

    if (.not. have_option(path)) then
      return
    end if

    X => extract_vector_field(state, "Coordinate")
    ! We can't use the external mesh in the extruded case -- these have to go on the
    ! CoordinateMesh.
    !mesh => get_external_mesh((/state/))
    mesh => X%mesh

    if (.not. have_option(path // "/tensor_field::MinimumEdgeLengths")) then
      ewrite(-1,*) "Warning: adaptivity turned on, but no edge length limits available?"
      return
    end if

    is_constant = (have_option(path // "/tensor_field::MinimumEdgeLengths/anisotropic_symmetric/constant"))
    if (is_constant) then
      call allocate(min_edge, mesh, "MinimumEdgeLengths", field_type=FIELD_TYPE_CONSTANT)
      call initialise_field(min_edge, path // "/tensor_field::MinimumEdgeLengths", X)
      call allocate(max_eigen, mesh, "MaxMetricEigenbound", field_type=FIELD_TYPE_CONSTANT)
      call set(max_eigen, eigenvalue_from_edge_length(node_val(min_edge, 1)))
    else
      call allocate(min_edge, mesh, "MinimumEdgeLengths")
      call initialise_field(min_edge, path // "/tensor_field::MinimumEdgeLengths", X)
      call allocate(max_eigen, mesh, "MaxMetricEigenbound")
      do node=1,node_count(mesh)
        call set(max_eigen, node, eigenvalue_from_edge_length(node_val(min_edge, node)))
      end do
    end if
      
    call insert(state, max_eigen, "MaxMetricEigenbound")
    call deallocate(min_edge)
    call deallocate(max_eigen)

    is_constant = (have_option(path // "/tensor_field::MaximumEdgeLengths/anisotropic_symmetric/constant"))
    if (is_constant) then
      call allocate(max_edge, mesh, "MaximumEdgeLengths", field_type=FIELD_TYPE_CONSTANT)
      call initialise_field(max_edge, path // "/tensor_field::MaximumEdgeLengths", X)
      call allocate(min_eigen, mesh, "MinMetricEigenbound", field_type=FIELD_TYPE_CONSTANT)
      call set(min_eigen, eigenvalue_from_edge_length(node_val(max_edge, 1)))
    else
      call allocate(max_edge, mesh, "MaximumEdgeLengths")
      call initialise_field(max_edge, path // "/tensor_field::MaximumEdgeLengths", X)
      call allocate(min_eigen, mesh, "MinMetricEigenbound")
      do node=1,node_count(mesh)
        call set(min_eigen, node, eigenvalue_from_edge_length(node_val(max_edge, node)))
      end do
    end if
      
    call insert(state, min_eigen, "MinMetricEigenbound")
    call deallocate(max_edge)
    call deallocate(min_eigen)

  end subroutine allocate_metric_limits

  function get_quad_family() result(quad_family)
    character(len=OPTION_PATH_LEN) :: quad_family_str
    integer :: quad_family

    if (have_option("/geometry/quadrature/quadrature_family")) then
      call get_option("/geometry/quadrature/quadrature_family", quad_family_str)
      select case (quad_family_str)
        case("family_cools")
          quad_family = FAMILY_COOLS
        case("family_grundmann_moeller")
          quad_family = FAMILY_GM
        case ("family_wandzura")
          quad_family = FAMILY_WANDZURA
      end select
    else
      quad_family = FAMILY_COOLS
    end if
  end function get_quad_family

end module populate_state_module
