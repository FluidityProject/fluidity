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
module populate_state_module
  use elements
  use state_module
  use FLDebug
  use spud
  use mesh_files
  use vtk_cache_module
  use global_parameters, only: OPTION_PATH_LEN, is_active_process, pi, &
    no_active_processes, topology_mesh_name, adaptivity_mesh_name, &
    periodic_boundary_option_path, domain_bbox, domain_volume, surface_radius
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
  use hadapt_extrude_radially
  use initialise_fields_module
  use transform_elements
  use parallel_tools
  use boundary_conditions_from_options
  use nemo_states_module
  use data_structures
  use fields_halos
  use read_triangle
  use initialise_ocean_forcing_module

  implicit none

  private

  public populate_state
  public populate_state_module_check_options
  public insert_external_mesh, insert_derived_meshes, &
       allocate_field_as_constant, allocate_and_insert_fields, &
       initialise_prognostic_fields, set_prescribed_field_values, &
       alias_fields, mesh_name, &
       allocate_and_insert_auxilliary_fields, &
       initialise_field, allocate_metric_limits, &
       make_mesh_periodic_from_options, make_mesh_unperiodic_from_options, &
       compute_domain_statistics

  interface allocate_field_as_constant
    
    module procedure allocate_field_as_constant_scalar, allocate_field_as_constant_vector, &
          allocate_field_as_constant_tensor

  end interface allocate_field_as_constant
    
  !! A list of locations in which additional scalar/vector/tensor fields
  !! are to be found. These are absolute paths in the schema.
  character(len=OPTION_PATH_LEN), dimension(9) :: additional_fields_absolute=&
       (/ &
       "/ocean_biology/pznd                                                                                                   ", &
       "/ocean_biology/six_component                                                                                          ", &
       "/ocean_forcing/iceshelf_meltrate/Holland08                                                                            ", &
       "/ocean_forcing/bulk_formulae/output_fluxes_diagnostics                                                                ", &
       "/porous_media                                                                                                         ", &
       "/porous_media_dual                                                                                                    ", &
       "/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/les_model/dynamic_les ", &
       "/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/les_model/second_order", &
       "/material_phase[0]/sediment/                                                                                          " &
       /)
       
  !! A list of relative paths under /material_phase[i]
  !! that are searched for additional fields to be added.
  character(len=OPTION_PATH_LEN), dimension(13) :: additional_fields_relative=&
       (/ &
       "/subgridscale_parameterisations/Mellor_Yamada                                                       ", &
       "/subgridscale_parameterisations/prescribed_diffusivity                                              ", &
       "/subgridscale_parameterisations/GLS                                                                 ", &
       "/subgridscale_parameterisations/k-epsilon                                                           ", &
       "/subgridscale_parameterisations/k-epsilon/debugging_options/source_term_output_fields               ", &
       "/subgridscale_parameterisations/k-epsilon/debugging_options/prescribed_source_terms                 ", &
       "/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/les_model/second_order", &
       "/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/les_model/fourth_order", &
       "/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/les_model/wale        ", &
       "/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/les_model/dynamic_les ", &
       "/vector_field::Velocity/prognostic/equation::ShallowWater                                           ", &
       "/vector_field::Velocity/prognostic/equation::ShallowWater/bottom_drag                               ", &
       "/vector_field::BedShearStress/diagnostic/calculation_method/velocity_gradient                       " &
       /)

  !! Relative paths under a field that are searched for grandchildren
  !! (moved here because of extremely obscure intel ICE -Stephan)
  character(len=OPTION_PATH_LEN), dimension(1):: &
         grandchild_paths = (/&
         &    "/spatial_discretisation/inner_element" &
         /)

contains


  subroutine populate_state(states)
    use Profiler
    type(state_type), pointer, dimension(:) :: states

    integer :: nstates ! number of states
    integer :: i

    ewrite(1,*) "In populate_state"
    call profiler_tic("I/O")
    call tictoc_clear(TICTOC_ID_IO_READ)

    ! Find out how many states there are
    nstates=option_count("/material_phase")
    allocate(states(1:nstates))
    do i = 1, nstates
       call nullify(states(i))
       call set_option_path(states(i), "/material_phase["//int2str(i-1)//"]")
    end do

    call initialise_ocean_forcing_readers

    call insert_external_mesh(states, save_vtk_cache = .true.)

    call insert_derived_meshes(states)

    !If any meshes have constraints, allocate an appropriate trace mesh
    call insert_trace_meshes(states)

    call compute_domain_statistics(states)

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
    call profiler_toc("I/O")
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
         & mesh_file_format, from_file_path
    integer, dimension(:), pointer :: coplanar_ids
    integer, dimension(3) :: mesh_dims
    integer :: i, j, nmeshes, nstates, quad_degree, stat
    type(element_type), pointer :: shape
    type(quadrature_type), pointer :: quad
    logical :: from_file, extruded
    integer :: dim, mdim, loc, column_ids
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

       from_file_path = trim(mesh_path) // "/from_file"
       from_file = have_option(from_file_path)
       if (.not. from_file) then
         from_file_path = trim(mesh_path) // "/from_mesh/extrude/checkpoint_from_file"
         extruded = have_option(from_file_path)
       else
         extruded = .false.
       end if

       if(from_file .or. extruded) then

          ! Get file format
          ! Can remove stat test when mesh format data backwards compatibility is removed
          call get_option(trim(from_file_path)//"/format/name", mesh_file_format, stat)
          ! Can remove following when mesh format data backwards compatibility is removed
          if(stat /= 0) then
             ewrite(0, *) "Warning: Mesh format name attribute missing for mesh " // trim(mesh_path)
             call get_option(trim(from_file_path)//"/format", mesh_file_format)
          end if

          ! Get filename for mesh, and other options
          call get_option(trim(from_file_path)//"/file_name", mesh_file_name)
          call get_option("/geometry/quadrature/degree", quad_degree)
          quad_family = get_quad_family()

          ! to make sure that the dimension is set even if MPI is not being used
          call get_option('/geometry/dimension', dim)

          if (is_active_process) then
            select case (mesh_file_format)
            case ("triangle", "gmsh", "exodusii")
              ! Get mesh dimension if present
              call get_option(trim(mesh_path)//"/from_file/dimension", mdim, stat)
              ! Read mesh
              if(stat==0) then
                 position=read_mesh_files(trim(mesh_file_name), &
                      quad_degree=quad_degree, &
                      quad_family=quad_family, mdim=mdim, &
                      format=mesh_file_format)
              else
                 position=read_mesh_files(trim(mesh_file_name), &
                      quad_degree=quad_degree, &
                      quad_family=quad_family, &
                      format=mesh_file_format)
              end if
              ! After successfully reading in an ExodusII mesh, change the option
              ! mesh file format to "gmsh", as the write routines for ExodusII are currently
              ! not implemented. Thus, checkpoints etc are dumped as gmsh mesh files
              if (trim(mesh_file_format)=="exodusii") then
                mesh_file_format = "gmsh"
                call set_option_attribute(trim(from_file_path)//"/format/name", trim(mesh_file_format), stat=stat)
                if (stat /= SPUD_NO_ERROR) then
                  FLAbort("Failed to set the mesh format to gmsh (required for checkpointing). Spud error code is: "//int2str(stat))
                end if
              end if
              mesh=position%mesh
            case ("vtu")
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
            case default
              ewrite(-1,*) trim(mesh_file_format), " is not a valid format for a mesh file"
              FLAbort("Invalid format for mesh file")
            end select
         end if

          if (no_active_processes /= getnprocs()) then
            ! not all processes are active, they need to be told the mesh dimensions

            ! receive the mesh dimension from rank 0
            if (getrank()==0) then
              if (is_active_process) then
                ! normally rank 0 should always be active, so it knows the dimensions
                mesh_dims(1)=mesh_dim(mesh)
                mesh_dims(2)=ele_loc(mesh,1)
                if (associated(mesh%columns)) then
                  mesh_dims(3)=1
                else
                  mesh_dims(3)=0
                end if
                ! The coordinate dimension is not the same as the mesh dimension
                ! in the case of spherical shells, and needs to be broadcast as
                ! well.  And this needs to be here to allow for the special case
                ! below
                dim=position%dim
              else
                ! this is a special case for a unit test with 1 inactive process
                call get_option('/geometry/dimension', mesh_dims(1))
                mesh_dims(2)=mesh_dims(1)+1
                mesh_dims(3)=0
                dim = mesh_dims(1)
              end if
            end if
            call MPI_bcast(mesh_dims, 3, getpinteger(), 0, MPI_COMM_FEMTOOLS, stat)
            call MPI_bcast(dim, 1, getpinteger(), 0, MPI_COMM_FEMTOOLS, stat)
          end if

          
          if (.not. is_active_process) then
            ! is_active_process records whether we have data on disk or not
            ! see the comment in Global_Parameters. In this block, 
            ! we want to allocate an empty mesh and positions.

            mdim=mesh_dims(1)
            loc=mesh_dims(2)
            column_ids=mesh_dims(3)

            allocate(quad)
            allocate(shape)
            quad = make_quadrature(loc, mdim, degree=quad_degree, family=quad_family)
            shape=make_element_shape(loc, mdim, 1, quad)
            call allocate(mesh, nodes=0, elements=0, shape=shape, name="EmptyMesh")
            call allocate(position, dim, mesh, "EmptyCoordinate")
            call add_faces(mesh)
            if (column_ids>0) then
              ! the association status of mesh%columns should be collective
              allocate(mesh%columns(1:0))
            end if

            ! Reference counting cleanups.
            call deallocate(mesh)
            call deallocate(quad)
            call deallocate(shape)
            
            deallocate(quad)
            deallocate(shape)

          end if

          ! if there is a derived mesh which specifies periodic bcs 
          ! to be *removed*, we assume the external mesh is periodic
          mesh%periodic = option_count("/geometry/mesh/from_mesh/&
             &periodic_boundary_conditions/remove_periodicity")>0

          ! Get mesh name. This must be done after the mesh file has
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
            if (no_active_processes == 1) then
              call create_empty_halo(position)
            else
              call read_halos(mesh_file_name, position)              
            end if
            ! Local element ordering needs to be consistent between processes, otherwise
            ! code in Halos_Repair (used in halo construction of derived meshes) will fail
            if (.not. verify_consistent_local_element_numbering(position%mesh)) then
              ewrite(-1,*) "The local element ordering is not the same between processes"
              ewrite(-1,*) "that see the same element. This is a necessary condition on the"
              ewrite(-1,*) "decomposed input meshes for fluidity. The fact that you've"
              ewrite(-1,*) "obtained such meshes is likely a bug in fldecomp or the"
              ewrite(-1,*) "checkpointing code. Please report to the fluidity mailing"
              ewrite(-1,*) "list and state exactly how you've obtained your input files."
              FLAbort("Inconsistent local element ordering")
            end if
            mesh = position%mesh
          end if
          
          ! coplanar ids are create here already and stored on the mesh, 
          ! so its derived meshes get the same coplanar ids
          ! (must be done after halo registration)
          if (.not. mesh_periodic(mesh)) then
            ! for periodic meshes, we postpone till we've derived the non-periodic mesh
            call get_coplanar_ids(mesh, position, coplanar_ids)
          end if
          
          if (.not. have_option(trim(mesh_path)//'/exclude_from_mesh_adaptivity')) then
            ! We register this as the topology mesh
            ! this is the mesh used by adaptivity for error measures and such
            ! (it may gets replaced if adding periodicity or extrusion)
            topology_mesh_name = mesh%name
            ! same for the mesh to be handled by adapt_state()
            ! (this gets replaced in case adding periodicity but not by extrusion)
            adaptivity_mesh_name = mesh%name
          end if
                    
          call surface_id_stats(mesh, position)
          
        end if
        
        if (from_file) then
          
          ! Insert mesh and position field into states(1) and
          ! alias it to all the others
          call insert(states, mesh, mesh%name)
          call insert(states, position, position%name)
          call deallocate(position)
          
        else if (extruded) then
          
          ! This will be picked up by insert_derived_meshes and changed
          ! appropriately
          call insert(states, position, "AdaptedExtrudedPositions")
          call deallocate(position)
          
        end if

    end do external_mesh_loop
    
    if(.not. present_and_true(save_vtk_cache)) then
       ! Flush the cache
       call vtk_cache_finalise()
    end if

    call toc(TICTOC_ID_IO_READ)

  end subroutine insert_external_mesh
    
  subroutine insert_derived_meshes(states, skip_extrusion)
    ! Insert derived meshes in state
    type(state_type), intent(inout), dimension(:) :: states
    ! if present and true: skip extrusion of meshes, and insert 0 node dummy meshes 
    ! instead (will have correct shape and dimension)
    logical, optional, intent(in):: skip_extrusion
    
    character(len=FIELD_NAME_LEN) :: mesh_name    
    character(len=OPTION_PATH_LEN) :: mesh_path
    logical :: incomplete, updated
    integer :: i
    integer :: nmeshes
    
    ! Get number of meshes
    nmeshes=option_count("/geometry/mesh")
    periodic_boundary_option_path=""

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

          call insert_derived_mesh(trim(mesh_path), &
                                   trim(mesh_name), &
                                   incomplete, &
                                   updated, &
                                   states, &
                                   skip_extrusion = skip_extrusion)
          
       end do derived_mesh_loop

       ! If we didn't skip any fields then we are done.
       if (.not.incomplete) exit outer_loop

       ! If we did skip fields and didn't update any fields this pass, then
       ! we have unresolvable dependencies.
       if (.not.updated) then
          FLExit("Unresolvable mesh dependencies")
       end if

    end do outer_loop

    ! not really a derived mesh but this is a relatively clean place to set the transform_to_physical
    ! spherical flag so that the main Coordinate field is interpretted as being spherical at the gauss
    ! points
    if (have_option('/geometry/spherical_earth/analytical_mapping/')) then
      call set_analytical_spherical_mapping()
    end if

  end subroutine insert_derived_meshes
           
  subroutine insert_derived_mesh(mesh_path, mesh_name, incomplete, updated, states, skip_extrusion)
          
    ! Insert one derived mesh given by mesh path and mesh_name
       
    character(len=*), intent(in) :: mesh_path
    character(len=*), intent(in) :: mesh_name    
    logical, intent(inout) :: incomplete
    logical, intent(inout) :: updated
    type(state_type), intent(inout), dimension(:) :: states
    ! if present and true: skip extrusion of meshes, and insert 0 node dummy meshes 
    ! instead (will have correct shape and dimension)
    logical, optional, intent(in):: skip_extrusion
   
    type(mesh_type) :: mesh, model_mesh
    type(vector_field), pointer :: position, modelposition
    type(vector_field) :: periodic_position, nonperiodic_position, extrudedposition, coordinateposition
    type(element_type) :: full_shape
    type(quadrature_type) :: quad

    character(len=FIELD_NAME_LEN) :: model_mesh_name
    character(len=OPTION_PATH_LEN) :: shape_type, cont
    logical :: new_cont, extrusion, periodic, remove_periodicity
    logical :: new_shape_type, new_degree, from_shape, make_new_mesh
    integer :: from_degree, from_shape_type, from_cont, j, stat
    integer :: quadrature_degree, h_dim
    logical :: exclude_from_mesh_adaptivity

    if (has_mesh(states(1), mesh_name)) then
       ! We already did this one.
       return
    end if
          
    if(have_option(trim(mesh_path)//"/from_mesh")) then
             
       ! Get model mesh name
       call get_option(trim(mesh_path)//"/from_mesh/mesh[0]/name", model_mesh_name)
             
       ! Extract model mesh
       model_mesh=extract_mesh(states(1), trim(model_mesh_name), stat=stat)
       if (stat/=0) then
          ! The mesh from which this mesh is derived is not yet
          ! present.
          incomplete=.true.
          return
       end if
             
       ! Find out if the new mesh is different from the old mesh and if
       ! so, find out how it differs - in the options check
       ! we've made sure only one of those (or both new_shape and new_cont) are .true.
       ! If there are no differences, do not create new mesh.
       from_shape=have_option(trim(mesh_path)//"/from_mesh/mesh_shape")

       ! 1. If mesh shape options are specified, check if they are different to the model mesh.
       if (from_shape) then
         ! 1.1. Check polynomial_degree option
         call get_option(trim(mesh_path)//"/from_mesh/mesh_shape/polynomial_degree", &
              from_degree, stat)
         if(stat==0) then
           ! Is polynomial_degree the same as model mesh?
           if(from_degree==model_mesh%shape%degree) then
             new_degree=.false.
           else
             new_degree=.true.
           end if
         ! If degree is not specified, use the model mesh degree.
         else
           new_degree=.false.
         end if

         ! 1.2. Check element_type option
         call get_option(trim(mesh_path)//"/from_mesh/mesh_shape/element_type", &
              shape_type, stat)
         if(stat==0) then
           ! Set comparison variable from_shape_type
           if(trim(shape_type)=="lagrangian") then
             from_shape_type=ELEMENT_LAGRANGIAN
           else if(trim(shape_type)=="bubble") then
             from_shape_type=ELEMENT_BUBBLE
          else if(trim(shape_type)=="trace") then
             from_shape_type=ELEMENT_TRACE
           end if
           ! If new_shape_type does not match model mesh shape type, make new mesh.
           if(from_shape_type == model_mesh%shape%numbering%type) then
             new_shape_type=.false.
           else
             new_shape_type=.true.
           end if
         ! If no element_type is specified, assume it is the same as model mesh
         ! and do not create new mesh.
         else
           new_shape_type=.false.
         end if
       ! Else if no mesh shape options are set, do not make new mesh.
       else
         new_degree=.false.; new_shape_type=.false.
       end if

       ! 2. If mesh_continuity is specified, check if it is different to the model mesh.
       call get_option(trim(mesh_path)//"/from_mesh/mesh_continuity", cont, stat)
       if(stat==0) then
         if(trim(cont)=="discontinuous") then
           from_cont=-1
         else if(trim(cont)=="continuous") then
           from_cont=0
         end if
         ! 2.1. If continuity is not the same as model mesh, create new mesh.
         if(from_cont==model_mesh%continuity) then
           new_cont=.false.
         else
           new_cont=.true.
         end if
       ! If no continuity is specified, assume it is the same as model mesh,
       ! and do not create a new mesh.
       else
         new_cont=.false.
       end if

       ! 3. If any of the above are true, make new mesh.
       make_new_mesh = new_shape_type .or. new_degree .or. new_cont

       extrusion=have_option(trim(mesh_path)//"/from_mesh/extrude")
       periodic=have_option(trim(mesh_path)//"/from_mesh/periodic_boundary_conditions")
       exclude_from_mesh_adaptivity=have_option(trim(mesh_path)//"/exclude_from_mesh_adaptivity")

       if (periodic) then
         ! there is an options check to guarantee that all periodic bcs have remove_periodicity
         remove_periodicity=option_count(trim(mesh_path)//"/from_mesh/periodic_boundary_conditions/remove_periodicity")>0
         if (remove_periodicity) then
           if (.not. mesh_periodic(model_mesh)) then
              ewrite(0,*) "In derivation of mesh ", trim(mesh_name), " from ", trim(model_mesh_name)
              FLExit("Trying to remove periodic bcs from non-periodic mesh.")
           end if
         end if
       end if
             
       ! We added at least one mesh on this pass.
       updated=.true.
                
       if (extrusion) then

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
                  
            modelposition => extract_vector_field(states(1), trim(model_mesh_name)//"Coordinate")
             
            if (present_and_true(skip_extrusion)) then

              ! the dummy mesh does need a shape of the right dimension
              h_dim = mesh_dim(modelposition)
              call get_option("/geometry/quadrature/degree", quadrature_degree)
              quad = make_quadrature(vertices=h_dim + 2, dim=h_dim + 1, degree=quadrature_degree)
              full_shape = make_element_shape(vertices=h_dim + 2, dim=h_dim + 1, degree=1, quad=quad)
              call deallocate(quad)

              call allocate(mesh, nodes=0, elements=0, shape=full_shape, name=mesh_name)
              call deallocate(full_shape)
              allocate(mesh%columns(1:0))
              call add_faces(mesh)
              mesh%periodic=modelposition%mesh%periodic
              call allocate(extrudedposition, h_dim+1, mesh, "EmptyCoordinate") ! name is fixed below
              call deallocate(mesh)
              if (IsParallel()) call create_empty_halo(extrudedposition)
            else if (have_option('/geometry/spherical_earth/')) then
              call extrude_radially(modelposition, mesh_path, extrudedposition)
            else
              call extrude(modelposition, mesh_path, extrudedposition)
            end if
                
         end if
               
         mesh = extrudedposition%mesh
               
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
               
       else if (make_new_mesh) then

         mesh = make_mesh_from_options(model_mesh, mesh_path)
             
       else if (periodic) then
                
         if (remove_periodicity) then
           ! model mesh can't be the CoordinateMesh:
           periodic_position=extract_vector_field(states(1), trim(model_mesh_name)//"Coordinate")
           nonperiodic_position = make_mesh_unperiodic_from_options( &
             periodic_position, mesh_path)
                
           ! the positions of this mesh have to be stored now 
           ! as it cannot be interpolated later.
           if (mesh_name=="CoordinateMesh") then
              nonperiodic_position%name = "Coordinate"
           else
              nonperiodic_position%name = trim(mesh_name)//"Coordinate"
           end if
           call insert(states, nonperiodic_position, nonperiodic_position%name)
           call deallocate(nonperiodic_position)
                 
           mesh=nonperiodic_position%mesh
           call incref(mesh)
                 
         else
           ! this means we can only periodise a mesh with an associated position field
           if (trim(model_mesh_name) == "CoordinateMesh") then
             position => extract_vector_field(states(1), 'Coordinate')
           else
             position => extract_vector_field(states(1), trim(model_mesh_name)//'Coordinate')
           end if
           periodic_position = make_mesh_periodic_from_options(position, mesh_path)
           ! Ensure the name and option path are set on the original
           ! mesh descriptor.
           periodic_position%mesh%name = mesh_name
           periodic_position%mesh%option_path = trim(mesh_path)

           mesh = periodic_position%mesh
           call incref(mesh)
           call insert(states, periodic_position, trim(periodic_position%name))
           call deallocate(periodic_position)
         end if
               
       else
          ! copy mesh unchanged, new reference
          mesh=model_mesh                
          call incref(mesh)
                
       end if
             
       mesh%name = mesh_name
                
       ! Set mesh option path.
       mesh%option_path = trim(mesh_path)
             
       ! if this is the coordinate mesh then we should insert the coordinate field
       ! also meshes excluded from adaptivity all have their own coordinate field
       ! for extrusion and periodic: the coordinate field has already been inserted above
       if ((trim(mesh_name)=="CoordinateMesh" .or. exclude_from_mesh_adaptivity) &
          .and. .not. (extrusion .or. periodic)) then

          if (model_mesh_name=="CoordinateMesh") then
             modelposition => extract_vector_field(states(1), "Coordinate")
          else
             modelposition => extract_vector_field(states(1), trim(model_mesh_name)//"Coordinate")
          end if
                
          if (mesh_name=="CoordinateMesh") then
            call allocate(coordinateposition, modelposition%dim, mesh, "Coordinate")
          else
            call allocate(coordinateposition, modelposition%dim, mesh, trim(mesh_name)//"Coordinate")
          end if
                
          ! remap the external mesh positions onto the CoordinateMesh... this requires that the space
          ! of the coordinates spans that of the external mesh
          call remap_field(from_field=modelposition, to_field=coordinateposition)
                
          if (mesh_name=="CoordinateMesh" .and. have_option('/geometry/spherical_earth/')) then

            if (have_option('/geometry/spherical_earth/superparametric_mapping/')) then
              call higher_order_sphere_projection(modelposition, coordinateposition)
            end if

          endif

          ! insert into states(1) and alias to all others
          call insert(states, coordinateposition, coordinateposition%name)
          ! drop reference to the local copy of the Coordinate field
          call deallocate(coordinateposition)
       end if
       if (trim(mesh_name)=="CoordinateMesh" .and. mesh_periodic(mesh)) then
         FLExit("CoordinateMesh may not be periodic")
       end if             

       ! Insert mesh into all states
       call insert(states, mesh, mesh%name)

       if (.not. have_option(trim(mesh_path)//'/exclude_from_mesh_adaptivity')) then
         ! update info for adaptivity/error metric code:
               
         if (extrusion .or. (periodic .and. .not. remove_periodicity)) then
           ! this is the name of the mesh to be used by the error metric for adaptivity
           topology_mesh_name=mesh%name
         end if
               
         if ((extrusion.and..not.have_option('/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity')) &
             .or.(periodic .and. .not. remove_periodicity)) then
           ! this is the name of the mesh to be adapted by adaptivity
           adaptivity_mesh_name=mesh%name
         end if

         if (periodic .and. trim(periodic_boundary_option_path(mesh%shape%dim)) == "") then
           periodic_boundary_option_path(mesh%shape%dim) = trim(mesh_path)
         end if
       end if
             
       call deallocate(mesh)
             
    end if
       
  end subroutine insert_derived_mesh

  subroutine insert_trace_meshes(states)
    !If any meshes have constraints, allocate an appropriate trace mesh
    type(state_type), dimension(:), intent(inout) :: states
    !
    type(mesh_type) :: from_mesh, model_mesh
    type(mesh_type) :: trace_mesh
    type(quadrature_type) :: quad
    type(element_type) :: trace_shape
    integer :: mesh_no, trace_degree, dim, loc, constraint_choice, &
         &quad_degree
    logical :: allocate_trace_mesh
    character(len=FIELD_NAME_LEN) :: model_mesh_name

    do mesh_no = 1, mesh_count(states(1))
       allocate_trace_mesh = .false.
       from_mesh = extract_mesh(states(1),mesh_no)
       if(associated(from_mesh%shape%constraints)) then
          constraint_choice = from_mesh%shape%constraints%type
          if(constraint_choice.ne.CONSTRAINT_NONE) then
             select case(constraint_choice)
             case (CONSTRAINT_BDM)
                trace_degree = from_mesh%shape%degree
             case (CONSTRAINT_RT)
                trace_degree = from_mesh%shape%degree-1
             case (CONSTRAINT_BDFM)
                trace_degree = from_mesh%shape%degree-1
             case default
                FLAbort('Constraint type not supported')
             end select
             dim = from_mesh%shape%dim
             loc=from_mesh%shape%quadrature%vertices

             ! Get model mesh name
             call get_option("/geometry/mesh["//int2str(mesh_no)//&
                  &"]/from_mesh/mesh[0]/name",&
                  &model_mesh_name)
             
             ! Extract model mesh
             model_mesh=extract_mesh(states(1), trim(model_mesh_name))
             
             !Make quadrature
             call get_option("/geometry/quadrature/degree",&
                  & quad_degree)
             quad=make_quadrature(loc, dim, &
                  degree=quad_degree, family=get_quad_family())
             !allocate shape
             trace_shape=make_element_shape(loc, dim, trace_degree, &
                  &quad,type=ELEMENT_TRACE)
             !deallocate quadrature (just drop a reference)
             call deallocate(quad)
             !allocate mesh
             trace_mesh=make_mesh(model_mesh, trace_shape, continuity=-1,&
                  name=trim(from_mesh%name)//"Trace")
             !deallocate shape (just drop a reference)
             call deallocate(trace_shape)
             !insert into states
             call insert(states,trace_mesh,trace_mesh%name)
             !deallocate mesh (just drop a reference)
             call deallocate(trace_mesh)
          end if
       end if
    end do

  end subroutine insert_trace_meshes

  function make_mesh_from_options(from_mesh, mesh_path) result (mesh)
    ! make new mesh changing shape or continuity of from_mesh
    type(mesh_type):: mesh
    type(mesh_type), intent(in):: from_mesh
    character(len=*), intent(in):: mesh_path
    
    character(len=FIELD_NAME_LEN) :: mesh_name
    character(len=OPTION_PATH_LEN) :: continuity_option, element_option, constraint_option_string
    type(quadrature_type):: quad
    type(element_type):: shape
    integer :: constraint_choice
    integer:: loc, dim, poly_degree, continuity, new_shape_type, quad_degree, stat
    logical :: new_shape
    
    ! Get new mesh shape information
    
    new_shape = have_option(trim(mesh_path)//"/from_mesh/mesh_shape")
    if(new_shape) then
      ! Get new mesh element type
      call get_option(trim(mesh_path)//"/from_mesh/mesh_shape/element_type", &
                      element_option, stat)
      if(stat==0) then
        if(trim(element_option)=="lagrangian") then
           new_shape_type=ELEMENT_LAGRANGIAN
        else if(trim(element_option)=="bubble") then
           new_shape_type=ELEMENT_BUBBLE
        else if(trim(element_option)=="trace") then
           new_shape_type=ELEMENT_TRACE
        end if
      else
        new_shape_type=from_mesh%shape%numbering%type
      end if
      
      ! degree is the degree of the Lagrange polynomials (even if you add in a bubble function)
      call get_option(trim(mesh_path)//"/from_mesh/mesh_shape/polynomial_degree", &
                      poly_degree, default=from_mesh%shape%degree)
    
      ! loc is the number of vertices of the element
      loc=from_mesh%shape%loc
      ! dim is the dimension
      dim=from_mesh%shape%dim
      ! Make quadrature
      call get_option("/geometry/quadrature/degree",&
           & quad_degree)
      quad=make_quadrature(loc, dim, degree=quad_degree, family=get_quad_family())
      ! Get element constraints
      call get_option(trim(mesh_path)//"/from_mesh/constraint_type",&
           constraint_option_string, stat)
      if(stat==0) then
         if(trim(constraint_option_string)=="BDFM") then
            constraint_choice=CONSTRAINT_BDFM
         else if(trim(constraint_option_string)=="RT") then
            constraint_choice=CONSTRAINT_RT
         else if(trim(constraint_option_string)=="BDM") then
            constraint_choice=CONSTRAINT_BDM
         else if(trim(constraint_option_string)=="none") then
            constraint_choice=CONSTRAINT_NONE
         end if
      else
         constraint_choice = CONSTRAINT_NONE
      end if
      
      ! Make new mesh shape
      shape=make_element_shape(loc, dim, poly_degree, quad,&
           &type=new_shape_type,constraint_type_choice=constraint_choice)
      call deallocate(quad) ! Really just drop a reference.
    else
      shape=from_mesh%shape
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
      continuity=from_mesh%continuity
    end if

    ! Get mesh name.
    call get_option(trim(mesh_path)//"/name", mesh_name)

    ! Make new mesh
    mesh=make_mesh(from_mesh, shape, continuity, mesh_name)

    ! Set mesh option path
    mesh%option_path = trim(mesh_path)

    ! Drop one reference to shape
    call deallocate(shape)
               
  end function make_mesh_from_options
  
  function make_mesh_periodic_from_options(position, mesh_path) result (position_out)
    ! make a periodic mesh as specified by options
    type(vector_field):: position_out
    type(vector_field), intent(in):: position
    character(len=*), intent(in):: mesh_path

    
    type(vector_field):: from_position
    type(integer_hash_table):: periodic_face_map
    character(len=FIELD_NAME_LEN):: bc_name, mesh_name
    character(len=OPTION_PATH_LEN) :: periodic_mapping_python
    integer, dimension(:), allocatable :: physical_boundary_ids, aliased_boundary_ids
    integer, dimension(2) :: shape_option
    integer:: n_periodic_bcs
    integer:: j
    logical :: fiddled_with_faces
    
    assert(has_faces(position%mesh))
    
    from_position=position
    
    ! builds up a map from aliased to physical faces
    call allocate(periodic_face_map)
    
    n_periodic_bcs=option_count(trim(mesh_path)//"/from_mesh/periodic_boundary_conditions")
    ewrite(2,*) "n_periodic_bcs=", n_periodic_bcs
    call incref(from_position)
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
       

       fiddled_with_faces = .false.
       if (.not. has_faces(from_position%mesh)) then
         from_position%mesh%faces => position%mesh%faces
         fiddled_with_faces = .true.
       end if
       position_out=make_mesh_periodic(from_position,&
          physical_boundary_ids,aliased_boundary_ids, &
          periodic_mapping_python, periodic_face_map=periodic_face_map)
       if (fiddled_with_faces) then
         from_position%mesh%faces => null()
       end if
       call deallocate(from_position)
       from_position=position_out
       
       deallocate( physical_boundary_ids, aliased_boundary_ids )
    end do
    
    call add_faces(position_out%mesh, model=position%mesh, periodic_face_map=periodic_face_map)
    
    call deallocate(periodic_face_map)
    
    ! finally fix the name of the produced mesh and its coordinate field
    call get_option(trim(mesh_path)//'/name', mesh_name)
    position_out%mesh%name=mesh_name
    if (mesh_name=="CoordinateMesh") then
      position_out%name="Coordinate"
    else
      position_out%name=trim(mesh_name)//"Coordinate"
    end if
    
  end function make_mesh_periodic_from_options

  function make_mesh_unperiodic_from_options(from_position, mesh_path, aliased_to_new_node_number, stat) result (position)
    ! make a periodic mesh as specified by options
    type(vector_field):: position
    type(vector_field), intent(in):: from_position
    character(len=*), intent(in):: mesh_path
    integer, intent(out), optional :: stat
    type(integer_hash_table), optional, intent(out) :: aliased_to_new_node_number

    type(vector_field):: lfrom_position, nonperiodic_position
    character(len=FIELD_NAME_LEN):: bc_name, mesh_name
    character(len=OPTION_PATH_LEN) :: periodic_mapping_python
    integer, dimension(:), allocatable :: physical_boundary_ids, aliased_boundary_ids
    integer, dimension(2) :: shape_option
    integer:: n_periodic_bcs
    integer:: j
    type(integer_hash_table) :: laliased_to_new_node_number
    type(integer_set) :: all_periodic_bc_ids
    logical :: fiddled_with_faces

    if (present(stat)) then
      stat = 0
    end if

    ! Get mesh name.
    call get_option(trim(mesh_path)//"/name", mesh_name)
    
    ! get our own reference of from_position, that we can throw away again
    lfrom_position=from_position
    call incref(lfrom_position)
    
    n_periodic_bcs=option_count(trim(mesh_path)//"/from_mesh/periodic_boundary_conditions")
    ewrite(2,*) "n_periodic_bcs=", n_periodic_bcs
    if (n_periodic_bcs == 0) then
      ewrite(-1,*) "You almost certainly didn't mean to pass in this option path."
      ewrite(-1,*) "trim(mesh_path): ", trim(mesh_path)
      ewrite(-1,*) "mesh_name: ", trim(mesh_name)
      FLAbort("No periodic boundary conditions to unwrap!")
    end if

    call allocate(all_periodic_bc_ids)
    do j=0, n_periodic_bcs-1
       shape_option = option_shape(trim(mesh_path)//"/from_mesh/periodic_boundary_conditions["//int2str(j)//"]/physical_boundary_ids")
       allocate( physical_boundary_ids(shape_option(1)) )
       call get_option(trim(mesh_path)//"/from_mesh/periodic_boundary_conditions["//int2str(j)//"]/physical_boundary_ids",physical_boundary_ids)
       call insert(all_periodic_bc_ids, physical_boundary_ids)
       deallocate(physical_boundary_ids)

       shape_option = option_shape(trim(mesh_path)//"/from_mesh/periodic_boundary_conditions["//int2str(j)//"]/aliased_boundary_ids")
       allocate( aliased_boundary_ids(shape_option(1)) )
       call get_option(trim(mesh_path)//"/from_mesh/periodic_boundary_conditions["//int2str(j)//"]/aliased_boundary_ids",aliased_boundary_ids)
       call insert(all_periodic_bc_ids, aliased_boundary_ids)
       deallocate(aliased_boundary_ids)
    end do

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
       
       ewrite(2,*) 'Removing periodicity from mesh'
       
       fiddled_with_faces = .false.
       if (.not. has_faces(lfrom_position%mesh)) then
         lfrom_position%mesh%faces => from_position%mesh%faces
         fiddled_with_faces = .true.
       end if

       nonperiodic_position=make_mesh_unperiodic(lfrom_position,&
          physical_boundary_ids,aliased_boundary_ids, &
          periodic_mapping_python, mesh_name, all_periodic_bc_ids, laliased_to_new_node_number)

       if (fiddled_with_faces) then
         lfrom_position%mesh%faces => null()
       end if

       if (associated(lfrom_position%mesh%halos)) then
         assert(associated(lfrom_position%mesh%element_halos))
         call derive_nonperiodic_halos_from_periodic_halos(nonperiodic_position, lfrom_position, laliased_to_new_node_number)
       end if
       call deallocate(lfrom_position)
       if (present(aliased_to_new_node_number)) then
         aliased_to_new_node_number = laliased_to_new_node_number
       else
         call deallocate(laliased_to_new_node_number)
       end if
       lfrom_position=nonperiodic_position
       
       deallocate( physical_boundary_ids, aliased_boundary_ids )
    end do
      
    ! assumes all periodic bcs have been removed
    ! this is checked for in add_faces
    ! this flag needs setting before the call to add_faces
    nonperiodic_position%mesh%periodic=.false.

    assert(associated(nonperiodic_position%mesh%shape%numbering))
    
    if (has_faces(from_position%mesh)) then
      call add_faces(nonperiodic_position%mesh, model=from_position%mesh, stat=stat)
    end if
    
    position=nonperiodic_position

    call deallocate(all_periodic_bc_ids)
    
  end function make_mesh_unperiodic_from_options
  
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

    character(len=OPTION_PATH_LEN) :: field_name, absolute_path
    integer :: i, istate ! counters
    integer :: nstates ! number of states
    character(len=255) :: tmp ! temporary string to make life a little easier
    type(scalar_field), pointer :: fshistory_sfield
    integer :: fshistory_levels 
    
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
            states(1), dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    end if

    ! Field that controls the weighting of partitions:
    if (have_option('/flredecomp/field_weighted_partitions')) then
       call allocate_and_insert_scalar_field('/flredecomp/field_weighted_partitions/scalar_field::FieldWeightedPartitionValues', states(1))
    end if
    
    if (have_option('/mesh_adaptivity/hr_adaptivity/zoltan_options/field_weighted_partitions')) then
       call allocate_and_insert_scalar_field('/mesh_adaptivity/hr_adaptivity/zoltan_options/field_weighted_partitions/scalar_field::FieldWeightedPartitionValues', states(1))
    end if
    
    ! grid velocity
    if (have_option('/mesh_adaptivity/mesh_movement/vector_field::GridVelocity')) then
       call allocate_and_insert_vector_field('/mesh_adaptivity/mesh_movement/vector_field::GridVelocity', &
            states(1), dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    end if

    ! solar irradiance submodel (hyperlight)
    if (have_option("/ocean_biology/lagrangian_ensemble/hyperlight")) then 
       call allocate_and_insert_irradiance(states(1))
    end if

    ! insert electrical property fields
    do i=1,nstates
      tmp = '/material_phase['//int2str(i-1)//']/electrical_properties/coupling_coefficients/'
      ! Electrokinetic coupling coefficient scalar field
      if (have_option(trim(tmp)//'scalar_field::Electrokinetic')) then
        call allocate_and_insert_scalar_field(trim(tmp)//'scalar_field::Electrokinetic', &
                                              states(i), &
                                              field_name='Electrokinetic')
      end if
      ! Thermoelectric coupling coefficient scalar field
      if (have_option(trim(tmp)//'scalar_field::Thermoelectric')) then
        call allocate_and_insert_scalar_field(trim(tmp)//'scalar_field::Thermoelectric', &
                                              states(i), &
                                              field_name='Thermoelectric')
      end if
      ! Electrochemical coupling coefficient scalar field
      if (have_option(trim(tmp)//'scalar_field::Electrochemical')) then
        call allocate_and_insert_scalar_field(trim(tmp)//'scalar_field::Electrochemical', &
                                              states(i), &
                                              field_name='Electrochemical')
      end if
    end do

    ! Harmonic Analysis History fields
    if (has_scalar_field(states(1),'FreeSurfaceHistory') ) then
      fshistory_sfield => extract_scalar_field(states(1), 'FreeSurfaceHistory')
      ! levels: the number of levels which will be saved. Too old levels will be overwritten by new ones.
      if (have_option(trim(complete_field_path(fshistory_sfield%option_path)) // "/algorithm/levels")) then
        call get_option(trim(complete_field_path(fshistory_sfield%option_path)) // "/algorithm/levels", fshistory_levels)
        fshistory_levels=max(fshistory_levels,0)
      else
        fshistory_levels=50
      end if
      do i=1,fshistory_levels
          call allocate_and_insert_scalar_field('', states(1), parent_mesh='PressureMesh', field_name='harmonic'//int2str(i))
      end do
    end if

    ! insert miscellaneous scalar fields
    do i=1, size(additional_fields_absolute)
       if (have_option(trim(additional_fields_absolute(i)))) then

          call allocate_and_insert_one_phase(additional_fields_absolute(i), states(1), &
             dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
          
       end if
    end do
    
    do i=1, size(additional_fields_relative)
       do istate = 1, size(states)
         absolute_path = "/material_phase["//int2str(istate-1)//"]/"//trim(additional_fields_relative(i))
         if (have_option(absolute_path)) then

            call allocate_and_insert_one_phase(absolute_path, states(istate), &
               dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
            
         end if
       end do
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

    end subroutine allocate_and_insert_one_phase

    subroutine allocate_and_insert_irradiance(state)
      ! Allocate irradiance fields for 36 wavebands in PAR
      type(state_type), intent(inout) :: state
      integer :: j
      real :: lambda
      character(len=OPTION_PATH_LEN) :: light_path, field_name

      ! Replicate irradiance template field for all wavebands
      light_path = "/ocean_biology/lagrangian_ensemble/hyperlight"      
      frequency_field_loop: do j=0,35         
         lambda = 350.0 + (j * 10.0)
         field_name="Irradiance_"//int2str(NINT(lambda))
         call allocate_and_insert_scalar_field(&
                  trim(light_path)&
                  //"/scalar_field::IrradianceTemplate", &
                  state, field_name=trim(field_name), &
                  dont_allocate_prognostic_value_spaces&
                  =dont_allocate_prognostic_value_spaces)
      end do frequency_field_loop

      ! Create PAR irradiance field
      if (have_option("/ocean_biology/lagrangian_ensemble/hyperlight/scalar_field::IrradiancePAR")) then 
         call allocate_and_insert_scalar_field(&
                  trim(light_path)&
                  //"/scalar_field::IrradiancePAR", &
                  state, field_name="IrradiancePAR", &
                  dont_allocate_prognostic_value_spaces&
                  =dont_allocate_prognostic_value_spaces)
      end if
    end subroutine allocate_and_insert_irradiance

  end subroutine allocate_and_insert_fields

  subroutine alias_fields(states)
    type(state_type), dimension(:), intent(inout) :: states

    character(len=OPTION_PATH_LEN) :: path
    character(len=OPTION_PATH_LEN) :: state_name, aliased_field_name, field_name
    integer :: stat
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
    ! distance to top and bottom
    if (have_option('/geometry/ocean_boundaries')) then
    
       sfield = extract_scalar_field(states(1), 'DistanceToTop')
       sfield%aliased = .true.
       do i = 1,nstates-1
          call insert(states(i+1), sfield, 'DistanceToTop')
       end do
       
       sfield = extract_scalar_field(states(1), 'DistanceToBottom')
       sfield%aliased = .true.
       do i = 1,nstates-1
          call insert(states(i+1), sfield, 'DistanceToBottom')
       end do

    end if

    ! direction of gravity
    if (have_option('/physical_parameters/gravity/vector_field::GravityDirection')) then
       vfield=extract_vector_field(states(1), 'GravityDirection')
       vfield%aliased = .true.
       do i = 1,nstates-1

          call insert(states(i+1), vfield, 'GravityDirection')

       end do
    end if

    ! grid velocity
    if (have_option('/mesh_adaptivity/mesh_movement/vector_field::GridVelocity')) then
    
       ! Save path to field
       path="/mesh_adaptivity/mesh_movement/vector_field::GridVelocity"

       ! If field is aliased, find which field it is aliased to, extract that field from the correct state and insert into state(1)
       is_aliased=have_option(trim(path)//"/aliased")
       if(is_aliased) then
           call get_option(trim(path)//"/name", field_name)
           call get_option(trim(path)//"/aliased/material_phase_name", state_name)
           call get_option(trim(path)//"/aliased/field_name", aliased_field_name)

           k=get_state_index(states, trim(state_name))
           vfield=extract_vector_field(states(k), trim(aliased_field_name))
           vfield%name = trim(field_name)  ! this seems to be necessary to preserve the aliased field's original name
           vfield%aliased = .true.
           call insert(states(1), vfield, trim(field_name))

       end if
    
       vfield=extract_vector_field(states(1), 'GridVelocity')
       vfield%aliased = .true.
       do i = 1,nstates-1

          call insert(states(i+1), vfield, 'GridVelocity')

       end do
    end if

    ! Deal with subgridscale parameterisations.
    call alias_diffusivity(states)
    
    ! Porous media fields
    have_porous_media: if (have_option('/porous_media')) then
       
       ! alias the Porosity field
       sfield=extract_scalar_field(states(1), 'Porosity')
       sfield%aliased = .true.
       do i = 1,nstates-1
          call insert(states(i+1), sfield, 'Porosity')
       end do
       
       ! alias the AbsolutePermeability field which may be 
       ! either scalar or vector (if present)
       
       sfield=extract_scalar_field(states(1), 'AbsolutePermeability', stat = stat)
       if (stat == 0) then       
          sfield%aliased = .true.
          do i = 1,nstates-1
             call insert(states(i+1), sfield, 'AbsolutePermeability')
          end do       
       end if
       
       vfield=extract_vector_field(states(1), 'AbsolutePermeability', stat = stat)
       if (stat == 0) then       
          vfield%aliased = .true.
          do i = 1,nstates-1
             call insert(states(i+1), vfield, 'AbsolutePermeability')
          end do       
       end if
    
    end if have_porous_media

    ! Porous media dual fields
    have_porous_media_dual: if (have_option('/porous_media_dual')) then
       
       ! alias the PorosityDual field
       sfield=extract_scalar_field(states(1), 'PorosityDual')
       sfield%aliased = .true.
       do i = 1,nstates-1
          call insert(states(i+1), sfield, 'PorosityDual')
       end do
       
       ! alias the AbsolutePermeabilityDual field if present
       sfield=extract_scalar_field(states(1), 'AbsolutePermeabilityDual', stat=stat)
       if (stat == 0) then
          sfield%aliased = .true.
          do i = 1,nstates-1
             call insert(states(i+1), sfield, 'AbsolutePermeabilityDual')
          end do       
       end if
       
    end if have_porous_media_dual

  end subroutine alias_fields

  subroutine alias_diffusivity(states)
    !!< Where fields get their diffusivity from a subgridscale
    !!< parameterisation, it is necessary to alias their diffusivity to the
    !!< diffusivity provided by the parameterisation.
    !!<
    !!< At this stage only prescribed diffusivity, the Generic Length Scale ocean model
    !!< and the K-Epsilon turbulence model are handled via this route.
    !!< Mellor-Yamada is pending a rewrite.
    type(state_type), dimension(:), intent(inout) :: states
    type(scalar_field), pointer :: sfield
    type(tensor_field) :: tfield

    integer :: i, s, stat

    ! Prescribed diffusivity
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

    ! Eddy diffusivity from Generic Length Scale Ocean model
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

  function allocate_scalar_field_as_constant(option_path) result(is_constant)
    !!< Return whether the supplied option path signals a constant
    !!< field

    character(len = *), intent(in) :: option_path

    logical :: is_constant

    is_constant = .false.
    if(option_count(trim(option_path) // "/prescribed/value") == 1) then
       is_constant = have_option(trim(option_path) // "/prescribed/value[0]/constant")
    end if

  end function allocate_scalar_field_as_constant
  
  function allocate_vector_field_as_constant(option_path) result(is_constant)
    !!< Return whether the supplied option path signals a constant
    !!< field

    character(len = *), intent(in) :: option_path

    logical :: is_constant

    is_constant = .false.
    if(option_count(trim(option_path) // "/prescribed/value") == 1) then
       is_constant = have_option(trim(option_path) // "/prescribed/value[0]/constant")
    end if

  end function allocate_vector_field_as_constant
  
  function allocate_tensor_field_as_constant(option_path) result(is_constant)
    !!< Return whether the supplied option path signals a constant
    !!< field

    character(len = *), intent(in) :: option_path

    logical :: is_constant

    if(option_count(trim(option_path) // "/prescribed/value") == 1) then
      is_constant = have_option(trim(option_path) // "/prescribed/value/isotropic/constant") .or. &
        & have_option(trim(option_path) // "/prescribed/value/anisotropic_symmetric/constant") .or. &
        & have_option(trim(option_path) // "/prescribed/value/anisotropic_asymmetric/constant")
    else
      is_constant = .false.
    end if

  end function allocate_tensor_field_as_constant

  function allocate_field_as_constant_scalar(s_field) result(is_constant)
    !!< Return whether the options tree defines the supplied scalar field to
    !!< be constant

    type(scalar_field), intent(in) :: s_field

    logical :: is_constant

    is_constant = allocate_scalar_field_as_constant(s_field%option_path)

  end function allocate_field_as_constant_scalar

  function allocate_field_as_constant_vector(v_field) result(is_constant)
    !!< Return whether the options tree defines the supplied vector field to
    !!< be constant

    type(vector_field), intent(in) :: v_field

    logical :: is_constant

    is_constant = allocate_vector_field_as_constant(v_field%option_path)

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
      is_constant = allocate_tensor_field_as_constant(t_field%option_path)
    end if

  end function allocate_field_as_constant_tensor

  recursive subroutine allocate_and_insert_scalar_field(option_path, state, &
    parent_mesh, parent_name, field_name, &
    dont_allocate_prognostic_value_spaces)

    character(len=*), intent(in) :: option_path
    type(state_type), intent(inout) :: state
    character(len=*), intent(in), optional :: parent_mesh
    character(len=*), intent(in), optional :: parent_name
    character(len=*), optional, intent(in):: field_name
    logical, optional, intent(in):: dont_allocate_prognostic_value_spaces

    logical :: is_prognostic, is_prescribed, is_diagnostic, is_aliased
    ! paths for options and child fields
    character(len=OPTION_PATH_LEN) :: path, adapt_path
    ! Strings for names
    character(len=OPTION_PATH_LEN) :: lfield_name, mesh_name
    type(scalar_field) :: field
    type(mesh_type), pointer :: mesh
    logical :: backward_compatibility, is_constant

    is_aliased=have_option(trim(option_path)//"/aliased")
    if(is_aliased) return

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
    ! Any fields that require backward compatibility are badly behaved, as they
    ! modify constant fields. *Do not add to this list!* Construct an
    ! appropriate diagnostic algorithm instead (possibly an internal).
    backward_compatibility = .false.

    ! Find out what kind of field we have
    is_prognostic=have_option(trim(path)//"/prognostic")
    is_prescribed=have_option(trim(path)//"/prescribed")
    is_diagnostic=have_option(trim(path)//"/diagnostic")

    is_constant=allocate_tensor_field_as_constant(path)

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
    else if(is_constant .and. .not. backward_compatibility) then
         
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
       call allocate_and_insert_scalar_field(adapt_path, state, parent_mesh=topology_mesh_name, &
          parent_name=lfield_name, &
          dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    else if(have_option(trim(adapt_path)//"/relative_measure")) then
       adapt_path=trim(adapt_path)//"/relative_measure/scalar_field::InterpolationErrorBound"
       call allocate_and_insert_scalar_field(adapt_path, state, parent_mesh=topology_mesh_name, &
          parent_name=lfield_name, &
          dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    end if

  end subroutine allocate_and_insert_scalar_field

  recursive subroutine allocate_and_insert_vector_field(option_path, state, parent_mesh, parent_name, &
                                                        field_name, dont_allocate_prognostic_value_spaces)

    character(len=*), intent(in) :: option_path
    type(state_type), intent(inout) :: state
    character(len=*), intent(in), optional :: parent_mesh
    character(len=*), intent(in), optional :: parent_name
    character(len=*), optional, intent(in):: field_name
    logical, intent(in), optional :: dont_allocate_prognostic_value_spaces
    
    integer :: dim
    logical :: is_prognostic, is_prescribed, is_diagnostic, is_aliased
    ! paths for options and child fields
    character(len=OPTION_PATH_LEN) :: path, adapt_path
    ! strings for names
    character(len=OPTION_PATH_LEN) :: lfield_name, mesh_name
    type(mesh_type), pointer :: mesh
    type(vector_field) :: field
    logical :: backward_compatibility, is_constant

    is_aliased=have_option(trim(option_path)//"/aliased")
    if(is_aliased) return
    
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
    ewrite(1,*) "In allocate_and_insert_vector_field, field is: ", trim(lfield_name)

    ! Do we need backward compatibility?
    ! If we need backward compatibility, then no matter how the field
    ! is described in XML, a value space will be allocated, for old-style
    ! code to use.
    ! If we do not need backward compatibility, we can make big savings
    ! on constant fields.
    ! Any fields that require backward compatibility are badly behaved, as they
    ! modify constant fields. *Do not add to this list!* Construct an
    ! appropriate diagnostic algorithm instead (possibly an internal).
    backward_compatibility = .false.

    ! Find out what kind of field we have
    is_prognostic=have_option(trim(path)//"/prognostic")
    is_prescribed=have_option(trim(path)//"/prescribed")
    is_diagnostic=have_option(trim(path)//"/diagnostic")

    is_constant=allocate_vector_field_as_constant(path)

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
       call allocate(field, dim, mesh, name=trim(lfield_name), &
          field_type=FIELD_TYPE_DEFERRED)
    else if(is_constant .and. .not. backward_compatibility) then
         
       ! Allocate as constant field if possible (and we don't need backward compatibility)
       call allocate(field, dim, mesh, name=trim(lfield_name), &
          field_type=FIELD_TYPE_CONSTANT)
       call zero(field)
       
    else
       ! If we have to keep backward compatibility, then
       ! just allocate the value space as normal,
       ! and don't try any funny tricks to save memory.

       ! Allocate field
       call allocate(field, dim, mesh, trim(lfield_name))
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
    call allocate_and_insert_grandchildren(path, state, mesh_name, lfield_name, &
        dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)

    ! Check for adaptivity weights associated with this field:
    adapt_path=trim(path)//"/adaptivity_options"
    if(have_option(trim(adapt_path)//"/absolute_measure")) then
       adapt_path=trim(adapt_path)//"/absolute_measure/vector_field::InterpolationErrorBound"
       call allocate_and_insert_vector_field(adapt_path, state, topology_mesh_name, lfield_name, &
          dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    else if(have_option(trim(adapt_path)//"/relative_measure")) then
       adapt_path=trim(adapt_path)//"/relative_measure/vector_field::InterpolationErrorBound"
       call allocate_and_insert_vector_field(adapt_path, state, topology_mesh_name, lfield_name, &
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

    logical :: backward_compatibility, is_prescribed, is_diagnostic, is_constant, is_aliased
    ! paths for options and child fields
    character(len=OPTION_PATH_LEN) :: path, adapt_path
    character(len=OPTION_PATH_LEN) :: field_name, mesh_name
    type(tensor_field) :: field
    type(mesh_type), pointer:: mesh

    is_aliased=have_option(trim(option_path)//"/aliased")
    if(is_aliased) return

    ! Save option_path
    path=trim(option_path)

    call get_option(trim(path)//"/name", field_name)
    if(present(parent_name)) then
       if(trim(field_name)/="Viscosity") then
          field_name=trim(parent_name)//trim(field_name)
       end if
    end if
    ewrite(1,*) "In allocate_and_insert_tensor_field, field is: ", trim(field_name)

    ! Do we need backward compatibility?
    ! If we need backward compatibility, then no matter how the field
    ! is described in XML, a value space will be allocated, for old-style
    ! code to use.
    ! If we do not need backward compatibility, we can make big savings
    ! on constant fields.
    ! Any fields that require backward compatibility are badly behaved, as they
    ! modify constant fields. *Do not add to this list!* Construct an
    ! appropriate diagnostic algorithm instead (possibly an internal).
    backward_compatibility = any(field_name == (/"ElectricalPotentialDiffusivity      "/))

    ! Find out what kind of field we have
    is_prescribed=have_option(trim(path)//"/prescribed")
    is_diagnostic=have_option(trim(path)//"/diagnostic")
    is_constant=allocate_tensor_field_as_constant(path)

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
       
    else if(is_constant .and. .not. backward_compatibility) then
         
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
       call allocate_and_insert_tensor_field(adapt_path, state, topology_mesh_name, field_name, &
            dont_allocate_prognostic_value_spaces=dont_allocate_prognostic_value_spaces)
    else if(have_option(trim(adapt_path)//"/relative_measure")) then
       adapt_path=trim(adapt_path)//"/relative_measure/tensor_field::InterpolationErrorBound"
       call allocate_and_insert_tensor_field(adapt_path, state, topology_mesh_name, field_name, &
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
    
    ewrite(2,*) "    Inserting children of: ",trim(path)
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
    exclude_interpolated, exclude_nonreprescribed, initial_mesh, time)

    type(state_type), dimension(:), intent(in):: states
    !! don't prescribe the fields with interpolation options
    logical, intent(in), optional :: exclude_interpolated
    !! do not prescribe the fields that have requested not to be represcribed
    logical, intent(in), optional :: exclude_nonreprescribed
    !! indicates whether we're prescribing on the initial mesh, if not (default)
    !! the fields with needs_initial_mesh(field) are left untouched, they have to
    !! be interpolated (somewhere else)
    logical, intent(in), optional:: initial_mesh
    !! current time if not using that in the options tree
    real, intent(in), optional :: time

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
                position, time=time)
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
                position, time=time)
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
                position, time=time)
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
                position, phase_path=trim(phase_path))
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
                position, phase_path=trim(phase_path))
          end if
       end do

    end do
      
    if (.not. present_and_true(save_vtk_cache)) then
       ! flush the cache
       call vtk_cache_finalise()
    end if

  end subroutine initialise_prognostic_fields

  subroutine allocate_and_insert_auxilliary_fields(states, force_prescribed_diagnositc_allocate_old_iterated)
    ! Set up some auxilliary fields to the prognostic fields.
    ! i.e. old and iterated fields depending on options set
    type(state_type), dimension(:), intent(inout):: states
    ! If present and true then force the allocation of Old and 
    ! Iterated fields for prescribed and diagnostic fields
    ! rather than aliasing them.
    logical, intent(in), optional :: force_prescribed_diagnositc_allocate_old_iterated
    
    ! Local variables
    logical :: l_force_prescribed_diagnositc_allocate_old_iterated
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield

    type(scalar_field) :: aux_sfield
    type(vector_field) :: aux_vfield
    type(tensor_field) :: aux_tfield

    integer :: iterations
    logical :: steady_state_global, prognostic, prescribed, diagnostic, gravity

    character(len=FIELD_NAME_LEN) :: field_name
    character(len=OPTION_PATH_LEN) :: state_path, field_path

    integer :: nsfields, nvfields, ntfields, p, f, p2, stat
    real :: current_time

    type(mesh_type), pointer :: x_mesh

    ewrite(1,*) "In allocate_and_insert_auxilliary_fields"
    
    if (present(force_prescribed_diagnositc_allocate_old_iterated)) then
       l_force_prescribed_diagnositc_allocate_old_iterated = force_prescribed_diagnositc_allocate_old_iterated
    else
       l_force_prescribed_diagnositc_allocate_old_iterated = .false.
    end if
    
    call get_option("/timestepping/nonlinear_iterations", iterations, default=1)
    steady_state_global = have_option("/timestepping/steady_state")

    ! old and iterated fields
    do p = 1, size(states)

      ! Get number of scalar fields that are children of this state
      nsfields=scalar_field_count(states(p))

      ! Loop over scalar fields
      sfields_loop: do f=1, nsfields

        sfield => extract_scalar_field(states(p), f)

        ! Save path to field
        field_path=trim(sfield%option_path)

        ! Get field name - this checks if the field has an option_path
        call get_option(trim(field_path)//"/name", field_name, stat)
                
        if((stat==0).and.(.not.aliased(sfield))) then

          prognostic=have_option(trim(sfield%option_path)//"/prognostic")
          prescribed=have_option(trim(sfield%option_path)//"/prescribed")
          diagnostic=have_option(trim(sfield%option_path)//"/diagnostic")

          if( ((prognostic.or.diagnostic) .and. &
               ((steady_state_global.and.steady_state_field(sfield)) .or. &
                (iterations>1) .or. &
                have_option(trim(sfield%option_path) // '/prognostic/spatial_discretisation/discontinuous_galerkin/slope_limiter::FPN') .or. &
                l_force_prescribed_diagnositc_allocate_old_iterated)) .or. &
              (prescribed .and. l_force_prescribed_diagnositc_allocate_old_iterated) ) then

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

          if( ((prognostic.or.diagnostic) .and. &
               ((convergence_field(sfield).and.(iterations>1)) .or. &
                l_force_prescribed_diagnositc_allocate_old_iterated)) .or. &
              (prescribed .and. l_force_prescribed_diagnositc_allocate_old_iterated) ) then

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

      end do sfields_loop

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

          if( ((prognostic.or.diagnostic) .and. &
               ((steady_state_global.and.steady_state_field(vfield)) .or. &
                (iterations>1) .or. &
                l_force_prescribed_diagnositc_allocate_old_iterated)) .or. &
              (prescribed .and. l_force_prescribed_diagnositc_allocate_old_iterated) ) then

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

          if( ((prognostic.or.diagnostic) .and. &
               ((convergence_field(vfield).and.(iterations>1)).or. &
                l_force_prescribed_diagnositc_allocate_old_iterated)) .or. &
              (prescribed .and. l_force_prescribed_diagnositc_allocate_old_iterated)) then

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

      ! Get number of tensor fields that are children of this state
      ntfields=tensor_field_count(states(p))

      ! Loop over tensor fields
      do f=1, ntfields

        tfield => extract_tensor_field(states(p), f)

        ! Save path to field
        field_path=trim(tfield%option_path)

        ! Get field name - this checks if the field has an option_path
        call get_option(trim(field_path)//"/name", field_name, stat)
                
        if((stat==0).and.(.not.aliased(tfield))) then

          prognostic=have_option(trim(tfield%option_path)//"/prognostic")
          prescribed=have_option(trim(tfield%option_path)//"/prescribed")
          diagnostic=have_option(trim(tfield%option_path)//"/diagnostic")

          if( ((prognostic.or.diagnostic) .and. &
                l_force_prescribed_diagnositc_allocate_old_iterated) .or. &
              (prescribed .and. l_force_prescribed_diagnositc_allocate_old_iterated) ) then

            call allocate(aux_tfield, tfield%mesh, "Old"//trim(tfield%name))
            call zero(aux_tfield)
            call insert(states(p), aux_tfield, trim(aux_tfield%name))
            call deallocate(aux_tfield)

          else

            aux_tfield = extract_tensor_field(states(p), trim(tfield%name))
            aux_tfield%name = "Old"//trim(tfield%name)
            aux_tfield%option_path=""  ! blank the option path so that it 
                                       ! doesn't get picked up in the next 
                                       ! aliased field loop
            aux_tfield%aliased=.true.
            call insert(states(p), aux_tfield, trim(aux_tfield%name))

          end if

          if( ((prognostic.or.diagnostic) .and. &
                l_force_prescribed_diagnositc_allocate_old_iterated) .or. &
              (prescribed .and. l_force_prescribed_diagnositc_allocate_old_iterated) ) then

            call allocate(aux_tfield, tfield%mesh, "Iterated"//trim(tfield%name))
            call zero(aux_tfield)
            call insert(states(p), aux_tfield, trim(aux_tfield%name))
            call deallocate(aux_tfield)

          else

            aux_tfield = extract_tensor_field(states(p), trim(tfield%name))
            aux_tfield%name = "Iterated"//trim(tfield%name)
            aux_tfield%option_path=""  ! blank the option path so that it 
                                       ! doesn't get picked up in the next 
                                       ! aliased field loop
            aux_tfield%aliased=.true.
            call insert(states(p), aux_tfield, trim(aux_tfield%name))

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

          if(prognostic.or.prescribed.or.diagnostic) then

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

          if(prognostic.or.prescribed.or.diagnostic) then

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

      ! Get number of tensor fields that are children of this state
      ntfields=tensor_field_count(states(p))

      ! Loop over tensor fields
      do f=1, ntfields

        tfield => extract_tensor_field(states(p), f)

        ! Save path to field
        field_path=trim(tfield%option_path)

        ! Get field name - this checks if the field has an option_path
        ! but if it's aliased the name that it gets from the option path will be of the field it's aliased to!
        call get_option(trim(field_path)//"/name", field_name, stat)

        if((stat==0).and.aliased(tfield).and.(tfield%option_path(:15)=="/material_phase")) then

          prognostic=have_option(trim(tfield%option_path)//"/prognostic")
          prescribed=have_option(trim(tfield%option_path)//"/prescribed")
          diagnostic=have_option(trim(tfield%option_path)//"/diagnostic")

          if(prognostic.or.prescribed.or.diagnostic) then

            do p2 = 1, size(states)
              write(state_path, '(a,i0,a)') "/material_phase[",p2-1,"]"
              if(starts_with(trim(field_path), trim(state_path))) exit
            end do

            if(p2==size(states)+1) then
              FLAbort("tensor_field aliased but could not find to which material_phase")
            end if

            aux_tfield=extract_tensor_field(states(p2), "Old"//trim(field_name))
            aux_tfield%name = "Old"//trim(tfield%name)
            aux_tfield%aliased = .true.
            aux_tfield%option_path = ""  ! blank the option path for consistency
            call insert(states(p), aux_tfield, trim(aux_tfield%name))

            aux_tfield=extract_tensor_field(states(p2), "Iterated"//trim(field_name))
            aux_tfield%name = "Iterated"//trim(tfield%name)
            aux_tfield%aliased = .true.
            aux_tfield%option_path = ""  ! blank the option path for consistency
            call insert(states(p), aux_tfield, trim(aux_tfield%name))

          end if

        end if

      end do

    end do      
    
    ! for mesh movement we need a "OriginalCoordinate", 
    ! "OldCoordinate" and "IteratedCoordinate" fields
    ! inserted in each state similar to "Coordinate"
    if (have_option('/mesh_adaptivity/mesh_movement')) then 
       vfield => extract_vector_field(states(1), name="Coordinate")
       
      ! first original coordinate field:
      call allocate(aux_vfield, vfield%dim, vfield%mesh, &
        name="Original"//trim(vfield%name))
      call set(aux_vfield, vfield)
      aux_vfield%option_path=""
      ! insert into states(1) and alias it to all other states.
      call insert(states, aux_vfield, trim(aux_vfield%name))
      call deallocate(aux_vfield)
       
      ! exactly the same for old coordinate field:
      call allocate(aux_vfield, vfield%dim, vfield%mesh, &
        name="Old"//trim(vfield%name))
      call set(aux_vfield, vfield)
      aux_vfield%option_path=""
      ! insert into states(1) and alias it to all other states.
      call insert(states, aux_vfield, trim(aux_vfield%name))
      call deallocate(aux_vfield)
        
      ! and again for the iterated coordinate field (the most up to date one):
      call allocate(aux_vfield, vfield%dim, vfield%mesh, &
        name="Iterated"//trim(vfield%name))
      call set(aux_vfield, vfield)
      aux_vfield%option_path=""
      ! insert into states(1) and alias it to all other states.
      call insert(states, aux_vfield, trim(aux_vfield%name))
      call deallocate(aux_vfield)
       
    else
    
      aux_vfield=extract_vector_field(states(1), name="Coordinate")
      aux_vfield%name = "Original"//trim(aux_vfield%name)
      aux_vfield%aliased = .true.
      aux_vfield%option_path = ""
      ! insert into states(1) and alias it to all other states.
      call insert(states, aux_vfield, trim(aux_vfield%name))
        
      aux_vfield=extract_vector_field(states(1), name="Coordinate")
      aux_vfield%name = "Old"//trim(aux_vfield%name)
      aux_vfield%aliased = .true.
      aux_vfield%option_path = ""
      ! insert into states(1) and alias it to all other states.
      call insert(states, aux_vfield, trim(aux_vfield%name))
        
      aux_vfield=extract_vector_field(states(1), name="Coordinate")
      aux_vfield%name = "Iterated"//trim(aux_vfield%name)
      aux_vfield%aliased = .true.
      aux_vfield%option_path = ""
      ! insert into states(1) and alias it to all other states.
      call insert(states, aux_vfield, trim(aux_vfield%name))
        
    end if

    x_mesh => extract_mesh(states(1), "CoordinateMesh")
    ! need a GridVelocity even if we're not moving the mesh
    if (.not.have_option('/mesh_adaptivity/mesh_movement')) then
      call allocate(aux_vfield, mesh_dim(x_mesh), x_mesh, "GridVelocity", field_type = FIELD_TYPE_CONSTANT)
      call zero(aux_vfield)
      aux_vfield%option_path = ""
      call insert(states, aux_vfield, trim(aux_vfield%name))
      call deallocate(aux_vfield)
    end if

    ! Disgusting and vomitous hack to ensure that time is output in
    ! vtu files.
    call allocate(aux_sfield, x_mesh, "Time", field_type=FIELD_TYPE_CONSTANT)
    call get_option("/timestepping/current_time", current_time)
    call set(aux_sfield, current_time)
    aux_sfield%option_path = ""
    call insert(states, aux_sfield, trim(aux_sfield%name))
    call deallocate(aux_sfield)
    
    ! Porous media fields - insert alias of fields in first state into all others
    have_porous_media: if (have_option('/porous_media')) then
       
       ! alias the OldPorosity field
       aux_sfield=extract_scalar_field(states(1), 'OldPorosity')
       aux_sfield%aliased = .true.
       aux_sfield%option_path = ""
       do p = 1,size(states)-1
          call insert(states(p+1), aux_sfield, 'OldPorosity')
       end do
       
       ! alias the OldAbsolutePermeability field which may be 
       ! either scalar or vector (if present)
       
       aux_sfield=extract_scalar_field(states(1), 'OldAbsolutePermeability', stat = stat)
       if (stat == 0) then       
          aux_sfield%aliased = .true.
          aux_sfield%option_path = ""
          do p = 1,size(states)-1
             call insert(states(p+1), aux_sfield, 'OldAbsolutePermeability')
          end do       
       end if
       
       aux_vfield=extract_vector_field(states(1), 'OldAbsolutePermeability', stat = stat)
       if (stat == 0) then       
          aux_vfield%aliased = .true.
          aux_vfield%option_path = ""
          do p = 1,size(states)-1
             call insert(states(p+1), aux_vfield, 'OldAbsolutePermeability')
          end do       
       end if

       ! alias the IteratedPorosity field
       aux_sfield=extract_scalar_field(states(1), 'IteratedPorosity')
       aux_sfield%aliased = .true.
       aux_sfield%option_path = ""
       do p = 1,size(states)-1
          call insert(states(p+1), aux_sfield, 'IteratedPorosity')
       end do
       
       ! alias the IteratedAbsolutePermeability field which may be 
       ! either scalar or vector (if present)
       
       aux_sfield=extract_scalar_field(states(1), 'IteratedAbsolutePermeability', stat = stat)
       if (stat == 0) then       
          aux_sfield%aliased = .true.
          aux_sfield%option_path = ""
          do p = 1,size(states)-1
             call insert(states(p+1), aux_sfield, 'IteratedAbsolutePermeability')
          end do       
       end if
       
       aux_vfield=extract_vector_field(states(1), 'IteratedAbsolutePermeability', stat = stat)
       if (stat == 0) then       
          aux_vfield%aliased = .true.
          aux_vfield%option_path = ""
          do p = 1,size(states)-1
             call insert(states(p+1), aux_vfield, 'IteratedAbsolutePermeability')
          end do       
       end if
    
    end if have_porous_media

    ! Porous media dual fields - insert alias of fields in first state into all others
    have_porous_media_dual: if (have_option('/porous_media_dual')) then
       
       ! alias the OldPorosityDual field
       aux_sfield=extract_scalar_field(states(1), 'OldPorosityDual')
       aux_sfield%aliased = .true.
       aux_sfield%option_path = ""
       do p = 1,size(states)-1
          call insert(states(p+1), aux_sfield, 'OldPorosityDual')
       end do
       
       ! alias the OldAbsolutePermeabilityDual field if present
       aux_sfield=extract_scalar_field(states(1), 'OldAbsolutePermeabilityDual', stat=stat)
       if (stat == 0) then
          aux_sfield%aliased = .true.
          aux_sfield%option_path = ""
          do p = 1,size(states)-1
             call insert(states(p+1), aux_sfield, 'OldAbsolutePermeabilityDual')
          end do       
       end if
       
       ! alias the IteratedPorosityDual field
       aux_sfield=extract_scalar_field(states(1), 'IteratedPorosityDual')
       aux_sfield%aliased = .true.
       aux_sfield%option_path = ""
       do p = 1,size(states)-1
          call insert(states(p+1), aux_sfield, 'IteratedPorosityDual')
       end do
       
       ! alias the IteratedAbsolutePermeabilityDual field if present
       aux_sfield=extract_scalar_field(states(1), 'IteratedAbsolutePermeabilityDual', stat=stat)
       if (stat == 0) then
          aux_sfield%aliased = .true.
          aux_sfield%option_path = ""
          do p = 1,size(states)-1
             call insert(states(p+1), aux_sfield, 'IteratedAbsolutePermeabilityDual')
          end do       
       end if
       
    end if have_porous_media_dual
    
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
             FLExit("No mesh for field "//trim(field_path))
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
    
  subroutine create_empty_halo(position)
    !!< Auxilary subroutine that creates node and element halos for position with no sends or receives
    type(vector_field), intent(inout):: position
    
    integer:: nprocs, j
    
    nprocs = getnprocs()
    allocate(position%mesh%halos(2))
    allocate(position%mesh%element_halos(2))
    do j=1,2
      ! Nodal halo
      call allocate(position%mesh%halos(j), nprocs = nprocs, nreceives = spread(0, 1, nprocs), &
                    nsends = spread(0, 1, nprocs), &
                    data_type=HALO_TYPE_CG_NODE, ordering_scheme=HALO_ORDER_TRAILING_RECEIVES, &
                    nowned_nodes = node_count(position),&
                    name="EmptyHalo")
      assert(trailing_receives_consistent(position%mesh%halos(j)))
      call create_global_to_universal_numbering(position%mesh%halos(j))
      call create_ownership(position%mesh%halos(j))

      ! Element halo
      call allocate(position%mesh%element_halos(j), nprocs = nprocs, nreceives = spread(0, 1, nprocs), &
                    nsends = spread(0, 1, nprocs), &
                    data_type=HALO_TYPE_ELEMENT, ordering_scheme=HALO_ORDER_TRAILING_RECEIVES, &
                    nowned_nodes = ele_count(position), &
                    name="EmptyHalo")
      assert(trailing_receives_consistent(position%mesh%element_halos(j)))
      call create_global_to_universal_numbering(position%mesh%element_halos(j))
      call create_ownership(position%mesh%element_halos(j))
    end do
    
  end subroutine create_empty_halo


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
    mesh => extract_mesh(state, trim(topology_mesh_name))

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

  subroutine compute_domain_statistics(states)
    type(state_type), dimension(:), intent(in) :: states
    integer :: dim
    type(vector_field), pointer :: positions
    integer :: ele
    real :: vol
    type(scalar_field) :: temp_s_field

    positions => extract_vector_field(states(1), "Coordinate")
    if (allocated(domain_bbox)) then
      deallocate(domain_bbox)
    end if
    allocate(domain_bbox(positions%dim, 2))
    domain_bbox = 0.0

    do dim=1,positions%dim
      domain_bbox(dim, 1) = minval(positions%val(dim,:))
      domain_bbox(dim, 2) = maxval(positions%val(dim,:))
      ewrite(2,*) "domain_bbox - dim, range =", dim, domain_bbox(dim,:)
    end do

    vol = 0.0
    do ele=1,ele_count(positions)
      vol = vol + element_volume(positions, ele)
    end do

    domain_volume = vol
    ewrite(2,*) "domain_volume =", domain_volume

    !If on-the-sphere, calculate the radius of the sphere.
    if (have_option("/geometry/spherical_earth/")) then
      temp_s_field = magnitude(positions)
      surface_radius = maxval(temp_s_field)
      call allmax(surface_radius)
      ! Need to deallocate the magnitude field create, or we get a leak
      call deallocate(temp_s_field)
    end if


  end subroutine compute_domain_statistics
  
  subroutine populate_state_module_check_options

    character(len=OPTION_PATH_LEN) :: problem_type

    ! Check mesh options
    call check_mesh_options

    call check_geometry_options

    ! check problem specific options:
    call get_option("/problem_type", problem_type)
    select case (problem_type)
    case ("fluids")
    case ("oceans")
       call check_ocean_options
    case ("large_scale_ocean_options")
       call check_large_scale_ocean_options
    case ("multimaterial")
       call check_multimaterial_options
    case ("stokes")
       call check_stokes_options
    case ("foams")
       call check_foams_options
    case ("multiphase")
       call check_multiphase_options
    case default
       ewrite(0,*) "Problem type:", trim(problem_type)
       FLAbort("Error unknown problem_type")
    end select
    ewrite(2,*) 'Done with problem type choice'

  end subroutine populate_state_module_check_options

  subroutine check_geometry_options

    logical :: on_sphere
    integer :: i, nstates

    on_sphere = have_option("/geometry/spherical_earth")

    if (on_sphere) then
      nstates=option_count("/material_phase")

      state_loop: do i=0, nstates-1
        if (have_option("/material_phase[" // int2str(i) // "]/vector_field::Velocity/prognostic")) then
          if (.not. (have_option("/material_phase[" // int2str(i) // "]/vector_field::Velocity/prognostic" // &
                                 "/spatial_discretisation/continuous_galerkin/buoyancy" // &
                                 "/radial_gravity_direction_at_gauss_points") .or. &
                     have_option("/material_phase[" // int2str(i) // "]/vector_field::Velocity/prognostic" // &
                                 "/spatial_discretisation/discontinuous_galerkin/buoyancy" // &
                                 "/radial_gravity_direction_at_gauss_points"))) then
            ewrite(0,*) "WARNING: the /geometry/spherical_earth option no long automatically makes the buoyancy radial."
            ewrite(0,*) "To recreate the previous behaviour it is now necessary to turn on the "
            ewrite(0,*) "buoyancy/radial_gravity_direction_at_gauss_points underneath the Velocity spatial_discretisation."
          end if
        end if
      end do state_loop
    end if

  end subroutine check_geometry_options

  subroutine check_mesh_options

    character(len=OPTION_PATH_LEN) :: path
    character(len=OPTION_PATH_LEN) :: field_name, mesh_name, from_mesh_name, phase_name
    integer :: i, j ! counters
    integer :: nstates ! number of states
    integer :: nfields ! number of fields
    integer :: nmeshes ! number of meshes
    integer :: n_external_meshes ! number of meshes from file
    integer :: n_external_meshes_excluded_from_mesh_adaptivity
    integer :: periodic_mesh_count ! number of meshes with periodic_boundary_conition options
    ! logicals to find out if we have certain options
    logical :: is_aliased

    ! Get number of meshes
    nmeshes=option_count("/geometry/mesh")

    ewrite(2,*) "Checking mesh options."
    ewrite(2,*) "There are", nmeshes, "meshes."

    n_external_meshes=0
    n_external_meshes_excluded_from_mesh_adaptivity=0

    mesh_loop1: do i=0, nmeshes-1

       ! Save mesh path
       path="/geometry/mesh["//int2str(i)//"]"

       if(have_option(trim(path)//"/from_file")) then

          n_external_meshes=n_external_meshes+1
          
          if (have_option(trim(path)//"/exclude_from_mesh_adaptivity")) then
            n_external_meshes_excluded_from_mesh_adaptivity=n_external_meshes_excluded_from_mesh_adaptivity+1
          end if

       else if (.not. have_option(trim(path)//"/from_mesh")) then

          call get_option(trim(path)//"/name", mesh_name)
          ewrite(-1,*) "In options for /geometry/mesh ("//trim(mesh_name)//"):"
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
    if(n_external_meshes-n_external_meshes_excluded_from_mesh_adaptivity>1) then
      ewrite(-1,*) "With multiple external (from_file) meshes"
      FLExit("Only one external mesh may leave out the exclude_from_mesh_adaptivity option.")
    end if

    ! Check that dimension of mesh is the same as the dimension defined in the options file
    ! ...that's not so easy: let's just do this in insert_external_mesh() with a nice FLEXit

    periodic_mesh_count = 0
    ! Check that the meshes required to make other meshes are present.
    mesh_loop2: do i=0, nmeshes-1

       ! Save mesh path
       path="/geometry/mesh["//int2str(i)//"]"
       call get_option(trim(path)//"/name", mesh_name)

       if (have_option(trim(path)//"/from_mesh")) then

          call get_option(trim(path)//"/from_mesh/mesh[0]/name", from_mesh_name)
          if (.not. have_option("/geometry/mesh::"//trim(from_mesh_name))) then

             ewrite(-1,*) "Unknown mesh: ", trim(from_mesh_name)
             ewrite(-1,*) "Specified as source (from_mesh) for ", trim(mesh_name)
             FLExit("Error in /geometry/mesh: unknown mesh.")

          end if
          
          if (have_option("/geometry/mesh::"//trim(from_mesh_name)//&
             "/exclude_from_mesh_adaptivity") .and. .not. &
             have_option(trim(path)//"/exclude_from_mesh_adaptivity")) then
             ! if the from_mesh is excluded, the mesh itself also needs to be
             ewrite(-1,*) "In derivation of mesh ", trim(mesh_name), " from ", trim(from_mesh_name)
             ewrite(-1,*) "A mesh derived from a mesh with exclude_from_mesh_adaptivity needs to have this options as well."
             FLExit("Missing exclude_from_mesh_adaptivity option")
          end if
          
          if (have_option(trim(path)//"/from_mesh/extrude") .and. ( &
             
             have_option(trim(path)//"/from_mesh/mesh_shape") .or. &
             have_option(trim(path)//"/from_mesh/mesh_continuity") .or. &
             have_option(trim(path)//"/from_mesh/periodic_boundary_conditions") &
             ) ) then
             
             ewrite(-1,*) "In derivation of mesh ", trim(mesh_name), " from ", trim(from_mesh_name)
             ewrite(-1,*) "When extruding a mesh, you cannot at the same time"
             ewrite(-1,*) "change its shape, continuity or add periodic bcs."
             ewrite(-1,*) "Need to do this in seperate step (derivation)."
             FLExit("Error in /geometry/mesh with extrude option")
             
          end if

          if (have_option(trim(path)//"/from_mesh/periodic_boundary_conditions")) then
            
             ! can't combine with anything else
             if ( &
             
                have_option(trim(path)//"/from_mesh/mesh_shape") .or. &
                have_option(trim(path)//"/from_mesh/mesh_continuity") .or. &
                have_option(trim(path)//"/from_mesh/extrude") &
             ) then
             
                ewrite(-1,*) "In derivation of mesh ", trim(mesh_name), " from ", trim(from_mesh_name)
                ewrite(-1,*) "When adding or removing periodicity to a mesh, you cannot at the same time"
                ewrite(-1,*) "change its shape, continuity or extrude a mesh."
                ewrite(-1,*) "Need to do this in seperate step (derivation)."
                FLExit("Error in /geometry/mesh with extrude option")
                
             end if
             
             if (have_option(trim(path)//"/from_mesh/periodic_boundary_conditions/remove_periodicity")) then
               
                ! check to see the from_mesh is not non-periodic is done above
                
                ! check that all periodic bcs have remove_periodicity
                if (option_count(trim(path)//"/from_mesh/periodic_boundary_conditions")/= &
                        option_count(trim(path)//"/from_mesh/periodic_boundary_conditions/remove_periodicity")) then
                   ewrite(-1,*) "In derivation of mesh ", trim(mesh_name), " from ", trim(from_mesh_name)
                   FLExit("All or none of the periodic_boundary_conditions need to have the option remove_periodicity.")
                end if
                
             else
             
                ! really periodic
                
                if (mesh_name=="CoordinateMesh") then
                  ewrite(-1,*) "In derivation of mesh ", trim(mesh_name), " from ", trim(from_mesh_name)
                  FLExit("CoordinateMesh may not be made periodic.")
                end if
                
                periodic_mesh_count=periodic_mesh_count+1
                
                if (periodic_mesh_count>1) then
                  ewrite(-1,*) "In the derivation of periodic meshes, all periodic boundary conditions"
                  ewrite(-1,*) "have to be applied at once. Thus only one mesh may have periodic_boundary_conditions"
                  ewrite(-1,*) "specified under /geometry/mesh::PeriodicMesh/from_mesh and all other periodic meshes"
                  ewrite(-1,*) "should be derived from this mesh."
                  FLExit("More than one mesh with periodic_boundary_conditions")
                end if
                
                if (.not. have_option("/geometry/mesh::"//trim(from_mesh_name)// &
                  "/from_file")) then
                  ewrite(-1,*) "In derivation of mesh ", trim(mesh_name), " from ", trim(from_mesh_name)
                  ewrite(-1,*) "In the derivation of periodic meshes, the first periodic mesh,"
                  ewrite(-1,*) "which has the periodic_boundary_conditions specified, must be derived"
                  ewrite(-1,*) "directly from the external (from_file) mesh."
                  FLExit("Periodic mesh not from from_file mesh")
                end if
             end if
             
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
             call get_option(trim(complete_field_path(path))//"/mesh[0]/name", mesh_name)

             if (.not. have_option("/geometry/mesh::"//trim(mesh_name))) then

                ewrite(-1,*) "Unknown mesh: ", trim(mesh_name)
                ewrite(-1,*) "Specified as mesh for scalar_field ", trim(field_name)
                ewrite(-1,*) "In material_phase ", trim(phase_name)
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
             call get_option(trim(complete_field_path(path))//"/mesh[0]/name", mesh_name)

             if (.not. have_option("/geometry/mesh::"//trim(mesh_name))) then

                ewrite(-1,*) "Unknown mesh: ", trim(mesh_name)
                ewrite(-1,*) "Specified as mesh for vector_field ", trim(field_name)
                ewrite(-1,*) "In material_phase ", trim(phase_name)
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
             call get_option(trim(complete_field_path(path))//"/mesh[0]/name", mesh_name)

             if (.not. have_option("/geometry/mesh::"//trim(mesh_name))) then

                ewrite(-1,*) "Unknown mesh: ", trim(mesh_name)
                ewrite(-1,*) "Specified as mesh for tensor_field ", trim(field_name)
                ewrite(-1,*) "In material_phase ", trim(phase_name)
                FLExit("Error: unknown mesh.")

             end if

          end if

       end do tensor_field_loop

    end do state_loop

  end subroutine check_mesh_options
  
  subroutine check_ocean_options

    character(len=OPTION_PATH_LEN) str, velocity_path, pressure_path, tmpstring
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

  subroutine check_large_scale_ocean_options
  
    character(len=OPTION_PATH_LEN) str, velocity_path, pressure_path, tmpstring, temperature_path, salinity_path,continuity2, continuity1, velmesh, pressuremesh, preconditioner
    logical on_sphere, constant_gravity
    integer iterations, poly
    if (option_count('/material_phase')/=1) then
       FLExit("The checks for problem_type oceans only work for single phase.")
    endif
    
    ! Velocity options checks   
    velocity_path="/material_phase[0]/vector_field::Velocity/prognostic"
    if (have_option(trim(velocity_path))) then
       str=trim(velocity_path)//'/spatial_discretisation/continuous_galerkin'
       if (have_option(trim(str))) then
          FLExit("For large scale ocean problems you need discontinuous galerkin velocity.")
       end if
       if (.not. have_option(trim(velocity_path)//'/equation::Boussinesq')) then
          FLExit("Wrong Velocity equation type - should be Boussinesq")
       end if  
       if(.not.(have_option(trim(velocity_path)//"/spatial_discretisation/discontinuous_galerkin/advection_scheme/upwind")).and. &
         (.not.(have_option("timestepping/steady_state")))) then
          ewrite(0,*)("WARNING: You should probably have advection_scheme/upwind under velocity")
       end if  
       if(.not.(have_option(trim(velocity_path)//"/spatial_discretisation/discontinuous_galerkin/advection_scheme/integrate_advection_by_parts/twice"))) then 
          FLExit("Should have Velocity/spatial_discretisation/advection_scheme/integrate_advection_by_parts/twice")
       end if
       if(have_option(trim(velocity_path)//"/spatial_discretisation/discontinuous_galerkin/mass_terms/lump_mass_matrix")) then
         FLExit("Should not lump mass matrix in large-scale ocean simulations")
       end if
       if (.not.have_option(trim(velocity_path)//"/spatial_discretisation/discontinuous_galerkin/viscosity_scheme/bassi_rebay").and. .not.have_option(trim(velocity_path)//"/spatial_discretisation/discontinuous_galerkin/viscosity_scheme/compact_discontinuous_galerkin")) then
         FLExit("Should have Bassi Rebay or compact discontinuous galerkin Viscosity scheme (under Velocity)")
       end if
       if (have_option(trim(velocity_path)//"/spatial_discretisation/discontinuous_galerkin/viscosity_scheme/compact_discontinuous_galerkin") ) then
          call get_option(trim(velocity_path)//"/solver/preconditioner/name",preconditioner)
          if (preconditioner .ne. "sor") then
            FLExit("You need sor preconditioner for velocity with compact discontinuous galerkin viscosity.")
          end if
       end if
       if (.not.have_option(trim(velocity_path)//"/temporal_discretisation/discontinuous_galerkin/maximum_courant_number_per_subcycle")) then
        ewrite(0,*)  ("WARNING: You may wish to switch on velocity/prognostic/temporal_discretisation/discontinuous_galerkin/maximum_courant_number_per_subcycle ")
       end if
       if (.not.have_option(trim(velocity_path)//"/spatial_discretisation/discontinuous_galerkin/advection_scheme/project_velocity_to_continuous")) then
         FLExit("You need to switch on velocity/prognostic/spatial_discretisation/discontinuous_galerkin/advection_scheme/project_velocity_to_continuous ")
       end if
  
  end if
  
!Timestepping options
  if(.not.(have_option("/timestepping/nonlinear_iterations"))) then
     FLExit("You should turn on timestepping/nonlinear_iterations and set to a number greater than 1")
  end if
  if((have_option("/timestepping/nonlinear_iterations"))) then
     call get_option(("/timestepping/nonlinear_iterations"), iterations)
       if(iterations .lt. 2 ) then
         FLExit("timestepping/nonlinear_iterations should be set to a number greater than 1")
       end if
  end if

! Subtract out hydrostatic level option
  if(.not.(have_option("material_phase/equation_of_state/fluids/linear/subtract_out_hydrostatic_level"))) then
     FLExit("You should switch on material_phase/equation_of_state/subtract_out_hydrostatic_level")
  end if

! Geometry ocean boundaries
  if(.not.(have_option("/geometry/ocean_boundaries"))) then
     FLExit("You need to switch on geometry/ocean_boundaries")
  end if

!Pressure options checks
  pressure_path="/material_phase[0]/scalar_field::Pressure/prognostic"
    if(have_option("/material_phase[0]/scalar_field::Pressure/prognostic")) then
       if (.not.have_option(trim(pressure_path)//"/scheme/use_projection_method")) then
          FLExit("For ocean problems you should use the projection method under scheme for pressure")
       end if
       call get_option(trim(pressure_path)//"/scheme/poisson_pressure_solution", tmpstring)
       select case (tmpstring)
       case ("never")
          ewrite(0,*) ("WARNING: Poisson pressure solution is set to never.")
       case ("only first timestep")
         ewrite(0,*)("WARNING: Poisson pressure solution is set to only first time step")
       end select
       if (.not.have_option(trim(pressure_path)//"/spatial_discretisation/continuous_galerkin")) then
          FLExit("For ocean problems you should use continuous galerkin pressure")
       end if
       if (.not.have_option(trim(pressure_path)//"/solver/preconditioner/vertical_lumping")) then
          ewrite(0,*)("WARNING: Vertical lumping not used during pressure solve.  Consider switching on pressure/vertical_lumping.")
       end if
       if (.not.have_option(trim(pressure_path)//"/spatial_discretisation/continuous_galerkin/remove_stabilisation_term")) then
          FLExit("Use remove stabilisation term under pressure")
       end if
       if (.not.have_option(trim(pressure_path)//"/spatial_discretisation/continuous_galerkin/integrate_continuity_by_parts")) then
          FLExit("Use integrate continuity by parts under pressure")
       end if
    end if

    ! Salinity options checks
    salinity_path="/material_phase[0]/scalar_field::Salinity/prognostic"
    if(have_option("/material_phase[0]/scalar_field::Salinity")&
         .and.(.not.(have_option("/material_phase[0]/equation_of_state/fluids/linear/salinity_dependency")))) then
       ewrite(0,*) "WARNING: You have a salinity field but it will not affect the density of the fluid."
    end if

   ! Temperature options checks
   temperature_path="/material_phase[0]/scalar_field::Temperature/prognostic"
    if(have_option("/material_phase[0]/scalar_field::Temperature")&
         .and.(.not.(have_option("/material_phase[0]/equation_of_state/fluids/linear/temperature_dependency") ))) then
       ewrite(0,*) "WARNING: You have a temperature field but it will not affect the density of the fluid."
    end if
   
    ! Check that the gravity field is not constant for spherical problems
    on_sphere=have_option("/geometry/spherical_earth")
    constant_gravity=have_option("/physical_parameters/gravity/vector_field::GravityDirection/prescribed/value[0]/constant")
    if(on_sphere .and. constant_gravity) then
       FLExit("GravityDirection set incorrectly for spherical geometry.")
    end if

! Check velocity mesh continuity
  call get_option("/material_phase[0]/vector_field::Velocity/prognostic/mesh/name",velmesh)
  call get_option("/geometry/mesh::"//trim(velmesh)//"/from_mesh/mesh_continuity",continuity2)
    
  if (trim(continuity2).ne."discontinuous") then
    FLExit("The velocity mesh is not discontinuous")
  end if

  ! Check pressure mesh continuity 
  call get_option("/material_phase[0]/scalar_field::Pressure/prognostic/mesh/name",pressuremesh)
  if (have_option("/geometry/mesh::"//trim(pressuremesh)//"/from_mesh/mesh_continuity"))then
    call get_option("/geometry/mesh::"//trim(pressuremesh)//"/from_mesh/mesh_continuity",continuity1)
    if (trim(continuity1).ne."continuous")then
      FLExit ("Pressure mesh is not continuous")
    end if
  end if
  ! Check pressure mesh polynomial order
  if (.not.have_option("/geometry/mesh::"//trim(pressuremesh)//"/from_mesh/mesh_shape"))then
    ewrite (0,*)"WARNING: You should have the pressure mesh shape set to polynomial order 2"
  end if
  if (have_option("/geometry/mesh::"//trim(pressuremesh)//"/from_mesh/mesh_shape/polynomial_degree"))then
    call get_option("/geometry/mesh::"//trim(pressuremesh)//"/from_mesh/mesh_shape/polynomial_degree",poly)
    if (poly.ne.2) then
      ewrite (0,*)"WARNING: You should have the pressure mesh shape set to polynomial order 2"
    end if 
  end if

    !Check for viscosity field
  if (.not.have_option("/material_phase[0]/vector_field::Velocity/prognostic/tensor_field::Viscosity"))then
    ewrite(0,*)"WARNING: You have no viscosity field"
  end if
  ! Check for absorption term
if (.not.have_option("/material_phase[0]/vector_field::Velocity/prognostic/vector_field::Absorption"))then
    ewrite(0,*)"WARNING: you may wish to add an absorption term under velocity"
  end if
  !Check for temperature diffusivity
  if (have_option("/material_phase[0]/scalar_field::Temperature/prognostic")) then
    if (.not. have_option("/material_phase[0]/scalar_field::Temperature/prognostic/tensor_field::Diffusivity")) then
    ewrite(0,*)"WARNING: you have a prognostic temperature field but no diffusivity"
    end if
  end if
  !Check for salinity diffusivity
  if (have_option("/material_phase[0]/scalar_field::Salinity/prognostic")) then
    if (.not. have_option("/material_phase[0]/scalar_field::Salinity/prognostic/tensor_field::Diffusivity")) then
    ewrite(0,*)"WARNING: you have a prognostic salinity field but no diffusivity"
    end if
  end if

  end subroutine check_large_scale_ocean_options

  subroutine check_multimaterial_options

    integer :: neos, nmat, i
    logical :: have_vfrac, have_dens
    
    integer :: diagnosticvolumefraction_count, density_count, &
               viscosity_count, surfacetension_count

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
    
    diagnosticvolumefraction_count = option_count(&
                &'/material_phase/scalar_field::MaterialVolumeFraction/diagnostic')
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

    viscosity_count = option_count('/material_phase/&
                &vector_field::Velocity/prognostic/&
                &tensor_field::Viscosity/diagnostic')
    if(viscosity_count>1) then
      ewrite(-1,*) viscosity_count, 'diagnostic bulk Viscosities.'
      FLExit("Only 1 diagnostic bulk Viscosity is allowed")
    end if

    surfacetension_count = option_count('/material_phase/&
                &vector_field::Velocity/prognostic&
                &/tensor_field::SurfaceTension/diagnostic')
    if(surfacetension_count>1) then
      ewrite(-1,*) surfacetension_count, 'diagnostic surface tensions.'
      FLExit("Only 1 diagnostic surface tension is allows")
    end if

  end subroutine check_multimaterial_options

  subroutine check_stokes_options

    ! Check options for Stokes flow simulations.

    integer :: i, nmat
    character(len=OPTION_PATH_LEN) :: velocity_path, pressure_path, schur_path
    character(len=FIELD_NAME_LEN)  :: schur_preconditioner, inner_matrix, pc_type
    logical :: exclude_mass, exclude_advection
    real :: theta

    nmat = option_count("/material_phase")

    do i = 0, nmat-1
       velocity_path="/material_phase["//int2str(i)//"]/vector_field::Velocity/prognostic"

       if (have_option(trim(velocity_path))) then
         
          ! Check that mass and advective terms are excluded:
          exclude_mass = have_option(trim(velocity_path)//&
               "/spatial_discretisation/continuous_galerkin/mass_terms"//&
               &"/exclude_mass_terms").or.&
               have_option(trim(velocity_path)//&
               "/spatial_discretisation/discontinuous_galerkin/mass_terms"//&
               &"/exclude_mass_terms")
          
          exclude_advection = have_option(trim(velocity_path)//&
               "/spatial_discretisation/continuous_galerkin/advection_terms"//&
               &"/exclude_advection_terms").or.&
               have_option(trim(velocity_path)//&
               "/spatial_discretisation/discontinuous_galerkin/advection_scheme/none") 
          
          if(.not.(exclude_mass) .OR. .not.(exclude_advection)) then
             FLExit("For Stokes problems you need to exclude the mass and advection terms.")
          end if

          ! Check that theta = 1 (we must be implicit as we have no time term!)
          call get_option(trim(velocity_path)//'/temporal_discretisation/theta/', theta)
          if(theta /= 1.) then
             FLExit("For Stokes problems, theta (under velocity) must = 1")
          end if

       end if

       pressure_path="/material_phase["//int2str(i)//"]/scalar_field::Pressure/prognostic"

       if (have_option(trim(pressure_path))) then  

          schur_path = "/material_phase["//int2str(i)//"]/scalar_field::Pressure/prognostic/"//&
               &"scheme/use_projection_method/full_schur_complement"

          if(have_option(trim(schur_path))) then
             
             call get_option(trim(schur_path)//"/preconditioner_matrix[0]/name", schur_preconditioner)
             
             select case(schur_preconditioner)
             case("ScaledPressureMassMatrix")
                ! Check pressure_mass_matrix preconditioner is compatible with viscosity tensor:
                if(have_option(trim(velocity_path)//&
                     &"/tensor_field::Viscosity/prescribed/value"//&
                     &"/anisotropic_symmetric").or.&
                     have_option(trim(velocity_path)//&
                     &"/tensor_field::Viscosity/prescribed/value"//&
                     &"/anisotropic_asymmetric")) then
                   ewrite(-1,*) "WARNING - At present, the viscosity scaling for the pressure mass matrix is"
                   ewrite(-1,*) "taken from the 1st component of the viscosity tensor. Such a scaling"
                   ewrite(-1,*) "is only valid when all components of each viscosity tensor are constant."
                end if
             case("NoPreconditionerMatrix")
                ! Check no preconditioner is selected when no preconditioner matrix is desired:
                call get_option("/material_phase["//int2str(i)//&
                     "]/scalar_field::Pressure/prognostic/solver/preconditioner/name", pc_type)
                if(pc_type /= 'none') FLExit("If no preconditioner is desired, set pctype='none'.")
             end select
             
             ! Check inner matrix is valid for Stokes - must have full viscous terms
             ! included. Stokes does not have a mass matrix.
             call get_option(trim(schur_path)//"/inner_matrix[0]/name", inner_matrix)
             
             if(trim(inner_matrix)/="FullMomentumMatrix") then
                ewrite(-1,*) "For Stokes problems, FullMomentumMatrix must be specified under:"
                ewrite(-1,*) "scalar_field::Pressure/prognostic/scheme/use_projection_method& "
                ewrite(-1,*) "&/full_schur_complement/inner_matrix"
                FLExit("For Stokes problems, change --> FullMomentumMatrix")
             end if

          end if
          
       end if
       
    end do

  end subroutine check_stokes_options

  subroutine check_implicit_solids_options

    integer :: nmat, i
    logical :: have_scon, have_spha, have_oneway, have_twoway

    nmat = option_count("/material_phase")

    do i = 0, nmat-1
       have_scon = have_option("/material_phase["//int2str(i)//&
            "]/scalar_field::SolidConcentration")
       have_spha = have_option("/material_phase["//int2str(i)//&
            "]/scalar_field::SolidPhase")
       if((.not.have_scon).or.(.not.have_spha)) then
          FLExit("An implicit solid needs a SolidConcentration and a SolidPhase.")
       end if
    end do
    
    have_oneway = have_option("/material_phase/one_way_coupling")
    have_twoway = have_option("/material_phase/two_way_coupling")

    if((.not.have_oneway).or.(.not.have_twoway)) then
       FLExit("Implicit_solids should be run with either a one-way coupling or a two-way coupling.")
    end if
       
  end subroutine check_implicit_solids_options

  subroutine check_foams_options
    ! Check options for liquid drainage in foam simulations.

    character(len=OPTION_PATH_LEN) :: velocity_path, pressure_path, drainage_lambda_path, compressible_eos_path, foam_velocity_path
    logical :: exclude_mass, equation_drainage, compressible_projection, prescribed_lambda, foam_eos, foam_velocity, Drainage_K1, Drainage_K2, source, absorption

    ! Check that the local length of Plateau borders per unit volume (lambda) is provided.
    compressible_eos_path="/material_phase[0]/equation_of_state/compressible"
    foam_eos = have_option(trim(compressible_eos_path)//&
            "/foam")
    if(.not.(foam_eos)) then
       FLExit("The first material_phase in a foam problem must have a foam equation of state.")
    end if


    pressure_path="/material_phase[0]/scalar_field::Pressure/prognostic"
    if (have_option(trim(pressure_path))) then

       ! Check that compressible projection method is used:
       compressible_projection = have_option("/material_phase[0]"//&
            "/equation_of_state/compressible")

       if(.not.(compressible_projection)) then
          FLExit("For foam problems you need to use a compressible eos.")
       end if
    end if


    velocity_path="/material_phase[0]/vector_field::Velocity/prognostic"
    if (have_option(trim(velocity_path))) then

       ! Check that the equation type for drainage of liquid in foams is selected:
       equation_drainage = have_option(trim(velocity_path)//&
            "/equation::Drainage")

       if(.not.(equation_drainage)) then
          FLExit("For foam problems you need to select Drainage as your equation type.")
       end if

       ! Check that the mass term is excluded:
       exclude_mass = have_option(trim(velocity_path)//&
            "/spatial_discretisation&
            &/continuous_galerkin/mass_terms&
            &/exclude_mass_terms").or.&
                      have_option(trim(velocity_path)//&
            "/spatial_discretisation&
            &/discontinuous_galerkin/mass_terms&
            &/exclude_mass_terms")

       if(.not.(exclude_mass)) then
          FLExit("For foam problems you need to exclude the mass term.")
       end if

       ! Check that source and absorption are provided:
       source = have_option(trim(velocity_path)//&
             "/vector_field::Source")

       if(.not.(source)) then
          FLExit("You need a velocity source term for foam simulations.")
       end if

       absorption = have_option(trim(velocity_path)//&
             "/vector_field::Absorption")

       if(.not.(absorption)) then
          FLExit("You need a velocity absorption term for foam simulations.")
       end if

       ! Check that K1 and K2 fields are provided:
       Drainage_K1 = have_option(trim(velocity_path)//&
             "/vector_field::DrainageK1")

       if(.not.(Drainage_K1)) then
          FLExit("You need DrainageK1 vector field for foam simulations.")
       end if

       Drainage_K2 = have_option(trim(velocity_path)//&
             "/scalar_field::DrainageK2")

       if(.not.(Drainage_K2)) then
          FLExit("You need DrainageK2 scalar field for foam simulations.")
       end if

    end if

    ! Check that there is a Foam Velocity field.
    foam_velocity_path="/material_phase[0]/vector_field::FoamVelocity"
    foam_velocity = have_option(trim(foam_velocity_path)//&
            "/prescribed").or.&
                    have_option(trim(foam_velocity_path)//&
            "/diagnostic")
    if(.not.(foam_velocity)) then
       FLExit("For foam simulations you need either a prescribed or a diagnostic Foam Velocity field.")
    end if

    ! Check that the local length of Plateau borders per unit volume (lambda) is provided.
    drainage_lambda_path="/material_phase[0]/scalar_field::DrainageLambda"
    prescribed_lambda = have_option(trim(drainage_lambda_path)//&
            "/prescribed")
    if(.not.(prescribed_lambda)) then
       FLExit("For foam simulations you need a DrainageLambda field which at the moment must be prescribed.")
    end if

  end subroutine check_foams_options
  
  subroutine check_multiphase_options
    !!< Options checking for multi-phase flow simulations.
    ! This currently assumes that all phases have prognostic velocity fields;
    ! we will deal with prescribed velocities later.

    integer :: nmat, i
    logical :: have_vfrac, prognostic_velocity
    
    integer :: diagnostic_vfrac_count

    nmat = option_count("/material_phase")

    do i = 0, nmat-1
       have_vfrac = have_option("/material_phase["//int2str(i)//&
            "]/scalar_field::PhaseVolumeFraction")
       prognostic_velocity = have_option("/material_phase["//int2str(i)//&
            "]/vector_field::Velocity/prognostic")
       if(prognostic_velocity .and. .not.have_vfrac) then
          FLExit("All phases need a PhaseVolumeFraction.")
       end if
    end do
    
    diagnostic_vfrac_count = option_count(&
                &'/material_phase/scalar_field::PhaseVolumeFraction/diagnostic')
    if(diagnostic_vfrac_count > 1) then
      ewrite(-1,*) diagnostic_vfrac_count, 'diagnostic PhaseVolumeFractions.'
      FLExit("Only 1 diagnostic PhaseVolumeFraction is allowed")
    end if

  end subroutine check_multiphase_options

end module populate_state_module
