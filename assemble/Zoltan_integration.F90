#include "fdebug.h"
#include "confdefs.h"

module zoltan_integration

#ifdef HAVE_ZOLTAN
  use spud
  use fldebug
  use global_parameters, only: real_size, OPTION_PATH_LEN, topology_mesh_name,&
      FIELD_NAME_LEN
  use futils, only: int2str, present_and_true
  use quadrature
  use element_numbering, only: ele_local_num
  use elements
  use mpi_interfaces
  use data_structures
  use parallel_tools
  use memory_diagnostics
  use sparse_tools
  use linked_lists
  use halos_ownership
  use parallel_fields
  use metric_tools
  use fields
  use state_module
  use field_options
  use vtk_interfaces
  use zoltan
  use halos_derivation
  use halos
  use sparsity_patterns_meshes
  use reserve_state_module
  use boundary_conditions
  use c_interfaces
  use detector_data_types
  use boundary_conditions_from_options
  use pickers
  use detector_tools
  use detector_parallel
  use hadapt_advancing_front
  use fields_halos
  use populate_state_module
  use surface_id_interleaving
  use adapt_integration
  use zoltan_global_variables
  use zoltan_detectors
  use zoltan_callbacks

  implicit none

  integer, save :: max_coplanar_id 

  public :: zoltan_drive
  private

  contains

  subroutine zoltan_drive(states, final_adapt_iteration, global_min_quality, metric, full_metric, initialise_fields, &
    skip_extrusion_after, skip_extruded_mesh_migration, flredecomping, input_procs, target_procs)

    type(state_type), dimension(:), intent(inout), target :: states
    logical, intent(in) :: final_adapt_iteration
    ! returns the minimum element quality. When using libadapitivity (instead of mba2d/3d),
    ! it is based on the minimum nodal quality where the nodal quality is computed from the
    ! *maximum* quality of the adjacent elements at each node - as this is closer to libadaptivity's
    ! termination criterion
    real, intent(out), optional :: global_min_quality
    ! the metric is the metric we base the quality functions on
    type(tensor_field), intent(inout), optional :: metric
    ! the full_metric is the metric we need to interpolate
    type(tensor_field), intent(inout), optional :: full_metric
    ! if present and true: don't bother redistributing fields that can be reinitialised
    logical, intent(in), optional :: initialise_fields
    ! if present and true: don't extrude meshes after decomposition
    logical, intent(in), optional :: skip_extrusion_after
    ! if present and true: only decompose and migrate the horizontal mesh
    logical, intent(in), optional :: skip_extruded_mesh_migration
    ! Are we flredecomping? If so, this should be true
    logical, intent(in), optional :: flredecomping
    ! If flredecomping then these values should be provided
    integer, intent(in), optional :: input_procs, target_procs

    type(zoltan_struct), pointer :: zz

    logical :: changes
    integer(zoltan_int) :: num_gid_entries, num_lid_entries
    integer(zoltan_int), dimension(:), pointer :: p1_export_global_ids => null()
    integer(zoltan_int), dimension(:), pointer :: p1_export_local_ids => null()
    integer(zoltan_int), dimension(:), pointer :: p1_export_procs => null()
    integer(zoltan_int) :: p1_num_import, p1_num_export
    integer(zoltan_int), dimension(:), pointer :: p1_import_global_ids => null()
    integer(zoltan_int), dimension(:), pointer :: p1_import_local_ids => null()
    integer(zoltan_int), dimension(:), pointer :: p1_import_procs => null()
    integer, save :: dumpno = 0

    type(tensor_field) :: new_metric
    type(mesh_type), pointer :: full_mesh
    
    integer(zoltan_int), dimension(:), pointer :: p1_export_local_ids_full => null()
    integer(zoltan_int), dimension(:), pointer :: p1_export_procs_full => null()
    integer(zoltan_int) :: p1_num_export_full
    type(vector_field) :: zoltan_global_new_positions_m1d
    real :: load_imbalance_tolerance
    logical :: flredecomp
    real :: minimum_quality
    integer :: flredecomp_input_procs = -1, flredecomp_target_procs = -1

    ewrite(1,*) "In zoltan_drive"

    if (.not. present(flredecomping)) then
        flredecomp = .false.
    else
        flredecomp = flredecomping
    end if

    if (flredecomp) then
       ! check for required optional arguments
       if (present(input_procs)) then
          flredecomp_input_procs = input_procs
       else
          FLAbort("input_procs must be supplied when flredecomping.")
       end if
       if (present(target_procs)) then
          flredecomp_target_procs = target_procs
       else
          FLAbort("target_procs must be supplied when flredecomping.")
       end if

       zoltan_global_base_option_path = '/flredecomp'
    else
       ! check invalid optional arguments haven't been supplied
       if (present(input_procs)) then
          FLAbort("input_procs should only be provided when flredecomping.")
       end if
       if (present(target_procs)) then
          FLAbort("target_procs should only be provided when flredecomping.")
       end if
       
       zoltan_global_base_option_path = '/mesh_adaptivity/hr_adaptivity/zoltan_options'
    end if

    zoltan_global_field_weighted_partitions = &
     have_option(trim(zoltan_global_base_option_path) // "/field_weighted_partitions")

    call setup_module_variables(states, final_adapt_iteration, zz, flredecomp)

    call setup_quality_module_variables(metric, minimum_quality) ! this needs to be called after setup_module_variables
                                        ! (but only on the 2d mesh with 2+1d adaptivity)
    if (present(global_min_quality)) then
       if (.NOT. final_adapt_iteration) then
          global_min_quality = minimum_quality
       else
          ! On final iteration we do not calculate the minimum element quality
          global_min_quality = 1.0
       end if
    end if

    full_mesh => extract_mesh(states(1), trim(topology_mesh_name))
    zoltan_global_migrate_extruded_mesh = (mesh_dim(full_mesh) /= mesh_dim(zoltan_global_zz_mesh)) .and. &
      .not. present_and_true(skip_extruded_mesh_migration)

    if(zoltan_global_migrate_extruded_mesh .AND. zoltan_global_field_weighted_partitions) then
        ewrite(-1,*) "Cannot weight mesh partitions based upon extruded columns"// &
                     "and a prescribed field. Select one option only or fix the code."
        FLExit("Use Weighted mesh partitions for EITHER extruded meshes or prescribed fields")
    end if
    
    if(zoltan_global_migrate_extruded_mesh) then
       call create_columns_sparsity(zoltan_global_columns_sparsity, full_mesh)
    end if

    load_imbalance_tolerance = get_load_imbalance_tolerance(final_adapt_iteration)
    call set_zoltan_parameters(final_adapt_iteration, flredecomp, flredecomp_target_procs, load_imbalance_tolerance, zz)

    call zoltan_load_balance(zz, changes, num_gid_entries, num_lid_entries, &
       & p1_num_import, p1_import_global_ids, p1_import_local_ids, p1_import_procs, & 
       & p1_num_export, p1_export_global_ids, p1_export_local_ids, p1_export_procs, &
       & load_imbalance_tolerance, flredecomp, flredecomp_input_procs, flredecomp_target_procs)

    if (.not. changes) then
      ewrite(1,*) "Zoltan decided no change was necessary, exiting"
      call deallocate_zoltan_lists(p1_import_global_ids, p1_import_local_ids, p1_import_procs, &
           & p1_export_global_ids, p1_export_local_ids, p1_export_procs)
      call cleanup_basic_module_variables(zz)
      call cleanup_quality_module_variables
      dumpno = dumpno + 1

      if (final_adapt_iteration) then
        ! interpolation does not interpolate in the halo regions, so we need a halo update afterwards
        ! normally this happens automatically due to the subsequent zoltan migration process, however
        ! if zoltan decides to not do anything we need to do it manually. We only need it in the final one
        ! because interpolation does halo update the old fields before the adapt.
        call halo_update(states)
      end if
      return
    end if

    if(zoltan_global_migrate_extruded_mesh) then
      call derive_full_export_lists(states, p1_num_export, p1_export_local_ids, p1_export_procs, &
           & p1_num_export_full, p1_export_local_ids_full, p1_export_procs_full)
    end if

    ! The general plan:
    ! Just send the nodes you own, along with a note of their dependencies
    ! The receiving process loops through all its receives and records whom
    ! it needs to receive from (the OLD owner)

    ! It builds an import list from that, then migrates again

    call are_we_keeping_or_sending_nodes(p1_num_export, p1_export_local_ids, p1_export_procs)

    ! Migrate here
    ! for nodes I am going to own
    call zoltan_migration_phase_one(zz, & 
         & p1_num_import, p1_import_global_ids, p1_import_local_ids, p1_import_procs, &
         & p1_num_export, p1_export_global_ids, p1_export_local_ids, p1_export_procs)
    call deallocate_zoltan_lists(p1_import_global_ids, p1_import_local_ids, p1_import_procs, &
         & p1_export_global_ids, p1_export_local_ids, p1_export_procs)
    ! deal with reconstructing new mesh, positions, etc. for processes who are only exporting
    call deal_with_exporters
    ! for halo nodes those nodes depend on
    call zoltan_migration_phase_two(zz)
    call deallocate_my_lists

    call reconstruct_enlist
    call reconstruct_senlist
    call reconstruct_halo(zz)
    
    if(zoltan_global_migrate_extruded_mesh) then
      zoltan_global_new_positions_m1d = zoltan_global_new_positions ! save a reference to the horizontal mesh you've just load balanced
      call copy(zoltan_global_universal_to_new_local_numbering_m1d, zoltan_global_universal_to_new_local_numbering)
      
      call cleanup_basic_module_variables(zz)
      ! don't clean up the quality variables now 
      ! (they'll be deallocated later but we don't need to use them in the vertically_structured section
      ! so we don't need to reallocate them either)
      call cleanup_other_module_variables
      
      call setup_module_variables(states, final_adapt_iteration, zz, flredecomp, mesh_name = topology_mesh_name)


      load_imbalance_tolerance = get_load_imbalance_tolerance(final_adapt_iteration)
      call set_zoltan_parameters(final_adapt_iteration, flredecomp, flredecomp_target_procs, load_imbalance_tolerance, zz)

      call reset_zoltan_lists_full(zz, &
       & p1_num_export_full, p1_export_local_ids_full, p1_export_procs_full, &
       & p1_num_import, p1_import_global_ids, p1_import_local_ids, p1_import_procs, &
       & p1_num_export, p1_export_global_ids, p1_export_local_ids, p1_export_procs)
      
      ! It builds an import list from that, then migrates again

      call are_we_keeping_or_sending_nodes(p1_num_export, p1_export_local_ids, p1_export_procs)

      ! Migrate here
      ! for nodes I am going to own
      call zoltan_migration_phase_one(zz, & 
       & p1_num_import, p1_import_global_ids, p1_import_local_ids, p1_import_procs, &
       & p1_num_export, p1_export_global_ids, p1_export_local_ids, p1_export_procs)

      call deallocate_zoltan_lists(p1_import_global_ids, p1_import_local_ids, p1_import_procs, &
           & p1_export_global_ids, p1_export_local_ids, p1_export_procs)

      call deal_with_exporters
      ! for halo nodes those nodes depend on
      call zoltan_migration_phase_two(zz)
      call deallocate_my_lists

      call reconstruct_enlist
      call reconstruct_senlist
      call reconstruct_halo(zz)

      if (.not. verify_consistent_local_element_numbering(zoltan_global_new_positions%mesh) ) then
        ewrite(-1,*) "For the extruded mesh, the local element numbering of elements in the halo region" // &
                     "is not consistent with that of the element owner. This is likely" // & 
                     "due to a bug in zoltan. Please report" // &
                     "to the fluidity mailing list"
        FLExit("Need a consistent local element ordering in parallel")
      end if
      
      deallocate(zoltan_global_universal_columns)
      call deallocate(zoltan_global_universal_to_new_local_numbering_m1d)
    end if

    ! At this point, we now have the balanced linear external mesh.
    ! Get populate_state to allocate the fields and such on this new
    ! mesh.

    call initialise_transfer(zz, states, zoltan_global_new_positions_m1d, metric, full_metric, new_metric, initialise_fields, skip_extrusion_after)

    ! And now transfer the field data around.
    call transfer_fields(zz)

    call deallocate(zoltan_global_new_positions)
    if(zoltan_global_migrate_extruded_mesh) then
      call deallocate(zoltan_global_new_positions_m1d)
    end if

    call finalise_transfer(states, metric, full_metric, new_metric)

    call cleanup_basic_module_variables(zz)
    call cleanup_quality_module_variables
    call cleanup_other_module_variables
    
    dumpno = dumpno + 1

    ewrite(1,*) "Exiting zoltan_drive"

  end subroutine zoltan_drive

  subroutine setup_module_variables(states, final_adapt_iteration, zz, flredecomp, mesh_name)
    type(state_type), dimension(:), intent(inout), target :: states
    logical, intent(in) :: final_adapt_iteration
    logical, intent(in) :: flredecomp
    type(zoltan_struct), pointer, intent(out) :: zz
    
    type(mesh_type), pointer :: mesh_ptr
    character(len=*), optional :: mesh_name
    integer :: nhalos, stat
    integer, dimension(:), allocatable :: owned_nodes
    integer :: i, j, floc, eloc
    integer, dimension(:), allocatable :: face_nodes
    integer :: old_element_number, universal_element_number, face_number, universal_surface_element_number
    integer, dimension(:), allocatable :: interleaved_surface_ids

    if (final_adapt_iteration) then
       zoltan_global_calculate_edge_weights = .false.
    else
       zoltan_global_calculate_edge_weights = .true.
    end if
    
    zoltan_global_max_edge_weight_on_node => extract_scalar_field(states(1), "MaxEdgeWeightOnNodes", stat) 
    if (stat == 0) then
       zoltan_global_output_edge_weights = .true.
    end if
    
    ! set quality_tolerance
    if (have_option(trim(zoltan_global_base_option_path) // "/element_quality_cutoff")) then
       call get_option(trim(zoltan_global_base_option_path) // "/element_quality_cutoff", zoltan_global_quality_tolerance)
       ! check that the value is reasonable
       if (zoltan_global_quality_tolerance < 0. .or. zoltan_global_quality_tolerance > 1.) then
          FLExit("element_quality_cutoff should be between 0 and 1. Default is 0.6")
       end if
    else
       zoltan_global_quality_tolerance = 0.6
    end if

    if(present(mesh_name)) then
       zoltan_global_zz_mesh = extract_mesh(states(1), trim(mesh_name))
    else if (flredecomp) then
       zoltan_global_zz_mesh = get_external_mesh(states)
    else
       ! This should actually be using get_external_mesh(), i.e. zoltan redistributes the mesh
       ! that all other meshes are derived from. However, in the case that we extrude and use
       ! generic adaptivity (not 2+1), the external horizontal mesh becomes seperated from the
       ! other meshes and it's actually the 3d adapted mesh that we derive everything from. We
       ! should probably actually completely remove the external horizontal mesh from the options
       ! tree, and make the adapted 3d mesh the external mesh. For now we keep it around and leave
       ! it in its old decomposition.
       call find_mesh_to_adapt(states(1), mesh_ptr)
       zoltan_global_zz_mesh = mesh_ptr
    end if
    call incref(zoltan_global_zz_mesh)
    if (zoltan_global_zz_mesh%name=="CoordinateMesh") then
       zoltan_global_zz_positions = extract_vector_field(states, "Coordinate")
    else
       zoltan_global_zz_positions = extract_vector_field(states, trim(zoltan_global_zz_mesh%name)//"Coordinate")
    end if
    call incref(zoltan_global_zz_positions)
    
    zoltan_global_zz_nelist => extract_nelist(zoltan_global_zz_mesh)
    
    zz => Zoltan_Create(halo_communicator(zoltan_global_zz_mesh))
    
    nhalos = halo_count(zoltan_global_zz_mesh)
    assert(nhalos == 2)
    zoltan_global_zz_halo => zoltan_global_zz_mesh%halos(nhalos)
    
    nhalos = element_halo_count(zoltan_global_zz_mesh)
    assert(nhalos >= 1)
    zoltan_global_zz_ele_halo => zoltan_global_zz_mesh%element_halos(nhalos)
    
    zoltan_global_zz_sparsity_one => get_csr_sparsity_firstorder(states, zoltan_global_zz_mesh, zoltan_global_zz_mesh)
    zoltan_global_zz_sparsity_two => get_csr_sparsity_secondorder(states, zoltan_global_zz_mesh, zoltan_global_zz_mesh)
    
    allocate(owned_nodes(halo_nowned_nodes(zoltan_global_zz_halo)))
    call allocate(zoltan_global_universal_to_old_local_numbering)
    call get_owned_nodes(zoltan_global_zz_halo, owned_nodes)
    do i=1,size(owned_nodes)
       call insert(zoltan_global_universal_to_old_local_numbering, halo_universal_number(zoltan_global_zz_halo, owned_nodes(i)), owned_nodes(i))
    end do
    deallocate(owned_nodes)

    call allocate(zoltan_global_uen_to_old_local_numbering)
    call allocate(zoltan_global_old_local_numbering_to_uen)
    do i=1,ele_count(zoltan_global_zz_positions)
       call insert(zoltan_global_uen_to_old_local_numbering, halo_universal_number(zoltan_global_zz_ele_halo, i), i)
       call insert(zoltan_global_old_local_numbering_to_uen, i, halo_universal_number(zoltan_global_zz_ele_halo, i))
    end do
    
    allocate(zoltan_global_receives(halo_proc_count(zoltan_global_zz_halo)))
    do i=1,size(zoltan_global_receives)
       call allocate(zoltan_global_receives(i))
    end do

    ! set up zoltan_global_old_snelist
    allocate(zoltan_global_old_snelist(node_count(zoltan_global_zz_positions)))
    call allocate(zoltan_global_old_snelist)
    call allocate(zoltan_global_universal_surface_number_to_surface_id)
    call allocate(zoltan_global_universal_surface_number_to_element_owner)
    allocate(interleaved_surface_ids(surface_element_count(zoltan_global_zz_positions)))
    call interleave_surface_ids(zoltan_global_zz_mesh, interleaved_surface_ids, max_coplanar_id)

    ! this is another thing that needs to be generalised for mixed meshes
    floc = face_loc(zoltan_global_zz_positions, 1)
    eloc = ele_loc(zoltan_global_zz_positions, 1)
    allocate(face_nodes(1:floc))
    
    do i=1, surface_element_count(zoltan_global_zz_positions)
       old_element_number = face_ele(zoltan_global_zz_positions, i)
       universal_element_number = halo_universal_number(zoltan_global_zz_ele_halo, old_element_number)
       face_number = local_face_number(zoltan_global_zz_positions, i)
       universal_surface_element_number = (universal_element_number-1)*eloc + face_number
       
       call insert(zoltan_global_universal_surface_number_to_surface_id, universal_surface_element_number, interleaved_surface_ids(i))
       call insert(zoltan_global_universal_surface_number_to_element_owner, universal_surface_element_number, universal_element_number)

       face_nodes = face_global_nodes(zoltan_global_zz_mesh, i)
       do j=1, floc
          call insert(zoltan_global_old_snelist(face_nodes(j)), universal_surface_element_number)
       end do
    end do
    
    deallocate(interleaved_surface_ids)
    deallocate(face_nodes)
    
    zoltan_global_preserve_mesh_regions = associated(zoltan_global_zz_mesh%region_ids)
    ! this deals with the case where some processors have no elements
    ! (i.e. when used to flredecomp from 1 to many processors)
    call allor(zoltan_global_preserve_mesh_regions) 
    if(zoltan_global_preserve_mesh_regions) then
       call allocate(zoltan_global_universal_element_number_to_region_id)
       do i = 1, element_count(zoltan_global_zz_positions)
          universal_element_number = halo_universal_number(zoltan_global_zz_ele_halo, i)
          call insert(zoltan_global_universal_element_number_to_region_id, universal_element_number, zoltan_global_zz_positions%mesh%region_ids(i))
       end do
    end if

    if(zoltan_global_field_weighted_partitions) then
       zoltan_global_field_weighted_partition_values = extract_scalar_field(states, "FieldWeightedPartitionValues") 
       assert(zoltan_global_field_weighted_partition_values%mesh == zoltan_global_zz_mesh)

       if(zoltan_global_field_weighted_partition_values%mesh%name /= zoltan_global_zz_mesh%name) then
          ewrite(-1,*) "FieldWeightedPartitionValues and Zoltan Global ZZ Mesh must be on the " // &
                       "same mesh. 99.9% of the time, this means that FieldWeightedPartitionValues " // & 
                       "must be on the external mesh."
          FLExit("FieldWeightedPartitionValues must be on the external mesh")
       end if

       call incref(zoltan_global_field_weighted_partition_values)

    end if

  end subroutine setup_module_variables

  subroutine setup_quality_module_variables(metric, minimum_quality)
    ! setups the field zoltan_global_element_quality (used to determine edge weights)
    ! and returns minimum_quality (to be used as zoltan iteration termination criterion)
    ! the metric is the metric we base the quality functions on
    type(tensor_field), intent(in), optional :: metric
    ! returns the minimum element quality. When using libadapitivity (instead of mba2d/3d),
    ! it is based on the minimum nodal quality where the nodal quality is computed from the
    ! *maximum* quality of the adjacent elements at each node - as this is closer to libadaptivity's
    ! termination criterion
    real, intent(out):: minimum_quality
    
    type(mesh_type) :: pwc_mesh
    integer :: node
    integer, dimension(:), pointer :: elements
    logical :: use_pain_functional
    
    ! And the element quality measure
    use_pain_functional = present(metric) .and. mesh_dim(zoltan_global_zz_mesh)==3 .and. &
          .not. have_option("/mesh_adaptivity/hr_adaptivity/adaptivity_library/libmba3d")
    if (use_pain_functional) then
       ! with libadaptivity use the Pain functional
       call element_quality_pain_p0(zoltan_global_zz_positions, metric, zoltan_global_element_quality)
       ! the rest of the zoltan wrappers have been written assuming the lipnikov 
       ! functional where q=0 is bad and q=1 is perfect. With the pain functional, 
       ! q'=0 is perfect and q'=\infty is bad. Therefore map q' -> q=1/(q'+1), so
       ! that we get the same behaviour
       call addto(zoltan_global_element_quality, 1.0)
       call invert(zoltan_global_element_quality)
    else if (present(metric)) then
       ! with mba2d or mba3d use the lipnikov functional:
       call element_quality_p0(zoltan_global_zz_positions, metric, zoltan_global_element_quality)
    else
       pwc_mesh = piecewise_constant_mesh(zoltan_global_zz_mesh, "PWCMesh")
       call allocate(zoltan_global_element_quality, pwc_mesh, "ElementQuality", field_type=FIELD_TYPE_CONSTANT)
       call set(zoltan_global_element_quality, 1.0)
       call deallocate(pwc_mesh)
    end if

    minimum_quality = minval(zoltan_global_element_quality)
    ewrite(1,*) "local minimum element quality = ", minimum_quality
    call allmin(minimum_quality)
    ewrite(1,*) "global minimum element quality = ", minimum_quality

    if (use_pain_functional) then
       ! libadaptivity terminates if for all possible operations, any of the affected
       ! elements have a quality that is above a threshold. This means that if an element
       ! that hasn't reached the threshold yet can only be improved via operations that
       ! affect a neighbour that is already good enough - it will be kept as it is.
       ! Therefore we compute the best quality element adjacent to each node, and then
       ! take the minimum over all nodes. The zoltan iterations termination criterion isbased
       ! on this minimum, saying that if all nodes should have at least a good enough element adjacent 
       ! to it - as we can't guarantee that elements of worse quality adjacent to a node will
       ! ever be changed by libadaptivity
       minimum_quality = 1.0
       do node=1, size(zoltan_global_zz_nelist, 1)
          elements => row_m_ptr(zoltan_global_zz_nelist, node)
          minimum_quality = min(minimum_quality, maxval(node_val(zoltan_global_element_quality, elements)))
       end do
       ewrite(1,*) "local minimum achievable quality = ", minimum_quality
       call allmin(minimum_quality)
       ewrite(1,*) "global minimum achievable quality = ", minimum_quality
    end if
    
  end subroutine setup_quality_module_variables

  function get_load_imbalance_tolerance(final_adapt_iteration) result(load_imbalance_tolerance)
    logical, intent(in) :: final_adapt_iteration    
 
    real, parameter :: default_load_imbalance_tolerance = 1.05
    real, parameter :: final_iteration_load_imbalance_tolerance = 1.02
    real :: load_imbalance_tolerance

    if (.NOT. final_adapt_iteration) then
       ! if user has passed us the option then use the load imbalance tolerance they supplied,
       ! else use the default load imbalance tolerance
       call get_option(trim(zoltan_global_base_option_path) // "/load_imbalance_tolerance", load_imbalance_tolerance, &
          & default = default_load_imbalance_tolerance)
       
       ! check the value is reasonable
       if (load_imbalance_tolerance < 1.0) then
          FLExit("load_imbalance_tolerance should be greater than or equal to 1. Default is 1.5")
       end if
    else
       load_imbalance_tolerance = final_iteration_load_imbalance_tolerance
    end if

  end function get_load_imbalance_tolerance

  subroutine set_zoltan_parameters(final_adapt_iteration, flredecomp, target_procs, &
     & load_imbalance_tolerance, zz)
    logical, intent(in) :: final_adapt_iteration
    logical, intent(in) :: flredecomp
    integer, intent(in) :: target_procs
    real, intent(in) :: load_imbalance_tolerance
    type(zoltan_struct), pointer, intent(in) :: zz    

    integer(zoltan_int) :: ierr
    character (len = FIELD_NAME_LEN) :: method, graph_checking_level
    character (len = 10) :: string_load_imbalance_tolerance

    if (debug_level()>1) then
       ierr = Zoltan_Set_Param(zz, "DEBUG_LEVEL", "1"); assert(ierr == ZOLTAN_OK)
    else         
       ierr = Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0"); assert(ierr == ZOLTAN_OK)
    end if

    ! convert load_imbalance_tolerance to a string for setting the option in Zoltan
    write(string_load_imbalance_tolerance, '(f6.3)' ) load_imbalance_tolerance
    ierr = Zoltan_Set_Param(zz, "IMBALANCE_TOL", string_load_imbalance_tolerance); assert(ierr == ZOLTAN_OK)
    ewrite(2,*) 'Initial load_imbalance_tolerance set to ', load_imbalance_tolerance

    ! For flredecomp if we are not an active process, then let's set the number of local parts to be zero
    if (flredecomp) then
       if (getprocno() > target_procs) then
          ierr = Zoltan_Set_Param(zz, "NUM_LOCAL_PARTS", "0"); assert(ierr == ZOLTAN_OK)
       else
          ierr = Zoltan_Set_Param(zz, "NUM_LOCAL_PARTS", "1"); assert(ierr == ZOLTAN_OK)
       end if
       ierr = Zoltan_set_Param(zz, "NUM_GLOBAL_PARTS", int2str(target_procs)); assert(ierr == ZOLTAN_OK)
    end if
    
    if (.NOT. final_adapt_iteration) then
       if (have_option(trim(zoltan_global_base_option_path) // "/partitioner")) then
          if (have_option(trim(zoltan_global_base_option_path) // "/partitioner/metis"))  then
             ierr = Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH"); assert(ierr == ZOLTAN_OK)
             ierr = Zoltan_Set_Param(zz, "GRAPH_PACKAGE", "PARMETIS"); assert(ierr == ZOLTAN_OK)
             ewrite(3,*) "Setting the partitioner to be ParMETIS."
             ! turn off graph checking unless debugging, this was filling the error file with Zoltan warnings
             if (have_option(trim(zoltan_global_base_option_path) // "/zoltan_debug/graph_checking")) then
                call get_option(trim(zoltan_global_base_option_path) // "/zoltan_debug/graph_checking", graph_checking_level)
                ierr = Zoltan_Set_Param(zz, "CHECK_GRAPH", trim(graph_checking_level)); assert(ierr == ZOLTAN_OK)
             else
                ierr = Zoltan_Set_Param(zz, "CHECK_GRAPH", "0"); assert(ierr == ZOLTAN_OK)
             end if
          end if
          
          if (have_option(trim(zoltan_global_base_option_path) // "/partitioner/zoltan")) then
             
             call get_option(trim(zoltan_global_base_option_path) // "/partitioner/zoltan/method", method)
             
             if (trim(method) == "graph") then
                ierr = Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH"); assert(ierr == ZOLTAN_OK)
                ierr = Zoltan_Set_Param(zz, "GRAPH_PACKAGE", "PHG"); assert(ierr == ZOLTAN_OK)
                ewrite(3,*) "Setting the partitioner to be Zoltan-Graph."
             else if (trim(method) == "hypergraph") then
                ierr = Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH"); assert(ierr == ZOLTAN_OK)
                ierr = Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG"); assert(ierr == ZOLTAN_OK)
                ewrite(3,*) "Setting the partitioner to be Zoltan-Hypergraph."
             end if
             
          end if
       
          if (have_option(trim(zoltan_global_base_option_path) // "/partitioner/scotch")) then
             ierr = Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH"); assert(ierr == ZOLTAN_OK)
             ierr = Zoltan_Set_Param(zz, "GRAPH_PACKAGE", "SCOTCH"); assert(ierr == ZOLTAN_OK)
             ewrite(3,*) "Setting the partitioner to be Scotch."
             ! Probably not going to want graph checking unless debugging
             if (have_option(trim(zoltan_global_base_option_path) // "/zoltan_debug/graph_checking")) then
                call get_option(trim(zoltan_global_base_option_path) // "/zoltan_debug/graph_checking", graph_checking_level)
                ierr = Zoltan_Set_Param(zz, "CHECK_GRAPH", trim(graph_checking_level)); assert(ierr == ZOLTAN_OK)
             else
                ierr = Zoltan_Set_Param(zz, "CHECK_GRAPH", "0"); assert(ierr == ZOLTAN_OK)
             end if
          end if
         
       else
          ! Use the Zoltan graph partitioner by default
          ierr = Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH"); assert(ierr == ZOLTAN_OK)
          ierr = Zoltan_Set_Param(zz, "GRAPH_PACKAGE", "PHG"); assert(ierr == ZOLTAN_OK)
          ewrite(3,*) "No partitioner option set, defaulting to using Zoltan-Graph."
       end if

    else

       if (have_option(trim(zoltan_global_base_option_path) // "/final_partitioner")) then
          
          if (have_option(trim(zoltan_global_base_option_path) // "/final_partitioner/metis"))  then
             ierr = Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH"); assert(ierr == ZOLTAN_OK)
             ierr = Zoltan_Set_Param(zz, "GRAPH_PACKAGE", "PARMETIS"); assert(ierr == ZOLTAN_OK)
             ewrite(3,*) "Setting the final partitioner to be ParMETIS."
             ! turn off graph checking unless debugging, this was filling the error file with Zoltan warnings
             if (have_option(trim(zoltan_global_base_option_path) // "/zoltan_debug/graph_checking")) then
                call get_option(trim(zoltan_global_base_option_path) // "/zoltan_debug/graph_checking", graph_checking_level)
                ierr = Zoltan_Set_Param(zz, "CHECK_GRAPH", trim(graph_checking_level)); assert(ierr == ZOLTAN_OK)
             else
                ierr = Zoltan_Set_Param(zz, "CHECK_GRAPH", "0"); assert(ierr == ZOLTAN_OK)
             end if
          end if
          
          if (have_option(trim(zoltan_global_base_option_path) // "/final_partitioner/zoltan")) then
             
             call get_option(trim(zoltan_global_base_option_path) // "/final_partitioner/zoltan/method", method)
             
             if (trim(method) == "graph") then
                ierr = Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH"); assert(ierr == ZOLTAN_OK)
                ierr = Zoltan_Set_Param(zz, "GRAPH_PACKAGE", "PHG"); assert(ierr == ZOLTAN_OK)
                ewrite(3,*) "Setting the final partitioner to be Zoltan-Graph."
             else if (trim(method) == "hypergraph") then
                ierr = Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH"); assert(ierr == ZOLTAN_OK)
                ierr = Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG"); assert(ierr == ZOLTAN_OK)
                ewrite(3,*) "Setting the final partitioner to be Zoltan-Hypergraph."
             end if
             
          end if
       
          if (have_option(trim(zoltan_global_base_option_path) // "/final_partitioner/scotch")) then
             ierr = Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH"); assert(ierr == ZOLTAN_OK)
             ierr = Zoltan_Set_Param(zz, "GRAPH_PACKAGE", "SCOTCH"); assert(ierr == ZOLTAN_OK)
                ewrite(3,*) "Setting the final partitioner to be Scotch."
             ! Probably not going to want graph checking unless debugging
             if (have_option(trim(zoltan_global_base_option_path) // "/zoltan_debug/graph_checking")) then
                call get_option(trim(zoltan_global_base_option_path) // "/zoltan_debug/graph_checking", graph_checking_level)
                ierr = Zoltan_Set_Param(zz, "CHECK_GRAPH", trim(graph_checking_level)); assert(ierr == ZOLTAN_OK)
             else
                ierr = Zoltan_Set_Param(zz, "CHECK_GRAPH", "0"); assert(ierr == ZOLTAN_OK)
             end if
          end if
          
       else
          ! Use ParMETIS by default on the final adapt iteration
          ierr = Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH"); assert(ierr == ZOLTAN_OK)
          ierr = Zoltan_Set_Param(zz, "GRAPH_PACKAGE", "PARMETIS"); assert(ierr == ZOLTAN_OK)
          ewrite(3,*) "No final partitioner option set, defaulting to using ParMETIS."
       end if

    end if

    ! Choose the appropriate partitioning method based on the current adapt iteration.
    ! The default is currently to do a clean partition on all adapt iterations to produce a 
    ! load balanced partitioning and to limit the required number of adapts. In certain cases,
    ! repartitioning or refining may lead to improved performance, and this is optional for all
    ! but the final adapt iteration.
    if (final_adapt_iteration) then

       ierr = Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION"); assert(ierr == ZOLTAN_OK)
       ewrite(3,*) "Setting partitioning approach to PARTITION."
       if (have_option(trim(zoltan_global_base_option_path) // "/final_partitioner/metis") .OR. &
          & (.NOT.(have_option(trim(zoltan_global_base_option_path) // "/final_partitioner")))) then
          ! chosen to match what Sam uses
          ierr = Zoltan_Set_Param(zz, "PARMETIS_METHOD", "PartKway"); assert(ierr == ZOLTAN_OK)
          ewrite(3,*) "Setting ParMETIS method to PartKway."
       end if

    else

       if (have_option(trim(zoltan_global_base_option_path) // "/load_balancing_approach/")) then
          if (have_option(trim(zoltan_global_base_option_path) // "/load_balancing_approach/partition"))  then
             ! Partition from scratch, not taking the current data distribution into account:
             ierr = Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION"); assert(ierr == ZOLTAN_OK)
             ewrite(3,*) "Setting partitioning approach to PARTITION."
          else if (have_option(trim(zoltan_global_base_option_path) // "/load_balancing_approach/repartition"))  then
             ! Partition but try to stay close to the curruent partition/distribution:
             ierr = Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION"); assert(ierr == ZOLTAN_OK)
             ewrite(3,*) "Setting partitioning approach to REPARTITION."
          else if (have_option(trim(zoltan_global_base_option_path) // "/load_balancing_approach/refine"))  then
             ! Refine the current partition/distribution; assumes only small changes:
             ierr = Zoltan_Set_Param(zz, "LB_APPROACH", "REFINE"); assert(ierr == ZOLTAN_OK)
             ewrite(3,*) "Setting partitioning approach to REFINE."
          end if
       else
          ierr = Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION"); assert(ierr == ZOLTAN_OK)
          ewrite(3,*) "Setting partitioning approach to PARTITION."
       end if

       if (have_option(trim(zoltan_global_base_option_path) // "/partitioner/metis"))  then
          ! chosen to match what Sam uses
          ierr = Zoltan_Set_Param(zz, "PARMETIS_METHOD", "AdaptiveRepart"); assert(ierr == ZOLTAN_OK)
          ewrite(3,*) "Setting ParMETIS method to AdaptiveRepart."
          ierr = Zoltan_Set_Param(zz, "PARMETIS_ITR", "100000.0"); assert(ierr == ZOLTAN_OK)
       end if

    end if

    ierr = Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); assert(ierr == ZOLTAN_OK)
    ierr = Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1"); assert(ierr == ZOLTAN_OK)
    ierr = Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1"); assert(ierr == ZOLTAN_OK)
    ierr = Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "1"); assert(ierr == ZOLTAN_OK)
    ierr = Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); assert(ierr == ZOLTAN_OK)
    ! scotch and parmetis do not behave correctly when one, or more input partitions are empty
    ! setting this to 2, means doing a simple scatter of contiguous chunks of vertices (ignoring
    ! any weighting) before going into scotch/parmetis
    ! 0 = never scatter, 1 (default) = only scatter if there is just one non-empty partition, 3 = always scatter
    ! I believe this is ignored with Zoltan Hypergraph
    ! empty *input* partitions should only occur with flredecomp, as we don't handle empty output partitions
    ! and adaptivity never reduces the number of vertices to 0
    ierr = Zoltan_Set_Param(zz, "SCATTER_GRAPH", "2"); assert(ierr == ZOLTAN_OK)
    
    ierr = Zoltan_Set_Fn(zz, ZOLTAN_NUM_OBJ_FN_TYPE, zoltan_cb_owned_node_count);      assert(ierr == ZOLTAN_OK)
    ierr = Zoltan_Set_Fn(zz, ZOLTAN_OBJ_LIST_FN_TYPE, zoltan_cb_get_owned_nodes);      assert(ierr == ZOLTAN_OK)
    ierr = Zoltan_Set_Fn(zz, ZOLTAN_NUM_EDGES_MULTI_FN_TYPE, zoltan_cb_get_num_edges); assert(ierr == ZOLTAN_OK)
    ierr = Zoltan_Set_Fn(zz, ZOLTAN_EDGE_LIST_MULTI_FN_TYPE, zoltan_cb_get_edge_list); assert(ierr == ZOLTAN_OK)
    ierr = Zoltan_Set_Fn(zz, ZOLTAN_OBJ_SIZE_MULTI_FN_TYPE, zoltan_cb_pack_node_sizes); assert(ierr == ZOLTAN_OK)
    ierr = Zoltan_Set_Fn(zz, ZOLTAN_PACK_OBJ_MULTI_FN_TYPE, zoltan_cb_pack_nodes); assert(ierr == ZOLTAN_OK)
    ierr = Zoltan_Set_Fn(zz, ZOLTAN_UNPACK_OBJ_MULTI_FN_TYPE, zoltan_cb_unpack_nodes); assert(ierr == ZOLTAN_OK)
    
  end subroutine set_zoltan_parameters

  subroutine deallocate_zoltan_lists(p1_import_global_ids, p1_import_local_ids, p1_import_procs, &
       & p1_export_global_ids, p1_export_local_ids, p1_export_procs)
    
    integer(zoltan_int), dimension(:), pointer, intent(inout) :: p1_import_global_ids
    integer(zoltan_int), dimension(:), pointer, intent(inout) :: p1_import_local_ids
    integer(zoltan_int), dimension(:), pointer, intent(inout) :: p1_import_procs    
    integer(zoltan_int), dimension(:), pointer, intent(inout) :: p1_export_global_ids
    integer(zoltan_int), dimension(:), pointer, intent(inout) :: p1_export_local_ids
    integer(zoltan_int), dimension(:), pointer, intent(inout) :: p1_export_procs

    integer(zoltan_int) :: ierr

    ! deallocates the memory which was allocated in Zoltan
    ierr = Zoltan_LB_Free_Data(p1_import_global_ids, p1_import_local_ids, p1_import_procs, &
         & p1_export_global_ids, p1_export_local_ids, p1_export_procs)
    
    assert(ierr == ZOLTAN_OK)
  end subroutine deallocate_zoltan_lists

  subroutine cleanup_basic_module_variables(zz)
    ! This routine deallocates everything that is guaranteed to be allocated
    ! regardless of whether Zoltan actually wants to change anything or not.
    ! (except for the quality functions)
    type(zoltan_struct), pointer, intent(inout) :: zz    

    call deallocate(zoltan_global_universal_surface_number_to_surface_id)
    call deallocate(zoltan_global_universal_surface_number_to_element_owner)
    call deallocate(zoltan_global_old_snelist)
    deallocate(zoltan_global_old_snelist)
    call deallocate(zoltan_global_receives)
    deallocate(zoltan_global_receives)
    
    call deallocate(zoltan_global_universal_to_old_local_numbering)
    call deallocate(zoltan_global_uen_to_old_local_numbering)
    call deallocate(zoltan_global_old_local_numbering_to_uen)
    
    zoltan_global_new_positions%refcount => null()
    
    call deallocate(zoltan_global_zz_mesh)
    call deallocate(zoltan_global_zz_positions)
    zoltan_global_zz_sparsity_one => null()
    zoltan_global_zz_sparsity_two => null()
    zoltan_global_zz_halo => null()
    zoltan_global_zz_ele_halo => null()
    call Zoltan_Destroy(zz)
    
    zoltan_global_preserve_columns = .false.
  end subroutine cleanup_basic_module_variables
  
  subroutine cleanup_quality_module_variables
    ! This routine deallocates the module quality fields.
    call deallocate(zoltan_global_element_quality)
    if(zoltan_global_migrate_extruded_mesh) then
       call deallocate(zoltan_global_columns_sparsity)
    end if
    if(zoltan_global_field_weighted_partitions) then
       call deallocate(zoltan_global_field_weighted_partition_values)
    end if

  end subroutine cleanup_quality_module_variables

  subroutine cleanup_other_module_variables
    call deallocate(zoltan_global_nodes_we_are_sending)
    call deallocate(zoltan_global_nodes_we_are_keeping)
    call deallocate(zoltan_global_universal_to_new_local_numbering)
    call deallocate(zoltan_global_new_nodes)
    call deallocate(zoltan_global_uen_to_new_local_numbering)
  end subroutine cleanup_other_module_variables

  subroutine zoltan_load_balance(zz, changes, num_gid_entries, num_lid_entries, &
     & p1_num_import, p1_import_global_ids, p1_import_local_ids, p1_import_procs, &
     & p1_num_export, p1_export_global_ids, p1_export_local_ids, p1_export_procs, &
     load_imbalance_tolerance, flredecomp, input_procs, target_procs)
    
    type(zoltan_struct), pointer, intent(in) :: zz    
    logical, intent(out) :: changes    
    
    integer(zoltan_int), intent(out) :: num_gid_entries, num_lid_entries
    integer(zoltan_int), intent(out) :: p1_num_import
    integer(zoltan_int), dimension(:), pointer, intent(out) :: p1_import_global_ids
    integer(zoltan_int), dimension(:), pointer, intent(out) :: p1_import_local_ids
    integer(zoltan_int), dimension(:), pointer, intent(out) :: p1_import_procs
    integer(zoltan_int), intent(out) :: p1_num_export
    integer(zoltan_int), dimension(:), pointer, intent(out) :: p1_export_global_ids 
    integer(zoltan_int), dimension(:), pointer, intent(out) :: p1_export_local_ids
    integer(zoltan_int), dimension(:), pointer, intent(out) :: p1_export_procs
    real, intent(inout) :: load_imbalance_tolerance
    logical, intent(in) :: flredecomp
    integer, intent(in) :: input_procs, target_procs

    ! These variables are needed when flredecomping as we then use Zoltan_LB_Partition
    integer(zoltan_int), dimension(:), pointer :: import_to_part
    integer(zoltan_int), dimension(:), pointer :: export_to_part
    integer(zoltan_int), dimension(:), pointer :: null_pointer => null()

    integer(zoltan_int) :: ierr
    integer :: i, node
    integer :: num_nodes, num_nodes_after_balance
    integer :: min_num_nodes_after_balance, total_num_nodes_before_balance, total_num_nodes_after_balance
    integer :: num_empty_partitions, empty_partition
    character (len = 10) :: string_load_imbalance_tolerance

    ewrite(1,*) 'in zoltan_load_balance'

    num_nodes = zoltan_global_zz_halo%nowned_nodes

    ! Special case when flredecomping - don't check for empty partitions
    if (flredecomp) then

       ! calculate total number of owned nodes before the load balance
       call mpi_allreduce(num_nodes, total_num_nodes_before_balance, 1, getPINTEGER(), &
          & MPI_SUM, MPI_COMM_FEMTOOLS, ierr)

       ! Need to use Zoltan_LB_Partition when flredecomping as NUM_LOCAL_PART and NUM_GLOBAL_PART are
       ! meant to be invalid for Zoltan_LB_Balance (actually appear to be valid even then but better
       ! to follow the doc)
       ierr = Zoltan_LB_Partition(zz, changes, num_gid_entries, num_lid_entries, p1_num_import, p1_import_global_ids, &
          & p1_import_local_ids, p1_import_procs, import_to_part, p1_num_export, p1_export_global_ids,  &
          & p1_export_local_ids, p1_export_procs, export_to_part)
       assert(ierr == ZOLTAN_OK)

       ! calculate how many owned nodes we'd have after doing the planned load balancing
       num_nodes_after_balance = num_nodes + p1_num_import - p1_num_export

       ! calculate total number of owned nodes after the load balance
       call mpi_allreduce(num_nodes_after_balance, total_num_nodes_after_balance, 1, getPINTEGER(), &
          & MPI_SUM, MPI_COMM_FEMTOOLS, ierr)

       if (total_num_nodes_before_balance .NE. total_num_nodes_after_balance) then
          FLAbort("The total number of nodes before load balancing does not equal the total number of nodes after the load balancing.")
       end if

       if (target_procs < input_procs) then
          ! We're expecting some processes to have empty partitions when using flredecomp to
          ! reduce the number of active processes. The plan is to calculate the number of 
          ! processes with an empty partition after the load balance and check this is how
          ! many we'd expect to be empty
          if (num_nodes_after_balance > 0) then
             empty_partition = 0
          else
             empty_partition = 1
          end if
          call mpi_allreduce(empty_partition, num_empty_partitions, 1, getPINTEGER(), &
             & MPI_SUM, MPI_COMM_FEMTOOLS, ierr)
          
          if (num_empty_partitions /= (input_procs - target_procs)) then
             FLAbort("The correct number of processes did not have empty partitons after the load balancing.")
          end if
       else
          ! If using flredecomp to increase the number of active processes then no process
          ! should have an empty partition after the load balance
          if (num_nodes_after_balance == 0) then
             FLAbort("After load balancing process would have an empty partition.")
          end if
       end if

       if (p1_num_import>0) then
         ! It appears that with gcc5 this routine crashes if p1_num_import==0
         ! not entirely sure whether this is a bug in zoltan with gcc5 or
         ! whether we are indeed not suppposed to deallocate this if there are no imports
         ierr = Zoltan_LB_Free_Part(null_pointer, null_pointer, null_pointer, import_to_part); assert(ierr == ZOLTAN_OK)
       end if
       if (p1_num_export>0) then
         ! see comment above, p1_num_import -> p1_num_export
         ierr = Zoltan_LB_Free_Part(null_pointer, null_pointer, null_pointer, export_to_part); assert(ierr == ZOLTAN_OK)
       end if

    else

       min_num_nodes_after_balance = 0
       do while (min_num_nodes_after_balance == 0)
          
          ierr = Zoltan_LB_Balance(zz, changes, num_gid_entries, num_lid_entries, p1_num_import, p1_import_global_ids, &
             &    p1_import_local_ids, p1_import_procs, p1_num_export, p1_export_global_ids, p1_export_local_ids, p1_export_procs)
          assert(ierr == ZOLTAN_OK)
          
          ! calculate how many owned nodes we'd have after doing the planned load balancing
          num_nodes_after_balance = num_nodes + p1_num_import - p1_num_export
          
          ! find the minimum number of owned nodes any process would have after doing the planned load balancing
          call mpi_allreduce(num_nodes_after_balance, min_num_nodes_after_balance, 1, getPINTEGER(), &
             & MPI_MIN, MPI_COMM_FEMTOOLS, ierr)
          assert(ierr == MPI_SUCCESS)
          
          if (min_num_nodes_after_balance == 0) then
             ewrite(2,*) 'Empty partion would be created with load_imbalance_tolerance of', load_imbalance_tolerance
             load_imbalance_tolerance = 0.95 * load_imbalance_tolerance
             if (load_imbalance_tolerance < 1.075) then

                ewrite(1,*) 'Could not prevent empty partions by tightening load_imbalance_tolerance.'
                ewrite(1,*) 'Attempting to load balance with no edge-weights.'
                
                ! Reset the load_imbalance_tolerance
                ierr = Zoltan_Set_Param(zz, "IMBALANCE_TOL", "1.075"); assert(ierr == ZOLTAN_OK)
                ! Turn off the edge-weight calculation
                zoltan_global_calculate_edge_weights = .false.
                
                ierr = Zoltan_LB_Balance(zz, changes, num_gid_entries, num_lid_entries, p1_num_import, p1_import_global_ids, &
                   &    p1_import_local_ids, p1_import_procs, p1_num_export, p1_export_global_ids, p1_export_local_ids, p1_export_procs)
                assert(ierr == ZOLTAN_OK)
                
                ! calculate how many owned nodes we'd have after doing the planned load balancing
                num_nodes_after_balance = num_nodes + p1_num_import - p1_num_export
                
                ! find the minimum number of owned nodes any process would have after doing the planned load balancing
                call mpi_allreduce(num_nodes_after_balance, min_num_nodes_after_balance, 1, getPINTEGER(), &
                   & MPI_MIN, MPI_COMM_FEMTOOLS, ierr)
                assert(ierr == MPI_SUCCESS)
                
                if (min_num_nodes_after_balance == 0) then
                   FLAbort("Could not stop Zoltan creating empty partitions.")
                else
                   ewrite(-1,*) 'Load balancing was carried out without edge-weighting being applied. Mesh may not be of expected quality.'
                end if
             else
                ! convert load_imbalance_tolerance to a string for setting the option in Zoltan
                write(string_load_imbalance_tolerance, '(f6.3)' ) load_imbalance_tolerance
                ierr = Zoltan_Set_Param(zz, "IMBALANCE_TOL", string_load_imbalance_tolerance); assert(ierr == ZOLTAN_OK)
                
                ewrite(2,*) 'Tightened load_imbalance_tolerance to ', load_imbalance_tolerance
             end if
          end if
       end do
    end if
   
    do i=1,p1_num_export
       node = p1_export_local_ids(i)
       assert(node_owned(zoltan_global_zz_halo, node))
    end do

    ewrite(1,*) 'exiting zoltan_load_balance'

  end subroutine zoltan_load_balance

  subroutine derive_full_export_lists(states, p1_num_export, p1_export_local_ids, p1_export_procs, &
       & p1_num_export_full, p1_export_local_ids_full, p1_export_procs_full)
    type(state_type), dimension(:), intent(inout), target :: states
    integer(zoltan_int), intent(in) :: p1_num_export
    integer(zoltan_int), dimension(:), pointer, intent(in) :: p1_export_local_ids
    integer(zoltan_int), dimension(:), pointer, intent(in) :: p1_export_procs

    integer(zoltan_int), intent(out) :: p1_num_export_full    
    integer(zoltan_int), dimension(:), pointer, intent(out) :: p1_export_local_ids_full
    integer(zoltan_int), dimension(:), pointer, intent(out) :: p1_export_procs_full

    type(mesh_type), pointer :: full_mesh
    integer :: i, column
    integer, dimension(:), pointer :: column_nodes
    
    integer :: last_full_node
    logical :: lpreserve_columns
    
    ewrite(1,*) 'in derive_full_export_lists'
    
    p1_num_export_full = 0
    do i = 1, p1_num_export
       column = p1_export_local_ids(i)
       p1_num_export_full = p1_num_export_full + row_length(zoltan_global_columns_sparsity, column)
    end do
      
    allocate(p1_export_local_ids_full(p1_num_export_full))
    p1_export_local_ids_full = 0
    allocate(p1_export_procs_full(p1_num_export_full))
    p1_export_procs_full = -1
    
    last_full_node = 1
    do i = 1, p1_num_export
       column = p1_export_local_ids(i)
       column_nodes => row_m_ptr(zoltan_global_columns_sparsity, column)
       
       p1_export_local_ids_full(last_full_node:last_full_node+size(column_nodes)-1) = column_nodes
       p1_export_procs_full(last_full_node:last_full_node+size(column_nodes)-1) = p1_export_procs(i)
       
       last_full_node = last_full_node + size(column_nodes)
    end do
    assert(last_full_node-1==p1_num_export_full)
    assert(all(p1_export_local_ids_full>0))
    assert(all(p1_export_procs_full>-1))
    
    full_mesh => extract_mesh(states(1), trim(topology_mesh_name))
    lpreserve_columns = associated(full_mesh%columns)
    call allor(lpreserve_columns)
    if(lpreserve_columns) then
       allocate(zoltan_global_universal_columns(node_count(full_mesh)))
       zoltan_global_universal_columns = halo_universal_numbers(zoltan_global_zz_halo, full_mesh%columns)
    end if
    
    ewrite(1,*) 'exiting derive_full_export_lists'
    
  end subroutine derive_full_export_lists

  subroutine reset_zoltan_lists_full(zz, &
       & p1_num_export_full, p1_export_local_ids_full, p1_export_procs_full, &
       & p1_num_import, p1_import_global_ids, p1_import_local_ids, p1_import_procs, &
       & p1_num_export, p1_export_global_ids, p1_export_local_ids, p1_export_procs)

    type(zoltan_struct), pointer, intent(in) :: zz    

    integer(zoltan_int), intent(in) :: p1_num_export_full
    integer(zoltan_int), dimension(:), pointer, intent(in) :: p1_export_local_ids_full
    integer(zoltan_int), dimension(:), pointer, intent(in) :: p1_export_procs_full

    integer(zoltan_int), intent(out) :: p1_num_import
    integer(zoltan_int), dimension(:), pointer, intent(out) :: p1_import_global_ids
    integer(zoltan_int), dimension(:), pointer, intent(out) :: p1_import_local_ids
    integer(zoltan_int), dimension(:), pointer, intent(out) :: p1_import_procs
    integer(zoltan_int), intent(out) :: p1_num_export
    integer(zoltan_int), dimension(:), pointer, intent(out) :: p1_export_global_ids
    integer(zoltan_int), dimension(:), pointer, intent(out) :: p1_export_local_ids
    integer(zoltan_int), dimension(:), pointer, intent(out) :: p1_export_procs

    integer :: i, ierr
    
    ewrite(1,*) 'in reset_zoltan_lists_full'
    
    p1_num_export = p1_num_export_full
    
    p1_export_local_ids => p1_export_local_ids_full
    p1_export_procs => p1_export_procs_full
    
    allocate(p1_export_global_ids(p1_num_export))
    p1_export_global_ids = 0
    
    do i=1,p1_num_export
       p1_export_global_ids(i) = halo_universal_number(zoltan_global_zz_halo, p1_export_local_ids(i))
    end do
    assert(all(p1_export_global_ids>0))
    
    p1_num_import = 0
    p1_import_local_ids => null()
    p1_import_procs => null()
    p1_import_global_ids => null()
    
    ierr = Zoltan_Compute_Destinations(zz, &
         & p1_num_export, p1_export_global_ids, p1_export_local_ids, p1_export_procs, &
         & p1_num_import, p1_import_global_ids, p1_import_local_ids, p1_import_procs)
    assert(ierr == ZOLTAN_OK)
    
    zoltan_global_preserve_columns = associated(zoltan_global_zz_positions%mesh%columns)
    call allor(zoltan_global_preserve_columns)
    
    ewrite(1,*) 'exiting reset_zoltan_lists_full'
    
  end subroutine reset_zoltan_lists_full


  subroutine are_we_keeping_or_sending_nodes(p1_num_export, p1_export_local_ids, p1_export_procs)
    integer(zoltan_int), intent(in) :: p1_num_export
    integer(zoltan_int), dimension(:), pointer, intent(in) :: p1_export_local_ids
    integer(zoltan_int), dimension(:), pointer, intent(in) :: p1_export_procs

    integer :: node, i
    integer, dimension(halo_nowned_nodes(zoltan_global_zz_halo)) :: owned_nodes
    
    call allocate(zoltan_global_nodes_we_are_sending)
    call allocate(zoltan_global_nodes_we_are_keeping)
    
    do i=1,p1_num_export
       call insert(zoltan_global_nodes_we_are_sending, p1_export_local_ids(i), p1_export_procs(i))
    end do
    
    call get_owned_nodes(zoltan_global_zz_halo, owned_nodes)
    do i=1,size(owned_nodes)
       node = owned_nodes(i)
       if (.not. has_key(zoltan_global_nodes_we_are_sending, node)) then
          call insert(zoltan_global_nodes_we_are_keeping, node)
       end if
    end do
  end subroutine are_we_keeping_or_sending_nodes

  subroutine zoltan_migration_phase_one(zz, &
       & p1_num_import, p1_import_global_ids, p1_import_local_ids, p1_import_procs, &
       & p1_num_export, p1_export_global_ids, p1_export_local_ids, p1_export_procs)

    type(zoltan_struct), pointer, intent(in) :: zz    

    integer(zoltan_int), intent(in) :: p1_num_import
    integer(zoltan_int), dimension(:), pointer, intent(in) :: p1_import_global_ids
    integer(zoltan_int), dimension(:), pointer, intent(in) :: p1_import_local_ids
    integer(zoltan_int), dimension(:), pointer, intent(in) :: p1_import_procs
    integer(zoltan_int), intent(in) :: p1_num_export
    integer(zoltan_int), dimension(:), pointer, intent(in) :: p1_export_global_ids
    integer(zoltan_int), dimension(:), pointer, intent(in) :: p1_export_local_ids
    integer(zoltan_int), dimension(:), pointer, intent(in) :: p1_export_procs


    integer(zoltan_int) :: ierr
    integer(zoltan_int), dimension(:), pointer :: import_to_part, export_to_part
    
    ewrite(1,*) "In zoltan_migration_phase_one; objects to import: ", p1_num_import
    ewrite(1,*) "In zoltan_migration_phase_one; objects to export: ", p1_num_export
    import_to_part => null()
    export_to_part => null()

    ierr = Zoltan_Migrate(zz, p1_num_import, p1_import_global_ids, p1_import_local_ids, p1_import_procs, &
         & import_to_part, p1_num_export, p1_export_global_ids, p1_export_local_ids, p1_export_procs, export_to_part) 
    assert(ierr == ZOLTAN_OK)
  end subroutine zoltan_migration_phase_one

  subroutine deal_with_exporters
    ! The unpack routine is the ones that do most of the heavy lifting
    ! in reconstructing the new mesh, positions, etc.
    ! But they don't get called if a process is only an exporter!
    ! So we need to reconstruct the positions and nelist for the
    ! nodes we own, and mark the request for the halos in the same way.
    ! (The owner of one of our halo nodes might have changed, and we
    ! don't know that, or who owns it, so we have to ask the old owner,
    ! i.e. participate in the phase_two migration).
    integer :: i, j, old_local_number, new_local_number
    integer, dimension(:), pointer :: neighbours
    integer :: universal_number
    type(integer_set) :: halo_nodes_we_need_to_know_about
    type(integer_hash_table) :: universal_number_to_old_owner
    logical :: changed
    integer :: old_owner, new_owner
    type(mesh_type) :: new_mesh
    type(integer_set) :: halo_nodes_we_currently_own
    integer :: rank
    
    if (associated(zoltan_global_new_positions%refcount)) return
    ewrite(1,*) "In deal_with_exports"

    call allocate(zoltan_global_new_nodes)
    
    ! Need to allocate zoltan_global_new_nodes, zoltan_global_new_elements, zoltan_global_new_positions, zoltan_global_new_nelist, zoltan_global_universal_to_new_local_numbering
    do i=1,key_count(zoltan_global_nodes_we_are_keeping)
       old_local_number = fetch(zoltan_global_nodes_we_are_keeping, i)
       call insert(zoltan_global_new_nodes, halo_universal_number(zoltan_global_zz_halo, old_local_number), changed=changed)
    end do
    
    call allocate(halo_nodes_we_need_to_know_about)
    call allocate(halo_nodes_we_currently_own)
    call allocate(universal_number_to_old_owner)
    
    rank = getrank()
    
    do i=1,key_count(zoltan_global_nodes_we_are_keeping)
       old_local_number = fetch(zoltan_global_nodes_we_are_keeping, i)
       neighbours => row_m_ptr(zoltan_global_zz_sparsity_two, old_local_number)
       do j=1,size(neighbours)
          universal_number = halo_universal_number(zoltan_global_zz_halo, neighbours(j))
          call insert(zoltan_global_new_nodes, universal_number, changed=changed)
          if (changed) then
             old_owner = halo_node_owner(zoltan_global_zz_halo, neighbours(j)) - 1
             if (old_owner == rank) then
                call insert(halo_nodes_we_currently_own, neighbours(j))
             else
                call insert(halo_nodes_we_need_to_know_about, universal_number)
             end if
             call insert(universal_number_to_old_owner, universal_number, old_owner)
          end if
       end do
    end do

    ! We need to process these too -- nodes we own now, but will
    ! not own later, and will become our halo nodes.
    do i=1,key_count(halo_nodes_we_currently_own)
       universal_number = halo_universal_number(zoltan_global_zz_halo, fetch(halo_nodes_we_currently_own, i))
       call insert(zoltan_global_new_nodes, universal_number)
    end do

    call invert_set(zoltan_global_new_nodes, zoltan_global_universal_to_new_local_numbering)
    call allocate(new_mesh, key_count(zoltan_global_new_nodes), 0, zoltan_global_zz_mesh%shape, trim(zoltan_global_zz_mesh%name))
    new_mesh%option_path = zoltan_global_zz_mesh%option_path
    if(zoltan_global_preserve_columns) then
       allocate(new_mesh%columns(key_count(zoltan_global_new_nodes)))
    end if
    call allocate(zoltan_global_new_positions, zoltan_global_zz_positions%dim, new_mesh, trim(zoltan_global_zz_positions%name))
    zoltan_global_new_positions%option_path = zoltan_global_zz_positions%option_path
    call deallocate(new_mesh)
    allocate(zoltan_global_new_snelist(key_count(zoltan_global_new_nodes)))
    call allocate(zoltan_global_new_snelist)
    call allocate(zoltan_global_new_surface_elements)
    
    allocate(zoltan_global_new_nelist(key_count(zoltan_global_new_nodes)))
    do i=1,key_count(zoltan_global_new_nodes)
       call allocate(zoltan_global_new_nelist(i))
    end do
    call allocate(zoltan_global_new_elements)
    
    do i=1,key_count(zoltan_global_nodes_we_are_keeping)
       old_local_number = fetch(zoltan_global_nodes_we_are_keeping, i)
       universal_number = halo_universal_number(zoltan_global_zz_halo, old_local_number)
       new_local_number = fetch(zoltan_global_universal_to_new_local_numbering, universal_number)
       call set(zoltan_global_new_positions, new_local_number, node_val(zoltan_global_zz_positions, old_local_number))
       
       if(zoltan_global_preserve_columns) then
          zoltan_global_new_positions%mesh%columns(new_local_number) = fetch(zoltan_global_universal_to_new_local_numbering_m1d, &
               zoltan_global_universal_columns(old_local_number))
       end if
       
       neighbours => row_m_ptr(zoltan_global_zz_nelist, old_local_number)
       do j=1,size(neighbours)
          call insert(zoltan_global_new_nelist(new_local_number), halo_universal_number(zoltan_global_zz_ele_halo, neighbours(j)))
          call insert(zoltan_global_new_elements, halo_universal_number(zoltan_global_zz_ele_halo, neighbours(j)))
          ! don't need to add anything to zoltan_global_universal_element_number_to_region_id because we already have it
       end do

       ! and record the snelist information
       do j=1,key_count(zoltan_global_old_snelist(old_local_number))
          call insert(zoltan_global_new_snelist(new_local_number), fetch(zoltan_global_old_snelist(old_local_number), j))
          call insert(zoltan_global_new_surface_elements, fetch(zoltan_global_old_snelist(old_local_number), j))
          ! we don't need to add anything to the zoltan_global_universal_surface_number_to_surface_id because we already have it
       end do
    end do

    ! Set the positions and nelist of halo_nodes_we_currently_own
    do i=1,key_count(halo_nodes_we_currently_own)
       old_local_number = fetch(halo_nodes_we_currently_own, i)
       universal_number = halo_universal_number(zoltan_global_zz_halo, old_local_number)
       new_local_number = fetch(zoltan_global_universal_to_new_local_numbering, universal_number)
       call set(zoltan_global_new_positions, new_local_number, node_val(zoltan_global_zz_positions, old_local_number))
       
       if(zoltan_global_preserve_columns) then
          zoltan_global_new_positions%mesh%columns(new_local_number) = fetch(zoltan_global_universal_to_new_local_numbering_m1d, &
               zoltan_global_universal_columns(old_local_number))
       end if

       neighbours => row_m_ptr(zoltan_global_zz_nelist, old_local_number)
       do j=1,size(neighbours)
          call insert(zoltan_global_new_nelist(new_local_number), halo_universal_number(zoltan_global_zz_ele_halo, neighbours(j)))
          call insert(zoltan_global_new_elements, halo_universal_number(zoltan_global_zz_ele_halo, neighbours(j)))
          ! don't need to add anything to zoltan_global_universal_element_number_to_region_id because we already have it
       end do

       new_owner = fetch(zoltan_global_nodes_we_are_sending, old_local_number)
       call insert(zoltan_global_receives(new_owner+1), universal_number)
       
       ! and record the snelist information
       do j=1,key_count(zoltan_global_old_snelist(old_local_number))
          call insert(zoltan_global_new_snelist(new_local_number), fetch(zoltan_global_old_snelist(old_local_number), j))
          call insert(zoltan_global_new_surface_elements, fetch(zoltan_global_old_snelist(old_local_number), j))
       end do
    end do
    call deallocate(halo_nodes_we_currently_own)
    
    zoltan_global_my_num_import = key_count(halo_nodes_we_need_to_know_about)
    allocate(zoltan_global_my_import_procs(zoltan_global_my_num_import))
    allocate(zoltan_global_my_import_global_ids(zoltan_global_my_num_import))
    do i=1,zoltan_global_my_num_import
       universal_number = fetch(halo_nodes_we_need_to_know_about, i)
       zoltan_global_my_import_global_ids(i) = universal_number
       zoltan_global_my_import_procs(i) = fetch(universal_number_to_old_owner, universal_number)
       assert(zoltan_global_my_import_procs(i) /= getrank())
    end do

    call deallocate(halo_nodes_we_need_to_know_about)
    call deallocate(universal_number_to_old_owner)

  end subroutine deal_with_exporters

  subroutine zoltan_migration_phase_two(zz)
    type(zoltan_struct), pointer, intent(in) :: zz    

    integer(zoltan_int) :: ierr
    integer(zoltan_int), dimension(:), pointer :: import_to_part, export_to_part, export_procs
    integer(zoltan_int), dimension(:), pointer :: export_global_ids, export_local_ids, import_local_ids
    integer(zoltan_int) :: num_export

    ! Register the new callback functions for packing and unpacking
    ierr = Zoltan_Set_Fn(zz, ZOLTAN_OBJ_SIZE_MULTI_FN_TYPE, zoltan_cb_pack_halo_node_sizes); assert(ierr == ZOLTAN_OK)
    ierr = Zoltan_Set_Fn(zz, ZOLTAN_PACK_OBJ_MULTI_FN_TYPE, zoltan_cb_pack_halo_nodes); assert(ierr == ZOLTAN_OK)
    ierr = Zoltan_Set_Fn(zz, ZOLTAN_UNPACK_OBJ_MULTI_FN_TYPE, zoltan_cb_unpack_halo_nodes); assert(ierr == ZOLTAN_OK)

    import_to_part => null()
    export_to_part => null()
    export_global_ids => null()
    export_local_ids => null()
    export_procs => null()
    num_export = -1

    if (.not. associated(zoltan_global_my_import_procs)) then
       ! We have nothing else to receive, but we still need to take part in
       ! the communication. Set num_import to 0
       zoltan_global_my_num_import = 0
       allocate(zoltan_global_my_import_global_ids(0))
       allocate(zoltan_global_my_import_procs(0))
    end if
    ewrite(1,*) "In zoltan_migration_phase_two; objects to import: ", zoltan_global_my_num_import
    
    ! We should be able to do:
    ! import_local_ids => null()
    ! from my reading of the Zoltan docs. But it actually doesn't appear to be the case.
    allocate(import_local_ids(zoltan_global_my_num_import))
    if (zoltan_global_my_num_import > 0) then
       import_local_ids = 666
    end if

    assert(associated(zoltan_global_my_import_global_ids))
    assert(associated(import_local_ids))
    assert(associated(zoltan_global_my_import_procs))
    assert(all(zoltan_global_my_import_procs >= 0))
    assert(all(zoltan_global_my_import_procs < getnprocs()))
    ierr = Zoltan_Compute_Destinations(zz, zoltan_global_my_num_import, zoltan_global_my_import_global_ids, import_local_ids, zoltan_global_my_import_procs, &
         & num_export, export_global_ids, export_local_ids, export_procs)
    assert(ierr == ZOLTAN_OK)
    ewrite(1,*) "In zoltan_migration_phase_two; objects to export: ", num_export

    ierr = Zoltan_Migrate(zz, zoltan_global_my_num_import, zoltan_global_my_import_global_ids, import_local_ids, zoltan_global_my_import_procs, &
         & import_to_part, num_export, export_global_ids, export_local_ids, export_procs, export_to_part) 
    assert(ierr == ZOLTAN_OK)

    ierr = Zoltan_LB_Free_Part(export_global_ids, export_local_ids, export_procs, export_to_part)
    assert(ierr == ZOLTAN_OK)
    
    deallocate(import_local_ids)
  end subroutine zoltan_migration_phase_two

  subroutine deallocate_my_lists
    deallocate(zoltan_global_my_import_global_ids)
    deallocate(zoltan_global_my_import_procs)
  end subroutine deallocate_my_lists

  subroutine reconstruct_enlist
    type(csr_sparsity):: eelist
    type(mesh_type):: temporary_mesh
    type(integer_set), dimension(key_count(zoltan_global_new_elements)) :: enlists
    integer :: i, j, k, expected_loc, full_elements, connected_elements
    integer :: universal_number, new_local_number
    type(integer_set) :: new_elements_we_actually_have
    integer, dimension(:), pointer:: neigh
    
    ewrite(1,*) "In reconstruct_enlist"
    
    ! zoltan_global_new_elements currently contains the universal numbers of elements
    ! we don't fully have and won't be in the final mesh.
    ! So uen_to_new_local_numbering here is just temporary -- we will
    ! construct the proper version later.
    call invert_set(zoltan_global_new_elements, zoltan_global_uen_to_new_local_numbering)
    
    do i=1,key_count(zoltan_global_new_elements)
       call allocate(enlists(i))
    end do
    
    ! Invert the nelists to give the enlists
    
    do i=1,key_count(zoltan_global_new_nodes)
       do j=1,key_count(zoltan_global_new_nelist(i))
          universal_number = fetch(zoltan_global_new_nelist(i), j)
          new_local_number = fetch(zoltan_global_uen_to_new_local_numbering, universal_number)
          call insert(enlists(new_local_number), i)
       end do
    end do

    call deallocate(zoltan_global_uen_to_new_local_numbering)
    
    ! Now, some of these will be degenerate, because the halo nodes will refer
    ! to elements we don't know about. We can tell these apart because they
    ! are incomplete.
    
    full_elements = 0
    ! For mixed meshes, this should be the loc of the positions mesh
    ! note we know the universal number -- fetch(zoltan_global_new_elements, i)
    ! However, it's constant for now
    expected_loc = zoltan_global_zz_mesh%shape%loc
    
    ! First, count how many we have
    do i=1,key_count(zoltan_global_new_elements)
       if (key_count(enlists(i)) == expected_loc) then
          full_elements = full_elements + 1
          !        else
          !          write(0,*) "Element ", fetch(zoltan_global_new_elements, i), " is degenerate. Dropping .."
       end if
    end do

    ewrite(2,*) "Found ", key_count(zoltan_global_new_elements), " possible new elements."
    ewrite(2,*) "Of these, ", full_elements, " are non-degenerate."

    ! Now we construct a temporary mesh of full elements
    ! This mesh is temporary because we also want to drop elements that are not connected
    ! to any other elements
      
    call allocate(temporary_mesh, node_count(zoltan_global_new_positions), full_elements, zoltan_global_new_positions%mesh%shape, &
         name="TemporaryZoltanMesh")
    j = 1
    do i=1,key_count(zoltan_global_new_elements)
       if (key_count(enlists(i)) == expected_loc) then
          call set_ele_nodes(temporary_mesh, j, set2vector(enlists(i)))
          j = j + 1
       end if
    end do
    
    call add_nelist(temporary_mesh)
    call extract_lists(temporary_mesh, eelist=eelist)
    
    connected_elements=0
    do i=1, element_count(temporary_mesh)
       neigh => row_m_ptr(eelist, i)
       if (any(neigh>0)) connected_elements = connected_elements + 1
    end do
    
    ewrite(2,*) "Of the ", full_elements, " full elements, ", connected_elements, " are connected."
    
    call allocate(new_elements_we_actually_have)
    call allocate(zoltan_global_uen_to_new_local_numbering)
    zoltan_global_new_positions%mesh%elements = connected_elements
    deallocate(zoltan_global_new_positions%mesh%ndglno)
    allocate(zoltan_global_new_positions%mesh%ndglno(connected_elements * expected_loc))
#ifdef HAVE_MEMORY_STATS
    call register_allocation("mesh_type", "integer", connected_elements * expected_loc, name=zoltan_global_new_positions%mesh%name)
#endif
    if(zoltan_global_preserve_mesh_regions) then
       allocate(zoltan_global_new_positions%mesh%region_ids(connected_elements))
    end if
      
    j = 1 ! index connected full elements (new local element numbering)
    k = 1 ! indexes full elements
    do i=1, key_count(zoltan_global_new_elements)
       if (key_count(enlists(i)) == expected_loc) then
          ! only for full elements
          neigh => row_m_ptr(eelist, k)
          if (any(neigh>0)) then
             ! of these only the connected elements
             universal_number = fetch(zoltan_global_new_elements, i)
             call set_ele_nodes(zoltan_global_new_positions%mesh, j, set2vector(enlists(i)))
             call insert(new_elements_we_actually_have, universal_number)
             call insert(zoltan_global_uen_to_new_local_numbering, universal_number, j)
             if(zoltan_global_preserve_mesh_regions) then
                zoltan_global_new_positions%mesh%region_ids(j) = fetch(zoltan_global_universal_element_number_to_region_id, universal_number)
             end if
             j = j + 1
          end if
          k = k + 1
       end if
    end do
        
    assert( k==full_elements+1 )
    assert( j==connected_elements+1 )

    do i=1,size(enlists)
       call deallocate(enlists(i))
    end do
        
    call deallocate(temporary_mesh)
    call deallocate(zoltan_global_new_nelist)
    deallocate(zoltan_global_new_nelist)
      
    ! New elements is no longer valid, as we have lost the degenerate elements
    call deallocate(zoltan_global_new_elements)
    zoltan_global_new_elements = new_elements_we_actually_have
    
    call deallocate(zoltan_global_universal_element_number_to_region_id)
    
    ! Bingo! Our mesh has an enlist.
    ewrite(1,*) "Exiting reconstruct_enlist"
    
  end subroutine reconstruct_enlist
    
  subroutine reconstruct_senlist
    type(integer_set), dimension(key_count(zoltan_global_new_surface_elements)) :: senlists
    integer :: i, j, expected_loc, full_elements
    integer :: universal_number, new_local_number
    type(integer_hash_table) :: universal_surface_element_to_local_numbering
    integer, dimension(:), allocatable, target :: surface_ids, element_owners
    type(csr_sparsity), pointer :: nnlist
    
    logical, dimension(key_count(zoltan_global_new_surface_elements)) :: keep_surface_element
    integer, dimension(:), allocatable :: sndgln
    integer :: universal_element_number
    
    ewrite(1,*) "In reconstruct_senlist"
    
    ! zoltan_global_new_surface_elements currently contains the universal numbers of surface elements
    ! we don't fully have and won't be in the final mesh.
    ! So universal_surface_element_to_local_numbering here is just temporary.
    call invert_set(zoltan_global_new_surface_elements, universal_surface_element_to_local_numbering)
    
    do i=1,key_count(zoltan_global_new_surface_elements)
       call allocate(senlists(i))
    end do
    
    ! Invert the snelists to give the senlists
    
    do i=1,key_count(zoltan_global_new_nodes)
       do j=1,key_count(zoltan_global_new_snelist(i))
          universal_number = fetch(zoltan_global_new_snelist(i), j)
          new_local_number = fetch(universal_surface_element_to_local_numbering, universal_number)
          call insert(senlists(new_local_number), i)
       end do
    end do
    
    call deallocate(universal_surface_element_to_local_numbering)
    
    ! Now, some of these will be degenerate, because the halo nodes will refer
    ! to elements we don't know about. We can tell these apart because they
    ! are incomplete.
    
    full_elements = 0
    ! For mixed meshes, this should taken from the loc of the positions mesh
    ! However, it's constant for now
    expected_loc = face_loc(zoltan_global_zz_mesh, 1)
    
    ! First, count how many we have
    nnlist => extract_nnlist(zoltan_global_new_positions%mesh)
    do i=1,key_count(zoltan_global_new_surface_elements)
       j = key_count(senlists(i))
       assert(j <= expected_loc)
       if (j == expected_loc) then
          ! We also need to check if we have the parent volume element --
          ! it's possible to get all the information for a face, without having the
          ! corresponding element!
          universal_number = fetch(zoltan_global_new_surface_elements, i)
          universal_element_number = fetch(zoltan_global_universal_surface_number_to_element_owner, universal_number)
          keep_surface_element(i) = has_key(zoltan_global_uen_to_new_local_numbering, universal_element_number)
          if (keep_surface_element(i)) full_elements = full_elements + 1
       else
          keep_surface_element(i) = .false.
          ! write(0,*) "Surface element ", fetch(zoltan_global_new_surface_elements, i), " is degenerate. Dropping .."
          ! write(0,*) "Local nodes: ", set2vector(senlists(i))
       end if
    end do

    ewrite(2,*) "Found ", key_count(zoltan_global_new_surface_elements), " possible new surface elements."
    ewrite(2,*) "Of these, ", full_elements, " are non-degenerate."

    ! And now fill in the non-degenerate ones
    
    allocate(sndgln(full_elements * expected_loc))
    allocate(surface_ids(full_elements))
    allocate(element_owners(full_elements))
    
    j = 1
    do i=1,key_count(zoltan_global_new_surface_elements)
       if (keep_surface_element(i)) then
          universal_number = fetch(zoltan_global_new_surface_elements, i)
          sndgln( (j-1)*expected_loc+1 : j*expected_loc ) = set2vector(senlists(i))
          surface_ids(j) = fetch(zoltan_global_universal_surface_number_to_surface_id, universal_number)
          universal_element_number = fetch(zoltan_global_universal_surface_number_to_element_owner, universal_number)
          element_owners(j) = fetch(zoltan_global_uen_to_new_local_numbering, universal_element_number)
          j = j + 1
       end if
    end do
    assert(j == full_elements + 1)

    if (zoltan_global_zz_mesh%faces%has_discontinuous_internal_boundaries) then
      ! for internal facet pairs, the surface ids are not necessarily the same (this is used in periodic meshes)
      ! we need to tell add_faces which facets is on which side by supplying element ownership info
      call add_faces(zoltan_global_new_positions%mesh, sndgln=sndgln, boundary_ids=surface_ids, element_owner=element_owners)
    else
      ! surface ids on facets pairs are assumed consistent - add_faces will copy the first of the pair it encounters
      ! in sndgln on either side - the next copy in sndgln is ignored (only checked that its surface id is consistent)
      call add_faces(zoltan_global_new_positions%mesh, sndgln=sndgln, boundary_ids=surface_ids, &
        allow_duplicate_internal_facets=.true.)
    end if
    
    do i=1,size(senlists)
       call deallocate(senlists(i))
    end do
    
    ! New elements is no longer valid, as we have lost the degenerate elements
    call deallocate(zoltan_global_new_surface_elements)
    call deallocate(zoltan_global_universal_surface_number_to_surface_id)
    call deallocate(zoltan_global_universal_surface_number_to_element_owner)
    call deallocate(zoltan_global_new_snelist)
    deallocate(zoltan_global_new_snelist)
    
    deallocate(sndgln)
    deallocate(surface_ids)
    deallocate(element_owners)
    
    call deinterleave_surface_ids(zoltan_global_new_positions%mesh, max_coplanar_id)
    
    ! Bingo! Our mesh has an senlist.
    ewrite(1,*) "Exiting reconstruct_senlist"
  end subroutine reconstruct_senlist

  subroutine reconstruct_halo(zz)
    ! At this point, the receives sets have been populated with all
    ! the universal node numbers we need to receive from each process.
    ! So, we are going to use zoltan to invert this to compute
    ! the send list for each process too.
    ! Then we will allocate the l2n halo and set it.
    ! Then we will chop it down to form the l1n halo, the l1e halo, and the
    ! l2e halo. 
    ! Supply the peeps with jeeps, brick apiece, capiche?

    type(zoltan_struct), pointer, intent(in) :: zz    

    integer :: num_import, num_export
    integer, dimension(:), pointer :: import_global_ids, import_local_ids, import_procs
    integer, dimension(:), pointer :: export_global_ids, export_local_ids, export_procs, export_to_part
    integer :: ierr, i, head
    type(integer_set), dimension(size(zoltan_global_receives)) :: sends
    integer, dimension(size(zoltan_global_receives)) :: nreceives, nsends
    integer, dimension(ele_count(zoltan_global_new_positions)) :: ele_renumber_permutation
    integer, dimension(node_count(zoltan_global_new_positions)) :: node_renumber_permutation
    integer :: universal_element_number, old_new_local_element_number, new_new_local_element_number
    integer :: universal_node_number, old_new_local_node_number, new_new_local_node_number
    
    integer, dimension(ele_count(zoltan_global_new_positions)) :: old_new_region_ids
    
    ewrite(1,*) "In reconstruct_halo"
    
    num_import = 0
    do i=1,size(zoltan_global_receives)
       nreceives(i) = key_count(zoltan_global_receives(i))
       num_import = num_import + nreceives(i)
    end do
    
    allocate(import_global_ids(num_import))
    allocate(import_local_ids(num_import))
    allocate(import_procs(num_import))
    
    import_local_ids = 666
    head = 1
    do i=1,size(zoltan_global_receives)
       import_global_ids(head:head + nreceives(i) - 1) = set2vector(zoltan_global_receives(i))
       import_procs(head:head + nreceives(i) - 1) = i - 1
       head = head + nreceives(i)
    end do
    
    export_global_ids => null()
    export_local_ids => null()
    export_procs => null()
    export_to_part => null()
    
    ierr = Zoltan_Compute_Destinations(zz, &
         & num_import, import_global_ids, import_local_ids, import_procs, &
         & num_export, export_global_ids, export_local_ids, export_procs)
    assert(ierr == ZOLTAN_OK)
    
    ! Now we know the sends too! Thanks, Zoltan!
    
    deallocate(import_global_ids)
    deallocate(import_local_ids)
    deallocate(import_procs)
    
    ! Now create the sends sets .. easier than pulling it straight out of Zoltan's data structures
    ! as Zoltan does NOT explicitly guarantee that the sends are organised such that export_procs looks
    ! like 
    ! [ sends to process 0 | sends to process 1 | sends to process 2 ... ]
    ! If such a guarantee were available, then it would be just as easy to use that, but there
    ! is no such guarantee given in the documentation ...
    ! Is an appetite for destruction, slap a murder rap on this production
    
    do i=1,size(sends)
       call allocate(sends(i))
    end do
    
    do i=1,num_export
       call insert(sends(export_procs(i)+1), export_global_ids(i))
    end do
    
    do i=1,size(sends)
       nsends(i) = key_count(sends(i))
    end do
    
    ! Allocate the halo and such
    ! We had to grow dreads to change our description, two cops is on a milkbox, missin'
    
    allocate(zoltan_global_new_positions%mesh%halos(2))
    call allocate(zoltan_global_new_positions%mesh%halos(2), &
         nsends = nsends, &
         nreceives = nreceives, &
         name = halo_name(zoltan_global_zz_halo), &
         communicator = halo_communicator(zoltan_global_zz_halo), &
         nowned_nodes = key_count(zoltan_global_new_nodes) - num_import, &
         data_type = halo_data_type(zoltan_global_zz_halo))

    do i=1,size(zoltan_global_receives)
       call set_halo_sends(zoltan_global_new_positions%mesh%halos(2), i, fetch(zoltan_global_universal_to_new_local_numbering, set2vector(sends(i))))
       call set_halo_receives(zoltan_global_new_positions%mesh%halos(2), i, fetch(zoltan_global_universal_to_new_local_numbering, set2vector(zoltan_global_receives(i))))
    end do

    ! Now derive all the other halos ...
    ! And me teevee's always off, cause I see something truly black then
    
    call derive_l1_from_l2_halo(zoltan_global_new_positions%mesh, ordering_scheme = HALO_ORDER_GENERAL, create_caches = .false.)
    call renumber_positions_trailing_receives(zoltan_global_new_positions, permutation=node_renumber_permutation)
    assert(has_ownership(zoltan_global_new_positions%mesh%halos(2)))
    assert(has_ownership(zoltan_global_new_positions%mesh%halos(1)))
    allocate(zoltan_global_new_positions%mesh%element_halos(2))
    call derive_element_halo_from_node_halo(zoltan_global_new_positions%mesh, create_caches = .false.)
    call renumber_positions_elements_trailing_receives(zoltan_global_new_positions, permutation=ele_renumber_permutation)
    assert(has_ownership(zoltan_global_new_positions%mesh%halos(2)))
    assert(has_ownership(zoltan_global_new_positions%mesh%halos(1)))
    
    if(zoltan_global_preserve_mesh_regions) then
       old_new_region_ids = zoltan_global_new_positions%mesh%region_ids
    end if
    
    ! The previous routine has renumbered all the local elements to put the element halos in
    ! trailing receives status. However, we need the universal element number -> local element number
    ! later in the field transfer. So we need to update our records now.
    do i=1,ele_count(zoltan_global_new_positions)
       universal_element_number = fetch(zoltan_global_new_elements, i)
       old_new_local_element_number = i
       new_new_local_element_number = ele_renumber_permutation(old_new_local_element_number)
       call insert(zoltan_global_uen_to_new_local_numbering, universal_element_number, new_new_local_element_number)
       if(zoltan_global_preserve_mesh_regions) then
          zoltan_global_new_positions%mesh%region_ids(new_new_local_element_number) = old_new_region_ids(old_new_local_element_number)
       end if
    end do
    ! We're also going to need the universal to new local numbering for 2+1d adaptivity
    do i=1,node_count(zoltan_global_new_positions)
       universal_node_number = fetch(zoltan_global_new_nodes, i)
       old_new_local_node_number = i
       new_new_local_node_number = node_renumber_permutation(old_new_local_node_number)
       call insert(zoltan_global_universal_to_new_local_numbering, universal_node_number, new_new_local_node_number)
    end do
    
    do i=1,halo_count(zoltan_global_new_positions)
       assert(halo_verifies(zoltan_global_new_positions%mesh%halos(i), zoltan_global_new_positions))
    end do
    
    ! Now cleanup
    
    ierr = Zoltan_LB_Free_Part(export_global_ids, export_local_ids, export_procs, export_to_part)
    assert(ierr == ZOLTAN_OK)

    do i=1,size(sends)
       call deallocate(sends(i))
    end do
      
    call reorder_element_numbering(zoltan_global_new_positions)
    
    ewrite(1,*) "Exiting reconstruct_halo"
    
  end subroutine reconstruct_halo
    
  subroutine initialise_transfer(zz, states, zoltan_global_new_positions_m1d, metric, full_metric, new_metric, initialise_fields, skip_extrusion_after)
    type(zoltan_struct), pointer, intent(in) :: zz
    type(state_type), dimension(:), intent(inout), target :: states
    type(vector_field), intent(inout) :: zoltan_global_new_positions_m1d
    type(tensor_field), intent(inout), optional :: metric
    type(tensor_field), intent(inout), optional :: full_metric
    type(tensor_field), intent(out) :: new_metric
    logical, intent(in), optional :: initialise_fields
    logical, intent(in), optional :: skip_extrusion_after

    integer :: i
    type(state_type), dimension(size(states)) :: interpolate_states
    integer(zoltan_int) :: ierr
    character(len=FIELD_NAME_LEN), dimension(:), allocatable :: mesh_names
    type(mesh_type), pointer :: mesh
    integer :: no_meshes
    
    ewrite(1,*) 'in initialise_transfer'
    
    ! Set up zoltan_global_source_states
    do i=1,size(states)
       call select_fields_to_interpolate(states(i), interpolate_states(i), no_positions=.true., &
            first_time_step=initialise_fields)
       ! Remove the current state as we've copied the bits we need
       call deallocate(states(i))
    end do
    
    ! Interpolate the metric, too
    if (present(full_metric)) then
       call insert(interpolate_states(1), full_metric, "ErrorMetric")
       call deallocate(full_metric)
    else if (present(metric)) then
       call insert(interpolate_states(1), metric, "ErrorMetric")
       call deallocate(metric)
    end if
    
    allocate( mesh_names(1:mesh_count(interpolate_states(1))) )
    no_meshes = 0
    do i=1, mesh_count(interpolate_states(1))
       mesh => extract_mesh(interpolate_states(1), i)
       if (zoltan_global_migrate_extruded_mesh .and. mesh_dim(mesh)/=mesh_dim(zoltan_global_new_positions)) cycle
       no_meshes = no_meshes + 1
       mesh_names(no_meshes) = mesh%name
    end do
    
    allocate(zoltan_global_source_states(no_meshes))
    call halo_update(interpolate_states, level=1)
    ! Place the fields we've picked out to interpolate onto the correct meshes of zoltan_global_source_states
    call collect_fields_by_mesh(interpolate_states, mesh_names(1:no_meshes), zoltan_global_source_states)
    
    ! Finished with interpolate_states for setting up zoltan_global_source_states
    do i=1,size(interpolate_states)
       call deallocate(interpolate_states(i))
    end do
    
    if (mesh_periodic(zoltan_global_zz_mesh)) then
       zoltan_global_new_positions%mesh%periodic = .true.
    end if
    
    ! Start setting up states so that it can be populated with migrated fields data
    
    ! Put the new positions mesh into states
    
    if(zoltan_global_migrate_extruded_mesh) then
       if (mesh_periodic(zoltan_global_zz_mesh)) then
          zoltan_global_new_positions_m1d%mesh%periodic = .true.
       end if
       call insert(states, zoltan_global_new_positions_m1d%mesh, name = zoltan_global_new_positions_m1d%mesh%name)
       call insert(states, zoltan_global_new_positions_m1d, name = zoltan_global_new_positions_m1d%name)
    end if
    
    call insert(states, zoltan_global_new_positions%mesh, name = zoltan_global_new_positions%mesh%name)
    call insert(states, zoltan_global_new_positions, name = zoltan_global_new_positions%name)

    ! Check the number of halos in our new mesh
    ! Used in a few places within the Zoltan callback and detector routines
    zoltan_global_new_positions_mesh_nhalos = halo_count(zoltan_global_new_positions%mesh)
    assert(zoltan_global_new_positions_mesh_nhalos == 2)
    
    ! Allocate a new metric field on the new positions mesh and zero it
    if (present(metric).or.present(full_metric)) then
       call allocate(new_metric, zoltan_global_new_positions%mesh, "ErrorMetric")
       call set(new_metric,spread(spread(666.0, 1, new_metric%dim(1)), 2, new_metric%dim(2)))
    end if

    ! Setup meshes and fields on states
    call restore_reserved_meshes(states)
    call insert_derived_meshes(states, skip_extrusion=skip_extrusion_after)
    call allocate_and_insert_fields(states)
    call restore_reserved_fields(states)
    
    ! And set up zoltan_global_target_states based on states
    do i=1,size(states)
       call select_fields_to_interpolate(states(i), interpolate_states(i), no_positions=.true., &
            first_time_step=initialise_fields)
    end do

    ! Metric will be interpolated too, it is 666.0 (for debugging purposes) at this point
    if (present(metric).or.present(full_metric)) then
       call insert(interpolate_states(1), new_metric, "ErrorMetric")
    end if
    
    allocate(zoltan_global_target_states(no_meshes))
    call collect_fields_by_mesh(interpolate_states, mesh_names(1:no_meshes), zoltan_global_target_states)
      
    ! Finished with interpolate states for setting up zoltan_global_target_states
    do i=1,size(interpolate_states)
       call deallocate(interpolate_states(i))
    end do

    ! Tell Zoltan which callback functions to use for the migration
    ierr = Zoltan_Set_Fn(zz, ZOLTAN_OBJ_SIZE_MULTI_FN_TYPE, zoltan_cb_pack_field_sizes); assert(ierr == ZOLTAN_OK)
    ierr = Zoltan_Set_Fn(zz, ZOLTAN_PACK_OBJ_MULTI_FN_TYPE, zoltan_cb_pack_fields); assert(ierr == ZOLTAN_OK)
    ierr = Zoltan_Set_Fn(zz, ZOLTAN_UNPACK_OBJ_MULTI_FN_TYPE, zoltan_cb_unpack_fields); assert(ierr == ZOLTAN_OK)
    
    ewrite(1,*) 'exiting initialise_transfer'
    
  end subroutine initialise_transfer

  subroutine update_detector_list_element(detector_list_array)
    ! Update the detector%element field for every detector left in our list
    ! and check that we did not miss any in the first send
    ! broadcast them if we did
    type(detector_list_ptr), dimension(:), intent(inout) :: detector_list_array

    type(detector_linked_list), pointer :: detector_list => null()
    type(detector_linked_list) ::  detector_send_list
    type(detector_type), pointer :: detector => null(), send_detector => null()
    integer :: i, j, send_count, ierr
    integer :: old_local_element_number, new_local_element_number, old_universal_element_number
    integer, allocatable :: ndets_being_sent(:)
    real, allocatable :: send_buff(:,:), recv_buff(:,:)
    logical do_broadcast
    type(element_type), pointer :: shape

    ewrite(1,*) "In update_detector_list_element"

    send_count=0

    !Loop for detectors or particles with no attributes

    do j = 1, size(detector_list_array)
       detector_list => detector_list_array(j)%ptr
       ewrite(2,*) "Length of detector list to be updated: ", detector_list%length

       detector => detector_list%first
       !Cycle if particle has attributes
       if (associated(detector)) then
          if (size(detector%attributes)>=1) then
             cycle
          end if
       end if
       do while (associated(detector))

          old_local_element_number = detector%element

          if (.not. has_key(zoltan_global_old_local_numbering_to_uen, old_local_element_number)) then
             ewrite(-1,*) "Zoltan can't find old element number for detector ", detector%id_number
             FLAbort('Trying to update unknown detector in Zoltan')
          end if
          old_universal_element_number = fetch(zoltan_global_old_local_numbering_to_uen, old_local_element_number)

          if(has_key(zoltan_global_uen_to_new_local_numbering, old_universal_element_number)) then
             ! Update the element number for the detector
             detector%element = fetch(zoltan_global_uen_to_new_local_numbering, old_universal_element_number)
             detector => detector%next
          else
             ! We no longer own the element containing this detector, and cannot establish its new
             ! owner from the halo, because the boundary has moved too far.
             ! Since we have no way of determining the new owner we are going to broadcast the detector
             ! to all procs, so move it to the send list and count how many detectors we're sending.
             ewrite(2,*) "Found non-local detector, initialising broadcast..."
             send_count = send_count + 1
          
             ! Store the old universal element number for unpacking to new local at the receive
             detector%element = old_universal_element_number
             detector%list_id=detector_list%id

             ! Remove detector from detector list
             send_detector => detector
             detector => detector%next
             call move(send_detector, detector_list, detector_send_list)
          end if
       end do
    end do

    ! Find out how many detectors each process wants to broadcast
    allocate(ndets_being_sent(getnprocs()))
    call mpi_allgather(send_count, 1, getPINTEGER(), ndets_being_sent, 1 , getPINTEGER(), MPI_COMM_FEMTOOLS, ierr)
    assert(ierr == MPI_SUCCESS)

    ! Check whether we have to perform broadcast, if not return
    do_broadcast=.false.
    if (any(ndets_being_sent > 0)) then
       do_broadcast=.true.
    end if
    if (.not. do_broadcast) then
       ewrite(1,*) "Exiting update_detector_list_element"
       return
    end if
    ewrite(2,*) "Broadcast required, initialising..."

    ! Allocate memory for all the detectors you're going to send
    allocate(send_buff(send_count,zoltan_global_ndata_per_det))
      
    detector => detector_send_list%first
    do i=1,send_count
       ! Pack the detector information and delete from send_list (delete advances detector to detector%next)
       call pack_detector(detector, send_buff(i, 1:zoltan_global_ndata_per_det), zoltan_global_ndims)
       call delete(detector, detector_send_list)
    end do

    ! Broadcast detectors whose new owner we can't identify
    do i=1,getnprocs()
       if (ndets_being_sent(i) > 0) then

          if (i == getprocno()) then
             ! Broadcast the detectors you want to send
             ewrite(2,*) "Broadcasting ", send_count, " detectors"
             call mpi_bcast(send_buff,send_count*zoltan_global_ndata_per_det, getPREAL(), i-1, MPI_COMM_FEMTOOLS, ierr)
             assert(ierr == MPI_SUCCESS)
          else
             ! Allocate memory to receive into
             allocate(recv_buff(ndets_being_sent(i),zoltan_global_ndata_per_det))
             
             ! Receive broadcast
             ewrite(2,*) "Receiving ", ndets_being_sent(i), " detectors from process ", i
             call mpi_bcast(recv_buff,ndets_being_sent(i)*zoltan_global_ndata_per_det, getPREAL(), i-1, MPI_COMM_FEMTOOLS, ierr)
             assert(ierr == MPI_SUCCESS)

             ! Unpack detector if you own it
             do j=1,ndets_being_sent(i)

                ! Allocate and unpack the detector
                shape=>ele_shape(zoltan_global_new_positions,1)                     
                call allocate(detector, zoltan_global_ndims, local_coord_count(shape))
                call unpack_detector(detector, recv_buff(j, 1:zoltan_global_ndata_per_det), zoltan_global_ndims)

                if (has_key(zoltan_global_uen_to_new_local_numbering, detector%element)) then 
                   new_local_element_number = fetch(zoltan_global_uen_to_new_local_numbering, detector%element)
                   if (element_owned(zoltan_global_new_positions%mesh, new_local_element_number)) then
                      detector%element = new_local_element_number
                      call insert(detector, detector_list_array(detector%list_id)%ptr)
                      detector => null()
                   else
                      call delete(detector)
                   end if
                else
                   call delete(detector)
                end if
             end do

             deallocate(recv_buff)
          end if
       end if
    end do
    
    deallocate(ndets_being_sent)
    deallocate(send_buff)
    ewrite(1,*) "Exiting update_detector_list_element"

  end subroutine update_detector_list_element
    
  subroutine update_particle_list_element(detector_list_array)
    ! Update the detector%element field for every particle left in our list
    type(detector_list_ptr), dimension(:), intent(inout) :: detector_list_array

    type(detector_linked_list), pointer :: detector_list => null()
    type(detector_linked_list) ::  detector_send_list
    type(detector_type), pointer :: detector => null(), send_detector => null()
    integer :: i, j, send_count, ierr,k
    integer :: old_local_element_number, new_local_element_number, old_universal_element_number
    integer, allocatable :: ndets_being_sent(:)
    real, allocatable :: send_buff(:,:), recv_buff(:,:)
    logical :: do_broadcast, sent
    integer, dimension(3) :: attribute_size
    integer :: total_attributes
    type(element_type), pointer :: shape

    ewrite(1,*) "In update_particle_list_element"

    !Loop for particles with attributes

    do j = 1, size(detector_list_array)
       send_count=0
       detector_list => detector_list_array(j)%ptr
       ewrite(2,*) "Length of detector list to be updated: ", detector_list%length

       detector => detector_list%first
       if (associated(detector)) then
          if (size(detector%attributes)<1) then
             detector => null()
          else
             attribute_size(1)=size(detector%attributes)
             attribute_size(2)=size(detector%old_attributes)
             attribute_size(3)=size(detector%old_fields)
          end if
       end if
       do while (associated(detector))

          old_local_element_number = detector%element

          if (.not. has_key(zoltan_global_old_local_numbering_to_uen, old_local_element_number)) then
             ewrite(-1,*) "Zoltan can't find old element number for particle ", detector%id_number
             FLAbort('Trying to update unknown particle in Zoltan')
          end if
          old_universal_element_number = fetch(zoltan_global_old_local_numbering_to_uen, old_local_element_number)

          if(has_key(zoltan_global_uen_to_new_local_numbering, old_universal_element_number)) then
             ! Update the element number for the particle
             detector%element = fetch(zoltan_global_uen_to_new_local_numbering, old_universal_element_number)
             detector => detector%next
          else
             ! We no longer own the element containing this particle, and cannot establish its new
             ! owner from the halo, because the boundary has moved too far.
             ! Since we have no way of determining the new owner we are going to broadcast the particle
             ! to all procs, so move it to the send list and count how many particles we're sending.
             ewrite(2,*) "Found non-local particle, initialising broadcast..."
             send_count = send_count + 1
          
             ! Store the old universal element number for unpacking to new local at the receive
             detector%element = old_universal_element_number
             detector%list_id=detector_list%id

             ! Remove detector from detector list
             send_detector => detector
             detector => detector%next
             call move(send_detector, detector_list, detector_send_list)
          end if
       end do
       ! Find out how many particles each process wants to broadcast
       allocate(ndets_being_sent(getnprocs()))
       call mpi_allgather(send_count, 1, getPINTEGER(), ndets_being_sent, 1 , getPINTEGER(), MPI_COMM_FEMTOOLS, ierr)
       assert(ierr == MPI_SUCCESS)
       ! Check whether we have to perform broadcast, if not return
       do_broadcast=.false.
       if (any(ndets_being_sent > 0)) then
          do_broadcast=.true.
       end if
       if (.not. do_broadcast) then
          deallocate(ndets_being_sent)
          cycle
       end if
       i=1
       sent=.false.
       do while (sent.eqv..false.)
          if (ndets_being_sent(i)>0) then
             call mpi_bcast(attribute_size, 3, getPINTEGER(), i-1, MPI_COMM_FEMTOOLS, ierr)
             assert(ierr == MPI_SUCCESS)
             total_attributes=sum(attribute_size)
             sent=.true.
          end if
          i=i+1
       end do
       ewrite(2,*) "Broadcast required, initialising..."
       
       ! Allocate memory for all the particles you're going to send
       allocate(send_buff(send_count,zoltan_global_ndata_per_det+total_attributes))
       
       detector => detector_send_list%first
       do i=1,send_count
          ! Pack the particle information and delete from send_list (delete advances particle to detector%next)
          call pack_detector(detector, send_buff(i,1:zoltan_global_ndata_per_det+total_attributes), zoltan_global_ndims, attribute_size=attribute_size)
          call delete(detector, detector_send_list)
       end do

       ! Broadcast particles whose new owner we can't identify
       do i=1,getnprocs()
          if (ndets_being_sent(i) > 0) then
             
             if (i == getprocno()) then
                ! Broadcast the particles you want to send
                ewrite(2,*) "Broadcasting ", send_count, " particles"
                call mpi_bcast(send_buff,send_count*(zoltan_global_ndata_per_det+total_attributes), getPREAL(), i-1, MPI_COMM_FEMTOOLS, ierr)
                assert(ierr == MPI_SUCCESS)
             else
                ! Allocate memory to receive into
                allocate(recv_buff(ndets_being_sent(i),zoltan_global_ndata_per_det+total_attributes))
                
                ! Receive broadcast
                ewrite(2,*) "Receiving ", ndets_being_sent(i), " particles from process ", i
                call mpi_bcast(recv_buff,ndets_being_sent(i)*(zoltan_global_ndata_per_det+total_attributes), getPREAL(), i-1, MPI_COMM_FEMTOOLS, ierr)
                assert(ierr == MPI_SUCCESS)
                
                ! Unpack particle if you own it
                do k=1,ndets_being_sent(i)
                   
                   ! Allocate and unpack the particle
                   shape=>ele_shape(zoltan_global_new_positions,1)                     
                   call allocate(detector, zoltan_global_ndims, local_coord_count(shape), attribute_size)
                   call unpack_detector(detector, recv_buff(k, 1:zoltan_global_ndata_per_det+total_attributes), zoltan_global_ndims, attribute_size=attribute_size)
                   
                   if (has_key(zoltan_global_uen_to_new_local_numbering, detector%element)) then 
                      new_local_element_number = fetch(zoltan_global_uen_to_new_local_numbering, detector%element)
                      if (element_owned(zoltan_global_new_positions%mesh, new_local_element_number)) then
                         detector%element = new_local_element_number
                         call picker_inquire(zoltan_global_new_positions, detector%position, detector%element, detector%local_coords, global=.false.)
                         call insert(detector, detector_list_array(detector%list_id)%ptr)
                         detector => null()
                      else
                         call delete(detector)
                      end if
                   else
                      call delete(detector)
                   end if
                end do
                
                deallocate(recv_buff)
             end if
          end if
       end do
       deallocate(send_buff)
       deallocate(ndets_being_sent)
    end do

    ewrite(1,*) "Exiting update_particle_list_element"

  end subroutine update_particle_list_element

  subroutine transfer_fields(zz)
    ! OK! So, here is how this is going to work. We are going to 
    ! loop through every element in which you own at least one node, and note
    ! that its information needs to be sent to the owners of its vertices.
    ! We also have to take special care of self-sends, since they don't
    ! get taken care of in the zoltan communication.

    type(zoltan_struct), pointer, intent(in) :: zz    
    
    integer :: old_ele
    integer, dimension(:), pointer :: old_local_nodes, nodes
    type(element_type), pointer :: eshape
    type(integer_set), dimension(halo_proc_count(zoltan_global_zz_halo)) :: sends
    integer :: i, j, new_owner, universal_element_number
    type(integer_set) :: self_sends
    integer :: num_import, num_export
    integer :: original_zoltan_global_unpacked_detectors_list_length
    integer, dimension(:,:), allocatable :: vertex_order
    integer, dimension(:), pointer :: import_global_ids, import_local_ids, import_procs, import_to_part
    integer, dimension(:), pointer :: export_global_ids, export_local_ids, export_procs, export_to_part
    integer :: head
    integer(zoltan_int) :: ierr
    
    integer :: old_universal_element_number, new_local_element_number, old_local_element_number
    integer :: state_no, field_no
    type(scalar_field), pointer :: source_sfield, target_sfield
    type(vector_field), pointer :: source_vfield, target_vfield
    type(tensor_field), pointer :: source_tfield, target_tfield

    type(detector_list_ptr), dimension(:), pointer :: detector_list_array => null()
    type(detector_type), pointer :: detector => null(), add_detector => null()

    ewrite(1,*) 'in transfer_fields'
    
    do i=1,size(sends)
       call allocate(sends(i))
    end do
    call allocate(self_sends)
    
    
    do old_ele=1,ele_count(zoltan_global_zz_positions)
       universal_element_number = halo_universal_number(zoltan_global_zz_ele_halo, old_ele)
       old_local_nodes => ele_nodes(zoltan_global_zz_positions, old_ele)
       if (.not. any(nodes_owned(zoltan_global_zz_halo, old_local_nodes))) cycle
       do i=1,size(old_local_nodes)
          if (has_value(zoltan_global_nodes_we_are_keeping, old_local_nodes(i))) then
             assert(node_owned(zoltan_global_zz_halo, old_local_nodes(i)))
             call insert(self_sends, universal_element_number)
          else if (has_key(zoltan_global_nodes_we_are_sending,  old_local_nodes(i))) then
             new_owner = fetch(zoltan_global_nodes_we_are_sending, old_local_nodes(i))
             call insert(sends(new_owner+1), universal_element_number)
          end if
       end do
    end do
    
    num_export = sum(key_count(sends))
    allocate(export_global_ids(num_export))
    allocate(export_procs(num_export))

    ! allocate array for storing the number of detectors in each of the elements to be transferred
    allocate(zoltan_global_ndets_in_ele(num_export))
    zoltan_global_ndets_in_ele(:) = 0

    ! calculate the amount of data to be transferred per detector
    zoltan_global_ndims = zoltan_global_zz_positions%dim
    zoltan_global_ndata_per_det = detector_buffer_size(zoltan_global_ndims, .false.)
    ewrite(2,*) "Amount of data to be transferred per detector: ", zoltan_global_ndata_per_det
    
    head = 1
    do i=1,size(sends)
       export_global_ids(head:head + key_count(sends(i)) - 1) = set2vector(sends(i))
       export_procs(head:head + key_count(sends(i)) - 1) = i - 1
       head = head + key_count(sends(i))
    end do
    
    allocate(export_local_ids(num_export))
    export_local_ids = 666
    
    import_global_ids => null()
    import_local_ids => null()
    import_procs => null()
    import_to_part => null()
    export_to_part => null()
    
    ierr = Zoltan_Compute_Destinations(zz, &
         & num_export, export_global_ids, export_local_ids, export_procs, &
         & num_import, import_global_ids, import_local_ids, import_procs)
    assert(ierr == ZOLTAN_OK)

    ! Get all detector lists
    call get_registered_detector_lists(detector_list_array)

    ! Log list lengths
    if (get_num_detector_lists()>0) then
       ewrite(2,*) "Before migrate, we have", get_num_detector_lists(), "detector lists:"
       do j = 1, size(detector_list_array)
          ewrite(2,*) "Detector list", j, "has", detector_list_array(j)%ptr%length, "local and ", detector_list_array(j)%ptr%total_num_det, "global detectors"
       end do
    end if

    ierr = Zoltan_Migrate(zz, num_import, import_global_ids, import_local_ids, import_procs, &
         & import_to_part, num_export, export_global_ids, export_local_ids, export_procs, export_to_part)
    
    assert(ierr == ZOLTAN_OK)
    
    deallocate(export_local_ids)
    deallocate(export_procs)
    deallocate(export_global_ids)

    deallocate(zoltan_global_ndets_in_ele)
    
    ! update the local detectors and make sure we didn't miss any in the first send
    if (get_num_detector_lists()>0) then
       call update_detector_list_element(detector_list_array)
       call update_particle_list_element(detector_list_array)
    end if

    ! Merge in any detectors we received as part of the transfer to our detector list
    detector => zoltan_global_unpacked_detectors_list%first
    original_zoltan_global_unpacked_detectors_list_length = zoltan_global_unpacked_detectors_list%length

    do j=1, original_zoltan_global_unpacked_detectors_list_length
       add_detector => detector
       detector => detector%next

       add_detector%name=int2str(add_detector%id_number)!Temporary change to fix particle spawning issues in parallel

       ! move detector to the correct list
       call move(add_detector, zoltan_global_unpacked_detectors_list, detector_list_array(add_detector%list_id)%ptr)
    end do

    assert(zoltan_global_unpacked_detectors_list%length==0)
    ewrite(2,*) "Merged", original_zoltan_global_unpacked_detectors_list_length, "detectors with local detector lists"

    ! Log list lengths
    if (get_num_detector_lists()>0) then
       call get_registered_detector_lists(detector_list_array)
       ewrite(2,*) "After migrate and merge, we have", get_num_detector_lists(), "detector lists:"
       do j = 1, size(detector_list_array)
          ewrite(2,*) "Detector list", j, "has", detector_list_array(j)%ptr%length, "local and ", detector_list_array(j)%ptr%total_num_det, "global detectors"
       end do
    end if

    ierr = Zoltan_LB_Free_Part(import_global_ids, import_local_ids, import_procs, import_to_part)
    assert(ierr == ZOLTAN_OK)
    
    ! for all self-send elements, establish the vertex order in the new element such
    ! that it matches the local element ordering of the old element
    allocate(vertex_order(1:ele_loc(zoltan_global_new_positions,1), key_count(self_sends)))
    do i=1, key_count(self_sends)
      old_universal_element_number = fetch(self_sends, i)
      new_local_element_number = fetch(zoltan_global_uen_to_new_local_numbering, old_universal_element_number)
      old_local_element_number = fetch(zoltan_global_uen_to_old_local_numbering, old_universal_element_number)
      ! this function takes the universal node numbers of the old element and the new global (over the local domain)
      ! node numbers, to compute the vertex order
      vertex_order(:,i) = local_vertex_order( &
         halo_universal_number(zoltan_global_zz_halo, ele_nodes(zoltan_global_zz_positions, old_local_element_number)), &
         ele_nodes(zoltan_global_new_positions, new_local_element_number))
    end do
    
    do state_no=1,size(zoltan_global_source_states)
       assert(scalar_field_count(zoltan_global_source_states(state_no)) == scalar_field_count(zoltan_global_target_states(state_no)))
       assert(vector_field_count(zoltan_global_source_states(state_no)) == vector_field_count(zoltan_global_target_states(state_no)))
       assert(tensor_field_count(zoltan_global_source_states(state_no)) == tensor_field_count(zoltan_global_target_states(state_no)))
       
       do field_no=1,scalar_field_count(zoltan_global_source_states(state_no))
          source_sfield => extract_scalar_field(zoltan_global_source_states(state_no), field_no)
          target_sfield => extract_scalar_field(zoltan_global_target_states(state_no), field_no)
          assert(trim(source_sfield%name) == trim(target_sfield%name))
          
          do i=1,key_count(self_sends)
             old_universal_element_number = fetch(self_sends, i)
             new_local_element_number = fetch(zoltan_global_uen_to_new_local_numbering, old_universal_element_number)
             old_local_element_number = fetch(zoltan_global_uen_to_old_local_numbering, old_universal_element_number)
             eshape => ele_shape(target_sfield, new_local_element_number)
             nodes => ele_nodes(target_sfield, new_local_element_number)
             call set(target_sfield, nodes(ele_local_num(vertex_order(:,i), eshape%numbering)), &
                               ele_val(source_sfield, old_local_element_number))
          end do
       end do
       
       do field_no=1,vector_field_count(zoltan_global_source_states(state_no))
          source_vfield => extract_vector_field(zoltan_global_source_states(state_no), field_no)
          target_vfield => extract_vector_field(zoltan_global_target_states(state_no), field_no)
          assert(trim(source_vfield%name) == trim(target_vfield%name))
          if (source_vfield%name == zoltan_global_new_positions%name) cycle
          
          do i=1,key_count(self_sends)
             old_universal_element_number = fetch(self_sends, i)
             new_local_element_number = fetch(zoltan_global_uen_to_new_local_numbering, old_universal_element_number)
             old_local_element_number = fetch(zoltan_global_uen_to_old_local_numbering, old_universal_element_number)
             eshape => ele_shape(target_vfield, new_local_element_number)
             nodes => ele_nodes(target_vfield, new_local_element_number)
             call set(target_vfield, nodes(ele_local_num(vertex_order(:,i), eshape%numbering)), &
                               ele_val(source_vfield, old_local_element_number))
          end do
       end do
       
       do field_no=1,tensor_field_count(zoltan_global_source_states(state_no))
          source_tfield => extract_tensor_field(zoltan_global_source_states(state_no), field_no)
          target_tfield => extract_tensor_field(zoltan_global_target_states(state_no), field_no)
          assert(trim(source_tfield%name) == trim(target_tfield%name))
          
          do i=1,key_count(self_sends)
             old_universal_element_number = fetch(self_sends, i)
             new_local_element_number = fetch(zoltan_global_uen_to_new_local_numbering, old_universal_element_number)
             old_local_element_number = fetch(zoltan_global_uen_to_old_local_numbering, old_universal_element_number)
             eshape => ele_shape(target_tfield, new_local_element_number)
             nodes => ele_nodes(target_tfield, new_local_element_number)
             call set(target_tfield, nodes(ele_local_num(vertex_order(:,i), eshape%numbering)), &
                               ele_val(source_tfield, old_local_element_number))
          end do
       end do
    end do
    
    call deallocate(self_sends)
    call deallocate(sends)
    deallocate(vertex_order)    
    
    call halo_update(zoltan_global_target_states)
    
    ewrite(1,*) 'exiting transfer_fields'
    
  end subroutine transfer_fields

  subroutine finalise_transfer(states, metric, full_metric, new_metric)
    type(state_type), dimension(:), intent(inout), target :: states

    type(tensor_field), intent(inout), optional :: metric
    type(tensor_field), intent(inout), optional :: full_metric
    type(tensor_field), intent(in) :: new_metric

    integer :: i
    call set_prescribed_field_values(states, exclude_interpolated = .true.)
    call populate_boundary_conditions(states)
    call set_boundary_conditions_values(states)
    call set_dirichlet_consistent(states)
    call alias_fields(states)
    
    if (present(full_metric)) then
       full_metric = new_metric
       call halo_update(full_metric)
    else if (present(metric)) then
       metric = new_metric
       call halo_update(metric)
    end if
    
    do i=1,size(zoltan_global_source_states)
       call deallocate(zoltan_global_source_states(i))
       call deallocate(zoltan_global_target_states(i))
    end do
    deallocate(zoltan_global_source_states)
    deallocate(zoltan_global_target_states)
  end subroutine finalise_transfer

  subroutine dump_linear_mesh
    type(scalar_field) :: sends, receives, unn
    integer :: i, proc
    
    assert(associated(zoltan_global_new_positions%refcount))
    assert(zoltan_global_new_positions%refcount%count == 1)
    assert(associated(zoltan_global_new_positions%mesh%refcount))
    assert(zoltan_global_new_positions%mesh%refcount%count == 1)
    
    call allocate(sends, zoltan_global_new_positions%mesh, "Sends")
    call zero(sends)
    call allocate(receives, zoltan_global_new_positions%mesh, "Receives")
    call zero(receives)
    call allocate(unn, zoltan_global_new_positions%mesh, "NewUniversalNodeNumber")
    call zero(unn)
    
    do proc=1,halo_proc_count(zoltan_global_new_positions%mesh%halos(2))
       do i=1,size(zoltan_global_new_positions%mesh%halos(2)%sends(proc)%ptr)
          call set(sends, zoltan_global_new_positions%mesh%halos(2)%sends(proc)%ptr(i), 1.0)
       end do
       do i=1,size(zoltan_global_new_positions%mesh%halos(2)%receives(proc)%ptr)
          call set(receives, zoltan_global_new_positions%mesh%halos(2)%receives(proc)%ptr(i), 1.0)
       end do
    end do
    
    do i=1,node_count(zoltan_global_new_positions)
       call set(unn, i, float(halo_universal_number(zoltan_global_new_positions%mesh%halos(2), i)))
    end do
    
    call deallocate(sends)
    call deallocate(receives)
    call deallocate(unn)
  end subroutine dump_linear_mesh

  subroutine dump_suggested_owner(states, p1_num_export, p1_export_local_ids, p1_export_procs)
    type(state_type), dimension(:), intent(inout), target :: states
    integer(zoltan_int), intent(in) :: p1_num_export
    integer, dimension(:), pointer, intent(in) :: p1_export_local_ids, p1_export_procs

    integer :: rank, i
    type(scalar_field) :: suggested_owner, unn
    type(vector_field) :: positions
    
    rank = getrank()
    call allocate(suggested_owner, zoltan_global_zz_mesh, "SuggestedOwner")
    
    call set(suggested_owner, float(rank))
    do i=1,p1_num_export
       call set(suggested_owner, p1_export_local_ids(i), float(p1_export_procs(i)))
    end do
    
    call allocate(unn, zoltan_global_zz_mesh, "OldUniversalNodeNumber")
    do i=1,node_count(unn)
       call set(unn, i, float(halo_universal_number(zoltan_global_zz_halo, i)))
    end do
    
    positions = get_coordinate_field(states(1), zoltan_global_zz_mesh)
    call halo_update(suggested_owner)
    call deallocate(positions)
    call deallocate(suggested_owner)
    call deallocate(unn)
    
  end subroutine dump_suggested_owner

#endif
end module zoltan_integration
