!    Copyrigh (C) 2006 Imperial College London and others.
!
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineeringp
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

module checkpoint
  !!< Checkpointing routines

  use fields
  use fields_data_types
  use field_options
  use fldebug
  use futils
  use global_parameters, only : FIELD_NAME_LEN, OPTION_PATH_LEN, simulation_start_time
  use halos
  use parallel_tools
  use spud
  use state_module
  use vtk_interfaces
  use mesh_files
  use detector_data_types
  use diagnostic_variables
  use mpi_interfaces

  implicit none

  private

  public :: do_checkpoint_simulation, checkpoint_simulation, checkpoint_detectors, checkpoint_check_options

  integer, parameter :: PREFIX_LEN = OPTION_PATH_LEN

  integer, save :: fhdet=0

contains

  function do_checkpoint_simulation(dump_no) result(do_checkpoint)
    !!< Determine whether a checkpoint should be written.

    integer, intent(in) :: dump_no

    logical :: do_checkpoint

    integer :: cp_dump_period
    real :: current_time

    ! Is checkpointing enabled?
    if(.not. have_option("/io/checkpointing")) then
      do_checkpoint = .false.
    else
      call get_option("/io/checkpointing/checkpoint_period_in_dumps", cp_dump_period)
      call get_option("/timestepping/current_time", current_time)

      ! Permit a checkpointing period of zero
      cp_dump_period = max(cp_dump_period, 1)

      ! Checkpoint if cp_no is divisible by cp_dump_period
      if(mod(dump_no, cp_dump_period) /= 0) then
        do_checkpoint = .false.
      ! unless this is at simulation start and checkpointing is not enabled at simulation start
      else if(current_time == simulation_start_time .and. .not. have_option("/io/checkpointing/checkpoint_at_start")) then
        do_checkpoint = .false.
      else
        do_checkpoint = .true.
      end if
    end if

    if(do_checkpoint) then
      ewrite(1, *) "do_checkpoint returning .true."
    else
      ewrite(1, *) "do_checkpoint returning .false."
    end if

  end function do_checkpoint_simulation

  subroutine checkpoint_simulation(state, prefix, postfix, cp_no, protect_simulation_name, &
    keep_initial_data, ignore_detectors, number_of_partitions)
    !!< Checkpoint the whole simulation
    
    type(state_type), dimension(:), intent(in) :: state
    !! Default value simulation_name
    character(len = *), optional, intent(in) :: prefix
    !! Default value "checkpoint"
    character(len = *), optional, intent(in) :: postfix
    integer, optional, intent(in) :: cp_no
    !! If present and .false., do not protect the simulation_name when
    !! checkpointing the options tree
    logical, optional, intent(in) :: protect_simulation_name
    !! If present and .true.: do not checkpoint fields that can be reinitialsed and do not 
    !! checkpoint extruded meshes if the extrusion can be repeated using the initial sizing_function, 
    !! i.e. if this run has not been started with a checkpointed extruded mesh (extrude/checkpoint_from_file)
    logical, optional, intent(in) :: keep_initial_data
    !! When using flredecomp to re-partition the domain, the detectors
    !! are ignored since they are kept in one header and one data file,
    !! not in a set of per-process files.  So flredecomp does not want to
    !! checkpoint detectors.
    logical, optional, intent(in) :: ignore_detectors
    !! If present, only write for processes 1:number_of_partitions (assumes the other partitions are empty)
    integer, optional, intent(in):: number_of_partitions

    character(len = PREFIX_LEN) :: lpostfix, lprefix

    ewrite(1, *) "Checkpointing simulation"

    if(present(prefix)) then
      lprefix = prefix
    else
      call get_option("/simulation_name", lprefix)
    end if
    if(len_trim(lprefix) == 0) then
      FLAbort("Checkpoint prefix cannot have length zero")
    end if
    
    if(present(postfix)) then
      lpostfix = postfix
    else
      lpostfix = "checkpoint"
    end if

    call checkpoint_state(state, lprefix, postfix = lpostfix, cp_no = cp_no, &
      keep_initial_data = keep_initial_data, number_of_partitions=number_of_partitions)
    if(have_option("/io/detectors") &
         .and. .not.present_and_true(ignore_detectors)) then
      call checkpoint_detectors(state, lprefix, postfix = lpostfix, cp_no = cp_no)
    end if
    if(getrank() == 0) then
      ! Only rank zero should write out the options tree in parallel
      call checkpoint_options(lprefix, postfix = lpostfix, cp_no = cp_no, &
        & protect_simulation_name = .not. present_and_false(protect_simulation_name))
    end if
    
  end subroutine checkpoint_simulation

  subroutine checkpoint_detectors(state, prefix, postfix, cp_no)
    !!< Checkpoint detectors

    type(state_type), dimension(:), intent(in) :: state
    character(len = *), intent(in) :: prefix
    !! Default value "checkpoint"
    character(len = *), optional, intent(in) :: postfix
    integer, optional, intent(in) :: cp_no

    type(detector_type), pointer :: node
    character(len = OPTION_PATH_LEN) :: detectors_cp_filename, path
    character(len = PREFIX_LEN) :: lpostfix
    integer :: det_unit, i, j, IERROR
    integer :: static_dete, lagrangian_dete, det_arrays, array_dete_lag, &
      & total_dete_lag
    type(vector_field), pointer :: vfield

    integer(KIND=MPI_OFFSET_KIND) :: location_to_write, offset
    integer, ALLOCATABLE, DIMENSION(:) :: status
    integer :: nints, realsize, dimen, total_num_det, number_total_columns
    real, dimension(:), allocatable :: buffer

    ewrite(1, *) "Checkpointing detectors"

    assert(len_trim(prefix) > 0)
    
    if(present(postfix)) then
      lpostfix = postfix
    else
      lpostfix = "checkpoint"
    end if

    static_dete = option_count("/io/detectors/static_detector")
    lagrangian_dete = option_count("/io/detectors/lagrangian_detector")
    det_arrays = option_count("/io/detectors/detector_array")
    array_dete_lag = 0
 
    do i = 1, det_arrays
      path = "/io/detectors/detector_array[" // int2str(i - 1) // "]"
      if(have_option(trim(path) // "/lagrangian")) then
        call get_option(trim(path) // "/number_of_detectors", j)
        array_dete_lag = array_dete_lag + j
     end if
    end do
    total_dete_lag = lagrangian_dete + array_dete_lag
    if(total_dete_lag == 0) then
      ewrite(1, *) "No Lagrangian detectors - not checkpointing detectors"
      ewrite(1, *) "Exiting checkpoint_detectors"
      return
    end if

    ! Construct a new detector checkpoint filename
    detectors_cp_filename = trim(prefix)
    if(present(cp_no)) detectors_cp_filename = trim(detectors_cp_filename) // "_" // int2str(cp_no)
    detectors_cp_filename = trim(detectors_cp_filename) // "_" // trim(lpostfix)

!!!!! Writing of position detectors before checkpointing in serial  !!!!!

    if(getprocno() == 1) then
      det_unit = free_unit()

      ! Write the detectors positions of the Lagrangian detectors into the
      ! checkpoint detector file in binary format

      ! Before writing the positions, we write a header that contains the names
      ! of the groups of detectors in the order that they were read

      open(unit = det_unit, &
        & file = trim(detectors_cp_filename) // '_det.groups', &
        & action = "write")
      do i = 1, size(default_stat%detector_group_names) 
        write(det_unit,'(a,i0)') &
          & default_stat%detector_group_names(i), default_stat%number_det_in_each_group(i)
      end do
      close(det_unit)
    
    end if

    !!< Writes detector last position into detectors file using MPI output 
    ! commands so that when running in parallel all processors can write at the same time information into the file at the right location.
    
    call MPI_FILE_OPEN(MPI_COMM_FEMTOOLS, trim(detectors_cp_filename) // '_det.positions.dat', MPI_MODE_CREATE + MPI_MODE_RDWR, MPI_INFO_NULL, fhdet, IERROR)

    ewrite(1,*) "after openning the IERROR is:", IERROR

    allocate( status(MPI_STATUS_SIZE) )

    call MPI_TYPE_EXTENT(getpreal(), realsize, ierror)

    vfield => extract_vector_field(state(1),"Velocity")

    dimen=vfield%dim

    total_num_det = 0

    do i =1, size(default_stat%number_det_in_each_group)

      total_num_det=total_num_det+default_stat%number_det_in_each_group(i)

    end do

    number_total_columns=total_num_det*dimen

    node => default_stat%detector_list%first

    location_to_write=0

    positionloop_cp: do i=1, default_stat%detector_list%length
      offset = location_to_write+(node%id_number-1)*size(node%position)*realsize
      ewrite(1,*) "after file set view position IERROR is:", IERROR

      allocate(buffer(size(node%position)))
      buffer=node%position
      nints=size(node%position)

      call MPI_FILE_WRITE_AT(fhdet,offset,buffer,nints,getpreal(),status,IERROR)

      ewrite(1,*) "after sync position IERROR is:", IERROR
      deallocate(buffer)
      node => node%next
    end do positionloop_cp

    call update_detectors_options(trim(detectors_cp_filename) // "_det", "binary")

    if (fhdet/=0) then
       call MPI_FILE_CLOSE(fhdet, IERROR) 
       if (IERROR/=0) then
          ewrite(0,*) "Warning: failed to close .detector checkpoint file open with mpi_file_open"
       end if
    end if
    
    ewrite(1, *) "Exiting detectors"

  end subroutine checkpoint_detectors  

  subroutine update_detectors_options(filename,format)

    !!< Updates the initial options of the detectors (options tree in diamond)

    character(len = *), intent(in) :: filename
    character(len = *), intent(in) :: format

    integer :: stat, i, python_or_file_dete, python_functions_or_files, static_dete, lagrangian_dete
    character(len = FIELD_NAME_LEN), dimension(:), allocatable :: type_detectors
    character(len = 254) :: temp_string

    static_dete = option_count("/io/detectors/static_detector")
    lagrangian_dete = option_count("/io/detectors/lagrangian_detector")
    python_functions_or_files = option_count("/io/detectors/detector_array")
    python_or_file_dete = 0

    do i = 0, static_dete-1  
       call delete_option("/io/detectors/static_detector[" // int2str(0) // "]")
    end do

    do i = 0, static_dete-1
       temp_string=default_stat%detector_group_names(i+1)

        ewrite(1,*) 'In update_detectors_options static det loop'
        ewrite(1,*) temp_string

       call set_option_attribute("/io/detectors/static_detector::" // trim(temp_string) // "/from_checkpoint_file/file_name", trim(filename), stat)

        if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
            FLAbort("Failed to set detectors options filename when checkpointing detectors with option path " // "/io/detectors/static_detector")
        end if

        call set_option("/io/detectors/static_detector::" // trim(temp_string) // "/from_checkpoint_file/format/", trim(format), stat)

        if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING) then
            FLAbort("Failed to set detectors options format when checkpointing detectors with option path " // "/io/detectors/static_detector")
        end if  
    end do

    do i = 0, lagrangian_dete-1  
       call delete_option("/io/detectors/lagrangian_detector[" // int2str(0) // "]")
    end do

    do i = 0, lagrangian_dete-1

       temp_string=default_stat%detector_group_names(i+1+static_dete)
        
       call set_option_attribute("/io/detectors/lagrangian_detector::" // trim(temp_string) // "/from_checkpoint_file/file_name", trim(filename), stat)

        if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
            FLAbort("Failed to set detectors options filename when checkpointing detectors with option path " // "/io/detectors/lagrangian_detector")
        end if

        call set_option("/io/detectors/lagrangian_detector::" // trim(temp_string) // "/from_checkpoint_file/format/", trim(format), stat)

        if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING) then
            FLAbort("Failed to set detectors options format when checkpointing detectors with option path " // "/io/detectors/lagrangian_detector")
        end if  
    end do

    allocate(type_detectors(python_functions_or_files))

    do i = 0, python_functions_or_files-1  

         if (have_option("/io/detectors/detector_array[" // int2str(0) // "]"//"/lagrangian")) then
             type_detectors(i+1)='LAGRANGIAN'
         else
             type_detectors(i+1)='STATIC'
         end if

         call delete_option("/io/detectors/detector_array[" // int2str(0) // "]")

    end do

    do i = 0, python_functions_or_files-1  
        temp_string=default_stat%detector_group_names(i+1+static_dete+lagrangian_dete)

        ewrite(1,*) 'In update_detectors_options'
        ewrite(1,*) temp_string

        ewrite(1,*) 'In update_detectors_options'
        ewrite(1,*) temp_string

        call set_option("/io/detectors/detector_array::" // trim(temp_string) // "/number_of_detectors/", &
                & default_stat%number_det_in_each_group(i+1+static_dete+lagrangian_dete), stat = stat) 

        ewrite(1,*) 'In update_detectors_options'
        ewrite(1,*) default_stat%number_det_in_each_group(i+1+static_dete+lagrangian_dete)

        assert(any(stat == (/SPUD_NO_ERROR, SPUD_NEW_KEY_WARNING/)))

        ewrite(1,*) 'In update_detectors_options'
        ewrite(1,*) default_stat%number_det_in_each_group(i+1+static_dete+lagrangian_dete)
        
        if (type_detectors(i+1)=='LAGRANGIAN') then
             call add_option("/io/detectors/detector_array::" // trim(temp_string) // "/lagrangian", stat = stat) 
             assert(any(stat == (/SPUD_NO_ERROR, SPUD_NEW_KEY_WARNING/)))
        else
             call add_option("/io/detectors/detector_array::" // trim(temp_string) // "/static", stat = stat)
             assert(any(stat == (/SPUD_NO_ERROR, SPUD_NEW_KEY_WARNING/)))
        end if

        call set_option_attribute("/io/detectors/detector_array::" // trim(temp_string) // "/from_checkpoint_file/file_name", trim(filename), stat)

        if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
            FLAbort("Failed to set detectors options filename when checkpointing detectors with option path " // "io/detectors/detector_array")
        end if

        call set_option("/io/detectors/detector_array::" // trim(temp_string) // "/from_checkpoint_file/format/", trim(format), stat)

        if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING) then
            FLAbort("Failed to set detectors options format when checkpointing detectors with option path " // "/io/detectors/detector_array")
        end if  
    end do

  deallocate(type_detectors)

  end subroutine update_detectors_options

  subroutine checkpoint_state(state, prefix, postfix, cp_no, keep_initial_data, number_of_partitions)
    !!< Checkpoint state.

    type(state_type), dimension(:), intent(in) :: state
    character(len = *), intent(in) :: prefix
    character(len = *), optional, intent(in) :: postfix
    integer, optional, intent(in) :: cp_no
    !! If present and .true.: do not checkpoint fields that can be reinitialsed and do not 
    !! checkpoint extruded meshes if the extrusion can be repeated using the initial sizing_function, 
    !! i.e. if this run has not been started with a checkpointed extruded mesh (extrude/checkpoint_from_file)
    logical, optional, intent(in) :: keep_initial_data
    !! If present, only write for processes 1:number_of_partitions (assumes the other partitions are empty)
    integer, optional, intent(in):: number_of_partitions

    call checkpoint_meshes(state, prefix, postfix, cp_no, keep_initial_data=keep_initial_data, number_of_partitions=number_of_partitions)
    call checkpoint_fields(state, prefix, postfix, cp_no, keep_initial_data=keep_initial_data, number_of_partitions=number_of_partitions)

  end subroutine checkpoint_state

  subroutine checkpoint_meshes(state, prefix, postfix, cp_no, keep_initial_data, number_of_partitions)
    !!< Checkpoint the meshes in state. Outputs to mesh files with names:
    !!<   [prefix]_[mesh_name][_cp_no][_postfix][_process].[extention]
    !!< where cp_no is optional and the process number is added in parallel.
    !!< Also outputs a .halo file if running in parallel.

    type(state_type), dimension(:), intent(in) :: state
    character(len = *), intent(in) :: prefix
    character(len = *), optional, intent(in) :: postfix
    integer, optional, intent(in) :: cp_no
    ! if present and true: do not checkpoint extruded meshes that can be re-extruded
    logical, optional, intent(in) :: keep_initial_data
    !! If present, only write for processes 1:number_of_partitions (assumes the other partitions are empty)
    integer, optional, intent(in):: number_of_partitions

    type(vector_field), pointer:: position
    character(len = FIELD_NAME_LEN) :: mesh_name, mesh_format
    character(len = OPTION_PATH_LEN) :: mesh_path, mesh_filename
    integer :: i, n_meshes, stat1, stat2, nparts
    type(mesh_type), pointer :: mesh, external_mesh
    logical :: from_file, extruded
    
    assert(len_trim(prefix) > 0)

    if (present(number_of_partitions)) then
      nparts = number_of_partitions
    else
      nparts = getnprocs()
    end if

    n_meshes = option_count("/geometry/mesh")
    do i = 0, n_meshes - 1
      ! Dump each mesh listed under /geometry/mesh that is from_file
      mesh_path = "/geometry/mesh[" // int2str(i) // "]"
      
      from_file = have_option(trim(mesh_path) // "/from_file")
      extruded = have_option(trim(mesh_path) // "/from_mesh/extrude") .and. &
        .not. (present_and_true(keep_initial_data) .and. &
               .not. have_option(trim(mesh_path) // "/from_mesh/extrude/checkpoint_from_file"))
        
      if(from_file .or. extruded) then
        ! Find the mesh (looking in first state)
        call get_option(trim(mesh_path) // "/name", mesh_name)
        mesh => extract_mesh(state(1), trim(mesh_name))

        ewrite(2, *) "Checkpointing mesh " // trim(mesh_name)

        ! Construct a new mesh filename
        mesh_filename = trim(prefix) // "_" // trim(mesh%name)
        if(present(cp_no)) mesh_filename = trim(mesh_filename) // "_" // int2str(cp_no)
        if(present_and_nonempty(postfix)) mesh_filename = trim(mesh_filename) // "_" // trim(postfix)

        ! Update the options tree (required for options tree checkpointing)
        if (from_file) then
          call set_option_attribute(trim(mesh_path) // "/from_file/file_name", trim(mesh_filename))
          call get_option(trim(mesh_path) // "/from_file/format/name", mesh_format)
        else if (extruded) then

          ! the mesh format is determined from the external mesh
          external_mesh => get_external_mesh(state)
          call get_option(trim(external_mesh%option_path) // "/from_file/format/name", mesh_format)

          call set_option_attribute(trim(mesh_path) // "/from_mesh/extrude/checkpoint_from_file/format/name", trim(mesh_format), stat=stat1)
          call set_option_attribute(trim(mesh_path) // "/from_mesh/extrude/checkpoint_from_file/file_name", trim(mesh_filename), stat=stat2)
          if ((stat1/=SPUD_NO_ERROR .and. stat1/=SPUD_NEW_KEY_WARNING) .or. &
             & (stat2/=SPUD_NO_ERROR .and. stat2/=SPUD_NEW_KEY_WARNING)) then
            FLAbort("Failed to modify extrude options for checkpointing.")
          end if
        end if

        ! Write out the mesh using a suitable coordinate field
        if (mesh%name=="CoordinateMesh") then
          if (have_option("/mesh_adaptivity/mesh_movement/free_surface")) then
            ! we don't want/need to checkpoint the moved mesh, as the mesh movement will again be applied after the restart
            ! based on the checkpointed FreeSurface field.
            position => extract_vector_field(state(1), "OriginalCoordinate", stat=stat1)
            if (stat1/=0) then
              ! some cases (e.g. flredecomp) OriginalCoordinate doesn't exist
              position => extract_vector_field(state(1), "Coordinate")
            end if
          else
            position => extract_vector_field(state(1), "Coordinate")
          end if
        else
          position => extract_vector_field(state(1), trim(mesh%name)//"Coordinate")
        end if

        if(nparts > 1) then
          call write_mesh_files(parallel_filename(mesh_filename), mesh_format, position, number_of_partitions=number_of_partitions)
          ! Write out the halos
          ewrite(2, *) "Checkpointing halos"
          call write_halos(mesh_filename, mesh, number_of_partitions=number_of_partitions)
        else
          ! Write out the mesh
          call write_mesh_files(mesh_filename, mesh_format, position, number_of_partitions=number_of_partitions)
        end if

      end if

   end do

  end subroutine checkpoint_meshes

  subroutine checkpoint_fields(state, prefix, postfix, cp_no, keep_initial_data, number_of_partitions)
    !!< Checkpoint the fields in state. Outputs to vtu files with names:
    !!<   [prefix]_[_state name]_[mesh_name][_cp_no][_postfix][_process].vtu
    !!< where the state name is added if multiple states are passed, cp_no is
    !!< optional and the process number is added in parallel.

    type(state_type), dimension(:), intent(in) :: state
    character(len = *), intent(in) :: prefix
    character(len = *), optional, intent(in) :: postfix
    integer, optional, intent(in) :: cp_no
    ! if present and true: do not checkpoint fields that can be reinitialised
    logical, optional, intent(in) :: keep_initial_data
    !! If present, only write for processes 1:number_of_partitions (assumes the other partitions are empty)
    integer, optional, intent(in):: number_of_partitions

    character(len = OPTION_PATH_LEN) :: vtu_filename
    integer :: i, j, k, nparts, n_ps_fields_on_mesh, n_pv_fields_on_mesh, n_pt_fields_on_mesh
    type(mesh_type), pointer :: mesh
    type(scalar_field), pointer :: s_field
    type(scalar_field), dimension(:), allocatable :: ps_fields_on_mesh
    type(vector_field), pointer :: positions, v_field
    type(vector_field), dimension(:), allocatable :: pv_fields_on_mesh
    type(tensor_field), pointer :: t_field
    type(tensor_field), dimension(:), allocatable :: pt_fields_on_mesh

    integer :: stat

    assert(len_trim(prefix) > 0)

    if (present(number_of_partitions)) then
      nparts = number_of_partitions
    else
      nparts = getnprocs()
    end if

    do i = 1, size(state)
      if (have_option("/mesh_adaptivity/mesh_movement/free_surface")) then
        ! we don't want/need to checkpoint the moved mesh, as the mesh movement will again be applied after the restart
        ! based on the checkpointed FreeSurface field.
        positions => extract_vector_field(state(i), "OriginalCoordinate", stat=stat)
        if (stat/=0) then
          ! some cases (e.g. flredecomp) OriginalCoordinate doesn't exist
          positions => extract_vector_field(state(1), "Coordinate")
        end if
      else
        if (have_option("/mesh_adaptivity/mesh_movement")) then
          ! for other mesh movement schemes we should probably write out the "moved" mesh to retain that information
          ! However if the checkpointed mesh is not the "CoordinateMesh", it will have its own MeshNameCoordinate field
          ! that hasn't moved; leading to a failure to restart from the checkpoint. This is only one of the possible problems
          ! generally this functionality is untested.
          ewrite(0,*) "WARNING: using mesh_movement with checkpointing is untested and likely broken."
        end if
        positions => extract_vector_field(state(i), "Coordinate")
      end if
      do j = 1, size(state(i)%meshes)
        mesh => state(i)%meshes(j)%ptr

        ! Construct a new field checkpoint filename
        vtu_filename = trim(prefix)
        if(size(state) > 1) vtu_filename = trim(vtu_filename) // "_" // trim(state(i)%name)
        vtu_filename = trim(vtu_filename) // "_" // trim(mesh%name)
        if(present(cp_no)) vtu_filename = trim(vtu_filename) // "_" // int2str(cp_no)
        if(present_and_nonempty(postfix)) vtu_filename = trim(vtu_filename) // "_" // trim(postfix)
        if(nparts > 1) then
          vtu_filename = trim(vtu_filename) // ".pvtu"
        else
          vtu_filename = trim(vtu_filename) // ".vtu"
        end if

        if(associated(state(i)%scalar_fields)) then
          allocate(ps_fields_on_mesh(size(state(i)%scalar_fields)))
        else
          allocate(ps_fields_on_mesh(0))
        end if
        if(associated(state(i)%vector_fields)) then
          allocate(pv_fields_on_mesh(size(state(i)%vector_fields)))
        else
          allocate(pv_fields_on_mesh(0))
        end if
        if(associated(state(i)%tensor_fields)) then
          allocate(pt_fields_on_mesh(size(state(i)%tensor_fields)))
        else
          allocate(pt_fields_on_mesh(0))
        end if
        n_ps_fields_on_mesh = 0
        n_pv_fields_on_mesh = 0
        n_pt_fields_on_mesh = 0

        if(associated(state(i)%scalar_fields)) then
          do k = 1, size(state(i)%scalar_fields)
            s_field => state(i)%scalar_fields(k)%ptr
            ! If the mesh names match
            if(trim(s_field%mesh%name) == trim(mesh%name)) then
              ! and the meshes are the same
              if(s_field%mesh == mesh) then
                ! and either the field is prognostic, prescribed and interpolated, or diagnostic and checkpointed
                if(have_option(trim(s_field%option_path) // "/prognostic") &
                    & .or. (have_option(trim(s_field%option_path) // "/prescribed") &
                    & .and. interpolate_field(s_field)) & 
                    & .or. have_option(trim(s_field%option_path) // "/diagnostic/output/checkpoint")) then
                  ! but not aliased
                  if(.not. aliased(s_field)) then

                    if(have_option(trim(complete_field_path(s_field%option_path)) // "/exclude_from_checkpointing")) cycle
                    ! needs_initial_mesh indicates the field is from_file (i.e. we're dealing with a checkpoint)
                    if(present_and_true(keep_initial_data) .and. (.not. needs_initial_mesh(s_field) .and. .not. have_option(trim(s_field%option_path) // "/diagnostic/output/checkpoint"))) cycle

                    ewrite(2, *) "Checkpointing scalar field " // trim(s_field%name) // " in state " // trim(state(i)%name) &
                        & // "on the " // trim(mesh%name)

                    n_ps_fields_on_mesh = n_ps_fields_on_mesh + 1
                    ps_fields_on_mesh(n_ps_fields_on_mesh) = s_field

                    if(have_option(trim(s_field%option_path) // "/prognostic")) then
                      ewrite(2, *) "Updating initial conditions for " // trim(s_field%name)
                      call update_initial_condition_options(trim(s_field%option_path), trim(vtu_filename), "vtu")
                    else if (have_option(trim(s_field%option_path) // "/prescribed").and. &
                        & interpolate_field(s_field)) then
                      ewrite(2, *) "Updating values for " // trim(s_field%name)
                      call update_value_options(trim(s_field%option_path), trim(vtu_filename), "vtu")
                    else if (have_option(trim(s_field%option_path) // "/diagnostic/output/checkpoint").and. &
	                      & interpolate_field(s_field)) then
                      ewrite(2, *) "... diagnostic field"
                    else
                      FLAbort("Can only checkpoint prognostic or prescribed (with interpolation options) fields.")
                    end if
                  end if
                end if
              end if
            end if
          end do
        end if

        if(associated(state(i)%vector_fields)) then
          do k = 1, size(state(i)%vector_fields)
            v_field => state(i)%vector_fields(k)%ptr
            ! If the mesh names match
            if(trim(v_field%mesh%name) == trim(mesh%name)) then
              ! and the meshes are the same
              if(v_field%mesh == mesh) then
                ! and either the field is prognostic, prescribed and interpolated, or diagnostic and checkpointed
                if(have_option(trim(v_field%option_path) // "/prognostic") &
                    & .or. (have_option(trim(v_field%option_path) // "/prescribed") &
                    & .and. interpolate_field(v_field)) &
                    & .or. have_option(trim(v_field%option_path) // "/diagnostic/output/checkpoint")) then
                  ! but not aliased
                  if(.not. aliased(v_field)) then

                    if(have_option(trim(complete_field_path(v_field%option_path)) // "/exclude_from_checkpointing")) cycle
                    ! needs_initial_mesh indicates the field is from_file (i.e. we're dealing with a checkpoint)
                    if(present_and_true(keep_initial_data) .and. (.not. needs_initial_mesh(v_field) .and. .not. have_option(trim(v_field%option_path) // "/diagnostic/output/checkpoint"))) cycle

                    ewrite(2, *) "Checkpointing vector field " // trim(v_field%name) // " in state " // trim(state(i)%name) &
                        & // "on the " // trim(mesh%name)

                    n_pv_fields_on_mesh = n_pv_fields_on_mesh + 1
                    pv_fields_on_mesh(n_pv_fields_on_mesh) = v_field

                    if(have_option(trim(v_field%option_path) // "/prognostic")) then
                      ewrite(2, *) "Updating initial conditions for " // trim(v_field%name)
                      call update_initial_condition_options(trim(v_field%option_path), trim(vtu_filename), "vtu")
                    else if (have_option(trim(v_field%option_path) // "/prescribed").and. &
                        & interpolate_field(v_field)) then
                      ewrite(2, *) "Updating values for " // trim(v_field%name)
                      call update_value_options(trim(v_field%option_path), trim(vtu_filename), "vtu")
                    else if (have_option(trim(v_field%option_path) // "/diagnostic/output/checkpoint").and. &
	                      & interpolate_field(v_field)) then
                      ewrite(2, *) "... diagnostic field"
                    else
                      FLAbort("Can only checkpoint prognostic or prescribed (with interpolation options) fields.")
                    end if
                  end if
                end if
              end if
            end if
          end do
        end if

        if(associated(state(i)%tensor_fields)) then
          do k = 1, size(state(i)%tensor_fields)
            t_field => state(i)%tensor_fields(k)%ptr
            ! If the mesh names match
            if(trim(t_field%mesh%name) == trim(mesh%name)) then
              ! and the meshes are the same
              if(t_field%mesh == mesh) then
                ! and either the field is prognostic, prescribed and interpolated, or diagnostic and checkpointed
                if(have_option(trim(t_field%option_path) // "/prognostic") &
                    & .or. (have_option(trim(t_field%option_path) // "/prescribed") &
                    & .and. interpolate_field(t_field)) &
                    & .or. have_option(trim(t_field%option_path) // "/diagnostic/output/checkpoint")) then
                  ! but not aliased
                  if(.not. aliased(t_field)) then

                    if(have_option(trim(complete_field_path(t_field%option_path)) // "/exclude_from_checkpointing")) cycle
                    ! needs_initial_mesh indicates the field is from_file (i.e. we're dealing with a checkpoint)
                    if(present_and_true(keep_initial_data) .and. (.not. needs_initial_mesh(t_field) .and. .not. have_option(trim(t_field%option_path) // "/diagnostic/output/checkpoint"))) cycle

                    ewrite(2, *) "Checkpointing tensor field " // trim(t_field%name) // " in state " // trim(state(i)%name) &
                        & // "on the " // trim(mesh%name)

                    n_pt_fields_on_mesh = n_pt_fields_on_mesh + 1
                    pt_fields_on_mesh(n_pt_fields_on_mesh) = t_field

                    if(have_option(trim(t_field%option_path) // "/prognostic")) then
                      ewrite(2, *) "Updating initial conditions for " // trim(t_field%name)
                      call update_initial_condition_options(trim(t_field%option_path), trim(vtu_filename), "vtu")
                    else if (have_option(trim(t_field%option_path) // "/prescribed").and. &
                        & interpolate_field(t_field)) then
                      ewrite(2, *) "Updating values for " // trim(t_field%name)
                      call update_value_options(trim(t_field%option_path), trim(vtu_filename), "vtu")
                    else if (have_option(trim(t_field%option_path) // "/diagnostic/output/checkpoint").and. &
                        & interpolate_field(t_field)) then
                      ewrite(2, *) "... diagnostic field"
                    else
                      FLAbort("Can only checkpoint prognostic or prescribed (with interpolation options) fields.")
                    end if
                  end if
                end if
              end if
            end if
          end do
        end if
        if(n_ps_fields_on_mesh + n_pv_fields_on_mesh + n_pt_fields_on_mesh > 0) then
          call vtk_write_fields(vtu_filename, position = positions, model = mesh, &
            & sfields = ps_fields_on_mesh(:n_ps_fields_on_mesh), vfields = pv_fields_on_mesh(:n_pv_fields_on_mesh), &
            & tfields = pt_fields_on_mesh(:n_pt_fields_on_mesh), number_of_partitions=number_of_partitions, stat=stat)
        end if

        deallocate(ps_fields_on_mesh)
        deallocate(pv_fields_on_mesh)
        deallocate(pt_fields_on_mesh)
      end do
    end do

  contains

    subroutine update_initial_condition_options(path, filename, format)
      !!< Updates the initial condition options for a prognostic field with
      !!< options path path

      character(len = *), intent(in) :: path
      character(len = *), intent(in) :: filename
      character(len = *), intent(in) :: format

      integer :: stat, ic, nics

      nics = option_count(trim(path) // "/prognostic/initial_condition")
      do ic = 0, nics-1  ! do while seemed to break, don't know why
        call delete_option(trim(path) // "/prognostic/initial_condition[" // int2str(0) // "]")
      end do
      call set_option_attribute(trim(path) // "/prognostic/initial_condition::WholeMesh/from_file/file_name", trim(filename), stat)
      if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
        FLAbort("Failed to set initial condition filename when checkpointing field with option path " // trim(path))
      end if
      call set_option_attribute(trim(path) // "/prognostic/initial_condition::WholeMesh/from_file/format/name", trim(format), stat)
      if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING) then
        FLAbort("Failed to set initial condition format when checkpointing field with option path " // trim(path))
      end if

    end subroutine update_initial_condition_options

    subroutine update_value_options(path, filename, format)
      !!< Updates the value options for a prescribed field with
      !!< options path path

      character(len = *), intent(in) :: path
      character(len = *), intent(in) :: filename
      character(len = *), intent(in) :: format

      integer :: stat, value, nvalues

      nvalues = option_count(trim(path) // "/prescribed/value")
      do value = 0, nvalues-1
        call delete_option(trim(path) // "/prescribed/value[" // int2str(0) // "]")
      end do
      call set_option_attribute(trim(path) // "/prescribed/value::WholeMesh/from_file/file_name", trim(filename), stat)
      if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
        FLAbort("Failed to set value filename when checkpointing field with option path " // trim(path))
      end if
      call set_option_attribute(trim(path) // "/prescribed/value::WholeMesh/from_file/format/name", trim(format), stat)
      if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING) then
        FLAbort("Failed to set value format when checkpointing field with option path " // trim(path))
      end if

    end subroutine update_value_options

  end subroutine checkpoint_fields
  
  subroutine checkpoint_options(prefix, postfix, cp_no, protect_simulation_name)
    !!< Checkpoint the entire options tree. Outputs to an FLML file with name:
    !!<   [prefix][_cp_no][_postfix].flml
    !!< where cp_no is optional.
    !!< Note that the simulation name in the checkpointed FLML file has, by
    !!< default, "_checkpoint" appended to it.

    character(len = *), intent(in) :: prefix
    character(len = *), intent(in) :: postfix
    integer, optional, intent(in) :: cp_no
    !! If present and .false., do not protect the simulation_name when
    !! checkpointing the options tree
    logical, optional, intent(in) :: protect_simulation_name

    character(len = OPTION_PATH_LEN) :: simulation_name, options_file_filename
    logical :: lprotect_simulation_name

    ewrite(2, *) "Checkpointing options tree"

    assert(len_trim(prefix) > 0)

    lprotect_simulation_name = .not. present_and_false(protect_simulation_name)

    if(lprotect_simulation_name) then
      call get_option("/simulation_name", simulation_name)
      
      ! Temporarily change the simulation name. This is used to protect
      ! checkpointed files (as these will necessarily themselves have
      ! checkpointing enabled).
      call set_option("/simulation_name", trim(simulation_name) // "_checkpoint")
    end if

    ! Construct a new options file filename
    options_file_filename = trim(prefix)
    if(present(cp_no)) options_file_filename = trim(options_file_filename) // "_" // int2str(cp_no)
    if(present_and_nonempty(postfix)) options_file_filename = trim(options_file_filename) // "_" // trim(postfix)
    options_file_filename = trim(options_file_filename) //  ".flml"

    call write_options(options_file_filename)

    if(lprotect_simulation_name) then
      ! Revert the simulation name
      call set_option("/simulation_name", trim(simulation_name))
    end if

  end subroutine checkpoint_options

  subroutine checkpoint_check_options
    !!< Check checkpointing related options

    integer :: cp_period, stat

    if(.not. have_option("/io/checkpointing")) then
      ! Nothing to check
      return
    end if

    ewrite(2, *) "Checking checkpointing options"

#ifndef HAVE_VTK
    ewrite(0, *) "Warning: Checkpointing is enabled, but Fluidity has been compiled without VTK support"
#endif

    call get_option("/io/checkpointing/checkpoint_period_in_dumps", cp_period, stat)
    if(stat /= 0) then
      FLExit("Checkpoint period (in dumps) required for checkpointing")
    end if
    if(cp_period < 0) then
      FLExit("Checkpoint period (in dumps) cannot be negative")
    end if

    ewrite(2, *) "Finished checking checkpointing options"

  end subroutine checkpoint_check_options

end module checkpoint
