#include "fdebug.h"

module interpolation_ensemble_state_on_supermesh
  use fields
  use state_module
  use vtk_interfaces
  use pseudo_supermesh
  use interpolation_manager
  use spud
  use unify_meshes_module
  use mesh_files
  use populate_state_module
  use populate_sub_state_module
  use Profiler
  use elements
  use state_module
  use FLDebug
!  use spud
  use mesh_files
  use vtk_cache_module
  use global_parameters, only: OPTION_PATH_LEN, is_active_process, pi, &
    no_active_processes, topology_mesh_name, adaptivity_mesh_name, &
    periodic_boundary_option_path, domain_bbox, domain_volume
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
  !  use hadapt_extrude_radially
  use limit_metric_module
  use initialise_fields_module
  use transform_elements
  use parallel_tools
  use boundary_conditions_from_options
  use nemo_states_module
  use data_structures
  use fields_halos
  use read_triangle
  !  use sediment, only: get_nSediments, get_sediment_name
  use futils, only: int2str
  use conformity_measurement
  use interpolation_module
  use merge_tensors
  use adapt_state_module
 implicit none

contains

  subroutine interpolation_ensembles_on_supermesh(ensemble_state_new, ensemble_state_old, no)

  type(state_type), dimension(:), intent(in) :: ensemble_state_old
  type(state_type), dimension(:), intent(inout) :: ensemble_state_new
  integer, intent(in) :: no

  type(state_type) :: tmp_state_1, tmp_state_2

  type(vector_field) :: positions_old,positions_new,out_positions
  !type(vector_field),pointer :: out_positions
  character(len=255) :: filename
  character(len=255), dimension(:), allocatable :: files
  integer :: argc
  integer :: i, status,ifield,nfield,nfield_old,jfield
  type(state_type) :: initial_state
  integer :: ierr
  integer :: stat, mxnods
  integer :: nrens
  integer :: field_count
  type(vector_field) :: velocity
  type(scalar_field) :: pressure, free_surface,time
  type(element_type) :: shape
  integer :: quad_degree

  type(vector_field),pointer :: field
   type(vector_field):: vfield2
   type(scalar_field):: sfield2
   type(tensor_field):: tfield2

   type(vector_field), pointer :: vfield,vfield_old,vfield_new,vfield_tmp
   type(scalar_field), pointer :: sfield,sfield_old,sfield_new
   type(tensor_field), pointer :: tfield,tfield_old,tfield_new


    print*,'In interpolation_ensembles_on_supermesh'
    !nrens = size(ensemble_state_old)    
!    call mpi_init(ierr)
    call set_option('/mesh_adaptivity/hr_adaptivity/maximum_number_of_nodes', mxnods, stat=stat)
    call get_option('/geometry/quadrature/degree', quad_degree)     
  
    print*,'linear_interpolation'
!    do i=1, nrens
!       print*,'****iren=',i
       print*,'interpolation of scalar'
       nfield = scalar_field_count( ensemble_state_old(1) )
       do ifield = 1, nfield 
          sfield_old => extract_scalar_field( ensemble_state_old(1), ifield )
          !    call allocate( sfield_old,  sfield%mesh, trim(sfield%name)   )
          print*,trim(sfield_old%name)
          !    call set(sfield_old, sfield)
          positions_old = extract_vector_field(ensemble_state_old(1), "Coordinate")
          ! new field
          sfield_new => extract_scalar_field( ensemble_state_new(no), ifield )
          !       call allocate( sfield_new,  sfield%mesh, trim(sfield%name)   )
          print*,trim(sfield_new%name)
          !       call set(sfield_new, sfield)
          positions_new = extract_vector_field(ensemble_state_new(no), "Coordinate")
          call linear_interpolation(sfield_old, positions_old, sfield_new,positions_new)
       end do
       
       
       print*,'interpolation of vector'
       nfield = vector_field_count( ensemble_state_new(no) )
       do ifield = 1, nfield-1 
          print*,'ifield',ifield
          vfield_new => extract_vector_field( ensemble_state_new(no), ifield )
          do jfield = 1,vector_field_count( ensemble_state_old(1) )
             vfield_tmp => extract_vector_field( ensemble_state_old(1), jfield )
             print*,trim(vfield_tmp%name)
             print*,'11tmp'
             if(trim(vfield_tmp%name)==trim(vfield_new%name)) then
                vfield_old => extract_vector_field( ensemble_state_old(1), jfield )
                exit
             endif
          enddo
          !          call allocate( vfield_old, vfield%dim, vfield%mesh, trim(vfield%name)   )
          print*,trim(vfield_old%name)
          print*,'11'
          !          call set(vfield_old, vfield)
          positions_old = extract_vector_field(ensemble_state_old(1), "Coordinate")
          ! new field
          !  call allocate( vfield_new, vfield%dim,vfield%mesh, trim(vfield%name)   )
          print*,trim(vfield_new%name)
          print*,'22'
          !  call set(vfield_new, vfield)
          positions_new = extract_vector_field(ensemble_state_new(no), "Coordinate")
          call linear_interpolation(vfield_old, positions_old, vfield_new,positions_new)
       end do
       
       print*,'interpolation of tensor'
       nfield = tensor_field_count( ensemble_state_new(no) )
       do ifield = 1, nfield 
          tfield_old => extract_tensor_field( ensemble_state_old(1), ifield )
          !    call allocate( tfield_old,  tfield%mesh, trim(tfield%name)   )
          !    call set(tfield_old, tfield)
          positions_old = extract_vector_field(ensemble_state_old(1), "Coordinate")
          ! new field
          tfield_new => extract_tensor_field( ensemble_state_new(no), ifield )
          !   call allocate( tfield_new,  tfield%mesh, trim(tfield%name)   )
          !   call set(tfield_new, tfield)
          positions_new = extract_vector_field(ensemble_state_new(no), "Coordinate")
          call linear_interpolation(tfield_old, positions_old, tfield_new,positions_new)
          
       end do
       print*,'%%%%%%%%%%%%%%%%%%%%%%%%%i=',i
       !       call linear_interpolation(ensemble_state_old(i), ensemble_state_new(i))
       pressure = extract_scalar_field(ensemble_state_new(no), "Pressure")
       print*,pressure%val(1:20)
       print*,'pressurenew',node_count(pressure%mesh)
       velocity = extract_vector_field(ensemble_state_new(no), "Velocity")
       print*,'velocitynew',velocity%val(1,1:20)
 !   enddo
    

  end subroutine interpolation_ensembles_on_supermesh


    !!femtools/Fields_Calculations.F90
    !function norm2_difference_single(fieldA, positionsA, fieldB, positionsB) result(norm)
 

  subroutine compute_pseudo_supermesh_from_mesh(super_positions, no_its, mxnods)
    !!< snapshots is a list of VTUs containing the meshes we want
    !!< to merge.
    !!< starting_positions is the initial mesh + positions to interpolate the
    !!< metric tensor fields describing the snapshot meshes onto.
    !!< super_positions is the output -- a positions field on a mesh.
!    character(len=255), dimension(:), intent(in) :: snapshots
    type(vector_field), dimension(:), allocatable  :: starting_positions
!    type(vector_field), dimension(:), intent(in) :: starting_positions
    type(vector_field), intent(out) :: super_positions
    integer, intent(in), optional :: no_its, mxnods

    integer :: lno_its
    integer :: it, i

    type(mesh_type) :: current_mesh, vtk_mesh
    type(vector_field) :: current_pos, vtk_pos
    type(state_type) :: vtk_state, temp_state
    type(state_type) :: interpolation_input, interpolation_output

    type(tensor_field) :: merged_metric, interpolated_metric
    type(tensor_field) :: vtk_metric
    integer :: quad_degree
    integer :: nrens
    character(len = OPTION_PATH_LEN) :: simulation_name
    character(len=1024) :: filename

      call get_option('/geometry/quadrature/degree', quad_degree)
      call get_option('/simulation_name',simulation_name)

      nrens=2
      allocate(starting_positions(1:nrens))
      do i=1,nrens
         write(filename, '(a, i0, a)') trim(simulation_name)//'_CoordinateMesh_', i,'_checkpoint'
         starting_positions(i)=read_mesh_files(filename=filename, quad_degree=quad_degree,format="gmsh")
      enddo

    if (present(no_its)) then
      lno_its = no_its
    else
      lno_its = 3
    end if

    current_pos = starting_positions(1)
    call incref(starting_positions(1))
    current_mesh = starting_positions(1)%mesh
    call incref(current_mesh)

    do it=1,lno_its
      call allocate(merged_metric, current_mesh, "MergedMetric")
      call zero(merged_metric)
      call allocate(interpolated_metric, current_mesh, "InterpolatedMetric")
      call zero(interpolated_metric)
      call insert(interpolation_output, interpolated_metric, "InterpolatedMetric")
      call insert(interpolation_output, current_mesh, "CoordinateMesh")
      call insert(interpolation_output, current_pos, "Coordinate")
      

      do i=1,size(starting_positions)
        call zero(interpolated_metric)
        !!call vtk_read_state(trim(snapshots(i)), vtk_state)
        !!vtk_mesh = extract_mesh(vtk_state, "Mesh")
        !!vtk_pos  = extract_vector_field(vtk_state, "Coordinate")
        vtk_mesh = starting_positions(i)%mesh
        vtk_pos  = starting_positions(i)

        call allocate(vtk_metric, vtk_mesh, "MeshMetric")
        print*,'before compute_mesh_metric'
        call compute_mesh_metric(vtk_pos, vtk_metric)
        print*,'after compute_mesh_metric'
        call insert(interpolation_input, vtk_metric, "InterpolatedMetric")
        call insert(interpolation_input, vtk_mesh, "CoordinateMesh")
        call insert(interpolation_input, vtk_pos, "Coordinate")
        print*,'vtk_write_state("interpolation_input")'
        !!call vtk_write_state("interpolation_input", i, state=(/interpolation_input/))
        call linear_interpolation(interpolation_input, interpolation_output)
        print*,'vtk_write_state("interpolation_output")'
        !!call vtk_write_state("interpolation_output", i, state=(/interpolation_output/))
        call merge_tensor_fields(merged_metric, interpolated_metric)
        call deallocate(vtk_metric)
        call deallocate(interpolation_input)
        call deallocate(vtk_state)
      end do

      call deallocate(interpolated_metric)
      call deallocate(interpolation_output)

      call insert(temp_state, current_mesh, "CoordinateMesh")
      call insert(temp_state, current_pos, "Coordinate")
      ! Assuming current_mesh had a refcount of one,
      ! it now has a refcount of two.
!      call get_edge_lengths(merged_metric, edgelen)
!      call vtk_write_fields("supermesh_before_adapt", it, current_pos, current_mesh, sfields=(/edgelen/), tfields=(/merged_metric/))
      if (present(mxnods)) then
        call limit_metric(current_pos, merged_metric, min_nodes=1, max_nodes=mxnods)
      end if
      print*,'before adapt_state'
      call adapt_state(temp_state, merged_metric)
!      call vtk_write_state("supermesh_after_adapt", it, state=(/temp_state/))
      call deallocate(merged_metric)
!      call deallocate(edgelen)
      ! Now it has a refcount of one, as adapt_state
      ! has destroyed the old one and created a new mesh
      ! with refcount one.

      ! We're finished with the current_mesh, so let it be
      ! deallocated if no one else is using it.
      call deallocate(current_mesh)
      call deallocate(current_pos)

      current_mesh = extract_mesh(temp_state, "CoordinateMesh")
      current_pos  = extract_vector_field(temp_state, "Coordinate")
      call incref(current_mesh)
      call incref(current_pos)
      call deallocate(temp_state)
    end do

    call deallocate(current_mesh)
    super_positions = current_pos
  end subroutine compute_pseudo_supermesh_from_mesh

  subroutine allocate_new_ensemble_state(out_positions, ensemble_state_new)
    type(vector_field), intent(in) :: out_positions
    type(state_type), dimension(:), intent(inout) :: ensemble_state_new

    integer :: nrens,nfield, i,ifield
    type(vector_field), pointer:: vfield
    type(scalar_field), pointer:: sfield
    type(tensor_field), pointer:: tfield
    type(vector_field) :: vfield2
    type(scalar_field) :: sfield2
    type(tensor_field) :: tfield2

    nrens = size(ensemble_state_new)
    
    do i = 1, nrens
       call nullify(ensemble_state_new(i)) 
       call set_option_path(ensemble_state_new(i), "/material_phase["//int2str(i-1)//"]")
    end do

  
    do i=1, size(ensemble_state_new)       
       call insert(ensemble_state_new(i), out_positions, "Coordinate")
       call insert(ensemble_state_new(i), out_positions%mesh, "CoordinateMesh")
    end do
    
   
    call insert_derived_meshes( ensemble_state_new)
    call allocate_and_insert_fields(ensemble_state_new)
    
    
    do i = 2, size(ensemble_state_new)
       nfield = scalar_field_count( ensemble_state_new(1) )
       do ifield = 1, nfield 
          sfield => extract_scalar_field( ensemble_state_new(1), ifield )
          call allocate( sfield2,  sfield%mesh, trim(sfield%name)   )
          call set(sfield2, sfield)
          call insert(ensemble_state_new(i), sfield2, trim(sfield%name))
       end do
       
       
       nfield = vector_field_count( ensemble_state_new(1) )
       print*,'&&&&vector',nfield
       do ifield = 1, nfield 
          vfield => extract_vector_field( ensemble_state_new(1), ifield )
          call allocate( vfield2,  vfield%dim, vfield%mesh, trim(vfield%name)   )
          call set(vfield2, vfield)
          call insert(ensemble_state_new(i), vfield2, trim(vfield%name))
       end do
       
       nfield = tensor_field_count( ensemble_state_new(1) )
       do ifield = 1, nfield 
          tfield => extract_tensor_field( ensemble_state_new(1), ifield )
          call allocate( tfield2,  tfield%mesh, trim(tfield%name)    )
          call set(tfield2, tfield)
          call insert(ensemble_state_new(i), tfield2, trim(tfield%name))
       end do
    end do
  end subroutine allocate_new_ensemble_state
 
end module Interpolation_ensemble_state_on_supermesh
