#include "fdebug.h"

program pseudo_supermesh_parallel

  use fields
  use state_module
  use vtk_interfaces
  use spud
  use parallel_tools
  use mpi
  use interpolation_module
  use adapt_state_unittest_module, adapt_state => adapt_state_unittest
  use conformity_measurement
  use merge_tensors
  use limit_metric_module
  use unittest_tools
  use metric_tools
  use global_parameters, only: current_debug_level
  use read_triangle
  use reference_counting
  implicit none

  character(len=255) :: initial_file
  character(len=255), dimension(:), allocatable :: files
  integer :: i
  integer :: ierr
  integer :: mxnods

  ! nfiles x dim x 2 x nprocs
  real, dimension(:, :, :, :), allocatable :: bboxes
  type(vector_field) :: positions
  integer :: dim, nprocs, nfiles
  integer :: adapt_iteration, no_adapt_iterations
  real, dimension(:, :), allocatable :: my_bbox
  integer, dimension(:), allocatable :: intersecting_subdomains
  integer :: intersecting_subdomain, file
  type(tensor_field) :: merged_metric, interpolated_metric, vtk_metric
  type(state_type) :: vtk_state, interpolation_input, interpolation_output
  type(vector_field), pointer :: vtk_positions
  type(state_type) :: adapting_state
  integer :: rank
  integer :: xpct

  call mpi_init(ierr)
  call python_init

  no_adapt_iterations = 3
  rank = getrank()

  current_debug_level = 0
  nprocs = getnprocs()

  call parse_args(initial_file, mxnods, files)
  positions = read_triangle_files(trim(initial_file), quad_degree=4)
  call halo_business(nprocs, positions, initial_file)

  dim = positions%dim
  nfiles = size(files)

  allocate(bboxes(nfiles, dim, 2, nprocs))
  call compute_bboxes(files, bboxes)
  allocate(my_bbox(dim, 2))
  my_bbox = compute_bbox(positions)
!  call print_references(-1)

  do adapt_iteration=1,no_adapt_iterations
    ewrite(0,*) "rank: ", rank, ": starting adapt iteration ", adapt_iteration
    call allocate(merged_metric, positions%mesh, "MergedMetric")
    call allocate(interpolated_metric, positions%mesh, "Metric")
    call zero(merged_metric)

    do file=1,nfiles
      ewrite(0,*) "  rank: ", rank, ": processing file " // trim(files(file))
      call get_intersecting_subdomains(bboxes(file, :, :, :), my_bbox, intersecting_subdomains)
      ewrite(0,*) " rank: ", rank, ": intersecting subdomains: ", intersecting_subdomains
      do i=1,size(intersecting_subdomains)
        call zero(interpolated_metric)
        intersecting_subdomain = intersecting_subdomains(i)
        ewrite(0,*) "    rank: ", rank, ": processing subdomain ", intersecting_subdomain, "; intersection ", i, " of ", size(intersecting_subdomains)

        ! So, this subdomain might have some info for us .. 
        ! Compute its metric, interpolate, and merge.
        call vtk_read_state(trim(subdomain_name(files(file), intersecting_subdomain-1)), vtk_state)
        vtk_positions => extract_vector_field(vtk_state, "Coordinate")

        call insert(interpolation_input, extract_mesh(vtk_state, "Mesh"), "Mesh")
        call insert(interpolation_input, vtk_positions, "Coordinate")
        ewrite(0,*) "    rank: ", rank, ": computing subdomain mesh metric"
        call compute_mesh_metric(vtk_positions, vtk_metric)
        call insert(interpolation_input, vtk_metric, "Metric")
        call deallocate(vtk_metric)
        call deallocate(vtk_state)

        call insert(interpolation_output, positions%mesh, "Mesh")
        call insert(interpolation_output, positions, "Coordinate")
        call insert(interpolation_output, interpolated_metric, "Metric")

        ewrite(0,*) "    rank: ", rank, ": interpolating subdomain mesh metric"
        call linear_interpolation(interpolation_input, interpolation_output, different_domains=.true.)
        ewrite(0,*) "    rank: ", rank, ": merging subdomain mesh metric"
        call merge_tensor_fields(merged_metric, interpolated_metric)
        xpct = expected_elements(positions, merged_metric)
        ewrite(0,*) "    rank: ", rank, ": expected number of elements: ", xpct
        call deallocate(interpolation_input)
        call deallocate(interpolation_output)
      end do

      deallocate(intersecting_subdomains)
    end do

    call deallocate(interpolated_metric)

    ! Limit metric to mxnods
    ewrite(0,*) "rank: ", rank, ": limiting maximum number of nodes"
    call limit_metric(positions, merged_metric, min_nodes=1, max_nodes=mxnods)

    ! Adapt here
    call insert(adapting_state, positions%mesh, "CoordinateMesh")
    call insert(adapting_state, positions, "Coordinate")
    xpct = expected_elements(positions, merged_metric)
    ewrite(0,*) "rank: ", rank, ": expected number of elements: ", xpct
    !call vtk_write_fields("before_adapt", adapt_iteration, position=positions, model=positions%mesh, tfields=(/merged_metric/))
    ewrite(0,*) "rank: ", rank, ": entering adapt"

    current_debug_level = 3
    call adapt_state(adapting_state, merged_metric)
    current_debug_level = 0
    call deallocate(merged_metric)
    call deallocate(positions)

    positions = extract_vector_field(adapting_state, "Coordinate")
    call incref(positions)
    call deallocate(adapting_state)

    !call vtk_write_fields("after_adapt", adapt_iteration, position=positions, model=positions%mesh)
    my_bbox = compute_bbox(positions)
!    call print_references(-1)
  end do

  call vtk_write_fields("pseudo_supermesh", position=positions, model=positions%mesh)

  call deallocate(positions)
  deallocate(bboxes)
  deallocate(my_bbox)
  deallocate(files)
  call print_references(-1)
  call mpi_finalize(ierr)

  contains

  subroutine halo_business(nprocs, positions, filename)
    use halos
    use global_parameters, only: halo_tag, halo_tag_p
    use flcomms_io
    integer, intent(in) :: nprocs
    type(vector_field), intent(inout) :: positions
    character(len=255), intent(in) :: filename

    if (nprocs > 1) then
      call read_halos(filename(1:len_trim(filename)))

      allocate(positions%mesh%halos(2))
      call import_halo(halo_tag, positions%mesh%halos(1))
      call import_halo(halo_tag_p, positions%mesh%halos(2))
      allocate(positions%mesh%element_halos(2))
      call derive_element_halo_from_node_halo(positions%mesh, &
        & ordering_scheme = HALO_ORDER_TRAILING_RECEIVES)       
    end if
  end subroutine halo_business

  function subdomain_name(pvtu_name, subdomain) result(name)
    character(len=255), intent(in) :: pvtu_name
    integer, intent(in) :: subdomain
    character(len=255) :: name

    name = pvtu_name(1:len_trim(pvtu_name) - len(".pvtu")) // "_" // int2str(subdomain) // ".vtu"
  end function subdomain_name

  subroutine set_metric_to_value(metric, value)
    type(tensor_field), intent(inout) :: metric
    real, intent(in) :: value
    integer :: node, i
    real :: eigval
    real, dimension(metric%dim, metric%dim) :: met

    eigval = 1/value**2
    met = 0.0
    do i=1,metric%dim
      met(i, i) = eigval
    end do

    do node=1,node_count(metric)
      call set(metric, node, met)
    end do
  end subroutine set_metric_to_value

  subroutine get_intersecting_subdomains(bboxes, my_bbox, intersecting_subdomains)
    use intersection_finder_module
    real, dimension(:, :, :), intent(in) :: bboxes
    real, dimension(:, :), intent(in) :: my_bbox
    integer, dimension(:), allocatable, intent(out) :: intersecting_subdomains

    integer :: nprocs
    integer :: i
    integer, dimension(:), allocatable :: tmp_array
    integer :: array_sz

    array_sz = 0
    nprocs = size(bboxes, 3)
    allocate(tmp_array(nprocs))

    do i=1,nprocs
      if (bbox_predicate(bboxes(:, :, i), my_bbox)) then
        array_sz = array_sz + 1
        tmp_array(array_sz) = i
      end if
    end do

    allocate(intersecting_subdomains(array_sz))
    intersecting_subdomains(1:array_sz) = tmp_array(1:array_sz)
    deallocate(tmp_array)
  end subroutine get_intersecting_subdomains

  subroutine compute_bboxes(files, bboxes)
    character(len=255), dimension(:), intent(in) :: files
    real, dimension(:, :, :, :), intent(inout) :: bboxes
    real, dimension(size(bboxes, 1), size(bboxes, 2), size(bboxes, 3)) :: my_bboxes
    integer :: rank
    integer :: i, nfiles
    type(state_type) :: state
    type(vector_field), pointer :: positions
    integer :: ierr
    character(len=255) :: subname

    nfiles = size(files)
    rank = getrank()

    do i=1,nfiles
      subname = subdomain_name(files(i), rank)
      call vtk_read_state(trim(subname), state)
      positions => extract_vector_field(state, "Coordinate")
      my_bboxes(i, :, :) = compute_bbox(positions)
      call deallocate(state)
    end do

    call mpi_gather(my_bboxes, size(my_bboxes), getpreal(), &
                    bboxes, size(my_bboxes), getpreal(), &
                    0, MPI_COMM_WORLD, ierr)
    assert(ierr == 0)

    call mpi_bcast(bboxes, size(bboxes), getpreal(), 0, MPI_COMM_WORLD, ierr)
    assert(ierr == 0)
  end subroutine compute_bboxes
  
  subroutine parse_args(initial_file, mxnods, files)
    character(len=255), intent(out) :: initial_file
    character(len=255), dimension(:), allocatable, intent(out) :: files
    character(len=255) :: mxnods_buffer
    integer, intent(out) :: mxnods
    integer :: argc, status, i

    argc = command_argument_count()
    if (argc < 4) then
      write(0,*) "Usage: pseudo_supermesh_parallel initial_file max_nodes input_mesh [input_mesh ...]"
      stop
    end if

    call get_command_argument(1, value=initial_file, status=status)
    select case(status)
    case(1:)
       write(0,*) "Initial VTU filename not found"
       stop
    case(:-1)
       write(0,*) "Warning: truncating filename"
    end select

    call get_command_argument(2, value=mxnods_buffer, status=status)
    read(mxnods_buffer, *, iostat=status) mxnods
    if (status /= 0) then
      write(0,*) "Could not parse maximum number of nodes"
      stop
    end if

    allocate(files(argc-2))
    do i=3,argc
      call get_command_argument(i, value=files(i-2), status=status)
      select case(status)
      case(1:)
         write(0,*) "Initial VTU filename not found"
         stop
      case(:-1)
         write(0,*) "Warning: truncating filename"
      end select
    end do
  end subroutine parse_args

  function compute_bbox(positions) result(bbox)
    type(vector_field), intent(in) :: positions
    real, dimension(positions%dim, 2) :: bbox

    integer :: i

    do i=1,positions%dim
      bbox(i, 1) = minval(positions%val(i)%ptr)
      bbox(i, 2) = maxval(positions%val(i)%ptr)
    end do
  end function compute_bbox
    
end program pseudo_supermesh_parallel
