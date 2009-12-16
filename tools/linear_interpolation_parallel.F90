#include "fdebug.h"

program linear_interpolation_parallel

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
  type(vector_field), pointer :: positions
  integer :: dim, nprocs, nfiles
  integer :: adapt_iteration, no_adapt_iterations
  real, dimension(:, :), allocatable :: my_bbox
  integer, dimension(:), allocatable :: intersecting_subdomains
  integer :: intersecting_subdomain, file
  type(state_type) :: my_state, vtk_state, interpolation_input, interpolation_output
  type(vector_field), pointer :: vtk_positions
  type(state_type) :: adapting_state
  integer :: rank
  integer :: xpct
  type(scalar_field) :: target_sfield, source_sfield
  type(vector_field) :: target_vfield, source_vfield
  integer :: field, mesh
  character(len=255) :: rank_str
  integer :: rank_digits

  call mpi_init(ierr)
  call python_init

  rank = getrank()
  ! printing out stuff is infuriating in fortran. I don't want 30 spaces ..
  rank_digits = max(floor(log10(float(rank))) + 1, 1)
  rank_str(1:rank_digits) = int2str(rank)

  current_debug_level = 0
  nprocs = getnprocs()

  call parse_args(initial_file, files)
  call vtk_read_state(subdomain_name(initial_file, rank), my_state)
  positions => extract_vector_field(my_state, "Coordinate")

  dim = positions%dim
  nfiles = size(files)

  allocate(bboxes(nfiles, dim, 2, nprocs))
  call compute_bboxes(files, bboxes)
  allocate(my_bbox(dim, 2))
  my_bbox = compute_bbox(positions)

  do file=1,nfiles
    ewrite(0,*) "  rank " // rank_str(1:rank_digits) // ": processing file " // trim(files(file))
    call get_intersecting_subdomains(bboxes(file, :, :, :), my_bbox, intersecting_subdomains)
    ewrite(0,*) " rank " // rank_str(1:rank_digits) // ": intersecting subdomains: ", intersecting_subdomains
    do i=1,size(intersecting_subdomains)
      intersecting_subdomain = intersecting_subdomains(i)
      ewrite(0,*) "    rank " // rank_str(1:rank_digits) // ": processing subdomain ", intersecting_subdomain, "; intersection ", i, " of ", size(intersecting_subdomains)

      ! So, this subdomain might have some info for us .. 
      call vtk_read_state(trim(subdomain_name(files(file), intersecting_subdomain-1)), vtk_state)
      vtk_positions => extract_vector_field(vtk_state, "Coordinate")
      call insert(interpolation_input, extract_mesh(vtk_state, "Mesh"), "Mesh")
      call insert(interpolation_input, vtk_positions, "Coordinate")

      ! first time through the loop, need to allocate the target fields
      do field=1,scalar_field_count(vtk_state)
        source_sfield = extract_scalar_field(vtk_state, field)
        call insert(interpolation_input, source_sfield, trim(source_sfield%name))
        if (i == 1) then
          call allocate(target_sfield, positions%mesh, trim(source_sfield%name))
          call zero(target_sfield)
          call insert(interpolation_output, target_sfield, trim(source_sfield%name))
          call deallocate(target_sfield)
        end if
      end do

      do field=1,vector_field_count(vtk_state)
        source_vfield = extract_vector_field(vtk_state, field)
        if (trim(source_vfield%name) == "Coordinate") then
          cycle
        end if
        call insert(interpolation_input, source_vfield, trim(source_vfield%name))

        if (i == 1) then
          call allocate(target_vfield, source_vfield%dim, positions%mesh, trim(source_vfield%name))
          call zero(target_vfield)
          call insert(interpolation_output, target_vfield, trim(source_vfield%name))
          call deallocate(target_vfield)
        end if
      end do

      call deallocate(vtk_state)

      if (i == 1) then
        call insert(interpolation_output, positions%mesh, "Mesh")
        call insert(interpolation_output, positions, "Coordinate")
      end if

      ewrite(0,*) "    rank " // rank_str(1:rank_digits) // ": interpolating"
      call linear_interpolation(interpolation_input, interpolation_output, different_domains=.true.)
      call deallocate(interpolation_input)

      !call vtk_write_state("interpolation_output", i, state=(/interpolation_output/))

    end do

    do field=1,scalar_field_count(interpolation_output)
      target_sfield = extract_scalar_field(interpolation_output, field)
      target_sfield%name = trim(target_sfield%name) // int2str(file)
      call insert(my_state, target_sfield, trim(target_sfield%name))
    end do

    do field=1,vector_field_count(interpolation_output)
      target_vfield = extract_vector_field(interpolation_output, field)
      if (trim(target_vfield%name) == "Coordinate") then
        cycle
      end if
      target_vfield%name = trim(target_vfield%name) // int2str(file)
      call insert(my_state, target_vfield, trim(target_vfield%name))
    end do

    call deallocate(interpolation_output)
    deallocate(intersecting_subdomains)
  end do

  call insert(my_state, positions%mesh, "Mesh")
  call vtk_write_state("linear_interpolation_parallel", model="Mesh", state=(/my_state/))
  call deallocate(my_state)

  deallocate(bboxes)
  deallocate(my_bbox)
  deallocate(files)
  call print_references(-1)
  call mpi_finalize(ierr)

  contains

  function subdomain_name(pvtu_name, subdomain) result(name)
    character(len=255), intent(in) :: pvtu_name
    integer, intent(in) :: subdomain
    character(len=255) :: name

    name = pvtu_name(1:len_trim(pvtu_name) - len(".pvtu")) // "_" // int2str(subdomain) // ".vtu"
  end function subdomain_name

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

    nfiles = size(files)
    rank = getrank()

    do i=1,nfiles
      call vtk_read_state(trim(subdomain_name(files(i), rank)), state)
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
  
  subroutine parse_args(initial_file, files)
    character(len=255), intent(out) :: initial_file
    character(len=255), dimension(:), allocatable, intent(out) :: files
    integer :: argc, status, i

    argc = command_argument_count()
    if (argc < 2) then
      write(0,*) "Usage: linear_interpolation_parallel initial_file input_mesh [input_mesh ...]"
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

    allocate(files(argc-1))
    do i=2,argc
      call get_command_argument(i, value=files(i-1), status=status)
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
    
end program linear_interpolation_parallel
