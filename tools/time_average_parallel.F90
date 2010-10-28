#include "fdebug.h"

program time_average_parallel

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
  use mesh_files
  use reference_counting
  implicit none

  character(len=255) :: initial_file
  character(len=255) :: fieldname
  integer :: i, j, rank
  integer :: ierr

  type(vector_field), pointer :: positions
  integer :: dim, nprocs, nfiles
  type(state_type) :: state
  type(scalar_field), pointer :: sfield
  type(scalar_field) :: savg
  type(vector_field), pointer :: vfield
  type(vector_field) :: vavg
  integer :: stat, argc, status

  call mpi_init(ierr)
  call python_init
  rank = getrank()

  current_debug_level = 0
  nprocs = getnprocs()
  write(0,*) "process count: ", nprocs
  write(0,*) "rank: ", rank

  argc = command_argument_count()
  if (argc == 0) then
    write(0,*) "Usage: time_average_parallel initial_file FieldName"
    stop
  end if

  write(0,*) "argument count: ", argc
  call get_command_argument(1, value=initial_file, status=status)
  write(0,*) "initial_file: ", initial_file

  call vtk_read_state(subdomain_name(initial_file, rank), state)
  positions => extract_vector_field(state, "Coordinate")
  !call parse_args(initial_file, fieldname)

  ! Modification to parse multiple fields
  do i = 2, argc
    call get_command_argument(i, value=fieldname, status=status)
    write(0,*) "fieldname: ", fieldname

    sfield => extract_scalar_field(state, trim(fieldname) // "1", stat=stat)
    if (stat == 0) then
      call allocate(savg, sfield%mesh, trim(fieldname) // "Average")
      call zero(savg)
      write(0,*) "field: ", trim(fieldname)

      j = 1
      do 
        sfield => extract_scalar_field(state, trim(fieldname) // int2str(i), stat=stat)
        if (stat /= 0) then
          exit
        end if
        call addto(savg, sfield)
        j = j + 1
      end do

      call scale(savg, 1.0/(i-1))
      call insert(state, savg, trim(fieldname) // "Average")
      call deallocate(savg)
    else
      vfield => extract_vector_field(state, trim(fieldname) // "1")
      call allocate(vavg, vfield%dim, vfield%mesh, trim(fieldname) // "Average")
      call zero(vavg)
      write(0,*) "field, dim: ", trim(fieldname), vfield%dim

      j = 1
      do 
        vfield => extract_vector_field(state, trim(fieldname) // int2str(i), stat=stat)
        if (stat /= 0) then
          exit
        end if
        call addto(vavg, vfield)
        j = j + 1
      end do

      call scale(vavg, 1.0/(i-1))
      call insert(state, vavg, trim(fieldname) // "Average")
      call deallocate(vavg)
    end if

    call insert(state, positions%mesh, "Mesh")
    call vtk_write_state("time_average_parallel", model="Mesh", state=(/state/))
    call deallocate(state)

  end do

  call print_references(-1)
  call mpi_finalize(ierr)

  contains

  function subdomain_name(pvtu_name, subdomain) result(name)
    character(len=255), intent(in) :: pvtu_name
    integer, intent(in) :: subdomain
    character(len=255) :: name

    name = pvtu_name(1:len_trim(pvtu_name) - len(".pvtu")) // "_" // int2str(subdomain) // ".vtu"
  end function subdomain_name

  ! No longer needed
  subroutine parse_args(initial_file, fieldname)
    character(len=255), intent(out) :: initial_file
    character(len=255), intent(out) :: fieldname
    integer :: argc, status

    argc = command_argument_count()
    write(0,*) "argument count: ", argc
    call get_command_argument(1, value=initial_file, status=status)
    write(0,*) "initial_file: ", initial_file
    call get_command_argument(2, value=fieldname, status=status)
    write(0,*) "fieldname: ", fieldname

    if (argc /= 2) then
      write(0,*) "Usage: time_average_parallel initial_file FieldName"
      stop
    end if

    call get_command_argument(1, value=initial_file, status=status)
    select case(status)
    case(:-1)
       write(0,*) "Warning: truncating filename"
    end select

    call get_command_argument(2, value=fieldname, status=status)
    select case(status)
    case(:-1)
       write(0,*) "Warning: truncating filename"
    end select
  end subroutine parse_args

end program time_average_parallel
