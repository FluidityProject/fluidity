!    Copyright (C) 2006 Imperial College London and others.
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

module particles

  use fldebug
  use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, &
       & PYTHON_FUNC_LEN, integer_size, real_size
  use state_module
  use elements
  use detector_move_lagrangian
  use detector_data_types
  use detector_tools
  use detector_parallel
  use fields
  use parallel_tools
  use pickers
  use mpi_interfaces
  use futils, only: int2str, free_unit
  use spud
  use parallel_fields
  use diagnostic_variables, only: field_tag, initialise_constant_diagnostics

  implicit none

  private

  public :: initialise_particles, move_particles, write_particles_loop, destroy_particles, &
            update_particle_attributes

  type(detector_linked_list), allocatable, dimension(:), save :: particle_lists

contains

  subroutine initialise_particles(filename,state)
    !!Initialise particles and set up particle file headers (per particle array)
    character(len = *), intent(in) :: filename
    type(state_type), dimension(:), intent(in) :: state

    character(len=FIELD_NAME_LEN) ::particle_name, funcnam
    character(len=PYTHON_FUNC_LEN) :: func
    character(len = OPTION_PATH_LEN) :: particle_file_filename, particles_cp_filename, name
    character(len = FIELD_NAME_LEN) :: buffer, fmt

    real, allocatable, dimension(:,:) :: coords
    real, allocatable, dimension(:) :: attribute_vals
    real, allocatable, dimension(:) :: packed_buff
    real, allocatable, dimension(:) :: particle_location
    type(vector_field), pointer :: xfield
    real:: current_time

    integer :: str_size
    integer :: i, j, m, n, l
    integer :: python_particles_func, dim, attribute_dims
    integer :: particle_file_unit=0, particle_checkpoint_unit=0
    integer :: column, IERROR, totaldet_global


    logical :: from_checkpoint

    ewrite(2,*), "In initialise_particles"

    !Check whether there are any particles.
    python_particles_func = option_count("/particles/particle_array")

    if (python_particles_func==0) return

    !Set up particle_lists
    allocate(particle_lists(python_particles_func))

    !Allocate parameters
    xfield=>extract_vector_field(state(1), "Coordinate")
    call get_option("/geometry/dimension",dim)
    call get_option("/timestepping/current_time", current_time)
    allocate(particle_location(dim))

    ! If the option
    ! "from_checkpoint_file" exists, it means we are continuing the simulation
    ! after checkpointing and the reading of the particle positions must be
    ! done from a file
    if (have_option("/particles/particle_array/from_checkpoint_file")) then
       from_checkpoint=.true.
    else
       from_checkpoint=.false.
    end if

    do i = 1,python_particles_func
       write(buffer, "(a,i0,a)") "/particles/particle_array[",i-1,"]"
       call get_option(trim(buffer)//"/number_of_particles", j)
       particle_lists(i)%total_num_det=j
       allocate(particle_lists(i)%detector_names(j))
       !Register this I/O list with a global list of detectors/particles
       call register_detector_list(particle_lists(i))

       !Find number of attributes 
       attribute_dims=0
       if (have_option('/particles/particle_array['//int2str(i-1)//']/attributes/attribute')) then
          attribute_dims=option_count('/particles/particle_array['//int2str(i-1)//']/attributes/attribute')
       end if

       ! Enable particles to drift with the mesh
       if (have_option("/particles/move_with_mesh")) then
          particle_lists(i)%move_with_mesh=.true.
       end if
       
       ! Set flag for NaN particle output
       if (have_option("/particles/write_nan_outside_domain")) then
          particle_lists(i)%write_nan_outside=.true.
       end if

       
       !Read particles from options
       if (.not.from_checkpoint) then
          ewrite(2,*) "Reading particles from options"

          call get_option(trim(buffer)//"/name", funcnam)
          
          str_size=len_trim(int2str(j))
          fmt="(a,I"//int2str(str_size)//"."//int2str(str_size)//")"

          if (.not.have_option(trim(buffer)//"/initial_position/from_file")) then
             ! Reading particles from a python function
             call get_option(trim(buffer)//"/initial_position/python", func)
             allocate(coords(dim,j))
             call set_detector_coords_from_python(coords, j, func, current_time)
             
             do l=1,j
                write(particle_name, fmt) trim(funcnam)//"_", l
                particle_lists(i)%detector_names(l)=trim(particle_name)
                call create_single_particle(particle_lists(i), xfield, coords(:,l), &
                     attribute_dims, l, LAGRANGIAN_DETECTOR, trim(particle_name))
             end do
             if (attribute_dims.ne.0) then
                call set_particle_attributes(state, attribute_dims, xfield, dim, current_time, i)
             end if
             
             deallocate(coords)
          else
             
             ! Reading from a binary file where the user has placed the particle information
             particle_file_unit=free_unit()
             call get_option("/particles/particle_array/initial_position/from_file/file_name",particle_file_filename)
             
#ifdef STREAM_IO
             open(unit = particle_file_unit, file = trim(particle_file_filename), &
                  & action = "read", access = "stream", form = "unformatted")
#else
             FLAbort("No stream I/O support")
#endif

             do l=1,j
                write(particle_name, fmt) trim(funcnam)//"_", l
                particle_lists(i)%detector_names(l)=trim(particle_name)
                read(particle_file_unit) particle_location
                call create_single_particle(particle_lists(i), xfield, particle_location, &
                     attribute_dims, l, LAGRANGIAN_DETECTOR, trim(particle_name), attribute_vals)          
             end do
             deallocate(attribute_vals)  
          end if
       else
          ewrite(2,*) "Reading particles from checkpoint"
          ! If reading from checkpoint file:
          ! Particles checkpoint file names end in _par, with.groups appended for the header file
          ! and .attributes.dat appended for the binary data file that holds the positions and attributes

          allocate(attribute_vals(attribute_dims))
          
          particle_checkpoint_unit=free_unit()
          call get_option("/particles/particle_array/from_checkpoint_file/file_name",particles_cp_filename)
          
#ifdef STREAM_IO
          open(unit = particle_checkpoint_unit, file = trim(particles_cp_filename) // '.attributes.dat', &
               & action = "read", access = "stream", form = "unformatted")
#else
          FLAbort("No stream I/O support")
#endif
          
          !Read in particle locations from checkpoint file    
          call get_option(trim(buffer)//"/name", funcnam)     
          allocate(packed_buff(dim+attribute_dims))
          str_size=len_trim(int2str(j))
          fmt="(a,I"//int2str(str_size)//"."//int2str(str_size)//")"

          do m=1,j
             write(particle_name, fmt) trim(funcnam)//"_", m
             read(particle_checkpoint_unit) packed_buff
             particle_location=packed_buff(1:dim)
             if (attribute_dims.NE.0) then
                attribute_vals=packed_buff(dim+1:dim+attribute_dims)
             end if
             call create_single_particle(particle_lists(i), xfield, &
                  particle_location, attribute_dims, m, LAGRANGIAN_DETECTOR, trim(particle_name), attribute_vals)
             particle_lists(i)%detector_names(m)=trim(particle_name)
          end do
          deallocate(packed_buff)
          deallocate(attribute_vals)
       end if !from checkpoint
       
       !Set type of output file
       
       particle_lists(i)%binary_output = .true.
       if (have_option("/particles/ascii_output")) then
          particle_lists(i)%binary_output= .false.
          if(isparallel()) then
             FLAbort("No support for ascii detector output in parallel. Please use binary output.")
          end if
       end if

       ! Only the first process should write the header file
       if (getprocno() == 1) then
          particle_lists(i)%output_unit=free_unit()
          open(unit=particle_lists(i)%output_unit, file=trim(filename)//'.particles.'//trim(funcnam), action="write")
          
          write(particle_lists(i)%output_unit, '(a)') "<header>"
          call initialise_constant_diagnostics(particle_lists(i)%output_unit, &
               binary_format = particle_lists(i)%binary_output)
          
          ! Initial columns are elapsed time and dt.
          buffer=field_tag(name="ElapsedTime", column=1, statistic="value")
          write(particle_lists(i)%output_unit, '(a)') trim(buffer)
          buffer=field_tag(name="dt", column=2, statistic="value")
          write(particle_lists(i)%output_unit, '(a)') trim(buffer)
          
          ! Next columns contain the positions of all the particles
          column=2
          positionloop: do m=1,j
             buffer=field_tag(name=particle_lists(i)%detector_names(m), column=column+1,&
                  statistic="position", components=xfield%dim)
             write(particle_lists(i)%output_unit, '(a)') trim(buffer)
             column=column+xfield%dim 
          end do positionloop
          
          ! Next columns contain attributes of particles
          attributeloop: do m=1,j
             do n=1,attribute_dims
                call get_option ('/particles/particle_array['&
                     //int2str(i-1)//']/attributes/attribute['//int2str(n-1)//']/name', name)
                buffer=field_tag(name=particle_lists(i)%detector_names(m), column=column+1,&
                     statistic=name, components=1)
                write(particle_lists(i)%output_unit, '(a)')trim(buffer)
                column=column+1
             end do
          end do attributeloop            
       end if

       if (getprocno() == 1) then
          write(particle_lists(i)%output_unit, '(a)') "</header>"
          flush(particle_lists(i)%output_unit)
          ! when using mpi_subroutines to write into the particles file we need to close the file since 
          ! filename.particles.dat needs to be open now with MPI_OPEN
          if ((isparallel()).or.(particle_lists(i)%binary_output)) then
             close(particle_lists(i)%output_unit)
          end if
       end if


       if ((isparallel()).or.(particle_lists(i)%binary_output)) then
          ! bit of hack to delete any existing .particles.dat file
          ! if we don't delete the existing .particles.dat would simply be opened for random access and 
          ! gradually overwritten, mixing particle output from the current with that of a previous run
          call MPI_FILE_OPEN(MPI_COMM_FEMTOOLS, trim(filename) // '.particles.'//trim(funcnam)//'.dat', MPI_MODE_CREATE + MPI_MODE_RDWR + MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, particle_lists(i)%mpi_fh, IERROR)
          ewrite(2,*) "1", IERROR
          call MPI_FILE_CLOSE(particle_lists(i)%mpi_fh, IERROR)
          ewrite(2,*) "2", IERROR
          call MPI_FILE_OPEN(MPI_COMM_FEMTOOLS, trim(filename) // '.particles.'//trim(funcnam)//'.dat', MPI_MODE_CREATE + MPI_MODE_RDWR, MPI_INFO_NULL, particle_lists(i)%mpi_fh, IERROR)
          ewrite(2,*) "3", IERROR, MPI_SUCCESS
          assert(ierror == MPI_SUCCESS)
       end if
       !Get options for lagrangian particle movement
       call read_detector_move_options(particle_lists(i), "/particles")    
    end do

    do i = 1,python_particles_func
       ! And finally some sanity checks
       write(buffer, "(a,i0,a)") "/particles/particle_array[",i-1,"]"
       call get_option(trim(buffer)//"/name",name)
       totaldet_global=particle_lists(i)%length
       call allsum(totaldet_global)
       ewrite(2,*) "Found", particle_lists(i)%length, "local and ", totaldet_global, "global particles for particle array ", trim(name)

       assert(totaldet_global==particle_lists(i)%total_num_det)
    end do


  end subroutine initialise_particles

  subroutine create_single_particle(detector_list, xfield, position, attribute_dims, id, type, name, attribute_vals)
    ! Allocate a single particle, populate and insert it into the given list
    ! In parallel, first check if the particle would be local and only allocate if it is
    type(detector_linked_list), intent(inout) :: detector_list
    type(vector_field), pointer :: xfield
    real, dimension(xfield%dim), intent(in) :: position
    integer, intent(in) :: id, type
    character(len=*), intent(in) :: name

    type(detector_type), pointer :: detector
    type(element_type), pointer :: shape
    real, dimension(xfield%dim+1) :: lcoords
    integer :: element
    real, dimension(attribute_dims), intent(in), optional :: attribute_vals
    integer, intent(in) :: attribute_dims

    real ::  dt
    
    shape=>ele_shape(xfield,1)
    assert(xfield%dim+1==local_coord_count(shape))
    ! Determine element and local_coords from position
    ! In parallel, global=.false. can often work because there will be
    ! a halo of non-owned elements in your process and so you can work out
    ! ownership without communication.  But in general it won't work.
    call picker_inquire(xfield,position,element,local_coord=lcoords,global=.true.)
    call get_option("/timestepping/timestep", dt)
    ! If we're in parallel and don't own the element, skip this particle
    if (isparallel()) then
       if (element<0) return
       if (.not.element_owned(xfield,element)) return
    else
       ! In serial make sure the particle is in the domain
       ! unless we have the write_nan_outside override
       if (element<0 .and. .not.detector_list%write_nan_outside) then
          ewrite(-1,*) "Dealing with particle ", id, " named: ", trim(name)
          FLExit("Trying to initialise particle outside of computational domain")
       end if
    end if
    ! Otherwise, allocate and insert particle
    allocate(detector)
    allocate(detector%position(xfield%dim))
    allocate(detector%local_coords(local_coord_count(shape)))
!    allocate(detector%attributes(attribute_dims))
    call insert(detector, detector_list)
    ! Populate particle
    detector%name=name
    detector%position=position
    detector%element=element
    detector%local_coords=lcoords
    detector%type=type
    detector%id_number=id

    if (present(attribute_vals)) then
       allocate(detector%attributes(attribute_dims))
       detector%attributes = attribute_vals
    end if
  end subroutine create_single_particle

  subroutine move_particles(state, dt, timestep)
    !!Routine to loop over particle arrays and call move_lagrangian_detectors
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: dt
    integer, intent(in) :: timestep

    integer :: particle_arrays, attribute_dims
    integer :: i

    !Check whether there are any particles.
    particle_arrays = option_count("/particles/particle_array")
    if (particle_arrays==0) return

    ewrite(2,*), "In move_particles"
    do i = 1, particle_arrays
       attribute_dims = option_count("/particles/particle_array["//int2str(i-1)// &
            "]/attributes/attribute")
       call move_lagrangian_detectors(state, particle_lists(i), dt, timestep, attribute_dims)
    end do

  end subroutine move_particles

  subroutine set_particle_attributes(state, attribute_dims, xfield, dim, time, i)
    !!Routine to set particle attributes
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: time
    type(vector_field), pointer, intent(in) :: xfield
    character(len=PYTHON_FUNC_LEN) :: func
    character(len=FIELD_NAME_LEN), allocatable, dimension(:) :: field_name
    integer, intent(in) :: i, dim, attribute_dims

    type(detector_type), pointer :: particle

    real, allocatable, dimension(:,:) :: position
    real, allocatable, dimension(:,:) :: attribute_array

    real :: constant
    integer :: k, j, n, m, nparticles
    real, allocatable, dimension(:,:) :: lcoords
    integer, allocatable, dimension(:) :: ele

    nparticles = particle_lists(i)%length

    if (nparticles.eq.0) then
       return
    end if
    allocate(position(dim,nparticles))
    allocate(attribute_array(attribute_dims,nparticles))
    particle => particle_lists(i)%first
    allocate(lcoords(size(particle%local_coords),nparticles))
    allocate(ele(nparticles))
    do j = 1,nparticles
       position(:,j) = particle%position
       lcoords(:,j) = particle%local_coords
       ele(j) = particle%element
       particle => particle%next
    end do
    
    do k = 0, attribute_dims-1
       if (have_option('/particles/particle_array['// &
            int2str(i-1)//']/attributes/attribute['//int2str(k)//']/constant')) then
          call get_option('/particles/particle_array['// &
               int2str(i-1)//']/attributes/attribute['//int2str(k)//']/constant', constant)
          call set_particle_constant_from_options(attribute_array(k+1,:), nparticles, constant)
       else if (have_option('/particles/particle_array['// &
            int2str(i-1)//']/attributes/attribute['//int2str(k)//']/python')) then
          call get_option('/particles/particle_array['// &
               int2str(i-1)//']/attributes/attribute['//int2str(k)//']/python', func)
          call set_particle_attribute_from_python(attribute_array(k+1,:), position(:,:), nparticles, func, time)
       else if (have_option('/particles/particle_array['// &
            int2str(i-1)//']/attributes/attribute['//int2str(k)//']/python_fields')) then
          call get_option('/particles/particle_array['// &
               int2str(i-1)//']/attributes/attribute['//int2str(k)//']/python_fields', func)
          m=option_count('/particles/particle_array['// &
               int2str(i-1)//']/attributes/attribute['//int2str(k)//']/python_fields/field_name')
          allocate(field_name(m))
          do n=0,m-1
             call get_option('/particles/particle_array['// &
                  int2str(i-1)//']/attributes/attribute['//int2str(k)// &
                  ']/python_fields/field_name['//int2str(n)//']/name', field_name(n+1))
          end do
          call set_particle_fields_from_python(state, xfield, dim, position(:,:), nparticles, ele(:), lcoords(:,:), attribute_array(k+1,:), func, time, field_name)
          deallocate(field_name)
       end if
    end do

    particle => particle_lists(i)%first
    do j = 1,nparticles
       if (.not. allocated(particle%attributes)) then
          allocate(particle%attributes(attribute_dims))
       end if
       particle%attributes = attribute_array(:,j)
       particle => particle%next
    end do

    deallocate(position)
    deallocate(attribute_array)
    deallocate(lcoords)
    deallocate(ele)

  end subroutine set_particle_attributes

  subroutine update_particle_attributes(state, time)
    !!Routine to loop over particle arrays and update particle attributes
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: time
    type(vector_field), pointer :: xfield

    integer :: i
    integer :: particle_arrays, dim, attribute_dims

    !Check whether there are any particles.
    particle_arrays = option_count("/particles/particle_array")

    if (particle_arrays==0) return

    ewrite(2,*), "In update_particle_attributes"

    !Allocate parameters
    xfield=>extract_vector_field(state(1), "Coordinate")
    call get_option("/geometry/dimension",dim)

    !Update particle attributes by array
    do i = 1,particle_arrays
       attribute_dims = option_count('/particles/particle_array['//int2str(i-1)//']/attributes/attribute')
       call set_particle_attributes(state, attribute_dims, xfield, dim, time, i)
    end do

  end subroutine update_particle_attributes

  subroutine write_particles_loop(state, time, dt)
    !!Subroutine to loop over particle_lists and call write_particles for each list
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: time, dt

    integer :: particle_arrays, attribute_dims
    integer :: i

    !Check whether there are any particles.
    particle_arrays = option_count("/particles/particle_array")
    if (particle_arrays==0) return

    do i = 1, particle_arrays
       attribute_dims = option_count("/particles/particle_array["//int2str(i-1)// &
            "]/attributes/attribute")
       call write_particles(state, particle_lists(i), attribute_dims, time, dt)
    end do
    

  end subroutine write_particles_loop

  subroutine write_particles(state, detector_list, attribute_dims, time, dt)
    !!< Write values of particles to the previously opened particles file.
    type(state_type), dimension(:), intent(in) :: state
    type(detector_linked_list), intent(inout) :: detector_list
    real, intent(in) :: time, dt

    character(len=10) :: format_buffer
    integer :: i, check_no_det, totaldet_global
    integer, intent(in) :: attribute_dims
    type(detector_type), pointer :: detector

    ewrite(1,*) "In write_particles"

    !Computing the global number of particles. This is to prevent hanging
    !when there are no particles on any processor
    check_no_det=1
    if (detector_list%length==0) then
       check_no_det=0
    end if
    call allmax(check_no_det)
    if (check_no_det==0) then
       return
    end if

    ! If isparallel() or binary output use this:
    if ((isparallel()).or.(detector_list%binary_output)) then    
       call write_mpi_out_particles(state,detector_list, attribute_dims, time,dt)
    else ! This is only for single processor with ascii output
       if(getprocno() == 1) then
          if(detector_list%binary_output) then
             write(detector_list%output_unit) time
             write(detector_list%output_unit) dt
          else
             format_buffer=reals_format(1)
             write(detector_list%output_unit, format_buffer, advance="no") time
             write(detector_list%output_unit, format_buffer, advance="no") dt
          end if
       end if

       ! Next columns contain the positions of all the particles.
       detector => detector_list%first
       positionloop: do i=1, detector_list%length
          if(detector_list%binary_output) then
             write(detector_list%output_unit) detector%position
          else
             format_buffer=reals_format(size(detector%position))
             write(detector_list%output_unit, format_buffer, advance="no") &
                  detector%position
          end if
          detector => detector%next
       end do positionloop

       ! Next columns contain the attributes of particles
       detector => detector_list%first
       attributeloop: do i=1,detector_list%length
          if (attribute_dims.ne.0) then
             if(detector_list%binary_output) then
                write(detector_list%output_unit) detector%attributes
             else
                format_buffer=reals_format(attribute_dims)
                write(detector_list%output_unit, format_buffer, advance="no") &
                     detector%attributes
             end if
          end if
          detector => detector%next
       end do attributeloop
       
       ! Output end of line
       if (.not. detector_list%binary_output) then
          ! Output end of line
          write(detector_list%output_unit,'(a)') ""
       end if
       flush(detector_list%output_unit)
    end if

    totaldet_global=detector_list%length
    call allsum(totaldet_global)
    ewrite(2,*) "Found", detector_list%length, "local and", totaldet_global, "global detectors"
    
    if (totaldet_global/=detector_list%total_num_det) then
       ewrite(2,*) "We have either duplication or have lost some det"
       ewrite(2,*) "totaldet_global", totaldet_global
       ewrite(2,*) "total_num_det", detector_list%total_num_det
    end if
    
    ewrite(1,*) "Exiting write_particles"
    
  contains
    
    function reals_format(reals)
      character(len=10) :: reals_format
      integer :: reals
      
      write(reals_format, '(a,i0,a)') '(',reals,'e15.6e3)'
      
    end function reals_format
    
  end subroutine write_particles

  subroutine write_mpi_out_particles(state,detector_list, attribute_dims, time,dt)
    !!< Writes particle information (position, attributes, etc.) into particle file using MPI output 
    ! commands so that when running in parallel all processors can write at the same time information into the file at the right location.       

    type(state_type), dimension(:), intent(in) :: state
    type(detector_linked_list), intent(inout) :: detector_list
    real, intent(in) :: time, dt
    integer, intent(in) :: attribute_dims

    integer :: i, ierror, realsize, dim, procno
    integer(KIND = MPI_OFFSET_KIND) :: location_to_write, offset
    integer :: number_total_columns

    type(vector_field), pointer :: vfield
    type(detector_type), pointer :: node

    ewrite(2, *) "In write_mpi_out_particles"

    detector_list%mpi_write_count = detector_list%mpi_write_count + 1
    ewrite(2, *) "Writing particle output ", detector_list%mpi_write_count
    
    procno = getprocno()


    call mpi_type_extent(getpreal(), realsize, ierror)
    assert(ierror == MPI_SUCCESS)

    vfield => extract_vector_field(state, "Coordinate")
    dim = vfield%dim
                           ! Time data
    number_total_columns = 2 + &
                           ! Particle coordinates
                         & detector_list%total_num_det * dim + &
                           ! Particle attributes
                         & detector_list%total_num_det * attribute_dims

    ! raise kind of one of the variables (each individually is a 4 byte-integer) such that the calculation is coerced to be of MPI_OFFSET_KIND (typically 8 bytes)
    ! this is necessary for files bigger than 2GB
    location_to_write = (int(detector_list%mpi_write_count, kind=MPI_OFFSET_KIND) - 1) * number_total_columns * realsize

    if(procno == 1) then
      ! Output time data
      call mpi_file_write_at(detector_list%mpi_fh, location_to_write, time, 1, getpreal(), MPI_STATUS_IGNORE, ierror)
      assert(ierror == MPI_SUCCESS)
        
      call mpi_file_write_at(detector_list%mpi_fh, location_to_write + realsize, dt, 1, getpreal(), MPI_STATUS_IGNORE, ierror)
      assert(ierror == MPI_SUCCESS)
    end if
    location_to_write = location_to_write + 2 * realsize

    node => detector_list%first
    position_loop: do i = 1, detector_list%length
      ! Output detector coordinates
      assert(size(node%position) == dim)  
    
      offset = location_to_write + (node%id_number - 1) * dim * realsize

      call mpi_file_write_at(detector_list%mpi_fh, offset, node%position, dim, getpreal(), MPI_STATUS_IGNORE, ierror)
      assert(ierror == MPI_SUCCESS)
      node => node%next
    end do position_loop
    assert(.not. associated(node))
    location_to_write = location_to_write + detector_list%total_num_det * dim * realsize

    if (attribute_dims.ne.0) then
       node => detector_list%first
       attribute_loop: do i = 1, detector_list%length
          assert(size(node%attributes) == attribute_dims)
          
          offset = location_to_write + (node%id_number - 1) * attribute_dims * realsize
          call mpi_file_write_at(detector_list%mpi_fh, offset, node%attributes, attribute_dims, getpreal(), MPI_STATUS_IGNORE, ierror)
          assert(ierror == MPI_SUCCESS)
          node => node%next
       end do attribute_loop
       assert(.not. associated(node))
    end if
    
    call mpi_file_sync(detector_list%mpi_fh, ierror)
    assert(ierror == MPI_SUCCESS)

!    ! The following was used when debugging to check some of the data written
!    ! into the file
!    ! Left here in case someone would like to use the mpi_file_read_at for
!    ! debugging or checking
!   
!    number_total_columns = 2 + total_num_det * dim
!    allocate(buffer(number_total_columns))
!
!    call mpi_file_read_at(fh, 0, buffer, size(buffer), getpreal(), status, ierror)
!    call mpi_get_count(status, getpreal(), count,  ierror)
!    assert(ierror == MPI_SUCCESS)
!
!    call mpi_barrier(MPI_COMM_FEMTOOLS, ierror)
!    assert(ierror == MPI_SUCCESS)
!    
!    deallocate(buffer)
!    ewrite(2, "(a,i0,a)") "Read ", count, " reals"

    ewrite(2, *) "Exiting write_mpi_out"
   
  end subroutine write_mpi_out_particles

  subroutine destroy_particles()
    !Deallocate all particle arrays (detector lists)

    deallocate(particle_lists)
    
  end subroutine destroy_particles

end module particles

  
    
    
    
