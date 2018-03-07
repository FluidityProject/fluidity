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
  use iso_c_binding, only: C_NULL_CHAR

  implicit none

  private

  public :: initialise_particles, move_particles, write_particles_loop, destroy_particles, &
            update_particle_attributes, checkpoint_particles_loop

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
    real, allocatable, dimension(:,:) :: attribute_vals
    real, allocatable, dimension(:) :: packed_buff
    real, allocatable, dimension(:) :: particle_location
    type(vector_field), pointer :: xfield
    real:: current_time

    integer, dimension(3) :: attributes_buffer

    integer :: str_size
    integer :: i, j, m, n, l
    integer :: python_particles_func, dim
    integer :: particle_file_unit=0, particle_checkpoint_unit=0
    integer :: column, IERROR, totaldet_global
    integer :: nfields, nprescribed, ndiagnostic, nprognostic


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

    do i = 1,python_particles_func

       ! If the option
       ! "from_checkpoint_file" exists, it means we are continuing the simulation
       ! after checkpointing and the reading of the particle positions must be
       ! done from a file
       
       if (have_option("/particles/particle_array["//int2str(i-1)//"]/initial_position/from_checkpoint_file")) then
          from_checkpoint=.true.
       else
          from_checkpoint=.false.
       end if
       write(buffer, "(a,i0,a)") "/particles/particle_array[",i-1,"]"
       call get_option(trim(buffer)//"/number_of_particles", j)
       call get_option(trim(buffer)//"/name", funcnam)
       particle_lists(i)%total_num_det=j
       allocate(particle_lists(i)%detector_names(j))
       !Register this I/O list with a global list of detectors/particles
       call register_detector_list(particle_lists(i))

       !Find number of attributes 
       attributes_buffer(1)=0
       if (have_option('/particles/particle_array['//int2str(i-1)//']/attributes/attribute')) then
          attributes_buffer(1)=option_count('/particles/particle_array['//int2str(i-1)//']/attributes/attribute')
       end if
       attributes_buffer(2)=0
       attributes_buffer(3)=0
       do m = 1,attributes_buffer(1)
          if (have_option('/particles/particle_array['//int2str(i-1)//']/attributes/attribute['//int2str(m-1)//']/python_fields')) then
             nprognostic = option_count('/material_phase/scalar_field/prognostic/particles/include_in_particles/store_old_field')
             nprescribed = option_count('/material_phase/scalar_field/prescribed/particles/include_in_particles/store_old_field')
             ndiagnostic = option_count('/material_phase/scalar_field/diagnostic/particles/include_in_particles/store_old_field')
             attributes_buffer(3) = ndiagnostic + nprescribed + nprognostic
             if (have_option('/particles/particle_array['//int2str(i-1)//']/attributes/attribute['//int2str(m-1)//']/python_fields/store_old_attribute')) then
                attributes_buffer(2)=attributes_buffer(2)+1
             end if
          end if
       end do

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
          
          str_size=len_trim(int2str(j))
          fmt="(a,I"//int2str(str_size)//"."//int2str(str_size)//")"

          if (.not.have_option(trim(buffer)//"/initial_position/from_file")) then
             ! Reading particles from a python function
             call get_option(trim(buffer)//"/initial_position/python", func)
             allocate(coords(dim,j))
             call set_detector_coords_from_python(coords, j, func, current_time)
             
             do l=1,j
                write(particle_name, fmt) trim(funcnam)//"_", l
                call create_single_particle(particle_lists(i), xfield, coords(:,l), &
                     l, LAGRANGIAN_DETECTOR, trim(particle_name))
             end do
             if (attributes_buffer(1).ne.0) then
                call set_particle_attributes(state, attributes_buffer, xfield, dim, current_time, i)
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
             
             !do l=1,j
             !   write(particle_name, fmt) trim(funcnam)//"_", l
             !   read(particle_file_unit) particle_location
             !   call create_single_particle(particle_lists(i), xfield, particle_location, &
             !        attribute_dims, l, LAGRANGIAN_DETECTOR, trim(particle_name), attribute_vals)          
             !end do
             !deallocate(attribute_vals)  !!!this hasn't been allocated
          end if
       else
          ewrite(2,*) "Reading particles from checkpoint"
          ! If reading from checkpoint file:
          ! Particles checkpoint file names end in _par, with.groups appended for the header file
          ! and .attributes.dat appended for the binary data file that holds the positions and attributes

          allocate(attribute_vals(3,maxval(attributes_buffer)))
          
          particle_checkpoint_unit=free_unit()
          call get_option("/particles/particle_array["//int2str(i-1)//"]/initial_position/from_checkpoint_file/file_name",particles_cp_filename)
          
#ifdef STREAM_IO
          open(unit = particle_checkpoint_unit, file = trim(particles_cp_filename) //'.attributes.dat', &
               & action = "read", access = "stream", form = "unformatted")
#else
          FLAbort("No stream I/O support")
#endif
          
          !Read in particle locations from checkpoint file       
          allocate(packed_buff(dim+attributes_buffer(1)+attributes_buffer(2)+attributes_buffer(3)))
          str_size=len_trim(int2str(j))
          fmt="(a,I"//int2str(str_size)//"."//int2str(str_size)//")"

          do m=1,j
             write(particle_name, fmt) trim(funcnam)//"_", m
             read(particle_checkpoint_unit) packed_buff
             particle_location=packed_buff(1:dim)
             if (attributes_buffer(1).NE.0) then
                attribute_vals(1,1:attributes_buffer(1))=packed_buff(dim+1:dim+attributes_buffer(1))
                attribute_vals(2,1:attributes_buffer(2))=packed_buff(dim+1+attributes_buffer(1):dim+attributes_buffer(1)+attributes_buffer(2))
                attribute_vals(3,1:attributes_buffer(3))=packed_buff(dim+1+attributes_buffer(1)+attributes_buffer(2):dim+attributes_buffer(1)+ &
                     attributes_buffer(2)+attributes_buffer(3))
             end if
             call create_single_particle(particle_lists(i), xfield, &
                  particle_location, m, LAGRANGIAN_DETECTOR, trim(particle_name),attributes_buffer=attributes_buffer,attribute_vals= attribute_vals)
          end do
          deallocate(packed_buff)
          deallocate(attribute_vals)
       end if !from checkpoint
       
       !Set type of output file
       
       particle_lists(i)%binary_output = .true.
       if (have_option("/particles/ascii_output")) then
          particle_lists(i)%binary_output= .false.
          if(isparallel()) then
             FLExit("Error: No support for ascii detector output in parallel. Please use binary output by turning off the ascii_output option.")
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
             do n=1,attributes_buffer(1)
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
          if (particle_lists(i)%binary_output) then
             close(particle_lists(i)%output_unit)
          end if
       end if


       if (particle_lists(i)%binary_output) then
          ! bit of hack to delete any existing .particles.dat file
          ! if we don't delete the existing .particles.dat would simply be opened for random access and 
          ! gradually overwritten, mixing particle output from the current with that of a previous run
          call MPI_FILE_OPEN(MPI_COMM_FEMTOOLS, trim(filename) // '.particles.'//trim(funcnam)//'.dat', MPI_MODE_CREATE + MPI_MODE_RDWR + MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, particle_lists(i)%mpi_fh, IERROR)
          call MPI_FILE_CLOSE(particle_lists(i)%mpi_fh, IERROR)
          call MPI_FILE_OPEN(MPI_COMM_FEMTOOLS, trim(filename) // '.particles.'//trim(funcnam)//'.dat', MPI_MODE_CREATE + MPI_MODE_RDWR, MPI_INFO_NULL, particle_lists(i)%mpi_fh, IERROR)
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

  subroutine create_single_particle(detector_list, xfield, position, id, type, name, attributes_buffer, attribute_vals)
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
    real, dimension(:,:), intent(in), optional :: attribute_vals
    integer, dimension(3), intent(in), optional :: attributes_buffer

    real ::  dt
    
    shape=>ele_shape(xfield,1)
    assert(xfield%dim+1==local_coord_count(shape))
    detector_list%detector_names(id)=name
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
    call insert(detector, detector_list)
    ! Populate particle
    detector%name=name
    detector%position=position
    detector%element=element
    detector%local_coords=lcoords
    detector%type=type
    detector%id_number=id

    ewrite(2,*) "147258", id

    if (present(attribute_vals)) then
       allocate(detector%attributes(attributes_buffer(1)))
       detector%attributes = attribute_vals(1,1:attributes_buffer(1))
       allocate(detector%old_attributes(attributes_buffer(2)))
       detector%old_attributes = attribute_vals(2,1:attributes_buffer(2))
       allocate(detector%old_fields(attributes_buffer(3)))
       detector%old_fields = attribute_vals(3,1:attributes_buffer(3))
    end if
  end subroutine create_single_particle

  subroutine move_particles(state, dt, timestep)
    !!Routine to loop over particle arrays and call move_lagrangian_detectors
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: dt
    integer, intent(in) :: timestep

    integer, dimension(3) :: attributes_buffer
    integer :: particle_arrays
    integer :: i, m
    integer :: nprescribed, ndiagnostic, nprognostic

    !Check whether there are any particles.
    particle_arrays = option_count("/particles/particle_array")
    if (particle_arrays==0) return

    ewrite(2,*), "In move_particles"
    do i = 1, particle_arrays
       attributes_buffer(1)=0
       if (have_option('/particles/particle_array['//int2str(i-1)//']/attributes/attribute')) then
          attributes_buffer(1)=option_count('/particles/particle_array['//int2str(i-1)//']/attributes/attribute')
       end if
       attributes_buffer(2)=0
       attributes_buffer(3)=0
       do m = 1,attributes_buffer(1)
          if (have_option('/particles/particle_array['//int2str(i-1)//']/attributes/attribute['//int2str(m-1)//']/python_fields')) then
             nprognostic = option_count('/material_phase/scalar_field/prognostic/particles/include_in_particles/store_old_field')
             nprescribed = option_count('/material_phase/scalar_field/prescribed/particles/include_in_particles/store_old_field')
             ndiagnostic = option_count('/material_phase/scalar_field/diagnostic/particles/include_in_particles/store_old_field')
             attributes_buffer(3) = ndiagnostic + nprescribed + nprognostic
             if (have_option('/particles/particle_array['//int2str(i-1)//']/attributes/attribute['//int2str(m-1)//']/python_fields/store_old_attribute')) then
                attributes_buffer(2)=attributes_buffer(2)+1
             end if
          end if
       end do
       call move_lagrangian_detectors(state, particle_lists(i), dt, timestep, attributes_buffer)
    end do

  end subroutine move_particles

  subroutine set_particle_attributes(state, attributes_buffer, xfield, dim, time, i)
    !!Routine to set particle attributes
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: time
    type(vector_field), pointer, intent(in) :: xfield
    character(len=PYTHON_FUNC_LEN) :: func
    character(len=FIELD_NAME_LEN), allocatable, dimension(:) :: field_name
    integer, intent(in) :: i, dim
    integer, dimension(3), intent(in) :: attributes_buffer

    type(detector_type), pointer :: particle

    real, allocatable, dimension(:,:) :: positions
    real, allocatable, dimension(:,:) :: attribute_array
    real, allocatable, dimension(:,:) :: old_attributes
    character, allocatable, dimension(:,:) :: old_att_names
    real, allocatable, dimension(:) :: zeros_dum
    character(len = OPTION_PATH_LEN) :: old_name

    real :: constant
    integer :: k, j, nparticles, l, m
    real, allocatable, dimension(:,:) :: lcoords
    integer, allocatable, dimension(:) :: ele

    nparticles = particle_lists(i)%length

    if (nparticles.eq.0) then
       return
    end if
    
    !Set parameters to calculate attributes
    
    allocate(positions(dim,nparticles))
    allocate(attribute_array(attributes_buffer(1),nparticles))
    particle => particle_lists(i)%first
    allocate(lcoords(size(particle%local_coords),nparticles))
    allocate(ele(nparticles))
    allocate(old_attributes(attributes_buffer(2),nparticles))
    allocate(zeros_dum(attributes_buffer(2)))
    do k = 1,attributes_buffer(2)
       zeros_dum(k) = 0
    end do
    do j = 1,nparticles
       positions(:,j) = particle%position
       lcoords(:,j) = particle%local_coords
       ele(j) = particle%element
       if (allocated(particle%old_attributes)) then
          old_attributes(:,j) = particle%old_attributes
       else
          old_attributes(:,j) = zeros_dum(:)
       end if
       particle => particle%next
    end do

    allocate(old_att_names(FIELD_NAME_LEN,attributes_buffer(2)))
    l=1
    do k = 1,attributes_buffer(1)
       if (have_option('/particles/particle_array['// &
            int2str(i-1)//']/attributes/attribute['//int2str(k-1)//']/python_fields/store_old_attribute')) then
          call get_option('/particles/particle_array['// &
               int2str(i-1)//']/attributes/attribute['//int2str(k-1)//']/name', old_name)
          old_att_names(1,l) = 'O'
          old_att_names(2,l) = 'l'
          old_att_names(3,l) = 'd'
          do j = 4,len_trim(old_name)+3
             old_att_names(j,l)=old_name(j-3:j-3)
          end do
          old_att_names(j,l) = C_NULL_CHAR
          l=l+1
       end if
    end do
    
    do k = 0, attributes_buffer(1)-1
       if (have_option('/particles/particle_array['// &
            int2str(i-1)//']/attributes/attribute['//int2str(k)//']/constant')) then
          particle => particle_lists(i)%first
          if (allocated(particle%attributes)) then
             do j = 1,nparticles   
                attribute_array(k+1,j) = particle%attributes(k+1)
                particle => particle%next
             end do
          else
             call get_option('/particles/particle_array['// &
                  int2str(i-1)//']/attributes/attribute['//int2str(k)//']/constant', constant)
             call set_particle_constant_from_options(attribute_array(k+1,:), nparticles, constant)
          end if
       else if (have_option('/particles/particle_array['// &
            int2str(i-1)//']/attributes/attribute['//int2str(k)//']/python')) then
          call get_option('/particles/particle_array['// &
               int2str(i-1)//']/attributes/attribute['//int2str(k)//']/python', func)
          call set_particle_attribute_from_python(attribute_array(k+1,:), positions(:,:), nparticles, func, time)
       else if (have_option('/particles/particle_array['// &
            int2str(i-1)//']/attributes/attribute['//int2str(k)//']/python_fields')) then
          call get_option('/particles/particle_array['// &
               int2str(i-1)//']/attributes/attribute['//int2str(k)//']/python_fields', func)
          call set_particle_fields_from_python(particle_lists(i), attributes_buffer, state, xfield, dim, positions(:,:), lcoords(:,:), ele(:), nparticles, attribute_array(k+1,:), old_att_names, old_attributes, func, time)
       else if (have_option('/particles/particle_array['// &
            int2str(i-1)//']/attributes/attribute['//int2str(k)//']/from_checkpoint_file')) then
          particle => particle_lists(i)%first
          do j = 1,nparticles   
             attribute_array(k+1,j) = particle%attributes(k+1)
             particle => particle%next
          end do
       end if
    end do
    
    !Set attribute values and old_attribute values
    particle => particle_lists(i)%first
    do j = 1,nparticles 
       if (.not.allocated(particle%attributes)) then
          allocate(particle%attributes(attributes_buffer(1)))
       end if
       if (.not.allocated(particle%old_attributes)) then
          allocate(particle%old_attributes(attributes_buffer(2)))
       end if
       particle%attributes = attribute_array(:,j)
       if (attributes_buffer(2).ne.0) then
          m = 1
          do k = 1,attributes_buffer(1)
             if (have_option('/particles/particle_array['// &
                  int2str(i-1)//']/attributes/attribute['//int2str(k-1)//']/python_fields/store_old_attribute')) then
                particle%old_attributes(m) = particle%attributes(k)
                m=m+1
             end if
          end do
       end if
       particle => particle%next
    end do

    if (attributes_buffer(3).ne.0) then
       call update_particle_fields(state, attributes_buffer, ele, lcoords, nparticles, i)
    end if
    
    deallocate(positions)
    deallocate(lcoords)
    deallocate(ele)
    deallocate(attribute_array)
    deallocate(old_attributes)
    deallocate(zeros_dum)
    deallocate(old_att_names)

  end subroutine set_particle_attributes

  subroutine update_particle_attributes(state, time)
    !!Routine to loop over particle arrays and update particle attributes
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: time
    type(vector_field), pointer :: xfield

    integer, dimension(3) :: attributes_buffer
    integer :: i, nprescribed, ndiagnostic, nprognostic, m
    integer :: particle_arrays, dim

    !Check whether there are any particles.
    particle_arrays = option_count("/particles/particle_array")

    if (particle_arrays==0) return

    ewrite(2,*), "In update_particle_attributes"

    !Allocate parameters
    xfield=>extract_vector_field(state(1), "Coordinate")
    call get_option("/geometry/dimension",dim)

    !Update particle attributes by array
    do i = 1,particle_arrays
       attributes_buffer(1)=0
       if (have_option('/particles/particle_array['//int2str(i-1)//']/attributes/attribute')) then
          attributes_buffer(1)=option_count('/particles/particle_array['//int2str(i-1)//']/attributes/attribute')
       end if
       attributes_buffer(2)=0
       attributes_buffer(3)=0
       do m = 1,attributes_buffer(1)
          if (have_option('/particles/particle_array['//int2str(i-1)//']/attributes/attribute['//int2str(m-1)//']/python_fields')) then
             nprognostic = option_count('/material_phase/scalar_field/prognostic/particles/include_in_particles/store_old_field')
             nprescribed = option_count('/material_phase/scalar_field/prescribed/particles/include_in_particles/store_old_field')
             ndiagnostic = option_count('/material_phase/scalar_field/diagnostic/particles/include_in_particles/store_old_field')
             attributes_buffer(3) = ndiagnostic + nprescribed + nprognostic
             if (have_option('/particles/particle_array['//int2str(i-1)//']/attributes/attribute['//int2str(m-1)//']/python_fields/store_old_attribute')) then
                attributes_buffer(2)=attributes_buffer(2)+1
             end if
          end if
       end do
       if (attributes_buffer(1).ne.0) then
          call set_particle_attributes(state, attributes_buffer, xfield, dim, time, i)
       end if
    end do

  end subroutine update_particle_attributes

  subroutine update_particle_fields(state, attributes_buffer, ele, lcoords, nparticles, i)
    
    type(state_type), dimension(:), intent(in) :: state
    integer, dimension(3), intent(in) :: attributes_buffer
    real, dimension(:,:), intent(in) :: lcoords
    integer, dimension(:), intent(in) :: ele
    integer, intent(in) :: nparticles
    integer, intent(in) :: i

    character(len = OPTION_PATH_LEN) :: name
    type(scalar_field), pointer :: sfield
    real, allocatable, dimension(:,:) :: old_field_vals
    real :: value
    type(detector_type), pointer :: particle
    integer :: phase, f, num_fields, l, j

    allocate(old_field_vals(attributes_buffer(3),nparticles))
    l=1
    do phase=1,size(state)
       num_fields = option_count('/material_phase[' &
            //int2str(phase-1)//']/scalar_field')
       do f = 1, num_fields
          call get_option('material_phase['//int2str(phase-1)//']/scalar_field['//int2str(f-1)//']/name', name)
          sfield => extract_scalar_field(state(phase),name)  
          if (have_option(trim(sfield%option_path)//"/prescribed/particles/include_in_particles/store_old_field").or. &
               have_option(trim(sfield%option_path)//"/diagnostic/particles/include_in_particles/store_old_field").or. &
               have_option(trim(sfield%option_path)//"/prognostic/particles/include_in_particles/store_old_field")) then
             do j = 1,nparticles 
                value = eval_field(ele(j), sfield, lcoords(:,j))
                old_field_vals(l,j)=value
             end do
             l=l+1
          end if
       end do
    end do

    particle => particle_lists(i)%first
    do j = 1,nparticles
       if (.not.allocated(particle%old_fields)) then
          allocate(particle%old_fields(attributes_buffer(3)))
       end if
       particle%old_fields=old_field_vals(:,j)
       particle=>particle%next
    end do

    deallocate(old_field_vals)

  end subroutine update_particle_fields

  subroutine write_particles_loop(state, time, dt)
    !!Subroutine to loop over particle_lists and call write_particles for each list
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: time, dt

    integer :: attribute_dims
    integer :: particle_arrays
    integer :: i

    !Check whether there are any particles.
    particle_arrays = option_count("/particles/particle_array")
    if (particle_arrays==0) return

    do i = 1, particle_arrays
       attribute_dims=option_count('/particles/particle_array['//int2str(i-1)//']/attributes/attribute')
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

    ! If binary output use this:
    if (detector_list%binary_output) then    
       call write_mpi_out_particles(state,detector_list, attribute_dims, time,dt)
    else ! This is only for single processor with ascii output
       if(getprocno() == 1) then
          format_buffer=reals_format(1)
          write(detector_list%output_unit, format_buffer, advance="no") time
          write(detector_list%output_unit, format_buffer, advance="no") dt
       end if

       ! Next columns contain the positions of all the particles.
       detector => detector_list%first
       positionloop: do i=1, detector_list%length
          format_buffer=reals_format(size(detector%position))
          write(detector_list%output_unit, format_buffer, advance="no") &
               detector%position
          detector => detector%next
       end do positionloop

       ! Next columns contain the attributes of particles
       detector => detector_list%first
       attributeloop: do i=1,detector_list%length
          if (attribute_dims.ne.0) then
             format_buffer=reals_format(attribute_dims)
             write(detector_list%output_unit, format_buffer, advance="no") &
                  detector%attributes
          end if
          detector => detector%next
       end do attributeloop
       
       ! Output end of line
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

  subroutine write_mpi_out_particles(state, detector_list, attribute_dims, time, dt)
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

  subroutine checkpoint_particles_loop(state,prefix,postfix,cp_no)
    !!Subroutine to loop over particle_lists and call checkpoint_particles for each list
    integer, parameter :: PREFIX_LEN = OPTION_PATH_LEN
    type(state_type), dimension(:), intent(in) :: state
    character(len = *), intent(in) :: prefix
    !! Default value "checkpoint"
    character(len = *), intent(in) :: postfix
    integer, optional, intent(in) :: cp_no
    character(len = PREFIX_LEN) :: lpostfix
    character(len = FIELD_NAME_LEN) :: name

    integer, dimension(3) :: attributes_buffer
    integer :: nprescribed, ndiagnostic, nprognostic
    integer :: particle_arrays
    integer :: i, m
    
    !Check whether there are any particles.
    particle_arrays = option_count("/particles/particle_array")
    if (particle_arrays==0) return

    ewrite(1, *) "Checkpointing particles"

    assert(len_trim(prefix) > 0)

    lpostfix = postfix

    do i = 1, particle_arrays
       attributes_buffer(1)=0
       if (have_option('/particles/particle_array['//int2str(i-1)//']/attributes/attribute')) then
          attributes_buffer(1)=option_count('/particles/particle_array['//int2str(i-1)//']/attributes/attribute')
       end if
       attributes_buffer(2)=0
       attributes_buffer(3)=0
       do m = 1,attributes_buffer(1)
          if (have_option('/particles/particle_array['//int2str(i-1)//']/attributes/attribute['//int2str(m-1)//']/python_fields')) then
             nprognostic = option_count('/material_phase/scalar_field/prognostic/particles/include_in_particles/store_old_field')
             nprescribed = option_count('/material_phase/scalar_field/prescribed/particles/include_in_particles/store_old_field')
             ndiagnostic = option_count('/material_phase/scalar_field/diagnostic/particles/include_in_particles/store_old_field')
             attributes_buffer(3) = ndiagnostic + nprescribed + nprognostic
             if (have_option('/particles/particle_array['//int2str(i-1)//']/attributes/attribute['//int2str(m-1)//']/python_fields/store_old_attribute')) then
                attributes_buffer(2)=attributes_buffer(2)+1
             end if
          end if
       end do
       call get_option("/particles/particle_array["//int2str(i-1)//"]/name", name)
       call checkpoint_particles(state,prefix,lpostfix,cp_no,particle_lists(i),attributes_buffer,name)
    end do

  end subroutine checkpoint_particles_loop

  subroutine checkpoint_particles(state,prefix,lpostfix,cp_no,particle_list,attributes_buffer,name)
    !!<Checkpoint Particles

    type(state_type), dimension(:), intent(in) :: state
    character(len = *), intent(in) :: prefix
    !! Default value "checkpoint"
    character(len = *), intent(in) :: lpostfix
    integer, optional, intent(in) :: cp_no
    character(len = *), intent(in) :: name

    type(detector_linked_list), intent(inout) :: particle_list
    integer, dimension(3), intent(in) :: attributes_buffer

    integer(KIND=MPI_OFFSET_KIND) :: location_to_write, offset
    
    type(detector_type), pointer :: node
    character(len = OPTION_PATH_LEN) :: particles_cp_filename
    type(vector_field), pointer :: vfield

    integer, ALLOCATABLE, DIMENSION(:) :: status
    real, dimension(:), allocatable :: buffer
    
    integer, save :: fhdet=0
    integer :: i, IERROR
    integer :: nints, realsize, dimen, num_particles, number_total_columns

    num_particles = particle_list%total_num_det

    ! Construct a new particle checkpoint filename
    !!!get name of particle array here to construct the output file
    particles_cp_filename = trim(prefix)
    if(present(cp_no)) particles_cp_filename = trim(particles_cp_filename) // "_" // int2str(cp_no)
    particles_cp_filename = trim(particles_cp_filename) // "_" // trim(lpostfix)

!!!!! Writing of position particles before checkpointing in serial  !!!!!

    !!< Writes particle last position into particles file using MPI output 
    ! commands so that when running in parallel all processors can write at the same time information into the file at the right location.
    
    call MPI_FILE_OPEN(MPI_COMM_FEMTOOLS, trim(particles_cp_filename) // '_particles.' // trim(name) // '.attributes.dat', MPI_MODE_CREATE + MPI_MODE_RDWR, MPI_INFO_NULL, fhdet, IERROR)

    ewrite(1,*) "after opening the IERROR is:", IERROR

    allocate( status(MPI_STATUS_SIZE) )

    call MPI_TYPE_EXTENT(getpreal(), realsize, ierror)

    vfield => extract_vector_field(state(1),"Velocity")

    dimen=vfield%dim

    number_total_columns=num_particles*(dimen+attributes_buffer(1)+attributes_buffer(2)+attributes_buffer(3)) !!!!or just *dimen?

    node => particle_list%first

    location_to_write=0

    positionloop_cp: do i=1, particle_list%length
      offset = location_to_write+(node%id_number-1)*(size(node%position)+attributes_buffer(1)+attributes_buffer(2)+attributes_buffer(3))*realsize
      ewrite(1,*) "after file set view position IERROR is:", IERROR

      allocate(buffer(size(node%position)+attributes_buffer(1)+attributes_buffer(2)+attributes_buffer(3)))
      buffer(1:size(node%position))=node%position
      if (attributes_buffer(1).NE.0) then
         buffer(1+size(node%position):size(node%position)+attributes_buffer(1))=node%attributes
      end if
      if (attributes_buffer(2).NE.0) then
         buffer(1+size(node%position)+attributes_buffer(1):size(node%position)+attributes_buffer(1) &
              +attributes_buffer(2))=node%old_attributes
      end if
      if (attributes_buffer(3).NE.0) then
         buffer(1+size(node%position)+attributes_buffer(1)+attributes_buffer(2):size(node%position)+attributes_buffer(1) &
              +attributes_buffer(2)+attributes_buffer(3))=node%old_fields
      end if
      nints=size(node%position)+attributes_buffer(1)+attributes_buffer(2)+attributes_buffer(3)
      
      call MPI_FILE_WRITE_AT(fhdet,offset,buffer,nints,getpreal(),status,IERROR)

      ewrite(1,*) "after sync position IERROR is:", IERROR
      deallocate(buffer)
      node => node%next
    end do positionloop_cp

    call update_particle_options(trim(particles_cp_filename) // "_particles", "binary", particle_list, name, attributes_buffer(1))

    if (fhdet/=0) then
       call MPI_FILE_CLOSE(fhdet, IERROR) 
       if (IERROR/=0) then
          ewrite(0,*) "Warning: failed to close .particles checkpoint file open with mpi_file_open"
       end if
    end if
    
  end subroutine checkpoint_particles

  subroutine update_particle_options(filename, format, particle_list, name, attribute_dims)
    !!< Updates the initial options of the detectors (options tree in diamond)

    character(len = *), intent(in) :: filename
    character(len = *), intent(in) :: format
    character(len = *), intent(in) :: name

    type(detector_linked_list), intent(inout) :: particle_list
    integer, intent(in) :: attribute_dims

    integer :: num_particles, i, stat
    logical :: particles_c

    character(len = 254) :: temp_string

    num_particles = particle_list%total_num_det

    temp_string=name

    ewrite(1,*) 'In update_particles_options'
    ewrite(1,*) temp_string
    
    call delete_option("/particles/particle_array::" // trim(temp_string) // "/number_of_particles")
    call delete_option("/particles/particle_array::" // trim(temp_string) // "/initial_position")
    
    call set_option("/particles/particle_array::" // trim(temp_string) // "/number_of_particles/", &
         & num_particles, stat = stat)
    
    ewrite(1,*) 'In update_particles_options'
    ewrite(1,*) num_particles
    
    assert(any(stat == (/SPUD_NO_ERROR, SPUD_NEW_KEY_WARNING/)))
    
    call set_option_attribute("/particles/particle_array::" // trim(temp_string) // "/initial_position/from_checkpoint_file/file_name", trim(filename)// "." // trim(temp_string), stat)
    
    if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
       FLAbort("Failed to set particles options filename when checkpointing particles with option path " // "/particles/particle_array::" // trim(temp_string))
    end if
    
    call set_option("/particles/particle_array::" // trim(temp_string) // "/initial_position/from_checkpoint_file/format/", trim(format), stat)
    
    if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING) then
       FLAbort("Failed to set particles options format when checkpointing particles with option path " // "/particles/particle_array")
    end if

    do i = 1, attribute_dims     
       particles_c = have_option("/particles/particle_array::" // trim(temp_string) // "/attributes/attribute["//int2str(i-1)//"]/constant")
       if (particles_c) then
          call delete_option("/particles/particle_array::" // trim(temp_string) // "/attributes/attribute["//int2str(i-1)//"]/constant")
          call set_option_attribute("/particles/particle_array::" // trim(temp_string) // "/attributes/attribute["//int2str(i-1)// &
               "]/from_checkpoint_file/file_name", trim(filename) // "." // trim(temp_string), stat)
          if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
           !  FLAbort("Failed to set scalar field particles filename when checkpointing with option path /particles/particle_array::" // &
            !      trim(temp_string) // "/attributes/attribute["//int2str(i-1)//"]/constant")
          end if
          call set_option("/particles/particle_array::" // trim(temp_string) // "/attributes/attribute["//int2str(i-1)// &
               "]/from_checkpoint_file/format/", trim(format), stat)
          if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING) then
            ! FLAbort("Failed to set scalar field particles options format when checkpointing with option path /particles/particle_array::" // &
            !      trim(temp_string) // "/attributes/attribute["//int2str(i-1)//"]/constant")
          end if
       end if
    end do

  end subroutine update_particle_options

  subroutine destroy_particles()
    !Deallocate all particle arrays (detector lists)
    if (allocated(particle_lists)) then
       deallocate(particle_lists)
    end if
    
  end subroutine destroy_particles

end module particles

  
    
    
    
