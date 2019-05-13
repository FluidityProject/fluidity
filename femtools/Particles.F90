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
!    You should have received a copy of the GNU Lesser General PublicS
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"


module particles
  use fldebug
  use iso_c_binding, only: C_NULL_CHAR
  use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, &
       & PYTHON_FUNC_LEN, integer_size, real_size
  use futils, only: int2str, free_unit
  use elements
  use mpi_interfaces
  use parallel_tools
  use spud
  use parallel_fields
  use fields
  use state_module
  use field_options
  use detector_data_types
  use pickers
  use detector_tools
  use detector_parallel
  use detector_move_lagrangian
  use diagnostic_variables, only: field_tag, initialise_constant_diagnostics

  use H5hut

  implicit none

  private

  public :: initialise_particles, move_particles, write_particles_loop, destroy_particles, &
            update_particle_attributes_and_fields, checkpoint_particles_loop

  type(detector_linked_list), allocatable, dimension(:), save :: particle_lists !!Particle lists with dimension equal to the number of particle subgroups

contains

  subroutine initialise_particles(filename,state)
    !!Initialise particles and set up particle file headers (per particle array)
    character(len = *), intent(in) :: filename
    type(state_type), dimension(:), intent(in) :: state

    character(len=FIELD_NAME_LEN) :: subname
    character(len = OPTION_PATH_LEN) :: group_path, subgroup_path

    type(vector_field), pointer :: xfield
    real :: current_time

    integer, dimension(3) :: attribute_size
    integer :: sub_particles
    integer :: i, k
    integer :: dim, particle_groups, total_arrays, list_counter
    integer, dimension(:), allocatable :: particle_arrays
    integer :: totaldet_global
    integer :: nfields, nprescribed, ndiagnostic, nprognostic

    logical :: from_file

    ewrite(2,*), "In initialise_particles"

    !Check whether there are any particles.
    particle_groups = option_count("/particles/particle_group")

    if (particle_groups==0) return

    !Set up particle_lists
    allocate(particle_arrays(particle_groups))
    total_arrays = 0
    do i = 1,particle_groups
       particle_arrays(i) = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup")
       total_arrays = total_arrays + particle_arrays(i)
    end do
    allocate(particle_lists(total_arrays))

    !Allocate parameters
    xfield=>extract_vector_field(state(1), "Coordinate")
    call get_option("/geometry/dimension",dim)
    call get_option("/timestepping/current_time", current_time)

    list_counter=1

    !Number of old_fields stored on particles
    nprognostic = option_count('/material_phase/scalar_field/prognostic/particles/include_in_particles/store_old_field')
    nprescribed = option_count('/material_phase/scalar_field/prescribed/particles/include_in_particles/store_old_field')
    ndiagnostic = option_count('/material_phase/scalar_field/diagnostic/particles/include_in_particles/store_old_field')
    
    do i = 1,particle_groups
       !Set group_path
       group_path = "/particles/particle_group["//int2str(i-1)//"]"
       do k = 1, particle_arrays(i)
          !subgroup_path
          subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
          ! If the option
          ! "from_file" exists, it means we are continuing the simulation
          ! after checkpointing and the reading of the particle positions must be
          ! done from a file

          from_file=have_option(trim(subgroup_path)//"/initial_position/from_file")
          call get_option(trim(subgroup_path)//"/number_of_particles", sub_particles) !number of particles in subgroup
          call get_option(trim(subgroup_path)//"/name", subname) !Name of particle subgroup
          particle_lists(list_counter)%total_num_det=sub_particles
          allocate(particle_lists(list_counter)%detector_names(sub_particles))
          !Register this I/O list with a global list of detectors/particles
          call register_detector_list(particle_lists(list_counter))
          
          !Find number of attributes 
          attribute_size(:)=0
          if (have_option(trim(subgroup_path) // "/attributes/attribute")) then
             attribute_size(1)=option_count(trim(subgroup_path) // "/attributes/attribute")
          end if

          if (option_count(trim(subgroup_path) // "/attributes/attribute/python_fields")>0) then
             attribute_size(3) = ndiagnostic+nprescribed+nprognostic
             attribute_size(2) = option_count(trim(subgroup_path) //"/attributes/attribute/python_fields/store_old_attribute")
          end if
          
          ! Enable particles to drift with the mesh
          if (have_option("/particles/move_with_mesh")) then
             particle_lists(list_counter)%move_with_mesh=.true.
          end if
          
          ! Set flag for NaN particle output
          if (have_option("/particles/write_nan_outside_domain")) then
             particle_lists(list_counter)%write_nan_outside=.true.
          end if
          
          !Read particles from options
          if (.not.from_file) then
             
             call read_particles_from_python(sub_particles, subname, current_time, state, attribute_size, xfield, dim, subgroup_path, particle_lists(list_counter))
             
          else !Read particles from file
             
             call read_particles_from_file(sub_particles, subname, attribute_size, state, xfield, dim, subgroup_path, particle_lists(list_counter))
             
          end if !Read particles

          call set_particle_output_file(sub_particles, subname, filename, attribute_size, xfield, subgroup_path, particle_lists(list_counter))
          !Get options for lagrangian particle movement
          call read_detector_move_options(particle_lists(list_counter), "/particles")
          list_counter = list_counter + 1
       end do
    end do

    ! And finally some sanity checks
    list_counter=1
    do i = 1,particle_groups
       group_path = "/particles/particle_group["//int2str(i-1)//"]"
       do k = 1,particle_arrays(i)
          subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
          call get_option(trim(subgroup_path)//"/name",subname)
          totaldet_global=particle_lists(list_counter)%length
          call allsum(totaldet_global)
          ewrite(2,*) "Found", particle_lists(list_counter)%length, "local and ", totaldet_global, "global particles for particle array ", trim(subname)
          
          assert(totaldet_global==particle_lists(list_counter)%total_num_det)
          list_counter = list_counter + 1
       end do
    end do

    deallocate(particle_arrays)
    
  end subroutine initialise_particles

  subroutine read_particles_from_python(sub_particles, subname, current_time, state, attribute_size, xfield, dim, subgroup_path, p_list)
    ! Reading particles from a python function

    type(state_type), dimension(:), intent(in) :: state
    type(vector_field), pointer, intent(in) :: xfield
    type(detector_linked_list), intent(inout) :: p_list
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path
    character(len=FIELD_NAME_LEN), intent(in) :: subname
    real, intent(in) :: current_time
    integer, intent(in) :: sub_particles, dim
    integer, dimension(3), intent(in) :: attribute_size

    character(len=PYTHON_FUNC_LEN) :: func
    character(len = FIELD_NAME_LEN) :: fmt
    character(len=FIELD_NAME_LEN) :: particle_name
    real, allocatable, dimension(:,:) :: coords !array to hold coordinates of particles for initialisation
    integer :: l, str_size

    ewrite(2,*) "Reading particles from options"
    call get_option(trim(subgroup_path)//"/initial_position/python", func)
    allocate(coords(dim,sub_particles))
    call set_detector_coords_from_python(coords, sub_particles, func, current_time)

    str_size=len_trim(int2str(sub_particles))
    fmt="(a,I"//int2str(str_size)//"."//int2str(str_size)//")"

    do l=1,sub_particles
       write(particle_name, fmt) trim(subname)//"_", l
       call create_single_particle(p_list, xfield, coords(:,l), &
            l, LAGRANGIAN_DETECTOR, trim(particle_name), attribute_size)
    end do
    if (attribute_size(1)/=0) then
       call update_particle_subgroup_attributes_and_fields(state, xfield, current_time, subgroup_path, p_list)
    end if
    
    deallocate(coords)

  end subroutine read_particles_from_python

  subroutine read_particles_from_file(sub_particles, subname, attribute_size, state, xfield, dim, subgroup_path, p_list)
    ! If reading from file:
    ! Particles checkpoint file names end in _par, with.groups appended for the header file
    ! and .attributes.dat appended for the binary data file that holds the positions and attributes

    type(state_type), dimension(:), intent(in) :: state
    type(vector_field), pointer, intent(in) :: xfield
    type(detector_linked_list), intent(inout) :: p_list
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path
    character(len=FIELD_NAME_LEN), intent(in) :: subname
    integer, intent(in) :: sub_particles, dim
    integer, dimension(3), intent(in) :: attribute_size

    real, allocatable, dimension(:,:) :: attribute_vals !array to hold particle attribute values for initialisation
    real, allocatable, dimension(:) :: positions !array to hold particle coordinates if from checkpoint file

    integer :: h5_ierror, m, i, j, str_size, old_attrib, old_field
    integer(kind=8) :: h5_id
    character(len=OPTION_PATH_LEN) :: particles_cp_filename
    character(len=FIELD_NAME_LEN) :: attname, particle_name, fmt
    type(scalar_field), pointer :: sfield

    ewrite(2,*) "Reading particles from file"

    allocate(positions(dim))
    allocate(attribute_vals(3,maxval(attribute_size)))

    str_size = len_trim(int2str(sub_particles))
    fmt = "(a,I"//int2str(str_size)//"."//int2str(str_size)//")"

    call get_option(trim(subgroup_path) // "/initial_position/from_file/file_name", particles_cp_filename)

    h5_id = h5_openfile(trim(particles_cp_filename), H5_O_RDONLY, H5_PROP_DEFAULT)
    h5_ierror = h5_setstep(h5_id, int(1, 8))

    do m = 1, sub_particles
      write(particle_name, fmt) trim(subname)//"_", m

      ! set view for this particle
      h5_ierror = h5pt_setview(h5_id, int(m, 8), int(m, 8))
      if (dim >= 1) &
           h5_ierror = h5pt_readdata_r8(h5_id, "x", positions(1))
      if (dim >= 2) &
           h5_ierror = h5pt_readdata_r8(h5_id, "y", positions(2))
      if (dim >= 3) &
           h5_ierror = h5pt_readdata_r8(h5_id, "z", positions(3))

      old_attrib = 0
      ! read out attributes by name
      if (attribute_size(1) /= 0) then
        do i = 1, attribute_size(1)
          call get_option(trim(subgroup_path) // "/attributes/attribute["//int2str(i-1)//"]/name", attname)
          h5_ierror = h5pt_readdata_r8(h5_id, trim(attname), attribute_vals(1, i))

          ! read old attribute if required
          if (have_option(trim(subgroup_path) // "/attributes/attribute["//int2str(i-1)//"]/python_fields/store_old_attribute")) then
            old_attrib = old_attrib + 1
            h5_ierror = h5pt_readdata_r8(h5_id, "old"//trim(attname), attribute_vals(2, old_attrib))
          end if
        end do
      end if
      assert(old_attrib == attribute_size(2))

      ! read in fields
      old_field = 0
      if (attribute_size(3) /= 0) then
        do i = 1, size(state)
          do j = 1, size(state(i)%scalar_names)
            sfield => extract_scalar_field(state(i), state(i)%scalar_names(j))
            if (sfield%option_path == "" .or. aliased(sfield)) then
              cycle
            else if (have_option(trim(complete_field_path(sfield%option_path)) // "/particles/include_in_particles/store_old_field")) then
              old_field = old_field + 1
              h5_ierror = h5pt_readdata_r8(h5_id, "field"//trim(state(i)%scalar_names(j)), attribute_vals(3, old_field))
            end if
          end do
        end do
      end if
      assert(old_field == attribute_size(3))

      call create_single_particle(p_list, xfield, &
           positions, m, LAGRANGIAN_DETECTOR, trim(particle_name), &
           attribute_size, attribute_vals=attribute_vals)
     end do

     deallocate(attribute_vals)
     deallocate(positions)
  end subroutine read_particles_from_file

  subroutine set_particle_output_file(sub_particles, subname, filename, attribute_size, xfield, subgroup_path, p_list)
    !Sets up the particle output file type and output file.

    type(vector_field), pointer, intent(in) :: xfield
    type(detector_linked_list), intent(inout) :: p_list
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path
    character(len = *), intent(in) :: filename
    integer, intent(in) :: sub_particles
    integer, dimension(3), intent(in) :: attribute_size
    character(len=FIELD_NAME_LEN), intent(in) :: subname

    integer :: file_id
    p_list%h5_id = h5_openfile(trim(filename) // '.particles.' // trim(subname) // '.h5part', H5_O_WRONLY, H5_PROP_DEFAULT)

    ! optionally set any file attributes here?

  end subroutine set_particle_output_file

  subroutine create_single_particle(detector_list, xfield, position, id, type, name, attribute_size, attribute_vals)
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
    integer, dimension(3), intent(in) :: attribute_size

    real ::  dt

    shape=>ele_shape(xfield,1)
    assert(xfield%dim+1==local_coord_count(shape))
    detector_list%detector_names(id)=name
    ! Determine element and local_coords from position
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

    allocate(detector%attributes(attribute_size(1)))
    allocate(detector%old_attributes(attribute_size(2)))
    allocate(detector%old_fields(attribute_size(3)))

    if (present(attribute_vals)) then
       detector%attributes = attribute_vals(1,1:attribute_size(1))
       detector%old_attributes = attribute_vals(2,1:attribute_size(2))
       detector%old_fields = attribute_vals(3,1:attribute_size(3))
    else
       detector%attributes(:) = 0
       detector%old_attributes(:) = 0
       detector%old_fields(:) = 0
    end if
    
  end subroutine create_single_particle

  subroutine move_particles(state, dt, timestep)
    !!Routine to loop over particle arrays and call move_lagrangian_detectors
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: dt
    integer, intent(in) :: timestep

    integer, dimension(3) :: attribute_size
    integer :: particle_groups, list_counter
    integer, dimension(:), allocatable :: particle_arrays
    integer :: i, m, k
    integer :: nprescribed, ndiagnostic, nprognostic
    character(len = OPTION_PATH_LEN) :: group_path, subgroup_path

    !Check whether there are any particles.

    particle_groups = option_count("/particles/particle_group")
    if (particle_groups==0) return

    !Number of old_fields stored on particles
    nprognostic = option_count('/material_phase/scalar_field/prognostic/particles/include_in_particles/store_old_field')
    nprescribed = option_count('/material_phase/scalar_field/prescribed/particles/include_in_particles/store_old_field')
    ndiagnostic = option_count('/material_phase/scalar_field/diagnostic/particles/include_in_particles/store_old_field')

    !Set up particle_lists
    allocate(particle_arrays(particle_groups))
    particle_arrays(:) = 0
    do i = 1,particle_groups
       particle_arrays(i) = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup")
    end do
    ewrite(2,*), "In move_particles"
    list_counter = 1
    do i = 1, particle_groups
       group_path = "/particles/particle_group["//int2str(i-1)//"]"
       do k = 1, particle_arrays(i)
          subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
          attribute_size(1)=0
          if (have_option(trim(subgroup_path) // "/attributes/attribute")) then
             attribute_size(1)=option_count(trim(subgroup_path) // "/attributes/attribute")
          end if
          attribute_size(2)=0
          attribute_size(3)=0
          do m = 1,attribute_size(1)
             if (have_option(trim(subgroup_path) // "/attributes/attribute["//int2str(m-1)//"]/python_fields")) then
                attribute_size(3) = ndiagnostic + nprescribed + nprognostic
                if (have_option(trim(subgroup_path) // "/attributes/attribute["//int2str(m-1)//"]/python_fields/store_old_attribute")) then
                   attribute_size(2)=attribute_size(2)+1
                end if
             end if
          end do
          call move_lagrangian_detectors(state, particle_lists(list_counter), dt, timestep, attribute_size)
          list_counter = list_counter + 1
       end do
    end do

  end subroutine move_particles

  subroutine update_particle_subgroup_attributes_and_fields(state, xfield, time, subgroup_path, p_list)
    !!Routine to set particle attributes 
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: time
    type(vector_field), pointer, intent(in) :: xfield
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path
    type(detector_linked_list), intent(in) :: p_list
    
    character(len=PYTHON_FUNC_LEN) :: func
    character(len=FIELD_NAME_LEN), allocatable, dimension(:) :: field_name
    type(detector_type), pointer :: particle

    real, allocatable, dimension(:,:) :: positions
    real, allocatable, dimension(:,:) :: attribute_array
    real, allocatable, dimension(:,:) :: old_attributes
    character, allocatable, dimension(:,:) :: old_att_names
    character(len = OPTION_PATH_LEN) :: old_name

    real :: constant
    integer :: j, nparticles, l, m, n
    integer, allocatable, dimension(:) :: store_old_att
    real, allocatable, dimension(:,:) :: lcoords
    integer, allocatable, dimension(:) :: ele

    !Check if this processor contains particles
    nparticles = p_list%length

    if (nparticles==0) then
       return
    end if

    !Set parameters to calculate attributes
    particle => p_list%first
    allocate(positions(size(particle%position),nparticles))
    allocate(attribute_array(size(particle%attributes),nparticles))
    allocate(lcoords(size(particle%local_coords),nparticles))
    allocate(ele(nparticles))
    allocate(old_attributes(size(particle%old_attributes),nparticles))
    allocate(old_att_names(FIELD_NAME_LEN,size(particle%old_attributes)))

    particle => p_list%first
    do j = 1,nparticles
       positions(:,j) = particle%position
       lcoords(:,j) = particle%local_coords
       ele(j) = particle%element
       old_attributes(:,j) = particle%old_attributes
       particle => particle%next
    end do

    l=1
    particle => p_list%first
    do n = 1,size(particle%attributes)
       if (have_option(trim(subgroup_path) // '/attributes/attribute['//int2str(n-1)//']/python_fields/store_old_attribute')) then
          call get_option(trim(subgroup_path) // '/attributes/attribute['//int2str(n-1)//']/name', old_name)
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

    do n = 0,size(particle%attributes)-1
       if (have_option(trim(subgroup_path) // '/attributes/attribute['//int2str(n)//']/constant')) then
          call get_option(trim(subgroup_path) // '/attributes/attribute['//int2str(n)//']/constant', constant)
          attribute_array(n+1,:) = constant
       else if (have_option(trim(subgroup_path) // '/attributes/attribute['//int2str(n)//']/python')) then
          call get_option(trim(subgroup_path) // '/attributes/attribute['//int2str(n)//']/python', func)
          call set_particle_attribute_from_python(attribute_array(n+1,:), positions(:,:), nparticles, func, time)
       else if (have_option(trim(subgroup_path) // '/attributes/attribute['//int2str(n)//']/python_fields')) then
          call get_option(trim(subgroup_path) // '/attributes/attribute['//int2str(n)//']/python_fields', func)
          call set_particle_attribute_from_python_fields(p_list, state, xfield, positions(:,:), lcoords(:,:), ele(:), nparticles, &
               & attribute_array(n+1,:), old_att_names, old_attributes, func, time)
       else if (have_option(trim(subgroup_path) // '/attributes/attribute['//int2str(n)//']/from_checkpoint_file')) then
          particle => p_list%first
          attribute_array(n+1,:) = particle%attributes(n+1)
       end if
    end do
    
    !Set attribute values and old_attribute values
    particle => p_list%first
    if (size(particle%old_attributes)==0) then
       do j = 1,nparticles
          particle%attributes = attribute_array(:,j)
          particle => particle%next
       end do
    else
       allocate(store_old_att(size(particle%attributes)))
       do n=1, size(particle%attributes)
          if (have_option(trim(subgroup_path) // '/attributes/attribute['//int2str(n-1)//']/python_fields/store_old_attribute')) then
             store_old_att(n)=1
          else
             store_old_att(n)=0
          end if
       end do
       do j = 1,nparticles
          particle%attributes = attribute_array(:,j)
          m=1
          do n = 1,size(particle%attributes)
             if (store_old_att(n)==1) then
                particle%old_attributes(m) = particle%attributes(n)
                m=m+1
             end if
          end do
          particle => particle%next
       end do
    end if

    particle => p_list%first
    if (size(particle%old_fields)/=0) then
       call update_particle_subgroup_fields(state, ele, lcoords, p_list)
    end if

    deallocate(positions)
    deallocate(lcoords)
    deallocate(ele)
    deallocate(attribute_array)
    deallocate(old_attributes)
    deallocate(old_att_names)

  end subroutine update_particle_subgroup_attributes_and_fields

  subroutine update_particle_attributes_and_fields(state, time)
    !!Routine to loop over particle arrays and update particle attributes
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: time
    type(vector_field), pointer :: xfield
    type(detector_type), pointer :: particle
    character(len = OPTION_PATH_LEN) :: group_path, subgroup_path

    integer :: i, k
    integer :: particle_groups, list_counter
    integer, dimension(:), allocatable :: particle_arrays

    !Check whether there are any particles.

    particle_groups = option_count("/particles/particle_group")
    if (particle_groups==0) return

    !Set up particle_lists
    allocate(particle_arrays(particle_groups))
    particle_arrays(:) = 0
    do i = 1,particle_groups
       particle_arrays(i) = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup")
    end do

    ewrite(2,*), "In update_particle_attributes_and_fields"

    !Allocate parameters
    xfield=>extract_vector_field(state(1), "Coordinate")
    list_counter = 1

    !Update particle attributes by array
    do i = 1,particle_groups
       group_path = "/particles/particle_group["//int2str(i-1)//"]"
       do k = 1, particle_arrays(i)
          subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
          if (particle_lists(list_counter)%length==0) then
             list_counter = list_counter + 1
             cycle
          end if
          particle => particle_lists(list_counter)%first
          if (size(particle%attributes)/=0) then
             call update_particle_subgroup_attributes_and_fields(state, xfield, time, subgroup_path, particle_lists(list_counter))
          end if
          list_counter = list_counter + 1
       end do
    end do
  end subroutine update_particle_attributes_and_fields

  subroutine update_particle_subgroup_fields(state, ele, lcoords, p_list)
    
    type(state_type), dimension(:), intent(in) :: state
    real, dimension(:,:), intent(in) :: lcoords
    integer, dimension(:), intent(in) :: ele
    type(detector_linked_list), intent(in) :: p_list

    character(len = OPTION_PATH_LEN) :: name
    type(scalar_field), pointer :: sfield
    real, allocatable, dimension(:,:) :: old_field_vals
    real :: value
    type(detector_type), pointer :: particle
    integer :: phase, f, l, j

    particle => p_list%first
    allocate(old_field_vals(size(particle%old_fields),size(lcoords(1,:))))
    l=1
    do phase=1,size(state)
       do f = 1, size(state(phase)%scalar_names)
          sfield => extract_scalar_field(state(phase),state(phase)%scalar_names(f))
          if (sfield%option_path=="".or.aliased(sfield)) then
             cycle
          else if (have_option(trim(complete_field_path(sfield%option_path)) // "/particles/include_in_particles/store_old_field")) then
             do j = 1,size(lcoords(1,:))
                value = eval_field(ele(j), sfield, lcoords(:,j))
                old_field_vals(l,j)=value
             end do
             l=l+1
          end if
       end do
    end do

    particle => p_list%first
    do j = 1,size(lcoords(1,:))
       particle%old_fields=old_field_vals(:,j)
       particle=>particle%next
    end do

    deallocate(old_field_vals)

  end subroutine update_particle_subgroup_fields

  subroutine write_particles_loop(state, timestep, time)
    !!Subroutine to loop over particle_lists and call write_particles for each list
    type(state_type), dimension(:), intent(in) :: state
    integer, intent(in) :: timestep
    real, intent(in) :: time

    integer :: attribute_dims
    integer :: i, k
    integer :: particle_groups, particle_subgroups, list_counter
    character(len=OPTION_PATH_LEN) :: group_path, subgroup_path

    !Check whether there are any particles.

    particle_groups = option_count("/particles/particle_group")
    if (particle_groups==0) return

    ewrite(1,*) "In write_particles_loop"
    
    list_counter = 1
    do i = 1, particle_groups
      group_path = "/particles/particle_group["//int2str(i-1)//"]"
      particle_subgroups = option_count(trim(group_path) // "/particle_subgroup")
      do k = 1, particle_subgroups
        subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
        attribute_dims=option_count(trim(subgroup_path) // '/attributes/attribute')
        call write_particles_subgroup(state, particle_lists(list_counter), attribute_dims, timestep, time, subgroup_path)
        list_counter = list_counter + 1
      end do
    end do
  end subroutine write_particles_loop

  subroutine write_particles_subgroup(state, detector_list, attribute_dims, timestep, time, subgroup_path)
    !!< Write values of particles to the previously opened particles file.
    type(state_type), dimension(:), intent(in) :: state
    type(detector_linked_list), intent(inout) :: detector_list
    integer, intent(in) :: timestep
    real, intent(in) :: time
    integer, intent(in) :: attribute_dims !dimensions of particles attribute information carried (attributes at current timestep, field values and attribute values at previous timestep)
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path

    integer :: h5_ierror, dim, i
    real, dimension(:,:), allocatable :: positions, attrib_data
    type(vector_field), pointer :: vfield
    type(detector_type), pointer :: node
    character(len=FIELD_NAME_LEN) :: attname

    ewrite(1,*) "In write_particles"

    ! create new step
    h5_ierror = h5_setstep(detector_list%h5_id, int(timestep, 8))

    ! write time as a step attribute
    h5_ierror = h5_writestepattrib_r8(detector_list%h5_id, "time", [time], int(1, 8))
    
    ! set the number of particles this process is going to write
    h5_ierror = h5pt_setnpoints(detector_list%h5_id, int(detector_list%length, 8))

    ! set up arrays to hold all node data (this won't work with large numbers of particles)
    vfield => extract_vector_field(state, "Coordinate")
    dim = vfield%dim
    allocate(positions(detector_list%length, 3))
    allocate(attrib_data(detector_list%length, attribute_dims))

    node => detector_list%first
    position_loop: do i = 1, detector_list%length
      assert(size(node%position) == dim)
      assert(size(node%attributes) == attribute_dims)

      positions(i,1:dim) = node%position(:)
      attrib_data(i,:) = node%attributes(:)

      node => node%next
    end do position_loop

    ! write out position
    if (dim >= 1) &
         h5_ierror = h5pt_writedata_r8(detector_list%h5_id, "x", positions(:,1))
    if (dim >= 2) &
         h5_ierror = h5pt_writedata_r8(detector_list%h5_id, "y", positions(:,2))
    if (dim >= 3) then
      h5_ierror = h5pt_writedata_r8(detector_list%h5_id, "z", positions(:,3))
    else
      positions(:,3) = 0.
      h5_ierror = h5pt_writedata_r8(detector_list%h5_id, "z", positions(:,3))
    end if

    ! write out attributes
    attribute_loop: do i = 1, attribute_dims
      call get_option(trim(subgroup_path) // "/attributes/attribute["//int2str(i-1)//"]/name", attname)
      h5_ierror = h5pt_writedata_r8(detector_list%h5_id, trim(attname), attrib_data(:,i))
    end do attribute_loop

    deallocate(attrib_data)
    deallocate(positions)
  end subroutine write_particles_subgroup

  subroutine checkpoint_particles_loop(state,prefix,postfix,cp_no)
    !!Subroutine to loop over particle_lists and call checkpoint_particles for each list
    type(state_type), dimension(:), intent(in) :: state
    character(len = *), intent(in) :: prefix
    character(len = *), intent(in) :: postfix
    integer, optional, intent(in) :: cp_no !Checkpoint number of the simulation
    
    character(len = OPTION_PATH_LEN) :: lpostfix
    character(len=OPTION_PATH_LEN) :: group_path, subgroup_path, subgroup_path_name, name

    integer, dimension(3) :: attribute_size
    integer :: nprescribed, ndiagnostic, nprognostic
    integer, dimension(:), allocatable :: particle_arrays
    integer :: i, m, k, particle_groups, list_counter

    !Check whether there are any particles.

    particle_groups = option_count("/particles/particle_group")
    if (particle_groups==0) return

    !Set up particle_lists
    allocate(particle_arrays(particle_groups))
    particle_arrays(:) = 0
    do i = 1,particle_groups
       particle_arrays(i) = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup")
    end do

    ewrite(1, *) "Checkpointing particles"

    assert(len_trim(prefix) > 0)

    lpostfix = postfix

    !Number of old_fields stored on particles
    nprognostic = option_count('/material_phase/scalar_field/prognostic/particles/include_in_particles/store_old_field')
    nprescribed = option_count('/material_phase/scalar_field/prescribed/particles/include_in_particles/store_old_field')
    ndiagnostic = option_count('/material_phase/scalar_field/diagnostic/particles/include_in_particles/store_old_field')

    list_counter = 1
    do i = 1, particle_groups
       group_path = "/particles/particle_group["//int2str(i-1)//"]"
       do k = 1, particle_arrays(i)
          subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
          subgroup_path_name = trim(group_path) // "/particle_subgroup::"!!option path used in update_particle_options

          ! count number of attributes in this subgroup
          attribute_size(1)=0
          if (have_option(trim(subgroup_path) // '/attributes/attribute')) then
             attribute_size(1)=option_count(trim(subgroup_path) // '/attributes/attribute')
           end if

           ! count number of required old attributes (for fields)
           ! if there are any fields, save the fields included in particles as well
          attribute_size(2)=0
          attribute_size(3)=0
          do m = 1,attribute_size(1)
             if (have_option(trim(subgroup_path) // '/attributes/attribute['//int2str(m-1)//']/python_fields')) then
                attribute_size(3) = ndiagnostic + nprescribed + nprognostic
                if (have_option(trim(subgroup_path) // '/attributes/attribute['//int2str(m-1)//']/python_fields/store_old_attribute')) then
                   attribute_size(2)=attribute_size(2)+1
                end if
             end if
          end do
          call get_option(trim(subgroup_path) // "/name", name)
          call checkpoint_particles_subgroup(state,prefix,lpostfix,cp_no,particle_lists(list_counter),attribute_size,name,subgroup_path,subgroup_path_name)
          list_counter = list_counter + 1
       end do
    end do

  end subroutine checkpoint_particles_loop

  subroutine checkpoint_particles_subgroup(state,prefix,lpostfix,cp_no,particle_list,attribute_size,name,subgroup_path,subgroup_path_name)
    !!<Checkpoint Particles

    type(state_type), dimension(:), intent(in) :: state
    character(len = *), intent(in) :: prefix
    character(len = *), intent(in) :: lpostfix
    integer, optional, intent(in) :: cp_no !Checkpoint number of the simulation
    character(len = *), intent(in) :: name
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path, subgroup_path_name

    type(detector_linked_list), intent(inout) :: particle_list
    integer, dimension(3), intent(in) :: attribute_size

    character(len=OPTION_PATH_LEN) :: particles_cp_filename
    character(len=FIELD_NAME_LEN) :: attname
    integer :: h5_ierror, i, j, dim, old_attrib, old_field
    integer(kind=8) :: h5_id
    real, dimension(:,:), allocatable :: positions, attrib_data, old_attrib_data, old_field_data
    type(vector_field), pointer :: vfield
    type(scalar_field), pointer :: sfield
    type(detector_type), pointer :: node

    ! construct a new particle checkpoint filename
    particles_cp_filename = trim(prefix)
    if(present(cp_no)) particles_cp_filename = trim(particles_cp_filename) // "_" // int2str(cp_no)
    particles_cp_filename = trim(particles_cp_filename) // "_" // trim(lpostfix)
    particles_cp_filename = trim(particles_cp_filename) // "_particles." // trim(name) // ".h5part"

    ! open output file
    h5_id = h5_openfile(trim(particles_cp_filename), H5_O_WRONLY, H5_PROP_DEFAULT)
    h5_ierror = h5_setstep(h5_id, int(1, 8))
    h5_ierror = h5pt_setnpoints(h5_id, int(particle_list%length, 8))

    ! get dimension of particle positions
    vfield => extract_vector_field(state(1), "Coordinate")
    dim = vfield%dim

    ! allocate arrays for node data
    allocate(positions(particle_list%length, dim))
    allocate(attrib_data(particle_list%length, attribute_size(1)))
    allocate(old_attrib_data(particle_list%length, attribute_size(2)))
    allocate(old_field_data(particle_list%length, attribute_size(3)))

    node => particle_list%first
    positionloop_cp: do i = 1, particle_list%length
      ! collect positions
      assert(size(node%position) == dim)

      positions(i,:) = node%position(:)
      if (attribute_size(1) /= 0) &
           attrib_data(i,:) = node%attributes(:)
      if (attribute_size(2) /= 0) &
           old_attrib_data(i,:) = node%old_attributes(:)
      if (attribute_size(3) /= 0) &
           old_field_data(i,:) = node%old_fields(:)

      node => node%next
    end do positionloop_cp

    ! write out position
    if (dim >= 1) &
         h5_ierror = h5pt_writedata_r8(h5_id, "x", positions(:,1))
    if (dim >= 2) &
         h5_ierror = h5pt_writedata_r8(h5_id, "y", positions(:,2))
    if (dim >= 3) &
         h5_ierror = h5pt_writedata_r8(h5_id, "z", positions(:,3))

    old_attrib = 0

    if (attribute_size(1) /= 0) then
      attribute_loop: do i = 1, attribute_size(1)
        call get_option(trim(subgroup_path) // "/attributes/attribute["//int2str(i-1)//"]/name", attname)
        h5_ierror = h5pt_writedata_r8(h5_id, trim(attname), attrib_data(:,i))

        ! collapse in old attribute loop
        if (have_option(trim(subgroup_path) // "/attributes/attribute["//int2str(i-1)//"]/python_fields/store_old_attribute")) then
          old_attrib = old_attrib + 1
          h5_ierror = h5pt_writedata_r8(h5_id, "old"//trim(attname), old_attrib_data(:,old_attrib))
        end if
      end do attribute_loop
    end if
    assert(old_attrib == attribute_size(2))

    old_field = 0
    if (attribute_size(3) /= 0) then
      do i = 1, size(state)
        do j = 1, size(state(i)%scalar_names)
          sfield => extract_scalar_field(state(i), state(i)%scalar_names(j))
          if (sfield%option_path == "" .or. aliased(sfield)) then
            cycle
          else if (have_option(trim(complete_field_path(sfield%option_path)) // "/particles/include_in_particles/store_old_field")) then
            old_field = old_field + 1
            h5_ierror = h5pt_writedata_r8(h5_id, "field"//trim(state(i)%scalar_names(j)), old_field_data(:,old_field))
          end if
        end do
      end do
    end if
    assert(old_field == attribute_size(3))

    ! update schema file to read this subgroup from the checkpoint file
    call update_particle_subgroup_options(trim(particles_cp_filename), "binary", particle_list, name, attribute_size(1), subgroup_path_name)

    deallocate(old_field_data)
    deallocate(old_attrib_data)
    deallocate(attrib_data)
    deallocate(positions)

    h5_ierror = h5_closefile(h5_id)

  end subroutine checkpoint_particles_subgroup

  subroutine update_particle_subgroup_options(filename, format, particle_list, name, attribute_dims, subgroup_path_name)
    !! Updates the initial options of particles in the schema file for reinitialization after checkpointing.
    !! Updates schema options for the initial number of particles and their initial positions. 

    character(len = *), intent(in) :: filename
    character(len = *), intent(in) :: format
    character(len = *), intent(in) :: name

    type(detector_linked_list), intent(inout) :: particle_list
    integer, intent(in) :: attribute_dims
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path_name

    integer :: num_particles, j, stat
    logical :: particles_c

    character(len = 254) :: temp_string

    num_particles = particle_list%total_num_det

    temp_string=name

    ewrite(1,*) 'In update_particles_options'
    ewrite(1,*) temp_string
    
    call delete_option(trim(subgroup_path_name) // trim(temp_string) // "/number_of_particles")
    call delete_option(trim(subgroup_path_name) // trim(temp_string) // "/initial_position")
    
    call set_option(trim(subgroup_path_name) // trim(temp_string) // "/number_of_particles/", &
         & num_particles, stat = stat)
    
    assert(any(stat == (/SPUD_NO_ERROR, SPUD_NEW_KEY_WARNING/)))
    
    call set_option_attribute(trim(subgroup_path_name) // trim(temp_string) // "/initial_position/from_file/file_name", trim(filename), stat)
    
    if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
       FLAbort("Failed to set particles options filename when checkpointing particles with option path " // "/particles/particle_array::" // trim(temp_string))
    end if
    
    call set_option(trim(subgroup_path_name) // trim(temp_string) // "/initial_position/from_file/format/", trim(format), stat)
    
    if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING) then
       FLAbort("Failed to set particles options format when checkpointing particles with option path " // "/particles/particle_group")
    end if
    
    do j = 1, attribute_dims     
       particles_c = have_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/attribute["//int2str(j-1)//"]/constant")
       if (particles_c) then
          call delete_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/attribute["//int2str(j-1)//"]/constant")
          call set_option_attribute(trim(subgroup_path_name) // trim(temp_string) // "/attributes/attribute["//int2str(j-1)// &
               "]/from_checkpoint_file/file_name", trim(filename) // "." // trim(temp_string), stat)
          if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
             FLAbort("Failed to set scalar field particles filename when checkpointing")
          end if
          call set_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/attribute["//int2str(j-1)// &
               "]/from_checkpoint_file/format/", trim(format), stat)
          if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING) then
             FLAbort("Failed to set scalar field particles options format when checkpointing")
          end if
       end if
    end do
    
  end subroutine update_particle_subgroup_options

  subroutine destroy_particles()
    integer :: i, particle_groups, h5_ierror

    if (allocated(particle_lists)) then
      ! gracefully clean up output files
      particle_groups = size(particle_lists)
      do i = 1, particle_groups
        h5_ierror = h5_closefile(particle_lists(i)%h5_id)
      enddo

      !Deallocate all particle arrays (detector lists)
      deallocate(particle_lists)
    end if
    
  end subroutine destroy_particles

end module particles

  
    
    
    
