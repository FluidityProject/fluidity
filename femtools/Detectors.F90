!    Copyright (C) 2006 Imperial College London and others.
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

module detectors
  !!< Wrapper module for all things detector related

  use detector_data_types
  use detector_tools
  use detector_parallel

  !! Stuff we still need...
  use state_module
  use fields
  use field_options
  use spud
  use global_parameters, only: OPTION_PATH_LEN, PYTHON_FUNC_LEN
  use diagnostic_variables
  use integer_hash_table_module

  implicit none

  private

  public :: initialise_detectors, read_detector_move_options

  interface detector_field
     module procedure detector_field_scalar, detector_field_vector
  end interface

  character(len=OPTION_PATH_LEN), parameter :: rk_gs_path="/explicit_runge_kutta_guided_search"

contains

  function detector_field_scalar(sfield)
    !!< Return whether the supplied field should be included in the .detector file
    logical :: detector_field_scalar
    type(scalar_field), target, intent(in) :: sfield
    
    if (sfield%option_path=="".or.aliased(sfield)) then
       detector_field_scalar=.false.
    else
       detector_field_scalar = have_option(&
            trim(complete_field_path(sfield%option_path)) // &
            "/detectors/include_in_detectors")
    end if

  end function detector_field_scalar

  function detector_field_vector(vfield)
    !!< Return whether the supplied field should be included in the .detector file
    logical :: detector_field_vector
    type(vector_field), target, intent(in) :: vfield

    if (vfield%option_path=="".or.aliased(vfield)) then
       detector_field_vector=.false.
    else
       detector_field_vector = have_option(&
            trim(complete_field_path(vfield%option_path)) // &
            "/detectors/include_in_detectors")
    end if

  end function detector_field_vector

  subroutine read_detector_move_options(detector_list, options_path)
    ! Subroutine to allocate the detector parameters, 
    ! including RK stages and update vector
    type(detector_linked_list), intent(inout) :: detector_list
    character(len=*), intent(in) :: options_path

    integer :: i, j, k
    real, allocatable, dimension(:) :: stage_weights
    integer, dimension(2) :: option_rank

    if (have_option(trim(options_path))) then

       call get_option(trim(options_path)//"/subcycles",detector_list%n_subcycles)
       call get_option(trim(options_path)//"/search_tolerance",detector_list%search_tolerance)

       ! Forward Euler options
       if (have_option(trim(options_path)//"/forward_euler_guided_search")) then
          detector_list%velocity_advection = .true.
          detector_list%n_stages = 1
          allocate(detector_list%timestep_weights(detector_list%n_stages))
          detector_list%timestep_weights = 1.0
       end if

       ! Parameters for classical Runge-Kutta
       if (have_option(trim(options_path)//"/rk4_guided_search")) then
          detector_list%velocity_advection = .true.
          detector_list%n_stages = 4
          allocate(stage_weights(detector_list%n_stages*(detector_list%n_stages-1)/2))
          stage_weights = (/0.5, 0., 0.5, 0., 0., 1./)
          allocate(detector_list%stage_matrix(detector_list%n_stages,detector_list%n_stages))
          detector_list%stage_matrix = 0.
          k = 0
          do i = 1, detector_list%n_stages
             do j = 1, detector_list%n_stages
                if(i>j) then
                   k = k + 1
                   detector_list%stage_matrix(i,j) = stage_weights(k)
                end if
             end do
          end do
          allocate(detector_list%timestep_weights(detector_list%n_stages))
          detector_list%timestep_weights = (/ 1./6., 1./3., 1./3., 1./6. /)
       end if

       ! Generic Runge-Kutta options
       if (have_option(trim(options_path)//trim(rk_gs_path))) then
          detector_list%velocity_advection = .true.
          call get_option(trim(options_path)//trim(rk_gs_path)//"/n_stages",detector_list%n_stages)

          ! Allocate and read stage_matrix from options
          allocate(stage_weights(detector_list%n_stages*(detector_list%n_stages-1)/2))
          option_rank = option_shape(trim(options_path)//trim(rk_gs_path)//"/stage_weights")
          if (option_rank(2).ne.-1) then
             FLExit('Stage Array wrong rank')
          end if
          if (option_rank(1).ne.size(stage_weights)) then
             ewrite(-1,*) 'size expected was', size(stage_weights)
             ewrite(-1,*) 'size actually was', option_rank(1)
             FLExit('Stage Array wrong size')
          end if
          call get_option(trim(options_path)//trim(rk_gs_path)//"/stage_weights",stage_weights)
          allocate(detector_list%stage_matrix(detector_list%n_stages,detector_list%n_stages))
          detector_list%stage_matrix = 0.
          k = 0
          do i = 1, detector_list%n_stages
             do j = 1, detector_list%n_stages
                if(i>j) then
                   k = k + 1
                   detector_list%stage_matrix(i,j) = stage_weights(k)
                end if
             end do
          end do

          ! Allocate and read timestep_weights from options
          allocate(detector_list%timestep_weights(detector_list%n_stages))
          option_rank = option_shape(trim(options_path)//trim(rk_gs_path)//"/timestep_weights")
          if (option_rank(2).ne.-1) then
             FLExit('Timestep Array wrong rank')
          end if
          if (option_rank(1).ne.size(detector_list%timestep_weights)) then
             FLExit('Timestep Array wrong size')
          end if
          call get_option(trim(options_path)//trim(rk_gs_path)//"/timestep_weights",detector_list%timestep_weights)
       end if

       ! Boundary reflection
       if (have_option(trim(options_path)//trim("/reflect_on_boundary"))) then
          detector_list%reflect_on_boundary=.true.
       end if

       if (have_option(trim(options_path)//trim("/parametric_guided_search"))) then
          detector_list%tracking_method = GUIDED_SEARCH_TRACKING
       elseif (have_option(trim(options_path)//trim("/geometric_tracking"))) then
          detector_list%tracking_method = GEOMETRIC_TRACKING
       else
          if (check_any_lagrangian(detector_list)) then
             ewrite(-1,*) "Found lagrangian detectors, but no tracking options"
             FLExit('No lagrangian particle tracking method specified')
          end if
       end if

    else
       if (check_any_lagrangian(detector_list)) then
          ewrite(-1,*) "Found lagrangian detectors, but no timstepping options"
          FLExit('No lagrangian timestepping specified')
       end if
    end if

  end subroutine read_detector_move_options

  subroutine initialise_detectors(filename, state)
    !!< Set up the detector file headers. This has the same syntax as the
    !!< .stat file
    character(len = *), intent(in) :: filename
    type(state_type), dimension(:), intent(in) :: state

    character(len=FIELD_NAME_LEN) ::funcnam, temp_name, detector_name
    character(len=PYTHON_FUNC_LEN) :: func

    integer :: column, i, j, k, phase, m, IERROR, field_count, totaldet_global
    integer :: static_dete, python_functions_or_files, total_dete, total_dete_groups, lagrangian_dete
    integer :: python_dete, ndete, dim, str_size, type_det
    integer, dimension(2) :: shape_option
    character(len = 254) :: buffer, material_phase_name, fmt
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield, xfield
    real, allocatable, dimension(:,:) :: coords
    real, allocatable, dimension(:) :: detector_location
    real:: current_time
    character(len = OPTION_PATH_LEN) :: detectors_cp_filename, detector_file_filename

    type(detector_type), pointer :: detector
    type(element_type), pointer :: shape

    ! Idempotency check
    if (default_stat%detectors_initialised) return
    default_stat%detectors_initialised=.true.

    ewrite(2,*) "In initialise_detectors"

    ! Check whether there are actually any detectors.
    static_dete = option_count("/io/detectors/static_detector")
    lagrangian_dete = option_count("/io/detectors/lagrangian_detector")
    python_functions_or_files = option_count("/io/detectors/detector_array")
    python_dete = 0
 
    do i=1,python_functions_or_files
       write(buffer, "(a,i0,a)") "/io/detectors/detector_array[",i-1,"]"
       call get_option(trim(buffer)//"/number_of_detectors", j)
       python_dete=python_dete+j
    end do
   
    total_dete=static_dete+lagrangian_dete+python_dete
    default_stat%detector_list%total_num_det=total_dete

    total_dete_groups=static_dete+lagrangian_dete+python_functions_or_files
 
    ! Register this I/O detector list with a global list of detector lists
    call register_detector_list(default_stat%detector_list)
    default_stat%detector_list%name = trim(filename)

    allocate(default_stat%detector_group_names(total_dete_groups))
    allocate(default_stat%number_det_in_each_group(total_dete_groups))
    allocate(default_stat%detector_list%detector_names(total_dete))
    
    if (total_dete==0) return

    xfield=>extract_vector_field(state(1), "Coordinate")
    shape=>ele_shape(xfield,1)
    call get_option("/geometry/dimension",dim)
    call get_option("/timestepping/current_time", current_time)
    allocate(detector_location(dim))

    ! Enable detectors to drift with the mesh
    if (have_option("/io/detectors/move_with_mesh")) then
       default_stat%detector_list%move_with_mesh=.true.
    end if

    ! Set flag for NaN detector output
    if (have_option("/io/detectors/write_nan_outside_domain")) then
       default_stat%detector_list%write_nan_outside=.true.
    end if

    ! Always use binary format in parallel
    if (have_option("/io/detectors/binary_output") .or. isparallel()) then
       default_stat%detector_list%binary_output=.true.
    end if
    
    ! Retrieve the position of each detector. If the option
    ! "from_checkpoint_file" exists, it means we are continuing the simulation
    ! after checkpointing and the reading of the detector positions must be
    ! done from a file
    if (have_option("/io/detectors/static_detector/from_checkpoint_file").or. & 
& have_option("/io/detectors/lagrangian_detector/from_checkpoint_file").or. &
& have_option("/io/detectors/detector_array/from_checkpoint_file")) then
       default_stat%from_checkpoint=.true.
    else
       default_stat%from_checkpoint=.false.
    end if

    ! Read detectors from options
    if (.not.default_stat%from_checkpoint) then
       ewrite(2,*) "Reading detectors from options"

       ! Read all single static detector from options
       do i=1,static_dete
          write(buffer, "(a,i0,a)") "/io/detectors/static_detector[",i-1,"]"

          shape_option=option_shape(trim(buffer)//"/location")
          assert(xfield%dim==shape_option(1))
          call get_option(trim(buffer)//"/location", detector_location)

          ! The arrays below contain information about the order in which detector
          ! groups are read and how many detectors there are in each group. This is
          ! used when checkpointing detectors. In particular, when continuing a
          ! simulation from a checkpoint, with these arrays we make sure we read
          ! back the detectors from the file in the same order than at the beginning
          ! of the simulation for consistency. All the .detectors files with
          ! detector data (position, value of variables at those positions, etc.) 
          ! will have the information in the same order.
          call get_option(trim(buffer)//"/name", detector_name)
          default_stat%detector_group_names(i)=detector_name
          default_stat%number_det_in_each_group(i)=1.0
          default_stat%detector_list%detector_names(i)=detector_name

          call create_single_detector(default_stat%detector_list, xfield, &
                detector_location, STATIC_DETECTOR, trim(detector_name), i)
       end do

       ! Read all single lagrangian detector from options
       do i=1,lagrangian_dete
          write(buffer, "(a,i0,a)") "/io/detectors/lagrangian_detector[",i-1,"]"

          shape_option=option_shape(trim(buffer)//"/location")
          assert(xfield%dim==shape_option(1))
          call get_option(trim(buffer)//"/location", detector_location)

          call get_option(trim(buffer)//"/name", detector_name)
          default_stat%detector_group_names(static_dete+i)=detector_name
          default_stat%number_det_in_each_group(static_dete+i)=1.0
          default_stat%detector_list%detector_names(static_dete+i)=detector_name

          call create_single_detector(default_stat%detector_list, xfield, &
                detector_location, LAGRANGIAN_DETECTOR, trim(detector_name), static_dete+1)
       end do

       k=static_dete+lagrangian_dete+1

       do i=1,python_functions_or_files
          write(buffer, "(a,i0,a)") "/io/detectors/detector_array[",i-1,"]"

          call get_option(trim(buffer)//"/name", funcnam)
          call get_option(trim(buffer)//"/number_of_detectors", ndete)
          str_size=len_trim(int2str(ndete))
          fmt="(a,I"//int2str(str_size)//"."//int2str(str_size)//")"

          if (have_option(trim(buffer)//"/lagrangian")) then
             type_det=LAGRANGIAN_DETECTOR
          else
             type_det=STATIC_DETECTOR
          end if

          default_stat%detector_group_names(i+static_dete+lagrangian_dete)=trim(funcnam)
          default_stat%number_det_in_each_group(i+static_dete+lagrangian_dete)=ndete

          if (.not.have_option(trim(buffer)//"/from_file")) then

             ! Reading detectors from a python function
             call get_option(trim(buffer)//"/python", func)
             allocate(coords(dim,ndete))
             call set_detector_coords_from_python(coords, ndete, func, current_time)
          
             do j=1,ndete
                write(detector_name, fmt) trim(funcnam)//"_", j
                default_stat%detector_list%detector_names(k)=trim(detector_name)

                call create_single_detector(default_stat%detector_list, xfield, &
                       coords(:,j), type_det, trim(detector_name), k)
                k=k+1           
             end do
             deallocate(coords)

          else

             ! Reading from a binary file where the user has placed the detector positions
             if (getprocno() == 1) then
                 default_stat%detector_file_unit=free_unit()
                 call get_option("/io/detectors/detector_array/from_file/file_name",detector_file_filename)

#ifdef STREAM_IO
                 open(unit = default_stat%detector_file_unit, file = trim(detector_file_filename), &
                      & action = "read", access = "stream", form = "unformatted")
#else
                 FLAbort("No stream I/O support")
#endif
       
                 do j=1,ndete
                    write(detector_name, fmt) trim(funcnam)//"_", j
                    default_stat%detector_list%detector_names(k)=trim(detector_name)
                    read(default_stat%detector_file_unit) detector_location
                    call create_single_detector(default_stat%detector_list, xfield, &
                          detector_location, type_det, trim(detector_name), k)
                    k=k+1          
                 end do
              end if                
          end if
       end do      
    else 
       ewrite(2,*) "Reading detectors from checkpoint"

       ! If reading from checkpoint file:
       ! Detector checkpoint file names end in _det, with.groups appended for the header file
       ! and .positions.dat appended for the binary data file that holds the positions

       default_stat%detector_checkpoint_unit=free_unit()
       if (have_option("/io/detectors/static_detector")) then
          call get_option("/io/detectors/static_detector/from_checkpoint_file/file_name",detectors_cp_filename)  
       elseif (have_option("/io/detectors/lagrangian_detector")) then 
          call get_option("/io/detectors/lagrangian_detector/from_checkpoint_file/file_name",detectors_cp_filename)  
       else 
          call get_option("/io/detectors/detector_array/from_checkpoint_file/file_name",detectors_cp_filename)  
       end if 

       open(unit=default_stat%detector_checkpoint_unit, file=trim(detectors_cp_filename) // '.groups', action="read") 

       ! First we read the header of checkpoint_file to get the order in which the detectors were read initialliy
       do i=1, total_dete_groups 
          read(default_stat%detector_checkpoint_unit,'(a,i10)') default_stat%detector_group_names(i), default_stat%number_det_in_each_group(i)
       end do

       close(default_stat%detector_checkpoint_unit)

#ifdef STREAM_IO
       open(unit = default_stat%detector_checkpoint_unit, file = trim(detectors_cp_filename) // '.positions.dat', &
             & action = "read", access = "stream", form = "unformatted")
#else
       FLAbort("No stream I/O support")
#endif
 
       ! Read in order the last positions of the detectors from the binary file.

       do j=1,size(default_stat%detector_group_names)
          do i=1,static_dete
             write(buffer, "(a,i0,a)") "/io/detectors/static_detector[",i-1,"]"
             call get_option(trim(buffer)//"/name", temp_name)
       
             if (default_stat%detector_group_names(j)==temp_name) then
                read(default_stat%detector_checkpoint_unit) detector_location
                call create_single_detector(default_stat%detector_list, xfield, &
                      detector_location, STATIC_DETECTOR, trim(temp_name), i)                  
             else
                cycle
             end if
          end do
       end do

       do j=1,size(default_stat%detector_group_names)
          do i=1,lagrangian_dete
             write(buffer, "(a,i0,a)") "/io/detectors/lagrangian_detector[",i-1,"]"
             call get_option(trim(buffer)//"/name", temp_name)

             if (default_stat%detector_group_names(j)==temp_name) then
                read(default_stat%detector_checkpoint_unit) detector_location
                call create_single_detector(default_stat%detector_list, xfield, &
                      detector_location, LAGRANGIAN_DETECTOR, trim(temp_name), static_dete+1) 
             else
                cycle
             end if
          end do
       end do

       k=static_dete+lagrangian_dete+1

       do j=1,size(default_stat%detector_group_names) 
          do i=1,python_functions_or_files
             write(buffer, "(a,i0,a)") "/io/detectors/detector_array[",i-1,"]"      
             call get_option(trim(buffer)//"/name", temp_name)

             if (default_stat%detector_group_names(j)==temp_name) then
                call get_option(trim(buffer)//"/number_of_detectors", ndete)
                str_size=len_trim(int2str(ndete))
                fmt="(a,I"//int2str(str_size)//"."//int2str(str_size)//")"

                if (have_option(trim(buffer)//"/lagrangian")) then
                   type_det=LAGRANGIAN_DETECTOR
                else
                   type_det=STATIC_DETECTOR
                end if

                do m=1,default_stat%number_det_in_each_group(j)
                   write(detector_name, fmt) trim(temp_name)//"_", m
                   read(default_stat%detector_checkpoint_unit) detector_location
                   call create_single_detector(default_stat%detector_list, xfield, &
                          detector_location, type_det, trim(detector_name), k) 
                   k=k+1           
                end do
             else                     
                cycle                   
             end if             
          end do
       end do

    end if  ! from_checkpoint

     ! Loop over all fields in state and record the ones we want to output
     allocate (default_stat%detector_list%sfield_list(size(state)))
     allocate (default_stat%detector_list%vfield_list(size(state)))
     phaseloop: do phase=1,size(state)
        material_phase_name=trim(state(phase)%name)

        ! Count the scalar fields to include in detectors
        field_count = 0
        do i = 1, size(state(phase)%scalar_names)
           sfield => extract_scalar_field(state(phase),state(phase)%scalar_names(i))   
           if (detector_field(sfield)) field_count = field_count + 1
        end do 
        allocate(default_stat%detector_list%sfield_list(phase)%ptr(field_count))
        default_stat%detector_list%num_sfields=default_stat%detector_list%num_sfields + field_count

        ! Loop over scalar fields again to store names and create header lines
        field_count = 1
        do i=1, size(state(phase)%scalar_names)

           sfield => extract_scalar_field(state(phase),state(phase)%scalar_names(i))
           if(.not. detector_field(sfield)) then
              cycle
           end if

           ! Store name of included scalar field
           default_stat%detector_list%sfield_list(phase)%ptr(field_count)=state(phase)%scalar_names(i)
           field_count = field_count + 1
        end do

        ! Count the vector fields to include in detectors
        field_count = 0
        do i = 1, size(state(phase)%vector_names)
           vfield => extract_vector_field(state(phase),state(phase)%vector_names(i))   
           if (detector_field(vfield)) field_count = field_count + 1
        end do 
        allocate(default_stat%detector_list%vfield_list(phase)%ptr(field_count))
        default_stat%detector_list%num_vfields=default_stat%detector_list%num_vfields + field_count

        ! Loop over vector fields again to store names and create header lines
        field_count = 1
        do i=1, size(state(phase)%vector_names)

           vfield => extract_vector_field(state(phase),state(phase)%vector_names(i))
           if(.not. detector_field(vfield)) then
              cycle
           end if

           ! Store name of included vector field
           default_stat%detector_list%vfield_list(phase)%ptr(field_count)=state(phase)%vector_names(i)
           field_count = field_count + 1
        end do

     end do phaseloop

    !Get options for lagrangian detector movement
    if (check_any_lagrangian(default_stat%detector_list)) then
       call read_detector_move_options(default_stat%detector_list, "/io/detectors/lagrangian_timestepping")
    end if

    ! Write the header information into the .detectors file
    call write_detector_header(state, default_stat%detector_list)

    ! And finally some sanity checks
    totaldet_global=default_stat%detector_list%length
    call allsum(totaldet_global)
    ewrite(2,*) "Found", default_stat%detector_list%length, "local and ", totaldet_global, "global detectors"

    assert(totaldet_global==default_stat%detector_list%total_num_det)

  end subroutine initialise_detectors

  subroutine list_det_into_csr_sparsity(detector_list,ihash_sparsity,list_into_array,element_detector_list,count)
!! This subroutine creates a hash table called ihash_sparsity and a csr_sparsity matrix called element_detector_list that we use to find out 
!! how many detectors a given element has and we also obtain the location (row index) of those detectors in an array called list_into_array. 
!! This array contains the information of detector_list but in an array format, each row of the array contains the information of a detector. 
!! By accessing the array at that/those row indexes we can extract the information (position, id_number, type, etc.) of each detector present 
!! in the element (each row index corresponds to a detector)

    type(detector_linked_list), intent(inout) :: detector_list
    type(integer_hash_table), intent(inout) :: ihash_sparsity
    type(csr_sparsity), intent(inout) :: element_detector_list
    real, dimension(:,:), allocatable, intent(inout) :: list_into_array
    integer, intent(in) :: count
    
    integer, dimension(:), allocatable:: detector_count ! detectors per element
    integer, dimension(:), pointer:: detectors
    type(detector_type), pointer :: node
    type(vector_field), pointer :: vfield, xfield
    integer :: dim, i, ele, no_rows, entries, row, pos

    if (detector_list%length/=0) then
   
       node => detector_list%first

       dim=size(node%position)

       do i=1, detector_list%length

             list_into_array(i,1:dim)=node%position
             list_into_array(i,dim+1)=node%element
             list_into_array(i,dim+2)=node%id_number
             list_into_array(i,dim+3)=node%type
             list_into_array(i,dim+4)=0.0

             node => node%next
    
       end do

       ! create map between element and rows, where each row corresponds to an element 
       ! with one or more detectors
    
       ! loop over detectors: 
       ! detector is in element ele

       no_rows=count

      ! count number of detectors per row
      allocate(detector_count(1:no_rows))
      detector_count=0
      ! loop over detectors:
      node => detector_list%first

      do i=1, detector_list%length

         ele=node%element
         if  (has_key(ihash_sparsity, ele)) then
            row=fetch(ihash_sparsity, ele)
            detector_count(row)=detector_count(row)+1
         end if
         node => node%next

      end do 

      ! set up %findrm, the beginning of each row in memory
      pos=1 ! position in colm
      do row=1, no_rows
         element_detector_list%findrm(row)=pos
         pos=pos+detector_count(row)
      end do
      element_detector_list%findrm(row)=pos

      ! fill up the rows with the rom_number of the detectors in the list_into_array
      detector_count=0
      ! loop over detectors:

      do i=1, detector_list%length
      
         ele=list_into_array(i,dim+1)  
         if (has_key(ihash_sparsity, ele)) then
            row=fetch(ihash_sparsity, ele)
            detectors => row_m_ptr(element_detector_list, row)
            detector_count(row)=detector_count(row)+1
            detectors(detector_count(row))=i
         end if

      end do

      deallocate(detector_count)

    end if

  end subroutine list_det_into_csr_sparsity

end module detectors
