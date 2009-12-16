! This module loads NEMO data into the states specified in the flml
module nemo_states_module

use NEMO_load_fields_vars
use spud
use fields
use state_module
use boundary_conditions
use global_parameters, only: OPTION_PATH_LEN, pi, current_debug_level
use coordinates
use Field_Options

implicit none

character(len=FIELD_NAME_LEN), dimension(:), pointer, save :: nemo_scalar_field_names
character(len=FIELD_NAME_LEN), dimension(:), pointer, save :: nemo_vector_field_names
integer, save :: no_nemo_scalar_fields=0
integer, save :: no_nemo_vector_fields=0

contains

subroutine insert_nemo_scalar_field(field)

  type(scalar_field), intent(in) :: field

  character(len=FIELD_NAME_LEN), dimension(:), pointer:: prev_field_names

  if (no_nemo_scalar_fields>0) then
    prev_field_names => nemo_scalar_field_names
  end if
  allocate(nemo_scalar_field_names(1:no_nemo_scalar_fields+1))
  nemo_scalar_field_names(no_nemo_scalar_fields+1)=field%name
  if (no_nemo_scalar_fields>0) then
    nemo_scalar_field_names(1:no_nemo_scalar_fields)=prev_field_names
    deallocate(prev_field_names)
  end if

  no_nemo_scalar_fields = no_nemo_scalar_fields+1
  
end subroutine

subroutine insert_nemo_vector_field(field)

  type(vector_field), intent(in) :: field

  character(len=FIELD_NAME_LEN), dimension(:), pointer:: prev_field_names

  if (no_nemo_vector_fields>0) then
    prev_field_names => nemo_vector_field_names
  end if
  allocate(nemo_vector_field_names(1:no_nemo_vector_fields+1))
  nemo_vector_field_names(no_nemo_vector_fields+1)=field%name
  if (no_nemo_vector_fields>0) then
    nemo_vector_field_names(1:no_nemo_vector_fields)=prev_field_names
  end if

  no_nemo_vector_fields = no_nemo_vector_fields+1
  
end subroutine

subroutine set_nemo_fields(state)

  type(state_type), intent(in) :: state

  logical, save:: first_time=.true.
  integer :: i
  type(scalar_field), pointer:: sfield
  type(vector_field), pointer:: vfield
  character(len=OPTION_PATH_LEN) :: format

  call load_nemo_values(state)

  if (first_time) then
    do i=1, no_nemo_scalar_fields

      sfield => extract_scalar_field(state, nemo_scalar_field_names(i))

      if (have_option(trim(sfield%option_path)//"/prognostic")) then

        call get_option(trim(sfield%option_path) // "/prognostic/initial_condition/NEMO_data/format", format)

        select case (format)
          case ("Temperature")
            call remap_field(temperature_t, sfield)
          case ("Salinity")
            call remap_field(salinity_t, sfield)
          case ("Free-surface height")
            call remap_field(pressure_t, sfield)
        end select
      
      endif
    enddo
    do i=1, no_nemo_vector_fields

      vfield => extract_vector_field(state, nemo_vector_field_names(i))

      if (have_option(trim(vfield%option_path)//"/prognostic")) then

        call get_option(trim(vfield%option_path) // "/prognostic/initial_condition/NEMO_data/format", format)

        select case (format)
          case ("Velocity")
            call remap_field(velocity_t, vfield)
        end select
      
      endif
    enddo
    first_time=.false.
  end if

  do i=1, no_nemo_scalar_fields

    sfield => extract_scalar_field(state, nemo_scalar_field_names(i))

    if (have_option(trim(sfield%option_path)//"/prescribed")) then
      
      call get_option(trim(sfield%option_path) // "/prescribed/value/NEMO_data/format", format)

      select case (format)
        case ("Temperature")
          call remap_field(temperature_t, sfield)
        case ("Salinity")
          call remap_field(salinity_t, sfield)
        case ("Free-surface height")
          call remap_field(pressure_t, sfield)
      end select

    endif
  enddo
  do i=1, no_nemo_vector_fields
    vfield => extract_vector_field(state, nemo_vector_field_names(i))
    if (have_option(trim(vfield%option_path)//"/prescribed")) then
      call get_option(trim(vfield%option_path) // "/prescribed/value/NEMO_data/format", format)
        select case (format)
          case ("Velocity")
            call remap_field(velocity_t, vfield)
        end select
    endif
  enddo

  call deallocate_temp_fields

end subroutine

subroutine load_nemo_values(state)

  type(state_type), intent(in) :: state

  type(mesh_type) :: input_mesh
  character(len=FIELD_NAME_LEN) input_mesh_name

  real :: current_time
  logical*1 :: on_sphere

  ! Temporary arrays to store the data read from the netCDF
  real, dimension(3) :: temp_vector_3D
  real, dimension(:), allocatable :: X, Y, Z, Temperature, Salinity
  real, dimension(:), allocatable :: U, V, W, SSH, Pressure
  integer :: NNodes, i
  
  call get_option('/ocean_forcing/mesh_choice/mesh/name', input_mesh_name)
  input_mesh = extract_mesh(state, input_mesh_name)
  position = get_coordinate_field(state, input_mesh)

  NNodes=node_count(input_mesh)

  allocate(X(NNodes), Y(NNodes), Z(NNodes), Temperature(NNodes), Salinity(NNodes), &
           U(NNodes), V(NNodes), W(NNodes), SSH(NNodes), Pressure(NNodes))

  do i=1,NNodes
      temp_vector_3D = node_val(position,i)
      X(i) = temp_vector_3D(1) 
      Y(i) = temp_vector_3D(2)
      Z(i) = temp_vector_3D(3)
  end do

  call get_option("/timestepping/current_time",current_time)
  on_sphere=have_option('/geometry/spherical_earth')

  call get_nemo_variables(current_time, X, Y, Z, Temperature, Salinity, U, V, W, &
                          SSH, NNodes)

  call allocate(temperature_t, input_mesh, name="temperature")
  call allocate(salinity_t, input_mesh, name="salinity")
  call allocate(pressure_t, input_mesh, name="pressure")
  call allocate(velocity_t, 3, input_mesh, name="velocity")

  do i=1,NNodes
     call set(temperature_t,i,Temperature(i))
     call set(salinity_t,i,Salinity(i))
     call set(pressure_t,i,9.81*SSH(i))
  enddo

  do i=1,NNodes
     temp_vector_3D(1)=U(i)
     temp_vector_3D(2)=V(i)
     temp_vector_3D(3)=W(i)
     call set(velocity_t,i,temp_vector_3D)
  enddo

  deallocate(X, Y, Z, Temperature, Salinity, &
             U, V, W, SSH, Pressure)

end subroutine

end module nemo_states_module