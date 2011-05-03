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
program test_coupler
  use spud
  use fields
  use field_options
  use fields_allocates
  use fields_manipulation
  use state_module
  use FLDebug
  use populate_state_module
  use timeloop_utilities
  use sparsity_patterns_meshes
  use diagnostic_fields_wrapper
  use interpolation_manager
  use interpolation_module
  use python_state
  use vtk_interfaces
  use quadrature
  use shape_functions
  use prism
  use mpi_interfaces
  use parallel_tools
  implicit none

  type(state_type), dimension(:), pointer :: state
  integer :: ierr, stat, tracer_id

  Type(PRISM_time_struct) :: model_time
  Type(PRISM_time_struct),dimension(2) :: model_time_bounds

  !Coupler specific declarations
  integer :: comp_id, localComm

  interface !to stop warnings during compile
     subroutine set_global_debug_level(n)
       integer, intent(in) :: n
     end subroutine set_global_debug_level
  end interface

  call python_init()
  call read_command_line()

  if (have_option("/coupling")) then

     !Initialise PRISM. If coupled, do not mpi_init.
     call prism_init('source', ierr)
     call prism_init_comp(comp_id,'source', ierr)

     !Get communicator from OASIS and set as the communicator for femtools.
     call prism_get_localcomm(comp_id, localComm, ierr)
     call set_communicator(localComm)

     call populate_state(state)
     call gen_def_grid(state, comp_id, tracer_id)

     !Initialise OASIS4 time variable
     model_time = PRISM_jobstart_date
     model_time_bounds = PRISM_jobstart_date
     call PRISM_calc_newdate(model_time_bounds(1), -43200.0, ierr)
     call PRISM_calc_newdate(model_time_bounds(2), +43200.0, ierr)

  else

     call mpi_init(ierr)
     call populate_state(state)

  end if

  !Necessary despite no time dependence, suppresses error/warning
  call set_option("/timestepping/current_time", 0.0, stat=stat)

  call interpolate_to_grid(state)
  call execute_coupling(state, tracer_id, model_time, model_time_bounds)

  call prism_terminate(ierr)

contains

  subroutine gen_def_grid(state, comp_id, tracer_id)
    type(state_type), dimension(:), intent(inout) :: state
    integer, intent(in) :: comp_id
    integer, intent(out) :: tracer_id

    integer :: ierr
    integer :: grid_id, method_id, mask_id
    integer, dimension(2,3) :: local_shape

    !The tracer field on the target mesh
    type(scalar_field) :: T_tgt

    type(vector_field) :: X_tgt_Q1, X_tgt_Q0
    type(mesh_type) :: Mesh_tgt_Q1, Mesh_tgt_Q0

    integer, dimension(3) :: struct_sizes

    real, dimension(:), allocatable :: lon_data, lat_data, vrt_data
    real, dimension(:), allocatable :: points_lon, points_lat, points_vrt
    real, dimension(:,:), allocatable :: corners_lon, corners_lat, corners_vrt
    logical, dimension(:,:,:), allocatable :: mask

    integer :: num_nodes, num_elements, i1, i2, i3, i4, i5
    type(quadrature_type) :: tgt_quad
    type(element_type) :: tgt_shape

    call create_structured_grid(lon_data, lat_data, vrt_data)

    !Put structured grid in the OASIS format

    call data_to_points(lon_data,lat_data,vrt_data,points_lon,points_lat,points_vrt)

    call data_to_corners(lon_data,lat_data,vrt_data,corners_lon,corners_lat,corners_vrt)

    !Define grid to OASIS

    call define_local_shape(local_shape)

    call define_mask(mask)

    call prism_def_grid( grid_id, "GRID", comp_id, local_shape, &
         PRISM_reglonlatvrt, ierr)

    call prism_set_points (method_id, "points", grid_id, local_shape, &
         points_lon, points_lat, points_vrt, .true., ierr)

    call prism_set_corners (grid_id, 8, local_shape, &
         corners_lon, corners_lat, corners_vrt, ierr)

    call prism_set_mask(mask_id, grid_id, local_shape, mask, .true., ierr)

    call prism_def_var (tracer_id, "field_ou", grid_id, method_id, mask_id, (/3,0/), &
         local_shape, PRISM_Double_Precision, ierr)

    call prism_enddef(ierr)

    !-----------------------------------------------------
    !- Here we define the mesh_type for a Q1 mesh
    !- this is then used to generate Q0 mesh.
    tgt_quad = make_quadrature(4,2,degree=1,stat=stat)

    tgt_shape = make_element_shape(4,2,1,tgt_quad, stat=stat)

    call grid_counts(num_nodes, num_elements)

    call allocate(Mesh_tgt_Q1, num_nodes, num_elements, tgt_shape, "StructuredMeshQ1" )

    !READ struct_sizes FROM DIAMOND HERE                                                                                                                                                                                                   
    struct_sizes = (/ 5 , 5 , 1 /)

    do i1 = 1, num_elements
       i2 = i1 + floor(float(i1-1)/(struct_sizes(1)))
       call set_ele_nodes(Mesh_tgt_Q1, i1, (/ i2, i2+1, i2+struct_sizes(1)+1, i2+struct_sizes(1)+2  /) )
    end do
    
    call allocate(X_tgt_Q1, 2, Mesh_tgt_Q1, name="StructuredMeshQ1Coordinate")

    do i3 = 1, num_nodes

       i1 = mod(i3,(struct_sizes(1)+1))
       i4 = 1
       if (i1 == 0) then 
          i1 = struct_sizes(1)
          i4 = 2
       end if

       i2 = floor((float(i3)-1) / (struct_sizes(1)+1) ) + 1
       i5 = 1
       if (i2 == struct_sizes(2)+1) then
          i2 = struct_sizes(2)
          i5 = 2
       end if

       call set(X_tgt_Q1,i3,(/ corners_lon(i1,i4), corners_lat(i2,i5) /) )
    end do

    !- Done Q1 mesh, SructuredMesh with StructuredMeshCoordinate
    !-----------------------------------------------------
    !- Define Q0 mesh from Q1 mesh

    !Shape of Q0 elements
    tgt_shape = make_element_shape(4,2,0,tgt_quad, stat=stat)
    !Q0 mesh derived from Q1.
    Mesh_tgt_Q0 = make_mesh(Mesh_tgt_Q1 , shape=tgt_shape, continuity=-1, name="StructuredMeshQ0")

    !Allocate X_tgt_Q0 to store node coordinates
    call allocate(X_tgt_Q0, 2, Mesh_tgt_Q0, name="StructuredMeshQ0Coordinate")

    !Set nodal coordinates for each node
    !(note that num_nodes = num_elements for Q0 mesh)
    ! The center of each element.
    do i3 = 1, num_elements
       i1 = mod(i3,struct_sizes(1))
       if (i1 == 0) i1 = struct_sizes(1)

       i2 = floor((float(i3)-1) / struct_sizes(1) ) + 1
       call set(X_tgt_Q0,i3,(/ points_lon(i1), points_lat(i2) /) )
    end do

    !- Done Q0 mesh, StructuredMeshQ0 with StructuredMeshQ0Coordinates
    !---------------------------------------------------------

    call allocate(T_tgt, Mesh_tgt_Q0, name="GridTracer")

    call insert(state,X_tgt_Q0,name="Coordinate_tgt")
    call insert(state,T_tgt,name="Tracer_tgt")

  end subroutine gen_def_grid

  subroutine interpolate_to_grid(state)
    type(state_type), dimension(:), intent(inout) :: state

    !The tracer we are coupling
    type(scalar_field), pointer :: T_src

    !The tracer field on the target mesh
    type(scalar_field), pointer :: T_tgt

    !Coordinate fields for source and target
    type(vector_field), pointer :: X_src
    type(vector_field), pointer :: X_tgt_Q0

    T_src=>extract_scalar_field(state,"Tracer")
    X_src=>extract_vector_field(state,"Coordinate")

    T_tgt=>extract_scalar_field(state,"Tracer_tgt")
    X_tgt_Q0=>extract_vector_field(state,"Coordinate_tgt")

    !Perform linear interpolation to structured mesh
    call linear_interpolation(T_src, X_src, T_tgt, X_tgt_Q0)

  end subroutine interpolate_to_grid

  subroutine execute_coupling(state, tracer_id, model_time, model_time_bounds)
    type(state_type), dimension(:), intent(inout) :: state

    integer, intent(in) :: tracer_id
    Type(PRISM_time_struct), intent(in) :: model_time
    Type(PRISM_time_struct), dimension(2), intent(in) :: model_time_bounds

    integer :: ierr, info

    !The tracer field on the target mesh
    type(scalar_field), pointer :: T_tgt

    real, dimension(:,:,:), allocatable :: Tracer_coupler

    T_tgt=>extract_scalar_field(state,"Tracer_tgt")

    call set_coupler_variable_from_field(Tracer_coupler,T_tgt)

    call prism_put(tracer_id, model_time, model_time_bounds, &
         Tracer_coupler, info, ierr)

  end subroutine execute_coupling

  subroutine read_command_line()

    character(len=1024) :: argument
    integer :: status, argn

    call set_global_debug_level(0)

    argn=1
    do 
       call get_command_argument(argn, value=argument, status=status)
       argn=argn+1

       if (status/=0) then
          call usage
          stop
       end if

       exit
    end do

    call load_options(argument)

    return

    call usage
    stop

  end subroutine read_command_line

  subroutine usage
    write (0,*) "usage: test_coupler <options_file>"
    write (0,*) ""
  end subroutine usage

  subroutine create_structured_grid(lon_data,lat_data,vrt_data)
    real, dimension(:), allocatable, intent(out) :: lon_data, lat_data, vrt_data

    integer, dimension(3) :: struct_sizes
    real :: pos_start, pos_end
    real :: disp
    integer :: i2

    !READ struct_sizes FROM DIAMOND HERE
    struct_sizes = (/ 5 , 5 , 1 /)

    allocate(lon_data(struct_sizes(1)))
    allocate(lat_data(struct_sizes(2)))
    allocate(vrt_data(struct_sizes(3)))

    !READ pos_end AND pos_start FOR LONGITUDE FROM DIAMOND HERE
    pos_end = 120.0
    pos_start = 0.0

    disp = (pos_end - pos_start) / struct_sizes(1)
    do i2 = 1, struct_sizes(1)
       lon_data(i2) = pos_start + disp * (i2-1)
    enddo

    !READ pos_end AND pos_start FOR LATITUDE FROM DIAMOND HERE
    pos_end = 90.0
    pos_start = 0.0

    disp = (pos_end - pos_start) / struct_sizes(2)
    do i2 = 1, struct_sizes(2)
       lat_data(i2) = pos_start + disp * (i2-1)
    enddo

    !READ pos_end AND pos_start FOR VERTICAL FROM DIAMOND HERE
    pos_end = 001.0
    pos_start = 0.0

    disp = (pos_end - pos_start) / struct_sizes(3)
    do i2 = 1, struct_sizes(3)
       vrt_data(i2) = pos_start + disp * (i2-1)
    enddo

  end subroutine create_structured_grid

  subroutine data_to_points(lon_data,lat_data,vrt_data,points_lon,points_lat,points_vrt)
    real, dimension(:), intent(in) :: lon_data, lat_data, vrt_data
    real, dimension(:), allocatable, intent(out) :: points_lon, points_lat, points_vrt

    integer, dimension(3) :: struct_sizes
    real :: pos_end, pos_start
    real :: disp

    !READ struct_sizes FROM DIAMOND HERE                                                                                                                                                                                                   
    struct_sizes = (/ 5 , 5 , 1 /)

    allocate(points_lon(struct_sizes(1)))
    allocate(points_lat(struct_sizes(2)))
    allocate(points_vrt(struct_sizes(3)))
    
    !READ pos_end AND pos_start FOR LONGITUDE FROM DIAMOND HERE
    pos_end = 120.0
    pos_start = 0.0

    disp = (pos_end - pos_start) / struct_sizes(1)

    points_lon = lon_data + ( disp / 2.0 )

    !READ pos_end AND pos_start FOR LATITUDE FROM DIAMOND HERE
    pos_end = 090.0
    pos_start = 0.0

    disp = (pos_end - pos_start) / struct_sizes(2)

    points_lat = lat_data + ( disp / 2.0 )
    
    !READ pos_end AND pos_start FOR VERTICAL FROM DIAMOND HERE
    pos_end = 001.0
    pos_start = 0.0
    
    disp = (pos_end - pos_start) / struct_sizes(3)

    points_vrt = vrt_data + ( disp / 2.0 )

  end subroutine data_to_points

subroutine data_to_corners(lon_data,lat_data,vrt_data,corners_lon,corners_lat,corners_vrt)
    real, dimension(:), intent(in) :: lon_data, lat_data, vrt_data
    real, dimension(:,:), allocatable, intent(out) :: corners_lon, corners_lat, corners_vrt

    integer, dimension(3) :: struct_sizes
    real :: pos_end, pos_start
    real :: disp

    !READ struct_sizes FROM DIAMOND HERE                                                                                                                                                                                                   
    struct_sizes = (/ 5 , 5 , 1 /)

    allocate(corners_lon(struct_sizes(1),2))
    allocate(corners_lat(struct_sizes(2),2))
    allocate(corners_vrt(struct_sizes(3),2))
    
    !READ pos_end AND pos_start FOR LONGITUDE FROM DIAMOND HERE
    pos_end = 120.0
    pos_start = 0.0

    disp = (pos_end - pos_start) / struct_sizes(1)

    corners_lon(:,1) = lon_data 
    corners_lon(:,2) = lon_data + disp

    !READ pos_end AND pos_start FOR LATITUDE FROM DIAMOND HERE
    pos_end = 90.0
    pos_start = 0.0

    disp = (pos_end - pos_start) / struct_sizes(2)

    corners_lat(:,1) = lat_data
    corners_lat(:,2) = lat_data + disp

    !READ pos_end AND pos_start FOR VERTICAL FROM DIAMOND HERE
    pos_end = 001.0
    pos_start = 0.0
    
    disp = (pos_end - pos_start) / struct_sizes(3)

    corners_vrt(:,1) = vrt_data
    corners_vrt(:,2) = vrt_data + disp

  end subroutine data_to_corners

  subroutine define_local_shape(local_shape)
    integer, dimension(2,3), intent(inout) :: local_shape
    integer, dimension(3) :: struct_sizes

    !READ struct_sizes FROM DIAMOND HERE                                                                                                                                                                                                   
    struct_sizes = (/ 5 , 5 , 1 /)

    local_shape(1,1) = 1
    local_shape(1,2) = 1
    local_shape(1,3) = 1
    local_shape(2,1) = struct_sizes(1)
    local_shape(2,2) = struct_sizes(2)
    local_shape(2,3) = struct_sizes(3)
  end subroutine define_local_shape

  subroutine define_mask(mask)
    logical, dimension(:,:,:), allocatable, intent(out) :: mask

    integer, dimension(3) :: struct_sizes

    !READ struct_sizes FROM DIAMOND HERE
    struct_sizes = (/ 5 , 5 , 1 /)

    allocate(mask(struct_sizes(1),struct_sizes(2),struct_sizes(3)))

    mask = .true.
  end subroutine define_mask

  subroutine set_coupler_variable_from_field(Tracer_coupler, T_tgt)
    real, dimension(:,:,:), allocatable, intent(out) :: Tracer_coupler
    type(scalar_field), intent(in) :: T_tgt

    integer, dimension(3) :: struct_sizes

    !READ struct_sizes FROM DIAMOND HERE
    struct_sizes = (/ 5 , 5 , 1 /)

    allocate(Tracer_coupler(struct_sizes(1),struct_sizes(2),struct_sizes(3)))

    Tracer_coupler = 2.3

    Tracer_coupler = reshape(T_tgt%val, (/ struct_sizes(1), struct_sizes(2), struct_sizes(3) /) )

  end subroutine set_coupler_variable_from_field

  subroutine grid_counts(nodes,elements)
    integer, intent(out) :: nodes, elements

    integer, dimension(3) :: struct_sizes

    !READ struct_sizes FROM DIAMOND HERE
    struct_sizes = (/ 5, 5 , 1 /)

    nodes = (struct_sizes(1)+1) * (struct_sizes(2)+1)
    elements = struct_sizes(1) * struct_sizes(2)
  end subroutine grid_counts

end program test_coupler
