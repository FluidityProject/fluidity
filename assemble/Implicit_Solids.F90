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
!    C.Pain@Imperial.ac.uk
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

module implicit_solids
! these 5 need to be on top and in this order, 
! so as not to confuse silly old intel compiler 
  use quadrature
  use elements
  use sparse_tools
  use fields
  use state_module
!
  use vtk_interfaces
  use linked_lists
  use intersection_finder_module
  use tetrahedron_intersection_module
  use unify_meshes_module
  use unittest_tools
  use global_parameters, only: FIELD_NAME_LEN, PYTHON_FUNC_LEN, &
       dt, timestep, current_time
  use spud
  use timeloop_utilities
  use fefields, only: compute_lumped_mass
  use parallel_tools
  use interpolation_module
  use diagnostic_variables
  use qmesh_module
  use mesh_files
  use read_triangle
  use fields_manipulation

  implicit none

#ifdef USING_FEMDEM
  interface 
     subroutine y3d_allocate_femdem(string, nodes, elements, edges)
       character(len=*), intent(in) :: string
       integer, intent(out) :: nodes, elements, edges
     end subroutine y3d_allocate_femdem
  end interface

  interface
     subroutine y3d_populate_femdem(ele1, ele2 ,ele3, ele4, &
          face1, face2, face3)
       integer, dimension(*), intent(out) :: ele1, ele2, ele3, ele4
       integer, dimension(*), intent(out) :: face1, face2, face3
     end subroutine y3d_populate_femdem
  end interface

  interface 
     subroutine y3dfemdem(ext_mesh_name, dt, rho, &
          xs, ys, zs, us, vs, ws, uf, vf, wf)
       real, intent(in) :: dt, rho
       character(len=*), intent(in) :: ext_mesh_name
       real, dimension(*), intent(out) :: xs, ys, zs
       real, dimension(*), intent(out) :: us, vs, ws
       real, dimension(*), intent(in) :: uf, vf, wf
     end subroutine y3dfemdem
  end interface
#endif

  type(vector_field), target, save :: ext_pos_solid_vel, ext_pos_fluid_vel
  type(vector_field), target, save :: fl_pos_solid_vel

  type(scalar_field), save :: solid_local
  type(vector_field), save :: external_positions
  real, save :: beta
  real, dimension(:,:), allocatable, save :: translation_coordinates
  integer, save :: number_of_solids
  logical, save :: one_way_coupling, two_way_coupling, multiple_solids

  private
  public:: solids

contains

  subroutine solids(state, its)

    type(state_type), intent(inout) :: state
    type(vector_field), pointer :: positions
    integer, intent(in) :: its
    type(scalar_field), pointer :: solid
    integer, save :: itinoi
    logical, save :: print_drag, do_calculate_volume_fraction
    logical, save :: init=.false.
    integer :: dat_unit, stat
    character(len=PYTHON_FUNC_LEN) :: python_function
    character(len=field_name_len) :: external_mesh_name
    integer :: i, quad_degree
    integer, save :: dim

    ewrite(2, *) "inside implicit_solids"

    solid => extract_scalar_field(state, "SolidConcentration")
    call zero(solid)

    if (.not. init) then

       call get_option("/implicit_solids/beta", beta, default=1.)

       call get_option("/timestepping/nonlinear_iterations", itinoi)
       one_way_coupling = have_option("/implicit_solids/one_way_coupling/")
       two_way_coupling = have_option("/implicit_solids/two_way_coupling/")

       call get_option("/geometry/dimension", dim)

       if (one_way_coupling) then

          ! check for mutiple solids and get translation coordinates
          call get_option( &
               "/implicit_solids/one_way_coupling/number_of_solids", &
               number_of_solids, stat)

          if (stat /= 0) number_of_solids = 1
          multiple_solids = (number_of_solids>1)

          allocate(translation_coordinates(dim, number_of_solids))
          if (have_option("/implicit_solids/one_way_coupling/python")) then
             call get_option(&
                  "/implicit_solids/one_way_coupling/python", &
                  python_function)
             call set_detector_coords_from_python(translation_coordinates, &
                  number_of_solids, python_function, current_time)
          else
             translation_coordinates = 0.
          end if

          ! get external mesh and compare meshes dimensions
          call get_option("/implicit_solids/one_way_coupling/mesh_name", &
               external_mesh_name)

          call get_option("/geometry/quadrature/degree", quad_degree)

          external_positions = &
               read_triangle_serial(trim(external_mesh_name), &
               quad_degree=quad_degree)

          positions => extract_vector_field(state, "Coordinate")

          assert(positions%dim >= 2)
          assert(positions%dim == external_positions%dim)

          ! figure out if we want to print out the drag and initialise drag file
          print_drag = have_option("/implicit_solids/one_way_coupling/print_drag")

          if (GetRank() == 0 .and. print_drag) then
             dat_unit = free_unit()
             open(dat_unit, file="drag_force", status="replace")
             close(dat_unit)
          end if

       else if (two_way_coupling)  then
          call femdem_two_way_initialise(state)
       else
          FLExit("implicit_solids: Don't know what to do...")
       end if

       do_calculate_volume_fraction = .true.

       init=.true.
    end if

    if (one_way_coupling) then

       if (do_calculate_volume_fraction) then

          call allocate(solid_local, solid%mesh, "SolidConcentrationLocal")
          call zero(solid_local)

          do i = 1, number_of_solids
             call calculate_volume_fraction(state, i)
          end do
       end if

       call set(solid, solid_local)
       ewrite_minmax(solid)

       call set_absorption_coefficient(state)

       if (print_drag .and. its==itinoi) call print_drag_force(state)

       do_calculate_volume_fraction = .false.
       if (do_adapt_mesh(current_time, timestep) .and. its==itinoi) then
          call deallocate(solid_local)
          do_calculate_volume_fraction = .true.
       end if

    else if (two_way_coupling) then

       if (do_calculate_volume_fraction) then
          call allocate(solid_local, &
               solid%mesh, "SolidConcentrationLocal")
          call allocate(fl_pos_solid_vel, dim, &
               solid%mesh, "SolidVelocity")
       end if

       if (its == 1) then
          ! update          : 1. external_positions
          !                   2. ext_pos_solid_vel
          ! return to femdem: 1. ext_pos_fluid_vel
          call femdem_two_way_update

          ! interpolate the fluid velocity from
          ! the fluidity mesh to the femdem mesh
          call femdem_interpolation(state, "out")

          ! interpolate the solid velocity from
          ! the femdem mesh to the fluidity mesh
          call zero(fl_pos_solid_vel)
          call femdem_interpolation(state, "in")

          call zero(solid_local)
          call calculate_volume_fraction(state)
       end if

       call set(solid, solid_local)
       ewrite_minmax(solid)

       call set_source(state)
       call set_absorption_coefficient(state)

       do_calculate_volume_fraction = .false.
       if (do_adapt_mesh(current_time, timestep) .and. its==itinoi) then
          call deallocate(solid_local)
          call deallocate(fl_pos_solid_vel)
          do_calculate_volume_fraction = .true.
       end if

    end if

    ewrite(2, *) "leaving implicit_solids"

  end subroutine solids

  !----------------------------------------------------------------------------

  subroutine calculate_volume_fraction(state, solid_number)

    type(state_type),intent(inout) :: state
    integer, intent(in), optional :: solid_number
    type(vector_field), pointer :: positions
    type(vector_field) :: external_positions_local
    integer :: ele_A, ele_B, ele_C
    type(tet_type) :: tet_A, tet_B
    type(plane_type), dimension(:), allocatable :: planes_A
    integer :: stat, nintersections, i, j, k, ntests
    integer, dimension(:), pointer :: ele_A_nodes
    type(vector_field) :: intersection
    real, dimension(:,:), allocatable :: pos_A
    real, dimension(:), allocatable :: detwei
    real :: vol, ele_A_vol

    ewrite(2, *) "inside calculate_volume_fraction"

    positions => extract_vector_field(state, "Coordinate")

    external_positions_local = external_positions

    ! translate coordinates for multiple solids...
    if (one_way_coupling .and. multiple_solids) then
       do i = 1, node_count(external_positions)
          external_positions_local%val(1)%ptr(i) = &
               external_positions%val(1)%ptr(i) + &
               translation_coordinates(1, solid_number)
          external_positions_local%val(2)%ptr(i) = &
               external_positions%val(2)%ptr(i) + &
               translation_coordinates(2, solid_number)
          if (external_positions%dim == 3) then
             external_positions_local%val(3)%ptr(i) = &
                  external_positions%val(3)%ptr(i) + &
                  translation_coordinates(3, solid_number)
          end if
       end do
    end if

    call rtree_intersection_finder_set_input(positions)

    do ele_B = 1, ele_count(external_positions_local)

       call rtree_intersection_finder_find(external_positions_local, ele_B)
       call rtree_intersection_finder_query_output(nintersections)

       if (positions%dim == 3) then
          tet_B%v = ele_val(external_positions_local, ele_B)
       else
          call intersector_set_dimension(positions%dim)
       end if

       do j = 1, nintersections

          call rtree_intersection_finder_get_output(ele_A, j)

          if (positions%dim == 3) then

             if (ele_loc(positions, ele_A)==4) then
                ! tets
                allocate(planes_A(4))
                tet_A%v = ele_val(positions, ele_A)
                planes_A = get_planes(tet_A)
             else
                ! hexes
                allocate(planes_A(6))
                planes_A = get_planes(positions, ele_A)
             end if

             call intersect_tets(tet_B, planes_A, &
                  ele_shape(external_positions_local, ele_B), &
                  stat=stat, output=intersection)

             deallocate(planes_A)

          else

             allocate(pos_A(positions%dim, ele_loc(positions, ele_A)))
             pos_A = ele_val(positions, ele_A)
             intersection = intersect_elements(external_positions_local, &
                  ele_B, pos_A, ele_shape(external_positions_local, ele_B))
             deallocate(pos_A)
             stat = 0

          end if

          if (stat == 1) cycle

          vol = 0.
          do ele_C = 1, ele_count(intersection)
             vol = vol + abs(simplex_volume(intersection, ele_C))
          end do

          allocate(detwei(ele_ngi(positions, ele_A)))
          call transform_to_physical(positions, ele_A, detwei=detwei)
          ele_A_vol = sum(detwei)
          deallocate(detwei)

          ele_A_nodes => ele_nodes(positions, ele_A)
          do k = 1, size(ele_A_nodes)
             call addto(solid_local, ele_A_nodes(k), vol/ele_A_vol)
          end do

          call deallocate(intersection)

       end do

    end do

    do i = 1, node_count(solid_local)
       call set(solid_local, i, max(0., min(1., node_val(solid_local, i))))
    end do

    ewrite_minmax(solid_local)

    call finalise_tet_intersector
    call rtree_intersection_finder_reset(ntests)

    ewrite(2, *) "leaving calculate_volume_fraction"

  end subroutine calculate_volume_fraction

  !----------------------------------------------------------------------------

  subroutine set_absorption_coefficient(state)

    type(state_type),intent(inout) :: state
    type(vector_field), pointer :: absorption
    integer :: i, j
    real :: sigma

    ewrite(2, *) "inside set_absorption_coefficient"

    absorption => extract_vector_field(state, "VelocityAbsorption")
    call zero(absorption)

    do i = 1, node_count(absorption)
       sigma = node_val(solid_local, i)*beta/dt
       do j = 1, absorption%dim
          call set(absorption, j, i, sigma)
       end do
    end do

    ewrite(2, *) "leaving set_absorption_coefficient"

  end subroutine set_absorption_coefficient

  !----------------------------------------------------------------------------

  subroutine set_source(state)

    type(state_type), intent(inout) :: state
    type(vector_field), pointer :: source
    integer :: i, j
    real :: sigma

    ewrite(2, *) "inside set_source"

    source => extract_vector_field(state, "VelocitySource")
    call zero(source)

    do i = 1, node_count(source)
       sigma = node_val(solid_local, i)*beta/dt
       do j = 1, source%dim
          call set(source, j, i, sigma * node_val(fl_pos_solid_vel, j, i))
       end do
    end do

    ewrite(2, *) "leaving set_source"

  end subroutine set_source

  !----------------------------------------------------------------------------

  subroutine print_drag_force(state)

    type(state_type),intent(inout) :: state
    type(vector_field), pointer :: velocity, positions, absorption
    type(scalar_field) :: lumped_mass
    real, dimension(:), allocatable :: drag
    integer :: i, j

    ewrite(2, *) "inside print_drag_force"

    velocity => extract_vector_field(state, "Velocity")
    absorption => extract_vector_field(state, "VelocityAbsorption")
    positions => extract_vector_field(state, "Coordinate")

    call allocate(lumped_mass, positions%mesh, "Lumped mass")
    call compute_lumped_mass(positions, lumped_mass)    

    allocate(drag(positions%dim))
    drag = 0.

    do i = 1, nowned_nodes(positions)
       do j = 1, positions%dim
          drag(j) = drag(j) + &
               node_val(absorption, 1, i) * node_val(velocity, j, i) * &
               node_val(lumped_mass, i)
       end do
    end do

    call allsumv(drag)

    if (GetRank() == 0) then
       open(1453, file="drag_force", position="append")
       write(1453, *) &
            (drag(i), i = 1, positions%dim)
       close(1453)
    end if

    deallocate(drag)
    call deallocate(lumped_mass)

    ewrite(2, *) "leaving print_drag_force"

  end subroutine print_drag_force

  !----------------------------------------------------------------------------

  subroutine femdem_two_way_initialise(state)

    character(len=field_name_len) :: external_mesh_name
    type(vector_field), pointer :: positions
    type(state_type) :: state
    integer :: quad_degree
    type(mesh_type) :: mesh

    integer :: i, j, loc, sloc
    integer :: dim, nodes, elements, edges

    integer, dimension(:), allocatable :: ele1, ele2, ele3, ele4
    integer, dimension(:), allocatable :: face1, face2, face3

    type(quadrature_type) :: quad
    type(element_type) :: shape
    integer, dimension(:), allocatable :: sndglno, boundary_ids
    integer :: boundaries

    ewrite(2, *) "inside femdem_two_way_initialise"

    call get_option("/implicit_solids/two_way_coupling/mesh_name", &
         external_mesh_name)
    call get_option("/geometry/quadrature/degree", quad_degree)

    ! femdem only supports tets
    loc = 4
    sloc = 3

#ifdef USING_FEMDEM
    call y3d_allocate_femdem(trim(external_mesh_name)//char(0), &
         nodes, elements, edges)

    allocate(ele1(elements)); allocate(ele2(elements))
    allocate(ele3(elements)); allocate(ele4(elements))
    allocate(face1(edges)); allocate(face2(edges))
    allocate(face3(edges))

    call y3d_populate_femdem(ele1, ele2, ele3, ele4,&
         face1, face2, face3)
#endif

    positions => extract_vector_field(state, "Coordinate")

    call get_option("/geometry/dimension", dim)
    assert(dim == 3)

    quad = make_quadrature(loc, dim, degree=quad_degree)
    shape = make_element_shape(loc, dim, 1, quad)

    call allocate(mesh, nodes, elements, shape, name="ExternalCoordinateMesh")
    call allocate(external_positions, dim, mesh, name="ExternalCoordinate")

    ! initialise solid mesh coordinates
    do i = 1, nodes
       forall (j=1:dim)
          external_positions%val(j)%ptr(i) = -66.6
       end forall
    end do

    do i = 1, elements
       external_positions%mesh%ndglno((i-1)*loc+1:i*loc) = &
            (/ele1(i)+1, ele2(i)+1, ele3(i)+1, ele4(i)+1/)
    end do

    boundaries = 1
    sloc = loc - 1

    allocate(sndglno(edges*sloc))
    sndglno = 0

    allocate(boundary_ids(edges))
    boundary_ids = 67

    do i = 1, edges
       sndglno((i-1)*sloc+1:i*sloc)= &
            (/face1(i)+1, face2(i)+1, face3(i)+1/)
    end do

    call add_faces(external_positions%mesh, &
         sndgln=sndglno, &
         boundary_ids=boundary_ids)

    deallocate(sndglno)
    deallocate(boundary_ids)

    call deallocate(mesh)
    call deallocate_element(shape)
    call deallocate(quad)

    external_positions%dim=3

    ! this is the solid velocity on the solid mesh
    call allocate(ext_pos_solid_vel, external_positions%dim, &
         external_positions%mesh, name="femdem_solid_velocity")
    call zero(ext_pos_solid_vel)

    ! this is the interpolated fluid velocity
    ! on the solid mesh
    call allocate(ext_pos_fluid_vel, external_positions%dim, &
         external_positions%mesh, name="femdem_fluid_velocity")
    call zero(ext_pos_fluid_vel)

    assert(node_count(external_positions) == node_count(ext_pos_solid_vel))
    assert(node_count(ext_pos_fluid_vel) == node_count(ext_pos_solid_vel))

#ifdef USING_FEMDEM
    deallocate(ele1, ele2, ele3, ele4)
    deallocate(face1, face2, face3)
#endif

  end subroutine femdem_two_way_initialise

  !----------------------------------------------------------------------------

  subroutine femdem_two_way_update

    real :: rho, dt_f
    character(len=field_name_len) :: external_mesh_name

    rho = 1.
    dt_f = dt

    call get_option("/implicit_solids/two_way_coupling/mesh_name", &
         external_mesh_name)

#ifdef USING_FEMDEM
    ! in  :: ext_pos_fluid_vel
    ! out :: updated external_positions and ext_pos_solid_vel
    call y3dfemdem(trim(external_mesh_name)//char(0), dt_f, rho, &
         external_positions%val(1)%ptr, external_positions%val(2)%ptr, &
         external_positions%val(3)%ptr, &
         ext_pos_solid_vel%val(1)%ptr, ext_pos_solid_vel%val(2)%ptr, &
         ext_pos_solid_vel%val(3)%ptr, &
         ext_pos_fluid_vel%val(1)%ptr, ext_pos_fluid_vel%val(2)%ptr, &
         ext_pos_fluid_vel%val(3)%ptr)
#endif

  end subroutine femdem_two_way_update

  !----------------------------------------------------------------------------

  subroutine femdem_interpolation(state, operation)
 
    character(len=*), intent(in) :: operation
    type(state_type) :: state
    type(state_type) :: alg_ext, alg_fl

    type(mesh_type), pointer :: fl_mesh
    type(vector_field) :: fl_positions
    type(vector_field), pointer :: field_ext, field_fl

    ! this subroutine interpolates velocities
    ! between fluidity and femdem meshes
    ! if operation == "in"  : interpolate femdem into fluidity
    ! if operation == "out" : interpolate fluidity into femdem

    ! read in femdem data into alg_ext state
    ! all data is on the femdem mesh
    call insert(alg_ext, external_positions%mesh, "Mesh")
    call insert(alg_ext, external_positions, "Coordinate")

    if (operation == "in") then

       ! this is the solid velocity on the solid mesh
       field_ext => ext_pos_solid_vel
       call insert(alg_ext, field_ext, "SolidVelocity")

    else if (operation == "out") then

       ! this is the fluid velocity on the solid mesh
       ! this will be returned to femdem...
       field_ext => ext_pos_fluid_vel
       call zero(field_ext)
       call insert(alg_ext, field_ext, "Velocity")

    end if

    ! read in fluidity data into alg_new state
    ! all data is on the fluidity mesh
    fl_mesh => extract_mesh(state, "CoordinateMesh")
    fl_positions = extract_vector_field(state, "Coordinate")
    call insert(alg_fl, fl_mesh, "Mesh")
    call insert(alg_fl, fl_positions, "Coordinate")

    if (operation == "in") then

       ! this is the solid velocity on the fluidity mesh
       ! this will be used to set the source term...
       field_fl => fl_pos_solid_vel
       call zero(field_fl)
       call insert(alg_fl, field_fl, "SolidVelocity")

    else if (operation == "out") then

       ! this is the fluid velocity on the fluidity mesh
       field_fl => extract_vector_field(state, "Velocity")
       call insert(alg_fl, field_fl, "Velocity")

    end if

    ! perform interpolations
    if (operation == "in") then

       ! interpolate the solid velocity from
       ! the solid mesh to the fluidity mesh
       call linear_interpolation(alg_ext, alg_fl, different_domains=.true.)

       ewrite_minmax(ext_pos_solid_vel%val(1)%ptr)
       ewrite_minmax(ext_pos_solid_vel%val(2)%ptr)
       ewrite_minmax(ext_pos_solid_vel%val(3)%ptr)

       ewrite_minmax(field_fl%val(1)%ptr)
       ewrite_minmax(field_fl%val(2)%ptr)
       ewrite_minmax(field_fl%val(3)%ptr)

    else if (operation == "out") then

       ! interpolate the fluid velocity from
       ! the fluid mesh to the solid mesh
       call linear_interpolation(alg_fl, alg_ext, different_domains=.true.)

       ewrite_minmax(ext_pos_fluid_vel%val(1)%ptr)
       ewrite_minmax(ext_pos_fluid_vel%val(2)%ptr)
       ewrite_minmax(ext_pos_fluid_vel%val(3)%ptr)

       ewrite_minmax(field_fl%val(1)%ptr)
       ewrite_minmax(field_fl%val(2)%ptr)
       ewrite_minmax(field_fl%val(3)%ptr)

    else
       FLExit("Don't know what to interpolate...")
    end if
    
    call deallocate(alg_ext)
    call deallocate(alg_fl)
    
  end subroutine femdem_interpolation

  !----------------------------------------------------------------------------

end module implicit_solids
