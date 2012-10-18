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
#define INLINE_MATMUL

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
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN, &
       PYTHON_FUNC_LEN, dt, timestep, current_time
  use spud
  use timeloop_utilities
  use fefields, only: compute_lumped_mass
  use parallel_tools
  use diagnostic_variables
  use qmesh_module
  use mesh_files
  use read_triangle
  use fields_manipulation
  use solvers
  use pickers_inquire
  use transform_elements
  use field_derivatives
  use FLDebug
  use supermesh_construction
  use futils
  use meshdiagnostics
  use sparsity_patterns
  use vector_tools
  use tensors
  use fetools
  use interpolation_module
  use adjacency_lists
  use sparse_matrices_fields
  use bound_field_module
  use halos
  use diagnostic_fields
  use boundary_conditions
  use data_structures
  use edge_length_module
  use detector_tools, only: set_detector_coords_from_python

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
          face1, face2, face3, xs, ys, zs)
       integer, dimension(*), intent(out) :: ele1, ele2, ele3, ele4
       integer, dimension(*), intent(out) :: face1, face2, face3
       real, dimension(*), intent(out) :: xs, ys, zs
     end subroutine y3d_populate_femdem
  end interface

  interface 
     subroutine y3dfemdem(ext_mesh_name, dt, rho_f, ibulk, &
          xs, ys, zs, fxs, fys, fzs, uf, vf, wf, us, vs, ws, af)
       real, intent(in) :: dt, rho_f
       integer, intent(in) :: ibulk
       character(len=*), intent(in) :: ext_mesh_name
       real, dimension(*), intent(in) :: uf, vf, wf, af
       real, dimension(*), intent(out) :: xs, ys, zs
       real, dimension(*), intent(out) :: fxs, fys, fzs
       real, dimension(*), intent(out) :: us, vs, ws

     end subroutine y3dfemdem
  end interface
#endif

  interface interpolation_galerkin_femdem
     module procedure interpolation_galerkin_scalars, &
          interpolation_galerkin_single_state_femdem, &
          interpolation_galerkin_multiple_states_femdem
  end interface

  type(vector_field), target, save :: ext_pos_fluid_vel, ext_pos_solid_vel
  type(vector_field), target, save :: ext_pos_solid_force, fl_pos_solid_force
  type(scalar_field), target, save :: ext_pos_solid

  type(scalar_field), save :: solid_local, old_solid_local, interface_local
  type(vector_field), save :: external_positions

  type(tensor_field), save :: metric, edge_lengths

  real, save :: beta, radius, source_intensity, rho_f
  real, save :: solid_peclet_number, fluid_peclet_number

  real, dimension(:), allocatable, save :: pressure_gradient
  real, dimension(:, :), allocatable, save :: translation_coordinates

  integer, save :: number_of_solids
  logical, save :: one_way_coupling, two_way_coupling, multiple_solids
  logical, save :: have_temperature, have_viscosity, do_print_multiple_solids_diagnostics
  logical, save :: do_print_diagnostics, do_calculate_volume_fraction
  logical, save :: have_fixed_temperature, have_fixed_temperature_source, have_radius
  logical, save :: use_bulk_velocity, use_fluid_velocity, have_pressure_gradient

  type(ilist), dimension(:), allocatable, save :: node_to_particle
  type(ilist), save :: surface_nodes, surface_faces

  private
  public:: solids, implicit_solids_nonlinear_iteration_converged, &
       &   remove_dummy_field, add_mass_source_absorption, &
       &   implicit_solids_register_diagnostic, implicit_solids_update, &
       &   implicit_solids_check_options

contains

  subroutine solids(state, its, itinoi)

    type(state_type), intent(inout) :: state
    integer, intent(in) :: its, itinoi

    type(state_type), dimension(1) :: states
    type(tensor_field), pointer :: viscosity
    type(vector_field), pointer :: x
    type(scalar_field), pointer :: solid, old_solid, interface
    integer :: i, stat

    integer, save :: dim
    logical, save :: init=.false.

    ewrite(2, *) "inside implicit_solids"

    ! SolidConcentration will be called /alpha in the comments from now on
    ! Furthermore, the superscript will denote on which mesh the corresponding
    ! variable will live, and the subscript will distinguish between the phase,
    ! i.e. u_f^s will be the fluid velocity on the solid mesh.
    ! solid actually is /alpha_s^f:
    solid => extract_scalar_field(state, "SolidConcentration")

    if (.not. init) then

       call get_option("/geometry/dimension", dim)
       ! Scaling factor beta
       call get_option("/implicit_solids/beta", beta, default=1.)
       call get_option("/implicit_solids/source_intensity", &
            source_intensity, default=0.)
       call get_option("/material_phase::"//trim(state%name)// &
            "/equation_of_state/fluids/linear/reference_density", rho_f, default = 1.0)

       one_way_coupling = have_option("/implicit_solids/one_way_coupling/")
       two_way_coupling = have_option("/implicit_solids/two_way_coupling/")

       viscosity => extract_tensor_field(state, "Viscosity", stat)
       have_viscosity = stat == 0
       ! Initialisation according to 1-way or 2-way coupling:
       if (one_way_coupling) then
          call one_way_initialise(state)
       else if (two_way_coupling)  then
          call femdem_initialise(state)
       else
          FLAbort("implicit_solids: Don't know what to do...")
       end if
       ! Make sure /alpha_s^f is computed at first timestep:
       do_calculate_volume_fraction = .true.
       ! At this stage, everything is initialised:
       init=.true.
    end if

    ! 1-WAY COUPLING
    if (one_way_coupling) then
       ! Computation of /alpha_s^f and the SolidPhase
       ! at the first timestep and after each adapt:
       if (do_calculate_volume_fraction) then
          call allocate(solid_local, solid%mesh, "SolidConcentration")
          call zero(solid_local)

          ! The SolidPhase represents the surface of the immersed body
          ! on the fluids mesh (0 < /alpha_s^f < 1):
          call allocate(interface_local, solid%mesh, "SolidPhase")
          call zero(interface_local)

          allocate(node_to_particle(node_count(solid)))

          ! Computing /alpha_s^f:
          do i = 1, number_of_solids
             ewrite(2, *) "  calculating volume fraction for solid", i
             call calculate_volume_fraction(state, i)
          end do

          ! Compute the SolidPhase:
          call calculate_solid_fluid_interface(state)

          ! 'x' will be the pointer to the coordinate field of the fluids mesh
          x => extract_vector_field(state, "Coordinate")
          ! Allocating variables for computing the absorption term
          call allocate(metric, x%mesh, "ErrorMetric")
          call zero(metric)
          call allocate(edge_lengths, metric%mesh, "EdgeLengths")
          call zero(edge_lengths)

          ! calculate metric and edge lengths (used for the absorption term later on)
          states = (/state/)
          if (have_viscosity) then
             call qmesh(states, metric) ! metric only needed to get the edge_lengths
             call get_edge_lengths(metric, edge_lengths)
          end if
       end if ! end if do_calculate_volume_fraction

       interface => extract_scalar_field(state, "SolidPhase")
       call zero(interface)
       call set(interface, interface_local)

       ! solid being /alpha^f_s before an adapt,
       ! solid_local being the new /alpha^f_s, after an adapt
       call set(solid, solid_local)
       ewrite_minmax(solid)

       ! Set absorption term /sigma
       call set_absorption_coefficient(state)
       ! Set source term (only if temperature, pressure gradient or 2-way coupling)
       call set_source(state)

       if (have_fixed_temperature_source) &
            call calculate_temperature_diffusivity(state)

       ! Check if the mesh is adapted at end of this timestep
       do_calculate_volume_fraction = .false.
       if (do_adapt_mesh(current_time, timestep) .and. its==itinoi) then
          call deallocate(solid_local)
          call deallocate(interface_local)
          call deallocate(edge_lengths)
          call deallocate(metric)
          do_calculate_volume_fraction = .true.
       end if

    ! 2-WAY COUPLING
    else if (two_way_coupling) then

       old_solid => extract_scalar_field(state, "OldSolidConcentration")

       ! this logical is true after an adapt but the volume
       ! fraction is updated-calculated EVERY time step
       if (do_calculate_volume_fraction) then
          call allocate(solid_local, solid%mesh, "SolidConcentration")
          ! set the local field to the interpolated field after an adapt
          call set(solid_local, solid)

          call allocate(old_solid_local, solid%mesh, "OldSolidConcentration")
          call zero(old_solid_local)

          call allocate(fl_pos_solid_force, dim, solid%mesh, "SolidForce")
          call zero(fl_pos_solid_force)

          x => extract_vector_field(state, "Coordinate")
          call allocate(metric, x%mesh, "ErrorMetric")
          call allocate(edge_lengths, metric%mesh, "EdgeLengths")
       end if

       if (its == 1) then

          ! store previous time step volume fraction
          call set(old_solid_local, solid_local)

          ! interpolate the fluid or bulk velocity from
          ! the fluidity mesh to the femdem mesh
          call zero(ext_pos_fluid_vel)
          call set(ext_pos_solid, 1.)
          call femdem_interpolation(state, "out")

          ! update          : 1. external_positions
          !                   2. ext_pos_solid_force (F_s)
          ! return to femdem: 1. ext_pos_fluid_vel (fluid or bulk velocity)
          !                   2. ext_pos_fluid (\alpha_f)
          call femdem_update(state)
 
          ! interpolate the solid force from
          ! the femdem mesh to the fluidity mesh
          ! and update solid_local
          call zero(fl_pos_solid_force)
          call zero(solid_local)
          call femdem_interpolation(state, "in")

          ! bound solid to [0, 1]
          call bound_concentration()

          ! calculate metric and edge lengths
          call zero(metric)
          call zero(edge_lengths)
          states = (/state/)
          call qmesh(states, metric)
          call get_edge_lengths(metric, edge_lengths)

       end if ! end if (its==1)

       ! previous time step volume fraction
       call set(old_solid, old_solid_local)
       ! at the first time step there is no previous time
       ! step, so assume it's the same as the current
       if (timestep==1) call set(old_solid, solid_local)

       ! current time step volume fraction
       call set(solid, solid_local)

       ! Set absorption term /sigma
       call set_absorption_coefficient(state)
       ! set source
       call set_source(state)

       ! Check if the mesh is adapted at end of this timestep
       do_calculate_volume_fraction = .false.
       if (do_adapt_mesh(current_time, timestep) .and. its==itinoi) then
          call deallocate(solid_local)
          call deallocate(old_solid_local)
          call deallocate(fl_pos_solid_force)
          call deallocate(edge_lengths)
          call deallocate(metric)
          do_calculate_volume_fraction = .true.
       end if

    end if

    ewrite(2, *) "leaving implicit_solids"

  end subroutine solids

  !----------------------------------------------------------------------------

  subroutine calculate_volume_fraction(state, solid_number)

    type(state_type), intent(inout) :: state
    integer, intent(in) :: solid_number

    type(vector_field), pointer :: positions
    type(vector_field) :: external_positions_local
    integer :: ele_A, ele_B, ele_C
    type(tet_type) :: tet_A, tet_B
    type(plane_type), dimension(:), allocatable :: planes_A
    integer :: stat, nintersections, i, j, k, ntests
    integer, dimension(:), pointer :: ele_A_nodes
    type(vector_field) :: intersection
    real, dimension(:, :), allocatable :: pos_A
    real, dimension(:), allocatable :: detwei
    real :: vol, ele_A_vol

    ! positions: coordinates of fluid mesh:
    positions => extract_vector_field(state, "Coordinate")

    ! external_positions: coordinates of solid mesh (lives in the whole module)
    ! external_positions_local: coordinates of solid mesh (lives in this subroutine)
    call allocate(external_positions_local, external_positions%dim, &
         external_positions%mesh, name="LocalCoordinates")
    ! Set the solid coordinates of 'external_positions_local' (lives in subroutine)
    ! to the solid coordinates of 'external_positions'  (lives in module)
    call set(external_positions_local, external_positions)

    ! translate coordinates for multiple solids...
    ! IF multiple solids are derived from one solid,
    ! thus all solids are just a copy of the original one
    if (one_way_coupling .and. multiple_solids) then
       do i = 1, node_count(external_positions)
          external_positions_local%val(1, i) = external_positions%val(1, i) + translation_coordinates(1, solid_number)
          external_positions_local%val(2, i) = external_positions%val(2, i) + translation_coordinates(2, solid_number)
          if (external_positions%dim == 3) then
             external_positions_local%val(3, i) = external_positions%val(3, i) + translation_coordinates(3, solid_number)
          end if
       end do
    end if

    ! Set the input of the RTREE finder as the coordinates of the fluids mesh:
    call rtree_intersection_finder_set_input(positions)

    ! For all elements of the solids mesh:
    do ele_B = 1, ele_count(external_positions_local)

       ! Via RTREE, find intersection of solid element 'ele_B' with the
       ! input mesh (fluid coordinate mesh 'positions'):
       call rtree_intersection_finder_find(external_positions_local, ele_B)
       ! Fetch output, the number of intersections of this solid element:
       call rtree_intersection_finder_query_output(nintersections)


       if (positions%dim == 3) then
          ! Get (solid) element value (coordinate), form tetrahedra:
          tet_B%v = ele_val(external_positions_local, ele_B)
       else
          call intersector_set_dimension(positions%dim)
       end if

       ! 1st inner-loop
       ! For all intersections of solid element ele_B with fluid mesh:
       do j = 1, nintersections

          ! Get the donor (fluid) element which intersects with ele_B
          call rtree_intersection_finder_get_output(ele_A, j)

          ! If 3D
          if (positions%dim == 3) then
             ! Get the global coordinates of fluid element ele_A:
             if (ele_loc(positions, ele_A)==4) then
                ! if fluid element is a tetrahedra
                allocate(planes_A(4))
                tet_A%v = ele_val(positions, ele_A)
                planes_A = get_planes(tet_A)
             else
                ! if fluid element is not a tetrahedra, 
                ! assumed to be a hexahedra:
                allocate(planes_A(6))
                planes_A = get_planes(positions, ele_A)
             end if
             ! Get the coordinates of nodes of the intersection
             ! between of both elements, store in intersection
             call intersect_tets(tet_B, planes_A, &
                  ele_shape(external_positions_local, ele_B), &
                  stat=stat, output=intersection)

             deallocate(planes_A)

          else ! 2D

             allocate(pos_A(positions%dim, ele_loc(positions, ele_A)))
             pos_A = ele_val(positions, ele_A)
             intersection = intersect_elements(external_positions_local, &
                  ele_B, pos_A, ele_shape(external_positions_local, ele_B))
             deallocate(pos_A)
             stat = 0

          end if ! end of dim==3

          ! No intersection, cycle:
          if (stat == 1) cycle

          ! Compute intersection volume:
          vol = 0.
          do ele_C = 1, ele_count(intersection)
             vol = vol + abs(simplex_volume(intersection, ele_C))
          end do

          ! Compute the volume of the fluid element:
          allocate(detwei(ele_ngi(positions, ele_A)))
          call transform_to_physical(positions, ele_A, detwei=detwei)
          ele_A_vol = sum(detwei)
          deallocate(detwei)

          ! Compute the volume fraction:
          ! ele_A_nodes: pointer to global node numbers of
          ! fluid element ele_A of the coordinate mesh
          ele_A_nodes => ele_nodes(positions, ele_A)
          do k = 1, size(ele_A_nodes)
             ! Volume fraction by grandy projection
             call addto(solid_local, ele_A_nodes(k), vol/ele_A_vol)
             call insert_ascending(node_to_particle(ele_A_nodes(k)), solid_number)
          end do

          call deallocate(intersection)

       end do ! ele_A, j, nintersections

    end do ! ele_B, ele_count(external_positions_local)

    ! bound solid to [0, 1]
    call bound_concentration()

    call finalise_tet_intersector
    call rtree_intersection_finder_reset(ntests)

    call deallocate(external_positions_local)

  end subroutine calculate_volume_fraction

  !----------------------------------------------------------------------------

  subroutine bound_concentration()

    integer :: i

    do i = 1, node_count(solid_local)
       call set(solid_local, i, max(0., min(1., node_val(solid_local, i))))
    end do

  end subroutine bound_concentration

  !----------------------------------------------------------------------------

  subroutine set_absorption_coefficient(state)

    type(state_type), intent(inout) :: state

    type(tensor_field), pointer :: viscosity
    type(vector_field), pointer :: absorption
    type(scalar_field), pointer :: Tabsorption
    integer :: i, j
    real :: sigma, sigma_1, sigma_2

    ewrite(2, *) "inside set_absorption_coefficient"

    absorption => extract_vector_field(state, "VelocityAbsorption")
    call zero(absorption)

    if (have_viscosity) then

       viscosity => extract_tensor_field(state, "Viscosity")

       do i = 1, node_count(absorption)

          sigma_1 = node_val(solid_local, i) / dt
          sigma_2 = node_val(solid_local, i) * &
               node_val(viscosity, 1, 1, i) / maxval(node_val(edge_lengths, i))**2
          sigma = max(sigma_1, sigma_2) * beta

          do j = 1, absorption%dim
             call set(absorption, j, i, sigma)
          end do

       end do

    else

       do i = 1, node_count(absorption)

          sigma = node_val(solid_local, i) * beta / dt

          do j = 1, absorption%dim
             call set(absorption, j, i, sigma)
          end do
          
       end do

    end if

    if (have_fixed_temperature) then

       ewrite(3, *) "  set absorption for temperature"

       Tabsorption => extract_scalar_field(state, "TemperatureAbsorption")
       call zero(Tabsorption)
       call set(Tabsorption, absorption, 1)

    end if

    ewrite(2, *) "leaving set_absorption_coefficient"

  end subroutine set_absorption_coefficient

  !----------------------------------------------------------------------------

  subroutine set_source(state)

    type(state_type), intent(inout) :: state

    type(vector_field), pointer :: source, positions
    type(vector_field), pointer :: velocity, absorption
    type(scalar_field), pointer :: Tsource
    integer :: i, particle
    real :: x0, y0, z0, x, y, z, sigma

    ewrite(2, *) "inside set_source"

    if (two_way_coupling) then

       velocity => extract_vector_field(state, "Velocity")
       absorption => extract_vector_field(state, "VelocityAbsorption")

       source => extract_vector_field(state, "VelocitySource")
       call zero(source)

       ! (\rho_f \alpha_s / \Delta t) (\hat{u_f} or u) - F_s/ \Delta t
       do i = 1, node_count(source)
          call set(source, i, &
               node_val(absorption, i) * node_val(velocity, i) &
               - node_val(fl_pos_solid_force, i)/dt)
       end do

    end if

    if (have_pressure_gradient) then

       source => extract_vector_field(state, "VelocitySource")
       call zero(source)
       call set(source, pressure_gradient)
       ! remove the source from the solids
       do i = 1, source%dim
          source%val(i,:) = source%val(i,:) * (1. - solid_local%val)
       end do

    end if

    if (have_temperature) then

       ewrite(3, *) "  set source for temperature"

       Tsource => extract_scalar_field(state, "TemperatureSource")
       call zero(Tsource)

       if (have_fixed_temperature_source) then

          ewrite(3, *) "  source intensity", source_intensity

          if (have_radius) then

             ewrite(3, *) "  using an inner zone", radius

             positions => extract_vector_field(state, "Coordinate")

             do i = 1, node_count(Tsource)
                sigma = node_val(solid_local, i)

                if (associated(node_to_particle(i)%firstnode)) then

                   particle = node_to_particle(i)%firstnode%value

                   x0 = translation_coordinates(1, particle); x = node_val(positions, 1, i)
                   y0 = translation_coordinates(2, particle); y = node_val(positions, 2, i)
                   z0 = translation_coordinates(3, particle); z = node_val(positions, 3, i)

                   if ((x-x0)**2 + (y-y0)**2 + (z-z0)**2 > radius**2) sigma = 0.
                end if

                call set(Tsource, i, sigma * source_intensity)
             end do

          else

             do i = 1, node_count(Tsource)
                sigma = node_val(solid_local, i)
                call set(Tsource, i, sigma * source_intensity)
             end do

          end if

       else if (have_fixed_temperature) then

          ewrite(3, *) "  solid temperature", source_intensity

          do i = 1, node_count(Tsource)
             sigma = node_val(solid_local, i)*beta/dt
             call set(Tsource, i, sigma * source_intensity)
          end do

       end if
    end if

    ewrite(2, *) "leaving set_source"

  end subroutine set_source

  !----------------------------------------------------------------------------

  subroutine calculate_temperature_diffusivity(state)

    type(state_type), intent(inout) :: state

    type(tensor_field), pointer :: diffusivity
    integer :: i, j, k
    real :: l

    ewrite(2, *) "inside calculate_temperature_diffusivity"

    diffusivity => extract_tensor_field(state, "TemperatureDiffusivity")
    call zero(diffusivity)

    do i = 1, node_count(diffusivity)

       l = node_val(solid_local, i)/solid_peclet_number + &
            (1.-node_val(solid_local, i))/fluid_peclet_number

       do j = 1, diffusivity%dim(1)
          do k = 1, diffusivity%dim(2)
             if (j==k) call set (diffusivity, j ,k, i, l)
          end do
       end do

    end do

    ewrite(2, *) "leaving calculate_temperature_diffusivity"

  end subroutine calculate_temperature_diffusivity

  !----------------------------------------------------------------------------

  subroutine one_way_initialise(state)

    type(state_type), intent(in) :: state

    type(vector_field), pointer :: positions
    integer :: stat, quad_degree
    character(len=FIELD_NAME_LEN) :: external_mesh_name
    character(len=PYTHON_FUNC_LEN) :: python_function

    ! pointer to vector field of coordinates of fluids mesh:
    positions => extract_vector_field(state, "Coordinate")
    ! figure out if we want to print out diagnostics and initialise files  
    do_print_diagnostics = have_option("/implicit_solids/one_way_coupling/print_diagnostics")
    do_print_multiple_solids_diagnostics = &
       & have_option("/implicit_solids/one_way_coupling/multiple_solids/print_diagnostics")
    ! check for mutiple solids and get translation coordinates
    number_of_solids = 1
    multiple_solids = have_option("/implicit_solids/one_way_coupling/multiple_solids")
    if (multiple_solids) then
       call get_option("/implicit_solids/one_way_coupling/multiple_solids/number_of_solids", &
         & number_of_solids)
    end if

    ! In case of multiple immersed bodies, translate their coordinates:
    allocate(translation_coordinates(positions%dim, number_of_solids))
    if (multiple_solids) then
       call get_option(&
            "/implicit_solids/one_way_coupling/multiple_solids/python", &
            python_function)
       call set_detector_coords_from_python(translation_coordinates, &
            number_of_solids, python_function, current_time)
    else
       translation_coordinates = 0.
    end if

    ! get external mesh and compare meshes dimensions
    call get_option("/implicit_solids/one_way_coupling/mesh/file_name", &
         external_mesh_name)

    call get_option("/geometry/quadrature/degree", quad_degree)

    ! Read in the serial mesh of solid body
    ! external_positions lives in the whole module:
    external_positions = read_triangle_serial(trim(external_mesh_name), &
         quad_degree=quad_degree)

    ! 1D set-ups not supported:
    assert(positions%dim >= 2)
    ! Make sure dimensions of coordinate meshes of the
    ! fluids mesh and solids mesh are equal:
    assert(positions%dim == external_positions%dim)

    ! check temperature related options
    if (have_temperature .and. &
         .not. have_option("/implicit_solids/source")) then
       ewrite(-1, *) "WARNING: Implicit solids are not emitting!"
    end if
    if (.not. have_temperature .and. &
         have_option("/implicit_solids/source")) then
       FLExit("You need to use a Temperature field if you want to &
            have emitting solids.")
    end if

    have_fixed_temperature_source = have_option("/implicit_solids/source/temperature_source")
    have_fixed_temperature = have_option("/implicit_solids/source/temperature")

    if (have_fixed_temperature_source) then
       call get_option("/implicit_solids/source/temperature_source/solid_peclet_number", &
            solid_peclet_number)
       call get_option("/implicit_solids/source/temperature_source/fluid_peclet_number", &
            fluid_peclet_number)
       call get_option("/implicit_solids/source/temperature_source/source_intensity", &
            source_intensity)
       call get_option("/implicit_solids/source/temperature_source/zone_radius", radius, stat)
       have_radius = stat == 0
    end if

    if (have_fixed_temperature) then
       call get_option("/implicit_solids/source/temperature/temperature", source_intensity)
    end if

    have_pressure_gradient = have_option("/implicit_solids/pressure_gradient")
    if (have_pressure_gradient) then
       allocate(pressure_gradient(positions%dim))
       call get_option("/implicit_solids/pressure_gradient/", pressure_gradient)
    end if

  end subroutine one_way_initialise

  !----------------------------------------------------------------------------

  subroutine femdem_initialise(state)

    type(state_type), intent(in) :: state

    character(len=FIELD_NAME_LEN) :: external_mesh_name
    type(vector_field), pointer :: positions
    integer :: quad_degree
    type(mesh_type) :: mesh
    integer :: i, loc, sloc
    integer :: dim, nodes, elements, edges
    integer, dimension(:), allocatable :: ele1, ele2, ele3, ele4
    integer, dimension(:), allocatable :: face1, face2, face3
    real, dimension(:), allocatable :: xs, ys, zs
    type(quadrature_type) :: quad
    type(element_type) :: shape
    integer, dimension(:), allocatable :: sndglno, boundary_ids
    integer :: boundaries

    ewrite(2, *) "inside femdem_initialise"

    call get_option("/implicit_solids/two_way_coupling/mesh/file_name", &
         external_mesh_name)
    call get_option("/geometry/quadrature/degree", quad_degree)
    use_bulk_velocity = &
         have_option("/implicit_solids/two_way_coupling/fluids_scheme/use_bulk_velocity")
    use_fluid_velocity = &
         have_option("/implicit_solids/two_way_coupling/fluids_scheme/use_fluid_velocity")

    ! femdem only supports tets
    loc = 4
    sloc= 3

#ifdef USING_FEMDEM
    call y3d_allocate_femdem(trim(external_mesh_name)//char(0), &
         nodes, elements, edges)

    allocate(ele1(elements)); allocate(ele2(elements))
    allocate(ele3(elements)); allocate(ele4(elements))
    allocate(face1(edges)); allocate(face2(edges))
    allocate(face3(edges))

    allocate(xs(nodes)); allocate(ys(nodes))
    allocate(zs(nodes))

    call y3d_populate_femdem(ele1, ele2, ele3, ele4,&
         face1, face2, face3, xs, ys, zs)
#endif

    positions => extract_vector_field(state, "Coordinate")

    call get_option("/geometry/dimension", dim)
    assert(dim == 3)

    quad = make_quadrature(loc, dim, degree=quad_degree)
    shape = make_element_shape(loc, dim, 1, quad)

    call allocate(mesh, nodes, elements, shape, name="CoordinateMesh")
    call allocate(external_positions, dim, mesh, name="Coordinate")

    ! initialise solid mesh coordinates
    external_positions%val(1,:) = xs
    external_positions%val(2,:) = ys
    external_positions%val(3,:) = zs

    do i = 1, elements
       external_positions%mesh%ndglno((i-1)*loc+1:i*loc) = &
            (/ele1(i)+1, ele2(i)+1, ele3(i)+1, ele4(i)+1/)
    end do

    boundaries = 1
    sloc = loc - 1

    allocate(sndglno(edges*sloc)); sndglno = 0
    allocate(boundary_ids(edges)); boundary_ids = 66

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

    ! this is the solid force
    ! on the solid mesh
    call allocate(ext_pos_solid_force, external_positions%dim, &
         external_positions%mesh, name="SolidForce")
    call zero(ext_pos_solid_force)

    ! this is the interpolated fluid or bulk velocity
    ! on the solid mesh
    call allocate(ext_pos_fluid_vel, external_positions%dim, &
         external_positions%mesh, name="Velocity")
    call zero(ext_pos_fluid_vel)

    ! this is the interpolated solid velocity
    ! on the solid mesh
    call allocate(ext_pos_solid_vel, external_positions%dim, &
         external_positions%mesh, name="Velocity")
    call zero(ext_pos_solid_vel)

    ! this is the interpolated fluid volume fraction
    ! on the solid mesh
    call allocate(ext_pos_solid, &
         external_positions%mesh, name="SolidConcentration")
    call zero(ext_pos_solid)

    assert(node_count(external_positions) == node_count(ext_pos_solid_force))
    assert(node_count(ext_pos_fluid_vel) == node_count(ext_pos_solid_force))

#ifdef USING_FEMDEM
    deallocate(ele1, ele2, ele3, ele4)
    deallocate(face1, face2, face3)
    deallocate(xs, ys, zs)
#endif

  end subroutine femdem_initialise

  !----------------------------------------------------------------------------

  subroutine remove_dummy_field(state)

    type(state_type), intent(inout) :: state

    type(scalar_field), pointer :: dummy
    logical :: have_dummy
    integer :: stat

    ! remove dummy field used for
    ! adapt_at_first_timestep
    dummy => extract_scalar_field(state, "FirstAdaptDummy", stat)
    have_dummy = stat == 0

    if (have_dummy) then
       call remove_scalar_field(state, "FirstAdaptDummy")
       call remove_scalar_field(state, "FirstAdaptDummyInterpolationErrorBound")

       call remove_scalar_field(state, "OldFirstAdaptDummy")
       call remove_scalar_field(state, "OldFirstAdaptDummyInterpolationErrorBound")

       call remove_scalar_field(state, "IteratedFirstAdaptDummy")
       call remove_scalar_field(state, "IteratedFirstAdaptDummyInterpolationErrorBound")

       call delete_option('/material_phase[0]/scalar_field::FirstAdaptDummy')
    else
       ewrite(-1, *) "You should be using a FirstAdaptDummy field &
            & with adapt_at_first_timestep."
    end if

  end subroutine remove_dummy_field

  !----------------------------------------------------------------------------

  subroutine femdem_update(state)

    type(state_type), intent(in) :: state

    character(len=FIELD_NAME_LEN), save :: external_mesh_name
    integer, save :: use_bulk
    logical, save :: init=.false.
    !integer :: i

    ewrite(2, *) "inside femdem_update"

    if (.not. init) then
       call get_option("/implicit_solids/two_way_coupling/mesh/file_name", &
            external_mesh_name)

       if (use_bulk_velocity) then
          use_bulk = 1
       elseif (use_fluid_velocity) then
          use_bulk = 0
       else
          FLAbort("Don't recognise fluids scheme...")
       end if

       init=.true.
    end if

    if (use_fluid_velocity) then
       ewrite(3, *) "calculating the bulk velocity"
       ! ext_pos_solid_vel is u_s, not \hat{u}_s, so multiply
       ! by \alpha_s and then add to ext_pos_fluid_vel to 
       ! form the bulk velocity (\hat{u}_f+\hat{u}_s)
       call scale(ext_pos_solid_vel, ext_pos_solid)
       ! this removes fluid velocity from the volume of the solid
       !do i = 1, ext_pos_fluid_vel%dim
       !   ext_pos_fluid_vel%val(i,:)=ext_pos_fluid_vel%val(i,:) * (1. - ext_pos_solid%val)
       !end do
       call addto(ext_pos_fluid_vel, ext_pos_solid_vel)
    end if

    call zero(external_positions)
    call zero(ext_pos_solid_vel)

#ifdef USING_FEMDEM

    ewrite(2, *) "about to call femdem"

    ! in  :: ext_pos_fluid_vel, ext_pos_fluid = 1. - ext_pos_solid 
    ! out :: updated external_positions, ext_pos_solid_force
    call y3dfemdem(trim(external_mesh_name)//char(0), dt, rho_f, use_bulk, &
         external_positions%val(1,:), external_positions%val(2,:), &
         external_positions%val(3,:), &
         ext_pos_solid_force%val(1,:), ext_pos_solid_force%val(2,:), &
         ext_pos_solid_force%val(3,:), &
         ext_pos_fluid_vel%val(1,:), ext_pos_fluid_vel%val(2,:), &
         ext_pos_fluid_vel%val(3,:), &
         ext_pos_solid_vel%val(1,:), ext_pos_solid_vel%val(2,:), &
         ext_pos_solid_vel%val(3,:), &
         1.-ext_pos_solid%val(:))

    ewrite(2, *) "leaving femdem"

#endif

    ewrite(2, *) "leaving femdem_update"

  end subroutine femdem_update

  !----------------------------------------------------------------------------

  subroutine femdem_interpolation(state, operation)
 
    character(len=*), intent(in) :: operation
    type(state_type), intent(in) :: state

    type(state_type) :: alg_ext_v, alg_fl_v, alg_ext_s, alg_fl_s
    type(mesh_type), pointer :: fl_mesh
    type(vector_field) :: fl_positions
    type(vector_field), pointer :: field_ext_v, field_fl_v
    type(scalar_field), pointer :: field_ext_s, field_fl_s
    character(len=OPTION_PATH_LEN) :: &
         path = "/tmp/galerkin_projection/continuous"
    !integer :: stat

    ! this subroutine interpolates forces and velocities
    ! between fluid and solid mesh
    ! if operation == "in"  : interpolate from femdem to fluidity
    ! if operation == "out" : interpolate from fluidity to femdem

    ! read in femdem data into alg_ext state
    ! all data is on the femdem mesh
    call insert(alg_ext_v, external_positions%mesh, "Mesh")
    call insert(alg_ext_v, external_positions, "Coordinate")
    call insert(alg_ext_s, external_positions%mesh, "Mesh")
    call insert(alg_ext_s, external_positions, "Coordinate")

    ! read in fluidity data into alg_new state
    ! all data is on the fluidity mesh
    fl_mesh => extract_mesh(state, "CoordinateMesh")
    fl_positions = extract_vector_field(state, "Coordinate")
    call insert(alg_fl_v, fl_mesh, "Mesh")
    call insert(alg_fl_v, fl_positions, "Coordinate")
    call insert(alg_fl_s, fl_mesh, "Mesh")
    call insert(alg_fl_s, fl_positions, "Coordinate")

    !call set_solver_options(path, &
    !        ksptype = "cg", &
    !        pctype = "hypre", &
    !        rtol = 1.e-10, &
    !        atol = 0., &
    !        max_its = 10000)
    !
    !call add_option( &
    !     trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", stat)
    !call set_option( &
    !     trim(path)//"/solver/preconditioner[0]/hypre_type[0]/name", "boomeramg")

    call set_solver_options(path, &
            ksptype = "cg", &
            pctype = "mg", &
            rtol = 1.e-10, &
            atol = 0., &
            max_its = 10000)

    path = "/tmp"

    if (operation == "in") then

       fl_pos_solid_force%option_path = path

       ! this is the solid velocity on the fluidity mesh
       ! this will be used to set the source term...
       field_fl_v => fl_pos_solid_force
       call zero(field_fl_v)
       call insert(alg_fl_v, field_fl_v, "SolidForce")

       ! this is the solid velocity on the solid mesh
       field_ext_v => ext_pos_solid_force
       call insert(alg_ext_v, field_ext_v, "SolidForce")

    else if (operation == "out") then

       ext_pos_fluid_vel%option_path = path

       ! this is the fluid or bulk velocity on the solid mesh
       ! this will be returned to femdem...
       field_ext_v => ext_pos_fluid_vel
       call zero(field_ext_v)
       call insert(alg_ext_v, field_ext_v, "Velocity")

       ! this is the fluid or bulk velocity on the fluidity mesh
       field_fl_v => extract_vector_field(state, "Velocity")
       call insert(alg_fl_v, field_fl_v, "Velocity")

       ! this is the solid concentration on the solid mesh
       ! this will be returned to femdem...
       field_ext_s => ext_pos_solid
       if (use_fluid_velocity) call zero(field_ext_s)
       call insert(alg_ext_s, field_ext_s, "SolidConcentration")

       ! this is the solid concentration on the fluidity mesh
       field_fl_s => extract_scalar_field(state, "SolidConcentration")
       call insert(alg_fl_s, field_fl_s, "SolidConcentration")

    end if

    ! perform interpolations
    if (operation == "in") then
  
       ! interpolate the solid force from
       ! the solid mesh to the fluidity mesh
       call interpolation_galerkin_femdem(alg_ext_v, alg_fl_v, field=solid_local)

       ewrite_minmax(ext_pos_solid_force)
       ewrite_minmax(field_fl_v)

    else if (operation == "out") then

       ! interpolate the fluid or bulk velocity and the solid concentration from
       ! the fluid mesh to the solid mesh
       call interpolation_galerkin_femdem(alg_fl_v, alg_ext_v, femdem_out=.true.)

       if (use_fluid_velocity) &
          call linear_interpolation(alg_fl_s, alg_ext_s, different_domains=.true.)

       ewrite_minmax(field_fl_v)
       ewrite_minmax(ext_pos_fluid_vel)

       ewrite_minmax(field_fl_s)
       ewrite_minmax(field_ext_s)

    else
       FLAbort("Don't know what to interpolate...")
    end if

    call deallocate(alg_ext_v); call deallocate(alg_ext_s)
    call deallocate(alg_fl_v); call deallocate(alg_fl_s)

  end subroutine femdem_interpolation

  !----------------------------------------------------------------------------

  subroutine implicit_solids_nonlinear_iteration_converged()

    if (do_adapt_mesh(current_time, timestep)) then
       if (one_way_coupling) then
          call deallocate(solid_local)
          call deallocate(interface_local)
          call deallocate(edge_lengths)
          call deallocate(metric)
       else if (two_way_coupling) then
          call deallocate(solid_local)
          call deallocate(old_solid_local)
       end if
        do_calculate_volume_fraction = .true.
    end if

  end subroutine implicit_solids_nonlinear_iteration_converged

  !----------------------------------------------------------------------------

  subroutine implicit_solids_force_computation(state, force, particle_force)

    type(state_type), intent(in) :: state
    real, dimension(:), allocatable, intent(out) :: force
    real, dimension(:, :), allocatable, intent(out) :: particle_force

    type(vector_field), pointer :: velocity, positions, absorption
    type(scalar_field) :: lumped_mass, lumped_mass_velocity_mesh
    integer :: i, j, particle
    type(inode), pointer :: node1

    ewrite(2, *) "inside implicit_solids_force_computation"

    ! Update the computation of the force on the solids

    velocity => extract_vector_field(state, "Velocity")
    absorption => extract_vector_field(state, "VelocityAbsorption")

    positions => extract_vector_field(state, "Coordinate")

    call allocate(lumped_mass, positions%mesh, "LumpedMass")
    call compute_lumped_mass(positions, lumped_mass)    

    call allocate(lumped_mass_velocity_mesh, velocity%mesh, "LumpedMassVelocityMesh")
    call remap_field(lumped_mass, lumped_mass_velocity_mesh)

    allocate(particle_force(number_of_solids, velocity%dim)); particle_force = 0.
    allocate(force(velocity%dim)); force = 0.

    do i = 1, nowned_nodes(velocity)
       node1 => node_to_particle(i)%firstnode
       do while (associated(node1))
          particle = node1%value
          do j = 1, velocity%dim
             particle_force(particle, j) = particle_force(particle, j) + &
                  node_val(absorption, j, i) * node_val(velocity, j, i) * &
                  node_val(lumped_mass_velocity_mesh, i)
          end do
          node1 => node1%next  
       end do
    end do
    
    do i = 1, number_of_solids
       call allsum(particle_force(i, :))
    end do
    do i = 1, velocity%dim
       do j = 1, number_of_solids
          force(i) = force(i) + particle_force(j, i)
       end do
    end do

    call deallocate(lumped_mass)
    call deallocate(lumped_mass_velocity_mesh)

    ewrite(2, *) "leaving implicit_solids_force_computation"

  end subroutine implicit_solids_force_computation

  !----------------------------------------------------------------------------

  subroutine implicit_solids_temperature_computation(state, T_w_avg, q_avg, wall_temperature, q)

    type(state_type), intent(in) :: state
    real, dimension(:), allocatable, intent(out) :: wall_temperature, q
    real, intent(out) :: T_w_avg, q_avg

    type(vector_field), pointer :: positions
    type(vector_field) :: temperature_grad
    type(scalar_field), pointer :: temperature
    type(tensor_field), pointer :: temperature_conductivity
    type(scalar_field) :: lumped_mass    
    real, dimension(:), allocatable :: solid_mass
    integer, dimension(:), allocatable :: face_nodes
    integer :: i, gi, particle
    type(inode), pointer :: node1, node2
    real, dimension(:, :), allocatable :: temperature_grad_at_quad, normal
    real, dimension(:), allocatable :: detwei, T_grad_dot_n_at_quad
    type(element_type), pointer :: T_f_shape  
    real :: k_f

    ewrite(2, *) "inside implicit_solids_temperature_computation"

    positions => extract_vector_field(state, "Coordinate")

    call allocate(lumped_mass, positions%mesh, "LumpedMass")
    call compute_lumped_mass(positions, lumped_mass)    

    allocate(face_nodes(face_loc(positions, 1))); face_nodes = 0
    allocate(wall_temperature(number_of_solids)); wall_temperature = 0.
    allocate(solid_mass(number_of_solids)); solid_mass = 0.
    allocate(q(number_of_solids)); q = 0.

    temperature => extract_scalar_field(state, "Temperature")

    ! calculate the solid-fluid interface temperature
    node1 => surface_nodes%firstnode
    do while (associated(node1))

       i = node1%value ! node number
       node2 => node_to_particle(i)%firstnode

       do while (associated(node2))

          particle = node2%value

          wall_temperature(particle) = wall_temperature(particle) + &
               node_val(temperature, i) * node_val(lumped_mass, i)
          solid_mass(particle) = solid_mass(particle) + node_val(lumped_mass, i)

          node2 => node2%next
       end do

       node1 => node1%next
    end do

    call allsum(wall_temperature)
    call allsum(solid_mass)

    wall_temperature = wall_temperature / solid_mass

    ! calculate the solid-fluid interface heat transfer rate
    call allocate(temperature_grad, positions%dim, &
         temperature%mesh, "TemperatureGradient")
    call zero(temperature_grad)
    call grad(temperature, positions, temperature_grad)

    temperature_conductivity => extract_tensor_field(state, "TemperatureDiffusivity")
    if (have_fixed_temperature_source) then
       k_f = minval(temperature_conductivity%val(1, 1, :))
    else
       k_f = node_val(temperature_conductivity, 1, 1, 1)
    end if

    allocate(temperature_grad_at_quad(positions%dim, face_ngi(temperature, 1)))
    allocate(detwei(face_ngi(temperature, 1)))
    allocate(normal(positions%dim, face_ngi(temperature, 1)))
    allocate(T_grad_dot_n_at_quad(face_ngi(temperature, 1)))
    node1 => surface_faces%firstnode
    do while (associated(node1))

       i = node1%value ! face number
       face_nodes = face_global_nodes(positions, i)
       particle = node_to_particle(face_nodes(1))%firstnode%value
       temperature_grad_at_quad = face_val_at_quad(temperature_grad, i)
       call transform_facet_to_physical(positions, i, &
            detwei_f=detwei, normal=normal)
       T_f_shape => face_shape(temperature, i)
       T_grad_dot_n_at_quad = 0.
       do gi = 1, face_ngi(temperature, i)
          T_grad_dot_n_at_quad(gi) = &
               dot_product(temperature_grad_at_quad(:, gi), normal(:, gi))
       end do

       q(particle) = q(particle) + &
            sum(shape_rhs(T_f_shape, detwei*T_grad_dot_n_at_quad*k_f))
       node1 => node1%next
    end do

    call allsum(q)

    T_w_avg = 0.; q_avg = 0.
    do i = 1, number_of_solids
       T_w_avg = T_w_avg + wall_temperature(i)
       q_avg = q_avg + q(i) 
    end do
    T_w_avg = T_w_avg / number_of_solids
    q_avg = q_avg / number_of_solids

    deallocate(temperature_grad_at_quad, detwei, normal, T_grad_dot_n_at_quad)
    call deallocate(temperature_grad)

    deallocate(face_nodes, solid_mass)
    call deallocate(lumped_mass)

    ewrite(2, *) "leaving implicit_solids_temperature_computation"

  end subroutine implicit_solids_temperature_computation

  !----------------------------------------------------------------------------

  subroutine implicit_solids_update(state)

    type(state_type), intent(in) :: state
    type(vector_field), pointer :: positions
    type(scalar_field) :: lumped_mass

    real, dimension(:), allocatable :: force
    real, dimension(:, :), allocatable :: particle_force
    real, dimension(:), allocatable :: wall_temperature, q
    real :: T_w_avg, q_avg
    integer :: i, j, str_size
    character(len=254) :: fmt, buffer

    ewrite(2, *) "inside implicit_solids_update"
    
    ! Update the computation of the diagnostics
    ! Only one-way coupling for now
    if (one_way_coupling .and. do_print_diagnostics) then

       call implicit_solids_force_computation(state, force, particle_force)

       str_size=len_trim(int2str(number_of_solids))
       fmt="(I"//int2str(str_size)//"."//int2str(str_size)//")"

       ! Register the force on a solid body
       call set_diagnostic(name="Force", statistic="Value", value=(/ force /))
       
       if (do_print_multiple_solids_diagnostics) then
          do j = 1, number_of_solids
             write(buffer, fmt) j
             call set_diagnostic(name="ForceOnSolid"//buffer, statistic="Value", value=(/ particle_force(j, i) /))
          end do
       end if

       if (have_temperature) then
          call implicit_solids_temperature_computation(state, T_w_avg, q_avg, wall_temperature, q)

          ! Register the diagnostics
          call set_diagnostic(name="WallTemperature", statistic="Value", value=(/ T_w_avg /))
          call set_diagnostic(name="HeatTransfer", statistic="Value", value=(/ q_avg /))
          
          if (multiple_solids .and. do_print_multiple_solids_diagnostics) then
             do i = 1, number_of_solids
                write(buffer, fmt) i
                call set_diagnostic(name="WallTemperatureOnSolid"//buffer, statistic="Value", value=(/ wall_temperature(i) /))
                call set_diagnostic(name="HeatTransferAtSolid"//buffer, statistic="Value", value=(/ q(i) /))
             end do
          end if

          deallocate(wall_temperature, q)
       end if

       deallocate(force, particle_force)

       if (do_adapt_mesh(current_time, timestep)) then
          call deallocate(node_to_particle)
          deallocate(node_to_particle)
          call deallocate(surface_nodes)
          call deallocate(surface_faces)
       end if
    end if

    ewrite(2, *) "leaving implicit_solids_update"

  end subroutine implicit_solids_update

  !----------------------------------------------------------------------------

  subroutine calculate_solid_fluid_interface(state)

    type(state_type), intent(in) :: state

    type(vector_field), pointer :: x
    integer, dimension(:), allocatable :: face_nodes
    integer, dimension(:), pointer :: faces, nodes
    integer :: ele, face, nface, node, nnode, cnt

    ewrite(2, *) "  mapping the solid-fluid interface"

    x => extract_vector_field(state, "Coordinate")

    do ele = 1, ele_count(x)

       ! element with one face on the solid
       if (all(ele_val(solid_local, ele) == 0.)) cycle
       if (all(ele_val(solid_local, ele) == 1.)) cycle
       if (all(ele_val(solid_local, ele)  > 0.)) cycle
       if (sum(ele_val(solid_local, ele)) > face_loc(x, 1)) cycle
       nodes => ele_nodes(x, ele)
       cnt=0
       do nnode = 1, size(nodes)
          node = nodes(nnode)
          if (node_val(solid_local, node) > 0.) cnt=cnt+1
       end do
       if (cnt > face_loc(x, 1)) cycle

       faces => ele_faces(x, ele)

       do nface = 1, size(faces)
          face = faces(nface)

          ! this is the solid face 
          if (all(face_val(solid_local, face) /= 0.)) then

             allocate(face_nodes(face_loc(x, face)))
             face_nodes = face_global_nodes(x, face)

             ! interface node list
             do node = 1, size(face_nodes)
                if (face_nodes(node) < nowned_nodes(x)) &
                     call insert_ascending(surface_nodes, face_nodes(node))
             end do

             ! this is for printing the interface to the vtu
             call set(interface_local, face_nodes, spread(1., 1, size(face_nodes)))

             ! interface face list
             if (all(face_nodes < nowned_nodes(x))) &
                  call insert_ascending(surface_faces, face)

             deallocate(face_nodes)
             exit

          end if

       end do ! face
    end do ! ele

  end subroutine calculate_solid_fluid_interface

  !----------------------------------------------------------------------------

  subroutine add_mass_source_absorption(ct_rhs, state)

    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: ct_rhs

    type(scalar_field), pointer :: solid, old_solid, p
    type(vector_field), pointer :: x
    type(element_type), pointer :: p_shape 
    type(element_type) :: test_function
    integer :: ele

    ! if solving for the bulk
    ! velocity nothing to be done
    if (use_bulk_velocity) return

    ewrite(2, *) "inside add_mass_source_absorption"

    solid => extract_scalar_field(state, "SolidConcentration")
    old_solid => extract_scalar_field(state, "OldSolidConcentration")
    x => extract_vector_field(state, "Coordinate")
    p => extract_scalar_field(state, "Pressure")

    do ele = 1, element_count(p)
       p_shape => ele_shape(p, ele)
       test_function = p_shape

       call add_ct_rhs_element_cg(ele, test_function, &
            p_shape, x, p, solid, old_solid, ct_rhs)
    end do

    ewrite(2, *) "leaving add_mass_source_absorption"

  contains

    subroutine add_ct_rhs_element_cg(ele, test_function, &
         p_shape, x, p, solid, old_solid, ct_rhs)

      type(scalar_field), intent(inout) :: ct_rhs
      integer, intent(in) :: ele
      type(element_type), intent(in) :: test_function, p_shape
      type(vector_field), intent(in) :: x
      type(scalar_field), intent(in) :: solid, old_solid, p

      real, dimension(ele_ngi(p, ele)) :: detwei
      real, dimension(ele_loc(p, ele), ele_ngi(p, ele), x%dim) :: dp_t
      integer, dimension(:), pointer :: p_ele
      real, dimension(ele_loc(x, ele)) :: ds, mat

      call transform_to_physical(x, ele, &
           p_shape, dshape=dp_t, detwei=detwei)

      mat = shape_rhs(test_function, detwei)
      ds = ele_val(solid, ele) - ele_val(old_solid, ele)
      p_ele => ele_nodes(p, ele)

      call addto(ct_rhs, p_ele, mat*ds/dt)

    end subroutine add_ct_rhs_element_cg

  end subroutine add_mass_source_absorption

  !----------------------------------------------------------------------------

  subroutine interpolation_galerkin_single_state_femdem(old_state, &
       new_state, field, femdem_out)

    type(state_type), intent(inout) :: old_state, new_state
    type(scalar_field), intent(inout), optional :: field
    logical, intent(in), optional :: femdem_out

    type(state_type), dimension(1) :: old_states, new_states
    
    old_states = (/old_state/)
    new_states = (/new_state/)
    call interpolation_galerkin_femdem(old_states, &
         new_states, field=field, femdem_out=femdem_out)
    old_state = old_states(1)
    new_state = new_states(1)
    
  end subroutine interpolation_galerkin_single_state_femdem

  subroutine interpolation_galerkin_multiple_states_femdem(old_states, &
       new_states, field, femdem_out)

    type(state_type), dimension(:), intent(inout) :: old_states, new_states
    type(scalar_field), intent(inout), optional :: field
    logical, intent(in), optional :: femdem_out

    type(state_type), dimension(size(old_states)) :: old_fields_state
    type(state_type), dimension(size(old_states)) :: new_fields_state
    type(vector_field), pointer :: old_position, new_position
    integer :: i

    ewrite(1, *) "In interpolation_galerkin_multiple_states_femdem"

    call collapse_fields_in_state(old_states, old_fields_state)
    call collapse_fields_in_state(new_states, new_fields_state)
    call derive_collapsed_bcs(new_states, new_fields_state, bctype = "dirichlet")

    old_position => extract_vector_field(old_states(1), "Coordinate")
    new_position => extract_vector_field(new_states(1), "Coordinate")

    call interpolation_galerkin_scalars(old_fields_state, old_position, &
         new_fields_state, new_position, solid=field, femdem_out=femdem_out)

    do i = 1, size(old_fields_state)
      call deallocate(old_fields_state(i))
      call deallocate(new_fields_state(i))
    end do

    ewrite(1, *) "Exiting interpolation_galerkin_multiple_states_femdem"

  end subroutine interpolation_galerkin_multiple_states_femdem

  !----------------------------------------------------------------------------

  subroutine galerkin_projection_inner_loop_femdem(ele_B, little_mass_matrix, &
       detJ, local_rhs, conservation_tolerance, stat, field_counts, old_fields, &
       old_position, new_fields, new_position, inversion_matrices_A, supermesh_shape, &
       femdem_out, solid)

    integer, intent(in) :: ele_B
    real, dimension(:,:,:), intent(inout) :: little_mass_matrix
    real, dimension(:), intent(out) :: detJ
    real, dimension(:,:,:), intent(inout) :: local_rhs
    real, intent(in) :: conservation_tolerance
    integer, intent(out) :: stat

    integer, dimension(:), intent(in) :: field_counts
    logical, intent(in) :: femdem_out

    type(scalar_field), dimension(:,:), intent(in) :: old_fields
    type(vector_field), intent(in) :: old_position

    type(scalar_field), dimension(:,:), intent(inout) :: new_fields
    type(vector_field), intent(in) :: new_position
    real, dimension(:, :, :) :: inversion_matrices_A
    type(element_type), intent(inout) :: supermesh_shape

    type(scalar_field), intent(inout) :: solid
    integer, dimension(:), pointer :: ele_B_nodes

    real, dimension(ele_loc(new_position, ele_B), ele_loc(new_position, ele_B)) :: inversion_matrix_B, inversion_matrix_A
    real, dimension(ele_ngi(new_position, ele_B)) :: detwei_B
    real, dimension(supermesh_shape%ngi) :: detwei_C

    real :: vol_B, vols_C, vols_Cc, all_vols_C
    integer :: ele_A, ele_C, nloc, dim, j, k, l, n, loc, field, mesh, mesh_count, nintersections
    type(vector_field) :: intersection
    type(element_type), pointer :: B_shape

    real, dimension(new_position%dim+1, supermesh_shape%ngi) :: pos_at_quad_B, pos_at_quad_A, tmp_pos_at_quad
    real, dimension(size(local_rhs, 3), supermesh_shape%ngi) :: basis_at_quad_B, basis_at_quad_A
    real, dimension(size(local_rhs, 3),size(local_rhs, 3)) :: mat, mat_int

    real, dimension(new_position%dim, supermesh_shape%ngi) :: intersection_val_at_quad
    real, dimension(new_position%dim, new_position%dim, ele_ngi(new_position, 1)) :: invJ
    real, dimension(new_position%dim, ele_loc(new_position, ele_B)) :: pos_B

    type(plane_type), dimension(4) :: planes_B
    type(tet_type) :: tet_A, tet_B
    integer :: lstat

    real, dimension(size(local_rhs, 3)) :: tmp_local_rhs, tmp_ele_val

    all_vols_C = 0.

    local_rhs = 0.0

    mesh_count = size(field_counts)
    dim = mesh_dim(new_position)

    if (dim == 3) then
      tet_B%V = ele_val(new_position, ele_B)
      planes_B = get_planes(tet_B)
    else
      pos_B = ele_val(new_position, ele_B)
    end if

    ! First thing: assemble and invert the inversion matrix.
    call local_coords_matrix(new_position, ele_B, inversion_matrix_B)
    inversion_matrix_B = transpose(inversion_matrix_B)

    ! Second thing: assemble the mass matrix of B on the left.
    call compute_inverse_jacobian(new_position, ele_B, invJ=invJ, detJ=detJ, detwei=detwei_B)

    do mesh = 1, mesh_count
      if(field_counts(mesh)>0) then
        B_shape => ele_shape(new_fields(mesh,1),1)
          nloc = B_shape%loc
          little_mass_matrix(mesh, :nloc, :nloc) = shape_shape(B_shape, B_shape, detwei_B)
      end if
    end do

    vol_B = sum(detwei_B)
    vols_C = 0.0

    call rtree_intersection_finder_find(new_position, ele_B)
    call rtree_intersection_finder_query_output(nintersections)

    ! loop over the intersecting elements
    do n = 1, nintersections

      call rtree_intersection_finder_get_output(ele_A, n)

      ! but we only need that mapping for this ele_B now, so just compute it now
      if (dim == 3 .and. (intersector_exactness .eqv. .false.)) then
        tet_A%V = ele_val(old_position, ele_A)
        call intersect_tets(tet_A, planes_B, supermesh_shape, stat=lstat, output=intersection)
        if (lstat == 1) cycle
      else
        intersection = intersect_elements(old_position, ele_A, pos_B, supermesh_shape)
      end if

#ifdef DUMP_SUPERMESH_INTERSECTIONS
      if (ele_count(intersection) /= 0) then
        call vtk_write_fields("intersection", dump_idx, intersection, intersection%mesh)
        dump_idx = dump_idx + 1
      end if
#endif

      vols_Cc=0. ! inside the nintersection loop

      ! Loop over the supermesh elements, evaluate the basis functions at the
      ! quadrature points and integrate.
      do ele_C=1,ele_count(intersection)
        intersection_val_at_quad = ele_val_at_quad(intersection, ele_C)
        ! Compute the local coordinates in ele_B of the quadrature points of ele_C:
        tmp_pos_at_quad(1:dim, :) = intersection_val_at_quad
        tmp_pos_at_quad(dim+1, :) = 1.0
#ifdef INLINE_MATMUL
        forall (j=1:dim+1)
          forall (k=1:supermesh_shape%ngi)
            pos_at_quad_B(j, k) = sum(inversion_matrix_B(:, j) * tmp_pos_at_quad(:, k))
          end forall
        end forall
#else
        pos_at_quad_B = matmul(inversion_matrix_B, tmp_pos_at_quad)
#endif

        ! Compute the local coordinates in ele_A of the quadrature points of ele_C:
        tmp_pos_at_quad(1:dim, :) = intersection_val_at_quad
        tmp_pos_at_quad(dim+1, :) = 1.0
#ifdef INLINE_MATMUL
        inversion_matrix_A = transpose(inversion_matrices_A(:, :, ele_A))
        forall (j=1:dim+1)
          forall (k=1:supermesh_shape%ngi)
            pos_at_quad_A(j, k) = sum(inversion_matrix_A(:, j) * tmp_pos_at_quad(:, k))
          end forall
        end forall
#else
        pos_at_quad_A = matmul(inversion_matrices_A(:, :, ele_A), tmp_pos_at_quad)
#endif

        call transform_to_physical(intersection, ele_C, detwei_C)

        vols_C = vols_C + sum(detwei_C)
        vols_Cc = vols_Cc + sum(detwei_C)

        do mesh = 1, mesh_count
          if(field_counts(mesh)>0) then
            B_shape => ele_shape(new_fields(mesh,1),1)
            nloc = B_shape%loc
            ! This is an inlined eval_shape, optimised for P0 and P1
            ! Evaluate the basis functions at the local coordinates
            basis_at_quad_A = 0.0
            basis_at_quad_B = 0.0
            if (element_degree(new_fields(mesh,1),ele_B)==0) then
              basis_at_quad_A(:nloc,:) = 1.0
              basis_at_quad_B(:nloc,:) = 1.0
            elseif (element_degree(new_fields(mesh,1),ele_B)==1) then
              basis_at_quad_A(:nloc,:) = pos_at_quad_A 
              basis_at_quad_B(:nloc,:) = pos_at_quad_B 
            else
              do loc=1,nloc
                do j=1,ele_ngi(intersection, ele_C)
                  basis_at_quad_A(loc, j) = eval_shape(B_shape, loc, pos_at_quad_A(:, j))
                  basis_at_quad_B(loc, j) = eval_shape(B_shape, loc, pos_at_quad_B(:, j))
                end do
              end do
            end if

            ! Combined outer_product and tensormul_3_1 to see if it is faster.
            ! This is sort of like a mixed shape_shape.
            ! Here we assemble a little local part of the mixed mass matrix.
            mat = 0.0
            mat_int = 0.0
            do j=1,ele_ngi(intersection, ele_C)
              forall (k=1:nloc,l=1:nloc)
                mat(k, l) = mat(k, l) + detwei_C(j) * basis_at_quad_B(k, j) * basis_at_quad_A(l, j)
              end forall
            end do

            ! And now we apply that to the field to give the RHS contribution to the Galerkin
            ! projection.
            do field=1,field_counts(mesh)
#ifdef INLINE_MATMUL
              tmp_ele_val(:nloc) = ele_val(old_fields(mesh,field), ele_A)
              forall (j=1:nloc)
                tmp_local_rhs(j) = sum(mat(j, :nloc) * tmp_ele_val(:nloc))
              end forall
              local_rhs(mesh,field,:nloc) = local_rhs(mesh,field,:nloc) + tmp_local_rhs(:nloc)
#else
              local_rhs(mesh,field,:nloc) = local_rhs(mesh,field,:nloc) +&
                                    matmul(mat(:nloc,:nloc), ele_val(old_fields(mesh,field), ele_A))
#endif
            end do
          end if
        end do

      end do  ! intersection loop, i.e. ele_C loop

      all_vols_C = all_vols_C + vols_Cc

      call deallocate(intersection)

   end do ! nintersection loop, i.e. ele_A loop 

   if (.not.femdem_out) then
     ele_B_nodes => ele_nodes(new_position, ele_B)
     do loc = 1, size(ele_B_nodes)
       call addto(solid, ele_B_nodes(loc), all_vols_C / vol_B)
     end do
   end if

   if (femdem_out) then
     ! Check for supermeshing failures.
     if (abs(vol_B - vols_C)/vol_B > conservation_tolerance .and. & 
#ifdef DOUBLEP
       & abs(vol_B - vols_C) > 100.0 * 1.0e-12) then
#else
       & abs(vol_B - vols_C) > 100.0 * epsilon(0.0)) then
#endif
       ewrite(0,*) 'sum(detwei_B) = ', vol_B, ', all sum(detwei_C) = ', vols_C
       stat = 1
     else
       stat = 0
     end if
   end if

  end subroutine galerkin_projection_inner_loop_femdem

  !----------------------------------------------------------------------------

  subroutine interpolation_galerkin_scalars(old_fields_state, old_position, new_fields_state, new_position, map_BA, force_bounded, solid, femdem_out)
    type(state_type), dimension(:), intent(in) :: old_fields_state
    type(vector_field), intent(in) :: old_position

    type(state_type), dimension(:), intent(inout) :: new_fields_state
    type(vector_field), intent(in) :: new_position
    type(ilist), dimension(:), intent(in), optional, target :: map_BA
    logical, intent(in), optional :: force_bounded

    integer :: ele_B
    integer :: ele_A
    integer :: name, no_names, priority, f, field, field2, max_field_count
    
    type(scalar_field), dimension(:,:), allocatable :: old_fields, new_fields
    integer, dimension(size(old_fields_state)) :: field_counts
    
    type(scalar_field), dimension(:,:), allocatable :: named_fields, named_rhs
    character(len=FIELD_NAME_LEN), dimension(:), allocatable :: field_names
    integer, dimension(:), allocatable :: named_counts, priorities, named_indices
    integer, dimension(:,:), allocatable :: tmp_named_indices

    ! We want to compute the mixed mass matrix M^{BA}.
    ! But that's huge. So, we compute a part of a row (not even the whole row)
    ! and multiply it by a part of the solution on the old mesh A
    ! to get a component of the RHS of the matrix we want to solve.
    real, dimension(:,:,:), allocatable :: local_rhs
    real, dimension(:,:), allocatable :: little_rhs
    type(scalar_field), dimension(:,:), allocatable :: rhs
    ! For each element in B, we will need to identify the local coordinates in B
    ! of the positions of the gauss points of all its children C elements.
    ! So we'll need to assemble and invert that matrix (the global-to-local inversion matrix):
    real, dimension(ele_loc(new_position, 1), ele_loc(new_position, 1), ele_count(old_position)) :: inversion_matrices_A
    real, dimension(:,:,:), allocatable :: little_mass_matrix
    real, dimension(:,:,:), allocatable :: little_inverse_mass_matrix, little_inverse_mass_matrix_copy

    integer :: dim
    type(ilist), dimension(:), pointer :: lmap_BA
    type(quadrature_type) :: supermesh_quad
    type(element_type) :: supermesh_shape
    real :: int_old, int_new, cons_err, current_time
    logical, dimension(size(old_fields_state)) :: dg

    type(csr_matrix), dimension(:), allocatable :: M_B
    type(csr_sparsity) :: M_B_sparsity
    type(scalar_field), dimension(:), allocatable :: M_B_L
    type(scalar_field) :: inverse_M_B_L
    
    ! Boundedness stuff
    logical, dimension(:,:), allocatable :: bounded, lumped
    logical, dimension(:), allocatable :: coupled
    type(scalar_field) :: bounded_soln, max_bound, min_bound
    type(csr_sparsity), pointer :: nnlist
    integer :: node_B
    integer, dimension(:), pointer :: patch

    real :: upper_bound, lower_bound
    integer, dimension(:), pointer :: ele_nodes_B
    integer :: stat, statp
    logical :: l_apply_globally, u_apply_globally
    
    logical :: l_force_bounded
    
    integer :: max_loc, max_degree, nloc
    integer :: mesh, mesh_count
    
    real :: conservation_tolerance, tmp_tol

    type(element_type), pointer :: shape_B
    real, dimension(ele_ngi(new_position, 1)) :: detJ
    integer :: j

    logical :: new_positions_simplicial
    type(integer_set), dimension(:,:), allocatable :: bc_nodes
    character(len=FIELD_NAME_LEN) :: bctype
    integer, dimension(:), pointer :: surface_node_list
    logical, dimension(:, :), allocatable :: force_bc
    integer :: bc

    ! femdem
    logical :: femdem_in
    logical, intent(in), optional :: femdem_out
    integer :: ntests
    type(scalar_field), intent(inout), optional :: solid

    ewrite(1, *) "In interpolation_galerkin_scalars"

    stat = 0
    if(present(force_bounded)) then
      l_force_bounded = force_bounded
    else
      l_force_bounded = .false.
    end if
    
    ! Linear positions -- definitely linear positions.
    assert(old_position%mesh%shape%degree == 1)
    assert(continuity(old_position) >= 0)
    assert(continuity(new_position) >= 0)
    
    mesh_count = size(old_fields_state)
    max_field_count = 0
    field_counts = 0
    do mesh = 1, mesh_count
      field_counts(mesh) = scalar_field_count(old_fields_state(mesh))
      max_field_count = max(max_field_count, scalar_field_count(old_fields_state(mesh)))
    end do
    allocate(bounded(mesh_count, max_field_count))
    bounded = .false.
    allocate(lumped(mesh_count, max_field_count))
    lumped = .false.
    allocate(old_fields(mesh_count, max_field_count))
    allocate(new_fields(mesh_count, max_field_count))
    allocate(force_bc(mesh_count, max_field_count))
    allocate(bc_nodes(mesh_count, max_field_count))
    
    shape_B => ele_shape(new_position, 1)
    new_positions_simplicial = (shape_B%numbering%family == FAMILY_SIMPLEX)
    
    dim = mesh_dim(new_position)

    dg = .false.
    max_degree = 0
    max_loc = 0
    conservation_tolerance = 1.0
    do mesh = 1, size(old_fields_state)
      if(field_counts(mesh)>0) then
      
        do field = 1, field_counts(mesh)
          old_fields(mesh, field) = extract_scalar_field(old_fields_state(mesh), field)
          new_fields(mesh, field) = extract_scalar_field(new_fields_state(mesh), field)
          call zero(new_fields(mesh, field))
          bounded(mesh, field) = l_force_bounded.or.&
                          have_option(trim(complete_field_path(new_fields(mesh,field)%option_path, stat=statp)) // &
                                                "/galerkin_projection/continuous/bounded[0]")
          lumped(mesh, field) = have_option(trim(complete_field_path(new_fields(mesh,field)%option_path, stat=statp)) // &
                                                "/galerkin_projection/continuous/lump_mass_matrix")
          call get_option(trim(complete_field_path(new_fields(mesh,field)%option_path, stat=statp)) // &
                                                "/galerkin_projection/supermesh_conservation/tolerance", tmp_tol, default = 0.001)
          ! Let's check for a relative area/volume loss of 0.1% if none is specified
          conservation_tolerance = min(conservation_tolerance, tmp_tol)

          force_bc(mesh, field) = have_option(trim(complete_field_path(new_fields(mesh,field)%option_path, stat=statp)) &
            & // "/galerkin_projection/honour_strong_boundary_conditions")
          if (force_bc(mesh, field)) then

            if (.not. has_boundary_condition(new_fields(mesh, field), "dirichlet")) then
              ewrite(0, *) "Warning: For field: " // trim(new_fields(mesh, field)%name)
              ewrite(0, *) "Warning: Asked to honour strong boundary conditions through the Galerkin projection without any such BCs"
            end if

            call set_dirichlet_consistent(new_fields(mesh, field))

            call allocate(bc_nodes(mesh, field))
            do bc=1, get_boundary_condition_count(new_fields(mesh, field))
              call get_boundary_condition(new_fields(mesh, field), bc, type=bctype, surface_node_list=surface_node_list)
              if (trim(bctype) == "dirichlet") then
                call insert(bc_nodes(mesh, field), surface_node_list)
              end if
            end do
          end if
        end do
        
        dg(mesh) = (continuity(new_fields(mesh,1)) < 0)
        if(dg(mesh)) then
          bounded(mesh,:) = .false. ! not possible to have a bounded or lumped dg field 
          lumped(mesh,:) = .false.  ! so just to make sure set it to false
        end if
        
        max_degree = max(max_degree, element_degree(new_fields(mesh,1), 1))
        max_loc = max(max_loc, ele_loc(new_fields(mesh,1), 1))
      
      end if
    end do
    
    allocate(local_rhs(mesh_count, max_field_count, max_loc))
    allocate(little_mass_matrix(mesh_count, max_loc, max_loc))

    if (any(dg).and.new_positions_simplicial) then
      allocate(little_inverse_mass_matrix(mesh_count, max_loc, max_loc))
      allocate(little_inverse_mass_matrix_copy(mesh_count, max_loc, max_loc))
      little_inverse_mass_matrix = 0.0
      do mesh=1,mesh_count
        if((field_counts(mesh)>0).and.dg(mesh)) then
          shape_B => ele_shape(new_fields(mesh, 1), 1)
          nloc = ele_loc(new_fields(mesh, 1), 1)
          little_inverse_mass_matrix(mesh, :nloc, :nloc) = shape_shape(shape_B, shape_B, shape_B%quadrature%weight)
          call invert(little_inverse_mass_matrix(mesh, :nloc, :nloc))
        end if
      end do
    end if
    
    allocate(little_rhs(max_loc, max_field_count))

    if(any(.not.dg)) then
      ! if any meshes are not dg then we need a lhs matrix and a global rhs
      
      allocate(rhs(mesh_count, max_field_count))
      allocate(M_B(mesh_count))
      allocate(M_B_L(mesh_count))
      
      do mesh = 1, mesh_count
        if(.not.dg(mesh)) then
          if(field_counts(mesh)>0) then
            do field = 1, field_counts(mesh)
              call allocate(rhs(mesh,field), new_fields(mesh,field)%mesh, name = trim(new_fields(mesh,field)%name)//"RHS")
              call zero(rhs(mesh,field))
            end do
      
            if(.not.all(lumped(mesh,1:field_counts(mesh)))) then
              M_B_sparsity = make_sparsity(new_fields(mesh,1)%mesh, new_fields(mesh,1)%mesh, name="MassMatrixBSparsity")
            
              call allocate(M_B(mesh), M_B_sparsity, &
                            name=trim(new_fields(mesh,1)%mesh%name)//"MassMatrixB")
              call zero(M_B(mesh))
              
              call deallocate(M_B_sparsity)
            end if
            
            if(any(bounded(mesh,:)).or.any(lumped(mesh,:))) then
              call allocate(M_B_L(mesh), new_fields(mesh,1)%mesh, &
                            name=trim(new_fields(mesh,1)%mesh%name)//"LumpedMassMatrixB")
              call zero(M_B_L(mesh))
            end if
          end if
        end if
      end do
      
    end if
    
    supermesh_quad = make_quadrature(vertices=ele_loc(new_position, 1), dim=dim, degree=max(max_degree+max_degree, 1))
    supermesh_shape = make_element_shape(vertices=ele_loc(new_position, 1), dim=dim, degree=1, quad=supermesh_quad)

    ! figure out if this is a femdem interpolation
    femdem_in = .false.
    if (present(solid)) then
       call zero(solid)
       femdem_in = .true.
    end if
    if (femdem_in .or. present_and_true(femdem_out)) call rtree_intersection_finder_set_input(old_position)

    call intersector_set_dimension(dim)
    if (present(map_BA)) then
      lmap_BA => map_BA
    else if (.not.(femdem_in .or. present_and_true(femdem_out))) then
      allocate(lmap_BA(ele_count(new_position)))
      lmap_BA = intersection_finder(new_position, old_position)
    end if

    do ele_A=1,ele_count(old_position)
      call local_coords_matrix(old_position, ele_A, inversion_matrices_A(:, :, ele_A))
    end do

#ifdef DUMP_SUPERMESH_INTERSECTIONS
    call system("rm intersection*.vtu")
    dump_idx = 0
#endif

    ewrite(1, *) "Entering supermeshing loop"

      ewrite(1, *) "   ...for femdem"
      do ele_B=1,ele_count(new_position)

        call galerkin_projection_inner_loop_femdem(ele_B, little_mass_matrix, detJ, local_rhs, conservation_tolerance, stat, &
                                                     field_counts, old_fields, old_position, new_fields, new_position, &
                                                     inversion_matrices_A, supermesh_shape, present_and_true(femdem_out), solid)

        if (stat /= 0) then
          ! Uhoh! We haven't found all the mass for ele_B :-/
          ! The intersector has missed something (almost certainly due to
          ! finite precision arithmetic). Geometry is hard!
          ! So let's go all arbitrary precision on its ass.
          ! Data, Warp 0!
#ifdef HAVE_LIBCGAL
          ewrite(0,*) "Using CGAL to try to fix conservation error"
          call intersector_set_exactness(.true.)
          call galerkin_projection_inner_loop(ele_B, little_mass_matrix, detJ, local_rhs, conservation_tolerance, stat, &
               field_counts, old_fields, old_position, new_fields, new_position, &
               lmap_BA, inversion_matrices_A, supermesh_shape)
          if(stat/=0) then
             ewrite(0,*) "Sorry, CGAL failed to fix conservation error."
          end if
          call intersector_set_exactness(.false.)
#else
          ewrite(0,*) "Warning: it appears a supermesh intersection wasn't found resulting in a conservation error."
          ewrite(0,*) "Recompile with CGAL if you want to try to fix it."
#endif
        end if

        do mesh = 1, mesh_count
          if(field_counts(mesh)>0) then
            nloc = ele_loc(new_fields(mesh,1),1)
            ele_nodes_B => ele_nodes(new_fields(mesh,1), ele_B)
            if(dg(mesh)) then
              little_rhs = 0.0
              do field=1,field_counts(mesh)
                little_rhs(:nloc, field) = local_rhs(mesh,field,:nloc)
              end do

              if (any(force_bc(mesh,1:field_counts(mesh)))) then
                little_inverse_mass_matrix_copy=little_inverse_mass_matrix
              end if
            
              if (new_positions_simplicial) then
                do field=1,field_counts(mesh)
                  if (force_bc(mesh,field)) then
                    if (any(has_value(bc_nodes(mesh,field), ele_nodes_B))) then
                      local_rhs(mesh, field, :nloc)=local_rhs(mesh, field, :nloc)-matmul( little_mass_matrix(mesh,:nloc,:nloc)*abs(detJ(1)), ele_val(new_fields(mesh,field), ele_B) )
                      little_inverse_mass_matrix = little_inverse_mass_matrix_copy
                      do j=1, nloc
                        if (has_value(bc_nodes(mesh,field), ele_nodes_B(j))) then
                          little_inverse_mass_matrix(mesh, j,:)=0.0
                          little_inverse_mass_matrix(mesh, :,j)=0.0
                          little_inverse_mass_matrix(mesh, j,j)=1.0
                          local_rhs(mesh, field, j)=node_val(new_fields(mesh,field), ele_nodes_B(j))
                        end if
                      end do
                    end if
                  end if
#ifdef INLINE_MATMUL
                  forall (j=1:nloc)
                    little_rhs(j, field) = sum(little_inverse_mass_matrix(mesh, j, :nloc) * local_rhs(mesh, field, :nloc))
                  end forall
                  little_rhs(:nloc, field) = little_rhs(:nloc, field) / abs(detJ(1))
#else
                  little_rhs(:nloc, field) = matmul(little_inverse_mass_matrix(mesh, :nloc, :nloc) / abs(detJ(1)), little_rhs(:nloc, field))
#endif
                end do
              else
                call solve(little_mass_matrix(mesh,:nloc,:nloc), little_rhs(:nloc,:field_counts(mesh)))
              end if
            
              if (any(force_bc(mesh,1:field_counts(mesh)))) then
                little_inverse_mass_matrix=little_inverse_mass_matrix_copy
              end if

              do field = 1, field_counts(mesh)
                call set(new_fields(mesh,field), ele_nodes_B, little_rhs(:nloc, field))
              end do

            else
    
              do field=1,field_counts(mesh)
                call addto(rhs(mesh,field), ele_nodes_B, local_rhs(mesh,field,:nloc))
              end do
          
              if(.not.all(lumped(mesh,1:field_counts(mesh)))) then
                call addto(M_B(mesh), ele_nodes_B, ele_nodes_B, little_mass_matrix(mesh,:nloc,:nloc))
              end if
          
              if(any(bounded(mesh,:)).or.any(lumped(mesh,:))) then
                call addto(M_B_L(mesh), ele_nodes_B, sum(little_mass_matrix(mesh,:nloc,:nloc), 2))
              end if
            end if
          end if

        end do

      end do

      ewrite(1, *) "Supermeshing complete"

    do mesh = 1, mesh_count
      if(field_counts(mesh)>0) then
        if(.not.dg(mesh)) then
        
          if(any(bounded(mesh,:)).or.any(lumped(mesh,:))) then
            call allocate(inverse_M_B_L, M_B_L(mesh)%mesh, "InverseLumpedMass")
            call invert(M_B_L(mesh), inverse_M_B_L)
          end if
          
          do field=1,field_counts(mesh)
            if(lumped(mesh,field)) then
              call set(new_fields(mesh, field), rhs(mesh, field))
              call scale(new_fields(mesh,field), inverse_M_B_L)
              call halo_update(new_fields(mesh,field))
            else
              if (force_bc(mesh, field))  then
                  call apply_dirichlet_conditions(M_B(mesh), rhs(mesh, field), new_fields(mesh, field))
              end if
              call petsc_solve(new_fields(mesh, field), M_B(mesh), rhs(mesh, field), &
                & option_path=trim(complete_field_path(new_fields(mesh,field)%option_path, stat=statp)) &
                & // "/galerkin_projection/continuous")
              if (force_bc(mesh, field))  then
                ! clean up the rows made inactive for the strong bcs
                call reset_inactive(M_B(mesh))
              end if
            end if
          end do

          if(any(bounded(mesh,:))) then
          
            nnlist => extract_nnlist(new_fields(mesh,1))
            
            ! Ok. All that above was more or less the same as Galerkin projection. Here is 
            ! where we bound.
            
            ! to be able to couple the fields together we first need to group the fields by name
            ! and order them by priority
            ! so... let's get the priorities
            allocate(priorities(field_counts(mesh)))
            priorities = 0
            do field = 1, field_counts(mesh)
              call get_option(trim(new_fields(mesh,field)%option_path)//"/prognostic/priority", priorities(field), default=0)
            end do
              
            ! let's allocate some space (too much in fact but it's our best guess) for the counts of each name
            allocate(named_counts(field_counts(mesh)))
            named_counts = 0
            ! the names themselves
            allocate(field_names(field_counts(mesh)))
            field_names = ""
            ! the indices of each name
            allocate(tmp_named_indices(field_counts(mesh), field_counts(mesh)))
            tmp_named_indices = 0
            
            ! now loop through the fields collecting the actual number of fields with the same
            ! names and where they are located (their indices) in the current lists
            f = 0
            do field=1,field_counts(mesh)
              if(bounded(mesh,field)) then
                if(any(new_fields(mesh,field)%name==field_names(:sum(named_counts)))) cycle
                f = f + 1
                field_names(f) = trim(new_fields(mesh,field)%name)
                named_counts(f) = 1
                tmp_named_indices(f,named_counts(f)) = field
                do field2=1,field_counts(mesh)
                  if(field==field2) cycle
                  if(trim(new_fields(mesh,field2)%name)==field_names(f)) then
                    named_counts(f) = named_counts(f) + 1
                    tmp_named_indices(f,named_counts(f)) = field2
                  end if
                enddo
              end if
            end do
            no_names = f
            
            ! allocate the real space for them (still too much to avoid ragged arrays)
            allocate(named_fields(no_names, maxval(named_counts)))
            allocate(named_rhs(no_names, maxval(named_counts)))
            allocate(named_indices(maxval(named_counts)))
            
            do name = 1, no_names
              ! sort out their indices in order of priority
              f = 0
              named_indices = 0
              do priority = maxval(priorities), minval(priorities), -1
                do field = 1, named_counts(name)
                  if(priorities(tmp_named_indices(name,field))==priority) then
                    f = f + 1
                    named_indices(f) = tmp_named_indices(name, field)
                  end if
                end do
              end do
              
              ! and finally put them into a new list of fields sorted by name
              do field = 1, named_counts(name)
                named_fields(name, field) = new_fields(mesh, named_indices(field))
                named_rhs(name, field) = rhs(mesh, named_indices(field))
              end do
            end do
            
            do name = 1, no_names
            
              allocate(coupled(named_counts(name)))
              coupled = .false.
              do field = 1, named_counts(name)
                coupled(field) = have_option(trim(complete_field_path(named_fields(name,field)%option_path, stat=statp))// &
                  & "/galerkin_projection/continuous/bounded[0]/bounds/upper_bound/coupled")
              end do
          
              do field=1,named_counts(name)
              
                ewrite(2,*) 'Bounding field:', trim(named_fields(name,field)%name)
            
                ! Step 0. Compute bounds
                call get_option(trim(complete_field_path(named_fields(name,field)%option_path, stat=statp))// &
                  & "/galerkin_projection/continuous/bounded[0]/bounds/upper_bound", &
                  & upper_bound, default=huge(0.0)*epsilon(0.0))
                  
                u_apply_globally = have_option(trim(complete_field_path(named_fields(name,field)%option_path, stat=statp))// &
                  & "/galerkin_projection/continuous/bounded[0]/bounds/upper_bound/apply_globally")
                  
                call get_option(trim(complete_field_path(named_fields(name,field)%option_path, stat=statp))// &
                  & "/galerkin_projection/continuous/bounded[0]/bounds/lower_bound", &
                  & lower_bound, default=-huge(0.0)*epsilon(0.0))
                  
                l_apply_globally = have_option(trim(complete_field_path(named_fields(name,field)%option_path, stat=statp))// &
                  & "/galerkin_projection/continuous/bounded[0]/bounds/lower_bound/apply_globally")
                
                if((.not.u_apply_globally).or.(coupled(field))) then
                  call allocate(max_bound, named_fields(name,1)%mesh, "MaxBound")
                else
                  call allocate(max_bound, named_fields(name,1)%mesh, "MaxBound", field_type=FIELD_TYPE_CONSTANT)
                end if
                
                if(.not.l_apply_globally) then
                  call allocate(min_bound, named_fields(name,1)%mesh, "MinBound")
                else
                  call allocate(min_bound, named_fields(name,1)%mesh, "MinBound", field_type=FIELD_TYPE_CONSTANT)
                end if
                
                call set(max_bound, upper_bound)
                if(coupled(field)) then
                  do field2 = 1, field-1
                    if(coupled(field2)) call addto(max_bound, named_fields(name,field2), -1.0)
                  end do
                end if
                
                call set(min_bound, lower_bound)
                
                call allocate(bounded_soln, named_fields(name,1)%mesh, "BoundedSolution")
                call set(bounded_soln, named_rhs(name,field))
                call scale(bounded_soln, inverse_M_B_L)
                call halo_update(bounded_soln)
                
                do node_B=1,node_count(named_fields(name,1)%mesh)
                  patch => row_m_ptr(nnlist, node_B)
                  if(.not.u_apply_globally) then
                    call set(max_bound, node_B, max(min(maxval(bounded_soln%val(patch)), &
                                                        node_val(max_bound, node_B)), &
                                                    lower_bound))
                  end if
                  if(.not.l_apply_globally) then
                    call set(min_bound, node_B, max(min(minval(bounded_soln%val(patch)), &
                                                        node_val(max_bound, node_B)), &
                                                    lower_bound))
                  end if
                end do

                call halo_update(max_bound)
                ewrite_minmax(max_bound)

                call halo_update(min_bound)
                ewrite_minmax(min_bound)
                
                call bound_field(named_fields(name, field), max_bound, min_bound, &
                                 M_B(mesh), M_B_L(mesh), inverse_M_B_L, bounded_soln, &
                                 new_position)

                
                call deallocate(max_bound)
                call deallocate(min_bound)
                call deallocate(bounded_soln)
                
              end do
              
              deallocate(coupled)
              
            end do
            
            deallocate(priorities)
            deallocate(named_counts)
            deallocate(field_names)
            deallocate(tmp_named_indices)
            deallocate(named_fields)
            deallocate(named_rhs)
            deallocate(named_indices)

          end if
          
          if(any(bounded(mesh,:)).or.any(lumped(mesh,:))) then
            call deallocate(inverse_M_B_L)
            call deallocate(M_B_L(mesh))
          end if

        end if
      
        do field = 1, field_counts(mesh)
          if(have_option(trim(complete_field_path(new_fields(mesh,field)%option_path, stat=statp)) // &
                                                "/galerkin_projection/supermesh_conservation/print_field_integral")) then
            int_old = field_integral(old_fields(mesh,field), old_position)
            int_new = field_integral(new_fields(mesh,field), new_position)
            cons_err = abs(int_old-int_new)/abs(int_old)
            ewrite(2,*) "relative change in field integral: ", cons_err, " for field ", trim(new_fields(mesh,field)%name)
            call get_option(trim(complete_field_path(new_fields(mesh,field)%option_path, stat=statp)) // &
                                                  "/galerkin_projection/supermesh_conservation/print_field_integral/tolerance", &
                                                  tmp_tol)
            if (cons_err > tmp_tol) then
              call get_option("/timestepping/current_time", current_time)
              ewrite(0,*) "Warning: relative conservation error: ", cons_err, " for field ", trim(old_fields(mesh,field)%name), " at time: ", current_time
              call vtk_write_fields(trim(new_fields(mesh,field)%name)//"_conservation_error", 0, old_position, old_fields(mesh,field)%mesh, sfields=(/old_fields(mesh,field)/))
              call vtk_write_fields(trim(new_fields(mesh,field)%name)//"_conservation_error", 1, new_position, new_fields(mesh,field)%mesh, sfields=(/new_fields(mesh,field)/))
            end if
          end if
        end do

      end if
      
    end do

    call deallocate(supermesh_shape)
    call deallocate(supermesh_quad)

    do mesh = 1, mesh_count
      if(field_counts(mesh)>0) then
        if(.not.dg(mesh)) then
          if(.not.all(lumped(mesh,1:field_counts(mesh)))) then
            call deallocate(M_B(mesh))
          end if
          do field = 1, field_counts(mesh)
            call deallocate(rhs(mesh,field))
            if (force_bc(mesh, field)) then
              call deallocate(bc_nodes(mesh, field))
            end if
          end do
        end if
      end if
    end do
    if(any(.not.dg).and.(max_field_count>0)) then
      deallocate(M_B)
      deallocate(M_B_L)
      deallocate(rhs)
    end if
    deallocate(bounded)
    deallocate(old_fields)
    deallocate(new_fields)
    deallocate(local_rhs)
    deallocate(little_mass_matrix)
    deallocate(force_bc)
    deallocate(bc_nodes)
    if(any(dg).and.new_positions_simplicial) then
      deallocate(little_inverse_mass_matrix)
      deallocate(little_inverse_mass_matrix_copy)
    end if
    deallocate(little_rhs)

    call finalise_tet_intersector
    call rtree_intersection_finder_reset(ntests)

    ewrite(1, *) "Exiting interpolation_galerkin_scalars"
    
  end subroutine interpolation_galerkin_scalars

  !----------------------------------------------------------------------------

  subroutine implicit_solids_register_diagnostic
    
    integer :: i, str_size, ndim
    character(len=254) :: fmt, buffer

    if (.not. have_option("/implicit_solids")) return
    
    ! figure out if we want to print out diagnostics and initialise files
    do_print_diagnostics = &
         have_option("/implicit_solids/one_way_coupling/print_diagnostics")
    do_print_multiple_solids_diagnostics = &
         have_option("/implicit_solids/one_way_coupling/multiple_solids/print_diagnostics")

    ! check for mutiple solids and get translation coordinates
    number_of_solids = 1
    multiple_solids = &
         have_option("/implicit_solids/one_way_coupling/multiple_solids")
    if (multiple_solids) call get_option( &
         "/implicit_solids/one_way_coupling/multiple_solids/number_of_solids", &
         number_of_solids)

    str_size=len_trim(int2str(number_of_solids))
    fmt="(I"//int2str(str_size)//"."//int2str(str_size)//")"

    have_temperature = &
         have_option("/material_phase[0]/scalar_field::Temperature")

    if (do_print_diagnostics) then
      call get_option("/geometry/dimension", ndim)
      call register_diagnostic(dim=ndim, name="Force", statistic="Value")

      if (do_print_multiple_solids_diagnostics) then
          do i = 1, number_of_solids
             write(buffer, fmt) i
             call register_diagnostic(dim=ndim, name="ForceOnSolid"//buffer, statistic="Value")
          end do
       end if

       if (have_temperature) then
          call register_diagnostic(dim=1, name="WallTemperature", statistic="Value")
          call register_diagnostic(dim=1, name="HeatTransfer", statistic="Value")
 
          if (do_print_multiple_solids_diagnostics) then
             do i = 1, number_of_solids
                write(buffer, fmt) i
                call register_diagnostic(dim=1, name="WallTemperatureOnSolid"//buffer, statistic="Value")
                call register_diagnostic(dim=1, name="HeatTransferAtSolid"//buffer, statistic="Value")
             end do
          end if 
 
       end if
    end if

  end subroutine implicit_solids_register_diagnostic

  !----------------------------------------------------------------------------

  subroutine implicit_solids_check_options
     integer :: ndim
     
     ! Get dimension:
     call get_option("/geometry/dimension", ndim)
     ! Check options for Implicit Solids:
     if (have_option("/implicit_solids/one_way_coupling") .and. ndim==1) then
        ewrite(-1,*) "Error: The 1-way Fluid-Structure Interactions are not supported for 1D simulations via Implicit Solids"
        FLExit("Use a 2D or 3D set-up when using a 1-way coupled simulation via implicit solids")
     else if (have_option("/implicit_solids/two_way_coupling") .and. (.not. ndim==3)) then
        ewrite(-1,*) "Error: The 2-way coupling of Fluidity/FEMDEM via Implicit Solids"
        ewrite(-1,*) "is only supported for 3D simulations."
        FLExit("Use 3D when using a 2-way coupled simulation via implicit solids")
     end if
  end subroutine implicit_solids_check_options

  !----------------------------------------------------------------------------

end module implicit_solids
