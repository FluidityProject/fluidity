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
  use global_parameters, only: FIELD_NAME_LEN, current_time, dt
  use spud
  use timeloop_utilities
  use fefields, only: compute_lumped_mass
  use parallel_tools
  use interpolation_module
  use mesh_files

  implicit none

  interface 
     subroutine y3allocate_femdem(string)          
       character(len=*) :: string
     end subroutine y3allocate_femdem
  end interface

  interface 
     subroutine y3dfemdem(flag, dt, &
          xs, ys, zs, us, vs, ws, uf, vf, wf)
       integer, intent(in) :: flag
       real, intent(in) :: dt
       real, dimension(:) :: xs, ys, zs
       real, dimension(:) :: us, vs, ws, uf, vf, wf
     end subroutine y3dfemdem
  end interface

  private
  public:: solids

contains

  subroutine solids(state, its)

    type(state_type),intent(inout) :: state
    character(len=field_name_len) :: external_mesh_name
    type(vector_field), pointer :: positions, velocity
    type(vector_field), pointer :: absorption
    type(vector_field), save :: external_positions
    type(scalar_field), pointer :: solid
    type(scalar_field) :: lumped_mass
    integer :: ele_A, ele_B, ele_C
    type(tet_type) :: tet_A, tet_B
    type(plane_type), dimension(:), allocatable :: planes_A
    integer :: stat, nintersections, i, j, k
    integer :: ntests, dat_unit, its
    integer, dimension(:), pointer :: ele_A_nodes
    type(vector_field) :: intersection
    real, dimension(:,:), allocatable :: pos_A
    real, dimension(:), allocatable :: detwei, drag
    real :: vol, ele_A_vol, sigma
    real, save :: beta
    integer, save :: itinoi
    logical, save :: init=.false.

    ewrite(3, *) "inside femdem"

    solid => extract_scalar_field(state, "SolidConcentration")
    call zero(solid)

    positions => extract_vector_field(state, "Coordinate")

    if (.not. init) then
       call get_option("/implicit_solids/one_way_coupling/mesh_name", &
            external_mesh_name)
       external_positions = extract_vector_field(state, &
            trim(external_mesh_name)//"Coordinate")

       call get_option("/implicit_solids/beta", beta, default=1.)

       assert(positions%dim >= 2)
       assert(positions%dim == external_positions%dim)

       if (GetRank() == 0) then
          dat_unit = free_unit()
          open(dat_unit, file="drag_force", status="replace")
          close(dat_unit)                  
       end if

       call get_option("/timestepping/nonlinear_iterations", itinoi)

       init=.true.
    end if

    call rtree_intersection_finder_set_input(positions)

    do ele_B = 1, ele_count(external_positions)

       call rtree_intersection_finder_find(external_positions, ele_B)
       call rtree_intersection_finder_query_output(nintersections)

       if (positions%dim == 3) then
          tet_B%v = ele_val(external_positions, ele_B)
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
                  ele_shape(external_positions, ele_B), &
                  stat=stat, output=intersection)

             deallocate(planes_A)

          else

             allocate(pos_A(positions%dim, ele_loc(positions, ele_A)))
             pos_A = ele_val(positions, ele_A)
             intersection = intersect_elements(external_positions, &
                  ele_B, pos_A, ele_shape(external_positions, ele_B))
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
             call addto(solid, ele_A_nodes(k), vol/ele_A_vol)
          end do

          call deallocate(intersection)

       end do

    end do

    do i = 1, node_count(solid)
       call set(solid, i, max(0., min(1., node_val(solid, i))))
    end do

    ewrite_minmax(solid)

    absorption => extract_vector_field(state, "VelocityAbsorption")
    call zero(absorption)

    do i = 1, node_count(absorption)
       sigma = node_val(solid, i)*beta/dt
       do j = 1, absorption%dim
          call set(absorption, j, i, sigma)
       end do
    end do

    ! total solid drag force...
    velocity => extract_vector_field(state, "Velocity")

    call allocate(lumped_mass, positions%mesh, "Lumped mass")
    call compute_lumped_mass(positions, lumped_mass)    

    allocate(drag(positions%dim))
    drag = 0.

    do i = 1, node_count(positions)
       do j = 1, positions%dim
          drag(j) = drag(j) + &
               node_val(absorption, 1, i) * node_val(velocity, j, i) * &
               node_val(lumped_mass, i)
       end do
    end do

    call deallocate(lumped_mass)

    call allsumv(drag)

    if (GetRank() == 0 .and. its == itinoi) then
       open(1453, file="drag_force", position="append")
       write(1453, *) &
            (drag(i), i = 1, positions%dim)
       close(1453)
    end if

    deallocate(drag)

    if(simulation_completed(current_time)) call deallocate(external_positions)
    call finalise_tet_intersector
    call rtree_intersection_finder_reset(ntests)

    ewrite(3, *) "leaving femdem"

  end subroutine solids

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  subroutine femdem_two_way_initialise(state)

    integer :: flag = 1
    character(len=field_name_len) :: external_mesh_name
    type(vector_field), pointer :: external_positions
    type(vector_field), pointer :: ext_pos_solid_vel, ext_pos_fluid_vel
    type(state_type) :: state


    call get_option("/implicit_solids/two_way_coupling/mesh_name", &
         external_mesh_name)

    !call y3allocate_femdem(trim(external_mesh_name)//char(0))

    ! get external positions...
    external_positions => extract_vector_field(state, &
         trim(external_mesh_name)//"Coordinate")

    ! this is the solid velocity on the solid mesh
    ext_pos_solid_vel => extract_vector_field(state, "ParticleVector")
    call zero(ext_pos_solid_vel)

    ! this is the interpolated fluid velocity
    ! on the solid mesh
    ext_pos_fluid_vel => extract_vector_field(state, "ParticleForce")
    call zero(ext_pos_fluid_vel)

    assert(node_count(external_positions) == node_count(ext_pos_solid_vel))
    assert(node_count(ext_pos_fluid_vel) == node_count(ext_pos_solid_vel))

    ! out :: ext_pos_solid_vel (i.e. initial velocity source)
    !call y3dfemdem(flag, dt, &
    !     external_positions%val(1)%ptr, external_positions%val(2)%ptr, external_positions%val(3)%ptr, &
    !     ext_pos_solid_vel%val(1)%ptr, ext_pos_solid_vel%val(2)%ptr, ext_pos_solid_vel%val(3)%ptr, &
    !     ext_pos_fluid_vel%val(1)%ptr, ext_pos_fluid_vel%val(2)%ptr, ext_pos_fluid_vel%val(3)%ptr)
    
  end subroutine femdem_two_way_initialise

  !----------------------------------------------------------------------------

!!$  subroutine femdem_two_way_update(state)
!!$
!!$    integer :: flag = 0
!!$    type(vector_field), pointer :: external_positions
!!$    type(vector_field), pointer :: ext_pos_solid_vel, ext_pos_fluid_vel
!!$    type(state_type) :: state
!!$
!!$    ! in  :: ext_pos_fluid_vel
!!$    ! out :: updated external_positions and ext_pos_solid_vel
!!$    call y3dfemdem(flag, dt, &
!!$         external_positions%val(1)%ptr, external_positions%val(2)%ptr, external_positions%val(3)%ptr, &
!!$         ext_pos_solid_vel%val(1)%ptr, ext_pos_solid_vel%val(2)%ptr, ext_pos_solid_vel%val(3)%ptr, &
!!$         ext_pos_fluid_vel%val(1)%ptr, ext_pos_fluid_vel%val(2)%ptr, ext_pos_fluid_vel%val(3)%ptr)
!!$
!!$  end subroutine femdem_two_way_update

  !----------------------------------------------------------------------------

  subroutine femdem_interpolation(state, operation)
 
    character(len=field_name_len) :: external_mesh_name
    character(len=field_name_len), intent(in) :: operation
    type(state_type) :: state
    type(state_type) :: alg_ext, alg_fl

    type(mesh_type), pointer :: ext_mesh, fl_mesh
    type(vector_field) :: ext_positions, fl_positions
    type(vector_field), pointer :: field_ext, field_fl

    ! this subroutine interpolates velocities
    ! between fluidity and femdem meshes
    ! if operation == "in"  : interpolate femdem into fluidity
    ! if operation == "out" : interpolate fluidity into femdem

    ! read in femdem data into alg_old state
    ! all data is on the solid mesh
    call get_option("/implicit_solids/two_way_coupling/mesh_name", &
         external_mesh_name)
    ext_mesh => extract_mesh(state, &
         trim(external_mesh_name)//"Coordinate")
    ext_positions = extract_vector_field(state, &
         trim(external_mesh_name)//"Coordinate")
    call insert(alg_ext, ext_mesh, "Mesh")
    call insert(alg_ext, ext_positions, "Coordinate")

    if (operation == "in") then

       ! this is the solid velocity on the solid mesh
       field_ext => extract_vector_field(state, "ParticleVector")
       ! just rename "ParticleVector" to "SolidVelocity"
       call insert(alg_ext, field_ext, "SolidVelocity")

    else if (operation == "out") then

       ! this is the fluid velocity on the solid mesh
       ! this will be returned to femdem...
       field_ext => extract_vector_field(state, "ParticleForce")
       call zero(field_ext)
       ! just rename "ParticleForce" to "Velocity"
       call insert(alg_ext, field_ext, "Velocity")

    end if

    ! read in fluidity data into alg_new state
    ! all data is on the fluidity mesh
    fl_mesh => extract_mesh(state, "Coordinate")
    fl_positions = extract_vector_field(state, "Coordinate")
    call insert(alg_fl, fl_mesh, "Mesh")
    call insert(alg_fl, fl_positions, "Coordinate")

    if (operation == "in") then

       ! this is the solid velocity on the fluidity mesh
       ! this will be used to set the source term...
       field_fl => extract_vector_field(state, "SolidVelocity")
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

    else if (operation == "out") then

       ! interpolate the fluid velocity from
       ! the fluid mesh to the solid mesh
       call linear_interpolation(alg_fl, alg_ext, different_domains=.true.)

    else
       FLExit("Don't know what to interpolate...")
    end if
    
    call deallocate(alg_ext)
    call deallocate(alg_fl)
    
  end subroutine femdem_interpolation

!!$  !----------------------------------------------------------------------------
!!$!  subroutine femdem_source_and_absorption(state)
!!$!  end subroutine femdem_source_and_absorption
!!$  !----------------------------------------------------------------------------

end module implicit_solids
