!  Copyright (C) 2006 Imperial College London and others.
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
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
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

module fsi_projections

  use supermesh_construction
  use solvers
  use global_parameters, only: OPTION_PATH_LEN
  use linked_lists
  use bound_field_module
  use fefields, only: compute_lumped_mass
  use tetrahedron_intersection_module

  implicit none

  public :: fsi_one_way_galerkin_projection, fsi_one_way_grandy_interpolation

contains

  subroutine fsi_one_way_galerkin_projection(fieldF, positionsF, positionsS, alpha_sf)
    !! Return the volume fraction scalar field by projecting unity from the supermesh to the fluid mesh
    !! Since positionsF and positionsS are different, we need to supermesh!
    !! positionsF and positionsS are the coordinate fields of the fluid and solid mesh respectively

    ! F stands for fluid, S for solid, C for the supermesh in this case:
    type(vector_field), intent(inout) :: fieldF
    type(vector_field), intent(inout) :: positionsF, positionsS
    type(scalar_field), intent(inout) :: alpha_sf

    ! this is for using the advancing front algorithm, therefore commented for now:
    !type(ilist), dimension(ele_count(positionsS)) :: map_SF 
    integer :: ele_F, ele_S

    type(quadrature_type) :: supermesh_quad
    type(element_type) :: supermesh_positions_shape, supermesh_field_shape

    type(vector_field) :: supermesh
    type(mesh_type) :: supermesh_field_mesh

    ! Scalar fields for alpha:
    type(scalar_field) :: alpha_sf_on_solid
    type(scalar_field) :: alpha_sf_on_supermesh

    type(csr_sparsity) :: sparsity_fluid
    type(csr_matrix) :: mass_matrix_fluid
    type(scalar_field) :: rhs_alpha_sf

    real, dimension(ele_loc(positionsS, 1), ele_loc(positionsS, 1)) :: inversion_matrix_S
    real, dimension(ele_loc(positionsS, 1), ele_loc(positionsS, 1), ele_count(positionsF)) :: inversion_matrices_F
    integer :: dim, max_degree

    ! Variables for ilist when using rtree intersection finder:
    type(ilist) :: map_SF_rtree
    integer :: nele_fs, j, maplen, ntests

    character(len=OPTION_PATH_LEN) :: tmp

    integer :: stat

    ! As the projections are bounded, we need to define lumped versions of the mass matrices:
    !type(scalar_field) :: lumped_mass_matrix_solid, lumped_mass_matrix_fluid
    !type(scalar_field) :: lumped_inverse_mass_matrix_solid, lumped_inverse_mass_matrix_fluid
    type(scalar_field) :: lumped_mass_matrix_fluid
    type(scalar_field) :: lumped_inverse_mass_matrix_fluid


    ewrite(2,*) "inside one_way_unity_projection"

    dim = mesh_dim(positionsS)
    call intersector_set_dimension(dim)

    ! Define polynomial degree and quadrature of the supermesh:
    max_degree = max(element_degree(fieldF, 1), element_degree(positionsS, 1))
    supermesh_quad = make_quadrature(vertices=ele_loc(positionsS, 1), dim=dim, degree=2*max_degree)
    ! Define the shape function of the supermesh:
    supermesh_positions_shape = make_element_shape(vertices=ele_loc(positionsS, 1), dim=dim, degree=1, quad=supermesh_quad)
    supermesh_field_shape = make_element_shape(vertices=ele_loc(positionsS, 1), dim=dim, degree=max_degree, quad=supermesh_quad)
    ! For each element in S get the list of intersecting elements in F
    ! Replace advancing_front_intersection_finder for now
    ! (more work needs to be done to make this work in parallel for fsi-modelling)
    ! with the more robust rtree_intersection_finder:
    ! advancing_front:
    !map_SF = intersection_finder(positionsS, positionsF)
    ! rtree:
    call rtree_intersection_finder_set_input(positionsF)


    sparsity_fluid = make_sparsity(fieldF%mesh, fieldF%mesh, "FluidMassMatrixSparsity")
    call allocate(mass_matrix_fluid, sparsity_fluid, name="FluidMassMatrix")
    call zero(mass_matrix_fluid)
    call deallocate(sparsity_fluid)

    ! Allocate lumped versions of the mass matrices:
    call allocate(lumped_mass_matrix_fluid, fieldF%mesh, "FluidLumpedMassMatrix")
    call zero(lumped_mass_matrix_fluid)
    call compute_lumped_mass(positionsF, lumped_mass_matrix_fluid)
    call allocate(lumped_inverse_mass_matrix_fluid, fieldF%mesh, "FluidLumpedInverseMass")
    call invert(lumped_mass_matrix_fluid, lumped_inverse_mass_matrix_fluid)

    ! get the matrix for the coordinates of (fluid) mesh F:
    do ele_F=1,ele_count(positionsF)
      call local_coords_matrix(positionsF, ele_F, inversion_matrices_F(:, :, ele_F))
      call assemble_mass_matrix(mass_matrix_fluid, fieldF, positionsF, ele_F)
    end do

    ! Fields for Alpha:
    call allocate(alpha_sf_on_solid, positionsS%mesh, "AlphaSSSolidMesh")
    ! alpha_sf on the solid mesh is unity by definition:
    call set(alpha_sf_on_solid, 1.0)
    call allocate(rhs_alpha_sf, fieldF%mesh, "AlphaSFOnFluidRHS")
    call zero(rhs_alpha_sf)

    ! Looping over all the elements of (solid) mesh S:
    do ele_S=1,ele_count(positionsS)

      ! =====================================================================================================
      ! The following is obviously only required when using the rtree intersection finder:
      ! This will be replaced by code using the advancing front algorithm, once this has been modified to
      ! work in parallel for fluid-solid interaction modelling
      ! Via RTREE, find intersection of solid element 'ele_S' with the
      ! input mesh (=fluid mesh 'positionsF'):
      call rtree_intersection_finder_find(positionsS, ele_S)
      ! Fetch output, the number of intersections of this solid element
      ! with the fluid mesh (= positionsF):
      call rtree_intersection_finder_query_output(nele_fs)
      ! Generate an ilist of elements in positionsF that intersect with ele_S:
      do j=1, nele_fs
        ! Get the donor (fluid) element which intersects with ele_S
        call rtree_intersection_finder_get_output(ele_F, j)
        ! insert_ascending works, but maybe is not necessary
        !call insert_ascending(map_SF_rtree, ele_F)
        call insert(map_SF_rtree, ele_F)
      end do
      ! =====================================================================================================

      ! get the matrix for the coordinates of (solid) mesh S:
      call local_coords_matrix(positionsS, ele_S, inversion_matrix_S)
      ! Construct the supermesh associated with ele_S. (For advancing front algorithm, uncomment the following line:)
      !call construct_supermesh(positionsS, ele_S, positionsF, map_SF(ele_S), supermesh_positions_shape, supermesh)
      ! =====================================================================================================
      ! When using the rtree intersection finder:
      if (map_SF_rtree%length > 0) then
        stat = 0
        call construct_supermesh(positionsS, ele_S, positionsF, map_SF_rtree, supermesh_positions_shape, supermesh, stat=stat)
        ! if no intersection for proc x was found, no supermesh was created, then stat/=0
        if (stat /= 0) then
          cycle
        end if
        ! =====================================================================================================

        ! At this point, a portion of the supermesh is constructed, which refines the dimensional space
        ! of the intersection of ele_S with mesh F.
        ! This portion of the supermesh is stored in 'supermesh'

        ! 1st step: Compute the project unity from the solid mesh to the supermesh:

        ! Allocate field for alpha on the supermesh:
        supermesh_field_mesh = make_mesh(supermesh%mesh, supermesh_field_shape, -1, "SupermeshFieldMesh")
        call allocate(alpha_sf_on_supermesh, supermesh_field_mesh, "AlphaSFOnSupermesh")
        call zero(alpha_sf_on_supermesh)

        ! Get alpha for this portion of the supermesh:
        call compute_alpha_on_supermesh(supermesh_field_shape, positionsS, ele_S, ele_val(alpha_sf_on_solid, ele_S), &
                                        supermesh, inversion_matrix_S, alpha_sf_on_supermesh)

        ! 2nd step: Compute the RHS of the Galerkin projection:
        ! Elemental Values of alpha of this part of the supermesh are now computed.
        ! Compute the RHS of 
        ! M_f * alpha_f = M_{f,sup} * alpha_{sup}
        ! and afterwards, the supermesh can be deleted as it is no longer needed,
        ! remember, (M_{f,sup} * alpha_{sup}) lives on the fluid mesh, not the supermesh.
        call compute_rhs_galerkin_projection_oneway_fsi(ele_S, rhs_alpha_sf, positionsF, positionsS, alpha_sf_on_supermesh, &
                                                      supermesh, inversion_matrices_F)

        ! Portion of supermesh no longer needed
        call deallocate(supermesh)
        call deallocate(supermesh_field_mesh)
        call deallocate(alpha_sf_on_supermesh)
      end if
      ! Flush ilist of intersecting elements (when using the rtree intersection finder):
      !call flush_list(map_SF_rtree)
      call deallocate(map_SF_rtree)
    end do

    call deallocate(supermesh_quad)
    call deallocate(supermesh_positions_shape)
    call deallocate(supermesh_field_shape)
    ! Because of using the rtree intersection finder:
    call rtree_intersection_finder_reset(ntests)
    !call deallocate(map_SF_rtree)
    ! =====================================================================================================
    ! the following is commented because of using rtree instead of the advancing front algorithm:
    !do ele_S=1,ele_count(positionsS)
    !  call deallocate(map_SF(ele_S))
    !end do
    ! =====================================================================================================

    ! 3rd step: Project alpha from the supermesh to the fluid and solid mesh:
    ! loop over fluid and solid elements and solve the last equation, which will
    ! project alpha from the supermesh to the fluid mesh

    ! Set solver options for the interpolations:
    tmp = alpha_sf%option_path
    alpha_sf%option_path = "/embedded_models/fsi_model/one_way_coupling/inter_mesh_projection/galerkin_projection/continuous"

    ! Project alpha to the fluid mesh:
    call petsc_solve(alpha_sf, mass_matrix_fluid, rhs_alpha_sf)

    ! Resetting option path for alpha_sf:
    alpha_sf%option_path = trim(tmp)

    ! Get options from option tree if projection is bounded or not:
    if (have_option('/embedded_models/fsi_model/one_way_coupling/inter_mesh_projection/galerkin_projection/continuous/bounded[0]')) then
    
       call bound_projection(alpha_sf, rhs_alpha_sf, mass_matrix_fluid, lumped_mass_matrix_fluid, lumped_inverse_mass_matrix_fluid, positionsF)
    end if

    call deallocate(mass_matrix_fluid)
    call deallocate(rhs_alpha_sf)
    call deallocate(lumped_mass_matrix_fluid)
    call deallocate(lumped_inverse_mass_matrix_fluid)

    ewrite(2,*) "leaving one_way_unity_projection"

  end subroutine fsi_one_way_galerkin_projection

  subroutine compute_alpha_on_supermesh(supermesh_field_shape, solid_positions, ele_B, alpha_sf_ele_value, &
                                        supermesh, inversion_matrix_B, alpha_sf_on_supermesh)
    ! F stands for fluid, S for solid, C for the supermesh in this case:
    type(vector_field), intent(in) :: solid_positions, supermesh
    type(element_type), intent(in), target :: supermesh_field_shape
    integer, intent(in) :: ele_B
    real, dimension(:), intent(in) :: alpha_sf_ele_value
    real, dimension(:, :), intent(in) :: inversion_matrix_B

    type(scalar_field), intent(inout) :: alpha_sf_on_supermesh

    type(mesh_type) :: supermesh_field_mesh

    integer :: ele_C
    integer, dimension(:), pointer :: nodes
    integer :: node_C, i, j, k

    integer :: ele_A
    real, dimension(ele_loc(solid_positions, ele_B)) :: local_coords
    integer :: dim
    real :: val
    type(vector_field) :: supermesh_positions_remapped

    dim = solid_positions%dim

    supermesh_field_mesh = make_mesh(supermesh%mesh, supermesh_field_shape, -1, "SupermeshFieldMesh")
    call allocate(supermesh_positions_remapped, dim, supermesh_field_mesh, "SupermeshPositionsRemapped")
    call remap_field(supermesh, supermesh_positions_remapped)

    ! Looping over all the elements of the portion of the supermesh:
    do ele_C=1,ele_count(supermesh)
      ! Get the global node numbers of element 'ele_C' of the supermesh:
      nodes => ele_nodes(alpha_sf_on_supermesh, ele_C)
      ! ele_A is then the region (number) of the supermesh, which describes the intersection of 
      ! one element of mesh A with one element of mesh B; the supermesh being the supermesh of all 
      ! the intersections of element 'ele_B' with mesh A.
      ! store the inverted matrix with the coordinates of mesh region 'ele_A' in inversion_matrix_A:

      ! Loop over the nodes of ele_C:
      do i=1,size(nodes)
        ! node_C = global node number of the i-th node of element ele_C:
        node_C = nodes(i)

        ! set local_coords to the values of node 'node_C' on the element 'ele_C' of the supermesh
        local_coords(1:dim) = node_val(supermesh_positions_remapped, node_C); local_coords(dim+1) = 1.0
        ! Compute the matrix multiplication of the coordinate matrix
        ! for region ele_F and local_coords, then store the result in local_coords:
        local_coords = matmul(inversion_matrix_B, local_coords)
        val = 0.0
        do j=1,supermesh_field_shape%loc
          val = val + eval_shape(supermesh_field_shape, j, local_coords) * alpha_sf_ele_value(j)
        end do
        call set(alpha_sf_on_supermesh, node_C, val)

      end do ! end of looping over the nodes of element 'ele_C' of the supermesh

    end do ! end of looping over elements of the supermesh

    ! supermesh_positions_remapped no longer needed:
    call deallocate(supermesh_positions_remapped)
    call deallocate(supermesh_field_mesh)

  end subroutine compute_alpha_on_supermesh

  subroutine compute_rhs_galerkin_projection_oneway_fsi(ele_B, rhs_alpha_sf, positionsA, positionsB, alpha_sf_on_supermesh, &
                                                        supermesh, inversion_matrices_A)

    ! FSI: fluid-solid interactions
    integer, intent(in) :: ele_B
    type(scalar_field), intent(inout) :: rhs_alpha_sf
    type(vector_field), intent(in) :: positionsA, positionsB
    type(scalar_field), intent(in) :: alpha_sf_on_supermesh
    type(vector_field), intent(in) :: supermesh
    real, dimension(:, :, :), intent(in) :: inversion_matrices_A

    integer :: ele_A, ele_C

    do ele_C=1,ele_count(supermesh)
      ele_A = ele_region_id(supermesh, ele_C) ! get the parent fluid element
      call compute_rhs_galerkin_projection_oneway_fsi_ele(rhs_alpha_sf, positionsA, positionsB, alpha_sf_on_supermesh, &
                                                          supermesh, inversion_matrices_A(:, :, ele_A), ele_A, ele_B, ele_C)
    end do
  end subroutine compute_rhs_galerkin_projection_oneway_fsi

  subroutine compute_rhs_galerkin_projection_oneway_fsi_ele(rhs_alpha_sf, positionsA, positionsB, alpha_sf_on_supermesh, &
                                                            supermesh, inversion_matrix_A, ele_A, ele_B, ele_C)
    integer, intent(in) :: ele_A, ele_B, ele_C
    type(scalar_field), intent(inout) :: rhs_alpha_sf
    type(vector_field), intent(in) :: positionsA, positionsB
    type(scalar_field), intent(in) :: alpha_sf_on_supermesh
    type(vector_field), intent(in) :: supermesh
    real, dimension(:, :), intent(in) :: inversion_matrix_A

    real, dimension(ele_ngi(supermesh, ele_C)) :: detwei_C
    real, dimension(positionsA%dim+1, ele_ngi(supermesh, ele_C)) :: local_coord_at_quad_fluid
    real, dimension(positionsB%dim+1, ele_ngi(supermesh, ele_C)) :: local_coord_at_quad_solid
    real, dimension(supermesh%dim+1, ele_ngi(supermesh, ele_C)) :: global_coord_at_quad

    real, dimension(ele_loc(rhs_alpha_sf, ele_A), ele_ngi(supermesh, ele_C)) :: basis_at_quad_fluid
    real, dimension(ele_loc(rhs_alpha_sf, ele_C), ele_ngi(supermesh, ele_C)) :: basis_at_quad_supermesh

    real, dimension(ele_loc(rhs_alpha_sf, ele_A), ele_loc(alpha_sf_on_supermesh, ele_C)) :: mat_fluid

    real, dimension(ele_loc(supermesh, ele_C)) :: supermesh_alpha_sf_val

    integer :: j, k, l, dim

    ! Compute the local coordinates of the fluid and solid meshes.
    global_coord_at_quad(1:supermesh%dim, :) = ele_val_at_quad(supermesh, ele_C)
    global_coord_at_quad(supermesh%dim+1, :) = 1.0
    local_coord_at_quad_fluid = matmul(inversion_matrix_A, global_coord_at_quad)

    ! Compute the basis functions of the fluid, solid and supermesh force fields at these local coordinates.
    basis_at_quad_fluid = basis_at_quad(ele_shape(rhs_alpha_sf, ele_A), local_coord_at_quad_fluid)
    basis_at_quad_supermesh = basis_at_quad(ele_shape(alpha_sf_on_supermesh, ele_C), alpha_sf_on_supermesh%mesh%shape%n)

    ! We need to integrate over the supermesh element, so get detwei_C
    call transform_to_physical(supermesh, ele_C, detwei=detwei_C)

    ! Compute the little mixed mass matrices
    mat_fluid = 0.0
    !mat_solid = 0.0
    do j=1,ele_ngi(supermesh, ele_C)
      forall (k=1:ele_loc(rhs_alpha_sf, ele_A),l=1:ele_loc(alpha_sf_on_supermesh, ele_C))
        mat_fluid(k, l) = mat_fluid(k, l) + detwei_C(j) * basis_at_quad_fluid(k, j) * basis_at_quad_supermesh(l, j)
      end forall
    end do

    ! Now compute the rhs contribution for alpha:
    supermesh_alpha_sf_val = ele_val(alpha_sf_on_supermesh, ele_C)
    call addto(rhs_alpha_sf, ele_nodes(rhs_alpha_sf, ele_A), matmul(mat_fluid, supermesh_alpha_sf_val))

  end subroutine compute_rhs_galerkin_projection_oneway_fsi_ele

  subroutine assemble_mass_matrix(mass, field, positions, ele)
    type(csr_matrix), intent(inout) :: mass
    type(vector_field), intent(in) :: field
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: ele

    real, dimension(ele_ngi(field, ele)) :: detwei
    real, dimension(ele_loc(field, ele), ele_loc(field, ele)) :: little_mass

    call transform_to_physical(positions, ele, detwei=detwei)
    little_mass = shape_shape(ele_shape(field, ele), ele_shape(field, ele), detwei)
    call addto(mass, ele_nodes(field, ele), ele_nodes(field, ele), little_mass)
  end subroutine assemble_mass_matrix

  function basis_at_quad(shape, local_coords) result(basis)
    type(element_type), intent(in) :: shape
    real, dimension(:, :), intent(in) :: local_coords
    real, dimension(shape%loc, size(local_coords, 2)) :: basis
    integer :: loc, gi

    if (shape%degree == 0) then
      basis = 1.0
    elseif (shape%degree == 1) then
      basis = local_coords
    else
      do loc=1,shape%loc
        do gi=1,size(local_coords, 2)
          basis(loc, gi) = eval_shape(shape, loc, local_coords(:, gi))
        end do
      end do
    end if
  end function basis_at_quad

  subroutine bound_projection(field, rhs_field, mass_matrix, mass_matrix_lumped, inverse_mass_matrix_lumped, positions)  
    ! This subroutine bounds the projection of a scalar field

    type(scalar_field), intent(inout) :: field, mass_matrix_lumped
    type(scalar_field), intent(in) :: rhs_field, inverse_mass_matrix_lumped
    type(csr_matrix), intent(in) :: mass_matrix
    type(vector_field), intent(in) :: positions

    type(scalar_field) :: bounded_soln, max_bound, min_bound
    real :: upper_bound, lower_bound
    type(csr_sparsity), pointer :: nnlist
    integer :: node_B
    integer, dimension(:), pointer :: patch

    upper_bound = huge(0.0)*epsilon(0.0) ! default value
    lower_bound = -huge(0.0)*epsilon(0.0) ! default value

    ewrite(2,*) "BOUNDS==========================================================================="
    ewrite(2,*) "upper_bound = ", upper_bound
    ewrite(2,*) "lower_bound = ", lower_bound
    ewrite(2,*) "BOUNDS==========================================================================="

    ! Enable the bounds to vary locally
    call allocate(max_bound, field%mesh, "MaxBound")
    call allocate(min_bound, field%mesh, "MinBound")
    call set(max_bound, upper_bound)
    call set(min_bound, lower_bound)
    call allocate(bounded_soln, field%mesh, "BoundedSolution")

    ewrite(2,*) "***********BEFORE********************"
    ewrite_minmax(field%val(:))

    call set(bounded_soln, rhs_field)
    call scale(bounded_soln, inverse_mass_matrix_lumped)
    call halo_update(bounded_soln)

    nnlist => extract_nnlist(field)

    do node_B = 1,node_count(field%mesh)
      patch => row_m_ptr(nnlist, node_B)
      call set(max_bound, node_B, max(min(maxval(bounded_soln%val(patch)), node_val(max_bound, node_B)), lower_bound))
      call set(min_bound, node_B, max(min(minval(bounded_soln%val(patch)), node_val(max_bound, node_B)), lower_bound))
    end do

    call halo_update(max_bound)
    ewrite_minmax(max_bound)
    call halo_update(min_bound)
    ewrite_minmax(min_bound)

    call bound_field(field, max_bound, min_bound, mass_matrix, mass_matrix_lumped, &
                        & inverse_mass_matrix_lumped, bounded_soln, positions)

    ewrite(2,*) "***********AFTER********************"
    ewrite_minmax(field%val(:))

    call deallocate(max_bound)
    call deallocate(min_bound)
    call deallocate(bounded_soln)

   end subroutine bound_projection

   subroutine fsi_one_way_grandy_interpolation(fluid_position, solid_position, alpha)

       type(vector_field), pointer, intent(in) :: fluid_position
       type(vector_field), intent(in) :: solid_position
       type(scalar_field), intent(inout) :: alpha

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

       ! Set the input of the RTREE finder as the coordinates of the fluids mesh:
       call rtree_intersection_finder_set_input(fluid_position)

       ! For all elements of the solids mesh:
       do ele_B = 1, ele_count(solid_position)

          ! Via RTREE, find intersection of solid element 'ele_B' with the
          ! input mesh (fluid coordinate mesh 'positions'):
          call rtree_intersection_finder_find(solid_position, ele_B)
          ! Fetch output, the number of intersections of this solid element:
          call rtree_intersection_finder_query_output(nintersections)

          if (fluid_position%dim == 3) then
             ! Get (solid) element value (coordinate), form tetrahedra:
             tet_B%v = ele_val(solid_position, ele_B)
          else
             call intersector_set_dimension(fluid_position%dim)
          end if

          ! 1st inner-loop
          ! For all intersections of solid element ele_B with fluid mesh:
          do j = 1, nintersections
             ! Get the donor (fluid) element which intersects with ele_B
             call rtree_intersection_finder_get_output(ele_A, j)
             ! If 3D
             if (fluid_position%dim == 3) then
                ! Get the global coordinates of fluid element ele_A:
                if (ele_loc(fluid_position, ele_A)==4) then
                   ! if fluid element is a tetrahedra
                   allocate(planes_A(4))
                   tet_A%v = ele_val(fluid_position, ele_A)
                   planes_A = get_planes(tet_A)
                else
                   ! if fluid element is not a tetrahedra, 
                   ! assumed to be a hexahedra:
                   allocate(planes_A(6))
                   planes_A = get_planes(fluid_position, ele_A)
                end if
                ! Get the coordinates of nodes of the intersection
                ! between of both elements, store in intersection
                call intersect_tets(tet_B, planes_A, &
                     ele_shape(solid_position, ele_B), &
                     stat=stat, output=intersection)
                deallocate(planes_A)
             else ! 2D
                allocate(pos_A(fluid_position%dim, ele_loc(fluid_position, ele_A)))
                pos_A = ele_val(fluid_position, ele_A)
                intersection = intersect_elements(solid_position, ele_B, pos_A, ele_shape(solid_position, ele_B))
                deallocate(pos_A)
                stat = 0
             end if ! end of dim==3

             ! No intersection, cycle:
             if (stat == 1) cycle

             ! Compute intersection volume:
             vol = 0.0
             do ele_C = 1, ele_count(intersection)
                vol = vol + abs(simplex_volume(intersection, ele_C))
             end do

             ! Compute the volume of the fluid element:
             allocate(detwei(ele_ngi(fluid_position, ele_A)))
             call transform_to_physical(fluid_position, ele_A, detwei=detwei)
             ele_A_vol = sum(detwei)
             deallocate(detwei)

             ! Compute the volume fraction:
             ! ele_A_nodes: pointer to global node numbers of
             ! fluid element ele_A of the coordinate mesh
             ele_A_nodes => ele_nodes(fluid_position, ele_A)
             do k = 1, size(ele_A_nodes)
                ! Volume fraction by grandy projection
                call addto(alpha, ele_A_nodes(k), vol/ele_A_vol)
             end do
             call deallocate(intersection)
          end do
       end do

       call finalise_tet_intersector
       call rtree_intersection_finder_reset(ntests)

       ! Bound field:
       call bound_scalar_field(alpha, 0.0, 1.0)

    end subroutine fsi_one_way_grandy_interpolation

    subroutine bound_scalar_field(sfield, minvalue, maxvalue)
        type(scalar_field), intent(inout) :: sfield
        real, intent(in) :: minvalue, maxvalue
        integer :: i

        do i = 1, node_count(sfield)
           call set(sfield, i, max(minvalue, min(maxvalue, node_val(sfield, i))))
        end do

    end subroutine bound_scalar_field

end module fsi_projections

