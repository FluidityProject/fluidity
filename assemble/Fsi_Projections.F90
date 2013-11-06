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
  use fefields, only: compute_lumped_mass, project_field
  use tetrahedron_intersection_module

  implicit none

  interface bound_projection
    module procedure bound_projection_scalars, bound_projection_vectors
  end interface

  public :: fsi_one_way_galerkin_projection, fsi_one_way_grandy_interpolation, fsi_get_interface

  contains

  subroutine fsi_one_way_galerkin_projection(fieldF, positionsF, positionsS, alpha_sf, solver_option_path, solid_velocity_on_solid, solid_velocity_sf)
    !! Return the volume fraction scalar field by projecting unity from the supermesh to the fluid mesh
    !! and, iff given, the solid velocity from the solid (via the supermesh) to the fluid mesh.
    !! Since positionsF and positionsS are different, we need to supermesh!
    !! positionsF and positionsS are the coordinate fields of the fluid and solid mesh respectively

    ! F stands for fluid, S for solid, C for the supermesh in this case:
    type(vector_field), intent(inout) :: fieldF
    type(vector_field), intent(inout) :: positionsF, positionsS
    type(scalar_field), intent(inout) :: alpha_sf
    character(len=OPTION_PATH_LEN), intent(in) :: solver_option_path
    type(vector_field), intent(in), optional :: solid_velocity_on_solid
    type(vector_field), intent(inout), optional :: solid_velocity_sf

    type(ilist), dimension(ele_count(positionsS)) :: map_SF
    integer :: ele_F, ele_S

    type(quadrature_type) :: supermesh_quad
    type(element_type) :: supermesh_positions_shape, supermesh_field_shape

    type(vector_field) :: supermesh
    type(mesh_type) :: supermesh_field_mesh

    ! Scalar fields for alpha:
    type(scalar_field) :: alpha_sf_on_solid
    type(scalar_field) :: alpha_sf_on_supermesh
    ! Vector fields for solid velocity:
    type(vector_field) :: solid_velocity_on_supermesh

    type(csr_sparsity) :: sparsity_fluid
    type(csr_matrix) :: mass_matrix_fluid
    type(scalar_field) :: rhs_alpha_sf
    type(vector_field) :: rhs_fluid

    real, dimension(ele_loc(positionsS, 1), ele_loc(positionsS, 1)) :: inversion_matrix_S
    real, dimension(ele_loc(positionsS, 1), ele_loc(positionsS, 1), ele_count(positionsF)) :: inversion_matrices_F
    integer :: dim, max_degree

    character(len=OPTION_PATH_LEN) :: tmp

    integer :: stat

    ! As the projections are bounded, we need to define lumped versions of the mass matrices:
    !type(scalar_field) :: lumped_mass_matrix_solid, lumped_mass_matrix_fluid
    !type(scalar_field) :: lumped_inverse_mass_matrix_solid, lumped_inverse_mass_matrix_fluid
    type(scalar_field) :: lumped_mass_matrix_fluid
    type(scalar_field) :: lumped_inverse_mass_matrix_fluid


    ewrite(2,*) "inside fsi_one_way_galerkin_projection"

    dim = mesh_dim(positionsS)
    call intersector_set_dimension(dim)

    ! Define polynomial degree and quadrature of the supermesh:
    max_degree = max(element_degree(fieldF, 1), element_degree(positionsS, 1))
    supermesh_quad = make_quadrature(vertices=ele_loc(positionsS, 1), dim=dim, degree=2*max_degree)
    ! Define the shape function of the supermesh:
    supermesh_positions_shape = make_element_shape(vertices=ele_loc(positionsS, 1), dim=dim, degree=1, quad=supermesh_quad)
    supermesh_field_shape = make_element_shape(vertices=ele_loc(positionsS, 1), dim=dim, degree=max_degree, quad=supermesh_quad)
    ! For each element in S get the list of intersecting elements in F
    map_SF = rtree_intersection_finder(positionsS, positionsF)

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
    ! and for the solid velocity, IFF corresponding files were given:
    if (present(solid_velocity_on_solid) .and. present(solid_velocity_sf)) then
        call allocate(rhs_fluid, fieldF%dim, fieldF%mesh, "FluidRHS")
        call zero(rhs_fluid)
    end if 

    ! Looping over all the elements of (solid) mesh S:
    do ele_S=1,ele_count(positionsS)

      ! Only continue if the element ele_B has an intersection with the other mesh:
      if (.not. map_SF(ele_S)%length > 0) then
        cycle ! no intersection found for this ele_B
      end if

      ! get the matrix for the coordinates of (solid) mesh S:
      call local_coords_matrix(positionsS, ele_S, inversion_matrix_S)
      ! Construct the supermesh associated with ele_S. (For advancing front algorithm, uncomment the following line:)
      call construct_supermesh(positionsS, ele_S, positionsF, map_SF(ele_S), supermesh_positions_shape, supermesh, stat=stat)
      ! if no intersection for proc x was found, no supermesh was created, then stat/=0
      if (stat /= 0) then ! should not happen, as we catch that event above, but doesn't hurt double checking
        cycle
      end if

      ! At this point, a portion of the supermesh is constructed, which refines the dimensional space
      ! of the intersection of ele_S with mesh F.
      ! This portion of the supermesh is stored in 'supermesh'

      ! 1st step: Compute the project unity from the solid mesh to the supermesh:

      ! Allocate field for alpha on the supermesh:
      supermesh_field_mesh = make_mesh(supermesh%mesh, supermesh_field_shape, -1, "SupermeshFieldMesh")
      call allocate(alpha_sf_on_supermesh, supermesh_field_mesh, "AlphaSFOnSupermesh")
      call zero(alpha_sf_on_supermesh)
      ! Same for solid velocity (IFF it was given):
      if (present(solid_velocity_on_solid) .and. present(solid_velocity_sf)) then
          call allocate(solid_velocity_on_supermesh, solid_velocity_on_solid%dim, supermesh_field_mesh, "SolidVelocityOnSupermesh")
          call zero(solid_velocity_on_supermesh)
      end if

      ! Get alpha for this portion of the supermesh:
      if (present(solid_velocity_on_solid) .and. present(solid_velocity_sf)) then
        call compute_alpha_on_supermesh(supermesh_field_shape, positionsS, ele_S, ele_val(alpha_sf_on_solid, ele_S), &
                                        supermesh, inversion_matrix_S, alpha_sf_on_supermesh, &
                                        ele_val(solid_velocity_on_solid, ele_S), solid_velocity_on_supermesh)
      else
        call compute_alpha_on_supermesh(supermesh_field_shape, positionsS, ele_S, ele_val(alpha_sf_on_solid, ele_S), &
                                        supermesh, inversion_matrix_S, alpha_sf_on_supermesh)
      end if

      ! 2nd step: Compute the RHS of the Galerkin projection:
      ! Elemental Values of alpha of this part of the supermesh are now computed.
      ! Compute the RHS of 
      ! M_f * alpha_f = M_{f,sup} * alpha_{sup}
      ! and afterwards, the supermesh can be deleted as it is no longer needed,
      ! remember, (M_{f,sup} * alpha_{sup}) lives on the fluid mesh, not the supermesh.
      if (present(solid_velocity_on_solid) .and. present(solid_velocity_sf)) then
        call compute_rhs_galerkin_projection_oneway_fsi(ele_S, positionsF, positionsS, supermesh, inversion_matrices_F, &
                                                        rhs_alpha_sf, alpha_sf_on_supermesh, &
                                                        rhs_fluid, solid_velocity_on_supermesh)
      else
        call compute_rhs_galerkin_projection_oneway_fsi(ele_S, positionsF, positionsS, supermesh, inversion_matrices_F, &
                                                        rhs_alpha_sf, alpha_sf_on_supermesh)
      end if

      ! Portion of supermesh no longer needed
      call deallocate(supermesh)
      call deallocate(supermesh_field_mesh)
      call deallocate(alpha_sf_on_supermesh)
      if (present(solid_velocity_on_solid) .and. present(solid_velocity_sf)) then
        call deallocate(solid_velocity_on_supermesh)
      end if
    end do

    call deallocate(supermesh_quad)
    call deallocate(supermesh_positions_shape)
    call deallocate(supermesh_field_shape)
    do ele_S=1,ele_count(positionsS)
      call deallocate(map_SF(ele_S))
    end do

    ! 3rd step: Project alpha from the supermesh to the fluid and solid mesh:
    ! loop over fluid and solid elements and solve the last equation, which will
    ! project alpha from the supermesh to the fluid mesh

    ! Project alpha to the fluid mesh:
    call petsc_solve(alpha_sf, mass_matrix_fluid, rhs_alpha_sf, option_path=trim(solver_option_path))

    ! Same for solid velocity, IFF it was given:
    if (present(solid_velocity_on_solid) .and. present(solid_velocity_sf)) then
        ! Project solid velocity to the fluid mesh:
        call petsc_solve(solid_velocity_sf, mass_matrix_fluid, rhs_fluid, option_path=trim(solver_option_path))
    end if

    ! Get options from option tree if projection is bounded or not:
    if (have_option(trim(solver_option_path)//'/bounded[0]')) then
        call bound_projection(alpha_sf, rhs_alpha_sf, mass_matrix_fluid, &
                              lumped_mass_matrix_fluid, lumped_inverse_mass_matrix_fluid, &
                              positionsF)
        ! IFF solid velocity was given, bound that vector field as well:
        if (present(solid_velocity_on_solid) .and. present(solid_velocity_sf)) then
            call bound_projection(solid_velocity_sf, rhs_fluid, mass_matrix_fluid, &
                                  lumped_mass_matrix_fluid, lumped_inverse_mass_matrix_fluid, &
                                  positionsF)
        end if
    end if

    call deallocate(mass_matrix_fluid)
    call deallocate(alpha_sf_on_solid)
    call deallocate(rhs_alpha_sf)
    if (present(solid_velocity_on_solid) .and. present(solid_velocity_sf)) then
      call deallocate(rhs_fluid)
    end if
    call deallocate(lumped_mass_matrix_fluid)
    call deallocate(lumped_inverse_mass_matrix_fluid)

    ewrite(2,*) "leaving fsi_one_way_galerkin_projection"

  end subroutine fsi_one_way_galerkin_projection

  subroutine compute_alpha_on_supermesh(supermesh_field_shape, solid_positions, ele_B, alpha_sf_ele_value, &
                                        supermesh, inversion_matrix_B, alpha_sf_on_supermesh, solid_velocity_ele_value, solid_velocity_on_supermesh)
    ! F stands for fluid, S for solid, C for the supermesh in this case:
    type(vector_field), intent(in) :: solid_positions, supermesh
    type(element_type), intent(in), target :: supermesh_field_shape
    integer, intent(in) :: ele_B
    real, dimension(:), intent(in) :: alpha_sf_ele_value
    real, dimension(:, :), intent(in) :: inversion_matrix_B

    type(scalar_field), intent(inout) :: alpha_sf_on_supermesh
    real, dimension(:, :), intent(in), optional :: solid_velocity_ele_value
    type(vector_field), intent(inout), optional :: solid_velocity_on_supermesh

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
        
        ! Now the same for the solid velocity (IFF it is given):
        if (present(solid_velocity_ele_value) .and. present(solid_velocity_on_supermesh)) then
          do k=1,solid_positions%dim
            val = 0.0
            do j=1,supermesh_field_shape%loc
              val = val + eval_shape(supermesh_field_shape, j, local_coords) * solid_velocity_ele_value(k, j)
            end do
            call set(solid_velocity_on_supermesh, k, node_C, val)
          end do
        end if

      end do ! end of looping over the nodes of element 'ele_C' of the supermesh

    end do ! end of looping over elements of the supermesh

    ! supermesh_positions_remapped no longer needed:
    call deallocate(supermesh_positions_remapped)
    call deallocate(supermesh_field_mesh)

  end subroutine compute_alpha_on_supermesh

  subroutine compute_rhs_galerkin_projection_oneway_fsi(ele_B, positionsA, positionsB, supermesh, inversion_matrices_A, &
                                                         rhs_alpha_sf, alpha_sf_on_supermesh, &
                                                         rhs_fluid, solid_velocity_on_supermesh)

    ! FSI: fluid-solid interactions
    integer, intent(in) :: ele_B
    type(vector_field), intent(in) :: positionsA, positionsB
    type(vector_field), intent(in) :: supermesh
    real, dimension(:, :, :), intent(in) :: inversion_matrices_A
    type(scalar_field), intent(inout) :: rhs_alpha_sf
    type(scalar_field), intent(in) :: alpha_sf_on_supermesh
    type(vector_field), intent(inout), optional :: rhs_fluid
    type(vector_field), intent(in), optional :: solid_velocity_on_supermesh

    integer :: ele_A, ele_C

    if (present(rhs_fluid) .and. present(solid_velocity_on_supermesh)) then
      do ele_C=1,ele_count(supermesh)
        ele_A = ele_region_id(supermesh, ele_C) ! get the parent fluid element
        call compute_rhs_galerkin_projection_oneway_fsi_ele(positionsA, positionsB, supermesh, &
                                                            inversion_matrices_A(:, :, ele_A), ele_A, ele_B, ele_C, &
                                                            rhs_alpha_sf, alpha_sf_on_supermesh, &
                                                            rhs_fluid, solid_velocity_on_supermesh)
      end do
    else ! without solid velocity projection:
      do ele_C=1,ele_count(supermesh)
        ele_A = ele_region_id(supermesh, ele_C) ! get the parent fluid element
        call compute_rhs_galerkin_projection_oneway_fsi_ele(positionsA, positionsB, supermesh, &
                                                            inversion_matrices_A(:, :, ele_A), ele_A, ele_B, ele_C, &
                                                            rhs_alpha_sf, alpha_sf_on_supermesh)
      end do
    end if
  end subroutine compute_rhs_galerkin_projection_oneway_fsi

  subroutine compute_rhs_galerkin_projection_oneway_fsi_ele(positionsA, positionsB, supermesh, &
                                                             inversion_matrix_A, ele_A, ele_B, ele_C, &
                                                             rhs_alpha_sf, alpha_sf_on_supermesh, &
                                                             rhs_fluid, solid_velocity_on_supermesh)

    type(vector_field), intent(in) :: positionsA, positionsB
    type(vector_field), intent(in) :: supermesh
    real, dimension(:, :), intent(in) :: inversion_matrix_A
    integer, intent(in) :: ele_A, ele_B, ele_C
    type(scalar_field), intent(inout) :: rhs_alpha_sf
    type(scalar_field), intent(in) :: alpha_sf_on_supermesh
    type(vector_field), intent(inout), optional :: rhs_fluid
    type(vector_field), intent(in), optional :: solid_velocity_on_supermesh

    real, dimension(ele_ngi(supermesh, ele_C)) :: detwei_C
    real, dimension(positionsA%dim+1, ele_ngi(supermesh, ele_C)) :: local_coord_at_quad_fluid
    real, dimension(positionsB%dim+1, ele_ngi(supermesh, ele_C)) :: local_coord_at_quad_solid
    real, dimension(supermesh%dim+1, ele_ngi(supermesh, ele_C)) :: global_coord_at_quad

    real, dimension(ele_loc(rhs_alpha_sf, ele_A), ele_ngi(supermesh, ele_C)) :: basis_at_quad_fluid
    real, dimension(ele_loc(rhs_alpha_sf, ele_C), ele_ngi(supermesh, ele_C)) :: basis_at_quad_supermesh

    real, dimension(ele_loc(rhs_alpha_sf, ele_A), ele_loc(alpha_sf_on_supermesh, ele_C)) :: mat_fluid

    real, dimension(ele_loc(supermesh, ele_C)) :: supermesh_alpha_sf_val
    real, dimension(positionsB%dim, ele_loc(supermesh, ele_C)) :: supermesh_solid_velocity_val

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

    ! Now compute the rhs contribution for a vector field:
    if (present(rhs_fluid) .and. present(solid_velocity_on_supermesh)) then
      supermesh_solid_velocity_val = ele_val(solid_velocity_on_supermesh, ele_C)
      do dim=1,solid_velocity_on_supermesh%dim
        call addto(rhs_fluid, dim, ele_nodes(rhs_fluid, ele_A), matmul(mat_fluid, supermesh_solid_velocity_val(dim, :)))
      end do
    end if
    ! Now compute the rhs contribution for alpha (scalar field):
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

  subroutine bound_projection_vectors(vfield, rhs_vfield, mass_matrix, lumped_mass_matrix, lumped_inverse_mass_matrix, positions)
    ! Bounding a vector field
    type(vector_field), intent(inout) :: vfield
    type(vector_field), intent(in) :: rhs_vfield
    type(csr_matrix), intent(in) :: mass_matrix
    type(scalar_field), intent(inout) :: lumped_mass_matrix
    type(scalar_field), intent(in) :: lumped_inverse_mass_matrix
    type(vector_field), intent(in) :: positions

    type(scalar_field) :: sfield, rhs_sfield
    integer :: d

    ewrite(2,*) "Inside bound_projection_vectors"

    call allocate(sfield, vfield%mesh, "VectorComponent")
    call zero(sfield)
    call allocate(rhs_sfield, vfield%mesh, "RHSVectorComponent")
    call zero(rhs_sfield)

    ! For dimension of given vector field, give the interface 'bound_projection' each vector component at a time
    do d=1, positions%dim
      sfield%val(:) = vfield%val(d,:)
      rhs_sfield%val(:) = rhs_vfield%val(d,:)
      call bound_projection(sfield, rhs_sfield, mass_matrix, lumped_mass_matrix, lumped_inverse_mass_matrix, positions)
      vfield%val(d,:) = sfield%val(:)
    end do

    call deallocate(sfield)
    call deallocate(rhs_sfield)

    ewrite(2,*) "Leaving bound_projection_vectors"

  end subroutine bound_projection_vectors

  subroutine bound_projection_scalars(sfield, rhs_field, mass_matrix, lumped_mass_matrix, inverse_mass_matrix_lumped, positions)
    ! Bounding a field (only scalar) after an inter-mesh projection.
    ! Similar to subroutine 'interpolation_galerkin_scalars' in module 'conservative_interpolation_module'
    type(scalar_field), intent(inout) :: sfield, lumped_mass_matrix
    type(scalar_field), intent(in) :: rhs_field, inverse_mass_matrix_lumped
    type(csr_matrix), intent(in) :: mass_matrix
    type(vector_field), intent(in) :: positions

    type(scalar_field) :: bounded_soln, max_bound, min_bound
    real :: upper_bound, lower_bound
    type(csr_sparsity), pointer :: nnlist
    integer :: node_B
    integer, dimension(:), pointer :: patch
    
    ewrite(2,*) "Inside bound_projection_scalars"

    ! Step 0: Computing the bounds,
    ! here we hard-code the default bounds:
    upper_bound = huge(0.0)*epsilon(0.0)
    lower_bound = -huge(0.0)*epsilon(0.0)

    ! Set the default bounds to the whole mesh:
    call allocate(max_bound, sfield%mesh, "MaxBound")
    call allocate(min_bound, sfield%mesh, "MinBound")
    call set(max_bound, upper_bound)
    call set(min_bound, lower_bound)
    call allocate(bounded_soln, sfield%mesh, "BoundedSolution")

    ! Preparing the solution:
    call set(bounded_soln, rhs_field)
    call scale(bounded_soln, inverse_mass_matrix_lumped)
    call halo_update(bounded_soln)

    nnlist => extract_nnlist(sfield)

    do node_B = 1,node_count(sfield%mesh)
      patch => row_m_ptr(nnlist, node_B)
      call set(max_bound, node_B, max(min(maxval(bounded_soln%val(patch)), &
                                                   node_val(max_bound, node_B)), &
                                                   lower_bound))
      call set(min_bound, node_B, max(min(minval(bounded_soln%val(patch)), &
                                                   node_val(max_bound, node_B)), &
                                                   lower_bound))
    end do

    call halo_update(max_bound)
    call halo_update(min_bound)

    call bound_field(sfield, max_bound, min_bound, mass_matrix, lumped_mass_matrix, &
                     inverse_mass_matrix_lumped, bounded_soln, positions)

    call deallocate(max_bound)
    call deallocate(min_bound)
    call deallocate(bounded_soln)

    ewrite(2,*) "Leaving bound_projection_scalars"

  end subroutine bound_projection_scalars

  subroutine fsi_one_way_grandy_interpolation(fluid_position, solid_position, alpha)

     type(vector_field), pointer, intent(inout) :: fluid_position
     type(vector_field), intent(in) :: solid_position
     type(scalar_field), intent(inout) :: alpha
     type(scalar_field) :: alpha_coordinatemesh

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

     call allocate(alpha_coordinatemesh, fluid_position%mesh, 'SolidConcentrationCoordinateMesh')
     call zero(alpha_coordinatemesh)

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
              call addto(alpha_coordinatemesh, ele_A_nodes(k), vol/ele_A_vol)
           end do
           call deallocate(intersection)
        end do
     end do

     call finalise_tet_intersector
     call rtree_intersection_finder_reset(ntests)

     ! Bound field:
     call bound_scalar_field(alpha_coordinatemesh, 0.0, 1.0)

     ! Remap to alpha:
     call remap_field(alpha_coordinatemesh, alpha)
     call deallocate(alpha_coordinatemesh)

  end subroutine fsi_one_way_grandy_interpolation

  subroutine bound_scalar_field(sfield, minvalue, maxvalue)
      type(scalar_field), intent(inout) :: sfield
      real, intent(in) :: minvalue, maxvalue
      integer :: i, ele
      integer, dimension(:), pointer :: nodes

      do ele = 1, ele_count(sfield%mesh)
          nodes => ele_nodes(sfield%mesh, ele)
          do i = 1, size(nodes)
              call set(sfield, nodes(i), max(minvalue, min(maxvalue, node_val(sfield, nodes(i)))))
          end do
      end do

  end subroutine bound_scalar_field

  subroutine fsi_get_interface(fsi_interface, positionsF, positionsS)
  !! This subroutine returns a scalar field, that is 1 only at the elements which intersect
  !! with the surface elements of a solid mesh
      type(scalar_field), intent(inout) :: fsi_interface
      type(vector_field), pointer, intent(in) :: positionsF, positionsS

      integer :: ele_F, ele_S
      integer :: num_intersections
      integer :: i, j, n, dim, node
      
      integer, dimension(:), pointer :: solid_faces
      integer, dimension(:), pointer :: nodes
      integer, dimension(face_loc(positionsS, 1)) :: solid_face_nodes
      integer :: face_number

      ewrite(2,*) "inside fsi_get_interface"

      ! Set the dimension for the intersection finder:
      dim = mesh_dim(positionsS)
      call intersector_set_dimension(dim) 
      ! Use rtree intersection finder because it is more robust, e.g.
      ! it does not fail when running in parallel and no solid element is found in the 
      ! composition of the fluid mesh:
      ! Set input for rtree intersection finder:
      call rtree_intersection_finder_set_input(positionsF)

      do ele_S=1,ele_count(positionsS)
          solid_faces => ele_faces(positionsS, ele_S)
          do i=1, ele_face_count(positionsS,ele_S)
              !face_number = 0.0
              face_number = face_neigh(positionsS, solid_faces(i))
              if (face_number == solid_faces(i)) then
                  solid_face_nodes = 0.0
                  ! This means, the face 'face_number' has no neighbour face, meaning
                  ! this face is at the solid's boundary, 
                  ! so let's get intersecting fluid elements for this solid element:

                  ! Via rtree, find intersection of solid element with fluid elements:
                  call rtree_intersection_finder_find(positionsS, ele_S)
                  ! Fetch output, the number of intersections
                  call rtree_intersection_finder_query_output(num_intersections)
                  ! Set scalar field for intersection with the solid boundary
                  if (num_intersections .ne. 0) then
                      do j=1, num_intersections
                          ! Get the donor (fluid) element which intersects with ele_S
                          call rtree_intersection_finder_get_output(ele_F, j)
                          nodes => ele_nodes(positionsF, ele_F)
                          do n = 1, size(nodes)
                              node = nodes(n)
                              call set(fsi_interface, node, 1.0)
                          end do
                      end do
                  end if
              end if
          end do
      end do

      ewrite(2,*) "leaving fsi_get_interface"

  end subroutine fsi_get_interface

end module fsi_projections

