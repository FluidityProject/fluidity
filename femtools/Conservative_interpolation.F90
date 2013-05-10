#include "fdebug.h"

#define INLINE_MATMUL

module conservative_interpolation_module

  use FLDebug
  use quadrature
  use elements
  use fields
  use sparse_tools
  use supermesh_construction
  use futils
  use transform_elements
  use meshdiagnostics
  use sparsity_patterns
  use vector_tools
  use tensors
  use fetools
  use sparse_tools
  use interpolation_module
  use solvers
  use adjacency_lists
  use vtk_interfaces
  use unittest_tools
  use spud
  use global_parameters, only : FIELD_NAME_LEN, OPTION_PATH_LEN
  use intersection_finder_module
  use linked_lists
  use sparse_matrices_fields
  use bound_field_module
  use halos
  use diagnostic_fields
  use tetrahedron_intersection_module
  use boundary_conditions
  use data_structures
  implicit none

  interface interpolation_galerkin
    module procedure interpolation_galerkin_scalars, interpolation_galerkin_single_state, &
                     interpolation_galerkin_multiple_states
  end interface

  interface grandy_projection
    module procedure grandy_projection_scalars, grandy_projection_multiple_states
  end interface

  public :: interpolation_galerkin, grandy_projection
  
  private
  
#ifdef DUMP_SUPERMESH_INTERSECTIONS
  integer :: dump_idx
#endif

  contains

  subroutine galerkin_projection_inner_loop(ele_B, little_mass_matrix, detJ, local_rhs, conservation_tolerance, stat, &
                                            field_counts, old_fields, old_position, new_fields, new_position, &
                                            map_BA, inversion_matrices_A, supermesh_shape)
  
    integer, intent(in) :: ele_B
    real, dimension(:,:,:), intent(inout) :: little_mass_matrix
    real, dimension(:), intent(out) :: detJ
    real, dimension(:,:,:), intent(inout) :: local_rhs
    real, intent(in) :: conservation_tolerance
    integer, intent(out) :: stat
    
    integer, dimension(:), intent(in) :: field_counts
    
    type(scalar_field), dimension(:,:), intent(in) :: old_fields
    type(vector_field), intent(in) :: old_position

    type(scalar_field), dimension(:,:), intent(inout) :: new_fields
    type(vector_field), intent(in) :: new_position

    type(ilist), dimension(:), intent(in) :: map_BA
    real, dimension(:, :, :) :: inversion_matrices_A
    type(element_type), intent(inout) :: supermesh_shape
   
    real, dimension(ele_loc(new_position, ele_B), ele_loc(new_position, ele_B)) :: inversion_matrix_B, inversion_matrix_A
    real, dimension(ele_ngi(new_position, ele_B)) :: detwei_B
    real, dimension(supermesh_shape%ngi) :: detwei_C
    type(inode), pointer :: llnode
    
    real :: vol_B, vols_C
    integer :: ele_A, ele_C, nloc, dim, j, k, l, loc, field, mesh, mesh_count
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

    ! llnode is looping over the intersecting elements for this ele_B
    llnode => map_BA(ele_B)%firstnode

    vol_B = sum(detwei_B)
    vols_C = 0.0

    ! loop over the intersecting elements
    do while (associated(llnode))
      ele_A = llnode%value
      ! but we only need that mapping for this ele_B now, so just compute it now
      if (dim == 3 .and. (intersector_exactness .eqv. .false.)) then
        tet_A%V = ele_val(old_position, ele_A)
        call intersect_tets(tet_A, planes_B, supermesh_shape, stat=lstat, output=intersection)
        if (lstat == 1) then
          llnode => llnode%next
          cycle
        end if
      else
        intersection = intersect_elements(old_position, ele_A, pos_B, supermesh_shape)
      end if

#ifdef DUMP_SUPERMESH_INTERSECTIONS
      if (ele_count(intersection) /= 0) then
        call vtk_write_fields("intersection", dump_idx, intersection, intersection%mesh)
        dump_idx = dump_idx + 1
      end if
#endif

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
      end do

      llnode => llnode%next
      call deallocate(intersection)
    end do

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

  end subroutine galerkin_projection_inner_loop

  subroutine interpolation_galerkin_scalars(old_fields_state, old_position, new_fields_state, new_position, map_BA, force_bounded)
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

    call intersector_set_dimension(dim)
    if (present(map_BA)) then
      lmap_BA => map_BA
    else
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

    do ele_B=1,ele_count(new_position)

       call galerkin_projection_inner_loop(ele_B, little_mass_matrix, detJ, local_rhs, conservation_tolerance, stat, &
                                           field_counts, old_fields, old_position, new_fields, new_position, &
                                           lmap_BA, inversion_matrices_A, supermesh_shape)

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

      if (.not. present(map_BA)) then
        do ele_B=1,ele_count(new_position)
          call deallocate(lmap_BA(ele_B))
        end do
        deallocate(lmap_BA)
      end if

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

    ewrite(1, *) "Exiting interpolation_galerkin_scalars"
    
  end subroutine interpolation_galerkin_scalars

  subroutine interpolation_galerkin_single_state(old_state, new_state, map_BA)
    type(state_type), intent(inout) :: old_state, new_state
    type(ilist), dimension(:), intent(in), optional :: map_BA

    type(state_type), dimension(1) :: old_states, new_states
    
    old_states = (/old_state/)
    new_states = (/new_state/)
    call interpolation_galerkin(old_states, new_states, map_BA=map_BA)
    old_state = old_states(1)
    new_state = new_states(1)
    
  end subroutine interpolation_galerkin_single_state

  subroutine interpolation_galerkin_multiple_states(old_states, new_states, map_BA)
    type(state_type), dimension(:), intent(inout) :: old_states, new_states
    type(ilist), dimension(:), intent(in), optional :: map_BA

    type(state_type), dimension(size(old_states)) :: old_fields_state, new_fields_state
    type(vector_field), pointer :: old_position, new_position
    integer :: i

    ewrite(1, *) "In interpolation_galerkin_multiple_states"

    call collapse_fields_in_state(old_states, old_fields_state)
    call collapse_fields_in_state(new_states, new_fields_state)
    call derive_collapsed_bcs(new_states, new_fields_state, bctype = "dirichlet")

    old_position => extract_vector_field(old_states(1), "Coordinate")
    new_position => extract_vector_field(new_states(1), "Coordinate")

    call interpolation_galerkin_scalars(old_fields_state, old_position, new_fields_state, new_position, map_BA=map_BA)

    do i = 1, size(old_fields_state)
      call deallocate(old_fields_state(i))
      call deallocate(new_fields_state(i))
    end do

    ewrite(1, *) "Exiting interpolation_galerkin_multiple_states"

  end subroutine interpolation_galerkin_multiple_states

  subroutine grandy_projection_multiple_states(old_states, new_states, map_BA)
    type(state_type), dimension(:), intent(inout) :: old_states, new_states
    type(ilist), dimension(:), intent(in), optional :: map_BA

    type(state_type), dimension(size(old_states)) :: old_fields_state, new_fields_state
    type(scalar_field), dimension(:), pointer :: old_fields, new_fields
    type(vector_field), pointer :: old_position, new_position
    integer :: i

    ewrite(1, *) "In grandy_projection_multiple_states"

    call collapse_fields_in_state(old_states, old_fields_state)
    call collapse_fields_in_state(new_states, new_fields_state)
    call collapse_state(old_fields_state, old_fields)
    call collapse_state(new_fields_state, new_fields)

    old_position => extract_vector_field(old_states(1), "Coordinate")
    new_position => extract_vector_field(new_states(1), "Coordinate")

    call grandy_projection_scalars(old_fields, old_position, new_fields, new_position, map_BA=map_BA)

    do i = 1, size(old_fields_state)
      call deallocate(old_fields_state(i))
      call deallocate(new_fields_state(i))
    end do

    deallocate(old_fields)
    deallocate(new_fields)

    ewrite(1, *) "Exiting grandy_projection_multiple_states"

  end subroutine grandy_projection_multiple_states

  subroutine grandy_projection_scalars(old_fields, old_position, new_fields, new_position, map_BA)
    !!< Grandy, 1999.
    !!< 10.1006/jcph.1998.6125
    type(scalar_field), dimension(:), intent(in) :: old_fields
    type(vector_field), intent(in) :: old_position
    type(vector_field), intent(in) :: new_position
    type(scalar_field), dimension(:), intent(inout) :: new_fields

    integer :: ele_A, ele_B, ele_C
    real :: vol_A, vol_B, vol_C
    integer :: dim

    type(ilist), dimension(:), intent(in), optional, target :: map_BA
    type(ilist), dimension(:), pointer :: lmap_BA
    type(inode), pointer :: llnode

    type(vector_field) :: intersection
    type(quadrature_type) :: supermesh_quad
    type(element_type) :: supermesh_shape

    real, dimension(size(old_fields)) :: integral_A, integral_B
    real, dimension(ele_ngi(old_fields(1), 1)) :: detwei_A

    type(scalar_field), dimension(size(old_fields)) :: pwc_B

    integer :: field, field_cnt
    type(state_type) :: projection_state
    character(len=OPTION_PATH_LEN) :: old_path
    integer :: stat
    real, dimension(new_position%dim, ele_loc(new_position, 1)) :: pos_B

    field_cnt = size(old_fields)
    dim = mesh_dim(new_position)


    ! Linear positions -- definitely linear positions.
    assert(old_position%mesh%shape%degree == 1)
    assert(continuity(old_position) >= 0)
    assert(continuity(new_position) >= 0)
    do field=1,field_cnt
      pwc_B(field) = piecewise_constant_field(new_position%mesh, trim(old_fields(field)%name) // "PWC")
      call zero(pwc_B(field))
    end do

    supermesh_quad = make_quadrature(vertices=ele_loc(new_position, 1), dim=dim, degree=1)
    supermesh_shape = make_element_shape(vertices=ele_loc(new_position, 1), dim=dim, degree=1, quad=supermesh_quad)

    call intersector_set_dimension(dim)
    if (present(map_BA)) then
      lmap_BA => map_BA
    else
      allocate(lmap_BA(ele_count(new_position)))
      lmap_BA = intersection_finder(new_position, old_position)
    end if

    do ele_B=1,ele_count(new_position)
      llnode => lmap_BA(ele_B)%firstnode
      integral_B = 0.0
      vol_B = simplex_volume(new_position, ele_B)
      pos_B = ele_val(new_position, ele_B)

      do while(associated(llnode))
        ele_A = llnode%value
        vol_A = simplex_volume(old_position, ele_A)
        integral_A = 0.0

        call transform_to_physical(old_position, ele_A, detwei_A)

        do field=1,field_cnt
          integral_A(field) = dot_product(ele_val_at_quad(old_fields(field), ele_A), detwei_A)
        end do

        intersection = intersect_elements(old_position, ele_A, pos_B, supermesh_shape)
        do ele_C=1,ele_count(intersection)
          vol_C = simplex_volume(intersection, ele_C)
          do field=1,field_cnt
            integral_B(field) = integral_B(field) + integral_A(field)*(vol_C/vol_A)
          end do
        end do

        do field=1,field_cnt
          call set(pwc_B(field), ele_B, integral_B(field) / vol_B)
        end do

        call deallocate(intersection)

        llnode => llnode%next
      end do
    end do

    call deallocate(supermesh_shape)
    call deallocate(supermesh_quad)

    ! Now call the Galerkin projection routines to translate from P0 to Pn.

    do field=1,field_cnt
      call insert(projection_state, pwc_B(field), trim(new_fields(field)%name))
      call insert(projection_state, new_position, "Coordinate")

      old_path = new_fields(field)%option_path
      new_fields(field)%option_path = "/temporary"
      call set_option("/temporary/diagnostic/source_field_name", trim(new_fields(field)%name), stat=stat)
      call set_solver_options("/temporary/diagnostic", ksptype='cg', pctype='eisenstat', rtol=1.0e-10, max_its=20000)
      call zero(new_fields(field))
      call calculate_galerkin_projection(projection_state, new_fields(field))
      call delete_option("/temporary")
      new_fields(field)%option_path = old_path

      call deallocate(projection_state)
      call deallocate(pwc_B(field))
    end do
    
    if (.not. present(map_BA)) then
      do ele_B=1,ele_count(new_position)
        call deallocate(lmap_BA(ele_B))
      end do
      deallocate(lmap_BA)
    end if

  end subroutine grandy_projection_scalars

end module conservative_interpolation_module
