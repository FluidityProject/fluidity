#include "fdebug.h"

module field_derivatives
    !!< This module contains code to compute the derivatives of
    !!< scalar fields. It uses superconvergent patch recovery, see
    !!< Zienkiewicz & Zhu, Int. J. Numer. Methods Eng, 33, 1331-1364 (1992)
    !!< At present it only computes up to the second derivative.
    !!< (Since the k-th derivative of a scalar field is a rank-k tensor,
    !!< anything with k > 2 rapidly becomes far too big to store in memory.)

    use elements
    use fetools, only: shape_shape, shape_dshape, dshape_outer_dshape
    use fields
    use halos
    use eventcounter
    use transform_elements
    use vector_tools
    use vector_set
    use node_boundary
    use surfacelabels
    use vtk_interfaces
    use superconvergence
    use state_module
    use boundary_conditions, only: get_entire_boundary_condition
    implicit none

    interface compute_hessian_real
      module procedure compute_hessian_var
    end interface

    interface differentiate_field_lumped
      module procedure differentiate_field_lumped_single, differentiate_field_lumped_multiple, &
          differentiate_field_lumped_vector
    end interface
    
    interface u_dot_nabla
      module procedure u_dot_nabla_scalar, &
        & u_dot_nabla_vector
    end interface u_dot_nabla

    interface grad
      module procedure grad_scalar, grad_vector, grad_vector_tensor
    end interface grad

    private

    public :: strain_rate, differentiate_field, grad, compute_hessian, &
      domain_is_2d, patch_type, get_patch_ele, get_patch_node, get_quadratic_fit_qf, curl, &
      get_quadratic_fit_eqf, div, u_dot_nabla, get_cubic_fit_cf, differentiate_field_lumped

    public :: compute_hessian_qf, compute_hessian_eqf, compute_hessian_var
    
    contains

    subroutine differentiate_field_spr(infield, positions, derivatives, outfields, accuracy_at_cost)
      !!< This subroutine takes in a scalar field infield,
      !!< an array of allocated scalar fields outfields, and an array of logicals
      !!< telling it what derivatives to take.
      !!< If derivatives(1), (2) and (3) are all .true., then
      !!< outfields will have all three directional derivatives of infield.
      !!< If derivatives(1) is .false. and (2) and (3) are true, then
      !!< outfields will contain the y and z derivatives, etc.

      type(scalar_field), intent(inout) :: infield
      type(vector_field), intent(in) :: positions
      logical, dimension(:), intent(in) :: derivatives
      type(scalar_field), dimension(:), intent(inout) :: outfields
      logical, intent(in), optional :: accuracy_at_cost !!< Don't check for duplicate
                                                        !!< superconvergent points.
                                                        !!< Setting this improves
                                                        !!< accuracy but greatly increases
                                                        !!< runtime.

      ! (n%loc x n%superconvergence%nsp x dim)
      real, dimension(:,:,:), allocatable ::  dnsp_t ! Holds the output of transform_superconvergent_to_physical
      ! superconvergent positions: positions%dim x shape%superconvergence%nsp
      real, dimension(:, :), allocatable :: superconvergent_positions, A, b ! Ax = b
      real, dimension(:), allocatable :: node_position

      real, dimension(MATRIX_SIZE_SPR) :: P

      integer :: i, node, sp, cnt, ele, j, stat
      real :: diffval ! the value of the derivative
      type(patch_type) :: patch
      type(element_type) :: shape, x_shape
      logical :: already_processed
      type(mesh_type) :: mesh
      
      integer :: level
      integer :: vset

      i = 0; node = 0; sp = 0; cnt = 0; ele = 0; j = 0
      diffval = 0.0

      mesh = infield%mesh
      call add_nelist(mesh)

      assert(size(derivatives) .le. 3)

      cnt = 0
      do i=1,size(derivatives)
        if (derivatives(i)) cnt = cnt + 1
      end do 

      assert(size(outfields) .ge. cnt)

      call initialise_boundcount(infield%mesh, positions) 

      ! FIXME: these lines assume (for efficiency reasons, since it's currently true)
      ! that each element of the mesh is the same.
      shape = infield%mesh%shape
      x_shape = positions%mesh%shape
      allocate(dnsp_t(shape%loc, shape%superconvergence%nsp, mesh_dim(infield%mesh)))
      allocate(superconvergent_positions(positions%dim, shape%superconvergence%nsp))
      allocate(A(MATRIX_SIZE_SPR, MATRIX_SIZE_SPR), b(MATRIX_SIZE_SPR, cnt))
      allocate(node_position(positions%dim))

      ! This visits all nodes, and sets each once.
      ! FIXME: generalise this to meshes that aren't linear tets.
      ! You'll need to visit each /vertex/ node once,
      ! and average the values computed for each midpoint node
      ! (as values for those are computed more than once)

      call vecset_create(vset)

      do node=1,infield%mesh%nodes                               ! loop over nodes requested
        A = 0.0; b = 0.0                                         ! clear the linear system
        if (node_lies_on_boundary(node)) then
          level = 2
        else
          level = 1
        end if
        patch = get_patch_ele(infield%mesh, node, level=level)    ! form patch

        node_position = node_val(positions, node)

        do ele=1,patch%count                                   ! loop over elements around node
          !shape = ele_shape(infield%mesh, patch%elements(ele)) ! get the element type

          ! get the derivatives of the basis functions at the superconvergent points
          assert(associated(shape%superconvergence))
          call transform_superconvergent_to_physical(ele_val(positions, patch%elements(ele)), x_shape, shape, dnsp_t)

          ! get the positions of the superconvergent points
          superconvergent_positions = ele_val_at_superconvergent(positions, patch%elements(ele))

          do sp=1,shape%superconvergence%nsp
            if (.not. present(accuracy_at_cost)) then
              call vecset_is_present(vset, superconvergent_positions(:, sp), already_processed)
              if (already_processed) cycle
            end if

            ! construct the matrix.
            A = A + compute_matrix_contribution_spr(superconvergent_positions(:, sp), shape)

            ! construct the rhs.
            cnt = 0
            do i=1,size(derivatives)

              if (derivatives(i)) then
                cnt = cnt+1
                diffval = dot_product(ele_val(infield, patch%elements(ele)), dnsp_t(:, sp, i)) ! compute the direct derivative at the point
                b(:, cnt) = b(:, cnt) + compute_rhs_contribution_spr(superconvergent_positions(:, sp), shape, diffval)
              end if
            end do ! loop over derivatives
          end do ! loop over superconvergent points
        end do ! loop over elements

        ! OK. Now we have the linear system to be solved.
        ! First solve it:

        call solve(A, b, stat)
        assert(stat == 0)

        ! The solutions are now in the memory of b.

        ! Now we need to multiply P * the solution to get the recovered derivative at the point.
        do i=1,positions%dim
          node_position(i) = positions%val(i,node)
        end do
        ! So now get P:
        P = getP_spr(node_position, shape)

        ! Now compute the derivative. (At last!)
        cnt = 0
        do i=1,size(derivatives)
          if (derivatives(i)) then
            cnt = cnt + 1
            diffval = dot_product(P, b(:, cnt))
            outfields(cnt)%val(node) = diffval
          end if
        end do

        deallocate(patch%elements) ! Have to do it here; otherwise memory leak
        call vecset_clear(vset)
      end do ! loop over nodes
      deallocate(dnsp_t, superconvergent_positions, A, b, node_position)
      call vecset_destroy(vset)

      ! Boundary values aren't trustworthy on 2d domains, so set them from internal values
      call differentiate_boundary_correction(outfields, positions, shape, count(derivatives))
      
      cnt = 0
      do i=1,size(derivatives)
         if (derivatives(i)) then
            cnt = cnt + 1
            call halo_update(outfields(cnt))
         end if
      end do

    end subroutine differentiate_field_spr

    subroutine grad_scalar(infield, positions, gradient)
      !!< This routine computes the gradient of a field.
      !!< For a continuous gradient this lumps the mass matrix
      !!< in the Galerkin projection.
      type(scalar_field), intent(in) :: infield
      type(vector_field), intent(in) :: positions
      type(vector_field), intent(inout) :: gradient
      type(scalar_field), dimension(gradient%dim) :: pardiff
      logical, dimension(gradient%dim) :: derivatives
      integer :: i, dim

      !! Field over the entire surface mesh containing bc values:
      type(scalar_field) :: bc_value
      !! Integer array of all surface elements indicating bc type::
      integer, dimension(:), allocatable :: bc_type

      dim = gradient%dim
      do i=1,dim
        pardiff(i) = extract_scalar_field(gradient, i)
      end do

      ! we need all derivatives
      derivatives = .true.

      if (infield%mesh%continuity<0) then
        !! required for dg gradient calculation
        allocate(bc_type(1:surface_element_count(infield)))
        call get_entire_boundary_condition(infield, (/"weakdirichlet"/), bc_value, bc_type)
        
        call differentiate_field(infield, positions, derivatives, pardiff, bc_value, bc_type)
        
        call deallocate(bc_value)
        deallocate(bc_type)
      else
        call differentiate_field(infield, positions, derivatives, pardiff)
      end if

    end subroutine grad_scalar

    subroutine grad_vector(infield, positions, gradient)
      !!< This routine computes the gradient of a field.
      type(vector_field), intent(in) :: infield
      type(vector_field), intent(in) :: positions
      type(vector_field), dimension(infield%dim), intent(inout) :: gradient
      type(scalar_field), dimension(gradient(1)%dim) :: pardiff
      type(scalar_field) :: component
      logical, dimension(gradient(1)%dim) :: derivatives
      integer :: i, j, dim

      !! Field over the entire surface mesh containing bc values:
      type(vector_field) :: bc_value
      type(scalar_field) :: bc_component_value
      !! Integer array of all surface elements indicating bc type::
      integer, dimension(:,:), allocatable :: bc_type
      integer, dimension(:), allocatable :: bc_component_type

      if (infield%mesh%continuity<0) then
        !! required for dg gradient calculation
        allocate(bc_type(infield%dim, 1:surface_element_count(infield)))
        allocate(bc_component_type(1:surface_element_count(infield)))
        call get_entire_boundary_condition(infield, (/"weakdirichlet"/), bc_value, bc_type)
      end if

      dim = gradient(1)%dim

      do j=1,infield%dim

        component = extract_scalar_field(infield, j)

        do i=1,dim
          pardiff(i) = extract_scalar_field(gradient(j), i)
        end do

        derivatives = .true.

        if (infield%mesh%continuity<0) then
          bc_component_value = extract_scalar_field(bc_value, j)
          bc_component_type = bc_type(j,:)
          call differentiate_field(component, positions, derivatives, pardiff, bc_component_value, bc_component_type)
        else
          call differentiate_field(component, positions, derivatives, pardiff)
        end if
      end do

      if (infield%mesh%continuity<0) then
        call deallocate(bc_value)
        deallocate(bc_type, bc_component_type)
      end if

    end subroutine grad_vector

    subroutine grad_vector_tensor(infield,positions,t_field)
      !!< This routine computes the full (tensor) grad of an infield vector field
      type(vector_field), intent(in) :: infield
      type(vector_field), intent(in) :: positions
      type(tensor_field), intent(inout) :: t_field

      type(scalar_field), dimension(infield%dim) :: pardiff
      type(scalar_field) :: component

      real, dimension(t_field%dim(1),t_field%dim(2)) :: t
      logical, dimension(infield%dim) :: derivatives
      integer :: i, j
      integer :: node

      !! Field over the entire surface mesh containing bc values:
      type(vector_field) :: bc_value
      type(scalar_field) :: bc_component_value
      !! Integer array of all surface elements indicating bc type::
      integer, dimension(:,:), allocatable :: bc_type
      integer, dimension(:), allocatable :: bc_component_type

      if (infield%mesh%continuity<0) then
        !! required for dg gradient calculation
        allocate(bc_type(infield%dim, 1:surface_element_count(infield)))
        allocate(bc_component_type(1:surface_element_count(infield)))
        call get_entire_boundary_condition(infield, (/"weakdirichlet"/), bc_value, bc_type)
      end if

      do j=1,infield%dim

        component = extract_scalar_field(infield, j)

        do i=1,infield%dim
          pardiff(i) = extract_scalar_field(t_field,i,j)
        end do

        derivatives = .true.

        if (infield%mesh%continuity<0) then
          bc_component_value = extract_scalar_field(bc_value, j)
          bc_component_type = bc_type(j,:)
          call differentiate_field(component, positions, derivatives, pardiff, bc_component_value, bc_component_type)
        else
          call differentiate_field(component, positions, derivatives, pardiff)
        end if

      end do

      if (infield%mesh%continuity<0) then
        call deallocate(bc_value)
        deallocate(bc_type, bc_component_type)
      end if
 
    end subroutine grad_vector_tensor

    subroutine strain_rate(infield,positions,t_field)
      !!< This routine computes the strain rate of an infield vector field
      type(vector_field), intent(in) :: infield
      type(vector_field), intent(in) :: positions
      type(tensor_field), intent(inout) :: t_field

      type(scalar_field), dimension(infield%dim) :: pardiff
      type(scalar_field) :: component

      real, dimension(t_field%dim(1),t_field%dim(2)) :: t
      logical, dimension(infield%dim) :: derivatives
      integer :: i, j
      integer :: node

      do j=1,infield%dim

        component = extract_scalar_field(infield, j)

        do i=1,infield%dim
          pardiff(i) = extract_scalar_field(t_field,i,j)
        end do

        derivatives = .true.

        call differentiate_field(component, positions, derivatives, pardiff)

      end do
      
      ! Computing the final strain rate tensor      
      do node=1,node_count(t_field)
           t=node_val(t_field, node)
           call set(t_field, node, (t+transpose(t))/2) 
      end do 
 
    end subroutine strain_rate

    subroutine differentiate_field_qf(infield, positions, derivatives, pardiff)
    !!< This routine computes the derivative using the QF (quadratic-fit)
    !!< approach described in:
    !!< M.G. Vallet et al., Numerical comparison of some Hessian recovery techniques 
    !!< Int. J. Numer. Meth. Engng., in press
      type(scalar_field), intent(in), target :: infield
      type(vector_field), intent(in) :: positions
      logical, dimension(:), intent(in) :: derivatives
      type(scalar_field), dimension(:), intent(inout) :: pardiff

      type(mesh_type), pointer :: mesh
      type(element_type) :: t_shape, x_shape

      integer :: node
      integer :: i, j

      real :: x, y, z
      real, dimension(3) :: diffvals
      real, dimension(MATRIX_SIZE_QF) :: b

      mesh => infield%mesh
      t_shape = ele_shape(mesh, 1)
      x_shape = ele_shape(positions, 1)
      do i=1,count(derivatives .eqv. .true.)
        call zero(pardiff(i))
      end do

      if (maxval(infield%val) == minval(infield%val)) then
        ewrite(2,*) "+++: Field constant; returning 0.0"
        return
      end if

      call add_nelist(mesh)
      call initialise_boundcount(mesh, positions)

      do node=1,node_count(infield)
        x = node_val(positions, 1, node); y = node_val(positions, 2, node); z = node_val(positions, 3, node)
        b = get_quadratic_fit_qf(infield, positions, node)

        diffvals(1) = b(2) + 2 * b(5) * x + b(8) * y + b(9) * z + b(11) * y * z
        diffvals(2) = b(3) + 2 * b(6) * y + b(8) * x + b(10) * z + b(11) * x * z
        diffvals(3) = b(4) + 2 * b(7) * z + b(9) * x + b(10) * y + b(11) * x * y

        j = 1
        do i=1,3
          if (derivatives(i)) then
            pardiff(j)%val(node) = diffvals(i)
            j = j + 1
          end if
        end do
      end do

      call differentiate_boundary_correction(pardiff, positions, t_shape, count(derivatives .eqv. .true.))
      call differentiate_squash_pseudo2d(pardiff)

      do i=1,count(derivatives .eqv. .true.)
         call halo_update(pardiff(i))
      end do
    end subroutine differentiate_field_qf

    subroutine compute_hessian_qf(infield, positions, hessian)
    !!< This routine computes the hessian using the QF (quadratic-fit)
    !!< approach described in:
    !!< M.G. Vallet et al., Numerical comparison of some Hessian recovery techniques 
    !!< Int. J. Numer. Meth. Engng., in press
      type(scalar_field), intent(in), target :: infield
      type(vector_field), intent(in) :: positions
      type(tensor_field), intent(inout) :: hessian

      type(mesh_type), pointer :: mesh
      type(element_type) :: t_shape, x_shape

      integer :: node

      real :: b(MATRIX_SIZE_QF)
      real :: x, y, z

#ifdef QF_DEBUG
      type(scalar_field) :: variance
      type(patch_type) :: patch
      integer :: i, nnode
#endif

      mesh => infield%mesh

#ifdef QF_DEBUG
      call allocate(variance, mesh, "Root-Variance")
      call zero(variance)
#endif

      t_shape = ele_shape(mesh, 1)
      x_shape = ele_shape(positions, 1)
      call zero(hessian)
      if (maxval(infield%val) == minval(infield%val)) then
        ewrite(2,*) "+++: Field constant; returning 0.0"
        return
      end if

      call add_nelist(mesh)
      call initialise_boundcount(mesh, positions)

      do node=1,node_count(infield)
        x = node_val(positions, 1, node); y = node_val(positions, 2, node); z = node_val(positions, 3, node);
        b = get_quadratic_fit_qf(infield, positions, node)

#ifdef QF_DEBUG
        patch = get_patch_node(mesh, node, level=2, min_nodes=MATRIX_SIZE_QF)
        do i=1,patch%count
          nnode = patch%elements(i)
          call addto(variance, node, (node_val(infield, nnode) - evaluate_qf(b, node_val(positions, nnode)))**2)
        end do
        variance%val(node) = sqrt(variance%val(node))
#endif

        hessian%val(1, 1, node) = 2 * b(5)
        hessian%val(1, 2, node) = b(8) + z * b(11)
        hessian%val(1, 3, node) = b(9) + y * b(11)
        hessian%val(2, 1, node) = b(8) + z * b(11)
        hessian%val(2, 2, node) = 2 * b(6)
        hessian%val(2, 3, node) = b(10) + x * b(11)
        hessian%val(3, 1, node) = b(9) + y * b(11)
        hessian%val(3, 2, node) = b(10) + x * b(11)
        hessian%val(3, 3, node) = 2 * b(7)
      end do

      call hessian_boundary_correction(hessian, positions, t_shape)
      call hessian_squash_pseudo2d(hessian)

      call halo_update(hessian)

#ifdef QF_DEBUG
      call vtk_write_fields("qf_debug", 0, positions, mesh, sfields=(/infield, variance/))
      call deallocate(variance)
#endif
    end subroutine compute_hessian_qf

    recursive function get_quadratic_fit_qf(infield, positions, node, level, maxlevel) result(b)
      !!< Fit a quadratic function to infield around the node,
      !!< with a least squares approach.
      type(scalar_field), intent(in) :: infield
      type(vector_field), intent(in) :: positions
      integer :: node
      integer, optional :: level, maxlevel
      real, dimension(MATRIX_SIZE_QF) :: b

      type(patch_type) :: node_patch
      real, dimension(MATRIX_SIZE_QF, MATRIX_SIZE_QF) :: A
      real, dimension(MATRIX_SIZE_QF, 1) :: b_tmp
      real, dimension(MATRIX_SIZE_QF_2D, MATRIX_SIZE_QF_2D) :: A_2D
      real, dimension(MATRIX_SIZE_QF_2D, 1) :: b_tmp_2D
      integer :: i, nnode, stat, llevel, lmaxlevel
      type(mesh_type) :: mesh

      if (present(level)) then
        llevel = level
      else
        llevel = 2
      end if

      if (present(maxlevel)) then
        lmaxlevel = maxlevel
      else
        lmaxlevel = llevel + 1
      end if

      mesh = infield%mesh
      node_patch = get_patch_node(mesh, node, level=llevel, min_nodes=MATRIX_SIZE_QF)
      A = 0.0; b_tmp = 0.0
      do i=1,node_patch%count
        nnode = node_patch%elements(i)

        A = A + compute_matrix_contribution_qf(node_val(positions, nnode))
        b_tmp(:, 1) = b_tmp(:, 1) + compute_rhs_contribution_qf(node_val(positions, nnode), node_val(infield, nnode))
      end do

      if (pseudo2d_coord /= 0) then
        if (pseudo2d_coord == 1) then
          A_2D = A(QF_2D_X, QF_2D_X)
          b_tmp_2D(:, 1) = b_tmp(QF_2D_X, 1)
          call solve(A_2D, b_tmp_2D, stat)
          b = 0.0; b(QF_2D_X) = b_tmp_2D(:, 1)
        else if (pseudo2d_coord == 2) then
          A_2D = A(QF_2D_Y, QF_2D_Y)
          b_tmp_2D(:, 1) = b_tmp(QF_2D_Y, 1)
          call solve(A_2D, b_tmp_2D, stat)
          b = 0.0; b(QF_2D_Y) = b_tmp_2D(:, 1)
        else if (pseudo2d_coord == 3) then
          A_2D = A(QF_2D_Z, QF_2D_Z)
          b_tmp_2D(:, 1) = b_tmp(QF_2D_Z, 1)
          call solve(A_2D, b_tmp_2D, stat)
          b = 0.0; b(QF_2D_Z) = b_tmp_2D(:, 1)
        end if
      else
        call solve(A, b_tmp, stat)
        b = b_tmp(:, 1)
      end if


      if (llevel < lmaxlevel) then
        if (stat /= 0) then
          ! If the solver fails, go one more level deep
          ! to get enough equations to do the fitting
          b = get_quadratic_fit_qf(infield, positions, node, level=llevel+1)
        end if
      end if

      deallocate(node_patch%elements) ! Have to do it here; otherwise memory leak
    end function get_quadratic_fit_qf

    subroutine compute_hessian_int(infield, positions, hessian)
    !!< This routine computes the hessian using integration by parts.
    !!< See Buscaglia and Dari, Int. J. Numer. Meth. Engng., 40, 4119-4136 (1997)
      type(scalar_field), intent(inout) :: infield
      type(vector_field), intent(in) :: positions
      type(tensor_field), intent(inout) :: hessian

      ! For now, assume only one element type in the mesh.

      real, dimension(ele_ngi(infield, 1)) :: detwei
      real, dimension(ele_loc(infield, 1), ele_ngi(infield, 1), mesh_dim(infield)) :: dt_t
      type(element_type), pointer :: t_shape
      real, dimension(mesh_dim(infield), mesh_dim(infield), ele_loc(infield, 1), ele_loc(infield, 1)) :: r
      real, dimension(mesh_dim(infield), mesh_dim(infield), ele_loc(infield, 1)) :: r_ele
      type(scalar_field), target  :: lumped_mass_matrix
      real, dimension(ele_loc(infield, 1), ele_loc(infield, 1)) :: mass_matrix
      type(mesh_type) :: mesh

      integer :: ele, node,i,j


      call zero(hessian)
      if (maxval(infield%val) == minval(infield%val)) then
        ewrite(2,*) "+++: Field constant; returning 0.0"
        return
      end if

      mesh = infield%mesh
      call add_nelist(mesh)
      call initialise_boundcount(infield%mesh, positions)

      call allocate(lumped_mass_matrix, infield%mesh, "Lumped mass matrix")
      call zero(lumped_mass_matrix)

      t_shape => ele_shape(infield, 1)

      do ele=1,element_count(infield)
        ! Compute detwei.
        call transform_to_physical(positions, ele, t_shape, dshape=dt_t, detwei=detwei)

        ! Compute the tensor representing grad(N) grad(N)
        r = dshape_outer_dshape(dt_t, dt_t, detwei)
        !r_ele = 0.5 * (tensormul(r, ele_val(infield, ele), 3) + tensormul(r, ele_val(infield, ele), 4))
        !r_ele = tensormul(r, ele_val(infield, ele), 4)
        
        r_ele = 0.
        do i = 1,size(r,1)
           do j = 1,size(r,2)
              r_ele(i,j,:) = r_ele(i,j,:) + &
                   matmul(r(i,j,:,:),ele_val(infield,ele))
           end do
        end do
        call addto(hessian, ele_nodes(infield, ele), r_ele)

        ! Lump the mass matrix
        mass_matrix = shape_shape(t_shape, t_shape, detwei)
        call addto(lumped_mass_matrix, ele_nodes(infield, ele), sum(mass_matrix, 2))
      end do

      do node=1,node_count(infield)
        hessian%val(:, :, node) = (-1.0 / node_val(lumped_mass_matrix, node)) * hessian%val(:, :, node)
        !hessian%val(:, :, node) = (-1) * hessian%val(:, :, node)
      end do

      call hessian_boundary_correction(hessian, positions, t_shape)
      call deallocate(lumped_mass_matrix)
    end subroutine compute_hessian_int

    subroutine differentiate_boundary_correction(pardiff, positions, t_shape, count)
      !!< Implement the boundary correction routine for first derivatives.
      type(scalar_field), dimension(:), intent(inout), target :: pardiff
      type(vector_field), intent(in) :: positions
      type(element_type), intent(in) :: t_shape
      integer, intent(in) :: count

      type(scalar_field) :: node_weights
      type(mesh_type), pointer :: mesh

      integer :: i, j, k, dim, node, nnode, ele
      integer, dimension(:), pointer :: neighbour_elements, neighbour_nodes

      real :: sum_weights
      real, dimension(ele_ngi(pardiff(1), 1)) :: detwei
      real, dimension(ele_loc(pardiff(1), 1), ele_loc(pardiff(1), 1)) :: mass_matrix
      real, dimension(mesh_dim(pardiff(1))) :: old_val
      
      logical :: has_neighbouring_interior_node
      type(patch_type) :: node_patch
      type(csr_sparsity), pointer :: nelist

      !type(scalar_field) :: boundcounts, node_numbers, patch_nodes

      mesh => pardiff(1)%mesh
      nelist => extract_nelist(mesh)
      call allocate(node_weights, mesh, "NodeWeights")

      dim = mesh_dim(mesh)

      call initialise_boundcount(mesh, positions)

      do i=1,dim
        do node=1,node_count(mesh)
          if (node_boundary_count(node) >= get_expected_boundcount() + i) then
            call zero(node_weights)
            has_neighbouring_interior_node = .false.
            ! First we need to compute the weights for each neighbouring node.
            neighbour_elements => row_m_ptr(nelist, node)
            do j=1,size(neighbour_elements)
              ele = neighbour_elements(j)
              neighbour_nodes => ele_nodes(mesh, ele)
              call transform_to_physical(positions, ele, detwei=detwei)
              mass_matrix = shape_shape(t_shape, t_shape, detwei)
              ! In words: find the row of the mass matrix corresponding to the node we're interested in,
              ! and stuff the integral of the shape functions into node_weights.
              call addto(node_weights, neighbour_nodes, mass_matrix(:, find(neighbour_nodes, node)))
    
              ! Also: find out if the node has /any/ neighbouring interior nodes.
              if (.not. has_neighbouring_interior_node) then
                do k=1,size(neighbour_nodes)
                  nnode = neighbour_nodes(k)
                  if (.not. node_lies_on_boundary(nnode)) then
                    has_neighbouring_interior_node = .true.
                    exit
                  end if
                end do
              end if
            end do
    
            ! Now that we have the weights, let us use them.
    
            node_patch = get_patch_node(mesh, node)
            sum_weights = 0.0
            forall (j=1:count)
              old_val(j) = pardiff(j)%val(node)
              pardiff(j)%val(node) = 0.0
            end forall
            do j=1,node_patch%count
              nnode = node_patch%elements(j)
              ! If it's on the boundary, no ...
              if (&
              (has_neighbouring_interior_node .and. &
              (.not. node_lies_on_boundary(nnode))) &
              .or. ((.not. has_neighbouring_interior_node) .and. node_boundary_count(nnode) < node_boundary_count(node))) then
                sum_weights = sum_weights + node_val(node_weights, nnode)
                do k=1,count
                  pardiff(k)%val(node) = pardiff(k)%val(node) +  &
                                         pardiff(k)%val(nnode) * node_val(node_weights, nnode)
                end do
              end if
            end do
    
            if (sum_weights == 0.0) then
              forall (k=1:count)
                pardiff(k)%val(node) = old_val(k)
              end forall
            else
              forall (k=1:count)
                pardiff(k)%val(node) = pardiff(k)%val(node) / sum_weights
              end forall
            end if
            deallocate(node_patch%elements)
          end if
        end do
      end do
      call deallocate(node_weights)
    end subroutine differentiate_boundary_correction

    subroutine hessian_boundary_correction(hessian, positions, t_shape)
      !!< Implement the hessian boundary correction routine.
      type(tensor_field), intent(inout), target :: hessian
      type(vector_field), intent(in) :: positions
      type(element_type), intent(in) :: t_shape

      type(scalar_field) :: node_weights
      type(mesh_type), pointer :: mesh

      integer :: i, j, k, dim, node, nnode, ele
      integer, dimension(:), pointer :: neighbour_elements, neighbour_nodes

      real :: sum_weights
      real, dimension(ele_ngi(hessian, 1)) :: detwei
      real, dimension(ele_loc(hessian, 1), ele_loc(hessian, 1)) :: mass_matrix
      real, dimension(hessian%dim(1), hessian%dim(2)) :: old_val
      
      logical :: has_neighbouring_interior_node
      type(patch_type) :: node_patch
      type(csr_sparsity), pointer :: nelist

      assert(hessian%dim(1)==hessian%dim(2))

      mesh => hessian%mesh
      nelist => extract_nelist(mesh)
      call allocate(node_weights, mesh, "Node weights")
      dim = hessian%dim(1)

      call initialise_boundcount(hessian%mesh, positions)

      do i=1,dim
        do node=1,node_count(hessian)
          if (node_boundary_count(node) >= get_expected_boundcount() + i) then
            call zero(node_weights)
            has_neighbouring_interior_node = .false.
            ! First we need to compute the weights for each neighbouring node.
            neighbour_elements => row_m_ptr(nelist, node)
            do j=1,size(neighbour_elements)
              ele = neighbour_elements(j)
              neighbour_nodes => ele_nodes(hessian, ele)
              call transform_to_physical(positions, ele, detwei=detwei)
              mass_matrix = shape_shape(t_shape, t_shape, detwei)
              ! In words: find the row of the mass matrix corresponding to the node we're interested in,
              ! and stuff the integral of the shape functions into node_weights.
              call addto(node_weights, neighbour_nodes, mass_matrix(:, find(neighbour_nodes, node)))
    
              ! Also: find out if the node has /any/ neighbouring interior nodes.
              if (.not. has_neighbouring_interior_node) then
                do k=1,size(neighbour_nodes)
                  nnode = neighbour_nodes(k)
                  if (.not. node_lies_on_boundary(nnode)) then
                    has_neighbouring_interior_node = .true.
                    exit
                  end if
                end do
              end if
            end do
    
            ! Now that we have the weights, let us use them.
    
            node_patch = get_patch_node(mesh, node)
            sum_weights = 0.0
            old_val = hessian%val(:, :, node)
            hessian%val(:, :, node) = 0.0
            do j=1,node_patch%count
              nnode = node_patch%elements(j)
              ! If it's on the boundary, no ...
              if (&
              (has_neighbouring_interior_node .and. &
              (.not. node_lies_on_boundary(nnode))) &
              .or. ((.not. has_neighbouring_interior_node) .and. node_boundary_count(nnode) < node_boundary_count(node))) then
                sum_weights = sum_weights + node_val(node_weights, nnode)
                hessian%val(:, :, node) = hessian%val(:, :, node) +  &
                                          hessian%val(:, :, nnode) * node_val(node_weights, nnode)
              end if
            end do
    
            if (sum_weights == 0.0) then
              hessian%val(:, :, node) = old_val
            else
              hessian%val(:, :, node) = hessian%val(:, :, node) / sum_weights
            end if
            deallocate(node_patch%elements)
          end if
        end do
      end do

      call deallocate(node_weights)
    end subroutine hessian_boundary_correction

    subroutine hessian_squash_pseudo2d(hessian)
      !!< Squash derivatives in directions where no dynamics occur.
      type(tensor_field), intent(inout) :: hessian

      integer :: node

      if (pseudo2d_coord == 0) then
        return
      end if

      do node=1,node_count(hessian)
        hessian%val(:, pseudo2d_coord, node) = 0.0
        hessian%val(pseudo2d_coord, :, node) = 0.0
      end do
    end subroutine hessian_squash_pseudo2d

    subroutine differentiate_squash_pseudo2d(pardiff)
      !!< Squash derivatives in directions where no dynamics occur.
      type(scalar_field), dimension(:), intent(inout) :: pardiff

      if (pseudo2d_coord == 0) then
        return
      end if

      call zero(pardiff(pseudo2d_coord))
    end subroutine differentiate_squash_pseudo2d

    function find(array, val) result(loc)
      !!< Find the first instance of val in array.
      integer, intent(in), dimension(:) :: array
      integer, intent(in) :: val
      integer :: i, loc

      loc = -1
      do i=1,size(array)
        if (array(i) == val) then
          loc = i
          return
        end if
      end do
    end function find

    subroutine compute_hessian_spr(infield, positions, hessian, accuracy_at_cost)
    !!< This routine computes the hessian, applying differentiate_field_spr multiple times.
      type(scalar_field), intent(inout) :: infield
      type(vector_field), intent(in) :: positions
      type(tensor_field), intent(inout) :: hessian
      logical, intent(in), optional :: accuracy_at_cost

      type(scalar_field) :: pardiff(positions%dim), temp_fields(positions%dim) ! temp_fields not allocated
      logical :: derivatives(positions%dim)
      integer :: i, dim
!      integer, dimension(:), pointer :: elements
!      real, dimension(positions%dim, positions%dim) :: org_evectors, evectors
!      real, dimension(positions%dim) :: org_evalues, evalues
!      integer :: j, k, neigh_count, node

      if (maxval(infield%val) == minval(infield%val)) then
        hessian%val = 0.0
        ewrite(2,*) "+++: Field constant; returning 0.0"
        return
      end if

      call initialise_boundcount(infield%mesh, positions)

      dim = positions%dim
      do i=1,dim
        call allocate(pardiff(i), infield%mesh) ! allocate the partial derivatives
      end do

      ! First get the first derivatives.

      if (dim == 2) then
        derivatives(1) = .true. ; derivatives(2) = .true.
        call differentiate_field_spr(infield, positions, derivatives, pardiff, accuracy_at_cost)

        temp_fields(1) = extract_scalar_field_from_tensor_field(hessian, 1, 1) 
        temp_fields(2) = extract_scalar_field_from_tensor_field(hessian, 1, 2)
        derivatives(1) = .true. ; derivatives(2) = .false.
        call differentiate_field_spr(pardiff(1), positions, derivatives, temp_fields, accuracy_at_cost)
        
        ! Now the state of hessian is
        ! [ d^2f/dx^2 d^2f/dydx ]
        ! [    ??         ??    ]

        temp_fields(1) = extract_scalar_field_from_tensor_field(hessian, 2, 1)
        temp_fields(1)%val(:) = temp_fields(2)%val

        ! Now the state of hessian is
        ! [ d^2f/dx^2 d^2f/dydx ]
        ! [ d^2f/dxdy     ??    ]

        temp_fields(1) = extract_scalar_field_from_tensor_field(hessian, 2, 2)
        derivatives(1) = .false. ; derivatives(2) = .true.
        call differentiate_field_spr(pardiff(2), positions, derivatives, temp_fields, accuracy_at_cost)

        ! Now the state of hessian is
        ! [ d^2f/dx^2 d^2f/dydx ]
        ! [ d^2f/dydy d^2f/dy^2 ]

      end if

      if (dim == 3) then

        ! if the domain is pseudo-2d, we don't want to take z derivatives
        if (.not. domain_is_2d()) then
          derivatives(1) = .true. ; derivatives(2) = .true. ; derivatives(3) = .true.
          call differentiate_field_spr(infield, positions, derivatives, pardiff, accuracy_at_cost)
        else if (domain_is_2d_x()) then
          derivatives(1) = .false. ; derivatives(2) = .true. ; derivatives(3) = .true.
          call differentiate_field_spr(infield, positions, derivatives, pardiff, accuracy_at_cost)
          pardiff(1)%val = 0.0
        else if (domain_is_2d_y()) then
          derivatives(1) = .true. ; derivatives(2) = .false. ; derivatives(3) = .true.
          call differentiate_field_spr(infield, positions, derivatives, pardiff, accuracy_at_cost)
          pardiff(2)%val = 0.0
        else if (domain_is_2d_z()) then
          derivatives(1) = .true. ; derivatives(2) = .true. ; derivatives(3) = .false.
          call differentiate_field_spr(infield, positions, derivatives, pardiff, accuracy_at_cost)
          pardiff(3)%val = 0.0
        end if
        ewrite(2,*) "+++: Gradient computed"

        temp_fields(1) = extract_scalar_field_from_tensor_field(hessian, 1, 1)
        temp_fields(2) = extract_scalar_field_from_tensor_field(hessian, 1, 2)
        temp_fields(3) = extract_scalar_field_from_tensor_field(hessian, 1, 3)
        if (.not. domain_is_2d()) then
          derivatives(1) = .true. ; derivatives(2) = .true.; derivatives(3) = .true.
          call differentiate_field_spr(pardiff(1), positions, derivatives, temp_fields, accuracy_at_cost)
        else if (domain_is_2d_x()) then
          temp_fields(1)%val = 0.0
          temp_fields(2)%val = 0.0
          temp_fields(3)%val = 0.0
        else if (domain_is_2d_y()) then
          derivatives(1) = .true. ; derivatives(2) = .false.; derivatives(3) = .true.
          call differentiate_field_spr(pardiff(1), positions, derivatives, temp_fields, accuracy_at_cost)
          temp_fields(2)%val = 0.0
        else if (domain_is_2d_z()) then
          derivatives(1) = .true. ; derivatives(2) = .true.; derivatives(3) = .false.
          call differentiate_field_spr(pardiff(1), positions, derivatives, temp_fields, accuracy_at_cost)
          temp_fields(3)%val = 0.0
        end if

        ! Now the state of hessian is
        ! [ d^2f/dx^2 d^2f/dydx d^2f/dzdx ]
        ! [     ??        ??        ??    ]
        ! [     ??        ??        ??    ]

        temp_fields(1) = extract_scalar_field_from_tensor_field(hessian, 2, 1)
        temp_fields(1)%val(:) = temp_fields(2)%val
        temp_fields(1) = extract_scalar_field_from_tensor_field(hessian, 3, 1)
        temp_fields(1)%val(:) = temp_fields(3)%val

        ewrite(2,*) "+++: 2nd-order x derivatives computed"

        ! Now the state of hessian is
        ! [ d^2f/dx^2 d^2f/dydx d^2f/dzdx ]
        ! [ d^2f/dydx     ??        ??    ]
        ! [ d^2f/dzdx     ??        ??    ]

        temp_fields(1) = extract_scalar_field_from_tensor_field(hessian, 2, 2)
        temp_fields(2) = extract_scalar_field_from_tensor_field(hessian, 2, 3)
        if (.not. domain_is_2d()) then
          derivatives(1) = .false.; derivatives(2) = .true. ; derivatives(3) = .true. ! only want the y,z derivatives
          call differentiate_field_spr(pardiff(2), positions, derivatives, temp_fields, accuracy_at_cost)
        else if (domain_is_2d_x()) then
          derivatives(1) = .false.; derivatives(2) = .true. ; derivatives(3) = .true.
          call differentiate_field_spr(pardiff(2), positions, derivatives, temp_fields, accuracy_at_cost)
        else if (domain_is_2d_y()) then
          temp_fields(1)%val = 0.0
          temp_fields(2)%val = 0.0
        else if (domain_is_2d_z()) then
          derivatives(1) = .false.; derivatives(2) = .true. ; derivatives(3) = .false.
          call differentiate_field_spr(pardiff(2), positions, derivatives, temp_fields, accuracy_at_cost)
          temp_fields(2)%val = 0.0
        end if

        ! Now the state of hessian is
        ! [ d^2f/dx^2 d^2f/dydx d^2f/dzdx ]
        ! [ d^2f/dydx d^2f/dy^2 d^2f/dydz ]
        ! [ d^2f/dzdx     ??        ??    ]

        ewrite(2,*) "+++: 2nd-order y derivatives computed"

        temp_fields(1) = extract_scalar_field_from_tensor_field(hessian, 3, 2)
        temp_fields(1)%val(:) = temp_fields(2)%val

        ! Now the state of hessian is
        ! [ d^2f/dx^2 d^2f/dydx d^2f/dzdx ]
        ! [ d^2f/dydx d^2f/dy^2 d^2f/dydz ]
        ! [ d^2f/dzdx d^2f/dydx     ??    ]

        temp_fields(1) = extract_scalar_field_from_tensor_field(hessian, 3, 3)
        if (.not. domain_is_2d_z()) then
          derivatives(1) = .false.; derivatives(2) = .false. ; derivatives(3) = .true. ! only want z derivative
          call differentiate_field_spr(pardiff(3), positions, derivatives, temp_fields, accuracy_at_cost)
        else
          temp_fields(1)%val = 0.0
        end if

        ! Now the state of hessian is
        ! [ d^2f/dx^2 d^2f/dydx d^2f/dzdx ]
        ! [ d^2f/dydx d^2f/dy^2 d^2f/dydz ]
        ! [ d^2f/dzdx d^2f/dydx d^2f/dz^2 ]

        ewrite(2,*) "+++: 2nd-order z derivatives computed"
      end if
      
      do i=1,dim
        deallocate(pardiff(i)%val)
      end do
    end subroutine compute_hessian_spr

    subroutine compute_hessian_var(infield, positions, hessian)
    !!< This routine computes the hessian using a weak finite element formulation.
      type(scalar_field), intent(in) :: infield
      type(vector_field), intent(in) :: positions
      type(tensor_field), intent(inout), target :: hessian

      type(vector_field), target :: gradient
      type(mesh_type), pointer :: mesh

      real, dimension(ele_ngi(positions, 1)) :: detwei
      real, dimension(ele_loc(infield, 1), ele_ngi(infield, 1), mesh_dim(infield)) :: dt_t
      real, dimension(ele_loc(hessian, 1), ele_ngi(hessian, 1), mesh_dim(hessian)) :: dh_t
      type(element_type), pointer :: t_shape, h_shape
      real, dimension(mesh_dim(infield), ele_loc(hessian, 1), ele_loc(infield, 1)) :: r
      real, dimension(mesh_dim(infield), ele_loc(hessian, 1)) :: r_grad_ele
      real, dimension(mesh_dim(hessian), ele_loc(hessian, 1), ele_loc(hessian, 1)) :: r_hess
      real, dimension(mesh_dim(hessian), mesh_dim(hessian), ele_loc(hessian, 1)) :: r_hess_ele
      type(scalar_field)  :: lumped_mass_matrix
      real, dimension(ele_loc(hessian, 1), ele_loc(hessian, 1)) :: mass_matrix
      integer :: dim, i, j
      integer :: node, ele
      real, dimension(:, :), pointer :: hess_ptr

      mesh => hessian%mesh
      dim = mesh_dim(mesh)

      call zero(hessian)
      if (maxval(infield%val) == minval(infield%val)) then
        ewrite(2,*) "+++: Field constant; returning 0.0"
        return
      end if

      call allocate(lumped_mass_matrix, mesh, "Lumped mass matrix")
      call allocate(gradient, dim, mesh, "Gradient")

      call add_nelist(mesh)
      call initialise_boundcount(mesh, positions)

      call zero(lumped_mass_matrix)
      call zero(gradient)

      t_shape => ele_shape(infield, 1)
      h_shape => ele_shape(hessian, 1)

      ! First, compute gradient and mass matrix.
      do ele=1,element_count(infield)
        ! Compute detwei.
        call transform_to_physical(positions, ele, t_shape, dshape=dt_t, detwei=detwei)

        r = shape_dshape(h_shape, dt_t, detwei)
        r_grad_ele = tensormul(r, ele_val(infield, ele), 3)

        call addto(gradient, ele_nodes(gradient, ele), r_grad_ele)

        ! Lump the mass matrix
        mass_matrix = shape_shape(h_shape, h_shape, detwei)
        call addto(lumped_mass_matrix, ele_nodes(lumped_mass_matrix, ele), sum(mass_matrix, 2))
      end do

      do node=1,node_count(gradient)
        do i=1,dim
          gradient%val(i,node) = gradient%val(i,node) / node_val(lumped_mass_matrix, node)
        end do
      end do

      ! Testing: does this cause the lock exchange result to fail?
      !do i=1,dim
      !  grad_components(i) = extract_scalar_field(gradient, i)
      !end do
      !call differentiate_boundary_correction(grad_components, positions, x_shape, t_shape, dim)

      do ele=1,element_count(infield)
        call transform_to_physical(positions, ele, h_shape, dshape=dh_t, detwei=detwei)
        r_hess = shape_dshape(h_shape, dh_t, detwei)
        do i=1,dim
          r_hess_ele(i, :, :) = tensormul(r_hess, ele_val(gradient, i, ele), 3)
        end do
        call addto(hessian, ele_nodes(hessian, ele), r_hess_ele)
      end do

      do node=1,node_count(hessian)
        hess_ptr => hessian%val(:, :, node)
        hess_ptr = hess_ptr / node_val(lumped_mass_matrix, node)
         
        ! Now we need to make it symmetric, see?
        do i=1,dim
          do j=i+1,dim
            hess_ptr(i, j) = (hess_ptr(i, j) + hess_ptr(j, i)) / 2.0
            hess_ptr(j, i) = hess_ptr(i, j)
          end do
        end do
      end do

      call hessian_boundary_correction(hessian, positions, h_shape)

      call deallocate(lumped_mass_matrix)
      call deallocate(gradient)
    end subroutine compute_hessian_var

    subroutine differentiate_field_lumped_multiple(infields, positions, derivatives, pardiff)
    !!< This routine computes the first derivatives using a weak finite element formulation.
      type(scalar_field), dimension(:), intent(in) :: infields
      type(vector_field), intent(in) :: positions
      logical, dimension(:), intent(in) :: derivatives
      type(scalar_field), dimension(:,:), target, intent(inout) :: pardiff

      type(vector_field), dimension(size(infields)), target :: gradient
      type(mesh_type), pointer :: mesh

      type(scalar_field)  :: lumped_mass_matrix, inverse_lumped_mass
      logical, dimension( mesh_dim(infields(1)) ):: compute
      integer :: i, j, k
      integer :: ele

      mesh => pardiff(1, 1)%mesh

      do i=1, size(infields)
        ! don't compute if the field is constant
        compute(i)= (maxval(infields(i)%val) /= minval(infields(i)%val))
        ! check the infield is continuous!!!!
        if (infields(i)%mesh%continuity<0) then
          ewrite(0,*) "If the following error is directly due to user input"
          ewrite(0,*) "a check and a more helpful error message should be inserted in"
          ewrite(0,*) "the calling routine (outside field_derivatives) - please mantis this:"
          ewrite(0,*) "Error has occured in differentiate_field_lumped_multiple, with field, ", trim(infields(i)%name)
          FLAbort("The field_derivatives code cannot take the derivative of a discontinuous field")
        end if
      end do

      call allocate(lumped_mass_matrix, mesh, "LumpedMassMatrix")
      call zero(lumped_mass_matrix)
      
      do i=1, size(infields)
        if (compute(i)) then
          call allocate(gradient(i), positions%dim, mesh, "Gradient")
          call zero(gradient(i))
        end if
      end do


      ! First, compute gradient and mass matrix.
      do ele=1, element_count(mesh)
        call differentiate_field_ele(ele)
      end do
        
      do i=1, size(infields)
        if (compute(i)) then
          k=0
          do j=1, positions%dim
            if (derivatives(j)) then
              k=k+1
              call set( pardiff(k,i), gradient(i), dim=j )
            end if
          end do
        else
          do k=1, size(pardiff,1)
            call zero(pardiff(k,i))
          end do
        end if
      end do

      ! invert the lumped mass matrix
      call allocate(inverse_lumped_mass, mesh, "InverseLumpedMassMatrix")
      call invert(lumped_mass_matrix, inverse_lumped_mass)
      
      ! compute pardiff=M^-1*pardiff
      do i=1, size(infields)
        do k=1, size(pardiff,1)
          call scale(pardiff(k,i), inverse_lumped_mass)
        end do
      end do

      do i=1, size(infields)
        if (compute(i)) then
          call deallocate(gradient(i))
        end if
      end do
      call deallocate(lumped_mass_matrix)
      call deallocate(inverse_lumped_mass)

      contains
      
      subroutine differentiate_field_ele(ele)
        integer, intent(in):: ele
      
        real, dimension(mesh_dim(mesh), ele_loc(mesh, ele), ele_loc(infields(1), ele)) :: r
        real, dimension(ele_ngi(mesh, ele)) :: detwei
        real, dimension(ele_loc(infields(1), ele), ele_ngi(infields(1), ele), mesh_dim(infields(1))) :: dt_t
        real, dimension(ele_loc(mesh, ele), ele_loc(mesh, ele)) :: mass_matrix
        
        integer i
        
        ! Compute detwei.
        call transform_to_physical(positions, ele, &
           ele_shape(infields(1), ele), dshape=dt_t, detwei=detwei)

        r = shape_dshape(ele_shape(mesh, ele), dt_t, detwei)
        do i=1, size(infields)
          
          if (compute(i)) then
            call addto(gradient(i), ele_nodes(mesh, ele), &
               tensormul(r, ele_val(infields(i), ele), 3) )
          end if

        end do
        
        ! Lump the mass matrix
        mass_matrix = shape_shape(ele_shape(mesh, ele), ele_shape(mesh, ele), detwei)
        call addto(lumped_mass_matrix, ele_nodes(mesh, ele), sum(mass_matrix, 2))
        
      end subroutine differentiate_field_ele

    end subroutine differentiate_field_lumped_multiple

    subroutine differentiate_field_lumped_single(infield, positions, derivatives, pardiff)
    !!< This routine computes the first derivatives using a weak finite element formulation.
      type(scalar_field), intent(in), target :: infield
      type(vector_field), intent(in) :: positions
      logical, dimension(:), intent(in) :: derivatives
      type(scalar_field), dimension(:), intent(inout) :: pardiff
        
      type(scalar_field), dimension(size(pardiff),1) :: pardiffs
        
      pardiffs(:,1)=pardiff
      
      call differentiate_field_lumped_multiple( (/ infield /), positions, derivatives, pardiffs)
      
    end subroutine differentiate_field_lumped_single

    subroutine differentiate_field_lumped_vector(infield, positions, outfield)
    !!< This routine computes the derivatives of a vector field returning a tensor field
      type(vector_field), intent(in), target :: infield
      type(vector_field), intent(in) :: positions
      type(tensor_field), intent(inout) :: outfield
        
      logical, dimension( positions%dim ):: derivatives
      type(scalar_field), dimension( infield%dim ):: infields
      type(scalar_field), dimension( positions%dim, infield%dim ):: pardiffs
      integer i, j
        
      derivatives=.true.
      do i=1, infield%dim
        infields(i)=extract_scalar_field(infield, i)
        do j=1, positions%dim
          pardiffs(j, i)=extract_scalar_field(outfield, j, i)
        end do
      end do
        
      call differentiate_field_lumped_multiple( infields, positions, derivatives, pardiffs)
      
    end subroutine differentiate_field_lumped_vector
      
    subroutine differentiate_field(infield, positions, derivatives, pardiff, bc_value, bc_type)
      type(scalar_field), intent(in), target :: infield
      type(vector_field), intent(in) :: positions
      logical, dimension(:), intent(in) :: derivatives
      type(scalar_field), dimension(:), intent(inout) :: pardiff

      ! weak bc's are rquired to calculated gradient of dg fields
      type(scalar_field), intent(in), optional :: bc_value
      integer, dimension(:), intent(in), optional :: bc_type

      integer :: i
      type(mesh_type), pointer :: mesh

      if (infield%field_type == FIELD_TYPE_CONSTANT) then
        do i=1,count(derivatives)
          call zero(pardiff(i))
        end do
        return
      end if

      if (continuity(infield)<0) then
        call differentiate_discontinuous_field(infield, positions, derivatives, pardiff, bc_value, bc_type)
        return
      end if

      if (continuity(pardiff(1))<0) then
        call differentiate_field_discontinuous(infield, positions, derivatives, pardiff)
        return
      end if

      mesh => infield%mesh
      call add_nelist(mesh)
      call differentiate_field_lumped_single(infield, positions, derivatives, pardiff)

      if (pseudo2d_coord /= 0) then
        if (derivatives(pseudo2d_coord)) then
        ! which pardiff corresponds to dimension pseudo2d_coord?
        i = count(derivatives(1:pseudo2d_coord))
        call zero(pardiff(i))
        end if
      end if
      
    end subroutine differentiate_field

    subroutine differentiate_discontinuous_field(infield, positions, derivatives, pardiff, bc_value, bc_type)
      ! calculated using:
      ! N_i N_j grad_u = N_i delta u_h - 
      !                  ({N_i} (u_h^-n^- + u_h^+n^+)) on internal faces -
      !                  (N_i (u_h - u_b) n) on weak dirichlet boundaries
      ! where: {x} = average of x over face
      !        u_h = value of u in element
      !        u_b = dirichlet boundary value
      ! (see Bassi et. al. 2005 - Discontinuous Galerkin solution of the Reynolds-averaged
      ! Navier–Stokes and k–x turbulence model equations, pg. 517

      type(scalar_field), intent(in), target :: infield
      type(vector_field), intent(in) :: positions
      logical, dimension(:), intent(in) :: derivatives
      type(scalar_field), dimension(:), intent(inout) :: pardiff

      type(scalar_field), intent(in) :: bc_value
      integer, dimension(:), intent(in) :: bc_type
      
      integer :: ele, i

      if (infield%field_type == FIELD_TYPE_CONSTANT) then
        do i=1,count(derivatives)
          if (derivatives(i)) then
             call zero(pardiff(i))
          end if
        end do
        return
      end if

      ! only works if all pardiff fields are discontinuous:
      do i=1, count(derivatives)
        assert(pardiff(i)%mesh%continuity<0)
      end do
      
      ! calculate gradient
      do ele = 1, ele_count(infield)
        call calculate_grad_ele_dg(infield, positions, derivatives, pardiff, ele, bc_value, bc_type)
      end do

    end subroutine differentiate_discontinuous_field    

    subroutine calculate_grad_ele_dg(infield, positions, derivatives, pardiff, ele, bc_value, bc_type)
      type(scalar_field), intent(in), target :: infield
      type(vector_field), intent(in) :: positions
      logical, dimension(:), intent(in) :: derivatives
      type(scalar_field), dimension(:), intent(inout) :: pardiff  
      integer, intent(in) :: ele
      type(scalar_field), intent(in) :: bc_value
      integer, dimension(:), intent(in) :: bc_type    
      
      ! variables for interior integral
      type(element_type), pointer :: shape
      real, dimension(ele_loc(infield, ele), ele_ngi(infield, ele), positions%dim) :: dshape
      real, dimension(ele_ngi(infield, ele)) :: detwei
      real, dimension(positions%dim, ele_ngi(infield, ele)) :: grad_h_gi
      real, dimension(positions%dim, ele_loc(infield, ele)) :: rhs

      ! variables for surface integral
      integer :: ni, ele_2, face, face_2, i
      integer, dimension(:), pointer :: neigh

      ! inverse mass
      real, dimension(ele_loc(infield, ele), ele_loc(infield, ele)) :: inv_mass

      ! In parallel, we only construct the equations on elements we own, or
      ! those in the L1 halo.
      if (.not.(element_owned(infield, ele).or.element_neighbour_owned(infield, ele))) then
        return
      end if

      shape => ele_shape(infield, ele) 
      call transform_to_physical(positions, ele, shape, dshape, detwei)

      ! Calculate grad within the element
      grad_h_gi = ele_grad_at_quad(infield, ele, dshape)

      ! Assemble interior contributions to rhs
      rhs = shape_vector_rhs(shape, grad_h_gi, detwei)

      ! Interface integrals
      neigh=>ele_neigh(infield, ele)
      do ni=1,size(neigh)
        ! Find the relevant faces.
        ele_2 = neigh(ni)
        face = ele_face(infield, ele, ele_2)

        if (ele_2>0) then
          ! Internal faces.
          face_2=ele_face(infield, ele_2, ele)
        else
          ! External face.
          face_2=face
        end if

        call calculate_grad_ele_dg_interface(ele, face, face_2, ni, &
             & rhs, positions, infield, bc_value, bc_type)
      end do

      ! multiply by inverse of mass matrix
      inv_mass = inverse(shape_shape(shape, shape, detwei))
      do i = 1, positions%dim
        rhs(i,:) = matmul(inv_mass, rhs(i,:))
        if (derivatives(i)) then
          call set(pardiff(i), ele_nodes(pardiff(i), ele), rhs(i,:))
        end if
      end do

    end subroutine calculate_grad_ele_dg

    subroutine calculate_grad_ele_dg_interface(ele, face, face_2, &
         ni, rhs, positions, infield, bc_value, bc_type)

      !!< Construct the DG element boundary integrals on the ni-th face of
      !!< element ele.
      integer, intent(in) :: ele, face, face_2, ni
      type(scalar_field), intent(in) :: bc_value, infield
      type(vector_field), intent(in) :: positions
      integer, dimension(:), intent(in) :: bc_type
      real, dimension(positions%dim, ele_loc(infield, ele)), intent(inout) :: rhs

      ! Face objects and numberings.
      type(element_type), pointer :: shape
      real, dimension(positions%dim, face_ngi(infield, face)) :: normal
      real, dimension(face_ngi(infield, face)) :: detwei, in_q, in_q_2, in_bc_q
      real, dimension(positions%dim, face_ngi(infield, face)) :: vector
      real, dimension(positions%dim, face_loc(infield, face)) :: face_rhs
      real, dimension(ele_loc(infield, ele)) :: elenodes
      real, dimension(face_loc(infield, face)) :: facenodes

      integer :: i, j

      face_rhs = 0.0
      vector = 0.0

      ! shape and detwei are the same for both faces, normal+ = - normal-
      shape => face_shape(infield, face)
      call transform_facet_to_physical(positions, face, detwei_f=detwei, normal=normal)

      if (face==face_2) then  
        ! boundary faces - need to apply weak dirichlet bc's
        ! = - int_ v_h \cdot (u - u^b) n 
        ! first check for weak-dirichlet bc
        if (bc_type(face) == 1) then  
          in_q = face_val_at_quad(infield, face)
          in_bc_q = ele_val_at_quad(bc_value, face)

          do i=1, mesh_dim(infield)
            vector(i,:) = -1.0*(in_q(:) - in_bc_q(:))*normal(i,:)
          end do
          face_rhs = shape_vector_rhs(shape, vector, detwei) 
        end if
      else    
        ! internal face
        ! = int_ {v_h} \cdot J(x)  
        in_q = face_val_at_quad(infield, face)
        in_q_2 = face_val_at_quad(infield, face_2)

        do i=1, mesh_dim(infield)
          vector(i,:) = -0.5*(in_q(:) - in_q_2(:))*normal(i,:)
        end do
        face_rhs = shape_vector_rhs(shape, vector, detwei) 
      end if

      elenodes = ele_nodes(infield, ele)
      facenodes = face_global_nodes(infield, face)
      do i=1, face_loc(infield, face)
        do j=1, ele_loc(infield, face)
          if (facenodes(i) == elenodes(j)) then
            rhs(:, j) = rhs(:, j) + face_rhs(:, i)
          end if
        end do
      end do

    end subroutine calculate_grad_ele_dg_interface

    subroutine differentiate_field_discontinuous(infield, positions, derivatives, pardiff)
      type(scalar_field), intent(in), target :: infield
      type(vector_field), intent(in) :: positions
      logical, dimension(:), intent(in) :: derivatives
      type(scalar_field), dimension(:), intent(inout) :: pardiff
        
      type(element_type) xshape, inshape, dershape
      real, dimension(ele_loc(infield,1), ele_ngi(infield,1), size(derivatives)):: dinshape
      real, dimension(size(derivatives), size(derivatives), ele_ngi(infield,1)):: invJ
      real, dimension(ele_loc(pardiff(1),1)):: r
      real, dimension(size(r), size(r)):: M
      real, dimension(size(derivatives), size(r), ele_loc(infield,1)):: Q
      real, dimension(ele_ngi(infield,1)):: detwei
      integer ele, gi, i, j
      
      if (infield%field_type == FIELD_TYPE_CONSTANT) then
        do i=1,count(derivatives)
          if (derivatives(i)) then
             call zero(pardiff(i))
          end if
        end do
        return
      end if
      
      ! only works if all pardiff fields are discontinuous:
      do i=1, count(derivatives)
        assert(pardiff(i)%mesh%continuity<0)
      end do
      ! and the infield is continuous!!!!
      if (infield%mesh%continuity<0) then
        ewrite(0,*) "If the following error is directly due to user input"
        ewrite(0,*) "a check and a more helpful error message should be inserted in"
        ewrite(0,*) "the calling routine (outside field_derivatives) - please mantis this:"
        ewrite(0,*) "Error has occured in differentiate_field_discontinuous, with field, ", trim(infield%name)
        FLAbort("Shouldn't get here?")
      end if
      
      xshape=ele_shape(positions, 1)
      inshape=ele_shape(infield, 1)
      dershape=ele_shape(pardiff(1), 1)
      
      do ele=1, element_count(infield)
        
         ! calculate the transformed derivative of the shape function
         call compute_inverse_jacobian( positions, ele, invJ, detwei=detwei)
         do gi=1, inshape%ngi
            do i=1, inshape%loc
               dinshape(i,gi,:)=matmul(invJ(:,:,gi), inshape%dn(i,gi,:))
            end do
         end do
           
         M=shape_shape(dershape, dershape, detwei)
         Q=shape_dshape(dershape, dinshape, detwei)
         call invert(M)
         
         ! apply Galerkin projection M^{-1} Q \phi
         j=0
         do i=1, size(derivatives)
            if (derivatives(i)) then
               j=j+1
               r=matmul(M, matmul(Q(i,:,:), ele_val(infield, ele)))
               call set(pardiff(j), ele_nodes(pardiff(j), ele), r)
            end if
        end do
          
      end do
      
    end subroutine differentiate_field_discontinuous
      
    subroutine compute_hessian(infield, positions, hessian)
      type(scalar_field), intent(inout) :: infield
      type(vector_field), intent(in) :: positions
      type(tensor_field), intent(inout) :: hessian

      integer :: node

      if (infield%field_type == FIELD_TYPE_CONSTANT) then
        call zero(hessian)
        return
      end if

      call add_nelist(infield%mesh)
      call compute_hessian_real(infield, positions, hessian)

      if (pseudo2d_coord /= 0) then
        do node=1,node_count(hessian)
          hessian%val(pseudo2d_coord, :, node) = 0.0
          hessian%val(:, pseudo2d_coord, node) = 0.0
        end do
      end if
    end subroutine compute_hessian

    subroutine curl(infield, positions, curl_norm, curl_field)
      type(vector_field), intent(in) :: positions, infield
      type(scalar_field), intent(inout), optional :: curl_norm ! norm of curl_field
      type(vector_field), intent(inout), optional :: curl_field

      type(vector_field), dimension(positions%dim) :: grad_v
      integer :: i
      real :: w, a, b, c
      type(mesh_type) :: mesh
      
      assert(positions%dim == 3)

      mesh = infield%mesh
      call add_nelist(mesh)

      do i=1,positions%dim
        call allocate(grad_v(i), positions%dim, infield%mesh, "Grad V")
        call grad(extract_scalar_field(infield, i), positions, grad_v(i))
      end do

      if (present(curl_field)) then
        call zero(curl_field)
      end if

      if (present(curl_norm)) then
        call zero(curl_norm)
      end if

      do i=1,node_count(infield)
        a = grad_v(3)%val(2,i) - grad_v(2)%val(3,i) ! dw/dy - dv/dz
        b = grad_v(1)%val(3,i) - grad_v(3)%val(1,i) ! du/dz - dw/dx
        c = grad_v(2)%val(1,i) - grad_v(1)%val(2,i) ! dv/dx - du/dy
        if (present(curl_norm)) then
          w = sqrt(a**2 + b**2 + c**2)
          call addto(curl_norm, i, w)
        end if
        if (present(curl_field)) then
          call addto(curl_field, i, (/a, b, c/))
        end if
      end do

      do i=1,positions%dim
        call deallocate(grad_v(i))
      end do

    end subroutine curl
    
  subroutine u_dot_nabla_scalar(v_field, in_field, positions, out_field)
    !!< Calculates (u dot nabla) in_field for scalar fields
      
    type(vector_field), intent(in) :: v_field
    type(scalar_field), intent(in) :: in_field
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(inout) :: out_field
    
    integer :: i
    real, dimension(positions%dim) :: grad_val_at_node, &
      & v_field_val_at_node
    type(vector_field) :: gradient
    
    call allocate(gradient, positions%dim, in_field%mesh, "Gradient")
    
    call grad(in_field, positions, gradient)
    
    call zero(out_field)
    do i = 1, node_count(out_field)
      grad_val_at_node = node_val(gradient, i)
      v_field_val_at_node = node_val(v_field, i)
      call set(out_field, i, &
        & dot_product(v_field_val_at_node, grad_val_at_node))
    end do
    
    call deallocate(gradient)
    
  end subroutine u_dot_nabla_scalar
  
  subroutine u_dot_nabla_vector(v_field, in_field, positions, out_field)
    !!< Calculates (u dot nabla) in_field for vector fields
    
    type(vector_field), intent(in) :: v_field
    type(vector_field), intent(in) :: in_field
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(inout) :: out_field

    integer :: i
    type(scalar_field) :: out_field_comp

    do i = 1, v_field%dim
      out_field_comp = extract_scalar_field(out_field, i)
      call u_dot_nabla(v_field, &
        & extract_scalar_field(in_field, i), positions, &
        & out_field_comp)
    end do
    
  end subroutine u_dot_nabla_vector

!    subroutine compute_hessian_eqf(infield, positions, hessian)
!      !!< This routine computes the hessian using the method suggested
!      !!< to me by Chris Pain. Not published yet. Let's see if it works, first.
!      !!< The idea is, for a given element, construct a quadratic polynomial
!      !!< expansion of the field over that element by constraining it as follows:
!      !!< eqns. 1-4: the value at the node is the same
!      !!< eqns. 5-11: the recovered gradient in the direction of the element centroid is the same
!      type(scalar_field), intent(in) :: infield
!      type(vector_field), intent(in) :: positions
!      type(tensor_field), intent(inout) :: hessian
!
!      !! Linear tets only for now. Sorry, David.
!
!      integer :: i, j
!      real, dimension(MATRIX_SIZE_QF, MATRIX_SIZE_QF) :: A
!      real, dimension(MATRIX_SIZE_QF) :: b
!    end subroutine compute_hessian_eqf

    function get_quadratic_fit_eqf(infield, positions, ele, gradient) result(b)
      !!< Implement Chris' idea for computing the Hessian.
      !!< The idea is, for a given element, construct a quadratic polynomial
      !!< expansion of the field over that element by constraining it as follows:
      !!< eqns. 1-4: the value at the node is the same
      !!< eqns. 5-11: for each edge, the value of the gradient in the direction of the centroid is the same
      !!< Linear tets only.
      type(scalar_field), intent(in) :: infield
      type(vector_field), intent(in) :: positions
      integer, intent(in) :: ele
      real, dimension(ele_loc(infield, ele), positions%dim) :: gradient
      real, dimension(MATRIX_SIZE_QF) :: b

      real, dimension(MATRIX_SIZE_QF, MATRIX_SIZE_QF) :: A
      integer :: i, j, k
      integer, dimension(:), pointer :: nodelist
      real, dimension(positions%dim) :: centroid, edge_centre, dir

      ! i is the index of the equation we're writing down.
      i = 1
      nodelist => ele_nodes(infield, ele)
      centroid = 0.0
      A = 0.0

!      do j=1,4
!        centroid = centroid + node_val(positions, nodelist(j))
!      end do
!      centroid = centroid / 4.0
      centroid = insphere_tet(ele_val(positions, ele))

      do j=1,4
        A(i, :) = getP_qf(node_val(positions, nodelist(j)))
        b(i) = node_val(infield, nodelist(j))
        i = i + 1
      end do

      do j=1,4
        do k=j+1,4
          edge_centre = (node_val(positions, nodelist(j)) + node_val(positions, nodelist(k))) / 2.0
          dir = centroid - edge_centre
          !dir = node_val(positions, nodelist(k)) - node_val(positions, nodelist(j))
          A(i, 2) = dir(1); A(i, 3) = dir(2); A(i, 4) = dir(3)
          A(i, 5) = 2 * dir(1) * edge_centre(1); A(i, 6) = 2 * dir(2) * edge_centre(2); A(i, 7) = 2 * dir(3) * edge_centre(3);
          A(i, 8) = dir(1) * edge_centre(2) + dir(2) * edge_centre(1)
          A(i, 9) = dir(3) * edge_centre(1) + dir(1) * edge_centre(3)
          A(i, 10) = dir(2) * edge_centre(3) + dir(3) * edge_centre(2)
          b(i) = dot_product((gradient(j, :) + gradient(k, :)) / 2.0, dir)
          i = i + 1
        end do
      end do
      A(11, 11) = 1.0; b(11) = 0.0

      call solve(A, b)
    end function get_quadratic_fit_eqf
    
    subroutine compute_hessian_eqf(infield, positions, hessian)
      type(scalar_field), intent(in) :: infield
      type(vector_field), intent(in) :: positions
      type(tensor_field), intent(inout) :: hessian

      integer :: ele, j, node
      integer, dimension(:), pointer :: nodelist
      type(vector_field) :: gradient
      integer, dimension(node_count(hessian)) :: touched
      real, dimension(MATRIX_SIZE_QF) :: fit
      real, dimension(positions%dim, positions%dim) :: tmp_hessian
      real :: x, y, z
      type(element_type) :: t_shape, x_shape

      touched = 0
      call allocate(gradient, positions%dim, infield%mesh, "Gradient")
      call grad(infield, positions, gradient)
      call zero(hessian)
      t_shape = ele_shape(infield, 1)
      x_shape = ele_shape(positions, 1)

      do ele=1,ele_count(infield)
        nodelist => ele_nodes(infield, ele)
        fit = get_quadratic_fit_eqf(infield, positions, ele, transpose(ele_val(gradient, ele)))
        do j=1,ele_loc(infield, ele)
          node = nodelist(j)
          x = node_val(positions, 1, node); y = node_val(positions, 2, node); z = node_val(positions, 3, node);
          tmp_hessian(1, 1) = 2 * fit(5)
          tmp_hessian(1, 2) = fit(8) + z * fit(11)
          tmp_hessian(1, 3) = fit(9) + y * fit(11)
          tmp_hessian(2, 1) = fit(8) + z * fit(11)
          tmp_hessian(2, 2) = 2 * fit(6)
          tmp_hessian(2, 3) = fit(10) + x * fit(11)
          tmp_hessian(3, 1) = fit(9) + y * fit(11)
          tmp_hessian(3, 2) = fit(10) + x * fit(11)
          tmp_hessian(3, 3) = 2 * fit(7)
          hessian%val(:, :, node) = hessian%val(:, :, node) + tmp_hessian
          touched(node) = touched(node) + 1
        end do
      end do

      do node=1,node_count(hessian)
        hessian%val(:, :, node) = hessian%val(:, :, node) / touched(node)
      end do

      call hessian_boundary_correction(hessian, positions, t_shape)
      call hessian_squash_pseudo2d(hessian)

      call deallocate(gradient)
    end subroutine compute_hessian_eqf

    subroutine div(infield, positions, divergence)
      !! Implement div() operator.
      type(vector_field), intent(in):: infield, positions
      type(scalar_field), intent(inout), target  :: divergence

      type(scalar_field) :: component
      type(scalar_field), dimension(1) :: derivative
      type(mesh_type), pointer :: mesh
      logical, dimension(mesh_dim(infield)) :: derivatives
      integer :: i

      mesh => divergence%mesh
      call allocate(derivative(1), mesh, "Derivative")

      call zero(divergence)
      derivatives = .false.

      do i=1,mesh_dim(infield)
        derivatives(i) = .true.
        component = extract_scalar_field(infield, i)
        call differentiate_field(component, positions, derivatives, derivative)
        call addto(divergence, derivative(1))
        derivatives = .false.
      end do

      call deallocate(derivative(1))
    end subroutine div
    
    function insphere_tet(positions) result(centre)
      !! dim x loc
      real, dimension(3, 4), intent(in) :: positions
      real, dimension(size(positions, 1)) :: centre

      real, dimension(size(positions, 1)) :: u, v, w, p, q, r, O1, O2, y, s
      real :: t

      u = positions(:, 2) - positions(:, 1)
      v = positions(:, 3) - positions(:, 1)
      w = positions(:, 4) - positions(:, 1)
      p = cross_product(u, v); p = p / norm2(p)
      q = cross_product(v, w); q = q / norm2(q)
      r = cross_product(w, u); r = r / norm2(r)

      O1 = p -q
      O2 = q - r
      y = cross_product(O1, O2)

      O1 = u - w
      O2 = v - w
      s = cross_product(O1, O2); s = -1 * (s / norm2(s))

      O1 = s - p
      t = dot_product(w, s) / dot_product(y, O1)
      centre = positions(:, 1) + t * y
    end function insphere_tet

    recursive function get_cubic_fit_cf(infield, positions, node, level, maxlevel) result(b)
      !!< Fit a cubic function to infield around the node,
      !!< with a least squares approach.
      type(scalar_field), intent(in) :: infield
      type(vector_field), intent(in) :: positions
      integer :: node
      integer, optional :: level, maxlevel
      real, dimension(MATRIX_SIZE_CF) :: b

      type(patch_type) :: node_patch
      real, dimension(MATRIX_SIZE_CF, MATRIX_SIZE_CF) :: A
      real, dimension(MATRIX_SIZE_CF, 1) :: b_tmp
      real, dimension(MATRIX_SIZE_CF_2D, MATRIX_SIZE_CF_2D) :: A_2D
      real, dimension(MATRIX_SIZE_CF_2D, 1) :: b_tmp_2D
      integer :: i, nnode, stat, llevel, lmaxlevel
      type(mesh_type) :: mesh

      if (present(level)) then
        llevel = level
      else
        llevel = 2
      end if

      if (present(maxlevel)) then
        lmaxlevel = maxlevel
      else
        lmaxlevel = llevel + 1
      end if

      mesh = infield%mesh
      node_patch = get_patch_node(mesh, node, level=llevel, min_nodes=min(node_count(mesh), int(2.5 * MATRIX_SIZE_CF)))
      A = 0.0; b_tmp = 0.0
      do i=1,node_patch%count
        nnode = node_patch%elements(i)

        A = A + compute_matrix_contribution_cf(node_val(positions, nnode))
        b_tmp(:, 1) = b_tmp(:, 1) + compute_rhs_contribution_cf(node_val(positions, nnode), node_val(infield, nnode))
      end do

      if (pseudo2d_coord /= 0) then
        if (pseudo2d_coord == 1) then
          A_2D = A(CF_2D_X, CF_2D_X)
          b_tmp_2D(:, 1) = b_tmp(CF_2D_X, 1)
          call solve(A_2D, b_tmp_2D, stat)
          b = 0.0; b(CF_2D_X) = b_tmp_2D(:, 1)
        else if (pseudo2d_coord == 2) then
          A_2D = A(CF_2D_Y, CF_2D_Y)
          b_tmp_2D(:, 1) = b_tmp(CF_2D_Y, 1)
          call solve(A_2D, b_tmp_2D, stat)
          b = 0.0; b(CF_2D_Y) = b_tmp_2D(:, 1)
        else if (pseudo2d_coord == 3) then
          A_2D = A(CF_2D_Z, CF_2D_Z)
          b_tmp_2D(:, 1) = b_tmp(CF_2D_Z, 1)
          call solve(A_2D, b_tmp_2D, stat)
          b = 0.0; b(CF_2D_Z) = b_tmp_2D(:, 1)
        end if
      else
        call solve(A, b_tmp, stat)
        b = b_tmp(:, 1)
      end if


      if (llevel < lmaxlevel) then
        if (stat /= 0) then
          ! If the solver fails, go one more level deep
          ! to get enough equations to do the fitting
          b = get_cubic_fit_cf(infield, positions, node, level=llevel+1)
        end if
      end if

      deallocate(node_patch%elements) ! Have to do it here; otherwise memory leak
    end function get_cubic_fit_cf
end module field_derivatives
