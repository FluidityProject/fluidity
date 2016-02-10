#include "fdebug.h"
!#define EXTRA_SPECIAL_GRADATION_DEBUGGING

module gradation_metric
!!< This module implements a gradation algorithm
!!< to ensure smooth mesh transitions from small
!!< to large edge lengths.
!!< The algorithm implemented is a modification of
!!< "Anisotropic mesh gradation control", Li et. al,
!!< 13th International Meshing Roundtable, 2004

  use unittest_tools
  use adjacency_lists
  use linked_lists
  use metric_tools
  use fields
  use vtk_interfaces
  use node_boundary
  use edge_length_module
  use field_derivatives

  implicit none

  ! These are the DEFAULTS ONLY IF YOU DON'T CALL
  ! initialise_gradation_metric.
  ! initialise_gradation_metric changes them for real code.
  ! Bottom line: if you want to change whether gradation is used
  ! or the gradation constant, CHANGE THE VALUES in INITIALISE_GRADATION_METRIC
  logical :: use_gradation_metric = .false.
  logical :: gradation_initialised = .false.
  real :: gamma0 = 1.5 !!< Gamma is a measure of the smoothness of the transition
                       !!< an edge. Gamma0 is the maximum allowed value for gamma.

  real :: theta1 = 0.01 !!< Theta1 is the maximum angle we deem to be noise in the 
                        !!< misalignment of the metric tensors at distance 1.0.
  real :: theta0 = 0.35 !!< Theta0 is the maximum angle we deem to be noise in the
                        !!< misalignment of the metric tensors at distance 0.0.
  integer :: max_rot_its = 8
  real :: domain_scale

  contains

  subroutine initialise_gradation_metric
  
    use_gradation_metric=have_option("/mesh_adaptivity/hr_adaptivity/enable_gradation")

    if (have_option("/mesh_adaptivity/hr_adaptivity/enable_gradation/gradation_parameter")) then
      call get_option("/mesh_adaptivity/hr_adaptivity/enable_gradation/gradation_parameter", gamma0)
    else
      gamma0 = 1.5
    end if

    ewrite(2,*) 'gradation: ', use_gradation_metric, gamma0
    
  end subroutine initialise_gradation_metric

  subroutine form_gradation_metric(positions, error_metric, noits)
    type(tensor_field), intent(inout) :: error_metric !!< The metric formed so far
    type(vector_field), intent(in) :: positions
    integer, intent(out), optional :: noits

    type(csr_matrix) :: nnlist !!< Node-node adjacency list
    type(csr_sparsity), pointer :: nn_sparsity !!< Node-node adjacency list sparsity
    type(elist) :: edgelist    !!< Linked list of edges
    type(mesh_type) :: mesh
    integer :: p, q !! the nodes
    real, dimension(error_metric%dim(1), error_metric%dim(2)) :: vec_P, vec_Q ! eigenvectors
    real, dimension(error_metric%dim(1)) :: val_P, val_Q ! the eigenvalues
    logical :: vals_changed_P, vals_changed_Q !!< have P or Q changed? If so need to reform
    logical :: vecs_changed_P, vecs_changed_Q !!< the metric and update any surrounding nodes.

    integer :: boundcount_P, boundcount_Q, expected_boundcount
    logical :: do_warp_directions

    integer :: dim, count

    integer :: global_its, end_marker ! count how many sweeps this involves

    type(scalar_field) :: edgelen
    integer, save :: adaptcnt = 0

#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
    integer :: stepcount
    type(scalar_field) :: nodefield
#endif
    logical :: debug_metric

    dim = error_metric%dim(1)
    debug_metric = have_option("/mesh_adaptivity/hr_adaptivity/debug/write_metric_stages")

    ewrite(2,*) "++: Applying gradation"

    mesh = error_metric%mesh
    domain_scale = domain_length_scale(positions)
    call initialise_boundcount(mesh, positions)
#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
    write(0,*) "domain_scale == ", domain_scale
#endif

    if (domain_is_2d()) then
      expected_boundcount = 1
    else
      expected_boundcount = 0
    end if

    !! Here I describe a convention I use in the nnlist.
    !! Normally the val array of the CSR matrix doesn't exist
    !! (as the value doesn't really matter). Here I set
    !! (i,j) to < 0 if (i,j) is NOT in the linked list of edges,
    !! and > 0 if it IS in the linked list of edges.
    !! |(i,j)| is the number of times it's been checked + 1.

    nn_sparsity => extract_nnlist(mesh)
    call allocate(nnlist, nn_sparsity, type=CSR_INTEGER)
    nnlist%ival = -1

    call construct_edge_list(mesh, nnlist, edgelist) 

    if (debug_metric) then
      call allocate(edgelen, error_metric%mesh, "Desired edge lengths")
    end if

    ! OK. So now we have the edge list.
#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
    stepcount = 0
    call allocate(nodefield, error_metric%mesh, "Node number")
    call get_node_field(error_metric%mesh, nodefield)
#endif

    end_marker = edgelist%length
    global_its = 0

    do while (edgelist%length /= 0)

      ! Count the number of sweeps through the mesh.
      end_marker = end_marker - 1
      if (end_marker == 0) then
        global_its = global_its + 1
        end_marker = edgelist%length
      end if

      call wrap_pop(nnlist, edgelist, p, q, count)                               ! fetch the nodes
      
#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
      write(0,*) "----------------------------------------------------"
      write(0,*) "stepcount == ", stepcount
      write(0,*) "(p, q) == (", p, ", ", q, ")" 
#endif

      call eigendecomposition_symmetric(node_val(error_metric, p), vec_P, val_P) ! decompose
      call eigendecomposition_symmetric(node_val(error_metric, q), vec_Q, val_Q)

#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
      call check_basis(vec_P)
      call check_basis(vec_Q)
      write(0,*) "Input P:"
      call write_matrix(node_val(error_metric, p), "P")
      call write_matrix(vec_P, "vec_P")
      call write_vector(val_P, "val_P")
      write(0,*) "Input Q:"
      call write_matrix(node_val(error_metric, q), "Q")
      call write_matrix(vec_Q, "vec_Q")
      call write_vector(val_Q, "val_Q")
#endif

      vecs_changed_P = .false. ; vecs_changed_Q = .false.

      boundcount_P = node_boundary_count(p)
      boundcount_Q = node_boundary_count(q)

       do_warp_directions = .false.

#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
      write(0,*) "boundcount_P == ", boundcount_P
      write(0,*) "boundcount_Q == ", boundcount_Q
      write(0,*) "do_warp_directions == ", do_warp_directions
#endif

      if (do_warp_directions) then
        call warp_directions(vec_P, val_P, vecs_changed_P, vec_Q, val_Q, vecs_changed_Q, distance(positions, p, q))
      end if

#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
      write(0,*) "P after warping:"
      call write_matrix(vec_P, "vec_P")
      write(0,*) "Q after warping:"
      call write_matrix(vec_Q, "vec_Q")
      write(0,*) "vecs_changed_P == ", vecs_changed_P, "vecs_changed_Q == ", vecs_changed_Q
      call check_basis(vec_P)
      call check_basis(vec_Q)
#endif


      vals_changed_P = .false. ; vals_changed_Q = .false.
      call reduce_edgelen(node_val(error_metric, p), vec_P, val_P, vals_changed_P, &
                          node_val(error_metric, q), vec_Q, val_Q, vals_changed_Q, &
                          positions, p, q)

#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
      write(0,*) "P after reducing:"
      call write_vector(val_P, "val_P")
      write(0,*) "Q after reducing:"
      call write_vector(val_Q, "val_Q")
      write(0,*) "vals_changed_P == ", vals_changed_P, "vals_changed_Q == ", vals_changed_Q
#endif

      if (vals_changed_P) then
        call eigenrecomposition(error_metric%val(:, :, p), vec_P, val_P)
        call tag_edges(nnlist, edgelist, p, q, count)
      end if
      if (vals_changed_Q) then
        call eigenrecomposition(error_metric%val(:, :, q), vec_Q, val_Q)
        call tag_edges(nnlist, edgelist, q, p, count)
      end if
      if (count <= max_rot_its) then ! honour directional changes for the first 4 sweeps
        if (vecs_changed_P) then
          call eigenrecomposition(error_metric%val(:, :, p), vec_P, val_P)
          call tag_edges(nnlist, edgelist, p, q, count)
        end if
        if (vecs_changed_Q) then
          call eigenrecomposition(error_metric%val(:, :, q), vec_Q, val_Q)
          call tag_edges(nnlist, edgelist, q, p, count)
        end if
      end if

#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
      call check_metric(error_metric)
      call get_edge_lengths(error_metric, edgelen)
      call vtk_write_fields("data/gradation_debug", stepcount, positions, positions%mesh, &
                            sfields=(/nodefield, edgelen/), tfields=(/error_metric/))
      stepcount = stepcount + 1
#endif
    end do

    call deallocate(nnlist)

    if (debug_metric) then
      call get_edge_lengths(error_metric, edgelen)
      call vtk_write_fields(trim("gradation_metric"), adaptcnt, positions, positions%mesh, &
                             sfields=(/edgelen/), tfields=(/error_metric/))
      call deallocate(edgelen)
      adaptcnt = adaptcnt + 1
    endif
#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
    call deallocate(nodefield)
#endif

  ewrite(2,*) "Finished gradation algorithm: global iterations == ", global_its
  if (present(noits)) noits = global_its
  end subroutine form_gradation_metric

  subroutine wrap_pop(nnlist, edgelist, p, q, count)
    !!< Wrap the pop() routine of the linked list. I want
    !!< to tag in nnlist that that edge is not in the list anymore.
    type(csr_matrix), intent(inout) :: nnlist
    type(elist), intent(inout) :: edgelist
    integer, intent(out) :: p, q, count

    call spop(edgelist, p, q)
    count = ival(nnlist, p, q)
    call set(nnlist, p, q, -1 * count)
    call set(nnlist, q, p, -1 * count)
  end subroutine

  subroutine tag_edges(nnlist, edgelist, p, q, count)
    !!< Insert all nodes connected to p into the edge list,
    !!< except the one you just checked,
    !!< if they're not already there.
    type(csr_matrix), intent(inout) :: nnlist
    type(elist), intent(inout) :: edgelist
    integer, intent(in) :: p, q, count

    integer :: j
    integer, dimension(:), pointer :: cola

    cola => row_m_ptr(nnlist, p)
    !write(0,*) "Tagging all edges connected to ", p, "that are not", q, "."

    do j=1,row_length(nnlist,p)
      if (cola(j) == q) cycle
      call add_to_edge_list(nnlist, edgelist, p, cola(j), count)
    end do
  end subroutine tag_edges

  subroutine construct_edge_list(mesh, nnlist, edgelist)
   !! From the node-node adjacency list, construct a linked
   !! list of edges.
   type(mesh_type), intent(in) :: mesh
   type(csr_matrix), intent(inout) :: nnlist
   type(elist), intent(out) :: edgelist

   integer :: i,j, rowlen
   integer, dimension(:), pointer :: cola

   do i=1,mesh%nodes
     rowlen = row_length(nnlist, i)
     cola => row_m_ptr(nnlist, i)
     do j=1,rowlen
       call add_to_edge_list(nnlist, edgelist, i, cola(j), 0)
     end do
   end do

 end subroutine construct_edge_list

 subroutine add_to_edge_list(nnlist, edgelist, i, j, count)
   type(csr_matrix), intent(inout) :: nnlist
   type(elist), intent(inout) :: edgelist
   integer, intent(in) :: i, j, count
   
   if (i == j) return
   if (ival(nnlist, i, j) < 0) then
     call insert(edgelist, i, j)
     call set(nnlist, i, j, count + 1)
     call set(nnlist, j, i, count + 1)
   end if
  end subroutine add_to_edge_list

  subroutine compute_omega(vec_P, perm_P, vec_Q, perm_Q, dist, idx, omega, angle)
    !!< Compute omega, the maximum angle we deem to be noise between two
    !!< nodes separated by dist.
    real, dimension(:, :), intent(in) :: vec_P, vec_Q
    integer, dimension(size(vec_P, 1)), intent(in) :: perm_P, perm_Q
    real, intent(in) :: dist
    integer, intent(out) :: idx
    real, intent(out) :: omega
    real, intent(out) :: angle

    integer :: i, dim
    real :: curangle, maxangle

    idx = 0
    maxangle = 0.0
    dim = size(perm_P)

    do i=1,dim
      curangle = get_angle(vec_P(:, perm_Q(i)), vec_Q(:, perm_P(i)))
      if (curangle > maxangle) then
        maxangle = curangle
        idx = i
      end if
    end do

    angle = maxangle

    omega = theta0 + (theta1 - theta0) * dist / domain_scale
    if (omega < 0.0) omega = 0.0
  end subroutine compute_omega

  subroutine warp_directions(vec_P, val_P, changed_P, vec_Q, val_Q, changed_Q, dist)
    !!< Given two metric tensors,
    !!< modify the directions of their eigenvectors to be more aligned,
    !!< while respecting anisotropism.
 
    real, dimension(:, :), intent(inout) :: vec_P, vec_Q
    real, dimension(size(vec_P, 1)), intent(in) :: val_P, val_Q
    logical, intent(inout) :: changed_P, changed_Q
    real, intent(in) :: dist
 
    integer, dimension(size(vec_P, 1)), target :: perm_P, perm_Q     ! contains the permutation matching up vectors
                                                                     ! which old ones correspond to the new ones
    real :: omega, angle, angle_P, angle_Q
    integer :: dim, idx, i

    dim = size(vec_P, 1)
  
    !write(0,*) "Starting warp_directions"

    if (metric_isotropic(val_P) .and. metric_isotropic(val_Q)) then
      if (val_P(1) >= val_Q(1)) then 
        vec_Q = vec_P
      else 
        vec_P = vec_Q
      end if
      return
    end if

    if (metric_isotropic(val_P)) then
      vec_P = vec_Q
      return
    end if

    if (metric_isotropic(val_Q)) then
      vec_Q = vec_P
      return
    end if

    if (vec_P .feq. vec_Q) return

    do i=1,dim
      call match_up_ellipsoids(vec_P, val_P, perm_P, vec_Q, val_Q, perm_Q)
      call compute_omega(vec_P, perm_P, vec_Q, perm_Q, dist, idx, omega, angle)
#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
      call write_matrix(vec_P, "vec_P")
      call write_matrix(vec_Q, "vec_Q")
      write(0,*) "perm_P == ", perm_P
      write(0,*) "perm_Q == ", perm_Q
      write(0,*) "omega == ", omega
      write(0,*) "angle == ", angle
      write(0,*) "dist == ", dist
#endif
      if (angle < omega) then
        if (angle < 0.01) cycle
        angle_P = (aspect_ratio(val_Q) / (aspect_ratio(val_Q) + aspect_ratio(val_P)))  * angle
        angle_Q = (aspect_ratio(val_P) / (aspect_ratio(val_Q) + aspect_ratio(val_P)))  * angle
        call rotate_vec(vec_P, perm_P, vec_Q, perm_Q, idx, angle_P)
        call rotate_vec(vec_Q, perm_Q, vec_P, perm_P, idx, angle_Q)
        changed_P = .true. ; changed_Q = .true.
      end if
    end do
  end subroutine warp_directions

  subroutine rotate_vec(vec_A, perm_A, vec_B, perm_B, idx, angle)
    real, dimension(:, :), intent(inout) :: vec_A, vec_B
    integer, dimension(size(vec_A, 1)), intent(in) :: perm_A, perm_B
    real, intent(in) :: angle
    integer, intent(in) :: idx
    real, dimension(size(vec_A, 1), size(vec_A, 1)) :: mat
    real, dimension(size(vec_A, 1)) :: cross, tmpvec
    integer :: dim
    real :: angle_2d, in_angle, cur_angle

    dim = size(vec_A, 1)

    if (dim == 2) then
      angle_2d = get_angle_2d(vec_A(:, perm_B(idx)), vec_B(:, perm_A(idx)))
      if (angle_2d < 0.0) then
        angle_2d = -1 * angle
      else
        angle_2d = angle
      end if
      mat(1, 1) = cos(angle_2d) ; mat(1, 2) = -sin(angle_2d)
      mat(2, 1) = sin(angle_2d) ; mat(2, 2) = cos(angle_2d)
      vec_A = matmul(mat, vec_A)
    else if (dim == 3) then
      in_angle = get_angle(vec_A(:, perm_B(idx)), vec_B(:, perm_A(idx)))
      cross = cross_product(vec_A(:, perm_B(idx)), vec_B(:, perm_A(idx))) ; cross = cross / norm(cross)
      mat = get_rotation_matrix_cross(cross, angle)
      tmpvec = matmul(mat, vec_A(:, perm_B(idx)))
      cur_angle = get_angle(tmpvec, vec_B(:, perm_A(idx)))
      if (cur_angle > in_angle) then ! got it the wrong way round
        mat = get_rotation_matrix_cross(-1 * cross, angle)
      end if
      vec_A = matmul(mat, vec_A)
    end if
  end subroutine rotate_vec

  subroutine match_up_vectors(vec_P, permutation_P, vec_Q, permutation_Q)
    !!< We match up vectors by matching up each angle in turn.
    real, dimension(:, :), intent(in) :: vec_P, vec_Q
    integer, dimension(:), intent(out)   :: permutation_P, permutation_Q
 
    integer, dimension(2) :: idx
    real, dimension(size(vec_P, 2), (size(vec_Q, 2))) :: angle_matrix
 
    integer :: i, j, k, count_P, count_Q, dim, stat
 
    dim = size(vec_P, 1)
    count_P = size(vec_P, 2)
    count_Q = size(vec_Q, 2)

    permutation_P = 0
    permutation_Q = 0
 
    ! construct a record of the angles between each of the eigenvectors.

    do i=1,count_P
      do j=1,count_Q
        angle_matrix(i, j) = get_angle(vec_P(:, i), vec_Q(:, j))
      end do
    end do

    ! now loop over and find the closest pair.

    do k=1,count_P
      idx = minloc(angle_matrix); i = idx(1); j = idx(2)
      angle_matrix(i, :) = huge(0.0); angle_matrix(:, j) = huge(0.0)
      permutation_P(k) = j; permutation_Q(k) = i
    end do

    stat = 0
    if (size(permutation_P) == dim .and. size(permutation_Q) == dim) then
      call check_perm(permutation_P, stat)
      call check_perm(permutation_Q, stat)
    end if

    if (stat /= 0) then
      call write_matrix(vec_P, "vec_P")
      call write_matrix(vec_Q, "vec_Q")
      call write_vector(permutation_P, "perm_P")
      call write_vector(permutation_Q, "perm_Q")
      FLAbort("Permutation not correct!")
    end if
  end subroutine match_up_vectors

  subroutine match_up_ellipsoids(vec_P, val_P, perm_P, vec_Q, val_Q, perm_Q)
    !!< Match up the index of the biggest eigenvalue with the biggest eigenvalue,
    !!< the next biggest with the next biggest, etc.
    real, dimension(:, :), intent(in) :: vec_P, vec_Q
    real, dimension(:), intent(in) :: val_P, val_Q
    integer, dimension(size(val_P)), intent(out) :: perm_P, perm_Q
    real, dimension(size(val_P)) :: lval_P, lval_Q
    real, dimension(size(val_P), size(val_P) - 1) :: equatorial_basis_P, equatorial_basis_Q
    integer, dimension(size(val_P) - 1) :: eq_perm_P, eq_perm_Q
    integer :: i, dim, biggest_loc(1), j
    integer :: stat

    dim = size(val_P)

    if (metric_spheroid(val_P) .and. metric_spheroid(val_Q)) then
      perm_Q(1) = get_polar_index(val_P)
      perm_P(1) = get_polar_index(val_Q)

      j = 1
      do i=1,dim
        if (i /= perm_Q(1)) then
          equatorial_basis_P(:, j) = vec_P(:, i)
          eq_perm_P(j) = i
          j = j + 1
        end if
      end do

      j = 1
      do i=1,dim
        if (i /= perm_P(1)) then
          equatorial_basis_Q(:, j) = vec_Q(:, i)
          eq_perm_Q(j) = i
          j = j + 1
        end if
      end do

      call match_up_vectors(equatorial_basis_P, perm_P(2:), equatorial_basis_Q, perm_Q(2:))

      do i=1,dim-1
        perm_P(i+1) = eq_perm_Q(perm_P(i+1))
        perm_Q(i+1) = eq_perm_P(perm_Q(i+1))
      end do

      return
    end if

    lval_P = val_P; lval_Q = val_Q

    do i=1,dim
      biggest_loc = maxloc(lval_P)
      perm_Q(i) = biggest_loc(1)
      lval_P(biggest_loc(1)) = 0.0

      biggest_loc = maxloc(lval_Q)
      perm_P(i) = biggest_loc(1)
      lval_Q(biggest_loc(1)) = 0.0
    end do

    stat = 0
    call check_perm(perm_P, stat)
    call check_perm(perm_Q, stat)
    if (stat /= 0) then
      call write_vector(val_P, "val_P")
      call write_vector(val_Q, "val_Q")
      call write_matrix(vec_P, "vec_P")
      call write_matrix(vec_Q, "vec_Q")
      call write_vector(perm_P, "perm_P")
      call write_vector(perm_Q, "perm_Q")
      FLAbort("Permutation not correct!")
    end if
  end subroutine match_up_ellipsoids

  subroutine reduce_edgelen(mat_P, vec_P, val_P, changed_P, mat_Q, vec_Q, val_Q, changed_Q, positions, p, q)
    !!< Reduce the edge lengths to smoothen the mesh, according to section
    !!< 3.3 of the referenced paper (with modifications).
    real, dimension(:, :), intent(in) :: mat_P, mat_Q
    real, dimension(size(mat_P, 1), size(mat_P, 1)), intent(in) :: vec_P, vec_Q
    real, dimension(size(mat_P, 1)), intent(inout) :: val_P, val_Q
    logical, intent(inout) :: changed_P, changed_Q
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: p, q

    real, dimension(size(mat_Q, 1)) :: edgelen_P, edgelen_Q
    integer, dimension(size(mat_Q, 1)) :: perm_P, perm_Q
    integer :: i, dim
    real :: gamma, dist, tmp

    dim = size(mat_P, 1)

    do i=1,dim
      edgelen_P(i) = edge_length_from_eigenvalue(val_P(i))
      edgelen_Q(i) = edge_length_from_eigenvalue(val_Q(i))
    end do

    call match_up_vectors(vec_P, perm_P, vec_Q, perm_Q)
    dist = distance(positions, p, q)

#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
    write(0,*) "perm_P == ", perm_P
    write(0,*) "perm_Q == ", perm_Q
    call write_vector(edgelen_P, "edgelen_P")
    call write_vector(edgelen_Q, "edgelen_Q")
#endif

    do i=1,dim
      gamma = get_gamma(edgelen_P(perm_Q(i)), vec_P(:, perm_Q(i)), edgelen_Q(perm_P(i)), vec_Q(:, perm_P(i)), dist)
#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
      write(0,*) "gamma == ", gamma, "; gamma0 == ", gamma0
      call write_vector(edgelen_P, "edgelen_P")
      call write_vector(edgelen_Q, "edgelen_Q")
#endif
      if (gamma .fgt. gamma0) then
        if (edgelen_P(perm_Q(i)) > edgelen_Q(perm_P(i))) then
#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
          write(0,*) "reducing edge length ", perm_Q(i), " of P." 
          write(0,*) "old value == ", edgelen_P(perm_Q(i))
#endif
          tmp = edgelen_P(perm_Q(i))
          edgelen_P(perm_Q(i)) = dist * log(gamma0) + edgelen_Q(perm_P(i))
          if (tmp .fne. edgelen_P(perm_Q(i))) changed_P = .true.
#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
          write(0,*) "new value == ", edgelen_P(perm_Q(i))
#endif
        else
#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
          write(0,*) "reducing edge length ", perm_P(i), " of Q." 
          write(0,*) "old value == ", edgelen_Q(perm_P(i))
#endif
          tmp = edgelen_Q(perm_P(i))
          edgelen_Q(perm_P(i)) = distance(positions, p, q) * log(gamma0) + edgelen_P(perm_Q(i))
          if (tmp .fne. edgelen_Q(perm_P(i))) changed_Q = .true.
#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
          write(0,*) "new value == ", edgelen_Q(perm_P(i))
#endif
        end if
      end if
    end do

    do i=1,dim
      val_P(i) = eigenvalue_from_edge_length(edgelen_P(i))
      val_Q(i) = eigenvalue_from_edge_length(edgelen_Q(i))
    end do

#ifdef EXTRA_SPECIAL_GRADATION_DEBUGGING
    call write_vector(val_P, "P output eigenvalues")
    call write_vector(val_Q, "Q output eigenvalues")
#endif
  end subroutine reduce_edgelen

  function get_gamma(h_P, v_P, h_Q, v_Q, dist) result(gamma)
    !!< Gamma is the mesh size gradation measure,
    !!< a measure of how quickly the desired edge length changes
    !!< when going from node p to node q.
    !!< See section 2 of the referenced paper.
    real, intent(in) :: h_P, h_Q, dist
    real, dimension(:), intent(in) :: v_P, v_Q
    real :: gamma

    gamma = exp(abs(h_P - h_Q) / dist)
    !gamma = gamma + (1 - gamma) * (1 - dot_product(v_P, v_Q))
  end function get_gamma
end module gradation_metric
