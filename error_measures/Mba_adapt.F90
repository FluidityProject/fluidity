#include "fdebug.h"

module mba_adapt_module
  use fields
  use state_module
  use interpolation_module
  use meshdiagnostics
  use vtk_interfaces
  use eventcounter
  use node_boundary
  use metric_tools
  use limit_metric_module
#ifdef HAVE_MBA_2D
  use mba2d_module
#endif
  implicit none

  contains

  subroutine mba_adapt(state, metric)
    type(state_type), intent(inout) :: state
    type(tensor_field), intent(in) :: metric

#ifdef HAVE_MBA_2D

    type(mesh_type), pointer :: xmesh
    type(vector_field), pointer :: positions

    integer :: nonods, mxnods, stotel, mxface, totele, maxele
    real, dimension(:, :), allocatable :: pos
    integer, dimension(:, :), allocatable :: ipf
    integer, dimension(:, :), allocatable :: ipe
    real, dimension(:, :), allocatable :: parcrv
    integer, dimension(:), allocatable :: ipv
    integer, dimension(:), allocatable :: iFnc
    integer, dimension(:), allocatable :: lbE
    real, dimension(:, :), allocatable :: tmp_metric
    integer :: i, j
    real :: rQuality
    integer :: ierr, maxWr, maxWi
    real, dimension(:), allocatable :: rW
    integer, dimension(:), allocatable :: iW
    integer :: status
    type(state_type) :: new_state
    type(mesh_type) :: new_mesh
    type(scalar_field), dimension(scalar_field_count(state)) :: new_scalar
    type(vector_field), dimension(vector_field_count(state)-1) :: new_vector
    type(tensor_field), dimension(tensor_field_count(state)) :: new_tensor
    type(vector_field) :: new_positions
    type(scalar_field), pointer :: old_scalar
    type(vector_field), pointer :: old_vector
    type(tensor_field), pointer :: old_tensor
    integer :: npv
    integer :: xpctel

    assert(metric%dim == 2)

    positions => extract_vector_field(state, "Coordinate")
    xmesh => extract_mesh(state, "Mesh")
    call initialise_boundcount(xmesh, positions)

    if (.not. has_faces(xmesh)) then
      call add_faces(xmesh)
    end if

    ! Patrick, is this right? -Stephan
    mxnods = 50000

    nonods = node_count(xmesh)
    totele = ele_count(xmesh)
    stotel = surface_element_count(xmesh)
    mxface = int(max((float(mxnods) / float(nonods)) * stotel * 1.2, 10000.0))
    maxele = int(max((float(mxnods) / float(nonods)) * totele * 1.2, 10000.0))
    xpctel = max(expected_elements(positions, metric), 5)

    allocate(pos(2, mxnods))
    pos = 0.0
    do i=1,2
      pos(i, 1:nonods) = positions%val(i)%ptr
    end do

    allocate(ipf(4, mxface))
    ipf = 0
    do i=1,stotel
      ipf(1:2, i) = face_global_nodes(xmesh, i)
      ipf(3, i) = 0
    end do
    ipf(4, 1:stotel) = xmesh%faces%boundary_ids

    allocate(ipe(3, maxele))
    ipe = 0
    do i=1,totele
      ipe(:, i) = ele_nodes(xmesh, i)
    end do

    npv = count(boundcount > 1)
    allocate(ipv(npv))
    j = 1
    do i=1,nonods
      if (boundcount(i) > 1) then
        ipv(j) = i
        j = j + 1
      end if
    end do

    allocate(parcrv(2, mxface))
    parcrv = 0

    allocate(iFnc(mxface))
    iFnc = 0

    allocate(lbE(maxele))
    lbE = 1

    allocate(tmp_metric(3, mxnods))
    tmp_metric = 0.0
    do i=1,nonods
      tmp_metric(1, i) = metric%val(1, 1, i)
      tmp_metric(2, i) = metric%val(2, 2, i)
      tmp_metric(3, i) = metric%val(1, 2, i)
    end do

    maxWr = (4 * mxnods + 10 * nonods + mxface + maxele)  * 1.2
    maxWi = (6 * mxnods + 10 * nonods + 19 * mxface + 11 * maxele + 12 * totele) * 1.2
    allocate(rW(maxWr))
    allocate(iW(maxWi))

    status = 1

    CALL mbaNodal(                                   &
         nonods, mxnods, stotel, mxface, totele, maxele, npv, &
         pos, ipf, ipe, ipv, &
         CrvFunction_ani, parcrv, iFnc, &
         xpctel, &
         0, 0, ipv, ipv, lbE, &
         .true., status, &
         100, 30000, &
         tmp_metric, 0.9, rQuality, &
         maxWr, maxWi, rW, iW, &
         10, ierr)

    call incrementeventcounter(EVENT_ADAPTIVITY)
    call incrementeventcounter(EVENT_MESH_MOVEMENT)

    ! Hooray! You didn't crash. Congratulations. Now let's assemble the output and interpolate.

    call allocate(new_mesh, nonods, totele, ele_shape(xmesh, 1), "Mesh")
    new_mesh%ndglno = reshape(IPE(:, 1:totele), (/size(new_mesh%ndglno)/))
    call insert(new_state, new_mesh, "Mesh")
    call deallocate(new_mesh)

    call allocate(new_positions, 2, new_mesh, "Coordinate")
    call set_all(new_positions, pos(:, 1:nonods))
    call snap_positions(positions, new_positions)
    call insert(new_state, new_positions, "Coordinate")
    call deallocate(new_positions)

    do i=1,scalar_field_count(state)
      old_scalar => extract_scalar_field(state, i)
      call allocate(new_scalar(i), new_mesh, trim(old_scalar%name))
      call zero(new_scalar(i))
      call insert(new_state, new_scalar(i), trim(old_scalar%name))
      call deallocate(new_scalar(i))
    end do

    j = 1
    do i=1,vector_field_count(state)
      old_vector => extract_vector_field(state, i)
      if (trim(old_vector%name) == "Coordinate") then
        cycle
      end if
      call allocate(new_vector(j), 2, new_mesh, trim(old_vector%name))
      call zero(new_vector(j))
      call insert(new_state, new_vector(j), trim(old_vector%name))
      call deallocate(new_vector(j))
      j = j + 1
    end do

    do i=1,tensor_field_count(state)
      old_tensor => extract_tensor_field(state, i)
      call allocate(new_tensor(i), new_mesh, trim(old_tensor%name))
      call zero(new_tensor(i))
      call insert(new_state, new_tensor(i), trim(old_tensor%name))
      call deallocate(new_tensor(i))
    end do

    call vtk_write_state("new_state", 0, state=(/new_state/))
    call linear_interpolation(state, new_state)
    
    deallocate(pos)
    deallocate(ipf)
    deallocate(ipe)
    deallocate(ipv)
    deallocate(parcrv)
    deallocate(iFnc)
    deallocate(lbE)
    deallocate(rW)
    deallocate(iW)
    deallocate(tmp_metric)

    call deallocate(state)
    state = new_state
#else
    FLAbort("You called mba_adapt without the mba library.")
#endif
  end subroutine mba_adapt

  subroutine CrvFunction_ani(tc, xyc, iFnc)
  real  tc, xyc(2)
  integer :: iFnc

  return
  end subroutine

  ! Try to snap the positions of new nodes to the positions
  ! of old nodes. libmba, for example, maps domains to the
  ! square [0.1, 0.9]^2. This mapping back and forth causes small
  ! errors in the positions that screws up various things.
  ! Let's put them back ..
  ! The inputs to this are the linear topological meshes
  ! created by the adaptivity library
  subroutine snap_positions(old_positions, new_positions)
    type(vector_field), intent(inout) :: old_positions
    type(vector_field), intent(inout) :: new_positions

    real, dimension(node_count(new_positions)) :: map
    integer :: node_new, ele_old, node_old, loc, i, j, dim
    integer, dimension(:), pointer :: ele_nodes_old, ele_faces_old
    real :: snap_dist, ele_dist
    logical, dimension(node_count(old_positions)) :: used_old
    integer :: face_old
    real, dimension(old_positions%dim, old_positions%dim-1) :: basis
    real, dimension(old_positions%dim) :: proj

    used_old = .false.
    loc = ele_loc(old_positions, 1)
    map = get_element_mapping(old_positions, new_positions)
    dim = new_positions%dim
    if (.not. has_faces(old_positions%mesh)) then
      call add_faces(old_positions%mesh)
    end if

    nodeloop: do node_new=1,node_count(new_positions)
      ele_old = map(node_new)
      ele_nodes_old => ele_nodes(old_positions, ele_old)
      ele_faces_old => ele_faces(old_positions, ele_old)
      ele_dist = huge(0.0)
      do i=1,loc
        do j=i+1,loc
          ele_dist = min(ele_dist, norm2(node_val(old_positions, ele_nodes_old(i)) - node_val(old_positions, ele_nodes_old(j))))
        end do
      end do

      snap_dist = ele_dist / 1.0e3

      ! Snap to nodes
      do i=1,loc
        node_old = ele_nodes_old(i)
        if (norm2(node_val(new_positions, node_new) - node_val(old_positions, node_old)) < snap_dist .and. .not. used_old(node_old)) then
          used_old(node_old) = .true.
          do j=1,dim
            new_positions%val(j)%ptr(node_new) = old_positions%val(j)%ptr(node_old)
          end do
          cycle nodeloop
        end if

      end do

      ! Snap to faces
     do i=1,loc
       face_old = ele_faces_old(i)
       basis = face_basis(old_positions, face_old)
       proj = project_to_subspace(node_val(new_positions, node_new), basis)
       if (norm2(node_val(new_positions, node_new) - proj) < snap_dist) then
         do j=1,dim
           new_positions%val(j)%ptr(node_new) = proj(j)
         end do
         cycle nodeloop
       end if
     end do
    end do nodeloop
  end subroutine snap_positions

  ! Get the basis of the subspace containing the face of an element -- only works
  ! for linear positions.
  function face_basis(positions, face) result(basis)
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: face
    real, dimension(positions%dim, positions%dim-1) :: basis
    
    real, dimension(positions%dim, face_loc(positions, face)) :: pos
    integer :: i, j
    real, dimension(positions%dim) :: tmp_vec
    integer :: dim

    dim = positions%dim

    pos = face_val(positions, face)

    ! Get edge vectors and normalise
    do i=1,dim-1
      basis(:, i) = pos(:, i) - pos(:, i+1)
      basis(:, i) = basis(:, i) / norm2(basis(:, i))
    end do

    ! Orthogonalise
    if (dim > 2) then
      do i=2,dim
        tmp_vec = basis(:, i)
        do j=1,i-1
          tmp_vec = tmp_vec - (dot_product(basis(:, j), tmp_vec) * basis(:, j))
        end do
        basis(:, i) = tmp_vec / norm2(tmp_vec)
      end do
    end if

    call check_basis(basis)
  end function face_basis

end module mba_adapt_module
