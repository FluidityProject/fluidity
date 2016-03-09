#include "fdebug.h"

module goals
!!< A repository of goals and their gradients,
!!< suitable for use with goal-based optimisation.

  use fldebug
  use vector_tools
  use elements
  use spud
  use tensors, only: tensormul
  use fetools
  use unittest_tools, only: get_matrix_identity
  use fields
  use state_module
  use fefields
  use field_derivatives, only: grad, curl

  implicit none

  public :: goal_temp, goal_temp_grad, &
          & goal_enstrophy, goal_enstrophy_grad, &
          & goal_les_velocity, goal_les_velocity_grad, &
          & goal_les_velocity_squared, goal_les_velocity_squared_grad, &
          & goal_les_velocity_4th, goal_les_velocity_4th_grad, &
          & compute_goals
  private

  contains

  subroutine compute_goals(state)
    ! Check to see what goals are requested
    ! and print them out per timestep
    type(state_type), dimension(:), intent(in) :: state
    real :: acctim
    real :: goal_val

    call get_option("/timestepping/current_time", acctim)

    if (have_option("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity/enstrophy_goal")) then
      goal_val = goal_enstrophy(state)
      ewrite(1,*) "enstrophy_goal at time ", acctim, ": ", goal_val
    end if

    if (have_option("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity/temperature_gradient_goal")) then
      goal_val = goal_temp(state)
      ewrite(1,*) "temperature_gradient_goal at time ", acctim, ": ", goal_val
    end if

    if (have_option("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity/les_goal")) then
      goal_val = goal_les_velocity(state)
      ewrite(1,*) "les_goal at time ", acctim, ": ", goal_val
    end if

    if (have_option("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity/higher_order_les_goal")) then
      goal_val = goal_les_velocity_4th(state)
      ewrite(1,*) "higher_order_les_goal at time ", acctim, ": ", goal_val
    end if

  end subroutine compute_goals

  function goal_temp(state)
    ! Basic sample goal:
    ! int( transpose(grad(T)) . kappa . grad(T) ) dV.
    ! kappa is a matrix that weights gradients
    ! in certain directions more than others.
    type(state_type), dimension(:), intent(in) :: state
    real :: goal_temp

    integer :: ele, ngi
    real, dimension(:), allocatable :: detwei, result
    real, dimension(:, :, :), allocatable :: dn_t
    type(scalar_field), pointer :: temp
    type(vector_field), pointer :: positions
    type(element_type), pointer :: t_shape
    real, dimension(:, :), allocatable :: kappa
    real, dimension(:), allocatable :: vec

    temp => extract_scalar_field(state(1), "Temperature")
    positions => extract_vector_field(state(1), "Coordinate")

    allocate(detwei(ele_ngi(temp, 1)), result(ele_ngi(temp, 1)))
    allocate(dn_t(ele_loc(temp, 1), ele_ngi(temp, 1), positions%dim))
    allocate(kappa(positions%dim, positions%dim))
    allocate(vec(positions%dim))
    t_shape => ele_shape(temp, 1)

    goal_temp = 0.0

    kappa = get_matrix_identity(size(kappa, 1))

    do ele=1,ele_count(temp)
      call transform_to_physical(positions, ele, t_shape, dshape=dn_t, detwei=detwei)

      do ngi=1,ele_ngi(temp, 1)
        vec = matmul(ele_val(temp, ele), dn_t(:, ngi, :)) ! grad(T)
        result(ngi) = dot_product(vec, matmul(kappa, vec))
      end do

      goal_temp = goal_temp + dot_product(result, detwei)
    end do

    deallocate(detwei)
    deallocate(dn_t)
    deallocate(kappa)
    deallocate(vec)
  end function goal_temp

  subroutine goal_temp_grad(state, dep, adj)
    ! The derivative of the goal with respect to
    ! whatever variables it depends on.
    ! Here it only depends on "Temperature", so 
    ! other variables are ignored.
    type(state_type), dimension(:), intent(in) :: state
    character(len=*), intent(in) :: dep
    type(scalar_field), intent(inout) :: adj

    type(vector_field), pointer :: positions
    type(element_type), pointer :: t_shape
    type(scalar_field), pointer :: temp
    integer :: ele, loc, ngi
    integer, dimension(:), pointer :: nodelist
    real, dimension(ele_ngi(adj, 1)) :: detwei
    real, dimension(ele_loc(adj, 1), ele_ngi(adj, 1), mesh_dim(adj)) :: dn_t
    real, dimension(ele_loc(adj, 1), ele_ngi(adj, 1)) :: tmpA
    real, dimension(mesh_dim(adj), mesh_dim(adj)) :: kappa

    call zero(adj)

    if (dep(1:11) == "Temperature") then
      temp => extract_scalar_field(state(1), "Temperature")
    else if (dep(1:14) == "OldTemperature") then
      temp => extract_scalar_field(state(1), "OldTemperature")
    else
      return
    end if

    positions => extract_vector_field(state(1), "Coordinate")
    t_shape => ele_shape(adj, 1)

    kappa = get_matrix_identity(size(kappa, 1))

    do ele=1,ele_count(adj)
      nodelist => ele_nodes(adj, ele)
      call transform_to_physical(positions, ele, t_shape, dshape=dn_t, detwei=detwei)
      do loc=1,ele_loc(adj, ele)
        do ngi=1,ele_ngi(adj,ele)
          tmpA(loc, ngi) = dot_product(dn_t(loc, ngi, :), matmul(kappa, dn_t(loc, ngi, :)))
        end do
      end do

      ! tmpA now represents grad(N_i) . kappa . grad(N_i)
      
      do loc=1,ele_loc(adj, ele)
        tmpA(loc, :) = tmpA(loc, :) * node_val(temp, nodelist(loc)) 
      end do

      ! tmpA now represents T_i . (grad(N_i) . kappa . grad(N_i))

      call addto(adj, ele_nodes(adj, ele), 2.0 * matmul(tmpA, detwei))
    end do

  end subroutine goal_temp_grad

  function goal_enstrophy(state)
    ! Compute the enstrophy of the system:
    ! 1/2 * int(|curl(velocity)|**2) dV.
    type(state_type), dimension(:), intent(in) :: state
    real :: goal_enstrophy

    integer :: ele, node
    real, dimension(:), allocatable :: detwei
    type(scalar_field) :: vorticity
    type(vector_field), pointer :: positions, velocity
    type(element_type), pointer :: t_shape

    positions => extract_vector_field(state(1), "Coordinate")
    velocity  => extract_vector_field(state(1), "Velocity")
    call allocate(vorticity, velocity%mesh, "|Vorticity|**2")

    call curl(velocity, positions, curl_norm=vorticity)
    do node=1,node_count(vorticity)
      vorticity%val(node) = vorticity%val(node)**2
    end do

    allocate(detwei(ele_ngi(vorticity, 1)))
    t_shape => ele_shape(vorticity, 1)

    goal_enstrophy = 0.0

    do ele=1,ele_count(vorticity)
      call transform_to_physical(positions, ele, detwei=detwei)
      goal_enstrophy = goal_enstrophy + dot_product(ele_val_at_quad(vorticity, ele), detwei)
    end do

    goal_enstrophy = goal_enstrophy * 0.5

    deallocate(detwei)
    call deallocate(vorticity)
  end function goal_enstrophy

  subroutine goal_enstrophy_grad(state, dep, adj)
    ! The derivative of the goal with respect to
    ! whatever variables it depends on.
    ! Here it depends on "Velocity[123]".
    type(state_type), dimension(:), intent(in) :: state
    character(len=*), intent(in) :: dep
    type(scalar_field), intent(inout) :: adj

    type(vector_field), pointer :: positions, velocity
    type(element_type), pointer :: t_shape
    integer :: ele, node, j, i
    real, dimension(ele_ngi(adj, 1)) :: detwei
    real, dimension(ele_loc(adj, 1), ele_ngi(adj, 1), mesh_dim(adj)) :: dn_t
    type(mesh_type) :: mesh
    ! curdim: which component of velocity are we processing
    ! odima, odimb: the other components of velocity
    ! curloc: the index of the current node in the list of nodes
    ! associated with a particular element, i.e. where does it appear
    ! in the local list of nodes
    integer :: curdim, curloc, odima, odimb
    integer, dimension(:), pointer :: e_n ! ele_nodes(mesh, ele)
    type(patch_type) :: patch
    real, dimension(ele_ngi(adj, 1)) :: r

    call zero(adj)
    if (dep(1:8) /= "Velocity") then
      return
    end if

    mesh = adj%mesh
    call add_nelist(mesh)

    read(dep(10:10), *) curdim ! Velocity1 -> 1 etc.
    if (curdim == 1) then
      odima = 2
      odimb = 3
    else if (curdim == 2) then
      odima = 1
      odimb = 3
    else
      odima = 1
      odimb = 2
    end if

    positions => extract_vector_field(state(1), "Coordinate")
    velocity  => extract_vector_field(state(1), "Velocity")
    t_shape => ele_shape(adj, 1)

    do node=1,node_count(adj)
      patch = get_patch_ele(mesh, node, level=1)
      do j=1,patch%count
        ele = patch%elements(j)
        e_n => ele_nodes(mesh, ele)
        do i=1,ele_loc(mesh, ele)
          if (node == e_n(i)) then
            curloc = i
            exit
          end if
        end do

        call transform_to_physical(positions, ele, t_shape, dshape=dn_t, detwei=detwei)
        r = 2 * node_val(velocity, curdim, node) * dn_t(curloc, :, odima)**2
        r = r + 2 * node_val(velocity, curdim, node) * dn_t(curloc, :, odimb)**2
        r = r - node_val(velocity, odima, node) * (dn_t(curloc, :, curdim) * dn_t(curloc, :, odima))
        r = r - node_val(velocity, odimb, node) * (dn_t(curloc, :, curdim) * dn_t(curloc, :, odimb))
        r = 0.5 * r
        call addto(adj, node, dot_product(r, detwei))
      end do
    end do

  end subroutine goal_enstrophy_grad

  !! LES goal-based error measures
  function goal_les_temp(state)
    ! int( transpose(grad(T)) . kappa . grad(T) ) dV.
    ! kappa is the LES matrix.
    type(state_type), dimension(:), intent(in) :: state
    real :: goal_les_temp

    integer :: ele, ngi
    real, dimension(:), allocatable :: detwei, result
    real, dimension(:, :, :), allocatable :: dn_t
    type(scalar_field), pointer :: temp
    type(vector_field), pointer :: positions
    type(element_type), pointer :: t_shape
    real, dimension(:, :), allocatable :: kappa
    real, dimension(:), allocatable :: vec
    type(vector_field), pointer :: nvelocity, gvelocity
    type(tensor_field), pointer :: mesh_metric
    real, dimension(:, :, :), allocatable :: mesh_size_tensor_ele
    real, dimension(:, :), allocatable :: mesh_size_tensor
    integer :: stat, i

    temp => extract_scalar_field(state(1), "Temperature")
    positions => extract_vector_field(state(1), "Coordinate")
    nvelocity => extract_vector_field(state(1), "NonlinearVelocity")
    gvelocity => extract_vector_field(state(1), "GridVelocity")
    mesh_metric => extract_tensor_field(state(1), "MeshSizingMetric", stat=stat)

    allocate(detwei(ele_ngi(temp, 1)), result(ele_ngi(temp, 1)))
    allocate(dn_t(ele_loc(temp, 1), ele_ngi(temp, 1), positions%dim))
    allocate(kappa(positions%dim, positions%dim))
    allocate(vec(positions%dim))
    t_shape => ele_shape(temp, 1)

    goal_les_temp = 0.0

    if (stat /= 0) then ! MeshSizingMetric not found
      do ele=1,ele_count(temp)
        call transform_to_physical(positions, ele, t_shape, dshape=dn_t, detwei=detwei)

        do ngi=1,ele_ngi(temp, 1)
          kappa = les_tensor(nvelocity, gvelocity, dn_t, ele, ngi)
          vec = matmul(ele_val(temp, ele), dn_t(:, ngi, :)) ! grad(T)
          result(ngi) = dot_product(vec, matmul(kappa, vec))
        end do

        goal_les_temp = goal_les_temp + dot_product(result, detwei)
      end do
    else
      allocate(mesh_size_tensor_ele(positions%dim, positions%dim, ele_loc(temp, 1)))
      allocate(mesh_size_tensor(positions%dim, positions%dim))
      do ele=1,ele_count(temp)
        call transform_to_physical(positions, ele, t_shape, dshape=dn_t, detwei=detwei)
        mesh_size_tensor_ele = ele_val(mesh_metric, ele)
        mesh_size_tensor = 0.0
        do i=1,ele_loc(temp, 1)
          mesh_size_tensor = mesh_size_tensor + mesh_size_tensor_ele(:, :, i)
        end do
        mesh_size_tensor = mesh_size_tensor / ele_loc(temp, 1)
        call invert(mesh_size_tensor)

        do ngi=1,ele_ngi(temp, 1)
          kappa = les_tensor(nvelocity, gvelocity, dn_t, ele, ngi, size_tensor=mesh_size_tensor)
          vec = matmul(ele_val(temp, ele), dn_t(:, ngi, :)) ! grad(T)
          result(ngi) = dot_product(vec, matmul(kappa, vec))
        end do

        goal_les_temp = goal_les_temp + dot_product(result, detwei)
      end do
      deallocate(mesh_size_tensor)
      deallocate(mesh_size_tensor_ele)
    end if

    goal_les_temp = goal_les_temp * 0.5

    deallocate(detwei)
    deallocate(dn_t)
    deallocate(kappa)
    deallocate(vec)
  end function goal_les_temp

  function goal_les_velocity(state)
    ! 0.5 * int( transpose(grad(T)) . kappa . grad(T) ) dV.
    ! kappa is the LES matrix.
    type(state_type), dimension(:), intent(in) :: state
    real :: goal_les_velocity

    integer :: ele, ngi
    real, dimension(:), allocatable :: detwei, result
    real, dimension(:, :, :), allocatable :: dn_t
    type(vector_field), pointer :: positions
    type(element_type), pointer :: t_shape
    real, dimension(:, :), allocatable :: kappa
    real, dimension(:), allocatable :: vec
    type(vector_field), pointer :: nvelocity, gvelocity
    integer :: i
    type(tensor_field), pointer :: mesh_metric
    real, dimension(:, :, :), allocatable :: mesh_size_tensor_ele
    real, dimension(:, :), allocatable :: mesh_size_tensor
    integer :: stat

    positions => extract_vector_field(state(1), "Coordinate")
    nvelocity => extract_vector_field(state(1), "NonlinearVelocity")
    gvelocity => extract_vector_field(state(1), "GridVelocity")
    mesh_metric => extract_tensor_field(state(1), "MeshSizingMetric", stat=stat)

    allocate(detwei(ele_ngi(nvelocity, 1)), result(ele_ngi(nvelocity, 1)))
    allocate(dn_t(ele_loc(nvelocity, 1), ele_ngi(nvelocity, 1), positions%dim))
    allocate(kappa(positions%dim, positions%dim))
    allocate(vec(positions%dim))
    t_shape => ele_shape(nvelocity, 1)

    goal_les_velocity = 0.0

    if (stat /= 0) then ! MeshSizingMetric not found
      do ele=1,ele_count(nvelocity)
        call transform_to_physical(positions, ele, t_shape, dshape=dn_t, detwei=detwei)


        do ngi=1,ele_ngi(nvelocity, 1)
          kappa = les_tensor(nvelocity, gvelocity, dn_t, ele, ngi)
          result(ngi) = 0.0
          do i=1,nvelocity%dim
            vec = matmul(ele_val(nvelocity, i, ele), dn_t(:, ngi, :))
            result(ngi) = result(ngi) + dot_product(vec, matmul(kappa, vec))
          end do
        end do

        goal_les_velocity = goal_les_velocity + dot_product(result, detwei)
      end do
    else
      allocate(mesh_size_tensor_ele(positions%dim, positions%dim, ele_loc(positions, 1)))
      allocate(mesh_size_tensor(positions%dim, positions%dim))
      do ele=1,ele_count(nvelocity)
        call transform_to_physical(positions, ele, t_shape, dshape=dn_t, detwei=detwei)

        mesh_size_tensor_ele = ele_val(mesh_metric, ele)
        mesh_size_tensor = 0.0
        do i=1,ele_loc(positions, 1)
          mesh_size_tensor = mesh_size_tensor + mesh_size_tensor_ele(:, :, i)
        end do
        mesh_size_tensor = mesh_size_tensor / ele_loc(positions, 1)
        call invert(mesh_size_tensor)

        do ngi=1,ele_ngi(nvelocity, 1)
          kappa = les_tensor(nvelocity, gvelocity, dn_t, ele, ngi, size_tensor=mesh_size_tensor)
          result(ngi) = 0.0
          do i=1,nvelocity%dim
            vec = matmul(ele_val(nvelocity, i, ele), dn_t(:, ngi, :))
            result(ngi) = result(ngi) + dot_product(vec, matmul(kappa, vec))
          end do
        end do

        goal_les_velocity = goal_les_velocity + dot_product(result, detwei)
      end do
      deallocate(mesh_size_tensor)
      deallocate(mesh_size_tensor_ele)
    end if

    goal_les_velocity = goal_les_velocity * 0.5

    deallocate(detwei)
    deallocate(result)
    deallocate(dn_t)
    deallocate(kappa)
    deallocate(vec)
  end function goal_les_velocity

  function les_tensor(nvelocity, gvelocity, dn_t, ele, gi, size_tensor) result(kappa)
    !! Get the LES tensor at gauss point gi for element ele.
    type(vector_field), intent(in) :: gvelocity, nvelocity
    real, dimension(:, :, :), intent(in) :: dn_t
    real, dimension(3, 3) :: kappa
    real, dimension(3, 3), optional, intent(in) :: size_tensor
    integer, intent(in) :: ele, gi

    real, dimension(ele_ngi(nvelocity, ele)) :: mxx, mxy, mxz, myy, myz, mzz ! ugh

    if (present(size_tensor)) then
      mxx(gi) = size_tensor(1, 1)
      mxy(gi) = size_tensor(1, 2)
      mxz(gi) = size_tensor(1, 3)
      myy(gi) = size_tensor(2, 2)
      myz(gi) = size_tensor(2, 3)
      mzz(gi) = size_tensor(3, 3)
    else
      call SIZEGIELETENS(dn_t(:, :, 1), dn_t(:, :, 2), dn_t(:, :, 3), ele_loc(nvelocity, ele), &
                    &  ele_ngi(nvelocity, ele), gi, MXX(gi), MXY(gi), MXZ(gi), MYY(gi), MYZ(gi), MZZ(gi))
    end if

    call LESVIS(node_count(nvelocity), ele_count(nvelocity), ele_loc(nvelocity, ele), &
              & nvelocity%mesh%ndglno, ele_ngi(nvelocity, ele), ele, gi, &
              & dn_t(:, :, 1), dn_t(:, :, 2), dn_t(:, :, 3), &
              & nvelocity%val(1,:), nvelocity%val(2,:), nvelocity%val(3,:), & 
              & gvelocity%val(1,:), gvelocity%val(2,:), gvelocity%val(3,:), &
              & MXX, MXY, MXZ, MYY, MYZ, MZZ) 

    kappa(1, 1) = mxx(gi)
    kappa(1, 2) = mxy(gi)
    kappa(2, 1) = mxy(gi)
    kappa(1, 3) = mxz(gi)
    kappa(3, 1) = mxz(gi)
    kappa(2, 2) = myy(gi)
    kappa(2, 3) = myz(gi)
    kappa(3, 2) = myz(gi)
    kappa(3, 3) = mzz(gi)
  end function les_tensor

  function les_tensor_d(nvelocity, gvelocity, dn_t, ele, gi, nodes, dim, size_tensor) result(kappa)
    !! Get the derivative of the LES tensor at gauss point gi for element ele,
    !! with respect to dimension dim.
    type(vector_field), intent(in) :: gvelocity, nvelocity
    real, dimension(:, :, :), intent(in) :: dn_t
    real, dimension(3, 3) :: kappa
    real, dimension(3, 3), optional, intent(in) :: size_tensor
    integer, dimension(:), intent(in) :: nodes
    integer, intent(in) :: ele, gi, dim
    real, dimension(node_count(nvelocity)) :: nud, nvd, nwd

    real :: mxx, mxy, mxz, myy, myz, mzz ! ugh
    real :: mxxd, mxyd, mxzd, myyd, myzd, mzzd ! ugh
    real :: lenxx, lenxy, lenxz, lenyy, lenyz, lenzz ! ugh

    if (present(size_tensor)) then
      lenxx = size_tensor(1, 1)
      lenxy = size_tensor(1, 2)
      lenxz = size_tensor(1, 3)
      lenyy = size_tensor(2, 2)
      lenyz = size_tensor(2, 3)
      lenzz = size_tensor(3, 3)
    else
      call sizegieletens(dn_t(:, :, 1), dn_t(:, :, 2), dn_t(:, :, 3), ele_loc(nvelocity, ele), &
                    &  ele_ngi(nvelocity, ele), gi, lenxx, lenxy, lenxz, lenyy, lenyz, lenzz)
    end if

    ! set nud, nvd, nwd
    ! these tell the derivative routine what to differentiate against

    nud = 0.0
    nvd = 0.0
    nwd = 0.0

    select case(dim)
    case(1)
      nud(nodes) = 1.0
    case(2)
      nvd(nodes) = 1.0
    case(3)
      nvd(nodes) = 1.0
    end select

    call lesvis_d(node_count(nvelocity), ele_count(nvelocity), ele_loc(nvelocity, ele), &
              & nvelocity%mesh%ndglno, ele_ngi(nvelocity, ele), ele, gi, &
              & dn_t(:, :, 1), dn_t(:, :, 2), dn_t(:, :, 3), &
              & nvelocity%val(1,:), nud, nvelocity%val(2,:), nvd, nvelocity%val(3,:), nwd, & 
              & gvelocity%val(1,:), gvelocity%val(2,:), gvelocity%val(3,:), &
              & lenxx, lenxy, lenxz, lenyy, lenyz, lenzz, &
              & mxx, mxxd, mxy, mxyd, mxz, mxzd, myy, myyd, myz, myzd, mzz, mzzd) 

    kappa(1, 1) = mxxd
    kappa(1, 2) = mxyd
    kappa(2, 1) = mxyd
    kappa(1, 3) = mxzd
    kappa(3, 1) = mxzd
    kappa(2, 2) = myyd
    kappa(2, 3) = myzd
    kappa(3, 2) = myzd
    kappa(3, 3) = mzzd
  end function les_tensor_d

  subroutine goal_les_velocity_grad(state, dep, adj)
    ! The derivative of the goal with respect to
    ! whatever variables it depends on.
    ! Here it only depends on "NonlinearVelocity", so 
    ! other variables are ignored.
    type(state_type), dimension(:), intent(in) :: state
    character(len=*), intent(in) :: dep
    type(scalar_field), intent(inout) :: adj

    type(vector_field), pointer :: positions
    type(element_type), pointer :: t_shape
    type(scalar_field), pointer :: field
    integer :: ele, loc, ngi
    integer, dimension(:), pointer :: nodelist
    real, dimension(ele_ngi(adj, 1)) :: detwei
    real, dimension(ele_loc(adj, 1), ele_ngi(adj, 1), mesh_dim(adj)) :: dn_t
    real, dimension(ele_loc(adj, 1), ele_ngi(adj, 1)) :: tmpA, tmpB
    real, dimension(mesh_dim(adj), mesh_dim(adj)) :: kappa, kappa_d
    type(vector_field), pointer :: nvelocity, gvelocity
    integer :: idx, dim
    logical :: allocated
    type(tensor_field), pointer :: mesh_metric
    real, dimension(:, :, :), allocatable :: mesh_size_tensor_ele
    real, dimension(:, :), allocatable :: mesh_size_tensor
    integer :: stat, i

    call zero(adj)
    ! String parsing in Fortran
    ! is such a pain
    ! it's unreal
    idx = index(dep, "%")
    if (idx == 0) then
      return
    end if
    read(dep(idx+1:idx+1), *) dim

    if (dep(1:idx-1) == "NonlinearVelocity" .or. dep(1:idx-1) == "OldNonlinearVelocity") then
      field => extract_scalar_field(state(1), trim(dep), allocated=allocated)
    else
      return
    end if

    positions => extract_vector_field(state(1), "Coordinate")
    t_shape => ele_shape(adj, 1)
    nvelocity => extract_vector_field(state(1), "NonlinearVelocity")
    gvelocity => extract_vector_field(state(1), "GridVelocity")
    mesh_metric => extract_tensor_field(state(1), "MeshSizingMetric", stat=stat)

    allocate(mesh_size_tensor_ele(positions%dim, positions%dim, ele_loc(positions, 1)))
    allocate(mesh_size_tensor(positions%dim, positions%dim))
    do ele=1,ele_count(adj)
      nodelist => ele_nodes(adj, ele)
      call transform_to_physical(positions, ele, t_shape, dshape=dn_t, detwei=detwei)

      if (stat /= 0) then
        kappa_d = les_tensor_d(nvelocity, gvelocity, dn_t, ele, 1, nodelist, dim)
        kappa = les_tensor(nvelocity, gvelocity, dn_t, ele, 1)
      else
        mesh_size_tensor_ele = ele_val(mesh_metric, ele)
        mesh_size_tensor = 0.0
        do i=1,ele_loc(positions, 1)
          mesh_size_tensor = mesh_size_tensor + mesh_size_tensor_ele(:, :, i)
        end do
        mesh_size_tensor = mesh_size_tensor / ele_loc(positions, 1)
        call invert(mesh_size_tensor)

        kappa_d = les_tensor_d(nvelocity, gvelocity, dn_t, ele, 1, nodelist, dim, size_tensor=mesh_size_tensor)
        kappa = les_tensor(nvelocity, gvelocity, dn_t, ele, 1, size_tensor=mesh_size_tensor)
      end if

      do ngi=1,ele_ngi(adj,ele)
        do loc=1,ele_loc(adj, ele)
          tmpA(loc, ngi) = dot_product(dn_t(loc, ngi, :), matmul(kappa, dn_t(loc, ngi, :)))
          tmpB(loc, ngi) = dot_product(dn_t(loc, ngi, :), matmul(kappa_d, dn_t(loc, ngi, :)))
        end do
      end do

      ! tmpA now represents grad(N_i) . kappa . grad(N_i)
      ! tmpB now represents grad(N_i) . kappa_d . grad(N_i)
      
      do loc=1,ele_loc(adj, ele)
        tmpA(loc, :) = tmpA(loc, :) * node_val(field, nodelist(loc)) 
        tmpB(loc, :) = tmpB(loc, :) * 0.5 * (node_val(field, nodelist(loc)))**2
      end do

      ! tmpA now represents T_i . (grad(N_i) . kappa . grad(N_i))
      ! tmpB now represents 1/2 * grad(U_i) . kappa_d . grad(U_i)

      call addto(adj, ele_nodes(adj, ele), matmul(tmpA, detwei))
      call addto(adj, ele_nodes(adj, ele), matmul(tmpB, detwei))
    end do
    deallocate(mesh_size_tensor)
    deallocate(mesh_size_tensor_ele)

    if (allocated) then
      deallocate(field)
    end if
  end subroutine goal_les_velocity_grad

  subroutine goal_les_temp_grad(state, dep, adj)
    ! The derivative of the goal with respect to
    ! whatever variables it depends on.
    ! Here it only depends on "Temperature", so 
    ! other variables are ignored.
    type(state_type), dimension(:), intent(in) :: state
    character(len=*), intent(in) :: dep
    type(scalar_field), intent(inout) :: adj

    type(vector_field), pointer :: positions
    type(element_type), pointer :: t_shape
    type(scalar_field), pointer :: temp
    integer :: ele, loc, ngi
    integer, dimension(:), pointer :: nodelist
    real, dimension(ele_ngi(adj, 1)) :: detwei
    real, dimension(ele_loc(adj, 1), ele_ngi(adj, 1), mesh_dim(adj)) :: dn_t
    real, dimension(ele_loc(adj, 1), ele_ngi(adj, 1)) :: tmpA
    real, dimension(mesh_dim(adj), mesh_dim(adj)) :: kappa
    type(vector_field), pointer :: nvelocity, gvelocity
    type(tensor_field), pointer :: mesh_metric
    real, dimension(:, :, :), allocatable :: mesh_size_tensor_ele
    real, dimension(:, :), allocatable :: mesh_size_tensor
    integer :: stat, i

    call zero(adj)
    if (trim(dep) == "Temperature") then
      temp => extract_scalar_field(state(1), "Temperature")
    else if (trim(dep) == "OldTemperature") then
      temp => extract_scalar_field(state(1), "OldTemperature")
    else
      return
    end if

    positions => extract_vector_field(state(1), "Coordinate")
    t_shape => ele_shape(adj, 1)
    nvelocity => extract_vector_field(state(1), "NonlinearVelocity")
    gvelocity => extract_vector_field(state(1), "GridVelocity")
    mesh_metric => extract_tensor_field(state(1), "MeshSizingMetric", stat=stat)


    if (stat /= 0) then ! MeshSizingMetric not there
      do ele=1,ele_count(adj)
        nodelist => ele_nodes(adj, ele)
        call transform_to_physical(positions, ele, t_shape, dshape=dn_t, detwei=detwei)
        do ngi=1,ele_ngi(adj,ele)

          kappa = les_tensor(nvelocity, gvelocity, dn_t, ele, ngi)

          do loc=1,ele_loc(adj, ele)
            tmpA(loc, ngi) = dot_product(dn_t(loc, ngi, :), matmul(kappa, dn_t(loc, ngi, :)))
          end do
        end do

        ! tmpA now represents grad(N_i) . kappa . grad(N_i)
        
        do loc=1,ele_loc(adj, ele)
          tmpA(loc, :) = tmpA(loc, :) * node_val(temp, nodelist(loc)) 
        end do

        ! tmpA now represents T_i . (grad(N_i) . kappa . grad(N_i))

        call addto(adj, ele_nodes(adj, ele), matmul(tmpA, detwei))
      end do
    else
      allocate(mesh_size_tensor_ele(positions%dim, positions%dim, ele_loc(positions, 1)))
      allocate(mesh_size_tensor(positions%dim, positions%dim))
      do ele=1,ele_count(adj)
        nodelist => ele_nodes(adj, ele)
        call transform_to_physical(positions, ele, t_shape, dshape=dn_t, detwei=detwei)

        mesh_size_tensor_ele = ele_val(mesh_metric, ele)
        mesh_size_tensor = 0.0
        do i=1,ele_loc(positions, 1)
          mesh_size_tensor = mesh_size_tensor + mesh_size_tensor_ele(:, :, i)
        end do
        mesh_size_tensor = mesh_size_tensor / ele_loc(positions, 1)
        call invert(mesh_size_tensor)

        do ngi=1,ele_ngi(adj,ele)

          kappa = les_tensor(nvelocity, gvelocity, dn_t, ele, ngi, size_tensor=mesh_size_tensor)

          do loc=1,ele_loc(adj, ele)
            tmpA(loc, ngi) = dot_product(dn_t(loc, ngi, :), matmul(kappa, dn_t(loc, ngi, :)))
          end do
        end do

        ! tmpA now represents grad(N_i) . kappa . grad(N_i)
        
        do loc=1,ele_loc(adj, ele)
          tmpA(loc, :) = tmpA(loc, :) * node_val(temp, nodelist(loc)) 
        end do

        ! tmpA now represents T_i . (grad(N_i) . kappa . grad(N_i))

        call addto(adj, ele_nodes(adj, ele), matmul(tmpA, detwei))
      end do
      deallocate(mesh_size_tensor)
      deallocate(mesh_size_tensor_ele)
    end if

  end subroutine goal_les_temp_grad

  subroutine sizegieletens(nx,ny,nz,nloc,ngi,gi, &
 &     tensxx,tensxy,tensxz,tensyy,tensyz,tenszz )
  integer nloc,ngi
  integer gi
  real tensxx,tensxy,tensxz,tensyy,tensyz,tenszz
  real nx(nloc,ngi),ny(nloc,ngi),nz(nloc,ngi)
  
  real rn
  real udl,vdl,wdl
  integer iloc

  tensxx=0.0
  tensxy=0.0
  tensxz=0.0
  tensyy=0.0
  tensyz=0.0
  tenszz=0.0
     do 350 iloc=1,nloc
              rn=nx(iloc,gi)**2+ny(iloc,gi)**2+nz(iloc,gi)**2
              udl=1.*nx(iloc,gi)/rn
              vdl=1.*ny(iloc,gi)/rn 
              wdl=1.*nz(iloc,gi)/rn
              
              tensxx=tensxx + udl*udl
              tensxy=tensxy + udl*vdl
              tensxz=tensxz + udl*wdl
              tensyy=tensyy + vdl*vdl
              tensyz=tensyz + vdl*wdl
              tenszz=tenszz + wdl*wdl
350     continue
  return 
  end subroutine

     subroutine lesvis(nonods,totele,nloc,vondgl, ngi, ele,gi, &
 &             nx,ny,nz, &
 &             nu,nv,nw, ug,vg,wg, &
 &             mxxdtw,mxydtw,mxzdtw, & 
 &             myydtw,myzdtw,mzzdtw) 
     integer nonods,totele,nloc,ngi, ele,gi
     integer vondgl(totele*nloc)
     real nx(nloc,ngi),ny(nloc,ngi),nz(nloc,ngi)
     real nu(nonods),nv(nonods),nw(nonods)
     real ug(nonods),vg(nonods),wg(nonods)
     real mxxdtw(ngi),mxydtw(ngi),mxzdtw(ngi),myydtw(ngi)
     real myzdtw(ngi),mzzdtw(ngi)
     real sxxgi,sxygi,sxzgi,syygi,syzgi,szzgi,syxgi,szxgi,szygi
     real vis,cs2,fourcs
     real velud,velvd,velwd
     integer l,iglv
       sxxgi=0.
       sxygi=0.
       sxzgi=0.
       syygi=0.
       syzgi=0.
       szzgi=0.
       do 374 l=1,nloc
       iglv=vondgl((ele-1)*nloc+l)
         velud=nu(iglv)-ug(iglv)
         velvd=nv(iglv)-vg(iglv)
         velwd=nw(iglv)-wg(iglv)
         sxxgi = sxxgi + nx(l,gi)*velud
         sxygi = sxygi + 0.5*(ny(l,gi)*velud+nx(l,gi)*velvd)
         sxzgi = sxzgi + 0.5*(nz(l,gi)*velud+nx(l,gi)*velwd)
         syygi = syygi + ny(l,gi)*velvd
         syzgi = syzgi + 0.5*(nz(l,gi)*velvd+ny(l,gi)*velwd)
         szzgi = szzgi + nz(l,gi)*velwd
374        continue
       syxgi=sxygi
       szxgi=sxzgi
       szygi=syzgi
       vis=sqrt(2.* (sxxgi*sxxgi + sxygi*sxygi + sxzgi*sxzgi &
 &                +  syxgi*syxgi + syygi*syygi + syzgi*syzgi &
 &                +  szxgi*szxgi + szygi*szygi + szzgi*szzgi))
       cs2=0.1**2
       fourcs=4.*cs2
       mxxdtw(gi) = fourcs*vis*mxxdtw(gi)
       mxydtw(gi) = fourcs*vis*mxydtw(gi)
       mxzdtw(gi) = fourcs*vis*mxzdtw(gi)
       myydtw(gi) = fourcs*vis*myydtw(gi)
       myzdtw(gi) = fourcs*vis*myzdtw(gi)
       mzzdtw(gi) = fourcs*vis*mzzdtw(gi)
     return 
     end subroutine

!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 2.2.3 (r2178) - 11/14/2007 15:07
!  
!  Differentiation of lesvis in forward (tangent) mode:
!   variations  of output variables: mxzdtw mxydtw mxxdtw myzdtw
!                myydtw mzzdtw
!   with respect to input variables: nu nv nw
subroutine lesvis_d(nonods, totele, nloc, vondgl, ngi, ele, gi, nx, ny, &
&  nz, nu, nud, nv, nvd, nw, nwd, ug, vg, wg, lenxx, lenxy, lenxz, lenyy&
&  , lenyz, lenzz, mxxdtw, mxxdtwd, mxydtw, mxydtwd, mxzdtw, mxzdtwd, &
&  myydtw, myydtwd, myzdtw, myzdtwd, mzzdtw, mzzdtwd)
  implicit none
  integer, intent(in) :: nonods, totele, nloc, ngi, ele, gi
  integer, intent(in) :: vondgl(totele*nloc)
  real, intent(in) :: nx(nloc, ngi), ny(nloc, ngi), nz(nloc, ngi)
  real, intent(in) :: nu(nonods), nv(nonods), nw(nonods)
  real, intent(in) :: nud(nonods), nvd(nonods), nwd(nonods)
  real, intent(in) :: ug(nonods), vg(nonods), wg(nonods)
  real, intent(in) :: lenxx, lenxy, lenxz, lenyy
  real, intent(in) :: lenyz, lenzz
  real, intent(out) :: mxxdtw, mxydtw, mxzdtw, myydtw
  real, intent(out) :: mxxdtwd, mxydtwd, mxzdtwd, myydtwd
  real, intent(out) :: myzdtw, mzzdtw
  real, intent(out) :: myzdtwd, mzzdtwd
  real :: sxxgi, sxygi, sxzgi, syygi, syzgi, szzgi, syxgi, szxgi, szygi
  real :: sxxgid, sxygid, sxzgid, syygid, syzgid, szzgid, syxgid, szxgid&
&  , szygid
  real :: vis, cs2, fourcs
  real :: visd
  real :: velud, velvd, velwd
  real :: veludd, velvdd, velwdd
  integer :: l, iglv
  real :: arg1
  real :: arg1d
  intrinsic sqrt
  sxxgi = 0.
  sxygi = 0.
  sxzgi = 0.
  syygi = 0.
  syzgi = 0.
  szzgi = 0.
  syygid = 0.0
  sxygid = 0.0
  sxxgid = 0.0
  szzgid = 0.0
  syzgid = 0.0
  sxzgid = 0.0
  do l=1,nloc
    iglv = vondgl((ele-1)*nloc+l)
    veludd = nud(iglv)
    velud = nu(iglv) - ug(iglv)
    velvdd = nvd(iglv)
    velvd = nv(iglv) - vg(iglv)
    velwdd = nwd(iglv)
    velwd = nw(iglv) - wg(iglv)
    sxxgid = sxxgid + nx(l, gi)*veludd
    sxxgi = sxxgi + nx(l, gi)*velud
    sxygid = sxygid + 0.5*(ny(l, gi)*veludd+nx(l, gi)*velvdd)
    sxygi = sxygi + 0.5*(ny(l, gi)*velud+nx(l, gi)*velvd)
    sxzgid = sxzgid + 0.5*(nz(l, gi)*veludd+nx(l, gi)*velwdd)
    sxzgi = sxzgi + 0.5*(nz(l, gi)*velud+nx(l, gi)*velwd)
    syygid = syygid + ny(l, gi)*velvdd
    syygi = syygi + ny(l, gi)*velvd
    syzgid = syzgid + 0.5*(nz(l, gi)*velvdd+ny(l, gi)*velwdd)
    syzgi = syzgi + 0.5*(nz(l, gi)*velvd+ny(l, gi)*velwd)
    szzgid = szzgid + nz(l, gi)*velwdd
    szzgi = szzgi + nz(l, gi)*velwd
  end do
  syxgid = sxygid
  syxgi = sxygi
  szxgid = sxzgid
  szxgi = sxzgi
  szygid = syzgid
  szygi = syzgi
  arg1d = 2.*(sxxgid*sxxgi+sxxgi*sxxgid+sxygid*sxygi+sxygi*sxygid+sxzgid&
&    *sxzgi+sxzgi*sxzgid+syxgid*syxgi+syxgi*syxgid+syygid*syygi+syygi*&
&    syygid+syzgid*syzgi+syzgi*syzgid+szxgid*szxgi+szxgi*szxgid+szygid*&
&    szygi+szygi*szygid+szzgid*szzgi+szzgi*szzgid)
  arg1 = 2.*(sxxgi*sxxgi+sxygi*sxygi+sxzgi*sxzgi+syxgi*syxgi+syygi*syygi&
&    +syzgi*syzgi+szxgi*szxgi+szygi*szygi+szzgi*szzgi)
  if (arg1 .eq. 0.0) then
    visd = 0.0
  else
    visd = arg1d/(2.0*sqrt(arg1))
  end if
  vis = sqrt(arg1)
  cs2 = 0.1**2
  fourcs = 4.*cs2
  mxxdtwd = fourcs*lenxx*visd
  mxxdtw = fourcs*vis*lenxx
  mxydtwd = fourcs*lenxy*visd
  mxydtw = fourcs*vis*lenxy
  mxzdtwd = fourcs*lenxz*visd
  mxzdtw = fourcs*vis*lenxz
  myydtwd = fourcs*lenyy*visd
  myydtw = fourcs*vis*lenyy
  myzdtwd = fourcs*lenyz*visd
  myzdtw = fourcs*vis*lenyz
  mzzdtwd = fourcs*lenzz*visd
  mzzdtw = fourcs*vis*lenzz
  return
end subroutine lesvis_d

  function goal_les_velocity_4th(state)
    ! 0.5 * int( transpose(grad(T))   . kappa . grad(T) dv) -  
    ! 0.5 * int( transpose(grad(T)_h) . kappa . grad(T)_h dV).
    ! kappa is the LES matrix.
    type(state_type), dimension(:), intent(in) :: state
    real :: goal_les_velocity_4th

    integer :: ele, ngi
    real, dimension(:), allocatable :: detwei, result
    real, dimension(:, :, :), allocatable :: derivs
    real, dimension(:, :, :), allocatable :: dn_t
    type(vector_field), pointer :: positions
    type(element_type), pointer ::  t_shape
    real, dimension(:, :), allocatable :: kappa
    real, dimension(:), allocatable :: vec
    type(vector_field), pointer :: nvelocity, gvelocity
    integer :: i
    type(tensor_field), pointer :: mesh_metric
    real, dimension(:, :, :), allocatable :: mesh_size_tensor_ele
    real, dimension(:, :), allocatable :: mesh_size_tensor
    type(vector_field), dimension(3) :: grad_h_vel
    integer :: stat

    positions => extract_vector_field(state(1), "Coordinate")
    nvelocity => extract_vector_field(state(1), "NonlinearVelocity")
    gvelocity => extract_vector_field(state(1), "GridVelocity")
    mesh_metric => extract_tensor_field(state(1), "MeshSizingMetric", stat=stat)

    allocate(detwei(ele_ngi(nvelocity, 1)), result(ele_ngi(nvelocity, 1)), derivs(positions%dim, positions%dim, ele_ngi(nvelocity, 1)))
    allocate(dn_t(ele_loc(nvelocity, 1), ele_ngi(nvelocity, 1), positions%dim))
    allocate(kappa(positions%dim, positions%dim))
    allocate(vec(positions%dim))

    do i=1,nvelocity%dim
      call allocate(grad_h_vel(i), nvelocity%dim, nvelocity%mesh, "GradientNonlinearVelocity")
    end do
    call grad(nvelocity, positions, grad_h_vel)

    t_shape => ele_shape(nvelocity, 1)

    goal_les_velocity_4th = 0.0

    if (stat /= 0) then ! MeshSizingMetric not found
      do ele=1,ele_count(nvelocity)
        call transform_to_physical(positions, ele, t_shape, dshape=dn_t, detwei=detwei)

        do i=1,positions%dim
          derivs(i, :, :) = ele_val_at_quad(grad_h_vel(i), ele)
        end do

        do ngi=1,ele_ngi(nvelocity, 1)
          kappa = les_tensor(nvelocity, gvelocity, dn_t, ele, ngi)
          result(ngi) = 0.0
          do i=1,nvelocity%dim
            vec = matmul(ele_val(nvelocity, i, ele), dn_t(:, ngi, :))
            result(ngi) = result(ngi) + dot_product(vec, matmul(kappa, vec))
          end do

          do i=1,nvelocity%dim
            vec = derivs(i, :, ngi)
            result(ngi) = result(ngi) - dot_product(vec, matmul(kappa, vec))
          end do
        end do

        goal_les_velocity_4th = goal_les_velocity_4th + dot_product(result, detwei)
      end do
    else
      allocate(mesh_size_tensor_ele(positions%dim, positions%dim, ele_loc(positions, 1)))
      allocate(mesh_size_tensor(positions%dim, positions%dim))
      do ele=1,ele_count(nvelocity)
        call transform_to_physical(positions, ele, t_shape, dshape=dn_t, detwei=detwei)

        mesh_size_tensor_ele = ele_val(mesh_metric, ele)
        mesh_size_tensor = 0.0
        do i=1,ele_loc(positions, 1)
          mesh_size_tensor = mesh_size_tensor + mesh_size_tensor_ele(:, :, i)
        end do
        mesh_size_tensor = mesh_size_tensor / ele_loc(positions, 1)
        call invert(mesh_size_tensor)

        do i=1,positions%dim
          derivs(i, :, :) = ele_val_at_quad(grad_h_vel(i), ele)
        end do

        do ngi=1,ele_ngi(nvelocity, 1)
          kappa = les_tensor(nvelocity, gvelocity, dn_t, ele, ngi, size_tensor=mesh_size_tensor)
          result(ngi) = 0.0
          do i=1,nvelocity%dim
            vec = matmul(ele_val(nvelocity, i, ele), dn_t(:, ngi, :))
            result(ngi) = result(ngi) + dot_product(vec, matmul(kappa, vec))
          end do

          do i=1,nvelocity%dim
            vec = derivs(i, :, ngi)
            result(ngi) = result(ngi) - dot_product(vec, matmul(kappa, vec))
          end do
        end do

        goal_les_velocity_4th = goal_les_velocity_4th + dot_product(result, detwei)
      end do
      deallocate(mesh_size_tensor)
      deallocate(mesh_size_tensor_ele)
    end if

    goal_les_velocity_4th = goal_les_velocity_4th * 0.5

    deallocate(detwei)
    deallocate(dn_t)
    deallocate(kappa)
    deallocate(result)
    deallocate(vec)
    deallocate(derivs)
    do i=1,nvelocity%dim
      call deallocate(grad_h_vel(i))
    end do
  end function goal_les_velocity_4th

  subroutine goal_les_velocity_4th_grad(state, dep, adj)
    ! The derivative of the goal with respect to
    ! whatever variables it depends on.
    ! Here it only depends on "NonlinearVelocity", so 
    ! other variables are ignored.
    type(state_type), dimension(:), intent(in) :: state
    character(len=*), intent(in) :: dep
    type(scalar_field), intent(inout) :: adj

    type(vector_field), pointer :: positions
    type(element_type), pointer ::  t_shape
    type(scalar_field), pointer :: field
    integer :: ele, loc, ngi
    integer, dimension(:), pointer :: nodelist
    real, dimension(ele_ngi(adj, 1)) :: detwei
    real, dimension(ele_loc(adj, 1), ele_ngi(adj, 1), mesh_dim(adj)) :: dn_t
    real, dimension(ele_loc(adj, 1), ele_ngi(adj, 1)) :: tmpA, tmpB
    real, dimension(mesh_dim(adj), mesh_dim(adj)) :: kappa, kappa_d
    real, dimension(ele_loc(adj, 1), ele_loc(adj, 1)) :: ss
    real, dimension(ele_loc(adj, 1)) :: result, ones
    type(vector_field), pointer :: nvelocity, gvelocity
    integer :: idx, dim
    logical :: allocated
    type(tensor_field), pointer :: mesh_metric
    real, dimension(:, :, :), allocatable :: mesh_size_tensor_ele
    real, dimension(:, :), allocatable :: mesh_size_tensor
    type(vector_field) :: grad_h_vel
    integer :: stat, i
    real, dimension(mesh_dim(adj), ele_loc(adj, 1), ele_loc(adj, 1)) :: r
    real, dimension(mesh_dim(adj), ele_loc(adj, 1)) :: r_grad_ele
    type(scalar_field) :: lumped_mass


    call zero(adj)
    ! String parsing in Fortran
    ! is such a pain
    ! it's unreal
    idx = index(dep, "%")
    if (idx == 0) then
      return
    end if
    read(dep(idx+1:len(dep)), *) dim

    if (dep(1:idx-1) == "NonlinearVelocity" .or. dep(1:idx-1) == "OldNonlinearVelocity") then
      field => extract_scalar_field(state(1), trim(dep), allocated=allocated)
    else
      return
    end if

    positions => extract_vector_field(state(1), "Coordinate")
    t_shape => ele_shape(adj, 1)
    nvelocity => extract_vector_field(state(1), "NonlinearVelocity")
    gvelocity => extract_vector_field(state(1), "GridVelocity")
    mesh_metric => extract_tensor_field(state(1), "MeshSizingMetric", stat=stat)
    dim = nvelocity%dim
    ones = 1.0

    call allocate(grad_h_vel, nvelocity%dim, nvelocity%mesh, "GradientNonlinearVelocity")
    call allocate(lumped_mass, nvelocity%mesh, "LumpedMass")
    call grad(field, positions, grad_h_vel)
    call compute_lumped_mass(positions, lumped_mass)

    if (stat /= 0) then ! MeshSizingMetric not there
      do ele=1,ele_count(adj)
        nodelist => ele_nodes(adj, ele)
        call transform_to_physical(positions, ele, t_shape, dshape=dn_t, detwei=detwei)
        kappa_d = les_tensor_d(nvelocity, gvelocity, dn_t, ele, 1, nodelist, dim)
        kappa = les_tensor(nvelocity, gvelocity, dn_t, ele, 1)
        ss = shape_shape(t_shape, t_shape, detwei)

        do ngi=1,ele_ngi(adj,ele)
          do loc=1,ele_loc(adj, ele)
            tmpA(loc, ngi) = dot_product(dn_t(loc, ngi, :), matmul(kappa, dn_t(loc, ngi, :)))
            tmpB(loc, ngi) = dot_product(dn_t(loc, ngi, :), matmul(kappa_d, dn_t(loc, ngi, :)))
          end do
        end do

        ! tmpA now represents grad(N_i) . kappa . grad(N_i)
        ! tmpB now represents grad(N_i) . kappa_d . grad(N_i)
        
        do loc=1,ele_loc(adj, ele)
          tmpA(loc, :) = tmpA(loc, :) * node_val(field, nodelist(loc)) 
          tmpB(loc, :) = tmpB(loc, :) * 0.5 * (node_val(field, nodelist(loc)))**2
        end do

        ! tmpA now represents T_i . (grad(N_i) . kappa . grad(N_i))
        ! tmpB now represents 1/2 * grad(U_i) . kappa_d . grad(U_i)

        call addto(adj, ele_nodes(adj, ele), matmul(tmpA, detwei))
        call addto(adj, ele_nodes(adj, ele), matmul(tmpB, detwei))

        ! And now the terms that give the higher order.
        do loc=1,ele_loc(adj, ele)
          result(loc) = ss(loc, loc) * dot_product(node_val(grad_h_vel, nodelist(loc)), &
                                         & matmul(kappa_d, node_val(grad_h_vel, nodelist(loc))))
        end do

        call addto(adj, ele_nodes(adj, ele), -0.5 * result)

        r = shape_dshape(t_shape, dn_t, detwei)
        r_grad_ele = tensormul(r, ones, 3)

        do i=1,size(r_grad_ele, 1)
          r_grad_ele(i, :) = r_grad_ele(i, :) / ele_val(lumped_mass, ele)
        end do

        do loc=1,ele_loc(adj, ele)
          result(loc) = ss(loc, loc) * dot_product(node_val(grad_h_vel, nodelist(loc)), &
                                         & matmul(kappa, r_grad_ele(:, loc)))
        end do

        call addto(adj, ele_nodes(adj, ele), -1.0 * result)
      end do
    else
      allocate(mesh_size_tensor_ele(positions%dim, positions%dim, ele_loc(positions, 1)))
      allocate(mesh_size_tensor(positions%dim, positions%dim))
      do ele=1,ele_count(adj)
        nodelist => ele_nodes(adj, ele)
        call transform_to_physical(positions, ele, t_shape, dshape=dn_t, detwei=detwei)

        mesh_size_tensor_ele = ele_val(mesh_metric, ele)
        mesh_size_tensor = 0.0
        do i=1,ele_loc(positions, 1)
          mesh_size_tensor = mesh_size_tensor + mesh_size_tensor_ele(:, :, i)
        end do
        mesh_size_tensor = mesh_size_tensor / ele_loc(positions, 1)
        call invert(mesh_size_tensor)

        kappa_d = les_tensor_d(nvelocity, gvelocity, dn_t, ele, 1, nodelist, dim, size_tensor=mesh_size_tensor)
        kappa = les_tensor(nvelocity, gvelocity, dn_t, ele, 1, size_tensor=mesh_size_tensor)
        ss = shape_shape(t_shape, t_shape, detwei)

        do ngi=1,ele_ngi(adj,ele)
          do loc=1,ele_loc(adj, ele)
            tmpA(loc, ngi) = dot_product(dn_t(loc, ngi, :), matmul(kappa, dn_t(loc, ngi, :)))
            tmpB(loc, ngi) = dot_product(dn_t(loc, ngi, :), matmul(kappa_d, dn_t(loc, ngi, :)))
          end do
        end do

        ! tmpA now represents grad(N_i) . kappa . grad(N_i)
        ! tmpB now represents grad(N_i) . kappa_d . grad(N_i)
        
        do loc=1,ele_loc(adj, ele)
          tmpA(loc, :) = tmpA(loc, :) * node_val(field, nodelist(loc)) 
          tmpB(loc, :) = tmpB(loc, :) * 0.5 * (node_val(field, nodelist(loc)))**2
        end do

        ! tmpA now represents T_i . (grad(N_i) . kappa . grad(N_i))
        ! tmpB now represents 1/2 * grad(U_i) . kappa_d . grad(U_i)

        call addto(adj, ele_nodes(adj, ele), matmul(tmpA, detwei))
        call addto(adj, ele_nodes(adj, ele), matmul(tmpB, detwei))

        ! And now the terms that give the higher order.
        do loc=1,ele_loc(adj, ele)
          result(loc) = ss(loc, loc) * dot_product(node_val(grad_h_vel, nodelist(loc)), &
                                         & matmul(kappa_d, node_val(grad_h_vel, nodelist(loc))))
        end do

        call addto(adj, ele_nodes(adj, ele), -0.5 * result)

        r = shape_dshape(t_shape, dn_t, detwei)
        r_grad_ele = tensormul(r, ones, 3)

        do i=1,size(r_grad_ele, 1)
          r_grad_ele(i, :) = r_grad_ele(i, :) / ele_val(lumped_mass, ele)
        end do

        do loc=1,ele_loc(adj, ele)
          result(loc) = ss(loc, loc) * dot_product(node_val(grad_h_vel, nodelist(loc)), &
                                         & matmul(kappa, r_grad_ele(:, loc)))
        end do

        call addto(adj, ele_nodes(adj, ele), -1.0 * result)
      end do
      deallocate(mesh_size_tensor)
      deallocate(mesh_size_tensor_ele)
    end if

    call deallocate(grad_h_vel)
    call deallocate(lumped_mass)
    if (allocated) then
      deallocate(field)
    end if
  end subroutine goal_les_velocity_4th_grad

  function goal_les_velocity_squared(state)
    ! 0.5 * int( |kappa . grad(T)|^2 ) dV.
    ! kappa is the LES matrix.
    type(state_type), dimension(:), intent(in) :: state
    real :: goal_les_velocity_squared

    integer :: ele, ngi
    real, dimension(:), allocatable :: detwei, result
    real, dimension(:, :, :), allocatable :: dn_t
    type(vector_field), pointer :: positions
    type(element_type), pointer ::  t_shape
    real, dimension(:, :), allocatable :: kappa
    real, dimension(:), allocatable :: vec
    type(vector_field), pointer :: nvelocity, gvelocity
    integer :: i
    type(tensor_field), pointer :: mesh_metric
    real, dimension(:, :, :), allocatable :: mesh_size_tensor_ele
    real, dimension(:, :), allocatable :: mesh_size_tensor
    integer :: stat

    positions => extract_vector_field(state(1), "Coordinate")
    nvelocity => extract_vector_field(state(1), "NonlinearVelocity")
    gvelocity => extract_vector_field(state(1), "GridVelocity")
    mesh_metric => extract_tensor_field(state(1), "MeshSizingMetric", stat=stat)

    allocate(detwei(ele_ngi(nvelocity, 1)), result(ele_ngi(nvelocity, 1)))
    allocate(dn_t(ele_loc(nvelocity, 1), ele_ngi(nvelocity, 1), positions%dim))
    allocate(kappa(positions%dim, positions%dim))
    allocate(vec(positions%dim))
    t_shape => ele_shape(nvelocity, 1)

    goal_les_velocity_squared = 0.0

    if (stat /= 0) then ! MeshSizingMetric not found
      do ele=1,ele_count(nvelocity)
        call transform_to_physical(positions, ele, t_shape, dshape=dn_t, detwei=detwei)


        do ngi=1,ele_ngi(nvelocity, 1)
          kappa = les_tensor(nvelocity, gvelocity, dn_t, ele, ngi)
          result(ngi) = 0.0
          do i=1,nvelocity%dim
            vec = matmul(ele_val(nvelocity, i, ele), dn_t(:, ngi, :))
            vec = matmul(kappa, vec)
            result(ngi) = result(ngi) + norm2(vec)**2
          end do
        end do

        goal_les_velocity_squared = goal_les_velocity_squared + dot_product(result, detwei)
      end do
    else
      allocate(mesh_size_tensor_ele(positions%dim, positions%dim, ele_loc(positions, 1)))
      allocate(mesh_size_tensor(positions%dim, positions%dim))
      do ele=1,ele_count(nvelocity)
        call transform_to_physical(positions, ele, t_shape, dshape=dn_t, detwei=detwei)

        mesh_size_tensor_ele = ele_val(mesh_metric, ele)
        mesh_size_tensor = 0.0
        do i=1,ele_loc(positions, 1)
          mesh_size_tensor = mesh_size_tensor + mesh_size_tensor_ele(:, :, i)
        end do
        mesh_size_tensor = mesh_size_tensor / ele_loc(positions, 1)
        call invert(mesh_size_tensor)

        do ngi=1,ele_ngi(nvelocity, 1)
          kappa = les_tensor(nvelocity, gvelocity, dn_t, ele, ngi, size_tensor=mesh_size_tensor)
          result(ngi) = 0.0
          do i=1,nvelocity%dim
            vec = matmul(ele_val(nvelocity, i, ele), dn_t(:, ngi, :))
            vec = matmul(kappa, vec)
            result(ngi) = result(ngi) + norm2(vec)**2
          end do
        end do

        goal_les_velocity_squared = goal_les_velocity_squared + dot_product(result, detwei)
      end do
      deallocate(mesh_size_tensor)
      deallocate(mesh_size_tensor_ele)
    end if

    goal_les_velocity_squared = goal_les_velocity_squared * 0.5

    deallocate(detwei)
    deallocate(result)
    deallocate(dn_t)
    deallocate(kappa)
    deallocate(vec)
  end function goal_les_velocity_squared

  subroutine goal_les_velocity_squared_grad(state, dep, adj)
    ! The derivative of the goal with respect to
    ! whatever variables it depends on.
    ! Here it only depends on "NonlinearVelocity", so 
    ! other variables are ignored.
    type(state_type), dimension(:), intent(in) :: state
    character(len=*), intent(in) :: dep
    type(scalar_field), intent(inout) :: adj

    type(vector_field), pointer :: positions
    type(element_type), pointer :: t_shape
    type(scalar_field), pointer :: field
    integer :: ele, loc, ngi
    integer, dimension(:), pointer :: nodelist
    real, dimension(ele_ngi(adj, 1)) :: detwei
    real, dimension(ele_loc(adj, 1), ele_ngi(adj, 1), mesh_dim(adj)) :: dn_t
    real, dimension(ele_loc(adj, 1), ele_ngi(adj, 1)) :: tmpA, tmpB
    real, dimension(mesh_dim(adj), mesh_dim(adj)) :: kappa, kappa_d
    real, dimension(mesh_dim(adj)) :: tmp_vector
    type(vector_field), pointer :: nvelocity, gvelocity
    integer :: idx, dim
    logical :: allocated
    type(tensor_field), pointer :: mesh_metric
    real, dimension(:, :, :), allocatable :: mesh_size_tensor_ele
    real, dimension(:, :), allocatable :: mesh_size_tensor
    integer :: stat, i

    call zero(adj)
    ! String parsing in Fortran
    ! is such a pain
    ! it's unreal
    idx = index(dep, "%")
    if (idx == 0) then
      return
    end if
    read(dep(idx+1:idx+1), *) dim

    if (dep(1:idx-1) == "NonlinearVelocity" .or. dep(1:idx-1) == "OldNonlinearVelocity") then
      field => extract_scalar_field(state(1), trim(dep), allocated=allocated)
    else
      return
    end if

    positions => extract_vector_field(state(1), "Coordinate")
    t_shape => ele_shape(adj, 1)
    nvelocity => extract_vector_field(state(1), "NonlinearVelocity")
    gvelocity => extract_vector_field(state(1), "GridVelocity")
    mesh_metric => extract_tensor_field(state(1), "MeshSizingMetric", stat=stat)

    allocate(mesh_size_tensor_ele(positions%dim, positions%dim, ele_loc(positions, 1)))
    allocate(mesh_size_tensor(positions%dim, positions%dim))
    do ele=1,ele_count(adj)
      nodelist => ele_nodes(adj, ele)
      call transform_to_physical(positions, ele, t_shape, dshape=dn_t, detwei=detwei)

      if (stat /= 0) then
        kappa_d = les_tensor_d(nvelocity, gvelocity, dn_t, ele, 1, nodelist, dim)
        kappa = les_tensor(nvelocity, gvelocity, dn_t, ele, 1)
      else
        mesh_size_tensor_ele = ele_val(mesh_metric, ele)
        mesh_size_tensor = 0.0
        do i=1,ele_loc(positions, 1)
          mesh_size_tensor = mesh_size_tensor + mesh_size_tensor_ele(:, :, i)
        end do
        mesh_size_tensor = mesh_size_tensor / ele_loc(positions, 1)
        call invert(mesh_size_tensor)

        kappa_d = les_tensor_d(nvelocity, gvelocity, dn_t, ele, 1, nodelist, dim, size_tensor=mesh_size_tensor)
        kappa = les_tensor(nvelocity, gvelocity, dn_t, ele, 1, size_tensor=mesh_size_tensor)
      end if

      do ngi=1,ele_ngi(adj,ele)
        do loc=1,ele_loc(adj, ele)
          tmp_vector = matmul(kappa, dn_t(loc, ngi, :))
          tmpA(loc, ngi) = dot_product(tmp_vector, tmp_vector)
          tmpB(loc, ngi) = dot_product(tmp_vector, matmul(kappa_d, dn_t(loc, ngi, :)))
        end do
      end do

      ! tmpA now represents (kappa x gradN) dot (kappa x gradN)
      ! tmpA now represents (kappa x gradN) dot (kappa_d x gradN)

      do loc=1,ele_loc(adj, ele)
        tmpA(loc, :) = tmpA(loc, :) * node_val(field, nodelist(loc)) * 2
        tmpB(loc, :) = tmpB(loc, :) * 2 * (node_val(field, nodelist(loc)))**2
      end do

      ! tmpA now represents 2 * u_i * tmpA
      ! tmpB now represents 2 * u_i**2 * tmpB

      call addto(adj, ele_nodes(adj, ele), matmul(tmpA, detwei))
      call addto(adj, ele_nodes(adj, ele), matmul(tmpB, detwei))
    end do
    deallocate(mesh_size_tensor)
    deallocate(mesh_size_tensor_ele)

    if (allocated) then
      deallocate(field)
    end if
  end subroutine goal_les_velocity_squared_grad
end module goals 
