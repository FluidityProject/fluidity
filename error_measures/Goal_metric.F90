!#define GOAL_DEBUG
#include "fdebug.h"

module goal_metric
!!< Simple goal-based metric formation.

  use global_parameters, only: dt, FIELD_NAME_LEN, OPTION_PATH_LEN
  use futils, only: count_chars, multiindex
  use elements
  use spud
  use transform_elements, only: transform_to_physical
  use fetools, only: shape_shape
  use metric_tools, only: error_bound_name
  use fields
  use edge_length_module, only: get_edge_lengths
  use state_module, only: state_type, extract_scalar_field, extract_vector_field, insert
  use merge_tensors, only: merge_tensor_fields
  use vtk_interfaces, only: vtk_write_fields
  use field_derivatives, only: compute_hessian, patch_type, get_patch_node
  use form_metric_field, only: form_metric
  use goals
  use gradation_metric, only: form_gradation_metric

  implicit none

  !! Use goal-based metric method, or
  !! the standard user-defined weights?
  logical :: use_goal_metric = .false.

  !! How much error should we tolerate in the goal?
  !! If goal_rel_tolerance /= 0.0, then use that as a relative
  !! tolerance, otherwise use goal_tolerance.
  real :: goal_tolerance
  real :: goal_rel_tolerance = 0.0
  character(len=OPTION_PATH_LEN) :: goal_name

  !! What state variables does the goal depend on?
  !! The gradient of the goal with respect to
  !! each dependency must be available from goal_grad.
  character(len=FIELD_NAME_LEN), dimension(:), pointer :: goal_deps => null()

  interface form_goal_metric
    module procedure form_goal_metric_generic, form_goal_metric_specific
  end interface

  contains

  subroutine initialise_goal_metric
    integer, dimension(:), allocatable :: indices
    character(len=OPTION_PATH_LEN) :: name, deps
    integer :: i, nchild, idx, c

    if (have_option("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity")) then
      ! We have the option. Let's set up goal-based adaptivity. Exciting!
      use_goal_metric = .true.

      ! Get the tolerance, either relative or absolute.
      if (have_option("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity/relative_tolerance")) then
        call get_option("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity/relative_tolerance", goal_rel_tolerance)
      else if (have_option("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity/absolute_tolerance")) then
        call get_option("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity/absolute_tolerance", goal_tolerance)
      else
        FLExit("Must specify either an absolute tolerance or a relative tolerance for goal based adaptivity.")
      end if

      ! We need to find out which goal we actually want ...
      call get_number_of_children("/mesh_adaptivity/hr_adaptivity/goal_based&
           &_adaptivity",nchild)
      do i=0,nchild-1
        call get_child_name("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity", i, name)
        idx = index(trim(name), "tolerance")
        if (idx == 0) then
          goal_name = trim(name)
          exit
        end if
      end do

      ! Once we have the goal, find out what fields in state it depends on
      if (.not. associated(goal_deps)) then
        call get_option("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity/" // trim(goal_name) &
                      & // "/dependencies", deps)

        ! deps == "NonlinearVelocity%1 NonlinearVelocity%2 NonlinearVelocity%3" for example

        c = count_chars(trim(deps), " ")
        allocate(goal_deps(c + 1))
        allocate(indices(c + 2))
        
        if (c == 0) then
          goal_deps(1) = trim(deps)
        else
          indices(1) = 0
          indices(2:c+1) = multiindex(trim(deps), " ")
          indices(c+2) = len_trim(deps) + 1

          do i=1,c+1
            goal_deps(i) = deps(indices(i) + 1:indices(i+1) - 1)
          end do
        end if
      end if
    end if

    if (allocated(indices)) then
      deallocate(indices)
    end if
  end subroutine initialise_goal_metric

  subroutine form_goal_metric_generic(state, metric)
    type(state_type), dimension(:), intent(inout) :: state
    type(tensor_field), intent(inout) :: metric

    ! number of nonlinear iterations to perform, when applicable
    integer :: iters, i
    type(tensor_field) :: tmp_metric
    type(vector_field), pointer :: positions
 
    select case(goal_name)
    case("enstrophy_goal")
      call form_goal_metric_specific(state, metric, goal_enstrophy, goal_enstrophy_grad)
    case("temperature_gradient_goal")
      call form_goal_metric_specific(state, metric, goal_temp, goal_temp_grad)
    case("les_goal")
      call get_option("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity/les_goal/nonlinear_iterations", iters, default=3)
      call form_goal_metric_specific(state, metric, goal_les_velocity, goal_les_velocity_grad)
      if (iters > 0) then
        positions => extract_vector_field(state(1), "Coordinate")
        call form_gradation_metric(positions, metric)
        call allocate(tmp_metric, metric%mesh, "TmpMetric")
        do i=1,iters
          call insert(state(1), metric, "MeshSizingMetric")
          call form_goal_metric(state, tmp_metric, goal_les_velocity, goal_les_velocity_grad)
          call form_gradation_metric(positions, tmp_metric)
          metric%val = tmp_metric%val
        end do
        call deallocate(tmp_metric)
      end if
    case("les_squared_goal")
      call get_option("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity/les_squared_goal/nonlinear_iterations", iters, default=5)
      call form_goal_metric_specific(state, metric, goal_les_velocity_squared, goal_les_velocity_squared_grad)
      if (iters > 0) then
        positions => extract_vector_field(state(1), "Coordinate")
        call form_gradation_metric(positions, metric)
        call allocate(tmp_metric, metric%mesh, "TmpMetric")
        do i=1,iters
          call insert(state(1), metric, "MeshSizingMetric")
          call form_goal_metric(state, tmp_metric, goal_les_velocity_squared, goal_les_velocity_squared_grad)
          call form_gradation_metric(positions, tmp_metric)
          metric%val = tmp_metric%val
        end do
        call deallocate(tmp_metric)
      end if
    case("higher_order_les_goal")
      call get_option("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity/les_goal/nonlinear_iterations", iters, default=3)
      call form_goal_metric_specific(state, metric, goal_les_velocity_4th, goal_les_velocity_4th_grad)
      if (iters > 0) then
        positions => extract_vector_field(state(1), "Coordinate")
        call form_gradation_metric(positions, metric)
        call allocate(tmp_metric, metric%mesh, "TmpMetric")
        do i=1,iters
          call insert(state(1), metric, "MeshSizingMetric")
          call form_goal_metric(state, tmp_metric, goal_les_velocity_4th, goal_les_velocity_4th_grad)
          call form_gradation_metric(positions, tmp_metric)
          metric%val = tmp_metric%val
        end do
        call deallocate(tmp_metric)
      end if
    end select
  end subroutine form_goal_metric_generic

  subroutine form_goal_metric_specific(state, metric, goal, goal_grad)
    type(state_type), dimension(:), intent(inout) :: state
    type(tensor_field), intent(inout) :: metric
    !! This function computes the goal, given the state
    !! of the system.
    interface 
      function goal(state)
        use state_module, only:state_type
        real :: goal
        type(state_type), dimension(:), intent(in) :: state
      end function goal
    end interface

    !! This function computes the gradient of the goal
    !! with respect to a particular dependency, at a particular point.
    interface
      subroutine goal_grad(state, dep, adj)
        use fields_data_types, only:scalar_field
        use state_module, only: state_type
        type(state_type), dimension(:), intent(in) :: state
        character(len=*), intent(in) :: dep
        type(scalar_field), intent(inout) :: adj
      end subroutine goal_grad
    end interface

    character(len=FIELD_NAME_LEN) :: dep
    integer :: i
    type(tensor_field) :: dep_hessian, adj_hessian, tmp_metric
    type(scalar_field) :: dep_err_field, adj_err_field, adj_field
    type(scalar_field), pointer :: dep_field
    type(mesh_type), pointer :: mesh
    type(vector_field), pointer :: positions, velocity
    integer :: stat
    logical :: allocated
    type(scalar_field) :: edgelen
    integer, save :: adaptcnt = 0
    character(len=20) :: buf

    ! At the moment all the fields have to be on the same mesh.
    positions => extract_vector_field(state(1), "Coordinate")
    velocity => extract_vector_field(state(1), "Velocity", stat=stat)
    if (stat == 0) then
      mesh => velocity%mesh
    else
      mesh => positions%mesh
    end if
    call allocate(dep_hessian, mesh, "Hessian")
    call allocate(adj_hessian, mesh, "Adjoint Hessian")
    call allocate(tmp_metric, mesh, "Accumulated Adjoint Hessian")
    call allocate(adj_field, mesh, "Adjoint")
    call allocate(dep_err_field, mesh, "Error field")
    call allocate(adj_err_field, mesh, "Adjoint error field")
    call allocate(edgelen, mesh, "Desired edge lengths")
    call zero(metric)
    call zero(tmp_metric)

    write(buf, '(i0)') adaptcnt

    if (goal_rel_tolerance /= 0.0) then
      goal_tolerance = goal_rel_tolerance * goal(state)
    end if

    ! Loop over dependent variables of goal.
    do i=1,size(goal_deps)
      dep = goal_deps(i)
      dep_field => extract_scalar_field(state(1), dep, allocated=allocated)

      ! Compute "adjoint" solution.
      call compute_adjoint(goal_grad, state, dep, adj_field)
#ifdef GOAL_DEBUG
      call vtk_write_fields("goal_adjoint_" // trim(buf), i, positions, mesh, sfields=(/dep_field, adj_field/))
#endif

      ! Compute Hessians.
      call compute_hessian(dep_field, positions, dep_hessian)
#ifdef GOAL_DEBUG
      call get_edge_lengths(dep_hessian, edgelen)
      call vtk_write_fields("goal_fhessian_" // trim(buf), i, positions, mesh, sfields=(/dep_field, edgelen/), tfields=(/dep_hessian/))
#endif
      call compute_hessian(adj_field, positions, adj_hessian)
#ifdef GOAL_DEBUG
      call get_edge_lengths(adj_hessian, edgelen)
      call vtk_write_fields("goal_bhessian_" // trim(buf), i, positions, mesh, sfields=(/adj_field, edgelen/), tfields=(/adj_hessian/))
#endif

      ! Compute forward error field.
      call compute_err_field(dep_field, dep_hessian, positions, adj_err_field)
#ifdef GOAL_DEBUG
      call vtk_write_fields("goal_berrfield_" // trim(buf), i, positions, mesh, sfields=(/adj_err_field/))
#endif
      call compute_err_field(adj_field, adj_hessian, positions, dep_err_field)
#ifdef GOAL_DEBUG
      call vtk_write_fields("goal_ferrfield_" // trim(buf), i, positions, mesh, sfields=(/dep_err_field/))
#endif

      ! Form metric.
      call insert(state(1), dep_err_field, error_bound_name(trim(dep_field%name)))
      call form_metric(state(1), dep_hessian, dep_field)
#ifdef GOAL_DEBUG
      call get_edge_lengths(dep_hessian, edgelen)
      call vtk_write_fields("goal_fmetric_" // trim(buf), i, positions, mesh, sfields=(/dep_field, edgelen/), tfields=(/dep_hessian/))
#endif
      call insert(state(1), adj_err_field, error_bound_name(trim(adj_field%name)))
      call form_metric(state(1), adj_hessian, adj_field)
#ifdef GOAL_DEBUG
      call get_edge_lengths(adj_hessian, edgelen)
      call vtk_write_fields("goal_bmetric_" // trim(buf), i, positions, mesh, sfields=(/adj_field, edgelen/), tfields=(/adj_hessian/))
#endif

      ! Merge.
      call merge_tensor_fields(metric, dep_hessian)
#ifdef GOAL_DEBUG
      call get_edge_lengths(metric, edgelen)
      call vtk_write_fields("final_goal_fmetric_" // trim(buf), i, positions, mesh, sfields=(/edgelen/), tfields=(/metric/))
#endif
      call merge_tensor_fields(tmp_metric, adj_hessian)
#ifdef GOAL_DEBUG
      call get_edge_lengths(tmp_metric, edgelen)
      call vtk_write_fields("final_goal_bmetric_" // trim(buf), i, positions, mesh, sfields=(/edgelen/), tfields=(/tmp_metric/))
#endif

      if (allocated) deallocate(dep_field)
    end do
    adaptcnt = adaptcnt + 1

    call merge_tensor_fields(metric, tmp_metric, aniso_min=.true.)
#ifdef GOAL_DEBUG
      call get_edge_lengths(metric, edgelen)
      call vtk_write_fields("final_goal_metric_" // trim(buf), i, positions, mesh, sfields=(/edgelen/), tfields=(/metric/))
#endif

    call deallocate(dep_hessian)
    call deallocate(adj_hessian)
    call deallocate(tmp_metric)
    call deallocate(adj_field)
    call deallocate(dep_err_field)
    call deallocate(adj_err_field)
    call deallocate(edgelen)
  end subroutine form_goal_metric_specific

  subroutine compute_adjoint(goal_grad, state, dep, adj_field)
    type(state_type), dimension(:), intent(in) :: state
    character(len=FIELD_NAME_LEN), intent(in) :: dep
    type(scalar_field), intent(inout) :: adj_field
    !! This function computes the gradient of the goal
    !! with respect to a particular dependency, at a particular point.
    interface
      subroutine goal_grad(state, dep, adj)
        use fields_data_types, only: scalar_field
        use state_module, only: state_type
        type(state_type), dimension(:), intent(in) :: state
        character(len=*), intent(in) :: dep
        type(scalar_field), intent(inout) :: adj
      end subroutine goal_grad
    end interface

    call goal_grad(state, dep, adj_field)
    adj_field%val = adj_field%val * dt
  end subroutine compute_adjoint

  subroutine compute_err_field(field, hessian, positions, err_field)
    type(scalar_field), intent(in) :: field
    type(tensor_field), intent(in) :: hessian
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(inout) :: err_field

    type(mesh_type) :: mesh
    type(patch_type) :: patch
    integer :: node, nnode, i
    type(element_type), pointer :: t_shape, x_shape
    real, dimension(hessian%dim(1), hessian%dim(2)) :: avg_hessian, evecs
    real, dimension(hessian%dim(1)) :: edge, evals

    real :: err

    mesh = field%mesh
    call add_nelist(mesh)
    t_shape => ele_shape(field, 1)
    x_shape => ele_shape(positions, 1)

    do node=1,node_count(mesh)
      patch = get_patch_node(mesh, node, level=1)
      err = 0.0

      ! estimate residual.
      do i=1,patch%count
        nnode = patch%elements(i)
        edge = node_val(positions, node) - node_val(positions, nnode)
        avg_hessian = (node_val(hessian, node) + node_val(hessian, nnode)) / 2.0
        call eigendecomposition_symmetric(avg_hessian, evecs, evals)
        evals = abs(evals)
        call eigenrecomposition(avg_hessian, evecs, evals)
        err = max(err, dot_product(edge, matmul(avg_hessian, edge)))
      end do

      if (err == 0.0) then
        err_field%val(node) = 1.0e10
      else 
        err_field%val(node) = goal_tolerance / (node_count(mesh) * err)
      end if
      deallocate(patch%elements)
    end do
  end subroutine compute_err_field

end module goal_metric
