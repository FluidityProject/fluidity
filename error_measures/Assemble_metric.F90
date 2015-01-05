#include "fdebug.h"

module metric_assemble

  use VTK_interfaces
  use surfacelabels
!  use fields
  use field_derivatives
  use form_metric_field
  use merge_tensors
  use edge_length_module
  use interpolation_metric
  use goals
  use goal_metric
  use gradation_metric
  use bounding_box_metric
  use boundary_metric
  use geometric_constraints_metric
  use limit_metric_module
  use metric_tools
  use state_module
  use halos
  use spud
  use parallel_tools
  use metric_advection
  use anisotropic_gradation
  use richardson_metric_module
  use anisotropic_zz_module
  use aspect_ratios_module
  use reference_meshes
  use hadapt_metric_based_extrude, only: get_1d_mesh, recombine_metric, get_1d_tensor
  
  implicit none

  private
  public :: assemble_metric, apply_vertical_gradation, apply_horizontal_gradation
  
  contains

  subroutine assemble_metric(state, error_metric)
    !!< This routine drives the metric assembly logic, which is contained in
    !!< other modules.

    type(state_type), dimension(:), intent(inout) :: state
    type(tensor_field), intent(inout) :: error_metric

    integer :: i, stat

    type(vector_field), pointer :: positions
    integer, save :: adaptcnt = 0
    character(len=20) :: buf
    logical :: debug_metric, vertically_structured_adaptivity, split_gradation
    type(tensor_field), pointer :: max_tensor

    ewrite(2,*) "+: Assembling metric"
    debug_metric = have_option("/mesh_adaptivity/hr_adaptivity/debug/write_metric_stages")
    ! is this metric going to be collapsed in the vertical to do horizontal adaptivity with it?
    vertically_structured_adaptivity = have_option("/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity")
    ! are we waiting until later to apply the gradation?
    split_gradation = have_option("/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity/split_gradation")

    do i=1,size(state)
      positions => extract_vector_field(state(i), "Coordinate", stat=stat)
      if (stat == 0) exit
    end do

    max_tensor => extract_tensor_field(state(1), "MinMetricEigenbound")

    if (debug_metric) then
      do i=1,size(state)
        write(buf, '(i0)') i 

        call vtk_write_state(trim("metric_input_") // trim(buf),&
             & adaptcnt, state=(/state(i)/)) 
      end do
      adaptcnt = adaptcnt + 1
    end if

    call zero(error_metric)

    call initialise_boundcount(error_metric%mesh, positions)

    call initialise_interpolation_metric
    call initialise_goal_metric
    call initialise_gradation_metric
    call initialise_anisotropic_gradation
    call initialise_bounding_box_metric
    call initialise_boundary_metric
    call initialise_geometric_constraints_metric
    call initialise_metric_advection
    call initialise_richardson_number_metric

    if (use_goal_metric) then
!     ***** Note that the following interface causes issues with the Intel compiler. Changed for now as a workaround ******
!      call form_goal_metric(state, error_metric)
      call form_goal_metric_generic(state, error_metric)
!     ***** End of Intel compiler workaround *****
      call halo_update(error_metric)
    else if (use_interpolation_metric) then
      call form_interpolation_metric(state, error_metric)
      call halo_update(error_metric)
      call form_anisotropic_zz_metric(state, error_metric)
    end if
    
    if (use_richardson_number_metric) then
      call form_richardson_number_metric(state, error_metric)
      call halo_update(error_metric)
    end if 
    
    if (use_boundary_metric) then
      call form_boundary_metric(error_metric, positions)
      call halo_update(error_metric)
    end if

    call enforce_reference_meshes(state, positions, error_metric)

    if (use_geometric_constraints_metric) then
      call form_geometric_constraints_metric(positions, error_metric, state(1))
      call halo_update(error_metric)
    end if

    ! gradation might be left until the metric is split into horizontal and vertical components
    if(.not.split_gradation) then
    
      call apply_gradation(error_metric, positions, state(1))

    end if

    if (use_metric_advection) then
      call form_advection_metric(error_metric, state(1))
      call halo_update(error_metric)
    end if

    if (use_bounding_box_metric) then
      call form_bounding_box_metric(positions, error_metric, max_tensor)
      call halo_update(error_metric)
    end if

    call bound_metric_aspect_ratio(error_metric)
    ! for vertically structured, the limiting should happen after the vertical collapsing
    if (.not. vertically_structured_adaptivity) then
      call limit_metric(positions, error_metric)
    end if
    call halo_update(error_metric)
    
  end subroutine assemble_metric
  
  subroutine apply_gradation(error_metric, positions, state, gamma)
    type(tensor_field), intent(inout) :: error_metric
    type(vector_field), intent(in) :: positions
    type(state_type), intent(in) :: state
    type(tensor_field), intent(in), optional :: gamma
    
    integer :: noits, grad_count

#ifdef HAVE_MPI
    include 'mpif.h'
    integer::noits_max, ierr
#endif
  
    if (use_anisotropic_gradation) then
      call form_anisotropic_gradation_metric(error_metric, positions, state, gamma_field=gamma)
    else if (use_gradation_metric) then
      if (.not. isparallel()) then
        call form_gradation_metric(positions, error_metric)
      else
        noits = 2
        grad_count = 0
        do while ((noits.gt.1).and.(grad_count.le.3))
           call form_gradation_metric(positions, error_metric, noits)
           call halo_update(error_metric)
           grad_count = grad_count + 1
#ifdef HAVE_MPI
           CALL MPI_Allreduce (noits, noits_max, 1, getpinteger(), MPI_MAX, MPI_COMM_FEMTOOLS, ierr)
           assert(ierr == MPI_SUCCESS)
           noits = noits_max
#endif
        end do
      end if
    end if
  
  end subroutine apply_gradation

  subroutine apply_horizontal_gradation(state, horizontal_metric, full_metric, horizontal_positions)
    type(state_type), intent(in) :: state
    type(tensor_field), intent(inout) :: horizontal_metric
    type(tensor_field), intent(in) :: full_metric
    type(vector_field), intent(in) :: horizontal_positions

    logical :: is_constant
    type(tensor_field) :: full_gamma, horizontal_gamma
    character(len=*), parameter :: path = "/mesh_adaptivity/hr_adaptivity/anisotropic_gradation/tensor_field::Gamma"
    type(vector_field) :: full_positions

    type(tensor_field), pointer :: min_bound, max_bound
    type(tensor_field) :: horizontal_min_bound, horizontal_max_bound

    if(.not.(use_anisotropic_gradation.or.use_gradation_metric)) return
    
    ewrite(1,*) 'in apply_horizontal_gradation'

    if (use_anisotropic_gradation) then
      is_constant = have_option(path//"/anisotropic_symmetric/constant")
      if (is_constant) then
        call allocate(full_gamma, full_metric%mesh, "Gamma", field_type=FIELD_TYPE_CONSTANT)
      else
        call allocate(full_gamma, full_metric%mesh, "Gamma")
      end if

      full_positions = get_coordinate_field(state, full_metric%mesh)
      call initialise_field(full_gamma, path, full_positions)

      call project_metric_to_surface(full_gamma, horizontal_positions, horizontal_gamma)

      call apply_gradation(horizontal_metric, horizontal_positions, state, gamma=horizontal_gamma)

      call deallocate(full_gamma)
      call deallocate(horizontal_gamma)
      call deallocate(full_positions)

    else
      call apply_gradation(horizontal_metric, horizontal_positions, state)
    end if

    if(use_anisotropic_gradation) then
      ! bounding the metric won't have happened yet because the min and max
      ! eigenbounds are the wrong size... so do it now!
      ! (min and max refer to edge lengths, not eigenvalues)
      min_bound => extract_tensor_field(state, "MaxMetricEigenbound")
      max_bound => extract_tensor_field(state, "MinMetricEigenbound")

      call project_metric_to_surface(min_bound, horizontal_positions, horizontal_min_bound)
      call project_metric_to_surface(max_bound, horizontal_positions, horizontal_max_bound)

      call bound_metric(horizontal_metric, horizontal_min_bound, horizontal_max_bound)

      call deallocate(horizontal_min_bound)
      call deallocate(horizontal_max_bound)
    end if

  end subroutine apply_horizontal_gradation

  subroutine apply_vertical_gradation(state, full_metric, full_positions, horizontal_positions)
    type(state_type), intent(in) :: state
    type(vector_field), intent(in) :: full_positions
    type(tensor_field), intent(inout) :: full_metric
    type(vector_field), intent(in) :: horizontal_positions
    
    integer :: column, node, i
    type(csr_sparsity) :: back_columns
    type(element_type) :: oned_shape
    type(quadrature_type) :: oned_quad
    integer :: quadrature_degree
    integer, parameter :: loc=2
    type(vector_field) :: oned_positions
    type(tensor_field) :: oned_metric

    logical :: is_constant
    type(tensor_field) :: full_gamma, oned_gamma
    character(len=*), parameter :: gamma_path = "/mesh_adaptivity/hr_adaptivity/anisotropic_gradation/tensor_field::Gamma"
    character(len=*), parameter :: path = "/mesh_adaptivity/hr_adaptivity/"

    type(tensor_field) :: min_edge, max_edge
    type(tensor_field) :: min_bound, max_bound
    type(tensor_field) :: oned_min_bound, oned_max_bound
    
    if(.not.(use_anisotropic_gradation.or.use_gradation_metric)) return

    ewrite(1,*) 'in apply_vertical_gradation'
    
    call get_option("/geometry/quadrature/degree", quadrature_degree)
    oned_quad = make_quadrature(vertices=loc, dim=1, degree=quadrature_degree)
    oned_shape = make_element_shape(vertices=loc, dim=1, degree=1, quad=oned_quad)
    call deallocate(oned_quad)

    call create_columns_sparsity(back_columns, full_positions%mesh)

    if (use_anisotropic_gradation) then
      ! unfortunately there are a few things that we need that would normally be in state
      ! but aren't at this stage... so allocate them and initialise them ourselves now
      is_constant = have_option(gamma_path//"/anisotropic_symmetric/constant")
      if (is_constant) then
        call allocate(full_gamma, full_metric%mesh, "Gamma", field_type=FIELD_TYPE_CONSTANT)
      else
        call allocate(full_gamma, full_metric%mesh, "Gamma")
      end if

      call initialise_field(full_gamma, gamma_path, full_positions)
      
      is_constant = (have_option(path // "/tensor_field::MinimumEdgeLengths/anisotropic_symmetric/constant"))
      if (is_constant) then
        call allocate(min_edge, full_metric%mesh, "MinimumEdgeLengths", field_type=FIELD_TYPE_CONSTANT)
        call initialise_field(min_edge, path // "/tensor_field::MinimumEdgeLengths", full_positions)
        call allocate(min_bound, full_metric%mesh, "MaxMetricEigenbound", field_type=FIELD_TYPE_CONSTANT)
        call set(min_bound, eigenvalue_from_edge_length(node_val(min_edge, 1)))
      else
        call allocate(min_edge, full_metric%mesh, "MinimumEdgeLengths")
        call initialise_field(min_edge, path // "/tensor_field::MinimumEdgeLengths", full_positions)
        call allocate(min_bound, full_metric%mesh, "MaxMetricEigenbound")
        do node=1,node_count(full_metric%mesh)
          call set(min_bound, node, eigenvalue_from_edge_length(node_val(min_edge, node)))
        end do
      end if
        
      call deallocate(min_edge)

      is_constant = (have_option(path // "/tensor_field::MaximumEdgeLengths/anisotropic_symmetric/constant"))
      if (is_constant) then
        call allocate(max_edge, full_metric%mesh, "MaximumEdgeLengths", field_type=FIELD_TYPE_CONSTANT)
        call initialise_field(max_edge, path // "/tensor_field::MaximumEdgeLengths", full_positions)
          call allocate(max_bound, full_metric%mesh, "MinMetricEigenbound", field_type=FIELD_TYPE_CONSTANT)
        call set(max_bound, eigenvalue_from_edge_length(node_val(max_edge, 1)))
      else
        call allocate(max_edge, full_metric%mesh, "MaximumEdgeLengths")
        call initialise_field(max_edge, path // "/tensor_field::MaximumEdgeLengths", full_positions)
        call allocate(max_bound, full_metric%mesh, "MinMetricEigenbound")
        do node=1,node_count(full_metric%mesh)
          call set(max_bound, node, eigenvalue_from_edge_length(node_val(max_edge, node)))
        end do
      end if
      
      if (mesh_periodic(full_metric)) then
        do i=1, mesh_dim(full_metric)
          if (minval(max_edge%val(i,i,:))<0.33*(domain_bbox(i,2)-domain_bbox(i,1))) then
            ewrite(0,*) "WARNING: Your MaximumEdgeLengths size is bigger than a third of the domain size."
            ewrite(0,*) "With periodic adaptivity this is probably not safe."
          end if
        end do
      end if
        
      call deallocate(max_edge)
      
    end if

    do column=1,node_count(horizontal_positions)
      if(use_anisotropic_gradation) then
        call get_1d_mesh(column, full_positions, back_columns, full_metric, oned_shape, oned_positions)
        
        call allocate(oned_metric, oned_positions%mesh, "1DMetric")
        call allocate(oned_gamma, oned_positions%mesh, "1DGamma", field_type=full_gamma%field_type)
        call get_1d_tensor(column, full_metric, oned_metric, back_columns)
        call get_1d_tensor(column, full_gamma, oned_gamma, back_columns)
        
        call apply_gradation(oned_metric, oned_positions, state, gamma=oned_gamma)
        
        call deallocate(oned_gamma)
        
        ! bounding the metric won't have happened yet because the min and max
        ! eigenbounds are the wrong size (and they weren't in state)... so do it now!
        ! (min and max refer to edge lengths, not eigenvalues)
        call allocate(oned_min_bound, oned_positions%mesh, "1DMaxMetricEigenbound", field_type=min_bound%field_type)
        call allocate(oned_max_bound, oned_positions%mesh, "1DMinMetricEigenbound", field_type=max_bound%field_type)
        call get_1d_tensor(column, min_bound, oned_min_bound, back_columns)
        call get_1d_tensor(column, max_bound, oned_max_bound, back_columns)

        call bound_metric(oned_metric, oned_min_bound, oned_max_bound)

        call deallocate(oned_min_bound)
        call deallocate(oned_max_bound)
      else
        call get_1d_mesh(column, full_positions, back_columns, full_metric, oned_shape, oned_positions)
        
        call allocate(oned_metric, oned_positions%mesh, "1DMetric")
        
        call get_1d_tensor(column, full_metric, oned_metric, back_columns)
        
        call apply_gradation(oned_metric, oned_positions, state)
      end if

      call recombine_metric(full_metric, column, oned_metric, back_columns)

      call deallocate(oned_positions)
      call deallocate(oned_metric)
    end do

    if (use_anisotropic_gradation) then
      call deallocate(full_gamma)
      call deallocate(min_bound)
      call deallocate(max_bound)
    end if
    call deallocate(back_columns)
    call deallocate(oned_shape)

  end subroutine apply_vertical_gradation

end module metric_assemble
