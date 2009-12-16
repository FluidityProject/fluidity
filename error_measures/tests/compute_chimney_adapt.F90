
subroutine compute_chimney_adapt
#define TEMP_ERROR 0.2
#define TEMP_REL .true.
#define TEMP_MIN 0.001
#define TEMP_SQ  .false.

#define VEL_ERROR (/0.5, 0.5, 0.5/)
#define VEL_REL .false.
#define VEL_MIN (/0.001, 0.001, 0.001/)
#define VEL_SQ .false.
#define NADAPT 1

  use global_parameters, only: current_debug_level, pseudo2d_coord
  use metric_assemble
  use adapt_state_unittest_module, only : adapt_state => adapt_state_unittest
  use form_metric_field
  use state_module
  use vtk_interfaces
  use vector_tools
  use unittest_tools
  use edge_length_module
  use gradation_metric
  use mpi
  use interpolation_error
  use field_options
  use smoothing_module
  implicit none
  
  type(state_type) :: state, dummy(1)
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions
  type(scalar_field), pointer :: porosity, permeability
  type(scalar_field), pointer :: porosity_old, permeability_old
  type(tensor_field), pointer :: metric
  type(tensor_field) :: tmp_metric
  type(scalar_field) :: edgelen
  type(metric_options) :: opts
  type(anisotropic_edge_options), target :: edge_opts
  real, dimension(4), target :: xhsamp, yhsamp, zhsamp
  real, dimension(4), target :: hminxx, hminxy, hminxz, hminyy, hminyz, hminzz
  real, dimension(4), target :: hmaxxx, hmaxxy, hmaxxz, hmaxyy, hmaxyz, hmaxzz
  integer :: nhsamp
  integer :: i, stat
  real, dimension(3, 3)::alpha=(/0.001, 0.0, 0.0/, /0.0, 0.001, 0.0/,/0.0, 0.0, 0.001/)

  interface
    function solution(pos)
      real, dimension(:) :: pos
      real :: solution
    end function
  end interface
  interface
    function gradsoln(pos)
      real, dimension(:) :: pos
      real, dimension(size(pos)) :: gradsoln
    end function
  end interface

  call vtk_read_state("/home/gormo/reservoir.vtu", state)

  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  porosity => extract_scalar_field(state, "porosity")
  permeability => extract_scalar_field(state, "permeability")

  call adaptivity_options(state, porosity, TEMP_ERROR, TEMP_REL, TEMP_MIN)
  porosity_old = porosity
  call smooth_scalar(porosity_old, positions, porosity, alpha)

  call adaptivity_options(state, permeability, TEMP_ERROR, TEMP_REL, TEMP_MIN)
  permeability_old = permeability
  call smooth_scalar(permeability_old, positions, permeability, alpha)
  
  call allocate(tmp_metric, mesh, "Metric")
  call insert(state, tmp_metric, "Metric")
  call deallocate(tmp_metric)
  metric => extract_tensor_field(state, "Metric")

  opts%min_edge_length = 0.002
  opts%max_edge_length = 10.0
  opts%use_anisotropic_edge_length = .false.
  xhsamp(1) = 0.0; yhsamp(1) = 0.0; zhsamp(1) = 0.0
  xhsamp(2) = 1.0; yhsamp(2) = 0.0; zhsamp(2) = 1.0
  xhsamp(3) = 1.0; yhsamp(3) = 1.0; zhsamp(3) = 0.0
  xhsamp(4) = 0.0; yhsamp(4) = 1.0; zhsamp(4) = 1.0
  hminxx = 0.00001; hminxy = 0.0; hminxz = 0.0; hminyy = 0.00001; hminyz = 0.0; hminzz = 0.00001
  hmaxxx = 0.5; hmaxxy = 0.0; hmaxxz = 0.0; hmaxyy = 0.5; hmaxyz = 0.0; hmaxzz = 0.4
  nhsamp = 4
  edge_opts%no_samp = NHSAMP
  edge_opts%x => XHSAMP(1:NHSAMP)
  edge_opts%y => YHSAMP(1:NHSAMP)
  edge_opts%z => ZHSAMP(1:NHSAMP)
  edge_opts%hminxx => HMINXX(1:NHSAMP); edge_opts%hmaxxx => HMAXXX(1:NHSAMP)
  edge_opts%hminxy => HMINXY(1:NHSAMP); edge_opts%hmaxxy => HMAXXY(1:NHSAMP)
  edge_opts%hminxz => HMINXZ(1:NHSAMP); edge_opts%hmaxxz => HMAXXZ(1:NHSAMP)
  edge_opts%hminyy => HMINYY(1:NHSAMP); edge_opts%hmaxyy => HMAXYY(1:NHSAMP)
  edge_opts%hminyz => HMINYZ(1:NHSAMP); edge_opts%hmaxyz => HMAXYZ(1:NHSAMP)
  edge_opts%hminzz => HMINZZ(1:NHSAMP); edge_opts%hmaxzz => HMAXZZ(1:NHSAMP)
  opts%anisotropic_edge_opts => edge_opts

  call set_option("/mesh_adaptivity/hr_adaptivity/geometric_constraints", 1.3, stat=stat)

  dummy(1) = state
  call assemble_metric(dummy, metric, opts)
  state = dummy(1)
  call allocate(edgelen, mesh, "Edge lengths")
  call get_edge_lengths(metric, edgelen)
  call insert(state, edgelen, "Edge lengths")
  call deallocate(edgelen)
  call vtk_write_state("data/chimney_adapt", 0, state=(/state/))
  call adapt_state(state, metric)

  do i=1,NADAPT-1
    mesh => extract_mesh(state, "Mesh")
    positions => extract_vector_field(state, "Coordinate")
    porosity => extract_scalar_field(state, "porosity")
    permeability => extract_scalar_field(state, "permeability")

    call adaptivity_options(state, porosity, TEMP_ERROR, TEMP_REL, TEMP_MIN)
    call adaptivity_options(state, permeability, TEMP_ERROR, TEMP_REL, TEMP_MIN)
    
    call deallocate(metric); call allocate(metric, mesh, "Metric")
    dummy(1) = state
    call assemble_metric(dummy, metric, opts)
    state = dummy(1)
    call vtk_write_state("/tmp/reservoir_adapt", i, state=(/state/))
    call adapt_state(state, metric) 
  end do

  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")

  porosity => extract_scalar_field(state, "porosity")
  permeability => extract_scalar_field(state, "permeability")

  call vtk_write_state("/tmp/reservoir_adapt", NADAPT, state=(/state/))
end subroutine compute_chimney_adapt
