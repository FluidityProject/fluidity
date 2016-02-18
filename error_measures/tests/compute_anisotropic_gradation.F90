subroutine compute_anisotropic_gradation

  use mesh_files
  use anisotropic_gradation
  use mba_adapt_module
  use state_module
  use vtk_interfaces

  type(vector_field) :: positions
  type(tensor_field) :: metric, gamma
  real, dimension(2, 2) :: id, nid
  type(state_type) :: state

  interface
    function set_gamma(pos)
      real, dimension(:) :: pos
      real, dimension(size(pos), size(pos)) :: set_gamma
    end function
  end interface

  positions = read_mesh_files("data/laplacian_grid.2", quad_degree=4, format="gmsh")
  call insert(state, positions, "Coordinate")
  call insert(state, positions%mesh, "Mesh")
  call allocate(metric, positions%mesh, "Metric")

  id = reshape((/1.0, 0.0, 0.0, 1.0/), (/2, 2/))
  nid = reshape((/0.1, 0.0, 0.0, 5.0/), (/2, 2/))

  call set(metric, id) ! an isotropic edge length of 0.01
  call set(metric, 1, id * 1000000)  ! an isotropic edge length of 0.001

  call allocate(gamma, positions%mesh, "Gamma", FIELD_TYPE_NORMAL)
  call set_from_function(gamma, set_gamma, positions) 

  call form_anisotropic_gradation_metric(metric, positions, gamma_field=gamma)
  call vtk_write_state("data/anisotropic_gradation", 0, state=(/state/))
  call mba_adapt(state, metric)
  call vtk_write_state("data/anisotropic_gradation", 1, state=(/state/))

end subroutine compute_anisotropic_gradation

function set_gamma(pos)
  real, dimension(:) :: pos
  real, dimension(size(pos), size(pos)) :: set_gamma
  real :: x, y

  x = pos(1); y = pos(2)

  set_gamma = 0.0
  set_gamma(1, 1) = max(x, 0.01)
  !set_gamma(2, 2) = max(1 - x, 0.01)
  set_gamma(2, 2) = 1.0
end function
