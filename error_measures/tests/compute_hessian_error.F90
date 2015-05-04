#define DIMENSION 2

subroutine compute_hessian_error

  use mesh_files
  use field_derivatives
  use fields
  use matrix_norms
  use vtk_interfaces
  use global_parameters
  use state_module
  implicit none

!  type(quadrature_type) :: quad
!  type(element_type) :: x_shape
  type(vector_field), pointer :: positions
  type(scalar_field) :: field
  type(tensor_field) :: hessian, ex_hessian
  type(state_type) :: state
  integer :: i, node

  real :: h

  interface
    function exact_field(pos) result(val)
      real, dimension(:), intent(in) :: pos
      real :: val
    end function exact_field
  end interface

  interface
    function exact_hessian(pos) result(hess)
      real, dimension(:), intent(in) :: pos
      real, dimension(size(pos), size(pos)) :: hess
    end function exact_hessian
  end interface

!  quad = make_quadrature(vertices = DIMENSION+1, dim =DIMENSION, degree=8)
!  x_shape = make_element_shape(vertices = DIMENSION+1, dim =DIMENSION, degree=1, quad=quad)
!  positions = read_mesh_files("data/square.1", x_shape, format="gmsh")

  pseudo2d_coord = 3

  call vtk_read_state("data/stretchedjack-itv20.vtu", state)
  positions => extract_vector_field(state, "Coordinate")

  call allocate(field, positions%mesh, "Field")
  call allocate(hessian, positions%mesh, "Hessian")
  call allocate(ex_hessian, positions%mesh, "Exact Hessian")

  call set_from_function(field, exact_field, positions)
  call set_from_function(ex_hessian, exact_hessian, positions)

  call compute_hessian_qf(field, positions, hessian)
  if (pseudo2d_coord /= 0) then
    do node=1,node_count(hessian)
      hessian%val(pseudo2d_coord, :, node) = 0.0
      hessian%val(:, pseudo2d_coord, node) = 0.0
    end do
  end if

  do i=1,node_count(positions)
    hessian%val(:, :, i) = hessian%val(:, :, i) - ex_hessian%val(:, :, i)
  end do

  h = sqrt(1.0 / ele_count(positions)) ! Area / elements

  write(0,*) "h: ", h
  write(0,*) "-------------- QF ---------------"
  write(0,*) "one: ", one_norm(hessian)
  write(0,*) "two: ", two_norm(hessian)
  write(0,*) "inf: ", inf_norm(hessian)

  call compute_hessian_eqf(field, positions, hessian)
  if (pseudo2d_coord /= 0) then
    do node=1,node_count(hessian)
      hessian%val(pseudo2d_coord, :, node) = 0.0
      hessian%val(:, pseudo2d_coord, node) = 0.0
    end do
  end if

  do i=1,node_count(positions)
    hessian%val(:, :, i) = hessian%val(:, :, i) - ex_hessian%val(:, :, i)
  end do

  write(0,*) "-------------- EQF ---------------"
  write(0,*) "one: ", one_norm(hessian)
  write(0,*) "two: ", two_norm(hessian)
  write(0,*) "inf: ", inf_norm(hessian)

  call compute_hessian_var(field, positions, hessian)
  if (pseudo2d_coord /= 0) then
    do node=1,node_count(hessian)
      hessian%val(pseudo2d_coord, :, node) = 0.0
      hessian%val(:, pseudo2d_coord, node) = 0.0
    end do
  end if

  do i=1,node_count(positions)
    hessian%val(:, :, i) = hessian%val(:, :, i) - ex_hessian%val(:, :, i)
  end do

  write(0,*) "-------------- VAR ---------------"
  write(0,*) "one: ", one_norm(hessian)
  write(0,*) "two: ", two_norm(hessian)
  write(0,*) "inf: ", inf_norm(hessian)
end subroutine compute_hessian_error

!function exact_field(pos) result(val)
!  real, dimension(:), intent(in) :: pos
!  real :: val
!  real :: x, y
!
!  x = pos(1); y = pos(2)
!
!  val = exp(-25 * x) + exp(-25 * y)
!end function exact_field
!
!function exact_hessian(pos) result(hess)
!  real, dimension(:), intent(in) :: pos
!  real, dimension(size(pos), size(pos)) :: hess
!  real :: x, y
!  x = pos(1); y = pos(2)
!
!  hess = 0.0
!  hess(1, 1) = 625*exp(-(25*x))
!  hess(1, 2) = 0
!  hess(2, 1) = hess(1, 2)
!  hess(2, 2) = 625*exp(-(25*x))
!end function exact_hessian

function exact_field(pos) result(val)
  real, dimension(:), intent(in) :: pos
  real :: val
  real :: x, y, z

  x = pos(1); y = pos(2); z = pos(3)
  val = exp(-100 * ((x - 0.5)**2 + (y - 0.5)**2))
end function exact_field

function exact_hessian(pos) result(hess)
  real, dimension(:), intent(in) :: pos
  real, dimension(size(pos), size(pos)) :: hess
  real :: x, y, z
  x = pos(1); y = pos(2); z = pos(3)

  hess = 0.0
  hess(1, 1) = 40000*(x - 0.5)**2*exp(-100*((y - 0.5)**2 + (x - 0.5)**2)) - 200*exp(-100*((y - 0.5)**2 + (x - 0.5)**2))
  hess(1, 2) = 40000*(x - 0.5)*(y - 0.5)*exp(-100*((y - 0.5)**2 + (x - 0.5)**2))
  hess(2, 1) = hess(1, 2)
  hess(2, 2) = 40000*(y - 0.5)**2*exp(-100*((y - 0.5)**2 + (x - 0.5)**2)) - 200*exp(-100*((y - 0.5)**2 + (x - 0.5)**2))
end function exact_hessian
