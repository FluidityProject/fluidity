module simple_advection
implicit none

contains

subroutine advection_action(x, u, c, Ac)
  real, dimension(:), intent(in) :: x
  real, dimension(:), intent(in) :: u
  real, dimension(:), intent(in) :: c
  real, dimension(:), intent(out) :: Ac

  integer :: ele, ele_count, node_count
  integer, dimension(2) :: ele_nodes
  real, dimension(2) :: ele_tmp
  node_count = size(x)
  if (size(c) /= node_count .or. size(u) /= node_count) then
    write(0,*) "Huh? Everything has to be consistent"
    stop
  end if

  ele_count = node_count-1 ! 1D only, baby
  Ac = 0.0

  do ele=1,ele_count
    ele_nodes = (/ele, ele+1/)
    call ele_advection_action(ele, ele_nodes, x, u, c, Ac)
  end do
end subroutine advection_action

subroutine ele_advection_action(ele, ele_nodes, x, u, c, Ac)
  integer, intent(in) :: ele
  integer, dimension(2), intent(in) :: ele_nodes
  real, dimension(:), intent(in) :: x, u, c
  real, dimension(:), intent(inout) :: Ac

  real, dimension(2, 2) :: A

               ! loc x ngi
  real, dimension(2, 2) :: shape_n
               ! log x ngi x dim
  real, dimension(2, 2, 1) :: dshape_n
  real :: h
  real, dimension(2) :: detwei

  integer :: i, j
  real, dimension(2) :: u_at_quad

  shape_n(1, :) = (/0.78867513459481298, 0.21132486540518702/) ! values of basis functions at quad points
  shape_n(2, :) = (/0.21132486540518702, 0.78867513459481298/)
  dshape_n(1, :, 1) = (/-1, -1/) ! values of derivatives of basis functions at quad points
  dshape_n(2, :, 1) = (/+1, +1/)

  h = x(ele_nodes(2)) - x(ele_nodes(1)) ! step size
  dshape_n = 1.0/h * dshape_n ! transform_to_physical
  detwei = (/0.5, 0.5/) * h

  ! Replacement matmul
  do i=1,2
    u_at_quad(i) = sum(u(ele_nodes) * shape_n(:, i))
  end do

  do i=1,2
    do j=1,2
      ! Replacement dot_product
      A(i,j) = sum((shape_n(i,:) * (dshape_n(j, :, 1) * u_at_quad)) * detwei)
    end do
  end do

  ! Enforce dirichlet BCs
  if (ele_nodes(1) == 1) then
    A(1,:) = (/1.0, 0.0/)
  end if
  if (ele_nodes(2) == size(x)) then
    A(2,:) = (/0.0, 1.0/)
  end if

  ! Replacement matmul
  do i=1,2
    Ac(ele_nodes(i)) = Ac(ele_nodes(i)) + sum(A(i,:) * c(ele_nodes))
  end do
end subroutine ele_advection_action

end module simple_advection
