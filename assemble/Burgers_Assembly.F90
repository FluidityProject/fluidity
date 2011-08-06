  !    Copyright (C) 2006 Imperial College London and others.
  !
  !    Please see the AUTHORS file in the main source directory for a full list
  !    of copyright holders.
  !
  !    Prof. C Pain
  !    Applied Modelling and Computation Group
  !    Department of Earth Science and Engineering
  !    Imperial College London
  !
  !    amcgsoftware@imperial.ac.uk
  !
  !    This library is free software; you can redistribute it and/or
  !    modify it under the terms of the GNU Lesser General Public
  !    License as published by the Free Software Foundation,
  !    version 2.1 of the License.
  !
  !    This library is distributed in the hope that it will be useful,
  !    but WITHOUT ANY WARRANTY; without even the implied warranty of
  !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  !    Lesser General Public License for more details.
  !
  !    You should have received a copy of the GNU Lesser General Public
  !    License along with this library; if not, write to the Free Software
  !    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
  !    USA
module burgers_assembly
  use fields
  use sparse_tools
  use spud
  use transform_elements
  use fetools
  implicit none

  private
  public :: assemble_advection_matrix


  contains

  subroutine assemble_advection_matrix(advection_matrix, x, u_left, u_right)
    type(csr_matrix), intent(inout) :: advection_matrix
    type(vector_field), intent(in) :: x
    type(scalar_field), intent(in), target :: u_left, u_right

    type(scalar_field) :: nu
    type(mesh_type), pointer :: mesh
    integer :: ele
    real :: itheta

    call get_option(trim(u_left%option_path) // "/prognostic/temporal_discretisation/relaxation", itheta, default=0.5)
    mesh => u_left%mesh
    call allocate(nu, mesh, "NonlinearVelocity")
    call set(nu, u_left)
    call scale(nu, (1.0 - itheta))
    call addto(nu, u_right, scale=itheta)
    nu%option_path = u_left%option_path

    call zero(advection_matrix)
    if (.not. have_option(trim(nu%option_path) // "/prognostic/remove_advection_term")) then
      do ele=1,ele_count(nu)
        call assemble_advection_matrix_ele(advection_matrix, x, nu, ele)
      end do
    end if

    call deallocate(nu)
  end subroutine assemble_advection_matrix

  subroutine assemble_advection_matrix_ele(advection_matrix, x, nu, ele)
    type(csr_matrix), intent(inout) :: advection_matrix
    type(vector_field), intent(in) :: x
    type(scalar_field), intent(in) :: nu
    integer, intent(in) :: ele

    real, dimension(ele_loc(nu, ele), ele_loc(nu, ele)) :: little_advection_matrix
    real, dimension(ele_ngi(nu, ele)) :: detwei
    real, dimension(ele_loc(nu, ele), ele_ngi(nu, ele), x%dim) :: du_t

    call transform_to_physical(x, ele, ele_shape(nu, ele), detwei=detwei, dshape=du_t)
    little_advection_matrix = shape_vector_dot_dshape(ele_shape(nu, ele), ele_val_at_quad(nu, ele), du_t, detwei)
    call addto(advection_matrix, ele_nodes(nu, ele), ele_nodes(nu, ele), little_advection_matrix)
  end subroutine assemble_advection_matrix_ele

end module burgers_assembly
