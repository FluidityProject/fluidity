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

#include "fdebug.h"

module gradient_matrix_cg

  use global_parameters, only: OPTION_PATH_LEN
  use fldebug
  use quadrature
  use futils
  use spud
  use transform_elements
  use fetools
  use fields
  use state_module
  use boundary_conditions
  use field_derivatives
  use field_options, only: complete_field_path

  implicit none

  private
  public :: assemble_gradient_matrix_cg

contains

    subroutine assemble_gradient_matrix_cg(C_m, state, c_rhs, & 
                                           test_mesh, field, option_path, &
                                           grad_mass, div_mass)

      ! inputs/outputs
      ! bucket full of fields
      type(state_type), intent(inout) :: state

      ! the pressure gradient and compressible gradient matrices
      type(block_csr_matrix), intent(inout) :: C_m

      type(vector_field), intent(inout), optional :: c_rhs

      type(mesh_type), intent(in) :: test_mesh

      type(scalar_field), intent(inout) :: field

      character(len=*), intent(in), optional :: option_path

      type(csr_matrix), intent(inout), optional :: grad_mass, div_mass

      ! local

      integer, dimension(:), pointer :: test_nodes, field_nodes
      integer, dimension(:), allocatable :: test_nodes_bdy, field_nodes_bdy

      real, dimension(:,:,:), allocatable :: ele_mat, ele_mat_bdy
      type(element_type), pointer :: field_shape, test_shape
      real, dimension(:,:,:), allocatable :: dfield_t
      real, dimension(:,:,:), allocatable :: dtest_t
      real, dimension(:), allocatable :: detwei, detwei_bdy
      real, dimension(:,:), allocatable :: normal_bdy

      ! loop integers
      integer :: ele, sele, dim

      ! pointer to coordinates
      type(vector_field), pointer :: x
      character(len=OPTION_PATH_LEN) :: l_option_path

      ! integrate by parts
      logical :: integrate_by_parts

      integer, dimension(:), allocatable :: field_bc_type
      type(scalar_field) :: field_bc

      real, dimension(:,:), allocatable :: div_mass_mat, grad_mass_mat
      
      integer :: stat

      ! =============================================================
      ! Subroutine to construct the matrix CT_m (a.k.a. C1/2/3T).
      ! =============================================================

      ewrite(2,*) 'In assemble_divergence_matrix_cg'

      if(present(option_path)) then
        l_option_path = trim(option_path)
      else
        l_option_path = trim(field%option_path)
      end if

      x=>extract_vector_field(state, "Coordinate")

      integrate_by_parts=have_option(trim(complete_field_path(l_option_path, stat=stat))//&
          &"/integrate_gradient_by_parts")

      ! Clear memory of arrays being designed
      call zero(C_m)
      if(present(c_rhs)) call zero(c_rhs)
      if(present(div_mass)) call zero(div_mass)
      if(present(grad_mass)) call zero(grad_mass)

      allocate(dfield_t(ele_loc(field, 1), ele_ngi(field, 1), mesh_dim(field)), &
               dtest_t(ele_loc(test_mesh, 1), ele_ngi(test_mesh, 1), mesh_dim(field)), &
               ele_mat(mesh_dim(test_mesh), ele_loc(test_mesh, 1), ele_loc(field, 1)), &
               detwei(ele_ngi(field, 1)), &
               grad_mass_mat(ele_loc(test_mesh, 1), ele_loc(test_mesh, 1)), &
               div_mass_mat(ele_loc(field, 1), ele_loc(field, 1)))

      do ele=1, element_count(test_mesh)
        test_nodes=>ele_nodes(test_mesh, ele)
        field_nodes=>ele_nodes(field, ele)

        test_shape=>ele_shape(test_mesh, ele)
        field_shape=>ele_shape(field, ele)

        if(integrate_by_parts) then
          ! transform the pressure derivatives into physical space
          ! (and get detwei)
          call transform_to_physical(X, ele, test_shape, dshape=dtest_t,&
               & detwei=detwei)

          ele_mat = -dshape_shape(dtest_t, field_shape, detwei)
        else
          ! transform the velociy derivatives into physical space
          ! (and get detwei)
          call transform_to_physical(X, ele, field_shape, dshape=dfield_t,&
               & detwei=detwei) 

          ele_mat = shape_dshape(test_shape, dfield_t, detwei)
        end if

        do dim = 1, mesh_dim(test_mesh)
          call addto(c_m, dim, 1, test_nodes, field_nodes, ele_mat(dim,:,:))
        end do

        if(present(div_mass)) then

          div_mass_mat = shape_shape(field_shape, field_shape, detwei)
          call addto(div_mass, field_nodes, field_nodes, div_mass_mat)

        end if

        if(present(grad_mass)) then

          grad_mass_mat = shape_shape(test_shape, test_shape, detwei)
          call addto(grad_mass, test_nodes, test_nodes, grad_mass_mat)

        end if

      end do

      if(integrate_by_parts) then

        allocate(detwei_bdy(face_ngi(field, 1)), &
                normal_bdy(mesh_dim(field), face_ngi(field, 1)))
        allocate(field_nodes_bdy(field%mesh%faces%shape%loc))
        allocate(test_nodes_bdy(test_mesh%faces%shape%loc))
        allocate(ele_mat_bdy(mesh_dim(field), face_loc(test_mesh, 1), face_loc(field, 1)))

        assert(surface_element_count(test_mesh)==surface_element_count(field))
        allocate(field_bc_type(surface_element_count(field)))
        call get_entire_boundary_condition(field, (/"weakdirichlet"/), field_bc, field_bc_type)

        do sele = 1, surface_element_count(test_mesh)

          test_shape=>face_shape(test_mesh, sele)
          field_shape=>face_shape(field, sele)

          test_nodes_bdy=face_global_nodes(test_mesh, sele)
          field_nodes_bdy=face_global_nodes(field, sele)

          call transform_facet_to_physical(X, sele, &
              &                          detwei_f=detwei_bdy,&
              &                          normal=normal_bdy) 

          ele_mat_bdy = shape_shape_vector(test_shape, field_shape, detwei_bdy, normal_bdy)

          do dim = 1, mesh_dim(field)
            if((field_bc_type(sele)==1).and.present(c_rhs)) then
              call addto(c_rhs, dim, test_nodes_bdy, &
                          -matmul(ele_mat_bdy(dim,:,:), &
                          ele_val(field_bc, sele)))
            else
              call addto(c_m, dim, 1, test_nodes_bdy, field_nodes_bdy, &
                          ele_mat_bdy(dim,:,:))
            end if
          end do

        end do

        call deallocate(field_bc)
        deallocate(field_bc_type)
        deallocate(detwei_bdy, normal_bdy)
        deallocate(test_nodes_bdy, field_nodes_bdy)

      end if

      deallocate(detwei)

    end subroutine assemble_gradient_matrix_cg

end module gradient_matrix_cg

