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

module shallow_water_diagnostics
  use diagnostic_source_fields
  use field_options
  use fields_manipulation
  use initialise_fields_module
  use fields
  use fldebug
  use global_parameters, only : timestep, OPTION_PATH_LEN
  use spud
  use state_fields_module
  use state_module
  use vtk_cache_module

  implicit none

  private

  public :: calculate_discontinuity_detector, calculate_manifold_divergence
  
contains

  ! Compute the discontinuity detector for 
  ! a scalar CG field (only works for P2+bubble currently) 
  subroutine calculate_discontinuity_detector(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    type(scalar_field), pointer :: source_field

    type(vector_field), pointer :: X
    !! Reconstructed scalar field
    type(scalar_field) :: rs_field
    real, allocatable, dimension(:) :: transfer
    real, allocatable, dimension(:) :: lumped_mass
    integer :: ele

    X => extract_vector_field(state, "Coordinate")
    source_field => scalar_source_field(state, s_field)
    call allocate(rs_field,source_field%mesh,"ReconstructedField")
    call zero(rs_field)
    
    !! Select transfer operator
    if((source_field%mesh%shape%numbering%type==ELEMENT_BUBBLE).and. &
         & source_field%mesh%shape%ndof==7) then
       allocate(transfer(7))
       transfer = (/ 1./20., 2./15., 1./20., &
            & 2./15., 2./15., 1./20., 9./20. /)/2
    else
       FLAbort('Unsupported element type.')
    end if

    !! Project to P0, and then project back to higher-order field
    !! using diagonal quadrature
    !! If it was a nodal basis 
    !! <gamma, s> = <gamma, S>
    !! LHS = \sum_e \sum_i w_i s_(e,i)
    !! RHS = \sum_e \sum_i w_i S_e

    !! Compute element means
    allocate(lumped_mass(node_count(rs_field)))

    lumped_mass = 0.

    do ele = 1, ele_count(rs_field)
       call assemble_reconstruction_ele(ele)
    end do    

    rs_field%val = rs_field%val/lumped_mass

    !! Need to take care of fact that not a nodal basis
    !! In each element, have computed value at s_(e,7)
    !! need to subtract off contributions from P2 component
    if(rs_field%mesh%shape%numbering%type==ELEMENT_BUBBLE) then
       do ele = 1, ele_count(rs_field)
          call fix_bubble_component(ele)
       end do
    end if

    !! Compute the discontinuity indicator
    !! Assumes a DG field for s_field
    call zero(s_field)
    do ele = 1, ele_count(s_field)
       call assemble_discontinuity_detector_ele(ele)
    end do

    call deallocate(rs_field)

  contains 
    
    subroutine assemble_discontinuity_detector_ele(ele)
      integer, intent(in) :: ele
      !
      real, dimension(ele_ngi(source_field,ele)) :: source_q, rs_q
      real, dimension(ele_loc(s_field,ele)) :: detector_vals

      source_q = ele_val_at_quad(source_field,ele)
      rs_q = ele_val_at_quad(rs_field,ele)

      detector_vals = maxval(abs(source_q-rs_q))
      call set(s_field,ele_nodes(s_field,ele),detector_vals)

    end subroutine assemble_discontinuity_detector_ele

    subroutine assemble_reconstruction_ele(ele)
      integer, intent(in) :: ele
      !
      real, dimension(ele_ngi(rs_field,ele)) :: s_gi, detwei
      real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
      real :: area, element_mean
      integer, dimension(:), pointer :: s_nodes

      !Uses compute_jacobian to be compatible with manifold stuff
      call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), &
           detwei=detwei,J=J)
      
      s_gi = ele_val_at_quad(source_field,ele)

      element_mean = sum(s_gi*detwei)
      area = sum(detwei)      

      s_nodes => ele_nodes(rs_field,ele)

      lumped_mass(s_nodes) = lumped_mass(s_nodes) + &
           & transfer*area

      call addto(rs_field,s_nodes,transfer*element_mean)
      
    end subroutine assemble_reconstruction_ele

    subroutine fix_bubble_component(ele)
      integer, intent(in) :: ele
      !
      integer :: nloc
      real, dimension(ele_loc(rs_field,ele)) :: N_vals
      integer, pointer, dimension(:) :: f_nodes
      real, dimension(ele_loc(rs_field,ele)) :: f_vals

      assert(rs_field%mesh%shape%numbering%type==ELEMENT_BUBBLE)      
      nloc = ele_loc(rs_field,ele)

      !Basis functions evaluated at bubble node.
      N_vals = eval_shape(ele_shape(rs_field,1), (/1.0/3.0,1.0/3.0,1.0/3.0/))

      f_nodes => ele_nodes(rs_field,ele)
      f_vals = ele_val(rs_field,ele)
      f_vals(nloc) = f_vals(nloc)-sum(f_vals(1:(nloc-1))*N_vals(1:(nloc-1)))
      f_vals(nloc) = f_vals(nloc)/N_vals(nloc)
      
      call set(rs_field,f_nodes(nloc),f_vals(nloc))

    end subroutine fix_bubble_component
    
  end subroutine calculate_discontinuity_detector

  subroutine calculate_manifold_divergence(state, s_field)
    implicit none
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field

    type(vector_field), pointer :: U, X
    type(scalar_field), pointer :: D
    integer :: ele

    U => extract_vector_field(state, "LocalVelocity")
    X => extract_vector_field(state, "Coordinate")
    D => extract_scalar_field(state, "LayerThickness")

    do ele = 1, ele_count(s_field)
       call assemble_manifold_divergence_ele(U, X, D, s_field, ele)
    end do


  end subroutine calculate_manifold_divergence

  subroutine assemble_manifold_divergence_ele(U, X, D, s_field, ele)
    implicit none
    type(vector_field), intent(in) :: U, X
    type(scalar_field), intent(in) :: D
    type(scalar_field), intent(inout) :: s_field
    integer, intent(in) :: ele

    real, dimension(mesh_dim(U),ele_loc(U,ele)) ::&
         & U_vals
    real, dimension(mesh_dim(U),ele_loc(U,ele),ele_loc(D,ele)) ::&
         & l_div_mat
      real, dimension(ele_loc(D,ele),ele_loc(D,ele)) :: &
           & d_mass_mat
    real, dimension(ele_loc(D,ele)) :: div_loc
    real, dimension(ele_ngi(s_field,ele)) :: detwei
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    type(element_type), pointer :: u_shape,d_shape
    integer :: dim1

    U_vals = ele_val(U,ele)
    u_shape => ele_shape(U,ele)
    d_shape => ele_shape(D,ele)

    l_div_mat = dshape_shape(u_shape%dn,d_shape,D_shape%quadrature%weight)

    div_loc = 0.
    do dim1 = 1, U%dim
       div_loc = div_loc + matmul(transpose(l_div_mat(dim1,:,:))&
            &,U_vals(dim1,:))
    end do

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), &
         detwei=detwei,J=J)
    d_mass_mat = shape_shape(d_shape,d_shape,detwei)

    call solve(d_mass_mat, div_loc)

    call set(s_field, ele_nodes(s_field,ele), div_loc)

  end subroutine assemble_manifold_divergence_ele

end module shallow_water_diagnostics
