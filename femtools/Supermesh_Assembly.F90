!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineeringp
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

module supermesh_assembly

! these 5 need to be on top and in this order, so as not to confuse silly old intel compiler 
  use quadrature
  use elements
  use sparse_tools
  use fields
  use state_module
!
  use adaptive_interpolation_module
  use fldebug
  use field_options
  use interpolation_module
  use intersection_finder_module
  use linked_lists
  use state_fields_module
  use solvers
  use supermesh_construction
  use transform_elements
  
  implicit none
  
  private
  
  public :: project_donor_shape_to_supermesh, &
    & project_target_shape_to_supermesh, construct_supermesh_ele, &
    & extruded_shape_function, generate_supermesh_node_ownership, &
    & project_donor_field_to_supermesh, project_target_field_to_supermesh, &
    & galerkin_projection_scalars, compute_inner_product_sa
    
  interface generate_supermesh_local_coords
    module procedure generate_supermesh_local_coords_ele,  &
      & generate_supermesh_local_coords_eles
  end interface generate_supermesh_local_coords
    
  interface project_donor_shape_to_supermesh
    module procedure project_donor_shape_to_supermesh_mesh, &
      & project_donor_shape_to_supermesh_shape
  end interface project_donor_shape_to_supermesh
  
  interface project_target_shape_to_supermesh
    module procedure project_target_shape_to_supermesh_mesh, &
      & project_target_shape_to_supermesh_shape
  end interface project_target_shape_to_supermesh
  
  interface construct_supermesh_dn
    module procedure construct_supermesh_dn_ele, construct_supermesh_dn_eles, &
      & construct_supermesh_dn_ele_ele_c
  end interface construct_supermesh_dn
  
  interface project_donor_field_to_supermesh
    module procedure project_donor_field_to_supermesh_scalar
  end interface project_donor_field_to_supermesh
  
  interface project_target_field_to_supermesh
    module procedure project_target_field_to_supermesh_scalar
  end interface project_target_field_to_supermesh
    
  interface construct_supermesh_ele
    module procedure construct_supermesh_ele_single_state, &
      & construct_supermesh_ele_multiple_states
  end interface construct_supermesh_ele
  
contains

  subroutine generate_supermesh_local_coords_ele(ele, positions, positions_c, base_shape_c, &
    & l_coords)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: positions_c
    type(element_type), intent(in) :: base_shape_c
    
    real, dimension(ele_loc(positions, ele), ele_ngi(positions, ele), ele_count(positions_c)) :: l_coords
    
    integer :: ele_c
    type(mesh_type) :: positions_c_remap_mesh
    type(vector_field) :: positions_c_remap
    
    if(ele_shape(positions_c, ele) == base_shape_c) then
      positions_c_remap = positions_c
      call incref(positions_c_remap)
    else
      positions_c_remap_mesh = make_mesh(positions_c%mesh, base_shape_c, continuity = -1, name = "CoordinateRemapMesh")
      call allocate(positions_c_remap, positions_c%dim, positions_c_remap_mesh, name = "CoordinateRemap")
      call deallocate(positions_c_remap_mesh)
      call remap_field(positions_c, positions_c_remap)
    end if

    do ele_c = 1, size(l_coords, 3)
      l_coords(:, :, ele_c) = local_coords(positions, ele, ele_val_at_quad(positions_c_remap, ele_c))
    end do

    call deallocate(positions_c_remap)
  
  end subroutine generate_supermesh_local_coords_ele

  subroutine generate_supermesh_local_coords_eles(eles, positions, positions_c, base_shape_c, &
    & l_coords)
    type(vector_field), intent(in) :: positions_c
    integer, dimension(ele_count(positions_c)), intent(in) :: eles
    type(vector_field), intent(in) :: positions
    type(element_type), intent(in) :: base_shape_c
    
    real, dimension(ele_loc(positions, 1), ele_ngi(positions, 1), ele_count(positions_c)) :: l_coords
    
    integer :: ele_c
    type(mesh_type) :: positions_c_remap_mesh
    type(vector_field) :: positions_c_remap
    
    assert(ele_count(positions_c) > 0)
    if(ele_shape(positions_c, 1) == base_shape_c) then
      positions_c_remap = positions_c
      call incref(positions_c_remap)
    else
      positions_c_remap_mesh = make_mesh(positions_c%mesh, base_shape_c, continuity = -1, name = "CoordinateRemapMesh")
      call allocate(positions_c_remap, positions_c%dim, positions_c_remap_mesh, name = "CoordinateRemap")
      call deallocate(positions_c_remap_mesh)
      call remap_field(positions_c, positions_c_remap)
    end if

    do ele_c = 1, size(l_coords, 3)
      l_coords(:, :, ele_c) = local_coords(positions, eles(ele_c), ele_val_at_quad(positions_c_remap, ele_c))
    end do

    call deallocate(positions_c_remap)
  
  end subroutine generate_supermesh_local_coords_eles

  subroutine project_donor_shape_to_supermesh_mesh(positions_a, shape_mesh, positions_c, &
    & shapes_c, form_dn)
    type(vector_field), intent(in) :: positions_a
    type(mesh_type), intent(in) :: shape_mesh
    type(vector_field), intent(in) :: positions_c
    type(element_type), dimension(:), allocatable, intent(out) :: shapes_c
    ! If present and .false., do not form the shape function derivatives
    logical, optional, intent(in) :: form_dn
        
    assert(ele_count(shape_mesh) > 0)
    call project_donor_shape_to_supermesh(positions_a, ele_shape(shape_mesh, 1), positions_c, &
      & shapes_c, form_dn = form_dn)
    
  end subroutine project_donor_shape_to_supermesh_mesh

  subroutine project_donor_shape_to_supermesh_shape(positions_a, base_shape_c, positions_c, &
    & shapes_c, form_dn)
    type(vector_field), intent(in) :: positions_a
    type(element_type), target, intent(in) :: base_shape_c
    type(vector_field), intent(in) :: positions_c
    type(element_type), dimension(:), allocatable, intent(out) :: shapes_c
    ! If present and .false., do not form the shape function derivatives
    logical, optional, intent(in) :: form_dn
    
    integer :: dim, degree, coords, i, j, loc, ngi
    integer, dimension(:), pointer :: eles_a
    logical :: lform_dn
    real, dimension(ele_loc(positions_a, 1), ele_ngi(positions_a, 1), ele_count(positions_c)) :: l_coords
    type(quadrature_type), pointer :: quad
    type(ele_numbering_type), pointer :: ele_num
    
    lform_dn = .not. present_and_false(form_dn)
    
    eles_a => ele_region_ids(positions_c)
    
    quad => base_shape_c%quadrature
    
    dim = base_shape_c%dim
    loc = base_shape_c%loc
    ngi = quad%ngi
    coords = local_coord_count(base_shape_c)
    degree = base_shape_c%degree

    if(base_shape_c%degree > 0 .or. lform_dn) then
      call generate_supermesh_local_coords(eles_a, positions_a, positions_c, base_shape_c, &
        & l_coords)
    end if
    
    allocate(shapes_c(ele_count(positions_c)))
    do i = 1, size(shapes_c)
       ele_num => find_element_numbering(&
            &vertices = base_shape_c%numbering%vertices, &
            &dimension = dim, degree =&
            & degree)    
      call allocate(shapes_c(i), ele_num=ele_num, ngi = ngi)
      
      shapes_c(i)%degree = degree
      shapes_c(i)%numbering => find_element_numbering(vertices = loc, dimension = dim, degree = degree)
      shapes_c(i)%quadrature = quad
      call incref(quad)
      
      shapes_c(i)%dn = huge(0.0) 
      assert(.not. associated(shapes_c(i)%dn_s))
      assert(.not. associated(shapes_c(i)%n_s))
      deallocate(shapes_c(i)%spoly)
      nullify(shapes_c(i)%spoly)
      deallocate(shapes_c(i)%dspoly)
      nullify(shapes_c(i)%dspoly)
      
      select case(base_shape_c%degree)
        case(0)
          shapes_c(i)%n = 1.0
        case(1)
          if(ele_numbering_family(base_shape_c) == FAMILY_SIMPLEX) then
            shapes_c(i)%n = l_coords(:, :, i)   
          else
            do j = 1, ngi
              shapes_c(i)%n(:, j) = eval_shape(base_shape_c, l_coords(:, j, i))
            end do   
          end if
        case default
          do j = 1, ngi
            shapes_c(i)%n(:, j) = eval_shape(base_shape_c, l_coords(:, j, i))
          end do    
      end select
    end do

    if(lform_dn) then
      call construct_supermesh_dn(eles_a, positions_a, positions_c, l_coords, base_shape_c, shapes_c)
    end if
    
  end subroutine project_donor_shape_to_supermesh_shape

  subroutine project_target_shape_to_supermesh_mesh(ele_b, &
    & positions_b, shape_mesh, positions_c, &
    & shapes_c, form_dn)
    integer, intent(in) :: ele_b
    type(vector_field), intent(in) :: positions_b
    type(mesh_type), intent(in) :: shape_mesh
    type(vector_field), intent(in) :: positions_c
    type(element_type), dimension(:), allocatable, intent(out) :: shapes_c
    ! If present and .false., do not form the shape function derivatives
    logical, optional, intent(in) :: form_dn
        
    call project_target_shape_to_supermesh(ele_b, &
      & positions_b, ele_shape(shape_mesh, ele_b), positions_c, &
      & shapes_c, form_dn = form_dn)
    
  end subroutine project_target_shape_to_supermesh_mesh

  subroutine project_target_shape_to_supermesh_shape(ele_b, &
    & positions_b, base_shape_c, positions_c, &
    & shapes_c, form_dn)
    integer, intent(in) :: ele_b
    type(vector_field), intent(in) :: positions_b
    type(element_type), target, intent(in) :: base_shape_c
    type(vector_field), intent(in) :: positions_c
    type(element_type), dimension(:), allocatable, intent(out) :: shapes_c
    ! If present and .false., do not form the shape function derivatives
    logical, optional, intent(in) :: form_dn
    
    integer :: dim, degree, coords, i, j, loc, ngi
    logical :: lform_dn
    real, dimension(ele_loc(positions_b, ele_b), ele_ngi(positions_b, ele_b), ele_count(positions_c)) :: l_coords
    type(quadrature_type), pointer :: quad
    type(ele_numbering_type), pointer :: ele_num
    
    lform_dn = .not. present_and_false(form_dn)
    
    quad => base_shape_c%quadrature
    
    dim = base_shape_c%dim
    loc = base_shape_c%loc
    ngi = quad%ngi
    coords = local_coord_count(base_shape_c)
    degree = base_shape_c%degree
    
    if(base_shape_c%degree > 0 .or. lform_dn) then
      call generate_supermesh_local_coords(ele_b, positions_b, positions_c, base_shape_c, &
        & l_coords)
    end if
    
    allocate(shapes_c(ele_count(positions_c)))
    do i = 1, size(shapes_c)
       ele_num => find_element_numbering(&
            vertices = base_shape_c%numbering%vertices, dimension = dim, degree =&
            & degree)    
      call allocate(shapes_c(i), ele_num, ngi = ngi)
      
      shapes_c(i)%degree = degree
      shapes_c(i)%numbering => find_element_numbering(vertices = loc, dimension = dim, degree = degree)
      shapes_c(i)%quadrature = quad
      call incref(quad)
      
      shapes_c(i)%dn = huge(0.0) 
      assert(.not. associated(shapes_c(i)%dn_s))
      assert(.not. associated(shapes_c(i)%n_s))
      deallocate(shapes_c(i)%spoly)
      nullify(shapes_c(i)%spoly)
      deallocate(shapes_c(i)%dspoly)
      nullify(shapes_c(i)%dspoly)
      
      select case(base_shape_c%degree)
        case(0)
          shapes_c(i)%n = 1.0
        case(1)
          if(ele_numbering_family(base_shape_c) == FAMILY_SIMPLEX) then
            shapes_c(i)%n = l_coords(:, :, i)   
          else
            do j = 1, ngi
              shapes_c(i)%n(:, j) = eval_shape(base_shape_c, l_coords(:, j, i))
            end do   
          end if
        case default
          do j = 1, ngi
            shapes_c(i)%n(:, j) = eval_shape(base_shape_c, l_coords(:, j, i))
          end do    
      end select
    end do

    if(lform_dn) then
      call construct_supermesh_dn(ele_b, positions_b, positions_c, l_coords, base_shape_c, shapes_c)
    end if
    
  end subroutine project_target_shape_to_supermesh_shape
    
  subroutine construct_supermesh_dn_ele(ele, positions, positions_c, l_coords, base_shape, shapes_c)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: positions_c
    real, dimension(ele_loc(positions, ele), ele_ngi(positions, ele), ele_count(positions_c)), intent(in) :: l_coords
    type(element_type), intent(in) :: base_shape
    type(element_type), dimension(ele_count(positions_c)), intent(inout) :: shapes_c
    
    integer :: i, j, k
    real, dimension(positions%dim, positions%dim, ele_ngi(positions, ele)) :: invj
    real, dimension(positions_c%dim, positions_c%dim, ele_ngi(positions, ele)) :: j_c
    
    if(base_shape%degree == 0) then
      ! This case is nice and easy
      do i = 1, size(shapes_c)
        shapes_c(i)%dn = 0.0
      end do
      return
    end if
    
    ! We need to form dn such that a transform_to_physical gives us the
    ! transformed shape function derivatives at the quadrature points of the
    ! supermesh element. A simple eval_dshape(...) isn't going to cut it, so we:
    
    call compute_inverse_jacobian(positions, ele, invj)
    
    do i = 1, size(shapes_c)
      assert(ele_ngi(positions, ele) == ele_ngi(positions_c, i))

      ! First evaluate the transformed shape function derivatives at the
      ! quadrature points of the supermesh element (what we want a
      ! transform_to_physical to give us)
      do j = 1, size(shapes_c(i)%dn, 2)
        shapes_c(i)%dn(:, j, :) = eval_dshape_transformed(base_shape, l_coords(:, j, i), invj)
      end do
    
      ! Then apply the inverse transform on the supermesh element
      call compute_jacobian(positions_c, i, j_c)
      forall(j = 1:size(shapes_c(i)%dn, 1), k = 1:size(shapes_c(i)%dn, 2))
        shapes_c(i)%dn(j, k, :) = matmul(j_c(:, :, k), shapes_c(i)%dn(j, k, :))
      end forall
    end do
    
  end subroutine construct_supermesh_dn_ele
  
  subroutine construct_supermesh_dn_eles(eles, positions, positions_c, l_coords, base_shape, shapes_c)
    type(vector_field), intent(in) :: positions_c
    integer, dimension(ele_count(positions_c)), intent(in) :: eles
    type(vector_field), intent(in) :: positions
    real, dimension(ele_loc(positions, 1), ele_ngi(positions, 1), ele_count(positions_c)), intent(in) :: l_coords
    type(element_type), intent(in) :: base_shape
    type(element_type), dimension(ele_count(positions_c)), intent(inout) :: shapes_c
    
    integer :: ele_c
    
    do ele_c = 1, size(shapes_c)
      call construct_supermesh_dn(eles(ele_c), ele_c, positions, positions_c, l_coords(:, :, ele_c), base_shape, shapes_c(ele_c))
    end do
    
  end subroutine construct_supermesh_dn_eles
  
  subroutine construct_supermesh_dn_ele_ele_c(ele, ele_c, positions, positions_c, l_coords, base_shape, shape_c)
    integer, intent(in) :: ele
    integer, intent(in) :: ele_c
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: positions_c
    real, dimension(ele_loc(positions, ele), ele_ngi(positions, ele)), intent(in) :: l_coords
    type(element_type), intent(in) :: base_shape
    type(element_type), intent(inout) :: shape_c
    
    integer :: i, j
    real, dimension(positions%dim, positions%dim, ele_ngi(positions, ele)) :: invj
    real, dimension(positions_c%dim, positions_c%dim, ele_ngi(positions, ele)) :: j_c
    
    assert(ele_ngi(positions, ele) == ele_ngi(positions_c, ele_c))
    
    if(base_shape%degree == 0) then
      ! This case is nice and easy
      shape_c%dn = 0.0
      return
    end if
    
    ! We need to form dn such that a transform_to_physical gives us the
    ! transformed shape function derivatives at the quadrature points of the
    ! supermesh element. A simple eval_dshape(...) isn't going to cut it, so we:
    
    call compute_inverse_jacobian(positions, ele, invj)
    
    ! First evaluate the transformed shape function derivatives at the
    ! quadrature points of the supermesh element (what we want a
    ! transform_to_physical to give us)
    do i = 1, size(shape_c%dn, 2)
      shape_c%dn(:, i, :) = eval_dshape_transformed(base_shape, l_coords(:, i), invj)
    end do
  
    ! Then apply the inverse transform on the supermesh element
    call compute_jacobian(positions_c, ele_c, j_c)
    forall(i = 1:size(shape_c%dn, 1), j = 1:size(shape_c%dn, 2))
      shape_c%dn(i, j, :) = matmul(j_c(:, :, j), shape_c%dn(i, j, :))
    end forall
    
  end subroutine construct_supermesh_dn_ele_ele_c
    
  function extruded_shape_function(ele_surf, ele_vol, positions_surf, positions_vol, shape_surf, shape_vol, &
    & form_dn) result(shape_surf_ext)
    !!< Extrude a surface shape function
    
    integer, intent(in) :: ele_surf
    integer, intent(in) :: ele_vol
    type(vector_field), intent(in) :: positions_surf
    type(vector_field), intent(in) :: positions_vol
    type(element_type), intent(in) :: shape_surf
    type(element_type), target, intent(in) :: shape_vol
    ! If present and .false., do not form the shape function derivatives
    logical, optional, intent(in) :: form_dn
    type(ele_numbering_type), pointer :: ele_num
    type(element_type) :: shape_surf_ext
    
    integer :: coords, degree, dim, i, loc, ngi
    real, dimension(positions_vol%dim - 1, ele_ngi(positions_vol, ele_vol)) :: positions_gi_vol
    real, dimension(ele_loc(positions_surf, ele_surf), ele_ngi(positions_vol, ele_vol)) :: l_coords
    logical :: lform_dn
    type(quadrature_type), pointer :: quad
    
    lform_dn = .not. present_and_false(form_dn)
    
    quad => shape_vol%quadrature
    
    dim = positions_vol%dim
    loc = shape_surf%loc
    ngi = quad%ngi
    coords = local_coord_count(shape_vol)
    degree = shape_surf%degree
    ele_num => &
         &find_element_numbering(vertices = shape_surf%numbering%vertices, &
         &dimension = dim - 1, degree = degree)
    ! Note that the extruded surface mesh shape function takes its number of
    ! quadrature points from the volume shape function
    call allocate_element(shape_surf_ext, ele_num=ele_num, ngi = ngi)
    shape_surf_ext%degree = degree
    shape_surf_ext%numbering => find_element_numbering(vertices = loc, dimension = dim - 1, degree = degree)
    shape_surf_ext%quadrature = quad
    call incref(quad)

    shape_surf_ext%dn = huge(0.0) 
    assert(.not. associated(shape_surf_ext%dn_s))
    assert(.not. associated(shape_surf_ext%n_s))
    deallocate(shape_surf_ext%spoly)
    nullify(shape_surf_ext%spoly)
    deallocate(shape_surf_ext%dspoly)
    nullify(shape_surf_ext%dspoly)

    select case(degree)
      case(0)
        shape_surf_ext%n = 1.0
      case(1)
        do i = 1, dim - 1
          positions_gi_vol(i, :) = ele_val_at_quad(positions_vol, ele_vol, i)
        end do
        l_coords = local_coords(positions_surf, ele_surf, positions_gi_vol)
        
        if(ele_numbering_family(shape_surf) == FAMILY_SIMPLEX) then
          shape_surf_ext%n = l_coords
        else
          do i = 1, ngi
            shape_surf_ext%n(:, i) = eval_shape(shape_surf, l_coords(:, i))
          end do
        end if
      case default
        do i = 1, dim - 1
          positions_gi_vol(i, :) = ele_val_at_quad(positions_vol, ele_vol, i)
        end do  
        l_coords = local_coords(positions_surf, ele_surf, positions_gi_vol)      
        
        do i = 1, ngi
          shape_surf_ext%n(:, i) = eval_shape(shape_surf, l_coords(:, i))
        end do
    end select

    if(lform_dn) then
      FLAbort("Shape function derivative extrude not yet available")
    end if
  
  end function extruded_shape_function
  
  subroutine generate_supermesh_node_ownership(positions_c, mesh_c, map)
    type(vector_field), intent(in) :: positions_c
    type(mesh_type), intent(in) :: mesh_c
    integer, dimension(:), allocatable, intent(out) :: map
    
    integer :: i
    
    assert(ele_count(mesh_c) > 0)
    allocate(map(ele_count(mesh_c) * ele_loc(mesh_c, 1)))
    do i = 1, ele_count(mesh_c)
      map(ele_nodes(mesh_c, i)) = ele_region_id(positions_c, i)
    end do
    
  end subroutine generate_supermesh_node_ownership
    
  function project_donor_field_to_supermesh_scalar(positions_a, positions_c, field_a) result(field_a_c)
    !!< Project a donor field onto the supermesh
    
    type(vector_field), intent(in) :: positions_a
    type(vector_field), intent(in) :: positions_c
    type(scalar_field), intent(in) :: field_a
    
    type(scalar_field) :: field_a_c
    
    integer, dimension(:), allocatable :: map
    type(mesh_type) :: mesh_a_c
    type(vector_field) :: positions_c_remap
    
    ! Allocate the supermesh field
    assert(ele_count(field_a) > 0)
    mesh_a_c = make_mesh(positions_c%mesh, ele_shape(field_a, 1), continuity = -1)
    call allocate(field_a_c, mesh_a_c, name = trim(field_a%name) // "Supermesh")
    call deallocate(mesh_a_c)
    
    ! Generate the map from nodes in the supermesh field to elements in the
    ! donor field
    call generate_supermesh_node_ownership(positions_c, mesh_a_c, map)
    
    ! We need the "target" positions handed to linear_interpolation to share its
    ! mesh with the supermesh field
    if(positions_c%mesh == mesh_a_c) then
      positions_c_remap = positions_c
      call incref(positions_c_remap)
    else
      call allocate(positions_c_remap, positions_c%dim, mesh_a_c, name = "CoordinateRemap")
      call remap_field(positions_c, positions_c_remap)
    end if
    
    ! Project - consistent interpolation onto the supermesh is lossless
    call linear_interpolation(field_a, positions_a, field_a_c, positions_c_remap, map = map)
    
    ! Cleanup
    deallocate(map)
    call deallocate(positions_c_remap)
    
  end function project_donor_field_to_supermesh_scalar
  
  function project_target_field_to_supermesh_scalar(ele_b, positions_b, positions_c, field_b) result(field_b_c)
    !!< Project a target field onto the supermesh
    
    integer, intent(in) :: ele_b
    type(vector_field), intent(in) :: positions_b
    type(vector_field), intent(in) :: positions_c
    type(scalar_field), intent(in) :: field_b
    
    type(scalar_field) :: field_b_c
    
    type(mesh_type) :: mesh_b_c
    type(vector_field) :: positions_c_remap
    
    ! Allocate the supermesh field
    assert(ele_count(field_b) > 0)
    mesh_b_c = make_mesh(positions_c%mesh, ele_shape(field_b, 1), continuity = -1)
    call allocate(field_b_c, mesh_b_c, name = trim(field_b%name) // "Supermesh")
    call deallocate(mesh_b_c)
        
    ! We need the "target" positions handed to linear_interpolation to share its
    ! mesh with the supermesh field
    if(positions_c%mesh == mesh_b_c) then
      positions_c_remap = positions_c
      call incref(positions_c_remap)
    else
      call allocate(positions_c_remap, positions_c%dim, mesh_b_c, name = "CoordinateRemap")
      call remap_field(positions_c, positions_c_remap)
    end if
    
    ! Project - consistent interpolation onto the supermesh is lossless
    assert(ele_count(mesh_b_c) > 0)
    call linear_interpolation(field_b, positions_b, field_b_c, positions_c_remap, map = spread(ele_b, 1, ele_count(mesh_b_c) * ele_loc(mesh_b_c, 1)))
    
    ! Cleanup
    call deallocate(positions_c_remap)
    
  end function project_target_field_to_supermesh_scalar
  
  subroutine construct_supermesh_ele_single_state(ele_b, positions_a, positions_b, map_ba, &
    & state_a, shape_mesh_b, &
    & state_c, positions_c, shapes_c, &
    & form_dn, single_mesh_state)
    integer, intent(in) :: ele_b
    type(vector_field), intent(in) :: positions_a
    type(vector_field), intent(in) :: positions_b
    type(ilist), intent(in) :: map_ba
    type(state_type), intent(in) :: state_a
    type(mesh_type), intent(in) :: shape_mesh_b
    type(state_type), intent(out) :: state_c
    type(vector_field), intent(out) :: positions_c
    type(element_type), dimension(:), allocatable, intent(out) :: shapes_c
    ! If present and .false., do not form the shape function derivatives
    logical, optional, intent(in) :: form_dn
    ! If present and .true., assume state_a contains fields all on the same mesh
    logical, optional, intent(in) :: single_mesh_state
    
    type(state_type), dimension(1) :: states_a, states_c
    
    states_a = (/state_a/)
    call construct_supermesh_ele(ele_b, positions_a, positions_b, map_ba, &
      & states_a, shape_mesh_b, &
      & states_c, positions_c, shapes_c, &
      & form_dn = form_dn, mesh_sorted_states = single_mesh_state)
    state_c = states_c(1)
  
  end subroutine construct_supermesh_ele_single_state
  
  subroutine construct_supermesh_ele_multiple_states(ele_b, positions_a, positions_b, map_ba, &
    & states_a, shape_mesh_b, &
    & states_c, positions_c, shapes_c, &
    & form_dn, mesh_sorted_states)
    integer, intent(in) :: ele_b
    type(vector_field), intent(in) :: positions_a
    type(vector_field), intent(in) :: positions_b
    type(ilist), intent(in) :: map_ba
    type(state_type), dimension(:), intent(in) :: states_a
    type(mesh_type), intent(in) :: shape_mesh_b
    type(state_type), dimension(size(states_a)), intent(out) :: states_c
    type(vector_field), intent(out) :: positions_c
    type(element_type), dimension(:), allocatable, intent(out) :: shapes_c
    ! If present and .false., do not form the shape function derivatives
    logical, optional, intent(in) :: form_dn
    ! If present and .true., assume states_a is sorted by meshes
    logical, optional, intent(in) :: mesh_sorted_states
    
    integer :: i, j, stat
    integer, dimension(:), allocatable :: map
    type(element_type), pointer :: shape_c
    type(mesh_type) :: mesh_c
    type(mesh_type), pointer :: mesh_a
    type(scalar_field), pointer :: s_field_a
    type(scalar_field) :: s_field_c
    type(state_type), dimension(:), allocatable :: sorted_states_a, sorted_states_c
    type(tensor_field), pointer :: t_field_a
    type(tensor_field) :: t_field_c
    type(vector_field), pointer :: v_field_a
    type(vector_field) :: v_field_c
    
    ! Supermesh
    shape_c => ele_shape(positions_b, ele_b)
    call construct_supermesh(positions_b, ele_b, positions_a, map_ba, shape_c, positions_c)
    call insert(states_c, positions_c, "Coordinate")
    call insert(states_c, positions_c%mesh, "CoordinateMesh")
    
    ! Generate the supermesh shape functions. These are the shape functions of
    ! the target mesh projected onto the supermesh.
    call project_target_shape_to_supermesh(ele_b, &
      & positions_b, shape_mesh_b, positions_c, &
      & shapes_c, form_dn = form_dn)
    
    ! Generate the supermesh fields. These are the fields of the donor mesh
    ! projected onto the supermesh. 
    if(present_and_true(mesh_sorted_states)) then
      do i = 1, size(states_a)
        mesh_a => single_state_mesh(states_a(i), stat = stat)
        if(stat /= 0) cycle
        
        assert(ele_count(mesh_a) > 0)
        shape_c => ele_shape(mesh_a, 1)
        mesh_c = make_mesh(positions_c%mesh, shape_c, continuity = -1, name = mesh_a%name)
        call insert(states_c(i), mesh_c, mesh_c%name)
      
        do j = 1, scalar_field_count(states_a(i))
          s_field_a => extract_scalar_field(states_a(i), j)
          ! We set all fields to have type FIELD_TYPE_NORMAL to keep the
          ! interpolation routines happy. Alternatively, we could modify the
          ! linear interpolation code to handle FIELD_TYPE_CONSTANT fields.
          call allocate(s_field_c, mesh_c, s_field_a%name, field_type = FIELD_TYPE_NORMAL)
          call insert(states_c(i), s_field_c, s_field_c%name)
          call deallocate(s_field_c)
        end do
        
        do j = 1, vector_field_count(states_a(i))
          v_field_a => extract_vector_field(states_a(i), j)
          if(trim(v_field_a%name) == "Coordinate") cycle
          call allocate(v_field_c, v_field_a%dim, mesh_c, v_field_a%name, field_type = FIELD_TYPE_NORMAL)
          call insert(states_c(i), v_field_c, v_field_c%name)
          call deallocate(v_field_c)
        end do
        
        do j = 1, tensor_field_count(states_a(i))
          t_field_a => extract_tensor_field(states_a(i), j)
          call allocate(t_field_c, mesh_c, t_field_a%name, field_type = FIELD_TYPE_NORMAL)
          call insert(states_c(i), t_field_c, t_field_c%name)
          call deallocate(t_field_c)
        end do
        
        call generate_supermesh_node_ownership(positions_c, mesh_c, map)
        call linear_interpolation(states_a(i), states_c(i), map = map)
        deallocate(map)
        
        call deallocate(mesh_c)
      end do
    else
      do i = 1, size(states_a)      
        do j = 1, mesh_count(states_a(i))
          mesh_a => extract_mesh(states_a(i), j)
          assert(ele_count(mesh_a) > 0)
          shape_c => ele_shape(mesh_a, 1)
          mesh_c = make_mesh(positions_c%mesh, shape_c, continuity = -1, name = mesh_a%name)
          call insert(states_c(i), mesh_c, mesh_c%name)
          call deallocate(mesh_c)
        end do
      
        do j = 1, scalar_field_count(states_a(i))
          s_field_a => extract_scalar_field(states_a(i), j)
          mesh_c = extract_mesh(states_c(i), s_field_a%mesh%name)
          call allocate(s_field_c, mesh_c, s_field_a%name, field_type = FIELD_TYPE_NORMAL)
          call insert(states_c(i), s_field_c, s_field_c%name)
          call deallocate(s_field_c)
        end do
        
        do j = 1, vector_field_count(states_a(i))
          v_field_a => extract_vector_field(states_a(i), j)
          if(trim(v_field_a%name) == "Coordinate") cycle
          mesh_c = extract_mesh(states_c(i), v_field_a%mesh%name)
          call allocate(v_field_c, v_field_a%dim, mesh_c, v_field_a%name, field_type = FIELD_TYPE_NORMAL)
          call insert(states_c(i), v_field_c, v_field_c%name)
          call deallocate(v_field_c)
        end do
        
        do j = 1, tensor_field_count(states_a(i))
          t_field_a => extract_tensor_field(states_a(i), j)
          mesh_c = extract_mesh(states_c(i), t_field_a%mesh%name)
          call allocate(t_field_c, mesh_c, t_field_a%name, field_type = FIELD_TYPE_NORMAL)
          call insert(states_c(i), t_field_c, t_field_c%name)
          call deallocate(t_field_c)
        end do
      end do
           
      call sort_states_by_mesh(states_a, sorted_states_a)
      call sort_states_by_mesh(states_c, sorted_states_c)
      
      do i = 1, size(sorted_states_c)
        mesh_c = single_state_mesh(sorted_states_c(i), stat = stat)          
        if(stat == 0) then
          call generate_supermesh_node_ownership(positions_c, mesh_c, map)
          call linear_interpolation(sorted_states_a(i), sorted_states_c(i), map = map)
          deallocate(map)
        end if
        
        call deallocate(sorted_states_a(i))
        call deallocate(sorted_states_c(i))
      end do
      
      deallocate(sorted_states_a)
      deallocate(sorted_states_c)
    end if
    
  end subroutine construct_supermesh_ele_multiple_states
  
  function single_state_mesh(state, stat) result(mesh)
    type(state_type), intent(in) :: state
    integer, optional, intent(out) :: stat
    
    type(mesh_type), pointer :: mesh
    
    type(scalar_field), pointer :: s_field
    type(vector_field), pointer :: v_field
    type(tensor_field), pointer :: t_field
    
    if(present(stat)) stat = 0
    
    if(scalar_field_count(state) > 0) then
      s_field => extract_scalar_field(state, 1)
      mesh => s_field%mesh
    else if(vector_field_count(state) > 0) then
      v_field => extract_vector_field(state, 1)
      mesh => v_field%mesh
    else if(tensor_field_count(state) > 0) then
      t_field => extract_tensor_field(state, 1)
      mesh => t_field%mesh
    else
      if(present(stat)) then
        stat = 1
        return
      else
        ewrite(-1, *) "For state " // trim(state%name)
        FLAbort("No mesh found")
      end if
    end if
    
  end function single_state_mesh
  
  subroutine galerkin_projection_scalars(states_a, positions_a, states_b, positions_b)
    type(state_type), dimension(:), intent(in) :: states_a
    type(vector_field), intent(in) :: positions_a
    type(state_type), dimension(size(states_a)), intent(inout) :: states_b
    type(vector_field), intent(in) :: positions_b
        
    integer :: ele_b, ele_c, field_count, i, j
    type(csr_matrix), pointer :: mass_matrix
    type(element_type), dimension(:), allocatable :: shapes_c
    type(ilist), dimension(ele_count(positions_b)) :: map_ba
    type(mesh_type), pointer :: mesh_b
    type(scalar_field), dimension(:), allocatable :: rhs
    type(scalar_field), pointer :: s_field_b
    type(state_type) :: state_c
    type(vector_field) :: positions_c

    call intersector_set_dimension(positions_b%dim)
    
    map_ba = intersection_finder(positions_b, positions_a)

    do i = 1, size(states_b)
      field_count = scalar_field_count(states_b(i))
      if(field_count == 0) cycle
    
      s_field_b => extract_scalar_field(states_b(i), 1)
      mesh_b => s_field_b%mesh
    
      select case(mesh_b%continuity)
        case(0)
          mass_matrix => get_mass_matrix(states_b(i), mesh_b)
          allocate(rhs(scalar_field_count(states_b(i))))
          do j = 1, field_count
            call allocate(rhs(j), mesh_b, "GalerkinProjectionRHS" // int2str(j))
            call zero(rhs(j))
          end do
        
          do ele_b = 1, ele_count(positions_b)
            call construct_supermesh_ele(ele_b, positions_a, positions_b, map_ba(ele_b), &
              & states_a(i), mesh_b, &
              & state_c, positions_c, shapes_c, &
              & form_dn = .false., single_mesh_state = .true.)

            do ele_c = 1, ele_count(positions_c)
              call assemble_galerkin_projection_scalars_ele(ele_c, ele_b, positions_c, state_c, shapes_c(ele_c), rhs)
            end do

            call deallocate(state_c)
            call deallocate(positions_c)
            do j = 1, size(shapes_c)
              call deallocate(shapes_c(j))
            end do
            deallocate(shapes_c)
          end do

          do j = 1, field_count
            s_field_b => extract_scalar_field(states_b(i), j)
            call petsc_solve(s_field_b, mass_matrix, rhs(j), &
              & option_path = trim(complete_field_path(s_field_b%option_path)) // "/galerkin_projection/continuous")

            call deallocate(rhs(j))
          end do
          deallocate(rhs)
        case(-1)
          do ele_b = 1, ele_count(positions_b)
            call construct_supermesh_ele(ele_b, positions_a, positions_b, map_ba(ele_b), &
              & states_a(i), mesh_b, &
              & state_c, positions_c, shapes_c, &
              & form_dn = .false., single_mesh_state = .true.)
              
            call solve_galerkin_projection_scalars_dg_ele(ele_b, positions_b, positions_c, mesh_b, states_b(i), state_c, shapes_c)
              
            call deallocate(state_c)
            call deallocate(positions_c)
            do j = 1, size(shapes_c)
              call deallocate(shapes_c(j))
            end do
            deallocate(shapes_c)
          end do
        case default
          ewrite(-1, "(a,i0)") "For mesh continuity ", mesh_b%continuity
          FLAbort("Unrecognised mesh continuity")
      end select
    end do
    
    do i = 1, size(map_ba)
      call deallocate(map_ba(i))
    end do

  end subroutine galerkin_projection_scalars  
  
  subroutine assemble_galerkin_projection_scalars_ele(ele, ele_out, positions, state, shape, rhs)
    integer, intent(in) :: ele
    integer, intent(in) :: ele_out
    type(vector_field), intent(in) :: positions
    type(state_type), intent(in) :: state
    type(element_type), intent(in) :: shape
    type(scalar_field), intent(inout), dimension(:) :: rhs
    
    integer :: field
    real, dimension(ele_ngi(positions, ele)) :: detwei
    type(scalar_field), pointer :: s_field
    
    assert(size(rhs) == scalar_field_count(state))

    call transform_to_physical(positions, ele, detwei = detwei)
    do field = 1, size(rhs)
      s_field => extract_scalar_field(state, field)
      call addto(rhs(field), ele_nodes(rhs(field), ele_out), &
        & shape_rhs(shape, detwei * ele_val_at_quad(s_field, ele)))
    end do
    
  end subroutine assemble_galerkin_projection_scalars_ele
  
  subroutine solve_galerkin_projection_scalars_dg_ele(ele_b, positions_b, positions_c, mesh_b, state_b, state_c, shapes_c)
    integer, intent(in) :: ele_b
    type(vector_field), intent(in) :: positions_b
    type(vector_field), intent(in) :: positions_c
    type(mesh_type), intent(in) :: mesh_b
    type(state_type), intent(in) :: state_b
    type(state_type), intent(in) :: state_c
    type(element_type), dimension(ele_count(positions_c)), intent(in) :: shapes_c
    
    integer :: i, j
    real, dimension(ele_loc(mesh_b, ele_b), scalar_field_count(state_c)) :: little_rhs
    real, dimension(ele_ngi(positions_b, ele_b)) :: detwei
    real, dimension(ele_loc(mesh_b, ele_b), ele_loc(mesh_b, ele_b)) :: little_mass
    type(scalar_field), pointer :: s_field_b
    
    call transform_to_physical(positions_b, ele_b, detwei = detwei)
      
    little_mass = shape_shape(ele_shape(mesh_b, ele_b), ele_shape(mesh_b, ele_b), detwei)
    little_rhs = 0.0
    do i = 1, scalar_field_count(state_b)
      do j = 1, ele_count(positions_c)
        call assemble_galerkin_projection_scalars_dg_ele(j, positions_c, state_c, shapes_c(j), little_rhs)
      end do
    end do

    call solve(little_mass, little_rhs)

    do i = 1, scalar_field_count(state_b)
      s_field_b => extract_scalar_field(state_b, i)
      call set(s_field_b, ele_nodes(s_field_b, ele_b), little_rhs(:, i))
    end do
    
  end subroutine solve_galerkin_projection_scalars_dg_ele
  
  subroutine assemble_galerkin_projection_scalars_dg_ele(ele, positions, state, shape, little_rhs)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(state_type), intent(in) :: state
    type(element_type), intent(in) :: shape
    real, dimension(shape%loc, scalar_field_count(state)), intent(inout) :: little_rhs
    
    integer :: i
    real, dimension(ele_ngi(positions, ele)) :: detwei
    type(scalar_field), pointer :: s_field
    
    call transform_to_physical(positions, ele, detwei = detwei)
      
    do i = 1, scalar_field_count(state)
      s_field => extract_scalar_field(state, i)
      little_rhs(:, i) = little_rhs(:, i) + shape_rhs(shape, detwei * ele_val_at_quad(s_field, ele))
    end do
    
  end subroutine assemble_galerkin_projection_scalars_dg_ele
  
  function compute_inner_product_sa(positions_a, positions_b, a, b) result(val)
    type(vector_field), intent(in) :: positions_a
    type(vector_field), intent(in) :: positions_b
    type(scalar_field), intent(in) :: a
    type(scalar_field), intent(in) :: b
    
    real :: val
        
    integer :: ele_b, ele_c
    type(ilist), dimension(ele_count(positions_b)) :: map_ba
    type(scalar_field) :: a_c, b_c
    type(vector_field) :: positions_c
    
    val = 0.0
    
    call intersector_set_dimension(positions_a%dim)
    
    map_ba = intersection_finder(positions_b, positions_a)
    
    do ele_b = 1, ele_count(positions_b)
      ! Supermesh
      call construct_supermesh(positions_b, ele_b, positions_a, map_ba(ele_b), ele_shape(positions_b, ele_b), positions_c)
      if(ele_count(positions_c) == 0) then
        call deallocate(positions_c)
        cycle
      end if
      
      ! Project a onto the supermesh
      a_c = project_donor_field_to_supermesh(positions_a, positions_c, a)
      ! Project b onto the supermesh
      b_c = project_target_field_to_supermesh(ele_b, positions_b, positions_c, b)
        
      do ele_c = 1, ele_count(positions_c)        
        ! Compute the contribution to the inner product
        call add_inner_product_ele(ele_c, positions_c, a_c, b_c, val)
      end do
        
      ! Cleanup
      call deallocate(positions_c)
      call deallocate(a_c)
      call deallocate(b_c)
    end do
    
    do ele_b = 1, ele_count(positions_b)
      call deallocate(map_ba(ele_b))
    end do
    
  end function compute_inner_product_sa
  
  subroutine add_inner_product_ele(ele, positions, a, b, val)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: a
    type(scalar_field), intent(in) :: b
    real, intent(inout) :: val
    
    real, dimension(ele_ngi(positions, ele)) :: detwei
    
    call transform_to_physical(positions, ele, detwei = detwei)
          
    val = val + dot_product(ele_val(a, ele), matmul(&
        &  shape_shape(ele_shape(a, ele), ele_shape(b, ele), detwei), ele_val(b, ele)))

  end subroutine add_inner_product_ele

end module supermesh_assembly
