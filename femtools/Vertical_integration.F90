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
!    C.Pain@Imperial.ac.uk
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

module vertical_integration
  !!< routines for doing vertical integration
  use spud
  use Fields
  use Sparse_tools
  use Sparsity_patterns
  use elements
  use boundary_conditions
  use solvers
  implicit none

  private

  type(csr_matrix), save :: surface_mass, surface_laplacian, &
       volume2surface, &
       CT_vert, vertical_laplacian, volume_mass, &
       horizontal_laplacian

  type(scalar_field), save :: height
  type(mesh_type), pointer :: v_mesh, s_mesh

  logical :: initialised = .false.

contains

  subroutine initialise_vertical_operators(pressure,positions)
    type(scalar_field), intent(inout) :: pressure
    type(vector_field), intent(inout) :: positions
    !local variables
    type(scalar_field) :: unity

    if(.not.initialised) then
       !set up operators
       call get_vertical_operators(pressure,positions)

       !get height
       call allocate(height,pressure%mesh)
       call allocate(unity,pressure%mesh,field_type=FIELD_TYPE_CONSTANT)
       unity%val = 1.0
       call get_vertical_integral(unity,height)

       initialised = .true.
    end if
  end subroutine initialise_vertical_operators

  subroutine deallocate_vertical_operators()
    if(initialised) then
       call deallocate(height)
       call deallocate(surface_mass)
       call deallocate(surface_laplacian)
       call deallocate(volume2surface)
       call deallocate(ct_vert)
       call deallocate(vertical_laplacian)
       call deallocate(horizontal_laplacian)
       call deallocate(volume_mass)
       initialised = .false.
    end if
  end subroutine deallocate_vertical_operators

  subroutine apply_vertical_preconditioner(X,rhs)
    real, dimension(:), intent(in) :: rhs
    real, dimension(:), intent(out) :: X
    !local variables
    !---------------------
    type(scalar_field) :: hlap_rhs, hlap_rhs_v_int, g, g_vol, rhs_vol, &
         X_vol

    call allocate(hlap_rhs,v_mesh)
    call allocate(hlap_rhs_v_int,s_mesh)
    call allocate(g,s_mesh)
    call allocate(g_vol,v_mesh)
    call allocate(rhs_vol,v_mesh)
    call allocate(X_vol,v_mesh)

    !get horizontal laplacian of rhs -- Hlap_rhs
    call mult(hlap_rhs%val,horizontal_laplacian,rhs)

    !combine with rhs
    hlap_rhs%val = rhs - hlap_rhs%val

    !vertically-average Hlap_rhs -- hlap_rhs_v_int
    call get_vertical_integral(hlap_rhs, hlap_rhs_v_int)

    !divide Hlap_rhs by height    
    hlap_rhs_v_int%val = hlap_rhs_v_int%val/height%val

    !solve horizontal equation for g
    g%options%abs_error = 1.0e-8
    g%options%max_its = 10000
    call zero(g)
    call petsc_solve(g, surface_laplacian, hlap_rhs_v_int)

    !extrapolate g to volume g_vol
    call apply_vertical_extrapolation(g_vol, g)

    !take horizontal laplacian of g_vol and put in rhs_vol
    call mult(rhs_vol%val,horizontal_laplacian,g_vol%val)
    rhs_vol%val = rhs_vol%val + rhs
    
    !solve vertical equation for X_vol
    x_vol%options%abs_error = 1.0e-8
    x_vol%options%max_its = 10000
    call zero(x_vol)
    call petsc_solve(x_vol, vertical_laplacian, x_vol)

    !copy X_vol into X
    X = x_vol%val

    call deallocate(hlap_rhs)
    call deallocate(hlap_rhs_v_int)
    call deallocate(g_vol)
    call deallocate(rhs_vol)
    call deallocate(X_vol)

  end subroutine apply_vertical_preconditioner

  subroutine get_vertical_integral(field, integral_field)
    !compute vertical integral
    !assumes field has been multiplied by the mass matrix as that is
    !what we use in operators

    type(scalar_field), intent(inout) :: field
    type(scalar_field), intent(inout) :: integral_field

    !local variables
    type(scalar_field) :: integral_field_rhs, integral_field_vol, &
         field2

    call allocate(integral_field_rhs,field%mesh)
    call allocate(integral_field_vol,field%mesh)
    call allocate(field2,field%mesh)

    !divide by mass matrix - field --> field2

    field2%options%abs_error = 1.0e-8
    field2%options%max_its = 10000
    call zero(field2)
    call petsc_solve(field2, volume_mass, field)

    !multiply by C - field2 --> integral_field_rhs

    call mult_T(integral_field_rhs%val, CT_vert, integral_field%val)
    
    !divide by vertical Laplacian - integral_field_rhs --> integral_field_vol

    integral_field_vol%options%abs_error = 1.0e-8
    integral_field_vol%options%max_its = 10000
    call zero(integral_field_vol)
    call petsc_solve(integral_field_vol, vertical_laplacian, &
         integral_field_rhs)

    !push onto surface mesh - integral_field_vol --> integral_field_vol

    call mult(integral_field%val, volume2surface, &
         integral_field_vol%val)

    call deallocate(integral_field_rhs)
    call deallocate(integral_field_vol)
    call deallocate(field2)

  end subroutine get_vertical_integral

  subroutine apply_vertical_extrapolation(field, surface_field)
    
    type(scalar_field), intent(inout) :: field
    type(scalar_field), intent(inout) :: surface_field

    !local variables
    type(scalar_field) :: field2, field3, field4, surface_field2

    call allocate(field2,field%mesh)
    call allocate(field3,field%mesh)
    call allocate(field4,field%mesh)
    call allocate(surface_field2,surface_field%mesh)

    !multiply by surface mass - surface_field --> surface_field2

    call mult_T(surface_field2%val, surface_mass, surface_field%val)

    !put into volume field - surface_field2 --> field4

    call mult_T(field4%val, volume2surface, &
         surface_field2%val)

    !divide by vertical Laplacian - field4 --> field3

    call petsc_solve(field3, vertical_laplacian, field4)

    !compute vertical gradient - field3 --> field2

    call mult(field2%val, CT_vert, field3%val)

    !divide by volume mass - field2 --> field

    ewrite(-1,*) 'warning, setting field options in extrapolation'
    field%options%abs_error = 1.0e-8
    field%options%max_its = 10000
    call zero(field)

    call petsc_solve(field, volume_mass, field2)

    call deallocate(field4)
    call deallocate(field3)
    call deallocate(field2)
    
  end subroutine apply_vertical_extrapolation

  subroutine get_vertical_operators(&
       & pressure, positions)

    type(scalar_field), intent(inout), target :: pressure
    type(vector_field), intent(in) :: positions
    !local variables
    integer :: ele !counting variable for elements
    integer :: node !counting variable for nodes
    type(csr_sparsity) :: CT_vert_sparsity, volume2surface_sparsity, &
         & surface_mass_sparsity, vertical_laplacian_sparsity
    type(scalar_field) :: field
    !stuff for storing top surface mesh
    integer, dimension(:), pointer :: s_e_list, s_n_list
    type(scalar_field) :: field_s
    type(vector_field) :: positions_s
    !variables for getting surface ids
    integer, dimension(:), allocatable:: surface_ids
    integer shape_option(2)

    call allocate(field,pressure%mesh)

    if(field%mesh%continuity<0) then
       FLAbort('Only works for continuous fields')
    end if

    ! Get vector of surface ids at top
    if(have_option('/geometry/ocean_boundaries/top_surface_ids')) then
       shape_option=option_shape('/geometry/ocean_boundaries/top_surface_ids')
       allocate(surface_ids(1:shape_option(1)))
       call get_option('/geometry/ocean_boundaries/top_surface_ids', &
            &surface_ids)
       ! Add boundary condition that marks the top of the domain
       call add_boundary_condition(field, name="TopSurface", &
            type="neumann", boundary_ids=surface_ids)
       deallocate(surface_ids)
    else
       ewrite(-1,*) "No top surface ids so using id 1"
       call add_boundary_condition(field, name="TopSurface", &
            type="neumann", boundary_ids=(/1/))
    end if
    v_mesh => pressure%mesh
    call get_boundary_condition(field, "TopSurface", &
         surface_element_list=s_e_list, surface_mesh=s_mesh, &
         surface_node_list=s_n_list)

    call allocate(positions_s,3,s_mesh)
    call allocate(field_s,s_mesh)

    call remap_vector_field_to_surface(&
         positions, positions_s, s_e_list)

    ! Calculate sparsities based on connectivity of
    ! field
    vertical_laplacian_sparsity= &
         make_sparsity(field%mesh, field%mesh, &
         name='VerticalLaplacianSparsity')
    call allocate(vertical_laplacian,Vertical_laplacian_Sparsity)
    call allocate(volume_mass,Vertical_laplacian_Sparsity)
    call zero(Vertical_laplacian)
    call zero(Volume_mass)

    CT_vert_sparsity=make_sparsity(field%mesh, field%mesh, &
         name='CT_vert_sparsity')
    call allocate(CT_vert,CT_vert_sparsity)
    call zero(CT_vert)

    surface_mass_sparsity= &
         make_sparsity(s_mesh, s_mesh, name ='surface_mesh_sparsity')
    call allocate(surface_mass,surface_mass_sparsity)
    call allocate(surface_laplacian, surface_mass_sparsity)
    call zero(surface_mass)
    call zero(surface_laplacian)
    
    volume2surface_sparsity = make_sparsity(&
         s_mesh, field%mesh, 'volume2surface_sparsity')
    call allocate(volume2surface, volume2surface_sparsity)
    call zero(volume2surface)

    do ele = 1, element_count(field) 
       call assemble_vertical_operators_element(&
            positions,field,ele)
    end do

    do ele = 1, element_count(s_mesh)
       call assemble_vertical_operators_surface_element(&
            positions_s, field_s,ele)
    end do

    do node = 1, node_count(s_mesh)
       call set(volume2surface,node,s_n_list(node),1.0)
    end do

    call deallocate(positions_s)
    call deallocate(field_s)

    ! Get vector of surface ids at bottom
    if(have_option('/geometry/ocean_boundaries/top_surface_ids')) then
       shape_option=option_shape('/geometry/ocean_boundaries/bottom_surface_ids')
       allocate(surface_ids(1:shape_option(1)))
       call get_option('/geometry/ocean_boundaries/bottom_surface_ids', &
            &surface_ids)
       ! Add boundary condition that marks the top of the domain
       call add_boundary_condition(field, name="BottomSurface", &
            type="dirichlet", boundary_ids=surface_ids)
       deallocate(surface_ids)
    else
       ewrite(-1,*) "No top surface ids so using id 2"
       call add_boundary_condition(field, name="BottomSurface", &
            type="dirichlet", boundary_ids=(/2/))
    end if

    call get_boundary_condition(field, "BottomSurface", &
         surface_element_list=s_e_list, surface_mesh=s_mesh, &
         surface_node_list=s_n_list)

    do node = 1, size(s_n_list)
       call set(vertical_laplacian, s_n_list(node), &
            s_n_list(node), INFINITY)
    end do
    
    call set(surface_laplacian, 1, 1, INFINITY)
    
    call deallocate(field)

  end subroutine get_vertical_operators

  subroutine assemble_vertical_operators_element(&
            positions,field,ele)
    type(scalar_field), intent(in) :: field
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: ele
    ! Locations of quadrature points
    real, dimension(positions%dim,ele_ngi(positions,ele)) :: X_quad
    ! Derivatives of shape function:
    real, dimension(ele_loc(field,ele), &
         ele_ngi(field,ele), positions%dim) :: dshape_field, dshape_v_field
    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,ele)) :: detwei    
    ! Node numbers of field element.
    integer, dimension(:), pointer :: ele_field
    ! Shape functions.
    type(element_type), pointer :: shape_field
    ! Local vertical laplacian matrix 
    real, dimension(ele_loc(field, ele), ele_loc(field, ele)) :: zz_mat
    real, dimension(ele_loc(field, ele), ele_loc(field, ele)) :: mass_mat
    real, dimension(ele_loc(field, ele), ele_loc(field, ele)) :: ct_mat
    real, dimension(ele_loc(field, ele), ele_loc(field, ele)) :: laph_mat
    ! normal vector
    real,dimension(positions%dim) :: n
    integer :: i,j

    !--------------------------
    
    if(positions%dim==3) n = (/0.0,0.0,1.0/)
    if(positions%dim==2) n = (/0.0,1.0/)

    ele_field=>ele_nodes(field, ele)
    shape_field=>ele_shape(field, ele)

    ! Locations of quadrature points.
    X_quad=ele_val_at_quad(positions, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, shape_field, &
         dshape=dshape_field, detwei=detwei)

    ! Extract component in "up" direction
    forall(i=1:ele_loc(field,ele),j=1:ele_ngi(field,ele))
       dshape_v_field(i,j,:) = dshape_field(i,j,:)*n
    end forall

    ! Local assembly:
    !it's actually minus the laplacian
    zz_mat= dshape_dot_dshape(dshape_v_field, dshape_v_field, detwei)
    !it's actually minus the gradient
    ct_mat=sum(dshape_shape(dshape_v_field,shape_field,detwei),1)
    !horizontal laplacian
    laph_mat= dshape_dot_dshape(dshape_field,dshape_field, detwei) - zz_mat
    !mass
    mass_mat=shape_shape(shape_field,shape_field,detwei)

    ! Global assembly: 
    call addto(vertical_laplacian, ele_field, &
         ele_field, zz_mat) 
    call addto(horizontal_laplacian, ele_field, &
         ele_field, laph_mat) 
    call addto(ct_vert, ele_field, ele_field, &
         ct_mat)
    call addto(volume_mass, ele_field, &
         ele_field, mass_mat) 

  end subroutine assemble_vertical_operators_element

  subroutine assemble_vertical_operators_surface_element(&
       positions, field, ele)
    type(scalar_field), intent(in) :: field
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: ele
    ! Locations of quadrature points
    real, dimension(positions%dim,ele_ngi(positions,ele)) :: X_quad
    ! Derivatives of shape function:
    real, dimension(ele_loc(field,ele), &
         ele_ngi(field,ele), positions%dim) :: dshape_field
    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,ele)) :: detwei    
    ! Node numbers of field element.
    integer, dimension(:), pointer :: ele_field
    ! Shape functions.
    type(element_type), pointer :: shape_field
    ! Local vertical laplacian matrix 
    real, dimension(ele_loc(field, ele), ele_loc(field, ele)) :: lap_mat
    real, dimension(ele_loc(field, ele), ele_loc(field, ele)) :: mass_mat

    !--------------------------
    
    ele_field=>ele_nodes(field, ele)
    shape_field=>ele_shape(field, ele)

    ! Locations of quadrature points.
    X_quad=ele_val_at_quad(positions, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(positions, ele, shape_field, &
         dshape=dshape_field, detwei=detwei)

    ! Local assembly:
    lap_mat=-dshape_dot_dshape(dshape_field, dshape_field, detwei)
    mass_mat=shape_shape(shape_field,shape_field,detwei)
    
    ! Global assembly:
    call addto(surface_laplacian, ele_field, ele_field, lap_mat)
    call addto(surface_mass, ele_field, ele_field, mass_mat)

  end subroutine assemble_vertical_operators_surface_element
  
end module vertical_integration
