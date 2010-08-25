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

module divergence_matrix_cg

  use quadrature
  use fields
  use field_derivatives
  use state_module
  use futils
  use fetools
  use spud
  use boundary_conditions
  use global_parameters, only: OPTION_PATH_LEN
  use transform_elements
  use fldebug
  use field_options, only: complete_field_path
  use upwind_stabilisation
  use equation_of_state

  implicit none

  private
  public :: assemble_divergence_matrix_cg, assemble_compressible_divergence_matrix_cg

  ! Stabilisation schemes
  integer, parameter :: STABILISATION_NONE = 0, &
    & STABILISATION_STREAMLINE_UPWIND = 1, STABILISATION_SUPG = 2
  ! Stabilisation scheme
  integer :: stabilisation_scheme
  integer :: nu_bar_scheme
  real :: nu_bar_scale

contains

    subroutine assemble_divergence_matrix_cg(CT_m, state, ct_rhs, & 
                                             test_mesh, field, option_path, &
                                             div_mass, grad_mass, &
                                             div_mass_lumped,&
                                             grad_mass_lumped, get_ct) 

      ! inputs/outputs
      ! bucket full of fields
      type(state_type), intent(inout) :: state

      ! the velocity divergence gradient matrices
      type(block_csr_matrix), intent(inout) :: CT_m

      type(scalar_field), intent(inout), optional :: ct_rhs

      type(mesh_type), intent(in) :: test_mesh

      type(vector_field), intent(inout) :: field

      character(len=*), intent(in), optional :: option_path

      type(csr_matrix), intent(inout), optional :: div_mass, grad_mass
      
      type(scalar_field), intent(inout), optional :: div_mass_lumped, grad_mass_lumped

      logical, intent(in), optional :: get_ct

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
      type(vector_field), pointer :: coordinate
      character(len=OPTION_PATH_LEN) :: l_option_path

      ! integrate by parts
      logical :: integrate_by_parts

      integer, dimension(:,:), allocatable :: field_bc_type
      type(vector_field) :: field_bc

      real, dimension(:,:), allocatable :: div_mass_mat, grad_mass_mat
      
      integer :: stat

      logical :: l_get_ct

      ! =============================================================
      ! Subroutine to construct the matrix CT_m (a.k.a. C1/2/3T).
      ! =============================================================

      ewrite(2,*) 'In assemble_divergence_matrix_cg'

      if(present(get_ct)) then
        l_get_ct = get_ct
      else
        l_get_ct = .true.
      end if

      if(present(option_path)) then
        l_option_path = trim(option_path)
      else
        l_option_path = trim(field%option_path)
      end if

      coordinate=>extract_vector_field(state, "Coordinate")

      integrate_by_parts=have_option(trim(complete_field_path(l_option_path, stat))//&
          &"/spatial_discretisation/continuous_galerkin/integrate_continuity_by_parts")&
          .or. have_option(trim(complete_field_path(l_option_path, stat))//&
          &"/integrate_divergence_by_parts")&
          .or. have_option(trim(complete_field_path(l_option_path, stat))//&
          &"/spatial_discretisation/continuous_galerkin/integrate_divergence_by_parts")&
          .or. have_option(trim(complete_field_path(l_option_path, stat))//&
          &"/spatial_discretisation/discontinuous_galerkin")

      ewrite(2,*) "Divergence is integrated by parts: ", integrate_by_parts

      if(present(ct_rhs)) call zero(ct_rhs)

      if (l_get_ct) then

         ! Clear memory of arrays being designed
         call zero(CT_m)
         if(present(div_mass)) call zero(div_mass)
         if(present(grad_mass)) call zero(grad_mass)
         if(present(div_mass_lumped)) call zero(div_mass_lumped)
         if(present(grad_mass_lumped)) call zero(grad_mass_lumped)

         allocate(dfield_t(ele_loc(field, 1), ele_ngi(field, 1), field%dim), &
              dtest_t(ele_loc(test_mesh, 1), ele_ngi(test_mesh, 1), field%dim), &
              ele_mat(field%dim, ele_loc(test_mesh, 1), ele_loc(field, 1)), &
              detwei(ele_ngi(field, 1)), &
              grad_mass_mat(ele_loc(field, 1), ele_loc(field, 1)), &
              div_mass_mat(ele_loc(test_mesh, 1), ele_loc(test_mesh, 1)))
         
         do ele=1, element_count(test_mesh)

            test_nodes=>ele_nodes(test_mesh, ele)
            field_nodes=>ele_nodes(field, ele)

            test_shape=>ele_shape(test_mesh, ele)
            field_shape=>ele_shape(field, ele)
            
            if(integrate_by_parts) then
               ! transform the pressure derivatives into physical space
               ! (and get detwei)
               call transform_to_physical(coordinate, ele, test_shape, &
                    dshape=dtest_t, detwei=detwei)

               ele_mat = -dshape_shape(dtest_t, field_shape, detwei)
            else
               ! transform the velociy derivatives into physical space
               ! (and get detwei)
               call transform_to_physical(coordinate, ele, field_shape, &
                    dshape=dfield_t, detwei=detwei)
               
               ele_mat = shape_dshape(test_shape, dfield_t, detwei)
            end if
            
            do dim = 1, field%dim
               call addto(ct_m, 1, dim, test_nodes, field_nodes, ele_mat(dim,:,:))
            end do
            
            if(present(div_mass).or.present(div_mass_lumped)) then
               
               div_mass_mat = shape_shape(test_shape, test_shape, detwei)
               
               if(present(div_mass)) then
                  
                  call addto(div_mass, test_nodes, test_nodes, div_mass_mat)
                  
               end if
               
               if(present(div_mass_lumped)) then
                  
                  call addto(div_mass_lumped, test_nodes, sum(div_mass_mat, 2))
                  
               end if
               
            end if
            
            if(present(grad_mass).or.present(grad_mass_lumped)) then
               
               grad_mass_mat = shape_shape(field_shape, field_shape, detwei)
               
               if(present(grad_mass)) then
                  
                  call addto(grad_mass, field_nodes, field_nodes, grad_mass_mat)
                  
               end if
          
               if(present(grad_mass_lumped)) then
                  
                  call addto(grad_mass_lumped, field_nodes, sum(grad_mass_mat, 2))
                  
               end if
               
            end if
        
         end do

      end if

      if(integrate_by_parts) then

        allocate(detwei_bdy(face_ngi(field, 1)), &
                normal_bdy(field%dim, face_ngi(field, 1)))
        allocate(field_nodes_bdy(field%mesh%faces%shape%loc))
        allocate(test_nodes_bdy(test_mesh%faces%shape%loc))
        allocate(ele_mat_bdy(field%dim, face_loc(test_mesh, 1), face_loc(field, 1)))

        assert(surface_element_count(test_mesh)==surface_element_count(field))
        allocate(field_bc_type(field%dim, surface_element_count(field)))
        call get_entire_boundary_condition(field, (/ &
          "weakdirichlet ", &
          "no_normal_flow", &
          "internal      ", &
          "free_surface  "/), field_bc, field_bc_type)

        do sele = 1, surface_element_count(test_mesh)

          if(any(field_bc_type(:,sele)==2)&
               .or.any(field_bc_type(:,sele)==3)&
               .or.any(field_bc_type(:,sele)==4)) cycle
          
          test_shape=>face_shape(test_mesh, sele)
          field_shape=>face_shape(field, sele)

          test_nodes_bdy=face_global_nodes(test_mesh, sele)
          field_nodes_bdy=face_global_nodes(field, sele)

          call transform_facet_to_physical(coordinate, sele, &
              &                          detwei_f=detwei_bdy,&
              &                          normal=normal_bdy) 

          ele_mat_bdy = shape_shape_vector(test_shape, field_shape, detwei_bdy, normal_bdy)

          do dim = 1, field%dim
            if((field_bc_type(dim, sele)==1).and.present(ct_rhs)) then
              call addto(ct_rhs, test_nodes_bdy, &
                          -matmul(ele_mat_bdy(dim,:,:), &
                          ele_val(field_bc, dim, sele)))
            else
               if (l_get_ct) then
                  call addto(ct_m, 1, dim, test_nodes_bdy, field_nodes_bdy, &
                       ele_mat_bdy(dim,:,:))
               end if
            end if
          end do

        end do

        call deallocate(field_bc)
        deallocate(field_bc_type)
        deallocate(detwei_bdy, normal_bdy)
        deallocate(test_nodes_bdy, field_nodes_bdy)

      end if

    end subroutine assemble_divergence_matrix_cg

    subroutine assemble_compressible_divergence_matrix_cg(ctp_m, state, ct_rhs)
      
      ! inputs/outputs
      ! bucket full of fields
      type(state_type), dimension(:), intent(inout) :: state

      ! the compressible divergence matrix
      type(block_csr_matrix), intent(inout) :: ctp_m

      type(scalar_field), intent(inout), optional :: ct_rhs

      if((size(state)==1).and.(.not.has_scalar_field(state(1), "MaterialVolumeFraction"))) then
      
        call assemble_1mat_compressible_divergence_matrix_cg(ctp_m, state(1), ct_rhs)
        
      else
      
        FLExit("Multimaterial compressible continuous_galerkin pressure not possible.")
        
      end if
    
    end subroutine assemble_compressible_divergence_matrix_cg

    subroutine assemble_1mat_compressible_divergence_matrix_cg(ctp_m, state, ct_rhs) 

      ! inputs/outputs
      ! bucket full of fields
      type(state_type), intent(inout) :: state

      ! the compressible divergence matrix
      type(block_csr_matrix), intent(inout) :: ctp_m

      type(scalar_field), intent(inout), optional :: ct_rhs

      ! local
      type(mesh_type), pointer :: test_mesh

      type(vector_field), pointer :: field

      integer, dimension(:), pointer :: test_nodes, field_nodes
      integer, dimension(:), allocatable :: test_nodes_bdy, field_nodes_bdy

      real, dimension(:,:,:), allocatable :: ele_mat, ele_mat_bdy
      type(element_type), pointer :: field_shape, test_shape_ptr, density_shape
      type(element_type) :: test_shape
      real, dimension(:,:,:), allocatable :: dfield_t, dtest_t, ddensity_t
      real, dimension(:), allocatable :: detwei, detwei_bdy, density_bdy
      real, dimension(:,:,:), allocatable :: j_mat
      real, dimension(:,:), allocatable :: normal_bdy
      
      real, dimension(:), allocatable :: density_at_quad, olddensity_at_quad
      real, dimension(:,:), allocatable :: density_grad_at_quad, nlvelocity_at_quad

      ! loop integers
      integer :: ele, sele, dim

      ! pointer to coordinates
      type(vector_field), pointer :: coordinate, nonlinearvelocity, velocity
      type(scalar_field), pointer :: pressure, density, olddensity
      real :: theta, dt

      ! integrate by parts
      logical :: integrate_by_parts

      integer, dimension(:,:), allocatable :: field_bc_type
      type(vector_field) :: field_bc

      integer, dimension(:), allocatable :: density_bc_type
      type(scalar_field) :: density_bc

      real, dimension(:,:), allocatable :: mass_p_mat
      
      integer :: stat

      ! =============================================================
      ! Subroutine to construct the matrix ctp_m (a.k.a. C1/2/3TP).
      ! =============================================================

      ewrite(2,*) 'In assemble_1mat_compressible_divergence_matrix_cg'

      coordinate=> extract_vector_field(state, "Coordinate")
      
      density => extract_scalar_field(state, "Density")
      olddensity => extract_scalar_field(state, "OldDensity")
      
      pressure => extract_scalar_field(state, "Pressure")
      
      velocity=>extract_vector_field(state, "Velocity")
      nonlinearvelocity=>extract_vector_field(state, "NonlinearVelocity") ! maybe this should be updated after the velocity solve?
      
      integrate_by_parts=have_option(trim(complete_field_path(density%option_path, stat))//&
          &"/spatial_discretisation/continuous_galerkin/advection_terms/integrate_advection_by_parts")&
          .or. have_option(trim(complete_field_path(velocity%option_path, stat))//&
          &"/spatial_discretisation/discontinuous_galerkin")

      ewrite(2,*) "Compressible divergence is integrated by parts: ", integrate_by_parts

      if(have_option(trim(density%option_path) // "/prognostic/spatial_discretisation/&
                      &continuous_galerkin/stabilisation/streamline_upwind")) then
        ewrite(2, *) "Streamline upwind stabilisation"
        FLExit("SU stabilisation broken with continuity at the moment.")
        stabilisation_scheme = STABILISATION_STREAMLINE_UPWIND
        call get_upwind_options(trim(density%option_path) // & 
                                "/prognostic/spatial_discretisation/continuous_galerkin/&
                                &stabilisation/streamline_upwind", &
                                & nu_bar_scheme, nu_bar_scale)
      else if(have_option(trim(density%option_path) // &
                          "/prognostic/spatial_discretisation/continuous_galerkin/&
                          &stabilisation/streamline_upwind_petrov_galerkin")) then
        ewrite(2, *) "SUPG stabilisation"
        stabilisation_scheme = STABILISATION_SUPG
        call get_upwind_options(trim(density%option_path) // & 
                                "/prognostic/spatial_discretisation/continuous_galerkin/&
                                &stabilisation/streamline_upwind_petrov_galerkin", &
                                & nu_bar_scheme, nu_bar_scale)
      else
        ewrite(2, *) "No stabilisation"
        stabilisation_scheme = STABILISATION_NONE
      end if

      call get_option(trim(complete_field_path(density%option_path, stat))//&
          &"/temporal_discretisation/theta", theta)
      call get_option("/timestepping/timestep", dt)
      
      test_mesh => pressure%mesh
      field => velocity

      if(present(ct_rhs)) call zero(ct_rhs)

      ! Clear memory of arrays being designed
      call zero(ctp_m)

      allocate(dfield_t(ele_loc(field, 1), ele_ngi(field, 1), field%dim), &
               dtest_t(ele_loc(test_mesh, 1), ele_ngi(test_mesh, 1), field%dim), &
               ddensity_t(ele_loc(density, 1), ele_ngi(density, 1), field%dim), &
               ele_mat(field%dim, ele_loc(test_mesh, 1), ele_loc(field, 1)), &
               detwei(ele_ngi(field, 1)), &
               density_at_quad(ele_ngi(density, 1)), &
               olddensity_at_quad(ele_ngi(density, 1)), &
               nlvelocity_at_quad(nonlinearvelocity%dim, ele_ngi(nonlinearvelocity, 1)), &
               density_grad_at_quad(ele_ngi(density,1), field%dim), &
               j_mat(field%dim, field%dim, ele_ngi(density, 1)))
      
      do ele=1, element_count(test_mesh)

        test_nodes=>ele_nodes(test_mesh, ele)
        field_nodes=>ele_nodes(field, ele)

        test_shape_ptr => ele_shape(test_mesh, ele)
        field_shape=>ele_shape(field, ele)
        density_shape => ele_shape(density, ele)
        
        density_at_quad = ele_val_at_quad(density, ele)
        olddensity_at_quad = ele_val_at_quad(olddensity, ele)
        
        nlvelocity_at_quad = ele_val_at_quad(nonlinearvelocity, ele)
        
        if(any(stabilisation_scheme == (/STABILISATION_STREAMLINE_UPWIND, STABILISATION_SUPG/))) then
          call transform_to_physical(coordinate, ele, test_shape_ptr, dshape = dtest_t, &
                                     detwei = detwei, j = j_mat)
        else
          call transform_to_physical(coordinate, ele, test_shape_ptr, dshape = dtest_t, detwei=detwei)
        end if
        
        if(.not.integrate_by_parts) then
          ! transform the field (velocity) derivatives into physical space
          call transform_to_physical(coordinate, ele, field_shape, dshape=dfield_t)
          
          if(test_shape_ptr==density_shape) then
            ddensity_t = dtest_t
          else
            call transform_to_physical(coordinate, ele, density_shape, dshape = ddensity_t)
          end if
        else
          dfield_t = 0.0
          ddensity_t = 0.0
        end if
        
        select case(stabilisation_scheme)
          case(STABILISATION_SUPG)
            test_shape = make_supg_shape(test_shape_ptr, dtest_t, nlvelocity_at_quad, j_mat, &
              & nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
          case default
            test_shape = test_shape_ptr
            call incref(test_shape)
        end select
        ! Important note: with SUPG the test function derivatives have not been
        ! modified - i.e. dtest_t is currently used everywhere. This is fine for P1,
        ! but is not consistent for P>1.


        if(integrate_by_parts) then
            ! if SUPG is fixed for P>1 then this dtest_t should be updated
            ele_mat = -dshape_shape(dtest_t, field_shape, &
                       detwei*(theta*density_at_quad + (1-theta)*olddensity_at_quad))
        else
            density_grad_at_quad = theta*(ele_grad_at_quad(density, ele, ddensity_t))+&
                                   (1-theta)*(ele_grad_at_quad(olddensity, ele, ddensity_t))
            
            ele_mat = shape_dshape(test_shape, dfield_t, &
                                   detwei*(theta*density_at_quad + (1-theta)*olddensity_at_quad)) + &
                      shape_shape_vector(test_shape, field_shape, detwei, &
                                         transpose(density_grad_at_quad))
        end if
        
        ! Stabilisation does not return the right shape for this operator!
        select case(stabilisation_scheme)
          case(STABILISATION_STREAMLINE_UPWIND)
!             ele_mat = ele_mat + &
!               & element_upwind_stabilisation_div(field_shape, ddensity_t, nlvelocity_at_quad, j_mat, detwei, &
!               & nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
        end select
        
        do dim = 1, field%dim
            call addto(ctp_m, 1, dim, test_nodes, field_nodes, ele_mat(dim,:,:))
        end do
        
        call deallocate(test_shape)
        
      end do

      if(integrate_by_parts) then

        allocate(detwei_bdy(face_ngi(field, 1)), &
                 density_bdy(face_ngi(density, 1)), &
                 normal_bdy(field%dim, face_ngi(field, 1)))
        allocate(field_nodes_bdy(field%mesh%faces%shape%loc))
        allocate(test_nodes_bdy(test_mesh%faces%shape%loc))
        allocate(ele_mat_bdy(field%dim, face_loc(test_mesh, 1), face_loc(field, 1)))

        assert(surface_element_count(test_mesh)==surface_element_count(field))
        allocate(field_bc_type(field%dim, surface_element_count(field)))
        call get_entire_boundary_condition(field, (/ &
          "weakdirichlet ", &
          "no_normal_flow", &
          "internal      ", &
          "free_surface  "/), field_bc, field_bc_type)

        allocate(density_bc_type(surface_element_count(density)))
        call get_entire_boundary_condition(density, (/ &
          "weakdirichlet"/), density_bc, density_bc_type)

        do sele = 1, surface_element_count(test_mesh)

          if(any(field_bc_type(:,sele)==2)&
               .or.any(field_bc_type(:,sele)==3)&
               .or.any(field_bc_type(:,sele)==4)) cycle
          
          test_shape_ptr=>face_shape(test_mesh, sele)
          field_shape=>face_shape(field, sele)

          test_nodes_bdy=face_global_nodes(test_mesh, sele)
          field_nodes_bdy=face_global_nodes(field, sele)
          
          if(density_bc_type(sele)==1) then
            ! not considering time varying bc yet!
            density_bdy = ele_val_at_quad(density_bc, sele)
          else
            density_bdy = theta*face_val_at_quad(density, sele) + &
                                  (1-theta)*face_val_at_quad(olddensity, sele)
          end if

          call transform_facet_to_physical(coordinate, sele, &
              &                          detwei_f=detwei_bdy,&
              &                          normal=normal_bdy) 

          ele_mat_bdy = shape_shape_vector(test_shape_ptr, field_shape, &
                                           detwei_bdy*density_bdy, normal_bdy)

          do dim = 1, field%dim
            if((field_bc_type(dim, sele)==1).and.present(ct_rhs)) then
              call addto(ct_rhs, test_nodes_bdy, &
                          -matmul(ele_mat_bdy(dim,:,:), &
                          ele_val(field_bc, dim, sele)))
            else
              call addto(ctp_m, 1, dim, test_nodes_bdy, field_nodes_bdy, &
                    ele_mat_bdy(dim,:,:))
            end if
          end do

        end do

        call deallocate(field_bc)
        call deallocate(density_bc)
        deallocate(field_bc_type, density_bc_type)
        deallocate(detwei_bdy, normal_bdy, density_bdy)
        deallocate(test_nodes_bdy, field_nodes_bdy)

      end if
      
    end subroutine assemble_1mat_compressible_divergence_matrix_cg

end module divergence_matrix_cg

