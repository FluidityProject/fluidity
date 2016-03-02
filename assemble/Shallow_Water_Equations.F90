!    Copyright (C) 2006 Imperial College London and others.
!
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
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

! This module implements the 2D shallow water equations within the normal
! fluidity binary configured from an ordinary .flml, going through the 
! usual momentum_equation code path. It is
! a different implementation than the one in main/Shallow_Water.F90
! that has its own binary and schema.
module shallow_water_equations
  use fldebug
  use global_parameters, only: OPTION_PATH_LEN
  use spud
  use transform_elements
  use sparse_tools
  use fields
  use fetools
  use state_module
  use boundary_conditions

  implicit none

  private

  public :: assemble_shallow_water_projection,  assemble_swe_divergence_matrix_cg,&
       shallow_water_equations_check_options

  contains

  ! Assemble the shallow water continuity equation by adding the time derivative (\eta^{n+1}-\eta^n)/dt
  ! This is based on assemble_1mat_compressible_projection_cg, but a lot simpler:
  !  The density \rho is here the total water depth (bottom depth+fs elevation). Its time-derivative however 
  !  is the same as the free surface time derivative. Pressure is g*\eta, so drhodp is simply 1/g

  subroutine assemble_shallow_water_projection(state, cmc, rhs, dt, theta_pg, theta_divergence, reassemble_cmc_m)

    ! This only works for single material_phase, so just pass me the first state:
    type(state_type), intent(inout) :: state
    ! the lhs and rhs to add into:
    type(csr_matrix), intent(inout) :: cmc
    type(scalar_field), intent(inout) :: rhs

    real, intent(in) :: dt
    real, intent(in) :: theta_pg, theta_divergence
    logical, intent(in) :: reassemble_cmc_m

    integer, dimension(:), pointer :: test_nodes
    type(element_type), pointer :: test_shape
    real, dimension(:), allocatable :: detwei
    real, dimension(:), allocatable :: delta_p_at_quad
    real :: rho0, g
    integer :: ele

    type(vector_field), pointer :: coordinate
    type(scalar_field), pointer :: pressure, old_pressure
    
    ewrite(1,*) 'Entering assemble_shallow_water_projection'
    
    call zero(rhs)

    ! only do all this if we need to make cmc (otherwise we'd be adding repeatedly)
    if(reassemble_cmc_m) then
      coordinate=> extract_vector_field(state, "Coordinate")
      
      pressure => extract_scalar_field(state, "Pressure")
      old_pressure => extract_scalar_field(state, "OldPressure")

      call get_option("/physical_parameters/gravity/magnitude", g)
      ! this is the density in front of the du/dt time-derivative which at the 
      ! moment is hard-coded to 1.0
      rho0 = 1.0
  
      allocate(detwei(ele_ngi(pressure, 1)), &
               delta_p_at_quad(ele_ngi(pressure, 1)))
      
      do ele=1, element_count(pressure)
      
        test_nodes=>ele_nodes(pressure, ele)
  
        test_shape => ele_shape(pressure, ele)
        
        delta_p_at_quad = ele_val_at_quad(pressure, ele) - ele_val_at_quad(old_pressure, ele)
                          
        call transform_to_physical(coordinate, ele, detwei=detwei)
  
        ! Time derivative: pressure correction \Phi that we are solving for is: \Phi=theta_div*theta_pg*(p^{n+1}-p^*)*dt
        ! Thus time derivative lhs is (\eta^{n+1}-\eta^*)/dt = \Phi/(g*dt*dt/theta_div/theta_pg)
        call addto(cmc, test_nodes, test_nodes, &
           shape_shape(test_shape, test_shape, detwei)/(rho0*g*dt*dt*theta_divergence*theta_pg))
        ! Time derivative rhs: -(\eta^*-\eta^n)/dt = (p^*-p^n)/(g*dt)
        call addto(rhs, test_nodes, &
           -shape_rhs(test_shape, detwei*delta_p_at_quad)/(rho0*g*dt))
        
      end do

      deallocate(detwei, delta_p_at_quad)
  
    end if

  end subroutine assemble_shallow_water_projection

  subroutine assemble_swe_divergence_matrix_cg(ctp_m, state, ct_rhs)

    ! This only works for single material_phase, so just pass me the first state:
    type(state_type), intent(inout) :: state

    ! the divergence matrix and rhs to be assembled
    type(block_csr_matrix), pointer :: ctp_m
    type(scalar_field), intent(inout), optional :: ct_rhs

    ! local
    type(mesh_type), pointer :: test_mesh

    type(vector_field), pointer :: field

    integer, dimension(:), pointer :: test_nodes, field_nodes

    real, dimension(:,:,:), allocatable :: ele_mat, ele_mat_bdy
    type(element_type), pointer :: field_shape, test_shape
    real, dimension(:,:,:), allocatable :: dfield_t, dtest_t, dbottom_t
    real, dimension(:), allocatable :: detwei
    real, dimension(:,:), allocatable :: normal_bdy
    
    real, dimension(:), allocatable :: depth_at_quad
    real, dimension(:,:), allocatable :: depth_grad_at_quad

    integer :: ele, sele, dim, xdim, ngi

    type(vector_field), pointer :: coordinate, velocity
    type(scalar_field), pointer :: pressure, old_pressure, bottom_depth
    real :: theta, dt, g, rho0

    ! integrate by parts
    logical :: integrate_by_parts

    integer, dimension(:,:), allocatable :: field_bc_type
    type(vector_field) :: field_bc

    ewrite(1,*) 'In assemble_swe_divergence_matrix_cg'

    coordinate=> extract_vector_field(state, "Coordinate")
    
    pressure => extract_scalar_field(state, "Pressure")
    old_pressure => extract_scalar_field(state, "OldPressure")
    bottom_depth => extract_scalar_field(state, "BottomDepth")
    
    velocity=>extract_vector_field(state, "Velocity")
    
    ! note that unlike for compressible we adhere to the option under pressure 
    ! (unless DG for which we don't have a choice)
    integrate_by_parts=have_option(trim(pressure%option_path)// &
           &"/prognostic/spatial_discretisation/integrate_continuity_by_parts") &
       &.or. have_option(trim(velocity%option_path)// &
           &"/prognostic/spatial_discretisation/discontinuous_galerkin")

    ewrite(2,*) "SWE divergence is integrated by parts: ", integrate_by_parts

    ! note that this is different then compressible, as we don't have
    ! a prognostic density equivalent, just treating it as a nonlinear relaxation here
    call get_option(trim(velocity%option_path)// &
                       &"/prognostic/temporal_discretisation/relaxation", theta)
    call get_option("/timestepping/timestep", dt)
    call get_option("/physical_parameters/gravity/magnitude", g)
    ! this is the density in front of the du/dt time-derivative which at the 
    ! moment is hard-coded to 1.0
    rho0 = 1.0
    
    test_mesh => pressure%mesh
    ! the field that ctp_m is multiplied with:
    field => velocity

    if(present(ct_rhs)) call zero(ct_rhs)

    ! Clear memory of arrays being designed
    call zero(ctp_m)

    ngi = ele_ngi(test_mesh, 1)
    xdim = coordinate%dim
    allocate(dtest_t(ele_loc(test_mesh, 1), ngi, xdim), &
             dfield_t(ele_loc(field, 1), ngi, xdim), &
             dbottom_t(ele_loc(bottom_depth, 1), ngi, xdim), &
             ele_mat(field%dim, ele_loc(test_mesh, 1), ele_loc(field, 1)), &
             detwei(ngi), &
             depth_at_quad(ngi), &
             depth_grad_at_quad(xdim, ngi))
    
    do ele=1, element_count(test_mesh)

      test_nodes => ele_nodes(test_mesh, ele)
      field_nodes => ele_nodes(field, ele)

      test_shape => ele_shape(test_mesh, ele)
      field_shape => ele_shape(field, ele)
      call transform_to_physical(coordinate, ele, test_shape, dshape = dtest_t, detwei=detwei)

      depth_at_quad = (theta*ele_val_at_quad(pressure, ele) + &
             (1-theta)*ele_val_at_quad(old_pressure, ele))/(rho0*g) + &
             ele_val_at_quad(bottom_depth, ele)
      
      if(integrate_by_parts) then

        ele_mat = -dshape_shape(dtest_t, field_shape, detwei*depth_at_quad)

      else
          
        ! transform the field (velocity) derivatives into physical space
        call transform_to_physical(coordinate, ele, field_shape, dshape=dfield_t)
        if (bottom_depth%mesh%shape==field_shape) then
          dbottom_t = dfield_t
        else if (bottom_depth%mesh%shape==test_shape) then
          dbottom_t = dtest_t
        else
          call transform_to_physical(coordinate, ele, ele_shape(bottom_depth, ele), &
            dshape=dbottom_t)
        end if
        
        assert( test_shape==pressure%mesh%shape )
        depth_grad_at_quad = (theta*ele_grad_at_quad(pressure, ele, dtest_t) + &
             (1-theta)*ele_grad_at_quad(old_pressure, ele, dtest_t))/(rho0*g) + &
             ele_grad_at_quad(bottom_depth, ele, dbottom_t)

        ele_mat = shape_dshape(test_shape, dfield_t, detwei*depth_at_quad) + &
                  shape_shape_vector(test_shape, field_shape, detwei, depth_grad_at_quad)

      end if
      
      do dim = 1, field%dim
        call addto(ctp_m, 1, dim, test_nodes, field_nodes, ele_mat(dim,:,:))
      end do

    end do
    deallocate(dtest_t, dfield_t, dbottom_t, &
             ele_mat, detwei, depth_at_quad, &
             depth_grad_at_quad)

    if(integrate_by_parts) then

      ngi = face_ngi(field,1)
      xdim = coordinate%dim
      allocate(detwei(ngi), &
               depth_at_quad(ngi), &
               normal_bdy(xdim, ngi))
      allocate(ele_mat_bdy(field%dim, face_loc(test_mesh, 1), face_loc(field, 1)))

      assert(surface_element_count(test_mesh)==surface_element_count(field))
      allocate(field_bc_type(field%dim, surface_element_count(field)))
      call get_entire_boundary_condition(field, (/ &
        "weakdirichlet ", &
        "no_normal_flow", &
        "internal      "/), &
        field_bc, field_bc_type)

      do sele = 1, surface_element_count(test_mesh)

        if(any(field_bc_type(:,sele)==2)&
             .or.any(field_bc_type(:,sele)==3)) cycle
        
        test_shape => face_shape(test_mesh, sele)
        field_shape => face_shape(field, sele)

        call transform_facet_to_physical(coordinate, sele, &
            &                          detwei_f=detwei, &
            &                          normal=normal_bdy) 

        depth_at_quad = (theta*face_val_at_quad(pressure, sele) + &
             (1-theta)*face_val_at_quad(old_pressure, sele))/(rho0*g) + &
             face_val_at_quad(bottom_depth, sele)

        ele_mat_bdy = shape_shape_vector(test_shape, field_shape, &
                                         detwei*depth_at_quad, normal_bdy)

        do dim = 1, field%dim
          if((field_bc_type(dim, sele)==1).and.present(ct_rhs)) then
            call addto(ct_rhs, face_global_nodes(test_mesh, sele), &
                        -matmul(ele_mat_bdy(dim,:,:), &
                        ele_val(field_bc, dim, sele)))
          else
            call addto(ctp_m, 1, dim, face_global_nodes(test_mesh, sele), &
                face_global_nodes(field, sele), ele_mat_bdy(dim,:,:))
          end if
        end do

      end do

      call deallocate(field_bc)
      deallocate(field_bc_type, depth_at_quad, detwei, ele_mat_bdy)
      deallocate(normal_bdy)

    end if
    
  end subroutine assemble_swe_divergence_matrix_cg

  subroutine shallow_water_equations_check_options

    character(len=*), dimension(1:3), parameter:: forbidden_bc_types = (/ &
      "free_surface ", "bulk_formulae", &
      "wind_forcing "  /)
    character(len=OPTION_PATH_LEN):: velocity_option_path, pressure_option_path
    real:: beta
    integer:: i, dim

    if (.not. have_option("/material_phase/vector_field::Velocity/prognostic/equation::ShallowWater")) return

    ewrite(2,*) "Checking shallow water options"

    call get_option("/geometry/dimension", dim)
    if (dim==3) then
      FLExit("With equation type ShallowWater you need a 2D mesh and configuration")
    end if
    if (have_option("/geometry/spherical_earth")) then
      FLExit("Equation type ShallowWater is not implemented for spherical_earth")
    end if
    if (have_option("/geometry/ocean_boundaries")) then
      FLExit("Do not specify /geometry/ocean_boundaries with equation type ShallowWater")
    end if
    if (.not. have_option("/physical_parameters/gravity")) then
      ewrite(0,*) "Missing option /physical_parameters/gravity " // &
         & "(you may ignore /physical_parameters/vector_field::GravityDirection)"
      FLExit("With equation type Shallow water you need to specify gravity")
    end if

    ! Options under /material_phase
    if (option_count("/material_phase")/=1) then
      FLExit("Equation type ShallowWater only works with a single material_phase")
    end if

    ! These don't make sense - so let's just forbid them to avoid confusion
    if (have_option("/material_phase/equation_of_state") .or. &
      have_option("/material_phase/scalar_field::Density")) then
      FLExit("With equation type ShallowWater you're not allowed an equation_of_state or Density field")
    end if

    ! Pressure options

    pressure_option_path = "/material_phase/scalar_field::Pressure/prognostic/"
    if (.not. have_option(pressure_option_path)) then
      FLExit("With equation type ShallowWater you need a prognostic Pressure field")
    end if

    if (have_option(trim(pressure_option_path)//"solver/iterative_method::cg")) then
      ewrite(0,*) "For shallow water equations the pressure matrix is asymmetric"
      ewrite(0,*) 'Therefore you should not use "cg" as the linear solver'
      ewrite(0,*) 'Use "gmres" instead'
      FLExit('Cannot use "cg" as linear solver for Pressure with ShallowWater')
    end if

    if (.not. have_option(trim(pressure_option_path)// &
      "/spatial_discretisation/continuous_galerkin")) then
      ! it might also work with dg, but someone should test that - definitely won't work with cv
      FLExit('Equation type ShallowWater only work with a continuous galerkin Pressure')
    end if

    ! Velocity options

    velocity_option_path = "/material_phase/vector_field::Velocity/prognostic/"
    if (.not. have_option(velocity_option_path)) then
      FLExit("With equation type ShallowWater you need a prognostic Velocity field")
    end if

    call get_option(trim(velocity_option_path)//"/spatial_discretisation/conservative_advection", beta)
    if (beta>0.0) then
      ewrite(0,*) "The shallow water equations are implemented in non-conservative form"
      ewrite(0,*) "To be consistent with that the velocity advection term should be in non-conservative form"
      FLExit("Velocity option spatial_discretisation/conservative_advection should be set to 0.0")
    end if

    if (have_option(trim(velocity_option_path)//"/vertical_stabilization")) then
      FLExit("With equation type ShallowWater you cannot use the option vertical_stabilization")
    end if

    ! boundary conditions that don't make sense:
    do i=1, size(forbidden_bc_types)
      if (have_option(trim(velocity_option_path)//"/boundary_conditions/type::"// &
        trim(forbidden_bc_types(i)))) then
        FLExit("Can't have "//trim(forbidden_bc_types(i))//" boundary condition with ShallowWater")
      end if
    end do


  end subroutine shallow_water_equations_check_options

end module shallow_water_equations

