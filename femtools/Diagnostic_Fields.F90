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

module diagnostic_fields
  !!< A module to calculate diagnostic fields.

  use fldebug
  use global_parameters, only:FIELD_NAME_LEN, current_time, OPTION_PATH_LEN
  use futils
  use spud
  use Vector_Tools
  use parallel_tools
  use quicksort
  use sparse_tools
  use CV_Faces
  use fetools
  use unittest_tools
  use fields
  use state_module
  use halos
  use boundary_conditions
  use field_derivatives
  use field_options
  use sparse_matrices_fields
  use fefields, only: compute_lumped_mass, compute_cv_mass
  use MeshDiagnostics
  use CV_Shape_Functions, only: make_cv_element_shape, make_cvbdy_element_shape
  use CVTools
  use cv_options
  use CV_Upwind_Values
  use CV_Face_Values, only: evaluate_face_val, theta_val
  use sparsity_patterns
  use sparsity_patterns_meshes
  use solvers
  use state_fields_module
  use interpolation_module
  use streamfunction

  implicit none

  private

  public :: insert_diagnostic_field, calculate_diagnostic_variable
  public :: calculate_cfl_number, calculate_galerkin_projection
  
  interface calculate_diagnostic_variable
     module procedure calculate_scalar_diagnostic_variable_single_state, &
          & calculate_scalar_diagnostic_variable_multiple_states, &
          & calculate_vector_diagnostic_variable_single_state, &
          & calculate_vector_diagnostic_variable_multiple_states, &
          & calculate_tensor_diagnostic_variable_single_state, &
          & calculate_tensor_diagnostic_variable_multiple_states
  end interface

  interface calculate_absolute_difference
    module procedure calculate_absolute_difference_scalar, calculate_absolute_difference_vector
  end interface

  interface calculate_galerkin_projection
    module procedure calculate_galerkin_projection_scalar, calculate_galerkin_projection_vector, &
                     calculate_galerkin_projection_tensor
  end interface
  

contains
  
  subroutine insert_diagnostic_field(state, d_field_name, &
    & d_field_mesh, d_field_rank, stat)
    !!< Insert a new diagnostic field of specified rank into state

    type(state_type), intent(inout) :: state
    character(len = *), intent(in) :: d_field_name
    type(mesh_type), intent(inout) :: d_field_mesh
    integer, intent(in) :: d_field_rank
    integer, intent(out), optional :: stat

    type(scalar_field), pointer :: s_field
    type(tensor_field), pointer :: t_field
    type(vector_field), pointer :: v_field

    select case(d_field_rank)
      case(0)
        allocate(s_field)
        call allocate(s_field, d_field_mesh, d_field_name)
        call calculate_diagnostic_variable(state, d_field_name, s_field, &
          & stat)
        if(.not. present_and_zero(stat)) then
          call insert(state, s_field, d_field_name)
        end if
        call deallocate(s_field)
        deallocate(s_field)
      case(1)
        allocate(v_field)
        call allocate(v_field, mesh_dim(d_field_mesh), d_field_mesh, &
          & d_field_name)
        call calculate_diagnostic_variable(state, d_field_name, v_field, &
          & stat)
        if(.not. present_and_zero(stat)) then
          call insert(state, v_field, d_field_name)
        end if
        call deallocate(v_field)
        deallocate(v_field)
      case(2)
        allocate(t_field)
        call allocate(t_field, d_field_mesh, d_field_name)
        call calculate_diagnostic_variable(state, d_field_name, t_field, &
          & stat)
        if(.not. present_and_zero(stat)) then
          call insert(state, t_field, d_field_name)
        end if
        call deallocate(t_field)
        deallocate(t_field)
      case default
        if(present(stat)) then
          stat = 1
          return
        else
          FLExit("Invalid diagnostic field rank")
        end if
    end select

  end subroutine insert_diagnostic_field

  subroutine calculate_scalar_diagnostic_variable_single_state(state, d_field_name, d_field, stat, dt, option_path)
    !!< Calculate the specified scalar diagnostic field d_field_name from state
    !!< and return the field in d_field.

    type(state_type), intent(inout) :: state
    character(len = *), intent(in) :: d_field_name
    type(scalar_field), intent(inout) ::d_field
    integer, optional, intent(out) :: stat
    real, intent(in), optional :: dt
    character(len=*), intent(in), optional :: option_path

    select case(d_field_name)

      case("CFLNumber")
        call calculate_cfl_number(state, d_field, dt=dt)

      case("ControlVolumeCFLNumber")
        call calculate_courant_number_cv(state, d_field, dt=dt)

      case("DG_CourantNumber")
        call calculate_courant_number_DG(state, d_field, dt=dt)
        
      case("CVMaterialDensityCFLNumber")
        call calculate_matdens_courant_number_cv(state, d_field, dt=dt)

      case("GridReynoldsNumber")
        call calculate_grid_reynolds_number(state, d_field)

      case("GridPecletNumber")
        call calculate_grid_peclet_number(state, d_field)

      case("HorizontalVelocityDivergence")
        call calculate_horizontal_velocity_divergence(state, d_field, stat)

      case("KineticEnergyDensity")
        call calculate_ke_density(state, d_field, stat)

      case("GravitationalPotentialEnergyDensity")
        call calculate_pe_density(state, d_field, stat)

      case("IsopycnalCoordinate")
        call calculate_isopycnal_coordinate(state, d_field, stat)
      
      case("BackgroundPotentialEnergyDensity")
        call calculate_back_pe_density(state, d_field, stat)
      
      case("HorizontalStreamFunction")
        call calculate_horizontal_streamfunction(state, d_field)
        
      case("StreamFunction")
        call calculate_stream_function_2d(state, d_field, stat)
        
      CASE("MultiplyConnectedStreamFunction")
        call calculate_stream_function_multipath_2d(state, d_field)

      case("Time")
        call set(d_field, current_time)
        
      case("VelocityDivergence")
        call calculate_velocity_divergence(state, d_field, stat)

      case("Speed")
        call calculate_speed(state, d_field, stat)
     
      case("DiffusiveDissipation")
        call calculate_diffusive_dissipation(state, d_field, stat)
      
      case("RichardsonNumber")
        call calculate_richardson_number_new(state, d_field)

      case("AbsoluteDifference")
        call calculate_absolute_difference(state, d_field)

      case("GalerkinProjection")
        call calculate_galerkin_projection(state, d_field)
       
      case("UniversalNumber")
        call calculate_universal_number(d_field)

      case("NodeOwner")
        call calculate_node_owner(d_field)

      case default
        if(present(stat)) then
          stat = 1
          return
        else
          FLExit("Invalid diagnostic scalar field name supplied")
        end if

    end select

  end subroutine calculate_scalar_diagnostic_variable_single_state
  
  subroutine calculate_scalar_diagnostic_variable_multiple_states(state, d_field_name, d_field, stat)
    !!< Calculate the specified scalar diagnostic field d_field_name from
    !!< the supplied states and return the field in d_field.
    
    type(state_type), dimension(:), intent(in) :: state
    character(len = *), intent(in) :: d_field_name
    type(scalar_field), intent(inout) :: d_field
    integer, optional, intent(out) :: stat
    
    select case(d_field_name)

      case default
        if(present(stat)) then
          stat = 1
          return
        else
          FLExit("Invalid diagnostic scalar field name supplied")
        end if

    end select
    
  end subroutine calculate_scalar_diagnostic_variable_multiple_states

  subroutine calculate_vector_diagnostic_variable_single_state(state, d_field_name, &
    & d_field, stat)
    !!< Calculate the specified vector diagnostic field d_field_name from
    !!< state and return the field in d_field.

    type(state_type), intent(inout) :: state
    character(len = *), intent(in) :: d_field_name
    type(vector_field), intent(inout) :: d_field
    integer, optional, intent(out) :: stat

    select case(d_field_name)

      ! Inner element fields

      case("InnerElementFullVelocity")
        call calculate_sgs_full_velocity(state, d_field, stat)

      case("InnerElementFullVorticity")
        call calculate_sgs_full_vorticity(state, d_field, stat)

      case("InnerElementVorticity")
        call calculate_sgs_vorticity(state, d_field, stat)

      case("DgMappedVelocity")
        call calculate_dg_mapped_cg_velocity(state, d_field, stat)

      case("DgMappedVorticity")
        call calculate_dg_mapped_cg_vorticity(state, d_field, stat)

      ! Inner element fields end

      case("LinearMomentum")
        call calculate_linear_momentum(state, d_field)

      case("AbsoluteDifference")
        call calculate_absolute_difference(state, d_field)

      case("BedShearStress")
        call calculate_bed_shear_stress(state, d_field)

      case("MaxBedShearStress")
        call calculate_max_bed_shear_stress(state, d_field)

      case("GalerkinProjection")
        call calculate_galerkin_projection(state, d_field)
        
      case("DiagnosticCoordinate")
        call calculate_diagnostic_coordinate_field(state, d_field)
        
      case default
        if(present(stat)) then
          stat = 1
          return
        else
          FLExit("Invalid diagnostic vector field name supplied")
        end if

    end select

  end subroutine calculate_vector_diagnostic_variable_single_state
  
  subroutine calculate_vector_diagnostic_variable_multiple_states(state, d_field_name, d_field, stat)
    !!< Calculate the specified vector diagnostic field d_field_name from
    !!< the supplied states and return the field in d_field.
    
    type(state_type), dimension(:), intent(in) :: state
    character(len = *), intent(in) :: d_field_name
    type(vector_field), intent(inout) ::d_field
    integer, optional, intent(out) :: stat
    
    select case(d_field_name)

      case default
        if(present(stat)) then
          stat = 1
          return
        else
          FLExit("Invalid diagnostic vector field name supplied")
        end if

    end select
    
  end subroutine calculate_vector_diagnostic_variable_multiple_states

  subroutine calculate_tensor_diagnostic_variable_single_state(state, d_field_name, &
    & d_field, stat)
    !!< Calculate the specified tensor diagnostic field d_field_name from
    !!< state and return the field in d_field.

    type(state_type), intent(in) :: state
    character(len = *), intent(in) :: d_field_name
    type(tensor_field), intent(inout) :: d_field
    integer, optional, intent(out) :: stat

    select case(d_field_name)

      case default
        if(present(stat)) then
          stat = 1
          return
        else
          FLExit("Invalid diagnostic tensor field name supplied")
        end if

    end select

  end subroutine calculate_tensor_diagnostic_variable_single_state
  
  subroutine calculate_tensor_diagnostic_variable_multiple_states(state, d_field_name, d_field, stat)
    !!< Calculate the specified tensor diagnostic field d_field_name from
    !!< the supplied states and return the field in d_field.
    
    type(state_type), dimension(:), intent(in) :: state
    character(len = *), intent(in) :: d_field_name
    type(tensor_field), intent(inout) ::d_field
    integer, optional, intent(out) :: stat
    
    select case(d_field_name)

      case default
        if(present(stat)) then
          stat = 1
          return
        else
          FLExit("Invalid diagnostic tensor field name supplied")
        end if

    end select
    
  end subroutine calculate_tensor_diagnostic_variable_multiple_states

  subroutine calculate_CFL_number(State, CFL, dt)
    !! Calculate the CFL number as a field.
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: cfl
    real, intent(in), optional :: dt

    type(vector_field), pointer :: U, X
    real :: l_dt
    integer :: ele, gi
    ! Transformed quadrature weights.
    real, dimension(ele_ngi(CFL, 1)) :: detwei
    ! Inverse of the local coordinate change matrix.
    real, dimension(mesh_dim(CFL), mesh_dim(CFL), ele_ngi(CFL, 1)) :: invJ
    ! velocity/dx at each quad point.
    real, dimension(mesh_dim(CFL), ele_ngi(CFL, 1)) :: CFL_q
    ! current element global node numbers.
    integer, dimension(:), pointer :: ele_cfl
    ! local cfl matrix on the current element.
    real, dimension(ele_loc(CFL, 1),ele_loc(CFL, 1)) :: CFL_mat
    ! current CFL element shape
    type(element_type), pointer :: CFL_shape

    U=>extract_vector_field(state, "Velocity")
    X=>extract_vector_field(state, "Coordinate")
    
    if(present(dt)) then
      l_dt = dt
    else
      call get_option("/timestepping/timestep",l_dt)
    end if
    assert(allfequals(l_dt))

    call zero(cfl)

    do ele=1, element_count(CFL)
       ele_CFL=>ele_nodes(CFL, ele)
       CFL_shape=>ele_shape(CFL, ele)

       call compute_inverse_jacobian(X, ele, detwei=detwei, invJ=invJ)

       ! Calculate the CFL number at each quadrature point.
       ! The matmul is the transpose of what I originally thought it should
       ! be. I don't understand why it's this way round but the results
       ! appear correct. -dham
       CFL_q=ele_val_at_quad(U, ele)
       do gi=1, size(detwei)
          CFL_q(:,gi)=l_dt*matmul(CFL_q(:,gi), invJ(:,:,gi))
       end do

       ! Project onto the basis functions to recover CFL at each node.
       CFL_mat=matmul(inverse(shape_shape(CFL_shape, CFL_shape, detwei)), &
            shape_shape(CFL_shape, CFL_shape, &
            &             detwei*maxval(abs(CFL_q),1)))

       ! CFL is inherently discontinuous. In the case where a continuous
       ! mesh is provided for CFL, the following takes the safest option
       ! of taking the maximum value at a node.
       CFL%val(ele_CFL)=max(CFL%val(ele_CFL), sum(CFL_mat,2))

    end do
    
    !call halo_max(cfl)

  end subroutine calculate_CFL_number

  subroutine calculate_grid_reynolds_number(State, GRN)
    !! Calculate the grid reynolds number as a field.
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: grn

    type(vector_field), pointer :: U, X
    integer :: ele, gi, stat, a, b
    ! Transformed quadrature weights.
    real, dimension(ele_ngi(GRN, 1)) :: detwei
    ! Inverse of the local coordinate change matrix.
    real, dimension(mesh_dim(GRN), mesh_dim(GRN), ele_ngi(GRN, 1)) :: J
    ! velocity/dx at each quad point.
    real, dimension(mesh_dim(GRN), ele_ngi(GRN, 1)) :: GRN_q
    ! viscosity at each quad point
    real, dimension(mesh_dim(GRN), mesh_dim(GRN), ele_ngi(GRN,1)) :: vis_q
    ! density at each quad point
    real, dimension(ele_ngi(GRN,1)) :: den_q    
    ! current element global node numbers.
    integer, dimension(:), pointer :: ele_grn
    ! local grn matrix on the current element.
    real, dimension(ele_loc(GRN, 1),ele_loc(GRN, 1)) :: GRN_mat
    ! current GRN element shape
    type(element_type), pointer :: GRN_shape
    type(tensor_field), pointer :: viscosity
    type(scalar_field), pointer :: density
    logical :: include_density_field, use_stress_form
    
    U=>extract_vector_field(state, "Velocity")
    X=>extract_vector_field(state, "Coordinate")

    call zero(grn)

    viscosity => extract_tensor_field(state,'Viscosity')
    
    include_density_field = have_option(trim(GRN%option_path)//'/diagnostic/include_density_field')
    
    if (include_density_field) then
       density => extract_scalar_field(state,'Density', stat = stat)
       if (stat /= 0) then
          FLExit('To include the Density field in the Grid Reynolds number calculation Density must exist in the material_phase state')
       end if
    end if

    if (have_option(trim(U%option_path)//&
            &"/prognostic/spatial_discretisation/continuous_galerkin"//&
            &"/stress_terms/stress_form") .or. &
            have_option(trim(U%option_path)//&
            &"/prognostic/spatial_discretisation/continuous_galerkin"//&
            &"/stress_terms/partial_stress_form")) then
       use_stress_form = .true.
    else
       use_stress_form = .false.
    end if
    
    do ele=1, element_count(GRN)
       ele_GRN=>ele_nodes(GRN, ele)
       GRN_shape=>ele_shape(GRN, ele)

       call compute_jacobian(X, ele, J=J, detwei=detwei)

       ! Calculate the GRN number at each quadrature point.
       ! The matmul is as given by dham
       GRN_q=ele_val_at_quad(U, ele)
       vis_q=ele_val_at_quad(viscosity, ele)

       ! for full and partial stress form we need to set the off diagonal terms of the viscosity tensor to zero
       ! to be able to invert it 
       if (use_stress_form) then
          do a=1,size(vis_q,1)
             do b=1,size(vis_q,2)
                if(a.eq.b) cycle
                vis_q(a,b,:) = 0.0
             end do
          end do
       end if

       do gi=1, size(detwei)
          GRN_q(:,gi)=matmul(GRN_q(:,gi), J(:,:,gi))
          GRN_q(:,gi)=matmul(inverse(vis_q(:,:,gi)), GRN_q(:,gi))
       end do
       
       ! include the density field if required also at the quad point
       if (include_density_field) then
          den_q=ele_val_at_quad(density, ele)
          do gi=1,size(detwei)
              GRN_q(:,gi)=den_q(gi)*GRN_q(:,gi)
          end do          
       end if
       
       ! Project onto the basis functions to recover GRN at each node.
       GRN_mat=matmul(inverse(shape_shape(GRN_shape, GRN_shape, detwei)), &
            shape_shape(GRN_shape, GRN_shape, &
            &             detwei*maxval(abs(GRN_q),1)))

       ! GRN is inherently discontinuous. In the case where a continuous
       ! mesh is provided for GRN, the following takes the safest option
       ! of taking the maximum value at a node.
       GRN%val(ele_GRN)=max(GRN%val(ele_GRN), sum(GRN_mat,2))

    end do

  end subroutine calculate_grid_reynolds_number

  subroutine calculate_grid_peclet_number(State, GPN)
    !! Calculate the grid peclet number as a field.
    !! Basically a rehash of the grid reynolds number calculation
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: gpn

    type(vector_field), pointer :: U, X
    integer :: ele, gi, stat
    ! Transformed quadrature weights.
    real, dimension(ele_ngi(GPN, 1)) :: detwei
    ! Inverse of the local coordinate change matrix.
    real, dimension(mesh_dim(GPN), mesh_dim(GPN), ele_ngi(GPN, 1)) :: J
    ! velocity/dx at each quad point.
    real, dimension(mesh_dim(GPN), ele_ngi(GPN, 1)) :: GPN_q
    real, dimension(mesh_dim(GPN), mesh_dim(GPN), ele_ngi(GPN,1)) :: diffus_q
    ! current element global node numbers.
    integer, dimension(:), pointer :: ele_gpn
    ! local grn matrix on the current element.
    real, dimension(ele_loc(GPN, 1),ele_loc(GPN, 1)) :: GPN_mat
    ! current GPN element shape
    type(element_type), pointer :: GPN_shape
    type(tensor_field), pointer :: diffusivity
    character(len=FIELD_NAME_LEN) :: field_name

    U=>extract_vector_field(state, "Velocity")
    X=>extract_vector_field(state, "Coordinate")

    call zero(gpn)

    call get_option(trim(GPN%option_path)//"/diagnostic/field_name", field_name)

    diffusivity => extract_tensor_field(state,trim(field_name)//'Diffusivity&
         &',stat=stat)
    
    if(stat/=0) then

      FLExit("Can't calculate Peclet number, no diffusivity")
    else

        do ele=1, element_count(GPN)
           ele_GPN=>ele_nodes(GPN, ele)
           GPN_shape=>ele_shape(GPN, ele)

           call compute_jacobian(X, ele, J=J, detwei=detwei)

           ! Calculate the GPN number at each quadrature point.
           ! The matmul is as given by dham
           GPN_q=ele_val_at_quad(U, ele)
           diffus_q=ele_val_at_quad(diffusivity, ele)
           do gi=1, size(detwei)
              GPN_q(:,gi)=matmul(GPN_q(:,gi), J(:,:,gi))
              GPN_q(:,gi)=matmul(inverse(diffus_q(:,:,gi)), GPN_q(:,gi))
           end do

           ! Project onto the basis functions to recover GPN at each node.
           GPN_mat=matmul(inverse(shape_shape(GPN_shape, GPN_shape, detwei)), &
            shape_shape(GPN_shape, GPN_shape, &
            &             detwei*maxval(abs(GPN_q),1)))

           ! GRN is inherently discontinuous. In the case where a continuous
           ! mesh is provided for GRN, the following takes the safest option
           ! of taking the maximum value at a node.
           GPN%val(ele_GPN)=max(GPN%val(ele_GPN), sum(GPN_mat,2))

        end do

    end if

  end subroutine calculate_grid_peclet_number

  subroutine calculate_horizontal_velocity_divergence(state, hveld_field, stat)
    !!< Calculate the horizontal velocity divergence field

    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: hveld_field
    integer, intent(out), optional :: stat

    integer :: i
    type(vector_field) :: hvel_field
    type(vector_field), pointer :: g_direction_field, positions, vel_field

    do i = 1, 3
      select case(i)
        case(1)
          g_direction_field => extract_vector_field(state, "GravityDirection", &
            & stat)
        case(2)
          positions => extract_vector_field(state, "Coordinate", stat)
        case(3)
          vel_field => extract_vector_field(state, "Velocity", stat)
        case default
          FLAbort("Invalid loop index")
      end select
      if(present_and_nonzero(stat)) then
        return
      end if
    end do

    call allocate(hvel_field, mesh_dim(vel_field%mesh), vel_field%mesh, &
      & "HorizontalVelocity")
      
    if (continuity(vel_field)<0 .or. &
        element_degree(vel_field,1)/=element_degree(g_direction_field,1)) then
      FLExit("HorizontalVelocityDivergence does not work for discontinuous or higher order fields.")
    end if

    do i = 1, node_count(hvel_field)
      call set(hvel_field, i, node_val(vel_field, i) - &
        & dot_product(node_val(vel_field, i), &
        & node_val(g_direction_field, i)) * node_val(g_direction_field, i))
    end do

    call div(hvel_field, positions, hveld_field)

    call deallocate(hvel_field)

  end subroutine calculate_horizontal_velocity_divergence

  subroutine calculate_ke_density(state, ke_density_field, stat)
    !!< Calculate the kinetic energy density field
    !!< Beware what your "Density" is!

    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: ke_density_field
    integer, intent(out), optional :: stat

    integer :: i
    type(scalar_field), pointer :: rho_field
    type(vector_field), pointer :: vel_field

    do i = 1, 2
      select case(i)
        case(1)
          rho_field => extract_scalar_field(state, "Density", stat)
        case(2)
          vel_field => extract_vector_field(state, "Velocity", stat)
        case default
          FLAbort("Invalid loop index")
      end select
      if(present_and_nonzero(stat)) then
        return
      end if
    end do

    if(present(stat) .and. (.not. rho_field%mesh == ke_density_field%mesh &
      & .or. .not. vel_field%mesh == ke_density_field%mesh)) then
      stat = 1
      return
    else
      assert(rho_field%mesh == ke_density_field%mesh)
      assert(vel_field%mesh == ke_density_field%mesh)
    end if

    call zero(ke_density_field)
    do i = 1, node_count(ke_density_field)
      call set(ke_density_field, i, &
        & 0.5 * node_val(rho_field, i) * norm2(node_val(vel_field, i))**2)
    end do

  end subroutine calculate_ke_density

  subroutine calculate_pe_density(state, pe_density_field, stat)
    !!< Calculate the gravitational potential energy density field
    !!< Currently assumes a constant gravity field
    !!< Beware what your "Density" is!

    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: pe_density_field
    integer, intent(out), optional :: stat

    integer :: i
    real :: g
    real, dimension(mesh_dim(pe_density_field)) :: g_direction, zero_point
    type(scalar_field), pointer :: rho_field
    type(vector_field), pointer :: positions, positions_remap

    do i = 1, 5
      select case(i)
        case(1)
          call get_option("/physical_parameters/gravity/magnitude", g, stat)
        case(2)
          call get_option( &
            &"/physical_parameters/gravity/" // &
            & "vector_field::GravityDirection/prescribed/value[0]/constant", &
            & g_direction, stat) ! assuming only 1 g_direction as this
                                 ! subroutine isn't set up to support a varying
                                 ! gravity direction over the mesh - need to
                                 ! modify this to get it working with
                                 ! region_ids but assuming this isn't a
                                 ! priority as it's assumed constant already
        case(3)
          call get_option( &
            & trim(pe_density_field%option_path) // "/diagnostic/zero_point", &
            & zero_point, stat)
        case(4)
          positions => extract_vector_field(state, "Coordinate", stat)
        case(5)
          rho_field => extract_scalar_field(state, "Density", stat)
        case default
          FLAbort("Invalid loop index")
      end select
      if(present_and_nonzero(stat)) then
        return
      end if
    end do

    if(present(stat) .and. (.not. rho_field%mesh == pe_density_field%mesh)) then
      stat = 1
      return
    else
      assert(rho_field%mesh == pe_density_field%mesh)
    end if

    if(positions%mesh == pe_density_field%mesh) then
      positions_remap => positions
    else
      allocate(positions_remap)
      call allocate(positions_remap, mesh_dim(pe_density_field%mesh), &
        & pe_density_field%mesh, "Coordinate")
      call remap_field(positions, positions_remap)
    end if

    g_direction = g_direction / norm2(g_direction)

    call zero(pe_density_field)
    do i = 1, node_count(pe_density_field)
      call set(pe_density_field, i, node_val(rho_field, i) * (-1.0) * &
        g * dot_product(g_direction, node_val(positions_remap, i) - zero_point))
    end do

    if(.not. positions%mesh == pe_density_field%mesh) then
      call deallocate(positions_remap)
      deallocate(positions_remap)
    end if

  end subroutine calculate_pe_density

  subroutine calculate_isopycnal_coordinate(state, isopycnal_coordinate, stat)
    !!< Calculate the isopycnal coordinate
    !!< You must be using control volumes for temperature
    !!< You need to set up a prescribed diagnostic scalar field called Width
    !!< which describes the width of your domain as a function of height
    !!< Assumes gravity is in y-direction in 2D, z-direction in 3D

    type(state_type), intent(in) :: state    
    type(scalar_field), intent(inout) :: isopycnal_coordinate
    integer, intent(out), optional :: stat

    integer, dimension(:), allocatable :: index, index2
    integer :: i, j, k, lstat
    real :: z_star, volume, volume_k, volume_rho
    type(scalar_field), pointer :: rho_field
    type(scalar_field) :: lumped_mass, lumped_mass_depth
    type(vector_field), pointer :: Xfield
    type(vector_field) :: Xfield_depth
    type(mesh_type), pointer :: mesh
    character(len = FIELD_NAME_LEN) :: fine_mesh_name

    rho_field => extract_scalar_field(state, "Density", lstat)
    
    Xfield => extract_vector_field(state, "Coordinate")
    call get_option(trim(isopycnal_coordinate%option_path)//"/diagnostic/fine_mesh/name",fine_mesh_name)
    mesh => extract_mesh(state,fine_mesh_name)
    Xfield_depth = get_coordinate_field(state, mesh)
    
    allocate(index(1:node_count(rho_field)))    
    allocate(index2(1:node_count(Xfield_depth)))
    
    ! reorder density from smallest to largest
    call qsort(rho_field%val, index)
    
    ! reorder vertical coordinate
    call qsort(Xfield_depth%val(Xfield_depth%dim,:), index2)
    
    call allocate(lumped_mass, rho_field%mesh, name="LumpedMass")
    call allocate(lumped_mass_depth, mesh, name="LumpedMassDepth")
    
    call zero(lumped_mass)
    call zero(lumped_mass_depth)

    call compute_lumped_mass(Xfield, lumped_mass)
    call compute_lumped_mass(Xfield_depth, lumped_mass_depth)
    

    j = 1
    volume = 0.0
    z_star = 0.0
    do i=node_count(rho_field),1,-1
      volume_rho = node_val(lumped_mass,index(i))
      do k = j, node_count(Xfield_depth)
        if (volume > volume_rho) exit
        volume_k = node_val(lumped_mass_depth,index2(k))
        volume = volume + volume_k
        z_star = z_star + node_val(Xfield_depth,Xfield_depth%dim,index2(k))*volume_k
      end do
      call set(isopycnal_coordinate,index(i),z_star/volume)
      j = min(k, node_count(Xfield_depth))
      volume = volume - volume_rho                      ! left over volume from redistribution
      z_star  = node_val(Xfield_depth,Xfield_depth%dim,index2(j))*volume
    end do

    call deallocate(lumped_mass)
    call deallocate(lumped_mass_depth)
    call deallocate(Xfield_depth)
    deallocate(index)
    deallocate(index2)
    
  end subroutine calculate_isopycnal_coordinate

  subroutine calculate_back_pe_density(state, back_pe_density_field, stat)
    !!< Calculate background potential energy density
    !!< You must have isopycnal_coordinate
    !!< You must be using control volumes for temperature

    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: back_pe_density_field
    integer, intent(out), optional :: stat
    integer :: lstat, i

    real :: g
    type(scalar_field), pointer :: pert_rho_field, rho_field
    type(scalar_field), pointer :: isopycnal_coordinate
    
    rho_field => extract_scalar_field(state, "Density", lstat)
    if (lstat /= 0) then
      if (present(stat)) then
        stat = lstat
      else
        FLExit("Need density")
      end if
    end if

    isopycnal_coordinate => extract_scalar_field(state, "IsopycnalCoordinate", lstat)
    if (lstat /= 0) then
      if (present(stat)) then
        stat = lstat
      else
        FLExit("Need isopycnal coordinate")
      end if
    end if
    
    call get_option("/physical_parameters/gravity/magnitude", g, lstat)
    if (lstat /= 0) then
      if (present(stat)) then
        stat = lstat
      else
        FLExit("Need gravity")
      end if
    end if
    
    call zero(back_pe_density_field)
    do i = 1, node_count(back_pe_density_field)
      call set(back_pe_density_field, i, node_val(rho_field, i) * &
        & g * node_val(isopycnal_coordinate,i))
    end do

  end subroutine calculate_back_pe_density
    
  subroutine calculate_horizontal_streamfunction(state, psi)
    !!< Calculate the horizontal stream function psi where:
    !!<   \partial_x \psi = -v
    !!<   \partial_y \psi = u
    
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: psi
    
    integer :: i
    type(csr_matrix) :: matrix
    type(csr_sparsity), pointer :: sparsity
    type(scalar_field) :: rhs
    type(vector_field), pointer :: gravity_direction, positions, velocity
    
    ewrite(1, *) "In calculate_horizontal_streamfunction"
    ewrite(2, *) "Computing horizontal stream function for state " // trim(state%name)

    if(psi%mesh%continuity /= 0) then
      FLExit("HorizontalStreamFunction requires a continuous mesh")
    end if
    if(mesh_dim(psi%mesh) /= 3) then
      FLExit("HorizontalStreamFunction only works in 3D")
    end if

    ! Extract the Coordinate field
    positions => extract_vector_field(state, "Coordinate")
    assert(positions%dim == mesh_dim(psi))
    assert(ele_count(positions) == ele_count(psi))

    ! Extract velocity    
    velocity => extract_vector_field(state, "Velocity")
    assert(velocity%dim == mesh_dim(psi))
    assert(ele_count(velocity) == ele_count(psi))
    ewrite_minmax(velocity)
    
    ! Extract gravity direction
    gravity_direction => extract_vector_field(state, "GravityDirection")
    assert(gravity_direction%dim == mesh_dim(psi))
    assert(ele_count(gravity_direction) == ele_count(psi))
    
    ! Allocate / extract from state
    sparsity => get_csr_sparsity_firstorder(state, psi%mesh, psi%mesh)
    call allocate(matrix, sparsity, name = trim(psi%name) // "Matrix")
    call allocate(rhs, psi%mesh, trim(psi%name) // "Rhs")
    
    ! Assemble
    call zero(matrix)
    call zero(rhs)
    do i = 1, ele_count(rhs)
      call assemble_horizontal_streamfunction_element(i, psi, matrix, rhs, positions, velocity, gravity_direction)
    end do
    ewrite_minmax(rhs)
    
    ! Boundary conditions - apply strong Dirichlet boundary condition of zero
    ! on all surfaces for now
    do i = 1, surface_element_count(rhs)
      call addto_diag(matrix, face_global_nodes(rhs, i), spread(INFINITY, 1, face_loc(rhs, i)))
    end do
    
    ! Solve
    call petsc_solve(psi, matrix, rhs)
    ewrite_minmax(psi)
    
    ! Deallocate
    call deallocate(matrix)
    call deallocate(rhs)
    
    ewrite(1, *) "Exiting calculate_horizontal_streamfunction"
    
  end subroutine calculate_horizontal_streamfunction
  
  subroutine assemble_horizontal_streamfunction_element(ele, psi, matrix, rhs, positions, velocity, gravity_direction)
    integer, intent(in) :: ele
    type(scalar_field), intent(in) :: psi
    type(csr_matrix), intent(inout) :: matrix
    type(scalar_field), intent(inout) :: rhs
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: velocity
    type(vector_field), intent(in) :: gravity_direction
    
    integer :: i, j
    integer, dimension(:), pointer :: element_nodes
    real, dimension(ele_ngi(psi, ele)) :: detwei, vorticity_h_gi
    real, dimension(mesh_dim(psi), ele_ngi(psi, ele)) :: gravity_direction_gi, vorticity_gi
    real, dimension(ele_loc(psi, ele), ele_ngi(psi, ele), mesh_dim(psi)) :: dn_t_h
    real, dimension(ele_loc(velocity, ele), ele_ngi(velocity, ele), mesh_dim(psi)) :: du_t
    type(element_type), pointer ::  psi_shape, velocity_shape
    
    assert(ele_ngi(velocity, ele) == ele_ngi(psi, ele))
    assert(ele_ngi(gravity_direction, ele) == ele_ngi(psi, ele))
    
    psi_shape => ele_shape(psi, ele)
    velocity_shape => ele_shape(velocity, ele)
    
    call transform_to_physical(positions, ele, psi_shape, dshape = dn_t_h, detwei = detwei)
    if(psi_shape == velocity_shape) then
      du_t = dn_t_h
    else
      call transform_to_physical(positions, ele, velocity_shape, dshape = du_t)
    end if
    
    gravity_direction_gi = ele_val_at_quad(gravity_direction, ele)
    
    forall(i = 1:size(dn_t_h, 1), j = 1:size(dn_t_h, 2))
      dn_t_h(i, j, :) = dn_t_h(i, j, :) - dot_product(dn_t_h(i, j, :), gravity_direction_gi(:, j)) * gravity_direction_gi(:, j)
    end forall
    
    vorticity_gi = ele_curl_at_quad(velocity, ele, du_t)
    do i = 1, size(vorticity_h_gi)
      vorticity_h_gi(i) = -dot_product(vorticity_gi(:, i), gravity_direction_gi(:, i))
    end do
    
    element_nodes => ele_nodes(psi, ele)
    
    call addto(matrix, element_nodes, element_nodes, dshape_dot_dshape(dn_t_h, dn_t_h, detwei))
    call addto(rhs, element_nodes, shape_rhs(psi_shape, detwei * vorticity_h_gi))
  
  end subroutine assemble_horizontal_streamfunction_element
  
  subroutine calculate_sgs_full_velocity(state, sgs_full, stat)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: sgs_full
    integer, intent(out), optional :: stat

    integer :: i
    type(vector_field), pointer :: v_field, sgs_component

    do i = 1, 2
      select case(i)
        case(1)
          sgs_component => extract_vector_field(state, "VelocityInnerElement", stat)
        case(2)
          v_field => extract_vector_field(state, "Velocity", stat)
        case default
          FLAbort("Invalid loop index")
      end select
      if(present_and_nonzero(stat)) then
        return
      end if
    end do

    call remap_field(from_field = v_field , to_field = sgs_full)
    call addto(sgs_full,sgs_component)

  end subroutine calculate_sgs_full_velocity

  subroutine calculate_sgs_full_vorticity(state, vort_field, stat)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: vort_field
    integer, intent(out), optional :: stat

    integer :: i
    type(vector_field), pointer :: positions, v_field

    do i = 1, 2
      select case(i)
        case(1)
          positions => extract_vector_field(state, "Coordinate", stat)
        case(2)
          v_field => extract_vector_field(state, "InnerElementFullVelocity", stat)
        case default
          FLAbort("Invalid loop index")
      end select
      if(present_and_nonzero(stat)) then
        return
      end if
    end do

    call curl(v_field, positions, curl_field = vort_field)

  end subroutine calculate_sgs_full_vorticity

  subroutine calculate_sgs_vorticity(state, vort_field, stat)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: vort_field
    integer, intent(out), optional :: stat

    integer :: i
    type(vector_field), pointer :: positions, v_field

    do i = 1, 2
      select case(i)
        case(1)
          positions => extract_vector_field(state, "Coordinate", stat)
        case(2)
          v_field => extract_vector_field(state, "VelocityInnerElement", stat)
        case default
          FLAbort("Invalid loop index")
      end select
      if(present_and_nonzero(stat)) then
        return
      end if
    end do

    call curl(v_field, positions, curl_field = vort_field)

  end subroutine calculate_sgs_vorticity

  subroutine calculate_dg_mapped_cg_velocity(state, dg_field, stat)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: dg_field
    integer, intent(out), optional :: stat

    integer :: i
    type(vector_field), pointer :: positions, v_field

    do i = 1, 2
      select case(i)
        case(1)
          positions => extract_vector_field(state, "Coordinate", stat)
        case(2)
          v_field => extract_vector_field(state, "Velocity", stat)
        case default
          FLAbort("Invalid loop index")
      end select
      if(present_and_nonzero(stat)) then
        return
      end if
    end do

    call remap_field(from_field = v_field , to_field = dg_field)

  end subroutine calculate_dg_mapped_cg_velocity

  subroutine calculate_dg_mapped_cg_vorticity(state, vort_field, stat)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: vort_field
    integer, intent(out), optional :: stat

    integer :: i
    type(vector_field), pointer :: positions, v_field

    do i = 1, 2
      select case(i)
        case(1)
          positions => extract_vector_field(state, "Coordinate", stat)
        case(2)
          v_field => extract_vector_field(state, "DgMappedVelocity", stat)
        case default
          FLAbort("Invalid loop index")
      end select
      if(present_and_nonzero(stat)) then
        return
      end if
    end do

    call curl(v_field, positions, curl_field = vort_field)

  end subroutine calculate_dg_mapped_cg_vorticity

  subroutine calculate_speed(state, speed_field, stat)
    !!< Calculate the speed field

    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: speed_field
    integer, intent(out), optional :: stat
    integer :: lstat

    integer :: i
    type(vector_field), pointer :: vel_field
    
    vel_field => extract_vector_field(state, "Velocity", lstat)
    if (lstat /= 0) then
      if (present(stat)) then
        stat = lstat
      else
        FLExit("Need Velocity")
      end if
    end if
    
    call zero(speed_field)
    do i = 1, node_count(speed_field)
      call set(speed_field, i, norm2(node_val(vel_field, i)))
    end do

  end subroutine calculate_speed

  subroutine calculate_diffusive_dissipation(state, diffusive_dissipation_field, stat)
    !!< Calculate -2*kappa*g*drho_dy 
    !!< this can be used to calculate diffusive dissipation
    !!< 2D at the moment
    !!< it probably should be generalised 
    !!< currently assumes a constant gravity field
    !!< also assumes an isotropic diffusivity
    
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: diffusive_dissipation_field
    integer, intent(out), optional :: stat
    
    integer :: i
    real :: g
    type(scalar_field), pointer :: rho_field
    type(vector_field), pointer :: positions
    type(scalar_field), dimension(1) :: drho_dy
    
    rho_field => extract_scalar_field(state, "Density", stat)
    positions => extract_vector_field(state, "Coordinate", stat)
      
    if(present_and_nonzero(stat)) then
      return
    end if
    
    call get_option("/physical_parameters/gravity/magnitude", g, stat)

    call allocate(drho_dy(1), rho_field%mesh, "DRhoDy")
    
    call differentiate_field(rho_field, &
     & positions, (/.false., .true./), drho_dy)
    
    call zero(diffusive_dissipation_field)
    do i = 1, node_count(diffusive_dissipation_field)
      call set(diffusive_dissipation_field, i, -g * node_val(drho_dy(1),i))
    end do  
    
    call deallocate(drho_dy(1))

  end subroutine calculate_diffusive_dissipation

  subroutine calculate_richardson_number_old(state, richardson_number_field)
    !!< Calculate the Richardson number field
    !!< Defined in Turner, Buoyancy Effects in Fluids, p.12 as
    !!< Ri = \frac{N^2}{(\frac{\partial u}{\partial z})^2+(\frac{\partial v}{\partial z})^2}
    !!< with N^2 = -\frac{g}{\rho_0}\frac{\partial \rho}{\partial z}
    !!< currently assumes a constant gravity field
    
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: richardson_number_field
    
#ifdef DDEBUG
    integer :: stat
    type(vector_field), pointer :: gravity_direction
#endif
    integer :: i
    real :: g
    type(scalar_field), pointer :: pert_rho_field
    type(vector_field), pointer :: positions, vel_field
    type(scalar_field), dimension(1) :: du_dz
    type(scalar_field), dimension(1) :: dv_dz
    type(scalar_field), dimension(1) :: drho_dz
    
    positions => extract_vector_field(state, "Coordinate")
    vel_field => extract_vector_field(state, "Velocity")
    pert_rho_field => extract_scalar_field(state, "PerturbationDensity")
        
    call get_option("/physical_parameters/gravity/magnitude", g)
#ifdef DDEBUG
    gravity_direction => extract_vector_field(state, "GravityDirection", stat)
    if(stat == 0) then
      select case(gravity_direction%dim)
        case(3)
          assert(all(abs(gravity_direction%val(1,:)) < epsilon(0.0)))
          assert(all(abs(gravity_direction%val(2,:)) < epsilon(0.0)))
          assert(all(abs(gravity_direction%val(3,:) + 1.0) < epsilon(0.0)))
        case(2)
          assert(all(abs(gravity_direction%val(1,:)) < epsilon(0.0)))
          assert(all(abs(gravity_direction%val(2,:) + 1.0) < epsilon(0.0)))
        case default
          FLAbort("Invalid dimension")
      end select
    end if
#endif
    
    assert(positions%mesh == richardson_number_field%mesh)
    assert(vel_field%mesh == richardson_number_field%mesh)
    assert(pert_rho_field%mesh == richardson_number_field%mesh)
    
    select case(positions%dim)
      case(3)
        call allocate(du_dz(1), vel_field%mesh, "DuDz")
        call allocate(dv_dz(1), vel_field%mesh, "DvDz")
        call allocate(drho_dz(1), pert_rho_field%mesh, "DRhoDz")
        
        call differentiate_field(extract_scalar_field(vel_field, 1), &
         & positions, (/.false., .false., .true./), du_dz)
        call differentiate_field(extract_scalar_field(vel_field, 2), &
         & positions, (/.false., .false., .true./), dv_dz)
        call differentiate_field(pert_rho_field, &
         & positions, (/.false., .false., .true./), drho_dz)
        
        call zero(richardson_number_field)
        do i = 1, node_count(richardson_number_field)
          call set(richardson_number_field, i, -g * node_val(drho_dz(1),i) /  &
           & (node_val(du_dz(1),i) ** 2 + node_val(dv_dz(1),i) ** 2))
        end do  
        
        call deallocate(du_dz(1))
        call deallocate(dv_dz(1))
        call deallocate(drho_dz(1))
      case(2)
        ! Actually dy
        call allocate(du_dz(1), vel_field%mesh, "DuDz")
        call allocate(drho_dz(1), pert_rho_field%mesh, "DRhoDz")
        call differentiate_field(extract_scalar_field(vel_field, 1), &
         & positions, (/.false., .true./), du_dz)
        call differentiate_field(pert_rho_field, &
         & positions, (/.false., .true./), drho_dz)
       
        call zero(richardson_number_field)
        do i = 1, node_count(richardson_number_field)
          call set(richardson_number_field, i, -g * node_val(drho_dz(1),i) /  &
           & (node_val(du_dz(1),i) ** 2))
        end do  
        
        call deallocate(du_dz(1))
        call deallocate(drho_dz(1))
      case default
        FLAbort("Invalid dimension")
    end select
    
  end subroutine calculate_richardson_number_old
  
  subroutine calculate_richardson_number_new(state, ri)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: ri
    
    integer :: i
    real :: g
    type(scalar_field), pointer :: masslump, perturbation_density
    type(vector_field), pointer :: gravity_direction, positions, velocity
    
    ewrite(1, *) "In calculate_richardson_number"
    ewrite(2, *) "Computing shear Richardson number for state " // trim(state%name)

    if(ri%mesh%continuity /= 0) then
      FLExit("RichardsonNumber requires a continuous mesh")
    end if

    ! Extract the Coordinate field
    positions => extract_vector_field(state, "Coordinate")
    assert(positions%dim == mesh_dim(ri))
    assert(ele_count(positions) == ele_count(ri))

    ! Extract velocity    
    velocity => extract_vector_field(state, "Velocity")
    assert(velocity%dim == mesh_dim(ri))
    assert(ele_count(velocity) == ele_count(ri))
    ewrite_minmax(velocity)
    
    ! Extract gravity
    gravity_direction => extract_vector_field(state, "GravityDirection")
    assert(gravity_direction%dim == mesh_dim(gravity_direction))
    assert(ele_count(gravity_direction) == ele_count(gravity_direction))
    call get_option("/physical_parameters/gravity/magnitude", g)
    
    ! Extract perturbation density
    perturbation_density => extract_scalar_field(state, "PerturbationDensity")
    ewrite_minmax(perturbation_density)
        
    ! Assemble
    call zero(ri)
    do i = 1, ele_count(ri)
      call assemble_richardson_number_element(i, ri, positions, velocity, g, perturbation_density)
    end do
    ewrite_minmax(ri)
    
    masslump => get_lumped_mass(state, ri%mesh)
    
    ! Solve (somewhat trivial)
    ri%val = ri%val / masslump%val
    ewrite_minmax(ri)
  
  end subroutine calculate_richardson_number_new
  
  subroutine assemble_richardson_number_element(ele, ri, positions, velocity, g, perturbation_density)
    integer, intent(in) :: ele
    type(scalar_field), intent(inout) :: ri
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: velocity
    real, intent(in) :: g
    type(scalar_field), intent(in) :: perturbation_density
    
    integer :: dim, i
    integer, dimension(:), pointer :: element_nodes
    real, dimension(ele_ngi(ri, ele)) :: denomenator_gi, detwei
    real, dimension(mesh_dim(ri), ele_ngi(ri, ele)) :: grad_theta_gi
    real, dimension(ele_loc(ri, ele), ele_ngi(ri, ele), mesh_dim(ri)) :: dn_t
    real, dimension(ele_loc(perturbation_density, ele), ele_ngi(perturbation_density, ele), mesh_dim(ri)) :: dtheta_t
    real, dimension(ele_loc(velocity, ele), ele_ngi(velocity, ele), mesh_dim(ri)) :: du_t
    type(element_type), pointer :: theta_shape, ri_shape, velocity_shape
    
    assert(ele_ngi(velocity, ele) == ele_ngi(ri, ele))
    assert(ele_ngi(perturbation_density, ele) == ele_ngi(ri, ele))
    
    dim = mesh_dim(ri)
    
    ri_shape => ele_shape(ri, ele)
    velocity_shape => ele_shape(velocity, ele)
    theta_shape => ele_shape(perturbation_density, ele)
    
    call transform_to_physical(positions, ele, ri_shape, &
      & dshape = dn_t, detwei = detwei)
    if(ri_shape == velocity_shape) then
      du_t = dn_t
    else
      call transform_to_physical(positions, ele, velocity_shape, dshape = du_t)
    end if
    if(theta_shape == velocity_shape) then
      dtheta_t = dn_t
    else
      call transform_to_physical(positions, ele, theta_shape, dshape = dtheta_t)
    end if
    
    grad_theta_gi = ele_grad_at_quad(perturbation_density, ele, dtheta_t)

    denomenator_gi = 0.0
    do i = 1, dim - 1
      denomenator_gi = denomenator_gi + (matmul(ele_val(velocity, i, ele), du_t(:, :, dim)) ** 2)
    end do
      
    element_nodes => ele_nodes(ri, ele)
    
    call addto(ri, element_nodes, &
      ! Note well: if \frac{d\theta}{dz} and
      ! ((\frac{du}{dz})^2 + \frac{dv}{dz})^2 are zero, then we pick up a value
      ! of NaN here
      & shape_rhs(ri_shape, -detwei * g * grad_theta_gi(dim, :) / denomenator_gi) &
      & )
    
  end subroutine assemble_richardson_number_element

  subroutine calculate_stream_function_2d(state, streamfunc, stat)
    !!< Calculate the stream function for a 
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: streamfunc
    integer, intent(out), optional :: stat
    
    integer :: i, lstat, ele
    type(vector_field), pointer :: X, U
    type(csr_sparsity) :: psi_sparsity
    type(csr_matrix) :: psi_mat
    type(scalar_field) :: rhs

    do i = 1, 2
       select case(i)
       case(1)
          X => extract_vector_field(state, "Coordinate", stat)
       case(2)
          U => extract_vector_field(state, "Velocity", stat)
       case default
          FLAbort("Invalid loop index")
       end select
       if(present_and_nonzero(stat)) then
          return
       end if
    end do
    
    assert(X%dim==2)
    ! No discontinuous stream functions.
    assert(continuity(streamfunc)>=0)
    
    psi_sparsity = extract_csr_sparsity(state, &
               &                      "StreamFunctionSparsity", lstat)
    if (lstat/=0) then
       psi_sparsity = make_sparsity(streamfunc%mesh, streamfunc%mesh, &
            "StreamFunctionSparsity")
    else
       call incref(psi_sparsity)
    end if
    
    call allocate(psi_mat, psi_sparsity, name="StreamFunctionMatrix")

    call zero(psi_mat)
    call allocate(rhs, streamfunc%mesh, "StreamFunctionRHS")
    call zero(rhs)

    do ele=1, element_count(streamfunc)
       
       call calculate_streamfunc_ele(psi_mat, rhs, ele, X, U)

    end do
    
    call petsc_solve(streamfunc, psi_mat, rhs)

    call deallocate(rhs)
    call deallocate(psi_mat)
    call deallocate(psi_sparsity)

  contains
    
    subroutine calculate_streamfunc_ele(psi_mat, rhs, ele, X, U)
      type(csr_matrix), intent(inout) :: psi_mat
      type(scalar_field), intent(inout) :: rhs
      type(vector_field), intent(in) :: X,U
      integer, intent(in) :: ele

      ! Transformed gradient function for velocity.
      real, dimension(ele_loc(U, ele), ele_ngi(U, ele), mesh_dim(U)) :: du_t
      ! Ditto for the stream function, psi
      real, dimension(ele_loc(rhs, ele), ele_ngi(rhs, ele), mesh_dim(rhs))&
           & :: dpsi_t 

      ! Local vorticity_matrix
      real, dimension(2, ele_loc(rhs, ele), ele_loc(U, ele)) ::&
           & lvorticity_mat
      ! Local vorticity
      real, dimension(ele_loc(rhs, ele)) :: lvorticity

      ! Variable transform times quadrature weights.
      real, dimension(ele_ngi(U,ele)) :: detwei
      
      type(element_type), pointer :: U_shape, psi_shape
      integer, dimension(:), pointer :: psi_ele, neigh
      integer :: i, ni, face

      U_shape=> ele_shape(U, ele)
      psi_shape=> ele_shape(rhs, ele)
      psi_ele=>ele_nodes(rhs, ele)
      
      ! Transform U derivatives and weights into physical space.
      call transform_to_physical(X, ele, U_shape, dshape=du_t, detwei=detwei)
      ! Ditto psi.
      call transform_to_physical(X, ele, psi_shape, dshape=dpsi_t)

      call addto(psi_mat, psi_ele, psi_ele, &
           -dshape_dot_dshape(dpsi_t, dpsi_t, detwei))
      
      lvorticity_mat=shape_curl_shape_2d(psi_shape, du_t, detwei)
      
      lvorticity=0.0
      do i=1,2
         lvorticity=lvorticity &
              +matmul(lvorticity_mat(i,:,:), ele_val(U, i, ele))
      end do
      
      call addto(rhs, psi_ele, lvorticity)
      
      neigh=>ele_neigh(U, ele)
      
      neighbourloop: do ni=1,size(neigh)
         ! Find boundaries.
         if (neigh(ni)<=0) then

            face=ele_face(rhs, ele, neigh(ni))

            ! Strong dirichlet condition (currently the only thing supported)
            call addto_diag(psi_mat, &
                 face_global_nodes(rhs, face), &
                 spread(INFINITY, 1, face_loc(rhs,face)))
         end if
      end do neighbourloop

    end subroutine calculate_streamfunc_ele

  end subroutine calculate_stream_function_2d

  subroutine calculate_velocity_divergence(state, div_u, stat)
    !!< Calculate div u for diagnosing continuity issues.
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: div_u
    integer, intent(out), optional :: stat

    integer :: i
    type(vector_field), pointer :: positions, v_field

    do i = 1, 2
      select case(i)
        case(1)
          positions => extract_vector_field(state, "Coordinate", stat)
        case(2)
          v_field => extract_vector_field(state, "Velocity", stat)
        case default
          FLAbort("Invalid loop index")
      end select
      if(present_and_nonzero(stat)) then
        return
      end if
    end do

    call div(v_field, positions, div_u)

  end subroutine calculate_velocity_divergence

  subroutine calculate_courant_number_dg(state, courant, dt)
    !!< Calculate courant number for DG velocity fields
    !!< == positive fluxes of unit function into element
    !!< *dt/volume of element

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: courant
    real, intent(in), optional :: dt
    !
    type(vector_field), pointer :: u, x
    real :: l_dt
    integer :: ele, stat

    u=>extract_vector_field(state, "NonlinearVelocity",stat)
    if(stat.ne.0) then    
       u=>extract_vector_field(state, "Velocity",stat)
       if(stat.ne.0) then
          FLExit('Missing velocity field!')
       end if
    end if
    x=>extract_vector_field(state, "Coordinate")

    if(present(dt)) then
       l_dt = dt
    else
       call get_option("/timestepping/timestep",l_dt)
    end if    
    
    call zero(courant)
    
    do ele = 1, element_count(courant)
       call calculate_courant_number_dg_ele(courant,x,u,ele,l_dt)
    end do
        
    ! the courant values at the edge of the halo are going to be incorrect
    ! this matters when computing the max courant number
    call halo_update(courant)

  end subroutine calculate_courant_number_dg

  subroutine calculate_courant_number_dg_ele(courant,x,u,ele,dt)
    type(vector_field), intent(in) :: x, u
    type(scalar_field), intent(inout) :: courant
    real, intent(in) :: dt
    integer, intent(in) :: ele
    !
    real :: Vol
    real :: Flux
    integer :: ni, ele_2, face, face_2
    integer, dimension(:), pointer :: neigh
    real, dimension(ele_ngi(u,ele)) :: detwei
    real, dimension(face_ngi(u,1)) :: detwei_f
    real, dimension(U%dim, face_ngi(U, 1)) :: normal, U_f_quad
    real, dimension(face_ngi(U,1)) :: flux_quad
    integer, dimension(:), pointer :: u_ele
    real :: val
    real, dimension(ele_loc(u,ele)) :: Vals
    !
    !Get element volume
    call transform_to_physical(X, ele, detwei=detwei)
    Vol = sum(detwei)
    
    !Get fluxes
    Flux = 0.0
    neigh=>ele_neigh(U, ele)
    do ni = 1, size(neigh)
       ele_2=neigh(ni)
       face=ele_face(U, ele, ele_2)
       if(ele_2<0.0) then
          face_2 = face
       else
          face_2=ele_face(U, ele_2, ele)
       end if
       
       U_f_quad =0.5*(face_val_at_quad(U, face)&
            &     +face_val_at_quad(U, face_2))

       call transform_facet_to_physical(X, face, &
            &                          detwei_f=detwei_f,&
            &                          normal=normal) 

       Flux_quad = -sum(U_f_quad*normal,1)
       Flux_quad = max(Flux_quad,0.0)

       Flux = Flux + sum(Flux_quad*detwei_f)
    end do

    u_ele => ele_nodes(U,ele)

    Val = Flux/Vol*dt
    Vals = Val
    call set(Courant,U_ele,Vals)

  end subroutine calculate_courant_number_dg_ele

   subroutine calculate_courant_number_cv(state, courant, dt)

      type(state_type), intent(inout) :: state
      type(scalar_field), intent(inout) :: courant
      real, intent(in), optional :: dt

      type(vector_field), pointer :: u, ug, x
      real :: l_dt
      integer :: ele, iloc, oloc, gi, ggi, sele, face, ni, face_2

      integer :: quaddegree
      type(cv_faces_type) :: cvfaces
      type(element_type) :: x_cvshape, x_cvbdyshape
      type(element_type) :: u_cvshape, u_cvbdyshape
      type(element_type) :: ug_cvshape, ug_cvbdyshape
      type(scalar_field), pointer :: cvmass
      real, dimension(:,:), allocatable :: x_ele, x_ele_bdy
      real, dimension(:,:), allocatable :: x_f, u_f, ug_f, u_bdy_f, ug_bdy_f
      real, dimension(:,:), allocatable :: normal, normal_bdy
      real, dimension(:), allocatable :: detwei, detwei_bdy
      real, dimension(:), allocatable :: normgi
      integer, dimension(:), pointer :: nodes, x_nodes, neigh
      integer, dimension(:), allocatable :: nodes_bdy
      real :: udotn, income
      ! logical array indicating if a face has already been visited by the opposing node
      logical, dimension(:), allocatable :: notvisited
            
      integer, dimension(:), allocatable :: courant_bc_type
      type(scalar_field) :: courant_bc

      type(vector_field) :: x_courant ! coordinates on courant mesh
      
      logical :: move_mesh

      ewrite(1,*) 'in calculate_courant_number_cv'
      
      move_mesh = have_option("/mesh_adaptivity/mesh_movement")

      udotn=0.0 ! to stop valgrind complaining about it being unitialised

      u=>extract_vector_field(state, "NonlinearVelocity")
      x=>extract_vector_field(state, "Coordinate")

      if(move_mesh) ug=>extract_vector_field(state, "GridVelocity")

      if(present(dt)) then
        l_dt = dt
      else
        call get_option("/timestepping/timestep",l_dt)
      end if

      call zero(courant)

      x_courant=get_coordinate_field(state, courant%mesh)
      
      ! determine the cv mass matrix to use for the length scale
      cvmass => get_cv_mass(state, courant%mesh)
      ewrite_minmax(cvmass)    
      
      if(courant%mesh%shape%degree /= 0) then

        call get_option("/geometry/quadrature/controlvolume_surface_degree", &
                      quaddegree, default=1)

        cvfaces=find_cv_faces(vertices=ele_vertices(courant,1), &
                              dimension=mesh_dim(courant), &
                              polydegree=courant%mesh%shape%degree, &
                              quaddegree=quaddegree)

        u_cvshape=make_cv_element_shape(cvfaces, u%mesh%shape)
        x_cvshape=make_cv_element_shape(cvfaces, x%mesh%shape)

        allocate(x_ele(x%dim,ele_loc(x,1)), &
                x_f(x%dim, x_cvshape%ngi), &
                u_f(u%dim, u_cvshape%ngi), &
                detwei(x_cvshape%ngi), &
                normal(x%dim, x_cvshape%ngi), &
                normgi(x%dim))
        allocate(notvisited(x_cvshape%ngi))

        if(move_mesh) then
          ug_cvshape=make_cv_element_shape(cvfaces, ug%mesh%shape)
          allocate(ug_f(ug%dim, ug_cvshape%ngi))
        end if

        do ele=1, element_count(courant)
          x_ele=ele_val(x, ele)
          x_f=ele_val_at_quad(x, ele, x_cvshape)
          u_f=ele_val_at_quad(u, ele, u_cvshape)
          if(move_mesh) ug_f = ele_val_at_quad(ug, ele, ug_cvshape)
          nodes=>ele_nodes(courant, ele)
          x_nodes=>ele_nodes(x_courant, ele)

          call transform_cvsurf_to_physical(x_ele, x_cvshape, &
                                            detwei, normal, cvfaces)

          notvisited=.true.

          do iloc = 1, courant%mesh%shape%loc

            do face = 1, cvfaces%faces

              if(cvfaces%neiloc(iloc, face) /= 0) then
                oloc = cvfaces%neiloc(iloc, face)

                do gi = 1, cvfaces%shape%ngi

                  ggi = (face-1)*cvfaces%shape%ngi + gi

                  ! have we been here before?
                  if(notvisited(ggi)) then
                    notvisited(ggi)=.false.

                    normgi=orientate_cvsurf_normgi(node_val(x_courant, x_nodes(iloc)),x_f(:,ggi),normal(:,ggi))

                    if(move_mesh) then
                      udotn=dot_product((u_f(:,ggi)-ug_f(:,ggi)), normgi)
                    else
                      udotn=dot_product(u_f(:,ggi), normgi)
                    end if

                    if(udotn>0.0) then
                      income=0.0
                    else
                      income=1.0
                    end if

                    call addto(courant, nodes(iloc), abs(udotn)*(1.-income)*detwei(ggi))
                    call addto(courant, nodes(oloc), abs(udotn)*income*detwei(ggi)) ! notvisited

                  end if ! notvisited

                end do

              end if
            end do
          end do
        end do

        u_cvbdyshape=make_cvbdy_element_shape(cvfaces, u%mesh%faces%shape)
        x_cvbdyshape=make_cvbdy_element_shape(cvfaces, x%mesh%faces%shape)

        allocate(x_ele_bdy(x%dim,face_loc(x,1)), &
                u_bdy_f(u%dim, u_cvbdyshape%ngi), &
                detwei_bdy(x_cvbdyshape%ngi), &
                normal_bdy(x%dim, x_cvbdyshape%ngi))
        allocate(nodes_bdy(face_loc(courant, 1)))
        allocate(courant_bc_type(surface_element_count(courant)))
        
        if(move_mesh) then
          ug_cvbdyshape=make_cvbdy_element_shape(cvfaces, ug%mesh%faces%shape)
          allocate(ug_bdy_f(ug%dim, ug_cvbdyshape%ngi))
        end if

        ! get the fields over the surface containing the bcs
        call get_entire_boundary_condition(courant, (/"internal"/), courant_bc, courant_bc_type)
        
        do sele=1,surface_element_count(courant)
        
          if(courant_bc_type(sele)==1) cycle

          ele = face_ele(x, sele)
          x_ele = ele_val(x, ele)
          x_ele_bdy = face_val(x, sele)
          nodes_bdy=face_global_nodes(courant, sele)

          call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, &
                                x_cvbdyshape, normal_bdy, detwei_bdy)

          u_bdy_f=face_val_at_quad(u, sele, u_cvbdyshape)
          if(move_mesh) ug_bdy_f=face_val_at_quad(ug, sele, ug_cvbdyshape)

          do iloc = 1, courant%mesh%faces%shape%loc

            do face = 1, cvfaces%sfaces

              if(cvfaces%sneiloc(iloc,face)/=0) then

                do gi = 1, cvfaces%shape%ngi

                  ggi = (face-1)*cvfaces%shape%ngi + gi
                  
                    if(move_mesh) then
                      udotn=dot_product((u_bdy_f(:,ggi)-ug_bdy_f(:,ggi)), normal_bdy(:,ggi))
                    else
                      udotn=dot_product(u_bdy_f(:,ggi), normal_bdy(:,ggi))
                    end if

                    if(udotn>0.0) then
                      income=0.0
                    else
                      income=1.0
                    end if

                    call addto(courant, nodes_bdy(iloc), abs(udotn)*(1.0-income)*detwei_bdy(ggi))

                end do

              end if

            end do

          end do

        end do

        deallocate(x_ele, x_f, u_f, detwei, normal, normgi)
        deallocate(x_ele_bdy, u_bdy_f, detwei_bdy, normal_bdy)
        deallocate(nodes_bdy)
        deallocate(notvisited)
        call deallocate(x_cvbdyshape)
        call deallocate(u_cvbdyshape)
        call deallocate(x_cvshape)
        call deallocate(u_cvshape)
        call deallocate(cvfaces)
        call deallocate(courant_bc)
        
        if(move_mesh) then
          call deallocate(ug_cvshape)
          call deallocate(ug_cvbdyshape)
          deallocate(ug_f, ug_bdy_f)
        end if

      else

        allocate(detwei(face_ngi(courant, 1)), &
                 u_f(u%dim, face_ngi(u, 1)), &
                 normal(x%dim, face_ngi(courant, 1)))

        do ele = 1, element_count(courant)

          nodes=>ele_nodes(courant, ele)
          assert(size(nodes)==1)

          neigh=>ele_neigh(courant, ele)

          do ni= 1, size(neigh)

            face = ele_face(courant, ele, neigh(ni))

            if(neigh(ni)>0) then
              ! internal face
              face_2=ele_face(courant, neigh(ni), ele)
            else
              ! external face
              face_2 = face
            end if

            call transform_facet_to_physical(x, face, detwei_f=detwei, normal=normal)

            ! if velocity is dg then use a trapezoidal rule (otherwise this will
            ! all cancel out to give the face value)
            u_f = 0.5*(face_val_at_quad(u, face) + face_val_at_quad(u, face_2))
            if(move_mesh) then
              u_f = u_f - face_val_at_quad(ug, face)
            end if

            call addto(courant, nodes(1), &
                 sum(sum(u_f*normal,1)*merge(1.0,0.0,.not.(sum(u_f*normal,1)<0.0))*detwei))

          end do


        end do

      end if

      courant%val = courant%val*l_dt/cvmass%val

      call deallocate(x_courant)
      
      call halo_update(courant)

   end subroutine calculate_courant_number_cv

   subroutine calculate_matdens_courant_number_cv(state, courant, dt)

      type(state_type), intent(inout) :: state
      type(scalar_field), intent(inout) :: courant
      real, intent(in), optional :: dt

      type(vector_field), pointer :: u, ug, x
      type(vector_field) :: x_courant
      type(scalar_field), pointer :: matdens, oldmatdens
      real :: l_dt
      integer :: ele, iloc, oloc, gi, ggi, sele, face, stat

      integer :: quaddegree
      type(cv_faces_type) :: cvfaces
      type(element_type) :: u_cvshape, u_cvbdyshape
      type(element_type) :: ug_cvshape, ug_cvbdyshape
      type(element_type) :: t_cvshape, t_cvbdyshape
      type(element_type) :: x_cvshape, x_cvbdyshape
      type(scalar_field) :: cvmass
      real, dimension(:,:), allocatable :: x_ele, x_ele_bdy
      real, dimension(:,:), allocatable :: x_f, u_f, u_bdy_f, ug_f, ug_bdy_f
      real, dimension(:,:), allocatable :: normal, normal_bdy
      real, dimension(:), allocatable :: detwei, detwei_bdy
      real, dimension(:), allocatable :: normgi
      integer, dimension(:), pointer :: nodes, x_nodes
      integer, dimension(:), allocatable :: nodes_bdy
      real :: udotn, income
      logical :: inflow

      ! options for the material density field
      type(cv_options_type) :: matdens_options
      real :: matdens_face_val, oldmatdens_face_val, matdens_theta_val
      ! type of courant number we want to use (if any)
      character(len=FIELD_NAME_LEN) :: cfl_type
      ! the courant number field
      type(scalar_field) :: cfl_no
      type(csr_sparsity) :: mesh_sparsity
      type(csr_matrix) :: matdens_upwind, oldmatdens_upwind
      real, dimension(:), allocatable :: matdens_ele, oldmatdens_ele, &
                                         cfl_ele
      integer, dimension(:), allocatable :: matdens_bc_type
      type(scalar_field) :: matdens_bc
      real, dimension(:), allocatable :: matdens_ele_bdy, oldmatdens_ele_bdy, &
              ghost_matdens_ele_bdy, ghost_oldmatdens_ele_bdy
      type(state_type), dimension(1) :: state_array

      ! logical array indicating if a face has already been visited by the opposing node
      logical, dimension(:), allocatable :: notvisited
      
      logical :: move_mesh
      
      move_mesh = have_option("/mesh_adaptivity/mesh_movement")

      udotn=0.0 ! to stop valgrind complaining about it being unitialised

      u=>extract_vector_field(state, "IteratedVelocity")
      x=>extract_vector_field(state, "Coordinate")
      x_courant = get_coordinate_field(state, courant%mesh)

      if(move_mesh) ug=>extract_vector_field(state, "GridVelocity")

      ! extract the material density and get all the relevent options:
      matdens=>extract_scalar_field(state, "MaterialDensity")
      oldmatdens=>extract_scalar_field(state, "OldMaterialDensity", stat)
      if(stat/=0) then
        oldmatdens=>matdens
      end if
      matdens_options=get_cv_options(matdens%option_path, matdens%mesh%shape%numbering%family, mesh_dim(matdens))

      ! hmmm, do we need a cfl no.?
      ! the density discretisation might need it
      ! clearly this can't use a CVMaterialDensityCFLNumber as that's what we're
      ! currently trying to find out!
      cfl_type=""
      call allocate(cfl_no, matdens%mesh, "CourantNumber")

      call get_option(trim(complete_cv_field_path(matdens%option_path))//&
                             "/face_value[0]/courant_number[0]/name", &
                            cfl_type, stat)
      if(stat==0) then
         select case(trim(cfl_type))
         case("CVMaterialDenstiyCFLNumber")
            FLAbort("You can't use the field you're in the process of creating!")
         case default
            ! otherwise we want to calculate a node centred field of the cfl number
            call calculate_diagnostic_variable(state, trim(cfl_type), cfl_no)
         end select
      else
         ! if we don't need a cfl number then just set it all to 1
         call set(cfl_no, 1.0)
      end if

      ! allocate upwind value matrices
      mesh_sparsity=make_sparsity(matdens%mesh, matdens%mesh, "MaterialDensitySparsity")

      call allocate(matdens_upwind, mesh_sparsity, name="MaterialDensityUpwindValues")
      call allocate(oldmatdens_upwind, mesh_sparsity, name="OldMaterialDensityUpwindValues")

      ! does the density field need upwind values?
      if(need_upwind_values(matdens_options)) then

        state_array(1) = state  ! a hack to let find_upwind_values accept a single state

        call find_upwind_values(state_array, x_courant, matdens, matdens_upwind, &
                                oldmatdens, oldmatdens_upwind &
                                )

      else

        call zero(matdens_upwind)
        call zero(oldmatdens_upwind)

      end if

      if(present(dt)) then
        l_dt = dt
      else
        call get_option("/timestepping/timestep",l_dt)
      end if

      call get_option("/geometry/quadrature/controlvolume_surface_degree", &
                     quaddegree, default=1)

      call zero(courant)

      cvfaces=find_cv_faces(vertices=ele_vertices(courant,1), &
                            dimension=mesh_dim(courant), &
                            polydegree=courant%mesh%shape%degree, &
                            quaddegree=quaddegree)

      u_cvshape=make_cv_element_shape(cvfaces, u%mesh%shape)
      x_cvshape=make_cv_element_shape(cvfaces, x%mesh%shape)
      t_cvshape=make_cv_element_shape(cvfaces, matdens%mesh%shape)

      allocate(x_ele(x%dim, ele_loc(x,1)), &
               x_f(x%dim, x_cvshape%ngi), &
               u_f(u%dim, u_cvshape%ngi), &
               detwei(x_cvshape%ngi), &
               normal(x%dim, x_cvshape%ngi), &
               normgi(x%dim))
      allocate(cfl_ele(ele_loc(cfl_no, 1)), &
               matdens_ele(ele_loc(matdens, 1)), &
               oldmatdens_ele(ele_loc(matdens, 1)))
      allocate(notvisited(x_cvshape%ngi))
      
      if(move_mesh) then
        ug_cvshape = make_cv_element_shape(cvfaces, ug%mesh%shape)
        allocate(ug_f(ug%dim, ug_cvshape%ngi))
      end if

      call allocate(cvmass, courant%mesh, "CV mass")
      call compute_cv_mass(x, cvmass)
      cvmass%val = cvmass%val*(matdens_options%theta*matdens%val+(1.0-matdens_options%theta)*oldmatdens%val)

      do ele=1, element_count(courant)
        x_ele=ele_val(x, ele)
        x_f=ele_val_at_quad(x, ele, x_cvshape)
        u_f=ele_val_at_quad(u, ele, u_cvshape)
        if(move_mesh) ug_f = ele_val_at_quad(ug, ele, ug_cvshape)
        nodes=>ele_nodes(courant, ele)
        x_nodes=>ele_nodes(x_courant, ele)

        call transform_cvsurf_to_physical(x_ele, x_cvshape, &
                                          detwei, normal, cvfaces)

        matdens_ele = ele_val(matdens, ele)
        oldmatdens_ele = ele_val(oldmatdens, ele)

        cfl_ele = ele_val(cfl_no, ele)

        notvisited=.true.

        do iloc = 1, courant%mesh%shape%loc

          do face = 1, cvfaces%faces

            if(cvfaces%neiloc(iloc, face) /= 0) then
              oloc = cvfaces%neiloc(iloc, face)

              do gi = 1, cvfaces%shape%ngi

                ggi = (face-1)*cvfaces%shape%ngi + gi

                ! have we been here before?
                if(notvisited(ggi)) then
                  notvisited(ggi)=.false.

                  normgi=orientate_cvsurf_normgi(node_val(x_courant, x_nodes(iloc)),x_f(:,ggi),normal(:,ggi))

                  if(move_mesh) then
                    udotn=dot_product((u_f(:,ggi)-ug_f(:,ggi)), normgi)
                  else
                    udotn=dot_product(u_f(:,ggi), normgi)
                  end if

                  inflow = (udotn<=0.0)

                  income = merge(1.0,0.0,inflow)

                  select case (matdens%field_type)
                  case(FIELD_TYPE_CONSTANT)

                      matdens_face_val = matdens_ele(iloc)
                      oldmatdens_face_val = oldmatdens_ele(iloc)

                  case default

                      call evaluate_face_val(matdens_face_val, oldmatdens_face_val, &
                                             iloc, oloc, ggi, nodes, &
                                             t_cvshape,&
                                             matdens_ele, oldmatdens_ele, &
                                             matdens_upwind, oldmatdens_upwind, &
                                             inflow, cfl_ele, &
                                             matdens_options)

                  end select

                  matdens_theta_val=theta_val(iloc, oloc, &
                                       matdens_face_val, &
                                       oldmatdens_face_val, &
                                       matdens_options%theta, l_dt, udotn, &
                                       x_ele, matdens_options%limit_theta, &
                                       matdens_ele, oldmatdens_ele)


                  call addto(courant, nodes(iloc), abs(udotn)*(1.-income)*detwei(ggi)*matdens_theta_val)
                  call addto(courant, nodes(oloc), abs(udotn)*income*detwei(ggi)*matdens_theta_val) ! notvisited

                end if ! notvisited

              end do

            end if
          end do
        end do
      end do

      u_cvbdyshape=make_cvbdy_element_shape(cvfaces, u%mesh%faces%shape)
      x_cvbdyshape=make_cvbdy_element_shape(cvfaces, x%mesh%faces%shape)
      t_cvbdyshape=make_cvbdy_element_shape(cvfaces, matdens%mesh%faces%shape)

      allocate(x_ele_bdy(x%dim,face_loc(x,1)), &
              u_bdy_f(u%dim, u_cvbdyshape%ngi), &
              detwei_bdy(x_cvbdyshape%ngi), &
              normal_bdy(x%dim, x_cvbdyshape%ngi), &
              matdens_ele_bdy(face_loc(matdens,1)), &
              oldmatdens_ele_bdy(face_loc(oldmatdens,1)), &
              ghost_matdens_ele_bdy(face_loc(matdens,1)), &
              ghost_oldmatdens_ele_bdy(face_loc(oldmatdens,1)))
      allocate(matdens_bc_type(surface_element_count(matdens)), &
               nodes_bdy(face_loc(courant,1)))
               
      if(move_mesh) then
        ug_cvbdyshape=make_cvbdy_element_shape(cvfaces, ug%mesh%faces%shape)
        allocate(ug_bdy_f(ug%dim, ug_cvbdyshape%ngi))
      end if

      ! get the fields over the surface containing the bcs
      call get_entire_boundary_condition(matdens, (/"weakdirichlet", &
                                                    "internal     "/), matdens_bc, matdens_bc_type)

      do sele=1,surface_element_count(courant)
      
        if(matdens_bc_type(sele)==2) cycle

        ele = face_ele(x, sele)
        x_ele = ele_val(x, ele)
        x_ele_bdy = face_val(x, sele)
        nodes_bdy=face_global_nodes(courant, sele)

        call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, &
                              x_cvbdyshape, normal_bdy, detwei_bdy)

        u_bdy_f=face_val_at_quad(u, sele, u_cvbdyshape)
        if(move_mesh) ug_bdy_f=face_val_at_quad(ug, sele, ug_cvbdyshape)

        ! deal with bcs for tdensity
        if(matdens_bc_type(sele)==1) then
          ghost_matdens_ele_bdy=ele_val(matdens_bc, sele)
        else
          ghost_matdens_ele_bdy=face_val(matdens, sele)
        end if

        if(matdens_bc_type(sele)==1) then
          ghost_oldmatdens_ele_bdy=ele_val(matdens_bc, sele) ! not considering time varying bcs yet
        else
          ghost_oldmatdens_ele_bdy=face_val(oldmatdens, sele)
        end if

        matdens_ele_bdy=face_val(matdens, sele)
        oldmatdens_ele_bdy=face_val(oldmatdens, sele)

        do iloc = 1, courant%mesh%faces%shape%loc

          do face = 1, cvfaces%sfaces

            if(cvfaces%sneiloc(iloc,face)/=0) then

              do gi = 1, cvfaces%shape%ngi

                ggi = (face-1)*cvfaces%shape%ngi + gi

                  if(move_mesh) then
                    udotn=dot_product((u_bdy_f(:,ggi)-ug_bdy_f(:,ggi)), normal_bdy(:,ggi))
                  else
                    udotn=dot_product(u_bdy_f(:,ggi), normal_bdy(:,ggi))
                  end if

                  if(udotn>0.0) then
                     income=0.0
                  else
                     income=1.0
                  end if

                  ! as we're on the boundary it's not possible to use high order methods so just
                  ! default to the pivotted solution method (first order upwinding)
                  ! if the flow is incoming then use the bc ghost values
                  ! if the flow is outgoing then use the surface nodes value

                  matdens_face_val = income*ghost_matdens_ele_bdy(iloc) + (1.-income)*matdens_ele_bdy(iloc)
                  oldmatdens_face_val = income*ghost_oldmatdens_ele_bdy(iloc) + (1.-income)*oldmatdens_ele_bdy(iloc)

                  matdens_theta_val = matdens_options%theta*matdens_face_val + (1.-matdens_options%theta)*oldmatdens_face_val

                  call addto(courant, nodes_bdy(iloc), abs(udotn)*(1.0-income)*detwei_bdy(ggi)*matdens_theta_val)

               end do

            end if

          end do

        end do

      end do

      courant%val = courant%val*l_dt/cvmass%val

      deallocate(x_ele, x_f, u_f, detwei, normal, normgi)
      deallocate(x_ele_bdy, u_bdy_f, detwei_bdy, normal_bdy)
      deallocate(nodes_bdy)
      deallocate(notvisited)
      call deallocate(u_cvbdyshape)
      call deallocate(x_cvbdyshape)
      call deallocate(t_cvbdyshape)
      call deallocate(u_cvshape)
      call deallocate(x_cvshape)
      call deallocate(t_cvshape)
      call deallocate(cvfaces)
      call deallocate(cvmass)
      call deallocate(cfl_no)
      call deallocate(matdens_upwind)
      call deallocate(oldmatdens_upwind)
      call deallocate(matdens_bc)
      call deallocate(mesh_sparsity)
      
      if(move_mesh) then
        call deallocate(ug_cvshape)
        call deallocate(ug_cvbdyshape)
        deallocate(ug_f, ug_bdy_f)
      end if

   end subroutine calculate_matdens_courant_number_cv

   subroutine calculate_linear_momentum(state, momentum)

      type(state_type), intent(in) :: state
      type(vector_field), intent(inout) :: momentum

      type(scalar_field), pointer :: density
      type(vector_field), pointer :: velocity
      
      type(scalar_field), pointer :: tmpdensity

      density => extract_scalar_field(state, "Density")
      if(.not.density%mesh==momentum%mesh) then
        allocate(tmpdensity)
        call allocate(tmpdensity, momentum%mesh, "TmpDensity")
        call remap_field(density, tmpdensity)
      else
        tmpdensity => density
      end if
      
      velocity => extract_vector_field(state, "Velocity")
      
      call remap_field(velocity, momentum)
      call scale(momentum, tmpdensity)
      
      if(.not.density%mesh==momentum%mesh) then
        call deallocate(tmpdensity)
        deallocate(tmpdensity)
      end if

   end subroutine calculate_linear_momentum

   subroutine calculate_absolute_difference_scalar(state, difference)

      type(state_type), intent(in) :: state
      type(scalar_field), intent(inout) :: difference

      type(scalar_field), pointer :: field_a, field_b
      
      type(scalar_field), pointer :: l_field_a, l_field_b
      logical :: remap_field_a, remap_field_b
      type(vector_field), pointer :: field_coordinate, difference_coordinate
      
      character(len=FIELD_NAME_LEN) :: field_name_a, field_name_b

      real :: av_diff, max_a, max_b, min_a, min_b, av_a, av_b

      call get_option(trim(difference%option_path)//"/diagnostic/field_name_a", field_name_a)
      call get_option(trim(difference%option_path)//"/diagnostic/field_name_b", field_name_b)

      field_a => extract_scalar_field(state, trim(field_name_a))
      if(mesh_compatible(field_a%mesh, difference%mesh)) then
        remap_field_a=.false.
        l_field_a=>field_a
      else
        remap_field_a=.true.
        allocate(l_field_a)
        call allocate(l_field_a, difference%mesh, trim(field_a%name))
        call zero(l_field_a)
        
        field_coordinate => get_external_coordinate_field(state, field_a%mesh)
        difference_coordinate => get_external_coordinate_field(state, difference%mesh)
        
        call linear_interpolation(field_a, field_coordinate, l_field_a, difference_coordinate)
      end if
      
      field_b => extract_scalar_field(state, trim(field_name_b))
      if(mesh_compatible(field_b%mesh, difference%mesh)) then
        remap_field_b=.false.
        l_field_b=>field_b
      else
        remap_field_b=.true.
        allocate(l_field_b)
        call allocate(l_field_b, difference%mesh, trim(field_b%name))
        call zero(l_field_b)
        
        field_coordinate => get_external_coordinate_field(state, field_b%mesh)
        difference_coordinate => get_external_coordinate_field(state, difference%mesh)
        
        call linear_interpolation(field_b, field_coordinate, l_field_b, difference_coordinate)
      end if

      if (have_option(trim(difference%option_path)//"/diagnostic/relative_to_average")) then
        call field_stats(l_field_a, max=max_a)
        call field_stats(l_field_a, min=min_a)
        av_a = (max_a+min_a)/2.0
        call field_stats(l_field_b, max=max_b)
        call field_stats(l_field_b, min=min_b)
        av_b = (max_b+min_b)/2.0
        av_diff = av_a-av_b
      else
        av_diff = 0.0
      end if
      
      call set(difference, l_field_a)
      call addto(difference, l_field_b, -1.0)
      call addto(difference, -av_diff)

      difference%val = abs(difference%val)

      if (have_option(trim(difference%option_path)//"/diagnostic/ignore_boundaries")) then
        if (associated(difference%mesh%faces%surface_node_list)) then
          if (size(difference%mesh%faces%surface_node_list)>0) then
            call set(difference, difference%mesh%faces%surface_node_list, &
                    spread(0.0, 1, size(difference%mesh%faces%surface_node_list)))
          end if
        end if
      end if
      
      if(remap_field_a) then
        call deallocate(l_field_a)
        deallocate(l_field_a)
      end if
      if(remap_field_b) then
        call deallocate(l_field_b)
        deallocate(l_field_b)
      end if

   end subroutine calculate_absolute_difference_scalar

   subroutine calculate_absolute_difference_vector(state, difference)

      type(state_type), intent(in) :: state
      type(vector_field), intent(inout) :: difference

      type(vector_field), pointer :: field_a, field_b
      
      type(vector_field), pointer :: l_field_a, l_field_b
      logical :: remap_field_a, remap_field_b
      type(vector_field), pointer :: field_coordinate, difference_coordinate
      
      character(len=FIELD_NAME_LEN) :: field_name_a, field_name_b
      integer :: i

      real, dimension(difference%dim) :: av_diff
      real :: max_a, max_b, min_a, min_b, av_a, av_b
      type(scalar_field) :: field_comp

      call get_option(trim(difference%option_path)//"/diagnostic/field_name_a", field_name_a)
      call get_option(trim(difference%option_path)//"/diagnostic/field_name_b", field_name_b)

      field_a => extract_vector_field(state, trim(field_name_a))
      if(mesh_compatible(field_a%mesh, difference%mesh)) then
        remap_field_a=.false.
        l_field_a=>field_a
      else
        remap_field_a=.true.
        allocate(l_field_a)
        call allocate(l_field_a, field_a%dim, difference%mesh, trim(field_a%name))
        call zero(l_field_a)
        
        field_coordinate => get_external_coordinate_field(state, field_a%mesh)
        difference_coordinate => get_external_coordinate_field(state, difference%mesh)
        
        call linear_interpolation(field_a, field_coordinate, l_field_a, difference_coordinate)
      end if
      
      field_b => extract_vector_field(state, trim(field_name_b))
      if(mesh_compatible(field_b%mesh, difference%mesh)) then
        remap_field_b=.false.
        l_field_b=>field_b
      else
        remap_field_b=.true.
        allocate(l_field_b)
        call allocate(l_field_b, field_b%dim, difference%mesh, trim(field_b%name))
        call zero(l_field_b)
        
        field_coordinate => get_external_coordinate_field(state, field_b%mesh)
        difference_coordinate => get_external_coordinate_field(state, difference%mesh)
        
        call linear_interpolation(field_b, field_coordinate, l_field_b, difference_coordinate)
      end if

      av_diff = 0.0
      if (have_option(trim(difference%option_path)//"/diagnostic/relative_to_average")) then
        do i = 1, difference%dim
          field_comp = extract_scalar_field(l_field_a, i)
          call field_stats(field_comp, max=max_a)
          call field_stats(field_comp, min=min_a)
          av_a = (max_a+min_a)/2.0
          field_comp = extract_scalar_field(l_field_b, i)
          call field_stats(field_comp, max=max_b)
          call field_stats(field_comp, min=min_b)
          av_b = (max_b+min_b)/2.0
          av_diff(i) = av_a-av_b
        end do
      end if
      
      call set(difference, l_field_a)
      call addto(difference, l_field_b, -1.0)
      call addto(difference, -av_diff)

      do i = 1, difference%dim
        difference%val(i,:) = abs(difference%val(i,:))
      end do
      
      if(remap_field_a) then
        call deallocate(l_field_a)
        deallocate(l_field_a)
      end if
      if(remap_field_b) then
        call deallocate(l_field_b)
        deallocate(l_field_b)
      end if

   end subroutine calculate_absolute_difference_vector

   subroutine calculate_bed_shear_stress(state, bed_shear_stress)
!
      type(state_type), intent(inout) :: state
      type(vector_field), intent(inout) :: bed_shear_stress
      type(scalar_field) :: masslump
      type(vector_field), pointer :: U, X
      type(tensor_field), pointer :: visc
      integer, dimension(:), allocatable :: faceglobalnodes
      integer :: i,j,snloc,ele,sele,globnod,face,node,stat
      real :: speed,density,drag_coefficient

      !! for DG
      !! Field that holds the gradient of velocity in boundary elements
      type(tensor_field), target :: dummy_visc
      integer :: grad_u_stat, visc_stat
      !! surface mesh, element and node list
      type(mesh_type), pointer :: surface_mesh
      integer, dimension(:), allocatable :: surface_element_list
      !! surface fields
      type(vector_field) :: bed_shear_stress_surface
      type(tensor_field) :: grad_U, visc_surface, grad_u_surface
      
      ewrite(2,*) 'in calculate bed_shear_stress'

      if (have_option(trim(bed_shear_stress%option_path)//"/prescribed")) then
        ewrite(2,*) 'prescribed bed_shear_stress - not calculating'
        return
      end if

      ! assumes constant density
      call get_option(trim(bed_shear_stress%option_path)//"/diagnostic/density", density)

      ! calculate using drag coefficient
      if (have_option(trim(bed_shear_stress%option_path)//&
           &"/diagnostic/calculation_method/drag_coefficient")) then

         call zero(bed_shear_stress) 

         call get_option(trim(bed_shear_stress%option_path)//&
              & "/diagnostic/calculation_method/drag_coefficient",&
              & drag_coefficient)

         U => extract_vector_field(state, "Velocity")
         snloc = face_loc(U, 1)
         allocate( faceglobalnodes(1:snloc) )
         do sele=1,surface_element_count(U)
            ele = face_ele(U, sele)
            faceglobalnodes = face_global_nodes(U, sele)
            do j = 1,snloc
               globnod = faceglobalnodes(j)
               speed = norm2(node_val(U, globnod))
               call set(bed_shear_stress, globnod, density*drag_coefficient*speed * node_val(U, globnod))
            end do
         end do
         deallocate( faceglobalnodes )
         
      ! calculate using velocity gradient
      else if (have_option(trim(bed_shear_stress%option_path)//&
           &"/diagnostic/calculation_method/velocity_gradient")) then

         call zero(bed_shear_stress) 

         visc => extract_tensor_field(state, "Viscosity", visc_stat)
         if (visc_stat /= 0.0) then
            ewrite(0,*) 'Warning: No viscosity specified - assumed to be 1.0 for bed shear calculation'
            call allocate(dummy_visc, bed_shear_stress%mesh, 'dummy_visc')
            call zero(dummy_visc)
            do i = 1, dummy_visc%dim(1)
               call set(dummy_visc, i, i, 1.0)
            end do
            visc => dummy_visc            
         end if
         U    => extract_vector_field(state, "Velocity")
         X    => extract_vector_field(state, "Coordinate")

         ! Check velociy and bed shear stress meshes are consistent
         if (continuity(bed_shear_stress) /= continuity(U) .or. &
             element_degree(bed_shear_stress, 1) /= element_degree(U, 1)) then
            FLAbort('Bed shear stress and velocity mesh must have the same continuity and degree')
         end if  
         
         if(continuity(bed_shear_stress)>=0) then       
            ! We need to calculate a global lumped mass over the surface elements
            call allocate(masslump, bed_shear_stress%mesh, 'Masslump')
            call zero(masslump)

            do face = 1, surface_element_count(bed_shear_stress)
               call calculate_bed_shear_stress_ele_cg(bed_shear_stress, masslump, face, X, U,&
                    & visc, density)
            end do

            where (masslump%val/=0.0)
               masslump%val=1./masslump%val
            end where
            call scale(bed_shear_stress, masslump)
            call deallocate(masslump)
         else
            ! We do DG differently. First the gradient of the velocity field is calculated   
            ! using the field_derivatives code.
            ! Then we use this field to determine the bed shear stress using:
            ! N_i N_j tau_b = N_i nu grad_u . |n|

            ! create a field to store the gradient on
            call allocate(grad_U, bed_shear_stress%mesh, 'grad_U')
            call zero(grad_U)
            ! calculate gradient of velocity
            call grad(U, X, grad_U)

            allocate(surface_element_list(surface_element_count(bed_shear_stress)))
            ! generate list of surface elements
            do i=1, surface_element_count(bed_shear_stress)
               surface_element_list(i)=i
            end do

            ! create surface field
            surface_mesh => get_dg_surface_mesh(bed_shear_stress%mesh)
            call allocate(bed_shear_stress_surface, bed_shear_stress%dim, surface_mesh)

            ! remap required fields to the boundary surfaces
            call allocate(grad_u_surface, surface_mesh, dim=grad_u%dim)
            call remap_field_to_surface(grad_u, grad_u_surface, surface_element_list)
            call allocate(visc_surface, surface_mesh, dim=visc%dim)
            call remap_field_to_surface(visc, visc_surface, surface_element_list)

            ! calculate bed shear stress
            do face = 1, ele_count(bed_shear_stress_surface)
               call calculate_bed_shear_stress_ele_dg(bed_shear_stress_surface, face, X, grad_u_surface,&
                    & visc_surface, density)

               ! copy values to volume field - can be done element by element as the surface is generated
               ! as we are in DG
               call set(bed_shear_stress, &
                    face_global_nodes(bed_shear_stress, face), &
                    ele_val(bed_shear_stress_surface, face))
            end do

            call deallocate(bed_shear_stress_surface)
            call deallocate(grad_u)
            call deallocate(grad_u_surface)
            call deallocate(visc_surface)
         end if

         if (visc_stat /= 0) then
            call deallocate(dummy_visc)
         end if
      else
         FLAbort('Unknown bed shear stress calculation method')
      end if

   end subroutine calculate_bed_shear_stress

   subroutine calculate_bed_shear_stress_ele_cg(bed_shear_stress, masslump, face, X, U, visc&
        &, density)

     type(vector_field), intent(inout) :: bed_shear_stress
     type(scalar_field), intent(inout) :: masslump
     type(vector_field), intent(in), pointer :: X, U
     type(tensor_field), intent(in), pointer :: visc
     integer, intent(in) :: face
     real, intent(in) :: density

     integer :: i, j, i_gi, ele, dim
     type(element_type), pointer :: f_shape, shape, X_f_shape, X_shape
     real, dimension(face_ngi(X, face)) :: detwei
     real, dimension(X%dim, face_ngi(X, face)) :: normal, normal_shear_at_quad, X_ele
     real, dimension(X%dim) :: abs_normal
     real, dimension(ele_loc(X, face_ele(X, face)), face_ngi(X, face), X%dim) :: ele_dshape_at_face_quad
     real, dimension(X%dim, X%dim, face_ngi(X, face)) :: grad_U_at_quad, visc_at_quad, shear_at_quad  
     real, dimension(X%dim, face_loc(U, face)) :: normal_shear_at_loc
     real, dimension(face_loc(X, face), face_loc(U, face)) :: mass

     ele    = face_ele(X, face) ! ele number for volume mesh
     dim    = mesh_dim(bed_shear_stress) ! field dimension 

     ! get shape functions
     f_shape => face_shape(U, face)     
     shape   => ele_shape(U, ele)     
     
     call transform_facet_to_physical(X, face, shape, ele_dshape_at_face_quad, &
                                      detwei_f = detwei, normal = normal)
    
     ! Calculate grad U at the surface element quadrature points 
     do i=1, dim
        do j=1, dim
           grad_U_at_quad(i, j, :) = &
                & matmul(ele_val(U, j, ele), ele_dshape_at_face_quad(:,:,i))
        end do
     end do

     visc_at_quad = face_val_at_quad(visc, face)
     X_ele = face_val_at_quad(X, face)
     do i_gi = 1, face_ngi(X, face)
        ! determine shear ( nu*(grad_u + grad_u.T) )   
        shear_at_quad(:,:,i_gi) = matmul(grad_U_at_quad(:,:,i_gi) + transpose(grad_U_at_quad(:,:,i_gi)), visc_at_quad(:,:,i_gi))

        ! Get absolute of normal vector
        do i = 1,dim
           abs_normal(i) = abs(normal(i,i_gi))
        end do

        ! Multiply by surface normal (dim,sgi) to obtain shear in direction normal
        ! to surface (not sure why it is transpose(shear) but this gives the
        ! correct answer?? sp911)
        normal_shear_at_quad(:,i_gi) = matmul(transpose(shear_at_quad(:,:,i_gi)), abs_normal) 
     end do  

     normal_shear_at_loc = shape_vector_rhs(f_shape, normal_shear_at_quad, density *&
          & detwei)

     ! for CG we need to calculate a global lumped mass
     mass = shape_shape(f_shape, f_shape, detwei)
     call addto(masslump, face_global_nodes(bed_shear_stress,face), sum(mass,1))

     ! add to bed_shear_stress field
     call addto(bed_shear_stress, face_global_nodes(bed_shear_stress,face), normal_shear_at_loc)

   end subroutine calculate_bed_shear_stress_ele_cg

   subroutine calculate_bed_shear_stress_ele_dg(bss, ele, X, grad_U, visc, density)

     type(vector_field), intent(inout) :: bss
     type(vector_field), intent(in), pointer :: X
     type(tensor_field), intent(in) :: visc
     type(tensor_field), intent(in) :: grad_U
     integer, intent(in) :: ele
     real, intent(in) :: density

     integer :: i, j, i_gi
     type(element_type), pointer :: shape
     real, dimension(ele_ngi(bss, ele)) :: detwei
     real, dimension(X%dim, ele_ngi(bss, ele)) :: normal, normal_shear_at_quad, X_at_quad
     real, dimension(X%dim) :: abs_normal
     real, dimension(X%dim, X%dim, ele_ngi(grad_U, ele)) :: grad_U_at_quad, visc_at_quad, shear_at_quad  
     real, dimension(X%dim, ele_loc(bss, ele)) :: rhs
     real, dimension(ele_loc(bss, ele), ele_loc(bss, ele)) :: inv_mass

     ! get shape functions
     shape => ele_shape(bss, ele)      
      
     call transform_facet_to_physical(X, ele, detwei_f = detwei, normal = normal)

     visc_at_quad = ele_val_at_quad(visc, ele)
     grad_U_at_quad = ele_val_at_quad(grad_U, ele)

     do i_gi = 1, ele_ngi(bss, ele)
        ! determine shear ( nu*(grad_u + grad_u.T ) )   
        shear_at_quad(:,:,i_gi) = density * matmul(grad_U_at_quad(:,:,i_gi) + transpose(grad_U_at_quad(:,:,i_gi)), visc_at_quad(:,:,i_gi))

        ! Get absolute of normal vector
        do i = 1, bss%dim
           abs_normal(i) = abs(normal(i,i_gi))
        end do

        ! Multiply by surface normal (dim,sgi) to obtain shear in direction normal
        ! to surface (not sure why it is transpose(shear) but this gives the
        ! correct answer?? sp911)
        normal_shear_at_quad(:,i_gi) = matmul(transpose(shear_at_quad(:,:,i_gi)), abs_normal) 
     end do  
     
     ! project on to basis functions to recover value at nodes
     rhs = shape_vector_rhs(shape, normal_shear_at_quad, detwei)
     inv_mass = inverse(shape_shape(shape, shape, detwei))
     do i = 1, X%dim
        rhs(i, :) = matmul(inv_mass, rhs(i, :))
     end do

     ! add to bss field
     call addto(bss, ele_nodes(bss,ele), rhs)
          
   end subroutine calculate_bed_shear_stress_ele_dg

   subroutine calculate_max_bed_shear_stress(state, max_bed_shear_stress)
!
      type(state_type), intent(in) :: state
      type(vector_field), intent(inout) :: max_bed_shear_stress

      type(vector_field), pointer :: bed_shear_stress
      type(scalar_field) :: magnitude_max_bss, magnitude_bss
      real :: current_time, spin_up_time
      integer stat, i

      call get_option(trim(max_bed_shear_stress%option_path)//"/spin_up_time", spin_up_time)
      call get_option("/timestepping/current_time", current_time)

      if(current_time>=spin_up_time) then
         
         ! Use the value already calculated previously
         bed_shear_stress => extract_vector_field(state, "BedShearStress", stat)  
         if(stat /= 0) then  
            ewrite(-1,*) "You need BedShearStress turned on to calculate MaxBedShearStress."
            FLExit("Turn on BedShearStress")
         end if

         ! We actually care about the vector that causes the maximum magnitude
         ! of bed shear stress, so check the magnitude and store if higher than
         ! what we already have.
         magnitude_max_bss = magnitude(max_bed_shear_stress)
         magnitude_bss = magnitude(bed_shear_stress)

         do i=1,node_count(magnitude_bss)
            if (node_val(magnitude_bss,i) .gt. node_val(magnitude_max_bss,i)) then
               call set(max_bed_shear_stress,i,node_val(bed_shear_stress,i))
            end if
         end do

         call deallocate(magnitude_max_bss)
         call deallocate(magnitude_bss)
      else
        call zero(max_bed_shear_stress)
      end if

   end subroutine calculate_max_bed_shear_stress

   subroutine calculate_galerkin_projection_scalar(state, field)
     type(state_type), intent(in) :: state
     type(scalar_field), intent(inout) :: field
       
     character(len=len_trim(field%option_path)) :: path
     character(len=FIELD_NAME_LEN) :: field_name
     type(scalar_field), pointer :: projected_field
     type(vector_field), pointer :: positions
     type(csr_sparsity) :: mass_sparsity
     type(csr_matrix) :: mass
     type(scalar_field) :: rhs, mass_lumped, inverse_mass_lumped
     logical :: dg
     logical :: check_integrals
     logical :: apply_bcs, lump_mass

     integer :: ele

     dg = (continuity(field) < 0)

     check_integrals = .false.
     apply_bcs = .true.

     path = field%option_path
     call get_option(path // "/diagnostic/source_field_name", field_name)
     lump_mass=have_option(path // "/diagnostic/lump_mass")
     projected_field => extract_scalar_field(state, trim(field_name))
     positions => extract_vector_field(state, "Coordinate")

     ! Assuming they're on the same quadrature
     assert(ele_ngi(field, 1) == ele_ngi(projected_field, 1))

     if (.not. dg .and. .not. lump_mass) then
       mass_sparsity = make_sparsity(field%mesh, field%mesh, name="MassMatrixSparsity")
       call allocate(mass, mass_sparsity, name="MassMatrix")
       call zero(mass)
     else if (lump_mass) then
       call allocate(mass_lumped, field%mesh, name="GalerkinProjectionMassLumped")
       call zero(mass_lumped)
     end if
     
     if (lump_mass .or. .not. dg) then
       call allocate(rhs, field%mesh, name="GalerkinProjectionRHS")
       call zero(rhs)
     end if

     do ele=1,ele_count(field)
       call assemble_galerkin_projection(field, projected_field, positions, &
                                      &  mass, rhs, ele, dg)
     end do

     if (lump_mass) then
       call allocate(inverse_mass_lumped, field%mesh, &
          name="GalerkinProjectionInverseMassLumped")
       call invert(mass_lumped, inverse_mass_lumped)
       call set(field, rhs)
       call scale(field, inverse_mass_lumped)
       call deallocate(mass_lumped)
       call deallocate(inverse_mass_lumped)
       call deallocate(rhs)
     else if (.not. dg) then
       call petsc_solve(field, mass, rhs, option_path=path // "/diagnostic")
       call deallocate(mass)
       call deallocate(mass_sparsity)
       call deallocate(rhs)
     end if

     if (check_integrals) then
       assert(field_integral(field, positions) .feq. field_integral(projected_field, positions))
     end if

     contains
     
       subroutine assemble_galerkin_projection(field, projected_field, positions, mass, rhs, ele, dg)
         type(scalar_field), intent(inout) :: field
         type(scalar_field), intent(in) :: projected_field
         type(vector_field), intent(in) :: positions
         type(csr_matrix), intent(inout) :: mass
         type(scalar_field), intent(inout) :: rhs
         integer, intent(in) :: ele
         logical, intent(in) :: dg

         type(element_type), pointer :: field_shape, proj_field_shape

         real, dimension(ele_loc(field, ele)) :: little_rhs
         real, dimension(ele_loc(field, ele), ele_loc(field, ele)) :: little_mass
         real, dimension(ele_loc(field, ele), ele_loc(projected_field, ele)) :: little_mba
         real, dimension(ele_loc(field, ele), ele_loc(projected_field, ele)) :: little_mba_int
         real, dimension(ele_ngi(field, ele)) :: detwei
         real, dimension(ele_loc(projected_field, ele)) :: proj_field_val

         integer :: i, j, k

         field_shape => ele_shape(field, ele)
         proj_field_shape => ele_shape(projected_field, ele)

         call transform_to_physical(positions, ele, detwei=detwei)

         little_mass = shape_shape(field_shape, field_shape, detwei)

         ! And compute the product of the basis functions
         little_mba = 0
         do i=1,ele_ngi(field, ele)
           forall(j=1:ele_loc(field, ele), k=1:ele_loc(projected_field, ele))
             little_mba_int(j, k) = field_shape%n(j, i) * proj_field_shape%n(k, i)
           end forall
           little_mba = little_mba + little_mba_int * detwei(i)
         end do

         proj_field_val = ele_val(projected_field, ele)
         little_rhs = matmul(little_mba, proj_field_val)
         
         if (lump_mass) then
           call addto(mass_lumped, ele_nodes(field, ele), &
             sum(little_mass,2))
           call addto(rhs, ele_nodes(field, ele), little_rhs)
         else if (dg) then
           call solve(little_mass, little_rhs)
           call set(field, ele_nodes(field, ele), little_rhs)
         else
           call addto(mass, ele_nodes(field, ele), ele_nodes(field, ele), little_mass)
           call addto(rhs, ele_nodes(field, ele), little_rhs)
         end if
         
       end subroutine assemble_galerkin_projection
      
   end subroutine calculate_galerkin_projection_scalar

   subroutine calculate_galerkin_projection_vector(state, field)
     type(state_type), intent(in) :: state
     type(vector_field), intent(inout) :: field
     character(len=len_trim(field%option_path)) :: path
     character(len=FIELD_NAME_LEN) :: field_name
     type(vector_field), pointer :: projected_field
     type(vector_field), pointer :: positions
     type(csr_sparsity) :: mass_sparsity
     type(csr_matrix) :: mass
     type(vector_field) :: rhs
     type(scalar_field) :: mass_lumped, inverse_mass_lumped
     logical :: dg
     logical :: check_integrals
     logical :: apply_bcs, lump_mass

     integer :: ele

     dg = (continuity(field) < 0)
     check_integrals = .false.
     apply_bcs = .true.

     path = field%option_path
     call get_option(path // "/diagnostic/source_field_name", field_name)
     lump_mass=have_option(path // "/diagnostic/lump_mass")
     projected_field => extract_vector_field(state, trim(field_name))
     positions => extract_vector_field(state, "Coordinate")

     ! Assuming they're on the same quadrature
     assert(ele_ngi(field, 1) == ele_ngi(projected_field, 1))

     if (.not. dg .and. .not. lump_mass) then
       mass_sparsity = make_sparsity(field%mesh, field%mesh, name="MassMatrixSparsity")
       call allocate(mass, mass_sparsity, name="MassMatrix")
       call zero(mass)
     else if (lump_mass) then
       call allocate(mass_lumped, field%mesh, name="GalerkinProjectionMassLumped")
       call zero(mass_lumped)
     end if
     
     if (lump_mass .or. .not. dg) then
       call allocate(rhs, field%dim, field%mesh, name="GalerkinProjectionRHS")
       call zero(rhs)
     end if

     do ele=1,ele_count(field)
       call assemble_galerkin_projection(field, projected_field, positions, &
                                      &  mass, rhs, ele, dg)
     end do

     if (lump_mass) then
       call allocate(inverse_mass_lumped, field%mesh, &
          name="GalerkinProjectionInverseMassLumped")
       call invert(mass_lumped, inverse_mass_lumped)
       call set(field, rhs)
       call scale(field, inverse_mass_lumped)
       call deallocate(mass_lumped)
       call deallocate(inverse_mass_lumped)
       call deallocate(rhs)
     else if (.not. dg) then
       call petsc_solve(field, mass, rhs, option_path=path // "/diagnostic")
       call deallocate(mass)
       call deallocate(mass_sparsity)
       call deallocate(rhs)
     end if

     if (check_integrals) then
       assert(all(abs(field_integral(field, positions) - field_integral(projected_field, positions)) < epsilon(0.0_4)))
     end if

     contains
     
       subroutine assemble_galerkin_projection(field, projected_field, positions, mass, rhs, ele, dg)
         type(vector_field), intent(inout) :: field
         type(vector_field), intent(in) :: projected_field
         type(vector_field), intent(in) :: positions
         type(csr_matrix), intent(inout) :: mass
         type(vector_field), intent(inout) :: rhs
         integer, intent(in) :: ele
         logical, intent(in) :: dg

         type(element_type), pointer :: field_shape, proj_field_shape

         real, dimension(ele_loc(field, ele), field%dim) :: little_rhs
         real, dimension(ele_loc(field, ele), ele_loc(field, ele)) :: little_mass
         real, dimension(ele_loc(field, ele), ele_loc(projected_field, ele)) :: little_mba
         real, dimension(ele_loc(field, ele), ele_loc(projected_field, ele)) :: little_mba_int
         real, dimension(ele_ngi(field, ele)) :: detwei
         real, dimension(field%dim, ele_loc(projected_field, ele)) :: proj_field_val

         integer :: i, j, k

         field_shape => ele_shape(field, ele)
         proj_field_shape => ele_shape(projected_field, ele)

         call transform_to_physical(positions, ele, detwei=detwei)

         little_mass = shape_shape(field_shape, field_shape, detwei)

         ! And compute the product of the basis functions
         little_mba = 0
         do i=1,ele_ngi(field, ele)
           forall(j=1:ele_loc(field, ele), k=1:ele_loc(projected_field, ele))
             little_mba_int(j, k) = field_shape%n(j, i) * proj_field_shape%n(k, i)
           end forall
           little_mba = little_mba + little_mba_int * detwei(i)
         end do

         proj_field_val = ele_val(projected_field, ele)
         do i=1,field%dim
           little_rhs(:, i) = matmul(little_mba, proj_field_val(i, :))
         end do

         if (lump_mass) then
           call addto(mass_lumped, ele_nodes(field, ele), &
             sum(little_mass,2))
           call addto(rhs, ele_nodes(field, ele), transpose(little_rhs))
         else if (dg) then
           call solve(little_mass, little_rhs)
           call set(field, ele_nodes(field, ele), transpose(little_rhs))
         else
           call addto(mass, ele_nodes(field, ele), ele_nodes(field, ele), little_mass)
           call addto(rhs, ele_nodes(field, ele), transpose(little_rhs))
         end if
         
       end subroutine assemble_galerkin_projection
         
   end subroutine calculate_galerkin_projection_vector

   subroutine calculate_universal_number(field)
     !!< Output the universal numbering associated with field. Clearly this
     !!< is primarily of interest for debugging.
     type(scalar_field) :: field
     
     integer i
     type(halo_type) :: halo

     if(.not.isparallel()) then
        do i=1, node_count(field)
           call set(field, i, real(i))
        end do

     else
        halo=field%mesh%halos(2)


        do i=1, node_count(field)
           
           call set(field, i, real(halo_universal_number(halo, i)))
           
        end do
     
     end if

   end subroutine calculate_universal_number

   subroutine calculate_node_owner(field)
     !!< Output the process owning each node in field. Clearly this
     !!< is primarily of interest for debugging.
     type(scalar_field) :: field
     
     integer i
     type(halo_type) :: halo

     if(.not.isparallel()) then
        do i=1, node_count(field)
           call set(field, i, real(i))
        end do

     else
        halo=field%mesh%halos(2)


        do i=1, node_count(field)
           
           call set(field, i, real(halo_node_owner(halo, i)))
           
        end do
     
     end if

   end subroutine calculate_node_owner
   
   subroutine calculate_galerkin_projection_tensor(state, field, solver_path)
     type(state_type), intent(in) :: state
     type(tensor_field), intent(inout) :: field
     character(len=*), intent(in), optional :: solver_path
     character(len=len_trim(field%option_path)) :: path
     character(len=FIELD_NAME_LEN) :: field_name
     type(tensor_field), pointer :: projected_field
     type(vector_field), pointer :: positions
     type(csr_sparsity) :: mass_sparsity
     type(csr_matrix) :: mass
     type(tensor_field) :: rhs
     logical :: dg
     logical :: check_integrals
     logical :: apply_bcs

     integer :: ele

     dg = (continuity(field) < 0)
     check_integrals = .false.
     apply_bcs = .true.

     path = field%option_path
     call get_option(path // "/diagnostic/source_field_name", field_name)
     projected_field => extract_tensor_field(state, trim(field_name))
     positions => extract_vector_field(state, "Coordinate")

     ! Assuming they're on the same quadrature
     assert(ele_ngi(field, 1) == ele_ngi(projected_field, 1))

     if (.not. dg) then
       mass_sparsity = make_sparsity(field%mesh, field%mesh, name="MassMatrixSparsity")
       call allocate(mass, mass_sparsity, name="MassMatrix")
       call zero(mass)
       call allocate(rhs, field%mesh, name="GalerkinProjectionRHS")
       call zero(rhs)
     end if

     do ele=1,ele_count(field)
       call assemble_galerkin_projection(field, projected_field, positions, &
                                      &  mass, rhs, ele, dg)
     end do

     if (.not. dg) then
       if (present(solver_path)) then
         call petsc_solve(field, mass, rhs, option_path=trim(solver_path))
       else
         call petsc_solve(field, mass, rhs, option_path=path // "/diagnostic")
       endif
       call deallocate(mass)
       call deallocate(mass_sparsity)
       call deallocate(rhs)
     end if

!     if (check_integrals) then
!       assert(all(abs(field_integral(field, positions) - field_integral(projected_field, positions)) < epsilon(0.0_4)))
!     end if

     contains
       subroutine assemble_galerkin_projection(field, projected_field, positions, mass, rhs, ele, dg)
         type(tensor_field), intent(inout) :: field
         type(tensor_field), intent(in) :: projected_field
         type(vector_field), intent(in) :: positions
         type(csr_matrix), intent(inout) :: mass
         type(tensor_field), intent(inout) :: rhs
         integer, intent(in) :: ele
         logical, intent(in) :: dg

         type(element_type), pointer :: field_shape, proj_field_shape

         real, dimension(field%dim(1), field%dim(2), ele_loc(field, ele)) :: little_rhs
         real, dimension(ele_loc(field, ele), ele_loc(field, ele)) :: little_mass
         real, dimension(ele_loc(field, ele), ele_loc(projected_field, ele)) :: little_mba
         real, dimension(ele_loc(field, ele), ele_loc(projected_field, ele)) :: little_mba_int
         real, dimension(ele_ngi(field, ele)) :: detwei
         real, dimension(field%dim(1), field%dim(2), ele_loc(projected_field, ele)) :: proj_field_val

         integer :: i, j, k

         field_shape => ele_shape(field, ele)
         proj_field_shape => ele_shape(projected_field, ele)

         call transform_to_physical(positions, ele, detwei=detwei)

         little_mass = shape_shape(field_shape, field_shape, detwei)

         ! And compute the product of the basis functions
         little_mba = 0
         do i=1,ele_ngi(field, ele)
           forall(j=1:ele_loc(field, ele), k=1:ele_loc(projected_field, ele))
             little_mba_int(j, k) = field_shape%n(j, i) * proj_field_shape%n(k, i)
           end forall
           little_mba = little_mba + little_mba_int * detwei(i)
         end do

         proj_field_val = ele_val(projected_field, ele)
         do i=1,field%dim(1)
           do j=1,field%dim(2)
             little_rhs(i, j, :) = matmul(little_mba, proj_field_val(i, j, :))
           end do
         end do

         if (dg) then
           FLAbort("You just need to write the appropriate solve interface.")
           !call solve(little_mass, little_rhs)
           call set(field, ele_nodes(field, ele), little_rhs)
         else
           call addto(mass, ele_nodes(field, ele), ele_nodes(field, ele), little_mass)
           call addto(rhs, ele_nodes(field, ele), little_rhs)
         end if
       end subroutine assemble_galerkin_projection
   end subroutine calculate_galerkin_projection_tensor
   
   subroutine calculate_diagnostic_coordinate_field(state, field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: field
    
    type(vector_field) :: coordinate_field
      
    coordinate_field = get_nodal_coordinate_field(state, field%mesh)
    call set(field, coordinate_field)
    call deallocate(coordinate_field)
   
   end subroutine calculate_diagnostic_coordinate_field
   
end module diagnostic_fields
