#include "fdebug.h"
module pv_inversion
  use state_module
  use boundary_conditions_from_options
  use fields
  use fields_base
  use fields_manipulation
  use SPUD
  use FLdebug
  use fetools
  use sparsity_patterns
  use state_module
  use solvers
  use vector_tools
  use boundary_conditions
  use transform_elements
  use global_parameters, only: OPTION_PATH_LEN
  implicit none
  
  private

  public :: solve_streamfunction_qg, streamfunction2velocity

contains

  subroutine streamfunction2velocity(state)
    !!< Compute the velocity from the streamfunction
    !!< It is done pointwise so velocity is DG
    !! State variable
    type(state_type), intent(in) :: state
    ! local variables
    !! streamfunction
    type(scalar_field), pointer :: streamfunction
    !! velocity and coordinates
    type(vector_field), pointer :: velocity, bg_velocity, X
    !! component of background velocity
    type(scalar_field) :: bg_comp
    !! element index
    integer :: ele, stat, i

    ! Get streamfunction, coordinate, velocity
    streamfunction => extract_scalar_field(state,'Streamfunction')
    X => extract_vector_field(state,'Coordinate')    
    velocity => extract_vector_field(state,'GeostrophicVelocity')

    ! wipe out velocity field
    call zero(velocity)

    ! loop over elements
    element_loop: do ele = 1, element_count(streamfunction)
       ! construct velocity on this element
       call streamfunction2velocity_ele(streamfunction, X, velocity, ele)
    end do element_loop

    bg_velocity => extract_vector_field(state,'BackgroundVelocity', stat)
    if(stat==0) then
       do i=1,mesh_dim(X)-1
          bg_comp = extract_scalar_field(bg_velocity, i, stat)
          call addto(velocity, i, bg_comp)
       end do
    end if

  end subroutine streamfunction2velocity

  subroutine streamfunction2velocity_ele(streamfunction, X, &
    velocity, ele)
    !!< Compute the velocity from the streamfunction for a single element
    !! streamfunction to be converted to velocity
    type(scalar_field), intent(in) :: streamfunction
    !! velocity to be converted from streamfunction
    type(vector_field), intent(inout) :: velocity
    !! coordinates field
    type(vector_field), intent(in) :: X
    !! index of element
    integer, intent(in) :: ele
    !
    ! Variable transform times quadrature weights.
    real, dimension(ele_ngi(streamfunction,ele)) :: detwei
    ! Transformed gradient function for streamfunction
    real, dimension(ele_loc(streamfunction, ele), & 
         ele_ngi(streamfunction,ele), &
         mesh_dim(streamfunction)) :: dstreamfunction
    ! Bilinear forms.
    real, dimension(ele_loc(velocity,ele), &
         ele_loc(velocity,ele)) :: vmass_mat_loc
    ! node and shape pointers
    integer, dimension(:), pointer :: streamfunction_ele, v_ele
    type(element_type), pointer :: streamfunction_shape, v_shape
    real, dimension(velocity%dim, &
         ele_loc(velocity,ele)) :: V_loc
    real, dimension(ele_ngi(streamfunction,ele),mesh_dim(velocity)) :: &
         sfn_grad_quad
 
    !----------------------------------------------------------------------
    ! Establish local node lists
    !----------------------------------------------------------------------

    streamfunction_ele=>ele_nodes(streamfunction,ele)
    v_ele=>ele_nodes(velocity,ele)

    !----------------------------------------------------------------------
    ! Establish local shape functions
    !----------------------------------------------------------------------

    streamfunction_shape=>ele_shape(streamfunction, ele)
    v_shape=>ele_shape(velocity, ele)

    ! Transform streamfunction derivatives and weights into physical space.
    call transform_to_physical(X, ele,&
         & streamfunction_shape , dshape=dstreamfunction, detwei=detwei)

    ! mass matrix for DG velocity, inverted
    vmass_mat_loc = shape_shape(v_shape,v_shape,detwei)
    call invert(vmass_mat_loc)
    ! gradients of stream function
    sfn_grad_quad=ele_grad_at_quad(streamfunction, ele,dstreamfunction)
    ! z-component of velocity is zero. If we're 3d then this is component 3.
    ! If we're 2d, check for presence of N to see if we're (x,y) or (x,z)
    if(velocity%dim==2) then
       V_loc(1,:) = -shape_rhs(v_shape,detwei*sfn_grad_quad(:,2))
       V_loc(1,:) = matmul(vmass_mat_loc,V_loc(1,:))
       V_loc(2,:) =  shape_rhs(v_shape,detwei*sfn_grad_quad(:,1))
       V_loc(2,:) = matmul(vmass_mat_loc,V_loc(2,:))
    else if(have_option("/physical_parameters/buoyancy_frequency")) then
       V_loc(1,:) = 0.
       V_loc(2,:) =  shape_rhs(v_shape,detwei*sfn_grad_quad(:,1))
       V_loc(2,:) = matmul(vmass_mat_loc,V_loc(2,:))
       V_loc(3,:) = 0.
    else
       V_loc(1,:) = -shape_rhs(v_shape,detwei*sfn_grad_quad(:,2))
       V_loc(1,:) = matmul(vmass_mat_loc,V_loc(1,:))
       V_loc(2,:) =  shape_rhs(v_shape,detwei*sfn_grad_quad(:,1))
       V_loc(2,:) = matmul(vmass_mat_loc,V_loc(2,:))
       V_loc(3,:) = 0.
    end if
    !adding velocity values into velocity field
    !we can do this as velocity is discontinuous so nodes are not
    !shared between elements
    call set(velocity,v_ele,V_loc)

  end subroutine streamfunction2velocity_ele

  subroutine solve_streamfunction_qg(state)
    !!< Compute the streamfunction from the PV
    !! state variable
    type(state_type), intent(inout) :: state

    type(scalar_field), pointer :: streamfunction
    type(scalar_field) :: RHS
    type(vector_field), pointer :: positions
    type(csr_sparsity) :: streamfunction_sparsity
    type(csr_matrix) :: streamfunction_mat

    !allocate rhs and matrix
    streamfunction => extract_scalar_field(state, "Streamfunction")
    call allocate(RHS,streamfunction%mesh, "StreamfunctionRHS")
    RHS%option_path = streamfunction%option_path
    streamfunction_sparsity = make_sparsity(streamfunction%mesh, &
         streamfunction%mesh, "StreamfunctionSparisity")
    call allocate(streamfunction_mat,streamfunction_sparsity)

    call zero(RHS)

    call construct_streamfunction_qg(streamfunction_mat, RHS, state)

    call apply_dirichlet_conditions(streamfunction_mat, RHS, &
         streamfunction)

    positions=>extract_vector_field(state, "Coordinate")
    call apply_buoyancy_boundary_conditions(RHS, streamfunction, positions)

    !zero streamfunction and solve for it
    call zero(streamfunction)

    ewrite(1,*) "Entering streamfunction solve"
    call petsc_solve(streamfunction, streamfunction_mat, RHS)

    call deallocate(streamfunction_mat)
    call deallocate(streamfunction_sparsity)
    call deallocate(RHS)

  end subroutine solve_streamfunction_qg

  subroutine construct_streamfunction_qg(streamfunction_mat, RHS, state)
    !!< construct the matrix-vector equation for PV inversion
    !! matrix for PV inversion operator
    type(csr_matrix), intent(inout) :: streamfunction_mat
    !! right-hand side for pv inversion equation
    type(scalar_field), intent(inout) :: RHS
    !! state variable
    type(state_type), intent(inout) :: state
    !
    integer :: ele
    type(scalar_field), pointer :: streamfunction, PV
    type(vector_field), pointer :: X, beta

    PV => extract_scalar_field(state, 'PotentialVorticity')
    streamfunction => extract_scalar_field(state, 'Streamfunction')
    X => extract_vector_field(state,'Coordinate')
    beta => extract_vector_field(state, "Beta")

    element_loop: do ele = 1, element_count(PV)
       
       call construct_streamfunction_qg_element(ele,streamfunction_mat, &
            RHS, PV, streamfunction, X, beta)

    end do element_loop

  end subroutine construct_streamfunction_qg

  subroutine construct_streamfunction_qg_element(ele, streamfunction_mat, &
       RHS, PV, streamfunction, X, beta)
    !!< construct the matrix-vector equation for PV inversion for one element
    !! matrix for PV inversion equation
    type(csr_matrix), intent(inout) :: streamfunction_mat
    !! field for right-hand side of PV inversion equation
    type(scalar_field), intent(inout) :: RHS
    !! potential vorticity field
    type(scalar_field), intent(in) :: PV, streamfunction
    !! coordinates field
    type(vector_field), intent(in) :: X, beta
    !! index of element
    integer, intent(in) :: ele
    !
    ! Variable transform times quadrature weights.
    real, dimension(ele_ngi(streamfunction,ele)) :: detwei, pv_quad    
    real, dimension(ele_ngi(streamfunction,ele)) :: coriolis_quad
    real, dimension(X%dim,ele_ngi(streamfunction,ele)) :: X_quad, beta_quad
    ! Transformed gradient function for streamfunction
    real, dimension(ele_loc(streamfunction, ele), & 
         ele_ngi(streamfunction,ele), &
         mesh_dim(streamfunction)) :: dstreamfunction
    ! Bilinear forms.
    real, dimension(ele_loc(streamfunction,ele), &
         ele_loc(streamfunction,ele)) :: &
         streamfunction_mat_loc
    real, dimension(mesh_dim(streamfunction),mesh_dim(streamfunction), &
         ele_ngi(streamfunction,ele)) :: lengthscale_mat_quad
    ! node and shape pointers
    integer, dimension(:), pointer :: streamfunction_ele
    type(element_type), pointer :: streamfunction_shape, pv_shape
    real, dimension(ele_loc(streamfunction, ele)) :: lrhs
    !! coriolis parameter f
    real :: f
    !! Brunt-Vaisaila frequency
    real :: N
    integer :: i, stat

    !----------------------------------------------------------------------
    ! Establish local node lists
    !----------------------------------------------------------------------

    streamfunction_ele=>ele_nodes(streamfunction,ele)

    !----------------------------------------------------------------------
    ! Establish local shape functions
    !----------------------------------------------------------------------

    streamfunction_shape=>ele_shape(streamfunction, ele)
    pv_shape=>ele_shape(pv, ele)

    ! Transform streamfunction derivatives and weights into physical space.
    call transform_to_physical(X, ele,&
         & streamfunction_shape , dshape=dstreamfunction, detwei=detwei)

    call get_option("/physical_parameters/coriolis/f_plane/f",f, stat)
    if(stat/=0) then
       call get_option("/physical_parameters/coriolis/beta_plane/f_0",f)
    end if

    !lengthscales in operator
    lengthscale_mat_quad = 0.
    lengthscale_mat_quad(1,1,:) = 1.
    call get_option("/physical_parameters/buoyancy_frequency",N, stat)
    ! If N present, check whether we're 2d x-z or 3d. If N not present then
    ! we're 2d x-y
    if(stat==0 .and. mesh_dim(streamfunction)==2) then
       lengthscale_mat_quad(2,2,:) = f*f/N/N
    else if(stat==0 .and. mesh_dim(streamfunction)==3) then
       lengthscale_mat_quad(2,2,:) = 1.
       lengthscale_mat_quad(3,3,:) = f*f/N/N
    else
       lengthscale_mat_quad(2,2,:) = 1.
    end if

    streamfunction_mat_loc = dshape_tensor_dshape(dstreamfunction, &
         lengthscale_mat_quad, &
         dstreamfunction, detwei)

    pv_quad=ele_val_at_quad(pv,ele)
    beta_quad=ele_val_at_quad(beta,ele)
    X_quad=ele_val_at_quad(X,ele)
    do i=1,size(coriolis_quad)
       coriolis_quad(i) = dot_product(beta_quad(:,i), X_quad(:,i))
    end do

    call addto(streamfunction_mat, &
         streamfunction_ele, &
         streamfunction_ele, streamfunction_mat_loc)

    lrhs=-shape_rhs(streamfunction_shape,(pv_quad-coriolis_quad)*detwei)
    call addto(RHS, streamfunction_ele, lrhs)

  end subroutine construct_streamfunction_qg_element

  subroutine apply_buoyancy_boundary_conditions(RHS, streamfunction, positions)

    type(scalar_field), intent(inout) :: RHS
    type(scalar_field), intent(in) :: streamfunction
    type(vector_field), intent(in) :: positions

    type(scalar_field), pointer :: buoyancy
    integer, dimension(:), pointer :: surface_element_list
    type(element_type), pointer:: streamfn_face_shape
    ! note that we assume all shapes to be the same in each element
    real, dimension(face_ngi(positions,1)) :: detwei_face
    real, dimension(face_loc(streamfunction,1)) :: bc_val
    real, dimension(face_loc(streamfunction,1), face_loc(streamfunction,1)) :: face_mat
    real :: f
    integer, dimension(face_loc(streamfunction,1)) :: streamfn_face_nodes
    integer :: ele, face, i, stat
    character(len=FIELD_NAME_LEN) :: bc_name, bc_type

    ! If there are no buoyancy boundary conditions, return.
    if (.not.(has_boundary_condition(streamfunction, "buoyancy"))) return

    call get_option("/physical_parameters/coriolis/f_plane/f",f, stat)
    if(stat/=0) then
       call get_option("/physical_parameters/coriolis/beta_plane/f_0",f)
    end if

    bcloop: do i=1, get_boundary_condition_count(streamfunction) 

       ! Get bc info
       call get_boundary_condition(streamfunction, i, name=bc_name,&
            type=bc_type, surface_element_list=surface_element_list)

       if (bc_type/="buoyancy") cycle bcloop

       ewrite(1,*) "applying buoyancy boundary condition: ", trim(bc_name)

       buoyancy=>extract_surface_field(streamfunction, bc_name, &
            "value")

       ! Loop over elements in surface mesh
       do ele=1, ele_count(buoyancy)

          bc_val=ele_val(buoyancy, ele)

          face=surface_element_list(ele)

          call transform_facet_to_physical(positions, face, detwei_face)

          streamfn_face_shape=>face_shape(streamfunction, face)
          face_mat=shape_shape(streamfn_face_shape, streamfn_face_shape,&
               detwei_face)

          ! global node numbers of face
          streamfn_face_nodes=face_global_nodes(streamfunction, face)

          ! insert in RHS
          call addto(RHS, streamfn_face_nodes, matmul(face_mat, &
               bc_val))

       end do

    end do bcloop

  end subroutine apply_buoyancy_boundary_conditions

end module pv_inversion
