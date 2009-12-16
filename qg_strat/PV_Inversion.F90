#include "fdebug.h"
module pv_inversion
  use state_module
  use boundary_conditions_from_options
  use fields
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
  implicit none
  
  private

  !! coriolis parameter f
  real :: f
  !! coriolis parameter beta
  real :: beta
  !! Brunt-Vaisaila frequency
  real :: N
  !! initialisation flag
  logical :: initialised = .false.
  !! sparsity for streamfunction elliptic operator matrix
  type(csr_sparsity), save :: streamfunction_sparsity
  !! streamfunction elliptic operator matrix
  type(csr_matrix), save :: streamfunction_mat
  !! streamfunction field
  type(scalar_field), pointer :: streamfunction
  !! right-hand side field for PV inversion equation
  type(scalar_field), save :: RHS

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
    type(vector_field), pointer :: velocity, X
    !! element index
    integer :: ele

    ! Get streamfunction, coordinate, velocity
    streamfunction => extract_scalar_field(state,'Streamfunction')
    X => extract_vector_field(state,'Coordinate')    
    velocity => extract_vector_field(state,'NonlinearVelocity')

    ! wipe out velocity field
    call zero(velocity)

    ! loop over elements
    element_loop: do ele = 1, element_count(streamfunction)
       ! construct velocity on this element
       call streamfunction2velocity_ele(state,streamfunction,X,velocity,ele)
    end do element_loop

  end subroutine streamfunction2velocity

  subroutine streamfunction2velocity_ele(state,streamfunction,X,velocity,ele)
    !!< Compute the velocity from the streamfunction for a single element
    !! state variable
    type(state_type), intent(in) :: state
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
    real, dimension(mesh_dim(velocity), &
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
    V_loc(1,:) = -shape_rhs(v_shape,detwei*sfn_grad_quad(:,2))
    V_loc(2,:) =  shape_rhs(v_shape,detwei*sfn_grad_quad(:,1))
    V_loc(3,:) = 0.
    V_loc(1,:) = matmul(vmass_mat_loc,V_loc(1,:))
    V_loc(2,:) = matmul(vmass_mat_loc,V_loc(2,:))
    V_loc(3,:) = matmul(vmass_mat_loc,V_loc(3,:))
    !adding velocity values into velocity field
    !we can do this as velocity is discontinuous so nodes are not
    !shared between elements
    call set(velocity,v_ele,V_loc)

  end subroutine streamfunction2velocity_ele

  subroutine solve_streamfunction_qg(state, PV)
    !!< Compute the streamfunction from the PV
    !! state variable
    type(state_type), intent(inout) :: state
    !! potential vorticity
    type(scalar_field), intent(in) :: PV

    type(vector_field), pointer :: positions
    character(len=FIELD_NAME_LEN) :: bc_path, bc_name
    integer :: i, n_bcs

    if(.not.initialised) then
       call initialise_pv_inversion(state)       
    end if

    !allocate rhs and matrix
    streamfunction => extract_scalar_field(state, "Streamfunction")
    call allocate(RHS,streamfunction%mesh, "StreamfunctionRHS")
    RHS%option_path = streamfunction%option_path
    streamfunction_sparsity = make_sparsity(streamfunction%mesh, &
         streamfunction%mesh, "StreamfunctionSparisity")
    call allocate(streamfunction_mat,streamfunction_sparsity)

    call zero(RHS)

    call construct_streamfunction_qg(streamfunction_mat, RHS, PV, state)

    call apply_dirichlet_conditions(streamfunction_mat,RHS, &
         streamfunction, -1.0)

    n_bcs=option_count("/material_phase["//int2str(0)//&
         "]/scalar_field::Streamfunction/prognostic/boundary_conditions/&
         type::buoyancy")

    ! Find out whether we need to apply buoyancy boundary conditions
    do i=0, n_bcs-1
       ewrite(1,*) "there are buoyancy boundary conditions - good luck!"
       bc_path="/material_phase["//int2str(0)//&
            "]/scalar_field::Streamfunction/prognostic/boundary_conditions["&
            //int2str(i)//"]"
       if(have_option(trim(bc_path)//"/type::buoyancy")) then
          call get_option(trim(bc_path)//"/name", bc_name)
          positions=>extract_vector_field(state, "Coordinate")
          call apply_buoyancy_boundary_condition(RHS, streamfunction, &
               positions, bc_name)
       end if
    end do
    ewrite(1,*) "Out of buoyancy boundary condition loop"

    !zero streamfunction and solve for it
    call zero(streamfunction)

    ewrite(1,*) "Entering streamfunction solve"
    call petsc_solve(streamfunction,streamfunction_mat,RHS)

    call deallocate(streamfunction_mat)
    call deallocate(streamfunction_sparsity)
    call deallocate(RHS)

  end subroutine solve_streamfunction_qg

  subroutine construct_streamfunction_qg(streamfunction_mat, RHS, PV, state)
    !!< construct the matrix-vector equation for PV inversion
    !! matrix for PV inversion operator
    type(csr_matrix), intent(inout) :: streamfunction_mat
    !! right-hand side for pv inversion equation
    type(scalar_field), intent(inout) :: RHS
    !! potential vorticity
    type(scalar_field), intent(in) :: PV
    !! state variable
    type(state_type), intent(inout) :: state
    !
    integer :: ele
    type(vector_field), pointer :: X

    X => extract_vector_field(state,'Coordinate')

    element_loop: do ele = 1, element_count(PV)
       
       call construct_streamfunction_qg_element(ele,streamfunction_mat, &
            RHS, PV, state, X)

    end do element_loop

  end subroutine construct_streamfunction_qg

  subroutine construct_streamfunction_qg_element(ele,streamfunction_mat, &
       RHS, PV, state, X)
    !!< construct the matrix-vector equation for PV inversion for one element
    !! matrix for PV inversion equation
    type(csr_matrix), intent(inout) :: streamfunction_mat
    !! field for right-hand side of PV inversion equation
    type(scalar_field), intent(inout) :: RHS
    !! potential vorticity field
    type(scalar_field), intent(in) :: PV
    !! state variable
    type(state_type), intent(inout) :: state
    !! coordinates field
    type(vector_field), intent(in) :: X
    !! index of element
    integer, intent(in) :: ele
    !
    ! Variable transform times quadrature weights.
    real, dimension(ele_ngi(streamfunction,ele)) :: detwei, pv_quad    
    real, dimension(X%dim,ele_ngi(streamfunction,ele)) :: X_quad
    ! Transformed gradient function for streamfunction
    real, dimension(ele_loc(streamfunction, ele), & 
         ele_ngi(streamfunction,ele), &
         mesh_dim(streamfunction)) :: dstreamfunction
    ! Bilinear forms.
    real, dimension(ele_loc(streamfunction,ele), &
         ele_loc(streamfunction,ele)) :: &
         streamfunction_mat_loc, mass_mat_loc
    real, dimension(mesh_dim(streamfunction),mesh_dim(streamfunction), &
         ele_ngi(streamfunction,ele)) :: lengthscale_mat_quad
    ! node and shape pointers
    integer, dimension(:), pointer :: streamfunction_ele
    type(element_type), pointer :: streamfunction_shape, pv_shape

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

    !lengthscales in operator
    lengthscale_mat_quad = 0.
    lengthscale_mat_quad(1,1,:) = 1.
    lengthscale_mat_quad(2,2,:) = 1.
    lengthscale_mat_quad(3,3,:) = f*f/N/N

    streamfunction_mat_loc = -dshape_tensor_dshape(dstreamfunction, &
         lengthscale_mat_quad, &
         dstreamfunction, detwei)
    mass_mat_loc = shape_shape(streamfunction_shape,streamfunction_shape, &
         detwei)

    pv_quad=ele_val_at_quad(pv,ele)
    X_quad = ele_val_at_quad(X,ele)

    call addto(streamfunction_mat, &
         streamfunction_ele, &
         streamfunction_ele, streamfunction_mat_loc)

    call addto(RHS, streamfunction_ele, &
         shape_rhs(streamfunction_shape,(pv_quad-beta*X_quad(2,:))*detwei))

  end subroutine construct_streamfunction_qg_element

  subroutine initialise_pv_inversion(state)
    !!< set up PV inversion module
    !!state variable
    type(state_type), intent(inout) :: state
    real, dimension(3) :: beta_vector
    integer :: stat
    !these should be scalar fields 
    if(initialised) then
       !we will want to reinitialise after an adapt
       FLAbort('tried to initialise pv solve but already initialised')
    end if
    beta = 0.0
    call get_option("/physical_parameters/coriolis/f_plane/f_0",f, stat)
    if(stat/=0) then
       call get_option("/physical_parameters/coriolis/beta_plane/f_0",f)
       call get_option("/physical_parameters/coriolis/beta_plane/beta_vector",&
            & beta_vector)
       beta = beta_vector(2)
    end if
    call get_option("/physical_parameters/buoyancy_frequency",N)
    initialised = .true.
  end subroutine initialise_pv_inversion

  subroutine apply_buoyancy_boundary_condition(RHS, streamfunction, positions, bc_name)

    type(scalar_field), intent(inout) :: RHS
    type(scalar_field), intent(in) :: streamfunction
    type(vector_field), intent(in) :: positions
    character(len=*), intent(in) :: bc_name

    type(scalar_field), pointer :: buoyancy
    integer, dimension(:), pointer :: surface_element_list
    type(element_type), pointer:: streamfn_face_shape
    ! note that we assume all shapes to be the same in each element
    real, dimension(face_ngi(positions,1)) :: detwei_face
    real, dimension(face_loc(streamfunction,1)) :: buoyancy_val
    real, dimension(face_loc(streamfunction,1), face_loc(streamfunction,1)) :: face_mat
    integer, dimension(face_loc(streamfunction,1)) :: streamfn_face_nodes
    integer :: ele, face

    ewrite(1,*) "applying buoyancy boundary condition: ", trim(bc_name)

    buoyancy=>extract_surface_field(streamfunction, bc_name, &
         "value")

    ! Get list of surface elements on which to apply bc
    call get_boundary_condition(streamfunction, bc_name,&
         surface_element_list=surface_element_list)

    ! Loop over elements in surface mesh
    do ele=1, ele_count(buoyancy)

       buoyancy_val=ele_val(buoyancy, ele)

       face=surface_element_list(ele)

       call transform_facet_to_physical(positions, face, detwei_face)

       streamfn_face_shape=>face_shape(streamfunction, face)
       face_mat=shape_shape(streamfn_face_shape, streamfn_face_shape,&
            detwei_face)

       ! global node numbers of face
       streamfn_face_nodes=face_global_nodes(streamfunction, face)

       ! insert in RHS
       call addto(RHS, streamfn_face_nodes, matmul(face_mat, &
            buoyancy_val))

    end do

  end subroutine apply_buoyancy_boundary_condition

end module pv_inversion
