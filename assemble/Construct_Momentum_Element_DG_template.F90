! Basic template for DG assembly routines. This is used to generate optimised versions.
! For simplicity, removed following functionality:
!   -   Arbitrary upwind viscosity scheme

#define ELEMENT_CONFIG

#include "element_type_defs.h"

!------------------ Start of template code ------------------------------

subroutine construct_momentum_elements_dg_ELEMENT_CONFIG( ele, big_m, rhs, &
    &X, U, U_nl, U_mesh, X_old, X_new, &
    & u_shape, p_shape, q_shape, &
    & Source, Buoyancy, hb_density, hb_pressure, gravity, Abs, &
    &Viscosity, swe_bottom_drag, swe_u_nl, P, old_pressure, Rho, surfacetension, q_mesh, &
    &velocity_bc, velocity_bc_type, &
    &pressure_bc, pressure_bc_type, &
    &turbine_conn_mesh, on_sphere, depth, have_wd_abs, alpha_u_field, Abs_wd, &
    &vvr_sf, ib_min_grad, nvfrac, &
    &inverse_mass, inverse_masslump, mass, subcycle_m, partial_stress, &
    have_les, have_isotropic_les, &
    smagorinsky_coefficient, eddy_visc, tensor_eddy_visc, &
    prescribed_filter_width, distance_to_wall, &
    y_plus_debug, les_filter_width_debug )

    !!< Construct the momentum equation for discontinuous elements in
    !!< acceleration form.
    implicit none

    !    type(integer_set), dimension(:), intent(in), pointer :: colours

    !! Index of current element
    integer, intent(in) :: ele
    !! Main momentum matrix.
    type(petsc_csr_matrix), intent(inout) :: big_m
    !! Momentum right hand side vector for each point.
    type(vector_field), intent(inout) :: rhs
    !! Auxiliary variable mesh
    type(mesh_type), intent(in) :: q_mesh
    type(mesh_type), intent(in) :: turbine_conn_mesh
    !!
    type(block_csr_matrix), intent(inout), optional :: subcycle_m

    !! Position, velocity and source fields.
    type(scalar_field), intent(in) :: buoyancy
    type(vector_field), intent(in) :: X, U, U_nl, Source, gravity, Abs
    type(vector_field), pointer :: U_mesh, X_old, X_new
    !! Viscosity
    type(tensor_field) :: Viscosity
    type(scalar_field) :: P, Rho
    type(scalar_field), intent(in) :: hb_density, hb_pressure
    !! surfacetension
    type(tensor_field) :: surfacetension
    !! field containing the bc values of velocity
    type(vector_field), intent(in) :: velocity_bc
    !! array of the type of bc (see get_entire_boundary_condition call above)
    integer, dimension(:,:), intent(in) :: velocity_bc_type
    !! same for pressure
    type(scalar_field), intent(in) :: pressure_bc
    integer, dimension(:), intent(in) :: pressure_bc_type
    !! fields only used for swe bottom drag (otherwise unitialised)
    type(scalar_field), intent(in) :: swe_bottom_drag, old_pressure
    type(vector_field), intent(in) :: swe_u_nl

    !! Inverse mass matrix
    type(block_csr_matrix), intent(inout), optional :: inverse_mass
    !! Mass lumping for each point
    type(vector_field), intent(inout), optional :: inverse_masslump
    !! Optional separate mass matrix.
    type(csr_matrix), intent(inout), optional :: mass
    logical, intent(in) :: have_wd_abs !! Wetting and drying switch, if TRUE, alpha_u_field must be passed as well
    type(scalar_field), intent(in) :: alpha_u_field
    type(vector_field), intent(in) :: Abs_wd

    ! Bilinear forms.

    real, dimension(NLOC, NLOC) :: Coriolis_mat, rho_mat, rho_move_mat, mass_mat, inverse_mass_mat, Advection_mat, Source_mat

    real, dimension(NDIM, NLOC, NLOC) :: ele2grad_mat, Abs_mat
    real, dimension(NDIM, NDIM, NLOC, NLOC) :: Abs_mat_sphere

    real, dimension(NDIM, NLOC) :: Abs_lump
    real, dimension(NDIM, NDIM, NLOC) :: Abs_lump_sphere

    real, dimension(NLOC) :: source_lump

    real, dimension(NLOC, NLOC) :: Q_inv
    real, dimension(NDIM, NLOC, EFLOC) :: Grad_u_mat_q, Div_u_mat_q
    real, dimension(NDIM, NDIM, EFLOC, EFLOC) :: Viscosity_mat
    real, dimension(NDIM) :: node_stress_diag, resid_stress_term
    real :: visc_dot_prod
    real, dimension(NDIM, NDIM) :: visc_grad_dot_u
    real, dimension(NDIM, NDIM, NLOC) :: Viscosity_ele
    real, dimension(NDIM, NDIM, NGI) :: visc_ele_quad

    real, dimension(NDIM, NLOC) :: x_val, x_val_2, u_val

    ! \Int_{ele} N_i kappa N_j dV, used for CDG fluxes
    real, dimension(NDIM, NDIM, NLOC, NLOC) :: kappa_mat

    ! Local assembly matrices.
    real, dimension(NLOC) :: l_MassLump, l_move_masslump

    ! Local node number map for 2nd order element.
    integer, dimension(EFLOC) :: local_glno

    ! Local variables.

    ! Neighbour element, face, neighbour face, no. internal element nodes
    integer :: ele_2, ele_2_X, face, face_2, loc
    ! Count variable for loops over dimension.
    integer :: dim, dim1, dim2, dim3, dim4
    ! Loops over faces.
    integer :: ni
    ! Array bounds for faces of the 2nd order element.
    integer :: start, finish

    ! Variable transform times quadrature weights.
    real, dimension(NGI) :: detwei, detwei_old, detwei_new, coefficient_detwei
    ! Transformed gradient function for velocity.
    ! Transformed gradient function for grid velocity.
    ! Transformed gradient function for auxiliary variable.
    real, dimension(NLOC, NGI, NDIM) :: du_t, dug_t, dq_t

    ! Density at quadrature points.
    ! Coriolis magnitude and sign at quadrature points.

    real, dimension(NGI) :: Rho_q, Coriolis_q
    ! Different velocities at quad points.
    real, dimension(NDIM, NGI) :: u_nl_q
    real, dimension(NGI) :: u_nl_div_q

    ! surface tension terms
    real, dimension(NDIM, NDIM, NGI) :: tension
    real, dimension(NDIM, NGI) :: dtensiondj

    ! Node and shape pointers.

    ! JRM: Can you avoid using pointers?
    ! integer, dimension(:), pointer :: u_ele, p_ele
    integer, dimension(NLOC ):: u_ele, p_ele, x_ele
    type(element_type), intent(in), pointer :: u_shape, p_shape, q_shape
    !    type(element_type), pointer :: u_face_shape, p_face_shape, q_face_shape

    ! Neighbours of this element.
    integer, dimension(:), pointer :: neigh, X_neigh
    ! Whether the velocity field is continuous and if it is piecewise constant.
    logical :: dg, p0
    integer :: i, gi
    logical :: boundary_element, turbine_face

    ! What we will be adding to the matrix and RHS - assemble these as we
    ! go, so that we only do the calculations we really need
    real, dimension(NDIM, EFLOC) :: big_m_diag_addto, rhs_addto

    real, dimension(NDIM, NDIM, EFLOC, EFLOC) :: big_m_tensor_addto
    logical, dimension(NDIM, NDIM) :: diagonal_block_mask, off_diagonal_block_mask
    ! Addto matrices for when subcycling is performed
    real, dimension(NDIM, NDIM, EFLOC, EFLOC) :: subcycle_m_tensor_addto

    !Switch to select if we are assembling the primal or dual form
    logical :: primal

    ! In parallel, we assemble terms on elements we own, and those in
    ! the L1 element halo
    logical :: assemble_element

    ! If on the sphere evaluate gravity direction at the gauss points
    logical :: on_sphere

    ! Absorption matrices
    real, dimension(NDIM, NGI) :: absorption_gi
    real, dimension(NDIM, NDIM, NGI) :: tensor_absorption_gi

    ! Add vertical velocity relaxation to the absorption if present
    real, intent(in) :: vvr_sf
    real, dimension(NDIM, NDIM, NGI) :: vvr_abs
    real, dimension(NDIM, NGI) :: vvr_abs_diag
    real, dimension(NGI) :: depth_at_quads
    type(scalar_field), intent(in) :: depth

    ! Add implicit buoyancy to the absorption if present
    real, intent(in) :: ib_min_grad
    real, dimension(NDIM, NDIM, NGI) :: ib_abs
    real, dimension(NDIM, NGI) :: ib_abs_diag
    real, dimension(NLOC, NGI, NDIM) :: dt_rho
    real, dimension(NDIM, NGI) :: grav_at_quads, grad_rho
    real, dimension(NDIM, NLOC) :: ele_grav_val
    real, dimension(NGI) :: drho_dz

!    ! Non-linear approximation to the PhaseVolumeFraction field
    type(scalar_field), intent(in) :: nvfrac
!    type(element_type), pointer :: nvfrac_shape
!    ! Transformed gradient function for the non-linear PhaseVolumeFraction.
!    real, dimension(NLOC, NGI, NDIM) :: dnvfrac_t
!    ! real, dimension(,:,:), allocatable :: dnvfrac_t
!    ! nvfrac at quadrature points.
!    real, dimension(NGI) :: nvfrac_gi, u_nl_dot_grad_nvfrac_gi
!    real, dimension(NDIM, NGI) :: grad_nvfrac_gi

    ! Moving mesh
    real, dimension(NDIM, NGI) :: ele_u_mesh_quad

    ! element centre and neighbour centre
    ! for IP parameters

    real, dimension(NDIM) :: ele_centre, neigh_centre, face_centre, face_centre_2
    real :: turbine_fluxfac

    real, dimension(NGI) :: alpha_u_quad

    ! added for partial stress form (sp911)
    logical, intent(in) :: partial_stress

    ! LES - sp911
    logical, intent(inout) :: have_les, have_isotropic_les
    real, intent(in) :: smagorinsky_coefficient
    type(scalar_field), pointer, intent(inout) :: eddy_visc, y_plus_debug, &
        & les_filter_width_debug
    type(tensor_field), pointer, intent(inout) :: tensor_eddy_visc
    type(scalar_field), pointer, intent(in) :: prescribed_filter_width, distance_to_wall
    real, dimension(NDIM, NDIM, FLOC) :: tmp_face_tensor

    integer :: iloc, jloc, idim, jdim

    integer :: nmat_d1,nmat_d2



     ! ========== INTERNAL INTERFACE VARIABLES ==========

    logical :: CDG_switch_in

    ! Matrix for assembling primal fluxes
    ! Note that this assumes same order polys in each element
    ! Code will need reorganising for p-refinement
    real, dimension(2, FLOC, NLOC) :: face_primal_fluxes_mat
    real, dimension(FLOC, FLOC) :: face_shape_shape_work

    ! Matrix for assembling penalty fluxes
    ! Note that this assumes same order polys in each element
    ! Code will need reorganising for p-refinement
    real, dimension(2, FLOC, FLOC) :: face_penalty_fluxes_mat

    ! \Int_{s_ele} N_iN_j n ds, used for CDG fluxes
    !real, dimension(mesh_dim(U),ele_loc(U,ele),ele_loc(U,ele)) :: &
    !    & normal_mat
    ! I think the above is wrong, and the dimensions below are correct.
    real, dimension(NDIM, FLOC, FLOC) :: face_normal_mat

    ! \Int_{s_ele} N_iN_j kappa.n ds, used for CDG fluxes
    ! Note that this assumes same order polys in each element
    ! Code will need reorganising for p-refinement
    !    real, dimension(mesh_dim(U),face_loc(U,face),face_loc(U,face)) :: &
    !        & kappa_normal_mat
    real, dimension(NDIM, FLOC, FLOC) :: face_kappa_normal_mat

    ! Face objects and numberings.
    ! type(element_type), intent(in), pointer :: u_shape, u_shape_2, p_shape, q_shape
    type(element_type) :: face_u_shape, face_p_shape, face_q_shape
    integer, dimension(FLOC) :: u_face_l, u_mesh_glno, u_face_glno_1, u_face_glno_2, Rho_face_glno_1

    integer, dimension(FLOC) :: q_face_l

    ! Note that both sides of the face can be assumed to have the same
    ! number of quadrature points.
    real, dimension(FNGI) :: face_Rho_q !, face_nvfrac_gi
    real, dimension(FLOC) :: face_Rho_val
    real, dimension(NDIM, FNGI) :: face_normal, face_u_nl_q, face_u_f_q, face_u_f2_q, face_div_u_f_q
    real, dimension(NDIM, FNGI) :: face_u_mesh_quad

    logical, dimension(FNGI) :: face_inflow
    real, dimension(FNGI) :: face_u_nl_q_dotn, face_income
    ! Variable transform times quadrature weights.
    real, dimension(FNGI) :: face_detwei, face_detwei_work
    real, dimension(FNGI) :: face_inner_advection_integral, face_outer_advection_integral

    ! Bilinear forms
    !real, dimension(face_loc(U,face),face_loc(U,face)) :: nnAdvection_out
    !real, dimension(face_loc(U,face),face_loc(U,face_2)) :: nnAdvection_in
    ! real, dimension(1, mesh_dim(U), face_loc(P,face),face_loc(U,face)) :: mnCT
    real, dimension(FLOC, FLOC) :: face_nnAdvection_out
    real, dimension(FLOC, FLOC) :: face_nnAdvection_in

    ! This is almost good enough.
    real, dimension(1, NDIM, P_FLOC, FLOC) :: face_mnCT


    ! Viscosity values on face (used for CDG and IP fluxes)
    real, dimension(NDIM, NDIM, FNGI) :: face_kappa_gi
    real, dimension(NDIM, NDIM, FLOC) :: face_visc_val

    ! surfacetension stuff
    real, dimension(NDIM, NDIM, FNGI) :: face_tension_q

    integer :: face_dim, face_start, face_finish
    logical :: face_boundary, face_free_surface, face_no_normal_flow, face_l_have_pressure_bc
    logical, dimension(NDIM) :: face_dirichlet

    integer :: face_d1, face_d2

    ! ========== END OF INTERFACE VARIABLES ==========

    have_les=.false.

    dg=.true.
    p0=.false.

#if defined (SCHEME_CDG)
    !    select case (viscosity_scheme)
    !        case (CDG)
            primal = .true.

#elif defined (SCHEME_IP)
    !        case (IP)
            primal = .true.
#else
    !        case default
    primal = .false.
    !    end select
#endif


    !    u_face_shape => U%mesh%faces%shape
    !    p_face_shape => P%mesh%faces%shape
    !    q_face_shape => q_mesh%faces%shape

    ! In parallel, we construct terms on elements we own and those in
    ! the L1 element halo.
    ! Note that element_neighbour_owned(U, ele) may return .false. if
    ! ele is owned.  For example, if ele is the only owned element on
    ! this process.  Hence we have to check for element ownership
    ! directly as well.
    assemble_element = .not.dg.or.element_neighbour_owned(U, ele).or.element_owned(U, ele)

    big_m_diag_addto = 0.0
    big_m_tensor_addto = 0.0
    if(subcycle) then
        subcycle_m_tensor_addto = 0.0
    end if

    rhs_addto = 0.0

    !----------------------------------------------------------------------
    ! Establish local node lists
    !----------------------------------------------------------------------

    u_ele=ele_nodes(U,ele)  ! Velocity
    p_ele=ele_nodes(P,ele)  ! Pressure
    x_ele=ele_nodes(X,ele) ! Coords

    local_glno(:NLOC)=u_ele ! Viscosity node list


    !print*, "u_shape%n(:,:):", u_shape%n(:,:)

    ! x_val = ele_val(X,ele)
    do concurrent (i=1:NDIM)
        x_val(i, :) = X%val(i, x_ele)
    end do

    ! Transform U derivatives and weights into physical space.
    ! if(.not.p0) then
    call transform_to_physical(X, ele,&
        & u_shape , dshape=du_t, detwei=detwei)
    !    else
    !        call transform_to_physical(X, ele, &
    !            & detwei=detwei)
    !        du_t = 0.0
    !    end if

    !print*, "ele, detwei:", ele, detwei


    if(move_mesh) then
        ele_u_mesh_quad=matmul(U_mesh%val(:, x_ele), u_shape%n)

        call transform_to_physical(X_old, ele, &
            & detwei=detwei_old)
        call transform_to_physical(X_new, ele, &
            & detwei=detwei_new)
        if(have_advection.and..not.integrate_by_parts_once) then
            call transform_to_physical(X, ele, &
                & U_shape, dshape = dug_t)
        end if
    end if

    if(have_viscosity.and.(.not.(q_mesh==u%mesh))) then
        ! Transform q derivatives into physical space.
        call transform_to_physical(X, ele,&
            & q_shape , dshape=dq_t)
    else
        dq_t=du_t
    end if

    !----------------------------------------------------------------------
    ! Construct element-wise quantities.
    !----------------------------------------------------------------------

                ! JRM: Start inlining?
    Rho_q=ele_val_at_quad(Rho, ele)

!    if(multiphase) then
!
!        ! If the Velocity and PhaseVolumeFraction meshes are different, then we need to
!        ! compute the derivatives of the PhaseVolumeFraction shape functions.
!        if(.not.(nvfrac%mesh == u%mesh)) then
!            nvfrac_shape => ele_shape(nvfrac%mesh, ele)
!            call transform_to_physical(X, ele, nvfrac_shape, dshape=dnvfrac_t)
!        else
!            dnvfrac_t = du_t
!        end if
!
!        nvfrac_gi = ele_val_at_quad(nvfrac, ele)
!        grad_nvfrac_gi = ele_grad_at_quad(nvfrac, ele, dnvfrac_t)
!    end if

    ! Viscosity matrices
    if(Viscosity%field_type==FIELD_TYPE_CONSTANT) then
        do concurrent(iloc=1:NLOC)
            Viscosity_ele(:, :, iloc)=Viscosity%val(:, :, 1)
        end do
    else
        ! Field type is FIELD_TYPE_NORMAL
        Viscosity_ele=Viscosity%val(:,:,u_ele)
    end if

    if (assemble_element) then
        ! Isotropic / scalar LES
        if(have_les) then
            if(have_isotropic_les) then
                do concurrent(dim1=1:NDIM)
                    Viscosity_ele(dim1, dim1, :) = &
                        Viscosity_ele(dim1, dim1, :)+eddy_visc%val(u_ele)
                end do
            else
                ! Anisotropic / tensor LES
                Viscosity_ele = Viscosity_ele + tensor_eddy_visc%val(:,:, u_ele)
             end if
        end if
        u_val = u%val(:, u_ele)
    end if
    ! visc_ele_quad=ele_val_at_quad(Viscosity,ele)
    visc_ele_quad = tensormul(Viscosity_ele, u_shape%n)


    !----------------------------------------------------------------------
    ! Construct bilinear forms.
    !----------------------------------------------------------------------

    ! Element density matrix.
    ! (compute for first component only at first, others are copied
    !  when necessary)
    !
    ! POTENTIAL OPTIMISATION.
    ! This doesn't have to be done _at_all_ unless the mesh has changed.
    ! Cache values?

    if (move_mesh) then
        ! this rho_mat (and l_masslump) is only used in the actual mass term in big_m
        ! (and its derivative inverse_mass or inverse_mass_lump)
        ! so should be evaluated at t+dt
        ! rho_mat = shape_shape(u_shape, u_shape, detwei_new*Rho_q)
        do concurrent(iloc=1:NLOC, jloc=1:NLOC)
            ! Main mass matrix.
            rho_mat(iloc,jloc)=dot_product(u_shape%n(iloc,:)*u_shape%n(jloc,:), detwei_new*Rho_q)
        end do
    else

!        if(multiphase) then
!            rho_mat = shape_shape(u_shape, u_shape, detwei*Rho_q*nvfrac_gi)
!        else
            rho_mat = shape_shape(u_shape, u_shape, detwei*Rho_q)
!        end if

    end if
    l_masslump= sum(rho_mat,2)

    if(present(mass)) then
        ! Return mass separately.
        ! NOTE: this doesn't deal with mesh movement
        call addto(mass, u_ele, u_ele, Rho_mat)
    else
        if(have_mass.and.assemble_element) then
            if(lump_mass) then
                do dim = 1, NDIM
                    big_m_diag_addto(dim, :NLOC) = big_m_diag_addto(dim, :NLOC) + l_masslump
                end do
            else
                do dim = 1, NDIM
                    big_m_tensor_addto(dim, dim, :NLOC, :NLOC) = big_m_tensor_addto(dim, dim, :NLOC, :NLOC) + rho_mat
                end do
            end if
        end if
        if (move_mesh.and.assemble_element) then
            ! In the unaccelerated form we solve:
            !  /
            !  |  N^{n+1} u^{n+1}/dt - N^{n} u^n/dt + ... = f
            !  /
            ! so in accelerated form:
            !  /
            !  |  N^{n+1} du + (N^{n+1}- N^{n}) u^n/dt + ... = f
            !  /
            ! where du=(u^{n+1}-u^{n})/dt is the acceleration.
            ! Put the (N^{n+1}-N^{n}) u^n term on the rhs
            rho_move_mat = shape_shape(u_shape, u_shape, (detwei_new-detwei_old)*Rho_q)
            if(lump_mass) then
                l_move_masslump= sum(rho_move_mat,2)
                !                do dim = 1, NDIM
                !                    rhs_addto(dim,:NLOC) = rhs_addto(dim,:NLOC) - l_move_masslump*u_val(dim,:)/dt
                !                end do
                rhs_addto(1,:NLOC) = rhs_addto(1,:NLOC) - l_move_masslump*u_val(1,:)/dt
                rhs_addto(2,:NLOC) = rhs_addto(2,:NLOC) - l_move_masslump*u_val(2,:)/dt
                rhs_addto(3,:NLOC) = rhs_addto(3,:NLOC) - l_move_masslump*u_val(3,:)/dt

            else
                !                do dim = 1, NDIM
                !                    rhs_addto(dim,:NLOC) = rhs_addto(dim,:NLOC) - matmul(rho_move_mat, u_val(dim,:))/dt
                !                end do
                rhs_addto(1,:NLOC) = rhs_addto(1,:NLOC) - matmul(rho_move_mat, u_val(1,:))/dt
                rhs_addto(2,:NLOC) = rhs_addto(2,:NLOC) - matmul(rho_move_mat, u_val(2,:))/dt
                rhs_addto(3,:NLOC) = rhs_addto(3,:NLOC) - matmul(rho_move_mat, u_val(3,:))/dt
            end if
        end if
    end if

    if(have_coriolis.and.(rhs%dim>1).and.assemble_element) then
        Coriolis_q=coriolis(ele_val_at_quad(X,ele))

        ! Element Coriolis parameter matrix.
        Coriolis_mat = shape_shape(u_shape, u_shape, Rho_q*Coriolis_q*detwei)

        ! cross terms in U_ and V_ for coriolis
        big_m_tensor_addto(U_, V_, :NLOC, :NLOC) = big_m_tensor_addto(U_, V_, :NLOC, :NLOC) - dt*theta*coriolis_mat
        big_m_tensor_addto(V_, U_, :NLOC, :NLOC) = big_m_tensor_addto(V_, U_, :NLOC, :NLOC) + dt*theta*coriolis_mat

        if(acceleration)then
            rhs_addto(U_, :NLOC) = rhs_addto(U_, :NLOC) + matmul(coriolis_mat, u_val(V_,:))
            rhs_addto(V_, :NLOC) = rhs_addto(V_, :NLOC) - matmul(coriolis_mat, u_val(U_,:))
        end if
    end if

    if(have_advection.and.assemble_element) then
        ! Advecting velocity at quadrature points.

        ! U_nl_q=ele_val_at_quad(U_nl,ele)
        U_nl_q = matmul(U_nl%val(:, u_ele), u_shape%n)

        if(integrate_conservation_term_by_parts) then

!            if(multiphase) then
!                ! Element advection matrix
!                !         /                                                /
!                !  - beta | (grad T dot U_nl) T Rho vfrac dV + (1. - beta) | T (vfrac U_nl dot grad T) Rho dV
!                !         /                                                /
!                Advection_mat = -beta*dshape_dot_vector_shape(du_t, U_nl_q, u_shape, detwei*Rho_q*nvfrac_gi) &
!                    + (1.-beta)*shape_vector_dot_dshape(u_shape, U_nl_q, du_t, detwei*Rho_q*nvfrac_gi)
!            else
                ! Element advection matrix
                !         /                                          /
                !  - beta | (grad T dot U_nl) T Rho dV + (1. - beta) | T (U_nl dot grad T) Rho dV
                !         /                                          /
                Advection_mat = -beta*dshape_dot_vector_shape(du_t, U_nl_q, u_shape, detwei*Rho_q) &
                    + (1.-beta)*shape_vector_dot_dshape(u_shape, U_nl_q, du_t, detwei*Rho_q)
!            end if

            if(move_mesh) then
                if(integrate_by_parts_once) then
!                    Advection_mat = Advection_mat &
!                        + dshape_dot_vector_shape(du_t, ele_val_at_quad(U_mesh,ele), u_shape, detwei * Rho_q)
                    Advection_mat = Advection_mat &
                        + dshape_dot_vector_shape(du_t, ele_u_mesh_quad, &
                        u_shape, detwei * Rho_q)
                else
!                    Advection_mat = Advection_mat &
!                        - shape_vector_dot_dshape(u_shape, ele_val_at_quad(U_mesh,ele), du_t, detwei * Rho_q) &
!                        - shape_shape(u_shape, u_shape, ele_div_at_quad(U_mesh, ele, dug_t) * detwei * Rho_q)
                    Advection_mat = Advection_mat &
                        - shape_vector_dot_dshape(u_shape, ele_u_mesh_quad, &
                            du_t, detwei * Rho_q) &
                        - shape_shape(u_shape, u_shape, ele_div_at_quad(U_mesh, ele, dug_t) * detwei * Rho_q)
                end if
            end if
        else
            ! Introduce grid velocities
            if (move_mesh) then
                ! NOTE: this modifies the velocities stored at the gauss pts.
                U_nl_q = U_nl_q - ele_u_mesh_quad
            end if
            U_nl_div_q=ele_div_at_quad(U_nl, ele, du_t)

            if(integrate_by_parts_once) then

!                if(multiphase) then
!                    ! Element advection matrix
!                    !    /                                                /
!                    !  - | (grad T dot U_nl vfrac) T Rho dV - (1. - beta) | T ( div(U_nl vfrac) ) T Rho dV
!                    !    /                                                /
!
!                    ! We need to compute \int{T div(u_nl vfrac) T},
!                    ! so split up the div using the product rule and compute
!                    ! \int{T vfrac div(u_nl) T} + \int{T u_nl grad(vfrac) T}
!                    do i = 1, NGI
!                        u_nl_dot_grad_nvfrac_gi(i) = dot_product(U_nl_q(:,i), grad_nvfrac_gi(:,i))
!                    end do
!                    Advection_mat = -dshape_dot_vector_shape(du_t, U_nl_q, u_shape, detwei*Rho_q*nvfrac_gi) &
!                        - (1.-beta) * (shape_shape(u_shape, u_shape, U_nl_div_q*detwei*Rho_q*nvfrac_gi) + &
!                        shape_shape(u_shape, u_shape, detwei*Rho_q*u_nl_dot_grad_nvfrac_gi))
!                else
                    ! Element advection matrix
                    !    /                                          /
                    !  - | (grad T dot U_nl) T Rho dV - (1. - beta) | T ( div U_nl ) T Rho dV
                    !    /                                          /
                    Advection_mat = - dshape_dot_vector_shape(du_t, U_nl_q, u_shape, detwei*Rho_q) &
                        - (1.-beta) * shape_shape(u_shape, u_shape, U_nl_div_q*detwei*Rho_q)
!                end if

            else

!                if(multiphase) then
!                    ! Element advection matrix
!                    !  /                                         /
!                    !  | T (vfrac U_nl dot grad T) Rho dV + beta | T ( div (vfrac U_nl) ) T Rho dV
!                    !  /                                         /
!
!                    ! We need to compute \int{T div(vfrac u_nl) T},
!                    ! so split up the div using the product rule and compute
!                    ! \int{T vfrac div(u_nl) T} + \int{T u_nl grad(vfrac) T}
!                    do i = 1, NGI
!                        u_nl_dot_grad_nvfrac_gi(i) = dot_product(U_nl_q(:,i), grad_nvfrac_gi(:,i))
!                    end do
!                    Advection_mat = shape_vector_dot_dshape(u_shape, U_nl_q, du_t, detwei*Rho_q*nvfrac_gi) &
!                        + beta * (shape_shape(u_shape, u_shape, U_nl_div_q*detwei*Rho_q*nvfrac_gi) + &
!                        shape_shape(u_shape, u_shape, detwei*Rho_q*u_nl_dot_grad_nvfrac_gi))
!                else
                    ! Element advection matrix
                    !  /                                   /
                    !  | T (U_nl dot grad T) Rho dV + beta | T ( div U_nl ) T Rho dV
                    !  /                                   /
                    Advection_mat = shape_vector_dot_dshape(u_shape, U_nl_q, du_t, detwei*Rho_q) &
                        + beta * shape_shape(u_shape, u_shape, U_nl_div_q * detwei*Rho_q)
!                end if

                if(move_mesh) then
                    Advection_mat = Advection_mat &
                        - shape_shape(u_shape, u_shape, ele_div_at_quad(U_mesh, ele, dug_t) * detwei * Rho_q)
                end if
            end if
        end if

        do dim = 1, NDIM
            if(subcycle) then
                subcycle_m_tensor_addto(dim, dim, :NLOC, :NLOC) &
                    &= subcycle_m_tensor_addto(dim, dim, :NLOC, :NLOC) &
                    &+ advection_mat
            else
                big_m_tensor_addto(dim, dim, :NLOC, :NLOC) &
                    &= big_m_tensor_addto(dim, dim, :NLOC, :NLOC) &
                    &+ dt*theta*advection_mat
            end if
            if(acceleration.and..not.subcycle) then
                rhs_addto(dim, :NLOC) = rhs_addto(dim, :NLOC) - matmul(advection_mat, u_val(dim,:))
            end if
        end do

    end if

    if(have_source.and.acceleration.and.assemble_element) then
        ! Momentum source matrix.
        Source_mat = shape_shape(U_shape, ele_shape(Source,ele), detwei*Rho_q)
        if(lump_source) then
            source_lump = sum(source_mat, 2)
            do dim = 1, NDIM
                ! lumped source
                rhs_addto(dim, :NLOC) = rhs_addto(dim, :NLOC) + source_lump*(ele_val(source, dim, ele))
            end do
        else
            do dim = 1, NDIM
                ! nonlumped source
                rhs_addto(dim, :NLOC) = rhs_addto(dim, :NLOC) + matmul(source_mat, ele_val(source, dim, ele))
            end do
        end if
    end if

    if(have_gravity.and.acceleration.and.assemble_element) then
        ! buoyancy
        if(subtract_out_reference_profile) then
            ! coefficient_detwei = detwei*gravity_magnitude*&
            !   (ele_val_at_quad(buoyancy, ele)-ele_val_at_quad(hb_density, ele))
            coefficient_detwei = detwei*gravity_magnitude* &
                (matmul(buoyancy%val(x_ele), u_shape%n) &
                    - matmul(hb_density%val(x_ele), u_shape%n))
        else
            ! coefficient_detwei = detwei*gravity_magnitude*ele_val_at_quad(buoyancy, ele)
            coefficient_detwei = detwei*gravity_magnitude* &
                matmul(buoyancy%val(x_ele), u_shape%n)
        end if

        if (on_sphere) then
            ! If were on a spherical Earth evaluate the direction of the gravity vector
            ! exactly at quadrature points.
            rhs_addto(:, :NLOC) = rhs_addto(:, :NLOC) + shape_vector_rhs(u_shape, &
                sphere_inward_normal_at_quad_ele(X, ele), &
                coefficient_detwei)
        else

!            if(multiphase) then
!                rhs_addto(:, :NLOC) = rhs_addto(:, :NLOC) + shape_vector_rhs(u_shape, &
!                    ele_val_at_quad(gravity, ele), &
!                    coefficient_detwei*nvfrac_gi)
!            else
                rhs_addto(:, :NLOC) = rhs_addto(:, :NLOC) + shape_vector_rhs(u_shape, &
                    ele_val_at_quad(gravity, ele), &
                    coefficient_detwei)
!            end if

        end if
    end if

    if((have_absorption.or.have_vertical_stabilization.or.have_wd_abs .or. have_swe_bottom_drag) .and. &
        (assemble_element .or. pressure_corrected_absorption)) then

        absorption_gi=0.0
        tensor_absorption_gi=0.0
        absorption_gi = ele_val_at_quad(Abs, ele)
        if (on_sphere.and.have_absorption) then ! Rotate the absorption
            tensor_absorption_gi=rotate_diagonal_to_sphere_gi(X, ele, absorption_gi)
        end if

        vvr_abs_diag=0.0
        vvr_abs=0.0
        ib_abs=0.0
        ib_abs_diag=0.0

        if (have_vertical_velocity_relaxation) then

            ! Form the vertical velocity relaxation absorption term
            ! depth_at_quads=ele_val_at_quad(depth, ele)
            depth_at_quads=matmul(depth%val(x_ele), u_shape%n)
            if (.not.on_sphere) then
                ! grav_at_quads=ele_val_at_quad(gravity, ele)

                if(gravity%field_type==FIELD_TYPE_CONSTANT) then
                    do concurrent (iloc=1:NLOC)
                        ele_grav_val(:,iloc) = gravity%val(:,1)
                    end do
                    grav_at_quads=matmul(ele_grav_val, u_shape%n)
                else
                    ! FIELD_TYPE_NORMAL
                    grav_at_quads=matmul(gravity%val(:,x_ele), u_shape%n)
                end if

                do concurrent (i=1:NGI)
                    vvr_abs_diag(:,i)=vvr_sf*gravity_magnitude*dt*grav_at_quads(:,i)*rho_q(i)/depth_at_quads(i)
                end do
            else
                ! on_sphere
                do concurrent(i=1:NGI)
                    vvr_abs_diag(3,i)=-vvr_sf*gravity_magnitude*dt*rho_q(i)/depth_at_quads(i)
                end do
                vvr_abs=rotate_diagonal_to_sphere_gi(X, ele, vvr_abs_diag)
            end if

        end if

        if (have_implicit_buoyancy) then

            call transform_to_physical(X, ele, ele_shape(buoyancy,ele), dshape=dt_rho)
            grad_rho=ele_grad_at_quad(buoyancy, ele, dt_rho)

            ! Calculate the gradient in the direction of gravity
            if (on_sphere) then
                grav_at_quads=sphere_inward_normal_at_quad_ele(X, ele)
            else
!                grav_at_quads=ele_val_at_quad(gravity, ele)
                if(gravity%field_type==FIELD_TYPE_CONSTANT) then
                        do concurrent (iloc=1:NLOC)
                            ele_grav_val(:,iloc) = gravity%val(:,1)
                        end do
                        grav_at_quads=matmul(ele_grav_val, u_shape%n)
                else
                        ! FIELD_TYPE_NORMAL
                        grav_at_quads=matmul(gravity%val(:,x_ele), u_shape%n)
                end if
           end if

            do i=1, NGI
                drho_dz(i)=dot_product(grad_rho(:,i),grav_at_quads(:,i)) ! Divide this by rho_0 for non-Boussinesq?
                if (drho_dz(i) < ib_min_grad) drho_dz(i)=ib_min_grad ! Default ib_min_grad=0.0
            end do

            ! Form the implicit buoyancy absorption terms
            if (on_sphere) then
                do i=1, NGI
                    ib_abs_diag(3,i)=-theta*dt*gravity_magnitude*drho_dz(i)
                end do
                ib_abs=rotate_diagonal_to_sphere_gi(X, ele, ib_abs_diag)
            else
                do i=1, NGI
                    ib_abs_diag(:,i)=theta*dt*gravity_magnitude*drho_dz(i)*grav_at_quads(:,i)
                end do
            end if

        end if

        ! Add any vertical stabilization to the absorption term
        if (on_sphere) then
            tensor_absorption_gi=tensor_absorption_gi-vvr_abs-ib_abs
            absorption_gi=absorption_gi-vvr_abs_diag-ib_abs_diag
        else
            absorption_gi=absorption_gi-vvr_abs_diag-ib_abs_diag
        end if

        if (have_swe_bottom_drag) then
            ! first compute total water depth H
            depth_at_quads = ele_val_at_quad(depth, ele) + (theta_nl*ele_val_at_quad(p, ele) + (1.0-theta_nl)*ele_val_at_quad(old_pressure, ele))/gravity_magnitude
            ! now reuse depth_at_quads to be the absorption coefficient: C_D*|u|/H
            depth_at_quads = (ele_val_at_quad(swe_bottom_drag, ele)*sqrt(sum(ele_val_at_quad(swe_u_nl, ele)**2, dim=1)))/depth_at_quads
            do i=1, NDIM
                absorption_gi(i,:) = absorption_gi(i,:) + depth_at_quads
            end do

        end if

        ! If on the sphere then use 'tensor' absorption. Note that using tensor absorption means that, currently,
        ! the absorption cannot be used in the pressure correction.
        if (on_sphere) then

            Abs_mat_sphere = shape_shape_tensor(U_shape, U_shape, detwei*rho_q, tensor_absorption_gi)
            Abs_mat = shape_shape_vector(U_shape, U_shape, detwei*rho_q, absorption_gi)
            if (have_wd_abs) then
                FLExit("Wetting and drying absorption does currently not work on the sphere.")
            end if

            if(lump_abs) then

                Abs_lump_sphere = sum(Abs_mat_sphere, 4)
                if (assemble_element) then
                    do dim = 1, NDIM
                        do dim2 = 1, NDIM
                            do i = 1, NLOC
                                big_m_tensor_addto(dim, dim2, i, i) = big_m_tensor_addto(dim, dim2, i, i) + &
                                    & dt*theta*Abs_lump_sphere(dim,dim2,i)
                            end do
                        end do
                        if (acceleration) then
                            rhs_addto(dim, :NLOC) = rhs_addto(dim, :NLOC) - Abs_lump_sphere(dim,dim,:)*u_val(dim,:)
                            ! off block diagonal absorption terms
                            do dim2 = 1, NDIM
                                if (dim==dim2) cycle ! The dim=dim2 terms were done above
                                rhs_addto(dim, :NLOC) = rhs_addto(dim, :NLOC) - Abs_lump_sphere(dim,dim2,:)*u_val(dim2,:)
                            end do
                        end if
                    end do
                end if

                if (present(inverse_masslump) .and. pressure_corrected_absorption) then
                    ! assert(lump_mass)
                    abs_lump = sum(Abs_mat, 3)
                    if(have_mass) then
                        do dim = 1, NDIM
                            call set( inverse_masslump, dim, u_ele, &
                                1.0/(l_masslump+dt*theta*abs_lump(dim,:)) )
                        end do
                    else
                        do dim = 1, NDIM
                            call set( inverse_masslump, dim, u_ele, &
                                1.0/(dt*theta*abs_lump(dim,:)) )
                        end do
                    end if
                end if

            else

                if (assemble_element) then
                    do dim = 1, NDIM
                        do dim2 = 1, NDIM
                            big_m_tensor_addto(dim, dim2, :NLOC, :NLOC) = &
                                big_m_tensor_addto(dim, dim2, :NLOC, :NLOC) + &
                                & dt*theta*Abs_mat_sphere(dim,dim2,:,:)
                        end do
                        if (acceleration) then
                            rhs_addto(dim, :NLOC) = rhs_addto(dim, :NLOC) - matmul(Abs_mat_sphere(dim,dim,:,:), u_val(dim,:))
                            ! off block diagonal absorption terms
                            do dim2 = 1, NDIM
                                if (dim==dim2) cycle ! The dim=dim2 terms were done above
                                rhs_addto(dim, :NLOC) = rhs_addto(dim, :NLOC) - matmul(Abs_mat_sphere(dim,dim2,:,:), u_val(dim2,:))
                            end do
                        end if
                    end do
                end if
                Abs_lump_sphere = 0.0
                if (present(inverse_mass) .and. pressure_corrected_absorption) then
                    ! assert(.not. lump_mass)
                    if(have_mass) then
                        do dim = 1, NDIM
                            call set(inverse_mass, dim, dim, u_ele, u_ele, &
                                inverse(rho_mat + dt*theta*Abs_mat(dim,:,:)))
                        end do
                    else
                        do dim = 1, NDIM
                            call set(inverse_mass, dim, dim, u_ele, u_ele, &
                                inverse(dt*theta*Abs_mat(dim,:,:)))
                        end do
                    end if
                end if

            end if

        else

            Abs_mat = shape_shape_vector(U_shape, U_shape, detwei*rho_q, absorption_gi)

            if (have_wd_abs) then
                !! Wetting and drying absorption becomes active when water level reaches d_0
                alpha_u_quad=ele_val_at_quad(alpha_u_field, ele)
                Abs_mat = Abs_mat + shape_shape_vector(U_shape, U_shape, alpha_u_quad*detwei*rho_q, &
                    &                                 ele_val_at_quad(Abs_wd,ele))
            end if

            if(lump_abs) then
                abs_lump = sum(Abs_mat, 3)

                if (assemble_element) then
                    do dim = 1, NDIM
                        big_m_diag_addto(dim, :NLOC) = big_m_diag_addto(dim, :NLOC) + dt*theta*abs_lump(dim,:)
                    end do
                    if(acceleration) then
                        do dim = 1, NDIM
                            rhs_addto(dim, :NLOC) = rhs_addto(dim, :NLOC) - abs_lump(dim,:)*u_val(dim,:)
                        end do
                    end if
                end if
                if (present(inverse_masslump) .and. pressure_corrected_absorption) then
                    !assert(lump_mass)
                    if(have_mass) then
                        do dim = 1, NDIM
                            call set( inverse_masslump, dim, u_ele, &
                                1.0/(l_masslump+dt*theta*abs_lump(dim,:)) )
                        end do
                    else
                        do dim = 1, NDIM
                            call set( inverse_masslump, dim, u_ele, &
                                1.0/(dt*theta*abs_lump(dim,:)) )
                        end do
                    end if
                end if

            else

                if (assemble_element) then
                    do dim = 1, NDIM
                        big_m_tensor_addto(dim, dim, :NLOC, :NLOC) = big_m_tensor_addto(dim, dim, :NLOC, :NLOC) + &
                            & dt*theta*Abs_mat(dim,:,:)
                    end do
                    if(acceleration) then
                        do dim = 1, NDIM
                            rhs_addto(dim, :NLOC) = rhs_addto(dim, :NLOC) - matmul(Abs_mat(dim,:,:), u_val(dim,:))
                        end do
                    end if
                end if
                if (present(inverse_mass) .and. pressure_corrected_absorption) then
                    !assert(.not. lump_mass)
                    if(have_mass) then
                        do dim = 1, NDIM
                            call set(inverse_mass, dim, dim, u_ele, u_ele, &
                                inverse(rho_mat + dt*theta*Abs_mat(dim,:,:)))
                        end do
                    else
                        do dim = 1, NDIM
                            call set(inverse_mass, dim, dim, u_ele, u_ele, &
                                inverse(dt*theta*Abs_mat(dim,:,:)))
                        end do
                    end if
                end if

            end if

        end if

    end if

    if ((((.not.have_absorption).and.(.not.have_vertical_stabilization).and.(.not.have_wd_abs)) .or. (.not.pressure_corrected_absorption)).and.(have_mass)) then
        ! no absorption: all mass matrix components are the same
        if (present(inverse_mass) .and. .not. lump_mass) then
            inverse_mass_mat=inverse(rho_mat)
            call set(inverse_mass, 1, 1, u_ele, u_ele, inverse_mass_mat)
            if (.not. inverse_mass%equal_diagonal_blocks) then
                ! after the strong dirichlet bcs have been applied, the diagonal
                ! blocks will be different. So for now we just copy:
                do dim = 2, NDIM
                    call set(inverse_mass, dim, dim, u_ele, u_ele, inverse_mass_mat)
                end do
            end if
        end if
        if (present(inverse_masslump) .and. lump_mass) then
            do dim = 1, NDIM
                call set(inverse_masslump, dim, u_ele, 1.0/l_masslump)
            end do
        end if

    end if

    ! Viscosity.
    Viscosity_mat=0
    if(assemble_element) then

        if (primal) then

            if(partial_stress) then
                ! This is where to stick the partial stress stuff for LES
                ! - Angus

                ! Highly experimental at the moment..
                ! - Angus (26/03/2016)
                do gi=1, NGI
                    do iloc=1, NLOC
                        node_stress_diag = matmul(du_t(iloc,gi,:), visc_ele_quad(:,:, gi))
                        visc_grad_dot_u = visc_ele_quad(:,:, gi) * sum(du_t(iloc,gi, :))

                        do jloc=1, NLOC
                            visc_dot_prod = dot_product(node_stress_diag, du_t(jloc,gi,:))*detwei(gi)
                            resid_stress_term = matmul(visc_grad_dot_u, du_t(jloc,gi,:))*detwei(gi)

                            do concurrent(idim=1:NDIM)
                                Viscosity_mat(idim,idim,iloc,jloc) = &
                                    Viscosity_mat(idim,idim,iloc,jloc) &
                                    + visc_dot_prod + resid_stress_term(idim)
                            end do
                        end do
                    end do
                end do


            else
                ! Tensor form
                do gi=1, NGI
                    do iloc=1, NLOC
                        node_stress_diag = matmul(du_t(iloc,gi,:), visc_ele_quad(:,:, gi))

                        ! Is this faster than that below?
                        do jloc=1, NLOC
                            visc_dot_prod = dot_product(node_stress_diag, du_t(jloc,gi,:))*detwei(gi)

                            do concurrent(idim=1:NDIM)
                                Viscosity_mat(idim,idim,iloc,jloc) = &
                                    Viscosity_mat(idim,idim,iloc,jloc) &
                                    + visc_dot_prod
                            end do
                        end do

                        !   Or is this better optimised?
!                        do concurrent(idim=1:NDIM,jloc=1:NLOC)
!                                Viscosity_mat(idim,idim,iloc,jloc) = &
!                                    Viscosity_mat(idim,idim,iloc,jloc) &
!                                    + dot_product(node_stress_diag, du_t(jloc,gi,:))*detwei(gi)
!                        end do

                    end do
                end do

            end if


#if defined (SCHEME_CDG) || defined (SCHEME_IP)
                !Compute a matrix which maps ele vals to ele grad vals
                !This works since the gradient of the shape function
                !lives in the original polynomial space -- cjc

                Mass_mat = shape_shape(u_shape, u_shape, detwei)
                inverse_mass_mat = mass_mat
                call invert(inverse_mass_mat)
#endif

            ! Get kappa mat for CDG

#if defined (SCHEME_CDG)
        !            if(viscosity_scheme==CDG) then
!                if(multiphase) then
!                    ! kappa = mu*vfrac for multiphase
!                    kappa_mat = shape_shape_tensor(u_shape,u_shape,detwei*nvfrac_gi, &
!                        & visc_ele_quad)
!                else
                    kappa_mat = shape_shape_tensor(u_shape,u_shape,detwei, &
                        & visc_ele_quad )
!                end if
        !            end if
#endif

        else
            ! Tau Q = grad(u)
!            if(multiphase) then
!                ! We define the auxiliary variable as vfrac*q = vfrac*div(u)
!                ! to obtain the correct form of the grad_u_mat_q matrix. This way,
!                ! transpose(grad_u_mat_q) gives the correct form of the viscosity term.
!                Q_inv= shape_shape(q_shape, q_shape, detwei*nvfrac_gi)
!            else
                Q_inv= shape_shape(q_shape, q_shape, detwei)
!            end if

            call invert(Q_inv)
            call cholesky_factor(Q_inv)

            Grad_U_mat_q=0.0
            Div_U_mat_q=0.0
            if(.not.p0) then

!                if(multiphase) then
!                    ! Split up -\int{grad(N_A vfrac) N_B} using the product rule
!                    ! and compute -\int{grad(N_A) vfrac N_B} - \int{N_A grad(vfrac) N_B}
!                    Grad_U_mat_q(:, :, :NLOC) = -dshape_shape(dq_t, u_shape, detwei*nvfrac_gi) - &
!                        & shape_shape_vector(q_shape, u_shape, detwei, grad_nvfrac_gi)
!                else
                    Grad_U_mat_q(:, :, :NLOC) = -dshape_shape(dq_t, u_shape, detwei)
!                end if

            end if
        end if
    end if

    if(have_surfacetension.and.(.not.p0).and.assemble_element) then
        if(integrate_surfacetension_by_parts) then
            tension = ele_val_at_quad(surfacetension, ele)

            rhs_addto(:,:NLOC) = rhs_addto(:,:NLOC) - &
                &dshape_dot_tensor_rhs(du_t, tension, detwei)
        else
            dtensiondj = ele_div_at_quad_tensor(surfacetension, ele, du_t)

            rhs_addto(:,:NLOC) = rhs_addto(:,:NLOC) + &
                & shape_vector_rhs(u_shape,dtensiondj,detwei)
        end if
    end if

    !-------------------------------------------------------------------
    ! Interface integrals
    !-------------------------------------------------------------------

    if(assemble_element) then
        neigh=>ele_neigh(U, ele)

        ! x_neigh/=t_neigh only on periodic boundaries.
        x_neigh=>ele_neigh(X, ele)
        ! x_neigh = neigh


        ! Local node map counter.
        start=NLOC+1
        ! Flag for whether this is a boundary element.
        boundary_element=.false.

        ! Stuff that can be used with faces, but only called once.

        ! Face shapes are the same, so just set them here.
        face_u_shape=face_shape(U, 1)
        face_p_shape=face_shape(P, 1)
        face_q_shape=face_shape(q_mesh, 1)
       
        ! If constant density, just do this once 
        if(Rho%field_type==FIELD_TYPE_CONSTANT) then
           face_Rho_val = Rho%val(1)
           face_Rho_q = matmul(face_Rho_val, face_u_shape%n)
        end if

        ! If constant viscosity, just do this once 
        if(Viscosity%field_TYPE==FIELD_TYPE_CONSTANT) then
           do concurrent (iloc=1:FLOC)
              face_visc_val(:,:,iloc) = Viscosity%val(:,:,1)
           end do
           face_kappa_gi = tensormul( face_visc_val, &
                face_u_shape%n )
        end if

        
        neighbourloop: do ni=1, NFACES

            !----------------------------------------------------------------------
            ! Find the relevant faces.
            !----------------------------------------------------------------------
            turbine_face=.false.
            ! These finding routines are outside the inner loop so as to allow
            ! for local stack variables of the right size in
            ! construct_momentum_interface_dg.

            ele_2=neigh(ni)

            ! Note that although face is calculated on field U, it is in fact
            ! applicable to any field which shares the same mesh topology.
            face=ele_face(U, ele, ele_2)

            if (ele_2>0) then
                ! Internal faces.
                face_2=ele_face(U, ele_2, ele)
               ! Check if face is turbine face (note: get_entire_boundary_condition only returns "applied" boundaries and we reset the apply status in each timestep)
            !            elseif (velocity_bc_type(1,face)==4 .or. velocity_bc_type(1,face)==5) then
            !                face_2=face_neigh(turbine_conn_mesh, face)
            !                turbine_face=.true.
            else
                ! External face.
                face_2=face
                boundary_element=.true.
            end if

            !Compute distance between cell centre and neighbouring cell centre
            !This is for Interior Penalty Method -- cjc
            !--------------
#if defined (SCHEME_IP)
            !            if(viscosity_scheme==IP) then

                if(edge_length_option==USE_ELEMENT_CENTRES) then
                    ele_2_X = x_neigh(ni)
                    ele_centre = sum(X_val,2)/size(X_val,2)
                    face_centre = sum(face_val(X,face),2)/size(face_val(X,face),2)

                    if(face==face_2) then
                        ! Boundary case. We compute 2x the distance to the face centre
                        h0 = 2*sqrt( sum(ele_centre - face_centre)**2 )
                    else if (ele_2/=x_neigh(ni)) then
                        ! Periodic boundary case. We have to cook up the coordinate by
                        ! adding vectors to the face from each side.
                        x_val_2 = ele_val(X,ele_2_X)
                        neigh_centre = sum(X_val_2,2)/size(X_val_2,2)
                        face_centre_2 = &
                            sum(face_val(X,face_2),2)/size(face_val(X,face_2),2)
                        h0 = sqrt ( sum(ele_centre - face_centre)**2 )
                        h0 = h0 + sqrt( sum(neigh_centre - face_centre_2)**2 )
                    else
                        x_val_2 = ele_val(X,ele_2_X)
                        neigh_centre = sum(X_val_2,2)/size(X_val_2,2)
                        h0 = sqrt ( sum(ele_centre - neigh_centre)**2 )
                    end if
                end if
            !            end if
#endif
            !--------------

            finish=start+FLOC-1

            local_glno(start:finish)=face_global_nodes(U, face_2)

                        ! Always primal now, and no turbines
            !                call construct_momentum_interface_dg_ELEMENT_CONFIG(ele, face, face_2, ni,&
            !                    & big_m_tensor_addto, &
            !                    & rhs_addto, Grad_U_mat_q, Div_U_mat_q, X,&
            !                    & Rho, U, U_nl, U_mesh, P, q_mesh, &
            !                    & surfacetension, &
            !                    & velocity_bc, velocity_bc_type, &
            !                    & pressure_bc, pressure_bc_type, hb_pressure, &
            !                    & subcycle_m_tensor_addto, nvfrac, &
            !                    & ele2grad_mat=ele2grad_mat, kappa_mat=kappa_mat, &
            !                    & inverse_mass_mat=inverse_mass_mat, &
            !                    & viscosity=viscosity, viscosity_mat=viscosity_mat, &
            !                    & tensor_eddy_visc=tensor_eddy_visc)



            ! ********** START OF INLINED CODE **********

            face_start=NLOC+(ni-1)*FLOC+1
            face_finish=face_start+FLOC-1


            ! u_face_l=face_local_nodes(U, face)
            u_face_l = U%mesh%faces%face_lno( FLOC*(face-1)+1 : FLOC*face )
            if(move_mesh) u_mesh_glno = face_global_nodes(U_mesh, face)

            u_face_glno_1 = face_global_nodes(U_nl, face)
            u_face_glno_2 = face_global_nodes(U_nl, face_2)
            rho_face_glno_1 = face_global_nodes(Rho, face)

            ! face_u_shape_2=>face_shape(U, face_2)

            q_face_l = q_mesh%faces%face_lno( FLOC*(face-1)+1 : FLOC*face )

            ! Get Density and (non-linear) PhaseVolumeFraction values
            ! at the Gauss points on the current face.
            ! face_Rho_q=face_val_at_quad(Rho, face)
            ! Only called if density varies

            if(Rho%field_type==FIELD_TYPE_NORMAL) then
               face_Rho_q = matmul(Rho%val(Rho_face_glno_1), face_u_shape%n)
            end if

            !    if(multiphase) then
            !        face_nvfrac_gi = face_val_at_quad(nvfrac, face)
            !    end if

            ! allocate( face_kappa_gi(Viscosity%dim(1), Viscosity%dim(2), &
            !    face_ngi(Viscosity,face)) )

! Only for Primal schemes (CDG, IP, etc.). Not needed for Bassi-Rebay
! (Am I sure about this?)
#ifndef SCHEME_BASSI
            ! If viscosity varies, calculate matrix at every iteration
            if(Viscosity%field_TYPE==FIELD_TYPE_NORMAL) then
               face_kappa_gi = tensormul( Viscosity%val(:,:,u_face_glno_1), &
                    face_u_shape%n )
            end if

            ! For Large-Eddy Simulation

!            if(associated(tensor_eddy_visc)) then
!                face_kappa_gi = face_kappa_gi + &
!                     tensormul( tensor_eddy_visc%val(:,:,u_face_glno_1), &
!                     face_u_shape%n )
!            end if

            if(have_les) then
                if(have_isotropic_les) then
                    tmp_face_tensor=0.
                    do concurrent(face_d1=1:NDIM)
                        tmp_face_tensor(face_d1, face_d1,:) = eddy_visc%val(u_face_glno_1)
                    end do
                    face_kappa_gi = face_kappa_gi + &
                         tensormul( tmp_face_tensor, &
                         face_u_shape%n )
                else
                    face_kappa_gi = face_kappa_gi + &
                         tensormul( tensor_eddy_visc%val(:,:,u_face_glno_1), &
                         face_u_shape%n )
                end if
             end if
#endif

!            if(multiphase) then
!                ! Multiply the viscosity tensor by the PhaseVolumeFraction
!                ! since kappa = viscosity*vfrac for multiphase flow simulations
!                face_nvfrac_gi = face_val_at_quad(nvfrac, face)
!                do face_d1=1, NDIM
!                    do face_d2=1, NDIM
!                        face_kappa_gi(face_d1,face_d2,:) = face_kappa_gi(face_d1,face_d2,:)*face_nvfrac_gi
!                    end do
!                end do
!            end if

            ! Boundary nodes have both faces the same.
            face_boundary=(face==face_2)
            face_dirichlet=.false.
            face_free_surface=.false.
            face_no_normal_flow=.false.
            face_l_have_pressure_bc=.false.
            if (face_boundary) then
                do face_dim=1, NDIM
                    if (velocity_bc_type(face_dim,face)==1) then
                        face_dirichlet(face_dim)=.true.
                    end if
                end do
                ! free surface b.c. is set for the 1st (normal) component
                if (velocity_bc_type(1,face)==2) then
                    face_free_surface=.true.
                end if
                ! no normal flow b.c. is set for the 1st (normal) component
                if (velocity_bc_type(1,face)==3) then
                    ! No normal flow is implemented here by switching off the
                    ! advection boundary integral.
                    face_no_normal_flow=.true.
                end if
                face_l_have_pressure_bc = pressure_bc_type(face) > 0
            end if

            !----------------------------------------------------------------------
            ! Change of coordinates on face.
            !----------------------------------------------------------------------
            call transform_facet_to_physical(X, face,&
                &                          detwei_f=face_detwei,&
                &                          normal=face_normal)

            !----------------------------------------------------------------------
            ! Construct bilinear forms.
            !----------------------------------------------------------------------

            if(have_advection.and..not.face_no_normal_flow) then
                ! Advecting velocity at quadrature points.

                !print*, "face_global_nodes:", face_global_nodes(U,face)
                !print*, "u_face_l:", u_face_l
                ! shape=>face_shape(field,face_number)
                ! quad_val=matmul(face_val(field, face_number), shape%n)

                ! face_u_f_q = face_val_at_quad(U_nl, face)
                face_u_f_q = matmul(U_nl%val(:,u_face_glno_1), face_u_shape%n)
                ! face_u_f2_q = face_val_at_quad(U_nl, face_2)
                face_u_f2_q = matmul(U_nl%val(:,u_face_glno_2), face_u_shape%n)

                face_u_nl_q=0.5*(face_u_f_q+face_u_f2_q)

                face_div_u_f_q = face_u_f_q


                ! Mesh velocity at quadrature points.
                if(move_mesh) then
                    ! here we assume that U_mesh at face is the same as U_mesh at face_2
                    ! if it isn't then you're in trouble because your mesh will tear
                    ! itself apart
                    face_u_mesh_quad = matmul( U_mesh%val(:, u_mesh_glno), &
                         face_u_shape%n )
                    ! face_u_nl_q=face_u_nl_q - face_val_at_quad(U_mesh, face)
                    face_u_nl_q=face_u_nl_q - face_u_mesh_quad
                    ! the velocity on the internal face isn't used again so 
                    ! we can modify it directly here...
                    face_u_f_q = face_u_f_q - face_u_mesh_quad
                end if

                face_u_nl_q_dotn = sum(face_u_nl_q*face_normal,1)

                ! Inflow is true if the flow at this gauss point is directed
                ! into this element.
                face_inflow= face_u_nl_q_dotn<0.0
                face_income = merge(1.0,0.0,face_inflow)

                ! Calculate outflow boundary integral.
                ! can anyone think of a way of optimising this more to avoid
                ! superfluous operations (i.e. multiplying things by 0 or 1)?

                ! first the integral around the inside of the element
                ! (this is the flux *out* of the element)
                face_inner_advection_integral = (1.-face_income)*face_u_nl_q_dotn
                if(.not.integrate_by_parts_once) then
                    ! i.e. if we're integrating by parts twice
                    face_inner_advection_integral = face_inner_advection_integral &
                        - sum(face_u_f_q*face_normal,1)
                end if
                if(integrate_conservation_term_by_parts) then
                    if(integrate_by_parts_once) then
                        face_inner_advection_integral = face_inner_advection_integral &
                            - (1.-beta)*sum(face_div_u_f_q*face_normal,1)
                    else
                        ! i.e. integrating by parts twice
                        face_inner_advection_integral = face_inner_advection_integral &
                            + beta*sum(face_div_u_f_q*face_normal,1)
                    end if
                end if

!                if(multiphase) then
!                    face_nnAdvection_out=shape_shape(face_u_shape, face_u_shape,  &
!                        &                        face_inner_advection_integral * face_detwei * face_Rho_q * face_nvfrac_gi)
!                else
                    face_nnAdvection_out=shape_shape(face_u_shape, face_u_shape,  &
                        &                        face_inner_advection_integral * face_detwei * face_Rho_q)
!                end if

                ! now the integral around the outside of the element
                ! (this is the flux *in* to the element)
                face_outer_advection_integral = face_income * face_u_nl_q_dotn
!                if(multiphase) then
!                    !            face_nnAdvection_in=shape_shape(face_u_shape, face_u_shape_2, &
!                    !                &                       face_outer_advection_integral * face_detwei * face_Rho_q * face_nvfrac_gi)
!                    face_nnAdvection_in=shape_shape(face_u_shape, face_u_shape, &
!                        &                       face_outer_advection_integral * face_detwei * face_Rho_q * face_nvfrac_gi)
!                else
                    !            face_nnAdvection_in=shape_shape(face_u_shape, face_u_shape_2, &
                    !                &                       face_outer_advection_integral * face_detwei * face_Rho_q)
                    face_nnAdvection_in=shape_shape(face_u_shape, face_u_shape, &
                        &                       face_outer_advection_integral * face_detwei * face_Rho_q)
!                end if

                do face_dim = 1, NDIM

                    ! Insert advection in matrix.

                    ! Outflow boundary integral.
                    if(subcycle) then
                        subcycle_m_tensor_addto(face_dim, face_dim, u_face_l, u_face_l) = &
                            &subcycle_m_tensor_addto(face_dim, face_dim, u_face_l, u_face_l) + &
                            &face_nnAdvection_out

                        if (.not.face_dirichlet(face_dim)) then
                            subcycle_m_tensor_addto(face_dim, face_dim, u_face_l, face_start:face_finish) = &
                                &subcycle_m_tensor_addto(face_dim, face_dim, u_face_l, face_start:face_finish)&
                                &+face_nnAdvection_in
                        end if
                    else
                        big_m_tensor_addto(face_dim, face_dim, u_face_l, u_face_l) = &
                            big_m_tensor_addto(face_dim, face_dim, u_face_l, u_face_l) + &
                            face_nnAdvection_out*dt*theta

                        if (.not.face_dirichlet(face_dim)) then
                            big_m_tensor_addto(face_dim, face_dim, u_face_l, face_start:face_finish) = &
                                big_m_tensor_addto(face_dim, face_dim, u_face_l, face_start:face_finish) + &
                                face_nnAdvection_in*dt*theta
                        end if

                        if (.not.face_dirichlet(face_dim)) then
                            ! For interior interfaces this is the upwinding term. For a
                            ! Neumann boundary it's necessary to apply downwinding here
                            ! to maintain the surface integral. Fortunately, since
                            ! face_2==face for a boundary this is automagic.

                            if (acceleration) then
                                rhs_addto(face_dim,u_face_l) = rhs_addto(face_dim,u_face_l) &
                                    ! Outflow boundary integral.
                                    -matmul(face_nnAdvection_out,face_val(U,face_dim,face))&
                                    ! Inflow boundary integral.
                                    -matmul(face_nnAdvection_in,face_val(U,face_dim,face_2))
                            end if

                        else

                            rhs_addto(face_dim,u_face_l) = rhs_addto(face_dim,u_face_l) &
                                ! Outflow boundary integral.
                                -matmul(face_nnAdvection_out,face_val(U,face_dim,face))&
                                ! Inflow boundary integral.
                                -matmul(face_nnAdvection_in,ele_val(velocity_bc,face_dim,face))
                        end if
                    end if
                end do

            end if

            ! Boundary term in grad_U.
            !   /
            !   | q, u, normal dx
            !   /

#if defined (SCHEME_IP)
            !    select case (viscosity_scheme)
            !        case (IP)
            face_primal_fluxes_mat = 0.0
            face_penalty_fluxes_mat = 0.0
            call primal_fluxes
            call interior_penalty
            call local_assembly_primal_face
            call local_assembly_ip_face
#elif defined (SCHEME_CDG)
            !        case (CDG)
            face_primal_fluxes_mat = 0.0
            face_penalty_fluxes_mat = 0.0
            call primal_fluxes
            if(.not.remove_penalty_fluxes) call interior_penalty
            call get_normal_mat
            call local_assembly_primal_face
            call local_assembly_cdg_face
            call local_assembly_ip_face
#elif defined (SCHEME_BASSI)
            call bassi_rebay_viscosity
#endif

            !    end select

            if(have_surfacetension.and.integrate_surfacetension_by_parts) then
                face_tension_q = 0.5*face_val_at_quad(surfacetension,face)+0.5*face_val_at_quad(surfacetension,face_2)
                rhs_addto(:,u_face_l) = rhs_addto(:,u_face_l) + shape_tensor_dot_vector_rhs(face_u_shape, face_tension_q, face_normal, face_detwei)
            end if


            !----------------------------------------------------------------------
            ! Perform global assembly.
            !----------------------------------------------------------------------

            ! Insert pressure boundary integral.
            if (l_include_pressure_bcs .and. face_boundary .and. face_l_have_pressure_bc) then

!                if(multiphase) then
!                    !            face_mnCT(1,:,:,:) = shape_shape_vector(face_p_shape, face_u_shape_2, face_detwei*face_nvfrac_gi, face_normal)
!                    face_mnCT(1,:,:,:) = shape_shape_vector(face_p_shape, face_u_shape, face_detwei*face_nvfrac_gi, face_normal)
!                else
                    !            face_mnCT(1,:,:,:) = shape_shape_vector(face_p_shape, face_u_shape_2, face_detwei, face_normal)
                    face_mnCT(1,:,:,:) = shape_shape_vector(face_p_shape, face_u_shape, face_detwei, face_normal)
!                end if
                ! for both weak and strong pressure dirichlet bcs:
                !      /
                ! add -|  N_i M_j \vec n p_j, where p_j are the prescribed bc values
                !      /
                if(subtract_out_reference_profile) then
                    do face_dim = 1, NDIM
                        rhs_addto(face_dim,u_face_l) = rhs_addto(face_dim,u_face_l) - &
                            matmul( ele_val(pressure_bc, face) - face_val(hb_pressure, face), face_mnCT(1,face_dim,:,:) )
                    end do
                else
                    do face_dim = 1, NDIM
                        rhs_addto(face_dim,u_face_l) = rhs_addto(face_dim,u_face_l) - &
                            matmul( ele_val(pressure_bc, face), face_mnCT(1,face_dim,:,:) )
                    end do
                end if
            end if

            ! ********** END OF INLINED CODE **********

            start=start+FLOC

        end do neighbourloop

        !----------------------------------------------------------------------
        ! Construct local diffusivity operator for DG.
        !----------------------------------------------------------------------

#if defined(SCHEME_BASSI)
        if (partial_stress) then
            call local_assembly_bassi_rebay_stress_form
        else
            call local_assembly_bassi_rebay
        end if
#endif

        ! if(have_viscosity) then

        if (boundary_element) then

            ! Weak application of dirichlet conditions on viscosity term.

            weak_dirichlet_loop: do i=1,2
                ! this is done in 2 passes
                ! iteration 1: wipe the rows corresponding to weak dirichlet boundary faces
                ! iteration 2: for columns corresponding to weak dirichlet boundary faces,
                !               move this coefficient multiplied with the bc value to the rhs
                !               then wipe the column
                ! The 2 iterations are necessary for elements with more than one weak dirichlet boundary face
                ! as we should not try to move the coefficient in columns corresponding to boundary face 1
                ! in rows correspoding to face 2 to the rhs, i.e. we need to wipe *all* boundary rows first.

                do dim=1, NDIM

                    ! Local node map counter.
                    start=NLOC+1

                    boundary_neighbourloop: do ni=1, NFACES
                        ele_2=neigh(ni)

                        ! Note that although face is calculated on field U, it is in fact
                        ! applicable to any field which shares the same mesh topology.
                        if (ele_2>0) then
                            ! Interior face - we need the neighbouring face to
                            ! calculate the new start
                            face=ele_face(U, ele_2, ele)
                        else
                            ! Boundary face
                            face=ele_face(U, ele, ele_2)
                            if (velocity_bc_type(dim,face)==1) then

                                ! Dirichlet condition.

                                finish=start+FLOC-1

                                if (i==1) then
                                    ! Wipe out boundary condition's coupling to itself.
                                    Viscosity_mat(:,dim,start:finish,:)=0.0
                                else
                                    ! Add BC into RHS
                                    !
                                    do dim1=1, NDIM
                                        rhs_addto(dim1,:) = rhs_addto(dim1,:) &
                                            & -matmul(Viscosity_mat(dim1,dim,:,start:finish), &
                                            & ele_val(velocity_bc,dim,face))
                                    end do
                                    ! Ensure it is not used again.
                                    Viscosity_mat(:,dim,:,start:finish)=0.0
                                end if
                               ! Check if face is turbine face (note: get_entire_boundary_condition only returns
                               ! "applied" boundaries and we reset the apply status in each timestep)
                            elseif (velocity_bc_type(dim,face)==4 .or. velocity_bc_type(dim,face)==5) then
                                face=face_neigh(turbine_conn_mesh, face)
                            end if
                        end if
                        start=start+FLOC

                    end do boundary_neighbourloop

                end do

            end do weak_dirichlet_loop

        end if

        ! Insert viscosity in matrix.
        big_m_tensor_addto = big_m_tensor_addto + Viscosity_mat*theta*dt

        if (acceleration) then
            do dim1=1, NDIM
                do dim2=1, NDIM
                    rhs_addto(dim1, :) = rhs_addto(dim1, :) &
                        - matmul(Viscosity_mat(dim1,dim2,:,:), &
                        node_val(U, dim2, local_glno))
                end do
            end do
        end if

        ! end if !have_viscosity

    end if !dg.and.(have_viscosity.or.have_advection)

    !----------------------------------------------------------------------
    ! Perform global assembly.
    !----------------------------------------------------------------------

    if (assemble_element) then

        ! add lumped terms to the diagonal of the matrix
        !        call add_diagonal_to_tensor(big_m_diag_addto, big_m_tensor_addto)

        do loc=1, EFLOC
            do dim=1, NDIM
                big_m_tensor_addto(dim, dim, loc, loc) &
                    = big_m_tensor_addto(dim, dim, loc, loc) + big_m_diag_addto(dim, loc)
            end do
        end do


        ! We always have DG, viscosity and advection for accelerated code.
        !        if(dg.and.(have_viscosity.or.have_advection)) then

        ! first the diagonal blocks, i.e. the coupling within the element
        ! and neighbouring face nodes but with the same component
        !            if(have_viscosity) then
        if(partial_stress) then
            call addto(big_m, local_glno, local_glno, &
                big_m_tensor_addto)
        else
            ! add to the matrix
            call addto(big_m, local_glno, local_glno, big_m_tensor_addto, &
                block_mask=diagonal_block_mask)
        end if
        ! add to the rhs
        call addto(rhs, local_glno, rhs_addto)

        !            else
        !                ! add to the matrix
        !                call addto(big_m, u_ele, local_glno, big_m_tensor_addto(:,:,:NLOC,:), &
        !                    block_mask=diagonal_block_mask)
        !                ! add to the rhs
        !                call addto(rhs, u_ele, rhs_addto(:,:NLOC))
        !            end if

        if(subcycle) then
            call addto(subcycle_m, u_ele, local_glno,&
                &subcycle_m_tensor_addto(:,:,:NLOC,:), &
                &block_mask=diagonal_block_mask)
        end if
        if(.not. partial_stress .and. have_coriolis) then
            ! add in coupling between different components, but only within the element
            call addto(big_m, u_ele, u_ele, &
                big_m_tensor_addto(:,:,:NLOC,:NLOC), block_mask&
                &=off_diagonal_block_mask)
        end if
    !        else
    !            ! in this case we only have coupling between nodes within the element
    !            if (have_coriolis) then
    !                call addto(big_m, u_ele, u_ele, big_m_tensor_addto(:,:,:NLOC,:NLOC))
    !            else
    !                ! add to the matrix
    !                call addto(big_m, u_ele, u_ele, big_m_tensor_addto(:,:,:NLOC,:NLOC), &
    !                    block_mask=diagonal_block_mask)
    !            end if
    !            ! add to the rhs
    !            call addto(rhs, u_ele, rhs_addto(:,:NLOC))
    !        end if

    end if

!contains
!
!    subroutine add_diagonal_to_tensor(big_m_diag_addto, big_m_tensor_addto)
!        real, dimension(NDIM, EFLOC), intent(in) :: big_m_diag_addto
!        real, dimension(NDIM, NDIM, EFLOC, EFLOC), intent(inout) :: big_m_tensor_addto
!
!        integer :: dim, loc
!
!        do loc=1, EFLOC
!            do dim=1, NDIM
!                big_m_tensor_addto(dim, dim, loc, loc) = big_m_tensor_addto(dim, dim, loc, loc) + big_m_diag_addto(dim, loc)
!            end do
!        end do
!
!    end subroutine add_diagonal_to_tensor

contains


    subroutine bassi_rebay_viscosity
        integer :: idim
        real, dimension(FNGI) :: coefficient_detwei

        do idim=1,NDIM

            coefficient_detwei = face_detwei*face_normal(idim,:)
! No multiphase for now.
!
!            if(multiphase) then
!                coefficient_detwei = coefficient_detwei*face_nvfrac_gi
!            end if

            if(.not. face_boundary) then
                ! Internal face.
                Grad_U_mat_q(idim, q_face_l, U_face_l)=&
                    Grad_U_mat_q(idim, q_face_l, U_face_l) &
                    +0.5*shape_shape(face_q_shape, face_U_shape, coefficient_detwei)

                ! External face.
                Grad_U_mat_q(idim, q_face_l, start:finish)=&
                    +0.5*shape_shape(face_q_shape, face_U_shape, coefficient_detwei)
            else
                ! Boundary case. Put the whole integral in the external bit.

                ! External face.
                Grad_U_mat_q(idim, q_face_l, face_start:face_finish)=&
                    +shape_shape(face_q_shape, face_U_shape, coefficient_detwei)
            end if
        end do

    end subroutine bassi_rebay_viscosity


    subroutine get_normal_mat
        !!< We assemble
        !!< \int_e N_i N_j n dS
        !!< where n is the normal
        !!< indices are (dim1, loc1, loc2)

        integer :: nmat_d1,nmat_d2

!        face_normal_mat = shape_shape_vector(face_u_shape,face_u_shape,face_detwei,face_normal)
        ! inlined version of above
        do concurrent(iloc=1:FLOC,jloc=1:FLOC)
           ! Main mass matrix.
           face_normal_mat(:,iloc,jloc)=&
                matmul(face_normal*spread(face_u_shape%n(iloc,:)*face_u_shape%n(jloc,:),1,NDIM),face_detwei)
        end do



        !!< We assemble
        !!< \int_e N_i N_j kappa.n dS
        !!< where n is the normal
        !!< indices are (dim1, loc1, loc2)

        face_kappa_normal_mat = 0
        do nmat_d1 = 1, NDIM
            do nmat_d2 = 1, NDIM
               face_detwei_work=face_detwei* &
                    & face_kappa_gi(nmat_d1,nmat_d2,:)*face_normal(nmat_d2,:)

!                face_kappa_normal_mat(nmat_d1,:,:) = face_kappa_normal_mat(nmat_d1,:,:) + &
!                    & shape_shape(face_u_shape,face_u_shape,face_detwei* &
!                    & face_kappa_gi(nmat_d1,nmat_d2,:)*face_normal(nmat_d2,:))

               ! inlined version of above
               do concurrent(iloc=1:FLOC,jloc=1:FLOC)
                  face_shape_shape_work(iloc,jloc)= &
                       dot_product(face_u_shape%n(iloc,:)*face_u_shape%n(jloc,:),face_detwei_work)
               end do

               face_kappa_normal_mat(nmat_d1,:,:) = &
                    face_kappa_normal_mat(nmat_d1,:,:) + face_shape_shape_work
            end do
        end do

    end subroutine get_normal_mat

    subroutine primal_fluxes

        !!< Notes for primal fluxes which are present in the interior penalty
        !!< and CDG methods (and, I believe, the LDG method when written in
        !! primal form)

        !!< We assemble

        !!< -Int_e [u]{kappa grad v} + [v]{kappa grad u}

        !!< = -Int_e 1/2(u^+n^+ + u^-n^-).(kappa^+ grad v^+ + kappa^- grad v^-)
        !!<  -Int_e 1/2(v^+n^+ + v^-n^-).(kappa^+ grad u^+ + kappa^- grad u^-)

        !!< Where + is the ele side, and - is the ele_2 side, and e is the edge

        !!<Computing grad u (and v) requires a element transform to physical
        !!<so we only assemble the + parts here, and the minus parts of the grad
        !!<will be assembled when we visit that element

        !!<So we assemble

        !!<  -Int_e 1/2 (u^+ - u^-)n^+.kappa^+ grad v^+
        !!< -Int_e 1/2 (v^+ - v^-)n^+.kappa^+ grad u^+

        !!<Actually we won't even do that, we'll just assemble the second
        !!<line, and apply the transpose operator

        !!<Note that grad v is obtained in element ele from ele2grad_mat

        !!<On the (Dirichlet) boundary we are assembling

        !!< -Int_e (v n. kappa grad u + u n. kappa grad v)

        !!< In practise we'll assemble it everywhere and only
        !!< add it on if we have a Dirichlet boundary

        !!< face_primal_fluxes_mat(1,:,:) maps from ele degrees of freedom
        !!<                                 to internal face dof
        !!< face_primal_fluxes_mat(2,:,:) maps from ele degrees of freedom
        !!<                                 to external face dof
        !!<                                 or face boundary conditions

        !!< For the extra CDG term, we assemble

        !!< -Int_e (C_{12}.[u][kappa grad v] + C_{12}.[v][kappa grad u]

        !!<=-Int C_{12}.(u^+n^+ + u^-n^-)((kappa^+ grad v^+).n^+ +(kappa^- grad v^-).n^-)
        !!<=-Int C_{12}.(v^+n^+ + v^-n^-)((kappa^+ grad u^+).n^+ +(kappa^- grad u^-).n^-)
        !!< Where + is the ele side, and - is the ele_2 side, and e is the
        !! edge

        !!< C_{12} = either (1/2)n^+ or (1/2)n^-
        !!< Take (1/2)n^+ if switch_g . n^+>

        !!<Computing grad u (and v) requires a element transform to physical
        !!<so we only assemble the + parts here, and the minus parts of the grad
        !!<will be assembled when we visit that element

        !!<So we assemble

        !!< - or + Int_e 1/2 (u^+ - u^-) (kappa^+ grad v^+).n^+
        !!< - or + Int_e 1/2 (v^+ - v^-) (kappa^+ grad u^+).n^+

        !! Compare with the primal flux term

        !!<  -Int_e 1/2 (u^+ - u^-)n^+.kappa^+ grad v^+
        !!< -Int_e 1/2 (v^+ - v^-)n^+.kappa^+ grad u^+

        !!< where we take the minus if switch_g.n^+>0 and plus otherwise

        !!< Note that this means that it cancels the primal term if
        !!<switch_g.n^+<0 and doubles it otherwise

        integer :: primal_d1, primal_d2
        real :: primal_flux_factor

#if defined (SCHEME_CDG)
        !        if(viscosity_scheme==CDG) then
            primal_flux_factor = 0.0
            CDG_switch_in = (sum(switch_g(1:NDIM)*sum(face_normal,2)/size(face_normal,2))>0)
            if(CDG_switch_in) primal_flux_factor = 1.0
#else
        !        else
        primal_flux_factor = 0.5
        CDG_switch_in = .true.
        !        end if
#endif

        do primal_d1 = 1, NDIM
            do primal_d2 = 1, NDIM

               ! These get used throughout this section.
               do concurrent (iloc=1:FLOC,jloc=1:FLOC)
                  face_detwei_work = face_detwei * &
                       face_normal(primal_d1,:) * &
                       face_kappa_gi(primal_d1,primal_d2,:)
                  face_shape_shape_work(iloc,jloc)= &
                       dot_product(face_u_shape%n(iloc,:) * &
                       face_u_shape%n(jloc,:), &
                       face_detwei_work)
               end do

                !  -Int_e 1/2 (u^+ - u^-)n^+.kappa^+ grad v^+
                if(.not.face_boundary) then
                    ! Internal face.
                    if(CDG_switch_in) then

!                        face_primal_fluxes_mat(1,:,:) =&
!                            face_primal_fluxes_mat(1,:,:)&
!                            -primal_flux_factor*matmul( &
!                            shape_shape(face_u_shape,face_u_shape, &
!                            & face_detwei * face_normal(primal_d1,:) * face_kappa_gi(primal_d1,primal_d2,:)), &
!                            ele2grad_mat(primal_d2,U_face_l,:))

                       face_primal_fluxes_mat(1,:,:) = &
                            face_primal_fluxes_mat(1,:,:) &
                            -primal_flux_factor*matmul( &
                            face_shape_shape_work, ele2grad_mat(primal_d2,U_face_l,:))

                        ! External face.
!                        face_primal_fluxes_mat(2,:,:) =&
!                            face_primal_fluxes_mat(2,:,:)&
!                            +primal_flux_factor*matmul( &
!                            shape_shape(face_u_shape,face_u_shape, &
!                            & face_detwei * face_normal(primal_d1,:) * face_kappa_gi(primal_d1,primal_d2,:)), &
!                            ele2grad_mat(primal_d2,U_face_l,:))
                        face_primal_fluxes_mat(2,:,:) =&
                            face_primal_fluxes_mat(2,:,:)&
                            +primal_flux_factor*matmul( &
                            face_shape_shape_work, ele2grad_mat(primal_d2,U_face_l,:))
                    end if
                else
                    !If a Dirichlet boundary, we add these terms, otherwise not.

                    !we do the entire integral on the inside face
                    face_primal_fluxes_mat(1,:,:) =&
                        face_primal_fluxes_mat(1,:,:)&
                        -matmul( &
                        face_shape_shape_work, &
                        ele2grad_mat(primal_d2,U_face_l,:))

                    !There is also a corresponding boundary condition integral
                    !on the RHS
!                    face_primal_fluxes_mat(2,:,:) =&
!                        face_primal_fluxes_mat(2,:,:)&
!                        +matmul( &
!                        shape_shape(face_u_shape,face_u_shape, &
!                        & face_detwei * face_normal(primal_d1,:) * face_kappa_gi(primal_d1,primal_d2,:)), &
!                        ele2grad_mat(primal_d2,U_face_l,:))

                    face_primal_fluxes_mat(2,:,:) =&
                        face_primal_fluxes_mat(2,:,:)&
                        +matmul( &
                        face_shape_shape_work, &
                        ele2grad_mat(primal_d2,U_face_l,:))
                end if
            end do
        end do

    end subroutine primal_fluxes



    ! Bassi-Rebay tensor form - ie. homogenous viscosity, incompressible
    subroutine local_assembly_bassi_rebay

        integer :: d3
        real, dimension(NDIM, NDIM, NLOC) :: tensor_visc

        do dim1=1, NDIM
            do dim2=1,NDIM
                do d3 = 1, NDIM

                    ! Div U * G^U * Viscosity * G * Grad U
                    ! Where G^U*G = inverse(Q_mass)
                    Viscosity_mat(d3,d3,:,:)=Viscosity_mat(d3,d3,:,:)&
                        +matmul(matmul(transpose(grad_U_mat_q(dim1,:,:))&
                        &         ,mat_diag_mat(Q_inv, tensor_visc(dim1,dim2,:)))&
                        &     ,grad_U_mat_q(dim2,:,:))

                end do
            end do
        end do

    end subroutine local_assembly_bassi_rebay


    ! Bassi-Rebay tensor form - the one you want for incompressible LES.
    subroutine local_assembly_bassi_rebay_stress_form

        ! Instead of:
        !   M_v = G^T_m (\nu Q^{-1})_mn G_n
        ! We construct:
        !   M_v_rs = G^T_m A_rmsn Q^{-1} G_n
        ! where A is a dim x dim x dim x dim linear operator:
        !   A_rmsn = \partial ( \nu ( u_{r,m} + u_{m,r} ) ) / \partial u_{s,n}
        !   where a_{b,c} = \partial a_b / \partial x_c
        ! off diagonal terms define the coupling between the velocity components

        real, dimension(NLOC, NLOC) :: Q_visc
        real, dimension(NLOC) :: isotropic_visc
        integer :: adim1, adim2, adim3, adim4

        isotropic_visc = Viscosity_ele(1,1,:)

        ! isotropic viscosity (just take the first component as scalar value)
        ! Q_visc = mat_diag_mat(Q_inv, Viscosity_ele(1,1,:))
        ! ** BUG FIX: must use isotropic_visc here, so effect of LES not ignored
        Q_visc = mat_diag_mat( Q_inv, isotropic_visc )

        do adim1=1,NDIM
            do adim2=1,NDIM
                do adim3=1,NDIM
                    do adim4=1,NDIM
                        if (adim1==adim2 .and. adim2==adim3 .and. adim3==adim4) then
                            Viscosity_mat(adim1,adim3,:,:) = Viscosity_mat(adim1,adim3,:,:) &
                                + 2.0 * matmul(matmul(transpose(grad_U_mat_q(adim2,:,:)),Q_visc),grad_U_mat_q(adim4,:,:))
                        else if  ((adim1==adim3 .and. adim2==adim4) .or. (adim2==adim3 .and. adim1==adim4)) then
                            Viscosity_mat(adim1,adim3,:,:) = Viscosity_mat(adim1,adim3,:,:) &
                                + matmul(matmul(transpose(grad_U_mat_q(adim2,:,:)),Q_visc),grad_U_mat_q(adim4,:,:))
                        end if
                    end do
                end do
            end do
        end do

    end subroutine local_assembly_bassi_rebay_stress_form


    subroutine interior_penalty

        !! Ripped from Advection_Diffusion_DG.F90 == cjc

        !! We assemble

        !! Int_e [u][v]

        !! = Int_e C(u^+n^+ + u^-n^-).(v^+n^+ + v^-n^-)

        !! Where + is the ele side, and - is the ele_2 side, and e is the edge
        !! and C is the penalty parameter

        !! We are only storing trial functions from this element, so
        !! will assemble the u^+ parts only, the u^- parts will be done
        !! from the other side

        !!So we assemble

        !!  Int_e C u^+ (v^+ - v^-)

        !!On the (Dirichlet) boundary we are assembling

        !! Int_e C uv

        !! In practise we'll assemble it everywhere and only
        !! add it on if we have a Dirichlet boundary

        !! face_penalty_fluxes_mat(1,:,:) maps from internal face dof
        !!                                 to internal face dof
        !! face_penalty_fluxes_mat(2,:,:) maps from internal face dof
        !!                                 to external face dof
        !!                                 or face boundary conditions

        ! Penalty parameter is C_0/h where h is the distance between the
        ! cell centre and the neighbours cell centre

        real :: C_h
        integer :: ip_d1, ip_d2

        real, dimension(FNGI) :: kappa_n

        ! FIXME: JRM HACK
        !face_kappa_gi = face_kappa_gi + ele_val_at_quad(tensor_eddy_visc, ele)

        kappa_n = 0.0
        do ip_d1 = 1, NDIM
            !            do ip_d2 = 1, NDIM
            !                kappa_n = kappa_n + &
            !                    face_normal(ip_d1,:)*face_kappa_gi(ip_d1,ip_d2,:)*face_normal(ip_d2,:)
            !            end do
            kappa_n = kappa_n + &
                face_normal(ip_d1,:)*face_kappa_gi(ip_d1,1,:)*face_normal(1,:)
            kappa_n = kappa_n + &
                face_normal(ip_d1,:)*face_kappa_gi(ip_d1,2,:)*face_normal(2,:)
            kappa_n = kappa_n + &
                face_normal(ip_d1,:)*face_kappa_gi(ip_d1,3,:)*face_normal(3,:)

        end do

        ! FIXME: JRM HACK
        ! face_kappa_gi = face_kappa_gi - ele_val_at_quad(tensor_eddy_visc, ele)

        if(EDGE_LENGTH_OPTION==USE_FACE_INTEGRALS) then
            h0 = sum(face_detwei)
            h0 = sqrt(h0)
        end if

        if(cdg_penalty) then
            C_h = Interior_Penalty_Parameter
        else
            C_h = Interior_Penalty_Parameter*(h0**edge_length_power)
        end if

        !If a dirichlet boundary then we add these terms, otherwise not

        face_penalty_fluxes_mat(1,:,:) =&
            face_penalty_fluxes_mat(1,:,:)+&
            !     C_h*shape_shape(face_u_shape,face_u_shape,face_detwei)
            C_h*shape_shape(face_u_shape,face_u_shape,face_detwei*kappa_n)

        face_penalty_fluxes_mat(2,:,:) =&
            face_penalty_fluxes_mat(2,:,:)-&
                 !C_h*shape_shape(face_u_shape,face_u_shape,face_detwei)
            C_h*shape_shape(face_u_shape,face_u_shape,face_detwei*kappa_n)

    end subroutine interior_penalty

    subroutine local_assembly_ip_face
        implicit none

        integer :: laip_d
        ! Don't need these I think
        !        integer :: nfele, nele
        !integer, dimension(FLOC) :: laip_U_face_loc

        !laip_U_face_loc=face_local_nodes(U, face)



        if (face_boundary) then
            do laip_d=1, NDIM
                if(face_dirichlet(laip_d)) then
                    !!These terms are not included on Neumann integrals

                    !! Internal Degrees of Freedom

                    !penalty flux

                    Viscosity_mat(laip_d,laip_d,u_face_l,u_face_l) = &
                        Viscosity_mat(laip_d,laip_d,u_face_l,u_face_l) + &
                        face_penalty_fluxes_mat(1,:,:)

                    !! External Degrees of Freedom

                    !!penalty fluxes

                    Viscosity_mat(laip_d,laip_d,u_face_l,face_start:face_finish) = &
                        Viscosity_mat(laip_d,laip_d,u_face_l,face_start:face_finish) + &
                        face_penalty_fluxes_mat(2,:,:)

                end if
            end do
        else
            do laip_d=1, NDIM
                !! Internal Degrees of Freedom

                !penalty flux

                Viscosity_mat(laip_d,laip_d,u_face_l,u_face_l) = &
                    Viscosity_mat(laip_d,laip_d,u_face_l,u_face_l) + &
                    face_penalty_fluxes_mat(1,:,:)

                !! External Degrees of Freedom

                !!penalty fluxes

                Viscosity_mat(laip_d,laip_d,u_face_l,face_start:face_finish) = &
                    Viscosity_mat(laip_d,laip_d,u_face_l,face_start:face_finish) + &
                    face_penalty_fluxes_mat(2,:,:)

            end do
        end if

    end subroutine local_assembly_ip_face




    subroutine local_assembly_primal_face
        implicit none

        integer :: pface_j, pface_dim
        ! integer, dimension(FLOC) :: pface_U_face_loc

        ! pface_U_face_loc=face_local_nodes(U, face)


        if (face_boundary) then
            do pface_dim=1, NDIM
                if(face_dirichlet(pface_dim)) then
                    !!These terms are not included on Neumann integrals

                    !! Internal Degrees of Freedom

                    !primal fluxes

                    Viscosity_mat(pface_dim,pface_dim,u_face_l,1:NLOC) = &
                        Viscosity_mat(pface_dim,pface_dim,u_face_l,1:NLOC) + &
                        face_primal_fluxes_mat(1,:,:)

                    do pface_j = 1, FLOC
                        Viscosity_mat(pface_dim,pface_dim,1:NLOC,u_face_l(pface_j)) = &
                            Viscosity_mat(pface_dim,pface_dim,1:NLOC,u_face_l(pface_j)) + &
                            face_primal_fluxes_mat(1,pface_j,:)
                    end do

                    !primal fluxes

                    Viscosity_mat(pface_dim,pface_dim,1:NLOC,face_start:face_finish) = &
                        Viscosity_mat(pface_dim,pface_dim,1:NLOC,face_start:face_finish) + &
                        transpose(face_primal_fluxes_mat(2,:,:))

                end if
            end do
        else
            do pface_dim=1, NDIM
                !! Internal Degrees of Freedom

                !primal fluxes

                Viscosity_mat(pface_dim,pface_dim,u_face_l,1:NLOC) = &
                    Viscosity_mat(pface_dim,pface_dim,u_face_l,1:NLOC) + &
                    face_primal_fluxes_mat(1,:,:)

                do pface_j = 1, FLOC
                    Viscosity_mat(pface_dim,pface_dim,1:NLOC,u_face_l(pface_j)) = &
                        Viscosity_mat(pface_dim,pface_dim,1:NLOC,u_face_l(pface_j)) + &
                        face_primal_fluxes_mat(1,pface_j,:)
                end do

                !! External Degrees of Freedom

                !primal fluxes

                Viscosity_mat(pface_dim,pface_dim,face_start:face_finish,1:NLOC) = &
                    Viscosity_mat(pface_dim,pface_dim,face_start:face_finish,1:NLOC) + &
                    face_primal_fluxes_mat(2,:,:)

                Viscosity_mat(pface_dim,pface_dim,1:NLOC,face_start:face_finish) = &
                    Viscosity_mat(pface_dim,pface_dim,1:NLOC,face_start:face_finish) + &
                    transpose(face_primal_fluxes_mat(2,:,:))

            end do
        end if

    end subroutine local_assembly_primal_face

    subroutine local_assembly_cdg_face
        implicit none
        !!< This code assembles the cdg fluxes involving the r_e and l_e lifting
        !!< operators.

        !!< We assemble the operator
        !!< \int (r^e([v]) + l^e(C_{12}.[v]) + r^e_D(v).\kappa.
        !!< (r^e([u]) + l^e(C_{12}.[u]) + r^e_D(u))dV (*)
        !!< This is done by forming the operator R:
        !!< \int v R(u)dV  = \int v (r^e([u]) + l^e(C_{12}.[u]) + r^e_D(u)) dV
        !!< and then constructing
        !!< \int R(v).\kappa.R(u) dV

        !!< The lifting operator r^e is defined by
        !!< \int_E \tau . r^e([u]) dV = - \int_e {\tau}.[u] dS
        !!< = -\frac{1}{2} \int_e {\tau^+ + \tau^-}.(u^+n^+ + u^-n^-) dS
        !!< = -\frac{1}{2} \int_e {\tau^+ + \tau^-}.n^+(u^+ - u^-) dS

        !!< Where + is the ele side, and - is the ele_2 side, and e is the edge

        !!< The lifting operator l^e is defined by
        !!< \int_E \tau . l^e(C_{12}.[u])dV = - \int_e C_{12}.[u][\tau] dS
        !!< = -\int C_{12}.(u^+n^+ + u^-n^-)(\tau^+.n^+ +\tau^-n^-) dS

        !!< C_{12} = either (1/2)n^+ or (1/2)n^-
        !!< Take (1/2)n^+ if switch_g . n^+> 0

        !!becomes
        !!< = \int_e (- or +)(u^+ - u^-)n^+.(\tau^+ - \tau^-) dS
        !!< with minus sign if switch_g  n^+ > 0

        !!< So adding r^e and l^e gives

        !!< = -\frac{1}{2} \int_e {\tau^+ + \tau^-}.n^+(u^+ - u^-) dS
        !!<     + \int_e (- or +)(u^+ - u^-)n^+.(\tau^+ - \tau^-) dS

        !!< = -\int_e \tau^+.n^+(u^+ - u^-) dS if switch_g n^+ > 0
        !!< = -\int_e \tau^-.n^+(u^+ - u^-) dS otherwise

        !!< so definition of r^e+l^e operator is
        !!< \int_E \tau.R(u) dV = -\int_e \tau^+.n^+(u^+ - u^-) dS if switch > 0
        !!< \int_E \tau.R(u) dV = -\int_e \tau^-.n^+(u^+ - u^-) dS if switch < 0

        !!< we are doing DG so the basis functions which are non-zero in E are
        !!< zero outside E, so \tau^- vanishes in this formula, so we get
        !!< \int_E \tau.R(u) dV = -\int_e \tau.n^+(u^+ - u^-) dS if switch > 0
        !!< and R(u) = 0 otherwise.

        !!< finally the boundary lifting operator r^e_D
        !!< \int_E \tau.r^e_D(u) dV =  -\int_e u\tau.n dS

        !!< We assemble the binary form (*) locally with
        !!< B(u,v) = p^TR^T.K.Rq, where p is the vector of coefficients of u in
        !!< element E plus the coefficients of u on the face e on the other side
        !!< K is the matrix obtained from the bilinear form
        !!< \int_E N_i \kappa N_j dV where \kappa is the viscosity tensor and
        !! N_i are the basis functions with support in element E

        !!< The matrix R maps from the coefficients of a scalar field  on both sides of face e
        !!< to the coefficients of a vector field with inside element E
        !!< i.e. size (dim x loc(E),2 x loc(e))
        !!< because of symmetry we just store (dim x loc(E), loc(e)) values
        !!< The matrix K maps from vector fields inside element E to vector
        !!< fields inside element E
        !!< i.e. size (dim x loc(E), dim x loc(E))
        !!< Hence, R^TKR maps from the coefficients of a scalar field on both
        !!< sides of face e to themselves
        !!< i.e. size (2 x loc(E), 2 x
        !!< It can be thus interpreted as a fancy penalty term for
        !!< discontinuities, a useful one because it is scale invariant

        !!< The matrix R can be formed by constructing the bilinear form matrix
        !!< for r^e, l^e and r^e_D, and then dividing by the elemental mass
        !!<  matrix on E

        !!< we place R^TKR into Viscosity_mat which maps from u
        !!< coefficients in element E plus those on the other side of face e
        !!< to themselves, hence it has size (loc(E) + loc(e), loc(E) + loc(e))

        !!< R^TKR is stored in add_mat which has size(2 x loc(e), 2 x loc(e))

        !!< we are using a few other pre-assembled local matrices
        !!< normal_mat is \int_e \tau.(un) dS (has size (dim x loc(e),loc(e))
        !!< normal_kappa_mat is \int_e \tau.\kappa.(un) dS
        !!< has size (dim x loc(e), loc(e))
        !!< inverse_mass_mat is the inverse mass in E


        ! Internal variables
        integer :: cdg_i, cdg_j, cdg_dim1, cdg_dim2, cdg_face1, cdg_face2, cdg_outer_dim
        ! integer, dimension(FLOC) :: cdg_U_face_loc
        real, dimension(NDIM, NLOC, FLOC) :: cdg_R_mat
        real, dimension(2,2, FLOC, FLOC) :: cdg_add_mat

        ! cdg_U_face_loc=face_local_nodes(U, face)

        cdg_R_mat = 0.
        do concurrent (cdg_dim1=1:NDIM, cdg_i=1:NLOC, cdg_j=1:FLOC)
            cdg_R_mat(cdg_dim1,cdg_i,cdg_j) = &
                &sum(inverse_mass_mat(cdg_i,u_face_l)*face_normal_mat(cdg_dim1, :, cdg_j))
        end do

        do cdg_outer_dim=1, NDIM

            cdg_add_mat = 0.0
            if(face_boundary) then
                if (face_dirichlet(cdg_outer_dim)) then
                    !Boundary case
                    ! R(/tau,u) = -\int_e \tau.n u  dS
                    !do cdg_dim1 = 1, mesh_dim(U)
                    !   do cdg_dim2 = 1, mesh_dim(U)
                    !      cdg_add_mat(1,1,:,:) = cdg_add_mat(1,1,:,:) + &
                    !           matmul(transpose(cdg_R_mat(cdg_dim1,:,:)), &
                    !           &matmul(kappa_mat(cdg_dim1,cdg_dim2,:,:),cdg_R_mat(cdg_dim2,:,:)))
                    !      cdg_add_mat(2,2,:,:) = cdg_add_mat(2,2,:,:) + &
                    !           matmul(transpose(cdg_R_mat(cdg_dim1,:,:)), &
                    !           &matmul(kappa_mat(cdg_dim1,cdg_dim2,:,:),cdg_R_mat(cdg_dim2,:,:)))
                    !   end do
                    !end do

                    do cdg_face1 = 1, 2
                        do cdg_face2 = 1, 2
                            do cdg_dim1 = 1, NDIM
                                do cdg_dim2 = 1, NDIM
                                    cdg_add_mat(cdg_face1,cdg_face2,:,:) = cdg_add_mat(cdg_face1,cdg_face2,:,:) + &
                                        &(-1.)**(cdg_face1+cdg_face2)*matmul(transpose(cdg_R_mat(cdg_dim1,:,:)), &
                                        &matmul(kappa_mat(cdg_dim1,cdg_dim2,:,:),cdg_R_mat(cdg_dim2,:,:)))
                                end do
                            end do
                        end do
                    end do

                end if
            else if(CDG_switch_in) then
                ! interior case
                ! R(\tau,u) = -\int_e \tau.n^+(u^+ - u^-) dS
                do cdg_face1 = 1, 2
                    do cdg_face2 = 1, 2
                        do cdg_dim1 = 1, NDIM
                            do cdg_dim2 = 1, NDIM
                                cdg_add_mat(cdg_face1,cdg_face2,:,:) = cdg_add_mat(cdg_face1,cdg_face2,:,:) + &
                                    &(-1.)**(cdg_face1+cdg_face2)*matmul(transpose(cdg_R_mat(cdg_dim1,:,:)), &
                                    &matmul(kappa_mat(cdg_dim1,cdg_dim2,:,:),cdg_R_mat(cdg_dim2,:,:)))
                            end do
                        end do
                    end do
                end do
            end if

            !cdg_face1 = 1, cdg_face2 = 1

            Viscosity_mat(cdg_outer_dim,cdg_outer_dim, u_face_l,u_face_l) = &
                &Viscosity_mat(cdg_outer_dim, cdg_outer_dim, u_face_l,u_face_l) + &
                &cdg_add_mat(1,1,:,:)

            !cdg_face1 = 1, cdg_face2 = 2

            Viscosity_mat(cdg_outer_dim, cdg_outer_dim, u_face_l,face_start:face_finish) = &
                &Viscosity_mat(cdg_outer_dim, cdg_outer_dim, u_face_l,face_start:face_finish) + &
                &cdg_add_mat(1,2,:,:)

            !cdg_face1 = 2, cdg_face2 = 1

            Viscosity_mat(cdg_outer_dim,cdg_outer_dim,face_start:face_finish,u_face_l) = &
                Viscosity_mat(cdg_outer_dim,cdg_outer_dim,face_start:face_finish,u_face_l) + &
                &cdg_add_mat(2,1,:,:)

            !cdg_face1 = 2, cdg_face2 = 2

            Viscosity_mat(cdg_outer_dim,cdg_outer_dim,face_start:face_finish,face_start:face_finish) = &
                &Viscosity_mat(cdg_outer_dim,cdg_outer_dim,face_start:face_finish,face_start:face_finish) + &
                &cdg_add_mat(2,2,:,:)
        end do

    end subroutine local_assembly_cdg_face



end subroutine construct_momentum_elements_dg_ELEMENT_CONFIG
