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

module advection_local_DG
  !!< This module contains the Discontinuous Galerkin form of the advection
  !!< -diffusion equation for scalars.
  use elements
  use sparse_tools
  use fetools
  use dgtools
  use fields
  use fefields
  use state_module
  use shape_functions
  use global_numbering
  use transform_elements
  use vector_tools
  use fldebug
  use vtk_interfaces
  use Coordinates
  use petsc_solve_state_module
  use boundary_conditions
  use boundary_conditions_from_options
  use spud
  use upwind_stabilisation
  use slope_limiters_dg
  use sparsity_patterns
  use sparse_matrices_fields
  use sparsity_patterns_meshes
  use manifold_projections
  use diagnostic_fields, only: calculate_diagnostic_variable
  use global_parameters, only : FIELD_NAME_LEN

  implicit none

  ! Buffer for output messages.
  character(len=255), private :: message

  private
  public solve_advection_dg_subcycle, solve_vector_advection_dg_subcycle

  ! Local private control parameters. These are module-global parameters
  ! because it would be expensive and/or inconvenient to re-evaluate them
  ! on a per-element or per-face basis
  real :: dt, theta

  ! Whether to include various terms
  logical :: include_advection
  ! Discretisation to use for advective flux.
  integer :: flux_scheme
  integer, parameter :: UPWIND_FLUX=1
  integer, parameter :: LAX_FRIEDRICHS_FLUX=2
  
  ! Boundary condition types:
  ! (the numbers should match up with the order in the 
  !  get_entire_boundary_condition call)
  integer :: BCTYPE_WEAKDIRICHLET=1, BCTYPE_DIRICHLET=2

  logical :: include_mass

  ! Stabilisation schemes.
  integer :: stabilisation_scheme
  integer, parameter :: NONE=0
  integer, parameter :: UPWIND=1

contains

  subroutine solve_advection_dg_subcycle(field_name, state, velocity_name)
    !!< Construct and solve the advection equation for the given
    !!< field using discontinuous elements.
    
    !! Name of the field to be solved for.
    character(len=*), intent(in) :: field_name
    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state
    !! Optional velocity name
    character(len = *), intent(in) :: velocity_name

    !! Tracer to be solved for.
    type(scalar_field), pointer :: T, T_old, s_field

    !! Coordinate and velocity fields
    type(vector_field), pointer :: X, U_nl

    !! Velocity field in Cartesian coordinates for slope limiter
    type(vector_field) :: U_nl_cartesian

    !! Change in T over one timestep.
    type(scalar_field) :: delta_T

    !! Sparsity of advection matrix.    
    type(csr_sparsity), pointer :: sparsity
    
    !! System matrix.
    type(csr_matrix) :: matrix, mass, inv_mass

    !! Sparsity of mass matrix.
    type(csr_sparsity) :: mass_sparsity

    !! Right hand side vector.
    type(scalar_field) :: rhs

    !! Whether to invoke the slope limiter
    logical :: limit_slope
    !! Which limiter to use
    integer :: limiter

    !! Number of advection subcycles.
    integer :: subcycles
    real :: max_courant_number

    character(len=FIELD_NAME_LEN) :: limiter_name
    integer :: i

    T=>extract_scalar_field(state, field_name)
    T_old=>extract_scalar_field(state, "Old"//field_name)
    X=>extract_vector_field(state, "Coordinate")

    ! Reset T to value at the beginning of the timestep.
    call set(T, T_old)

    sparsity => get_csr_sparsity_firstorder(state, T%mesh, T%mesh)

    call allocate(matrix, sparsity) ! Add data space to the sparsity
    ! pattern.

    mass_sparsity=make_sparsity_dg_mass(T%mesh)
    call allocate(mass, mass_sparsity)
    call allocate(inv_mass, mass_sparsity)

    ! Ensure delta_T inherits options from T.
    call allocate(delta_T, T%mesh, "delta_T")
    delta_T%option_path = T%option_path
    call allocate(rhs, T%mesh, trim(field_name)//" RHS")
   
    call construct_advection_dg(matrix, rhs, field_name, state, &
         mass, velocity_name=velocity_name)

    call get_dg_inverse_mass_matrix(inv_mass, mass)
    
    ! Note that since theta and dt are module global, these lines have to
    ! come after construct_advection_diffusion_dg.
    call get_option("/timestepping/timestep", dt)
    
    if(have_option(trim(T%option_path)//&
         &"/prognostic/temporal_discretisation/discontinuous_galerkin/&
         &/number_advection_subcycles")) then
       call get_option(trim(T%option_path)//&
            &"/prognostic/temporal_discretisation/discontinuous_galerkin/&
            &/number_advection_subcycles", subcycles)
    else
       call get_option(trim(T%option_path)//&
            &"/prognostic/temporal_discretisation/discontinuous_galerkin/&
            &/maximum_courant_number_per_subcycle", Max_Courant_number)
       
       s_field => extract_scalar_field(state, "DG_CourantNumber")
       call calculate_diagnostic_variable(state, "DG_CourantNumber_Local", &
            & s_field)
       
       subcycles = ceiling( maxval(s_field%val)/Max_Courant_number)
       call allmax(subcycles)
       ewrite(2,*) 'Number of subcycles for tracer eqn: ', subcycles
    end if

    limit_slope=.false.
    if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/slope_limiter")) then
       limit_slope=.true.
       
       ! Note unsafe for mixed element meshes
       if (element_degree(T,1)==0) then
          FLExit("Slope limiters make no sense for degree 0 fields")
       end if

       call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter/name",limiter_name)

       select case(trim(limiter_name))
       case("Cockburn_Shu")
          limiter=LIMITER_COCKBURN
       case("Hermite_Weno")
          limiter=LIMITER_HERMITE_WENO
       case("minimal")
          limiter=LIMITER_MINIMAL
       case("FPN")
          limiter=LIMITER_FPN
       case("Vertex_Based")
          limiter=LIMITER_VB
       case default
          FLAbort('No such limiter')
       end select
       
    end if

    U_nl=>extract_vector_field(state, velocity_name)
    ! Allocate and zero cartesian velocity field
    call allocate(U_nl_cartesian, U_nl%dim, U_nl%mesh, "U_nl_cartesian")
    call zero(U_nl_cartesian)

    do i=1, subcycles

       ! dT = Advection * T
       call mult(delta_T, matrix, T)
       ! dT = dT + RHS
       call addto(delta_T, RHS, -1.0)
       ! dT = M^(-1) dT
       call dg_apply_mass(inv_mass, delta_T)
       
       ! T = T + dt/s * dT
       call addto(T, delta_T, scale=-dt/subcycles)
       call halo_update(T)
       ewrite_minmax(delta_T)
       if (limit_slope) then
          ! Filter wiggles from T
          call limit_slope_dg(T, U_nl, X, state, limiter)
       end if

    end do

    call deallocate(delta_T)
    call deallocate(matrix)
    call deallocate(mass)
    call deallocate(inv_mass)
    call deallocate(mass_sparsity)
    call deallocate(rhs)

  end subroutine solve_advection_dg_subcycle

  subroutine construct_advection_dg(big_m, rhs, field_name,&
       & state, mass, velocity_name) 
    !!< Construct the advection equation for discontinuous elements in
    !!< acceleration form.
    !!< 
    !!< If mass is provided then the mass matrix is not added into big_m or
    !!< rhs. It is instead returned as mass. This may be useful for testing
    !!< or for solving equations otherwise than in acceleration form.

    !! Main advection matrix.    
    type(csr_matrix), intent(inout) :: big_m
    !! Right hand side vector.
    type(scalar_field), intent(inout) :: rhs
    
    !! Name of the field to be advected.
    character(len=*), intent(in) :: field_name
    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state
    !! Optional separate mass matrix.
    type(csr_matrix), intent(inout), optional :: mass
    !! Optional velocity name
    character(len = *), intent(in), optional :: velocity_name

    !! Position, and velocity fields.
    type(vector_field) :: X, U, U_nl
    !! Tracer to be solved for.
    type(scalar_field) :: T

    !! Local velocity name
    character(len = FIELD_NAME_LEN) :: lvelocity_name

    !! Element index
    integer :: ele

    !! Status variable for field extraction.
    integer :: stat

    !! Field over the entire surface mesh containing bc values:
    type(scalar_field) :: bc_value
    !! Integer array of all surface elements indicating bc type
    !! (see below call to get_entire_boundary_condition):
    integer, dimension(:), allocatable :: bc_type

    ewrite(1,*) "Writing advection equation for "&
         &//trim(field_name)

    ! These names are based on the CGNS SIDS.
    T=extract_scalar_field(state, field_name)
    X=extract_vector_field(state, "Coordinate")

    if(present(velocity_name)) then
      lvelocity_name = velocity_name
    else
      lvelocity_name = "NonlinearVelocity"
    end if

    if (.not.have_option(trim(T%option_path)//"/prognostic&
         &/spatial_discretisation/discontinuous_galerkin&
         &/advection_scheme/none")) then
       U_nl=extract_vector_field(state, lvelocity_name)
       call incref(U_nl)
       include_advection=.true.
    else
       ! Forcing a zero NonlinearVelocity will disable advection.
       U=extract_vector_field(state, "Velocity", stat)
       if (stat/=0) then 
          FLExit("Oh dear, no velocity field. A velocity field is required for advection!")
       end if
       call allocate(U_nl, U%dim, U%mesh, "LocalNonlinearVelocity")
       call zero(U_nl)
       include_advection=.false.
    end if

    flux_scheme=UPWIND_FLUX
    if (have_option(trim(T%option_path)//"/prognostic&
         &/spatial_discretisation/discontinuous_galerkin&
         &/advection_scheme/lax_friedrichs")) then
       flux_scheme=LAX_FRIEDRICHS_FLUX
    end if

    include_mass = .not. have_option(trim(T%option_path)//&
           "/prognostic/spatial_discretisation/discontinuous_galerkin/mass_terms/exclude_mass_terms")
           
    ! Switch on upwind stabilisation if requested.
    if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/upwind_stabilisation")) then
       stabilisation_scheme=UPWIND
    else
       stabilisation_scheme=NONE
    end if

    assert(has_faces(X%mesh))
    assert(has_faces(T%mesh))
    
    ! Enquire about boundary conditions we're interested in
    ! Returns an integer array bc_type over the surface elements
    ! that indicates the bc type (in the order we specified, i.e.
    ! BCTYPE_WEAKDIRICHLET=1)
    allocate( bc_type(1:surface_element_count(T)) )
    call get_entire_boundary_condition(T, &
       & (/"weakdirichlet"/), &
       & bc_value, bc_type)

    call zero(big_m)
    call zero(RHS)
    if (present(mass)) call zero(mass)

    element_loop: do ele=1,element_count(T)
       
       call construct_adv_element_dg(ele, big_m, rhs,&
            & X, T, U_nl, &
            bc_value, bc_type, mass)
       
    end do element_loop
    
    ! Drop any extra field references.

    call deallocate(U_nl)
    call deallocate(bc_value)

  end subroutine construct_advection_dg

  subroutine construct_adv_element_dg(ele, big_m, rhs,&
       & X, T, U_nl, &
       & bc_value, bc_type, &
       & mass)
    !!< Construct the advection_diffusion equation for discontinuous elements in
    !!< acceleration form.
    implicit none
    !! Index of current element
    integer :: ele
    !! Main advection matrix.
    type(csr_matrix), intent(inout) :: big_m
    !! Right hand side vector.
    type(scalar_field), intent(inout) :: rhs
    !! Field over the entire surface mesh containing bc values:
    type(scalar_field), intent(in):: bc_value
    !! Integer array of all surface elements indicating bc type
    !! (see above call to get_entire_boundary_condition):
    integer, dimension(:), intent(in):: bc_type
    !! Optional separate mass matrix.
    type(csr_matrix), intent(inout), optional :: mass
    
    !! Position and velocity.
    type(vector_field), intent(in) :: X, U_nl

    type(scalar_field), intent(in) :: T

    ! Bilinear forms.
    real, dimension(ele_loc(T,ele), ele_loc(T,ele)) :: &
         mass_mat
    real, dimension(ele_loc(T,ele), ele_loc(T,ele)) :: &
         Advection_mat, Ad_mat1, Ad_mat2

    ! Local assembly matrices.
    real, dimension(ele_loc(T,ele), ele_loc(T,ele)) :: l_T_mat
    real, dimension(ele_loc(T,ele)) :: l_T_rhs

    ! Local variables.
    
    ! Neighbour element, face and neighbour face.
    integer :: ele_2, face, face_2
    ! Loops over faces.
    integer :: ni
    
    ! Transform from local to physical coordinates.
    real, dimension(ele_ngi(X,ele)) :: detwei
    real, dimension(mesh_dim(U_nl), X%dim, ele_ngi(X,ele)) :: J_mat 

    ! Different velocities at quad points.
    real, dimension(U_nl%dim, ele_ngi(U_nl, ele)) :: U_quad
    real, dimension(U_nl%dim, ele_ngi(U_nl, ele)) :: U_nl_q
    real, dimension(ele_ngi(U_nl, ele)) :: U_nl_div_q

    ! Node and shape pointers.
    integer, dimension(:), pointer :: T_ele
    type(element_type), pointer :: T_shape, U_shape
    ! Neighbours of this element.
    integer, dimension(:), pointer :: neigh, x_neigh
    ! Whether the tracer field is continuous.
    logical :: dg

    integer :: gi,i,j

    logical :: boundary_element

    dg=continuity(T)<0

    !----------------------------------------------------------------------
    ! Establish local node lists
    !----------------------------------------------------------------------
    
    T_ele=>ele_nodes(T,ele)  ! Tracer node numbers

    !----------------------------------------------------------------------
    ! Establish local shape functions
    !----------------------------------------------------------------------

    T_shape=>ele_shape(T, ele)
    U_shape=>ele_shape(U_nl, ele)

    !==========================
    ! Coordinates
    !==========================

    ! Get J_mat
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), detwei=detwei, J=J_mat)
    
    !----------------------------------------------------------------------
    ! Construct bilinear forms.
    !----------------------------------------------------------------------

    ! Element density matrix.
    !  /
    !  | T T  dV
    !  / 
    mass_mat = shape_shape(T_shape, T_shape, detwei)

    if (include_advection) then

      ! Advecting velocity at quadrature points.
      U_quad=ele_val_at_quad(U_nl,ele)
      U_nl_q=U_quad

      U_nl_div_q=ele_div_at_quad(U_nl, ele, U_shape%dn)

      ! Element advection matrix
      !    /                           /
      !  - | (grad T dot U_nl) T dV -  | T ( div U_nl ) T dV
      !    /                           /

      Ad_mat1 = -dshape_dot_vector_shape(T_shape%dn, U_nl_q, T_shape, T_shape%quadrature%weight)
      Ad_mat2 = -shape_shape(T_shape, T_shape, U_nl_div_q * T_shape%quadrature%weight)
      Advection_mat = -dshape_dot_vector_shape(T_shape%dn, U_nl_q, T_shape, T_shape%quadrature%weight)  &
           -shape_shape(T_shape, T_shape, U_nl_div_q * T_shape%quadrature%weight)
   else
      Advection_mat=0.0
   end if

   !----------------------------------------------------------------------
   ! Perform global assembly.
   !----------------------------------------------------------------------

   l_T_rhs=0.0

   ! Assemble matrix.
    
   ! Advection.
   l_T_mat= Advection_mat

   if (present(mass)) then
      ! Return mass separately.
      call addto(mass, t_ele, t_ele, mass_mat)
   else
      if(include_mass) then
         ! Put mass in the matrix.
         l_T_mat=l_T_mat+mass_mat
      end if
   end if

   call addto(big_m, t_ele, t_ele, l_T_mat)

   !-------------------------------------------------------------------
   ! Interface integrals
   !-------------------------------------------------------------------
    
   neigh=>ele_neigh(T, ele)
   ! x_neigh/=t_neigh only on periodic boundaries.
   x_neigh=>ele_neigh(X, ele)

   ! Flag for whether this is a boundary element.
   boundary_element=.false.

   neighbourloop: do ni=1,size(neigh)

      !----------------------------------------------------------------------
      ! Find the relevant faces.
      !----------------------------------------------------------------------
       
      ! These finding routines are outside the inner loop so as to allow
      ! for local stack variables of the right size in
      ! construct_add_diff_interface_dg.

      ele_2=neigh(ni)
       
      ! Note that although face is calculated on field U, it is in fact
      ! applicable to any field which shares the same mesh topology.
      face=ele_face(T, ele, ele_2)
    
      if (ele_2>0) then
         ! Internal faces.
         face_2=ele_face(T, ele_2, ele)
      else
         ! External face.
         face_2=face
         boundary_element=.true.
      end if

      call construct_adv_interface_dg(ele, face, face_2,&
           & big_m, rhs, X, T, U_nl,&
           & bc_value, bc_type)

   end do neighbourloop
    
 end subroutine construct_adv_element_dg
  
  subroutine construct_adv_interface_dg(ele, face, face_2, &
       big_m, rhs, &
       & X, T, U_nl,&
       & bc_value, bc_type)

    !!< Construct the DG element boundary integrals on the ni-th face of
    !!< element ele.
    implicit none

    integer, intent(in) :: ele, face, face_2
    type(csr_matrix), intent(inout) :: big_m
    type(scalar_field), intent(inout) :: rhs
    ! We pass these additional fields to save on state lookups.
    type(vector_field), intent(in) :: X, U_nl
    type(scalar_field), intent(in) :: T
   !! Field over the entire surface mesh containing bc values:
    type(scalar_field), intent(in):: bc_value
    !! Integer array of all surface elements indicating bc type
    !! (see above call to get_entire_boundary_condition):
    integer, dimension(:), intent(in):: bc_type

    ! Face objects and numberings.
    type(element_type), pointer :: T_shape, T_shape_2
    integer, dimension(face_loc(T,face)) :: T_face, T_face_l
    integer, dimension(face_loc(T,face_2)) :: T_face_2
    integer, dimension(face_loc(U_nl,face)) :: U_face
    integer, dimension(face_loc(U_nl,face_2)) :: U_face_2

    ! Note that both sides of the face can be assumed to have the same
    ! number of quadrature points.
    real, dimension(U_nl%dim, face_ngi(U_nl, face)) :: n1, n2, normal, U_nl_q,&
         & u_f_q, u_f2_q, div_u_f_q
    logical, dimension(face_ngi(U_nl, face)) :: inflow
    real, dimension(face_ngi(U_nl, face)) :: U_nl_q_dotn1, U_nl_q_dotn2, U_nl_q_dotn, income
    ! Variable transform times quadrature weights.
    real, dimension(face_ngi(T,face)) :: detwei
    real, dimension(U_nl%dim, X%dim, face_ngi(T,face)) :: J
    real, dimension(face_ngi(T,face)) :: inner_advection_integral, outer_advection_integral

    ! Bilinear forms
    real, dimension(face_loc(T,face),face_loc(T,face)) :: nnAdvection_out
    real, dimension(face_loc(T,face),face_loc(T,face_2)) :: nnAdvection_in

    logical :: boundary, dirichlet

    ! Lax-Friedrichs flux parameter
    real :: C

    integer :: i

    T_face=face_global_nodes(T, face)
    T_face_l=face_local_nodes(T, face)
    T_shape=>face_shape(T, face)

    T_face_2=face_global_nodes(T, face_2)
    T_shape_2=>face_shape(T, face_2)
    
    ! Boundary nodes have both faces the same.
    boundary=(face==face_2)
    dirichlet=.false.
    if (boundary) then
       if (bc_type(face)==BCTYPE_WEAKDIRICHLET) then
         dirichlet=.true.
       end if
    end if

    !Unambiguously calculate the normal using the face with the higher
    !face number. This is so that the normal is identical on both sides.
    ! Jemma: need to be more careful here - actually have to calculate normal for face2

    n1=get_normal(local_face_number(T%mesh,face))
    n2=get_normal(local_face_number(T%mesh,face_2))
    
    !----------------------------------------------------------------------
    ! Construct element-wise quantities.
    !----------------------------------------------------------------------

    if (include_advection.and.(flux_scheme==UPWIND_FLUX)) then
       
       ! Advecting velocity at quadrature points.
       U_f_q = face_val_at_quad(U_nl, face)
       U_f2_q = face_val_at_quad(U_nl, face_2)

       if(local_face_number(T%mesh, face)==3) then
          U_f_q=sqrt(2.)*U_f_q
       end if
       if(local_face_number(T%mesh, face_2)==3) then
          U_f2_q=sqrt(2.)*U_f2_q
       end if

       u_nl_q_dotn1 = sum(U_f_q*n1,1)
       u_nl_q_dotn2 = -sum(U_f2_q*n2,1)
       U_nl_q_dotn=0.5*(u_nl_q_dotn1+u_nl_q_dotn2)

       ! Inflow is true if the flow at this gauss point is directed
       ! into this element.
       inflow = u_nl_q_dotn<0.0
       income = merge(1.0,0.0,inflow)

       !----------------------------------------------------------------------
       ! Construct bilinear forms.
       !----------------------------------------------------------------------
       
       ! Calculate outflow boundary integral.
       ! can anyone think of a way of optimising this more to avoid
       ! superfluous operations (i.e. multiplying things by 0 or 1)?

       ! first the integral around the inside of the element
       ! (this is the flux *out* of the element)
       inner_advection_integral = (1.-income)*u_nl_q_dotn
       nnAdvection_out=shape_shape(T_shape, T_shape,  &
            &                     inner_advection_integral*T_shape%quadrature%weight)
!       print*, 'js:nnAdvection_out:'
!       print*, nnAdvection_out
       ! now the integral around the outside of the element
       ! (this is the flux *in* to the element)
       outer_advection_integral = income*u_nl_q_dotn
       nnAdvection_in=shape_shape(T_shape, T_shape_2, &
            &                       outer_advection_integral*T_shape%quadrature%weight)
!       print*, 'js:nnAdvection_in:'
!       print*, nnAdvection_in
       
    else if (include_advection.and.(flux_scheme==LAX_FRIEDRICHS_FLUX)) then

       FLExit("Haven't worked out Lax-Friedrichs yet")

    end if

    !----------------------------------------------------------------------
    ! Perform global assembly.
    !----------------------------------------------------------------------

    ! Insert advection in matrix.
    if (include_advection) then
    
       ! Outflow boundary integral.
       call addto(big_M, T_face, T_face,&
            nnAdvection_out)
       
       if (.not.dirichlet) then
          ! Inflow boundary integral.
          call addto(big_M, T_face, T_face_2,&
               nnAdvection_in)
       end if

       ! Insert advection in RHS.
       
       if (dirichlet) then

          ! Inflow and outflow of Dirichlet value.
          call addto(RHS, T_face, &
               -matmul(nnAdvection_in,&
               ele_val(bc_value, face)))
          
       end if
       
    end if

    contains

      function get_normal(face) result(norm)

        integer, intent(in) :: face
        real, dimension(U_nl%dim, face_ngi(U_nl,face)) :: norm

        integer :: i

        if(U_nl%dim==1) then
           if(face==1) then
              forall(i=1:face_ngi(U_nl,face)) norm(1,i)=1.
           else if(face==2) then
              forall(i=1:face_ngi(U_nl,face)) norm(1,i)=-1.
           end if
        else if(U_nl%dim==2) then
           if(face==1) then
              forall(i=1:face_ngi(U_nl,face)) norm(1:2,i)=(/-1.,0./)
           else if(face==2) then
              forall(i=1:face_ngi(U_nl,face)) norm(1:2,i)=(/0.,-1./)
           else if(face==3) then
              forall(i=1:face_ngi(U_nl,face)) norm(1:2,i)=(/1/sqrt(2.),1/sqrt(2.)/)
           else
              FLAbort('Oh dear oh dear')
           end if
        end if

      end function get_normal

  end subroutine construct_adv_interface_dg

  subroutine solve_vector_advection_dg_subcycle(field_name, state, velocity_name)
    !!< Construct and solve the advection equation for the given
    !!< field using discontinuous elements.
    
    !! Name of the field to be solved for.
    character(len=*), intent(in) :: field_name
    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state
    !! Optional velocity name
    character(len = *), intent(in) :: velocity_name

    !! Velocity to be solved for.
    type(vector_field), pointer :: U, U_old

    !! Coordinate and advecting velocity fields
    type(vector_field), pointer :: X, U_nl

    !! Velocity field in Cartesian coordinates for slope limiter
    type(vector_field) :: U_cartesian
    type(scalar_field) :: U_component

    !! Change in U over one timestep.
    type(vector_field) :: delta_U

    !! DG Courant number field
    type(scalar_field), pointer :: s_field

    !! Sparsity of advection matrix.    
    type(csr_sparsity), pointer :: sparsity
    
    !! System matrix.
    type(block_csr_matrix) :: matrix, mass, inv_mass

    !! Sparsity of mass matrix.
    type(csr_sparsity) :: mass_sparsity

    !! Right hand side vector.
    type(vector_field) :: rhs

    !! Whether to invoke the slope limiter
    logical :: limit_slope
    !! Which limiter to use
    integer :: limiter

    !! Number of advection subcycles.
    integer :: subcycles
    real :: max_courant_number

    character(len=FIELD_NAME_LEN) :: limiter_name
    integer :: i, j, dim

    U=>extract_vector_field(state, field_name)
    U_old=>extract_vector_field(state, "Old"//field_name)
    X=>extract_vector_field(state, "Coordinate")

    dim=mesh_dim(U)

    ! Reset U to value at the beginning of the timestep.
    call set(U, U_old)

    sparsity => get_csr_sparsity_firstorder(state, U%mesh, U%mesh)

    ! Add data space to the sparsity pattern.
    call allocate(matrix, sparsity, (/dim,dim/))

    mass_sparsity=make_sparsity_dg_mass(U%mesh)
    call allocate(mass, mass_sparsity, (/dim,dim/))
    call allocate(inv_mass, mass_sparsity, (/dim,dim/))
    call zero(mass)
    call zero(inv_mass)

    ! Ensure delta_U inherits options from U.
    call allocate(delta_U, U%dim, U%mesh, "delta_U")
    delta_U%option_path = U%option_path
    call allocate(rhs, U%dim, U%mesh, trim(field_name)//" RHS")
   
    call construct_vector_advection_dg(matrix, rhs, field_name, state, &
         mass, inv_mass, velocity_name=velocity_name)

!    call get_dg_inverse_mass_matrix(inv_mass, mass)
    
    ! Note that since theta and dt are module global, these lines have to
    ! come after construct_advection_diffusion_dg.
    call get_option("/timestepping/timestep", dt)
    
    if(have_option(trim(U%option_path)//&
         &"/prognostic/temporal_discretisation/discontinuous_galerkin/&
         &/number_advection_subcycles")) then
       call get_option(trim(U%option_path)//&
            &"/prognostic/temporal_discretisation/discontinuous_galerkin/&
            &/number_advection_subcycles", subcycles)
    else
       call get_option(trim(U%option_path)//&
            &"/prognostic/temporal_discretisation/discontinuous_galerkin/&
            &/maximum_courant_number_per_subcycle", Max_Courant_number)
       
       s_field => extract_scalar_field(state, "DG_CourantNumber")
       call calculate_diagnostic_variable(state, "DG_CourantNumber_Local", &
            & s_field)
       ewrite_minmax(s_field)
       
       subcycles = ceiling( maxval(s_field%val)/Max_Courant_number)
       call allmax(subcycles)
       ewrite(2,*) 'Number of subcycles for velocity eqn: ', subcycles
    end if

    limit_slope=.false.
    if (have_option(trim(U%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/slope_limiter")) then
       limit_slope=.true.
       
       ! Note unsafe for mixed element meshes
       if (element_degree(U,1)==0) then
          FLExit("Slope limiters make no sense for degree 0 fields")
       end if

       call get_option(trim(U%option_path)//"/prognostic/spatial_discretisation/&
            &discontinuous_galerkin/slope_limiter/name",limiter_name)

       select case(trim(limiter_name))
       case("Cockburn_Shu")
          limiter=LIMITER_COCKBURN
       case("Hermite_Weno")
          limiter=LIMITER_HERMITE_WENO
       case("minimal")
          limiter=LIMITER_MINIMAL
       case("FPN")
          limiter=LIMITER_FPN
       case("Vertex_Based")
          limiter=LIMITER_VB
       case default
          FLAbort('No such limiter')
       end select
       
    end if

    U_nl=>extract_vector_field(state, velocity_name)
    ! Allocate and zero cartesian velocity field
    call allocate(U_cartesian, X%dim, U_nl%mesh, "U_cartesian")
    call zero(U_cartesian)

    do i=1, subcycles

       print*, 'subcycle=', i

       ! dU = Advection * U
       call mult(delta_U, matrix, U)
       ewrite_minmax(delta_U)
       ! dU = dU + RHS
       call addto(delta_U, RHS, -1.0)
       ewrite_minmax(delta_U)
       ! dU = M^(-1) dU
       call dg_apply_mass(inv_mass, delta_U)
       ewrite_minmax(delta_U)
       ewrite_minmax(inv_mass)
       
       ! U = U + dt/s * dU
       call addto(U, delta_U, scale=-dt/subcycles)
       call halo_update(U)
       ewrite_minmax(delta_U)
!       limit_slope=.true.
!       limiter=LIMITER_VB
       if (limit_slope) then
          call project_local_to_cartesian(X, U, U_cartesian)
          ! Filter wiggles from U
          do j=1,U_cartesian%dim
             U_component=extract_scalar_field(U_cartesian,j)
             call limit_slope_dg(U_component, U_nl, X, state, limiter)
          end do
          call project_cartesian_to_local(X, U_cartesian, U)
       end if

    end do

    call deallocate(delta_U)
    call deallocate(matrix)
    call deallocate(mass)
    call deallocate(inv_mass)
    call deallocate(mass_sparsity)
    call deallocate(rhs)

  end subroutine solve_vector_advection_dg_subcycle

  subroutine construct_vector_advection_dg(big_m, rhs, field_name,&
       & state, mass, inv_mass, velocity_name) 
    !!< Construct the advection equation for discontinuous elements in
    !!< acceleration form.
    !!< 
    !!< If mass is provided then the mass matrix is not added into big_m or
    !!< rhs. It is instead returned as mass. This may be useful for testing
    !!< or for solving equations otherwise than in acceleration form.

    !! Main advection matrix.    
    type(block_csr_matrix), intent(inout) :: big_m
    !! Right hand side vector.
    type(vector_field), intent(inout) :: rhs
    
    !! Name of the field to be advected.
    character(len=*), intent(in) :: field_name
    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state
    !! Optional separate mass matrix.
    type(block_csr_matrix), intent(inout), optional :: mass, inv_mass
    !! Optional velocity name
    character(len = *), intent(in), optional :: velocity_name

    !! Position, velocity and advecting velocity fields.
    type(vector_field) :: X, U, U_nl

    !! Local velocity name
    character(len = FIELD_NAME_LEN) :: lvelocity_name

    !! Element index
    integer :: ele

    !! Status variable for field extraction.
    integer :: stat

    !! Field over the entire surface mesh containing bc values:
    type(vector_field) :: bc_value
    !! Integer array of all surface elements indicating bc type
    !! (see below call to get_entire_boundary_condition):
    integer, dimension(:,:), allocatable :: bc_type

    ewrite(1,*) "Writing advection equation for "&
         &//trim(field_name)

    ! These names are based on the CGNS SIDS.
    U=extract_vector_field(state, field_name)
    X=extract_vector_field(state, "Coordinate")

    if (.not.U%mesh%shape%degree==1) then
       FLAbort('This bit only works for linear elements at the moment')
    end if

    if(present(velocity_name)) then
      lvelocity_name = velocity_name
    else
      lvelocity_name = "NonlinearVelocity"
    end if

    U_nl=extract_vector_field(state, lvelocity_name)
    call incref(U_nl)

    flux_scheme=UPWIND_FLUX

    include_mass = .not. have_option(trim(U%option_path)//&
           "/prognostic/spatial_discretisation/discontinuous_galerkin/mass_terms/exclude_mass_terms")
           
    assert(has_faces(X%mesh))
    assert(has_faces(U%mesh))
    
    ! Enquire about boundary conditions we're interested in
    ! Returns an integer array bc_type over the surface elements
    ! that indicates the bc type (in the order we specified, i.e.
    ! BCTYPE_WEAKDIRICHLET=1)
    allocate( bc_type(1:U%dim,1:surface_element_count(U)) )
    call get_entire_boundary_condition(U, &
       & (/"weakdirichlet"/), &
       & bc_value, bc_type)

    call zero(big_m)
    call zero(RHS)
    if (present(mass)) call zero(mass)

    element_loop: do ele=1,element_count(U)
       
       call construct_vector_adv_element_dg(ele, big_m, rhs,&
            & X, U, U_nl, &
            bc_value, bc_type, mass, inv_mass)
       
    end do element_loop
    
    ! Drop any extra field references.

    call deallocate(U_nl)
    call deallocate(bc_value)

  end subroutine construct_vector_advection_dg

  subroutine construct_vector_adv_element_dg(ele, big_m, rhs,&
       & X, U, U_nl, &
       & bc_value, bc_type, &
       & mass, inv_mass)
    !!< Construct the advection_diffusion equation for discontinuous elements in
    !!< acceleration form.
    implicit none
    !! Index of current element
    integer :: ele
    !! Main advection matrix.
    type(block_csr_matrix), intent(inout) :: big_m
    !! Right hand side vector.
    type(vector_field), intent(inout) :: rhs
    !! Field over the entire surface mesh containing bc values:
    type(vector_field), intent(in):: bc_value
    !! Integer array of all surface elements indicating bc type
    !! (see above call to get_entire_boundary_condition):
    integer, dimension(:,:), intent(in):: bc_type
    !! Optional separate mass matrix.
    type(block_csr_matrix), intent(inout), optional :: mass, inv_mass
    
    !! Position, velocity and advecting velocity.
    type(vector_field), intent(in) :: X, U, U_nl

    ! Bilinear forms.
    real, dimension(mesh_dim(U),mesh_dim(U),ele_loc(U,ele),ele_loc(U,ele)) :: mass_mat
    real, dimension(mesh_dim(U),mesh_dim(U),ele_loc(U,ele),ele_loc(U,ele)) :: &
         Advection_mat, Ad_mat1, Ad_mat2

    ! Local assembly matrices.
    real, dimension(mesh_dim(U),mesh_dim(U),ele_loc(U,ele), ele_loc(U,ele)) :: l_U_mat
    real, dimension(mesh_dim(U)*ele_loc(U,ele), mesh_dim(U)*ele_loc(U,ele)) :: l_mass
    real, dimension(ele_loc(U,ele)) :: l_U_rhs

    ! Local variables.
    
    ! Neighbour element, face and neighbour face.
    integer :: ele_2, face, face_2
    ! Loops over faces.
    integer :: ni
    
    ! Transform from local to physical coordinates.
    real, dimension(ele_ngi(X,ele)) :: detwei, detJ
    real, dimension(mesh_dim(U), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(mesh_dim(U), mesh_dim(U), ele_ngi(X,ele)) :: G

    ! Different velocities at quad points.
    real, dimension(U_nl%dim, ele_ngi(U_nl, ele)) :: U_quad
    real, dimension(U_nl%dim, ele_ngi(U_nl, ele)) :: U_nl_q
    real, dimension(ele_ngi(U_nl, ele)) :: U_nl_div_q

    ! Node and shape pointers.
    integer, dimension(:), pointer :: U_ele
    type(element_type), pointer :: U_shape, U_nl_shape
    ! Neighbours of this element.
    integer, dimension(:), pointer :: neigh, x_neigh

    integer :: i, gi, dim, dim1, dim2, nloc

    logical :: boundary_element

    dim=mesh_dim(U)

    !----------------------------------------------------------------------
    ! Establish local node lists
    !----------------------------------------------------------------------
    
    U_ele=>ele_nodes(U,ele)  ! Velocity node numbers

    !----------------------------------------------------------------------
    ! Establish local shape functions
    !----------------------------------------------------------------------

    U_shape=>ele_shape(U, ele)
    U_nl_shape=>ele_shape(U_nl, ele)

    !==========================
    ! Coordinates
    !==========================

    ! Get J and detJ
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), detwei=detwei, J=J, detJ=detJ)

    !----------------------------------------------------------------------
    ! Construct bilinear forms.
    !----------------------------------------------------------------------

    ! Element mass matrix.
    !  /
    !  | W G U  dV
    !  / 
    do gi=1,ele_ngi(X,ele)
       G(:,:,gi)=matmul(J(:,:,gi), transpose(J(:,:,gi)))/detJ(gi)
    end do

    mass_mat=shape_shape_tensor(u_shape, u_shape, &
         u_shape%quadrature%weight, G)

    ! Advecting velocity at quadrature points.
    U_quad=ele_val_at_quad(U_nl,ele)
    U_nl_q=U_quad
    U_nl_div_q=ele_div_at_quad(U_nl, ele, U_shape%dn)

    ! Element advection matrix
    !    /                           /
    !  - | (grad W dot U_nl) T dV -  | W ( div U_nl ) G W dV
    !    /                           /

    Ad_mat1 = -dshape_dot_vector_tensor_shape(U_shape%dn, U_nl_q, G, U_shape, U_shape%quadrature%weight)
!    print*, 'U_shape%quadrature%weight'
!    print*, U_shape%quadrature%weight
!    print*, 'U_nl_q'
!    print*, U_nl_q
!    print*, 'Ad_mat1'
!    print*, Ad_mat1
    Ad_mat2 = -shape_shape_tensor(U_shape, U_shape, U_nl_div_q * U_shape%quadrature%weight, G)
!    print*, 'U_nl_div_q'
!    print*, U_nl_div_q
!    print*, 'Ad_mat2'
!    print*, Ad_mat2
    Advection_mat = -dshape_dot_vector_tensor_shape(U_shape%dn, U_nl_q, G, U_shape, U_shape%quadrature%weight)  &
           -shape_shape_tensor(U_shape, U_shape, U_nl_div_q * U_shape%quadrature%weight, G)

   !----------------------------------------------------------------------
   ! Perform global assembly.
   !----------------------------------------------------------------------

   l_U_rhs=0.0

   ! Assemble matrix.
    
   ! Advection.
   l_U_mat= Advection_mat

   if (present(mass)) then
      ! Return mass separately.
      do dim1 = 1,dim
         do dim2 = 1,dim
            call addto(mass, dim1, dim2, U_ele, U_ele, mass_mat(dim1,dim2,:,:))
         end do
      end do
      ! calculate inverse_mass here too for now
      nloc=ele_loc(U,ele)
      do dim1 = 1,dim
         do dim2 = 1,dim
            l_mass(nloc*(dim1-1)+1:nloc*dim1, nloc*(dim2-1)+1:nloc*dim2)=&
                 mass_mat(dim1,dim2,:,:)
         end do
      end do
!      print*, 'l_mass'
!      print*, l_mass
      call invert(l_mass)
!      print*, 'inverse l_mass'
!      print*, l_mass
      do dim1 = 1,dim
         do dim2 = 1,dim
            call addto(inv_mass, dim1, dim2, U_ele, U_ele, &
                 l_mass(nloc*(dim1-1)+1:nloc*dim1,nloc*(dim2-1)+1:nloc*dim2))
         end do
      end do
      ewrite_minmax(inv_mass)
   else
      if(include_mass) then
         ! Put mass in the matrix.
         l_U_mat=l_U_mat+mass_mat
      end if
   end if

   do dim1 = 1,dim
      do dim2 = 1,dim
         call addto(big_m, dim1, dim2, U_ele, U_ele, l_U_mat(dim1,dim2,:,:))
      end do
   end do

   !-------------------------------------------------------------------
   ! Interface integrals
   !-------------------------------------------------------------------
    
   neigh=>ele_neigh(U, ele)
   ! x_neigh/=t_neigh only on periodic boundaries.
   x_neigh=>ele_neigh(X, ele)

   ! Flag for whether this is a boundary element.
   boundary_element=.false.

   neighbourloop: do ni=1,size(neigh)

      !----------------------------------------------------------------------
      ! Find the relevant faces.
      !----------------------------------------------------------------------
       
      ! These finding routines are outside the inner loop so as to allow
      ! for local stack variables of the right size in
      ! construct_add_diff_interface_dg.

      ele_2=neigh(ni)
       
      ! Note that although face is calculated on field U, it is in fact
      ! applicable to any field which shares the same mesh topology.
      face=ele_face(U, ele, ele_2)
    
      if (ele_2>0) then
         ! Internal faces.
         face_2=ele_face(U, ele_2, ele)
      else
         ! External face.
         face_2=face
         boundary_element=.true.
      end if

      call construct_vector_adv_interface_dg(ele, face, face_2,&
           & big_m, rhs, X, U, U_nl,&
           & bc_value, bc_type)

   end do neighbourloop
    
 end subroutine construct_vector_adv_element_dg
  
  subroutine construct_vector_adv_interface_dg(ele, face, face_2, &
       big_m, rhs, &
       & X, U, U_nl,&
       & bc_value, bc_type)

    !!< Construct the DG element boundary integrals on the ni-th face of
    !!< element ele.
    implicit none

    integer, intent(in) :: ele, face, face_2
    type(block_csr_matrix), intent(inout) :: big_m
    type(vector_field), intent(inout) :: rhs
    ! We pass these additional fields to save on state lookups.
    type(vector_field), intent(in) :: X, U, U_nl
   !! Field over the entire surface mesh containing bc values:
    type(vector_field), intent(in):: bc_value
    !! Integer array of all surface elements indicating bc type
    !! (see above call to get_entire_boundary_condition):
    integer, dimension(:,:), intent(in):: bc_type

    ! Face objects and numberings.
    type(element_type), pointer :: U_shape, U_shape_2
    integer, dimension(face_loc(U,face)) :: U_face, U_face_l
    integer, dimension(face_loc(U,face_2)) :: U_face_2
    integer, dimension(face_loc(U_nl,face)) :: U_nl_face
    integer, dimension(face_loc(U_nl,face_2)) :: U_nl_face_2

    ! Note that both sides of the face can be assumed to have the same
    ! number of quadrature points.
    real, dimension(U_nl%dim, face_ngi(U_nl, face)) :: n1, n2, normal, U_nl_q,&
         & U_nl_f_q, U_nl_f2_q
    logical, dimension(face_ngi(U_nl, face)) :: inflow
    real, dimension(face_ngi(U_nl, face)) :: U_nl_q_dotn1, U_nl_q_dotn2, U_nl_q_dotn, income
    ! Variable transform times quadrature weights.
    real, dimension(ele_ngi(X,ele)) :: detwei, detJ
    real, dimension(face_ngi(X,face)) :: detwei_f
    real :: detJ_f
    real, dimension(mesh_dim(U), X%dim) :: J_f
    real, dimension(mesh_dim(U), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(mesh_dim(U), mesh_dim(U), face_ngi(X,face)) :: G
    real, dimension(face_ngi(U,face)) :: inner_advection_integral, outer_advection_integral

    ! Bilinear forms
    real, dimension(mesh_dim(U),mesh_dim(U),face_loc(U,face),face_loc(U,face)) :: nnAdvection_out
    real, dimension(mesh_dim(U),mesh_dim(U),face_loc(U,face),face_loc(U,face_2)) :: nnAdvection_in

    integer :: dim, dim1, dim2, gi, i
    integer, parameter :: cyc3(1:5)=(/ 1, 2, 3, 1, 2 /)

    dim=mesh_dim(U)

    U_face=face_global_nodes(U, face)
    U_face_l=face_local_nodes(U, face)
    U_shape=>face_shape(U, face)

    U_face_2=face_global_nodes(U, face_2)
    U_shape_2=>face_shape(U, face_2)
    
    !Unambiguously calculate the normal using the face with the higher
    !face number. This is so that the normal is identical on both sides.
    ! Jemma: need to be more careful here - actually have to calculate normal for face2

    n1=get_normal(local_face_number(U%mesh,face))
    n2=get_normal(local_face_number(U%mesh,face_2))
    
    !----------------------------------------------------------------------
    ! Construct element-wise quantities.
    !----------------------------------------------------------------------

    ! Advecting velocity at quadrature points.
    U_nl_f_q = face_val_at_quad(U_nl, face)
    U_nl_f2_q = face_val_at_quad(U_nl, face_2)

    if(local_face_number(U%mesh, face)==3) then
       U_nl_f_q=sqrt(2.)*U_nl_f_q
    end if
    if(local_face_number(U%mesh, face_2)==3) then
       U_nl_f2_q=sqrt(2.)*U_nl_f2_q
    end if

    u_nl_q_dotn1 = sum(U_nl_f_q*n1,1)
    u_nl_q_dotn2 = -sum(U_nl_f2_q*n2,1)
    U_nl_q_dotn=0.5*(u_nl_q_dotn1+u_nl_q_dotn2)
!    print*, 'U_nl_q_dotn'
!    print*, U_nl_q_dotn

    ! Inflow is true if the flow at this gauss point is directed
    ! into this element.
    inflow = u_nl_q_dotn<0.0
    income = merge(1.0,0.0,inflow)

    ! Calculate tensor G on face
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), detwei=detwei, J=J, detJ=detJ)
    J_f(:,:)=J(:,:,1)
    detJ_f=0.0
    select case (mesh_dim(U))
    case(1)
       detJ_f=1.0
    case(2)
       detJ_f = sqrt(J_f(1,1)**2 + J_f(2,1)**2)
    end select
    do gi=1,face_ngi(X,ele)
       G(:,:,gi)=matmul(J_f(:,:), transpose(J_f(:,:)))/detJ_f
    end do
!    print*, 'G'
!    print*, G

    !----------------------------------------------------------------------
    ! Construct bilinear forms.
    !----------------------------------------------------------------------
       
    ! Calculate outflow boundary integral.
    ! can anyone think of a way of optimising this more to avoid
    ! superfluous operations (i.e. multiplying things by 0 or 1)?

    ! first the integral around the inside of the element
    ! (this is the flux *out* of the element)
    inner_advection_integral = (1.-income)*u_nl_q_dotn
    nnAdvection_out=shape_shape_tensor(U_shape, U_shape,  &
         &            inner_advection_integral*U_shape%quadrature%weight, G)
!    print*, 'js:nnAdvection_out:'
!    print*, nnAdvection_out
    ! now the integral around the outside of the element
    ! (this is the flux *in* to the element)
    outer_advection_integral = income*u_nl_q_dotn
    nnAdvection_in=shape_shape_tensor(U_shape, U_shape_2, &
         &            outer_advection_integral*U_shape%quadrature%weight, G)
!    print*, 'js:nnAdvection_in:'
!    print*, nnAdvection_in
       
    !----------------------------------------------------------------------
    ! Perform global assembly.
    !----------------------------------------------------------------------

    ! Insert advection in matrix.
    
    ! Outflow boundary integral.
    do dim1 = 1, dim
       do dim2 = 1, dim
          call addto(big_M, dim1, dim2, U_face, U_face, nnAdvection_out(dim1,dim2,:,:))
       end do
    end do

    ! Inflow boundary integral.
    do dim1 = 1, dim
       do dim2 = 1, dim
          call addto(big_M, dim1, dim2, U_face, U_face_2, nnAdvection_in(dim1,dim2,:,:))
       end do
    end do

    contains

      function get_normal(face) result(norm)

        integer, intent(in) :: face
        real, dimension(U_nl%dim, face_ngi(U_nl,face)) :: norm

        integer :: i

        if(U_nl%dim==1) then
           if(face==1) then
              forall(i=1:face_ngi(U_nl,face)) norm(1,i)=1.
           else if(face==2) then
              forall(i=1:face_ngi(U_nl,face)) norm(1,i)=-1.
           end if
        else if(U_nl%dim==2) then
           if(face==1) then
              forall(i=1:face_ngi(U_nl,face)) norm(1:2,i)=(/-1.,0./)
           else if(face==2) then
              forall(i=1:face_ngi(U_nl,face)) norm(1:2,i)=(/0.,-1./)
           else if(face==3) then
              forall(i=1:face_ngi(U_nl,face)) norm(1:2,i)=(/1/sqrt(2.),1/sqrt(2.)/)
           else
              FLAbort('Oh dear oh dear')
           end if
        end if

      end function get_normal

    end subroutine construct_vector_adv_interface_dg

end module advection_local_DG
