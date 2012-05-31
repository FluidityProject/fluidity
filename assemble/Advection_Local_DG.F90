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
  !!< equation when vector fields are represented in the local coordinates
  !!< from the Piola transformation.
  use elements
  use global_parameters, only:current_debug_level, OPTION_PATH_LEN
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
  use Tensors
  use petsc_solve_state_module
  use spud
  use slope_limiters_dg
  use sparsity_patterns
  use sparse_matrices_fields
  use sparsity_patterns_meshes
  use Vector_Tools
  use manifold_tools
  use bubble_tools
  use diagnostic_fields, only: calculate_diagnostic_variable
  use global_parameters, only : FIELD_NAME_LEN

  implicit none

  ! Buffer for output messages.
  character(len=255), private :: message

  private
  public solve_advection_dg_subcycle, solve_vector_advection_dg_subcycle,&
       & solve_advection_cg_tracer, compute_U_residual, get_pv

  !parameters specifying vector upwind options
  integer, parameter :: VECTOR_UPWIND_EDGE=1, VECTOR_UPWIND_SPHERE=2

  ! Local private control parameters. These are module-global parameters
  ! because it would be expensive and/or inconvenient to re-evaluate them
  ! on a per-element or per-face basis
  real :: dt

  ! Whether to include various terms
  logical :: include_advection
  contains

    subroutine solve_advection_dg_subcycle(field_name, state, velocity_name,&
         &continuity, flux)
    !!< Construct and solve the advection equation for the given
    !!< field using discontinuous elements.
    
    !! Name of the field to be solved for.
    character(len=*), intent(in) :: field_name
    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state
    !! Optional velocity name
    character(len = *), intent(in) :: velocity_name
    !! Solve continuity equation as opposed to advection equation
    logical, intent(in), optional :: continuity
    !! Return a flux variable such that T-T_old=div(Flux)
    type(vector_field), intent(inout), optional :: Flux

    !! Tracer to be solved for.
    type(scalar_field), pointer :: T, T_old, s_field

    !! Coordinate and velocity fields
    type(vector_field), pointer :: X, U_nl

    !! Velocity field in Cartesian coordinates for slope limiter
    type(vector_field) :: U_nl_cartesian

    !! Change in T over one timestep.
    type(scalar_field) :: delta_T, delta_T_total

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

    !! Flux reconstruction stuff:
    !! Upwind fluxes trace variable
    type(scalar_field) :: UpwindFlux
    !! Mesh where UpwindFlux lives
    type(mesh_type), pointer :: lambda_mesh
    !! Upwind Flux matrix
    type(csr_matrix) :: UpwindFluxMatrix

    character(len=FIELD_NAME_LEN) :: limiter_name
    integer :: i, ele

    ewrite(1,*) 'subroutine solve_advection_dg_subcycle'


    T=>extract_scalar_field(state, field_name)
    T_old=>extract_scalar_field(state, "Old"//field_name)
    X=>extract_vector_field(state, "Coordinate")
    U_nl=>extract_vector_field(state, velocity_name)

    ewrite(2,*) 'UVALS in dg', maxval(abs(U_nl%val))

    ! Reset T to value at the beginning of the timestep.
    call set(T, T_old)
    if(present(Flux)) then
       call zero(Flux)
    end if

    sparsity => get_csr_sparsity_firstorder(state, T%mesh, T%mesh)
    call allocate(matrix, sparsity) ! Add data space to the sparsity
    ! pattern.
    call zero(matrix)
    
    !Flux reconstruction allocations
    if(present(Flux)) then
       !Allocate trace variable for upwind fluxes
       lambda_mesh=>extract_mesh(state, "VelocityMeshTrace")
       call allocate(UpwindFlux,lambda_mesh, "UpwindFlux")
       call zero(UpwindFlux)
       !Allocate matrix that maps T values to upwind fluxes
       sparsity => get_csr_sparsity_firstorder(state,lambda_mesh,T%mesh)
       call allocate(UpwindFluxMatrix,sparsity)
       call zero(UpwindFluxMatrix)
    end if

    mass_sparsity=make_sparsity_dg_mass(T%mesh)
    call allocate(mass, mass_sparsity)
    call zero(mass)
    call allocate(inv_mass, mass_sparsity)
    call zero(inv_mass)

    ! Ensure delta_T inherits options from T.
    call allocate(delta_T, T%mesh, "delta_T")
    call zero(delta_T)
    delta_T%option_path = T%option_path
    if(present(Flux)) then
       call allocate(delta_T_total, T%mesh, "deltaT_total")
       call zero(delta_T_total)
    end if
    call allocate(rhs, T%mesh, trim(field_name)//" RHS")
    call zero(rhs)

    if(present(Flux)) then
       call construct_advection_dg(matrix, rhs, field_name, state, &
            mass, velocity_name=velocity_name,continuity=continuity, &
            UpwindFlux=upwindflux, UpwindFluxMatrix=UpwindFluxMatrix)
    else
       call construct_advection_dg(matrix, rhs, field_name, state, &
            mass, velocity_name=velocity_name,continuity=continuity)
    end if

    call get_dg_inverse_mass_matrix(inv_mass, mass)
    
    ! Note that since dt is module global, these lines have to
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
            &maximum_courant_number_per_subcycle", Max_Courant_number)
       
       s_field => extract_scalar_field(state, "DG_CourantNumber")
       call calculate_diagnostic_variable(state, "DG_CourantNumber_Local", &
            & s_field)
       ewrite_minmax(s_field)
       subcycles = ceiling( maxval(s_field%val)/Max_Courant_number)
       call allmax(subcycles)
       ewrite(2,*) 'Number of subcycles for tracer eqn: ', subcycles
    end if

    limit_slope=.false.
    if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/slope_limiter")) then
       limit_slope=.true.

       call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation/discontinuous_galerkin/slope_limiter/name",limiter_name)

       select case(trim(limiter_name))
       case("Vertex_Based")
          limiter=LIMITER_VB
       case default
          FLAbort('No such limiter')
       end select
       
       ! Note unsafe for mixed element meshes
       if (element_degree(T,1)==0) then
          FLExit("Slope limiters make no sense for degree 0 fields")
       end if

    end if

    if (limit_slope) then
       ! Filter wiggles from T
       call limit_slope_dg(T, U_nl, X, state, limiter, delta_T)
       if(present(flux)) then
          call zero(upwindflux)
          call mult(delta_T_total, mass, delta_T)
          call update_flux(Flux, Delta_T_total, UpwindFlux)
       end if
   end if

    do i=1, subcycles
       ewrite(1,*) 'SUBCYCLE', i

       ! dT = Advection * T
       call mult(delta_T, matrix, T)
       ! dT = dT + RHS
       call addto(delta_T, RHS, -1.0)
       if(present(Flux)) then
          if(maxval(abs(RHS%val))>1.0e-8) then
             FLAbort('Flux reconstruction doesn''t work if diffusion/bcs present at the moment')
          end if
       end if
       call scale(delta_T, -dt/subcycles)
       ! dT = M^(-1) dT
       if(present(flux)) then
          call mult(UpwindFlux, upwindfluxmatrix, T)
          call scale(UpwindFlux,-dt/subcycles)
          call update_flux(Flux, Delta_T, UpwindFlux)
       end if
       ! T = T + dt/s * dT
       call dg_apply_mass(inv_mass, delta_T)
       ewrite(2,*) 'Delta_T', maxval(abs(delta_T%val))
       call addto(T, delta_T)
       !Probably need to halo_update(Flux) as well       
       call halo_update(T)

       if (limit_slope) then
          ! Filter wiggles from T
          call limit_slope_dg(T, U_nl, X, state, limiter, delta_T)
          if(present(flux)) then
             call mult(delta_t_total, mass, delta_t)
             call zero(UpwindFlux)
             call update_flux(Flux, Delta_T_total, UpwindFlux)
          end if
       end if
    end do

    if(present(Flux)) then
       if(have_option('/material_phase::Fluid/scalar_field::LayerThickness/p&
            &rognostic/spatial_discretisation/debug')) then
          do ele = 1, ele_count(Flux)
             call check_flux(Flux,T,T_old,X,ele)
          end do
       end if
       call deallocate(UpwindFlux)
       call deallocate(UpwindFluxMatrix)
       call deallocate(delta_T_total)
    end if

    call deallocate(delta_T)
    call deallocate(matrix)
    call deallocate(mass)
    call deallocate(inv_mass)
    call deallocate(mass_sparsity)
    call deallocate(rhs)

  end subroutine solve_advection_dg_subcycle

  subroutine check_flux(Flux,T,T_old,X,ele)
    type(vector_field), intent(in) :: Flux, X
    type(scalar_field), intent(in) :: T,T_old
    integer, intent(in) :: ele
    !
    real, dimension(Flux%dim,ele_loc(Flux,ele)) :: Flux_vals
    real, dimension(ele_ngi(T,ele)) :: Div_Flux_gi
    real, dimension(ele_ngi(T,ele)) :: Delta_T_gi
    real, dimension(ele_loc(T,ele)) :: Delta_T_rhs, Div_Flux_rhs
    integer :: dim1, loc
    real, dimension(ele_ngi(X,ele)) :: detwei
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J_mat

    Delta_T_gi = ele_val_at_quad(T,ele)-ele_val_at_quad(T_old,ele)
    Flux_vals = ele_val(Flux,ele)

    Div_Flux_gi = 0.
    !dn is loc x ngi x dim    
    do dim1 = 1, Flux%dim
       do loc = 1, ele_loc(Flux,ele)
          Div_Flux_gi = Div_Flux_gi + Flux%mesh%shape%dn(loc,:,dim1)*&
               &Flux_vals(dim1,loc)
       end do
    end do

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele),&
         detwei=detwei,J=J_mat)
    
    Delta_T_rhs = shape_rhs(ele_shape(T,ele),Delta_T_gi*detwei)
    Div_Flux_rhs = shape_rhs(ele_shape(T,ele),Div_Flux_gi*&
         Flux%mesh%shape%quadrature%weight)

    assert(maxval(abs(Delta_T_rhs-Div_Flux_rhs))<1.0e-10)

  end subroutine check_flux

  subroutine update_flux(Flux, Delta_T, UpwindFlux)
    !! Subroutine to compute div-conforming Flux such that
    !! div Flux = Delta_T
    !! with Flux.n = UpwindFlux on element boundaries
    !! and then to add to Flux
    type(vector_field), intent(inout) :: Flux
    type(scalar_field), intent(in) :: Delta_T
    type(scalar_field), intent(in) :: UpwindFlux
    !
    integer :: ele
    !Checking the Flux is in local representation
    assert(Flux%dim==mesh_dim(Flux))

    do ele = 1, ele_count(Flux)
       call update_flux_ele(Flux, Delta_T, UpwindFlux,ele)
    end do
    
  end subroutine update_flux

  subroutine update_flux_ele(Flux, Delta_T, UpwindFlux,ele)
    !! Subroutine to compute div-conforming Flux such that
    !! div Flux = Delta_T
    !! with Flux.n = UpwindFlux on element boundaries
    type(vector_field), intent(inout) :: Flux
    type(scalar_field), intent(in) :: Delta_T
    type(scalar_field), intent(in) :: UpwindFlux
    integer, intent(in) :: ele
    !
    real, dimension(Flux%dim,ele_loc(Flux,ele)) :: Flux_vals
    real, dimension(ele_loc(Delta_T,ele)) :: Delta_T_vals
    real, dimension(Flux%dim*ele_loc(Flux,ele),&
         Flux%dim*ele_loc(Flux,ele)) :: Flux_mat
    real, dimension(Flux%dim*ele_loc(Flux,ele)) :: Flux_rhs
    integer :: row, ni, ele2, face, face2,floc, i, dim1, dim2
    integer, dimension(:), pointer :: neigh
    type(constraints_type), pointer :: flux_constraint
    type(element_type), pointer :: flux_shape, T_shape
    real, dimension(Flux%dim,ele_loc(Delta_T,ele),&
         &ele_loc(Flux,ele)) :: grad_mat
    real, dimension(ele_loc(Delta_t,ele)) :: delta_t_val
    real, dimension(ele_loc(flux,ele),ele_loc(flux,ele)) &
         &:: flux_mass_mat
    real, dimension(ele_loc(Delta_t,ele)) :: div_flux_val
    real, dimension(ele_ngi(Delta_t,ele)) :: div_flux_gi
    real :: total_flux, residual
    !! Need to set up and solve equation for flux DOFs
    !! Rows are as follows:
    !! All of the normal fluxes through edges, numbered according to
    !! UpwindFlux numbering
    !! Inner product of Flux with grad basis equaling gradient of delta_T
    !! plus boundary conditions
    !! Inner product of solution with curl basis equals zero
    !! inner product of solution with constraint basis equals zero

    !! Debugging check:
    !! \int \Delta T dV = \int div Flux dV = \int Flux.n dS
    !!                                     = \int \hat{Flux}.\hat{n dS}
    

    neigh => ele_neigh(Flux,ele)
    total_flux = 0.
    do ni = 1, size(neigh)
       ele2 = neigh(ni)
       face = ele_face(Flux,ele,ele2)
       if(ele2>0) then
          face2 = ele_face(Flux,ele2,ele)
       else
          face2 = face
       end if
       if(face>face2) then
          total_flux=total_flux+sum(face_val(UpwindFlux,face))
       else
          total_flux=total_flux-sum(face_val(UpwindFlux,face))
       end if
    end do
    residual = abs(total_flux-sum(ele_val(Delta_T,ele)))&
         &/max(1.0,maxval(ele_val(Delta_T,ele)))

    assert(residual<1.0e-10)

    Flux_rhs = 0.
    Flux_mat = 0.
    flux_shape => ele_shape(flux,ele)
    T_shape => ele_shape(delta_t,ele)
    flux_constraint => flux%mesh%shape%constraints

    row = 1
    !first do the face DOFs
    neigh => ele_neigh(Flux,ele)
    do ni = 1, size(neigh)
       ele2 = neigh(ni)
       face = ele_face(Flux,ele,ele2)
       if(ele2>0) then
          face2 = ele_face(Flux,ele2,ele)
       else
          face2 = face
       end if
       floc = face_loc(upwindflux,face)
       
       call update_flux_face(&
            Flux,UpwindFlux,face,face2,&
            ele,ele2,Flux_mat(row:row+floc-1,:),&
            Flux_rhs(row:row+floc-1))
       row = row + floc
    end do
    
    !Next the gradient basis
    !Start with \nabla\cdot F = \Delta T
    !Integrate with test function over element
    ! <\phi,Div F>_E = <\phi, \Delta T>_E

    grad_mat = shape_dshape(T_shape, flux_shape%dn,&
           & T_shape%quadrature%weight)
    delta_t_val = ele_val(delta_t,ele)
    do i = 1, ele_loc(Delta_T,ele)-1
       do dim1 = 1, flux%dim
          Flux_mat(row,(dim1-1)*ele_loc(flux,ele)+1:dim1*ele_loc(flux,ele))&
               &=grad_mat(dim1,i,:)
       end do
       flux_rhs(row) = delta_t_val(i)
       row = row + 1
    end do
    
    !Then the curl basis
    !Adding in the curl terms
    flux_mass_mat = shape_shape(flux_shape,flux_shape,T_shape%quadrature%weight)
    do i = 1, flux_constraint%n_curl_basis
       flux_mat(row,:) = 0.
       flux_rhs(row) = 0.
       do dim1 = 1, mesh_dim(flux)
          do dim2 = 1, mesh_dim(flux)
             flux_mat(&
                  row,(dim2-1)*ele_loc(flux,ele)+1:dim2*ele_loc(flux,ele)) =&
                  flux_mat(&
                  row,(dim2-1)*ele_loc(flux,ele)+1:dim2*ele_loc(flux,ele)) +&
                  matmul(flux_constraint%curl_basis(i,:,dim1),&
                  flux_mass_mat)
          end do
       end do
       row = row + 1
    end do

    !Then the constraints
    do i = 1, flux_constraint%n_constraints
       do dim1 = 1, mesh_dim(flux)
          flux_mat(&
               row,(dim1-1)*ele_loc(flux,ele)+1:&
               dim1*ele_loc(flux,ele)) =&
               flux_constraint%orthogonal(i,:,dim1)
       end do
       row = row + 1
    end do

    call solve(flux_mat,flux_rhs)

    do dim1 = 1, mesh_dim(flux)
       flux_vals(dim1,:) = &
            flux_rhs((dim1-1)*ele_loc(flux,ele)+1:dim1*ele_loc(flux,ele))
    end do

    !A debugging check.
    !Computing < \phi, div F> in local coordinates
    !and comparing with <\phi, delta_T> in global coordinates
    !(works because of pullback formula)

    !dn is loc x ngi x dim
    !Flux is dim x loc
    div_flux_gi = 0.
    do dim1 =1, size(flux_vals,1)
       do i=1,size(flux_vals,2)
          div_flux_gi = div_flux_gi + flux_vals(dim1,i)*&
               &flux_shape%dn(i,:,dim1)
       end do
    end do
    div_flux_val = shape_rhs(T_shape,div_flux_gi*T_shape%quadrature%weight)
    residual = maxval(abs(div_flux_val-delta_T_val))/max(1.0&
         &,maxval(abs(delta_T_val)))
    assert(residual<1.0e-10)
    
    !Add the flux in
    call addto(flux,ele_nodes(flux,ele),flux_vals)
  end subroutine update_flux_ele

  subroutine update_flux_face(&
       Flux,UpwindFlux,face,face2,&
       ele,ele2,Flux_mat_rows,Flux_rhs_rows)
    type(vector_field), intent(in) :: Flux
    type(scalar_Field), intent(in) :: UpwindFlux
    integer, intent(in) :: face,face2,ele,ele2
    real, dimension(face_loc(upwindflux,ele),&
         &mesh_dim(flux)*ele_loc(flux,ele)), &
         &intent(inout) :: flux_mat_rows
    real, dimension(flux%mesh%shape%constraints%n_face_basis), intent(inout)&
         &:: flux_rhs_rows
    !
    integer :: dim1, row, flux_dim1
    real, dimension(mesh_dim(flux), face_ngi(flux, face)) :: n_local
    real, dimension(face_ngi(flux,face)) :: detwei_f
    real :: weight
    type(element_type), pointer :: flux_face_shape, upwindflux_face_shape
    integer, dimension(face_loc(flux,face)) :: flux_face_nodes
    real, dimension(mesh_dim(flux),face_loc(upwindflux,face),&
         face_loc(flux,face)) :: face_mat
    !
    flux_face_shape=>face_shape(flux, face)
    upwindflux_face_shape=>face_shape(upwindflux, face)
    flux_face_nodes=face_local_nodes(flux, face)

    !Get normal in local coordinates
    call get_local_normal(n_local,weight,flux,&
         &local_face_number(flux%mesh,face))

    detwei_f = weight*flux_face_shape%quadrature%weight
    face_mat = shape_shape_vector(&
         upwindflux_face_shape,flux_face_shape,detwei_f,n_local)
    !Equation is:
    ! <\phi,Flux.n> = <\phi,UpwindFlux> for all trace test functions \phi.

    do row = 1, size(flux_mat_rows,1)
       flux_mat_rows(row,:) = 0.
       do dim1 = 1,mesh_dim(flux)
          flux_dim1 = (dim1-1)*ele_loc(flux,ele)
          flux_mat_rows(row,flux_dim1+flux_face_nodes)&
               &=face_mat(dim1,row,:)
       end do
    end do
    !Need to adopt a sign convention:
    !If face>face_2
    !We store flux with normal pointing from face into face_2
    !Otherwise 
    !We store flux with normal pointing from face_2 into face
    ! Outflow boundary integral.
    if(face>face2) then
       flux_rhs_rows = face_val(UpwindFlux,face)
    else
       flux_rhs_rows = -face_val(UpwindFlux,face)
    end if
  end subroutine update_flux_face

  subroutine construct_advection_dg(big_m, rhs, field_name,&
       & state, mass, velocity_name, continuity,&
       & UpwindFlux, UpwindFluxMatrix)
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
    !! solve the continuity equation
    logical, intent(in), optional :: continuity
    !! Upwind flux field
    type(scalar_field), intent(in), optional :: UpwindFlux
    !! Upwind flux matrix
    type(csr_matrix), intent(inout), optional :: UpwindFluxMatrix

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
 
    assert(has_faces(X%mesh))
    assert(has_faces(T%mesh))
    
    call zero(big_m)
    call zero(RHS)
    call zero(mass)

    element_loop: do ele=1,element_count(T)
       
       call construct_adv_element_dg(ele, big_m, rhs,&
            & X, T, U_nl, mass, continuity=continuity,&
            & upwindflux=upwindflux, upwindfluxmatrix=upwindfluxmatrix)
       
    end do element_loop
    
    ! Drop any extra field references.

    call deallocate(U_nl)

  end subroutine construct_advection_dg

  subroutine construct_adv_element_dg(ele, big_m, rhs,&
       & X, T, U_nl, mass, continuity, upwindflux, upwindfluxmatrix)
    !!< Construct the advection_diffusion equation for discontinuous elements in
    !!< acceleration form.
    implicit none
    !! Index of current element
    integer :: ele
    !! Main advection matrix.
    type(csr_matrix), intent(inout) :: big_m
    !! Right hand side vector.
    type(scalar_field), intent(inout) :: rhs
    !! Optional separate mass matrix.
    type(csr_matrix), intent(inout) :: mass
    
    !! Position and velocity.
    type(vector_field), intent(in) :: X, U_nl

    type(scalar_field), intent(in) :: T

    !! Upwind flux field
    type(scalar_field), intent(in), optional :: UpwindFlux
    !! Upwind flux matrix
    type(csr_matrix), intent(inout), optional :: UpwindFluxMatrix


    !! Solve the continuity equation rather than advection
    logical, intent(in), optional :: continuity

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
    integer, dimension(:), pointer :: neigh

    integer :: gi,i,j

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
      !  - | (grad phi dot U_nl) T dV -| phi ( div U_nl ) T dV
      !    /                           /
      Advection_mat = -dshape_dot_vector_shape(T_shape%dn, U_nl_q, T_shape,&
           & T_shape%quadrature%weight)
      if(.not.present_and_true(continuity)) then
         Advection_mat = Advection_mat&
              &-shape_shape(T_shape, T_shape, U_nl_div_q * T_shape&
              &%quadrature%weight)
      end if
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

   call addto(mass, t_ele, t_ele, mass_mat)

   call addto(big_m, t_ele, t_ele, l_T_mat)

   !-------------------------------------------------------------------
   ! Interface integrals
   !-------------------------------------------------------------------
    
   neigh=>ele_neigh(T, ele)

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
      end if

      call construct_adv_interface_dg(ele, face, face_2,&
           & big_m, rhs, X, T, U_nl, upwindflux, upwindfluxmatrix)

   end do neighbourloop
    
 end subroutine construct_adv_element_dg
  
  subroutine construct_adv_interface_dg(ele, face, face_2, &
       big_m, rhs, X, T, U_nl,upwindflux, upwindfluxmatrix)

    !!< Construct the DG element boundary integrals on the ni-th face of
    !!< element ele.
    implicit none

    integer, intent(in) :: ele, face, face_2
    type(csr_matrix), intent(inout) :: big_m
    type(scalar_field), intent(inout) :: rhs
    ! We pass these additional fields to save on state lookups.
    type(vector_field), intent(in) :: X, U_nl
    type(scalar_field), intent(in) :: T
    !! Upwind flux field
    type(scalar_field), intent(in), optional :: UpwindFlux
    !! Upwind flux matrix
    type(csr_matrix), intent(inout), optional :: UpwindFluxMatrix

    ! Face objects and numberings.
    type(element_type), pointer :: T_shape, T_shape_2
    integer, dimension(face_loc(T,face)) :: T_face
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

    ! Normal weights
    real :: w1, w2
    integer :: i

    T_face=face_global_nodes(T, face)
    T_shape=>face_shape(T, face)

    T_face_2=face_global_nodes(T, face_2)
    T_shape_2=>face_shape(T, face_2)
    
    !Unambiguously calculate the normal using the face with the higher
    !face number. This is so that the normal is identical on both sides.
    ! Jemma: need to be more careful here - actually have to calculate normal for face2

    call get_local_normal(n1, w1, U_nl, local_face_number(U_nl%mesh,face))
    call get_local_normal(n2, w2, U_nl, local_face_number(U_nl%mesh,face_2))
    
    !----------------------------------------------------------------------
    ! Construct element-wise quantities.
    !----------------------------------------------------------------------

    if (include_advection) then
       
       ! Advecting velocity at quadrature points.
       U_f_q = face_val_at_quad(U_nl, face)
       U_f2_q = face_val_at_quad(U_nl, face_2)

       u_nl_q_dotn1 = sum(U_f_q*w1*n1,1)
       u_nl_q_dotn2 = -sum(U_f2_q*w2*n2,1)
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

       ! now the integral around the outside of the element
       ! (this is the flux *in* to the element)
       outer_advection_integral = income*u_nl_q_dotn
       nnAdvection_in=shape_shape(T_shape, T_shape_2, &
            &                       outer_advection_integral*T_shape%quadrature%weight)

    end if

    !----------------------------------------------------------------------
    ! Perform global assembly.
    !----------------------------------------------------------------------

    ! Insert advection in matrix.
    if (include_advection) then
    
       ! Outflow boundary integral.
       call addto(big_M, T_face, T_face,&
            nnAdvection_out)
       
       ! Inflow boundary integral.
       call addto(big_M, T_face, T_face_2,&
            nnAdvection_in)

       if(present(UpwindFlux).and.present(UpwindFluxMatrix)) then
          call construct_upwindflux_interface(upwindflux,upwindfluxmatrix)
       end if
    end if

    contains 
      subroutine construct_upwindflux_interface(upwindflux,upwindfluxmatrix)
        type(scalar_field), intent(in) :: upwindflux
        type(csr_matrix), intent(inout) :: upwindfluxmatrix
        !
        ! Bilinear forms
        real, dimension(face_loc(upwindflux,face),&
             face_loc(T,face)) :: upwindflux_mat_in,upwindflux_mat_out
        type(element_type), pointer :: flux_shape
        integer, dimension(face_loc(upwindflux,face)) :: flux_face
        integer :: sign

        flux_shape=>face_shape(upwindflux, face)
        flux_face=face_global_nodes(upwindflux, face)

       upwindflux_mat_in=shape_shape(flux_shape, T_shape,  &
            income*u_nl_q_dotn*T_shape%quadrature%weight)
       upwindflux_mat_out=shape_shape(flux_shape, T_shape,  &
            (1.0-income)*u_nl_q_dotn*T_shape%quadrature%weight)

       !Need to adopt a sign convention:
       !If face>face_2
       !We store flux with normal pointing from face into face_2
       !Otherwise 
       !We store flux with normal pointing from face_2 into face
       ! Outflow boundary integral.

       !inflow = u_nl_q_dotn<0.0
       !income = 1 if(inflow) and 0 otherwise
       if(face>face_2) then
          sign = 1
       else
          sign = -1
       end if

       if(face.ne.face_2) then
          call addto(upwindfluxmatrix, flux_face, T_face_2,&
               0.5*sign*upwindflux_mat_in)
       end if
       call addto(upwindfluxmatrix, flux_face, T_face,&
            0.5*sign*upwindflux_mat_out)

      end subroutine construct_upwindflux_interface

  end subroutine construct_adv_interface_dg  

  subroutine solve_advection_cg_tracer(Q,D,D_old,Flux,QF,state)
    !!< Solve the continuity equation for D*Q with Flux
    !!< d/dt (D*Q) + div(Flux*Q) = diffusion terms.
    !!< Done on the vorticity mesh
    !!< Return a PV flux Q defined on the velocity mesh
    !!< which is the projection of Flux*Q into the velocity space
    !!< Note that Flux and Q contain a factor of dt for convenience.
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(in),target :: D,D_old
    type(scalar_field), intent(inout),target :: Q
    type(vector_field), intent(in) :: Flux
    type(vector_field), intent(inout) :: QF
    !
    type(csr_sparsity), pointer :: Q_sparsity
    type(csr_matrix) :: Adv_mat
    type(scalar_field) :: Q_rhs, Qtest1,Qtest2
    type(vector_field), pointer :: X, down
    type(scalar_field), pointer :: Q_old
    integer :: ele
    real :: dt, t_theta

    ewrite(1,*) '  subroutine solve_advection_cg_tracer('

    X=>extract_vector_field(state, "Coordinate")
    down=>extract_vector_field(state, "GravityDirection")
    Q_old=>extract_scalar_field(state, "Old"//trim(Q%name))
    call get_option('/timestepping/timestep',dt)
    call get_option('/timestepping/theta',t_theta)
    
    !set up matrix and rhs
    Q_sparsity => get_csr_sparsity_firstorder(state, Q%mesh, Q%mesh)
    call allocate(adv_mat, Q_sparsity) ! Add data space to the sparsity
    ! pattern.
    call zero(adv_mat)
    call allocate(Q_rhs,Q%mesh,trim(Q%name)//"RHS")
    call zero(Q_rhs)
    call zero(QF)

    do ele = 1, ele_count(Q)
       call construct_advection_cg_tracer_ele(Q_rhs,adv_mat,Q,D,D_old,Flux&
            &,X,down,dt,t_theta,ele)
    end do

    ewrite(2,*) 'Q_RHS', maxval(abs(Q_rhs%val))

    call petsc_solve(Q,adv_mat,Q_rhs)

    !! Compute the PV flux to pass to velocity equation
    do ele = 1, ele_count(Q)
       call construct_pv_flux_ele(QF,Q,Q_old,D,D_old,Flux,X,down,t_theta,ele)
    end do

    if(have_option('/material_phase::Fluid/scalar_field::PotentialVorticity/&
         &prognostic/debug')) then
       call allocate(Qtest1,Q%mesh,trim(Q%name)//"test1")
       call zero(Qtest1)
       call allocate(Qtest2,Q%mesh,trim(Q%name)//"test2")
       call zero(Qtest2)
       do ele = 1, ele_count(Q)
          call test_pv_flux_ele(Qtest1,Qtest2,QF,Q,Q_old,D,D_old,&
               Flux,X,down,t_theta,ele)
       end do
       ewrite(2,*) 'Error = ', maxval(abs(Qtest1%val-Qtest2%val))
       assert(maxval(abs(Qtest1%val-Qtest2%val))<1.0e-10)
       ewrite(2,*) 'test passed'
       call deallocate(Qtest1)
       call deallocate(Qtest2)
    end if

    call deallocate(adv_mat)
    call deallocate(Q_rhs)

    ewrite(1,*) 'END  subroutine solve_advection_cg_tracer('

  end subroutine solve_advection_cg_tracer
  
  subroutine construct_advection_cg_tracer_ele(Q_rhs,adv_mat,Q,D,D_old,Flux,&
       & X,down,dt,t_theta,ele)
    type(scalar_field), intent(in) :: D,D_old,Q
    type(scalar_field), intent(inout) :: Q_rhs
    type(csr_matrix), intent(inout) :: Adv_mat
    type(vector_field), intent(in) :: X, Flux, down
    integer, intent(in) :: ele
    real, intent(in) :: dt, t_theta
    !
    real, dimension(ele_loc(Q,ele),ele_loc(Q,ele)) :: l_adv_mat
    real, dimension(ele_loc(Q,ele)) :: l_rhs
    real, dimension(Flux%dim,ele_ngi(Flux,ele)) :: Flux_gi
    real, dimension(ele_ngi(X,ele)) :: detwei, Q_gi, D_gi, D_old_gi
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    type(element_type), pointer :: Q_shape
    integer :: loc

    Q_shape => ele_shape(Q,ele)
    D_gi = ele_val_at_quad(D,ele)
    D_old_gi = ele_val_at_quad(D_old,ele)
    Q_gi = ele_val_at_quad(Q,ele)
    Flux_gi = ele_val_at_quad(Flux,ele)

    !Equations
    ! D^{n+1} = D^n + div Flux --- POINTWISE
    ! So (integrating by parts)
    ! <\gamma, D^{n+1}> = <\gamma, D^n> - <\nabla\gamma, F>
    ! and a consistent theta-method for Q is
    ! <\gamma, D^{n+1}Q^{n+1}> = <\gamma, D^nQ^n> - 
    !                                     \theta<\nabla\gamma,FQ^{n+1}> 
    !                                  -(1-\theta)<\nabla\gamma,FQ^n>


    ! Get J and detwei
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), detwei&
         &=detwei, J=J)
    !Advection terms
    l_adv_mat = t_theta*dshape_dot_vector_shape(&
         Q_shape%dn,Flux_gi,Q_shape,Q_shape%quadrature%weight)
    l_rhs = - (1-t_theta)*dshape_dot_vector_rhs(&
         Q_shape%dn,Flux_gi,Q_gi*Q_shape%quadrature%weight)
    !Mass terms
    l_adv_mat = l_adv_mat + &
         shape_shape(ele_shape(Q,ele),ele_shape(Q,ele),D_gi*detwei)
    l_rhs = l_rhs + shape_rhs(ele_shape(Q,ele),Q_gi*D_old_gi*detwei)

    call addto(Q_rhs,ele_nodes(Q_rhs,ele),l_rhs)
    call addto(adv_mat,ele_nodes(Q_rhs,ele),ele_nodes(Q_rhs,ele),&
         l_adv_mat)    

  end subroutine construct_advection_cg_tracer_ele

  subroutine construct_pv_flux_ele(QFlux,Q,Q_old,D,D_old,&
       Flux,X,down,t_theta,ele)
    type(scalar_field), intent(in) :: Q,Q_old,D,D_old
    type(vector_field), intent(inout) :: QFlux
    type(vector_field), intent(in) :: X, Flux, Down
    integer, intent(in) :: ele
    real, intent(in) :: t_theta
    !
    real, dimension(mesh_dim(X),ele_loc(QFlux,ele)) :: QFlux_perp_rhs
    real, dimension(Flux%dim,ele_ngi(Flux,ele)) :: Flux_gi, Flux_perp_gi
    real, dimension(ele_ngi(X,ele)) :: Q_gi,Q_old_gi,D_gi,D_old_gi
    real, dimension(X%dim, ele_ngi(Q, ele)) :: up_gi
    real, dimension(mesh_dim(QFlux), mesh_dim(QFlux), ele_ngi(QFlux,ele)) &
         :: Metric, Metricf
    real, dimension(mesh_dim(QFlux),mesh_dim(QFlux),ele_loc(QFlux,ele)&
         &,ele_loc(QFlux,ele)) :: l_u_mat
    real, dimension(mesh_dim(QFlux)*ele_loc(QFlux,ele),mesh_dim(QFlux)&
         &*ele_loc(QFlux,ele)) :: solve_mat
    real, dimension(ele_loc(Q,ele)) :: Q_test1_rhs, Q_test2_rhs
    real, dimension(mesh_dim(Qflux)*ele_loc(QFlux,ele)) :: solve_rhs
    real, dimension(mesh_dim(Qflux),ele_ngi(Qflux,ele)) :: &
         contravariant_velocity_gi, velocity_perp_gi
    type(element_type), pointer :: Q_shape, QFlux_shape
    integer :: loc, orientation,gi, dim1, dim2
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(ele_ngi(X,ele)) :: detwei, detJ
    real, dimension(X%dim, X%dim, ele_ngi(X,ele)) :: rot

    Q_shape => ele_shape(Q,ele)
    QFlux_shape => ele_shape(QFlux,ele)
    Q_gi = ele_val_at_quad(Q,ele)
    Q_old_gi = ele_val_at_quad(Q_old,ele)
    D_gi = ele_val_at_quad(D,ele)
    D_old_gi = ele_val_at_quad(D_old,ele)
    Flux_gi = ele_val_at_quad(Flux,ele)

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J, detwei, detJ)

    up_gi = -ele_val_at_quad(down,ele)
    call get_up_gi(X,ele,up_gi,orientation)

    !Equations
    ! <\nabla \gamma, (-v,u)> = <(\gamma_x, \gamma_y),(-v,u)>
    !                         = <(\gamma_y,-\gamma_x),( u,v)>
    !                         = <-\nabla^\perp\gamma,(u,v)>

    ! and a consistent theta-method for Q is
    ! <\gamma, D^{n+1}Q^{n+1}> = <\gamma, D^nQ^n> - 
    !                                     \theta<\nabla\gamma,FQ^{n+1}> 
    !                                  -(1-\theta)<\nabla\gamma,FQ^n>
    ! So velocity update is 
    ! <\gamma, \zeta^{n+1}> = <-\nabla^\perp\gamma,u^{n+1}>
    != <-\nabla^\perp\gamma,u^n> 
    !  +\theta<-\nabla^\perp\gamma,(FQ^{n+1})^\perp>
    !  +(1-\theta)<-\nabla^\perp\gamma,(FQ^n)^\perp>
    ! taking w = -\nabla^\perp\gamma,
    ! <w,u^{n+1}>=<w,u^n> + <w,F(\theta Q^{n+1}+(1-\theta)Q^n)^\perp>
    ! We actually store the inner product with test function (FQ_rhs)
    ! (avoids having to do a solve)

    do gi  = 1, ele_ngi(QFlux,ele)
       Flux_gi(:,gi) = Flux_gi(:,gi)*(t_theta*Q_gi(gi) + &
         (1-t_theta)*Q_old_gi(gi))
    end do

    do gi=1, ele_ngi(QFlux,ele)
       rot(1,:,gi)=(/0.,-up_gi(3,gi),up_gi(2,gi)/)
       rot(2,:,gi)=(/up_gi(3,gi),0.,-up_gi(1,gi)/)
       rot(3,:,gi)=(/-up_gi(2,gi),up_gi(1,gi),0./)
    end do
    do gi=1,ele_ngi(QFlux,ele)
       Metric(:,:,gi)=matmul(J(:,:,gi), transpose(J(:,:,gi)))/detJ(gi)
       Metricf(:,:,gi)=matmul(J(:,:,gi), &
            matmul(rot(:,:,gi), transpose(J(:,:,gi))))/detJ(gi)
    end do

    do gi = 1, ele_ngi(QFlux,ele)
       Flux_perp_gi(:,gi) = matmul(Metricf(:,:,gi),Flux_gi(:,gi))
    end do
    
    QFlux_perp_rhs = shape_vector_rhs(QFlux_shape,Flux_perp_gi,&
         & QFlux_shape%quadrature%weight)
    call set(QFlux,ele_nodes(QFlux,ele),QFlux_perp_rhs)

  end subroutine construct_pv_flux_ele

  subroutine test_pv_flux_ele(Qtest1,Qtest2,QFlux,Q,Q_old,D,D_old,&
       Flux,X,down,t_theta,ele)
    type(scalar_field), intent(in) :: Q,Q_old,D,D_old
    type(scalar_field), intent(inout) :: Qtest1, Qtest2
    type(vector_field), intent(in) :: X, Flux, Down,QFlux
    integer, intent(in) :: ele
    real, intent(in) :: t_theta
    !
    real, dimension(mesh_dim(X),ele_loc(QFlux,ele)) :: QFlux_perp_rhs
    real, dimension(ele_ngi(X,ele)) :: Q_gi,Q_old_gi,D_gi,D_old_gi
    real, dimension(Flux%dim,ele_ngi(Flux,ele)) :: Flux_gi, Flux_perp_gi
    real, dimension(X%dim, ele_ngi(Q, ele)) :: up_gi
    real, dimension(mesh_dim(QFlux), mesh_dim(QFlux), ele_ngi(QFlux,ele)) &
         :: Metric, Metricf
    real, dimension(mesh_dim(QFlux),mesh_dim(QFlux),ele_loc(QFlux,ele)&
         &,ele_loc(QFlux,ele)) :: l_u_mat
    real, dimension(mesh_dim(QFlux)*ele_loc(QFlux,ele),mesh_dim(QFlux)&
         &*ele_loc(QFlux,ele)) :: solve_mat
    real, dimension(ele_loc(Q,ele)) :: Q_test1_rhs, Q_test2_rhs
    real, dimension(mesh_dim(Qflux)*ele_loc(QFlux,ele)) :: solve_rhs
    real, dimension(mesh_dim(Qflux),ele_ngi(Qflux,ele)) :: &
         contravariant_velocity_gi, velocity_perp_gi
    type(element_type), pointer :: Q_shape, QFlux_shape
    integer :: loc, orientation,gi, dim1, dim2
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(ele_ngi(X,ele)) :: detwei, detJ
    real, dimension(X%dim, X%dim, ele_ngi(X,ele)) :: rot

    Q_shape => ele_shape(Q,ele)
    QFlux_shape => ele_shape(QFlux,ele)
    Q_gi = ele_val_at_quad(Q,ele)
    Q_old_gi = ele_val_at_quad(Q_old,ele)
    D_gi = ele_val_at_quad(D,ele)
    D_old_gi = ele_val_at_quad(D_old,ele)
    Flux_gi = ele_val_at_quad(Flux,ele)
    Qflux_perp_rhs = ele_val(Qflux,ele)

    do gi  = 1, ele_ngi(QFlux,ele)
       Flux_gi(:,gi) = Flux_gi(:,gi)*(t_theta*Q_gi(gi) + &
         (1-t_theta)*Q_old_gi(gi))
    end do

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J, detwei, detJ)

    up_gi = -ele_val_at_quad(down,ele)
    call get_up_gi(X,ele,up_gi,orientation)

    do gi=1, ele_ngi(QFlux,ele)
       rot(1,:,gi)=(/0.,-up_gi(3,gi),up_gi(2,gi)/)
       rot(2,:,gi)=(/up_gi(3,gi),0.,-up_gi(1,gi)/)
       rot(3,:,gi)=(/-up_gi(2,gi),up_gi(1,gi),0./)
    end do
    do gi=1,ele_ngi(QFlux,ele)
       Metric(:,:,gi)=matmul(J(:,:,gi), transpose(J(:,:,gi)))/detJ(gi)
       Metricf(:,:,gi)=matmul(J(:,:,gi), &
            matmul(rot(:,:,gi), transpose(J(:,:,gi))))/detJ(gi)
    end do

    l_u_mat = shape_shape_tensor(qflux_shape, qflux_shape, &
         qflux_shape%quadrature%weight,Metric)
    solve_mat = 0.
    solve_rhs = 0.
    do dim1 = 1, mesh_dim(qflux)
       do dim2 = 1, mesh_dim(qflux)
          solve_mat( (dim1-1)*ele_loc(qflux,ele)+1:dim1*ele_loc(qflux,ele),&
               (dim2-1)*ele_loc(qflux,ele)+1:dim2*ele_loc(qflux,ele)) = &
               l_u_mat(dim1,dim2,:,:)
       end do
       solve_rhs( (dim1-1)*ele_loc(qflux,ele)+1:dim1*ele_loc(qflux,ele))&
            = QFlux_perp_rhs(dim1,:)
    end do
    call solve(solve_mat,solve_rhs)
    !Don't need to include constraints since u eqn is
    ! <w,Q^\perp> + <<[w],\lambda>> + \sum_i C_i\cdot w \Gamma_i = RHS
    ! but if w=-\nabla^\perp\gamma then [w]=0 and C_i\cdot w=0.

    do dim1 = 1, mesh_dim(qflux)
       flux_perp_gi(dim1,:) = matmul(transpose(QFlux_shape%n),&
            solve_rhs((dim1-1)*ele_loc(qflux,ele)+1:dim1*ele_loc(qflux&
            &,ele)))
    end do
    !Now compute vorticity
    do gi=1,ele_ngi(X,ele)
       Metric(:,:,gi)=matmul(J(:,:,gi), transpose(J(:,:,gi)))/detJ(gi)
       contravariant_velocity_gi(:,gi) = &
            matmul(Metric(:,:,gi),Flux_perp_gi(:,gi))
    end do

    select case(mesh_dim(X))
    case (2)
       velocity_perp_gi(1,:) = -contravariant_velocity_gi(2,:)
       velocity_perp_gi(2,:) =  contravariant_velocity_gi(1,:)
       Q_test1_rhs = dshape_dot_vector_rhs(Q%mesh%shape%dn, &
            velocity_perp_gi,Q%mesh%shape%quadrature%weight)
       Q_test1_rhs = q_test1_rhs*orientation
    case default
       FLAbort('Exterior derivative not implemented for given mesh dimension' )
    end select

    Q_test2_rhs = shape_rhs(ele_shape(Q,ele),(Q_gi*D_gi-Q_old_gi&
         &*D_old_gi)*detwei)

    call addto(Qtest1,ele_nodes(Q,ele),Q_test1_rhs)
    call addto(Qtest2,ele_nodes(Q,ele),Q_test2_rhs)

  end subroutine test_pv_flux_ele

  subroutine get_PV(state,PV,velocity,D,path)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(in) :: velocity
    type(scalar_field), intent(inout) :: PV
    type(scalar_field), intent(in) :: D
    character(len=OPTION_PATH_LEN), optional :: path
    !
    type(vector_field), pointer :: X, down
    type(scalar_field), pointer :: Coriolis
    type(csr_matrix) :: pv_mass_matrix
    type(csr_sparsity), pointer :: pv_mass_sparsity, curl_sparsity
    type(scalar_field) :: pv_rhs
    integer :: ele

    ewrite(1,*) 'subroutine get_pv'

    pv_mass_sparsity => &
         get_csr_sparsity_firstorder(state, pv%mesh, pv%mesh)
    call allocate(pv_mass_matrix,pv_mass_sparsity)
    call zero(pv_mass_matrix)

    X=>extract_vector_field(state, "Coordinate")
    down=>extract_vector_field(state, "GravityDirection")
    Coriolis=>extract_scalar_field(state, "Coriolis")

    call allocate(pv_rhs, pv%mesh, 'PvRHS')
    call zero(pv_rhs)
    
    do ele = 1, ele_count(pv)
       call assemble_pv_ele(pv_mass_matrix,pv_rhs,velocity,D,Coriolis,X&
            &,down,ele)
    end do

    if(present(path)) then
       call petsc_solve(pv,pv_mass_matrix,pv_rhs,option_path=path)
    else
       call petsc_solve(pv,pv_mass_matrix,pv_rhs)
    end if

    call deallocate(pv_mass_matrix)
    call deallocate(pv_rhs)

  end subroutine get_pv

  subroutine assemble_PV_ele(pv_mass_matrix,pv_rhs,velocity,D,Coriolis,X&
       &,down,ele)
    !!Subroutine to compute the pv from a *local velocity* field
    type(vector_field), intent(in) :: X, down, velocity
    type(scalar_field), intent(inout) :: pv_rhs
    type(scalar_field), intent(in) :: D, Coriolis
    type(csr_matrix), intent(inout) :: pv_mass_matrix
    integer, intent(in) :: ele
    !
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(ele_ngi(X,ele)) :: detwei, detJ, l_f, l_d
    real, dimension(ele_loc(pv_rhs,ele),ele_loc(pv_rhs,ele))&
         :: l_mass_mat
    real, dimension(ele_loc(pv_rhs,ele))&
         :: l_rhs
    real, dimension(mesh_dim(velocity),ele_ngi(velocity,ele)) :: velocity_gi,&
         contravariant_velocity_gi, velocity_perp_gi
    real, dimension(X%dim, ele_ngi(X,ele)) :: up_gi
    integer :: orientation, gi
    real, dimension(mesh_dim(X), mesh_dim(X), ele_ngi(X,ele)) :: Metric
    !

    up_gi = -ele_val_at_quad(down,ele)
    call get_up_gi(X,ele,up_gi,orientation)

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J, detwei, detJ)    
    l_d = ele_val_at_quad(D,ele)
    l_mass_mat = shape_shape(ele_shape(pv_rhs,ele),&
         ele_shape(pv_rhs,ele),detwei*l_d)

    velocity_gi = ele_val_at_quad(velocity,ele)
    do gi=1,ele_ngi(X,ele)
       Metric(:,:,gi)=matmul(J(:,:,gi), transpose(J(:,:,gi)))/detJ(gi)
       contravariant_velocity_gi(:,gi) = &
            matmul(Metric(:,:,gi),velocity_gi(:,gi))
    end do

    !Relative vorticity

    ! pv = \nabla^\perp\cdot u = v_x - u_y
    ! < \gamma, pv > = <\gamma, v_x - u_y> = <-\gamma_x,v>+<\gamma_y,u>
     ! < -\nabla^\perp \gamma, J^TJ u/detJ> in local coordinates
    ! < \nabla \gamma, (J^TJ u)^\perp/detJ> in local coordinates
    ! requires us to know the orientation of the manifold

    select case(mesh_dim(X))
    case (2)
       velocity_perp_gi(1,:) = -contravariant_velocity_gi(2,:)
       velocity_perp_gi(2,:) =  contravariant_velocity_gi(1,:)
       l_rhs = dshape_dot_vector_rhs(pv_rhs%mesh%shape%dn, &
            velocity_perp_gi,X%mesh%shape%quadrature%weight)
       l_rhs = l_rhs*orientation
    case default
       FLAbort('Exterior derivative not implemented for given mesh dimension')
    end select

    ! Coriolis term
    l_f = ele_val_at_quad(Coriolis,ele)
    l_rhs = l_rhs + shape_rhs(ele_shape(pv_rhs,ele),l_f*detwei)

    call addto(pv_mass_matrix,ele_nodes(pv_rhs,ele),&
         ele_nodes(pv_rhs,ele),l_mass_mat)
    call addto(pv_rhs,ele_nodes(pv_rhs,ele),l_rhs)

  end subroutine assemble_pv_ele

  subroutine compute_U_residual(UResidual,oldU,oldD,newU,newD,PVFlux,state)
    !!< Compute the residual in the U equation, given PVFlux computed
    !!< from the PV advection equation. Note that the PVFlux has 
    !!< already been perped, and multiplied by test functions.
    !!< Residual will be multiplied by test functions
    type(vector_field), intent(inout) :: UResidual
    type(vector_field), intent(in) :: oldU,newU,PVFlux
    type(scalar_field), intent(in) :: oldD,newD
    type(state_type), intent(inout) :: state
    !
    integer :: ele
    type(vector_field), pointer :: X
    real :: dt, theta,g 

    ewrite(1,*) '  subroutine compute_U_residual('

    X=>extract_vector_field(state, "Coordinate")
    call get_option("/timestepping/timestep", dt)
    call get_option("/timestepping/theta",theta)
    call get_option("/physical_parameters/gravity/magnitude", g)

    call zero(UResidual)
    do ele = 1, ele_count(UResidual)
       call compute_U_residual_ele(UResidual,oldU,oldD,newU,newD,PVFlux,X,theta&
            &,dt,g,ele)
    end do
    
    ewrite(1,*) 'END subroutine compute_U_residual('

  end subroutine compute_U_residual

  subroutine compute_U_residual_ele(UResidual,oldU,oldD,newU,newD,PVFlux,X,theta&
       &,dt,g,ele)
    type(vector_field), intent(inout) :: UResidual
    type(vector_field), intent(in) :: oldU,newU,PVFlux,X
    type(scalar_field), intent(in) :: oldD,newD
    real, intent(in) :: dt,theta,g
    integer, intent(in) :: ele
    !
    real, dimension(ele_ngi(X,ele)) :: detwei, detJ
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(mesh_dim(oldU), mesh_dim(oldU)) :: Metric
    real, dimension(mesh_dim(oldU), ele_ngi(X,ele)) :: U_rhs, newU_rhs
    real, dimension(ele_ngi(X,ele)) :: D_gi, newD_gi, D_bar_gi, K_bar_gi
    real, dimension(mesh_dim(X), ele_ngi(X,ele)) :: U_gi,newU_gi
    real, dimension(X%dim, ele_ngi(X,ele)) :: U_cart_gi,newU_cart_gi
    real, dimension(mesh_dim(oldU),ele_loc(UResidual,ele)) :: UR_rhs
    integer :: gi
    type(element_type), pointer :: U_shape

    !r[w] = <w,u^{n+1}-u^n> + <w,FQ^\perp> 
    !           - dt*<div w, g\bar{h} + \bar{|u|^2/2}>

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), &
         detwei=detwei, J=J, detJ=detJ)
    U_shape=>ele_shape(oldU, ele)

    !Get all the variables at the quadrature points
    D_gi = ele_val_at_quad(oldD,ele)
    newD_gi = ele_val_at_quad(newD,ele)
    U_gi = ele_val_at_quad(oldU,ele)
    newU_gi = ele_val_at_quad(newU,ele)
    U_rhs = 0.
    newU_rhs = 0.
    do gi=1,ele_ngi(oldU,ele)
       Metric=matmul(J(:,:,gi), transpose(J(:,:,gi)))/detJ(gi)
       U_rhs(:,gi) = matmul(Metric,U_gi(:,gi))
       newU_rhs(:,gi) = matmul(Metric,newU_gi(:,gi))
    end do
    U_cart_gi = 0.
    newU_cart_gi = 0.
    do gi = 1, ele_ngi(oldU,ele)
       U_cart_gi(:,gi) = matmul(transpose(J(:,:,gi)),U_gi(:,gi))/detJ(gi)
       newU_cart_gi(:,gi) = matmul(transpose(J(:,:,gi)),newU_gi(:,gi))/detJ(gi)
    end do

    !First the PV Flux (perped)
    UR_rhs = -ele_val(PVFlux,ele)
    !Now the time derivative term
    UR_rhs = UR_rhs + shape_vector_rhs(U_shape,newU_rhs-U_rhs,&
         U_shape%quadrature%weight)
    !Now the gradient terms (done in local coordinates)
    D_bar_gi = theta*newD_gi + (1-theta)*D_gi
    K_bar_gi = 0.5*(theta*sum(newU_cart_gi,1)+(1-theta)*sum(U_cart_gi,1))

    UR_rhs = UR_rhs - dt * dshape_rhs(U_shape%dn,&
         (g*D_bar_gi + K_bar_gi)*U_shape%quadrature%weight)

    call set(UResidual,ele_nodes(UResidual,ele),UR_rhs)

  end subroutine compute_U_residual_ele

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
    type(vector_field), pointer :: U, U_old, U_cartesian

    !! Coordinate and advecting velocity fields
    type(vector_field), pointer :: X, U_nl, down

    !! Temporary velocity fields
    type(vector_field) :: U_tmp, U_cartesian_tmp

    !! Velocity components Cartesian coordinates for slope limiter
    type(scalar_field) :: U_component

    !! Change in U over one timestep.
    type(vector_field) :: delta_U, delta_U_tmp

    !! DG Courant number field
    type(scalar_field), pointer :: s_field

    !! Sparsity of advection matrix.    
    type(csr_sparsity), pointer :: sparsity
    
    !! System matrices.
    type(block_csr_matrix) :: A, L, mass_local, inv_mass_local, &
         inv_mass_cartesian

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
    integer :: i, j, dim, ele, dim1
    integer :: upwinding_option

    if(have_option('/material_phase::Fluid/vector_field::Velocity/prognostic&
         &/spatial_discretisation/discontinuous_galerkin/advection_scheme/ed&
         &ge_coordinates_upwind')) then
       upwinding_option = VECTOR_UPWIND_EDGE
    else
       if(have_option('/material_phase::Fluid/vector_field::Velocity/prognostic&
            &/spatial_discretisation/discontinuous_galerkin/advection_scheme/sphere_coordinates_upwind')) then
          upwinding_option = VECTOR_UPWIND_SPHERE
       else
          FLAbort('Unknown upwinding option')
       end if
    end if

    U=>extract_vector_field(state, field_name)
    U_old=>extract_vector_field(state, "Old"//field_name)
    X=>extract_vector_field(state, "Coordinate")
    U_cartesian=>extract_vector_field(state, "Velocity")
    down=>extract_vector_field(state, "GravityDirection")

    dim=mesh_dim(U)

    ! Reset U to value at the beginning of the timestep.
    call set(U, U_old)

    sparsity => get_csr_sparsity_firstorder(state, U%mesh, U%mesh)

    ! Add data space to the sparsity pattern.
    call allocate(A, sparsity, (/dim,X%dim/))
    call zero(A)
    call allocate(L, sparsity, (/dim,X%dim/))
    call zero(L)

    mass_sparsity=make_sparsity_dg_mass(U%mesh)
    call allocate(mass_local, mass_sparsity, (/dim,dim/))
    call zero(mass_local)
    call allocate(inv_mass_local, mass_sparsity, (/dim,dim/))
    call zero(inv_mass_local)
    call allocate(inv_mass_cartesian, mass_sparsity, (/X%dim,X%dim/))
    call zero(inv_mass_cartesian)

    ! Ensure delta_U inherits options from U.
    call allocate(delta_U, U%dim, U%mesh, "delta_U")
    call zero(delta_U)
    call allocate(delta_U_tmp, U%dim, U%mesh, "delta_U")
    call zero(delta_U_tmp)
    delta_U%option_path = U_cartesian%option_path
    call allocate(rhs, U%dim, U%mesh, trim(field_name)//" RHS")
    call zero(rhs)

    ! allocate temp velocity fields
    call allocate(U_tmp, U%dim, U%mesh, 'tmpU')
    call zero(U_tmp)
    call allocate(U_cartesian_tmp, U_cartesian%dim, U_cartesian%mesh, 'tmpUcartesian')
    call zero(U_cartesian_tmp)
   
    call construct_vector_advection_dg(A, L, mass_local, &
         inv_mass_local, inv_mass_cartesian,&
         rhs, field_name, state, upwinding_option, down, &
         velocity_name=velocity_name)
    
    ! Note that since dt is module global, these lines have to
    ! come after construct_advection_diffusion_dg.
    call get_option("/timestepping/timestep", dt)
    
    if(have_option(trim(U_cartesian%option_path)//&
         &"/prognostic/temporal_discretisation/discontinuous_galerkin/&
         &/number_advection_subcycles")) then
       call get_option(trim(U_cartesian%option_path)//&
            &"/prognostic/temporal_discretisation/discontinuous_galerkin/&
            &/number_advection_subcycles", subcycles)
    else
       call get_option(trim(U_cartesian%option_path)//&
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
    if (have_option(trim(U_cartesian%option_path)//&
         "/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/slope_limiter")) then
       limit_slope=.true.
       
       ! Note unsafe for mixed element meshes
       if (element_degree(U,1)==0) then
          FLExit("Slope limiters make no sense for degree 0 fields")
       end if

    end if

    do i=1, subcycles
       call mult_t(U_cartesian_tmp, L, U)
       call mult(U_cartesian, inv_mass_cartesian, U_cartesian_tmp)
        ! dU = Advection * U
       ! A maps from cartesian to local
       call mult(delta_U_tmp, A, U_cartesian)
       ! dU = dU + RHS
       !------------------------------------------
       !Is there anything in RHS?
       call addto(delta_U_tmp, RHS, -1.0)
       !------------------------------------------
       ! dU = M^(-1) dU
       call mult(delta_U, inv_mass_local, delta_U_tmp)
       !------------------------------------------
       !Do we need the below?
       call mult_t(U_cartesian_tmp, L, delta_U)
       call mult(U_cartesian, inv_mass_cartesian, U_cartesian_tmp)
       !------------------------------------------

       ! U = U + dt/s * dU
       call addto(U, delta_U, scale=-dt/subcycles)
       call halo_update(U)

       if (limit_slope) then

          call mult_t(U_cartesian_tmp, L, U)
          call mult(U_cartesian, inv_mass_cartesian, U_cartesian_tmp)

          !Really crap limiter
          !Just limit each cartesian component
          do dim1 = 1, U_cartesian%dim
             u_component = extract_scalar_field(U_cartesian, dim1)
             call limit_vb(state, U_component)
          end do
          !limiter=VECTOR_LIMITER_VB
          !call limit_slope_dg(U_cartesian, state, limiter)

          call mult(U_tmp, L, U_cartesian)
          call mult(U, inv_mass_local, U_tmp)
       end if

    end do

    call deallocate(delta_U)
    call deallocate(A)
    call deallocate(L)
    call deallocate(mass_local)
    call deallocate(inv_mass_local)
    call deallocate(inv_mass_cartesian)
    call deallocate(mass_sparsity)
    call deallocate(rhs)

  end subroutine solve_vector_advection_dg_subcycle

  subroutine construct_vector_advection_dg(A, L, mass_local, inv_mass_local,&
       & inv_mass_cartesian, rhs, field_name, &
       & state, upwinding_option, down, velocity_name) 
    !!< Construct the advection equation for discontinuous elements in
    !!< acceleration form.
    !!< 
    !!< If mass is provided then the mass matrix is not added into big_m or
    !!< rhs. It is instead returned as mass. This may be useful for testing
    !!< or for solving equations otherwise than in acceleration form.

    !! Main advection matrix.    
    type(block_csr_matrix), intent(inout) :: A, L
    !! Mass matrices.
    type(block_csr_matrix), intent(inout) :: mass_local, inv_mass_local, &
         inv_mass_cartesian
    !! Right hand side vector.
    type(vector_field), intent(inout) :: rhs
    type(vector_field), intent(in) :: down
    !! Name of the field to be advected.
    character(len=*), intent(in) :: field_name
    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state
    integer, intent(in) :: upwinding_option

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

    ewrite(1,*) "Writing advection equation for "&
         &//trim(field_name)

    ! These names are based on the CGNS SIDS.
    U=extract_vector_field(state, field_name)
    X=extract_vector_field(state, "Coordinate")

    if(present(velocity_name)) then
      lvelocity_name = velocity_name
    else
      lvelocity_name = "NonlinearVelocity"
    end if

    U_nl=extract_vector_field(state, lvelocity_name)
    call incref(U_nl)

    assert(has_faces(X%mesh))
    assert(has_faces(U%mesh))
    
    call zero(A)
    call zero(RHS)
    call zero(mass_local)

    element_loop: do ele=1,element_count(U)
       
       call construct_vector_adv_element_dg(ele, A, L,&
            & mass_local, inv_mass_local, inv_mass_cartesian, rhs,&
            & X, U, U_nl, upwinding_option, down)
       
    end do element_loop
    
    ! Drop any extra field references.

    call deallocate(U_nl)

  end subroutine construct_vector_advection_dg

  subroutine construct_vector_adv_element_dg(ele, A, L,&
       & mass_local, inv_mass_local, inv_mass_cartesian, rhs,&
       & X, U, U_nl, upwinding_option, down)
    !!< Construct the advection_diffusion equation for discontinuous elements in
    !!< acceleration form.
    implicit none
    !! Index of current element
    integer :: ele
    !! Main advection matrix.
    type(block_csr_matrix), intent(inout) :: A
    !! Transformation matrix
    type(block_csr_matrix), intent(inout) :: L
    !! Mass matrices.
    type(block_csr_matrix), intent(inout) :: mass_local, inv_mass_local, &
         inv_mass_cartesian
    !! Right hand side vector.
    type(vector_field), intent(inout) :: rhs
    
    !! Position, velocity and advecting velocity.
    type(vector_field), intent(in) :: X, U, U_nl, down

    !! Upwinding option
    integer, intent(in) :: upwinding_option

    ! Bilinear forms.
    real, dimension(mesh_dim(U),mesh_dim(U),ele_loc(U,ele),ele_loc(U,ele)) :: local_mass_mat
    real, dimension(mesh_dim(U),X%dim,ele_loc(U,ele),ele_loc(U,ele)) :: l_mat
    real, dimension(ele_loc(U,ele),ele_loc(U,ele)) :: cartesian_mass_mat
    real, dimension(mesh_dim(U),X%dim,ele_loc(U,ele),ele_loc(U,ele)) :: &
         Advection_mat

    ! Local assembly matrices.
    real, dimension(mesh_dim(U)*ele_loc(U,ele), mesh_dim(U)*ele_loc(U,ele)) :: l_mass, l_mass_cartesian

    ! Local variables.
    
    ! Neighbour element, face and neighbour face.
    integer :: ele_2, face, face_2
    ! Loops over faces.
    integer :: ni
    
    ! Transform from local to physical coordinates.
    real, dimension(ele_ngi(X,ele)) :: detwei, detJ
    real, dimension(mesh_dim(U), X%dim, ele_ngi(X,ele)) :: J, J_scaled
    real, dimension(X%dim, mesh_dim(U), ele_ngi(X,ele)) :: pinvJ
    real, dimension(ele_loc(U,ele), ele_ngi(X,ele), X%dim) :: dshape
    real, dimension(mesh_dim(U), mesh_dim(U), ele_ngi(X,ele)) :: G

    ! Different velocities at quad points.
    real, dimension(U_nl%dim, ele_ngi(U_nl, ele)) :: U_nl_q
    real, dimension(X%dim, ele_ngi(U_nl, ele)) :: U_cartesian_q
    real, dimension(ele_ngi(U_nl, ele)) :: U_nl_div_q

    ! Node and shape pointers.
    integer, dimension(:), pointer :: U_ele
    type(element_type), pointer :: U_shape, U_nl_shape
    ! Neighbours of this element.
    integer, dimension(:), pointer :: neigh

    integer :: i, gi, dim, dim1, dim2, nloc, k

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
       G(:,:,gi)=matmul(J(:,:,gi),transpose(J(:,:,gi)))/detJ(gi)
       J_scaled(:,:,gi)=J(:,:,gi)/detJ(gi)
    end do

    local_mass_mat=shape_shape_tensor(u_shape, u_shape, &
         u_shape%quadrature%weight, G)

    cartesian_mass_mat=shape_shape(u_shape, u_shape, detwei)

    ! Transformation matrix
    l_mat=shape_shape_tensor(u_shape, u_shape, u_shape%quadrature%weight, J)

    ! Advecting velocity at quadrature points.
    U_nl_q=ele_val_at_quad(U_nl,ele)
    U_nl_div_q=ele_div_at_quad(U_nl, ele, U_shape%dn)

    ! WE HAVE INTEGRATED BY PARTS TWICE
    ! Element advection matrix
    !    / 
    !    | w . div (\bar{U} \tensor u )dV
    !    / 
    
    ! becomes
    
    !    /
    !    | \phi_i (dx/d\xi)^T.\div_L (\bar{U}_L \tensor \phi_j)dV_L
    !    /
    ! underscore U indicates local coordinates

    ! split into two terms

    !    /
    !    | \phi_i (dx/d\xi)^T(\div_L \bar{U}_L)\phi_j dV_L
    !    /

    ! and 

    !    /
    !    | \phi_i (dx/d\xi)^T\bar{U}_L . grad_L \phi_j dV_L
    !    /


    U_cartesian_q=0.
    do gi=1,ele_ngi(X,ele)
       U_cartesian_q(:,gi)=matmul(transpose(J(:,:,gi)),U_nl_q(:,gi))
       U_cartesian_q(:,gi)=U_cartesian_q(:,gi)/detJ(gi)
    end do
    Advection_mat = shape_vector_dot_dshape_tensor(U_shape, U_nl_q, U_shape&
         &%dn, J_scaled, U_shape%quadrature%weight)

   !----------------------------------------------------------------------
   ! Perform global assembly.
   !----------------------------------------------------------------------

   ! Assemble matrices.

   ! Return local mass separately.
   do dim1 = 1,dim
      do dim2 = 1,dim
         call addto(mass_local, dim1, dim2, U_ele, U_ele, &
              local_mass_mat(dim1,dim2,:,:))
      end do
   end do

   ! calculate inverse_mass_local
   nloc=ele_loc(U,ele)
   do dim1 = 1,dim
      do dim2 = 1,dim
         l_mass(nloc*(dim1-1)+1:nloc*dim1, nloc*(dim2-1)+1:nloc*dim2)=&
              local_mass_mat(dim1,dim2,:,:)
      end do
   end do

   call invert(l_mass)

   do dim1 = 1,dim
      do dim2 = 1,dim
         call addto(inv_mass_local, dim1, dim2, U_ele, U_ele, &
              l_mass(nloc*(dim1-1)+1:nloc*dim1,nloc*(dim2-1)+1:nloc*dim2))
      end do
   end do

   ! calculate inverse_mass_cartesian
   call invert(cartesian_mass_mat)

   do dim1 = 1,X%dim
         call addto(inv_mass_cartesian, dim1, dim1, U_ele, U_ele, cartesian_mass_mat)
   end do

   ! advection matrix
   do dim1 = 1,dim
      do dim2 = 1,X%dim
         call addto(A, dim1, dim2, U_ele, U_ele, Advection_mat(dim1,dim2,:,:))
      end do
   end do

   ! transformation matrix
   do dim1 = 1,dim
      do dim2 = 1,X%dim
         call addto(L, dim1, dim2, U_ele, U_ele, l_mat(dim1,dim2,:,:))
      end do
   end do

   !-------------------------------------------------------------------
   ! Interface integrals
   !-------------------------------------------------------------------
    
   neigh=>ele_neigh(U, ele)

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
      end if

      call construct_vector_adv_interface_dg(ele, ele_2, face, face_2,&
           & A, rhs, X, U, U_nl, upwinding_option, down)

   end do neighbourloop
    
 end subroutine construct_vector_adv_element_dg
  
  subroutine construct_vector_adv_interface_dg(ele, ele_2, face, face_2, &
       & A, rhs, X, U, U_nl, upwinding_option, down)

    !!< Construct the DG element boundary integrals on the ni-th face of
    !!< element ele.
    implicit none

    integer, intent(in) :: ele, ele_2, face, face_2
    type(block_csr_matrix), intent(inout) :: A
    type(vector_field), intent(inout) :: rhs
    ! We pass these additional fields to save on state lookups.
    type(vector_field), intent(in) :: X, U, U_nl, down
    integer, intent(in) :: upwinding_option

    ! Face objects and numberings.
    type(element_type), pointer :: U_shape, U_shape_2
    integer, dimension(face_loc(U,face)) :: U_face, U_face_l
    integer, dimension(face_loc(U,face_2)) :: U_face_2
    integer, dimension(face_loc(U_nl,face)) :: U_nl_face
    integer, dimension(face_loc(U_nl,face_2)) :: U_nl_face_2

    ! Note that both sides of the face can be assumed to have the same
    ! number of quadrature points.
    real, dimension(U_nl%dim, face_ngi(U_nl, face)) :: n1_l, n2_l, U_nl_q,&
         & U_nl_f_q, U_nl_f2_q
    logical, dimension(face_ngi(U_nl, face)) :: inflow
    real, dimension(face_ngi(U_nl, face)) :: U_nl_q_dotn1_l, U_nl_q_dotn2_l, U_nl_q_dotn_l, income
    real, dimension(X%dim, face_ngi(X,face)) :: n1, n2, t11, t12, t21,t22,&
         & pole_axis
    real, dimension(X%dim, ele_ngi(X,ele)) :: up_gi, up_gi2
    real, dimension(X%dim, face_ngI(X,face)) :: X_quad
    real, dimension(X%dim, X%dim, face_ngi(X,face)) :: Btmp
    real, dimension(mesh_dim(U), X%dim, face_ngi(X,face)) :: B
    ! Variable transform times quadrature weights.
    real, dimension(ele_ngi(X,ele)) :: detwei, detJ
    real, dimension(face_ngi(X,face)) :: detwei_f, detJ_f
    real, dimension(mesh_dim(U), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(mesh_dim(U), X%dim, face_ngi(X,face)) :: J_scaled
    real, dimension(face_ngi(X,face)) :: inner_advection_integral, outer_advection_integral

    ! Bilinear forms
    real, dimension(mesh_dim(U),X%dim,face_loc(U,face),face_loc(U,face)) :: nnAdvection_out
    real, dimension(mesh_dim(U),X%dim,face_loc(U,face),face_loc(U,face_2)) :: nnAdvection_in

    ! normal weights
    real :: w1, w2
    integer :: dim, dim1, dim2, gi, i, k

    dim=mesh_dim(U)

    U_face=face_global_nodes(U, face)
    U_shape=>face_shape(U, face)
    X_quad = face_val_at_quad(X, face)

    U_face_2=face_global_nodes(U, face_2)
    U_shape_2=>face_shape(U, face_2)
    
    !Unambiguously calculate the normal using the face with the higher
    !face number. This is so that the normal is identical on both sides.

    call get_local_normal(n1_l, w1, U, local_face_number(U%mesh,face))
    call get_local_normal(n2_l, w2, U, local_face_number(U%mesh,face_2))
    
    !----------------------------------------------------------------------
    ! Construct element-wise quantities.
    !----------------------------------------------------------------------

    ! Advecting velocity at quadrature points.
    U_nl_f_q = face_val_at_quad(U_nl, face)
    U_nl_f2_q = face_val_at_quad(U_nl, face_2)

    u_nl_q_dotn1_l = sum(U_nl_f_q*w1*n1_l,1)
    u_nl_q_dotn2_l = -sum(U_nl_f2_q*w2*n2_l,1)
    U_nl_q_dotn_l=0.5*(u_nl_q_dotn1_l+u_nl_q_dotn2_l)

    ! Inflow is true if the flow at this gauss point is directed
    ! into this element.
    inflow = u_nl_q_dotn_l<0.0
    income = merge(1.0,0.0,inflow)

    ! Calculate tensor on face
    if(ele_loc(X,ele).ne.3) then
       ewrite(1,*) 'Nloc=',ele_loc(X,ele)
       FLAbort('Hard coded for linear elements at the moment')
    end if
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), detwei=detwei, J=J, detJ=detJ)
    forall(gi=1:face_ngi(X,face))
       J_scaled(:,:,gi)=J(:,:,1)/detJ(1)
    end forall

    select case(upwinding_option)
    case (VECTOR_UPWIND_EDGE)
       ! Get outward pointing normals in physical space for face1
       ! inward pointing normals in physical space for face2
       n1=get_face_normal_manifold(X, ele, face)
       !needs checking
       if(ele_2<0) then
          ! external boundary
          n2=n1
       else
          n2=-get_face_normal_manifold(X, ele_2, face_2)
       end if
       ! Form 'bending' tensor
       ! This rotates the normal to face 2 into -normal from face 1
       ! u_b = t(t.u_b) + n_b(n_b.u_b)
       ! u_a = t(t.u_b) + n_a(n_b.u_b)
       !     = u_b - n_b(n_b.u_b) + n_a(n_b.u_b)
       !     = (I + (n_a-n_b)n_b^T)u_b == B u_b
       forall(gi=1:face_ngi(X,face))
          Btmp(:,:,gi)=outer_product(n1(:,gi)-n2(:,gi),n2(:,gi))
          forall(i=1:size(Btmp,1)) Btmp(i,i,gi)=Btmp(i,i,gi)+1.0
          B(:,:,gi)=matmul(J_scaled(:,:,gi),Btmp(:,:,gi))
       end forall
    case (VECTOR_UPWIND_SPHERE)
       if(mesh_dim(X).ne.2) then
          FLExit('Assumes 2d surface. Surface of the sphere, in fact.')
       end if
       !Get element normal
       up_gi = -ele_val_at_quad(down,ele)
       call get_up_gi(X,ele,up_gi)
       up_gi2 = -ele_val_at_quad(down,ele_2)
       call get_up_gi(X,ele_2,up_gi2)
       ! Form 'bending' tensor
       ! Coordinate system is:
       ! t1: normal to local normal and (0,0,1)
       ! t2: normal to t1 and local normal
       ! u2 = t21(t21.u2) + t22(t22.u2)
       ! u1 = t11(t21.u2) + t12(t22.u2)
       !    = (t11 t21^T + t12 t22^T)u2

       pole_axis = 0.
       pole_axis(3,:) = 1.
       btmp = 0.
       do gi = 1, face_ngi(X,face)
          t11(:,gi) = cross_product(pole_axis(:,gi),up_gi(:,1))
          t11(:,gi) = t11(:,gi)/sqrt(sum(t11(:,gi)**2))
          t12(:,gi) = cross_product(up_gi(:,1),t11(:,gi))
          t21(:,gi) = cross_product(pole_axis(:,gi),up_gi2(:,1))
          t21(:,gi) = t21(:,gi)/sqrt(sum(t21(:,gi)**2))
          t22(:,gi) = cross_product(up_gi2(:,1),t21(:,gi))
          
          forall(i=1:3,k=1:3)
             Btmp(i,k,gi) = t11(i,gi)*t21(k,gi) + t12(i,gi)*t22(k,gi)
          end forall
          B(:,:,gi)=matmul(J_scaled(:,:,gi),Btmp(:,:,gi))
       end do
    case default
       FLAbort('Unknown vector upwinding option')
    end select

    !----------------------------------------------------------------------
    ! Construct bilinear forms.
    !----------------------------------------------------------------------
       
    ! Calculate outflow boundary integral.
    ! can anyone think of a way of optimising this more to avoid
    ! superfluous operations (i.e. multiplying things by 0 or 1)?

    ! first the integral around the inside of the element
    ! (this is the flux *out* of the element)
    inner_advection_integral = -income*u_nl_q_dotn_l
    nnAdvection_out=shape_shape_tensor(U_shape, U_shape,  &
         &inner_advection_integral*U_shape%quadrature%weight, J_scaled)
    
    ! now the integral around the outside of the element
    ! (this is the flux *in* to the element)
    outer_advection_integral = income*u_nl_q_dotn_l
    nnAdvection_in=shape_shape_tensor(U_shape, U_shape_2, &
         &outer_advection_integral*U_shape%quadrature%weight, B)
       
    !----------------------------------------------------------------------
    ! Perform global assembly.
    !----------------------------------------------------------------------

    ! Insert advection in matrix.

    ! Outflow boundary integral.
    do dim1 = 1, dim
       do dim2 = 1, X%dim
          call addto(A, dim1, dim2, U_face, U_face, nnAdvection_out(dim1,dim2,:,:))
       end do
    end do

    ! Inflow boundary integral.
    do dim1 = 1, dim
       do dim2 = 1, X%dim
          call addto(A, dim1, dim2, U_face, U_face_2, nnAdvection_in(dim1,dim2,:,:))
        end do
    end do

    contains

      ! copy of outer_product in vector tools, but now as function
      pure function outer_product(x, y) result (M)
        !!< Give two column vectors x, y
        !!< compute the matrix xy*

        real, dimension(:), intent(in) :: x, y
        real, dimension(size(x), size(y)) :: M
        integer :: i, j

        forall (i=1:size(x))
           forall (j=1:size(y))
              M(i, j) = x(i) * y(j)
           end forall
        end forall
      
      end function outer_product
    
  end subroutine construct_vector_adv_interface_dg

end module advection_local_DG
