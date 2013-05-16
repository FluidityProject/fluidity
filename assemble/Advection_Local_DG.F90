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
  use shallow_water_diagnostics

  implicit none

  ! Buffer for output messages.
  character(len=255), private :: message

  private
  public solve_advection_dg_subcycle, &
       & solve_advection_cg_tracer, compute_U_residual, get_pv

  ! Local private control parameters. These are module-global parameters
  ! because it would be expensive and/or inconvenient to re-evaluate them
  ! on a per-element or per-face basis
  real :: dt

  ! Whether to include various terms
  logical :: include_advection

  ! parameters for SSPRK schemes
  real, dimension(1), target :: alpha1 = (/ 0.0 /)
  real, dimension(2), target :: alpha2 = (/ 0.0, 0.5 /)
  real, dimension(3), target :: alpha3 = (/ 0.0, 0.75, 0.3333333333333333 /)

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
    type(scalar_field) :: delta_T, delta_T_total, Tstage

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
    type(vector_field) :: tmpflux
    !! Mesh where UpwindFlux lives
    type(mesh_type), pointer :: lambda_mesh
    !! Upwind Flux matrix
    type(csr_matrix) :: UpwindFluxMatrix

    character(len=FIELD_NAME_LEN) :: limiter_name
    integer :: i, k, ele, nstages
    real :: meancheck, ldt
    real, dimension(:), pointer :: alpha

    ewrite(1,*) 'subroutine solve_advection_dg_subcycle'


    T=>extract_scalar_field(state, field_name)
    call allocate(Tstage,T%mesh,"Stage"//field_name)
    T_old=>extract_scalar_field(state, "Old"//field_name)
    X=>extract_vector_field(state, "Coordinate")
    U_nl=>extract_vector_field(state, velocity_name)

    ewrite(2,*) 'UVALS in dg', maxval(abs(U_nl%val))

    ! Reset T to value at the beginning of the timestep.
    call set(T, T_old)
    if(present(Flux)) then
       call zero(Flux)
       call allocate(tmpflux,flux%dim,flux%mesh, "TmpFlux")
       call zero(tmpflux)
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
    !/material_phase::Fluid/scalar_field::LayerThickness/prognostic/spatial_discretisation/discontinuous_galerkin/slope_limiter::Vertex_Based
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

    !! Tmpflux comes from the initial limiting step and gets added on
    !! at the end.
    if (limit_slope) then
       ! Filter wiggles from T
       call limit_vb_manifold(state,T,delta_T)
       if(present(flux)) then
          call zero(upwindflux)
          call mult(delta_T_total, mass, delta_T)
          call update_flux(flux, Delta_T_total, UpwindFlux, Mass)
       end if
   end if

   call get_option(trim(T%option_path)//'/prognostic/temporal_discretisation&
        &/discontinuous_galerkin/SSPRK_order',nstages)
   select case(nstages)
   case (1)
      alpha => alpha1
   case (2)
      alpha => alpha2
   case (3)
      alpha => alpha3
   case default
      FLAbort('Number of stages not supported.')
   end select

   ldt = dt/subcycles
   subcycle_loop: do i=1, subcycles
      ewrite(1,*) 'SUBCYCLE', i

      if(present(flux)) then
         call set(tmpflux,flux)
      end if

      call set(Tstage,T)

      RKstage_loop: do k = 1, nstages

         ! dT = Advection * T
         call mult(delta_T, matrix, Tstage)
         ! dT = dT + RHS
         call addto(delta_T, RHS, -1.0)
         if(present(Flux)) then
            if(maxval(abs(RHS%val))>1.0e-8) then
               FLAbort('Flux reconstruction doesn''t work if diffusion/bcs pres&
                    &ent at the moment')
            end if
         end if
         call scale(delta_T, -ldt)
         ! dT = M^(-1) dT
         if(present(flux)) then
            call mult(UpwindFlux, upwindfluxmatrix, Tstage)
            call scale(UpwindFlux,-ldt)
            call update_flux(Flux, Delta_T, UpwindFlux, Mass)
            call scale(Flux, 1-alpha(k))
            call addto(Flux, tmpflux, alpha(k))
         end if
         ! T = T + dt/s * dT
         call dg_apply_mass(inv_mass, delta_T)
         ewrite(2,*) 'Delta_T', maxval(abs(delta_T%val))
         call addto(Tstage, delta_T)
         call scale(Tstage, 1-alpha(k))
         call addto(Tstage,T,alpha(k))
         !Probably need to halo_update(Flux) as well       
         call halo_update(Tstage)
         
         if (limit_slope) then
            ! Filter wiggles from T
            call limit_vb_manifold(state,Tstage,delta_T)
            if(present(flux)) then
               call mult(delta_t_total, mass, delta_t)
               call zero(UpwindFlux)
               call update_flux(Flux, Delta_T_total, UpwindFlux, Mass)
            end if
         end if
      end do RKstage_loop

      call set(T,Tstage)
      if(present(flux)) then
         !! Add on the initial flux due to the first slope limiter
         call addto(Flux,tmpflux)
      end if

   end do subcycle_loop

   if(present(Flux)) then
      !! Debugging checks
      if(have_option('/material_phase::Fluid/scalar_field::LayerThickness/p&
           &rognostic/spatial_discretisation/debug')) then
         do ele = 1, ele_count(Flux)
            call check_flux(Flux,T,T_old,X,mass,ele)
         end do
      end if
      call deallocate(UpwindFlux)
      call deallocate(UpwindFluxMatrix)
      call deallocate(delta_T_total)
      call deallocate(tmpflux)
   end if
   
   call deallocate(tstage)
   call deallocate(delta_T)
   call deallocate(matrix)
   call deallocate(mass)
   call deallocate(inv_mass)
   call deallocate(mass_sparsity)
   call deallocate(rhs)
   
 end subroutine solve_advection_dg_subcycle

  subroutine check_flux(Flux,T,T_old,X,mass,ele)
    type(vector_field), intent(in) :: Flux, X
    type(scalar_field), intent(in) :: T,T_old
    integer, intent(in) :: ele
    type(csr_matrix), intent(in) :: mass
    !
    real, dimension(Flux%dim,ele_loc(Flux,ele)) :: Flux_vals
    real, dimension(ele_ngi(T,ele)) :: Div_Flux_gi
    real, dimension(ele_ngi(T,ele)) :: Delta_T_gi, T_gi
    real, dimension(ele_loc(T,ele)) :: Delta_T_rhs, Div_Flux_rhs,T_rhs
    integer :: dim1, loc
    real, dimension(ele_ngi(X,ele)) :: detwei
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J_mat
    real :: residual
    !

    !! Check that T-T_{old} = div(F)

    Delta_T_gi = ele_val_at_quad(T,ele)-ele_val_at_quad(T_old,ele)
    T_gi = ele_val_at_quad(T,ele)
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
    T_rhs = shape_rhs(ele_shape(T,ele),T_gi*detwei)
    Div_Flux_rhs = shape_rhs(ele_shape(T,ele),Div_Flux_gi*&
         Flux%mesh%shape%quadrature%weight)
    
    residual = &
         maxval(abs(Delta_T_rhs-Div_Flux_rhs))/max(1.0&
         &,maxval(abs(T_rhs)))
    if(residual>1.0e-6) then
       ewrite(2,*) residual, maxval(abs(T_rhs))
       FLExit('Flux residual error.')
    end if

  end subroutine check_flux

  subroutine update_flux(Flux, Delta_T, UpwindFlux, Mass)
    !! Subroutine to compute div-conforming Flux such that
    !! div Flux = Delta_T
    !! with Flux.n = UpwindFlux on element boundaries
    !! and then to add to Flux
    type(vector_field), intent(inout) :: Flux
    type(scalar_field), intent(in) :: Delta_T
    type(scalar_field), intent(in) :: UpwindFlux
    type(csr_matrix), intent(in) :: Mass
    !
    integer :: ele, i
    real :: Area
    integer, dimension(:), pointer :: T_ele

    !Checking the Flux is in local representation
    assert(Flux%dim==mesh_dim(Flux))

    do ele = 1, ele_count(Flux)
       T_ele => ele_nodes(delta_T,ele)
       area = 0.
       do i = 1, size(T_ele)
          area = area + &
               sum(row_val(Mass,t_ele(i)))
       end do
       call update_flux_ele(Flux, Delta_T, UpwindFlux,ele,area)
    end do

  end subroutine update_flux

  subroutine update_flux_ele(Flux, Delta_T, UpwindFlux,ele,area)
    !! Subroutine to compute div-conforming Flux such that
    !! div Flux = Delta_T
    !! with Flux.n = UpwindFlux on element boundaries
    type(vector_field), intent(inout) :: Flux
    type(scalar_field), intent(in) :: Delta_T
    type(scalar_field), intent(in) :: UpwindFlux
    integer, intent(in) :: ele
    real, intent(in) :: area
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
         &/max(1.0,maxval(ele_val(Delta_T,ele)))/area
    if(residual>1.0e-10) then
       ewrite(0,*) 'Residual = ', residual
       FLExit('Bad residual')
    end if

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

    if(flux_constraint%n_grad_basis>0) then
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
    end if

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
    print*, row, size(flux_mat,1), size(flux_mat,2)
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

    delta_t_val = ele_val(delta_t,ele)
    div_flux_val = shape_rhs(T_shape,div_flux_gi*T_shape%quadrature%weight)

    residual = maxval(abs(div_flux_val-delta_T_val))/max(1.0&
         &,maxval(abs(delta_T_val)))/area
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

    !Get normals
    call get_local_normal(n1, w1, U_nl, local_face_number(U_nl%mesh,face))
    call get_local_normal(n2, w2, U_nl, local_face_number(U_nl%mesh,face_2))    

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
       
       !Factor of 0.5 is because we visit from both sides
       if(face.ne.face_2) then
          call addto(upwindfluxmatrix, flux_face, T_face_2,&
               0.5*sign*upwindflux_mat_in)
       end if
       call addto(upwindfluxmatrix, flux_face, T_face,&
            0.5*sign*upwindflux_mat_out)

      end subroutine construct_upwindflux_interface

  end subroutine construct_adv_interface_dg  

  subroutine solve_advection_cg_tracer(Q,D,D_old,Flux,U_nl,QF,state)
    !!< Solve the continuity equation for D*Q with Flux
    !!< d/dt (D*Q) + div(Flux*Q) = diffusion terms.
    !!< Done on the vorticity mesh
    !!< Return a PV flux Q defined on the velocity mesh
    !!< which is the projection of Flux*Q into the velocity space
    !!< Note that Flux and Q contain a factor of dt for convenience.
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(in),target :: D,D_old
    type(scalar_field), intent(inout),target :: Q
    type(vector_field), intent(in) :: Flux,U_nl
    type(vector_field), intent(inout) :: QF
    !
    type(csr_sparsity), pointer :: Q_sparsity
    type(csr_matrix) :: Adv_mat
    type(scalar_field) :: Q_rhs, Qtest1,Qtest2
    type(vector_field), pointer :: X, down
    type(scalar_field), pointer :: Q_old
    integer :: ele
    real :: dt, t_theta, residual, disc_max, c1
    type(scalar_field), pointer :: discontinuity_detector_field
    character(len = OPTION_PATH_LEN) :: discontinuity_detector_name
    type(scalar_field), dimension(:), pointer :: Qstages

    ewrite(1,*) '  subroutine solve_advection_cg_tracer('

    X=>extract_vector_field(state, "Coordinate")
    down=>extract_vector_field(state, "GravityDirection")
    Q_old=>extract_scalar_field(state, "Old"//trim(Q%name))
    call get_option('/timestepping/timestep',dt)
    call get_option(trim(Q%option_path)//'/prognostic/timestepping/theta',t_theta)
    
    !set up matrix and rhs
    Q_sparsity => get_csr_sparsity_firstorder(state, Q%mesh, Q%mesh)
    call allocate(adv_mat, Q_sparsity) ! Add data space to the sparsity
    ! pattern.
    call zero(adv_mat)
    call allocate(Q_rhs,Q%mesh,trim(Q%name)//"RHS")
    call zero(Q_rhs)
    call zero(QF)

    if(have_option(trim(Q%option_path)//'/prognostic/spatial_discretisation/&
         &continuous_galerkin/discontinuity_capturing')) then
       call get_option(trim(Q%option_path)//'/prognostic/spatial_discretisati&
            &on/continuous_galerkin/discontinuity_capturing/discontinuity_in&
            &dicator_name',discontinuity_detector_name)

       discontinuity_detector_field => extract_scalar_field(&
            &state,trim(discontinuity_detector_name))
       disc_max = maxval(discontinuity_detector_field%val)
    end if

    if(have_option(trim(Q%option_path)//'/prognostic/spatial_discretisation/&
         &continuous_galerkin/theta_method')) then
       !! Use theta method to integrate q
       do ele = 1, ele_count(Q)
          call construct_advection_cg_tracer_theta_ele(Q_rhs,adv_mat,Q,D,D_old,&
               &Discontinuity_detector_field,Flux&
               &,X,down,U_nl,dt,t_theta,disc_max,ele)
       end do
       ewrite(2,*) 'Q_RHS', maxval(abs(Q_rhs%val))
       call petsc_solve(Q,adv_mat,Q_rhs)
       !! Compute the PV flux to pass to velocity equation
       do ele = 1, ele_count(Q)
          call construct_pv_flux_theta_ele(QF,Q,Q_old,D,D_old,&
               &Discontinuity_detector_field,Flux,&
               &X,down,U_nl,t_theta,disc_max,ele)
       end do
    else if (have_option(trim(Q%option_path)//'/prognostic/spatial_discretis&
         &ation/continuous_galerkin/taylor_galerkin')) then

       call get_option(trim(Q%option_path)//'/prognostic/spatial_discretis&
            &ation/continuous_galerkin/taylor_galerkin/eta',eta)
       
       !Set up stage coefficients
       if (have_option(trim(Q%option_path)//'/prognostic/spatial_discretis&
         &ation/continuous_galerkin/taylor_galerkin/Tbar2_3_scheme')) then
          !! Use the \bar{T}(2,3) scheme with positive sign
          n_stages = 2
          
          c1 = 0.5*(1 + (-1./3.+8*eta)**0.5)
          mcoeffs(1,:) = (/ c1, 0 /)
          mcoeffs(2,:) = (/ 0.5*(3-1./c1), 0.5*(1./c1-1) /)
          ncoeffs(1,:) = (/ 0.5*c1**2-eta, 0.0 /)
          ncoeffs(2,:) = (/ 0.25*(3*c1-1)-eta, 0.25*(1-c1) /)
       else
          FLAbort('Unknown choice of TG scheme')
       end if
       call construct_taylor_galerkin_stage_ele(&
            & Q_rhs,adv_mat,Q_stages,D,D_old,Flux,X,&
            & dt,ele)
       call petsc_solve(Q_stages
       FLAbort('this is where TG option goes')
    end if

    if(have_option('/material_phase::Fluid/scalar_field::PotentialVorticity/&
         &prognostic/debug')) then
       call allocate(Qtest1,Q%mesh,trim(Q%name)//"test1")
       call zero(Qtest1)
       call allocate(Qtest2,Q%mesh,trim(Q%name)//"test2")
       call zero(Qtest2)
       do ele = 1, ele_count(Q)
          call test_pv_flux_ele(Qtest1,Qtest2,QF,Q,Q_old,D,D_old,&
               Flux,X,t_theta,ele)
       end do
       ewrite(2,*) 'Error = ', maxval(abs(Qtest2%val)), maxval(abs(Qtest1%val))
       ewrite(2,*) 'Error = ', maxval(abs(Qtest1%val-Qtest2%val)), maxval(Qtest1%val)
       residual = maxval(abs(Qtest1%val-Qtest2%val))/max(1.0&
            &,maxval(abs(Qtest1%val)))
       ewrite(2,*) 'residual from qtest', residual
       assert(residual<1.0e-6)
       ewrite(2,*) 'test passed'
       call deallocate(Qtest1)
       call deallocate(Qtest2)
    end if

    call deallocate(adv_mat)
    call deallocate(Q_rhs)

    if(have_option(trim(Q%option_path)//'/prognostic/spatial_discretisation/&
         &continuous_galerkin/discontinuity_capturing')) then
       call calculate_discontinuity_detector(state,&
            discontinuity_detector_field)
    end if

    ewrite(1,*) 'END  subroutine solve_advection_cg_tracer('

  end subroutine solve_advection_cg_tracer
  
  subroutine construct_advection_cg_tracer_theta_ele(Q_rhs,adv_mat,Q,D,D_old,&
       Discontinuity_detector_field,Flux,&
       & X,down,U_nl,dt,t_theta,disc_max,ele)
    type(scalar_field), intent(in) :: D,D_old,Q, Discontinuity_detector_field
    type(scalar_field), intent(inout) :: Q_rhs
    type(csr_matrix), intent(inout) :: Adv_mat
    type(vector_field), intent(in) :: X, Flux, down,U_nl
    integer, intent(in) :: ele
    real, intent(in) :: dt, t_theta,disc_max
    !
    real, dimension(ele_loc(Q,ele),ele_loc(Q,ele)) :: l_adv_mat,tmp_mat
    real, dimension(ele_loc(Q,ele)) :: l_rhs
    real, dimension(U_nl%dim,ele_ngi(U_nl,ele)) :: U_nl_gi
    real, dimension(X%dim, ele_ngi(U_nl,ele)) :: U_cart_gi
    real, dimension(Flux%dim,ele_ngi(Flux,ele)) :: Flux_gi
    real, dimension(Flux%dim,FLux%dim,ele_ngi(Flux,ele)) :: Grad_Flux_gi,&
         & tensor_gi
    real, dimension(ele_ngi(X,ele)) :: detwei, Q_gi, D_gi, &
         & D_old_gi,div_flux_gi, detwei_l, detJ, DI_gi
    real, dimension(ele_loc(D,ele)) :: D_val
    real, dimension(ele_loc(Q,ele)) :: Q_val
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    type(element_type), pointer :: Q_shape, Flux_shape
    integer :: loc,dim1,dim2,gi
    real :: tau, alpha, tol, area, h,c_sc,ratio
    real, dimension(mesh_dim(U_nl), mesh_dim(U_nl), ele_ngi(Q,ele)) &
         :: MetricT
    real, dimension(mesh_dim(X),ele_ngi(X,ele)) :: gradQ
    type(element_type), pointer :: D_shape
    !

    ! Get J and detwei
    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), detwei&
         &=detwei, J=J, detJ=detJ)

    D_shape => ele_shape(D,ele)
    Q_shape => ele_shape(Q,ele)
    Flux_shape => ele_shape(Flux,ele)
    D_val = invert_pi_ele(ele_val(D,ele),D_shape,detwei)
    D_gi = matmul(transpose(D_shape%n),D_val)
    D_val = invert_pi_ele(ele_val(D_old,ele),D_shape,detwei)
    D_old_gi = matmul(transpose(D_shape%n),D_val)
    Q_gi = ele_val_at_quad(Q,ele)
    Flux_gi = ele_val_at_quad(Flux,ele)
    U_nl_gi = ele_val_at_quad(U_nl,ele)
    Q_val = ele_val(Q,ele)

    !Equations
    ! D^{n+1} = D^n + \pi div Flux 
    ! If define \int_K\phi\hat{D}/J J dS = \int_K\phi D J dS, then
    ! <\phi,\pi(\hat{D}/J)> = <\phi, D> for all \phi,
    ! => \pi \hat{D} = D
    ! \hat{D}^{n+1}/J = \hat{D}^n/J + div Flux

    ! So (integrating by parts)
    ! <\gamma, \hat{D}^{n+1}/J> = 
    !             <\gamma, D^n/J> - <\nabla\gamma, F>

    ! and a consistent theta-method for Q is
    ! <\gamma,\hat{D}^{n+1}Q^{n+1}/J> = 
    !                  <\gamma, \hat{D}^nQ^n/J> - 
    !                                     \theta<\nabla\gamma,FQ^{n+1}> 
    !                                  -(1-\theta)<\nabla\gamma,FQ^n>

    detwei_l = Q_shape%quadrature%weight
    !Advection terms
    !! <-grad gamma, FQ >
    !! Can be evaluated locally using pullback
    tmp_mat = dshape_dot_vector_shape(&
         Q_shape%dn,Flux_gi,Q_shape,detwei_l)
    l_adv_mat = t_theta*tmp_mat
    l_rhs = -(1-t_theta)*matmul(tmp_mat,Q_val)

    !Mass terms
    !! <gamma, QD>
    !! Requires transformed integral
    tmp_mat = shape_shape(ele_shape(Q,ele),ele_shape(Q,ele),D_gi*detwei_l)
    l_adv_mat = l_adv_mat + tmp_mat
    tmp_mat = shape_shape(ele_shape(Q,ele),ele_shape(Q,ele),D_old_gi*detwei_l)
    l_rhs = l_rhs + matmul(tmp_mat,Q_val)

    !Laplacian filter options
    if(have_option(trim(Q%option_path)//'/prognostic/spatial_discretisation/&
         &continuous_galerkin/theta_method/laplacian_filter')) then
       FLExit('Laplacian filter not coded yet.')
    end if
    !Streamline upwinding options
    if(have_option(trim(Q%option_path)//'/prognostic/spatial_discretisation/&
         &continuous_galerkin/theta_method/streamline_upwinding')) then
       !! Replace test function \gamma with \gamma + \tau u\cdot\nabla\gamma
       call get_option(trim(Q%option_path)//'/prognostic/spatial_discretisat&
            &ion/continuous_galerkin/theta_method/streamline_upwinding/alpha',alpha)
       call get_option(trim(Q%option_path)//'/prognostic/spatial_discretisat&
            &ion/continuous_galerkin/theta_method/streamline_upwinding/tol',tol)

       !Gradients of fluxes
       Grad_flux_gi = ele_grad_at_quad(Flux,ele,Flux_shape%dn)
       div_flux_gi = 0.
       tensor_gi = 0.
       do dim1 = 1, Flux%dim
          div_flux_gi = div_flux_gi + grad_flux_gi(dim1,dim1,:)
          do dim2 = 1, Flux%dim
             tensor_gi(dim1,dim2,:) = U_nl_gi(dim1,:)*Flux_gi(dim2,:)
          end do
       end do
       !real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
       do gi = 1, ele_ngi(X,ele)
          U_cart_gi(:,gi) = matmul(transpose(J(:,:,gi)),U_nl_gi(:,gi))/detJ(gi)
       end do

       if(maxval(sqrt(sum(U_cart_gi**2,1)))<tol) then
          tau = 0.
       else
          area = sum(detwei)
          h = sqrt(4*area/sqrt(3.))
          tau = h*alpha/maxval(sqrt(sum(U_cart_gi**2,1)))
       end if
       !! Mass terms
       !! tau < u . grad gamma, QD>
       !! = tau < grad gamma, u QD>
       !! can do locally because of pullback formula
       tmp_mat = dshape_dot_vector_shape(&
            Q_shape%dn,U_nl_gi,Q_shape,D_gi*detwei_l/detJ)
       l_adv_mat = l_adv_mat + tau*tmp_mat
       tmp_mat = dshape_dot_vector_shape(&
            Q_shape%dn,U_nl_gi,Q_shape,D_old_gi*detwei_l/detJ)
       l_rhs = l_rhs + tau*matmul(tmp_mat,Q_val)
       !! Advection terms
       !! tau < u. grad gamma, div (F Q)>
       !! = tau <grad gamma, u (F.grad Q + Q div F)>
       !! can do locally because of pullback formula
       !! tensor = u outer F (F already contains factor of dt)
       tmp_mat = dshape_dot_vector_shape(&
            Q_shape%dn,U_nl_gi,Q_shape,div_flux_gi*detwei_l/detJ) &
            + dshape_tensor_dshape(&
            Q_shape%dn,tensor_gi,Q_shape%dn,detwei_l/detJ)
       l_adv_mat = l_adv_mat - tau*t_theta*tmp_mat
       l_rhs = l_rhs + tau*(1-t_theta)*matmul(tmp_mat,Q_val)
    end if
    !Discontinuity capturing options
    if(have_option(trim(Q%option_path)//'/prognostic/spatial_discretisation/&
         &continuous_galerkin/discontinuity_capturing')) then

       call get_option(trim(Q%option_path)//'/prognostic/&
            &spatial_discretisation/&
            &continuous_galerkin/discontinuity_capturing/&
            &scaling_coefficient',c_sc)
       call get_option(trim(Q%option_path)//'/prognostic/&
            &spatial_discretisation/&
            &continuous_galerkin/discontinuity_capturing/&
            &filter_ratio',ratio)

       !need to work out diffusion term in physical coordinates

       ! grad perp psi = J\hat{grad}^\perp\hat{\psi}/det J

       ! grad perp gamma . grad perp psi det J
       ! = (-gamma_2)^T(M_11 M_12)(-psi_2)
       ! = ( gamma_1)  (M_21 M_22)( psi_1)
       !
       ! = (gamma_2)^T( M_11 -M_12)( psi_2)
       ! = (gamma_1)  (-M_21  M_22)( psi_1)

       ! = (gamma_1)^T( M_22 -M_21)( psi_1)
       ! = (gamma_2)  (-M_12  M_11)( psi_2)
       
       do gi=1,ele_ngi(Q,ele)
          MetricT(:,:,gi)=matmul(J(:,:,gi), transpose(J(:,:,gi)))/detJ(gi)**2
       end do
       MetricT(1,1,:) = MetricT(2,2,:)
       MetricT(2,2,:) = MetricT(1,1,:)
       MetricT(1,2,:) =-MetricT(2,1,:)
       MetricT(2,1,:) =-MetricT(1,2,:)

       area = sum(detwei)
       h = sqrt(4*area/sqrt(3.))
       !D_gi = 0.5*(ele_val_at_quad(D_old,ele) + &
       !     ele_val_at_quad(D,ele))
       DI_gi = ele_val_at_quad(discontinuity_detector_field,ele)
       if(DI_gi(1)<ratio*disc_max) then
          DI_gi = 0.
       end if
       tmp_mat = dshape_tensor_dshape(Q_shape%dn, &
            MetricT, Q_shape%dn,&
            h*D_gi*DI_gi*Q_shape%quadrature%weight)
       l_adv_mat = l_adv_mat + dt*t_theta*tmp_mat
       l_rhs = l_rhs - (1-t_theta)*dt*matmul(tmp_mat,Q_val)
    end if

    call addto(Q_rhs,ele_nodes(Q_rhs,ele),l_rhs)
    call addto(adv_mat,ele_nodes(Q_rhs,ele),ele_nodes(Q_rhs,ele),&
         l_adv_mat)    

  end subroutine construct_advection_cg_tracer_theta_ele

  subroutine construct_pv_flux_theta_ele(QFlux,Q,Q_old,D,D_old,&
       Discontinuity_detector_field,Flux,X,down,U_nl,t_theta,disc_max,ele)
    type(scalar_field), intent(in) :: Q,Q_old,D,D_old,&
         & Discontinuity_detector_field
    type(vector_field), intent(inout) :: QFlux
    type(vector_field), intent(in) :: X, Flux, Down, U_nl
    integer, intent(in) :: ele
    real, intent(in) :: t_theta,disc_max
    !
    real, dimension(mesh_dim(X),ele_loc(QFlux,ele)) :: QFlux_perp_rhs
    real, dimension(Flux%dim,ele_ngi(Flux,ele)) :: Flux_gi, Flux_perp_gi,&
         QFlux_gi
    real, dimension(ele_ngi(X,ele)) :: Q_gi,Q_old_gi,D_gi,D_old_gi,DI_gi,Qbar_gi
    real, dimension(ele_loc(D,ele)) :: D_val
    real, dimension(X%dim, ele_ngi(Q, ele)) :: up_gi
    real, dimension(mesh_dim(QFlux),mesh_dim(QFlux),ele_loc(QFlux,ele)&
         &,ele_loc(QFlux,ele)) :: l_u_mat
    real, dimension(mesh_dim(QFlux)*ele_loc(QFlux,ele),mesh_dim(QFlux)&
         &*ele_loc(QFlux,ele)) :: solve_mat
    real, dimension(ele_loc(Q,ele)) :: Q_test1_rhs, Q_test2_rhs
    real, dimension(mesh_dim(Qflux)*ele_loc(QFlux,ele)) :: solve_rhs
    real, dimension(mesh_dim(Qflux),ele_ngi(Qflux,ele)) :: &
         contravariant_velocity_gi, velocity_perp_gi
    type(element_type), pointer :: Q_shape, QFlux_shape, Flux_shape
    integer :: loc, orientation,gi, dim1, dim2
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(ele_ngi(X,ele)) :: detwei, detJ, detwei_l, div_flux_gi
    real, dimension(mesh_dim(Q),ele_ngi(Q,ele)) :: gradQ_gi
    real, dimension(Flux%dim,FLux%dim,ele_ngi(Flux,ele)) :: Grad_Flux_gi
    real, dimension(U_nl%dim,ele_ngi(U_nl,ele)) :: U_nl_gi
    real, dimension(X%dim,ele_ngi(U_nl,ele)) :: U_cart_gi
    real, dimension(mesh_dim(U_nl), mesh_dim(U_nl), ele_ngi(Q,ele)) &
         :: MetricT
    real :: residual, alpha, tol, tau, area, h, c_sc, ratio
    type(element_type), pointer :: D_shape
    !

    !! How this works is the following. We want to find F in the DG space
    !! such that 
    !! 
    !! d/dt<-grad gamma, u> = <grad gamma, F>
    !! = \sum_E \int_E grad gamma.F dx 
    !! = \sum_E \int_{\hat{E}} \hat{grad}\hat{gamma}\cdot\hat{F}dx
    !! = \sum_E \sum_g (\hat{grad}\hat{\gamma})_g\cdot\hat{F}_g w_g
    !! for all test functions \gamma.
    !! All of these are exact identities by the pullback formula
    
    !! So, if we have a flux term of the form
    !! = \sum_E \sum_g (\hat{grad}\hat{\gamma})_g\cdot N_g w_g
    !! where N_g are arbitrary flux values defined at quadrature points
    !! then we can define \hat{F}_g via
    !! \sum_E  \sum_g v_g\cdot\hat{F}_g w_g = \sum_E\sum_g v_g\cdot N_g w_g
    !! For all dg test functions v (and therefore v=\hat{\grad}\hat{\gamma}
    !! as a special case).

    D_shape => ele_shape(D,ele)

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J, detwei, detJ)

    Q_shape => ele_shape(Q,ele)
    Flux_shape => ele_shape(Flux,ele)
    QFlux_shape => ele_shape(QFlux,ele)
    Q_gi = ele_val_at_quad(Q,ele)
    Q_old_gi = ele_val_at_quad(Q_old,ele)
    Qbar_gi = t_theta*Q_gi + (1-t_theta)*Q_old_gi
    D_val = invert_pi_ele(ele_val(D,ele),D_shape,detwei)
    D_gi = matmul(transpose(D_shape%n),D_val)
    D_val = invert_pi_ele(ele_val(D_old,ele),D_shape,detwei)
    D_old_gi = matmul(transpose(D_shape%n),D_val)
    Flux_gi = ele_val_at_quad(Flux,ele)
    U_nl_gi = ele_val_at_quad(U_nl,ele)

    detwei_l = Q_shape%quadrature%weight
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
    ! So vorticity update is 
    ! <\gamma, \zeta^{n+1}> = <-\nabla^\perp\gamma,u^{n+1}>
    != <-\nabla^\perp\gamma,u^n> 
    !  +\theta<-\nabla^\perp\gamma,(FQ^{n+1})^\perp>
    !  +(1-\theta)<-\nabla^\perp\gamma,(FQ^n)^\perp>
    ! taking w = -\nabla^\perp\gamma,
    ! <w,u^{n+1}>=<w,u^n> + <w,F(\theta Q^{n+1}+(1-\theta)Q^n)^\perp>
    ! <w,u^{n+1}>=<w,u^n> + <w,\bar{Q}^\perp>
    ! We actually store the inner product with test function (FQ_rhs)
    ! (avoids having to do a solve)
    
    do gi  = 1, ele_ngi(QFlux,ele)
       QFlux_gi(:,gi) = Flux_gi(:,gi)*Qbar_gi(gi)
    end do

    !------------------------------------------------

    !! Additional dissipative terms.
    !! They all take the form
    !! <\nabla gamma, G> 
    !! Which can be evaluated in local coordinates due to pullback magic

    !Laplacian filter options
    if(have_option(trim(Q%option_path)//'/prognostic/spatial_discretisation/&
         &continuous_galerkin/theta_method/laplacian_filter')) then
       FLExit('Laplacian filter not coded yet.')
    end if
    !Streamline upwinding options
    if(have_option(trim(Q%option_path)//'/prognostic/spatial_discretisation/&
         &continuous_galerkin/theta_method/streamline_upwinding')) then
       !! Replace test function \gamma with \gamma + \tau u\cdot\nabla\gamma
       call get_option(trim(Q%option_path)//'/prognostic/spatial_discretisat&
            &ion/continuous_galerkin/theta_method/streamline_upwinding/alpha',alpha)
       call get_option(trim(Q%option_path)//'/prognostic/spatial_discretisat&
            &ion/continuous_galerkin/theta_method/streamline_upwinding/tol',tol)

       !Gradients of fluxes
       Grad_flux_gi = ele_grad_at_quad(Flux,ele,Flux_shape%dn)
       div_flux_gi = 0.
       do dim1 = 1, Flux%dim
          div_flux_gi = div_flux_gi + grad_flux_gi(dim1,dim1,:)
       end do
       !real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
       do gi = 1, ele_ngi(X,ele)
          U_cart_gi(:,gi) = matmul(transpose(J(:,:,gi)),U_nl_gi(:,gi))/detJ(gi)
       end do

       if(maxval(sqrt(sum(U_cart_gi**2,1)))<tol) then
          tau = 0.
       else
          area = sum(detwei)
          h = sqrt(4*area/sqrt(3.))
          tau = h*alpha/maxval(sqrt(sum(U_cart_gi**2,1)))
       end if
       !! Mass terms
       !! tau < u . grad gamma, -(QD)^{n+1}+(QD)^n>
       !! = tau < grad gamma, u(-(QD)^{n+1}+(QD)^n)>
       !! can do locally because of pullback formula
       !! Extra flux is tau u\Delta(QD)
       !! minus sign from definition of vorticity
       !! tau
       do gi = 1, ele_ngi(Q,ele)
          QFLux_gi(:,gi) = QFlux_gi(:,gi) + &
               & tau*U_nl_gi(:,gi)*(Q_gi(gi)*D_gi(gi)-Q_old_gi(gi)&
               &*D_old_gi(gi))/detJ(gi)

       end do
       !! Advection terms
       !! tau < u. grad gamma, div (F\bar{Q})>
       !! can do locally because of pullback formula
       !! = <grad gamma, tau u (F.grad\bar{Q} + \bar{Q} div F)/det J>_{Ehat}
       !! minus sign because of definition of vorticity
       !! Contains factor of dt
       
       GradQ_gi = t_theta*ele_grad_at_quad(Q,ele,Q_shape%dn) + &
            (1-t_theta)*ele_grad_at_quad(Q_old,ele,Q_shape%dn)

       do gi = 1, ele_ngi(Q,ele)
          QFlux_gi(:,gi) = QFlux_gi(:,gi) - tau*u_nl_gi(:,gi)*(&
               & sum(Flux_gi(:,gi)*gradQ_gi(:,gi)) + &
               & div_flux_gi(gi)*Qbar_gi(gi))/detJ(gi)
       end do
    end if

    !------------------------------------------------

    !Discontinuity capturing options
    if(have_option(trim(Q%option_path)//'/prognostic/spatial_discretisation/&
         &continuous_galerkin/discontinuity_capturing')) then

       call get_option(trim(Q%option_path)//'/prognostic/&
            &spatial_discretisation/&
            &continuous_galerkin/discontinuity_capturing/&
            &scaling_coefficient',c_sc)
       call get_option(trim(Q%option_path)//'/prognostic/&
            &spatial_discretisation/&
            &continuous_galerkin/discontinuity_capturing/&
            &filter_ratio',ratio)
              
       do gi=1,ele_ngi(Q,ele)
          MetricT(:,:,gi)=matmul(J(:,:,gi), transpose(J(:,:,gi)))/detJ(gi)**2
       end do
       MetricT(1,1,:) = MetricT(2,2,:)
       MetricT(2,2,:) = MetricT(1,1,:)
       MetricT(1,2,:) =-MetricT(2,1,:)
       MetricT(2,1,:) =-MetricT(1,2,:)

       area = sum(detwei)
       h = sqrt(4*area/sqrt(3.))
       DI_gi = ele_val_at_quad(discontinuity_detector_field,ele)

       if(DI_gi(1)<ratio*disc_max) then
          DI_gi = 0.
       end if
       GradQ_gi = t_theta*ele_grad_at_quad(Q,ele,Q_shape%dn) + &
            (1-t_theta)*ele_grad_at_quad(Q_old,ele,Q_shape%dn)

       !! Some notes to sort out the sign:
       !! Vorticity update is
       !! <gamma, zeta^{n+1}> = <gamma, zeta^n> 
       !!                        - dt*<grad gamma,eta D grad q>
       !! so
       !! <-grad^\perp gamma,Delta u> = dt*<-grad^\perp gamma, eta D grad^
       !! \perp q>
       !! <w,Delta u> = dt*<w,eta D grad^perp q>
       
       !tmp_mat = dshape_tensor_dshape(Q_shape%dn, &
       !     MetricT, Q_shape%dn,&
       !     h*D_gi*DI_gi*Q_shape%quadrature%weight)

       do gi = 1, ele_ngi(Q,ele)
          QFlux_gi(:,gi) = QFlux_gi(:,gi) + dt*h*D_gi(gi)*DI_gi(gi)* &
               matmul(MetricT(:,:,gi),GradQ_gi(:,gi))
       end do       
       
    end if

    !------------------------------------------------

    !! Putting it all together:

    !! Evaluate 
    !! < w, F^\perp > in local coordinates

    Flux_perp_gi(1,:) = -orientation*QFlux_gi(2,:)
    Flux_perp_gi(2,:) =  orientation*QFlux_gi(1,:)

    QFlux_perp_rhs = shape_vector_rhs(QFlux_shape,Flux_perp_gi,&
         & QFlux_shape%quadrature%weight)

    call set(QFlux,ele_nodes(QFlux,ele),QFlux_perp_rhs)

  end subroutine construct_pv_flux_theta_ele

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
    type(element_type), pointer :: Q_shape, QFlux_shape,D_shape
    integer :: loc, orientation,gi, dim1, dim2
    real, dimension(mesh_dim(X), X%dim, ele_ngi(X,ele)) :: J
    real, dimension(ele_ngi(X,ele)) :: detwei, detJ
    real, dimension(X%dim, X%dim, ele_ngi(X,ele)) :: rot
    real, dimension(ele_loc(D,ele)) :: D_val

    D_shape => ele_shape(D,ele)
    Q_shape => ele_shape(Q,ele)
    QFlux_shape => ele_shape(QFlux,ele)
    Q_gi = ele_val_at_quad(Q,ele)
    Q_old_gi = ele_val_at_quad(Q_old,ele)

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J, detwei, detJ)
    D_val = invert_pi_ele(ele_val(D,ele),D_shape,detwei)
    D_gi = matmul(transpose(D_shape%n),D_val)
    D_val = invert_pi_ele(ele_val(D_old,ele),D_shape,detwei)
    D_old_gi = matmul(transpose(D_shape%n),D_val)

    Flux_gi = ele_val_at_quad(Flux,ele)
    Qflux_perp_rhs = ele_val(Qflux,ele)

    do gi  = 1, ele_ngi(QFlux,ele)
       Flux_gi(:,gi) = Flux_gi(:,gi)*(t_theta*Q_gi(gi) + &
         (1-t_theta)*Q_old_gi(gi))
    end do

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

    !!Q_test1_rhs contains -<\nabla\gamma,Q>

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

    !!Q_test2_rhs contains <\gamma,Q^{n+1}\tilde{D}^{n+1}-
    !!                             Q^n\tilde{D}^n>

    Q_test2_rhs = shape_rhs(ele_shape(Q,ele),(Q_gi*D_gi-Q_old_gi&
         &*D_old_gi)*Q_shape%quadrature%weight)

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
    !!DEbugging
    integer ::dim1
    real, dimension(3,3) :: X_val

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
    real, dimension(ele_loc(D,ele)) :: D_val
    real, dimension(ele_loc(pv_rhs,ele),ele_loc(pv_rhs,ele))&
         :: l_mass_mat
    real, dimension(ele_loc(pv_rhs,ele))&
         :: l_rhs
    real, dimension(mesh_dim(velocity),ele_ngi(velocity,ele)) :: velocity_gi,&
         contravariant_velocity_gi, velocity_perp_gi
    real, dimension(X%dim, ele_ngi(X,ele)) :: up_gi
    integer :: orientation, gi
    real, dimension(mesh_dim(X), mesh_dim(X), ele_ngi(X,ele)) :: Metric
    type(element_type), pointer :: D_shape
    !

    D_shape => ele_shape(D,ele)

    up_gi = -ele_val_at_quad(down,ele)
    call get_up_gi(X,ele,up_gi,orientation)

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), J, detwei, detJ)    
    D_val = ele_val(D,ele)
    D_val = invert_pi_ele(D_val,D_shape,detwei)
    l_d = matmul(transpose(D_shape%n),D_val)

    l_mass_mat = shape_shape(ele_shape(pv_rhs,ele),&
         ele_shape(pv_rhs,ele),l_d*d_shape%quadrature%weight)

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

  function invert_pi_ele(D_val_in,D_shape,detwei) result (D_val_out)
    real, intent(in), dimension(:) :: D_val_in
    type(element_type), intent(in) :: D_shape
    real, intent(in), dimension(:) :: detwei
    real, dimension(size(D_val_in)) :: D_val_out !The result
    !
    real, dimension(size(D_val_in),size(D_val_in)) :: proj_mat, &
         & mass_mat

    !! \sum_i \phi_i \hat{D}_i w_i = \sum_i \phi_i D_i J_i w_i.

    proj_mat = shape_shape(D_shape,D_shape,D_shape%quadrature%weight)
    mass_mat = shape_shape(D_shape,D_shape,detwei)
    D_val_out = matmul(mass_mat,D_val_in)
    call solve(proj_mat,D_val_out)

  end function invert_pi_ele

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
    integer :: ele, stat
    type(vector_field), pointer :: X
    type(scalar_field), pointer :: Orography
    real :: dt, theta,g 
    type(scalar_field), target :: l_orography
    logical :: have_orography

    ewrite(1,*) '  subroutine compute_U_residual('

    X=>extract_vector_field(state, "Coordinate")
    call get_option("/timestepping/timestep", dt)
    call get_option("/timestepping/theta",theta)
    call get_option("/physical_parameters/gravity/magnitude", g)

    Orography=>extract_scalar_field(state, "Orography", stat)
    have_orography = (stat==0)
    if(.not.have_orography) then
       call allocate(l_orography,X%mesh,"L_Orography")
       orography => l_orography
       call zero(l_orography)
    end if

    call zero(UResidual)
    do ele = 1, ele_count(UResidual)
       call compute_U_residual_ele(UResidual,oldU,oldD,newU,newD,PVFlux,&
            orography,X,theta&
            &,dt,g,ele)
    end do
    
    if(.not.have_orography) then
       call deallocate(l_orography)
    end if

    ewrite(1,*) 'END subroutine compute_U_residual('

  end subroutine compute_U_residual

  subroutine compute_U_residual_ele(UResidual,oldU,oldD,newU,newD,PVFlux,&
       orography,X,theta&
       &,dt,g,ele)
    type(vector_field), intent(inout) :: UResidual
    type(vector_field), intent(in) :: oldU,newU,PVFlux,X
    type(scalar_field), intent(in) :: oldD,newD,orography
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
    real, dimension(ele_ngi(X,ele)) :: orography_gi

    !r[w] = <w,u^{n+1}-u^n> - <w,FQ^\perp> 
    !           - dt*<div w, g\bar{h} + \bar{|u|^2/2}>

    call compute_jacobian(ele_val(X,ele), ele_shape(X,ele), &
         detwei=detwei, J=J, detJ=detJ)
    U_shape=>ele_shape(oldU, ele)

    !Get all the variables at the quadrature points
    D_gi = ele_val_at_quad(oldD,ele)
    orography_gi = ele_val_at_quad(orography,ele)
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

    !!Residual is
    !! U_new - (U_old + dt*Q^\perp - dt grad(g\bar{D} + \bar{K}))
    !First the PV Flux (perped, contains factor of dt)
    UR_rhs = -ele_val(PVFlux,ele)
    !Now the time derivative term
    UR_rhs = UR_rhs + shape_vector_rhs(U_shape,newU_rhs-U_rhs,&
         U_shape%quadrature%weight)
    !Now the gradient terms (done in local coordinates)
    D_bar_gi = theta*newD_gi + (1-theta)*D_gi + orography_gi
    K_bar_gi = 0.5*(theta*sum(newU_cart_gi**2,1)+(1-theta)*sum(U_cart_gi**2,1))
    !integration by parts, so minus sign
    UR_rhs = UR_rhs - dt * dshape_rhs(U_shape%dn,&
         (g*D_bar_gi + K_bar_gi)*U_shape%quadrature%weight)

    call set(UResidual,ele_nodes(UResidual,ele),UR_rhs)

  end subroutine compute_U_residual_ele

end module advection_local_DG
