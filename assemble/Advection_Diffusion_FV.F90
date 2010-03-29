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

module advection_diffusion_fv
  !!< This module contains the finite volume form of the advection
  !!< -diffusion equation for scalars.
  use elements
  use sparse_tools
  use fetools
  use fields
  use state_module
  use shape_functions
  use global_numbering
  use transform_elements
  use vector_tools
  use fldebug
  use vtk_interfaces
  use solvers
  use boundary_conditions
  use boundary_conditions_from_options
  use spud
  use upwind_stabilisation
  use sparsity_patterns, only: make_sparsity, make_sparsity_transpose
  use sparse_matrices_fields
  use field_options

  implicit none

  private
  public solve_advection_diffusion_fv, construct_advection_diffusion_fv

  ! Local private control parameters. These are module-global parameters
  ! because it would be expensive and/or inconvenient to re-evaluate them
  ! on a per-element or per-face basis
  real :: dt, theta

  ! Weight between conservative and non-conservative forms of the advection
  ! equation. 
  ! 1 is for conservative 0 is for non-conservative.
  real :: beta

  ! do we have a diffusivity or not?
  integer :: diff_stat

  ! Discretisation to use for diffusion term.
  integer :: diffusion_scheme
  integer, parameter :: BASSIREBAY=1
  integer, parameter :: LDG=2
  ! Boundary condition types:
  ! (the numbers should match up with the order in the 
  !  get_entire_boundary_condition call)
  integer :: BCTYPE_WEAKDIRICHLET=1

  real :: penalty
  logical :: penalise

contains

  subroutine solve_advection_diffusion_fv(field_name, state)
    !!< Construct and solve the advection-diffusion equation for the given
    !!< field unsing element centred finite volumes.
    
    !! Name of the field to be solved for.
    character(len=*), intent(in) :: field_name
    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state

    !! Tracer to be solved for.
    type(scalar_field), pointer :: T, T_old
    !! Change in T over one timestep.
    type(scalar_field) :: delta_T
    !! System matrix.
    type(csr_matrix) :: matrix
    !! Right hand side vector.
    type(scalar_field) :: rhs
    !! Sparsity of advection_diffusion matrix
    type(csr_sparsity) :: sparsity

    type(tensor_field), pointer :: diffusivity

    T=>extract_scalar_field(state, field_name)
    T_old=>extract_scalar_field(state, "Old"//field_name)

    ! Reset T to value at the beginning of the timestep.
    call set(T, T_old)

    if(has_tensor_field(state, trim(field_name)//"Diffusivity")) then
      diffusivity=>extract_tensor_field(state, trim(field_name)//"Diffusivity")
      sparsity=make_sparsity_transpose(t%mesh, diffusivity%mesh, trim(field_name)//"AdvectionDiffusionSparsity")
    else
      sparsity=make_sparsity(t%mesh, t%mesh, trim(field_name)//"AdvectionSparsity")
    end if

    call allocate(matrix, sparsity, name=trim(field_name)//"AdvectionDiffusionMatrix") ! Add data space to the sparsity

    ! Ensure delta_T inherits options from T.
    call allocate(delta_T, T%mesh, trim(field_name)//"Change")
    delta_T%option_path = T%option_path
    call allocate(rhs, T%mesh, trim(field_name)//"RHS")
    
    call construct_advection_diffusion_fv(matrix, rhs, field_name, state) 

    ! Apply strong dirichlet boundary conditions.
    ! This is for big spring boundary conditions.
    call apply_dirichlet_conditions(matrix, rhs, T, dt)

    call zero(delta_T) ! Impose zero initial guess.
    ! Solve for the change in T.
    call petsc_solve(delta_T, matrix, rhs)

    ! Add the change in T to T.
    call addto(T, delta_T, dt)

    call deallocate(delta_T)
    call deallocate(matrix)
    call deallocate(rhs)
    call deallocate(sparsity)

  end subroutine solve_advection_diffusion_fv

  subroutine construct_advection_diffusion_fv(big_m, rhs, field_name,&
       & state, mass) 
    !!< Construct the advection_diffusion equation for finite volumes in
    !!< acceleration form.
    !!< 
    !!< If mass is provided then the mass matrix is not added into big_m or
    !!< rhs. It is instead returned as mass. This may be useful for testing
    !!< or for solving equations otherwise than in accelleration form.
    
    !! Main advection_diffusion matrix.    
    type(csr_matrix), intent(inout) :: big_m
    !! Right hand side vector.
    type(scalar_field), intent(inout) :: rhs
    
    !! Name of the field to be advected.
    character(len=*), intent(in) :: field_name
    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state
    !! Optional separate mass matrix.
    type(csr_matrix), intent(inout), optional :: mass

    !! Position, and velocity fields.
    type(vector_field) :: X, U, U_nl, U_mesh
    !! Tracer to be solved for.
    type(scalar_field) :: T
    !! Diffusivity
    type(tensor_field) :: Diffusivity

    !! Source and absorption
    type(scalar_field) :: Source, Abs

    !! Element index
    integer :: ele

    !! Status variable for field extraction.
    integer :: stat

    !! Shape function for the auxiliary variable in the second order
    !! operator. 
    type(element_type), pointer :: q_shape
    type(element_type), pointer :: T_shape

    !! Mesh for auxiliary variable
    type(mesh_type) :: q_mesh

    type(csr_sparsity) :: grad_sparsity, div_sparsity, sparsity_q
    type(block_csr_matrix) :: div_m, grad_t_m
    type(vector_field) :: grad_rhs
    type(scalar_field) :: diff_rhs
    type(csr_matrix) :: inverse_mass, D_m
      
    !! Field over the entire surface mesh containing bc values:
    type(scalar_field) :: bc_value
    !! Integer array of all surface elements indicating bc type
    !! (see below call to get_entire_boundary_condition):
    integer, dimension(:), allocatable :: bc_type

    integer :: quaddegree

    ewrite(1,*) "Writing finite volume advection-diffusion equation for "&
         &//trim(field_name)

    ! These names are based on the CGNS SIDS.
    T=extract_scalar_field(state, field_name)
    X=extract_vector_field(state, "Coordinate")
    U=extract_vector_field(state, "Velocity")

    U_nl=extract_vector_field(state, "NonlinearVelocity")

    U_mesh=extract_vector_field(state, "GridVelocity")
!    Abs=>extract_scalar_field(state, "Absorption")

    Diffusivity=extract_tensor_field(state, trim(field_name)//"Diffusivity"&
         &, stat=diff_stat)
    if (diff_stat/=0) then
       call allocate(Diffusivity, T%mesh, trim(field_name)//"Diffusivity",&
            FIELD_TYPE_CONSTANT)
       call zero(Diffusivity)
    else
       ! Grab an extra reference to cause the deallocate below to be safe.
       call incref(Diffusivity)
    end if

    ! Determine the scheme to use to discretise diffusivity.
    if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
         &finite_volume/diffusion_scheme/bassi_rebay")) then
       diffusion_scheme=BASSIREBAY
    else if (have_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
         &finite_volume/diffusion_scheme/ldg")) then
       diffusion_scheme=LDG
    else
       FLAbort("Unknown diffusion scheme.")
    end if

    Source=extract_scalar_field(state, trim(field_name)//"Source"&
         &, stat=stat)
    if (stat/=0) then
       call allocate(Source, T%mesh, trim(field_name)//"Source",&
            FIELD_TYPE_CONSTANT)
       call zero(Source)
    else
       ! Grab an extra reference to cause the deallocate below to be safe.
       call incref(Source)
    end if

    Abs=extract_scalar_field(state, trim(field_name)//"Absorption"&
         &, stat=stat)
    if (stat/=0) then
       call allocate(Abs, T%mesh, trim(field_name)//"Absorption",&
            FIELD_TYPE_CONSTANT)
       call zero(Abs)
    else
       ! Grab an extra reference to cause the deallocate below to be safe.
       call incref(Abs)
    end if

    ! Retrieve scalar options from the options dictionary.
    call get_option(trim(T%option_path)//&
         &"/prognostic/temporal_discretisation/theta", theta)
    call get_option("/timestepping/timestep", dt)
    call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
         &conservative_advection", beta)

    penalise = .false.
    call get_option(trim(T%option_path)//"/prognostic/spatial_discretisation/&
           &finite_volume/diffusion_scheme/bassi_rebay/penalty", penalty, default = 0.0, stat=stat)
    if(stat==0) penalise = .true.

    q_mesh=diffusivity%mesh

    div_sparsity = make_sparsity(t%mesh, q_mesh, trim(field_name)//"GradientSparsity")
    grad_sparsity = make_sparsity(t%mesh, q_mesh, trim(field_name)//"GradientSparsity")
    sparsity_q = make_sparsity(q_mesh, q_mesh, trim(field_name)//"DiffusivityMassSparsity")
    call allocate(div_m, sparsity=grad_sparsity, blocks=(/1, mesh_dim(T)/), &
                  name=trim(field_name)//"AuxiliaryDivergenceMatrix")
    call zero(div_m)
    call allocate(grad_t_m, sparsity=grad_sparsity, blocks=(/1, mesh_dim(T)/), &
                  name=trim(field_name)//"AuxiliaryGradientMatrixTransposed")
    call zero(grad_t_m)
    call allocate(grad_rhs, mesh_dim(T), q_mesh, trim(field_name)//"AuxiliaryGradientRHS")
    call zero(grad_rhs)
    call allocate(diff_rhs, T%mesh, trim(field_name)//"DiffusionRHS")
    call zero(diff_rhs)
    call allocate(inverse_mass, sparsity_q, name=trim(field_name)//"InverseMass")
    call zero(inverse_mass)
    call allocate(D_m, big_m%sparsity, name=trim(field_name)//"DiffusivityMatrix")
    call zero(D_m)

    assert(has_faces(X%mesh))
    assert(has_faces(T%mesh))

    call zero(big_m)
    call zero(RHS)
    if (present(mass)) call zero(mass)
    
    ! Enquire about boundary conditions we're interested in
    ! Returns an integer array bc_type over the surface elements
    ! that indicates the bc type (in the order we specified, i.e.
    ! BCTYPE_WEAKDIRICHLET=1)
    allocate( bc_type(1:surface_element_count(T)) )
    call get_entire_boundary_condition(T, &
       & (/"weakdirichlet"/), &
       & bc_value, bc_type)    

    element_loop: do ele=1,element_count(T)

       call construct_adv_diff_element_fv(ele, big_m, rhs, &
            & state, X, T, U_nl, U_mesh, Source, Abs,&
            & Diffusivity, bc_value, bc_type, &
            & q_mesh, div_m, grad_t_m, grad_rhs, inverse_mass, mass, &
            & D_m, diff_rhs)

    end do element_loop

    if(diff_stat==0) then
      call assemble_diffusion_m_fv(big_m, rhs, T, D_m, diff_rhs, &
                                  div_m, grad_t_m, grad_rhs, &
                                  diffusivity, inverse_mass, dt, theta)
    end if

    ! Drop any extra field references.
    call deallocate(Diffusivity)
    call deallocate(Source)
    call deallocate(Abs)
    call deallocate(div_m)
    call deallocate(grad_t_m)
    call deallocate(sparsity_q)
    call deallocate(grad_sparsity)
    call deallocate(div_sparsity)
    call deallocate(diff_rhs)
    call deallocate(grad_rhs)
    call deallocate(D_m)
    call deallocate(inverse_mass)

  end subroutine construct_advection_diffusion_fv

  subroutine assemble_diffusion_m_fv(big_m, rhs, T, D_m, diff_rhs, &
                                     div_m, grad_t_m, grad_rhs, &
                                     diffusivity, inverse_mass, dt, theta)

    type(csr_matrix), intent(inout) :: big_m
    type(scalar_field), intent(inout) :: rhs
    type(scalar_field), intent(inout) :: T

    type(csr_matrix), intent(inout) :: D_m
    type(scalar_field), intent(inout) :: diff_rhs

    type(block_csr_matrix), intent(inout) :: div_m, grad_t_m
    type(vector_field), intent(inout) :: grad_rhs

    type(tensor_field), intent(in) :: diffusivity
    type(csr_matrix), intent(in) :: inverse_mass

    real, intent(in) :: dt, theta

    type(csr_matrix) :: grad_m_t_block, ldiv_diff, D_m_temp, DM_m_temp
    type(scalar_field) :: diffusivity_dd, grad_rhs_dim, temprhs

    integer :: dim1, dim2, row
    real, dimension(:), pointer :: row_val
    integer, dimension(:), pointer :: row_indices

    logical :: isotropic

    type(scalar_field) :: MT_old

    ewrite(1,*) 'in assemble_diffusion_m_fv'

    call allocate(MT_old, rhs%mesh, name="MT_oldProduct" )
    call zero(MT_old)

    call allocate(ldiv_diff, grad_t_m%sparsity, name = "LocalDivergenceByLumpedMass")
    call zero(ldiv_diff)

    call allocate(temprhs, diff_rhs%mesh, name= "AuxilliaryTempRHS")
    call zero(temprhs)

    ! an optimisation that reduces the number of matrix multiplies if we're isotropic
    isotropic = isotropic_field(diffusivity)

    do dim1 = 1, diffusivity%dim
      do dim2 = 1, diffusivity%dim

        if((isotropic).and.(dim1/=dim2)) cycle

        diffusivity_dd = extract_scalar_field(diffusivity, dim1, dim2)

        ! multiply the divergence matrix by the diffusivity and the inverse mass matrix
        do row = 1, size(div_m, 1)
          row_indices=>row_m_ptr(div_m, row)
          row_val=>row_val_ptr(div_m, 1, dim1, row)
          call set(ldiv_diff, (/row/), row_indices, &
                  spread((row_val*node_val(diffusivity_dd, row_indices)), 1, 1))
        end do

        DM_m_temp = matmul_T(ldiv_diff, inverse_mass, model=grad_t_m%sparsity)

        ! multiply the divergence and mass by the gradient matrix
        grad_m_t_block = block(grad_t_m, 1, dim2)
        D_m_temp = matmul_T(DM_m_temp, grad_m_t_block, model=D_m%sparsity)

        call addto(D_m, D_m_temp)

        ! multiply the divergence matrix by the gradient of any boundary terms
        grad_rhs_dim = extract_scalar_field(grad_rhs, dim2)
        ewrite_minmax(grad_rhs_dim%val)
        call mult(temprhs, DM_m_temp, grad_rhs_dim)

        call addto(diff_rhs, temprhs)

        call deallocate(D_m_temp)
        call deallocate(DM_m_temp)

      end do
    end do

    call mult(MT_old, D_m, T)
    call addto(rhs, MT_old, -1.0)
    call addto(rhs, diff_rhs, -1.0)

    call addto(big_m, D_m, theta*dt)

    call deallocate(temprhs)
    call deallocate(ldiv_diff)
    call deallocate(MT_old)

  end subroutine assemble_diffusion_m_fv

  subroutine construct_adv_diff_element_fv(ele, big_m,rhs,&
       & state, X, T, U_nl, U_mesh, Source, Abs, Diffusivity,&
       & bc_value, bc_type, &
       & q_mesh, div_m, grad_t_m, grad_rhs, inverse_mass, mass, &
       & D_m, diff_rhs)
    !!< Construct the advection_diffusion equation for finite volumes in
    !!< acceleration form.
    implicit none
    !! Index of current element
    integer :: ele
    !! Main advection_diffusion matrix.
    type(csr_matrix), intent(inout) :: big_m
    !! Right hand side vector.
    type(scalar_field), intent(inout) :: rhs
    !! Field over the entire surface mesh containing bc values:
    type(scalar_field), intent(in) :: bc_value
    !! Integer array of all surface elements indicating bc type
    !! (see below call to get_entire_boundary_condition):
    integer, dimension(:), intent(in) :: bc_type

    !! Auxiliary variable mesh
    type(mesh_type), intent(in) :: q_mesh

    type(block_csr_matrix), intent(inout) :: div_m, grad_t_m
    type(vector_field), intent(inout) :: grad_rhs
    type(csr_matrix), intent(inout) :: inverse_mass

    !! Optional separate mass matrix.
    type(csr_matrix), intent(inout), optional :: mass

    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state   

    type(csr_matrix), intent(inout) :: D_m
    type(scalar_field), intent(inout) :: diff_rhs

    !! Position and velocity.
    type(vector_field) :: X, U_nl, U_mesh

    type(scalar_field) :: T, Source, Abs
    !! Diffusivity
    type(tensor_field) :: Diffusivity

    
    ! Bilinear forms.
    real, dimension(ele_loc(T,ele), ele_loc(T,ele)) :: &
         mass_mat
    real, dimension(ele_loc(T,ele), ele_loc(T,ele)) ::  Abs_mat
    real, dimension(ele_loc(q_mesh,ele), ele_loc(q_mesh,ele)) :: Q_inv 
!     real, dimension(mesh_dim(T), ele_loc(q_mesh,ele), &
!          ele_and_faces_loc(T,ele)) :: Div_T_mat
    real, dimension(mesh_dim(T), ele_loc(t,ele), &
         ele_and_faces_loc(q_mesh,ele)) :: Div_mat
    real, dimension(mesh_dim(T), ele_and_faces_loc(T,ele), &
         ele_loc(q_mesh,ele)) :: Grad_T_mat
    logical, dimension(ele_and_faces_loc(T,ele), ele_loc(q_mesh,ele)) :: grad_t_mask
    real, dimension(mesh_dim(T), ele_loc(T,ele)) :: grad_rhs_local

    real, dimension(ele_loc(T,ele), ele_and_faces_loc(T,ele)) :: penalty_mat
    real, dimension(ele_loc(T,ele)) :: diff_rhs_local

    ! Local assembly matrices.
    real, dimension(ele_loc(T,ele), ele_loc(T,ele)) :: l_T_mat

    ! Local node number map for 2nd order element.
    integer, dimension(ele_and_faces_loc(q_mesh,ele)) :: q_local_glno
    integer, dimension(ele_and_faces_loc(T,ele)) :: t_local_glno

    ! Local variables.
    
    ! Neighbour element, face and neighbour face.
    integer :: ele_2, face, face_2
    ! Count variable for loops over dimension.
    integer :: dim1, dim2
    ! Loops over faces.
    integer :: ni
    ! Array bounds for faces of the 2nd order element.
    integer :: q_start, t_start, q_finish, t_finish
    
    ! Variable transform times quadrature weights.
    real, dimension(ele_ngi(T,ele)) :: detwei
    ! Transform from local to physical coordinates.
    real, dimension(U_nl%dim, U_nl%dim, ele_ngi(T,ele)) :: invJ 
    ! Transformed gradient function for tracer.
    real, dimension(ele_loc(T, ele), ele_ngi(T, ele), mesh_dim(T)) :: dt_t
    ! Transformed gradient function for velocity.
    real, dimension(ele_loc(U_nl, ele), ele_ngi(U_nl, ele), mesh_dim(T)) ::&
         & du_t 
    ! Transformed gradient function for auxiliary variable.
    real, dimension(ele_loc(q_mesh,ele), ele_ngi(q_mesh,ele), mesh_dim(T)) :: dq_t
    ! Different velocities at quad points.
    real, dimension(U_nl%dim, ele_ngi(U_nl, ele)) :: u_nl_q, u_mesh_q
    real, dimension(ele_ngi(U_nl, ele)) :: u_div_q, u_nl_div_q

    ! Node and shape pointers.
    integer, dimension(:), pointer :: t_ele, q_ele
    type(element_type), pointer :: t_shape, u_shape, q_shape
    ! Neighbours of this element.
    integer, dimension(:), pointer :: neigh
    ! Whether the tracer field is continuous.
    logical :: dg

    logical :: boundary_element

    integer :: iloc, oloc, row, col

    dg=continuity(T)<0

    !----------------------------------------------------------------------
    ! Establish local node lists
    !----------------------------------------------------------------------

    T_ele=>ele_nodes(T,ele)  ! field
    q_ele=>ele_nodes(q_mesh, ele) ! diffusivity

    t_local_glno=0
    t_local_glno(:size(t_ele))=t_ele ! Diffusivity node list.

    q_local_glno=0
    q_local_glno(:size(q_ele))=q_ele ! Diffusivity node list.

    !----------------------------------------------------------------------
    ! Establish local shape functions
    !----------------------------------------------------------------------

    t_shape=>ele_shape(T, ele)
    u_shape=>ele_shape(U_nl, ele)
    q_shape=>ele_shape(q_mesh, ele)

    ! Transform Tracer derivatives and weights into physical space. If
    ! necessary, grab invJ as well.
    call transform_to_physical(X, ele,&
        & t_shape , dshape=dt_t, detwei=detwei)

    ! Transform U_nl derivatives and weights into physical space.
    call transform_to_physical(X, ele,&
         & u_shape , dshape=du_t)

    ! Transform q derivatives into physical space.
    call transform_to_physical(X, ele,&
         & q_shape , dshape=dq_t)
        
    !----------------------------------------------------------------------
    ! Construct element-wise quantities.
    !----------------------------------------------------------------------

    ! Introduce grid velocities in non-linear terms. 
    ! Advecting velocity at quadrature points.
    U_nl_q=ele_val_at_quad(U_nl,ele)
    ! Divergence of advecting velocity.
    U_nl_div_q=ele_div_at_quad(U_nl, ele, du_t)

    ! Mesh velocity at quadrature points.
    U_mesh_q=ele_val_at_quad(U_mesh,ele)
    ! Divergence of mesh movement.
    U_div_q=ele_div_at_quad(U_mesh, ele, du_t)
    
    !----------------------------------------------------------------------
    ! Construct bilinear forms.
    !----------------------------------------------------------------------

    ! Element density matrix.
    !  /
    !  | T T  dV
    !  / 
    mass_mat = shape_shape(T_shape, T_shape, detwei)

    ! Mesh movement terms still to be added.

    ! Absorption matrix.
    Abs_mat = shape_shape(T_shape, T_shape, detwei*ele_val_at_quad(Abs,ele))

    ! Tau Q = grad(u)
    Q_inv= shape_shape(q_shape, q_shape, detwei)
    call invert(Q_inv)

    call addto(inverse_mass, q_ele, q_ele, Q_inv)

    grad_t_mask=.false.
    Grad_T_mat=0.0
    Div_mat=0.0
    penalty_mat = 0.0
    diff_rhs_local = 0.0
    Grad_T_mat(:, :size(T_ele), :) = shape_dshape(T_shape, dq_t, detwei)
    Div_mat(:, :, :size(q_ele)) = -dshape_shape(dt_t, q_shape, detwei)
    grad_t_mask(:size(t_ele),:) = .true.

    grad_rhs_local = 0.0

    !----------------------------------------------------------------------
    ! Perform global assembly.
    !----------------------------------------------------------------------

    ! Right hand side field.
    call addto(RHS, t_ele, &
         ! Source term
          sum(detwei)*ele_val(Source,ele) &
    )
    
    ! Assemble matrix.
    
    ! Advection.
    l_T_mat= 0.0

    if (present(mass)) then
       ! Return mass separately.
       call addto(mass, t_ele, t_ele, mass_mat)
    else
       ! Put mass in the matrix.
       l_T_mat=l_T_mat+mass_mat
    end if

    call addto(big_m, t_ele, t_ele, l_T_mat)

    !-------------------------------------------------------------------
    ! Interface integrals
    !-------------------------------------------------------------------
    
    neigh=>ele_neigh(T, ele)
    ! Local node map counter.
    q_start=size(q_ele)+1
    t_start=size(t_ele)+1
    ! Flag for whether this is a boundary element.
    boundary_element=.false.

    neighbourloop: do ni=1,size(neigh)

       !----------------------------------------------------------------------
       ! Find the relevant faces.
       !----------------------------------------------------------------------
       
       ! These finding routines are outside the inner loop so as to allow
       ! for local stack variables of the right size in
       ! construct_add_diff_interface_fv.

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

       if (dg) then
          q_finish=q_start+face_loc(q_mesh, face_2)-1
          q_local_glno(q_start:q_finish)=face_global_nodes(q_mesh, face_2)

          t_finish=t_start+face_loc(t, face_2)-1
          t_local_glno(t_start:t_finish)=face_global_nodes(t, face_2)
       end if

       call construct_adv_diff_interface_fv(ele, ele_2, face, face_2, ni,&
            & big_m, rhs, Grad_T_mat, Div_mat, grad_t_mask, state, X, T, U_nl,&
            & bc_value, bc_type, &
            & U_mesh, q_mesh,  &
            & grad_rhs_local, penalty_mat, diff_rhs_local)
     
       if (dg) then
          q_start=q_start+face_loc(q_mesh, face_2)
          t_start=t_start+face_loc(t, face_2)
       end if
       
    end do neighbourloop

    do dim1 = 1, mesh_dim(T)
      call addto(div_m, 1, dim1, t_ele, q_local_glno, div_mat(dim1,:,:))
      call addto(grad_t_m, 1, dim1, t_local_glno, q_ele, grad_t_mat(dim1,:,:), mask=grad_t_mask)
      call addto(grad_rhs, dim1, T_ele, grad_rhs_local(dim1,:))
    end do
    call addto(D_m, t_ele, t_local_glno, penalty_mat)
    call addto(diff_rhs, t_ele, diff_rhs_local)

  end subroutine construct_adv_diff_element_fv

  subroutine construct_adv_diff_interface_fv(ele, ele_2, face, face_2, &
       ni, big_m, rhs, Grad_T_mat, Div_mat, grad_t_mask, state, X, T, U_nl,&
       & bc_value, bc_type, &
       & U_mesh, q_mesh, &
       & grad_rhs_local, penalty_mat, diff_rhs_local)
    !!< Construct the FV element boundary integrals on the ni-th face of
    !!< element ele.
    implicit none

    integer, intent(in) :: ele, ele_2, face, face_2, ni
    type(csr_matrix), intent(inout) :: big_m
    type(scalar_field), intent(inout) :: rhs
    real, dimension(:,:,:), intent(inout) :: Grad_T_mat, Div_mat
    logical, dimension(:,:), intent(inout) :: grad_t_mask
    real, dimension(:,:), intent(inout) :: grad_rhs_local
    type(state_type), intent(in) :: state
    ! We pass these additional fields to save on state lookups.
    type(vector_field), intent(in) :: X, U_nl, U_mesh
    type(scalar_field), intent(in) :: T
    !! Mesh of the auxiliary variable in the second order operator.
    type(mesh_type), intent(in) :: q_mesh
    real, dimension(:,:), intent(inout) :: penalty_mat
    real, dimension(:), intent(inout) :: diff_rhs_local
    !! Field over the entire surface mesh containing bc values:
    type(scalar_field), intent(in) :: bc_value
    !! Integer array of all surface elements indicating bc type
    !! (see below call to get_entire_boundary_condition):
    integer, dimension(:), intent(in):: bc_type

    ! Face objects and numberings.
    type(element_type), pointer :: T_shape, T_shape_2, q_shape, q_shape_2
    integer, dimension(face_loc(T,face)) :: T_face, T_face_l
    integer, dimension(face_loc(T,face_2)) :: T_face_2
    ! This has to be a pointer to work around a stupid gcc bug.
    integer, dimension(:), pointer :: q_face_l

    ! Note that both sides of the face can be assumed to have the same
    ! number of quadrature points.
    real, dimension(U_nl%dim, face_ngi(U_nl, face)) :: normal, u_nl_q,&
         & u_mesh_q, u_div_q
    logical, dimension(face_ngi(U_nl, face)) :: inflow
    ! Variable transform times quadrature weights.
    real, dimension(face_ngi(T,face)) :: detwei

    ! Bilinear forms
    real, dimension(face_loc(T,face),face_loc(T,face)) :: nnAdvection_out
    real, dimension(face_loc(T,face),face_loc(T,face_2)) :: nnAdvection_in
    
    integer :: dim, q_start, t_start, q_finish, t_finish
    logical :: boundary, dirichlet

    t_face=face_global_nodes(T, face)

    t_face_l=face_local_nodes(T, face)
    t_shape=>face_shape(T, face)

    t_face_2=face_global_nodes(T, face_2)
    t_shape_2=>face_shape(T, face_2)
    
    q_face_l=>face_local_nodes(q_mesh, face)
    q_shape=>face_shape(q_mesh, face)

    q_shape_2=>face_shape(q_mesh, face_2)

    ! Boundary nodes have both faces the same.
    boundary=(face==face_2)
    dirichlet=.false.
    if (boundary) then
       if (bc_type(face)==BCTYPE_WEAKDIRICHLET) then
          dirichlet=.true.
       end if
    end if

    !----------------------------------------------------------------------
    ! Change of coordinates on face.
    !----------------------------------------------------------------------

    call transform_facet_to_physical(X, face, &
         &                          detwei_f=detwei,&
         &                          normal=normal)

    !----------------------------------------------------------------------
    ! Construct element-wise quantities.
    !----------------------------------------------------------------------

    U_nl_q=0.0
    U_mesh_q=0.0
    U_div_q=0.0

    ! Introduce grid velocities in non-linear terms. 

    ! Advecting velocity at quadrature points.
    U_nl_q=0.5*(face_val_at_quad(U_nl, face)+face_val_at_quad(U_nl, face_2))
    
    ! Mesh velocity at quadrature points.
    U_mesh_q=face_val_at_quad(U_mesh, face)
    ! Divergence of mesh movement.
    !       U_div_q=ele_face_div_at_quad(U_mesh, ele, ni)
    
    ! Inflow is true if the flow at this gauss point is directed
    ! into this element.
    inflow= sum((U_nl_q-U_mesh_q)*normal,1)<0.0

    !----------------------------------------------------------------------
    ! Construct bilinear forms.
    !----------------------------------------------------------------------

    ! Calculate outflow boundary integral.
    nnAdvection_out=shape_shape(T_shape, T_shape,  &
         &                        (merge(1.0,0.0,.not.inflow) &
         ! Jump condition interior term.
         &                               -(1.0-beta)) &
         &                          *sum((U_nl_q)*normal,1) &
         &                          *detwei) 
    
    nnAdvection_in=shape_shape(T_shape, T_shape_2, &
         &                       merge(1.0,0.0,inflow) &
         &                         *sum((U_nl_q)*normal,1) &
         &                         *detwei) 

    ! Boundary term in grad_U.
    !   /
    !   | q, u, normal dx
    !   /
    q_start=ele_loc(q_mesh, ele)+(ni-1)*face_loc(q_mesh, face_2)+1
    q_finish=q_start+face_loc(q_mesh, face_2)-1

    t_start=ele_loc(t, ele)+(ni-1)*face_loc(t, face_2)+1
    t_finish=t_start+face_loc(t, face_2)-1

    select case(diffusion_scheme)
    case(BASSIREBAY)
      call bassi_rebay_diffusion
    case(LDG)
      call ldg_diffusion
    end select

    !----------------------------------------------------------------------
    ! Perform global assembly.
    !----------------------------------------------------------------------

    ! Insert advection in matrix.

    ! Outflow boundary integral.
    call addto(big_M, T_face, T_face,&
            nnAdvection_out*dt*theta)
       
    if (.not.dirichlet) then
       ! Inflow boundary integral.
       call addto(big_M, T_face, T_face_2,&
            nnAdvection_in*dt*theta)
    end if

    ! Insert advection in RHS.
       
    if (.not.dirichlet) then
       ! For interior interfaces this is the upwinding term. For a Neumann
       ! boundary it's necessary to apply downwinding here to maintain the
       ! surface integral. Fortunately, since face_2==face for a boundary
       ! this is automagic.

       call addto(RHS, T_face, &
            ! Outflow boundary integral.
            -matmul(nnAdvection_out,face_val(T,face))&
            ! Inflow boundary integral.
            -matmul(nnAdvection_in,face_val(T,face_2)))

    else

       ! Inflow and outflow of Dirichlet value.
       call addto(RHS, T_face, &
            -matmul(nnAdvection_in,&
            ele_val(bc_value, face)))
       
       ! The interior integral is still interior!
       call addto(RHS, T_face, &
            -matmul(nnAdvection_out,face_val(T,face)))

    end if

  contains

    subroutine bassi_rebay_diffusion
      
      do dim=1,mesh_dim(T)

         if(.not.boundary) then
            ! Internal face.
            div_mat(dim, T_face_l, q_face_l)=&
                 div_mat(dim, T_face_l, q_face_l) &
!                  +0.5*sum(detwei*normal(dim,:))/size(q_face_l)
                 +0.5*shape_shape(t_shape, q_shape, detwei*normal(dim,:))

            grad_t_mat(dim, t_face_l, q_face_l)=&
                 grad_t_mat(dim, t_face_l, q_face_l) &
                 -0.5*shape_shape(t_shape, q_shape, detwei*normal(dim,:))
            grad_t_mask(t_face_l, q_face_l) = .true.

            penalty_mat(t_face_l, t_face_l) = penalty_mat(t_face_l, t_face_l) &
                  +penalty*shape_shape(t_shape, t_shape, detwei)

            ! External face.
            div_mat(dim, T_face_l, q_start:q_finish)=&
                 div_mat(dim, T_face_l, q_start:q_finish) &
!                  +0.5*sum(detwei*normal(dim,:))/size(q_face_l)
                 +0.5*shape_shape(t_shape, q_shape_2, detwei*normal(dim,:))

            grad_t_mat(dim, t_start:t_finish, q_face_l)=&
                 grad_t_mat(dim, t_start:t_finish, q_face_l) &
                 -0.5*shape_shape(t_shape_2, q_shape, detwei*normal(dim,:))
            grad_t_mask(t_start:t_finish, q_face_l) = .true.

            penalty_mat(t_face_l, t_start:t_finish) = penalty_mat(t_face_l, t_start:t_finish) &
                  -penalty*shape_shape(t_shape, t_shape_2, detwei)

         else
!             ! Boundary case. Put the whole integral in the external bit.
            if(dirichlet) then
              grad_rhs_local(dim, T_face_l) = grad_rhs_local(dim, T_face_l) &
                 -sum(detwei*normal(dim,:))*ele_val(bc_value,face)

              div_mat(dim, T_face_l, q_face_l)=&
                  div_mat(dim, T_face_l, q_face_l) &
                  +sum(detwei*normal(dim,:))

              penalty_mat(t_face_l, t_face_l) = penalty_mat(t_face_l, t_face_l) &
                  +penalty*shape_shape(t_shape, t_shape, detwei)

              diff_rhs_local(t_face_l) = diff_rhs_local(t_face_l) &
                  -penalty*sum(detwei)*ele_val(bc_value, face)

            else

              ! External face.
              grad_t_mat(dim, t_face_l, q_face_l)=&
                 grad_t_mat(dim, t_face_l, q_face_l) &
                 -shape_shape(t_shape, q_shape, detwei*normal(dim,:))
              grad_t_mask(t_face_l, q_face_l) = .true.

            end if


         end if
      end do

    end subroutine bassi_rebay_diffusion

    subroutine ldg_diffusion
      
      do dim=1,mesh_dim(T)

         if(.not.boundary) then
            ! Internal face.
            div_mat(dim, T_face_l, q_face_l)=&
                 div_mat(dim, T_face_l, q_face_l) &
                 +shape_shape(t_shape, q_shape, detwei*normal(dim,:))

!             grad_t_mat(dim, t_face_l, q_face_l)=&
!                  grad_t_mat(dim, t_face_l, q_face_l) &
!                  -shape_shape(t_shape, q_shape, detwei*normal(dim,:))
!             grad_t_mask(t_face_l, q_face_l) = .true.

            ! External face.
!             div_mat(dim, T_face_l, q_start:q_finish)=&
!                  div_mat(dim, T_face_l, q_start:q_finish) &
!                  +shape_shape(t_shape, q_shape_2, detwei*normal(dim,:))

            grad_t_mat(dim, t_start:t_finish, q_face_l)=&
                 grad_t_mat(dim, t_start:t_finish, q_face_l) &
                 -shape_shape(t_shape_2, q_shape, detwei*normal(dim,:))
            grad_t_mask(t_start:t_finish, q_face_l) = .true.

         else
!             ! Boundary case. Put the whole integral in the external bit.
            if(dirichlet) then
              grad_rhs_local(dim, T_face_l) = grad_rhs_local(dim, T_face_l) &
                 -sum(detwei*normal(dim,:))*ele_val(bc_value,face)

              div_mat(dim, T_face_l, q_face_l)=&
                  div_mat(dim, T_face_l, q_face_l) &
                  +sum(detwei*normal(dim,:))

            else

              ! External face.
              grad_t_mat(dim, t_face_l, q_face_l)=&
                 grad_t_mat(dim, t_face_l, q_face_l) &
                 -shape_shape(t_shape, q_shape, detwei*normal(dim,:))
              grad_t_mask(t_face_l, q_face_l) = .true.

            end if


         end if
      end do

    end subroutine ldg_diffusion

  end subroutine construct_adv_diff_interface_fv

end module advection_diffusion_fv
