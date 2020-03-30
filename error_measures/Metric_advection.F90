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

module metric_advection

  use fldebug
  use vector_tools
  use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
  use futils, only: int2str
  use elements
  use spud
  use parallel_tools
  use sparse_tools
  use shape_functions
  use cv_faces
  use transform_elements
  use fetools
  use unittest_tools
  use fields
  use state_module
  use vtk_interfaces
  use sparse_matrices_fields
  use solvers
  use boundary_conditions
  use merge_tensors
  use edge_length_module
  use sparsity_patterns
  use cv_shape_functions
  use field_options, only: get_coordinate_field
  use cvtools
  use cv_options
  use cv_upwind_values
  use cv_face_values, only: theta_val, evaluate_face_val
  use sparsity_patterns_meshes
  use fefields, only: compute_cv_mass
  use state_fields_module, only: get_cv_mass
  use diagnostic_variables, only: field_tag
  use diagnostic_fields, only: calculate_diagnostic_variable
  use form_metric_field
  use populate_state_module
  use adaptive_timestepping

  implicit none

  private
  public :: form_advection_metric, initialise_metric_advection, use_metric_advection

  logical :: use_metric_advection

contains

  subroutine initialise_metric_advection
    if (have_option("/mesh_adaptivity/hr_adaptivity/metric_advection")) then
      use_metric_advection = .true.
    else
      use_metric_advection = .false.
    end if
  end subroutine initialise_metric_advection

  function expand_idx(i, j, dim) result(k)
    integer, intent(in) :: i, j, dim
    integer :: k

    k = i + dim * (j - 1)
  end function expand_idx

  subroutine form_advection_metric(tfield, state)
    !! Name of the field to be solved for.
    type(tensor_field), intent(inout) :: tfield
    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state

    ! Coordinate field
    type(vector_field), pointer :: x
    type(vector_field) :: x_tfield


    ! LHS equation matrix.
    type(csr_matrix) :: M, A_m
    ! sparsity structure to construct the matrices with
    type(csr_sparsity), pointer :: mesh_sparsity

    ! Change in tfield over one timestep.
    ! Right hand side vector, cv mass matrix, 
    ! locally iterated field (for advection iterations) 
    ! and local old field (for subcycling)
    type(tensor_field) :: rhs, advit_tfield, delta_tfield, accum_tfield
    type(tensor_field) :: l_tfield, tmp_tfield
    type(scalar_field), pointer :: t_cvmass
    type(scalar_field) :: cvmass
      
    ! local copy of option_path for solution field
    character(len=OPTION_PATH_LEN) :: option_path

    ! number of advection iterations and subcycles
    integer :: rk_iterations, adv_iterations, no_subcycles
    ! iterators
    integer :: rk_it, adv_it, i, sub, j, period_in_timesteps
    ! time (to output to file), timestep, iterations tolerance, subcycling timestep
    real :: actual_dt, adapt_dt, error, sub_dt
    real :: max_cfl, scale_adv

    ! degree of quadrature to use on each control volume face
    integer :: quaddegree
    ! control volume face information
    type(cv_faces_type) :: cvfaces
    ! control volume shape function for volume and boundary
    type(element_type) :: u_cvshape, u_cvbdyshape
    type(element_type) :: x_cvshape, x_cvbdyshape
    type(element_type) :: t_cvshape, t_cvbdyshape

    ! options wrappers for tfield
    type(cv_options_type) :: tfield_options

    ! success indicators?
    integer :: stat, cfl_stat
    ! type of courant number we want to use
    character(len=FIELD_NAME_LEN) :: cfl_type
    ! the courant number field
    type(scalar_field) :: cfl_no
    ! nonlinear and grid velocities
    type(vector_field), pointer :: nu, ug
    ! relative velocity
    type(vector_field) :: relu
    ! assume explicitness?
    logical :: explicit
    ! if we're subcycling how fast can we go?
    real :: max_sub_cfl
    ! construct the matrices
    logical :: getmat
    
    logical :: output_subcycle_vtus, output_final_vtus
    type(scalar_field) :: edgelen
    integer, save :: adaptcnt = 0
    
    ewrite(1,*) 'in metric advection'

    ! extract lots of fields:
    ! the actual thing we're trying to solve for
    option_path="/mesh_adaptivity/hr_adaptivity/metric_advection"

    ! now we can get the options for these fields
    ! handily wrapped in a new type...
    tfield_options=get_cv_options(option_path, tfield%mesh%shape%numbering%family, mesh_dim(tfield))

    output_subcycle_vtus = have_option(trim(option_path)//"/output/output_subcycle_vtus")
    output_final_vtus = have_option(trim(option_path)//"/output/output_final_vtus")
    
    ! extract fields from state
    nu=>extract_vector_field(state, "NonlinearVelocity")
    ug=>extract_vector_field(state, "GridVelocity")
    x=>extract_vector_field(state, "Coordinate")
    x_tfield = get_coordinate_field(state, tfield%mesh)

    ! find relative velocity
    call allocate(relu, nu%dim, nu%mesh, "RelativeVelocity")
    call set(relu, nu)
    call addto(relu, ug, -1.0)

    if(output_subcycle_vtus.or.output_final_vtus) then
      call allocate(edgelen, tfield%mesh, "Edge lengths")
    end if

    ! create control volume shape functions
    call get_option("/geometry/quadrature/controlvolume_surface_degree", &
                   quaddegree, default=1)
    cvfaces=find_cv_faces(vertices=ele_vertices(tfield, 1), &
                          dimension=tfield%mesh%shape%dim, &
                          polydegree=tfield%mesh%shape%degree, &
                          quaddegree=quaddegree)
    u_cvshape=make_cv_element_shape(cvfaces, nu%mesh%shape)
    x_cvshape=make_cv_element_shape(cvfaces, x%mesh%shape)
    t_cvshape=make_cv_element_shape(cvfaces, tfield%mesh%shape)
    u_cvbdyshape=make_cvbdy_element_shape(cvfaces, nu%mesh%faces%shape)
    x_cvbdyshape=make_cvbdy_element_shape(cvfaces, x%mesh%faces%shape)
    t_cvbdyshape=make_cvbdy_element_shape(cvfaces, tfield%mesh%faces%shape)
    
    ! get the mesh sparsity for the matrices
    mesh_sparsity=>get_csr_sparsity_firstorder(state, tfield%mesh, tfield%mesh)

    ! allocate the lhs matrix
    call allocate(M, mesh_sparsity, name="MetricMatrix")
    call allocate(A_m, mesh_sparsity, name="MetricAdvectionMatrix")

    ! allocate the rhs of the equation
    call allocate(rhs, tfield%mesh, name="MetricAdvectionRHS")

    ! allocate a field to store the subcycling iterations in
    call allocate(l_tfield, tfield%mesh, name="LocalMetric")
    call set(l_tfield, tfield)
    ! allocate a field to use to merge with the main field
    call allocate(tmp_tfield, tfield%mesh, name="TempMetric")
    call zero(tmp_tfield)
    ! allocate a field to store the locally iterated values in
    call allocate(advit_tfield, tfield%mesh, name="AdvIteratedMetric")
    call set(advit_tfield, tfield)
    ! allocate a field to store the accumulated Runge-Kutta slope in
    call allocate(accum_tfield, tfield%mesh, name="AccumAdvectedMetric")
    call zero(accum_tfield)
    ! allocate a field to store the change between the old and new values
    call allocate(delta_tfield, tfield%mesh, name="Delta_Metric")
    call zero(delta_tfield) ! Impose zero initial guess.
    ! Ensure delta_tfield inherits options from mesh_adaptivity
    delta_tfield%option_path = option_path

    explicit = have_option(trim(option_path)//"/explicit")

    ! find out how many iterations we'll be doing
    rk_iterations =1
!     ! Runge-Kutta 4
!     rk_iterations = 4
    
    call get_option(trim(option_path)//"/temporal_discretisation&
                    &/control_volumes/number_advection_iterations", &
                    adv_iterations, default=1)

    ! get the timestep information
    call get_option("/timestepping/timestep", actual_dt)
    call get_option("/mesh_adaptivity/hr_adaptivity/period", adapt_dt, stat)
    if (stat /= 0) then
      call get_option("/mesh_adaptivity/hr_adaptivity/period_in_timesteps", period_in_timesteps)
      adapt_dt = period_in_timesteps * actual_dt
    end if
    
    call get_option(trim(option_path)//"/temporal_discretisation&
                    &/scale_advection_time", scale_adv, default=1.1)
    adapt_dt = adapt_dt*scale_adv
      
    no_subcycles = 1
    ! are we subcycling?    
    call get_option(trim(option_path)//"/temporal_discretisation&    
                    &/number_advection_subcycles", no_subcycles, stat=stat)
    if(stat/=0) then
      ! have not specified a number of subcycles but perhaps we're using a 
      ! courant number definition?
      call get_option(trim(option_path)//"/temporal_discretisation/maximum_courant_number_per_subcycle", &
                      max_sub_cfl)
      call get_option(trim(option_path)//"/temporal_discretisation/maximum_courant_number_per_subcycle/courant_number[0]/name", cfl_type)
      call allocate(cfl_no, tfield%mesh, "CourantNumber")
      call calculate_diagnostic_variable(state, trim(cfl_type), cfl_no, dt=adapt_dt, &
         & option_path=trim(option_path)//"/temporal_discretisation/maximum_courant_number_per_subcycle/courant_number[0]")
      max_cfl = maxval(cfl_no%val)
      call allmax(max_cfl)
      ewrite(2,*) "max_cfl = ", max_cfl
      call deallocate(cfl_no)
      
      no_subcycles=ceiling(max_cfl/max_sub_cfl)
    end if
    
    sub_dt=adapt_dt/real(no_subcycles)

    ! find the cv mass that is used for the time term derivative
    t_cvmass => get_cv_mass(state, tfield%mesh)
    ewrite_minmax(t_cvmass)    
    
    call allocate(cvmass, tfield%mesh, "LocalCVMass")
    call set(cvmass, t_cvmass)
    ewrite_minmax(cvmass)    

    ewrite(2,*) 'no_subcycles = ', no_subcycles
    ewrite(2,*) 'rk_iterations = ', rk_iterations
    ewrite(2,*) 'adv_iterations = ', adv_iterations
    do sub=1,no_subcycles
      
      call zero(accum_tfield)
      
      do rk_it = 1, rk_iterations

        do adv_it = 1, adv_iterations
          getmat=(adv_it==1).and.(sub==1).and.(rk_it==1)
  
          ! record the value of tfield since the previous iteration
  
          if(explicit) then
            call set(cvmass, t_cvmass)
          else
            call zero(M)
            call addto_diag(M, t_cvmass)
          end if
          call zero(rhs)
  
          ! If we've passed the first iteration/subcycle so we have A_m don't enter the next step.
          ! Also if we're using first order upwinding (and not using a spatially varying theta 
          ! - limit_theta) so there's no need to assemble the
          ! nonlinear rhs (assuming we've enforced pivot theta = theta as cv_options should do)
          ! then just multiply things out
          if(getmat.or.(tfield_options%facevalue/=CV_FACEVALUE_FIRSTORDERUPWIND).or.(tfield_options%limit_theta)) then
            ! we need the matrix (probably the first iteration/subcycle) or we need
            ! the nonlinear rhs (if we're not using first order upwinding) or we're
            ! using a spatially varying theta that has to be multiplied into the
            ! assembly...
            ! so go into assembly
            call assemble_advection_m_cv(A_m, rhs, &
                                        advit_tfield, l_tfield, tfield_options, &
                                        cvfaces, x_cvshape, x_cvbdyshape, &
                                        u_cvshape, u_cvbdyshape, t_cvshape, t_cvbdyshape, &
                                        state, relu, x, x_tfield, &
                                        getmat, sub_dt, rk_it, delta_tfield)
          end if
  
          ! assemble it all into a coherent equation
          call assemble_field_eqn_cv(M, A_m, cvmass, rhs, &
                                    advit_tfield, l_tfield, &
                                    sub_dt, explicit, tfield_options)
  
          if(explicit) then
            do i = 1, delta_tfield%dim(1)
              do j = 1, delta_tfield%dim(2)
                delta_tfield%val(i,j,:) = rhs%val(i,j,:)/cvmass%val(:)
              end do
            end do
          else
            call zero(delta_tfield) ! Impose zero initial guess.
            ! Solve for the change in T.
            call petsc_solve(delta_tfield, M, rhs, symmetric=.true.) 
          end if
          
          call set(advit_tfield, l_tfield)  
          call addto(advit_tfield, delta_tfield, sub_dt)
          
        end do ! advection iterations
          
        if(rk_iterations > 1) then
          call addto(accum_tfield, delta_tfield, rk_coeff(adv_it))
        end if
      
      end do ! runge-kutta iterations
      
      if(rk_iterations>1) then
        call addto(l_tfield, accum_tfield, sub_dt / 6.0)
      else
        call set(l_tfield, advit_tfield)
      end if
      
      call set(tmp_tfield, l_tfield)
      call merge_tensor_fields(tfield, tmp_tfield)
      
      if(output_subcycle_vtus) then
        call get_edge_lengths(l_tfield, edgelen)
        call vtk_write_fields(trim("advected_metric_subcycle_")//int2str(adaptcnt), (sub-1), x, x%mesh, &
                                  sfields=(/edgelen/), tfields=(/l_tfield/))
      end if

    end do ! subcycle loop

    call bound_metric(tfield, state)
    
    if(output_final_vtus) then
      call get_edge_lengths(tfield, edgelen)
      call vtk_write_fields(trim("advected_metric_final"), adaptcnt, x, x%mesh, &
                                sfields=(/edgelen/), tfields=(/tfield/))
    end if

    if(output_subcycle_vtus.or.output_final_vtus) then
      adaptcnt = adaptcnt + 1
      call deallocate(edgelen)
    end if

    call deallocate(tmp_tfield)
    call deallocate(l_tfield)
    call deallocate(delta_tfield)
    call deallocate(advit_tfield)
    call deallocate(accum_tfield)
    call deallocate(rhs)
    call deallocate(M)
    call deallocate(A_m)
    call deallocate(u_cvbdyshape)
    call deallocate(x_cvbdyshape)
    call deallocate(t_cvbdyshape)
    call deallocate(u_cvshape)
    call deallocate(x_cvshape)
    call deallocate(t_cvshape)
    call deallocate(cvfaces)
    call deallocate(relu)
    call deallocate(x_tfield)
    call deallocate(cvmass)

  end subroutine form_advection_metric
  ! end of solution wrapping subroutines
  !************************************************************************
  !************************************************************************
  ! equation wrapping subroutines
  subroutine assemble_field_eqn_cv(M, A_m, cvmass, rhs, &
                                  tfield, oldtfield, &
                                  dt, explicit, tfield_options)

    ! This subroutine assembles the equation
    ! M(T^{n+1}-T^{n})/\Delta T = rhs
    ! for control volumes.
    ! By the time you get here M should already contain the mass
    ! components (and if back compatible the diffusional components)
    ! of the equation.

    ! inputs/outputs:
    ! lhs matrix
    type(csr_matrix), intent(inout) :: M
    ! matrix containing advective terms - to be incorporated
    ! into M during this subroutine
    type(csr_matrix), intent(inout) :: A_m
    type(scalar_field), intent(inout) :: cvmass
    ! rhs of equation
    type(tensor_field), intent(inout) :: rhs
    ! the field we are solving for
    type(tensor_field), intent(inout) :: tfield, oldtfield
    ! the timestep
    real, intent(in) :: dt
    ! we're explicit!
    logical, intent(in) :: explicit
    ! options wrappers for tfield
    type(cv_options_type), intent(in) :: tfield_options

    ! local memory:
    ! for all equation types:
    ! product of A_m and oldtfield
    type(scalar_field) :: A_mT_old

    type(scalar_field) :: oldtfield_scomp

    integer :: i, j

    ! allocate some memory for assembly
    call allocate(A_mT_old, rhs%mesh, name="A_mT_oldProduct" )
    
    do i=1,oldtfield%dim(1)
      do j=i,oldtfield%dim(2)
        ! construct rhs
        oldtfield_scomp = extract_scalar_field(oldtfield, i, j)
        call mult(A_mT_old, A_m, oldtfield_scomp)
        call addto(rhs, i, j, A_mT_old, -1.0)
        if (j /= i) then
          call addto(rhs, j, i, A_mT_old, -1.0)
        end if
      end do
    end do

    if((tfield_options%facevalue==CV_FACEVALUE_FIRSTORDERUPWIND).and.&
       (.not.tfield_options%limit_theta))then
       ! in this case the pivot solution is first order upwinding so
       ! A_m is all we need (i.e. we assume the rhs has been zeroed)
       ! also we're not using a spatially varying theta so A_m has been
       ! assembled excluding it (so it needs to be multiplied in now)
    
      ! [M + dt*theta*A_m](T^{n+1}-T^{n})/dt = rhs - A_m*T^{n}
      
      ! construct M
      if(.not.explicit) then
        call addto(M, A_m, tfield_options%theta*dt)
      end if
      
    else

      ! [M + dt*A_m](T^{n+1}-T^{n})/dt = rhs - A_m*T^{n}

      ! construct M
      if(.not.explicit) then
        call addto(M, A_m, dt)
      end if
      
    end if

    call deallocate(A_mT_old)

  end subroutine assemble_field_eqn_cv

  !************************************************************************
  !************************************************************************
  ! assembly subroutines 
  subroutine assemble_advection_m_cv(A_m, rhs, &
                                     tfield, oldtfield, tfield_options, &
                                     cvfaces, x_cvshape, x_cvbdyshape, &
                                     u_cvshape, u_cvbdyshape, t_cvshape, t_cvbdyshape, &
                                     state, relu, x, x_tfield, getmat, dt, rk_it, delta_tfield)

    ! This subroutine assembles the advection matrix and rhs for
    ! control volume field equations such that:
    ! A_m = div(\rho u T) - (1-beta)*T*div(\rho u)

    ! inputs/outputs:
    ! the advection matrix
    type(csr_matrix), intent(inout) :: A_m
    ! the rhs of the control volume field eqn
    type(tensor_field), intent(inout) :: rhs

    ! the field being solved for
    type(tensor_field), intent(inout), target :: tfield
    ! previous time level of the field being solved for
    type(tensor_field), intent(inout) :: oldtfield, delta_tfield
    ! a type containing all the tfield options
    type(cv_options_type), intent(in) :: tfield_options

    ! information about cv faces
    type(cv_faces_type), intent(in) :: cvfaces
    ! shape functions for region and surface
    type(element_type), intent(in) :: x_cvshape, x_cvbdyshape
    type(element_type), intent(in) :: u_cvshape, u_cvbdyshape
    type(element_type), intent(in) :: t_cvshape, t_cvbdyshape
    ! bucket full of fields
    type(state_type), intent(inout) :: state
    ! the relative velocity
    type(vector_field), intent(in) :: relu
    ! the coordinates
    type(vector_field), intent(inout) :: x, x_tfield
    ! logical indicating if the matrix should be constructed
    ! or if it exists already from a previous iteration
    logical, intent(in) :: getmat
    ! timestep
    real, intent(in) :: dt
    ! which runge kutta iteration are we on
    integer, intent(in) :: rk_it

    ! local memory:
    ! allocatable memory for coordinates, velocity, normals, determinants, nodes
    ! and the cfl number at the gauss pts and nodes
    real, dimension(:,:), allocatable :: x_ele, x_ele_bdy
    real, dimension(:,:), allocatable :: x_f, u_f, u_bdy_f
    real, dimension(:,:), allocatable :: normal, normal_bdy
    real, dimension(:), allocatable :: detwei, detwei_bdy
    real, dimension(:), allocatable :: normgi
    integer, dimension(:), pointer :: nodes, x_nodes, upwind_nodes
    integer, dimension(:), allocatable :: nodes_bdy
    real, dimension(:), allocatable :: cfl_ele

    ! allocatable memory for the values of the field at the nodes
    ! and on the boundary and for ghost values outside the boundary
    real, dimension(:), allocatable :: tfield_ele, oldtfield_ele
    real, dimension(:), allocatable :: tfield_ele_bdy, oldtfield_ele_bdy
    real, dimension(:), allocatable :: ghost_tfield_ele_bdy, ghost_oldtfield_ele_bdy

    ! some memory used in assembly of the face values
    real :: tfield_theta_val, tfield_pivot_val
    real :: tfield_face_val, oldtfield_face_val

    ! logical array indicating if a face has already been visited by the opposing node
    logical, dimension(:), allocatable :: notvisited

    ! loop integers
    integer :: ele, sele, iloc, oloc, face, gi, ggi


    ! mesh sparsity for upwind value matrices
    type(csr_sparsity), pointer :: mesh_sparsity
    ! upwind value matrices for the fields and densities
    type(csr_matrix)  :: tfield_upwind, &
          oldtfield_upwind

    ! incoming or outgoing flow
    real :: udotn, income
    logical :: inflow
    ! time and face discretisation
    real :: ptheta, ftheta, beta

    ! the type of the bc if integrating over domain boundaries
    integer, dimension(:), allocatable :: tfield_bc_type
    ! fields for the bcs over the entire surface mesh
    type(scalar_field) :: tfield_bc

    ! Bilinear form for the tensor twisting terms.
    real, dimension(tfield%dim(1), tfield%dim(2), ele_loc(tfield,1), ele_loc(tfield,1)) :: &
         Twist_mat
    ! Twist terms applied to explicit T_guess values.
    real, dimension(ele_loc(tfield,1)) :: Twist_rhs
    ! Transformed gradient function for velocity.
    real, dimension(ele_loc(relu, 1), ele_ngi(relu, 1), mesh_dim(tfield)) ::&
      & du_t
    type(element_type), pointer :: u_shape, t_shape

    type(scalar_field) :: tfield_scomp, oldtfield_scomp
    type(tensor_field), pointer :: U_nl_J
    type(tensor_field) :: bounded_U_nl_J

    real, dimension(:,:), allocatable :: mat_local
    real, dimension(:,:,:), allocatable :: rhs_local
    
    integer :: i, j, k, dimi, dimj

    real, dimension(relu%dim, relu%dim) :: id, value, strain, rotate
    type(tensor_field), pointer :: max_eigenbound
    real :: frob_norm
    integer :: node

    ewrite(2,*) 'assemble_advection_m_cv'

    ! allocate upwind value matrices
    mesh_sparsity=>get_csr_sparsity_firstorder(state, tfield%mesh, tfield%mesh)

    call allocate(tfield_upwind, mesh_sparsity, name="TFieldUpwindValues")
    call allocate(oldtfield_upwind, mesh_sparsity, name="OldTFieldUpwindValues")

    ! allocate memory for assembly
    allocate(x_ele(x%dim,ele_loc(x,1)), &
             x_f(x%dim, x_cvshape%ngi), &
             u_f(relu%dim, u_cvshape%ngi), &
             detwei(x_cvshape%ngi), &
             normal(x%dim, x_cvshape%ngi), &
             normgi(x%dim))
    allocate(cfl_ele(ele_loc(tfield,1)), &
             tfield_ele(ele_loc(tfield,1)), &
             oldtfield_ele(ele_loc(oldtfield,1)))
    allocate(notvisited(x_cvshape%ngi))
    allocate(mat_local(ele_loc(tfield,1), ele_loc(tfield,1)), &
             rhs_local(tfield%dim(1), tfield%dim(2), ele_loc(tfield,1)))

    ! Clear memory of arrays being designed
    if(getmat) call zero(A_m)

    if((tfield_options%facevalue==CV_FACEVALUE_FIRSTORDERUPWIND).and.&
        (.not.tfield_options%limit_theta))then
      dimi = 1
      dimj = 1
    else
      dimi = tfield%dim(1)
      dimj = tfield%dim(2)
    end if
    
    do i=1,dimi
      do j=i,dimj

        tfield_scomp = extract_scalar_field(tfield, i, j)
        oldtfield_scomp = extract_scalar_field(oldtfield, i, j)

        ! does the field need upwind values
        if(need_upwind_values(tfield_options)) then

          call find_upwind_values(state, x_tfield, tfield_scomp, tfield_upwind, &
                                  oldtfield_scomp, oldtfield_upwind, &
                                  option_path=trim(tfield%option_path))

        else

          call zero(tfield_upwind)
          call zero(oldtfield_upwind)

        end if

        ! some temporal discretisation options for clarity
        ptheta = tfield_options%ptheta
        beta = tfield_options%beta

        ! loop over elements
        do ele=1, element_count(tfield)
          x_ele=ele_val(x, ele)
          x_f=ele_val_at_quad(x, ele, x_cvshape)
          u_f=ele_val_at_quad(relu, ele, u_cvshape)
          nodes=>ele_nodes(tfield, ele)
          x_nodes=>ele_nodes(x_tfield, ele)
          if((tfield_options%upwind_scheme==CV_UPWINDVALUE_PROJECT_POINT).or.&
            (tfield_options%upwind_scheme==CV_UPWINDVALUE_PROJECT_GRAD)) then
            upwind_nodes=>x_nodes
          else
            upwind_nodes=>nodes
          end if

          ! find determinant and unorientated normal
          call transform_cvsurf_to_physical(x_ele, x_cvshape, &
                                            detwei, normal, cvfaces)

          cfl_ele = 1.0

          tfield_ele = ele_val(tfield_scomp, ele)
          oldtfield_ele = ele_val(oldtfield_scomp, ele)

          notvisited=.true.
          
          mat_local = 0.0
          rhs_local = 0.0

          ! loop over nodes within this element
          do iloc = 1, tfield%mesh%shape%loc

            ! loop over cv faces internal to this element
            do face = 1, cvfaces%faces

              ! is this a face neighbouring iloc?
              if(cvfaces%neiloc(iloc, face) /= 0) then
                oloc = cvfaces%neiloc(iloc, face)

                ! loop over gauss points on face
                do gi = 1, cvfaces%shape%ngi

                  ! global gauss pt index
                  ggi = (face-1)*cvfaces%shape%ngi + gi

                  ! have we been here before?
                  if(notvisited(ggi)) then
                    notvisited(ggi)=.false.

                    ! correct the orientation of the normal so it points away from iloc
                    normgi=orientate_cvsurf_normgi(node_val(x_tfield, x_nodes(iloc)),x_f(:,ggi),normal(:,ggi))

                    ! calculate u.n
                    udotn=dot_product(u_f(:,ggi), normgi(:))

                    inflow = (udotn<=0.0)

                    income = merge(1.0,0.0,inflow)

                    if((tfield_options%facevalue==CV_FACEVALUE_FIRSTORDERUPWIND).and.&
                       (.not.tfield_options%limit_theta))then
                      ! if we need the matrix then assemble it now
                      assert((i==1).and.(j==1).and.getmat)
                    
                      mat_local(iloc, oloc) = mat_local(iloc, oloc) &
                                            + detwei(ggi)*udotn*income
                      mat_local(oloc, iloc) = mat_local(oloc, iloc) &
                                            + detwei(ggi)*(-udotn)*(1.-income)
                      mat_local(iloc, iloc) = mat_local(iloc, iloc) &
                                            + detwei(ggi)*udotn*(1.0-income) &
                                            - (1.-beta)*detwei(ggi)*udotn
                      mat_local(oloc, oloc) = mat_local(oloc, oloc) &
                                            + detwei(ggi)*(-udotn)*income &
                                            - (1.-beta)*detwei(ggi)*(-udotn)

                    else
                      ! calculate the iterated pivot value (so far only does first order upwind)
                      ! which will be subtracted out from the rhs such that with an increasing number
                      ! of iterations the true implicit lhs pivot is cancelled out (if it converges!)
                      tfield_pivot_val = income*tfield_ele(oloc) + (1.-income)*tfield_ele(iloc)
                      
                      ! evaluate the nonlinear face value that will go into the rhs
                      ! this is the value that you choose the discretisation for and
                      ! that will become the dominant term once convergence is achieved
                      call evaluate_face_val(tfield_face_val, oldtfield_face_val, & 
                                            iloc, oloc, ggi, upwind_nodes, &
                                            t_cvshape, &
                                            tfield_ele, oldtfield_ele, &
                                            tfield_upwind, oldtfield_upwind, &
                                            inflow, cfl_ele, &
                                            tfield_options)

                      ! perform the time discretisation
                      tfield_theta_val=theta_val(iloc, oloc, &
                                          tfield_face_val, &
                                          oldtfield_face_val, &
                                          tfield_options%theta, dt, udotn, &
                                          x_ele, tfield_options%limit_theta, &
                                          tfield_ele, oldtfield_ele, &
                                          ftheta=ftheta)

                      rhs_local(i, j, iloc) = rhs_local(i, j, iloc) &
                                      + ptheta*udotn*detwei(ggi)*tfield_pivot_val &
                                      - udotn*detwei(ggi)*tfield_theta_val &
                                      + (1.-ftheta)*(1.-beta)*detwei(ggi)*udotn*oldtfield_ele(iloc)
                        
                      if(j/=i) then
                        rhs_local(j, i, iloc) = rhs_local(j, i, iloc) &
                                        + ptheta*udotn*detwei(ggi)*tfield_pivot_val &
                                        - udotn*detwei(ggi)*tfield_theta_val &
                                        + (1.-ftheta)*(1.-beta)*detwei(ggi)*udotn*oldtfield_ele(iloc)
                      end if
                        
                      rhs_local(i, j, oloc) = rhs_local(i, j, oloc) &
                                      + ptheta*(-udotn)*detwei(ggi)*tfield_pivot_val &
                                      - (-udotn)*detwei(ggi)*tfield_theta_val &
                                      + (1.-ftheta)*(1.-beta)*detwei(ggi)*(-udotn)*oldtfield_ele(oloc)
      
                      if(j/=i) then
                        rhs_local(j, i, oloc) = rhs_local(j, i, oloc) &
                                        + ptheta*(-udotn)*detwei(ggi)*tfield_pivot_val &
                                        - (-udotn)*detwei(ggi)*tfield_theta_val &
                                        + (1.-ftheta)*(1.-beta)*detwei(ggi)*(-udotn)*oldtfield_ele(oloc)
                      end if
                      
                      ! if we need the matrix then assemble it now
                      if(getmat .and. i == 1 .and. j == 1) then
                        mat_local(iloc, oloc) = mat_local(iloc, oloc) &
                                              + ptheta*detwei(ggi)*udotn*income
                        mat_local(oloc, iloc) = mat_local(oloc, iloc) &
                                              + ptheta*detwei(ggi)*(-udotn)*(1.-income)
                        mat_local(iloc, iloc) = mat_local(iloc, iloc) &
                                              + ptheta*detwei(ggi)*udotn*(1.0-income) &
                                              - ftheta*(1.-beta)*detwei(ggi)*udotn
                        mat_local(oloc, oloc) = mat_local(oloc, oloc) &
                                              + ptheta*detwei(ggi)*(-udotn)*income &
                                              - ftheta*(1.-beta)*detwei(ggi)*(-udotn)

                      end if
                      
                    end if

                  end if ! notvisited
                end do ! gi
              end if ! neiloc
            end do ! face
          end do ! iloc
          
          if(getmat.and.(i==1).and.(j==1)) then
            call addto(A_m, nodes, nodes, mat_local)
          end if
          
          if((tfield_options%facevalue/=CV_FACEVALUE_FIRSTORDERUPWIND).or.&
              (tfield_options%limit_theta))then
            call addto(rhs, nodes, rhs_local)
          end if
        end do ! ele
      end do ! j
    end do ! i

!     deallocate(detwei)
!     u_shape=>ele_shape(relu, 1)
!     t_shape=>ele_shape(tfield, 1)
!     allocate(detwei(ele_ngi(relu, 1)))
! 
!   ! Twisting terms.
!   ! Added in here separately to not confuse the CV formulation of the advection equation
!
!   ! Bound the velocity jacobian
!   U_nl_J => extract_tensor_field(state, "VelocityJacobian")
!   call allocate(bounded_U_nl_J, U_nl_J%mesh, "BoundedVelocityJacobian")
!   id = get_matrix_identity(U_nl%dim)
!   max_eigenbound => extract_tensor_field(state, "MaxMetricEigenbound")
!
!   do node=1,node_count(bounded_U_nl_J)
!     value = node_val(U_nl_J, node)
!
!     ! Divide into symmetric and asymmetric parts
!     strain = 0.5 * (value + transpose(value))
!     rotate = value - strain
!
!     frob_norm = frob(node_val(oldtfield, node) - node_val(max_eigenbound, node))
!     strain = anisotropic_min_abs(strain, id * (1.0 + frob_norm/100.0))
!     write(0,*) "factor == ", 1.0 + frob_norm/100.0
!     !strain = id
!     value = strain + rotate
!
!     call set(bounded_U_nl_J, node, value)
!   end do
!
!    do ele=1,ele_count(tfield)
!      ! Transform U_nl derivatives and weights into physical space.
!      call transform_to_physical(ele_val(X,ele), ele_shape(X,ele), u_shape , dm_t=du_t, detwei=detwei)
!      U_nl_J_q=ele_val_at_quad(bounded_U_nl_J, ele)
!      Twist_mat = shape_shape_tensor(T_shape, T_shape, detwei, U_nl_J_q)
!      do i=1,tfield%dim
!        do j=i,tfield%dim
!          twist_rhs=0.0
!          do k=1,tfield%dim
!             twist_rhs = twist_rhs + matmul(Twist_mat(i, k, :, :), &
!                  &      ele_val(oldtfield, k, j, ele) + (dt/rk_coeff(rk_it)) * ele_val(delta_tfield, k, j, ele))
!
!             twist_rhs = twist_rhs + matmul(Twist_mat(k, j, :, :), &
!                  &      ele_val(oldtfield, i, k, ele) + (dt/rk_coeff(rk_it)) * ele_val(delta_tfield, i, k, ele))
!          end do
!         
!          call addto(rhs, i, j, ele_nodes(tfield, ele), -1.0 * twist_rhs)
!          if (j /= i) then
!            call addto(rhs, j, i, ele_nodes(tfield, ele), -1.0 * twist_rhs)
!          end if
!
!        end do ! j
!      end do ! i
!    end do ! ele
!
!    call deallocate(bounded_U_nl_J)

    call deallocate(tfield_upwind)
    call deallocate(oldtfield_upwind)

    deallocate(x_ele, x_f, detwei, normal, normgi, u_f)
    deallocate(cfl_ele, tfield_ele, oldtfield_ele)
    deallocate(notvisited)

    contains

    function frob(R)
      real, dimension(:, :), intent(in) :: R
      real :: frob

      real, dimension(size(R, 1), size(R, 2)) :: X
      integer :: i
      real :: trace

      X = matmul(transpose(R), R)
      trace = 0.0
      do i=1,size(R, 1)
        trace = trace + X(i, i)
      end do
      frob = sqrt(trace)
    end function frob

    function anisotropic_min_abs(tensor1, tensor2) result(tensor3)
      real, dimension(:, :), intent(in) :: tensor1, tensor2
      real, dimension(size(tensor1, 1), size(tensor1, 1)) :: tensor3, F, T, Finv, evecs
      real, dimension(size(tensor1, 1)) :: evals
      integer :: i, dim
      real :: e_sign, e_abs

      dim = size(tensor1, 1)

      if (all(tensor2 == 0.0)) then
        call eigendecomposition_symmetric(tensor1, evecs, evals)
        do i=1,dim
          evals(i) = min(evals(i), 0.0)
        end do
        call eigenrecomposition(tensor3, evecs, evals)
        return
      end if

      ! So we are dealing with the non-degenerate case.

      call eigendecomposition_symmetric(tensor2, evecs, evals)
      F = get_deformation_matrix(tensor2, evecs, evals)
      Finv = inverse(F)

      T = transpose(Finv)
      tensor3 = matmul(matmul(T, tensor1), transpose(T))
      call eigendecomposition_symmetric(tensor3, evecs, evals)

      do i=1,dim
        e_sign = sign(1.0, evals(i))
        e_abs = abs(evals(i))
        evals(i) = min(e_abs, 1.0) * e_sign
      end do
      call eigenrecomposition(tensor3, evecs, evals)
      T = F
      tensor3 = matmul(matmul(transpose(T), tensor3), T)
    end function


  end subroutine assemble_advection_m_cv
  ! end of assembly subroutines
  !************************************************************************

  function rk_coeff(i)
  ! RK4 coefficients
    integer, intent(in) :: i
    real :: rk_coeff

    select case(i)
      case(1)
        rk_coeff = 1.0
      case(2)
        rk_coeff = 2.0
      case(3)
        rk_coeff = 2.0
      case(4)
        rk_coeff = 1.0
    end select

  end function

  subroutine insert_bounded_velocity_jacobian(state)
    type(state_type), intent(inout) :: state

    type(tensor_field) :: U_nl_J, rhs
    type(vector_field), pointer :: U_nl, X
    real, dimension(:), allocatable :: detwei
    real, dimension(:, :, :), allocatable :: U_nl_J_q, du_t, tensor_rhs
    real, dimension(:, :), allocatable :: little_mass_matrix
    integer :: ele, i, j, dim, loc, ngi
    type(scalar_field) :: lumped_mass
    type(element_type), pointer :: u_shape

    U_nl => extract_vector_field(state, "NonlinearVelocity")
    X => extract_vector_field(state, "Coordinate")
    u_shape => ele_shape(U_nl, 1)
    dim = U_nl%dim
    loc = ele_loc(U_nl, 1)
    ngi = ele_ngi(U_nl, 1)
    call allocate(U_nl_J, U_nl%mesh, "VelocityJacobian")
    call zero(U_nl_J)
    call allocate(rhs, U_nl%mesh, "VelocityJacobianRHS")
    call zero(rhs)
    call allocate(lumped_mass, U_nl%mesh, "LumpedMass")
    call zero(lumped_mass)

    allocate(detwei(ngi))
    allocate(U_nl_J_q(dim, dim, ngi))
    allocate(du_t(loc, ngi, dim))
    allocate(little_mass_matrix(loc, loc))
    allocate(tensor_rhs(dim, dim, loc))

    do ele=1,ele_count(U_nl)
      call transform_to_physical(X, ele, u_shape , dshape=du_t, detwei=detwei)
      little_mass_matrix = shape_shape(u_shape, u_shape, detwei)
      call addto(lumped_mass, ele_nodes(U_nl, ele), sum(little_mass_matrix, 2))

      U_nl_J_q = ele_jacobian_at_quad(U_nl, ele, du_t)
      tensor_rhs = shape_tensor_rhs(u_shape, U_nl_J_q, detwei)
      call addto(rhs, ele_nodes(U_nl, ele), tensor_rhs)
    end do

    if (U_nl%mesh%shape%degree /= 1) then
      ewrite(-1,*) "You need to write the code to do the full Galerkin projection here."
      ewrite(-1,*) "It's easy. Instead of lumping the mass, just form the full mass matrix"
      ewrite(-1,*) "and call petsc_solve."
      FLExit("bounded_velocity_jacobian not available for non P1 velocities")
    end if

    do i=1,dim
      do j=1,dim
        U_nl_J%val(i, j, :) = rhs%val(i, j, :) * (1.0 / lumped_mass%val)
      end do
    end do
    call deallocate(rhs)

    call insert(state, U_nl_J, "VelocityJacobian")
    call deallocate(U_nl_J)

    call deallocate(lumped_mass)
    
    deallocate(detwei)
    deallocate(U_nl_J_q)
    deallocate(du_t)
    deallocate(little_mass_matrix)
    deallocate(tensor_rhs)

  end subroutine insert_bounded_velocity_jacobian
end module metric_advection
