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

  module momentum_equation_reduced_adjoint
    use fields
    use state_module
    use spud
    use fldebug
    !      use momentum_cg
    use divergence_matrix_cv
    use divergence_matrix_cg
    use momentum_dg
    !      use assemble_cmc
    use field_priority_lists
    use momentum_diagnostic_fields, only: calculate_momentum_diagnostics
    use field_options
    !      use compressible_projection
    use boundary_conditions
    use boundary_conditions_from_options
    use sparse_matrices_fields
    use sparse_tools
    use sparse_tools_petsc
    !      use free_surface_module
    use solvers
    !      use full_projection
    !      use petsc_solve_state_module
    use Profiler
    !      use geostrophic_pressure
    !      use hydrostatic_pressure
    !      use vertical_balance_pressure
    !      use foam_drainage, only: calculate_drainage_source_absor
    !      use oceansurfaceforcing
    !      use drag_module
    use parallel_tools
    use linked_lists
    use sparsity_patterns_meshes
    use state_matrices_module
    use vtk_interfaces
    !      use rotated_boundary_conditions
    !      use Weak_BCs
    use reduced_model_runtime
    use state_fields_module
    !      use Tidal_module
    use Coordinates
    use diagnostic_fields, only: calculate_diagnostic_variable
    use dgtools, only: dg_apply_mass
    !      use slope_limiters_dg
    !      use implicit_solids
    use multiphase_module
    !      use pressure_dirichlet_bcs_cv
    use reduced_projection
    use sparse_tools_petsc 
    use global_parameters, only:theta
    
    
    implicit none

    private
    
    public :: solve_momentum_reduced_adjoint
    ! The timestep
    ! The timestep
    real :: dt
    
    ! Are we going to form the Diagonal Schur complement preconditioner?
    logical :: get_diag_schur
    ! Do we need the scaled pressure mass matrix?
    logical :: get_scaled_pressure_mass_matrix
    ! Do we need an auxiliary matrix for full_projection solve?
    logical :: assemble_schur_auxiliary_matrix
    
    ! Do we want to use the compressible projection method?
    logical :: use_compressible_projection
    ! Are we doing a full Schur solve?
    logical :: full_schur
    ! Are we lumping mass or assuming consistent mass?
    logical, dimension(:), allocatable :: lump_mass
    ! are we using a cv pressure 
    logical :: cv_pressure 
    ! for a CG pressure are we testing the continuity with cv
    logical :: cg_pressure_cv_test_continuity
    
    ! Do we need to reassemble the C^T or CMC matrices?
    logical :: reassemble_all_ct_m, reassemble_all_cmc_m
    
    ! Do we want to apply a theta weighting to the pressure gradient term?
    logical :: use_theta_pg
    
    ! Is a theta-weighting term present in the velocity-divergence?
    logical :: use_theta_divergence
    
    ! Are we using a discontinuous Galerkin discretisation?
    logical, dimension(:), allocatable :: dg
    ! True if advection-subcycling is performed
    logical, dimension(:), allocatable :: subcycle
    
    ! Apply KMK stabilisation?
    logical :: apply_kmk
    
    logical :: diagonal_big_m
    logical :: pressure_debugging_vtus
    
    ! Add viscous terms to inverse_masslump for low Re which is only used for pressure correction
    logical :: low_re_p_correction_fix
    
    ! Increased each call to momentum equation, used as index for pressure debugging vtus
    integer, save :: pdv_count = -1
    
    logical, dimension(:), allocatable :: sphere_absorption
    
    ! Are we running a multi-phase simulation?
    logical :: multiphase
    
    !! True if the momentum equation should be solved with the reduced model.
    logical :: reduced_model,DEIM
    integer :: ii,jj,deim_number,udim
    type(vector_field) :: deim_rhs_u
    type(scalar_field) :: deim_rhs_p
    ! real, dimension(:), allocatable :: pod_coef_deim
    type(state_type), dimension(:), pointer :: deim_state => null()   !output state from full model
    type(state_type), dimension(:), pointer :: deim_state_res => null()   !output state from reduced_model
    ! type(state_type), dimension(:), pointer :: deim_state_resl => null() 
    type(state_type), dimension(:), allocatable :: deim_state_resl
    

  contains
    
    subroutine solve_momentum_reduced_adjoint(state, at_first_timestep, timestep, POD_state, snapmean, eps, its)
      !!< Construct and solve the momentum and continuity equations
      !!< using Chorin's projection method (Chorin, 1968)
       
      ! An array of buckets full of fields
      ! The whole array is needed for the sake of multimaterial assembly
      type(state_type), dimension(:), intent(inout) :: state
     
      logical, intent(in) :: at_first_timestep
      integer, intent(in) :: timestep
      type(state_type), dimension(:,:,:), intent(inout) :: POD_state 
      
      ! Counter iterating over each state
      integer :: istate 
      
      ! The pressure projection matrix (extracted from state)
      type(csr_matrix), pointer :: cmc_m
      
      ! logical to indicate whether ct_m and cmc_m need reassembling
      ! (used for each state within the assembly loop)
      logical :: reassemble_ct_m, reassemble_cmc_m
      ! is there a pressure in state?
      logical :: have_pressure
      ! Are we solving a Poisson pressure equation?
      logical :: poisson_p
      
      ! Matrix sparsity patterns for the matrices we allocate locally
      type(csr_sparsity), pointer :: u_sparsity
      
      !! Locally allocated matrices:
      ! Matrix for split explicit advection
      type(petsc_csr_matrix), dimension(:), allocatable, target  :: big_m_tmp
      type(block_csr_matrix), dimension(:), allocatable :: subcycle_m
      ! Pointer to matrix for full projection solve:
      type(petsc_csr_matrix_pointer), dimension(:), allocatable :: inner_m
      ! Pointer to preconditioner matrix for full projection solve:
      type(csr_matrix), pointer :: full_projection_preconditioner
      ! Auxiliary matrix for full_projection solve
      type(csr_sparsity), pointer :: schur_auxiliary_matrix_sparsity
      type(csr_matrix) :: schur_auxiliary_matrix
      ! Scaled pressure mass matrix - used for preconditioning full projection solve:
      type(csr_matrix), target :: scaled_pressure_mass_matrix
      type(csr_sparsity), pointer :: scaled_pressure_mass_matrix_sparsity
      ! Left hand matrix of CMC. For incompressibe flow this points to ct_m as they are identical, 
      ! unless for CG pressure with CV tested continuity case when this matrix will be the 
      ! CV divergence tested matrix and ct_m the CG divergence tested matrix (right hand matrix of CMC).
      ! For compressible flow this differs to ct_m in that it will contain the variable density.
      type(block_csr_matrix_pointer), dimension(:), allocatable :: ctp_m
      ! The lumped mass matrix (may vary per component as absorption could be included)
      type(vector_field), dimension(1:size(state)) :: inverse_masslump, visc_inverse_masslump
      ! Mass matrix
      type(petsc_csr_matrix), dimension(1:size(state)), target :: mass
      ! For DG:
      type(block_csr_matrix), dimension(1:size(state)):: inverse_mass
       
      ! Momentum RHS
      type(vector_field), dimension(1:size(state)):: rhs_deim, rhs_advec,rhs_deim_res
      
      ! Projection RHS
      type(scalar_field) :: projec_rhs
      
      ! Do we want to assemble the KMK stabilisation matrix?
      logical :: assemble_kmk
      
      ! Change in pressure
      type(scalar_field) :: delta_p
      ! Change in velocity
      type(vector_field) :: delta_u
      
      ! Dummy fields
      type(scalar_field), pointer :: dummyscalar, dummydensity, dummypressure
      
      ! Pressure and density
      type(scalar_field), pointer :: p, density
      type(mesh_type), pointer :: p_mesh
      ! Velocity and space
      type(vector_field), pointer :: u, x
      
      ! with free-surface or compressible pressure projection pressures 
      ! are at integer time levels and we apply a theta weighting to the
      ! pressure gradient term
      real :: theta_pg
      ! in this case p_theta=theta_pg*p+(1-theta_pg)*old_p
      type(scalar_field), pointer :: old_p, p_theta
      ! With free-surface or compressible-projection the velocity divergence is
      ! calculated at time n+theta_divergence instead of at the end of the timestep
      real :: theta_divergence
      type(vector_field), pointer :: old_u
      ! all of this only applies if use_theta_pg .eqv. .true.
      ! without a free surface, or with a free surface and theta==1
      ! use_theta_pg .eqv. .false. and p_theta => p
      
      ! What is the equation type?
      character(len=FIELD_NAME_LEN) :: equation_type, poisson_scheme, schur_scheme, pressure_pmat
      
      integer :: stat
      real ::  theta_pp
      
      ! The list of stiff nodes
      ! This is saved because the list is only formed when cmc is assembled, which
      ! isn't necessarily every time this subroutine is called but the list is
      ! still needed to fix the rhs (applying the fix to cmc itself wipes out the
      ! information that would be required to recompile the list)
      type(ilist), save :: stiff_nodes_list
      
      
      !! Variables for multi-phase flow model
      integer :: prognostic_count
      ! Do we have a prognostic pressure field to solve for?
      logical :: prognostic_p = .false.
      ! Prognostic pressure field's state index (if present)
      integer :: prognostic_p_istate 
      ! The 'global' CMC matrix (the sum of all individual phase CMC matrices)
      type(csr_matrix), pointer :: cmc_global
      ! An array of submaterials of the current phase in state(istate).
      type(state_type), dimension(:), pointer :: submaterials
      ! The index of the current phase (i.e. state(istate)) in the submaterials array
      integer :: submaterials_istate
      ! Do we have fluid-particle drag between phases?
      logical :: have_fp_drag
      real :: scale
      
      
      !!for reduced model
      type(vector_field), pointer :: snapmean_velocity
      type(scalar_field), pointer :: snapmean_pressure
      type(vector_field), pointer :: POD_velocity, POD_velocity_deim, velocity_deim,velocity_deim_snapmean
      type(scalar_field), pointer :: POD_pressure
      
      type(pod_matrix_type) :: pod_matrix, pod_matrix_mass, pod_matrix_adv,pod_matrix_B
      type(pod_rhs_type) :: pod_rhs
      real, dimension(:), allocatable :: pod_coef,pod_coef_dt
      real, dimension(:,:), allocatable :: pod_sol_velocity, pod_ct_m
      real, dimension(:), allocatable :: pod_sol_pressure
      integer :: d, i, j, k,dim,jd
      real, intent(in) :: eps
      logical, intent(in) :: snapmean
      logical :: timestep_check
      type(scalar_field) :: u_cpt
      
      !  real, dimension(:,:), allocatable :: A_deim ,mom_rhs_deim  !() name change
      !  real, dimension(:), allocatable :: b_deim,Ny
      real, dimension(:,:), allocatable :: P_mn,Temp_pu,Temp_vu,Temp_vupu,mom_rhs_deim !A_deim,
      real, dimension(:,:,:), allocatable :: V_kn,U_nm, U_mn
      real, dimension(:), allocatable :: Ny ,b_deim
      real, dimension(:,:), allocatable :: A_deim    !() name change
      integer :: unodes,d1, AI,AJ,AM,AN,d2   
      
      !petro
      type(scalar_field), pointer :: POD_u_scalar  !petro
      type(pod_rhs_type)::pod_rhs_old
      type(vector_field), pointer :: POD_u
      real, dimension(:), allocatable :: theta_pet,smean_gi,KB_pod,b_pod,PSI_GI,PSI_OLD_GI,psi_old,GRADX,GRADT,DIFFGI_pod,PSI
      real, dimension(:,:), allocatable ::leftsvd,leftsvd_gi,leftsvd_x_gi, AMAT_pod, KMAT_pod
      
      !real, dimension(:,:), allocatable :: snapmean_velocity
      type(element_type), pointer :: x_shape,xf_shape
      real, dimension(:,:),allocatable ::X_val !for getting the result of detwei
      real, dimension(:), allocatable :: detwei(:),res(:)
      real ,dimension(:,:),allocatable :: JJ   !(X%dim, mesh_dim(X)),
      real :: det ! for getting the result of detwei
      real, dimension(:,:,:), allocatable ::  dshape
      
      ! real, dimension(POD_u%dim,ele_loc(POD_u,ele)) :: X_val !for getting the result of detwei
      ! new added for petro-galerkin
      integer :: TOTELE,NLOC,NGI,nsvd !,stat
      integer ::  POD_num
      type(mesh_type), pointer :: pod_umesh
      REAL PSI_THETA
      INTEGER ELE,globi,globj,isvd,jsvd,jloc,GI 
      REAL A_STAR_COEF,AX_STAR,P_STAR_POD
      REAL DT1
      real noloc,nogas
      logical :: petrov
      real, dimension(:,:), allocatable :: pod_matrix_snapmean,  pod_matrix_adv_mean
      real, dimension(:), allocatable :: pod_rhs_snapmean
      real, dimension(:,:), allocatable :: pod_matrix_perturbed,pod_matrix_adv_perturbed
      real, dimension(:), allocatable :: pod_rhs_perturbed,pod_ct_rhs,pod_rhs_adv_perturbed
      
      !nonlinear_iteration_loop
      integer :: nonlinear_iterations
      integer, intent(in) :: its
      
      !free surface matrix
      type(csr_matrix) :: fs_m
      logical :: on_sphere, have_absorption, have_vertical_stabilization
      type(vector_field), pointer :: dummy_absorption
      logical :: lump_mass_form
      real :: finish_time,current_time
      integer :: total_timestep
      ! Adjoint model
      logical adjoint_reduced
      type(pod_rhs_type) :: pod_rhs_adjoint
      type(pod_matrix_type) :: adjoint_pod_A,adjoint_pod_B,adjoint_A_extra,adjoint_A
      real, dimension(:,:), allocatable :: pod_coef_all,pod_coef_adjoint ! solution of adjoint

      ! gradient:
      real, dimension(:), allocatable :: g,ds,ds_tmp
      
      !  else !do adjoint_reduced
      call get_option("/timestepping/current_time", current_time)
      call get_option("/timestepping/finish_time", finish_time)       
      call get_option("/timestepping/timestep", dt)
      call get_option(trim(u%option_path)//"/prognostic/temporal_discretisation/theta", &
           theta)
      total_timestep=(finish_time-current_time)/dt
      call allocate(adjoint_pod_A, POD_velocity, POD_pressure)
      call allocate(adjoint_A_extra, POD_velocity, POD_pressure) 
      call allocate(adjoint_A, POD_velocity, POD_pressure)
      call allocate(adjoint_pod_B, POD_velocity, POD_pressure)
      allocate(pod_coef_all(total_timestep,((u%dim+1)*size(POD_state,1))))
      allocate(pod_coef_adjoint(total_timestep,((u%dim+1)*size(POD_state,1))))
      read(100,*)((pod_coef_all(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=0,total_timestep)
      open(unit=20,file='pod_matrix_snapmean')
      read(20,*)((pod_matrix_snapmean(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
      close(20)
      
      if(timestep.eq.1) then  !!! do we run the adjoint from timestep = 1 or total_timestep??
         allocate(g((u%dim+1)*size(POD_state,1))) !!  where we deallocate it?
      endif
      allocate(ds((u%dim+1)*size(POD_state,1)))
      allocate(ds_tmp((u%dim+1)*size(POD_state,1)))

      g = 0.0
      ds = 0.0
      ds_tmp = 0.0

      pod_matrix%val=pod_matrix_snapmean(:,:)
      pod_matrix_B%val=0.0
      
      open(1,file='coef_pod_all_obv')
      read(1,*)((pod_coef_adjoint(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=0,total_timestep)
      ! pod_coef_adjoint() save the coef_pod_all_obv temporarily
      close(1)
      pod_coef_adjoint(:,:)=pod_coef_adjoint(:,:)-pod_coef_all(:,:)
      open(30,file='advection_matrix_perturbed')
      adjoint_pod_A%val =0.0
      adjoint_pod_B%val =0.0
      do k=1,size(POD_state,1)
         do d=1, u%dim
            !adjoint_pod_A%val size same as pod_matrix. read pod_matrix_perturbed here including theta*dt
            read(30,*)((pod_matrix_perturbed(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
            ! delete theta*dt from pod_matrix_perturbed
            pod_matrix_perturbed = pod_matrix_perturbed/(theta*dt)
            pod_matrix%val=pod_matrix%val+  &
                 theta*pod_coef_all(timestep,k+(d-1)*size(POD_state,1))*pod_matrix_perturbed(:,:)
            pod_matrix_B%val = pod_matrix_B%val + &
                 (1.-theta)*pod_coef_all(timestep,k+(d-1)*size(POD_state,1))*pod_matrix_perturbed(:,:)
            adjoint_pod_A%val((d-1)*size(POD_state,1)+k,:)= &
                 theta*matmul(pod_matrix_perturbed(:,:),pod_coef_all(timestep+1, :))
            ! theta*matmul(pod_matrix_perturbed(:,:),pod_coef_all(timestep+1, k+(d-1)*size(POD_state,1)))
            adjoint_pod_B%val((d-1)*size(POD_state,1)+k,:)= &
                 (1.-theta)*matmul(pod_matrix_perturbed(:,:),pod_coef_all(timestep, :))
            !(1.-theta)*matmul(pod_matrix_perturbed(:,:),pod_coef_all(timestep, k+(d-1)*size(POD_state,1)))
            !   pod_matrix%val=pod_matrix%val+pod_coef_all(timestep, k+(d-1)*size(POD_state,1))*pod_matrix_perturbed(:,:)
            
         enddo
         
         read(30,*)((pod_matrix_perturbed(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))

         pod_matrix%val=pod_matrix%val+theta*pod_coef_all(timestep,k+u%dim*size(POD_state,1))*pod_matrix_perturbed(:,:)
         pod_matrix_B%val = pod_matrix_B%val + (1.-theta)*pod_coef_all(timestep,k+u%dim*size(POD_state,1))*pod_matrix_perturbed(:,:)

         adjoint_pod_A%val((d-1)*size(POD_state,1)+k,:)=&
              theta*matmul(pod_matrix_perturbed(:,:),pod_coef_all(timestep+1, :))
         !theta*matmul(pod_matrix_perturbed(:,:),pod_coef_all(timestep+1, k+(d-1)*size(POD_state,1)))
         adjoint_pod_B%val((d-1)*size(POD_state,1)+k,:)=&
              (1.-theta)*matmul(pod_matrix_perturbed(:,:),pod_coef_all(timestep, :))
         !(1.-theta)*matmul(pod_matrix_perturbed(:,:),pod_coef_all(timestep, k+(d-1)*size(POD_state,1)))
         !   pod_matrix%val=pod_matrix%val+pod_coef_all(timestep, k+(d-1)*size(POD_state,1))*pod_matrix_perturbed(:,:)
      enddo
      
      close(30)
      
      if(timestep.eq.total_timestep) then
         adjoint_A_extra%val=adjoint_pod_A%val      
         adjoint_A_extra%val=adjoint_pod_A%val + pod_matrix_B%val
         adjoint_A_extra%val=transpose(adjoint_A_extra%val)
         adjoint_A%val=transpose(pod_matrix%val)     
         pod_rhs_adjoint%val=matmul(adjoint_A_extra%val,pod_rhs_adjoint%val)
      else
         adjoint_A_extra%val=adjoint_pod_A%val+ adjoint_pod_B%val        
         adjoint_A_extra%val=adjoint_pod_A%val+ adjoint_pod_B%val + pod_matrix_B%val
         adjoint_A_extra%val=transpose(adjoint_A_extra%val)
         adjoint_A%val=transpose(pod_matrix%val)     
         pod_rhs_adjoint%val=matmul(adjoint_A_extra%val,pod_rhs_adjoint%val)    
      endif
      call solve(adjoint_A%val,pod_rhs_adjoint%val)
      pod_coef_adjoint(timestep,:)=pod_rhs_adjoint%val(:)
      if(timestep.eq.total_timestep) then
         ds(:) = matmul(adjoint_pod_B%val, pod_coef_all(1, :) )  !!!Is pod_coef_all(1, :) the initial one??
         g = ds
         do i = 1,size(g)
            ds = 0.0
            ds(i) = 1.0
            ds_tmp = matmul(adjoint_pod_B%val,ds)
            g(i)= g(i) + dot_product(pod_rhs_adjoint%val,ds_tmp)
         enddo
      endif
      
      call deallocate(adjoint_pod_A)
      call deallocate(adjoint_A_extra)
      call deallocate(adjoint_A)
      call deallocate(adjoint_pod_B)
      deallocate(pod_coef_all)
      deallocate(pod_coef_adjoint)
      deallocate(ds)
      deallocate(ds_tmp)
      ! endif !end if(.not.adjoint_reduced)
      
    end subroutine solve_momentum_reduced_adjoint
    
  end module momentum_equation_reduced_adjoint
