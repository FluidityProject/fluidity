#include "fdebug.h"

   module momentum_equation_reduced

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
      use global_parameters, only: option_path_len,current_time
      use write_state_module
     ! use steam_nrs_module
      use Bp_Neutral_Net 
      use ann_bp
      implicit none

      private
      public :: solve_momentum_reduced, deim_state, deim_state_res,deim_state_resl,deim_number,pod_coef_all_o

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
      character(len = OPTION_PATH_LEN) :: simulation_name
      real :: finish_time
      integer :: pod_coef_size
      real, dimension(:,:), allocatable ::pod_coef_all_o,pod_coef_adjoint,pod_coef_all_temp
     ! real, dimension(:,:), allocatable ::pod_coef_all_diff, pod_coef_all_ann
      real :: coef,coefdiff,output
   contains

     subroutine solve_momentum_reduced(state, u,p,big_m, ct_m, mom_rhs, ct_rhs, inverse_masslump, &
          at_first_timestep, timestep, POD_state, POD_state_deim, snapmean, eps, its,total_timestep,if_optimal,fs_m)
       
       !!< Construct and solve the momentum and continuity equations
       !!< using Chorin's projection method (Chorin, 1968)
       
       ! An array of buckets full of fields
       ! The whole array is needed for the sake of multimaterial assembly
       type(state_type), dimension(:), intent(inout) :: state
       ! Momentum LHS
       type(petsc_csr_matrix), dimension(:), intent(inout), target :: big_m
       ! The pressure gradient matrix (extracted from state)
       type(block_csr_matrix_pointer), dimension(:), intent(inout) :: ct_m
       type(vector_field), dimension(1:size(state)), intent(inout) :: mom_rhs
       type(scalar_field), dimension(1:size(state)), intent(inout) :: ct_rhs
       logical, intent(in) :: at_first_timestep
       integer, intent(in) :: timestep
       type(state_type), dimension(:,:,:), intent(inout) :: POD_state
       type(state_type), dimension(:), intent(inout) :: POD_state_deim
       logical, optional, intent(in) :: if_optimal

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
       
       type(pod_matrix_type) :: pod_matrix, pod_matrix_mass, pod_matrix_adv 
       type(pod_rhs_type) :: pod_rhs
       real, dimension(:), allocatable :: pod_coef,pod_coef_dt,pod_coef0
       real, dimension(:,:), allocatable :: pod_sol_velocity, pod_ct_m
       real, dimension(:), allocatable :: pod_sol_pressure
       integer :: d, i, j, k,dim,jd
       real, intent(in) :: eps
       logical, intent(in) :: snapmean
       logical :: timestep_check
       type(scalar_field) :: u_cpt
       
       !real, dimension(:,:), allocatable :: A_deim ,mom_rhs_deim  !() name change
       !real, dimension(:), allocatable :: b_deim,Ny
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
       integer, intent(in) :: its, total_timestep
       
       !free surface matrix
       !type(csr_matrix) :: fs_m
       type(csr_matrix), optional, intent(inout) :: fs_m
       logical :: on_sphere, have_absorption, have_vertical_stabilization
       type(vector_field), pointer :: dummy_absorption
       logical :: lump_mass_form
       
       ! Adjoint model
       logical :: reduced_adjoint
       type(pod_rhs_type) :: pod_rhs_adjoint
       integer :: timestep_tmp, dump_no_tmp

       ! optimal_observation
       logical :: adaptive_observation
       integer :: no_optimal_obv
       real,dimension(:), allocatable :: u_tmp,optimal_location
       integer :: ctrl_timestep

       dump_no_tmp = 0

       ewrite(1,*) 'Entering solve_momentum_reduced'
       ewrite(1,*) '*********************************'
       call get_option('/timestepping/nonlinear_iterations', nonlinear_iterations, default=1)
       !!print*,'its = ', its, ' nonlinear_iterations = ', nonlinear_iterations
       nonlinear_iterations=1

       u => extract_vector_field(state(1), "Velocity", stat)
       p => extract_scalar_field(state(1), "Pressure")

        
       !enter reduced model process after the last nonliear_iteration
       ! if(its==nonlinear_iterations)then
       call profiler_tic("reduced_model_loop")
       call get_option(&
            '/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
       deim=have_option("/reduced_model/discrete_empirical_interpolation_method")
       reduced_adjoint= have_option("/reduced_model/adjoint")

       petrov=have_option("/reduced_model/execute_reduced_model/petrov_galerkin")  
       call get_option('/reduced_model/pod_basis_formation/pod_basis_count', deim_number)
       
       reduced_model_loop: do istate = 1, size(state)
          ! Get the velocity


          lump_mass_form = have_option(trim(u%option_path)//&
               &"/prognostic/spatial_discretisation"//&
               &"/continuous_galerkin/mass_terms/lump_mass_matrix")
          call get_option(trim(u%option_path)//"/prognostic/temporal_discretisation/theta", &
               theta)
          ! If there's no velocity then cycle
          if(stat/=0) cycle
          ! If this is an aliased velocity then cycle
          if(aliased(u)) cycle
          ! If the velocity isn't prognostic then cycle
          if(.not.have_option(trim(u%option_path)//"/prognostic")) cycle
          
          if(have_option(trim(u%option_path)//"/prognostic/spatial_discretisation"//&
               &"/continuous_galerkin").or.&
               have_option(trim(u%option_path)//"/prognostic/spatial_discretisation"//&
               &"/discontinuous_galerkin")) then
             
             !!speedup version of reduced model
             
             ! Allocate the change in pressure field
             call allocate(delta_p, p%mesh, "DeltaP")
             delta_p%option_path = trim(p%option_path)
             call zero(delta_p)
             
             call allocate(delta_u, u%dim, u%mesh, "DeltaU")
             delta_u%option_path = trim(u%option_path)
             call zero(delta_u)
             
            
             call get_option("/timestepping/timestep", dt)
             
             POD_velocity=>extract_vector_field(POD_state(1,1,istate), "Velocity")
             POD_pressure=>extract_scalar_field(POD_state(1,2,istate), "Pressure")    
             call allocate(pod_matrix, POD_velocity, POD_pressure)
             call allocate(pod_matrix_mass, POD_velocity, POD_pressure)
             call allocate(pod_matrix_adv, POD_velocity, POD_pressure)             
             call allocate(pod_rhs, POD_velocity, POD_pressure)
             
             allocate(pod_matrix_snapmean((u%dim+1)*size(POD_state,1),(u%dim+1)*size(POD_state,1)))
             allocate(pod_matrix_adv_mean((u%dim+1)*size(POD_state,1),(u%dim+1)*size(POD_state,1)))
             allocate(pod_rhs_snapmean((u%dim+1)*size(POD_state,1)))
             allocate(pod_matrix_perturbed((u%dim+1)*size(POD_state,1),(u%dim+1)*size(POD_state,1)))
             allocate(pod_matrix_adv_perturbed((u%dim+1)*size(POD_state,1),(u%dim+1)*size(POD_state,1)))
             allocate(pod_rhs_perturbed((u%dim+1)*size(POD_state,1)))
             allocate(pod_rhs_adv_perturbed((u%dim+1)*size(POD_state,1)))
             allocate(pod_coef((u%dim+1)*size(POD_state,1)))
             allocate(pod_coef_dt((u%dim+1)*size(POD_state,1)))
             allocate(pod_ct_m(size(POD_state,1),u%dim*size(POD_state,1)))
             if(timestep.eq.total_timestep) then
                allocate(pod_coef0((u%dim+1)*size(POD_state,1)))
             endif

             pod_coef_size=(u%dim+1)*size(POD_state,1)

             if(have_option("/reduced_model/adjoint")) then
  		!========================================================================
		! ADJOINT MODEL
		!========================================================================
                call get_option("/timestepping/current_time", current_time)
                call get_option("/timestepping/finish_time", finish_time)       
                call get_option("/timestepping/timestep", dt)
!                total_timestep=int((finish_time-current_time)/dt)+1
                allocate(pod_coef_adjoint(0:total_timestep,(u%dim+1)*size(POD_state,1)))
                pod_coef_adjoint =0.0
                call  solve_reduced_adjoint(state, (total_timestep-timestep), POD_state,total_timestep,pod_coef_adjoint)
           
                deallocate(pod_coef_adjoint)
             else
                if(present(if_optimal)) then
                   timestep_tmp = timestep+1
                else
                   if(timestep==1.and. eps.eq.0.0 .and. .not.snapmean)then
                      !get the initial velocity
                      u=>extract_vector_field(state(istate), "Velocity")
                      !get the initial pressure
                      p=>extract_scalar_field(state(istate), "Pressure")

                      adaptive_observation = .false.
                      if(adaptive_observation) then
                         if(.false.) then
                            open(10,file = 'observation_u_ini.dat')
                            write(10,*) (u%val(1,i),i=1,size(u%val,2))
                            close(10)
                            open(10,file = 'observation_v_ini.dat')
                            write(10,*) (u%val(2,i),i=1,size(u%val,2))
                            close(10)
                            open(10,file = 'observation_p_ini.dat')
                            write(10,*) (p%val(i),i=1,size(u%val,2))
                            close(10)
!                            stop 24
                         else
                            no_optimal_obv = 100
                            allocate(u_tmp(size(u%val,2)))
                            allocate(optimal_location(size(u%val,2)))
                            open(10,file='optimal_sensor_location.dat')
                            read(10,*) (optimal_location(i),i=1,size(optimal_location))
                            close(10)  
                            ! u component
                            open(10,file = 'observation_u_ini.dat')
                            read(10,*) (u_tmp(i),i=1,size(u%val,2))
                            close(10)
                            do i=1,no_optimal_obv
                               j=optimal_location(i)
                               print*,u%val(1,j),u_tmp(j)
                      !         u%val(1,j)=u_tmp(j)
                               u%val(1,size(optimal_location)-i)=u_tmp(size(optimal_location)-i)
                            enddo
                            ! v component
                            open(10,file = 'observation_v_ini.dat')
                            read(10,*) (u_tmp(i),i=1,size(u%val,2))
                            close(10)
                            do i=1,no_optimal_obv
                               j=optimal_location(i)
                               print*,u%val(2,j),u_tmp(j)
                             !  u%val(2,j)=u_tmp(j)
                               u%val(2,size(optimal_location)-i)=u_tmp(size(optimal_location)-i)
                            enddo
                           ! p 
                            open(10,file = 'observation_p_ini.dat')
                            read(10,*) (u_tmp(i),i=1,size(p%val))
                            close(10)
                            do i=1,no_optimal_obv
                               j=optimal_location(i)
                     !          p%val(j)=u_tmp(j)
                               p%val(size(optimal_location)-i)=u_tmp(size(optimal_location)-i)
                            enddo
                            deallocate(u_tmp)
                            deallocate(optimal_location)
                         endif
                      endif
                   !   stop 66
                      call project_from_full_to_pod(istate,  pod_state, state, pod_coef)
                      ! save the initial coeficient
                      open(101,file='coef_pod_all')
                      write(101,*)(pod_coef(i),i=1,(u%dim+1)*size(POD_state,1))
                      close(101)
                      open(101,file='coef_pod_all0')
                      write(101,*)(pod_coef(i),i=1,(u%dim+1)*size(POD_state,1))
                      close(101)
                      timestep_tmp = timestep+1
                   else
                      timestep_tmp = timestep
                   endif
                endif
                !========================================================================
		! FORWARD MODEL
		!========================================================================
                if (timestep_tmp==1)then
                   ! add -ct_m^T*p to mom_rhs
                   !----------------------
                   !                     call advance_velocity(state, istate, x, u, p_theta, big_m, ct_m, &
                   !                                        mom_rhs, subcycle_m, inverse_mass)
                   
                   !!the ct_rhs for reduced model
                   !                     ct_rhs(istate)%val=ct_rhs(istate)%val/dt/theta_divergence
                   !                ct_rhs(istate)%val=ct_rhs(istate)%val/theta_divergence
                   ct_rhs(istate)%val=ct_rhs(istate)%val
                   
                   
                   print*,has_boundary_condition(u,'free_surface')
                   !project to reduced matrix
                   if(has_boundary_condition(u,'free_surface'))then
                    !  ct_rhs(istate)%val=0
                      call project_reduced(big_m(istate), mom_rhs(istate), ct_m(istate)%ptr, ct_rhs(istate), pod_matrix, &
                           pod_rhs, dt, POD_state(:,:,istate), fs_m=fs_m)
                   else
                      call project_reduced(big_m(istate), mom_rhs(istate), ct_m(istate)%ptr, ct_rhs(istate), pod_matrix, &
                           pod_rhs, dt, POD_state(:,:,istate))
                   endif
                   !save the pod_matrix and pod_rhs for snapmean state
                   if(snapmean)then
                      ewrite(1,*)'snapmean'
                      open(unit=20,file='pod_matrix_snapmean')
                      write(20,*)((pod_matrix%val(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
                      write(20,*)(pod_rhs%val(i),i=1,(u%dim+1)*size(POD_state,1))
                      close(20)
                      allocate(dg(size(state)))
                      dg(istate) = have_option(trim(u%option_path)//&
                           &"/prognostic/spatial_discretisation"//&
                           &"/discontinuous_galerkin")
                      
                      allocate(big_m_tmp(size(state)))
                      
                      if(dg(istate)) then
                         call allocate_big_m_dg(state(istate), big_m_tmp(istate), u)
                      else
                         ! Create a sparsity if necessary or pull it from state:
                         u_sparsity => get_csr_sparsity_firstorder(state, u%mesh, u%mesh)
                         !diagonal_big_m = .true.

                         call allocate(big_m_tmp(istate), u_sparsity, (/u%dim, u%dim/), &
                              diagonal=diagonal_big_m, name="BIG_m")
                      end if
                      call zero(big_m_tmp(istate))
                      pod_rhs%val=0.0
                      !  pod_matrix%val=0.0
                      pod_matrix_mass%val=0.0
                      pod_matrix_adv%val=0.0
                      if(lump_mass_form) then
                         do dim = 1, u%dim !inverse_masslump(istate)%dim
                            do i = 1,node_count(u)
                               !! inverse_masslump = masslump in reduced order modelling 
                               !! we comment invert(inverse_masslump) in Momentum_CG
                               !! here, big_m_tmp is tmp matrix
                               call addto_diag(big_m_tmp(istate), dim, dim, i, inverse_masslump(istate)%val(dim,i) )
                                                      
                            enddo
                         enddo
                     
                          if(has_boundary_condition(u,'free_surface'))then
                          !print *,'ct_rhs(istate)%val', ct_rhs(istate)%val
                          !ct_rhs(istate)%val=0
                          call project_reduced(big_m_tmp(istate), mom_rhs(istate), ct_m(istate)%ptr, ct_rhs(istate), pod_matrix_mass, &
                               pod_rhs, dt, POD_state(:,:,istate),fs_m=fs_m)
                        else
                          call project_reduced(big_m_tmp(istate), mom_rhs(istate), ct_m(istate)%ptr, ct_rhs(istate), pod_matrix_mass, &
                               pod_rhs, dt, POD_state(:,:,istate))
                         endif
                      
                            
                         
                      else
                         ! print*,'mass'
                         !call project_reduced(mass(istate), mom_rhs(istate), ct_m(istate)%ptr, ct_rhs(istate), pod_matrix, &
                         !   pod_ct_m, pod_rhs, 1.0, POD_state(:,:,istate))
                      endif

                           pod_matrix_mass%val(u%dim*size(POD_state,1)+1:(u%dim+1)*size(POD_state,1),:) =0.0
                           pod_matrix_mass%val(1:u%dim*size(POD_state,1),u%dim*size(POD_state,1)+1:(u%dim+1)*size(POD_state,1)) =0.0
                      
                      !save the mass_matrix 
                      open(unit=20,file='mass_matrix')
                      write(20,*)((pod_matrix_mass%val(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
                      close(20)

                      
                      !save the advection_matrix 
                      pod_matrix_adv%val(:,:)=pod_matrix%val(:,:)-pod_matrix_mass%val(:,:)
                      !! clear the rest part in advection_matrix
      !                pod_matrix_adv%val(u%dim*size(POD_state,1)+1:(u%dim+1)*size(POD_state,1),:) =0.0
      !                pod_matrix_adv%val(1:u%dim*size(POD_state,1), u%dim*size(POD_state,1)+1:(u%dim+1)*size(POD_state,1)) = 0.0
                      open(unit=21,file='advection_matrix_snapmean')
                      write(21,*)((pod_matrix_adv%val(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
                      close(21)
                      !  open(unit=21,file='big_m')
                      !  write(21,*)((pod_matrix%val(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
                      !   close(21)
                      
                      !  open(unit=20,file='ct_m_matrix')
                      ! write(20,*) ((pod_ct_m(i,j),j=1,u%dim*size(POD_state,1)),i=1,size(POD_state,1))
                      !  close(20)
                      deallocate(big_m_tmp)
                      deallocate(dg)
                      
                   endif


                
                   if(eps.ne.0.0)then
                      pod_matrix_mass%val=0.0
                      pod_matrix_adv%val=0.0
                      pod_matrix_adv_perturbed=0.0
                      open(unit=20,file='mass_matrix')
                      read(20,*)((pod_matrix_mass%val(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
                      close(20)
                      pod_matrix_adv%val(:,:)=pod_matrix%val(:,:)-pod_matrix_mass%val(:,:)
                      !! clear the rest part in advection_matrix
       !               pod_matrix_adv%val(u%dim*size(POD_state,1)+1:(u%dim+1)*size(POD_state,1),:) =0.0
       !               pod_matrix_adv%val(1:u%dim*size(POD_state,1), u%dim*size(POD_state,1)+1:(u%dim+1)*size(POD_state,1)) = 0.0
                      
                      open(unit=20,file='advection_matrix_snapmean')
                      read(20,*)((pod_matrix_adv_mean(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
                      close(20)
                      pod_matrix_adv_perturbed(:,:)=(pod_matrix_adv%val(:,:)-pod_matrix_adv_mean(:,:))/eps
                      !save the advection_matrix 
                      ! open(unit=21,file='advection_matrix_perturbed')
                      write(60,*)((pod_matrix_adv_perturbed(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
                      !  close(21)
                      
                      ewrite(1,*)'perturbed'
                      ! print*,'perturbed'
                      open(unit=20,file='pod_matrix_snapmean')
                      read(20,*)((pod_matrix_snapmean(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
                      read(20,*)(pod_rhs_snapmean(i),i=1,(u%dim+1)*size(POD_state,1))
                      close(20)
                      pod_matrix%val(:,:)=(pod_matrix%val(:,:)-pod_matrix_snapmean(:,:))/eps
                      pod_rhs%val(:)=(pod_rhs%val(:)-pod_rhs_snapmean(:))/eps
                      
                      
                      call project_from_full_to_pod(istate,  pod_state, state, pod_coef)

                      call Matrix_vector_multiplication(size(pod_coef),pod_rhs_adv_perturbed, &
                           pod_matrix_adv_perturbed,pod_coef)
                      
                      pod_rhs%val = pod_rhs%val + pod_rhs_adv_perturbed/(theta*dt) ! delete advection term in rhs
                   
                      !save pod_matrix and pod_rhs for perturbed state (related to POD_velocity%dim, size(POD_state,1))
                      write(30,*)((pod_matrix%val(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
                      write(50,*)(pod_rhs%val(i),i=1,(u%dim+1)*size(POD_state,1))! delete advection term in rhs
                      
                   elseif(eps.eq.0.0 .and. .not.snapmean)then
                      !need to exclude snapmean
                      ewrite(1,*)'reduced'
                      call solve(pod_matrix%val, pod_rhs%val)
                      pod_coef_dt(:)=pod_rhs%val
                      pod_rhs%val=0.0
                      
                      !get the initial velocity
                      u=>extract_vector_field(state(istate), "Velocity")
                      !get the initial pressure
                      p=>extract_scalar_field(state(istate), "Pressure")

                      call project_from_full_to_pod(istate,  pod_state, state, pod_coef)
                      ! save the initial coeficient
                      open(101,file='coef_pod_all')
                      write(101,*)(pod_coef(i),i=1,(u%dim+1)*size(POD_state,1))
                      close(101)
                      
                      !!print*,'before add pod_coef_dt*dt to pod_coef_initial'
                      !!print*,pod_coef(2),pod_coef_dt(2)
                      !divided by the number of nonlinear_iterations???
                      ! for  velocity, here we solved du/dt
                      pod_coef(1:size(POD_state,1)*u%dim)=pod_coef(1:size(POD_state,1)*u%dim)+dt*pod_coef_dt(1:size(POD_state,1)*u%dim)
                      ! for pressure, here we solve dp
                      pod_coef(size(POD_state,1)*u%dim+1:size(POD_state,1)*(u%dim+1))=  &
                           pod_coef(size(POD_state,1)*u%dim+1:size(POD_state,1)*(u%dim+1))+ &
                           pod_coef_dt(size(POD_state,1)*u%dim+1:size(POD_state,1)*(u%dim+1))

                      ! for pressure, here we solve p^n+1    
                     ! pod_coef(size(POD_state,1)*u%dim+1:size(POD_state,1)*(u%dim+1))=  &
                      !     - pod_coef_dt(size(POD_state,1)*u%dim+1:size(POD_state,1)*(u%dim+1))


                      !pod_coef(:)=pod_coef_dt(:)
                      !!print*,pod_coef(2)
                      rmsestep=timestep
                      !save pod_coef for timestep 2
                      write(40,*)(pod_coef(i),i=1,(u%dim+1)*size(POD_state,1))
                      ! save pod_coef for all the time levels
                     
                      ! save the pod_coef at ttimestep =1
                      open(101,file='coef_pod_all', position='append',ACTION='WRITE')
                      write(101,*)(pod_coef(i),i=1,(u%dim+1)*size(POD_state,1))
                      close(101)

                      call project_full(delta_u, delta_p, pod_sol_velocity, pod_sol_pressure,POD_state(:,:,istate), pod_coef)
                      !save u1 
                      
                      ! print*,'pod_coef at timestep1'
                      !  print*,pod_coef
                   endif
                   
                else !timestep.gt.1
                   ewrite(1,*)'timestep=',timestep
                    
                   POD_velocity=>extract_vector_field(POD_state(1,1,istate), "Velocity")
                   POD_pressure=>extract_scalar_field(POD_state(1,2,istate), "Pressure")
                   
                   !pod_coef is from the previous timestep, argument?
                   if(timestep==1) then
                      open(40,file='coef_pod_all0')
                   else
                      open(40,file='pod_coef')
                   endif
                   !open(40,file='pod_coef')
                   read(40,*)(pod_coef(i),i=1,(u%dim+1)*size(POD_state,1))
                   close(40)

                   if(timestep==1) then
                      open(40,file='coef_pod_all')
                      write(40,*)(pod_coef(i),i=1,(u%dim+1)*size(POD_state,1))
                      close(40)
                   endif
                   !open(40,file='pod_coef')


                   open(unit=20,file='pod_matrix_snapmean')
                   read(20,*)((pod_matrix_snapmean(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
                   read(20,*)(pod_rhs_snapmean(i),i=1,(u%dim+1)*size(POD_state,1))
                   close(20)
                   
                   pod_matrix%val=pod_matrix_snapmean(:,:)
                   pod_rhs%val=pod_rhs_snapmean(:)
                   pod_matrix_adv_perturbed =0.0
                   pod_rhs_adv_perturbed =0.0
                   
                   do k=1,size(POD_state,1)
                      do d=1, u%dim
                         
                         read(30,*)((pod_matrix_perturbed(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
                         read(50,*)(pod_rhs_perturbed(i),i=1,(u%dim+1)*size(POD_state,1))
                         
                         pod_matrix%val=pod_matrix%val+pod_coef(k+(d-1)*size(POD_state,1))*pod_matrix_perturbed(:,:)
                         pod_rhs%val=pod_rhs%val+pod_coef(k+(d-1)*size(POD_state,1))*pod_rhs_perturbed(:)
                         
                         read(60,*)((pod_matrix_adv_perturbed(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
                         !                        pod_matrix_adv_perturbed= pod_coef(k+(d-1)*size(POD_state,1))*pod_matrix_adv_perturbed(:,:)
                         pod_matrix_perturbed= pod_coef(k+(d-1)*size(POD_state,1))*pod_matrix_perturbed(:,:)
                         !                        call Matrix_vector_multiplication(size(pod_coef),pod_rhs_adv_perturbed, &
                         !                             pod_matrix_adv_perturbed,pod_coef)
                         call Matrix_vector_multiplication(size(pod_coef),pod_rhs_adv_perturbed, &
                              pod_matrix_perturbed,pod_coef)
                         pod_rhs%val = pod_rhs%val-pod_rhs_adv_perturbed/(theta*dt) ! add advection term in rhs
                         
                      enddo
                   
                      read(30,*)((pod_matrix_perturbed(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
                      read(50,*)(pod_rhs_perturbed(i),i=1,(u%dim+1)*size(POD_state,1))
                      
                      pod_matrix%val=pod_matrix%val+pod_coef(k+u%dim*size(POD_state,1))*pod_matrix_perturbed(:,:)
                      pod_rhs%val=pod_rhs%val+pod_coef(k+u%dim*size(POD_state,1))*pod_rhs_perturbed(:)
                   enddo
                   
                      
                   call solve(pod_matrix%val,pod_rhs%val)  
                   pod_coef_dt=0.0
                   pod_coef_dt=pod_rhs%val
                   pod_rhs%val=0.0
                   ! for  velocity, here we solved du/dt
                   pod_coef(1:size(POD_state,1)*u%dim)=pod_coef(1:size(POD_state,1)*u%dim)+dt*pod_coef_dt(1:size(POD_state,1)*u%dim)
                   ! for pressure, here we solve dp  
                   pod_coef(size(POD_state,1)*u%dim+1:size(POD_state,1)*(u%dim+1))=  &
                        pod_coef(size(POD_state,1)*u%dim+1:size(POD_state,1)*(u%dim+1))- &
                        pod_coef_dt(size(POD_state,1)*u%dim+1:size(POD_state,1)*(u%dim+1))
                     ! for pressure, here we solve p^n+1    
                    ! pod_coef(size(POD_state,1)*u%dim+1:size(POD_state,1)*(u%dim+1))=  &
                      !   -pod_coef_dt(size(POD_state,1)*u%dim+1:size(POD_state,1)*(u%dim+1))
                   
                   !pod_coef(:)=pod_coef_dt(:)
                   
                   if(.true.) then
                   if(timestep.eq.total_timestep) then
                      
                      open(101,file='coef_pod_all0')
                      read(101,*)(pod_coef0(i),i=1,(u%dim+1)*size(POD_state,1))
                      close(101)
                      call project_full(delta_u, delta_p, pod_sol_velocity, pod_sol_pressure, POD_state(:,:,istate), pod_coef0)
                      snapmean_velocity=>extract_vector_field(POD_state(1,1,istate),"SnapmeanVelocity")
                      snapmean_pressure=>extract_Scalar_field(POD_state(1,2,istate),"SnapmeanPressure")
                      !call addto(u, delta_u, dt)
                      u%val=snapmean_velocity%val
                      call addto(u, delta_u)
                      p%val=snapmean_pressure%val
                      call addto(p, delta_p)
                      dump_no_tmp=total_timestep
                      call write_state(dump_no_tmp, state)                      
                   endif
                   endif

                  if (have_option("/reduced_model/Non_intrusive")) then
                   if (have_option("/reduced_model/training")) then
                          ctrl_timestep=1499
                   else 
                          ctrl_timestep=19
                   endif
                    
                 
                   if(timestep>ctrl_timestep) then !
                      
                        ! i= (u%dim)*size(POD_state,1) 
                        ! call  ann_bp_main(total_timestep,timestep,pod_coef(i),i) 
                        do i=1,  (u%dim)*size(POD_state,1) 
                          output=0.0  
                          pod_coef(i)=0           
                          call  ann_bp_main(total_timestep,timestep,pod_coef(i),i) !recalculate the pod_coef(i) using BP_ANN one by one
                        enddo 
                    endif 
                  !  if (timestep.eq.260) stop 300
                   !save pod_coef, rewrite every timestep 
                  endif ! non-intrusive
                   open(40,file='pod_coef')
                   write(40,*)(pod_coef(i),i=1,(u%dim+1)*size(POD_state,1)) ! ignore for ann
                   close(40)
                   

                   ! save pod_coef for all the time levels
                   
                   open(101,file='coef_pod_all', position='append',ACTION='WRITE')
                   write(101,*)(pod_coef(i),i=1,(u%dim+1)*size(POD_state,1))
                   close(101)
                   print * ,'pod_coef(i)', pod_coef((u%dim)*size(POD_state,1))

                   allocate(pod_coef_all_temp(total_timestep,((u%dim+1)*size(POD_state,1))))
                   open(11,file='coef_pod_all_obv')
                  !open(11,file='coef_pod_all_obv',position='append',ACTION='READ')
                  read(11,*)((pod_coef_all_temp(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,total_timestep)                  
                  close(11)
                    !pod_coef(:)=pod_coef_all_temp(timestep,:) ! test the rom using full model's results

                   rmsestep=timestep
                   call project_full(delta_u, delta_p, pod_sol_velocity, pod_sol_pressure, POD_state(:,:,istate), pod_coef)
                   
                   deallocate(pod_coef_all_temp)
                endif ! timestep          
                
                !1000 continue
                
                snapmean_velocity=>extract_vector_field(POD_state(1,1,istate),"SnapmeanVelocity")
                snapmean_pressure=>extract_Scalar_field(POD_state(1,2,istate),"SnapmeanPressure")
                !call addto(u, delta_u, dt)
                u%val=snapmean_velocity%val
                call addto(u, delta_u)
                p%val=snapmean_pressure%val
                call addto(p, delta_p)

             endif ! adjoint             
             call deallocate(delta_p)
             call deallocate(delta_u)
             deallocate(pod_matrix_snapmean)
             deallocate(pod_rhs_snapmean)
             deallocate(pod_matrix_perturbed)
             deallocate(pod_matrix_adv_perturbed)
             deallocate(pod_rhs_perturbed)
             deallocate(pod_rhs_adv_perturbed)
             deallocate(pod_ct_m)
             call deallocate(pod_matrix)
             call deallocate(pod_matrix_mass)
             call deallocate(pod_rhs)
             deallocate(pod_coef)

             if(timestep.eq.total_timestep) then
                deallocate(pod_coef0)
             endif
           
          endif ! prognostic velocity
       enddo reduced_model_loop
       
        end subroutine solve_momentum_reduced

        SUBROUTINE Matrix_vector_multiplication(size_vector,A_phi,Mat_A,phi)
          IMPLICIT NONE
          integer, intent(in):: size_vector
          real, dimension(:,:), intent(in):: Mat_A
          real, dimension(:), intent(in):: phi
          real, intent(out):: A_phi(size_vector)
          
          ! local
          integer ii,jj,kk,i,j,k,nod,row,column
          
          A_phi=0.0
          
          do row= 1,size(Mat_A(:,1))
             do column = 1,size_vector
                A_phi(row)=A_phi(row)+Mat_A(row,column)*phi(column)
             enddo
          enddo
          
        END SUBROUTINE Matrix_vector_multiplication


        function inv(A) result(Ainv)
          ! Returns the inverse of a matrix calculated by finding the LU
          ! decomposition.  Depends on LAPACK.
          real, dimension(:,:), intent(in) :: A
          real, dimension(size(A,1),size(A,2)) :: Ainv
          
          !real, dimension(size(A,1)) :: work  ! work array for LAPACK
          real, dimension(:),allocatable :: work
          !  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
          integer, dimension(:),allocatable :: ipiv      
          integer :: n, info
          
          ! External procedures defined in LAPACK
          external DGETRF
          external DGETRI
          allocate(work(size(A,1)))
          allocate(ipiv(size(A,1)))
          ! Store A in Ainv to prevent it from being overwritten by LAPACK
          Ainv = A
          n = size(A,1)
          
          ! DGETRF computes an LU factorization of a general M-by-N matrix A
          ! using partial pivoting with row interchanges.
          call DGETRF(n, n, Ainv, n, ipiv, info)
          
          if (info /= 0) then
             stop 'Matrix is numerically singular!'
          end if
          
          ! DGETRI computes the inverse of a matrix using the LU factorization
          ! computed by DGETRF.
          call DGETRI(n, Ainv, n, ipiv, work, n, info)
          
          if (info /= 0) then
             stop 'Matrix inversion failed!'
          end if
          deallocate(work)
          deallocate(ipiv)
        end function inv

        subroutine solve_reduced_adjoint(state, timestep, POD_state,total_timestep,pod_coef_adjoint)
        type(state_type), dimension(:), intent(inout) :: state
        type(state_type), dimension(:,:,:), intent(inout) :: POD_state
        integer,intent(in) :: total_timestep
        real, dimension(:,:),intent(inout) ::pod_coef_adjoint
        integer :: timestep
        real :: dt
        integer :: d,i,k,j
        ! Adjoint model
      	logical :: reduced_adjoint
      	type(pod_rhs_type) :: pod_rhs_adjoint
      	type(pod_matrix_type) :: adjoint_pod_A, adjoint_pod_B,adjoint_A_extra,adjoint_A ,pod_matrix,pod_matrix_B
      	type(pod_matrix_type) :: adjoint_pod_A0, adjoint_pod_B0, pod_matrix_B0
      	real, dimension(:,:), allocatable :: pod_coef_all,pod_coef_all_obv,pod_matrix_perturbed,&
                     pod_matrix_snapmean ! solution of adjoint
        real, dimension(:,:), allocatable :: pod_coef_adjoint_tmp
        type(vector_field), pointer :: u
        type(scalar_field), pointer :: p
        type(vector_field), pointer :: POD_velocity,snapmean_velocity
        type(scalar_field), pointer :: POD_pressure,snapmean_pressure
        type(pod_matrix_type) :: pod_matrix_mass
        ! Change in pressure
        type(scalar_field) :: delta_p
        ! Change in velocity
        type(vector_field) :: delta_u
        ! optimal sensor location
        integer, dimension(:), allocatable :: optimal_location
        real, dimension(:), allocatable :: u_total
        
      	! gradient:
      	real, dimension(:), allocatable :: g,ds,ds_tmp

       logical :: adaptive_observation

       ewrite(1,*) 'Adjoint model, timestep',timestep,total_timestep
       ewrite(1,*) '================================================'
!      	call get_option("/timestepping/current_time", current_time)
!      	call get_option("/timestepping/finish_time", finish_time)       
      	call get_option("/timestepping/timestep", dt)
!        total_timestep=int((finish_time-current_time)/dt)+1
        u=>extract_vector_field(state(1), "Velocity")
        p=>extract_scalar_field(state(1), "Pressure")
        POD_velocity=>extract_vector_field(POD_state(1,1,1), "Velocity")
        POD_pressure=>extract_scalar_field(POD_state(1,2,1), "Pressure")
      	call get_option(trim(u%option_path)//"/prognostic/temporal_discretisation/theta",theta)
      	call allocate(adjoint_pod_A, POD_velocity, POD_pressure)
      	call allocate(adjoint_A_extra, POD_velocity, POD_pressure) 
      	call allocate(adjoint_A, POD_velocity, POD_pressure)
      	call allocate(adjoint_pod_B, POD_velocity, POD_pressure)
        call allocate(pod_matrix, POD_velocity, POD_pressure)
        call allocate(pod_matrix_B, POD_velocity, POD_pressure)
        call allocate(pod_rhs_adjoint, POD_velocity, POD_pressure)
        call allocate(pod_matrix_mass, POD_velocity, POD_pressure)
      	call allocate(adjoint_pod_B0, POD_velocity, POD_pressure)
      	call allocate(adjoint_pod_A0, POD_velocity, POD_pressure)
        call allocate(pod_matrix_B0, POD_velocity, POD_pressure)

        allocate(pod_coef_adjoint_tmp(0:total_timestep,(u%dim+1)*size(POD_state,1)))
        pod_coef_adjoint_tmp = 0.0
      	allocate(pod_coef_all(0:total_timestep,(u%dim+1)*size(POD_state,1)))
      	allocate(pod_coef_all_obv(0:total_timestep,(u%dim+1)*size(POD_state,1)))
        allocate(g ((u%dim+1)*size(POD_state,1))) 
        allocate(ds((u%dim+1)*size(POD_state,1)))
      	allocate(ds_tmp((u%dim+1)*size(POD_state,1)))
        allocate(pod_matrix_perturbed((u%dim+1)*size(POD_state,1),(u%dim+1)*size(POD_state,1)))
        !allocate(pod_matrix((u%dim+1)*size(POD_state,1),(u%dim+1)*size(POD_state,1)))
        !allocate(pod_matrix_B((u%dim+1)*size(POD_state,1),(u%dim+1)*size(POD_state,1)))
        allocate(pod_matrix_snapmean((u%dim+1)*size(POD_state,1),(u%dim+1)*size(POD_state,1)))
!        allocate(pod_rhs_snapmean((u%dim+1)*size(POD_state,1)))
        open(unit=61,file='coef_pod_all')
        read(61,*)((pod_coef_all(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=0,total_timestep)
        close(61)
!        open(unit=20,file='pod_matrix_snapmean')
!        read(20,*)((pod_matrix_snapmean(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
!        close(20)
        open(unit=20,file='mass_matrix')
        read(20,*)((pod_matrix_mass%val(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
        close(20)

        open(unit=21,file='advection_matrix_snapmean')
        read(21,*)((pod_matrix_snapmean(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
        close(21)

        ! Allocate the change in pressure field
        call allocate(delta_p, p%mesh, "DeltaP")
        delta_p%option_path = trim(p%option_path)
        call zero(delta_p)
        
        call allocate(delta_u, u%dim, u%mesh, "DeltaU")
        delta_u%option_path = trim(u%option_path)
        call zero(delta_u)

      	g = 0.0
      	ds = 0.0
      	ds_tmp = 0.0

      	pod_matrix%val=pod_matrix_snapmean(:,:)/(theta*dt)
      	pod_matrix_B%val=pod_matrix_snapmean(:,:)/(theta*dt)
      	pod_matrix_B0%val=pod_matrix_snapmean(:,:)/(theta*dt)
      
      	open(1,file='coef_pod_all_obv')
      	read(1,*)((pod_coef_all_obv(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,total_timestep)
      	! pod_coef_adjoint() save the coef_pod_all_obv temporarily
      	close(1)
      	open(30,file='advection_matrix_perturbed')
      	adjoint_pod_A%val =0.0
      	adjoint_pod_B%val =0.0
      	do k=1,size(POD_state,1)
           do d=1, u%dim
            !adjoint_pod_A%val size same as pod_matrix. read pod_matrix_perturbed here including theta*dt
            read(30,*)((pod_matrix_perturbed(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
            !read(30,*)((pod_matrix_perturbed(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))
            ! delete theta*dt from pod_matrix_perturbed
            pod_matrix_perturbed = pod_matrix_perturbed/(theta*dt)
            pod_matrix%val=pod_matrix%val+  &
                 pod_coef_all(timestep-1,k+(d-1)*size(POD_state,1))*pod_matrix_perturbed(:,:)
            pod_matrix_B%val=pod_matrix_B%val+  &
                 pod_coef_all(timestep,k+(d-1)*size(POD_state,1))*pod_matrix_perturbed(:,:)
            adjoint_pod_A%val((d-1)*size(POD_state,1)+k,:)= &
                 theta*matmul(pod_matrix_perturbed(:,:),pod_coef_all(timestep+1, :))
            adjoint_pod_B%val((d-1)*size(POD_state,1)+k,:)= &
                 (1.-theta)*matmul(pod_matrix_perturbed(:,:),pod_coef_all(timestep, :))

            !if(timestep.eq.1) then
               adjoint_pod_A0%val((d-1)*size(POD_state,1)+k,:)= &
                    theta*matmul(pod_matrix_perturbed(:,:),pod_coef_all(timestep, :))
               adjoint_pod_B0%val((d-1)*size(POD_state,1)+k,:)= &
                    (1.-theta)*matmul(pod_matrix_perturbed(:,:),pod_coef_all(timestep-1, :))
               pod_matrix_B0%val=pod_matrix_B0%val+  &
                    pod_coef_all(timestep-1,k+(d-1)*size(POD_state,1))*pod_matrix_perturbed(:,:)
            !endif
            
         enddo
         
         read(30,*)((pod_matrix_perturbed(i,j),j=1,(u%dim+1)*size(POD_state,1)),i=1,(u%dim+1)*size(POD_state,1))

         pod_matrix%val=pod_matrix%val+pod_coef_all(timestep-1,k+u%dim*size(POD_state,1))*pod_matrix_perturbed(:,:)
         pod_matrix_B%val=pod_matrix_B%val+pod_coef_all(timestep,k+u%dim*size(POD_state,1))*pod_matrix_perturbed(:,:)

         adjoint_pod_A%val((d-1)*size(POD_state,1)+k,:)=&
              theta*matmul(pod_matrix_perturbed(:,:),pod_coef_all(timestep+1, :))
         adjoint_pod_B%val((d-1)*size(POD_state,1)+k,:)=&
              (1.-theta)*matmul(pod_matrix_perturbed(:,:),pod_coef_all(timestep, :))

         !if(timestep.eq.1) then
            adjoint_pod_A0%val((d-1)*size(POD_state,1)+k,:)=&
                 theta*matmul(pod_matrix_perturbed(:,:),pod_coef_all(timestep, :))
            adjoint_pod_B0%val((d-1)*size(POD_state,1)+k,:)=&
                 (1.-theta)*matmul(pod_matrix_perturbed(:,:),pod_coef_all(timestep-1, :))
            pod_matrix_B0%val=pod_matrix_B0%val+pod_coef_all(timestep-1,k+u%dim*size(POD_state,1))*pod_matrix_perturbed(:,:)
        !endif
      enddo
      
      	close(30)

       pod_matrix_B%val = (1.-theta)*pod_matrix_B%val - (pod_matrix_mass%val)/dt
       pod_matrix_B0%val = (1.-theta)*pod_matrix_B0%val - (pod_matrix_mass%val)/dt
       pod_matrix%val = theta*pod_matrix%val + (pod_matrix_mass%val)/dt


       if(timestep.eq.total_timestep-1) then
          pod_coef_adjoint_tmp(timestep+1,:) = 0.0
       else
          open(10,file='pod_coef_adj.dat')
          read(10,*) (pod_coef_adjoint_tmp(timestep+1,i),i=1,(u%dim+1)*size(POD_state,1))
          close(10)
       endif

       
       !    adjoint_A_extra%val=adjoint_pod_A%val+ adjoint_pod_B%val        
       adjoint_A_extra%val=adjoint_pod_A%val + adjoint_pod_B%val + pod_matrix_B%val
       adjoint_A_extra%val=transpose(adjoint_A_extra%val)
       adjoint_A%val=transpose(pod_matrix%val)     
       pod_rhs_adjoint%val=-matmul(adjoint_A_extra%val,pod_coef_adjoint_tmp(timestep+1,:))    

!       if(2*int(timestep/2).eq.timestep) then
       if(timestep.eq.total_timestep-1) then
 !      if(timestep.eq.total_timestep-1) then
          pod_rhs_adjoint%val(1:u%dim*size(POD_state,1))= pod_rhs_adjoint%val(1:u%dim*size(POD_state,1))+ &
               pod_coef_all_obv(timestep,1:u%dim*size(POD_state,1))-pod_coef_all(timestep,1:u%dim*size(POD_state,1))
       endif

       call solve(adjoint_A%val,pod_rhs_adjoint%val)
       pod_coef_adjoint_tmp(timestep,:)=pod_rhs_adjoint%val(:)
       print*,'222',pod_rhs_adjoint%val(:)
       open(10,file='pod_coef_adj.dat')
       write(10,*) (pod_coef_adjoint_tmp(timestep,i),i=1,(u%dim+1)*size(POD_state,1))
       close(10)

       !if(timestep.eq.1) then !! this is only for initial conditions.... we are now looking for the sensitivity w.r.t. results at each time level.
         !ds(:) = matmul(-adjoint_pod_A0%val-adjoint_pod_B0%val-pod_matrix_B0%val, pod_coef_all(0, :))  !!!Is pod_coef_all(1, :) the initial one??
          g = 0.0
          do i = 1,size(g)
             ds = 0.0
             ds(i) = 1.0
             ds_tmp = matmul(-adjoint_pod_A0%val-adjoint_pod_B0%val-pod_matrix_B0%val,ds)
           ! print*,'ds_tmp', ds_tmp
             g(i)= g(i) + dot_product(pod_rhs_adjoint%val,ds_tmp)
          enddo
          
          call project_adjoint_gradient(delta_u, delta_p, POD_state(:,:,1), g)
          
                snapmean_velocity=>extract_vector_field(POD_state(1,1,1),"SnapmeanVelocity")
                snapmean_pressure=>extract_Scalar_field(POD_state(1,2,1),"SnapmeanPressure")
                !call addto(u, delta_u, dt)
                !if(timestep==total_timestep) then
                   u%val=0.0
                   call addto(u, delta_u)
                   p%val=0.0
                   call addto(p, delta_p)
                !endif
          call addto(u, delta_u)
          call addto(p, delta_p)

          ! sort out the order of magniture of delta_u -- important map of sensor location
          adaptive_observation = .true.
          if(adaptive_observation) then
             if(timestep.eq.1) then
                allocate(optimal_location(size(delta_u%val,2)))
                allocate(u_total(size(delta_u%val,2)))
                u_total=0.0
                do i =1,size(delta_u%val,2)
                   do d =1,u%dim
                      u_total(i) =  u_total(i)+delta_u%val(d,i)**2
                   enddo
                   u_total(i) = sqrt(u_total(i))
                enddo
             endif
          endif
          

          print*,'iiiiiii',g
          if(timestep.eq.1) then
             open(unit=2,file='adjoint_g')
             write(2,*)(g(i),i=1, (u%dim+1)*size(POD_state,1))
             close(2)
             
          endif
     
          if(timestep.eq.total_timestep-1) then
              open(unit=2,file='adjoint_g_alltime')
              write(2,*)(g(i),i=1, (u%dim+1)*size(POD_state,1))
              close(2)
          else
              open(unit=2,file='adjoint_g_alltime',ACTION='WRITE',position='append')
              write(2,*)(g(i),i=1, (u%dim+1)*size(POD_state,1))
              close(2)
          endif

       !endif

       pod_coef_adjoint = pod_coef_adjoint_tmp

       call deallocate(adjoint_pod_A)
       call deallocate(adjoint_A_extra)
       call deallocate(adjoint_A)
       call deallocate(adjoint_pod_B)
       deallocate(pod_coef_adjoint_tmp)
       deallocate(pod_coef_all)
       deallocate(ds)
       deallocate(ds_tmp)
       deallocate(pod_matrix_perturbed)
       call deallocate(pod_matrix)
       call deallocate(pod_matrix_B)
       deallocate(pod_matrix_snapmean)
       call deallocate(pod_matrix_mass)

       deallocate(optimal_location)
       deallocate(u_total)
     end subroutine solve_reduced_adjoint

     subroutine max_sort(n,a,counter)
       implicit none
       integer :: n
       real, intent(in) :: a(n)
       real :: b(n)
       integer, intent(out) :: counter(n)
       integer i,j  ! loop counter
       real max  ! find a maximum value of one loop
       real temp ! 
       integer temp2 ! 
       !sort

       b=a
       do i=1,n
          counter(i) = i
       enddo
       do i=1,n-1
          max=b(i)     ! 
          do j=i+1,n
!             print*,max,b(j)
             if ( max < b(j) ) then      
!             print*,'2222',max,b(j)
                 temp=b(j)          
                 temp2=counter(j)          
                b(j)=b(i)
                b(i)=temp
                counter(j)= counter(i)
                counter(i)= temp2
                max=b(i)
             end if
          end do
       end do
     end subroutine max_sort


     
          
  SUBROUTINE NEUEQNSTA(TEMP,PRES,DENSTY,EQNSTA)

       IMPLICIT NONE
       REAL TEMP,PRES!,DENSTY
       REAL, intent(inout) :: DENSTY
       INTEGER EQNSTA
! This sub finds the density DENSTY from TEMP and PRES
!       INTEGER :: NLAYER1,NLAFIN,NTRAIN,NONODS,MXNCOLM,NLAYERS,NITS
! NLAYER1,NLAFIN no of nodes in the 1st (input layer) and final (output layer)
! EQNSTA=7 lime EoS
! EQNSTA=8 Silica EoS
! EQNSTA=9 granite EoS
       INTEGER, PARAMETER :: NLAYER1=2,NLAFIN=1,NONODS=15
!       INTEGER, PARAMETER :: NLAYER1=2,NLAFIN=1,NONODS=17,NTRAIN=10000
       INTEGER, PARAMETER :: MXNCOLM=NONODS*NONODS,NLAYERS=4
       REAL :: WEIGHT(MXNCOLM),WEINEW(MXNCOLM)
       REAL :: MINT,MAXT,MINP,MAXP,MIND,MAXD
       REAL :: FUNC,FUNNEW,PERTERB
       INTEGER :: COLM(MXNCOLM),FINDRM(NONODS+1)
       INTEGER :: NLAY(NLAYERS)
       INTEGER :: ITNGOOD,ITS,COUNT,ILAYER,I,NCOLM, IRED,ITRAIN,I1
       INTEGER :: LOPT,NOREXP
!
       NLAY(1)=2
       NLAY(2)=6
       NLAY(3)=6
       NLAY(4)=1
!
       LOPT=EQNSTA-6
!
! Define the sparcity of the matrix and set WEIGHT
       CALL SPAWEI(WEIGHT,MXNCOLM,NONODS,  &
                   NCOLM,FINDRM,COLM,NLAY,NLAYERS)
!
       CALL DEFNWEI(WEIGHT,NCOLM,LOPT,NOREXP,  & 
         MINT,MAXT,MINP,MAXP,MIND,MAXD)
!
       CALL EOS_NEU(TEMP,PRES,DENSTY, &
           MINT,MAXT,MINP,MAXP,MIND,MAXD, &
           WEIGHT,NONODS, MXNCOLM, &
           NCOLM,FINDRM,COLM,NLAYER1,NLAFIN,NOREXP)
!
       RETURN
       END SUBROUTINE NEUEQNSTA
	
   SUBROUTINE SPAWEI(WEIGHT,MXNCOLM,NONODS, &
                   NCOLM,FINDRM,COLM,NLAY,NLAYERS)
! Define the sparcity of the matrix and set WEIGHT

       IMPLICIT NONE
       INTEGER :: MXNCOLM,NONODS,NCOLM,NLAYERS
       REAL :: WEIGHT(MXNCOLM)
       !!INTEGER :: NLAY(NLAYERS),FINDRM(NONODS+1),COLM(NCOLM)
       INTEGER :: NLAY(NLAYERS),FINDRM(NONODS+1),COLM(MXNCOLM)
! Local variables
       INTEGER :: NOD,COUNT,ILAYER,II,JJ,INLAYS
!
       NOD=0
       COUNT=0
       INLAYS=0
       DO ILAYER=1,NLAYERS
         DO II=1,NLAY(ILAYER)
           NOD=NOD+1
           FINDRM(NOD)=COUNT+1
           IF(ILAYER.EQ.1) THEN
           ELSE
             DO JJ=1,NLAY(ILAYER-1)
               COUNT=COUNT+1
               COLM(COUNT)=INLAYS-NLAY(ILAYER-1)+JJ
!               ewrite(3,*) 'COUNT,COLM(COUNT):',COUNT,COLM(COUNT)
!               ewrite(3,*) 'INLAYS-NLAY(ILAYER-1)+JJ:', &
!                        INLAYS,NLAY(ILAYER-1),JJ
!               STOP 7
             END DO
           ENDIF
         END DO
         INLAYS=INLAYS+NLAY(ILAYER) 
       END DO
       IF(NONODS.NE.NOD) THEN
          ewrite(3,*) 'NONODS WAS NOT SET TO CORRECT VALUE'
          ewrite(3,*) 'NOD,NONODS:',NOD,NONODS
          STOP 77
       ENDIF
       NCOLM=COUNT
       FINDRM(NONODS+1)=NCOLM+1
!       ewrite(3,*) 'NONODS,NCOLM=',NONODS,NCOLM
!        stop 4
!
!       DO NOD=1,NONODS
!         ewrite(3,*) 'NOD,FINDRM(NOD),FINDRM(NOD+1)-1:', &
!                  NOD,FINDRM(NOD),FINDRM(NOD+1)-1
!         ewrite(3,*) 'COLM:', &
!                (COLM(COUNT),COUNT=FINDRM(NOD),FINDRM(NOD+1)-1)
!       END DO
!       STOP 33
!
       DO COUNT=1,NCOLM
         WEIGHT(COUNT)=0.
       END DO
       RETURN
       END SUBROUTINE SPAWEI

       SUBROUTINE GETNEUVALS(NEUVAL,WEIGHT,NONODS, MXNCOLM, &
               NCOLM,FINDRM,COLM,NLAYER1,NLAFIN,NOREXP)
! NLAYER1,NLAFIN no of nodes in the 1st (input layer) and final (output layer)
! This sub calculates the neuron values. 
! If NOREXP=1 the dont find exponent of output. 
       IMPLICIT NONE
       LOGICAL OUTEXP
!       PARAMETER(OUTEXP=.FALSE.)
!       PARAMETER(OUTEXP=.true.)
       INTEGER :: NLAYER1,NLAFIN,NONODS,NCOLM, MXNCOLM
       REAL :: NEUVAL(NONODS),WEIGHT(NCOLM)
       INTEGER :: FINDRM(NONODS+1),COLM(MXNCOLM)
       INTEGER :: NOREXP
! LOCAL VARIABLES...
       REAL :: SUM
       INTEGER :: NOD,COUNT
       OUTEXP=(NOREXP.EQ.0)
       DO NOD=NLAYER1+1,NONODS
         SUM=0.
         DO COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
           SUM=SUM+WEIGHT(COUNT)*NEUVAL(COLM(COUNT))
         END DO
         IF(NOD.EQ.NONODS) THEN
           IF(OUTEXP) THEN
              NEUVAL(NOD)=1./(1.+EXP(-SUM))
           ELSE
              NEUVAL(NOD)=SUM
           ENDIF
         ELSE
           NEUVAL(NOD)=1./(1.+EXP(-SUM))
         ENDIF
       END DO
       RETURN
       END SUBROUTINE GETNEUVALS

      
       SUBROUTINE DEFNWEI(WEIGHT,NCOLM,LOPT,NOREXP,  & 
         MINT,MAXT,MINP,MAXP,MIND,MAXD)

       IMPLICIT NONE
       INTEGER NCOLM,LOPT,NOREXP
       REAL WEIGHT(NCOLM)
       REAL MINT,MAXT,MINP,MAXP,MIND,MAXD
!
       IF(LOPT.EQ.1) THEN
! Eqn of state for lime...
       CALL DEFNWEILIM(WEIGHT,NCOLM,NOREXP,  & 
         MINT,MAXT,MINP,MAXP,MIND,MAXD)
       ELSE IF(LOPT.EQ.2) THEN
! Eqn of state for Silica...
       CALL DEFNWEISIL(WEIGHT,NCOLM,NOREXP,  & 
         MINT,MAXT,MINP,MAXP,MIND,MAXD)
       ELSE IF(LOPT.EQ.3) THEN
! Eqn of state for Granite...
       CALL DEFNWEIGRA(WEIGHT,NCOLM,NOREXP,  & 
         MINT,MAXT,MINP,MAXP,MIND,MAXD)
       ENDIF
!
        RETURN
        END SUBROUTINE DEFNWEI
!
!
!
!
       SUBROUTINE DEFNWEILIM(WEIGHT,NCOLM,NOREXP,  & 
         MINT,MAXT,MINP,MAXP,MIND,MAXD)

       IMPLICIT NONE
       INTEGER NCOLM,NOREXP
       REAL WEIGHT(NCOLM)
       REAL MINT,MAXT,MINP,MAXP,MIND,MAXD
!
! Eqn of state for lime...
        NOREXP=0 
!
        WEIGHT(1)=      5.97685132589016     
        WEIGHT(2)=      1.10400778174330     
        WEIGHT(3)=      -21.2825376649738     
        WEIGHT(4)=      0.169090329919894     
        WEIGHT(5)=       20.1788017885642     
        WEIGHT(6)=      0.131574759151893     
        WEIGHT(7)=       12.3755801537776     
        WEIGHT(8)=      0.766587143843660     
        WEIGHT(9)=      -24.4043122338779     
        WEIGHT(10)=      0.232113849596108     
        WEIGHT(11)=       4.85706354480158     
        WEIGHT(12)=      -1.26789757374056     
        WEIGHT(13)=     -0.255388642114362     
        WEIGHT(14)=       12.0554633479762     
        WEIGHT(15)=      -6.94915549914577     
        WEIGHT(16)=      -2.09264002488505     
        WEIGHT(17)=       15.8594342730431     
        WEIGHT(18)=       1.72679479920685     
        WEIGHT(19)=      -2.06771382299692     
        WEIGHT(20)=      0.395823993317512     
        WEIGHT(21)=      -1.12455553498827     
        WEIGHT(22)=      0.750060594828952     
        WEIGHT(23)=     -0.473464061154712     
        WEIGHT(24)=      -1.09279477569545     
        WEIGHT(25)=       4.10419973811868     
        WEIGHT(26)=      -15.2397251897475     
        WEIGHT(27)=       14.8805451851581     
        WEIGHT(28)=       9.99889550046277     
        WEIGHT(29)=      -20.6890317411573     
        WEIGHT(30)=       5.36436208014150     
        WEIGHT(31)=      -1.24892426108564     
        WEIGHT(32)=       1.51633332142468     
        WEIGHT(33)=      -3.14208255247642     
        WEIGHT(34)=      -2.56368066528082     
        WEIGHT(35)=       9.65149011111481D-002
        WEIGHT(36)=     -0.507349715070754     
        WEIGHT(37)=      0.423474661102974     
        WEIGHT(38)=      -29.6165259232046     
        WEIGHT(39)=       12.3201603741371     
        WEIGHT(40)=       4.00375892664630     
        WEIGHT(41)=      -33.2267781535008     
        WEIGHT(42)=      -3.26392133049933     
        WEIGHT(43)=      -1.24302679976193     
        WEIGHT(44)=       4.59939718630899     
        WEIGHT(45)=      -5.25970543903552     
        WEIGHT(46)=      -3.00135185108571     
        WEIGHT(47)=       6.34139347514880     
        WEIGHT(48)=      -1.84842529815566     
        WEIGHT(49)=      -18.5019808067015     
        WEIGHT(50)=      1.88379314700385     
        WEIGHT(51)=      17.5824625711418     
        WEIGHT(52)=    -0.527333559763088     
        WEIGHT(53)=     -18.5389144969952     
        WEIGHT(54)=      6.24702737659637 
!    MINT,MAXT,MINP,MAXP,MIND,MAXD...
        MINT=298.150000000000        
        MAXT=37340.0952636910        
        MINP=1013250.00000000     
        MAXP=126898713.821683        
        MIND=1.31000000000000        
        MAXD=3.34412549863268 
        RETURN
        END SUBROUTINE DEFNWEILIM
!
!
!
!
       SUBROUTINE DEFNWEISIL(WEIGHT,NCOLM,NOREXP,  & 
         MINT,MAXT,MINP,MAXP,MIND,MAXD)

       IMPLICIT NONE
       INTEGER NCOLM,NOREXP
       REAL WEIGHT(NCOLM)
       REAL MINT,MAXT,MINP,MAXP,MIND,MAXD
!
! Eqn of state for Silica...
        NOREXP=0 
!
        WEIGHT(1)=   6.12626077572350     
        WEIGHT(2)=    0.272564908604740     
        WEIGHT(3)=    -9.75616050871453     
        WEIGHT(4)=    0.246912396109541     
        WEIGHT(5)=    -6.53574055759190     
        WEIGHT(6)=   -0.384054877796388     
        WEIGHT(7)=     -8.44375031552827     
        WEIGHT(8)=    -3.329123698799069E-002
        WEIGHT(9)=     -8.27742549631889     
        WEIGHT(10)=    -9.890057731467744E-002
        WEIGHT(11)=      9.73715337228838     
        WEIGHT(12)=    -0.240771751819873     
        WEIGHT(13)=     -4.49184124614962     
        WEIGHT(14)=      5.95516268582194     
        WEIGHT(15)=      2.82549172297567     
        WEIGHT(16)=    -1.894941028182234E-002
        WEIGHT(17)=     0.601560527201315     
        WEIGHT(18)=     -4.47260922589827     
        WEIGHT(19)=      2.09913243736491     
        WEIGHT(20)=     -20.5263926531534     
        WEIGHT(21)=     -12.2964085736977     
        WEIGHT(22)=     -15.5042442607703     
        WEIGHT(23)=     -17.9007234974105     
        WEIGHT(24)=      14.4442370983336     
        WEIGHT(25)=      1.31102712743980     
        WEIGHT(26)=      2.24711186705534     
        WEIGHT(27)=     0.862146059851389     
        WEIGHT(28)=      7.08865488987401     
        WEIGHT(29)=      4.18955246476026     
        WEIGHT(30)=    -0.989849174284737     
        WEIGHT(31)=     1.447487508202793E-002
        WEIGHT(32)=      3.48014767000874     
        WEIGHT(33)=      2.13069663254350     
        WEIGHT(34)=     -4.73426230781354     
        WEIGHT(35)=      4.49227690172869     
        WEIGHT(36)=      3.47395096346542     
        WEIGHT(37)=     -1.19391789986029     
        WEIGHT(38)=     -4.64835965448362     
        WEIGHT(39)=    -0.635582999915386     
        WEIGHT(40)=      2.15717035380372     
        WEIGHT(41)=     0.489452839441452     
        WEIGHT(42)=     -2.36516255020337     
        WEIGHT(43)=    -0.765529516551553     
        WEIGHT(44)=     -4.07263693916288     
        WEIGHT(45)=      2.99555734116876     
        WEIGHT(46)=     -3.86497423995451     
        WEIGHT(47)=     -2.35841013369974     
        WEIGHT(48)=     -1.43724880804879     
        WEIGHT(49)=      3.73463484818312     
        WEIGHT(50)=    -21.2831183705293     
        WEIGHT(51)=    0.384125180599244     
        WEIGHT(52)=     2.24339682643455     
        WEIGHT(53)=    -2.24290011447767     
        WEIGHT(54)=    0.146262916094321  
!    MINT,MAXT,MINP,MAXP,MIND,MAXD...   
   MINT=298.150000000000        
   MAXT=37340.0952636910        
   MINP=2026500.00000000     
   MAXP=127911963.821683       
   MIND=1.911301743049916E-002   
   MAXD=2.70021028300145 
   RETURN
   END SUBROUTINE DEFNWEISIL
!
!    
!
!
       SUBROUTINE DEFNWEIGRA(WEIGHT,NCOLM,NOREXP,  & 
         MINT,MAXT,MINP,MAXP,MIND,MAXD)

       IMPLICIT NONE
       INTEGER NCOLM,NOREXP
       REAL WEIGHT(NCOLM)
       REAL MINT,MAXT,MINP,MAXP,MIND,MAXD
!
! Eqn of state for Granite...
        NOREXP=1
!
        WEIGHT(1)=     -1.25210991901872     
        WEIGHT(2)=      1.84286851967128D-002
        WEIGHT(3)=     -1.14594954358655     
        WEIGHT(4)=      8.40714881584379D-002
        WEIGHT(5)=      2.65597303171804D-002
        WEIGHT(6)=      1.55508856931543D-002
        WEIGHT(7)=     0.836002693041384     
        WEIGHT(8)=     -9.05909371896358D-002
        WEIGHT(9)=    -0.268794055455616     
        WEIGHT(10)=     -7.17034776815134D-002
        WEIGHT(11)=      1.24561087811569     
        WEIGHT(12)=     0.115202509920321     
        WEIGHT(13)=     0.592933066624559     
        WEIGHT(14)=     0.295443685199725     
        WEIGHT(15)=     0.170527010853434     
        WEIGHT(16)=     -1.60166187085996D-002
        WEIGHT(17)=     0.116064480970310     
        WEIGHT(18)=    -0.256250630064082     
        WEIGHT(19)=      1.67727052032752     
        WEIGHT(20)=      1.19632898578538     
        WEIGHT(21)=      1.82714810873958D-002
        WEIGHT(22)=    -0.791636710562007     
        WEIGHT(23)=     0.656387670046782     
        WEIGHT(24)=    -0.925300133084808     
        WEIGHT(25)=    -0.356175762842861     
        WEIGHT(26)=    -0.550524110962181     
        WEIGHT(27)=    -0.298790079514946     
        WEIGHT(28)=     0.322634049945166     
        WEIGHT(29)=     -2.40837140811772D-002
        WEIGHT(30)=      1.03851390805661     
        WEIGHT(31)=     0.616274464012693     
        WEIGHT(32)=     0.443745520350746     
        WEIGHT(33)=    -0.190649055120316     
        WEIGHT(34)=    -0.388445592379031     
        WEIGHT(35)=      9.01767363969453D-003
        WEIGHT(36)=    -0.385043937917201     
        WEIGHT(37)=     0.545501292960052     
        WEIGHT(38)=     0.472350823673588     
        WEIGHT(39)=    -0.232523211547267     
        WEIGHT(40)=    -0.392064393082278     
        WEIGHT(41)=    -0.116452739687640     
        WEIGHT(42)=    -0.516573848801829     
        WEIGHT(43)=     0.225262483882030     
        WEIGHT(44)=    -0.101955071126241     
        WEIGHT(45)=      5.47774259524076D-002
        WEIGHT(46)=     -9.01064058715637D-002
        WEIGHT(47)=    -0.101044266061742     
        WEIGHT(48)=     -9.15587454508020D-002
        WEIGHT(49)=     -8.02903612304057D-002
        WEIGHT(50)=     2.16295456236704     
        WEIGHT(51)=    -1.68201258314670     
        WEIGHT(52)=    0.415871026552132     
        WEIGHT(53)=    0.591087476532526     
        WEIGHT(54)=   -0.257449667903680     
!    MINT,MAXT,MINP,MAXP,MIND,MAXD... 
   MINT=298.150000000000        
   MAXT=37340.0952636910        
   MINP=2026500.00000000     
   MAXP=127911963.821683        
!   MIND=2276.12204764802        
!   MAXD=2675.75066083717        
   MIND=2.27612204764802        
   MAXD=2.67575066083717     
   RETURN
   END SUBROUTINE DEFNWEIGRA

      SUBROUTINE EOS_NEU(TEMP,PRES,DENSTY, &
           MINT,MAXT,MINP,MAXP,MIND,MAXD, &
           WEIGHT,NONODS, MXNCOLM, &
           NCOLM,FINDRM,COLM,NLAYER1,NLAFIN,NOREXP)
! t=(max-min)*scaled+min
! This is the neural network EoS calculates density DENSTY 
! given the temp TEMP and pressure PRES. 
! MINT=min value of temp in data base used for training.
       IMPLICIT NONE
       INTEGER :: MXNODS
       PARAMETER(MXNODS=500)
       REAL :: NEUVAL(MXNODS)  
       INTEGER :: NLAYER1,NLAFIN,NONODS,NCOLM, MXNCOLM
       INTEGER :: FINDRM(NONODS+1),COLM(MXNCOLM)
       INTEGER :: NOREXP
       REAL :: TEMP,PRES!,DENSTY
       REAL, intent(inout) :: DENSTY
       REAL :: MINT,MAXT,MINP,MAXP,MIND,MAXD
       REAL :: WEIGHT(NCOLM)
! Local variables...
       REAL :: SCALT,SCALP,SCALD
         SCALT=(TEMP-MINT)/(MAXT-MINT)  
         SCALP=(PRES-MINP)/(MAXP-MINP) 
         NEUVAL(1)=SCALT
         NEUVAL(2)=SCALP
         CALL GETNEUVALS(NEUVAL,WEIGHT,NONODS, MXNCOLM, &
               NCOLM,FINDRM,COLM,NLAYER1,NLAFIN,NOREXP)
         SCALD=NEUVAL(NONODS)
! Unscale density SCALD to produce DENSTY
         DENSTY=(MAXD-MIND)*SCALD+MIND
         RETURN
         END SUBROUTINE EOS_NEU
!

   end module momentum_equation_reduced
