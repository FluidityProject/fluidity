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
      real, dimension(:,:), allocatable ::pod_coef_all_o,pod_coef_adjoint
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
                      x => extract_vector_field(state(istate), "Coordinate") 
                      ! get the coordinate values
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
                            stop 24
                         else
                            no_optimal_obv = 600
                            ! save the optimal coordinate values
                            allocate(u_tmp(size(u%val,2)))
                            allocate(optimal_location(size(u%val,2)))
                            open(10,file='optimal_sensor_location.dat')
                            read(10,*) (optimal_location(i),i=1,size(optimal_location))
                            close(10)  
if(.false.) then
                open(10,file='optimal_sensor_coordinate_50.dat')
                do i=1,50 ! only save 200 optimal location
                   j=optimal_location(i)
                   write(10,*) x%val(1,j),x%val(2,j)
                enddo
                close(10)
                open(10,file='optimal_sensor_coordinate_51_100.dat')
                do i=51,100 ! only save 200 optimal location
                   j=optimal_location(i)
                   write(10,*) x%val(1,j),x%val(2,j)
                enddo
                close(10)
                open(10,file='optimal_sensor_coordinate_100_200.dat')
                do i=101,200 ! only save 200 optimal location
                   j=optimal_location(i)
                   write(10,*) x%val(1,j),x%val(2,j)
                enddo
                close(10)
                open(10,file='optimal_sensor_coordinate_200_300.dat')
                do i=201,300 ! only save 200 optimal location
                   j=optimal_location(i)
                   write(10,*) x%val(1,j),x%val(2,j)
                enddo
                close(10)
    
stop 36
endif
                            ! u component
                            open(10,file = 'observation_u_ini.dat')
                            read(10,*) (u_tmp(i),i=1,size(u%val,2))
                            close(10)
                            do i=1,no_optimal_obv
                               j=optimal_location(i)
                               print*,u%val(1,j),u_tmp(j)
                  !             u%val(1,j)=u_tmp(j)
                               u%val(1,2000-i)=u_tmp(1000-i)
                  !             u%val(1,i)=u_tmp(i)
                            enddo
                            ! v component
                            open(10,file = 'observation_v_ini.dat')
                            read(10,*) (u_tmp(i),i=1,size(u%val,2))
                            close(10)
                            do i=1,no_optimal_obv
                               j=optimal_location(i)
                               print*,u%val(2,j),u_tmp(j)
                               ! use the optimal observational data
                  !             u%val(2,j)=u_tmp(j)
                               ! use the random observational data
                               u%val(2,2000-i)=u_tmp(1000-i)
                           !    u%val(2,i)=u_tmp(i)
                            enddo
                           ! p 
                            open(10,file = 'observation_p_ini.dat')
                            read(10,*) (u_tmp(i),i=1,size(p%val))
                            close(10)
                            do i=1,no_optimal_obv
                               j=optimal_location(i)
                               ! use the optimal observational data
                   !            p%val(j)=u_tmp(j)
                               ! use the random observational data
                               p%val(2000-i)=u_tmp(1000-i)
                    !           p%val(i)=u_tmp(i)
                            enddo
                            deallocate(u_tmp)
                            deallocate(optimal_location)
                         endif
                      endif

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
                      pod_coef(size(POD_state,1)*u%dim+1:size(POD_state,1)*(u%dim+1))=  &
                           - pod_coef_dt(size(POD_state,1)*u%dim+1:size(POD_state,1)*(u%dim+1))


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
                     pod_coef(size(POD_state,1)*u%dim+1:size(POD_state,1)*(u%dim+1))=  &
                         -pod_coef_dt(size(POD_state,1)*u%dim+1:size(POD_state,1)*(u%dim+1))
                   
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
                !      dump_no_tmp = total_timestep
                      call write_state(dump_no_tmp, state)
                   endif
                   endif
                   !save pod_coef, rewrite every timestep
                   open(40,file='pod_coef')
                   write(40,*)(pod_coef(i),i=1,(u%dim+1)*size(POD_state,1))
                   close(40)
                   ! save pod_coef for all the time levels

                   open(101,file='coef_pod_all', position='append',ACTION='WRITE')
                   write(101,*)(pod_coef(i),i=1,(u%dim+1)*size(POD_state,1))
                   close(101)

                   rmsestep=timestep
                   call project_full(delta_u, delta_p, pod_sol_velocity, pod_sol_pressure, POD_state(:,:,istate), pod_coef)
                   
                   
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

             if(timestep.eq.total_timestep-1) then
                
                open(10, file = 'num_last.dat')
                do i=1,u%dim
                   write(10,*) (u%val(i,j),j=1,node_count(u))
                enddo
                write(10,*) (p%val(j),j=1,node_count(p))
                close(10)
             endif

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
        type(vector_field), pointer :: u,x
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
       real,dimension(:), allocatable :: misfit_POD_rhs

       ewrite(1,*) 'Adjoint model, timestep',timestep,total_timestep
       ewrite(1,*) '================================================'
!      	call get_option("/timestepping/current_time", current_time)
!      	call get_option("/timestepping/finish_time", finish_time)       
      	call get_option("/timestepping/timestep", dt)
!        total_timestep=int((finish_time-current_time)/dt)+1
        u=>extract_vector_field(state(1), "Velocity")
        p=>extract_scalar_field(state(1), "Pressure")
        x => extract_vector_field(state(1), "Coordinate") 
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
          ! calculate the misfit between the data and numerical results, then project onto the reduced space
          allocate(misfit_POD_rhs((u%dim+1)*size(POD_state,1)))
          call misfit_POD(POD_state,misfit_POD_rhs)
!print*,misfit_POD_rhs
!print*,'****************'
!print*,pod_coef_all_obv(timestep,1:u%dim*size(POD_state,1))-pod_coef_all(timestep,1:u%dim*size(POD_state,1))
!stop 99
          pod_rhs_adjoint%val(1:u%dim*size(POD_state,1))= pod_rhs_adjoint%val(1:u%dim*size(POD_state,1))+ &
!               pod_coef_all_obv(timestep,1:u%dim*size(POD_state,1))-pod_coef_all(timestep,1:u%dim*size(POD_state,1))
               misfit_POD_rhs(1:u%dim*size(POD_state,1))

          deallocate(misfit_POD_rhs)
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
                call max_sort(size(optimal_location),u_total,optimal_location)
                open(10,file='optimal_sensor_location.dat')
                write(10,*) (optimal_location(i),i=1,size(optimal_location))
                close(10)  
                open(10,file='optimal_sensor_coordinate.dat')
                do i=1,2000 ! only save 200 optimal location
                   j=optimal_location(i)
                   write(10,*) x%val(1,j),x%val(2,j)
                enddo
                close(10)
                deallocate(optimal_location)
                deallocate(u_total)          
                stop 34
             endif
          endif

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
       deallocate(pod_coef_all_obv)

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

     subroutine misfit_POD(POD_state,misfit_POD_rhs)
       implicit none
       type(state_type), dimension(:,:,:), intent(inout) :: POD_state
       integer:: n
       real, dimension(:),intent(inout) :: misfit_POD_rhs
       real, allocatable, dimension(:) :: num, obs,diff
       type(vector_field), pointer :: pod_velocity,x
       type(scalar_field), pointer :: pod_pressure
       integer :: istate,nonods
       integer :: i,j,k,node
       type(vector_field), pointer :: snapmean_velocity
       type(scalar_field), pointer :: snapmean_pressure
       print*,'misfit_POD'
       istate =1
       pod_velocity=>extract_vector_field(POD_state(1,1,istate), "Velocity")
       x => extract_vector_field(POD_state(1,1,istate), "Coordinate") 
       nonods=node_count(pod_velocity)
       allocate(num( (pod_velocity%dim+1)*nonods))
       allocate(obs( (pod_velocity%dim+1)*nonods))
       allocate(diff( (pod_velocity%dim+1)*nonods))
       num=0.0
       obs=0.0
       diff=0.0

       open(10, file = 'num_last.dat')
       do i=1,pod_velocity%dim
          read(10,*) (num((i-1)*nonods+j),j=1,nonods)
       enddo
       read(10,*) (num(pod_velocity%dim*nonods+j),j=1,nonods)
       close(10) 

       open(20, file = 'obs_last.dat')
       do i=1,pod_velocity%dim
          read(20,*) (obs((i-1)*nonods+j),j=1,nonods)
       enddo
       read(20,*) (obs(pod_velocity%dim*nonods+j),j=1,nonods)
       close(20) 

       diff = obs-num
       do node = 1, nonods
!          if(abs(x%val(1,node)).gt.0.3.or.abs(x%val(1,node)).lt.0.05.or.  &
!               abs(x%val(2,node)).gt.0.9.or.abs(x%val(2,node)).lt.0.7) then
          if(abs(x%val(1,node)).gt.1.05.or.abs(x%val(1,node)).lt.0.79.or.  &
               abs(x%val(2,node)).gt.0.65.or.abs(x%val(2,node)).lt.0.47) then
      !       data(node) = 0.0
      !       num(node) = 0.0
             diff(node) = 0.0
          endif
       enddo
    do i=1,size(POD_state,1)
       pod_velocity=>extract_vector_field(POD_state(i,1,istate), "Velocity")
       snapmean_velocity=>extract_vector_field(POD_state(i,1,istate),"SnapmeanVelocity")
       snapmean_pressure=>extract_Scalar_field(POD_state(i,2,istate),"SnapmeanPressure") 
       do j=1,pod_velocity%dim
          misfit_POD_rhs(i+size(POD_state,1)*(j-1))=   &
               dot_product(POD_velocity%val(j,:), diff((j-1)*nonods+1:j*nonods) )
       enddo
       POD_pressure=>extract_scalar_field(POD_state(i,2,istate), "Pressure")
       misfit_POD_rhs(i+size(POD_state,1)*pod_velocity%dim)= &
            dot_product(POD_pressure%val(:), diff(pod_velocity%dim*nonods+1:(pod_velocity%dim+1)*nonods) )   
    enddo
    deallocate(num)
    deallocate(obs)
    deallocate(diff)
       
  end subroutine misfit_POD
  
end module momentum_equation_reduced
