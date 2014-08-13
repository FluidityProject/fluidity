
!The typs originally declared in the darcy impes assemble are noved to this seperat file 
!to avoid circular dependency between modules
!lcai 18 July 2014

module darcy_impes_assemble_type   

#include "fdebug.h"

  use spud
  use fields
  use state_module
  use fldebug
  use sparse_tools_petsc
  use cv_faces
  use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
  
  !!*****03 July 2014 LCai ************use Leaching Chemical model****!
  use darcy_impes_leaching_types


  implicit none
  private
  
  public :: darcy_impes_type, &
            darcy_impes_cv_options_type,&
            darcy_impes_relperm_corr_options_type,&
            darcy_impes_subcycle_options_type,&
            darcy_impes_adaptive_dt_options_type,&
            darcy_impes_eos_options_type,&
            cached_face_value_type,&
            darcy_impes_generic_prog_sfield_type

               
   
   ! Options associated with the relative permeability correlation
   type darcy_impes_relperm_corr_options_type
      ! The correlation type to use
      integer :: type
      ! The exponent for the power law correlation for each phase
      real, dimension(:), pointer :: exponents
      ! The residual saturation for each phase
      real, dimension(:), pointer :: residual_saturations
      ! The cut off value for saturation when calculating relperm
      real, dimension(:), pointer :: cutoff_saturations
      ! The scaling_coefficients when calculating relperm
      real, dimension(:), pointer :: scaling_coefficients      
   end type darcy_impes_relperm_corr_options_type
   
   ! Options associated with the CV discretisation for the darcy impes solver
   type darcy_impes_cv_options_type
      ! what CV_FACEVALUE_ to use
      integer :: facevalue
      ! the number of face value iterations
      integer :: number_face_value_iteration
      ! whether to limit the face value spatially
      logical :: limit_facevalue
      ! what CV_LIMITER_ to use
      integer :: limiter
      ! the slopes of the limiter being used (if Sweby)
      real, dimension(2) :: limiter_slopes
      ! what upwind scheme are we using (if any)
      integer :: upwind_scheme   
   end type darcy_impes_cv_options_type
   
   ! Options associated with explicit advection subcycling for CV scalar field solver
   type darcy_impes_subcycle_options_type
      logical :: have
      integer :: number
      logical :: consistent
   end type darcy_impes_subcycle_options_type
   
   ! Options associated with adaptive time stepping
   type darcy_impes_adaptive_dt_options_type
      logical :: have
      real    :: requested_cfl
      real    :: min_dt
      real    :: max_dt
      real    :: increase_tolerance
      logical :: min_dt_terminate_if_reached
      logical :: at_first_dt
   end type darcy_impes_adaptive_dt_options_type
   
   ! Options associated with EoS for a phase
   type darcy_impes_eos_options_type
      logical :: have_fluids_linear
      real    :: fluids_linear_reference_density
   end type darcy_impes_eos_options_type
   
   ! Data associated with the cached CV face values
   type cached_face_value_type
      real, dimension(:,:,:,:), pointer :: relperm      ! nsub, ngi_vele, vele, nphase
      real, dimension(:,:,:,:), pointer :: relperm_bdy  ! nsub, ngi_sele, sele, nphase
      real, dimension(:,:,:),   pointer :: den          !       ngi_vele, vele, nphase
      real, dimension(:,:,:),   pointer :: den_bdy      !       ngi_sele, sele, nphase
      real, dimension(:,:),     pointer :: detwei       !       ngi_vele, vele
      real, dimension(:,:),     pointer :: detwei_bdy   !       ngi_sele, sele
      real, dimension(:,:,:),   pointer :: normal       ! ndim, ngi_vele, vele
      real, dimension(:,:,:),   pointer :: normal_bdy   ! ndim, ngi_sele, sele
      real, dimension(:,:,:,:), pointer :: p_dshape     ! nloc, ngi_vele, ndim, vele
!!!!! *** THIS IS NOT POSSIBLE YET ***
!!!!!      real, dimension(:,:,:,:), pointer :: p_dshape_bdy ! nloc, ngi_sele, ndim, sele
      logical                           :: cached_detwei_normal
      logical                           :: cached_p_dshape
   end type cached_face_value_type
   
   ! Data associated with a generic prognostic sfield, some of which will point to fields in state
   type darcy_impes_generic_prog_sfield_type
      ! The phase this sfield is associated with
      integer :: phase
      ! LCai 08 Aug 2013***********************
      character(len=FIELD_NAME_LEN) :: source_name !the name of the source immobile filed
      !****************Finish*******************
      type(scalar_field), pointer :: sfield
      type(scalar_field), pointer :: old_sfield
      type(scalar_field), pointer :: sfield_abs
      type(scalar_field), pointer :: sfield_src
      type(tensor_field), pointer :: sfield_diff
      logical :: have_diff
      logical :: have_abs
      logical :: have_src
      logical :: have_adv
      logical :: have_MIM_source !LCai 08 Aug 2013
      type(darcy_impes_cv_options_type) :: sfield_cv_options
      !******For the Leaching chemical source term***LCai** 30 June 2014*****!
      type(leach_chemical_prog_sfield_src) :: lc_src
      !*****Finish*******

   end type darcy_impes_generic_prog_sfield_type
  


   !***********************************LCai 08 Aug 2013**************
   type immobile_prog_sfield_type
      integer :: phase !the phase of this field
      character(len=FIELD_NAME_LEN) :: source_name !the name of the source mobile filed
      type(scalar_field), pointer :: sfield
      type(scalar_field), pointer :: old_sfield
   end type immobile_prog_sfield_type
   
   
   !***********************************LCai 23 July & 08 & 16 Aug 2013**************
   !the Mobile-Immobile model options
   type darcy_impes_MIM_options_type
      type(scalar_field_pointer), dimension(:), pointer :: immobile_saturation
      type(scalar_field_pointer), dimension(:), pointer :: old_immobile_saturation
      type(scalar_field_pointer), dimension(:), pointer :: mobile_saturation
      type(scalar_field_pointer), dimension(:), pointer :: old_mobile_saturation
      type(scalar_field_pointer), dimension(:), pointer :: mass_trans_coef
      type(scalar_field_pointer), dimension(:), pointer :: old_mass_trans_coef
      type(immobile_prog_sfield_type), dimension(:), pointer :: immobile_prog_sfield
      type(scalar_field) :: MIM_src !the source term of MIM added to the matrix 
      type(scalar_field) :: MIM_src_s !the second source term of MIM added to the matrix
      ! *** Flag for Whether there is Mobile-Immobile model
      logical, dimension(:), pointer :: have_MIM
      ! *** Flag to check wether the MIM exist in at least one phase
      logical :: have_MIM_phase
      ! *** Flag for Whether there is Mass transfer coefficient
      logical, dimension(:), pointer :: have_mass_trans_coef
      ! *** Flag for wether there are immobile prognostic fields
      !***NOT USED**logical, dimension(:), pointer :: have_immobile_prog_sfield
   end type darcy_impes_MIM_options_type
   !***********Finish******************LCai ***********************************
   
   
   type darcy_impes_type
      ! *** Pointers to fields from state that have array length of number of phases ***
      type(vector_field_pointer), dimension(:), pointer :: darcy_velocity
      type(scalar_field_pointer), dimension(:), pointer :: mobility      
      type(scalar_field_pointer), dimension(:), pointer :: fractional_flow
      type(scalar_field_pointer), dimension(:), pointer :: saturation
      type(scalar_field_pointer), dimension(:), pointer :: old_saturation
      type(scalar_field_pointer), dimension(:), pointer :: saturation_source
      type(scalar_field_pointer), dimension(:), pointer :: relative_permeability
      type(scalar_field_pointer), dimension(:), pointer :: old_relative_permeability
      type(scalar_field_pointer), dimension(:), pointer :: viscosity
      type(scalar_field_pointer), dimension(:), pointer :: cfl
      type(scalar_field_pointer), dimension(:), pointer :: pressure
      type(scalar_field_pointer), dimension(:), pointer :: capilliary_pressure
      type(scalar_field_pointer), dimension(:), pointer :: density
      type(scalar_field_pointer), dimension(:), pointer :: old_density

      ! *** Pointers to fields from state that are NOT phase dependent ***
      type(mesh_type),    pointer :: pressure_mesh
      type(mesh_type),    pointer :: elementwise_mesh
      type(scalar_field), pointer :: average_pressure
      type(scalar_field), pointer :: porosity
      type(scalar_field), pointer :: old_porosity
      type(scalar_field), pointer :: absolute_permeability
      type(vector_field), pointer :: positions
      type(vector_field), pointer :: total_darcy_velocity
      type(scalar_field), pointer :: total_mobility
      type(scalar_field), pointer :: sum_saturation
      type(scalar_field), pointer :: div_total_darcy_velocity
      type(vector_field), pointer :: gravity_direction
      type(vector_field), pointer :: bulk_darcy_velocity
      ! ***** Pointers to DUAL fields from state that have array length of number of phases *****
      type(scalar_field_pointer), dimension(:), pointer :: pressure_other_porous_media
      type(scalar_field_pointer), dimension(:), pointer :: transmissibility_lambda_dual
      ! *** Pointer to the pressure mesh - pressure mesh sparsity, used for pressure matrix and finding CV upwind values ***
      type(csr_sparsity), pointer :: sparsity_pmesh_pmesh
      ! *** Data associated with generic prognostic scalar fields, some of which points to fields in state ***
      type(darcy_impes_generic_prog_sfield_type), dimension(:), pointer :: generic_prog_sfield
      ! *** Data associated with GradientPressure fields ***
      type(element_type)                                :: gradient_pressure_shape
      type(mesh_type)                                   :: gradient_pressure_mesh
      type(vector_field_pointer), dimension(:), pointer :: gradient_pressure
      type(vector_field_pointer), dimension(:), pointer :: iterated_gradient_pressure
      ! *** Data associated with gravity ***
      logical            :: have_gravity
      real               :: gravity_magnitude
      type(vector_field) :: gravity
      ! *** Fields allocated here used in assemble algorithm ***
      type(vector_field)     :: positions_pressure_mesh
      type(csr_matrix)       :: matrix
      type(petsc_csr_matrix) :: dual_block_pressure_matrix
      type(csr_matrix)       :: pressure_matrix
      type(scalar_field), pointer :: pressure_rhs
      type(scalar_field)     :: lhs
      type(scalar_field)     :: rhs
      type(scalar_field)     :: rhs_full
      type(scalar_field)     :: rhs_high_resolution
      type(scalar_field)     :: inverse_cv_mass_pressure_mesh
      type(scalar_field)     :: cv_mass_pressure_mesh_with_source
      type(scalar_field)     :: cv_mass_pressure_mesh_with_porosity   
      type(scalar_field)     :: cv_mass_pressure_mesh_with_old_porosity 
      type(scalar_field)     :: inverse_cv_sa_pressure_mesh
      type(scalar_field)     :: work_array_of_size_pressure_mesh
      type(scalar_field)     :: modified_relative_permeability
      type(scalar_field), pointer :: constant_zero_sfield_pmesh
      type(scalar_field), pointer :: constant_zero_sfield_elementwisemesh
      type(scalar_field_pointer), dimension(:), pointer :: old_saturation_subcycle
      ! ***** coupling field associated with DUAL model *****
      type(scalar_field) :: cv_mass_pressure_mesh_with_lambda_dual
      ! ***** rhs src field associated with DUAL model *****
      type(scalar_field) :: rhs_dual
      ! *** Data associated with v, pressure and sfield BC allocated here ***
      type(mesh_type)                          :: bc_surface_mesh
      type(scalar_field)                       :: v_bc_value
      integer,           dimension(:), pointer :: v_bc_flag
      type(scalar_field)                       :: sfield_bc_value
      integer,           dimension(:), pointer :: sfield_bc_flag
      type(scalar_field)                       :: pressure_bc_value
      integer,           dimension(:), pointer :: pressure_bc_flag
      type(scalar_field)                       :: inverse_characteristic_length
      real                                     :: weak_pressure_bc_coeff
      ! *** The number of phase and the CV surface quadrature degree to use ***
      integer :: number_phase, cv_surface_quaddegree   
      ! *** Data specifically associated with the CV discretisation allocated here ***
      type(darcy_impes_cv_options_type)       :: saturation_cv_options
      type(darcy_impes_cv_options_type)       :: relperm_cv_options
      type(darcy_impes_cv_options_type)       :: density_cv_options
      type(cv_faces_type)                     :: cvfaces
      type(element_type)                      :: x_cvshape_full
      type(element_type)                      :: p_cvshape_full
      type(element_type)                      :: gradp_cvshape_full
      type(element_type)                      :: x_cvshape
      type(element_type)                      :: p_cvshape
      type(element_type)                      :: gradp_cvshape
      type(element_type)                      :: x_cvbdyshape
      type(element_type)                      :: p_cvbdyshape      
      type(element_type)                      :: gradp_cvbdyshape
      type(darcy_impes_subcycle_options_type) :: subcy_opt_sat
      type(csr_matrix)                        :: relperm_upwind
      type(csr_matrix)                        :: sfield_upwind
      ! *** The minimum value of saturation to use in the denominator of modrelperm face value ***
      real :: minimum_denominator_saturation_value
      ! *** A flag for whether saturation face values need to be determined
      logical :: determine_saturation_face_values
      ! *** Data associate with the relperm correlation options ***
      type(darcy_impes_relperm_corr_options_type) :: relperm_corr_options
      
      !***********************************LCai 23 & 27 July & 22 Aug 2013***********************
      type (darcy_impes_MIM_options_type) :: MIM_options
      type(scalar_field), pointer  :: sat_ADE !saturation used for adv-diff equation for solving prog sfield
      type(scalar_field), pointer  :: old_sat_ADE 

      type(scalar_field) :: porosity_pmesh ! the porosity based on pressure mesh
      type(scalar_field) :: old_porosity_pmesh
      real, dimension(1) :: porosity_cnt !constant value for porosity as a constant
      real, dimension(1) :: old_porosity_cnt
      logical :: prt_is_constant !is the porosity a constant
     
      !***********Finish******************LCai **********************************************
      
      ! *** Flag for each phase for whether there is capilliary pressure ***
      logical, dimension(:), pointer :: have_capilliary_pressure
      ! *** Flag for each phase for whether there is a saturation source ***
      logical, dimension(:), pointer :: have_saturation_source
      ! *** Flag for whether the first phase saturation is diagnostic, else it is prognostic ***
      logical :: phase_one_saturation_diagnostic      
      ! *** Flag for whether the first phase pressure is prognostic, else it is prescribed *** 
      logical :: first_phase_pressure_prognostic

      ! *** The cached face values for all phases ***
      type(cached_face_value_type) :: cached_face_value
      ! *** The advection subcycle timestep size ***
      real    :: dt_subcycle
      ! *** Time time step size, also stored here for convenience ***
      real :: dt
      ! *** The current time, also stored here for convenience ***
      real :: current_time 
      ! *** Non linear iteration, also stored here for convencience ***
      integer :: nonlinear_iter
      ! *** Max Non linear iteration for this timestep, also stored here for convencience ***
      integer :: max_nonlinear_iter_this_timestep
      ! *** The number of volume, surface elements and nodes, stored here for convenience ***
      integer :: number_vele
      integer :: number_sele
      integer :: number_pmesh_node
      ! *** Geometric dimension, also stored here for convenience ***
      integer :: ndim
      ! *** Options data associated with adaptive time stepping stored here ***
      type(darcy_impes_adaptive_dt_options_type) :: adaptive_dt_options
      ! *** Options data associated with each phase EoS stored here *** 
      type(darcy_impes_eos_options_type), dimension(:), pointer :: eos_options
      ! *** Pointer to main state array, for convenience ***
      type(state_type), dimension(:), pointer :: state
      !*******03 July 2014 Lcai*****Leaching chemical model*************!
      type(leach_chemical_type) :: lc
      !********03 July 2014 Lcai*****Finish****************************!

   end type darcy_impes_type

end module
