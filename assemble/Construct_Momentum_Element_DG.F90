#include "fdebug.h"  
#undef ASSERT
#ifdef NDEBUG
#define ASSERT(X)
#else
#ifdef FORTRAN_DISALLOWS_LONG_LINES
#define ASSERT(X) IF(.NOT.(X)) FLAbort("Failed assertion ")
#else
#define ASSERT(X) IF(.NOT.(X)) FLAbort("Failed assertion X")
#endif
#endif

module momentum_element_dg

  use spud
  use fldebug
  use vector_tools
  use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN, COLOURING_DG2, &
       COLOURING_DG0
#ifdef _OPENMP
  use omp_lib
#endif
  use integer_set_module
  use parallel_tools
  use sparse_tools
  use shape_functions
  use transform_elements
  use fetools
  use parallel_fields
  use fields
  use profiler
  use petsc_tools
  use sparse_tools_petsc
  use sparse_matrices_fields
  use state_module
  use vtk_interfaces
  use halos
  use field_options
  use fefields
  use boundary_conditions, only: has_boundary_condition, get_entire_boundary_condition
  use field_derivatives
  use coordinates
  use solvers
  use sparsity_patterns
  use dgtools
  use smoothing_module
  use sparsity_patterns_meshes
  use boundary_conditions_from_options
  use coriolis_module, only : coriolis, set_coriolis_parameters
  use turbine
  use diagnostic_fields
  use slope_limiters_dg
  use colouring
  use multiphase_module

  ! Module private variables for model options. This prevents us having to
  ! do dictionary lookups for every element (or element face!)
  real :: dt, theta, theta_nl
  logical :: lump_mass, lump_abs, lump_source, subcycle

  ! Whether the advection term is only integrated by parts once.
  logical :: integrate_by_parts_once=.false.
  ! Whether the conservation term is integrated by parts or not
  logical :: integrate_conservation_term_by_parts=.false.
  ! Whether or not to integrate the surface tension term by parts
  logical :: integrate_surfacetension_by_parts

  ! Weight between conservative and non-conservative forms of the advection
  ! equation. 
  ! 1 is for conservative 0 is for non-conservative.
  real :: beta

  ! Discretisation to use for viscosity term.
  integer :: viscosity_scheme
  integer, parameter :: ARBITRARY_UPWIND=1
  integer, parameter :: BASSI_REBAY=2
  integer, parameter :: CDG=3
  integer, parameter :: IP=4

  ! Method for getting h0 in IP
  integer :: edge_length_option
  integer, parameter :: USE_FACE_INTEGRALS=1
  integer, parameter :: USE_ELEMENT_CENTRES=2

  ! Parameters for interior penalty method
  real :: Interior_Penalty_Parameter, edge_length_power, h0

  ! Flag indicating whether to include pressure bcs (not for cv pressure)
  logical :: l_include_pressure_bcs
  
  ! which terms do we have?
  logical :: have_mass
  logical :: have_source
  logical :: have_gravity
  logical :: on_sphere, radial_gravity
  logical :: have_absorption
  logical :: have_vertical_stabilization
  logical :: have_implicit_buoyancy
  logical :: have_vertical_velocity_relaxation
  logical :: have_swe_bottom_drag
  ! implicit absorption is corrected by the pressure correction
  ! by combining the implicit part of absorption with the mass term of u^{n+1}
  logical :: pressure_corrected_absorption
  logical :: have_viscosity
  logical :: have_surfacetension
  logical :: have_coriolis
  logical :: have_advection
  logical :: move_mesh
  logical :: have_pressure_bc
  logical :: subtract_out_reference_profile
  logical :: have_les
  
  real :: gravity_magnitude

  ! CDG stuff
  real, dimension(3) :: switch_g
  logical :: CDG_penalty
  logical :: remove_penalty_fluxes

  ! Are we running a multi-phase flow simulation?
  logical :: multiphase

interface
     subroutine construct_momentum_element_dg_interface(ele, big_m, rhs, &
       &X, U, U_nl, U_mesh, X_old, X_new, Source, Buoyancy, hb_density, hb_pressure, gravity, Abs, &
       &Viscosity, swe_bottom_drag, swe_u_nl, P, old_pressure, Rho, surfacetension, q_mesh, &
       &velocity_bc, velocity_bc_type, &
       &pressure_bc, pressure_bc_type, &
       &turbine_conn_mesh, depth, have_wd_abs, alpha_u_field, Abs_wd, &
       &vvr_sf, ib_min_grad, nvfrac, &
       &inverse_mass, inverse_masslump, mass, subcycle_m, subcycle_rhs, partial_stress, &
       &smagorinsky_coefficient, eddy_visc, prescribed_filter_width, distance_to_wall, &
       &y_plus_debug, les_filter_width_debug)
       use sparse_tools
       use fields
       use sparse_tools_petsc
    !!< Construct the momentum equation for discontinuous elements in
    !!< acceleration form.
       implicit none
       !! Index of current element
       integer :: ele
       !! Main momentum matrix.
       type(petsc_csr_matrix), intent(inout) :: big_m
       !! Momentum right hand side vector for each point.
       type(vector_field), intent(inout) :: rhs
       !! Auxiliary variable mesh
       type(mesh_type), intent(in) :: q_mesh
       type(mesh_type), intent(in) :: turbine_conn_mesh
       type(scalar_field), intent(in) :: depth
       !! 

       real, intent(in) :: vvr_sf, ib_min_grad

       type(block_csr_matrix), intent(inout), optional :: subcycle_m
       type(vector_field), intent(inout), optional :: subcycle_rhs
       
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

       type(scalar_field), intent(in) :: nvfrac

       ! added for partial stress form (sp911)
       logical, intent(in) :: partial_stress

       ! LES - sp911
       real, intent(in) :: smagorinsky_coefficient
       type(scalar_field), pointer, intent(inout) :: eddy_visc, y_plus_debug, &
            & les_filter_width_debug
       type(scalar_field), pointer, intent(in) :: prescribed_filter_width, distance_to_wall
       

     end subroutine construct_momentum_element_dg_interface
  end interface

  interface
     subroutine construct_momentum_interface_dg_interface(ele, face, face_2, ni, &
       & big_m_tensor_addto, &
       & rhs_addto, Grad_U_mat, Div_U_mat, X, Rho, U,&
       & U_nl, U_mesh, P, q_mesh, surfacetension, &
       & velocity_bc, velocity_bc_type, &
       & pressure_bc, pressure_bc_type, hb_pressure, &
       & subcycle_m_tensor_addto, subcycle_rhs_addto, nvfrac, &
       & ele2grad_mat, kappa_mat, inverse_mass_mat, &
       & viscosity, viscosity_mat, partial_stress )
       use fields, only: scalar_field, vector_field, tensor_field, mesh_type
       implicit none
       integer, intent(in) :: ele, face, face_2, ni
       real, dimension(:,:,:,:), intent(inout) :: big_m_tensor_addto
       real, dimension(:,:,:,:), intent(inout) :: subcycle_m_tensor_addto
       real, dimension(:,:), intent(inout) :: rhs_addto, subcycle_rhs_addto
       real, dimension(:,:,:), intent(inout) :: Grad_U_mat, Div_U_mat
       ! We pass these additional fields to save on state lookups.
       type(vector_field), intent(in) :: X, U, U_nl
       type(vector_field), pointer :: U_mesh
       type(scalar_field), intent(in) :: Rho, P
       type(scalar_field), intent(in) :: nvfrac
       !! Mesh of the auxiliary variable in the second order operator.
       type(mesh_type), intent(in) :: q_mesh
       !! surfacetension
       type(tensor_field), intent(in) :: surfacetension
       !! Boundary conditions associated with this interface (if any).
       type(vector_field), intent(in) :: velocity_bc
       integer, dimension(:,:), intent(in) :: velocity_bc_type
       type(scalar_field), intent(in) :: pressure_bc
       integer, dimension(:), intent(in) :: pressure_bc_type
       type(scalar_field), intent(in) :: hb_pressure
       !! Computation of primal fluxes and penalty fluxes
       real, intent(in), optional, dimension(:,:,:) :: ele2grad_mat
       !! \Int_{ele} N_i kappa N_j dV, used for CDG fluxes
       real, dimension(:,:,:,:), intent(in), optional :: kappa_mat
       !! Inverse element mass matrix.
       real, dimension(:,:), intent(in), optional :: inverse_mass_mat
       type(tensor_field), intent(in), optional :: viscosity
       !! Local viscosity matrix for assembly.
       real, intent(inout), dimension(:,:,:,:), optional :: viscosity_mat
       logical, intent(in) :: partial_stress

     end subroutine construct_momentum_interface_dg_interface
  end interface

  procedure(construct_momentum_element_dg_interface), pointer :: construct_momentum_element_dg
  procedure(construct_momentum_interface_dg_interface), pointer :: construct_momentum_interface_dg

  contains



!!! Two dimensional

#define NDIM(variable) 2
#define MDIM(variable) 2
#define NLOC_X 3

!!! P1DGP2

#define NLOC(variable,ele) 3
#define FLOC(variable,face) 2
#define EFLOC(variable,face) 9
#define P_FLOC 3

#define NGI(variable,ele) 6
#define FNGI(variable,face) 4
#define WRAP_NAME(name) name##_P1DGP2_2D_6GI
#include "Construct_Momentum_Element_DG_template.F90"
#undef NGI
#undef FNGI
#undef WRAP_NAME
#define NGI(variable,ele) 16
#define FNGI(variable,face) 8
#define WRAP_NAME(name) name##_P1DGP2_2D_16GI
#include "Construct_Momentum_Element_DG_template.F90"
#undef NGI
#undef FNGI
#undef WRAP_NAME


#undef NLOC
#undef FLOC
#undef EFLOC
#undef P_FLOC




#undef NDIM
#undef MDIM
#undef NDIM_X
#undef NLOC_X

!!! Three dimensional


#define NDIM(variable) 3
#define MDIM(variable) 3
#define NLOC_X 4

!!! P1DGP2

#define NLOC(variable,ele) 4
#define FLOC(variable,face) 3

#define EFLOC(variable,face) 16
#define P_FLOC 6
#define NGI(variable,ele) 11
#define FNGI(variable,face) 6

#define WRAP_NAME(name) name##_P1DGP2_3D_11GI
#include "Construct_Momentum_Element_DG_template.F90"

#undef NGI
#undef FNGI
#undef WRAP_NAME

#undef NLOC
#undef FLOC
#undef EFLOC
#undef P_FLOC

#undef NDIM
#undef MDIM
#undef NDIM_X
#undef NLOC_X

!!! generic

#define GENERIC
#define NDIM(variable) field_dim(variable)
#define MDIM(variable) mesh_dim(variable)
#define NLOC(variable,ele) ele_loc(variable,ele)
#define NGI(variable,ele) ele_ngi(U,ele)
#define FNGI(variable,face) face_ngi(variable,face) 
#define FLOC(variable,face) face_loc(U,face)
#define EFLOC(variable,face) ele_and_faces_loc(U,ele)
#define P_FLOC face_loc(P,face)
#define NLOC_X ele_loc(X,ele)
#define WRAP_NAME(name) name##_generic
#include "Construct_Momentum_Element_DG_template.F90"
#undef NDIM
#undef NDIM_X
#undef NLOC
#undef NLOC_X
#undef NGI
#undef FNGI
#undef FLOC
#undef EFLOC
#undef P_FLOC
#undef GENERIC

end module momentum_element_dg
