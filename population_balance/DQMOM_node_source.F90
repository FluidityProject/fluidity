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
!    but WITHOUT ANY WARRANTY; without even the implied arranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module dqmom

  use fldebug
  use vector_tools
  use spud
  use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
  use futils, only: int2str
  use sparse_tools
  use fields
  use state_module
  use state_fields_module
  use initialise_fields_module
  use solvers

  implicit none

  public dqmom_init, dqmom_calculate_source_terms, dqmom_calculate_abscissa,&
       & dqmom_check_options, dqmom_calculate_moments, dqmom_calculate_statistics, dqmom_apply_min_weight

  private

  real, save     :: fields_min = 1.0e-11
  !! TODO:
  !! 1. Make the algorithm work for multiple phases and make sure it works for
  !!    several pop_balances - currently works for only one population balance 
  !! 2. Make check_options run
  !! 3. Check all prognostic fields are identical (possible bar initial conditions)
  !! 4. Check prognostic Source terms are set to diagnostic, Internal
  !! 5. Check diagnostic fields are set to Internal

contains

  subroutine get_pop_field(state, i_pop, i_field, type, item, stat, iterated)

    ! Returns the required population balance field

    type(state_type), intent(in) :: state
    integer, intent(in) :: i_pop, i_field
    character(len=FIELD_NAME_LEN) :: type ! moments, weights, ascissa, or weighted_abscissa
    type(scalar_field), pointer, intent(out) :: item
    integer, intent(out), optional :: stat
    logical, intent(in), optional :: iterated
    
    character(len=FIELD_NAME_LEN) :: name
    character(len=OPTION_PATH_LEN) :: option_path 

    call get_pop_option_path(state, i_pop, option_path)
    call get_option(trim(option_path)//'/'//trim(type)//'/scalar_field['//int2str(i_field -&
         1)//']/name', name) 
    if (present(iterated)) then
       if (iterated) then
          name = 'Iterated'//trim(name)
       end if
    end if
    item => extract_scalar_field(state, trim(name), stat)
  end subroutine get_pop_field

  subroutine get_pop_option_path(state, i_pop, option_path)

    type(state_type), intent(in) :: state
    integer, intent(in) :: i_pop
    character(len=OPTION_PATH_LEN), intent(out) :: option_path 

    option_path = trim(state%option_path)//'/population_balance'
    if (option_count(option_path) > 1) then
       option_path = trim(state%option_path)//'/population_balance['&
            &//int2str(i_pop-1)//']'
    end if

  end subroutine get_pop_option_path

  subroutine dqmom_init(states)

    type(state_type), intent(in), dimension(:) :: states

    integer :: i_state, i_pop
    character(len=OPTION_PATH_LEN) :: option_path 

    do i_state = 1, option_count("/material_phase")
       do i_pop = 1, option_count(trim(states(i_state)%option_path)//&
            '/population_balance')
          call get_pop_option_path(states(i_state), i_pop, option_path)
          if (have_option(trim(option_path)// &
               '/calculate_initial_conditions_from_moments')) then
             call dqmom_PD_algorithm(states(i_state), i_pop)
          end if
       end do
    end do

    call dqmom_calculate_abscissa(states)

  end subroutine dqmom_init

  subroutine dqmom_PD_algorithm(state, i_pop)

    ! Algorithm for calculating abscissa and weights from the moment of a distribution
    ! See Gordon 1968 and Mcgraw 1997
   
    type(state_type), intent(in) :: state
    integer, intent(in) :: i_pop

    type(scalar_field_pointer), dimension(:), allocatable :: moments
    type(scalar_field_pointer), dimension(:), allocatable :: weighted_abscissa
    type(scalar_field_pointer), dimension(:), allocatable :: weights
    type(vector_field), pointer:: position

    ! J - Jacobi matrix
    ! P and alpha are required working arrays
    ! e_values - eigenvalues
    ! e_vectors - corresponding eigenvectors
    real, dimension(:,:), allocatable :: P, Jac, e_vectors
    real, dimension(:), allocatable :: alpha, e_values
    character(len=OPTION_PATH_LEN) :: option_path
    character(len=FIELD_NAME_LEN) :: type
    integer :: i_field, n_abscissa, i_node, i, j, stat

    call get_pop_option_path(state, i_pop, option_path)
    n_abscissa = option_count(trim(option_path)//'/abscissa/scalar_field')

    allocate(moments(2*n_abscissa))
    allocate(weighted_abscissa(n_abscissa))
    allocate(weights(n_abscissa))
    allocate(P(2*n_abscissa + 1,2*n_abscissa + 1))
    allocate(Jac(n_abscissa,n_abscissa))
    allocate(alpha(2*n_abscissa))
    allocate(e_vectors(n_abscissa,n_abscissa))
    allocate(e_values(n_abscissa))

    position => extract_vector_field(state, "Coordinate")
    do i_field = 1, 2*n_abscissa
       ! collect moment fields and initialise from initial conditions
       ! will not be done automatically as this is a diagnostic field
       type = 'moments'
       call get_pop_field(state, i_pop, i_field, type, moments(i_field)%ptr)
       call zero(moments(i_field)%ptr)
       call initialise_field_over_regions(moments(i_field)%ptr, &
            trim(moments(i_field)%ptr%option_path)//'/diagnostic/initial_condition', position)
    end do
    
    do i_field = 1, n_abscissa
       ! collect weighted_abscissa and weight fields and zero
       type = 'weights'
       call get_pop_field(state, i_pop, i_field, type, weights(i_field)%ptr)
       type = 'weighted_abscissa'
       call get_pop_field(state, i_pop, i_field, type, &
            weighted_abscissa(i_field)%ptr)
       call zero(weights(i_field)%ptr)
       call zero(weighted_abscissa(i_field)%ptr)
    end do

    do i_node = 1, node_count(moments(1)%ptr)
       ! zero working arrays
       P = 0.0
       Jac = 0.0
       alpha = 0.0
       e_vectors = 0.0
       e_values = 0.0
       
       ! Construct P matrix
       P(1,1) = 1.0
       ! set the zero'th moment to 1.0 for the purposes of finding the abscissa
       ! weights will be multiplied by the zero'th moment later to obtain their correct
       ! values
       P(1,2) = 1.0       
       do i = 2, 2*n_abscissa
          P(i,2) = (-1)**(i-1)*(node_val(moments(i)%ptr, i_node)/node_val(moments(1)%ptr, i_node))
       end do
       do j = 3, 2*n_abscissa + 1
          do i = 1, (2*n_abscissa + 1) - (j - 1)
             P(i,j) = P(1,j-1)*P(i+1,j-2) - P(1,j-2)*P(i+1,j-1)
          end do
       end do
       ! Construct alpha array
       alpha(1) = 0.0
       do i = 2, 2*n_abscissa
          alpha(i) = P(1,i+1)/(P(1,i)*P(1,i-1))
       end do
       ! Construct jacobi matrix
       do i = 1, n_abscissa
          Jac(i,i) = alpha(2*i) + alpha(2*i - 1)
       end do
       do i = 1, n_abscissa - 1
          Jac(i, i+1) = ((alpha(2*i + 1)*alpha(2*i))**2.0)**0.25
          Jac(i+1, i) = ((alpha(2*i + 1)*alpha(2*i))**2.0)**0.25
       end do
       
       ! Calculate eigenvalues and eigenvectors
       call eigendecomposition_symmetric(Jac, e_vectors, e_values, stat)
       if (stat /= 0) then
          FLExit('Cannot compute abscissa and weights using PD algorithm')
       end if

       do i_field = 1, n_abscissa
          ! set prognostic field values
          call set(weights(i_field)%ptr, i_node, &
               node_val(moments(1)%ptr, i_node) * e_vectors(1,i_field)**2)
          call set(weighted_abscissa(i_field)%ptr, i_node, &
               node_val(weights(i_field)%ptr, i_node) * e_values(i_field))
       end do
       
    end do

    deallocate(moments, weighted_abscissa, weights, P, Jac, alpha, e_vectors, e_values)
    
  end subroutine dqmom_PD_algorithm

  subroutine dqmom_calculate_abscissa(states)

    type(state_type), dimension(:), intent(in) :: states

    type(scalar_field), pointer :: abscissa, weight, weighted_abscissa
    type(scalar_field) :: inv_weight
    integer :: i_state, i_pop, i_abscissa
    character(len=OPTION_PATH_LEN) :: option_path 
    character(len=FIELD_NAME_LEN) :: type
    
    ewrite(1, *) "In dqmom_calculate_abscissa" 
    do i_state = 1, option_count("/material_phase")
       do i_pop = 1, option_count(trim(states(i_state)%option_path)//'/population_balance')
          call get_pop_option_path(states(i_state), i_pop, option_path)
          do i_abscissa = 1, option_count(trim(option_path)//'/abscissa/scalar_field')
             ! get required fields
             type = 'abscissa'
             call get_pop_field(states(i_state), i_pop, i_abscissa, type, abscissa)
             type = 'weights'
             call get_pop_field(states(i_state), i_pop, i_abscissa, type, weight)
             type = 'weighted_abscissa'
             call get_pop_field(states(i_state), i_pop, i_abscissa, type, weighted_abscissa)

             ! calculate the inverse of the weights
             call allocate(inv_weight, weight%mesh, 'InverseWeight')
             call set(inv_weight, weight)
             where (inv_weight%val/=0.0)
                inv_weight%val=1./inv_weight%val
             end where

             ! calculate abscissa
             call set(abscissa, weighted_abscissa)
             call scale(abscissa, inv_weight)

             call deallocate(inv_weight)
          end do
       end do
    end do
    ewrite(1, *) "Exiting dqmom_calculate_abscissa"
  end subroutine dqmom_calculate_abscissa

  subroutine dqmom_apply_min_weight(states)

    type(state_type), dimension(:), intent(in) :: states

    type(scalar_field), pointer :: weight, weighted_abscissa
    integer :: i_state, i_pop, i_abscissa
    real :: minvalue
    character(len=OPTION_PATH_LEN) :: option_path
    character(len=FIELD_NAME_LEN) :: type

    ewrite(1, *) "In dqmom_apply_min_weight"
    do i_state = 1, option_count("/material_phase")
       do i_pop = 1, option_count(trim(states(i_state)%option_path)//'/population_balance')
          call get_pop_option_path(states(i_state), i_pop, option_path)
          do i_abscissa = 1, option_count(trim(option_path)//'/abscissa/scalar_field')
             ! get required fields
             type = 'weights'
             call get_pop_field(states(i_state), i_pop, i_abscissa, type, weight)
             type = 'weighted_abscissa'
             call get_pop_field(states(i_state), i_pop, i_abscissa, type, weighted_abscissa)

             ! Get minimum weight value from the diamond file
             call get_option(trim(option_path)//'/minimum_weight', minvalue)
             where (weight%val.le.minvalue)
                weight%val=minvalue
             end where

             if (have_option(trim(option_path)//'/minimum_weighted_abscissa')) then
                call get_option(trim(option_path)//'/minimum_weighted_abscissa', minvalue)
                where (weighted_abscissa%val.le.minvalue)
                   weighted_abscissa%val=minvalue
                end where
             endif
          end do
       end do
    end do
    ewrite(1, *) "Exiting dqmom_apply_min_weight"
  end subroutine dqmom_apply_min_weight

  subroutine dqmom_calculate_source_terms(states, it)

    type(state_type), dimension(:), intent(inout) :: states
    integer, intent(in) :: it    

    integer :: i_state, i_pop, cont_state = -1
    character(len=FIELD_NAME_LEN) :: cont_state_name

    ! find the continuous state number
    if (have_option("/population_balance_continuous_phase_name")) then
       call get_option("/population_balance_continuous_phase_name", cont_state_name)
       do i_state = 1, option_count("/material_phase")
          if (states(i_state)%name == cont_state_name) then
             cont_state = i_state
             exit
          end if
       end do
       if (cont_state <0) then
          FLAbort("Continuous state name you mentioned could not be located in the population balance calculations.")
       end if
    end if

    do i_state = 1, option_count("/material_phase")
       do i_pop = 1, option_count(trim(states(i_state)%option_path)//'/population_balance')
          call dqmom_calculate_source_term_pop(states(i_state), it, i_pop, states, cont_state)
       end do
    end do

  end subroutine dqmom_calculate_source_terms

  subroutine dqmom_calculate_source_term_pop(state, it, i_pop, states, cont_state)

    type(state_type), intent(inout) :: state
    type(state_type), dimension(:), intent(inout) :: states
    integer, intent(in) :: it
    integer, intent(in) :: cont_state

    type(scalar_field_pointer), dimension(:), allocatable :: abscissa,&
         weight, it_abscissa, it_weight, weighted_abscissa, s_weighted_abscissa, s_weight, a_weighted_abscissa, a_weight
    type(scalar_field), pointer :: lumped_mass, turbulent_dissipation, sponge_field
    type(tensor_field), pointer :: viscosity_continuous
    type(csr_matrix), pointer :: mass_matrix
    type(scalar_field), dimension(:), allocatable :: r_abscissa, r_weight
    type(tensor_field), pointer :: D
    type(vector_field), pointer :: X
    type(scalar_field) :: dummy_scalar
    real :: theta, cond, growth_r, internal_dispersion_coeff, aggregation_freq_const, breakage_freq_const, breakage_freq_degree, perturb_val, C5
    integer :: i_pop, N, i, j, stat, i_node
    character(len=OPTION_PATH_LEN) :: option_path 
    character(len=FIELD_NAME_LEN) :: type, field_name, growth_type, aggregation_freq_type, breakage_freq_type, breakage_dist_type, singular_option
    logical :: have_D = .false.
    logical :: have_growth = .FALSE.
    logical :: have_internal_dispersion = .FALSE.
    logical :: have_aggregation = .FALSE.
    logical :: have_breakage = .FALSE.

    call get_pop_option_path(state, i_pop, option_path)
    N = option_count(trim(option_path)//'/abscissa/scalar_field')
    allocate(abscissa(N), weight(N), it_abscissa(N), it_weight(N), &
         r_abscissa(N), r_weight(N), weighted_abscissa(N), s_weighted_abscissa(N), s_weight(N), a_weighted_abscissa(N), a_weight(N))

    do i = 1, N
       ! collect abscissa and weight fields
       type = 'weights'
       call get_pop_field(state, i_pop, i, type, weight(i)%ptr)
       call get_pop_field(state, i_pop, i, type, it_weight(i)%ptr, iterated=.true.)
       type = 'abscissa'
       call get_pop_field(state, i_pop, i, type, abscissa(i)%ptr)
       call get_pop_field(state, i_pop, i, type, it_abscissa(i)%ptr, iterated=.true.)
       type = 'weighted_abscissa'
       call get_pop_field(state, i_pop, i, type, weighted_abscissa(i)%ptr)

       ! get source fields (note this is the weighted abscissa source not the abscissa source)
       s_weight(i)%ptr => extract_scalar_field(state, trim(weight(i)%ptr%name)//'Source')
       call get_option(trim(option_path)//'/weighted_abscissa/scalar_field['// &
            int2str(i - 1)//']/name', field_name) 
       s_weighted_abscissa(i)%ptr => extract_scalar_field(state, trim(field_name)//'Source')
       call zero(s_weight(i)%ptr)
       call zero(s_weighted_abscissa(i)%ptr)

       ! relax non-linear values
       ! all temporal relaxations must be the same
       call get_option(trim(weight(i)%ptr%option_path)// &
            '/prognostic/temporal_discretisation/theta', theta)
       ! do not recalculate source terms if theta = 0.0 after first non-linear iteration
       if ((theta == 0.0) .and. (it /= 1)) then 
          return
       end if
       call allocate(r_abscissa(i), abscissa(i)%ptr%mesh, &
            "RelaxedAbscissa"//int2str(i))
       call allocate(r_weight(i), weight(i)%ptr%mesh, "RelaxedWeight"//int2str(i))
       call set(r_abscissa(i), it_abscissa(i)%ptr)
       call scale(r_abscissa(i), theta)
       call addto(r_abscissa(i), abscissa(i)%ptr, 1.0-theta)
       call set(r_weight(i), it_weight(i)%ptr)
       call scale(r_weight(i), theta)
       call addto(r_weight(i), weight(i)%ptr, 1.0-theta)       
    end do    

    call allocate(dummy_scalar, r_abscissa(1)%mesh, name="DummyScalar")

    ! get diffusion
    D => extract_tensor_field(state, trim(weight(1)%ptr%name)//'Diffusivity',stat)
    if (stat == 0) then
       have_D = .true.
    end if

    ! check for growth term
    if (have_option(trim(option_path)//'/population_balance_source_terms/growth')) then
       have_growth = .TRUE.
       if (have_option(trim(option_path)//'/population_balance_source_terms/growth/power_law_growth')) then
          growth_type = 'power_law_growth';
          call get_option(trim(option_path)//'/population_balance_source_terms/growth/power_law_growth', growth_r)
       end if
    else
       have_growth = .FALSE.
    end if

    ! check for internal dispersion term
    if (have_option(trim(option_path)//'/population_balance_source_terms/internal_dispersion')) then
       have_internal_dispersion = .TRUE.
       call get_option(trim(option_path)//'/population_balance_source_terms/internal_dispersion', internal_dispersion_coeff)
    else
       have_internal_dispersion = .FALSE.
    end if

    ! check for aggregation term
    if (have_option(trim(option_path)//'/population_balance_source_terms/aggregation')) then
       have_aggregation = .TRUE.
       if (have_option(trim(option_path)//'/population_balance_source_terms/aggregation/aggregation_frequency/constant_aggregation')) then
          aggregation_freq_type = 'constant_aggregation'
          call get_option(trim(option_path)//'/population_balance_source_terms/aggregation/aggregation_frequency/constant_aggregation', aggregation_freq_const)
       else if (have_option(trim(option_path)//'/population_balance_source_terms/aggregation/aggregation_frequency/hydrodynamic_aggregation')) then
          aggregation_freq_type = 'hydrodynamic_aggregation'
       else if (have_option(trim(option_path)//'/population_balance_source_terms/aggregation/aggregation_frequency/sum_aggregation')) then
          aggregation_freq_type = 'sum_aggregation'
          call get_option(trim(option_path)//'/population_balance_source_terms/aggregation/aggregation_frequency/sum_aggregation', aggregation_freq_const)
       else if (have_option(trim(option_path)//'/population_balance_source_terms/aggregation/aggregation_frequency/laakkonen_2007_aggregation')) then
          aggregation_freq_type = 'laakkonen_2007_aggregation'
          call get_option(trim(option_path)//'/population_balance_source_terms/aggregation/aggregation_frequency/laakkonen_2007_aggregation/C5', C5, default = 0.88)
          if (.not. have_option("/population_balance_continuous_phase_name")) then
             FLAbort("Enable the option population_balance_continuous_phase_name and provide a name for the continuous phase&
                      as it is needed for extracting the turbulence dissipation needed in laakkonen_2007_aggregation kernel")
          end if
          turbulent_dissipation => extract_scalar_field(states(cont_state), "TurbulentDissipation", stat=stat)
          if (stat/=0) then
             FLAbort("I can't find the Turbulent Dissipation field of continuous phase for population balance aggregation term calculations.")
          end if
          if (have_option(trim(states(cont_state)%option_path)//'/subgridscale_parameterisations/k-epsilon')) then
             viscosity_continuous => extract_tensor_field(states(cont_state), "BackgroundViscosity", stat=stat)  
             if (stat/=0) then
                FLAbort("I can't find the Background Viscosity field in k-epsilon for continuous phase for population balance aggregation term calculations.")
             end if
          else 
             viscosity_continuous => extract_tensor_field(states(cont_state), "Viscosity", stat=stat)
             if (stat/=0) then
                FLAbort("I can't find the Viscosity field for continuous phase for population balance aggregation term calculations.")
             end if
          end if 
       end if
    else
       have_aggregation = .FALSE.
    end if

    ! check for breakage term
    if (have_option(trim(option_path)//'/population_balance_source_terms/breakage')) then
       have_breakage = .TRUE.
       if (have_option(trim(option_path)//'/population_balance_source_terms/breakage/breakage_frequency/constant_breakage')) then
          breakage_freq_type = 'constant_breakage'
          call get_option(trim(option_path)//'/population_balance_source_terms/breakage/breakage_frequency/constant_breakage', breakage_freq_const)
       else if (have_option(trim(option_path)//'/population_balance_source_terms/breakage/breakage_frequency/power_law_breakage')) then
          breakage_freq_type = 'power_law_breakage'
          call get_option(trim(option_path)//'/population_balance_source_terms/breakage/breakage_frequency/power_law_breakage/coefficient', breakage_freq_const)
          call get_option(trim(option_path)//'/population_balance_source_terms/breakage/breakage_frequency/power_law_breakage/degree', breakage_freq_degree)
       else if (have_option(trim(option_path)//'/population_balance_source_terms/breakage/breakage_frequency/laakkonen_breakage')) then
          breakage_freq_type = 'laakkonen_frequency'
          if (aggregation_freq_type /= 'laakkonen_2007_aggregation') then
             if (.not. have_option("/population_balance_continuous_phase_name")) then
                FLAbort("Enable the option population_balance_continuous_phase_name and provide a name for the continuous phase &
                         as it is needed for extracting the turbulence dissipation needed in laakkonen_frequency kernel")
             end if
             turbulent_dissipation => extract_scalar_field(states(cont_state), "TurbulentDissipation", stat=stat)
             if (stat/=0) then
                FLAbort("I can't find the Turbulent Dissipation field of continuous phase for population balance breakage term calculations.")
             end if
             if (have_option(trim(states(cont_state)%option_path)//'/subgridscale_parameterisations/k-epsilon')) then
                viscosity_continuous => extract_tensor_field(states(cont_state), "BackgroundViscosity", stat=stat)
                if (stat/=0) then
                   FLAbort("I can't find the Background Viscosity field in k-epsilon for continuous phase for population balance aggregation term calculations.")
                end if
             else 
                viscosity_continuous => extract_tensor_field(states(cont_state), "Viscosity", stat=stat)
                if (stat/=0) then
                   FLAbort("I can't find the Viscosity field for continuous phase for population balance aggregation term calculations.")
                end if
             end if
          end if
       end if

       if (have_option(trim(option_path)//'/population_balance_source_terms/breakage/distribution_function/symmetric_fragmentation')) then
          breakage_dist_type = 'symmetric_fragmentation'
       else if (have_option(trim(option_path)//'/population_balance_source_terms/breakage/distribution_function/mcCoy_madras_2003')) then
          breakage_dist_type = 'mcCoy_madras_2003'
       else if (have_option(trim(option_path)//'/population_balance_source_terms/breakage/distribution_function/laakkonen_2007')) then
          breakage_dist_type = 'laakkonen_2007'
       end if
   
    else
       have_breakage = .FALSE.
    end if

    ! get ill-conditioned matrix settings
    call get_option(trim(option_path)//'/ill_conditioned_matrices/required_condition_number', cond)
    if (have_option(trim(option_path)//'/ill_conditioned_matrices/set_source_to_zero')) then
       singular_option = 'set_source_to_zero';
    else if (have_option(trim(option_path)//'/ill_conditioned_matrices/perturbate')) then 
       singular_option = 'perturbate';
       call get_option(trim(option_path)//'/ill_conditioned_matrices/perturbate/perturbation', perturb_val)
    else if (have_option(trim(option_path)//'/ill_conditioned_matrices/do_nothing')) then
       singular_option = 'do_nothing';
    end if

    X => extract_vector_field(state, 'Coordinate')

    ! assembly loop
    do i = 1, node_count(r_abscissa(1))
       call dqmom_calculate_source_term_node(r_abscissa, r_weight, s_weighted_abscissa, s_weight, &
                &D, have_D, have_growth, growth_type, growth_r, have_internal_dispersion, internal_dispersion_coeff, &
                &have_aggregation, aggregation_freq_type, aggregation_freq_const, C5, &
                &have_breakage, breakage_freq_type, breakage_freq_const, breakage_freq_degree, breakage_dist_type, &
                &turbulent_dissipation, viscosity_continuous, X, singular_option, perturb_val, cond, i)       
    end do
    
    ! Checking if the source terms need to be implemented as absorption
    if(have_option(trim(option_path)//'/apply_source_as_absorption')) then
       do i =1, N
          a_weight(i)%ptr => extract_scalar_field(state, trim(weight(i)%ptr%name)//'Absorption', stat)
          if (stat/=0) then
             FLAbort("Absorption scalar field could not be extracted for population balance weights. How can I apply the source as absorption now!")
          end if
          where (weight(i)%ptr%val >= fields_min)
              a_weight(i)%ptr%val=-1./weight(i)%ptr%val
          elsewhere
              a_weight(i)%ptr%val=-1./fields_min
          end where
          call scale(a_weight(i)%ptr, s_weight(i)%ptr)

          call get_option(trim(option_path)//'/weighted_abscissa/scalar_field['// &
            int2str(i - 1)//']/name', field_name)
          a_weighted_abscissa(i)%ptr => extract_scalar_field(state, trim(field_name)//'Absorption', stat)    
          if (stat/=0) then
             FLAbort("Absorption scalar field could not be extracted for population balance weighted_abscissa. How can I apply the source as absorption now!")
          end if
          ! set dummy_scalar = weight * abscissa
          call set(dummy_scalar, abscissa(i)%ptr)
          call scale(dummy_scalar, weight(i)%ptr)
          where (dummy_scalar%val >= fields_min)
              a_weighted_abscissa(i)%ptr%val=-1./dummy_scalar%val
          elsewhere
              a_weighted_abscissa(i)%ptr%val=-1./fields_min
          end where
          call scale(a_weighted_abscissa(i)%ptr, s_weighted_abscissa(i)%ptr)

          ! make source terms zero now to prevent applying them twice
          call zero(s_weight(i)%ptr)
          call zero(s_weighted_abscissa(i)%ptr)
       end do
       if(have_option(trim(option_path)//'/apply_source_as_absorption/include_sponge_region')) then
          ! extract name of the sponge field
          call get_option(trim(option_path)//'/apply_source_as_absorption/include_sponge_region/sponge_scalar_field_name', field_name)
          sponge_field => extract_scalar_field(state, field_name, stat)
          if (stat/=0) then
             FLAbort("Scalar sponge field could not be located in the state.")
          end if
          ! define temp_field - use dummy_scalar
          call zero(dummy_scalar)
          where (sponge_field%val<0.001)
              dummy_scalar%val = 1.0
          end where
          ! abs_new = abs_old*temp_field + sponge_field... make sure the sponge field stays the same
          do i=1, N
             call scale(a_weight(i)%ptr, dummy_scalar)
             call addto(a_weight(i)%ptr, sponge_field)
             call scale(a_weighted_abscissa(i)%ptr, dummy_scalar)
             call addto(a_weighted_abscissa(i)%ptr, sponge_field)
          end do
       end if

    end if

    ! S = S_c + S_p phi_p. If S is positive, S=S_c otherwise S=S_p phi_p. This makes sure that the scalar remains non-negative. 
    ! See Pg 145 Numerical Heat Transfer and Fluid Flow by Suhas V. Patankar
    if(have_option(trim(option_path)//'/apply_source_as_absorption_for_negative_source_only')) then
       do i =1, N
          a_weight(i)%ptr => extract_scalar_field(state, trim(weight(i)%ptr%name)//'Absorption', stat)
          if (stat/=0) then
             FLAbort("Absorption scalar field could not be extracted for population balance weights. How can I apply the source as absorption now!")
          end if
          call zero(a_weight(i)%ptr)
          do i_node=1, node_count(weight(i)%ptr)
             if (node_val(s_weight(i)%ptr, i_node)<0.0) then
                if (node_val(weight(i)%ptr, i_node)>fields_min) then
                   call set(a_weight(i)%ptr, i_node, -1.0*node_val(s_weight(i)%ptr, i_node)*(1./node_val(weight(i)%ptr, i_node)))
                else
                   call set(a_weight(i)%ptr, i_node, -1.0*node_val(s_weight(i)%ptr, i_node)*(1./fields_min))
                end if
                call set(s_weight(i)%ptr, i_node, 0.0)
             end if
          end do

          call get_option(trim(option_path)//'/weighted_abscissa/scalar_field['// &
            int2str(i - 1)//']/name', field_name)
          a_weighted_abscissa(i)%ptr => extract_scalar_field(state, trim(field_name)//'Absorption', stat)
          if (stat/=0) then
             FLAbort("Absorption scalar field could not be extracted for population balance weighted_abscissa. How can I apply the source as absorption now!")
          end if
          call zero(a_weighted_abscissa(i)%ptr)
          ! set dummy_scalar = weight * abscissa
          call set(dummy_scalar, abscissa(i)%ptr)
          call scale(dummy_scalar, weight(i)%ptr)
          do i_node=1, node_count(weight(i)%ptr)
             if (node_val(s_weighted_abscissa(i)%ptr, i_node)<0.0) then
                if (node_val(dummy_scalar, i_node)>fields_min) then
                   call set(a_weighted_abscissa(i)%ptr, i_node, -1.0*node_val(s_weighted_abscissa(i)%ptr, i_node)*(1./node_val(dummy_scalar, i_node)))
                else
                   call set(a_weighted_abscissa(i)%ptr, i_node, -1.0*node_val(s_weighted_abscissa(i)%ptr, i_node)*(1./fields_min))
                end if 
                call set(s_weighted_abscissa(i)%ptr, i_node, 0.0)
             end if
          end do
       end do

       if(have_option(trim(option_path)//'/apply_source_as_absorption_for_negative_source_only/include_sponge_region')) then
          ! extract name of the sponge field
          call get_option(trim(option_path)//'/apply_source_as_absorption_for_negative_source_only/include_sponge_region/sponge_scalar_field_name', field_name)
          sponge_field => extract_scalar_field(state, field_name, stat)
          if (stat/=0) then
             FLAbort("Scalar sponge field could not be located in the state.")
          end if
          ! define temp_field - use dummy_scalar
          call zero(dummy_scalar)
          where (sponge_field%val<0.001)
              dummy_scalar%val = 1.0
          end where
          ! abs_new = abs_old*temp_field + sponge_field... make sure the sponge field stays the same
          do i=1, N
             call scale(a_weight(i)%ptr, dummy_scalar)
             call addto(a_weight(i)%ptr, sponge_field)
             call scale(a_weighted_abscissa(i)%ptr, dummy_scalar)
             call addto(a_weighted_abscissa(i)%ptr, sponge_field)
          end do
       end if

    end if


    call deallocate(dummy_scalar)
    do i = 1, N
       call deallocate(r_abscissa(i))
       call deallocate(r_weight(i))
    end do
    deallocate(abscissa, weight, it_abscissa, it_weight, &
         r_abscissa, r_weight, weighted_abscissa, s_weight, s_weighted_abscissa, a_weight, a_weighted_abscissa)

  end subroutine dqmom_calculate_source_term_pop

  subroutine dqmom_calculate_source_term_node(abscissa, weight, s_weighted_abscissa, s_weight, &
                 &D, have_D, have_growth, growth_type, growth_r, have_internal_dispersion, internal_dispersion_coeff, &
                 &have_aggregation, aggregation_freq_type, aggregation_freq_const, C5, &
                 &have_breakage, breakage_freq_type, breakage_freq_const, breakage_freq_degree, breakage_dist_type, &
                 &turbulent_dissipation, viscosity_continuous, &
                 &X, singular_option, perturb_val, cond, node)

    type(scalar_field), dimension(:), intent(in) :: abscissa, weight
    type(scalar_field), intent(in) :: turbulent_dissipation
    type(tensor_field), intent(in) :: viscosity_continuous
    type(scalar_field_pointer), dimension(:), intent(inout) :: s_weighted_abscissa, s_weight
    type(tensor_field), pointer, intent(in) :: D
    type(vector_field), pointer, intent(in) :: X
    integer, intent(in) :: node
    real, intent(in) :: cond, growth_r, internal_dispersion_coeff, aggregation_freq_const, breakage_freq_const, breakage_freq_degree, perturb_val, C5
    logical, intent(in) :: have_D, have_growth, have_internal_dispersion, have_aggregation, have_breakage
    character(len=FIELD_NAME_LEN), intent(in) :: growth_type, aggregation_freq_type, breakage_freq_type, breakage_dist_type, singular_option
    
    real, dimension(1, size(abscissa)) :: abscissa_val
    real, dimension(1, size(abscissa)*2, size(abscissa)*2) :: A
    real, dimension(1, size(abscissa)*2) :: S_rhs  ! source term (includes growth, breakage and coalescence term)
    real, dimension(1, size(abscissa)*2, size(abscissa)) :: moment_daughter_dist_func
    real, dimension(1, size(abscissa)) :: break_freq
    real, dimension(1, size(abscissa), size(abscissa)) :: aggregation_freq   ! at present it is not dependent on space coordinate, but can be dependent and will have to be a scalar field
    real, dimension(size(abscissa)*2, 1) :: b
    real, dimension(1, size(abscissa)) :: abscissa_S
    real, dimension(1, size(weight)) :: weight_S
    real :: eps_node
    real, dimension(:,:), allocatable :: visc_node
    real, dimension(size(abscissa)*2, size(abscissa)*2) :: svd_tmp1, svd_tmp2
    real, dimension(size(abscissa)*2) :: SV
    integer :: stat, N, i, j, k

    real :: sigma, density_continuous, density_dispersed

    N = size(abscissa)
    
    ! construct A matrices (lhs knowns)
    do i = 1, N
       abscissa_val(1,i) = node_val(abscissa(i), node)       
    end do
    A = A_matrix(abscissa_val)

    ! initialize dqmom source term to zero 
    S_rhs = 0.0

    ! construct S vector (rhs pt.3) for GROWTH term
    if (have_growth) then
       if (growth_type=='power_law_growth') then
          do i = 1, 2*N
             do j = 1, N
                S_rhs(1,i) = S_rhs(1,i) + (i-1)*node_val(weight(j), node)*abscissa_val(1,j)**(i-2+growth_r)
             end do
          end do
       end if
    end if

    ! construct S vector (rhs pt.3) for INTERNAL DISPERSION
    if (have_internal_dispersion) then
       do i = 1, 2*N
          do j = 1, N
             S_rhs(1,i) = S_rhs(1,i) + (i-1)*(i-2)*node_val(weight(j), node)*(abscissa_val(1,j)**(i-3))*internal_dispersion_coeff
          end do
       end do
    end if

    ! construct S vector for BREAKAGE
    
    if (have_breakage) then
       if (breakage_freq_type=='constant_breakage') then
          break_freq = breakage_freq_const
       else if (breakage_freq_type=='power_law_breakage') then
          do i = 1, N
             break_freq(1,i) = breakage_freq_const*abscissa_val(1,i)**breakage_freq_degree
          end do
       else if (breakage_freq_type=='laakkonen_frequency') then
          density_continuous = 998.2
          density_dispersed = 1.205
          sigma = 0.072
          eps_node = node_val(turbulent_dissipation,node)
          ! Assuming isotropic molecular viscosity here   
          allocate(visc_node(viscosity_continuous%dim(1), viscosity_continuous%dim(1)))
          visc_node = node_val(viscosity_continuous,node)
          do i = 1, N
             break_freq(1,i) = 6.0*eps_node**(1./3) * erfc(sqrt( 0.04*(sigma/density_continuous)*(1./(eps_node**(2./3) * abscissa_val(1,i)**(5./3))) + 0.01*(visc_node(1,1)/sqrt(density_continuous*density_dispersed))*(1./(eps_node**(1./3)*abscissa_val(1,i)**(4./3)))))
          end do
          deallocate(visc_node)
       end if

       if (breakage_dist_type=='symmetric_fragmentation') then
          do i = 1, 2*N
             do j = 1, N
                moment_daughter_dist_func(1,i,j) = (2.0**(((3-(i-1))/3.0)))*(abscissa_val(1,j)**(i-1))   
             end do
          end do
       else if (breakage_dist_type=='mcCoy_madras_2003') then
          do i = 1, 2*N
             do j = 1, N
                moment_daughter_dist_func(1,i,j) = (6.0/((i-1)+3.0))*(abscissa_val(1,j)**(i-1))
             end do
          end do
       else if (breakage_dist_type=='laakkonen_2007') then
          do i = 1, 2*N
             do j = 1, N
                moment_daughter_dist_func(1,i,j) = 180.0*(abscissa_val(1,j)**(i-1))*(1./((i-1)+15) - 2./((i-1)+12) + 1./((i-1)+9))
             end do
          end do
       end if

       do i = 1, 2*N 
          do j = 1, N
             ! birth term due to breakage
             S_rhs(1,i) = S_rhs(1,i) + break_freq(1,j)*node_val(weight(j), node)*moment_daughter_dist_func(1,i,j)   ! daughter distribution function already includes the factor for number of particles formed after breakage
             ! death term due to breakage
             S_rhs(1,i) = S_rhs(1,i) - break_freq(1,j)*node_val(weight(j), node)*(abscissa_val(1,j)**(i-1))
          end do
       end do
    endif

    
    !!! construct S vector for AGGREGATION
    if (have_aggregation) then
       if (aggregation_freq_type=='constant_aggregation') then
          aggregation_freq = aggregation_freq_const
       else if (aggregation_freq_type=='hydrodynamic_aggregation') then
          do i = 1, N
             do j = 1, N
                aggregation_freq(1,i,j) = abscissa_val(1,i)**3 + abscissa_val(1,j)**3
             end do
          end do
       else if (aggregation_freq_type=='sum_aggregation') then
          do i = 1, N
             do j = 1, N
                aggregation_freq(1,i,j) = (abscissa_val(1,i) + abscissa_val(1,j))*aggregation_freq_const
             end do
          end do
       else if (aggregation_freq_type=='laakkonen_2007_aggregation') then
          density_continuous = 998.2
          sigma = 0.072
          eps_node = node_val(turbulent_dissipation, node)
          ! Assuming isotropic molecular viscosity here   
          allocate(visc_node(viscosity_continuous%dim(1), viscosity_continuous%dim(1)))
          visc_node=node_val(viscosity_continuous, node)
          do i = 1, N
             do j = 1, N
                aggregation_freq(1,i,j) = C5 * eps_node**(1./3) * (abscissa_val(1,i) + abscissa_val(1,j))**2 * (abscissa_val(1,i)**(2./3) + abscissa_val(1,j)**(2./3))**(1./2) * exp(-6.0E9*((visc_node(1,1)*density_continuous)/sigma**2)*eps_node*((abscissa_val(1,i)*abscissa_val(1,j))/(abscissa_val(1,i)+abscissa_val(1,j)))**4)
             end do
          end do
          deallocate(visc_node)
       end if

       do i = 1, 2*N
          do j = 1, N
             do k = 1, N
                ! birth term due to aggregation
                S_rhs(1,i) = S_rhs(1,i) + 0.5 * aggregation_freq(1,j,k) * node_val(weight(j), node) * node_val(weight(k), node) * &
                                           ((abs(abscissa_val(1,j)**3 + abscissa_val(1,k)**3)**(1.0/3.0)) * &
                                           sign(1.0,(abscissa_val(1,j)**3 + abscissa_val(1,k)**3)))**(i-1)
                ! death term due to aggregation
                S_rhs(1,i) = S_rhs(1,i) - aggregation_freq(1,j,k) * node_val(weight(j), node) * node_val(weight(k), node) * abscissa_val(1,j)**(i-1)
             end do
          end do
       end do
    endif

    ! check for ill-conditioned matrices
    if (singular_option=='set_source_to_zero') then
       call svd(A(1,:,:), svd_tmp1, SV, svd_tmp2)
       if (SV(size(SV))/SV(1) < cond) then
          ewrite(2,*) 'ill-conditioned matrix found', SV(size(SV))/SV(1)
          S_rhs(1,:)=0.0
          A(1,:,:) = 0.0
          do j = 1, 2*N
             A(1,j,j) = 1.0
          end do
       end if
    else if (singular_option=='do_nothing') then
       call svd(A(1,:,:), svd_tmp1, SV, svd_tmp2)
       if (SV(size(SV))/SV(1) < cond) then
          ewrite(2,*) 'ill-conditioned matrix found but doing nothing about it', SV(size(SV))/SV(1)
       end if
    else if(singular_option=='perturbate') then
       call svd(A(1,:,:), svd_tmp1, SV, svd_tmp2)
       if (SV(size(SV))/SV(1) < cond) then
          ewrite(2,*) 'ill-conditioned matrix found and perturbating', SV(size(SV))/SV(1)
          do i=1,N-1
             abscissa_val(1,i)=abscissa_val(1,i)-perturb_val
          end do
          A=A_matrix(abscissa_val)
          call svd(A(1,:,:), svd_tmp1, SV, svd_tmp2)
          ewrite(2,*) 'Condition number after perturbating', SV(size(SV))/SV(1)
       end if
    end if

    ! solve linear system to find source values for weights and weighted-abscissa equations
    b(:,1) = S_rhs(1,:)   
    call dqmom_solve(A(1,:,:), b, stat)
    weight_S(1,:) = b(:N,1)
    abscissa_S(1,:) = b(N+1:,1)

    ! Add to source fields
    do i = 1, N
       call addto(s_weight(i)%ptr, node, weight_S(1,i))
       call addto(s_weighted_abscissa(i)%ptr, node, abscissa_S(1,i))
    end do

  end subroutine dqmom_calculate_source_term_node

  function A_matrix(abscissa)

    real, dimension(:,:), intent(in) :: abscissa

    real, dimension(size(abscissa,1), size(abscissa,2)*2, size(abscissa,2)*2) :: A_matrix
    integer :: i, j, N

    N = size(abscissa,2)        
    do i = 1, 2*N
       do j = 1, N
          A_matrix(:,i,j) = (2-i)*abscissa(:,j)**(i-1)
          A_matrix(:,i,j+N) = (i-1)*abscissa(:,j)**(i-2)
       end do
    end do

  end function A_matrix

  subroutine dqmom_solve(A, b, stat)
    !!< Solve Ax=b for right hand sides b putting the result in b.
    !!<
    !!< This is simply a wrapper for lapack.
    real, dimension(:,:), intent(in) :: A
    real, dimension(:,:), intent(inout) :: b
    integer, optional, intent(out) :: stat

    real, dimension(size(A,1), size(A,2)) :: Atmp
    integer, dimension(size(A,1)) :: ipiv
    integer :: info

    interface 
#ifdef DOUBLEP
       SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
         INTEGER :: INFO, LDA, LDB, N, NRHS
         INTEGER :: IPIV( * )
         REAL ::  A( LDA, * ), B( LDB, * )
       END SUBROUTINE DGESV
#else
       SUBROUTINE SGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
         INTEGER :: INFO, LDA, LDB, N, NRHS
         INTEGER :: IPIV( * )
         REAL ::  A( LDA, * ), B( LDB, * )
       END SUBROUTINE SGESV
#endif
    end interface

    if (present(stat)) stat = 0

    Atmp=A

#ifdef DOUBLEP
    call dgesv(&
#else
    call sgesv(&
#endif
    size(A,1), size(b,2), Atmp, size(A,1), ipiv, b, size(b,1), info) 
    
    if (present(stat)) then
       stat = info
    end if

  end subroutine dqmom_solve

  subroutine dqmom_calculate_moments(state)

    type(state_type), intent(in) :: state

    type(scalar_field), pointer :: abscissa, weight, moment
    type(scalar_field) :: work
    integer :: i_pop, N, i_moment, i_abscissa, i
    character(len=OPTION_PATH_LEN) :: option_path 
    character(len=FIELD_NAME_LEN) :: type
    
    ewrite(1, *) "In dqmom_calculate_moments"
    do i_pop = 1, option_count(trim(state%option_path)//'/population_balance')
       call get_pop_option_path(state, i_pop, option_path)
       N = option_count(trim(option_path)//'/abscissa/scalar_field')
       do i_moment = 1, 2*N
          type = 'moments'
          call get_pop_field(state, i_pop, i_moment, type, moment)
          call zero(moment)
          call allocate(work, moment%mesh, 'work')
          do i_abscissa = 1, N
             ! get required fields
             type = 'abscissa'
             call get_pop_field(state, i_pop, i_abscissa, type, abscissa)
             type = 'weights'
             call get_pop_field(state, i_pop, i_abscissa, type, weight)

             ! calculate moment -> m(i) = sum_j(w(j)*q(j)**i)
             call zero(work)
             call set(work, 1.0)
             do i = 1, i_moment - 1
                call scale(work, abscissa)
             end do
             call scale(work, weight)
             call addto(moment, work)
          end do
          call deallocate(work)
       end do
    end do
    ewrite(1, *) "Exiting dqmom_calculate_moments"
  end subroutine dqmom_calculate_moments

  subroutine dqmom_calculate_statistics(state)

    type(state_type), intent(in) :: state

    type(scalar_field_pointer), dimension(:), allocatable :: moments
    type(scalar_field), pointer :: stats
    integer :: i_pop, N, i, j, i_stat
    real :: mean, std, scaling_factor_Dia
    character(len=OPTION_PATH_LEN) :: option_path 
    character(len=FIELD_NAME_LEN) :: type
    
    ewrite(1, *) "In dqmom_calculate_statistics"
    do i_pop = 1, option_count(trim(state%option_path)//'/population_balance')
       call get_pop_option_path(state, i_pop, option_path)

       if (option_count(trim(option_path)//'/statistics/scalar_field') > 0) then
          
          N = option_count(trim(option_path)//'/moments/scalar_field')
          allocate(moments(N))
          do i = 1, N
             ! get required fields
             type = 'moments'
             call get_pop_field(state, i_pop, i, type, moments(i)%ptr)
          end do

          type = 'statistics'
          do i_stat = 1, option_count(trim(option_path)//'/statistics/scalar_field')
             call get_pop_field(state, i_pop, i_stat, type, stats)
             if (trim(stats%name) == "Mean") then
                call zero(stats)
                do j = 1, node_count(stats)
                   call set(stats, j, node_val(moments(2)%ptr,j)/node_val(moments(1)%ptr,j))
                end do
             end if
             if (trim(stats%name) == "StandardDeviation") then
                call zero(stats)
                do j = 1, node_count(stats)
                   mean = node_val(moments(2)%ptr,j)/node_val(moments(1)%ptr,j)
                   call set(stats, j, (node_val(moments(3)%ptr,j)/node_val(moments(1)%ptr,j) -&
                        mean**2)**0.5)
                end do
             end if
             if (trim(stats%name) == "Skew") then
                call zero(stats)
                do j = 1, node_count(stats)
                   mean = node_val(moments(2)%ptr,j)/node_val(moments(1)%ptr,j)
                   std = (node_val(moments(3)%ptr,j)/node_val(moments(1)%ptr,j) - mean**2)**0.5
                   call set(stats, j, &
                        (node_val(moments(4)%ptr,j)/node_val(moments(1)%ptr,j) - &
                        3*mean*std**2.0 - mean**3.0) / std**3.0)
                end do
             end if
             if (trim(stats%name) == "SauterMeanDia") then
                call zero(stats)
                if (have_option(trim(option_path)//'/scaling_factor_Dia')) then
                   call get_option(trim(option_path)//'/scaling_factor_Dia', scaling_factor_Dia)
                else
                   scaling_factor_Dia=1.0 
                end if
                do j = 1, node_count(stats)
                   call set(stats, j, scaling_factor_Dia*(node_val(moments(4)%ptr,j)/node_val(moments(3)%ptr,j)) )
                end do
             end if
             if (trim(stats%name) == "MeanDia10") then
                call zero(stats)
                if (have_option(trim(option_path)//'/scaling_factor_Dia')) then
                   call get_option(trim(option_path)//'/scaling_factor_Dia', scaling_factor_Dia)
                else
                   scaling_factor_Dia=1.0
                end if
                do j = 1, node_count(stats)
                   call set(stats, j, scaling_factor_Dia*(node_val(moments(2)%ptr,j)/node_val(moments(1)%ptr,j)) )
                end do
             end if

          end do

          deallocate(moments)

       end if
    end do
    ewrite(1, *) "Exiting dqmom_calculate_statistics"


  end subroutine dqmom_calculate_statistics

  subroutine dqmom_check_options

!    type(state_type), intent(in) :: state
    
    integer :: i_pop, i_state, i_field, stat
    integer :: n_abscissa, n_weights, n_weighted_abscissa, n_moments, n_statistics
    character(len=OPTION_PATH_LEN) :: option_path 
    character(len=FIELD_NAME_LEN)  :: old_msh, amsh, wmsh, wamsh, mmsh, smsh, type

    write(*,*) 'in dqmom_check_options'

    do i_state = 1, option_count("/material_phase")
       do i_pop = 1, option_count('material_phase['//int2str(i_state-1)//']/population_balance')
          option_path = 'material_phase['//int2str(i_state-1)//']/population_balance['//int2str(i_pop-1)//']'

          ! Check there are the same number of abscissa, weight and weighted abscissa fields
          n_abscissa = option_count(trim(option_path)//'/abscissa/scalar_field')
          n_weights = option_count(trim(option_path)//'/weights/scalar_field')
          n_weighted_abscissa = option_count(trim(option_path)// &
               '/weighted_abscissa/scalar_field')
          if((n_weights /= n_abscissa) .or. (n_weights /= n_weighted_abscissa)) then
             FLExit("The number of weights, abscissas and weighted abscissa scalar fields must be the same in the population balance solver")
          end if

          ! Check there are sufficient abscissas to calculate requested moments
          n_moments = option_count(trim(option_path)//'/moments/scalar_field')
          if((n_moments / 2) /= n_abscissa) then
             FLExit("The number of moments must be twice the number of abscissas in the population balance solver")
          end if       

          ! Need to check that all fields are on the same mesh
    
          type='abscissa'
          do i_field = 1, n_abscissa
             call get_option(trim(option_path)//'/'//trim(type)//'/scalar_field['//int2str(i_field - 1)//']/diagnostic/mesh/name', amsh, stat)
             if (stat /= 0) then
                FLExit("Abscissa scalar field must be diagnostic in population balance solver")
             else
                if (i_field==1) then 
                   old_msh=amsh
                else 
                   if (trim(amsh)/=trim(old_msh)) then
                      FLExit("All abscissas must be on the same mesh")
                   else
                      old_msh=amsh
                   end if
                end if
             end if  
          end do

          type='weights'
          do i_field = 1, n_abscissa
             call get_option(trim(option_path)//'/'//trim(type)//'/scalar_field['//int2str(i_field - 1)//']/prognostic/mesh/name', wmsh, stat)
             if (stat /= 0) then
                FLExit("Weight scalar field must be prognostic in population balance solver")
             else
                if (i_field==1) then
                   old_msh=wmsh
                else                
                   if (trim(wmsh)/=trim(old_msh)) then
                      FLExit("All weights must be on the same mesh")
                   else
                      old_msh=wmsh
                   end if
                end if
             end if             
          end do

          type='weighted_abscissa'
          do i_field = 1, n_abscissa
             call get_option(trim(option_path)//'/'//trim(type)//'/scalar_field['//int2str(i_field - 1)//']/prognostic/mesh/name', wamsh, stat)
             if (stat /= 0) then
                FLExit("Weighted abscissa scalar field must be prognostic in population balance solver")
             else
                if (i_field==1) then
                   old_msh=wamsh
                else
                   if (trim(wamsh)/=trim(old_msh)) then
                      FLExit("All weighted abscissas must be on the same mesh")
                   else
                      old_msh=wamsh
                   end if
                end if
             end if
          end do

          type='moments'
          do i_field = 1, n_moments
             call get_option(trim(option_path)//'/'//trim(type)//'/scalar_field['//int2str(i_field - 1)//']/diagnostic/mesh/name', mmsh, stat)
             if (stat /= 0) then
                FLExit("Moment scalar field must be diagnostic in population balance solver")
             else
                if (i_field==1) then
                   old_msh=mmsh
                else
                   if (trim(mmsh)/=trim(old_msh)) then
                      FLExit("All moments must be on the same mesh")
                   else
                      old_msh=mmsh
                   end if
                end if
             end if
          end do          

          n_statistics = option_count(trim(option_path)//'/statistics/scalar_field')
          type='statistics'
          do i_field = 1, n_statistics
             call get_option(trim(option_path)//'/'//trim(type)//'/scalar_field['//int2str(i_field - 1)//']/diagnostic/mesh/name', smsh, stat)
             if (stat /= 0) then
                FLExit("Statitics scalar field must be diagnostic in population balance solver")
             else
                if (i_field==1) then
                   old_msh=smsh
                else
                   if (trim(smsh)/=trim(old_msh)) then
                      FLExit("All statistics must be on the same mesh")
                   else
                      old_msh=smsh
                   end if
                end if
             end if
          end do        

          if (n_statistics == 0) then
             smsh = amsh
          end if

          if ((amsh/=wmsh) .or. (amsh/=wamsh) .or. (amsh/=mmsh) .or. (amsh/=smsh)) then
             FLExit("Abscissas, weights, weighted abscissas, moments and statistics - all must be on the same mesh in population balance solver")
          end if

       end do
    end do

    write(*,*) 'Finished dqmom_check_options'
  end subroutine dqmom_check_options

end module dqmom
