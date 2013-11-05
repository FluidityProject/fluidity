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

  use vector_tools
  use spud
  use fields
  use state_module
  use state_fields_module
  use initialise_fields_module
  use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
!  use sparsity_patterns_meshes
  use sparse_tools
  use solvers

  implicit none

  public dqmom_init, dqmom_calculate_source_terms, dqmom_calculate_abscissa,&
       & dqmom_check_options, dqmom_calculate_moments, dqmom_calculate_statistics

  private

  !! TODO:
  !! 1. Make field_finding algorithm work for multiple phases and make sure it works for
  !!    several pop_balances
  !! 3. Make check_options run
  !! 4. Check all prognostic fields are identical (possible bar initial conditions)
  !! 5. Check prognostic Source terms are set to diagnostic, Internal
  !! 6. Minimum weight
  !! 7. Check diagnostic fields are set to Internal
  !! 8. Entire set up of fields is a horrible hack - could do with sorting out - suggest 
  !!    new blueprint for how this could work - will only work for one pop_balance per phase. IMPORTANT
  !! 9. Make changes to function dqmom_calculate_abscissa so that it works for a state insteaad of all states. 
  !!    Make necessary changes to the place this is called in the code. 
  !! 10.In Fluids.F90 use conditions on presence of population balance tree and call the functions inside the if condition.
  !!    Loop on states inside Fluids.F90 (or any other file where DQMOM functions are called) instead of doing it inside DQMOM functions.  

  !! Changes made by gb812
  !! 1. in 'Construct C matrix (rhs pt. 2)' -> took grad_D out or squared term
  !! 2. in 'Construct C matrix (rhs pt. 2)' -> while summing grad_D in the next line, changed from sum(grad_D) to sum(grad_D,1) 
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
            &//int2str(i_pop)//']'
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

  subroutine dqmom_calculate_source_terms(states, it)

    type(state_type), dimension(:), intent(inout) :: states
    integer, intent(in) :: it    

    integer :: i_state, i_pop

    do i_state = 1, option_count("/material_phase")
       do i_pop = 1, option_count(trim(states(i_state)%option_path)//'/population_balance')
          call dqmom_calculate_source_term_pop(states(i_state), it, i_pop)
       end do
    end do

  end subroutine dqmom_calculate_source_terms

  subroutine dqmom_calculate_source_term_pop(state, it, i_pop)

    type(state_type), intent(inout) :: state
    integer, intent(in) :: it

    type(scalar_field_pointer), dimension(:), allocatable :: abscissa,&
         weight, it_abscissa, it_weight, s_weighted_abscissa, s_weight
    type(scalar_field), pointer :: lumped_mass
    type(csr_matrix), pointer :: mass_matrix
    type(scalar_field), dimension(:), allocatable :: r_abscissa, r_weight
    type(tensor_field), pointer :: D
    type(vector_field), pointer :: X
    real :: theta, cond, growth_r, internal_dispersion_coeff, aggregation_freq_const, breakage_freq_const, breakage_freq_degree, perturb_val
    integer :: i_pop, N, i, j, stat
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
         r_abscissa(N), r_weight(N), s_weighted_abscissa(N), s_weight(N))

    do i = 1, N
       ! collect abscissa and weight fields
       type = 'weights'
       call get_pop_field(state, i_pop, i, type, weight(i)%ptr)
       call get_pop_field(state, i_pop, i, type, it_weight(i)%ptr, iterated=.true.)
       type = 'abscissa'
       call get_pop_field(state, i_pop, i, type, abscissa(i)%ptr)
       call get_pop_field(state, i_pop, i, type, it_abscissa(i)%ptr, iterated=.true.)

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
       end if

       if (have_option(trim(option_path)//'/population_balance_source_terms/breakage/distribution_function/symmetric_fragmentation')) then
          breakage_dist_type = 'symmetric_fragmentation'
       else if (have_option(trim(option_path)//'/population_balance_source_terms/breakage/distribution_function/mcCoy_madras_2003')) then
          breakage_dist_type = 'mcCoy_madras_2003'
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
!    else if (have_option(trim(option_path)//'/ill_conditioned_matrices/average_surrounding_nodes')) then        
!       singular_option = 'average_surrounding_nodes';
    end if

    X => extract_vector_field(state, 'Coordinate')
    ! assembly loop
    do i = 1, ele_count(r_abscissa(1))
       call dqmom_calculate_source_term_ele(r_abscissa, r_weight, s_weighted_abscissa, s_weight, &
                &D, have_D, have_growth, growth_type, growth_r, have_internal_dispersion, internal_dispersion_coeff, &
                &have_aggregation, aggregation_freq_type, aggregation_freq_const, &
                &have_breakage, breakage_freq_type, breakage_freq_const, breakage_freq_degree, breakage_dist_type, X, singular_option, perturb_val, cond, i)       
    end do

    ! for non-DG we apply inverse mass globally
    if(continuity(r_abscissa(1))>=0) then
       if(have_option(trim(option_path)//'/adv_diff_source_term_interpolation/use_full_mass_matrix')) then
          mass_matrix => get_mass_matrix(state, r_abscissa(1)%mesh)
          ! r_abscissa(1) is used as a dummy scalar field here since it is not needed anymore. 
          do j = 1, N
             call zero(r_abscissa(1))
             call petsc_solve(r_abscissa(1), mass_matrix, s_weight(j)%ptr, trim(option_path)//'/adv_diff_source_term_interpolation/use_full_mass_matrix')
             call set(s_weight(j)%ptr, r_abscissa(1))
             call zero(r_abscissa(1))
             call petsc_solve(r_abscissa(1), mass_matrix, s_weighted_abscissa(j)%ptr, trim(option_path)//'/adv_diff_source_term_interpolation/use_full_mass_matrix')
             call set(s_weighted_abscissa(j)%ptr, r_abscissa(1))
          end do
       else if(have_option(trim(option_path)//'/adv_diff_source_term_interpolation/use_mass_lumping')) then
          lumped_mass => get_lumped_mass(state, r_abscissa(1)%mesh)
          do j = 1, N
             do i = 1, node_count(r_abscissa(1))
                call set(s_weighted_abscissa(j)%ptr, i, node_val(s_weighted_abscissa(j)%ptr,i)&
                     &/node_val(lumped_mass,i))
                call set(s_weight(j)%ptr, i, node_val(s_weight(j)%ptr,i)&
                     &/node_val(lumped_mass,i))
             end do
          end do
       else 
          FLAbort("Check the .flml file. You must specify an option under 'population_balance/adv_diff_source_term_interpolation'")
       end if
    end if

    do i = 1, N
       call deallocate(r_abscissa(i))
       call deallocate(r_weight(i))
    end do
    deallocate(abscissa, weight, it_abscissa, it_weight, &
         r_abscissa, r_weight)

  end subroutine dqmom_calculate_source_term_pop

  subroutine dqmom_calculate_source_term_ele(abscissa, weight, s_weighted_abscissa, s_weight, &
                 &D, have_D, have_growth, growth_type, growth_r, have_internal_dispersion, internal_dispersion_coeff, &
                 &have_aggregation, aggregation_freq_type, aggregation_freq_const, &
                 &have_breakage, breakage_freq_type, breakage_freq_const, breakage_freq_degree, breakage_dist_type, X, singular_option, perturb_val, cond, ele)

    type(scalar_field), dimension(:), intent(in) :: abscissa, weight
    type(scalar_field_pointer), dimension(:), intent(inout) :: s_weighted_abscissa, s_weight
    type(tensor_field), pointer, intent(in) :: D
    type(vector_field), pointer, intent(in) :: X
    integer, intent(in) :: ele
    real, intent(in) :: cond, growth_r, internal_dispersion_coeff, aggregation_freq_const, breakage_freq_const, breakage_freq_degree, perturb_val
    logical, intent(in) :: have_D, have_growth, have_internal_dispersion, have_aggregation, have_breakage
    character(len=FIELD_NAME_LEN), intent(in) :: growth_type, aggregation_freq_type, breakage_freq_type, breakage_dist_type, singular_option

    real, dimension(ele_ngi(abscissa(1), ele), size(abscissa)) :: abscissa_val_at_quad
    real, dimension(ele_ngi(abscissa(1), ele), size(abscissa)*2, size(abscissa)*2) :: A
    real, dimension(ele_ngi(abscissa(1), ele), size(abscissa)*2, size(abscissa)) :: A_3
    real, dimension(ele_ngi(abscissa(1), ele), size(abscissa)) :: C
    real, dimension(ele_ngi(abscissa(1), ele), size(abscissa)*2) :: S_rhs  ! source term (includes growth, breakage and coalescence term): gb 15-11-2012
    real, dimension(ele_ngi(abscissa(1), ele), size(abscissa)*2, size(abscissa)) :: moment_daughter_dist_func
    real, dimension(ele_ngi(abscissa(1), ele), size(abscissa)) :: break_freq
    real, dimension(ele_ngi(abscissa(1), ele), size(abscissa), size(abscissa)) :: aggregation_freq   ! at present it is not dependent on space coordinate, but can be dependent and will have to be a scalar fields
    real, dimension(size(abscissa)*2, 1) :: b
    real, dimension(ele_ngi(abscissa(1), ele), size(abscissa)) :: abscissa_S_at_quad
    real, dimension(ele_ngi(abscissa(1), ele), size(weight)) :: weight_S_at_quad
    real, dimension(ele_loc(abscissa(1), ele)) :: abscissa_S_at_nodes
    real, dimension(ele_loc(abscissa(1), ele)) :: weight_S_at_nodes
    real, dimension(ele_loc(abscissa(1), ele), ele_loc(abscissa(1), ele)) :: invmass
    real, dimension(ele_ngi(abscissa(1), ele)) :: detwei
    real, dimension(X%dim, ele_ngi(abscissa(1), ele)) :: grad_D
    real, dimension(X%dim, X%dim, ele_ngi(abscissa(1), ele)) :: D_at_quad
    type(element_type), pointer :: shape
    integer, dimension(:), pointer :: nodes
    real, dimension(ele_loc(abscissa(1), ele), ele_ngi(abscissa(1), ele), X%dim) :: dshape
    real, dimension(size(abscissa)*2, size(abscissa)*2) :: svd_tmp1, svd_tmp2
    real, dimension(size(abscissa)*2) :: SV
    integer :: stat, N, i, j, k
    real :: curr_time
    real, dimension(size(abscissa)) :: perturb
    integer :: p_num
    real, dimension(size(abscissa)) :: absc_init

    N = size(abscissa)
    
    nodes => ele_nodes(abscissa(1), ele)
    shape => ele_shape(abscissa(1), ele)

    call transform_to_physical(X, ele, shape, dshape=dshape, detwei=detwei)

    ! construct A matrices (lhs knowns)
    do i = 1, N
       abscissa_val_at_quad(:,i) = ele_val_at_quad(abscissa(i), ele)       
!!       print*, "abscissa = ", abscissa_val_at_quad(1,i)
    end do
    A = A_matrix(abscissa_val_at_quad)

    ! construct A_3 matrix (rhs pt.1)
    do i = 1, 2*N
       do j = 1, N
          A_3(:,i,j) = (i-1)*(i-2)*ele_val_at_quad(abscissa(j), ele)**(i-3)
       end do
    end do

    ! construct C matrix (rhs pt.2)
    if (have_D) then
       do i = 1, N
          D_at_quad = ele_val_at_quad(D, ele)
          do j = 1, X%dim
             grad_D(j,:) = D_at_quad(j,j,:)
          end do
          grad_D = ((ele_grad_at_quad(abscissa(i), ele, dshape))**2)*grad_D
          C(:,i) = ele_val_at_quad(weight(i), ele)*sum(grad_D,1)
       end do
    else
       C = 0.0
    end if


    !! gb nov-2012---
    ! initialize dqmom source term to zero 
    S_rhs = 0.0

    ! construct S vector (rhs pt.3) for GROWTH term
    if (have_growth) then
       if (growth_type=='power_law_growth') then
          do i = 1, 2*N
             do j = 1, N
                S_rhs(:,i) = S_rhs(:,i) + (i-1)*ele_val_at_quad(weight(j), ele)*(ele_val_at_quad(abscissa(j), ele)**(i-2))*(ele_val_at_quad(abscissa(j), ele)**growth_r)
             end do 
          end do 
       end if
    end if
    
    if (have_internal_dispersion) then
       do i = 1, 2*N
          do j = 1, N
!             S_rhs(:,i) = S_rhs(:,i) + (i-1)*(i-2)*ele_val_at_quad(weight(j), ele)*(ele_val_at_quad(abscissa(j), ele)**(i-3))*Diffusion_internal(j)
             S_rhs(:,i) = S_rhs(:,i) + (i-1)*(i-2)*ele_val_at_quad(weight(j), ele)*(abscissa_val_at_quad(:,j)**(i-3))*internal_dispersion_coeff   ! abscissa_val_at_quad() takes care of the perturbed abscissas
          end do
       end do
    end if

    !!! construct S vector for BREAKAGE
!    do i = 1, N
!       break_freq(:,i) = 0.02*abscissa_val_at_quad(:,i)**3
!    end do
    
    if (have_breakage) then
       if (breakage_freq_type=='constant_breakage') then
          break_freq = breakage_freq_const
       else if (breakage_freq_type=='power_law_breakage') then
          do i = 1, N
             break_freq(:,i) = breakage_freq_const*abscissa_val_at_quad(:,i)**breakage_freq_degree
          end do
       end if

       if (breakage_dist_type=='symmetric_fragmentation') then
          do i = 1, 2*N
             do j = 1, N
                moment_daughter_dist_func(:,i,j) = (2.0**(((3-(i-1))/3.0)))*(abscissa_val_at_quad(:,j)**(i-1))   
!                moment_daughter_dist_func(:,i,j) = (2.0**(((3-(i-1))/3.0)))*(ele_val_at_quad(abscissa(j), ele)**(i-1))
             end do
          end do
       else if (breakage_dist_type=='mcCoy_madras_2003') then
          do i = 1, 2*N
             do j = 1, N
                moment_daughter_dist_func(:,i,j) = (6.0/((i-1)+3.0))*(abscissa_val_at_quad(:,j)**(i-1))
!                moment_daughter_dist_func(:,i,j) = (6.0/((i-1)+3.0))*(ele_val_at_quad(abscissa(j), ele)**(i-1))
             end do
          end do
       end if

       do i = 1, 2*N 
          do j = 1, N
             ! birth term due to breakage
             S_rhs(:,i) = S_rhs(:,i) + break_freq(:,j)*ele_val_at_quad(weight(j), ele)*moment_daughter_dist_func(:,i,j)   ! daughter distribution function already includes the factor for number of particles formed after breakage
             ! death term due to breakage
             S_rhs(:,i) = S_rhs(:,i) - break_freq(:,j)*ele_val_at_quad(weight(j), ele)*(abscissa_val_at_quad(:,j)**(i-1))
!             S_rhs(:,i) = S_rhs(:,i) - break_freq(:,j)*ele_val_at_quad(weight(j), ele)*(ele_val_at_quad(abscissa(j), ele)**(i-1))
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
                aggregation_freq(:,i,j) = abscissa_val_at_quad(:,i)**3 + abscissa_val_at_quad(:,j)**3
             end do
          end do
       end if

       do i = 1, 2*N
          do j = 1, N
             do k = 1, N
                ! birth term due to aggregation
                S_rhs(:,i) = S_rhs(:,i) + 0.5 * aggregation_freq(:,j,k) * ele_val_at_quad(weight(j), ele) * ele_val_at_quad(weight(k), ele) * &
                                           ((abs(abscissa_val_at_quad(:,j)**3 + abscissa_val_at_quad(:,k)**3)**(1.0/3.0)) * &
                                           sign(1.0,(abscissa_val_at_quad(:,j)**3 + abscissa_val_at_quad(:,k)**3)))**(i-1)
!                S_rhs(:,i) = S_rhs(:,i) + 0.5 * aggregation_freq(:,j,k) * ele_val_at_quad(weight(j), ele) * ele_val_at_quad(weight(k), ele) * &
!                                           &(ele_val_at_quad(abscissa(j), ele)**3 + ele_val_at_quad(abscissa(k), ele)**3)**((i-1)/3.0)
                ! death term due to aggregation
                S_rhs(:,i) = S_rhs(:,i) - aggregation_freq(:,j,k) * ele_val_at_quad(weight(j), ele) * ele_val_at_quad(weight(k), ele) * abscissa_val_at_quad(:,j)**(i-1)
!                S_rhs(:,i) = S_rhs(:,i) - aggregation_freq(:,j,k) * ele_val_at_quad(weight(j), ele)*ele_val_at_quad(weight(k), ele)*ele_val_at_quad(abscissa(j), ele)**(i-1)
             end do
          end do
       end do
    endif

    !! ----------------

!    perturb(1) = 0.00
!    perturb(2) = 0.01
!    perturb(3) = -0.01
    ! check for ill-conditioned matrices
    if (singular_option=='perturbate') then
       do i = 1, ele_ngi(abscissa(1), ele)
          call svd(A(i,:,:), svd_tmp1, SV, svd_tmp2)
          p_num=1
          do while (SV(size(SV))/SV(1) < cond)
          !if (SV(size(SV))/SV(1) < cond) then
             print*, "SINGULAR MATRIX"
             print*, "condition number=", (SV(size(SV))/SV(1))
             do j = 1, N
                !! perturbation coefficient should be taken from diamond options with a check for default value
                call random_number(perturb(1))
                perturb(1)=-1.0+2.0*perturb(1)
                print*, "perturb(1)=", N, perturb(1)
                abscissa_val_at_quad(i,j) = abscissa_val_at_quad(i,j)*(1.0+perturb(1)*perturb_val)
                print*, "abscissa", j, " = ", abscissa_val_at_quad(i,j)
             end do
             A = A_matrix(abscissa_val_at_quad)
             call svd(A(i,:,:), svd_tmp1, SV, svd_tmp2)
             print*, "perturbed condition number=", (SV(size(SV))/SV(1))
             p_num=p_num+1
             print*, "p_num=", p_num
!             if (SV(size(SV))/SV(1) < cond) then
!                print*, "perturbing does not help"
!             end if
!          end if
          end do
       end do
    else if (singular_option=='set_source_to_zero') then
       do i = 1, ele_ngi(abscissa(1), ele)
          call svd(A(i,:,:), svd_tmp1, SV, svd_tmp2)
          if (SV(size(SV))/SV(1) < cond) then
             ewrite(2,*) 'ill-conditioned matrix found', SV(size(SV))/SV(1)
             A(i,:,:) = 0.0
             A_3(i,:,:) = 0.0
             C(i,:) = 0.0
             do j = 1, 2*N
                A(i,j,j) = 1.0
             end do
          end if
       end do
    end if
!!       p_num = 1


!!       do j = 1, N
!!          absc_init(j) = abscissa_val_at_quad(i,j)
!!       end do

!!       do while (SV(size(SV))/SV(1) < cond)
!          call get_option("/timestepping/current_time", curr_time)   ! get the current simulation time
          !! a better check for the first timestep instead of small value of the current time
 !         if (curr_time < 1e-7) then   ! this block perturbates the abscissas in A matrix if current time is 0
!!             print*, "perturbing %g \n", p_num
!!             do j = 1, N
                !! perturbation coefficient should be taken from diamond options with a check for default value
!!                call random_number(perturb(1))
!!                abscissa_val_at_quad(i,j) = absc_init(j) + perturb(1)
!                abscissa_val_at_quad(:,j) = abscissa_val_at_quad(:,j) + (2-j)*(0.1) 
!!                print*, "abscissa", j, " = ", abscissa_val_at_quad(i,j)
!!             end do
!!             A = A_matrix(abscissa_val_at_quad)
!!             call svd(A(i,:,:), svd_tmp1, SV, svd_tmp2)
!!             print*, "condition number=", (SV(size(SV))/SV(1))
!!             p_num = p_num+1
!             print*, "perturbed condition number=", (SV(size(SV))/SV(1))
!             if (SV(size(SV))/SV(1) < cond) then
!                print*, "perturbing does not help"
!             end if

!          else          
!             ewrite(2,*) 'ill-conditioned matrix found'
!             A(i,:,:) = 0.0
!             A_3(i,:,:) = 0.0
!             C(i,:) = 0.0
!             do j = 1, 2*N
!                A(i,j,j) = 1.0
!             end do
!          end if
!       end if
!!       end do
!    end do
           



    ! solve linear system to find source values
    do i = 1, ele_ngi(abscissa(1), ele)
       b(:,1) = matmul(A_3(i,:,:), C(i,:)) + S_rhs(i,:)   !! gb 15-11-2012! added S_rhs term
       call dqmom_solve(A(i,:,:), b, stat)
       weight_S_at_quad(i,:) = b(:N,1)
       abscissa_S_at_quad(i,:) = b(N+1:,1)
    end do    

    ! In the DG case we apply the inverse mass locally.
    invmass = inverse(shape_shape(shape, shape, detwei))

    ! integrate and add to source fields
    do i = 1, N
       weight_S_at_nodes = shape_rhs(shape, detwei* weight_S_at_quad(:,i))
       abscissa_S_at_nodes = shape_rhs(shape, detwei* abscissa_S_at_quad(:,i))
       if(continuity(abscissa(1))<0) then
          weight_S_at_nodes = matmul(weight_S_at_nodes, invmass)
          abscissa_S_at_nodes = matmul(abscissa_S_at_nodes, invmass)
       end if
       call addto(s_weight(i)%ptr, nodes, weight_S_at_nodes)
       call addto(s_weighted_abscissa(i)%ptr, nodes, abscissa_S_at_nodes)
    end do

  end subroutine dqmom_calculate_source_term_ele

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
    real :: mean, std
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
                do j = 1, node_count(stats)
                   call set(stats, j, node_val(moments(4)%ptr,j)/node_val(moments(3)%ptr,j))
                end do
             end if

          end do

          deallocate(moments)

       end if
    end do
    ewrite(1, *) "Exiting dqmom_calculate_statistics"


  end subroutine dqmom_calculate_statistics

  subroutine dqmom_check_options(state)

    type(state_type), intent(in) :: state
    
    integer :: i_pop
    integer :: n_abscissa, n_weights, n_weighted_abscissa, n_moments
    character(len=OPTION_PATH_LEN) :: option_path 

    write(*,*) 'in dqmom_check_options'

    do i_pop = 1, option_count(trim(state%option_path)//'/population_balance')
       option_path = trim(state%option_path)//'/population_balance['//int2str(i_pop)//']'

       ! Check there are the same number of abscissa, weight and weighted abscissa fields
       n_abscissa = option_count(trim(option_path)//'/abscissa/scalar_field')
       n_weights = option_count(trim(option_path)//'/weights/scalar_field')
       n_weighted_abscissa = option_count(trim(option_path)// &
            '/weighted_abscissa/scalar_field')
       assert(n_weights == n_abscissa) 
       assert(n_weights == n_weighted_abscissa)

       ! Check there are sufficient abscissa to calculate requested moments
       n_moments = option_count(trim(option_path)//'/moments/scalar_field')
       assert((n_moments / 2) == n_abscissa)       

       ! Need to check that all fields are on the same mesh
    end do

  end subroutine dqmom_check_options

end module dqmom
