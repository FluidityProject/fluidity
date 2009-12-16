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

module bound_field_module

  use FLDebug
  use fields
  use spud
  use global_parameters, only : FIELD_NAME_LEN, real_8
  use field_options
  use sparse_matrices_fields
  use quicksort
  use sparse_tools
  use sparsity_patterns
  use transform_elements
  use fetools
  use halos
  use parallel_fields
  use parallel_tools
  implicit none

  type(scalar_field), save :: func_target_field, func_lumped_mass
  real, save :: integral_value
  type(csr_matrix), pointer, save :: func_mass_matrix
  real, dimension(:, :), allocatable, save :: func_detwei
  
  integer, save :: functional_type
  integer, parameter, public :: FUNCTIONAL_VEC_L2=0, &
                                FUNCTIONAL_LUMPED_VEC_L2=1, &
                                FUNCTIONAL_FUNC_L2=2

  interface
    subroutine algencan(epsfeas, epsopt, iprint, ncomp, n, x, l, u, m, &
      & lambda, equatn, linear, coded, checkder, f, cnorm, snorm, nlpsupn, inform)
      use global_parameters, only : real_8
      implicit none
      real(kind = real_8) :: epsfeas
      real(kind = real_8) :: epsopt
      integer :: iprint
      integer :: ncomp
      integer :: n
      real(kind = real_8), dimension(n) :: x
      real(kind = real_8), dimension(n) :: l
      real(kind = real_8), dimension(n) :: u
      integer :: m
      real(kind = real_8), dimension(m) :: lambda
      logical, dimension(m) :: equatn
      logical, dimension(m) :: linear
      logical, dimension(10) :: coded
      logical :: checkder
      real(kind = real_8) :: f
      real(kind = real_8) :: cnorm
      real(kind = real_8) :: snorm
      real(kind = real_8) :: nlpsupn
      integer :: inform
    end subroutine algencan
  end interface

  public :: bound_field, bound_field_algencan

  contains
  
  subroutine bound_field(field, max_bound, min_bound, &
                         mass, lumpedmass, inverse_lumpedmass, &
                         bounded_soln, positions)
  
    type(scalar_field), intent(inout) :: field
    type(scalar_field), intent(in) :: max_bound, min_bound
    type(csr_matrix), intent(in) :: mass
    type(scalar_field), intent(inout) :: lumpedmass
    type(scalar_field), intent(in) :: inverse_lumpedmass, bounded_soln
    type(vector_field), intent(in) :: positions
    
    character(len=FIELD_NAME_LEN) :: algorithm_type
    integer :: statp
  
    call get_option(trim(complete_field_path(field%option_path, statp))// &
                  & "/galerkin_projection/continuous/bounded[0]/name", &
                    algorithm_type, default="Diffuse")
                    ! have to give it a default as this is unit tested
  
    select case(algorithm_type)
    case("Diffuse")
      call bound_field_diffuse(field, max_bound, min_bound, &
                         mass, lumpedmass, inverse_lumpedmass)
    case("Algencan")
      call bound_field_algencan(field, positions, max_bound, min_bound, lumpedmass, &
                          bounded_soln=bounded_soln, mass_matrix=mass)
    case default
      FLAbort("Unknown bounding algorithm selected")
    end select
  
  end subroutine bound_field

  subroutine bound_field_diffuse(field, max_bound, min_bound, &
                         mass, lumpedmass, inverse_lumpedmass)
  
  type(scalar_field), intent(inout) :: field
  type(scalar_field), intent(in) :: max_bound, min_bound
  type(csr_matrix), intent(in) :: mass
  type(scalar_field), intent(in) :: lumpedmass, inverse_lumpedmass  
  
  integer :: j, k, iters, node, statp
  real :: bound_tol, repair_tol
  real :: field_node_val, max_bound_node_val, min_bound_node_val
  real :: lumpedmass_node_val, lump_ratio, target_change
  integer :: target_node
  type(scalar_field) :: deviation, diffused, increase_capacity, decrease_capacity
  integer, dimension(:), allocatable :: increase_capacity_idx, decrease_capacity_idx, &
                                        deviation_idx
  real :: maxdeviation, mindeviation
  real :: local_decrease_cap, local_increase_cap, local_deviation
  real :: total_decrease_cap, total_increase_cap, total_deviation

  ewrite(1,*) 'In bound_field_diffuse'

  call allocate(deviation, field%mesh, "BoundedDeviation")
  call allocate(diffused, field%mesh, "DiffusedDeviation")

  call get_option(trim(complete_field_path(field%option_path, statp))// &
                  & "/galerkin_projection/continuous/bounded[0]/boundedness_iterations", &
                  & iters, default=1000)
          
  call get_option(trim(complete_field_path(field%option_path, statp))// &
                  & "/galerkin_projection/continuous/bounded[0]/boundedness_iterations/tolerance", &
                  & bound_tol, default=epsilon(0.0))
          
  ! Iterate!
  iterloop: do k=1,iters

    ! Step 1. Compute deviation
    call zero(deviation)
    do node=1,node_count(field)
      field_node_val = node_val(field, node)
      max_bound_node_val = node_val(max_bound, node)
      min_bound_node_val = node_val(min_bound, node)
      if (field_node_val > max_bound_node_val) then
        call set(deviation, node, field_node_val - max_bound_node_val)
      else if (field_node_val < min_bound_node_val) then
        call set(deviation, node, field_node_val - min_bound_node_val)
      end if
    end do

    maxdeviation = maxval(deviation)
    call allmax(maxdeviation)
    mindeviation = minval(deviation)
    call allmin(mindeviation)

    if (maxdeviation < bound_tol .and. mindeviation > -bound_tol) then
      exit iterloop
    end if

    ! Step 2. Diffuse!
    call mult(diffused, mass, deviation)
    call scale(diffused, inverse_lumpedmass)

    ! Step 3. Update!
    call addto(diffused, deviation, -1.0)
    call addto(field, diffused, 2.0)

    call halo_update(field, 1)
  end do iterloop
  
  ewrite(2,*) "Bounded interpolation for field ", trim(field%name), " took ", k, " iterations"
  ewrite(2,*) "Before final diffusion:"
  ewrite(2,*) "maxval(deviation): ", maxval(deviation), "; minval(deviation): ", minval(deviation)
  
  if(have_option(trim(complete_field_path(field%option_path, statp))// &
    & "/galerkin_projection/continuous/bounded[0]/repair_deviations")) then
                
      call get_option(trim(complete_field_path(field%option_path, statp))// &
                  & "/galerkin_projection/continuous/bounded[0]/repair_deviations/tolerance", &
                  & repair_tol, default=epsilon(0.0))

      ! If we reached the maximum number of iterations or we still have deviations, then the time has come
      ! to perform surgery. We shuffle the remaining deviation (hopefully small!)
      ! around without regard to physical locality. Hopefully the user has
      ! specified enough iterations so that this is 1e-10 range ..
      if ((k>iters) .or. (maxdeviation > repair_tol) .or. (mindeviation < -repair_tol)) then
        ! Here, increase_capacity is the capacity to absorb increases,
        ! and decrease_capacity is the capacity to absorb decreases.
        ewrite(2,*) 'Repairing deviations'
        
        call allocate(increase_capacity, field%mesh, "IncreaseAbsorptionCapacity")
        call allocate(decrease_capacity, field%mesh, "DecreaseAbsorptionCapacity")
        allocate(increase_capacity_idx(node_count(field)))
        allocate(decrease_capacity_idx(node_count(field)))
        allocate(deviation_idx(node_count(field)))
      
        call zero(deviation)
        call zero(increase_capacity)
        call zero(decrease_capacity)
        do node=1,node_count(field)
          if(.not.node_owned(field,node)) cycle
          field_node_val = node_val(field, node)
          max_bound_node_val = node_val(max_bound, node)
          min_bound_node_val = node_val(min_bound, node)
          if (field_node_val > max_bound_node_val) then
            call set(deviation, node, field_node_val - max_bound_node_val)
            ! here we intentionally ignore the full decrease_capacity of this node 
            ! (which would be = -min_bound_node_val + field_node_val) in order
            ! to get the sharpest interpolation possible... i.e. overshoots are only
            ! taken down to the max_bound rather than allowing their full capacity down
            ! to the min_bound to be used (this could lead to a solution that's more diffuse
            ! than the lumped bounded solution)
            call set(decrease_capacity, node, field_node_val - max_bound_node_val)
          else if (field_node_val < min_bound_node_val) then
            call set(deviation, node, field_node_val - min_bound_node_val)
            ! here we intentionally ignore the full increase_capacity of this node 
            ! (which would be = max_bound_node_val - field_node_val) in order
            ! to get the sharpest interpolation possible... i.e. overshoots are only
            ! taken down to the max_bound rather than allowing their full capacity down
            ! to the min_bound to be used (this could lead to a solution that's more diffuse
            ! than the lumped bounded solution)
            call set(increase_capacity, node, -(field_node_val - min_bound_node_val))
          else
            call set(increase_capacity, node, max_bound_node_val - field_node_val)
            call set(decrease_capacity, node,-min_bound_node_val + field_node_val)
          end if
        end do
        
#ifdef DDEBUG
        local_increase_cap = sum(increase_capacity%val*lumpedmass%val)
        local_decrease_cap = sum(decrease_capacity%val*lumpedmass%val)
        local_deviation = sum(abs(deviation%val)*lumpedmass%val)
        ewrite(2,*) 'local_deviation = ', local_deviation
        ewrite(2,*) 'local_decrease_cap+local_increase_cap = ', local_decrease_cap+local_increase_cap

        if(isparallel()) then
          total_increase_cap = local_increase_cap 
          call allsum(total_increase_cap)
          total_decrease_cap = local_decrease_cap 
          call allsum(total_decrease_cap)
          total_deviation = local_deviation
          call allsum(total_deviation)

          ewrite(2,*) 'total_deviation = ', total_deviation
          ewrite(2,*) 'total_decrease_cap+total_increase_cap = ', total_decrease_cap+total_increase_cap
        end if

        if(local_deviation>local_decrease_cap+local_increase_cap) then
          ewrite(0,*) 'local_deviation = ', local_deviation
          ewrite(0,*) 'local_decrease_cap+local_increase_cap = ', local_decrease_cap+local_increase_cap

          if(isparallel()) then 
            ewrite(0,*) 'total_deviation = ', total_deviation
            ewrite(0,*) 'total_decrease_cap+total_increase_cap = ', total_decrease_cap+total_increase_cap

            if(total_deviation<=total_increase_cap+total_decrease_cap) then
              ewrite(0,*) 'Enough capacity globally... should parallelise the repairing of deviations'
            else
              ewrite(0,*) 'Not enough capacity globally either... oh dear.'
            end if
          end if

          ewrite(0,*) "Warning: insufficient capacity to repair deviations.  Will try my best."
        end if
#endif

        ! sort into increasing order
        call qsort(increase_capacity%val, increase_capacity_idx)
        call qsort(decrease_capacity%val, decrease_capacity_idx)
        call qsort(abs(deviation%val), deviation_idx)
        
        ! Now we can distribute stuff around.

        ! starting with the largest deviation
        deviation_loop: do j=node_count(field),1,-1
          node = deviation_idx(j)
          if(.not.node_owned(field, node)) cycle
          lumpedmass_node_val = node_val(lumpedmass, node)
          k = node_count(field)+1
          do while((node_val(deviation, node) > repair_tol).and.&
                   (k>1))
            k = k - 1
            ! Choose the node with the largest capacity
            ! (at least at the start since we don't reorder once we get going)...
            target_node = increase_capacity_idx(k)
            if(.not.node_owned(field,target_node)) cycle
            if(node_val(increase_capacity, target_node)<epsilon(0.0)) cycle
            ! If we take x from a node here and put it at a node there,
            ! we break conservation. Instead, we need to multiply by the ratio
            ! of the lumped mass .. ugly .. I know .. I'm sorry.
            lump_ratio = node_val(lumpedmass, target_node) / lumpedmass_node_val

            ! Obviously, we would /like/ to take the whole deviation and place it
            ! at target_node. But we might not be able to. 
            target_change = min(node_val(increase_capacity, target_node), node_val(deviation, node) / lump_ratio)
            
            call addto(field, target_node, target_change)
            call addto(field, node, -target_change * lump_ratio)
            
            call addto(deviation, node, -target_change * lump_ratio)
            call addto(decrease_capacity, node, -target_change * lump_ratio)
            
            field_node_val = node_val(field, target_node)
            max_bound_node_val = node_val(max_bound, target_node)
            min_bound_node_val = node_val(min_bound, target_node)
            if (field_node_val > max_bound_node_val) then
              call set(deviation, target_node, field_node_val - max_bound_node_val)
              call set(increase_capacity, target_node, 0.0)
              ! see note above about setting of the decrease_capacity
              call set(decrease_capacity, node, field_node_val - max_bound_node_val)
            else if (field_node_val < min_bound_node_val) then
              call set(deviation, target_node, field_node_val - min_bound_node_val)
              ! see note above about setting of the increase_capacity
              call set(increase_capacity, node, -(field_node_val - min_bound_node_val))
              call set(decrease_capacity, target_node, 0.0)
            else
              call set(deviation, target_node, 0.0)
              call set(increase_capacity, target_node, max_bound_node_val - field_node_val)
              call set(decrease_capacity, target_node,-min_bound_node_val + field_node_val)
            end if
          end do
          
          k = node_count(field)
          do while((node_val(deviation, node) < -repair_tol).and.&
                   (k>1))
            k = k - 1
            target_node = decrease_capacity_idx(k)
            if(.not.node_owned(field, target_node)) cycle
            if(node_val(decrease_capacity, target_node)<epsilon(0.0)) cycle

            lump_ratio = node_val(lumpedmass, target_node) / lumpedmass_node_val

            ! target_change is positive -- we will be subtracting it from target_node
            target_change = min(node_val(decrease_capacity, target_node), -node_val(deviation, node) / lump_ratio)

            call addto(field, target_node, -target_change)
            call addto(field, node, target_change * lump_ratio)
            
            call addto(deviation, node, target_change * lump_ratio)
            call addto(increase_capacity, node, -target_change * lump_ratio)
            
            field_node_val = node_val(field, target_node)
            max_bound_node_val = node_val(max_bound, target_node)
            min_bound_node_val = node_val(min_bound, target_node)
            if (field_node_val > max_bound_node_val) then
              call set(deviation, target_node, field_node_val - max_bound_node_val)
              call set(increase_capacity, target_node, 0.0)
              ! see note above about setting of the decrease_capacity
              call set(decrease_capacity, node, field_node_val - max_bound_node_val)
            else if (field_node_val < min_bound_node_val) then
              call set(deviation, target_node, field_node_val - min_bound_node_val)
              ! see note above about setting of the increase_capacity
              call set(increase_capacity, node, -(field_node_val - min_bound_node_val))
              call set(decrease_capacity, target_node, 0.0)
            else
              call set(deviation, target_node, 0.0)
              call set(increase_capacity, target_node, max_bound_node_val - field_node_val)
              call set(decrease_capacity, target_node,-min_bound_node_val + field_node_val)
            end if
          end do

        end do deviation_loop
      
        call deallocate(increase_capacity)
        call deallocate(decrease_capacity)
        deallocate(increase_capacity_idx)
        deallocate(decrease_capacity_idx)
        deallocate(deviation_idx)
    
      end if
    
    end if
    
    call deallocate(diffused)
    call deallocate(deviation)

    ! have only done 1st order halo updates until now, so finish with a full update
    call halo_update(field)

  end subroutine bound_field_diffuse

  subroutine bound_field_algencan(field, positions, ub, lb, lumped_mass, bounded_soln, mass_matrix)
    type(scalar_field), intent(inout) :: field, lumped_mass
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: ub, lb
    type(scalar_field), intent(in), optional :: bounded_soln
    type(csr_matrix), intent(in), optional, target :: mass_matrix

    integer :: prob_size

    ! We only allow one constraint - conservation
    integer, parameter :: no_constraints=1
    real(kind = real_8), dimension(no_constraints) :: lambda
    logical, dimension(no_constraints) :: equatn, linear
    logical, dimension(10) :: coded
    logical :: checkder
    real(kind = real_8) :: epsfeas, epsopt
    integer :: iprint, ncomp
    real(kind = real_8) :: cnorm, snorm, nlpsupn, f
    integer :: inform
    
    integer :: ele
    type(csr_sparsity) :: mass_matrix_sparsity

    character(len=FIELD_NAME_LEN) :: functional_name
    real :: weight
    integer :: statp

#ifndef DOUBLEP
    real(kind = real_8), dimension(:), allocatable :: x
#endif

#ifdef HAVE_ALGENCAN

    ewrite(1,*) 'In bound_field_algencan'
    
    call get_option(trim(complete_field_path(field%option_path, statp))// &
                  & "/galerkin_projection/continuous/bounded[0]/functional[0]/name", &
                    functional_name, default="L2")
                    
    select case(functional_name)
    case("L2")
      functional_type = FUNCTIONAL_VEC_L2
      ewrite(2,*) "using functional L2"
    case("LumpedMassL2")
      functional_type = FUNCTIONAL_LUMPED_VEC_L2
      ewrite(2,*) "using functional LumpedMassL2"
    case("IntegralL2")
      functional_type = FUNCTIONAL_FUNC_L2
      ewrite(2,*) "using functional IntegralL2"
    case default
      FLAbort("Unrecognised functional type.")
    end select

    call get_option(trim(complete_field_path(field%option_path, statp))// &
                  & "/galerkin_projection/continuous/bounded[0]/functional[0]/weight", &
                    weight, default=1.0)
    
    assert(.not. associated(func_target_field%val))
    call allocate(func_target_field, field%mesh, "TemporaryTargetingField")
    call allocate(func_lumped_mass, lumped_mass%mesh, "TemporaryLumpedMass")
    call set(func_target_field, field)
    call scale(func_target_field, weight)
    call set(func_lumped_mass, lumped_mass)
    integral_value = dot_product(lumped_mass%val, field%val)

    if(functional_type==FUNCTIONAL_FUNC_L2) then
      allocate(func_detwei(ele_count(func_target_field), ele_ngi(func_target_field, 1)))
      if (present(mass_matrix)) then
        func_mass_matrix => mass_matrix
      else
        mass_matrix_sparsity = make_sparsity(func_target_field%mesh, func_target_field%mesh, "MassMatrixSparsityTemporary")
        allocate(func_mass_matrix)
        call allocate(func_mass_matrix, mass_matrix_sparsity, name="MassMatrixTemporary")
        call zero(func_mass_matrix)
      end if
      do ele=1,ele_count(func_target_field)
        call transform_to_physical(positions, ele, &
                                &  detwei=func_detwei(ele, :))
        if (.not. present(mass_matrix)) then
          call addto(func_mass_matrix, ele_nodes(func_target_field, ele), ele_nodes(func_target_field, ele), &
                    shape_shape(ele_shape(func_target_field, ele), ele_shape(func_target_field, ele), func_detwei(ele, :)))
        end if
  
      end do
    end if

    if (present(bounded_soln)) then
      call addto(func_target_field, bounded_soln)
      call scale(func_target_field, 0.5)
    end if

    !.. algencan stuff

    ! Initial guess for Lagrange multiplier
    lambda = 0.0

    ! Are the constraints equality or less-than relations?
    ! .true. => equality
    ! .false. => less-than
    ! We want to enforce mass conservation, which is an equality constraint
    equatn(1) = .true.

    ! Are the constraints linear?
    ! Mass conservation is indeed a linear constraint
    linear(1) = .true.

    ! What routines do we have?
    coded( 1) = .true.  ! evalf    
    coded( 2) = .true.  ! evalg   
    coded( 3) = .true.  ! evalh                                                                                                                            
    coded( 4) = .true.  ! evalc   
    coded( 5) = .true.  ! evaljac
    coded( 6) = .false.  ! evalhc  
    coded( 7) = .false. ! evalfc
    coded( 8) = .false. ! evalgjac
    coded( 9) = .false. ! evalhl
    coded(10) = .false. ! evalhlp

    checkder = .false.
    prob_size = node_count(field)

    ! constraint tolerance
    epsfeas = 1.0d-12
    ! gradient tolerance
    epsopt = 1.0d-05
    iprint = 0
    ncomp = 5

#ifdef DOUBLEP
    call algencan(epsfeas, epsopt, iprint, ncomp, prob_size, field%val, lb%val, ub%val, no_constraints, &
         & lambda, equatn, linear, coded, checkder, f, cnorm, snorm, nlpsupn, inform)
#else
    allocate(x(node_count(field)))
    x = field%val
    call algencan(epsfeas, epsopt, iprint, ncomp, prob_size, x, real(lb%val, kind = real_8), real(ub%val, kind = real_8), no_constraints, &
         & lambda, equatn, linear, coded, checkder, f, cnorm, snorm, nlpsupn, inform)
    field%val = x
    deallocate(x)
#endif

    call deallocate(func_lumped_mass)
    call deallocate(func_target_field)

    if(functional_type==FUNCTIONAL_FUNC_L2) then
      deallocate(func_detwei)
      if (.not. present(mass_matrix)) then
        call deallocate(mass_matrix_sparsity)
        call deallocate(func_mass_matrix)
        deallocate(func_mass_matrix)
      end if
    endif

#else
    FLAbort("You need algencan to bound a field with algencan.")
#endif
  end subroutine bound_field_algencan

end module bound_field_module

! deliberately external so algencan can find it
subroutine evalf(n, x, f, flag)
  use fields
  use bound_field_module
  use vector_tools
  use fetools
  implicit none
  integer, intent(in) :: n
  real, dimension(n), intent(in) :: x
  real, intent(out) :: f
  integer, intent(out) :: flag
  integer :: i
  
  integer :: ele
  integer, dimension(:), pointer :: nodes
  real, dimension(:), allocatable :: integrand

  flag = 0

  select case(functional_type)
  case(FUNCTIONAL_VEC_L2)
    f = norm2(func_target_field%val - x)**2
  case(FUNCTIONAL_LUMPED_VEC_L2)
    f = 0
    do i=1,n
      f = f + (func_lumped_mass%val(i) * (func_target_field%val(i) - x(i)))**2
    end do
  case(FUNCTIONAL_FUNC_L2)
    allocate(integrand(ele_ngi(func_target_field, 1)))
    f = 0
    do ele=1,ele_count(func_target_field)
      nodes => ele_nodes(func_target_field, ele)
      integrand = matmul((ele_val(func_target_field, ele) - x(nodes))**2, func_target_field%mesh%shape%n)
      f = f + dot_product(integrand, func_detwei(ele, :))
    end do
    deallocate(integrand)
  end select

end subroutine evalf

subroutine evalg(n, x, g, flag)
  use fields
  use bound_field_module
  use fetools
  implicit none
  integer, intent(in) :: n
  real, dimension(n), intent(in) :: x
  real, dimension(*) :: g
  integer, intent(out) :: flag
  
  integer :: ele
  integer, dimension(:), pointer :: nodes
  real, dimension(:), allocatable :: integrand
  
  flag = 0
  select case(functional_type)
  case(FUNCTIONAL_VEC_L2)
    g(1:n) = -2 * (func_target_field%val - x)
  case(FUNCTIONAL_LUMPED_VEC_L2)
    g(1:n) = -2 * func_lumped_mass%val * (func_target_field%val - x)
  case(FUNCTIONAL_FUNC_L2)
    allocate(integrand(ele_ngi(func_target_field, 1)))
    g(1:n) = 0
    do ele=1,ele_count(func_target_field)
      nodes => ele_nodes(func_target_field, ele)
      integrand = matmul(-2 * (ele_val(func_target_field, ele) - x(nodes)), func_target_field%mesh%shape%n)
      g(nodes) = g(nodes) + shape_rhs(ele_shape(func_target_field, ele), integrand * func_detwei(ele, :))
    end do
    deallocate(integrand)
  end select

end subroutine evalg

subroutine evalh(n, x, hlin, hcol, hval, hnnz, flag)
  use bound_field_module
  implicit none
  integer, intent(in) :: n
  real, dimension(n), intent(in) :: x
  real, dimension(*) :: hval
  integer, dimension(*) :: hlin, hcol
  integer, intent(out) :: hnnz, flag
  integer :: i, j

  integer :: tmp_nnz
  real, dimension(size(func_mass_matrix%val)) :: tmp_hval
  integer, dimension(size(func_mass_matrix%val)) :: tmp_hcol, tmp_hlin

  flag = 0

  select case(functional_type)
  case(FUNCTIONAL_VEC_L2)
    hnnz = n
    hval(1:n) = 2.0
    hlin(1:n) = (/ (i,i=1,n) /)
    hcol(1:n) = hlin(1:n)
  case(FUNCTIONAL_LUMPED_VEC_L2)
    hnnz = n
    hval(1:n) = 2.0 * func_lumped_mass%val
    hlin(1:n) = (/ (i,i=1,n) /)
    hcol(1:n) = hlin(1:n)
  case(FUNCTIONAL_FUNC_L2)
    tmp_nnz = size(func_mass_matrix%val)
    tmp_hval(1:tmp_nnz) = func_mass_matrix%val
    do i=1,size(func_mass_matrix%sparsity%findrm)-1
      tmp_hlin(func_mass_matrix%sparsity%findrm(i):func_mass_matrix%sparsity%findrm(i)-1) = i
    end do
    tmp_hcol(1:tmp_nnz) = func_mass_matrix%sparsity%colm
  
    hnnz = (tmp_nnz + n)/2
    j = 1
    do i=1,tmp_nnz
      if (tmp_hlin(i) .ge. tmp_hcol(i)) then
        hval(j) = tmp_hval(i)
        hlin(j) = tmp_hlin(i)
        hcol(j) = tmp_hcol(i)
        j = j + 1
      end if
    end do
  end select

end subroutine evalh

subroutine evalc(n, x, ind, c, flag)
  use fields
  use bound_field_module
  implicit none
  integer, intent(in) :: n, ind
  real, dimension(n), intent(in) :: x
  real, intent(out) :: c
  integer, intent(out) :: flag

  assert(ind == 1)
  flag = 0

  c = (dot_product(func_lumped_mass%val, x) - integral_value)/integral_value
end subroutine evalc

subroutine evaljac(n, x, ind, jcvar, jcval, jcnnz, flag)
  use fields
  use bound_field_module
  implicit none
  integer, intent(in) :: n, ind
  real, dimension(n), intent(in) :: x
  integer, intent(out) :: flag, jcnnz
  integer, dimension(*) :: jcvar
  real, dimension(*) :: jcval
  integer :: i

  assert(ind == 1)
  flag = 0

  jcnnz = n
  jcvar(1:n) = (/ (i,i=1,n) /)
  jcval(1:n) = func_lumped_mass%val/integral_value
end subroutine evaljac

subroutine evalhc(n, x, ind, hclin, hccol, hcval, hcnnz, flag)
  implicit none
  integer, intent(in) :: n, ind
  real, dimension(n), intent(in) :: x
  integer, intent(out) :: flag, hcnnz
  integer, dimension(*) :: hclin, hccol
  real, dimension(*) :: hcval

  flag = 0
  hcnnz = 0
end subroutine evalhc

!     ******************************************************************
!     ******************************************************************

subroutine evalfc(n,x,f,m,c,flag)

implicit none

!     SCALAR ARGUMENTS
integer flag,m,n
double precision f

!     ARRAY ARGUMENTS
double precision c(m),x(n)

flag = - 1

end subroutine

!     ******************************************************************
!     ******************************************************************

subroutine evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,flag)

implicit none

!     SCALAR ARGUMENTS
integer flag,jcnnz,m,n

!     ARRAY ARGUMENTS
integer jcfun(*),jcvar(*)
double precision g(n),jcval(*),x(n)

flag = - 1

end subroutine

!     ******************************************************************
!     ******************************************************************

subroutine evalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

implicit none

!     SCALAR ARGUMENTS
logical goth
integer flag,m,n
double precision sf

!     ARRAY ARGUMENTS
double precision hp(n),lambda(m),p(n),sc(m),x(n)

!     This subroutine might compute the product of the Hessian of the
!     Lagrangian times vector p (just the Hessian of the objective 
!     function in the unconstrained or bound-constrained case). 
!     
!     Parameters of the subroutine:
!
!     On Entry:
!
!     n        integer,
!              number of variables,
!
!     x        double precision x(n),
!              current point,
!
!     m        integer,
!              number of constraints,
!
!     lambda   double precision lambda(m),
!              vector of Lagrange multipliers,
!
!     p        double precision p(n),
!              vector of the matrix-vector product,
!
!     goth     logical,
!              can be used to indicate if the Hessian matrices were
!              computed at the current point. It is set to .false.
!              by the optimization method every time the current
!              point is modified. Sugestion: if its value is .false. 
!              then compute the Hessians, save them in a common 
!              structure and set goth to .true.. Otherwise, just use 
!              the Hessians saved in the common block structure,
!
!     On Return
!
!     hp       double precision hp(n),
!              Hessian-vector product,
!
!     goth     logical,
!              see above,
!              
!     flag     integer,
!              You must set it to any number different of 0 (zero) if 
!              some error ocurred during the evaluation of the 
!              Hessian-vector product. (For example, trying to compute 
!              the square root of a negative number, dividing by zero 
!              or a very small number, etc.) If everything was o.k. you 
!              must set it equal to zero.

flag = - 1

end subroutine

subroutine evalhl(n,x,m,lambda,scalef,scalec,hllin,hlcol,hlval, &
hlnnz,flag)

implicit none

!     SCALAR ARGUMENTS
integer flag,hlnnz,m,n
double precision scalef

!     ARRAY ARGUMENTS
integer hlcol(*),hllin(*)
double precision hlval(*),lambda(m),scalec(m),x(n)

flag = - 1

end subroutine
