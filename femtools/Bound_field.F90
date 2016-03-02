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

module bound_field_module

  use FLDebug
  use spud
  use global_parameters, only : FIELD_NAME_LEN, real_8
  use quicksort
  use parallel_tools
  use sparse_tools
  use transform_elements
  use fetools
  use parallel_fields
  use fields
  use field_options
  use sparse_matrices_fields
  use sparsity_patterns
  use halos
  implicit none

  type(scalar_field), save :: func_target_field, func_lumped_mass
  real, save :: integral_value
  type(csr_matrix), pointer, save :: func_mass_matrix
  real, dimension(:, :), allocatable, save :: func_detwei
  
  integer, save :: functional_type
  integer, parameter, public :: FUNCTIONAL_VEC_L2=0, &
                                FUNCTIONAL_LUMPED_VEC_L2=1, &
                                FUNCTIONAL_FUNC_L2=2

  private 

  public :: bound_field, bound_field_diffuse

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
  
    call get_option(trim(complete_field_path(field%option_path, stat=statp))// &
                  & "/galerkin_projection/continuous/bounded[0]/name", &
                    algorithm_type, default="Diffuse")
                    ! have to give it a default as this is unit tested
  
    select case(algorithm_type)
    case("Diffuse")
      call bound_field_diffuse(field, max_bound, min_bound, &
                         mass, lumpedmass, inverse_lumpedmass)
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

  call get_option(trim(complete_field_path(field%option_path, stat=statp))// &
                  & "/galerkin_projection/continuous/bounded[0]/boundedness_iterations", &
                  & iters, default=1000)
          
  call get_option(trim(complete_field_path(field%option_path, stat=statp))// &
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
  
  if(have_option(trim(complete_field_path(field%option_path, stat=statp))// &
    & "/galerkin_projection/continuous/bounded[0]/repair_deviations")) then
                
      call get_option(trim(complete_field_path(field%option_path, stat=statp))// &
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

end module bound_field_module
