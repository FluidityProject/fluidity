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

module multimaterial_module
  !! This module contains the options and material properties used
  !! when running a multimaterial simulation.
  use fldebug
  use state_module
  use fields
  use spud
  use fefields, only: compute_cv_mass
  use global_parameters, only: OPTION_PATH_LEN
  use field_priority_lists
  use field_options
  use diagnostic_fields_matrices
  use cv_upwind_values
  use equation_of_state, only: compressible_material_eos
  implicit none

  interface calculate_bulk_property
    module procedure calculate_bulk_scalar_property, calculate_bulk_vector_property, calculate_bulk_tensor_property
  end interface

  interface add_scaled_material_property
    module procedure add_scaled_material_property_scalar, add_scaled_material_property_vector, add_scaled_material_property_tensor
  end interface

  private
  public :: initialise_diagnostic_material_properties, &
            calculate_material_mass, calculate_bulk_material_pressure, &
            calculate_sum_material_volume_fractions, calculate_material_volume, &
            calculate_bulk_property, add_scaled_material_property, calculate_surfacetension, &
            calculate_diagnostic_material_volume_fraction, order_states_priority

contains

  subroutine calculate_surfacetension(state, surfacetension)
    ! calculates the surface tension in tensor form
    type(state_type), dimension(:), intent(inout) :: state
    type(tensor_field), intent(inout) :: surfacetension
    
    type(scalar_field), pointer :: volumefraction
    integer :: i, node, dimi, dimj, stat
    logical :: prognostic
    type(scalar_field) :: grad_mag, grad_mag2
    type(vector_field) :: gradient
    
    type(vector_field) :: normals
    logical, dimension(:), allocatable :: on_boundary
    type(vector_field), pointer :: x
    integer, dimension(2) :: shape_option
    integer, dimension(:), allocatable :: surface_ids
    
    real :: coeff, eq_angle
    real, dimension(surfacetension%dim(1), surfacetension%dim(2)) :: tensor
    
    if(size(state)==1) then
      FLExit("Don't know how to calculate a surface tension with only one material_phase.")
    end if
    
    x => extract_vector_field(state(1), "Coordinate")
    
    call allocate(gradient, surfacetension%dim(1), surfacetension%mesh, "Gradient")
    gradient%option_path = surfacetension%option_path
    
    call zero(surfacetension)
    
    do i = 1, size(state)
      volumefraction => extract_scalar_field(state(i), "MaterialVolumeFraction")
      
      prognostic = have_option(trim(volumefraction%option_path)//"/prognostic")
      
      if(prognostic.and.(.not.aliased(volumefraction))) then
        call get_option(trim(volumefraction%option_path)//&
                        "/prognostic/surface_tension/surface_tension_coefficient", &
                        coeff, default=0.0)
      
        call zero(gradient)
        
        call calculate_div_t_cv(state(i), gradient)
        
        ! the magnitude of the field gradient is a regularisation of the delta function
        ! indicating where the interface is
        grad_mag = magnitude(gradient)
        
        ! normalise the gradient
        do node = 1, node_count(surfacetension)
          if(node_val(grad_mag, node)>epsilon(0.0)) then
            call set(gradient, node, node_val(gradient, node)/node_val(grad_mag, node))
          else
            call set(gradient, node, spread(0.0, 1, gradient%dim))
            call set(grad_mag, node, 0.0)
          end if
        end do
        
        ! if we have an equilibrium contact angle then modify the gradient near the requested walls
        ! (note that grad_mag is unchanged)
        call get_option(trim(volumefraction%option_path)//&
                       "/prognostic/surface_tension/equilibrium_contact_angle", eq_angle, stat)
        if(stat==0) then
          call allocate(normals, mesh_dim(surfacetension), surfacetension%mesh, "NormalsToBoundary")
          call zero(normals)
          
          allocate(on_boundary(node_count(surfacetension)))
          on_boundary = .false.
        
          shape_option=option_shape(trim(volumefraction%option_path) // &
                 & "/prognostic/surface_tension/equilibrium_contact_angle/surface_ids")           
          allocate(surface_ids(1:shape_option(1)))
          call get_option(trim(volumefraction%option_path)//&
                 &"/prognostic/surface_tension/equilibrium_contact_angle/surface_ids", surface_ids)

          call calculate_boundary_normals(surfacetension%mesh, x, &
                                          normals, on_boundary, &
                                          surface_ids = surface_ids)
          
          do node = 1, node_count(surfacetension)
            if(on_boundary(node)) then
              call set(gradient, node, &
                      (node_val(normals, node)*cos(eq_angle)+node_val(gradient, node)*sin(eq_angle)))
            end if
          end do
          
          grad_mag2 = magnitude(gradient)
        
          ! renormalise the gradient
          do node = 1, node_count(surfacetension)
            if(node_val(grad_mag2, node)>epsilon(0.0)) then
              call set(gradient, node, node_val(gradient, node)/node_val(grad_mag2, node))
            else
              call set(gradient, node, spread(0.0, 1, gradient%dim))
            end if
          end do
          
          call deallocate(grad_mag2)
          
          deallocate(on_boundary)
          call deallocate(normals)
          deallocate(surface_ids)
          
        end if
                
        do node = 1, node_count(surfacetension)
          tensor = 0.0
          do dimi = 1, size(tensor,1)
            do dimj = 1, size(tensor,2)
              if(dimi==dimj) tensor(dimi,dimj) = coeff*node_val(grad_mag, node)
              tensor(dimi,dimj) = tensor(dimi,dimj) - &
                                  coeff*node_val(gradient, dimi, node)*node_val(gradient, dimj, node)*&
                                  node_val(grad_mag, node)
            end do
          end do
          
          call addto(surfacetension, node, tensor)
          
        end do
        
        call deallocate(grad_mag)
        
      end if
    
    end do
    
    call deallocate(gradient)
    
  end subroutine calculate_surfacetension

  subroutine initialise_diagnostic_material_properties(state)

    type(state_type), dimension(:), intent(inout) :: state

    !locals  
    integer :: stat, i
    type(scalar_field), pointer :: sfield
    logical :: prognostic

    do i = 1, size(state)

      sfield=>extract_scalar_field(state(i),'MaterialDensity',stat)
      if(stat==0) then
        prognostic=(have_option(trim(sfield%option_path)//'/prognostic'))
        if((.not.aliased(sfield)).and. prognostic) then
          call compressible_material_eos(state(i),materialdensity=sfield)
        end if
      end if

    end do

  end subroutine initialise_diagnostic_material_properties
  
  subroutine calculate_diagnostic_material_volume_fraction(state)

    type(state_type), dimension(:), intent(inout) :: state
    
    !locals
    type(scalar_field), pointer :: materialvolumefraction
    integer :: i, stat, diagnostic_count, diagnostic_state_index
    type(scalar_field) :: sumvolumefractions
    type(scalar_field), pointer :: sfield
    logical :: diagnostic

    ! How many diagnostic MaterialVolumeFraction fields do we have in state?
    ! Note that state contains all the submaterials of the current phase, including the phase itself.
    ! Therefore, if the only material is the phase itself, diagnostic_count should be 0. Otherwise,
    ! it should be 1.
    diagnostic_count = 0
    do i = 1, size(state)
       if(have_option(trim(state(i)%option_path)//"/scalar_field::MaterialVolumeFraction/diagnostic")) then
          diagnostic_count = diagnostic_count + 1
          ! Record the index of the state containing the diagnostic MaterialVolumeFraction field
          diagnostic_state_index = i 
       end if
    end do

    if(diagnostic_count>1) then
      ewrite(0,*) diagnostic_count, ' diagnostic MaterialVolumeFractions'
      FLExit("Only one diagnostic MaterialVolumeFraction permitted.")
    end if

    if(diagnostic_count==1) then
      ! Extract the diagnostic volume fraction
      materialvolumefraction => extract_scalar_field(state(diagnostic_state_index), 'MaterialVolumeFraction')
      
      call allocate(sumvolumefractions, materialvolumefraction%mesh, 'Sum of volume fractions')
      call zero(sumvolumefractions)
      
      do i = 1,size(state)
        sfield=>extract_scalar_field(state(i),'MaterialVolumeFraction',stat)
        diagnostic=(have_option(trim(sfield%option_path)//'/diagnostic'))
        if ( (stat==0).and.(.not. aliased(sfield)).and.(.not.diagnostic)) then
          call addto(sumvolumefractions, sfield)
        end if
      end do
      
      call set(materialvolumefraction, 1.0)
      call addto(materialvolumefraction, sumvolumefractions, -1.0)
      call deallocate(sumvolumefractions)
    end if

  end subroutine calculate_diagnostic_material_volume_fraction
  
  subroutine order_states_priority(state, state_order)
    type(state_type), dimension(:), intent(inout) :: state
    integer, dimension(:), intent(inout) :: state_order
    
    type(scalar_field), pointer :: volumefraction
    logical, dimension(size(state)) :: priority_states
    integer, dimension(size(state)) :: state_priorities

    integer :: i, p, f

    assert(size(state_order)==size(state))

    priority_states = .false.
    state_priorities = 0
    do i = 1, size(state)
      volumefraction=>extract_scalar_field(state(i), "MaterialVolumeFraction")

      if(have_option(trim(volumefraction%option_path)//"/prognostic/priority")) then
        call get_option(trim(volumefraction%option_path)//"/prognostic/priority", state_priorities(i))
        priority_states(i) = .true.
      end if
    end do

    do i = 1, size(state)
      if(.not.priority_states(i)) then
        state_priorities(i) = minval(state_priorities)-1
      end if
    end do

    ! now work out the right order
    f = 0
    state_order = 0
    do p = maxval(state_priorities), minval(state_priorities), -1
      do i=1, size(state)
        if(state_priorities(i)==p) then
          f = f + 1
          state_order(f) = i
        end if
      end do
    end do
    assert(all(state_order>0))

  end subroutine order_states_priority

  subroutine calculate_bulk_scalar_property(state,bulkfield,materialname,mean_type,momentum_diagnostic)

    type(state_type), dimension(:), intent(inout) :: state
    type(scalar_field), intent(inout) :: bulkfield
    character(len=*), intent(in) :: materialname
    character(len=*), intent(in), optional :: mean_type
    logical, intent(in), optional ::momentum_diagnostic

    !locals
    integer :: i, stat
    character(len=OPTION_PATH_LEN) :: l_mean_type
    type(scalar_field), pointer :: sfield
 
    integer, dimension(size(state)) :: state_order
    type(scalar_field) :: sumvolumefractionsbound

    ewrite(1,*) 'In calculate_bulk_scalar_property:', trim(bulkfield%name)

    if (present(mean_type)) then
      l_mean_type = mean_type
    else
      l_mean_type = "arithmetic"
    end if

    select case(l_mean_type)
      case("arithmetic")
        call zero(bulkfield)
      case("harmonic")
        call zero(bulkfield)
      case("geometric")
        call set(bulkfield, 1.0)
      case default
        FLExit("Invalid mean_type in calculate_bulk_property")
    end select

    call order_states_priority(state, state_order)

    call allocate(sumvolumefractionsbound, bulkfield%mesh, "SumMaterialVolumeFractionsBound")
    call set(sumvolumefractionsbound, 1.0)
    
    do i = 1, size(state)
  
      ewrite(2,*) 'Considering state: ', state(state_order(i))%name
      sfield => extract_scalar_field(state(state_order(i)), trim(materialname), stat)
      if(stat==0) then
        call add_scaled_material_property(state(state_order(i)), bulkfield, sfield, &
                                          sumvolumefractionsbound=sumvolumefractionsbound, &
                                          mean_type=l_mean_type, momentum_diagnostic=momentum_diagnostic)
      end if

    end do

    select case(l_mean_type)
      case("harmonic")
        call invert(bulkfield, tolerance=tiny(0.0))
    end select
    
    call deallocate(sumvolumefractionsbound)
    
  end subroutine calculate_bulk_scalar_property

  subroutine calculate_bulk_vector_property(state,bulkfield,materialname,mean_type,momentum_diagnostic)

    type(state_type), dimension(:), intent(inout) :: state
    type(vector_field), intent(inout) :: bulkfield
    character(len=*), intent(in) :: materialname
    character(len=*), intent(in), optional :: mean_type
    logical, intent(in), optional :: momentum_diagnostic

    !locals
    integer :: i, stat
    character(len=OPTION_PATH_LEN) :: l_mean_type
    type(vector_field), pointer :: vfield
 
    integer, dimension(size(state)) :: state_order
    type(scalar_field) :: sumvolumefractionsbound

    ewrite(1,*) 'In calculate_bulk_vector_property:', trim(bulkfield%name)

    if (present(mean_type)) then
      l_mean_type = mean_type
    else
      l_mean_type = "arithmetic"
    end if

    select case(l_mean_type)
      case("arithmetic")
        call zero(bulkfield)
      case("harmonic")
        call zero(bulkfield)
      case("geometric")
        do i = 1, bulkfield%dim
          call set(bulkfield, i, 1.0)
        end do
      case default
        FLExit("Invalid mean_type in calculate_bulk_property")
    end select

    call order_states_priority(state, state_order)
    
    call allocate(sumvolumefractionsbound, bulkfield%mesh, "SumMaterialVolumeFractionsBound")
    call set(sumvolumefractionsbound, 1.0)

    do i = 1, size(state)
  
      ewrite(2,*) 'Considering state: ', state(state_order(i))%name
      vfield => extract_vector_field(state(state_order(i)), trim(materialname), stat)
      if(stat==0) then
        call add_scaled_material_property(state(state_order(i)), bulkfield, vfield, &
                                          sumvolumefractionsbound=sumvolumefractionsbound, &
                                          mean_type=l_mean_type, momentum_diagnostic=momentum_diagnostic)
      end if

    end do
    
    select case(l_mean_type)
      case("harmonic")
        call invert(bulkfield, tolerance=tiny(0.0))
    end select
    
    call deallocate(sumvolumefractionsbound)

  end subroutine calculate_bulk_vector_property

  subroutine calculate_bulk_tensor_property(state,bulkfield,materialname,mean_type,momentum_diagnostic)

    type(state_type), dimension(:), intent(inout) :: state
    type(tensor_field), intent(inout) :: bulkfield
    character(len=*), intent(in) :: materialname
    character(len=*), intent(in), optional :: mean_type
    logical, intent(in), optional :: momentum_diagnostic

    !locals
    integer :: i, j, stat
    character(len=OPTION_PATH_LEN) :: l_mean_type
    type(tensor_field), pointer :: tfield
 
    integer, dimension(size(state)) :: state_order
    type(scalar_field) :: sumvolumefractionsbound

    ewrite(1,*) 'In calculate_bulk_tensor_property:', trim(bulkfield%name)

    if (present(mean_type)) then
      l_mean_type = mean_type
    else
      l_mean_type = "arithmetic"
    end if

    select case(l_mean_type)
      case("arithmetic")
        call zero(bulkfield)
      case("harmonic")
        call zero(bulkfield)
      case("geometric")
        do i = 1, bulkfield%dim(1)
          do j = 1, bulkfield%dim(2)
            call set(bulkfield, i, j, 1.0)
          end do
        end do
      case default
        FLExit("Invalid mean_type in calculate_bulk_property")
    end select

    call order_states_priority(state, state_order)
    
    call allocate(sumvolumefractionsbound, bulkfield%mesh, "SumMaterialVolumeFractionsBound")
    call set(sumvolumefractionsbound, 1.0)

    do i = 1, size(state)
  
      ewrite(2,*) 'Considering state: ', state(state_order(i))%name
      tfield => extract_tensor_field(state(state_order(i)), trim(materialname), stat)
      if(stat==0) then
        call add_scaled_material_property(state(state_order(i)), bulkfield, tfield, &
                                          sumvolumefractionsbound=sumvolumefractionsbound, &
                                          mean_type=l_mean_type, momentum_diagnostic=momentum_diagnostic)
      end if
    
    end do

    select case(l_mean_type)
      case("harmonic")
        call invert(bulkfield, tolerance=tiny(0.0))
    end select
    
    call deallocate(sumvolumefractionsbound)
    
  end subroutine calculate_bulk_tensor_property

  subroutine get_scalable_volume_fraction(scaledvfrac, state, sumvolumefractionsbound, momentum_diagnostic)

    type(scalar_field), intent(inout) :: scaledvfrac
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout), optional :: sumvolumefractionsbound
    logical, intent(in), optional :: momentum_diagnostic

    type(scalar_field), pointer :: volumefraction, oldvolumefraction
    type(vector_field), pointer :: velocity
    type(scalar_field) :: remapvfrac

    integer :: stat
    real :: theta
    
    logical :: cap
    real:: u_cap_val, l_cap_val

    volumefraction => extract_scalar_field(state, 'MaterialVolumeFraction')
    
    call remap_field(volumefraction, scaledvfrac)
          
    if(present_and_true(momentum_diagnostic)) then
      velocity => extract_vector_field(state, 'Velocity', stat=stat)
      if(stat==0) then
        call get_option(trim(velocity%option_path)//'/prognostic/temporal_discretisation/theta', &
                        theta, stat)
        if(stat==0) then
          call allocate(remapvfrac, scaledvfrac%mesh, "RemappedMaterialVolumeFraction")
          
          oldvolumefraction => extract_scalar_field(state, 'OldMaterialVolumeFraction')
          call remap_field(oldvolumefraction, remapvfrac)
          
          call scale(scaledvfrac, theta)
          call addto(scaledvfrac, remapvfrac, (1.-theta))
          
          call deallocate(remapvfrac)
        end if
      end if
    end if
          
    cap = (have_option(trim(complete_field_path(volumefraction%option_path))//"/cap_values"))
          
    if(cap) then
      ! this capping takes care of under or overshoots in this volume fraction individually
      ! these will have typically occurred during advection
            
      call get_option(trim(complete_field_path(volumefraction%option_path))//"/cap_values/upper_cap", &
                      u_cap_val, default=huge(0.0)*epsilon(0.0))
      call get_option(trim(complete_field_path(volumefraction%option_path))//"/cap_values/lower_cap", &
                      l_cap_val, default=-huge(0.0)*epsilon(0.0))
            
      call bound(scaledvfrac, l_cap_val, u_cap_val)
          
      if(present(sumvolumefractionsbound)) then
        assert(sumvolumefractionsbound%mesh==scaledvfrac%mesh)
        ! this capping takes care of overlapping volume fractions
        call bound(scaledvfrac, upper_bound=sumvolumefractionsbound) 
        call addto(sumvolumefractionsbound, scaledvfrac, scale=-1.0)
        ewrite_minmax(sumvolumefractionsbound)
      end if
          
    end if
    ewrite_minmax(scaledvfrac)
 
  end subroutine get_scalable_volume_fraction

  subroutine add_scaled_material_property_scalar(state,bulkfield,field,sumvolumefractionsbound,mean_type,momentum_diagnostic)

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: bulkfield, field
    type(scalar_field), intent(inout), optional :: sumvolumefractionsbound
    logical, intent(in), optional :: momentum_diagnostic
    character(len=*), intent(in), optional :: mean_type

    !locals
    character(len=OPTION_PATH_LEN) :: l_mean_type
    type(scalar_field) :: scaledvfrac
    type(scalar_field) :: tempfield
    
    if (present(mean_type)) then
      l_mean_type = mean_type
    else
      l_mean_type = "arithmetic"
    end if

    call allocate(tempfield, bulkfield%mesh, "Temp"//trim(bulkfield%name))
    call allocate(scaledvfrac, bulkfield%mesh, "ScaledMaterialVolumeFraction")

    call get_scalable_volume_fraction(scaledvfrac, state, &
                                      sumvolumefractionsbound=sumvolumefractionsbound, &
                                      momentum_diagnostic=momentum_diagnostic)

    call remap_field(field, tempfield)
    select case(l_mean_type)
      case("arithmetic")
        call scale(tempfield, scaledvfrac)
        call addto(bulkfield, tempfield)
      case("harmonic")
        call invert(tempfield, tolerance=tiny(0.0))
        call scale(tempfield, scaledvfrac)
        call addto(bulkfield, tempfield)
      case("geometric")
        call power(tempfield, scaledvfrac)
        call scale(bulkfield, tempfield)
      case default
        FLExit("Invalid mean_type in add_scaled_material_property")
    end select

    call deallocate(tempfield)
    call deallocate(scaledvfrac)

  end subroutine add_scaled_material_property_scalar

  subroutine add_scaled_material_property_vector(state,bulkfield,field,sumvolumefractionsbound,mean_type,momentum_diagnostic)

    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: bulkfield, field
    type(scalar_field), intent(inout), optional :: sumvolumefractionsbound
    logical, intent(in), optional :: momentum_diagnostic
    character(len=*), intent(in), optional :: mean_type

    !locals
    character(len=OPTION_PATH_LEN) :: l_mean_type
    type(scalar_field) :: scaledvfrac
    type(vector_field) :: tempfield

    if (present(mean_type)) then
      l_mean_type = mean_type
    else
      l_mean_type = "arithmetic"
    end if

    call allocate(tempfield, bulkfield%dim, bulkfield%mesh, "Temp"//trim(bulkfield%name))
    call allocate(scaledvfrac, bulkfield%mesh, "ScaledMaterialVolumeFraction")
    
    call get_scalable_volume_fraction(scaledvfrac, state, &
                                      sumvolumefractionsbound=sumvolumefractionsbound, &
                                      momentum_diagnostic=momentum_diagnostic)
          
    call remap_field(field, tempfield)
    select case(l_mean_type)
      case("arithmetic")
        call scale(tempfield, scaledvfrac)
        call addto(bulkfield, tempfield)
      case("harmonic")
        call invert(tempfield, tolerance=tiny(0.0))
        call scale(tempfield, scaledvfrac)
        call addto(bulkfield, tempfield)
      case("geometric")
        call power(tempfield, scaledvfrac)
        call scale(bulkfield, tempfield)
      case default
        FLExit("Invalid mean_type in add_scaled_material_property")
    end select
          
    call deallocate(tempfield)
    call deallocate(scaledvfrac)

  end subroutine add_scaled_material_property_vector

  subroutine add_scaled_material_property_tensor(state,bulkfield,field,sumvolumefractionsbound,mean_type,momentum_diagnostic)

    type(state_type), intent(inout) :: state
    type(tensor_field), intent(inout) :: bulkfield, field
    type(scalar_field), intent(inout), optional :: sumvolumefractionsbound
    logical, intent(in), optional :: momentum_diagnostic
    character(len=*), intent(in), optional :: mean_type

    !locals
    character(len=OPTION_PATH_LEN) :: l_mean_type
    type(scalar_field) :: scaledvfrac
    type(tensor_field) :: tempfield

    if (present(mean_type)) then
      l_mean_type = mean_type
    else
      l_mean_type = "arithmetic"
    end if

    call allocate(tempfield, bulkfield%mesh, "Temp"//trim(bulkfield%name))
    call allocate(scaledvfrac, bulkfield%mesh, "ScaledMaterialVolumeFraction")

    call get_scalable_volume_fraction(scaledvfrac, state, &
                                      sumvolumefractionsbound=sumvolumefractionsbound, &
                                      momentum_diagnostic=momentum_diagnostic)
    
    call remap_field(field, tempfield)
    select case(l_mean_type)
      case("arithmetic")
        call scale(tempfield, scaledvfrac)
        call addto(bulkfield, tempfield)
      case("harmonic")
        call invert(tempfield, tolerance=tiny(0.0))
        call scale(tempfield, scaledvfrac)
        call addto(bulkfield, tempfield)
      case("geometric")
        call power(tempfield, scaledvfrac)
        call scale(bulkfield, tempfield)
      case default
        FLExit("Invalid mean_type in add_scaled_material_property")
    end select
          
    call deallocate(tempfield)
    call deallocate(scaledvfrac)

  end subroutine add_scaled_material_property_tensor

  subroutine calculate_material_volume(state, materialvolume)

    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: materialvolume
  
    ! local
    type(scalar_field) :: cvmass
    type(scalar_field), pointer :: volumefraction
    type(vector_field), pointer :: coordinates
  
    coordinates=>extract_vector_field(state, "Coordinate")
  
    call allocate(cvmass, materialvolume%mesh, "CV mass")
    call zero(cvmass)
  
    call compute_cv_mass(coordinates, cvmass)
  
    volumefraction=>extract_scalar_field(state,"MaterialVolumeFraction")
  
    materialvolume%val=volumefraction%val*cvmass%val
  
    call deallocate(cvmass)

  end subroutine calculate_material_volume

  subroutine calculate_material_mass(state, materialmass)

    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: materialmass
  
    ! local
    integer :: stat
    type(scalar_field) :: cvmass
    type(scalar_field), pointer :: volumefraction, materialdensity
    type(vector_field), pointer :: coordinates
    real :: rho_0
  
    coordinates=>extract_vector_field(state, "Coordinate")
  
    call allocate(cvmass, materialmass%mesh, "CV mass")
    call zero(cvmass)
  
    call compute_cv_mass(coordinates, cvmass)

    volumefraction=>extract_scalar_field(state,"MaterialVolumeFraction")
    call set(materialmass, volumefraction)
    call scale(materialmass, cvmass)
  
    materialdensity=>extract_scalar_field(state,"MaterialDensity", stat=stat)
    if(stat==0) then
      call scale(materialmass, materialdensity)
    else
      call get_option("/material_phase::"//trim(state%name)&
                      //"/equation_of_state/fluids/linear/reference_density", rho_0)
      call scale(materialmass, rho_0)
    end if
  
    call deallocate(cvmass)

  end subroutine calculate_material_mass

  subroutine calculate_bulk_material_pressure(state,bulkmaterialpressure)

    type(state_type), dimension(:), intent(inout) :: state
    type(scalar_field) :: bulkmaterialpressure

    !locals
    integer :: i, stat
    type(scalar_field), pointer :: volumefraction
    type(scalar_field) :: materialpressure

    call zero(bulkmaterialpressure)

    call allocate(materialpressure, bulkmaterialpressure%mesh, "TempBulkPressure")

    do i = 1, size(state)

      volumefraction => extract_scalar_field(state(i), 'MaterialVolumeFraction', stat)

      if (stat==0) then

         call compressible_material_eos(state(i), materialpressure=materialpressure)

         bulkmaterialpressure%val=bulkmaterialpressure%val+volumefraction%val*materialpressure%val

      end if

    end do

    call deallocate(materialpressure)

  end subroutine calculate_bulk_material_pressure
  
  subroutine calculate_sum_material_volume_fractions(state,sumvolumefractions)

    type(state_type), dimension(:), intent(inout) :: state
    type(scalar_field), intent(inout) :: sumvolumefractions

    !locals
    integer :: i, stat
    type(scalar_field), pointer :: sfield
    logical :: prognostic, diagnostic, prescribed

    diagnostic = have_option(trim(sumvolumefractions%option_path)//"/diagnostic")
    if(.not.diagnostic) return

    call zero(sumvolumefractions)

    do i = 1,size(state)
      sfield=>extract_scalar_field(state(i),'MaterialVolumeFraction',stat)
      if(stat==0) then
        prognostic = have_option(trim(sfield%option_path)//"/prognostic")
        prescribed = have_option(trim(sfield%option_path)//"/prescribed")
        if ((.not.aliased(sfield)).and.(prognostic.or.prescribed)) then
          call addto(sumvolumefractions, sfield)
        end if
      end if
    end do

  end subroutine calculate_sum_material_volume_fractions

end module multimaterial_module
