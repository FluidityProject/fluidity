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

module porous_media
  !! This module constructs terms for solving darcy's law:
  !! 1. It adapts the advection source to account for porous media
  !! 2. It puts the viscosity-permeability term into the momentum equation
  !! 
  !! It also deals with electrokinetic coupling and the additional terms
  !! required for multiphase flow in porous media (i.e. capillary pressure
  !! and relative permeabilities etc
  use fldebug
  use state_module
  use fields
  use spud
  implicit none

  public :: porous_media_advection, &
            porous_media_momentum, &
            capillary_pressure, &
            relative_permeability, &
            normalised_saturation, &
            brine_salinity, &
            calculate_porous_media_absorption

  contains
  
  subroutine porous_media_advection(state)
    ! specify interface variables
    type(state_type), dimension(:), intent(inout) :: state

    ! specify local variables
    integer :: nstates, i

    ewrite(3,*) 'In porous_media_advection'

    ! counting phases -- is it multiphase?
    nstates = option_count("/material_phase")

    do i=1,nstates
      ewrite(3,*) 'Considering phase',i,'of',nstates
    end do
  end subroutine porous_media_advection

  subroutine porous_media_momentum(state)
    ! declare interface variables
    type(state_type), dimension(:), intent(inout) :: state
    ! declare local variables
    integer :: nstates, stat, i, j
    real, pointer :: kr, kr_expo, cp_A, cp_expo
    real, dimension(:), pointer :: relperm
    type(scalar_field), pointer :: porosity, phase_volume_fraction
    type(scalar_field) :: viscosityxx
    type(vector_field), pointer :: velocity, absorption
    type(tensor_field), pointer :: viscosity

    ewrite(3,*) 'In porous_media_momentum'

    ! Solving equation: q = - k/eta grad P
    ! Calculate absorption term: volume_fraction * eta / k

    ! Counting phases - is multiphase?
    nstates = option_count("/material_phase")

    ! Get porosity and permeability out of state
    porosity => extract_scalar_field(state(1),"Porosity",stat)
    if (stat==0) then
      ewrite(3,*) 'Porosity present and correct'
      ewrite_minmax(porosity)
    else
      FLAbort("Problem reading in porosity")
    end if


    if (nstates>1) then
      ! MULTIPHASE
      ewrite(3,*) 'Multiphase'
      do i=1,size(state)
        ! read in volume fraction
        phase_volume_fraction => extract_scalar_field(state(i),"PhaseVolumeFraction",stat)

        ! call to capillary pressure subroutine
        if (have_option("/porous_media/multiphase_parameters/cp_A")) then
          !call capillary_pressure()
        end if

        ! add gravity
        if (have_option("/physical_parameters/gravity")) then
          !add some density related effect
          ! Is this necessary - or is it taken care of elsewhere?
        end if

        ! calculate absorption term for each phase
        ! read in viscosity
        viscosity => extract_tensor_field(state(i),"Viscosity",stat)
        if (stat==0) then
           ewrite(3,*) 'Found a viscosity field'
        else
           FLAbort('Need viscosity field')
        end if
        viscosityxx=extract_scalar_field_from_tensor_field(viscosity, 1, 1)
        ewrite_minmax(viscosityxx)
        ! read in rel perm options
        call get_option("/porous_media/multiphase_parameters/kr" // int2str(i) //"", kr)
        call get_option("/porous_media/multiphase_parameters/kr_expo" //int2str(i) //"", kr_expo)
        ! end if
        !       call allocate(relperm, phase_volume_fraction%mesh, name="Relperms")
        call relative_permeability(phase_volume_fraction,relperm,kr,kr_expo)
      end do
    else
      ! SINGLE PHASE
      ewrite(3,*) 'Single phase'

      ! check absorption term exists
      velocity => extract_vector_field(state(1),"Velocity",stat)
      absorption => extract_vector_field(state(1),"VelocityAbsorption",stat)
      if ((stat/=0).and.have_option(trim(velocity%option_path)//"/prognostic")) then
        FLAbort('Need to switch on velocity absorption')
      end if
    end if
  end subroutine porous_media_momentum

  subroutine capillary_pressure()
    !!< compute capillary pressure
  end subroutine capillary_pressure

  subroutine relative_permeability(volume_fraction,relperm,kr,expo)
    !!< calculating parameterisation of rel perms:
    !!< relperm = kr
    ! declare interface variables
    type(scalar_field) :: volume_fraction
    real, dimension(:), pointer :: relperm
    real :: kr, expo
  end subroutine relative_permeability

  function normalised_saturation(Sw,Swc,Sro)
    !!< calculate normalised saturation, given
    !!< Sw  = water saturation
    !!< Swc = connate water saturation
    !!< Sro = residual oil saturation
    real, intent(in) :: Sw, Swc, Sro
    real :: normalised_saturation, nominator
    nominator = Sw-Swc
    if (abs(nominator)>1.0e-3) then
      normalised_saturation = nominator/(1.0-Swc-Sro)
    else
      normalised_saturation = 0.0
    end if
  end function normalised_saturation

  function brine_salinity(sigma_w)
    !!< calculate brine salinity from specified electrical conductivity
    !!< see Worthington et al. (1990)
    real, intent(in) :: sigma_w
    real :: brine_salinity, logsw
    real, dimension(5) :: c
    c(1) = -1.03024
    c(2) = 1.06627
    c(3) = 2.41239e-2
    c(4) = 3.68102e-3
    c(5) = 1.46369e-4
    logsw = log10(sigma_w)
    brine_salinity = 10**(c(1)+c(2)*logsw+c(3)*logsw**2+c(4)*logsw**3+c(5)*logsw**4)
  end function brine_salinity

  subroutine calculate_porous_media_absorption(state, i, absorption, stat)
    !!< calculate the absorption term for a vector velocity absorption field
    ! specify interface variables
    type(state_type), dimension(:), intent(inout) :: state
    integer, intent(in) :: i
    type(vector_field), intent(inout) :: absorption
    integer, intent(out) :: stat

    ! specify local variables
    integer :: j, k
    type(scalar_field) :: s_viscosity, s_permeability
    type(vector_field), pointer :: v_permeability
    type(tensor_field), pointer :: t_viscosity, t_permeability

    stat = -1

    ! extract viscosity field as a scalar field
    t_viscosity => extract_tensor_field(state(i), "Viscosity", stat)
    if (stat/=0) then
       FLAbort('Phase '//int2str(i)//' is missing a viscosity field.')
    end if
    s_viscosity = extract_scalar_field_from_tensor_field(t_viscosity, 1, 1)
    ewrite_minmax(s_viscosity)

    ! check permeability field and appropriately compute absorption term
    if (have_option("/porous_media/scalar_field::Permeability")) then
      ! SCALAR PERMEABILITY
      s_permeability = extract_scalar_field(state(1), "Permeability", stat)
      do j=1,absorption%dim
        do k=1,node_count(absorption)
          call set(absorption, j, k, s_viscosity%val(k)/s_permeability%val(k))
        end do
      end do
    elseif (have_option("/porous_media/vector_field::Permeability")) then
      ! VECTOR PERMEABILITY
      v_permeability => extract_vector_field(state(1), "Permeability", stat)
      do j=1,absorption%dim
        s_permeability = extract_scalar_field_from_vector_field(v_permeability, j)
        do k=1,node_count(absorption)
          call set(absorption, j, k, s_viscosity%val(k)/s_permeability%val(k))
        end do
      end do
    elseif (have_option("/porous_media/tensor_field::Permeability")) then
      ! TENSOR PERMEABILITY
      t_permeability => extract_tensor_field(state(1), "Permeability", stat)
      FLAbort("Cannot currently solve porous media problems with a full tensor permeability.")
    end if
    if (stat/=0) then
      FLAbort("Porous media problems require a permeability field.")
    end if

    stat = 0
  end subroutine calculate_porous_media_absorption

end module porous_media

