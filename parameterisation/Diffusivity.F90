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

module diffusivity
  !!< A wrapper module which calls routines implementing different
  !!< diffusivity models. These routines are invoked by setting 
  !!< 50<DISOTT<60.
  !!<
  !!< The scheme to be used and the parameters are specified in the file
  !!< diffusivity_input.dat which must be placed in the same directory as the
  !!< fluidity input files. diffusivity_input.dat is a Fortran namelist input
  !!< file with the following format:
  !!<
  !!<  &diffusivity_input
  !!<  scheme = "scheme_name",
  !!<  <variable> = <value>,
  !!<  ...
  !!<  <variable> = <value>,
  !!<  /
  !!<
  !!< The following variables apply to the following schemes:
  !!<
  !!< ==Isotropic Diffusivity==
  !!<  scheme="Isotropic"
  !!<  real :: iso_diffusivity ! The diffusivity
  !!<
  !!< ==Anisotropic Diffusivity==
  !!<  scheme="Anisotropic"
  !!<  real :: horizontal_diffusivity, vertical_diffusivity
  !!<
  !!< ==Gent McWilliams==
  !!<  scheme="Gent McWilliams"
  !!<  real :: a_isoneutral ! The isoneutral diffusivity.
  !!<  real :: kappa        ! The Gent McWilliams diffusivity
  !!<  real :: R1           ! 1st baroclinic Rossby radius (R1=c/f)
  !!<  real :: S_max        ! maximum slope
  !!<  real :: S_d          ! width of transition region for exponential taper
  use Isotropic
  use FLDebug
  use Anisotropic
  use Gent_McWilliams
  use Elements
  use FETools
  use spud
  use global_parameters, only: new_options, OPTION_PATH_LEN, FIELD_NAME_LEN
  implicit none

  private

  public :: initialise_diffusivity, get_diffusivity, diffusivity_check_options

  !! a_i is isoneutral diffusion and kappa is the Gent McWilliams
  !! diffusivity.
  !! R1 is the first baroclinic Rossby radius
  !! S_max is the (user defined) maximum density slope above which tapering
  !! is required
  !! S_d controls the width of the exponential taper
  real, save :: a_isoneutral, kappa, R1, S_max, S_d, c, iso_diffusivity, & 
       & horizontal_diffusivity, vertical_diffusivity


  !! Diffusivity scheme in use.
  character(len=42), public, save :: diffusivity_scheme    

contains

  subroutine initialise_diffusivity
    !!< The namelist must list all the parameters it is possible to read
    !!< from diffusivity_input.dat

    ! Unit number to read diffusivity information on.
    integer :: unit

    ! check that the config file is read only once.
    logical, save :: initialised=.false.

    ! diffusivity scheme.
    character(len=42), save :: scheme    

    character(len=OPTION_PATH_LEN) :: GM_path

    namelist/diffusivity_input/scheme, a_isoneutral, kappa, R1, S_max, S_d, &
         & iso_diffusivity, horizontal_diffusivity, vertical_diffusivity

    if (.not. initialised) then

       if(new_options) then

          GM_path="/material_phase[0]/subgridscale_parameterisations/Gent_McWilliams"
          if(have_option(trim(GM_path))) then

             scheme="Gent McWilliams"
             call get_option(trim(GM_path)//"/isoneutral_diffusivity", a_isoneutral)
             call get_option(trim(GM_path)//"/GentMcWilliams_diffusivity", kappa)
             call get_option(trim(GM_path)//"/tapering/maximum_density_slope", S_max)
             call get_option(trim(GM_path)//"/tapering/first_baroclinic_Rossby_radius", R1)
             call get_option(trim(GM_path)//"/tapering/width_of_exponential_taper_region", S_d)
             call get_option(trim(GM_path)//"/tapering/constant_for_linear_taper", c)
          end if

       else

          ewrite(-1,*) 'Calling GM without new options is a bit broken'
          FLAbort('Calling GM without new options is a bit broken')
          unit=free_unit()
          open(unit=unit, file="diffusivity_input.dat", action="read",&
               & status="old", err=666)
          
          read(unit, nml=diffusivity_input)
       
          close(unit)

       end if

       ewrite(2, *) "Diffusivity scheme:", scheme
       initialised=.true.

       diffusivity_scheme=scheme

    end if

    return
       
666 FLAbort("Required file diffusivity_input.dat could not be opened")

  end subroutine initialise_diffusivity

  function get_diffusivity(tracer_shape, density, density_shape, &
       & density_dshape, topdis, botdis, d, gi) result (mu)
    !!< Calculates diffusivity (mu) for the field tracer at Gauss point gi

    !! Shape functions for input fields.
    type(element_type), intent(in):: tracer_shape, density_shape
    !! Value of input fields at nodes.
    real, intent(in):: density(density_shape%loc)
    !! Derivatives of input fields at nodes.
    real, intent(in):: &
         density_dshape(density_shape%loc, density_shape%ngi, density_shape%dim) 
    !! Distance from top/bottom
    real, intent(in):: topdis(density_shape%ngi),botdis(density_shape%ngi)
    !! some measure (don't ask) of distance from boundaries
    real, intent(in):: d(density_shape%ngi)
    !! Gauss point at which to calculate diffusivity
    integer, intent(in):: gi

    !! diffusivity tensor (3x3 in 3 dimensions)
    real :: mu(tracer_shape%dim, tracer_shape%dim)
       
    call initialise_diffusivity
 
    select case (diffusivity_scheme)

    case ("Isotropic")

       call isotropic_diffusivity(mu, iso_diffusivity)
       
    case ("Anisotropic")

       call anisotropic_diffusivity(mu, horizontal_diffusivity, &
            & vertical_diffusivity)
       
    case ("Gent McWilliams")

       call gent_mcwilliams_diffusivity(mu, density, density_shape, &
            & density_dshape, topdis, botdis, d, gi, a_isoneutral, &
            & kappa, R1, S_max, S_d, c)

    case default
       ! Can't happen
       FLAbort("Unknown diffusivity scheme")
       
    end select
        
  end function get_diffusivity

  subroutine diffusivity_check_options
    !! This check ensures that the subgridscale parameterisations requested
    !! for particular scalar fields have actually been switched on in Diamond.
    character(len=OPTION_PATH_LEN) :: field_path, phase_path
    character(len=FIELD_NAME_LEN) :: field_name, phase_name,&
         & parameterisation

    integer :: i, j, nstates, nfields

    nstates=option_count("/material_phase")
    
    state_loop: do i=0, nstates-1
       
       phase_path="/material_phase["//int2str(i)//"]"
       
       ! Get number of scalar fields that are children of this state
       nfields=option_count(trim(phase_path)//"/scalar_field")
       
       ! Loop over scalar fields
       scalar_field_loop: do j=0, nfields-1
          
          ! Save path to field
          field_path=trim(phase_path)//"/scalar_field["//int2str(j)//"]"

          if (have_option(trim(field_path)//&
               "/prognostic/subgridscale_parameterisation")) then
             
             call get_option(trim(field_path)//&
               "/prognostic/subgridscale_parameterisation/name", parameterisation)

             ! Ensure that the requested diffusivity scheme is available.
             if (.not.have_option(trim(phase_path)//&
                  "/subgridscale_parameterisations/"//trim(parameterisation)))&
                  & then

                call get_option(trim(phase_path)//"/name", phase_name)
                call get_option(trim(field_path)//"/name", field_name)
                
                ewrite(0,*) trim(phase_name)//"::"//trim(field_name)//&
                     " wants "//trim(parameterisation)//&
                     " but it is not enabled."
                FLExit("Error: inconsistent parameterisation options")

             end if

             ! Check for users specifying a parameterisation and
             ! prescribing a diffusivity.
             if (have_option(trim(field_path)//&
               "/prognostic/tensor_field::Diffusivity")) then

                call get_option(trim(phase_path)//"/name", phase_name)
                call get_option(trim(field_path)//"/name", field_name)

                ewrite(0,*) trim(phase_name)//"::"//trim(field_name)//&
                     " has both a Diffusivity and a subgridscale &
                     &parameterisation."
                ewrite(0,*) "This is not currently supported."
                FLExit("Error: inconsistent parameterisation options")

             end if

          end if

       end do scalar_field_loop

    end do state_loop

  end subroutine diffusivity_check_options

end module diffusivity
