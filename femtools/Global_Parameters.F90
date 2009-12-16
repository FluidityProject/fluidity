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

module global_parameters
  !!< This routine exists to save us all from argument list hell!
  !!<
  !!< All the global parameters which don't change while fluidity is running
  !!< should live here. I am building this up as I encounter more parameters
  !!< in the code. It would be great if others did the same.
  !!<
  !!< The correct syntax for accessing this module is:
  !!<
  !!< use global_parameters, only: parameter1, parameter2 ...
  !!<
  !!< Try to only use the parameters which are needed locally.
  
  ! Debug specific paramaters are contained in fldebug_parameters
  ! (to resolve build dependencies)
  use fldebug_parameters
  
  implicit none
    
  !------------------------------------------------------------------------
  ! As yet unidentified parameters.
  !------------------------------------------------------------------------
  
  !! Scaling factors?
  ! This is still used by fluids, so we'd better make sure it's initialised
  real, save :: SCFACTH0 = 1.0
    
  !------------------------------------------------------------------------
  ! Global flag for whether we are running off the new xml options file.
  ! Legacy variable - to be removed in the future
  !------------------------------------------------------------------------    
  logical, parameter :: new_options = .true.

  !------------------------------------------------------------------------
  ! Precision parameters
  !------------------------------------------------------------------------
  !! Number of digits past the decimal point for a real
  integer, parameter :: real_digits_10 = precision(0.0)

  !! Integer size in bytes
  integer, parameter :: integer_size=bit_size(0)/8
  ! The real_size depends on reals having the obvious ieee754 sizes. 
  !! Real size in bytes
#ifdef DOUBLEP
  integer, parameter :: real_size=8
#else
  integer, parameter :: real_size=4
#endif

  integer, parameter :: int_4 = selected_int_kind(4), &
                      & int_8 = selected_int_kind(8), &
                     & int_16 = selected_int_kind(16)
  ! real variable declarations of the form:
  !   real*4 :: real_var
  ! are not portable. Use these instead:
  !   real(real_4) :: real_val
  integer, parameter :: real_4 = selected_real_kind(4), &
                      & real_8 = selected_real_kind(8), &
                     & real_16 = selected_real_kind(16)

  !------------------------------------------------------------------------
  ! Parameters specifying scope bounds for fluidity.
  !------------------------------------------------------------------------    
  !! Maximum number of phases supported by fluidity. 
  integer, parameter :: MXNPHA=4

  !------------------------------------------------------------------------
  ! Parameters controlling the scheme used in the flow core.
  !------------------------------------------------------------------------
  
  !! The timestep.
  real, save, target ::  dt
  
  !! The simulation start time
  real, save :: simulation_start_time
  !! The simulation start CPU time
  real, save :: simulation_start_cpu_time
  !! The simulation start wall time
  real, save :: simulation_start_wall_time
  
  !! Accumulated system time.
  real, save, target :: ACCTIM

  !! The implicitness of the Crank-Nicolson scheme for the momentum
  !! equation.
  real, save, target :: THETA(MXNPHA)

  !! The discretisation of the momentum equaition is controlled through
  !! DISOPT(DISCRETIZATION OPTION).
  !!
  !! All methods use a balancing diffusion term and THETA time stepping. 
  !! The following are absolute values of DISOPT...
  !!  * DISOPT=1 - balancing diffusion based on (x,y) space.
  !!  * DISOPT=2 - Laxwendrof balancing diffusion.
  !!  * DISOPT=3 - (x,y,t) -balancing diffusion. 
  !!  * DISOPT=4 - No balancing diffusion.
  !!  * DISOPT=5 - nonlinear streamline and cross stream diffusion.
  !!  * DISOPT=6 - nonlinear upwind in steapest direction.
  !!  * DISOPT=7 - nonlinear streamline+ cross stream diffusion(but restricted)
  !!  * DISOPT=42- LES option using constant length scale.
  !!  * DISOPT=43- LES option using isotropic length scale.
  !!  * DISOPT=44- LES option which uses no balancing diffusion.
  !!  * DISOPT=45- LES option which uses no balancing diffusion.
  !!  * DISOPT=46- same as 45 but with 4th order dissipation.
  !!  * 70<=DISOPT<=85 - Discontinuous Galerkin schemes.
  !!  * DISOPT=125 - NO balancing diffusion(DISOPT=4)and take out non-linear terms.
  integer, save, target :: DISOPT(MXNPHA)

  !! Disopn is the discretisation option for nonlinear term.
  !!
  !! If NDISOP.ne.0 then set advection to zero and treat advection 
  !! using the high res method.
  integer, save, target :: ndisop(MXNPHA)
  !! NSUBVLOC is the no. of basis for sub-element modelling for vels.
  !! NSUBNVLOC is the no. of basis for sub-element modelling of advection.
  integer, save, target :: NSUBVLOC(MXNPHA)
  integer, save, target :: NSUBNVLOC(MXNPHA)

  !! If lump then lump the mass matrix in the momentum equation.
  logical, save, target :: LUMP(MXNPHA)  

  !! If abslum then lump the absorbtion.
  logical, save, target :: ABSLUM(MXNPHA)

  real, parameter:: pi = 3.1415926535897931

  !------------------------------------------------------------------------
  ! Parameters specifying the coordinate system.
  !------------------------------------------------------------------------

  !! Flag for 3 dimensions.
  logical, save :: D3
  !! Flag for cylindrical coordinates. 
  logical, save :: DCYL  

  !! Flag for pseudo2d domains.
  integer, save :: pseudo2d_coord = 0

  !! Flag for a cartesian sphere. 0 is a flat earth, other values are the
  !! different sphere options.
  integer, save :: isphere=0

  !------------------------------------------------------------------------
  ! Parameters controlling geophysical flow.
  !------------------------------------------------------------------------
  
  !! Parameter controlling the geostrophic balance pressure.
  integer, save :: GEOBAL
  real, save :: GRAVTY !gravity magnitude
  real, dimension(MXNPHA), save :: dengam
  
  !------------------------------------------------------------------------
  ! Parameters for parallel
  !------------------------------------------------------------------------
  
  !! Halo tags of the first and second halo
  integer, parameter :: halo_tag = 1, halo_tag_p = 2
  !! Halo tag for the quadratic tetrahedral meshes
  integer, parameter :: halo_tag_t10 = 3
  !! Halo tag for level 1 and 2 dg (ie face adjacency) halos
  !! Note that the -1 halo is currently broken.
  integer, parameter :: halo_tag_dg1 = -2, halo_tag_dg2 = -2
  
  !! When upscaling a problem (e.g. from 16 to 32 processors),
  !! we want to run sam on 32 processors to do the domain decomposition.
  !! But only 16 processors will have data on disk. However,
  !! all 32 processors still have to go through populate_state
  !! to make sure it goes through all the MPI calls and doesn't
  !! deadlock. So we record whether the process is an "active" process,
  !! one that has data on disk.
  logical :: is_active_process = .true.

  !------------------------------------------------------------------------
  ! Ident names and numbers.
  !------------------------------------------------------------------------

  ! Zeroing long strings is EXPENSIVE.
  ! (See the commit message for r11059)
  ! That is why we supply an empty_name and an empty_path
  ! as, e.g.,
  ! field%option_path=empty_path
  ! is much much quicker than
  ! field%option_path="" .
  ! This is probably a bug in gcc, but it is a bug in gcc
  ! that we have to live with.

  !! Field names are permitted to be as long as Fortran names.
  integer, parameter :: FIELD_NAME_LEN=101
  character(len=FIELD_NAME_LEN) :: empty_name=""

  !! Maximum length of an option path
  integer, parameter :: OPTION_PATH_LEN=8192
  character(len=OPTION_PATH_LEN) :: empty_path=""

  !! Maximum length of a python string representing a function
  integer, parameter :: PYTHON_FUNC_LEN=8192
  character(len=PYTHON_FUNC_LEN) :: empty_python_func=""

  !! indexes between states and phases
  integer, save, dimension(:), allocatable :: phase2state_index, state2phase_index

  type ident
     !!< Type for storing idents. These associate GEM ident numbers with
     !!< field names.
     character(len=FIELD_NAME_LEN) ::name
     integer :: ident
  end type ident

  !! list of legacy "idents" Don't add any new ones!!
  integer, parameter :: ident_count=11
  type(ident), dimension(ident_count), parameter :: ident_list = (/ &
       ident("Salinity", 42), &
       ident("Temperature", -1), &
       ident("Tracer", 666), & ! This is for totally passive tracer tests.
       ident("SecondFluid",56), &
       ident("DiffuseInterface",57), &
       ident("KineticEnergy", 101), &
       ident("TurbulentLengthScalexKineticEnergy", 102), &
       ident("VerticalViscosity", 103), &
       ident("VerticalDiffusivity", 104), &
       ident("FreeSurface", -29), &
       ident("BalancePressure", -2003) &
       /)


  private :: ident, ident_list

contains
  
  elemental function ident_name(ident)
    !!< Given an ident number return the field name.
    character(len=FIELD_NAME_LEN) :: ident_name
    integer, intent(in) :: ident

    integer :: loc

    loc=0
    ! Unfortunately the following line caused an internal compiler error in
    ! gfortran.
    !loc=minval(ident_list%ident, (ident_list%ident==ident))
    do loc=1,size(ident_list)
       if (ident_list(loc)%ident==ident) then
          goto 42
       end if
    end do
    
    ! Fallthrough case for unknown fields.
!     write(ident_name, '(a,i0)') "UnknownField", ident  ! can't use this as in case of
!                                                        ! negative idents it introduces hyphens
!                                                        ! which aren't allowed in the options dictionary
    write(ident_name, '(a,i0)') "UnknownField", abs(ident) ! very nonideal solution
  
    return
    
42  ident_name=ident_list(loc)%name
    
  end function ident_name

  pure function name_ident(name)
    !!< Given a field name return the ident number
    character(len=FIELD_NAME_LEN), intent(in) :: name
    integer :: name_ident

    integer :: loc

    loc=0
    ! Unfortunately the following line caused an internal compiler error in
    ! gfortran.
    !loc=minval(ident_list%ident, (ident_list%ident==ident))
    do loc=1,size(ident_list)
       if (ident_list(loc)%name==name) then
          goto 43
       end if
    end do
    
    ! Fallthrough case for unknown fields.
    name_ident = 0

    return
    
43  name_ident = ident_list(loc)%ident
    
  end function name_ident

  function domain_is_2d() result(bool)
    !!< Is the domain pseudo2d or not?
    logical :: bool

    bool = .false.
    if (pseudo2d_coord > 0 .and. pseudo2d_coord < 4) bool = .true.
  end function domain_is_2d

  function domain_is_2d_x() result(bool)
    !!< Is the domain pseudo2d in the x direction?
    logical :: bool
    bool = (pseudo2d_coord == 1)
  end function domain_is_2d_x

  function domain_is_2d_y() result(bool)
    !!< Is the domain pseudo2d in the y direction?
    logical :: bool
    bool = (pseudo2d_coord == 2)
  end function domain_is_2d_y
  
  function domain_is_2d_z() result(bool)
    !!< Is the domain pseudo2d in the z direction?
    logical :: bool
    bool = (pseudo2d_coord == 3)
  end function domain_is_2d_z
end module global_parameters
