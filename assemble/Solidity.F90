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
!    q
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

module solidity
  !! This module contains the options and material properties used
  !! when running FLUIDITY in the SOLIDITY mode
  use fldebug
  use futils
  use spud
  use global_parameters, only: OPTION_PATH_LEN
  implicit none

  private
  public :: GET_SOLIDITY_OPTIONS

contains

  SUBROUTINE GET_SOLIDITY_OPTIONS(SOLIDS, MKCOMP)
    ! This subroutine is now legacy and only acting as a back compatibility
    ! link to new options for solids.
    ! Description                                   Programmer      Date
    ! ==================================================================
    ! Original version ..............................  GSC    2006-03-29
    ! Added input of cohesion and angle of friction..  CRGW   2006-06-10
    ! Stripped down legacy version...................  CRGW   2008-04-14

    ! Outputs
    INTEGER, INTENT(OUT) :: SOLIDS       !! Solid material modeling on or off
    INTEGER, INTENT(OUT) :: MKCOMP       !! Compressibility option (set to zero to turn off)

    integer :: i, f
    character(len=OPTION_PATH_LEN) :: thismaterial_phase
      
    ! Set default (off) values
    SOLIDS       = 0
    MKCOMP       = 0

    do i = 1,option_count('/material_phase')
      thismaterial_phase = '/material_phase['//int2str(i-1)//']'

      if (have_option(trim(thismaterial_phase)//'/vector_field::Velocity&
          &/prognostic/tensor_field::Elasticity')) then
        SOLIDS = 1
        if (have_option(trim(thismaterial_phase)//'/vector_field::Velocity&
          &/prognostic/tensor_field::Viscosity')) then
          SOLIDS = 2
        end if
        if ((have_option(trim(thismaterial_phase)//'/vector_field::Velocity&
          &/prognostic/scalar_field::Cohesion')).or.&
            (have_option(trim(thismaterial_phase)//'/vector_field::Velocity&
          &/prognostic/scalar_field::FrictionAngle'))) then
          SOLIDS = 3
        end if
      end if

      if (have_option(trim(thismaterial_phase)//'/equation_of_state&
          &/multimaterial/miegrunneisen')) then
        MKCOMP = 3
      end if

    end do

  END SUBROUTINE GET_SOLIDITY_OPTIONS

end module solidity

