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

module copy_outof_into_state
  !! This module enables the multiphase prototype code to interact with state
  !! First copying everything required from state to the old variable space
  !! Second copying everything back to state for output 

  use FLDebug
  use state_module
  use populate_state_module
  use diagnostic_variables
  use diagnostic_fields_wrapper
  use global_parameters, only: option_path_len
  use diagnostic_fields_new, only : &
    & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
    & check_diagnostic_dependencies
  use spud
  implicit none
  
  private
  
  public :: copy_outof_state, copy_into_state
  
  contains
  
  subroutine copy_outof_state()
  
  end subroutine copy_outof_state
  
  subroutine copy_into_state()
  
  end subroutine copy_into_state

end module copy_outof_into_state