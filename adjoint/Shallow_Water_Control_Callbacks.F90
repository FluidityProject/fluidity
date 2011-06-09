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

module shallow_water_adjoint_controls
  use fields
  use global_parameters, only : PYTHON_FUNC_LEN, OPTION_PATH_LEN
  use spud
  use state_module
  use python_state

    implicit none

    private

    public :: shallow_water_adjoint_timestep_callback 

    contains 
    
    subroutine shallow_water_adjoint_timestep_callback(state, timestep, functional_name)
      use global_parameters, only: OPTION_PATH_LEN
      use state_module
      type(state_type), dimension(:), intent(inout) :: state
      integer, intent(in) :: timestep
      character(len=OPTION_PATH_LEN), intent(in) :: functional_name

      character(len=OPTION_PATH_LEN) :: field_name, control_type, control_deriv_name, field_deriv_name
      integer :: nb_controls
      integer :: i
      type(scalar_field), pointer :: sfield
      type(vector_field), pointer :: vfield
      type(tensor_field), pointer :: tfield
      
      nb_controls = option_count("/adjoint/controls/control")
      do i = 0, nb_controls-1
        call get_option("/adjoint/controls/control[" // int2str(i) //"]/name", control_deriv_name)
        control_deriv_name = trim(control_deriv_name) // "_Derivative"
        call get_option("/adjoint/controls/control[" // int2str(i) //"]/type/name", control_type)
        call get_option("/adjoint/controls/control[" // int2str(i) //"]/type::" // trim(control_type) // "/field_name", field_name)
      end do

    end subroutine shallow_water_adjoint_timestep_callback 

end module shallow_water_adjoint_controls
