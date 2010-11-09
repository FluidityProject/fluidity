!    Copyright (C) 2009 Imperial College London and others.
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

subroutine test_pade_equation_of_state
  !!< Test the nonlinear equation of state against the check values in the
  !!< original publication. The checks are from p740 of:
  !!< McDougall, T.J., D.R. Jackett, D.G. Wright, and R. Feistel, 2003:
  !!< Accurate and Computationally Efficient Algorithms for Potential
  !!< Temperature and Density of Seawater. J. Atmos. Oceanic Technol., 20,
  !!< 730-741. DOI: 10.1175/1520-0426(2003)20

  use equation_of_state
  use unittest_tools
  implicit none
  
  real :: rho

  call mcD_J_W_F2002(rho, 25., 35., invert_distance(2000.))

  call report_test("[Pade equation of state: check point 1.]", &
      & .not.fequals(1000*rho+1000, 1031.654229, 1.e-6), .false., &
      & "Equation of state returned wrong value")

  call mcD_J_W_F2002(rho, 20., 20., invert_distance(1000.))

  call report_test("[Pade equation of state: check point 2.]", &
      & .not.fequals(1000*rho+1000, 1017.726743, 1.e-6), .false., &
      & "Equation of state returned wrong value")

  call mcD_J_W_F2002(rho, 12., 40., invert_distance(8000.))

  call report_test("[Pade equation of state: check point 3.]", &
      & .not.fequals(1000*rho+1000, 1062.928258, 1.e-6), .false., &
      & "Equation of state returned wrong value")

contains
  
  function invert_distance(pressure)
    real :: invert_distance
    real, intent(in) :: pressure

    invert_distance=pressure/( 9.81*1000.0*1.0e-4 )

  end function invert_distance

end subroutine test_pade_equation_of_state
