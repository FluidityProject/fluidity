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

module detector_data_types

  use fldebug
  use global_parameters, only : FIELD_NAME_LEN
  
  implicit none
  
  private
  
  public :: detector_type, detector_linked_list, STATIC_DETECTOR, LAGRANGIAN_DETECTOR

  integer, parameter :: STATIC_DETECTOR=1, LAGRANGIAN_DETECTOR=2  

  !! Type for caching detector position and search information.
  type detector_type
     !! Physical location of the detector.
     real, dimension(:), allocatable :: position
     !! Name of the detector in input and output.
     character(len=FIELD_NAME_LEN) :: name 
     !! Whether the detector is on the current processor. 
     logical :: local
     !! Element number in which the detector lies.
     integer :: element
     !! Local coordinates of the detector in that element.
     real, dimension(:), allocatable :: local_coords
     !! Whether the detector is static or Lagrangian.
     integer :: type = STATIC_DETECTOR
     !! Bisected time step used when moving the Lagrangian detectors
     real :: dt
     !! Identification number indicating the order in which the detectors are read
     integer :: id_number
     !! Processor that owns initially the detector (for parallel)
     integer :: initial_owner= -1
     TYPE (detector_type), POINTER :: next=> null()
     TYPE (detector_type), POINTER :: previous=> null() 
  end type detector_type
  
  type detector_linked_list
     integer :: length=0
     TYPE (detector_type), POINTER :: firstnode => null()
     TYPE (detector_type), POINTER :: lastnode => null()
  end type detector_linked_list

end module detector_data_types
