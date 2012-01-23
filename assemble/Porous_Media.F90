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

module porous_media

  use fldebug
  use state_module
  use fields
  use spud
  implicit none
  
  private

  public :: calculate_porous_media_absorption

  contains

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
       FLExit('Phase '//int2str(i)//' is missing a viscosity field.')
    end if
    s_viscosity = extract_scalar_field_from_tensor_field(t_viscosity, 1, 1)
    ewrite_minmax(s_viscosity)

    ! check permeability field and appropriately compute absorption term
    if (have_option("/porous_media/scalar_field::Permeability")) then
      ! SCALAR PERMEABILITY
      s_permeability = extract_scalar_field(state(1), "Permeability", stat)
      do j=1,absorption%dim
        do k=1,node_count(absorption)
          call set(absorption, j, k, node_val(s_viscosity, k)/node_val(s_permeability, k))
        end do
      end do
    elseif (have_option("/porous_media/vector_field::Permeability")) then
      ! VECTOR PERMEABILITY
      v_permeability => extract_vector_field(state(1), "Permeability", stat)
      do j=1,absorption%dim
        s_permeability = extract_scalar_field_from_vector_field(v_permeability, j)
        do k=1,node_count(absorption)
          call set(absorption, j, k, node_val(s_viscosity, k)/node_val(s_permeability, k))
        end do
      end do
    elseif (have_option("/porous_media/tensor_field::Permeability")) then
      ! TENSOR PERMEABILITY
      t_permeability => extract_tensor_field(state(1), "Permeability", stat)
      FLExit("Cannot currently solve porous media problems with a full tensor permeability.")
    end if
    if (stat/=0) then
      FLExit("Porous media problems require a permeability field.")
    end if
    stat = 0
  end subroutine calculate_porous_media_absorption

end module porous_media

