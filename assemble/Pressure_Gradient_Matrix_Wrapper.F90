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

module pressure_gradient_matrix_wrapper
  use fldebug
  use state_module
  use sparse_tools
  use spud
  use divergence_matrix_cv
  use divergence_matrix_cg
  use global_parameters, only: OPTION_PATH_LEN
  use fields
  implicit none 

  private
  public :: assemble_pressure_gradient_matrix

contains

  !************************************************************************
  subroutine assemble_pressure_gradient_matrix(getc12, mkcomp, CT_m, CTP_m, &
                                               state, istate, pressure_option_path)
    ! =============================================================
    ! Subroutine to construct the matrices CMC, C1/2/3T and 
    ! (if MKCOMP = 3) C1/2/3PT.
    ! =============================================================
  
    ! =============================================================
    ! inputs 
    ! =============================================================

    ! bucket full of fields
    logical, intent(in) :: getc12
    integer, intent(in) :: mkcomp

    type(state_type), dimension(:), intent(inout) :: state

    integer, intent(in) :: istate

    character(len=OPTION_PATH_LEN), intent(in) :: pressure_option_path

    ! the pressure gradient and compressible gradient matrices
    type(block_csr_matrix), intent(inout) :: CT_m, CTP_m

    type(scalar_field), pointer :: pressure
    type(vector_field), pointer :: velocity

    integer :: i

    ewrite(2,*) 'In assemble_pressure_gradient_matrix'

    pressure=>extract_scalar_field(state(istate), "Pressure")
    velocity=>extract_vector_field(state(istate), "Velocity")

    ! do we want C1/2/3T (i.e. GETC12)?
    IF(GETC12) THEN

      if(have_option(trim(pressure_option_path)//&
            '/prognostic/spatial_discretisation/control_volumes')) then

        call assemble_divergence_matrix_cv(CT_m, state(istate), &
                                           test_mesh=pressure%mesh, field=velocity)

        if(mkcomp.eq.3) then
          call assemble_compressible_divergence_matrix_cv(CTP_m, state)
        end if

      else

        call assemble_divergence_matrix_cg(CT_m, state(istate), & 
                                          test_mesh=pressure%mesh, field=velocity, &
                                          option_path=trim(pressure_option_path))


      end if

      do i = 1, CT_m%blocks(2)
        ewrite_minmax(CT_m%val(1,i)%ptr(:))
      end do

      if(mkcomp.eq.3) then
        do i = 1, CTP_m%blocks(2)
          ewrite_minmax(CTP_m%val(1,i)%ptr(:))
        end do
      end if

    end if

  end subroutine assemble_pressure_gradient_matrix

end module pressure_gradient_matrix_wrapper

