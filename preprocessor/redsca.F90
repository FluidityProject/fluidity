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

module redsca_module

  use parallel_tools
  use state_module
  use fetools
  use spud
  use flcomms_module
  use fldebug
  use global_parameters, only : halo_tag, halo_tag_p

  implicit none

contains
  SUBROUTINE REDSCA(PROCNO,&
       &     NNODP, NNODPP, NDPSET, &
       &     OPTSOU,NSOUPT,MXNSOU,FIESOU)

    INTEGER NDPSET
    INTEGER PROCNO,NNODP,NNODPP
    ! This is for distributed sources for a field...
    INTEGER OPTSOU,NSOUPT,MXNSOU
    INTEGER FIESOU(MXNSOU)
    
    INTEGER I,J,NCARS
    
    !     state/options stuff -- cjc
    integer :: mat, n_pressure, tmpint, stat

    if(IsParallel()) then
       nnodp = get_nowned_nodes(halo_tag)
       nnodpp = get_nowned_nodes(halo_tag_p)
    else
       nnodp = 0
       nnodpp = 0
    end if

    NDPSET = 0 ! to avoid it being unitialised when there's no prognostic pressure
    n_pressure = option_count("/material_phase/scalar_field::Pressure/prognostic")
    if(n_pressure>1) then ! just a test to make sure we're not doing something we didn't mean to
       FLAbort("Multiple prognostic pressures... don't know what to do.")
    end if
    if(n_pressure==1) then  ! no point in looping if there's no prognostic pressure
       do mat = 0, option_count("/material_phase")-1
          call get_option("/material_phase["//int2str(mat)//"]/scalar_field::Pressure&
               &/prognostic/reference_node", tmpint, stat=stat)
          if(stat==0) then ! covers cases of there being no prognostic pressure in this
             NDPSET=tmpint  ! material_phase and of it not being set
          end if
       end do
    end if

    ! Ensure that only processor number 1 is the only one to set a
    ! reference pressure node.
    IF(ProcNo.NE.1) THEN
       NDPSET = 0
    END IF

  END SUBROUTINE REDSCA

end module redsca_module
