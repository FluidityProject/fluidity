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

module automatic_timestep

  use fldebug
  use parallel_tools

  implicit none

  private

  public :: auttim

contains

  SUBROUTINE AUTTIM(PARA,OLDU,NEWU,NONODS,NNODP,&
       &    OLDDT,NEWDT,AUTACC,AUTGOB,PLIMDT,LIMIDT,TAGAIN)
    INTEGER PARA, NONODS, NNODP
    REAL OLDU(NONODS),NEWU(NONODS)
    REAL OLDDT,NEWDT, AUTACC,AUTGOB,PLIMDT,LIMIDT
    REAL T1OVT2
    ! AUTACC is the automatic accuracy it must be greater than one 
    ! and less than 2. just above one is most accurat EG 1.05.
    ! If AUTGOB * (accuracy suggested by AUTACC) is not reached 
    ! then repeat time step.
    ! If TAGAIN repeat time step, else generate a new time step NEWDT
    ! TAGAIN is always false if AUTGOB.le.0.
    ! PLIMDT=the previous limit on the actual time step 
    ! needed to obtain min accuracy requirments. 
    REAL LANDA,TOLER
    LOGICAL TAGAIN
    INTEGER I, NOD
    REAL RNEWU, ROLDU, RRR
    ! This sub finds the next time step DT automatically

    TAGAIN=.FALSE.
    !
    TOLER=0.
    do  I=1,NNODP! Was loop 5
       TOLER =MAX(TOLER, ABS(OLDU(I)))
    end do ! Was loop 5
    IF(PARA.EQ.1) CALL ALLMAX(TOLER)
    IF(TOLER.LT.1.E-2) TOLER=1.E-2
    TOLER=0.01*TOLER
    !
    !
    T1OVT2=0.
    do  NOD=1,NNODP! Was loop 10
       RNEWU=ABS(NEWU(NOD))
       ROLDU=ABS(OLDU(NOD))
       T1OVT2=MAX(T1OVT2, &
            &          (ABS(RNEWU-ROLDU)+ROLDU)/MAX(ROLDU,TOLER)  )
    end do ! Was loop 10
    IF(PARA.EQ.1) CALL ALLMAX(T1OVT2)
    ! This stops it from using time steps more 
    ! then about 3 times that of previous step.
    RRR=1.+0.3*(AUTACC-1.)
    IF(T1OVT2.LT.RRR) T1OVT2=RRR
    !        ewrite(3,*)'AUTTIM: T1OVT2=',T1OVT2
    !
    LANDA=LOG(T1OVT2)/OLDDT
    !          ewrite(3,*)'AUTTIM: LANDA=',LANDA
    RRR=LOG(AUTACC)/LANDA
    NEWDT=MIN(NEWDT, RRR)
    IF(AUTGOB.GT.0) &
         &         LIMIDT=MIN(LIMIDT, &
         &         LOG(AUTGOB*(AUTACC-1.)+1.)/LANDA)
    ewrite(3,*)'TOLER,OLDDT,NEWDT,AUTACC,AUTGOB,LIMIDT:',&
         &             TOLER,OLDDT,NEWDT,AUTACC,AUTGOB,LIMIDT
    ewrite(3,*)'SUGGESTED NEW TIME STEP HERE IS:',RRR
    !
    IF(NEWDT.GT.PLIMDT) THEN
       TAGAIN=.TRUE.
    ELSE
       TAGAIN=.FALSE.
    ENDIF
    IF(AUTGOB.LE.0) TAGAIN=.FALSE.
    !
    RETURN
  end subroutine auttim

end module automatic_timestep
