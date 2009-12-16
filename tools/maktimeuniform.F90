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
subroutine fprint_backtrace
end subroutine fprint_backtrace

PROGRAM MAIN 
  use FLDebug
  IMPLICIT NONE 
  REAL MAXX,MAXY,MINY,MINX
  REAL DT,RVALUE,ACCTIM,ACCTIM2,VALUE,VALUE2,VALINT
  REAL TIME,ACCTI1,LTIME,W
  INTEGER NONODS,FREDOP,TOTELE,ELE,ILINE,ICOUNT
  INTEGER I,SOMTIM,J,NUMBER,IREP,NUMITS,NTIME,NT,IL
  LOGICAL MAKE2
  ! THIS PROGRAM PLOTS VELS FOR LINEAR MODEL

  !

  ewrite(3,*) 'enter time step we want'
  ewrite(3,*) 'IF YOU ENTER -VE THEN MAKE THE PRODUCED'
  ewrite(3,*) 'DATA A POWER OF 2'
  ewrite(3,*) 'BY SELECTING A SMALLER TIME STEP' 
  read(*,*) DT 
  !
  MAKE2=(DT.LT.0)
  DT=ABS(DT) 
  !
  !
  !
  !
  OPEN(2,FILE='forlog.data',STATUS= 'UNKNOWN')
  !
  ICOUNT=1
  read(2,*,end=3000) ACCTI1,VALUE
  do iline=1,100000000
     read(2,*,end=3000) ACCTIM,VALUE
     ICOUNT=ICOUNT+1
  end do
3000 CONTINUE
  ewrite(3,*) 'ACCTIM,RVALUE:',ACCTIM,VALUE
  !
  CLOSE(2)
  !
  ! MODIFY DT to make it a power of 2...
  ewrite(3,*) 'MAKE2=',MAKE2
  IF(MAKE2) THEN
     LTIME=ACCTIM-ACCTI1
     NTIME= LTIME/DT + 1.E-25 
     ! obtain an NTIME that is a power of two
     NT=1
     do I=1,10000
        NT=NT*2
        IF(NT.GE.NTIME) GOTO 1000
     END DO
1000 CONTINUE
     NTIME=NT
     DT=LTIME/(NTIME-1)
     ewrite(3,*) 'recalculated NTIME,DT:',NTIME,DT
     ewrite(3,*) 'the power of two=',I 
  ENDIF
  !
  !
  OPEN(2,FILE='forlog.data',    STATUS='UNKNOWN')
  OPEN(9,FILE='foruni.data',STATUS='UNKNOWN')
  ! 
  ewrite(3,*) 'this has produced a file called foruni.data'
  !
  icount=1
  read(2,*,end=4000) ACCTIM, VALUE
  !        write(9,*)         ACCTIM, VALUE
  write(9,*)         0.0, VALUE
  TIME=ACCTIM+DT
  do  iline=1,100000000! Was loop 10
     read(2,*,end=4000) ACCTIM2,VALUE2
     !        ewrite(3,*) 'time,ACCTIM2,VALUE2:',time,ACCTIM2,VALUE2
     !
     do IL=1,100000
        IF((TIME.GT.ACCTIM).AND.(TIME.LE.ACCTIM2)) THEN
           ! INTERPOLATE VALUE... 
           W=(TIME-ACCTIM)/(ACCTIM2-ACCTIM)
           VALINT=(1.-W)*VALUE + W*VALUE2
           icount=icount+1
           !             if(icount.le.ntime) write(9,*) TIME,VALINT
           if(icount.le.ntime) write(9,*) TIME-ACCTI1,VALINT
           !             ewrite(3,*) 'TIME-ACCTI1,VALINT:',TIME-ACCTI1,VALINT
           ! Next time value...
           TIME=TIME+DT
        ELSE
           IF(TIME.GE.ACCTIM2) GOTO 5000
        ENDIF
     END DO
5000 CONTINUE
     !
     ACCTIM=ACCTIM2
     VALUE =VALUE2 
  end do ! Was loop 10
  !
4000 CONTINUE
  ewrite(3,*) 'ACCTIM2,VALUE2:',ACCTIM2,VALUE2
  !
  STOP
END PROGRAM MAIN
!
!
!


!
