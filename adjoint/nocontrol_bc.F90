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
SUBROUTINE NOCONTROL_BC(X,Y,Z,                                   &
     NONODS,NOBCU,NOBCV,NOBCW,MAXGRAD,D3,ITSTIME,        &
     MAXREC,MRECORD,NOCVA,G,IMEM_CVA)
   use FLDebug
  use signals
  use spud
  IMPLICIT NONE
  REAL    ::MAXGRAD
  LOGICAL ::D3
  INTEGER ::NOBCU,NOBCV,NOBCW !passed dawn from the forward model
  INTEGER  ::NONODS,NOCVA,ITSTIME,MAXREC
  REAL    ::X(NONODS),Y(NONODS),Z(NONODS)
  REAL    ::G(NOCVA)    ! the vector of the gradient
  INTEGER    ::MRECORD(MAXREC)
  INTEGER    ::LSGRAD
  INTEGER   ::IMEM_CVA(NOCVA)

  ! Local....
  REAL,DIMENSION(:),ALLOCATABLE::G_IN,G_OUT
  INTEGER ::K1,K2        ! the records for writing to G
  INTEGER ::II,NOD
  REAL    ::TT
  REAL    ::X_in,X_out

  if(have_option(trim("/model/fluids/adjoint/bctest"))) then
     X_in = -7.5
  else 
     X_in = 0.0
  endif

  ! Local total size of G at each time step
  LSGRAD = NOBCU+NOBCV
  IF(D3) THEN
     LSGRAD = LSGRAD + NOBCW
  ENDIF

  ALLOCATE(G_IN(LSGRAD))
  ALLOCATE(G_OUT(LSGRAD))
  G_IN(1:LSGRAD) = 0.0
  G_OUT(1:LSGRAD) = 0.0


  ewrite(3,*) 'NOCVA in SUB_grad',NOCVA
  K1 = MRECORD(ITSTIME)
  K2 = K1+NOBCU-1
  ewrite(3,*) 'bcu'
  DO II=1,NOBCU
     NOD=IMEM_CVA(K1+II-1)
     ewrite(3,*) 'II,NOCVA,nod,nonods',II,NOCVA,nod,nonods
     IF( ABS(X(NOD)-X_in).GT.1.0E-6 ) THEN
!     IF( (ABS(X(NOD)-0.0).GT.1.0E-6).AND.(ABS(X(NOD)-30000.0).GT.1.0E-6) ) THEN
         G(K1+II-1)=0.0
!        IF(ITSTIME.NE.1) G(K1+II+LSGRAD-1)=0.0
     ENDIF
  if( abs(x(nod)-X_in).lt.1.0e-6) then
    G_IN(II) = G(K1+II-1)
  endif 
  if( abs(x(nod)-30000.0).lt.1.0e-6) then
    G_OUT(II) = G(K1+II-1)
  endif 

  ENDDO

  K1 = K2+1             
  K2 = K1+NOBCV-1
  ewrite(3,*) 'bcv'
  DO II=1,NOBCV
     NOD=IMEM_CVA(K1+II-1)
!     IF( ABS(X(NOD)-0.0).GT.1.0E-6 ) THEN
     IF( (ABS(X(NOD)-X_in).GT.1.0E-6).AND.(ABS(X(NOD)-30000.0).GT.1.0E-6) ) THEN
        G(K1+II-1)=0.0
!        IF(ITSTIME.NE.1) G(K1+II+LSGRAD-1)=0.0
     ENDIF
  if( abs(x(nod)-X_in).lt.1.0e-6) then
    G_IN(NOBCU+II) = G(K1+II-1)
  endif 
  if( abs(x(nod)-30000.0).lt.1.0e-6) then
    G_OUT(NOBCU+II) = G(K1+II-1)
  endif 

  ENDDO

  IF(D3) THEN
     K1 = K2+1           
     K2 = K1+NOBCW-1
     DO II=1,NOBCW
        NOD=IMEM_CVA(K1+II-1)
!     IF( ABS(X(NOD)-0.0).GT.1.0E-6 ) THEN
     IF( (ABS(X(NOD)-X_in).GT.1.0E-6).AND.(ABS(X(NOD)-30000.0).GT.1.0E-6) ) THEN
           !     IF(X(NOD).NE.0.0) THEN
           G(K1+II-1)=0.0
!           IF(ITSTIME.NE.1) G(K1+II+LSGRAD-1)=0.0
        ENDIF
  if( abs(x(nod)-X_in).lt.1.0e-6) then
    G_IN(NOBCU+NOBCV+II) = G(K1+II-1)
  endif 
  if( abs(x(nod)-30000.0).lt.1.0e-6) then
    G_OUT(NOBCU+NOBCV+II) = G(K1+II-1)
  endif 
    ENDDO
  ENDIF

  MAXGRAD = 0.0
  DO II =1,NOCVA
     MAXGRAD = MAX(MAXGRAD,ABS(G(II)))
     if(ABS(G(II)).gt.3.0) then
          ewrite(3,*) 'MAXGRAD in SUB_noncontrol_bc',ABS(G(II))
     endif
  ENDDO

  if(itstime.eq.1) then
    open(1,file='gradient_inlet.dat',status='replace')
    open(3,file='gradient_inlet1.dat',status='replace')
    open(2,file='gradient_outlet.dat',status='replace')
    write(1,*) itstime
    write(2,*) itstime
  else
    open(1,file='gradient_inlet.dat',position='append')
    open(3,file='gradient_inlet1.dat',position='append')
    open(2,file='gradient_outlet.dat',position='append')
    write(1,*) itstime
    write(2,*) itstime
  endif

  write(1,*) (G_IN(II),II=1,LSGRAD)
  write(2,*) (G_OUT(II),II=1,LSGRAD)
  write(3,*) itstime,G_OUT(1)
 close(1)
 close(2)
 close(3)

  ewrite(3,*) 'deallocate1...'
     DEALLOCATE(G_IN)
  ewrite(3,*) 'deallocate2...'
     DEALLOCATE(G_OUT)
  ewrite(3,*) 'deallocate3...'

END SUBROUTINE NOCONTROL_BC

!================================


SUBROUTINE NOCONTROL_BC_FS(X,Y,Z,                        &
     NONODS,NOBCH,MAXGRAD,D3,ITSTIME,NTIME,              &
     MAXREC,MRECORD,NOCVA,G,IMEM_CVA)
   use FLDebug
  IMPLICIT NONE
  REAL    ::MAXGRAD
  LOGICAL ::D3
  INTEGER ::NOBCH   !passed dawn from the forward model
  INTEGER ::NONODS,NOCVA,ITSTIME,MAXREC,NTIME
  REAL    ::X(NONODS),Y(NONODS),Z(NONODS)
  REAL    ::G(NOCVA)    ! the vector of the gradient
  INTEGER ::MRECORD(MAXREC)
  INTEGER ::LSGRAD
  INTEGER ::IMEM_CVA(NOCVA)


  ! Local....
  REAL,DIMENSION(:),ALLOCATABLE::G_IN,G_OUT
  INTEGER ::K1,K2        ! the records for writing to G
  INTEGER ::II,NOD
  REAL    ::TT

  ! Local total size of G at each time step
  LSGRAD = NOBCH

  ALLOCATE(G_IN(LSGRAD))
!  ALLOCATE(G_OUT(LSGRAD))
  G_IN(1:LSGRAD) = 0.0
!  CALL RCLEAR(G_OUT,LSGRAD)


  ewrite(3,*) 'NOCVA in SUB_grad_fs',NOCVA
  K1 = MRECORD(ITSTIME)
  K2 = K1+NOBCH-1
  ewrite(3,*) 'bch'
  DO II=1,NOBCH
     NOD=IMEM_CVA(K1+II-1)
     ewrite(3,*) 'II,NOCVA,nod,nonods',II,NOCVA,nod,nonods
     IF( ABS(X(NOD)-0.0).GT.1.0E-6 ) THEN
!     IF( (ABS(X(NOD)-0.0).GT.1.0E-6).AND.(ABS(X(NOD)-30000.0).GT.1.0E-6) ) THEN
         G(K1+II-1)=0.0
!        IF(ITSTIME.NE.1) G(K1+II+LSGRAD-1)=0.0
     ENDIF
  if( abs(x(nod)-0.0).lt.1.0e-6) then
    G_IN(II) = G(K1+II-1)
  endif 
!  if( abs(x(nod)-30000.0).lt.1.0e-6) then
!    G_OUT(II) = G(K1+II-1)
!  endif 

  ENDDO

  if(itstime.eq.1) then
    open(1,file='gradient_inlet.dat',status='replace')
    open(3,file='gradient_inlet1.dat',status='replace')
    write(1,*) itstime
    write(2,*) itstime
  else
    open(1,file='gradient_inlet.dat',position='append')
    open(3,file='gradient_inlet1.dat',position='append')
    write(1,*) itstime
  endif

  write(1,*) (G_IN(II),II=1,LSGRAD)
!  write(2,*) (G_OUT(II),II=1,LSGRAD)
  write(3,*) ntime-itstime+1,G_IN(1)
 close(1)
 close(2)
 close(3)
 

tt = (itstime-1)*200.0

!IF(tt.lt.20000.0.or.tt.gt.64000.0) THEN
!IF(tt.le.21600.0.or.tt.gt.64800.0) THEN
!  K1 = MRECORD(ITSTIME)
!  K2 = K1+NOBCH-1
!  DO II=1,NOBCH
!    G(K1+II-1)=0.0
!    G_IN(II) = G(K1+II-1)
!  ENDDO
!ENDIF


  if(itstime.eq.1) then
    open(3,file='gradient_inlet2.dat',status='replace')
  else
    open(3,file='gradient_inlet2.dat',position='append')
  endif

  write(3,*) ntime-itstime+1,G_IN(1)

  close(3)

 IF( (ITSTIME-NINT(NTIME/2.0)).LT.1.0E-6 ) THEN
  MAXGRAD = 0.0
  DO II =1,NOCVA
     MAXGRAD = MAX(MAXGRAD,ABS(G(II)))
  ENDDO
 ENDIF
  ewrite(3,*) 'deallocate1...'
     DEALLOCATE(G_IN)
  ewrite(3,*) 'deallocate2...'
 !    DEALLOCATE(G_OUT)
  ewrite(3,*) 'deallocate3...'


END SUBROUTINE NOCONTROL_BC_FS

