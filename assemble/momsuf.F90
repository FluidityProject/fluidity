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

module shear_stress

  use fldebug
  use position_in_matrix

  implicit none

  private

  public :: calculate_shear_stress

contains

  SUBROUTINE calculate_shear_stress(ACCTIM,U,V,W,SCFACTH0, &
       VECX,VECY,VECZ, NONODS,STOTEL,SNLOC,SNGI, &
       D3,DCYL,DSPH, &
       SNDGLN, &
       SN,SNLX,SNLY,SWEIGH, &
       X,Y,Z,R0,D0, &
       BIGM,NBIGM, &
       THETA,FINDRM,COLM,NCOLM,CENTRM,ITSITI,ITINOI, &
! make implicit in pressure as well
       DT,ML,SHEARSTRESSX,GOT_SHEARSTRESSX,SHEARSTRESSY,GOT_SHEARSTRESSY, &
       SHEARSTRESSZ,GOT_SHEARSTRESSZ,MEANSHEARSTRESSX, &
       GOT_MEANSHEARSTRESSX,MEANSHEARSTRESSY,GOT_MEANSHEARSTRESSY, &
       MEANSHEARSTRESSZ,GOT_MEANSHEARSTRESSZ)

! boundary elements - in which I,J=1(x), I,J=3(z). WE may use 
! its symmetry at a latter stage. 
! R0=the minimum distance from ocean bed to centre of earth.
! D0=H/r0 where H is the height of ocean above R0 (H not used here)
! For dimensional run D0=0.
! NB for cylindrical coords (r,z)           corresponds to (y,x)
! NB for spherical   coords (\psi,\theta,r) corresponds to (x,y,z)
    REAL    TWOPI,R0,SCFACTH0,ACCTIM
    INTEGER NONODS,STOTEL,SNLOC,SNGI
    REAL    U(NONODS),V(NONODS),W(NONODS)
    REAL    VECX(NONODS),VECY(NONODS),VECZ(NONODS)
    REAL    X(NONODS),Y(NONODS),Z(NONODS)
    LOGICAL D3,DCYL,DSPH
    INTEGER SNDGLN(STOTEL*SNLOC)
    REAL    SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
    REAL    SWEIGH(SNGI)
    REAL    DT,ML(NONODS)
    REAL    SHEARSTRESSX(NONODS),SHEARSTRESSY(NONODS),SHEARSTRESSZ(NONODS)
    LOGICAL GOT_SHEARSTRESSX,GOT_SHEARSTRESSY,GOT_SHEARSTRESSZ
    REAL    MEANSHEARSTRESSX(NONODS),MEANSHEARSTRESSY(NONODS),MEANSHEARSTRESSZ(NONODS)
    LOGICAL GOT_MEANSHEARSTRESSX,GOT_MEANSHEARSTRESSY,GOT_MEANSHEARSTRESSZ

    INTEGER SELE,L,GI,Q,P,GLOBI,GLOBJ,COUNT,LOWER,UPPER
    INTEGER LEGL,LIGL,LJGL,HIT,HIT2,NOD
    REAL    DXDLX,DXDLY,DYDLX,DYDLY,DZDLX,DZDLY
    REAL    D0,DD0
    REAL    RGI,A,B,C,RNN
    INTEGER IGL,ILOC,JLOC
    REAL    NN,NNDRG,NNX,NNY,NNZ
    REAL    DRAG
    INTEGER NBIGM
    REAL    BIGM(NBIGM),THETA,DRAG_THETA
    INTEGER DRAG_OPTION
    INTEGER NCOLM
    INTEGER FINDRM(NONODS+1),COLM(NCOLM),CENTRM(NONODS),ICENT
    integer itsiti,itinoi
    INTEGER IBL11,IBL12,IBL13,IBL21,IBL22,IBL23,IBL31,IBL32,IBL33
    LOGICAL BLKSYM,LUMPDRAG
    REAL    CD,SUM,lmpdrg
    REAL    spint,rub,tempssx,tempssy,tempssz       
    integer rubi      
    logical got_tidal,got_fsopts
    CHARACTER*4096 data_file
    real N_TIME_AV
    save N_TIME_AV
! *****************************************************************

    ewrite(1,*) 'in calculate_shear_stress'

    if(ITSITI.EQ.ITINOI) then
       CD=0.003 

! Compute the max bed shear stress
      ewrite(3,*) "in max bed shear stress calculation"
      ewrite(3,*) "GOT_SHEARSTRESSX:", GOT_SHEARSTRESSX
       IF(GOT_SHEARSTRESSX) THEN ! first find spint
          spint=0.0
          data_file = ' '
          data_file(1:13) = 'FSoptions.dat'
          INQUIRE(file=data_file, exist=got_fsopts)
          if(got_fsopts) then
998          FORMAT(I9)
999          FORMAT(1E15.7)
             OPEN(555,status='unknown',file=data_file)
             READ(555,998) rubi
             READ(555,999) rub
             READ(555,999) spint
             CLOSE(555)
          endif
       ENDIF
       ewrite(3,*) "GOT_SHEARSTRESSX,GOT_SHEARSTRESSY,GOT_SHEARSTRESSZ:", &
              GOT_SHEARSTRESSX,GOT_SHEARSTRESSY,GOT_SHEARSTRESSZ

       IF( GOT_SHEARSTRESSX.AND.GOT_SHEARSTRESSY.AND.GOT_SHEARSTRESSZ ) then
          do globi=1,nonods
             IF(acctim.lt.spint) THEN
                SHEARSTRESSX(GLOBI) = 0.
                SHEARSTRESSY(GLOBI) = 0.
                SHEARSTRESSZ(GLOBI) = 0.
             ENDIF
             tempssx = 1023.0*CD*SQRT(U(GLOBI)**2 +V(GLOBI)**2 + &
                  W(GLOBI)**2)*U(GLOBI)
             tempssy = 1023.0*CD*SQRT(U(GLOBI)**2 +V(GLOBI)**2 + &
                  W(GLOBI)**2)*V(GLOBI)
             tempssz = 1023.0*CD*SQRT(U(GLOBI)**2 +V(GLOBI)**2 + &
                  W(GLOBI)**2)*W(GLOBI)
             if( ( tempssx**2 + tempssy**2 + tempssz**2 ).gt. &
                  ( ( SHEARSTRESSX(GLOBI)**2 + SHEARSTRESSY(GLOBI)**2 + &
                  SHEARSTRESSZ(GLOBI)**2 )) ) THEN
                SHEARSTRESSX(GLOBI) = tempssx
                SHEARSTRESSY(GLOBI) = tempssy
                SHEARSTRESSZ(GLOBI) = tempssz
             endif
          enddo
       ELSEIF ( GOT_SHEARSTRESSX ) THEN ! If just got one component then just output the magnitude and not direction.
          do globi=1,nonods     
             IF(acctim.lt.spint) THEN
                SHEARSTRESSX(GLOBI) = 0.
             ENDIF
             SHEARSTRESSX(GLOBI) = MAX(SHEARSTRESSX(GLOBI), &
                  abs(1023.0*CD*(U(GLOBI)**2 + V(GLOBI)**2 + W(GLOBI)**2 )) )
          enddo
       ENDIF

! Compute the MEAN bed shear stress
       IF(GOT_MEANSHEARSTRESSX) THEN ! first find spint
          spint=0.0
          data_file = ' '
          data_file(1:13) = 'FSoptions.dat'
          INQUIRE(file=data_file, exist=got_fsopts)
          if(got_fsopts) then
             OPEN(556,status='unknown',file=data_file)    
             READ(556,998) rubi
             READ(556,999) rub
             READ(556,999) spint
             CLOSE(556)
          endif
       ENDIF
       IF(GOT_MEANSHEARSTRESSX.AND.GOT_MEANSHEARSTRESSY.AND. &
            GOT_MEANSHEARSTRESSZ ) then
          IF(acctim.lt.spint) N_TIME_AV = 0
          N_TIME_AV = N_TIME_AV + 1  
          do globi=1,nonods     
             IF(acctim.lt.spint) THEN
                MEANSHEARSTRESSX(GLOBI) = 0.
                MEANSHEARSTRESSY(GLOBI) = 0.
                MEANSHEARSTRESSZ(GLOBI) = 0.
             ENDIF

             tempssx = 1023.0*CD*SQRT(U(GLOBI)**2 +V(GLOBI)**2 +W(GLOBI) **2)

             MEANSHEARSTRESSX(GLOBI) = MEANSHEARSTRESSX(GLOBI)* &
                  (N_TIME_AV-1)/N_TIME_AV + tempssx*U(GLOBI)/N_TIME_AV
             MEANSHEARSTRESSY(GLOBI) = MEANSHEARSTRESSY(GLOBI)* &
                  (N_TIME_AV-1)/N_TIME_AV + tempssx*V(GLOBI)/N_TIME_AV
             MEANSHEARSTRESSZ(GLOBI) = MEANSHEARSTRESSZ(GLOBI)* &
                  (N_TIME_AV-1)/N_TIME_AV + tempssx*W(GLOBI)/N_TIME_AV

          enddo
       ELSEIF ( GOT_MEANSHEARSTRESSX ) THEN ! If just got one component then just output the magnitude and not direction.
          IF(acctim.lt.spint) N_TIME_AV = 0
          N_TIME_AV = N_TIME_AV + 1
          do globi=1,nonods     
             IF(acctim.lt.spint) MEANSHEARSTRESSX(GLOBI) = 0.
             MEANSHEARSTRESSX(GLOBI) = MEANSHEARSTRESSX(GLOBI)* &
                  (N_TIME_AV-1)/N_TIME_AV +1023.0*CD*(U(GLOBI)**2 + &
                  V(GLOBI)**2+W(GLOBI)**2)/N_TIME_AV
    ewrite(3,*) ' just left mean bed shear stess vector calculations'
          enddo
       ENDIF
    endif

    ewrite(1,*) 'end of calculate_shear_stress'

  END SUBROUTINE calculate_shear_stress

end module shear_stress
