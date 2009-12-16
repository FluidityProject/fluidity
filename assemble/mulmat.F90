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
!
module mulmat_module
use sml
use fldebug
use AllSorts
implicit none
  
contains

  SUBROUTINE MULMAT( APK, XK, A, FINA, COLA, &
       NCOLA, NBIGM, MATSTR, NONODS, LENGTH, NNODP,SYM)
    !     -------------------------------------------------------------
    ! - APK = BIGM*XK. 
    ! This works for matrices obtained from DIFF2D & DIFF3D as well. 
    ! IF SYM then it is assumed that only the upper diagonal part of 
    ! the matrix is stored. 
    ! NNODP=no of subdomain nodes excluding interface nodes
    ! NNODP=NONODS if no paral
    ! NONODS=no of subdomain nodes plus halo nodes.
    ! LENGTH=NONODS*IDIM
    ! The halo variables of APK contain nothing useful.
    ! MATSTR contains the matrix structure. 
    !
    INTEGER NBIGM
    INTEGER NCOLA, LENGTH, NONODS,NNODP, IDIM,MATSTR
    INTEGER FINA(NONODS+1), COLA(NCOLA)
    INTEGER COUNT,COL, COLPN, COLP2N,COLP3N,COLP4N
    INTEGER CL1,CL2,CL3,CL4
    LOGICAL SYM
    LOGICAL TRANSP,NTRANS,LD3
    LOGICAL BLKSYM,BLKDIA
    ! NB The normal is pointing out of the domain. 
    !  The rotation matrix in 3-D is R=  
    !   T1X   T1Y   T1Z
    !   T2X   T2Y   T2Z
    !   NX    NY    NZ
    !  The rotation matrix in 2-D is R=  
    !   T1X   T1Y   
    !   NX    NY    
    ! WHEN ROTATIONS ARE USED THIS SUB FORMS the vector...
    ! apk=(R^T A R + D) xk in which D contains the b.c's DM1,DM2,DM3.  
    ! 
    REAL    APK(LENGTH), XK(LENGTH), A(NBIGM)
    INTEGER I, IBL22, IBL33, I2, IC, I3, IPN, IP2N, NC
    INTEGER IBL11,IBL21,IBL31,IBL32,IBL41,IBL42,IBL43,IBL44
    INTEGER IBL12,IBL13,IBL14,IBL23,IBL24,IBL34,IP3N,I4
    INTEGER ICOUNT,IL,JL

    IDIM=LENGTH/NONODS

    IF(MATSTR.EQ.6) THEN
       !         print *,'in mulmat idum=',idim
       !         stop 3831
       CALL DG_MAT_SOL_OR_MUL_VEC_MOM(A,FINA,FINA,&
            COLA,NCOLA,NBIGM,&
            IDIM,NONODS,&
            APK,XK,0,.FALSE.)
       RETURN
    ENDIF
    !
    !       ewrite(3,*) 'in mulmat idim,ncola,nbigm,nonods,length,BLKSYM:',
    !     &          idim,ncola,nbigm,nonods,length,BLKSYM
    !
    BLKSYM=.TRUE.
    IF(NBIGM.EQ.IDIM*IDIM*NCOLA) BLKSYM=.FALSE.
    !
    !      ewrite(3,*)'IN MULMAT idim,NBIGM,ncola',idim,NBIGM,ncola
    !
    IF(IDIM.EQ.1) THEN
       IF(SYM) THEN
          do  I = 1, NNODP! Was loop 210
             APK(I)=A(FINA(I))*XK(I)
          end do ! Was loop 210
          do  I = 1, NNODP! Was loop 220
             do  COUNT = FINA(I)+1,FINA(I+1)-1! Was loop 220
                COL=COLA(COUNT)
                APK(I) = APK(I) + A(COUNT)*XK(COL)
                APK(COL) = APK(COL) + A(COUNT)*XK(I)
             end do ! Was loop 220
          end do ! Was loop 220
          do  I = NNODP+1,NONODS! Was loop 211
             APK(I)=0.
          end do ! Was loop 211
       ELSE
          !         ewrite(1,*) 'entering lookp nnodp=',nnodp
          do  I = 1, NNODP! Was loop 200
             APK(I)=0.
             !           ewrite(3,*) 'ncola,i,FINA(I),FINA(I+1)=',
             !     &              ncola,i,FINA(I),FINA(I+1)
             do  COUNT = FINA(I),FINA(I+1)-1! Was loop 200
                APK(I) = APK(I) + A(COUNT)*XK(COLA(COUNT))
             end do ! Was loop 200
          end do ! Was loop 200
          !          ewrite(3,*) 'finished loop'
       ENDIF
    ELSE
       !
       do  I=1,LENGTH! Was loop 5
          APK(I)=0.
       end do ! Was loop 5
       !
       IF(MATSTR.GE.1) THEN
          ! Assume the matrix structure is other way arrounf
          CALL MULOTH(APK, XK, A, FINA, COLA, &
               NCOLA, NBIGM,MATSTR, NONODS, LENGTH, NNODP)
          RETURN
       ENDIF
       !
       !
       BLKDIA=.FALSE.
       IF((IDIM.EQ.2).AND.(NBIGM.LE.2*NCOLA)) BLKDIA=.TRUE.
       IF((IDIM.EQ.3).AND.(NBIGM.LE.3*NCOLA)) BLKDIA=.TRUE.
       !
       !
       IF(BLKDIA) THEN
          IBL22=NCOLA
          IBL33=2*NCOLA
          IF(NBIGM.EQ.NCOLA) THEN
             IBL22=0
             IBL33=0
          ENDIF
          !
          IF(IDIM.EQ.2) THEN
             do  I = 1, NNODP! Was loop 8253
                I2=I+NONODS
                do  IC = FINA(I),FINA(I+1)-1! Was loop 8153
                   CL1=COLA(IC)
                   CL2=CL1+NONODS
                   APK(I) =APK(I) + A(IC)     *XK(CL1)          
                   APK(I2)=APK(I2)+A(IC+IBL22)*XK(CL2)
                end do ! Was loop 8153
             end do ! Was loop 8253
          ENDIF
          IF(IDIM.EQ.3) THEN
             do  I = 1, NNODP! Was loop 8252
                I2=I+NONODS
                I3=I+2*NONODS
                do  IC = FINA(I),FINA(I+1)-1! Was loop 8152
                   CL1=COLA(IC)
                   CL2=CL1+  NONODS
                   CL3=CL1+2*NONODS
                   APK(I) = APK(I) + A(IC)    *XK(CL1)          
                   APK(I2)=APK(I2)+A(IC+IBL22)*XK(CL2)
                   APK(I3)=APK(I3)+A(IC+IBL33)*XK(CL3)
                end do ! Was loop 8152
             end do ! Was loop 8252
          ENDIF
          IF(IDIM.GE.4) THEN
             ewrite(3,*)'INSIDE MULMAT NO OPTION'
             STOP
          ENDIF
          !
          ! ELSE(BLKDIA)
       ELSE
          IF(IDIM.EQ.2) THEN
             !
             IF(BLKSYM) THEN
                do  I = 1, NNODP! Was loop 720
                   !
                   IPN=I+NONODS
                   do  COUNT = FINA(I),FINA(I+1)-1! Was loop 710
                      !
                      COL  =COLA(COUNT)
                      COLPN=COL+NONODS
                      APK(I)  = APK(I)   + A(COUNT)      *XK(COL)
                      APK(COL)= APK(COL) + A(COUNT+NCOLA)*XK(IPN)
                      !
                      APK(IPN) = APK(IPN)+A(COUNT+NCOLA)  *XK(COL)&
                           +A(COUNT+2*NCOLA)*XK(COLPN)
                   end do ! Was loop 710
                   !
                end do ! Was loop 720
             ELSE
                do  I = 1, NNODP! Was loop 7201
                   !
                   IPN=I+NONODS
                   do  COUNT = FINA(I),FINA(I+1)-1! Was loop 7101
                      COL  =COLA(COUNT)
                      COLPN=COL+NONODS
                      APK(I) = APK(I) +A(COUNT)      *XK(COL) &
                           +A(COUNT+NCOLA)*XK(COLPN)
                      !
                      APK(IPN)=APK(IPN)+A(COUNT+2*NCOLA)*XK(COL)&
                           +A(COUNT+3*NCOLA)*XK(COLPN)
                   end do ! Was loop 7101
                   !
                end do ! Was loop 7201
             ENDIF
          ENDIF
          ! ENDOF IF IDIM=2
          !
          IF(IDIM.EQ.3) THEN
             !
             IF(BLKSYM) THEN
                do  I = 1, NNODP! Was loop 825
                   IPN =I+  NONODS
                   IP2N=I+2*NONODS
                   do  COUNT = FINA(I),FINA(I+1)-1! Was loop 815
                      !
                      COL   =COLA(COUNT)
                      COLPN =COL+  NONODS
                      COLP2N=COL+2*NONODS
                      !
                      APK(COL) = APK(COL) + A(COUNT+NCOLA)  *XK(IPN)&
                           + A(COUNT+3*NCOLA)*XK(IP2N)
                      !
                      APK(COLPN) = APK(COLPN) + A(COUNT+4*NCOLA)*XK(IP2N)
                      !
                      APK(I) = APK(I) + A(COUNT)*XK(COL)
                      !
                      APK(IPN) = APK(IPN) + A(COUNT+NCOLA)  *XK(COL) &
                           + A(COUNT+2*NCOLA)*XK(COLPN)
                      !
                      APK(IP2N) = APK(IP2N)+ A(COUNT+3*NCOLA)*XK(COL) &
                           + A(COUNT+4*NCOLA)*XK(COLPN) &
                           + A(COUNT+5*NCOLA)*XK(COLP2N)
                      !
                   end do ! Was loop 815
                end do ! Was loop 825
             ELSE
                NC=NCOLA
                do  I = 1, NNODP! Was loop 8251
                   I2=I+  NONODS
                   I3=I+2*NONODS
                   do  IC = FINA(I),FINA(I+1)-1! Was loop 8151
                      !
                      CL1=COLA(IC)
                      CL2=CL1+  NONODS
                      CL3=CL1+2*NONODS
                      !
                      APK(I) = APK(I)+A(IC)     *XK(CL1)+A(IC+NC)  *XK(CL2)&
                           +A(IC+2*NC)*XK(CL3)
                      APK(I2)=APK(I2)+A(IC+3*NC)*XK(CL1)+A(IC+4*NC)*XK(CL2)&
                           +A(IC+5*NC)*XK(CL3)
                      APK(I3)=APK(I3)+A(IC+6*NC)*XK(CL1)+A(IC+7*NC)*XK(CL2)&
                           +A(IC+8*NC)*XK(CL3)
                      !
                   end do ! Was loop 8151
                end do ! Was loop 8251
             ENDIF
             !
             ! END OF IF IDIM=3
          ENDIF
          !
          IF(IDIM.EQ.4) THEN
             IF(BLKSYM) THEN
                IBL11=0
                IBL21=NCOLA
                IBL22=2*NCOLA
                IBL31=3*NCOLA
                IBL32=4*NCOLA
                IBL33=5*NCOLA
                IBL41=6*NCOLA
                IBL42=7*NCOLA
                IBL43=8*NCOLA
                IBL44=9*NCOLA
                !
                IBL12=NCOLA
                IBL13=3*NCOLA
                IBL14=6*NCOLA
                IBL23=4*NCOLA
                IBL24=7*NCOLA
                IBL34=8*NCOLA
                do  I = 1, NNODP! Was loop 7875
                   IPN =I+  NONODS
                   IP2N=I+2*NONODS
                   IP3N=I+3*NONODS
                   do  COUNT = FINA(I),FINA(I+1)-1! Was loop 2875
                      !
                      COL   =COLA(COUNT)
                      COLPN =COL+  NONODS
                      COLP2N=COL+2*NONODS
                      COLP3N=COL+3*NONODS
                      !
                      APK(COL) = APK(COL) + A(COUNT+IBL12)*XK(IPN)&
                           + A(COUNT+IBL13)*XK(IP2N)&
                           + A(COUNT+IBL14)*XK(IP3N)
                      !
                      APK(COLPN) = APK(COLPN) + A(COUNT+IBL23)*XK(IP2N)&
                           + A(COUNT+IBL24)*XK(IP3N)
                      !
                      APK(COLP2N)= APK(COLP2N)+ A(COUNT+IBL34)*XK(IP3N)
                      !
                      APK(I) = APK(I) + A(COUNT+IBL11)*XK(COL)
                      !
                      APK(IPN) = APK(IPN) + A(COUNT+IBL21)  *XK(COL) &
                           + A(COUNT+IBL22)*XK(COLPN)
                      !
                      APK(IP2N) = APK(IP2N)+ A(COUNT+IBL31)*XK(COL) &
                           + A(COUNT+IBL32)*XK(COLPN) &
                           + A(COUNT+IBL33)*XK(COLP2N)
                      !
                      APK(IP3N) = APK(IP3N)+ A(COUNT+IBL41)*XK(COL) &
                           + A(COUNT+IBL42)*XK(COLPN) &
                           + A(COUNT+IBL43)*XK(COLP2N)&
                           + A(COUNT+IBL44)*XK(COLP3N)
                      !
                   end do ! Was loop 2875
                end do ! Was loop 7875
             ELSE
                NC=NCOLA
                IBL11=0
                IBL12=NCOLA
                IBL13=2*NCOLA
                IBL14=3*NCOLA
                !
                IBL21=4*NCOLA
                IBL22=5*NCOLA
                IBL23=6*NCOLA
                IBL24=7*NCOLA
                !
                IBL31=8*NCOLA
                IBL32=9*NCOLA
                IBL33=10*NCOLA
                IBL34=11*NCOLA
                !
                IBL41=12*NCOLA
                IBL42=13*NCOLA
                IBL43=14*NCOLA
                IBL44=15*NCOLA
                !
                !        ewrite(3,*) 'IBL11,A(IBL11+1):',IBL11,A(IBL11+1)
                !        ewrite(3,*) 'IBL12,A(IBL12+1):',IBL12,A(IBL12+1)
                !        ewrite(3,*) 'IBL13,A(IBL13+1):',IBL13,A(IBL13+1)
                !        ewrite(3,*) 'IBL14,A(IBL14+1):',IBL14,A(IBL14+1)
                !        ewrite(3,*) ' '
                !
                !        ewrite(3,*) 'IBL21,A(IBL21+1):',IBL21,A(IBL21+1)
                !        ewrite(3,*) 'IBL22,A(IBL22+1):',IBL22,A(IBL22+1)
                !        ewrite(3,*) 'IBL23,A(IBL23+1):',IBL23,A(IBL23+1)
                !        ewrite(3,*) 'IBL24,A(IBL24+1):',IBL24,A(IBL24+1)
                !        ewrite(3,*) ' '
                !
                !        ewrite(3,*) 'IBL31,A(IBL31+1):',IBL31,A(IBL31+1)
                !        ewrite(3,*) 'IBL32,A(IBL32+1):',IBL32,A(IBL32+1)
                !        ewrite(3,*) 'IBL33,A(IBL33+1):',IBL33,A(IBL33+1)
                !        ewrite(3,*) 'IBL34,A(IBL34+1):',IBL34,A(IBL34+1)
                !        ewrite(3,*) ' '
                !
                !        ewrite(3,*) 'IBL41,A(IBL41+1):',IBL41,A(IBL41+1)
                !        ewrite(3,*) 'IBL42,A(IBL42+1):',IBL42,A(IBL42+1)
                !        ewrite(3,*) 'IBL43,A(IBL43+1):',IBL43,A(IBL43+1)
                !        ewrite(3,*) 'IBL44,A(IBL44+1):',IBL44,A(IBL44+1)
                !        ewrite(3,*) ' '
                !
                !       ewrite(3,*) 'idim,ncola,nbigm,nonods,length,BLKSYM:',
                !     &          idim,ncola,nbigm,nonods,length,BLKSYM
                do  I = 1, NNODP! Was loop 2182
                   I2=I+  NONODS
                   I3=I+2*NONODS
                   I4=I+3*NONODS
                   do  COUNT = FINA(I),FINA(I+1)-1! Was loop 2181
                      !
                      CL1=COLA(COUNT)
                      CL2=CL1+  NONODS
                      CL3=CL1+2*NONODS
                      CL4=CL1+3*NONODS
                      !
                      APK(I) = APK(I)   + A(COUNT+IBL11)*XK(CL1)&
                           + A(COUNT+IBL12)*XK(CL2)&
                           + A(COUNT+IBL13)*XK(CL3)&
                           + A(COUNT+IBL14)*XK(CL4)
                      !
                      APK(I2) = APK(I2) + A(COUNT+IBL21)*XK(CL1) &
                           + A(COUNT+IBL22)*XK(CL2)&
                           + A(COUNT+IBL23)*XK(CL3)&
                           + A(COUNT+IBL24)*XK(CL4)
                      !
                      APK(I3) = APK(I3) + A(COUNT+IBL31)*XK(CL1) &
                           + A(COUNT+IBL32)*XK(CL2) &
                           + A(COUNT+IBL33)*XK(CL3)&
                           + A(COUNT+IBL34)*XK(CL4)
                      !
                      APK(I4) = APK(I4) + A(COUNT+IBL41)*XK(CL1) &
                           + A(COUNT+IBL42)*XK(CL2) &
                           + A(COUNT+IBL43)*XK(CL3)&
                           + A(COUNT+IBL44)*XK(CL4)
                      !
                   end do ! Was loop 2181
                end do ! Was loop 2182
             ENDIF
          ENDIF
          !
          IF(IDIM.GE.5) THEN
             IF(BLKSYM) THEN
                do  I = 1, NNODP! Was loop 2825
                   do  COUNT = FINA(I),FINA(I+1)-1! Was loop 2815
                      !
                      CL1=COLA(COUNT)
                      ICOUNT=0
                      do  IL=1,IDIM! Was loop 143
                         do  JL=1,IL-1! Was loop 144
                            !
                            APK(I+(IL-1)*NONODS) = APK(I+(IL-1)*NONODS)   &
                                 + A(COUNT+NCOLA*ICOUNT)*XK(CL1+(JL-1)*NONODS)
                            ! The upper part ...
                            APK(CL1+(JL-1)*NONODS) = APK(CL1+(JL-1)*NONODS)   &
                                 + A(COUNT+NCOLA*ICOUNT)*XK(I+(IL-1)*NONODS)
                            !
                            ICOUNT=ICOUNT+1
                            !
                         end do ! Was loop 144
                         !            JL=IL
                         APK(I+(IL-1)*NONODS) = APK(I+(IL-1)*NONODS)   &
                              + A(COUNT+NCOLA*ICOUNT)*XK(CL1+(IL-1)*NONODS)
                         ICOUNT=ICOUNT+1
                      end do ! Was loop 143
                      !
                   end do ! Was loop 2815
                end do ! Was loop 2825
             ELSE
                do  I = 1, NNODP! Was loop 2782
                   do  COUNT = FINA(I),FINA(I+1)-1! Was loop 7181
                      !
                      CL1=COLA(COUNT)
                      do  IL=1,IDIM! Was loop 133
                         do  JL=1,IDIM! Was loop 133
                            !
                            APK(I+(IL-1)*NONODS) = APK(I+(IL-1)*NONODS)   &
                                 + A(COUNT+NCOLA*( (IL-1)*IDIM+ JL-1))*XK(CL1+(JL-1)*NONODS)
                         end do ! Was loop 133
                      end do ! Was loop 133
                      !
                   end do ! Was loop 7181
                end do ! Was loop 2782
             ENDIF
             ! END OF IF (IDIM.GE.5) THEN
          ENDIF
          !
          !
          ! ENDIF ELSE (BLKDIA)
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE MULMAT
  !
  !
  !
  !
  SUBROUTINE MULOTH(APK, XK, A, FINA, COLA, &
       NCOLA, NBIGM, MATSTR, NONODS, LENGTH, NNODP) 
    IMPLICIT NONE
    ! Assume the matrix structure is other way arround. 
    ! I.E. moments first then nodes. 
    ! - APK = BIGM*XK. 
    ! This works for matrices obtained from DIFF2D & DIFF3D as well. 
    ! IF SYM then it is assumed that only the upper diagonal part of 
    ! the matrix is stored. 
    ! NNODP=no of subdomain nodes excluding interface nodes
    ! NNODP=NONODS if no paral
    ! NONODS=no of subdomain nodes plus halo nodes.
    ! LENGTH=NONODS*IDIM
    ! The halo variables of APK contain nothing useful.
    ! MATSTR=matrix stucture. 
    !
    INTEGER NCOLA, LENGTH, NONODS,NNODP, IDIM, IDIM2,MATSTR, NBIGM
    INTEGER FINA(NONODS+1), COLA(NCOLA)
    INTEGER COUNT,COL, COLPN, COLP2N,COLP3N,COLP4N
    INTEGER CL1,CL2,CL3,CL4
    LOGICAL BLKSYM
    ! The block structure of the matrix is:
    !   1 * 
    !   2 3   *-block not stored -for 2-D
    ! or 
    !  1 2
    !  3 4  -for 2-D 
    ! for 3-D ...
    !   1 * *
    !   2 3 *
    !   4 5 6    *-block not stored -for 3-D
    ! OR 
    !   1 2 3
    !   4 5 6
    !   7 8 9    -for 3-D
    !
    REAL    APK(LENGTH), XK(LENGTH), A(NBIGM)
    INTEGER I, ICOUNT, IL, JL


    IDIM=LENGTH/NONODS
    IDIM2=IDIM**2
    BLKSYM=.FALSE.
    IF(NBIGM.NE.IDIM2*NCOLA) BLKSYM=.TRUE.
    !      ewrite(3,*) 'IN MULOTH:idim,matstr,sym:',
    !     :              idim,matstr,blksym
    IF(BLKSYM) THEN
       do  I = 1, NNODP! Was loop 2825
          do  COUNT = FINA(I),FINA(I+1)-1! Was loop 2815
             !
             CL1=COLA(COUNT)
             ICOUNT=1
             do  IL=1,IDIM! Was loop 143
                do  JL=1,IL-1! Was loop 144
                   !
                   !          APK(I+(IL-1)*NONODS) = APK(I+(IL-1)*NONODS)   
                   !     &      + A(COUNT+NCOLA*ICOUNT)*XK(CL1+(JL-1)*NONODS)
                   APK( (I-1)*IDIM + IL ) = APK((I-1)*IDIM + IL)   &
                        + A( (COUNT-1)*IDIM2+ICOUNT)*XK( (CL1-1)*IDIM+JL)
                   ! The upper part ...
                   !          APK(CL1+(JL-1)*NONODS) = APK(CL1+(JL-1)*NONODS)   
                   !     &      + A(COUNT+NCOLA*ICOUNT)*XK(I+(IL-1)*NONODS)
                   APK( (CL1-1)*IDIM+JL) = APK( (CL1-1)*IDIM+JL)   &
                        + A( (COUNT-1)*IDIM2+ICOUNT)*XK((I-1)*IDIM+IL)
                   !
                   ICOUNT=ICOUNT+1
                   !
                end do ! Was loop 144
                !            JL=IL
                !           APK(I+(IL-1)*NONODS) = APK(I+(IL-1)*NONODS)   
                !     &      + A(COUNT+NCOLA*ICOUNT)*XK(CL1+(IL-1)*NONODS)
                APK( (I-1)*IDIM+IL) = APK((I-1)*IDIM+IL)   &
                     + A( (COUNT-1)*IDIM2+ICOUNT)*XK( (CL1-1)*IDIM+IL)
                ICOUNT=ICOUNT+1
             end do ! Was loop 143
             !
          end do ! Was loop 2815
       end do ! Was loop 2825
    ELSE
       IF(MATSTR.EQ.2) THEN
          ! Adrians ordering...
          do  I = 1, NNODP! Was loop 2782
             do  COUNT = FINA(I),FINA(I+1)-1! Was loop 7181
                !
                CL1=COLA(COUNT)
                do  IL=1,IDIM! Was loop 133
                   do  JL=1,IDIM! Was loop 133
                      !
                      !          APK(I+(IL-1)*NONODS) = APK(I+(IL-1)*NONODS)   
                      !     &      + A(COUNT+NCOLA*( (IL-1)*IDIM+ JL-1))*XK(CL1+(JL-1)*NONODS)
                      !
                      APK((I-1)*IDIM+IL) = APK((I-1)*IDIM+IL)   &
                           + A( (COUNT-1)*IDIM2 + (JL-1)*IDIM+ IL)&
                           *XK( (CL1-1)*IDIM+JL)
                      !
                   end do ! Was loop 133
                end do ! Was loop 133
                !
             end do ! Was loop 7181
          end do ! Was loop 2782
          !
       ELSE
          do  I = 1, NNODP! Was loop 3782
             do  COUNT = FINA(I),FINA(I+1)-1! Was loop 8181
                !
                CL1=COLA(COUNT)
                do  IL=1,IDIM! Was loop 333
                   do  JL=1,IDIM! Was loop 333
                      !
                      !          APK(I+(IL-1)*NONODS) = APK(I+(IL-1)*NONODS)   
                      !     &      + A(COUNT+NCOLA*( (IL-1)*IDIM+ JL-1))*XK(CL1+(JL-1)*NONODS)
                      !
                      APK((I-1)*IDIM+IL) = APK((I-1)*IDIM+IL)   &
                           + A( (COUNT-1)*IDIM2 + (IL-1)*IDIM+ JL)&
                           *XK( (CL1-1)*IDIM+JL)
                      !           if((I-1)*IDIM+IL.eq.4) print *,'a,XK( (CL1-1)*IDIM+JL):',
                      !     &      A( (COUNT-1)*IDIM2 + (IL-1)*IDIM+ JL),XK( (CL1-1)*IDIM+JL)
                      !
                   end do ! Was loop 333
                end do ! Was loop 333
                !
             end do ! Was loop 8181
          end do ! Was loop 3782
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE MULOTH
  !
  !
  !
  !       
  SUBROUTINE DG_MAT_SOL_OR_MUL_VEC_MOM(AMATELE,CENTRMELE,FINDRMELE,&
       COLMELE,NCOLMELE,NBIGMELE,&
       BIG_NLOC,TOTELE,&
       B,VEC,IMM,SOLVE)
    ! matrix vector multiplication for B= (R (AMATELE +C)R^T+D)  VEC
    ! IF SOLVE then solve A B=VEC for b
    IMPLICIT NONE
    INTEGER TOTELE
    INTEGER BIG_NLOC
    INTEGER CENTRMELE(TOTELE),FINDRMELE(TOTELE+1)
    INTEGER NBIGMELE,NCOLMELE,COLMELE(NCOLMELE)
    REAL AMATELE(NBIGMELE)
    REAL B(TOTELE*BIG_NLOC)
    REAL VEC(TOTELE*BIG_NLOC)
    LOGICAL SOLVE
    INTEGER IMM
    ! IMM is preconditioner.
    ! Corriolis...
    ! place into rotations matrices:
    ! DGT1X,DGT1Y,DGT1Z, 
    ! DGT2X,DGT2Y,DGT2Z, 
    ! DGNX,DGNY,DGNZ,
    ! and b.c's               DGDM1,DGDM2,DGDM3
    ! Local variables...
    ! Integer pointers...
    INTEGER COMAT11,COMAT12,COMAT13, &
         COMAT21,COMAT22,COMAT23, &
         COMAT31,COMAT32,COMAT33,  &
         DGT1X,DGT1Y,DGT1Z, &
         DGT2X,DGT2Y,DGT2Z, &
         DGNX,DGNY,DGNZ,&
         DGDM1,DGDM2,DGDM3
    INTEGER ELE,ILOC,NLOC,RPT
    REAL, ALLOCATABLE, DIMENSION(:)::BX
    REAL, ALLOCATABLE, DIMENSION(:)::BY
    REAL, ALLOCATABLE, DIMENSION(:)::BZ
    REAL, ALLOCATABLE, DIMENSION(:)::VECX
    REAL, ALLOCATABLE, DIMENSION(:)::VECY
    REAL, ALLOCATABLE, DIMENSION(:)::VECZ

    NLOC=4

    ALLOCATE(BX(TOTELE*NLOC))
    ALLOCATE(BY(TOTELE*NLOC))
    ALLOCATE(BZ(TOTELE*NLOC))

    ALLOCATE(VECX(TOTELE*NLOC))
    ALLOCATE(VECY(TOTELE*NLOC))
    ALLOCATE(VECZ(TOTELE*NLOC))

    do ELE=1,TOTELE
       do ILOC=1,NLOC
          VECX((ELE-1)*NLOC+ILOC)=VEC((ELE-1)*BIG_NLOC+ILOC)
          VECY((ELE-1)*NLOC+ILOC)=VEC((ELE-1)*BIG_NLOC+ILOC+NLOC)
          VECZ((ELE-1)*NLOC+ILOC)=VEC((ELE-1)*BIG_NLOC+ILOC+2*NLOC)
       END DO
    END DO

    RPT=NBIGMELE-21*TOTELE*NLOC+1
    CALL ALOMEM(COMAT11,   RPT, TOTELE*NLOC, NBIGMELE)
    CALL ALOMEM(COMAT12,   RPT, TOTELE*NLOC, NBIGMELE)
    CALL ALOMEM(COMAT13,   RPT, TOTELE*NLOC, NBIGMELE)

    CALL ALOMEM(COMAT21,   RPT, TOTELE*NLOC, NBIGMELE)
    CALL ALOMEM(COMAT22,   RPT, TOTELE*NLOC, NBIGMELE)
    CALL ALOMEM(COMAT23,   RPT, TOTELE*NLOC, NBIGMELE)

    CALL ALOMEM(COMAT31,   RPT, TOTELE*NLOC, NBIGMELE)
    CALL ALOMEM(COMAT32,   RPT, TOTELE*NLOC, NBIGMELE)
    CALL ALOMEM(COMAT33,   RPT, TOTELE*NLOC, NBIGMELE)

    CALL ALOMEM(DGT1X,   RPT, TOTELE*NLOC, NBIGMELE)
    CALL ALOMEM(DGT1Y,   RPT, TOTELE*NLOC, NBIGMELE)
    CALL ALOMEM(DGT1Z,   RPT, TOTELE*NLOC, NBIGMELE)

    CALL ALOMEM(DGT2X,   RPT, TOTELE*NLOC, NBIGMELE)
    CALL ALOMEM(DGT2Y,   RPT, TOTELE*NLOC, NBIGMELE)
    CALL ALOMEM(DGT2Z,   RPT, TOTELE*NLOC, NBIGMELE)

    CALL ALOMEM(DGNX,    RPT, TOTELE*NLOC, NBIGMELE)
    CALL ALOMEM(DGNY,    RPT, TOTELE*NLOC, NBIGMELE)
    CALL ALOMEM(DGNZ,    RPT, TOTELE*NLOC, NBIGMELE)

    CALL ALOMEM(DGDM1,   RPT, TOTELE*NLOC, NBIGMELE)
    CALL ALOMEM(DGDM2,   RPT, TOTELE*NLOC, NBIGMELE)
    CALL ALOMEM(DGDM3,   RPT, TOTELE*NLOC, NBIGMELE)

    IF(SOLVE) THEN
       CALL DG_SOLVE_PRE_MOM(AMATELE,CENTRMELE,FINDRMELE,&
            COLMELE,NCOLMELE,NBIGMELE,&
            NLOC,BIG_NLOC,TOTELE,&
            BX,BY,BZ,&
            VECX,VECY,VECZ,&
            AMATELE(COMAT11),AMATELE(COMAT12),AMATELE(COMAT13),AMATELE(COMAT21),&
            AMATELE(COMAT22),AMATELE(COMAT23),AMATELE(COMAT31),AMATELE(COMAT32),AMATELE(COMAT33),& 
            AMATELE(DGT1X),AMATELE(DGT1Y),AMATELE(DGT1Z), &
            AMATELE(DGT2X),AMATELE(DGT2Y),AMATELE(DGT2Z), &
            AMATELE(DGNX),AMATELE(DGNY),AMATELE(DGNZ),&
            AMATELE(DGDM1),AMATELE(DGDM2),AMATELE(DGDM3),IMM)
    ELSE
       !         
       CALL DG_MAT_VEC_MOM(AMATELE,CENTRMELE,FINDRMELE,&
            COLMELE,NCOLMELE,NBIGMELE,&
            NLOC,BIG_NLOC,TOTELE,&
            BX,BY,BZ,&
            VECX,VECY,VECZ,&
            AMATELE(COMAT11),AMATELE(COMAT12),AMATELE(COMAT13),AMATELE(COMAT21),&
            AMATELE(COMAT22),AMATELE(COMAT23),AMATELE(COMAT31),AMATELE(COMAT32),AMATELE(COMAT33),&
            AMATELE(DGT1X),AMATELE(DGT1Y),AMATELE(DGT1Z), &
            AMATELE(DGT2X),AMATELE(DGT2Y),AMATELE(DGT2Z), &
            AMATELE(DGNX),AMATELE(DGNY),AMATELE(DGNZ),&
            AMATELE(DGDM1),AMATELE(DGDM2),AMATELE(DGDM3),IMM)
    ENDIF

    do ELE=1,TOTELE
       do ILOC=1,NLOC
          B((ELE-1)*BIG_NLOC+ILOC)       =BX((ELE-1)*NLOC+ILOC)
          B((ELE-1)*BIG_NLOC+ILOC+NLOC)  =BY((ELE-1)*NLOC+ILOC)
          B((ELE-1)*BIG_NLOC+ILOC+2*NLOC)=BZ((ELE-1)*NLOC+ILOC)
       END DO
    END DO

    RETURN
  END SUBROUTINE DG_MAT_SOL_OR_MUL_VEC_MOM
  !
  !
  !
  !
  SUBROUTINE DG_MAT_VEC_MOM(AMATELE,CENTRMELE,FINDRMELE,&
       COLMELE,NCOLMELE,NBIGMELE,&
       NLOC,BIG_NLOC,TOTELE,&
       BX,BY,BZ,&
       VECX,VECY,VECZ,&
       COMAT11,COMAT12,COMAT13,COMAT21,&
       COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
       DGT1X,DGT1Y,DGT1Z, &
       DGT2X,DGT2Y,DGT2Z, &
       DGNX,DGNY,DGNZ,&
       DGDM1,DGDM2,DGDM3,IMM)
    ! matrix vector multiplication for B= (R (AMATELE +C) R^T +D)  VEC
    IMPLICIT NONE
    LOGICAL NO_ROTS
    PARAMETER(NO_ROTS=.TRUE.)
    INTEGER NLOC,TOTELE
    INTEGER BIG_NLOC
    INTEGER CENTRMELE(TOTELE),FINDRMELE(TOTELE+1)
    INTEGER NCOLMELE,NBIGMELE,COLMELE(NCOLMELE)
    REAL AMATELE(NBIGMELE)
    REAL BX(TOTELE*NLOC),BY(TOTELE*NLOC),BZ(TOTELE*NLOC)
    REAL VECX(TOTELE*NLOC),VECY(TOTELE*NLOC),VECZ(TOTELE*NLOC)
    ! Corriolis...
    REAL COMAT11(TOTELE,NLOC),COMAT12(TOTELE,NLOC),COMAT13(TOTELE,NLOC)
    REAL COMAT21(TOTELE,NLOC),COMAT22(TOTELE,NLOC),COMAT23(TOTELE,NLOC)
    REAL COMAT31(TOTELE,NLOC),COMAT32(TOTELE,NLOC),COMAT33(TOTELE,NLOC)
    ! place into rotations matrices:
    ! DGT1X,DGT1Y,DGT1Z, 
    ! DGT2X,DGT2Y,DGT2Z, 
    ! DGNX,DGNY,DGNZ,
    ! and b.c's               DGDM1,DGDM2,DGDM3
    REAL DGT1X(TOTELE*NLOC),DGT1Y(TOTELE*NLOC),DGT1Z(TOTELE*NLOC)
    REAL DGT2X(TOTELE*NLOC),DGT2Y(TOTELE*NLOC),DGT2Z(TOTELE*NLOC)
    REAL DGNX(TOTELE*NLOC),DGNY(TOTELE*NLOC),DGNZ(TOTELE*NLOC)
    REAL DGDM1(TOTELE*NLOC),DGDM2(TOTELE*NLOC),DGDM3(TOTELE*NLOC)
    INTEGER IMM
    ! Local variables...
    INTEGER ILOC,JLOC,ILOC3,JLOC3,ELE,ELE2
    INTEGER COUNT,IDG,JDG,GLOBI,NOD,ID,JD,II
    INTEGER BIG_ILOC,BIG_JLOC,IJDISP,NBLOCK
    REAL RAMAT
    REAL, ALLOCATABLE, DIMENSION(:)::RT_VECX
    REAL, ALLOCATABLE, DIMENSION(:)::RT_VECY
    REAL, ALLOCATABLE, DIMENSION(:)::RT_VECZ

    ALLOCATE(RT_VECX(TOTELE*NLOC))
    ALLOCATE(RT_VECY(TOTELE*NLOC))
    ALLOCATE(RT_VECZ(TOTELE*NLOC))
    RT_VECX=VECX
    RT_VECY=VECY
    RT_VECZ=VECZ

    ! Rotate RT_VEC=RT*VEC...
    if(.NOT.NO_ROTS) then
       CALL UVW_ROT_DG_ALWAYS(NLOC,BIG_NLOC,TOTELE,&
            .TRUE.,RT_VECX,RT_VECY,RT_VECZ,&
            DGT1X,DGT1Y,DGT1Z, &
            DGT2X,DGT2Y,DGT2Z, &
            DGNX,DGNY,DGNZ)
    endif
    BX=0.0
    BY=0.0
    BZ=0.0

    do ELE=1,TOTELE
       !         if(ele.eq.1) print *,'FINDRMELE(ELE),FINDRMELE(ELE+1)-1:',
       !     &                         FINDRMELE(ELE),FINDRMELE(ELE+1)-1
       do COUNT=FINDRMELE(ELE),FINDRMELE(ELE+1)-1
          ELE2=COLMELE(COUNT)
          ! Multiply block by R block R^T
          do ILOC=1,NLOC
             IDG=(ELE-1)*NLOC+ILOC
             do JLOC=1,NLOC
                JDG=(ELE2-1)*NLOC+JLOC
                IJDISP=(ILOC-1)*NLOC + JLOC
                RAMAT=AMATELE((COUNT-1)*NLOC*NLOC+IJDISP)
                BX(IDG)=BX(IDG)+RAMAT*RT_VECX(JDG)
                BY(IDG)=BY(IDG)+RAMAT*RT_VECY(JDG)
                BZ(IDG)=BZ(IDG)+RAMAT*RT_VECZ(JDG)
                !        if(idg.eq.4) print *,'RAMAT,RT_VECX(JDG):',RAMAT,RT_VECX(JDG)
             END DO
          END DO
       END DO
    END DO

    ! CORRIOLIS TERM...         
    do ELE=1,TOTELE
       ! Multiply block by R block R^T
       do ILOC=1,NLOC
          IDG=(ELE-1)*NLOC+ILOC
          BX(IDG)=BX(IDG) +COMAT11(ELE,ILOC)*RT_VECX(IDG)+COMAT12(ELE,ILOC)*RT_VECY(IDG) &
               +COMAT13(ELE,ILOC)*RT_VECZ(IDG)
          BY(IDG)=BY(IDG) +COMAT21(ELE,ILOC)*RT_VECX(IDG)+COMAT22(ELE,ILOC)*RT_VECY(IDG) &
               +COMAT23(ELE,ILOC)*RT_VECZ(IDG)
          BZ(IDG)=BZ(IDG) +COMAT31(ELE,ILOC)*RT_VECX(IDG)+COMAT32(ELE,ILOC)*RT_VECY(IDG) &
               +COMAT33(ELE,ILOC)*RT_VECZ(IDG)
       END DO
    END DO
    ! Rotate B=R B...
    !       if(.false.) then
    if(.NOT.NO_ROTS) then
       CALL UVW_ROT_DG_ALWAYS(NLOC,BIG_NLOC,TOTELE,&
            .FALSE.,BX,BY,BZ,&
            DGT1X,DGT1Y,DGT1Z, &
            DGT2X,DGT2Y,DGT2Z, &
            DGNX,DGNY,DGNZ)
    endif
    ! B=B+D*VEC
    do ELE=1,TOTELE
       do ILOC=1,NLOC
          IDG=(ELE-1)*NLOC+ILOC
          BX(IDG)=BX(IDG) +DGDM1(IDG)*VECX(IDG)
          BY(IDG)=BY(IDG) +DGDM2(IDG)*VECY(IDG)
          BZ(IDG)=BZ(IDG) +DGDM3(IDG)*VECZ(IDG)
       END DO
    END DO
    RETURN
  END SUBROUTINE DG_MAT_VEC_MOM
  !
  !
  !
  !
  SUBROUTINE DG_SOLVE_PRE_MOM(AMATELE,CENTRMELE,FINDRMELE,&
       COLMELE,NCOLMELE,NBIGMELE,&
       NLOC,BIG_NLOC,TOTELE,&
       BX,BY,BZ,&
       VECX,VECY,VECZ,&
       COMAT11,COMAT12,COMAT13,COMAT21,&
       COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
       DGT1X,DGT1Y,DGT1Z, &
       DGT2X,DGT2Y,DGT2Z, &
       DGNX,DGNY,DGNZ,&
       DGDM1,DGDM2,DGDM3,IMM)
    ! solve A B=VEC for b using 1 step BFBGS & A= (R (AMATELE +C)R^T +D)
    ! STARTING FROM 0. 
    INTEGER NLOC,TOTELE
    INTEGER BIG_NLOC
    INTEGER CENTRMELE(TOTELE),FINDRMELE(TOTELE+1)
    INTEGER NCOLMELE,NBIGMELE,COLMELE(NCOLMELE)
    REAL AMATELE(NBIGMELE)
    REAL BX(TOTELE*NLOC),BY(TOTELE*NLOC),BZ(TOTELE*NLOC)
    REAL VECX(TOTELE*NLOC),VECY(TOTELE*NLOC),VECZ(TOTELE*NLOC)
    INTEGER IMM
    ! IMM=preconditioning option. 
    ! Corriolis...
    REAL COMAT11(TOTELE,NLOC),COMAT12(TOTELE,NLOC),COMAT13(TOTELE,NLOC)
    REAL COMAT21(TOTELE,NLOC),COMAT22(TOTELE,NLOC),COMAT23(TOTELE,NLOC)
    REAL COMAT31(TOTELE,NLOC),COMAT32(TOTELE,NLOC),COMAT33(TOTELE,NLOC)
    ! place into rotations matrices:
    ! DGT1X,DGT1Y,DGT1Z, 
    ! DGT2X,DGT2Y,DGT2Z, 
    ! DGNX,DGNY,DGNZ,
    ! and b.c's               DGDM1,DGDM2,DGDM3
    REAL DGT1X(TOTELE*NLOC),DGT1Y(TOTELE*NLOC),DGT1Z(TOTELE*NLOC)
    REAL DGT2X(TOTELE*NLOC),DGT2Y(TOTELE*NLOC),DGT2Z(TOTELE*NLOC)
    REAL DGNX(TOTELE*NLOC),DGNY(TOTELE*NLOC),DGNZ(TOTELE*NLOC)
    REAL DGDM1(TOTELE*NLOC),DGDM2(TOTELE*NLOC),DGDM3(TOTELE*NLOC)
    ! Local variables...
    INTEGER ILOC,JLOC,ILOC3,JLOC3,ELE,ELE2
    INTEGER COUNT,IDG,JDG,GLOBI,NOD,ID,JD,I,J,ILOOP,NILOOP
    INTEGER BIG_ILOC,BIG_JLOC,IJDISP
    INTEGER ELESTART,ELEFINI,ELESTEP
    INTEGER COUNTSTART,COUNTFINI
    REAL RAMAT
    REAL RMATBLOKRT(BIG_NLOC,BIG_NLOC),MATBLOKRT(BIG_NLOC,BIG_NLOC)
    REAL MATBLOK(BIG_NLOC,BIG_NLOC),LOCVEC(BIG_NLOC)
    REAL RLOCWORK(BIG_NLOC),LOCWORK(BIG_NLOC)
    REAL RLOC(BIG_NLOC,BIG_NLOC),RTLOC(BIG_NLOC,BIG_NLOC)
    REAL LOCB(BIG_NLOC),RT_LOCB(BIG_NLOC)
    REAL, ALLOCATABLE, DIMENSION(:)::RT_BX
    REAL, ALLOCATABLE, DIMENSION(:)::RT_BY
    REAL, ALLOCATABLE, DIMENSION(:)::RT_BZ
    REAL, ALLOCATABLE, DIMENSION(:)::WORKX
    REAL, ALLOCATABLE, DIMENSION(:)::WORKY
    REAL, ALLOCATABLE, DIMENSION(:)::WORKZ

    if(IMM.EQ.0) then
       do ELE=1,TOTELE
          COUNT=CENTRMELE(ELE)
          do ILOC=1,NLOC
             IDG=(ELE-1)*NLOC+ILOC
             JLOC=ILOC
             IJDISP=(ILOC-1)*NLOC + JLOC
             RAMAT=AMATELE((COUNT-1)*NLOC*NLOC+IJDISP)
             BX(IDG)=VECX(IDG)/(RAMAT+DGDM1(IDG)+COMAT11(ELE,ILOC))
             BY(IDG)=VECY(IDG)/(RAMAT+DGDM2(IDG)+COMAT22(ELE,ILOC))
             BZ(IDG)=VECZ(IDG)/(RAMAT+DGDM3(IDG)+COMAT33(ELE,ILOC))
             !               print *,'(RAMAT+DGDM1(IDG)):',(RAMAT+DGDM1(IDG))
             !               print *,'(RAMAT+DGDM2(IDG)):',(RAMAT+DGDM2(IDG))
             !               print *,'(RAMAT+DGDM3(IDG)):',(RAMAT+DGDM3(IDG))
          END DO
       END DO
       !           stop 393
       RETURN
    endif


    ALLOCATE(RT_BX(TOTELE*NLOC))
    ALLOCATE(RT_BY(TOTELE*NLOC))
    ALLOCATE(RT_BZ(TOTELE*NLOC))
    ALLOCATE(WORKX(TOTELE*NLOC))
    ALLOCATE(WORKY(TOTELE*NLOC))
    ALLOCATE(WORKZ(TOTELE*NLOC))
    WORKX=0.0
    WORKY=0.0
    WORKZ=0.0
    NILOOP=2
    IF(IMM.EQ.-1) NILOOP=1

    do ILOOP=1,NILOOP
       IF(ILOOP.EQ.1) THEN
          ELESTART=1
          ELEFINI=TOTELE
          ELESTEP=1
       ELSE
          ELESTART=TOTELE
          ELEFINI=1
          ELESTEP=-1
       ENDIF
       do ELE=ELESTART,ELEFINI,  ELESTEP
          IF(ILOOP.EQ.1) THEN
             COUNTSTART=FINDRMELE(ELE)
             COUNTFINI=CENTRMELE(ELE)-1
          ELSE
             COUNTSTART=CENTRMELE(ELE)+1  
             COUNTFINI=FINDRMELE(ELE+1)-1
          ENDIF
          do COUNT=COUNTSTART,COUNTFINI
             ELE2=COLMELE(COUNT)
             do ILOC=1,NLOC
                IDG=(ELE-1)*NLOC+ILOC
                do JLOC=1,NLOC
                   JDG=(ELE2-1)*NLOC+JLOC
                   IJDISP=(ILOC-1)*NLOC + JLOC
                   RAMAT=AMATELE((COUNT-1)*NLOC*NLOC+IJDISP)
                   WORKX(IDG)=WORKX(IDG)-RAMAT*RT_BX(JDG)
                   WORKY(IDG)=WORKY(IDG)-RAMAT*RT_BY(JDG)
                   WORKZ(IDG)=WORKZ(IDG)-RAMAT*RT_BZ(JDG)
                END DO
             END DO
          END DO

          ! Calculate matrix block...
          MATBLOK=0.0
          COUNT=CENTRMELE(ELE)
          do ILOC=1,NLOC
             do JLOC=1,NLOC
                ID=1
                JD=ID
                BIG_ILOC=(ID-1)*NLOC+ILOC
                BIG_JLOC=(JD-1)*NLOC+JLOC
                IJDISP=(ILOC-1)*NLOC + JLOC
                RAMAT=AMATELE((COUNT-1)*NLOC*NLOC+IJDISP)
                MATBLOK(BIG_ILOC,BIG_JLOC)=RAMAT
                ID=2
                JD=ID
                BIG_ILOC=(ID-1)*NLOC+ILOC
                BIG_JLOC=(JD-1)*NLOC+JLOC
                IJDISP=(ILOC-1)*NLOC + JLOC
                MATBLOK(BIG_ILOC,BIG_JLOC)=RAMAT
                ID=3
                JD=ID
                BIG_ILOC=(ID-1)*NLOC+ILOC
                BIG_JLOC=(JD-1)*NLOC+JLOC
                IJDISP=(ILOC-1)*NLOC + JLOC
                MATBLOK(BIG_ILOC,BIG_JLOC)=RAMAT
             END DO
             JLOC=ILOC
             ID=1
             JD=1
             BIG_ILOC=(ID-1)*NLOC+ILOC
             BIG_JLOC=(JD-1)*NLOC+JLOC
             MATBLOK(BIG_ILOC,BIG_JLOC)=MATBLOK(BIG_ILOC,BIG_JLOC)+COMAT11(ELE,ILOC)
             ID=1
             JD=2
             BIG_ILOC=(ID-1)*NLOC+ILOC
             BIG_JLOC=(JD-1)*NLOC+JLOC
             MATBLOK(BIG_ILOC,BIG_JLOC)=MATBLOK(BIG_ILOC,BIG_JLOC)+COMAT12(ELE,ILOC)
             ID=1
             JD=3
             BIG_ILOC=(ID-1)*NLOC+ILOC
             BIG_JLOC=(JD-1)*NLOC+JLOC
             MATBLOK(BIG_ILOC,BIG_JLOC)=MATBLOK(BIG_ILOC,BIG_JLOC)+COMAT13(ELE,ILOC)
             ID=2
             JD=1
             BIG_ILOC=(ID-1)*NLOC+ILOC
             BIG_JLOC=(JD-1)*NLOC+JLOC
             MATBLOK(BIG_ILOC,BIG_JLOC)=MATBLOK(BIG_ILOC,BIG_JLOC)+COMAT21(ELE,ILOC)
             ID=2
             JD=2
             BIG_ILOC=(ID-1)*NLOC+ILOC
             BIG_JLOC=(JD-1)*NLOC+JLOC
             MATBLOK(BIG_ILOC,BIG_JLOC)=MATBLOK(BIG_ILOC,BIG_JLOC)+COMAT22(ELE,ILOC)
             ID=2
             JD=3
             BIG_ILOC=(ID-1)*NLOC+ILOC
             BIG_JLOC=(JD-1)*NLOC+JLOC
             MATBLOK(BIG_ILOC,BIG_JLOC)=MATBLOK(BIG_ILOC,BIG_JLOC)+COMAT23(ELE,ILOC)
             ID=3
             JD=1
             BIG_ILOC=(ID-1)*NLOC+ILOC
             BIG_JLOC=(JD-1)*NLOC+JLOC
             MATBLOK(BIG_ILOC,BIG_JLOC)=MATBLOK(BIG_ILOC,BIG_JLOC)+COMAT31(ELE,ILOC)
             ID=3
             JD=2
             BIG_ILOC=(ID-1)*NLOC+ILOC
             BIG_JLOC=(JD-1)*NLOC+JLOC
             MATBLOK(BIG_ILOC,BIG_JLOC)=MATBLOK(BIG_ILOC,BIG_JLOC)+COMAT32(ELE,ILOC)
             ID=3
             JD=3
             BIG_ILOC=(ID-1)*NLOC+ILOC
             BIG_JLOC=(JD-1)*NLOC+JLOC
             MATBLOK(BIG_ILOC,BIG_JLOC)=MATBLOK(BIG_ILOC,BIG_JLOC)+COMAT33(ELE,ILOC)
          END DO
          ! form R*MATBLOK*RT
          
          RMATBLOKRT=MATBLOK
          ! Put in BC's
          do ILOC=1,NLOC
             ID=1
             BIG_ILOC=(ID-1)*NLOC+ILOC
             RMATBLOKRT(BIG_ILOC,BIG_ILOC)=RMATBLOKRT(BIG_ILOC,BIG_ILOC)+DGDM1((ELE-1)*NLOC+ILOC)
             ID=2
             BIG_ILOC=(ID-1)*NLOC+ILOC
             RMATBLOKRT(BIG_ILOC,BIG_ILOC)=RMATBLOKRT(BIG_ILOC,BIG_ILOC)+DGDM2((ELE-1)*NLOC+ILOC)
             ID=3
             BIG_ILOC=(ID-1)*NLOC+ILOC
             RMATBLOKRT(BIG_ILOC,BIG_ILOC)=RMATBLOKRT(BIG_ILOC,BIG_ILOC)+DGDM3((ELE-1)*NLOC+ILOC)
          END DO
          ! Solve MAT X=B (NB MAT is overwritten).  
          do ILOC=1,NLOC
             IDG=(ELE-1)*NLOC+ILOC
             ID=1
             BIG_ILOC=(ID-1)*NLOC+ILOC
             LOCWORK(BIG_ILOC)=WORKX(IDG)
             ID=2
             BIG_ILOC=(ID-1)*NLOC+ILOC
             LOCWORK(BIG_ILOC)=WORKY(IDG)
             ID=3
             BIG_ILOC=(ID-1)*NLOC+ILOC
             LOCWORK(BIG_ILOC)=WORKZ(IDG)
          END DO
          ! RLOCWORK=R*LOCWORK
          RLOCWORK=LOCWORK

          do ILOC=1,NLOC
             IDG=(ELE-1)*NLOC+ILOC
             ID=1
             BIG_ILOC=(ID-1)*NLOC+ILOC
             RLOCWORK(BIG_ILOC)=RLOCWORK(BIG_ILOC)+VECX(IDG)
             ID=2
             BIG_ILOC=(ID-1)*NLOC+ILOC
             RLOCWORK(BIG_ILOC)=RLOCWORK(BIG_ILOC)+VECY(IDG)
             ID=3
             BIG_ILOC=(ID-1)*NLOC+ILOC
             RLOCWORK(BIG_ILOC)=RLOCWORK(BIG_ILOC)+VECZ(IDG)
          END DO
          ! SOLVE BLOCK SYSTEM...
          CALL SMLINN(RMATBLOKRT,LOCB,RLOCWORK,BIG_NLOC,BIG_NLOC)
          !            print *,'RLOCWORK:',RLOCWORK
          !            print *,'locb:',locb
          !            print *,'matrix:'
          !            DO I=1,BIG_NLOC
          !              print *,(RMATBLOKRT(i,j),J=1,BIG_NLOC)
          !            END DO
          !            stop 381
          ! RT_LOCB=RT*LOCB
          RT_LOCB=LOCB

          do ILOC=1,NLOC
             IDG=(ELE-1)*NLOC+ILOC
             ID=1
             BIG_ILOC=(ID-1)*NLOC+ILOC
             BX(IDG)   =LOCB(BIG_ILOC)
             RT_BX(IDG)=RT_LOCB(BIG_ILOC)
             ID=2
             BIG_ILOC=(ID-1)*NLOC+ILOC
             BY(IDG)   =LOCB(BIG_ILOC)
             RT_BY(IDG)=RT_LOCB(BIG_ILOC)
             ID=3
             BIG_ILOC=(ID-1)*NLOC+ILOC
             BZ(IDG)   =LOCB(BIG_ILOC)
             RT_BZ(IDG)=RT_LOCB(BIG_ILOC)
          END DO
          ! Next ELE...
       END DO
       ! Next ILOOP...
    END DO
    RETURN
  END SUBROUTINE DG_SOLVE_PRE_MOM
  !
  !
  !
  !
  SUBROUTINE UVW_ROT_DG_ALWAYS(NLOC,BIG_NLOC,TOTELE,&
       TRANSP,DGU,DGV,DGW,&
       DGT1X,DGT1Y,DGT1Z, &
       DGT2X,DGT2Y,DGT2Z, &
       DGNX,DGNY,DGNZ)
    ! apply U=R U or if TRANSP then U=R^T U
    IMPLICIT NONE
    INTEGER NLOC,TOTELE
    INTEGER BIG_NLOC
    LOGICAL TRANSP
    REAL DGU(TOTELE*NLOC),DGV(TOTELE*NLOC),DGW(TOTELE*NLOC)
    ! place into rotations matrices:
    ! DGT1X,DGT1Y,DGT1Z, 
    ! DGT2X,DGT2Y,DGT2Z, 
    ! DGNX,DGNY,DGNZ,
    ! and b.c's  DGDM1,DGDM2,DGDM3
    REAL DGT1X(TOTELE*NLOC),DGT1Y(TOTELE*NLOC),DGT1Z(TOTELE*NLOC)
    REAL DGT2X(TOTELE*NLOC),DGT2Y(TOTELE*NLOC),DGT2Z(TOTELE*NLOC)
    REAL DGNX(TOTELE*NLOC),DGNY(TOTELE*NLOC),DGNZ(TOTELE*NLOC)
    ! Local variables...
    INTEGER SILOC,ILOC,SGI,DGNODI,ILOC3,JLOC3,NODDG,ELE
    REAL RLOC(BIG_NLOC,BIG_NLOC),RTLOC(BIG_NLOC,BIG_NLOC)
    REAL MATROT(BIG_NLOC,BIG_NLOC)
    REAL SOUSUBU(NLOC),SOUSUBV(NLOC),SOUSUBW(NLOC)


    do ELE=1,TOTELE
       RLOC=0.0
       do ILOC=1,NLOC
          NODDG=(ELE-1)*NLOC+ILOC
          !  The rotation matrix in 3-D is R=  
          !   T1X   T1Y   T1Z
          !   T2X   T2Y   T2Z
          !   NX    NY    NZ
          RLOC(ILOC,ILOC)       =DGT1X(NODDG)
          RLOC(ILOC,ILOC+NLOC)  =DGT1Y(NODDG)
          RLOC(ILOC,ILOC+2*NLOC)=DGT1Z(NODDG)
          RLOC(ILOC+NLOC,ILOC)       =DGT2X(NODDG)
          RLOC(ILOC+NLOC,ILOC+NLOC)  =DGT2Y(NODDG)
          RLOC(ILOC+NLOC,ILOC+2*NLOC)=DGT2Z(NODDG)
          RLOC(ILOC+2*NLOC,ILOC)       =DGNX(NODDG)
          RLOC(ILOC+2*NLOC,ILOC+NLOC)  =DGNY(NODDG)
          RLOC(ILOC+2*NLOC,ILOC+2*NLOC)=DGNZ(NODDG)
       END DO

       IF(.NOT.TRANSP) THEN
          do ILOC3=1,BIG_NLOC
             do JLOC3=1,BIG_NLOC
                RTLOC(ILOC3,JLOC3)=RLOC(JLOC3,ILOC3)
             END DO
          END DO
          MATROT=RTLOC
       ELSE
          MATROT=RLOC
       ENDIF

       do ILOC=1,NLOC
          NODDG=(ELE-1)*NLOC+ILOC
          SOUSUBU(ILOC)=DGU(NODDG)
          SOUSUBV(ILOC)=DGV(NODDG)
          SOUSUBW(ILOC)=DGW(NODDG)

          ! Rotate rhs vecs and set to zero (for zero bc)...
          DGU(NODDG)=MATROT(ILOC,ILOC)*SOUSUBU(ILOC)&
               +MATROT(ILOC,ILOC+NLOC)*SOUSUBV(ILOC)&
               +MATROT(ILOC,ILOC+2*NLOC)*SOUSUBW(ILOC)   
          DGV(NODDG)=MATROT(ILOC+NLOC,ILOC)*SOUSUBU(ILOC)&
               +MATROT(ILOC+NLOC,ILOC+NLOC)*SOUSUBV(ILOC)&
               +MATROT(ILOC+NLOC,ILOC+2*NLOC)*SOUSUBW(ILOC)
          DGW(NODDG)=MATROT(ILOC+2*NLOC,ILOC)*SOUSUBU(ILOC)&
               +MATROT(ILOC+2*NLOC,ILOC+NLOC)*SOUSUBV(ILOC)&
               +MATROT(ILOC+2*NLOC,ILOC+2*NLOC)*SOUSUBW(ILOC)
       END DO
    END DO
    RETURN
  END SUBROUTINE UVW_ROT_DG_ALWAYS
end module mulmat_module
