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
!
!
!
module rotated_boundary_conditions_legacy
use fldebug
use flcomms_module
implicit none

contains
  
         SUBROUTINE GETROT( GETTAN, &
     &         NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
     &         NNODRO,D3)
! Find a default set of tangent vectors. 
! This subroutine works out the TANGENT vectors for rotation from 
! the commponents of the normal direction to the boundary(NX,NY,NZ)
! This sub is only called if GETTAN is true. 
! NODROT contains the nodes to rotate. 
! NB The normal is pointing out of the domain. 
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
!  The rotation matrix in 2-D is R=  
!   T1X   T1Y   
!   NX    NY    
! 
         REAL TOLER
         PARAMETER(TOLER=1.E-4)
       INTEGER NNODRO
       REAL NX(NNODRO),NY(NNODRO),NZ(NNODRO)
       REAL T1X(NNODRO),T1Y(NNODRO),T1Z(NNODRO)
       REAL T2X(NNODRO),T2Y(NNODRO),T2Z(NNODRO)

       INTEGER II, NOD
       REAL RN
       LOGICAL D3,GETTAN
         
          ewrite(1,*)'just inside  getrot'

       IF(GETTAN) THEN
! &&&&&&&&&&&&&&
       IF(.NOT.D3) THEN
      do  II=1,NNODRO! Was loop 10
            T1X(II)=-NY(II)
            T1Y(II)=NX(II)
!             ewrite(3,*) 'NX(II),NY(II):',NX(II),NY(II)
!             ewrite(3,*) 'T1X(II),T1Y(II):',T1X(II),T1Y(II)
      end do ! Was loop 10
!           stop 3831
       ENDIF
       IF(D3) THEN
      do  II=1,NNODRO! Was loop 20
!          ewrite(3,*)'ii,NX(II),NY(II),NZ(II):',
!     &                  ii,NX(II),NY(II),NZ(II)
            RN=SQRT(NX(II)**2 + NY(II)**2)
!              IF(RN.GT.TOLER) THEN
              IF(.false.) THEN
            T1X(II)=NY(II)/RN
            T1Y(II)=-NX(II)/RN
            T1Z(II)=0.
!
            T2X(II)=NX(II)*NZ(II)/RN
            T2Y(II)=NY(II)*NZ(II)/RN
            T2Z(II)=-RN
              ELSE
! Normal approx alighned with z-axis
              RN=SQRT(NY(II)**2 + NZ(II)**2)
!
!                IF(RN.GT.TOLER) THEN
                IF(.false.) THEN
                   T1X(II)=0.
                   T1Y(II)= NZ(II)/RN
                   T1Z(II)=-NY(II)/RN
!
                   T2X(II)=-RN
                   T2Y(II)=NX(II)*NY(II)/RN
                   T2Z(II)=NX(II)*NZ(II)/RN
                ELSE
! Normal approx alighned with x-axis as well
! the thing to use here...
                RN=SQRT(NX(II)**2 + NZ(II)**2)
                T1X(II)=-NZ(II)/RN
                T1Y(II)=0.
                T1Z(II)= NX(II)/RN

              T2X(II)= NX(II)*NY(II)/RN
              T2Y(II)=-RN
              T2Z(II)= NY(II)*NZ(II)/RN
                ENDIF
 
              ENDIF
      end do ! Was loop 20
       ENDIF
! &&&&&&&&&&&&&& -ENDOF IF(GETTAN) 
       ENDIF
          ewrite(3,*)'finished getrot'
       END subroutine

         SUBROUTINE ROTOCE( GETTAN, NODROT,&
     &         NOBCU,BCU2,NOBCV,BCV2,NOBCW,BCW2,&
     &         NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
     &         NNODRO,D3)
! Find a default set of tangent vectors. 
! This subroutine works out the TANGENT vectors for rotation from 
! the commponents of the normal direction to the boundary(NX,NY,NZ)
! This sub is only called if GETTAN is true. 
! NODROT contains the nodes to rotate. 
! NB The normal is pointing out of the domain. 
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
!  The rotation matrix in 2-D is R=  
!   T1X   T1Y   
!   NX    NY    
! This sub rotates for ocean modelling. 
       INTEGER II,NOD2,JJ,NOD
       REAL RN,TOLER
       PARAMETER(TOLER=1.E-4)
       INTEGER NNODRO
       INTEGER NODROT(NNODRO)
       INTEGER NOBCU,NOBCV,NOBCW
       INTEGER BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
       REAL NX(NNODRO),NY(NNODRO),NZ(NNODRO)
       REAL T1X(NNODRO),T1Y(NNODRO),T1Z(NNODRO)
       REAL T2X(NNODRO),T2Y(NNODRO),T2Z(NNODRO)

       LOGICAL D3,GETTAN
         
       ewrite(2,*)'just inside  rotoce NNODRO:',NNODRO
       ewrite(2,*)'NOBCU,NOBCV,NOBCW:',NOBCU,NOBCV,NOBCW

       if(GETTAN) THEN
          do II=1,NNODRO
                CALL GETROTSIM(D3, &
     &               T1X(II),T1Y(II),T1Z(II), &
     &               T2X(II),T2Y(II),T2Z(II), &
     &               NX(II),NY(II),NZ(II))
          END DO
       ENDIF
       
       END subroutine rotoce

       SUBROUTINE GETROTSIM(D3, &
     &         T1X,T1Y,T1Z, &
     &         T2X,T2Y,T2Z, &
     &         NX,NY,NZ)
! Find a default set of tangent vectors. 
! This subroutine works out the TANGENT vectors for rotation from 
! the commponents of the normal direction to the boundary(NX,NY,NZ)
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
!  The rotation matrix in 2-D is R=  
!   T1X   T1Y   
!   NX    NY    
! 
         REAL TOLER
         PARAMETER(TOLER=1.E-4)
         LOGICAL D3
         REAL T1X,T1Y,T1Z
         REAL T2X,T2Y,T2Z
         REAL NX,NY,NZ
! Local variables...
         REAL RNX,RNY,RNZ
         REAL RN
    
! &&&&&&&&&&&&&&
          IF(.NOT.D3) THEN
               T1X=-NY
               T1Y=NX
          ENDIF
          IF(D3) THEN
    
         RNZ=NX**2 + NY**2
         RNX=NY**2 + NZ**2
         RNY=NX**2 + NZ**2
! Choose the values of RN? that is largest...

         IF((RNZ.GT.RNX).AND.(RNZ.GT.RNY)) THEN
    
                 RN=SQRT(NX**2 + NY**2)
                 T1X=NY/RN
                 T1Y=-NX/RN
                 T1Z=0.

                 T2X=NX*NZ/RN
                 T2Y=NY*NZ/RN
                 T2Z=-RN

         ELSE IF(RNX.GT.RNY) THEN
              
! Normal approx alighned with z-axis
                 RN=SQRT(NY**2 + NZ**2)
                 T1X=0.
                 T1Y= NZ/RN
                 T1Z=-NY/RN

                 T2X=-RN
                 T2Y=NX*NY/RN
                 T2Z=NX*NZ/RN
         ELSE
                
! Normal approx alighned with x-axis as well
! the thing to use here...
                   RN=SQRT(NX**2 + NZ**2)
                   T1X=-NZ/RN
                   T1Y=0.
                   T1Z= NX/RN
                   T2X= NX*NY/RN
                   T2Y=-RN
                   T2Z= NY*NZ/RN
         ENDIF
          
          ENDIF
          
          END subroutine
            
SUBROUTINE GETROTGRAV(&
               T1X,T1Y,T1Z, &
               T2X,T2Y,T2Z, &
               NORMX,NORMY,NORMZ, &
               COGRAX,BSOUX,BSOUY,BSOUZ, &
               VGRAVX,VGRAVY,VGRAVZ, xnonod)
! calculate Rotation matrix R from taking 
! the normal direction to be �ve the direction of gravity. 
       INTEGER xnonod
! R=
       REAL, intent(out):: T1X(xnonod),T1Y(xnonod),T1Z(xnonod)
       REAL, intent(out):: T2X(xnonod),T2Y(xnonod),T2Z(xnonod) 
       REAL, intent(out):: NORMX(xnonod),NORMY(xnonod),NORMZ(xnonod) 
       LOGICAL, intent(in):: COGRAX
       REAL, intent(in):: BSOUX,BSOUY,BSOUZ
       REAL, intent(in):: VGRAVX(xnonod),VGRAVY(xnonod),VGRAVZ(xnonod)
! Local variables�
       REAL RN
       INTEGER IGL
       DO IGL=1,xnonod
! Calculate rotations...
         IF(COGRAX) THEN
           NORMX(IGL)=-BSOUX
           NORMY(IGL)=-BSOUY
           NORMZ(IGL)=-BSOUZ
         ELSE
           NORMX(IGL)=-VGRAVX(IGL)
           NORMY(IGL)=-VGRAVY(IGL)
           NORMZ(IGL)=-VGRAVZ(IGL)
!           ewrite(3,*) 'igl,3s:',VGRAVX(IGL),VGRAVY(IGL),VGRAVZ(IGL)
         ENDIF
         RN=SQRT(NORMX(IGL)**2+NORMY(IGL)**2+NORMZ(IGL)**2)
         NORMX(IGL)=NORMX(IGL)/RN
         NORMY(IGL)=NORMY(IGL)/RN
         NORMZ(IGL)=NORMZ(IGL)/RN
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NORMX NORMY NORMZ
! Calculate T1X,T1Y,T1Z, T2X,T2Y,T2Z from NORMX,NORMY,NORMZ

         CALL GETROTSIM(.TRUE., &
               T1X(IGL),T1Y(IGL),T1Z(IGL), &
               T2X(IGL),T2Y(IGL),T2Z(IGL), &
               NORMX(IGL),NORMY(IGL),NORMZ(IGL))
       END DO
       RETURN
       END subroutine
         
      SUBROUTINE GETNOR( &
     &     NX,NY,NZ, &
     &     NODROT,NNODRO,D3,&
     &     NONODS,FREDOP,&
     &     C1T,C2T,C3T,FINDCT,COLCT,NCT,&
     &     PARA,halo_tag,&
     &     NNODP)
!     because of boundary conditions we make normal point into domain. 
!     Find a default set of normal vectors. 
!     This sub is only called if GETTAN is true. 
!     NODROT contains the nodes to rotate. 
!     NB The normal is pointing out of the domain. 
!     The rotation matrix in 3-D is R=  
!     T1X   T1Y   T1Z
!     T2X   T2Y   T2Z
!     NX    NY    NZ
!     The rotation matrix in 2-D is R=  
!     T1X   T1Y   
!     NX    NY    
!     SHAPNX,SHAPNY,SHAPNZ are the working arrays.
      INTEGER NCT
      REAL TOLER
      PARAMETER(TOLER=1.E-4)
      INTEGER NNODRO,NONODS,FREDOP
      REAL NX(NNODRO),NY(NNODRO),NZ(NNODRO)
      REAL, dimension(:), allocatable:: SHAPNX,SHAPNY,SHAPNZ, spair, spair2
      INTEGER NODROT(NNODRO)
      REAL C1T(NCT),C2T(NCT),C3T(NCT)
      INTEGER COL
      INTEGER I,J,FINDCT(FREDOP+1),COLCT(NCT)
!     parallel stuff...
      INTEGER PARA,halo_tag
      INTEGER NNODP
      LOGICAL D3
      INTEGER NDIM
      REAL RR,RN
      INTEGER NOD,II
      
      allocate(shapnx(1:nonods), shapny(1:nonods), shapnz(1:nonods), &
        spair(3*nonods), spair2(3*nonods))
!     
!     Work out the shape functions at all the nodes.
!     SHAPNX(i) is the derivative of shape function i in X-direction.
      IF(.NOT.D3) THEN
      do NOD=1,NONODS
            SHAPNX(NOD)=0.
            SHAPNY(NOD)=0.
         END DO
      do J=1,FREDOP
      do I=FINDCT(J),FINDCT(J+1)-1
               COL=COLCT(I)
               SHAPNX(COL)=SHAPNX(COL)+C1T(I)
               SHAPNY(COL)=SHAPNY(COL)+C2T(I)
            END DO
         END DO
      ELSE
      do NOD=1,NONODS
            SHAPNX(NOD)=0.
            SHAPNY(NOD)=0.
            SHAPNZ(NOD)=0.
         END DO
      do J=1,FREDOP
      do I=FINDCT(J),FINDCT(J+1)-1
               COL=COLCT(I)
               SHAPNX(COL)=SHAPNX(COL)+C1T(I)
               SHAPNY(COL)=SHAPNY(COL)+C2T(I)
               SHAPNZ(COL)=SHAPNZ(COL)+C3T(I)
            END DO
         END DO
      ENDIF
!     
      IF(D3) THEN
      do II=1,NNODRO
            NOD=NODROT(II)
            NX(II)=SHAPNX(NOD)
            NY(II)=SHAPNY(NOD)
            NZ(II)=SHAPNZ(NOD)
            RN=SQRT(NX(II)**2+NY(II)**2+NZ(II)**2)
!     -values so normal pointing into domain
            NX(II)=-NX(II)/RN
            NY(II)=-NY(II)/RN
            NZ(II)=-NZ(II)/RN
         END DO
      ELSE
      do II=1,NNODRO
            NOD=NODROT(II)
            NX(II)=SHAPNX(NOD)
            NY(II)=SHAPNY(NOD)
            RN=SQRT(NX(II)**2+NY(II)**2)
!     -values so normal pointing into domain
            NX(II)=-NX(II)/RN
            NY(II)=-NY(II)/RN
         END DO
      ENDIF
      IF(PARA.EQ.1) THEN 
         NDIM=2
         IF(D3) NDIM=3
         SPAIR(1:NDIM*NONODS) = 0.0
      do II=1,NNODRO
            NOD=NODROT(II)
            SPAIR(NOD)         =NX(II)
            SPAIR(NOD+NONODS)  =NY(II)
            IF(D3) SPAIR(NOD+2*NONODS)=NZ(II)
         END DO
         CALL HALGET(SPAIR,NONODS,NDIM*NONODS,NNODP,&
     &        halo_tag)
      do II=1,NNODRO
            NOD=NODROT(II)
            RR=SPAIR(NOD)**2+SPAIR(NOD+NONODS)**2
            IF(D3) RR=RR+SPAIR(NOD+2*NONODS)**2
            IF(RR.GT.0.3) THEN
               NX(II)       =SPAIR(NOD)
               NY(II)       =SPAIR(NOD+NONODS)
               IF(D3) NZ(II)=SPAIR(NOD+2*NONODS)
            ENDIF
         END DO
      END IF
      
      deallocate(shapnx, shapny, shapnz, spair, spair2)
      
      END subroutine getnor

  SUBROUTINE ROTUVW(U,V,W,NONODS,&
     &         NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
     &         NODROT,NNODRO,D3,BACK)
! IF .NOT.BACK This sub rotates U,V,W alighned with the X,Y,Z coordinates
! to U,V,W ALIGHNED WITH THE BOUNDARY. 
!
! IF BACK then rotate u,v,w ALIGHNED WITH THE BOUNDARY 
! to U,V,W alighned with the X,Y,Z coordinates.
! NODROT contains the nodes to rotate. 
! NB The normal is pointing out of the domain. 
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
!  The rotation matrix in 2-D is R=  
!   T1X   T1Y   
!   NX    NY    
! 
       INTEGER NNODRO,NONODS
       REAL NX(NNODRO),NY(NNODRO),NZ(NNODRO)
       REAL T1X(NNODRO),T1Y(NNODRO),T1Z(NNODRO)
       REAL T2X(NNODRO),T2Y(NNODRO),T2Z(NNODRO)
       INTEGER NODROT(NNODRO)
!
       REAL U(NONODS),V(NONODS),W(NONODS)
       LOGICAL D3,BACK
       INTEGER II,NOD
       REAL RU,RV,RW
       
  IF(BACK) THEN
! perform R q. 
    IF(.NOT.D3) THEN
      do  II=1,NNODRO! Was loop 107
           NOD=NODROT(II)
           RU=U(NOD)
           RV=V(NOD)
           U(NOD)=RU*T1X(II) + RV*T1Y(II)
           V(NOD)=RU*NX(II) + RV*NY(II)
      end do ! Was loop 107
    ENDIF
    IF(D3) THEN
      do  II=1,NNODRO! Was loop 207
           NOD=NODROT(II)
           RU=U(NOD)
           RV=V(NOD)
           RW=W(NOD)
           U(NOD)=RU*T1X(II) + RV*T1Y(II) + RW*T1Z(II)
           V(NOD)=RU*T2X(II) + RV*T2Y(II) + RW*T2Z(II)
           W(NOD)=RU*NX(II)  + RV*NY(II)  + RW*NZ(II)
      end do ! Was loop 207
    ENDIF
  ENDIF

  IF(.NOT.BACK) THEN
! perform R^T q. 
    IF(.NOT.D3) THEN
      do  II=1,NNODRO! Was loop 10
           NOD=NODROT(II)
           RU=U(NOD)
           RV=V(NOD)
           U(NOD)=RU*T1X(II) + RV*NX(II)
           V(NOD)=RU*T1Y(II) + RV*NY(II)
      end do ! Was loop 10
    ENDIF
    IF(D3) THEN
      do  II=1,NNODRO! Was loop 20
           NOD=NODROT(II)
           RU=U(NOD)
           RV=V(NOD)
           RW=W(NOD)
           U(NOD)=RU*T1X(II) + RV*T2X(II) + RW*NX(II)
           V(NOD)=RU*T1Y(II) + RV*T2Y(II) + RW*NY(II)
           W(NOD)=RU*T1Z(II) + RV*T2Z(II) + RW*NZ(II)
      end do ! Was loop 20
    ENDIF
  ENDIF
  
  END subroutine rotuvw
       
  SUBROUTINE RBIGRT(BIGM1,CENTRM,FINDRM,COLM,&
     &         NCOLM,NBIGM,NONODS,&
     &         NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
     &         NODROT,NNODRO,D3)
         IMPLICIT NONE
! This subroutine works out R BIGM1 R^T where U,V and possibly W 
! are coupled through cross derivatives. 
! NODROT contains the nodes to rotate. 
! NB The normal is pointing out of the domain. 
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
!  The rotation matrix in 2-D is R=  
!   T1X   T1Y   
!   NX    NY    
! 
         INTEGER NONODS,NCOLM,NBIGM,NNODRO
         REAL NX(NNODRO),NY(NNODRO),NZ(NNODRO)
         REAL T1X(NNODRO),T1Y(NNODRO),T1Z(NNODRO)
         REAL T2X(NNODRO),T2Y(NNODRO),T2Z(NNODRO)
         INTEGER NODROT(NNODRO)

         REAL BIGM1(NBIGM)
         INTEGER CENTRM(NONODS),FINDRM(NONODS+1),COLM(NCOLM)
         INTEGER ROW,COUNT,COL,COUNT2
         LOGICAL D3

        REAL NXII,NYII,NZII, T1XII,T1YII,T1ZII, T2XII,T2YII,T2ZII
        REAL A11,A12,A13, A21,A22,A23, A31,A32,A33
        REAL SA11,SA12,SA13, SA21,SA22,SA23, SA31,SA32,SA33
        REAL A,B,C, D,E,F, G,H,K
        LOGICAL BLKSYM
        INTEGER II,NOD,ICENT,ICOUNT
       
 
      BLKSYM=.FALSE.
      IF((.NOT.D3).AND.(3*NCOLM.EQ.NBIGM))BLKSYM=.TRUE.
      IF(       D3.AND.(6*NCOLM.EQ.NBIGM))BLKSYM=.TRUE.

      IF(BLKSYM) THEN

        IF(.NOT.D3) THEN
          do  II=1,NNODRO! Was loop 60
             NOD=NODROT(II)
             ICENT=CENTRM(NOD)
             SA11=BIGM1(ICENT)
             SA12=BIGM1(ICENT+NCOLM)
             SA21=BIGM1(ICENT+NCOLM)
             SA22=BIGM1(ICENT+2*NCOLM)
             NXII=NX(II)
             NYII=NY(II)
             T1XII=T1X(II)
             T1YII=T1Y(II)
! In the matrix find the ROWS for which COLUMN NOD is non-zero.
             do  COUNT=FINDRM(NOD),FINDRM(NOD+1)-1! Was loop 70
                  ROW=COLM(COUNT)
! Find COUNT2 which is the transpose position of COUNT. 
! Can replace next 3 lines with the efficient serching algorithm in HART. 
                  do  ICOUNT=FINDRM(ROW),FINDRM(ROW+1)-1
                      IF(COLM(ICOUNT).EQ.NOD) COUNT2=ICOUNT
                  end do
! Now form BIGM1*R^T.
                  A11=BIGM1(COUNT2)
                  A12=BIGM1(COUNT+NCOLM)
                  A21=BIGM1(COUNT2+NCOLM)
                  A22=BIGM1(COUNT2+2*NCOLM)
! Form A*RT. 
                  BIGM1(COUNT2)=A11*T1XII&
     &                         +A12*T1YII
! Miss out the transpose 
                  BIGM1(COUNT2+NCOLM)=A11*NXII&
     &                         +A12*NYII
                  BIGM1(COUNT2+2*NCOLM)=A21*NXII&
     &                                         +A22*NYII
!
! Now form R * BIGM1.
                  A11=BIGM1(COUNT)
                  A12=BIGM1(COUNT2+NCOLM)
                  A21=BIGM1(COUNT+NCOLM)
                  A22=BIGM1(COUNT+2*NCOLM)
! Form R * (A*RT) where A*RT is contained in Aij.
                  BIGM1(COUNT)=A11*T1XII+A21*T1YII
! Miss out the transpose 
!                 BIGM1(COUNT2+NCOLM)=A12*T1XII+A22*T1YII
                  BIGM1(COUNT+NCOLM)=A11*NXII+A21*NYII
                  BIGM1(COUNT+2*NCOLM)=A12*NXII+A22*NYII
             end do
! Correct some values.(At the cross points).  R*(SA*RT)
             A=SA11*T1XII+SA12*T1YII
             B=SA11*NXII+SA12*NYII
             C=SA21*T1XII+SA22*T1YII
             D=SA21*NXII+ SA22*NYII
! 
             BIGM1(ICENT)=T1XII*A+T1YII*C
             BIGM1(ICENT+NCOLM)=NXII*A+NYII*B
             BIGM1(ICENT+2*NCOLM)=NXII*B+NYII*D
          end do ! Was loop 60
        ENDIF
        
        IF(D3) THEN
          do  II=1,NNODRO! Was loop 160
             NOD=NODROT(II)
             ICENT=CENTRM(NOD)
             SA11=BIGM1(ICENT)
             SA21=BIGM1(ICENT+NCOLM)
             SA22=BIGM1(ICENT+2*NCOLM)
!
             SA31=BIGM1(ICENT+3*NCOLM)
             SA32=BIGM1(ICENT+4*NCOLM)
             SA33=BIGM1(ICENT+5*NCOLM)
!
             SA12=SA21
             SA13=SA31
             SA23=SA32
             NXII=NX(II)
             NYII=NY(II)
             NZII=NZ(II)
!
             T1XII=T1X(II)
             T1YII=T1Y(II)
             T1ZII=T1Z(II)
!
             T2XII=T2X(II)
             T2YII=T2Y(II)
             T2ZII=T2Z(II)
! In the matrix find the ROWS for which COLUMN NOD is non-zero.
             do  COUNT=FINDRM(NOD),FINDRM(NOD+1)-1! Was loop 170
                  ROW=COLM(COUNT)
! Find COUNT2 which is the transpose position of COUNT. 
! Can replace next 3 lines with the efficient serching algorithm in HART. 
                  do  ICOUNT=FINDRM(ROW),FINDRM(ROW+1)-1! Was loop 180
                      IF(COLM(ICOUNT).EQ.NOD) COUNT2=ICOUNT
                  end do ! Was loop 180
! Now form BIGM1*R^T.
                  A11=BIGM1(COUNT2)
                  A12=BIGM1(COUNT+NCOLM)
                  A13=BIGM1(COUNT+3*NCOLM)
!
                  A21=BIGM1(COUNT2+NCOLM)
                  A22=BIGM1(COUNT2+2*NCOLM)
                  A23=BIGM1(COUNT+4*NCOLM)
!
                  A31=BIGM1(COUNT2+3*NCOLM)
                  A32=BIGM1(COUNT2+4*NCOLM)
                  A33=BIGM1(COUNT2+5*NCOLM)
! Form A*RT. 
                  BIGM1(COUNT2)        =A11*T1XII+A12*T1YII+A13*T1ZII
                  BIGM1(COUNT2+NCOLM)  =A21*T1XII+A22*T1YII+A23*T1ZII
                  BIGM1(COUNT2+3*NCOLM)=A31*T1XII+A32*T1YII+A33*T1ZII
! Miss out the transposeS
                  BIGM1(COUNT2+2*NCOLM)=A21*T2XII+A22*T2YII+A23*T2ZII
                  BIGM1(COUNT2+4*NCOLM)=A31*T2XII+A32*T2YII+A33*T2ZII
!
                  BIGM1(COUNT2+5*NCOLM)=A31*NXII+A32*NYII+A33*NZII
!
! Now form R * BIGM1.
                  A11=BIGM1(COUNT)
                  A12=BIGM1(COUNT2+NCOLM)
                  A13=BIGM1(COUNT2+3*NCOLM)
!
                  A21=BIGM1(COUNT+NCOLM)
                  A22=BIGM1(COUNT+2*NCOLM)
                  A23=BIGM1(COUNT2+4*NCOLM)
!
                  A31=BIGM1(COUNT+3*NCOLM)
                  A32=BIGM1(COUNT+4*NCOLM)
                  A33=BIGM1(COUNT+5*NCOLM)
! Form R * (A*RT) where A*RT is contained in Aij.
                  BIGM1(COUNT)        =A11*T1XII+A21*T1YII+A31*T1ZII
                  BIGM1(COUNT+NCOLM)  =A11*T2XII+A21*T2YII+A31*T2ZII
                  BIGM1(COUNT+3*NCOLM)=A11*NXII+A21*NYII+A31*NZII
! Miss out the transposeS
                  BIGM1(COUNT+2*NCOLM)=A12*T2XII+A22*T2YII+A32*T2ZII
                  BIGM1(COUNT+4*NCOLM)=A12*NXII+A22*NYII+A32*NZII
!
                  BIGM1(COUNT2+5*NCOLM)=A13*NXII+A23*NYII+A33*NZII
             end do ! Was loop 170
! Correct some values.(At the cross points).  R*(SA*RT)
             A=SA11*T1XII+SA12*T1YII+SA13*T1ZII
             B=SA11*T2XII+SA12*T2YII+SA13*T2ZII
             C=SA11*NXII+SA12*NYII+SA13*NZII
!
             D=SA21*T1XII+SA22*T1YII+SA23*T1ZII
             E=SA21*T2XII+SA22*T2YII+SA23*T2ZII
             F=SA21*NXII+SA22*NYII+SA23*NZII
!
             G=SA31*T1XII+SA32*T1YII+SA33*T1ZII
             H=SA31*T2XII+SA32*T2YII+SA33*T2ZII
             K=SA31*NXII+SA32*NYII+SA33*NZII
!
             BIGM1(ICENT)        =T1XII*A+T1YII*D+T1ZII*G
             BIGM1(ICENT+NCOLM)  =T2XII*A+T2YII*D+T2ZII*G
             BIGM1(ICENT+3*NCOLM)=NXII*A+NYII*D+NZII*G
!
             BIGM1(ICENT+2*NCOLM)=T2XII*B+T2YII*E+T2ZII*H
             BIGM1(ICENT+4*NCOLM)=NXII*B+NYII*E+NZII*H
!
             BIGM1(ICENT+5*NCOLM)=NXII*C+NYII*F+NZII*K
!
        end do ! Was loop 160
      ENDIF
! IF(BLKSYM)ELSE...
      ELSE
        IF(.NOT.D3) THEN
          do  II=1,NNODRO! Was loop 560
             NOD=NODROT(II)
             ICENT=CENTRM(NOD)
             SA11=BIGM1(ICENT)
             SA12=BIGM1(ICENT+NCOLM)
             SA21=BIGM1(ICENT+2*NCOLM)
             SA22=BIGM1(ICENT+3*NCOLM)
             NXII=NX(II)
             NYII=NY(II)
             T1XII=T1X(II)
             T1YII=T1Y(II)
! In the matrix find the ROWS for which COLUMN NOD is non-zero.
             do  COUNT=FINDRM(NOD),FINDRM(NOD+1)-1! Was loop 570
                  ROW=COLM(COUNT)
! Find COUNT2 which is the transpose position of COUNT. 
! Can replace next 3 lines with the efficient serching algorithm in HART. 
                  do  ICOUNT=FINDRM(ROW),FINDRM(ROW+1)-1! Was loop 580
                      IF(COLM(ICOUNT).EQ.NOD) COUNT2=ICOUNT
      end do ! Was loop 580
! Now form BIGM1*R^T.
                  A11=BIGM1(COUNT2)
                  A12=BIGM1(COUNT2+NCOLM)
                  A21=BIGM1(COUNT2+2*NCOLM)
                  A22=BIGM1(COUNT2+3*NCOLM)
! Form A*RT. 
                  BIGM1(COUNT2)        =A11*T1XII+A12*T1YII
                  BIGM1(COUNT2+NCOLM)  =A11*NXII +A12*NYII
                  BIGM1(COUNT2+2*NCOLM)=A21*T1XII+A22*T1YII
                  BIGM1(COUNT2+3*NCOLM)=A21*NXII +A22*NYII
!
! Now form R * BIGM1.
                  A11=BIGM1(COUNT)
                  A12=BIGM1(COUNT+NCOLM)
                  A21=BIGM1(COUNT+2*NCOLM)
                  A22=BIGM1(COUNT+3*NCOLM)
! Form R * (A*RT) where A*RT is contained in Aij.
                  BIGM1(COUNT)        =A11*T1XII+A21*T1YII
                  BIGM1(COUNT+NCOLM)  =A12*T1XII+A22*T1YII
                  BIGM1(COUNT+2*NCOLM)=A11*NXII +A21*NYII
                  BIGM1(COUNT+3*NCOLM)=A12*NXII +A22*NYII
      end do ! Was loop 570
! Correct some values.(At the cross points).  R*(SA*RT)
              A=SA11*T1XII+SA12*T1YII
              B=SA11*NXII +SA12*NYII
              C=SA21*T1XII+SA22*T1YII
              D=SA21*NXII +SA22*NYII
! 
              BIGM1(ICENT)        =T1XII*A+T1YII*C
              BIGM1(ICENT+NCOLM)  =T1XII*B+T1YII*D
              BIGM1(ICENT+2*NCOLM)=NXII *A+NYII*C
              BIGM1(ICENT+3*NCOLM)=NXII *B+NYII*D
      end do ! Was loop 560
         ENDIF
        
         IF(D3) THEN
      do  II=1,NNODRO! Was loop 5160
             NOD=NODROT(II)
             ICENT=CENTRM(NOD)
             SA11=BIGM1(ICENT)
             SA12=BIGM1(ICENT+NCOLM)
             SA13=BIGM1(ICENT+2*NCOLM)
!
             SA21=BIGM1(ICENT+3*NCOLM)
             SA22=BIGM1(ICENT+4*NCOLM)
             SA23=BIGM1(ICENT+5*NCOLM)
!
             SA31=BIGM1(ICENT+6*NCOLM)
             SA32=BIGM1(ICENT+7*NCOLM)
             SA33=BIGM1(ICENT+8*NCOLM)
!
             NXII=NX(II)
             NYII=NY(II)
             NZII=NZ(II)
!
             T1XII=T1X(II)
             T1YII=T1Y(II)
             T1ZII=T1Z(II)
!
             T2XII=T2X(II)
             T2YII=T2Y(II)
             T2ZII=T2Z(II)
! In the matrix find the ROWS for which COLUMN NOD is non-zero.
      do  COUNT=FINDRM(NOD),FINDRM(NOD+1)-1! Was loop 5170
                  ROW=COLM(COUNT)
! Find COUNT2 which is the transpose position of COUNT. 
! Can replace next 3 lines with the efficient serching algorithm in HART. 
      do  ICOUNT=FINDRM(ROW),FINDRM(ROW+1)-1! Was loop 5180
                      IF(COLM(ICOUNT).EQ.NOD) COUNT2=ICOUNT
      end do ! Was loop 5180
! Now form BIGM1*R^T.
                  A11=BIGM1(COUNT2        )
                  A12=BIGM1(COUNT2+  NCOLM)
                  A13=BIGM1(COUNT2+2*NCOLM)
!
                  A21=BIGM1(COUNT2+3*NCOLM)
                  A22=BIGM1(COUNT2+4*NCOLM)
                  A23=BIGM1(COUNT2+5*NCOLM)
!
                  A31=BIGM1(COUNT2+6*NCOLM)
                  A32=BIGM1(COUNT2+7*NCOLM)
                  A33=BIGM1(COUNT2+8*NCOLM)
! Form A*RT. 
                  BIGM1(COUNT2)        =A11*T1XII+A12*T1YII+A13*T1ZII
                  BIGM1(COUNT2+  NCOLM)=A11*T2XII+A12*T2YII+A13*T2ZII
                  BIGM1(COUNT2+2*NCOLM)=A11*NXII +A12*NYII +A13*NZII
!
                  BIGM1(COUNT2+3*NCOLM)=A21*T1XII+A22*T1YII+A23*T1ZII
                  BIGM1(COUNT2+4*NCOLM)=A21*T2XII+A22*T2YII+A23*T2ZII
                  BIGM1(COUNT2+5*NCOLM)=A21*NXII +A22*NYII +A23*NZII
!
                  BIGM1(COUNT2+6*NCOLM)=A31*T1XII+A32*T1YII+A33*T1ZII
                  BIGM1(COUNT2+7*NCOLM)=A31*T2XII+A32*T2YII+A33*T2ZII
                  BIGM1(COUNT2+8*NCOLM)=A31*NXII +A32*NYII +A33*NZII
!
! Now form R * BIGM1.
                  A11=BIGM1(COUNT        )
                  A12=BIGM1(COUNT+  NCOLM)
                  A13=BIGM1(COUNT+2*NCOLM)
!
                  A21=BIGM1(COUNT+3*NCOLM)
                  A22=BIGM1(COUNT+4*NCOLM)
                  A23=BIGM1(COUNT+5*NCOLM)
!
                  A31=BIGM1(COUNT+6*NCOLM)
                  A32=BIGM1(COUNT+7*NCOLM)
                  A33=BIGM1(COUNT+8*NCOLM)
! Form R * (A*RT) where A*RT is contained in Aij.
                  BIGM1(COUNT        )=A11*T1XII+A21*T1YII+A31*T1ZII
                  BIGM1(COUNT+NCOLM  )=A12*T1XII+A22*T1YII+A32*T1ZII
                  BIGM1(COUNT+2*NCOLM)=A13*T1XII+A23*T1YII+A33*T1ZII
!
                  BIGM1(COUNT+3*NCOLM)=A11*T2XII+A21*T2YII+A31*T2ZII
                  BIGM1(COUNT+4*NCOLM)=A12*T2XII+A22*T2YII+A32*T2ZII
                  BIGM1(COUNT+5*NCOLM)=A13*T2XII+A23*T2YII+A33*T2ZII
!
                  BIGM1(COUNT+6*NCOLM)=A11*NXII+A21*NYII+A31*NZII
                  BIGM1(COUNT+7*NCOLM)=A12*NXII+A22*NYII+A32*NZII
                  BIGM1(COUNT+8*NCOLM)=A13*NXII+A23*NYII+A33*NZII
      end do ! Was loop 5170
! Correct some values.(At the cross points).  R*(SA*RT)
              A=SA11*T1XII+SA12*T1YII+SA13*T1ZII
              B=SA11*T2XII+SA12*T2YII+SA13*T2ZII
              C=SA11*NXII +SA12*NYII +SA13*NZII
!
              D=SA21*T1XII+SA22*T1YII+SA23*T1ZII
              E=SA21*T2XII+SA22*T2YII+SA23*T2ZII
              F=SA21*NXII +SA22*NYII +SA23*NZII
!
              G=SA31*T1XII+SA32*T1YII+SA33*T1ZII
              H=SA31*T2XII+SA32*T2YII+SA33*T2ZII
              K=SA31*NXII +SA32*NYII +SA33*NZII
!
              BIGM1(ICENT)        =T1XII*A+T1YII*D+T1ZII*G
              BIGM1(ICENT+  NCOLM)=T1XII*B+T1YII*E+T1ZII*H
              BIGM1(ICENT+2*NCOLM)=T1XII*C+T1YII*F+T1ZII*K
!
              BIGM1(ICENT+3*NCOLM)=T2XII*A+T2YII*D+T2ZII*G
              BIGM1(ICENT+4*NCOLM)=T2XII*B+T2YII*E+T2ZII*H
              BIGM1(ICENT+5*NCOLM)=T2XII*C+T2YII*F+T2ZII*K
!
              BIGM1(ICENT+6*NCOLM)=NXII*A+NYII*D+NZII*G
              BIGM1(ICENT+7*NCOLM)=NXII*B+NYII*E+NZII*H
              BIGM1(ICENT+8*NCOLM)=NXII*C+NYII*F+NZII*K
!
      end do ! Was loop 5160
         ENDIF
         ENDIF
         END subroutine rbigrt
         
end module rotated_boundary_conditions_legacy
