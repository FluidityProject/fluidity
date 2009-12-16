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

module spaerr_module

  use FLDebug
  use tr2d_module

implicit none

private

public :: jacdia, jacim1, lesime, lesime2, lesime0, lesimec, eletens, get_xpctel

contains

  SUBROUTINE JACVOL(VOL,X,Y,Z,D3)
! This sub calculates the volume of an element. 
! X,Y,Z are the coords of the 4 nodes. 

    REAL VOL,X(4),Y(4),Z(4) 
    LOGICAL D3
! Local variables...
    REAL X12,X13,X14,Y12,Y13,Y14,Z12,Z13,Z14 

    IF(D3) THEN
       X12 = X(2) - X(1)
       X13 = X(3) - X(1)
       X14 = X(4) - X(1)
       Y12 = Y(2) - Y(1)
       Y13 = Y(3) - Y(1)
       Y14 = Y(4) - Y(1)
       Z12 = Z(2) - Z(1)
       Z13 = Z(3) - Z(1)
       Z14 = Z(4) - Z(1)

       VOL=ABS(X12*Y13*Z14 + X13*Y14*Z12 + X14*Y12*Z13&
     &          - X14*Y13*Z12 - X13*Y12*Z14 - X12*Y14*Z13)/6.
    ELSE
! Caculate area=VOL in 2-D. (half base x height) 
       VOL=0.5*( (X(2)-X(1))*(Y(3)-Y(1)) &
     &               -(X(3)-X(1))*(Y(2)-Y(1)) )
    ENDIF

  END SUBROUTINE JACVOL

  SUBROUTINE JACDET(DET,M,NDIM)
! This sub calculates the volume of an element. 
! X,Y,Z are the coords of the 4 nodes. 
    INTEGER NDIM
    REAL DET,M(NDIM,NDIM) 

    IF(NDIM.EQ.3) THEN
       DET = M(1,1)*( M(2,2)*M(3,3) - M(2,3)*M(3,2) )&
     &      + M(2,1)*( M(3,2)*M(1,3) - M(1,2)*M(3,3) )&
     &      + M(3,1)*( M(1,2)*M(2,3) - M(2,2)*M(1,3) )
    ELSE
       DET = M(1,1)*M(2,2)-M(1,2)*M(2,1) 
    ENDIF

  END SUBROUTINE JACDET

  SUBROUTINE JACDIA(AA,V,D,N, &
! Working arrays...
     &       A) 
! This sub performs Jacobi rotations of a symmetric matrix in order to 
! find the eigen-vectors V and the eigen values A so 
! that AA=V^T D V & D is diagonal. 
! It uses the algorithm of Matrix Computations 2nd edition, p196. 
    REAL TOLER,CONVEG
    PARAMETER(TOLER=1.E-14,CONVEG=1.E-7) 
    INTEGER N
    REAL AA(N,N),V(N,N),D(N), A(N,N)
! Local variables...
    REAL R,ABSA,MAXA,COSAL2,COSALF,SINAL2,SINALF,MAXEIG
    INTEGER ITS,NITS,Q,P,QQ,PP 

    NITS=9*(N*N-N)/2

    do P=1,N
       do Q=1,N
          V(P,Q)=0.
          A(P,Q)=AA(P,Q)
       END DO
    END DO
    
!     Check first whether matrix is diagonal
    IF(Q.EQ.0) THEN

       do PP=1,N
          D(PP) = A(PP,PP)
          do QQ=1,N
             IF(PP.EQ.QQ) THEN
                V(PP,QQ) = 1.0
             ELSE
                V(PP,QQ) = 0.0
             END IF
          END DO
       END DO
       RETURN
    END IF

    MAXEIG=0.
    do P=1,N
       V(P,P)=1.0
       MAXEIG=MAX(MAXEIG,ABS(A(P,P)))
    END DO
    IF(MAXEIG.LT.TOLER) THEN
       D(1:N) = 0.0
       GOTO 2000
    ENDIF

    do ITS=1,NITS
! Find maximum on upper diagonal of matrix. 
! QQ is the coln; PP is the row. 
       Q=0
       P=0
       MAXA=0.
       do PP=1,N-1
          do QQ=PP+1,N
             ABSA=ABS(A(PP,QQ)) 
             IF(ABSA.GT.MAXA) THEN
                MAXA=ABSA
                Q=QQ
                P=PP
             ENDIF
          END DO
       END DO

       IF(MAXA/MAXEIG.LT.CONVEG) GOTO 2000
! Rotate with (Q,P) postions.
       R=MAX(TOLER,SQRT( (A(P,P)-A(Q,Q))**2 + 4.*A(P,Q)**2 ) ) 
       IF(A(P,P).GT.A(Q,Q)) THEN
          COSAL2=0.5+0.5*(A(P,P)-A(Q,Q))/R
          COSALF=SQRT(COSAL2)
          IF(ABS(COSALF).LT.TOLER) COSALF=TOLER
          SINALF=A(Q,P)/(R*COSALF)
       ELSE
          SINAL2=0.5-0.5*(A(P,P)-A(Q,Q))/R
          SINALF=SQRT(SINAL2)
          IF(ABS(SINALF).LT.TOLER) SINALF=TOLER
          COSALF=A(Q,P)/(R*SINALF)
       ENDIF
! Pre and Post multiply of A=R^T A R  by rotation matrix. 
       CALL JACPRE(-SINALF,COSALF,P,Q,A,N)
       CALL JACPOS( SINALF,COSALF,P,Q,A,N)
! Accumulate rotations V=R^T V
       CALL JACPRE(-SINALF,COSALF,P,Q,V,N)
    end do

2000      CONTINUE 
! Put e-values in a vector...
    do Q=1,N
       D(Q)=A(Q,Q)
    END DO

  END SUBROUTINE JACDIA

  SUBROUTINE JACPRE(SINALF,COSALF,P,Q,A,N)
! This sub performs matrix-matrix multiplication A=R*A. 
! PRE-MULTIPLY matrix A by transpose of Rotation matrix 
! is realised by passing -SINALF down into SINALF. 
    INTEGER N
    REAL SINALF,COSALF,A(N,N)
    INTEGER P,Q
! Local variables...
    INTEGER I
    REAL P1I
!
! Premultiply by rotation matrix...
    do I=1,N
! Row P 1st...
       P1I   =COSALF*A(P,I)-SINALF*A(Q,I)
! Row 2nd put strait in A...
       A(Q,I)=SINALF*A(P,I)+COSALF*A(Q,I)
       A(P,I)=P1I 
    END DO

  END SUBROUTINE JACPRE

  SUBROUTINE JACPOS(SINALF,COSALF,P,Q,A,N)
! This sub performs matrix-matrix multiplication A=A*R. 
! POST-MULTIPLY matrix A by transpose of Rotation matrix 
! is realised by passing -SINALF down into SINALF. 
    INTEGER N
    REAL SINALF,COSALF,A(N,N)
    INTEGER P,Q
! Local variables...
    INTEGER I
    REAL IP1
!
! Post multiply by rotation matrix...
    do I=1,N
! Column P 1st...
       IP1   = COSALF*A(I,P)+SINALF*A(I,Q)
! column 2nd put strait in A...
       A(I,Q)=-SINALF*A(I,P)+COSALF*A(I,Q)
       A(I,P)=IP1
    END DO

  END SUBROUTINE JACPOS

  SUBROUTINE JACIM1(H1,&
     &      H2,HCOM,NDIM,&
! Working arrays...
     &      V1,V2,A,VV,DD,D1,D2,GETMIN)
! This sub combines the two metrics H1, H2 into HCOM. 
! If GETMIN then find max inner ellipsoid otherwise min-original.
    REAL SMALL
    PARAMETER(SMALL=1.E-20)
    INTEGER NDIM
    REAL H1(NDIM,NDIM),H2(NDIM,NDIM),HCOM(NDIM,NDIM)
    REAL V1(NDIM,NDIM),V2(NDIM,NDIM),A(NDIM,NDIM)
    REAL VV(NDIM,NDIM)
    REAL DD(NDIM),D1(NDIM),D2(NDIM)
    LOGICAL GETMIN
! Local variables...
    INTEGER I
    REAL RMIN1,RMAX1,RMIN2,RMAX2

! This sub superimposes 2 metrics H1 & H2 into HCOM.
    CALL JACDIA(H1,V1,D1,NDIM, &
! Working arrays...
     &       A) 
    CALL JACDIA(H2,V2,D2,NDIM, &
! Working arrays...
     &       A) 
! Choose metric which is least distorted to distort Euclidean space.
    RMIN1= 1.E+20
    RMAX1=-1.E+20
    RMIN2= 1.E+20
    RMAX2=-1.E+20

    do I=1,NDIM
       RMIN1=MIN(RMIN1,D1(I))
       RMAX1=MAX(RMAX1,D1(I))
       RMIN2=MIN(RMIN2,D2(I))
       RMAX2=MAX(RMAX2,D2(I))
    END DO
    RMAX1=MAX(RMAX1,SMALL)
    RMAX2=MAX(RMAX2,SMALL)

    IF(RMIN1/RMAX1.GT.RMIN2/RMAX2) THEN
       CALL JACIM2(V1,D1,H2,HCOM,NDIM,&
! Working arrays...
     &        A,VV,V2,DD,D2,GETMIN)
    ELSE
       CALL JACIM2(V2,D2,H1,HCOM,NDIM,&
! Working arrays...
     &        A,VV,V1,DD,D1,GETMIN)
    ENDIF

  END SUBROUTINE JACIM1

  SUBROUTINE JACIM2(V1,D1,H2,HCOM,NDIM,&
! Working arrays...
     &        A,VV,VVI,DD,DDI,GETMIN)
! This sub find e-vals & vecs of  (D^{-0.5} V M V^T D^{-0.5}). 
! HCOM is the combined Hessian. H2 is matrix M in the above equation.
    INTEGER NDIM
    REAL V1(NDIM,NDIM),D1(NDIM),H2(NDIM,NDIM),HCOM(NDIM,NDIM)
    REAL A(NDIM,NDIM),VV(NDIM,NDIM),VVI(NDIM,NDIM)
    REAL DD(NDIM),DDI(NDIM)
    LOGICAL GETMIN
! Local variables...
    INTEGER I,J,P
    REAL RSUM,RS

! DDI=D1^{-0.5}. DD=D1^{0.5}.  
    do I=1,NDIM
       RS=SQRT(D1(I))
       DDI(I)=1./RS
       DD(I) =RS
    END DO
! VVI=DDI^T V1.  V1=DD^T V1. 
    do I=1,NDIM
       do J=1,NDIM
          VVI(I,J)=DDI(I)*V1(I,J)
          V1(I,J) =DD(I) *V1(I,J)
       END DO
    END DO
! A=H2*VVI^T
    do I=1,NDIM
       do J=1,NDIM
          RSUM=0.0
          do P=1,NDIM
             RSUM=RSUM+H2(I,P)*VVI(J,P)
          END DO
          A(I,J)=RSUM
       END DO
    END DO
! H2=VVI*A
    do I=1,NDIM
       do J=1,NDIM
          RSUM=0.0
          do P=1,NDIM
             RSUM=RSUM+VVI(I,P)*A(P,J)
          END DO
          H2(I,J)=RSUM
       END DO
    END DO

! Find the e-vectors & values of H2.
    CALL JACDIA(H2,VV,DD,NDIM, &
! Working arrays...
     &                  A) 
    IF(GETMIN) THEN
       do I=1,NDIM
          DD(I)=MIN(DD(I),1.0)
       END DO
    ELSE
! The MAX was originally used here...
       do I=1,NDIM
          DD(I)=MAX(DD(I),1.0)
       END DO
    ENDIF
! H2=VV^T DD VV **********
! A=DD*VV
    do I=1,NDIM
       do J=1,NDIM
          A(I,J)=DD(I)*VV(I,J)
       END DO
    END DO
! H2=VV^T A
    do I=1,NDIM
       do J=1,NDIM
          RSUM=0.0
          do P=1,NDIM
             RSUM=RSUM+VV(P,I)*A(P,J)
          END DO
          H2(I,J)=RSUM
       END DO
    END DO

! HCOM=V1^T H2 V1 ************
! A=H2*V1 
    do I=1,NDIM
       do J=1,NDIM
          RSUM=0.0
          do P=1,NDIM
             RSUM=RSUM+H2(I,P)*V1(P,J)
          END DO
          A(I,J)=RSUM
       END DO
    END DO
! HCOM=V1^T A
    do I=1,NDIM
       do J=1,NDIM
          RSUM=0.0
          do P=1,NDIM
             RSUM=RSUM+V1(P,I)*A(P,J)
          END DO
          HCOM(I,J)=RSUM
       END DO
    END DO

  END SUBROUTINE JACIM2

  SUBROUTINE LESIME(NONODS,NDIM, FERROR,&
       &         MUPTXX,MUPTXY,MUPTXZ,&
       &         MUPTYY,MUPTYZ,MUPTZZ)
! This sub calculates the inverse of the metric contained in FERROR
! and puts it into 
! MUPTXX,MUPTYY =components of diffusivities. 
    INTEGER NONODS,NDIM
    REAL FERROR(NONODS*NDIM*NDIM)
    REAL MUPTXX(NONODS),MUPTXY(NONODS),MUPTXZ(NONODS)
    REAL MUPTYY(NONODS),MUPTYZ(NONODS),MUPTZZ(NONODS)
! Local variables...
    INTEGER NOD,IDIM,JDIM,K,NDIM2
    REAL AA(3,3),V(3,3),D(3), A(3,3)
    REAL RSUM

    NDIM2=NDIM**2

    do NOD=1,NONODS
! This sub performs Jacobi rotations of a symmetric matrix in order to 
! find the eigen-vectors V and the eigen values A so 
       CALL JACDIA(FERROR((NOD-1)*NDIM2+1),V,D,NDIM,&
! Working arrays...
              &       A)
!
! FORM V^T D^{-1} V
       do IDIM=1,NDIM
          D(IDIM)=1./D(IDIM)
       END DO
! FORM V^T D V
       do IDIM=1,NDIM
          do JDIM=1,NDIM
! AA=D V
             AA(IDIM,JDIM)=D(IDIM)*V(IDIM,JDIM)
          END DO
       END DO
! A=V^T*AA
       do IDIM=1,NDIM
          do JDIM=1,NDIM
             RSUM=0.
             do K=1,NDIM
                RSUM=RSUM+V(K,IDIM)*AA(K,JDIM)
             END DO
             A(IDIM,JDIM)=RSUM
          END DO
       END DO
       MUPTXX(NOD)=A(1,1)
       MUPTXY(NOD)=A(1,2)
       MUPTXZ(NOD)=A(1,3)
       MUPTYY(NOD)=A(2,2)
       MUPTYZ(NOD)=A(2,3)
       MUPTZZ(NOD)=A(3,3)

    end do

  END SUBROUTINE LESIME

  SUBROUTINE LESIME2(NONODS,NDIM, FERROR,&
     &         MUPTXX,MUPTXY,MUPTXZ,&
     &         MUPTYY,MUPTYZ,MUPTZZ)
! This sub calculates the inverse of the metric contained in FERROR
! and puts it into 
! MUPTXX,MUPTYY =components of diffusivities. 
! This just copies FERROR to MUPTXX etc.
    INTEGER NONODS,NDIM
    REAL FERROR(NONODS*NDIM*NDIM)
    REAL MUPTXX(NONODS),MUPTXY(NONODS),MUPTXZ(NONODS)
    REAL MUPTYY(NONODS),MUPTYZ(NONODS),MUPTZZ(NONODS)
! Local variables...
    INTEGER NOD,NDIM2

    NDIM2=NDIM**2

    do NOD=1,NONODS
       MUPTXX(NOD)=FERROR((NOD-1)*NDIM2+1)
       MUPTXY(NOD)=FERROR((NOD-1)*NDIM2+2)
       MUPTXZ(NOD)=FERROR((NOD-1)*NDIM2+3)
       MUPTYY(NOD)=FERROR((NOD-1)*NDIM2+5)
       MUPTYZ(NOD)=FERROR((NOD-1)*NDIM2+6)
       MUPTZZ(NOD)=FERROR((NOD-1)*NDIM2+9)
    end do

  END SUBROUTINE LESIME2

  SUBROUTINE LESIME0(NONODS,NDIM, FERROR,&
     &         MUPTXX,MUPTXY,MUPTXZ,&
     &         MUPTYY,MUPTYZ,MUPTZZ)
! This sub puts an inhomogenious isotropic length scale into the 
! LES model. 
! This sub calculates the inverse of the metric contained in FERROR
! and puts it into 
! MUPTXX,MUPTYY =components of diffusivities. 
    INTEGER NONODS,NDIM
    REAL FERROR(NONODS*NDIM*NDIM)
    REAL MUPTXX(NONODS),MUPTXY(NONODS),MUPTXZ(NONODS)
    REAL MUPTYY(NONODS),MUPTYZ(NONODS),MUPTZZ(NONODS)
! Local variables...
    INTEGER NOD,IDIM,NDIM2
    REAL V(3,3),D(3),A(3,3)
    LOGICAL PRISCR
    REAL DELSUM,ISOL

    NDIM2=NDIM**2

    do NOD=1,NONODS
! This sub performs Jacobi rotations of a symmetric matrix in order to 
! find the eigen-vectors V and the eigen values A so 
       CALL JACDIA(FERROR((NOD-1)*NDIM2+1),V,D,NDIM,&
! Working arrays...
              &       A)

! Find length scales in each direction
       do IDIM=1,NDIM
          D(IDIM)=1./D(IDIM)
       END DO

! Calculate a single characteristic length for the element
       DELSUM=0.
       do IDIM=1,NDIM
          DELSUM=DELSUM+D(IDIM)
       END DO
       ISOL=DELSUM/3.0

       MUPTXX(NOD)=ISOL
       MUPTXY(NOD)=ISOL
       MUPTXZ(NOD)=ISOL
       MUPTYY(NOD)=ISOL
       MUPTYZ(NOD)=ISOL
       MUPTZZ(NOD)=ISOL

    end do

  END SUBROUTINE LESIME0

  SUBROUTINE LESIMEC(NONODS,NDIM,&
     &         MUPTXX,MUPTXY,MUPTXZ,&
     &         MUPTYY,MUPTYZ,MUPTZZ)
! This sub assigns a fixed, constant length scale to MUPTXX etc 
! to be used in the LES model.

    INTEGER NONODS,NDIM
    REAL MUPTXX(NONODS),MUPTXY(NONODS),MUPTXZ(NONODS)
    REAL MUPTYY(NONODS),MUPTYZ(NONODS),MUPTZZ(NONODS)
! Local variables...
    INTEGER NOD
    REAL FIXEDL,FIXEDL2
      
    FIXEDL=0.02
    FIXEDL2=FIXEDL*FIXEDL

    do NOD=1,NONODS
       MUPTXX(NOD)=FIXEDL2
       MUPTXY(NOD)=FIXEDL2
       MUPTXZ(NOD)=FIXEDL2
       MUPTYY(NOD)=FIXEDL2
       MUPTYZ(NOD)=FIXEDL2
       MUPTZZ(NOD)=FIXEDL2
    end do

  END SUBROUTINE LESIMEC

!     Given a metric the function returns the number of tetrahedral
!     elements we would expect after adapting to the metric
  INTEGER FUNCTION GET_XPCTEL(Metric, X, Y, Z, NDGLNO, NNodes, &
     &     NTetra, NLOC)

    INTEGER NNodes
    REAL Metric(NNodes*9), X(NNodes), Y(NNodes), Z(NNodes)
    INTEGER NTetra, NLOC, NDGLNO(NTetra*NLOC)
      
!     VOLSCA=the quantity to scale the predicted no of elements 
    REAL VOLSCA
      
!     volume of unit tetrahedral(VOL1TE) =1./sqrt(72). (each side has length 1.)
    REAL VOL1EL
      
    PARAMETER(VOLSCA=1.2, VOL1EL=0.11785113)
    REAL VOLUME, VOL, DET, MEANM(9), XX(4), YY(4), ZZ(4), MXDET, SMDET
    INTEGER ELE, Node, I, ILOC
      
    VOLUME=0.0
    MXDET = -1E+30
    SMDET = 0.0
    do ELE=1,NTetra
       do I=1,9
          MEANM(I) = 0.0
       END DO
       do ILOC=1,4
          Node = NDGLNO((ELE-1)*NLOC+ILOC)
          XX(ILOC) = X(Node)
          YY(ILOC) = Y(Node)
          ZZ(ILOC) = Z(Node)
          do I=1,9
             MEANM(I) = MEANM(I) + Metric((Node-1)*9+I)*0.25
          END DO
       END DO
         
       CALL JACVOL(VOL, XX, YY, ZZ, .TRUE.)
       if (vol <= 0.0) then
          ewrite(-1,*) "Volume of element", ele, "not positive!"
          ewrite(-1,*) "Bad news! Your inputs to adaptivity are fecked."
       end if
       CALL JACDET(DET, MEANM, 3)
       MXDET = MAX(MXDET,DET)
       SMDET = SMDET + DET
         
       VOLUME = VOLUME + VOL*SQRT(MAX(0.0,DET))
    END DO
      
    I = INT(VOLSCA*VOLUME/VOL1EL)
    IF( I .LT. 2 ) THEN
       ewrite(0,*)&
            &    '+++ Warning: GET_XPCTEL looks wrong: ',I
       ewrite(0,*)&
            &    '  VOLSCA,VOLUME,VOL1EL: ',VOLSCA,VOLUME,VOL1EL
       ewrite(0,*)&
            &    '  Determinant max & sum: ',MXDET,SMDET
       ewrite(0,*)&
            &    '  NNOD,NELM,NLOC: ',NNodes,NTetra,NLOC
    END IF
    GET_XPCTEL = I
      
  END FUNCTION GET_XPCTEL

  subroutine eletens(x,y,z, &
    & tensxx,tensxy,tensxz,tensyy,tensyz,tenszz, &
    & fredop,nonods,xnonod,totele,nloc,disopt,ndglno)
!     This sub calculates the NODE-WISE TENSOR TENS?? THAT 
!     REPRESENTS THE SIZE AND SHAPE OF THE SURROUNDING ELEMENTS.
!     DISOPT=LES option.
    integer ngi
    parameter(ngi=4)
    integer fredop,nonods,xnonod,totele,nloc,disopt
    real tensxx(fredop),tensxy(fredop),tensxz(fredop)
    real tensyy(fredop),tensyz(fredop),tenszz(fredop)
    real x(xnonod),y(xnonod),z(xnonod)
    integer ndglno(totele*nloc)
    !     Local variables...
    real, dimension(:), allocatable :: n,nlx,nly,nlz,weight
    real, dimension(:), allocatable :: detwei,m
    real, dimension(:), allocatable :: udl,vdl,wdl,gamma,pml
    
    ! Dummy argument used to replace passing of rmem(rpt)
    real, dimension(0) :: unused_real_arg
    
    allocate(n(nloc * ngi))
    allocate(nlx(nloc * ngi))
    allocate(nly(nloc * ngi))
    allocate(nlz(nloc * ngi))
    allocate(weight(ngi))
    allocate(detwei(ngi))
    allocate(m(nloc * ngi))
      
    call tr3d(.false.,ngi,nloc,nloc,&
   &     m,weight,n,nlx,nly,nlz,&
   &     0,0,unused_real_arg,unused_real_arg,unused_real_arg,unused_real_arg)
    
    allocate(udl(nloc * nloc * ngi))
    allocate(vdl(nloc * nloc * ngi))
    allocate(wdl(nloc * nloc * ngi))
    allocate(gamma(nloc * nloc * ngi))
      
    allocate(pml(fredop))      
    
    call eletens2(x,y,z,&
   &     tensxx,tensxy,tensxz,tensyy,tensyz,tenszz, pml,&
   &     fredop,totele,disopt,&
   &     n,nlx,nly,nlz, &
   &     n,weight,ngi,nloc,nloc, &
   &     ndglno,ndglno,nonods,&
   &     detwei,&
   &     udl,vdl,wdl,&
   &     gamma)

    deallocate(pml)
    
    deallocate(udl)
    deallocate(vdl)
    deallocate(wdl)
    deallocate(gamma)

    deallocate(n)
    deallocate(nlx)
    deallocate(nly)
    deallocate(nlz)
    deallocate(weight)
    deallocate(detwei)
    deallocate(m)

  end subroutine eletens

  SUBROUTINE ELETENS2(X,Y,Z,&
     &     TENSXX,TENSXY,TENSXZ,TENSYY,TENSYZ,TENSZZ, PML,&
     &     FREDOP,TOTELE,DISOPT,&
     &     N,NLX,NLY,NLZ, &
     &     M,WEIGHT,NGI,NLOC,MLOC, &
     &     PNDGLN,XONDGL,XNONOD,&
     &     DETWEI,&
     &     UDL,VDL,WDL,&
     &     GAMMA  )
    use sml
!     This sub calculates the NODE-WISE TENSOR TENS?? THAT 
!     REPRESENTS THE SIZE AND SHAPE OF THE SURROUNDING ELEMENTS.
!     DISOPT=LES option.
    INTEGER XNONOD
    INTEGER MLOC,NLOC,NGI
    INTEGER TOTELE,FREDOP,DISOPT
!     
    REAL TENSXX(FREDOP),TENSXY(FREDOP),TENSXZ(FREDOP)
    REAL TENSYY(FREDOP),TENSYZ(FREDOP),TENSZZ(FREDOP)
    REAL PML(FREDOP)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    INTEGER PNDGLN(TOTELE*MLOC)
    INTEGER XONDGL(TOTELE*NLOC)
!     
    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL M(MLOC,NGI)
    REAL WEIGHT(NGI)
!     
    REAL DETWEI(NGI)
    REAL UDL(NLOC*NLOC,NGI),VDL(NLOC*NLOC,NGI)
    REAL WDL(NLOC*NLOC,NGI)
    REAL GAMMA(NLOC*NLOC,NGI)
!     REAL MX(MLOC,NGI),MY(MLOC,NGI),MZ(MLOC,NGI)
!     HX,HY-characteristic length scales in x,y directions.
!     Local variables...
    REAL RN
    INTEGER NDIM
    PARAMETER(NDIM=3)
    REAL AA(NDIM,NDIM),V(NDIM,NDIM),D(NDIM),A(NDIM,NDIM)
!     
    INTEGER ELE,ILOC,L,L1,L2,IGLX1,IGLX2,IGLX,ID,NID
    INTEGER GI,GLOBI
!     
    REAL HXGI,HYGI,HZGI,DETJ
    REAL HOVERQ
    REAL AGI,BGI,CGI,DGI
    REAL EGI,FGI,GGI,HGI,KGI
    REAL A11,A12,A13, A21,A22,A23, A31,A32,A33
    REAL RWIND
    REAL RFACT,RT1,RT2,RT3,D1,D2,D3
    ewrite(1, *) "SUBROUTINE ELETENS2()"

    RWIND =1./REAL(6)
    NID=NLOC*NLOC

    TENSXX(1:FREDOP) = 0.0
    TENSXY(1:FREDOP) = 0.0
    TENSXZ(1:FREDOP) = 0.0
    TENSYY(1:FREDOP) = 0.0
    TENSYZ(1:FREDOP) = 0.0
    TENSZZ(1:FREDOP) = 0.0
    PML(1:FREDOP) = 0.0

!     This subroutine forms a contabution to the Right Hand Side
!     of Poissons pressure equation, as well as  F1 & F2.
!     
    do ELE=1,TOTELE
       do GI=1,NGI

          AGI=0.
          BGI=0.
          CGI=0.

          DGI=0.
          EGI=0.
          FGI=0.

          GGI=0.
          HGI=0.
          KGI=0.

          do L=1,NLOC
             IGLX=XONDGL((ELE-1)*NLOC+L)
             AGI=AGI+NLX(L,GI)*X(IGLX) 
             BGI=BGI+NLX(L,GI)*Y(IGLX) 
             CGI=CGI+NLX(L,GI)*Z(IGLX) 

             DGI=DGI+NLY(L,GI)*X(IGLX) 
             EGI=EGI+NLY(L,GI)*Y(IGLX) 
             FGI=FGI+NLY(L,GI)*Z(IGLX) 

             GGI=GGI+NLZ(L,GI)*X(IGLX) 
             HGI=HGI+NLZ(L,GI)*Y(IGLX) 
             KGI=KGI+NLZ(L,GI)*Z(IGLX) 

          end do

          DETJ=AGI*(EGI*KGI-FGI*HGI)&
     &           -BGI*(DGI*KGI-FGI*GGI)&
     &           +CGI*(DGI*HGI-EGI*GGI)
          DETWEI(GI)=abs(DETJ)*WEIGHT(GI)

!     For coefficient in the inverse mat of the jacobian. 
          A11= (EGI*KGI-FGI*HGI) /DETJ
          A21=-(DGI*KGI-FGI*GGI) /DETJ
          A31= (DGI*HGI-EGI*GGI) /DETJ

          A12=-(BGI*KGI-CGI*HGI) /DETJ
          A22= (AGI*KGI-CGI*GGI) /DETJ
          A32=-(AGI*HGI-BGI*GGI) /DETJ

          A13= (BGI*FGI-CGI*EGI) /DETJ
          A23=-(AGI*FGI-CGI*DGI) /DETJ
          A33= (AGI*EGI-BGI*DGI) /DETJ

          RWIND =1./REAL(6)
          NID=NLOC*NLOC
!     **********calculate normalised velocitys across element...          
          ID=0
          do L1=1,NLOC
             IGLX1=XONDGL((ELE-1)*NLOC+L1) 
             do L2=1,NLOC
                IGLX2=XONDGL((ELE-1)*NLOC+L2) 
                ID=ID+1
                if(l1.eq.l2) then
                   UDL(ID,GI)=0.0
                   VDL(ID,GI)=0.0
                   WDL(ID,GI)=0.0
                   GAMMA(ID,GI)=0.0
                else
                   UDL(ID,GI)=X(IGLX1)-X(IGLX2)
                   VDL(ID,GI)=Y(IGLX1)-Y(IGLX2)
                   WDL(ID,GI)=Z(IGLX1)-Z(IGLX2)
!     Normalise 
                   RN=SQRT(UDL(ID,GI)**2+VDL(ID,GI)**2+WDL(ID,GI)**2)
                   UDL(ID,GI)=UDL(ID,GI)/RN
                   VDL(ID,GI)=VDL(ID,GI)/RN
                   WDL(ID,GI)=WDL(ID,GI)/RN

!     put the length scale in
                   HXGI=ABS(A11*UDL(ID,GI)+A12*VDL(ID,GI)+A13*WDL(ID,GI))
                   HYGI=ABS(A21*UDL(ID,GI)+A22*VDL(ID,GI)+A23*WDL(ID,GI))
                   HZGI=ABS(A31*UDL(ID,GI)+A32*VDL(ID,GI)+A33*WDL(ID,GI))
!     HX,HY are the characteristic length scales in x,y directions. 
                   HOVERQ=RN
                   GAMMA(ID,GI)=RWIND*HOVERQ 
                endif
             END DO
          END DO
!     **********calculate normalised velocitys across element...            

       end do

       do ILOC=1,MLOC
          GLOBI=PNDGLN((ELE-1)*MLOC+ILOC)

          do GI=1,NGI
             do ID=1,NID
!     TENS=GAMMA(ID,GI) Q Q^T
                RFACT=N(ILOC,GI)*DETWEI(GI)*GAMMA(ID,GI)
                  
                TENSXX(GLOBI)=TENSXX(GLOBI) + RFACT*UDL(ID,GI)*UDL(ID,GI)
                TENSXY(GLOBI)=TENSXY(GLOBI) + RFACT*UDL(ID,GI)*VDL(ID,GI)
                TENSXZ(GLOBI)=TENSXZ(GLOBI) + RFACT*UDL(ID,GI)*WDL(ID,GI)
                TENSYY(GLOBI)=TENSYY(GLOBI) + RFACT*VDL(ID,GI)*VDL(ID,GI)
                TENSYZ(GLOBI)=TENSYZ(GLOBI) + RFACT*VDL(ID,GI)*WDL(ID,GI)
                TENSZZ(GLOBI)=TENSZZ(GLOBI) + RFACT*WDL(ID,GI)*WDL(ID,GI)
                  
!     USE THE COMPONENT OF DIFLIN THE X,Y & Z-DIRECTIONS 
!     RESPECTIVELY FOR C1T,C2T,C3T.
             end do
             PML(GLOBI) =PML(GLOBI) +M(ILOC,GI)*DETWEI(GI) 
          end do
       end do
    end do

    do GLOBI=1,FREDOP
       TENSXX(GLOBI)=TENSXX(GLOBI)/PML(GLOBI)
       TENSXY(GLOBI)=TENSXY(GLOBI)/PML(GLOBI)
       TENSXZ(GLOBI)=TENSXZ(GLOBI)/PML(GLOBI)
       TENSYY(GLOBI)=TENSYY(GLOBI)/PML(GLOBI)
       TENSYZ(GLOBI)=TENSYZ(GLOBI)/PML(GLOBI)
       TENSZZ(GLOBI)=TENSZZ(GLOBI)/PML(GLOBI)

!     nb we want 1/L^2 - at the moment we have L on the diagonal.
!     Make sure the eigen-values are positive...
       AA(1,1)=TENSXX(GLOBI)
       AA(1,2)=TENSXY(GLOBI)
       AA(1,3)=TENSXZ(GLOBI)
         
       AA(2,1)=TENSXY(GLOBI)
       AA(2,2)=TENSYY(GLOBI)
       AA(2,3)=TENSYZ(GLOBI)
         
       AA(3,1)=TENSXZ(GLOBI)
       AA(3,2)=TENSYZ(GLOBI)
       AA(3,3)=TENSZZ(GLOBI)
       CALL JACDIA(AA,V,D,NDIM,A)
!     TENSOR=V^T 1/D**2 V

       IF((DISOPT.GE.45).AND.(DISOPT.LE.48)) THEN
!     SET to metric which has 1/h^2 in it...
          D1=1./MAX(1.E-16,D(1)**2)
          D2=1./MAX(1.E-16,D(2)**2)
          D3=1./MAX(1.E-16,D(3)**2)
       ELSE IF((DISOPT.EQ.43).OR.(DISOPT.EQ.44)) THEN
!     set to inverse of metric which is a multiple of the tensor
          D1=MAX(1.E-16,D(1)**2)
          D2=MAX(1.E-16,D(2)**2)
          D3=MAX(1.E-16,D(3)**2)
       ELSE
          FLAbort("NOT A VALID OPTION FOR LES ASSEMBLED EQNS")
       ENDIF
       CALL SMLMA3( &
     &        TENSXX(GLOBI),TENSXY(GLOBI),TENSXZ(GLOBI),&
     &        RT1,          TENSYY(GLOBI),TENSYZ(GLOBI),&
     &        RT2,          RT3,          TENSZZ(GLOBI), &
!     =
     &        V(1,1),   V(2,1),   V(3,1),  &
     &        V(1,2),   V(2,2),   V(3,2),  &
     &        V(1,3),   V(2,3),   V(3,3),             &
!     *          
     &        D1*V(1,1),   D1*V(1,2),   D1*V(1,3),  &
     &        D2*V(2,1),   D2*V(2,2),   D2*V(2,3),  &
     &        D3*V(3,1),   D3*V(3,2),   D3*V(3,3) ) 
    END DO

  end subroutine eletens2

end module spaerr_module
