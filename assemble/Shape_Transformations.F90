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

module Shape_Transformations
  ! This module contains all kinds of calculations related to the tranformation
  ! of shape functions to physical space.
  use Element_Numbering, only: ilink2, silink2
  use AllSorts
  use sml
  implicit none

contains

  SUBROUTINE CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MLOC,NGI, &
       N,NLX,NLY,NLZ, M,MLX,MLY,MLZ, WEIGHT, DETWEI,D3,DCYL, &
       NX,NY,NZ, MX,MY,MZ, &
       A11,A12,A13, A21,A22,A23, A31,A32,A33, &
       XD,YD,ZD, &
       ISPHERE) 
    INTEGER, intent(in):: ELE,TOTELE,XNONOD,NLOC,MLOC,NGI,ISPHERE
    REAL, intent(in):: X(XNONOD),Y(XNONOD),Z(XNONOD)
    INTEGER, intent(in):: XONDGL(TOTELE*NLOC)
    !     Volume shape functions...
    REAL, intent(in):: N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL, intent(in):: M(MLOC,NGI),MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
    REAL, intent(in):: WEIGHT(NGI)
    REAL, intent(out):: DETWEI(NGI)
    REAL, intent(out):: NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    REAL, intent(out):: MX(MLOC,NGI),MY(MLOC,NGI),MZ(MLOC,NGI)
    LOGICAL, intent(in):: D3,DCYL
    REAL, intent(out):: A11(NGI),A12(NGI),A13(NGI)
    REAL, intent(out):: A21(NGI),A22(NGI),A23(NGI)
    REAL, intent(out):: A31(NGI),A32(NGI),A33(NGI) 
    REAL, intent(out):: XD(NGI),YD(NGI),ZD(NGI)
    !     Local variables...
    INTEGER IGL,GI,ILOC
    !     Local coords...
    INTEGER NCLOC

    IF(ISPHERE.EQ.0) THEN
       CALL DETNLXMAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MLOC,NGI, &
            NLX,NLY,NLZ,MLX,MLY,MLZ, WEIGHT, DETWEI, &
            NX,NY,NZ, MX,MY,MZ, &
            A11,A12,A13, A21,A22,A23, A31,A32,A33)
       DO GI=1,NGI
          XD(GI)=0.0
          YD(GI)=0.0
          ZD(GI)=0.0
          DO ILOC=1,NLOC
             IGL=XONDGL((ELE-1)*NLOC+ILOC)
             XD(GI)=XD(GI)+N(ILOC,GI)*X(IGL)
             YD(GI)=YD(GI)+N(ILOC,GI)*Y(IGL)
             ZD(GI)=ZD(GI)+N(ILOC,GI)*Z(IGL)
          END DO
       END DO
    ELSE
       NCLOC=10

       CALL CALHIGHNOD2(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MLOC,NGI, &
            N,NLX,NLY,NLZ,MLX,MLY,MLZ, WEIGHT, DETWEI,D3,DCYL, &
            NX,NY,NZ, MX,MY,MZ, &
            NCLOC,M,MLX,MLY,MLZ,&
            A11,A12,A13, A21,A22,A23, A31,A32,A33,&
            XD,YD,ZD)
    ENDIF
  END SUBROUTINE CALHIGHNODAI

  SUBROUTINE CALHIGHNOD2(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MLOC,NGI, &
       N,NLX,NLY,NLZ,MLX,MLY,MLZ, WEIGHT, DETWEI,D3,DCYL, &
       NX,NY,NZ, MX,MY,MZ, &
       NCLOC,NC,NCLX,NCLY,NCLZ,&
       A11,A12,A13, A21,A22,A23, A31,A32,A33,&
       XD,YD,ZD)
    ! This sub calculates the nodes for high order tets of degree ISPHERE 
    INTEGER, intent(in):: ELE,TOTELE,XNONOD,NLOC,MLOC,NGI,NCLOC
    REAL, intent(in):: X(XNONOD),Y(XNONOD),Z(XNONOD)
    INTEGER, intent(in):: XONDGL(TOTELE*NLOC)
    !     Volume shape functions...
    REAL, intent(in)::N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL, intent(in)::MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
    REAL, intent(in):: NC(NCLOC,NGI),NCLX(NCLOC,NGI),NCLY(NCLOC,NGI),NCLZ(NCLOC,NGI)
    REAL, intent(in):: WEIGHT(NGI)
    REAL, intent(out):: DETWEI(NGI)
    REAL, intent(out):: NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    REAL, intent(out):: MX(MLOC,NGI),MY(MLOC,NGI),MZ(MLOC,NGI)
    REAL, intent(out):: A11(NGI),A12(NGI),A13(NGI)
    REAL, intent(out):: A21(NGI),A22(NGI),A23(NGI)
    REAL, intent(out):: A31(NGI),A32(NGI),A33(NGI) 
    REAL, intent(out):: XD(NGI),YD(NGI),ZD(NGI)
    LOGICAL, intent(in):: D3,DCYL
    !     Local variables...
    INTEGER ILOC,GI,JGLX, JLOC
    REAL MAT(NCLOC,NCLOC)
    REAL BX(NCLOC),BY(NCLOC),BZ(NCLOC),BR(NCLOC)
    REAL XL(NCLOC),YL(NCLOC),ZL(NCLOC),RL(NCLOC)
    REAL RAD,NCN,VOLUME

    IF(NCLOC.EQ.10) THEN
       !     Calculate initial XL,YL,ZL for a 10 node tet in a fast way.
       CALL UNWIND_XL_CALC(XL,YL,ZL,NCLOC,XNONOD,NLOC,TOTELE,XONDGL,&
            ELE,X,Y,Z)

    ELSE
       CALL DETNNN(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NGI, &
            N,NLX,NLY,NLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL) 
       ! Solve matrix eqn for intermedia coords.
       DO ILOC=1,NCLOC
          DO JLOC=1,NCLOC
             MAT(ILOC,JLOC)=0.0
             DO GI=1,NGI
                MAT(ILOC,JLOC)=MAT(ILOC,JLOC)+DETWEI(GI)*NC(ILOC,GI)*NC(JLOC,GI)
             END DO
          END DO

          BX(ILOC)=0.0
          BY(ILOC)=0.0
          BZ(ILOC)=0.0
          BR(ILOC)=0.0
          DO JLOC=1,NLOC
             JGLX=XONDGL((ELE-1)*NLOC+JLOC)
             NCN=0.0
             DO GI=1,NGI
                NCN=NCN+DETWEI(GI)*NC(ILOC,GI)*N(JLOC,GI)
             END DO
             BX(ILOC)=BX(ILOC)+NCN*X(JGLX)
             BY(ILOC)=BY(ILOC)+NCN*Y(JGLX)
             BZ(ILOC)=BZ(ILOC)+NCN*Z(JGLX)
             BR(ILOC)=BR(ILOC)+NCN*SQRT( X(JGLX)**2+Y(JGLX)**2+Z(JGLX)**2 )
          END DO
       END DO
       !         
       ! Solve matrix eqn for radii at all nodes from the radii at the 
       ! corder nodes. 
       ! Solve MAT XL=BX (NB BDIAG is overwritten)
       CALL SMLINNGOT(MAT,XL,BX,NCLOC,NCLOC,.FALSE.)
       CALL SMLINNGOT(MAT,YL,BY,NCLOC,NCLOC,.TRUE.)
       CALL SMLINNGOT(MAT,ZL,BZ,NCLOC,NCLOC,.TRUE.)
       CALL SMLINNGOT(MAT,RL,BR,NCLOC,NCLOC,.TRUE.)
       ! Adjust the coords so we have a linear variation in radius
       DO ILOC=1,NCLOC
          RAD=SQRT(XL(ILOC)**2+YL(ILOC)**2+ZL(ILOC)**2)
          XL(ILOC)=XL(ILOC)*RL(ILOC)/RAD
          YL(ILOC)=YL(ILOC)*RL(ILOC)/RAD
          ZL(ILOC)=ZL(ILOC)*RL(ILOC)/RAD
       END DO
    ENDIF


    CALL DETNLXMSUP(NLOC,MLOC,NGI, &
         NLX,NLY,NLZ, MLX,MLY,MLZ, WEIGHT, DETWEI, &
         NX,NY,NZ, MX,MY,MZ, &
         NCLOC,NCLX,NCLY,NCLZ,XL,YL,ZL,&
         A11,A12,A13, A21,A22,A23, A31,A32,A33)
    DO GI=1,NGI
       XD(GI)=0.0
       YD(GI)=0.0
       ZD(GI)=0.0
       DO ILOC=1,NCLOC
          XD(GI)=XD(GI)+NC(ILOC,GI)*XL(ILOC)
          YD(GI)=YD(GI)+NC(ILOC,GI)*YL(ILOC)
          ZD(GI)=ZD(GI)+NC(ILOC,GI)*ZL(ILOC)
       END DO
    END DO
    RETURN
  END SUBROUTINE CALHIGHNOD2

  SUBROUTINE CALHIGHNODAINN(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MLOC,NGI, &
       N,NLX,NLY,NLZ,MLX,MLY,MLZ, WEIGHT, DETWEI,D3,DCYL, &
       ISPHERE) 
    !     This sub calculates the nodes for high order tets of degree ISPHERE 
    INTEGER ELE,TOTELE,XNONOD,NLOC,MLOC,NGI,ISPHERE
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    INTEGER XONDGL(TOTELE*NLOC)
    !     Volume shape functions...
    REAL N(NLOC,NGI), NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
    REAL WEIGHT(NGI)
    REAL DETWEI(NGI)
    LOGICAL D3,DCYL
    !     Local variables...
    REAL VOLUME
    !     Local coords...
    INTEGER NCLOC

    IF(ISPHERE.EQ.0) THEN
       CALL DETNNN(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NGI,  &
            N,NLX,NLY,NLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL)
    ELSE
       NCLOC=10

       CALL CALHIGHNOD2NN(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NGI, &
            WEIGHT,DETWEI, &
            NCLOC,MLX,MLY,MLZ)
    ENDIF
  END SUBROUTINE CALHIGHNODAINN

  SUBROUTINE CALHIGHNOD2NN(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NGI, &
       WEIGHT, DETWEI, &
       NCLOC,NCLX,NCLY,NCLZ)
    ! This sub calculates the nodes for high order tets of degree ISPHERE 
    INTEGER ELE,TOTELE,XNONOD,NLOC,NGI,NCLOC
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    INTEGER XONDGL(TOTELE*NLOC)
    !     Volume shape functions...
    REAL NCLX(NCLOC,NGI),NCLY(NCLOC,NGI),NCLZ(NCLOC,NGI)
    REAL WEIGHT(NGI)
    REAL DETWEI(NGI)
    !     Local variables...
    REAL XL(NCLOC),YL(NCLOC),ZL(NCLOC)

    !     Calculate initial XL,YL,ZL for a 10 node tet in a fast way.
    CALL UNWIND_XL_CALC(XL,YL,ZL,NCLOC,XNONOD,NLOC,TOTELE,XONDGL,&
         ELE,X,Y,Z)

    CALL DETNLXMSUPNN(NGI, &
         WEIGHT,DETWEI, &
         NCLOC,NCLX,NCLY,NCLZ,XL,YL,ZL)

  END SUBROUTINE CALHIGHNOD2NN

  SUBROUTINE DETNLXMSUPNN(NGI, &
       WEIGHT,DETWEI, &
       NCLOC,NCLX,NCLY,NCLZ,XL,YL,ZL)
    ! Calculate DETWEI,NX,NY,NZ, MX,MY,MZ for element ELE
    ! For coefficient in the inverse mat of the jacobian. 
    ! This works for a superparametric fem NC,NCLX,NCLY,NCLZ are basis functions 
    ! associated with coordinates.
    ! (XL,YL,ZL) are th local coordinates of the nodes.
    INTEGER NGI,NCLOC
    ! Volume shape functions...
    REAL NCLX(NCLOC,NGI),NCLY(NCLOC,NGI),NCLZ(NCLOC,NGI)
    REAL WEIGHT(NGI)
    REAL DETWEI(NGI)
    REAL XL(NCLOC),YL(NCLOC),ZL(NCLOC)
    ! Local variables...
    REAL AGI,BGI,CGI, DGI,EGI,FGI, GGI,HGI,KGI
    REAL DETJ,VOLUME
    INTEGER GI,L

    VOLUME=0.
    DO GI=1,NGI

       AGI=0.
       BGI=0.
       CGI=0.

       DGI=0.
       EGI=0.
       FGI=0.

       GGI=0.
       HGI=0.
       KGI=0.

       DO L=1,NCLOC

          AGI=AGI+NCLX(L,GI)*XL(L) 
          BGI=BGI+NCLX(L,GI)*YL(L) 
          CGI=CGI+NCLX(L,GI)*ZL(L) 

          DGI=DGI+NCLY(L,GI)*XL(L) 
          EGI=EGI+NCLY(L,GI)*YL(L) 
          FGI=FGI+NCLY(L,GI)*ZL(L) 

          GGI=GGI+NCLZ(L,GI)*XL(L) 
          HGI=HGI+NCLZ(L,GI)*YL(L) 
          KGI=KGI+NCLZ(L,GI)*ZL(L) 

       end do

       DETJ=AGI*(EGI*KGI-FGI*HGI) &
            -BGI*(DGI*KGI-FGI*GGI) &
            +CGI*(DGI*HGI-EGI*GGI)
       DETWEI(GI)=abs(DETJ)*WEIGHT(GI)
       VOLUME=VOLUME+DETWEI(GI)

    end do

  END SUBROUTINE DETNLXMSUPNN



  SUBROUTINE UNWIND_XL_CALC(XL,YL,ZL,NCLOC,XNONOD,NLOC,TOTELE,XONDGL, &
       ELE,X,Y,Z)
    IMPLICIT NONE
    ! Calculate initial XL,YL,ZL for a 10 node tet in a fast way.
    ! If TEST print out the original outputs of this sub...
    logical, parameter:: test=.false.
    INTEGER, intent(in):: NCLOC,XNONOD,NLOC,TOTELE
    INTEGER, intent(in):: XONDGL(NLOC*TOTELE)
    REAL, intent(out):: XL(NCLOC),YL(NCLOC),ZL(NCLOC)
    INTEGER, intent(in):: ELE
    REAL, intent(in):: X(XNONOD),Y(XNONOD),Z(XNONOD)
    ! Local variables...
    INTEGER ILOC,JLOC,KLOC,IGL,JGL
    REAL RL,RAD
    REAL RADLOC(4)
    ! Function...
    DO ILOC=1,NLOC
       IGL=XONDGL((ELE-1)*NLOC+ILOC)
       RADLOC(ILOC)=SQRT(X(IGL)**2+Y(IGL)**2+Z(IGL)**2)
    END DO

    if(test) then
       DO ILOC=1,NLOC
          IGL=XONDGL((ELE-1)*NLOC+ILOC)
          ! the corner node...
          JLOC=ILOC
          KLOC=ILINK2(ILOC,JLOC)
          XL(KLOC)=X(IGL)
          YL(KLOC)=Y(IGL)
          ZL(KLOC)=Z(IGL)
          ! the mid side nodes...
          DO JLOC=ILOC+1,NLOC
             JGL=XONDGL((ELE-1)*NLOC+JLOC)
             KLOC=ILINK2(ILOC,JLOC)
             XL(KLOC)=0.5*(X(IGL)+X(JGL))
             YL(KLOC)=0.5*(Y(IGL)+Y(JGL))
             ZL(KLOC)=0.5*(Z(IGL)+Z(JGL))
             ! Adjust the coords so we have a linear variation in radius.
             RL=0.5*(RADLOC(ILOC)+RADLOC(JLOC))
             RAD=SQRT(XL(KLOC)**2+YL(KLOC)**2+ZL(KLOC)**2)
             XL(KLOC)=XL(KLOC)*RL/RAD
             YL(KLOC)=YL(KLOC)*RL/RAD
             ZL(KLOC)=ZL(KLOC)*RL/RAD
          END DO
       END DO
    endif
    if(.not.test) then
       iloc=1
       IGL=XONDGL((ELE-1)*NLOC+ILOC)
       ! The corner node iloc...
       kloc=1
       XL(KLOC)=X(IGL)
       YL(KLOC)=Y(IGL)
       ZL(KLOC)=Z(IGL)
       jloc=2
       JGL=XONDGL((ELE-1)*NLOC+JLOC)
       kloc=2
       XL(KLOC)=0.5*(X(IGL)+X(JGL))
       YL(KLOC)=0.5*(Y(IGL)+Y(JGL))
       ZL(KLOC)=0.5*(Z(IGL)+Z(JGL))
       ! Adjust the coords so we have a linear variation in radius.
       RL=0.5*(RADLOC(ILOC)+RADLOC(JLOC))
       RAD=SQRT(XL(KLOC)**2+YL(KLOC)**2+ZL(KLOC)**2)
       XL(KLOC)=XL(KLOC)*RL/RAD
       YL(KLOC)=YL(KLOC)*RL/RAD
       ZL(KLOC)=ZL(KLOC)*RL/RAD
       jloc=3
       JGL=XONDGL((ELE-1)*NLOC+JLOC)
       kloc=4
       XL(KLOC)=0.5*(X(IGL)+X(JGL))
       YL(KLOC)=0.5*(Y(IGL)+Y(JGL))
       ZL(KLOC)=0.5*(Z(IGL)+Z(JGL))
       ! Adjust the coords so we have a linear variation in radius.
       RL=0.5*(RADLOC(ILOC)+RADLOC(JLOC))
       RAD=SQRT(XL(KLOC)**2+YL(KLOC)**2+ZL(KLOC)**2)
       XL(KLOC)=XL(KLOC)*RL/RAD
       YL(KLOC)=YL(KLOC)*RL/RAD
       ZL(KLOC)=ZL(KLOC)*RL/RAD
       jloc=4
       JGL=XONDGL((ELE-1)*NLOC+JLOC)
       kloc=7
       XL(KLOC)=0.5*(X(IGL)+X(JGL))
       YL(KLOC)=0.5*(Y(IGL)+Y(JGL))
       ZL(KLOC)=0.5*(Z(IGL)+Z(JGL))
       ! Adjust the coords so we have a linear variation in radius.
       RL=0.5*(RADLOC(ILOC)+RADLOC(JLOC))
       RAD=SQRT(XL(KLOC)**2+YL(KLOC)**2+ZL(KLOC)**2)
       XL(KLOC)=XL(KLOC)*RL/RAD
       YL(KLOC)=YL(KLOC)*RL/RAD
       ZL(KLOC)=ZL(KLOC)*RL/RAD
       iloc=2
       IGL=XONDGL((ELE-1)*NLOC+ILOC)
       ! The corner node iloc...
       kloc=3
       XL(KLOC)=X(IGL)
       YL(KLOC)=Y(IGL)
       ZL(KLOC)=Z(IGL)
       jloc=3
       JGL=XONDGL((ELE-1)*NLOC+JLOC)
       kloc=5
       XL(KLOC)=0.5*(X(IGL)+X(JGL))
       YL(KLOC)=0.5*(Y(IGL)+Y(JGL))
       ZL(KLOC)=0.5*(Z(IGL)+Z(JGL))
       ! Adjust the coords so we have a linear variation in radius.
       RL=0.5*(RADLOC(ILOC)+RADLOC(JLOC))
       RAD=SQRT(XL(KLOC)**2+YL(KLOC)**2+ZL(KLOC)**2)
       XL(KLOC)=XL(KLOC)*RL/RAD
       YL(KLOC)=YL(KLOC)*RL/RAD
       ZL(KLOC)=ZL(KLOC)*RL/RAD
       jloc=4
       JGL=XONDGL((ELE-1)*NLOC+JLOC)
       kloc=8
       XL(KLOC)=0.5*(X(IGL)+X(JGL))
       YL(KLOC)=0.5*(Y(IGL)+Y(JGL))
       ZL(KLOC)=0.5*(Z(IGL)+Z(JGL))
       ! Adjust the coords so we have a linear variation in radius.
       RL=0.5*(RADLOC(ILOC)+RADLOC(JLOC))
       RAD=SQRT(XL(KLOC)**2+YL(KLOC)**2+ZL(KLOC)**2)
       XL(KLOC)=XL(KLOC)*RL/RAD
       YL(KLOC)=YL(KLOC)*RL/RAD
       ZL(KLOC)=ZL(KLOC)*RL/RAD
       iloc=3
       IGL=XONDGL((ELE-1)*NLOC+ILOC)
       ! The corner node iloc...
       kloc=6
       XL(KLOC)=X(IGL)
       YL(KLOC)=Y(IGL)
       ZL(KLOC)=Z(IGL)
       jloc=4
       JGL=XONDGL((ELE-1)*NLOC+JLOC)
       kloc=9
       XL(KLOC)=0.5*(X(IGL)+X(JGL))
       YL(KLOC)=0.5*(Y(IGL)+Y(JGL))
       ZL(KLOC)=0.5*(Z(IGL)+Z(JGL))
       ! Adjust the coords so we have a linear variation in radius.
       RL=0.5*(RADLOC(ILOC)+RADLOC(JLOC))
       RAD=SQRT(XL(KLOC)**2+YL(KLOC)**2+ZL(KLOC)**2)
       XL(KLOC)=XL(KLOC)*RL/RAD
       YL(KLOC)=YL(KLOC)*RL/RAD
       ZL(KLOC)=ZL(KLOC)*RL/RAD
       iloc=4
       IGL=XONDGL((ELE-1)*NLOC+ILOC)
       ! The corner node iloc...
       kloc=10
       XL(KLOC)=X(IGL)
       YL(KLOC)=Y(IGL)
       ZL(KLOC)=Z(IGL)

    endif
    RETURN
  END SUBROUTINE UNWIND_XL_CALC

  SUBROUTINE DETNLXMAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MLOC,NGI, &
       NLX,NLY,NLZ,MLX,MLY,MLZ, WEIGHT, DETWEI, &
       NX,NY,NZ, MX,MY,MZ,&
       A11,A12,A13, A21,A22,A23, A31,A32,A33)
    ! Calculate DETWEI,NX,NY,NZ, MX,MY,MZ for element ELE
    ! For coefficient in the inverse mat of the jacobian.
    INTEGER ELE,TOTELE,XNONOD,NLOC,MLOC,NGI
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    INTEGER XONDGL(TOTELE*NLOC)
    ! Volume shape functions...
    REAL NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
    REAL WEIGHT(NGI)
    REAL DETWEI(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    REAL MX(MLOC,NGI),MY(MLOC,NGI),MZ(MLOC,NGI)
    REAL A11(NGI),A12(NGI),A13(NGI)
    REAL A21(NGI),A22(NGI),A23(NGI)
    REAL A31(NGI),A32(NGI),A33(NGI)
    ! Local variables...
    REAL AGI,BGI,CGI, DGI,EGI,FGI, GGI,HGI,KGI
    REAL DETJ,VOLUME
    INTEGER GI,L,IGLX

    VOLUME=0.
    DO GI=1,NGI

       AGI=0.
       BGI=0.
       CGI=0.

       DGI=0.
       EGI=0.
       FGI=0.

       GGI=0.
       HGI=0.
       KGI=0.

       DO L=1,NLOC
          IGLX=XONDGL((ELE-1)*NLOC+L)
          ! NB R0 does not appear here although the z-coord might be Z+R0. 
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

       DETJ=AGI*(EGI*KGI-FGI*HGI) &
            -BGI*(DGI*KGI-FGI*GGI) &
            +CGI*(DGI*HGI-EGI*GGI)
       DETWEI(GI)=abs(DETJ)*WEIGHT(GI)
       VOLUME=VOLUME+DETWEI(GI)

       A11(GI)= (EGI*KGI-FGI*HGI) /DETJ
       A21(GI)=-(DGI*KGI-FGI*GGI) /DETJ
       A31(GI)= (DGI*HGI-EGI*GGI) /DETJ

       A12(GI)=-(BGI*KGI-CGI*HGI) /DETJ
       A22(GI)= (AGI*KGI-CGI*GGI) /DETJ
       A32(GI)=-(AGI*HGI-BGI*GGI) /DETJ

       A13(GI)= (BGI*FGI-CGI*EGI) /DETJ
       A23(GI)=-(AGI*FGI-CGI*DGI) /DETJ
       A33(GI)= (AGI*EGI-BGI*DGI) /DETJ

       DO L=1,NLOC
          NX(L,GI)= A11(GI)*NLX(L,GI)+A12(GI)*NLY(L,GI)+A13(GI)*NLZ(L,GI)
          NY(L,GI)= A21(GI)*NLX(L,GI)+A22(GI)*NLY(L,GI)+A23(GI)*NLZ(L,GI)
          NZ(L,GI)= A31(GI)*NLX(L,GI)+A32(GI)*NLY(L,GI)+A33(GI)*NLZ(L,GI)
       end do

       DO L=1,MLOC
          MX(L,GI)= A11(GI)*MLX(L,GI)+A12(GI)*MLY(L,GI)+A13(GI)*MLZ(L,GI)
          MY(L,GI)= A21(GI)*MLX(L,GI)+A22(GI)*MLY(L,GI)+A23(GI)*MLZ(L,GI)
          MZ(L,GI)= A31(GI)*MLX(L,GI)+A32(GI)*MLY(L,GI)+A33(GI)*MLZ(L,GI)
       ENDDO

    end do
  END SUBROUTINE DETNLXMAI

  subroutine detnlxmsup(nloc,mloc,ngi, &
       nlx,nly,nlz, mlx,mly,mlz, weight, detwei, &
       nx,ny,nz, mx,my,mz, &
       ncloc,nclx,ncly,nclz,xl,yl,zl, &
       a11,a12,a13, a21,a22,a23, a31,a32,a33)
    ! calculate detwei,nx,ny,nz, mx,my,mz for element ele
    ! for coefficient in the inverse mat of the jacobian. 
    ! this works for a superparametric fem nc,nclx,ncly,nclz are basis functions 
    ! associated with coordinates.
    ! (xl,yl,zl) are th local coordinates of the nodes.
    integer, intent(in)::nloc,mloc,ngi,ncloc
    real, intent(in):: xl(ncloc),yl(ncloc),zl(ncloc)
    real, intent(out):: a11(ngi),a12(ngi),a13(ngi)
    real, intent(out):: a21(ngi),a22(ngi),a23(ngi)
    real, intent(out):: a31(ngi),a32(ngi),a33(ngi) 
    ! volume shape functions...
    real, intent(in)::nlx(nloc,ngi),nly(nloc,ngi),nlz(nloc,ngi)
    real, intent(in)::mlx(mloc,ngi),mly(mloc,ngi),mlz(mloc,ngi)
    real, intent(in)::nclx(ncloc,ngi),ncly(ncloc,ngi),nclz(ncloc,ngi)
    real, intent(in)::weight(ngi)
    real, intent(out)::detwei(ngi)
    real, intent(out)::nx(nloc,ngi),ny(nloc,ngi),nz(nloc,ngi)
    real, intent(out)::mx(mloc,ngi),my(mloc,ngi),mz(mloc,ngi)
    ! local variables...
    real agi,bgi,cgi, dgi,egi,fgi, ggi,hgi,kgi
    real detj,volume
    integer gi,l

    volume=0.
    do gi=1,ngi

       agi=0.
       bgi=0.
       cgi=0.

       dgi=0.
       egi=0.
       fgi=0.

       ggi=0.
       hgi=0.
       kgi=0.

       do l=1,ncloc
          agi=agi+nclx(l,gi)*xl(l) 
          bgi=bgi+nclx(l,gi)*yl(l) 
          cgi=cgi+nclx(l,gi)*zl(l) 

          dgi=dgi+ncly(l,gi)*xl(l) 
          egi=egi+ncly(l,gi)*yl(l) 
          fgi=fgi+ncly(l,gi)*zl(l) 

          ggi=ggi+nclz(l,gi)*xl(l) 
          hgi=hgi+nclz(l,gi)*yl(l) 
          kgi=kgi+nclz(l,gi)*zl(l) 
       end do

       detj=agi*(egi*kgi-fgi*hgi) &
            -bgi*(dgi*kgi-fgi*ggi) &
            +cgi*(dgi*hgi-egi*ggi)
       detwei(gi)=abs(detj)*weight(gi)
       volume=volume+detwei(gi)

       a11(gi)= (egi*kgi-fgi*hgi) /detj
       a21(gi)=-(dgi*kgi-fgi*ggi) /detj
       a31(gi)= (dgi*hgi-egi*ggi) /detj

       a12(gi)=-(bgi*kgi-cgi*hgi) /detj
       a22(gi)= (agi*kgi-cgi*ggi) /detj
       a32(gi)=-(agi*hgi-bgi*ggi) /detj

       a13(gi)= (bgi*fgi-cgi*egi) /detj
       a23(gi)=-(agi*fgi-cgi*dgi) /detj
       a33(gi)= (agi*egi-bgi*dgi) /detj

       do l=1,nloc
          nx(l,gi)= a11(gi)*nlx(l,gi)+a12(gi)*nly(l,gi)+a13(gi)*nlz(l,gi)
          ny(l,gi)= a21(gi)*nlx(l,gi)+a22(gi)*nly(l,gi)+a23(gi)*nlz(l,gi)
          nz(l,gi)= a31(gi)*nlx(l,gi)+a32(gi)*nly(l,gi)+a33(gi)*nlz(l,gi)
       end do

       do l=1,mloc
          mx(l,gi)= a11(gi)*mlx(l,gi)+a12(gi)*mly(l,gi)+a13(gi)*mlz(l,gi)
          my(l,gi)= a21(gi)*mlx(l,gi)+a22(gi)*mly(l,gi)+a23(gi)*mlz(l,gi)
          mz(l,gi)= a31(gi)*mlx(l,gi)+a32(gi)*mly(l,gi)+a33(gi)*mlz(l,gi)
       enddo

    end do

  end subroutine detnlxmsup

  SUBROUTINE DGSDETNXLOC(ELE,SILOC2ILOC,XONDGL, &
       TOTELE,NLOC,XNONOD,SNLOC,SNGI, &
       X,Y,Z, &
       SN,SNLX,SNLY, SWEIGH, SDETWE,SAREA, D3,DCYL, &
       NORMXN,NORMYN,NORMZN, &
       NORMX,NORMY,NORMZ, &
       ISPHERE,SMLOC,SM,SMLX,SMLY)
    INTEGER TOTELE,NLOC,XNONOD,SNLOC,SNGI
    INTEGER ELE,SILOC2ILOC(SNLOC),XONDGL(TOTELE*NLOC)
    INTEGER SMLOC,ISPHERE
    ! ISPHERE if ge 1 represents the polynomial of the superparametric
    ! mapping of the spherical earth.
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
    REAL SM(SNLOC,SNGI),SMLX(SNLOC,SNGI),SMLY(SNLOC,SNGI)
    REAL SWEIGH(SNGI)
    REAL SDETWE(SNGI)
    REAL SAREA
    LOGICAL D3,DCYL
    REAL NORMXN(SNGI),NORMYN(SNGI),NORMZN(SNGI)
    REAL NORMX,NORMY,NORMZ
    ! Local variables...
    REAL XSL(SMLOC),YSL(SMLOC),ZSL(SMLOC)
    INTEGER SL,IGLX,ILOC,JLOC,SILOC,SJLOC,INOD,JNOD,KLOC,SKLOC
    REAL RADI,RADJ,RADP,RAD

    IF(ISPHERE.LE.1) THEN
       DO SL=1,SNLOC
          ILOC=SILOC2ILOC(SL)
          IGLX=XONDGL((ELE-1)*NLOC+ILOC)
          XSL(SL)=X(IGLX)
          YSL(SL)=Y(IGLX)
          IF(D3) ZSL(SL)=Z(IGLX)
       END DO

       CALL DGSDETNXLOC2(SNLOC,SNGI, &
            XSL,YSL,ZSL, &
            SN,SNLX,SNLY, SWEIGH, SDETWE,SAREA, D3,DCYL, &
            NORMXN,NORMYN,NORMZN, &
            NORMX,NORMY,NORMZ)
    ELSE
       DO SILOC=1,SNLOC
          ILOC=SILOC2ILOC(SILOC)
          INOD=XONDGL((ELE-1)*NLOC+ILOC)
          DO SJLOC=SILOC,SNLOC
             JLOC=SILOC2ILOC(SJLOC)
             JNOD=XONDGL((ELE-1)*NLOC+JLOC)
             KLOC=ILINK2(ILOC,JLOC)
             SKLOC=SILINK2(SILOC,SJLOC)
             XSL(SKLOC)=0.5*(X(INOD)+X(JNOD))
             YSL(SKLOC)=0.5*(Y(INOD)+Y(JNOD))
             ZSL(SKLOC)=0.5*(Z(INOD)+Z(JNOD))
             ! Assume the mid side nodes are on the sphere
             IF(ILOC.NE.JLOC) THEN
                RADI=SQRT(X(INOD)**2+Y(INOD)**2+Z(INOD)**2)
                RADJ=SQRT(X(JNOD)**2+Y(JNOD)**2+Z(JNOD)**2)
                RADP=0.5*(RADI+RADJ)
                RAD=SQRT(XSL(SKLOC)**2+YSL(SKLOC)**2+ZSL(SKLOC)**2)
                XSL(SKLOC)=XSL(SKLOC)*RADP/RAD
                YSL(SKLOC)=YSL(SKLOC)*RADP/RAD
                ZSL(SKLOC)=ZSL(SKLOC)*RADP/RAD
             ENDIF
          END DO
       END DO

       CALL DGSDETNXLOC2(SMLOC,SNGI, &
            XSL,YSL,ZSL, &
            SM,SMLX,SMLY, SWEIGH, SDETWE,SAREA, D3,DCYL, &
            NORMXN,NORMYN,NORMZN, &
            NORMX,NORMY,NORMZ)

    ENDIF

  end SUBROUTINE DGSDETNXLOC

  SUBROUTINE DGSDETNXLOC2(SNLOC,SNGI, &
       XSL,YSL,ZSL, &
       SN,SNLX,SNLY, SWEIGH, SDETWE,SAREA, D3,DCYL, &
       NORMXN,NORMYN,NORMZN, &
       NORMX,NORMY,NORMZ)
    INTEGER SNLOC,SNGI
    REAL XSL(SNLOC),YSL(SNLOC),ZSL(SNLOC)
    REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
    REAL SWEIGH(SNGI)
    REAL SDETWE(SNGI)
    REAL SAREA
    LOGICAL D3,DCYL
    REAL NORMXN(SNGI),NORMYN(SNGI),NORMZN(SNGI)
    REAL NORMX,NORMY,NORMZ
    real, parameter::pi=3.141592654
    INTEGER GI,SL,IGLX
    REAL DXDLX,DXDLY,DYDLX,DYDLY,DZDLX,DZDLY
    REAL A,B,C,DETJ,RGI,TWOPI

    SAREA=0.

    IF(D3) THEN
       DO GI=1,SNGI

          DXDLX=0.
          DXDLY=0.
          DYDLX=0.
          DYDLY=0.
          DZDLX=0.
          DZDLY=0.

          DO SL=1,SNLOC
             DXDLX=DXDLX + SNLX(SL,GI)*XSL(SL)
             DXDLY=DXDLY + SNLY(SL,GI)*XSL(SL)
             DYDLX=DYDLX + SNLX(SL,GI)*YSL(SL)
             DYDLY=DYDLY + SNLY(SL,GI)*YSL(SL)
             DZDLX=DZDLX + SNLX(SL,GI)*ZSL(SL)
             DZDLY=DZDLY + SNLY(SL,GI)*ZSL(SL)
          END DO
          A = DYDLX*DZDLY - DYDLY*DZDLX
          B = DXDLX*DZDLY - DXDLY*DZDLX
          C = DXDLX*DYDLY - DXDLY*DYDLX

          DETJ=SQRT( A**2 + B**2 + C**2)
          SDETWE(GI)=DETJ*SWEIGH(GI)
          SAREA=SAREA+SDETWE(GI)

          ! Calculate the normal at the Gauss pts...
          ! Perform x-product. N=T1 x T2
          CALL NORMGI(NORMXN(GI),NORMYN(GI),NORMZN(GI), &
               DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,DZDLY, &
               NORMX,NORMY,NORMZ)

       END DO
       ! IF(D3) THEN...
    ELSE
       ! ENDOF IF(D3) THEN ELSE...
       TWOPI=1.0
       IF(DCYL) TWOPI=2.*PI

       DO GI=1,SNGI
          RGI=0.
          DXDLX=0.
          DXDLY=0.
          DYDLX=0.
          DYDLY=0.
          DZDLX=0.
          ! DZDLY=1 is to calculate the normal.
          DZDLY=1.
          DO SL=1,SNLOC
             DXDLX=DXDLX + SNLX(SL,GI)*XSL(SL)
             DYDLX=DYDLX + SNLX(SL,GI)*YSL(SL)
             RGI=RGI+SN(SL,GI)*YSL(IGLX)
          END DO
          IF(.NOT.DCYL) RGI=1.0
          DETJ=SQRT( DXDLX**2 + DYDLX**2 )
          SDETWE(GI)=TWOPI*RGI*DETJ*SWEIGH(GI)
          SAREA=SAREA+SDETWE(GI)
          NORMZ=0.
          CALL NORMGI(NORMXN(GI),NORMYN(GI),NORMZN(GI), &
               DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,DZDLY, &
               NORMX,NORMY,NORMZ)
       END DO
    ENDIF

  end SUBROUTINE DGSDETNXLOC2

  SUBROUTINE SCVDETNX( ELE,      GI,      NODI,  &
                                !     - INTEGERS
       NLOC,     SVNGI,   TOTELE,   &
       XNDGLN,   XNONOD,&
                                !     - REALS
       CVDETWEI, CVNORMX, CVNORMY,  &
       CVNORMZ,  SVN,     SVNLX,    &
       SVNLY,    SVWEIGH, PNTX,     &
       PNTY,     PNTZ,    X,        &
       Y,        Z,  &
                                !     - LOGICALS
       D3,       DCYL,    GOTPNTS )
    !     --------------------------------------------------
    !     
    !     - this subroutine calculates the control volume (CV) 
    !     - CVNORMX, CVNORMY, CVNORMZ normals at the Gaussian 
    !     - integration points GI. NODI = the current global 
    !     - node number for the co-ordinates. 
    !     
    !     -------------------------------
    !     - date last modified : 15/03/2003
    !     -------------------------------
    INTEGER ELE,    GI
    INTEGER NODI,   NLOC
    INTEGER SVNGI,  TOTELE
    INTEGER XNONOD     

    INTEGER XNDGLN(TOTELE*NLOC)

    REAL CVDETWEI(SVNGI)   
    REAL CVNORMX(SVNGI),    CVNORMY(SVNGI)
    REAL CVNORMZ(SVNGI),    SVN(NLOC,SVNGI)
    REAL SVNLX(NLOC,SVNGI), SVNLY(NLOC,SVNGI)
    REAL SVWEIGH(SVNGI),    X(XNONOD)
    REAL Y(XNONOD),         Z(XNONOD)

    LOGICAL D3,DCYL, GOTPNTS

    !     - Local variables
    INTEGER NODJ,JLOC

    REAL A, B, C
    REAL DETJ
    REAL DXDLX, DXDLY, DYDLX
    REAL DYDLY, DZDLX, DZDLY

    REAL TWOPI
    real, parameter::pi=3.14159265

    REAL POSVGIX, POSVGIY, POSVGIZ
    REAL PNTX,    PNTY,    PNTZ
    REAL RGI

    IF( D3 ) THEN

       DXDLX = 0.0
       DXDLY = 0.0

       DYDLX = 0.0
       DYDLY = 0.0

       DZDLX = 0.0
       DZDLY = 0.0

       POSVGIX = 0.0
       POSVGIY = 0.0
       POSVGIZ = 0.0

       do  JLOC = 1, NLOC! Was loop 100

          NODJ = XNDGLN((ELE-1)*NLOC+JLOC)

          DXDLX = DXDLX + SVNLX(JLOC,GI)*X(NODJ)
          DXDLY = DXDLY + SVNLY(JLOC,GI)*X(NODJ) 
          DYDLX = DYDLX + SVNLX(JLOC,GI)*Y(NODJ) 
          DYDLY = DYDLY + SVNLY(JLOC,GI)*Y(NODJ) 
          DZDLX = DZDLX + SVNLX(JLOC,GI)*Z(NODJ) 
          DZDLY = DZDLY + SVNLY(JLOC,GI)*Z(NODJ) 

          POSVGIX = POSVGIX + SVN(JLOC,GI)*X(NODJ)
          POSVGIY = POSVGIY + SVN(JLOC,GI)*Y(NODJ)
          POSVGIZ = POSVGIZ + SVN(JLOC,GI)*Z(NODJ)
       end do ! Was loop 100

       !     - Note that POSVGIX,POSVGIY and POSVGIZ can be considered as the 
       !     - components of the Gauss pnt GI with the co-ordinate origin 
       !     - positioned at the current control volume NODI.

       IF( GOTPNTS ) THEN

          POSVGIX = POSVGIX - PNTX 
          POSVGIY = POSVGIY - PNTY 
          POSVGIZ = POSVGIZ - PNTZ 

       ELSE
          PNTX = POSVGIX
          PNTY = POSVGIY
          PNTZ = POSVGIZ

          POSVGIX = POSVGIX - X(NODI)
          POSVGIY = POSVGIY - Y(NODI)
          POSVGIZ = POSVGIZ - Z(NODI)

       END IF

       A = DYDLX*DZDLY - DYDLY*DZDLX
       B = DXDLX*DZDLY - DXDLY*DZDLX
       C = DXDLX*DYDLY - DXDLY*DYDLX
       !
       !     - Calculate the determinant of the Jacobian at Gauss pnt GI.
       !     
       DETJ = SQRT( A**2 + B**2 + C**2 )
       !     
       !     - Calculate the determinant times the surface weight at Gauss pnt GI.
       !     
       CVDETWEI(GI) = DETJ*SVWEIGH(GI)
       !     
       !     - Calculate the normal at the Gauss pts
       !     - TANX1 = DXDLX, TANY1 = DYDLX, TANZ1 = DZDLX,    
       !     - TANX2 = DXDLY, TANY2 = DYDLY, TANZ2 = DZDLY
       !     - Perform cross-product. N = T1 x T2
       !     
       CALL NORMGI( CVNORMX(GI), CVNORMY(GI), CVNORMZ(GI),&
            DXDLX,       DYDLX,       DZDLX, &
            DXDLY,       DYDLY,       DZDLY,&
            POSVGIX,     POSVGIY,     POSVGIZ ) 

    ELSE 

       TWOPI = 1.0

       IF( DCYL ) TWOPI = 2.0*PI

       RGI   = 0.0

       DXDLX = 0.0
       DXDLY = 0.0

       DYDLX = 0.0
       DYDLY = 0.0

       DZDLX = 0.0
       !     
       !     - Note that we set the derivative wrt to y of coordinate z to 1.0
       !     
       DZDLY = 1.0

       POSVGIX = 0.0
       POSVGIY = 0.0
       POSVGIZ = 0.0

       do  JLOC = 1, NLOC! Was loop 300

          NODJ = XNDGLN((ELE-1)*NLOC+JLOC)

          DXDLX = DXDLX + SVNLX(JLOC,GI)*X(NODJ) 
          DYDLX = DYDLX + SVNLX(JLOC,GI)*Y(NODJ) 

          POSVGIX = POSVGIX + SVN(JLOC,GI)*X(NODJ)
          POSVGIY = POSVGIY + SVN(JLOC,GI)*Y(NODJ)

          RGI = RGI + SVN(JLOC,GI)*Y(NODJ)

       end do ! Was loop 300
       !     
       !     - Note that POSVGIX and POSVGIY can be considered as the components 
       !     - of the Gauss pnt GI with the co-ordinate origin positioned at the
       !     - current control volume NODI.
       !     
       IF( GOTPNTS ) THEN

          POSVGIX = POSVGIX - PNTX 
          POSVGIY = POSVGIY - PNTY 

       ELSE

          POSVGIX = POSVGIX - X(NODI)
          POSVGIY = POSVGIY - Y(NODI)

       END IF

       IF( .NOT. DCYL ) RGI = 1.0 

       DETJ = SQRT( DXDLX**2 + DYDLX**2 )
       CVDETWEI(GI)  = TWOPI*RGI*DETJ*SVWEIGH(GI)
       !
       !     - Calculate the normal at the Gauss pts
       !     - TANX1 = DXDLX, TANY1 = DYDLX, TANZ1 = DZDLX,    
       !     - TANX2 = DXDLY, TANY2 = DYDLY, TANZ2 = DZDLY
       !     - Perform cross-product. N = T1 x T2
       !     
       CALL NORMGI( CVNORMX(GI), CVNORMY(GI), CVNORMZ(GI),&
            DXDLX,       DYDLX,       DZDLX, &
            DXDLY,       DYDLY,       DZDLY,&
            POSVGIX,     POSVGIY,     POSVGIZ )
       !     
       !     - End of GI loop
       !     
    ENDIF

  END SUBROUTINE SCVDETNX

  SUBROUTINE DGSIMPLNORM(ELE,SILOC2ILOC,TOTELE,NLOC,SNLOC,XONDGL, &
       X,Y,Z,XNONOD,NORMX,NORMY,NORMZ)
    ! Form approximate surface normal (NORMX,NORMY,NORMZ)
    INTEGER ELE,TOTELE,NLOC,SNLOC,XNONOD
    INTEGER SILOC2ILOC(SNLOC),XONDGL(TOTELE*NLOC)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    REAL NORMX,NORMY,NORMZ
    ! Local variables...
    REAL XC,YC,ZC,SXC,SYC,SZC,NORM,NDOTN
    REAL APRNORMX,APRNORMY,APRNORMZ
    INTEGER XNODI,ILOC,SILOC,IGL1,IGL2,IGL3
    XC=0.0
    YC=0.0
    ZC=0.0
    DO ILOC=1,NLOC
       XNODI=XONDGL((ELE-1)*NLOC+ILOC)
       XC=XC+X(XNODI)/REAL(NLOC)
       YC=YC+Y(XNODI)/REAL(NLOC)
       ZC=ZC+Z(XNODI)/REAL(NLOC)
    END DO

    SXC=0.0
    SYC=0.0
    SZC=0.0

    DO SILOC=1,SNLOC
       ILOC=SILOC2ILOC(SILOC)
       XNODI=XONDGL((ELE-1)*NLOC+ILOC)
       SXC=SXC+X(XNODI)/REAL(SNLOC)
       SYC=SYC+Y(XNODI)/REAL(SNLOC)
       SZC=SZC+Z(XNODI)/REAL(SNLOC)
    END DO
    APRNORMX=SXC-XC
    APRNORMY=SYC-YC
    APRNORMZ=SZC-ZC

    IGL1 =XONDGL((ELE-1)*NLOC+SILOC2ILOC(1))
    IGL2 =XONDGL((ELE-1)*NLOC+SILOC2ILOC(2))
    IGL3 =XONDGL((ELE-1)*NLOC+SILOC2ILOC(3))

    CALL XPROD(NORMX,NORMY,NORMZ, &
         X(IGL1)-X(IGL2),Y(IGL1)-Y(IGL2),Z(IGL1)-Z(IGL2), &
         X(IGL1)-X(IGL3),Y(IGL1)-Y(IGL3),Z(IGL1)-Z(IGL3) )

    NDOTN=APRNORMX*NORMX+APRNORMY*NORMY+APRNORMZ*NORMZ
    IF(NDOTN.LT.0.0) THEN
       NORMX=-NORMX
       NORMY=-NORMY
       NORMZ=-NORMZ
    ENDIF
    NORM=SQRT(NORMX**2+NORMY**2+NORMZ**2)
    NORMX=NORMX/NORM
    NORMY=NORMY/NORM
    NORMZ=NORMZ/NORM

  END SUBROUTINE DGSIMPLNORM

  SUBROUTINE NORMGI(NORMXN,NORMYN,NORMZN, &
       DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,DZDLY, &
       NORMX,NORMY,NORMZ) 
    !     Calculate the normal at the Gauss pts...
    !     Perform x-product. N=T1 x T2
    REAL NORMXN,NORMYN,NORMZN
    real     DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,DZDLY
    real     NORMX,NORMY,NORMZ
    !     Local variables...
    REAL RN,SIRN

    CALL XPROD(NORMXN,NORMYN,NORMZN, &
         DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,DZDLY)

    RN=SQRT(NORMXN**2+NORMYN**2+NORMZN**2)
    SIRN=SIGN(1.0/RN,NORMXN*NORMX+NORMYN*NORMY+NORMZN*NORMZ)
    NORMXN=SIRN*NORMXN
    NORMYN=SIRN*NORMYN
    NORMZN=SIRN*NORMZN

  END subroutine normgi
    
    subroutine transform_to_physical_fast_3d(X, dn, invJ, detJ)
    !! Column n of X is the position of the nth node. (3 x nloc)
    !! only need position of n nodes since Jacobian is only calculated once
    real, dimension(:,:), intent(in) :: X
    !! Derivative of shape function (nloc x ngi x dim)
    real, dimension(:,:,:), intent(in) :: dn
    !! Inverse Jacobian at each gauss point (3 x 3 x ngi)
    real, dimension(:,:,:), intent(out) :: invJ
    !! Determinant of Jacobian at each gauss point (ngi)
    real, dimension(:), intent(out) :: detJ
    
    integer, parameter:: cyc3(1:5)=(/ 1, 2, 3, 1, 2 /)
    real, dimension(3,3) :: J_T ! transpose of Jacobian
    integer gi, i, k
    
    do gi=1, size(dn,2)
      
      J_t=matmul(X(:,:), dn(:, gi, :))
      ! Calculate (scaled) inverse using recursive determinants.
      forall (i=1:3,k=1:3) 
        invJ(i,k,gi)=J_t(cyc3(i+1),cyc3(k+1))*J_t(cyc3(i+2),cyc3(k+2)) &
            -J_t(cyc3(i+2),cyc3(k+1))*J_t(cyc3(i+1),cyc3(k+2))
      end forall
      ! Form determinant by expanding minors.
      detJ(gi)=dot_product(J_t(:,1), invJ(:,1,gi))
      
      ! Scale inverse by determinant.
      invJ(:,:,gi)=invJ(:,:,gi)/detJ(gi)
      
    end do
  
  end subroutine transform_to_physical_fast_3d

end module Shape_Transformations
