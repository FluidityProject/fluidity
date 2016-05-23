!    Copyright (C) 2008 Imperial College London and others.
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


module alturbine_utils

    implicit none

    public QuatRot, cross

contains 
 
    !SUBROUTINE CalcBEGeom(BNum)

    !    implicit none

    !    integer :: BNum
    !    integer :: nbe, nei, FlipN, nej, j
    !    real :: sEM, tEM, nEM
    !    real :: PE(3), sE(3), tE(3), normE(3), P1(3), P2(3), P3(3), P4(3), V1(3), V2(3), V3(3), V4(3), A1(3), A2(3)

    !    ! Calculates element geometry from element end geometry

    !    ! JCM: Eventually, should just be able to loop through Blades(BNum) data structure
    !    ! While data is still held in arrays concatenated across blades, need to replicate
    !    ! nbe (stored in configr) from Blades(1).NElem
    !    nbe=Blades(1)%NElem

    !    FlipN=Blades(BNum)%FlipN

    !    nei=1+(BNum-1)*(nbe+1)

    !    do j=1,nbe
    !        nej=nei+j

    !        ! Element center locations
    !        xBC(nej)=(xBE(nej)+xBE(nej-1))/2.0
    !        yBC(nej)=(yBE(nej)+yBE(nej-1))/2.0
    !        zBC(nej)=(zBE(nej)+zBE(nej-1))/2.0

    !        ! Set spanwise and tangential vectors
	!	    sE=-(/xBE(nej)-xBE(nej-1),yBE(nej)-yBE(nej-1),zBE(nej)-zBE(nej-1)/) ! nominal element spanwise direction set opposite to QC line
	!	    sEM=sqrt(dot_product(sE,sE))
	!	    sE=sE/sEM
	!	    tE=(/txBE(nej)+txBE(nej-1),tyBE(nej)+tyBE(nej-1),tzBE(nej)+tzBE(nej-1)/)/2.0
	!	    ! Force tE normal to sE
	!	    tE=tE-dot_product(tE,sE)*sE
	!	    tEM=sqrt(dot_product(tE,tE))
	!	    tE=tE/tEM
	!	    sxBC(nej)=sE(1)
	!	    syBC(nej)=sE(2)
	!	    szBC(nej)=sE(3)
    !        txBC(nej)=tE(1)
    !        tyBC(nej)=tE(2)
    !        tzBC(nej)=tE(3)

	!	    ! Calc normal vector
	!	    Call cross(sE(1),sE(2),sE(3),tE(1),tE(2),tE(3),normE(1),normE(2),normE(3))
	!	    nEM=sqrt(dot_product(normE,normE))
	!	    normE=normE/nEM
	!	    nxBC(nej)=normE(1)
	!	    nyBC(nej)=normE(2)
	!	    nzBC(nej)=normE(3)

	!	    ! Flip normal direction if requested
	!	    CircSign(nej)=1.0
	!	    if (FlipN .eq. 1) then
	!            nxBC(nej)=-nxBC(nej)
	!            nyBC(nej)=-nyBC(nej)
	!            nzBC(nej)=-nzBC(nej)
    !            sxBC(nej)=-sxBC(nej)
    !            syBC(nej)=-syBC(nej)
    !            szBC(nej)=-szBC(nej)
    !            CircSign(nej)=-1.0
	!	    end if

	!	    ! Calc element area and chord
	!	    P1=(/xBE(nej-1)-0.25*CtoR(nej-1)*txBE(nej-1),yBE(nej-1)-0.25*CtoR(nej-1)*tyBE(nej-1),zBE(nej-1)-0.25*CtoR(nej-1)*tzBE(nej-1)/)
	!	    P2=(/xBE(nej-1)+0.75*CtoR(nej-1)*txBE(nej-1),yBE(nej-1)+0.75*CtoR(nej-1)*tyBE(nej-1),zBE(nej-1)+0.75*CtoR(nej-1)*tzBE(nej-1)/)
	!	    P3=(/xBE(nej)+0.75*CtoR(nej)*txBE(nej),yBE(nej)+0.75*CtoR(nej)*tyBE(nej),zBE(nej)+0.75*CtoR(nej)*tzBE(nej)/)
	!	    P4=(/xBE(nej)-0.25*CtoR(nej)*txBE(nej),yBE(nej)-0.25*CtoR(nej)*tyBE(nej),zBE(nej)-0.25*CtoR(nej)*tzBE(nej)/)
	!	    V1=P2-P1
	!	    V2=P3-P2
	!	    V3=P4-P3
	!	    V4=P1-P4
	!	    ! Calc quad area from two triangular facets
	!	    Call cross(V1(1),V1(2),V1(3),V2(1),V2(2),V2(3),A1(1),A1(2),A1(3))
	!	    A1=A1/2.0
    !        Call cross(V3(1),V3(2),V3(3),V4(1),V4(2),V4(3),A2(1),A2(2),A2(3))
    !        A2=A2/2.0
	!	    eArea(nej)=sqrt(dot_product(A1,A1))+sqrt(dot_product(A2,A2))
	!	    ! Calc average element chord from area and span
	!	    eChord(nej)=eArea(nej)/sEM

    !    end do
    !End SUBROUTINE CalcBEGeom
    

    SUBROUTINE cross(ax,ay,az,bx,by,bz,cx,cy,cz) 

        real ax,ay,az,bx,by,bz,cx,cy,cz 

        cx = ay*bz - az*by
        cy = az*bx - ax*bz
        cz = ax*by - ay*bx

    End SUBROUTINE cross


    subroutine QuatRot(vx,vy,vz,Theta,nRx,nRy,nRz,Ox,Oy,Oz,vRx,vRy,vRz)

        ! % Perform rotation of vector v around normal vector nR using the
        ! % quaternion machinery.
        ! % v: input vector
        ! % Theta: rotation angle (rad)
        ! % nR: normal vector around which to rotate
        ! % Origin: origin point of rotation
        ! %
        ! % vR: Rotated vector

        real :: vx,vy,vz,Theta,nRx,nRy,nRz,Ox,Oy,Oz,vRx,vRy,vRz     

        real :: p(4,1), pR(4,1), q(4), qbar(4), nRMag, vOx, vOy, vOz
        real :: QL(4,4), QbarR(4,4)

        ! Force normalize nR
        nRMag=sqrt(nRx**2+nRy**2+nRz**2)
        nRx=nRx/nRMag
        nRy=nRy/nRMag
        nRz=nRz/nRMag

        ! Quaternion form of v
        vOx=vx-Ox
        vOy=vy-Oy
        vOz=vz-Oz
        p=reshape((/0.0,vOx,vOy,vOz/),(/4,1/))

        ! Rotation quaternion and conjugate
        q=(/cos(Theta/2),nRx*sin(Theta/2),nRy*sin(Theta/2),nRz*sin(Theta/2)/)
        qbar=(/q(1),-q(2),-q(3),-q(4)/)

        QL=transpose(reshape((/q(1), -q(2), -q(3), -q(4), &
            q(2),  q(1), -q(4),  q(3), &
            q(3),  q(4),  q(1), -q(2), &
            q(4), -q(3),  q(2),  q(1)/),(/4,4/)))

        QbarR=transpose(reshape((/qbar(1), -qbar(2), -qbar(3), -qbar(4), &
            qbar(2),  qbar(1),  qbar(4), -qbar(3), &
            qbar(3), -qbar(4),  qbar(1),  qbar(2), &
            qbar(4),  qbar(3), -qbar(2),  qbar(1)/),(/4,4/)))

        ! Rotate p
        pR=matmul(matmul(QbarR,QL),p)
        vRx=pR(2,1)+Ox
        vRy=pR(3,1)+Oy
        vRz=pR(4,1)+Oz

    end subroutine QuatRot

end module alturbine_utils
