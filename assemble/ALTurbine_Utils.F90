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
   
    integer function FindMinimum(x,Start,End)
    implicit none
    integer, dimension(1:),intent(in) :: x
    integer, intent(in)               :: Start, End
    integer                           :: Minimum
    integer                           :: Location
    integer                           :: i

    minimum = x(start)
    Location = Start 
    do i=start+1,End
        if(x(i) < Minimum) then
            Minimum = x(i)
            Location = i 
        end if
    end do
    FindMinimum = Location
    end function FindMinimum

    subroutine swap(a,b)
    implicit none
    integer, intent(inout) :: a,b
    integer                :: Temp

    Temp = a
    a = b
    b = Temp
    end subroutine swap

    subroutine sort(x,size)
    implicit none
    integer, dimension(1:), intent(INOUT) :: x
    integer, intent(in)                   :: size
    integer                               :: i 
    integer                               :: Location

    do i=1,Size-1
        location=FindMinimum(x,i,size)
        call swap(x(i),x(Location))
    end do
    end subroutine


 
end module alturbine_utils
