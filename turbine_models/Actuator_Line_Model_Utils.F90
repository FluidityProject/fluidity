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
  
module actuator_line_model_utils

    use fldebug
    use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, PYTHON_FUNC_LEN, pi

    implicit none

    public QuatRot, cross, IsoKernel, AnIsoKernel

contains 
 
    SUBROUTINE cross(ax,ay,az,bx,by,bz,cx,cy,cz) 

        real ax,ay,az,bx,by,bz,cx,cy,cz 

        cx = ay*bz - az*by
        cy = az*bx - ax*bz
        cz = ax*by - ay*bx

    End SUBROUTINE cross

   subroutine QuatRot(vx,vy,vz,Theta,Rx,Ry,Rz,Ox,Oy,Oz,vRx,vRy,vRz)
   
   ! % Perform rotation of vector v around normal vector nR using the
   ! % quaternion machinery.
   ! % v: input vector
   ! % Theta: rotation angle (rad)
   ! % nR: normal vector around which to rotate
   ! % Origin: origin point of rotation
   ! %
   ! % vR: Rotated vector

   implicit none
   real,intent(in) :: vx,vy,vz,Theta,Rx,Ry,Rz,Ox,Oy,Oz
   real,intent(inout):: vRx,vRy,vRz     
   real :: nRx,nRy,nRz 
   real :: p(4,1), pR(4,1), q(4), qbar(4), RMag, vOx, vOy, vOz
   real :: QL(4,4), QbarR(4,4)
    
   ! Force normalize nR
   RMag=sqrt(Rx**2.0+Ry**2.0+Rz**2.0)
   nRx=Rx/RMag
   nRy=Ry/RMag
   nRz=Rz/RMag

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

    subroutine IDW(Ncol,Xcol,Ycol,Zcol,Fxcol,Fycol,Fzcol,p,Xmesh,Ymesh,Zmesh,Fxmesh,Fymesh,Fzmesh)
        implicit none
        integer, intent(in) :: Ncol
        real, dimension(Ncol),intent(in) :: Xcol,Ycol,Zcol,Fxcol,Fycol,Fzcol 
        real, intent(in) :: Xmesh,Ymesh,Zmesh
        integer,intent(in) :: p
        real, intent(inout) :: Fxmesh,Fymesh,Fzmesh
        
        real,dimension(Ncol) :: d(Ncol), w(Ncol)
        real ::  wsum
        integer :: i,imin

        wsum=0.0
        do i=1,Ncol     
        d(i)=sqrt((Xcol(i)-Xmesh)**2+(Ycol(i)-Ymesh)**2+(Zcol(i)-Zmesh)**2)
        w(i)=1/d(i)**p
        wsum=wsum+w(i)
        end do
        
        if (minval(d)<0.001) then
            imin=minloc(d,1)
            Fxmesh=Fxcol(imin)
            Fymesh=Fycol(imin)
            Fzmesh=Fzcol(imin)
        else
            Fxmesh=0.0
            Fymesh=0.0
            Fzmesh=0.0
            do i=1,Ncol
            Fxmesh=Fxmesh+w(i)*Fxcol(i)/wsum
            Fymesh=Fymesh+w(i)*Fycol(i)/wsum
            Fzmesh=Fzmesh+w(i)*Fzcol(i)/wsum
            enddo
        endif

    end subroutine IDW

    real function IsoKernel(dr,mesh_size,chord,dim)
    
        implicit none
        integer,intent(in) :: dim
        real,intent(in) ::dr, mesh_size,chord
        real :: epsilon_par,epsilon_threshold,epsilon_chord
        integer :: j
    
          
           epsilon_threshold = 2.0*mesh_size
           epsilon_chord = 0.4*chord
           
           if(epsilon_chord >= epsilon_threshold) then
              epsilon_par=epsilon_chord
           else
              epsilon_par=epsilon_threshold
           end if
               

            if(dim==2) then    
            IsoKernel = 1.0/(epsilon_par**2*pi)*exp(-(dr/epsilon_par)**2.0)
            elseif(dim==3) then
            IsoKernel = 1.0/(epsilon_par**3.0*pi**1.5)*exp(-(dr/epsilon_par)**2.0)
            else 
            FLAbort("1D source not implemented")
            endif
        
    end function IsoKernel
    
    real function AnIsoKernel(dx,dy,dz,nx,ny,nz,tx,ty,tz,sx,sy,sz,ec,et,es)
    
        implicit none
        real,intent(in) :: dx,dy,dz,nx,ny,nz,tx,ty,tz,sx,sy,sz,ec,et,es
        real :: n,t,s

        n=dx*nx+dy*ny+dz*nz ! normal projection
        t=dx*tx+dy*ty+dz*tz ! Chordwise projection
        s=dx*sx+dy*sy+dz*sz ! Spanwise projection

        if(abs(s)<=es) then
        AnIsoKernel = exp(-((n/et)**2.0+(t/ec)**2.0))/(ec*et*pi)
        else
        AnIsoKernel = 0.0
        endif
    
    end function AnIsoKernel

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

  function dirname(number)
    integer,intent(in)    :: number
    character(len=6)  :: dirname

    ! Cast the (rounded) number to string using 6 digits and
    ! leading zeros
    write (dirname, '(I6.1)')  number
    ! This is the same w/o leading zeros  
    !write (dirname, '(I6)')  nint(number)

    ! This is for one digit (no rounding)
    !write (dirname, '(F4.1)')  number
  end function

 
end module actuator_line_model_utils
