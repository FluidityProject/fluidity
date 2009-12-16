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

module CriticalTimeStep
  implicit none

  private
  public::GetCriticalDt
contains

  ! For a given mesh and velocity field, find the time step 
  ! that gives a maximum Courant number of unity. If the returned
  ! timestep is less than one then all velocities are zero and thus
  ! The time step can be infinity.
  SUBROUTINE GetCriticalDt(X, Y, Z, XONDGL, &
       Ux, Uy, Uz,                            &
       UG, VG, WG,                            &
       TOTELE, NONODS, NLOC, NGI,             &
       NLX, NLY, NLZ, N, WEIGHT,D3,             &
       CriticalDT,ISPHERE)
    use parallel_tools
    IMPLICIT NONE
    INTEGER TOTELE,NONODS,NLOC,NGI
    REAL X(NONODS), Y(NONODS), Z(NONODS)
    INTEGER XONDGL(TOTELE*NLOC)
    REAL Ux(NONODS),Uy(NONODS),Uz(NONODS)
    REAL UG(NONODS),VG(NONODS),WG(NONODS)
    
    ! Shape functions
    REAL NLX(NLOC,NGI), NLY(NLOC,NGI), NLZ(NLOC,NGI)
    
    ! Interpolation weights at quadrature points
    REAL WEIGHT(NGI)
    REAL N(NLOC,NGI)
    LOGICAL D3

    REAL CriticalDT
    INTEGER ISPHERE

    ! Jacobian matrix, J
    REAL AGI,BGI,CGI, DGI,EGI,FGI, GGI,HGI,KGI

    ! J^-1
    REAL A11,A12,A13, A21,A22,A23, A31,A32,A33
    REAL DETJ

    ! Velocity at a quadrature point
    REAL Vx, Vy, Vz, V

    ! Scale vector
    REAL Lx, Ly, Lz, L

    REAL DT, hmax, h2
    INTEGER N1,N2,I,J,GI,IGLX,ELE
    integer myrank, element_owner

    CriticalDT = -1.0
    myrank=GetRank()    

    ! Ensure that the helo is healty
    call flcomms_update(1, Ux, 1, 1, 0)
    call flcomms_update(1, Uy, 1, 1, 0)
    call flcomms_update(1, Uz, 1, 1, 0)

    call flcomms_update(1, UG, 1, 1, 0)
    call flcomms_update(1, VG, 1, 1, 0)
    call flcomms_update(1, WG, 1, 1, 0)
    
    DO ELE=1,TOTELE
       call flcomms_get_element_owner(2, ele, element_owner)         
       if(element_owner/=myrank) then
          cycle
       end if

       ! Find longest edge in the element squared.
       ! This is used to protect us from small V
       h2 = 0.0
       DO N1=1, NLOC
          I = XONDGL((ELE-1)*NLOC+N1)
          DO N2=1, NLOC
             J = XONDGL((ELE-1)*NLOC+N2)
             IF(D3) THEN
                h2 = MAX(h2, (X(I)-X(J))*(X(I)-X(J)) + (Y(I)-Y(J))*(Y(I)-Y(J)) + (Z(I)-Z(J))*(Z(I)-Z(J)))
             ELSE
                h2 = MAX(h2, (X(I)-X(J))*(X(I)-X(J)) + (Y(I)-Y(J))*(Y(I)-Y(J)))
             END IF
          END DO
       END DO
       hmax = SQRT(h2)

       DO GI=1,NGI
          IF(D3) THEN
             AGI=0.
             BGI=0.
             CGI=0.

             DGI=0.
             EGI=0.
             FGI=0.

             GGI=0.
             HGI=0.
             KGI=0.

             Vx=0.
             Vy=0.
             Vz=0.

             ! Calculate the Jacobian matrix for this quadriture point
             DO I=1,NLOC
                IGLX=XONDGL((ELE-1)*NLOC+I)

                ! NB R0 does not appear here although the z-coord might be Z+R0. 
                AGI=AGI+NLX(I,GI)*X(IGLX) 
                BGI=BGI+NLX(I,GI)*Y(IGLX) 
                CGI=CGI+NLX(I,GI)*Z(IGLX) 

                DGI=DGI+NLY(I,GI)*X(IGLX) 
                EGI=EGI+NLY(I,GI)*Y(IGLX) 
                FGI=FGI+NLY(I,GI)*Z(IGLX) 

                GGI=GGI+NLZ(I,GI)*X(IGLX) 
                HGI=HGI+NLZ(I,GI)*Y(IGLX) 
                KGI=KGI+NLZ(I,GI)*Z(IGLX)

                ! Interpolate velocity onto the quadrature point
                !              Vx = Vx + N(I,GI)*Ux(IGLX)
                !              Vy = Vy + N(I,GI)*Uy(IGLX)
                !              Vz = Vz + N(I,GI)*Uz(IGLX)
                Vx = Vx + N(I,GI)*(Ux(IGLX)-UG(IGLX))
                Vy = Vy + N(I,GI)*(Uy(IGLX)-VG(IGLX))
                Vz = Vz + N(I,GI)*(Uz(IGLX)-WG(IGLX))
             END DO
             V = SQRT(Vx*Vx + Vy*Vy + Vz*Vz)

             ! ewrite(3,*) ""
             ! ewrite(3,*) "Jacobian matrix"
             ! ewrite(3,*)  AGI, BGI, CGI
             ! ewrite(3,*)  DGI, EGI, FGI
             ! ewrite(3,*)  GGI, HGI, KGI
             ! ewrite(3,*) ""

             DETJ= AGI*(EGI*KGI-FGI*HGI) &
                  -BGI*(DGI*KGI-FGI*GGI) &
                  +CGI*(DGI*HGI-EGI*GGI)
             ! ewrite(3,*) "Jacobian =",DETJ

             ! For coefficient in the inverse matrix of the jacobian. 
             A11= (EGI*KGI-FGI*HGI)/DETJ
             A21=-(DGI*KGI-FGI*GGI)/DETJ
             A31= (DGI*HGI-EGI*GGI)/DETJ

             A12=-(BGI*KGI-CGI*HGI)/DETJ
             A22= (AGI*KGI-CGI*GGI)/DETJ
             A32=-(AGI*HGI-BGI*GGI)/DETJ

             A13= (BGI*FGI-CGI*EGI)/DETJ
             A23=-(AGI*FGI-CGI*DGI)/DETJ
             A33= (AGI*EGI-BGI*DGI)/DETJ

             ! ewrite(3,*) ""
             ! ewrite(3,*) "Jacobian matrix^-1"
             ! ewrite(3,*)  A11, A21, A31
             ! ewrite(3,*)  A12, A22, A32
             ! ewrite(3,*)  A13, A23, A33
             ! ewrite(3,*) ""

             Lx = A11*Vx + A12*Vy + A13*Vz
             Ly = A21*Vx + A22*Vy + A23*Vz
             Lz = A31*Vx + A32*Vy + A33*Vz
             ! Lx = AGI*Vx + BGI*Vy + CGI*Vz
             ! Ly = DGI*Vx + EGI*Vy + FGI*Vz
             ! Lz = GGI*Vx + HGI*Vy + KGI*Vz

             L  = SQRT(Lx*Lx + Ly*Ly + Lz*Lz)

          ELSE

             AGI=0.
             BGI=0.
             CGI=0.
             DGI=0.

             DO I=1,NLOC
                IGLX=XONDGL((ELE-1)*NLOC+I)

                AGI=AGI + NLX(I,GI)*X(IGLX) 
                BGI=BGI + NLX(I,GI)*Y(IGLX) 
                CGI=CGI + NLY(I,GI)*X(IGLX) 
                DGI=DGI + NLY(I,GI)*Y(IGLX) 

                ! Interpolate velocity onto the quadrature point
                !              Vx = Vx + N(I,GI)*Ux(IGLX)
                !              Vy = Vy + N(I,GI)*Uy(IGLX)
                Vx = Vx + N(I,GI)*(Ux(IGLX)-UG(IGLX))
                Vy = Vy + N(I,GI)*(Uy(IGLX)-VG(IGLX))
             END DO
             V = SQRT(Vx*Vx + Vy*Vy)

             DETJ= AGI*DGI-BGI*CGI

             A11 = DGI/DETJ
             A12 =-BGI/DETJ
             A21 =-CGI/DETJ
             A22 = AGI/DETJ

             Lx = A11*Vx + A12*Vy
             Ly = A21*Vx + A22*Vy
             L = SQRT(Lx*Lx + Ly*Ly)

          END IF

          IF(V.GT.0.0) THEN
             L = V/L
             L = MIN(L, hmax)

             ! DT that would give a Courant number of unity
             DT = L/V

             IF(CriticalDT.LT.0.0) THEN
                CriticalDT = DT
             ELSE
                CriticalDT = MIN(CriticalDT, DT)
             END IF

          END IF
       END DO ! DO GI=1, NGI
    END DO

    CALL ALLMIN(CriticalDT)

    RETURN
  END SUBROUTINE GetCriticalDt
  !
end module CriticalTimeStep
