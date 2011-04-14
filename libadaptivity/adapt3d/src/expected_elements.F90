! Copyright (C) 2006 Imperial College London and others.
!    
! Please see the AUTHORS file in the main source directory for a full list
! of copyright holders.
!
! Prof. C Pain
! Applied Modelling and Computation Group
! Department of Earth Science and Engineering
! Imperial College London
!
! C.Pain@Imperial.ac.uk
!    
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation,
! version 2.1 of the License.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
! USA

!< Given a metric tensor field and a mesh, calculate the number of
!< tetrahedral elements we would expect after adapting to that metric.
module predicted_elements
  implicit none
  private
  public::predicted_number_elements
  
contains
  
  integer function predicted_number_elements(Metric, X, Y, Z, NDGLNO,&
       NNodes, NTetra, nloc)    
    integer, intent(in)::NNodes
    real, intent(in)::Metric(NNodes*9), X(NNodes), Y(NNodes), Z(NNodes)
    integer, intent(in)::NTetra, nloc, NDGLNO(NTetra*NLOC)
    
    ! The quantity to scale the predicted no of elements 
    real, parameter::VOLSCA=1.2
    
    ! Volume of unit tetrahedral(VOL1TE) =1./sqrt(72). (each side has length 1.)
    real, parameter::VOL1EL=0.11785113
    
    real VOLUME, VOL, DET, MEANM(9), XX(4), YY(4), ZZ(4), MXDET, SMDET
    integer ELE, Node, I, ILOC
    
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
          write(0, *) "WARNING: Volume of element", ele, "not positive!"
          write(0, *) "WARNING: Bad news! Your inputs to adaptivity are fecked."
       end if
       CALL JACDET(DET, MEANM, 3)
       MXDET = MAX(MXDET,DET)
       SMDET = SMDET + DET
       
       VOLUME = VOLUME + VOL*SQRT(MAX(0.0,DET))
    END DO
    
    I = INT(VOLSCA*VOLUME/VOL1EL)
    IF( I .LT. 2 ) THEN
       write(0, *) 'WARNING: predicted number of elements looks wrong: ',I
       write(0, *) 'VOLSCA,VOLUME,VOL1EL = ',VOLSCA,VOLUME,VOL1EL
       write(0, *) 'Determinant max & sum = ',MXDET,SMDET
       write(0, *) 'NNOD,NELM,NLOC = ',NNodes,NTetra,NLOC
    END IF
    
    predicted_number_elements = I
  end function predicted_number_elements

  ! This sub calculates the volume of an element. 
  ! X,Y,Z are the coords of the 4 nodes.   
  subroutine jacvol(VOL, X, Y, Z, D3)
    real, intent(out)::vol
    real, intent(in)::X(4),Y(4),Z(4) 
    LOGICAL, intent(in)::D3

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
    
  end subroutine jacvol

  ! This sub calculates the volume of an element. 
  ! X,Y,Z are the coords of the 4 nodes.   
  subroutine jacdet(DET,M,NDIM)
    integer, intent(in)::NDIM
    real, intent(out)::DET
    real, intent(in)::M(NDIM,NDIM) 
    
    IF(NDIM.EQ.3) THEN
       DET = M(1,1)*( M(2,2)*M(3,3) - M(2,3)*M(3,2) )&
            &      + M(2,1)*( M(3,2)*M(1,3) - M(1,2)*M(3,3) )&
            &      + M(3,1)*( M(1,2)*M(2,3) - M(2,2)*M(1,3) )
    ELSE
       DET = M(1,1)*M(2,2)-M(1,2)*M(2,1) 
    ENDIF
    
  end subroutine jacdet
end module predicted_elements

! I don't feel up to hacking at at F90 name mangling for C
! interfacing, therefore:
integer function get_predicted_nelements(Metric, X, Y, Z, NDGLNO, &
     NNodes, NTetra, nloc)    
  use predicted_elements
  integer, intent(in)::NNodes
  real, intent(in)::Metric(NNodes*9), X(NNodes), Y(NNodes), Z(NNodes)
  integer, intent(in)::NTetra, nloc, NDGLNO(NTetra*NLOC)
  
  get_predicted_nelements = predicted_number_elements(Metric, X, Y, Z, NDGLNO, &
       NNodes, NTetra, nloc)    
  return
end function get_predicted_nelements
