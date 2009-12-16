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

module lagrot_module

  use rgptwe_module

  implicit none

  private

  public :: lagrot

contains

      SUBROUTINE LAGROT(WEIT,QUAPOS,NDGI,GETNDP)
!     This computes the weight and points for standard Gaussian quadrature.
!     IF(GETNDP) then get the POSITION OF THE NODES 
!     AND DONT BOTHER WITH THE WEITS.
      INTEGER NDGI
      REAL WEIT(NDGI),QUAPOS(NDGI)
      LOGICAL GETNDP
      LOGICAL WEIGHT
      INTEGER IG
!     
      IF(.NOT.GETNDP) THEN
         WEIGHT=.TRUE.
      do IG=1,NDGI
            WEIT(IG)=RGPTWE(IG,NDGI,WEIGHT)
         END DO
!     
         WEIGHT=.FALSE.
      do IG=1,NDGI
            QUAPOS(IG)=RGPTWE(IG,NDGI,WEIGHT)
         END DO
      ELSE
         IF(NDGI.EQ.1) THEN
            QUAPOS(1)=0.
         ELSE
      do IG=1,NDGI
               QUAPOS(IG)= -1+2.*REAL(IG-1)/REAL(NDGI-1)
            END DO
         ENDIF
      ENDIF

  end subroutine lagrot

end module lagrot_module
