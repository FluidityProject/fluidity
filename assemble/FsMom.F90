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

module fsmom_module

  implicit none
  
  private
  
  !public :: 
  
contains

        SUBROUTINE FsMom(NONODS,XNONOD,FREDOP,STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI, &
                       FSNDID,SOURCX,SOURCY,SOURCZ,SNDGLN,X,Y,Z,H,              &
                       SN,SNLX,SNLY,SWEIGH,D3)                      

      LOGICAL  D3
      INTEGER  STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI
      INTEGER  NONODS,XNONOD,FREDOP
      INTEGER  SNDGLN(STOTEL*SNLOC)
      REAL     FSNDID(NONODS)
      REAL     X(XNONOD),Y(XNONOD),Z(XNONOD),H(XNONOD)
      REAL     ZD(NGI)
      REAL     SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
      REAL,DIMENSION(:,:),ALLOCATABLE  :: SNX,SNY,SNZ
      REAL     SWEIGH(SNGI)
      REAL     SOURCX(NONODS),SOURCY(NONODS),SOURCZ(NONODS)
! Working arrays/local variables
      REAL     AGI,BGI,CGI,DGI
      REAL     RBIGM2X,RBIGM2Y,RBIGM1,RNN
      INTEGER  SELE,ILOC,JLOC,SILOC,SJLOC,L,IGL,IGLT,GLOBI,GLOBJ,GLOBIP,GLOBJP,GI
      INTEGER  HIT,NOD,NODI,NODJ,ITNOD,JTNOD,INUM,UPPER,LOWER,COUNT
      INTEGER  COUNTER
      REAL     DETWEI(NGI),DET
      LOGICAL  YES
      REAL,PARAMETER::  WSF=1.0  !0.1509   !WSF = g/H, g is gravity, H is average height
      REAL,PARAMETER::  H0 = 0.0 !, Average water depth = 65.0

      ALLOCATE( SNX(SNLOC,SNGI) )
      ALLOCATE( SNY(SNLOC,SNGI) )
      ALLOCATE( SNZ(SNLOC,SNGI) )
      
! Counter counts how many nodes are found on the free surface
       COUNTER = 0         
   DO  SELE=1,STOTEL
! NB. These are the surface elements
         HIT=0
         DO SILOC=1,SNLOC
           IGL = SNDGLN((SELE-1)*SNLOC + SILOC)
           IF(FSNDID(IGL).GT.0.1) HIT=HIT+1      
!           ewrite(3,*) FSNDID(IGL),X(IGL),Y(IGL),Z(IGL)    
         ENDDO
!
!C If HIT=SNLOC then we are on the free surface so continue
     IF(HIT.EQ.SNLOC) THEN
!C SILOC=1 OK here as this is a surfac rather than volume element.
           COUNTER = COUNTER +1

         DO GI = 1,SNGI 
            AGI = 0.
            BGI = 0.
            CGI = 0.
            DGI = 0.
            ZD(GI) = 0.
            DO L = 1,SNLOC
              IGL = SNDGLN((SELE-1)*SNLOC+L)

              AGI = AGI + SNLX(L,GI)*X(IGL)
              BGI = BGI + SNLX(L,GI)*Y(IGL)
              CGI = CGI + SNLY(L,GI)*X(IGL)
              DGI = DGI + SNLY(L,GI)*Y(IGL)

              ZD(GI) = ZD(GI) + SN(L,GI)*Z(IGL)
            ENDDO

            DET = AGI*DGI-BGI*CGI
            DETWEI(GI) = DET*SWEIGH(GI)
            DO L=1,SNLOC
              SNX(L,GI) = (DGI*SNLX(L,GI)-BGI*SNLY(L,GI))/DET
              SNY(L,GI) = (-CGI*SNLX(L,GI)+AGI*SNLY(L,GI))/DET
              SNZ(L,GI) = 0.0
            END DO 

         ENDDO !end gi           
!C----------------------------- 
         DO ILOC = 1,SNLOC
                 GLOBI = SNDGLN((SELE-1)*SNLOC+ILOC)
             DO JLOC = 1,SNLOC
                  GLOBJ = SNDGLN((SELE-1)*SNLOC+JLOC)
                  RBIGM2X = 0.
                  RBIGM2Y = 0.
                  RBIGM1 = 0.
                  DO GI = 1,SNGI
                      RNN    = SN(ILOC,GI)*SN(JLOC,GI)*DETWEI(GI)                    
! RBIGM2  cf  VLND in diff2d.f
                    RBIGM2X = RBIGM2X + DETWEI(GI)*ZD(GI)*SN(ILOC,GI)*SNX(JLOC,GI)
                    RBIGM2Y = RBIGM2Y + DETWEI(GI)*ZD(GI)*SN(ILOC,GI)*SNY(JLOC,GI)
                    RBIGM1 = RBIGM1 + RNN
                  ENDDO
!cc source term on RHS
               SOURCX(GLOBI) = SOURCX(GLOBI) +RBIGM2X*H(GLOBJ)
               SOURCY(GLOBI) = SOURCY(GLOBI) +RBIGM2Y*H(GLOBJ)
               SOURCZ(GLOBI) = SOURCZ(GLOBI) +RBIGM1*H(GLOBJ)
             ENDDO

         ENDDO
     ENDIF

  ENDDO  !end sele

      DEALLOCATE( SNX )
      DEALLOCATE( SNY )
      DEALLOCATE( SNZ )

! Take into accouting the dimension in Eqs. We have to multiply source by WSF

   DO NOD = 1, NONODS
               SOURCX(GLOBI) = SOURCX(NOD)*WSF
               SOURCY(GLOBI) = SOURCY(NOD)*WSF
               SOURCZ(GLOBI) = SOURCZ(NOD)*WSF
   ENDDO


       END SUBROUTINE FsMom

end module fsmom_module
