!     Copyright (C) 2006 Imperial College London and others.
!     
!     Please see the AUTHORS file in the main source directory for a full list
!     of copyright holders.
!     
!     Prof. C Pain
!     Applied Modelling and Computation Group
!     Department of Earth Science and Engineering
!     Imperial College London
!     
!     C.Pain@Imperial.ac.uk
!     
!     This library is free software; you can redistribute it and/or
!     modify it under the terms of the GNU Lesser General Public
!     License as published by the Free Software Foundation,
!     version 2.1 of the License.
!     
!     This library is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!     Lesser General Public License for more details.
!     
!     You should have received a copy of the GNU Lesser General Public
!     License along with this library; if not, write to the Free Software
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!     USA

#include "fdebug.h"

module diagnostic_body_drag_module

  use allsorts
  use fldebug
  use parallel_tools
  use surface_integration
  
  implicit none
  
  private
  
  !public :: 
  
contains
   
      SUBROUTINE diagnostic_body_drag(FX,FY,FZ,ELE_SURFACE_FLAG,X,Y,Z,U,V,W,P,VISC,&
     & NORMALX,NORMALY,NORMALZ,&
     & NONODS,TOTELE,STOTEL,NDGLNO,XONDGL,SNDGLN,&
     & NLOC,NGI,SNLOC,SNGI,&
     & N,NLX,NLY,NLZ,WEIGHT,DETWEI,&
     & NX,NY,NZ,&
     & SN,SNLX,SNLY,SWEIGH,&
     & FINDRM,COLM,NCOLM,D3,DCYL,&
     & TAUXX,TAUXY,TAUXZ,&
     & TAUYY,TAUYZ,TAUZZ)
      INTEGER  NONODS,STOTEL,NLOC,SNLOC,NGI,SNGI,TOTELE
      INTEGER  ELE_SURFACE_FLAG(STOTEL)
      REAL     X(NONODS),Y(NONODS),Z(NONODS)
      REAL     U(NONODS),V(NONODS),W(NONODS),P(NONODS),VISC,TWOVISC
      REAL     NORMALX(NONODS),NORMALY(NONODS),NORMALZ(NONODS)
      REAL     N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
      REAL     WEIGHT(NGI),DETWEI(NGI),SWEIGH(SNGI),VOLUME
      REAL     SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
      REAL     NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
      INTEGER  FINDRM(NONODS+1),NCOLM,COLM(NCOLM),NDGLNO(TOTELE*NLOC),XONDGL(TOTELE*NLOC),SNDGLN(STOTEL*SNLOC)
      REAL     TAUXX(TOTELE),TAUXY(TOTELE),TAUXZ(TOTELE),TAUYY(TOTELE),TAUYZ(TOTELE),TAUZZ(TOTELE)
      REAL     FX,FY,FZ,CD,CL
      LOGICAL  D3,DCYL
      INTEGER  SELE,ELE,VELE,ILOC,JLOC,SILOC,SJLOC,GLOBI,GLOBJ,SGLOBI,GI,COUNT,L,IGL,COUNTER,NOD1,NOD2,NOD3
      REAL     DXDLX,DXDLY,DYDLX,DYDLY,DZDLX,DZDLY,A,B,C,RNN
      REAL     NORMALXSELE,NORMALYSELE,NORMALZSELE,XA,XB,YA,YB,ZA,ZB,AREA1,ABSA,ABSB,SINTET,COSTET,TOTAL_AREA,total_area2
      REAL     norm12,norm23,norm31,sumlen,radcheck
      integer myrank, element_owner

      myrank=GetRank()    

! Calculate the element-wise entries of the stress tensor
      ewrite(2,*) 'In diagnostic_body_drag'
      IF(SNLOC.NE.3) THEN
        ewrite(-1,*) 'This subroutine makes an assumption that ',&
     &  'requires linear tets so cannot be used here.'
        STOP 2304
      ENDIF

      total_area = 0.0
! Make some assumptions
      TWOVISC = 2.0*VISC

!      ewrite(3,*) 'SWEIGH:',SWEIGH
!      ewrite(3,*) 'VISC:',VISC
!      ewrite(3,*) 'normalx',r2norm(normalx,nonods,0)
!      ewrite(3,*) 'u',r2norm(u,nonods,0)
!      ewrite(3,*) 'P',r2norm(p,nonods,0)
      
      TAUXX(1:TOTELE) = 0.0
      TAUXY(1:TOTELE) = 0.0
      TAUXZ(1:TOTELE) = 0.0
      TAUYY(1:TOTELE) = 0.0
      TAUYZ(1:TOTELE) = 0.0
      TAUZZ(1:TOTELE) = 0.0
    
      TOTAL_AREA = 0.0      
      total_area2 = 0.0      

! Calculate the drag and lift forces and coefficients.

      FX = 0.0
      FY = 0.0
      FZ = 0.0      
      counter = 0      
      do SELE = 1,STOTEL
        IF(ELE_SURFACE_FLAG(SELE).EQ.1) THEN ! On the surface of interest

! WHAT IS CORRESPONDING VOLUME ELE? ! THIS CAN BE DONE BETTER
          VELE = 0
      do ELE = 1,TOTELE
            COUNT = 0
      do SILOC = 1,SNLOC
            SGLOBI = SNDGLN((SELE-1)*SNLOC + SILOC)
      do ILOC = 1,NLOC
                GLOBI = NDGLNO((ELE-1)*NLOC + ILOC)
                IF(GLOBI.EQ.SGLOBI) COUNT = COUNT + 1
              END DO
            END DO 
            IF(COUNT.EQ.SNLOC) THEN
               VELE = ELE
               EXIT
            ENDIF             
          END DO
          IF(VELE.EQ.0) THEN
            STOP 5824    
          ENDIF
            ! Check if this processor should be integrating over this element - avoid double-counting
            call flcomms_get_element_owner(2, vele, element_owner)         
            if(element_owner/=myrank) then
               cycle
            end if
            counter = counter+1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Compute grad of viscous tensor on this volume element
!
        CALL DETNLX(vELE, X,Y,Z, XONDGL, TOTELE,NONODS,NLOC,NGI, &
     &        N,NLX,NLY,NLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL, &
     &        NX,NY,NZ)
      do ILOC = 1,NLOC
            GLOBI = NDGLNO((vELE-1)*NLOC + ILOC)
            TAUXX(vELE) = TAUXX(vELE) + U(GLOBI)*NX(ILOC,1)   ! Assuming linear tet basis funcions and hence constant derivs 
            TAUXY(vELE) = TAUXY(vELE) + U(GLOBI)*NY(ILOC,1) + V(GLOBI)*NX(ILOC,1)
            TAUXZ(vELE) = TAUXZ(vELE) + U(GLOBI)*NZ(ILOC,1) + W(GLOBI)*NX(ILOC,1)
            TAUYY(vELE) = TAUYY(vELE) + V(GLOBI)*NY(ILOC,1)
            TAUYZ(vELE) = TAUYZ(vELE) + V(GLOBI)*NZ(ILOC,1) + W(GLOBI)*NY(ILOC,1)
            TAUZZ(vELE) = TAUZZ(vELE) + W(GLOBI)*NZ(ILOC,1)
          END DO
          TAUXX(vELE) = TWOVISC*TAUXX(vELE)
          TAUYY(vELE) = TWOVISC*TAUYY(vELE)
          TAUZZ(vELE) = TWOVISC*TAUZZ(vELE)
          TAUXY(vELE) =    VISC*TAUXY(vELE)
          TAUXZ(vELE) =    VISC*TAUXZ(vELE)
          TAUYZ(vELE) =    VISC*TAUYZ(vELE)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      do GI = 1,SNGI
            DXDLX = 0.
            DXDLY = 0.
            DYDLX = 0.
            DYDLY = 0.
            DZDLX = 0.
            DZDLY = 0.       
      do L = 1,SNLOC
              IGL   = SNDGLN((SELE-1)*SNLOC+L)
              DXDLX = DXDLX + SNLX(L,GI)*X(IGL)
              DXDLY = DXDLY + SNLY(L,GI)*X(IGL) 
              DYDLX = DYDLX + SNLX(L,GI)*Y(IGL)
              DYDLY = DYDLY + SNLY(L,GI)*Y(IGL)
              DZDLX = DZDLX + SNLX(L,GI)*Z(IGL) 
              DZDLY = DZDLY + SNLY(L,GI)*Z(IGL)
            END DO
            A = DYDLX*DZDLY - DYDLY*DZDLX
            B = DXDLX*DZDLY - DXDLY*DZDLX
            C = DXDLX*DYDLY - DXDLY*DYDLX      
            DETWEI(GI) = SQRT( A**2 + B**2 + C**2 )*SWEIGH(GI)        
            total_area = total_area+DETWEI(GI)
          END DO
!          
! lets work out a surface element normal

!          NOD1 = SNDGLN((SELE-1)*SNLOC+1)
!          NOD2 = SNDGLN((SELE-1)*SNLOC+2)
!          NOD3 = SNDGLN((SELE-1)*SNLOC+3)

!          normalxsele = (normalx(NOD1)+normalx(NOD2)+normalx(NOD3))/3.
!          normalysele = (normaly(NOD1)+normaly(NOD2)+normaly(NOD3))/3.
!          normalzsele = (normalz(NOD1)+normalz(NOD2)+normalz(NOD3))/3.

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! and the area of this triangle - 
! need to uncomment if using the p/w constant on an element surface for P or not including P in integrals
!          XB=X(NOD2)-X(NOD1)
!          XA=X(NOD3)-X(NOD1)
!          YB=Y(NOD2)-Y(NOD1)
!          YA=Y(NOD3)-Y(NOD1)
!          ZB=Z(NOD2)-Z(NOD1)
!          ZA=Z(NOD3)-Z(NOD1)
!          ABSA=SQRT(XA**2+YA**2+ZA**2)
!          ABSB=SQRT(XB**2+YB**2+ZB**2)
!          COSTET=((XB*XA)+(YB*YA)+(ZB*ZA))/(ABSA*ABSB)
!          SINTET=SQRT(1-COSTET**2)
!          AREA1=0.5*ABSA*ABSB*SINTET
!          TOTAL_AREA = TOTAL_AREA + AREA1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Alternative method - more expensive than above?
!          norm12 = sqrt( (x(nod1)-x(nod2))**2 + (y(nod1)-y(nod2))**2 + (z(nod1)-z(nod2))**2 )
!          norm23 = sqrt( (x(nod2)-x(nod3))**2 + (y(nod2)-y(nod3))**2 + (z(nod2)-z(nod3))**2 )
!          norm31 = sqrt( (x(nod3)-x(nod1))**2 + (y(nod3)-y(nod1))**2 + (z(nod3)-z(nod1))**2 )
!          sumlen = (norm12 + norm23 + norm31)/2.0
!          area1 = sqrt(sumlen*(sumlen-norm12)*(sumlen-norm23)*(sumlen-norm31))
!          total_area2 = total_area2+area1
!


      do SILOC = 1,SNLOC
            GLOBI = SNDGLN((SELE-1)*SNLOC + SILOC)
      do SJLOC = 1,SNLOC
              GLOBJ = SNDGLN((SELE-1)*SNLOC + SJLOC)
              RNN = 0.0              
      do GI = 1,SNGI
                RNN = RNN + SN(SILOC,GI)*SN(SJLOC,GI)*DETWEI(GI)
              END DO

               FX = FX + (NORMALX(GLOBJ)*TAUXX(VELE) + NORMALY(GLOBJ)*TAUXY(VELE) + NORMALZ(GLOBJ)*TAUXZ(VELE) -&
     &                    NORMALX(GLOBJ)*P(GLOBJ))*RNN
               FY = FY + (NORMALX(GLOBJ)*TAUXY(VELE) + NORMALY(GLOBJ)*TAUYY(VELE) + NORMALZ(GLOBJ)*TAUYZ(VELE) - &
     &                    NORMALY(GLOBJ)*P(GLOBJ))*RNN
               FZ = FZ + (NORMALX(GLOBJ)*TAUXZ(VELE) + NORMALY(GLOBJ)*TAUYZ(VELE) + NORMALZ(GLOBJ)*TAUZZ(VELE) - &
     &                    NORMALZ(GLOBJ)*P(GLOBJ))*RNN

!               FX = FX + (NORMALXsele*TAUXX(VELE) + NORMALYsele*TAUXY(VELE) + NORMALZsele*TAUXZ(VELE) -
!     &                    NORMALXsele*P(GLOBJ))*RNN
!               FY = FY + (NORMALXsele*TAUXY(VELE) + NORMALYsele*TAUYY(VELE) + NORMALZsele*TAUYZ(VELE) - 
!     &                    NORMALYsele*P(GLOBJ))*RNN
!               FZ = FZ + (NORMALXsele*TAUXZ(VELE) + NORMALYsele*TAUYZ(VELE) + NORMALZsele*TAUZZ(VELE) - 
!     &                    NORMALZsele*P(GLOBJ))*RNN

!               FX = FX - (NORMALXsele*P(GLOBJ))*RNN
!               FY = FY - (NORMALYsele*P(GLOBJ))*RNN
!               FZ = FZ - (NORMALZsele*P(GLOBJ))*RNN
!               FX = FX - (NORMALX(globj)*P(GLOBJ))*RNN
!               FY = FY - (NORMALY(globj)*P(GLOBJ))*RNN
!               FZ = FZ - (NORMALZ(globj)*P(GLOBJ))*RNN

            END DO
          END DO
        ENDIF               
      END DO

      call allsum(total_area)
      call allsum(fx)
      call allsum(fy)
      call allsum(fz)    
      ewrite(2,*) 'integrated over this many surface elements:',counter
      ewrite(2,*) 'total area integrated over:',total_area
      ewrite(2,*) 'fx,fy,fz:',fx,fy,fz

  end subroutine diagnostic_body_drag
    
      SUBROUTINE ON_BODY(X,Y,Z,NOD_NORMAL_X,NOD_NORMAL_Y,NOD_NORMAL_Z,ELE_SURFACE_FLAG,&
     & NONODS,TOTELE,STOTEL,NLOC,SNLOC,SNGI,NDGLNO,SNDGLN,TSNDGL,FLAG_NODES,GOT_FLAG_NODES,&
     & FINDRM,COLM,NCOLM,SALPHE,SUFNOD,SN,SNLX,SNLY,SWEIGH,IDFIELD,&
     & NOBCT,BCT1,BCT2,NNODP,halo_tag)
      INTEGER NONODS,TOTELE,STOTEL,NLOC,SNLOC,SNGI,SUFNOD,ELE_SURFACE_FLAG(STOTEL)      
      REAL    X(NONODS),Y(NONODS),Z(NONODS),SALPHE(SUFNOD)
      REAL    NOD_NORMAL_X(NONODS),NOD_NORMAL_Y(NONODS),NOD_NORMAL_Z(NONODS), FLAG_NODES(NONODS)
      REAL    SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI),SWEIGH(SNGI)
      LOGICAL GOT_FLAG_NODES
      INTEGER NDGLNO(NLOC*TOTELE),SNDGLN(SNLOC*STOTEL),TSNDGL(SNLOC*STOTEL)
      REAL    RADIUS,XC,YC,ZC,NORM,TOL,RADSQ,VOLELE,AGI,BGI,CGI,DGI,RN,FX,FY,FZ,AX,AY,AZ,DETJ
      REAL     DXDLX,DXDLY,DYDLX,DYDLY,DZDLX,DZDLY,A,B,C
      INTEGER COUNT,COUNTER,GLOBI,SELE,SILOC,IGL,IGLT,HIT,HIT2,GI,IGLX,I1,I2,I3,NOD,L
      INTEGER NCOLM,FINDRM(NONODS+1),COLM(NCOLM)
      REAL IDFIELD(NONODS)
      INTEGER  NOBCT,BCT2(NOBCT),NNODP,halo_tag
      REAL     BCT1(NOBCT)      

      INTEGER, ALLOCATABLE, DIMENSION(:)::NOD_SURFACE_FLAG
      REAL, ALLOCATABLE, DIMENSION(:)::ELE_NORMAL_X
      REAL, ALLOCATABLE, DIMENSION(:)::ELE_NORMAL_Y
      REAL, ALLOCATABLE, DIMENSION(:)::ELE_NORMAL_Z                  
      REAL, ALLOCATABLE, DIMENSION(:)::ELE_AREA
      REAL, ALLOCATABLE, DIMENSION(:)::AREA_SUM

           
      ALLOCATE(NOD_SURFACE_FLAG(NONODS))       
      ALLOCATE(ELE_NORMAL_X(STOTEL))       
      ALLOCATE(ELE_NORMAL_Y(STOTEL))     
      ALLOCATE(ELE_NORMAL_Z(STOTEL))       
      ALLOCATE(ELE_AREA(STOTEL))
      ALLOCATE(AREA_SUM(NONODS))


      ewrite(3,*) 'In ON_BODY',SNGI,SWEIGH,GOT_FLAG_NODES
      IDFIELD(1:NONODS)=0.0

      ewrite(3,*) 'BCT',nobct
      do nod=1,nobct
        idfield(BCT2(NOD))=1.0 ! Now using this rather than SALPHE which was corrupt for gid mesh for some reason
      enddo
      
                         
      ELE_SURFACE_FLAG(1:STOTEL) = 0
      NOD_SURFACE_FLAG(1:NONODS) = 0
      ELE_NORMAL_X(1:STOTEL) = 0.0
      ELE_NORMAL_Y(1:STOTEL) = 0.0
      ELE_NORMAL_Z(1:STOTEL) = 0.0
      NOD_NORMAL_X(1:NONODS) = 0.0
      NOD_NORMAL_Y(1:NONODS) = 0.0
      NOD_NORMAL_Z(1:NONODS) = 0.0
      ELE_AREA(1:STOTEL) = 0.0
      AREA_SUM(1:NONODS) = 0.0
      IF(GOT_FLAG_NODES) FLAG_NODES(1:NONODS) = 0.0  
      COUNTER = 0
          
! calculate surface nodes on body from BC defined in gem
      do SELE=1,STOTEL
        HIT=0
        HIT2=0
      do SILOC=1,SNLOC
!c          IGLT = TSNDGL((SELE-1)*SNLOC + SILOC)  
!c          IF(SALPHE(IGLT).GT.0.9999)  HIT = HIT+1                         
           IGL = SNDGLN((SELE-1)*SNLOC + SILOC)  
           IF( (idfield(IGL).GT.0.999) )  HIT = HIT+1
        ENDDO
        IF(HIT.EQ.SNLOC) THEN
!          ewrite(3,*) 'SELE:',SELE
          ELE_SURFACE_FLAG(SELE) = 1
      do SILOC=1,SNLOC
            IGL = SNDGLN((SELE-1)*SNLOC + SILOC)  
            NOD_SURFACE_FLAG(IGL) = 1                   
          ENDDO          
! Normal to this surface
! First calculate an approx normal (normx,normy,normz)...
          CALL SIMPNORM(FX,FY,FZ,.true.,&
     &           SNDGLN,STOTEL,SNLOC,NONODS,NONODS,SELE,&
     &           X,Y,Z,NCOLM,FINDRM,COLM) 
          I1=SNDGLN((SELE-1)*SNLOC+1)
          I2=SNDGLN((SELE-1)*SNLOC+2)
          I3=SNDGLN((SELE-1)*SNLOC+3)
! normal to two vectors
          CALL ANOFAC(I1,I2,I3, X,Y,Z,1.0, &
     &          AX,AY,AZ,.TRUE.,NONODS,-FX,-FY,-FZ)
          RN=SQRT(AX**2+AY**2+AZ**2)
          ELE_NORMAL_X(SELE) = AX/RN
          ELE_NORMAL_Y(SELE) = AY/RN
          ELE_NORMAL_Z(SELE) = AZ/RN
!          ewrite(3,*) 'Face Normal:',ELE_NORMAL_X(SELE),ELE_NORMAL_Y(SELE),ELE_NORMAL_Z(SELE)
! Area of this surface element
                  VOLELE = 0.0
      do GI = 1,SNGI
                     DXDLX  = 0.
                     DXDLY  = 0.
                     DYDLX  = 0.
                     DYDLY  = 0.
                     DZDLX  = 0.
                     DZDLY  = 0.           
      do L = 1,SNLOC
                        IGL   = SNDGLN((SELE-1)*SNLOC+L)
                        DXDLX = DXDLX + SNLX(L,GI)*X(IGL)
                        DXDLY = DXDLY + SNLY(L,GI)*X(IGL) 
                        DYDLX = DYDLX + SNLX(L,GI)*Y(IGL) 
                        DYDLY = DYDLY + SNLY(L,GI)*Y(IGL) 
                        DZDLX = DZDLX + SNLX(L,GI)*Z(IGL) 
                        DZDLY = DZDLY + SNLY(L,GI)*Z(IGL)
                     ENDDO
                     A = DYDLX*DZDLY - DYDLY*DZDLX
                     B = DXDLX*DZDLY - DXDLY*DZDLX
                     C = DXDLX*DYDLY - DXDLY*DYDLX
                     DETJ = SQRT( A**2 + B**2 + C**2 )
                     VOLELE = VOLELE + DETJ*SWEIGH(GI)                     
!            ewrite(3,*) 'VOLELE,DETJ,SWEIGH(GI):',VOLELE,DETJ,SWEIGH(GI)
                   ENDDO
          ELE_AREA(SELE) = VOLELE
!          ewrite(3,*) 'SELE,ELE_AREA(SELE):',SELE,ELE_AREA(SELE)
        ENDIF        
      ENDDO
! Now project the surface normals to the nodes using the weighting given by the areas of the surface elements
      do SELE=1,STOTEL
        IF(ELE_SURFACE_FLAG(SELE).EQ.1) THEN
!          ewrite(3,*) 'SELE:',SELE
          COUNTER = COUNTER+1        
      do SILOC=1,SNLOC
            IGL = SNDGLN((SELE-1)*SNLOC + SILOC)  
            IF(GOT_FLAG_NODES) FLAG_NODES(IGL) = 1.0      
            NOD_NORMAL_X(IGL) = NOD_NORMAL_X(IGL) + ELE_NORMAL_X(SELE)*ELE_AREA(SELE)
            NOD_NORMAL_Y(IGL) = NOD_NORMAL_Y(IGL) + ELE_NORMAL_Y(SELE)*ELE_AREA(SELE)
            NOD_NORMAL_Z(IGL) = NOD_NORMAL_Z(IGL) + ELE_NORMAL_Z(SELE)*ELE_AREA(SELE)
            AREA_SUM(IGL) = AREA_SUM(IGL) + ELE_AREA(SELE)
          ENDDO
        ENDIF
      END DO
      do NOD=1,NNODP
        IF(NOD_SURFACE_FLAG(NOD).EQ.1) THEN
          NOD_NORMAL_X(NOD) = NOD_NORMAL_X(NOD)/AREA_SUM(NOD)    
          NOD_NORMAL_Y(NOD) = NOD_NORMAL_Y(NOD)/AREA_SUM(NOD)    
          NOD_NORMAL_Z(NOD) = NOD_NORMAL_Z(NOD)/AREA_SUM(NOD)
          RN = SQRT(NOD_NORMAL_X(NOD)**2 + NOD_NORMAL_Y(NOD)**2 + NOD_NORMAL_Z(NOD)**2)
          NOD_NORMAL_X(NOD) = NOD_NORMAL_X(NOD)/RN   
          NOD_NORMAL_Y(NOD) = NOD_NORMAL_Y(NOD)/RN  
          NOD_NORMAL_Z(NOD) = NOD_NORMAL_Z(NOD)/RN

!          ewrite(3,*) 'Coordinate:',X(NOD),Y(NOD),Z(NOD)
!          ewrite(3,*) 'Node Normal:',NOD_NORMAL_X(NOD),NOD_NORMAL_Y(NOD),NOD_NORMAL_Z(NOD)
        ENDIF 
      ENDDO
      CALL flcomms_update(halo_tag, NOD_NORMAL_X, 1, 1, 0)
      CALL flcomms_update(halo_tag, NOD_NORMAL_Y, 1, 1, 0)
      CALL flcomms_update(halo_tag, NOD_NORMAL_Z, 1, 1, 0)
      ewrite(3,*) 'FOUND THIS MANY SURFACE ELEMENTS ON THE CYLINDER:',COUNTER
      ewrite(3,*)&
     &     'sum(ELE_AREA),sum(ELE_SURFACE_FLAG:',sum(ELE_AREA),sum(ELE_SURFACE_FLAG)
      ewrite(3,*) 'sum(ELE_SURFACE_FLAG):',sum(ELE_SURFACE_FLAG)
      ewrite(3,*) 'sum(idfield(1:NNODP)):',sum(idfield(1:NNODP))

  end subroutine on_body

      SUBROUTINE ON_CYLINDER(X,Y,Z,NORMALX,NORMALY,NORMALZ,ELE_SURFACE_FLAG,&
     & NONODS,TOTELE,STOTEL,NLOC,SNLOC,NDGLNO,SNDGLN,FLAG_NODES,GOT_FLAG_NODES)
      INTEGER NONODS,TOTELE,STOTEL,NLOC,SNLOC,ELE_SURFACE_FLAG(STOTEL)      
      REAL    X(NONODS),Y(NONODS),Z(NONODS)
      REAL    NORMALX(NONODS),NORMALY(NONODS),NORMALZ(NONODS), FLAG_NODES(NONODS)
      LOGICAL GOT_FLAG_NODES
      INTEGER NDGLNO(NLOC*TOTELE),SNDGLN(SNLOC*STOTEL)
      REAL    RADIUS,XC,YC,ZC,NORM,TOL,RADSQ
      INTEGER COUNT,COUNTER,GLOBI,SELE,SILOC
      
      IF(GOT_FLAG_NODES) FLAG_NODES(1:NONODS) = 0.0

      XC = 0.5
      YC = 0.2
      ZC = 0.0
      RADIUS = 0.05
      RADSQ = RADIUS**2
      TOL = 0.0001
      
      COUNTER = 0
      ELE_SURFACE_FLAG(1:STOTEL) = 0
      do SELE = 1,STOTEL
        COUNT = 0
      do SILOC = 1,SNLOC
          GLOBI = SNDGLN((SELE-1)*SNLOC + SILOC)
          IF( ABS(( (X(GLOBI)-XC)**2 + (Y(GLOBI)-YC)**2 )-RADSQ).LT.TOL ) COUNT = COUNT + 1
        END DO
       
        IF(COUNT.EQ.SNLOC) THEN 
          COUNTER = COUNTER+1         
          ELE_SURFACE_FLAG(SELE) = 1
      do SILOC = 1,SNLOC
            GLOBI = SNDGLN((SELE-1)*SNLOC + SILOC)
            IF(GOT_FLAG_NODES) FLAG_NODES(GLOBI) = 1.0
            NORMALX(GLOBI) = X(GLOBI) - XC 
            NORMALY(GLOBI) = Y(GLOBI) - YC
            NORMALZ(GLOBI) = 0.0
            NORM = SQRT(NORMALX(GLOBI)**2 + NORMALY(GLOBI)**2 + NORMALZ(GLOBI)**2)
            NORMALX(GLOBI) = NORMALX(GLOBI)/NORM
            NORMALY(GLOBI) = NORMALY(GLOBI)/NORM
            NORMALZ(GLOBI) = NORMALZ(GLOBI)/NORM
          END DO
        END IF

      END DO
      ewrite(3,*) 'FOUND THIS MANY SURFACE ELEMENTS ON THE CYLINDER:',COUNTER

  end subroutine on_cylinder

      SUBROUTINE ON_CYLINDER2D(X,Y,Z,NORMALX,NORMALY,NORMALZ,ELE_SURFACE_FLAG,&
     & NONODS,TOTELE,STOTEL,NLOC,SNLOC,NDGLNO,SNDGLN,FLAG_NODES,GOT_FLAG_NODES)
      INTEGER NONODS,TOTELE,STOTEL,NLOC,SNLOC,ELE_SURFACE_FLAG(STOTEL)      
      REAL    X(NONODS),Y(NONODS),Z(NONODS)
      REAL    NORMALX(NONODS),NORMALY(NONODS),NORMALZ(NONODS), FLAG_NODES(NONODS)
      LOGICAL GOT_FLAG_NODES
      INTEGER NDGLNO(NLOC*TOTELE),SNDGLN(SNLOC*STOTEL)
      REAL    RADIUS,XC,YC,ZC,NORM,TOL,RADSQ
      INTEGER COUNT,COUNTER,GLOBI,SELE,SILOC
      
      IF(GOT_FLAG_NODES) FLAG_NODES(1:NONODS) = 0.0

      XC = 0.2
      YC = 0.2
      ZC = 0.0
      RADIUS = 0.05
      RADSQ = RADIUS**2
      TOL = 0.0001
      
      COUNTER = 0
      ELE_SURFACE_FLAG(1:STOTEL) = 0
      do SELE = 1,STOTEL
        COUNT = 0
      do SILOC = 1,SNLOC
          GLOBI = SNDGLN((SELE-1)*SNLOC + SILOC)
          IF( ABS(( (X(GLOBI)-XC)**2 + (Y(GLOBI)-YC)**2 )-RADSQ).LT.TOL ) COUNT = COUNT + 1
        END DO
       
        IF(COUNT.EQ.SNLOC) THEN 
          COUNTER = COUNTER+1         
          ELE_SURFACE_FLAG(SELE) = 1
      do SILOC = 1,SNLOC
            GLOBI = SNDGLN((SELE-1)*SNLOC + SILOC)
            IF(GOT_FLAG_NODES) FLAG_NODES(GLOBI) = 1.0
            NORMALX(GLOBI) = X(GLOBI) - XC 
            NORMALY(GLOBI) = Y(GLOBI) - YC
            NORMALZ(GLOBI) = 0.0
            NORM = SQRT(NORMALX(GLOBI)**2 + NORMALY(GLOBI)**2 + NORMALZ(GLOBI)**2)
            NORMALX(GLOBI) = NORMALX(GLOBI)/NORM
            NORMALY(GLOBI) = NORMALY(GLOBI)/NORM
            NORMALZ(GLOBI) = NORMALZ(GLOBI)/NORM
          END DO
        END IF

      END DO
      ewrite(3,*) 'FOUND THIS MANY SURFACE ELEMENTS ON THE CYLINDER:',COUNTER

  end subroutine on_cylinder2d
  
end module diagnostic_body_drag_module

