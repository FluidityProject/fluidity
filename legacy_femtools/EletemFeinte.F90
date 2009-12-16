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

module EletemFeinte_module

   use FLDebug

   contains

!----------------------------------------------------------------------------------------
  
SUBROUTINE FEINTE(SHORT,LONG,MSF,MLF,DISCRA,GELONG,&
                 X,Y,Z,&
                 M,N,NLX,NLY,NLZ,WEIGHT,MLOC,NLOC,NGI,D3,DCYL, &
                 NDGLNO,VNDGLN,PNDGLN,&
                 DETWEI,RGI,&
                 TOTELE,FREDOP,NONODS,LNONOD)
!!< This sub interpolates a vector LONG from the basis funs N 
!!< onto the basis function M- which can be discontinuous. 
!!< If DISCRA then we feed back into the vector SHORT then 
!!< discretised version. 
!!< If DISCRA then we do,nt solve a matrix equation
!!< If GELONG then get a long vector from a short one.

IMPLICIT NONE

REAL PIE,TWOPIE
PARAMETER(PIE=3.141592654)
INTEGER TOTELE,FREDOP,NONODS,LNONOD,MLOC,NLOC,NGI
REAL M(MLOC,NGI)
REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
REAL WEIGHT(NGI)
REAL DETJ,DETWEI(NGI),RGI(NGI)
REAL SHORT(FREDOP),LONG(LNONOD)
REAL MSF(FREDOP),MLF(NONODS)
REAL X(NONODS),Y(NONODS),Z(NONODS)
! NB in cylindrical coords Y is the radius. 
LOGICAL D3,DCYL
LOGICAL DISCRA,GELONG
INTEGER NDGLNO(NLOC*TOTELE),VNDGLN(NLOC*TOTELE)
INTEGER PNDGLN(MLOC*TOTELE)

REAL VLM,VLMN,VLNM
REAL AGI,BGI,CGI,DGI,EGI,FGI,GGI,HGI,KGI
INTEGER L,ELE,GLOBI,GLOBJ,GI
INTEGER IGL,ILOC,JLOC,I

ewrite(3,*) 'inside feinte'
         
ewrite(3,*)'1 FEINTE TOTELE,FREDOP,NONODS,D3,DCYL',&
            TOTELE,FREDOP,NONODS,D3,DCYL
ewrite(3,*)'d3,dcyl,nloc,mloc,ngi:',d3,dcyl,nloc,mloc,ngi

 IF(.NOT.GELONG) THEN 
   SHORT(1:FREDOP) = 0.0
 ELSE       
   LONG(1:LNONOD) = 0.0         
 ENDIF
       
ewrite(3,*)'=DISCRA,GELONG,FREDOP:',DISCRA,GELONG,FREDOP
 
 IF(.NOT.DISCRA) THEN
   IF(.NOT.GELONG) THEN
      MSF(1:FREDOP) = 0.0
   ELSE
      MLF(1:NONODS) = 0.0
   ENDIF
 ENDIF

 do  GI=1,NGI
   RGI(GI)=1.0
 end do 
 
 TWOPIE=1.0
 IF(DCYL) TWOPIE=2.0*PIE
 
! Finite element interpolation(but lumped)
 do  ELE=1,TOTELE

   IF(.NOT.D3) THEN
   
      do  GI=1,NGI
         AGI=0.
         BGI=0.
         CGI=0.
         DGI=0.

         IF(DCYL) RGI(GI)=0.

         do  L=1,NLOC

            IGL=NDGLNO((ELE-1)*NLOC+L)

            AGI=AGI + NLX(L,GI)*X(IGL) 
            BGI=BGI + NLX(L,GI)*Y(IGL) 
            CGI=CGI + NLY(L,GI)*X(IGL) 
            DGI=DGI + NLY(L,GI)*Y(IGL) 
            IF(DCYL)   RGI(GI)=RGI(GI)+N(L,GI)*Y(IGL)
            
         end do 

         DETJ= AGI*DGI-BGI*CGI 

         DETWEI(GI)=TWOPIE*DETJ*WEIGHT(GI)

      end do 
      
   ELSE
   
      do  GI=1,NGI

         AGI=0.
         BGI=0.
         CGI=0.

         DGI=0.
         EGI=0.
         FGI=0.

         GGI=0.
         HGI=0.
         KGI=0.

         do  L=1,NLOC
         
            IGL=NDGLNO((ELE-1)*NLOC+L)
            AGI=AGI+NLX(L,GI)*X(IGL) 
            BGI=BGI+NLX(L,GI)*Y(IGL) 
            CGI=CGI+NLX(L,GI)*Z(IGL) 

            DGI=DGI+NLY(L,GI)*X(IGL) 
            EGI=EGI+NLY(L,GI)*Y(IGL) 
            FGI=FGI+NLY(L,GI)*Z(IGL) 

            GGI=GGI+NLZ(L,GI)*X(IGL) 
            HGI=HGI+NLZ(L,GI)*Y(IGL) 
            KGI=KGI+NLZ(L,GI)*Z(IGL) 

         end do 

         DETJ=AGI*(EGI*KGI-FGI*HGI) &
                  -BGI*(DGI*KGI-FGI*GGI) &
                  +CGI*(DGI*HGI-EGI*GGI)
                  
         DETWEI(GI)=DETJ*WEIGHT(GI)
         
      end do 
      
   ENDIF

   IF(.NOT.GELONG) THEN

      do  ILOC=1,MLOC
         
         GLOBI=PNDGLN( (ELE-1)*MLOC +ILOC)

         IF(.NOT.DISCRA) THEN
         
            do  JLOC=1,MLOC
            
               VLM=0.
               
               do  GI=1,NGI
               
                  VLM=VLM+ M(ILOC,GI)*M(JLOC,GI)*DETWEI(GI)*RGI(GI)
                  
               end do 
               
               MSF(GLOBI)=MSF(GLOBI)+VLM
               
            end do 
            
         ENDIF

         do  JLOC=1,NLOC
            
            GLOBJ=VNDGLN((ELE-1)*NLOC+JLOC)
            VLMN=0.
            
            do  GI=1,NGI
            
               VLMN=VLMN+M(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)*RGI(GI)
               
            end do 
            
            SHORT(GLOBI)=SHORT(GLOBI)+VLMN*LONG(GLOBJ)
            
         end do 

      end do 

   ELSE

      do  ILOC=1,NLOC
         
         GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)

         IF(.NOT.DISCRA) THEN
         
            do  JLOC=1,NLOC
            
               VLM=0.
               
               do  GI=1,NGI
               
                  VLM=VLM+ N(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)*RGI(GI)
                  
               end do
               
               MLF(GLOBI)=MLF(GLOBI)+VLM
               
            end do 
           
         ENDIF

         do  JLOC=1,MLOC
         
            GLOBJ=PNDGLN( (ELE-1)*MLOC +JLOC)
             
            VLNM=0.
            
            do  GI=1,NGI
            
               VLNM=VLNM+N(ILOC,GI)*M(JLOC,GI)*DETWEI(GI)*RGI(GI)
             
            end do 
            
            LONG(GLOBI)=LONG(GLOBI)+VLNM*SHORT(GLOBJ)

         end do 
         
      end do 
         
! ENDOF IF(GELONG) THEN ELSE..
   ENDIF

 end do ! TOTELE LOOP

! now solve for SHORT using lumped approximation. 
ewrite(3,*)'obtaining short DISCRA:',DISCRA
ewrite(3,*) 'about to solve'

 IF(.NOT.DISCRA) THEN
 
   IF(.NOT.GELONG) THEN
   
      do  I=1,FREDOP
      
         SHORT(I)=SHORT(I)/MSF(I)
         
      end do 
      
   ELSE
   
      do  I=1,NONODS
      
         LONG(I)=LONG(I)/MLF(I)
         
      end do 
      
   ENDIF
   
 ENDIF
 
ewrite(3,*) 'finishing feinte'
RETURN
END SUBROUTINE FEINTE


end module EletemFeinte_module
