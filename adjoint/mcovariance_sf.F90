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
SUBROUTINE MCOVARIANCE_SF(GCOVARI,F_SMOOTHNESS,NOCVA,NONODS,XNONOD,NBIGM,        &
                              STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI,                         &
                              FINDRM,NCOLM,COLM,                                         &
                              FSNDID,SNDGLN,X,Y,Z,H,                                     &
                              SN,SNLX,SNLY,SWEIGH,D3,                                    &
                              ITSTIME,MAXREC,MRECORD,IMEM_CVA,CENTRM,LINITS,GLOITS,DCYL)                      

      use FLDebug
      use AdvectionDiffusion
      use position_in_matrix
      IMPLICIT NONE
      LOGICAL  D3,DCYL
      INTEGER  ITSTIME,NOCVA,MAXREC,NCOLM,LINITS,GLOITS
      INTEGER  STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI
      INTEGER  NONODS,XNONOD,NBIGM
      REAL     F_SMOOTHNESS
      REAL     GCOVARI(NOCVA)
      INTEGER  FINDRM(NONODS+1),COLM(NCOLM)
      INTEGER  SNDGLN(STOTEL*SNLOC)
      REAL     FSNDID(NONODS)
      REAL     X(XNONOD),Y(XNONOD),Z(XNONOD),H(XNONOD)
      REAL     SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
      INTEGER  MRECORD(MAXREC),IMEM_CVA(NOCVA)
      INTEGER  CENTRM(NONODS)
! the tensor matrix
      REAL,DIMENSION(:),ALLOCATABLE  :: kmupxx,kmupxy,kmupxz,kmupyx,kmupyy,kmupyz,kmupzx,kmupzy,kmupzz
      REAL,DIMENSION(:),ALLOCATABLE  :: GRAVE,KSF,MMASS,MML
      REAL,DIMENSION(:,:),ALLOCATABLE  :: SNX,SNY,SNZ
      REAL     SWEIGH(SNGI)
! Working arrays/local variables
      REAL     AGI,BGI,CGI,DGI
      REAL     AIJ,RNN
      INTEGER  SELE,ILOC,JLOC,SILOC,SJLOC,L,IGL,IGLT,GLOBI,GLOBJ,GLOBIP,GLOBJP,GI
      INTEGER  HIT,NOD,NODI,NODJ,ITNOD,JTNOD,INUM,UPPER,LOWER,COUNT
      INTEGER  COUNTER
      REAL     DETWEI(NGI),DET
      LOGICAL  YES
! local variables
       REAL NORMXN(SNGI),NORMYN(SNGI),NORMZN(SNGI)
       REAL NORMX,NORMY,NORMZ
      REAL     kmupxxD(NGI),kmupxyD(NGI),kmupxzD(NGI), &
               kmupyxD(NGI),kmupyyD(NGI),kmupyzD(NGI),kmupzxD(NGI),kmupzyD(NGI),kmupzzD(NGI)
      INTEGER  ICENT
      INTEGER  i, j, k,k1,k2
      REAL SAREA

      ALLOCATE( KSF(NCOLM) )
      ALLOCATE( GRAVE(NONODS) )
      ALLOCATE( MMASS(NCOLM) )
      ALLOCATE( MML(NCOLM) )
      ALLOCATE( SNX(SNLOC,SNGI) )
      ALLOCATE( SNY(SNLOC,SNGI) )
      ALLOCATE( SNZ(SNLOC,SNGI) )
      ALLOCATE( kmupxx(NONODS) )
      ALLOCATE( kmupxy(NONODS) )
      ALLOCATE( kmupxz(NONODS) )
      ALLOCATE( kmupyx(NONODS) )
      ALLOCATE( kmupyy(NONODS) )
      ALLOCATE( kmupyz(NONODS) )
      ALLOCATE( kmupzx(NONODS) )
      ALLOCATE( kmupzy(NONODS) )
      ALLOCATE( kmupzz(NONODS) )
 
      KSF(1:NCOLM) = 0.0
      GRAVE(1:NONODS) = 0.0
      MMASS(1:NCOLM) = 0.0
      MML(1:NCOLM) = 0.0
      kmupxx(1:NONODS) = 0.0
      kmupxy(1:NONODS) = 0.0
      kmupxz(1:NONODS) = 0.0
      kmupyx(1:NONODS) = 0.0
      kmupyy(1:NONODS) = 0.0
      kmupyz(1:NONODS) = 0.0
      kmupzx(1:NONODS) = 0.0
      kmupzy(1:NONODS) = 0.0
      kmupzz(1:NONODS) = 0.0

      DO SILOC =1,SNLOC
      DO GI =1,SNGI
        SNX(SILOC,GI)=0.0
        SNY(SILOC,GI)=0.0
        SNZ(SILOC,GI)=0.0
      ENDDO
      ENDDO


     IF(ITSTIME.EQ.1) THEN
       GCOVARI(1:NOCVA) = 0.0
     ENDIF

! give the tesor along the inlet---in future, it should be written in a general way.

     DO NOD =1,NONODS        
      IF( ABS(X(NOD)-0.0).LT.1.0E-6 ) THEN
        kmupyy(NOD) = 100.0  !1.0
      ELSE
        kmupyy(NOD) = 0.0
      ENDIF

     ENDDO
      

! Counter counts how many nodes are found on the free surface
       COUNTER = 0         
   DO  SELE=1,STOTEL
! NB. These are the surface elements
         HIT=0
         DO SILOC=1,SNLOC
           IGL = SNDGLN((SELE-1)*SNLOC + SILOC)
           IF(FSNDID(IGL).GT.0.1) HIT=HIT+1      
!           ewrite(3,*) FSNDID(IGL),X(IGL),Y(IGL),Z(IGL) 
!           ewrite(3,*)  'IGL*****************************',IGL  
         ENDDO
!
!C If HIT=SNLOC then we are on the free surface so continue
     IF(HIT.EQ.SNLOC) THEN
!C SILOC=1 OK here as this is a surfac rather than volume element.

         if(x(igl).lt.1000.0) then
           ewrite(3,*) FSNDID(IGL),X(IGL),Y(IGL),Z(IGL) 
           ewrite(3,*)  'IGL*****************************',IGL  
         endif
           COUNTER = COUNTER +1

         DO GI = 1,SNGI 
            AGI = 0.
            BGI = 0.
            CGI = 0.
            DGI = 0.
              kmupxxD(GI) = 0.
              kmupxyD(GI) = 0.
              kmupxzD(GI) = 0.
              kmupyxD(GI) = 0.
              kmupyyD(GI) = 0.
              kmupyzD(GI) = 0.
              kmupzxD(GI) = 0.
              kmupzyD(GI) = 0.
              kmupzzD(GI) = 0.
            DO L = 1,SNLOC
              IGL = SNDGLN((SELE-1)*SNLOC+L)

              AGI = AGI + SNLX(L,GI)*X(IGL)
              BGI = BGI + SNLX(L,GI)*Y(IGL)
              CGI = CGI + SNLY(L,GI)*X(IGL)
              DGI = DGI + SNLY(L,GI)*Y(IGL)
              kmupxxD(GI) = kmupxxD(GI) + SN(L,GI)*kmupxx(IGL)
              kmupxyD(GI) = kmupxyD(GI) + SN(L,GI)*kmupxy(IGL)
              kmupxzD(GI) = kmupxzD(GI) + SN(L,GI)*kmupxz(IGL)
              kmupyxD(GI) = kmupyxD(GI) + SN(L,GI)*kmupyx(IGL)
              kmupyyD(GI) = kmupyyD(GI) + SN(L,GI)*kmupyy(IGL)
              kmupyzD(GI) = kmupyzD(GI) + SN(L,GI)*kmupyz(IGL)
              kmupzxD(GI) = kmupzxD(GI) + SN(L,GI)*kmupzx(IGL)
              kmupzyD(GI) = kmupzyD(GI) + SN(L,GI)*kmupzy(IGL)
              kmupzzD(GI) = kmupzzD(GI) + SN(L,GI)*kmupzz(IGL)
            ENDDO

!            DET = AGI*DGI-BGI*CGI
!            DETWEI(GI) = DET*SWEIGH(GI)

         ENDDO !end gi      

       CALL SDETNX2(SELE,SNDGLN,                                     &
             STOTEL,XNONOD,SNLOC,SNGI,                              &
             X,Y,Z,                                                 & 
             SN,SNLX,SNLY, SWEIGH, DETWEI,SAREA, D3,DCYL,           &
             NORMXN,NORMYN,NORMZN,                                  &
             NORMX,NORMY,NORMZ) 

   DO ILOC = 1,SNLOC
           GLOBI = SNDGLN((SELE-1)*SNLOC+ILOC)
       DO JLOC = 1,SNLOC
            GLOBJ = SNDGLN((SELE-1)*SNLOC+JLOC)

              IF(X(GLOBI).LT.1.0E-6.AND.X(GLOBJ).LT.1.0E-6) THEN              
                  RNN  =0.0
            DO GI = 1,SNGI
!                RNN  = RNN+SN(ILOC,GI)*SN(JLOC,GI)*DETWEI(GI)
                RNN  = RNN+SN(ILOC,GI)*SN(JLOC,GI)*DETWEI(GI)*kmupyyD(GI)
                  ENDDO !enddo gi

! the mass matric MMASS
!--------------------
!....          Find COUNT. 
                CALL POSINMAT(COUNT,GLOBI,GLOBJ,NONODS,FINDRM,COLM,NCOLM)

                MMASS(COUNT)= MMASS(COUNT)+RNN


          ICENT=CENTRM(GLOBI)
 ! Lumped MAss Matrix
!--------------------
                MML(ICENT)=MML(ICENT)+RNN
              ENDIF
          ENDDO !enddo jloc

         ENDDO !enddo iloc

     ENDIF
  ENDDO  !end sele

! the anisotropic matrix
   open(1,file='mml.dat')
          DO NOD=1,NONODS
             DO COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
               KSF(COUNT) = -(MMASS(COUNT)-MML(COUNT))
       write(1,*) NOD,COUNT,KSF(COUNT),MML(COUNT),MMASS(COUNT)
             END DO
          ENDDO 
  close(1)


!....    Spatial regularisation...

! Gradient.....
   open(1,file='grave.dat')
   open(2,file='grave1.dat')
            DO NOD=1,NONODS
               DO COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
!                  GRAVE(NOD)=GRAVE(NOD)+KSF(COUNT)*H(COLM(COUNT))        
                  GRAVE(NOD)=GRAVE(NOD)+KSF(COUNT)*H(COLM(COUNT))       
       write(1,*) NOD,grave(nod),COUNT,KSF(COUNT),COLM(COUNT),H(COLM(COUNT)),X(COLM(COUNT))

               END DO
       write(1,*) '*************'
     if(grave(nod).gt.1.0e-6) then
       write(2,*) NOD,grave(nod)
     endif
            END DO
   close(1)
   close(2)
 

! for forward model
           DO I= MRECORD(ITSTIME),MRECORD(ITSTIME+1)-1
             NOD=IMEM_CVA(I)
             GCOVARI(I)=GRAVE(NOD)
           ENDDO
! cost function
             IF(ITSTIME.EQ.1) THEN
               F_SMOOTHNESS=0.0
             ENDIF
             DO NOD=1,NONODS
               F_SMOOTHNESS=F_SMOOTHNESS + 0.5*H(NOD)*GRAVE(NOD)
             END DO
    if(.true.) then
      open(1,file='gcovari.dat')
       write(1,*) 'nocva=',nocva,F_SMOOTHNESS
       write(1,*) 'gcovari inside grad_BC_FS'
       write(1,*) (gcovari(i),i=MRECORD(ITSTIME),MRECORD(ITSTIME+1)-1)
     close(1)
    endif

      DEALLOCATE( KSF )
      DEALLOCATE( GRAVE )
      DEALLOCATE( MMASS )
      DEALLOCATE( MML )
      DEALLOCATE( SNX )
      DEALLOCATE( SNY )
      DEALLOCATE( SNZ )
      DEALLOCATE( kmupxx )
      DEALLOCATE( kmupxy )
      DEALLOCATE( kmupxz )
      DEALLOCATE( kmupyx )
      DEALLOCATE( kmupyy )
      DEALLOCATE( kmupyz )
      DEALLOCATE( kmupzx )
      DEALLOCATE( kmupzy )
      DEALLOCATE( kmupzz )

END SUBROUTINE MCOVARIANCE_SF






SUBROUTINE MCOVARIANCE_SF1(GRAVE,GCOVARI,F_SMOOTHNESS,NOCVA,NONODS,XNONOD,NBIGM,        &
                              STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI,                         &
                              FINDRM,NCOLM,COLM,                                         &
                              FSNDID,SNDGLN,X,Y,Z,H,                                     &
                              SN,SNLX,SNLY,SWEIGH,D3,                                    &
                              ITSTIME,MAXREC,MRECORD,IMEM_CVA,CENTRM,LINITS,GLOITS,DCYL)                      

     use FLDebug
      use AdvectionDiffusion
      use position_in_matrix
      IMPLICIT NONE
      LOGICAL  D3,DCYL
      INTEGER  ITSTIME,NOCVA,MAXREC,NCOLM,LINITS,GLOITS
      INTEGER  STOTEL,TOTELE,NLOC,SNLOC,NGI,SNGI
      INTEGER  NONODS,XNONOD,NBIGM
      REAL     F_SMOOTHNESS
      REAL     GCOVARI(NOCVA),GRAVE(NONODS)
      INTEGER  FINDRM(NONODS+1),COLM(NCOLM)
      INTEGER  SNDGLN(STOTEL*SNLOC)
      REAL     FSNDID(NONODS)
      REAL     X(XNONOD),Y(XNONOD),Z(XNONOD),H(XNONOD)
      REAL     SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
      INTEGER  MRECORD(MAXREC),IMEM_CVA(NOCVA)
      INTEGER  CENTRM(NONODS)
! the tensor matrix
      REAL,DIMENSION(:),ALLOCATABLE  :: kmupxx,kmupxy,kmupxz,kmupyx,kmupyy,kmupyz,kmupzx,kmupzy,kmupzz
      REAL,DIMENSION(:),ALLOCATABLE  :: KSF,MMASS,MML
      REAL,DIMENSION(:,:),ALLOCATABLE  :: SNX,SNY,SNZ
      REAL     SWEIGH(SNGI)
! Working arrays/local variables
      REAL     AGI,BGI,CGI,DGI
      REAL     AIJ,RNN
      INTEGER  SELE,ILOC,JLOC,SILOC,SJLOC,L,IGL,IGLT,GLOBI,GLOBJ,GLOBIP,GLOBJP,GI
      INTEGER  HIT,NOD,NODI,NODJ,ITNOD,JTNOD,INUM,UPPER,LOWER,COUNT
      INTEGER  COUNTER
      REAL     DETWEI(NGI),DET
      LOGICAL  YES
! local variables
       REAL NORMXN(SNGI),NORMYN(SNGI),NORMZN(SNGI)
       REAL NORMX,NORMY,NORMZ
      REAL     kmupxxD(NGI),kmupxyD(NGI),kmupxzD(NGI), &
               kmupyxD(NGI),kmupyyD(NGI),kmupyzD(NGI),kmupzxD(NGI),kmupzyD(NGI),kmupzzD(NGI)
      INTEGER  ICENT
      INTEGER  i, j, k,k1,k2
      REAL SAREA

      ALLOCATE( KSF(NCOLM) )
      ALLOCATE( MMASS(NCOLM) )
      ALLOCATE( MML(NCOLM) )
      ALLOCATE( SNX(SNLOC,SNGI) )
      ALLOCATE( SNY(SNLOC,SNGI) )
      ALLOCATE( SNZ(SNLOC,SNGI) )
      ALLOCATE( kmupxx(NONODS) )
      ALLOCATE( kmupxy(NONODS) )
      ALLOCATE( kmupxz(NONODS) )
      ALLOCATE( kmupyx(NONODS) )
      ALLOCATE( kmupyy(NONODS) )
      ALLOCATE( kmupyz(NONODS) )
      ALLOCATE( kmupzx(NONODS) )
      ALLOCATE( kmupzy(NONODS) )
      ALLOCATE( kmupzz(NONODS) )
 
      KSF(1:NCOLM) = 0.0
      MMASS(1:NCOLM) = 0.0
      MML(1:NCOLM) = 0.0
      kmupxx(1:NONODS) = 0.0
      kmupxy(1:NONODS) = 0.0
      kmupxz(1:NONODS) = 0.0
      kmupyx(1:NONODS) = 0.0
      kmupyy(1:NONODS) = 0.0
      kmupyz(1:NONODS) = 0.0
      kmupzx(1:NONODS) = 0.0
      kmupzy(1:NONODS) = 0.0
      kmupzz(1:NONODS) = 0.0

      DO SILOC =1,SNLOC
      DO GI =1,SNGI
        SNX(SILOC,GI)=0.0
        SNY(SILOC,GI)=0.0
        SNZ(SILOC,GI)=0.0
      ENDDO
      ENDDO


     IF(ITSTIME.EQ.1) THEN
       GCOVARI(1:NOCVA) = 0.0
     ENDIF

! give the tesor along the inlet---in future, it should be written in a general way.

     DO NOD =1,NONODS        
      IF( ABS(X(NOD)-0.0).LT.1.0E-6 ) THEN
        kmupyy(NOD) = 100.0 !100
      ELSE
        kmupyy(NOD) = 0.0
      ENDIF

     ENDDO
      

! Counter counts how many nodes are found on the free surface
       COUNTER = 0         
   DO  SELE=1,STOTEL
! NB. These are the surface elements
         HIT=0
         DO SILOC=1,SNLOC
           IGL = SNDGLN((SELE-1)*SNLOC + SILOC)
           IF(FSNDID(IGL).GT.0.1) HIT=HIT+1      
!           ewrite(3,*) FSNDID(IGL),X(IGL),Y(IGL),Z(IGL) 
!           ewrite(3,*)  'IGL*****************************',IGL  
         ENDDO
!
!C If HIT=SNLOC then we are on the free surface so continue
     IF(HIT.EQ.SNLOC) THEN
!C SILOC=1 OK here as this is a surfac rather than volume element.

         if(x(igl).lt.1000.0) then
           ewrite(3,*) FSNDID(IGL),X(IGL),Y(IGL),Z(IGL) 
           ewrite(3,*)  'IGL*****************************',IGL  
         endif
           COUNTER = COUNTER +1

         DO GI = 1,SNGI 
            AGI = 0.
            BGI = 0.
            CGI = 0.
            DGI = 0.
              kmupxxD(GI) = 0.
              kmupxyD(GI) = 0.
              kmupxzD(GI) = 0.
              kmupyxD(GI) = 0.
              kmupyyD(GI) = 0.
              kmupyzD(GI) = 0.
              kmupzxD(GI) = 0.
              kmupzyD(GI) = 0.
              kmupzzD(GI) = 0.
            DO L = 1,SNLOC
              IGL = SNDGLN((SELE-1)*SNLOC+L)

              AGI = AGI + SNLX(L,GI)*X(IGL)
              BGI = BGI + SNLX(L,GI)*Y(IGL)
              CGI = CGI + SNLY(L,GI)*X(IGL)
              DGI = DGI + SNLY(L,GI)*Y(IGL)
              kmupxxD(GI) = kmupxxD(GI) + SN(L,GI)*kmupxx(IGL)
              kmupxyD(GI) = kmupxyD(GI) + SN(L,GI)*kmupxy(IGL)
              kmupxzD(GI) = kmupxzD(GI) + SN(L,GI)*kmupxz(IGL)
              kmupyxD(GI) = kmupyxD(GI) + SN(L,GI)*kmupyx(IGL)
              kmupyyD(GI) = kmupyyD(GI) + SN(L,GI)*kmupyy(IGL)
              kmupyzD(GI) = kmupyzD(GI) + SN(L,GI)*kmupyz(IGL)
              kmupzxD(GI) = kmupzxD(GI) + SN(L,GI)*kmupzx(IGL)
              kmupzyD(GI) = kmupzyD(GI) + SN(L,GI)*kmupzy(IGL)
              kmupzzD(GI) = kmupzzD(GI) + SN(L,GI)*kmupzz(IGL)
            ENDDO

!            DET = AGI*DGI-BGI*CGI
!            DETWEI(GI) = DET*SWEIGH(GI)

         ENDDO !end gi      

       CALL SDETNX2(SELE,SNDGLN,                                     &
             STOTEL,XNONOD,SNLOC,SNGI,                              &
             X,Y,Z,                                                 & 
             SN,SNLX,SNLY, SWEIGH, DETWEI,SAREA, D3,DCYL,           &
             NORMXN,NORMYN,NORMZN,                                  &
             NORMX,NORMY,NORMZ) 

   DO ILOC = 1,SNLOC
           GLOBI = SNDGLN((SELE-1)*SNLOC+ILOC)
       DO JLOC = 1,SNLOC
            GLOBJ = SNDGLN((SELE-1)*SNLOC+JLOC)

              IF(X(GLOBI).LT.1.0E-6.AND.X(GLOBJ).LT.1.0E-6) THEN              
                  RNN  =0.0
            DO GI = 1,SNGI
!                RNN  = RNN+SN(ILOC,GI)*SN(JLOC,GI)*DETWEI(GI)
                RNN  = RNN+SN(ILOC,GI)*SN(JLOC,GI)*DETWEI(GI)*kmupyyD(GI)
                  ENDDO !enddo gi

! the mass matric MMASS
!--------------------
!....          Find COUNT. 
                CALL POSINMAT(COUNT,GLOBI,GLOBJ,NONODS,FINDRM,COLM,NCOLM)

                MMASS(COUNT)= MMASS(COUNT)+RNN


          ICENT=CENTRM(GLOBI)
 ! Lumped MAss Matrix
!--------------------
                MML(ICENT)=MML(ICENT)+RNN
              ENDIF
          ENDDO !enddo jloc

         ENDDO !enddo iloc

     ENDIF
  ENDDO  !end sele

! the anisotropic matrix
   open(1,file='mml.dat')
          DO NOD=1,NONODS
             DO COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
               KSF(COUNT) = -(MMASS(COUNT)-MML(COUNT))
       write(1,*) NOD,COUNT,KSF(COUNT),MML(COUNT),MMASS(COUNT)
             END DO
          ENDDO 
  close(1)


!....    Spatial regularisation...

! Gradient.....
   open(1,file='grave.dat')
   open(2,file='grave1.dat')
            DO NOD=1,NONODS
               DO COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
!                  GRAVE(NOD)=GRAVE(NOD)+KSF(COUNT)*H(COLM(COUNT))        
                  GRAVE(NOD)=GRAVE(NOD)+KSF(COUNT)*H(COLM(COUNT))       
       write(1,*) NOD,grave(nod),COUNT,KSF(COUNT),COLM(COUNT),H(COLM(COUNT)),X(COLM(COUNT))
               END DO
       write(1,*) '*************'
     if(grave(nod).gt.1.0e-6) then
       write(2,*) NOD,grave(nod)
     endif
            END DO
   close(1)
   close(2)

! cost function
             IF(ITSTIME.EQ.1) THEN
               F_SMOOTHNESS=0.0
             ENDIF
             DO NOD=1,NONODS
               F_SMOOTHNESS=F_SMOOTHNESS + 0.5*H(NOD)*GRAVE(NOD)
             END DO
 

      DEALLOCATE( KSF )
      DEALLOCATE( MMASS )
      DEALLOCATE( MML )
      DEALLOCATE( SNX )
      DEALLOCATE( SNY )
      DEALLOCATE( SNZ )
      DEALLOCATE( kmupxx )
      DEALLOCATE( kmupxy )
      DEALLOCATE( kmupxz )
      DEALLOCATE( kmupyx )
      DEALLOCATE( kmupyy )
      DEALLOCATE( kmupyz )
      DEALLOCATE( kmupzx )
      DEALLOCATE( kmupzy )
      DEALLOCATE( kmupzz )

END SUBROUTINE MCOVARIANCE_SF1
