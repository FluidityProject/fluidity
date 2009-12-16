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

#include "confdefs.h"
#include "fdebug.h"

module verticaldg_module
  use FLDebug
  use mesh_connections
  use shape_transformations
  use coordinates
  use Transform_elements
  use position_in_matrix
  use tr2d_module
#ifdef ADJOINT
  use assnav_adjoint_module
  use diff3d_adjoint_module
#else
  use assnav_module
  use diff3d_module
#endif
  use coriolis_module
  use spud

  implicit none

  private
  
  public :: verticaldg, verticaldgstore, longalltopsufdg, alltopsufdg,  &
       ctdgsup, dgsdetnxloc, dgsdetnxloc2, dgsimplnorm

contains

  SUBROUTINE VERTICALDG(VECX,VECY,VECZ, &
       NSUBVLOC,VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD, &
       DGP,ELEORD,GETELEORD,NLEVEL, &
       GETFREES,FREES,EQETID,GRAVTY,  &
       PNORMX,PNORMY,PNORMZ, &
       PGRAVX,PGRAVY,PGRAVZ,NDENST, &
       POTGRA,USEPOTGRA,  &
       PGNODS, &
       X,Y,Z, U,V,W, &
       INUSUB,NSUBNVLOC,NUSUB,NVSUB,NWSUB, &
       NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB,COGRAX, &
       M,MLX,MLY,MLZ,&
       N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC, &
       TOTELE, &
       GEOBAL,&
       NDGLNO,PGNDGLNO, XONDGL, &
       NONODS,XNONOD, &
       scfacth0,ISPHERE, &
       FINDRM,COLM,NCOLM,&
       IDGBCMOM,IDGBOYONLY,PSMOOTH,COLELE4,PHSOLVE,wdgravity)
    ! This subroutine solves the DG eqns for balanced pressure DGP.
    ! Number of source nodes, advection nodes and node positions.
    ! Number of elements, velocity nodes, ???
    ! If GEOBAL=-23 then have quadratic variation in the vertical with a DG method.
    ! If GEOBAL=-24 then have quadratic variation in the vertical with a DG method but
    ! use no b.c's when putting into the momentum eqns.
    ! If GEOBAL=-25 same as -24 but do not solve for verticaldg pressure 
    ! but use a quadratic pressure i.e. P1_(SGS or DG)-P2. 
    ! If GEOBAL=-26 same as -24 but do not solve for verticaldg pressure 
    ! & solve for p and free surface together.
    ! If GEOBAL=-27 same as -23 but solve for p and free surface together.
    ! If PHSOLVE then solve the vertical DG eqns else use the previous 
    ! hydrostic pressure value. 
    integer, intent(in) ::  totele,nonods,xnonod,pgnods,NSUBVLOC
    integer, intent(in):: ngi,nloc,mloc
    real, intent(inout):: vecx(nonods),vecy(nonods),vecz(nonods)
    integer, intent(in):: nsubnvloc
    real, intent(inout):: VECX_SGSADD(TOTELE*NSUBVLOC),VECY_SGSADD(TOTELE*NSUBVLOC)
    real, intent(inout):: VECZ_SGSADD(TOTELE*NSUBVLOC)
    real, intent(inout):: dgp(totele*mloc)

    real, parameter:: tolercent=0.00001
    INTEGER DGSMOOTH
    PARAMETER(DGSMOOTH=1)
    ! DGSMOOTH=1 for smoothing across elements.
    ! DGSMOOTH=1 for smoothing within element.

    integer, intent(inout):: eleord(totele) ! in if geteleord==.false.
    logical, intent(in):: geteleord
    integer, intent(in):: nlevel
    real, intent(in):: scfacth0
    real, intent(in):: u(nonods),v(nonods),w(nonods)
    integer, intent(in):: inusub
    real, intent(in):: nusub(totele*nsubnvloc),nvsub(totele*nsubnvloc)
    real, intent(in):: nwsub(totele*nsubnvloc)
  INTEGER, intent(in):: NSOGRASUB
  REAL, intent(in):: SOXGRASUB(TOTELE*NSOGRASUB),SOYGRASUB(TOTELE*NSOGRASUB)
  REAL, intent(in):: SOZGRASUB(TOTELE*NSOGRASUB)
     LOGICAL, intent(in):: COGRAX
    real, intent(in):: x(xnonod),y(xnonod),z(xnonod)
    logical, intent(in):: getfrees
    real, intent(in):: frees(pgnods),eqetid(pgnods)
    ! global node number maps.
    integer, intent(in):: ndglno(totele*nloc),pgndglno(totele*mloc)
    integer, intent(in):: xondgl(totele*nloc)
    integer, intent(in):: findrm(nonods+1),ncolm,colm(ncolm)

    integer, intent(in):: usepotgra,isphere,idgbcmom,idgboyonly
    real, intent(in):: pnormx(pgnods),pnormy(pgnods),pnormz(pgnods)
    real, intent(in):: pgravx(pgnods),pgravy(pgnods),pgravz(pgnods), ndenst(nonods)
    real, intent(in):: potgra(pgnods*usepotgra)

    real, intent(in):: n(nloc,ngi),nlx(nloc,ngi),nly(nloc,ngi),nlz(nloc,ngi)
    real, intent(in):: weight(ngi)

    real, intent(in):: m(mloc,ngi),mlx(mloc,ngi),mly(mloc,ngi),mlz(mloc,ngi)
    integer, intent(in):: geobal
    real, intent(in):: gravty
    Real, intent(in):: wdgravity(totele)

    ! local variables...
    integer smloc,sngi,snloc
    parameter(snloc=3,smloc=6,sngi=7)
    real sn(snloc,sngi),snlx(snloc,sngi),snly(snloc,sngi),snlz(snloc,sngi)
    real sweigh(sngi)
    real sm(smloc,sngi),smlx(smloc,sngi),smly(smloc,sngi),smlz(smloc,sngi)
    real sl1(sngi), sl2(sngi), sl3(sngi), sl4(sngi)
    real ud(ngi),vd(ngi),wd(ngi)
    real xd(ngi),yd(ngi),zd(ngi)
    real detwei(ngi)
    real nx(nloc,ngi),ny(nloc,ngi),nz(nloc,ngi)
    real mx(mloc,ngi),my(mloc,ngi),mz(mloc,ngi)
    real sdetwe(sngi)
    real fsxgi(ngi),fsygi(ngi),fszgi(ngi)
    real fsxgic(ngi),fsygic(ngi),fszgic(ngi)
    real gipgravx(ngi),gipgravy(ngi),gipgravz(ngi)
    real a11(ngi),a12(ngi),a13(ngi)
    real a21(ngi),a22(ngi),a23(ngi)
    real a31(ngi),a32(ngi),a33(ngi)
    real dirx(ngi),diry(ngi),dirz(ngi)
    real ddxdx(ngi),ddydy(ngi),ddzdz(ngi)
    real ome2(ngi)
    real souxgi(ngi),souygi(ngi),souzgi(ngi)
    real rx(ngi),ry(ngi),rz(ngi)
    real rxc(ngi),ryc(ngi),rzc(ngi)

    real snormxn(sngi),snormyn(sngi),snormzn(sngi)
    real sdirx(sngi),sdiry(sngi),sdirz(sngi)
    real ndotq(sngi),income(sngi)

    real srtvx(sngi),srtvy(sngi),srtvz(sngi)
    real srestx(sngi),sresty(sngi),srestz(sngi)

    real rtrest11(ngi),rtrest12(ngi),rtrest13(ngi)
    real rtrest21(ngi),rtrest22(ngi),rtrest23(ngi)
    real rtrest31(ngi),rtrest32(ngi),rtrest33(ngi)
    real rxr(ngi),ryr(ngi),rzr(ngi)
    real rxrc(ngi),ryrc(ngi),rzrc(ngi)

    integer siloc2iloc(snloc),othiloc(nloc)
    integer mlocsiloc2iloc(smloc),mlocothiloc(mloc)

    integer ele,iloc,jloc
    integer gi,sgi,globi,globj,globj2
    integer count2
    logical topsuf,mombc

    logical psmooth,fsmooth
    integer colele4(totele*5)
    logical phsolve

    real normx,normy,normz
    real sarea
    integer ele2
    integer jloc2,igl, max_its,adjits
    integer siloc,sjloc,IDGNOD
    real rnn
    real rn,invlmx,invlmy,invlmz,bothvlmx,bothvlmy,bothvlmz,poth
    real rrvecx,rrvecy,rrvecz
    REAL RRSX,RRSY,RRSZ,RNRRS,GRA_SIG_CHAN
    REAL RDIRX,RDIRY,RDIRZ
    integer iface,nface
    parameter(nface=4)
    integer loclist(nface,3)
  
    INTEGER one_vdg_ser_it
    SAVE one_vdg_ser_it
    DATA one_vdg_ser_it /0/

    real, allocatable, dimension(:)::rhsvec
    real, allocatable, dimension(:)::vecx2
    real, allocatable, dimension(:)::vecy2
    real, allocatable, dimension(:)::vecz2

    real, allocatable, dimension(:)::vecxf
    real, allocatable, dimension(:)::vecyf
    real, allocatable, dimension(:)::veczf
    real, allocatable, dimension(:)::ml
    real, allocatable, dimension(:)::VECX_SGSADD2
    real, allocatable, dimension(:)::VECY_SGSADD2
    real, allocatable, dimension(:)::VECZ_SGSADD2
    real, allocatable, dimension(:)::MLELE

    ewrite(1,*) 'in VERTICALDG psmooth=',psmooth

    ALLOCATE(RHSVEC(TOTELE*MLOC))

    ALLOCATE(VECX2(NONODS))
    ALLOCATE(VECY2(NONODS))
    ALLOCATE(VECZ2(NONODS))
    ALLOCATE(VECXF(NONODS))
    ALLOCATE(VECYF(NONODS))
    ALLOCATE(VECZF(NONODS))
    ALLOCATE(ML(NONODS))
    ALLOCATE(VECX_SGSADD2(NSUBVLOC*TOTELE))
    ALLOCATE(VECY_SGSADD2(NSUBVLOC*TOTELE))
    ALLOCATE(VECZ_SGSADD2(NSUBVLOC*TOTELE))
    ALLOCATE(MLELE(NSUBVLOC*TOTELE))

    VECX2(1:NONODS) = 0.0
    VECY2(1:NONODS) = 0.0
    VECZ2(1:NONODS) = 0.0
    VECXF(1:NONODS) = 0.0
    VECYF(1:NONODS) = 0.0
    VECZF(1:NONODS) = 0.0
    ML(1:NONODS) = 0.0
    VECX_SGSADD2 = 0.0
    VECY_SGSADD2 = 0.0
    VECZ_SGSADD2 = 0.0
    MLELE = 0.0

    ! Include the DG b.c's in the momentum eqns?...
    IF(GEOBAL.LT.-30) THEN
       MOMBC=(IDGBCMOM.EQ.1)
    ELSE
       MOMBC=.TRUE.
       IF(GEOBAL.EQ.-24) MOMBC=.FALSE.
       IF(GEOBAL.EQ.-25) MOMBC=.FALSE.
       IF(GEOBAL.EQ.-26) MOMBC=.FALSE.
    ENDIF
    ! Smooth out the resulting pressure force...
    !     PSMOOTH=.TRUE.
    ! If FSMOOTH smooth out the free surface balance and Coriolis
    FSMOOTH=.FALSE.

    ! Calculate surface shape functions SN, SM etc...
    CALL TRIQUA(SL1, SL2, SL3, SL3, SWEIGH, .FALSE.,SNGI)
    CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
         SNLOC,SNGI, &
         SN,SNLX,SNLY,SNLZ)
    CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
         SMLOC,SNGI, &
         SM,SMLX,SMLY,SMLZ)

    ! The local cords are in anti-clockwise order
    IFACE=1
    LOCLIST(IFACE,1)=1
    LOCLIST(IFACE,2)=2
    LOCLIST(IFACE,3)=3
    IFACE=2
    LOCLIST(IFACE,1)=1
    LOCLIST(IFACE,2)=2
    LOCLIST(IFACE,3)=4
    IFACE=3
    LOCLIST(IFACE,1)=3
    LOCLIST(IFACE,2)=2
    LOCLIST(IFACE,3)=4
    IFACE=4
    LOCLIST(IFACE,1)=1
    LOCLIST(IFACE,2)=3
    LOCLIST(IFACE,3)=4

    ! CALCULATE ELEMENT TO ELEMENT LIST**************
    ! and put in order of faces

    ewrite(1,*) 'going into 1ST ELEMENT LOOP PHSOLVE=',PHSOLVE

    IF(PHSOLVE) THEN

       IF(GETELEORD) THEN

          ! *********FIND ELEMENT ORDERING******************************
          ewrite(2,*) 'REORDERING THE ELEMENTS'
          call element_order(eleord,nlevel,totele,snloc,nloc,smloc,mloc,sngi, &
               xondgl,ndglno, &
               colele4,loclist,pgnods,nonods,xnonod,x,y,z, &
               sn,snlx,snly, sweigh, &
               isphere,sm,smlx,smly, &
               findrm,colm,ncolm, &
               pgndglno,pnormx,pnormy,pnormz,&
               one_vdg_ser_it,adjits)
       ENDIF

       ! Calculate part of rhs vec RHSVEC********************.
       call verticaldg_rhs(rhsvec, &
            ndglno, pgndglno, xondgl, &
            m, mlx, mly, mlz, &
            n, nlx, nly, nlz, &
            weight, &
            GEOBAL,&
            x, y, z, &
            u, v, w, &
            inusub,nsubnvloc,nusub,nvsub,nwsub, &
            NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB,COGRAX, &
            pgravx, pgravy, pgravz, &
            pnormx, pnormy, pnormz, &
            potgra, ndenst, frees, eqetid, &
            gravty, &
            isphere, scfacth0, idgboyonly, usepotgra, getfrees,totele, wdgravity)

       ewrite(2,*) 'r2norm(rhsvec,mloc*totele,0):',r2norm(rhsvec,mloc*totele,0)

       if((nlevel>1).or.(one_vdg_ser_it.eq.1)) then
          max_its=1
       else
          max_its=10
       end if

       call verticaldg_solve(dgp, rhsvec, &
            tolercent, max_its, &
            eleord, colele4, &
            ndglno, pgndglno, xondgl, &
            m, mlx, mly, mlz, &
            n, nlx, nly, nlz, &
            weight, &
            x, y, z, &
            pnormx, pnormy, pnormz, &
            isphere, &
            sn, snlx, snly, &
            sm, smlx, smly, &
            sweigh)

    end if

    ewrite(2,*) 'passed ok'
    CALL PMINMX(DGP,TOTELE*MLOC,'******DGP  ')

    ewrite(2,*) 'r2norm(dgp,mloc*totele,0):',r2norm(dgp,mloc*totele,0)

    ewrite(1,*) '******* 2-NOW PUT SOLUTION PSI into VECX,VECY,VECZ'

    element_loop: DO ELE=1,TOTELE

       ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
       CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MLOC,NGI, &
            N,NLX,NLY,NLZ, M,MLX,MLY,MLZ, WEIGHT, DETWEI,.TRUE.,.FALSE., &
            NX,NY,NZ, MX,MY,MZ, &
            A11,A12,A13, A21,A22,A23, A31,A32,A33, &
            XD,YD,ZD, &
            ISPHERE)

       DO GI=1,NGI
          OME2(GI)=2.*FUNOME(XD(GI),YD(GI),ZD(GI))
          IF((GEOBAL.EQ.-12).or.(GEOBAL.EQ.-22)) OME2(GI)=0.
          IF(IDGBOYONLY.EQ.1) OME2(GI)=0.

          GIPGRAVX(GI)=0.0
          GIPGRAVY(GI)=0.0
          GIPGRAVZ(GI)=0.0
          DIRX(GI)=0.0
          DIRY(GI)=0.0
          DIRZ(GI)=0.0
          DDXDX(GI)=0.0
          DDYDY(GI)=0.0
          DDZDZ(GI)=0.0
          FSXGI(GI)=0.0
          FSYGI(GI)=0.0
          FSZGI(GI)=0.0
          FSXGIC(GI)=0.0
          FSYGIC(GI)=0.0
          FSZGIC(GI)=0.0
          SOUXGI(GI)=0.0
          SOUYGI(GI)=0.0
          SOUZGI(GI)=0.0
          UD(GI)=0.0
          VD(GI)=0.0
          WD(GI)=0.0
          RX(GI)=0.0
          RY(GI)=0.0
          RZ(GI)=0.0
          RXC(GI)=0.0
          RYC(GI)=0.0
          RZC(GI)=0.0

          DO ILOC=1,MLOC
             IGL=PGNDGLNO((ELE-1)*MLOC+ILOC)
             IF(USEPOTGRA.EQ.1) THEN
                GIPGRAVX(GI)=GIPGRAVX(GI) + MX(ILOC,GI)*POTGRA(IGL)
                GIPGRAVY(GI)=GIPGRAVY(GI) + MY(ILOC,GI)*POTGRA(IGL)
                GIPGRAVZ(GI)=GIPGRAVZ(GI) + MZ(ILOC,GI)*POTGRA(IGL)
             ELSE
                GIPGRAVX(GI)=GIPGRAVX(GI) + M(ILOC,GI)*PGRAVX(IGL)
                GIPGRAVY(GI)=GIPGRAVY(GI) + M(ILOC,GI)*PGRAVY(IGL)
                GIPGRAVZ(GI)=GIPGRAVZ(GI) + M(ILOC,GI)*PGRAVZ(IGL)
             ENDIF
             DIRX(GI)=DIRX(GI) - M(ILOC,GI)*PNORMX(IGL)
             DIRY(GI)=DIRY(GI) - M(ILOC,GI)*PNORMY(IGL)
             DIRZ(GI)=DIRZ(GI) - M(ILOC,GI)*PNORMZ(IGL)
             DDXDX(GI)=DDXDX(GI)-MX(ILOC,GI)*PNORMX(IGL)
             DDYDY(GI)=DDYDY(GI)-MY(ILOC,GI)*PNORMY(IGL)
             DDZDZ(GI)=DDZDZ(GI)-MZ(ILOC,GI)*PNORMZ(IGL)

             IF(GETFREES) THEN
                FSXGI(GI)=FSXGI(GI) + MX(ILOC,GI)*(FREES(IGL)-EQETID(IGL))*GRAVTY*wdgravity(ele)
                FSYGI(GI)=FSYGI(GI) + MY(ILOC,GI)*(FREES(IGL)-EQETID(IGL))*GRAVTY*wdgravity(ele)
                FSZGI(GI)=FSZGI(GI) + MZ(ILOC,GI)*(FREES(IGL)-EQETID(IGL))*GRAVTY*wdgravity(ele)
             ENDIF
          END DO

          RN=SQRT(DIRX(gi)**2+DIRY(gi)**2+DIRZ(GI)**2)
          DIRX(GI)=DIRX(GI)/RN
          DIRY(GI)=DIRY(GI)/RN
          DIRZ(GI)=DIRZ(GI)/RN
          !  The rotation matrix in 3-D is R=
          !   T1X   T1Y   T1Z
          !   T2X   T2Y   T2Z
          !   NORMX NORMY NORMZ
          ! Calculate T1X,T1Y,T1Z, T2X,T2Y,T2Z from NORMX,NORMY,NORMZ
         if(MOMBC) then

             RTREST11(GI)=1.0 - DIRX(GI)*DIRX(GI)
             RTREST12(GI)=    - DIRX(GI)*DIRY(GI)
             RTREST13(GI)=    - DIRX(GI)*DIRZ(GI)

             RTREST21(GI)=    - DIRY(GI)*DIRX(GI)
             RTREST22(GI)=1.0 - DIRY(GI)*DIRY(GI)
             RTREST23(GI)=    - DIRY(GI)*DIRZ(GI)

             RTREST31(GI)=    - DIRZ(GI)*DIRX(GI)
             RTREST32(GI)=    - DIRZ(GI)*DIRY(GI)
             RTREST33(GI)=1.0 - DIRZ(GI)*DIRZ(GI)
          else
             RTREST11(GI)=1.0
             RTREST12(GI)=0.0
             RTREST13(GI)=0.0

             RTREST21(GI)=0.0
             RTREST22(GI)=1.0
             RTREST23(GI)=0.0

             RTREST31(GI)=0.0
             RTREST32(GI)=0.0
             RTREST33(GI)=1.0
          endif

          DO ILOC=1,NLOC
             IGL=NDGLNO((ELE-1)*NLOC+ILOC)
             SOUXGI(GI)=SOUXGI(GI) + N(ILOC,GI)*GIPGRAVX(GI)*NDENST(IGL)
             SOUYGI(GI)=SOUYGI(GI) + N(ILOC,GI)*GIPGRAVY(GI)*NDENST(IGL)
             SOUZGI(GI)=SOUZGI(GI) + N(ILOC,GI)*GIPGRAVZ(GI)*NDENST(IGL)
             UD(GI)=UD(GI) + N(ILOC,GI)*U(IGL)
             VD(GI)=VD(GI) + N(ILOC,GI)*V(IGL)
             WD(GI)=WD(GI) + N(ILOC,GI)*W(IGL)
             IF((INUSUB.NE.0).AND.(NSUBNVLOC.NE.0)) THEN
                IDGNOD=(ELE-1)*NLOC+ILOC
                UD(GI)=UD(GI) + N(ILOC,GI)*NUSUB(IDGNOD)
                VD(GI)=VD(GI) + N(ILOC,GI)*NVSUB(IDGNOD)
                WD(GI)=WD(GI) + N(ILOC,GI)*NWSUB(IDGNOD)
             ENDIF
          END DO

             IF(NSOGRASUB.NE.0) THEN
               IF(COGRAX) THEN
                 RDIRX=0.0
                 RDIRY=0.0
                 RDIRZ=1.0
               ELSE
                 RDIRX=XD(GI)
                 RDIRY=YD(GI)
                 RDIRZ=ZD(GI)
                 RN=SQRT(RDIRX**2+RDIRY**2+RDIRZ**2)
                 RDIRX=RDIRX/RN 
                 RDIRY=RDIRY/RN
                 RDIRZ=RDIRZ/RN
               ENDIF
               
               RRSX=0.0
               RRSY=0.0
               RRSZ=0.0
               DO ILOC=1,NLOC
                 RRSX=RRSX+N(ILOC,GI)*SOXGRASUB((ELE-1)*NLOC+ILOC)
                 RRSY=RRSY+N(ILOC,GI)*SOYGRASUB((ELE-1)*NLOC+ILOC)
                 RRSZ=RRSZ+N(ILOC,GI)*SOZGRASUB((ELE-1)*NLOC+ILOC)
               END DO
               RNRRS=SQRT(RRSX**2+RRSY**2+RRSZ**2)
               GRA_SIG_CHAN=SIGN(1.0,RDIRX*RRSX+RDIRY*RRSY+RDIRZ*RRSZ)
               
               SOUXGI(GI)=SOUXGI(GI) + RDIRX*GRA_SIG_CHAN*RNRRS
               SOUYGI(GI)=SOUYGI(GI) + RDIRY*GRA_SIG_CHAN*RNRRS
               SOUZGI(GI)=SOUZGI(GI) + RDIRZ*GRA_SIG_CHAN*RNRRS
             ENDIF

          RX(GI)=(-FSXGI(GI)+VD(GI)*OME2(GI)+ SOUXGI(GI) )
          RY(GI)=(-FSYGI(GI)-UD(GI)*OME2(GI)+ SOUYGI(GI) )
          RZ(GI)=(-FSZGI(GI) + SOUZGI(GI) )
          RXR(GI)=RTREST11(GI)*RX(GI)+RTREST12(GI)*RY(GI)+RTREST13(GI)*RZ(GI)
          RYR(GI)=RTREST21(GI)*RX(GI)+RTREST22(GI)*RY(GI)+RTREST23(GI)*RZ(GI)
          RZR(GI)=RTREST31(GI)*RX(GI)+RTREST32(GI)*RY(GI)+RTREST33(GI)*RZ(GI)
          RXC(GI)=(-FSXGIC(GI)+VD(GI)*OME2(GI) )
          RYC(GI)=(-FSYGIC(GI)-UD(GI)*OME2(GI) )
          RZC(GI)=(-FSZGIC(GI)  )
          RXRC(GI)=RTREST11(GI)*RXC(GI)+RTREST12(GI)*RYC(GI)+RTREST13(GI)*RZC(GI)
          RYRC(GI)=RTREST21(GI)*RXC(GI)+RTREST22(GI)*RYC(GI)+RTREST23(GI)*RZC(GI)
          RZRC(GI)=RTREST31(GI)*RXC(GI)+RTREST32(GI)*RYC(GI)+RTREST33(GI)*RZC(GI)
       END DO

       iloc_loop: DO ILOC=1,NLOC
          GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
          jloc_loop: DO JLOC=1,MLOC
             GLOBJ=(ELE-1)*MLOC+JLOC

             RRVECX=0.0
             RRVECY=0.0
             RRVECZ=0.0

             RNN=0.0
             DO GI=1,NGI

                RRVECX=RRVECX-N(ILOC,GI)*DETWEI(GI)*DGP(GLOBJ) &
                     *(RTREST11(GI)*MX(JLOC,GI)+RTREST12(GI)*MY(JLOC,GI)+RTREST13(GI)*MZ(JLOC,GI))
                RRVECY=RRVECY-N(ILOC,GI)*DETWEI(GI)*DGP(GLOBJ) &
                     *(RTREST21(GI)*MX(JLOC,GI)+RTREST22(GI)*MY(JLOC,GI)+RTREST23(GI)*MZ(JLOC,GI))
                RRVECZ=RRVECZ-N(ILOC,GI)*DETWEI(GI)*DGP(GLOBJ) &
                     *(RTREST31(GI)*MX(JLOC,GI)+RTREST32(GI)*MY(JLOC,GI)&
                     &+RTREST33(GI)*MZ(JLOC,GI))

                RNN=RNN+N(ILOC,GI)*DETWEI(GI)*M(JLOC,GI)

             END DO
             VECX(GLOBI)=VECX(GLOBI) +RRVECX
             VECY(GLOBI)=VECY(GLOBI) +RRVECY
             VECZ(GLOBI)=VECZ(GLOBI) +RRVECZ

             VECX2(GLOBI)=VECX2(GLOBI) +RRVECX
             VECY2(GLOBI)=VECY2(GLOBI) +RRVECY
             VECZ2(GLOBI)=VECZ2(GLOBI) +RRVECZ

             IF(NSUBVLOC.NE.0) THEN
                IDGNOD=(ELE-1)*NSUBVLOC+ILOC
                VECX_SGSADD(IDGNOD)=VECX_SGSADD(IDGNOD)+RRVECX
                VECY_SGSADD(IDGNOD)=VECY_SGSADD(IDGNOD)+RRVECY
                VECZ_SGSADD(IDGNOD)=VECZ_SGSADD(IDGNOD)+RRVECZ

                VECX_SGSADD2(IDGNOD)=VECX_SGSADD2(IDGNOD)+RRVECX
                VECY_SGSADD2(IDGNOD)=VECY_SGSADD2(IDGNOD)+RRVECY
                VECZ_SGSADD2(IDGNOD)=VECZ_SGSADD2(IDGNOD)+RRVECZ
                MLELE(IDGNOD)=MLELE(IDGNOD)+RNN
             ENDIF

             ML(GLOBI)=ML(GLOBI)+RNN

          END DO jloc_loop

          DO GI=1,NGI
             VECX(GLOBI)=VECX(GLOBI)+N(ILOC,GI)*DETWEI(GI)*RXR(GI)
             VECY(GLOBI)=VECY(GLOBI)+N(ILOC,GI)*DETWEI(GI)*RYR(GI)
             VECZ(GLOBI)=VECZ(GLOBI)+N(ILOC,GI)*DETWEI(GI)*RZR(GI)

             VECXF(GLOBI)=VECXF(GLOBI)+N(ILOC,GI)*DETWEI(GI)*RXRC(GI)
             VECYF(GLOBI)=VECYF(GLOBI)+N(ILOC,GI)*DETWEI(GI)*RYRC(GI)
             VECZF(GLOBI)=VECZF(GLOBI)+N(ILOC,GI)*DETWEI(GI)*RZRC(GI)
             IF(NSUBVLOC.NE.0) THEN
                IDGNOD=(ELE-1)*NSUBVLOC+ILOC
                VECX_SGSADD(IDGNOD)=VECX_SGSADD(IDGNOD)+N(ILOC,GI)*DETWEI(GI)*RXR(GI)
                VECY_SGSADD(IDGNOD)=VECY_SGSADD(IDGNOD)+N(ILOC,GI)*DETWEI(GI)*RYR(GI)
                VECZ_SGSADD(IDGNOD)=VECZ_SGSADD(IDGNOD)+N(ILOC,GI)*DETWEI(GI)*RZR(GI)
             ENDIF
          END DO
       END DO iloc_loop

       ! Put surface integrals into matrix and rhs*****************************
       if(MOMBC) then
          face_loop: DO IFACE=1,4
             COUNT2=(ELE-1)*5+IFACE
             ! NB colele4 is ordered in terms of faces.
             ELE2=COLELE4(COUNT2)

             ! Surface eleemnt of domain...
             DO SILOC=1,SNLOC
                SILOC2ILOC(SILOC)=LOCLIST(IFACE,SILOC)
             END DO
             ! Form OTHILOC from SILOC2ILOC
             IF(ELE2.NE.0) THEN
                CALL ELE2ELEFACE(OTHILOC,SILOC2ILOC,SNLOC,NLOC,TOTELE, &
                     XONDGL,NDGLNO,ELE,ELE2)
             ELSE
                DO ILOC=1,NLOC
                   OTHILOC(ILOC)=ILOC
                END DO
             ENDIF

             ! Form MLOCOTHILOC and MLOCSILOC2ILOC
             CALL MLOCELE2ELEFACE(MLOCOTHILOC,MLOCSILOC2ILOC,SNLOC,NLOC, &
                  SMLOC,MLOC,OTHILOC,SILOC2ILOC)

             ! PERFORM SURFACE INTEGRATION         !
             ! Form approximate surface normal (NORMX,NORMY,NORMZ)
             CALL DGSIMPLNORM(ELE,SILOC2ILOC,TOTELE,NLOC,SNLOC,XONDGL, &
                  X,Y,Z,XNONOD,NORMX,NORMY,NORMZ)

             CALL DGSDETNXLOC(ELE,SILOC2ILOC,XONDGL, &
                  TOTELE,NLOC,XNONOD,SNLOC,SNGI, &
                  X,Y,Z, &
                  SN,SNLX,SNLY, SWEIGH, SDETWE,SAREA, .TRUE.,.FALSE., &
                  SNORMXN,SNORMYN,SNORMZN, &
                  NORMX,NORMY,NORMZ, &
                  ISPHERE,SMLOC,SM,SMLX,SMLY)

             TOPSUF=.FALSE.
             DO SGI=1,SNGI
                SDIRX(SGI)=0.0
                SDIRY(SGI)=0.0
                SDIRZ(SGI)=0.0
                DO SILOC=1,SMLOC
                   ILOC=MLOCSILOC2ILOC(SILOC)
                   IGL=PGNDGLNO((ELE-1)*MLOC+ILOC)
                   SDIRX(SGI)=SDIRX(SGI)-SM(SILOC,SGI)*PNORMX(IGL)
                   SDIRY(SGI)=SDIRY(SGI)-SM(SILOC,SGI)*PNORMY(IGL)
                   SDIRZ(SGI)=SDIRZ(SGI)-SM(SILOC,SGI)*PNORMZ(IGL)
                END DO
                RN=SQRT(SDIRX(SGI)**2+SDIRY(SGI)**2+SDIRZ(SGI)**2)
                SDIRX(SGI)=SDIRX(SGI)/RN
                SDIRY(SGI)=SDIRY(SGI)/RN
                SDIRZ(SGI)=SDIRZ(SGI)/RN

                SRTVX(SGI)=SDIRX(SGI)*SDIRX(SGI)*SNORMXN(SGI) &
                     +SDIRX(SGI)*SDIRY(SGI)*SNORMYN(SGI) &
                     +SDIRX(SGI)*SDIRZ(SGI)*SNORMZN(SGI)
                SRTVY(SGI)=SDIRY(SGI)*SDIRX(SGI)*SNORMXN(SGI) &
                     +SDIRY(SGI)*SDIRY(SGI)*SNORMYN(SGI) &
                     +SDIRY(SGI)*SDIRZ(SGI)*SNORMZN(SGI)
                SRTVZ(SGI)=SDIRZ(SGI)*SDIRX(SGI)*SNORMXN(SGI) &
                     +SDIRZ(SGI)*SDIRY(SGI)*SNORMYN(SGI) &
                     +SDIRZ(SGI)*SDIRZ(SGI)*SNORMZN(SGI)
                SRESTX(SGI)=SNORMXN(SGI) - SRTVX(SGI)
                SRESTY(SGI)=SNORMYN(SGI) - SRTVY(SGI)
                SRESTZ(SGI)=SNORMZN(SGI) - SRTVZ(SGI)
                NDOTQ(SGI)=SNORMXN(SGI)*SDIRX(SGI)+SNORMYN(SGI)*SDIRY(SGI) &
                     +SNORMZN(SGI)*SDIRZ(SGI)
                INCOME(SGI)=0.0
                IF(NDOTQ(SGI).LT.0.0) INCOME(SGI)=1.0
                IF(ELE2.EQ.0) THEN
                   IF(NDOTQ(SGI).LT.-0.8) TOPSUF=.TRUE.
                ENDIF
             END DO

             ! Perform surface integration...
             DO SILOC=1,SNLOC
                ILOC=SILOC2ILOC(SILOC)
                GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
                DO SJLOC=1,SMLOC
                   JLOC=MLOCSILOC2ILOC(SJLOC)
                   GLOBJ=(ELE-1)*MLOC+JLOC
                   JLOC2=MLOCOTHILOC(JLOC)
                   GLOBJ2=(ELE2-1)*MLOC+JLOC2

                   INVLMX=0.0
                   INVLMY=0.0
                   INVLMZ=0.0
                   BOTHVLMX=0.0
                   BOTHVLMY=0.0
                   BOTHVLMZ=0.0
                   DO SGI=1,SNGI
                      RNN=SN(SILOC,SGI)*SM(SJLOC,SGI)*SDETWE(SGI)
                      BOTHVLMX=BOTHVLMX+SRESTX(SGI)*RNN
                      BOTHVLMY=BOTHVLMY+SRESTY(SGI)*RNN
                      BOTHVLMZ=BOTHVLMZ+SRESTZ(SGI)*RNN
                   END DO
                   ! Assume a zero inlet b.c when we have an inlet into
                   ! the domain -no surface integral
                   ! Include surface eleement at bottom and sides of domain only...

                   IF(ELE2.EQ.0) THEN
                      POTH=DGP(GLOBJ)
                   ELSE
                      POTH=DGP(GLOBJ2)
                   ENDIF
                   VECX(GLOBI)=VECX(GLOBI)+INVLMX*(DGP(GLOBJ)-POTH) &
                        +0.5*BOTHVLMX*(DGP(GLOBJ)-POTH)
                   VECY(GLOBI)=VECY(GLOBI)+INVLMY*(DGP(GLOBJ)-POTH) &
                        +0.5*BOTHVLMY*(DGP(GLOBJ)-POTH)
                   VECZ(GLOBI)=VECZ(GLOBI)+INVLMZ*(DGP(GLOBJ)-POTH) &
                        +0.5*BOTHVLMZ*(DGP(GLOBJ)-POTH)
                   VECX2(GLOBI)=VECX2(GLOBI)+INVLMX*(DGP(GLOBJ)-POTH) &
                        +0.5*BOTHVLMX*(DGP(GLOBJ)-POTH)
                   VECY2(GLOBI)=VECY2(GLOBI)+INVLMY*(DGP(GLOBJ)-POTH) &
                        +0.5*BOTHVLMY*(DGP(GLOBJ)-POTH)
                   VECZ2(GLOBI)=VECZ2(GLOBI)+INVLMZ*(DGP(GLOBJ)-POTH) &
                        +0.5*BOTHVLMZ*(DGP(GLOBJ)-POTH)
                   IF(NSUBVLOC.NE.0) THEN
                      IDGNOD=(ELE-1)*NSUBVLOC+ILOC
                      VECX_SGSADD(IDGNOD)=VECX_SGSADD(IDGNOD)+INVLMX*(DGP(GLOBJ)-POTH) &
                           +0.5*BOTHVLMX*(DGP(GLOBJ)-POTH)
                      VECY_SGSADD(IDGNOD)=VECY_SGSADD(IDGNOD)+INVLMY*(DGP(GLOBJ)-POTH) &
                           +0.5*BOTHVLMY*(DGP(GLOBJ)-POTH)
                      VECZ_SGSADD(IDGNOD)=VECZ_SGSADD(IDGNOD)+INVLMZ*(DGP(GLOBJ)-POTH) &
                           +0.5*BOTHVLMZ*(DGP(GLOBJ)-POTH)

                      VECX_SGSADD2(IDGNOD)=VECX_SGSADD2(IDGNOD)+INVLMX*(DGP(GLOBJ)-POTH) &
                           +0.5*BOTHVLMX*(DGP(GLOBJ)-POTH)
                      VECY_SGSADD2(IDGNOD)=VECY_SGSADD2(IDGNOD)+INVLMY*(DGP(GLOBJ)-POTH) &
                           +0.5*BOTHVLMY*(DGP(GLOBJ)-POTH)
                      VECZ_SGSADD2(IDGNOD)=VECZ_SGSADD2(IDGNOD)+INVLMZ*(DGP(GLOBJ)-POTH) &
                           +0.5*BOTHVLMZ*(DGP(GLOBJ)-POTH)
                   ENDIF
                END DO
             END DO

          END DO face_loop
       ENDIF

    END DO element_loop

    ! Smooth out the resulting force...
!    print *,'((INUSUB.NE.0).AND.(NSUBNVLOC.NE.0).AND.(NSOGRASUB.NE.0)):', &
!             ((INUSUB.NE.0).AND.(NSUBNVLOC.NE.0).AND.(NSOGRASUB.NE.0))
!    print *,'INUSUB,NSUBNVLOC,NSOGRASUB:', INUSUB,NSUBNVLOC,NSOGRASUB
!    print *,'PSMOOTH,FSMOOTH:',PSMOOTH,FSMOOTH
!    stop 234

    IF((INUSUB.NE.0).AND.(NSUBNVLOC.NE.0).AND.(NSOGRASUB.NE.0)) THEN
    IF(PSMOOTH.OR.FSMOOTH) THEN
!    IF(.true.) THEN
!    IF(.true.) THEN
    
      CALL SMOOTH_FORCE(NONODS,NSUBVLOC,NLOC,TOTELE, &
        NDGLNO, PSMOOTH, FSMOOTH, DGSMOOTH, &
        VECX,VECY,VECZ, VECX2,VECY2,VECZ2, ML,MLELE, &
        VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD, &
        VECX_SGSADD2,VECY_SGSADD2,VECZ_SGSADD2, &
        VECXF,VECYF,VECZF, &
! For sub CALHIGHNODAINN...
        X,Y,Z, XONDGL, XNONOD,MLOC,NGI, &
        N,NLX,NLY,NLZ,MLX,MLY,MLZ, WEIGHT, DETWEI, &
        ISPHERE)
        
    ENDIF
    ENDIF

    ewrite(1,*) 'finished verticaldg'

  END SUBROUTINE VERTICALDG
!
!
!
!
     SUBROUTINE SMOOTH_FORCE(NONODS,NSUBVLOC,NLOC,TOTELE, &
        NDGLNO, PSMOOTH, FSMOOTH, DGSMOOTH, &
        VECX,VECY,VECZ, VECX2,VECY2,VECZ2, ML,MLELE, &
        VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD, &
        VECX_SGSADD2,VECY_SGSADD2,VECZ_SGSADD2, &
        VECXF,VECYF,VECZF, &
! For sub CALHIGHNODAINN...
        X,Y,Z, XONDGL, XNONOD,MLOC,NGI, &
        N,NLX,NLY,NLZ,MLX,MLY,MLZ, WEIGHT, DETWEI, &
        ISPHERE)
     INTEGER NONODS,NSUBVLOC,NLOC,TOTELE,XNONOD,MLOC,NGI
     INTEGER NDGLNO(TOTELE*NLOC) 
     LOGICAL PSMOOTH, FSMOOTH
     INTEGER DGSMOOTH
     REAL VECX(NONODS),VECY(NONODS),VECZ(NONODS)
     REAL VECX2(NONODS),VECY2(NONODS),VECZ2(NONODS)
     REAL ML(NONODS),MLELE(NSUBVLOC*TOTELE)
    real, intent(inout):: VECX_SGSADD(TOTELE*NSUBVLOC),VECY_SGSADD(TOTELE*NSUBVLOC)
    real, intent(inout):: VECZ_SGSADD(TOTELE*NSUBVLOC)
    real, intent(inout):: VECX_SGSADD2(TOTELE*NSUBVLOC),VECY_SGSADD2(TOTELE*NSUBVLOC)
    real, intent(inout):: VECZ_SGSADD2(TOTELE*NSUBVLOC)
     REAL VECXF(NONODS),VECYF(NONODS),VECZF(NONODS)
! For sub CALHIGHNODAINN...
    real X(XNONOD),Y(XNONOD),Z(XNONOD)
    real, intent(in):: n(nloc,ngi),nlx(nloc,ngi),nly(nloc,ngi),nlz(nloc,ngi)
    real, intent(in):: weight(ngi),detwei(ngi)
    real, intent(in):: mlx(mloc,ngi),mly(mloc,ngi),mlz(mloc,ngi)
    integer, intent(in):: xondgl(totele*nloc)
    integer ISPHERE
! local variables...
    INTEGER GLOBI,GLOBJ,IDGNOD,JDGNOD,ELE,ILOC,JLOC,GI
    REAL RNN
       IF(PSMOOTH) THEN
          DO GLOBI=1,NONODS
             VECX(GLOBI)=VECX(GLOBI)-VECX2(GLOBI)
             VECY(GLOBI)=VECY(GLOBI)-VECY2(GLOBI)
             VECZ(GLOBI)=VECZ(GLOBI)-VECZ2(GLOBI)
          END DO
          DO GLOBI=1,NONODS
             VECX2(GLOBI)=VECX2(GLOBI)/ML(GLOBI)
             VECY2(GLOBI)=VECY2(GLOBI)/ML(GLOBI)
             VECZ2(GLOBI)=VECZ2(GLOBI)/ML(GLOBI)
          END DO
          IF(NSUBVLOC.NE.0) THEN
             DO IDGNOD=1,NSUBVLOC*TOTELE
                VECX_SGSADD(IDGNOD)=VECX_SGSADD(IDGNOD)-VECX_SGSADD2(IDGNOD)
                VECY_SGSADD(IDGNOD)=VECY_SGSADD(IDGNOD)-VECY_SGSADD2(IDGNOD)
                VECZ_SGSADD(IDGNOD)=VECZ_SGSADD(IDGNOD)-VECZ_SGSADD2(IDGNOD)
             END DO
             DO IDGNOD=1,NSUBVLOC*TOTELE
                VECX_SGSADD2(IDGNOD)=VECX_SGSADD2(IDGNOD)/MLELE(IDGNOD)
                VECY_SGSADD2(IDGNOD)=VECY_SGSADD2(IDGNOD)/MLELE(IDGNOD)
                VECZ_SGSADD2(IDGNOD)=VECZ_SGSADD2(IDGNOD)/MLELE(IDGNOD)
             END DO
          ENDIF
       ENDIF

       IF(FSMOOTH) THEN
          DO GLOBI=1,NONODS
             VECX(GLOBI)=VECX(GLOBI)-VECXF(GLOBI)
             VECY(GLOBI)=VECY(GLOBI)-VECYF(GLOBI)
             VECZ(GLOBI)=VECZ(GLOBI)-VECZF(GLOBI)
          END DO
          DO GLOBI=1,NONODS
             VECXF(GLOBI)=VECXF(GLOBI)/ML(GLOBI)
             VECYF(GLOBI)=VECYF(GLOBI)/ML(GLOBI)
             VECZF(GLOBI)=VECZF(GLOBI)/ML(GLOBI)
          END DO
       ENDIF

       DO ELE=1,TOTELE

          ! Calculate DETWEI for element ELE -NEW VERSION
          CALL CALHIGHNODAINN(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MLOC,NGI, &
               N,NLX,NLY,NLZ,MLX,MLY,MLZ, WEIGHT, DETWEI,.TRUE.,.FALSE., &
               ISPHERE)

          DO ILOC=1,NLOC
             GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
             DO JLOC=1,NLOC
                GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)
                RNN=0.0
                DO GI=1,NGI
                   RNN=RNN+N(ILOC,GI)*DETWEI(GI)*N(JLOC,GI)
                END DO
                IF(PSMOOTH) THEN
                   VECX(GLOBI)=VECX(GLOBI)+RNN*VECX2(GLOBJ)
                   VECY(GLOBI)=VECY(GLOBI)+RNN*VECY2(GLOBJ)
                   VECZ(GLOBI)=VECZ(GLOBI)+RNN*VECZ2(GLOBJ)
                   IF(NSUBVLOC.NE.0) THEN
                      IDGNOD=(ELE-1)*NLOC+ILOC
                      JDGNOD=(ELE-1)*NLOC+JLOC
                      IF(DGSMOOTH.EQ.2) THEN
                         ! Smoothing local to an element...
                         VECX_SGSADD(IDGNOD)=VECX_SGSADD(IDGNOD)+RNN*VECX_SGSADD2(JDGNOD)
                         VECY_SGSADD(IDGNOD)=VECY_SGSADD(IDGNOD)+RNN*VECY_SGSADD2(JDGNOD)
                         VECZ_SGSADD(IDGNOD)=VECZ_SGSADD(IDGNOD)+RNN*VECZ_SGSADD2(JDGNOD)
                      ELSE
                         ! Smoothing across elements...
                         VECX_SGSADD(IDGNOD)=VECX_SGSADD(IDGNOD)+RNN*VECX2(GLOBJ)
                         VECY_SGSADD(IDGNOD)=VECY_SGSADD(IDGNOD)+RNN*VECY2(GLOBJ)
                         VECZ_SGSADD(IDGNOD)=VECZ_SGSADD(IDGNOD)+RNN*VECZ2(GLOBJ)
                      ENDIF
                   ENDIF
                ENDIF
                IF(FSMOOTH) THEN
                   VECX(GLOBI)=VECX(GLOBI)+RNN*VECXF(GLOBJ)
                   VECY(GLOBI)=VECY(GLOBI)+RNN*VECYF(GLOBJ)
                   VECZ(GLOBI)=VECZ(GLOBI)+RNN*VECZF(GLOBJ)
                ENDIF
             END DO
          END DO

       END DO
      END SUBROUTINE SMOOTH_FORCE
!
!
!
!
   subroutine verticaldg_rhs(rhsvec, &
       ndglno, pgndglno, xondgl, &
       m, mlx, mly, mlz, &
       n, nlx, nly, nlz, &
       weight, &
       GEOBAL,&
       x, y, z, &
       u, v, w, &
       inusub,nsubnvloc,nusub,nvsub,nwsub, &
       NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB,COGRAX, &
       pgravx, pgravy, pgravz, &
       pnormx, pnormy, pnormz, &
       potgra, ndenst, frees, eqetid, &
       gravty, &
       isphere, scfacth0, idgboyonly, usepotgra, getfrees,totele,wdgravity)

    real, intent(out):: rhsvec(:)
    integer, dimension(:), intent(in):: ndglno, pgndglno, xondgl
    real, dimension(:,:), intent(in):: m, mlx, mly, mlz
    real, dimension(:,:), intent(in):: n, nlx, nly, nlz
    real, dimension(:), intent(in):: weight
    integer, intent(in):: geobal
    real, dimension(:), intent(in):: x, y, z
    real, dimension(:), intent(in):: u, v, w
    integer, intent(in):: inusub,nsubnvloc,NSOGRASUB
    real, dimension(:), intent(in):: nusub, nvsub, nwsub
    real, dimension(:), intent(in):: SOXGRASUB,SOYGRASUB,SOZGRASUB
    real, dimension(:), intent(in):: pgravx, pgravy, pgravz
    real, dimension(:), intent(in):: pnormx, pnormy, pnormz
    real, dimension(:), intent(in):: potgra, ndenst, frees, eqetid
    integer, intent(in):: isphere, idgboyonly, usepotgra
    real, intent(in):: scfacth0, gravty
    logical, intent(in):: getfrees,COGRAX

    logical, parameter:: D3=.true., DCYL=.false.

    real, dimension(size(m,2)):: detwei
    real, dimension(size(m,2)):: xd, yd, zd
    real, dimension(size(m,2)):: a11,a12,a13, a21,a22,a23, a31,a32,a33
    real rn
     REAL RRSX,RRSY,RRSZ,RNRRS,GRA_SIG_CHAN
     REAL RDIRX,RDIRY,RDIRZ

    real, dimension(size(m,2), size(m,1)):: mt
    real, dimension(size(n,2), size(n,1)):: nt
    real, dimension(size(m,2), 3, size(m,1)):: dm_t
    real, dimension(size(m,1), size(m,2), 3):: dnc
    real, dimension(3, size(m,1), size(m,2)):: dm
    real, dimension(3, size(m, 1)):: xyz
    real, dimension(3, 3, size(m,2)):: invJ
    real, dimension(size(m,2)):: detJ
    real, dimension(size(m,1)):: potgraloc, freespotloc
    real, dimension(size(m,1), 3):: pnorm, pgravloc
    real, dimension(size(m,2),3):: gipgrav, dir
    real, dimension(size(n,2),3):: ud, r
    real, dimension(size(n,1),3):: uloc
    real, dimension(size(n,1)):: ndenstloc
    real, dimension(size(n,2)):: omegi
    real, dimension(size(m,1), size(m,2)):: mx, my, mz
    real, dimension(size(n,1), size(n,2)):: nx, ny, nz
    integer totele, xnonod, nloc, mloc, xoloc, ncloc, ngi
    integer ndglno_index, pgndglno_index, rhsa, rhsb
    integer ele, gi, igl, i, iloc, idgnod,iglx

    real, intent(in) :: wdgravity(totele)

    ! sizes of things:
    mloc=size(m,1)
    nloc=size(n,1)
    ngi=size(m,2)
    ASSERT(ngi==size(n,2))
    !     totele=size(rhsvec)/mloc
    ASSERT(totele==size(rhsvec)/mloc)
    xnonod=size(x)
    xoloc=size(xondgl)/totele

    ! need transpose of m (of shape ngi x mloc)
    ! so that multiply with mloc x dim gives ngi x dim
    mt=transpose(m)
    ! similar for n
    nt=transpose(n)
    if (isphere==2) then
       ! this assumes same element for pressure as for coordinates !!!
       ncloc=mloc
       ! derivative of shape function for coordinates
       dnc(:,:,1)=mlx
       dnc(:,:,2)=mly
       dnc(:,:,3)=mlz
       ! derivative of shape function for pressure 
       ! (same thing but in different order)
       dm(1,:,:)=mlx
       dm(2,:,:)=mly
       dm(3,:,:)=mlz
    end if

    ! indices in resp. rhsvec, ndglno and pgndglno
    rhsa=1
    rhsb=mloc
    ndglno_index=0
    pgndglno_index=0

    elements: DO ELE=1,TOTELE

       if (isphere==2) then

          call unwind_xl_calc(xyz(1,:), xyz(2,:), xyz(3,:), & 
               ncloc,xnonod,xoloc,totele,xondgl, ele, x, y, z)
          call transform_to_physical_fast_3d(xyz, dnc, invJ, detJ)
          detwei=detJ*weight
          forall (gi=1:ngi)
             dm_t(gi,:,:)=matmul( invJ(:,:,gi), dm(:,:,gi) )
          end forall

       else

          CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MLOC,NGI, &
               N,NLX,NLY,NLZ, M,MLX,MLY,MLZ, WEIGHT, DETWEI,D3,DCYL, &
               NX,NY,NZ, mx, my, mz, &
               A11,A12,A13, A21,A22,A23, A31,A32,A33, &
               XD,YD,ZD, &
               ISPHERE)
          dm_t(:,1,:)=transpose(mx)
          dm_t(:,2,:)=transpose(my)
          dm_t(:,3,:)=transpose(mz)
       end if

       ! make vectors of size mloc for var. on pressure mesh
       do iloc=1, mloc
          IGL=PGNDGLNO(pgndglno_index+ILOC)
          if (usepotgra==1) then
             potgraloc(iloc)=potgra(igl)
          else
             pgravloc(iloc, 1)=pgravx(igl)
             pgravloc(iloc, 2)=pgravy(igl)
             pgravloc(iloc, 3)=pgravz(igl)
          end if
          pnorm(iloc, 1)=pnormx(igl)
          pnorm(iloc, 2)=pnormy(igl)
          pnorm(iloc, 3)=pnormz(igl)
          if (getfrees) then
             freespotloc(iloc)=(frees(igl)-eqetid(igl))*gravty*wdgravity(ele)
          end if
       end do

       ! Establish gravitational direction:
       ! mt (ngi x mloc) * pnorm (mloc x 3) gives ngi x 3
       dir=matmul(mt, pnorm)
       ! normalise on each gauss point:
       do gi=1, ngi
          rn=sqrt(sum(dir(gi,:)**2))
          dir(gi,:)=dir(gi,:)/rn
       end do

       ! Gravity source:
       if (usepotgra==1) then
          ! potgraloc (mloc) * dm_t (mloc x ngi*3) gives ngi*3 object
          ! reshaped in ngi x 3 object:
          gipgrav=reshape( &
               matmul(  reshape( dm_t, (/ ngi*3, mloc /) ), potgraloc), &
               (/ ngi, 3 /) )
       else
          ! mt (ngi x mloc) * pgravloc (mloc x 3) gives ngi x 3
          gipgrav=matmul(mt, pgravloc)
       end if

       ! Bouyancy term:
       do iloc=1, nloc
          igl=ndglno(ndglno_index+iloc)
          ndenstloc(iloc)=ndenst(igl)
          IF((INUSUB.NE.0).AND.(NSUBNVLOC.NE.0)) THEN
             IDGNOD=(ELE-1)*NLOC+ILOC
             uloc(iloc,1)=u(igl)+NUSUB(IDGNOD)
             uloc(iloc,2)=v(igl)+NVSUB(IDGNOD)
             uloc(iloc,3)=w(igl)+NWSUB(IDGNOD)
          ELSE
             uloc(iloc,1)=u(igl)
             uloc(iloc,2)=v(igl)
             uloc(iloc,3)=w(igl)
          ENDIF
       end do
       do i=1, 3
          r(:,i)=matmul(ndenstloc, n)*gipgrav(:,i)
       end do

        
           IF(NSOGRASUB.NE.0) THEN
             XD=0.0
             YD=0.0
             ZD=0.0
             DO ILOC=1,NLOC
               iglx=XONDGL((ele-1)*nloc+iloc)
               do gi=1,ngi
                 XD(GI)=XD(GI)+N(ILOC,GI)*X(IGLX)
                 YD(GI)=YD(GI)+N(ILOC,GI)*Y(IGLX)
                 ZD(GI)=ZD(GI)+N(ILOC,GI)*Z(IGLX)
               END DO
             END DO
             do gi=1,ngi
               DO ILOC=1,NLOC
                 RRSX=RRSX+N(ILOC,GI)*SOXGRASUB((ELE-1)*NLOC+ILOC)
                 RRSY=RRSY+N(ILOC,GI)*SOYGRASUB((ELE-1)*NLOC+ILOC)
                 RRSZ=RRSZ+N(ILOC,GI)*SOZGRASUB((ELE-1)*NLOC+ILOC)
               END DO
               IF(COGRAX) THEN
                 RDIRX=0.0
                 RDIRY=0.0
                 RDIRZ=1.0
               ELSE
                 RDIRX=XD(GI)
                 RDIRY=YD(GI)
                 RDIRZ=ZD(GI)
                 RN=SQRT(RDIRX**2+RDIRY**2+RDIRZ**2)
                 RDIRX=RDIRX/RN 
                 RDIRY=RDIRY/RN
                 RDIRZ=RDIRZ/RN
               ENDIF
               RRSX=0.0
               RRSY=0.0
               RRSZ=0.0
               DO ILOC=1,NLOC
                 RRSX=RRSX+N(ILOC,GI)*SOXGRASUB((ELE-1)*NLOC+ILOC)
                 RRSY=RRSY+N(ILOC,GI)*SOYGRASUB((ELE-1)*NLOC+ILOC)
                 RRSZ=RRSZ+N(ILOC,GI)*SOZGRASUB((ELE-1)*NLOC+ILOC)
               END DO
               RNRRS=SQRT(RRSX**2+RRSY**2+RRSZ**2)
               GRA_SIG_CHAN=SIGN(1.0,RDIRX*RRSX+RDIRY*RRSY+RDIRZ*RRSZ)
               
               r(gi,1)=r(gi,1) + RDIRX*GRA_SIG_CHAN*RNRRS
               r(gi,2)=r(gi,2) + RDIRY*GRA_SIG_CHAN*RNRRS
               r(gi,3)=r(gi,3) + RDIRZ*GRA_SIG_CHAN*RNRRS
             end do
           ENDIF
        
       ! Free surface contribution:
       if (getfrees) then

          ! dm_t (ngi*3 x mloc) * freespotloc (mloc) gives ngi*3 object
          r=r-reshape( &
               matmul(  reshape( dm_t, (/ ngi*3, mloc /) ), freespotloc), &
               (/ ngi, 3 /) )
       end if

       ! Coriolis contribution:
       ! nt (ngi x nloc) * uloc (nloc x 3) gives ud (ngi x 3)
       ud=matmul(nt, uloc)

       ! Coriolis using other coordinates:
       do gi=1, ngi
          omegi(gi)=funome(xd(gi), yd(gi), zd(gi))
       end do
       r(:,1)=r(:,1)+ud(:,2)*omegi*2
       r(:,2)=r(:,2)-ud(:,1)*omegi*2

       ! inner product of grav. dir. normal vector dir (ngi x 3)
       ! and force vector, composed above, r (ngi x 3) gives ngi entries
       ! integrated with shape functions m (mloc x ngi) gives mloc contributions
       rhsvec( rhsa:rhsb )=-matmul(m, sum(dir*r,2)*detwei)

       rhsa=rhsb+1
       rhsb=rhsb+mloc
       ndglno_index=ndglno_index+nloc
       pgndglno_index=pgndglno_index+mloc

    end do elements

  end subroutine verticaldg_rhs

  subroutine verticaldg_solve(dgp, rhsvec, &
       tolercent, max_its, &
       eleord, colele4, &
       ndglno, pgndglno, xondgl, &
       m, mlx, mly, mlz, &
       n, nlx, nly, nlz, &
       weight, &
       x, y, z, &
       pnormx, pnormy, pnormz, &
       isphere, &
       sn, snlx, snly, &
       sm, smlx, smly, &
       sweigh)

    real, dimension(:), intent(inout):: dgp
    real, dimension(:), intent(in):: rhsvec
    real, intent(in):: tolercent
    integer, intent(in):: max_its
    integer, dimension(:), intent(in):: eleord, colele4
    integer, dimension(:), intent(in):: ndglno, pgndglno, xondgl
    real, dimension(:,:), intent(in):: m, mlx, mly, mlz
    real, dimension(:,:), intent(in):: n, nlx, nly, nlz
    real, dimension(:), intent(in):: weight
    real, dimension(:), intent(in):: x, y, z
    real, dimension(:), intent(in):: pnormx, pnormy, pnormz
    integer, intent(in):: isphere

    logical, parameter:: D3=.true., DCYL=.false.
    integer, parameter:: snloc=3, smloc=6, sngi=7

    real, dimension(snloc, sngi), intent(in):: sn, snlx, snly
    real, dimension(smloc, sngi), intent(in):: sm, smlx, smly
    real, dimension(sngi), intent(in):: sweigh

    real, dimension( size(m,1), size(m,1) ) :: mat
    real, dimension( size(m,1) ) :: rhsloc, locdgp
    real, dimension( 3, size(m,1) ) :: xyz
    real, dimension( size(m,1), size(m,2), 3 ) :: dnc
    real, dimension( 3, 3, size(m,2) ) :: invJ
    real, dimension( size(m,2) ) :: detJ, detwei
    real, dimension( 3, size(m,1), size(m,2) ) :: dm, dm_t
    real, dimension( size(m,1), 3 ) :: pnorm
    real, dimension( size(m,2), 3 ) :: dir
    real, dimension( size(m,2), size(m,1) ) :: dir_dm_t, mt
    real rn
    integer its, i, gi, ele, igl
    integer mloc, nloc, ncloc, xoloc,ngi
    integer totele, xnonod
    ! stuff needed for calhighnodai
    real, dimension(size(m,2)):: xd, yd, zd
    real, dimension(size(m,2)):: a11,a12,a13, a21,a22,a23, a31,a32,a33
    real, dimension(size(n,1), size(n,2)):: nx, ny, nz
    real, dimension(size(m,1), size(m,2)):: mx, my, mz

    real maxdiff, maxdgp
    real normx, normy, normz, sarea, invlm, poth
    integer, parameter:: loclist(4,snloc)=reshape( &
         (/ 1,1,3,1, 2,2,2,3, 3,4,4,4 /), (/ 4,3 /) )
    integer, dimension(snloc)::  siloc2iloc
    integer, dimension(size(ndglno)/size(eleord)):: othiloc
    integer, dimension(smloc):: mlocsiloc2iloc
    integer, dimension(size(m,1)):: mlocothiloc
    real, dimension(sngi):: sdetwe, snormxn, snormyn, snormzn
    real, dimension(sngi):: sdirx, sdiry, sdirz, ndotq, income
    integer iface, count2, ele2, siloc, sjloc, iloc, iloc2, globi, globi2
    integer jloc, jloc2, globj, globj2, sgi
    logical:: topsuf 

    ! sizes of things:
    totele=size(eleord)
    xnonod=size(x)
    mloc=size(m,1)
    nloc=size(ndglno)/totele
    xoloc=size(xondgl)/totele
    ncloc=mloc
    ngi=size(m,2)

    mt=transpose(m)
    ! derivative of shape function for coordinates
    dnc(:,:,1)=mlx
    dnc(:,:,2)=mly
    dnc(:,:,3)=mlz
    ! derivative of shape function for pressure 
    ! (same thing but in different order)
    dm(1,:,:)=mlx
    dm(2,:,:)=mly
    dm(3,:,:)=mlz

    iterations: do its=1, max_its
       ewrite(2,*) 'its=',its

       maxdgp=0.0
       maxdiff=0.0
       elements: DO i=1, totele
          ele=eleord(i)

          rhsloc(1:mloc) = rhsvec( (ele-1)*mloc+1:ele*mloc )

          if (isphere==2) then

             call unwind_xl_calc(xyz(1,:), xyz(2,:), xyz(3,:), & 
                  ncloc,xnonod,xoloc,totele,xondgl, ele, x, y, z)
             call transform_to_physical_fast_3d(xyz, dnc, invJ, detJ)
             detwei=detJ*weight
             forall (gi=1:ngi)
                dm_t(:,:,gi)=matmul( invJ(:,:,gi), dm(:,:,gi) )
             end forall

          else

             CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MLOC,NGI, &
                  N,NLX,NLY,NLZ, M,MLX,MLY,MLZ, WEIGHT, DETWEI,D3,DCYL, &
                  NX,NY,NZ, mx, my, mz, &
                  A11,A12,A13, A21,A22,A23, A31,A32,A33, &
                  XD,YD,ZD, &
                  ISPHERE)
             dm_t(1,:,:)=mx
             dm_t(2,:,:)=my
             dm_t(3,:,:)=mz
          end if

          ! make pnorm vector of size mloc
          do iloc=1, mloc
             igl=pgndglno((ele-1)*mloc+iloc)
             pnorm(iloc, 1)=pnormx(igl)
             pnorm(iloc, 2)=pnormy(igl)
             pnorm(iloc, 3)=pnormz(igl)
          end do

          ! Establish gravitational direction:
          ! mt (ngi x mloc) * pnorm (mloc x 3) gives ngi x 3
          dir=-matmul(mt, pnorm)
          ! normalise on each gauss point:
          do gi=1, ngi
             rn=sqrt(sum(dir(gi,:)**2))
             dir(gi,:)=dir(gi,:)*(detwei(gi)/rn)
             ! dir ([ngi x] 3) * dm_t (3 x mloc [x ngi]) gives [ngi x] mloc
             dir_dm_t(gi,:)=matmul( dir(gi, :), dm_t(:, :, gi) )
          end do

          ! m (mloc x ngi) * dir_dm_t (ngi x mloc) gives mat (mloc x mloc)
          mat=matmul( m, dir_dm_t )

          ! Put surface integrals into matrix and rhs*****************************
          faces: DO IFACE=1,4
             COUNT2=(ELE-1)*5+IFACE
             ! NB colele4 is ordered in terms of faces.
             ELE2=COLELE4(COUNT2)

             ! Surface element of domain...
             DO SILOC=1,SNLOC
                SILOC2ILOC(SILOC)=LOCLIST(IFACE,SILOC)
             END DO
             ! Form OTHILOC from SILOC2ILOC
             IF(ELE2.NE.0) THEN
                CALL ELE2ELEFACE(OTHILOC,SILOC2ILOC,SNLOC,NLOC,TOTELE, &
                     XONDGL,NDGLNO,ELE,ELE2)
             ELSE
                DO ILOC=1,NLOC
                   OTHILOC(ILOC)=ILOC
                END DO
             ENDIF

             ! Form MLOCOTHILOC and MLOCSILOC2ILOC
             CALL MLOCELE2ELEFACE(MLOCOTHILOC,MLOCSILOC2ILOC,SNLOC,NLOC, &
                  SMLOC,MLOC,OTHILOC,SILOC2ILOC)

             ! PERFORM SURFACE INTEGRATION         !
             ! Form approximate surface normal (NORMX,NORMY,NORMZ)
             CALL DGSIMPLNORM(ELE,SILOC2ILOC,TOTELE,NLOC,SNLOC,XONDGL, &
                  X,Y,Z,XNONOD,NORMX,NORMY,NORMZ)

             CALL DGSDETNXLOC(ELE,SILOC2ILOC,XONDGL, &
                  TOTELE,NLOC,XNONOD,SNLOC,SNGI, &
                  X,Y,Z, &
                  SN,SNLX,SNLY, SWEIGH, SDETWE,SAREA, .TRUE.,.FALSE., &
                  SNORMXN,SNORMYN,SNORMZN, &
                  NORMX,NORMY,NORMZ, &
                  ISPHERE,SMLOC,SM,SMLX,SMLY)

             TOPSUF=.FALSE.
             DO SGI=1,SNGI
                SDIRX(SGI)=0.0
                SDIRY(SGI)=0.0
                SDIRZ(SGI)=0.0
                DO SILOC=1,SMLOC
                   ILOC=MLOCSILOC2ILOC(SILOC)
                   IGL=PGNDGLNO((ELE-1)*MLOC+ILOC)
                   SDIRX(SGI)=SDIRX(SGI)-SM(SILOC,SGI)*PNORMX(IGL)
                   SDIRY(SGI)=SDIRY(SGI)-SM(SILOC,SGI)*PNORMY(IGL)
                   SDIRZ(SGI)=SDIRZ(SGI)-SM(SILOC,SGI)*PNORMZ(IGL)
                END DO
                RN=SQRT(SDIRX(SGI)**2+SDIRY(SGI)**2+SDIRZ(SGI)**2)
                SDIRX(SGI)=SDIRX(SGI)/RN
                SDIRY(SGI)=SDIRY(SGI)/RN
                SDIRZ(SGI)=SDIRZ(SGI)/RN
                NDOTQ(SGI)=SNORMXN(SGI)*SDIRX(SGI)+SNORMYN(SGI)*SDIRY(SGI) &
                     +SNORMZN(SGI)*SDIRZ(SGI)
                INCOME(SGI)=0.0
                IF(NDOTQ(SGI).LT.0.0) INCOME(SGI)=1.0
                IF(ELE2.EQ.0) THEN
                   IF(NDOTQ(SGI).LT.-0.8) TOPSUF=.TRUE.
                ENDIF
             END DO

             ! Perform surface integration...
             DO SILOC=1,SMLOC
                ILOC=MLOCSILOC2ILOC(SILOC)
                GLOBI=(ELE-1)*MLOC+ILOC
                ILOC2=MLOCOTHILOC(ILOC)
                GLOBI2=(ELE2-1)*MLOC+ILOC2
                DO SJLOC=1,SMLOC
                   JLOC=MLOCSILOC2ILOC(SJLOC)
                   GLOBJ=(ELE-1)*MLOC+JLOC
                   JLOC2=MLOCOTHILOC(JLOC)
                   GLOBJ2=(ELE2-1)*MLOC+JLOC2

                   INVLM =0.0
                   DO SGI=1,SNGI
                      INVLM= INVLM +INCOME(SGI)*NDOTQ(SGI)*SM(SILOC,SGI)*SM(SJLOC,SGI)*SDETWE(SGI)
                   END DO
                   IF(TOPSUF) THEN
                      POTH=0.0
                   ELSE IF(ELE2.EQ.0) THEN
                      POTH=DGP(GLOBJ)
                      INVLM=0.0
                   ELSE
                      POTH=DGP(GLOBJ2)
                   ENDIF
                   MAT(ILOC,JLOC)=MAT(ILOC,JLOC)-INVLM
                   rhsloc(ILOC)=rhsloc(ILOC)-INVLM*POTH

                END DO
             END DO

          END DO faces

          CALL SMLINN(MAT,LOCDGP,rhsloc,MLOC,MLOC)

          DO ILOC=1,MLOC
             MAXDIFF=MAX(ABS(DGP((ELE-1)*MLOC+ILOC)-LOCDGP(ILOC)),MAXDIFF)
             MAXDGP=MAX(ABS(LOCDGP(ILOC)),MAXDGP)
             DGP((ELE-1)*MLOC+ILOC)=LOCDGP(ILOC)
          END DO

       END DO elements
       ewrite(2,*) 'maxdiff,maxdgp,TOLERCENT*MAXDGP:',maxdiff,maxdgp,TOLERCENT*MAXDGP
       IF(MAXDIFF.LE.TOLERCENT*MAXDGP) exit
    end do iterations

  end subroutine verticaldg_solve

  SUBROUTINE VERTICALDGSTORE(VECX,VECY,VECZ, &
       NSUBVLOC,VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD, &
       DGP,ELEORD,GETELEORD,NLEVEL, &
       GETFREES,FREES,EQETID,GRAVTY,  &
       PNORMX,PNORMY,PNORMZ, &
       PGRAVX,PGRAVY,PGRAVZ,NDENST, &
       POTGRA,USEPOTGRA,  &
       PGNODS, &
       X,Y,Z, U,V,W, &
       INUSUB,NSUBNVLOC,NUSUB,NVSUB,NWSUB, &
       NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB,COGRAX, &
       M,MLX,MLY,MLZ,&
       N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC, &
       TOTELE, &
       GEOBAL,&
       NDGLNO,PGNDGLNO, XONDGL, &
       NONODS,XNONOD, &
       scfacth0,ISPHERE, &
       FINDRM,COLM,NCOLM, &
       IDGBCMOM,IDGBOYONLY,PSMOOTH,COLELE4,PHSOLVE,wdgravity )
    ! This subroutine solves the DG eqns for balanced pressure DGP.
    ! for speed this sub stores a large DG matrix.
    ! If PHSOLVE then solve the vertical DG eqns else use the previous 
    ! hydrostic pressure value. 
    IMPLICIT NONE
    ! Number of source nodes, advection nodes and node positions.
    ! Number of elements, velocity nodes, ???
    ! If GEOBAL=-23 then have quadratic variation in the vertical with a DG method.
    ! If GEOBAL=-24 then have quadratic variation in the vertical with a DG method but
    ! use no b.c's when putting into the momentum eqns.
    ! If GEOBAL=-25 same as -24 but do not solve for verticaldg pressure 
    ! but use a quadratic pressure i.e. P1_(SGS or DG)-P2. 
    ! If GEOBAL=-26 same as -24 but do not solve for verticaldg pressure 
    ! & solve for p and free surface together.
    ! If GEOBAL=-27 same as -23 but solve for p and free surface together.
    ! IF NMIS_FORM_MAT_DGP is the no of time this sub is 
    ! visited when PHSOLVE=.true. before we re-form the vertical dg matrix 
    ! - this matrix is stored between visits to this sub. IF=0,1 then re-form every visit. 
    REAL TOLERCENT
    PARAMETER(TOLERCENT=0.00001)
    INTEGER DGSMOOTH,NMIS_FORM_MAT_DGP
    PARAMETER(DGSMOOTH=1,NMIS_FORM_MAT_DGP=1)
    ! DGSMOOTH=1 for smoothing across elements.
    ! DGSMOOTH=1 for smoothing within element.
    INTEGER, intent(in) :: TOTELE,NONODS,XNONOD,PGNODS
    INTEGER, intent(in) :: NGI,NLOC,MLOC

    REAL DGP(TOTELE*MLOC)
    INTEGER ELEORD(TOTELE)
    LOGICAL GETELEORD
    INTEGER NLEVEL
    REAL SCFACTH0
    REAL VECX(NONODS),VECY(NONODS),VECZ(NONODS)
    INTEGER NSUBVLOC
    real, intent(inout):: VECX_SGSADD(TOTELE*NSUBVLOC),VECY_SGSADD(TOTELE*NSUBVLOC)
    real, intent(inout):: VECZ_SGSADD(TOTELE*NSUBVLOC)
    REAL U(NONODS),V(NONODS),W(NONODS)
    INTEGER INUSUB,NSUBNVLOC
    REAL NUSUB(TOTELE*NSUBNVLOC),NVSUB(TOTELE*NSUBNVLOC)
    REAL NWSUB(TOTELE*NSUBNVLOC)
  INTEGER NSOGRASUB
  REAL SOXGRASUB(TOTELE*NSOGRASUB),SOYGRASUB(TOTELE*NSOGRASUB)
  REAL SOZGRASUB(TOTELE*NSOGRASUB)
     LOGICAL COGRAX
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    LOGICAL GETFREES
    REAL FREES(PGNODS),EQETID(PGNODS)
    ! Global node number maps.
    INTEGER NDGLNO(TOTELE*NLOC),PGNDGLNO(TOTELE*MLOC)
    INTEGER XONDGL(TOTELE*NLOC)
    INTEGER FINDRM(NONODS+1),NCOLM,COLM(NCOLM)

    INTEGER USEPOTGRA,ISPHERE,IDGBCMOM,IDGBOYONLY
    REAL PNORMX(PGNODS),PNORMY(PGNODS),PNORMZ(PGNODS)
    REAL PGRAVX(PGNODS),PGRAVY(PGNODS),PGRAVZ(PGNODS), NDENST(NONODS)
    REAL POTGRA(PGNODS*USEPOTGRA)

    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL WEIGHT(NGI)

    REAL M(MLOC,NGI),MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
    INTEGER GEOBAL
    REAL GRAVTY

    Real, intent(in) :: wdgravity(totele)

    ! Local variables...
    INTEGER SMLOC,SNGI,SNLOC
    PARAMETER(SNLOC=3,SMLOC=6,SNGI=7)
    REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI),SNLZ(SNLOC,SNGI)
    REAL SWEIGH(SNGI)
    REAL SM(SMLOC,SNGI),SMLX(SMLOC,SNGI),SMLY(SMLOC,SNGI),SMLZ(SMLOC,SNGI)
    REAL SL1(SNGI), SL2(SNGI), SL3(SNGI), SL4(SNGI)
    REAL UD(NGI),VD(NGI),WD(NGI)
    REAL XD(NGI),YD(NGI),ZD(NGI)
    REAL DETWEI(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    REAL MX(MLOC,NGI),MY(MLOC,NGI),MZ(MLOC,NGI)
    REAL SDETWE(SNGI)
   
    REAL FSXGI(NGI),FSYGI(NGI),FSZGI(NGI)
    REAL FSXGIC(NGI),FSYGIC(NGI),FSZGIC(NGI)
    REAL GIPGRAVX(NGI),GIPGRAVY(NGI),GIPGRAVZ(NGI)
    REAL A11(NGI),A12(NGI),A13(NGI)
    REAL A21(NGI),A22(NGI),A23(NGI)
    REAL A31(NGI),A32(NGI),A33(NGI)
    REAL DIRX(NGI),DIRY(NGI),DIRZ(NGI)
    REAL OME2(NGI)
    REAL SOUXGI(NGI),SOUYGI(NGI),SOUZGI(NGI)
    REAL RX(NGI),RY(NGI),RZ(NGI)
    REAL RXC(NGI),RYC(NGI),RZC(NGI)

    REAL SNORMXN(SNGI),SNORMYN(SNGI),SNORMZN(SNGI)
    REAL SDIRX(SNGI),SDIRY(SNGI),SDIRZ(SNGI)
    REAL NDOTQ(SNGI),INCOME(SNGI)

    REAL SRTVX(SNGI),SRTVY(SNGI),SRTVZ(SNGI)
    REAL SRESTX(SNGI),SRESTY(SNGI),SRESTZ(SNGI)
    REAL RTREST11(NGI),RTREST12(NGI),RTREST13(NGI)
    REAL RTREST21(NGI),RTREST22(NGI),RTREST23(NGI)
    REAL RTREST31(NGI),RTREST32(NGI),RTREST33(NGI)
    REAL RXR(NGI),RYR(NGI),RZR(NGI)
    REAL RXRC(NGI),RYRC(NGI),RZRC(NGI)
    REAL LOCDGP(MLOC)

    INTEGER SILOC2ILOC(SNLOC),OTHILOC(NLOC)
    INTEGER MLOCSILOC2ILOC(SMLOC),MLOCOTHILOC(MLOC)

    INTEGER ELE,ILOC,JLOC
    INTEGER GI,SGI,GLOBI,GLOBJ,GLOBI2,GLOBJ2
    INTEGER COUNT2,IDGNOD
    LOGICAL TOPSUF,MOMBC

    LOGICAL PSMOOTH,FSMOOTH
    INTEGER COLELE4(TOTELE*5)
    LOGICAL PHSOLVE

    REAL NORMX,NORMY,NORMZ
    REAL SAREA
    INTEGER ELE2,ELE3,II,adjits
    INTEGER ILOC2,JLOC2,IGL,NITS,ITS 
    INTEGER SILOC,SJLOC
    LOGICAL GET_NEW_MAT
    REAL INVLM,RNN
    REAL RN,INVLMX,INVLMY,INVLMZ,BOTHVLMX,BOTHVLMY,BOTHVLMZ,POTH
    REAL RRVECX,RRVECY,RRVECZ
    REAL MAXDGP,MAXDIFF
    REAL RRSX,RRSY,RRSZ,RNRRS,GRA_SIG_CHAN
     REAL RDIRX,RDIRY,RDIRZ
    INTEGER IFACE,NFACE
    PARAMETER(NFACE=4)
    INTEGER LOCLIST(NFACE,3)
    ! Only need to sweep through the elements one 
    ! if one_vdg_ser_it=1 else iterate to convergence. 
    INTEGER one_vdg_ser_it,COUNT_MIS_FORM_MAT_DGP
    SAVE one_vdg_ser_it
    DATA one_vdg_ser_it /0/
    SAVE COUNT_MIS_FORM_MAT_DGP
    DATA COUNT_MIS_FORM_MAT_DGP /0/

    REAL MAT(MLOC,MLOC),VECLOC(MLOC)
    REAL, ALLOCATABLE, DIMENSION(:)::RHSVEC
    REAL, ALLOCATABLE, DIMENSION(:)::VECX2
    REAL, ALLOCATABLE, DIMENSION(:)::VECY2
    REAL, ALLOCATABLE, DIMENSION(:)::VECZ2
    REAL, ALLOCATABLE, DIMENSION(:)::VECXF
    REAL, ALLOCATABLE, DIMENSION(:)::VECYF
    REAL, ALLOCATABLE, DIMENSION(:)::VECZF
    REAL, ALLOCATABLE, DIMENSION(:)::ML
    real, allocatable, dimension(:)::VECX_SGSADD2
    real, allocatable, dimension(:)::VECY_SGSADD2
    real, allocatable, dimension(:)::VECZ_SGSADD2
    real, allocatable, dimension(:)::MLELE
    
    REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:)::MATELE
    REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)::MATFACE
    INTEGER, ALLOCATABLE, SAVE, DIMENSION(:,:,:)::ILOCMATFACE
    INTEGER, ALLOCATABLE, SAVE, DIMENSION(:,:,:)::JLOC2MATFACE

    ewrite(1,*) 'in VERTICALDGSTORE'

    ALLOCATE(RHSVEC(TOTELE*MLOC))
    ALLOCATE(VECX2(NONODS))
    ALLOCATE(VECY2(NONODS))
    ALLOCATE(VECZ2(NONODS))
    ALLOCATE(VECXF(NONODS))
    ALLOCATE(VECYF(NONODS))
    ALLOCATE(VECZF(NONODS))
    ALLOCATE(ML(NONODS))
    ALLOCATE(VECX_SGSADD2(NSUBVLOC*TOTELE))
    ALLOCATE(VECY_SGSADD2(NSUBVLOC*TOTELE))
    ALLOCATE(VECZ_SGSADD2(NSUBVLOC*TOTELE))
    ALLOCATE(MLELE(NSUBVLOC*TOTELE))

    VECX2(1:NONODS) = 0.0
    VECY2(1:NONODS) = 0.0
    VECZ2(1:NONODS) = 0.0
    VECXF(1:NONODS) = 0.0
    VECYF(1:NONODS) = 0.0
    VECZF(1:NONODS) = 0.0
    ML(1:NONODS) = 0.0
    VECX_SGSADD2 = 0.0
    VECY_SGSADD2 = 0.0
    VECZ_SGSADD2 = 0.0
    MLELE=0.0
    GET_NEW_MAT=.TRUE.

    IF(GEOBAL.LT.-30) THEN
       MOMBC=(IDGBCMOM.EQ.1)
    ELSE
       MOMBC=.TRUE.
       IF(GEOBAL.EQ.-24) MOMBC=.FALSE.
       IF(GEOBAL.EQ.-25) MOMBC=.FALSE.
       IF(GEOBAL.EQ.-26) MOMBC=.FALSE.
    ENDIF

    FSMOOTH=.FALSE.

    ! Calculate surface shape functions SN, SM etc...
    CALL TRIQUA(SL1, SL2, SL3, SL3, SWEIGH, .FALSE.,SNGI)
    CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
         SNLOC,SNGI, &
         SN,SNLX,SNLY,SNLZ)
    CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
         SMLOC,SNGI, &
         SM,SMLX,SMLY,SMLZ)

    ! The local cords are in anti-clockwise order
    IFACE=1
    LOCLIST(IFACE,1)=1
    LOCLIST(IFACE,2)=2
    LOCLIST(IFACE,3)=3
    IFACE=2
    LOCLIST(IFACE,1)=1
    LOCLIST(IFACE,2)=2
    LOCLIST(IFACE,3)=4
    IFACE=3
    LOCLIST(IFACE,1)=3
    LOCLIST(IFACE,2)=2
    LOCLIST(IFACE,3)=4
    IFACE=4
    LOCLIST(IFACE,1)=1
    LOCLIST(IFACE,2)=3
    LOCLIST(IFACE,3)=4

    ! CALCULATE ELEMENT TO ELEMENT LIST**************
    ! and put in order of faces

    ewrite(1,*) 'going into 1ST ELEMENT LOOP PHSOLVE=',PHSOLVE
    IF(PHSOLVE) THEN 

       ! Alocate matrix...
       ewrite(2,*) 'matface mem:',TOTELE*4*sMLOC*sMLOC

       IF(GETELEORD) THEN

          ! *********FIND ELEMENT ORDERING******************************

          CALL ELEMENT_ORDER(ELEORD,NLEVEL,TOTELE,SNLOC,NLOC,SMLOC,MLOC,SNGI, &
               XONDGL,NDGLNO, &
               COLELE4,LOCLIST,PGNODS,NONODS,XNONOD,X,Y,Z, &
               SN,SNLX,SNLY, SWEIGH, &
               ISPHERE,SM,SMLX,SMLY, &
               FINDRM,COLM,NCOLM, &
               PGNDGLNO,PNORMX,PNORMY,PNORMZ,&
               one_vdg_ser_it,adjits)
    ewrite(1,*) 'one_vdg_ser_it,adjits=',one_vdg_ser_it,adjits
          COUNT_MIS_FORM_MAT_DGP=10000
       END IF
    
       COUNT_MIS_FORM_MAT_DGP=COUNT_MIS_FORM_MAT_DGP+1
       IF(COUNT_MIS_FORM_MAT_DGP.GE.NMIS_FORM_MAT_DGP) THEN
         GET_NEW_MAT=.TRUE.
         COUNT_MIS_FORM_MAT_DGP=0
       ELSE
         GET_NEW_MAT=.FALSE.
       ENDIF
       
    ewrite(1,*) 'COUNT_MIS_FORM_MAT_DGP,NMIS_FORM_MAT_DGP,GET_NEW_MAT,AllOCATED(MATFACE):', &
                 COUNT_MIS_FORM_MAT_DGP,NMIS_FORM_MAT_DGP,GET_NEW_MAT,AllOCATED(MATFACE)
    
       IF(GET_NEW_MAT) THEN
          IF(AllOCATED(MATFACE)) DEALLOCATE(MATFACE)
          IF(ALLOCATED(ILOCMATFACE)) DEALLOCATE(ILOCMATFACE)
          IF(ALLOCATED(JLOC2MATFACE)) DEALLOCATE(JLOC2MATFACE)
          IF(ALLOCATED(MATELE)) DEALLOCATE(MATELE)
      
          ALLOCATE(MATFACE(TOTELE,4,SMLOC,SMLOC))
          ewrite(2,*) 'ILOCMATFACE mem:',TOTELE*4*sMLOC
          ALLOCATE(ILOCMATFACE(TOTELE,4,SMLOC))
          ewrite(2,*) 'JLOC2MATFACE mem:',TOTELE*4*sMLOC
          ALLOCATE(JLOC2MATFACE(TOTELE,4,SMLOC))
          ewrite(2,*) 'matele mem:',TOTELE*MLOC*MLOC
          ALLOCATE(MATELE(TOTELE,MLOC,MLOC))
          MATELE = 0.0
          MATFACE = 0.0
       ENDIF

       ! Calculate part of rhs vec RHSVEC***********************************.
       RHSVEC(1:TOTELE*MLOC) = 0.0

       ! FORM THE MATRIX TO BE SOLVED (MATGLO,MATFACE,RHSVEC) FOR DGP **
       second_element_loop: DO ELE=1,TOTELE

          ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
          CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MLOC,NGI, &
               N,NLX,NLY,NLZ, M,MLX,MLY,MLZ, WEIGHT, DETWEI,.TRUE.,.FALSE., &
               NX,NY,NZ, MX,MY,MZ, &
               A11,A12,A13, A21,A22,A23, A31,A32,A33, &
               XD,YD,ZD, &
               ISPHERE)

          DO GI=1,NGI
             OME2(GI)=2.*FUNOME(XD(GI),YD(GI),ZD(GI))
             IF((GEOBAL.EQ.-12).or.(GEOBAL.EQ.-22)) OME2(GI)=0.
             IF(IDGBOYONLY.EQ.1) OME2(GI)=0.

             GIPGRAVX(GI)=0.0
             GIPGRAVY(GI)=0.0
             GIPGRAVZ(GI)=0.0
             DIRX(GI)=0.0
             DIRY(GI)=0.0
             DIRZ(GI)=0.0
             FSXGI(GI)=0.0
             FSYGI(GI)=0.0
             FSZGI(GI)=0.0
             SOUXGI(GI)=0.0
             SOUYGI(GI)=0.0
             SOUZGI(GI)=0.0
             UD(GI)=0.0
             VD(GI)=0.0
             WD(GI)=0.0
             RX(GI)=0.0
             RY(GI)=0.0
             RZ(GI)=0.0
             DO ILOC=1,MLOC
                IGL=PGNDGLNO((ELE-1)*MLOC+ILOC)
                IF(USEPOTGRA.EQ.1) THEN
                   GIPGRAVX(GI)=GIPGRAVX(GI) + MX(ILOC,GI)*POTGRA(IGL)
                   GIPGRAVY(GI)=GIPGRAVY(GI) + MY(ILOC,GI)*POTGRA(IGL)
                   GIPGRAVZ(GI)=GIPGRAVZ(GI) + MZ(ILOC,GI)*POTGRA(IGL)
                ELSE
                   GIPGRAVX(GI)=GIPGRAVX(GI) + M(ILOC,GI)*PGRAVX(IGL)
                   GIPGRAVY(GI)=GIPGRAVY(GI) + M(ILOC,GI)*PGRAVY(IGL)
                   GIPGRAVZ(GI)=GIPGRAVZ(GI) + M(ILOC,GI)*PGRAVZ(IGL)
                ENDIF
                DIRX(GI)=DIRX(GI) - M(ILOC,GI)*PNORMX(IGL)
                DIRY(GI)=DIRY(GI) - M(ILOC,GI)*PNORMY(IGL)
                DIRZ(GI)=DIRZ(GI) - M(ILOC,GI)*PNORMZ(IGL)
                IF(GETFREES) THEN
                   FSXGI(GI)=FSXGI(GI) + MX(ILOC,GI)*(FREES(IGL)-EQETID(IGL))*GRAVTY*wdgravity(ele)
                   FSYGI(GI)=FSYGI(GI) + MY(ILOC,GI)*(FREES(IGL)-EQETID(IGL))*GRAVTY*wdgravity(ele)
                   FSZGI(GI)=FSZGI(GI) + MZ(ILOC,GI)*(FREES(IGL)-EQETID(IGL))*GRAVTY*wdgravity(ele)
                ENDIF
             END DO

             RN=SQRT(DIRX(gi)**2+DIRY(gi)**2+DIRZ(GI)**2)
             DIRX(GI)=DIRX(GI)/RN
             DIRY(GI)=DIRY(GI)/RN
             DIRZ(GI)=DIRZ(GI)/RN
             DO ILOC=1,NLOC
                IGL=NDGLNO((ELE-1)*NLOC+ILOC)
                SOUXGI(GI)=SOUXGI(GI) + N(ILOC,GI)*GIPGRAVX(GI)*NDENST(IGL)
                SOUYGI(GI)=SOUYGI(GI) + N(ILOC,GI)*GIPGRAVY(GI)*NDENST(IGL)
                SOUZGI(GI)=SOUZGI(GI) + N(ILOC,GI)*GIPGRAVZ(GI)*NDENST(IGL)
                UD(GI)=UD(GI) + N(ILOC,GI)*U(IGL)
                VD(GI)=VD(GI) + N(ILOC,GI)*V(IGL)
                WD(GI)=WD(GI) + N(ILOC,GI)*W(IGL)
                IF((INUSUB.NE.0).AND.(NSUBNVLOC.NE.0)) THEN
                   IDGNOD=(ELE-1)*NLOC+ILOC
                   UD(GI)=UD(GI) + N(ILOC,GI)*NUSUB(IDGNOD)
                   VD(GI)=VD(GI) + N(ILOC,GI)*NVSUB(IDGNOD)
                   WD(GI)=WD(GI) + N(ILOC,GI)*NWSUB(IDGNOD)
                ENDIF
             END DO

             IF(NSOGRASUB.NE.0) THEN
               
               IF(COGRAX) THEN
                 RDIRX=0.0
                 RDIRY=0.0
                 RDIRZ=1.0
               ELSE
                 RDIRX=XD(GI)
                 RDIRY=YD(GI)
                 RDIRZ=ZD(GI)
                 RN=SQRT(RDIRX**2+RDIRY**2+RDIRZ**2)
                 RDIRX=RDIRX/RN 
                 RDIRY=RDIRY/RN
                 RDIRZ=RDIRZ/RN
               ENDIF
               
               RRSX=0.0
               RRSY=0.0
               RRSZ=0.0
               DO ILOC=1,NLOC
                 RRSX=RRSX+N(ILOC,GI)*SOXGRASUB((ELE-1)*NLOC+ILOC)
                 RRSY=RRSY+N(ILOC,GI)*SOYGRASUB((ELE-1)*NLOC+ILOC)
                 RRSZ=RRSZ+N(ILOC,GI)*SOZGRASUB((ELE-1)*NLOC+ILOC)
               END DO
               RNRRS=SQRT(RRSX**2+RRSY**2+RRSZ**2)
               GRA_SIG_CHAN=SIGN(1.0,RDIRX*RRSX+RDIRY*RRSY+RDIRZ*RRSZ)
               
               SOUXGI(GI)=SOUXGI(GI) + RDIRX*GRA_SIG_CHAN*RNRRS
               SOUYGI(GI)=SOUYGI(GI) + RDIRY*GRA_SIG_CHAN*RNRRS
               SOUZGI(GI)=SOUZGI(GI) + RDIRZ*GRA_SIG_CHAN*RNRRS
             ENDIF

             RX(GI)=(-FSXGI(GI)+VD(GI)*OME2(GI)+ SOUXGI(GI) )
             RY(GI)=(-FSYGI(GI)-UD(GI)*OME2(GI)+ SOUYGI(GI) )
             RZ(GI)=(-FSZGI(GI) + SOUZGI(GI) )

          END DO

          DO ILOC=1,MLOC
             RNN=0.0
             DO GI=1,NGI
                RNN=RNN &
                     +M(ILOC,GI)*(DIRX(GI)*RX(GI)  &
                     +DIRY(GI)*RY(GI)+DIRZ(GI)*RZ(GI))*DETWEI(GI)
             END DO
             RHSVEC((ELE-1)*MLOC+ILOC)=RHSVEC((ELE-1)*MLOC+ILOC)  +RNN
          END DO

        IF(GET_NEW_MAT) THEN
          DO ILOC=1,MLOC
             DO JLOC=1,MLOC
                RNN=0.0
                DO GI=1,NGI
                   RNN=RNN &
                        +M(ILOC,GI)*(DIRX(GI)*MX(JLOC,GI)+DIRY(GI)*MY(JLOC,GI)  &
                        +DIRZ(GI)*MZ(JLOC,GI) )*DETWEI(GI)
                END DO
                MATELE(ELE,ILOC,JLOC)=MATELE(ELE,ILOC,JLOC) +RNN
             END DO
          END DO

          ! Put surface integrals into matrix and rhs*****************************
          second_face_loop: DO IFACE=1,4
             COUNT2=(ELE-1)*5+IFACE
             ! NB colele4 is ordered in terms of faces.
             ELE2=COLELE4(COUNT2)

             ! Surface element of domain...
             DO SILOC=1,SNLOC
                SILOC2ILOC(SILOC)=LOCLIST(IFACE,SILOC)
             END DO
             ! Form OTHILOC from SILOC2ILOC
             IF(ELE2.NE.0) THEN
                CALL ELE2ELEFACE(OTHILOC,SILOC2ILOC,SNLOC,NLOC,TOTELE, &
                     XONDGL,NDGLNO,ELE,ELE2)
             ELSE
                DO ILOC=1,NLOC
                   OTHILOC(ILOC)=ILOC
                END DO
             ENDIF

             ! Form MLOCOTHILOC and MLOCSILOC2ILOC
             CALL MLOCELE2ELEFACE(MLOCOTHILOC,MLOCSILOC2ILOC,SNLOC,NLOC, &
                  SMLOC,MLOC,OTHILOC,SILOC2ILOC)

             ! PERFORM SURFACE INTEGRATION         !
             ! Form approximate surface normal (NORMX,NORMY,NORMZ)
             CALL DGSIMPLNORM(ELE,SILOC2ILOC,TOTELE,NLOC,SNLOC,XONDGL, &
                  X,Y,Z,XNONOD,NORMX,NORMY,NORMZ)

             CALL DGSDETNXLOC(ELE,SILOC2ILOC,XONDGL, &
                  TOTELE,NLOC,XNONOD,SNLOC,SNGI, &
                  X,Y,Z, &
                  SN,SNLX,SNLY, SWEIGH, SDETWE,SAREA, .TRUE.,.FALSE., &
                  SNORMXN,SNORMYN,SNORMZN, &
                  NORMX,NORMY,NORMZ, &
                  ISPHERE,SMLOC,SM,SMLX,SMLY)

             TOPSUF=.FALSE.
             DO SGI=1,SNGI
                SDIRX(SGI)=0.0
                SDIRY(SGI)=0.0
                SDIRZ(SGI)=0.0
                DO SILOC=1,SMLOC
                   ILOC=MLOCSILOC2ILOC(SILOC)
                   IGL=PGNDGLNO((ELE-1)*MLOC+ILOC)
                   SDIRX(SGI)=SDIRX(SGI)-SM(SILOC,SGI)*PNORMX(IGL)
                   SDIRY(SGI)=SDIRY(SGI)-SM(SILOC,SGI)*PNORMY(IGL)
                   SDIRZ(SGI)=SDIRZ(SGI)-SM(SILOC,SGI)*PNORMZ(IGL)
                END DO
                RN=SQRT(SDIRX(SGI)**2+SDIRY(SGI)**2+SDIRZ(SGI)**2)
                SDIRX(SGI)=SDIRX(SGI)/RN
                SDIRY(SGI)=SDIRY(SGI)/RN
                SDIRZ(SGI)=SDIRZ(SGI)/RN
                NDOTQ(SGI)=SNORMXN(SGI)*SDIRX(SGI)+SNORMYN(SGI)*SDIRY(SGI) &
                     +SNORMZN(SGI)*SDIRZ(SGI)
                INCOME(SGI)=0.0
                IF(NDOTQ(SGI).LT.0.0) INCOME(SGI)=1.0
                IF(ELE2.EQ.0) THEN
                   IF(NDOTQ(SGI).LT.-0.8) TOPSUF=.TRUE.
                ENDIF
             END DO

             ! Perform surface integration...
             DO SILOC=1,SMLOC
                ILOC=MLOCSILOC2ILOC(SILOC)
                GLOBI=(ELE-1)*MLOC+ILOC
                ILOC2=MLOCOTHILOC(ILOC)
                GLOBI2=(ELE2-1)*MLOC+ILOC2
                DO SJLOC=1,SMLOC
                   JLOC=MLOCSILOC2ILOC(SJLOC)
                   GLOBJ=(ELE-1)*MLOC+JLOC
                   JLOC2=MLOCOTHILOC(JLOC)
                   GLOBJ2=(ELE2-1)*MLOC+JLOC2
                   INVLM =0.0
                   DO SGI=1,SNGI
                      INVLM= INVLM +INCOME(SGI)*NDOTQ(SGI)*SM(SILOC,SGI)*SM(SJLOC,SGI)*SDETWE(SGI)
                   END DO
                   MATELE(ELE,ILOC,JLOC)=MATELE(ELE,ILOC,JLOC)-INVLM
                   IF(TOPSUF) THEN
                      MATFACE(ELE,IFACE,SILOC,SJLOC)=MATFACE(ELE,IFACE,SILOC,SJLOC)-0.0
                   ELSE IF(ELE2.EQ.0) THEN
                      MATELE(ELE,ILOC,JLOC)=MATELE(ELE,ILOC,JLOC)+INVLM
                   ELSE
                      MATFACE(ELE,IFACE,SILOC,SJLOC)=MATFACE(ELE,IFACE,SILOC,SJLOC)-INVLM
                   ENDIF
                   ILOCMATFACE(ELE,IFACE,SILOC) =ILOC
                   JLOC2MATFACE(ELE,IFACE,SJLOC)=JLOC2

                END DO
             END DO

          END DO second_face_loop
! END OF IF(GET_MAT) THEN          
        ENDIF

       END DO second_element_loop
       ewrite(2,*) 'r2norm(rhsvec,mloc*totele,0):',r2norm(rhsvec,mloc*totele,0)

       ! Only 1 iteration is needed for a structured mesh...
       NITS=30
       IF((NLEVEL.GT.1).OR.(one_vdg_ser_it.EQ.1)) NITS=1
       interation_loop: DO ITS=1,NITS
          ewrite(2,*) 'its=',its
!          print *,'***********************************'
!          print *,'***********************************'
!          print *,'***vertical dg its=',its
!          print *,'***********************************'
!          print *,'***********************************' 
          MAXDGP=0.0
          MAXDIFF=0.0
          DO II=1,TOTELE
             ELE=ELEORD(II)
             DO ILOC=1,MLOC
                DO JLOC=1,MLOC
                   MAT(ILOC,JLOC)=MATELE(ELE,ILOC,JLOC)
                END DO
             END DO
             VECLOC(1:MLOC) = RHSVEC((ELE-1)*MLOC+1:(ELE-1)*MLOC+1 + MLOC - 1)

             DO IFACE=1,4
                COUNT2=(ELE-1)*5+IFACE
                ! NB colele4 is ordered in terms of faces.
                ELE2=COLELE4(COUNT2)
                ELE3=ELE2
                IF(ELE2.EQ.0) ELE3=ELE
                DO SILOC=1,SMLOC
                   ILOC=ILOCMATFACE(ELE,IFACE,SILOC)
                   DO SJLOC=1,SMLOC
                      JLOC2=JLOC2MATFACE(ELE,IFACE,SJLOC)
                      GLOBJ2=(ELE3-1)*MLOC+JLOC2
                      VECLOC(ILOC)=VECLOC(ILOC)+MATFACE(ELE,IFACE,SILOC,SJLOC)*DGP(GLOBJ2)
                   END DO
                END DO
             END DO

             ! Solve MAT DGP=VECLOC (Store the decomposition and use this after the 1st iteration)
             IF((ITS.EQ.1).AND.(GET_NEW_MAT)) THEN
                CALL SMLINNGOT(MAT,LOCDGP,VECLOC,MLOC,MLOC,.FALSE.)
                ! Store the decomposition...
                DO ILOC=1,MLOC
                   DO JLOC=1,MLOC
                      MATELE(ELE,ILOC,JLOC)=MAT(ILOC,JLOC)
                   END DO
                END DO
             ELSE
                CALL SMLINNGOT(MAT,LOCDGP,VECLOC,MLOC,MLOC,.TRUE.)
             ENDIF

             DO ILOC=1,MLOC
                MAXDIFF=MAX(ABS(DGP((ELE-1)*MLOC+ILOC)-LOCDGP(ILOC)),MAXDIFF)
                MAXDGP=MAX(ABS(LOCDGP(ILOC)),MAXDGP)
                DGP((ELE-1)*MLOC+ILOC)=LOCDGP(ILOC)
             END DO

          END DO
          ewrite(2,*) 'maxdiff,maxdgp,TOLERCENT*MAXDGP:',maxdiff,maxdgp,TOLERCENT*MAXDGP
          IF(MAXDIFF.LE.TOLERCENT*MAXDGP) exit
       END DO interation_loop

       ! ENDOF IF(PHSOLVE) THEN...      
    ENDIF

    ewrite(2,*) 'passed ok'
    CALL PMINMX(DGP,TOTELE*MLOC,'******DGP  ')

    ewrite(2,*) 'r2norm(dgp,mloc*totele,0):',r2norm(dgp,mloc*totele,0)
    ewrite(2,*) 'r2norm(SOXGRASUB,NSOGRASUB*totele,0):',r2norm(SOXGRASUB,NSOGRASUB*totele,0)
    ewrite(2,*) 'r2norm(SOyGRASUB,NSOGRASUB*totele,0):',r2norm(SOyGRASUB,NSOGRASUB*totele,0)
    ewrite(2,*) 'r2norm(SOzGRASUB,NSOGRASUB*totele,0):',r2norm(SOzGRASUB,NSOGRASUB*totele,0)

    ewrite(1,*) '******* 2-NOW PUT SOLUTION PSI into VECX,VECY,VECZ'
    ! NOW PUT SOLUTION PSI into VECX,VECY,VECZ:

    third_element_loop: DO ELE=1,TOTELE

       ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
       CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MLOC,NGI, &
            N,NLX,NLY,NLZ, M,MLX,MLY,MLZ, WEIGHT, DETWEI,.TRUE.,.FALSE., &
            NX,NY,NZ, MX,MY,MZ, &
            A11,A12,A13, A21,A22,A23, A31,A32,A33, &
            XD,YD,ZD, &
            ISPHERE)

       DO GI=1,NGI
          OME2(GI)=2.*FUNOME(XD(GI),YD(GI),ZD(GI))
          IF((GEOBAL.EQ.-12).or.(GEOBAL.EQ.-22)) OME2(GI)=0.
          IF(IDGBOYONLY.EQ.1) OME2(GI)=0.

          GIPGRAVX(GI)=0.0
          GIPGRAVY(GI)=0.0
          GIPGRAVZ(GI)=0.0
          DIRX(GI)=0.0
          DIRY(GI)=0.0
          DIRZ(GI)=0.0
          FSXGI(GI)=0.0
          FSYGI(GI)=0.0
          FSZGI(GI)=0.0
          FSXGIC(GI)=0.0
          FSYGIC(GI)=0.0
          FSZGIC(GI)=0.0
          SOUXGI(GI)=0.0
          SOUYGI(GI)=0.0
          SOUZGI(GI)=0.0
          UD(GI)=0.0
          VD(GI)=0.0
          WD(GI)=0.0
          RX(GI)=0.0
          RY(GI)=0.0
          RZ(GI)=0.0
          RXC(GI)=0.0
          RYC(GI)=0.0
          RZC(GI)=0.0
          DO ILOC=1,MLOC
             IGL=PGNDGLNO((ELE-1)*MLOC+ILOC)
             IF(USEPOTGRA.EQ.1) THEN
                GIPGRAVX(GI)=GIPGRAVX(GI) + MX(ILOC,GI)*POTGRA(IGL)
                GIPGRAVY(GI)=GIPGRAVY(GI) + MY(ILOC,GI)*POTGRA(IGL)
                GIPGRAVZ(GI)=GIPGRAVZ(GI) + MZ(ILOC,GI)*POTGRA(IGL)
             ELSE
                GIPGRAVX(GI)=GIPGRAVX(GI) + M(ILOC,GI)*PGRAVX(IGL)
                GIPGRAVY(GI)=GIPGRAVY(GI) + M(ILOC,GI)*PGRAVY(IGL)
                GIPGRAVZ(GI)=GIPGRAVZ(GI) + M(ILOC,GI)*PGRAVZ(IGL)
             ENDIF
             DIRX(GI)=DIRX(GI) - M(ILOC,GI)*PNORMX(IGL)
             DIRY(GI)=DIRY(GI) - M(ILOC,GI)*PNORMY(IGL)
             DIRZ(GI)=DIRZ(GI) - M(ILOC,GI)*PNORMZ(IGL)
             IF(GETFREES) THEN
                FSXGI(GI)=FSXGI(GI) + MX(ILOC,GI)*(FREES(IGL)-EQETID(IGL))*GRAVTY*wdgravity(ele)
                FSYGI(GI)=FSYGI(GI) + MY(ILOC,GI)*(FREES(IGL)-EQETID(IGL))*GRAVTY*wdgravity(ele)
                FSZGI(GI)=FSZGI(GI) + MZ(ILOC,GI)*(FREES(IGL)-EQETID(IGL))*GRAVTY*wdgravity(ele)
             ENDIF
          END DO

          RN=SQRT(DIRX(gi)**2+DIRY(gi)**2+DIRZ(GI)**2)
          DIRX(GI)=DIRX(GI)/RN
          DIRY(GI)=DIRY(GI)/RN
          DIRZ(GI)=DIRZ(GI)/RN
          !  The rotation matrix in 3-D is R=
          !   T1X   T1Y   T1Z
          !   T2X   T2Y   T2Z
          !   NORMX NORMY NORMZ
          ! Calculate T1X,T1Y,T1Z, T2X,T2Y,T2Z from NORMX,NORMY,NORMZ
          if(MOMBC) then

             RTREST11(GI)=1.0 - DIRX(GI)*DIRX(GI)
             RTREST12(GI)=    - DIRX(GI)*DIRY(GI)
             RTREST13(GI)=    - DIRX(GI)*DIRZ(GI)

             RTREST21(GI)=    - DIRY(GI)*DIRX(GI)
             RTREST22(GI)=1.0 - DIRY(GI)*DIRY(GI)
             RTREST23(GI)=    - DIRY(GI)*DIRZ(GI)

             RTREST31(GI)=    - DIRZ(GI)*DIRX(GI)
             RTREST32(GI)=    - DIRZ(GI)*DIRY(GI)
             RTREST33(GI)=1.0 - DIRZ(GI)*DIRZ(GI)
          else
             RTREST11(GI)=1.0
             RTREST12(GI)=0.0
             RTREST13(GI)=0.0

             RTREST21(GI)=0.0
             RTREST22(GI)=1.0
             RTREST23(GI)=0.0

             RTREST31(GI)=0.0
             RTREST32(GI)=0.0
             RTREST33(GI)=1.0
          endif

          DO ILOC=1,NLOC
             IGL=NDGLNO((ELE-1)*NLOC+ILOC)
             SOUXGI(GI)=SOUXGI(GI) + N(ILOC,GI)*GIPGRAVX(GI)*NDENST(IGL)
             SOUYGI(GI)=SOUYGI(GI) + N(ILOC,GI)*GIPGRAVY(GI)*NDENST(IGL)
             SOUZGI(GI)=SOUZGI(GI) + N(ILOC,GI)*GIPGRAVZ(GI)*NDENST(IGL)
             UD(GI)=UD(GI) + N(ILOC,GI)*U(IGL)
             VD(GI)=VD(GI) + N(ILOC,GI)*V(IGL)
             WD(GI)=WD(GI) + N(ILOC,GI)*W(IGL)
             IF((INUSUB.NE.0).AND.(NSUBNVLOC.NE.0)) THEN
                IDGNOD=(ELE-1)*NLOC+ILOC
                UD(GI)=UD(GI) + N(ILOC,GI)*NUSUB(IDGNOD)
                VD(GI)=VD(GI) + N(ILOC,GI)*NVSUB(IDGNOD)
                WD(GI)=WD(GI) + N(ILOC,GI)*NWSUB(IDGNOD)
             ENDIF
          END DO
          
             IF(NSOGRASUB.NE.0) THEN
               IF(COGRAX) THEN
                 RDIRX=0.0
                 RDIRY=0.0
                 RDIRZ=1.0
               ELSE
                 RDIRX=XD(GI)
                 RDIRY=YD(GI)
                 RDIRZ=ZD(GI)
                 RN=SQRT(RDIRX**2+RDIRY**2+RDIRZ**2)
                 RDIRX=RDIRX/RN 
                 RDIRY=RDIRY/RN
                 RDIRZ=RDIRZ/RN
               ENDIF
               
               RRSX=0.0
               RRSY=0.0
               RRSZ=0.0
               DO ILOC=1,NLOC
                 RRSX=RRSX+N(ILOC,GI)*SOXGRASUB((ELE-1)*NLOC+ILOC)
                 RRSY=RRSY+N(ILOC,GI)*SOYGRASUB((ELE-1)*NLOC+ILOC)
                 RRSZ=RRSZ+N(ILOC,GI)*SOZGRASUB((ELE-1)*NLOC+ILOC)
               END DO
               RNRRS=SQRT(RRSX**2+RRSY**2+RRSZ**2)
               GRA_SIG_CHAN=SIGN(1.0,RDIRX*RRSX+RDIRY*RRSY+RDIRZ*RRSZ)
               
               SOUXGI(GI)=SOUXGI(GI) + RDIRX*GRA_SIG_CHAN*RNRRS
               SOUYGI(GI)=SOUYGI(GI) + RDIRY*GRA_SIG_CHAN*RNRRS
               SOUZGI(GI)=SOUZGI(GI) + RDIRZ*GRA_SIG_CHAN*RNRRS
             ENDIF

          RX(GI)=-FSXGI(GI)+VD(GI)*OME2(GI)+ SOUXGI(GI)
          RY(GI)=-FSYGI(GI)-UD(GI)*OME2(GI)+ SOUYGI(GI)
          RZ(GI)=-FSZGI(GI) + SOUZGI(GI)
          RXR(GI)=RTREST11(GI)*RX(GI)+RTREST12(GI)*RY(GI)+RTREST13(GI)*RZ(GI)
          RYR(GI)=RTREST21(GI)*RX(GI)+RTREST22(GI)*RY(GI)+RTREST23(GI)*RZ(GI)
          RZR(GI)=RTREST31(GI)*RX(GI)+RTREST32(GI)*RY(GI)+RTREST33(GI)*RZ(GI)
          
          RXC(GI)=(-FSXGIC(GI)+VD(GI)*OME2(GI) )
          RYC(GI)=(-FSYGIC(GI)-UD(GI)*OME2(GI) )
          RZC(GI)=(-FSZGIC(GI)  )
          RXRC(GI)=RTREST11(GI)*RXC(GI)+RTREST12(GI)*RYC(GI)+RTREST13(GI)*RZC(GI)
          RYRC(GI)=RTREST21(GI)*RXC(GI)+RTREST22(GI)*RYC(GI)+RTREST23(GI)*RZC(GI)
          RZRC(GI)=RTREST31(GI)*RXC(GI)+RTREST32(GI)*RYC(GI)+RTREST33(GI)*RZC(GI)
       END DO

       DO ILOC=1,NLOC
          GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
          DO JLOC=1,MLOC
             GLOBJ=(ELE-1)*MLOC+JLOC

             RRVECX=0.0
             RRVECY=0.0
             RRVECZ=0.0

             RNN=0.0
             DO GI=1,NGI

                RRVECX=RRVECX-N(ILOC,GI)*DETWEI(GI)*DGP(GLOBJ) &
                     *(RTREST11(GI)*MX(JLOC,GI)+RTREST12(GI)*MY(JLOC,GI)+RTREST13(GI)*MZ(JLOC,GI))
                RRVECY=RRVECY-N(ILOC,GI)*DETWEI(GI)*DGP(GLOBJ) &
                     *(RTREST21(GI)*MX(JLOC,GI)+RTREST22(GI)*MY(JLOC,GI)+RTREST23(GI)*MZ(JLOC,GI))
                RRVECZ=RRVECZ-N(ILOC,GI)*DETWEI(GI)*DGP(GLOBJ) &
                     *(RTREST31(GI)*MX(JLOC,GI)+RTREST32(GI)*MY(JLOC,GI)+RTREST33(GI)*MZ(JLOC,GI))

                RNN=RNN+N(ILOC,GI)*DETWEI(GI)*M(JLOC,GI)

             END DO
             VECX(GLOBI)=VECX(GLOBI) +RRVECX
             VECY(GLOBI)=VECY(GLOBI) +RRVECY
             VECZ(GLOBI)=VECZ(GLOBI) +RRVECZ

             VECX2(GLOBI)=VECX2(GLOBI) +RRVECX
             VECY2(GLOBI)=VECY2(GLOBI) +RRVECY
             VECZ2(GLOBI)=VECZ2(GLOBI) +RRVECZ

             IF(NSUBVLOC.NE.0) THEN
                IDGNOD=(ELE-1)*NSUBVLOC+ILOC
                VECX_SGSADD(IDGNOD)=VECX_SGSADD(IDGNOD)+RRVECX
                VECY_SGSADD(IDGNOD)=VECY_SGSADD(IDGNOD)+RRVECY
                VECZ_SGSADD(IDGNOD)=VECZ_SGSADD(IDGNOD)+RRVECZ

                VECX_SGSADD2(IDGNOD)=VECX_SGSADD2(IDGNOD)+RRVECX
                VECY_SGSADD2(IDGNOD)=VECY_SGSADD2(IDGNOD)+RRVECY
                VECZ_SGSADD2(IDGNOD)=VECZ_SGSADD2(IDGNOD)+RRVECZ
                MLELE(IDGNOD)=MLELE(IDGNOD)+RNN
             ENDIF


             ML(GLOBI)=ML(GLOBI)+RNN

          END DO

          DO GI=1,NGI
             VECX(GLOBI)=VECX(GLOBI)+N(ILOC,GI)*DETWEI(GI)*RXR(GI)
             VECY(GLOBI)=VECY(GLOBI)+N(ILOC,GI)*DETWEI(GI)*RYR(GI)
             VECZ(GLOBI)=VECZ(GLOBI)+N(ILOC,GI)*DETWEI(GI)*RZR(GI)

             VECXF(GLOBI)=VECXF(GLOBI)+N(ILOC,GI)*DETWEI(GI)*RXRC(GI)
             VECYF(GLOBI)=VECYF(GLOBI)+N(ILOC,GI)*DETWEI(GI)*RYRC(GI)
             VECZF(GLOBI)=VECZF(GLOBI)+N(ILOC,GI)*DETWEI(GI)*RZRC(GI)
             IF(NSUBVLOC.NE.0) THEN
                IDGNOD=(ELE-1)*NSUBVLOC+ILOC
                VECX_SGSADD(IDGNOD)=VECX_SGSADD(IDGNOD)+N(ILOC,GI)*DETWEI(GI)*RXR(GI)
                VECY_SGSADD(IDGNOD)=VECY_SGSADD(IDGNOD)+N(ILOC,GI)*DETWEI(GI)*RYR(GI)
                VECZ_SGSADD(IDGNOD)=VECZ_SGSADD(IDGNOD)+N(ILOC,GI)*DETWEI(GI)*RZR(GI)
             ENDIF
          END DO
       END DO

       ! Put surface integrals into matrix and rhs*****************************
       if(MOMBC) then
          third_face_loop: DO IFACE=1,4
             COUNT2=(ELE-1)*5+IFACE
             ! NB colele4 is ordered in terms of faces.
             ELE2=COLELE4(COUNT2)

             ! Surface element of domain...
             DO SILOC=1,SNLOC
                SILOC2ILOC(SILOC)=LOCLIST(IFACE,SILOC)
             END DO
             ! Form OTHILOC from SILOC2ILOC
             IF(ELE2.NE.0) THEN
                CALL ELE2ELEFACE(OTHILOC,SILOC2ILOC,SNLOC,NLOC,TOTELE, &
                     XONDGL,NDGLNO,ELE,ELE2)
             ELSE
                DO ILOC=1,NLOC
                   OTHILOC(ILOC)=ILOC
                END DO
             ENDIF

             ! Form MLOCOTHILOC and MLOCSILOC2ILOC
             CALL MLOCELE2ELEFACE(MLOCOTHILOC,MLOCSILOC2ILOC,SNLOC,NLOC, &
                  SMLOC,MLOC,OTHILOC,SILOC2ILOC)

             ! PERFORM SURFACE INTEGRATION         !
             ! Form approximate surface normal (NORMX,NORMY,NORMZ)
             CALL DGSIMPLNORM(ELE,SILOC2ILOC,TOTELE,NLOC,SNLOC,XONDGL, &
                  X,Y,Z,XNONOD,NORMX,NORMY,NORMZ)

             CALL DGSDETNXLOC(ELE,SILOC2ILOC,XONDGL, &
                  TOTELE,NLOC,XNONOD,SNLOC,SNGI, &
                  X,Y,Z, &
                  SN,SNLX,SNLY, SWEIGH, SDETWE,SAREA, .TRUE.,.FALSE., &
                  SNORMXN,SNORMYN,SNORMZN, &
                  NORMX,NORMY,NORMZ, &
                  ISPHERE,SMLOC,SM,SMLX,SMLY)

             TOPSUF=.FALSE.
             DO SGI=1,SNGI
                SDIRX(SGI)=0.0
                SDIRY(SGI)=0.0
                SDIRZ(SGI)=0.0
                DO SILOC=1,SMLOC
                   ILOC=MLOCSILOC2ILOC(SILOC)
                   IGL=PGNDGLNO((ELE-1)*MLOC+ILOC)
                   SDIRX(SGI)=SDIRX(SGI)-SM(SILOC,SGI)*PNORMX(IGL)
                   SDIRY(SGI)=SDIRY(SGI)-SM(SILOC,SGI)*PNORMY(IGL)
                   SDIRZ(SGI)=SDIRZ(SGI)-SM(SILOC,SGI)*PNORMZ(IGL)
                END DO
                RN=SQRT(SDIRX(SGI)**2+SDIRY(SGI)**2+SDIRZ(SGI)**2)
                SDIRX(SGI)=SDIRX(SGI)/RN
                SDIRY(SGI)=SDIRY(SGI)/RN
                SDIRZ(SGI)=SDIRZ(SGI)/RN

                SRTVX(SGI)=SDIRX(SGI)*SDIRX(SGI)*SNORMXN(SGI) &
                     +SDIRX(SGI)*SDIRY(SGI)*SNORMYN(SGI) &
                     +SDIRX(SGI)*SDIRZ(SGI)*SNORMZN(SGI)
                SRTVY(SGI)=SDIRY(SGI)*SDIRX(SGI)*SNORMXN(SGI) &
                     +SDIRY(SGI)*SDIRY(SGI)*SNORMYN(SGI) &
                     +SDIRY(SGI)*SDIRZ(SGI)*SNORMZN(SGI)
                SRTVZ(SGI)=SDIRZ(SGI)*SDIRX(SGI)*SNORMXN(SGI) &
                     +SDIRZ(SGI)*SDIRY(SGI)*SNORMYN(SGI) &
                     +SDIRZ(SGI)*SDIRZ(SGI)*SNORMZN(SGI)
                SRESTX(SGI)=SNORMXN(SGI) - SRTVX(SGI)
                SRESTY(SGI)=SNORMYN(SGI) - SRTVY(SGI)
                SRESTZ(SGI)=SNORMZN(SGI) - SRTVZ(SGI)
                NDOTQ(SGI)=SNORMXN(SGI)*SDIRX(SGI)+SNORMYN(SGI)*SDIRY(SGI) &
                     +SNORMZN(SGI)*SDIRZ(SGI)
                INCOME(SGI)=0.0
                IF(NDOTQ(SGI).LT.0.0) INCOME(SGI)=1.0
                IF(ELE2.EQ.0) THEN
                   IF(NDOTQ(SGI).LT.-0.8) TOPSUF=.TRUE.
                ENDIF
             END DO

             ! Perform surface integration...
             siloc_loop: DO SILOC=1,SNLOC
                ILOC=SILOC2ILOC(SILOC)
                GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
                sjloc_loop: DO SJLOC=1,SMLOC
                   JLOC=MLOCSILOC2ILOC(SJLOC)
                   GLOBJ=(ELE-1)*MLOC+JLOC
                   JLOC2=MLOCOTHILOC(JLOC)
                   GLOBJ2=(ELE2-1)*MLOC+JLOC2

                   INVLMX=0.0
                   INVLMY=0.0
                   INVLMZ=0.0
                   BOTHVLMX=0.0
                   BOTHVLMY=0.0
                   BOTHVLMZ=0.0
                   DO SGI=1,SNGI
                      RNN=SN(SILOC,SGI)*SM(SJLOC,SGI)*SDETWE(SGI)
                      BOTHVLMX=BOTHVLMX+SRESTX(SGI)*RNN
                      BOTHVLMY=BOTHVLMY+SRESTY(SGI)*RNN
                      BOTHVLMZ=BOTHVLMZ+SRESTZ(SGI)*RNN
                   END DO
                   ! Assume a zero inlet b.c when we have an inlet into the 
                   ! domain -no surface integral
                   ! Include surface eleement at bottom and sides of domain only...
                   IF(ELE2.EQ.0) THEN
                      POTH=DGP(GLOBJ)
                   ELSE
                      POTH=DGP(GLOBJ2) 
                   ENDIF
                   VECX(GLOBI)=VECX(GLOBI)+INVLMX*(DGP(GLOBJ)-POTH) &
                        +0.5*BOTHVLMX*(DGP(GLOBJ)-POTH)
                   VECY(GLOBI)=VECY(GLOBI)+INVLMY*(DGP(GLOBJ)-POTH) &
                        +0.5*BOTHVLMY*(DGP(GLOBJ)-POTH)
                   VECZ(GLOBI)=VECZ(GLOBI)+INVLMZ*(DGP(GLOBJ)-POTH) &
                        +0.5*BOTHVLMZ*(DGP(GLOBJ)-POTH)

                   VECX2(GLOBI)=VECX2(GLOBI)+INVLMX*(DGP(GLOBJ)-POTH) &
                        +0.5*BOTHVLMX*(DGP(GLOBJ)-POTH)
                   VECY2(GLOBI)=VECY2(GLOBI)+INVLMY*(DGP(GLOBJ)-POTH) &
                        +0.5*BOTHVLMY*(DGP(GLOBJ)-POTH)
                   VECZ2(GLOBI)=VECZ2(GLOBI)+INVLMZ*(DGP(GLOBJ)-POTH) &
                        +0.5*BOTHVLMZ*(DGP(GLOBJ)-POTH)
                   IF(NSUBVLOC.NE.0) THEN
                      IDGNOD=(ELE-1)*NSUBVLOC+ILOC
                      VECX_SGSADD(IDGNOD)=VECX_SGSADD(IDGNOD)+INVLMX*(DGP(GLOBJ)-POTH) &
                           +0.5*BOTHVLMX*(DGP(GLOBJ)-POTH)
                      VECY_SGSADD(IDGNOD)=VECY_SGSADD(IDGNOD)+INVLMY*(DGP(GLOBJ)-POTH) &
                           +0.5*BOTHVLMY*(DGP(GLOBJ)-POTH)
                      VECZ_SGSADD(IDGNOD)=VECZ_SGSADD(IDGNOD)+INVLMZ*(DGP(GLOBJ)-POTH) &
                           +0.5*BOTHVLMZ*(DGP(GLOBJ)-POTH)

                      VECX_SGSADD2(IDGNOD)=VECX_SGSADD2(IDGNOD)+INVLMX*(DGP(GLOBJ)-POTH) &
                           +0.5*BOTHVLMX*(DGP(GLOBJ)-POTH)
                      VECY_SGSADD2(IDGNOD)=VECY_SGSADD2(IDGNOD)+INVLMY*(DGP(GLOBJ)-POTH) &
                           +0.5*BOTHVLMY*(DGP(GLOBJ)-POTH)
                      VECZ_SGSADD2(IDGNOD)=VECZ_SGSADD2(IDGNOD)+INVLMZ*(DGP(GLOBJ)-POTH) &
                           +0.5*BOTHVLMZ*(DGP(GLOBJ)-POTH)
                   ENDIF
                END DO sjloc_loop
             END DO siloc_loop

          END DO third_face_loop
       ENDIF

    END DO third_element_loop

    ! Smooth out the resulting force...

    IF((INUSUB.NE.0).AND.(NSUBNVLOC.NE.0).AND.(NSOGRASUB.NE.0)) THEN
    IF(PSMOOTH.OR.FSMOOTH) THEN 
     
      CALL SMOOTH_FORCE(NONODS,NSUBVLOC,NLOC,TOTELE, &
        NDGLNO, PSMOOTH, FSMOOTH, DGSMOOTH, &
        VECX,VECY,VECZ, VECX2,VECY2,VECZ2, ML,MLELE, &
        VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD, &
        VECX_SGSADD2,VECY_SGSADD2,VECZ_SGSADD2, &
        VECXF,VECYF,VECZF, &
! For sub CALHIGHNODAINN...
        X,Y,Z, XONDGL, XNONOD,MLOC,NGI, &
        N,NLX,NLY,NLZ,MLX,MLY,MLZ, WEIGHT, DETWEI, &
        ISPHERE)
        
    ENDIF
    ENDIF  
    
    IF(NMIS_FORM_MAT_DGP.LE.1) THEN
!    IF(COUNT_MIS_FORM_MAT_DGP+1.GE.NMIS_FORM_MAT_DGP)  THEN
! Dont storing matrix or de-allocate ready for next visit to this sub. 
          IF(AllOCATED(MATFACE)) DEALLOCATE(MATFACE)
          IF(ALLOCATED(ILOCMATFACE)) DEALLOCATE(ILOCMATFACE)
          IF(ALLOCATED(JLOC2MATFACE)) DEALLOCATE(JLOC2MATFACE)
          IF(ALLOCATED(MATELE)) DEALLOCATE(MATELE)
!    ENDIF
    ENDIF

    ewrite(1,*) 'finished verticaldgstore'

  END SUBROUTINE VERTICALDGSTORE

  subroutine longalltopsufdg(longdifmath,longdifvecx, &
       xorig,yorig,zorig, &
       topdis,botdis,atheta,absorb, &
       gravty,ctheta, &
       nloc, &
       totele, &
       geobal,&
       ndglno,xondgl, &
       nonods,xnonod, &
       scfacth0,isphere, &
       pgfindrm,pgcentrm,pgcolm,npgcolm,pgnods, &
       cograx,&
       u,v,w, &
       INUSUB,NSUBNVLOC,NUSUB,NVSUB,NWSUB, &
       nofsta,velmean,dt, &
       mloc,pgndglno,frees,freesold, &
       fstheta, rnonl_coefficient)
    ! This subroutine caluculates free surface mass and
    ! surface momentum flux contributions by NOT integrating over the
    ! top surface of the ocean.
    ! M is the basis function assocaited with normal pressure.
    ! the 4th order switch switches on the lumping.
    use FLDebug
    use assnav_module
    IMPLICIT NONE
    ! Number of source nodes, advection nodes and node positions.
    ! Number of elements
    ! IF SMALLSTEN then use a small stencil other wise use a more accurate
    ! enlarged stencil.
    LOGICAL SLUMMAT,SLUMRHS,SMALLSTEN,NONLIN
    PARAMETER(SLUMMAT=.FALSE.,SLUMRHS=.FALSE.,SMALLSTEN=.TRUE.)
    ! If NONLIN then use Nonlinear dissipation.
    PARAMETER(NONLIN=.true.)
    ! IF SLUMMAT then lump the free surface height mass matrix.
    ! IF SLUMRHS then lump the rhs contribution to matrix eqn.
    INTEGER TOTELE,NONODS,XNONOD,PGNODS
    INTEGER NLOC

    REAL XORIG(XNONOD),YORIG(XNONOD),ZORIG(XNONOD)
    ! Global node number maps.
    INTEGER NDGLNO(TOTELE*NLOC)
    INTEGER XONDGL(TOTELE*NLOC)
    REAL GRAVTY,CTHETA,ATHETA
    REAL TOPDIS(XNONOD),BOTDIS(XNONOD)
    REAL ABSORB(NONODS)
    real scfacth0
    INTEGER GEOBAL
    ! This is the THETA time stepping parameter for the free surface.
    REAL FSTHETA
    REAL LONGDIFMATH(NPGCOLM)
    REAL LONGDIFVECX(PGNODS)
    INTEGER NPGCOLM,PGFINDRM(PGNODS+1),PGCOLM(NPGCOLM),PGCENTRM(PGNODS)
    LOGICAL COGRAX
    REAL U(NONODS),V(NONODS),W(NONODS)
    INTEGER INUSUB,NSUBNVLOC
    REAL NUSUB(TOTELE*NSUBNVLOC),NVSUB(TOTELE*NSUBNVLOC)
    REAL NWSUB(TOTELE*NSUBNVLOC)

    INTEGER ISPHERE,VELMEAN
    INTEGER NOFSTA
    REAL DT
    INTEGER MLOC,PGNDGLNO(MLOC*TOTELE)
    REAL FREES(PGNODS),FREESOLD(PGNODS)

    real, intent(in):: rnonl_coefficient

    ! Local variables...
    INTEGER SNLOC,SNGI,SNCLOC,NCLOC
    PARAMETER(SNLOC=3,SNCLOC=6,SNGI=7,NCLOC=10)
    REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI),SNLZ(SNLOC,SNGI)
    
    REAL SWEIGH(SNGI)
    REAL SNC(SNCLOC,SNGI),SNCLX(SNCLOC,SNGI),SNCLY(SNCLOC,SNGI),SNCLZ(SNCLOC,SNGI)
    REAL SNCX(SNCLOC,SNGI),SNCY(SNCLOC,SNGI),SNCZ(SNCLOC,SNGI) 
    REAL XSL(SNCLOC),YSL(SNCLOC),ZSL(SNCLOC)
    REAL SUD(SNGI),SVD(SNGI),SWD(SNGI)
    REAL SABSGI(SNGI),SSPEED(SNGI),SDEPTH(SNGI),GIFREES(SNGI),GIFREESOLD(SNGI)
    REAL SL1(SNGI), SL2(SNGI), SL3(SNGI), SL4(SNGI)
    
    REAL SDETWE(SNGI)
    
    REAL A11,A12,A13
    REAL A21,A22,A23
    REAL A31,A32,A33
    REAL AGI,BGI,CGI
    REAL DGI,EGI,FGI
    REAL GGI,HGI,KGI
    REAL SA11(SNGI),SA12(SNGI),SA13(SNGI)
    REAL SA21(SNGI),SA22(SNGI),SA23(SNGI)
    REAL SA31(SNGI),SA32(SNGI),SA33(SNGI)
    REAL SI11(SNGI),SI12(SNGI),SI13(SNGI)
    REAL SI21(SNGI),SI22(SNGI),SI23(SNGI)
    REAL SI31(SNGI),SI32(SNGI),SI33(SNGI)
    REAL SKSTAB(SNGI),DXSQRTSTA(SNGI),SQRTSTA(SNGI)
    
    REAL NORGRAVX,NORGRAVY,NORGRAVZ
    REAL NGDOTNORM(SNGI)

    REAL SNORMXN(SNGI),SNORMYN(SNGI),SNORMZN(SNGI)
    REAL SDIRX(SNGI),SDIRY(SNGI),SDIRZ(SNGI)
    REAL LOCPNORMX(NCLOC),LOCPNORMY(NCLOC),LOCPNORMZ(NCLOC)

    INTEGER SILOC2ILOC(SNLOC),OTHILOC(NLOC)
   
    INTEGER MLOCSILOC2ILOC(SNCLOC),MLOCOTHILOC(NCLOC)
   
    INTEGER ELE,ILOC,JLOC
    INTEGER SGI,GLOBI,GLOBJ
    INTEGER COUNT2

    REAL NDOTQ(SNGI),INCOME(SNGI)
    REAL SAREA
    REAL OME2,DX
    INTEGER ELE2
    INTEGER IGL
    INTEGER SILOC,SJLOC,SKLOC,IGLX,JGLX
    INTEGER KLOC,INOD,JNOD
    INTEGER POSMAT
    REAL RADI,RADJ,RADP,RAD,XDGI,YDGI,ZDGI,RGI
    REAL RNN,OUTVLMX,OUTVLMY,OUTVLMZ,RN
    
    REAL NORMX,NORMY,NORMZ
    
    REAL NNINVD,RESID,RNONL,KMATH,KKS
    
    REAL RTREST11,RTREST12,RTREST13
    REAL RTREST21,RTREST22,RTREST23
    REAL RTREST31,RTREST32,RTREST33
    REAL TEMSNCX,TEMSNCY,TEMSNCZ
    
    LOGICAL TOPSUF
    INTEGER ICENT,GLOBK,IDGNOD
    INTEGER IFACE,NFACE
    PARAMETER(NFACE=4)
    INTEGER LOCLIST(NFACE,3)

    INTEGER, ALLOCATABLE, DIMENSION(:)::FINELE
    INTEGER, ALLOCATABLE, DIMENSION(:)::COLELE4
    
    REAL, ALLOCATABLE, DIMENSION(:)::ZER

    ALLOCATE(FINELE(TOTELE+1))
    ALLOCATE(COLELE4(TOTELE*5))
    ALLOCATE(ZER(PGNODS))

    LONGDIFVECX(1:PGNODS)=0.0
    LONGDIFMATH(1:NPGCOLM)=0.0

    ZER(1:PGNODS)=0.0

    ewrite(1,*) 'in LONGALLTOPSUFDG'

    ! Calculate surface shape functions SN, SM etc...
    CALL TRIQUA(SL1, SL2, SL3, SL3, SWEIGH, .FALSE.,SNGI)
    CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
         SNLOC,SNGI, &
         SN,SNLX,SNLY,SNLZ)
    CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
         SNCLOC,SNGI, &
         SNC,SNCLX,SNCLY,SNCLZ)

    ewrite(1,*) 'going into GETFINELE4'
    CALL GETFINELE4(TOTELE,NLOC,SNLOC,NDGLNO, COLELE4,NONODS)

    ! The local cords are in anti-clockwise order
    IFACE=1
    LOCLIST(IFACE,1)=1
    LOCLIST(IFACE,2)=2
    LOCLIST(IFACE,3)=3
    IFACE=2
    LOCLIST(IFACE,1)=1
    LOCLIST(IFACE,2)=2
    LOCLIST(IFACE,3)=4
    IFACE=3
    LOCLIST(IFACE,1)=3
    LOCLIST(IFACE,2)=2
    LOCLIST(IFACE,3)=4
    IFACE=4
    LOCLIST(IFACE,1)=1
    LOCLIST(IFACE,2)=3
    LOCLIST(IFACE,3)=4

    element_loop: DO ELE=1,TOTELE

       ! Put surface integrals into matrix and rhs*****************************
       face_loop: DO IFACE=1,4
          COUNT2=(ELE-1)*5+IFACE
          ! NB colele4 is ordered in terms of faces.
          ELE2=COLELE4(COUNT2)

          IF(ELE2.EQ.0) THEN

             ! This is for a surface element...
             IF(COGRAX) THEN
                DO KLOC=1,NCLOC
                   LOCPNORMX(KLOC)=0.0
                   LOCPNORMY(KLOC)=0.0
                   LOCPNORMZ(KLOC)=1.0
                END DO
             ELSE
                ! Assume its a spherical Earth...
                DO ILOC=1,NLOC
                   INOD=XONDGL((ELE-1)*NLOC+ILOC)
                   DO JLOC=ILOC,NLOC
                      JNOD=XONDGL((ELE-1)*NLOC+JLOC)
                      KLOC=ILINK2(ILOC,JLOC)
                      LOCPNORMX(KLOC)=0.5*(XORIG(INOD)+XORIG(JNOD))
                      LOCPNORMY(KLOC)=0.5*(YORIG(INOD)+YORIG(JNOD))
                      LOCPNORMZ(KLOC)=0.5*(ZORIG(INOD)+ZORIG(JNOD))
                   END DO
                END DO
                DO KLOC=1,NCLOC
                   RN=SQRT(LOCPNORMX(KLOC)**2+LOCPNORMY(KLOC)**2+LOCPNORMZ(KLOC)**2)
                   LOCPNORMX(KLOC)=LOCPNORMX(KLOC)/RN
                   LOCPNORMY(KLOC)=LOCPNORMY(KLOC)/RN
                   LOCPNORMZ(KLOC)=LOCPNORMZ(KLOC)/RN
                END DO
             ENDIF

             ! Surface element of domain...
             DO SILOC=1,SNLOC
                SILOC2ILOC(SILOC)=LOCLIST(IFACE,SILOC)
             END DO
             ! Form OTHILOC from SILOC2ILOC
             IF(ELE2.NE.0) THEN
                CALL ELE2ELEFACE(OTHILOC,SILOC2ILOC,SNLOC,NLOC,TOTELE, &
                     XONDGL,NDGLNO,ELE,ELE2)
             ELSE
                DO ILOC=1,NLOC
                   OTHILOC(ILOC)=ILOC
                END DO
             ENDIF

             ! Form MLOCOTHILOC and MLOCSILOC2ILOC
             CALL MLOCELE2ELEFACE(MLOCOTHILOC,MLOCSILOC2ILOC,SNLOC,NLOC, &
                  SNCLOC,NCLOC,OTHILOC,SILOC2ILOC)

             ! PERFORM SURFACE INTEGRATION         !
             ! Form approximate surface normal (NORMX,NORMY,NORMZ)
             CALL DGSIMPLNORM(ELE,SILOC2ILOC,TOTELE,NLOC,SNLOC,XONDGL, &
                  XORIG,YORIG,ZORIG,XNONOD,NORMX,NORMY,NORMZ)

             CALL DGSDETNXLOC(ELE,SILOC2ILOC,XONDGL, &
                  TOTELE,NLOC,XNONOD,SNLOC,SNGI, &
                  XORIG,YORIG,ZORIG, &
                  SN,SNLX,SNLY, SWEIGH, SDETWE,SAREA, .TRUE.,.FALSE., &
                  SNORMXN,SNORMYN,SNORMZN, &
                  NORMX,NORMY,NORMZ, &
                  ISPHERE,SNCLOC,SNC,SNCLX,SNCLY)

             TOPSUF=.FALSE.
             DO SGI=1,SNGI
                SDIRX(SGI)=0.0
                SDIRY(SGI)=0.0
                SDIRZ(SGI)=0.0
                DO SILOC=1,SNCLOC
                   ILOC=MLOCSILOC2ILOC(SILOC)
                   SDIRX(SGI)=SDIRX(SGI)-SNC(SILOC,SGI)*LOCPNORMX(ILOC)
                   SDIRY(SGI)=SDIRY(SGI)-SNC(SILOC,SGI)*LOCPNORMY(ILOC)
                   SDIRZ(SGI)=SDIRZ(SGI)-SNC(SILOC,SGI)*LOCPNORMZ(ILOC)
                END DO
                RN=SQRT(SDIRX(SGI)**2+SDIRY(SGI)**2+SDIRZ(SGI)**2)
                SDIRX(SGI)=SDIRX(SGI)/RN
                SDIRY(SGI)=SDIRY(SGI)/RN
                SDIRZ(SGI)=SDIRZ(SGI)/RN
                NDOTQ(SGI)=SNORMXN(SGI)*SDIRX(SGI)+SNORMYN(SGI)*SDIRY(SGI) &
                     +SNORMZN(SGI)*SDIRZ(SGI)
                INCOME(SGI)=0.0
                IF(NDOTQ(SGI).LT.0.0) INCOME(SGI)=1.0
                IF(ELE2.EQ.0) THEN
                   IF(NDOTQ(SGI).LT.-0.8) TOPSUF=.TRUE.
                ENDIF
             END DO

             IF(TOPSUF) THEN

                DO SILOC=1,SNLOC
                   ILOC=SILOC2ILOC(SILOC)
                   INOD=XONDGL((ELE-1)*NLOC+ILOC)
                   DO SJLOC=SILOC,SNLOC,1
                      JLOC=SILOC2ILOC(SJLOC)
                      JNOD=XONDGL((ELE-1)*NLOC+JLOC)
                      KLOC=ILINK2(ILOC,JLOC)
                      SKLOC=SILINK2(SILOC,SJLOC)
                      XSL(SKLOC)=0.5*(XORIG(INOD)+XORIG(JNOD))
                      YSL(SKLOC)=0.5*(YORIG(INOD)+YORIG(JNOD))
                      ZSL(SKLOC)=0.5*(ZORIG(INOD)+ZORIG(JNOD))
                      ! Assume the mid side nodes are on the sphere
                      IF(ISPHERE.GE.2) THEN
                         IF(ILOC.NE.JLOC) THEN
                            RADI=SQRT(XORIG(INOD)**2+YORIG(INOD)**2+ZORIG(INOD)**2)
                            RADJ=SQRT(XORIG(JNOD)**2+YORIG(JNOD)**2+ZORIG(JNOD)**2)
                            RADP=0.5*(RADI+RADJ)
                            RAD=SQRT(XSL(SKLOC)**2+YSL(SKLOC)**2+ZSL(SKLOC)**2)
                            XSL(SKLOC)=XSL(SKLOC)*RADP/RAD
                            YSL(SKLOC)=YSL(SKLOC)*RADP/RAD
                            ZSL(SKLOC)=ZSL(SKLOC)*RADP/RAD
                         ENDIF
                      ENDIF
                   END DO
                END DO
                ! Now amend (XSL,YSL,ZSL) in the direction of the free surface
                ! change...
                DO SKLOC=1,SNCLOC
                   KLOC=MLOCSILOC2ILOC(SKLOC)
                   GLOBK =PGNDGLNO((ELE-1)*MLOC+KLOC)
                   XSL(SKLOC)=XSL(SKLOC) + LOCPNORMX(KLOC)*FREES(GLOBK)
                   YSL(SKLOC)=YSL(SKLOC) + LOCPNORMY(KLOC)*FREES(GLOBK)
                   ZSL(SKLOC)=ZSL(SKLOC) + LOCPNORMZ(KLOC)*FREES(GLOBK)
                END DO

                ! Recalculate the normal...
                CALL DGSDETNXLOC2(SNCLOC,SNGI, &
                     XSL,YSL,ZSL, &
                     SNC,SNCLX,SNCLY, SWEIGH, SDETWE,SAREA,.TRUE.,.FALSE., &
                     SNORMXN,SNORMYN,SNORMZN, &
                     NORMX,NORMY,NORMZ)

                sgi_loop: DO SGI=1,SNGI

                   AGI=0.
                   BGI=0.
                   CGI=0.

                   DGI=0.
                   EGI=0.
                   FGI=0.

                   XDGI=0.0
                   YDGI=0.0
                   ZDGI=0.0
                   DO SILOC=1,SNCLOC
                      ! NB R0 does not appear here although the z-coord might be Z+R0.
                      AGI=AGI+SNCLX(SILOC,SGI)*XSL(SILOC)
                      BGI=BGI+SNCLX(SILOC,SGI)*YSL(SILOC)
                      CGI=CGI+SNCLX(SILOC,SGI)*ZSL(SILOC)

                      DGI=DGI+SNCLY(SILOC,SGI)*XSL(SILOC)
                      EGI=EGI+SNCLY(SILOC,SGI)*YSL(SILOC)
                      FGI=FGI+SNCLY(SILOC,SGI)*ZSL(SILOC)

                      XDGI=XDGI+SNC(SILOC,SGI)*XSL(SILOC)
                      YDGI=YDGI+SNC(SILOC,SGI)*YSL(SILOC)
                      ZDGI=ZDGI+SNC(SILOC,SGI)*ZSL(SILOC)
                   end do
                   IF(COGRAX) THEN
                      GGI=0.0
                      HGI=0.0
                      KGI=1.0
                      NORGRAVX=0.0
                      NORGRAVY=0.0
                      NORGRAVZ=1.0
                   ELSE
                      RGI=SQRT(XDGI**2+YDGI**2+ZDGI**2)
                      GGI=RGI/MAX(1.0E-10,SQRT(RGI**2-(YDGI**2+ZDGI**2)))
                      HGI=RGI/MAX(1.0E-10,SQRT(RGI**2-(XDGI**2+ZDGI**2)))
                      KGI=RGI/MAX(1.0E-10,SQRT(RGI**2-(XDGI**2+YDGI**2)))
                      NORGRAVX=XDGI/RGI
                      NORGRAVY=YDGI/RGI
                      NORGRAVZ=ZDGI/RGI
                   ENDIF
                   NGDOTNORM(SGI)=NORGRAVX*SNORMXN(SGI)+NORGRAVY*SNORMYN(SGI) &
                        +NORGRAVZ*SNORMZN(SGI)

                   CALL INV3X3( &
                        AGI,BGI,CGI,DGI,EGI,FGI,GGI,HGI,KGI, &
                        A11,A12,A13,A21,A22,A23,A31,A32,A33)

                   DO SILOC=1,SNCLOC
                      SNCX(SILOC,SGI)= A11*SNCLX(SILOC,SGI)+A12*SNCLY(SILOC,SGI)
                      SNCY(SILOC,SGI)= A21*SNCLX(SILOC,SGI)+A22*SNCLY(SILOC,SGI)
                      SNCZ(SILOC,SGI)= A31*SNCLX(SILOC,SGI)+A32*SNCLY(SILOC,SGI)
                   end do
                   ! rotate and rotate back qithout vertical (E.G. r on
                   !  spherical Earth) component:

                   RTREST11=1.0 - NORGRAVX*NORGRAVX
                   RTREST12=    - NORGRAVX*NORGRAVY
                   RTREST13=    - NORGRAVX*NORGRAVZ

                   RTREST21=    - NORGRAVY*NORGRAVX
                   RTREST22=1.0 - NORGRAVY*NORGRAVY
                   RTREST23=    - NORGRAVY*NORGRAVZ

                   RTREST31=    - NORGRAVZ*NORGRAVX
                   RTREST32=    - NORGRAVZ*NORGRAVY
                   RTREST33=1.0 - NORGRAVZ*NORGRAVZ
                   DO SILOC=1,SNCLOC
                      TEMSNCX=SNCX(SILOC,SGI)
                      TEMSNCY=SNCY(SILOC,SGI)
                      TEMSNCZ=SNCZ(SILOC,SGI)
                      SNCX(SILOC,SGI)= RTREST11*TEMSNCX+RTREST12*TEMSNCY+RTREST13*TEMSNCZ
                      SNCY(SILOC,SGI)= RTREST21*TEMSNCX+RTREST22*TEMSNCY+RTREST23*TEMSNCZ
                      SNCZ(SILOC,SGI)= RTREST31*TEMSNCX+RTREST32*TEMSNCY+RTREST33*TEMSNCZ
                   END DO

                   SUD(SGI)=0.0
                   SVD(SGI)=0.0
                   SWD(SGI)=0.0
                   SABSGI(SGI)=0.0
                   DO SILOC=1,SNLOC
                      ILOC=SILOC2ILOC(SILOC)
                      IGL =NDGLNO((ELE-1)*NLOC+ILOC)
                      SUD(SGI)=SUD(SGI) + SN(SILOC,SGI)*U(IGL)
                      SVD(SGI)=SVD(SGI) + SN(SILOC,SGI)*V(IGL)
                      SWD(SGI)=SWD(SGI) + SN(SILOC,SGI)*W(IGL)
                      IF((INUSUB.NE.0).AND.(NSUBNVLOC.NE.0)) THEN
                         IDGNOD=(ELE-1)*NLOC+ILOC
                         SUD(SGI)=SUD(SGI) + SN(SILOC,SGI)*NUSUB(IDGNOD)
                         SVD(SGI)=SVD(SGI) + SN(SILOC,SGI)*NVSUB(IDGNOD)
                         SWD(SGI)=SWD(SGI) + SN(SILOC,SGI)*NWSUB(IDGNOD)
                      ENDIF
                      SABSGI(SGI)=SABSGI(SGI)+SN(SILOC,SGI)*ABSORB(IGL)
                   END DO
                   SSPEED(SGI)=SQRT(SUD(SGI)**2+SVD(SGI)**2+SWD(SGI)**2)
                   SDEPTH(SGI)=0.0
                   DO SILOC=1,SNLOC
                      ILOC=SILOC2ILOC(SILOC)
                      IGLX=XONDGL((ELE-1)*NLOC+ILOC)
                      SDEPTH(SGI)=SDEPTH(SGI)+SN(SILOC,SGI)*(TOPDIS(IGLX)+BOTDIS(IGLX))
                   END DO
                   GIFREES(SGI)=0.0
                   GIFREESOLD(SGI)=0.0
                   DO SJLOC=1,SNCLOC
                      JLOC=MLOCSILOC2ILOC(SJLOC)
                      GLOBJ =PGNDGLNO((ELE-1)*MLOC+JLOC)
                      GIFREES(SGI)=GIFREES(SGI)+       SNC(SJLOC,SGI)*FREES(GLOBJ)
                      GIFREESOLD(SGI)=GIFREESOLD(SGI)+ SNC(SJLOC,SGI)*FREESOLD(GLOBJ)
                   END DO

                   OME2=2.*FUNOME(XDGI,YDGI,ZDGI)
                   IF((GEOBAL.EQ.-12).or.(GEOBAL.EQ.-22)) OME2=0.
                   SA11(SGI)=1.0         + DT*ATHETA*SABSGI(SGI)
                   SA12(SGI)=-DT*CTHETA*OME2
                   SA13(SGI)=0.0
                   SA21(SGI)=DT*CTHETA*OME2
                   SA22(SGI)=1.0         + DT*ATHETA*SABSGI(SGI)
                   SA23(SGI)=0.0
                   SA31(SGI)=0.0
                   SA32(SGI)=0.0
                   SA33(SGI)=1.0         + DT*ATHETA*SABSGI(SGI)
                   CALL INV3X3( &
                        SA11(SGI),SA12(SGI),SA13(SGI),SA21(SGI),SA22(SGI),SA23(SGI),SA31(SGI),SA32(SGI),SA33(SGI), &
                        SI11(SGI),SI12(SGI),SI13(SGI),SI21(SGI),SI22(SGI),SI23(SGI),SI31(SGI),SI32(SGI),SI33(SGI))
                END DO sgi_loop

                ! Calculate DX - the characteristic horizontal element length scale
                ! IF NOFSTA=1 switch off stabilisation =0 include it.
                DX=0.0
                DO SILOC=1,SNLOC
                   ILOC=SILOC2ILOC(SILOC)
                   IGLX=XONDGL((ELE-1)*NLOC+ILOC)
                   DO SJLOC=SILOC+1,SNLOC,1
                      JLOC=SILOC2ILOC(SJLOC)
                      JGLX=XONDGL((ELE-1)*NLOC+JLOC)
                      DX=MAX(DX, ABS(XORIG(IGLX)-XORIG(JGLX)), &
                           ABS(YORIG(IGLX)-YORIG(JGLX)), ABS(ZORIG(IGLX)-ZORIG(JGLX)) )
                   END DO
                END DO
                IF(NOFSTA.EQ.0) THEN
                   DO SGI=1,SNGI
                      SKSTAB(SGI)=(  (dx*sqrt(gravty*Sdepth(Sgi))/dt)   )/(1.+DT*ATHETA*SABSGI(SGI))
                      IF(NONLIN) THEN
                         RESID=-(GIFREES(SGI)-GIFREESOLD(SGI))/DT &
                              +SNORMXN(SGI)*SUD(SGI)+SNORMYN(SGI)*SVD(SGI)+SNORMZN(SGI)*SWD(SGI)
                         ! 0.1 is too dissipative for tsunami and 2nd option
                         !  below.
                         !            RNONL=0.01*(DX/Sdepth(Sgi))*ABS(RESID)/MAX(SSPEED(SGI),1.E-10)
                         rnonl=rnonl_coefficient*abs(resid)/ &
                              max(abs(snormxn(sgi)*sud(sgi)+snormyn(sgi)&
                              &*svd(sgi)+snormzn(sgi)*swd(sgi)),1.e-10)
                         SKSTAB(SGI)=SKSTAB(SGI)*MIN(1.0,RNONL)
                      ENDIF
                      SQRTSTA(SGI)=SQRT(SKSTAB(SGI))
                      DXSQRTSTA(SGI)=DX*SQRTSTA(SGI)
                   END DO
                ELSE
                   SKSTAB(1:SNGI)=0.0
                   SQRTSTA(1:SNGI)=0.0
                   DXSQRTSTA(1:SNGI)=0.0
                ENDIF

                ! Perform surface integration...
                siloc_loop: DO SILOC=1,SNCLOC
                   ILOC=MLOCSILOC2ILOC(SILOC)
                   GLOBI =PGNDGLNO((ELE-1)*MLOC+ILOC)
                   sjloc_loop: DO SJLOC=1,SNCLOC
                      JLOC=MLOCSILOC2ILOC(SJLOC)
                      GLOBJ =PGNDGLNO((ELE-1)*MLOC+JLOC)

                      NNINVD=0.0

                      KMATH=0.0
                      KKS=0.0

                      DO SGI=1,SNGI
                         RNN=SNC(SILOC,SGI)*SNC(SJLOC,SGI)*SDETWE(SGI)
                         NNINVD=NNINVD+1.0*NGDOTNORM(SGI)*RNN

                         KMATH=KMATH+FSTHETA*ngdotnorm(sgi)*GRAVTY*SDETWE(SGI)*SDEPTH(SGI)*( &
                              SNCX(SILOC,SGI)*(SI11(SGI)*SNCX(SJLOC,SGI)+SI12(SGI)*SNCY(SJLOC,SGI)+SI13(SGI)*SNCZ(SJLOC,SGI)) &
                              + SNCY(SILOC,SGI)*(SI21(SGI)*SNCX(SJLOC,SGI)+SI22(SGI)*SNCY(SJLOC,SGI)+SI23(SGI)*SNCZ(SJLOC,SGI)) &
                              + SNCZ(SILOC,SGI)*(SI31(SGI)*SNCX(SJLOC,SGI)&
                              &+SI32(SGI)*SNCY(SJLOC,SGI)+SI33(SGI)*SNCZ(SJLOC,SGI)))

                         IF(NOFSTA.EQ.0) THEN
                            KKS=KKS+SKSTAB(SGI)*SDETWE(SGI)*(SNCX(SILOC,SGI)*SNCX(SJLOC,SGI) &
                                 +SNCY(SILOC,SGI)*SNCY(SJLOC,SGI)+SNCZ(SILOC,SGI)*SNCZ(SJLOC,SGI))
                         ENDIF
                      END DO
                      ! Find POSMAT. ***************************************
                      CALL POSINMAT(POSMAT,GLOBI,GLOBJ, &
                           PGNODS,PGFINDRM,PGCOLM,NPGCOLM)
                      ! Found POSMAT ***************************************

                      ICENT=PGCENTRM(GLOBI)
                      IF(VELMEAN.EQ.1) THEN
                         LONGDIFMATH(POSMAT)=LONGDIFMATH(POSMAT)+KMATH +NNINVD/(DT**2)
                      ELSE
                         LONGDIFMATH(POSMAT)=LONGDIFMATH(POSMAT)+KMATH +NNINVD/(DT**2) + KKS
                      ENDIF

                   END DO sjloc_loop
                   ! Now for the vector....
                   DO SJLOC=1,SNCLOC
                      JLOC=MLOCSILOC2ILOC(SJLOC)
                      GLOBJ =PGNDGLNO((ELE-1)*MLOC+JLOC)

                      NNINVD=0.0
                      KKS=0.0
                      DO SGI=1,SNGI
                         NNINVD=NNINVD+NGDOTNORM(SGI)*SNC(SILOC,SGI)*SNC(SJLOC,SGI)*SDETWE(SGI)
                         IF(NOFSTA.EQ.0) THEN
                            KKS=KKS+SKSTAB(SGI)*SDETWE(SGI)*(SNCX(SILOC,SGI)*SNCX(SJLOC,SGI) &
                                 +SNCY(SILOC,SGI)*SNCY(SJLOC,SGI)+SNCZ(SILOC,SGI)*SNCZ(SJLOC,SGI))
                         ENDIF
                      END DO
                      ! The matrix eqn associated with velocity star
                      IF(VELMEAN.EQ.1) THEN
                         LONGDIFVECX(GLOBI)=LONGDIFVECX(GLOBI)-(NNINVD/(DT**2))*(FREES(GLOBJ)-FREESOLD(GLOBJ))
                      ELSE
                         LONGDIFVECX(GLOBI)=LONGDIFVECX(GLOBI)-(NNINVD/(DT**2))*(FREES(GLOBJ)-FREESOLD(GLOBJ))
                         LONGDIFVECX(GLOBI)=LONGDIFVECX(GLOBI)- KKS*FREES(GLOBJ)
                      ENDIF
                   END DO

                   ! surface flux term:
                   OUTVLMX=0.0
                   OUTVLMY=0.0
                   OUTVLMZ=0.0
                   DO SGI=1,SNGI
                      RNN=SNC(SILOC,SGI)*SDETWE(SGI)
                      OUTVLMX=OUTVLMX+RNN*SNORMXN(SGI)*SUD(SGI)
                      OUTVLMY=OUTVLMY+RNN*SNORMYN(SGI)*SVD(SGI)
                      OUTVLMZ=OUTVLMZ+RNN*SNORMZN(SGI)*SWD(SGI)
                   END DO
                   LONGDIFVECX(GLOBI)=LONGDIFVECX(GLOBI) &
                        +(OUTVLMX+OUTVLMY+OUTVLMZ)/DT

                END DO siloc_loop
             ENDIF
          ENDIF

       END DO face_loop

    END DO element_loop

  END SUBROUTINE LONGALLTOPSUFDG

  subroutine alltopsufdg(shortdifmath,shortdifvecx, &
       wetdrymat,sufml, &
       xorig,yorig,zorig, &
       topdis,botdis,atheta,absorb, &
       gravty,ctheta, &
       nloc, &
       totele, &
       geobal,&
       ndglno,xondgl, &
       nonods,xnonod, &
       scfacth0,isphere, &
       findrm,colm,ncolm, &
       fincmc,midcmc,colcmc,ncmc,fredop, &
       cograx, &
       u,v,w, &
       INUSUB,NSUBNVLOC,NUSUB,NVSUB,NWSUB, &
       nofsta,nonlin,velmean,dt, & 
       pgnods,mloc,pgndglno,frees,freesold, &
       fstheta, rnonl_coefficient, wetdry,colele4,wdgravity)
    ! This subroutine caluculates free surface mass and 
    ! surface momentum flux contributions by NOT integrating over the
    ! top surface of the ocean.
    ! M is the basis function assocaited with normal pressure.
    ! the 4th order switch switches on the lumping.
    use FLDebug
    use assnav_module
    use assnav_module

    IMPLICIT NONE
    ! Number of source nodes, advection nodes and node positions.
    ! Number of elements
    ! IF SMALLSTEN then use a small stencil other wise use a more accurate
    ! enlarged stencil.
    LOGICAL SLUMMAT,SLUMRHS,SMALLSTEN
    PARAMETER(SLUMMAT=.FALSE.,SLUMRHS=.FALSE.,SMALLSTEN=.TRUE.)
    ! IF SLUMMAT then lump the free surface height mass matrix.
    ! IF SLUMRHS then lump the rhs contribution to matrix eqn.
    INTEGER TOTELE,NONODS,XNONOD
    INTEGER NLOC

    REAL XORIG(XNONOD),YORIG(XNONOD),ZORIG(XNONOD)
    ! Global node number maps.
    INTEGER NDGLNO(TOTELE*NLOC)
    INTEGER XONDGL(TOTELE*NLOC)
    REAL GRAVTY,CTHETA,ATHETA
    REAL TOPDIS(XNONOD),BOTDIS(XNONOD)
    REAL ABSORB(NONODS)
    real scfacth0
    INTEGER GEOBAL
    ! This is the THETA time stepping parameter for the free surface.
    REAL FSTHETA
    real, intent(in):: rnonl_coefficient
    INTEGER WETDRY
    INTEGER COLELE4(TOTELE*5)
    REAL SHORTDIFMATH(NCMC)
    REAL SHORTDIFVECX(NONODS)
    INTEGER NCOLM
    INTEGER FINDRM(NONODS+1),COLM(NCOLM)
    INTEGER FREDOP,NCMC,FINCMC(FREDOP+1),COLCMC(NCMC),MIDCMC(FREDOP)
    LOGICAL COGRAX
    REAL U(NONODS),V(NONODS),W(NONODS)
    INTEGER INUSUB,NSUBNVLOC
    REAL NUSUB(TOTELE*NSUBNVLOC),NVSUB(TOTELE*NSUBNVLOC)
    REAL NWSUB(TOTELE*NSUBNVLOC)

    INTEGER ISPHERE,VELMEAN
    INTEGER NOFSTA
    LOGICAL NONLIN
    REAL DT
    INTEGER PGNODS,MLOC,PGNDGLNO(MLOC*TOTELE)
    REAL FREES(PGNODS),FREESOLD(PGNODS)

    REAL SUFML(NONODS),WETDRYMAT(NCMC*WETDRY)

    ! Local variables...
    INTEGER SNLOC,SNGI,SNCLOC,NCLOC
    PARAMETER(SNLOC=3,SNCLOC=6,SNGI=7,NCLOC=10)
    REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI),SNLZ(SNLOC,SNGI)
    REAL SNX(SNLOC,SNGI),SNY(SNLOC,SNGI),SNZ(SNLOC,SNGI)
   
    REAL SWEIGH(SNGI)
    REAL SNC(SNCLOC,SNGI),SNCLX(SNCLOC,SNGI),SNCLY(SNCLOC,SNGI),SNCLZ(SNCLOC,SNGI)
    REAL SNCX(SNCLOC,SNGI),SNCY(SNCLOC,SNGI),SNCZ(SNCLOC,SNGI)
    REAL XSL(SNCLOC),YSL(SNCLOC),ZSL(SNCLOC)
    REAL SUD(SNGI),SVD(SNGI),SWD(SNGI)
    REAL SABSGI(SNGI),SSPEED(SNGI),SDEPTH(SNGI),GIFREES(SNGI),GIFREESOLD(SNGI)
    REAL SL1(SNGI), SL2(SNGI), SL3(SNGI), SL4(SNGI)
 
    REAL SDETWE(SNGI)
   
    REAL A11,A12,A13
    REAL A21,A22,A23
    REAL A31,A32,A33
    REAL AGI,BGI,CGI
    REAL DGI,EGI,FGI
    REAL GGI,HGI,KGI
    REAL SA11(SNGI),SA12(SNGI),SA13(SNGI)
    REAL SA21(SNGI),SA22(SNGI),SA23(SNGI)
    REAL SA31(SNGI),SA32(SNGI),SA33(SNGI)
    REAL SI11(SNGI),SI12(SNGI),SI13(SNGI)
    REAL SI21(SNGI),SI22(SNGI),SI23(SNGI)
    REAL SI31(SNGI),SI32(SNGI),SI33(SNGI)
    REAL SKSTAB(SNGI),DXSQRTSTA(SNGI),SQRTSTA(SNGI)
    
    REAL NORGRAVX,NORGRAVY,NORGRAVZ
    REAL NGDOTNORM(SNGI)
    REAL SQRTDISPX(SNGI),SQRTDISPY(SNGI),SQRTDISPZ(SNGI)
    REAL SNORMXN(SNGI),SNORMYN(SNGI),SNORMZN(SNGI)
    REAL SDIRX(SNGI),SDIRY(SNGI),SDIRZ(SNGI)
    REAL LOCPNORMX(NCLOC),LOCPNORMY(NCLOC),LOCPNORMZ(NCLOC)

    INTEGER SILOC2ILOC(SNLOC),OTHILOC(NLOC)

    INTEGER MLOCSILOC2ILOC(SNCLOC),MLOCOTHILOC(NCLOC)
  
    INTEGER ELE,ILOC,JLOC
    INTEGER SGI,GLOBI,GLOBJ,GLOBK
    INTEGER COUNT2
    
    REAL NDOTQ(SNGI),INCOME(SNGI)
    REAL UTOT(NLOC),VTOT(NLOC),WTOT(NLOC)
    REAL SAREA
    REAL OME2,DX
    INTEGER ELE2
    INTEGER IGL,IDGNOD
    INTEGER SILOC,SJLOC,SKLOC,IGLX,JGLX
    INTEGER KLOC,INOD,JNOD
    INTEGER POSMAT,POSCT
    REAL RADI,RADJ,RADP,RAD,XDGI,YDGI,ZDGI,RGI
    REAL RNN,OUTVLMX,OUTVLMY,OUTVLMZ,RN
   
    REAL NORMX,NORMY,NORMZ
    
    REAL NNINVD,RESID,RNONL,KMATH,KKS
    REAL SNSNX,SNSNY,SNSNZ
    REAL SNSNXX,SNSNYY,SNSNZZ
    REAL NNMLX,NNMLY,NNMLZ,STANNML
    
    REAL RTREST11,RTREST12,RTREST13
    REAL RTREST21,RTREST22,RTREST23
    REAL RTREST31,RTREST32,RTREST33
    REAL TEMSNX,TEMSNY,TEMSNZ
    REAL TEMSNCX,TEMSNCY,TEMSNCZ
    REAL CONSUM,CONOUT
    LOGICAL TOPSUF,LUMMAT,LUMMAT2,LUMMAT3
    INTEGER ICENT
    INTEGER IFACE,NFACE
    PARAMETER(NFACE=4)
    INTEGER LOCLIST(NFACE,3)

    Real, intent(in) :: wdgravity(totele)

    REAL, ALLOCATABLE, DIMENSION(:)::vec
    INTEGER, ALLOCATABLE, DIMENSION(:)::FINELE
    REAL, ALLOCATABLE, DIMENSION(:)::CXT
    REAL, ALLOCATABLE, DIMENSION(:)::CYT
    REAL, ALLOCATABLE, DIMENSION(:)::CZT
    REAL, ALLOCATABLE, DIMENSION(:)::CXXT
    REAL, ALLOCATABLE, DIMENSION(:)::CYYT
    REAL, ALLOCATABLE, DIMENSION(:)::CZZT
    REAL, ALLOCATABLE, DIMENSION(:)::MLX
    REAL, ALLOCATABLE, DIMENSION(:)::MLY
    REAL, ALLOCATABLE, DIMENSION(:)::MLZ
    REAL, ALLOCATABLE, DIMENSION(:)::STAML
    REAL, ALLOCATABLE, DIMENSION(:)::ZER
    REAL, ALLOCATABLE, DIMENSION(:)::MSUFML

    if (.not. SMALLSTEN) then
      FLAbort("SMALLSTEN==.false. no longer works")
    end if

    ALLOCATE(vec(nonods))
    ALLOCATE(FINELE(TOTELE+1))
    ALLOCATE(CXT(NCOLM))
    ALLOCATE(CYT(NCOLM))
    ALLOCATE(CZT(NCOLM))
    ALLOCATE(CXXT(NCOLM))
    ALLOCATE(CYYT(NCOLM))
    ALLOCATE(CZZT(NCOLM))
    ALLOCATE(MLX(FREDOP))
    ALLOCATE(MLY(FREDOP))
    ALLOCATE(MLZ(FREDOP))
    ALLOCATE(STAML(FREDOP))
    ALLOCATE(ZER(FREDOP))
    ALLOCATE(MSUFML(PGNODS))
    MSUFML=0.0

    vec(1:nonods) = 0.0
    SHORTDIFVECX(1:NONODS) = 0.0
    SHORTDIFMATH(1:NCMC) = 0.0
    IF(WETDRY.EQ.1) WETDRYMAT(1:NCMC)=0.0

    CXT(1:NCOLM) = 0.0
    CYT(1:NCOLM) = 0.0
    CZT(1:NCOLM) = 0.0
    CXXT(1:NCOLM) = 0.0
    CYYT(1:NCOLM) = 0.0
    CZZT(1:NCOLM) = 0.0
    MLX(1:FREDOP) = 0.0
    MLY(1:FREDOP) = 0.0
    MLZ(1:FREDOP) = 0.0
    STAML(1:FREDOP) = 0.0
    ZER(1:FREDOP) = 0.0
    SUFML(1:NONODS) = 0.0
    CONOUT=0.0

    ! Calculate surface shape functions SN, SM etc...
    CALL TRIQUA(SL1, SL2, SL3, SL3, SWEIGH, .FALSE.,SNGI)
    CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
         SNLOC,SNGI, &
         SN,SNLX,SNLY,SNLZ)
    CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
         SNCLOC,SNGI, &
         SNC,SNCLX,SNCLY,SNCLZ)

    ! CALCULATE ELEMENT TO ELEMENT LIST**************
    ! and put in order of faces

    ! The local cords are in anti-clockwise order
    IFACE=1
    LOCLIST(IFACE,1)=1
    LOCLIST(IFACE,2)=2
    LOCLIST(IFACE,3)=3
    IFACE=2
    LOCLIST(IFACE,1)=1
    LOCLIST(IFACE,2)=2
    LOCLIST(IFACE,3)=4
    IFACE=3
    LOCLIST(IFACE,1)=3
    LOCLIST(IFACE,2)=2
    LOCLIST(IFACE,3)=4
    IFACE=4
    LOCLIST(IFACE,1)=1
    LOCLIST(IFACE,2)=3
    LOCLIST(IFACE,3)=4

    element_loop: DO ELE=1,TOTELE

       ! Put surface integrals into matrix and rhs*****************************
       face_loop: DO IFACE=1,4
          COUNT2=(ELE-1)*5+IFACE
          ! NB colele4 is ordered in terms of faces.
          ELE2=COLELE4(COUNT2)

          IF(ELE2.EQ.0) THEN

             ! This is for a surface element...
             IF(COGRAX) THEN
                DO KLOC=1,NCLOC
                   LOCPNORMX(KLOC)=0.0
                   LOCPNORMY(KLOC)=0.0
                   LOCPNORMZ(KLOC)=1.0
                END DO
             ELSE
                ! Assume its a spherical Earth...
                DO ILOC=1,NLOC
                   INOD=XONDGL((ELE-1)*NLOC+ILOC)
                   DO JLOC=ILOC,NLOC
                      JNOD=XONDGL((ELE-1)*NLOC+JLOC)
                      KLOC=ILINK2(ILOC,JLOC)
                      LOCPNORMX(KLOC)=0.5*(XORIG(INOD)+XORIG(JNOD))
                      LOCPNORMY(KLOC)=0.5*(YORIG(INOD)+YORIG(JNOD))
                      LOCPNORMZ(KLOC)=0.5*(ZORIG(INOD)+ZORIG(JNOD))
                   END DO
                END DO
                DO KLOC=1,NCLOC
                   RN=SQRT(LOCPNORMX(KLOC)**2+LOCPNORMY(KLOC)**2+LOCPNORMZ(KLOC)**2)
                   LOCPNORMX(KLOC)=LOCPNORMX(KLOC)/RN
                   LOCPNORMY(KLOC)=LOCPNORMY(KLOC)/RN
                   LOCPNORMZ(KLOC)=LOCPNORMZ(KLOC)/RN
                END DO
             ENDIF

             ! Surface element of domain...
             DO SILOC=1,SNLOC
                SILOC2ILOC(SILOC)=LOCLIST(IFACE,SILOC)
             END DO
             ! Form OTHILOC from SILOC2ILOC
             IF(ELE2.NE.0) THEN
                CALL ELE2ELEFACE(OTHILOC,SILOC2ILOC,SNLOC,NLOC,TOTELE, &
                     XONDGL,NDGLNO,ELE,ELE2)
             ELSE
                DO ILOC=1,NLOC
                   OTHILOC(ILOC)=ILOC
                END DO
             ENDIF

             ! Form MLOCOTHILOC and MLOCSILOC2ILOC
             CALL MLOCELE2ELEFACE(MLOCOTHILOC,MLOCSILOC2ILOC,SNLOC,NLOC, &
                  SNCLOC,NCLOC,OTHILOC,SILOC2ILOC)

             ! PERFORM SURFACE INTEGRATION         !
             ! Form approximate surface normal (NORMX,NORMY,NORMZ)
             CALL DGSIMPLNORM(ELE,SILOC2ILOC,TOTELE,NLOC,SNLOC,XONDGL, &
                  XORIG,YORIG,ZORIG,XNONOD,NORMX,NORMY,NORMZ)

             CALL DGSDETNXLOC(ELE,SILOC2ILOC,XONDGL, &
                  TOTELE,NLOC,XNONOD,SNLOC,SNGI, &
                  XORIG,YORIG,ZORIG, &
                  SN,SNLX,SNLY, SWEIGH, SDETWE,SAREA, .TRUE.,.FALSE., &
                  SNORMXN,SNORMYN,SNORMZN, &
                  NORMX,NORMY,NORMZ, &
                  ISPHERE,SNCLOC,SNC,SNCLX,SNCLY)

             TOPSUF=.FALSE.
             DO SGI=1,SNGI
                SDIRX(SGI)=0.0
                SDIRY(SGI)=0.0
                SDIRZ(SGI)=0.0
                DO SILOC=1,SNCLOC
                   ILOC=MLOCSILOC2ILOC(SILOC)
                   SDIRX(SGI)=SDIRX(SGI)-SNC(SILOC,SGI)*LOCPNORMX(ILOC)
                   SDIRY(SGI)=SDIRY(SGI)-SNC(SILOC,SGI)*LOCPNORMY(ILOC)
                   SDIRZ(SGI)=SDIRZ(SGI)-SNC(SILOC,SGI)*LOCPNORMZ(ILOC)
                END DO
                RN=SQRT(SDIRX(SGI)**2+SDIRY(SGI)**2+SDIRZ(SGI)**2)
                SDIRX(SGI)=SDIRX(SGI)/RN
                SDIRY(SGI)=SDIRY(SGI)/RN
                SDIRZ(SGI)=SDIRZ(SGI)/RN
                NDOTQ(SGI)=SNORMXN(SGI)*SDIRX(SGI)+SNORMYN(SGI)*SDIRY(SGI) &
                     +SNORMZN(SGI)*SDIRZ(SGI)
                INCOME(SGI)=0.0
                IF(NDOTQ(SGI).LT.0.0) INCOME(SGI)=1.0
                IF(ELE2.EQ.0) THEN
                   IF(NDOTQ(SGI).LT.-0.8) TOPSUF=.TRUE.
                ENDIF
             END DO

             IF(TOPSUF) THEN

                ! calculate SNX,SNY,SNZ***
                DO SILOC=1,SNLOC
                   ILOC=SILOC2ILOC(SILOC)
                   INOD=XONDGL((ELE-1)*NLOC+ILOC)
                   DO SJLOC=SILOC,SNLOC,1
                      JLOC=SILOC2ILOC(SJLOC)
                      JNOD=XONDGL((ELE-1)*NLOC+JLOC)
                      KLOC=ILINK2(ILOC,JLOC)
                      SKLOC=SILINK2(SILOC,SJLOC)
                      XSL(SKLOC)=0.5*(XORIG(INOD)+XORIG(JNOD))
                      YSL(SKLOC)=0.5*(YORIG(INOD)+YORIG(JNOD))
                      ZSL(SKLOC)=0.5*(ZORIG(INOD)+ZORIG(JNOD))
                      ! Assume the mid side nodes are on the sphere
                      IF(ISPHERE.GE.2) THEN
                         IF(ILOC.NE.JLOC) THEN
                            RADI=SQRT(XORIG(INOD)**2+YORIG(INOD)**2+ZORIG(INOD)**2)
                            RADJ=SQRT(XORIG(JNOD)**2+YORIG(JNOD)**2+ZORIG(JNOD)**2)
                            RADP=0.5*(RADI+RADJ)
                            RAD=SQRT(XSL(SKLOC)**2+YSL(SKLOC)**2+ZSL(SKLOC)**2)
                            XSL(SKLOC)=XSL(SKLOC)*RADP/RAD
                            YSL(SKLOC)=YSL(SKLOC)*RADP/RAD
                            ZSL(SKLOC)=ZSL(SKLOC)*RADP/RAD
                         ENDIF
                      ENDIF
                   END DO
                END DO
                ! Now amend (XSL,YSL,ZSL) in the direction of the free surface
                ! change...
                DO SKLOC=1,SNCLOC
                   KLOC=MLOCSILOC2ILOC(SKLOC)
                   GLOBK =PGNDGLNO((ELE-1)*MLOC+KLOC)
                   XSL(SKLOC)=XSL(SKLOC) + LOCPNORMX(KLOC)*FREES(GLOBK)
                   YSL(SKLOC)=YSL(SKLOC) + LOCPNORMY(KLOC)*FREES(GLOBK)
                   ZSL(SKLOC)=ZSL(SKLOC) + LOCPNORMZ(KLOC)*FREES(GLOBK)
                END DO

                ! Recalculate the normal...
                CALL DGSDETNXLOC2(SNCLOC,SNGI, &
                     XSL,YSL,ZSL, &
                     SNC,SNCLX,SNCLY, SWEIGH, SDETWE,SAREA,.TRUE.,.FALSE., &
                     SNORMXN,SNORMYN,SNORMZN, &
                     NORMX,NORMY,NORMZ)

                UTOT=0.0
                VTOT=0.0
                WTOT=0.0
                DO SILOC=1,SNLOC
                   ILOC=SILOC2ILOC(SILOC)
                   IGL =NDGLNO((ELE-1)*NLOC+ILOC)
                   UTOT(ILOC)=U(IGL)
                   VTOT(ILOC)=V(IGL)
                   WTOT(ILOC)=W(IGL)
                   IF((INUSUB.NE.0).AND.(NSUBNVLOC.NE.0)) THEN
                      IDGNOD=(ELE-1)*NLOC+ILOC
                      UTOT(ILOC)=UTOT(ILOC) + NUSUB(IDGNOD)
                      VTOT(ILOC)=VTOT(ILOC) + NVSUB(IDGNOD)
                      WTOT(ILOC)=WTOT(ILOC) + NWSUB(IDGNOD)
                   ENDIF
                END DO

                sgi_loop: DO SGI=1,SNGI

                   AGI=0.
                   BGI=0.
                   CGI=0.

                   DGI=0.
                   EGI=0.
                   FGI=0.

                   XDGI=0.0
                   YDGI=0.0
                   ZDGI=0.0
                   DO SILOC=1,SNCLOC
                      ! NB R0 does not appear here although the z-coord might be Z+R0.
                      AGI=AGI+SNCLX(SILOC,SGI)*XSL(SILOC)
                      BGI=BGI+SNCLX(SILOC,SGI)*YSL(SILOC)
                      CGI=CGI+SNCLX(SILOC,SGI)*ZSL(SILOC)

                      DGI=DGI+SNCLY(SILOC,SGI)*XSL(SILOC)
                      EGI=EGI+SNCLY(SILOC,SGI)*YSL(SILOC)
                      FGI=FGI+SNCLY(SILOC,SGI)*ZSL(SILOC)

                      XDGI=XDGI+SNC(SILOC,SGI)*XSL(SILOC)
                      YDGI=YDGI+SNC(SILOC,SGI)*YSL(SILOC)
                      ZDGI=ZDGI+SNC(SILOC,SGI)*ZSL(SILOC)
                   END DO
                   IF(COGRAX) THEN
                      GGI=0.0
                      HGI=0.0
                      KGI=1.0
                      NORGRAVX=0.0
                      NORGRAVY=0.0
                      NORGRAVZ=1.0
                   ELSE
                      RGI=SQRT(XDGI**2+YDGI**2+ZDGI**2)
                      GGI=RGI/MAX(1.0E-10,SQRT(RGI**2-(YDGI**2+ZDGI**2)))
                      HGI=RGI/MAX(1.0E-10,SQRT(RGI**2-(XDGI**2+ZDGI**2)))
                      KGI=RGI/MAX(1.0E-10,SQRT(RGI**2-(XDGI**2+YDGI**2)))
                      NORGRAVX=XDGI/RGI
                      NORGRAVY=YDGI/RGI
                      NORGRAVZ=ZDGI/RGI
                   ENDIF
                   NGDOTNORM(SGI)=NORGRAVX*SNORMXN(SGI)+NORGRAVY*SNORMYN(SGI) &
                        +NORGRAVZ*SNORMZN(SGI)

                   CALL INV3X3( &
                        AGI,BGI,CGI,DGI,EGI,FGI,GGI,HGI,KGI, &
                        A11,A12,A13,A21,A22,A23,A31,A32,A33)

                   DO SILOC=1,SNLOC
                      SNX(SILOC,SGI)= A11*SNLX(SILOC,SGI)+A12*SNLY(SILOC,SGI)
                      SNY(SILOC,SGI)= A21*SNLX(SILOC,SGI)+A22*SNLY(SILOC,SGI)
                      SNZ(SILOC,SGI)= A31*SNLX(SILOC,SGI)+A32*SNLY(SILOC,SGI)
                   END DO
                   DO SILOC=1,SNCLOC
                      SNCX(SILOC,SGI)= A11*SNCLX(SILOC,SGI)+A12*SNCLY(SILOC,SGI)
                      SNCY(SILOC,SGI)= A21*SNCLX(SILOC,SGI)+A22*SNCLY(SILOC,SGI)
                      SNCZ(SILOC,SGI)= A31*SNCLX(SILOC,SGI)+A32*SNCLY(SILOC,SGI)
                   END DO
                   ! rotate and rotate back without vertical (E.G. r on
                   !  spherical Earth)  component:

                   RTREST11=1.0 - NORGRAVX*NORGRAVX
                   RTREST12=    - NORGRAVX*NORGRAVY
                   RTREST13=    - NORGRAVX*NORGRAVZ

                   RTREST21=    - NORGRAVY*NORGRAVX
                   RTREST22=1.0 - NORGRAVY*NORGRAVY
                   RTREST23=    - NORGRAVY*NORGRAVZ

                   RTREST31=    - NORGRAVZ*NORGRAVX
                   RTREST32=    - NORGRAVZ*NORGRAVY
                   RTREST33=1.0 - NORGRAVZ*NORGRAVZ
                   DO SILOC=1,SNLOC
                      TEMSNX=SNX(SILOC,SGI)
                      TEMSNY=SNY(SILOC,SGI)
                      TEMSNZ=SNZ(SILOC,SGI)
                      SNX(SILOC,SGI)= RTREST11*TEMSNX+RTREST12*TEMSNY+RTREST13*TEMSNZ
                      SNY(SILOC,SGI)= RTREST21*TEMSNX+RTREST22*TEMSNY+RTREST23*TEMSNZ
                      SNZ(SILOC,SGI)= RTREST31*TEMSNX+RTREST32*TEMSNY+RTREST33*TEMSNZ
                   END DO
                   DO SILOC=1,SNCLOC
                      TEMSNCX=SNCX(SILOC,SGI)
                      TEMSNCY=SNCY(SILOC,SGI)
                      TEMSNCZ=SNCZ(SILOC,SGI)
                      SNCX(SILOC,SGI)= RTREST11*TEMSNCX+RTREST12*TEMSNCY+RTREST13*TEMSNCZ
                      SNCY(SILOC,SGI)= RTREST21*TEMSNCX+RTREST22*TEMSNCY+RTREST23*TEMSNCZ
                      SNCZ(SILOC,SGI)= RTREST31*TEMSNCX+RTREST32*TEMSNCY+RTREST33*TEMSNCZ
                   END DO

                   SUD(SGI)=0.0
                   SVD(SGI)=0.0
                   SWD(SGI)=0.0
                   SABSGI(SGI)=0.0
                   DO SILOC=1,SNLOC
                      ILOC=SILOC2ILOC(SILOC)
                      IGL =NDGLNO((ELE-1)*NLOC+ILOC)
                      SUD(SGI)=SUD(SGI) + SN(SILOC,SGI)*U(IGL)
                      SVD(SGI)=SVD(SGI) + SN(SILOC,SGI)*V(IGL)
                      SWD(SGI)=SWD(SGI) + SN(SILOC,SGI)*W(IGL)
                      IF((INUSUB.NE.0).AND.(NSUBNVLOC.NE.0)) THEN
                         IDGNOD=(ELE-1)*NLOC+ILOC
                         SUD(SGI)=SUD(SGI) + SN(SILOC,SGI)*NUSUB(IDGNOD)
                         SVD(SGI)=SVD(SGI) + SN(SILOC,SGI)*NVSUB(IDGNOD)
                         SWD(SGI)=SWD(SGI) + SN(SILOC,SGI)*NWSUB(IDGNOD)
                      ENDIF
                      SABSGI(SGI)=SABSGI(SGI)+SN(SILOC,SGI)*ABSORB(IGL)
                   END DO
                   SSPEED(SGI)=SQRT(SUD(SGI)**2+SVD(SGI)**2+SWD(SGI)**2)
                   SDEPTH(SGI)=0.0
                   DO SILOC=1,SNLOC
                      ILOC=SILOC2ILOC(SILOC)
                      IGLX=XONDGL((ELE-1)*NLOC+ILOC)
                      SDEPTH(SGI)=SDEPTH(SGI)+SN(SILOC,SGI)*(TOPDIS(IGLX)+BOTDIS(IGLX))
                   END DO
                   GIFREES(SGI)=0.0
                   GIFREESOLD(SGI)=0.0
                   DO SJLOC=1,SNCLOC
                      JLOC=MLOCSILOC2ILOC(SJLOC)
                      GLOBJ =PGNDGLNO((ELE-1)*MLOC+JLOC)
                      GIFREES(SGI)=GIFREES(SGI)+       SNC(SJLOC,SGI)*FREES(GLOBJ)
                      GIFREESOLD(SGI)=GIFREESOLD(SGI)+ SNC(SJLOC,SGI)*FREESOLD(GLOBJ)
                   END DO

                   OME2=2.*FUNOME(XDGI,YDGI,ZDGI)
                   IF((GEOBAL.EQ.-12).or.(GEOBAL.EQ.-22)) OME2=0.
                   SA11(SGI)=1.0         + DT*ATHETA*SABSGI(SGI)
                   SA12(SGI)=-DT*CTHETA*OME2
                   SA13(SGI)=0.0
                   SA21(SGI)=DT*CTHETA*OME2
                   SA22(SGI)=1.0         + DT*ATHETA*SABSGI(SGI)
                   SA23(SGI)=0.0
                   SA31(SGI)=0.0
                   SA32(SGI)=0.0
                   SA33(SGI)=1.0         + DT*ATHETA*SABSGI(SGI)
                   CALL INV3X3( &
                        SA11(SGI),SA12(SGI),SA13(SGI),SA21(SGI),SA22(SGI),SA23(SGI),SA31(SGI),SA32(SGI),SA33(SGI), &
                        SI11(SGI),SI12(SGI),SI13(SGI),SI21(SGI),SI22(SGI),SI23(SGI),SI31(SGI),SI32(SGI),SI33(SGI))

                   SQRTDISPX(SGI)=SQRT(FSTHETA*SI11(SGI)*SDEPTH(SGI)*GRAVTY)
                   SQRTDISPY(SGI)=SQRT(FSTHETA*SI22(SGI)*SDEPTH(SGI)*GRAVTY)
                   SQRTDISPZ(SGI)=SQRT(FSTHETA*SI33(SGI)*SDEPTH(SGI)*GRAVTY)
                END DO sgi_loop

                ! Calculate DX - the characteristic horizontal element length scale
                ! IF NOFSTA=1 switch off stabilisation =0 include it.
                DX=0.0
                DO SILOC=1,SNLOC
                   ILOC=SILOC2ILOC(SILOC)
                   IGLX=XONDGL((ELE-1)*NLOC+ILOC)
                   DO SJLOC=SILOC+1,SNLOC,1
                      JLOC=SILOC2ILOC(SJLOC)
                      JGLX=XONDGL((ELE-1)*NLOC+JLOC)
                      DX=MAX(DX, ABS(XORIG(IGLX)-XORIG(JGLX)), ABS(YORIG(IGLX)-YORIG(JGLX)), &
                           ABS(ZORIG(IGLX)-ZORIG(JGLX)) )
                   END DO
                END DO
                IF(NOFSTA.EQ.0) THEN
                   DO SGI=1,SNGI
                      SKSTAB(SGI)=(  (dx*sqrt(gravty*Sdepth(Sgi))/dt)   )/(1.+DT*ATHETA*SABSGI(SGI))
                      IF(NONLIN) THEN
                         RESID=-(GIFREES(SGI)-GIFREESOLD(SGI))/DT &
                              +(SNORMXN(SGI)*SUD(SGI)+SNORMYN(SGI)*SVD(SGI)+SNORMZN(SGI)*SWD(SGI))/NGDOTNORM(SGI)
                         rnonl=rnonl_coefficient*abs(resid)/ &
                              max(abs((snormxn(sgi)*sud(sgi)+snormyn(sgi)*svd(sgi)+snormzn(sgi)*swd(sgi)) &
                              /ngdotnorm(sgi)),1.e-10)
                         SKSTAB(SGI)=SKSTAB(SGI)*MIN(1.0,RNONL)
                      ENDIF
                      SQRTSTA(SGI)=SQRT(SKSTAB(SGI))
                      DXSQRTSTA(SGI)=DX*SQRTSTA(SGI)
                   END DO
                ELSE
                   SKSTAB(1:SNGI) = 0.0
                   SQRTSTA(1:SNGI) = 0.0
                   DXSQRTSTA(1:SNGI) = 0.0
                ENDIF

                ! Perform surface integration...
                siloc_loop: DO SILOC=1,SNLOC
                   ILOC=SILOC2ILOC(SILOC)
                   GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
                   vec(globi)=1.0
                   sjloc_loop: DO SJLOC=1,SNLOC
                      JLOC=SILOC2ILOC(SJLOC)
                      GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)

                      OUTVLMX=0.0
                      OUTVLMY=0.0
                      OUTVLMZ=0.0
                      NNINVD=0.0

                      KMATH=0.0
                      KKS=0.0

                      NNMLX=0.0
                      NNMLY=0.0
                      NNMLZ=0.0

                      STANNML=0.0

                      SNSNX=0.0
                      SNSNY=0.0
                      SNSNZ=0.0

                      SNSNXX=0.0
                      SNSNYY=0.0
                      SNSNZZ=0.0

                      second_sgi_loop: DO SGI=1,SNGI
                         RNN=SN(SILOC,SGI)*SN(SJLOC,SGI)*SDETWE(SGI)
                         OUTVLMX=OUTVLMX+RNN*SNORMXN(SGI)/NGDOTNORM(SGI)
                         OUTVLMY=OUTVLMY+RNN*SNORMYN(SGI)/NGDOTNORM(SGI)
                         OUTVLMZ=OUTVLMZ+RNN*SNORMZN(SGI)/NGDOTNORM(SGI)
                         NNINVD=NNINVD+1.0*RNN
                         NNMLX=NNMLX+RNN
                         NNMLY=NNMLY+RNN
                         NNMLZ=NNMLZ+RNN

                         STANNML=STANNML+RNN

                         IF(.NOT.SMALLSTEN) THEN
                            SNSNX=SNSNX+ SQRTDISPX(SGI)*SDETWE(SGI)*SNX(SILOC,SGI)*SN(SJLOC,SGI)
                            SNSNY=SNSNY+ SQRTDISPY(SGI)*SDETWE(SGI)*SNY(SILOC,SGI)*SN(SJLOC,SGI)
                            SNSNZ=SNSNZ+ SQRTDISPZ(SGI)*SDETWE(SGI)*SNZ(SILOC,SGI)*SN(SJLOC,SGI)
                         ENDIF

                         SNSNXX=SNSNXX+ SQRTSTA(SGI)*SDETWE(SGI)*SNX(SILOC,SGI)*SN(SJLOC,SGI)
                         SNSNYY=SNSNYY+ SQRTSTA(SGI)*SDETWE(SGI)*SNY(SILOC,SGI)*SN(SJLOC,SGI)
                         SNSNZZ=SNSNZZ+ SQRTSTA(SGI)*SDETWE(SGI)*SNZ(SILOC,SGI)*SN(SJLOC,SGI)

                         IF(SMALLSTEN) THEN
                            KMATH=KMATH+FSTHETA*wdgravity(ele)*GRAVTY*SDETWE(SGI)*SDEPTH(SGI)*( &
                                 SNX(SILOC,SGI)*(SI11(SGI)*SNX(SJLOC,SGI)+SI12(SGI)*SNY(SJLOC,SGI)+SI13(SGI)*SNZ(SJLOC,SGI)) &
                                 + SNY(SILOC,SGI)*(SI21(SGI)*SNX(SJLOC,SGI)&
                                 &+SI22(SGI)*SNY(SJLOC,SGI)+SI23(SGI)*SNZ(SJLOC,SGI)) &
                                 + SNZ(SILOC,SGI)*(SI31(SGI)*SNX(SJLOC,SGI)&
                                 &+SI32(SGI)*SNY(SJLOC,SGI)+SI33(SGI)*SNZ(SJLOC,SGI)))
                         ENDIF
                         IF(NOFSTA.EQ.0) THEN
                            KKS=KKS+SKSTAB(SGI)*SDETWE(SGI)*(SNX(SILOC,SGI)*SNX(SJLOC,SGI) &
                                 +SNY(SILOC,SGI)*SNY(SJLOC,SGI)+SNZ(SILOC,SGI)*SNZ(SJLOC,SGI))
                         ENDIF
                      END DO second_sgi_loop
                      ! Find POSMAT. ***************************************
                      CALL POSINMAT(POSMAT,GLOBI,GLOBJ, &
                           FREDOP,FINCMC,COLCMC,NCMC)
                      CALL POSINMAT(POSCT,GLOBI,GLOBJ, &
                           NONODS,FINDRM,COLM,NCOLM)
                      ! Found POSMAT ***************************************
                      CXT(POSCT)=CXT(POSCT)+SNSNX
                      CYT(POSCT)=CYT(POSCT)+SNSNY
                      CZT(POSCT)=CZT(POSCT)+SNSNZ

                      CXXT(POSCT)=CXXT(POSCT)+SNSNXX
                      CYYT(POSCT)=CYYT(POSCT)+SNSNYY
                      CZZT(POSCT)=CZZT(POSCT)+SNSNZZ

                      MLX(GLOBI)=MLX(GLOBI)+NNMLX
                      MLY(GLOBI)=MLY(GLOBI)+NNMLY
                      MLZ(GLOBI)=MLZ(GLOBI)+NNMLZ

                      STAML(GLOBI)=STAML(GLOBI)+STANNML 
                      SUFML(GLOBI)=SUFML(GLOBI)+NNINVD

                      LUMMAT=.false.
                      LUMMAT2=.false.
                      LUMMAT3=.false.
                      IF(WETDRY.EQ.1) THEN
                         WETDRYMAT(POSMAT)=WETDRYMAT(POSMAT)+NNINVD
                      ENDIF
                      ICENT=MIDCMC(GLOBI)
                      IF(VELMEAN.EQ.1) THEN
                         IF(LUMMAT) THEN
                            SHORTDIFMATH(ICENT)=SHORTDIFMATH(ICENT)+NNINVD/(DT**2)
                            SHORTDIFMATH(POSMAT)=SHORTDIFMATH(POSMAT)+KMATH
                         ELSE
                            SHORTDIFMATH(POSMAT)=SHORTDIFMATH(POSMAT)+KMATH +NNINVD/(DT**2)
                         ENDIF
                      ELSE
                         IF(LUMMAT) THEN
                            SHORTDIFMATH(ICENT)=SHORTDIFMATH(ICENT)+NNINVD/(DT**2)
                            SHORTDIFMATH(POSMAT)=SHORTDIFMATH(POSMAT)+KMATH + KKS
                         ELSE
                            SHORTDIFMATH(POSMAT)=SHORTDIFMATH(POSMAT)+KMATH +NNINVD/(DT**2) + KKS
                         ENDIF
                      ENDIF
                      ! surface flux term:
                      IF(LUMMAT3) THEN
                         SHORTDIFVECX(GLOBI)=SHORTDIFVECX(GLOBI) &
                              +(OUTVLMX*UTOT(ILOC)+OUTVLMY*VTOT(ILOC)+OUTVLMZ*WTOT(ILOC))/DT
                      ELSE
                         SHORTDIFVECX(GLOBI)=SHORTDIFVECX(GLOBI) &
                              +(OUTVLMX*UTOT(JLOC)+OUTVLMY*VTOT(JLOC)+OUTVLMZ*WTOT(JLOC))/DT
                         conout=conout+OUTVLMX*UTOT(JLOC)+OUTVLMY*VTOT(JLOC)+OUTVLMZ*WTOT(JLOC)
                      ENDIF

                   END DO sjloc_loop
                   ! Now for the vector....
                   DO SJLOC=1,SNCLOC
                      JLOC=MLOCSILOC2ILOC(SJLOC)
                      GLOBJ =PGNDGLNO((ELE-1)*MLOC+JLOC)

                      NNINVD=0.0
                      KKS=0.0
                      DO SGI=1,SNGI
                         NNINVD=NNINVD+SN(SILOC,SGI)*SNC(SJLOC,SGI)*SDETWE(SGI)
                         IF(NOFSTA.EQ.0) THEN
                            KKS=KKS+SKSTAB(SGI)*SDETWE(SGI)*(SNX(SILOC,SGI)*SNCX(SJLOC,SGI) &
                                 +SNY(SILOC,SGI)*SNCY(SJLOC,SGI)+SNZ(SILOC,SGI)*SNCZ(SJLOC,SGI))
                         ENDIF
                      END DO
                      MSUFML(GLOBJ)=MSUFML(GLOBJ)+NNINVD
                      ! The matrix eqn associated with velocity star
                      IF(VELMEAN.EQ.1) THEN
                         IF(LUMMAT2) THEN
                            SHORTDIFVECX(GLOBI)=SHORTDIFVECX(GLOBI)-(NNINVD/(DT**2))*(FREES(GLOBI)-FREESOLD(GLOBI))
                         ELSE
                            SHORTDIFVECX(GLOBI)=SHORTDIFVECX(GLOBI)-(NNINVD/(DT**2))*(FREES(GLOBJ)-FREESOLD(GLOBJ))
                         ENDIF
                      ELSE
                         IF(LUMMAT2) THEN
                            SHORTDIFVECX(GLOBI)=SHORTDIFVECX(GLOBI)-(NNINVD/(DT**2))*(FREES(GLOBI)-FREESOLD(GLOBI))
                         ELSE
                            SHORTDIFVECX(GLOBI)=SHORTDIFVECX(GLOBI)-(NNINVD/(DT**2))*(FREES(GLOBJ)-FREESOLD(GLOBJ))
                         ENDIF
                         SHORTDIFVECX(GLOBI)=SHORTDIFVECX(GLOBI)- KKS*FREES(GLOBJ)
                      ENDIF
                   END DO

                END DO siloc_loop
             ENDIF
          ENDIF

       END DO face_loop

    END DO element_loop

    ! Conservation check...
    CONSUM=0.0
    SAREA=0.0
    DO GLOBI=1,PGNODS
       CONSUM=CONSUM+MSUFML(GLOBI)*FREES(GLOBI)
       SAREA=SAREA+MSUFML(GLOBI)
    END DO
    EWRITE(1,*) 'Conservation check surface area*** (m^2) SAREA=',SAREA
    SAREA=0.0
    DO GLOBI=1,NONODS
       SAREA=SAREA+SUFML(GLOBI)
    END DO
    EWRITE(1,*) 'Conservation check surface area*** (m^2) SAREA=',SAREA
    EWRITE(1,*) 'Conservation check *** (m^3) CONSUM=',CONSUM
    EWRITE(1,*) 'This is the averaged error in conservation across the surface'
    EWRITE(1,*) 'Conservation check *** (m^3) CONSUM=',CONSUM
    EWRITE(1,*) 'CONSUM should not change and may be zero if'
    EWRITE(1,*) 'simulation started with a flat free surface and no inlet'  
    EWRITE(1,*) 'CONSUM is the volume error of the free surface'  
    EWRITE(1,*) 'Conservation check of the source integral *** (m^3/S) CONOUT=',CONOUT   
    EWRITE(1,*) 'Conservation check of the source integral *** (m/S) CONOUT/SAREA=',CONOUT/SAREA  

  END SUBROUTINE ALLTOPSUFDG

  SUBROUTINE CTDGSUP(C1T,C2T,C3T, &
       XORIG,YORIG,ZORIG, &
       N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI, &
       TOTELE,D3, &
       NDGLNO,XONDGL, &
       NONODS,XNONOD, &
       ISPHERE, &
       FINDCT,COLCT,NCT, &
       COGRAX, &
       NOBCFSimp,BCT2FSimp, &
       PGNODS,NCLOC2,FREES,PGNDGLNO, &
       TOPDIS,BOTDIS,COLELE4)
    ! This subroutine caluculates C1T,C2T,C3T by NOT integrating
    ! the pressure term by parts and integrating over the
    ! top surface of the ocean.
    ! M is the basis function assocaited with normal pressure.
    use FLDebug
    use assnav_module
    IMPLICIT NONE
    ! Number of source nodes, advection nodes and node positions.
    ! Number of elements
    LOGICAL SLUM
    PARAMETER(SLUM=.FALSE.)
    ! IF SLUM then lump the mass matrix in the surface integral.
    INTEGER TOTELE,NONODS,XNONOD
    INTEGER NGI,NLOC
    LOGICAL D3, COGRAX

    REAL XORIG(XNONOD),YORIG(XNONOD),ZORIG(XNONOD)
    ! Global node number maps.
    INTEGER NDGLNO(TOTELE*NLOC)
    INTEGER XONDGL(TOTELE*NLOC)
    INTEGER NCT
    INTEGER FINDCT(NONODS+1),COLCT(NCT)
    REAL C1T(NCT),C2T(NCT),C3T(NCT)

    INTEGER ISPHERE

    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL WEIGHT(NGI)
    ! boundary conditions for the free surface where we want to avoid
    ! applying n\cdot u=0 b.c.
    INTEGER NOBCFSimp
    INTEGER BCT2FSimp(NOBCFSimp)
    INTEGER PGNODS,NCLOC2,PGNDGLNO(TOTELE*NCLOC2)
    REAL FREES(PGNODS)
    REAL TOPDIS(XNONOD),BOTDIS(XNONOD)
    INTEGER COLELE4(TOTELE*5)

    ! Local variables...
    INTEGER SNLOC,SNGI,SNCLOC,NCLOC
    PARAMETER(SNLOC=3,SNCLOC=6,SNGI=7,NCLOC=10)
    REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI),SNLZ(SNLOC,SNGI)
    REAL SWEIGH(SNGI)
    REAL SNC(SNCLOC,SNGI),SNCLX(SNCLOC,SNGI),SNCLY(SNCLOC,SNGI),SNCLZ(SNCLOC,SNGI)
    REAL NC(NCLOC,NGI),NCLX(NCLOC,NGI),NCLY(NCLOC,NGI),NCLZ(NCLOC,NGI)
    REAL SL1(SNGI), SL2(SNGI), SL3(SNGI), SL4(SNGI)
    REAL L1(NGI), L2(NGI), L3(NGI), L4(NGI)
    REAL DETWEI(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    REAL NCX(NCLOC,NGI),NCY(NCLOC,NGI),NCZ(NCLOC,NGI)
    REAL SDETWE(SNGI)
    REAL PGWEIG(NGI)
    REAL A11(NGI),A12(NGI),A13(NGI)
    REAL A21(NGI),A22(NGI),A23(NGI)
    REAL A31(NGI),A32(NGI),A33(NGI)
    !
    REAL SNORMXN(SNGI),SNORMYN(SNGI),SNORMZN(SNGI)
    REAL SDIRX(SNGI),SDIRY(SNGI),SDIRZ(SNGI)
    REAL LOCPNORMX(NCLOC),LOCPNORMY(NCLOC),LOCPNORMZ(NCLOC)
    REAL NORSCALE(NCLOC)
    REAL XL(NCLOC),YL(NCLOC),ZL(NCLOC)
    REAL XSL(SNCLOC),YSL(SNCLOC),ZSL(SNCLOC)

    INTEGER SILOC2ILOC(SNLOC),OTHILOC(NLOC)
    INTEGER MLOCSILOC2ILOC(SNCLOC),MLOCOTHILOC(NCLOC)
    INTEGER ELE,ILOC,JLOC,SKLOC
    INTEGER GI,SGI,GLOBI,GLOBJ,GLOBJ2,GLOBK
    INTEGER COUNT2

    REAL NDOTQ(SNGI),INCOME(SNGI)
    REAL SAREA
    INTEGER ELE2
    INTEGER SILOC,SJLOC
    INTEGER KLOC,INOD,JNOD
    INTEGER POSMAT,COUNT
    REAL RNN,OUTVLMX,OUTVLMY,OUTVLMZ,RN
    REAL NXN,NYN,NZN

    REAL NORMX,NORMY,NORMZ
    REAL RADI,RADJ,RADP,RAD
    LOGICAL TOPSUF,BNDSUF,LUMMAT
    INTEGER ICENT,ICOUNT,II
    INTEGER IFACE,NFACE
    PARAMETER(NFACE=4)
    INTEGER LOCLIST(NFACE,3)

    REAL, ALLOCATABLE, DIMENSION(:)::vec
    INTEGER, ALLOCATABLE, DIMENSION(:)::FINELE
    INTEGER, ALLOCATABLE, DIMENSION(:)::BNDCON

    ALLOCATE(vec(nonods))
    ALLOCATE(FINELE(TOTELE+1))

    ! BNDCON marks the boundary condition nodes...
    ALLOCATE(BNDCON(NONODS))

    BNDCON(1:NONODS)=0
    DO II=1,NOBCFSimp
       GLOBI=BCT2FSimp(II)
       BNDCON(GLOBI)=1
    END DO

    ! Calculate surface shape functions SN, SM etc...
    CALL TRIQUA(SL1, SL2, SL3, SL3, SWEIGH, .FALSE.,SNGI)
    CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
         SNLOC,SNGI, &
         SN,SNLX,SNLY,SNLZ)
    CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
         SNCLOC,SNGI, &
         SNC,SNCLX,SNCLY,SNCLZ)

    ! Work out the shape functions and there derivatives...
    CALL TRIQUA(L1, L2, L3, L4, PGWEIG, D3,NGI)
    CALL SHATRI(L1, L2, L3, L4, PGWEIG, D3, &
         NCLOC,NGI,  &
         NC,NCLX,NCLY,NCLZ)

    C1T(1:NCT)=0.0
    C2T(1:NCT)=0.0
    C3T(1:NCT)=0.0

    first_element_loop: DO ELE=1,TOTELE

       ! This is for a surface element...
       IF(COGRAX) THEN
          DO KLOC=1,NCLOC
             LOCPNORMX(KLOC)=0.0
             LOCPNORMY(KLOC)=0.0
             LOCPNORMZ(KLOC)=1.0
          END DO
          DO ILOC=1,NLOC
             INOD=XONDGL((ELE-1)*NLOC+ILOC)
             DO JLOC=ILOC,NLOC
                JNOD=XONDGL((ELE-1)*NLOC+JLOC)
                KLOC=ILINK2(ILOC,JLOC)
                NORSCALE(KLOC)=0.5*(BOTDIS(INOD)+BOTDIS(JNOD)) &
                     /(0.5*(BOTDIS(INOD)+BOTDIS(JNOD)+TOPDIS(INOD)+TOPDIS(JNOD)))
             END DO
          END DO
       ELSE
          ! Assume its a spherical Earth...
          DO ILOC=1,NLOC
             INOD=XONDGL((ELE-1)*NLOC+ILOC)
             DO JLOC=ILOC,NLOC
                JNOD=XONDGL((ELE-1)*NLOC+JLOC)
                KLOC=ILINK2(ILOC,JLOC)
                LOCPNORMX(KLOC)=0.5*(XORIG(INOD)+XORIG(JNOD))
                LOCPNORMY(KLOC)=0.5*(YORIG(INOD)+YORIG(JNOD))
                LOCPNORMZ(KLOC)=0.5*(ZORIG(INOD)+ZORIG(JNOD))
                NORSCALE(KLOC)=0.5*(BOTDIS(INOD)+BOTDIS(JNOD)) &
                     /(0.5*(BOTDIS(INOD)+BOTDIS(JNOD)+TOPDIS(INOD)+TOPDIS(JNOD)))
             END DO
          END DO
          DO KLOC=1,NCLOC
             RN=SQRT(LOCPNORMX(KLOC)**2+LOCPNORMY(KLOC)**2+LOCPNORMZ(KLOC)**2)
             LOCPNORMX(KLOC)=LOCPNORMX(KLOC)/RN
             LOCPNORMY(KLOC)=LOCPNORMY(KLOC)/RN
             LOCPNORMZ(KLOC)=LOCPNORMZ(KLOC)/RN
          END DO
       ENDIF
       DO ILOC=1,NLOC
          INOD=XONDGL((ELE-1)*NLOC+ILOC)
          DO JLOC=ILOC,NLOC,1
             JNOD=XONDGL((ELE-1)*NLOC+JLOC)
             KLOC=ILINK2(ILOC,JLOC)
             XL(KLOC)=0.5*(XORIG(INOD)+XORIG(JNOD))
             YL(KLOC)=0.5*(YORIG(INOD)+YORIG(JNOD))
             ZL(KLOC)=0.5*(ZORIG(INOD)+ZORIG(JNOD))
             ! Assume the mid side nodes are on the sphere
             IF(ISPHERE.GE.2) THEN
                IF(ILOC.NE.JLOC) THEN
                   RADI=SQRT(XORIG(INOD)**2+YORIG(INOD)**2+ZORIG(INOD)**2)
                   RADJ=SQRT(XORIG(JNOD)**2+YORIG(JNOD)**2+ZORIG(JNOD)**2)
                   RADP=0.5*(RADI+RADJ)
                   RAD=SQRT(XL(KLOC)**2+YL(KLOC)**2+ZL(KLOC)**2)
                   XL(KLOC)=XL(KLOC)*RADP/RAD
                   YL(KLOC)=YL(KLOC)*RADP/RAD
                   ZL(KLOC)=ZL(KLOC)*RADP/RAD
                ENDIF
             ENDIF
          END DO
       END DO
       ! Now amend (XL,YL,ZL) in the direction of the free surface
       ! change...
       DO KLOC=1,NCLOC
          GLOBK =PGNDGLNO((ELE-1)*NCLOC+KLOC)
          XL(KLOC)=XL(KLOC) + LOCPNORMX(KLOC)*NORSCALE(KLOC)*FREES(GLOBK)
          YL(KLOC)=YL(KLOC) + LOCPNORMY(KLOC)*NORSCALE(KLOC)*FREES(GLOBK)
          ZL(KLOC)=ZL(KLOC) + LOCPNORMZ(KLOC)*NORSCALE(KLOC)*FREES(GLOBK)
       END DO

       ! Get determinant and derivatives...
       CALL DETNLXMSUP(NLOC,NCLOC,NGI, &
            NLX,NLY,NLZ, NCLX,NCLY,NCLZ, WEIGHT, DETWEI, &
            NX,NY,NZ, NCX,NCY,NCZ, &
            NCLOC,NCLX,NCLY,NCLZ,XL,YL,ZL, &
            A11,A12,A13, A21,A22,A23, A31,A32,A33)

       DO ILOC=1,NLOC
          GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
          DO JLOC=1,NLOC
             GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)
             NXN=0.0
             NYN=0.0
             NZN=0.0
             DO GI=1,NGI
                NXN=NXN+NX(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
                NYN=NYN+NY(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
                NZN=NZN+NZ(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
             END DO
             ! Find POSMAT. ***************************************
             CALL POSINMAT(POSMAT,GLOBI,GLOBJ, &
                  NONODS,FINDCT,COLCT,NCT)
             ! Found POSMAT ***************************************
             C1T(POSMAT)=C1T(POSMAT)-NXN
             C2T(POSMAT)=C2T(POSMAT)-NYN
             C3T(POSMAT)=C3T(POSMAT)-NZN
          END DO
       END DO

    END DO first_element_loop

    ! The local cords are in anti-clockwise order
    IFACE=1
    LOCLIST(IFACE,1)=1
    LOCLIST(IFACE,2)=2
    LOCLIST(IFACE,3)=3
    IFACE=2
    LOCLIST(IFACE,1)=1
    LOCLIST(IFACE,2)=2
    LOCLIST(IFACE,3)=4
    IFACE=3
    LOCLIST(IFACE,1)=3
    LOCLIST(IFACE,2)=2
    LOCLIST(IFACE,3)=4
    IFACE=4
    LOCLIST(IFACE,1)=1
    LOCLIST(IFACE,2)=3
    LOCLIST(IFACE,3)=4
    vec(1:nonods)=0.0

    second_element_loop:DO ELE=1,TOTELE

       ! Put surface integrals into matrix and rhs*****************************
       face_loop: DO IFACE=1,4
          COUNT2=(ELE-1)*5+IFACE
          ! NB colele4 is ordered in terms of faces.
          ELE2=COLELE4(COUNT2)

          IF(ELE2.EQ.0) THEN

             ! This is for a surface element...
             IF(COGRAX) THEN
                DO KLOC=1,NCLOC
                   LOCPNORMX(KLOC)=0.0
                   LOCPNORMY(KLOC)=0.0
                   LOCPNORMZ(KLOC)=1.0
                END DO
                DO ILOC=1,NLOC
                   INOD=XONDGL((ELE-1)*NLOC+ILOC)
                   DO JLOC=ILOC,NLOC
                      JNOD=XONDGL((ELE-1)*NLOC+JLOC)
                      KLOC=ILINK2(ILOC,JLOC)
                      NORSCALE(KLOC)=0.5*(BOTDIS(INOD)+BOTDIS(JNOD)) &
                           /(0.5*(BOTDIS(INOD)+BOTDIS(JNOD)+TOPDIS(INOD)&
                           &+TOPDIS(JNOD)))
                   ENDDO
                ENDDO
             ELSE
                ! Assume its a spherical Earth...
                DO ILOC=1,NLOC
                   INOD=XONDGL((ELE-1)*NLOC+ILOC)
                   DO JLOC=ILOC,NLOC
                      JNOD=XONDGL((ELE-1)*NLOC+JLOC)
                      KLOC=ILINK2(ILOC,JLOC)
                      LOCPNORMX(KLOC)=0.5*(XORIG(INOD)+XORIG(JNOD))
                      LOCPNORMY(KLOC)=0.5*(YORIG(INOD)+YORIG(JNOD))
                      LOCPNORMZ(KLOC)=0.5*(ZORIG(INOD)+ZORIG(JNOD))
                      NORSCALE(KLOC)=0.5*(BOTDIS(INOD)+BOTDIS(JNOD)) &
                           /(0.5*(BOTDIS(INOD)+BOTDIS(JNOD)+TOPDIS(INOD)&
                           &+TOPDIS(JNOD)))
                   END DO
                END DO
                DO KLOC=1,NCLOC
                   RN=SQRT(LOCPNORMX(KLOC)**2+LOCPNORMY(KLOC)**2+LOCPNORMZ(KLOC)**2)
                   LOCPNORMX(KLOC)=LOCPNORMX(KLOC)/RN
                   LOCPNORMY(KLOC)=LOCPNORMY(KLOC)/RN
                   LOCPNORMZ(KLOC)=LOCPNORMZ(KLOC)/RN
                END DO
             ENDIF

             ! Surface element of domain...
             DO SILOC=1,SNLOC
                SILOC2ILOC(SILOC)=LOCLIST(IFACE,SILOC)
             END DO
             ! Form OTHILOC from SILOC2ILOC
             IF(ELE2.NE.0) THEN
                CALL ELE2ELEFACE(OTHILOC,SILOC2ILOC,SNLOC,NLOC,TOTELE, &
                     XONDGL,NDGLNO,ELE,ELE2)
             ELSE
                DO ILOC=1,NLOC
                   OTHILOC(ILOC)=ILOC
                END DO
             ENDIF

             ! Form MLOCOTHILOC and MLOCSILOC2ILOC
             CALL MLOCELE2ELEFACE(MLOCOTHILOC,MLOCSILOC2ILOC,SNLOC,NLOC, &
                  SNCLOC,NCLOC,OTHILOC,SILOC2ILOC)

             ! PERFORM SURFACE INTEGRATION         !
             ! Form approximate surface normal (NORMX,NORMY,NORMZ)
             CALL DGSIMPLNORM(ELE,SILOC2ILOC,TOTELE,NLOC,SNLOC,XONDGL, &
                  XORIG,YORIG,ZORIG,XNONOD,NORMX,NORMY,NORMZ)

             CALL DGSDETNXLOC(ELE,SILOC2ILOC,XONDGL, &
                  TOTELE,NLOC,XNONOD,SNLOC,SNGI, &
                  XORIG,YORIG,ZORIG, &
                  SN,SNLX,SNLY, SWEIGH, SDETWE,SAREA, .TRUE.,.FALSE., &
                  SNORMXN,SNORMYN,SNORMZN, &
                  NORMX,NORMY,NORMZ, &
                  ISPHERE,SNCLOC,SNC,SNCLX,SNCLY)

             TOPSUF=.FALSE.
             DO SGI=1,SNGI
                SDIRX(SGI)=0.0
                SDIRY(SGI)=0.0
                SDIRZ(SGI)=0.0
                DO SILOC=1,SNCLOC
                   ILOC=MLOCSILOC2ILOC(SILOC)
                   SDIRX(SGI)=SDIRX(SGI)-SNC(SILOC,SGI)*LOCPNORMX(ILOC)
                   SDIRY(SGI)=SDIRY(SGI)-SNC(SILOC,SGI)*LOCPNORMY(ILOC)
                   SDIRZ(SGI)=SDIRZ(SGI)-SNC(SILOC,SGI)*LOCPNORMZ(ILOC)
                END DO
                RN=SQRT(SDIRX(SGI)**2+SDIRY(SGI)**2+SDIRZ(SGI)**2)
                SDIRX(SGI)=SDIRX(SGI)/RN
                SDIRY(SGI)=SDIRY(SGI)/RN
                SDIRZ(SGI)=SDIRZ(SGI)/RN
                NDOTQ(SGI)=SNORMXN(SGI)*SDIRX(SGI)+SNORMYN(SGI)*SDIRY(SGI) &
                     +SNORMZN(SGI)*SDIRZ(SGI)
                INCOME(SGI)=0.0
                IF(NDOTQ(SGI).LT.0.0) INCOME(SGI)=1.0
                IF(ELE2.EQ.0) THEN
                   IF(NDOTQ(SGI).LT.-0.8) TOPSUF=.TRUE.
                ENDIF
             END DO
             ! BNDSUF=.true. if we are on a boundary condition surface i.e. inlet-outlet
             ! for the free surface.
             ICOUNT=0
             DO SILOC=1,SNLOC
                ILOC=SILOC2ILOC(SILOC)
                GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
                IF(BNDCON(GLOBI).EQ.1) ICOUNT=ICOUNT+1
             ENDDO
             BNDSUF=(ICOUNT.EQ.SNLOC)

             IF(TOPSUF.OR.BNDSUF) THEN

                DO SILOC=1,SNLOC
                   ILOC=SILOC2ILOC(SILOC)
                   INOD=XONDGL((ELE-1)*NLOC+ILOC)
                   DO SJLOC=SILOC,SNLOC,1
                      JLOC=SILOC2ILOC(SJLOC)
                      JNOD=XONDGL((ELE-1)*NLOC+JLOC)
                      KLOC=ILINK2(ILOC,JLOC)
                      SKLOC=SILINK2(SILOC,SJLOC)
                      XSL(SKLOC)=0.5*(XORIG(INOD)+XORIG(JNOD))
                      YSL(SKLOC)=0.5*(YORIG(INOD)+YORIG(JNOD))
                      ZSL(SKLOC)=0.5*(ZORIG(INOD)+ZORIG(JNOD))
                      ! Assume the mid side nodes are on the sphere
                      IF(ISPHERE.GE.2) THEN
                         IF(ILOC.NE.JLOC) THEN
                            RADI=SQRT(XORIG(INOD)**2+YORIG(INOD)**2+ZORIG(INOD)**2)
                            RADJ=SQRT(XORIG(JNOD)**2+YORIG(JNOD)**2+ZORIG(JNOD)**2)
                            RADP=0.5*(RADI+RADJ)
                            RAD=SQRT(XSL(SKLOC)**2+YSL(SKLOC)**2+ZSL(SKLOC)**2)
                            XSL(SKLOC)=XSL(SKLOC)*RADP/RAD
                            YSL(SKLOC)=YSL(SKLOC)*RADP/RAD
                            ZSL(SKLOC)=ZSL(SKLOC)*RADP/RAD
                         ENDIF
                      ENDIF
                   END DO
                END DO
                ! Now amend (XSL,YSL,ZSL) in the direction of the free surface
                ! change...
                DO SKLOC=1,SNCLOC
                   KLOC=MLOCSILOC2ILOC(SKLOC)
                   GLOBK =PGNDGLNO((ELE-1)*NCLOC+KLOC)
                   XSL(SKLOC)=XSL(SKLOC) + LOCPNORMX(KLOC)*NORSCALE(KLOC)*FREES(GLOBK)
                   YSL(SKLOC)=YSL(SKLOC) + LOCPNORMY(KLOC)*NORSCALE(KLOC)*FREES(GLOBK)
                   ZSL(SKLOC)=ZSL(SKLOC) + LOCPNORMZ(KLOC)*NORSCALE(KLOC)*FREES(GLOBK)
                END DO

                ! Recalculate the normal...
                CALL DGSDETNXLOC2(SNCLOC,SNGI, &
                     XSL,YSL,ZSL, &
                     SNC,SNCLX,SNCLY, SWEIGH, SDETWE,SAREA,.TRUE.,.FALSE., &
                     SNORMXN,SNORMYN,SNORMZN, &
                     NORMX,NORMY,NORMZ)
                ! ********NEW INSERT*********
                ! Perform surface integration...
                DO SILOC=1,SNLOC
                   ILOC=SILOC2ILOC(SILOC)
                   GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
                   ! Find ICENT
                   LUMMAT=.false.
                   IF(LUMMAT) THEN
                      DO COUNT=FINDCT(GLOBI),FINDCT(GLOBI+1)-1
                         IF(COLCT(COUNT).EQ.GLOBI) ICENT=COUNT
                      END DO
                   ENDIF
                   vec(globi)=1.0
                   DO SJLOC=1,SNLOC
                      JLOC=SILOC2ILOC(SJLOC)
                      GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)

                      OUTVLMX=0.0
                      OUTVLMY=0.0
                      OUTVLMZ=0.0
                      DO SGI=1,SNGI
                         RNN=SN(SILOC,SGI)*SN(SJLOC,SGI)*SDETWE(SGI)
                         OUTVLMX=OUTVLMX+RNN*SNORMXN(SGI)
                         OUTVLMY=OUTVLMY+RNN*SNORMYN(SGI)
                         OUTVLMZ=OUTVLMZ+RNN*SNORMZN(SGI)
                      END DO
                      ! Find POSMAT. ***************************************
                      IF(SLUM) THEN
                         GLOBJ2=GLOBI
                      ELSE
                         GLOBJ2=GLOBJ
                      ENDIF
                      CALL POSINMAT(POSMAT,GLOBI,GLOBJ2, &
                           NONODS,FINDCT,COLCT,NCT)
                      ! Found POSMAT ***************************************
                      IF(LUMMAT) THEN
                         C1T(ICENT)=C1T(ICENT)+OUTVLMX
                         C2T(ICENT)=C2T(ICENT)+OUTVLMY
                         C3T(ICENT)=C3T(ICENT)+OUTVLMZ
                      ELSE
                         C1T(POSMAT)=C1T(POSMAT)+OUTVLMX
                         C2T(POSMAT)=C2T(POSMAT)+OUTVLMY
                         C3T(POSMAT)=C3T(POSMAT)+OUTVLMZ
                      ENDIF

                   END DO
                END DO
             ENDIF
          ENDIF

       END DO face_loop

    END DO second_element_loop

    ewrite(2,*) 'r2norm(c1t,nct,0):',r2norm(c1t,nct,0)
    ewrite(2,*) 'r2norm(c2t,nct,0):',r2norm(c2t,nct,0)
    ewrite(2,*) 'r2norm(c3t,nct,0):',r2norm(c3t,nct,0)

  END SUBROUTINE CTDGSUP

  SUBROUTINE ELEMENT_ORDER(ELEORD,NLEVEL,TOTELE,SNLOC,NLOC,SMLOC,MLOC,SNGI, &
       XONDGL,NDGLNO, &
       COLELE4,LOCLIST,PGNODS,NONODS,XNONOD,X,Y,Z, &
       SN,SNLX,SNLY, SWEIGH, &
       ISPHERE,SM,SMLX,SMLY, &
       FINDRM,COLM,NCOLM, &
       PGNDGLNO,PNORMX,PNORMY,PNORMZ, &
       one_vdg_ser_it,adjits)
    use FLDebug
    IMPLICIT NONE
    INTEGER NFACE,MXITS
    PARAMETER(NFACE=4,MXITS=50)
    ! MXITS=max no of amendment iterations to find a route through mesh
    ! so we dont need to iterate or iterate minimally. 
    INTEGER LOCLIST(NFACE,3)
    INTEGER NLEVEL,TOTELE,SNLOC,NLOC,SMLOC,MLOC,SNGI
    INTEGER one_vdg_ser_it,adjits
    INTEGER PGNODS,NONODS,XNONOD
    INTEGER ELEORD(TOTELE),COLELE4(5*TOTELE)
    ! Global node number maps.
    INTEGER NDGLNO(TOTELE*NLOC),PGNDGLNO(TOTELE*MLOC)
    INTEGER XONDGL(TOTELE*NLOC)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)

    INTEGER ISPHERE
    REAL PNORMX(PGNODS),PNORMY(PGNODS),PNORMZ(PGNODS)
    REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
    REAL SWEIGH(SNGI)
    REAL SM(SMLOC,SNGI),SMLX(SMLOC,SNGI),SMLY(SMLOC,SNGI)
    INTEGER NCOLM
    INTEGER FINDRM(NONODS+1),COLM(NCOLM)

    ! Local variables...
    REAL SDETWE(SNGI)
    REAL SNORMXN(SNGI),SNORMYN(SNGI),SNORMZN(SNGI)
    REAL SDIRX(SNGI),SDIRY(SNGI),SDIRZ(SNGI)
    REAL NDOTQ(SNGI),INCOME(SNGI)

    INTEGER SILOC2ILOC(SNLOC),OTHILOC(NLOC)
    INTEGER MLOCSILOC2ILOC(SMLOC),MLOCOTHILOC(MLOC)
    INTEGER ELE,ILOC
    INTEGER SGI
    INTEGER SILOC,IFACE,ELE2,II,IGL,ITS
    INTEGER COUNT2
    REAL SAREA,RN,SUFPROD,RSUM
    LOGICAL TOPSUF,CONVERG

    INTEGER NODI,NFRONT,IFR,COUNT,ICOL,ELENO,ICOUNT,INLEV,ELENO2
    LOGICAL FOUND

    REAL NORMX,NORMY,NORMZ

    INTEGER, ALLOCATABLE, DIMENSION(:)::IFRONT
    INTEGER, ALLOCATABLE, DIMENSION(:)::INVELEORD
    REAL, ALLOCATABLE, DIMENSION(:)::SUF_IN_FLX
    INTEGER, ALLOCATABLE, DIMENSION(:)::ISUF_IN_FLX

    ALLOCATE(IFRONT(NONODS))
    ALLOCATE(INVELEORD(TOTELE))
    IFRONT(1:NONODS)=0
    ALLOCATE(SUF_IN_FLX(5*TOTELE))
    ALLOCATE(ISUF_IN_FLX(5*TOTELE))
    SUF_IN_FLX=0.0
    ISUF_IN_FLX=0

    IF(NLEVEL.LE.1) THEN

       ! For an unstructured mesh...
       element_loop: DO ELE=1,TOTELE

          ! Put surface integrals into matrix and rhs*****************************
          face_loop: DO IFACE=1,4
             COUNT2=(ELE-1)*5+IFACE
             ! NB colele4 is ordered in terms of faces.
             ELE2=COLELE4(COUNT2)

             ! Surface eleemnt of domain...
             DO SILOC=1,SNLOC
                SILOC2ILOC(SILOC)=LOCLIST(IFACE,SILOC)
             END DO
             ! Form OTHILOC from SILOC2ILOC
             IF(ELE2.NE.0) THEN
                CALL ELE2ELEFACE(OTHILOC,SILOC2ILOC,SNLOC,NLOC,TOTELE, &
                     XONDGL,NDGLNO,ELE,ELE2)
             ELSE
                DO ILOC=1,NLOC
                   OTHILOC(ILOC)=ILOC
                END DO
             ENDIF

             ! Form MLOCOTHILOC and MLOCSILOC2ILOC
             CALL MLOCELE2ELEFACE(MLOCOTHILOC,MLOCSILOC2ILOC,SNLOC,NLOC, &
                  SMLOC,MLOC,OTHILOC,SILOC2ILOC)

             ! PERFORM SURFACE INTEGRATION         !
             ! Form approximate surface normal (NORMX,NORMY,NORMZ)
             CALL DGSIMPLNORM(ELE,SILOC2ILOC,TOTELE,NLOC,SNLOC,XONDGL, &
                  X,Y,Z,XNONOD,NORMX,NORMY,NORMZ)

             CALL DGSDETNXLOC(ELE,SILOC2ILOC,XONDGL, &
                  TOTELE,NLOC,XNONOD,SNLOC,SNGI, &
                  X,Y,Z, &
                  SN,SNLX,SNLY, SWEIGH, SDETWE,SAREA, .TRUE.,.FALSE., &
                  SNORMXN,SNORMYN,SNORMZN, &
                  NORMX,NORMY,NORMZ, &
                  ISPHERE,SMLOC,SM,SMLX,SMLY)

             TOPSUF=.FALSE.
             SUFPROD=0.0
             DO SGI=1,SNGI
                SDIRX(SGI)=0.0
                SDIRY(SGI)=0.0
                SDIRZ(SGI)=0.0
                DO SILOC=1,SMLOC
                   ILOC=MLOCSILOC2ILOC(SILOC)
                   IGL=PGNDGLNO((ELE-1)*MLOC+ILOC)
                   SDIRX(SGI)=SDIRX(SGI)-SM(SILOC,SGI)*PNORMX(IGL)
                   SDIRY(SGI)=SDIRY(SGI)-SM(SILOC,SGI)*PNORMY(IGL)
                   SDIRZ(SGI)=SDIRZ(SGI)-SM(SILOC,SGI)*PNORMZ(IGL)
                END DO
                RN=SQRT(SDIRX(SGI)**2+SDIRY(SGI)**2+SDIRZ(SGI)**2)
                SDIRX(SGI)=SDIRX(SGI)/RN
                SDIRY(SGI)=SDIRY(SGI)/RN
                SDIRZ(SGI)=SDIRZ(SGI)/RN
                NDOTQ(SGI)=SNORMXN(SGI)*SDIRX(SGI)+SNORMYN(SGI)*SDIRY(SGI) &
                     +SNORMZN(SGI)*SDIRZ(SGI)
                INCOME(SGI)=0.0
                IF(NDOTQ(SGI).LT.0.0) INCOME(SGI)=1.0
                SUFPROD=SUFPROD-INCOME(SGI)*NDOTQ(SGI)*SDETWE(SGI)
                IF(ELE2.EQ.0) THEN
                   IF(NDOTQ(SGI).LT.-0.8) TOPSUF=.TRUE.
                ENDIF
             END DO
             SUF_IN_FLX(COUNT2)=SUF_IN_FLX(COUNT2)+SUFPROD

             ! Perform surface integration...
             IF(TOPSUF) THEN
                DO SILOC=1,SNLOC
                   ILOC=SILOC2ILOC(SILOC)
                   NODI=NDGLNO((ELE-1)*NLOC+ILOC)
                   IFRONT(NODI)=1
                END DO
             ENDIF

          END DO face_loop

       END DO element_loop

       ! Grow the front from the front on the surface...
       NFRONT=1
       DO IFR=1,100000
          FOUND=.FALSE.
          DO NODI=1,NONODS
             IF(IFRONT(NODI).EQ.IFR) THEN
                DO COUNT=FINDRM(NODI),FINDRM(NODI+1)-1
                   ICOL=COLM(COUNT)
                   IF(IFRONT(ICOL).EQ.0) THEN
                      IFRONT(ICOL)=IFR+1
                      FOUND=.TRUE.
                   ENDIF
                END DO
             ENDIF
          END DO
          IF(.NOT.FOUND) exit
          NFRONT=NFRONT+1
       END DO
       ! ENDOF IF(NLEVEL.LE.1) THEN ...
    ENDIF

    IF(NLEVEL.GT.1) THEN
       ! Structured mesh...
       NFRONT=NLEVEL
       INLEV=NONODS/NLEVEL
       do nodi=1,nonods
          ifront(nodi)=(NLEVEL+1)-(INT((NODi-1)/INLEV)+1)
       end do
    ENDIF

    ! From fronts get ordering of elements starting from the top...
    INVELEORD(1:TOTELE)=0
    ELEORD(1:TOTELE)=0
    ELENO=0
    DO IFR=1,NFRONT
       DO II=4,1,-1
          DO ELE=1,TOTELE
             IF(INVELEORD(ELE).EQ.0) THEN
                ! Count no of nodes of element ELE in front...
                ICOUNT=0
                DO ILOC=1,NLOC
                   NODI=NDGLNO((ELE-1)*NLOC+ILOC)
                   IF(IFRONT(NODI).EQ.IFR) ICOUNT=ICOUNT+1
                END DO
                IF(ICOUNT.EQ.II) THEN
                   ELENO=ELENO+1
                   INVELEORD(ELE)=ELENO
                   ELEORD(ELENO)=ELE
                ENDIF
             ENDIF
          END DO
       END DO
    END DO
    ewrite(2,*) 'eleno,totele,NFRONT:',eleno,totele,NFRONT
    
    IF(NLEVEL.LE.1) THEN
    ! Amend ordering so that we try and get a route through
    ! by swopping over ordering if ordering is the opposite. 
      DO ELENO=1,TOTELE
        ELE=ELEORD(ELENO)
        RSUM=0.0
        DO IFACE=1,4
             COUNT2=(ELE-1)*5+IFACE
             ! NB colele4 is ordered in terms of faces.
             RSUM=RSUM+SUF_IN_FLX(COUNT2)
        END DO
        DO IFACE=1,4
             COUNT2=(ELE-1)*5+IFACE
             ! NB colele4 is ordered in terms of faces.
             ELE2=COLELE4(COUNT2)
             IF(ELE2.NE.0) THEN
             IF(RSUM.GT.1.E-20) THEN
                IF(SUF_IN_FLX(COUNT2)/RSUM.GT.0.001) THEN
                   ISUF_IN_FLX(COUNT2)=1
                ENDIF
             ENDIF
             ENDIF
        END DO
      END DO
      
    DO ITS=1,MXITS
      CONVERG=.TRUE.
      DO ELENO=1,TOTELE
        ELE=ELEORD(ELENO)
        DO IFACE=1,4
             COUNT2=(ELE-1)*5+IFACE
             ! NB colele4 is ordered in terms of faces.
             ELE2=COLELE4(COUNT2)
             IF(ELE2.NE.0) THEN
               ELENO2=INVELEORD(ELE2)
               IF(ELENO2.GT.ELENO) THEN
               IF(ISUF_IN_FLX(COUNT2).EQ.1) THEN
! nb: INVELEORD(ELE)=ELENO
! nb: ELEORD(ELENO)=ELE
! Swop over...
                  ELEORD(ELENO)=ELE2
                  ELEORD(ELENO2)=ELE
                  INVELEORD(ELE)=ELENO2
                  INVELEORD(ELE2)=ELENO 
                  CONVERG=.FALSE.
               ENDIF
               ENDIF
             ENDIF
        END DO
      END DO
      IF(CONVERG) THEN
          one_vdg_ser_it=1
          adjits=ITS
          EXIT
      ENDIF
    END DO
    ENDIF

    DEALLOCATE(IFRONT)
    DEALLOCATE(INVELEORD)
    DEALLOCATE(SUF_IN_FLX)
    DEALLOCATE(ISUF_IN_FLX)

  END SUBROUTINE ELEMENT_ORDER

end module verticaldg_module
