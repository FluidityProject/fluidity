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

module assnav_module

  use detnlxr_module
  use allsorts
  use Coordinates
  use Element_Numbering, only: ilink2
  use fields
  use FLdebug
  use global_parameters, only: new_options, OPTION_PATH_LEN
  use shape_transformations
  use solvers
  use sparse_tools
  use state_module
  use Pressure_Gradient_Matrix_Wrapper
  use divergence_matrix_cv
  use solid_assembly, only: solid3d
  use hart3d_allsorts
  use tr2d_module
  use sparse_tools
  use momentum_matrix_cg_2d
  use position_in_matrix
  use spud
  use coriolis_module
  
  use diff3d_module
  
  implicit none
    
  private
  
  public ::  assnav, amendmvmesh, gespsi, inv3x3, mapgravdis, midsnods, &
    & pgassem1p, pgmome1p, vispre
  
contains

  SUBROUTINE ASSNAV&
  &         (acctim,U,V,W,&
  !     --------------------------------------START OF ADD BY CRGW 31/03/06
  !                   PLASTIC STRAINS:
  &          STVPXX,STVPYY,STVPZZ,&
  &          STVPYZ,STVPXZ,STVPXY,&
  !     --------------------------------------END OF ADD BY CRGW 31/03/06
  &          NU,NV,NW,UG,VG,WG,&
  &          SOURCX,SOURCY,SOURCZ,X,Y,Z,R0,D0,&
  &     NBUOY,BOY_ML,SOXGRA,SOYGRA,SOZGRA, &
  &     NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB, &
  &          VECX,VECY,VECZ, &
  &          C1T,C2T, BIGM1,NBIGM, &
  &          NCMC,&
  &          GRAFOR,ML,&
  &          FINDCT,COLCT,NCT,FREDOP,TOTELE, &
  &          FINDRM,COLM,NCOLM,NONODS,CENTRM, &
  &          M,&
  &          N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC, &
  &          DISOPT,DISOPN,DT,THETA,BETA,&
  &          GEOBAL,&
  &          LUMP,&
  &          MAKSYM,&
  &          GETC12,&
  &          CONVIS, &
  &          ABSLUM,SLUMP,&
  &          NDGLNO,PNDGLN, &
  &          SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
  &          DENPT,&
  &          MUPTXX,MUPTXY,MUPTXZ, MUPTYY,MUPTYZ,MUPTZZ,&
  &          ML2MXX,ML2MXY,ML2MXZ, ML2MYY,ML2MYZ,ML2MZZ, &
  &          XABSOR,YABSOR,ZABSOR, CVMVEC,&
  &          D3,DCYL,DSPH,&
  &          STOTEL,SNLOC,SNGI,SNDGLN,  &
  &          SN,SNLX,SNLY,SWEIGH, &
  &          SCFACTH0,&
  &          GOTTRAF, &
  &          VERSIO, ISPHERE,&
  &          nnodp, para, halo_tag, rotat, nodrot,nnodro,&
  &          state, istate, &
  &          big_m, CT_m, CTP_m, CMC_m,&
  !     --------------------------------------START OF ADD BY CRGW 14/03/06
  !               INTEGER SOLIDS FLAG:
  &          SOLIDS, MKCOMP,&
  &          pressure_option_path &
  !     --------------------------------------END OF ADD BY CRGW 14/03/06
  &      )

    ! This subroutine assebles the Navier Stokes equations 
    ! FOR SOLUTION LATER. 
    !
    ! N.B.  When DTLAX=0.5 *DT we have something similar to lax-wendrof scheme. 
    ! The Petrof Galerkin weighting is controlled through DTLAX. 
    ! IF SYM then make BIGM symmetrix . 
    ! IF LUMP then lump the mass matrix else dont. 
    ! BETA controls the conservative discretisation 
    ! BETA=0. is the usual advection form, BETA=1. is divergence form. 
    ! BETA=0.5 is an average of the divergence and advection form.(recommended).
    !  IF LEASQR then solve using a least squares approach it 
    ! only works for non-Navier Stokes equations with MU=0. and DTLAX=DT.
    ! NU,NV,NW are for the non-linear terms ordeneraly NU=U,NV=V,NW=W.
    ! UG,VG,WG are the grid velocities. 
    ! IF(QUICK) then use reduced quadrature and assume diff matrix is in 
    ! If DCYL BUT NOT 2.5 D cylinderical coords then W=0.
    ! IF DSPH then 3-D spherical coords. 
    ! H is for Petrof-Galerkin
    ! NB SOURCE, DENPT, viscosity are SNONOD in length and have ele-pt SONDGL.
    ! CVMVEC=virtual mass vector (is 0 in non 2-phase flow). 
    LOGICAL INMAT
    INTEGER LESNOD
    ! FOURTH=.true. then have a 4th order LES
    ! if LESNOD=1 then have a node-wise material 
    ! for LES - nec. for proper 4th order LES.
    PARAMETER(INMAT=.TRUE.)
    INTEGER SNONOD,VNONOD,XNONOD,NCT
    INTEGER DISOPT,DISOPN,GEOBAL
    INTEGER MLOC,NLOC,NGI,NCOLM,NBIGM
    INTEGER STOTEL,SNLOC,SNGI
    REAL DT,THETA,BETA,acctim
    INTEGER TOTELE,NONODS,FREDOP
    REAL U(VNONOD),V(VNONOD),W(VNONOD)
    REAL UG(VNONOD),VG(VNONOD),WG(VNONOD)
    REAL NU(VNONOD),NV(VNONOD),NW(VNONOD)
    REAL VECX(NONODS),VECY(NONODS),VECZ(NONODS)
    REAL BIGM1(NBIGM)
    INTEGER NCMC
    REAL GRAFOR(NONODS)
    REAL SCFACTH0 
    ! NB NBIGM=NCOLM when VARDIS=.FALSE.
    REAL ML(NONODS)
    REAL C1T(NCT),C2T(NCT)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    REAL SOURCX(SNONOD),SOURCY(SNONOD),SOURCZ(SNONOD)
      INTEGER NBUOY
      REAL BOY_ML(NBUOY*6*NONODS)
! soxgra,soygra,sozgra is the direction of gravity force when geobal.le.-10.or.(HYDROS.NE.0)
      REAL SOXGRA(SNONOD),SOYGRA(SNONOD),SOZGRA(SNONOD)
      INTEGER NSOGRASUB
      REAL SOXGRASUB(TOTELE*NSOGRASUB),SOYGRASUB(TOTELE*NSOGRASUB)
      REAL SOZGRASUB(TOTELE*NSOGRASUB)
    !
    REAL DENPT(SNONOD)
    REAL MUPTXX(SNONOD),MUPTXY(SNONOD),MUPTXZ(SNONOD)
    REAL MUPTYY(SNONOD),MUPTYZ(SNONOD),MUPTZZ(SNONOD) 
    REAL ML2MXX(SNONOD),ML2MXY(SNONOD),ML2MXZ(SNONOD)
    REAL ML2MYY(SNONOD),ML2MYZ(SNONOD),ML2MZZ(SNONOD)
    REAL XABSOR(SNONOD),YABSOR(SNONOD),ZABSOR(SNONOD)
    !
    INTEGER NDGLNO(TOTELE*NLOC),PNDGLN(TOTELE*MLOC)
    INTEGER SONDGL(TOTELE*NLOC),VONDGL(TOTELE*NLOC)
    INTEGER XONDGL(TOTELE*NLOC)
    INTEGER FINDRM(NONODS+1),COLM(NCOLM)
    INTEGER CENTRM(NONODS)
    INTEGER FINDCT(FREDOP+1),COLCT(NCT)
    
    logical, intent(in) :: rotat
    integer, intent(in) :: nnodro
    integer, intent(in), dimension(nnodro) :: nodrot
    integer, intent(in) :: nnodp, halo_tag, para
    
    REAL M(MLOC,NGI)
    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL WEIGHT(NGI)
    LOGICAL GETC12
    LOGICAL LUMP,MAKSYM,CONVIS 
    LOGICAL D3,DCYL,DSPH
    LOGICAL ABSLUM,SLUMP
    INTEGER SNDGLN(STOTEL*SNLOC)
    REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
    REAL SWEIGH(SNGI)
    REAL D0,R0
    INTEGER LES4TH
    ! such an efficient use of memory:
    integer, dimension(:), allocatable:: zero_vector
    ! volume fraction of phase(usually=1.)
    REAL CVMVEC(NONODS)
    INTEGER MXNODS
    !
    INTEGER GOTTRAF
    INTEGER VERSIO,ISPHERE
    !    Encapsulated state information containing fields. -dham
    type(state_type), dimension(:), intent(inout) :: state
    ! the current phase index to state
    integer, intent(in) :: istate
    !    System matrices.
    type(block_csr_matrix), intent(inout) :: big_m
    ! Transpose of the pressure gradient matrix
    ! (currently wrapped around C1T, C1T, C3T
    type(block_csr_matrix), intent(inout) :: CT_m
    ! Transpose of the modified compressible pressure gradient matrix
    type(block_csr_matrix), intent(inout) :: CTP_m
    type(csr_matrix), intent(inout) :: CMC_m
    !    Wrappers for Right Hand Side and Mass lumping
!!$    type(vector_field) :: RHS, Velocity
!!$    type(scalar_field) :: masslump
    
    !     --------------------------------------START OF ADD BY CRGW 14/03/06
    !       INTEGER SOLIDS FLAG:
    INTEGER SOLIDS
    !     --------------------------------------END OF ADD BY CRGW 14/03/06
    !     --------------------------------------START OF ADD BY CRGW 31/03/06
    !       PLASTIC STRAINS:
    REAL STVPXX(:), STVPYY(:), STVPZZ(:)
    REAL STVPYZ(:), STVPXZ(:), STVPXY(:)
    !     --------------------------------------END OF ADD BY CRGW 31/03/06
    INTEGER MKCOMP

    character(len=OPTION_PATH_LEN) pressure_option_path

    !
    ! Get the pointers for the element work space
    ewrite(3,*) 'SUBROUTINE ASSNAV()'
    
    !     
    ewrite(3,*)'nonods,snonod,xnonod,vnonod',&
         &            nonods,snonod,xnonod,vnonod

    MXNODS=MAX(NONODS,SNONOD,VNONOD,XNONOD) 
    allocate(zero_vector(1:MXNODS))
    ewrite_minmax(denpt)
    !
    IF(D3) THEN
        IF(DISOPT.EQ.46) THEN
            ! New 4th order method...
            LESNOD=1
            LES4TH=1
          ELSE
            LESNOD=0
            LES4TH=0
          ENDIF
          
          IF(SOLIDS.GT.0) THEN
            
            CALL SOLID3D( VECX, VECY, VECZ, &
                  &               BIGM1, &
                  &               ML, &
                  &               NBIGM, NCOLM, &
                  &               NLOC, NGI, &
                  &               DT, THETA, &
                  &               LUMP, &
                  ! plastic strains:
                  &               STVPXX, STVPYY, STVPZZ, &
                  &               STVPYZ, STVPXZ, STVPXY, &
                  ! solid switch:
                  &               SOLIDS, &
                  &               state(istate) )
            
          ELSE
          
            IF((ISPHERE.EQ.2).OR.((DISOPT.ge.144).AND.(DISOPT.LE.154))) THEN
                ! New 4th order method...
                LESNOD=1
                LES4TH=1
                !
                CALL DIFF3DSIMP&
                    &         (acctim,U,V,W,NU,NV,NW,UG,VG,WG,&
                    &          SOURCX,SOURCY,SOURCZ,X,Y,Z, D0,&
                    &          NBUOY,BOY_ML,SOXGRA,SOYGRA,SOZGRA,&
                    &          VECX,VECY,VECZ, &
                    &          BIGM1,NBIGM,&
                    &          ML,&
                    &          TOTELE, &
                    &          FINDRM,COLM,NCOLM,NONODS,CENTRM, &
                    &          N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI, &
                    &          DISOPT,DISOPN,DT,THETA,BETA,&
                    &          GEOBAL,&
                    &          LUMP,CONVIS,&
                    &          ABSLUM,&
                    &          NDGLNO, &
                    &          SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
                    &          DENPT,&
                    &          MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ, &
                    &          ML2MXX,ML2MXY,ML2MXZ,ML2MYY,ML2MYZ,ML2MZZ, &
                    &          XABSOR,YABSOR,ZABSOR, &
                    &          CVMVEC,&
                    &          LES4TH,&
                    &          LESNOD,&
                    &          SCFACTH0,&
                    &          nnodp, para, halo_tag,&
                    &          ISPHERE)
                ewrite(3,*) 'just outside diff3d'
            ELSE

                CALL DIFF3D&
                &         (acctim,U,V,W,NU,NV,NW,UG,VG,WG,&
                &          SOURCX,SOURCY,SOURCZ,X,Y,Z, D0,&
                &          VECX,VECY,VECZ, &
                &          BIGM1,NBIGM, &
                &          ML,&
                &          TOTELE, &
                &          FINDRM,COLM,NCOLM,NONODS,CENTRM, &
                &          N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI, &
                &          DISOPT,DISOPN,DT,THETA,BETA,&
                &          GEOBAL,&
                &          LUMP,MAKSYM,CONVIS,&
                &          ABSLUM,SLUMP,&
                &          NDGLNO, &
                &          SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
                &          DENPT,&
                &          MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ, &
                &          ML2MXX,ML2MXY,ML2MXZ,ML2MYY,ML2MYZ,ML2MZZ, &
                &          XABSOR,YABSOR,ZABSOR, &
                &          CVMVEC,&
                &          LES4TH,&
                &          LESNOD,&
                &          SCFACTH0,&
                &          nnodp,para,halo_tag,&
                &          GOTTRAF,VERSIO,ISPHERE)
                ewrite(3,*) 'just outside diff3d'
            ENDIF
          ENDIF

          ! diff3d, diff3dsimp and solid3d no longer assembly CT_m as it was being
          ! done twice for mixed cg-cv velocity-pressure discretisations
          ! Now assemble_cts_m assembles both CT_m and CTP_m for such mixed formulations
          ! (entirely in new code) and just CT_m for p1p1 cg formulations.
          call assemble_pressure_gradient_matrix(getc12, mkcomp, CT_m, CTP_m, &
                                                state, istate, pressure_option_path)

      ENDIF
      !
      IF((.NOT.D3).AND.(.NOT.DCYL)) THEN
          !
          CALL DIFF2D(U,V,NU,NV,UG,VG,&
              &          SOURCX,SOURCY,X,Y,&
              &          VECX,VECY,&
              &          C1T,C2T, BIGM1,NBIGM, ML,&
              &          FINDCT,COLCT,NCT,FREDOP, TOTELE, &
              &          FINDRM,COLM,NCOLM,NONODS,CENTRM, &
              &          M,N,NLX,NLY, WEIGHT, NLOC,NGI,MLOC, &
              &          DISOPT,DISOPN,DT,THETA,BETA,LUMP,MAKSYM,CONVIS,&
              &          GETC12,&
              &          ABSLUM,SLUMP,&
              &          NDGLNO,PNDGLN, &
              &          SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
              &          DENPT,&
              &          MUPTXX,MUPTXY,MUPTYY, XABSOR,YABSOR, &
              &          CVMVEC)

    ENDIF
      
    ewrite(1,*)  "END SUBROUTINE ASSNAV()"

  END SUBROUTINE ASSNAV

  SUBROUTINE VISPRE(PSIPRE,PG,NONODS,PGNODS,&
       &               NLOC,MLOC,TOTELE,&
       &               NDGLNO,PGNDGLNO )
    ! This sub maps the quadratic pressure to linear pressure...
    INTEGER NONODS,PGNODS,NLOC,MLOC,TOTELE
    REAL PSIPRE(NONODS),PG(PGNODS)
    INTEGER NDGLNO(TOTELE*NLOC),PGNDGLNO(TOTELE*MLOC)
    ! Local variables...
    INTEGER PNOD,NOD,ILOC,KLOC,ELE
    ! Function..
    IF((MLOC.LE.5).OR.(MLOC.EQ.8)) THEN
       PSIPRE(1:NONODS) = PG(1:NONODS)
    ELSE
       do ELE=1,TOTELE! Was loop 
          do ILOC=1,NLOC! Was loop 
             KLOC=ILINK2(ILOC,ILOC)
             NOD =NDGLNO((ELE-1)*NLOC+ILOC)
             PNOD=PGNDGLNO((ELE-1)*MLOC+KLOC)
             PSIPRE(NOD)=PG(PNOD)
          END DO
       END DO
    ENDIF
    RETURN
  END SUBROUTINE VISPRE
  !
  !
  !
  !
  SUBROUTINE GESPSI(PSIPRE,PG,NONODS,PGNODS,&
       &               NLOC,MLOC,TOTELE,&
       &               NDGLNO,PGNDGLNO )
    ! THIS SUB MAPS LINEAR PRESSURE TO QUADRATIC FOR INITIAL GUESS...
    INTEGER NONODS,PGNODS,NLOC,MLOC,TOTELE
    REAL PSIPRE(NONODS),PG(PGNODS)
    INTEGER NDGLNO(TOTELE*NLOC),PGNDGLNO(TOTELE*MLOC)
    ! Local variables...
    INTEGER PNOD,INOD,JNOD,ILOC,JLOC,KLOC,ELE
    ! Function..
    ewrite(3,*) 'NONODS,PGNODS,NLOC,MLOC,TOTELE:',&
         &            NONODS,PGNODS,NLOC,MLOC,TOTELE
    IF((MLOC.LE.5).OR.(MLOC.EQ.8)) THEN
       PG(1:NONODS) = PSIPRE(1:NONODS)
       IF(MLOC.EQ.5) PG(NONODS+1:NONODS+1 + TOTELE - 1) = 0.0
    ELSE
       do ELE=1,TOTELE! Was loop 
          do ILOC=1,NLOC! Was loop 
             do JLOC=1,NLOC! Was loop 
                KLOC=ILINK2(ILOC,JLOC)
                INOD =NDGLNO((ELE-1)*NLOC+ILOC)
                JNOD =NDGLNO((ELE-1)*NLOC+JLOC)
                PNOD=PGNDGLNO((ELE-1)*MLOC+KLOC)
                PG(PNOD)=0.5*(PSIPRE(INOD)+PSIPRE(JNOD))
             END DO
          END DO
       END DO
       IF(MLOC.EQ.11) PG(PGNODS-TOTELE+1:PGNODS-TOTELE+1 + TOTELE - 1) = 0.0
    ENDIF
    RETURN
  END SUBROUTINE GESPSI

  SUBROUTINE PGASSEM1P(DIFMATH,&
       &          DIFVECX,DIFVECY,DIFVECZ, &
       &          FREES,EQETID,GETFREES,GRAVTY, &
       &          TAKVERT, FSTAKVERT,HBALANC,PNORMX,PNORMY,PNORMZ, &
       &          PGRAVX,PGRAVY,PGRAVZ,NDENST,&
       &          POTGRA,USEPOTGRA,  &
       &          PGFINDRM,PGCOLM,NPGCOLM,PGNODS,&
       &          X,Y,Z, U,V,W,&
       &          M,MLX,MLY,MLZ,&
       &          N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC, &
       &          TOTELE,D3,DCYL, &
       &          GEOBAL,&
       &          NDGLNO,PGNDGLNO, XONDGL,&
       &          NONODS,XNONOD,&
       &          scfacth0,ISPHERE, &
       &          IDGBAL,ELLBOYONLY,ELLCORONLY )
    ! Put the result into momentum eqns VECX=VEXC+N_i dP/dx etc
    INTEGER NONODS,PGNODS,XNONOD,USEPOTGRA,ISPHERE
    INTEGER TOTELE,NLOC,MLOC,NGI
    INTEGER NPGCOLM
    INTEGER PGFINDRM(PGNODS+1),PGCOLM(NPGCOLM)
    INTEGER IDGBAL,ELLBOYONLY,ELLCORONLY
    REAL FREES(PGNODS), EQETID(PGNODS)
    LOGICAL GETFREES
    REAL GRAVTY
    ! PSIx,psiy are the derivatives of hydrostatic pressure in x,y-directions
    REAL DIFMATH(NPGCOLM)
    REAL DIFVECX(PGNODS),DIFVECY(PGNODS),DIFVECZ(PGNODS)
    LOGICAL TAKVERT,FSTAKVERT,HBALANC
    ! If TAKVERT then take out the vertical component of pressure/forcing. 
    ! If FSTAKVERT then take out the vertical component of free surface forcing. 
    ! If HBALANC then consider only horizontal derivatives so exclude buoyancy. 
    REAL PNORMX(PGNODS),PNORMY(PGNODS),PNORMZ(PGNODS)
    REAL PGRAVX(PGNODS),PGRAVY(PGNODS),PGRAVZ(PGNODS), NDENST(NONODS)
    REAL POTGRA(PGNODS*USEPOTGRA)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    REAL U(NONODS),V(NONODS),W(NONODS)
    INTEGER NDGLNO(TOTELE*NLOC),PGNDGLNO(TOTELE*MLOC)
    INTEGER XONDGL(TOTELE*NLOC)

    REAL M(MLOC,NGI),MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL WEIGHT(NGI)
    LOGICAL D3,DCYL
    INTEGER GEOBAL
    real scfacth0
    ! Local variables...
    REAL DETWEI(NGI)
    REAL XD(NGI),YD(NGI),ZD(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    REAL MX(MLOC,NGI),MY(MLOC,NGI),MZ(MLOC,NGI)
    INTEGER ELE,ILOC,JLOC,GLOBI,GLOBJ,POSMAT,GI,L,IGL,IGLP
    REAL NORMGIX(NGI),NORMGIY(NGI),NORMGIZ(NGI)
    REAL SOUXGI(NGI),SOUYGI(NGI),SOUZGI(NGI)
    REAL UD(NGI),VD(NGI),WD(NGI)
    REAL PSIXGI(NGI),PSIYGI(NGI),PSIZGI(NGI)
    REAL PVVX(NGI),PVVY(NGI),PVVZ(NGI)
    REAL FSXGI(NGI),FSYGI(NGI),FSZGI(NGI) 
    REAL GIPGRAVX(NGI),GIPGRAVY(NGI),GIPGRAVZ(NGI)
    REAL A11(NGI),A12(NGI),A13(NGI)
    REAL A21(NGI),A22(NGI),A23(NGI)
    REAL A31(NGI),A32(NGI),A33(NGI) 
    REAL RRX(NGI),RRY(NGI),RRZ(NGI)
    REAL KMATH,KMATX,KMATY,KMATZ
    REAL RX,RY,RZ,FRX,FRY,FRZ,OME2
    REAL DIFX,DIFY,DIFZ

    ! THIS SUB ASSEMBLES THE ELIPTIC EQN FOR GEOSTROPHIC PRESSURE

    DIFMATH(1:NPGCOLM) = 0.0
    DIFVECX(1:PGNODS) = 0.0
    DIFVECY(1:PGNODS) = 0.0
    DIFVECZ(1:PGNODS) = 0.0
    !
    do  ELE=1,TOTELE! Was loop 340

       CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MLOC,NGI, &
            &               N,NLX,NLY,NLZ, M,MLX,MLY,MLZ, WEIGHT, DETWEI,D3,DCYL, &
            &               NX,NY,NZ, MX,MY,MZ,&
            &               A11,A12,A13, A21,A22,A23, A31,A32,A33,&
            &               XD,YD,ZD,&
            &               ISPHERE) 
       !         ewrite(3,*) '1ele=',ele
       !
       do  GI=1,NGI! Was loop 331
          GIPGRAVX(GI)=0.
          GIPGRAVY(GI)=0.
          GIPGRAVZ(GI)=0.
          NORMGIX(GI)=0.
          NORMGIY(GI)=0.
          NORMGIZ(GI)=0.
          PVVX(GI)=0.
          PVVY(GI)=0.
          PVVZ(GI)=0.
          FSXGI(GI)=0.
          FSYGI(GI)=0.
          FSZGI(GI)=0.
          do  L=1,MLOC! Was loop 379
             IGL =PGNDGLNO((ELE-1)*MLOC+L)
             IF(USEPOTGRA.EQ.1) THEN
                GIPGRAVX(GI)=GIPGRAVX(GI) + MX(L,GI)*POTGRA(IGL)
                GIPGRAVY(GI)=GIPGRAVY(GI) + MY(L,GI)*POTGRA(IGL)
                GIPGRAVZ(GI)=GIPGRAVZ(GI) + MZ(L,GI)*POTGRA(IGL)
             ELSE
                GIPGRAVX(GI)=GIPGRAVX(GI) + M(L,GI)*PGRAVX(IGL)
                GIPGRAVY(GI)=GIPGRAVY(GI) + M(L,GI)*PGRAVY(IGL)
                GIPGRAVZ(GI)=GIPGRAVZ(GI) + M(L,GI)*PGRAVZ(IGL)
             ENDIF
             NORMGIX(GI)=NORMGIX(GI) + M(L,GI)*PNORMX(IGL)
             NORMGIY(GI)=NORMGIY(GI) + M(L,GI)*PNORMY(IGL)
             NORMGIZ(GI)=NORMGIZ(GI) + M(L,GI)*PNORMZ(IGL)
             IF(GETFREES) THEN 
                FSXGI(GI)=FSXGI(GI) + MX(L,GI)*(FREES(IGL)-EQETID(IGL))*GRAVTY
                FSYGI(GI)=FSYGI(GI) + MY(L,GI)*(FREES(IGL)-EQETID(IGL))*GRAVTY
                FSZGI(GI)=FSZGI(GI) + MZ(L,GI)*(FREES(IGL)-EQETID(IGL))*GRAVTY
             ENDIF
          end do ! Was loop 379
          !           SX=0.0
          !           SY=0.0
          !           SZ=0.0
          SOUXGI(GI)=0.0
          SOUYGI(GI)=0.0
          SOUZGI(GI)=0.0
          UD(GI)=0.0
          VD(GI)=0.0
          WD(GI)=0.0
          PSIXGI(GI)=0.0
          PSIYGI(GI)=0.0
          PSIZGI(GI)=0.0
          do L=1,NLOC! Was loop 
             IGL =NDGLNO((ELE-1)*NLOC+L) 
             IGLP=(ELE-1)*NLOC+L

             IF(ABS(GEOBAL).NE.20) THEN
                IF(.NOT.HBALANC) THEN
                   SOUXGI(GI)=SOUXGI(GI) + N(L,GI)*GIPGRAVX(GI)*NDENST(IGL)
                   SOUYGI(GI)=SOUYGI(GI) + N(L,GI)*GIPGRAVY(GI)*NDENST(IGL)
                   SOUZGI(GI)=SOUZGI(GI) + N(L,GI)*GIPGRAVZ(GI)*NDENST(IGL)
                ENDIF
             ENDIF
             UD(GI)=UD(GI) + N(L,GI)*U(IGL)
             VD(GI)=VD(GI) + N(L,GI)*V(IGL)
             WD(GI)=WD(GI) + N(L,GI)*W(IGL)
          END DO


          ! calculate RRX,RRY,RRZ...
          OME2=2.*FUNOME(XD(GI),YD(GI),ZD(GI))
          IF((GEOBAL.EQ.-12).or.(GEOBAL.EQ.-22)) OME2=0.
          IF(ELLBOYONLY.EQ.1) OME2=0.
          IF(ELLCORONLY.EQ.1) THEN
             SOUXGI(GI)=0.0
             SOUYGI(GI)=0.0
             SOUZGI(GI)=0.0
          ENDIF

          RRX(GI)=0.0
          RRY(GI)=0.0
          RRZ(GI)=0.0
          !
          IF(.NOT.(TAKVERT.OR.(IDGBAL.EQ.1))) THEN
             RRX(GI)=RRX(GI) -PSIXGI(GI)          +VD(GI)*OME2+ SOUXGI(GI) 
             RRY(GI)=RRY(GI) -PSIYGI(GI)          -UD(GI)*OME2+ SOUYGI(GI) 
             RRZ(GI)=RRZ(GI)                                  + SOUZGI(GI) 
          ENDIF
          !
          IF(TAKVERT.OR.(IDGBAL.EQ.1)) THEN 
             ! Subtract out the horizontal pressure gradient

             ! We have not included SOUXGI(GI) etc because we assume that it is in gravity direction.  
             RX= -PVVX(GI)-PSIXGI(GI)+VD(GI)*OME2 + SOUXGI(GI)
             RY= -PVVY(GI)-PSIYGI(GI)-UD(GI)*OME2 + SOUYGI(GI) 
             RZ= -PVVZ(GI) + SOUZGI(GI)
             RRX(GI)=RRX(GI) +RX
             RRY(GI)=RRY(GI) +RY
             RRZ(GI)=RRZ(GI) +RZ
             ! This is for subtracting out the vertical component...
             ! The rotation matrix in 3-D is R=  
             !          T1X   T1Y   T1Z
             !          T2X   T2Y   T2Z
             !          NORMX NORMY NORMZ
             ! (DIFX,DIFY,DIFZ)=R (RX,RY,RZ)
             DIFX=0.0
             DIFY=0.0
             DIFZ=NORMGIX(GI)*RX+ NORMGIY(GI)*RY+ NORMGIZ(GI)*RZ
             ! Rotate it back by multiplying by RT to get (VX,VY,VZ).
             RRX(GI)=RRX(GI) - NORMGIX(GI)*DIFZ
             RRY(GI)=RRY(GI) - NORMGIY(GI)*DIFZ
             RRZ(GI)=RRZ(GI) - NORMGIZ(GI)*DIFZ
          ENDIF
          !
          IF(.NOT.FSTAKVERT) THEN
             RRX(GI)=RRX(GI) -FSXGI(GI)
             RRY(GI)=RRY(GI) -FSYGI(GI)
             RRZ(GI)=RRZ(GI) -FSZGI(GI)
          ENDIF
          !
          IF(FSTAKVERT) THEN 
             ! Subtract out the horizontal pressure gradient

             ! We have not included SOUXGI(GI) etc because we assume that it is in gravity direction.  
             FRX= -FSXGI(GI)
             FRY= -FSYGI(GI)
             FRZ= -FSZGI(GI)
             RRX(GI)=RRX(GI) +FRX
             RRY(GI)=RRY(GI) +FRY
             RRZ(GI)=RRZ(GI) +FRZ
             ! This is for subtracting out the vertical component...
             ! The rotation matrix in 3-D is R=  
             !          T1X   T1Y   T1Z
             !          T2X   T2Y   T2Z
             !          NORMX NORMY NORMZ
             ! (DIFX,DIFY,DIFZ)=R (FRX,FRY,FRZ)
             DIFX=0.0
             DIFY=0.0
             DIFZ=NORMGIX(GI)*FRX+ NORMGIY(GI)*FRY+ NORMGIZ(GI)*FRZ
             ! Rotate it back by multiplying by RT to get (VX,VY,VZ).
             RRX(GI)=RRX(GI) - NORMGIX(GI)*DIFZ
             RRY(GI)=RRY(GI) - NORMGIY(GI)*DIFZ
             RRZ(GI)=RRZ(GI) - NORMGIZ(GI)*DIFZ
          ENDIF


       end do ! Was loop 331
       !
       ! FORM THE MATRIX EQUATION...
       do  ILOC=1,MLOC! Was loop 350
          !            ewrite(3,*) 'iloc=',iloc
          GLOBI =PGNDGLNO((ELE-1)*MLOC+ILOC)
          !            ewrite(3,*) 'PGNODS,GLOBI=',PGNODS,GLOBI
          ! 
          !
          do  JLOC=1,MLOC! Was loop 360
             !
             GLOBJ =PGNDGLNO((ELE-1)*MLOC+JLOC)

             KMATH=0. 
             !
             do  GI=1,NGI! Was loop 458
                KMATH=KMATH+DETWEI(GI)*(MX(ILOC,GI)*MX(JLOC,GI)&
                     &             +MY(ILOC,GI)*MY(JLOC,GI)+ MZ(ILOC,GI)*MZ(JLOC,GI))
             end do ! Was loop 458

             !            if((globi.eq.1).and.(globj.eq.1)) then
             !     print *,'ele,kmath:',ele,kmath
             !     endif

             !
             ! Find POSMAT. ***************************************
             CALL POSINMAT(POSMAT,GLOBI,GLOBJ,&
                  &                   PGNODS,PGFINDRM,PGCOLM,NPGCOLM)
             ! Found POSMAT ***************************************
             DIFMATH(POSMAT)=DIFMATH(POSMAT)+KMATH
          end do ! Was loop 360
          !
          ! Now obtain r.h.s vector
          KMATX=0.
          KMATY=0.
          KMATZ=0.
          !
          do  GI=1,NGI! Was loop 4582
             KMATX=KMATX+DETWEI(GI)*MX(ILOC,GI)*RRX(GI)
             KMATY=KMATY+DETWEI(GI)*MY(ILOC,GI)*RRY(GI) 
             KMATZ=KMATZ+DETWEI(GI)*MZ(ILOC,GI)*RRZ(GI) 
          end do ! Was loop 4582
          !
          DIFVECX(GLOBI)=DIFVECX(GLOBI)+KMATX
          DIFVECY(GLOBI)=DIFVECY(GLOBI)+KMATY
          DIFVECZ(GLOBI)=DIFVECZ(GLOBI)+KMATZ
          !
       end do ! Was loop 350
    end do ! Was loop 340

    !       CALL PMINMX(DIFVECX,PGNODS,'******DIFVECX  ')
    !       CALL PMINMX(DIFVECY,PGNODS,'******DIFVECY  ')
    !       CALL PMINMX(DIFVECZ,PGNODS,'******DIFVECZ  ')
    !       stop 383

    ewrite(3,*) 'NPGCOLM=',NPGCOLM
    !        ewrite(3,*) 'r2norm(difmat,ncolm,0):',r2norm(difmat,NPGCOLM,0)
    !        ewrite(3,*) 'DIFVEC(i):',(DIFVEC(i),i=1,30)
    !         do globi=1,-pgnods
    !           ewrite(3,*) globi,(difmat(posmat),posmat=pgfindrm(globi),pgfindrm(globi+1)-1)
    !           kmat=0.
    !           do posmat=pgfindrm(globi),pgfindrm(globi+1)-1
    !               kmat=kmat+difmat(posmat)
    !            end do
    !             ewrite(3,*) 'globi,kmat,c:',globi,kmat,difmat(pgcentrm(globi))
    !         end do
    !       ewrite(3,*) 'r2norm(difvec,pgnods,0):',r2norm(difvec,pgnods,0)

    RETURN
  END SUBROUTINE PGASSEM1P

  SUBROUTINE MAPGRAVDIS(GRAVOPT,DIFMATH,P,&
       &          VECGRA1,VECGRA2,VECGRA3,&
       &          PGRAVX,PGRAVY,PGRAVZ,&
       &          PGFINDRM,PGCOLM,NPGCOLM,PGNODS, &
       &          X,Y,Z, &
       &          M,MLX,MLY,MLZ,&
       &          N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC, &
       &          TOTELE,D3,DCYL,  &
       &          PGNDGLNO, XONDGL,&
       &          XNONOD,ISPHERE)
    ! discretise the eqns for correcting gravity so it is consistent with mesh.
    ! GRAVOPT=1 Form diffusion eqn for pressure
    ! GRAVOPT=2 form mass eqns for GRAVITY components
    INTEGER GRAVOPT,ISPHERE
    INTEGER PGNODS,XNONOD
    INTEGER TOTELE,NLOC,MLOC,NGI
    INTEGER NPGCOLM
    INTEGER PGFINDRM(PGNODS+1),PGCOLM(NPGCOLM)
    REAL DIFMATH(NPGCOLM)
    REAL P(PGNODS)
    REAL VECGRA1(PGNODS),VECGRA2(PGNODS),VECGRA3(PGNODS)
    REAL PGRAVX(PGNODS),PGRAVY(PGNODS),PGRAVZ(PGNODS)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    INTEGER PGNDGLNO(TOTELE*MLOC)
    INTEGER XONDGL(TOTELE*NLOC)

    REAL M(MLOC,NGI),MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL WEIGHT(NGI)
    LOGICAL D3,DCYL

    ! Local variables...
    INTEGER ELE,ILOC,JLOC,GLOBI,GLOBJ,POSMAT,GI,L,IGL
    REAL DETWEI(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    REAL MX(MLOC,NGI),MY(MLOC,NGI),MZ(MLOC,NGI)
    REAL GRAVXD(NGI),GRAVYD(NGI),GRAVZD(NGI)
    REAL PX(NGI),PY(NGI),PZ(NGI)
    REAL A11(NGI),A12(NGI),A13(NGI)
    REAL A21(NGI),A22(NGI),A23(NGI)
    REAL A31(NGI),A32(NGI),A33(NGI) 
    REAL XD(NGI),YD(NGI),ZD(NGI)
    REAL KMAT,NN,NNXP,NNYP,NNZP,NNGRAV,RAD

    DIFMATH(1:NPGCOLM) = 0.0
    VECGRA1(1:PGNODS) = 0.0
    VECGRA2(1:PGNODS) = 0.0
    VECGRA3(1:PGNODS) = 0.0

    !
    do  ELE=1,TOTELE! Was loop 340

       CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MLOC,NGI, &
            &               N,NLX,NLY,NLZ, M,MLX,MLY,MLZ, WEIGHT, DETWEI,D3,DCYL, &
            &               NX,NY,NZ, MX,MY,MZ, &
            &               A11,A12,A13, A21,A22,A23, A31,A32,A33,&
            &               XD,YD,ZD,&
            &               ISPHERE) 

       do  GI=1,NGI! Was loop 331

          IF(ISPHERE.EQ.0) THEN
             !           IF(.true.) THEN
             GRAVXD(GI)=0.
             GRAVYD(GI)=0.
             GRAVZD(GI)=0.
             do  L=1,MLOC! Was loop 79
                IGL =PGNDGLNO((ELE-1)*MLOC+L) 
                GRAVXD(GI)=GRAVXD(GI) + M(L,GI)*PGRAVX(IGL)
                GRAVYD(GI)=GRAVYD(GI) + M(L,GI)*PGRAVY(IGL)
                GRAVZD(GI)=GRAVZD(GI) + M(L,GI)*PGRAVZ(IGL)
             end do ! Was loop 79
          ELSE
             RAD=SQRT(XD(GI)**2+YD(GI)**2+ZD(GI)**2)
             GRAVXD(GI)=-XD(GI)/RAD
             GRAVYD(GI)=-YD(GI)/RAD
             GRAVZD(GI)=-ZD(GI)/RAD
          ENDIF

          PX(GI)=0.0
          PY(GI)=0.0
          PZ(GI)=0.0
          IF(GRAVOPT.NE.1) THEN
             ! Form mass eqns for GRAVITY components
             do ILOC=1,MLOC! Was loop 
                IGL=PGNDGLNO((ELE-1)*MLOC+ILOC)
                PX(GI)=PX(GI) + MX(ILOC,GI)*P(IGL)
                PY(GI)=PY(GI) + MY(ILOC,GI)*P(IGL)
                PZ(GI)=PZ(GI) + MZ(ILOC,GI)*P(IGL)
             END DO
          ENDIF

       end do ! Was loop 331

       ! FORM THE MATRIX EQUATION...
       do  ILOC=1,MLOC! Was loop 350
          GLOBI =PGNDGLNO((ELE-1)*MLOC+ILOC)
          do  JLOC=1,MLOC! Was loop 360
             !
             GLOBJ =PGNDGLNO((ELE-1)*MLOC+JLOC)

             KMAT=0.0
             NN=0.0
             !
             do  GI=1,NGI! Was loop 458
                KMAT=KMAT+DETWEI(GI)*(MX(ILOC,GI)*MX(JLOC,GI)&
                     &              + MY(ILOC,GI)*MY(JLOC,GI)+MZ(ILOC,GI)*MZ(JLOC,GI) )

                NN=NN + DETWEI(GI)*M(ILOC,GI)*M(JLOC,GI)
             end do ! Was loop 458

             !
             ! Find POSMAT. ***************************************
             CALL POSINMAT(POSMAT,GLOBI,GLOBJ,&
                  &                   PGNODS,PGFINDRM,PGCOLM,NPGCOLM)
             ! Found POSMAT ***************************************
             IF(GRAVOPT.EQ.1) THEN
                ! For diffusion eqn for pressure
                DIFMATH(POSMAT)=DIFMATH(POSMAT)+KMAT
             ELSE
                ! For mass eqns for GRAVITY components
                DIFMATH(POSMAT)=DIFMATH(POSMAT)+NN
             ENDIF

          end do ! Was loop 360
          !
          ! Now obtain r.h.s vector
          NNGRAV=0.
          NNXP=0.
          NNYP=0.
          NNZP=0.
          do  GI=1,NGI! Was loop 4582
             NNGRAV=NNGRAV+DETWEI(GI)&
                  &       *(MX(ILOC,GI)*GRAVXD(GI)+MY(ILOC,GI)*GRAVYD(GI)+MZ(ILOC,GI)*GRAVZD(GI))

             NNXP=NNXP+DETWEI(GI)*M(ILOC,GI)*PX(GI)
             NNYP=NNYP+DETWEI(GI)*M(ILOC,GI)*PY(GI)
             NNZP=NNZP+DETWEI(GI)*M(ILOC,GI)*PZ(GI)
          end do ! Was loop 4582
          !
          IF(GRAVOPT.EQ.1) THEN
             ! For diffusion eqn for pressure
             VECGRA1(GLOBI)=VECGRA1(GLOBI)+NNGRAV
          ELSE
             ! For mass eqns for GRAVITY components
             VECGRA1(GLOBI)=VECGRA1(GLOBI)+NNXP
             VECGRA2(GLOBI)=VECGRA2(GLOBI)+NNYP
             VECGRA3(GLOBI)=VECGRA3(GLOBI)+NNZP
          ENDIF
       end do ! Was loop 350
    end do ! Was loop 340

    RETURN
  END SUBROUTINE MAPGRAVDIS
  !     
  !     
  !     
  !     
  !     
  SUBROUTINE AMENDMVMESH(OPTION,VECX,VECY,VECZ,&
       &     UOLD,VOLD,WOLD,DT,&
       &     X,Y,Z, XOLD,YOLD,ZOLD, &
       &     N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,&
       &     TOTELE,D3,DCYL,  &
       &     XONDGL,NDGLNO, NONODS,XNONOD)
    !     Amend the r.h.s vector for a matrix eqn with a moving mesh. 
    !     Use a lumped mass matrix approach.
    INTEGER OPTION
    INTEGER NONODS,XNONOD
    INTEGER TOTELE,NLOC,NGI
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    REAL XOLD(XNONOD),YOLD(XNONOD),ZOLD(XNONOD)
    REAL VECX(NONODS),VECY(NONODS),VECZ(NONODS)
    REAL UOLD(NONODS),VOLD(NONODS),WOLD(NONODS)
    REAL DT
    INTEGER XONDGL(TOTELE*NLOC),NDGLNO(TOTELE*NLOC)

    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL WEIGHT(NGI)
    LOGICAL D3,DCYL

    !     Local variables...
    INTEGER ELE,ILOC,GLOBI,GI
    REAL DETWEI(NGI),DETWEIOLD(NGI)
    REAL NNDIF,VOLUME
    !     
    do  ELE=1,TOTELE! Was loop 340
       !     
       !     Calculate DETWEI,RA,NX,NY,NZ for element ELE
       CALL DETNNN(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NGI, &
            &        N,NLX,NLY,NLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL) 

       !     Calculate DETWEI,RA,NX,NY,NZ for element ELE
       CALL DETNNN(ELE, XOLD,YOLD,ZOLD, XONDGL, TOTELE,XNONOD,NLOC,NGI, &
            &        N,NLX,NLY,NLZ, WEIGHT, DETWEIOLD,VOLUME, D3,DCYL) 
       !     
       !     FORM THE MATRIX EQUATION...
       do  ILOC=1,NLOC! Was loop 350
          GLOBI =NDGLNO((ELE-1)*NLOC+ILOC)
          NNDIF=0.0
          do  GI=1,NGI! Was loop 458
             NNDIF=NNDIF+(DETWEIOLD(GI)-DETWEI(GI))*N(ILOC,GI)
          end do ! Was loop 458
          VECX(GLOBI)=VECX(GLOBI)+NNDIF*UOLD(GLOBI)/DT
          IF(OPTION.GE.2) VECY(GLOBI)=VECY(GLOBI)+NNDIF*VOLD(GLOBI)/DT
          IF(D3) THEN
             IF(OPTION.GE.3) VECZ(GLOBI)=VECZ(GLOBI)+NNDIF*WOLD(GLOBI)/DT
          ENDIF
       end do ! Was loop 350
    end do ! Was loop 340
    RETURN
  END SUBROUTINE AMENDMVMESH

  SUBROUTINE PGMOME1P(VECX,VECY,VECZ, PG, &
       &          FREES,EQETID,GETFREES,GRAVTY,&
       &          TAKVERT,FSTAKVERT, PNORMX,PNORMY,PNORMZ, &
       &          PGRAVX,PGRAVY,PGRAVZ,NDENST,&
       &          POTGRA,USEPOTGRA,  &
       &          SOURCX,SOURCY,SOURCZ,&
       &          X,Y,Z, U,V,W,&
       &          M,MLX,MLY,MLZ,&
       &          N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC, &
       &          TOTELE,D3,DCYL, &
       &          GEOBAL,&
       &          NDGLNO,PGNDGLNO, XONDGL,&
       &          NONODS,PGNODS,XNONOD,&
       &          scfacth0,ISPHERE, &
       &          IDGBAL,ELLBOYONLY,ELLCORONLY )
    ! Put the result into momentum eqns VECX=VECX+(CORRIOLIS+N_i dP/dx) etc
    ! If GEOBAL=-11 or -21 then include bouyancy in forcing terms...
    ! If GEOBAL=-12 or -22 then ONLY include bouyancy in forcing terms...
    INTEGER NONODS,PGNODS,XNONOD,USEPOTGRA,ISPHERE
    INTEGER TOTELE,NLOC,MLOC,NGI  
    INTEGER IDGBAL,ELLBOYONLY,ELLCORONLY
    LOGICAL TAKVERT,FSTAKVERT,GETFREES
    ! If TAKVERT then take out the vertical component of pressure/forcing. 
    ! If FSTAKVERT take out the vertical free surface gradient. 
    REAL PNORMX(PGNODS),PNORMY(PGNODS),PNORMZ(PGNODS)
    REAL PGRAVX(PGNODS),PGRAVY(PGNODS),PGRAVZ(PGNODS), NDENST(NONODS)
    REAL SOURCX(NONODS),SOURCY(NONODS),SOURCZ(NONODS)
    REAL VECX(NONODS),VECY(NONODS),VECZ(NONODS)
    REAL PG(PGNODS),FREES(PGNODS),EQETID(PGNODS)
    REAL POTGRA(PGNODS*USEPOTGRA)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    REAL U(NONODS),V(NONODS),W(NONODS)
    INTEGER NDGLNO(TOTELE*NLOC),PGNDGLNO(TOTELE*MLOC)
    INTEGER XONDGL(TOTELE*NLOC)

    REAL M(MLOC,NGI),MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL WEIGHT(NGI)
    LOGICAL D3,DCYL
    INTEGER GEOBAL
    REAL GRAVTY
    real scfacth0
    !
    ! Local variables...
    REAL DETWEI(NGI)
    REAL A11(NGI),A12(NGI),A13(NGI)
    REAL A21(NGI),A22(NGI),A23(NGI)
    REAL A31(NGI),A32(NGI),A33(NGI) 
    REAL XD(NGI),YD(NGI),ZD(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    REAL MX(MLOC,NGI),MY(MLOC,NGI),MZ(MLOC,NGI)
    INTEGER ELE,ILOC,GLOBI,GI,L,IGLX,IGL,IGLP
    REAL KMATX,KMATY,KMATZ,DIFX,DIFY,DIFZ
    REAL RX,RY,RZ,FRX,FRY,FRZ
    REAL SX,SY,SZ
    REAL NORMGIX(NGI),NORMGIY(NGI),NORMGIZ(NGI)
    REAL GIPGRAVX(NGI),GIPGRAVY(NGI),GIPGRAVZ(NGI)
    REAL PX(NGI),PY(NGI),PZ(NGI) 
    REAL PVVX(NGI),PVVY(NGI),PVVZ(NGI)
    REAL UD(NGI),VD(NGI),WD(NGI),OME2
    REAL SOUXGI(NGI),SOUYGI(NGI),SOUZGI(NGI) 
    REAL PSIXGI(NGI),PSIYGI(NGI),PSIZGI(NGI) 
    REAL FSXGI(NGI),FSYGI(NGI),FSZGI(NGI) 
    REAL RRX(NGI),RRY(NGI),RRZ(NGI)


    do  ELE=1,TOTELE! Was loop 340
       !
       CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MLOC,NGI, &
            &               N,NLX,NLY,NLZ, M,MLX,MLY,MLZ, WEIGHT, DETWEI,D3,DCYL, &
            &               NX,NY,NZ, MX,MY,MZ,&
            &               A11,A12,A13, A21,A22,A23, A31,A32,A33,&
            &               XD,YD,ZD,&
            &               ISPHERE) 
       !
       do  GI=1,NGI! Was loop 331
          GIPGRAVX(GI)=0.
          GIPGRAVY(GI)=0.
          GIPGRAVZ(GI)=0.
          do L=1,MLOC! Was loop 
             IGL =PGNDGLNO((ELE-1)*MLOC+L)
             IF(USEPOTGRA.EQ.1) THEN
                GIPGRAVX(GI)=GIPGRAVX(GI) + MX(L,GI)*POTGRA(IGL)
                GIPGRAVY(GI)=GIPGRAVY(GI) + MY(L,GI)*POTGRA(IGL)
                GIPGRAVZ(GI)=GIPGRAVZ(GI) + MZ(L,GI)*POTGRA(IGL)
             ELSE
                GIPGRAVX(GI)=GIPGRAVX(GI) + M(L,GI)*PGRAVX(IGL)
                GIPGRAVY(GI)=GIPGRAVY(GI) + M(L,GI)*PGRAVY(IGL)
                GIPGRAVZ(GI)=GIPGRAVZ(GI) + M(L,GI)*PGRAVZ(IGL)
             ENDIF
          END DO
          SX=0.0
          SY=0.0
          SZ=0.0
          UD(GI)=0.
          VD(GI)=0.
          WD(GI)=0.
          SOUXGI(GI)=0.
          SOUYGI(GI)=0.
          SOUZGI(GI)=0. 
          PSIXGI(GI)=0.
          PSIYGI(GI)=0.
          PSIZGI(GI)=0.
          do  L=1,NLOC! Was loop 79
             IGLX=XONDGL((ELE-1)*NLOC+L)
             IGL =NDGLNO((ELE-1)*NLOC+L)
             IGLP=(ELE-1)*NLOC+L
             UD(GI)=UD(GI) + N(L,GI)*U(IGL)
             VD(GI)=VD(GI) + N(L,GI)*V(IGL)
             WD(GI)=WD(GI) + N(L,GI)*W(IGL)
             IF((GEOBAL.EQ.-11).or.(GEOBAL.EQ.-21)&
                  &            .OR.(GEOBAL.EQ.-12).or.(GEOBAL.EQ.-22)) THEN
                !                 SOUXGI(GI)=SOUXGI(GI) + N(L,GI)*SOURCX(IGL)
                !                 SOUYGI(GI)=SOUYGI(GI) + N(L,GI)*SOURCY(IGL)
                !                 SOUZGI(GI)=SOUZGI(GI) + N(L,GI)*SOURCZ(IGL)
                SX=SX + N(L,GI)*SOURCX(IGL)
                SY=SY + N(L,GI)*SOURCY(IGL)
                SZ=SZ + N(L,GI)*SOURCZ(IGL)
                SOUXGI(GI)=SOUXGI(GI) + N(L,GI)*GIPGRAVX(GI)*NDENST(IGL)
                SOUYGI(GI)=SOUYGI(GI) + N(L,GI)*GIPGRAVY(GI)*NDENST(IGL)
                SOUZGI(GI)=SOUZGI(GI) + N(L,GI)*GIPGRAVZ(GI)*NDENST(IGL)
             ENDIF

          end do ! Was loop 79

          NORMGIX(GI)=0.
          NORMGIY(GI)=0.
          NORMGIZ(GI)=0.
          PX(GI)=0.
          PY(GI)=0.
          PZ(GI)=0. 
          PVVX(GI)=0.
          PVVY(GI)=0.
          PVVZ(GI)=0.
          FSXGI(GI)=0.
          FSYGI(GI)=0.
          FSZGI(GI)=0.
          do  L=1,MLOC! Was loop 792
             IGL=PGNDGLNO((ELE-1)*MLOC+L) 
             NORMGIX(GI)=NORMGIX(GI) + M(L,GI)*PNORMX(IGL)
             NORMGIY(GI)=NORMGIY(GI) + M(L,GI)*PNORMY(IGL)
             NORMGIZ(GI)=NORMGIZ(GI) + M(L,GI)*PNORMZ(IGL)
             PX(GI)=PX(GI) + MX(L,GI)*PG(IGL)
             PY(GI)=PY(GI) + MY(L,GI)*PG(IGL)
             PZ(GI)=PZ(GI) + MZ(L,GI)*PG(IGL) 
             IF(GETFREES) THEN 
                FSXGI(GI)=FSXGI(GI) + MX(L,GI)*(FREES(IGL)-EQETID(IGL))*GRAVTY
                FSYGI(GI)=FSYGI(GI) + MY(L,GI)*(FREES(IGL)-EQETID(IGL))*GRAVTY
                FSZGI(GI)=FSZGI(GI) + MZ(L,GI)*(FREES(IGL)-EQETID(IGL))*GRAVTY
             ENDIF
          end do ! Was loop 792

          ! Calculate RRX,RRY,RRZ...
          OME2=2.*FUNOME(XD(GI),YD(GI),ZD(GI))
          IF((GEOBAL.EQ.-12).or.(GEOBAL.EQ.-22)) OME2=0.
          IF(ELLBOYONLY.EQ.1) OME2=0.
          IF(ELLCORONLY.EQ.1) THEN
             SOUXGI(GI)=0.0
             SOUYGI(GI)=0.0
             SOUZGI(GI)=0.0
          ENDIF


          RRX(GI)=0.0
          RRY(GI)=0.0
          RRZ(GI)=0.0
          !
          IF(IDGBAL.NE.0) THEN 
             RX= -PX(GI)
             RY= -PY(GI)
             RZ= -PZ(GI) 
             RRX(GI)=RRX(GI)+RX
             RRY(GI)=RRY(GI)+RY
             RRZ(GI)=RRZ(GI)+RZ
             if(.false.) then
                ! This is for subtracting out the vertical component
                ! The rotation matrix in 3-D is R=  
                !          T1X   T1Y   T1Z
                !          T2X   T2Y   T2Z
                !          NORMX NORMY NORMZ
                ! (DIFX,DIFY,DIFZ)=R (RX,RY,RZ)
                DIFX=0.0
                DIFY=0.0
                DIFZ=NORMGIX(GI)*RX+ NORMGIY(GI)*RY+ NORMGIZ(GI)*RZ
                ! Rotate it back by multiplying by RT. 
                RRX(GI)=RRX(GI)-NORMGIX(GI)*DIFZ
                RRY(GI)=RRY(GI)-NORMGIY(GI)*DIFZ
                RRZ(GI)=RRZ(GI)-NORMGIZ(GI)*DIFZ
             end if
          ENDIF
          !
          IF(IDGBAL.EQ.0) THEN
             !
             IF(.NOT.TAKVERT) THEN
                RRX(GI)=RRX(GI)-PX(GI)-PSIXGI(GI)+VD(GI)*OME2+ SOUXGI(GI)
                RRY(GI)=RRY(GI)-PY(GI)-PSIYGI(GI)-UD(GI)*OME2+ SOUYGI(GI)
                RRZ(GI)=RRZ(GI)-PZ(GI) + SOUZGI(GI)
             ELSE
                RRX(GI)=RRX(GI)-PX(GI)-PSIXGI(GI) 
                RRY(GI)=RRY(GI)-PY(GI)-PSIYGI(GI)
                RRZ(GI)=RRZ(GI)-PZ(GI)
             ENDIF
             !

             IF(TAKVERT) THEN 
                !         IF(.true.) THEN 
                ! Subtract out the horizontal pressure gradient

                ! We may not included SOUXGI(GI) etc because we assume that it is in gravity direction.  
                RX= -PX(GI)+VD(GI)*OME2 + SOUXGI(GI)
                RY= -PY(GI)-UD(GI)*OME2 + SOUYGI(GI)
                RZ= -PZ(GI) + SOUZGI(GI)
                RRX(GI)=RRX(GI)+RX
                RRY(GI)=RRY(GI)+RY
                RRZ(GI)=RRZ(GI)+RZ
                ! This is for subtracting out the vertical component
                ! The rotation matrix in 3-D is R=  
                !          T1X   T1Y   T1Z
                !          T2X   T2Y   T2Z
                !          NORMX NORMY NORMZ
                ! (DIFX,DIFY,DIFZ)=R (RX,RY,RZ)
                DIFX=0.0
                DIFY=0.0
                DIFZ=NORMGIX(GI)*RX+ NORMGIY(GI)*RY+ NORMGIZ(GI)*RZ
                ! Rotate it back by multiplying by RT. 
                RRX(GI)=RRX(GI)-NORMGIX(GI)*DIFZ
                RRY(GI)=RRY(GI)-NORMGIY(GI)*DIFZ
                RRZ(GI)=RRZ(GI)-NORMGIZ(GI)*DIFZ
                ! endof IF(TAKVERT) THEN
             ENDIF



             IF(.NOT.FSTAKVERT) THEN
                RRX(GI)=RRX(GI)-FSXGI(GI)
                RRY(GI)=RRY(GI)-FSYGI(GI)
                RRZ(GI)=RRZ(GI)-FSZGI(GI)
             ENDIF

             IF(FSTAKVERT) THEN 
                ! Subtract out the horizontal pressure gradient

                ! We may not included SOUXGI(GI) etc because we assume that it is in gravity direction.  
                FRX= -FSXGI(GI)
                FRY= -FSYGI(GI)
                FRZ= -FSZGI(GI)
                RRX(GI)=RRX(GI)+FRX
                RRY(GI)=RRY(GI)+FRY
                RRZ(GI)=RRZ(GI)+FRZ
                ! This is for subtracting out the vertical component 
                ! The rotation matrix in 3-D is R=  
                !          T1X   T1Y   T1Z
                !          T2X   T2Y   T2Z
                !          NORMX NORMY NORMZ
                ! (DIFX,DIFY,DIFZ)=R (FRX,FRY,FRZ)
                DIFX=0.0
                DIFY=0.0
                DIFZ=NORMGIX(GI)*FRX+ NORMGIY(GI)*FRY+ NORMGIZ(GI)*FRZ
                ! Rotate it back by multiplying by RT. 
                RRX(GI)=RRX(GI)-NORMGIX(GI)*DIFZ
                RRY(GI)=RRY(GI)-NORMGIY(GI)*DIFZ
                RRZ(GI)=RRZ(GI)-NORMGIZ(GI)*DIFZ
                ! endof IF(TAKVERT) THEN
             ENDIF

             ! ENDOF IF(IDGBAL.NE.0) THEN
          ENDIF


       end do ! Was loop 331

       !
       ! FORM THE MATRIX EQUATION...
       do  ILOC=1,NLOC! Was loop 350
          GLOBI =NDGLNO((ELE-1)*NLOC+ILOC)
          ! 
          KMATX=0.
          KMATY=0.
          KMATZ=0. 
          !
          do  GI=1,NGI! Was loop 458
             KMATX=KMATX+DETWEI(GI)*N(ILOC,GI)*RRX(GI)
             KMATY=KMATY+DETWEI(GI)*N(ILOC,GI)*RRY(GI)
             KMATZ=KMATZ+DETWEI(GI)*N(ILOC,GI)*RRZ(GI)
          end do ! Was loop 458
          !
          VECX(GLOBI)=VECX(GLOBI)+KMATX
          VECY(GLOBI)=VECY(GLOBI)+KMATY
          VECZ(GLOBI)=VECZ(GLOBI)+KMATZ

          !
       end do ! Was loop 350
    end do ! Was loop 340
    !          CALL PMINMX(PG,PGNODS,'******PG  ')
    RETURN
  END SUBROUTINE PGMOME1P
  !
  !
  !
  !
  SUBROUTINE INV3X3(&
       &         AGI,BGI,CGI, DGI,EGI,FGI, GGI,HGI,KGI,&
       &         A11,A12,A13, A21,A22,A23, A31,A32,A33)
    ! This sub finds the inverse of a 3x3 matrix
    REAL AGI,BGI,CGI, DGI,EGI,FGI, GGI,HGI,KGI
    REAL A11,A12,A13, A21,A22,A23, A31,A32,A33
    ! Local variables...
    REAL DETJ
    !
    DETJ=AGI*(EGI*KGI-FGI*HGI)&
         &             -BGI*(DGI*KGI-FGI*GGI)&
         &             +CGI*(DGI*HGI-EGI*GGI)
    !
    A11= (EGI*KGI-FGI*HGI) /DETJ
    A21=-(DGI*KGI-FGI*GGI) /DETJ
    A31= (DGI*HGI-EGI*GGI) /DETJ
    !
    A12=-(BGI*KGI-CGI*HGI) /DETJ
    A22= (AGI*KGI-CGI*GGI) /DETJ
    A32=-(AGI*HGI-BGI*GGI) /DETJ
    !
    A13= (BGI*FGI-CGI*EGI) /DETJ
    A23=-(AGI*FGI-CGI*DGI) /DETJ
    A33= (AGI*EGI-BGI*DGI) /DETJ
    RETURN
  END SUBROUTINE INV3X3

  SUBROUTINE MIDSNODS(FINDRM,NONODS,FREDOP,&
       &                CENTRM,COLM,NCOLM,&
       &                NDGLNO,PNDGLN,TOTELE,NLOC,MLOC)
    INTEGER NONODS,NCOLM,NLOC,MLOC,TOTELE,FREDOP
    INTEGER FINDRM(NONODS+1),CENTRM(NONODS)
    INTEGER COLM(NCOLM)
    INTEGER NDGLNO(TOTELE*NLOC),PNDGLN(TOTELE*MLOC)
    ! FUNCTION
    INTEGER COUNT,I,ELE,PNODS,ILOC,JLOC,JCOL,INOD,JNOD,KLOC
    ! WORK SPACE...
    INTEGER, ALLOCATABLE, DIMENSION(:)::DONMAT
    ALLOCATE(DONMAT(NCOLM))
    !
    do I=1,NCOLM! Was loop 
       DONMAT(I)=0
    END DO
    do ELE=1,TOTELE! Was loop 
       do KLOC=1,MLOC! Was loop 
          PNDGLN((ELE-1)*MLOC+KLOC)=-1
       END DO
    END DO
    !
    PNODS=0
    do ELE=1,TOTELE! Was loop 
       do ILOC=1,NLOC! Was loop 
          do JLOC=1,NLOC! Was loop 
             INOD=NDGLNO((ELE-1)*NLOC+ILOC)
             JNOD=NDGLNO((ELE-1)*NLOC+JLOC)
             KLOC=ILINK2(ILOC,JLOC)
             ! Find JNOD...(cONSIDER ONLY UPPER DIAGONAL PART OF DONMAT)
             do COUNT=CENTRM(INOD),FINDRM(INOD+1)-1! Was loop 
                JCOL=COLM(COUNT)
                IF(JCOL.EQ.JNOD) THEN 
                   IF(DONMAT(COUNT).EQ.0) THEN
                      PNODS=PNODS+1
                      PNDGLN((ELE-1)*MLOC+KLOC)=PNODS
                      DONMAT(COUNT)=PNODS
                   ELSE
                      PNDGLN((ELE-1)*MLOC+KLOC)=DONMAT(COUNT)
                   ENDIF
                ENDIF
             ENDDO
             !
          END DO
       END DO
    END DO
    !
    IF(MLOC.EQ.11) THEN
       ! Add a bubble function...
       do ELE=1,TOTELE! Was loop 
          PNDGLN((ELE-1)*MLOC+MLOC)=PNODS + ELE
       END DO
       PNODS=PNODS+TOTELE
    ENDIF
    !
    FREDOP=PNODS
    !
    !       Sanity checks
    do ELE=1,TOTELE! Was loop 
       do KLOC=1,MLOC! Was loop 
          IF((PNDGLN((ELE-1)*MLOC+KLOC).LE.0).OR.&
               &          (PNDGLN((ELE-1)*MLOC+KLOC).GT.FREDOP)) THEN
             ewrite(3,*) 'ELE,KLOC,PNDGLN((ELE-1)*MLOC+KLOC):',&
                  &                   ELE,KLOC,PNDGLN((ELE-1)*MLOC+KLOC)
             do ILOC=1,4! Was loop 
                INOD=NDGLNO((ELE-1)*NLOC+ILOC)
                ewrite(3,*) 'ILOC,INOD:',ILOC,INOD
             END DO
             FLAbort("STOP 2392")
          ENDIF
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE MIDSNODS

end module assnav_module
