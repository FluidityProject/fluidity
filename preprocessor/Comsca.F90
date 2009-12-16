
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
module comsca_module
  use FLDebug
  use quadrature
  use elements
  use Fields
  use State_module
  use spud
  use global_parameters, only : FIELD_NAME_LEN, ident_name, phase2state_index, &
       state2phase_index, name_ident
  use legacy_boundary_conditions
  use legacy_field_lists
  use populate_state_module
  IMPLICIT NONE
  
  private
  
  public comsca, get_ntsol, get_nphase, get_ndisot
  
  contains

SUBROUTINE COMSCA(PARA,PROCNO,NPROCS, TOFIL,MXNTSO,MXNPHA, &
                                ! For the INTEGERS ...
     & NPHASE,NTSOL, IRAD,                     NWICEL,NLOC  ,NGI   ,MLOC,&
     & ITINOI, &
     & PFRONT,NCOLOP,QFRONT,NCOLOQ,            SNLOC, SNGI,&
     & DUMSAV,BINARY,MATRPT,VERSIO,            NPRESS,NPROPT,PORMED,INVPOT,&
     & FBMODL,PCENT, RADISO,OUTSTA,            EVTIOP,EVSBCY,EVMXOP,SMLOC,&
     & GEOBAL,OPTSOU,            ISPHERE,SIGMOP,MISNIT,&
                                ! Phase momentum equations...
     & DISOPT,DISOP2,QMISOU,            MULPA, CGSOLQ,GMRESQ,&
     & CGSQ  ,MIXMAS,UZAWA,                    POISON,PROJEC,&
     & TMPERI,&
     & PHAMOD,                                 MPHAOP,&
     & COMLAW,ADOPTU,ADOPTV,ADOPTW,            EQNSTA,PREOPT,NDISOP,NSUBVLOC,&
     & NSUBNVLOC,&
                                ! Field equation...
     & DISOTT,TPHASE,TPHAOP,                   TMULPA,CGSOLT,&
     & GMREST,CGST,TNOIT1,TNOIT2,              TLEFT1,TLEFT2,TIMM,TMINT1,&
     & TMINT2,TMISOU,                          TFRONT,TTPERI,IDENT,&
     & SCATVE,TELEDI,ADOPTT,                   IDEFRA,IDEHCO,IDEPRE,TANOMU,&
     & TSIGMO,NSUBTLOC,NSUBBOYLOC,&
     & THOPT,NDISOT,MIXEOS,MIXOCP,             MIXOKA,MIXODI,MIXOLA,MIXOFO,&
     & MIXG0O,MIXCHO,&
                                ! Radiation...
     & NOITSG,NOITG1,NOITG2,PRETYP,            MINTER,GLEFTP,NGROUP,NPAT,&
     & NMAT,  NMATD, NSCAT, NDELAY,            NSMAT, NITEIG,&
     & RADTEM,RADBUB,RADAIR,RADDIV,            RMISOU,NDFREE,&
                                ! This is for LOGICALS...
     & LOWQAD,D3,                              DCYL  ,DCYL25,NAV   ,MVMESH,&
     & ZERQG,             CMCHAN,PERIOD,HFUNCT,&
     & GETTAN,ROTAT, DSPH,  BHOUT,             RAD,   UVWADP,&
     & COGRAX,COGRAY,COGRAZ,                   ADMESH,&
                                ! Phase momentum equations...
     & LUMP  ,MAKSYM,                   CONVIS,PRESYM,MOMSYM,&
     & MASSYM,                                 ZERVEL,ZSOX  ,ZSOY,&
     & ZSOZ  ,CONSOX,CONSOY,CONSOZ,            CONDEN,CONMU,&
     & CONDIV,VARDIV,                          &
     & CHADEN,                                 ONEMU, TWOMU, ALLMU,&
     & XABSZE,YABSZE,ZABSZE,                   XABSCO,YABSCO,ZABSCO,&
     & ABSLUM,ABSLUP,SLUMP,                    COMPRE,&
                                ! Field equation...
     & BOUSIN,DEFALT,TLUMP,             ADV,&
     & DIFF,                                   TSYM  ,CONTEM,&
     & ZSOT,CONSOT,TCONDE,TCONMU,              TONEMU,TALLMU,TEQADP,  &
     & SUFTEM,TABSZE,TABSCO,&
                                ! Radiation...
     & GETVPT,GETTAB,TIME,ORTHOG,              TNOLIN,GSCOUP,PSIADP,&
                                ! For the REALS ...
     & ACCTIM,LTIME ,DT,&
     & ITHETA,ITIERR,STEDER,AUTACC,            AUTGOB,&
     & CPULIM,                   &
     & R0,D0,                                  ROTDIR,&
     & GRAVTY,CRIUP, MINCH, MAXCH,             CPUMES,MESTP1,MESTP2,&
     & XMAPAK,STIMDU,MXCHAN,VOILIM,            RESCOE,RESCOW,MUWALL,GASWDR,& 
     & ALFSTA,ALFST2,SPRESS,                   MAXKE, MINKE, SPHERI,EVTIR1,&
     & EVTIR2,EVGDIN,WATIME,            &
                                ! Phase momentum equations...
     & THETA ,BETA  ,DENGAM,&
     & RONSOX,RONSOY,                          RONSOZ,RDENP ,&
     & RMUPXX,RMUPXY,RMUPXZ,RMUPYY,            RMUPYZ,RMUPZZ,&
     & TEMINI,DENINI,                          BSOUX,BSOUY,BSOUZ,&
     & XABS,  YABS,  ZABS,                     ADWEIU,ADWEIV,ADWEIW,&
     & TMPER1,TMPER2,                          MPHANO,&
     & GASCON,HEACAP,GAMDE2,GAMDE3,            ADATOU,ADATOV,ADATOW,&
     & MXVDRA,PERCO1,PERCO2,PERCO3,            PERCO4,&
                                ! Field equation...
     & TTHETA,TBETA,TEROR1,                    TEROR2,TRELAX,&
     & RONTEM,                                 RONSOT,TDEN,&
     & TMUXX,TMUXY,TMUXZ,                      TMUYY,TMUYZ,TMUZZ,&
     & VOLDEN,VOLMU,TABS,                      TTPER1,TTPER2,&
     & TPHANO,ADWEIT,                          RSCATV,&
     & ADATOT,                                 TPORD1,TPORD2,THPAR1,THPAR2,&
     & THPAR3,MIXCO1,MIXCO2,MIXCO3,            MIXCO4,&
     & MIXCP1,MIXCP2,MIXCP3,MIXCP4,            MIXKA1,MIXKA2,MIXKA3,MIXKA4,&
     & MIXDI1,MIXDI2,MIXDI3,MIXDI4,            MIXCP5,MIXCP6,MIXCP7,MIXCP8,&
     & MIXLA1,MIXLA2,MIXLA3,MIXLA4,            MIXFO1,MIXFO2,MIXFO3,MIXFO4,&
     & MIXG01,MIXG02,MIXG03,MIXG04,            MIXG05,MIXG06,MIXG07,MIXG08,&
     & MIXCH1,MIXCH2,MIXCH3,MIXCH4,  &
                                ! Radiation...
     & ERRORG,ERROG1,ERROG2,GRELAX,            ERREIG,DENMIX,BUBZER,&
     & XDIFFU,NUCLEA,PATMOS,                   TCV,   MULIQU,&
     & MUMIX, HEAPFI,CONPFI,CONBAK,            EXPCOE,STADEN,STATEM,&
     & have_state,state,size_state)
  ! This sub communicates the common variables from processor one to 
  ! all other processes. 
  ! If HFUNCT then solve for a height function. 
  INTEGER ONE
  PARAMETER(ONE=1)
  INTEGER MXNTSO,MXNPHA
  
  ! Can die when all the I/O has been stripped from comsca
  integer, parameter :: redcha = -666
  
  CHARACTER*240 NAME
  LOGICAL TOFIL
  CHARACTER*240 CHARA
  INTEGER PARA,PROCNO,NPROCS
  ! Blank spaces in file for future use...
  REAL RBLANK
  INTEGER IBLANK
  LOGICAL LBLANK
  !
  ! For the INTEGERS ...
  INTEGER &
       & NPHASE,NTSOL, IRAD,  &
       & NWICEL,NLOC  ,NGI   ,MLOC,&
       & ITINOI,&
       & PFRONT,NCOLOP,QFRONT,NCOLOQ,&
       & SNLOC, SNGI,&
       & DUMSAV,BINARY,MATRPT,VERSIO,&
       & NPRESS,NPROPT,PORMED,INVPOT,&
       & FBMODL,PCENT, RADISO,OUTSTA,&
       & EVTIOP,EVSBCY,EVMXOP,SMLOC,&
       & GEOBAL,OPTSOU,ISPHERE,SIGMOP,MISNIT,&
                                ! Phase momentum equations...&
       & DISOPT(MXNPHA),DISOP2(MXNPHA),QMISOU(MXNPHA),&
       & MULPA(MXNPHA), CGSOLQ(MXNPHA),GMRESQ(MXNPHA),&
       & CGSQ(MXNPHA)  ,MIXMAS(MXNPHA),UZAWA(MXNPHA),&
       & POISON(MXNPHA),PROJEC(MXNPHA),&
       & TMPERI(MXNPHA),&
       & PHAMOD(MXNPHA),&
       & MPHAOP(MXNPHA,MXNPHA),&
       & COMLAW(MXNPHA),ADOPTU(MXNPHA),ADOPTV(MXNPHA),ADOPTW(MXNPHA),&
       & EQNSTA(MXNPHA),PREOPT(MXNPHA),NDISOP(MXNPHA),NSUBVLOC(MXNPHA),&
       & NSUBNVLOC(MXNPHA),&
                                ! Field equation...
       & DISOTT(MXNTSO),TPHASE(MXNTSO),TPHAOP(MXNTSO),&
       & TMULPA(MXNTSO),CGSOLT(MXNTSO),&
       & GMREST(MXNTSO),CGST(MXNTSO),TNOIT1(MXNTSO),TNOIT2(MXNTSO),&
       & TLEFT1(MXNTSO),TLEFT2(MXNTSO),TIMM(MXNTSO),TMINT1(MXNTSO),&
       & TMINT2(MXNTSO),TMISOU(MXNTSO),&
       & TFRONT(MXNTSO),TTPERI(MXNTSO),IDENT(MXNTSO),&
       & SCATVE(MXNTSO),TELEDI(MXNTSO),ADOPTT(MXNTSO),&
       & IDEFRA(MXNTSO),IDEHCO(MXNTSO),IDEPRE(MXNTSO),TANOMU(MXNTSO),&
       & TSIGMO(MXNTSO),NSUBTLOC(MXNTSO),NSUBBOYLOC(MXNTSO),&
       & THOPT(MXNTSO),NDISOT(MXNTSO),MIXEOS(MXNTSO),MIXOCP(MXNTSO),&
       & MIXOKA(MXNTSO),MIXODI(MXNTSO),MIXOLA(MXNTSO),MIXOFO(MXNTSO),&
       & MIXG0O(MXNTSO),MIXCHO(MXNTSO),&
                                ! Radiation...
       & NOITSG,NOITG1,NOITG2,PRETYP,&
       & MINTER,GLEFTP,NGROUP,NPAT,&
       & NMAT,  NMATD, NSCAT, NDELAY,&
       & NSMAT, NITEIG,&
       & RADTEM,RADBUB,RADAIR,RADDIV,&
       & RMISOU,NDFREE
  !
  ! This is for LOGICALS...
  LOGICAL &
       & LOWQAD,D3,&
       & DCYL  ,DCYL25,NAV   ,MVMESH,&
       & ZERQG,&
       & CMCHAN,PERIOD,HFUNCT,&
       & GETTAN,ROTAT, DSPH,  BHOUT,&
       & RAD,   UVWADP,&
       & COGRAX,COGRAY,COGRAZ,&
       & ADMESH,&
                                ! Phase momentum equations...
       & LUMP(MXNPHA)  ,MAKSYM(MXNPHA), &
       & CONVIS(MXNPHA),PRESYM(MXNPHA),MOMSYM(MXNPHA),&
       & MASSYM(MXNPHA),&
       & ZERVEL(MXNPHA),ZSOX(MXNPHA)  ,ZSOY(MXNPHA),&
       & ZSOZ(MXNPHA)  ,CONSOX(MXNPHA),CONSOY(MXNPHA),CONSOZ(MXNPHA),&
       & CONDEN(MXNPHA),CONMU(MXNPHA),&
       & CONDIV(MXNPHA),VARDIV(MXNPHA),&
       & CHADEN(MXNPHA),&
       & ONEMU(MXNPHA), TWOMU(MXNPHA), ALLMU(MXNPHA),&
       & XABSZE(MXNPHA),YABSZE(MXNPHA),ZABSZE(MXNPHA),&
       & XABSCO(MXNPHA),YABSCO(MXNPHA),ZABSCO(MXNPHA),&
       & ABSLUM(MXNPHA),ABSLUP(MXNPHA),SLUMP(MXNPHA),&
       & COMPRE(MXNPHA),&
                                ! Field equation...
       & BOUSIN(MXNTSO),DEFALT(MXNTSO),TLUMP(MXNTSO),&
       & ADV(MXNTSO),&
       & DIFF(MXNTSO),&
       & TSYM(MXNTSO)  ,CONTEM(MXNTSO),&
       & ZSOT(MXNTSO),CONSOT(MXNTSO),TCONDE(MXNTSO),TCONMU(MXNTSO),&
       & TONEMU(MXNTSO),TALLMU(MXNTSO),TEQADP(MXNTSO),  &
       & SUFTEM(MXNTSO),TABSZE(MXNTSO),TABSCO(MXNTSO),&
                                ! Radiation...
       & GETVPT,GETTAB,TIME,ORTHOG,&
       & TNOLIN,GSCOUP,PSIADP
  !
  ! For the REALS ...
  REAL &
       & ACCTIM,LTIME ,DT,   &
       & ITHETA,ITIERR,STEDER,AUTACC,&
       & AUTGOB,&
       & CPULIM,&
       & R0,D0,&
       & ROTDIR(3),&
       & GRAVTY,CRIUP, MINCH, MAXCH,&
       & CPUMES,MESTP1,MESTP2,&
       & XMAPAK,STIMDU,MXCHAN,VOILIM,&
       & RESCOE,RESCOW,MUWALL,GASWDR, &
       & ALFSTA,ALFST2,SPRESS,&
       & MAXKE, MINKE, SPHERI,EVTIR1,&
       & EVTIR2,EVGDIN,WATIME,&
                                ! Phase momentum equations...&
       & THETA(MXNPHA) ,BETA(MXNPHA)  ,DENGAM(MXNPHA),&
       & RONSOX(MXNPHA),RONSOY(MXNPHA),&
       & RONSOZ(MXNPHA),RDENP(MXNPHA),&
       & RMUPXX(MXNPHA),RMUPXY(MXNPHA),RMUPXZ(MXNPHA),RMUPYY(MXNPHA),&
       & RMUPYZ(MXNPHA),RMUPZZ(MXNPHA),&
       & TEMINI(MXNPHA),DENINI(MXNPHA),&
       & BSOUX(MXNPHA),BSOUY(MXNPHA),BSOUZ(MXNPHA),&
       & XABS(MXNPHA),  YABS(MXNPHA),  ZABS(MXNPHA),&
       & ADWEIU(MXNPHA),ADWEIV(MXNPHA),ADWEIW(MXNPHA),&
       & TMPER1(MXNPHA),TMPER2(MXNPHA),&
       & MPHANO(MXNPHA,MXNPHA),&
       & GASCON(MXNPHA),HEACAP(MXNPHA),GAMDE2(MXNPHA),GAMDE3(MXNPHA),&
       & ADATOU(MXNPHA),ADATOV(MXNPHA),ADATOW(MXNPHA),&
       & MXVDRA(MXNPHA),PERCO1(MXNPHA),PERCO2(MXNPHA),PERCO3(MXNPHA),&
       & PERCO4(MXNPHA),&
                                ! Field equation...
       & TTHETA(MXNTSO),TBETA(MXNTSO),TEROR1(MXNTSO),&
       & TEROR2(MXNTSO),TRELAX(MXNTSO),&
       & RONTEM(MXNTSO),&
       & RONSOT(MXNTSO),TDEN(MXNTSO),&
       & TMUXX(MXNTSO),TMUXY(MXNTSO),TMUXZ(MXNTSO),&
       & TMUYY(MXNTSO),TMUYZ(MXNTSO),TMUZZ(MXNTSO),&
       & VOLDEN(MXNTSO),VOLMU(MXNTSO),TABS(MXNTSO),&
       & TTPER1(MXNTSO),TTPER2(MXNTSO),&
       & TPHANO(MXNTSO),ADWEIT(MXNTSO),&
       & RSCATV(MXNTSO),&
       & ADATOT(MXNTSO),&
       & TPORD1(MXNTSO),TPORD2(MXNTSO),THPAR1(MXNTSO),THPAR2(MXNTSO),&
       & THPAR3(MXNTSO),MIXCO1(MXNTSO),MIXCO2(MXNTSO),MIXCO3(MXNTSO),&
       & MIXCO4(MXNTSO),&
       & MIXCP1(MXNTSO),MIXCP2(MXNTSO),MIXCP3(MXNTSO),MIXCP4(MXNTSO),  &
       & MIXKA1(MXNTSO),MIXKA2(MXNTSO),MIXKA3(MXNTSO),MIXKA4(MXNTSO),&
       & MIXDI1(MXNTSO),MIXDI2(MXNTSO),MIXDI3(MXNTSO),MIXDI4(MXNTSO),&
       & MIXCP5(MXNTSO),MIXCP6(MXNTSO),MIXCP7(MXNTSO),MIXCP8(MXNTSO),&
       & MIXLA1(MXNTSO),MIXLA2(MXNTSO),MIXLA3(MXNTSO),MIXLA4(MXNTSO),&
       & MIXFO1(MXNTSO),MIXFO2(MXNTSO),MIXFO3(MXNTSO),MIXFO4(MXNTSO),&
       & MIXG01(MXNTSO),MIXG02(MXNTSO),MIXG03(MXNTSO),MIXG04(MXNTSO),            &
       & MIXG05(MXNTSO),MIXG06(MXNTSO),MIXG07(MXNTSO),MIXG08(MXNTSO),&
       & MIXCH1(MXNTSO),MIXCH2(MXNTSO),MIXCH3(MXNTSO),MIXCH4(MXNTSO),   &
                                ! Radiation...
       & ERRORG,ERROG1,ERROG2,GRELAX,&
       & ERREIG,DENMIX,BUBZER, &
       & XDIFFU,NUCLEA,PATMOS, &
       & TCV,   MULIQU, &
       & MUMIX, HEAPFI,CONPFI,CONBAK,&
       & EXPCOE,STADEN,STATEM
  INTEGER IP,IP2,IT,K,KK
  ! integer to save index of field for which bsines should be called
  integer bsines_it
  !
  !     state/options stuff -- cjc
  logical, intent(in) :: have_state
  integer, intent(in) :: size_state
  type(state_type), intent(in) :: state(size_state)
  type(mesh_type), pointer :: U_mesh,P_mesh
  integer :: stat ! Status variable for option routines.
  character(len=4000) :: thisphase, thisfield, tmpstring
  integer :: tmpint ! Temporary variable for grabing options for
  integer :: i, counter, tmpindex
  real, dimension(:), allocatable :: tmprealvector
  real :: tmpreal
  real :: omega
  ! processing. 
  logical :: have_velocity, have_pressure
  integer :: istate
  
  ewrite(1,*) 'SUBROUTINE COMSCA'

  RBLANK=0.
  IBLANK=0
  LBLANK=.FALSE.

  IF(tofil) then
     FLAbort("TOFIL must be false")
  end IF
  if(.not.have_state) then
     FLAbort("A state must be provided")
  end if
  
  ewrite(2,*) 'READING IN OPTIONS'
  ! For the integers ...
  if(have_state) then
     !FLAbort('need to uncomment commented options calls before continuing') 
     !'
     call get_nphase(nphase)
     call get_ntsol(ntsol)
     irad = 0         !hardwired pending fix from Jefferson
  else
     READ(REDCHA,111) NPHASE,NTSOL, IRAD
  end if
  if(have_state) then
     U_mesh => extract_mesh(state(1),'VelocityMesh')
     P_mesh => extract_mesh(state(1),'PressureMesh')
     call get_option("/geometry/dimension", tmpint)
     allocate(tmprealvector(tmpint))
     
     counter = option_count('/material_phase/&
          &vector_field::Velocity/prognostic')
     if (counter>0) then
        have_velocity=.true. 
     else
        have_velocity=.false. 
     end if
     
     counter = option_count('/material_phase/&
          &scalar_field::Pressure/prognostic')
     if (counter==1) then
        have_pressure=.true. 
     elseif (counter==0) then
        have_pressure=.false.
     elseif (counter>1) then
        FLExit('Not sure legacy compatibility can cope with multiple pressures.')
     end if
     
     NWICEL = 0
     !if(U_mesh%shape%dim.ne.3) then
     !   FLAbort('need to do something for non3d')
     !end if
     if((U_mesh%shape%loc==8).and.(P_mesh%shape%loc==1)) then
        nwicel = 1
     elseif((U_mesh%shape%loc==20).and.(P_mesh%shape%loc==8)) then
        nwicel = 2
     elseif((U_mesh%shape%loc==27).and.(P_mesh%shape%loc==8)) then
        nwicel = 3
     elseif((U_mesh%shape%loc==4).and.(P_mesh%shape%loc==1)) then  
        nwicel = 4
     elseif((U_mesh%shape%loc==4).and.(P_mesh%shape%loc==4)) then  
        nwicel = 5
     end if
     nloc = U_mesh%shape%loc
     ngi  = U_mesh%shape%quadrature%ngi
     mloc = P_mesh%shape%loc
  else
     READ(REDCHA,111) NWICEL,NLOC  ,NGI   ,MLOC
  end if

  if(have_state) then
     call get_option('/timestepping/nonlinear_iterations',ITINOI,&
          & default=1)
  else
     READ(REDCHA,111) ITINOI
  end if
  if(have_state) then
     PFRONT=0 ! Waiting for SCHEMA
     NCOLOP=0 ! Don't know if this is in an appropriate place but needs to be
     ! set to 0 so that successive later integers can be added to it
     if (have_option('/mesh_adaptivity/mesh_movement/full_ale')) then
        NCOLOP = NCOLOP + 10
     endif
     
     QFRONT=0 ! Waiting for SCHEMA
     NCOLOQ=0 ! Waiting for SCHEMA
  else
     READ(REDCHA,111) PFRONT,NCOLOP,QFRONT,NCOLOQ
  end if
  if(have_state) then
     snloc=U_mesh%faces%shape%loc
     sngi=U_mesh%faces%shape%quadrature%ngi
  else
     READ(REDCHA,111) SNLOC, SNGI
  end if
  if (have_state) then
     if (have_option("/io/max_dump_file_count")) then
        call get_option("/io/max_dump_file_count", DUMSAV)
     else
        DUMSAV=huge(0) ! Effectively no limit.
     end if
     BINARY=0         ! Default from datflu
     MATRPT=0         ! Default from datflu apparently not used
     VERSIO=66666     ! Should no longer be needed
  else
     READ(REDCHA,111) DUMSAV,BINARY,MATRPT,VERSIO
  end if
  if (have_state) then
     NPRESS=1         ! Default from datflu MULTIPHASE?
     NPROPT=1         ! Default from datflu MULTIPHASE?
     PORMED=0         ! Default from datflu MULTIPHASE?
     INVPOT=0         ! Default from datflu MULTIPHASE
  else
     READ(REDCHA,111) NPRESS,NPROPT,PORMED,INVPOT
  end if
  IF((VERSIO.GE.7).AND.(NPHASE.GT.1)) THEN
     if (have_state) then
        FBMODL = 1    ! Default from datflu MULTIPHASE?
        PCENT = 0     ! Default from datflu MULTIPHASE?
        RADISO = 0    ! Default from datflu MULTIPHASE?
        OUTSTA = 0    ! Default from datflu 
     else
        READ(REDCHA,111) FBMODL,PCENT, RADISO,OUTSTA
     end if
     if (have_state) then
        EVTIOP = 0    ! Default from datflu apparently not used
        EVSBCY = 0    ! Default from datflu apparently not used
        EVMXOP = 0    ! Default from datflu RADIATION
     else
        READ(REDCHA,111) EVTIOP,EVSBCY,EVMXOP,IBLANK
        READ(REDCHA,111) IBLANK,IBLANK,IBLANK,IBLANK
        READ(REDCHA,111) IBLANK,IBLANK,IBLANK,IBLANK
     end if
  ELSE
     FBMODL=1
  ENDIF
  SMLOC=P_mesh%faces%shape%loc
  
  tmpindex = 0
  do i = 1,size(state)
     if (has_scalar_field(state(1), 'BalancePressure')) then
        if (tmpindex==0) then
           tmpindex = i
        else
           FLExit("More than one material_phase with a BalancePressure.  Legacy compatability doesn't know how to deal with this.")
        end if
     end if
  end do
  
  if (tmpindex/=0) then
     call get_option('/material_phase['//int2str(tmpindex-1)//']/scalar_field::BalancePressure&
          &/prognostic/spatial_discretisation/&
          &geostrophic_pressure_option',&
          tmpstring)
     ! only quadratic balance pressure is known to work
     ! so geobal -20/-21    -Stephan
     if (tmpstring=='include_buoyancy') then
        GEOBAL=-21
     else if (tmpstring=='exclude_buoyancy') then
        GEOBAL=-20
     else
        FLAbort("Unknown geostrophic pressure option.")
     end if
  else if (has_scalar_field(state(1), 'FreeSurface')) then
     ! only known-to-work f.s. option at the moment -Stephan
     GEOBAL=-23
  else
     GEOBAL=0
  endif
  
  if(have_option('/traffic_model'))then
     optsou = 2
  else      
     optsou=0
  end if
  isphere=0
  if(have_option("/geometry/spherical_earth/&
      &linear_mapping")) then
    isphere=1
  else if(have_option("/geometry/spherical_earth/&
      &quadratic_superparametric_mapping")) then
    isphere=2
  end if
  sigmop = 0 ! Default from datflu believed deprecated.
  misnit = 0 ! Default from datflu
        
  DO IP=1,MAX(1,NPHASE)
     if (have_state) then
        istate = phase2state_index(ip)
        write(thisphase,'(a, i0, a)') "/material_phase[",istate-1,"]"
        if (have_option(trim(thisphase)//&
             "/vector_field::Velocity/prognostic/&
             &spatial_discretisation/&
             &legacy_continuous_galerkin")) then
           if (have_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/&
                &legacy_continuous_galerkin/&
                &balancing_diffusion_x")) then
              DISOPT(IP)=1
           else if (have_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/&
                &legacy_continuous_galerkin/&
                &laxwendroff_balancing_diffusion")) then
              DISOPT(IP)=2
           else if (have_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/&
                &legacy_continuous_galerkin/&
                &balancing_diffusion_x-t")) then
              DISOPT(IP)=3
           else if (have_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/&
                &legacy_continuous_galerkin/&
                &no_balancing_diffusion")) then
              DISOPT(IP)=4
           else if (have_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/&
                &legacy_continuous_galerkin/&
                &nonlinear_streamline_w_crossstream_diffusion")) then
              DISOPT(IP)=5
           else if (have_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/&
                &legacy_continuous_galerkin/&
                &nonlinear_upwind_steepest")) then
              DISOPT(IP)=6
           else if (have_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/&
                &legacy_continuous_galerkin/&
                &nonlinear_streamline_w_restricted_crossstream_diffusion")) then
              DISOPT(IP)=7
           else if (have_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/&
                &legacy_continuous_galerkin/&
                &les_constant_length_scale")) then
              DISOPT(IP)=42
           else if (have_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/&
                &legacy_continuous_galerkin/&
                &les_isotropic_length_scale")) then
              DISOPT(IP)=43
           else if (have_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/&
                &legacy_continuous_galerkin/&
                &les_no_balancing_diffusion")) then
              DISOPT(IP)=44
           else if (have_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/&
                &legacy_continuous_galerkin/&
                &les_no_balancing_diffusion_2")) then
              DISOPT(IP)=45
           else if (have_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/&
                &legacy_continuous_galerkin/&
                &les_no_balancing_diffusion_fourth_order_dissipation")) then
              DISOPT(IP)=46
           else if (have_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/&
                &legacy_continuous_galerkin/&
                &les_tensor_form")) then
              DISOPT(IP)=47
           else if (have_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/&
                &legacy_continuous_galerkin/&
                &les_fourth_order")) then
              DISOPT(IP)=48
           else if (have_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/&
                &legacy_continuous_galerkin/&
                &no_balancing_diffusion_remove_nonlinear_terms")) then
              DISOPT(IP)=125
           else
              FLExit("No known legacy_continuous_galerkin discretisation option given.")
           end if
        else if (have_option(trim(thisphase)//&
             "/vector_field::Velocity/prognostic/&
             &spatial_discretisation/&
             &legacy_mixed_cv_cg")) then
           call get_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/&
                &legacy_mixed_cv_cg/&
                &legacy_disott", DISOPT(IP), default=4)
        else
           call get_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/&
                &legacy_discretisation/&
                &legacy_disopt", DISOPT(IP), default=0)
        end if
        QMISOU(IP)=0 !Deprecated PETSC?
        DISOP2(IP)=0 !Deprecated PETSC?
     else
        QMISOU(IP)=0
        READ(REDCHA,111) DISOPT(IP),DISOP2(IP),QMISOU(IP)
     end if
     if (have_state) then
        MULPA(IP)=0 ! Deprecated PETSC
        CGSOLQ(IP)=0 ! Deprecated PETSC
        GMRESQ(IP)=1 ! Deprecated PETSC
     else
        READ(REDCHA,111) MULPA(IP), CGSOLQ(IP),GMRESQ(IP)
     end if
     if (have_state) then
        CGSQ(IP)=0 ! Deprecated PETSC
        MIXMAS(IP)=0 !  Legacy pressure solver?
        UZAWA(IP)=0 ! Legacy pressure solver?
     else
        READ(REDCHA,111) CGSQ(IP)  ,MIXMAS(IP),UZAWA(IP)
     end if
     if (have_state) then
        if (have_pressure) then
           call get_option(trim(thisphase)//&
                "/scalar_field::Pressure/prognostic/&
                &scheme/poisson_pressure_solution", tmpstring)
           select case (tmpstring)
           case ("never")
              POISON(IP)=0
           case("every timestep")
              POISON(IP)=1
           case ("only first timestep")
              POISON(IP)=-1
           case default
              FLExit(trim(tmpstring)//" is not a legal poisson_pressure_solution") 
           end select
        else
           POISON(IP)=-66602
        end if
        if ((have_option(trim(thisphase)//&
             "/scalar_field::Pressure/prognostic/&
             &scheme/use_projection_method")).or.&
             (have_option(trim(thisphase)//&
             "/scalar_field::Pressure/prognostic/&
             &scheme/use_compressible_projection_method"))) then
           PROJEC(IP)=1
        else
           PROJEC(IP)=0
        end if

     else
        READ(REDCHA,111) POISON(IP),PROJEC(IP)
     end if
     if (have_state) then
        TMPERI(IP) = -66612 ! Hard coded time periodic boundary conditions
        ! are a dumb idea.
        ! In the new options BCs are given by field value not field acceleration.
     else
        READ(REDCHA,111) TMPERI(IP)
     end if
     if (have_state) then
        PHAMOD = 0 ! Dunno what this is.
     else
        READ(REDCHA,111) PHAMOD(IP)
     end if
     if (have_state) then
        MPHAOP = 0 ! MULTIPHASE.
     else
        READ(REDCHA,111) (MPHAOP(IP,IP2),IP2=1,MAX(1,NPHASE))
     end if
     IF(VERSIO.GE.3) THEN
        if (have_state) then
           COMLAW(IP)=0 ! MULTIPHASE
           if (have_option(trim(thisphase)//& 
                "/vector_field::Velocity/prognostic/adaptivity_options&
                &/absolute_measure_velocity")) then
              ADOPTU(IP)=0
              ADOPTV(IP)=0
              ADOPTW(IP)=0
           end if
        else
           READ(REDCHA,111) COMLAW(IP),ADOPTU(IP),ADOPTV(IP),ADOPTW(IP)
        end if
     ELSE
        COMLAW(IP)=0
     ENDIF
     IF(VERSIO.GE.11) THEN
        if (have_state) then
           if (.not.have_option(trim(thisphase)//"/equation_of_state"))&
                & then
              EQNSTA(IP)=0
           else if (have_option(trim(thisphase)//"/equation_of_state/fluids"))&
                & then
              if (have_option(trim(thisphase)//"/equation_of_state/fluids/linear")) then
                 if(have_option(trim(thisphase)//"/equation_of_state/fluids/linear/salinity_dependency")) then
                    EQNSTA(IP)=1
                 else
                    EQNSTA(IP)=0
                 end if
              else if (have_option(trim(thisphase)//"/equation_of_state/&
                   &fluids/ocean_pade_approximation")) then
                 if(have_option(trim(thisphase)//"/equation_of_state/&
                      &fluids/ocean_pade_approximation/&
                      &include_depth_below_surface")) then
                    EQNSTA(IP)=2
                 else
                    EQNSTA(IP)=3
                 end if
              end if
           else if (have_option(trim(thisphase)//"/equation_of_state/&
                &multimaterial")) then
              ! Multimaterial does not use eqnsta.
              EQNSTA(IP)=0
           else 
              FLAbort("Missing equation of state for "//trim(thisphase))
           end if

           if (have_option(trim(thisphase)//&
                "/scalar_field::Pressure/prognostic/&
                &spatial_discretisation/legacy_continuous_galerkin/&
                &integrate_continuity_by_parts")) then
              PREOPT(IP)=1
           else
              PREOPT(IP)=0
           end if
           if (have_velocity) then

              call get_ndisot(trim(thisphase)//'vector_field::Velocity', ndisop(ip))
           else
              NDISOP(IP)=-66613
           end if
           if (have_option(trim(thisphase)//&
                "/vector_field::Velocity/prognostic/&
                &spatial_discretisation/inner_element")) then
              NSUBVLOC(IP)=4
              NSUBNVLOC(IP)=4
              if (have_option(trim(thisphase)//&
                   "/vector_field::Velocity/prognostic/&
                   &spatial_discretisation/inner_element/&
                   &use_quadratic_pressure")) then
                 PREOPT(IP)=PREOPT(IP)+1000
                 ewrite(1,*)'Inner Element model - using *quadratic* pressure (preopt+1e3)'
              else
                 ewrite(1,*)'Inner Element model - using *linear* pressure'
              end if
              if (have_option(trim(thisphase)//&
                   "/vector_field::Velocity/prognostic/&
                   &spatial_discretisation/inner_element/&
                   &apply_full_discontinuous_Galerkin")) then
                 PREOPT(IP)=PREOPT(IP)+10000
                 ewrite(1,*)'Inner Element model - using full DG (preopt+1e4)'
              else
                 ewrite(1,*)'Inner Element model - not using full DG'
              end if
           else
              NSUBVLOC(IP)=0
              NSUBNVLOC(IP)=0
           end if
        else
           READ(REDCHA,111) EQNSTA(IP),PREOPT(IP),NDISOP(IP),NSUBVLOC(IP)
        end if
        if (.not.have_state) then
           READ(REDCHA,111) NSUBNVLOC(IP),IBLANK,   IBLANK,   IBLANK
           READ(REDCHA,111) IBLANK    ,IBLANK,   IBLANK,   IBLANK
           READ(REDCHA,111) IBLANK    ,IBLANK,   IBLANK,   IBLANK
        end if
        !            ewrite(3,*) 'before eqnsta(ip)=',eqnsta(ip)
        IF(EQNSTA(IP).EQ.0) EQNSTA(IP)=COMLAW(IP) ! COMLAW appears
        ! not to be set - dham.
        !            ewrite(3,*) 'after eqnsta(ip)=',eqnsta(ip)
     ELSE
        EQNSTA(IP)=COMLAW(IP) 
        PREOPT(IP)=0
        NDISOP(IP)=0
        NSUBVLOC(IP)=0
        NSUBNVLOC(IP)=0
     ENDIF
  end do ! IP=1,MAX(1,NPHASE)

  DO IT=1,NTSOL
     !     INFORMATION FOR DAVID HAMBURGER
     !     The name of field(it) is stored in field_name_list(it)
     if (have_state) then
        write(thisphase,'(a, i0, a)') "/material_phase["&
             &,field_state_list(IT)-1,"]" 
        write (thisfield, '(a)') trim(thisphase)//"/scalar_field::"&
             &//trim(field_name_list(IT))
        if (have_option(trim(thisfield)//"/prognostic/&
             &spatial_discretisation/legacy_mixed_cv_cg")) then
           call get_option(trim(thisfield)//"/prognostic/spatial_discretisation/&
                &legacy_mixed_cv_cg/legacy_disott",DISOTT(IT), default=-7) ! Default
        else if (have_option(trim(thisfield)//"/prognostic/&
             &spatial_discretisation/legacy_continuous_galerkin")) then
           if (have_option(trim(thisfield)//"/prognostic/&
                &spatial_discretisation/legacy_continuous_galerkin/petrof_galerkin_xt")) then
              DISOTT(IT)=1
           else if (have_option(trim(thisfield)//"/prognostic/&
                &spatial_discretisation/legacy_continuous_galerkin/weighted_petrof_galerkin_xt")) then
              DISOTT(IT)=2
           else if (have_option(trim(thisfield)//"/prognostic/&
                &spatial_discretisation/legacy_continuous_galerkin/petrof_galerkin_x")) then
              DISOTT(IT)=3
           else if (have_option(trim(thisfield)//"/prognostic/&
                &spatial_discretisation/legacy_continuous_galerkin/least_squares_xt")) then
              DISOTT(IT)=4
           else if (have_option(trim(thisfield)//"/prognostic/&
                &spatial_discretisation/legacy_continuous_galerkin/least_squares_x")) then
              DISOTT(IT)=5
           else if (have_option(trim(thisfield)//"/prognostic/&
                &spatial_discretisation/legacy_continuous_galerkin/least_squares_xtheta")) then
              DISOTT(IT)=6
           else if (have_option(trim(thisfield)//"/prognostic/&
                &spatial_discretisation/legacy_continuous_galerkin/galerkin_theta")) then
              DISOTT(IT)=7
           else if (have_option(trim(thisfield)//"/prognostic/&
                &spatial_discretisation/legacy_continuous_galerkin/space_weighting")) then
              DISOTT(IT)=8
           end if
           if (have_option(trim(thisfield)//"/prognostic/&
                &spatial_discretisation/legacy_continuous_galerkin/optimally_upwind")) then
              DISOTT(IT)=-1*DISOTT(IT)
           end if
        else 
           call get_option(trim(thisfield)//"/prognostic/spatial_discretisation/&
                &legacy_discretisation/legacy_disott",DISOTT(IT), default=0) ! Default
           ! is for non-prognostic fields.
        end if
        tphase(it) = state2phase_index(field_state_list(it))
        TPHAOP(IT) = 0 ! Believed unused.

     else
        READ(REDCHA,111) DISOTT(IT),TPHASE(IT),TPHAOP(IT)
     end if
     if (have_state) then
        TMULPA(IT) = 0 ! Believed unused.
        CGSOLT(IT) = 0 ! Deprecated PETSC
     else
        READ(REDCHA,111) TMULPA(IT),CGSOLT(IT)
     end if
     if (have_state) then
        GMREST(IT) = 1 ! Deprecated PETSC
        CGST(IT) = 0 ! Deprecated PETSC
        call get_option(trim(thisfield)//"/prognostic/solver/legacy/&
             &max_iterations_outer",TNOIT1(IT), default=0) ! Default
        ! is for non-prognostic fields.
        call get_option(trim(thisfield)//"/prognostic/solver/legacy/&
             &max_iterations_inner",TNOIT2(IT), default=0) ! Default
        ! is for non-prognostic fields.
     else
        READ(REDCHA,111) GMREST(IT),CGST(IT),TNOIT1(IT),TNOIT2(IT)
     end if
     if (have_state) then
        TLEFT1(IT)=1 ! Deprecated PETSC
        TLEFT2(IT)=1 ! Deprecated PETSC
        TIMM(IT)=1   ! Deprecated PETSC
        call get_option(trim(thisfield)//"/prognostic/solver/legacy/&
             &max_nonrestarted_iterations_outer",TMINT1(IT), default=0)
        ! Default is for non-prognostic fields.             
     else
        READ(REDCHA,111) TLEFT1(IT),TLEFT2(IT),TIMM(IT),TMINT1(IT)
     end if
     if (have_state) then
        call get_option(trim(thisfield)//"/prognostic/solver/legacy/&
             &max_nonrestarted_iterations_inner",TMINT2(IT), default=0)
        ! Default is for non-prognostic fields.             
        call get_option(trim(thisfield)//"/prognostic/solver/legacy/&
             &field_no_outer_iterations_ddm",TMISOU(IT), default=0)
        ! Default is for non-prognostic fields.                          
     else
        READ(REDCHA,111) TMINT2(IT),TMISOU(IT)
     end if
     if(have_state) then
        TFRONT(IT) = 0 ! tfront -- legacy solver option PETSC
        TTPERI(IT) = -66614 ! ttperi -- time-dependent bc option 
        ! In the new options BCs are given by field value not field acceleration.
        ident(it) = name_ident(field_name_list(it))
     else
        READ(REDCHA,111) TFRONT(IT),TTPERI(IT),IDENT(IT)
     end if

     IF(VERSIO.GE.2) THEN
        if (have_state) then 
           SCATVE(IT) = 0 ! Waiting for SCHEMA
           TELEDI(IT) = 0 ! Waiting for someone to care!
           if (have_option(trim(thisfield)//"/prognostic/&
                &adaptivity_options/relative_measure_velocity")) then
              ADOPTT(IT) = 1
           else
              ADOPTT(IT) = 0
           end if
        else
           READ(REDCHA,111) SCATVE(IT),TELEDI(IT),ADOPTT(IT),  IBLANK
        end if
        !              IF((VERSIO.GE.8).AND.(NPHASE.GT.1)) THEN
        IF(((VERSIO.GE.8).AND.(NPHASE.GT.1)).OR.(VERSIO.GE.22)) THEN
           if (have_state) then
              IDEFRA(IT) = 0 ! Default from datflu MULTIPHASE
              IDEHCO(IT) = 0 ! Default from datflu MULTIPHASE
              IDEPRE(IT) = 0 ! Default from datflu MULTIPHASE
              counter = option_count(trim(thisfield)//"/prognostic&
                   &/tensor_field::Diffusivity/prescribed/value/anisotropic_symmetric")
              if(counter>0) then
                 TANOMU(IT) = 1 ! Use off diagonal diffusivities
              else
                 TANOMU(IT) = 0 ! Don't use off diagonal diffusivities
              end if
           else
              READ(REDCHA,111) IDEFRA(IT),IDEHCO(IT),IDEPRE(IT)&
                   &,TANOMU(IT)
           end if
           !idefra -- some kind of fraction id, default 0
           !idepre -- form or pressure term to add to compressible temp equation
           !          default 0
           !tanomu -- mystery option, always seems to be 0
           !         ewrite(3,*) 'it,TANOMU(IT):',it,TANOMU(IT)
           if (have_state) then
              if (have_option(trim(thisfield)//&
                   "/prognostic/&
                   &spatial_discretisation/inner_element")) then
                 NSUBTLOC(IT) = 4
              else
                 NSUBTLOC(IT) = 0 ! Default from datflu MULTIPHASE
              end if
              TSIGMO(IT) = 0 ! Default from datflu MULTIPHASE
              NSUBBOYLOC(IT) = 0  ! Default from datflu MULTIPHASE
           else
              READ(REDCHA,111) TSIGMO(IT),NSUBTLOC(IT),NSUBBOYLOC(IT)&
                   &,IBLANK
              !tsigmo -- mystery, defaults to 0
              !nsubtloc -- number of degrees of freedom for subscale in T
              !nsuboyloc -- mystery, defaults to 0
              READ(REDCHA,111) IBLANK    ,IBLANK,    IBLANK,       IBLANK
              READ(REDCHA,111) IBLANK    ,IBLANK,    IBLANK,       IBLANK
              READ(REDCHA,111) IBLANK    ,IBLANK,    IBLANK,       IBLANK
              READ(REDCHA,111) IBLANK    ,IBLANK,    IBLANK,       IBLANK
           end if
        ELSE
           IDEFRA(IT)=0                 
           IDEHCO(IT)=0
           IDEPRE(IT)=0
           TANOMU(IT)=0
           TSIGMO(IT)=0
           NSUBTLOC(IT)=0
           NSUBBOYLOC(IT)=0
        ENDIF
     ELSE
        SCATVE(IT)=0
        TELEDI(IT)=0
        IDEFRA(IT)=0
        IDEHCO(IT)=0
        IDEPRE(IT)=0
        TANOMU(IT)=0
        TSIGMO(IT)=0
        NSUBTLOC(IT)=0
        NSUBBOYLOC(IT)=0
     ENDIF
     IF(VERSIO.GE.12) THEN
        if (have_state) then
           THOPT(IT) = 0 ! Unknown option.
           NDISOT(IT)=0

           call get_ndisot(trim(thisfield), ndisot(it))

           MIXEOS(IT)=0
           MIXOCP(IT) = 0 ! Waiting for SCHEMA heat capacity.
        else
           READ(REDCHA,111) THOPT(IT) ,NDISOT(IT),MIXEOS(IT)&
                &,MIXOCP(IT) 
        end if
        if (have_state) then
           MIXOKA(IT) = 0 ! Believed unused
           MIXODI(IT) = 0 ! Believed unused
           MIXOLA(IT) = 0 ! MULTIPHASE
           MIXOFO(IT) = 0 ! Believed unused
        else
           READ(REDCHA,111) MIXOKA(IT),MIXODI(IT),MIXOLA(IT)&
                &,MIXOFO(IT) 
        end if
        if (have_state) then
           MIXG0O(IT) = 0 ! Believed unused
           MIXCHO(IT) = 0 ! Believed unused
        else
           READ(REDCHA,111) MIXG0O(IT),MIXCHO(IT),IBLANK,IBLANK
        end if
     ELSE
        THOPT(IT)=0
        NDISOT(IT)=0
        MIXEOS(IT)=0
        MIXOCP(IT)=0
        MIXOKA(IT)=0
        MIXODI(IT)=0
        MIXOLA(IT)=0
        MIXOFO(IT)=0
        MIXG0O(IT)=0
        MIXCHO(IT)=0
     ENDIF
  end do ! IT=1,NTSOL
  IF(IRAD.NE.0) THEN
     if (have_option("/radiation")) then
        FLAbort("Radiation is not currently supported")
     end if
     !radiation options, I suggest we just bung them in one place
     READ(REDCHA,111) NOITSG,NOITG1,NOITG2,PRETYP
     READ(REDCHA,111) MINTER,GLEFTP,NGROUP,NPAT
     READ(REDCHA,111) NMAT,  NMATD, NSCAT, NDELAY
     READ(REDCHA,111) NSMAT, NITEIG
     READ(REDCHA,111) RADTEM,RADBUB,RADAIR,RADDIV
     READ(REDCHA,111) RMISOU,NDFREE
  ENDIF
  ! This is for logicals
  if (have_state) then
     LOWQAD=.false. ! Believed deprecated
     call get_option("/geometry/dimension", tmpint)
     D3=(tmpint==3)
  else
     READ(REDCHA,333) LOWQAD,D3
     !priscr -- printing option, usually false I think
     !wriscr -- ditto
     !lowqad -- seems not used, false
  end if
  if (have_state) then
     DCYL=.false. ! Cylindrical coordinates not currently supported.
     DCYL25=.false.
!!$           if (.not.D3) then
!!$              FLAbort("2D simulations are not currently supported")
!!$           end if
     NAV=.true. ! Waiting for SCHEMA
     MVMESH=have_option("/mesh_adaptivity/mesh_movement")
  else        
     READ(REDCHA,333) DCYL  ,DCYL25,NAV   ,MVMESH
     !dcyl -- set to false and assert that dimension of mesh is 3
     !dcyl25 - ditto
     !nav -- not in schema, should it be?
  end if
  ZERQG=.true. ! Always start mesh movement from 0.
  if (have_state) then
     CMCHAN=.false.
     do i = 1, size(state)
        if(have_option('/material_phase['//int2str(i-1)//'/&
             &scalar_field::Pressure/prognostic/scheme/&
             &update_discretised_equation')) then
           CMCHAN=.true.
        end if
     end do
     PERIOD=.false. 
     if(have_option('/geometry/mesh::VelocityMesh//&
          &from_mesh/periodic_boundary_conditions')) then
        PERIOD=.true.
     end if
     HFUNCT=.false. ! Believed unused.
  else
     READ(REDCHA,333) CMCHAN,PERIOD,HFUNCT
     !cmchan -- schema
     !period -- schema? if not, set to false and I will fix, cjc
     !hfunct -- dead option, set false
  end if
  if (have_state) then
     ! Look for rotated boundary conditions:
     call set_rotat_gettan(rotat, gettan)
     DSPH=.false. ! Believed unused

     BHOUT=.false.
     do i = 1,size(state)
        if(have_option('/material_phase['//int2str(i-1)//']/equation_of_state/&
             &fluids/linear/subtract_out_hydrostatic_level').or. &
             have_option('/material_phase['//int2str(i-1)//']/equation_of_state/&
             &fluids/ocean_pade_approximation')) then
           BHOUT=.true.
        end if
     end do
  else
     READ(REDCHA,333) GETTAN,ROTAT, DSPH,  BHOUT
     !gettan -- option to tell fluidity to construct tangents (usually
     !true)
     !rotat -- schema
     !dsph -- 3d spherical coordinates, maybe not used
     !bhout -- schema
  end if
  if (have_state) then
     rad=have_option("/radiation")
     if (rad) FLAbort("Radiation is not currently supported")
     UVWADP=.false. ! not used for anything anymore?
  else               
     READ(REDCHA,333) RAD,   UVWADP
     !rad -- radiation something, set false
     !uvwadp -- might have to ask about that one, could be important
  end if
  if (have_state) then
     COGRAX=.false.
     if(have_option("/physical_parameters/gravity/&
          &vector_field::GravityDirection/prescribed/value/python")&
          &.or. option_count("/physical_parameters/&
          &gravity/vector_field::GravityDirection/value") > 1) then
        COGRAX=.false.
     else
        COGRAX=.true.
     end if
     COGRAY=COGRAX
     COGRAZ=COGRAX
     ewrite(3,*) "COGRAX:", COGRAX
     ewrite(3,*) have_option("/physical_parameters/gravity/&
          &vector_field::GravityDirection/prescribed/value/python")
     ewrite(3,*) option_count("/physical_parameters/&
          &gravity/vector_field::GravityDirection/value")
  else
     READ(REDCHA,333) COGRAX,COGRAY,COGRAZ
     !constant gravity field
  end if
  if (have_state) then
     ADMESH=have_option("/mesh_adaptivity/hr_adaptivity")
  else
     READ(REDCHA,333) ADMESH
     READ(REDCHA,333) LBLANK,LBLANK,LBLANK,LBLANK
     READ(REDCHA,333) LBLANK,LBLANK,LBLANK,LBLANK
  end if
  DO IP=1,MAX(1,NPHASE)
     if (have_state) then
        istate = phase2state_index(ip)
        write(thisphase,'(a, i0, a)') "/material_phase[",istate-1,"]"

        LUMP(IP)=((have_option(trim(thisphase)//& 
             "/vector_field::Velocity/prognostic/&
             &spatial_discretisation/legacy_continuous_galerkin/lump_mass_matrix")).or.&
             (have_option(trim(thisphase)//& 
             "/vector_field::Velocity/prognostic/&
             &spatial_discretisation/discontinuous_galerkin/lump_mass_matrix")).or.&
             (have_option(trim(thisphase)//& 
             "/vector_field::Velocity/prognostic/&
             &spatial_discretisation/continuous_galerkin/mass_terms/lump_mass_matrix")).or.&
             (have_option(trim(thisphase)//& 
             "/vector_field::Velocity/prognostic/&
             &spatial_discretisation/legacy_discretisation/legacy_mlump")))


        MAKSYM(IP)=.false. ! Waiting for SCHEMA
     else
        READ(REDCHA,333) LUMP(IP)  ,MAKSYM(IP)
        !lump is lump_mass_matrix in schema
        !maksym -- bigm is symmetric
     end if
     if (have_state) then
        CONVIS(IP)=.true. ! Waiting for SCHEMA
        PRESYM(IP)=.false. ! Waiting for SCHEMA
        PRESYM(IP)=.false. ! Waiting for SCHEMA
     else
        READ(REDCHA,333) CONVIS(IP),PRESYM(IP),MOMSYM(IP)
        !convis -- conservative form of viscosity
     end if
     if (have_state) then
        MASSYM(IP) = have_option(trim(thisphase)//& 
             "/vector_field::Velocity/prognostic/&
             &spatial_discretisation/symmetric_mass_matrix")
     else
        READ(REDCHA,333) MASSYM(IP)
     end if
     if (have_state) then
        ZERVEL(IP)=.false. ! Initial conditions now happen in
        ! populate_state
        zsox(IP)=(.not.(have_option(trim(thisphase)//& 
             "/vector_field::Velocity/prognostic/&
             &vector_field::Source")))
        zsoy(IP)=zsox(IP)
     else
        READ(REDCHA,333) ZERVEL(IP),ZSOX(IP)  ,ZSOY(IP)
     end if
     !zervel -- zero velocity initially
     !zsox/y/z -- schema
     if (have_state) then
        zsoz(IP)=zsox(IP)
        CONSOX(IP)=.false.
        CONSOY(IP)=consox(ip) ! Believed unused
        CONSOZ(IP)=consox(ip) ! Believed unused
     else
        READ(REDCHA,333) ZSOZ(IP)  ,CONSOX(IP),CONSOY(IP),CONSOZ(IP)
        !consox/y/z appears not to be used, but check as could be 
        !to do with consistent mass which is somewhat important
     end if
     if (have_state) then
        CONDEN(IP)=.false. ! Believed unused
        CONMU(IP)=.false.
     else
        READ(REDCHA,333) CONDEN(IP),CONMU(IP)
        !ditto conden
     end if
     if (have_state) then
        CONDIV(IP)=.false. ! Believed unused
        VARDIV(IP)=.false. ! Unknown.
     else
        READ(REDCHA,333) CONDIV(IP),VARDIV(IP)
        !condiv as conden
        !vardiv means divergence varies, have div field
     end if
     if (have_state) then
        CHADEN(IP) = have_option(trim(thisphase)//"/equation_of_state/&
             &fluids/linear/legacy_bousinesq_change_density")
     else
        READ(REDCHA,333) CHADEN(IP)
     end if
     if (have_state) then
        ONEMU(IP)=.false.
        TWOMU(IP)=.false.
        ALLMU(IP)=.false.
        if((have_option(trim(thisphase)//"/vector_field::Velocity/&
             &prognostic/tensor_field::Viscosity"))) then
           counter = (option_count(trim(thisphase)//"/vector_field::Velocity/&
                &prognostic/tensor_field::Viscosity/prescribed/value/isotropic"))
           if((counter>0).and.&
                (counter==option_count(trim(thisphase)//"/vector_field::Velocity/&
                &prognostic/tensor_field::Viscosity/prescribed/value"))) then
              ONEMU(IP)=.true.
           else
              ALLMU(IP)=.true.
           end if
        end if
        counter = option_count(trim(thisphase)//"/vector_field::Velocity/&
             &prognostic/tensor_field::Viscosity/prescribed/value/anisotropic_asymmetric")
        if (counter>0) then
           ewrite(0,*) "warning: unless your code is written in the new coding structure &
                &you'll lose asymmetric viscosity components"
        end if
     else
        READ(REDCHA,333) ONEMU(IP), TWOMU(IP), ALLMU(IP)
     end if
     if (have_state) then
        XABSZE(IP)=.not.have_option(trim(thisphase)//& 
             "/vector_field::Velocity/prognostic/&
             &vector_field::Absorption/")
        YABSZE(IP)=XABSZE(IP)
        ZABSZE(IP)=XABSZE(IP)
     else
        READ(REDCHA,333) XABSZE(IP),YABSZE(IP),ZABSZE(IP)
        !have zero u/v/w absorption
     end if
     if (have_state) then
        ! We don't do initial conditions like this any more
        XABSCO(IP)=.false.
        YABSCO(IP)=.false.
        ZABSCO(IP)=.false.
     else
        READ(REDCHA,333) XABSCO(IP),YABSCO(IP),ZABSCO(IP)
        !constant value to set u/v/w absorption
     end if
     if (have_state) then
        ABSLUM(IP)=have_option(trim(thisphase)//&
             "/vector_field::Velocity/prognostic/&
             &vector_field::Absorption/lump_absorption")
        ABSLUP(IP)=.false. ! Believed unused
        SLUMP(IP)=have_option(trim(thisphase)//&
             "/vector_field::Velocity/prognostic/&
             &vector_field::Source/lump_source")
     else
        READ(REDCHA,333) ABSLUM(IP),ABSLUP(IP),SLUMP(IP)
        !abslum schema, abslum appears not to be used,slump means
        ! lump RHS
     end if
     IF(VERSIO.GE.1) THEN
        if (have_state) then
           compre(IP) = .false.
        else
           READ(REDCHA,333) COMPRE(IP),LBLANK,LBLANK,LBLANK
           !schema
        end if
     ELSE
        COMPRE(IP)=.FALSE.
     ENDIF
  end do

  bsines_it=0
  ! Looping over fields to find temperature and/or salinity 
  ! and associated equation of state options
  DO IT=1,NTSOL
     bousin(it)=.false.
     if (have_state) then
        if (trim(field_name_list(IT))=="Temperature" .and. &
             have_option(trim(thisphase)//"/equation_of_state/&
             &fluids/linear/temperature_dependency")) then
           bsines_it=it
        else if(trim(field_name_list(IT))=="Temperature" .and. &
             &have_option(trim(thisphase)//"/equation_of_state/&
             &fluids/ocean_pade_approximation")) then
           bsines_it=it
        else if (trim(field_name_list(IT))=="Salinity" .and. &
             have_option(trim(thisphase)//"/equation_of_state/&
             &fluids/linear/salinity_dependency")) then
           bsines_it=it
        else if(trim(field_name_list(IT))=="Salinity" .and. &
             &have_option(trim(thisphase)//"/equation_of_state/&
             &fluids/ocean_pade_approximation")) then
           bsines_it=it
        end if
     end if
  end DO
  ! If temperature and/or salinity affect the equation of state,
  ! set BOUSIN to be true for EITHER temperature or salinity -
  ! whichever is solved for last.
  if(bsines_it>0) then
     ewrite(2,*) "Setting BOUSIN to true for field with IT=", bsines_it
     ewrite(2,*) "This field is: ", trim(field_name_list(bsines_it))
     bousin(bsines_it)=.true.
  end if

  DO IT=1,NTSOL
      thisfield=field_optionpath_list(it)

      DEFALT(IT)=.false.  ! We don't do defaults like this anymore
      TLUMP(IT)=((have_option(trim(thisfield)//"/prognostic/&
           &spatial_discretisation/legacy_discretisation/legacy_tlump")).or.&
           (have_option(trim(thisfield)//"/prognostic/&
           &spatial_discretisation/legacy_continuous_galerkin/lump_mass_matrix")).or.&
           (have_option(trim(thisfield)//"/prognostic/&
           &spatial_discretisation/legacy_mixed_cv_cg")))
     if (have_state) then
        ADV(IT)=.true. ! We should be smarter in the SCHEMA
     else
        READ(REDCHA,333) ADV(IT)   
        !adv -- is advected
     end if
     if (have_state) then
        DIFF(IT)=.false. ! We don't do defaults like this anymore
        ! to solve a pure diffusion problem just have the velocity for
        ! that material_phase not set or set to 0
     else
        READ(REDCHA,333) DIFF(IT)
        !diff -- is diffused
     end if
     if (have_state) then
        TSYM(IT)=.false. ! Believed unused.
        CONTEM(IT)=.false. ! We don't do initial conditions like
        ! this anymore
     else
        READ(REDCHA,333) TSYM(IT)  ,CONTEM(IT)
        !TSYM seems not to be used
        !CONTEM -- temp set from constant
     end if
     if (have_state) then
        ZSOT(IT)=(.not.have_option(trim(thisfield)//&
             "/prognostic/scalar_field::Source"))
        CONSOT(IT)=.false. ! We don't do initial conditions like
        ! this anymore
        TCONDE(IT)=.true. ! Don't understand this.
        TCONMU=.false. ! We don't do initial conditions like
        ! this anymore
     else
        READ(REDCHA,333) ZSOT(IT),CONSOT(IT),TCONDE(IT),TCONMU(IT)
        !zsot schema, consot set source form constant
        !tconde -- constant density for temp (????)
        !tconmu -- schema
     end if
     if (have_state) then
        counter=option_count(trim(thisfield)//"/prognostic&
             &/tensor_field::Diffusivity/prescribed/value/isotropic")
        if(counter>0) then
           TONEMU(IT)=.true.
        else
           TONEMU(IT)=.false.
        end if
        counter=(option_count(trim(thisfield)//"/prognostic&
             &/tensor_field::Diffusivity/prescribed/value/anisotropic_symmetric")+&
             option_count(trim(thisfield)//"/prognostic&
             &/tensor_field::Diffusivity/prescribed/value/anisotropic_asymmetric"))
        if(counter>0) then
           TALLMU(IT)=.true.
        else
           TALLMU(IT)=.false.
        end if
        if(TALLMU(IT).and.TONEMU(IT)) then
           TONEMU(IT)=.false.
        end if
        counter=(option_count(trim(thisfield)//"/prognostic&
             &/tensor_field::Diffusivity/prescribed/value/anisotropic_asymmetric"))
        if (counter>0) then
           ewrite(0,*) "warning: unless your code is written in the new coding structure &
                &you'll lose asymmetric diffusivity components"
        end if
        TEQADP(IT)=.false. ! Not used for anything anymore?
     else
        READ(REDCHA,333) TONEMU(IT),TALLMU(IT),TEQADP(IT)   
        !schema,schema,something with adaptive timestepping
     end if
     if (have_state) then
        SUFTEM(IT)=have_option(trim(thisfield)//"/prognostic/&
             &discretisation/legacy_discretisation/legacy_suftem")               
        TABSZE(IT)=(.not.have_option(trim(thisfield)//&
             "/prognostic/scalar_field::Absorption"))
        TABSCO(IT)=.false. !We don't do initial conditions like that
        ! anymore.
     else
        READ(REDCHA,333) SUFTEM(IT),TABSZE(IT),TABSCO(IT)
        !schema,zero for absorption in t, set constant absorption
     end if
  end do
  IF(IRAD.NE.0) THEN
     if (have_state) then
        FLAbort("Radiation not supported in new options yet")
     end if
     !radiation stuff, bundle into own options
     READ(REDCHA,333) GETVPT,GETTAB,TIME,ORTHOG
     READ(REDCHA,333) TNOLIN,GSCOUP,PSIADP
  ENDIF
  ! For the reals ...
  if (have_state) then
     call get_option("/timestepping/current_time", ACCTIM)
     call get_option("/timestepping/finish_time", LTIME)
     call get_option("/timestepping/timestep", DT)
  else
     READ(REDCHA,222) ACCTIM,LTIME ,DT    
  end if
  if (have_state) then

     tmpreal = 0.0
     tmpindex = 0
     do i = 1,size(state)
        call get_option('/material_phase['//int2str(i-1)//']/vector_field::Velocity/&
             &prognostic/temporal_discretisation/relaxation', tmpreal, stat=stat)
        if(stat==0) then
           if(tmpindex==0) then
              tmpindex = i
              ITHETA = tmpreal
           else
              if(tmpreal==ITHETA) then
                 ewrite(0,*) 'Warning: two values of ITHETA not supported in legacy code structure'
              else
                 FLExit('Two different values of ITHETA have been set in different material_phases.  This is not supported by legacy code structure.')
              end if
           end if
        else
          ITHETA = 1.0
        end if
     end do
     call get_option("/timestepping/nonlinear_iterations/&
          &tolerance", ITIERR,default=0.0)
     call get_option("/timestepping/steady_state/tolerance", &
          STEDER, default = -666.01)
     call get_option("/timestepping/adaptive_timestep&
          &/tolerance", autacc, default=0.0)
  else
     READ(REDCHA,222) ITHETA,ITIERR,STEDER,AUTACC
  end if
  if (have_state) then
     call get_option("/timestepping/adaptive_timestep&
          &/adaptive_timestep_repeat", AUTGOB, default=0.0)
  else
     READ(REDCHA,222) AUTGOB
  end if
  
  if (have_state) then
     R0=0.0 ! Waiting for SCHEMA
     D0=0.0 ! Waiting for SCHEMA
  else
     READ(REDCHA,222) R0,D0
     !radius of earth?, some kind of general constant overloaded everywhere
  end if
  if (have_state) then
     ROTDIR(1)=0.0 ! Waiting for SCHEMA
     ROTDIR(2)=0.0
     ROTDIR(3)=0.0
  else
     READ(REDCHA,222) ROTDIR(1),ROTDIR(2),ROTDIR(3)
     !rotation direction for a vector
  end if
  if (have_state) then
     call get_option("/physical_parameters/gravity/magnitude", GRAVTY, &
          & default=0.0)
     CRIUP=0.0 ! Waiting for SCHEMA
     call get_option("/mesh_adaptivity/hr_adaptivity/&
          &constant_size_constraint/minimum_edge_length", MINCH,&
          & default=0.0)
     call get_option("/mesh_adaptivity/hr_adaptivity/&
          &constant_size_constraint/maximum_edge_length", MAXCH,&
          & default=0.0)
  else
     READ(REDCHA,222) GRAVTY,CRIUP, MINCH, MAXCH
     !schema, lower limit of mesh refinement
  end if
  if (have_state) then
     call get_option("/mesh_adaptivity/hr_adaptivity/&
          &cpu_period", CPUMES, default=0.0)
     call get_option("/mesh_adaptivity/hr_adaptivity/&
          &functional_tolerance", MESTP1, default=0.0)
     call get_option("/mesh_adaptivity/hr_adaptivity/&
          &maximum_aspect_ratio", MESTP2, default=0.0)
  else
     READ(REDCHA,222) CPUMES,MESTP1,MESTP2
  end if
  if (have_state) then
     XMAPAK=0.0
     STIMDU=0.0
     MXCHAN=0.0
     VOILIM=0.6 ! WTF?
  else
     READ(REDCHA,222) XMAPAK,STIMDU,MXCHAN,VOILIM
     !xmapak something odd to do with event
     !stimdu, time between diagnostic dumps
     !mxchan, some momentum exchange calculation
     !voilim seems obselete
  end if
  IF(((NPHASE.GT.1).OR.(VERSIO.GE.20))&
       &        .AND.(VERSIO.GE.4)) THEN
     !bunch of multiphase stuff, defaults below
     if (have_state) then
        RESCOE=0.995 ! MULTIPHASE?
        RESCOW=0.75 ! MULTIPHASE?
        MUWALL=0.2 ! MULTIPHASE?
        GASWDR=0.025 ! MULTIPHASE?
        ALFSTA=0.62 ! MULTIPHASE?
        ALFST2=10.
        call get_option("/timestepping/wall_time_limit", WATIME,&
             & default=0.0)
     else
        READ(REDCHA,222) RESCOE,RESCOW,MUWALL,GASWDR
        READ(REDCHA,222) ALFSTA,ALFST2,SPRESS,PATMOS
        READ(REDCHA,222) MAXKE, MINKE, SPHERI,EVTIR1
        READ(REDCHA,222) EVTIR2,EVGDIN,WATIME
     end if
  ELSE
     RESCOE=0.995
     RESCOW=0.75
     MUWALL=0.2
     GASWDR=0.025
     ALFSTA=0.62
     ALFST2=0.62
     SPRESS=10.
     WATIME=0.0
  ENDIF
          
  IF(VERSIO.LE.8) VOILIM=0.6
  DO IP=1,MAX(1,NPHASE)
     if (have_state) then
        istate = phase2state_index(IP)
        write(thisphase,'(a, i0, a)') "/material_phase[",istate-1,"]"
        if (have_velocity) then
           call get_option(trim(thisphase)//"/vector_field::Velocity/&
                &prognostic/temporal_discretisation/theta",&
                & THETA(IP), default=0.5)
           call get_option(trim(thisphase)//"/vector_field::Velocity/&
                &prognostic/spatial_discretisation/&
                &conservative_advection", BETA(IP), default=0.0)
        else
           THETA(IP)=666.02
           BETA(IP)=666.03
        end if
        call get_option(trim(thisphase)//&
             '/scalar_field::Pressure/prognostic/atmospheric_pressure',&
             PATMOS, default=0.0)
     else
        READ(REDCHA,222) THETA(IP) ,BETA(IP) 
        !schema,schema,schema,mystery probably legacy solve option
     end if
     if (have_state) then
        call get_option(trim(thisphase)//"/equation_of_state/&
             &fluids/linear/temperature_dependency/&
             &thermal_expansion_coefficient", DENGAM(IP), &
             default=-666.11)
     else
        READ(REDCHA,222) DENGAM(IP)
     end if
     if (have_state) then
        RONSOX(IP)=0.0 ! We don't do initial conditions like that
        ! anymore.
        RONSOY(IP)=0.0
     else
        READ(REDCHA,222) RONSOX(IP),RONSOY(IP)
        !value to initalise velocity sources as a constant
     end if
     if (have_state) then
        RONSOZ(IP)=0.0
        RDENP(IP)=0.0  !SCHEMA?
     else
        READ(REDCHA,222) RONSOZ(IP),RDENP(IP)
        ! same for w, and denp
     end if
     !        ewrite(3,*) '57777 DENGAM(1),rdenp(1):',DENGAM(1),rdenp(1)
     !        RDENP(1)=1.0
     if (have_state) then
        RMUPXX(IP)=0.0
        RMUPXY(IP)=0.0
        RMUPXZ(IP)=0.0
        RMUPYY(IP)=0.0
        RMUPYZ(IP)=0.0
        RMUPZZ(IP)=0.0
     else
        READ(REDCHA,222) RMUPXX(IP),RMUPXY(IP),RMUPXZ(IP),RMUPYY(IP)
        !constant to set viscosity
        READ(REDCHA,222) RMUPYZ(IP),RMUPZZ(IP)
     end if
     if (have_state) then
        call get_option(trim(thisphase)//"/equation_of_state/&
             &fluids/linear/temperature_dependency/&
             &reference_temperature", TEMINI(IP), &
             default=-666.12)
        call get_option(trim(thisphase)//"/equation_of_state/&
             &fluids/linear/reference_density", DENINI(IP), &
             default=-666.13)
     else
        READ(REDCHA,222) TEMINI(IP),DENINI(IP)
     end if
     if (have_state) then
        call get_option("/physical_parameters/gravity/&
             &vector_field::GravityDirection/prescribed/value[0]/constant",&
             & tmprealvector, stat=stat) ! back compatibility with BSOUX, BSOUY, BSOUZ only
        ! supports one constant direction so assume here
        ! that this will be provided by the first value of
        ! the prescribed gravity field
        if (stat/=0) then
           tmprealvector=0.0
        end if
        BSOUX(IP)=tmprealvector(1)
        if (size(tmprealvector)>1) then
           BSOUY(IP)=tmprealvector(2)
        else
           BSOUY(IP)=0.0
        end if
        if (size(tmprealvector)>2) then
           BSOUZ(IP)=tmprealvector(3)
        else
           BSOUZ(IP)=0.0
        end if
     else
        READ(REDCHA,222) BSOUX(IP),BSOUY(IP),BSOUZ(IP)
     end if
     if (have_state) then
        XABS(IP)=0.0 ! We don't do initial conditions like this anymore.
        YABS(IP)=0.0
        ZABS(IP)=0.0
     else
        READ(REDCHA,222) XABS(IP),  YABS(IP),  ZABS(IP)
        !values for constant absorption
     end if
     if (have_state) then
        call get_option(trim(thisphase)//"vector_field::Velocity/&
             &adaptivity_options/absolute_measure_velocity/&
             &desired_error", tmprealvector, stat=stat)
        if (stat/=0) then
           call get_option(trim(thisphase)//"vector_field::Velocity/&
                &adaptivity_options/relative_measure_velocity/&
                &desired_error", tmprealvector, stat=stat)
        end if
        if (stat/=0) then
           call get_option(trim(thisphase)//"vector_field::Velocity/&
                &adaptivity_options/absolute_measure_speed/&
                &desired_error", tmprealvector, stat=stat)
        end if
        if (stat/=0) then
           call get_option(trim(thisphase)//"vector_field::Velocity/&
                &adaptivity_options/relative_measure_speed/&
                &desired_error", tmprealvector, stat=stat)
        end if
        if (stat/=0) then
           tmprealvector=0.0
        end if
        ADWEIU(IP)=tmprealvector(1)
        if (size(tmprealvector)>1) then
           ADWEIV(IP)=tmprealvector(2)
        end if
        if (size(tmprealvector)>2) then
           ADWEIW(IP)=tmprealvector(3)
        else
           ADWEIW(IP)=0.0
        end if
     else
        READ(REDCHA,222) ADWEIU(IP),ADWEIV(IP),ADWEIW(IP)
     end if
     if (have_state) then
        TMPER1(IP)=0.0
        TMPER2(IP)=0.0
     else
        READ(REDCHA,222) TMPER1(IP),TMPER2(IP)
     end if
     !schema
     if(have_state) then
        mphano(ip,:) = 1.0
     else
        READ(REDCHA,222) (MPHANO(IP,IP2),IP2=1,MAX(1,NPHASE))
     end if
     IF(VERSIO.GE.3) THEN
        if(have_state) then
           gascon(ip)=0.0
           heacap(ip) = 1.0 !legacy, set from default
           call get_option(trim(thisphase)//"/equation_of_state/&
                &fluids/linear/salinity_dependency/&
                &saline_contraction_coefficient", GAMDE2(IP),&
                &default=0.0)
           call get_option(trim(thisphase)//"/equation_of_state/&
                &fluids/linear/salinity_dependency/&
                &reference_salinity", GAMDE3(IP), default=0.0)
        else
           READ(REDCHA,222) GASCON(IP),HEACAP(IP),GAMDE2(IP),GAMDE3(IP)
        end if
        !schema,heat capacity, salinity constants (I think)
     ELSE
        GAMDE2(IP)=0.0
        GAMDE3(IP)=0.0
     ENDIF
     IF((VERSIO.GE.4).AND.ADMESH) THEN
        if (have_state) then
           call get_option(trim(thisphase)//"vector_field::Velocity/&
                &adaptivity_options/relative_measure_velocity/&
                &desired_error", tmprealvector, stat=stat)
           if (stat/=0) then
              call get_option(trim(thisphase)//"vector_field::Velocity/&
                   &adaptivity_options/absolute_measure_velocity/&
                   &desired_error", tmprealvector, stat=stat)
           end if
           if (stat/=0) then
              tmprealvector=0.0
           end if
           ADATOU(IP)=tmprealvector(1)
           if (size(tmprealvector)>1) then
              ADATOV(IP)=tmprealvector(2)
           else
              ADATOV(IP)=0
           end if
           if (size(tmprealvector)>2) then
              ADATOW(IP)=tmprealvector(3)
           else
              ADATOW(IP)=0.0
           end if
        else
           READ(REDCHA,222) ADATOU(IP),ADATOV(IP),ADATOW(IP),RBLANK
           READ(REDCHA,222) RBLANK,RBLANK,RBLANK,RBLANK
        end if
     ENDIF
     IF((VERSIO.GE.9).AND.(NPHASE.GT.1)) THEN
        if(have_state) then
           mxvdra(ip) = 1.0e+20
           perco1(ip) = 0.0
           perco2(ip) = 0.0
           perco3(ip) = 0.0
           perco4(ip) = 0.0
        else
           READ(REDCHA,222) MXVDRA(IP),PERCO1(IP),PERCO2(IP),PERCO3(IP)
           READ(REDCHA,222) PERCO4(IP),RBLANK,RBLANK,RBLANK
        end if
     ELSE
        MXVDRA(IP)=1.E+20
     ENDIF
  end do

  DO IT=1,NTSOL
     if(have_state) then
        write(thisphase,'(a, i0, a)') "/material_phase["&
             &,field_state_list(IT)-1,"]" 
        write (thisfield, '(a)') trim(thisphase)//"/scalar_field::"&
             &//trim(field_name_list(IT))

        call get_option(trim(thisfield)//'/prognostic/solver&
             &/legacy/tolerance_outer',teror1(it), default=0.0)
        if (field_name_list(it)=='BalancePressure' .or. &
             field_name_list(it)=='FreeSurface') then
           ttheta(it)=0.0
           tbeta(it)=0.0
        else
           call get_option(trim(thisfield)//'/prognostic/&
                &temporal_discretisation/theta', TTHETA(IT), default=1.0)
           call get_option(trim(thisfield)//'/prognostic/&
                &/spatial_discretisation/&
                &conservative_advection',tbeta(it), default=0.0)
        end if
     else
        READ(REDCHA,222) TTHETA(IT),TBETA(IT),TEROR1(IT)
     end if
     if(have_state) then
        teror2(it) = 1.e-10 !legacy solver option
        trelax(it) = 0.     !legacy solver option
     else
        READ(REDCHA,222) TEROR2(IT),TRELAX(IT)
     end if
     !schema,schema,legacy solver option
     if(have_state) then
        rontem(it) = 0.
     else
        READ(REDCHA,222) RONTEM(IT)
     end if
     !legacy solver option
     if(have_state) then
        ronsot(it) = 0.           !legacy solver option
        tden(it) = 1.0
     else
        READ(REDCHA,222) RONSOT(IT),TDEN(IT)
     end if
     if(have_state) then
        tmuxx(it) = 0.
        tmuyy(it) = 0.
        tmuzz(it) = 0.
     else
        READ(REDCHA,222) TMUXX(IT),TMUYY(IT),TMUZZ(IT)
     end if
     if(have_state) then
        volden(it) = 0.
        volmu(it) = 0.
        tabs(it) = 0.
     else
        READ(REDCHA,222) VOLDEN(IT),VOLMU(IT),TABS(IT)
     end if
     !volume fraction density constant value
     !another volume fraction option
     !constant value to initialise T absorption
     if(have_state) then              
        ttper1(it) = 0.
        ttper2(it) = 0.
     else
        READ(REDCHA,222) TTPER1(IT),TTPER2(IT)
     end if
     if(have_state) then
        TPHANO(IT) = 1.0 !legacy, default
        call get_option(trim(thisfield)//'/adaptivity_options/&
             &absolute_measure/desired_error',adweit(it), stat=stat)
        if (stat/=0) then
           call get_option(trim(thisfield)//'/adaptivity_options/&
                &relative_measure/desired_error',adweit(it), stat=stat)
        end if
     else
        READ(REDCHA,222) TPHANO(IT),ADWEIT(IT)
     end if
     IF(VERSIO.GE.2) THEN
        if(have_state) then
           RSCATV(IT) = 0. !not sure what that does
           TMUXY(IT)=0.0
           TMUXZ(IT)=0.0
           TMUYZ(IT)=0.0
        else
           READ(REDCHA,222) RSCATV(IT),TMUXY(IT),TMUXZ(IT),TMUYZ(IT)
        end if
     ELSE
        RSCATV(IT)=0.0
        TMUXY(IT)=0.0
        TMUXZ(IT)=0.0
        TMUYZ(IT)=0.0
     ENDIF
     IF((VERSIO.GE.4).AND.ADMESH) THEN
        if(.not. have_state) then
           READ(REDCHA,222) ADATOT(IT),RBLANK,RBLANK,RBLANK
           !schema
           READ(REDCHA,222) RBLANK,RBLANK,RBLANK,RBLANK
        end if
     ENDIF
     IF(VERSIO.GE.10) THEN
        if(have_state) then
        else
           READ(REDCHA,222) TPORD1(IT),TPORD2(IT),THPAR1(IT),THPAR2(IT)

           READ(REDCHA,222) THPAR3(IT),MIXCO1(IT),MIXCO2(IT),MIXCO3(IT)
           !mystery time dependent terms in heat equation
           READ(REDCHA,222) MIXCO4(IT),RBLANK,RBLANK,RBLANK
           !schema
        end if
     ENDIF
     IF(VERSIO.GE.19) THEN
        if(have_state) then
        else
           !these are all heat capacity coefficients
           READ(REDCHA,222) MIXCP1(IT),MIXCP2(IT),MIXCP3(IT),MIXCP4(IT)
           READ(REDCHA,222) MIXKA1(IT),MIXKA2(IT),MIXKA3(IT),MIXKA4(IT)
           READ(REDCHA,222) MIXDI1(IT),MIXDI2(IT),MIXDI3(IT),MIXDI4(IT)
           READ(REDCHA,222) MIXCP5(IT),MIXCP6(IT),MIXCP7(IT),MIXCP8(IT)
           READ(REDCHA,222) MIXLA1(IT),MIXLA2(IT),MIXLA3(IT),MIXLA4(IT)
           READ(REDCHA,222) MIXFO1(IT),MIXFO2(IT),MIXFO3(IT),MIXFO4(IT)
        end if
     ENDIF
     !        ewrite(3,*) 'just before new bit it=',it
     IF(VERSIO.GE.20) THEN
        !same here
        if(have_state) then
        else
           READ(REDCHA,222) MIXG01(IT),MIXG02(IT),MIXG03(IT),MIXG04(IT)
           READ(REDCHA,222) MIXG05(IT),MIXG06(IT),MIXG07(IT),MIXG08(IT)
           READ(REDCHA,222) MIXCH1(IT),MIXCH2(IT),MIXCH3(IT),MIXCH4(IT)
           READ(REDCHA,222) RBLANK,RBLANK,RBLANK,RBLANK
           READ(REDCHA,222) RBLANK,RBLANK,RBLANK,RBLANK
           READ(REDCHA,222) RBLANK,RBLANK,RBLANK,RBLANK
        ENDIF
     end if
  end do
  IF(IRAD.NE.0) THEN
     if(have_state) then
     else
        !stick in radiation options
        READ(REDCHA,222) ERRORG,ERROG1,ERROG2,GRELAX
        READ(REDCHA,222) ERREIG,DENMIX,BUBZER
        IF(VERSIO.GE.3) THEN
           READ(REDCHA,222) XDIFFU,NUCLEA,PATMOS
        ELSE
           XDIFFU=50.0
           NUCLEA=0.05
           PATMOS=1.01325E+6
        ENDIF
        !more radiation stuff, use defaults below
        READ(REDCHA,222) TCV,   MULIQU
        READ(REDCHA,222) MUMIX, HEAPFI,CONPFI,CONBAK
        IF(VERSIO.GE.2) THEN
           READ(REDCHA,222) EXPCOE,STADEN,  STATEM,   RBLANK
           IF(VERSIO.GE.6) THEN
              STADEN=DENINI(1)
              STATEM=TEMINI(1)
           ENDIF
        ELSE
           EXPCOE=DENGAM(1)
        ENDIF
     end if
     ! endof IF(IRAD.NE.0) THEN...
  ENDIF
  !      ewrite(3,*)'in comsca read in reals' 
  !
  ewrite(3,*) 'PATMOS=',PATMOS
  ewrite(3,*) 'RDENP(1):',RDENP(1)
  !         STOP 242
  !
  IF(XMAPAK.EQ.0.0) XMAPAK=1.0
  IF(NPRESS.EQ.0) NPRESS=1
  !
  ! *****CHECK CONSISTENCY OF OPTIONS ************
  NPHASE=MAX(1,NPHASE)

111 FORMAT(10I9)
222 FORMAT(4E15.7)

333 FORMAT(4L7)


   !-- Inner Element setup ---------------------------------------------------------

   if ( have_state .and. &
        (any(NSUBVLOC(1:nphase)/=0) .or. any(NSUBTLOC(1:ntsol)/=0)) ) then
     ewrite(3,*)'Inner Element model detected'
     if (any(NSUBVLOC(1:nphase)/=0)) then
       ewrite(1,*)'Inner Element model applied to momentum transport'
       ewrite(3,*)'Inner Element model  nsubvloc[', NSUBVLOC,']'
       ewrite(3,*)'Inner Element model nsubnvloc[', NSUBNVLOC,']'
       ewrite(3,*)'Inner Element model    preopt[', PREOPT,']'
     end if
     if (any(NSUBTLOC(1:ntsol)/=0)) then
       ewrite(1,*)'Inner Element model applied to tracer transport'
       ewrite(3,*)'Inner Element model  nsubtloc[', NSUBTLOC,']'
     end if
   end if

   if (have_state) then
      do ip=1,max(1,nphase)
         istate = phase2state_index(ip)
      end do
   end if

   !-- Inner Element setup ---------------------------------------------------------

END SUBROUTINE COMSCA

subroutine get_ntsol(ntsol)
  integer, intent(out) :: ntsol
  integer :: nphases,nfields,ncars,p,f
  character(len=FIELD_NAME_LEN) :: tmpstring
  logical :: aliased, pressure, density

  ntsol = 0

  nphases = option_count('/material_phase')  
  do p = 0, nphases-1
     nfields = option_count('/material_phase[' &
          //int2str(p)//']/scalar_field')
     do f = 0, nfields-1
        aliased = have_option('/material_phase['// &
             int2str(p)//']/scalar_field['//int2str(f)//']/aliased')
        call get_option('/material_phase['// &
             int2str(p)//']/scalar_field['//int2str(f)//']/name', tmpstring)
        pressure = (trim(tmpstring)=='Pressure')
        density = (trim(tmpstring)=='Density')

        if ((.not.aliased).and.(.not.pressure).and.(.not.density)) then
          ntsol = ntsol + 1
        end if
     end do
     ! prognostic scalar fields for Mellor Yamada:
     if (have_option('/material_phase[' &
          //int2str(p)//']/subgridscale_parameterisations/Mellor_Yamada/scalar_field::KineticEnergy/prognostic')) then
        ntsol=ntsol + 1
     end if
     if (have_option('/material_phase[' &
          //int2str(p)//']/subgridscale_parameterisations/Mellor_Yamada/scalar_field::TurbulentLengthScalexKineticEnergy/prognostic')) then
        ntsol=ntsol + 1
     end if
     if (have_option('/material_phase[' &
          //int2str(p)//']/subgridscale_parameterisations/GLS/scalar_field::GLSTurbulentKineticEnergy/prognostic')) then
        ntsol=ntsol + 1
     end if
     if (have_option('/material_phase[' &
          //int2str(p)//']/subgridscale_parameterisations/GLS/scalar_field::GLSGenericSecondQuantity/prognostic')) then
        ntsol=ntsol + 1
     end if
     ! Prognostic sediment fields.
     ntsol=ntsol+option_count('/material_phase[' &
          //int2str(p)//']/sediment/sediment_class')
  end do
  
  ! tracers for traffic modelling
  if(have_option('traffic_model'))then
     if(have_option('traffic_model/scalar_field::TrafficTracerTemplate'))then
        call get_option('/traffic_model/number_of_vehicles',ncars)
        ntsol=ntsol+ncars
     endif
  endif
  
end subroutine get_ntsol

subroutine get_nphase(nphase)
  integer, intent(out) :: nphase
  integer :: nmaterial_phases,p

  nphase = 0

  nmaterial_phases = option_count('/material_phase')  
  do p = 0, nmaterial_phases-1
     if (have_option('/material_phase['//int2str(p)//']/vector_field::Velocity')) then
        ! don't know if prescribed or diagnostic fields should be included in nphase but
        ! suspect that for things like traffic they should be
        ! definitely don't want aliased - crgw
        if (.not.have_option('/material_phase['//int2str(p)//']/vector_field::Velocity/aliased')) then
          nphase = nphase + 1
        end if
     end if
  end do
end subroutine get_nphase

  subroutine get_ndisot(thisfield, ndisot)
    character(len=*), intent(in) :: thisfield
    integer, intent(out) :: ndisot

    ! local
    integer :: tmpint
    character(len=FIELD_NAME_LEN) :: face_value

    face_value=""
    ndisot=0
    if (have_option(trim(thisfield)//"/prognostic/&
        &spatial_discretisation/legacy_mixed_cv_cg")) then
      call get_option(trim(thisfield)//"/prognostic/&
              &spatial_discretisation/legacy_mixed_cv_cg/&
              &face_value[0]/name", face_value)
      select case(face_value)
      case ( "FirstOrderUpwind" )
          if (have_option(trim(thisfield)//"/prognostic/&
            &temporal_discretisation/control_volumes/&
            &limit_theta")) then
            ndisot=1
          else
            ndisot=0
          end if
      case ( "Trapezoidal" )
          if (have_option(trim(thisfield)//"/prognostic/&
            &temporal_discretisation/control_volumes/&
            &limit_theta")) then
            ndisot=3
          else
            ndisot=2
          end if

      case ( "FiniteElement" )
          if (have_option(trim(thisfield)//"/prognostic/&
            &temporal_discretisation/control_volumes/&
            &limit_theta")) then
            if (have_option(trim(thisfield)//"/prognostic/&
                &spatial_discretisation/&
                &legacy_mixed_cv_cg/face_value[0]/limit_face_value")) then
                ndisot=5
            else
                ndisot=7
            end if
          else
            if (have_option(trim(thisfield)//"/prognostic/&
                &spatial_discretisation/&
                &legacy_mixed_cv_cg/face_value[0]/limit_face_value")) then
                ndisot=4
            else
                ndisot=6
            end if
          end if
      case ( "HyperC", "UltraC" )
          if (have_option(trim(thisfield)//"/prognostic/&
            &temporal_discretisation/control_volumes/&
            &limit_theta")) then
            ndisot=9
          else
            ndisot=8
          end if
      case default
        FLAbort( "Unknown control volume discretisation type." )
      end select

      call get_option(trim(thisfield)//"/prognostic/&
          &temporal_discretisation/&
          &control_volumes/number_advection_iterations",&
          tmpint, default=1)
      ndisot=ndisot+tmpint*10

      if (have_option(trim(thisfield)//"/prognostic/&
          &spatial_discretisation/&
          &legacy_mixed_cv_cg/store_upwind_elements")) then
          ndisot=ndisot+1000
      end if
    else
      call get_option(trim(thisfield)//"/prognostic/&
            &spatial_discretisation/legacy_discretisation/legacy_ndisot"&
            &,ndisot, default=0) ! Default 
                                 ! is for non-prognostic fields.
    end if

  end subroutine get_ndisot

end module comsca_module
