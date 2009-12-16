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

module redfil_module
  use AllSorts
  use FLDebug
  use state_module
  use fields
  use spud
  use flcomms_module
  use global_parameters, only : FIELD_NAME_LEN, phase2state_index, &
       state2phase_index
  use Legacy_Boundary_Conditions
  use parallel_tools
  use populate_state_module
  use comsca_module
  use tr2d_module
  use sml
  use redsca_module
  use legacy_field_lists
  
  implicit none

  private 

  public :: redfil, sndgln2tsndgl, &
     initialise_state_phase_lists_from_options, &
     initialise_state_phase_lists_from_file
     
  contains

SUBROUTINE REDFIL(&
     !     THE VARIABLES DEFINED IN COMSCA(1ST LINES IS USED IN)
     !     *****************************************************&
     MXNTSO,MXNPHA,&
     !     For the INTEGERS ...
     NPHASE,NTSOL,&
     NLOC  ,NGI   ,MLOC,&
     ITINOI,&
     NCOLOP,&
     SNLOC, SNGI,VERSIO,&
     NPRESS,NPROPT,&
     RADISO,&
     GEOBAL,OPTSOU,ISPHERE,&
     MISNIT,&
     !     Phase momentum equations...
     DISOPT,MULPA, CGSOLQ,GMRESQ,&
     MIXMAS,UZAWA,POISON,PROJEC,&
     EQNSTA,PREOPT,NDISOP,NSUBVLOC, NSUBNVLOC,&
     !     Field equation...
     DISOTT,TPHASE,CGSOLT,GMREST,&
     IDENT,&
     TELEDI,&
     NSUBTLOC,&
     NDISOT,&
     !     This is for LOGICALS...
     D3,DCYL  ,NAV   ,MVMESH,&
     CMCHAN,GETTAN,&
     ROTAT,DSPH,BHOUT,RAD,COGRAX,COGRAY,COGRAZ,ADMESH,&
     !     Phase momentum equations...
     LUMP  ,MAKSYM,CONVIS,&
     CHADEN, &
     ABSLUM,&
     SLUMP,COMPRE,&
     !     Field equation...
     BOUSIN,TLUMP,&
     SUFTEM,&
     !     For the REALS ...
     ACCTIM,LTIME ,DT,&
     ITHETA,ITIERR,STEDER,&
     R0,D0,GRAVTY,&
     ALFST2,SPRESS,&
     !     Phase momentum equations...
     THETA ,BETA  , &
     DENGAM,&
     TEMINI,DENINI,&
     BSOUX,BSOUY,BSOUZ,&
     GAMDE2,GAMDE3,&
     !     Field equation...
     TTHETA,TBETA,&
     !     THE FOLLOWING ARE USED IN REDSCA*********************&
     NONODS,XNONOD,TOTELE,FREDOP,&
     NOBCU,NOBCV,NOBCW,  NOBCT, NNODP,&
     NNODPP, NDPSET,&
     NNODRO, STOTEL,&
     NSOUPT,MXNSOU,FIESOU,&
                                !     THE REMAINING VARIABLES DEFINED IN THIS SUB***********&
     procno,halo_tag, halo_tag_p, &
     BCU1_mem,BCU2_mem,BCV1_mem,BCV2_mem,BCW1_mem,BCW2_mem, & 
     BCT1_mem,BCT2_mem,&
     nstate,state, uses_old_code_path)

  INTEGER MXNTSO,MXNPHA
  INTEGER OPTSOU,ISPHERE,SIGMOP,MISNIT

  !     For the INTEGERS ...
  INTEGER &
       NPHASE,NTSOL, IRAD,&
       NWICEL,NLOC  ,NGI   ,MLOC,&
       ITINOI,&
       PFRONT,NCOLOP,QFRONT,NCOLOQ,&
       SNLOC, SNGI,&
       DUMSAV,BINARY,MATRPT,VERSIO,&
       NPRESS,NPROPT,PORMED,INVPOT,&
       FBMODL,PCENT, RADISO,OUTSTA,&
       EVTIOP,EVSBCY,EVMXOP,SMLOC,&
       GEOBAL,&
                                !     Phase momentum equations...
       DISOPT(MXNPHA),DISOP2(MXNPHA),QMISOU(MXNPHA),&
       MULPA(MXNPHA), CGSOLQ(MXNPHA),GMRESQ(MXNPHA),&
       CGSQ(MXNPHA)  ,MIXMAS(MXNPHA),UZAWA(MXNPHA),&
       POISON(MXNPHA),PROJEC(MXNPHA),&
       TMPERI(MXNPHA),&
       PHAMOD(MXNPHA),&
       MPHAOP(MXNPHA,MXNPHA),&
       COMLAW(MXNPHA),ADOPTU(MXNPHA),ADOPTV(MXNPHA),ADOPTW(MXNPHA),&
       EQNSTA(MXNPHA),PREOPT(MXNPHA),NDISOP(MXNPHA),NSUBVLOC(MXNPHA),&
       NSUBNVLOC(MXNPHA),&
                                !     Field equation...
       DISOTT(MXNTSO),TPHASE(MXNTSO),TPHAOP(MXNTSO),&
       TMULPA(MXNTSO),CGSOLT(MXNTSO),&
       GMREST(MXNTSO),CGST(MXNTSO),TNOIT1(MXNTSO),TNOIT2(MXNTSO),&
       TLEFT1(MXNTSO),TLEFT2(MXNTSO),TIMM(MXNTSO),TMINT1(MXNTSO),&
       TMINT2(MXNTSO),TMISOU(MXNTSO),&
       TFRONT(MXNTSO),TTPERI(MXNTSO),IDENT(MXNTSO),&
       SCATVE(MXNTSO),TELEDI(MXNTSO),ADOPTT(MXNTSO),&
       IDEFRA(MXNTSO),IDEHCO(MXNTSO),IDEPRE(MXNTSO),TANOMU(MXNTSO),&
       TSIGMO(MXNTSO),NSUBTLOC(MXNTSO),NSUBBOYLOC(MXNTSO),&
       THOPT(MXNTSO),NDISOT(MXNTSO),MIXEOS(MXNTSO),MIXOCP(MXNTSO),&
       MIXOKA(MXNTSO),MIXODI(MXNTSO),MIXOLA(MXNTSO),MIXOFO(MXNTSO),&
       MIXG0O(MXNTSO),MIXCHO(MXNTSO),&
                                !     Radiation...
       NOITSG,NOITG1,NOITG2,PRETYP,&
       MINTER,GLEFTP,NGROUP,NPAT,&
       NMAT,  NMATD, NSCAT, NDELAY,&
       NSMAT, NITEIG,&
       RADTEM,RADBUB,RADAIR,RADDIV,&
       RMISOU,NDFREE
  !     
  !     This is for LOGICALS...
  LOGICAL &
       LOWQAD,D3,&
       DCYL  ,DCYL25,NAV   ,MVMESH,&
       ZERQG,&
       CMCHAN,PERIOD,HFUNCT,&
       GETTAN,ROTAT, DSPH,  BHOUT,&
       RAD,   UVWADP,&
       COGRAX,COGRAY,COGRAZ,&
       ADMESH,&
                                !     Phase momentum equations...
       LUMP(MXNPHA)  ,MAKSYM(MXNPHA),&
       CONVIS(MXNPHA),PRESYM(MXNPHA),MOMSYM(MXNPHA),&
       MASSYM(MXNPHA),&
       ZERVEL(MXNPHA),ZSOX(MXNPHA)  ,ZSOY(MXNPHA),&
       ZSOZ(MXNPHA)  ,CONSOX(MXNPHA),CONSOY(MXNPHA),CONSOZ(MXNPHA),&
       CONDEN(MXNPHA),CONMU(MXNPHA),&
       CONDIV(MXNPHA),VARDIV(MXNPHA),&
       CHADEN(MXNPHA),&
       ONEMU(MXNPHA), TWOMU(MXNPHA), ALLMU(MXNPHA),&
       XABSZE(MXNPHA),YABSZE(MXNPHA),ZABSZE(MXNPHA),&
       XABSCO(MXNPHA),YABSCO(MXNPHA),ZABSCO(MXNPHA),&
       ABSLUM(MXNPHA),ABSLUP(MXNPHA),SLUMP(MXNPHA),&
       COMPRE(MXNPHA),&
                                !     Field equation...
       BOUSIN(MXNTSO),DEFALT(MXNTSO),TLUMP(MXNTSO),&
       ADV(MXNTSO),&
       DIFF(MXNTSO),&
       TSYM(MXNTSO)  ,CONTEM(MXNTSO),&
       ZSOT(MXNTSO),CONSOT(MXNTSO),TCONDE(MXNTSO),TCONMU(MXNTSO),&
       TONEMU(MXNTSO),TALLMU(MXNTSO),TEQADP(MXNTSO),&
       SUFTEM(MXNTSO),TABSZE(MXNTSO),TABSCO(MXNTSO),&
                                !     Radiation...
       GETVPT,GETTAB,TIME,ORTHOG,&
       TNOLIN,GSCOUP,PSIADP
  !     
  !     For the REALS ...
  REAL &
       ACCTIM,LTIME ,DT,&
       ITHETA,ITIERR,STEDER,AUTACC,&
       AUTGOB,&
       CPULIM,&
       R0,D0,&
       ROTDIR(3),&
       GRAVTY,CRIUP, MINCH, MAXCH,&
       CPUMES,MESTP1,MESTP2,&
       XMAPAK,STIMDU,MXCHAN,VOILIM,&
       RESCOE,RESCOW,MUWALL,GASWDR,&
       ALFSTA,ALFST2,SPRESS,&
       MAXKE, MINKE, SPHERI,EVTIR1,&
       EVTIR2,EVGDIN,WATIME,&
                                !     Phase momentum equations...
       THETA(MXNPHA) ,BETA(MXNPHA)  ,&
       DENGAM(MXNPHA),&
       RONSOX(MXNPHA),RONSOY(MXNPHA),&
       RONSOZ(MXNPHA),RDENP (MXNPHA),&
       RMUPXX(MXNPHA),RMUPXY(MXNPHA),RMUPXZ(MXNPHA),RMUPYY(MXNPHA),&
       RMUPYZ(MXNPHA),RMUPZZ(MXNPHA),&
       TEMINI(MXNPHA),DENINI(MXNPHA),&
       BSOUX(MXNPHA),BSOUY(MXNPHA),BSOUZ(MXNPHA),&
       XABS(MXNPHA),  YABS(MXNPHA),  ZABS(MXNPHA),&
       ADWEIU(MXNPHA),ADWEIV(MXNPHA),ADWEIW(MXNPHA),&
       TMPER1(MXNPHA),TMPER2(MXNPHA),&
       MPHANO(MXNPHA,MXNPHA),&
       GASCON(MXNPHA),HEACAP(MXNPHA),GAMDE2(MXNPHA),GAMDE3(MXNPHA),&
       ADATOU(MXNPHA),ADATOV(MXNPHA),ADATOW(MXNPHA),&
       MXVDRA(MXNPHA),PERCO1(MXNPHA),PERCO2(MXNPHA),PERCO3(MXNPHA),&
       PERCO4(MXNPHA),&
                                !     Field equation...
       TTHETA(MXNTSO),TBETA(MXNTSO),TEROR1(MXNTSO),&
       TEROR2(MXNTSO),TRELAX(MXNTSO),&
       RONTEM(MXNTSO),&
       RONSOT(MXNTSO),TDEN(MXNTSO),&
       TMUXX(MXNTSO),TMUXY(MXNTSO),TMUXZ(MXNTSO),&
       TMUYY(MXNTSO),TMUYZ(MXNTSO),TMUZZ(MXNTSO),&
       VOLDEN(MXNTSO),VOLMU(MXNTSO),TABS(MXNTSO),&
       TTPER1(MXNTSO),TTPER2(MXNTSO),&
       TPHANO(MXNTSO),ADWEIT(MXNTSO),&
       RSCATV(MXNTSO),&
       ADATOT(MXNTSO),&
       TPORD1(MXNTSO),TPORD2(MXNTSO),THPAR1(MXNTSO),THPAR2(MXNTSO),&
       THPAR3(MXNTSO),MIXCO1(MXNTSO),MIXCO2(MXNTSO),MIXCO3(MXNTSO),&
       MIXCO4(MXNTSO),&
       MIXCP1(MXNTSO),MIXCP2(MXNTSO),MIXCP3(MXNTSO),MIXCP4(MXNTSO),& 
       MIXKA1(MXNTSO),MIXKA2(MXNTSO),MIXKA3(MXNTSO),MIXKA4(MXNTSO),&
       MIXDI1(MXNTSO),MIXDI2(MXNTSO),MIXDI3(MXNTSO),MIXDI4(MXNTSO),&
       MIXCP5(MXNTSO),MIXCP6(MXNTSO),MIXCP7(MXNTSO),MIXCP8(MXNTSO),&
       MIXLA1(MXNTSO),MIXLA2(MXNTSO),MIXLA3(MXNTSO),MIXLA4(MXNTSO),&
       MIXFO1(MXNTSO),MIXFO2(MXNTSO),MIXFO3(MXNTSO),MIXFO4(MXNTSO),&
       MIXG01(MXNTSO),MIXG02(MXNTSO),MIXG03(MXNTSO),MIXG04(MXNTSO),&
       MIXG05(MXNTSO),MIXG06(MXNTSO),MIXG07(MXNTSO),MIXG08(MXNTSO),&
       MIXCH1(MXNTSO),MIXCH2(MXNTSO),MIXCH3(MXNTSO),MIXCH4(MXNTSO),&
                                !     Radiation...
       ERRORG,ERROG1, ERROG2, GRELAX,&
       ERREIG,DENMIX,BUBZER,&
       XDIFFU,NUCLEA,PATMOS,&
       TCV,   MULIQU,&
       MUMIX, HEAPFI,CONPFI,CONBAK,&
       EXPCOE,STADEN,STATEM
  !     
  !     This is for distributed sources for a field...
  INTEGER NSOUPT,MXNSOU
  INTEGER FIESOU(MXNSOU)

  !     THE FOLLOWING ARE USED IN REDSCA*********************&
  INTEGER NONODS,XNONOD,TOTELE,FREDOP,&
       NOBCU(MXNPHA),NOBCV(MXNPHA),NOBCW(MXNPHA),&
       NOBCT(MXNTSO), NNODP,&
       NNODPP, NDPSET,&
       STOTEL, nnodro
  !     THE REMAINING VARIABLES DEFINED IN THIS SUB***********&
  INTEGER halo_tag,halo_tag_p
  type(real_vector), dimension(:), pointer :: bct1_mem
  type(integer_vector), dimension(:), pointer :: bct2_mem
  type(real_vector), dimension(:), pointer :: bcu1_mem
  type(integer_vector), dimension(:), pointer :: bcu2_mem
  type(real_vector), dimension(:), pointer :: bcv1_mem
  type(integer_vector), dimension(:), pointer :: bcv2_mem
  type(real_vector), dimension(:), pointer :: bcw1_mem
  type(integer_vector), dimension(:), pointer :: bcw2_mem

  !     THE FOLLOWING ARE ONLY USED IN THIS ROUTINE
  INTEGER IT,IP,NOPHAS,NOBCTS
  INTEGER NProcs,ProcNo,PARA

  integer, intent(in) :: nstate
  type(state_type), dimension(nstate),intent(in), target :: state
    
  logical, intent(in):: uses_old_code_path

  type(mesh_type), pointer :: U_mesh,X_mesh,P_mesh
  type(vector_field), pointer :: Velocity
  type(scalar_field), pointer :: Field_Buffer
  
  integer :: istate

  ! Temporary parameters, for additional sanity
  logical, parameter :: have_state = .true.
  logical, parameter :: remesh = .false.
  logical, parameter :: tofil = .false.
  
  ! The following additional assumptions have been made when stripping this
  ! routine (and are asserted to hold following comsca):
  !integer, parameter :: irad = 0
  !logical, parameter :: nav = .true.
  !logical, parameter :: rad = .false. 

  ! This sub extracts data from state and the options tree
  ! It also "allocates" memory in imem and rmem
  
  ewrite(1, *) "In redfil"

  
  call get_ntsol(ntsol)

  call initialise_field_lists_from_options(state, ntsol)
  call initialise_state_phase_lists_from_options()

  NProcs = GetNProcs()
  PARA=0
  IF(IsParallel()) THEN
     PARA = 1
  END IF
    
  ewrite(2, "(a,i0)") "Size of state is ", size(state)
  ewrite(1, *) "Calling comsca from redfil"
  CALL COMSCA(PARA,PROCNO,NPROCS,      TOFIL,MXNTSO,MXNPHA,&
                          !     For the INTEGERS ...
       NPHASE,NTSOL,IRAD,              NWICEL,NLOC,NGI,MLOC,&
       ITINOI, &
       PFRONT,NCOLOP,QFRONT,NCOLOQ,    SNLOC, SNGI,&
       DUMSAV,BINARY,MATRPT,VERSIO,    NPRESS,NPROPT,PORMED,INVPOT,&
       FBMODL,PCENT, RADISO,OUTSTA,    EVTIOP,EVSBCY,EVMXOP,SMLOC,&
       GEOBAL,OPTSOU,    ISPHERE,SIGMOP,MISNIT,&
                          !     Phase momentum equations...
       DISOPT,DISOP2,QMISOU,    MULPA, CGSOLQ,GMRESQ,&
       CGSQ  ,MIXMAS,UZAWA,            POISON,PROJEC,&
       TMPERI,&
       PHAMOD,                         MPHAOP,&
       COMLAW,ADOPTU,ADOPTV,ADOPTW,    EQNSTA,PREOPT,NDISOP,NSUBVLOC,&
       NSUBNVLOC,&
                          !     Field equation...
       DISOTT,TPHASE,TPHAOP,           TMULPA,CGSOLT,&
       GMREST,CGST,TNOIT1,TNOIT2,      TLEFT1,TLEFT2,TIMM,TMINT1,&
       TMINT2,TMISOU,                  TFRONT,TTPERI,IDENT,&
       SCATVE,TELEDI,ADOPTT,           IDEFRA,IDEHCO,IDEPRE,TANOMU,&
       TSIGMO,NSUBTLOC,NSUBBOYLOC, &
       THOPT,NDISOT,MIXEOS,MIXOCP,     MIXOKA,MIXODI,MIXOLA,MIXOFO,&
       MIXG0O,MIXCHO,&
                          !     Radiation...
       NOITSG,NOITG1,NOITG2,PRETYP,    MINTER,GLEFTP,NGROUP,NPAT,&
       NMAT,  NMATD, NSCAT, NDELAY,    NSMAT, NITEIG,&
       RADTEM,RADBUB,RADAIR,RADDIV,    RMISOU,NDFREE,&
                          !     This is for LOGICALS...
       LOWQAD,D3,&
       DCYL  ,DCYL25,NAV,MVMESH,       ZERQG,&
       CMCHAN,PERIOD,HFUNCT,           GETTAN,ROTAT, DSPH,  BHOUT,&
       RAD,   UVWADP,                  COGRAX,COGRAY,COGRAZ,&
       ADMESH,&
                          !     Phase momentum equations...
       LUMP  ,MAKSYM,           CONVIS,PRESYM,MOMSYM,&
       MASSYM,                         ZERVEL,ZSOX  ,ZSOY,&
       ZSOZ  ,CONSOX,CONSOY,CONSOZ,    CONDEN,CONMU,&
       CONDIV,VARDIV,                  &
       CHADEN,                         ONEMU,TWOMU,ALLMU,&
       XABSZE,YABSZE,ZABSZE,           XABSCO,YABSCO,ZABSCO,&
       ABSLUM,ABSLUP,SLUMP,            COMPRE,&
                          !     Field equation...
       BOUSIN,DEFALT,TLUMP,     ADV,&
       DIFF,                           TSYM  ,CONTEM,&
       ZSOT,CONSOT,TCONDE,TCONMU,      TONEMU,TALLMU,TEQADP,&
       SUFTEM,TABSZE,TABSCO,&
                          !     Radiation...
       GETVPT,GETTAB,TIME,ORTHOG,      TNOLIN,GSCOUP,PSIADP,&
                          !     For the REALS ...
       ACCTIM,LTIME ,DT,&
       ITHETA,ITIERR,STEDER,AUTACC,    AUTGOB,&
       CPULIM,           &
       R0,D0,                          ROTDIR,&
       GRAVTY,CRIUP, MINCH, MAXCH,     CPUMES,MESTP1,MESTP2,&
       XMAPAK,STIMDU,MXCHAN,VOILIM,    RESCOE,RESCOW,MUWALL,GASWDR,&
       ALFSTA,ALFST2,SPRESS,           MAXKE, MINKE, SPHERI,EVTIR1,&
       EVTIR2,EVGDIN,WATIME,    &
                          !     Phase momentum equations...
       THETA ,BETA  , DENGAM,&
       RONSOX,RONSOY,                  RONSOZ,RDENP ,&
       RMUPXX,RMUPXY,RMUPXZ,RMUPYY,    RMUPYZ,RMUPZZ,&
       TEMINI,DENINI,                  BSOUX,BSOUY,BSOUZ,&
       XABS,  YABS,  ZABS,             ADWEIU,ADWEIV,ADWEIW,&
       TMPER1,TMPER2,                  MPHANO,&
       GASCON,HEACAP,GAMDE2,GAMDE3,    ADATOU,ADATOV,ADATOW,&
       MXVDRA,PERCO1,PERCO2,PERCO3,    PERCO4,&
                          !     Field equation...
       TTHETA,TBETA,TEROR1,            TEROR2,TRELAX,&
       RONTEM,                         RONSOT,TDEN,&
       TMUXX,TMUXY,TMUXZ,              TMUYY,TMUYZ,TMUZZ,&
       VOLDEN,VOLMU,TABS,              TTPER1,TTPER2,&
       TPHANO,ADWEIT,                  RSCATV,&
       ADATOT,                         TPORD1,TPORD2,THPAR1,THPAR2,&
       THPAR3,MIXCO1,MIXCO2,MIXCO3,    MIXCO4,&
       MIXCP1,MIXCP2,MIXCP3,MIXCP4,    MIXKA1,MIXKA2,MIXKA3,MIXKA4,&
       MIXDI1,MIXDI2,MIXDI3,MIXDI4,    MIXCP5,MIXCP6,MIXCP7,MIXCP8,&
       MIXLA1,MIXLA2,MIXLA3,MIXLA4,    MIXFO1,MIXFO2,MIXFO3,MIXFO4,&
       MIXG01,MIXG02,MIXG03,MIXG04,    MIXG05,MIXG06,MIXG07,MIXG08,&
       MIXCH1,MIXCH2,MIXCH3,MIXCH4,&
                          !     Radiation...
       ERRORG,ERROG1,ERROG2,GRELAX,    ERREIG,DENMIX,BUBZER,&
       XDIFFU,NUCLEA,PATMOS,           TCV,MULIQU,&
       MUMIX, HEAPFI,CONPFI,CONBAK,    EXPCOE,STADEN,STATEM,&
       have_state,state,size(state))
  ewrite(1, *) "Exited comsca"
  
  assert(nav)
  assert(.not. rad)
  assert(irad == 0)

  NOPHAS=MAX(1,NPHASE)


  ewrite(1, *) "Calling redsca from redfil"
  CALL REDSCA(PROCNO,&
       NNODP, NNODPP, NDPSET,&
       OPTSOU,NSOUPT,MXNSOU,FIESOU)
  ewrite(1, *) "Exited redsca"

  U_mesh => extract_mesh(state(1),'VelocityMesh')
  nonods = node_count(U_mesh)
  nloc = U_mesh%shape%loc
  totele = element_count(U_mesh)
  X_mesh => extract_mesh(state(1),'CoordinateMesh')
  xnonod = node_count(X_mesh)
  stotel=surface_element_count(X_mesh)
  snloc=face_loc(X_mesh,1)
  sngi=face_ngi(X_mesh,1)
  P_mesh => extract_mesh(state(1),'PressureMesh')
  fredop = node_count(P_mesh)
  mloc = P_mesh%shape%loc

  if(isparallel()) then
    nnodp = get_nowned_nodes(halo_tag)
    nnodpp = get_nowned_nodes(halo_tag_p)
  else
    nnodp = nonods
    nnodpp = fredop
  end if

  TELEDI = 0

  !     Get the number of phases
  NOPHAS=MAX(1,NPHASE)

  ! rest of the subroutine is setting legacy boundary conditions:
  if (.not. uses_old_code_path) return

  if(associated( bcu1_mem )) deallocate( bcu1_mem )
  if(associated( bcu2_mem )) deallocate( bcu2_mem )
  if(nphase>0) then
     allocate( bcu1_mem(nphase) )
     allocate( bcu2_mem(nphase) )
     call nullify(bcu1_mem)
     call nullify(bcu2_mem)
  else
     allocate( bcu1_mem(1) )
     allocate( bcu2_mem(1) )
     allocate( bcu1_mem(1)%ptr(1) )
     allocate( bcu2_mem(1)%ptr(1) )
  end if

  if(associated( bcv1_mem )) deallocate( bcv1_mem )
  if(associated( bcv2_mem )) deallocate( bcv2_mem )
  if(nphase>0) then
     allocate( bcv1_mem(nphase) )
     allocate( bcv2_mem(nphase) )
     call nullify(bcv1_mem)
     call nullify(bcv2_mem)
  else
     allocate( bcv1_mem(1) )
     allocate( bcv2_mem(1) )
     allocate( bcv1_mem(1)%ptr(1) )
     allocate( bcv2_mem(1)%ptr(1) )
  end if

  if(associated( bcw1_mem )) deallocate( bcw1_mem )
  if(associated( bcw2_mem )) deallocate( bcw2_mem )
  if(nphase>0) then
     allocate( bcw1_mem(nphase) )
     allocate( bcw2_mem(nphase) )
     call nullify(bcw1_mem)
     call nullify(bcw2_mem)
  else
     allocate( bcw1_mem(1) )
     allocate( bcw2_mem(1) )
     allocate( bcw1_mem(1)%ptr(1) )
     allocate( bcw2_mem(1)%ptr(1) )
  end if

  DO IP=1,NOPHAS
     
     !     Extract boundary conditions
     istate = phase2state_index(ip)
     Velocity => extract_vector_field(state(istate),'Velocity')
     call getbcuvw_sizes(Velocity,nobcu(ip),nobcv(ip),nobcw(ip))

     if(associated( bcu1_mem(IP)%ptr )) deallocate( bcu1_mem(IP)%ptr )
     if(associated( bcu2_mem(IP)%ptr )) deallocate( bcu2_mem(IP)%ptr )
     allocate( bcu1_mem(IP)%ptr(nobcu(ip)) ) 
     allocate( bcu2_mem(IP)%ptr(nobcu(ip)) ) 

     if(associated( bcv1_mem(IP)%ptr )) deallocate( bcv1_mem(IP)%ptr )
     if(associated( bcv2_mem(IP)%ptr )) deallocate( bcv2_mem(IP)%ptr )
     allocate( bcv1_mem(IP)%ptr(nobcv(ip)) ) 
     allocate( bcv2_mem(IP)%ptr(nobcv(ip)) ) 

     if(associated( bcw1_mem(IP)%ptr )) deallocate( bcw1_mem(IP)%ptr )
     if(associated( bcw2_mem(IP)%ptr )) deallocate( bcw2_mem(IP)%ptr )
     allocate( bcw1_mem(IP)%ptr(nobcw(ip)) ) 
     allocate( bcw2_mem(IP)%ptr(nobcw(ip)) ) 
      
     call getbcuvw12(Velocity, &
          bcu1_mem(ip)%ptr, &
          bcv1_mem(ip)%ptr, &
          bcw1_mem(ip)%ptr, &
          bcu2_mem(ip)%ptr, &
          bcv2_mem(ip)%ptr, &
          bcw2_mem(ip)%ptr)

     CALL CHKBC(NOBCU(IP),bcu2_mem(ip)%ptr,NONODS)
     CALL FixBoundaryConditions(halo_tag_p, NONODS, NOBCU(IP),&
          BCU1_mem(IP)%ptr, bcu2_mem(ip)%ptr)

     CALL CHKBC(NOBCV(IP),BCV2_mem(IP)%ptr,NONODS)
     CALL FixBoundaryConditions(halo_tag_p, NONODS, NOBCV(IP),&
        BCV1_mem(IP)%ptr, BCV2_mem(IP)%ptr)

     if(d3) then
        CALL CHKBC(NOBCW(IP),BCW2_mem(IP)%ptr,NONODS) 
        CALL FixBoundaryConditions(halo_tag_p, NONODS, NOBCW(IP),&
             BCW1_mem(IP)%ptr, BCW2_mem(IP)%ptr)
     end if
           
  END DO

  !     Rotation? 
  IF(ROTAT) THEN  
     ! here we assume ROTAT and GETTAN have been determined from 
     ! new options inside COMSCA -Stephan 
     
     ! new framework allows all combinations of ROTAT and GETTAN 
     ! to be different for velocities in different phases and even  
     ! for different b.c.s of velocity in one phase 
     ! legacy fluidity does not allow it, we should check -Stephan 
     
     ! by this assumption we Extracting everything from phase 1  
     Velocity => extract_vector_field(state(1),'Velocity')  
     call getnnodro(Velocity, nnodro) 
  end IF
    
  !     Extract boundary conditions - temperature.

  if(associated( bct1_mem )) deallocate( bct1_mem )
  if(associated( bct2_mem )) deallocate( bct2_mem )
  if(ntsol>0) then
     allocate( bct1_mem(ntsol) )
     allocate( bct2_mem(ntsol) )
     call nullify(bct1_mem)
     call nullify(bct2_mem)
  else
     allocate( bct1_mem(1) )
     allocate( bct2_mem(1) )
     allocate( bct1_mem(1)%ptr(1) )
     allocate( bct2_mem(1)%ptr(1) )
  end if

  nobcts = 0
  DO IT=1,NTSOL
     Field_Buffer => extract_scalar_field( &
          state(field_state_list(it)), &
          trim(field_name_list(it)) )
     call getbct_size(Field_Buffer,nobct(it))
     if(associated( bct1_mem(IT)%ptr )) deallocate( bct1_mem(IT)%ptr )
     if(associated( bct2_mem(IT)%ptr )) deallocate( bct2_mem(IT)%ptr )
     allocate( bct1_mem(IT)%ptr(nobct(it)) ) 
     allocate( bct2_mem(IT)%ptr(nobct(it)) ) 

     call getbct(field_buffer, bct1_mem(IT)%ptr, &
          & bct2_mem(IT)%ptr)

     CALL CHKBC(NOBCT(IT),BCT2_mem(IT)%ptr,node_count(Field_Buffer))
     CALL FixBoundaryConditions(halo_tag_p, node_count(Field_Buffer), NOBCT(IT),&
          BCT1_mem(IT)%ptr, BCT2_mem(IT)%ptr)
     !     Count the no of bc's
     NOBCTS=NOBCTS+NOBCT(IT)
  END DO

  ewrite(1, *) "Exiting redfil"
  
end subroutine redfil

subroutine sndgln2tsndgl(sndgln, tsndgl, stotel_times_snloc, nonods, sufnod)
  integer, intent(in):: stotel_times_snloc, nonods
  integer, intent(in):: sndgln(1:stotel_times_snloc)
  integer, intent(out):: tsndgl(1:stotel_times_snloc)
  integer, intent(out):: sufnod
  
  integer, dimension(:), allocatable:: nod2sufnod
  integer nod, i
  
  allocate(nod2sufnod(1:nonods))
  nod2sufnod=0
  
  ! mark surface nodes with nod2sufnod(nod)==1
  do i=1, stotel_times_snloc
    nod=sndgln(i)
    if (nod2sufnod(nod)==0) then
      nod2sufnod(nod)=1
    end if
  end do
    
  ! create numbering in the same order as full nodal numbering:
  sufnod=0
  do i=1, nonods
    if (nod2sufnod(i)==1) then
      sufnod=sufnod+1
      nod2sufnod(i)=sufnod
    end if
  end do
    
  ! map sndgln to tsndgl in surface node numbering
  do i=1, stotel_times_snloc
    nod=sndgln(i)
    tsndgl(i)=nod2sufnod(nod)
  end do
      
  deallocate(nod2sufnod)
  
end subroutine sndgln2tsndgl


subroutine initialise_state_phase_lists_from_options()

    logical, save:: initialised=.false.
    integer :: nphase, counter, p, nmaterial_phases

    if (initialised) return

    nmaterial_phases = option_count('/material_phase')  
    allocate(state2phase_index(nmaterial_phases))
    state2phase_index = 0

    call get_nphase(nphase)
    allocate(phase2state_index(nphase))

    counter = 0
    do p = 0, nmaterial_phases-1
      if (have_option('/material_phase['//int2str(p)//']/vector_field::Velocity')) then
          ! don't know if prescribed or diagnostic fields should be included in nphase but
          ! suspect that for things like traffic they should be
          ! definitely don't want aliased - crgw
          if (.not.have_option('/material_phase['//int2str(p)//']/vector_field::Velocity/aliased')) then
            counter = counter + 1
            state2phase_index(p+1) = counter
            phase2state_index(counter) = p+1
          end if
      end if
    end do
  
    initialised = .true.

end subroutine initialise_state_phase_lists_from_options

subroutine initialise_state_phase_lists_from_file(nphase)

    logical, save:: initialised=.false.
    integer, intent(in) :: nphase
    integer :: p

    if (initialised) return

    allocate(state2phase_index(nphase))
    allocate(phase2state_index(nphase))

    ! just a one to one relationship I think
    do p = 1, nphase
      phase2state_index(p) = p
      state2phase_index(p) = p
    end do

    initialised = .true.

end subroutine initialise_state_phase_lists_from_file

end module redfil_module
