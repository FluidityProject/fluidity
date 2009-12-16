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

module AdvectionDiffusion
  use climatology
  use detnlxr_module
  use parallel_tools
  use shape_transformations
  use AllSorts
  use FLDebug
  use quadrature
  use elements
  use state_module
  use OceanSurfaceForcing
  use solvers
  use time_dependent_t_source
  use volumesource
  use spud
  use Legacy_Advection_Diffusion_CV, only: ASSEMBLE_FIELD_ADVECTION_CV
  use Field_Equations_CV, only: advection_diffusion_eqn_cv
  use global_parameters, only: halo_tag, halo_tag_p, name_ident, OPTION_PATH_LEN
  use hart3d_allsorts
  use advection_diffusion_cg_3d
  use legacy_boundary_conditions
  use shape_module
  use surface_integration
  use EletemFeinte_module
  use tr2d_module
  use assnav_module
  use solmom_module
  use diff3d_module
  use advection_diffusion_cg_2d
  use assemble_boundary_conditions
  use field_options
  use redfil_module
  use spaerr_module
  use legacy_tensor_fields
  use sparsity_patterns_meshes

  implicit none
  
  private
  
  public :: advdif, get_copied_field, normgi
  
  interface
    subroutine fluxes_settimeseconds(current_time)
      implicit none
      real, intent(in) :: current_time
    end subroutine fluxes_settimeseconds
  end interface

contains

  SUBROUTINE ADVDIF(NEWMES2,&
       NOBCT,BCT1,BCT2,&
       NONODS,TOTELE,r0,d0, &
       NLOC,NGI,&
       DISOPT,NDISOT2,THETA,BETA,LUMP,&
       SNONOD,VNONOD, &
       CGSOLQ,GMRESQ,&
       SUFTEM,SNLOC,SNGI,&
                                ! sub element modelling...
       NSUBTLOC,NSUBNVLOC,&
       VERSIO,ISPHERE,&
       metric_tensor, FREDOP&
       ,state, field_name)

    !     THIS SUB FORMS AND SOLVES A TIME DEPENDENT ADVECTION 
    !     DIFFUSION EQUATION. 
    !     The discretization scheme is controlled through DISOPT(DISCRETIZATION OPTION)
    !     NDISOT2.NE.0 switches on the nonlinear discretisation CV method. 
    !     DISOPT .LT. 0 then optimaly upwinded. 
    !     The following are absolute values of DISOPT(NOT SPHERICAL COORDS)...
    !     DISOPT=1 - Petrof Galerkin x-t.
    !     DISOPT=2 - Petrof Galerkin x-t(BUT weighted to get correct diffusion at steady state).
    !     DISOPT=3 - Petrof Galerkin (x,y).
    !     DISOPT=4 - Least squares (x,y,t). (symmetric)
    !     DISOPT=5 - Least squares (x,y).
    !     DISOPT=6 - Least squares THETA-method. (symmetric)
    !     DISOPT=7 - Galerkin THETA-method.(symmetric when MAKSYM=.true.)
    !     DISOPT=8 - Space weighting. 
    !     For DISOPT=7 ...
    !     If THETA=0.5 then Crank-Nicolson is used in Navier Stokes equ's.
    !     If THETA=1.0 then bakward Euler is used in NS equns. THETA=2/3 Galerkin. 
    !     FOR SPHERICAL COORDS...
    !     DISOPT=1 - balancing diffusion based on (x,y) space.
    !     DISOPT=2 - Laxwendrof balancing diffusion.
    !     DISOPT=3 - (x,y,t) -balancing diffusion. 
    !     DISOPT=4 - No balancing diffusion.
    !     Similarly for the temperature equation.
    !     SCATVE=.true. Then scale velocities according to RSCATV
    !     IF(GETMAT) then get mass current matrix(may change from time step to time step). 
    !     IF MAKSYM then make BIGM symmetrix I.E take advection terms out of matrix.
    !     IF LUMP then lump the mass matrix else dont. 
    !     IF MLCENT then put ML to the diagonal of BIGM matrix. 
    !     BETA controls the conservative discretisation 
    !     BETA=0. is the usual advection form, BETA=1. is divergence form. 
    !     BETA=0.5 is an average of the divergence and advection form.(recommended).
    !     ABSORB =absorbtion of energy etc. 
    !     MUPTXX,MUPTYY =components of diffusivities. 
    !     NU,NV,NW are for the non-linear terms ordeneraly NU=U,NV=V,NW=W.
    !     UG,VG,WG are the grid velocities. 
    !     NDGLNO=element pter for unknowns. 
    !     SONDGL=element pter for materials and sources. 
    !     VONDGL=element pter for advection NU,UG.
    !     XONDGL=element pter for coordinates(x,y,z)
    !     If INMAT then return matricie(s).  
    !     NB We can solve for stream function with clever manipulation 
    !     of the quantities. 
    !     NB SOURCE,DENPT, & diffusivities are of length SNONOD
    !     and use same pointer SONDGL.
    real, intent(in) :: r0
    INTEGER NCOLM,NLOC,NGI
    REAL TOLER
    PARAMETER(TOLER=0.0001)
    INTEGER DISOPT
    integer NDISOT2
    INTEGER TOTELE,NONODS,SNONOD,VNONOD
    INTEGER SNLOC,SNGI
    INTEGER UZAWA,CGSOLQ,GMRESQ
    !
    REAL THETA,BETA
    !     If PHI=0.5 then Crank-Nicolson is used in Navier Stokes equ's.
    !     If PHI=1.0 then bakward Euler is used in NS equns.
    !     Similarly for the temperature equation except width variable THETA.
    !     have gravity in -ve y direction.

    LOGICAL GETC12,GETMAT,MLCENT
    LOGICAL LUMP,SYM,MASSON,VLKON

    !     FORMULATE BOUNDARY CONDITIONS
    !     This subroutine adds in the DIRICHLET boundary conditions. 
    !     IDIM controls how many values to set. 
    !     BCU2,BCV2,BCW2 contain the nodes at which the boundary 
    !     conditions are to be applied. 
    !     BCU1 caontains  the value of the specified acceleration of U 
    !     that is BCU1=(Unew - Uold)/DT
    INTEGER NOBCT
    INTEGER BCT2(NOBCT)
    REAL BCT1(NOBCT)
    !     This sub solves an advection diffusion equation for T 
    !     For a time dependent problem.  
    !     If DEFAULT then use the default settings to solve the equation. 
    !     When DEFAULT some of the settings can change. 
    !     T is from the previous time step. 
    !     RUB is a working array.
    !     For parallel processing ...
    !     NNODP=no of subdomain nodes excluding halo nodes. 
    LOGICAL D2
    !     These are pointer for element working space...
    LOGICAL SUFTEM
    !     SNDGLN - points to nodes (global) from surface elements. 
    !     TSNDGL - points to nodes (global) from surface elements(SURFACE ELEMENT NODES). 
    !     SN,SNLX,SNLY,SWEIGH -surface element shape functions. 
    !     STOTEL,SNLOC,SNGI - no of surface elements, nodes per, Gauss pts per.
    real, intent(in) :: d0
    INTEGER NSUBTLOC,NSUBNVLOC
    LOGICAL IGUESS
    LOGICAL PUTMAT
    INTEGER NDISOT
    real, dimension(:), allocatable :: TOLD,TNEW
    ! scalar fields to wrap TOLD and TNEW
    type(scalar_field) :: l_told_field, l_tnew_field
    INTEGER NITS,NONITS,DISOPN
    INTEGER MXNODS,IDIM,ele
    INTEGER IDUM(1)
    INTEGER NOD1,NOD2,NOD3,NOD4
    REAL VOLELE
    REAL RMAXT,RMINT,RMAXNU,RMINNU,RMAXNV,RMINNV
    real, dimension(:), allocatable::VECXSTOR
    real, dimension(:), allocatable::BIGMSTOR
    real, dimension(:), allocatable::MATTSTOR
    real, dimension(:), allocatable::tempML
    real, dimension(:), allocatable::tempSOURCE
    real, dimension(:), allocatable::tempABSORB
    real, dimension(:), allocatable::dm1 ! Something to do with boundaries.

    REAL, ALLOCATABLE, DIMENSION(:,:,:)::DMATINVCSTORE
    REAL, ALLOCATABLE, DIMENSION(:,:)::DINVSOUSUBTSTORE

    INTEGER ISUB
    INTEGER VERSIO,ISPHERE
    LOGICAL NEWMES2
    !     Encapsulated state information containing fields. -dham
    ! state now has to be passed in as an array for compatibility with subroutines
    ! that are used elsewhere in the code however it is only of length 1
    type(state_type), dimension(1), intent(inout) :: state
    !     The name of the field were advecting/diffusing
    !     Used to retrieve the field from state
    character(len=*), intent(in):: field_name
    real :: source_factor = 1.0

    ! for LES -Jemma
    type(vector_field), pointer :: Coordinate
    type(scalar_field) :: X_coord, Y_coord, Z_coord
    type(tensor_field) :: metric_tensor
    REAL, dimension(:), allocatable:: TL2MXX,TL2MXY,TL2MXZ
    REAL, dimension(:), allocatable:: TL2MYY,TL2MYZ,TL2MZZ
    integer :: NDIM, FREDOP

    type(vector_field) :: U_nl, U_nl_backup, Gravity
    type(scalar_field) :: Sink

    ! Diffusivity tensor
    type(tensor_field) :: diffusivity
    real, dimension(:), allocatable :: MUPTXX,MUPTXY,MUPTXZ,&
         MUPTYY,MUPTYZ,MUPTZZ

    ! Local memory
    type(scalar_field), pointer:: tfield
    character(len=OPTION_PATH_LEN) option_path
    REAL, ALLOCATABLE, DIMENSION(:)::VECT_SGSADD
           INTEGER NGI2
           REAL, ALLOCATABLE, DIMENSION(:,:)::N2
           REAL, ALLOCATABLE, DIMENSION(:,:)::NLX2
           REAL, ALLOCATABLE, DIMENSION(:,:)::NLY2
           REAL, ALLOCATABLE, DIMENSION(:,:)::NLZ2
           REAL, ALLOCATABLE, DIMENSION(:)::WEIGHT2
           REAL, ALLOCATABLE, DIMENSION(:)::L1,L2,L3,L4

    ! Seem to be unused, but are passed to solmom
    real, dimension(:), allocatable :: f, savdia
    
    ! Matrix, right hand side vector and mass-lump vector.
    real, dimension(:), pointer ::  VECX,BIGM,ML
    ! The matrix to be solved for.
    type(csr_matrix) :: M
    ! a matrix type to wrap MATRIX within advdif
    type(csr_matrix) :: A_m
    ! sparsity for matrix types
    type(csr_sparsity), pointer :: mesh_sparsity
    ! wrapper around vecx
    type(scalar_field) :: rhs, masslump

    character(len=OPTION_PATH_LEN) :: velocity_option_path
   
    ! Used for data unpacking
    integer :: stat
    ! Array of zeros so that unpacked state variables always point at something
    real, dimension(:), allocatable, target :: zero_real_array
    type(csr_sparsity), pointer :: sparsity
    type(mesh_type), pointer :: mesh
    type(scalar_field), pointer :: s_field
    type(vector_field), pointer :: v_field

    ! Unpacked state data
    integer :: xnonod, stotel, sufnod
    integer, dimension(:), pointer :: ndglno, sondgl, vondgl, xondgl
    integer, dimension(:), allocatable :: sndgln, tsndgl
    integer, dimension(:), pointer :: centrm, colm, findrm
    real, dimension(:), pointer :: x, y, z
    real, dimension(:), pointer :: xold, yold, zold
    real, dimension(:), pointer :: ug, vg, wg
    real, dimension(:), pointer :: nu, nv, nw
    real, dimension(:), pointer :: t, ttnew
    real, dimension(:), pointer :: tsub, tnewsub
    real, dimension(:), pointer :: nusub, nvsub, nwsub
    real, dimension(:), pointer :: absorb, source
    real, dimension(:), pointer :: weight, sweigh
    real, dimension(:, :), allocatable :: n, nlx, nly, nlz
    real, dimension(:, :), allocatable :: sn, snlx, snly, snlz
    real, dimension(:), allocatable :: salphe, taire
    real, dimension(:), allocatable :: denpt
    
    ! Unpacked options tree data
    integer :: ident
    logical :: d3, mvmesh
    real :: acctim, dt

    ! Unpacked parallel and halo data
    integer :: nnodp, para
    
    ! Temporary parameters
    logical, parameter :: dsph = .false.
    logical, parameter :: dcyl = .false.
    logical, parameter :: inmat = .true.
    logical, parameter :: maksym = .false.
    
    ewrite(1, *) "In advdif"
    ewrite(2, *) "Solving for field " // trim(field_name) // " in state " // trim(state(1)%name)
    
    ewrite(2, *) "Unpacking state data"  
    
    v_field => extract_vector_field(state(1), "Coordinate")
    xnonod = node_count(v_field)
    
    allocate(zero_real_array(max(snonod, totele * nsubtloc, totele * nsubnvloc, vnonod, xnonod)))
    zero_real_array = 0.0
 
    v_field => extract_vector_field(state(1), "Coordinate")
    x => zero_real_array(1:xnonod)
    y => zero_real_array(1:xnonod)
    z => zero_real_array(1:xnonod)
    if(mesh_dim(v_field) >= 1) x => v_field%val(1)%ptr
    if(mesh_dim(v_field) >= 2) y => v_field%val(2)%ptr
    if(mesh_dim(v_field) >= 3) z => v_field%val(3)%ptr
    
    mesh => extract_mesh(state(1), "CoordinateMesh")
    xondgl => mesh%ndglno
    stotel = surface_element_count(v_field)
    allocate(sndgln(stotel * snloc))
    call getsndgln(mesh, sndgln)
    sufnod = stotel * snloc
    allocate(tsndgl(sufnod))
    allocate(salphe(sufnod))
    allocate(taire(sufnod))
    if(option_count("/material_phase/scalar_field/prognostic/discretisation/legacy_discretisation/legacy_suftem") > 0) then
      call sndgln2tsndgl(sndgln, tsndgl, stotel * snloc, nonods, sufnod)
      s_field => extract_scalar_field(state(1), field_name)
      call getsalphetaire(s_field, tsndgl, salphe, taire)
    else
      tsndgl = 0
      salphe = 0
      taire = 0
    end if
    d3 = (mesh_dim(mesh) == 3)
    
    v_field => extract_vector_field(state(1), "OldCoordinate", stat)
    xold => zero_real_array(1:xnonod)
    yold => zero_real_array(1:xnonod)
    zold => zero_real_array(1:xnonod)
    if(stat == 0) then
      if(mesh_dim(v_field) >= 1) xold => v_field%val(1)%ptr
      if(mesh_dim(v_field) >= 2) yold => v_field%val(2)%ptr
      if(mesh_dim(v_field) >= 3) zold => v_field%val(3)%ptr
    end if
    
    v_field => extract_vector_field(state(1), "Velocity")
    
    mesh => extract_mesh(state(1), "VelocityMesh")
    ndglno => mesh%ndglno
    sondgl => ndglno
    vondgl => ndglno
    
    allocate(n(ele_loc(mesh, 1), ele_ngi(mesh, 1)))
    allocate(nlx(ele_loc(mesh, 1), ele_ngi(mesh, 1)))
    allocate(nly(ele_loc(mesh, 1), ele_ngi(mesh, 1)))
    allocate(nlz(ele_loc(mesh, 1), ele_ngi(mesh, 1)))
    call extract_old_element(ele_shape(mesh, 1), n, nlx, nly, nlz)
    allocate(sn(ele_loc(mesh, 1), face_ngi(mesh, 1)))
    allocate(snlx(ele_loc(mesh, 1), face_ngi(mesh, 1)))
    allocate(snly(ele_loc(mesh, 1), face_ngi(mesh, 1)))
    allocate(snlz(ele_loc(mesh, 1), face_ngi(mesh, 1)))
    call extract_old_element(face_shape(mesh, 1), sn, snlx, snly)

    allocate(DM1(NONODS))
    
    v_field => extract_vector_field(state(1), "GridVelocity", stat)
    ug => zero_real_array(1:vnonod)
    vg => zero_real_array(1:vnonod)
    wg => zero_real_array(1:vnonod)
    if(stat == 0) then
      if(mesh_dim(v_field) >= 1) ug => v_field%val(1)%ptr
      if(mesh_dim(v_field) >= 2) vg => v_field%val(2)%ptr
      if(mesh_dim(v_field) >= 3) wg => v_field%val(3)%ptr
    end if
    
    U_nl=extract_vector_field(state(1), "NonlinearVelocity")
    Sink=extract_scalar_field(state(1), trim(field_name)//"SinkingVelocity"&
         &, stat=stat)
    if (stat==0) then
       Gravity=extract_vector_field(state(1), "GravityDirection")

       U_nl_backup=U_nl

       call allocate(U_nl, U_nl%dim, U_nl%mesh, "NonlinearVelocity")
       
       call set(U_nl, U_nl_backup)
       call addto(U_nl, Gravity, scale=Sink)

    else
       ! Grab an extra reference to cause the deallocate below to be safe.
       call incref(U_nl)
    end if
    nu => zero_real_array(1:vnonod)
    nv => zero_real_array(1:vnonod)
    nw => zero_real_array(1:vnonod)
    if(mesh_dim(U_nl) >= 1) nu => U_nl%val(1)%ptr
    if(mesh_dim(U_nl) >= 2) nv => U_nl%val(2)%ptr
    if(mesh_dim(U_nl) >= 3) nw => U_nl%val(3)%ptr
    
    v_field => extract_vector_field(state(1), "VelocityInnerElement", stat)
    nusub => zero_real_array(1:totele * nsubnvloc)
    nvsub => zero_real_array(1:totele * nsubnvloc)
    nwsub => zero_real_array(1:totele * nsubnvloc)
    if(stat == 0) then
      if(mesh_dim(v_field) >= 1) nusub => v_field%val(1)%ptr
      if(mesh_dim(v_field) >= 2) nvsub => v_field%val(2)%ptr
      if(mesh_dim(v_field) >= 3) nwsub => v_field%val(3)%ptr
    end if

    mesh => extract_mesh(state(1), "PressureMesh")
    weight => mesh%shape%quadrature%weight
    sweigh => mesh%faces%shape%quadrature%weight
    
    s_field => extract_scalar_field(state(1), field_name)      
    if(node_count(s_field) == nonods) then
      t => s_field%val
    else
      FLAbort("Called advdif for a non-backwards-compatible tracer field")
    end if
    
    s_field => extract_scalar_field(state(1), "Iterated" // field_name)
    ttnew => s_field%val

    s_field => extract_scalar_field(state(1), trim(field_name) // "InnerElement", stat)
    tsub => zero_real_array(1:totele * nsubtloc)
    if(stat == 0) then
      tsub => s_field%val
    end if

    s_field => extract_scalar_field(state(1), "Iterated" // trim(field_name) // "InnerElement", stat)
    tnewsub => zero_real_array(1:totele * nsubtloc)
    if(stat == 0) then
      tnewsub => s_field%val
    end if

    s_field => extract_scalar_field(state(1), trim(field_name) // "Absorption", stat)
    absorb => zero_real_array(1:snonod)
    if(stat == 0) then
      if(.not. aliased(s_field) .and. node_count(s_field) == nonods) then
        absorb => s_field%val
      end if
    end if

    s_field => extract_scalar_field(state(1), trim(field_name) // "Source", stat)
    source => zero_real_array(1:snonod)
    if(stat == 0) then
      if(.not. aliased(s_field) .and. node_count(s_field) == nonods) then
        source => s_field%val
      end if
    end if
    
    ! Diffusivity
    diffusivity= extract_tensor_field(state(1), &
         trim(field_name) // "Diffusivity", stat)

    allocate(MUPTXX(SNONOD),MUPTXY(SNONOD),MUPTXZ(SNONOD), &
         MUPTYY(SNONOD),MUPTYZ(SNONOD),MUPTZZ(SNONOD))
    if(stat==0) then
       call copy_tensor_to_legacy(diffusivity, MUPTXX, MUPTXY, MUPTXZ, &
              MUPTYY, MUPTYZ, MUPTZZ)
    else
       MUPTXX=0.0
       MUPTXY=0.0
       MUPTXZ=0.0
       MUPTYY=0.0
       MUPTYZ=0.0
       MUPTZZ=0.0
    end if
        
    mesh => extract_mesh(state(1), "VelocityMesh")
    sparsity => get_csr_sparsity_firstorder(state, mesh, mesh)
    centrm => sparsity%centrm
    colm => sparsity%colm
    findrm => sparsity%findrm
    
    ncolm=size(colm)

    ! This looks unused
    allocate(denpt(nonods))
    denpt = 1.0
    
    ident = name_ident(field_name)
    
    sparsity => null()
    mesh => null()
    s_field => null()  
    v_field => null()
    
    ewrite(2, *) "Unpacking options tree data"
    
    mvmesh = have_option("/mesh_adaptivity/mesh_movement")
    call get_option("/timestepping/current_time", acctim)
    call get_option("/timestepping/timestep", dt)
    
    ewrite(2, *) "Unpacking parallel data"
    
    if(isparallel()) then
      para = 1
      nnodp = get_nowned_nodes(halo_tag)
    else
      para = 0
      nnodp = nonods
    end if
    
       ! won't work with multi-phase/material
     velocity_option_path='/material_phase[0]/vector_field::Velocity'

    !     ---------------------------------------------- End solidity addition

    !     NB LEASQR can only be used for D3 and D2. 
    !     If ADV then assume pure advection. 
    !     If DIFF then assume pure diffusion. 
    !     NB HART3D can assble 3 matrix equations from 3 equations 
    !     that are the same except for the source terms.( not used

    ewrite(2, *) 'r2norm(T,nnodp,0):',r2norm(T,nnodp,0)
    option_path=""      
    nullify(tfield)

    tfield => extract_scalar_field(state(1), trim(field_name))
    !     options tree path to field for solver options (via solmom)
    option_path=tfield%option_path

#if defined(HAVE_NETCDF)&&defined(HAVE_LIBUDUNITS)&&defined(DOUBLEP)
    if(have_option("/environmental_data/ERA40")) then
       call fluxes_settimeseconds(acctim)
    end if
    if(have_option("/environmental_data/climatology")) then
       call climatology_SetTimeSeconds(acctim)
    end if
#endif


    ewrite(2, *) 'D3:',D3
    ewrite(2, *) 'r2norm(source,nnodp,0):',r2norm(source,nnodp,0)
    ewrite(2, *) 'r2norm(T,nnodp,0):',r2norm(T,nnodp,0)
    ewrite(2, *) 'NOBCT=',NOBCT

    SYM=.FALSE.
    GETC12=.FALSE.
    GETMAT=.FALSE.
    MLCENT=.FALSE.
    MASSON=.FALSE. 
    VLKON=.FALSE.

    D2=.TRUE.
    IF(D3) D2=.FALSE.
    NDISOT=MOD(NDISOT2,1000)
    ! IF ABS(NDISOT2).GT.1000 then do not repeatedly form element 
    ! matrices used for limiting. 

    UZAWA=0

    !     Get the pointers for the element work space
    MXNODS=MAX(NONODS,SNONOD,VNONOD)

    !     Non-linear iteration loop (NDISOT=non-linear discretisation)...
    NONITS=INT(ABS(NDISOT)/10)
    DISOPN=MOD(INT(ABS(NDISOT)/1),10)

    if(have_option(trim(option_path)//"/prognostic/spatial_discretisation/legacy_mixed_cv_cg")) then
      DENPT = 1.0  ! density terms are dealt with later from state now
    end if
    mesh=> extract_mesh(state(1), "VelocityMesh")
    mesh_sparsity=>get_csr_sparsity_firstorder(state, mesh, mesh)
    mesh => null()
    call allocate(M, mesh_sparsity, name="FieldMatrix")
    call allocate(A_m, mesh_sparsity, name="FieldAdvectionMatrix")
    call allocate(rhs, tfield%mesh, name="RHSFieldEqn")
    call allocate(masslump, tfield%mesh, name="FieldMassLump")
    bigm=>M%val
    vecx=>rhs%val
    ML=>masslump%val

    ewrite(3,*) 'NONITS,DISOPN,NDISOT:',&
         NONITS,DISOPN,NDISOT

    IF(ABS(NDISOT).GE.10) THEN
       allocate(told(nonods), tnew(nonods))
       TOLD = T
       TNEW = T
       ! wrap the iterative fields for back compatibility
       l_told_field=wrap_scalar_field(tfield%mesh, val=TOLD, name="NonlinearAdvectionTOld")
       l_tnew_field=wrap_scalar_field(tfield%mesh, val=TNEW, name="NonlinearAdvectionTNew")
    ENDIF
    

    IF(D2) THEN
       ewrite(1,*)'About to disappear into hart2d'
       CALL HART2D(T,NU,NV,UG,VG,&
            SOURCE,X,Y,&
            VECX,&
            BIGM,BIGM, BIGM, ML,&
            FINDRM,COLM,NCOLM,NONODS,TOTELE,&
            FINDRM,COLM,NCOLM,NONODS,CENTRM, &
            N,N,NLX,NLY, WEIGHT, NLOC,NGI,NLOC, &
            DISOPT,DISOPN,ABS(NDISOT),DT,THETA,BETA,LUMP,MAKSYM,&
            .FALSE.,.FALSE.,&
            GETC12,INMAT,MLCENT, &
            NDGLNO,NDGLNO,   &
            SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
            DENPT,MUPTXX,MUPTYY, ABSORB)

       ewrite(1,*)'have left hart2d'
    ENDIF

    IF(D3) THEN

       ewrite(1,*) 'in advdif, if(d3)'
       ewrite(1,*) 'about to call one of the versions of hart3d'
       ewrite(2,*) 'NONODS=',NONODS

       IF(NSUBTLOC.NE.0) THEN


          ISUB=1

            
          IF(NGI.EQ.1) THEN
            NGI2=4
          ELSE
            NGI2=NGI
          ENDIF

          ALLOCATE(N2(NLOC,NGI2))
          ALLOCATE(NLX2(NLOC,NGI2))
          ALLOCATE(NLY2(NLOC,NGI2))
          ALLOCATE(NLZ2(NLOC,NGI2))
          ALLOCATE(WEIGHT2(NGI2))
          ALLOCATE(L1(NGI2))
          ALLOCATE(L2(NGI2))
          ALLOCATE(L3(NGI2))
          ALLOCATE(L4(NGI2))
          CALL TRIQUA(L1, L2, L3, L4, weight2, .true.,NGI2)
          CALL SHATRI(L1, L2, L3, L4, weight2, .true., &
     &            nLOC,NGI2, &
     &            n2,nLX2,nLY2,nLZ2)
            

          IF(BETA.GT.5.0) THEN

             ! Calculate ML...
             ALLOCATE(tempML(NONODS))
             ALLOCATE(tempSOURCE(NONODS))
             ALLOCATE(tempABSORB(NONODS))

             ML=0.0
             do ELE=1,TOTELE
                NOD1=NDGLNO((ELE-1)*NLOC+1)
                NOD2=NDGLNO((ELE-1)*NLOC+2)
                NOD3=NDGLNO((ELE-1)*NLOC+3)
                NOD4=NDGLNO((ELE-1)*NLOC+4)
                VOLELE=ABS(TetVolume( &
                     x(NOD1), y(NOD1), z(NOD1),         &
                     x(NOD2), y(NOD2), z(NOD2),         &
                     x(NOD3), y(NOD3), z(NOD3),         &
                     x(NOD4), y(NOD4), z(NOD4)))
                ML(NOD1)=ML(NOD1)+0.25*VOLELE
                ML(NOD2)=ML(NOD2)+0.25*VOLELE
                ML(NOD3)=ML(NOD3)+0.25*VOLELE
                ML(NOD4)=ML(NOD4)+0.25*VOLELE
             END DO
             VECX=0.0
             tempML=0.0
             ! *************ocean forcing/relaxation
             if(have_option("/environmental_data/climatology")) then
                ! ml is not changed here but is used (SEND DOWN tempML when changed)...
                if(ident.eq.-1) then
                   call RelaxTemperatureToClimatology(NONODS, X, Y, Z, &
                        T, ML, tempML, dt, centrm, BIGM, 0, VECX)
                else if(ident.eq.42) then
                   call RelaxSalinityToClimatology(NONODS, X, Y, Z, T, &
                        ML, tempML, dt, centrm, BIGM, 0, VECX)
                end if
             end if

             if(have_option("/environmental/ERA40").and.(ident.eq.-1)) then
                CALL SurfaceHeatFlux(VECX,X,Y,Z,&
                     NONODS,SNLOC,SNGI,&
                     SN,SNLX,SNLY,sWEIGH)
             else if(have_option("/environmental/ERA40").and.(ident.eq.42)) then
                call SurfaceSalinityFlux(VECX,X,Y,Z,T,&
                     NONODS,SNLOC,SNGI,&
                     SN,SNLX,SNLY,SWEIGH)
             endif

             if(have_option("/environmental_data/climatology")) then
                if(ident.eq.-1) then
                   call SurfaceTemperatureRelaxation(X,Y,Z,&
                        snloc, sngi, sn, snlx, snly, sweigh,&
                        T,tempML,dt,centrm,BIGM,0,VECX)
                else if(ident.eq.42) then
                   call SurfaceSalinityRelaxation(X,Y,Z,&
                        snloc, sngi, sn, snlx, snly, sweigh,&
                        T,tempML,dt,centrm,BIGM,0,VECX)
                end if
             end if
             ! *************ocean forcing/relaxation
             tempSOURCE=SOURCE+VECX/ML
             tempABSORB=ABSORB+tempML/(ML*DT)


             ALLOCATE(VECT_SGSADD(TOTELE*NLOC))
             VECT_SGSADD=0.0
             ewrite(1, *) 'entering SOLV_FIELD_SGSDG'
             CALL SOLV_FIELD_SGSDG(T,TTNEW,&
                  NU,NV,NW,UG,VG,WG,&
                  tempSOURCE,X,Y,Z,D0,&
                  VECT_SGSADD,0,&
                  VECT_SGSADD,&
                  BIGM, ML,&
                  TOTELE, &
                  FINDRM,COLM,NCOLM,NONODS,CENTRM, &
                  N2,NLX2,NLY2,NLZ2, WEIGHT2, NLOC,NGI2,&
                  DISOPT,DISOPN,DT,THETA,(BETA-10.0),LUMP,0,&
                  .FALSE.,LUMP,&
                  NDGLNO,&
                  SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
                  DENPT,&
                  MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
                  tempABSORB,VERSIO,ISPHERE,&
                  NNODP,&
                                ! for the SGS...
                  NSUBNVLOC, &
                  NUSUB,NVSUB,NWSUB,1,&
                  TSUB,TNEWSUB,ISUB,&
                  NCOLM,0,IDUM,&
                  .FALSE.,&
                  NOBCT,BCT1,BCT2,&
                                ! GMRES solver...
                  PARA,halo_tag,&
                                ! form DGCMC
                  D3, velocity_option_path) 
            goto 42 
          ENDIF


          ALLOCATE(DMATINVCSTORE(TOTELE,NSUBTLOC,NLOC))
          ALLOCATE(DINVSOUSUBTSTORE(TOTELE,NSUBTLOC))

          ewrite(1,*) 'calling hart3dsgs'
          CALL HART3DSGS(T,&
               NU,NV,NW,UG,VG,WG,&
               SOURCE,X,Y,Z,&
               VECX, &
               BIGM,ML,&
               TOTELE,&
               FINDRM,COLM,NCOLM,NONODS, &
               N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI, &
               DISOPT,DISOPN,DT,THETA,BETA,LUMP,&
               (ABS(NDISOT).GE.10),(ABS(NDISOT).GE.10),&
               NDGLNO, &
               SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
               MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ, &
               ABSORB,ISPHERE,&
               nnodp,para,halo_tag,&
               ! for the SGS...
               NSUBTLOC,NSUBNVLOC,&
               NUSUB,NVSUB,NWSUB,&
               TSUB, &
               DMATINVCSTORE,&
               DINVSOUSUBTSTORE)
       ELSE IF((ISPHERE.EQ.2).OR.((DISOPT.GE.146).AND.(DISOPT.LE.154)))&
            THEN

          ewrite(1,*) 'calling hart3dsim'
          ewrite(2,*) 'isphere, disopt', isphere, disopt
          CALL HART3Dsim(T,&
               NU,NV,NW,UG,VG,WG,&
               SOURCE,X,Y,Z,&
               VECX, &
               BIGM, ML,&
               TOTELE,&
               FINDRM,COLM,NCOLM,NONODS,centrm,&
               N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI, &
               DISOPT,DISOPN,DT,THETA,LUMP,&
               (NDISOT.NE.0),(NDISOT.NE.0),&
               NDGLNO, &
               SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
               MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ, &
               ABSORB, ISPHERE,&
               nnodp,para,halo_tag)

       ELSE

          call get_time_dependent_temperature_source()
          source_factor=get_time_dependent_source_factor(ACCTIM)
          ewrite(2,*) 'source_factor', source_factor

          call Volumesource_initialise()

          if(have_volumesource) then
             call set_volumesource(X(1:xnonod),Y(1:xnonod)&
                  ,Z(1:xnonod),source(1:nonods),state(1))
          end if

          ! IF LES THEN PUT INVERSE OF METRIC INTO MUPTXX,MUPTXY ETC
          !
          ! If DISOPT=42, constant length scale
          ! If DISOPT=43, isotropic length scale
          ! If DISOPT=44, inverted metric used
          ! If DISOPT=45 or 46, modified inverted metric using rotation tensors

          allocate(tl2mxx(snonod), tl2mxy(snonod), tl2mxz(snonod))
          allocate(tl2myy(snonod), tl2myz(snonod), tl2mzz(snonod))
          tl2mxx=0.0
          tl2mxy=0.0
          tl2mxz=0.0
          tl2myy=0.0
          tl2myz=0.0
          tl2mzz=0.0

          if(disopt>=42 .and. disopt<=50) then

             Coordinate => extract_vector_field(state(1), "Coordinate")
             X_coord=extract_scalar_field(Coordinate, 1)
             Y_coord=extract_scalar_field(Coordinate, 2)
             Z_coord=extract_scalar_field(Coordinate, 3)

             if(D3) then
                NDIM=3
             else
                NDIM=2
             end if
             
             call eletens(X_coord%val, Y_coord%val, Z_coord%val, &
                  tl2mxx, tl2mxy, tl2mxz, tl2myy, tl2myz, tl2mzz,&
                  fredop, snonod, xnonod, totele, nloc, disopt, ndglno)

             IF(DISOPT.EQ.42) THEN
                CALL LESIMEC(NONODS,NDIM,&
                     & TL2MXX,TL2MXY,TL2MXZ,&
                     & TL2MYY,TL2MYZ,TL2MZZ)
             ENDIF
             IF(DISOPT.EQ.43) THEN
                CALL LESIME0(NONODS,NDIM, metric_tensor%val,&
                     & TL2MXX,TL2MXY,TL2MXZ,&
                     & TL2MYY,TL2MYZ,TL2MZZ)
             ENDIF
             IF(DISOPT.EQ.44) THEN
                CALL LESIME(NONODS,NDIM, metric_tensor%val,&
                     &TL2MXX,TL2MXY,TL2MXZ,&
                     &TL2MYY,TL2MYZ,TL2MZZ)
             ENDIF
             IF((DISOPT.EQ.45).OR.(DISOPT.EQ.46)) THEN
                CALL LESIME2(NONODS,NDIM, metric_tensor%val,&
                     & TL2MXX,TL2MXY,TL2MXZ,&
                     & TL2MYY,TL2MYZ,TL2MZZ)
             ENDIF
          end if

          ewrite(1,*) 'calling hart3d'
          CALL HART3D(T,NU,NV,NW,UG,VG,WG,&
               source_factor*SOURCE,&
               X,Y,Z,D0,&
               VECX, &
               BIGM,BIGM,BIGM, BIGM, ML,&
               FINDRM,COLM,NCOLM,NONODS,TOTELE,&
               FINDRM,COLM,NCOLM,NONODS,CENTRM, &
               N,N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,NLOC, &
               DISOPT,ABS(NDISOT),DT,THETA,BETA,LUMP,MAKSYM,&
               (NDISOT.NE.0),(NDISOT.NE.0),&
               GETC12,INMAT,MLCENT, &
               NDGLNO,NDGLNO,   &
               SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
               DENPT,&
               MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ, &
               TL2MXX,TL2MXY,TL2MXZ,TL2MYY,TL2MYZ,TL2MZZ, &
               ABSORB, &
               VERSIO, state(1), tfield)
       ENDIF
       CALL PMINMX(vecx,NONODS,'******vecx  ')
    ENDIF

    ewrite(2,*)'dt,totele,nonods:',dt,totele,nonods

    ewrite(2,*) 'R2NORM OF VECX 2:',R2NORM(VECX,NNODP,PARA)

    !     apply the boundary conditions
    ewrite(2,*) 'SUFTEM=',SUFTEM
    ewrite(2,*) 'NONODS,SUFNOD,STOTEL,SNLOC,SNGI:',&
         NONODS,SUFNOD,STOTEL,SNLOC,SNGI

    IF(SUFTEM) THEN
       PUTMAT=.TRUE.
       CALL TEMSUF( .TRUE., PUTMAT, DT,THETA, &
            VECX,T,NONODS,SUFNOD,STOTEL,SNLOC,SNGI, &
            D3,DCYL,DSPH, &
            SNDGLN,TSNDGL, SALPHE,TAIRE, &
            SN,SNLX,SNLY,SWEIGH, &
            BIGM,NCOLM,COLM,FINDRM,  CENTRM, &
            LUMP,&
            X,Y,Z,R0,D0, &
            NU,NV,NW)
    ENDIF

    IF(MVMESH) THEN
       !     Solidity NOTE ------------------------------------ CRGW 30/01/07
       !     Section MAY NEED TO BE MODIFIED to avoid conflict with Lagrangian solids
       !     Amend the r.h.s vector for a matrix eqn with a moving mesh. 
       CALL AMENDMVMESH(1,VECX,VECX,VECX,&
            T,T,T,DT,&
            X,Y,Z, XOLD,YOLD,ZOLD, &
            N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,&
            TOTELE,D3,DCYL,  &
            XONDGL,NDGLNO, NONODS,XNONOD)
    ENDIF

    if(have_option("/environmental_data/climatology")) then
       if(ident.eq.-1) then
          call RelaxTemperatureToClimatology(NONODS, X, Y, Z, &
               T, ML, ML, dt, centrm, BIGM, ncolm, VECX)
       else if(ident.eq.42) then
          call RelaxSalinityToClimatology(NONODS, X, Y, Z, T, &
               ML, ML, dt, centrm, BIGM, ncolm, VECX)
       end if
    end if

    if(have_option("/environmental/ERA40").and.(ident.eq.-1)) then
       CALL SurfaceHeatFlux(VECX,X,Y,Z,&
            NONODS,SNLOC,SNGI,&
            SN,SNLX,SNLY,sWEIGH)
    else if(have_option("/environmental_data/ERA40").and.(ident.eq.42)) then
       call SurfaceSalinityFlux(VECX,X,Y,Z,T,&
            NONODS,SNLOC,SNGI,&
            SN,SNLX,SNLY,SWEIGH)
    endif

    if(have_option("/environmental_data/climatology")) then
       if(ident.eq.-1) then
          call SurfaceTemperatureRelaxation(X,Y,Z,&
               snloc, sngi, sn, snlx, snly, sweigh,&
               T,ML,dt,centrm,BIGM,ncolm,VECX)
       else if(ident.eq.42) then
          call SurfaceSalinityRelaxation(X,Y,Z,&
               snloc, sngi, sn, snlx, snly, sweigh,&
               T,ML,dt,centrm,BIGM,ncolm,VECX)
       end if
    end if

    DO  NITS=1,MAX(1,NONITS)! Was loop 135
       !
       IF(ABS(NDISOT).GE.10) THEN
          IF(NITS.EQ.1) THEN
             TNEW = TTNEW
          ELSE
             TNEW = T
          ENDIF
          T = TOLD
       ENDIF

       !     Add in the 1st order spatial derivative -transport. ******
       IF(ABS(NDISOT).GE.10) THEN
          IF(NITS.EQ.1) THEN
             ALLOCATE(VECXSTOR(NONODS))
             ALLOCATE(MATTSTOR(NCOLM))
             ALLOCATE(BIGMSTOR(NCOLM))
             VECXSTOR(1:NONODS)=VECX(1:NONODS)
             BIGMSTOR(1:NCOLM)=BIGM(1:NCOLM)
             getmat=.true.
          ELSE
             VECX(1:NONODS)=VECXSTOR(1:NONODS)
             BIGM(1:NCOLM)=BIGMSTOR(1:NCOLM)
             getmat=.false.
          ENDIF
          !      Solidity addition ------------------------------- CRGW 30/01/07
          !          --------------------------------------------------commented out 22/03/07 crgw

          if(have_option(trim(option_path)//"/prognostic/spatial_discretisation&
                         &/legacy_mixed_cv_cg")) then
            ! if using new options then go into the new subroutine
            ! from here this still depends on hart3/2d to have created
            ! the lumped mass matrix (and possibly any diffusional components)
            call advection_diffusion_eqn_cv(M, A_m, rhs, &
                                       l_tnew_field, state, &
                                       getmat, dt, &
                                       reference_field=tfield)
          else
              CALL ASSEMBLE_FIELD_ADVECTION_CV(&
                    NONODS,VNONOD,XNONOD,TOTELE,NLOC,NGI,&
                    SNLOC, STOTEL,&
                    NDGLNO,XONDGL,VONDGL,&
                    SNDGLN,&
                    FINDRM,COLM,NCOLM,CENTRM, &
                    BIGM,VECX,&
                    X,Y,Z,&
                    N,NLX,NLY,NLZ, WEIGHT,&
                    NU,NV,NW, UG,VG,WG, &
                    TNEW,TOLD,DENPT,DENPT,&
                    DISOPN,DT,THETA,BETA, D3,DCYL,&
                    PARA,halo_tag,NNODP, .false.&
                    , NOBCT, BCT1, BCT2&
                    , NITS, NDISOT2, NEWMES2, MATTSTOR&
                    , NSUBNVLOC,NUSUB,NVSUB,NWSUB, option_path&
                    )
  
          end if
       ENDIF


        ! apply the dirichlet bcs
        IDIM=1
        ewrite(2,*) 'R2NORM OF VECX1:',R2NORM(VECX,NNODP,PARA)
  
        ewrite(1,*)'calling boucon'
        CALL BOUCON(IDIM, NONODS, &
              DM1,NOBCT,BCT1,BCT2,VECX,&
              DM1,NOBCT,BCT1,BCT2,VECX,&
              DM1,NOBCT,BCT1,BCT2,VECX    )
        ewrite(3,*) 'R2NORM OF DM1:',R2NORM(DM1,NNODP,PARA)
        ewrite(3,*) 'R2NORM OF VECX-:',R2NORM(VECX,NNODP,PARA)

       !     Solve the equation for new T. 
       RMAXT=maxval(T)
       CALL ALLMAX(RMAXT)

       RMINT=minval(T)
       CALL ALLMIN(RMINT)

       RMAXNU=maxval(NU)
       CALL ALLMAX(RMAXNU)

       RMINNU=minval(NU)
       CALL ALLMIN(RMINNU)

       RMAXNV=maxval(NV)
       CALL ALLMAX(RMAXNV)

       RMINNV=minval(NV)
       CALL ALLMIN(RMINNV)

       ewrite(2,*)'***** MAX,MIN T IS =====',RMAXT,RMINT
       ewrite(2,*)'RMAXNU,RMINNU,RMAXNV,RMINNV:',&
            RMAXNU,RMINNU,RMAXNV,RMINNV

       UZAWA=CGSOLQ

!       ewrite_minmax(vecx)
       ewrite_minmax(bigm)
       ewrite(2,*) "R2NORM(VECX)",PARA,NONODS,R2NORM(VECX,NNODP,PARA)
       ewrite(2,*) 'in advdif, SYM:',SYM
       IGUESS=.FALSE.

       ! f and savdia seem to be unused by advdif
       allocate(f(nonods))
       allocate(savdia(nonods))
       IF(NITS.GE.2) THEN     
          IGUESS=.TRUE.
          CALL SOLMOM(T,TNEW,IGUESS,VECX,f,DT,&
               DM1,BIGM,FINDRM,COLM,CENTRM,NCOLM, &
               savdia,&
               NONODS,&
               GMRESQ,&
               PARA,halo_tag, NNODP,&
               option_path)
          
       ELSE
          IGUESS=.FALSE.
          CALL SOLMOM(T,T,IGUESS,VECX,f,DT,&
               DM1,BIGM,FINDRM,COLM,CENTRM,NCOLM, &
               savdia,&
               NONODS,&
               GMRESQ,&
               PARA,halo_tag,NNODP,&
               option_path)
       ENDIF
       deallocate(f)
       deallocate(savdia)
       
       IF(NSUBTLOC.NE.0) THEN
          ! This sub calculates the inner element model TSUB from T
          CALL GETSUBSCALE(TSUB,T,DINVSOUSUBTSTORE,DMATINVCSTORE,&
               NONODS,TOTELE,NLOC,NDGLNO,NSUBTLOC)
          TNEWSUB=TSUB
          DEALLOCATE(DMATINVCSTORE)
          DEALLOCATE(DINVSOUSUBTSTORE)
       ENDIF

       !     Add in the 1st order spatial derivative -transport. ******
       IF(ABS(NDISOT).GE.10) THEN
          !     Use and extended halo in the high resolution method...
          IF(D3.AND.(NLOC.EQ.4)) THEN
             IF(PARA.EQ.1) THEN
                CALL HALGET(TNEW,NONODS,NONODS,NNODP&
                     ,halo_tag_p)
                CALL HALGET(TOLD,NONODS,NONODS,NNODP&
                     ,halo_tag_p)
                IF(NITS.EQ.1) THEN
                   CALL HALGET(NU,NONODS,NONODS,NNODP,halo_tag_p)
                   CALL HALGET(NV,NONODS,NONODS,NNODP,halo_tag_p)
                   IF(D3) CALL HALGET(NW,NONODS,NONODS,NNODP&
                        ,halo_tag_p)
                   CALL HALGET(UG,NONODS,NONODS,NNODP,halo_tag_p)
                   CALL HALGET(VG,NONODS,NONODS,NNODP,halo_tag_p)
                   IF(D3) CALL HALGET(WG,NONODS,NONODS,NNODP&
                        ,halo_tag_p)
                ENDIF
             ENDIF
          ENDIF
       ENDIF

    END DO
      
      TTNEW=T

42  deallocate(sndgln) 
    deallocate(tsndgl)
    deallocate(salphe)
    deallocate(taire)
    deallocate(n)
    deallocate(nlx)
    deallocate(nly)
    deallocate(nlz)
    deallocate(sn)
    deallocate(snlx)
    deallocate(snly)
    deallocate(snlz)
    deallocate(denpt)
    ! Drop the reference and deallocate (if sinking velocity present).
    call deallocate(U_nl)

    ! Copy viscosity back in case it's changed.
    if (has_tensor_field(state(1),trim(field_name)//"Diffusivity")) then
       call copy_tensor_from_legacy(diffusivity, MUPTXX, MUPTXY, MUPTXZ, &
            MUPTYY, MUPTYZ, MUPTZZ)
    end if
        
#ifdef DDEBUG
    if(any(abs(zero_real_array) > 0.0)) then
      ewrite(0, *) "The zero real array seems to have been written to!"
      ewrite(0, *) "Max val: ", maxval(zero_real_array)
      ewrite(0, *) "Min val: ", minval(zero_real_array)
    end if
#endif
    deallocate(zero_real_array)

    call deallocate(M)
    call deallocate(A_m)
    call deallocate(rhs)
    call deallocate(masslump)
    if(ABS(NDISOT).GE.10) then
      call deallocate(l_told_field)
      call deallocate(l_tnew_field)
    end if

    ewrite(1, *) "Exiting advdif"

  end subroutine advdif

  LOGICAL FUNCTION LCONVR(GLORER,RELATI,NEW,PREV,TOTELE,PARA)
    INTEGER TOTELE
    REAL GLORER
    LOGICAL RELATI
    REAL NEW(TOTELE),PREV(TOTELE)
    INTEGER PARA
    ! Local variables...
    REAL TOLER,INFINY
    PARAMETER(TOLER=1.E-10,INFINY=1.E+20) 
    REAL RMIN,RMAX,MAXDIF,RLCONV
    INTEGER I
    RMIN= INFINY
    RMAX=-INFINY
    MAXDIF=0.
    do I=1,TOTELE
       RMIN=MIN(RMIN,NEW(I))
       RMAX=MAX(RMAX,NEW(I))
       MAXDIF=MAX(MAXDIF,ABS(PREV(I)-NEW(I)))
    END DO
    IF(PARA.NE.0) THEN
       CALL ALLMIN(RMIN)
       CALL ALLMAX(RMAX)
       CALL ALLMAX(MAXDIF)
    ENDIF
    IF(RELATI) THEN
       ! Relative error...
       LCONVR=(MAXDIF/MAX(TOLER,1.001*RMAX-RMIN).LT.GLORER)
       ewrite(3,*) 'MAXDIF,1.001*RMAX-RMIN:',&
            MAXDIF,1.001*RMAX-RMIN
       ewrite(3,*) 'MAXDIF/MAX(TOLER,1.001*RMAX-RMIN):',&
            MAXDIF/MAX(TOLER,1.001*RMAX-RMIN)
    ELSE
       ! Absoute error...
       LCONVR=(MAXDIF.LT.GLORER)
    ENDIF
    !
    IF(PARA.NE.0) THEN
       RLCONV=0.
       IF(LCONVR) RLCONV=1.0
       CALL ALLMAX(RLCONV)
       IF(RLCONV.GT.0.5) THEN
          LCONVR=.TRUE.
       ELSE
          LCONVR=.FALSE.
       ENDIF
    ENDIF
    ewrite(3,*) '****MAXDIF,RMAX,RMIN,GLORER:',&
         MAXDIF,RMAX,RMIN,GLORER
    !
    RETURN
  END function LCONVR

  subroutine get_copied_field(fieldname, state)

    type(state_type), intent(in) :: state
    character(len=*), intent(in) :: fieldname

    type(scalar_field), pointer :: copiedfield
    type(scalar_field), pointer :: tmpfield
    character(len=OPTION_PATH_LEN) :: tmpstring

    if(trim(fieldname)=="CopiedField") then
      copiedfield=>extract_scalar_field(state, "CopiedField")
      call get_option(trim(copiedfield%option_path)//"/prognostic/copy_from_field", &
                              tmpstring)
      tmpfield=>extract_scalar_field(state, "Old"//trim(tmpstring))
      call set(copiedfield, tmpfield)
    end if

  end subroutine get_copied_field

end module AdvectionDiffusion
