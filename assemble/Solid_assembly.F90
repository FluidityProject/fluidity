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
!************************************************************************
module Solid_assembly
  use quadrature
  use elements
  use fldebug
  use state_module
  use fields
  use elements
  use sparsity_patterns_meshes
  use spud
  implicit none

  private
  public :: SOLID3D, STEVAL

  contains

  SUBROUTINE SOLID3D( VECX, VECY, VECZ, &
                      BIGM, &
                      ML,  &
                      NBIGM, NCOLM, & 
                      NLOC, NGI, &
                      DT, THETA, &
                      LUMP, &
                      ! plastic strains:
                      STVPXX, STVPYY, STVPZZ, &
                      STVPYZ, STVPXZ, STVPXY, &
                      ! solid switch:
                      SOLIDS, &
                      state )
  
    ! =============================================================
    ! Subroutine to construct the matrices BIGM1, C1/2/3T, ML
    ! for a variety of solid types.
    ! =============================================================
    ! =============================================================
    ! inputs 
    ! =============================================================
    type(state_type), intent(inout) :: state

    INTEGER, INTENT(IN) :: NLOC, NGI, NBIGM, NCOLM
    ! nloc  = number of displacement nodes per element
    ! ngi   = number of gauss points (for current element)
    ! ncolm = length of each block of bigm
        
    REAL, INTENT(IN) :: DT, THETA
    ! dt    = time step
    ! theta = parameter \in [0,1] to select between explicit and implicit 
    !         treatment of stiffness terms
    !         theta = 0.0 -> explicit
    !         theta = 0.5 -> crank-nicolson
    !         theta = 1.0 -> implicit

    LOGICAL, INTENT(IN) :: LUMP
    ! lump   = if true lump the mass matrix terms (mii) to the diagonal of bigm(nbigm)
        
    INTEGER, INTENT(IN) :: SOLIDS
    ! solids = dictates behaviour of material:
    !          1: elastic... this sub does nothing
    !          2: kelvin...  this sub does nothing
    !          3: bingham... this sub corrects plastic overstep
    
    ! --------------------------------------start of add by crgw 31/03/06
    ! plasticity stuff:
    ! plastic strains:
    REAL, INTENT(IN) ::  STVPXX(:), STVPYY(:), STVPZZ(:)
    ! stvpxx(totele*ngi) = plastic strain in the x direction 
    !                      on a plane perpendicular to the x axis
    ! stvpyy(totele*ngi) = plastic strain in the y direction 
    !                      on a plane perpendicular to the y axis
    ! stvpzz(totele*ngi) = plastic strain in the z direction 
    !                      on a plane perpendicular to the z axis
          
    REAL, INTENT(IN) ::  STVPYZ(:), STVPXZ(:), STVPXY(:)
    ! stvpyz(totele*ngi) = plastic strain in the z direction 
    !                      on a plane perpendicular to the y axis
    ! stvpxz(totele*ngi) = plastic strain in the z direction 
    !                      on a plane perpendicular to the x axis
    ! stvpxy(totele*ngi) = plastic strain in the y direction 
    !                      on a plane perpendicular to the x axis
    ! --------------------------------------end of add by crgw 31/03/06
  
    ! =============================================================
    ! outputs
    ! =============================================================
    REAL, INTENT(INOUT) :: VECX(:), VECY(:), VECZ(:)
    ! vecx(nonods) = upper third of r.h.s vector (x force components)
    ! vecy(nonods) = middle third of r.h.s vector (y force components)
    ! vecz(nonods) = lower third of r.h.s vector (z force components)
    
    REAL, INTENT(INOUT) :: BIGM(:), ML(:)
    ! bigm(nbigm)   = left hand side matrix (sparsely stored as a vector) in equation:
    !                 bigm*[u;v;w] = [vecx;vecy;vecz]
    ! ml(nonods)    = lumped mass matrix
    
    ! =============================================================
    ! local variables
    ! =============================================================
    ! fractions/surds/parameters etc.:
    REAL :: PI
    PARAMETER(PI=3.1415926538)
    REAL :: TWOTHI, FOURTHI
    REAL :: ONESIX, ONETHI, ONEHAL
    PARAMETER(TWOTHI = 2./3., FOURTHI = 4./3.)
    PARAMETER(ONESIX = 1./6., ONETHI = 1./3., ONEHAL = 1./2.)
    ! end of parameters
    
    REAL :: MUGIXX, MUGIZZ
    ! mugixx = muptxx(snonod) evaluated at a gauss pt
    ! mugizz = muptzz(snonod) evaluated at a gauss pt
          
    INTEGER :: ELE, ILOC, JLOC, GI, COUNT
    ! ele   = integer counter through the elements
    ! iloc  = integer counter through the nodes in each element (a)
    ! jloc  = integer counter through the nodes in each element (b)
    ! gi    = integer counter through the gauss points in each element
    ! count = integer used to cycle through rows of bigm(nbigm) 
    !         to find correct position
          
    REAL :: NX(NLOC,NGI), NY(NLOC,NGI), NZ(NLOC,NGI)
    ! nx(nloc,ngi) = derivative of displacement shape function in global x direction
    !                (evaluated at a gauss point around a node)
    ! ny(nloc,ngi) = derivative of displacement shape function in global y direction
    !                (evaluated at a gauss point around a node)
    ! nz(nloc,ngi) = derivative of displacement shape function in global z direction
    !                (evaluated at a gauss point around a node)
          
    REAL :: UD(NGI), VD(NGI), WD(NGI)
    ! ud(ngi) = difference between new iterate velocity and grid velocity in x direction
    !           evaluated at the gauss points
    ! vd(ngi) = difference between new iterate velocity and grid velocity in y direction
    !           evaluated at the gauss points
    ! wd(ngi) = difference between new iterate velocity and grid velocity in z direction
    !           evaluated at the gauss points
        
    ! streamline-upwind/petrov-galerkin terms:
    REAL :: HXGI, HYGI, HZGI, HOVERQ
    ! hxgi = abs(first entry of the inverse of the jacobian times [ud; vd; wd])
    ! hygi = abs(second entry of the inverse of the jacobian times [ud; vd; wd])
    ! hzgi = abs(third entry of the inverse of the jacobian times [ud; vd; wd])
          
    REAL :: R1,R2,UU, HL
    ! uu,r1,r2 = terms used in assembly of petrov-galerkin term
    ! hl       = petrov galerkin term
    
    REAL :: GAMMA(NGI)
    ! gamma(ngi) = measure of grid size for streamline upwind/petrov galerkin term
        
    REAL :: AGI,BGI,CGI
    REAL :: DGI,EGI,FGI
    REAL :: GGI,HGI,KGI
    ! [agi, bgi, cgi;
    !  dgi, egi, fgi;
    !  ggi, hgi, kgi] = entries of jacobi matrix
    
    REAL :: DETJ, INVDET, DETWEI(NGI), DDETWE(NGI)
    ! detj   = determinant of jacobi matrix (jacobian)
    ! invdet = 1./detj 
    ! detwei(ngi) = determinant of jacobian times weight at the gauss point
    !               (evaluated by subroutine detnlxr)
    ! ddetwe(ngi) = detwei(at a gauss pt)*density(at a gauss pt)
          
    REAL :: A11,A12,A13
    REAL :: A21,A22,A23
    REAL :: A31,A32,A33
    ! [a11, a12, a13;
    !  a21, a22, a23;
    !  a31, a32, a33] = entries of inverse jacobi matrix
    
    INTEGER :: IBL11,IBL12,IBL13
    INTEGER :: IBL21,IBL22,IBL23
    INTEGER :: IBL31,IBL32,IBL33
    ! ibl11/12/13/21/22/23/31/32/33 = indexing to "blocks" of "matrix" bigm(nbigm)
          
    LOGICAL :: BLKSYM
    ! blksym = dictates symmetry of bigm(nbigm):
    !          when blksym = .true.:
    !          [ibl11 ibl12 ibl13
    !           ibl12 ibl22 ibl23
    !           ibl13 ibl23 ibl33]
    !          when blksym = .false.:
    !          [ibl11 ibl12 ibl13
    !           ibl21 ibl22 ibl23
    !           ibl31 ibl32 ibl33]
          
    REAL :: RLUMP
    ! rlump = controls lumping of mass matrix (dependent on lump)
          
    ! indexing integers:
    INTEGER :: IGLS, IGLV, IGLX, ICENT
    INTEGER :: JGLS
    INTEGER :: SGI
    INTEGER :: INOD, JNOD, POSMAT
    ! end of indexing integers
          
    REAL :: TWOTHIMU, FOURTHIMU
    ! twothimu  = -2./3.*mu_(e/v)
    ! fourthimu =  4./3.*mu_(e/v)
          
    REAL :: D11(NGI), D12(NGI), D13(NGI), D14(NGI), D15(NGI), D16(NGI)
    REAL :: D22(NGI), D23(NGI), D24(NGI), D25(NGI), D26(NGI)
    REAL :: D33(NGI), D34(NGI), D35(NGI), D36(NGI)
    REAL :: D44(NGI), D45(NGI), D46(NGI)
    REAL :: D55(NGI), D56(NGI)
    REAL :: D66(NGI)
    ! d11/12/13/14/15/16/22/23/24/25/26/33/34/35/36/44/45/46/55/56/66
    !  = components of the symmetric matrix d relating stress to deviatoric strain
    !  - for a linear isotropic 3-d solid the matrix d will have the form:
    !    4./3.*mu_e     -2./3.mu_e     -2./3.mu_e    0.    0.   0.
    !    -2./3.mu_e     4./3.*mu_e     -2./3.mu_e    0.    0.   0.
    !    -2./3.mu_e     -2./3.mu_e     4./3.*mu_e    0.    0.   0.
    !        0.             0.            0.        mu_e   0.   0.
    !        0.             0.            0.         0.   mu_e  0.
    !        0.             0.            0.         0.    0.  mu_e
    
    REAL :: F11(NGI), F12(NGI), F13(NGI), F14(NGI), F15(NGI), F16(NGI)
    REAL :: F22(NGI), F23(NGI), F24(NGI), F25(NGI), F26(NGI)
    REAL :: F33(NGI), F34(NGI), F35(NGI), F36(NGI)
    REAL :: F44(NGI), F45(NGI), F46(NGI)
    REAL :: F55(NGI), F56(NGI)
    REAL :: F66(NGI)
    ! f11/12/13/14/15/16/22/23/24/25/26/33/34/35/36/44/45/46/55/56/66
    !  = components of the symmetric matrix f relating stress to deviatoric strain rate
    !  - for a linear isotropic 3-d material the matrix f will have the form:
    !    4./3.*mu_v     -2./3.mu_v     -2./3.mu_v    0.    0.   0.
    !    -2./3.mu_v     4./3.*mu_v     -2./3.mu_v    0.    0.   0.
    !    -2./3.mu_v     -2./3.mu_v     4./3.*mu_v    0.    0.   0.
    !        0.             0.            0.        mu_v   0.   0.
    !        0.             0.            0.         0.   mu_v  0.
    !        0.             0.            0.         0.    0.  mu_v
    
    REAL :: K11, K12, K13
    REAL :: K21, K22, K23
    REAL :: K31, K32, K33
    ! k11/12/13/21/22/23/31/32/33
    !  = components of the elastic stress-strain stiffness matrix
          
    REAL :: C11, C12, C13
    REAL :: C21, C22, C23
    REAL :: C31, C32, C33
    ! c11/12/13/21/22/23/31/32/33
    !  = components of the viscous stress-strain stiffness matrix
          
    ! conservation of momentum matrix components:
    REAL :: MII, NII, NN
    ! mii   = mass matrix diagonal matrix component (including density)
    ! nii   = non-linear advection diagonal matrix component for velocity nodes
    ! nn    = mass matrix diagonal matrix component (excluding density)
          
    REAL :: NDENGI(NGI)
    ! ndengi(ngi)   = new iterate density (to be) evaluated at gauss points
    
    ! parameters to set the model type:     
    ! nb: default behaviour is for both to be false, when the model simulates a purely
    !     linearly elastic material
    LOGICAL :: KELVIN, BINGHAM
    ! kelvin  = simulates a kelvin substance with elastic and viscous components
    !           in parallel
    !           this option needs to be set to true for multi-material modelling where
    !           one section of the mass fraction is a solid and the other a liquid
    ! bingham = simulates a bingham substance which initially behaves elastically until
    !           a yield stress is reached at which point it begins to flow
    !           nb: think this only works for lagrangian (one material) models at the moment
  
    ! --------------------------------------start of add by crgw 31/03/06
    ! plasticity stuff:
    REAL :: BD11GI, BD12GI, BD13GI
    REAL :: BD14GI, BD15GI, BD16GI
    REAL :: BD21GI, BD22GI, BD23GI
    REAL :: BD24GI, BD25GI, BD26GI
    REAL :: BD31GI, BD32GI, BD33GI
    REAL :: BD34GI, BD35GI, BD36GI
    ! bd11/12/13/14/15/16/21/22/23/24/25/26/31/32/33/34/35/36gi
    !  = components of 3x6 "unfinished stiffness" matrix for evaluating
    !    rhs stiffness terms arising from the plastic strains
          
    REAL :: SPCONX, SPCONY, SPCONZ
    ! spconx = rhs x contribution arising from the multiplications of bd__gi with
    !          the plastic strains
    ! spcony = rhs y contribution arising from the multiplications of bd__gi with
    !          the plastic strains
    ! spconz = rhs z contribution arising from the multiplications of bd__gi with
    !          the plastic strains
    ! --------------------------------------end of add by crgw 31/03/06
    

    type(vector_field), pointer :: velocity, nlvelocity, gridvelocity, velocitysource
    type(scalar_field) :: velocityx, velocityy, velocityz, &
                          nlvelocityx, nlvelocityy, nlvelocityz, &
                          gridvelocityx, gridvelocityy, gridvelocityz
    real, dimension(:), pointer :: u, v, w, nu, nv, nw, ug, vg, wg, sourcx, sourcy, sourcz

    type(vector_field), pointer :: displacement
    type(scalar_field) :: displacementx, displacementy, displacementz
    real, dimension(:), pointer :: xdisp, ydisp, zdisp
  
    type(vector_field), pointer :: coordinate
    type(scalar_field) :: coordinatex, coordinatey, coordinatez
    real, dimension(:), pointer :: x, y, z
  
    type(tensor_field), pointer :: elasticity, viscosity
    type(scalar_field) :: elasticityxx, viscosityxx
    real, dimension(:), pointer :: muptxx, muptzz
  
    type(scalar_field), pointer :: density
    real, dimension(:), pointer :: denpt

    integer :: totele, nonods, stat
    logical :: have_gravity
    real :: gravity_mag

    integer, dimension(:), pointer :: ndglno, sondgl, vondgl, xondgl

    type(element_type), pointer :: shape

    REAL :: N(NLOC,NGI), NLX(NLOC,NGI), NLY(NLOC,NGI), NLZ(NLOC,NGI)
    ! n(nloc,ngi)   = displacement shape function 
    ! nlx(nloc,ngi) = derivative of displacement shape function in local x direction
    !                 (evaluated at a gauss point around a node)
    ! nly(nloc,ngi) = derivative of displacement shape function in local y direction
    !                 (evaluated at a gauss point around a node)
    ! nlz(nloc,ngi) = derivative of displacement shape function in local z direction
    !                 (evaluated at a gauss point around a node)
    
    real, dimension(:), pointer :: weight
    ! weight(ngi) = weights for gaussian quadrature at each of the gauss points

    type(csr_sparsity), pointer :: velocity_sparsity
    integer, dimension(:), pointer :: findrm, colm, centrm

    real, dimension(:), allocatable, target :: zero ! can't believe I'm doing this but this will just be temporary

    ! =============================================================
    ! This subroutine works in two halves.  The first evaluates the stress
    ! state from the previous iterate and determines if it lies within the
    ! yield envelope or not.  The second only occurs if it does not, in which
    ! case a local (at a single Gauss point) integration loop is performed to
    ! correct the overstep back to the yield surface.
    ! 
    ! called from fluids.f
    ! 
    ! Description                                   Programmer      Date
    ! ==================================================================
    ! Original version............................. CRGW        2006-06-02
    ! =============================================================
  
    ! =============================================================
    ! 
    ! The subroutine operates almost entirely within a loop through the elements.
    ! 
    ! The first half of this loop goes through the Gauss points first then the
    ! nodes (of the current element) evaluating the values of various parameters 
    ! and functions at the Gauss points.
    ! 
    ! The second half consists of two nested loops through the nodes of the element
    ! and a third loop through the Gauss points evaluating the integrals of the 
    ! parameters and the relationships (mostly shape functions) between the node pairs.
    ! 
    ! These integrals are then used to form the matrices BIGM1(NBIGM) (called locally BIGM),
    ! C1/2/3T and ML
    ! 
    ! Called from ASSNAV
    ! 
    ! Description                                   Programmer      Date
    ! ==================================================================
    ! Original version for elastic, Kelvin and 
    !   Bingham solids............................... CRGW       14/03/06
    ! Added compressibility........................... CRGW       14/05/06
    ! Major correction to Bingham solid (now mostly
    !     done in STEVAL)............................. CRGW       02/06/06
    ! =============================================================
  
    ewrite(1,*) 'IN SOLID3D()'
  
    displacement=>extract_vector_field(state, "Displacement")

    displacementx=extract_scalar_field_from_vector_field(displacement, 1)
    xdisp=>displacementx%val

    displacementy=extract_scalar_field_from_vector_field(displacement, 2)
    ydisp=>displacementy%val

    displacementz=extract_scalar_field_from_vector_field(displacement, 3)
    zdisp=>displacementz%val


    velocity=>extract_vector_field(state, "Velocity")

    velocityx=extract_scalar_field_from_vector_field(velocity, 1)
    u=>velocityx%val

    velocityy=extract_scalar_field_from_vector_field(velocity, 2)
    v=>velocityy%val

    velocityz=extract_scalar_field_from_vector_field(velocity, 3)
    w=>velocityz%val

    totele = element_count(velocity)
    nonods = node_count(velocity)

    allocate(zero(nonods))
    zero = 0.0
  
    nlvelocity=>extract_vector_field(state, "NonlinearVelocity")

    nlvelocityx=extract_scalar_field_from_vector_field(nlvelocity, 1)
    nu=>nlvelocityx%val

    nlvelocityy=extract_scalar_field_from_vector_field(nlvelocity, 2)
    nv=>nlvelocityy%val

    nlvelocityz=extract_scalar_field_from_vector_field(nlvelocity, 3)
    nw=>nlvelocityz%val


    gridvelocity=>extract_vector_field(state, "GridVelocity")

    gridvelocityx=extract_scalar_field_from_vector_field(gridvelocity, 1)
    ug=>gridvelocityx%val

    gridvelocityy=extract_scalar_field_from_vector_field(gridvelocity, 2)
    vg=>gridvelocityy%val

    gridvelocityz=extract_scalar_field_from_vector_field(gridvelocity, 3)
    wg=>gridvelocityz%val

    density=>extract_scalar_field(state, "Density")
    denpt=>density%val

    call get_option("/physical_parameters/gravity/magnitude", gravity_mag, stat)
    have_gravity = stat==0

    if(have_gravity) then
      velocitysource=>extract_vector_field(state, "GravityDirection")
      
      allocate(sourcx(nonods))
      allocate(sourcy(nonods))
      allocate(sourcz(nonods))
      
      select case(velocitysource%field_type)
      case(FIELD_TYPE_CONSTANT)
        sourcx = velocitysource%val(1)%ptr(1)*gravity_mag*denpt
        sourcy = velocitysource%val(2)%ptr(1)*gravity_mag*denpt
        sourcz = velocitysource%val(3)%ptr(1)*gravity_mag*denpt
      case default
        sourcx = velocitysource%val(1)%ptr*gravity_mag*denpt
        sourcy = velocitysource%val(2)%ptr*gravity_mag*denpt
        sourcz = velocitysource%val(3)%ptr*gravity_mag*denpt
      end select
    
    else

      sourcx=>zero
      sourcy=>zero
      sourcz=>zero

    end if

    coordinate=>extract_vector_field(state, "Coordinate")

    coordinatex=extract_scalar_field_from_vector_field(coordinate, 1)
    x=>coordinatex%val

    coordinatey=extract_scalar_field_from_vector_field(coordinate, 2)
    y=>coordinatey%val

    coordinatez=extract_scalar_field_from_vector_field(coordinate, 3)
    z=>coordinatez%val


    elasticity=>extract_tensor_field(state, "Elasticity")

    elasticityxx=extract_scalar_field_from_tensor_field(elasticity, 1, 1)
    muptzz=>elasticityxx%val


    viscosity=>extract_tensor_field(state, "Viscosity", stat)

    if(stat==0) then
      viscosityxx=extract_scalar_field_from_tensor_field(viscosity, 1, 1)
      muptxx=>viscosityxx%val
    else
      muptxx=>zero
    end if


    ndglno => velocity%mesh%ndglno
    sondgl => velocity%mesh%ndglno
    vondgl => velocity%mesh%ndglno
    xondgl => coordinate%mesh%ndglno

    shape=>ele_shape(velocity,1)
    weight=>shape%quadrature%weight

    call extract_old_element(shape, N, NLX, NLY, NLZ)

    velocity_sparsity=>get_csr_sparsity_firstorder(state, velocity%mesh, velocity%mesh)
    findrm=>velocity_sparsity%findrm
    centrm=>velocity_sparsity%centrm
    colm=>velocity_sparsity%colm

    ewrite_minmax(MUPTXX)
    ewrite_minmax(MUPTZZ)
        
    KELVIN = .FALSE.
    BINGHAM = .FALSE.
    IF(SOLIDS.EQ.2) KELVIN=.TRUE.
    IF(SOLIDS.EQ.3) BINGHAM=.TRUE.
  
    ewrite(3,*)'SOLIDS, KELVIN, BINGHAM = ', &
            SOLIDS, KELVIN, BINGHAM
    ewrite(3,*) 'NLOC, NGI = ', NLOC, NGI
        
  ! matrix column locations
    IF(NBIGM.EQ.9*NCOLM) THEN
      BLKSYM=.FALSE.
    ELSE
      BLKSYM=.TRUE.
    ENDIF
    IF(BLKSYM) THEN
      IBL11=0
      IBL12=1*NCOLM
      IBL13=3*NCOLM
      IBL21=1*NCOLM
      IBL22=2*NCOLM
      IBL23=4*NCOLM
      IBL31=3*NCOLM
      IBL32=4*NCOLM
      IBL33=5*NCOLM
    ELSE
      IBL11=0
      IBL12=1*NCOLM
      IBL13=2*NCOLM
      IBL21=3*NCOLM
      IBL22=4*NCOLM
      IBL23=5*NCOLM
      IBL31=6*NCOLM
      IBL32=7*NCOLM
      IBL33=8*NCOLM
    ENDIF
    
    ! clear all vectors and matrices
    VECX  = 0.
    VECY  = 0.
    VECZ  = 0.
    BIGM  = 0.
    ML    = 0.
    GAMMA = 0.
        
    RLUMP = 0.
    ! if lumping the mass matrix terms in bigm
    IF(LUMP) RLUMP = 1.

    ewrite_minmax(xdisp)
    ewrite_minmax(ydisp)
    ewrite_minmax(zdisp)

    ewrite_minmax(u)
    ewrite_minmax(v)
    ewrite_minmax(w)

    ewrite_minmax(nu)
    ewrite_minmax(nv)
    ewrite_minmax(nw)

    ewrite_minmax(ug)
    ewrite_minmax(vg)
    ewrite_minmax(wg)

    ewrite_minmax(denpt)

    ewrite_minmax(stvpxx)
    ewrite_minmax(stvpyy)
    ewrite_minmax(stvpzz)
    ewrite_minmax(stvpyz)
    ewrite_minmax(stvpxz)
    ewrite_minmax(stvpxy)

    ewrite_minmax(sourcx)
    ewrite_minmax(sourcy)
    ewrite_minmax(sourcz)

    ! loop through elements
    DO ELE=1,TOTELE
          
      ! loop through gauss points evaluating parameters at them
      DO GI=1,NGI
          
        SGI = (ELE-1)*NGI + GI
  
        ! zero the jacobian components         
        AGI = 0.
        BGI = 0.
        CGI = 0.
            
        DGI = 0.
        EGI = 0.
        FGI = 0.
        
        GGI = 0.
        HGI = 0.
        KGI = 0.
        
        ! zero non-linear advection terms
        UD(GI) = 0.
        VD(GI) = 0.
        WD(GI) = 0.
            
        ! zero temporary material parameters (prior to summation)
        MUGIXX = 0.
        MUGIZZ = 0.
  
        NDENGI(GI) = 0.
            
        ! loop through element nodes
        DO ILOC=1,NLOC
          ! indexing
          IGLS = SONDGL((ELE-1)*NLOC + ILOC)
          IGLX = XONDGL((ELE-1)*NLOC + ILOC)
          IGLV = VONDGL((ELE-1)*NLOC + ILOC)
          
          ! sum the components of the jacobian           
          AGI = AGI + NLX(ILOC,GI)*X(IGLX) 
          BGI = BGI + NLX(ILOC,GI)*Y(IGLX) 
          CGI = CGI + NLX(ILOC,GI)*Z(IGLX) 
              
          DGI = DGI + NLY(ILOC,GI)*X(IGLX) 
          EGI = EGI + NLY(ILOC,GI)*Y(IGLX) 
          FGI = FGI + NLY(ILOC,GI)*Z(IGLX) 
          
          GGI = GGI + NLZ(ILOC,GI)*X(IGLX) 
          HGI = HGI + NLZ(ILOC,GI)*Y(IGLX) 
          KGI = KGI + NLZ(ILOC,GI)*Z(IGLX) 
              
          ! non-linear advection velocties 
          UD(GI) = UD(GI) + N(ILOC,GI)*(NU(IGLV)-UG(IGLV))
          VD(GI) = VD(GI) + N(ILOC,GI)*(NV(IGLV)-VG(IGLV))
          WD(GI) = WD(GI) + N(ILOC,GI)*(NW(IGLV)-WG(IGLV))
          
          ! sum parameters at nodes to get the value at the gauss point
          MUGIXX = MUGIXX + N(ILOC,GI)*MUPTXX(IGLS)
          MUGIZZ = MUGIZZ + N(ILOC,GI)*MUPTZZ(IGLS)
          
          ! density
          NDENGI(GI) = NDENGI(GI) + N(ILOC,GI)*DENPT(IGLS)
        
        ! iloc=1,nloc
        ENDDO
  
        ! determinant of jacobi matrix (determinant)          
        DETJ = AGI*(EGI*KGI - FGI*HGI) - &
              BGI*(DGI*KGI - FGI*GGI) + &
              CGI*(DGI*HGI - EGI*GGI)
        
        DETWEI(GI)=abs(DETJ)*WEIGHT(GI)
        DDETWE(GI)=NDENGI(GI)*DETWEI(GI)
        INVDET    =1./DETJ
        
        ! calculate the components of the jacobi inverse
        A11= (EGI*KGI - FGI*HGI)*INVDET
        A21=-(DGI*KGI - FGI*GGI)*INVDET
        A31= (DGI*HGI - EGI*GGI)*INVDET
        
        A12=-(BGI*KGI - CGI*HGI)*INVDET
        A22= (AGI*KGI - CGI*GGI)*INVDET
        A32=-(AGI*HGI - BGI*GGI)*INVDET
        
        A13= (BGI*FGI - CGI*EGI)*INVDET
        A23=-(AGI*FGI - CGI*DGI)*INVDET
        A33= (AGI*EGI - BGI*DGI)*INVDET
        
        ! evaluate the derivatives of the velocity shape functions at the gauss points
        DO ILOC = 1, NLOC
          NX(ILOC, GI) = A11*NLX(ILOC, GI) + A12*NLY(ILOC, GI) + A13*NLZ(ILOC, GI)
          NY(ILOC, GI) = A21*NLX(ILOC, GI) + A22*NLY(ILOC, GI) + A23*NLZ(ILOC, GI)
          NZ(ILOC, GI) = A31*NLX(ILOC, GI) + A32*NLY(ILOC, GI) + A33*NLZ(ILOC, GI)
        ! iloc=1,nloc
        ENDDO
            
        ! calculate gamma for use in petrov-galerkin term
        HXGI = ABS(A11*UD(GI) + A12*VD(GI) + A13*WD(GI))
        HYGI = ABS(A21*UD(GI) + A22*VD(GI) + A23*WD(GI))
        HZGI = ABS(A31*UD(GI) + A32*VD(GI) + A33*WD(GI))
                              
        HOVERQ = 2./MAX(HXGI, HYGI, HZGI, 1.E-7)
                        
        GAMMA(GI) = ONEHAL*HOVERQ 
        
        FOURTHIMU = FOURTHI*MUGIZZ 
        TWOTHIMU  = -TWOTHI*MUGIZZ
              
        D11(GI) = FOURTHIMU
        D12(GI) = TWOTHIMU
        D13(GI) = TWOTHIMU
        D14(GI) = 0.
        D15(GI) = 0.
        D16(GI) = 0.
        
        D22(GI) = FOURTHIMU
        D23(GI) = TWOTHIMU
        D24(GI) = 0.
        D25(GI) = 0.
        D26(GI) = 0.
          
        D33(GI) = FOURTHIMU
        D34(GI) = 0.
        D35(GI) = 0.
        D36(GI) = 0.
        
        D44(GI) = MUGIZZ
        D45(GI) = 0.
        D46(GI) = 0.
        
        D55(GI) = MUGIZZ
        D56(GI) = 0.
        
        D66(GI) = MUGIZZ
        
        IF(KELVIN) THEN
          ! viscosities
          FOURTHIMU = FOURTHI*MUGIXX
          TWOTHIMU  = -TWOTHI*MUGIXX
          
          F11(GI) = FOURTHIMU
          F12(GI) = TWOTHIMU
          F13(GI) = TWOTHIMU
          F14(GI) = 0.
          F15(GI) = 0.
          F16(GI) = 0.
          
          F22(GI) = FOURTHIMU
          F23(GI) = TWOTHIMU
          F24(GI) = 0.
          F25(GI) = 0.
          F26(GI) = 0.
          
          F33(GI) = FOURTHIMU
          F34(GI) = 0.
          F35(GI) = 0.
          F36(GI) = 0.
          
          F44(GI) = MUGIXX
          F45(GI) = 0.
          F46(GI) = 0.
          
          F55(GI) = MUGIXX
          F56(GI) = 0.
          
          F66(GI) = MUGIXX
          
          ! if(kelvin)
        ENDIF
            
      ! gi=1,ngi
      ENDDO
  
      ! loop through element nodes (A)
      DO ILOC=1,NLOC
        ! indexing
        INOD=NDGLNO((ELE-1)*NLOC+ILOC)
        IGLS=SONDGL((ELE-1)*NLOC+ILOC)

        ! loop through element nodes (B)
        DO JLOC=1,NLOC
          ! indexing
          JNOD=NDGLNO((ELE-1)*NLOC+JLOC)
          JGLS=SONDGL((ELE-1)*NLOC+JLOC)
          
          ! zero all matrix contributions before gaussian quadrature             
          K11=0.
          K12=0.
          K13=0.
          
          K21=0.
          K22=0.
          K23=0.
          
          K31=0.
          K32=0.
          K33=0.
          
          C11=0.
          C12=0.
          C13=0.
          
          C21=0.
          C22=0.
          C23=0.
          
          C31=0.
          C32=0.
          C33=0.
          
          MII=0.
          
          NII=0.
          
          NN=0.
          
          HL=0.
          
          ! loop through gauss points calculating integrals             
          DO GI=1,NGI
            ! stiffness matrix contributions               
            K11 = K11 &
                  +( NX(ILOC, GI)* &
                      (D11(GI)*NX(JLOC, GI) + D15(GI)*NZ(JLOC, GI) + D16(GI)*NY(JLOC, GI)) &
                  + NZ(ILOC, GI)* &
                      (D15(GI)*NX(JLOC, GI) + D55(GI)*NZ(JLOC, GI) + D56(GI)*NY(JLOC, GI)) &
                  + NY(ILOC, GI)* &
                      (D16(GI)*NX(JLOC, GI) + D56(GI)*NZ(JLOC, GI) + D66(GI)*NY(JLOC, GI)) )*DETWEI(GI)
            
            K12 = K12 &
                  +( NX(ILOC, GI)* &
                      (D12(GI)*NY(JLOC, GI) + D14(GI)*NZ(JLOC, GI) + D16(GI)*NX(JLOC, GI)) &
                  + NZ(ILOC, GI)* &
                      (D25(GI)*NY(JLOC, GI) + D45(GI)*NZ(JLOC, GI) + D56(GI)*NX(JLOC, GI)) &
                  + NY(ILOC, GI)* &
                      (D26(GI)*NY(JLOC, GI) + D46(GI)*NZ(JLOC, GI) + D66(GI)*NX(JLOC, GI)) )*DETWEI(GI)
            
            K13 = K13 &
                  +( NX(ILOC, GI)* &
                      (D13(GI)*NZ(JLOC, GI) + D14(GI)*NY(JLOC, GI) + D15(GI)*NX(JLOC, GI)) &
                  + NZ(ILOC, GI)* &
                      (D35(GI)*NZ(JLOC, GI) + D45(GI)*NY(JLOC, GI) + D55(GI)*NX(JLOC, GI)) &
                  + NY(ILOC, GI)* &
                      (D36(GI)*NZ(JLOC, GI) + D46(GI)*NY(JLOC, GI) + D56(GI)*NX(JLOC, GI)) )*DETWEI(GI)
            
            K21 = K21 &
                  +( NY(ILOC, GI)* &
                      (D12(GI)*NX(JLOC, GI) + D25(GI)*NZ(JLOC, GI) + D26(GI)*NY(JLOC, GI)) &
                  + NZ(ILOC, GI)* &
                      (D14(GI)*NX(JLOC, GI) + D45(GI)*NZ(JLOC, GI) + D46(GI)*NY(JLOC, GI)) &
                  + NX(ILOC, GI)* &
                      (D16(GI)*NX(JLOC, GI) + D56(GI)*NZ(JLOC, GI) + D66(GI)*NY(JLOC, GI)) )*DETWEI(GI)
            
            K22 = K22 &
                  +( NY(ILOC, GI)* &
                      (D22(GI)*NY(JLOC, GI) + D24(GI)*NZ(JLOC, GI) + D26(GI)*NX(JLOC, GI)) &
                  + NZ(ILOC, GI)* &
                      (D24(GI)*NY(JLOC, GI) + D44(GI)*NZ(JLOC, GI) + D46(GI)*NX(JLOC, GI)) &
                  + NX(ILOC, GI)* &
                      (D26(GI)*NY(JLOC, GI) + D46(GI)*NZ(JLOC, GI) + D66(GI)*NX(JLOC, GI)) )*DETWEI(GI)
            
            K23 = K23 &
                  +( NY(ILOC, GI)* &
                      (D23(GI)*NZ(JLOC, GI) + D24(GI)*NY(JLOC, GI) + D25(GI)*NX(JLOC, GI)) &
                  + NZ(ILOC, GI)* &
                      (D34(GI)*NZ(JLOC, GI) + D44(GI)*NY(JLOC, GI) + D45(GI)*NX(JLOC, GI)) &
                  + NX(ILOC, GI)* &
                      (D36(GI)*NZ(JLOC, GI) + D46(GI)*NY(JLOC, GI) + D56(GI)*NX(JLOC, GI)) )*DETWEI(GI)
            
            K31 = K31 &
                  +( NZ(ILOC, GI)* &
                      (D13(GI)*NX(JLOC, GI) + D35(GI)*NZ(JLOC, GI) + D36(GI)*NY(JLOC, GI)) &
                  + NY(ILOC, GI)* &
                      (D14(GI)*NX(JLOC, GI) + D45(GI)*NZ(JLOC, GI) + D46(GI)*NY(JLOC, GI)) &
                  + NX(ILOC, GI)* &
                      (D15(GI)*NX(JLOC, GI) + D55(GI)*NZ(JLOC, GI) + D56(GI)*NY(JLOC, GI)) )*DETWEI(GI)
            
            K32 = K32 &
                  +( NZ(ILOC,GI)* &
                    (D23(GI)*NY(JLOC, GI) + D35(GI)*NZ(JLOC, GI) + D36(GI)*NX(JLOC, GI)) &
                  + NY(ILOC,GI)* &
                    (D24(GI)*NY(JLOC, GI) + D44(GI)*NZ(JLOC, GI) + D46(GI)*NX(JLOC, GI)) &
                  + NX(ILOC,GI)* &
                    (D25(GI)*NY(JLOC, GI) + D45(GI)*NZ(JLOC, GI) + D56(GI)*NX(JLOC, GI)) )*DETWEI(GI)
            
            K33 = K33 &
                  +( NZ(ILOC, GI)* &
                    (D33(GI)*NZ(JLOC, GI) + D34(GI)*NY(JLOC, GI) + D35(GI)*NX(JLOC, GI)) &
                  + NY(ILOC, GI)* &
                    (D34(GI)*NZ(JLOC, GI) + D44(GI)*NY(JLOC, GI) + D45(GI)*NX(JLOC, GI)) &
                  + NX(ILOC, GI)* &
                    (D35(GI)*NZ(JLOC, GI) + D45(GI)*NY(JLOC, GI) + D55(GI)*NX(JLOC, GI)) )*DETWEI(GI)
            
              IF(KELVIN) THEN
                ! viscosity matrix contributions               
                C11 = C11 &
                      +( NX(ILOC, GI)* &
                          (F11(GI)*NX(JLOC, GI) + F15(GI)*NZ(JLOC, GI) + F16(GI)*NY(JLOC, GI)) &
                      + NZ(ILOC, GI)* &
                          (F15(GI)*NX(JLOC, GI) + F55(GI)*NZ(JLOC, GI) + F56(GI)*NY(JLOC, GI)) &
                      + NY(ILOC, GI)* &
                          (F16(GI)*NX(JLOC, GI) + F56(GI)*NZ(JLOC, GI) + F66(GI)*NY(JLOC, GI)) )*DETWEI(GI)
                
                C12 = C12 &
                      +( NX(ILOC, GI)* &
                          (F12(GI)*NY(JLOC, GI) + F14(GI)*NZ(JLOC, GI) + F16(GI)*NX(JLOC, GI)) &
                      + NZ(ILOC, GI)* &
                          (F25(GI)*NY(JLOC, GI) + F45(GI)*NZ(JLOC, GI) + F56(GI)*NX(JLOC, GI)) &
                      + NY(ILOC, GI)* &
                          (F26(GI)*NY(JLOC, GI) + F46(GI)*NZ(JLOC, GI) + F66(GI)*NX(JLOC, GI)) )*DETWEI(GI)
                
                C13 = C13 &
                      +( NX(ILOC, GI)* &
                          (F13(GI)*NZ(JLOC, GI) + F14(GI)*NY(JLOC, GI) + F15(GI)*NX(JLOC, GI)) &
                      + NZ(ILOC, GI)* &
                          (F35(GI)*NZ(JLOC, GI) + F45(GI)*NY(JLOC, GI) + F55(GI)*NX(JLOC, GI)) &
                      + NY(ILOC, GI)* &
                          (F36(GI)*NZ(JLOC, GI) + F46(GI)*NY(JLOC, GI) + F56(GI)*NX(JLOC, GI)) )*DETWEI(GI)
                
                C21 = C21 &
                      +( NY(ILOC, GI)* &
                          (F12(GI)*NX(JLOC, GI) + F25(GI)*NZ(JLOC, GI) + F26(GI)*NY(JLOC, GI)) &
                      + NZ(ILOC, GI)* &
                          (F14(GI)*NX(JLOC, GI) + F45(GI)*NZ(JLOC, GI) + F46(GI)*NY(JLOC, GI)) &
                      + NX(ILOC, GI)* &
                          (F16(GI)*NX(JLOC, GI) + F56(GI)*NZ(JLOC, GI) + F66(GI)*NY(JLOC, GI)) )*DETWEI(GI)
                
                C22 = C22 &
                      +( NY(ILOC, GI)* &
                          (F22(GI)*NY(JLOC, GI) + F24(GI)*NZ(JLOC, GI) + F26(GI)*NX(JLOC, GI)) &
                      + NZ(ILOC, GI)* &
                          (F24(GI)*NY(JLOC, GI) + F44(GI)*NZ(JLOC, GI) + F46(GI)*NX(JLOC, GI)) &
                      + NX(ILOC, GI)* &
                          (F26(GI)*NY(JLOC, GI) + F46(GI)*NZ(JLOC, GI) + F66(GI)*NX(JLOC, GI)) )*DETWEI(GI)
                
                C23 = C23 &
                      +( NY(ILOC, GI)* &
                          (F23(GI)*NZ(JLOC, GI) + F24(GI)*NY(JLOC, GI) + F25(GI)*NX(JLOC, GI)) &
                      + NZ(ILOC, GI)* &
                          (F34(GI)*NZ(JLOC, GI) + F44(GI)*NY(JLOC, GI) + F45(GI)*NX(JLOC, GI)) &
                      + NX(ILOC, GI)* &
                          (F36(GI)*NZ(JLOC, GI) + F46(GI)*NY(JLOC, GI) + F56(GI)*NX(JLOC, GI)) )*DETWEI(GI)
                
                C31 = C31 &
                      +( NZ(ILOC, GI)* &
                          (F13(GI)*NX(JLOC, GI) + F35(GI)*NZ(JLOC, GI) + F36(GI)*NY(JLOC, GI)) &
                      + NY(ILOC, GI)* &
                          (F14(GI)*NX(JLOC, GI) + F45(GI)*NZ(JLOC, GI) + F46(GI)*NY(JLOC, GI)) &
                      + NX(ILOC, GI)* &
                          (F15(GI)*NX(JLOC, GI) + F55(GI)*NZ(JLOC, GI) + F56(GI)*NY(JLOC, GI)) )*DETWEI(GI)
                
                C32 = C32 &
                      +( NZ(ILOC, GI)* &
                          (F23(GI)*NY(JLOC, GI) + F35(GI)*NZ(JLOC, GI) + F36(GI)*NX(JLOC, GI)) &
                      + NY(ILOC, GI)* &
                          (F24(GI)*NY(JLOC, GI) + F44(GI)*NZ(JLOC, GI) + F46(GI)*NX(JLOC, GI)) &
                      + NX(ILOC, GI)* &
                          (F25(GI)*NY(JLOC, GI) + F45(GI)*NZ(JLOC, GI) + F56(GI)*NX(JLOC, GI)) )*DETWEI(GI)
                
                C33 = C33 &
                      +( NZ(ILOC, GI)* &
                          (F33(GI)*NZ(JLOC, GI) + F34(GI)*NY(JLOC, GI) + F35(GI)*NX(JLOC, GI)) &
                      + NY(ILOC, GI)* &
                          (F34(GI)*NZ(JLOC, GI) + F44(GI)*NY(JLOC, GI) + F45(GI)*NX(JLOC, GI)) &
                      + NX(ILOC, GI)* &
                          (F35(GI)*NZ(JLOC, GI) + F45(GI)*NY(JLOC, GI) + F55(GI)*NX(JLOC, GI)) )*DETWEI(GI)
              ! if(kelvin)
              ENDIF
              
              ! mass contribution (incl. density)
              MII = MII + N(ILOC, GI)*DDETWE(GI)*N(JLOC, GI)
              
              ! mass contribution (excl. density)
              NN = NN + N(ILOC, GI)*N(JLOC, GI)*DETWEI(GI)
              
              ! petrov-galerkin and non-linear advection terms
              R1 =  UD(GI)*NX(ILOC, GI) &
                  + VD(GI)*NY(ILOC, GI) &
                  + WD(GI)*NZ(ILOC, GI)   
              R2 = ( UD(GI)*NX(JLOC, GI) &
                  + VD(GI)*NY(JLOC, GI) &
                  + WD(GI)*NZ(JLOC, GI) )*DETWEI(GI)
              
              UU = GAMMA(GI)*NDENGI(GI)*R1*R2
              
              HL = HL + UU
              
              ! non-linear advection contribution
              NII = NII + N(ILOC, GI)*NDENGI(GI)*R2
                
            ! gi=1,ngi
            ENDDO

            ! create r.h.s of equation
            VECX(INOD) = VECX(INOD) &
                        -( (C11 + THETA*DT*K11)*U(JNOD)+ &
                          (C12 + THETA*DT*K12)*V(JNOD)+ &
                          (C13 + THETA*DT*K13)*W(JNOD) ) &
  !                       -( K11*XDISP(JNOD)+ &
  !                          K12*YDISP(JNOD)+ &
  !                          K13*ZDISP(JNOD) ) &
                        - NII*U(JNOD) &
                        - HL*U(JNOD) &
                        + NN*SOURCX(JNOD)
            VECY(INOD) = VECY(INOD) &
                        -( (C21 + THETA*DT*K21)*U(JNOD)+ &
                          (C22 + THETA*DT*K22)*V(JNOD)+ &
                          (C23 + THETA*DT*K23)*W(JNOD) ) &
  !                       -( K21*XDISP(JNOD)+ &
  !                          K22*YDISP(JNOD)+ &
  !                          K23*ZDISP(JNOD) ) &
                        - NII*V(JNOD) &
                        - HL*V(JNOD) &
                        + NN*SOURCY(JNOD)
            VECZ(INOD) = VECZ(INOD) &
                        -( (C31 + THETA*DT*K31)*U(JNOD)+ &
                          (C32 + THETA*DT*K32)*V(JNOD)+ &
                          (C33 + THETA*DT*K33)*W(JNOD) ) &
  !                       -( K31*XDISP(JNOD)+ &
  !                          K32*YDISP(JNOD)+ &
  !                          K33*ZDISP(JNOD) ) &
                        - NII*W(JNOD) &
                        - HL*W(JNOD) &
                        + NN*SOURCZ(JNOD)
      
            VECX(INOD) = VECX(INOD) &
                        - (K11*XDISP(JNOD)+ &
                          K12*YDISP(JNOD)+ &
                          K13*ZDISP(JNOD))
            VECY(INOD) = VECY(INOD) &
                        - (K21*XDISP(JNOD)+ &
                          K22*YDISP(JNOD)+ &
                          K23*ZDISP(JNOD))
            VECZ(INOD) = VECZ(INOD) &
                        - (K31*XDISP(JNOD)+ &
                          K32*YDISP(JNOD)+ &
                          K33*ZDISP(JNOD))
      
            ! find position in matrix
            DO COUNT = FINDRM(INOD), FINDRM(INOD+1)-1
              IF(JNOD.EQ.COLM(COUNT)) POSMAT=COUNT
            ENDDO
            
            ICENT = CENTRM(INOD)
            
            ! lump the mass matrix
            ML(INOD) = ML(INOD) + MII
            
            BIGM(ICENT + IBL11) = BIGM(ICENT + IBL11) + MII*RLUMP
            BIGM(ICENT + IBL22) = BIGM(ICENT + IBL22) + MII*RLUMP
            BIGM(ICENT + IBL33) = BIGM(ICENT + IBL33) + MII*RLUMP
            
            ! place other terms in matrix
            BIGM(POSMAT + IBL11) = BIGM(POSMAT + IBL11) &
                                  + THETA*THETA*DT*DT*K11 &
                                  + THETA*DT*C11 &
                                  + MII*(1. - RLUMP) & 
                                  + THETA*DT*NII &
                                  + THETA*DT*HL
            BIGM(POSMAT + IBL12) = BIGM(POSMAT + IBL12) &
                                  + THETA*THETA*DT*DT*K12 & 
                                  + THETA*DT*C12
            BIGM(POSMAT + IBL13) = BIGM(POSMAT + IBL13) &
                                  + THETA*THETA*DT*DT*K13 & 
                                  + THETA*DT*C13
            
            BIGM(POSMAT + IBL21) = BIGM(POSMAT + IBL21) &
                                  + THETA*THETA*DT*DT*K21 & 
                                  + THETA*DT*C21
            BIGM(POSMAT + IBL22) = BIGM(POSMAT + IBL22) &
                                  + THETA*THETA*DT*DT*K22 &
                                  + THETA*DT*C22 &
                                  + MII*(1. - RLUMP) &
                                  + THETA*DT*NII &
                                  + THETA*DT*HL
            BIGM(POSMAT + IBL23) = BIGM(POSMAT + IBL23) &
                                  + THETA*THETA*DT*DT*K23 & 
                                  + THETA*DT*C23
            
            BIGM(POSMAT + IBL31) = BIGM(POSMAT + IBL31) &
                                  + THETA*THETA*DT*DT*K31 &
                                  + THETA*DT*C31
            BIGM(POSMAT + IBL32) = BIGM(POSMAT + IBL32) &
                                  + THETA*THETA*DT*DT*K32 & 
                                  + THETA*DT*C32
            BIGM(POSMAT + IBL33) = BIGM(POSMAT + IBL33) &
                                  + THETA*THETA*DT*DT*K33  &
                                  + THETA*DT*C33  &
                                  + MII*(1. - RLUMP) &
                                  + THETA*DT*NII &
                                  + THETA*DT*HL
  
          ! jloc=1,nloc
          ENDDO
  
        ! zero before gaussian quadrature
        SPCONX = 0.
        SPCONY = 0.
        SPCONZ = 0.
                        
        IF(BINGHAM) THEN
                        
          ! loop through gauss points evaluating the integral of
          ! the deviatoric stresses at the gauss points times
          ! the derivatives of the shape functions
          DO GI=1,NGI
            SGI = (ELE-1)*NGI + GI
            
            BD11GI = D11(GI)*NX(ILOC,GI) + D15(GI)*NZ(ILOC,GI) + D16(GI)*NY(ILOC,GI)
            BD12GI = D12(GI)*NX(ILOC,GI) + D25(GI)*NZ(ILOC,GI) + D26(GI)*NY(ILOC,GI)
            BD13GI = D13(GI)*NX(ILOC,GI) + D35(GI)*NZ(ILOC,GI) + D36(GI)*NY(ILOC,GI)
            BD14GI = D14(GI)*NX(ILOC,GI) + D45(GI)*NZ(ILOC,GI) + D46(GI)*NY(ILOC,GI)
            BD15GI = D15(GI)*NX(ILOC,GI) + D55(GI)*NZ(ILOC,GI) + D56(GI)*NY(ILOC,GI)
            BD16GI = D16(GI)*NX(ILOC,GI) + D56(GI)*NZ(ILOC,GI) + D66(GI)*NY(ILOC,GI)
            
            BD21GI = D12(GI)*NY(ILOC,GI) + D14(GI)*NZ(ILOC,GI) + D16(GI)*NX(ILOC,GI)
            BD22GI = D22(GI)*NY(ILOC,GI) + D24(GI)*NZ(ILOC,GI) + D26(GI)*NX(ILOC,GI)
            BD23GI = D23(GI)*NY(ILOC,GI) + D34(GI)*NZ(ILOC,GI) + D36(GI)*NX(ILOC,GI)
            BD24GI = D24(GI)*NY(ILOC,GI) + D44(GI)*NZ(ILOC,GI) + D46(GI)*NX(ILOC,GI)
            BD25GI = D25(GI)*NY(ILOC,GI) + D45(GI)*NZ(ILOC,GI) + D56(GI)*NX(ILOC,GI)
            BD26GI = D26(GI)*NY(ILOC,GI) + D46(GI)*NZ(ILOC,GI) + D66(GI)*NX(ILOC,GI)
            
            BD31GI = D13(GI)*NZ(ILOC,GI) + D14(GI)*NY(ILOC,GI) + D15(GI)*NX(ILOC,GI)
            BD32GI = D23(GI)*NZ(ILOC,GI) + D24(GI)*NY(ILOC,GI) + D25(GI)*NX(ILOC,GI)
            BD33GI = D33(GI)*NZ(ILOC,GI) + D34(GI)*NY(ILOC,GI) + D35(GI)*NX(ILOC,GI)
            BD34GI = D34(GI)*NZ(ILOC,GI) + D44(GI)*NY(ILOC,GI) + D45(GI)*NX(ILOC,GI)
            BD35GI = D35(GI)*NZ(ILOC,GI) + D45(GI)*NY(ILOC,GI) + D55(GI)*NX(ILOC,GI)
            BD36GI = D36(GI)*NZ(ILOC,GI) + D46(GI)*NY(ILOC,GI) + D56(GI)*NX(ILOC,GI)
  
            IF(BINGHAM) THEN
              SPCONX = SPCONX &
                      +( BD11GI*STVPXX(SGI) &
                      + BD12GI*STVPYY(SGI) &
                      + BD13GI*STVPZZ(SGI) &
                      + BD14GI*STVPYZ(SGI) &
                      + BD15GI*STVPXZ(SGI) &
                      + BD16GI*STVPXY(SGI) )*DETWEI(GI)
              SPCONY = SPCONY &
                      +( BD21GI*STVPXX(SGI) &
                      + BD22GI*STVPYY(SGI) &
                      + BD23GI*STVPZZ(SGI) &
                      + BD24GI*STVPZZ(SGI) &
                      + BD25GI*STVPXZ(SGI) &
                      + BD26GI*STVPXY(SGI) )*DETWEI(GI)
              SPCONZ = SPCONZ &
                      +( BD31GI*STVPXX(SGI) &
                      + BD32GI*STVPYY(SGI) &
                      + BD33GI*STVPZZ(SGI) &
                      + BD34GI*STVPYZ(SGI) &
                      + BD35GI*STVPXZ(SGI) &
                      + BD36GI*STVPXY(SGI) )*DETWEI(GI)
            ENDIF
  
          ! gi=1,ngi
          ENDDO
  
        ! if(bingham)
        ENDIF
            
        VECX(INOD) = VECX(INOD) + SPCONX
        VECY(INOD) = VECY(INOD) + SPCONY
        VECZ(INOD) = VECZ(INOD) + SPCONZ
  
      ! iloc=1,nloc
      ENDDO
    
    ! ele=1,totele
    ENDDO
        
    ewrite_minmax(BIGM)
    ewrite_minmax(ML)
    ewrite_minmax(VECX)
    ewrite_minmax(VECY)
    ewrite_minmax(VECZ)  
    if(have_gravity) then
      deallocate(sourcx)
      deallocate(sourcy)
      deallocate(sourcz)
    end if
          
    ewrite(3,*)  'EXITING SOLID3D()'
    RETURN
  END SUBROUTINE SOLID3D
  !************************************************************************
          
  

  !************************************************************************
  SUBROUTINE STEVAL( &
                    NLOC,NGI,MLOC, &
                    STVPXX,STVPYY,STVPZZ, &
                    STVPYZ,STVPXZ,STVPXY, &
                    F, FOLD, &
                    SOLIDS, &
                    OFFYIELDENV &
                    ,state)
    ! =============================================================
    ! Subroutine to check whether the current iterate's result lies inside the
    ! yield envelope (if the material is Bingham).  If it fails this criterion
    ! a plastic integration step is used to return the stress state to the yield
    ! surface F = 0.
    ! =============================================================
    ! =============================================================
    ! inputs
    ! =============================================================
    type(state_type), intent(inout) :: state

    INTEGER, INTENT(IN) :: NLOC, NGI
    ! nloc  = number of displacement nodes per element
    ! ngi   = number of gauss points (for current element)
        
    INTEGER, INTENT(IN) :: MLOC
    ! mloc = number of pressure nodes per element
    
    INTEGER, INTENT(IN) :: SOLIDS
    ! solids = dictates behaviour of material:
    !          1: elastic... this sub does nothing
    !          2: kelvin...  this sub does nothing
    !          3: bingham... this sub corrects plastic overstep
    
    ! =============================================================
    ! outputs
    ! =============================================================
  
    ! plastic strains:
    REAL, INTENT(INOUT) :: STVPXX(:), STVPYY(:), STVPZZ(:)
    ! stvpxx(totele*ngi) = plastic strain in the x direction 
    !                      on a plane perpendicular to the x axis
    ! stvpyy(totele*ngi) = plastic strain in the y direction 
    !                      on a plane perpendicular to the y axis
    ! stvpzz(totele*ngi) = plastic strain in the z direction 
    !                      on a plane perpendicular to the z axis
    
    REAL, INTENT(INOUT) :: STVPYZ(:), STVPXZ(:), STVPXY(:)
    ! stvpyz(totele*ngi) = plastic strain in the z direction 
    !                      on a plane perpendicular to the y axis
    ! stvpxz(totele*ngi) = plastic strain in the z direction 
    !                      on a plane perpendicular to the x axis
    ! stvpxy(totele*ngi) = plastic strain in the y direction 
    !                      on a plane perpendicular to the x axis
    
    ! no real need for these to be passed in and out - except for gid?:
    REAL, INTENT(INOUT) :: F(:), FOLD(:)
    ! f(totele*ngi)    = compliance to yield surface after plastic correction
    ! fold(totele*ngi) = compliance to yield surface before plastic correction
  
    LOGICAL, INTENT(INOUT) :: OFFYIELDENV
    ! offyieldenv = indicates whether the timestep requires another iteration 
    !               to achieve the specified tolerance tol
        
    ! =============================================================
    ! local variables
    ! =============================================================
  
    ! fractions/surds/parameters etc:
    REAL :: TWOTHI, FOURTHI, SQRTHR
    REAL :: ONESIX, ONETHI, ONEHAL
    PARAMETER(TWOTHI = 2./3., FOURTHI = 4./3.)
    PARAMETER(ONESIX = 1./6., ONETHI = 1./3., ONEHAL = 1./2.)
    REAL :: PI
    PARAMETER(PI=3.1415926538)
    ! end of parameters etc.
    
    ! indexing integers:
    INTEGER :: IGLS, IGLV, IGLX, IGLP, SGI, K
    ! end of indexing integers
    
    REAL :: DGAMMA, D2GAMMA, DNORM
    ! dgamma  = plastic pseudoviscosity
    ! d2gamma = change in dgamma
    ! dnorm   = denominator in plastic corrector step
    
    REAL :: MUGIXX
    ! mugixx = muptxx(snonod) evaluated at a gauss pt
    
    INTEGER :: ELE, ILOC, GI
    ! ele  = integer counter through the elements
    ! iloc = integer counter through the nodes in each element (a)
    ! gi   = integer counter through the gauss points in each element
    
    REAL :: AGI,BGI,CGI
    REAL :: DGI,EGI,FGI
    REAL :: GGI,HGI,KGI
    ! [agi, bgi, cgi;
    !  dgi, egi, fgi;
    !  ggi, hgi, kgi] = entries of jacobi matrix
    
    REAL :: DETJ,INVDET
    ! detj   = determinant of jacobi matrix (jacobian)
    ! invdet = 1./detj 
    
    REAL :: A11,A12,A13
    REAL :: A21,A22,A23
    REAL :: A31,A32,A33
    ! [a11, a12, a13;
    !  a21, a22, a23;
    !  a31, a32, a33] = entries of inverse jacobi matrix
    
    REAL :: NX(NLOC,NGI), NY(NLOC,NGI), NZ(NLOC,NGI)
    ! nx(nloc,ngi) = derivative of displacement shape function in global x direction
    !                (evaluated at a gauss point around a node)
    ! ny(nloc,ngi) = derivative of displacement shape function in global y direction
    !                (evaluated at a gauss point around a node)
    ! nz(nloc,ngi) = derivative of displacement shape function in global z direction
    !                (evaluated at a gauss point around a node)
    
    REAL :: D11(NGI), D12(NGI), D13(NGI), D14(NGI), D15(NGI), D16(NGI)
    REAL :: D22(NGI), D23(NGI), D24(NGI), D25(NGI), D26(NGI)
    REAL :: D33(NGI), D34(NGI), D35(NGI), D36(NGI)
    REAL :: D44(NGI), D45(NGI), D46(NGI)
    REAL :: D55(NGI), D56(NGI)
    REAL :: D66(NGI)
    ! d11/12/13/14/15/16/22/23/24/25/26/33/34/35/36/44/45/46/55/56/66
    !   = components of the symmetric matrix d relating stress to deviatoric strain
    !   - for a linear isotropic 3-d solid the matrix d will have the form:
    !    4./3.*mu_e     -2./3.mu_e     -2./3.mu_e    0.    0.   0.
    !    -2./3.mu_e     4./3.*mu_e     -2./3.mu_e    0.    0.   0.
    !    -2./3.mu_e     -2./3.mu_e     4./3.*mu_e    0.    0.   0.
    !        0.             0.            0.        mu_e   0.   0.
    !        0.             0.            0.         0.   mu_e  0.
    !        0.             0.            0.         0.    0.  mu_e
    
    REAL :: TWOTHIMU, FOURTHIMU
    ! twothimu  = -2./3.*mu_(e/v)
    ! fourthimu = 4./3.*mu_(e/v)
    
    ! parameters to set the model type:     
    ! nb: default behaviour is for both to be false, when the model simulates a purely
    !     linearly elastic material
    LOGICAL :: KELVIN, BINGHAM
    ! kelvin  = simulates a kelvin substance with elastic and viscous components
    !           in parallel
    !           this option needs to be set to true for multi-material modelling where
    !           one section of the mass fraction is a solid and the other a liquid
    ! bingham = simulates a bingham substance which initially behaves elastically until
    !           a yield stress is reached at which point it begins to flow
    !           nb: think this only works for lagrangian (one material) models at the moment
    
    REAL :: SINFAN, COSFAN, NDENOM
    ! sinfan = sin(frcang)
    ! cosfan = cos(frcang)
    ! ndenom = denominator for evaluation of stress invariants
    
    ! pressure:
    REAL :: PGI(NGI), INVPGI, DENGI(NGI)
    ! p(fredop*npress) = pressure
    ! pgi(ngi)         = pressure evaluated at the gauss points
    ! invpgi           = 1./pgi(gi)
    ! dengi(ngi)       = density (to be) evaluated at gauss points
    
    REAL :: NXDGIX, NYDGIX, NZDGIX
    ! nxdgix(ngi) = derivative w.r.t x of the previous iterate's x displacement 
    !               evaluated at the gauss points
    ! nydgix(ngi) = derivative w.r.t x of the previous iterate's y displacement 
    !               evaluated at the gauss points
    ! nzdgix(ngi) = derivative w.r.t x of the previous iterate's z displacement 
    !               evaluated at the gauss points
    
    REAL :: NXDGIY, NYDGIY, NZDGIY
    ! nxdgiy(ngi) = derivative w.r.t y of the previous iterate's x displacement 
    !               evaluated at the gauss points
    ! nydgiy(ngi) = derivative w.r.t y of the previous iterate's y displacement 
    !               evaluated at the gauss points
    ! nzdgiy(ngi) = derivative w.r.t y of the previous iterate's z displacement 
    !               evaluated at the gauss points
    
    REAL :: NXDGIZ, NYDGIZ, NZDGIZ
    ! nxdgiz(ngi) = derivative w.r.t z of the previous iterate's x displacement 
    !               evaluated at the gauss points
    ! nydgiz(ngi) = derivative w.r.t z of the previous iterate's y displacement 
    !               evaluated at the gauss points
    ! nzdgiz(ngi) = derivative w.r.t z of the previous iterate's z displacement 
    !               evaluated at the gauss points
    
    REAL :: TAXXGI, TAYYGI, TAZZGI
    ! tauxxgi = deviatoric stress evaluated at the gauss points
    !           acting in x direction on planes perpendicular to the x axis
    ! tauyygi = deviatoric stress evaluated at the gauss points
    !           acting in y direction on planes perpendicular to the y axis
    ! tauzzgi = deviatoric stress evaluated at the gauss points
    !           acting in z direction on planes perpendicular to the z axis
    
    REAL :: TAXYGI, TAXZGI, TAYZGI
    ! tauxygi = deviatoric stress evaluated at the gauss points
    !           acting in y direction on planes perpendicular to the x axis
    ! tauxzgi = deviatoric stress evaluated at the gauss points
    !           acting in z direction on planes perpendicular to the x axis
    ! tauyzgi = deviatoric stress evaluated at the gauss points
    !           acting in z direction on planes perpendicular to the y axis
    
    REAL :: SRXXGI, SRYYGI, SRZZGI
    ! srxxgi = total stress evaluated at the gauss points = tauxxgi - p(gi)
    !          acting in x direction on planes perpendicular to the x axis
    ! sryygi = total stress evaluated at the gauss points = tauyygi - p(gi)
    !          acting in y direction on planes perpendicular to the y axis
    ! srzzgi = total stress evaluated at the gauss points = tauzzgi - p(gi)
    !          acting in z direction on planes perpendicular to the z axis
    
    ! q = 0: flow surface off which normality is enforced
    REAL :: DQDSXX, DQDSYY, DQDSZZ
    ! dqdsxx = derivative of q wrt xx component of total stress
    ! dqdsyy = derivative of q wrt xx component of total stress
    ! dqdszz = derivative of q wrt xx component of total stress
    
    REAL :: DQDSYZ, DQDSXZ, DQDSXY
    ! dqdsyz = derivative of q wrt yz component of total stress
    ! dqdsxz = derivative of q wrt xz component of total stress
    ! dqdsxy = derivative of q wrt xy component of total stress
    
    REAL :: DQDP, DQDJ2, DQDJ3
    ! dqdp  = derivative of q wrt pressure
    ! dqdj2 = derivative of q wrt j2 stress invariant
    ! dqdj3 = derivative of q wrt j3 stress invariant
    
    REAL :: D2QDP2, D2QDJ22, D2QDJ32
    ! d2qdp2  = second derivative of q wrt pressure - not used yet
    ! d2qdj22 = second derivative of q wrt j2 stress invariant - not used yet
    ! d2qdj32 = second derivative of q wrt j3 stress invariant - not used yet
    
    REAL :: D2QDPDJ2, D2QDPDJ3, D2QDJ2DJ3
    ! d2qdpdj2  = cross derivative of q wrt pressure and j2 stress invariant - not used yet
    ! d2qdpdj3  = cross derivative of q wrt pressure and j3 stress invariant - not used yet
    ! d2qdj2dj3 = cross derivative of q wrt j2 and j3 stress invariants - not used yet
    
    ! f = 0: yield surface
    REAL :: DFDSXX, DFDSYY, DFDSZZ
    ! dqdsxx = derivative of f wrt xx component of total stress
    ! dqdsyy = derivative of f wrt xx component of total stress
    ! dqdszz = derivative of f wrt xx component of total stress
    
    REAL :: DFDSYZ, DFDSXZ, DFDSXY
    ! dqdsyz = derivative of f wrt yz component of total stress
    ! dqdsxz = derivative of f wrt xz component of total stress
    ! dqdsxy = derivative of f wrt xy component of total stress
    
    REAL :: J2, J3, LODETH, SQRTJ2, DILANG
    ! j2     = j2 invariant of stress
    ! j3     = j3 invariant of stress
    ! lodeth = lode angle invariant of stress - not used yet
    ! sqrtj2 = sqrt(j2)
    ! dilang = dilation angle - not used yet
    
    REAL :: DFDP, DFDJ2, DFDJ3
    ! dfdp  = derivative of f wrt pressure
    ! dfdj2 = derivative of f wrt j2 stress invariant
    ! dfdj3 = derivative of f wrt j3 stress invariant
    
    REAL :: D2FDP2, D2FDJ22, D2FDJ32
    ! d2fdp2  = second derivative of f wrt pressure - not used yet
    ! d2fdj22 = second derivative of f wrt j2 stress invariant - not used yet
    ! d2fdj32 = second derivative of f wrt j3 stress invariant - not used yet
    
    REAL :: D2FDPDJ2, D2FDPDJ3, D2FDJ2DJ3
    ! d2fdpdj2  = cross derivative of f wrt pressure and j2 stress invariant - not used yet
    ! d2fdpdj3  = cross derivative of f wrt pressure and j3 stress invariant - not used yet
    ! d2fdj2dj3 = cross derivative of f wrt j2 and j3 stress invariants - not used yet
    
    REAL :: DPDS11,DPDS12,DPDS13,DPDS14,DPDS15,DPDS16
    REAL :: DPDS22,DPDS23,DPDS24,DPDS25,DPDS26
    REAL :: DPDS33,DPDS34,DPDS35,DPDS36
    REAL :: DPDS44,DPDS45,DPDS46
    REAL :: DPDS55,DPDS56
    REAL :: DPDS66
    ! dpds11/12/13/14/15/16/22/23/24/25/26/33/34/35/36/44/45/46/55/56/66
    !   = components of symmetric matrix that once multiplied with stress
    !     gives the 6 components of the derivative of pressure wrt stress
    
    REAL :: D2PDS211,D2PDS212,D2PDS213,D2PDS214,D2PDS215,D2PDS216
    REAL :: D2PDS222,D2PDS223,D2PDS224,D2PDS225,D2PDS226
    REAL :: D2PDS233,D2PDS234,D2PDS235,D2PDS236
    REAL :: D2PDS244,D2PDS245,D2PDS246
    REAL :: D2PDS255,D2PDS256
    REAL :: D2PDS266
    ! d2pds211/12/13/14/15/16/22/23/24/25/26/33/34/35/36/44/45/46/55/56/66
    !   = components of symmetric matrix that once multiplied with stress
    !     gives the 6 components of the second derivative of pressure wrt stress
    !   - not used yet
    
    REAL :: DJ2D11,DJ2D12,DJ2D13,DJ2D14,DJ2D15,DJ2D16
    REAL :: DJ2D22,DJ2D23,DJ2D24,DJ2D25,DJ2D26
    REAL :: DJ2D33,DJ2D34,DJ2D35,DJ2D36
    REAL :: DJ2D44,DJ2D45,DJ2D46
    REAL :: DJ2D55,DJ2D56
    REAL :: DJ2D66
    ! dj2d11/12/13/14/15/16/22/23/24/25/26/33/34/35/36/44/45/46/55/56/66
    !   = components of symmetric matrix that once multiplied with stress
    !     gives the 6 components of the derivative of the j2 stress invariant wrt stress
    
    REAL :: D2J2D211,D2J2D212,D2J2D213,D2J2D214,D2J2D215,D2J2D216
    REAL :: D2J2D222,D2J2D223,D2J2D224,D2J2D225,D2J2D226
    REAL :: D2J2D233,D2J2D234,D2J2D235,D2J2D236
    REAL :: D2J2D244,D2J2D245,D2J2D246
    REAL :: D2J2D255,D2J2D256
    REAL :: D2J2D266
    ! d2j2d211/12/13/14/15/16/22/23/24/25/26/33/34/35/36/44/45/46/55/56/66
    !   = components of symmetric matrix that once multiplied with stress
    !     gives the 6 components of the second derivative of the j2 stress invariant wrt stress
    !   - not used yet
    
    REAL :: DJ3D11,DJ3D12,DJ3D13,DJ3D14,DJ3D15,DJ3D16
    REAL :: DJ3D22,DJ3D23,DJ3D24,DJ3D25,DJ3D26
    REAL :: DJ3D33,DJ3D34,DJ3D35,DJ3D36
    REAL :: DJ3D44,DJ3D45,DJ3D46
    REAL :: DJ3D55,DJ3D56
    REAL :: DJ3D66
    ! dj3d11/12/13/14/15/16/22/23/24/25/26/33/34/35/36/44/45/46/55/56/66
    !   = components of symmetric matrix that once multiplied with stress
    !     gives the 6 components of the derivative of the j3 stress invariant wrt stress
    
    REAL :: D2J3D211,D2J3D212,D2J3D213,D2J3D214,D2J3D215,D2J3D216
    REAL :: D2J3D222,D2J3D223,D2J3D224,D2J3D225,D2J3D226
    REAL :: D2J3D233,D2J3D234,D2J3D235,D2J3D236
    REAL :: D2J3D244,D2J3D245,D2J3D246
    REAL :: D2J3D255,D2J3D256
    REAL :: D2J3D266
    ! d2j3d211/12/13/14/15/16/22/23/24/25/26/33/34/35/36/44/45/46/55/56/66
    !   = components of symmetric matrix that once multiplied with stress
    !     gives the 6 components of the second derivative of the j3 stress invariant wrt stress
    !   - not used yet
    
    REAL :: TOL
    PARAMETER (TOL = 1.E-6)
    ! tol = parameter indicating the desired tolerance within which 
    !       the stress state should lie relative to the yield surface f = 0

    type(vector_field), pointer :: displacement
    type(scalar_field) :: displacementx, displacementy, displacementz
    real, dimension(:), pointer :: nxdisp, nydisp, nzdisp
  
    type(vector_field), pointer :: coordinate
    type(scalar_field) :: coordinatex, coordinatey, coordinatez
    real, dimension(:), pointer :: x, y, z
  
    type(tensor_field), pointer :: elasticity
    type(scalar_field) :: elasticityxx
    real, dimension(:), pointer :: muptxx
  
    type(scalar_field), pointer :: density
    real, dimension(:), pointer :: denpt

    type(scalar_field), pointer :: pressure
    real, dimension(:), pointer :: p

    type(scalar_field), pointer :: cohesion, frictionangle
    ! bingham material paramters:
    real :: frcang, cohesi
    ! frcang = angle of internal friction (radians)
    ! cohesi = cohesion of material
    
    integer :: totele

    integer, dimension(:), pointer :: sondgl, vondgl, xondgl, pndgln

    REAL :: N(NLOC,NGI), NLX(NLOC,NGI), NLY(NLOC,NGI), NLZ(NLOC,NGI)
    ! n(nloc,ngi)   = displacement shape function 
    ! nlx(nloc,ngi) = derivative of displacement shape function in local x direction
    !                 (evaluated at a gauss point around a node)
    ! nly(nloc,ngi) = derivative of displacement shape function in local y direction
    !                 (evaluated at a gauss point around a node)
    ! nlz(nloc,ngi) = derivative of displacement shape function in local z direction
    !                 (evaluated at a gauss point around a node)
    
    REAL :: M(MLOC,NGI), MLX(MLOC,NGI), MLY(MLOC,NGI), MLZ(MLOC,NGI)
    ! m(mloc,ngi)   = pressure shape function
    !                 (evaluated at a gauss point around a node)

    ! =============================================================
    ! This subroutine works in two halves.  The first evaluates the stress
    ! state from the previous iterate and determines if it lies within the
    ! yield envelope or not.  The second only occurs if it does not, in which
    ! case a local (at a single Gauss point) integration loop is performed to
    ! correct the overstep back to the yield surface.
    ! 
    ! called from fluids.f
    ! 
    ! Description                                   Programmer      Date
    ! ==================================================================
    ! Original version............................. CRGW        2006-06-02
    ! =============================================================
  
    ewrite(3,*) 'IN STEVAL()'

    displacement=>extract_vector_field(state, "Displacement")

    displacementx=extract_scalar_field_from_vector_field(displacement, 1)
    nxdisp=>displacementx%val

    displacementy=extract_scalar_field_from_vector_field(displacement, 2)
    nydisp=>displacementy%val

    displacementz=extract_scalar_field_from_vector_field(displacement, 3)
    nzdisp=>displacementz%val


    coordinate=>extract_vector_field(state, "Coordinate")

    coordinatex=extract_scalar_field_from_vector_field(coordinate, 1)
    x=>coordinatex%val

    coordinatey=extract_scalar_field_from_vector_field(coordinate, 2)
    y=>coordinatey%val

    coordinatez=extract_scalar_field_from_vector_field(coordinate, 3)
    z=>coordinatez%val


    elasticity=>extract_tensor_field(state, "Elasticity")

    elasticityxx=extract_scalar_field_from_tensor_field(elasticity, 1, 1)
    muptxx=>elasticityxx%val

    density=>extract_scalar_field(state, "Density")
    denpt=>density%val

    pressure=>extract_scalar_field(state, "Pressure")
    p=>pressure%val

    cohesion => extract_scalar_field(state,'VelocityCohesion')
    cohesi = maxval(cohesion%val)

    frictionangle => extract_scalar_field(state,'VelocityFrictionAngle')
    frcang = maxval(frictionangle%val)

    totele = element_count(displacement)
  
    sondgl => displacement%mesh%ndglno
    vondgl => displacement%mesh%ndglno
    xondgl => coordinate%mesh%ndglno
    pndgln => pressure%mesh%ndglno

    call extract_old_element(ele_shape(displacement,1), N, NLX, NLY, NLZ)
    call extract_old_element(ele_shape(pressure,1), M, MLX, MLY, MLZ)

    F = 0.
    FOLD = 0.
  
    OFFYIELDENV  = .FALSE.
  
    KELVIN  = .FALSE.
    BINGHAM = .FALSE.
    IF(SOLIDS.EQ.2) KELVIN  = .TRUE.
    IF(SOLIDS.EQ.3) BINGHAM = .TRUE.
    
    ewrite(3,*) 'SOLIDS, KELVIN, BINGHAM, COHESI, FRCANG = ', &
            SOLIDS, KELVIN, BINGHAM, COHESI, FRCANG
    
    IF(.NOT.BINGHAM) RETURN
        
    ! yield and flow surface parameters       
    SQRTHR = SQRT(3.)               ! square root of 3
    DILANG = 0.                     ! dilation angle - not used yet
    
    SINFAN = SIN(FRCANG)            ! sine of friction angle
    COSFAN = COS(FRCANG)            ! cosine of friction angle
    NDENOM = 1./(3.-SINFAN)         ! denominator of invariant calculations
          
    ! flow surface first derivatives
    DQDP  = 0.
    DQDJ2 = SQRTHR
    DQDJ3 = 0.
    
    ! flow surface second derivatives
    D2QDP2  = 0.
    D2QDJ22 = 0.
    D2QDJ32 = 0.
    
    ! flow surface cross derivatives
    D2QDPDJ2  = 0.
    D2QDPDJ3  = 0.
    D2QDJ2DJ3 = 0.
    
    ! yield surface first derivatives
    DFDP  = -6.*SINFAN*NDENOM
    DFDJ2 = SQRTHR
    DFDJ3 = 0.
    
    ! yield surface second derivatives
    D2FDP2  = 0.
    D2FDJ22 = 0.
    D2FDJ32 = 0.
          
    ! yield surface cross derivatives
    D2FDPDJ2  = 0.
    D2FDPDJ3  = 0.
    D2FDJ2DJ3 = 0.
          
          
    ! loop through elements
    DO ELE=1,TOTELE
      
      ! loop through gauss points evaluating parameters at them
      DO GI=1,NGI
  
        SGI = (ELE-1)*NGI + GI
        
        ! zero the jacobian components         
        AGI=0.
        BGI=0.
        CGI=0.
        
        DGI=0.
        EGI=0.
        FGI=0.
        
        GGI=0.
        HGI=0.
        KGI=0.
        
        ! zero previous iterate's differentiated displacements
        NXDGIX=0.
        NXDGIY=0.
        NXDGIZ=0.
        
        NYDGIX=0.
        NYDGIY=0.
        NYDGIZ=0.
        
        NZDGIX=0.
        NZDGIY=0.
        NZDGIZ=0.
        
        PGI(GI)=0.
        DENGI(GI)=0.
                              
        ! zero temporary material parameters (prior to summation)
        MUGIXX=0.
  
        ! loop through element nodes
        DO ILOC=1,NLOC
          ! indexing
          IGLS=SONDGL((ELE-1)*NLOC+ILOC)
          IGLX=XONDGL((ELE-1)*NLOC+ILOC)
          
          ! sum the components of the jacobian           
          AGI=AGI+NLX(ILOC,GI)*X(IGLX) 
          BGI=BGI+NLX(ILOC,GI)*Y(IGLX) 
          CGI=CGI+NLX(ILOC,GI)*Z(IGLX) 
          
          DGI=DGI+NLY(ILOC,GI)*X(IGLX) 
          EGI=EGI+NLY(ILOC,GI)*Y(IGLX) 
          FGI=FGI+NLY(ILOC,GI)*Z(IGLX) 
          
          GGI=GGI+NLZ(ILOC,GI)*X(IGLX) 
          HGI=HGI+NLZ(ILOC,GI)*Y(IGLX) 
          KGI=KGI+NLZ(ILOC,GI)*Z(IGLX) 
          
          ! sum parameters at nodes to get the value at the gauss point
          MUGIXX=MUGIXX+N(ILOC,GI)*MUPTXX(IGLS)
          
          ! density evaluated at the gauss points
          DENGI(GI)  = DENGI(GI) +N(ILOC,GI)*DENPT(IGLS)
          
        ! iloc=1,nloc
        ENDDO
  
        ! determinant of jacobi matrix (determinant)          
        DETJ = AGI*(EGI*KGI - FGI*HGI) - &
              BGI*(DGI*KGI - FGI*GGI) + &
              CGI*(DGI*HGI - EGI*GGI)
        
        INVDET    =1./DETJ
        
        ! calculate the components of the jacobi inverse
        A11= (EGI*KGI - FGI*HGI)*INVDET
        A21=-(DGI*KGI - FGI*GGI)*INVDET
        A31= (DGI*HGI - EGI*GGI)*INVDET
        
        A12=-(BGI*KGI - CGI*HGI)*INVDET
        A22= (AGI*KGI - CGI*GGI)*INVDET
        A32=-(AGI*HGI - BGI*GGI)*INVDET
        
        A13= (BGI*FGI - CGI*EGI)*INVDET
        A23=-(AGI*FGI - CGI*DGI)*INVDET
        A33= (AGI*EGI - BGI*DGI)*INVDET
      
        ! evaluate the derivatives of the velocity shape functions at the gauss points
        DO ILOC=1,NLOC
          NX(ILOC,GI)=A11*NLX(ILOC,GI)+A12*NLY(ILOC,GI)+A13*NLZ(ILOC,GI)
          NY(ILOC,GI)=A21*NLX(ILOC,GI)+A22*NLY(ILOC,GI)+A23*NLZ(ILOC,GI)
          NZ(ILOC,GI)=A31*NLX(ILOC,GI)+A32*NLY(ILOC,GI)+A33*NLZ(ILOC,GI)
        ! iloc=1,nloc
        ENDDO
      
        FOURTHIMU = FOURTHI*MUGIXX 
        TWOTHIMU  = -TWOTHI*MUGIXX
        
        D11(GI)=FOURTHIMU
        D12(GI)=TWOTHIMU
        D13(GI)=TWOTHIMU
        D14(GI)=0.
        D15(GI)=0.
        D16(GI)=0.
        
        D22(GI)=FOURTHIMU
        D23(GI)=TWOTHIMU
        D24(GI)=0.
        D25(GI)=0.
        D26(GI)=0.
            
        D33(GI)=FOURTHIMU
        D34(GI)=0.
        D35(GI)=0.
        D36(GI)=0.
        
        D44(GI)=MUGIXX
        D45(GI)=0.
        D46(GI)=0.
        
        D55(GI)=MUGIXX
        D56(GI)=0.
        
        D66(GI)=MUGIXX
      
        IF(BINGHAM) THEN
          DO ILOC=1,MLOC
            ! indexing               
            IGLP=PNDGLN((ELE-1)*MLOC +ILOC)
            
            PGI(GI) = PGI(GI)+M(ILOC,GI)*P(IGLP)
          
          ENDDO
        ENDIF
  
        IF(BINGHAM) THEN
          ! to work out stress and strain
          DO ILOC=1,NLOC
            ! indexing
            IGLV=VONDGL((ELE-1)*NLOC+ILOC)
            
            ! previous iterate displacement differentiated at gauss points
            NXDGIX=NXDGIX+NX(ILOC,GI)*NXDISP(IGLV)
            NXDGIY=NXDGIY+NY(ILOC,GI)*NXDISP(IGLV)
            NXDGIZ=NXDGIZ+NZ(ILOC,GI)*NXDISP(IGLV)
            
            NYDGIX=NYDGIX+NX(ILOC,GI)*NYDISP(IGLV)
            NYDGIY=NYDGIY+NY(ILOC,GI)*NYDISP(IGLV)
            NYDGIZ=NYDGIZ+NZ(ILOC,GI)*NYDISP(IGLV)
            
            NZDGIX=NZDGIX+NX(ILOC,GI)*NZDISP(IGLV)
            NZDGIY=NZDGIY+NY(ILOC,GI)*NZDISP(IGLV)
            NZDGIZ=NZDGIZ+NZ(ILOC,GI)*NZDISP(IGLV)
  
          ! iloc=1,nloc
          ENDDO
          
          K = 0
          DGAMMA = 0.
  
  4585   CONTINUE
  
          TAXXGI =  D11(GI)*NXDGIX + D12(GI)*NYDGIY + D13(GI)*NZDGIZ &
                  + D14(GI)*(NYDGIZ + NZDGIY) &
                  + D15(GI)*(NXDGIZ + NZDGIX) &
                  + D16(GI)*(NXDGIY + NYDGIX) &
                  - ( D11(GI)*STVPXX(SGI) + D12(GI)*STVPYY(SGI) &
                    + D13(GI)*STVPZZ(SGI) + D14(GI)*STVPYZ(SGI) &
                    + D15(GI)*STVPXZ(SGI) + D16(GI)*STVPXY(SGI) )
                
          TAYYGI =  D12(GI)*NXDGIX + D22(GI)*NYDGIY + D23(GI)*NZDGIZ &
                  + D24(GI)*(NYDGIZ + NZDGIY) &
                  + D25(GI)*(NXDGIZ + NZDGIX) &
                  + D26(GI)*(NXDGIY + NYDGIX) &
                  - ( D12(GI)*STVPXX(SGI) + D22(GI)*STVPYY(SGI) &
                    + D23(GI)*STVPZZ(SGI) + D24(GI)*STVPYZ(SGI) &
                    + D25(GI)*STVPXZ(SGI) + D26(GI)*STVPXY(SGI) )
                  
          TAZZGI =  D13(GI)*NXDGIX + D23(GI)*NYDGIY + D33(GI)*NZDGIZ &
                  + D34(GI)*(NYDGIZ + NZDGIY) &
                  + D35(GI)*(NXDGIZ + NZDGIX) &
                  + D36(GI)*(NXDGIY + NYDGIX) &
                  - ( D13(GI)*STVPXX(SGI) + D23(GI)*STVPYY(SGI) &
                    + D33(GI)*STVPZZ(SGI) + D34(GI)*STVPYZ(SGI) &
                    + D35(GI)*STVPXZ(SGI) + D36(GI)*STVPXY(SGI) )
                  
          TAYZGI =  D14(GI)*NXDGIX + D24(GI)*NYDGIY + D34(GI)*NZDGIZ &
                  + D44(GI)*(NYDGIZ + NZDGIY) &
                  + D45(GI)*(NXDGIZ + NZDGIX) &
                  + D46(GI)*(NXDGIY + NYDGIX) &
                  - ( D14(GI)*STVPXX(SGI) + D24(GI)*STVPYY(SGI) &
                    + D34(GI)*STVPZZ(SGI) + D44(GI)*STVPYZ(SGI) &
                    + D45(GI)*STVPXZ(SGI) + D46(GI)*STVPXY(SGI) )
                  
          TAXZGI =  D15(GI)*NXDGIX + D25(GI)*NYDGIY + D35(GI)*NZDGIZ &
                  + D45(GI)*(NYDGIZ + NZDGIY) &
                  + D55(GI)*(NXDGIZ + NZDGIX) &
                  + D56(GI)*(NXDGIY + NYDGIX) &
                  - ( D15(GI)*STVPXX(SGI) + D25(GI)*STVPYY(SGI) &
                    + D35(GI)*STVPZZ(SGI) + D45(GI)*STVPYZ(SGI) &
                    + D55(GI)*STVPXZ(SGI) + D56(GI)*STVPXY(SGI) )
                  
          TAXYGI = D16(GI)*NXDGIX + D26(GI)*NYDGIY + D36(GI)*NZDGIZ &
                  + D46(GI)*(NYDGIZ + NZDGIY) &
                  + D56(GI)*(NXDGIZ + NZDGIX) &
                  + D66(GI)*(NXDGIY + NYDGIX) &
                  - ( D16(GI)*STVPXX(SGI) + D26(GI)*STVPYY(SGI) &
                    + D36(GI)*STVPZZ(SGI) + D46(GI)*STVPYZ(SGI) &
                    + D56(GI)*STVPXZ(SGI) + D66(GI)*STVPXY(SGI) )
          
          ! ! to calculate the principal deviatoric stresses: tau1, tau2, tau3
          ! ! trace of tau over three
          ! TRTAO3 = ONETHI*(TAUXXGI+TAUYYGI+TAUZZGI)
          ! ! determinant of tau minus i*trtao3 over two
          ! DTMTO2 =  ONEHAL*(  (TAUXXGI - TRTAO3)* &
          !                   ( (TAUYYGI - TRTAO3)*(TAUZZGI - TRTAO3) - TAUYZGI*TAUYZGI ) &
          !                   - TAUXYGI*( TAUXYGI*(TAUZZGI - TRTAO3) - TAUYZGI*TAUXZGI ) &
          !                   + TAUXZGI*( TAUXYGI*TAUYZGI - (TAUYYGI - TRTAO3)*TAUXZGI)  )
          ! ! sum of squares of tau minus i*trtao3 over six
          ! STMTO6 = ONESIX*(  (TAUXXGI - TRTAO3)*(TAUXXGI - TRTAO3) &
          !                 + 2.*TAUXYGI*TAUXYGI + 2.*TAUXZGI*TAUXZGI &
          !                 + (TAUYYGI - TRTAO3)*(TAUYYGI - TRTAO3) &
          !                 + (TAUZZGI - TRTAO3)*(TAUZZGI - TRTAO3) &
          !                 + 2.*TAUYZGI*TAUYZ  )
          !     
          ! PHI = ONETHI*ATAN(  SQRT( STMTO6*STMTO6*STMTO6 - DTMTO2*DTMTO2 )/DTMTO2  )
          ! SINPHI = SIN(PHI)
          ! COSPHI = COS(PHI)
          ! RTSTMT = SQRT(STMTO6)
          ! TAU1 = TRTAO3 + 2.*RTSTMT*COSPHI
          ! TAU2 = TRTAO3 - RTSTMT*(COSPHI + SQRT(3.)*SINPHI)
          ! TAU3 = TRTAO3 - RTSTMT*(COSPHI - SQRT(3.)*SINPHI)
          ! 
          ! TAUMAX = MAX(TAU1, TAU2, TAU3)
          ! TAUMIN = MIN(TAU1, TAU2, TAU3)
          ! TAUTMP = TAU1 + TAU2 + TAU3 - (TAUMAX + TAUMIN)
          ! 
          ! TAU1 = TAUMAX
          ! TAU2 = TAUTMP
          ! TAU3 = TAUMIN
              
          J2 = ONEHAL*( TAXXGI*TAXXGI &
                      + TAYYGI*TAYYGI &
                      + TAZZGI*TAZZGI ) &
              + TAYZGI*TAYZGI &
              + TAXZGI*TAXZGI &
              + TAXYGI*TAXYGI
          
          J3 = ONETHI*( TAXXGI*TAXXGI*TAXXGI &
                      + TAYYGI*TAYYGI*TAYYGI &
                      + TAZZGI*TAZZGI*TAZZGI ) &
              + TAXXGI*( TAXYGI*TAXYGI &
                        + TAXZGI*TAXZGI ) &
              + TAYYGI*( TAXYGI*TAXYGI &
                        + TAYZGI*TAYZGI ) &
              + TAZZGI*( TAXZGI*TAXZGI &
                        + TAYZGI*TAYZGI ) &
              + 2.*TAYZGI*TAXZGI*TAXYGI
                    
          SQRTJ2 = SQRT(J2)
                      
          LODETH = ONETHI*ASIN((-1.5)*SQRTHR*J3/(J2*SQRTJ2))
                      
          ! COHESI=9.8*0.05*DENGI(GI)
  
  !         should this have a PGI in????????? ----->          
  !         F(SGI) = SQRTHR*SQRTJ2 - 6.*SINFAN*NDENOM - 6.*COHESI*COSFAN*NDENOM
          F(SGI) = SQRTHR*SQRTJ2 - 6.*SINFAN*NDENOM*PGI(GI) - 6.*COHESI*COSFAN*NDENOM
              
          IF(K.EQ.0) FOLD(SGI) = F(SGI)
              
          IF(F(SGI).GE.TOL) THEN
            OFFYIELDENV = .TRUE.
                          
            INVPGI = -1./(9.*PGI(GI))
                        
            ! M^{0} = matrix of coefficients for derivative of p wrt stress
            DPDS11 = INVPGI
            DPDS12 = INVPGI
            DPDS13 = INVPGI
            DPDS14 = 0.
            DPDS15 = 0.
            DPDS16 = 0.
            
            DPDS22 = INVPGI
            DPDS23 = INVPGI
            DPDS24 = 0.
            DPDS25 = 0.
            DPDS26 = 0.
            
            DPDS33 = INVPGI
            DPDS34 = 0.
            DPDS35 = 0.
            DPDS36 = 0.
            
            DPDS44 = 0.
            DPDS45 = 0.
            DPDS46 = 0.
            
            DPDS55 = 0.
            DPDS56 = 0.
            
            DPDS66 = 0.
            
            ! dM^{0}/ds = matrix of second derivatives of p wrt stress
            D2PDS211 = 0.
            D2PDS212 = 0.
            D2PDS213 = 0.
            D2PDS214 = 0.
            D2PDS215 = 0.
            D2PDS216 = 0.
            
            D2PDS222 = 0.
            D2PDS223 = 0.
            D2PDS224 = 0.
            D2PDS225 = 0.
            D2PDS226 = 0.
            
            D2PDS233 = 0.
            D2PDS234 = 0.
            D2PDS235 = 0.
            D2PDS236 = 0.
            
            D2PDS244 = 0.
            D2PDS245 = 0.
            D2PDS246 = 0.
            
            D2PDS255 = 0.
            D2PDS256 = 0.
            
            D2PDS266 = 0.
            
            ! M^{I} = matrix of coefficients for derivative of j2 wrt stress
            DJ2D11 = TWOTHI
            DJ2D12 = -ONETHI
            DJ2D13 = -ONETHI
            DJ2D14 = 0.
            DJ2D15 = 0.
            DJ2D16 = 0.
            
            DJ2D22 = TWOTHI
            DJ2D23 = -ONETHI
            DJ2D24 = 0.
            DJ2D25 = 0.
            DJ2D26 = 0.
            
            DJ2D33 = TWOTHI
            DJ2D34 = 0.
            DJ2D35 = 0.
            DJ2D36 = 0.
            
            DJ2D44 = 2.
            DJ2D45 = 0.
            DJ2D46 = 0.
            
            DJ2D55 = 2.
            DJ2D56 = 0.
            
            DJ2D66 = 2.
            
            ! dM^{I}/ds = matrix of second derivatives of j2 wrt stress
            D2J2D211 = TWOTHI
            D2J2D212 = -ONETHI
            D2J2D213 = -ONETHI
            D2J2D214 = 0.
            D2J2D215 = 0.
            D2J2D216 = 0.
            
            D2J2D222 = TWOTHI
            D2J2D223 = -ONETHI
            D2J2D224 = 0.
            D2J2D225 = 0.
            D2J2D226 = 0.
            
            D2J2D233 = TWOTHI
            D2J2D234 = 0.
            D2J2D235 = 0.
            D2J2D236 = 0.
            
            D2J2D244 = 2.
            D2J2D245 = 0.
            D2J2D246 = 0.
            
            D2J2D255 = 2.
            D2J2D256 = 0.
            
            D2J2D266 = 2.
            
            ! M^{II} = matrix of coefficients for derivative of j3 wrt stress
            DJ3D11 = ONETHI*TAXXGI
            DJ3D12 = ONETHI*TAZZGI
            DJ3D13 = ONETHI*TAYYGI
            DJ3D14 = -TWOTHI*TAYZGI
            DJ3D15 = ONETHI*TAXZGI
            DJ3D16 = ONETHI*TAXYGI
            
            DJ3D22 = ONETHI*TAYYGI
            DJ3D23 = ONETHI*TAXXGI
            DJ3D24 = ONETHI*TAYZGI
            DJ3D25 = -TWOTHI*TAXZGI
            DJ3D26 = ONETHI*TAXYGI
            
            DJ3D33 = ONETHI*TAZZGI
            DJ3D34 = ONETHI*TAYZGI
            DJ3D35 = ONETHI*TAXZGI
            DJ3D36 = -TWOTHI*TAXYGI
            
            DJ3D44 = -TAXXGI
            DJ3D45 = TAXYGI
            DJ3D46 = TAXZGI
            
            DJ3D55 = -TAYYGI
            DJ3D56 = TAYZGI
            
            DJ3D66 = -TAZZGI
            
            ! dM^{II}/ds = matrix of second derivatives of j3 wrt stress
            D2J3D211 = TWOTHI*TAXXGI
            D2J3D212 = TWOTHI*TAZZGI
            D2J3D213 = TWOTHI*TAYYGI
            D2J3D214 = -FOURTHI*TAYZGI
            D2J3D215 = TWOTHI*TAXZGI
            D2J3D216 = TWOTHI*TAXYGI
            
            D2J3D222 = TWOTHI*TAYYGI
            D2J3D223 = TWOTHI*TAXXGI
            D2J3D224= TWOTHI*TAYZGI
            D2J3D225 = -FOURTHI*TAXZGI
            D2J3D226 = TWOTHI*TAXYGI
            
            D2J3D233 = TWOTHI*TAZZGI
            D2J3D234 = TWOTHI*TAYZGI
            D2J3D235 = TWOTHI*TAXZGI
            D2J3D236 = -FOURTHI*TAXYGI
            
            D2J3D244 = -2.*TAXXGI
            D2J3D245 = 2.*TAXYGI
            D2J3D246 = 2.*TAXZGI
            
            D2J3D255 = -2.*TAYYGI
            D2J3D256 = 2.*TAYZGI
            
            D2J3D266 = -2.*TAZZGI
            
            ! total stresses (11,22,33)
            SRXXGI = TAXXGI - PGI(GI)
            SRYYGI = TAYYGI - PGI(GI)
            SRZZGI = TAZZGI - PGI(GI)
            
            ! dqds__ = vector components of flow surface derivative wrt stress
            DQDSXX =  (DQDP*DPDS11 + DQDJ2*DJ2D11 + DQDJ3*DJ3D11)*SRXXGI &
                    + (DQDP*DPDS12 + DQDJ2*DJ2D12 + DQDJ3*DJ3D12)*SRYYGI &
                    + (DQDP*DPDS13 + DQDJ2*DJ2D13 + DQDJ3*DJ3D13)*SRZZGI &
                    + (DQDP*DPDS14 + DQDJ2*DJ2D14 + DQDJ3*DJ3D14)*TAYZGI &
                    + (DQDP*DPDS15 + DQDJ2*DJ2D15 + DQDJ3*DJ3D15)*TAXZGI &
                    + (DQDP*DPDS16 + DQDJ2*DJ2D16 + DQDJ3*DJ3D16)*TAXYGI
                    
            DQDSYY =  (DQDP*DPDS12 + DQDJ2*DJ2D12 + DQDJ3*DJ3D12)*SRXXGI &
                    + (DQDP*DPDS22 + DQDJ2*DJ2D22 + DQDJ3*DJ3D22)*SRYYGI &
                    + (DQDP*DPDS23 + DQDJ2*DJ2D23 + DQDJ3*DJ3D23)*SRZZGI &
                    + (DQDP*DPDS24 + DQDJ2*DJ2D24 + DQDJ3*DJ3D24)*TAYZGI &
                    + (DQDP*DPDS25 + DQDJ2*DJ2D25 + DQDJ3*DJ3D25)*TAXZGI &
                    + (DQDP*DPDS26 + DQDJ2*DJ2D26 + DQDJ3*DJ3D26)*TAXYGI
                    
            DQDSZZ =  (DQDP*DPDS13 + DQDJ2*DJ2D13 + DQDJ3*DJ3D13)*SRXXGI &
                    + (DQDP*DPDS23 + DQDJ2*DJ2D23 + DQDJ3*DJ3D23)*SRYYGI &
                    + (DQDP*DPDS33 + DQDJ2*DJ2D33 + DQDJ3*DJ3D33)*SRZZGI &
                    + (DQDP*DPDS34 + DQDJ2*DJ2D34 + DQDJ3*DJ3D34)*TAYZGI &
                    + (DQDP*DPDS35 + DQDJ2*DJ2D35 + DQDJ3*DJ3D35)*TAXZGI &
                    + (DQDP*DPDS36 + DQDJ2*DJ2D36 + DQDJ3*DJ3D36)*TAXYGI
                    
            DQDSYZ =  (DQDP*DPDS14 + DQDJ2*DJ2D14 + DQDJ3*DJ3D14)*SRXXGI &
                    + (DQDP*DPDS24 + DQDJ2*DJ2D24 + DQDJ3*DJ3D24)*SRYYGI &
                    + (DQDP*DPDS34 + DQDJ2*DJ2D34 + DQDJ3*DJ3D34)*SRZZGI &
                    + (DQDP*DPDS44 + DQDJ2*DJ2D44 + DQDJ3*DJ3D44)*TAYZGI &
                    + (DQDP*DPDS45 + DQDJ2*DJ2D45 + DQDJ3*DJ3D45)*TAXZGI &
                    + (DQDP*DPDS46 + DQDJ2*DJ2D46 + DQDJ3*DJ3D46)*TAXYGI
                    
            DQDSXZ =  (DQDP*DPDS15 + DQDJ2*DJ2D15 + DQDJ3*DJ3D15)*SRXXGI &
                    + (DQDP*DPDS25 + DQDJ2*DJ2D25 + DQDJ3*DJ3D25)*SRYYGI &
                    + (DQDP*DPDS35 + DQDJ2*DJ2D35 + DQDJ3*DJ3D35)*SRZZGI &
                    + (DQDP*DPDS45 + DQDJ2*DJ2D45 + DQDJ3*DJ3D45)*TAYZGI &
                    + (DQDP*DPDS55 + DQDJ2*DJ2D55 + DQDJ3*DJ3D55)*TAXZGI &
                    + (DQDP*DPDS56 + DQDJ2*DJ2D56 + DQDJ3*DJ3D56)*TAXYGI
                    
            DQDSXY =  (DQDP*DPDS16 + DQDJ2*DJ2D16 + DQDJ3*DJ3D16)*SRXXGI &
                    + (DQDP*DPDS26 + DQDJ2*DJ2D26 + DQDJ3*DJ3D26)*SRYYGI &
                    + (DQDP*DPDS36 + DQDJ2*DJ2D36 + DQDJ3*DJ3D36)*SRZZGI &
                    + (DQDP*DPDS46 + DQDJ2*DJ2D46 + DQDJ3*DJ3D46)*TAYZGI &
                    + (DQDP*DPDS56 + DQDJ2*DJ2D56 + DQDJ3*DJ3D56)*TAXZGI &
                    + (DQDP*DPDS66 + DQDJ2*DJ2D66 + DQDJ3*DJ3D66)*TAXYGI
            
            ! dfds__ = vector components of yield surface derivative wrt stress
            DFDSXX =  (DFDP*DPDS11 + DFDJ2*DJ2D11 + DFDJ3*DJ3D11)*SRXXGI &
                    + (DFDP*DPDS12 + DFDJ2*DJ2D12 + DFDJ3*DJ3D12)*SRYYGI &
                    + (DFDP*DPDS13 + DFDJ2*DJ2D13 + DFDJ3*DJ3D13)*SRZZGI &
                    + (DFDP*DPDS14 + DFDJ2*DJ2D14 + DFDJ3*DJ3D14)*TAYZGI &
                    + (DFDP*DPDS15 + DFDJ2*DJ2D15 + DFDJ3*DJ3D15)*TAXZGI &
                    + (DFDP*DPDS16 + DFDJ2*DJ2D16 + DFDJ3*DJ3D16)*TAXYGI
                      
            DFDSYY =  (DFDP*DPDS12 + DFDJ2*DJ2D12 + DFDJ3*DJ3D12)*SRXXGI &
                    + (DFDP*DPDS22 + DFDJ2*DJ2D22 + DFDJ3*DJ3D22)*SRYYGI &
                    + (DFDP*DPDS23 + DFDJ2*DJ2D23 + DFDJ3*DJ3D23)*SRZZGI &
                    + (DFDP*DPDS24 + DFDJ2*DJ2D24 + DFDJ3*DJ3D24)*TAYZGI &
                    + (DFDP*DPDS25 + DFDJ2*DJ2D25 + DFDJ3*DJ3D25)*TAXZGI &
                    + (DFDP*DPDS26 + DFDJ2*DJ2D26 + DFDJ3*DJ3D26)*TAXYGI
                    
            DFDSZZ =  (DFDP*DPDS13 + DFDJ2*DJ2D13 + DFDJ3*DJ3D13)*SRXXGI &
                    + (DFDP*DPDS23 + DFDJ2*DJ2D23 + DFDJ3*DJ3D23)*SRYYGI &
                    + (DFDP*DPDS33 + DFDJ2*DJ2D33 + DFDJ3*DJ3D33)*SRZZGI &
                    + (DFDP*DPDS34 + DFDJ2*DJ2D34 + DFDJ3*DJ3D34)*TAYZGI &
                    + (DFDP*DPDS35 + DFDJ2*DJ2D35 + DFDJ3*DJ3D35)*TAXZGI &
                    + (DFDP*DPDS36 + DFDJ2*DJ2D36 + DFDJ3*DJ3D36)*TAXYGI
                    
            DFDSYZ =  (DFDP*DPDS14 + DFDJ2*DJ2D14 + DFDJ3*DJ3D14)*SRXXGI &
                    + (DFDP*DPDS24 + DFDJ2*DJ2D24 + DFDJ3*DJ3D24)*SRYYGI &
                    + (DFDP*DPDS34 + DFDJ2*DJ2D34 + DFDJ3*DJ3D34)*SRZZGI &
                    + (DFDP*DPDS44 + DFDJ2*DJ2D44 + DFDJ3*DJ3D44)*TAYZGI &
                    + (DFDP*DPDS45 + DFDJ2*DJ2D45 + DFDJ3*DJ3D45)*TAXZGI &
                    + (DFDP*DPDS46 + DFDJ2*DJ2D46 + DFDJ3*DJ3D46)*TAXYGI
                        
            DFDSXZ =  (DFDP*DPDS15 + DFDJ2*DJ2D15 + DFDJ3*DJ3D15)*SRXXGI &
                    + (DFDP*DPDS25 + DFDJ2*DJ2D25 + DFDJ3*DJ3D25)*SRYYGI &
                    + (DFDP*DPDS35 + DFDJ2*DJ2D35 + DFDJ3*DJ3D35)*SRZZGI &
                    + (DFDP*DPDS45 + DFDJ2*DJ2D45 + DFDJ3*DJ3D45)*TAYZGI &
                    + (DFDP*DPDS55 + DFDJ2*DJ2D55 + DFDJ3*DJ3D55)*TAXZGI &
                    + (DFDP*DPDS56 + DFDJ2*DJ2D56 + DFDJ3*DJ3D56)*TAXYGI
            
            DFDSXY =  (DFDP*DPDS16 + DFDJ2*DJ2D16 + DFDJ3*DJ3D16)*SRXXGI &
                    + (DFDP*DPDS26 + DFDJ2*DJ2D26 + DFDJ3*DJ3D26)*SRYYGI &
                    + (DFDP*DPDS36 + DFDJ2*DJ2D36 + DFDJ3*DJ3D36)*SRZZGI &
                    + (DFDP*DPDS46 + DFDJ2*DJ2D46 + DFDJ3*DJ3D46)*TAYZGI &
                    + (DFDP*DPDS56 + DFDJ2*DJ2D56 + DFDJ3*DJ3D56)*TAXZGI &
                    + (DFDP*DPDS66 + DFDJ2*DJ2D66 + DFDJ3*DJ3D66)*TAXYGI
                
            ! ! d2qds2__ = matrix components of flow surface second derivative wrt stress
            ! !            note that this is a simplified version only suitable for
            ! !            von mises or drucker prager surfaces
            ! D2QDS211 = DQDP*D2PDS211+DQDJ2*D2J2D211+DQDJ3*D2J3D211 &
            !            D2QDP2*... ! this bit involves matrix multiplication...
            ! D2QDS212 = DQDP*D2PDS212+DQDJ2*D2J2D212+DQDJ3*D2J3D212
            ! D2QDS213 = DQDP*D2PDS213+DQDJ2*D2J2D213+DQDJ3*D2J3D213
            ! D2QDS214 = DQDP*D2PDS214+DQDJ2*D2J2D214+DQDJ3*D2J3D214
            ! D2QDS215 = DQDP*D2PDS215+DQDJ2*D2J2D215+DQDJ3*D2J3D215
            ! D2QDS216 = DQDP*D2PDS216+DQDJ2*D2J2D216+DQDJ3*D2J3D216
            ! 
            ! D2QDS222 = DQDP*D2PDS222+DQDJ2*D2J2D222+DQDJ3*D2J3D222
            ! D2QDS223 = DQDP*D2PDS223+DQDJ2*D2J2D223+DQDJ3*D2J3D223
            ! D2QDS224 = DQDP*D2PDS224+DQDJ2*D2J2D224+DQDJ3*D2J3D224
            ! D2QDS225 = DQDP*D2PDS225+DQDJ2*D2J2D225+DQDJ3*D2J3D225
            ! D2QDS226 = DQDP*D2PDS226+DQDJ2*D2J2D226+DQDJ3*D2J3D226
            ! 
            ! D2QDS233 = DQDP*D2PDS233+DQDJ2*D2J2D233+DQDJ3*D2J3D233
            ! D2QDS234 = DQDP*D2PDS234+DQDJ2*D2J2D234+DQDJ3*D2J3D234
            ! D2QDS235 = DQDP*D2PDS235+DQDJ2*D2J2D235+DQDJ3*D2J3D235
            ! D2QDS236 = DQDP*D2PDS236+DQDJ2*D2J2D236+DQDJ3*D2J3D236
            ! 
            ! D2QDS244 = DQDP*D2PDS244+DQDJ2*D2J2D244+DQDJ3*D2J3D244
            ! D2QDS245 = DQDP*D2PDS245+DQDJ2*D2J2D245+DQDJ3*D2J3D245
            ! D2QDS246 = DQDP*D2PDS246+DQDJ2*D2J2D246+DQDJ3*D2J3D246
            ! 
            ! D2QDS255 = DQDP*D2PDS255+DQDJ2*D2J2D255+DQDJ3*D2J3D255
            ! D2QDS256 = DQDP*D2PDS256+DQDJ2*D2J2D256+DQDJ3*D2J3D256
            ! 
            ! D2QDS266 = DQDP*D2PDS266+DQDJ2*D2J2D266+DQDJ3*D2J3D266
            !               
            ! ! xi^{-1} matrix = hessian?
            ! XIINV11 = D11 + DGAMMA*D2QDS211
            ! XIINV12 = D12 + DGAMMA*D2QDS212
            ! XIINV13 = D13 + DGAMMA*D2QDS213
            ! XIINV14 = D14 + DGAMMA*D2QDS214
            ! XIINV15 = D15 + DGAMMA*D2QDS215
            ! XIINV16 = D16 + DGAMMA*D2QDS216
            ! 
            ! XIINV22 = D22 + DGAMMA*D2QDS222
            ! XIINV23 = D23 + DGAMMA*D2QDS223
            ! XIINV24 = D24 + DGAMMA*D2QDS224
            ! XIINV25 = D25 + DGAMMA*D2QDS225
            ! XIINV26 = D26 + DGAMMA*D2QDS226
            ! 
            ! XIINV33 = D33 + DGAMMA*D2QDS233
            ! XIINV34 = D34 + DGAMMA*D2QDS234
            ! XIINV35 = D35 + DGAMMA*D2QDS235
            ! XIINV36 = D36 + DGAMMA*D2QDS236
            ! 
            ! XIINV44 = D44 + DGAMMA*D2QDS244
            ! XIINV45 = D45 + DGAMMA*D2QDS245
            ! XIINV46 = D46 + DGAMMA*D2QDS246
            ! 
            ! XIINV55 = D55 + DGAMMA*D2QDS255
            ! XIINV56 = D56 + DGAMMA*D2QDS256
            ! 
            ! XIINV66 = D66 + DGAMMA*D2QDS266
            
                
            DNORM = DFDSXX*( D11(GI)*DQDSXX + D12(GI)*DQDSYY + D13(GI)*DQDSZZ &
                            + D14(GI)*DQDSYZ + D15(GI)*DQDSXZ + D16(GI)*DQDSXY ) &
                  + DFDSYY*( D12(GI)*DQDSXX + D22(GI)*DQDSYY + D23(GI)*DQDSZZ &
                            + D24(GI)*DQDSYZ + D25(GI)*DQDSXZ + D26(GI)*DQDSXY ) &
                  + DFDSZZ*( D13(GI)*DQDSXX + D23(GI)*DQDSYY + D33(GI)*DQDSZZ &
                            + D34(GI)*DQDSYZ + D35(GI)*DQDSXZ + D36(GI)*DQDSXY ) &
                  + DFDSYZ*( D14(GI)*DQDSXX + D24(GI)*DQDSYY + D34(GI)*DQDSZZ &
                            + D44(GI)*DQDSYZ + D45(GI)*DQDSXZ + D46(GI)*DQDSXY ) &
                  + DFDSXZ*( D15(GI)*DQDSXX + D25(GI)*DQDSYY + D35(GI)*DQDSZZ &
                            + D45(GI)*DQDSYZ + D55(GI)*DQDSXZ + D56(GI)*DQDSXY ) &
                  + DFDSXY*( D16(GI)*DQDSXX + D26(GI)*DQDSYY + D36(GI)*DQDSZZ &
                            + D46(GI)*DQDSYZ + D56(GI)*DQDSXZ + D66(GI)*DQDSXY )
                
            D2GAMMA = F(SGI)/DNORM
                
            STVPXX(SGI) = D2GAMMA*DQDSXX + STVPXX(SGI)
            STVPYY(SGI) = D2GAMMA*DQDSYY + STVPYY(SGI)
            STVPZZ(SGI) = D2GAMMA*DQDSZZ + STVPZZ(SGI)
            STVPYZ(SGI) = D2GAMMA*DQDSYZ + STVPYZ(SGI)
            STVPXZ(SGI) = D2GAMMA*DQDSXZ + STVPXZ(SGI)
            STVPXY(SGI) = D2GAMMA*DQDSXY + STVPXY(SGI)
                        
            DGAMMA = DGAMMA + D2GAMMA
            K = K + 1
            
            GOTO 4585
  
          ! if(f(sgi).ge.0.)  
          ENDIF
              
        ! if(bingham)
        ENDIF
            
      ! gi=1,ngi
      ENDDO
              
    ! ele=1,totele
    ENDDO
          
    ewrite_minmax(stvpxx)
    ewrite_minmax(stvpyy)
    ewrite_minmax(stvpzz)
    ewrite_minmax(stvpyz)
    ewrite_minmax(stvpxz)
    ewrite_minmax(stvpxy)
    
    ewrite_minmax(F)
    ewrite_minmax(FOLD)
    
    ewrite(3,*)  'OFFYIELDENV = ', OFFYIELDENV
    
    ewrite(3,*)  'EXITING STEVAL()'
    RETURN
  END SUBROUTINE STEVAL
  !************************************************************************
  
end module Solid_assembly

