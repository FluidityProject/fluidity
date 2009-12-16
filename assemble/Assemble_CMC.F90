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

module assemble_CMC

  use dgtools
  use fldebug
  use state_module
  use sparse_tools
  use sparse_tools_petsc
  use spud
  use fields
  use sparse_matrices_fields
  use multimaterial_module, only: extract_prognostic_pressure
  use global_parameters, only: OPTION_PATH_LEN
  use allsorts
  use boundary_conditions
  use hart3d_allsorts
  use position_in_matrix
  use shape_transformations
  use sml
  use tr2d_module
  use linked_lists
  
  implicit none 

  ! Buffer for output messages.
  character(len=255), private :: message

  private
  
  public :: assemble_cmc_dg, assemble_masslumped_cmc, cmc_wrapper, &
            assemble_masslumped_ctm, repair_stiff_nodes, zero_stiff_nodes, &
            assemble_diagonal_schur

contains


    !************************************************************************
    SUBROUTINE CMC_WRAPPER( &
                        NONODS, XNONOD, TOTELE &
                      , GETC12, NCT, FINDCT, COLCT &
                      , CMCGET, FREDOP, FINCMC, COLCMC, MIDCMC &
                      , NDGLNO, PNDGLN, XNDGLN, SONDGLN &
                      , NCMC, C1T, C2T, C3T, C1TP, C2TP, C3TP &
                      , PREOPT, X, Y, Z & 
                      , D3, MKCOMP &
                      , NBUOY,BOY_ML,NDWISP &
                      , NLOC, NGI, N, NLX, NLY, NLZ &
                      , MLOC, M, MLX, MLY, MLZ, WEIGHT &
                      , DM1, DM2, DM3 &
                      , NDPSET, ML &
                      , CMC, KCMC &
                      , ROTAT, NNODRO, NODROT &
                      , NORMX, NORMY, NORMZ, T1X, T1Y, T1Z, T2X, T2Y, T2Z &
                      , ISPHERE &
                      , NOFILT &
                      , state &
                      )
      ! =============================================================
      ! Subroutine to construct the matrices CMC, C1/2/3T and 
      ! (if MKCOMP = 3) C1/2/3PT.
      ! =============================================================
    
      ! =============================================================
      ! inputs 
      ! =============================================================

      ! bucket full of fields
      type(state_type), dimension(:), intent(inout) :: state

      INTEGER, INTENT(IN) :: NONODS, XNONOD, TOTELE
      ! nonods = total number of nodes
      ! xnonod = length of co-ordinate vectors
      ! totele = total number of elements
    
      LOGICAL, INTENT(IN) :: GETC12, CMCGET
      INTEGER, INTENT(IN) :: NBUOY
      REAL, INTENT(IN) :: BOY_ML(NBUOY*9*NONODS)
      ! getc12 = if true assemble the c1/2/3t(nct) "matrices"
      ! cmcget = if true assemble the cmc matrix
      
      INTEGER, INTENT(IN) :: NCT, FREDOP, NCMC
      ! nct    = length of c1/2/3t vectors
      ! fredop = number of pressure degrees of freedom
      ! ncmc   = length of colcmc
      
      INTEGER, INTENT(IN) :: MKCOMP
      ! mkcomp       = determines if material is compressible or not
      !                mkcomp.ge.1 -> compressible
      
      INTEGER, INTENT(IN) :: NLOC, MLOC, NGI
      ! nloc  = number of displacement nodes per element
      ! mloc  = number of pressure nodes per element
      ! ngi   = number of gauss points (for current element)
      
      LOGICAL, INTENT(IN) :: D3
      ! d3   = flag for three dimensional simulation
      
      INTEGER, INTENT(IN) :: FINDCT(FREDOP+1), COLCT(NCT)
      ! findct(fredop+1) = integer array indexing the first entry of each
      !                    matrix "row" in c1/2/3t(nct) corresponding to a particular node
      ! colct(nct)       = integer array indexing equivalent columns corresponding
      !                    to entries of c1/2/3t(nct)
      
      INTEGER, INTENT(IN) :: FINCMC(FREDOP+1), COLCMC(NCMC), MIDCMC(FREDOP)
      ! fincmc(fredop+1) = array indexing first entry of each row of cmc
      !                    used here to ensure mcompl has same sparsity structure as cmc
      ! colcmc(ncmc)     = array indexing equivalent columns of cmc matrix
      ! midcmc(fredop)   = array indexing diagonal entries of cmc matrix - not used here
      
      INTEGER, INTENT(IN) :: NDGLNO(TOTELE*NLOC), PNDGLN(TOTELE*MLOC)
      ! ndglno(totele*nloc) = element pointers for unknowns
      ! pndgln(totele*mloc) = element pointer for pressure terms
      
      INTEGER, INTENT(IN) :: XNDGLN(TOTELE*NLOC), SONDGLN(TOTELE*NLOC)
      ! xndgln(totele*nloc) = element pointer for co-ordinates
      ! sondgln(totele*nloc) = element pointers for materials and sources
      
      REAL, INTENT(IN) :: N(NLOC, NGI), M(MLOC, NGI)
      ! n(nloc,ngi) = displacement shape function 
      ! m(mloc,ngi) = pressure shape function
      !               (evaluated at a gauss point around a node)
      
      REAL, INTENT(IN) :: NLX(NLOC, NGI), NLY(NLOC, NGI), NLZ(NLOC, NGI)
      ! nlx(nloc,ngi) = derivative of displacement shape function in local x direction
      !                 (evaluated at a gauss point around a node)
      ! nly(nloc,ngi) = derivative of displacement shape function in local y direction
      !                 (evaluated at a gauss point around a node)
      ! nlz(nloc,ngi) = derivative of displacement shape function in local z direction
      !                 (evaluated at a gauss point around a node)
      
      REAL, INTENT(IN) :: MLX(NLOC, NGI), MLY(NLOC, NGI), MLZ(NLOC, NGI)
      ! mlx(mloc,ngi) = derivative of pressure shape function in local x direction
      !                 (evaluated at a gauss point around a node)
      ! mly(mloc,ngi) = derivative of pressure shape function in local y direction
      !                 (evaluated at a gauss point around a node)
      ! mlz(mloc,ngi) = derivative of pressure shape function in local z direction
      !                 (evaluated at a gauss point around a node)
      
      REAL, INTENT(IN) :: X(XNONOD), Y(XNONOD), Z(XNONOD)
      ! x(nonod) = x co-ordinates of the nodes
      ! y(nonod) = y co-ordinates of the nodes
      ! z(nonod) = z co-ordinates of the nodes
      
      INTEGER, INTENT(IN) :: PREOPT
      ! preopt = indicates whether or not to integrate the continuity eqn using green's
      !          preopt = 1 -> green's
      !          preopt = 0 -> no green's
      
      REAL, INTENT(IN) :: WEIGHT(NGI)
      ! weight(ngi) = weights for gaussian quadrature at each of the gauss points
      
      REAL, INTENT(IN) :: ML(NONODS)
      ! ml(nonods)    = lumped mass matrix (including density term)
      
      REAL, INTENT(IN) :: DM1(NONODS), DM2(NONODS), DM3(NONODS) 
      ! dm1(nonods) = dirichlet boundary conditions on the x-component of velocity
      ! dm2(nonods) = dirichlet boundary conditions on the y-component of velocity
      ! dm3(nonods) = dirichlet boundary conditions on the z-component of velocity
      
      LOGICAL, INTENT(IN) :: ROTAT, NDWISP
      ! rotat  = flag indicating whether to rotat the co-ordinate system or not
      !          rotat = .true.  -> rotate the co-ordinates
      !          rotat = .false. -> do not rotate co-ordinates
      ! ndwisp = flag indicating whether the velocity and pressure nodes coincide
      !          ndwisp = .true.  -> nodes coincide
      !          ndwisp = .false. -> nodes different
      
      INTEGER, INTENT(IN) :: NDPSET, NNODRO
      ! ndpset = number of node where pressure is set to zero, 
      !          if all velocity boundaries are dirichlet
      ! nnodro = lenght of rotation vectors (number of nodes to be rotated?)
      
      INTEGER, INTENT(IN) :: NODROT(NNODRO)
      ! nodrot = vector indexing nodes to be rotated?
      
      REAL, INTENT(IN) :: NORMX(NNODRO), NORMY(NNODRO), NORMZ(NNODRO)
      REAL, INTENT(IN) :: T1X(NNODRO), T1Y(NNODRO), T1Z(NNODRO)
      REAL, INTENT(IN) :: T2X(NNODRO), T2Y(NNODRO), T2Z(NNODRO)
      ! [t1x   t1y   t1z;
      !  t2x   t2y   t2z;
      !  nx    ny    nz ] = components of rotation matrix (in 3D)
      ! [t1x   t1y;
      !  nx    ny ] = components of rotation matrix (in 2D)
      
      INTEGER :: ISPHERE, NOFILT
    
      ! =============================================================
      ! outputs
      ! =============================================================
      REAL, INTENT(INOUT) :: C1T(NCT), C2T(NCT), C3T(NCT)
      ! c1t(nct) = discretisation of x component of divergence operator
      ! c2t(nct) = discretisation of y component of divergence operator
      ! c3t(nct) = discretisation of z component of divergence operator
      
      REAL, INTENT(INOUT) :: C1TP(NCT), C2TP(NCT), C3TP(NCT)
      ! c1tp(nct) = discretisation of x component of density modified 
      !             divergence operator
      ! c2tp(nct) = discretisation of y component of density modified 
      !             divergence operator
      ! c3tp(nct) = discretisation of z component of density modified 
      !             divergence operator
    
      REAL, INTENT(INOUT) :: CMC(NCMC), KCMC(NCMC)
      ! cmc(ncmc)  = cmc matrix constructed in getcmc
      ! kcmc(ncmc) = 4th order pressure stabilisation constructed in cmcfilt
    

      ! =============================================================
      ! local variables
      ! =============================================================
      REAL, ALLOCATABLE :: ROTATVEC(:)
      ! rotatspace(rotatlen) = vector space used by getcmc if rotat = .true.
      
      INTEGER :: ROTATLEN
      ! rotatlen = length of rotatspace
      
      LOGICAL :: DUMMYROTAT
      ! dummyrotat = dummy flag to be passed to getcmc
      PARAMETER(DUMMYROTAT = .FALSE.)
      
      integer :: stat
      character(len=OPTION_PATH_LEN) :: pressure_option_path

      type(scalar_field), pointer :: pressure

      ! Cause the pressure warnings to only happen once.
      logical, save :: pressure_warned=.false.

      ! =============================================================
      ! 
      ! If CMCGET = .TRUE. then this subroutine calls GETCMC to assemble
      ! the CMC or GMC (using C*TP) matrix.
      ! If we are using coninuous_galerkin tets then CMCFILT is called
      ! to create the stabilisation (fourth order pressure term) matrix KCMC
      ! and add it to CMC (so long as NOFILT is false).
      !
      ! Called from NAVSTO (just before SOLNAV)
      ! 
      ! Description                                   Programmer      Date
      ! ==================================================================
      ! Original version  .............................. CRGW       03/09/06
      ! =============================================================
    
      EWRITE(2,*) 'ENTERING CMC_WRAPPER()'

      pressure_option_path=""
      pressure=>extract_prognostic_pressure(state, stat=stat)
      if(stat==0) then
        pressure_option_path=trim(pressure%option_path)
      else if((stat==1).and.(.not.pressure_warned)) then
        ewrite(0,*) "Warning: No prognostic pressure found in state."
        ewrite(0,*) "Strange that you've got here without a prognostic press&
             &ure."
        pressure_warned=.true.
      else if((stat==2).and.(.not.pressure_warned)) then
        ewrite(0,*) "Warning: Multiple prognostic pressures found."
        ewrite(0,*) "Not sure if new options are compatible with this yet."
        pressure_warned=.true.
      end if

      ! if we need to make cmc then call getcmc
      IF(CMCGET) THEN 
          ROTATLEN = 1
          ALLOCATE( ROTATVEC(ROTATLEN) )
        
          CALL GETCMC( NDPSET, D3, CMC, C1T, C2T, C3T, &
                      ML, ML, ML, DM1, DM2, DM3, &
                      NBUOY,BOY_ML, &
                      NONODS, FREDOP, NCT, NCMC, &
                      FINCMC, COLCMC, MIDCMC, &
                      FINDCT, COLCT, &
                      ! rotations:
                      ROTATVEC, ROTATLEN, &
                      DUMMYROTAT, NNODRO, NODROT, &
                      NORMX, NORMY, NORMZ, T1X, T1Y, T1Z, T2X, T2Y, T2Z &
                      ! compressible:
                      ,MKCOMP, C1TP, C2TP, C3TP &
                      )
    
           DEALLOCATE( ROTATVEC )
          
          ! if our pressure and velocity nodes coincide
          IF(NDWISP) THEN
             
            IF((have_option(trim(pressure_option_path)//&
                '/prognostic/spatial_discretisation/control_volumes')).OR.&
               ( NOFILT.EQ.1 )) THEN
              KCMC = 0.
            ELSE

              CALL CMCFILT( X, Y, Z, &
                            KCMC, CMC, NDPSET, D3, &
                            FINCMC, COLCMC, MIDCMC, NCMC, FREDOP, TOTELE, &
                            N, NLX, NLY, NLZ, &
                            M, MLX, MLY, MLZ, WEIGHT, NGI, NLOC, MLOC, &
                            PNDGLN, XNDGLN, XNONOD, &
                            ISPHERE )
             ENDIF
    
          ! if(ndwisp)
          ENDIF
      
      ! if(cmcget)
      ENDIF

      ewrite_minmax(CMC)
      EWRITE(2,*) 'EXITING CMC_WRAPPER()'

      RETURN
    END SUBROUTINE CMC_WRAPPER  
    !************************************************************************

    subroutine assemble_cmc_dg(CMC, CTP, CT, inverse_mass)
      !!< Assemble the pressure matrix C^T M^{-1} C for a DG mesh.
      !!< This currently does not support rotations.
      type(csr_matrix), intent(inout) :: CMC
      type(block_csr_matrix), intent(in) :: CTP, CT
      type(block_csr_matrix), intent(in):: inverse_mass

      type(csr_matrix) :: CM ! Temporary half way matrix
      type(csr_matrix) :: CMC_tmp ! Temporary accumulator.
      type(csr_matrix) :: CT_block ! The current CT dimension
      type(csr_matrix) :: CTP_block ! The current CTP dimension

      integer :: dim, dim2

      call zero(CMC)
      
      if(inverse_mass%diagonal) then

         do dim=1, blocks(CT,2)
            
            ! This is a borrowed reference so no need to deallocate.
            CT_block=block(CT, 1, dim)
            CTP_block=block(CTP, 1, dim)
            
            CM=matmul_T(CTP_block, block(inverse_mass,dim,dim), CT%sparsity)
            
            CMC_tmp=matmul_T(CM, CT_block, model=CMC%sparsity)
            
            call addto(CMC, CMC_tmp)
            
            call deallocate(CM)
            call deallocate(CMC_tmp)

         end do

      else

         do dim=1, blocks(CT,2)
            do dim2 = 1, blocks(CT,2)
               
               ! This is a borrowed reference so no need to deallocate.
               CT_block=block(CT, 1, dim)
               CTP_block=block(CTP, 1, dim2)
               
               CM=matmul_T(CTP_block, block(inverse_mass,dim2,dim), CT%sparsity)
               
               CMC_tmp=matmul_T(CM, CT_block, model=CMC%sparsity)
               
               call addto(CMC, CMC_tmp)
               
               call deallocate(CM)
               call deallocate(CMC_tmp)
            end do
         end do

      end if

      ewrite_minmax(cmc%val)

    end subroutine assemble_cmc_dg

    subroutine assemble_masslumped_cmc(cmc_m, ctp_m, inverse_masslump, ct_m)
      !!< Assemble the pressure matrix C_P^T M_l^{-1} C.
      !!< This currently does not support rotations.

      ! this subroutine is designed to start to replace cmc_wrapper and getcmc

      type(csr_matrix), intent(inout) :: cmc_m
      type(block_csr_matrix), intent(in) :: ctp_m, ct_m
      type(vector_field), intent(in) :: inverse_masslump

      ewrite(1,*) 'Entering assemble_masslumped_cmc'

      call mult_div_vector_div_T(cmc_m, ctp_m, inverse_masslump, ct_m)
      
      ewrite_minmax(cmc_m%val)

    end subroutine assemble_masslumped_cmc

    subroutine assemble_diagonal_schur(schur_diagonal_matrix,u,inner_m,ctp_m,ct_m)
      !!< Assemble the matrix C_P^T * [(Big_m)_diagonal]^-1 * C. 
      !!< This is used as a preconditioner for the full projection solve
      !!< when using the full momentum matrix.

      ! Fluidity velocity vector:
      type(vector_field), intent(in):: u
      ! Divergence matrices:
      type(block_csr_matrix), intent(in) :: ctp_m, ct_m
      ! Inner matrix of Schur complement:
      type(petsc_csr_matrix), intent(inout) :: inner_m   
      ! Product matrix:
      type(csr_matrix), intent(inout):: schur_diagonal_matrix
      ! Diagonal of inner_m - required to set up preconditioner matrix for Stokes problems
      ! (i.e. DiagonalSchurComplement):
      type(vector_field) :: inner_m_diagonal   

      ewrite(1,*) 'Entering assemble_diagonal_schur'

      call zero(schur_diagonal_matrix)
      call allocate(inner_m_diagonal,u%dim,u%mesh,"Diagonal_inner_m")
      call zero(inner_m_diagonal)
      call extract_diagonal(inner_m,inner_m_diagonal)
      call mult_div_invvector_div_T(schur_diagonal_matrix, ctp_m, inner_m_diagonal, ct_m)
      ewrite_minmax(schur_diagonal_matrix%val)
      call deallocate(inner_m_diagonal)

    end subroutine assemble_diagonal_schur

    subroutine repair_stiff_nodes(cmc_m, stiff_nodes_list)
    
      type(csr_matrix), intent(inout) :: cmc_m
      type(ilist), intent(inout) :: stiff_nodes_list
      
      integer :: row, col
      real, pointer :: row_diag
      integer, dimension(:), pointer :: row_m
      real, dimension(:), pointer :: row_val
      real :: tolerance
      
      ewrite(1,*) 'in repair_stiff_nodes()'
      
      tolerance = maxval(cmc_m%val)*epsilon(0.0)
      ewrite(2,*) 'tolerance = ', tolerance
      
      call flush_list(stiff_nodes_list)
      
      do row = 1, size(cmc_m, 1)
        row_diag=>diag_val_ptr(cmc_m, row)
        row_m=>row_m_ptr(cmc_m, row)
        row_val=>row_val_ptr(cmc_m, row)
        if(row_diag<tolerance) then
          ewrite(2,*) 'before, row, row_diag, sum(row_val) = ', row, row_diag, sum(abs(row_val))
          where(cmc_m%sparsity%colm==row) cmc_m%val = 0.0
          call zero_row(cmc_m, row)
          call addto_diag(cmc_m, row, 1.0)
          call insert(stiff_nodes_list, row)
        end if
      end do
      
      call print_list(stiff_nodes_list, 2)

    end subroutine repair_stiff_nodes
      
    subroutine zero_stiff_nodes(rhs, stiff_nodes_list)
    
      type(scalar_field), intent(inout) :: rhs
      type(ilist), intent(in) :: stiff_nodes_list
      
      integer :: row
      real, pointer :: row_diag
      real, dimension(:), pointer :: row_val
      real :: tolerance
      
      ewrite(1,*) 'in zero_stiff_nodes()'
      
      if(stiff_nodes_list%length>0) then
        ewrite(2,*) 'before node_val = ', node_val(rhs, list2vector(stiff_nodes_list))
        call set(rhs, list2vector(stiff_nodes_list), spread(0.0, 1, stiff_nodes_list%length))
      end if

    end subroutine zero_stiff_nodes
    
    subroutine assemble_masslumped_ctm(ctm_m, ctp_m, masslump)
      !!< Assemble the matrix C_P^T M_l^{-1}
      !!< This currently does not support rotations.

      type(block_csr_matrix), intent(in) :: ctp_m
      type(scalar_field), intent(in) :: masslump
      type(block_csr_matrix), intent(inout) :: ctm_m

      type(csr_matrix) :: lctm_m_block

      integer :: dim, row
      real, dimension(:), pointer :: row_val
      integer, dimension(:), pointer :: row_indices

      ewrite(1,*) 'Entering assemble_masslumped_ctm'

      call zero(ctm_m)

      do dim = 1, ctm_m%blocks(2)

        lctm_m_block = block(ctm_m, 1, dim)

        do row = 1, size(ctp_m, 1)
          row_indices=>row_m_ptr(ctp_m, row)
          row_val=>row_val_ptr(ctp_m, 1, dim, row)
          call set(lctm_m_block, (/row/), row_indices, &
                  spread((row_val/node_val(masslump, row_indices)), 1, 1))
        end do
        
        ewrite_minmax(ctm_m%val(1,dim)%ptr)

      end do

    end subroutine assemble_masslumped_ctm
    
  SUBROUTINE GETCMC(NDPSET,D3,CMC,C1T,C2T,C3T,&
       ML1,ML2,ML3,DM1,DM2,DM3,&
       NBUOY,BOY_ML, &
       NONODS,FREDOP,NCT,NCMC,&
       FINCMC,COLCMC,MIDCMC,&
       FINDCT,COLCT, &
                                ! RMRTIN is only used if ROTAT it has length =3*NONODS in2-D 
                                ! and 6*NONODS in 3-D.  
       RMRTIN,NRMRT,       &
                                ! Th following is for rotations only ROTAT=.TRUE.
       ROTAT,NNODRO,  NODROT,&
       NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z&
                                !              compressibility---------------------------------------------
       ,MKCOMP, C1TP, C2TP, C3TP) 
    INTEGER NDPSET,NCT,NCMC,NRMRT
    INTEGER NONODS,FREDOP
    REAL CTC
    REAL ML1(NONODS),ML2(NONODS),ML3(NONODS)
    INTEGER NBUOY
    REAL BOY_ML(NBUOY*9*NONODS)
    ! If ROTAT then it is assumed that C1T etc already conbtain the 
    ! rotation matrices via C1T=R^T*C1T
    REAL C1T(NCT),C2T(NCT),C3T(NCT),CMC(NCMC)
    ! TEMCMC is working space only used when we are performing a 
    ! multi-phase simulation. 
    ! Put volume fraction of this phase into RMEM(1) & make a no not 
    ! too close to 0.0 if not multi-phase then set volume frac=1.
    REAL DM1(NONODS),DM2(NONODS),DM3(NONODS)
    REAL RMRTIN(NRMRT)
    REAL AGI,BGI,CGI,DGI,EGI,FGI,GGI,HGI,KGI
    REAL DETJ
    INTEGER ROW,START,FINISH,ST2,FI2
    INTEGER FINDCT(FREDOP+1),COLCT(NCT)
    INTEGER FINCMC(FREDOP+1),COLCMC(NCMC)
    INTEGER MIDCMC(FREDOP)
    INTEGER ACR,COL,ST1,FI1,I,J,IBEG,IBEG2
    LOGICAL D3,LD3
    REAL &
         A11, A12, A13,&
         A21, A22, A23,&
         A31, A32, A33
    ! NB The normal is pointing out of the domain. 
    !  The rotation matrix in 3-D is R=  
    !   T1X   T1Y   T1Z
    !   T2X   T2Y   T2Z
    !   NX    NY    NZ
    !  The rotation matrix in 2-D is R=  
    !   T1X   T1Y   
    !   NX    NY     
    ! 
    INTEGER NNODRO
    REAL NX(NNODRO),NY(NNODRO),NZ(NNODRO)
    REAL T1X(NNODRO),T1Y(NNODRO),T1Z(NNODRO)
    REAL T2X(NNODRO),T2Y(NNODRO),T2Z(NNODRO)
    INTEGER NODROT(NNODRO)
    LOGICAL ROTAT
    INTEGER IDIM,II,NOD
    REAL RA,RB,RC,RD,RDET
    REAL RINA11,RINA12,RINA21,RINA22
    REAL RINA31,RINA32,RINA13,RINA23,RINA33
    INTEGER ICOL,ICOUNT

    !        compressibility---------------------------------------------
    REAL C1TP(NCT), C2TP(NCT), C3TP(NCT)
    INTEGER MKCOMP
    !        --------------------------------------added by crgw 18/07/06
    real, allocatable, dimension(:)::BOY_INV

    ! A COLOURIN TECHNIQUE WORK BE MORE EFFICIENT 
    ! HERE(SEE BOTTON OF SUB). 
    ! The following forms the matrix C1T RML C1 + C2T RML C2.
    ! Often there is only one diagonal so that ML=ML1,ML2,ML3. 
    ! If D3 then suppose the flow is 3-D else suppose it is 2-D. 
    ! May have to alter this to make it more efficient. 
    ! If a multi-phase solution then accumulate the 
    ! pressure matrix

    ewrite(1,*)'JUST INSIDE GETCMC'  
    !        ewrite(3,*)'ML1=',ML1 
    !        ewrite(3,*)'DM1=',DM1
    !
    do  I=1,NCMC! Was loop 23
       CMC(I)=0.
    end do ! Was loop 23
    !
    !
    !
    IF(NBUOY.NE.0) THEN
       ALLOCATE(BOY_INV(9*NONODS))
       DO NOD=1,NONODS
          ! Form A^{-1} 
          ! solve for inverse and put in BOY_INV:
          CALL INVXXX( &
               BOY_ML(NOD)+DM1(NOD),BOY_ML(NOD+NONODS),           BOY_ML(NOD+2*NONODS), &
               BOY_ML(NOD+3*NONODS),BOY_ML(NOD+4*NONODS)+DM2(NOD),BOY_ML(NOD+5*NONODS), &
               BOY_ML(NOD+6*NONODS),BOY_ML(NOD+7*NONODS),         BOY_ML(NOD+8*NONODS)+DM3(NOD),&
               BOY_INV(NOD),         BOY_INV(NOD+NONODS  ),BOY_INV(NOD+2*NONODS), &
               BOY_INV(NOD+3*NONODS),BOY_INV(NOD+4*NONODS),BOY_INV(NOD+5*NONODS), &
               BOY_INV(NOD+6*NONODS),BOY_INV(NOD+7*NONODS),BOY_INV(NOD+8*NONODS)   )
       END DO
    ELSE IF(ROTAT) THEN
       ! Form (R M_L R^T+D)^{-1} = RMRTIN
       ! NB in 2-D it has block numbering ...
       ! 1 *
       ! 2 3
       ! and in 3-D ...
       ! 1 * *
       ! 2 3 *
       ! 4 5 6       *-not stored
       !
       IF(.NOT.D3) THEN
          ! For 2-D ...
          IDIM=2
          LD3=.FALSE.
          do  I=1,3*NONODS! Was loop 19
             RMRTIN(I) =0.
          end do ! Was loop 19
          do  I=1,NONODS! Was loop 21
             RMRTIN(I)         =1.0/(ML1(I)+DM1(I))
             RMRTIN(I+2*NONODS)=1.0/(ML2(I)+DM2(I))
          end do ! Was loop 21
          !
          do  II=1,NNODRO! Was loop 31
             NOD=NODROT(II)
             ! Form R*(M_L R^T)
             CALL SMLMA2(&
                  A11, A12, &
                  A21, A22, &
                                ! =
                  T1X (II),  T1Y(II) ,&
                  NX(II)  ,   NY(II),    &
                                ! *
                  T1X(II)*ML1(NOD), NX(II)*ML1(NOD),&
                  T1Y(II)*ML2(NOD), NY(II)*ML2(NOD)   )
             ! add in D
             A11=A11+DM1(NOD)
             A22=A22+DM2(NOD)
             !
             ! Form A^{-1} 
             RA=A11
             RB=A12
             RC=A21
             RD=A22
             RDET=RA*RD-RB*RC
             RINA11=RD/RDET
             RINA12=-RB/RDET
             RINA21=-RC/RDET
             RINA22=RA/RDET
             ! set values of (RMR^T +D)^{-1}
             RMRTIN(NOD)         =RINA11
             RMRTIN(NOD+NONODS  )=RINA21
             RMRTIN(NOD+2*NONODS)=RINA22
          end do ! Was loop 31
          !
          ! IF(D3.EQ.0)
       ELSE
          ! For 3-D ...

          IDIM=3
          LD3=.TRUE.
          do  I=1,6*NONODS! Was loop 191
             RMRTIN(I) =0.
          end do ! Was loop 191
          do  I=1,NONODS! Was loop 22
             RMRTIN(I)         =1.0/(ML1(I)+DM1(I))
             RMRTIN(I+2*NONODS)=1.0/(ML2(I)+DM2(I))
             RMRTIN(I+5*NONODS)=1.0/(ML3(I)+DM3(I))
          end do ! Was loop 22
          !
          do  II=1,NNODRO! Was loop 30
             NOD=NODROT(II)
             ! Form R *(M_L R^T)
             CALL SMLMA3(&
                  A11, A12, A13,&
                  A21, A22, A23,&
                  A31, A32, A33,&
                                ! =
                  T1X(II),   T1Y(II),   T1Z(II),&
                  T2X(II),   T2Y(II),   T2Z(II),&
                  NX(II),    NY(II),    NZ(II),&
                                ! *
                  T1X(II)*ML1(NOD), T2X(II)*ML1(NOD), NX(II)*ML1(NOD),&
                  T1Y(II)*ML2(NOD), T2Y(II)*ML2(NOD), NY(II)*ML2(NOD),&
                  T1Z(II)*ML3(NOD), T2Z(II)*ML3(NOD), NZ(II)*ML3(NOD)  )
             ! add in D
             A11=A11+DM1(NOD)
             A22=A22+DM2(NOD)
             A33=A33+DM3(NOD)
             !
             ! Form A^{-1} 
             AGI=A11
             BGI=A12
             CGI=A13
             !
             DGI=A21
             EGI=A22
             FGI=A23
             !
             GGI=A31
             HGI=A32
             KGI=A33

             DETJ=AGI*(EGI*KGI-FGI*HGI)&
                  -BGI*(DGI*KGI-FGI*GGI)&
                  +CGI*(DGI*HGI-EGI*GGI)
             RINA11= (EGI*KGI-FGI*HGI) /DETJ
             RINA21=-(DGI*KGI-FGI*GGI) /DETJ
             RINA31= (DGI*HGI-EGI*GGI) /DETJ
             !
             RINA12=-(BGI*KGI-CGI*HGI) /DETJ
             RINA22= (AGI*KGI-CGI*GGI) /DETJ
             RINA32=-(AGI*HGI-BGI*GGI) /DETJ
             !
             RINA13= (BGI*FGI-CGI*EGI) /DETJ
             RINA23=-(AGI*FGI-CGI*DGI) /DETJ
             RINA33= (AGI*EGI-BGI*DGI) /DETJ
             ! set values of (RMR^T +D)^{-1}
             RMRTIN(NOD)         =RINA11
             RMRTIN(NOD+NONODS  )=RINA21
             RMRTIN(NOD+2*NONODS)=RINA22
             RMRTIN(NOD+3*NONODS)=RINA31
             RMRTIN(NOD+4*NONODS)=RINA32
             RMRTIN(NOD+5*NONODS)=RINA33
          end do ! Was loop 30
          !
          ! ELSE IF(D3.EQ.0)
       ENDIF
       !
       ! ENOF IF(ROTAT)
    ENDIF
    !
    !
    ! NOW FORM CMC
    IF(.NOT.D3) THEN
       !
       do  ROW=1,FREDOP! Was loop 303
          !         ewrite(3,*)'row,ml1(row):',row,ml1(row)
          START=FINCMC(ROW)
          FINISH=FINCMC(ROW+1)-1
          ST2=FINDCT(ROW)
          FI2=FINDCT(ROW+1)-1
          do  ACR=START,FINISH! Was loop 313
             COL=COLCMC(ACR)
             CTC=0.
             ! We wish to mult vector column 'ROW' of C1 by vector
             ! column 'COL' of C1.
             ST1=FINDCT(COL)
             FI1=FINDCT(COL+1)-1
             do  I=ST1,FI1! Was loop 323
                do  J=ST2,FI2! Was loop 333
                   IF(COLCT(I).EQ.COLCT(J)) THEN
                      ICOL=COLCT(I)
                      !
                      IF(.NOT.ROTAT) THEN
                         !        compressibility---------------------------------------------
                         IF(MKCOMP.EQ.3) THEN
                            CTC=CTC+C1T(I)*C1TP(J)/(ML1(ICOL)+DM1(ICOL)  )&
                                 +C2T(I)*C2TP(J)/(ML2(ICOL)+DM2(ICOL)  )
                         ELSE                     
                            !        --------------------------------------added by crgw 18/07/06
                            CTC=CTC+C1T(I)*C1T(J)/(ML1(ICOL)+DM1(ICOL)  )&
                                 +C2T(I)*C2T(J)/(ML2(ICOL)+DM2(ICOL)  )
                         ENDIF
                      ELSE
                         ! NB RMRTIN contains (R M R^T +D)^{-1}
                         RINA11=RMRTIN(ICOL) 
                         RINA12=RMRTIN(ICOL+NONODS) 
                         RINA21=RMRTIN(ICOL+NONODS)
                         RINA22=RMRTIN(ICOL+2*NONODS)
                         !
                         !        compressibility---------------------------------------------
                         IF(MKCOMP.EQ.3) THEN
                            CTC=CTC+C1T(I)*(C1TP(J)*RINA11+C2TP(J)*RINA12)&
                                 +C2T(I)*(C1TP(J)*RINA21+C2TP(J)*RINA22)
                         ELSE                     
                            !        --------------------------------------added by crgw 18/07/06
                            CTC=CTC+C1T(I)*(C1T(J)*RINA11+C2T(J)*RINA12)&
                                 +C2T(I)*(C1T(J)*RINA21+C2T(J)*RINA22)
                         ENDIF
                      ENDIF
                   ENDIF
                end do ! Was loop 333
             end do ! Was loop 323
             CMC(ACR)=CTC
          end do ! Was loop 313
       end do ! Was loop 303
    ENDIF
    !
    ewrite(2,*) 'MKCOMP', MKCOMP
    IF(D3) THEN
       do  ROW=1,FREDOP! Was loop 301
          START=FINCMC(ROW)
          FINISH=FINCMC(ROW+1)-1
          ST2=FINDCT(ROW)
          FI2=FINDCT(ROW+1)-1
          do  ACR=START,FINISH! Was loop 311
             COL=COLCMC(ACR)
             CTC=0.
             ! We wish to mult vector column 'ROW' of C1 by vector
             ! column 'COL' of C1.
             ST1=FINDCT(COL)
             FI1=FINDCT(COL+1)-1
             IBEG2=0
             do  I=ST1,FI1! Was loop 321
                IBEG=IBEG2
                do  J=ST2+IBEG,FI2! Was loop 331
                   IF(COLCT(I).EQ.COLCT(J)) THEN
                      ICOL=COLCT(I)
                      IBEG2=J-ST2+1
                      !
                      IF(NBUOY.NE.0) THEN
                         !
                         RINA11=BOY_INV(ICOL)         
                         RINA12=BOY_INV(ICOL+NONODS  )
                         RINA13=BOY_INV(ICOL+2*NONODS)

                         RINA21=BOY_INV(ICOL+3*NONODS)
                         RINA22=BOY_INV(ICOL+4*NONODS)
                         RINA23=BOY_INV(ICOL+5*NONODS)

                         RINA31=BOY_INV(ICOL+6*NONODS)
                         RINA32=BOY_INV(ICOL+7*NONODS)
                         RINA33=BOY_INV(ICOL+8*NONODS)

                         CTC=CTC &
                              +C1T(I)*(C1T(J)*RINA11+C2T(J)*RINA12+C3T(J)*RINA13) &
                              +C2T(I)*(C1T(J)*RINA21+C2T(J)*RINA22+C3T(J)*RINA23) &
                              +C3T(I)*(C1T(J)*RINA31+C2T(J)*RINA32+C3T(J)*RINA33)
                      ELSE IF(.NOT.ROTAT) THEN
                         !        compressibility---------------------------------------------
                         IF(MKCOMP.EQ.3) THEN
                            CTC=CTC+C1T(I)*C1TP(J)/(ML1(ICOL)+DM1(ICOL)  )&
                                 +C2T(I)*C2TP(J)/(ML2(ICOL)+DM2(ICOL)  )&
                                 +C3T(I)*C3TP(J)/(ML3(ICOL)+DM3(ICOL)  )

                         ELSE                     
                            !        --------------------------------------added by crgw 18/07/06
                            CTC=CTC+C1T(I)*C1T(J)/(ML1(ICOL)+DM1(ICOL)  )&
                                 +C2T(I)*C2T(J)/(ML2(ICOL)+DM2(ICOL)  )&
                                 +C3T(I)*C3T(J)/(ML3(ICOL)+DM3(ICOL)  )
                         ENDIF
                      ELSE
                         ! NB RMRTIN contains (R M R^T +D)^{-1}
                         !
                         RINA11= RMRTIN(ICOL)
                         !
                         RINA21= RMRTIN(ICOL+NONODS)
                         RINA22= RMRTIN(ICOL+2*NONODS)
                         !
                         RINA31= RMRTIN(ICOL+3*NONODS)
                         RINA32= RMRTIN(ICOL+4*NONODS)
                         RINA33= RMRTIN(ICOL+5*NONODS)
                         !
                         RINA12=RINA21
                         RINA13=RINA31
                         !
                         RINA23= RINA32
                         !        compressibility---------------------------------------------
                         IF(MKCOMP.EQ.3) THEN
                            CTC=CTC&
                                 +C1T(I)*(C1TP(J)*RINA11+C2TP(J)*RINA12+C3TP(J)*RINA13)&
                                 +C2T(I)*(C1TP(J)*RINA21+C2TP(J)*RINA22+C3TP(J)*RINA23)&
                                 +C3T(I)*(C1TP(J)*RINA31+C2TP(J)*RINA32+C3TP(J)*RINA33)
                         ELSE
                            !        --------------------------------------added by crgw 18/07/06
                            CTC=CTC&
                                 +C1T(I)*(C1T(J)*RINA11+C2T(J)*RINA12+C3T(J)*RINA13)&
                                 +C2T(I)*(C1T(J)*RINA21+C2T(J)*RINA22+C3T(J)*RINA23)&
                                 +C3T(I)*(C1T(J)*RINA31+C2T(J)*RINA32+C3T(J)*RINA33)
                         ENDIF

                      ENDIF
                      exit
                   ENDIF
                end do ! Was loop 331
             end do ! Was loop 321
             CMC(ACR)=CTC
          end do ! Was loop 311
       end do ! Was loop 301
    ENDIF
    !
    IF((NDPSET.GE.1).AND.(NDPSET.LE.FREDOP)) THEN 
       ! SET PRSEEURE NODE NDPSET TO ZERO.
       !        ewrite(3,*)'CMC(MIDCMC(NDPSET))=',CMC(MIDCMC(NDPSET))
       CMC(MIDCMC(NDPSET))=INFINITY
    ENDIF
    !
    !         i=83
    !         do icount=findct(i),findct(i+1)-1
    !        ewrite(3,*) 'c1t,c2t,c3t:',c1t(icount),c2t(icount),c3t(icount)
    !         end do
    !
    do  I=1,-FREDOP! Was loop 314
       ewrite(3,*)'I,COLM:',I,(COLCMC(ICOUNT),&
            ICOUNT=FINCMC(I),FINCMC(I+1)-1)
       ewrite(3,*)'cmc(midcmc(i))',cmc(midcmc(i))
       ewrite(3,*)'I,CMC:',I,(CMC(ICOUNT),&
            ICOUNT=FINCMC(I),FINCMC(I+1)-1)
    end do ! Was loop 314
    !          stop

    !
    !
    ! THE COLOURING ALGORITHM IS AS FOLLOWS. ...
    !        DO 10 COLOR=1,MXCOLO
    !             DO 20 I=1,FREDOP
    !                 V(I)=0.
    !                 IF(NODCOL(I).EQ.COLOR) V(I)=1.
    !20          CONTINUE
    !C PERFORM MATRIX VECTOR MULTIPLICATION 
    !C IE   x= CT ML C v 
    !             DO 30 I=1,FREDOP
    !             DO 30 COUNT=FINCMC(I),FINCMC(I+1)-1
    !                 IF(NODCOL(COLCMC(COUNT).EQ.COLOR) CMC(COUNT)=X(I)
    !30          CONTINUE
    !10     CONTINUE
    !       ewrite(3,*)'ML1=',ML1
    ewrite(1,*)'exiting getcmc'
    RETURN
  END SUBROUTINE GETCMC
  
  ! This sub calculates KCMC which filters out the singularities in 
  ! the pressure eqns.
  ! Also add KCMC into CMC. 
  SUBROUTINE CMCFILT(X,Y,Z, &
       KCMC,CMC,NDPSET,D3,&
       FINCMC,COLCMC,MIDCMC,NCMC,FREDOP,TOTELE,&
       N,NLX,NLY,NLZ, &
       M,MLX,MLY,MLZ, WEIGHT,NGI,NLOC,MLOC, &
       PNDGLN,XONDGL,XNONOD,&
       ISPHERE)
    INTEGER XNONOD
    INTEGER MLOC,NLOC,NGI
    INTEGER TOTELE,FREDOP,NCMC
    INTEGER NDPSET,ISPHERE
    LOGICAL D3
    ! If PHI=0.5 then Crank-Nicolson is used in Navier Stokes equ's.
    ! If PHI=1.0 then bakward Euler is used in NS equns.
    ! Similarly for the temperature equation except width variable THETA.
    ! have gravity in -ve y direction.
    ! ISPHERE indicates the order of the superparamtertic FEM mapping
    !
    REAL KCMC(NCMC),CMC(NCMC)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    INTEGER PNDGLN(TOTELE*MLOC)
    INTEGER XONDGL(TOTELE*NLOC)
    INTEGER FINCMC(FREDOP+1),COLCMC(NCMC),MIDCMC(FREDOP)

    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL M(MLOC,NGI),MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
    REAL WEIGHT(NGI)

    ! Local variables...
    INTEGER NSCMC
    INTEGER, ALLOCATABLE, DIMENSION(:)::SFINCMC,SCOLCMC
    REAL :: UDL(NLOC*NLOC,NGI),VDL(NLOC*NLOC,NGI),WDL(NLOC*NLOC,NGI)
    REAL :: GAMMA(NLOC*NLOC,NGI), PML(FREDOP)
    REAL, ALLOCATABLE :: IZER(:),SCMC(:)
    INTEGER DumbInt(1)
    REAL DUMMYARRAY0(0)
    REAL :: DETWEI(NGI), MX(MLOC,NGI), MY(MLOC,NGI), MZ(MLOC,NGI)
    ewrite(3, *) "SUBROUTINE CMCFILT()"
    ! Calculate stripped down sparcity pattern for diffusion 
    ! matrix in pressure filter...
    ! Use KCMC as work space here...

    ALLOCATE( SFINCMC(FREDOP+1) )
    ALLOCATE( SCOLCMC(NCMC) )
    !
    ! Trime down CMC sparcity ie) abtain SFINCMC,SCOLCMC
    ewrite(1, *) "going into LIMITCMC"
    CALL LIMITCMC(SFINCMC, SCOLCMC, NSCMC, NCMC, &
         FINCMC,COLCMC,NCMC,FREDOP,&
         TOTELE,MLOC,PNDGLN, KCMC)

    ewrite(3, *) "FINISHED LIMITCMC"

    ewrite(2,*) 'NSCMC, MLOC, NLOC', NSCMC, MLOC, NLOC

    ALLOCATE(SCMC(NSCMC))

    ! Calculate diffusion matrix for filter
    if(.true.) then
  !  if(isphere.ne.2) then
    ! original method...
    CALL CTDIFF(X,Y,Z, &
         SCMC,PML,&
         SFINCMC, SCOLCMC,NSCMC,FREDOP,TOTELE, &
         N,NLX,NLY,NLZ, &
         M,MLX,MLY,MLZ, WEIGHT,NGI,NLOC,MLOC, &
         PNDGLN,XONDGL,XNONOD,&
         DETWEI,&
         UDL,VDL,WDL,&
         GAMMA,&
         MX,MY,MZ,ISPHERE    )
     else
     ! new method...
    CALL CTDIFFSPH(X,Y,Z, &
         SCMC,PML,&
         SFINCMC, SCOLCMC,NSCMC,FREDOP,TOTELE, &
         N,NLX,NLY,NLZ, &
         WEIGHT,NGI,NLOC,MLOC, &
         PNDGLN,XONDGL,XNONOD,&
         DETWEI,&
         ISPHERE    )
     endif

    ! P1 & Q1 basis...
    ALLOCATE(IZER(NSCMC))
    IZER = 0.0
    ewrite(1, *) "going inot GETCMC"
    ! Calculate the actual filter..
    ! RMRTIN is only used if ROTAT it has length =3*NONODS in2-D 
    ! and 6*NONODS in 3-D.  -NRTDR has the same length  
    ! rotations only ROTAT=.TRUE.
    CALL GETCMC(NDPSET,D3,KCMC,SCMC,IZER,IZER, &
         PML,PML,PML,IZER,IZER,IZER,&
         0,IZER,&
         FREDOP,FREDOP,NSCMC,NCMC,&
         FINCMC,COLCMC,MIDCMC,&
         SFINCMC,SCOLCMC,&
         DUMMYARRAY0,0,        &
         .FALSE.,0,  DumbInt, &
         DUMMYARRAY0, DUMMYARRAY0, DUMMYARRAY0, DUMMYARRAY0, DUMMYARRAY0, DUMMYARRAY0, & 
         DUMMYARRAY0, DUMMYARRAY0, DUMMYARRAY0 &
  !      compressibility-added for consistency--------------------------------------
         , 0, IZER, IZER, IZER &
  !      --------------------------------------added by crgw 18/07/06
         )

    ewrite(3, *) "JUSTR OUT SIDE GETCMC"
    ! CMC=CMC+KCMC
    !              call rclear(kcmc,ncmc)
    CALL RADDIN(CMC,KCMC,NCMC) 

    DEALLOCATE( SFINCMC )
    DEALLOCATE( SCOLCMC )
    DEALLOCATE( IZER, SCMC )
    ewrite(3, *) "END SUBROUTINE CMCFILT"
    RETURN
  END SUBROUTINE CMCFILT
  
  SUBROUTINE CTDIFFSPH(X,Y,Z,&
     &     SCMC,PML,&
     &     SFINCMC,SCOLCMC,NSCMC,FREDOP,TOTELE, &
     &     N,NLX,NLY,NLZ, &
     &     WEIGHT,NGI,NLOC,MLOC, &
     &     PNDGLN,XONDGL,XNONOD,&
     &     DETWEI,&
     &     ISPHERE     )
!     This sub calculates the stabalisation diffusion matrix SCMC. 
!     Eventuall KCMC=SCMC^T PML^-1 SCMC
!     PNDGLN=element pter for unknown pressures. 
!     XONDGL=element pter for coordinates(x,y,z)
!     -------------------------------------------------------------------------------------------
!     If we are solving for another variable like temperature 
!     or chemical species then NU,NV,NW will be the velocities 
!     and U=TEMPERATURE OR CHEMICAL SPECIES. 
!     -------------------------------------------------------------------------------------------
!     NB SOURCE, DENPT, viscosity are SNONOD in length and have ele-pt SONDGL.
! ISPHERE indicates the order of the superparamtertic FEM mapping
      INTEGER XNONOD
      INTEGER MLOC,NLOC,NGI
      INTEGER TOTELE,FREDOP,NSCMC,ISPHERE
!     
      REAL SCMC(NSCMC),PML(FREDOP)
      REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
      INTEGER PNDGLN(TOTELE*MLOC)
      INTEGER XONDGL(TOTELE*NLOC)
      INTEGER SFINCMC(FREDOP+1),SCOLCMC(NSCMC)
!     
      REAL, target:: N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
      REAL WEIGHT(NGI)
!     
      REAL DETWEI(NGI)
!     HX,HY-characteristic length scales in x,y directions.
!     Local variables...
!     
      INTEGER ELE,ILOC,JLOC
      INTEGER GI,GLOBI,GLOBJ
      INTEGER POSMAT
!     
      REAL MATRIX,PMASS
      REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
      REAL TENSXX(NGI),TENSXY(NGI),TENSXZ(NGI),TENSYY(NGI),TENSYZ(NGI),TENSZZ(NGI)
      REAL UDLGI,VDLGI,WDLGI,RN,SQRTRN,HDIST
       REAL A11(NGI),A12(NGI),A13(NGI)
       REAL A21(NGI),A22(NGI),A23(NGI)
       REAL A31(NGI),A32(NGI),A33(NGI) 
       REAL XD(NGI),YD(NGI),ZD(NGI)
! Local coords...
         INTEGER NCLOC
         REAL, ALLOCATABLE, DIMENSION(:,:)::NC
         REAL, ALLOCATABLE, DIMENSION(:,:)::NCLX
         REAL, ALLOCATABLE, DIMENSION(:,:)::NCLY
         REAL, ALLOCATABLE, DIMENSION(:,:)::NCLZ
         REAL, ALLOCATABLE, DIMENSION(:,:)::NCX
         REAL, ALLOCATABLE, DIMENSION(:,:)::NCY
         REAL, ALLOCATABLE, DIMENSION(:,:)::NCZ

      ewrite(1, *) "SUBROUTINE CTDIFFSPH()"
!     
!     
      SCMC(1:NSCMC) = 0.0
      PML(1:FREDOP) = 0.0
!     
!     This subroutine forms a contabution to the Right Hand Side
!     of Poissons pressure equation, as well as  F1 & F2.
!
      NCLOC=4
      IF(ISPHERE.EQ.2) NCLOC=10

      ALLOCATE(NC(NCLOC,NGI))
      ALLOCATE(NCLX(NCLOC,NGI))
      ALLOCATE(NCLY(NCLOC,NGI))
      ALLOCATE(NCLZ(NCLOC,NGI))
      ALLOCATE(NCX(NCLOC,NGI))
      ALLOCATE(NCY(NCLOC,NGI))
      ALLOCATE(NCZ(NCLOC,NGI))
      
      IF(ISPHERE.NE.0) THEN
!         CALL TRIQUA(L1, L2, L3, L4, PGWEIG, .TRUE.,NGI)
! Work out the shape functions and there derivatives...
!         CALL SHATRI(L1, L2, L3, L4, PGWEIG, .TRUE., 
!     &            NCLOC,NGI,                          
!     &            NC,NCLX,NCLY,NCLZ) 
          ASSERT(ncloc==10)
          
          call QuadTetShapes(nc, nclx, ncly, nclz, ncloc, ngi)
      ELSE
          NC = N
          NCLX = NLX
          NCLY = NLY
          NCLZ = NLZ
      ENDIF
      ewrite(1,*) 'entering element loop'
      ewrite(2,*) 'nloc,mloc,ncloc,ngi:',nloc,mloc,ncloc,ngi
!     
      do  ELE=1,TOTELE! Was loop 340
      
! Get determinant and derivatives...

         CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NCLOC,NGI, &
     &               N,NLX,NLY,NLZ, NC,NCLX,NCLY,NCLZ, WEIGHT, DETWEI,.TRUE.,.FALSE., &
     &               NX,NY,NZ, NCX,NCY,NCZ,&
     &               A11,A12,A13, A21,A22,A23, A31,A32,A33,&
     &               XD,YD,ZD,&
     &               ISPHERE) 

      do  GI=1,NGI! Was loop 331
           TENSXX(GI)=0.0
           TENSXY(GI)=0.0
           TENSXZ(GI)=0.0
           TENSYY(GI)=0.0
           TENSYZ(GI)=0.0
           TENSZZ(GI)=0.0
!     nb we want L^2 - at the moment we have L^2 on the diagonal.
      do  ILOC=1,NLOC! Was loop 354
!     TENS=Q Q^T
                  RN=NX(ILOC,GI)**2+NY(ILOC,GI)**2+NZ(ILOC,GI)**2
                  SQRTRN=SQRT(RN)
                  HDIST=1./SQRTRN
                  UDLGI=NX(ILOC,GI)/SQRTRN
                  VDLGI=NY(ILOC,GI)/SQRTRN
                  WDLGI=NZ(ILOC,GI)/SQRTRN
                  
                  TENSXX(GI)=TENSXX(GI) + HDIST*UDLGI*UDLGI
                  TENSXY(GI)=TENSXY(GI) + HDIST*UDLGI*VDLGI
                  TENSXZ(GI)=TENSXZ(GI) + HDIST*UDLGI*WDLGI
                  TENSYY(GI)=TENSYY(GI) + HDIST*VDLGI*VDLGI
                  TENSYZ(GI)=TENSYZ(GI) + HDIST*VDLGI*WDLGI
                  TENSZZ(GI)=TENSZZ(GI) + HDIST*WDLGI*WDLGI
      end do ! Was loop 354
      end do ! Was loop 331
!     
      do  ILOC=1,MLOC! Was loop 350
            GLOBI=PNDGLN((ELE-1)*MLOC+ILOC)
      do  JLOC=1,MLOC! Was loop 360
                    GLOBJ=PNDGLN((ELE-1)*MLOC+JLOC)
!     
                    MATRIX=0.
                    PMASS =0. 
!     
      do  GI=1,NGI! Was loop 458
           MATRIX=MATRIX+(NX(ILOC,GI)*(NX(JLOC,GI)*TENSXX(GI)&
     &                 +NY(JLOC,GI)*TENSXY(GI)+NZ(JLOC,GI)*TENSXZ(GI))&
     &               +NY(ILOC,GI)*(NX(JLOC,GI)*TENSXY(GI)&
     &                 +NY(JLOC,GI)*TENSYY(GI)+NZ(JLOC,GI)*TENSYZ(GI))&
     &               +NZ(ILOC,GI)*(NX(JLOC,GI)*TENSXZ(GI)&
     &                 +NY(JLOC,GI)*TENSYZ(GI)+NZ(JLOC,GI)*TENSZZ(GI)) )*DETWEI(GI)
!     USE THE COMPONENT OF DIFLIN THE X,Y & Z-DIRECTIONS 
!     RESPECTIVELY FOR C1T,C2T,C3T.
                       PMASS =PMASS +N(ILOC,GI)*N(JLOC,GI)*DETWEI(GI) 
      end do ! Was loop 458
!     
!     we will now place contributions into the matrices SCMC and BIGT
!     these are can be symm matrices.
!     Find POSMAT. ***************************************
!     
                    CALL POSINMAT(POSMAT,GLOBI,GLOBJ,&
     &                   FREDOP,SFINCMC,SCOLCMC,NSCMC)
!     Find POSMAT. ***************************************
!     
                    PML(GLOBI)=PML(GLOBI)+PMASS
                    SCMC(POSMAT)=SCMC(POSMAT)+MATRIX 
!     
!     
      end do ! Was loop 360
      end do ! Was loop 350
      end do ! Was loop 340
         ewrite(1, *) "about to exit CTDIFFSPH()"

  end subroutine ctdiffsph
  
  SUBROUTINE LIMITCMC(SFINCMC,SCOLCMC,NSCMC,MXNSCMC,&
       &         FINCMC,COLCMC,NCMC,FREDOP,&
       &         TOTELE,MLOC,PNDGLN, KCMC)
    ! Trime down CMC sparcity ie) abtain SFINCMC,SCOLCMC
    ! KCMC is used as tempory storage here.
    INTEGER NSCMC,NCMC,MXNSCMC
    INTEGER FREDOP,TOTELE,MLOC,PNDGLN(TOTELE*MLOC)
    INTEGER SFINCMC(FREDOP+1),SCOLCMC(MXNSCMC)
    INTEGER FINCMC(FREDOP+1),COLCMC(NCMC)
    REAL KCMC(NCMC)
    ! Local variables...
    INTEGER ILOC,JLOC,PINOD,PJNOD,POSMAT,ELE,COUNT,COUNT2
    !
    KCMC(1:NCMC) = 0.0
    do ELE=1,TOTELE! Was loop 
       do ILOC=1,MLOC! Was loop 
          PINOD=PNDGLN((ELE-1)*MLOC+ILOC)
          do JLOC=1,MLOC! Was loop 
             PJNOD=PNDGLN((ELE-1)*MLOC+JLOC)
             CALL POSINMAT(POSMAT,PINOD,PJNOD,&
                  &                  FREDOP,FINCMC,COLCMC,NCMC)
             KCMC(POSMAT)=1.0
          END DO
       END DO
    END DO
    !
    COUNT2=1
    do PINOD=1,FREDOP! Was loop 
       SFINCMC(PINOD)=COUNT2
       do COUNT=FINCMC(PINOD),FINCMC(PINOD+1)-1! Was loop 
          IF(KCMC(COUNT).NE.0.0) THEN
             SCOLCMC(COUNT2)=COLCMC(COUNT)
             COUNT2=COUNT2+1
          END IF
       END DO
    END DO
    SFINCMC(FREDOP+1)=COUNT2
    NSCMC=COUNT2-1
    !       ewrite(3,*) 'scolcmc:',(scolcmc(i),i=1,50)
    RETURN
  END SUBROUTINE LIMITCMC
  
end module assemble_cmc
