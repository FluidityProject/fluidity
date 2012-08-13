
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
!    amcgsoftware@imperial.ac.uk
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

module solvers_module

  use fldebug
  use sparse_tools_petsc
  use solvers
  use fields
  use global_parameters, only: OPTION_PATH_LEN
  use spud

  implicit none
    
  private
  
  public  :: solver, PRES_DG_MULTIGRID
  
  interface solver
     module procedure solve_via_copy_to_petsc_csr_matrix
  end interface solver
  
contains

! -----------------------------------------------------------------------------

  subroutine solve_via_copy_to_petsc_csr_matrix(A, &
                                                x, &
                                                b, &
                                                findfe, &
                                                colfe, &
                                                option_path)
  
     !!< Solve a matrix Ax = b system via copying over to the
     !!< petsc_csr_matrix type and calling the femtools solver
     !!< using the spud options given by the field options path
     
    integer, dimension(:), intent(in) :: findfe 
    integer, dimension(:), intent(in) :: colfe
    real, dimension(:), intent(in) :: a,b
    real, dimension(:), intent(inout) :: x
    character(len=*), intent(in) :: option_path
    
    ! local variables
    integer :: i,j,k
    integer :: rows
    integer, dimension(:), allocatable :: dnnz
    type(petsc_csr_matrix) :: matrix
    type(scalar_field) :: rhs
    type(scalar_field) :: solution
    
    rows = size(x)
    
    ewrite(3,*) size(x)+1, size(findfe)

    assert(size(x) == size(b))
    assert(size(a) == size(colfe))
    assert((size(x)+1) == size(findfe))
    
    ! find the number of non zeros per row
    allocate(dnnz(rows)) ; dnnz = 0
    do i = 1,rows
       dnnz(i) = findfe(i+1) - findfe(i)
    end do
    
    ! allocate the petsc_csr_matrix using nnz (pass dnnz also for onnz) 
    call allocate(matrix, rows, rows, dnnz, dnnz, (/1,1/), name = 'dummy')
    call zero(matrix)
    
    ! add in the entries to petsc matrix
    do i = 1, rows
      do j = findfe(i), findfe(i+1) - 1
          k = colfe(j)
          call addto(matrix, blocki = 1, blockj = 1, i = i, j = k, val = a(j))
       end do
    end do 

    call assemble(matrix)    
    
    ! Set up rhs and initial guess which are scalar field types.        
    allocate(rhs%val(rows))
    allocate(solution%val(rows))
        
    do i = 1,rows
       call set(rhs, i, b(i))
       call set(solution, i,  x(i))
    end do
        
    ! solve matrix
    call petsc_solve(solution, &
                     matrix, &
                     rhs, &
                     option_path = trim(option_path))
    
    ! copy solution back    
    do i = 1,rows
       x(i) =  node_val(solution, i)
    end do
    
    ! deallocate as needed
    deallocate(dnnz)    
    deallocate(rhs%val)
    deallocate(solution%val)
    call deallocate(matrix)
      
  end subroutine solve_via_copy_to_petsc_csr_matrix




  SUBROUTINE PRES_DG_MULTIGRID(CMC, P, RHS, &
       NCOLCMC, CV_NONODS, FINDCMC, COLCMC, MIDCMC, &
       totele, cv_nloc, x_nonods, cv_ndgln, x_ndgln )
    !
    ! Solve CMC * P = RHS for RHS.
    ! form a discontinuous pressure mesh for pressure...
    implicit none
    INTEGER, intent( in ) ::  NCOLCMC, CV_NONODS, totele, cv_nloc, x_nonods
    REAL, DIMENSION( NCOLCMC ), intent( in ) ::  CMC
    REAL, DIMENSION( CV_NONODS ), intent( inout ) ::  P
    REAL, DIMENSION( CV_NONODS ), intent( in ) :: RHS
    INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCMC
    INTEGER, DIMENSION( NCOLCMC ), intent( in ) :: COLCMC
    INTEGER, DIMENSION( CV_NONODS ), intent( in ) :: MIDCMC
    INTEGER, DIMENSION( cv_nloc*totele ), intent( in ) :: cv_ndgln, x_ndgln

    REAL ERROR, RELAX, RELAX_DIAABS, RELAX_DIA
    INTEGER N_LIN_ITS, NGL_ITS

    PARAMETER(ERROR=1.E-15, RELAX=0.1, RELAX_DIAABS=0.0)
    PARAMETER(RELAX_DIA=2.0, N_LIN_ITS=2)
    PARAMETER(NGL_ITS=50) ! maybe we need to increase this...

    !  PARAMETER(ERROR=1.E-15, RELAX=0.05, RELAX_DIAABS=0.0)
    !  PARAMETER(RELAX_DIA=2.0, N_LIN_ITS=2)
    !  PARAMETER(NGL_ITS=2000)

    ! RELAX: overall relaxation coeff; =1 for no relaxation. 
    ! RELAX_DIAABS: relaxation of the absolute values of the sum of the row of the matrix;
    !                      - recommend >=2 for hard problems, =0 for easy
    ! RELAX_DIA: relaxation of diagonal; =1 no relaxation (normally applied). 
    ! N_LIN_ITS = no of linear iterations
    ! ERROR= solver tolerence between 2 consecutive iterations
    ! NGL_ITS = no of global its

    integer, dimension( : ), allocatable :: findcmc_small, colcmc_small, midcmc_small, &
         MAP_DG2CTY
    real, dimension( : ), allocatable :: cmc_small, resid_dg, resid_cty, &
         nods_sourou, DP_DG, DP_SMALL
    integer :: ele,cv_iloc, dg_nod, cty_nod, jcolcmc, jcolcmc_small
    integer :: mx_ncmc_small, ncmc_small, count, count2, count3, GL_ITS

    character(len=OPTION_PATH_LEN) :: path = "/tmp/pressure"
    integer :: stat

    ! obtain sparcity of a new matrix 
    mx_ncmc_small = ncolcmc * 4
    allocate( FINDcmc_small(x_nonods+1) )
    allocate( colcmc_small(mx_ncmc_small) )
    allocate( midcmc_small(x_nonods) )
    allocate( MAP_DG2CTY(cv_nonods) )

    EWRITE(3,*)'BEFORE pousinmc'
    ! lump the pressure nodes to take away the discontinuity...
    DO ELE = 1, TOTELE
       DO CV_ILOC = 1, CV_NLOC
          dg_nod = (ele-1) * cv_nloc + cv_iloc
          cty_nod = x_ndgln( (ele-1) * cv_nloc + cv_iloc)
          MAP_DG2CTY(dg_nod) = cty_nod
       END DO
    END DO

    CALL GET_SPAR_CMC_SMALL(FINDCMC_SMALL, COLCMC_SMALL, MIDCMC_SMALL, &
         MX_NCMC_SMALL, NCMC_SMALL, CV_NONODS, X_NONODS, MAP_DG2CTY, &
         FINDCMC, COLCMC, NCOLCMC)

    ewrite(3,*)'FINDCMC_SMALL:', FINDCMC_SMALL
    ewrite(3,*)'COLCMC_SMALL(1:NCMC_SMALL):', COLCMC_SMALL(1:NCMC_SMALL)
    !stop 3292

    allocate( cmc_small(ncmc_small) )
    allocate( resid_dg(cv_nonods) )
    allocate( resid_cty(x_nonods) )
    allocate( dp_small(x_nonods) )
    allocate( dp_dg(cv_nonods) )
    allocate( nods_sourou(x_nonods) )

    ewrite(3,*)'***forming cmc_small:'
    CMC_SMALL = 0.
    DO dg_nod = 1, CV_NONODS
       cty_nod = MAP_DG2CTY(dg_nod)
       ! add row dg_nod to row cty_nod of cty mesh
       DO COUNT = FINDCMC(DG_NOD), FINDCMC(DG_NOD+1) - 1
          jcolcmc = COLCMC(COUNT)
          jcolcmc_small = MAP_DG2CTY(jcolcmc)
          count2 = 0
          DO COUNT3 = FINDCMC_small(cty_NOD), FINDCMC_small(cty_NOD+1) - 1
             !ewrite(3,*)'dg_nod,cty_nod,jcolcmc_small,colcmc_small(count3):', &
             !     dg_nod,cty_nod,jcolcmc_small,colcmc_small(count3)
             if(colcmc_small(count3)==jcolcmc_small) count2=count3
          end do
          if(count2==0) then
             ewrite(3,*)'could not find coln'
             stop 3282
          end if
          CMC_SMALL(COUNT2) = CMC_SMALL(COUNT2) + CMC(COUNT)              
       END DO
    END DO

    DO GL_ITS = 1, NGL_ITS

       EWRITE(3,*)'GL_ITS=', GL_ITS
       ! SSOR smoother for the multi-grid method...
       ewrite(3,*)'before solving:',p

       if (gl_its>2) then

          if (.true.) then
             call set_solver_options(path, &
                  !ksptype = "gmres", &
                  pctype = "jacobi", & ! use this for P1DGP1DG
                  !pctype = "none", &   ! use this for P1DGP2DG
                  rtol = 1.e-10, &
                  atol = 1.e-15, &
                  max_its = 25)
             ! ignore solver failures...
             call add_option( &
                  trim(path)//"/solver/ignore_all_solver_failures", stat)
             CALL SOLVER( CMC, P, RHS, &
                  FINDCMC, COLCMC, &
                  option_path = path )
          end if

          if (.false.) then
             CALL SIMPLE_SOLVER( CMC, P, RHS,  &
                  NCOLCMC, CV_NONODS, FINDCMC, COLCMC, MIDCMC,  &
                  ERROR, RELAX, RELAX_DIAABS, RELAX_DIA, N_LIN_ITS )
             ewrite(3,*)'after solving:',p
          end if

          if (.false.) then
             CALL SOLVER( CMC, P, RHS, &
                  FINDCMC, COLCMC, &
                  option_path = '/material_phase[0]/scalar_field::Pressure' )
          end if

       end if

       resid_dg = rhs
       do dg_nod = 1, cv_nonods
          DO COUNT = FINDCMC(dg_NOD), FINDCMC(dg_NOD+1) - 1
             resid_dg(dg_nod) = resid_dg(dg_nod) - cmc(count) * P(COLCMC(COUNT))
          END DO
       end do

       ! Map resid_dg to resid_cty as well as the solution:
       resid_cty=0.
       nods_sourou=0.
       do dg_nod = 1, cv_nonods
          cty_nod = MAP_DG2CTY(dg_nod)
          resid_cty(cty_nod) = resid_cty(cty_nod)+resid_dg(dg_nod)
          nods_sourou(cty_nod) = nods_sourou(cty_nod)+1.
       end do
       ! We have added the rows together so no need to normalize residual. 
       resid_cty= resid_cty/nods_sourou

       ! Course grid solver...
       DP_SMALL = 0.
       EWRITE(3,*)'SOLVER'
       CALL SOLVER( CMC_SMALL(1:NCMC_SMALL), DP_SMALL, resid_cty, &
            FINDCMC_SMALL, COLCMC_SMALL(1:NCMC_SMALL), &
            option_path = '/material_phase[0]/scalar_field::Pressure')
       EWRITE(3,*)'OUT OF SOLVER'

       ! Map the corrections DP_SMALL to dg:
       DO dg_nod = 1, cv_nonods
          cty_nod = MAP_DG2CTY(dg_nod)
          DP_DG(DG_NOD) = DP_SMALL(CTY_NOD)
       END DO

       P = P + DP_DG

    END DO

    deallocate( FINDcmc_small, &
         colcmc_small, &
         midcmc_small, &
         MAP_DG2CTY )

    deallocate( cmc_small, &
         resid_dg, &
         resid_cty, &
         dp_small, &
         dp_dg, &
         nods_sourou )

    RETURN
  END SUBROUTINE PRES_DG_MULTIGRID


  SUBROUTINE GET_SPAR_CMC_SMALL(FINDCMC_SMALL,COLCMC_SMALL,MIDCMC_SMALL, &
       MX_NCMC_SMALL,NCMC_SMALL, CV_NONODS,X_NONODS, MAP_DG2CTY, &
       FINDCMC,COLCMC,NCOLCMC)
    ! Form sparcity COLCMC_SMALL, FINDCMC_SMALL 
    ! from FINDCMC,COLCMC...
    ! It lumps the DG pressure matrix to a continuous pressure matrix...
    INTEGER, intent( in ) ::  MX_NCMC_SMALL,NCOLCMC, CV_NONODS,X_NONODS
    INTEGER, intent( inout ) ::  NCMC_SMALL
    INTEGER, DIMENSION( x_NONODS + 1 ), intent( inout ) :: FINDCMC_SMALL
    INTEGER, DIMENSION( MX_NCMC_SMALL ), intent( inout ) :: COLCMC_SMALL
    INTEGER, DIMENSION( x_nonods ), intent( inout ) :: MIDCMC_SMALL

    INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCMC
    INTEGER, DIMENSION( NCOLCMC ), intent( in ) :: COLCMC

    INTEGER, DIMENSION( CV_nonods ), intent( in ) :: MAP_DG2CTY
    ! Local variables
    integer, dimension( : ), allocatable :: MX_NODS_ROW_SMALL, NODS_ROW_SMALL, &
         FINDCMC_SMALL_mx
    INTEGER :: dg_nod,cty_nod,count,jcolcmc,jcolcmc_small,count2,count3

    allocate(MX_NODS_ROW_SMALL(x_nonods))
    allocate(NODS_ROW_SMALL(x_nonods))
    allocate(FINDCMC_SMALL_mx(x_nonods+1))

    MX_NODS_ROW_SMALL=0
    DO dg_nod=1,CV_NONODS
       cty_nod=MAP_DG2CTY(dg_nod)
       MX_NODS_ROW_SMALL(cty_nod)=min( MX_NODS_ROW_SMALL(cty_nod)  &
            +FINDCMC(DG_NOD+1)-FINDCMC(DG_NOD),   x_nonods)
    END DO
    FINDCMC_small_MX(1)=1
    do cty_nod=2,x_nonods+1
       FINDCMC_small_MX(cty_nod)=FINDCMC_small_MX(cty_nod-1)+mx_NODS_ROW_SMALL(cty_nod-1)
    end do
    ewrite(3,*)'MAP_DG2CTY:',MAP_DG2CTY
    ewrite(3,*)'mx_NODS_ROW_SMALL:',mx_NODS_ROW_SMALL
    ewrite(3,*)'FINDCMC_small_MX:',FINDCMC_small_MX

    NODS_ROW_SMALL=0
    COLCMC_SMALL=0
    DO dg_nod=1,CV_NONODS
       cty_nod=MAP_DG2CTY(dg_nod)
       ! add row dg_nod to row cty_nod of cty mesh
       DO COUNT=FINDCMC(DG_NOD),FINDCMC(DG_NOD+1)-1
          jcolcmc=COLCMC(COUNT)
          jcolcmc_small=MAP_DG2CTY(jcolcmc)

          count2=0
          DO COUNT3=FINDCMC_small_mx(CTY_NOD),FINDCMC_SMALL_mx(CTY_NOD)+NODS_ROW_SMALL(cty_nod)-1
             if(colcmc_small(count3)==jcolcmc_small) count2=count3
          end do
          if(count2==0) then ! then put coln in as we have not found it in row 
             NODS_ROW_SMALL(cty_nod)=NODS_ROW_SMALL(cty_nod)+1
             COLCMC_SMALL(FINDCMC_small_mx(CTY_NOD)+NODS_ROW_SMALL(cty_nod)-1)=jcolcmc_small
             if(cty_nod==1) then
                ewrite(3,*)'dg_nod,cty_nod,jcolcmc_small:',dg_nod,cty_nod,jcolcmc_small
             endif
          endif
       END DO
    END DO
    ewrite(3,*)'NODS_ROW_SMALL:',NODS_ROW_SMALL

    FINDCMC_small(1)=1
    do cty_nod=2,x_nonods+1
       FINDCMC_small(cty_nod)=FINDCMC_small(cty_nod-1)+NODS_ROW_SMALL(cty_nod-1)
    end do
    NCMC_SMALL=FINDCMC_small(x_nonods+1)-1

    ! Shrink up the pointer list COLCMC_SMALL: 
    COUNT=0
    do cty_nod=1,x_nonods
       DO COUNT2=FINDCMC_small_MX(cty_nod),FINDCMC_small_MX(cty_nod+1)-1
          IF(COLCMC_SMALL(COUNT2).NE.0) THEN
             COUNT=COUNT+1
             COLCMC_SMALL(COUNT)=COLCMC_SMALL(COUNT2)
          ENDIF
       END DO
    END DO
    ! Put in assending coln order in each row...
    do cty_nod=1,x_nonods
       ewrite(3,*)'cty_nod,FINDCMC_small(cty_nod),FINDCMC_small(cty_nod+1)-1:', &
            cty_nod,FINDCMC_small(cty_nod),FINDCMC_small(cty_nod+1)-1
       call ibubble2(COLCMC_SMALL(FINDCMC_small(cty_nod):FINDCMC_small(cty_nod+1)-1))
       ewrite(3,*)'COLCMC_SMALL(FINDCMC_small(cty_nod):FINDCMC_small(cty_nod+1)-1):', &
            COLCMC_SMALL(FINDCMC_small(cty_nod):FINDCMC_small(cty_nod+1)-1)
    end do
    ! Calculate MIDCMC_SMALL...
    do cty_nod=1,x_nonods
       DO COUNT=FINDCMC_small(cty_nod),FINDCMC_small(cty_nod+1)-1
          IF(COLCMC_SMALL(COUNT)==cty_nod) MIDCMC_SMALL(CTY_NOD)=COUNT
       END DO
    END DO

    RETURN
  END SUBROUTINE GET_SPAR_CMC_SMALL
         


    subroutine ibubble2(ivec)
      ! sort ivec in increasing order
      implicit none
      integer, dimension( : ), intent( inout ) :: ivec
      ! Local variables
      integer :: nvec, i, j, itemp

      nvec = size(ivec)

!        ewrite(3,*)'before ivec:',ivec

      do j = 1, nvec
         do i = 1, nvec - 1
            if ( ivec( i ) > ivec( i + 1 ) ) then
               itemp = ivec( i + 1 )
               ivec( i + 1 ) = ivec( i )
               ivec( i ) = itemp
            end if
         end do
      end do
!        ewrite(3,*)'after ivec:',ivec
      return
    end subroutine ibubble2


    SUBROUTINE SIMPLE_SOLVER( CMC, P, RHS,  &
         NCMC, NONODS, FINCMC, COLCMC, MIDCMC,  &
         ERROR, RELAX, RELAX_DIAABS, RELAX_DIA, N_LIN_ITS )
      !
      ! Solve CMC * P = RHS for RHS.
      ! RELAX: overall relaxation coeff; =1 for no relaxation. 
      ! RELAX_DIAABS: relaxation of the absolute values of the sum of the row of the matrix;
      !               - recommend >=2 for hard problems, =0 for easy
      ! RELAX_DIA: relaxation of diagonal; =1 no relaxation (normally applied). 
      ! N_LIN_ITS = no of linear iterations
      ! ERROR= solver tolerence between 2 consecutive iterations
      implicit none
      REAL, intent( in ) :: ERROR, RELAX, RELAX_DIAABS, RELAX_DIA
      INTEGER, intent( in ) ::  N_LIN_ITS, NCMC, NONODS
      REAL, DIMENSION( NCMC ), intent( in ) ::  CMC
      REAL, DIMENSION( NONODS ), intent( inout ) ::  P
      REAL, DIMENSION( NONODS ), intent( in ) :: RHS
      INTEGER, DIMENSION( NONODS + 1 ), intent( in ) :: FINCMC
      INTEGER, DIMENSION( NCMC ), intent( in ) :: COLCMC
      INTEGER, DIMENSION( NONODS ), intent( in ) :: MIDCMC
      ! Local variables
      INTEGER :: ITS, ILOOP, ISTART, IFINI, ISTEP, NOD, COUNT
      REAL :: R, SABS_DIAG, RTOP, RBOT, POLD, MAX_ERR
      LOGICAL :: jacobi

      ewrite(3,*) 'In Solver'

      jacobi = .false. !.true.

      if(jacobi) then

      Loop_Non_Linear_Iter1: DO ITS = 1, N_LIN_ITS

            Loop_Nods1: DO NOD = 1,NONODS
               R =  CMC( MIDCMC( NOD ))*P(NOD)+ RHS( NOD )
               DO COUNT = FINCMC( NOD ), FINCMC( NOD + 1 ) - 1
                  R = R - CMC( COUNT ) * P( COLCMC( COUNT ))
               END DO
               RTOP = R 
               RBOT = CMC( MIDCMC( NOD ))
               P( NOD ) = RELAX * ( RTOP / RBOT ) + ( 1.0 - RELAX ) * P( NOD )
            END DO Loop_Nods1

         END DO Loop_Non_Linear_Iter1
     else

      Loop_Non_Linear_Iter: DO ITS = 1, N_LIN_ITS

         MAX_ERR = 0.0
         Loop_Internal: DO ILOOP = 1, 2
            IF( ILOOP == 1 ) THEN
               ISTART = 1
               IFINI = NONODS
               ISTEP = 1
            ELSE
               ISTART = NONODS
               IFINI = 1
               ISTEP = -1
            ENDIF

            Loop_Nods: DO NOD = ISTART, IFINI, ISTEP
               R = RELAX_DIA * CMC( MIDCMC( NOD )) * P( NOD ) + RHS( NOD )
               SABS_DIAG = 0.0
               DO COUNT = FINCMC( NOD ), FINCMC( NOD + 1 ) - 1
                  R = R - CMC( COUNT ) * P( COLCMC( COUNT ))
                  SABS_DIAG = SABS_DIAG + ABS( CMC( COUNT ))
               END DO
               RTOP = R + RELAX_DIAABS * SABS_DIAG * P( NOD )
               RBOT = RELAX_DIAABS * SABS_DIAG + RELAX_DIA * CMC( MIDCMC( NOD ))
               POLD = P( NOD )
               P( NOD ) = RELAX * ( RTOP / RBOT ) + ( 1.0 - RELAX ) * P( NOD )
               MAX_ERR = MAX( MAX_ERR, ABS( POLD - P( NOD )))
            END DO Loop_Nods
         END DO Loop_Internal

         IF( MAX_ERR < ERROR ) CYCLE

      END DO Loop_Non_Linear_Iter
        endif

      ewrite(3,*) 'Leaving Solver'

      RETURN
    END SUBROUTINE SIMPLE_SOLVER



end module solvers_module


