
!!!=========================================!!!
!!!           SOLVERS SUBROUTINES           !!!
!!!=========================================!!!



module solvers_module

contains

  SUBROUTINE SOLVER( CMC, P, RHS,  &
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

    write(357,*) 'In Solver'

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

    write(357,*) 'Leaving Solver'

    RETURN
  END SUBROUTINE SOLVER


  SUBROUTINE GETCMC( CMC, CTP, DIAINV, CT, FREDOP, NONODS, &
       NCOLCT, FINDCT, COLCT, &
       NCMC, FINCMC, COLCMC )
    IMPLICIT NONE
    INTEGER, intent( in ) :: FREDOP, NONODS, NCOLCT, NCMC
    INTEGER, DIMENSION( NCMC ), intent( inout ) :: CMC
    INTEGER, DIMENSION( NCOLCT ), intent( in ) :: CTP
    INTEGER, DIMENSION( NONODS ), intent( in ) :: DIAINV
    INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT, CT
    INTEGER, DIMENSION( FREDOP + 1 ), intent( in ) :: FINDCT, FINCMC
    INTEGER, DIMENSION( NCMC ), intent( in ) :: COLCMC

    ! Local variables...
    INTEGER :: ROW, START, FINISH, ST2, FI2, COL, ST1, FI1, I, J, ACR, ICOL
    REAL :: CTC

    write(357,*) 'In GetCMC'

    !write(357,*)'ncmc,FREDOP,nonods=',ncmc,FREDOP,nonods

    CMC( 1 : NCMC ) = 0.0

    Loop_Row: DO ROW = 1, FREDOP
       START = FINCMC( ROW )
       FINISH = FINCMC( ROW + 1 ) - 1
       ST2 = FINDCT( ROW )
       FI2 = FINDCT( ROW + 1 ) - 1

       DO ACR = START, FINISH, 1
          COL = COLCMC( ACR )
          CTC = 0.
          ! Multiplying vector column 'ROW' of C1 by vector column 'COL' of C1.
          ST1 = FINDCT( COL )
          FI1 = FINDCT( COL + 1 ) - 1

          Loop_I: DO I = ST1, FI1
             Loop_J: DO J = ST2, FI2
                IF( COLCT( I ) == COLCT( J )) THEN
                   ICOL = COLCT( I )
                   CTC = CTC + CT( I ) * CTP( J ) * DIAINV( ICOL )
                END IF
             END DO Loop_J
          END DO Loop_I

          CMC( ACR ) = CTC

       END DO

    END DO Loop_Row

    write(357,*) 'Leaving GetCMC'

    RETURN
  END SUBROUTINE GETCMC


  subroutine getnc( totele, nonods, nloc, nc )
    implicit none
    integer, intent( in ) :: totele, nonods, nloc
    integer, intent( inout ) :: nc
    ! Local
    integer :: nod, ele, ele1, ele2, count

    write(357,*) 'In GetNC'

    count = 0
    Loop_cvnod: do nod = 1, nonods
       ele1 = 1 + ( nod - 2 ) / ( nloc - 1 ) 
       ele2 = 1 + ( nod - 1 ) / ( nloc - 1 ) 

       Loop_element: do ele = max( 1, ele1 ), min( totele, ele2 ), 1
          count = count + 1
       end do Loop_element

    end do Loop_cvnod

    nc = count

    write(357,*) 'Leaving GetNC'

    return
  end subroutine getnc

  subroutine getnele( band, nonods, nele )
    implicit none
    integer, intent( in ) :: band, nonods
    integer, intent( inout ) :: nele
    ! band is the semi-bandwidth. 
    ! Local
    integer :: nod, count, iband, col

    write(357,*) 'In GetNELE'

    count = 0
    col = 0
    Loop_cvnod: do nod = 1, nonods

       Loop_Band: do iband = -band, band, 1
          col = nod + iband
          if( ( col >= 1 ) .and. ( col <= nonods ) ) count = count + 1
       end do Loop_Band

    end do Loop_cvnod

    nele = count

    print *,'GetNELE'

    return
  end subroutine getnele


end module solvers_module


