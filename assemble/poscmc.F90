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
module poscmc_module
  use FLDebug
  implicit none
  private
  public::poscmc
contains 
! This subroutine forms the matrix operating on the pressure vector
! It is found from C1T ML C1 + C2T ML C2.
! In the first part of COLCMC  contains: given a node the 
! pressure nodes surrounding it.
SUBROUTINE POSCMC(FINDCT,COLCT,NCT, &
     FINCMC,COLCMC,NIMEM, &
     NCMC, &
     FREDOP,NONODS,&
     MIDCMC, &
     NOINOD ,PRESYM)

  INTEGER, INTENT(IN)::NONODS,FREDOP,NIMEM,NCT
  INTEGER, INTENT(IN)::FINDCT(FREDOP+1),COLCT(NCT)
  INTEGER, INTENT(OUT)::FINCMC(FREDOP+1),COLCMC(NIMEM),NCMC
  ! Call this sub after we get the CT pointers and before 
  ! we get the M pointers. 
  INTEGER, INTENT(OUT)::MIDCMC(FREDOP)
  INTEGER, INTENT(OUT)::NOINOD(NONODS)
  LOGICAL, INTENT(IN)::PRESYM

  INTEGER COUNT,COUNT2
  INTEGER GLOBI,GLOBJ,I,IROW,INOD,JROW,JNOD,PTR
  
  ! Define a linked list
  TYPE node
     INTEGER :: ID                ! id number of node
     TYPE (node), POINTER :: next ! next node
  END TYPE node
  
  TYPE row
     TYPE (node), POINTER :: row
  END TYPE row
  
  TYPE(row), DIMENSION(:), ALLOCATABLE::Matrix, Matrix2
  TYPE(node), POINTER::List, Current, Current2, Next

  ! Initalise the linked lists
  ALLOCATE( Matrix(NONODS) )
  DO I=1, NONODS
     ALLOCATE( List )
     List%ID = -1
     NULLIFY( List%next )
     
     Matrix(I)%row => List
     NULLIFY(List)
  END DO

  NOINOD(1:NONODS)=0

  ! Given a vel node find the pressure nodes surrounding it. 
  DO IROW=1,FREDOP
     GLOBJ=IROW
     DO COUNT=FINDCT(IROW),FINDCT(IROW+1)-1
        GLOBI=COLCT(COUNT)
        List => Matrix(GLOBI)%row

        ! Check if the list is initalised
        IF(List%ID.EQ.-1) THEN
           List%ID = GLOBJ
           CYCLE
        END IF
        
        IF(GLOBJ.LT.List%ID) THEN
           ! Insert at start of list
           ALLOCATE(Current)
           Current%ID = GLOBJ
           Current%next => List
           
           Matrix(GLOBJ)%row => Current
           List => Matrix(GLOBJ)%row
        ELSE
           Current => List
           DO WHILE ( ASSOCIATED(Current) )
              IF(GLOBJ.EQ.Current%ID) THEN
                 ! Already have this node
                 exit 
              ELSE IF(.NOT.ASSOCIATED(Current%next)) THEN
                 ! End of list - insert this node
                 ALLOCATE(Current%next)
                 NULLIFY(Current%next%next)
                 Current%next%ID = GLOBJ
                 exit
              ELSE IF(GLOBJ.LT.Current%next%ID) THEN
                 ! Insert new node here
                 ALLOCATE(Next)
                 Next%ID = GLOBJ
                 Next%next => Current%next
                 Current%Next => Next
                 exit
              END IF
              Current => Current%next
           END DO
        END IF
        NOINOD(GLOBI)=NOINOD(GLOBI)+1
     END DO
  END DO

  ! Initalise the linked lists
  ALLOCATE( Matrix2(FREDOP) )
  DO I=1, FREDOP
     ALLOCATE( List )
     List%ID = -1
     NULLIFY( List%next )
     
     Matrix2(I)%row => List
     NULLIFY(List)
  END DO
  
  DO IROW=1,FREDOP 
     ! Find the PRESSURE NODES surrounding pressure node IROW
     DO COUNT=FINDCT(IROW),FINDCT(IROW+1)-1
        INOD=COLCT(COUNT)

        ! Find the pressure nodes surrounding node INOD 
        ! these will be connected to pressure node IROW. 
        Current => Matrix(INOD)%row
        DO WHILE ( ASSOCIATED(Current) )
           JROW = Current%ID
           IF((.NOT.PRESYM).OR.(JROW.GE.IROW)) THEN
              List => Matrix2(IROW)%row

              ! Check if the list is initialised
              IF(List%ID.EQ.-1) THEN
                 List%ID = JROW
                 CYCLE
              END IF
              
              IF(JROW.LT.List%ID) THEN
                 ! Insert at start of list
                 ALLOCATE(Current2)
                 Current2%ID = JROW
                 Current2%next => List
                 
                 Matrix2(IROW)%row => Current2
                 List => Matrix2(IROW)%row
              ELSE
                 Current2 => List
                 DO WHILE ( ASSOCIATED(Current2) )
                    IF(JROW.EQ.Current2%ID) THEN
                       ! Already have this node
                       exit 
                    ELSE IF(.NOT.ASSOCIATED(Current2%next)) THEN
                       ! End of list - insert this node
                       ALLOCATE(Current2%next)
                       NULLIFY(Current2%next%next)
                       Current2%next%ID = JROW
                       
                       exit
                    ELSE IF(JROW.LT.Current2%next%ID) THEN
                       ! Insert new node here
                       ALLOCATE(Next)
                       Next%ID = JROW
                       Next%next => Current2%next
                       Current2%Next => Next
                       exit
                    END IF
                    Current2 => Current2%next
                 END DO
              END IF
           END IF
           Current => Current%next
        END DO
     END DO
  END DO

  ! Delete Matrix
  DO IROW=1,NONODS
     Current => Matrix(IROW)%row
     DO WHILE ( ASSOCIATED(Current) )
        Next => Current%next
        DEALLOCATE(Current)
        Current => Next
     END DO
  END DO
  DEALLOCATE(Matrix);

  ! From matrix write COLCMC, FINCMC and MIDCMC
  ! linked list as we go
  PTR = 1
  DO IROW=1,FREDOP
     MIDCMC(IROW) = -1
     FINCMC(IROW) = PTR
  !OVERWRITE HERE:  !!!!
     Current => Matrix2(IROW)%row
     
     DO WHILE ( ASSOCIATED(Current) )
        IF(PTR.GT.NIMEM) THEN
           ewrite(-1,*) 'POSCMC: NOT ENOUGH INTEGER SPACE'
           ewrite(-1,*) '  NIMEM, PTR: ',NIMEM,PTR
           ewrite(-1,*) '  FREDOP, IROW: ',FREDOP,IROW
           !FLAbort("Integer memory too small")
           STOP
        END IF
        
        COLCMC(PTR)    = Current%ID
        
        IF(Current%ID.EQ.IROW) THEN
           MIDCMC(IROW) = PTR
        END IF
        
        Next => Current%next
        DEALLOCATE(Current)
        Current => Next
        
        PTR = PTR + 1
     END DO
  END DO
  
  NCMC = PTR-1
  FINCMC(FREDOP+1)=NCMC+1
 
  ! Can use this to verify result against old solution
  !  DO IROW=1,FREDOP
  !     ewrite(3,*) "#################################################"
  !     ewrite(3,*)  IROW, COLCMC(MIDCMC(IROW)), (COLCMC(I), I=FINCMC(IROW),FINCMC(IROW+1)-1)
  !     ewrite(3,*) "#################################################"
  !  END DO
  !  STOP
  
END SUBROUTINE POSCMC

end module poscmc_module
