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

SUBROUTINE VTKOUT(outName, vtkTitle, indexName, NNodes, NElems, X, Y, Z, &
     FIELDS, FSTRIDE, NFIELDS,                                           &
     NDGLNO, sizeNDGLNO, elementTypes, elementSizes, NODDOM)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN)::outName, vtkTitle, indexName
  INTEGER, INTENT(IN)::NNodes, NElems, FSTRIDE, NFIELDS, sizeNDGLNO,      &
       NDGLNO(sizeNDGLNO), NODDOM(NElems), elementTypes(NElems), elementSizes(NElems)
  REAL, INTENT(IN):: X(NNodes), Y(NNodes), Z(NNodes), FIELDS(NNodes*NFIELDS)
  
  INTEGER I, J, IERR, position, cell(100)
  LOGICAL haveIndex, moreFields
  CHARACTER(LEN=2048) type, dataName
  INTEGER  type_len, dataName_len
  REAL*4 dataArray(9), dataArray2(9)
  REAL, ALLOCATABLE, DIMENSION(:)::RNODDOM
  INTEGER ids(100), ids_cnt
  ids_cnt = 100;

#ifdef DEBUG
  ewrite(1,*)  'Just inside VTKOUT()'
  ewrite(3,*)  'outName  = ', outName
  ewrite(3,*)  'vtkTitle = ', vtkTitle
#endif
  ! Open and initalize the file
  CALL vtkopen(outName, LEN_TRIM(outName), vtkTitle, LEN_TRIM(vtkTitle))
  
#ifdef DEBUG
  ewrite(3,*)  'Writting out mesh geometry.'
#endif
  CALL vtkwritemesh(NNodes, NElems, X, Y, Z, NDGLNO, elementTypes, elementSizes)

  ! Output all the field-data
  CALL FL_FLDopen(indexName, LEN_TRIM(indexName), IERR)
  IF(IERR==0) THEN
     haveIndex = .TRUE.
  ELSE
     haveIndex = .FALSE.     
  END IF

  CALL vtkstartN()
  
  IF( haveIndex ) THEN
     
     DO
        CALL FL_FLDread(type, type_len, dataName, dataName_len, ids, ids_cnt, IERR)
        IF(IERR .NE. 0) EXIT
        DO J=1, ids_cnt
           ids(J) = ids(J)-1
        END DO
        
        ! Scalar support
        IF(type(1:type_len) .EQ. 'SCALAR') THEN

           CALL vtkwriteFSN(FIELDS(ids(1)*FSTRIDE+1), dataName, dataName_len)
           
           ! Vector support
        ELSE IF(type(1:type_len) .EQ. 'VECTOR') THEN

           CALL vtkwriteFVN(                                    &
                FIELDS(ids(1)*FSTRIDE+1),                     &
                FIELDS(ids(2)*FSTRIDE+1),                     &
                FIELDS(ids(3)*FSTRIDE+1), dataName, dataName_len)
           
        ! Tensor support
        ELSE IF(type(1:type_len) .EQ. 'TENSOR') THEN

           CALL vtkwriteFTN(                                  &
                FIELDS(ids(1)*FSTRIDE+1),                     &
                FIELDS(ids(2)*FSTRIDE+1),                     &
                FIELDS(ids(3)*FSTRIDE+1),                     &
                FIELDS(ids(4)*FSTRIDE+1),                     &
                FIELDS(ids(5)*FSTRIDE+1),                     &
                FIELDS(ids(6)*FSTRIDE+1),                     &
                FIELDS(ids(7)*FSTRIDE+1),                     &
                FIELDS(ids(8)*FSTRIDE+1),                     &
                FIELDS(ids(9)*FSTRIDE+1),                     &
                dataName, dataName_len)
           
        ELSE IF(type(1:type_len) .EQ. 'CONDUCTIVITY') THEN
           
           CALL vtkwriteFTN2(                                   &
                FIELDS(ids(1)*FSTRIDE+1),                     &
                FIELDS(ids(2)*FSTRIDE+1),                     &
                FIELDS(ids(3)*FSTRIDE+1),                     &
                FIELDS(ids(4)*FSTRIDE+1),                     &
                FIELDS(ids(5)*FSTRIDE+1),                     &
                FIELDS(ids(6)*FSTRIDE+1),                     &
                dataName, dataName_len)
           
        END IF
     END DO
     
     CALL FL_FLDclose();
     
  ELSE
     
     ! Default action is to output all the fields as scalars
     DO I=0, NFIELDS-1
        CALL vtkwriteFSN(FIELDS(I*FSTRIDE+1), "", 0)
     END DO
     
  END IF

  CALL vtkstartC();  
  ALLOCATE( RNODDOM(NElems) )
  DO I=1, NElems
    RNODDOM(I) = NODDOM(I)
  END DO
  CALL vtkwriteFSC(RNODDOM, "GraphPartitioning", 17)
  DEALLOCATE(RNODDOM)

  ! Finalize
  CALL vtkclose(outName, LEN_TRIM(outName), vtkTitle, LEN_TRIM(vtkTitle))
  
  RETURN
  
END SUBROUTINE VTKOUT
















