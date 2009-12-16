SUBROUTINE ExtractData(filename, flength, nelm, nnod, node_groups, szenls, nfields, nproperties,  &
     ndimensions, maxnamelen, vol_tracer, vol_tracer_name, int_tracer, int_tracer_name, node_positions)

!--------------------------------------------------------------------------------------------------

! ExtractConnectivity is a subroutine that yoinks data out of a vtu file.
       
  IMPLICIT NONE 

! Variable declarations for main subroutine.
  INTEGER :: err, i
  INTEGER, INTENT(IN) :: flength, ndimensions, nelm, nnod, nproperties, nfields, szenls, maxnamelen         

  INTEGER :: enlbas(nelm+1), enlist(szenls)

  REAL ::  x(nnod), y(nnod), z(nnod)
  REAL ::  fields(nnod, nfields), properties(nelm, nproperties)
  REAL, INTENT(INOUT) :: int_tracer(nnod), node_groups(nelm,4), node_positions(nnod,3), vol_tracer(nnod)

  CHARACTER*(*), INTENT(IN) :: filename, int_tracer_name, vol_tracer_name
  CHARACTER(len = maxnamelen) :: names(nfields+nproperties)

  LOGICAL :: got_vol_tracer, got_int_tracer

!--------------------------------------------------------------------------------------------------

  INTERFACE
     INTEGER FUNCTION freadvtkfile( filename, namelen, nnod, nelm, szenls,     &
          nfields, nproperties, ndimensions, x, y, z, fields, properties,     &
          enlbas, enlist, names, maxnamelen)
       CHARACTER*(*) :: filename
       INTEGER :: namelen, nnod, nelm, szenls, nfields, nproperties, ndimensions
       INTEGER :: enlbas(nelm+1), enlist(szenls)
       REAL :: x(nnod), y(nnod), z(nnod), fields(nnod,nfields), properties(nelm,nproperties)  
       CHARACTER(len = maxnamelen) :: names(nfields+nproperties)
     END FUNCTION freadvtkfile
  END INTERFACE

!--------------------------------------------------------------------------------------------------
! Read in the vtk file.

  err = freadvtkfile(filename, flength, nnod, nelm, szenls, nfields,     &
       nproperties, ndimensions, x, y, z, fields, properties, enlbas,     &
       enlist, names, maxnamelen)

  IF( err /= 0 ) THEN
     PRINT *, 'Something went wrong with VTK read!'
     STOP
  END IF

!--------------------------------------------------------------------------------------------------

! Extract the global node numbers for each element into the node_groups array.

  DO i = 1,nelm
     node_groups(i,:) = enlist( enlbas(i)+1:enlbas(i+1) )
  ENDDO

!--------------------------------------------------------------------------------------------------
! Extract the node positions.

  node_positions(:,1) = x
  node_positions(:,2) = y
  node_positions(:,3) = z

!--------------------------------------------------------------------------------------------------
! Extract the tracer fields for the volume defining field and the tracer to be integrated.

! Initialise flags to indicate whether the fields were successfully found.
  got_vol_tracer = .FALSE.
  got_int_tracer = .FALSE.

  DO i = 1,nfields

     IF( TRIM(vol_tracer_name) == TRIM(names(i)) )THEN
        vol_tracer(:) = fields(:,i)
        got_vol_tracer = .TRUE.
        PRINT*,'Successfully loaded vol_tracer field.'
     END IF

     IF( TRIM(int_tracer_name) == TRIM(names(i)) )THEN
        int_tracer(:) = fields(:,i)
        got_int_tracer = .TRUE.
        PRINT*,'Successfully loaded int_tracer field.'
    END IF

  END DO

! Warn the user if the fields weren't there and STOP.
  IF( got_vol_tracer .NEQV. .TRUE. )THEN
     PRINT*,'ERROR: Cannot find field named ',TRIM(vol_tracer_name),'.'
     PRINT*,'The fields available to you are:', names(:)
     PRINT*,'Could not extract vol_tracer array.'
     PRINT*,'STOPPING.'
     STOP
  ELSEIF( got_int_tracer .NEQV. .TRUE. )THEN
     PRINT*,'ERROR: Cannot find field named ',TRIM(int_tracer_name),'.'
     PRINT*,'The fields available to you are:', names(:)
     PRINT*,'Could not extract int_tracer array.'
     PRINT*,'STOPPING.'
     STOP
  ELSE
  END IF

!--------------------------------------------------------------------------------------------------

END SUBROUTINE ExtractData

!--------------------------------------------------------------------------------------------------
