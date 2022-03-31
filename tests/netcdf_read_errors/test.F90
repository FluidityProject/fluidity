PROGRAM test_getenv
  CHARACTER(LEN=100) :: hostname
  INTEGER :: status_value = 0

  CALL GET_ENVIRONMENT_VARIABLE("HOSTNAME", hostname, STATUS=status_value)

  IF (status_value == 2) THEN
    WRITE(*,*) 'WARNING: Processor does not support environment variables - hostname is unknown.'
    hostname = 'Unknown'
  ELSE IF (status_value == -1) THEN
    WRITE(*,*) 'WARNING: Hostname is too long for character variable - hostname is truncated.'
  ELSE IF (status_value == 1) THEN
    WRITE(*,*) 'WARNING: $HOSTNAME environment variable does not exist - hostname is unknown.'
    hostname = 'Unknown'
  END IF
END PROGRAM
