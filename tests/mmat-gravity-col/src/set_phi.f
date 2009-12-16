C	Function to set the initial volume fraction

        IMPLICIT NONE

C	X, Y, Z,  coordinates from GEM
C	T; value returned to GEM
C	I; counter to read coordinates from GEM

        REAL X, Y, Z
        REAL T
        INTEGER I
C
     


        DO I=0, 1500000
                READ(*, END=100, FMT=*) X, Y, Z

                T = 0.0
 
                IF(Y.LE.0.5.AND.Y.GE.0.0) THEN
             
                  T = 1.0
        
                ENDIF
 
                WRITE(*,*) T
        END DO
100     CONTINUE
 
        STOP
        END

