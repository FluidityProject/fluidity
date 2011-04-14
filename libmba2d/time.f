C ==========================================================
      Real Function ANItime(tmdata)
C ==========================================================
C The wrapper for etime(tmdata). ETIME returns a real number 
C which is the total CPU time used by the program. ETIME 
C uses TMDATA to return the user time (tmdata(1)) and the 
C system time (tmdata(2)) separately. 
C
C *** Remarks:
C        1. the user may set tmdata(1) = tmdata(2) = 0.0 and
C           ANItime = 0.0 without breaking the code.
C ==========================================================
      Real etime, tmdata(2)

      ANItime = 0.0
      tmdata = 0.0

      Return
      End

