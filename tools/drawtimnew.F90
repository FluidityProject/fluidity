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
subroutine fprint_backtrace
end subroutine fprint_backtrace

program drawtimenew
  use FLDebug
  integer MXNOS
  PARAMETER(MXNOS=3000000)
  REAL X(MXNOS),Y(MXNOS),Z(MXNOS)
  REAL XYZ(MXNOS)
  REAL MAXX,MAXY,MINY,MINX
  REAL RSUB(70000)
  REAL DT
  INTEGER NONODS,FREDOP,TOTELE,ELE
  INTEGER I,SOMTIM,J,NUMBER,IREP,NUMITS,K,JNUM,IDUM1, IDUM2 
  integer ii 
  INTEGER ILOOP    
  INTEGER MAKINT
  CHARACTER LA*48
  CHARACTER*240 CSTRI
  CHARACTER*240 FILNAM
  CHARACTER*20 CNUMB
  LOGICAL YINT,ZINT
  ! THIS PROGRAM PLOTS VELS FOR LINEAR MODEL
  ! Calculate the number of columns ncolm...
  !
  !
  !        ewrite(3,*) 'enter FILNAM'
  READ '(A240)', FILNAM
  !        READ(*,*) iloop
  !        K=INDEX(FILNAM,' ') - 1 
  !        ewrite(3,*) 'K=',K
  !        ewrite(3,*) 'FILNAM=',FILNAM(1:K)
  !
  !
  !        OPEN(2,FILE='rubish',STATUS= 'UNKNOWN')
  OPEN(2,FILE=FILNAM,STATUS= 'UNKNOWN')
  !
  do  ILOOP=1,1000000! Was loop 5000
     !
     read(2,55) cstri
     !         ewrite(3,*) 'CSTRI=',CSTRI
     !         ewrite(3,*) CSTRI
     !         ewrite(3,'A100') CSTRI
     ewrite(3,55) CSTRI
     if(cstri(1:2).eq.'@@') then
        K=INDEX(CSTRI,' ')-1
        CNUMB=CSTRI(3:K)
        K=INDEX(CNUMB,' ')-1
        JNUM=0
        !           ewrite(3,*) 'k=',k
        i=0
        do Ii=K,1,-1
           i=i+1
           J=MAKINT(CNUMB(Ii:Ii))
           !             ewrite(3,*) 'i,j:',i,j
           JNUM=JNUM+J*( 10**(I-1)  )
        END DO
        !           ewrite(3,*) 'NUMBER OF COLUMNS=',JNUM
        NCOLM=JNUM
        GOTO 5500
     ENDIF
     !
     IF(CSTRi(1:1).NE.'@') THEN
        ewrite(3,*) '*****You are using an old file****'
        ewrite(3,*) 'to draw a graph type: drawtimnew'
        ewrite(3,*) 'the first input is the filenam - not prompted for'
        !           ewrite(3,*) 'enter 2 dummy variables IDUM1, IDUM2'
        !           READ(*,*) IDUM1,IDUM2
        !         ewrite(3,*) '--IDUM1,IDUM2=',IDUM1,IDUM2
        ewrite(3,*) 'ENTER THE NUMBER OF COLUMNS'
        READ(*,*) NCOLM
        !         ewrite(3,*) '--ncolm=',ncolm
        ewrite(3,*) 'enter X-COLUMN, Y-COLUMN'
        CLOSE(2) 
        OPEN(2,FILE=FILNAM,STATUS= 'UNKNOWN')
        GOTO 5500
     ENDIF

  end do ! Was loop 5000

  !
5500 CONTINUE
  ! 

  !         ewrite(3,*) 'enter no of colmns'
  !         read(*,*) ncolm
  !
  !         ewrite(3,*) 'if any columns are neg integrate over time'
  !         ewrite(3,*) 'enter X-COLUMN, Y-COLUMN'
  read(*,*) NXCOL,NYCOL
  !
  if(NXCOL.gt.NCOLM) THEN
     !           ewrite(3,*) 'X-COLUMN TOO LARGE'
     STOP
  ENDIF
  if(NYCOL.gt.NCOLM) THEN
     !           ewrite(3,*) 'Y-COLUMN TOO LARGE'
     STOP
  ENDIF
  OPEN(7,FILE='forlog.data',STATUS= 'UNKNOWN')
  NZCOL=0
  YINT=.FALSE.
  ZINT=.FALSE.
  rmean1=0.
  rmean2=0.
  if(NYCOL.lt.0) YINT=.TRUE.
  if(NZCOL.lt.0) ZINT=.TRUE.
  NYCOL=ABS(NYCOL)
  NZCOL=ABS(NZCOL)
  !       ewrite(3,*) ' NXCOL,NYCOL,NZCOL:', NXCOL,NYCOL,NZCOL
  !
  ICOUNT=0
1000 continue
  !
  !
  do iline=1,int( (ncolm-1)/12 ) +1
     read(2,222,end=3000) &
          &  (Rsub(IT),IT=12*(iline-1)+1,min(12*iline,Ncolm))
  end do

  !
  !
  !        read(2,*,end=3000) (Rsub(IT),IT=1,Ncolm)
  ICOUNT=ICOUNT+1
  X(ICOUNT)=RSUB(NXCOL)
  Y(ICOUNT)=RSUB(NYCOL)
  IF(NZCOL.GT.0) Z(ICOUNT)=RSUB(NZCOL)
  !         write(7,*) x(ICOUNT),'  ',y(ICOUNT)
  !         write(7,222) x(ICOUNT),y(ICOUNT)
  write(7,444) x(ICOUNT),y(ICOUNT)
  i=icount
  rmean1=rmean1+Y(ICOUNT)
  if(NZCOL.GT.0) rmean2=rmean2+z(ICOUNT)
  !            ewrite(3,*) 'i,x(i),y(i):',i,x(i),y(i)
  IF(ICOUNT.GT.MXNOS-200) THEN
     ewrite(3,*) 'LIMIT HAS BEEN REACHED'
     GOTO 3000
  ENDIF
  goto 1000
3000 continue
  CLOSE(2)
  CLOSE(7)
  rmean1=rmean1/real(ICOUNT)
  if(NZCOL.GT.0) rmean2=rmean2/real(ICOUNT)
  !         ewrite(3,*) 'rmean1,rmean2:',rmean1,rmean2
  !
  IF(YINT) THEN
     OPEN(8,FILE='forlog1',STATUS= 'UNKNOWN')
     ! Work out the integral with repect to time of Y
     ewrite(3,*) 'in hereY'
     XYZ(1)=0.
     do I=2,ICOUNT
        DT=X(I)-X(I-1)
        RADD=0.5*(Y(I-1)+Y(I))*DT*3.204e-11
        if(XYZ(i-1).gt.1.e+20) radd=0.
        XYZ(I)=XYZ(I-1)+RADD
        !                ewrite(3,*) 'i,x(i),Y(i):',i,x(i),Y(i)
     END DO
     !
     do I=1,ICOUNT
        Y(I)=XYZ(I)
        if(abs(x(i)).lt.1.e-30) x(i)=0.
        if(abs(y(i)).lt.1.e-30) y(i)=0.
        !                write(8,*) x(I),y(I)
        !                write(8,222) x(I),y(I)
        write(8,444) x(I),y(I)
     END DO
     CLOSE(8)
  ENDIF
  !
  IF(ZINT) THEN
     ! Work out the integral with repect to time of Y
     ewrite(3,*) 'in hereZ'
     ewrite(3,*) 'in hereY'
     XYZ(1)=0.
     do I=2,ICOUNT
        DT=X(I)-X(I-1)
        RADD=0.5*(Z(I-1)+Z(I))*DT*3.204e-11
        if(XYZ(i-1).gt.1.e+20) radd=0.
        XYZ(I)=XYZ(I-1)+RADD
        !                ewrite(3,*) 'i,x(i),Z(i):',i,x(i),Z(i)
     END DO
     !
     do I=1,ICOUNT
        Z(I)=XYZ(I)
     END DO
  ENDIF
  !
  ! work out X-coord.
  MAXX=-1.E+20
  MAXY=-1.E+20
  MINY=1.E+20
  MINX=1.E+20
  do  I=1,ICOUNT! Was loop 20
     MAXX=MAX(X(I),MAXX)
     MAXY=MAX(Y(I),MAXY)
     MINY=MIN(MINY,Y(I))
     MINX=MIN(MINX,X(I))
     IF(NZCOL.GT.0) THEN
        MAXY=MAX(Z(I),MAXY)
        MINY=MIN(MINY,Z(I))
     ENDIF
     !            ewrite(3,*) 'i,x(i),y(i):',i,x(i),y(i)
  end do ! Was loop 20
  ewrite(3,*) 'MAXX,MAXY,MINY,MINX:',MAXX,MAXY,MINY,MINX
  !C           ewrite(3,*) 'drawing Y line'
  !
222 FORMAT(21E13.6)
  !222    FORMAT(21F13.6)
  !222    FORMAT(23E14.6)
444 FORMAT(E13.6,2x,E13.6)
  !222    FORMAT(8E13.6)
  !222    FORMAT(9(1PE20.7))
55 FORMAT(A100)

end program drawtimenew

! ******THIS IS THE END OF MAIN PROG ******
!
!
INTEGER FUNCTION MAKINT(CSTR)
  CHARACTER*1 CSTR
  IF(CSTR.EQ.'0') MAKINT=0
  IF(CSTR.EQ.'1') MAKINT=1
  IF(CSTR.EQ.'2') MAKINT=2
  IF(CSTR.EQ.'3') MAKINT=3
  IF(CSTR.EQ.'4') MAKINT=4
  IF(CSTR.EQ.'5') MAKINT=5
  IF(CSTR.EQ.'6') MAKINT=6
  IF(CSTR.EQ.'7') MAKINT=7
  IF(CSTR.EQ.'8') MAKINT=8
  IF(CSTR.EQ.'9') MAKINT=9
  RETURN
END FUNCTION MAKINT
!
!
!
!
!
