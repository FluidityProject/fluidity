C     ******************************************************************
C     ******************************************************************
      subroutine fparam(epsfeas,epsopt,iprint,ncomp)

      implicit none

C     SCALAR ARGUMENTS
      integer ncomp,iprint
      double precision epsfeas,epsopt

#include "dim.par"
#include "algparam.com"
#include "outtyp.com"
#include "fixvar.com"
#include "slacks.com"
#include "scaling.com"

C     PARAMETERS
      integer nwords
      parameter ( nwords = 14 )

C     LOCAL SCALARS
      logical lss,scl,connected
      integer i,ifirst,ikey,ilast,inum,j,file20_unit
      double precision dnum

C     LOCAL ARRAYS
      character * 80 line
      character * 10 keyword
      character *  4 lsssub
      character *  4 sclsub

C     DATA BLOCKS
      character * 1  addinfo(nwords),lower(26),upper(26)
      character * 10 dictionary(nwords)
      character * 38 description(nwords)

      data lower /'a','b','c','d','e','f','g','h','i','j','k','l','m',
     +            'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      data upper /'A','B','C','D','E','F','G','H','I','J','K','L','M',
     +            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/

      data dictionary( 1) /'INCREMENTA'/
      data dictionary( 2) /'HESSIAN-AP'/
      data dictionary( 3) /'TRUE-HESSI'/
      data dictionary( 4) /'PERFORM-AC'/
      data dictionary( 5) /'DIRECT-SOL'/
      data dictionary( 6) /'SCALE-LINE'/
      data dictionary( 7) /'DO-NOT-REM'/
      data dictionary( 8) /'ADD-SLACKS'/
      data dictionary( 9) /'DO-NOT-SCA'/
      data dictionary(10) /'IGNORE-OBJ'/
      data dictionary(11) /'FEASIBILIT'/
      data dictionary(12) /'OPTIMALITY'/
      data dictionary(13) /'OUTPUT-DET'/
      data dictionary(14) /'NCOMP-ARRA'/

      data description( 1) /'INCREMENTAL-QUOTIENTS-IN-CG           '/
      data description( 2) /'HESSIAN-APPROXIMATION-IN-CG           '/
      data description( 3) /'TRUE-HESSIAN-PRODUCT-IN-CG            '/
      data description( 4) /'PERFORM-ACCELERATION-STEP             '/
      data description( 5) /'DIRECT-SOLVER                         '/
      data description( 6) /'SCALE-LINEAR-SYSTEMS                  '/
      data description( 7) /'DO-NOT-REMOVE-FIXED-VARIABLES         '/
      data description( 8) /'ADD-SLACKS                            '/
      data description( 9) /'DO-NOT-SCALE-OBJECTIVE-AND-CONSTRAINTS'/
      data description(10) /'IGNORE-OBJECTIVE-FUNCTION             '/
      data description(11) /'FEASIBILITY-TOLERANCE                 '/
      data description(12) /'OPTIMALITY-TOLERANCE                  '/
      data description(13) /'OUTPUT-DETAIL                         '/
      data description(14) /'NCOMP-ARRAY                           '/

      data addinfo( 1) /' '/
      data addinfo( 2) /' '/
      data addinfo( 3) /' '/
      data addinfo( 4) /' '/
      data addinfo( 5) /' '/
      data addinfo( 6) /' '/
      data addinfo( 7) /' '/
      data addinfo( 8) /' '/
      data addinfo( 9) /' '/
      data addinfo(10) /' '/
      data addinfo(11) /'D'/
      data addinfo(12) /'D'/
      data addinfo(13) /'I'/
      data addinfo(14) /'I'/

C     EXTERNAL FUNCTIONS
      external lss,scl

C     BANNER
      if ( iprintctl(1) ) then
          write(*, 8000)
          write(file10_unit,8000)
      end if

C     OPENING THE SPECIFICATION FILE
      file20_unit_loop: do file20_unit=10, 99
  
        inquire(unit=file20_unit, opened=connected)
  
        if (.not.connected) exit file20_unit_loop
  
      end do file20_unit_loop
      
      open(file20_unit,err=300,file='algencan.dat',status='old')

      if ( iprintctl(2) ) then
          write(*, 9005)
          write(file10_unit,9005)
      end if

C     MAIN LOOP

 100  continue

C     READING LINES
      read(file20_unit,fmt=1000,err=400,end=200) line

C     PROCESS LINES

C     Find first character
      i = 1
 110  if ( i .le. 80 .and. line(i:i) .eq. ' ' ) then
          i = i + 1
          go to 110
      end if

C     Skip blank lines
      if ( i .gt. 80 ) then
          go to 100
      end if

      ifirst = i

C     Skip comments
      if ( line(ifirst:ifirst) .eq. '*' .or.
     +     line(ifirst:ifirst) .eq. '#' ) then
          go to 100
      end if

C     Find the end of the keyword
      i = ifirst + 1
 120  if ( i .le. 80 .and. line(i:i) .ne. ' ' ) then
          i = i + 1
          go to 120
      end if

      ilast = i - 1

C     Obtain the first 10 characters and convert to upper-case
      keyword = '          '
      do i = 1,min( 10, ilast - ifirst + 1 )
          keyword(i:i) = line(ifirst+i-1:ifirst+i-1)
          do j = 1,26
              if ( keyword(i:i) .eq. lower(j) ) then
                  keyword(i:i) = upper(j)
              end if
          end do
      end do

C     Look up the keyword in the dictionary
      i = 1
 130  if ( i .le. nwords .and. keyword .ne. dictionary(i) ) then
          i = i + 1
          go to 130
      end if

C     Ignore unknown keywords
      if ( i .gt. nwords ) then
          if ( iprintctl(2) ) then
              write(*, 9020) line(ifirst:ilast)
              write(file10_unit,9020) line(ifirst:ilast)
          end if
          go to 100
      end if

      ikey = i

C     Read additional information if needed
      if ( addinfo(ikey) .ne. ' ' ) then

C         Skip blanks
          i = ilast + 1
 140      if ( i .le. 80 .and. line(i:i) .eq. ' ' ) then
              i = i + 1
              go to 140
          end if

C         Ignore keywords without the required information
          if ( i .gt. 80 ) then
              if ( iprintctl(2) ) then
                  write(*, 9030) description(ikey)
                  write(file10_unit,9030) description(ikey)
              end if
              go to 100
          end if

C         Read additional information
          if ( addinfo(ikey) .eq. 'I' ) then
              read(unit=line(i:80),fmt=2000) inum

          else if ( addinfo(ikey) .eq. 'D' ) then
              read(unit=line(i:80),fmt=3000) dnum
          end if

      end if

C     Process keyword
      if ( iprintctl(2) ) then
          if ( addinfo(ikey) .eq. ' ' ) then
              write(*, 9040) description(ikey)
              write(file10_unit,9040) description(ikey)
          else if ( addinfo(ikey) .eq. 'I' ) then
              write(*, 9041) description(ikey),inum
              write(file10_unit,9041) description(ikey),inum
          else if ( addinfo(ikey) .eq. 'D' ) then
              write(*, 9042) description(ikey),dnum
              write(file10_unit,9042) description(ikey),dnum
          end if
      end if

C     Set the corresponding algencan argument
      if ( ikey .eq.  1 ) then
          hptype = 'INCQUO'

      else if ( ikey .eq.  2 ) then
          hptype = 'HAPPRO'

      else if ( ikey .eq.  3 ) then
          if ( hlpcoded .or. truehl ) then
              hptype = 'TRUEHL'
          else
              if ( iprintctl(2) ) then
                  write(* ,9100) description(ikey)
                  write(file10_unit,9100) description(ikey)
              end if
          end if

      else if ( ikey .eq.  4 ) then
          if ( .not. truehl ) then
              if ( iprintctl(2) ) then
                  write(* ,9110) description(ikey)
                  write(file10_unit,9110) description(ikey)
              end if
          else if ( .not. lss(lsssub) ) then
              if ( iprintctl(2) ) then
                  write(* ,9120) description(ikey)
                  write(file10_unit,9120) description(ikey)
              end if
          else
              skipacc = .false.
              if ( iprintctl(2) ) then
                  write(* ,9060) lsssub
                  write(file10_unit,9060) lsssub
              end if
          end if

      else if ( ikey .eq.  5 ) then
          if ( .not. truehl ) then
              if ( iprintctl(2) ) then
                  write(* ,9110) description(ikey)
                  write(file10_unit,9110) description(ikey)
              end if
          else if ( .not. lss(lsssub) ) then
              if ( iprintctl(2) ) then
                  write(* ,9120) description(ikey)
                  write(file10_unit,9120) description(ikey)
              end if
          else
              avoidds = .false.
              if ( iprintctl(2) ) then
                  write(* ,9050) lsssub
                  write(file10_unit,9050) lsssub
              end if
          end if

      else if ( ikey .eq.  6 ) then
          if ( lss(lsssub) ) then
              if ( lsssub .eq. 'MA57' ) then
                  sclsys = .true.
                  if ( iprintctl(2) ) then
                      write(* ,9070) lsssub
                      write(file10_unit,9070) lsssub
                  end if
              else if ( scl(sclsub) ) then
                  sclsys = .true.
                  if ( iprintctl(2) ) then
                      write(* ,9080) sclsub
                      write(file10_unit,9080) sclsub
                  end if
              else
                  if ( iprintctl(2) ) then
                      write(* ,9130) description(ikey)
                      write(file10_unit,9130) description(ikey)
                  end if
              end if
          else
              if ( iprintctl(2) ) then
                  write(* ,9130) description(ikey)
                  write(file10_unit,9130) description(ikey)
              end if
          end if

      else if ( ikey .eq.  7 ) then
          rmfixv = .false.

      else if ( ikey .eq.  8 ) then
          slacks = .true.

      else if ( ikey .eq.  9 ) then
          scale = .false.

      else if ( ikey .eq. 10 ) then
          ignoref = .true.

      else if ( ikey .eq. 11 ) then
          epsfeas = dnum

      else if ( ikey .eq. 12 ) then
          epsopt = dnum

      else if ( ikey .eq. 13 ) then
          iprint = inum

      else if ( ikey .eq. 14 ) then
          ncomp = inum
      end if

C     IIERATE
      go to 100

C     END OF LOOP

C     TERMINATIONS

C     CLOSING SPECIFICATION FILE
 200  continue
      close(file20_unit)
      go to 500

C     NO SPECIFICATION FILE
 300  continue
      if ( iprintctl(2) ) then
          write(*, 9000)
          write(file10_unit,9000)
      end if
      go to 500

C     ERROR READING THE SPECIFICATION FILE
 400  continue
      if ( iprintctl(2) ) then
          write(*, 9010)
          write(file10_unit,9010)
      end if
      go to 500

C     PRINTING PARAMETERS VALUES
 500  continue
      if ( iprintctl(2) ) then
          write(* ,4000) hptype,avoidds,skipacc,sclsys,rmfixv,slacks,
     +                   scale,epsfeas,epsopt,iprint,ncomp
          write(file10_unit,4000) hptype,avoidds,skipacc,sclsys,rmfixv,slacks,
     +                   scale,epsfeas,epsopt,iprint,ncomp
      end if

C     NON-EXECUTABLE STATEMENTS

 1000 format(A80)
 2000 format(BN,I20)
 3000 format(BN,F24.0)
 4000 format(/,1X,'ALGENCAN PARAMETERS:',
     +       /,1X,'hptype   = ',     14X,A6,
     +       /,1X,'avoidds  = ',     19X,L1,
     +       /,1X,'skipacc  = ',     19X,L1,
     +       /,1X,'sclsys   = ',     19X,L1,
     +       /,1X,'rmfixv   = ',     19X,L1,
     +       /,1X,'slacks   = ',     19X,L1,
     +       /,1X,'scale    = ',     19X,L1,
     +       /,1X,'epsfeas  = ',8X,1P,D12.4,
     +       /,1X,'epsopt   = ',8X,1P,D12.4,
     +       /,1X,'iprint   = ',        I20,
     +       /,1X,'ncomp    = ',        I20)
 8000 format(/,1X,78('='),
     +       /,1X,'This is ALGENCAN 2.2.1.',
     +       /,1X,'ALGENCAN, an augmented Lagrangian method for ',
     +            'nonlinear programming, is part of',/,1X,'the TANGO ',
     +            'Project: Trustable Algorithms for Nonlinear ',
     +            'General Optimization.',/,1X,'See ',
     +            'http://www.ime.usp.br/~egbirgin/tango/ for details.',
     +       /,1X,78('='))
 9000 format(/,1X,'The optional specification file algencan.dat was ',
     +            'not found in the current',/,1X,'directory (this is ',
     +            'not a problem nor an error). The default values ',
     +            'for the',/,1X,'ALGENCAN parameters will be used.')
 9005 format(/,1X,'Specification file algencan.dat is being used.')
 9010 format(/,1X,'Error reading specification file algencan.dat.')
 9020 format(  1X,'Ignoring unknown keyword ',A38)
 9030 format(  1X,'Ignoring incomplete keyword ',A38)
 9040 format(1X,A38)
 9041 format(1X,A38,5X,I20)
 9042 format(1X,A38,1X,1P,D24.8)

 9050 format(1X,'(Subroutine ',A4,' from HSL will be used as a direct ',
     +          'solver for Newtonian',/,1X,'linear systems.)')
 9060 format(1X,'(Subroutine ',A4,' from HSL will be used as a direct ',
     +          'solver for the',/,1X,'acceleration-step linear ',
     +          'systems.)')
 9070 format(1X,'(Linear systems will be scaled using the embedded ',
     +          'scaling option of',/,1X,'subroutine ',A4,' from HSL.)')
 9080 format(1X,'(Subroutine ',A4,' from HSL will be used for scaling ',
     +          'linear systems.)')

 9100 format(/,1X,'Warning: Ignoring keyword ',A38,'.',
     +       /,1X,'This option requires subroutines EVALH and EVALHC, ',
     +            'or, alternatively,',/,1X,'subroutine EVALHLP, to ',
     +            'be coded by the user. If you already coded them,',
     +       /,1X,'set array CODED in subrutine INIP appropiately.',/)
 9110 format(/,1X,'Warning: Ignoring keyword ',A38,'.',
     +       /,1X,'This option requires subroutines EVALH and EVALHC, ',
     +            'or, alternatively,',/,1X,'subroutine EVALHL, to ',
     +            'be coded by the user. If you already coded them,',
     +       /,1X,'set array CODED in subrutine INIP appropiately.',/)
 9120 format(/,1X,'Warning: Ignoring keyword ',A38,'.',
     +       /,1X,'This option requires subroutine MA27 or MA57 from ',
     +            'HSL to be provided by the',/,1X,'user. If you ',
     +            'have any of them, see the compilation instructions ',
     +            'for details.',/)
 9130 format(/,1X,'Warning: Ignoring keyword ',A38,'.',
     +       /,1X,'This option requires subroutine MA57 (that has an ',
     +            'embedded scaling option)',/,1X,'or, to be used in ',
     +            'connection with subroutine MA27, subroutine MC30 ',
     +            'or MC77',/,1X,'from HSL to be provided by the ',
     +            'user. If you have any of them, see the',/,1X,
     +            'compilation instructions for details.',/)

      end
