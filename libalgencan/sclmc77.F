#include "mc77ad.f"

C     sclinfo:
C
C     0: Success.
C     1: Invalid row or column index.
C     2: Duplicate entry.
C     6: Insufficient space to store the input matrix.
C     7: Insufficient double precision working space.

C     ******************************************************************
C     ******************************************************************

      logical function scl(sclsub)

C     SCALAR ARGUMENTS
      character * 4 sclsub

      sclsub = 'MC77'
      scl = .true.
      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sclini()

      implicit none

#include "dim.par"
#include "mc77dat.com"

      call mc77id(icntl,cntl)

C     Suppress monitoring, warning and error messages
      icntl(1) = -1
      icntl(2) = -1

C     Avoid checking input matrix indices
      icntl(4) =  1

C     Input matrix is symmetric
      icntl(6) =  1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sclana(nsys,hnnz,hlin,hcol,hval,hdiag,sclinfo)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz,sclinfo

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz),hdiag(nsys)
      double precision hval(hnnz)

#include "dim.par"
#include "outtyp.com"
#include "mc77dat.com"

C     LOCAL SCALARS
      logical dupent
      integer i,j,k

C     LOCAL ARRAYS
      integer hlen(nsysmax),hsta(nsysmax)

      if ( nsys .gt. nsysmax ) then
        ! INSUFFICIENT SPACE TO STORE THE INPUT MATRIX

          if ( iprintctl(3) ) then
              write(* ,9000) nsysmax,nsys
              write(file10_unit,9000) nsysmax,nsys
          end if

          sclinfo = 6
          return

      end if

      do k = 1,hnnz
         i = hlin(k)
         j = hcol(k)

         if ( i .lt. 1 .or. i .gt. nsys .or.
     +        j .lt. 1 .or. j .gt. nsys .or.
     +        i .lt. j ) then
           ! INVALID ROW OR COLUMN INDEX

             if ( iprintctl(3) ) then
                 write(* ,9010) k,i,j
                 write(file10_unit,9010) k,i,j
             end if

             sclinfo = 1
             return

         end if
      end do

      call coo2csr(nsys,hnnz,hlin,hcol,hval,hlen,hsta)

      do j = 1,nsys
         iw(j) = 0
      end do

      dupent = .false.

      do i = 1,nsys
         do k = hsta(i),hsta(i) + hlen(i) - 1
            j = hcol(k)

            if ( iw(j) .eq. i ) then
                dupent = .true.

                if ( iprintctl(3) ) then
                   write(* ,8000) i,j
                   write(file10_unit,8000) i,j
                end if

            else
                iw(j) = i
            end if

            hlin(k) = i

            if ( i .eq. j ) then
                hdiag(i) = k
            end if
         end do
      end do


      if ( dupent ) then
        ! DUPLICATE ENTRY
          sclinfo = 2

          if ( iprintctl(3) ) then
              write(* ,9020)
              write(file10_unit,9020)
          end if

          return
      end if

C     NON-EXECUTABLE STATEMENTS

 8000 format(/,1X,'Duplicate entry with row index ',I8,' and column ',
     +            'index ',I8,'.')

 9000 format(/,1X,'SCLANA-MC77 WARNING: Insufficient space to store ',
     +            'input matrix. Increase',
     +       /,1X,'parameter nsysmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try scaling again.')
 9010 format(/,1X,'SCLANA-MC77 WARNING: Entry ',I8,' has got invalid ',
     +            'row index ',I8,' or column index ',I8,'.')
 9020 format(/,1X,'SCLANA-MC77 WARNING: Split entries are not allowed ',
     +            'by the scaling subroutine MC77.')

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sclsys(nsys,hnnz,hlin,hcol,hval,s,sclinfo)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz,sclinfo

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz)
      double precision hval(hnnz),s(nsys)

#include "dim.par"
#include "outtyp.com"
#include "mc77dat.com"

C     LOCAL SCALARS
      integer i

      call mc77bd(0,nsys,nsys,hnnz,hlin,hcol,hval,iw,nsysmax,dw,hnnzmax,
     +icntl,cntl,info,rinfo)

      if ( info(1) .eq. 0 ) then
        ! SUCCESS

          sclinfo = 0

          do i = 1,hnnz
              hval(i) = hval(i) / dw(hlin(i)) / dw(hcol(i))
          end do

          do i = 1,nsys
             s(i) = 1.0d0 / dw(i)
          end do

          return
      end if

      if ( info(1) .eq. -6 ) then
        ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

          sclinfo = 7

          if ( iprintctl(3) ) then
              write(* ,9010) hnnzmax,info(2)
              write(file10_unit,9010) hnnzmax,info(2)
          end if

          return

      end if

C     UNHANDLED ERROR

      if ( iprintctl(3) ) then
          write(* ,9020) info(1)
          write(file10_unit,9020) info(1)
      end if

      stop

C     NON-EXECUTABLE STATEMENTS

 9010 format(/,1X,'SCLSYS-MC77 WARNING: Insufficient double ',
     +            'precision working space. Increase',
     +       /,1X,'parameter hnnzmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try scaling again.')
 9020 format(/,1X,'SCLSYS-MC77 ERROR: Unhandled error ',I16,'.',
     +       /,1X,'See documentation for details.')

      end
