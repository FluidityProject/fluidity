#include "mc30ad.f"

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

      sclsub = 'MC30'
      scl = .true.
      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sclini()

      implicit none

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
#include "mc30dat.com"

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
         w(j) = 0
      end do

      dupent = .false.

      do i = 1,nsys
         do k = hsta(i),hsta(i) + hlen(i) - 1
            j = hcol(k)

            if ( w(j) .eq. i ) then
                dupent = .true.

                if ( iprintctl(3) ) then
                   write(* ,8000) i,j
                   write(file10_unit,8000) i,j
                end if

            else
                w(j) = i
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

 8000 format(/,' Duplicate entry with row index ',I8,' and column',
     +         ' index ',I8,'.')

 9000 format(/,1X,'SCLANA-MC30 WARNING: Insufficient space to store ',
     +            'input matrix. Increase',
     +       /,1X,'parameter nsysmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try scaling again.')
 9010 format(/,1X,'SCLANA-MC30 WARNING: Entry ',I8,' has got invalid ',
     +            'row index ',I8,' or column index ',I8,'.')
 9020 format(/,1X,'SCLANA-MC30 WARNING: Split entries are not allowed ',
     +            'by the scaling subroutine MC30.')

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
#include "mc30dat.com"

C     LOCAL SCALARS
      integer i

      call mc30ad(nsys,hnnz,hval,hlin,hcol,s,w,0,sclinfo)

      do i = 1,hnnz
          hval(i) = hval(i) * exp( s(hlin(i)) + s(hcol(i)) )
      end do

      do i = 1,nsys
         s(i) = exp( s(i) )
      end do

      end
