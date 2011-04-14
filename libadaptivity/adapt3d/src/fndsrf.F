C Copyright (C) 2006 Imperial College London and others.
C 
C Please see the AUTHORS file in the main source directory for a full list
C of copyright holders.
C 
C Adrian Umpleby
C Applied Modelling and Computation Group
C Department of Earth Science and Engineering
C Imperial College London
C 
C adrian@Imperial.ac.uk
C 
C This library is free software; you can redistribute it and/or
C modify it under the terms of the GNU Lesser General Public
C License as published by the Free Software Foundation; either
C version 2.1 of the License.
C 
C This library is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C Lesser General Public License for more details.
C 
C You should have received a copy of the GNU Lesser General Public
C License along with this library; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
C USA
      subroutine fndsrf( i1, i2, i3, i4, enlbas, enlist, nelbas, nelist,
     :                   nnod, nelm, nlist, nelst, nxtelm )
c=======================================================================
c
c - This subroutine returns an element which has nodes I1,I2,I3 and I4.
c - If NXTELM is sent with a non-zero value it will not return that element
c - but will search for another one with those nodes. If no element is found
c - NXTELM will be returned as -1.
c
c - This has been written to be used to search for 3d elements that have common
c - faces, so if only three nodes are required for a triangular face then I4
c - should be set equal to I3.
c
c - It is quite possible to use this for 2d elements (i.e. searching for
c - common lines), in which case I3 and I4 should be set equal to I2.
c
c - SERIOUS CONSEQUENCES WILL UNFOLD IF ANY OF I1,I2,I3 OR I4 ARE SET TO ZERO!!
c
c=======================================================================
      IMPLICIT NONE
C
      integer I1, I2, I3, I4, NNOD, NELM, NLIST, NELST, NXTELM,
     :        enlbas(nelm+1), enlist(nlist),
     :        nelbas(nnod+1), nelist(nelst)
C
      INTEGER IBAS1, N1, I, IELM, IBAS, NLN, ILN
c
c      print*,'Into FNDSRF: ',nxtelm
c      print*,'   connected nodes: ',i1,i2,i3,i4
      ibas1 = abs(nelbas(i1))
c
      n1 = abs(nelbas(i1+1)) - ibas1
c      print*,'   elements: ',(nelist(ibas1+i),i=1,n1)
c
      do 1000 i = 1, n1
c
         ielm = nelist(ibas1+i)
c
         if( ielm .eq. nxtelm ) goto 1000
c
         ibas = abs(enlbas(ielm))
         if( ibas .lt. 0 ) goto 1000
c
         nln  = abs(enlbas(ielm+1)) - ibas
c         print*,'   element ',ielm,(enlist(ibas+iln),iln=1,nln)
c
         do iln = 1, nln
            if( enlist(ibas+iln) .eq. i2 ) goto 10
         end do
c
         goto 1000
c
  10     do iln = 1, nln
            if( enlist(ibas+iln) .eq. i3 ) goto 20
         end do
c
         goto 1000
c
  20     if( i3 .eq. i4 ) goto 30
         do iln = 1, nln
            if( enlist(ibas+iln) .eq. i4 ) goto 30
         end do
c
         goto 1000
c
  30     nxtelm = ielm
         goto 9000
c
 1000 continue
c
      nxtelm = -1
c
 9000 continue
c      print*,'Leaving FNDSRF: ',nxtelm
c
      return
      end
c
