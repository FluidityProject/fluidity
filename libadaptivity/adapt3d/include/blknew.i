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
C
C - this common block holds the info for the dynamic list of element and edge info
C
      INTEGER MAXEDE
C
      PARAMETER( MAXEDE = 2000 )
C
      INTEGER MAXBIG, MXNODS, EMTBIG, STTBIG, ENDBIG, NELEMS,
     :        EMTNOD, STTNOD, ENDNOD, NUMNDS, NFRTND, STFRND, NEDGES,
     :        ELSLFT, LSTLFT, ANISOT, 
     :        TOPBIG, TOPNOD
C
      INTEGER*8 TELCHK, TEDCHK, TNDCHK, NELCHK, NELADD, NELSUB
C
      LOGICAL IS3DMS, EDGON, ADPTNG, USEQLY
C
      COMMON / BLKNEW / EMTBIG, STTBIG, ENDBIG, NELEMS, NEDGES,
     :                  EMTNOD, STTNOD, ENDNOD, NUMNDS, MAXBIG,
     :                  NFRTND, STFRND, IS3DMS, ELSLFT, LSTLFT,
     :                  EDGON,  ANISOT, ADPTNG, MXNODS, USEQLY,
     :                  TOPBIG, TOPNOD
C
      COMMON / BLKNW2 / TELCHK, TEDCHK, TNDCHK,
     :                  NELCHK, NELADD, NELSUB
C
