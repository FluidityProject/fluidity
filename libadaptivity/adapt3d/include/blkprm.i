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
C - this block holds the parameter info for the dynamic lists
C
      INTEGER SIZNOD, NXTNOD, LSTNOD, SIZBIG, NXTBIG, LSTBIG
C
      PARAMETER( SIZNOD = 17, NXTNOD = 5, LSTNOD = 4,
     :           SIZBIG =  7, NXTBIG = 6, LSTBIG = 5 )
C
C-----------------------------------------------------------------------
C - current map of edge block is as follows:
C
C  1,2: node pointers for this edge
C  3: not used (set to zero when edge is created)
C  4: flags (see setflg.f)
C  5: pointer to previous edge in list
C  6: pointer to next edge in list
C  7: not used (set to zero when edge is created)
C
C-----------------------------------------------------------------------
C - current map of element blocks is as follows:
C
C First block:
C
C  1,2,3: first three connected element pointers for this element
C  4: final connected element pointer, but times two plus one
C  5: pointer to previous block in list (last block for previous element)
C  6: pointer on to second block for this element
C  7: in-sphere radius times ten thousand as an integer
C
C Second block:
C
C  1,2,3,4: first four edges for this element
C  5: pointer back to first block for this element
C  6: pointer on to third block for this element
C  7: sum of edge functionals times ten thousand (as int) for edges of this element
C
C Third block:
C
C  1,2: fifth and sixth edges for this element
C  3: not used (set to zero when element is created)
C  4: flags (see setflg.f)
C  5: pointer back to second block for this element
C  6: pointer to next block in list  (first block for next element)
C  7: element region/material distinguishing value
C
C-----------------------------------------------------------------------
C - current map of node block is as follows:
C
C  1,2,3: X,Y,Z co-ordinates
C  4: pointer to previous node in list
C  5: pointer to next node in list
C  6: flags (see stndfl.f)
C  7-15: nodal metric
C  16: element number from original mesh (before adapting) that contains this node
C  17: reserved for node number when creating fixed mesh
C
C-----------------------------------------------------------------------
C
