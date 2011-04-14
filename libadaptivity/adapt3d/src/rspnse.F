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
      LOGICAL FUNCTION RSPNSE( ASKSTR, LENASK )
C-----------------------------------------------------------------------
C
C
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER LENASK
C
      CHARACTER*120 ASKSTR
C
      LOGICAL OK
C
      INTEGER I
C
      CHARACTER*10 REPLY
C
C
  6   print*,' '
      print*, ASKSTR(1:LENASK)
      read(*,'(a)') reply
C
      i = 0
      OK = .FALSE.
  7   if( i .lt. 120 ) then
         i = i + 1
         if( reply(i:i) .eq. ' ' ) then
            goto 7
         else if( reply(i:i+2) .eq. 'YES' ) then
            OK = .TRUE.
         else if( reply(i:i+2) .eq. 'yes' ) then
            OK = .TRUE.
         else if( reply(i:i+1) .eq. 'Y ' ) then
            OK = .TRUE.
         else if( reply(i:i+1) .eq. 'y ' ) then
            OK = .TRUE.
         else if( reply(i:i+1) .eq. 'NO' ) then
            OK = .FALSE.
         else if( reply(i:i+1) .eq. 'no' ) then
            OK = .FALSE.
         else if( reply(i:i+1) .eq. 'n ' ) then
            OK = .FALSE.
         else if( reply(i:i+1) .eq. 'N ' ) then
            OK = .FALSE.
         ELSE
            GOTO 7
         END IF
      ELSE
         PRINT*,'I did not understand the reply...'
         goto 6
      END IF
C
      RSPNSE = OK
C
      RETURN
      END
C
