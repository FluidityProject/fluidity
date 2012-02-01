! Copyright (C) 2006 Imperial College London and others.
! 
! Please see the AUTHORS file in the main source directory for a full list
! of copyright holders.
! 
! Adrian Umpleby
! Applied Modelling and Computation Group
! Department of Earth Science and Engineering
! Imperial College London
! 
! adrian@Imperial.ac.uk
! 
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License.
! 
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
! USA

#include "fdebug_adapt.h"

      module flag_handling
      
        use addebug

        implicit none
  
#include "blkerr.i"
#include "blkprm.i"

        private
  
        public :: get_flag, set_flag
        
        integer, parameter, public :: biglst_flag_index = 4
        integer, parameter, public :: nodlst_flag_index = 6
        integer, parameter, public :: max_biglst_flag_bit = 29
        integer, parameter, public :: max_nodlst_flag_bit = 17

        interface get_flag
          module procedure get_biglst_flag, get_nodlst_flag
        end interface get_flag
        
        interface set_flag
          module procedure set_biglst_flag, set_nodlst_flag
        end interface set_flag

      contains

        function get_biglst_flag(biglst, iptr, ibit) result(get_flag)
!-----------------------------------------------------------------------
!
! - This subroutine checks the flag for the element or edge in block IPTR
! - to see if bit IBIT (from 1 to 29) is set (.TRUE.) or unset (.FALSE.).
!
! - Bits 1 to 10 are reserved for internal use...
!
! - Edges make use of bits 1 to 10 as follows:
! - Bit 1: set if block is not empty (cannot be set using SETFLG)
! - Bit 2: used to show geometry edges
! - Bit 3: used to show internal edges (should be overidden by bit 2)
! - Bit 4: used temporarily (by CHKFAC) for finding edges through faces
! - Bit 4: used temporarily (by ADAPT)  for finding edges through elements
! - Bit 5: used temporarily (by NDCNCT) for finding edges round a node
! - Bit 6: used to flag edges not to be touched
! - Bit 7: used (by ADAPT1) to keep a tally of elements added/removed/checked
! - Bit 8: used temporarily (by SETGMY) for calculating geometry information
!
! - Elements make use of bits 1 to 10 as follows:
! - Bit 1: set if block is not empty (cannot be set using SETFLG)
! - Bit 2: used temporarily (by EDGELS) for finding elements round an edge
! - Bit 3: used temporarily (by ADAPT1)  for finding pairs of elements
! - Bit 4: used temporarily (by NDCNCT) for finding elements round a node
! - Bit 5: used to flag elements not to be touched
! - Bit 6: used temporarily (by EDGTST) to flag surrounding elements during edge collapse
!
!-----------------------------------------------------------------------
  
          integer, dimension(:, :), intent(in) :: biglst
          integer, intent(in) :: iptr
          integer, intent(in) :: ibit
    
          logical :: get_flag
    
          integer :: next_biglst_index, flag
    
          ! Common block variables used:
          ! ierr    blkerr.i
          ! NXTBIG  blkprm.i
    
          if(ierr /= 0) then
            get_flag = .false.
            return
          end if
          
          assert(ibit > 0 .and. ibit <= max_biglst_flag_bit)  ! Invalid bit number
          assert(iptr > 0)  ! Invalid "pointer"
    
          flag = biglst(biglst_flag_index, iptr)
    
          if(flag == 0) then
            ! This flag is unset. We are only allowed to query the first bit
            ! of unset flags (the first bit indicates "set")
          
            assert(ibit == 1)  ! Block is empty
      
            get_flag = .false.
          else if(flag < 0) then
            ! A negative biglst_flag_index entry indicates that this block
            ! relates to an edge
            
            get_flag =  btest(-flag, ibit - 1)
          else
            ! A positive biglst_flag_index entry indicates that this block
            ! relates to an element. The flag lives in the third element block,
            ! so we need to skip on a few entries in biglst.
            
            next_biglst_index = biglst(NXTBIG, iptr)
            assert(biglst(biglst_flag_index, next_biglst_index) > 0) 
              ! Second block is not an element
            
            next_biglst_index = biglst(NXTBIG, next_biglst_index)
            flag = biglst(biglst_flag_index, next_biglst_index)
            assert(flag > 0)
              ! Third block is not an element
                          
            get_flag =  btest(flag, ibit - 1)
          end if
          
        end function get_biglst_flag
        
        subroutine set_biglst_flag(biglst, iptr, ibit, value)
!-----------------------------------------------------------------------
!
! - This subroutine sets/unsets bit IBIT (from 2 to 29) of the flag for the
! - element or edge in block IPTR to VALUE (.TRUE. or .FALSE.)
!
! - Bits 1 to 10 are reserved for internal use...
!
! - Edges make use of bits 1 to 10 as follows:
! - Bit 1: set if block is not empty - not allowed using this subroutine
! - Bit 2: used to show geometry edges
! - Bit 3: used to show internal edges (should be overidden by bit 2)
! - Bit 4: used temporarily (by CHKFAC) for finding edges through faces
! - Bit 4: used temporarily (by ADAPT)  for finding edges through elements
! - Bit 5: used temporarily (by NDCNCT) for finding edges round a node
! - Bit 6: used to flag edges not to be touched
! - Bit 7: used (by ADAPT1) to keep a tally of elements added/removed/checked
! - Bit 8: used temporarily (by SETGMY) for calculating geometry information
!
! - Elements make use of bits 1 to 10 as follows:
! - Bit 1: set if block is not empty (cannot be set using SETFLG)
! - Bit 2: used temporarily (by EDGELS) for finding elements round an edge
! - Bit 3: used temporarily (by ADAPT1)  for finding pairs of elements
! - Bit 4: used temporarily (by NDCNCT) for finding elements round a node
! - Bit 5: used to flag elements not to be touched
! - Bit 6: used temporarily (by EDGTST) to flag surrounding elements during edge collapse
!
!-----------------------------------------------------------------------
          integer, dimension(:, :), intent(inout) :: biglst
          integer, intent(in) :: iptr
          integer, intent(in) :: ibit
          logical, intent(in) :: value
          
          integer :: flag, next_biglst_index
          logical :: current_value
          
          ! Common block variables used:
          ! ierr    blkerr.i
          ! NXTBIG  blkprm.i
          
          if(ierr /= 0) return
          
          assert(ibit > 0 .and. ibit <= max_biglst_flag_bit)  ! Invalid bit number
          assert(iptr > 0)  ! Invalid "pointer"
          
          flag = biglst(biglst_flag_index, iptr)
          
          assert(flag /= 0)  ! Block is empty
          
          if(flag < 0) then
            ! A negative biglst_flag_index entry indicates that this block
            ! relates to an edge
                        
            flag = -flag
            
            if(btest(flag, ibit - 1) .eqv. value) return
                        
            if(value) then
              flag = ibset(flag, ibit - 1)
            else
              flag = ibclr(flag, ibit - 1)
            end if
            
            biglst(biglst_flag_index, iptr) = -flag
          else
            ! A positive biglst_flag_index entry indicates that this block
            ! relates to an element. The flag lives in the third element block,
            ! so we need to skip on a few entries in biglst.
          
            next_biglst_index = biglst(NXTBIG, iptr)
            assert(biglst(biglst_flag_index, next_biglst_index) > 0)
              ! Second block is not an element
            
            next_biglst_index = biglst(NXTBIG, next_biglst_index)
            flag = biglst(biglst_flag_index, next_biglst_index)
            assert(flag > 0)
              ! Third block is not an element
              
            if(btest(flag, ibit - 1) .eqv. value) return
              
            if(value) then
              flag = ibset(flag, ibit - 1)
            else
              flag = ibclr(flag, ibit - 1)
            end if
            
            biglst(biglst_flag_index, next_biglst_index) = flag
          end if
          
        end subroutine set_biglst_flag
        
        function get_nodlst_flag(nodlst, iptr, ibit) result(get_flag)
!-----------------------------------------------------------------------
!
! - This subroutine checks the flag for the node in block IPTR
! - to see if bit IBIT (from 1 to 17) is set (.TRUE.) or unset (.FALSE.).
!
! - Bits 1 to 10 are reserved for internal use...
!
! - These bits are used in the following ways:
! - Bit 1: set if block is not empty (cannot be set using STNDFL)
! - Bit 2: set if the node is an essential geometry node
! - Bit 3: set if the node is internal (should be overridden by bit 2)
! - Bit 4: set if the node splits a previous geometry edge
! - Bit 5: used temporarily (by ??????) to loop through nodes using edges
! - Bit 6: used temporarily (by ??????) to loop through nodes using elements
! - Bit 7: used temporarily (by EDGTST) during edge shrink evaluation
! - Bit 8: used to flag nodes not to be touched during adapt
! - Bit 9: used temporarily (by MKADPT) to flag halo nodes in a parallel problem
!
!-----------------------------------------------------------------------
        
          real, dimension(:, :), intent(in) :: nodlst
          integer, intent(in) :: iptr
          integer, intent(in) :: ibit
          
          logical :: get_flag
          
          integer :: flag
          
          ! Common block variables used:
          ! ierr    blkerr.i
          
          if(ierr /= 0) then
            get_flag = .false.
            return
          end if
          
          assert(ibit > 0 .and. ibit <= max_nodlst_flag_bit)  ! Invalid bit number
          assert(iptr > 0)  ! Invalid "pointer"
          
          flag = nodlst(nodlst_flag_index, iptr)
          
          if(flag == 0) then
            assert(ibit == 1)  ! Node block is empty
            
            get_flag = .false.
          else
            assert(flag > 0)  ! Node block is invalid
            
            get_flag = btest(flag, ibit - 1)
          end if
        
        end function get_nodlst_flag
     
        subroutine set_nodlst_flag(nodlst, iptr, ibit, value)
!-----------------------------------------------------------------------
!
! - This subroutine sets/unsets bit IBIT (from 2 to 17) of the flag for
! - the node in block IPTR to VALUE (.TRUE. or .FALSE.)
!
! - Bits 1 to 10 are reserved for internal use...
!
! - These bits are used in the following ways:
! - Bit 1: set if block is not empty - not allowed by this subroutine
! - Bit 2: set if the node is an essential geometry node
! - Bit 3: set if the node is internal (should be overridden by bit 2)
! - Bit 4: set if the node splits a previous geometry edge
! - Bit 5: used temporarily (by ??????) to loop through nodes using edges
! - Bit 6: used temporarily (by ??????) to loop through nodes using elements
! - Bit 7: used temporarily (by EDGTST) during edge shrink evaluation
! - Bit 8: used to flag nodes not to be touched during adapt
! - Bit 9: used temporarily (by MKADPT) to flag halo nodes in a parallel problem
!
!-----------------------------------------------------------------------
        
          real, dimension(:, :), intent(inout) :: nodlst
          integer, intent(in) :: iptr
          integer, intent(in) :: ibit
          logical, intent(in) :: value
        
          integer :: flag
        
          ! Common block variables used:
          ! ierr    blkerr.i
          
          if(ierr /= 0) return
      
          assert(ibit > 0 .and. ibit <= max_nodlst_flag_bit)  ! Invalid bit number
          assert(iptr > 0)  ! Invalid "pointer"
          
          flag = nodlst(nodlst_flag_index, iptr)
          
          assert(flag /= 0)  ! Node block is empty
          assert(flag > 0)  ! Node block is invalid
          
          if(btest(flag, ibit - 1) .eqv. value) return
          
          if(value) then
            flag = ibset(flag, ibit - 1)
          else
            flag = ibclr(flag, ibit - 1)
          end if
          
          nodlst(nodlst_flag_index, iptr) = flag
          
        end subroutine set_nodlst_flag
        
      end module flag_handling
