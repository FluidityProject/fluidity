
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
!    amcgsoftware@imperial.ac.uk
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


module printout

  use fldebug

  implicit none

contains

  subroutine open_output_file(output_channel,output_name,output_name_length,file_format)

    integer,                           intent(in) :: output_channel,output_name_length
    character(len=12),                 intent(in) :: file_format    
    character(len=output_name_length), intent(in) :: output_name

    ! local variables
    integer :: ierror

    ewrite(3,*) 'In open_output_file'

    open(output_channel,file=trim(output_name),status='replace',form=trim(file_format),action='write',iostat=ierror) 

    if (ierror .ne. 0) then       
       write(*,*) 'Problem opening output file ',trim(output_name)  
       stop 446       
    end if ! if (ierror .ne. 0) then      

    ewrite(3,*) 'Leaving open_output_file'

  end subroutine open_output_file

  !----------------------------------------------------------------

  subroutine write_integer_to_string(integer_variable,string_variable,len_string_variable)

    ! write an integer to a string

    integer, intent(in) :: integer_variable,len_string_variable
    character(len=len_string_variable), intent(inout) :: string_variable

    ewrite(3,*) 'In write_integer_to_string'

    if (integer_variable .lt. 10) then

       if (len_string_variable .lt. 1) then

          write(*,*) 'Write integer to string error - string to short,len_string_variable,integer_variable', &
               len_string_variable,integer_variable 
          stop 29134

       end if ! if (len_string_variable .lt. 1) then

       write(unit=string_variable,fmt='(I1)') integer_variable

    else if ((integer_variable .ge. 10) .and. (integer_variable .lt. 100)) then

       if (len_string_variable .lt. 2) then

          write(*,*) 'Write integer to string error - string to short,len_string_variable,integer_variable', &
               len_string_variable,integer_variable 
          stop 29134

       end if ! if (len_string_variable .lt. 2) then

       write(unit=string_variable,fmt='(I2)') integer_variable

    else if ((integer_variable .ge. 100) .and. (integer_variable .lt. 1000)) then

       if (len_string_variable .lt. 3) then

          write(*,*) 'Write integer to string error - string to short,len_string_variable,integer_variable', &
               len_string_variable,integer_variable 
          stop 29134

       end if ! if (len_string_variable .lt. 3) then

       write(unit=string_variable,fmt='(I3)') integer_variable

    else if ((integer_variable .ge. 1000) .and. (integer_variable .lt. 10000)) then

       if (len_string_variable .lt. 4) then

          write(*,*) 'Write integer to string error - string to short,len_string_variable,integer_variable', &
               len_string_variable,integer_variable 
          stop 29134

       end if ! if (len_string_variable .lt. 4) then

       write(unit=string_variable,fmt='(I4)') integer_variable

    else if ((integer_variable .ge. 10000) .and. (integer_variable .lt. 100000)) then

       if (len_string_variable .lt. 5) then

          write(*,*) 'Write integer to string error - string to short,len_string_variable,integer_variable', &
               len_string_variable,integer_variable 
          stop 29134

       end if ! if (len_string_variable .lt. 5) then

       write(unit=string_variable,fmt='(I5)') integer_variable

    else if ((integer_variable .ge. 100000) .and. (integer_variable .lt. 1000000)) then

       if (len_string_variable .lt. 6) then

          write(*,*) 'Write integer to string error - string to short,len_string_variable,integer_variable', &
               len_string_variable,integer_variable 
          stop 29134

       end if ! if (len_string_variable .lt. 6) then

       write(unit=string_variable,fmt='(I6)') integer_variable

    else if ((integer_variable .ge. 1000000) .and. (integer_variable .lt. 10000000)) then

       if (len_string_variable .lt. 7) then

          write(*,*) 'Write integer to string error - string to short,len_string_variable,integer_variable', &
               len_string_variable,integer_variable 
          stop 29134

       end if ! if (len_string_variable .lt. 7) then

       write(unit=string_variable,fmt='(I7)') integer_variable

    else if ((integer_variable .ge. 10000000) .and. (integer_variable .lt. 100000000)) then

       if (len_string_variable .lt. 8) then

          write(*,*) 'Write integer to string error - string to short,len_string_variable,integer_variable', &
               len_string_variable,integer_variable 
          stop 29134

       end if ! if (len_string_variable .lt. 8) then

       write(unit=string_variable,fmt='(I8)') integer_variable

    else

       write(*,*) 'Write integer to string error - integer to big,integer_variable',integer_variable
       stop 29135

    end if ! if (integer_variable .lt. 10) then

    ewrite(3,*) 'Leaving write_integer_to_string'

  end subroutine write_integer_to_string

  ! -----------------------------------------------------------------------------------------

  subroutine output_fem_sol_of_cv( unit, totele, cv_nonods, x_nonods, nphase, &
       cv_nloc, x_nloc, cv_ndgln, x_ndgln, x, femt )

    implicit none

    integer, intent( in ) :: unit, totele, cv_nonods, x_nonods, nphase, cv_nloc, x_nloc
    integer, dimension( totele * cv_nloc ), intent( in ) :: cv_ndgln
    integer, dimension( totele * x_nloc ), intent( in ) :: x_ndgln
    real, dimension( x_nonods ), intent( in ) :: x
    real, dimension( cv_nonods * nphase ), intent( in ) :: femt

    ! Local variables
    integer :: ele, cv_iloc

    ewrite(3,*) 'In output_fem_sol_of_cv'

    do ele = 1, totele
       do cv_iloc = 1, cv_nloc
          write( unit ,* ) x( x_ndgln(( ele - 1 ) * x_nloc + cv_iloc )), &
               femt( cv_ndgln(( ele - 1 ) * cv_nloc + cv_iloc ))
       end do
    end do

    ewrite(3,*) 'Leaving output_fem_sol_of_cv'

  end subroutine output_fem_sol_of_cv

  ! -----------------------------------------------------------------------------------------

  subroutine generate_name_dump( itime, unit, field, field_no )
    implicit none
    integer, intent( in ) :: itime
    integer, intent( inout ) :: unit
    character( len = 50 ), intent( in ) :: field
    integer, intent( in ) :: field_no

    ! Local variables
    character( len = 50 ) :: file_name_in, file_name_out, dump
    integer :: iaux, k, k1, k2, k3

    ewrite(3,*) 'In generate_name_dump'

    iaux = 9997
    file_name_in = 'test_'
    k1 = index( file_name_in, ' ' ) - 1 
    k2 = index( field, ' ' ) - 1 
    dump = '.d.'
    k3 =  index( dump, ' ' ) - 1 
    file_name_in = file_name_in( 1 : k1 ) // field( 1 : k2 ) // dump( 1 : k3 )
    !file_name_in = trim( trim( file_name_in )//trim( field )//trim( dump )

    open( iaux, file = 'tempfile', status = 'unknown' )
    k = index( file_name_in, ' ' ) - 1
    write( iaux, 222 ) file_name_in( 1 : k ), itime
    !trim( file_name_in ), '.d.', itime

    rewind( unit = iaux )
    read( iaux, * ) file_name_out
    k = index( file_name_out, ' ' ) - 1
    !unit = itime + field_no
    unit = 100000 + field_no
    open( unit, file = trim( file_name_out ), status = 'unknown' )

    close( iaux )


222 format( a, i0 )
    !222 format( a, a3, i0 )

    ewrite(3,*) 'Leaving generate_name_dump'

  end subroutine generate_name_dump

  ! -----------------------------------------------------------------------------------------

  subroutine printing_field_array( unit, totele, &
       cv_nonods, x_nonods, x_nloc, x_ndgln, cv_nloc, cv_ndgln, &
       pos_x, field_length, field, iphase )
    implicit none
    integer, intent( in ) :: unit, totele, x_nloc, cv_nonods, x_nonods, cv_nloc, iphase
    integer, dimension( totele * x_nloc ), intent( in ) :: x_ndgln
    integer, dimension( totele * cv_nloc ), intent( in ) :: cv_ndgln
    real, dimension( x_nonods ), intent( in ) :: pos_x
    integer, intent( in ) :: field_length
    real, dimension( field_length ), intent( in ) :: field

    ! Local variables
    integer :: cv_iloc, ele, xi_nod, xi_nod_plus, xi_nod_minus, field_nod
    real :: x_coord

    ewrite(3,*) 'In printing_field_array'

    Loop_Elements: do ele = 1, totele

       cv_iloc = 1
       xi_nod =      x_ndgln(( ele - 1 ) * x_nloc  + cv_iloc )
       xi_nod_plus = x_ndgln(( ele - 1 ) * x_nloc  + cv_iloc + 1 )
       field_nod = cv_ndgln(( ele - 1 ) * cv_nloc + cv_iloc ) + ( iphase - 1 ) * cv_nonods
       x_coord = pos_x( xi_nod ) + pos_x( xi_nod_plus )
       x_coord =  0.5 * x_coord
       write( unit, * ) pos_x( xi_nod ), field( field_nod )
       write( unit, * ) x_coord, field( field_nod )

       do cv_iloc = 2, cv_nloc - 1

          xi_nod =       x_ndgln(( ele - 1 ) * x_nloc  + cv_iloc )
          xi_nod_plus  = x_ndgln(( ele - 1 ) * x_nloc  + cv_iloc + 1 )
          xi_nod_minus = x_ndgln(( ele - 1 ) * x_nloc  + cv_iloc - 1 )
          x_coord = pos_x( xi_nod ) + pos_x( xi_nod_minus )
          x_coord = 0.5 * x_coord
          field_nod = cv_ndgln(( ele - 1 ) * cv_nloc + cv_iloc ) + ( iphase - 1 ) * cv_nonods
          write( unit, * ) x_coord, field( field_nod )
          x_coord = pos_x( xi_nod ) + pos_x( xi_nod_plus )
          x_coord =  0.5 * x_coord
          write( unit, * ) x_coord, field( field_nod )

       end do

       cv_iloc = cv_nloc
       xi_nod =       x_ndgln(( ele - 1 ) * x_nloc  + cv_iloc )
       xi_nod_minus = x_ndgln(( ele - 1 ) * x_nloc  + cv_iloc - 1 )
       x_coord = pos_x( xi_nod ) + pos_x( xi_nod_minus )
       x_coord = 0.5 * x_coord
       field_nod = cv_ndgln(( ele - 1 ) * cv_nloc + cv_iloc ) + ( iphase - 1 ) * cv_nonods
       write( unit, * ) x_coord, field( field_nod )
       write( unit, * ) pos_x( xi_nod ), field( field_nod )

    end do Loop_Elements

    ewrite(3,*) 'Leaving printing_field_array'

  end subroutine printing_field_array


  ! -----------------------------------------------------------------------------------------

  subroutine printing_field_array_veloc( unit, totele, &
       cv_nonods, x_nonods, x_nloc, x_ndgln, cv_nloc, cv_ndgln, &
       pos_x,  field)
    implicit none
    integer, intent( in ) :: unit, totele, x_nloc, cv_nonods, x_nonods, cv_nloc
    integer, dimension( totele * x_nloc ), intent( in ) :: x_ndgln
    integer, dimension( totele * cv_nloc ), intent( in ) :: cv_ndgln
    real, dimension( x_nonods ), intent( in ) :: pos_x
    real, dimension(cv_nloc*totele), intent( in ) :: field

    ! Local variables
    integer :: cv_iloc, ele, xi_nod, xi_nod_plus, xi_nod_minus, field_nod
    real :: x_coord

    ewrite(3,*) 'In printing_field_array'

    Loop_Elements: do ele = 1, totele

       cv_iloc = 1
       xi_nod =      x_ndgln(( ele - 1 ) * x_nloc  + cv_iloc )
       xi_nod_plus = x_ndgln(( ele - 1 ) * x_nloc  + cv_iloc + 1 )
       field_nod = ( ele - 1 ) * cv_nloc + cv_iloc 
       x_coord = pos_x( xi_nod ) + pos_x( xi_nod_plus )
       x_coord =  0.5 * x_coord
       write( unit, * ) pos_x( xi_nod ), field( field_nod )
       write( unit, * ) x_coord, field( field_nod )

       do cv_iloc = 2, cv_nloc - 1

          xi_nod =       x_ndgln(( ele - 1 ) * x_nloc  + cv_iloc )
          xi_nod_plus  = x_ndgln(( ele - 1 ) * x_nloc  + cv_iloc + 1 )
          xi_nod_minus = x_ndgln(( ele - 1 ) * x_nloc  + cv_iloc - 1 )
          x_coord = pos_x( xi_nod ) + pos_x( xi_nod_minus )
          x_coord = 0.5 * x_coord
          field_nod = ( ele - 1 ) * cv_nloc + cv_iloc
          write( unit, * ) x_coord, field( field_nod )
          x_coord = pos_x( xi_nod ) + pos_x( xi_nod_plus )
          x_coord =  0.5 * x_coord
          write( unit, * ) x_coord, field( field_nod )

       end do

       cv_iloc = cv_nloc
       xi_nod =       x_ndgln(( ele - 1 ) * x_nloc  + cv_iloc )
       xi_nod_minus = x_ndgln(( ele - 1 ) * x_nloc  + cv_iloc - 1 )
       x_coord = pos_x( xi_nod ) + pos_x( xi_nod_minus )
       x_coord = 0.5 * x_coord
       field_nod = ( ele - 1 ) * cv_nloc + cv_iloc
       write( unit, * ) x_coord, field( field_nod )
       write( unit, * ) pos_x( xi_nod ), field( field_nod )

    end do Loop_Elements

    ewrite(3,*) 'Leaving printing_field_array'

  end subroutine printing_field_array_veloc
  ! -----------------------------------------------------------------------------------------

  subroutine printing_veloc_field( unit, totele, &
       xu_nonods, xu_nloc, xu_ndgln, u_nloc, u_ndgln, &
       pos_x, u_nonods, field_length, field, iphase )
    implicit none
    integer, intent( in ) :: unit, totele, xu_nonods, xu_nloc
    integer, dimension( totele * xu_nloc ), intent( in ) :: xu_ndgln
    integer, intent( in ) :: u_nloc
    integer, dimension( totele * u_nloc ), intent( in ) :: u_ndgln
    real, dimension( xu_nonods ), intent( in ) :: pos_x
    integer, intent( in ) :: u_nonods
    integer, intent( in ) :: field_length
    real, dimension( field_length ), intent( in ) :: field
    integer, intent( in ) :: iphase

    ! Local variables
    integer :: x_iloc, ele, xi_nod, xi_nod_plus, &
         xi_nod_minus, field_nod
    real :: x_coord

    ewrite(3,*) 'In printing_veloc_field'
    
    Loop_Elements: do ele = 1, totele

       x_iloc = 1
       xi_nod =      xu_ndgln(( ele - 1 ) * xu_nloc  + x_iloc )
       xi_nod_plus = xu_ndgln(( ele - 1 ) * xu_nloc  + x_iloc + 1 )
       field_nod = u_ndgln(( ele - 1 ) * u_nloc + x_iloc ) + ( iphase - 1 ) * u_nonods
       x_coord = pos_x( xi_nod ) + pos_x( xi_nod_plus )
       x_coord =  0.5 * x_coord
       write( unit, * ) pos_x( xi_nod ), field( field_nod )
       write( unit, * ) x_coord, field( field_nod )

       do x_iloc = 2, xu_nloc - 1

          xi_nod =       xu_ndgln(( ele - 1 ) * xu_nloc  + x_iloc )
          xi_nod_plus  = xu_ndgln(( ele - 1 ) * xu_nloc  + x_iloc + 1 )
          xi_nod_minus = xu_ndgln(( ele - 1 ) * xu_nloc  + x_iloc - 1 )
          x_coord = pos_x( xi_nod ) + pos_x( xi_nod_minus )
          x_coord = 0.5 * x_coord
          field_nod = u_ndgln(( ele - 1 ) * u_nloc + x_iloc ) + ( iphase - 1 ) * u_nonods ! Is this correct? check later!
          write( unit, * ) x_coord, field( field_nod )
          x_coord = pos_x( xi_nod ) + pos_x( xi_nod_plus )
          x_coord =  0.5 * x_coord
          write( unit, * ) x_coord, field( field_nod )

       end do

       x_iloc = xu_nloc
       xi_nod =       xu_ndgln(( ele - 1 ) * xu_nloc  + x_iloc )
       xi_nod_minus = xu_ndgln(( ele - 1 ) * xu_nloc  + x_iloc - 1 )
       x_coord = pos_x( xi_nod ) + pos_x( xi_nod_minus )
       x_coord = 0.5 * x_coord
       field_nod = u_ndgln(( ele - 1 ) * u_nloc + x_iloc ) + ( iphase - 1 ) * u_nonods
       write( unit, * ) x_coord, field( field_nod )
       write( unit, * ) pos_x( xi_nod ), field( field_nod )

    end do Loop_Elements

    ewrite(3,*) 'Leaving printing_veloc_field'

  end subroutine printing_veloc_field

  ! -----------------------------------------------------------------------------------------

  subroutine printing_fw_field( unit, totele, &
       cv_nonods, cv_nloc, cv_ndgln, &
       field_length, field )
    implicit none
    integer, intent( in ) :: unit, totele, cv_nonods, cv_nloc
    integer, dimension( totele * cv_nloc ), intent( in ) :: cv_ndgln
    integer, intent( in ) :: field_length
    real, dimension( field_length ), intent( in ) :: field

    ! Local variables
    integer :: cv_iloc, ele, cv_nod
    real :: visc1, visc2, s_gc, s_or, fw, kr1, kr2, sat

    ewrite(3,*) 'In printing_fw_field'

    Loop_ELE: DO ELE = 1, TOTELE

       Loop_CVNLOC: DO CV_ILOC = 1, CV_NLOC

          CV_NOD = CV_NDGLN(( ELE - 1) * CV_NLOC + CV_ILOC )

          VISC1 = 0.4E-2
          VISC2 = 0.5E-2
          S_GC = 0.2 ! here S_gc --> S_wi
          S_OR = 0.3

          SAT = max( 0.0, FIELD( CV_NOD + CV_NONODS )) ! Second phase
          KR2 = ( 1. - ( 1. - SAT - S_GC ) / MAX( 1.E-5, 1. - S_GC - S_OR )) **2
          KR2 = KR2 / REAL( CV_NLOC )

          SAT = MAX( 0.0, FIELD( CV_NOD )) ! First phase
          KR1 = 0.7 * (( SAT - S_GC ) / MAX( 1.E-5, 1. - S_GC - S_OR )) **2
          KR1 = KR1 / REAL( CV_NLOC )

          FW = FW + 1. / ( 1. + VISC1 / KR1 * KR2 / VISC2 )
          FW = FW / REAL( CV_NLOC )

       END DO Loop_CVNLOC

       write( unit, * ) SAT, KR1

    END DO Loop_ELE

    ewrite(3,*) 'Leaving printing_fw_field'

  end subroutine printing_fw_field

  ! -----------------------------------------------------------------------------------------


end module printout
