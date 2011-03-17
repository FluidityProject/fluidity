

module printout

  implicit none

contains

  subroutine open_output_file(output_channel,output_name,output_name_length,file_format)

    integer,                           intent(in) :: output_channel,output_name_length
    character(len=12),                 intent(in) :: file_format    
    character(len=output_name_length), intent(in) :: output_name

    ! local variables
    integer :: ierror

    write(357,*) 'In open_output_file'

    open(output_channel,file=trim(output_name),status='replace',form=trim(file_format),action='write',iostat=ierror) 

    if (ierror .ne. 0) then       
       write(*,*) 'Problem opening output file ',trim(output_name)  
       stop 446       
    end if ! if (ierror .ne. 0) then      

    write(357,*) 'Leaving open_output_file'

  end subroutine open_output_file

  !----------------------------------------------------------------

  subroutine write_integer_to_string(integer_variable,string_variable,len_string_variable)

    ! write an integer to a string

    integer, intent(in) :: integer_variable,len_string_variable
    character(len=len_string_variable), intent(inout) :: string_variable

    write(357,*) 'In write_integer_to_string'

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

    write(357,*) 'Leaving write_integer_to_string'

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

    write(357,*) 'In output_fem_sol_of_cv'

    do ele = 1, totele
       do cv_iloc = 1, cv_nloc
          write( unit ,* ) x( x_ndgln(( ele - 1 ) * x_nloc + cv_iloc )), &
               femt( cv_ndgln(( ele - 1 ) * cv_nloc + cv_iloc ))
       end do
    end do

    write(357,*) 'Leaving output_fem_sol_of_cv'

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

    write(357,*) 'In generate_name_dump'

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

    write(357,*) 'Leaving generate_name_dump'

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

    write(357,*) 'In printing_field_array'

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

    write(357,*) 'Leaving printing_field_array'

  end subroutine printing_field_array

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

    write(357,*) 'In printing_veloc_field'

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

    write(357,*) 'Leaving printing_veloc_field'

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

    write(357,*) 'In printing_fw_field'

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

    write(357,*) 'Leaving printing_fw_field'

  end subroutine printing_fw_field

  ! -----------------------------------------------------------------------------------------

  subroutine check_sparsity( &
       u_pha_nonods, cv_pha_nonods, &
       u_nonods, cv_nonods, totele, &
       mx_ncolacv, ncolacv, finacv, colacv, midacv, & ! CV multi-phase eqns (e.g. vol frac, temp)
       nlenmcy, mx_ncolmcy, ncolmcy, finmcy, colmcy, midmcy, & ! Force balance plus cty multi-phase eqns
       mxnele, ncolele, midele, finele, colele, & ! Element connectivity 
       mx_ncoldgm_pha, ncoldgm_pha, coldgm_pha, findgm_pha, middgm_pha, & ! Force balance sparsity  
       mx_nct, ncolct, findct, colct, & ! CT sparsity - global cty eqn
       mx_nc, ncolc, findc, colc, & ! C sparsity operating on pressure in force balance
       mx_ncolcmc, ncolcmc, findcmc, colcmc, midcmc, & ! pressure matrix for projection method
       mx_ncolm, ncolm, findm, colm, midm )
    use spact

    implicit none
    integer, intent( in ) :: u_pha_nonods, cv_pha_nonods, u_nonods, cv_nonods, totele
    integer, intent ( in ) :: mx_ncolacv, ncolacv
    integer, dimension( cv_pha_nonods + 1 ), intent (in ) :: finacv
    integer, dimension( mx_ncolacv ), intent (in ) :: colacv
    integer, dimension( cv_pha_nonods ), intent (in ) :: midacv
    integer, intent ( in ) :: nlenmcy, mx_ncolmcy, ncolmcy
    integer, dimension( nlenmcy + 1 ), intent (in ) :: finmcy
    integer, dimension( mx_ncolmcy ), intent (in ) :: colmcy
    integer, dimension( nlenmcy ), intent (in ) :: midmcy
    integer, intent ( in ) :: mxnele, ncolele
    integer, dimension( totele ), intent (in ) :: midele
    integer, dimension( totele + 1 ), intent (in ) :: finele
    integer, dimension( mxnele ), intent (in ) :: colele
    integer, intent ( in ) :: mx_ncoldgm_pha, ncoldgm_pha
    integer, dimension( mx_ncoldgm_pha ), intent (in ) :: coldgm_pha
    integer, dimension( u_pha_nonods + 1 ), intent (in ) :: findgm_pha
    integer, dimension( u_pha_nonods ), intent (in ) :: middgm_pha
    integer, intent ( in ) :: mx_nct, ncolct
    integer, dimension( cv_nonods + 1 ), intent (in ) :: findct
    integer, dimension( mx_nct ), intent (in ) :: colct
    integer, intent ( in ) :: mx_nc, ncolc
    integer, dimension( u_nonods + 1 ), intent (in ) :: findc
    integer, dimension( mx_nc ), intent (in ) :: colc
    integer, intent ( in ) :: mx_ncolcmc, ncolcmc
    integer, dimension( cv_nonods + 1 ), intent (in ) :: findcmc
    integer, dimension( mx_ncolcmc ), intent (in ) :: colcmc
    integer, dimension( cv_nonods ), intent (in ) :: midcmc
    integer, intent ( in ) :: mx_ncolm, ncolm
    integer, dimension( cv_nonods + 1 ), intent (in ) :: findm
    integer, dimension( ncolm ), intent (in ) :: colm
    integer, dimension( cv_nonods ), intent (in ) :: midm

    ! Local variables
    integer, dimension( : ), allocatable :: dummy

    write(357,*) 'In check_sparsity'

    open( 15, file = 'CheckSparsityMatrix.dat', status = 'unknown' )
    write( 15, * )'########## FINMCY, MIDMCY, COLMCY ##################'
    write(15, * )'NCOLMCY:', NCOLMCY
    call checksparsity( .true., 15, NCOLMCY, NLENMCY, MX_NCOLMCY, FINMCY, MIDMCY, COLMCY )

    write( 15, * )'########## FINACV, COLACV, MIDACV ##################'
    write(15, * )'NCOLACV:', NCOLACV
    call checksparsity( .true., 15, NCOLACV, CV_PHA_NONODS, MX_NCOLACV, FINACV, MIDACV, COLACV  )

    write( 15, * )'########## FINELE, MIDELE, COLELE  ##################'
    write(15, * )'NCOLELE:',NCOLELE 
    call checksparsity( .true., 15, NCOLELE, TOTELE, MXNELE, FINELE, MIDELE, COLELE )

    allocate( dummy( CV_NONODS ))
    write( 15, * )'########## FINDCT, COLCT ##################'
    write(15, * )'NCOLCT:', NCOLCT
    call checksparsity( .false., 15, NCOLCT, CV_NONODS, MX_NCT, FINDCT, dummy, COLCT  )
    deallocate( dummy )

    allocate( dummy( U_NONODS ))
    write( 15, * )'########## FINDC, COLC ##################'
    write(15, * )'NCOLC:', NCOLC
    call checksparsity( .false., 15, NCOLC, U_NONODS, MX_NC, FINDC, dummy, COLC )
    deallocate( dummy )

    write( 15, * )'########## FINDGM_PHA, MIDDGM_PHA, COLDGM_PHA ##################'
    write(15, * )'NCOLDGM_PHA:',NCOLDGM_PHA 
    call checksparsity( .true., 15, NCOLDGM_PHA, U_PHA_NONODS, MX_NCOLDGM_PHA, FINDGM_PHA, MIDDGM_PHA, COLDGM_PHA )

    write( 15, * )'########## FINDCMC, MIDCMC, COLCMC ##################'
    write(15, * )'NCOLCMC:',NCOLCMC 
    call checksparsity( .true., 15, NCOLCMC, CV_NONODS, MX_NCOLCMC, FINDCMC, MIDCMC, COLCMC )

    write( 15, * )'########## FINDM, MIDM, COLM ##################'
    write(15, * )'NCOLM:',NCOLM 
    call checksparsity( .true., 15, NCOLM, CV_NONODS, MX_NCOLM, FINDM, MIDM, COLM )

    close( 15 )

    write(357,*) 'Leaving check_sparsity'

    return

  end subroutine check_sparsity

  subroutine mirror_data( unit_debug, problem, nphase, ncomp, totele, ndim, nlev, &
       u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, &
       cv_snloc,  p_snloc, stotel, &
       ncoef, nuabs_coefs, &
       u_ele_type, p_ele_type, mat_ele_type, cv_ele_type, &
       cv_sele_type, u_sele_type, &
       ntime, nits, ndpset, &
       dt, patmos, p_ini, t_ini, &
       t_beta, v_beta, t_theta, v_theta, u_theta, &
       t_disopt, u_disopt, v_disopt, t_dg_vel_int_opt, &
       u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
       domain_length, u_snloc, mat_nloc, cv_nonods, u_nonods, &
       p_nonods, mat_nonods, ncp_coefs, x_nonods, xu_nonods, &
       nlenmcy, &
       nopt_vel_upwind_coefs, &
       u_ndgln, xu_ndgln, cv_ndgln, x_ndgln, p_ndgln, &
       mat_ndgln, u_sndgln, cv_sndgln, x_sndgln, p_sndgln, &
       wic_vol_bc, wic_d_bc, wic_u_bc, wic_p_bc, wic_t_bc, & 
       suf_vol_bc, suf_d_bc, suf_cpd_bc, suf_t_bc, suf_p_bc, &
       suf_u_bc, &
       suf_u_bc_rob1, suf_u_bc_rob2, &
       opt_vel_upwind_coefs, &
       x, xu, nu, ug, &
       uabs_option, u_abs_stab, u_absorb, &
       u_source, &
       u, &
       den, satura, comp, p, cv_p, volfra_pore, perm )

    implicit none
    integer, intent( in ) :: unit_debug,problem, nphase, ncomp, totele, ndim, nlev, &
         u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, &
         cv_snloc,  p_snloc, stotel, &
         ncoef, nuabs_coefs, &
         u_ele_type, p_ele_type, mat_ele_type, cv_ele_type, &
         cv_sele_type, u_sele_type, &
         ntime, nits, ndpset, &
         t_disopt, u_disopt, v_disopt, t_dg_vel_int_opt, &
         u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
         nopt_vel_upwind_coefs, &
         u_snloc, mat_nloc, cv_nonods, u_nonods, &
         p_nonods, mat_nonods, ncp_coefs, x_nonods, xu_nonods, &
         nlenmcy
    real, intent( in ) :: dt, patmos, p_ini, t_ini, t_beta, v_beta, &
         t_theta, v_theta, u_theta, domain_length
    integer, dimension( stotel * nphase ), intent( in ) :: wic_vol_bc, wic_d_bc, wic_u_bc, wic_p_bc, wic_t_bc
    real, dimension( stotel * cv_snloc * nphase ), intent( in ) :: suf_vol_bc, suf_d_bc, suf_cpd_bc, suf_t_bc
    real, dimension( stotel * p_snloc * nphase ), intent( in ) :: suf_p_bc
    real, dimension( stotel * u_snloc * nphase ), intent( in ) :: suf_u_bc
    real, dimension( stotel * u_snloc * nphase ), intent( in ) :: suf_u_bc_rob1, suf_u_bc_rob2
    integer, dimension( totele * u_nloc ), intent( in ) :: u_ndgln
    integer, dimension( totele * xu_nloc ), intent( in ) :: xu_ndgln
    integer, dimension( totele * cv_nloc ), intent( in ) :: cv_ndgln
    integer, dimension( totele * cv_nloc ), intent( in ) :: x_ndgln
    integer, dimension( totele * p_nloc ), intent( in ) :: p_ndgln
    integer, dimension( totele * mat_nloc ), intent( in ) :: mat_ndgln
    integer, dimension( stotel * u_snloc ), intent( in ) :: u_sndgln
    integer, dimension( stotel * cv_snloc ), intent( in ) :: cv_sndgln, x_sndgln
    integer, dimension( stotel * p_snloc ), intent( in ) :: p_sndgln
    real, dimension( nopt_vel_upwind_coefs ), intent( in ) :: opt_vel_upwind_coefs
    real, dimension( x_nonods ), intent( in ) :: x
    real, dimension( xu_nonods ), intent( in ) :: xu
    real, dimension( u_nonods * nphase ), intent( in ) ::  nu, ug
    integer, dimension( nphase ), intent( in ) :: uabs_option
    real, dimension( mat_nonods, ndim * nphase, ndim * nphase ), intent( in ) :: u_abs_stab, u_absorb
    real, dimension( u_nonods * nphase ), intent( in ) :: u_source
    real, dimension( u_nonods * nphase ), intent( in ) :: u
    real, dimension( cv_nonods * nphase ), intent( in ) :: den, satura
    real, dimension( cv_nonods * nphase * ncomp ), intent( in ) :: comp 
    real, dimension( cv_nonods ), intent( in ) :: p, cv_p
    real, dimension( totele ), intent( in ) :: volfra_pore
    real, dimension( totele , ndim, ndim ), intent( in ) :: perm
    ! Local variables
    integer :: i, j, k

    write(357,*) 'In mirror_data'

    write( unit_debug, 201 ) 'problem, nphase, ncomp, totele, ndim, nlev: ', &
         problem, nphase, ncomp, totele, ndim, nlev

    write( unit_debug, 201 ) 'u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc: ' , &
         u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc

    write( unit_debug, 201 ) 'ncoef, nuabs_coefs, u_ele_type, p_ele_type, mat_ele_type: ', &
         ncoef, nuabs_coefs, u_ele_type, p_ele_type, mat_ele_type

    write( unit_debug, 201 ) 'cv_ele_type, cv_sele_type, u_sele_type, ntime, nits: ', &
         cv_ele_type, cv_sele_type, u_sele_type, ntime, nits

    write( unit_debug, 201 ) 'ndpset, t_disopt, u_disopt, v_disopt, t_dg_vel_int_opt: ', &
         ndpset, t_disopt, u_disopt, v_disopt, t_dg_vel_int_opt

    write( unit_debug, 201 ) 'u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, u_snloc, mat_nloc: ', & 
         u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, u_snloc, mat_nloc

    write( unit_debug, 201 ) 'cv_nonods, u_nonods, p_nonods, mat_nonods, ncp_coefs: ', &
         cv_nonods, u_nonods, p_nonods, mat_nonods, ncp_coefs

    write( unit_debug, 202 ) 'x_nonods, xu_nonods, nlenmcy: ', &
         x_nonods, xu_nonods, nlenmcy

    write( unit_debug, 203 ) 'dt, patmos, p_ini, t_ini: ', dt, patmos, p_ini, t_ini

    write( unit_debug, 204 ) 't_beta, v_beta, t_theta, v_theta, u_theta: ', &
         t_beta, v_beta, t_theta, v_theta, u_theta

    write( unit_debug, * ) 'domain_length: ', domain_length

    write( unit_debug, * ) '############################################'

    write( unit_debug, * ) 'wic_vol_bc( stotel * nphase ):', ( wic_vol_bc( i ), i = 1, stotel * nphase )

    write( unit_debug, * ) 'wic_d_bc( stotel * nphase ):', ( wic_d_bc( i ), i = 1, stotel * nphase )

    write( unit_debug, * ) 'wic_u_bc( stotel * nphase ):', ( wic_u_bc( i ), i = 1, stotel * nphase )

    write( unit_debug, * ) 'wic_p_bc( stotel * nphase ):', ( wic_p_bc( i ), i = 1, stotel * nphase )

    write( unit_debug, * ) 'wic_t_bc( stotel * nphase ):', ( wic_t_bc( i ), i = 1, stotel * nphase )

    write( unit_debug, * ) '############################################'

    write( unit_debug, * ) 'suf_vol_bc( stotel * cv_snloc * nphase ):', &
         ( suf_vol_bc( i ), i = 1, stotel * cv_snloc * nphase )

    write( unit_debug, * ) 'suf_d_bc( stotel * cv_snloc * nphase ):', &
         ( suf_d_bc( i ), i = 1, stotel * cv_snloc * nphase )

    write( unit_debug, * ) 'suf_cpd_bc( stotel * cv_snloc * nphase ):', &
         ( suf_cpd_bc( i ), i = 1, stotel * cv_snloc * nphase )

    write( unit_debug, * ) 'suf_t_bc( stotel * cv_snloc * nphase ):', &
         ( suf_t_bc( i ), i = 1, stotel * cv_snloc * nphase )

    write( unit_debug, * ) 'suf_p_bc( stotel * p_snloc * nphase ):', &
         ( suf_p_bc( i ), i = 1, stotel * p_snloc * nphase )

    write( unit_debug, * ) 'suf_u_bc( stotel * u_snloc * nphase ):', &
         ( suf_u_bc( i ), i = 1, stotel * u_snloc * nphase )

    write( unit_debug, * ) 'suf_u_bc_rob1( stotel * u_snloc * nphase ):', &
         ( suf_u_bc_rob1( i ), i = 1, stotel * u_snloc * nphase )

    write( unit_debug, * ) 'suf_u_bc_rob2( stotel * u_snloc * nphase ):', &
         ( suf_u_bc_rob2( i ), i = 1, stotel * u_snloc * nphase  )

    write( unit_debug, * ) 'suf_t_bc_rob1( stotel * cv_snloc * nphase ):', &
         ( suf_u_bc_rob1( i ), i = 1, stotel * cv_snloc * nphase )

    write( unit_debug, * ) 'suf_t_bc_rob2( stotel * cv_snloc * nphase ):', &
         ( suf_u_bc_rob2( i ), i = 1, stotel * cv_snloc * nphase  )

    write( unit_debug, * ) '############################################'

    write( unit_debug, * ) 'opt_vel_upwind_coefs( nopt_vel_upwind_coefs ):', &
         ( opt_vel_upwind_coefs( i ), i = 1, nopt_vel_upwind_coefs )

    write( unit_debug, * ) 'x( x_nonods ):', ( x( i ), i = 1, x_nonods )

    write( unit_debug, * ) 'xu( xu_nonods ):', ( xu( i ), i = 1, xu_nonods )

    write( unit_debug, * ) 'nu( u_nonods * nphase ):', &
         ( nu( i ), i = 1, u_nonods * nphase )

    write( unit_debug, * ) 'ug( u_nonods * nphase ):', &
         ( ug( i ), i = 1, u_nonods * nphase )

    write( unit_debug, * ) 'uabs_option( nphase ):', &
         ( uabs_option( i ), i = 1, nphase )

    write( unit_debug, * ) '############################################'

    write( unit_debug, * ) 'u_ndgln', size( u_ndgln ), ( u_ndgln( i ), i = 1, totele * u_nloc )
    write( unit_debug, * ) 'xu_ndgln', size( xu_ndgln ), ( xu_ndgln( i ), i = 1, totele * xu_nloc )
    write( unit_debug, * ) 'cv_ndgln', size( cv_ndgln ), ( cv_ndgln( i ), i = 1, totele * cv_nloc )
    write( unit_debug, * ) 'x_ndgln', size( x_ndgln ), ( x_ndgln( i ), i = 1,  totele * cv_nloc )
    write( unit_debug, * ) 'p_ndgln',  size( p_ndgln ), ( p_ndgln( i ), i = 1, totele * p_nloc )
    write( unit_debug, * ) 'mat_ndgln', size( mat_ndgln ), ( mat_ndgln( i ), i = 1, totele * mat_nloc )
    write( unit_debug, * ) 'u_sndgln', size( u_sndgln), ( u_sndgln( i ), i = 1, stotel * u_snloc )
    write( unit_debug, * ) 'cv_sndgln', size( cv_sndgln ), ( cv_sndgln( i ), i = 1, stotel * cv_snloc )
    write( unit_debug, * ) 'x_sndgln',  size( x_sndgln ),( x_sndgln( i ), i = 1, stotel * cv_snloc )
    write( unit_debug, * ) 'p_sndgln', size( p_sndgln ), ( p_sndgln( i ), i = 1, stotel * p_snloc )

    write( unit_debug, * ) '############################################'

    write( unit_debug, * ) 'u_abs_stab( mat_nonods, ndim * nphase, ndim * nphase ):'
    do i = 1, mat_nonods
       do j = 1, ndim * nphase
          write( unit_debug, * ) i , j, ( u_abs_stab( i, j, k ), k = 1, ndim * nphase )
       end do
    end do

    write( unit_debug, * ) '############################################'

    write( unit_debug, * ) 'u_absorb( mat_nonods, ndim * nphase, ndim * nphase ):'
    do i = 1, mat_nonods
       do j = 1, ndim * nphase
          write( unit_debug, * ) i , j, ( u_absorb( i, j, k ), k = 1, ndim * nphase )
       end do
    end do

    write( unit_debug, * ) '############################################'

    write( unit_debug, * ) 'u_source( u_nonods * nphase ):', &
         ( u_source( i ), i = 1, u_nonods * nphase )

    write( unit_debug, * ) '############################################'


    write( unit_debug, * ) 'u( u_nonods * nphase ):', &
         ( u( i ), i = 1, u_nonods * nphase )

    write( unit_debug, * ) 'den( cv_nonods * nphase ):', &
         ( den( i ), i = 1, cv_nonods * nphase )

    write( unit_debug, * ) 'satura( cv_nonods * nphase ):', &
         ( satura( i ), i = 1, cv_nonods * nphase )

    write( unit_debug, * ) 'comp( cv_nonods * nphase * ncomp ):', &
         ( comp( i ), i = 1, cv_nonods * nphase * ncomp )

    write( unit_debug, * ) 'p( cv_nonods ):', &
         ( p( i ), i = 1, cv_nonods )

    write( unit_debug, * ) 'cv_p( cv_nonods ):', &
         ( cv_p( i ), i = 1, cv_nonods )

    write( unit_debug, * ) 'volfra_pore( totele ):', &
         ( volfra_pore( i ), i = 1, totele )

    write( unit_debug, * ) 'perm( totele, ndim, ndim ):'
    do i = 1, totele
       do j = 1, ndim
          write( unit_debug, * ) i , j, ( perm( i, j, k ), k = 1, ndim )
       end do
    end do

    !write( unit_debug, * ) '(  ):', &
    !     ( ( i ), i = 1,  )

    write( unit_debug, * ) '############################################'


201 format( a, 1x, 5i5 )
202 format( a, 1x, 3i5 )
203 format( a, 1x, 4g10.4 )
204 format( a, 1x, 5g10.4 )

    write(357,*) 'Leaving mirror_data'

    return
  end subroutine mirror_data

  subroutine mirror_array_int( unit_debug, name, ndim, array )
    implicit none
    integer, intent( in ) :: unit_debug
    character( len = 50 ), intent( in ) :: name
    integer, intent( in ) :: ndim
    integer, dimension( ndim ), intent( in ) :: array

    ! Local variables
    integer :: idim, k

    k = index( name, ' ' ) - 1
    write( unit_debug, * ) name( 1 : k ), ( array( idim ), idim = 1, ndim )

    return
  end subroutine mirror_array_int


  subroutine mirror_array_real( unit_debug, name, ndim, array )
    implicit none
    integer, intent( in ) :: unit_debug
    character( len = 50 ), intent( in ) :: name
    integer, intent( in ) :: ndim
    real, dimension( ndim ), intent( in ) :: array

    ! Local variables
    integer :: idim, k

    k = index( name, ' ' ) - 1
    write( unit_debug, * ) name( 1 : k ), ( array( idim ), idim = 1, ndim )

    return
  end subroutine mirror_array_real


  subroutine mirror_matrix_int( unit_debug, name, ndim1, ndim2, array )
    implicit none
    integer, intent( in ) :: unit_debug
    character( len = 50 ), intent( in ) :: name
    integer, intent( in ) :: ndim1, ndim2
    integer, dimension( ndim1, ndim2 ), intent( in ) :: array

    ! Local variables
    integer :: idim1, idim2, k

    k = index( name, ' ' ) - 1

    write( unit_debug, * ) name( 1 : k )
    do idim1 = 1, ndim1
       write( unit_debug, * ) ( array( idim1, idim2 ), idim2 = 1, ndim2 )
    end do

    return
  end subroutine mirror_matrix_int

  subroutine mirror_matrix_real( unit_debug, name, ndim1, ndim2, array )
    implicit none
    integer, intent( in ) :: unit_debug
    character( len = 50 ), intent( in ) :: name
    integer, intent( in ) :: ndim1, ndim2
    real, dimension( ndim1, ndim2 ), intent( in ) :: array

    ! Local variables
    integer :: idim1, idim2, k

    k = index( name, ' ' ) - 1

    write( unit_debug, * ) name( 1 : k )
    do idim1 = 1, ndim1
       write( unit_debug, * ) ( array( idim1, idim2 ), idim2 = 1, ndim2 )
    end do

    return
  end subroutine mirror_matrix_real

end module printout



module input_var

contains


  subroutine get_entry( unit, len_name, ior, ifile, fcn_name, value_real, value_bool )
    implicit none
    integer, intent( in ) :: unit, len_name
    integer, intent( inout ) :: ior ! ior=13 indicates a non-scalar, ie, array, matrix or tensor
    character( len = len_name ), intent( inout ) :: ifile, fcn_name
    real, intent( inout ) :: value_real
    logical, intent( inout ) :: value_bool
    ! Local variables
    character( len = len_name ) :: name
    integer, parameter :: nexit = 4
    character( len = nexit ), parameter :: exit_file = 'exit'
    integer :: k
    real :: dummy_real

    !write(357,*) 'In get_entry'

    ior = 1
    ifile = ' '

    value_real = 0.
    value_bool = .true.

    read( unit, * ) name
    write( 357, * ) name
    if( ( name( 1 : 1 ) == '#' ) .or. ( name( 1 : 1 ) == ' ' ) ) then
       ior = 0
    else
       k = index( name, ' ' ) - 1
       if( name( 1 : nexit ) == exit_file ) then 
          ior = -1000
          return
       elseif(( name( 1 : 9 ) == 'lump_eqns' ) .or. &
           ( name( 1 : 21 ) == 'volfra_use_theta_flux' ) .or. ( name( 1 : 21 ) == 'volfra_get_theta_flux' ) .or. &
           ( name( 1 : 19 ) == 'comp_use_theta_flux' ) .or. ( name( 1 : 19 ) == 'comp_get_theta_flux' )) then
          backspace( unit )
          read( unit, * ) name, value_bool
          k = index( name, ' ' ) - 1
          ifile( 1 : k ) = name( 1 : k )
          ior = 5
       else
          backspace( unit )
          read( unit, * ) name, value_real

          if( value_real < -1000. ) then
             backspace( unit )
             read( unit, * ) name, dummy_real, fcn_name
             k = index( name, ' ' ) - 1
             ifile( 1 : k ) = name( 1 : k )
             fcn_name = trim( fcn_name )
             ior = 13
          elseif( abs( value_real + 1000. ) < 1.e-6 ) then ! so value is -1000
 ! Real or integer array
             backspace( unit )
             read( unit, * ) name, dummy_real, fcn_name
             k = index( name, ' ' ) - 1
             ifile( 1 : k ) = name( 1 : k )
             fcn_name = trim( fcn_name )
             ior = 15
          elseif( abs( value_real + 999. ) < 1.e-6 ) then ! so value is -999
! 2 x 2 matrix
             backspace( unit )
             read( unit, * ) name, dummy_real, fcn_name
             k = index( name, ' ' ) - 1
             ifile( 1 : k ) = name( 1 : k )
             fcn_name = trim( fcn_name )
             ior = 16
          elseif( abs( value_real + 998. ) < 1.e-6 ) then ! so value is -998
! 3 x 3 matrix
             backspace( unit )
             read( unit, * ) name, dummy_real, fcn_name
             k = index( name, ' ' ) - 1
             ifile( 1 : k ) = name( 1 : k )
             fcn_name = trim( fcn_name )
             ior = 17
          elseif( abs( value_real + 997  ) < 1.e-6 ) then ! so value is -997
! 4 x 4 matrix
             backspace( unit )
             read( unit, * ) name, dummy_real, fcn_name
             k = index( name, ' ' ) - 1
             ifile( 1 : k ) = name( 1 : k )
             fcn_name = trim( fcn_name )
             ior = 18
          else
             k = index( name, ' ' ) - 1
             ifile( 1 : k ) = name( 1 : k )
             ior = 3
          end if

       end if
    end if

    return
  end subroutine get_entry

  subroutine Scalar_Int_Val( unit, len_name, fcn_name, scalar_int )
    implicit none
    integer, intent( in ) :: unit, len_name
    character( len = len_name ), intent( in ) :: fcn_name
    integer, intent ( inout ) :: scalar_int
    ! Local variables
    integer :: idim, k

    k = index( fcn_name, ' ' ) - 1
    call system( 'rm -f fvalues' )
    call system( './' // fcn_name( 1 : k ) // ' > fvalues ' )
    open( unit + 1, file = 'fvalues', status = 'unknown' )

    read( unit + 1, * ) scalar_int

    close( unit + 1 )

    return
  end subroutine Scalar_Int_Val

  subroutine Scalar_Real_Val( unit, len_name, fcn_name, scalar_real )
    implicit none
    integer, intent( in ) :: unit, len_name
    character( len = len_name ), intent( in ) :: fcn_name
    real, intent ( inout ) :: scalar_real
    ! Local variables
    integer :: idim, k

    k = index( fcn_name, ' ' ) - 1
    call system( 'rm -f fvalues' )
    call system( './' // fcn_name( 1 : k ) // ' > fvalues ' )
    open( unit + 1, file = 'fvalues', status = 'unknown' )

    read( unit + 1, * ) scalar_real

    close( unit + 1 )

    return
  end subroutine Scalar_Real_Val


  subroutine Array_Int_Set( unit, ior, ndim, len_name, fcn_name, array )
    implicit none
    integer, intent( in ) :: unit, ior, ndim, len_name
    character( len = len_name ), intent( in ) :: fcn_name
    integer, dimension( ndim ), intent( inout ) :: array
    ! Local variables
    integer :: idim, k

    if (ior == 15 ) then

       call system( 'rm -f filedim' )
       open( unit + 1, file = 'filedim', status = 'unknown' )
       write( unit + 1, * ) ndim
       close( unit + 1 )

    endif

    k = index( fcn_name, ' ' ) - 1
    call system( 'rm -f fvalues' )
    call system( './' // fcn_name( 1 : k ) // ' > fvalues ' )
    open( unit + 1, file = 'fvalues', status = 'unknown' )

    do idim = 1, ndim
       read( unit + 1, * ) array( idim )
    end do

    close( unit + 1 )

    return
  end subroutine Array_Int_Set


  subroutine Array_Real_Set( unit, ior, ndim, len_name, fcn_name, array )
    implicit none
    integer, intent( in ) :: unit, ior, ndim, len_name
    character( len = len_name ), intent( in ) :: fcn_name
    real, dimension( ndim ), intent( inout ) :: array
    ! Local variables
    integer :: idim, k

    print*, 'ndim:', ndim
    if (ior == 15 ) then

       call system( 'rm -f filedim' )
       open( unit + 699, file = 'filedim', status = 'unknown' )
       write( unit + 699, * ) ndim
       close( unit + 699 )

    endif

    k = index( fcn_name, ' ' ) - 1
    call system( 'rm -f fvalues' )
    call system( './' // fcn_name( 1 : k ) // ' > fvalues' )
    open( unit + 699, file = 'fvalues', status = 'unknown' )

    do idim = 1, ndim
       read( unit + 699, * ) array( idim )
    end do

    close( unit + 699 )

    return
  end subroutine Array_Real_Set


  subroutine Array_Matrix2_Set( unit, ior, ndim1, ndim2, len_name, fcn_name, matrix )
    implicit none
    integer, intent( in ) :: unit, ior, ndim1, ndim2, len_name
    character( len = len_name ), intent( in ) :: fcn_name
    real, dimension( ndim1, ndim2 ), intent( inout ) :: matrix
    ! Local variables
    integer :: idim1, idim2, k

    if (ior == 16 ) then

       open( unit + 563, file = 'filedim', status = 'unknown' )
       write( unit + 563, * ) ndim1, ndim2
       close( unit + 563 )

    endif

    k = index( fcn_name, ' ' ) - 1
    call system( 'rm -f fvalues' )
    call system( './' // fcn_name( 1 : k ) // ' > fvalues ' )
    open( unit + 563, file = 'fvalues', status = 'unknown' )

    do idim1 = 1, ndim1
       do idim2 = 1, ndim2
          read( unit + 563, * ) matrix( idim1, idim2 )
       end do
    end do

    close( unit + 563 )

    return
  end subroutine Array_Matrix2_Set



  subroutine Array_Matrix3_Set( unit, ior, ndim1, ndim2, ndim3, len_name, fcn_name, matrix )
    implicit none
    integer, intent( in ) :: unit, ior, ndim1, ndim2, ndim3, len_name
    character( len = len_name ), intent( in ) :: fcn_name
    real, dimension( ndim1, ndim2, ndim3 ), intent( inout ) :: matrix
    ! Local variables
    integer :: idim1, idim2, idim3, k

    if (ior == 17 ) then

       open( unit + 1, file = 'filedim', status = 'unknown' )
       write( unit + 1, * ) ndim1, ndim2, ndim3
       close( unit + 1 )
    endif

    k = index( fcn_name, ' ' ) - 1
    call system( 'rm -f fvalues' )
    call system( './' // fcn_name( 1 : k ) // ' > fvalues ' )
    open( unit + 1, file = 'fvalues', status = 'unknown' )

    do idim1 = 1, ndim1
       do idim2 = 1, ndim2
          do idim3 = 1, ndim3
             read( unit + 1, * ) matrix( idim1, idim2, idim3 )
          end do
       end do
    end do

    close( unit + 1 )

    return
  end subroutine Array_Matrix3_Set


  subroutine Array_Matrix4_Set( unit, ior, ndim1, ndim2, ndim3, ndim4, len_name, fcn_name, matrix )
    implicit none
    integer, intent( in ) :: unit, ior, ndim1, ndim2, ndim3, ndim4, len_name
    character( len = len_name ), intent( in ) :: fcn_name
    real, dimension( ndim1, ndim2, ndim3, ndim4 ), intent( inout ) :: matrix
    ! Local variables
    integer :: idim1, idim2, idim3, idim4, k

    if (ior == 18 ) then

       open( unit + 1, file = 'filedim', status = 'unknown' )
       write( unit + 1, * ) ndim1, ndim2, ndim3, ndim4
       close( unit + 1 )

    endif

    k = index( fcn_name, ' ' ) - 1
    call system( 'rm -f fvalues' )
    call system( './' // fcn_name( 1 : k ) // ' > fvalues ' )
    open( unit + 1, file = 'fvalues', status = 'unknown' )

    do idim1 = 1, ndim1
       do idim2 = 1, ndim2
          do idim3 = 1, ndim3
             do idim4 = 1, ndim4
                read( unit + 1, * ) matrix( idim1, idim2, idim3, idim4 )
             end do
          end do
       end do
    end do

    close( unit + 1 )

    return
  end subroutine Array_Matrix4_Set


  subroutine read_scalar( unit, option_debug, problem, nphase, ncomp, totele, ndim, nlev, &
       u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, &
       cv_snloc, u_snloc, p_snloc, x_snloc, stotel, &
       ncoef, nuabs_coefs, &
       u_ele_type, p_ele_type, mat_ele_type, cv_ele_type, &
       cv_sele_type, u_sele_type, ntime, nits, nits_internal, noit_dim, &
       nits_flux_lim_volfra, nits_flux_lim_comp, &
       ndpset, &
       dt, patmos, p_ini, t_ini, &
       t_beta, v_beta, t_theta, v_theta, u_theta, &
       t_disopt, u_disopt, v_disopt, t_dg_vel_int_opt, &
       u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
       lump_eqns, & 
       volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, comp_get_theta_flux, &
       capil_pres_opt, ncapil_pres_coef, comp_diffusion_opt, ncomp_diff_coef, &
       domain_length )
    implicit none

    integer, intent( in ) :: unit
    integer, intent( inout ) :: option_debug, problem, nphase, ncomp, totele, ndim, nlev, &
         u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, &
         cv_snloc, u_snloc, p_snloc, x_snloc, stotel, &
         ncoef, nuabs_coefs, &
         u_ele_type, p_ele_type, mat_ele_type, cv_ele_type, &
         cv_sele_type, u_sele_type, ntime, nits, nits_internal, noit_dim, &
         nits_flux_lim_volfra, nits_flux_lim_comp, &
         ndpset
    real, intent( inout ) :: dt, patmos, p_ini, t_ini, &
         t_beta, v_beta, t_theta, v_theta, u_theta
    integer, intent( inout ) :: t_disopt, u_disopt, v_disopt, t_dg_vel_int_opt, &
         u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt
    logical, intent( inout ) :: lump_eqns, &
       volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, comp_get_theta_flux
    integer, intent( inout ) :: capil_pres_opt, ncapil_pres_coef, &
         comp_diffusion_opt, ncomp_diff_coef
    real, intent( inout ) :: domain_length
    ! Local variables
    integer, parameter :: len_name = 50
    character( len = len_name ) :: ifile, fcn_name
    integer :: ior, k
    real :: value_real
    logical :: value_bool

    ior = 0
    ncomp = 0

    do while ( .not. ( ior == -1000 ))
       call get_entry( unit, len_name, ior, ifile, fcn_name, value_real, value_bool )

       Conditional_IOR: if ( ior > 0 ) then
          k = index( ifile, ' ') - 1

          Select Case ( ifile( 1 : k ))

          Case( 'option_debug' );
            option_debug = real( value_real )

          Case( 'problem' );
             problem = real( value_real )

          Case( 'nphase' );
             nphase = int( value_real )

          Case( 'ncomp' );
             ncomp = int( value_real )

          Case( 'totele' );
             totele = int( value_real )

          Case( 'ndim' );
             ndim = int( value_real )

          Case( 'nlev' );
             nlev = int( value_real )

          Case( 'u_nloc' );
             u_nloc = int( value_real )

          Case( 'xu_nloc' );
             xu_nloc = int( value_real )

          Case( 'cv_nloc' );
             cv_nloc = int( value_real )

          Case( 'x_nloc' );
             x_nloc = int( value_real )

          Case( 'p_nloc' );
             p_nloc = int( value_real )

          Case( 'cv_snloc' );
             cv_snloc = int( value_real )

          Case( 'u_snloc' );
             u_snloc = int( value_real )

          Case( 'p_snloc' );
             p_snloc = int( value_real )

          Case( 'x_snloc' );
             x_snloc = int( value_real )

          Case( 'stotel' );
             stotel = int( value_real )

          Case( 'ncoef' );
             ncoef = int( value_real )

          Case( 'nuabs_coefs' );
             nuabs_coefs = int( value_real )

          Case( 'u_ele_type' );
             u_ele_type = int( value_real )

          Case( 'p_ele_type' );
             p_ele_type = int( value_real )

          Case( 'mat_ele_type' );
             mat_ele_type = int( value_real )

          Case( 'cv_ele_type' );
             cv_ele_type = int( value_real )

          Case( 'cv_sele_type' );
             cv_sele_type = int( value_real )

          Case( 'u_sele_type' );
             u_sele_type = int( value_real )

          Case( 'ntime' );
             ntime = int( value_real )

          Case( 'nits' );
             nits = int( value_real )

          Case( 'nits_internal' );
             nits_internal = int( value_real )

          Case( 'noit_dim' );
             noit_dim = int( value_real )

          Case( 'nits_flux_lim_volfra' );
             nits_flux_lim_volfra = int( value_real )

          Case( 'nits_flux_lim_comp' );
             nits_flux_lim_comp = int( value_real )

          Case( 'ndpset' );
             ndpset = int( value_real )

          Case( 'v_disopt' );
             v_disopt = int( value_real )

          Case( 't_disopt' );
             t_disopt = int( value_real )

          Case( 'u_disopt' );
             u_disopt = int( value_real )

          Case( 't_dg_vel_int_opt' );
             t_dg_vel_int_opt = int( value_real )

          Case( 'u_dg_vel_int_opt' );
             u_dg_vel_int_opt = int( value_real )

          Case( 'v_dg_vel_int_opt' );
             v_dg_vel_int_opt = int( value_real )

          Case( 'w_dg_vel_int_opt' );
             w_dg_vel_int_opt = int( value_real )

          Case( 'lump_eqns' );
             lump_eqns = value_bool

          Case( 'volfra_use_theta_flux' );
             volfra_use_theta_flux = value_bool

          Case( 'volfra_get_theta_flux' );
             volfra_get_theta_flux = value_bool

          Case( 'comp_use_theta_flux' );
             comp_use_theta_flux = value_bool

          Case( 'comp_get_theta_flux' );
             comp_get_theta_flux = value_bool

          Case( 'capil_pres_opt' );
             capil_pres_opt = int( value_real )

          Case( 'ncapil_pres_coef' );
             ncapil_pres_coef = int( value_real )

          Case( 'comp_diffusion_opt' );
             comp_diffusion_opt = int( value_real )

          Case( 'ncomp_diff_coef' );
             ncomp_diff_coef = int( value_real )

          Case( 'domain_length' );
             domain_length = int( value_real )

          Case( 'dt' );
             dt = value_real

          Case( 'patmos' );
             patmos = value_real

          Case( 'p_ini' );
             p_ini = value_real

          Case( 't_ini' );
             t_ini = value_real

          Case( 't_beta' );
             t_beta = value_real

          Case( 'v_beta' );
             v_beta = value_real

          Case( 't_theta' );
             t_theta = value_real

          Case( 'v_theta' );
             v_theta = value_real

          Case( 'u_theta' );
             u_theta = value_real

          Case DEFAULT
             write( 357, * ) 'Option not found in read_scalar subrt.', ifile( 1 : k ) 
             STOP 1877

          end Select

       end if Conditional_IOR

    end do

    return

  end subroutine read_scalar


  subroutine read_all( unit, nphase, ncomp, totele, ndim, &
       cv_snloc, u_snloc, p_snloc, stotel, &
       ncoef, nuabs_coefs, ncp_coefs, & 
       cv_nonods, u_nonods, &
       mat_nonods, x_nonods, xu_nonods, &
       nopt_vel_upwind_coefs, &
       wic_vol_bc, wic_d_bc, wic_u_bc, wic_p_bc, wic_t_bc, wic_comp_bc, & 
       suf_vol_bc, suf_d_bc, suf_cpd_bc, suf_t_bc, suf_p_bc, &
       suf_u_bc, suf_v_bc, suf_w_bc, suf_one_bc, suf_comp_bc, &
       suf_u_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob1, suf_v_bc_rob2, &
       suf_w_bc_rob1, suf_w_bc_rob2, suf_t_bc_rob1, suf_t_bc_rob2, &
       suf_vol_bc_rob1, suf_vol_bc_rob2, &
       suf_comp_bc_rob1, suf_comp_bc_rob2, &
       opt_vel_upwind_coefs, &
       volfra_error, volfra_relax, volfra_relax_diag, volfra_relax_row, volfra_relax_number_iterations, & 
       scalar_error, scalar_relax, scalar_relax_diag, scalar_relax_row, scalar_relax_number_iterations, & 
       global_error, global_relax, global_relax_diag, global_relax_row, global_relax_number_iterations, & 
       velocity_error, velocity_relax, velocity_relax_diag, velocity_relax_row, velocity_relax_number_iterations, & 
       pressure_error, pressure_relax, pressure_relax_diag, pressure_relax_row, pressure_relax_number_iterations, & 
       mass_matrix_error, mass_matrix_relax, mass_matrix_relax_diag, mass_matrix_relax_row, mass_matrix_relax_number_iterations, &
       in_ele_upwind, dg_ele_upwind, & 
       x, y, z, xu, yu, zu, nu, nv, nw, ug, vg, wg, &
       uabs_option, uabs_coefs, u_abs_stab, &
       u_absorb, t_absorb, v_absorb, comp_absorb, &
       u_source, t_source, v_source, comp_source, udiffusion, tdiffusion, &
       ncomp_diff_coef, comp_diffusion, comp_diff_coef, &
       ncapil_pres_coef, capil_pres_coef, & 
       u, v, w, &
       den, satura, comp, volfra, t, cv_one, p, cv_p, volfra_pore, perm, &
       K_Comp, alpha_beta, &
       eos_option, cp_option, eos_coefs, cp_coefs )

    implicit none
    integer, intent( in ) :: unit, nphase, ncomp, totele, ndim, &
         cv_snloc, u_snloc, p_snloc, stotel, &
         ncoef, nuabs_coefs, ncp_coefs, & 
         cv_nonods, u_nonods, &
         mat_nonods, x_nonods, xu_nonods, &
         nopt_vel_upwind_coefs
    integer, dimension( stotel * nphase ), intent( inout ) :: wic_vol_bc, wic_d_bc, wic_u_bc, wic_p_bc, wic_t_bc, wic_comp_bc 
    real, dimension( stotel * cv_snloc * nphase ), intent( inout ) :: suf_vol_bc, suf_d_bc, suf_cpd_bc, suf_t_bc
    real, dimension( stotel * p_snloc * nphase ), intent( inout ) :: suf_p_bc
    real, dimension( stotel * u_snloc * nphase ), intent( inout ) :: suf_u_bc, suf_v_bc, suf_w_bc
    real, dimension( stotel * cv_snloc * nphase ), intent( inout ) :: suf_one_bc
    real, dimension( stotel * cv_snloc * nphase * ncomp ), intent( inout ) :: suf_comp_bc
    real, dimension( stotel * u_snloc * nphase ), intent( inout ) :: suf_u_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob1, suf_v_bc_rob2, &
         suf_w_bc_rob1, suf_w_bc_rob2
    real, dimension( stotel * cv_snloc * nphase ), intent( inout ) :: suf_t_bc_rob1, suf_t_bc_rob2
    real, dimension( stotel * cv_snloc * nphase ), intent( inout ) :: suf_vol_bc_rob1, suf_vol_bc_rob2
    real, dimension( stotel * cv_snloc * nphase ), intent( inout ) :: suf_comp_bc_rob1, suf_comp_bc_rob2
    real, dimension( nopt_vel_upwind_coefs ), intent( inout ) :: opt_vel_upwind_coefs
    real, intent( inout ) :: volfra_error, volfra_relax, volfra_relax_diag, volfra_relax_row,  & 
         scalar_error, scalar_relax, scalar_relax_diag, scalar_relax_row, & 
         global_error, global_relax, global_relax_diag, global_relax_row, & 
         velocity_error, velocity_relax, velocity_relax_diag, velocity_relax_row, & 
         pressure_error, pressure_relax, pressure_relax_diag, pressure_relax_row, &
         mass_matrix_error, mass_matrix_relax, mass_matrix_relax_diag, mass_matrix_relax_row 
    integer, intent( inout ) :: volfra_relax_number_iterations, scalar_relax_number_iterations, &
         global_relax_number_iterations,  velocity_relax_number_iterations, &
         pressure_relax_number_iterations, mass_matrix_relax_number_iterations, &
         in_ele_upwind, dg_ele_upwind
    real, dimension( x_nonods ), intent( inout ) :: x, y, z
    real, dimension( xu_nonods ), intent( inout ) :: xu, yu, zu
    real, dimension( u_nonods * nphase ), intent( inout ) ::  nu, nv, nw, ug, vg, wg
    integer, dimension( nphase ), intent( inout ) :: uabs_option
    real, dimension( nphase, nuabs_coefs ), intent( inout ) :: uabs_coefs
    real, dimension( mat_nonods, ndim * nphase, ndim * nphase ), intent( inout ) :: u_abs_stab, u_absorb
    real, dimension( cv_nonods, nphase, nphase ), intent( inout ) :: t_absorb, v_absorb
    real, dimension( cv_nonods * nphase, nphase, nphase ), intent( inout ) :: comp_absorb
    real, dimension( u_nonods * nphase ), intent( inout ) :: u_source
    real, dimension( cv_nonods * nphase ), intent( inout ) :: t_source, v_source, comp_source
    real, dimension( mat_nonods, ndim, ndim, nphase ), intent( inout ) :: udiffusion, tdiffusion

    integer, intent( in ) :: ncomp_diff_coef
    real, dimension( mat_nonods, ndim, ndim, nphase ), intent( inout ) :: comp_diffusion
    real, dimension( ncomp, ncomp_diff_coef, nphase ), intent( inout ) :: comp_diff_coef
    integer, intent( in ) :: ncapil_pres_coef
    real, dimension( ncapil_pres_coef, nphase, nphase ), intent( inout ) :: &
         capil_pres_coef

    real, dimension( u_nonods * nphase ), intent( inout ) :: u, v, w
    real, dimension( cv_nonods * nphase ), intent( inout ) :: den, satura
    real, dimension( cv_nonods * nphase * ncomp ), intent( inout ) :: comp
    real, dimension( cv_nonods * nphase ), intent( inout ) :: volfra, t, cv_one
    real, dimension( cv_nonods ), intent( inout ) :: p, cv_p
    real, dimension( totele ), intent( inout ) :: volfra_pore
    real, dimension( totele , ndim, ndim ), intent( inout ) :: perm
    real, dimension( ncomp, nphase, nphase ), intent( inout ) :: K_Comp
    real, intent( inout ) :: alpha_beta
    integer, dimension( nphase ), intent( inout ) :: eos_option, cp_option
    real, dimension( nphase, ncoef ), intent( inout ) :: eos_coefs
    real, dimension( nphase, ncp_coefs ), intent( inout ) :: cp_coefs


    ! Local variables
    integer, parameter :: len_name = 50
    character( len = len_name ) :: ifile, fcn_name
    integer :: ior, k
    real :: value_real
    logical :: value_bool

    write(357,*) 'In read_all'

    ior = 0
    do while ( .not. ( ior == -1000 ))


       call get_entry( unit, len_name, ior, ifile, fcn_name, value_real, value_bool )

       Conditional_IOR: if ( ior > 0 ) then
          k = index( ifile, ' ') - 1
          Select Case ( ifile( 1 : k ))

             ! Scalar:
          Case( 'alpha_beta' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, alpha_beta )
             else
                alpha_beta = value_real
             endif

          Case( 'volfra_error' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, volfra_error )
             else
                volfra_error = value_real
             endif

          Case( 'scalar_error' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, scalar_error )
             else
                scalar_error = value_real
             endif

          Case( 'global_error' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, global_error )
             else
                global_error = value_real
             endif

          Case( 'velocity_error' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, velocity_error )
             else
                velocity_error = value_real
             endif

          Case( 'pressure_error' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, pressure_error )
             else
                pressure_error = value_real
             endif

          Case( 'mass_matrix_error' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, mass_matrix_error )
             else
                mass_matrix_error = value_real
             endif

          Case( 'volfra_relax' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, volfra_relax )
             else
                volfra_relax = value_real
             endif

          Case( 'scalar_relax' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, scalar_relax )
             else
                scalar_relax = value_real
             endif

          Case( 'global_relax' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, global_relax )
             else
                global_relax = value_real
             endif

          Case( 'velocity_relax' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, velocity_relax )
             else
                velocity_relax = value_real
             endif

          Case( 'pressure_relax' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, pressure_relax )
             else
                pressure_relax = value_real
             endif

          Case( 'mass_matrix_relax' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, mass_matrix_relax )
             else
                mass_matrix_relax = value_real
             endif

          Case( 'volfra_relax_diag' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, volfra_relax_diag )
             else
                volfra_relax_diag = value_real
             endif

          Case( 'scalar_relax_diag' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, scalar_relax_diag )
             else
                scalar_relax_diag = value_real
             endif

          Case( 'global_relax_diag' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, global_relax_diag )
             else
                global_relax_diag = value_real
             endif

          Case( 'velocity_relax_diag' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, velocity_relax_diag )
             else
                velocity_relax_diag = value_real
             endif

          Case( 'pressure_relax_diag' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, pressure_relax_diag )
             else
                pressure_relax_diag = value_real
             endif

          Case( 'mass_matrix_relax_diag' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, mass_matrix_relax_diag )
             else
                mass_matrix_relax_diag = value_real
             endif

          Case( 'volfra_relax_row' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, volfra_relax_row )
             else
                volfra_relax_row = value_real
             endif

          Case( 'scalar_relax_row' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, scalar_relax_row )
             else
                scalar_relax_row = value_real
             endif

          Case( 'global_relax_row' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, global_relax_row )
             else
                global_relax_row = value_real
             endif

          Case( 'velocity_relax_row' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, velocity_relax_row )
             else
                velocity_relax_row = value_real
             endif

          Case( 'pressure_relax_row' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, pressure_relax_row )
             else
                pressure_relax_row = value_real
             endif

          Case( 'mass_matrix_relax_row' );
             if( ior == 13 ) then
                call Scalar_Real_Val( unit, len_name, fcn_name, mass_matrix_relax_row )
             else
                mass_matrix_relax_row = value_real
             endif

          Case( 'volfra_relax_number_iterations' );
             if( ior == 13 ) then
                call Scalar_Int_Val( unit, len_name, fcn_name, volfra_relax_number_iterations )
             else
                volfra_relax_number_iterations = int( value_real )
             endif

          Case( 'scalar_relax_number_iterations' );
             if( ior == 13 ) then
                call Scalar_Int_Val( unit, len_name, fcn_name, scalar_relax_number_iterations )
             else
                scalar_relax_number_iterations = int( value_real )
             endif

          Case( 'global_relax_number_iterations' );
             if( ior == 13 ) then
                call Scalar_Int_Val( unit, len_name, fcn_name, global_relax_number_iterations )
             else
                global_relax_number_iterations = int( value_real )
             endif

          Case( 'velocity_relax_number_iterations' );
             if( ior == 13 ) then
                call Scalar_Int_Val( unit, len_name, fcn_name, velocity_relax_number_iterations )
             else
                velocity_relax_number_iterations = int( value_real )
             endif

          Case( 'pressure_relax_number_iterations' );
             if( ior == 13 ) then
                call Scalar_Int_Val( unit, len_name, fcn_name, pressure_relax_number_iterations )
             else
                pressure_relax_number_iterations = int( value_real )
             endif

          Case( 'mass_matrix_relax_number_iterations' );
             if( ior == 13 ) then
                call Scalar_Int_Val( unit, len_name, fcn_name, mass_matrix_relax_number_iterations )
             else
                mass_matrix_relax_number_iterations = int( value_real )
             endif

          Case( 'in_ele_upwind' );
             if( ior == 13 ) then
                call Scalar_Int_Val( unit, len_name, fcn_name, in_ele_upwind )
             else
                in_ele_upwind = int( value_real )
             endif

          Case( 'dg_ele_upwind' );
             if( ior == 13 ) then
                call Scalar_Int_Val( unit, len_name, fcn_name, dg_ele_upwind )
             else
                dg_ele_upwind = int( value_real )
             endif


             ! Integer arrays:
          Case( 'wic_vol_bc' );

             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Int_Set( unit, ior, stotel * nphase, len_name, fcn_name, wic_vol_bc )
             else
                wic_vol_bc = int( value_real )
             endif

          Case( 'wic_d_bc' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Int_Set( unit, ior, stotel * nphase, len_name, fcn_name, wic_d_bc )
             else
                wic_d_bc = int( value_real )
             endif

          Case( 'wic_u_bc' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Int_Set( unit, ior, stotel * nphase, len_name, fcn_name, wic_u_bc )
             else
                wic_u_bc = int( value_real )
             endif

          Case( 'wic_p_bc' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Int_Set( unit, ior, stotel * nphase, len_name, fcn_name, wic_p_bc )
             else
                wic_p_bc = int( value_real )
             endif

          Case( 'wic_t_bc' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Int_Set( unit, ior, stotel * nphase, len_name, fcn_name, wic_t_bc )
             else
                wic_t_bc = int( value_real )
             endif

          Case( 'wic_comp_bc' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Int_Set( unit, ior, stotel * nphase, len_name, fcn_name, wic_comp_bc )
             else
                wic_comp_bc = int( value_real )
             endif

          Case( 'uabs_option' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Int_Set( unit, ior, nphase, len_name, fcn_name, uabs_option )
             else
                uabs_option = int( value_real )
             endif

          Case( 'eos_option' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Int_Set( unit, ior, nphase, len_name, fcn_name, eos_option )
             else
                eos_option = int( value_real )
             endif

          Case( 'cp_option' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Int_Set( unit, ior, nphase, len_name, fcn_name, cp_option )
             else
                cp_option = int( value_real )
             endif

             ! Real arrays:
          Case( 'suf_vol_bc' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * cv_snloc * nphase, len_name, fcn_name, suf_vol_bc )
             else
                suf_vol_bc = value_real
             endif

          Case( 'suf_d_bc' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * cv_snloc * nphase, len_name, fcn_name, suf_d_bc )
             else
                suf_d_bc = value_real
             endif

          Case( 'suf_cpd_bc' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * cv_snloc * nphase, len_name, fcn_name, suf_cpd_bc )
             else
                suf_cpd_bc = value_real
             endif

          Case( 'suf_t_bc' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * cv_snloc * nphase, len_name, fcn_name, suf_t_bc )
             else
                suf_t_bc = value_real
             endif

          Case( 'suf_p_bc' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * p_snloc * nphase, len_name, fcn_name, suf_p_bc )
             else
                suf_p_bc = value_real
             endif

          Case( 'suf_u_bc' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * u_snloc * nphase, len_name, fcn_name, suf_u_bc )
             else
                suf_u_bc = value_real
             endif

          Case( 'suf_v_bc' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * u_snloc * nphase, len_name, fcn_name, suf_v_bc )
             else
                suf_v_bc = value_real
             endif

          Case( 'suf_w_bc' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * u_snloc * nphase, len_name, fcn_name, suf_w_bc )
             else
                suf_w_bc = value_real
             endif

          Case( 'suf_one_bc' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * cv_snloc * nphase, len_name, fcn_name, suf_one_bc )
             else
                suf_one_bc = value_real
             endif

          Case( 'suf_comp_bc' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * cv_snloc * nphase * ncomp, len_name, fcn_name, suf_comp_bc )
             else
                suf_comp_bc = value_real
             endif

          Case( 'suf_u_bc_rob1' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * u_snloc * nphase, len_name, fcn_name, suf_u_bc_rob1 )
             else
                suf_u_bc_rob1 = value_real
             endif

          Case( 'suf_u_bc_rob2' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * u_snloc * nphase, len_name, fcn_name, suf_u_bc_rob2 )
             else
                suf_u_bc_rob2 = value_real
             endif

          Case( 'suf_v_bc_rob1' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * u_snloc * nphase, len_name, fcn_name, suf_v_bc_rob1 )
             else
                suf_v_bc_rob1 = value_real
             endif

          Case( 'suf_v_bc_rob2' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * u_snloc * nphase, len_name, fcn_name, suf_v_bc_rob2 )
             else
                suf_v_bc_rob2 = value_real
             endif

          Case( 'suf_w_bc_rob1' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * u_snloc * nphase, len_name, fcn_name, suf_w_bc_rob1 )
             else
                suf_w_bc_rob1 = value_real
             endif

          Case( 'suf_w_bc_rob2' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * u_snloc * nphase, len_name, fcn_name, suf_w_bc_rob2 )
             else
                suf_w_bc_rob2 = value_real
             endif

          Case( 'suf_t_bc_rob1' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * cv_snloc * nphase, len_name, fcn_name, suf_t_bc_rob1 )
             else
                suf_t_bc_rob1 = value_real
             endif

          Case( 'suf_t_bc_rob2' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * cv_snloc * nphase, len_name, fcn_name, suf_t_bc_rob2 )
             else
                suf_t_bc_rob2 = value_real
             endif

          Case( 'suf_vol_bc_rob1' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * cv_snloc * nphase, len_name, fcn_name, suf_vol_bc_rob1 )
             else
                suf_vol_bc_rob1 = value_real
             endif

          Case( 'suf_vol_bc_rob2' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * cv_snloc * nphase, len_name, fcn_name, suf_vol_bc_rob2 )
             else
                suf_vol_bc_rob2 = value_real
             endif

          Case( 'suf_comp_bc_rob1' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * cv_snloc * nphase, len_name, fcn_name, suf_comp_bc_rob1 )
             else
                suf_comp_bc_rob1 = value_real
             endif

          Case( 'suf_comp_bc_rob2' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, stotel * cv_snloc * nphase, len_name, fcn_name, suf_comp_bc_rob2 )
             else
                suf_comp_bc_rob2 = value_real
             endif

          Case( 'opt_vel_upwind_coefs' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, nopt_vel_upwind_coefs, len_name, fcn_name, opt_vel_upwind_coefs )
             else
                opt_vel_upwind_coefs = value_real
             endif

          Case( 'x' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, x_nonods, len_name, fcn_name, x )
             else
                x = value_real
             endif

          Case( 'y' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, x_nonods, len_name, fcn_name, y )
             else
                y = value_real
             endif

          Case( 'z' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, x_nonods, len_name, fcn_name, z )
             else
                z = value_real
             endif

          Case( 'xu' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, xu_nonods, len_name, fcn_name, xu )
             else
                xu = value_real
             endif

          Case( 'yu' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, xu_nonods, len_name, fcn_name, yu )
             else
                yu = value_real
             endif

          Case( 'zu' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, xu_nonods, len_name, fcn_name, zu )
             else
                zu = value_real
             endif

          Case( 'nu' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, u_nonods * nphase, len_name, fcn_name, nu )
             else
                nu = value_real
             endif

          Case( 'nv' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, u_nonods * nphase, len_name, fcn_name, nv )
             else
                nv = value_real
             endif

          Case( 'nw' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, u_nonods * nphase, len_name, fcn_name, nw )
             else
                nw = value_real
             endif

          Case( 'ug' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, u_nonods * nphase, len_name, fcn_name, ug )
             else
                ug = value_real
             endif

          Case( 'vg' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, u_nonods * nphase, len_name, fcn_name, vg )
             else
                vg = value_real
             endif

          Case( 'wg' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, u_nonods * nphase, len_name, fcn_name, wg )
             else
                wg = value_real
             endif

          Case( 'u_source' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, u_nonods * nphase, len_name, fcn_name, u_source )
             else
                u_source = value_real
             endif

          Case( 't_source' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, cv_nonods * nphase, len_name, fcn_name, t_source )
             else
                t_source = value_real
             endif

          Case( 'v_source' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, cv_nonods * nphase, len_name, fcn_name, v_source )
             else
                v_source = value_real
             endif

          Case( 'comp_source' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, cv_nonods * nphase, len_name, fcn_name, comp_source )
             else
                comp_source = value_real
             endif

          Case( 'u' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, u_nonods * nphase, len_name, fcn_name, u )
             else
                u  = value_real
             endif

          Case( 'v' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, u_nonods * nphase, len_name, fcn_name, v )
             else
                v = value_real
             endif

          Case( 'w' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, u_nonods * nphase, len_name, fcn_name, w )
             else
                w = value_real
             endif

          Case( 'den' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, cv_nonods * nphase, len_name, fcn_name, den )
             else
                den = value_real
             endif

          Case( 'satura' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, cv_nonods * nphase, len_name, fcn_name, satura )
             else
                satura = value_real
             endif

          Case( 'comp' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, cv_nonods * nphase * ncomp, len_name, fcn_name, comp )
             else
                comp = value_real
             endif

          Case( 'volfra' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, cv_nonods * nphase , len_name, fcn_name, volfra )
             else
                volfra = value_real
             endif

          Case( 't' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, cv_nonods * nphase , len_name, fcn_name, t )
             else
                t = value_real
             endif

          Case( 'cv_one' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, cv_nonods * nphase , len_name, fcn_name, cv_one )
             else
                cv_one = value_real
             endif

          Case( 'p' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, cv_nonods, len_name, fcn_name, p )
             else
                p = value_real
             endif

          Case( 'cv_p' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, cv_nonods, len_name, fcn_name, cv_p )
             else
                cv_p = value_real
             endif

          Case( 'volfra_pore' );
             if(( ior == 13 ) .or. ( ior == 15 )) then
                call Array_Real_Set( unit, ior, totele, len_name, fcn_name, volfra_pore )
             else
                volfra_pore = value_real
             endif

          Case( 'K_Comp' );
             if(( ior == 13 ) .or. ( ior == 17 )) then
                call Array_Matrix3_Set( unit, ior, ncomp, nphase, nphase, len_name, fcn_name, K_Comp )
             else
                K_Comp = value_real
             endif

             ! Matrices and Tensors ( x , y )
          Case( 'uabs_coefs' );
             if(( ior == 13 ) .or. ( ior == 16 )) then
                call Array_Matrix2_Set( unit, ior, nphase, nuabs_coefs, len_name, fcn_name, uabs_coefs )
             else
                uabs_coefs = value_real
             endif

          Case( 'eos_coefs' );
             if(( ior == 13 ) .or. ( ior == 16 )) then
                call Array_Matrix2_Set( unit, ior, nphase, ncoef, len_name, fcn_name, eos_coefs )
             else
                eos_coefs = value_real
             endif

          Case( 'cp_coefs' );
             if(( ior == 13 ) .or. ( ior == 16 )) then
                call Array_Matrix2_Set( unit, ior, nphase, ncp_coefs, len_name, fcn_name, cp_coefs )
             else
                cp_coefs = value_real
             endif

             ! Matrices and Tensors ( x , y , z)

          Case( 'comp_diff_coef' );
             if(( ior == 13 ) .or. ( ior == 17 )) then
                call Array_Matrix3_Set( unit, ior, ncomp, ncomp_diff_coef, nphase, len_name, fcn_name, comp_diff_coef )
             else
                comp_diff_coef = value_real
             endif

          Case( 'capil_pres_coef' );
             if(( ior == 13 ) .or. ( ior == 17 )) then
                call Array_Matrix3_Set( unit, ior, ncapil_pres_coef, nphase, nphase, len_name, fcn_name, capil_pres_coef )
             else
                capil_pres_coef = value_real
             endif

          Case( 'u_abs_stab' );
             if(( ior == 13 ) .or. ( ior == 17 )) then
                call Array_Matrix3_Set( unit, ior,  mat_nonods, ndim * nphase, ndim * nphase , len_name, fcn_name, u_abs_stab )
             else
                u_abs_stab = value_real
             endif

          Case( 'u_absorb' );
             if(( ior == 13 ) .or. ( ior == 17 )) then
                call Array_Matrix3_Set( unit, ior,  mat_nonods, ndim * nphase, ndim * nphase , len_name, fcn_name, u_absorb )
             else
                u_absorb = value_real
             endif

          Case( 't_absorb' );
             if(( ior == 13 ) .or. ( ior == 17 )) then
                call Array_Matrix3_Set( unit, ior,  cv_nonods, nphase, nphase, len_name, fcn_name, t_absorb )
             else
                t_absorb  = value_real
             endif

          Case( 'v_absorb' );
             if(( ior == 13 ) .or. ( ior == 17 )) then
                call Array_Matrix3_Set( unit, ior,  cv_nonods, nphase, nphase, len_name, fcn_name, v_absorb )
             else
                v_absorb  = value_real
             endif

          Case( 'perm' );
             if(( ior == 13 ) .or. ( ior == 17 )) then
                call Array_Matrix3_Set( unit, ior, totele, ndim, ndim, len_name, fcn_name, perm )
             else
                perm = value_real
             endif

             ! Matrices and Tensors ( x , y , z , w)
          Case( 'udiffusion' );
             if(( ior == 13 ) .or. ( ior == 18 )) then
                call Array_Matrix4_Set( unit, ior,  mat_nonods, ndim, ndim, nphase, len_name, fcn_name, udiffusion )
             else
                udiffusion = value_real
             endif

          Case( 'tdiffusion' );
             if(( ior == 13 ) .or. ( ior == 18 )) then
                call Array_Matrix4_Set( unit, ior,  mat_nonods, ndim, ndim, nphase, len_name, fcn_name, tdiffusion )
             else
                tdiffusion = value_real
             endif

          Case( 'comp_diffusion' );
             if(( ior == 13 ) .or. ( ior == 18 )) then
                call Array_Matrix4_Set( unit, ior, mat_nonods, ndim, ndim, nphase, len_name, fcn_name, comp_diffusion )
             else
                comp_diffusion = value_real
             endif

          Case DEFAULT;
             write( 375, * ) 'Option not found - read_all subrt., ', ifile( 1 : k ) 
             stop 788

          end Select

       end if Conditional_IOR

    end do

    write(357,*) 'Leaving read_all'

    return
  end subroutine read_all


  subroutine readin_bc( unit_input, nphase, stotel, cv_snloc, p_snloc, u_snloc, &
       wic_vol_bc, wic_d_bc, wic_u_bc, wic_p_bc, wic_t_bc, & 
       suf_vol_bc, suf_d_bc, suf_cpd_bc, suf_t_bc, suf_p_bc, &
       suf_u_bc, suf_v_bc, suf_w_bc, suf_one_bc, &
       suf_u_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob1, suf_v_bc_rob2, &
       suf_w_bc_rob1, suf_w_bc_rob2, suf_t_bc_rob1, suf_t_bc_rob2  )
    implicit none
    integer, intent( in ) :: unit_input
    integer, intent( in ) :: nphase, stotel, cv_snloc, p_snloc, u_snloc
    integer, dimension( stotel * nphase ), intent( inout ) :: wic_vol_bc, wic_d_bc, wic_u_bc, wic_p_bc, wic_t_bc
    real, dimension( stotel * cv_snloc * nphase ), intent( inout ) :: suf_vol_bc, suf_d_bc, suf_cpd_bc, suf_t_bc
    real, dimension( stotel * p_snloc * nphase ), intent( inout ) :: suf_p_bc
    real, dimension( stotel * u_snloc * nphase ), intent( inout ) :: suf_u_bc, suf_v_bc, suf_w_bc
    real, dimension( stotel * cv_snloc * nphase ), intent( inout ) :: suf_one_bc
    real, dimension( stotel * u_snloc * nphase ), intent( inout ) :: suf_u_bc_rob1, suf_u_bc_rob2, suf_v_bc_rob1, suf_v_bc_rob2, &
         suf_w_bc_rob1, suf_w_bc_rob2
    real, dimension( stotel * cv_snloc * nphase ), intent( inout ) :: suf_t_bc_rob1, suf_t_bc_rob2

    ! Local variables
    character( len = 150) :: title

    write(357,*) 'In readin_bc'

    read( unit_input, * ) title( 1 : 150 )
    read( unit_input, * ) title( 1 : 150 )
    read( unit_input, * ) title( 1 : 150 )
    read( unit_input, * ) title( 1 : 150 )
    call input_int( unit_input, stotel * nphase, wic_vol_bc )
    call input_int( unit_input, stotel * nphase, wic_d_bc )
    call input_int( unit_input, stotel * nphase, wic_u_bc )
    call input_int( unit_input, stotel * nphase, wic_p_bc )
    call input_int( unit_input, stotel * nphase, wic_t_bc )

    call input_real( unit_input, stotel * cv_snloc * nphase, suf_vol_bc )
    call input_real( unit_input, stotel * cv_snloc * nphase, suf_d_bc )
    call input_real( unit_input, stotel * cv_snloc * nphase, suf_cpd_bc )
    call input_real( unit_input, stotel * cv_snloc * nphase, suf_t_bc )

    call input_real( unit_input, stotel * p_snloc * nphase, suf_p_bc )

    call input_real( unit_input, stotel * u_snloc * nphase, suf_u_bc )
    call input_real( unit_input, stotel * u_snloc * nphase, suf_v_bc )
    call input_real( unit_input, stotel * u_snloc * nphase, suf_w_bc )

    call input_real( unit_input, stotel * cv_snloc * nphase, suf_one_bc )

    call input_real( unit_input, stotel * u_snloc * nphase, suf_u_bc_rob1 )
    call input_real( unit_input, stotel * u_snloc * nphase, suf_u_bc_rob2 )
    call input_real( unit_input, stotel * u_snloc * nphase, suf_v_bc_rob1 )
    call input_real( unit_input, stotel * u_snloc * nphase, suf_v_bc_rob2 )
    call input_real( unit_input, stotel * u_snloc * nphase, suf_w_bc_rob1 )
    call input_real( unit_input, stotel * u_snloc * nphase, suf_w_bc_rob2 )
    call input_real( unit_input, stotel * cv_snloc * nphase, suf_t_bc_rob1 )
    call input_real( unit_input, stotel * cv_snloc * nphase, suf_t_bc_rob2 )

    write(357,*) 'Leaving readin_bc'

    return
  end subroutine readin_bc

  subroutine reading_initial_position_velocities( unit_input, nphase, &
       x_nonods, xu_nonods, u_nonods, &
       x, y, z, xu, yu, zu, nu, nv, nw, ug, vg, wg )

    implicit none
    integer, intent( in ) :: unit_input, nphase, x_nonods, xu_nonods,u_nonods
    real, dimension( x_nonods ), intent( inout ) :: x, y, z
    real, dimension( xu_nonods ), intent( inout ) :: xu, yu, zu
    real, dimension( u_nonods * nphase ), intent( inout ) ::  nu, nv, nw, ug, vg, wg

    write(357,*) 'In reading_initial_position_velocities'

    call input_real( unit_input, x_nonods, x )
    call input_real( unit_input, x_nonods, y )
    call input_real( unit_input, x_nonods, z )

    call input_real( unit_input, xu_nonods, xu )
    call input_real( unit_input, xu_nonods, yu )
    call input_real( unit_input, xu_nonods, zu )

    call input_real( unit_input, u_nonods * nphase, nu )
    call input_real( unit_input, u_nonods * nphase, nv )
    call input_real( unit_input, u_nonods * nphase, nw )

    call input_real( unit_input, u_nonods * nphase, ug )
    call input_real( unit_input, u_nonods * nphase, vg )
    call input_real( unit_input, u_nonods * nphase, wg )

    write(357,*) 'Leaving reading_initial_position_velocities'

    return
  end subroutine reading_initial_position_velocities


  subroutine input_int( unit, n, array )
    implicit none
    integer, intent( in ) :: unit, n
    integer, dimension( n ), intent( inout ) :: array

    ! Local variables
    character( len = 100 ) :: description, fcn_name
    integer :: value, i, k, npos, length
    integer, dimension( : ), allocatable :: ind

    read( unit, 101 ) description
    read( unit, * ) value
    if( value < -1000 ) then ! Pre-defined function that can be used
       read( unit, * ) fcn_name
       k = index( fcn_name, ' ' ) - 1
       call system( 'rm -f fvalues' )
       call system( './' // fcn_name( 1 : k ) // ' > fvalues ' )
       open( unit + 1, file = 'fvalues', status = 'unknown' )
       read( unit + 1, * ) length
       do i = 1, length
          read( unit + 1 , * ) array( i )
       end do
       close( unit + 1 )
    else
       array( 1 : n ) = value
       ! For specific values in the array, specify number of index (if negative, it is assumed 
       ! the array is filled with one single value as described by the previous 'read' statement.
       ! syntax is:
       ! 2 
       ! 1 6 5 8 --> two index to be modified, 1 and 6, and array in these positions  assumed 
       ! values of 5 and 8
       read( unit, * ) npos
       if( npos > n )then
          print*, 'number of index is larger the size of the array'
          stop 119
       elseif( npos < 0 ) then
          return
       else
          allocate( ind( npos ))
          read( unit, * ) (ind( i ), i = 1, npos ), ( array ( ind( i )), i = 1, npos )
       end if
    end if

101 format( a100)

    return
  end subroutine input_int

  subroutine input_real( unit, n, array )
    implicit none
    integer, intent( in ) :: unit, n
    real, dimension( n ), intent( inout ) :: array

    ! Local variables
    character( len = 100 ) :: description, fcn_name
    real :: value
    integer :: i, npos, k, length
    integer, dimension( : ), allocatable :: ind

    read( unit, 101 ) description
    read( unit, * ) value
    if( value < -1000. ) then ! pre-defined function that can be used
       read( unit, * ) fcn_name
       k = index( fcn_name, ' ' ) - 1
       call system( 'rm -f fvalues' )
       call system( './' // fcn_name( 1 : k ) // ' > fvalues ' )
       open( unit + 1, file = 'fvalues', status = 'unknown' )
       read( unit + 1, * ) length
       do i = 1, length
          read( unit + 1 , * ) array( i )
       end do
       close( unit + 1 )
    else
       array( 1 : n ) = value
       ! For specific values in the array, specify number of index (if negative, it is assumed 
       ! the array is filled with one single value as described by the previous 'read' statement.
       ! syntax is:
       ! 2 
       ! 1 6 5 8 --> two index to be modified, 1 and 6, and array in these positions  assumed 
       ! values of 5 and 8
       read( unit, * ) npos
       if( npos > n )then
          print*, 'number of index is larger the size of the array'
          stop 119
       elseif( npos < 0 ) then
          return
       else
          allocate( ind( npos ))
          read( unit, * ) (ind( i ), i = 1, npos ), ( array ( ind( i )), i = 1, npos )
       end if
    end if

101 format( a100)

    return
  end subroutine input_real

  subroutine input_matrix( unit, ndim1, ndim2, array )
    ! Here, user can either provide a single value for all components of
    ! the tensor or use a external function to design the tensor with 
    ! input array( i, j, k ) 
    implicit none
    integer, intent( in ) :: unit, ndim1, ndim2
    real, dimension( ndim1, ndim2 ), intent( inout ) :: array
    ! Local variables
    character( len = 100 ) :: description, fcn_name
    integer :: i, j,  k, length1, length2 
    real :: value

    read( unit, 101 ) description
    read( unit, * ) value
    if( value < -1000. ) then ! pre-defined function that can be used
       read( unit, * ) fcn_name
       k = index( fcn_name, ' ' ) - 1
       call system( 'rm -f fvalues' )
       call system( './' // fcn_name( 1 : k ) // ' > fvalues ' )
       open( unit + 1, file = 'fvalues', status = 'unknown' )
       read( unit + 1, * ) length1, length2 
       if( ( length1 > ndim1 ) .or. ( length2 > ndim2 ))then
          print*, 'Matrix dimension exceeded'
          stop 910
       end if
       do i = 1, length1
          do j = 1, length2
             read( unit + 1 , * ) array( i, j )
          end do
       end do
       close( unit + 1 )
    else
       array( 1 : ndim1, 1 : ndim2 ) = value
    end if

101 format( a100)

    return
  end subroutine input_matrix


  subroutine input_tensor( unit, ndim1, ndim2, ndim3, array )
    ! Here, user can either provide a single value for all components of
    ! the tensor or use a external function to design the tensor with 
    ! input array( i, j, k ) 
    implicit none
    integer, intent( in ) :: unit, ndim1, ndim2, ndim3
    real, dimension( ndim1, ndim2, ndim3 ), intent( inout ) :: array
    ! Local variables
    character( len = 100 ) :: description, fcn_name
    integer :: i, j,  k, length1, length2, length3 
    real :: value

    read( unit, 102 ) description
    read( unit, * ) value
    if( value < -1000. ) then ! pre-defined function that can be used
       read( unit, * ) fcn_name
       k = index( fcn_name, ' ' ) - 1
       call system( 'rm -f fvalues' )
       call system( './' // fcn_name( 1 : k ) // ' > fvalues ' )
       open( unit + 1, file = 'fvalues', status = 'unknown' )
       read( unit + 1, * ) length1, length2, length3 
       if( ( length1 > ndim1 ) .or. ( length2 > ndim2 ) .or. &
            ( length3 > ndim3 ))then
          print*, 'Tensor dimension exceeded'
          stop 909
       end if
       do i = 1, length1
          do j = 1, length2
             do k = 1, length3
                read( unit + 1 , * ) array( i, j, k )
             end do
          end do
       end do
       close( unit + 1 )
    else
       array( 1 : ndim1, 1 : ndim2, 1 : ndim3  ) = value
    end if

102 format( a100 )

    return
  end subroutine input_tensor



  subroutine input_tensor2( unit, ndim1, ndim2, ndim3, ndim4, array )
    ! Here, user can either provide a single value for all components of
    ! the tensor or use a external function to design the tensor with 
    ! input array( i, j, k ) 
    implicit none
    integer, intent( in ) :: unit, ndim1, ndim2, ndim3, ndim4
    real, dimension( ndim1, ndim2, ndim3, ndim4 ), intent( inout ) :: array
    ! Local variables
    character( len = 100 ) :: description, fcn_name
    integer :: i, j,  k, l, length1, length2, length3, length4 
    real :: value

    read( unit, 101 ) description
    read( unit, * ) value
    if( value < -1000. ) then ! pre-defined function that can be used
       read( unit, * ) fcn_name
       k = index( fcn_name, ' ' ) - 1
       call system( 'rm -f fvalues' )
       call system( './' // fcn_name( 1 : k ) // ' > fvalues ' )
       open( unit + 1, file = 'fvalues', status = 'unknown' )
       read( unit + 1, * ) length1, length2, length3 
       if( ( length1 > ndim1 ) .or. ( length2 > ndim2 ) .or. &
            ( length3 > ndim3 ) .or. ( length4 > ndim4 ))then
          print*, 'Tensor dimension exceeded'
          stop 910
       end if
       do i = 1, length1
          do j = 1, length2
             do k = 1, length3
                do l = 1, length4
                   read( unit + 1 , * ) array( i, j, k, l )
                end do
             end do
          end do
       end do
       close( unit + 1 )
    else
       array( 1 : ndim1, 1 : ndim2, 1 : ndim3, ndim4  ) = value
    end if

101 format( a100)

    return
  end subroutine input_tensor2


  subroutine allocating_global_nodes( ndim, totele, domain_length, &
       u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, mat_nloc, &
       cv_snloc, u_snloc, p_snloc, stotel, &
       cv_nonods, u_nonods, x_nonods, xu_nonods, &
       u_ele_type, cv_ele_type, &
       x, xu, &
       u_ndgln, xu_ndgln, cv_ndgln, x_ndgln, p_ndgln, &
       mat_ndgln, u_sndgln, cv_sndgln, p_sndgln )
    use printout

    implicit none
    integer, intent( in ) :: ndim, totele
    real, intent( in ) :: domain_length
    integer, intent( in ) :: u_nloc, xu_nloc, cv_nloc, x_nloc, p_nloc, mat_nloc, &
         cv_snloc, u_snloc, p_snloc, stotel, &
         cv_nonods, u_nonods, x_nonods, xu_nonods, &
         u_ele_type, cv_ele_type
    real, dimension( x_nonods ), intent( inout )  :: x    
    real, dimension( xu_nonods ), intent( inout ) :: xu
    integer, dimension( totele * u_nloc ), intent( inout )  :: u_ndgln
    integer, dimension( totele * xu_nloc ), intent( inout )  :: xu_ndgln
    integer, dimension( totele * cv_nloc ), intent( inout )  :: cv_ndgln
    integer, dimension( totele * cv_nloc ), intent( inout )  :: x_ndgln
    integer, dimension( totele * p_nloc ), intent( inout )  :: p_ndgln
    integer, dimension( totele * mat_nloc ), intent( inout )  :: mat_ndgln
    integer, dimension( stotel * u_snloc ), intent( inout )  :: u_sndgln
    integer, dimension( stotel * cv_snloc ), intent( inout )  :: cv_sndgln
    integer, dimension( stotel * p_snloc ), intent( inout )  :: p_sndgln

    ! Local variables:
    integer :: u_nloc2, u_nod, xu_nod, cv_nod, x_nod, mat_nod, p_nod
    real :: dx
    integer :: ele, iloc, cv_iloc
    character( len = 100 ) :: name

    write(357,*) 'In allocating_global_nodes'

    Conditional_NDIM: if( ndim == 1 ) then ! This needs to be updated for 2-3D

       dx = domain_length / real( totele )
       cv_sndgln( 1 ) = 1
       cv_sndgln( 2 ) = cv_nonods
       p_sndgln( 1 ) = 1
       p_sndgln( 2 ) = cv_nonods

       if( cv_ele_type == 2 ) then
          u_nloc2 = u_nloc / cv_nloc
          do cv_iloc = 1, cv_nloc 
             u_sndgln( cv_iloc ) = 1 + ( cv_iloc -1 ) * u_nloc2
             u_sndgln( cv_iloc + cv_nloc ) = u_nonods - u_nloc + cv_iloc * u_nloc2
          end do
       else
          u_sndgln( 1 ) = 1
          u_sndgln( 2 ) = u_nonods
       end if
       name = '####u_sndgln####'
       call mirror_array_int( 357, name,  stotel * u_snloc, u_sndgln )

       u_nod = 0
       xu_nod = 0 
       cv_nod = 0
       x_nod = 0
       mat_nod = 0
       p_nod = 0

       Loop_Elements: do ele = 1, totele

          Loop_U: do iloc = 1, u_nloc ! storing velocity nodes
             u_nod = u_nod + 1
             u_ndgln( ( ele - 1 ) * u_nloc + iloc ) = u_nod
          end do Loop_U

          Loop_XU:do iloc = 1, xu_nloc
             xu_nod = xu_nod + 1
             xu_ndgln( ( ele - 1 ) * xu_nloc + iloc ) = xu_nod
             if ( xu_nloc == 1 ) then
                xu( xu_nod ) = ( real( ele - 1 ) + 0.5 ) * dx
             else
                if( u_ele_type == 2 ) then
                   xu( xu_nod ) = real( ele - 1 ) * dx + real ( mod( iloc, cv_nloc ) - 1 ) &
                        * dx / real ( u_nloc - 1 )
                else
                   xu( xu_nod ) = real( ele - 1 ) * dx + real ( iloc - 1 ) &
                        * dx / real ( u_nloc - 1 )
                end if
             end if
          end do Loop_XU
          if( xu_nloc /= 1 ) xu_nod = xu_nod - 1

          Loop_P: do iloc = 1, p_nloc
             p_nod = p_nod + 1
             p_ndgln( ( ele - 1 ) * p_nloc + iloc ) = p_nod
          end do Loop_P
          !if( problem == 2 ) p_nod = p_nod - 1
          if( cv_nonods /= totele * cv_nloc ) p_nod = p_nod - 1

          Loop_CV: do iloc = 1, cv_nloc
             cv_nod = cv_nod + 1
             cv_ndgln( ( ele - 1 ) * cv_nloc + iloc ) = cv_nod
          end do Loop_CV
          !if( problem == 2 ) cv_nod = cv_nod - 1
          if( cv_nonods /= totele * cv_nloc ) cv_nod = cv_nod - 1

          Loop_X: do iloc = 1, x_nloc
             x_nod = x_nod + 1
             x_ndgln( ( ele - 1 ) * x_nloc + iloc ) = x_nod
             if ( x_nloc == 1 ) then
                x( x_nod ) = ( real( ele - 1 ) + 0.5 ) * dx 
             else
                x( x_nod ) = real( ele - 1 ) * dx + real( iloc - 1 ) * dx / real ( x_nloc - 1 )
             end if
          end do Loop_X
          if( x_nloc /= 1 ) x_nod = x_nod - 1

          Loop_Mat: do iloc = 1, mat_nloc
             mat_nod = mat_nod + 1
             mat_ndgln( ( ele - 1 ) * mat_nloc + iloc ) = mat_nod
          end do Loop_Mat

       end do Loop_Elements

    end if Conditional_NDIM

    write(357,*) 'Leaving allocating_global_nodes'


    return
  end subroutine allocating_global_nodes


  ! Initialising T and Told:
  subroutine initialise_scalar_fields( &
       problem, ndim, nphase, totele, domain_length, &
       x_nloc, cv_nloc, x_nonods, cv_nonods,  &
       x_ndgln, cv_ndgln, &
       x, told, t )
    implicit none
    integer, intent( in ) :: problem, ndim, nphase, totele
    real, intent( in ) :: domain_length
    integer, intent( in ) ::x_nloc, cv_nloc, x_nonods, cv_nonods
    integer, dimension( totele * cv_nloc ) :: x_ndgln, cv_ndgln
    real, dimension( x_nonods ), intent( in )  :: x
    real, dimension( cv_nonods * nphase ), intent( inout )  :: told, t
    ! Local variables
    integer :: iloc, ele, x_nod, cv_nod
    real :: dx

    write(357,*) 'In initialise_scalar_fields'

    Conditional_NDIM: if ( ndim == 1 ) then ! This may change for 2-3D

       dx = domain_length / real( totele )
       t = 0.
       Loop_Element: do ele = 1, totele
          do iloc = 1, cv_nloc
             x_nod = x_ndgln( ( ele - 1 ) * x_nloc + iloc )
             cv_nod = cv_ndgln( ( ele - 1 ) * cv_nloc + iloc )

             Select case( problem )

             case( -1, 0 ); ! CV-Adv test-case
                if( ( x( x_nod ) > 0.199 ) .and. ( x( x_nod ) < 0.401 )) t( cv_nod ) = 1.
                if( problem == 0 ) t( cv_nod ) = exp( -(( x( x_nod ) - 0.6 ) / 0.2 )**2 )

             case( 1, 2 ); ! BL problem
                if( x( x_nod ) < dx ) t( cv_nod ) = 0.5

             case DEFAULT
                write(357,*) 'It is necessary to set up -1 < problem < 1'
                stop 1789

             end Select

          end do
       end do Loop_Element
       told = t

    end if Conditional_NDIM

    write(357,*) 'Leaving initialise_scalar_fields'

    return
  end subroutine initialise_scalar_fields


  integer function Combination( n, r )
    ! This function performs the combinatorial:
    ! C(n,r) = n! / ( (n-r)! r! )
    implicit none
    integer :: n, r

    Combination = Permut( n ) / ( Permut( n - r ) * Permut( r ))

    return
  end function Combination

  integer function Permut( n )
    ! This function performs probabilistic permutation:
    ! P(n) = n!
    implicit none
    integer :: n
    integer :: i

    permut = 1
    do i = 1, n
       permut = permut * i
    end do

    return
  end function permut

end module input_var

