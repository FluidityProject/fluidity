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

module tidal_diagnostics

  use fldebug
  use global_parameters, only : timestep, OPTION_PATH_LEN, current_time
  use spud
  use futils
  use vector_tools, only: solve
  use fields
  use state_module
  use field_options
  use diagnostic_source_fields
  use initialise_fields_module
  use state_fields_module
  use tidal_module
  use write_state_module, only: do_write_state

  implicit none

  private

  public :: calculate_free_surface_history, calculate_tidal_harmonics

  ! Module level variables - 'cos we have both the free surface history and
  ! harmonic analyses to get options from
  integer, save :: nLevels_
  integer, save :: when_to_calculate_
  logical, save :: calc_diag_at_

type harmonic_field
   type(scalar_field) , pointer :: s_field
   character(len=OPTION_PATH_LEN) :: name ! name of scalar field
   integer :: sigmaIndex  ! 0 == constant consituent C0, -1 == Residual
   character(len=OPTION_PATH_LEN) :: target ! 'Amplitude' or 'Phase'
end type harmonic_field


contains

function get_number_of_harmonic_fields(state) result(N)

    type(state_type) :: state

    integer :: N
    character(len = OPTION_PATH_LEN) :: lalgorithm, base_path
    type(scalar_field), pointer :: iter_field
    integer :: ii

    N = 0
    
    do ii=1,scalar_field_count(state)
       iter_field => extract_scalar_field(state,ii)
       if (trim(iter_field%option_path)=='')then
          cycle
       end if
       base_path = trim(complete_field_path(iter_field%option_path)) // '/algorithm/'
       if (have_option(trim(base_path) // 'name')) then
           call get_option(trim(base_path) // 'name', lalgorithm, default = "Internal")
           if (trim(lalgorithm)=='tidal_harmonics') then
              N = N + 1
           end if
       end if
    end do

end function

subroutine calculate_free_surface_history(state, s_field)

    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    type(scalar_field), pointer :: hist_fs_field
    type(scalar_field), pointer :: fs_field
    character(len=OPTION_PATH_LEN) :: base_path
    integer :: stride, new_snapshot_index, levels, stat, timestep_counter
    real :: spin_up_time, current_time, timestep
    real, dimension(:), allocatable :: saved_snapshots_times

    ewrite(3,*) 'in free_surface_history_diagnostics'

    fs_field => extract_scalar_field(state,"FreeSurface",stat)
    call halo_update(fs_field)

    if(stat /= 0) then
      FLExit('I do not have a FreeSurface field so can not calculate diagnostics on it. Please switch on the FreeSurface diagnostic.')
      return
    end if

! get history options
   base_path=trim(complete_field_path(s_field%option_path)) // "/algorithm/"
   
   ! levels: the number of levels which will be saved. Too old levels will be overwritten by new ones.
   call get_option(trim(base_path) // "levels", levels, default=50)
   levels=max(levels,0)
   ! Set it for checkpointing, etc
   call set_option(trim(base_path) // "levels", levels, stat)
   nLevels_ = levels
   assert(any(stat == (/SPUD_NO_ERROR, SPUD_NEW_KEY_WARNING/)))

   ! The internal timestep counter of calculate_free_surface_history.
   if (have_option(trim(base_path) // "timestep_counter")) then
       call get_option(trim(base_path) // "timestep_counter", timestep_counter)
       timestep_counter=timestep_counter+1
   else
        timestep_counter=0
   end if
   ! Lets save the current timestep_counter in the option tree.
   call set_option(trim(base_path) // "timestep_counter", timestep_counter,stat)
   assert(any(stat == (/SPUD_NO_ERROR, SPUD_NEW_KEY_WARNING/)))

   ! stride: Defines how many timesteps shall be skipped between two history snapshots.
   call get_option(trim(base_path) // "stride", stride, default=50)
   ! Set it for checkpointing in the future
   call set_option(trim(base_path) // "stride", stride, stat)
   assert(any(stat == (/SPUD_NO_ERROR, SPUD_NEW_KEY_WARNING/)))

   ! When to calculate
   if (have_option(trim(base_path)//"/calculation_period"))   then
      call get_option(trim(base_path)//"/calculation_period",when_to_calculate_)
      calc_diag_at_ = .true.
   else
      calc_diag_at_ = .false.
   end if


   call get_option(trim(base_path)//"spin_up_time", spin_up_time, default=0.0)
   call get_option("/timestepping/current_time", current_time)
   call get_option("/timestepping/timestep", timestep)
   ! Spin up time is measured from the start of the simulation, not the start
   ! time (hence no start time in here).
   ! Note this is after the options check above as we add options to the tree if
   ! they aren't there and they are needed in tidal_harmonics diagnostics
   if(current_time+timestep<spin_up_time) then
       ewrite(4,*) "Still spinning up."
       return
   endif
   
   ! check if we want to save the current timestep at all
   if (mod(timestep_counter,stride)/=0) then
       ewrite(4,*) "Ignoring this timestep"
       return
   end if

   allocate(saved_snapshots_times(levels))
   if (have_option(trim(base_path) // "saved_snapshots_times")) then
       call get_option(trim(base_path) // "saved_snapshots_times", saved_snapshots_times)
   end if

   ! get index where we want to save the new snapshot.
   new_snapshot_index=mod(timestep_counter/stride,levels)+1 

  ! Save current free surface field in the history.
  ! Note: Since diagnositcs are executed after the solving step, we actually save the fields at current_time+timestep
  saved_snapshots_times(new_snapshot_index)=current_time+timestep
  call set_option(trim(base_path) // "saved_snapshots_times", saved_snapshots_times, stat)
  assert(any(stat == (/SPUD_NO_ERROR, SPUD_NEW_KEY_WARNING/)))
  ewrite(4,*) 'Filling history level: ', min(timestep_counter/stride+1,levels), '/', levels
 
  ! lets copy a snapshot of freesurface to s_field(new_snapshot_index)
  hist_fs_field => extract_scalar_field(state,'harmonic'//int2str(new_snapshot_index))
  call set(hist_fs_field,fs_field)
  deallocate(saved_snapshots_times)
end subroutine calculate_free_surface_history


subroutine calculate_tidal_harmonics(state, s_field)
   type(state_type), intent(in) :: state
   type(scalar_field), intent(in) :: s_field
   type(harmonic_field), dimension(:), allocatable :: harmonic_fields
   real, dimension(:), allocatable :: sigma
   integer, save :: last_update=-1, nohfs=-1, M=-1
   logical :: ignoretimestep
   real, dimension(:), allocatable :: saved_snapshots_times
   integer :: i, current_snapshot_index
   integer :: when_to_calculate
   

   ! Check dump period - if we're about to dump output, calculate, regardless of
   ! other options
   if (.not. do_write_state(current_time, timestep+1)) then 
       ! Note: diagnostics are done at the end of the timestemp, dumps at the
       ! begining. Hence the +1 on the timestep number above - we're
       ! anticipating a dump at the start of the next timestep
       ! Now check if the user wants a timestep
       if (calc_diag_at_) then
          if (.not. mod(timestep,when_to_calculate_+1) == 0) then
             return ! it's not time to calculate
          end if
       else
           return
       end if
   end if

   ! Only if Harmonics weren't already calculated in this timestep
   if (last_update/=timestep) then 
      ewrite(3,*) "In tidal_harmonics"
      allocate(harmonic_fields(get_number_of_harmonic_fields(state)))
      allocate(sigma(nLevels_))

      last_update=timestep
      call getFreeSurfaceHistoryData(state, ignoretimestep, saved_snapshots_times, current_snapshot_index)
      if (.not. ignoretimestep) then
         ! Initialize the harmonic fields and frequencies if not already done.
         call getHarmonicFields(state, harmonic_fields, nohfs, sigma, M)
         ewrite(4,*) 'Frequencies to analyse:'
         do i=1,M
            ewrite(4,*) sigma(i)
         end do
         ! Calculate harmonics and update (all!) harmonic fields
         call update_harmonic_fields(state, saved_snapshots_times, size(saved_snapshots_times), current_snapshot_index, sigma, M, harmonic_fields, nohfs)
      end if
      deallocate(harmonic_fields)
      deallocate(sigma)
   end if

   if (allocated(saved_snapshots_times)) then
     deallocate(saved_snapshots_times)
   end if


end subroutine calculate_tidal_harmonics

subroutine getFreeSurfaceHistoryData(state, ignoretimestep, saved_snapshots_times, current_snapshot_index)
    type(state_type), intent(in) :: state
    real, dimension(:), allocatable, intent(out) :: saved_snapshots_times
    logical, intent(out) :: ignoretimestep
    integer, intent(out) :: current_snapshot_index
    character(len = OPTION_PATH_LEN) :: free_surface_history_path
    integer :: timestep_counter, stride, levels, stat
    type(scalar_field), pointer :: fshistory_sfield

    ignoretimestep=.false.
    ! Find free_surface_history diagonstic field
    fshistory_sfield => extract_scalar_field(state, 'FreeSurfaceHistory')
    free_surface_history_path = trim(complete_field_path(fshistory_sfield%option_path)) // '/algorithm/'

   ! get information from free_surface_history diagnostic and check if the harmonic analysis needs to be calculated at this timestep
   ! These options should be available here, because we ran calculate_free_surface_history() as a dependency before
   call get_option(trim(free_surface_history_path) // "timestep_counter", timestep_counter)
   call get_option(trim(free_surface_history_path) // "stride", stride)
   call get_option(trim(free_surface_history_path) // "levels", levels)

   if( mod(timestep_counter,stride)/=0) then
        ewrite(4,*) 'Do nothing in this timestep.'
        ignoretimestep=.true.
        return
   end if
   if(timestep_counter/stride+1 .lt. levels) then
        ewrite(4,*) 'Do nothing until levels are filled up.'
        ignoretimestep=.true.
        return
   end if

   allocate(saved_snapshots_times(levels))
   if (have_option(trim(free_surface_history_path) // "saved_snapshots_times")) then
       call get_option(trim(free_surface_history_path) // "saved_snapshots_times", saved_snapshots_times)
   else
       call set_option(trim(free_surface_history_path) //"saved_snapshots_times", saved_snapshots_times,stat)
   end if
   current_snapshot_index=mod(timestep_counter/stride,levels)+1 
end subroutine getFreeSurfaceHistoryData


subroutine getHarmonicFields(state, harmonic_fields, nohfs, sigma, M)
   type(state_type), intent(in) :: state
   type(harmonic_field), dimension(:), intent(inout) :: harmonic_fields
   real, dimension(:), intent(inout) :: sigma
   integer, intent(inout) :: nohfs, M

   character(len = OPTION_PATH_LEN) :: lalgorithm, base_path, constituent_name, target
   integer :: i, ii
   real :: freq
   type(scalar_field), pointer :: iter_field

    nohfs=0  ! number of harmonic_fields
    M=0 ! number of sigmas
    ! Get desired constituents from the optione tree
    s_field_loop: do ii=1,scalar_field_count(state)
       iter_field => extract_scalar_field(state,ii)
       if (trim(iter_field%option_path)=='')then
          cycle
       end if
       base_path = trim(complete_field_path(iter_field%option_path)) // '/algorithm/'
       if (have_option(trim(base_path) // 'name')) then
           call get_option(trim(base_path) // 'name', lalgorithm, default = "Internal")
           if (trim(lalgorithm)=='tidal_harmonics') then
               nohfs=nohfs+1
               if (nohfs>size(harmonic_fields)) then
                  FLAbort('We found more tidal_harmonic fields than the space we allocated. Please report this as a bug')
               end if
               call get_option(trim(base_path) // 'constituent/name', constituent_name)
               harmonic_fields(nohfs)%s_field=>iter_field
               harmonic_fields(nohfs)%name=constituent_name

               if (trim(constituent_name)=='C0') then
                  harmonic_fields(nohfs)%sigmaIndex=0
                  harmonic_fields(nohfs)%target=''

               elseif (trim(constituent_name)=='Residual') then
                  harmonic_fields(nohfs)%sigmaIndex=-1
                  harmonic_fields(nohfs)%target=''
               ! Includes prescribed and custom frequencies:
               else
                  call get_option(trim(base_path) // 'target', target)
                  harmonic_fields(nohfs)%target=target
                  ! Check if we have this constituent frequency already in sigma and if not, save it
                  
                  call get_option(trim(base_path) // '/constituent', freq)
                  do i=1,M
                    if (abs(sigma(i)-freq) <= 1.0E-14) then
                      harmonic_fields(nohfs)%sigmaIndex=i
                      cycle s_field_loop
                    end if
                  end do
                  ! Frequency was not found, so lets add it to sigma
                  M=M+1
                  if (M>size(sigma)) then
                    FLExit('More frequencies than there are saved levels. Increase the number of levels (at least 2 times the number of frequencies)')
                  end if
                  sigma(M) = freq
                  harmonic_fields(nohfs)%sigmaIndex=M
               end if
           end if
       end if
    end do s_field_loop
    ewrite(4,*) 'Found ',  nohfs, ' constituents to analyse.' 
    ewrite(4,*) 'Found ', M, ' frequencies to analyse.'

   if ( M .le. 0 ) then
      FLExit("Internal error in calculate_tidal_harmonics(). No harmonic constituents were found in option tree.")
   end if
! check which constituent we want to calculate
!   call update_harmonic_analysis(state, levels, current_snapshot_index, saved_snapshots_times, sigma, M, s_field, myconstituent_id, target)
end subroutine getHarmonicFields

subroutine  update_harmonic_fields(state, snapshots_times, N, current_snapshot_index, sigma, M, harmonic_fields, nohfs)
    type(state_type), intent(in) :: state
    integer, intent(in) :: nohfs, M, N, current_snapshot_index ! M = size(sigma), N = size(snapshots_times)
    type(harmonic_field), dimension(:), intent(in) :: harmonic_fields
    real, dimension(:), intent(in) :: sigma, snapshots_times
    real, dimension(:,:), allocatable :: harmonic_A  ! for solving Ax=b system
    real, dimension(:), allocatable :: harmonic_x, harmonic_b, harmonic_time_series_vals_at_node
    integer :: i, ii, stat, node, nonodes
    logical :: forceC0toZero
    real :: residual
    type(scalar_field), pointer :: harmonic_current
!
    allocate(harmonic_A(2*M+1,2*M+1))
    allocate(harmonic_x(2*M+1))
    allocate(harmonic_b(2*M+1))
    allocate(harmonic_time_series_vals_at_node(N))

   ! Set this to .True. if you want to force C0 to 0
   forceC0toZero = .False.

   ! Loop over all nodes of the mesh
   harmonic_current => extract_scalar_field(state, 'harmonic1', stat)
   nonodes=node_count(harmonic_current)
   do node = 1,nonodes
      !Extract the free surface elevations at the current node
      do i = 1,N
       harmonic_current => extract_scalar_field(state,'harmonic'//int2str(i),stat)
       harmonic_time_series_vals_at_node(i) = node_val(harmonic_current,node)
      end do

      !Form and invert the least squares system (Unsorted version)
      call harmonic_analysis_at_single_node(N,snapshots_times,harmonic_time_series_vals_at_node,M,sigma,&
                                                  harmonic_A,harmonic_x,harmonic_b, forceC0toZero)
      ! Calculate residual
      harmonic_current => extract_scalar_field(state,'harmonic'//int2str(current_snapshot_index),stat)
      residual=node_val(harmonic_current,node)
      residual=residual-harmonic_x(1)
      do ii=1,M
           residual=residual-harmonic_x(1+ii)*cos(2*3.141592654*sigma(ii)*snapshots_times(current_snapshot_index))
           residual=residual-harmonic_x(1+M+ii)*sin(2*3.141592654*sigma(ii)*snapshots_times(current_snapshot_index))
      end do

      call save_harmonic_x_in_fields(harmonic_x, residual, M, harmonic_fields, nohfs, node)

    end do ! node loop
    deallocate(harmonic_A)
    deallocate(harmonic_x)
    deallocate(harmonic_b)
    deallocate(harmonic_time_series_vals_at_node)

 end subroutine update_harmonic_fields


subroutine save_harmonic_x_in_fields(harmonic_x, residual, M, harmonic_fields, nohfs, node)
    real, dimension(:), intent(in) :: harmonic_x
    real, intent(in) :: residual
    type(harmonic_field), dimension(:), intent(in) :: harmonic_fields
    integer, intent(in) :: nohfs, M, node
    integer :: i, MM
    real :: result
    type(scalar_field), pointer :: harmonic_current
    ! Loop over harmonic diagnostic fields
    do i=1,nohfs
      harmonic_current=>harmonic_fields(i)%s_field
      MM=harmonic_fields(i)%sigmaIndex

      ! Check if we want the C0 constituent
      if (MM==0) then
        call set(harmonic_current, node, harmonic_x(1))

      ! Check if we want the residual
      elseif (MM==-1) then
        call set(harmonic_current, node, residual)
      end if

      !stick the amplitude and phase into something that will be output
      if (harmonic_fields(i)%target=='Amplitude') then
         call set( harmonic_current, node, sqrt( harmonic_x(MM+1)**2 + harmonic_x(MM+1+M)**2 ) )
      elseif (harmonic_fields(i)%target=='Phase') then
         result = atan2(harmonic_x(MM+1+M),harmonic_x(MM+1))
         !*180.0/pi
         !if (phase < 0.0) phase = phase + 360.0
         call set( harmonic_current, node, result )
      end if
    end do
end subroutine save_harmonic_x_in_fields
!
!!!!!!!
 subroutine harmonic_analysis_at_single_node(N,harmonic_times_reordered,harmonic_time_series_vals_at_node,M,sigma,&
                                                harmonic_A,harmonic_x,harmonic_b, forceC0toZero)
    real, intent(in) :: harmonic_times_reordered(:),harmonic_time_series_vals_at_node(:),sigma(:)
    real, intent(inout) :: harmonic_A(:,:),harmonic_x(:),harmonic_b(:)
    integer, intent(in) :: M,N
    logical, optional, intent(in) :: forceC0toZero
    real :: C_k, S_k, CC_jk, SS_jk, SC_jk, CS_kj
    real :: pi = 3.141592654
    integer :: i, j, k, stat

! For the least squares system
       harmonic_A(1,1) = N

! Need the C_k and S_k for first row/column
       j = 1
       do k = 1,M
          C_k = 0.0
          S_k = 0.0
          do i = 1,N
             C_k = C_k + cos(2.*pi*sigma(k)*harmonic_times_reordered(i))
             S_k = S_k + sin(2.*pi*sigma(k)*harmonic_times_reordered(i))
          end do
          harmonic_A(1,k+1)   = C_k
          harmonic_A(1,k+1+M) = S_k
          harmonic_A(k+1,1)   = C_k
          harmonic_A(k+1+M,1) = S_k
       end do
! rest of the rows and columns of the matrix
       do j = 1,M
         do k = 1,M
          CC_jk = 0.0
          SS_jk = 0.0
          SC_jk = 0.0
          do i = 1,N
             CC_jk = CC_jk + cos(2.*pi*sigma(k)*harmonic_times_reordered(i))*cos(2.*pi*sigma(j)*harmonic_times_reordered(i))
             SS_jk = SS_jk + sin(2.*pi*sigma(k)*harmonic_times_reordered(i))*sin(2.*pi*sigma(j)*harmonic_times_reordered(i))
             SC_jk = SC_jk + cos(2.*pi*sigma(k)*harmonic_times_reordered(i))*sin(2.*pi*sigma(j)*harmonic_times_reordered(i))
          end do
          CS_kj = SC_jk
          harmonic_A(j+1,k+1)     = CC_jk ! top left quadrant
          harmonic_A(j+1+M,k+1+M) = SS_jk ! bottom right quadrant
          harmonic_A(j+1+M,k+1)   = SC_jk ! bottom left quadrant
          harmonic_A(k+1,j+1+M)   = SC_jk ! top right quadrant  (swap order we fill up matrix as CS_kj.eq.SC_jk, CS_kj.ne.SC_kj)
         end do
       end do
       if(present_and_true(forceC0toZero)) then
           harmonic_A(1,2:2*M+1)=0.0
           harmonic_A(2:2*M+1,1)=0.0
       end if
! now the rhs vector
       harmonic_b(1:2*M+1) = 0.0
       if(.not. present_and_true(forceC0toZero)) then
         do i = 1,N
            harmonic_b(1) = harmonic_b(1) + harmonic_time_series_vals_at_node(i)
         end do
       end if
       do j = 1,M
          do i = 1,N
             harmonic_b(j+1)   = harmonic_b(j+1)   + harmonic_time_series_vals_at_node(i)*cos(2.*pi*sigma(j)*harmonic_times_reordered(i))
             harmonic_b(j+1+M) = harmonic_b(j+1+M) + harmonic_time_series_vals_at_node(i)*sin(2.*pi*sigma(j)*harmonic_times_reordered(i))
          end do
       end do
! solve the system
       call solve(harmonic_A, harmonic_b, stat) ! Solve Ax=b, note that b will be overwritten
       harmonic_x = harmonic_b

  end subroutine harmonic_analysis_at_single_node

end module tidal_diagnostics

