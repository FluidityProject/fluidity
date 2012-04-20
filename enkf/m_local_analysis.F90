module m_local_analysis
  use spud
  use write_state_module
  use timeloop_utilities
  use global_parameters, only: option_path_len, current_time, dt, FIELD_NAME_LEN
  use FLDebug
  use vtk_interfaces
  use memory_diagnostics
  use populate_state_module
  use read_triangle
  use Field_Options
  use fields
  use state_module
  use analysis_module

contains 
subroutine local_analysis(nrens,nrobs,ndim,nvar,A,state,measurement_state,D,E,S)
! Computes the EnKF analysis
!   use mod_states
!   use mod_measurement
!   use m_getD
!   use m_active_obs
   
!In variables (variables allocated within a subroutine)
!Ensemble matrix
   implicit none

   integer, intent(in) :: nrens            ! Size of ensemble
   integer, intent(in) :: nrobs            ! Total nb of obs
   integer, intent(in) :: ndim             ! dimension of model state
   integer, intent(in) :: nvar             ! no of variables 
   type(state_type), dimension(:), intent(in) :: state, measurement_state
   real, intent(inout) :: A(ndim,nrens)    ! ensemble matrix
   real, intent(in)    :: S(nrobs,nrens)   ! matrix holding HA' 
   real, intent(in)    :: D(nrobs)         ! vector holding d-HA
   real, intent(in)    :: E(nrobs,nrens)   ! matrix holding observation perturbations
!   type(states),                 intent(inout) :: A(nrens)        ! Ensemble matrix
!   type(measurement),            intent(in)    :: obs(nrobs)      ! measurements
!   real,                         intent(in)    :: radius

! Variables for Local Analysis
   type(vector_field) :: surface_position
   real                   subA(nvar,nrens)      ! local Ensemble matrix
   logical                lobs(nvar)           ! which measurements are active
   real, allocatable, dimension(:,:) :: subS,subE,subD
   integer nobs, nobs_save                      ! Number of active measurements

   integer m,i,j,nods_surface
   logical test      ! test analysis.F90 at a certain node
   real radius
   real, dimension(:), allocatable :: x,y
   integer :: quad_degree

   print *, ' Entering m_local_analysis'
   open(11,file='local_assimilation.in',status='old')
   read(11,*) radius
   print*, ' Radius = ', radius, ' meters'
   read(11,*) nobs_save
   print*, nobs_save, ' observations in neighbourhood (max)'
   close(11)
   if (2*ndim*nobs_save < nrens*(nobs_save+ndim)) then
      print*, '(local_analysis) For ',ndim*nrens/(2*ndim-nrens), &
            ' obs and less: Reps computed'
   else
      print*, '(local_analysis) if ',nobs_save,'>',ndim*nrens/(2*ndim-nrens), &
            ' obs: X5 computed'
   endif
   
   

   call get_option("/geometry/quadrature/degree", quad_degree)
   surface_position = read_triangle_serial('HorizontalMesh', quad_degree=quad_degree)

   nods_surface = node_count(surface_position)

   do i=1,nods_surface
!!need to modify
!         test = (i==nx/2) .and. (j==ny/2)
!!         test = ((i==216) .and. (j==426)) .or. ((i==nx/2) .and. (j==ny/2))
!!         if (depths(i,j) <= 0.0) cycle

        do m=1,nrens            
           subA(:,m)=getA(A(:,m),ndim,nvar,i) ! ensemble for grid point i
        enddo

         nobs=nobs_save
         call active_obs(i,state,measurement_state,nrobs,lobs,nobs,radius) 
                        ! Returns the array lobs which is
                        ! true for the nobs observations
         if (nobs < 2) cycle
         if(test) then
!!            print*, ' Testing at node (+sigma.dat and Kriging weights)'
!!            print*, 'i, nobs, count(lobs), modlon(i,j),modlat(i,j)'
!!            write(*,'(4i5,2f8.2)') i,nobs, count(lobs),modlon(i,j),modlat(i,j)
         endif
         allocate(subD(nobs,nrens))
         allocate(subE(nobs,nrens))
         allocate(subS(nobs,nrens))
         call getD(D,subD,nrobs,nrens,lobs,nobs) ! the innovations to use 
         call getD(E,subE,nrobs,nrens,lobs,nobs) ! the observation errors to use
         call getD(S,subS,nrobs,nrens,lobs,nobs) ! the HA' to use
!         call analysisESSL(subA, subD, subE, subS,  ndim, nrens, nobs)
         call analysis(subA, subD, subE, subS,  ndim, nrens, nobs, test)
!         call analysis2(subA, subD, R, subS, ndim, nrens, nobs, test)
!         call analysis4c(subA, R, subS, meanD, Rot, ndim, nrens, nobs, test)
!         call analysis5c(subA, R, subS, meanD, Rot, ndim, nrens, nobs, test)

         if (test) print*, ' test done '
         do m=1,nrens
            call putA(subA(:,m),A(:,m),ndim,nvar,i)          ! ensemble for grid point i,j
         enddo

         deallocate(subD, subE, subS)
   enddo
   print *, ' Leaving m_local_analysis'

end subroutine local_analysis

      function getA(A,ndim,nvar,nod) result(subA)
      implicit none
      integer, intent(in) :: ndim             ! dimension of model state
      integer, intent(in) :: nvar             ! no of variables 
      real :: subA(nvar)
      real, dimension(:), intent(in) :: A    ! ensemble matrix
      integer, intent(in) :: nod
      !local 
      integer i,j,k
      do k =1, nvar
         subA(:)=A((k-1)*nvar+nod)
      enddo
   end function getA

   subroutine putA(subA,A,ndim,nvar,nod)
      implicit none
      integer, intent(in) :: ndim             ! dimension of model state
      integer, intent(in) :: nvar             ! no of variables 
      real, dimension(:), intent(in) :: subA
      real, dimension(:), intent(inout) :: A    ! ensemble matrix
      integer, intent(in) :: nod
      !local 
      integer i,j,k

      do k =1,nvar
         A((k-1)*nvar+nod)=subA(k)  
      enddo
   end subroutine putA

subroutine getD(D,subD,nrobs,nrens,lobs,nobs)
! Returns the subD matrix corresponding to active measurements
   implicit none
   integer, intent(in)  :: nrobs
   integer, intent(in)  :: nrens
   integer, intent(in)  :: nobs
   real,    intent(in)  :: D(nrobs,nrens)
   logical, intent(in)  :: lobs(nrobs)
   real,    intent(out) :: subD(nobs,nrens)

   integer j,m

   j=0
   do m=1,nrobs
      if (lobs(m)) then
         j=j+1
         subD(j,:)=D(m,:)
      endif
   enddo

end subroutine getD

subroutine active_obs(i,state,measurement_state,nrobs,lobs,nobs,radius)

! Calculate observations within range 'radius' of a model p-cell midpoint.
! This is calculated for all points on grid.
! Select only the nobs nearest ones, updates nobs if there are not that many
! observations within the radius
! Next improvement : impose nobs in each angular sector around the observation

!  use m_dist
!  use mod_dimensions
!  use mod_measurement
!  use mod_angles

  implicit none
  type(state_type), dimension(:), intent(in) :: state,measurement_state
  type(vector_field) :: surface_position,measurement_position
  integer, intent(in)           :: i
  integer, intent(in)           :: nrobs
  real, intent(in)              :: radius
!  real, intent(in)              :: modlon(nx,ny),modlat(nx,ny)
!  type(measurement), intent(in) :: obs(nrobs)
  logical, intent(out)          :: lobs(nrobs)  ! keep data or not
  integer, intent(inout)        :: nobs         ! nb nearest obs within radius

  real tmpdist    ! List of qualified oservations
  real shortest(nobs)    ! List of qualified oservations
  integer indexes(nobs)           ! their respective index
  real, dimension(:), allocatable :: x_mod,y_mod  ! positions on model surface mesh
  real, dimension(:), allocatable :: x_obs,y_obs   ! positions on measurement mesh
  real lat0,lon0
  integer m,k,ntemp, nend
  real obslon(nobs),obslat(nobs)
  integer :: quad_degree
  integer :: j

!  real,allocatable,dimension(:) :: obslon,obslat
!  print * , nobs ,' obs retained (max)'
  shortest(:)=radius+1.
  indexes(:)=0
  lobs(:)=.false.
  ntemp=0


   call get_option("/geometry/quadrature/degree", quad_degree)

   !Coordinate x1,y1 on Surface Mesh
   !--------------------------------
   surface_position = read_triangle_serial('HorizontalMesh', quad_degree=quad_degree)

   allocate(x_mod(node_count(surface_position)))
   allocate(y_mod(node_count(surface_position)))
   do k = 1, node_count(surface_position)
      x_mod(k)=node_val(surface_position, k, 1)
      y_mod(k)=node_val(surface_position, k, 2)
   end do

   !Coordinate x1,y1 on Measurement Mesh
   !--------------------------------------
   measurement_position = read_triangle_serial('MeasurementMesh', quad_degree=quad_degree)
   if(nrobs.ne.node_count(measurement_position)) then
      print*, 'nrobs.nq.node_count(measurement_position)' 
      stop 1
   endif

   allocate(x_obs(node_count(measurement_position)))
   allocate(y_obs(node_count(measurement_position)))
   do k = 1, node_count(measurement_position)
      x_obs(k)=node_val(measurement_position, k, 1)
      y_obs(k)=node_val(measurement_position, k, 2)
   end do

  do m=1,node_count(measurement_position)
!ok     if ( abs(obs(m)%ipiv - i) + abs(obs(m)%jpiv - j) > 4*nobs ) cycle
!  ***** delete the following line **** may fix it later
!     if ( (obs(m)%ipiv - i)**2 + (obs(m)%jpiv - j)**2 > (radius/20000.)**2 ) cycle
!     lon0 = ang180(obs(m)%lon+0.001)         ! Add small number to avoid
!     lat0 = obs(m)%lat+0.001         ! singularity in spherdist. Hrmph.

     tmpdist=dist(x_mod(i),y_mod(i),x_obs(m),y_obs(m))
!     write(*,*) tmpdist
     if (tmpdist <= radius) then
!       if(mod(i*j,100)==0) then
!          print '(2f6.1,f9.1,a2,f9.1)', lon0, lat0, tmpdist/1000.,'<',radius/1000.
!       endif
       ntemp = ntemp+1
       do k=1,min(nobs,ntemp)
             if (tmpdist<shortest(k)) then
                nend=min(nobs,max(ntemp,2))
                shortest(k+1:nend)=shortest(k:nend-1) ! shift toplist back
                shortest(k)=tmpdist                   ! to insert new chart
                if (indexes(nobs)>0) lobs(indexes(nobs))=.false. ! Disqualify
                indexes(k+1:nend)=indexes(k:nend-1) 
                indexes(k)=m
                lobs(m)=.true.
                exit
             endif
       enddo
     endif
  enddo

  if (ntemp<nobs) nobs=ntemp
!  allocate(obslon(1:nobs)) 
!  allocate(obslat(1:nobs)) 
!  obslon=obs(indexes(1:nobs))%lon 
!  obslat=obs(indexes(1:nobs))%lat 
!  if(mod(i*j,100)==0) then
!  if(mod(i*j,10500)==0) then
   if( nobs< 2) then 
!!     open(7,file='target.dat')
!!     write(7,'(a,2f7.2)') 'Target point Lon-Lat: ', modlon(i,j),modlat(i,j)
!!     write(7,*) nobs, ' active obs, distances (km), lon, lat :'
!!     write(7, '(10f6.1)')  shortest(1:nobs)/1000
     !write(7, '(10f6.1)')  obs(indexes(1:nobs))%lon
     !write(7, '(10f6.1)')  obs(indexes(1:nobs))%lat
!     write(7, '(10f6.1)')  obslon
!     write(7, '(10f6.1)')  obslat
   endif
end subroutine active_obs

real function dist(x1,y1,x2,y2)
! --- -----------------------------------------
! --- Computes the distance between geo. pos.
! --- lon1,lat1 and lon2,lat2. 
! --- INPUT is in degrees.
! --- -----------------------------------------

   implicit none

   real, intent(in) :: x1,y1
   real, intent(in) :: x2,y2

!   real, parameter :: invradian=0.017453292
   real, parameter :: rearth=6371001.0     ! Radius of earth

   real  dx,dy,dz,dr                       ! Cartesian distances

   dist=SQRT( (x1-x2)**2+(y1-y2)**2)       

end function dist

end module m_local_analysis
