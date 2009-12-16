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
!    C.Pain@Imperial.ac.uk
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

! #define DEBUG_OCEAN_SURFACE_FORCING 1

module OceanSurfaceForcing
  use climatology
  use elements
  use vertical_extrapolation_module
  use state_module
  use fields
  use fetools
  use transform_elements
  use global_parameters, only : OPTION_PATH_LEN, pi
  use boundary_conditions
  use spud
  use allsorts
  use coordinates
  use position_in_matrix
  
  implicit none

  private
  public::initialise, surface_stress, &
       SurfaceHeatFlux, SurfaceSalinityFlux, &
       SurfaceTemperatureRelaxation, SurfaceSalinityRelaxation, &
       RelaxSalinityToClimatology, RelaxTemperatureToClimatology, &
       wind_forcing

  interface initialise
    module procedure initialise_from_state, initialise_from_raw_data
  end interface initialise

  integer, save::NSElements
  integer, save, pointer::SENList(:)=>null(), SurfaceNodes(:)=>null()

  !density_air 1.3 kg m^-1 from Stewart  
  real, parameter::density_air=1.3         ! kg m^-1
  real, parameter::density_seawater=1023.0 ! kg m^-1
  real, parameter::Cp=3985                 ! Units J kg^-1 K^-1      
  
contains

  subroutine initialise_from_state(state)
    type(state_type), intent(in) :: state

    integer :: sloc, snelements
    integer, dimension(:), allocatable :: sndgln
    integer, dimension(:), pointer :: boundary_ids => null()
    type(mesh_type), pointer :: mesh => null()

    mesh => extract_mesh(state, "CoordinateMesh")

    snelements = surface_element_count(mesh)
    if(snelements > 0) then
      sloc = face_loc(mesh, 1)
      assert(associated(mesh%faces))
      assert(associated(mesh%faces%boundary_ids))
      boundary_ids => mesh%faces%boundary_ids
    else
      sloc = 0
      allocate(boundary_ids(0))
    end if
    allocate(sndgln(snelements * sloc))
    call getsndgln(mesh, sndgln)

    call initialise(sndgln, boundary_ids)

    deallocate(sndgln)
    if(snelements == 0) then
      deallocate(boundary_ids)
    end if

  end subroutine initialise_from_state

  subroutine initialise_from_raw_data(sndgln, boundary_ids)
    integer, intent(in)::sndgln(:), boundary_ids(:)
   
    logical, allocatable::on_top(:)
    integer i, j, loc, nid
    integer::snloc=0

    ewrite(3, *) "OceanSurfaceForcing.F90: subroutine initialize(...)"

    ! Find out how many top surface elements we've got
    NSElements = 0
    do i=1, size(boundary_ids)
       if(boundary_ids(i)==TOP_BOUNDARY_ID) then
          NSElements = NSElements + 1
       end if
    end do
    
    ewrite(3, *) "Number of free-surface elements: ", NSElements, " of ", size(boundary_ids)

    ! Number of nodes per element
    if(size(boundary_ids)>0) then
       snloc = size(sndgln)/size(boundary_ids)
    end if

    ! Local compressed copy of free-surface surface elements
    if(associated(SENList)) then
       deallocate(SENList)
    end if
    allocate(SENList(snloc*NSElements))
    allocate(on_top(maxval(SNDGLN)))
    on_top = .false.
    loc = 0
    do i=1, size(boundary_ids)
       if(boundary_ids(i)==TOP_BOUNDARY_ID) then
          do j=1, snloc
             nid = SNDGLN((i-1)*snloc+j)
             SENList(loc*snloc+j) = nid
             on_top(nid) = .true.
          end do
          loc = loc + 1
       end if
    end do
    assert(loc==NSElements)

    ! Get the list of nodes on the surface
    loc=0
    do i=1, size(on_top)
       if(on_top(i)) then
          loc = loc + 1
       end if
    end do
    if(associated(SurfaceNodes)) then
       deallocate(SurfaceNodes)
    end if
    allocate(SurfaceNodes(loc))
    loc=1
    do i=1, size(on_top)
       if(on_top(i)) then
          SurfaceNodes(loc) = i
          loc = loc + 1
       end if
    end do
    
  end subroutine initialise_from_raw_data
  
  subroutine get_sml(snloc, sngi, sn, snlx, snly, sweigh, x, y, z, sml)
    integer, intent(in)::snloc, sngi
    real, intent(in)::sn(snloc, sngi), snlx(snloc, sngi), snly(snloc, sngi)
    REAL, intent(in)::sweigh(sngi) 
    real, intent(in), dimension(:)::x, y, z
    real, intent(out)::sml(:)
    
    integer i, j, nid, sgi, XNONOD
    real NN, SDETWE(SNGI), SAREA
    
    XNONOD=size(x)
    sml = 0.0
    
    do i=1, NSElements
       CALL SDETNN(SENList((i-1)*snloc+1:i*snloc), &
            XNONOD,SNLOC,SNGI, &
            X,Y,Z, &
            SN, SNLX, SNLY, SWEIGH, SDETWE, SAREA, .true., .false.)
       do j=1,snloc
          nid = SENList((i-1)*snloc+j)
          ! Have a surface integral on outlet boundary...
          NN=0.0
          do sgi=1, sngi
             NN=NN+SDETWE(sgi)*SN(j, sgi)
          end do
          sml(nid) = sml(nid) + NN
       end do
    end do
  end subroutine get_sml

  ! Expects the speed to be given in meters-per-second
  function get_sea_surface_drag(Speed) result(c_D)
    implicit none
    real, intent(in)::Speed
    real c_D
    ! Ref: Numerical modelling of ocean dynamics, Z. Kowalik,
    ! T.S. Murty, ISBN 9-8102-1334-4, page 27 equation 1.129
    ! Still requires improvement
    c_D = 0.001*(0.75+0.067*Speed)
    
    return
  end function get_sea_surface_drag

  ! This subroutine applies ERA-40 wind stress data through surface
  ! integrals.
  subroutine surface_stress(U,V,W,NU,NV,NW,X,Y,Z, &
       snloc, sngi, sn, snlx, snly, sweigh, &
       theta,dt,centrm,ML,NCOLM, &
       BIGM, VECX, VECY, VECZ)
    real, intent(in), dimension(:)::U, V, W, NU, NV, NW, X, Y, Z
    integer, intent(in)::snloc, sngi
    real, intent(in)::sn(snloc, sngi), snlx(snloc, sngi), snly(snloc, sngi)
    real, intent(in)::sweigh(sngi) 
    real, intent(in)::theta,dt
    integer, intent(in)::centrm(:)
    real, intent(inout)::ML(:)
    integer, intent(in)::ncolm
    real, intent(inout), dimension(:)::BIGM,VECX,VECY,VECZ
    
    logical BLKSYM
    real tmpa
    integer i, nid
    real longitude, latitude, p10uv(2), Uair, Vair, Wair, Wind
    real CD
    integer IBL11, IBL12, IBL13, IBL21, IBL22, IBL23, IBL31, IBL32, IBL33
    
    integer xnonod
    real, allocatable, dimension(:)::sml
    
#ifdef DEBUG_OCEAN_SURFACE_FORCING
    integer, allocatable, dimension(:)::elementTypes, elementSizes
    real, allocatable, dimension(:)::windx, windy, windz,p10u,p10v
    integer, save::NCalls=0
    character(len=1024) str
#endif
    
    ewrite(3, *) "subroutine surface_stress()"
    
    xnonod = size(x)
    allocate(sml(xnonod))
    call get_sml(snloc, sngi, sn, snlx, snly, sweigh, x, y, z, sml)
    
#ifndef DOUBLEP
    FLAbort("ERA-40 interface does not support single-precision.")
#endif
    
    ewrite(3, *) "Read ECMWF ERA40 data"
    call fluxes_clearfields()
    call fluxes_addfieldofinterest("p10u", 4)
    call fluxes_addfieldofinterest("p10v", 4)
    
    ! Establish of the offsets to the block matricies
    IF(size(bigm).EQ.9*NCOLM) THEN
       BLKSYM=.FALSE.
    ELSE
       BLKSYM=.TRUE.
    ENDIF
    IF(BLKSYM) THEN
       IBL11=0
       IBL12=1*NCOLM
       IBL13=3*NCOLM
       IBL21=1*NCOLM
       IBL22=2*NCOLM
       IBL23=4*NCOLM
       IBL31=3*NCOLM
       IBL32=4*NCOLM
       IBL33=5*NCOLM
    ELSE
       IBL11=0
       IBL12=1*NCOLM
       IBL13=2*NCOLM
       IBL21=3*NCOLM
       IBL22=4*NCOLM
       IBL23=5*NCOLM
       IBL31=6*NCOLM
       IBL32=7*NCOLM
       IBL33=8*NCOLM
    ENDIF

#ifdef DEBUG_OCEAN_SURFACE_FORCING
    allocate(windx(size(x)), windy(size(x)), windz(size(x)), &
         p10u(size(x)), p10v(size(x)))
    windx = 0.0
    windy = 0.0
    windz = 0.0
    p10u = 0.0
    p10v = 0.0
#endif
    
    do i=1, size(SurfaceNodes)
       nid = SurfaceNodes(i)
       call LongitudeLatitude( (/ x(nid), y(nid), z(nid) /), longitude, latitude)
       call fluxes_getscalars(longitude, latitude, p10uv)
       
       call ll2r3_rotate(longitude, latitude, p10uv(1), p10uv(2), &
            Uair, Vair, Wair)
       
#ifdef DEBUG_OCEAN_SURFACE_FORCING
       ewrite(3, *) "wind ", p10uv(1), p10uv(2), Uair, Vair, Wair
       windx(nid) = Uair
       windy(nid) = Vair
       windz(nid) = Wair
       p10u(nid) = p10uv(1)
       p10v(nid) = p10uv(2)
#endif
       
       Wind = sqrt((Uair-NU(nid))**2 + (Vair-NV(nid))**2 + (Wair-NW(nid))**2)
       
       CD = get_sea_surface_drag(Wind)

       tmpa = sml(nid)*CD*density_air*Wind/density_seawater

       vecx(nid) = vecx(nid) + tmpa*(Uair-U(nid))
       vecy(nid) = vecy(nid) + tmpa*(Vair-v(nid))
       vecz(nid) = vecz(nid) + tmpa*(Wair-W(nid))

       if(ncolm.ne.0) then
          BIGM(centrm(nid) + IBL11) = BIGM(centrm(nid) + IBL11) + dt*theta*tmpa
          BIGM(centrm(nid) + IBL22) = BIGM(centrm(nid) + IBL22) + dt*theta*tmpa
          BIGM(centrm(nid) + IBL33) = BIGM(centrm(nid) + IBL33) + dt*theta*tmpa
       end if

       ML(nid) = ML(nid) + dt*theta*tmpa
    end do

#ifdef DEBUG_OCEAN_SURFACE_FORCING
    str = ' '
    write(str, '("global_stress",I6.6,".vtu")') NCalls
    NCalls = NCalls + 1

    call VTKOPEN(str, len_trim(str), "wind stress", 11)
    allocate(elementTypes(NSElements), elementSizes(NSElements))
    elementTypes = 5
    elementSizes = 3
    call VTKWRITEMESHD(size(x), NSElements, x, y, z, SENList, elementTypes, elementSizes) 
    call VTKWRITEDVN(windx, windy, windz, "wind", 4)
    call VTKWRITEDVN(vecx, vecy, vecz, "left hand side", 14)
    call VTKWRITEDSN(p10u, "p10u", 4)
    call VTKWRITEDSN(p10v, "p10v", 4)
    call VTKCLOSE()
    deallocate(elementTypes, elementSizes, windx, windy, windz)
#endif

    deallocate(sml)

    RETURN
  end subroutine surface_stress

  ! Calculates the salinity flux at the surface
  subroutine  SurfaceSalinityFlux(VEC,X,Y,Z,salinity, &
       NONODS,SNLOC,SNGI, &
       SN,SNLX,SNLY,SWEIGH)
    integer, intent(in)::NONODS,SNLOC,SNGI
    real, intent(inout)::VEC(NONODS)
    real, intent(in), dimension(nonods)::X,Y,Z,salinity
    real, intent(in)::SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI),SWEIGH(SNGI)

    REAL scalars(4)
    integer i, nid
    
    REAL, ALLOCATABLE, DIMENSION(:)::salinity_flux

    real outward_water_flux
    real longitude, latitude

#ifdef DEBUG_OCEAN_SURFACE_FORCING
    integer, allocatable, dimension(:)::elementTypes, elementSizes
    integer, save::NCalls=0
    real tmp_f(nonods), tmp_e(nonods), tmp_p(nonods)
    character(len=1024) str

    tmp_f = 0.0
    tmp_e = 0.0
    tmp_p = 0.0
#endif

    ewrite(3, *) "subroutine SurfaceSalinityFlux()"
    return

#ifndef DOUBLEP
    FLAbort("ERA-40 interface does not support single-precision.")
#endif
    
    allocate(salinity_flux(nonods))
    salinity_flux = 0.0

    ! ... salinity flux as defined by Adrian E. Gill,
    ! "Atmosphere-Ocean Dynamics", page 34
    call fluxes_clearfields()
    call fluxes_addfieldofinterest("e",    1)
    call fluxes_addfieldofinterest("tp",   2)
    
    ! Get coordinates in (longitude, latitude) and interpolate.
    do i=1, size(SurfaceNodes)
       nid = SurfaceNodes(i)
       call LongitudeLatitude( (/ x(nid), y(nid), z(nid) /), longitude, latitude)
       call fluxes_getscalars(longitude, latitude, scalars)
          
       ! Take care of the sign of evaporation!
       outward_water_flux = (-scalars(1)-scalars(2))/21600.0
          
#ifdef DEBUG_OCEAN_SURFACE_FORCING
       tmp_f(nid) = outward_water_flux
       tmp_e(nid) = scalars(1)
       tmp_p(nid) = scalars(2)
#endif
          
       salinity_flux(nid) = outward_water_flux*salinity(nid)/density_seawater
    end do
    
    ! Integrate over the top surface
    call integrate_over_surface(NONODS, NSElements, snloc, &
         SENList, X, Y, Z, sngi, SN, SNLX, SNLY, SWEIGH, &
         salinity_flux, vec)

#ifdef DEBUG_OCEAN_SURFACE_FORCING
    str = ' '
    write(str, '("salinity",I6.6,".vtu")') NCalls
    NCalls = NCalls + 1

    call VTKOPEN(str, len_trim(str), "E - P", 5)
    allocate(elementTypes(NSElements), elementSizes(NSElements))
    elementTypes = 5
    elementSizes = 3
    call VTKWRITEMESHD(nonods, NSElements, x, y, z, SENList, elementTypes, elementSizes)
    call VTKWRITEDSN(salinity_flux, "E - P", 5)
    call VTKWRITEDSN(vec, "vec", 3)
    call VTKWRITEDSN(tmp_f, "tmp_f", 5)
    call VTKWRITEDSN(tmp_e, "tmp_e", 5)
    call VTKWRITEDSN(tmp_p, "tmp_p", 5)
    call VTKCLOSE()
    deallocate(elementTypes, elementSizes)
#endif

    deallocate(salinity_flux)

    RETURN

  end subroutine SurfaceSalinityFlux

  ! This subroutine applies the surface heat flux.
  subroutine SurfaceHeatFlux(VEC,X,Y,Z,&
       NONODS,SNLOC,SNGI, &
       SN,SNLX,SNLY,SWEIGH)
    integer, intent(in)::NONODS,SNLOC,SNGI
    real, intent(in), dimension(nonods)::X,Y,Z
    real, intent(in)::SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI),SWEIGH(SNGI)
    real, intent(inout)::VEC(NONODS)

    ! More study is needed to properly model the thermal radiation
    ! (e.g. cloud effects). So using re-analysis data for the moment
    ! real, parameter::stephan_coefficient=5.7E-8 !W/m^2K^4
    ! real, parameter::emissivity=0.985

    real longitude, latitude, Q, scalars(4)
    integer i, nid

    REAL, ALLOCATABLE, DIMENSION(:)::tflux

#ifdef DEBUG_OCEAN_SURFACE_FORCING
    integer, allocatable, dimension(:)::elementTypes, elementSizes
    integer, save::NCalls=0
    real tmp_f(nonods), tmp_qb(nonods), tmp_ql(nonods), &
                    tmp_qa(nonods), tmp_qr(nonods)
    character(len=1024) str

    tmp_f = 0.0
    tmp_qb = 0.0
    tmp_ql = 0.0
    tmp_qa = 0.0
    tmp_qr = 0.0
#endif

    ewrite(3, *) "subroutine SurfaceHeatFlux()"
    return

#ifndef DOUBLEP
    FLAbort("ERA-40 interface does not support single-precision.")
#endif
    
    allocate(TFLUX(nonods))
    tflux = 0.0

    ! Upward heat flux (W m**-2) as defined by Adrian E. Gill,
    ! "Atmosphere-Ocean Dynamics", page 34
    call fluxes_clearfields()

    ! For more details see:
    ! http://tigge.ecmwf.int/tigge/d/show_archive/table=parameters/

    ! Surface latent heat flux (W/m^2s)
    ! Exchange of latent heat (due to phase transitions: evaporation,
    ! condensation) with the surface through turbulent diffusion.
    ! Validity: accumulated from the beginning of the forecast
    ! The flux sign convention is positive downwards.
    call fluxes_addfieldofinterest("slhf", 4)

    ! Surface sensible heat flux (W/m^2s)
    ! Exchange of heat (no phase transitions) with the surface through turbulent
    ! diffusion.
    ! Validity: accumulated from the beginning of the forecast
    ! The flux sign convention is positive downwards.
    call fluxes_addfieldofinterest("sshf", 4)

    ! Surface solar radiation (W/m^2s)
    call fluxes_addfieldofinterest("ssr",  3)
    
    ! Surface thermal radiation (W/m^2s)
    ! Net thermal (short wave) radiation at the surface.
    call fluxes_addfieldofinterest("str",  3)

    ! Get coordinates in (longitude, latitude) and interpolate.
    do i=1, size(SurfaceNodes)
       nid = SurfaceNodes(i)
       call LongitudeLatitude( (/ x(nid), y(nid), z(nid) /), longitude, latitude)
       call fluxes_getscalars(longitude, latitude, scalars)
       
       ! These are accumulated over 6 hours (see units above)
       ! therefore we have to divide by 6*60*60.
       scalars = scalars/21600.0
       
       ! Surface thermal radiation (black-body/long wavelength)
       ! This has to be changed to use the CORE bulk formula
       ! Qb = emissivity*stephan_coefficient*(Temperature(nid)+273.15)**4
       
       ! Total heat flux into ocean
       Q = scalars(1)+scalars(2)+scalars(3)+scalars(4)

#ifdef DEBUG_OCEAN_SURFACE_FORCING
       tmp_f(nid) = Q
       tmp_ql(nid) = scalars(1)
       tmp_qa(nid) = scalars(2)
       tmp_qr(nid) = scalars(3)
       tmp_qb(nid) = scalars(4)
#endif
          
       TFLUX(nid) = Q/(density_seawater*Cp)
    end do
    
    ! Integrate over the top surface
    call integrate_over_surface(NONODS, NSElements, snloc, &
         SENList, X, Y, Z, sngi, SN, SNLX, SNLY, SWEIGH, &
         TFLUX, vec)
    
#ifdef DEBUG_OCEAN_SURFACE_FORCING
    str = ' '
    write(str, '("tflux",I6.6,".vtu")') NCalls
    NCalls = NCalls + 1
    
    call VTKOPEN(str, len_trim(str), "tflux", 5)
    allocate(elementTypes(NSElements), elementSizes(NSElements))
    elementTypes = 5
    elementSizes = 3
    call VTKWRITEMESHD(nonods, NSElements, x, y, z, SENList, elementTypes, elementSizes)
    call VTKWRITEDSN(tflux, "tflux", 5)
    call VTKWRITEDSN(vec, "vec", 3)
    call VTKWRITEDSN(tmp_f, "tmp_f", 5)
    call VTKWRITEDSN(tmp_qb, "qb", 2)
    call VTKWRITEDSN(tmp_ql, "slhf", 4)
    call VTKWRITEDSN(tmp_qa, "sshf", 4)
    call VTKWRITEDSN(tmp_qr, "ssr", 3)
    call VTKCLOSE()
    deallocate(elementTypes, elementSizes)
#endif

    deallocate(tflux)

    RETURN

  end subroutine SurfaceHeatFlux

  ! Calculates the salinity relaxation at the surface
  subroutine  SurfaceSalinityRelaxation(X,Y,Z,&
       snloc, sngi, sn, snlx, snly, sweigh, &
       salinity,ML,dt,centrm,BIGM, ncolm, VEC)
    real, intent(in), dimension(:)::X,Y,Z
    integer, intent(in)::snloc, sngi
    real, intent(in)::sn(snloc, sngi), snlx(snloc, sngi), snly(snloc, sngi)
    real, intent(in)::sweigh(sngi) 
    real, intent(in)::salinity(:)
    real, intent(in)::dt
    integer, intent(in)::centrm(:), ncolm
    real, intent(inout)::ml(:), BIGM(:), VEC(:)

    integer i, nid
    
    real s_climate, relax

    integer xnonod
    real, allocatable, dimension(:)::sml

#ifdef DEBUG_OCEAN_SURFACE_FORCING
    integer, allocatable, dimension(:)::elementTypes, elementSizes
    integer, save::NCalls=0
    real, allocatable, dimension(:)::climate
    character(len=1024) str
#endif

    ewrite(3, *) "subroutine SurfaceSalinityRelaxation()"

    xnonod = size(x)
    allocate(sml(xnonod))
    call get_sml(snloc, sngi, sn, snlx, snly, sweigh, x, y, z, sml)
#ifdef DEBUG_OCEAN_SURFACE_FORCING
    allocate(climate(xnonod))
    climate = 0.0
#endif

#ifndef DOUBLEP
    FLAbort("Climatology interface does not support single-precision.")
#endif

    ! relax = 50.0/(15.0*24*3600)  ! 15 days relaxation time scale
    relax = 50.0/(3.0*24*3600)  ! 3 days relaxation time scale
    
    ! Get coordinates in (longitude, latitude) and interpolate.
    do i=1, size(SurfaceNodes)
       nid = SurfaceNodes(i)
       call climatology_GetSurfaceValue("salinity", (/ x(nid), y(nid), z(nid) /), s_climate)
#ifdef DEBUG_OCEAN_SURFACE_FORCING
       climate(nid) = s_climate
#endif
       VEC(nid) = VEC(nid) + sml(nid)*relax*(s_climate - salinity(nid))
       if(ncolm.ne.0) then
          BIGM(centrm(nid)) = BIGM(centrm(nid)) + dt*sml(nid)*relax
       end if
       ml(nid) = ml(nid) + dt*sml(nid)*relax
    end do
    
#ifdef DEBUG_OCEAN_SURFACE_FORCING
    str = ' '
    write(str, '("salinity_climate",I6.6,".vtu")') NCalls
    NCalls = NCalls + 1

    call VTKOPEN(str, len_trim(str), "salinity", 8)
    allocate(elementTypes(NSElements), elementSizes(NSElements))
    elementTypes = 5
    elementSizes = 3
    call VTKWRITEMESHD(xnonod, NSElements, x, y, z, &
         SENList, elementTypes, elementSizes)
    call VTKWRITEDSN(climate, "salinity", 8)
    call VTKWRITEDSN(sml, "sml", 3)
    call VTKCLOSE()
    deallocate(elementTypes, elementSizes, climate)
#endif

    deallocate(sml)
    return
  end subroutine SurfaceSalinityRelaxation
  
  ! This subroutine applies the surface temperature relaxation.
  subroutine SurfaceTemperatureRelaxation(X,Y,Z, &
       snloc, sngi, sn, snlx, snly, sweigh, &
       Temperature, ML, dt, centrm,&
       BIGM, ncolm, VEC)
    real, intent(in), dimension(:)::X,Y,Z
    integer, intent(in)::snloc, sngi
    real, intent(in)::sn(snloc, sngi), snlx(snloc, sngi), snly(snloc, sngi)
    real, intent(in)::sweigh(sngi)
    real, intent(in)::Temperature(:)
    real, intent(in)::dt
    integer, intent(in)::centrm(:), ncolm
    real, intent(inout)::ML(:), BIGM(:), VEC(:)
    
    integer i, nid
    
    real t_climate, relax

    integer xnonod
    real, allocatable, dimension(:)::sml
    
#ifdef DEBUG_OCEAN_SURFACE_FORCING
    integer, allocatable, dimension(:)::elementTypes, elementSizes
    integer, save::NCalls=0
    real, allocatable, dimension(:)::climate
    character(len=1024) str

#endif

    ewrite(3, *) "subroutine SurfaceTemperatureRelaxation()"
    
#ifndef DOUBLEP
    FLAbort("Climatology interface does not support single-precision.")
#endif
    
    xnonod = size(x)
    allocate(sml(xnonod))
    call get_sml(snloc, sngi, sn, snlx, snly, sweigh, x, y, z, sml)
#ifdef DEBUG_OCEAN_SURFACE_FORCING
    allocate(climate(xnonod))
    climate = 0.0
#endif
    
    ! relax = 50.0/(15.0*24*3600) ! 15 days relaxation time scale
    relax = 50.0/(3.0*24*3600) ! 3 days relaxation time scale
    
    do i=1, size(SurfaceNodes)
       nid = SurfaceNodes(i)
       call climatology_GetSurfaceValue("temperature", (/ x(nid), y(nid), z(nid) /), t_climate)
#ifdef DEBUG_OCEAN_SURFACE_FORCING
       climate(nid) = t_climate
#endif
       VEC(nid) = VEC(nid) + sml(nid)*relax*(t_climate - Temperature(nid))
       if(ncolm.ne.0) then
          BIGM(centrm(nid)) = BIGM(centrm(nid)) + dt*sml(nid)*relax
       end if
       ml(nid) = ml(nid) + dt*sml(nid)*relax
    end do
    
#ifdef DEBUG_OCEAN_SURFACE_FORCING
    str = ' '
    write(str, '("temperature_climate",I6.6,".vtu")') NCalls
    NCalls = NCalls + 1

    call VTKOPEN(str, len_trim(str), "temperature", 11)
    allocate(elementTypes(NSElements), elementSizes(NSElements))
    elementTypes = 5
    elementSizes = 3
    call VTKWRITEMESHD(xnonod, NSElements, x, y, z, &
         SENList, elementTypes, elementSizes)
    call VTKWRITEDSN(climate, "temperature", 11)
    call VTKWRITEDSN(sml, "sml", 3)
    call VTKCLOSE()
    deallocate(elementTypes, elementSizes, climate)
#endif

    deallocate(sml)    
    
    return
  end subroutine SurfaceTemperatureRelaxation

  real function sponge_invtau(x, y, z)
    real, intent(in)::x, y, z
    
    real longitude, latitude, alpha
    
    real::invtau_max=1.0/(5.0*24*3600)  ! 5 days
    
    real, parameter::lat_inside_n=60, lat_outside_n=80
    real, parameter::lat_inside_s=10, lat_outside_s=0
    
    call LongitudeLatitude( (/ x, y, z /), longitude, latitude)
    
    if((lat_inside_s.lt.latitude).and.(latitude.lt.lat_inside_n)) then
       sponge_invtau = -1.0
    else if(latitude.ge.lat_inside_n) then
       ! Linear profile
       ! sponge_invtau = invtau_max*(latitude-lat_inside_n)/(lat_outside_n-lat_inside_n)
     
       ! Exponential profile. alpha = log(2.0)/d sets the distance
       ! from the boundary for which the value of the strength of the sponge
       ! has halved.
       alpha = log(2.0)/3.0
       sponge_invtau = invtau_max*exp(-alpha*(lat_outside_n-latitude))
    else ! if(latitude.le.lat_inside_s) then
       ! Linear profile
       ! sponge_invtau = invtau_max*(lat_inside_s-latitude)/(lat_inside_s-lat_outside_s)

       ! Exponential profile.
       alpha = log(2.0)/1.5
       sponge_invtau = invtau_max*exp(-alpha*(latitude-lat_outside_s))
    end if

    return
  end function sponge_invtau

  subroutine RelaxSalinityToClimatology(NONODS, X, Y, Z, Salinity, &
       ML, tempML, dt, centrm, BIGM, ncolm, VEC)
    integer, intent(in)::NONODS
    real, intent(in), dimension(nonods)::X, Y, Z, Salinity, ML
    real, intent(inout)::tempML(nonods)
    real, intent(in)::dt
    integer, intent(in)::centrm(nonods), ncolm
    real, intent(inout)::BIGM(:), VEC(nonods)
    
    INTEGER I
    
    ! Relaxation parameter
    real relax
    
    ! Salinity from climatology
    real s_climate
    
    ! Get coordinates in (longitude, latitude) and interpolate.
    do i=1, nonods
       relax = sponge_invtau(x(i), y(i), z(i))
       if(relax.gt.0.0) then
          call climatology_GetValue("salinity", x(i), y(i), z(i), s_climate)
          
          VEC(i) = VEC(i) + ML(i)*relax*(s_climate - Salinity(i))
          if(ncolm.ne.0) then
             BIGM(centrm(i)) = BIGM(centrm(i)) + dt*ML(i)*relax
          else
             tempML(i)=tempML(i)+dt*ML(i)*relax
          end if
       end if
    end do

    return
  end subroutine RelaxSalinityToClimatology

  subroutine RelaxTemperatureToClimatology(NONODS, X, Y, Z, Temperature, &
       ML, tempML, dt, centrm, BIGM, ncolm, VEC)
    integer, intent(in)::NONODS
    real, intent(in), dimension(nonods)::X, Y, Z, Temperature, ML
    real, intent(inout)::tempML(nonods)
    real, intent(in)::dt
    integer, intent(in)::centrm(nonods), ncolm
    real, intent(inout)::BIGM(:), VEC(nonods)
    
    INTEGER I
    
    ! Relaxation parameter
    real relax
    
    ! Temperature from climatology
    real t_climate

    ! Get coordinates in (longitude, latitude) and interpolate.
    do i=1, nonods
       relax = sponge_invtau(x(i), y(i), z(i))
       if(relax.gt.0.0) then
          call climatology_GetValue("temperature", x(i), y(i), z(i), t_climate)
          
          VEC(i) = VEC(i) + ML(i)*relax*(t_climate - Temperature(i))
          if(ncolm.ne.0) then
             BIGM(centrm(i)) = BIGM(centrm(i)) + dt*ML(i)*relax
          else
             tempML(i)=tempML(i)+dt*ML(i)*relax
          end if
       end if
    end do

    return

  end subroutine RelaxTemperatureToClimatology
  
  subroutine wind_forcing(state, rhs)
    !!< Implements wind forcing from new options.
    type(state_type), intent(in):: state
    type(vector_field), intent(inout):: rhs
    
    type(vector_field), pointer:: velocity, positions, wind_surface_field
    type(scalar_field), pointer:: wind_drag_coefficient
    type(element_type) faceshape
    character(len=OPTION_PATH_LEN) bctype
    real, dimension(:,:), allocatable:: llpos, llpos_at_quads, ll_wind_at_nodes, ll_wind_at_quads
    real, dimension(:,:), allocatable:: Q, wind_at_quads, wind_at_nodes
    real, dimension(:), allocatable:: detwei, C_D, unorm
    real rho_air
    logical apply_wind_formula, on_sphere
    integer, dimension(:), allocatable:: face_nodes
    integer, dimension(:), pointer:: surface_element_list
    integer snloc, sngi, wdim, nobcs
    integer i, j, k, sele
    logical:: parallel_dg
    
    ewrite(1,*) 'Inside wind_forcing'
    ewrite_minmax(rhs%val(1)%ptr)
    ewrite_minmax(rhs%val(2)%ptr)
   
    velocity => extract_vector_field(state, "Velocity")
    positions => extract_vector_field(state, "Coordinate")
   
    parallel_dg=continuity(velocity)<0 .and. IsParallel()
    on_sphere=have_option('/geometry/spherical_earth')
    faceshape=face_shape(velocity, 1)
    snloc=face_loc(velocity, 1)
    sngi=face_ngi(velocity, 1)
    if (on_sphere) then
      wdim=3 ! dimension of the wind
      assert(velocity%dim==3)
    else
      wdim=velocity%dim-1 ! dimension of the wind
    end if
    
    allocate(detwei(1:sngi), Q(1:snloc,1:snloc), &
      face_nodes(1:snloc), C_D(1:sngi), unorm(1:sngi), &
      wind_at_quads(1:wdim,1:sngi), wind_at_nodes(1:wdim,1:snloc))
    if (on_sphere) then
      allocate( llpos(1:2,1:snloc), llpos_at_quads(1:2,1:sngi), &
        ll_wind_at_nodes(1:2,1:snloc), ll_wind_at_quads(1:2,1:sngi) )
    end if
    
    nobcs = get_boundary_condition_count(velocity)
    do i=1, nobcs
      call get_boundary_condition(velocity, i, type=bctype, &
          surface_element_list=surface_element_list)
      if (bctype=='wind_forcing') then  
        wind_surface_field => extract_surface_field(velocity, i, "WindSurfaceField")
        apply_wind_formula=has_scalar_surface_field(velocity, i, "WindDragCoefficient")
        if (apply_wind_formula) then
           wind_drag_coefficient => extract_scalar_surface_field(velocity, &
              i, "WindDragCoefficient")
           call get_option(trim(velocity%option_path)// &
              '/prognostic/boundary_conditions['//int2str(i-1)//']&
              &/type[0]/wind_velocity/density_air', rho_air)
        end if
              
        do j=1, size(surface_element_list) 
          sele=surface_element_list(j) ! face/surface element nr.
          
          if (parallel_dg) then
            if (.not. element_owned(velocity, face_ele(velocity,sele))) cycle
          end if
          
          ! compute integration weights detwei
          call transform_facet_to_physical(positions, sele, detwei)
          
          ! compute wind_at_nodes: specified wind stress or wind velocity
          if (on_sphere) then
            call LongitudeLatitude(face_val(positions,sele), &
               llpos(1,:), llpos(2,:))
            ll_wind_at_nodes=ele_val(wind_surface_field, j)
            call ll2r3_rotate(llpos(1,:), llpos(2,:), &
                ll_wind_at_nodes(1,:), ll_wind_at_nodes(2,:), &
                wind_at_nodes(1,:), wind_at_nodes(2,:), wind_at_nodes(3,:))
          else
            wind_at_nodes=ele_val(wind_surface_field, j)
          end if
          
          if (apply_wind_formula) then
            ! make wind_at nodes relative wind velocity:
            do k=1, wdim
               wind_at_nodes(k,:)=wind_at_nodes(k,:)-face_val(velocity, k, sele)
            end do

            ! drag coefficient:
            C_D=ele_val_at_quad(wind_drag_coefficient, j)
            
            ! compute wind_at_quads: wind velocity at quadr. points
            if (on_sphere) then
              call LongitudeLatitude(face_val_at_quad(positions,sele), &
                 llpos_at_quads(1,:), llpos_at_quads(2,:))
              ll_wind_at_quads=ele_val_at_quad(wind_surface_field, j)
              call ll2r3_rotate( &
                  llpos_at_quads(1,:), llpos_at_quads(2,:), &
                  ll_wind_at_quads(1,:), ll_wind_at_quads(2,:), &
                  wind_at_quads(1,:), wind_at_quads(2,:), wind_at_quads(3,:))
            else
            
              wind_at_quads=ele_val_at_quad(wind_surface_field, j)
            end if

            ! make wind_at_quads relative velocity (specified wind-sea surface velocity)
            ! at quadrature points
            do k=1, wdim
              wind_at_quads(k,:)=wind_at_quads(k,:)- &
                  face_val_at_quad(velocity, sele, dim=k)
            end do

            ! compute its norm, sum over dim=1 is sum over components
            unorm=sqrt(sum(wind_at_quads**2, dim=1))

            ! multiply at each gauss point:
            detwei=detwei*C_D*unorm*rho_air
          end if
          Q = shape_shape(faceshape, faceshape, detwei)
            
          ! add surface forcing in rhs of momentum equation:
          face_nodes=face_global_nodes(velocity, sele)
          do k=1, wdim
            call addto(rhs, k, face_nodes, matmul(Q, wind_at_nodes(k, :)) )
          end do
        end do
      end if
    end do
      
    deallocate(detwei, Q, face_nodes, C_D, unorm, &
       wind_at_quads, wind_at_nodes)
    if (on_sphere) then
      deallocate(llpos, llpos_at_quads, ll_wind_at_nodes, ll_wind_at_quads)
    end if
    
    ewrite_minmax(rhs%val(1)%ptr)
    ewrite_minmax(rhs%val(2)%ptr)
   
  end subroutine wind_forcing

  ! This is a generic subroutine to calculate some scalar value over a
  ! surface. This subroutine assumes that multiple terms may be
  ! involved in this integration, therefore, the array "integral" is
  ! not zeroed.
  subroutine integrate_over_surface(NNodes, NSElements, snloc, &
       SENList, X, Y, Z, sngi, SN, SNLX, SNLY, SWEIGH, scalar, &
       integral)
    integer, intent(in)::NNodes, NSElements, snloc, &
         SENList(NSElements*snloc), sngi
    real, intent(in)::X(NNodes), Y(NNodes), Z(NNodes), &
         SN(snloc,sngi), &
         SNLX(snloc,sngi),SNLY(snloc,sngi), &
         SWEIGH(sngi), scalar(NNodes)
    real, intent(inout)::integral(NNodes)
    
    integer ielem, igl, L, iloc, jloc, nid_i, nid_j, gi
    real DXDLX, DXDLY, DYDLX, &
         DYDLY, DZDLX, DZDLY
    real A, B, C
    real DETWEI(sngi), rnn

    ewrite(3, *) "subroutine integrate_over_surface( ... )"

    !     Integrate over surface
    do ielem = 1,NSElements
       do GI = 1,sngi
          DXDLX = 0.
          DYDLX = 0.
          DZDLX = 0.

          DXDLY = 0.
          DYDLY = 0.
          DZDLY = 0.

          do L = 1,snloc
             igl   = SENList((ielem-1)*snloc+L)
             DXDLX = DXDLX + SNLX(L,GI)*X(igl)
             DYDLX = DYDLX + SNLX(L,GI)*Y(igl)
             DZDLX = DZDLX + SNLX(L,GI)*Z(igl) 

             DXDLY = DXDLY + SNLY(L,GI)*X(igl)                
             DYDLY = DYDLY + SNLY(L,GI)*Y(igl)   
             DZDLY = DZDLY + SNLY(L,GI)*Z(igl)
          end do

          A = DYDLX*DZDLY - DYDLY*DZDLX
          B = DXDLX*DZDLY - DXDLY*DZDLX
          C = DXDLX*DYDLY - DXDLY*DYDLX
          DETWEI(GI) = SQRT(A*A + B*B + C*C)*SWEIGH(GI)
       end do

       do ILOC=1,snloc
          nid_i=SENList((ielem-1)*snloc+ILOC)
          do JLOC=1,snloc
             nid_j = SENList((ielem-1)*snloc+JLOC)
             rnn = 0.
             do GI=1,sngi
                rnn=rnn+SN(ILOC,GI)*SN(JLOC,GI)*DETWEI(GI)
             end do
             integral(nid_i) = integral(nid_i) + rnn*scalar(nid_j)
          end do
       end do
    end do

    return
  end subroutine integrate_over_surface

end module OceanSurfaceForcing
