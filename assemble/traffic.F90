!    Copyright (C) 2006 Imprial College London and others.
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
  
module traffic

  use pickers
  use FLDebug
  use parallel_tools
  use spud
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN, current_time, dt
  use state_module
  use boundary_conditions
  use generic_interface
  use fields
  use elements
  use parallel_tools
  use qmesh_module, only: do_adapt_mesh
  use vector_tools
  use field_options

  implicit none
  
  private
  public :: traffic_tracer, traffic_source, traffic_density_update, traffic_check_options
  
contains
  
  !============================================================================
  !============================== TRAFFIC TRACER ==============================
  !============================================================================
  subroutine traffic_tracer(field_name,state,timestep)
    use FLDebug
    implicit none
    
    character(len=*), intent(in) :: field_name
    integer,intent(in):: timestep
    logical,save:: lfiem_t
    type(state_type), intent(inout):: state
    type(scalar_field), pointer:: s_field
    character(len=option_path_len):: name

    real, save :: interxsou,interysou,interzsou
    real, save :: interxsoue,interysoue,interzsoue
    real, save :: interxxsou,interxysou,interxzsou
    real, save :: interyxsou,interyysou,interyzsou
    real, save :: interzxsou,interzysou,interzzsou
    real, save :: interxxsoue,interxysoue,interxzsoue
    real, save :: interyxsoue,interyysoue,interyzsoue
    real, save :: interzxsoue,interzysoue,interzzsoue
    real, save :: interintsou1,interintsou2,interintsou3
    real, save :: interintsou4,interintsou5,interintsou6
    real, save :: intert0sou,intertwidso

    integer:: isou,totele,nloc,nnodp,nonods,cars

    ewrite(3,*) "setting the traffic tracer source"

    lfiem_t=.true.

    call get_option("/traffic_model/number_of_vehicles",cars)

    name=trim(field_name)//"Source"
    s_field=>extract_scalar_field(state,name)

    nloc=ele_loc(s_field,1)
    totele=ele_count(s_field)
    
    nonods = s_field%mesh%nodes
    nnodp = nowned_nodes(s_field)

    do isou=1,cars

          call trafmodel(isou,                         &
               interxsou,interysou,interzsou,          &
               interxxsou,interxysou,interxzsou,       &
               interyxsou,interyysou,interyzsou,       &
               interzxsou,interzysou,interzzsou,       &
               interxsoue,interysoue,interzsoue,       &
               interxxsoue,interxysoue,interxzsoue,    &
               interyxsoue,interyysoue,interyzsoue,    &
               interzxsoue,interzysoue,interzzsoue,    &
               interintsou1,interintsou2,interintsou3, &
               interintsou4,interintsou5,interintsou6, &
               intert0sou,intertwidso,                 &
               state,timestep)

          if ( abs(intert0sou-current_time) <= 1.e-10 ) then

             call source_distribution(                    &
                  interxsoue,interysoue,interzsoue,       &
                  interxxsoue,interxysoue,interxzsoue,    &
                  interyxsoue,interyysoue,interyzsoue,    &
                  interzxsoue,interzysoue,interzzsoue,    &
                  interintsou1,interintsou2,interintsou3, &
                  interintsou4,interintsou5,interintsou6, &
                  intert0sou,intertwidso,                 &
                  nnodp,totele,nloc,state,lfiem_t,        &
                  s_field)
             
             lfiem_t=.false.

          end if
    enddo

    return
  end subroutine traffic_tracer

  !============================================================================
  !============================== TRAFFIC SOURCE ==============================
  !============================================================================
  subroutine traffic_source(state,timestep)
    use fldebug
    implicit none

    logical,save:: lfiem_s
    type(state_type),intent(inout):: state
    type(vector_field),pointer:: source,absorp,VelSource
    type(scalar_field),pointer:: denpt,cdenpt
    type(scalar_field):: sourceX,sourceY,sourceZ
    type(scalar_field):: absorpX
    type(scalar_field),pointer:: SPhase

    real, save :: interxsou,interysou,interzsou
    real, save :: interxsoue,interysoue,interzsoue
    real, save :: interxxsou,interxysou,interxzsou
    real, save :: interyxsou,interyysou,interyzsou
    real, save :: interzxsou,interzysou,interzzsou
    real, save :: interxxsoue,interxysoue,interxzsoue
    real, save :: interyxsoue,interyysoue,interyzsoue
    real, save :: interzxsoue,interzysoue,interzzsoue
    real, save :: interintsou1,interintsou2,interintsou3
    real, save :: interintsou4,interintsou5,interintsou6
    real, save :: intert0sou,intertwidso

    integer:: isou,totele,nloc,nonods,nnodp,nod,cars,timestep

    ewrite(3,*) "setting the traffic source"

    call get_option("/traffic_model/number_of_vehicles",cars)

    source  => extract_vector_field(state, "VelocitySource")
    if(.not. have_option(trim(source%option_path) // "/diagnostic/algorithm::Internal")) then
      FLAbort("Velocity source must be an internal diagnostic when using the traffic model")
    end if
    sourceX = extract_scalar_field(source, 1)
    sourceY = extract_scalar_field(source, 2)
    sourceZ = extract_scalar_field(source, 3)

    nonods = source%mesh%nodes
    nnodp = nowned_nodes(source)
    
    lfiem_s=.true.
    
    VelSource => extract_vector_field(state, "SolidVelocity")

    absorp  => extract_vector_field(state, "VelocityAbsorption")
    if(.not. have_option(trim(absorp%option_path) // "/diagnostic/algorithm::Internal")) then
      FLAbort("Velocity absorption must be an internal diagnostic when using the traffic model")
    end if
    absorpX = extract_scalar_field(absorp, 1)

    SPhase => extract_scalar_field(state,"SolidPhase")

    denpt  => extract_scalar_field(state, "Density")
    cdenpt => extract_scalar_field(state, "CopyofDensity")

    nloc=ele_loc(source,1)
    totele=ele_count(source)

    do isou=1,cars

       call trafmodel(isou,                         &
            interxsou,interysou,interzsou,          &
            interxxsou,interxysou,interxzsou,       &
            interyxsou,interyysou,interyzsou,       &
            interzxsou,interzysou,interzzsou,       &
            interxsoue,interysoue,interzsoue,       &
            interxxsoue,interxysoue,interxzsoue,    &
            interyxsoue,interyysoue,interyzsoue,    &
            interzxsoue,interzysoue,interzzsoue,    &
            interintsou1,interintsou2,interintsou3, &
            interintsou4,interintsou5,interintsou6, &
            intert0sou,intertwidso,                 &
            state,timestep)

       if ( abs(intert0sou-current_time) <= 1.e-10 ) then
          
          call source_distribution(                   &
               interxsou,interysou,interzsou,         &
               interxxsou,interxysou,interxzsou,      &
               interyxsou,interyysou,interyzsou,      &
               interzxsou,interzysou,interzzsou,      &
               interintsou1,interintsou2,interintsou3,&
               interintsou4,interintsou5,interintsou6,&
               intert0sou,intertwidso,                &
               nnodp,totele,nloc,state,lfiem_s,       &
               sourceX,sourceY,sourceZ,absorpX,SPhase)

          source%val(1)%ptr(:)=SourceX%val(:)
          source%val(2)%ptr(:)=SourceY%val(:)
          source%val(3)%ptr(:)=SourceZ%val(:)
          absorp%val(1)%ptr(:)=absorpX%val(:)
          absorp%val(2)%ptr(:)=absorpX%val(:)
          absorp%val(3)%ptr(:)=absorpX%val(:)
          
          lfiem_s=.false.

       end if
    enddo

    ! define velocities...
    do nod=1,nonods       
       VelSource%val(1)%ptr(nod)=node_val(source,1,nod)/max(1.e-10,node_val(absorp,1,nod))
       VelSource%val(2)%ptr(nod)=node_val(source,2,nod)/max(1.e-10,node_val(absorp,2,nod))
       VelSource%val(3)%ptr(nod)=node_val(source,3,nod)/max(1.e-10,node_val(absorp,3,nod))
       
       call set(Sphase,nod,min(0.999,node_val(SPhase,nod)))
       call set(cdenpt,nod,denpt%val(nod))
       call set(denpt,nod,(1.-node_val(Sphase,nod))*node_val(cdenpt,nod))
    end do

    ewrite(3,*) "done setting the traffic source for",cars,"vehicles"

    return
  end subroutine traffic_source

  subroutine source_distribution(             &
       interxsou,interysou,interzsou,         &
       interxxsou,interxysou,interxzsou,      &
       interyxsou,interyysou,interyzsou,      &
       interzxsou,interzysou,interzzsou,      &
       interintsou1,interintsou2,interintsou3,&
       interintsou4,interintsou5,interintsou6,&
       intert0sou,intertwidso,                &
       nnodp,totele,nloc,state,lfiem,         &
       source1,source2,source3,source4,source5)
    implicit none

    type(scalar_field)::source1
    logical:: lfiem
    type(scalar_field),optional::source2,source3,source4,source5

    type(state_type),intent(inout)::state
    integer:: nonods,nnodp

    real,allocatable,dimension(:):: mlad
    integer:: totele,nloc

    real:: interxsou,interysou,interzsou
    real:: interxxsou,interxysou,interxzsou
    real:: interyxsou,interyysou,interyzsou
    real:: interzxsou,interzysou,interzzsou
    real:: intert0sou,intertwidso
    real:: interintsou1,interintsou2,interintsou3
    real:: interintsou4,interintsou5,interintsou6

    real,dimension(:),pointer    :: x=>null(),y=>null(),z=>null()
    type(vector_field),pointer   :: v_field=>null()
    type(mesh_type),pointer      :: mesh=>null()
    integer,dimension(:), pointer:: ndglno

    mesh => extract_velocity_mesh(state)
    ndglno => mesh%ndglno
    
    v_field => extract_vector_field(state, "Coordinate")
    
    x => v_field%val(1)%ptr
    y => v_field%val(2)%ptr
    z => v_field%val(3)%ptr
    nonods=size(x)
    
    allocate(mlad(nonods))
    call getml_ad(mlad,nonods,totele,nloc,ndglno,x,y,z)
    
    if(lfiem)then 
       call set(source1,0.0)
       if(present(source2))then
          call set(source2,0.0)
          call set(source3,0.0)
          call set(source4,0.0)
          call set(source5,0.0)
       end if
    endif
    
    call addinsour(mlad,                        &
         interxsou,interysou,interzsou,         &
         interxxsou,interxysou,interxzsou,      &
         interyxsou,interyysou,interyzsou,      &
         interzxsou,interzysou,interzzsou,      &
         interintsou1,interintsou2,interintsou3,&
         interintsou4,interintsou5,interintsou6,&
         nonods,nnodp,x,y,z,                    &
         current_time,intert0sou,intertwidso,   &
         source1,source2,source3,source4,source5)
    
    deallocate(mlad)
    
    return
  end subroutine source_distribution

  subroutine addinsour(ml,        &
       xsou,ysou,zsou,            &
       xxsou,xysou,xzsou,         &
       yxsou,yysou,yzsou,         &
       zxsou,zysou,zzsou,         &
       intsou1,intsou2,intsou3,   &
       intsou4,intsou5,intsou6,   &
       nonods,nnodp,x,y,z,        &
       current_time,t0sou,twidso, &
       source1,source2,source3,   &
       source4,source5)
    implicit none

    type(scalar_field)::source1
    type(scalar_field),optional::source2,source3,source4,source5

    real,dimension(nonods):: x,y,z
    integer:: nonods,nnodp
    real   :: ml(nonods)
    real   :: xsou,ysou,zsou
    real   :: xxsou,xysou,xzsou
    real   :: yxsou,yysou,yzsou
    real   :: zxsou,zysou,zzsou
    real   :: intsou1,intsou2,intsou3
    real   :: intsou4,intsou5,intsou6
    real   :: current_time,t0sou,twidso

    real   :: vrot(3,3),drot(3)
    real   :: accsou,s1,s2,s3,s4,s5,f
    integer:: nod
    real,parameter::toler=1.e-10

    ! add the distribution into source

    ! find e-value decomposition...
    call rect_eig(xxsou,xysou,xzsou, &
         yxsou,yysou,yzsou,          &
         zxsou,zysou,zzsou,          &
         vrot,drot)
    
    ! working out integral of source - effectively source strength
    accsou=0.
    do nod=1,nnodp
       accsou=accsou+ml(nod)*                    &
            funsoup(x(nod),y(nod),z(nod),current_time, &
            xsou,ysou,zsou,t0sou,twidso,vrot,drot)
    enddo
    call allsum(accsou)
    
    if(abs(accsou)>toler) then
       ! working out source term for each node
       do nod=1,nonods

          f = funsoup(x(nod),y(nod),z(nod),current_time, &
               xsou,ysou,zsou,t0sou,twidso,vrot,drot)

          s1=0.;s2=0.;s3=0.;s4=0.;s5=0.
          if (present(source2)) then
             s1=node_val(source1,nod)+(intsou2/accsou) * f
             s2=node_val(source2,nod)+(intsou3/accsou) * f
             s3=node_val(source3,nod)+(intsou4/accsou) * f
             s4=node_val(source4,nod)+(intsou5/accsou) * f
             s5=node_val(source5,nod)+(intsou6/accsou) * f
         
             call set(source1,nod,s1)
             call set(source2,nod,s2)
             call set(source3,nod,s3)
             call set(source4,nod,s4)
             call set(source5,nod,s5)
          else
             
             s1=node_val(source1,nod)+(intsou1/accsou) * f
             call set(source1,nod,s1)
          end if
          
       end do
    endif
    
    return
  end subroutine addinsour

  subroutine rect_eig(m11,m12,m13, &
       m21,m22,m23,m31,m32,m33,    &
       vrot,drot)
    implicit none
    
    real, intent(in):: m11,m12,m13,m21,m22,m23,m31,m32,m33
    real   :: vrot(3,3),drot(3),a(3,3)
    
    integer,parameter::ndim=3
    
    ! vrot is the rotation matrix and drot the e-values

    a(1,1) = m11
    a(1,2) = m12
    a(1,3) = m13
    a(2,1) = m21
    a(2,2) = m22
    a(2,3) = m23
    a(3,1) = m31
    a(3,2) = m32
    a(3,3) = m33

    call eigendecomposition_symmetric(a, vrot, drot)
    
    return
  end subroutine rect_eig

  real function funsoup(x,y,z,t, &
       x0,y0,z0,t0,twind,vrot,drot)
    implicit none

    real:: x,y,z,t,x0,y0,z0
    real:: t0,twind
    real:: vrot(3,3),drot(3)
    real:: xx(3),xx0(3)

    ! note: in this subroutine the location of the source is denoted by x0,y0,z0
    ! x, y, z denote the co-ordinates of each node
    ! vrot is the rotation matrix and drot the e-values

    xx(1) = x
    xx(2) = y
    xx(3) = z
    xx0(1)= x0
    xx0(2)= y0
    xx0(3)= z0

    funsoup=funsou(xx,t,xx0,t0,twind,vrot,drot)

    return
  end function funsoup

  real function funsou(x,t,x0,t0,twind,vrot,drot)
    implicit none

    real   :: x(3),t,x0(3),t0,twind
    real   :: vrot(3,3),drot(3)

    ! finding out if the point (x,t) falling in the
    ! specified source region during t0 to t0+dt
    !
    ! (x,t) is the position and time of the computational domain
    ! (x0,t0) is the center position of the source at time t0
    ! twind is the duration that the source exist from t0 to t0+twind
    ! nb a=vrot^t drot vrot
    ! vrot is the rotation matrix and drot the e-values

    funsou = 0.0

    ! outside active time - nothing to do
    if(abs(t-t0)> 0.5*twind) then
       funsou=0.0
       return
    end if

    if (dist_rect(x,x0,vrot,drot)<=1.0) then
       funsou=1.
    else
       funsou=0.
    end if

  end function funsou

    real function dist_rect(x,x0,vrot,drot)
      implicit none

      real   :: x(3),x0(3)
      real   :: vrot(3,3),drot(3)
      real   :: p,vx(3)
      integer:: i,j

      ! calculate the value of quadratic 
      ! functions for a box shape
      ! nb a=v^t d v
      ! v: roation matrix; d: the e-values.

      do i = 1, 3
         p = 0.
         do j = 1, 3
            p= p + vrot(i,j)*(x(j)-x0(j))
         end do
         vx(i) = p
      end do

      dist_rect = max((vx(1)**2)*drot(1),(vx(2)**2)*drot(2), &
           (vx(3)**2)*drot(3) )

    end function dist_rect

    subroutine getml_ad(mlad,nonods,totele,nloc,ndglno,x,y,z)
      implicit none

      integer:: nonods,totele,nloc,ndglno(totele*nloc)
      real   :: mlad(nonods)
      real   :: x(nonods),y(nonods),z(nonods)

      real   :: x0(4),y0(4),z0(4)
      real   :: vol
      integer:: ele,iloc,nod

      ! this sub calculates the lumped mass matrix mlad

      mlad(1:nonods) = 0.0
      do ele=1,totele

         do iloc=1,nloc
            nod=ndglno((ele-1)*4+iloc)
            x0(iloc) = x(nod)
            y0(iloc) = y(nod)
            z0(iloc) = z(nod)
         end do
         vol = abs( tetvol(x0, y0, z0) )

         do iloc=1,nloc
            nod=ndglno((ele-1)*nloc+iloc)
            mlad(nod) = mlad(nod)+vol/4.0
         enddo
      enddo

    end subroutine getml_ad
    
    
    subroutine traffic_density_update(state)
      implicit none
      
      type(state_type)::state
      type(scalar_field),pointer:: denpt,cdenpt
      
      denpt => extract_scalar_field(state, "Density")
      cdenpt=> extract_scalar_field(state, "CopyofDensity")
      
      call set(denpt,cdenpt)
      
      return
    end subroutine traffic_density_update
    
    !============================================================================
    !============================= VISSIM INTERFACE =============================
    !============================================================================
    subroutine trafmodel(isou,                    &
         
         interxsou,interysou,interzsou,           &
         interxxsou,interxysou,interxzsou,        &
         interyxsou,interyysou,interyzsou,        &
         interzxsou,interzysou,interzzsou,        &
         
         interxsoue,interysoue,interzsoue,        &
         interxxsoue,interxysoue,interxzsoue,     &
         interyxsoue,interyysoue,interyzsoue,     &
         interzxsoue,interzysoue,interzzsoue,     &
         
         interintsou1,interintsou2,interintsou3,  &
         interintsou4,interintsou5,interintsou6,  & 
         
         intert0sou,intertwidso,                  &
         
         state,timestep)

      use FLDebug
      use parallel_tools
      use generic_interface
      implicit none

      integer,intent(in):: isou,timestep

      real,intent(out):: interxsou,interysou,interzsou
      real,intent(out):: interxxsou,interxysou,interxzsou
      real,intent(out):: interyxsou,interyysou,interyzsou
      real,intent(out):: interzxsou,interzysou,interzzsou

      real,intent(out):: interxsoue,interysoue,interzsoue
      real,intent(out):: interxxsoue,interxysoue,interxzsoue
      real,intent(out):: interyxsoue,interyysoue,interyzsoue
      real,intent(out):: interzxsoue,interzysoue,interzzsoue

      real,intent(out):: interintsou1,interintsou2,interintsou3
      real,intent(out):: interintsou4,interintsou5,interintsou6

      real,intent(out):: intert0sou,intertwidso

      type(state_type),intent(inout):: state

      real,allocatable,dimension(:,:),save :: trafficd
      real,allocatable,dimension(:)        :: traffic_data
      integer,save                         :: carno,ef,roadbc
      real,save                            :: vehicle_scale,exs,network_scale,lag
      real,save                            :: model_scale,rwsp,mwsp
      logical,save                         :: slfr
      character(len=option_path_len)       :: vehicle_activity_file
      character(len=option_path_len)       :: petrol_car,diesel_car,large_vehicle
      integer                              :: pc_stat,dc_stat,lv_stat
      logical,save                         :: petrol_car_db_exists
      logical,save                         :: diesel_car_db_exists
      logical,save                         :: large_vehicle_db_exists
      logical,save                         :: initz=.false.,have_bounding_box=.false.
      integer                              :: tmpncar

      real,save                            :: xmin,xmax,ymin,ymax
      real,allocatable,dimension(:),save   :: clg,cat
      real,allocatable,dimension(:),save   :: cax,cay,caz
      real,allocatable,dimension(:),save   :: cau,caa
      integer,allocatable,dimension(:),save:: cid,ctp,line
      real,allocatable,dimension(:,:),save :: petrol_car_emissions_factors
      real,allocatable,dimension(:,:),save :: diesel_car_emissions_factors
      real,allocatable,dimension(:,:),save :: large_vehicle_emissions_factors

      integer:: i,j,stat
      integer:: io,icount
      integer:: lnn,ex,ey,ez,carid
      integer:: length_unf

      real:: number_of_reads, relax
      real, dimension(:, :), allocatable :: emissions_factors_db
      real:: vx(3,3),vy(3,3),vz(3,3),d(3,3),rot1(3,3),rot2(3,3)
      real:: nxtx,nxty,nxtz,prvx,prvy,prvz
      real:: tanany,tananz,angy,angz
      real:: cx,cy,cz,crx,cry,crz
      real:: a1,a2,a3,b1,b2,b3,exd
      real,parameter:: pi = 3.1415926535897931

      logical:: carr

      type(vector_field),pointer:: position

      ! this sub is the traffic model interface
      ! VISSIM vehicle activivity output file format ::
      ! VehNr; Type; Length; t; WorldX; WorldY; vMS; a
      ! make sure you have not included pedestrians

      allocate(emissions_factors_db(162, 152))

      if (.not.initz) then

         call get_option("/traffic_model/number_of_vehicles",carno)
         call get_option("/traffic_model/vehicle_activity_file/file_name",vehicle_activity_file)
         call get_option("/traffic_model/traffic_model_time_lag",lag)
         call get_option("/traffic_model/scale/vehicle_scale",vehicle_scale)
         call get_option("/traffic_model/scale/exhaust_diameter",exs)
         call get_option("/traffic_model/scale/road_network_scale",network_scale)
         call get_option("/traffic_model/scale/model_scale",model_scale)
         call get_option("/traffic_model/scale/reference_velocity_field",rwsp)
         call get_option("/traffic_model/scale/reference_velocity_model",mwsp)
         call get_option("/traffic_model/sloping_floor",roadbc,stat)

         if (have_option("/traffic_model/bounding_box")) then
            call get_option("/traffic_model/bounding_box/x_min",xmin)
            call get_option("/traffic_model/bounding_box/x_max",xmax)
            call get_option("/traffic_model/bounding_box/y_min",ymin)
            call get_option("/traffic_model/bounding_box/y_max",ymax)
            have_bounding_box=.true.
         end if

         if (have_option("/traffic_model/emissions_factors/type::modelled_emissions_factors")) then
            call get_option("/traffic_model/emissions_factors/type::modelled_emissions_factors/passenger_vehicle_petrol",petrol_car,pc_stat)
            call get_option("/traffic_model/emissions_factors/type::modelled_emissions_factors/passenger_vehicle_diesel",diesel_car,dc_stat)
            call get_option("/traffic_model/emissions_factors/type::modelled_emissions_factors/learge_vehicle",large_vehicle,lv_stat)

            if (pc_stat==0) then        
               inquire(file = trim(petrol_car), exist = petrol_car_db_exists)
               if (petrol_car_db_exists) then
                  allocate(petrol_car_emissions_factors(162,152))
                  petrol_car_emissions_factors=0.
                  open(unit = 910, file = trim(petrol_car), action = "read")
                  do i=1,162
                     read(910,*) petrol_car_emissions_factors(i,1:152)
                  end do
                  close(910)
               else
                  FLExit("Could not locate the petrol vehicles' emissions factors database files")                  
               end if
            endif

            if (dc_stat==0)then 
               inquire(file = trim(diesel_car), exist = diesel_car_db_exists)
               if (diesel_car_db_exists) then
                  allocate(diesel_car_emissions_factors(162,152))
                  diesel_car_emissions_factors=0.
                  open(unit = 911, file = trim(diesel_car), action = "read")
                  do i=1,162
                     read(911,*) diesel_car_emissions_factors(i,1:152)
                  end do
                  close(911)
               else
                  FLExit("Could not locate the diesel veicles' emissions factors database file")
               end if
            end if

            if (lv_stat==0)then
               inquire(file = trim(large_vehicle), exist = large_vehicle_db_exists)
               if (large_vehicle_db_exists)then
                  allocate(large_vehicle_emissions_factors(162,152))
                  large_vehicle_emissions_factors=0.
                  open(unit = 912, file = trim(large_vehicle), action = "read")
                  do i=1,162
                     read(912,*) large_vehicle_emissions_factors(i,1:152)
                  end do
                  close(912)
               else
                  FLExit("Could not locate the large vehicles' emissions factors database file")
               end if
            end if
         end if
         
         if (stat==0) then
            slfr=.true.
         else
            slfr=.false.
         end if
         
         crx=-66.6; nxtx=-66.6
         cry=-66.6; nxty=-66.6
         crz=-66.6; nxtz=-66.6
         
         io=0
         j =0

         inquire(iolength=length_unf) number_of_reads
         open(unit = 913, file = vehicle_activity_file, status='old', access ='direct', recl=length_unf, form='unformatted')
         read(913,rec=1)number_of_reads
         close(913)

         allocate( traffic_data(int(number_of_reads)+1) )
         traffic_data=-666.66
         inquire(iolength=length_unf) traffic_data
         open(unit = 914, file = vehicle_activity_file, status='old', access ='direct', recl=length_unf, form='unformatted')
         read(914,rec=1) traffic_data

         allocate(trafficd(int(number_of_reads)/8,8))
         do i=1,int(number_of_reads),8
            trafficd(j+1,1:8) = traffic_data(i+1:i+8)
            j=j+1
         end do
         close(914)

         tmpncar= maxval(trafficd(:,1))
         ef     = j-1

         if (carno/=tmpncar) then
            FLExit("Check your vehicle activity file, number of cars don't match")
         endif
         
         allocate(cid(ef),ctp(ef))
         allocate(cax(ef),cay(ef),caz(ef))
         allocate(cat(ef))
         allocate(cau(ef),caa(ef),clg(ef)) 
         
         caz=-66.6
         
         do i=1,ef
            cid(i)=int(trafficd(i,1))                        ! vehicle id
            ctp(i)=int(trafficd(i,2))                        ! vehicle type
            clg(i)=trafficd(i,3)/vehicle_scale               ! vehicle length

            cat(i)=trafficd(i,4)*rwsp/(model_scale*mwsp)+lag ! simulation time

            cax(i)=trafficd(i,5)/network_scale               ! x
            cay(i)=trafficd(i,6)/network_scale               ! y
                        
            cau(i)=trafficd(i,7)                             ! u
            caa(i)=trafficd(i,8)                             ! a
         enddo
         
         if (slfr) then
            call get_floor(state,roadbc,cax,cay,ef,caz)
         endif
         
         initz=.true.
      endif
      
      if (slfr) then   
         if (do_adapt_mesh(current_time,timestep)) then
            call get_floor(state,roadbc,cax,cay,ef,caz)
         endif
      else
         position => extract_vector_field(state,"Coordinate")
         
         crz  = minval(position%val(3)%ptr)
         nxtz = minval(position%val(3)%ptr)
         prvz = minval(position%val(3)%ptr)
         
         call allmin(crz)
         call allmin(nxtz)
         call allmin(prvz)
      endif
      
      allocate(line(ef))
      line   = 0
      icount = 1

      do i=1,ef
         if (cid(i) == isou) then
            line(icount) = i
            icount = icount + 1
         endif
      enddo

      carr = .false.
      do i=1,icount-2
         if ( current_time >= cat(line(i)) .and. current_time < cat(line(i+1)) ) then
            lnn=i
            carr=.true.
         endif
      enddo

      if (carr) then

         relax = (current_time-cat(line(lnn)))/(cat(line(lnn+1))-cat(line(lnn)))

         ! current location
         crx = (1-relax) * cax(line(lnn))  + relax * cax(line(lnn+1))
         cry = (1-relax) * cay(line(lnn))  + relax * cay(line(lnn+1))

         ! check for bounding box
         if (have_bounding_box) then
            if (crx<xmin .or. crx>xmax .or. &
                cry<ymin .or. cry>ymax) then

               !interxsou   =-66; interysou   =-66; interzsou   =-66
               !interxsoue  =-66; interysoue  =-66; interzsoue  =-66
               !interxxsou  =-66; interxysou  =-66; interxzsou  =-66
               !interyxsou  =-66; interyysou  =-66; interyzsou  =-66
               !interzxsou  =-66; interzysou  =-66; interzzsou  =-66
               !interintsou1=-66; interintsou2=-66; interintsou3=-66
               !interintsou4=-66; interintsou5=-66; interintsou6=-66
               intert0sou=-1.e+20
               deallocate(line)
               deallocate(emissions_factors_db)
               return

            end if
         end if

         if (slfr) then
            crz = caz(line(lnn))
         endif

         ! next location
         nxtx = cax(line(lnn+1))
         nxty = cay(line(lnn+1))
         if (slfr) then
            nxtz = caz(line(lnn+1))
         end if
         
         ! previous location
         prvx = cax(line(lnn))
         prvy = cay(line(lnn))
         if (slfr) then
            prvz = caz(line(lnn))
         end if
         
         ! work out direction
         if (crx==nxtx .and. cry==nxty) then
            i=lnn
            do while (nxtx==crx .and. nxty==cry .and. line(i+2)/=0)
               i=i+1
               nxtx = cax(line(i+1))
               nxty = cay(line(i+1))
            end do
 
            if (slfr) then
               nxtz = caz(line(i+1))
            end if
         end if

         tananz = (nxty-cry) / (nxtx-crx)
         tanany = (nxtz-crz) / (sqrt((nxty-cry)**2+(nxtx-crx)**2))

         ! if not first time vehicle appears in sim get previous location
         ! if same as current
         if (crx==prvx .and. cry==prvy .and. lnn>1) then
            prvx = cax(line(lnn-1))
            prvy = cay(line(lnn-1))
         end if

         if (crx/=prvx .and. cry/=prvy) then
            tananz = (tananz + ( (cry-prvy) / (crx-prvx) ) ) / 2.
            tanany = (tanany + ( (crz-prvz) / (sqrt((cry-prvy)**2+(crx-prvx)**2)))  ) / 2.
         end if

         if (crx==nxtx .and. cry==nxty .and. crx/=prvx .and. cry/=prvy ) then
            tananz = (cry-prvy) / (crx-prvx) 
            tanany = (crz-prvz) / (sqrt((cry-prvy)**2+(crx-prvx)**2))
         end if

         ! if there is a bug in the traffic model and this vehicle is
         ! permanently stationaty then ignore it
         if (crx==nxtx .and. cry==nxty .and. crx==prvx .and. cry==prvy ) then
            intert0sou=-1.e+20
            deallocate(line)
            deallocate(emissions_factors_db)
            return
         end if

         ! set interXsou
         interxsou = crx
         interysou = cry

         ! get angles
         angz=atan(tananz)
         if (slfr) then
            angy=atan(tanany)
         else
            angy=0
         end if

         if (crx<nxtx) then
            ex=+1
         else
            ex=-1
         end if
         
         if (cry<nxty) then
            ey=+1
         else
            ey=-1
         end if
         
         if (crz<nxtz) then
            ez=+1
         else
            ez=-1
         end if
         
         if (abs(crx-nxtx)<=1.e-10 .and. cry<nxty) then
            ex=+1
            ey=+1
         else if (abs(crx-nxtx)<=1.e-10 .and. nxty<cry) then
            ex=-1
            ey=-1
         else if (abs(cry-nxty)<=1.e-10 .and. crx<nxtx) then
            ex=+1
            ey=+1
         else if (abs(cry-nxty)<=1.e-10 .and. nxtx<crx) then
            ex=-1
            ey=-1
         end if
         
         ! x, y, dimensions, exhaust & emissions...
         emissions_factors_db=0.
         carid=0
         if (ctp(line(lnn))==100 .or. ctp(line(lnn))==11 .or. &      ! passenger car, taxi, van
              ctp(line(lnn))==9 .or. ctp(line(lnn))==24) then

            if (ctp(line(lnn))==100 .or. ctp(line(lnn))==11 .or. &
                 ctp(line(lnn))==9) then
               ! here taxis are assumed to use petrol...needs to be changed...
               cx = clg(line(lnn))/vehicle_scale
               cy = 1.5           /vehicle_scale
               cz = 1.5           /vehicle_scale
               if (petrol_car_db_exists) then
                  emissions_factors_db = petrol_car_emissions_factors
               end if
            else
               cx = 5.31          /vehicle_scale
               cy = 1.6           /vehicle_scale
               cz = 1.8           /vehicle_scale
               if (diesel_car_db_exists) then
                  emissions_factors_db = diesel_car_emissions_factors
               end if
            end if
            
            a1 = abs((cx/2)*cos(angz))
            a2 = abs((cy/2-0.01*cy)*sin(angz))
            b1 = abs((cx/2)*sin(angz))
            b2 = abs((cy/2-0.01*cy)*cos(angz))

            if (ex==+1 .and. ey==+1) then
               interxsou=-a1+interxsou
               interysou=-b1+interysou

               interxsoue=-a1-a2+interxsou
               interysoue=-b1+b2+interysou

            else if (ex==-1 .and. ey==-1) then
               interxsou=+a1+interxsou
               interysou=+b1+interysou

               interxsoue=+a1+a2+interxsou
               interysoue=+b1-b2+interysou

            else if (ex==+1 .and. ey==-1) then
               interxsou=+a1+interxsou
               interysou=+b1+interysou

               interxsoue=-a1+a2+interxsou
               interysoue=+b1+b2+interysou

            else if (ex==-1 .and. ey==+1) then
               interxsou=-a1+interxsou
               interysou=-b1+interysou

               interxsoue=+a1-a2+interxsou
               interysoue=-b1-b2+interysou
            end if

         else if (ctp(line(lnn))==200 .or. ctp(line(lnn))==15 .or. &  ! bus, ar. bus
              ctp(line(lnn))==16) then

            if (ctp(line(lnn))==200 .or. ctp(line(lnn))==15) then
               cx = 11.541/vehicle_scale
               cy =  2.55 /vehicle_scale
               cz =  4.4  /vehicle_scale
            else
               cx = 18.59 /vehicle_scale
               cy =  2.55 /vehicle_scale
               cz =  2.1  /vehicle_scale
            end if
            if (large_vehicle_db_exists) then
               emissions_factors_db = large_vehicle_emissions_factors
            end if
            
            a1 = abs((cx/2)*cos(angz))
            b1 = abs((cx/2)*sin(angz))

            if (ex==+1 .and. ey==+1) then
               interxsou=-a1+interxsou
               interysou=-b1+interysou

               interxsoue=-a1+interxsou
               interxsoue=-b1+interysou

            else if (ex==-1 .and. ey==-1) then
               interxsou=+a1+interxsou
               interysou=+b1+interysou

               interxsoue=+a1+interxsou
               interxsoue=+b1+interysou

            else if (ex==+1 .and. ey==-1) then
               interxsou=+a1+interxsou
               interysou=+b1+interysou

               interxsoue=-a1+interxsou
               interxsoue=+b1+interysou

            else if (ex==-1 .and. ey==+1) then
               interxsou=-a1+interxsou
               interysou=-b1+interysou

               interxsoue=+a1+interxsou
               interxsoue=-b1+interysou
            end if

         else if (ctp(line(lnn))==300  .or. ctp(line(lnn))==40 .or. & ! HGV, LGV exhaust
              ctp(line(lnn))==201) then
            cx = 10.215/vehicle_scale 
            cy =  2.5  /vehicle_scale
            cz =  4.4  /vehicle_scale
            if (large_vehicle_db_exists) then
               emissions_factors_db = large_vehicle_emissions_factors
            end if

            a1 = abs((cx/2)*cos(angz))
            a2 = abs((cy/2)*sin(angz))
            a3 = abs((cx/5)*cos(angz))
            b1 = abs((cx/2)*sin(angz))
            b2 = abs((cy/2)*cos(angz))
            b3 = abs((cx/5)*sin(angz))

            if(ex==+1 .and. ey==+1)then
               interxsou=-a1+interxsou
               interysou=-b1+interysou

               interxsoue=+a1-a2-a3+interxsou
               interysoue=+b1+b2-b3+interysou

            else if (ex==-1 .and. ey==-1) then
               interxsou=+a1+interxsou
               interysou=+b1+interysou

               interxsoue=-a1+a2+a3+interxsou
               interysoue=-b1-b2+b3+interysou

            else if (ex==+1 .and. ey==-1) then
               interxsou=+a1+interxsou
               interysou=+b1+interysou

               interxsoue=+a1+a2-a3+interxsou
               interysoue=-b1+b2-b3+interysou

            else if (ex==-1 .and. ey==+1) then
               interxsou=-a1+interxsou
               interysou=-b1+interysou

               interxsoue=-a1-a2+a3+interxsou
               interysoue=+b1-b2-b3+interysou
            end if
         else
            intert0sou=-1.e+20
            ewrite(1,*) "Vehicle class not found. Check you vehicle activity file..."
            deallocate(line)
            deallocate(emissions_factors_db)
            return            
         end if

         interzsou  = crz + 2.0/model_scale + cz/2
         interzsoue = crz + 2.0/model_scale 

         if (ez==-1) then
            interxsoue = interzsoue*sin(angy) + interxsoue*cos(angy)
            interzsoue = interzsoue*cos(angy) - interxsoue*sin(angy)
         else if (ez==+1) then
            interxsoue =-interzsoue*sin(angy) + interxsoue*cos(angy)
            interzsoue = interzsoue*cos(angy) + interxsoue*sin(angy)
         end if

         interxxsou = 4/cx**2
         interxysou = 0.
         interxzsou = 0.
         interyxsou = 0.
         interyysou = 4/cy**2
         interyzsou = 0.
         interzxsou = 0.
         interzysou = 0.
         interzzsou = 4/cz**2  

         d(1,1) = interxxsou 
         d(1,2) = interxysou
         d(1,3) = interxzsou
         d(2,1) = interyxsou
         d(2,2) = interyysou
         d(2,3) = interyzsou
         d(3,1) = interzxsou
         d(3,2) = interzysou
         d(3,3) = interzzsou

         ! for rotations about x
         vx(1,1) = 1.
         vx(1,2) = 0.
         vx(1,3) = 0.
         vx(2,1) = 0.
         vx(2,2) = cos(0.)
         vx(2,3) = sin(0.)
         vx(3,1) = 0.
         vx(3,2) =-sin(0.)
         vx(3,3) = cos(0.)

         ! for rotations about y
         vy(1,1) = cos(angy)
         vy(1,2) = 0.
         vy(1,3) =-sin(angy)
         vy(2,1) = 0.
         vy(2,2) = 1.
         vy(2,3) = 0.
         vy(3,1) = sin(angy)
         vy(3,2) = 0.
         vy(3,3) = cos(angy)

         ! for rotations about z
         vz(1,1) = cos(angz)
         vz(1,2) = sin(angz)
         vz(1,3) = 0.
         vz(2,1) =-sin(angz)
         vz(2,2) = cos(angz)
         vz(2,3) = 0.
         vz(3,1) = 0.
         vz(3,2) = 0.
         vz(3,3) = 1.

         rot1 = matmul(matmul(vz,vy),vx)
         rot2 = matmul(matmul(transpose(rot1),d),rot1)

         interxxsou = rot2(1,1)
         interxysou = rot2(1,2)
         interxzsou = rot2(1,3)
         interyxsou = rot2(2,1)
         interyysou = rot2(2,2)
         interyzsou = rot2(2,3)
         interzxsou = rot2(3,1)
         interzysou = rot2(3,2)
         interzzsou = rot2(3,3)

         intertwidso = 10000.
         intert0sou  = current_time

         ! exhaust diameter - fixed for all vehicles...
         exd = exs/model_scale

         interxxsoue = 4/(exd**2)
         interxysoue = 0.
         interxzsoue = 0.
         interyxsoue = 0.
         interyysoue = 4/(exd**2)
         interyzsoue = 0.
         interzxsoue = 0.
         interzysoue = 0.
         interzzsoue = 4/(exd**2)
 
         if (have_option("/traffic_model/emissions_factors/") ) then
            call emfac(interintsou1,cau(line(lnn)),caa(line(lnn)),emissions_factors_db)
         else
            interintsou1=1.
         endif

         interintsou2 = ex*abs(cos(angz)*cos(angy))*((cau(line(lnn))*cx*cy*cz)/(dt*rwsp)) 
         interintsou3 = ey*abs(sin(angz)*cos(angy))*((cau(line(lnn))*cx*cy*cz)/(dt*rwsp))
         interintsou4 = ez*abs(          sin(angy))*((cau(line(lnn))*cx*cy*cz)/(dt*rwsp))
         interintsou5 =(cx*cy*cz)/dt
         interintsou6 = cx*cy*cz

         !ewrite(3,*) "...TRAFFIC MODEL INTERFACE PRINTS..."
         !ewrite(3,*) "time",acctim,acctim/rwsp*model_scale
         !ewrite(3,*) "vehicle type & id:",ctp(line(lnn)),isou
         !ewrite(3,*) "vehicle dimensions:",cx,cy,cz
         !ewrite(3,*) "angles (tanan,ang,sin,cos)"
         !ewrite(3,*) "y::",tanany,angy,":",sin(angy),cos(angy)
         !ewrite(3,*) "z::",tananz,angz,":",sin(angz),cos(angz)
         !ewrite(3,*) "previous,current,next"
         !ewrite(3,*) "x::",prvx,crx,nxtx
         !ewrite(3,*) "y::",prvy,cry,nxty
         !ewrite(3,*) "z::",prvz,crz,nxtz
         !ewrite(3,*) "vehicle coords, inters"
         !ewrite(3,*) interxsou,interysou,interzsou
         !ewrite(3,*) "exhaust coords"
         !ewrite(3,*) interxsoue,interysoue,interzsoue
         !ewrite(3,*) "rotated"
         !ewrite(3,*) rot2(1,1), rot2(1,2), rot2(1,3)
         !ewrite(3,*) rot2(2,1), rot2(2,2), rot2(2,3)
         !ewrite(3,*) rot2(3,1), rot2(3,2), rot2(3,3)
         !ewrite(3,*) "intensities"
         !ewrite(3,*) interintsou1
         !ewrite(3,*) interintsou2,interintsou3,interintsou4
         !ewrite(3,*) interintsou5,interintsou6
         !ewrite(3,*) "signs",ex,ey,ez
         !ewrite(3,*) "...END OF TRAFFIC PRINTS..."
      else
         intert0sou=-1.e+20
      end if

      deallocate(line)
      deallocate(emissions_factors_db)

      return
    end subroutine trafmodel

    !----------------------------------------------------------
    !----------------------------------------------------------

    subroutine emfac(interintsou,speed,accel,db)

      use FLDebug
      implicit none

      real,intent(out)                   :: interintsou
      real,intent(inout)                 :: speed,accel
      real,dimension(162,152),intent(in) :: db

      real   :: spac
      real   :: dec_speed,dec_spac
      integer:: u,ua
      real   :: speed_mat(2,1),spac_mat(1,2),ef(2,2)

      if (have_option("/traffic_model/emissions_factors/type::constant_emissions_factors"))then

         call get_option("/traffic_model/emissions_factors/type::constant_emissions_factors/emissions_factor",interintsou)

      elseif (have_option("/traffic_model/emissions_factors/type::modelled_emissions_factors"))then

         spac=speed*accel
         
         if (speed>160.) speed= 160.
         if (speed<0.)   speed=   0.         
         if (spac>70.)   spac =  70.
         if (spac<-80.)  spac = -80.
         
         do u=2,161
            if (speed>=db(u,1) .and. speed<=db(u+1,1))then
               do ua=2,151
                  if (spac>=db(1,ua) .and. spac<=db(1,ua+1))then

                     dec_speed=speed-int(speed)
                     dec_spac =abs(spac-int(spac))
                     if (spac<0.) dec_spac=1-dec_spac

                     spac_mat(1,1)=1-dec_spac
                     spac_mat(1,2)=dec_spac

                     speed_mat(1,1)=1-dec_speed
                     speed_mat(2,1)=dec_speed
                     
                     ef(1,1)=db(u,ua)
                     ef(1,2)=db(u,ua+1)
                     ef(2,1)=db(u+1,ua)
                     ef(2,2)=db(u+1,ua+1)

                     interintsou=maxval(matmul(matmul(spac_mat,ef),speed_mat))
                  endif
               enddo
            endif
         enddo     
      else
         interintsou=1.
      endif

      return
    end subroutine emfac
    
    !----------------------------------------------------------

    subroutine get_floor(state,roadbc,PointX,PointY,Points,PointZ)

      use FLDebug
      use parallel_tools
      use generic_interface 
      
      implicit none

      type(state_type),intent(inout)    :: state
      integer,intent(in)                :: roadbc
      integer,intent(in)                :: Points
      real,dimension(Points),intent(in) :: PointX,PointY
      real,dimension(Points),intent(out):: PointZ

      character(len=field_name_len)     :: name
      type(mesh_type),pointer           :: xmesh
      type(scalar_field),pointer        :: dist
      real,dimension(:),pointer         :: x,y,z
      integer                           :: xnonod,snloc,sele,ntri
      integer,pointer,dimension(:)      :: surface_element_list
      type(scalar_field)                :: botdis_field
      type(vector_field),pointer        :: positions
      integer,dimension(:),allocatable  :: senlist
      integer                           :: i,l1,l2,l3,pp1
      real                              :: zz

      integer,dimension(:),allocatable  :: eid
      real,dimension(:, :),allocatable     :: eshape
      type(vector_field)                :: bc_position, bc_position_2d
      type(mesh_type), pointer          :: surface_mesh

#ifdef HAVE_MPI
      include 'mpif.h'
      integer:: ierr
#endif

      !ewrite(3,*) "inside r-tree interface",roadbc

      name="DistanceToBottom"
      xmesh => extract_mesh(state,"CoordinateMesh")

      positions => extract_vector_field(state,"Coordinate")
      x => positions%val(1)%ptr
      y => positions%val(2)%ptr
      z => positions%val(3)%ptr

      xnonod=size(x)

      call allocate(botdis_field, xmesh,"DistanceToBottom")
      call add_boundary_condition(botdis_field,"bottom","surface", &
           (/ roadbc /))
      call insert(state,botdis_field,"DistanceToBottom")
      call deallocate(botdis_field)

      dist => extract_scalar_field(state,name)
      call get_boundary_condition(dist,1,surface_mesh=surface_mesh, &
           surface_element_list=surface_element_list)

      ! get the number of surface elements
      ntri = size(surface_element_list)
      ! get the number of nodes for each element
      snloc = face_loc(xmesh,1)

      ! surface element list
      allocate(senlist(ntri*snloc))

      do i=1, ntri
         ! get the surface element number
         sele=surface_element_list(i)
         ! get its global node numbers
         senlist((i-1)*snloc+1:i*snloc) = face_global_nodes(xmesh,sele)
      end do

      call allocate(bc_position,positions%dim,surface_mesh)
      call remap_field_to_surface(positions,bc_position,surface_element_list)
      
      allocate(eid(points), eshape(snloc, points))
      bc_position_2d = wrap_vector_field(bc_position%mesh, bc_position%val(1)%ptr, bc_position%val(2)%ptr, name = trim(bc_position%name) // "2d")
      call picker_inquire(bc_position_2d, coordsx = pointX, coordsy = pointY, &
        & eles = eid, local_coords = eshape, global = .true.)
      
      do i=1, points         
         if (eid(i) > 0)then
            
            pp1 = (eid(i)-1)*snloc+1
            
            l1  = senlist(pp1)
            l2  = senlist(pp1+1)
            l3  = senlist(pp1+2)
            
            zz = z(l1)*eshape(1, i)+ &
                 z(l2)*eshape(2, i)+ &
                 z(l3)*eshape(3, i)            
         endif
         
#ifdef HAVE_MPI
         call mpi_allreduce(zz,PointZ(i),1,GetPReal(),MPI_MAX,MPI_COMM_WORLD,ierr);
#else
         PointZ(i)=zz
#endif
      enddo

      call deallocate(bc_position_2d)
      call deallocate(bc_position)
      deallocate(eid,eshape)
      deallocate(senlist)
        
      return
    end subroutine get_floor
    
    !----------------------------------------------------------
    
    subroutine traffic_check_options()
      
      implicit none
      
      if (have_option("/traffic_model")) then
         
         if (option_count("/material_phase")/=1) then
            FLExit("We cannot model traffic underwater yet...")
         endif
         
         if (.not.have_option("/material_phase[0]/scalar_field::SolidPhase")) then
            FLExit("You haven't switched on the Volume of the vehicles")
         endif
         
         if (.not.have_option("/material_phase[0]/vector_field::Velocity/prognostic/vector_field::Source")) then
            FLExit("You haven't switched on the Source in the Velocity field")
         endif
         if (have_option("/material_phase[0]/vector_field::Velocity/prognostic/vector_field::Source/prescribed")) then
            FLExit("The Source here should always be diagnostic")
         endif
         
         if (.not.have_option("/material_phase[0]/vector_field::Velocity/prognostic/vector_field::Absorption")) then
            FLExit("You haven't switched on the Absorption in the Velocity field")
         endif
         if (have_option("/material_phase[0]/vector_field::Velocity/prognostic/vector_field::Absorption/prescribed")) then
            FLExit("The Absorption here should always be diagnostic")
         endif
         
         if (.not.have_option("/material_phase[0]/scalar_field::CopyofDensity")) then
            FLExit("You haven't switched on the CopyofDensity in the Velocity field")
         endif
         
         if (.not.have_option("/material_phase[0]/vector_field::SolidVelocity")) then
            FLExit("You haven't switched on the SolidVelocity in the Velocity field")
         endif
         
         if (have_option("/traffic_model/scalar_field::TrafficTracerTemplate")) then
            if(.not.have_option("/traffic_model/scalar_field::TrafficTracerTemplate/prognostic/equation::AdvectionDiffusion")) then
               FLExit("TrafficTracer::equation should be set to AdvectionDiffusion")
            endif
            if(.not.have_option("/traffic_model/scalar_field::TrafficTracerTemplate/prognostic/scalar_field::Source/diagnostic")) then
               FLExit("TrafficTracer::you should turn on the source and set it to be diagnostic")
            endif
         endif
         
         if (have_option("/traffic_model/emissions_factors/") &
              .and. .not.have_option("/traffic_model/scalar_field::TrafficTracerTemplate/")) then
            FLExit("You cannot use emissions factors if you haven't switched on the traffic tracers")
         endif
         
         if(.not.have_option("/mesh_adaptivity")) then
            FLExit("You should turn on adaptivity before you use this model")
         endif
         
      endif
      
      return
    end subroutine traffic_check_options
end module traffic
