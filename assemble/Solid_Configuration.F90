#include "confdefs.h"
#include "fdebug.h"

module solidconfiguration
  use fldebug
  use fields
  use fefields,only: compute_lumped_mass
  use fetools, only: X_, Y_, Z_
  use spud
  use state_module
  use transform_elements
  use boundary_conditions
  use global_parameters, only: OPTION_PATH_LEN, dt, current_time
  use write_triangle
  use parallel_tools
  use vector_tools, only: det
  use vtk_interfaces, only: vtk_write_fields
  implicit none

  private
  !for particle positions and radii ...
  integer              ,save :: nparticles,quadnodes
  real                 ,save :: cylinder_width
  real,    allocatable ,save :: centerpx(:),centerpy(:),centerpz(:)
  real,    allocatable ,save :: radii(:)

  !for particle translational and angular velocities, and forces.
  real,    allocatable ,save :: pvelx(:),pvely(:),pvelz(:)
  real,    allocatable ,save :: pavelx(:),pavely(:),pavelz(:)
  real,    allocatable ,save :: forcex(:),forcey(:),forcez(:)
  real,    allocatable ,save :: torquex(:),torquey(:),torquez(:)
  real,    allocatable ,save :: particle_mass(:)
  integer, allocatable       :: nelist(:,:),nelistlgth(:),eelist(:,:)

  !for particle tracking
  integer, allocatable ,save :: node_to_particle(:)

  !for external mesh representation of a particle
  integer              ,save :: ext_nonods,ext_totele
  real,    allocatable ,save :: ext_x(:,:),ext_y(:,:),ext_z(:,:)
  real,    allocatable ,save :: ext_u(:,:),ext_v(:,:),ext_w(:,:)
  integer, allocatable ,save :: quad2lin(:),lin2quad(:)
  real                 ,save :: vcheck_tol,acheck_tol
  integer, allocatable ,save :: bin(:,:,:),ind(:,:)

  !forces on the external mesh
  real,    allocatable, save :: ext_forcex(:,:),ext_forcey(:,:),ext_forcez(:,:)
  real,    allocatable, save :: ext_torquex(:,:),ext_torquey(:,:),ext_torquez(:,:)

  !general parameters for calculating forces
  real                     ,save :: solid_absorption_factor,solid_density,k
  real                     ,save :: solid_concentration_max
  logical                  ,save :: multimaterial,oneway,viscon
  logical                  ,save :: drag,output_drag,output_particle_vtus

  !for physical calculations
  real,                 save :: gravity_x,gravity_y,gravity_z,gravity_m,bottom
  real,    dimension(3),save :: gravity_vector
  logical, save              :: buoyancy
  real,save                  :: buoyancy_factor

  !for simple dynamics box size
  real,                 save :: xmin,ymin,zmin,xmax,ymax,zmax

  !option paths
  character(len=OPTION_PATH_LEN), save :: solid_path,dynamic_path
  character(len=OPTION_PATH_LEN), save :: solid_type,dynamic_type,mapping_type
  character(len=OPTION_PATH_LEN), save :: input_file_name,y3d_file_name,femdem3d_filename,quad2lin_filename

  !iteration counter
  integer, save:: iteration
  data iteration /0/

  !parameters
  integer nloc
  real pie
  parameter(nloc=4,pie=3.1415926535897931)

#ifdef USING_FEMDEM
  !  interface
  !     subroutine y3allocate_femdem(string)
  !       character(len=*) :: string
  !     end subroutine y3allocate_femdem
  !  end interface
  !  interface
  !     subroutine y3dfemdem(flag,timestep,posx,posy,posz,velx,vely,velz,forx,fory,forz,&
  !          linearnodes,lin2quad1,quad2lin1,string)
  !       integer,   intent(in) :: flag,linearnodes
  !       real,      intent(in) :: timestep
  !       real,    dimension(:) :: posx,posy,posz,velx,vely,velz,forx,fory,forz
  !       integer, dimension(:) :: lin2quad1,quad2lin1
  !       character(len=*)      :: string
  !     end subroutine y3dfemdem
  !  end interface
#endif

  public solid_configuration, smooth_out_solid_velocities,solid_drag_calculation,mesh2mesh_3d,&
       & drag_on_surface,calculate_particle_mass,solid_data_update,element_volume_s

contains
  subroutine solid_configuration(state,its,itinoi)
    implicit none
    !variables in the call
    type(state_type), dimension(:), intent(inout) :: state
    integer, intent(in)    :: its,itinoi  !this is just to know where in the non linear loop the subroutine is called

    !general integers
    integer  inod,ielem,node,stat,i,j
    integer  iloc,countinside,ipart

    !for tetrahedral local coordinate handling
    real xx(4),yy(4),zz(4),pvl(4),xx1,yy1,zz1,vcheck1
    real rx,ry,rz,x,y,z

    !options
    character(len=OPTION_PATH_LEN)       :: tmp_path

    !python scripts
    character(len=OPTION_PATH_LEN), save :: python_position,python_velocity,python_radius
    character(len=OPTION_PATH_LEN), save :: python_angular_velocity

    !for buffer zone
    !real profile1,profile2,profile3
    !real dist,const_profile
    !real buffer_conc,buffer

    !new variable types.
    type(scalar_field) visualize,matvolfrac,visualizesolid,particle_scalar,solid_concentration
    type(vector_field) svelocity,position1,position2,particle_vector,particle_force
    type(vector_field) velocity,velocityplotforsolids

    !logicals
    logical inside,got_particles,got_quad2lin

    !files
    integer dat_unit

    !reals
    real particle_volume,volume

    !switch
    integer, save:: start
    data start /1/

    !switches

    !Extract fields that will always be used here.
    position1=extract_vector_field(state(1),"Coordinate")
    svelocity=extract_vector_field(state(1),"SolidVelocity")
    visualizesolid=extract_scalar_field(state(1),"VisualizeSolid")
    solid_concentration=extract_scalar_field(state(1),"SolidConcentration")


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !section 1 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !read input parameters

    if (start.eq.1) then
       !Section1 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !read input parameters
       !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ewrite(0,*) '----------------------------------------------------------------'
       ewrite(0,*) 'Solid_Configuration : reading input...'
       solid_path="/imported_solids"
       call get_option(trim(solid_path)//"/solid_type/name",solid_type)
       dynamic_path=trim(solid_path)//"/dynamic_type"
       call get_option(trim(dynamic_path)//"/name",dynamic_type)
       !read common options to all solid_types.
       call get_option(trim(solid_path)//"/solid_concentration_max",&
            solid_concentration_max)
       call get_option(trim(solid_path)//"/solid_absorption_factor",&
            solid_absorption_factor)
       call get_option(trim(solid_path)//"/number_of_particles",nparticles)
       if(nparticles.le.0) then
          ewrite(-1,*) "Your options have: ", trim(solid_path)//"/number_of_particles" , "set to: ",nparticles
          FLExit('The number of particles cannot be negative or zero')
       end if

       oneway=.false.
       buoyancy=.false.
       viscon=.false.
       multimaterial=.false.
       drag=.false.
       output_drag=.false.
       output_particle_vtus=.false.

       !check if problem involves using multiple-fluids.
       if (have_option(trim(solid_path)//"/use_multimaterials")) then
          multimaterial=.true.
       end if
       !one way coupling. If true, solid is not affected by solid forces
       if (have_option(trim(solid_path)//"/oneway")) then
          oneway=.true.
       end if
       !turn on "artificial" buoyancy
       if (have_option(trim(solid_path)//"/buoyancy")) then
          buoyancy=.true.
       end if
       !visualize both fluids and solid parts in one scalar field
       if (have_option(trim(solid_path)//"/visualize_solidfluid")) then
          viscon=.true.
       end if
       if (have_option(trim(solid_path)//"/calculate_drag")) then
          drag=.true.
       end if
       if (have_option(trim(solid_path)//"/output_drag")) then
          output_drag=.true.
          call initialise_output_files()
       end if
       if (have_option(trim(solid_path)//"/output_particle_vtus")) then
          output_particle_vtus=.true.
       end if
       !read solid density
       call get_option('/imported_solids/solid_density',solid_density)

       !Read gravity info
       call get_option('/physical_parameters/gravity/magnitude',gravity_m)
       call get_option("/physical_parameters/gravity/&
            &vector_field::GravityDirection/prescribed/value[0]/constant",&
            & gravity_vector)
       gravity_x=-gravity_m*gravity_vector(1)
       gravity_y=-gravity_m*gravity_vector(2)
       gravity_z=-gravity_m*gravity_vector(3)

       !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !Allocate memory for geometry, position, and velocity information
       !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       allocate(centerpx(nparticles))
       allocate(centerpy(nparticles))
       allocate(centerpz(nparticles))
       allocate(pvelx(nparticles))
       allocate(pvely(nparticles))
       allocate(pvelz(nparticles))
       allocate(pavelx(nparticles))
       allocate(pavely(nparticles))
       allocate(pavelz(nparticles))
       allocate(particle_mass(nparticles))

       select case(trim(solid_type))
       case('spheres')
          allocate(radii(nparticles))
       case('small_spheres')
          allocate(radii(nparticles))
       case('cylinders')
          allocate(radii(nparticles))
          call get_option('/imported_solids/solid_type/cylinder_width',cylinder_width)
       case('external_2D_mesh')
          position2=extract_vector_field(state(1),"ParticleMeshCoordinate")
          allocate(ext_x(node_count(position2),nparticles))
          allocate(ext_y(node_count(position2),nparticles))
          allocate(ext_z(node_count(position2),nparticles))
          allocate(ext_u(node_count(position2),nparticles))
          allocate(ext_v(node_count(position2),nparticles))
          allocate(ext_w(node_count(position2),nparticles))
          call get_option(trim(solid_path)//"/solid_type/mapping_type/name", mapping_type)
          call get_option(trim(solid_path)//"/solid_type/volume_checking_tol", vcheck_tol)
          !construct the rest of the particles (they will share the connectivity
          do ipart=1,nparticles
             do inod=1,node_count(position2)
                ext_x(inod,ipart)=position2%val(X_,inod)
                ext_y(inod,ipart)=position2%val(Y_,inod)
                ext_z(inod,ipart)=position2%val(Z_,inod)
             end do
          end do
       case('external_3D_mesh')
          allocate(radii(nparticles))
          position2=extract_vector_field(state(1),"ParticleMeshCoordinate")
          allocate(ext_x(node_count(position2),nparticles))
          allocate(ext_y(node_count(position2),nparticles))
          allocate(ext_z(node_count(position2),nparticles))
          allocate(ext_u(node_count(position2),nparticles))
          allocate(ext_v(node_count(position2),nparticles))
          allocate(ext_w(node_count(position2),nparticles))
          call get_option(trim(solid_path)//"/solid_type/mapping_type/name", mapping_type)
          call get_option(trim(solid_path)//"/solid_type/volume_checking_tol", vcheck_tol)
          !construct the rest of the particles (they will share the connectivity
          do ipart=1,nparticles
             do inod=1,node_count(position2)
                ext_x(inod,ipart)=position2%val(X_,inod)
                ext_y(inod,ipart)=position2%val(Y_,inod)
                ext_z(inod,ipart)=position2%val(Z_,inod)
             end do
          end do
       end select

       !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !Set initial conditions for solid particles
       !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       select case(trim(dynamic_type))
       case("use_simple_dynamics")
          call get_option(trim(dynamic_path)//"/set_bottom",bottom)
          call get_option(trim(dynamic_path)//"/set_xmin",xmin)
          call get_option(trim(dynamic_path)//"/set_ymin",ymin)
          call get_option(trim(dynamic_path)//"/set_zmin",zmin)
          call get_option(trim(dynamic_path)//"/set_xmax",xmax)
          call get_option(trim(dynamic_path)//"/set_ymax",ymax)
          call get_option(trim(dynamic_path)//"/set_zmax",zmax)
          call get_option(trim(solid_path)//"/position_script",python_position)
          call get_option(trim(solid_path)//"/translation_velocity_script",python_velocity)
          call get_option(trim(solid_path)//"/angular_velocity_script",python_angular_velocity)
          if(trim(solid_type)=="cylinders".or.trim(solid_type)=="spheres".or.trim(solid_type)=="small_spheres") then
             call get_option(trim(solid_path)//"/radius_script",python_radius)
             !reading particle radius.
             call set_particle_sfield_from_python(python_radius,len(python_radius),&
                  nparticles, current_time, radii, stat)
             do ipart=1,nparticles
                ewrite(0,*)'Finihed reading radii: ',radii(ipart)
             end do
          end if
          !note: initial angular positions have not been implemented
          !external meshed particles start without any rotation
          !read in particle translational velocities and angular velocities.
          call set_particle_vfield_from_python(python_position,len(python_position), &
               nparticles, current_time, centerpx, centerpy, centerpz,stat)
          call set_particle_vfield_from_python(python_velocity,len(python_velocity), &
               nparticles, current_time, pvelx, pvely, pvelz, stat)
          call set_particle_vfield_from_python(python_angular_velocity,len(python_angular_velocity), &
               nparticles, current_time, pavelx, pavely, pavelz, stat)
          do ipart=1,nparticles
             ewrite(0,*)'Finihed reading pos: ',centerpx(ipart),centerpy(ipart),centerpz(ipart)
             ewrite(0,*)'Finihed reading pos: ',pvelx(ipart),pvely(ipart),pvelz(ipart)
          end do
          if (trim(solid_type)=='external_3D_mesh'.or.trim(solid_type)=='external_2D_mesh') then
             !position each particle...
             do ipart=1,nparticles
                do inod=1,node_count(position2)
                   ext_x(inod,ipart)= ext_x(inod,ipart) + centerpx(ipart)
                   ext_y(inod,ipart)= ext_y(inod,ipart) + centerpy(ipart)
                   ext_z(inod,ipart)= ext_z(inod,ipart) + centerpz(ipart)
                end do
             end do
          end if
       case("python_script")
          call get_option(trim(solid_path)//"/position_script",python_position)
          call get_option(trim(solid_path)//"/translation_velocity_script",python_velocity)
          call get_option(trim(solid_path)//"/angular_velocity_script",python_angular_velocity)
          if(trim(solid_type)=="cylinders".or.trim(solid_type)=="spheres".or.trim(solid_type)=="small_spheres") then
             call get_option(trim(solid_path)//"/radius_script",python_radius)
             !reading particle radius.
             call set_particle_sfield_from_python(python_radius,len(python_radius),&
                  nparticles, current_time, radii, stat)
          end if
          !note: initial angular positions have not been implemented
          !external meshed particles start without any rotation
          !read in particle translational velocities and angular velocities.
          call set_particle_vfield_from_python(python_position,len(python_position), &
               nparticles, current_time, centerpx, centerpy, centerpz,stat)
          call set_particle_vfield_from_python(python_velocity,len(python_velocity), &
               nparticles, current_time, pvelx, pvely, pvelz, stat)
          call set_particle_vfield_from_python(python_angular_velocity,len(python_angular_velocity), &
               nparticles, current_time, pavelx, pavely, pavelz, stat)
          if (trim(solid_type)=='external_3D_mesh'.or.trim(solid_type)=='external_2D_mesh') then
             !position each particle...
             do ipart=1,nparticles
                do inod=1,node_count(position2)
                   ext_x(inod,ipart)= ext_x(inod,ipart) + centerpx(ipart)
                   ext_y(inod,ipart)= ext_y(inod,ipart) + centerpy(ipart)
                   ext_z(inod,ipart)= ext_z(inod,ipart) + centerpz(ipart)
                end do
             end do
          end if
       case("from_input_file")
          tmp_path=trim(dynamic_path)//"/file_name"
          if (have_option(trim(tmp_path))) then
             call get_option(trim(tmp_path),input_file_name)
          end if
          inquire(file=trim(input_file_name),exist=got_particles)
          !Read input variables and sanity check
          if(.not.got_particles) then
             ewrite(-1,*) "Solids configuration"
             ewrite(-1,*) "Looking for file: ",trim(input_file_name)
             FLExit("Unfortunately, this file was not found.")
          end if
          dat_unit=free_unit()
          open(dat_unit,file=trim(input_file_name),status='OLD')
          !skip four lines
          do ipart=1,4
             read(dat_unit,*)
          end do
          !Reading particle reference center, radius, traslational velocity, and angular velocity
          do ipart=1,nparticles
             read(dat_unit,*) centerpx(ipart),centerpy(ipart),centerpz(ipart),&
                  radii(ipart),pvelx(ipart),pvely(ipart),pvelz(ipart),&
                  pavelx(ipart),pavely(ipart),pavelz(ipart)
          end do
          close(dat_unit)
          ewrite(0,*) "Finished reading particles from input file"
       case("use_Y3D")
       case("use_2Dfemdem")
       case("use_3Dfemdem")
#ifdef USING_FEMDEM
          !Initialize 3D femdem code
          !Read in quadratic-to-linear information
          call get_option(trim(dynamic_path)//'/quad2lin/file_name',quad2lin_filename)
          inquire(file=trim(quad2lin_filename),exist=got_quad2lin)
          !Read input variables and sanity check
          if(.not.got_quad2lin) then
             ewrite(-1,*) "Solids configuration: reading quadratic-to-linear info"
             ewrite(-1,*) "Looking for file: ",trim(input_file_name)
             FLExit("Unfortunately, this file was not found.")
          end if
          dat_unit=free_unit()
          open(dat_unit,file=trim(quad2lin_filename),status='OLD')
          allocate(lin2quad(node_count(position2)))
          read(dat_unit,*)
          read(dat_unit,*) quadnodes
          read(dat_unit,*)
          allocate(quad2lin(quadnodes))
          ewrite(0,*) quadnodes, node_count(position2)
          do inod=1,node_count(position2)
             read(dat_unit,*) lin2quad(inod)
          end do
          read(dat_unit,*)
          do inod=1,quadnodes
             read(dat_unit,*) quad2lin(inod)
          end do
          close(dat_unit)
          call get_option(trim(dynamic_path)//"/file_name",femdem3d_filename)
          ewrite(0,*) 'right before allocating y3dfemdem: ',trim(femdem3d_filename)
          !allocate 3dfemdem memory
          !call y3allocate_femdem(trim(femdem3d_filename)//char(0))
          ewrite(0,*) 'after allocating y3dfemdem'
          particle_vector=extract_vector_field(state(1),"ParticleVector")
          particle_force=extract_vector_field(state(1),"ParticleForce")
          particle_force%val(X_,:)=0.0
          particle_force%val(Y_,:)=0.0
          particle_force%val(Z_,:)=0.0
          ewrite(0,*) 'right before first call to y3dfemdem',node_count(position2),node_count(particle_vector)
          !initialize 3dfemdem code
          call y3dfemdem(1,DT,position2%val(X_,:),position2%val(Y_,:),position2%val(Z_,:), &
               particle_vector%val(X_,:),particle_vector%val(Y_,:),particle_vector%val(Z_,:), &
               particle_force%val(X_,:),particle_force%val(Y_,:),particle_force%val(Z_,:), &
               node_count(position2),lin2quad,quad2lin,&
               trim(femdem3d_filename)//char(0))
          ewrite(0,*) 'right after first call to y3dfemdem',node_count(position2),node_count(particle_vector)
          !Initializing velocity buffer
          do ipart=1,nparticles
             do inod=1,node_count(position2)
                ext_u(inod,ipart)=particle_vector%val(X_,inod)
                ext_v(inod,ipart)=particle_vector%val(Y_,inod)
                ext_w(inod,ipart)=particle_vector%val(Z_,inod)
             end do
          end do
#else
          ewrite(-1,*) "Error with your fluidity build"
          FLExit('You need to configure --with-femdem=/PATH_TO_FEMDEM3D_DIR')
#endif
       end select
       ewrite(0,*) 'finished reading solid_configuration parameters ...'
    end if !if(start.eq.1)

    !If particle motion is ruled by a python script, then recalculate its velocity
    if (trim(dynamic_type)=='python_script'.and.start/=1) then
       call set_particle_vfield_from_python(python_velocity,len(python_velocity), &
            nparticles, current_time, pvelx, pvely, pvelz, stat)
       call set_particle_vfield_from_python(python_angular_velocity,len(python_angular_velocity), &
            nparticles, current_time, pavelx, pavely, pavelz, stat)
    end if

    !If particle object is a 3D mesh, and not using 3Dfemdem, then new particle velocities
    !(which were calculated on the previous timestep) must be passed on to each particle mesh.
    if (trim(solid_type)=='external_3D_mesh') then
       if (trim(dynamic_type)=='from_input_file'.or. &
            trim(dynamic_type)=='python_script'.or. &
            trim(dynamic_type)=='use_simple_dynamics') then
          do ipart=1,nparticles
             do inod=1,node_count(position2)
                rx=ext_x(inod,ipart)-centerpx(ipart)
                ry=ext_y(inod,ipart)-centerpy(ipart)
                rz=ext_z(inod,ipart)-centerpz(ipart)
                ext_u(inod,ipart)= pvelx(ipart) - pavelz(ipart)*ry + pavely(ipart)*rz
                ext_v(inod,ipart)= pvely(ipart) - pavelx(ipart)*rz + pavelz(ipart)*rx
                ext_w(inod,ipart)= pvelz(ipart) - pavely(ipart)*rx + pavelx(ipart)*ry
             end do
          end do
       end if
    end if

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !section 2 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !mapping

    !start counters
    countinside=0

    !allocate space for particle tracking (node to particle list)
    allocate(node_to_particle(node_count(position1)))
    node_to_particle=0

    !Ensure that vectors are zero before mapping.
    solid_concentration%val=0.0
    visualizesolid%val=0.0
    do i=1,3
       svelocity%val(i,:)=0.0
    end do

    select case(trim(solid_type))
    case("cylinders")
       do inod=1,node_count(position1)
          x=position1%val(X_,inod)
          y=position1%val(Y_,inod)
          z=position1%val(Z_,inod)
          inside=.false.
          do ipart=1,nparticles
             !note: a cylinder is mapped with boundcheckc and a sphere with boundcheck
             call boundcheckc(x,y,radii(ipart),centerpx(ipart),centerpy(ipart),inside)
             if (inside) then
                countinside=countinside+1
                solid_concentration%val(inod)=1.0
                visualizesolid%val(inod)=1.0
                node_to_particle(inod)=ipart
                exit
             end if
          end do
          !given that the z coordinate is not needed. svelocity%val(Z_,:) doesn't need to be
          !calculated.
          if(inside) then
             rx=x-centerpx(ipart)
             ry=y-centerpy(ipart)
             svelocity%val(X_,inod)= pvelx(ipart) - pavelz(ipart)*ry
             svelocity%val(Y_,inod)= pvely(ipart) + pavelz(ipart)*rx
          end if
       end do
    case("spheres")
       do inod=1,node_count(position1)
          x=position1%val(X_,inod)
          y=position1%val(Y_,inod)
          z=position1%val(Z_,inod)
          inside=.false.
          do ipart=1,nparticles
             call boundcheck(x,y,z,radii(ipart),centerpx(ipart),centerpy(ipart),centerpz(ipart),inside)
             if (inside) then
                countinside=countinside+1
                solid_concentration%val(inod)=1.0
                visualizesolid%val(inod)=1.0
                node_to_particle(inod)=ipart
                exit
             end if
          end do
          if (inside) then
             rx=x-centerpx(ipart)
             ry=y-centerpy(ipart)
             rz=z-centerpz(ipart)
             svelocity%val(X_,inod)= pvelx(ipart) - pavelz(ipart)*ry + pavely(ipart)*rz
             svelocity%val(Y_,inod)= pvely(ipart) - pavelx(ipart)*rz + pavelz(ipart)*rx
             svelocity%val(Z_,inod)= pvelz(ipart) - pavely(ipart)*rx + pavelx(ipart)*ry
          end if
       end do
       visualizesolid%val=solid_concentration%val
    case("external_2D_mesh")
       ewrite(0,*) 'no mapping here yet for 2d mesh yet'
       FLExit('no mapping posibility for 2d meshes yet')
    case("external_3D_mesh")
       ewrite(0,*) 'Mapping external mesh with method: ',trim(mapping_type)
       position2=extract_vector_field(state(1),"ParticleMeshCoordinate")
       particle_vector=extract_vector_field(state(1),"ParticleVector")
       particle_scalar=extract_scalar_field(state(1),"ParticleScalar")
       particle_scalar%val=1.0
       do ipart=1,nparticles
          do inod=1,node_count(position2)
             position2%val(X_,inod)=ext_x(inod,ipart)
             position2%val(Y_,inod)=ext_y(inod,ipart)
             position2%val(Z_,inod)=ext_z(inod,ipart)
             particle_vector%val(X_,inod)=ext_u(inod,ipart)
             particle_vector%val(Y_,inod)=ext_v(inod,ipart)
             particle_vector%val(Z_,inod)=ext_w(inod,ipart)
          end do
          call  mesh2mesh_3d(state,ipart,node_to_particle, &
               "SolidConcentration","ParticleScalar",  &
               "SolidVelocity","ParticleVector")
       end do
       visualizesolid%val=solid_concentration%val
    case("small_spheres")
       !first check which element each particle lands. Brute force for now.
       ewrite(0,*) 'Mapping small spheres'
       do ipart=1,nparticles
          do ielem=1,element_count(position1)
             do iloc=1,nloc
                node = position1%mesh%ndglno((ielem-1)*nloc+iloc)
                xx(iloc) = position1%val(X_,node)
                yy(iloc) = position1%val(Y_,node)!original coordinates are stored in a 4 size vector.
                zz(iloc) = position1%val(Z_,node)
             end do
             volume=element_volume_s(nloc,xx,yy,zz)
             do iloc=1,nloc
                xx1=xx(iloc)
                yy1=yy(iloc)
                zz1=zz(iloc)
                xx(iloc)=centerpx(ipart)
                yy(iloc)=centerpy(ipart)
                zz(iloc)=centerpz(ipart)
                pvl(iloc)=element_volume_s(nloc,xx,yy,zz)
                xx(iloc)=xx1
                yy(iloc)=yy1
                zz(iloc)=zz1
             end do
             vcheck1=1.0-(abs(pvl(1))+abs(pvl(2))  &
                  +abs(pvl(3))+abs(pvl(4)))/abs(volume)
             if (abs(vcheck1).lt.1.e-8) then
                !calculate particle volume first
                particle_volume=4*pie*radii(ipart)**3.0/3.0
                !send concentration to the nodes
                do iloc=1,nloc
                   node = position1%mesh%ndglno((ielem-1)*nloc+iloc)
                   countinside=countinside+1
                   solid_concentration%val(node)=solid_concentration%val(node)+ abs(pvl(iloc)/volume)*particle_volume
                   rx=position1%val(1,node)-centerpx(ipart)
                   ry=position1%val(2,node)-centerpy(ipart)
                   rz=position1%val(3,node)-centerpz(ipart)
                   svelocity%val(X_,node)= pvelx(ipart) - pavelz(ipart)*ry + pavely(ipart)*rz
                   svelocity%val(Y_,node)= pvely(ipart) - pavelx(ipart)*rz + pavelz(ipart)*rx
                   svelocity%val(Z_,node)= pvelz(ipart) - pavely(ipart)*rx + pavelx(ipart)*ry
                   visualizesolid%val(node)=visualizesolid%val(node)+ abs(pvl(iloc)/volume)*particle_volume
                   node_to_particle(node)=ipart
                end do
                exit
             end if
          end do
       end do
       visualizesolid%val=solid_concentration%val
    end select

    !Calculate total mass of solid
    call calculate_particle_mass(state)

    !Make velocity field that does not include velocity inside the solid
    velocity=extract_vector_field(state(1),"Velocity")
    velocityplotforsolids=extract_vector_field(state(1),"VelocityPlotForSolids")
    do i=1,3
       velocityplotforsolids%val(i,:)=(1.-visualizesolid%val)*velocity%val(i,:)
    end do

    !if using multimaterials, then this option helps visualize both solid concentration and
    if (viscon) then
       visualize=extract_scalar_field(state(1),"VisualizeSolidFluid")
       matvolfrac=extract_scalar_field(state(1),"MaterialVolumeFraction")
       visualize%val=matvolfrac%val + visualizesolid%val
    end if
    start=0 !this switches of the parameter-reading side of the subroutine.
    !deallocate(node_to_particle)
    ewrite(0,*) 'Finished setting up Solid fields...'
  end subroutine solid_configuration

  subroutine initialise_output_files()
    implicit none
    integer dat_unit,myrank
    myrank=GetRank()
    !Initialize files for output to external files.
    !These files contain the total drag for each individual particle.
    if(myrank==0) then
       dat_unit=free_unit()
       open(dat_unit,file='particle_position',status='replace')
       close(dat_unit)
       open(dat_unit,file='particle_velocity',status='replace')
       close(dat_unit)
       open(dat_unit,file='particle_angular_velocity',status='replace')
       close(dat_unit)
       open(dat_unit,file='particle_force',status='replace')
       close(dat_unit)
       open(dat_unit,file='particle_torque',status='replace')
       close(dat_unit)
    end if
  end subroutine initialise_output_files

  subroutine smooth_out_solid_velocities(state,itinoi)
    implicit none
    type(state_type), dimension(:), intent(inout) :: state
    integer, intent(in) :: itinoi
    type(scalar_field) visualizesolid
    type(vector_field) svelocity,flow_velocity
    integer i

    visualizesolid=extract_scalar_field(state(1),"VisualizeSolid")
    svelocity=extract_vector_field(state(1),"SolidVelocity")
    flow_velocity=extract_vector_field(state(1),"NonlinearVelocity")

    !could also try using a buffer, and using solid_concentration (which goes up to 0.999), instead.
    do i=1,3
       flow_velocity%val(i,:)=(1.0-visualizesolid%val)*flow_velocity%val(i,:) + &
            visualizesolid%val*svelocity%val(i,:)
    end do
  end subroutine smooth_out_solid_velocities

  subroutine solid_drag_calculation(state,its,itinoi)
    implicit none
    type(state_type), dimension(:), intent(inout) :: state
    integer, intent(in) ::  its,itinoi
    integer stat,inod,i,ipart

    !new variable types.
    type(scalar_field) density,copyofdensity,visualize,matvolfrac,solid_concentration,particle_scalar
    type(vector_field) svelocity,absorption,velocity_source,flow_velocity,solid_force_vector
    type(vector_field) position2,particle_vector,particle_force
    !dummy
    integer, allocatable          :: dummy(:)
    real a,cdx,cdy,cdz,slipvelx,slipvely,slipvelz

    density=extract_scalar_field(state(1),"Density")
    absorption=extract_vector_field(state(1),"VelocityAbsorption")
    copyofdensity=extract_scalar_field(state(1),"CopyofDensity")
    velocity_source=extract_vector_field(state(1),"VelocitySource")
    flow_velocity=extract_vector_field(state(1),"Velocity")
    solid_force_vector=extract_vector_field(state(1),"SolidForce")
    svelocity=extract_vector_field(state(1),"SolidVelocity")
    solid_concentration=extract_scalar_field(state(1),"SolidConcentration")

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !Section 3 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !Absorption and source term calculation
    !just to make name shorter.
    k=solid_absorption_factor

    if(multimaterial) then
       do i=1,3
          !note: when using also multi materials, density is probably not worked out yet at this stage
          !try to work this out by averaging the fluid densities? This could avoid problems
          absorption%val(i,:)=solid_concentration%val*solid_density*k/dt
          velocity_source%val(i,:)=(1.0-solid_concentration%val)*velocity_source%val(i,:) + &
               solid_concentration%val*absorption%val(i,:)*svelocity%val(i,:)
       end do
    else
       !calculate force coefficients and source terms in each direction
       if (solid_type=='small_spheres') then
          !first the friction coefficients
          do i=1,node_count(solid_concentration)
             if (node_to_particle(i)/=0) then
                a=solid_concentration%val(i)
                ipart=node_to_particle(i)
!                slipvelx=abs(flow_velocity%val(1,i)-pvelx(ipart))
!                slipvely=abs(flow_velocity%val(2,i)-pvely(ipart))
!                slipvelz=abs(flow_velocity%val(3,i)-pvelz(ipart))

                slipvelx=abs(pvelx(ipart))
                slipvely=abs(pvely(ipart))
                slipvelz=abs(pvelz(ipart))
                cdx=drag_coefficient(a,slipvelx,density%val(i),2*radii(ipart),0.01)
                cdy=drag_coefficient(a,slipvely,density%val(i),2*radii(ipart),0.01)
                cdz=drag_coefficient(a,slipvelz,density%val(i),2*radii(ipart),0.01)
                if (slipvelx==0.0) cdx=0.0
                if (slipvely==0.0) cdy=0.0
                if (slipvelz==0.0) cdz=0.0
!                ewrite(0,*) 'data:  ',i,ipart,slipvelx,slipvely,slipvelz,a
!                ewrite(0,*) 'data1: ',radii(ipart),pvelx(ipart),pvely(ipart),pvelz(ipart)
!                ewrite(0,*) 'data2: ',cdx,cdy,cdz

                if(a>0.2) then
                   !Use ergun equation if solid volume fraction is above 0.2
                   absorption%val(1,i)=150.0*0.001*a**2.0/((1-a)*2*radii(ipart))**2.0 &
                        + 1.75*a*density%val(i)*slipvelx/(2*radii(ipart))
                   absorption%val(2,i)=150.0*0.001*a**2.0/((1-a)*2*radii(ipart))**2.0 &
                        + 1.75*a*density%val(i)*slipvely/(2*radii(ipart))
                   absorption%val(3,i)=150.0*0.001*a**2.0/((1-a)*2*radii(ipart))**2.0 &
                        + 1.75*a*density%val(i)*slipvelz/(2*radii(ipart))
!                   ewrite(0,*) 'Ergun : ',absorption%val(1,i),absorption%val(2,i),absorption%val(3,i)
                elseif (a<=0.2) then
                   !Use Wen and Yu's method if solid volume fraction is below 0.2
                   absorption%val(1,i)=0.75*cdx*slipvelx*density%val(i)*a/(2*radii(ipart)*(1-a)**2.65)
                   absorption%val(2,i)=0.75*cdy*slipvely*density%val(i)**2.0*a/(2*radii(ipart)*(1-a)**2.65)
                   absorption%val(3,i)=0.75*cdz*slipvelz*density%val(i)*a/(2*radii(ipart)*(1-a)**2.65)
!                   ewrite(0,*) 'Wen : ',absorption%val(1,i),absorption%val(2,i),absorption%val(3,i)
                endif
             end if
          end do
          !Now the source term.
          do i=1,3
             velocity_source%val(i,:)=(1.0-solid_concentration%val)*velocity_source%val(i,:) + &
                  absorption%val(i,:)*svelocity%val(i,:)
          end do
       else
          do i=1,3
             absorption%val(i,:)=solid_concentration%val*density%val*k/dt
             velocity_source%val(i,:)=(1.0-solid_concentration%val)*velocity_source%val(i,:) + &
                  solid_concentration%val*absorption%val(i,:)*svelocity%val(i,:)
          end do
       end if
    end if

    !calculate SolidForce fields
    call calculate_solid_force_vector(state)

    !tweak volume fraction (to prevent singularities)
    do inod=1,node_count(density)
       solid_concentration%val(inod) = min(solid_concentration_max, &
            solid_concentration%val(inod))
    end do

    !make copy of density and alter it to include solid_density.
    copyofdensity%val=density%val
    density%val=(1.0-solid_concentration%val)*density%val + solid_concentration%val*solid_density

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !section 4 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !calculating drag for each particle, depending on the particle type.
    allocate(forcex(nparticles))
    allocate(forcey(nparticles))
    allocate(forcez(nparticles))
    allocate(torquex(nparticles))
    allocate(torquey(nparticles))
    allocate(torquez(nparticles))
    !initialize forces on solids
    forcex=  0.0
    forcey=  0.0
    forcez=  0.0
    torquex = 0.0
    torquey = 0.0
    torquez = 0.0
    if(its==1) then
       ewrite(0,*) 'oneway :',oneway,' output drag forces:',output_drag
    end if

    !Calculate individual particle drag (This drag will vary with non-linear iterations)
    do ipart=1,nparticles
       call calculate_particle_drag(state,ipart)
       if (its==itinoi) then
          ewrite(0,'(A21,I0,A2,3F15.5)') ' forces for particle ',ipart,' :',forcex(ipart),forcey(ipart),forcez(ipart)
       end if
    end do

    !output individual particle drag
    if(its.eq.itinoi.and.output_drag) then
       call output_particle_info(nparticles,       &
            centerpx,centerpy,centerpz,            &
            pvelx,pvely,pvelz,                     &
            pavelx,pavely,pavelz,                  &
            forcex,forcey,forcez,                  &
            torquex,torquey,torquez)
    end if

    !If using external meshes, map the calculated solid force vector
    !back to them.
    if ((trim(solid_type)=='external_3D_mesh')) then
       !allocate memory for 3d external mesh forces.
       position2=extract_vector_field(state(1),"ParticleMeshCoordinate")
       allocate(ext_forcex(node_count(position2),nparticles))
       allocate(ext_forcey(node_count(position2),nparticles))
       allocate(ext_forcez(node_count(position2),nparticles))
       allocate(ext_torquex(node_count(position2),nparticles))
       allocate(ext_torquey(node_count(position2),nparticles))
       allocate(ext_torquez(node_count(position2),nparticles))
       !initialize forces on solids
       ext_forcex = 0.0
       ext_forcey = 0.0
       ext_forcez = 0.0
       ext_torquex = 0.0
       ext_torquey = 0.0
       ext_torquez = 0.0
       particle_force=extract_vector_field(state(1),"ParticleForce")
       particle_vector=extract_vector_field(state(1),"ParticleVector")
       particle_scalar=extract_scalar_field(state(1),"ParticleScalar")
       if (.not.oneway) then
          ewrite(0,*) 'Starting mapback of forces with method: ',trim(mapping_type)
          allocate(dummy(node_count(position2)))
          !Mapping Forces to external mesh
          do ipart=1,nparticles
             !put the coordinates for the corresponding particle into position2
             do inod=1,node_count(position2)
                position2%val(X_,inod)=ext_x(inod,ipart)
                position2%val(Y_,inod)=ext_y(inod,ipart)
                position2%val(Z_,inod)=ext_z(inod,ipart)
             end do
             !interpolate SolidForce values to the ParticleVector
             !(NOTE: dummy is used to simulate node_to_particle but is actually not used for anything at all)
             call  mesh2mesh_3d(state,ipart,dummy,       &
                  "ParticleScalar","SolidConcentration",  &
                  "ParticleForce","SolidForce")
             !save particle_force values in the respective external mesh force vector
             do inod=1,node_count(position2)
                ext_forcex(inod,ipart)=particle_force%val(X_,inod)
                ext_forcey(inod,ipart)=particle_force%val(Y_,inod)
                ext_forcez(inod,ipart)=particle_force%val(Z_,inod)
             end do
          end do
          ewrite(0,*) 'Finished mapback of forces.'
          deallocate(dummy)
       else
          ewrite(0,*) 'No mapback will be done. Particle forces are zero.'
          particle_force%val(X_,:)=0.0
          particle_force%val(Y_,:)=0.0
          particle_force%val(Z_,:)=0.0
       end if
       if (its==itinoi.and.output_particle_vtus) then
          call vtk_write_fields("particle", iteration, position2, position2%mesh,sfields=(/particle_scalar/),&
               &vfields=(/particle_force/))
       end if
    end if
  end subroutine solid_drag_calculation

  subroutine solid_data_update(state,its,itinoi)
    implicit none
    type(state_type), dimension(:), intent(inout) :: state
    integer, intent(in) :: its,itinoi
    type(scalar_field) :: copyofdensity,density
    type(vector_field) :: particle_vector,particle_force,position2
    integer inod,ipart

    !first recover saved density, and smooth out velocities if using non-linear iterations
    copyofdensity=extract_scalar_field(state(1),"CopyofDensity")
    density=extract_scalar_field(state(1),"Density")
    if(itinoi>1) then
       call smooth_out_solid_velocities(state,itinoi)
    end if
    density%val=copyofdensity%val

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !section 5 cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !updating particle position and velocity so they can be mapped the next timestep.

    if (its==itinoi) then
       ewrite(0,*) 'Dynamic type: ',trim(dynamic_type)
       select case (trim(dynamic_type))
       case('python_script')
          if (trim(solid_type)=='cylinders'.or. &
               trim(solid_type)=='spheres') then
             do ipart=1,nparticles
                centerpx(ipart)=centerpx(ipart) + pvelx(ipart)*dt
                centerpy(ipart)=centerpy(ipart) + pvely(ipart)*dt
                centerpz(ipart)=centerpz(ipart) + pvelz(ipart)*dt
             end do
          elseif (trim(solid_type)=='external_3D_mesh') then
             ewrite(0,*) 'updating movement of external mesh'
             position2=extract_vector_field(state(1),"ParticleMeshCoordinate")
             do ipart=1,nparticles
                do inod=1,node_count(position2)
                   ext_x(inod,ipart)= ext_x(inod,ipart) + ext_u(inod,ipart)*dt
                   ext_y(inod,ipart)= ext_y(inod,ipart) + ext_v(inod,ipart)*dt
                   ext_z(inod,ipart)= ext_z(inod,ipart) + ext_w(inod,ipart)*dt
                end do
             end do
             ewrite(0,*) 'Finished updating movement of external mesh'
          end if
       case('from_input_file')
          if (trim(solid_type)=='cylinders'.or. &
               trim(solid_type)=='spheres') then
             do ipart=1,nparticles
                centerpx(ipart)=centerpx(ipart) + pvelx(ipart)*dt
                centerpy(ipart)=centerpy(ipart) + pvely(ipart)*dt
                centerpz(ipart)=centerpz(ipart) + pvelz(ipart)*dt
             end do
          elseif (trim(solid_type)=='external_3D_mesh') then
             ewrite(0,*) 'updating movement of external mesh'
             position2=extract_vector_field(state(1),"ParticleMeshCoordinate")
             do ipart=1,nparticles
                do inod=1,node_count(position2)
                   ext_x(inod,ipart)= ext_x(inod,ipart) + ext_u(inod,ipart)*dt
                   ext_y(inod,ipart)= ext_y(inod,ipart) + ext_v(inod,ipart)*dt
                   ext_z(inod,ipart)= ext_z(inod,ipart) + ext_w(inod,ipart)*dt
                end do
             end do
             ewrite(0,*) 'Finished updating movement of external mesh'
          end if
       case('use_simple_dynamics')
          call hookforce(state)
       case('use_y3D')
          FLExit('DEM not implemented yet')
       case('use_3Dfemdem')
#ifdef USING_FEMDEM
          position2=extract_vector_field(state(1),"ParticleMeshCoordinate")
          particle_vector=extract_vector_field(state(1),"ParticleVector")
          particle_force=extract_vector_field(state(1),"ParticleForce")
          ewrite(0,*) 'right before second call to y3dfemdem'
          !inputs are in particle_force, while outputs are the new nodal
          !positions and velocities (position2, and particle_vector,
          !respectively)
          call y3dfemdem(0,DT,position2%val(X_,:),position2%val(Y_,:),position2%val(Z_,:), &
               particle_vector%val(X_,:),particle_vector%val(Y_,:),particle_vector%val(Z_,:), &
               particle_force%val(X_,:),particle_force%val(Y_,:),particle_force%val(Z_,:), &
               node_count(position2),lin2quad,quad2lin,&
               trim(femdem3d_filename)//char(0))
          ewrite(0,*) 'right after second call to y3dfemdem'
          !updating position and velocity buffers for the next timestep.
          do ipart=1,nparticles
             do inod=1,node_count(position2)
                ext_x(inod,ipart)=position2%val(X_,inod)
                ext_y(inod,ipart)=position2%val(Y_,inod)
                ext_z(inod,ipart)=position2%val(Z_,inod)
                ext_u(inod,ipart)=particle_vector%val(X_,inod)
                ext_v(inod,ipart)=particle_vector%val(Y_,inod)
                ext_w(inod,ipart)=particle_vector%val(Z_,inod)
             end do
          end do
#else
          ewrite(-1,*) "Error with your fluidity build"
          FLExit('You need to configure --with-femdem=/PATH_TO_FEMDEM3D_DIR')
#endif
       end select
    end if

    !forces may be deallocated as they are re-calculated every timestep
    deallocate(forcex)
    deallocate(forcey)
    deallocate(forcez)
    deallocate(torquex)
    deallocate(torquey)
    deallocate(torquez)
    ewrite(0,*) 'Deallocating'
    if ((trim(solid_type)=='external_2D_mesh').or. &
         (trim(solid_type)=='external_3D_mesh')) then
       deallocate(ext_forcex)
       deallocate(ext_forcey)
       deallocate(ext_forcez)
       deallocate(ext_torquex)
       deallocate(ext_torquey)
       deallocate(ext_torquez)
    end if
    deallocate(node_to_particle)
    if (its==itinoi) then
       iteration=iteration+1
    end if

  end subroutine solid_data_update

  subroutine boundcheck(x,y,z,bradius,x0,y0,z0,inside)
    implicit none
    real, intent(in) :: x,y,z,bradius,x0,y0,z0
    logical, intent(inout) :: inside
    !local variables
    real rr
    rr=((x-x0)**2.0+(y-y0)**2.0+(z-z0)**2.0)**0.5
    !note:  bounding radius should be bigger than the radius of the sphere!!
    if(rr.le.bradius) then
       inside=.true.
    else
       inside=.false.
    end if
  end subroutine boundcheck

  subroutine boundcheckc(x,y,radius,x0,y0,inside)
    implicit none
    real, intent(in) :: x,y,radius,x0,y0
    logical, intent(inout) :: inside
    !local variables
    real rr
    rr=((x-x0)**2.0+(y-y0)**2.0)**0.5
    !note:  bounding radius should be bigger than the radius of the sphere!!
    if(rr.le.radius) then
       inside=.true.
    else
       inside=.false.
    end if
  end subroutine boundcheckc

  !subroutine that calculates the volume of a tetrahedral element.
  !notice difference in convention with fluidity leads
  !opposite signs in the results.
  function element_volume_s(nloc,xx,yy,zz) result(volume)

    implicit none
    integer, intent(in) :: nloc
    real,    intent(in) :: xx(nloc),yy(nloc),zz(nloc)
    real :: volume
    integer i,ii,node(nloc)
    real matrix(3,3),a,b,c,d

    node(1)=2
    node(2)=3
    node(3)=4

    do i=1,3
       matrix(i,1)=xx(node(i))
       matrix(i,2)=yy(node(i))
       matrix(i,3)=zz(node(i))
    end do
    a=det(matrix)

    do i=1,3
       matrix(i,1)=1.0
    end do
    b=-det(matrix)

    do i=1,3
       matrix(i,2)=xx(node(i))
    end do
    c=det(matrix)

    do i=1,3
       matrix(i,2)=xx(node(i))
       matrix(i,3)=yy(node(i))
    end do
    d=-det(matrix)

    volume=(a+b*xx(1)+c*yy(1)+d*zz(1))/6.0
    return
  end function element_volume_s

  function element_derivatives(nloc,xx,yy,zz) result(dn)
    implicit none
    integer, intent(in) :: nloc
    real,    intent(in) :: xx(nloc),yy(nloc),zz(nloc)
    real :: dn(3)
    real :: volume6
    integer i,ii,node(nloc)
    real matrix(3,3),a,b,c,d

    node(1)=2
    node(2)=3
    node(3)=4

    do i=1,3
       matrix(i,1)=xx(node(i))
       matrix(i,2)=yy(node(i))
       matrix(i,3)=zz(node(i))
    end do
    a=det(matrix)

    do i=1,3
       matrix(i,1)=1.0
    end do
    b=-det(matrix)

    do i=1,3
       matrix(i,2)=xx(node(i))
    end do
    c=det(matrix)

    do i=1,3
       matrix(i,2)=xx(node(i))
       matrix(i,3)=yy(node(i))
    end do
    d=-det(matrix)

    volume6=(a+b*xx(1)+c*yy(1)+d*zz(1))
    dn(1)=b/volume6
    dn(2)=c/volume6
    dn(3)=d/volume6
    return
  end function element_derivatives

  subroutine mesh2mesh_3d(state,ipart,node_to_particle1, &
       field1name,field2name,vfield1name,vfield2name)
    !this subroutine maps values of field2(mesh2) to the mesh from field1 (mesh1)
    !obviously, field1 and field2 are on different meshes.
    implicit none
    type(state_type), intent(inout),dimension(:) :: state
    character(len=*) ,intent(in)          :: field1name,field2name
    character(len=*) ,intent(in)          :: vfield1name,vfield2name
    integer, intent(in)                   :: ipart
    integer, intent(inout),dimension(:)   :: node_to_particle1

    type(scalar_field) :: field1,field2
    type(vector_field) :: vfield1,vfield2
    type(vector_field) :: position1,position2
    character(len=OPTION_PATH_LEN) :: position_name1,position_name2

    integer nonods1,nonods2,totele1,totele2
    integer number_of_bins_x,number_of_bins_y,number_of_bins_z
    integer nfastest,nbrute,ntry,inod,ielem,iele,i,j,k,imin,new_i,new_j,new_k
    real ffmin

    real bin_x_max,bin_y_max,bin_z_max
    real bin_x_min,bin_y_min,bin_z_min,dx,dy,dz
    real xx(4),yy(4),zz(4),xx1,yy1,zz1,pvl(4),volume
    real sec1,sec2,shapef,vcheck1

    integer max_number_of_bins_x,max_number_of_bins_y,max_number_of_bins_z
    integer ii,iloc,nloc,itry,node

    logical inside
    call cpu_time(sec1)
    nloc=4

    field1=extract_scalar_field(state(1),trim(field1name))
    field2=extract_scalar_field(state(1),trim(field2name))
    vfield1=extract_vector_field(state(1),trim(vfield1name))
    vfield2=extract_vector_field(state(1),trim(vfield2name))

    !Obtain the respective mesh names and Coordinate names to extract nodal positions
    if (field1%mesh%name/="CoordinateMesh") then
       position_name1=trim(field1%mesh%name)//"Coordinate"
    else
       position_name1="Coordinate"
    end if
    if (field2%mesh%name/="CoordinateMesh") then
       position_name2=trim(field2%mesh%name)//"Coordinate"
    else
       position_name2="Coordinate"
    end if

    !mapping method derived from lohner's book (fastest neighbor to neighbor search)
    !with fallback on brute force
    call get_option("/imported_solids/solid_type/volume_checking_tol", vcheck_tol)
    if (trim(mapping_type)=='fastest_n_to_n') then
       call get_option("/imported_solids/solid_type/mapping_type/ntry", ntry)
    end  if

    nonods1=node_count(field1)
    nonods2=node_count(field2)
    totele1=element_count(field1)
    totele2=element_count(field2)

    !generate node to element list for mesh 2
    call generate_nelist(field2%mesh)
    !generate element to element list for mesh 2
    call generate_eelist(field2%mesh)

    call get_bin_size(state,position_name1,position_name2,ntry, &
         number_of_bins_x,number_of_bins_y,number_of_bins_z,          &
         bin_x_min,bin_y_min,bin_z_min,                               &
         bin_x_max,bin_y_max,bin_z_max,                               &
         dx,dy,dz)

    call construct_bin(                                      &
         state,                                              &
         number_of_bins_x,number_of_bins_y,number_of_bins_z, &
         bin_x_max,bin_y_max,bin_z_max,                      &
         bin_x_min,bin_y_min,bin_z_min,                      &
         position_name2,dx,dy,dz)

    !-------------------------------------------------------------------
    !performing fastest neighbor to neighbor search.
    nfastest=0
    nbrute=0
    ntry=20

    position1=extract_vector_field(state(1),trim(position_name1))
    position2=extract_vector_field(state(1),trim(position_name2))

    !search the field1 nodes: loop over the nodes skipping those outside of the bin's range
    do inod=1,nonods1
       inside=.false.
       call check_if_in_bin(inside,bin_x_max,bin_y_max,bin_z_max, &
            bin_x_min,bin_y_min,bin_z_min,                        &
            position1%val(X_,inod),position1%val(Y_,inod),position1%val(Z_,inod))
       if (inside) then
          !find which bin the node "inod" lands in
          i=min(max(int(abs((position1%val(X_,inod)-bin_x_min)/dx))+1,1),number_of_bins_x)
          j=min(max(int(abs((position1%val(Y_,inod)-bin_y_min)/dy))+1,1),number_of_bins_y)
          k=min(max(int(abs((position1%val(Z_,inod)-bin_z_min)/dz))+1,1),number_of_bins_z)

          !search which fluid2 mesh element associated with the bin where the inod lands.
          ielem=bin(i,j,k)
          if(ielem==-1) then
             !if no element is associated with the bin, search surrounding bins for a non -1 value of id
             do ii=1,27
                new_i = min(max(i + ind(ii,1),1),number_of_bins_x)
                new_j = min(max(j + ind(ii,2),1),number_of_bins_y)
                new_k = min(max(k + ind(ii,3),1),number_of_bins_z)
                ielem=bin(new_i,new_j,new_k)
                if(ielem/=-1) exit
             end do
          end if

          if (ielem/=-1) then
             !once an element is found, begin the neighbor element jumps.
             do itry=1,ntry
                do iloc=1,nloc
                   node = field2%mesh%ndglno((ielem-1)*nloc+iloc)
                   xx(iloc) = position2%val(X_,node)
                   yy(iloc) = position2%val(Y_,node) !original coordinates are stored in a 4 size vector.
                   zz(iloc) = position2%val(Z_,node)
                end do

                volume=element_volume_s(nloc,xx,yy,zz)

                ffmin=huge(0.0)
                do iloc=1,nloc
                   xx1=xx(iloc)
                   yy1=yy(iloc)
                   zz1=zz(iloc)
                   xx(iloc)=position1%val(X_,inod)
                   yy(iloc)=position1%val(Y_,inod)
                   zz(iloc)=position1%val(Z_,inod)
                   pvl(iloc)=element_volume_s(nloc,xx,yy,zz)
                   xx(iloc)=xx1
                   yy(iloc)=yy1
                   zz(iloc)=zz1
                   !calculate minimum shape function value.
                   if(pvl(iloc)/volume.lt.ffmin) then
                      imin=iloc
                      ffmin=pvl(iloc)/volume
                   end if
                end do

                vcheck1=1.0-(abs(pvl(1))+abs(pvl(2))+abs(pvl(3))+abs(pvl(4)))/abs(volume)
                if (abs(vcheck1).lt.vcheck_tol) then
                   field1%val(inod)=0.0
                   vfield1%val(X_,inod)=0.0
                   vfield1%val(Y_,inod)=0.0
                   vfield1%val(Z_,inod)=0.0
                   do iloc=1,4
                      node=field2%mesh%ndglno((ielem-1)*4+iloc)
                      shapef=pvl(iloc)/volume
                      field1%val(inod)=field1%val(inod) + field2%val(node)*shapef
                      vfield1%val(1,inod)= vfield1%val(1,inod) + &
                           vfield2%val(1,node)*shapef
                      vfield1%val(2,inod)= vfield1%val(2,inod) + &
                           vfield2%val(2,node)*shapef
                      vfield1%val(3,inod)= vfield1%val(3,inod) + &
                           vfield2%val(3,node)*shapef
                   end do
                   node_to_particle1(inod)=ipart
                   nfastest=nfastest+1
                   exit
                else
                   ielem=eelist(ielem,imin)!new node based on minimum value of form function
                   if(ielem.eq.-1) exit
                end if
             end do
          end if

          !if the search does not succeed, go on with brute force search.
          if (node_to_particle1(inod)==0) then
             do ielem=1,totele2
                do iloc=1,nloc
                   node = field2%mesh%ndglno((ielem-1)*nloc+iloc)
                   xx(iloc) = position2%val(X_,node)
                   yy(iloc) = position2%val(Y_,node)!original coordinates are stored in a 4 size vector.
                   zz(iloc) = position2%val(Z_,node)
                end do
                volume=element_volume_s(nloc,xx,yy,zz)
                do iloc=1,nloc
                   xx1=xx(iloc)
                   yy1=yy(iloc)
                   zz1=zz(iloc)
                   xx(iloc)=position1%val(X_,inod)
                   yy(iloc)=position1%val(Y_,inod)
                   zz(iloc)=position1%val(Z_,inod)
                   pvl(iloc)=element_volume_s(nloc,xx,yy,zz)
                   xx(iloc)=xx1
                   yy(iloc)=yy1
                   zz(iloc)=zz1
                end do
                vcheck1=1.0-(abs(pvl(1))+abs(pvl(2))  &
                     +abs(pvl(3))+abs(pvl(4)))/abs(volume)
                if (abs(vcheck1).lt.vcheck_tol) then
                   field1%val(inod)=0.0
                   do iloc=1,4
                      node=field2%mesh%ndglno((ielem-1)*4+iloc)
                      shapef=pvl(iloc)/volume
                      field1%val(inod)=field1%val(inod) + field2%val(node)*shapef
                   end do
                   vfield1%val(X_,inod)=0.0
                   vfield1%val(Y_,inod)=0.0
                   vfield1%val(Z_,inod)=0.0
                   do iloc=1,4
                      node=field2%mesh%ndglno((ielem-1)*4+iloc)
                      shapef=pvl(iloc)/volume
                      vfield1%val(X_,inod)= vfield1%val(X_,inod) + &
                           vfield2%val(X_,node)*shapef
                      vfield1%val(Y_,inod)= vfield1%val(Y_,inod) + &
                           vfield2%val(Y_,node)*shapef
                      vfield1%val(Z_,inod)= vfield1%val(Z_,inod) + &
                           vfield2%val(Z_,node)*shapef
                   end do
                   node_to_particle1(inod)=ipart
                   nbrute=nbrute+1
                   exit
                end if
             end do
          end if
       end if!(if inside the bin)
    end do

    ewrite(0,*) 'nodes mapped:...by nn search  :',nfastest,', and by brute force:',nbrute

    deallocate(eelist)
    deallocate(nelist)
    deallocate(nelistlgth)
    deallocate(bin)
    deallocate(ind)
    call cpu_time(sec2)
    ewrite(0,*) 'time spent on nn search: ',sec2-sec1
  end subroutine mesh2mesh_3d

  subroutine check_if_in_bin(inside,bin_x_max,bin_y_max,bin_z_max, &
       bin_x_min,bin_y_min,bin_z_min,                              &
       x,y,z)
    implicit none
    logical, intent(inout) :: inside
    real, intent(in) :: bin_x_max,bin_y_max,bin_z_max
    real, intent(in) :: bin_x_min,bin_y_min,bin_z_min
    real, intent(in) :: x,y,z

    if ((x<=bin_x_max.and.x>=bin_x_min).and.&
         (y<=bin_y_max.and.y>=bin_y_min).and.&
         (z<=bin_z_max.and.z>=bin_z_min)) then
       inside=.true.
    end if

  end subroutine check_if_in_bin

  subroutine get_bin_size(state,position_name1,position_name2,ntry, &
       number_of_bins_x,number_of_bins_y,number_of_bins_z,          &
       bin_x_min,bin_y_min,bin_z_min,                               &
       bin_x_max,bin_y_max,bin_z_max,                               &
       dx,dy,dz)

    implicit none
    type(state_type), dimension(:), intent(inout) :: state
    character(len=*) ,intent(in) :: position_name1,position_name2
    integer, intent(in)    :: ntry
    integer, intent(inout) :: number_of_bins_x,number_of_bins_y,number_of_bins_z
    real,    intent(inout) :: bin_x_min,bin_y_min,bin_z_min
    real,    intent(inout) :: bin_x_max,bin_y_max,bin_z_max
    real,    intent(inout) :: dx,dy,dz

    type (vector_field) :: position1,position2
    type (scalar_field) :: lumped_mass
    integer :: nonods1,nonods2
    real :: bin_x_min1,bin_y_min1,bin_z_min1
    real :: bin_x_max1,bin_y_max1,bin_z_max1
    real :: bin_x_min2,bin_y_min2,bin_z_min2
    real :: bin_x_max2,bin_y_max2,bin_z_max2
    real :: min_el_size,deltax,deltay,deltaz
    integer :: max_number_of_bins_x,max_number_of_bins_y,max_number_of_bins_z

    character(len=OPTION_PATH_LEN) :: tmpstring

    tmpstring="/imported_solids/solid_type/mapping_type"
    call get_option(trim(tmpstring)//"/max_number_of_bins_x",max_number_of_bins_x)
    call get_option(trim(tmpstring)//"/max_number_of_bins_y",max_number_of_bins_y)
    call get_option(trim(tmpstring)//"/max_number_of_bins_z",max_number_of_bins_z)

    nonods1=node_count(position1)
    nonods2=node_count(position2)
    !-----------------------------------------------------
    !extract coordinates corresponding to the mesh of field1
    position1=extract_vector_field(state(1),trim(position_name1))
    bin_x_max1=maxval(position1%val(X_,:),nonods1)
    bin_y_max1=maxval(position1%val(Y_,:),nonods1)
    bin_z_max1=maxval(position1%val(Z_,:),nonods1)
    bin_x_min1=minval(position1%val(X_,:),nonods1)
    bin_y_min1=minval(position1%val(Y_,:),nonods1)
    bin_z_min1=minval(position1%val(Z_,:),nonods1)

    position2=extract_vector_field(state(1),trim(position_name2))
    bin_x_max2=maxval(position2%val(X_,:),nonods2)
    bin_y_max2=maxval(position2%val(Y_,:),nonods2)
    bin_z_max2=maxval(position2%val(Z_,:),nonods2)
    bin_x_min2=minval(position2%val(X_,:),nonods2)
    bin_y_min2=minval(position2%val(Y_,:),nonods2)
    bin_z_min2=minval(position2%val(Z_,:),nonods2)

    !first get the maximum size of the interpolating box (bin size)
    bin_x_max=min(bin_x_max2,bin_x_max1)
    bin_y_max=min(bin_y_max2,bin_y_max1)
    bin_z_max=min(bin_z_max2,bin_z_max1)

    bin_x_min=max(bin_x_min2,bin_x_min1)
    bin_y_min=max(bin_y_min2,bin_y_min1)
    bin_z_min=max(bin_z_min2,bin_z_min1)

    deltax = abs(bin_x_max-bin_x_min)*0.05
    deltay = abs(bin_y_max-bin_y_min)*0.05
    deltaz = abs(bin_z_max-bin_z_min)*0.05

    bin_x_max = bin_x_max+deltax
    bin_y_max = bin_y_max+deltay
    bin_z_max = bin_z_max+deltaz

    bin_x_min = bin_x_min-deltax
    bin_y_min = bin_y_min-deltay
    bin_z_min = bin_z_min-deltaz

    !next, establish the number of bins by

    call allocate(lumped_mass, position2%mesh, "Lumped mass")
    !generate lumped mass matrix
    call compute_lumped_mass(position2, lumped_mass)
    min_el_size=minval(lumped_mass%val**(1./3.),nonods2)
    call deallocate(lumped_mass)

    number_of_bins_x=max(min(int(abs(bin_x_max-bin_x_min)/(3*min_el_size)),max_number_of_bins_x),1)
    number_of_bins_y=max(min(int(abs(bin_y_max-bin_y_min)/(3*min_el_size)),max_number_of_bins_y),1)
    number_of_bins_z=max(min(int(abs(bin_z_max-bin_z_min)/(3*min_el_size)),max_number_of_bins_z),1)

    dx=abs(bin_x_max-bin_x_min)/number_of_bins_x
    dy=abs(bin_y_max-bin_y_min)/number_of_bins_y
    dz=abs(bin_z_max-bin_z_min)/number_of_bins_z

    ewrite(0,*) 'number of bins: ',number_of_bins_x,number_of_bins_y,number_of_bins_z
    ewrite(0,*) 'bin_x_max1,bin_x_max2',bin_x_max1,bin_x_max2
    ewrite(0,*) 'bin_y_max1,bin_y_max2',bin_y_max1,bin_y_max2
    ewrite(0,*) 'bin_z_max1,bin_z_max2',bin_z_max1,bin_z_max2
    ewrite(0,*) 'bin_x_min1,bin_x_min2',bin_x_min1,bin_x_min2
    ewrite(0,*) 'bin_y_min1,bin_y_min2',bin_y_min1,bin_y_min2
    ewrite(0,*) 'bin_z_min1,bin_z_min2',bin_z_min1,bin_z_min2


    ewrite(0,*) 'dx,dy,dz : ',dx,dy,dz
    ewrite(0,*) 'binx : ',abs(bin_x_max-bin_x_min)
    ewrite(0,*) 'biny : ',abs(bin_y_max-bin_y_min)
    ewrite(0,*) 'binz : ',abs(bin_z_max-bin_z_min)
    ewrite(0,*) 'min_el_size: ',min_el_size

    !    DO A CHECK ON ALL NODES WHICH ARE LEFT IN AND OUT OF THE BIN BOUNDING BOX.
    !       MAKE SURE EVERYTHING IS RIGHT THERE. MAKE COUNTERS AND PRINT AMOUNT OF NODES>
  end subroutine get_bin_size

  subroutine construct_bin(                                &
       state,                                              &
       number_of_bins_x,number_of_bins_y,number_of_bins_z, &
       bin_x_max,bin_y_max,bin_z_max,                      &
       bin_x_min,bin_y_min,bin_z_min,                      &
       position_name,dx,dy,dz)

    implicit none
    type(state_type), dimension(:), intent(inout) :: state
    integer, intent(in)    :: number_of_bins_x,number_of_bins_y,number_of_bins_z
    real, intent(in)    :: bin_x_max,bin_y_max,bin_z_max
    real, intent(in)    :: bin_x_min,bin_y_min,bin_z_min
    character(len=OPTION_PATH_LEN), intent(in) :: position_name
    real,    intent(in) :: dx,dy,dz

    type (vector_field) :: position
    integer nloc,ielem,i,j,k,ii,iii,node(4)
    real avg(3)

    position=extract_vector_field(state(1),trim(position_name))
    nloc=4
    allocate(bin(number_of_bins_x,number_of_bins_y,number_of_bins_z))
    allocate(ind(27,3))
    call create_ind(ind)
    bin=-1
    !create bin: a loop over the external mesh's elements calculates each element's
    !volumetric center (xavg,yavg,zavg) and then assings it the closest integer valued
    !coordinates, which correspond to the bin.
    !note: try a nodal search instead of element search to find the elements
    !use the nelist!!!!

    do ielem=1,element_count(position)
       do iii=1,nloc
          node(iii)=position%mesh%ndglno((ielem-1)*nloc + iii)
       end do
       do ii=1,3
          avg(ii)=0
          do iii=1,nloc
             avg(ii)=avg(ii) + position%val(ii,node(iii))
          end do
          avg(ii)=avg(ii)/nloc
       end do

       i=int(abs((avg(1)-bin_x_min)/dx))+1 !normally this would have a +1
       j=int(abs((avg(2)-bin_y_min)/dy))+1 !but since there is a -1 layer surrounding the bin
       k=int(abs((avg(3)-bin_z_min)/dz))+1 !it needs to be skipped

       if(bin(i,j,k)==-1) then
          bin(i,j,k)=ielem
       end if
    end do
  end subroutine construct_bin

  subroutine calculate_particle_drag(state,ipart)
    implicit none
    type(state_type), dimension(:), intent(inout) :: state
    integer, intent(in)    :: ipart

    !local vars
    integer inod,counter
    real rx,ry,rz,forcexs,forceys,forcezs
    type(scalar_field) lumped_mass,visualizesolid
    type(vector_field) svelocity,absorption,velocity_source,velocity
    type(vector_field) , pointer :: positions
    integer myrank

    myrank=GetRank()

    absorption=extract_vector_field(state(1),"VelocityAbsorption")
    velocity_source=extract_vector_field(state(1),"VelocitySource")
    velocity=extract_vector_field(state(1),"Velocity")
    svelocity=extract_vector_field(state(1),"SolidVelocity")
    positions=>extract_vector_field(state(1),"Coordinate")

    call allocate(lumped_mass, positions%mesh, "Lumped mass")
    !generate lumped mass matrix
    call compute_lumped_mass(positions, lumped_mass)

    forcex(ipart)=0.0
    forcey(ipart)=0.0
    forcez(ipart)=0.0
    torquex(ipart)=0.0
    torquey(ipart)=0.0
    torquez(ipart)=0.0
    counter=0

    visualizesolid=extract_scalar_field(state(1),"VisualizeSolid")

    do inod=1,node_count(positions)
       if (node_to_particle(inod)==ipart) then

          forcexs = absorption%val(X_,inod)* &
               (velocity%val(X_,inod)-svelocity%val(X_,inod))*lumped_mass%val(inod)
          forceys = absorption%val(Y_,inod)* &
               (velocity%val(Y_,inod)-svelocity%val(Y_,inod))*lumped_mass%val(inod)
          forcezs = absorption%val(Z_,inod)* &
               (velocity%val(Z_,inod)-svelocity%val(Z_,inod))*lumped_mass%val(inod)

          forcex(ipart) = forcex(ipart) + forcexs
          forcey(ipart) = forcey(ipart) + forceys
          forcez(ipart) = forcez(ipart) + forcezs

          rx = positions%val(X_,inod)-centerpx(ipart)
          ry = positions%val(Y_,inod)-centerpy(ipart)
          rz = positions%val(Z_,inod)-centerpz(ipart)

          torquex(ipart) = torquex(ipart) - rz*forceys + ry*forcezs
          torquey(ipart) = torquey(ipart) - rx*forcezs + rz*forcexs
          torquez(ipart) = torquez(ipart) - ry*forcexs + rx*forceys

          counter=counter+1
       endif
    end do
    ewrite(0,*) 'Amount of nodes inside particle: ',ipart,' is: ',counter
    call deallocate(lumped_mass)
    call allsum(forcex(ipart))
    call allsum(forcey(ipart))
    call allsum(forcez(ipart))
    call allsum(torquex(ipart))
    call allsum(torquey(ipart))
    call allsum(torquez(ipart))

  end subroutine calculate_particle_drag
  subroutine calculate_solid_force_vector(state)

    implicit none
    type(state_type), dimension(:), intent(inout) :: state

    !local vars
    type(scalar_field) lumped_mass,visualizesolid
    type(vector_field) svelocity,absorption,velocity,solid_force
    type(vector_field) , pointer :: positions

    absorption=extract_vector_field(state(1),"VelocityAbsorption")
    velocity=extract_vector_field(state(1),"Velocity")
    svelocity=extract_vector_field(state(1),"SolidVelocity")
    visualizesolid=extract_scalar_field(state(1),"VisualizeSolid")
    positions=>extract_vector_field(state(1),"Coordinate")
    solid_force=extract_vector_field(state(1),"SolidForce")
    call allocate(lumped_mass, positions%mesh, "Lumped mass")

    !generate lumped mass matrix
    call compute_lumped_mass(positions, lumped_mass)

    solid_force%val(X_,:) = absorption%val(X_,:)*visualizesolid%val*(velocity%val(X_,:) - &
         svelocity%val(X_,:))*lumped_mass%val
    solid_force%val(Y_,:) = absorption%val(Y_,:)*visualizesolid%val*(velocity%val(Y_,:) - &
         svelocity%val(Y_,:))*lumped_mass%val
    solid_force%val(Z_,:) = absorption%val(Z_,:)*visualizesolid%val*(velocity%val(Z_,:) - &
         svelocity%val(Z_,:))*lumped_mass%val

    call deallocate(lumped_mass)

  end subroutine calculate_solid_force_vector

  subroutine calculate_particle_mass(state)
    implicit none
    type(state_type), dimension(:), intent(inout) :: state
    !local vars
    integer inod,myrank,nprivatenodes,ipart
    type(vector_field) positions
    type(scalar_field) lumped_mass,solid_concentration
    real :: total_volume
    myrank=GetRank()

    solid_concentration=extract_scalar_field(state(1),"SolidConcentration")
    positions=extract_vector_field(state(1),"Coordinate")
    call allocate(lumped_mass, positions%mesh, "Lumped mass")
    !generate lumped mass matrix
    call compute_lumped_mass(positions, lumped_mass)
    particle_mass=0.0
    total_volume=0.0

    nprivatenodes = nowned_nodes(positions)

    ewrite(0,*)'Private Nodes: ',nprivatenodes
    do ipart=1,nparticles
       do inod=1, nprivatenodes
          if(node_to_particle(inod)>0) then
             particle_mass(node_to_particle(inod))=particle_mass(node_to_particle(inod)) + &
                  solid_concentration%val(inod)*lumped_mass%val(inod)*solid_density
          end if
          total_volume=total_volume + solid_concentration%val(inod)*lumped_mass%val(inod)
       end do
    end do
    call allsum(total_volume)
    call allsumv(particle_mass)
    call deallocate(lumped_mass)
    if (myrank==0) then
       ewrite(0,*) 'Total Solid Volume: ',total_volume
       ewrite(0,*) 'Total Solid Mass  : ',total_volume*solid_density
    end if
    do ipart=1,nparticles
       ewrite(0,*) 'Particle : ',ipart,' mass = ',particle_mass(ipart)
    end do

  end subroutine calculate_particle_mass

  subroutine output_particle_info(nparticles,    &
       centerpx,centerpy,centerpz,               &
       pvelx,pvely,pvelz,                        &
       pavelx,pavely,pavelz,                     &
       forcex,forcey,forcez,                     &
       torquex,torquey,torquez)

    implicit none
    integer, intent(in) :: nparticles
    real,    intent(in), dimension(nparticles) :: centerpx,centerpy,centerpz
    real,    intent(in), dimension(nparticles) :: pvelx,pvely,pvelz
    real,    intent(in), dimension(nparticles) :: pavelx,pavely,pavelz
    real,    intent(in), dimension(nparticles) :: forcex,forcey,forcez
    real,    intent(in), dimension(nparticles) :: torquex,torquey,torquez
    character*5 efilnm
    !local variables
    integer ipart,myrank
    myrank=GetRank()
    if(myrank==0) then
       open(1266,file='particle_position',position='append')
       open(1267,file='particle_velocity',position='append')
       open(1268,file='particle_angular_velocity',position='append')
       open(1269,file='particle_force',position='append')
       open(1270,file='particle_torque',position='append')

       call i4_to_s_left (nparticles*3, efilnm )
       if (nparticles.le.5) then
          if(nparticles*3.le.9) then
             write(1266,'(i0,'//efilnm(1:1)//'e20.10)') iteration,           &
                  (centerpx(ipart),ipart=1,nparticles),                      &
                  (centerpy(ipart),ipart=1,nparticles),                      &
                  (centerpz(ipart),ipart=1,nparticles)
             write(1267,'(i0,'//efilnm(1:1)//'e20.10)') iteration,           &
                  (pvelx(ipart),ipart=1,nparticles),                         &
                  (pvely(ipart),ipart=1,nparticles),                         &
                  (pvelz(ipart),ipart=1,nparticles)

             write(1268,'(i0,'//efilnm(1:1)//'e20.10)') iteration,           &
                  (pavelx(ipart),ipart=1,nparticles),                        &
                  (pavely(ipart),ipart=1,nparticles),                        &
                  (pavelz(ipart),ipart=1,nparticles)

             write(1269,'(i0,'//efilnm(1:1)//'e20.10)') iteration,           &
                  (forcex(ipart),ipart=1,nparticles),                        &
                  (forcey(ipart),ipart=1,nparticles),                        &
                  (forcez(ipart),ipart=1,nparticles)
             write(1270,'(i0,'//efilnm(1:1)//'e20.10)') iteration,           &
                  (torquex(ipart),ipart=1,nparticles),                       &
                  (torquey(ipart),ipart=1,nparticles),                       &
                  (torquez(ipart),ipart=1,nparticles)
          else
             write(1266,'(i0,'//efilnm(1:2)//'e20.10)') iteration,           &
                  (centerpx(ipart),ipart=1,nparticles),                      &
                  (centerpy(ipart),ipart=1,nparticles),                      &
                  (centerpz(ipart),ipart=1,nparticles)
             write(1267,'(i0,'//efilnm(1:2)//'e20.10)') iteration,           &
                  (pvelx(ipart),ipart=1,nparticles),                         &
                  (pvely(ipart),ipart=1,nparticles),                         &
                  (pvelz(ipart),ipart=1,nparticles)

             write(1268,'(i0,'//efilnm(1:2)//'e20.10)') iteration,           &
                  (pavelx(ipart),ipart=1,nparticles),                        &
                  (pavely(ipart),ipart=1,nparticles),                        &
                  (pavelz(ipart),ipart=1,nparticles)
             write(1269,'(i0,'//efilnm(1:2)//'e20.10)') iteration,           &
                  (forcex(ipart),ipart=1,nparticles),                        &
                  (forcey(ipart),ipart=1,nparticles),                        &
                  (forcez(ipart),ipart=1,nparticles)
             write(1270,'(i0,'//efilnm(1:2)//'e20.10)') iteration,           &
                  (torquex(ipart),ipart=1,nparticles),                       &
                  (torquey(ipart),ipart=1,nparticles),                       &
                  (torquez(ipart),ipart=1,nparticles)
          end if
       else
          ewrite(0,*)'particle information output is only available for up to 5 particles'
          FLExit('dying... turn off output drag, or choose different number of particles')
       end if
       close(1266)
       close(1267)
       close(1268)
       close(1269)
       close(1270)
    end if
  end subroutine output_particle_info

  subroutine create_ind(ind)
    implicit none
    integer, intent(inout) :: ind(27,3)

    ind(1,1)=0
    ind(2,1)=0
    ind(3,1)=1
    ind(4,1)=1
    ind(5,1)=-1
    ind(6,1)=-1
    ind(7,1)=-1
    ind(8,1)=0
    ind(9,1)=1
    ind(10,1)=0
    ind(11,1)=0
    ind(12,1)=1
    ind(13,1)=1
    ind(14,1)=-1
    ind(15,1)=-1
    ind(16,1)=-1
    ind(17,1)=0
    ind(18,1)=1
    ind(19,1)=0
    ind(20,1)=0
    ind(21,1)=1
    ind(22,1)=1
    ind(23,1)=-1
    ind(24,1)=-1
    ind(25,1)=-1
    ind(26,1)=0
    ind(27,1)=1

    ind(1,2)=0
    ind(2,2)=1
    ind(3,2)=0
    ind(4,2)=1
    ind(5,2)=1
    ind(6,2)=0
    ind(7,2)=-1
    ind(8,2)=-1
    ind(9,2)=-1
    ind(10,2)=0
    ind(11,2)=1
    ind(12,2)=0
    ind(13,2)=1
    ind(14,2)=1
    ind(15,2)=0
    ind(16,2)=-1
    ind(17,2)=-1
    ind(18,2)=-1
    ind(19,2)=0
    ind(20,2)=1
    ind(21,2)=0
    ind(22,2)=1
    ind(23,2)=1
    ind(24,2)=0
    ind(25,2)=-1
    ind(26,2)=-1
    ind(27,2)=-1

    ind(1,3)=-1
    ind(2,3)=-1
    ind(3,3)=-1
    ind(4,3)=-1
    ind(5,3)=-1
    ind(6,3)=-1
    ind(7,3)=-1
    ind(8,3)=-1
    ind(9,3)=-1
    ind(10,3)=0
    ind(11,3)=0
    ind(12,3)=0
    ind(13,3)=0
    ind(14,3)=0
    ind(15,3)=0
    ind(16,3)=0
    ind(17,3)=0
    ind(18,3)=0
    ind(19,3)=1
    ind(20,3)=1
    ind(21,3)=1
    ind(22,3)=1
    ind(23,3)=1
    ind(24,3)=1
    ind(25,3)=1
    ind(26,3)=1
    ind(27,3)=1
  end subroutine create_ind

  subroutine i4_to_s_left ( intval, s )

    !*****************************************************************************80
    !
    !! i4_to_s_left converts an i4 to a left-justified string.
    !
    !  examples:
    !
    !    assume that s is 6 characters long:
    !
    !    intval  s
    !
    !         1  1
    !        -1  -1
    !         0  0
    !      1952  1952
    !    123456  123456
    !   1234567  ******  <-- not enough room!
    !
    !  modified:
    !
    !    28 july 2000
    !
    !  author:
    !
    !    john burkardt
    !
    !  parameters:
    !
    !    input, integer intval, an integer to be converted.
    !
    !    output, character ( len = * ) s, the representation of the integer.
    !    the integer will be left-justified.  if there is not enough space,
    !    the string will be filled with stars.
    !
    implicit none

    character c
    integer i
    integer idig
    integer ihi
    integer ilo
    integer intval
    integer ipos
    integer ival
    character*5 s

    s = ' '

    ilo = 1
    ihi = len ( s )

    if ( ihi <= 0 ) then
       return
    end if
    !
    !  make a copy of the integer.
    !
    ival = intval
    !
    !  handle the negative sign.
    !
    if ( ival < 0 ) then

       if ( ihi <= 1 ) then
          s(1:1) = '*'
          return
       end if

       ival = -ival
       s(1:1) = '-'
       ilo = 2

    end if
    !
    !  the absolute value of the integer goes into s(ilo:ihi).
    !
    ipos = ihi
    !
    !  find the last digit of ival, strip it off, and stick it into the string.
    !
    do

       idig = mod ( ival, 10 )
       ival = ival / 10

       if ( ipos < ilo ) then
          do i = 1, ihi
             s(i:i) = '*'
          end do
          return
       end if

       call digit_to_ch ( idig, c )

       s(ipos:ipos) = c
       ipos = ipos - 1

       if ( ival == 0 ) then
          exit
       end if

    end do
    !
    !  shift the string to the left.
    !
    s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
    s(ilo+ihi-ipos:ihi) = ' '

    return
  end subroutine i4_to_s_left

  subroutine digit_to_ch ( digit, ch )

    !*****************************************************************************80
    !
    !! digit_to_ch returns the character representation of a decimal digit.
    !
    !  discussion:
    !
    !    instead of char, we now use the achar function, which
    !    guarantees the ascii collating sequence.
    !
    !  example:
    !
    !    digit   ch
    !    -----  ---
    !      0    '0'
    !      1    '1'
    !    ...    ...
    !      9    '9'
    !     17    '*'
    !
    !  modified:
    !
    !    04 august 1999
    !
    !  author:
    !
    !    john burkardt
    !
    !  parameters:
    !
    !    input, integer digit, the digit value between 0 and 9.
    !
    !    output, character ch, the corresponding character.
    !
    implicit none

    character ch
    integer digit

    if ( 0 <= digit .and. digit <= 9 ) then

       ch = achar ( digit + 48 )

    else

       ch = '*'

    end if

    return
  end subroutine digit_to_ch

  subroutine generate_nelist(mesh)
    implicit none
    integer ielem,iloc,node
    integer maxconec
    type (mesh_type),intent(in) :: mesh
    parameter(maxconec=400)
    !    ewrite(0,*) 'nloc :',nloc,node_count(mesh),node_count(mesh)
    allocate(nelist(node_count(mesh),maxconec))
    allocate(nelistlgth(node_count(mesh)))
    nelistlgth=0
    nelist=0
    do ielem=1,element_count(mesh)
       do iloc=1,nloc
          node=mesh%ndglno(nloc*(ielem-1)+iloc)
          nelistlgth(node)=nelistlgth(node)+1
          if(nelistlgth(node).gt.maxconec) then
             ewrite(0,*)'too many connectivities'
          end if
          nelist(node,nelistlgth(node))=ielem
       end do
    end do
  end subroutine generate_nelist

  subroutine generate_eelist(mesh)
    implicit none
    integer ielem,iele,iloc,jloc,icon
    integer node,node1,node2,node3,node4,nface
    integer, allocatable :: index1(:,:)
    type (mesh_type),intent(in) :: mesh

    allocate(eelist(element_count(mesh),nloc))
    eelist=-1
    !-----------------------------------------------------
    !generate element to element list for mesh
    allocate(index1(4,3))
    index1(1,1)=2
    index1(1,2)=3
    index1(1,3)=4
    index1(2,1)=3
    index1(2,2)=4
    index1(2,3)=1
    index1(3,1)=4
    index1(3,2)=1
    index1(3,3)=2
    index1(4,1)=1
    index1(4,2)=2
    index1(4,3)=3

    do ielem=1,element_count(mesh)
       do iloc=1,nloc
          !Node1 is the node opposing the face being checked.
          !this face is consisting of node2,node3, and node4
          !elements surrounding one of these three nodes are checked.(node2)
          !when an element matches node2,node3,and node4, then there is a face match
          node1 = mesh%ndglno((ielem-1)*nloc + iloc)
          node2 = mesh%ndglno((ielem-1)*nloc + index1(iloc,1))
          node3 = mesh%ndglno((ielem-1)*nloc + index1(iloc,2))
          node4 = mesh%ndglno((ielem-1)*nloc + index1(iloc,3))

          do icon=1,nelistlgth(node2)
             iele=nelist(node2,icon)
             nface=0
             if(iele/=ielem) then
                do jloc=1,nloc
                   node=mesh%ndglno((iele-1)*nloc + jloc)
                   if(node==node2.or.node==node3.or.node==node4) then
                      nface=nface+1
                   end if
                end do
                if(nface==3) then
                   eelist(ielem,iloc)=iele
                   exit
                end if
             end if
          end do
       end do
    end do
    deallocate(index1)

  end subroutine generate_eelist

  subroutine check_volumes( &
       state,         &
       nonods,totele, &
       ndglno,nloc,   &
       x,y,z)
    implicit none
    type(state_type),dimension(:),intent(in) :: state
    integer, intent(in)    :: nonods
    integer, intent(in)    :: totele
    integer, intent(in)    :: ndglno(4*totele),nloc
    real,    intent(in)    :: x(nonods),y(nonods),z(nonods)

    integer iloc,ielem,node
    real xx(nloc),yy(nloc),zz(nloc),vl
    type(scalar_field) :: density,matvolfrac
    type(vector_field) :: position,velocity
    type(mesh_type)    :: mesh
    integer, allocatable , dimension(:) :: snodes
    integer sele,stat,sngi
    real, allocatable, dimension(:) :: face_detwei
    character(len=OPTION_PATH_LEN) :: tmpstring
    integer nsbc,nvbc,ph
    integer state_to_use

    phaseloop: do ph = 1, size(state)
       write(tmpstring, '(a,i0,a)') "/material_phase[",ph-1,"]"
       if(have_option(trim(tmpstring)//"/scalar_field::MaterialVolumeFraction/prognostic")) then
          state_to_use=ph
       end if
    end do phaseloop
    ewrite(2,*) 'element volumes--- jem'
    do ielem=1,totele
       do iloc=1,nloc
          node = ndglno((ielem-1)*nloc+iloc)
          xx(iloc) = x(node)
          yy(iloc) = y(node)
          zz(iloc) = z(node)
       end do
       vl=element_volume_s(nloc,xx,yy,zz)
       ewrite(2,*) 'volume: ',vl
    end do
    ewrite(2,*) 'end element volumes--- jem'
    density=extract_scalar_field(state(state_to_use),"Density",stat)
    position = extract_vector_field(state(state_to_use), "Coordinate")
    matvolfrac=extract_scalar_field(state(state_to_use),"MaterialVolumeFraction",stat)
    if(stat/=0) then
       ewrite(-1,*) 'Could not extract MaterialVolumeFraction'
       FLExit("You need to add the MaterialVolumeFraction scalar field to your options")
    end if
    velocity=extract_vector_field(state(state_to_use),"Velocity",stat)
    if(stat/=0) then
       ewrite(-1,*)('could not extract Velocity')
       FLExit("You need to add the Velocity vector field to your options")
    end if

    mesh=extract_mesh(state(state_to_use),'CoordinateMesh')
    sngi=face_ngi(density, 1)
    allocate(face_detwei(1:sngi))
    allocate(snodes(3))
    ewrite(2,*) 'surface elements--- jem'
    ewrite(2,*)'number of face elements :',surface_element_count(mesh)
    do sele=1, surface_element_count(mesh)
       snodes=face_global_nodes(mesh, sele)
       call transform_facet_to_physical(position, sele, face_detwei)
       ewrite(2,*) snodes(1),snodes(2),snodes(3),sum(face_detwei),surface_element_id(density,sele)
    end do
    ewrite(2,*) 'surface elements--- jem'
    deallocate(snodes)
    deallocate(face_detwei)

    nsbc=get_boundary_condition_count(matvolfrac)
    nvbc=get_boundary_condition_count(velocity)
    ewrite(2,*) 'boundary conditions matvolfrac', get_boundary_condition_count(matvolfrac),&
         get_boundary_condition_count(matvolfrac)
    ewrite(2,*) 'boundary conditions velocity', get_boundary_condition_count(velocity),&
         get_boundary_condition_count(velocity)
    if(nsbc.ne.0) then
       ewrite(2,*) 'matvolfrac boundary condition count is :',nsbc,'  --jem'
    else
       ewrite(2,*) 'matvolfrac boundary condition count is zero  --jem'
    end if


    if(nvbc.ne.0) then
       ewrite(2,*) 'velocity boundary condition count is :',nvbc,'  --jem'
    else
       ewrite(2,*) 'velocity boundary condition count is :',nvbc,'  --jem'
    end if
    call write_triangle_files("from_gid",state(state_to_use),position%mesh)


  end subroutine check_volumes

  !  subroutine modify_itinoi(itinoi)
  !    implicit none
  !    integer, intent(inout) :: itinoi
  !    integer             :: period
  !    integer, save       :: iteration
  !    data iteration /0/
  !
  !    call get_option("/mesh_adaptivity/hr_adaptivity/period_in_timesteps",period)
  !
  !    if (iteration==period) then
  !      itinoi=20
  !    elseif (iteration==period+1) then
  !      itinoi=15
  !    elseif (iteration==period+2) then
  !      itinoi=10
  !    elseif (iteration==period+3) then
  !      itinoi=5
  !      iteration=0
  !    else
  !      itinoi=1
  !    end if
  !    ewrite(0,*) 'itinoi modified to: ',itinoi
  !    iteration=iteration+1
  !
  !  end subroutine modify_itinoi


  subroutine gid2triangle(state)
    implicit none
    type(state_type),dimension(:),intent(in) :: state
    type(vector_field) :: position
    character(len=OPTION_PATH_LEN) :: tmpstring,filename
    integer ph,state_to_use

    ewrite(0,*) '----------------------------------------------------------------'
    if (have_option("/geometry/gid2triangle")) then
       phaseloop: do ph = 1, size(state)
          write(tmpstring, '(a,i0,a)') "/material_phase[",ph-1,"]"
          if((have_option(trim(tmpstring)//"/scalar_field::MaterialVolumeFraction/prognostic")).or. &
               (have_option(trim(tmpstring)//"/scalar_field::Pressure/prognostic")).or. &
               (have_option(trim(tmpstring)//"/scalar_field::Velocity/prognostic"))) then
             state_to_use=ph
             exit phaseloop
          end if
       end do phaseloop
       call get_option("/geometry/gid2triangle/file_name",filename)
       position = extract_vector_field(state(state_to_use), "Coordinate")
       call write_triangle_files(trim(filename),state(state_to_use),position%mesh)
       ewrite(0,*) '-finished translating gid mesh into triangle format'
       ewrite(0,*) '----------------------------------------------------------------'
       ewrite(-1,*)' now that you have output the gid mesh in triangle format,'
       ewrite(-1,*)' simply reload your flml file and change the import format'
       ewrite(-1,*)' to "triangle". remember to change the file name as well.'
       ewrite(-1,*)' note: this works safely for single mesh problems only!!'
       FLExit('exiting. re-run problem.')
    else
       return
    end if

  end subroutine gid2triangle


  subroutine drag_on_surface(state)
    implicit none
    type(state_type),dimension(:),intent(inout) :: state

    !Local Memory
    integer i,count,node,iloc,ielem,nod1,nod2,nod3,nloc
    type(scalar_field) pressure
    type(mesh_type) mesh
    type(vector_field) positions,flow_velocity
    real, allocatable, dimension(:)::snodes
    integer surface_number
    real :: fx,fy,fz,total_area,dn(3),normalx,normaly,normalz,dnormal
    real :: edge1(3),edge2(3),tauxx,tauxy,tauxz,tauyy,tauyz,tauzz,area
    real :: xx(3),yy(3),zz(3),pavg,visc
    mesh=extract_mesh(state(1),"CoordinateMesh")
    call get_option("/imported_solids/calculate_drag_on_surface/surface_id",surface_number)
    flow_velocity=extract_vector_field(state(1),"Velocity")
    positions=extract_vector_field(state(1),"Coordinate")
    pressure=extract_scalar_field(state(1),"Pressure")
    nloc=4
    fx=0.0
    fy=0.0
    fz=0.0
    total_area=0.0
    count=0
    visc=0.005
    allocate(snodes(3))
    do i=1,surface_element_count(mesh)
       if (mesh%faces%boundary_ids(i)==surface_number) then
          ewrite(0,*) 'face ',i, 'boundary_id ',mesh%faces%boundary_ids(i)
          snodes=face_global_nodes(mesh, i)
          ewrite(0,*) 'surface node list' ,snodes(1),snodes(2),snodes(3)
          count=count+1

          ielem=positions%mesh%faces%face_element_list(i)
          ewrite(0,*) 'surface element: ielem: ',ielem
          do iloc=1,nloc
             node=mesh%ndglno(nloc*(ielem-1) + iloc)
             xx(iloc)=positions%val(X_,node)
             yy(iloc)=positions%val(Y_,node)
             zz(iloc)=positions%val(Z_,node)
          end do

          dn=element_derivatives(nloc,xx,yy,zz)
          !Now calculating the components of the stress tensor
          !(Piecewise constant inside the element)

          TAUXX = 0.0
          TAUXY = 0.0
          TAUXZ = 0.0
          TAUYY = 0.0
          TAUYZ = 0.0
          TAUZZ = 0.0

          !      EWRITE(0,*) '--------------------------------------------------'
          DO ILOC = 1,NLOC
             Node = positions%mesh%ndglno((ielem-1)*nloc + iloc )
             TAUXX = TAUXX + Flow_velocity%val(X_,Node)*dn(1)   ! Assuming linear tet basis funcions and hence constant derivs
             TAUXY = TAUXY + Flow_velocity%val(X_,Node)*dn(2) + Flow_velocity%val(Y_,Node)*dn(1)
             TAUXZ = TAUXZ + Flow_velocity%val(X_,Node)*dn(3) + Flow_velocity%val(Z_,Node)*dn(1)
             TAUYY = TAUYY + Flow_velocity%val(Y_,Node)*dn(2)
             TAUYZ = TAUYZ + Flow_velocity%val(Y_,Node)*dn(3) + Flow_velocity%val(Z_,Node)*dn(2)
             TAUZZ = TAUZZ + Flow_velocity%val(Z_,Node)*dn(3)
             !          EWRITE(0,'(A,7E12.5)') 'VELOCITIES:',U(Node),V(Node),W(Node),VISC,dn(1),dn(2),dn(3)
          END DO

          TAUXX = 2.0*VISC*TAUXX
          TAUYY = 2.0*VISC*TAUYY
          TAUZZ = 2.0*VISC*TAUZZ
          TAUXY =     VISC*TAUXY
          TAUYZ =     VISC*TAUYZ
          TAUZZ =     VISC*TAUZZ
          !Finished calculating stress tensor

          !Now going to calculate forces
          NOD1=snodes(1)
          NOD2=snodes(2)
          NOD3=snodes(3)
          !Calculate area of surface element.

          !Now going to formulate the two vectors that define the plane that contains the surface element
          !Vector 1 goes from node 1 to node 2
          EDGE1(1)=positions%val(X_,NOD2)-positions%val(X_,NOD1)
          EDGE1(2)=positions%val(Y_,NOD2)-positions%val(Y_,NOD1)
          EDGE1(3)=positions%val(Z_,NOD2)-positions%val(Z_,NOD1)

          !Vector 2 goes from node 1 to node 3
          EDGE1(1)=positions%val(X_,NOD3)-positions%val(X_,NOD1)
          EDGE1(2)=positions%val(Y_,NOD3)-positions%val(Y_,NOD1)
          EDGE1(3)=positions%val(Z_,NOD3)-positions%val(Z_,NOD1)

          !Now calculating the cross product of the two EDGE vectors.
          NORMALX=EDGE1(2)*EDGE2(3)-EDGE1(3)*EDGE2(2)
          NORMALY=-(EDGE1(1)*EDGE2(3)-EDGE1(3)*EDGE2(1))
          NORMALZ=EDGE1(1)*EDGE2(2)-EDGE1(2)*EDGE2(1)

          !Suing the cross product to calculate the area (0.5*(EDGE1 X EDGE2))
          AREA= 0.5*(NORMALX**2.0+NORMALY**2.0+NORMALZ**2.0)**0.5
          TOTAL_AREA = TOTAL_AREA + AREA
          !       EWRITE(0,*) 'AREA: ',AREA

          !Calculating the norm of the results to normalize
          DNORMAL=(NORMALX**2.0 + NORMALY**2.0 + NORMALZ**2.0)**0.5

          !Creating normalized vector
          NORMALX=NORMALX/DNORMAL
          NORMALY=NORMALY/DNORMAL
          NORMALZ=NORMALZ/DNORMAL

          !      IF(FIRST) THEN
          !         EWRITE(0,*) NOD1,NOD2,NOD3
          !         EWRITE(0,'(3E12.5)') NORMALX,NORMALY,NORMALZ
          !         EWRITE(0,'(3E12.5)') TAUXX,TAUXY,TAUXZ
          !         EWRITE(0,'(2E12.5)') TAUYY,TAUYZ
          !         EWRITE(0,'(1E12.5)') TAUZZ
          !         DO GI=1,NGI
          !            EWRITE(0,'(4E12.5)') NX(1,GI),NX(2,GI),NX(3,GI),NX(4,GI)
          !         END DO
          !       END IF
          !      EWRITE(0,*) '=================================================='
          IF(NORMALX**2.0+NORMALY**2.0+NORMALZ**2.0.GT.1.00001) THEN
             EWRITE(0,*) 'NORMAL VECTOR is not normalized for element',ielem
             STOP
          END IF

          !Now integrating the force in each direction acting on the surface in question.
          PAVG =(pressure%val(NOD1)+pressure%val(NOD2)+pressure%val(NOD3))/3.0
          FX = FX + (NORMALX*TAUXX + NORMALY*TAUXY + NORMALZ*TAUXZ - NORMALX*PAVG)*AREA
          FY = FY + (NORMALX*TAUXY + NORMALY*TAUYY + NORMALZ*TAUYZ - NORMALY*PAVG)*AREA
          FZ = FZ + (NORMALX*TAUXZ + NORMALY*TAUYZ + NORMALZ*TAUZZ - NORMALZ*PAVG)*AREA
       end if
    end do
    deallocate(snodes)

  end subroutine drag_on_surface

  subroutine hookforce(state)
    implicit none
    type(state_type), dimension(:), intent(inout) :: state
    type(scalar_field) density
    real fluid_density,dxx,dyy,dzz,dx2,dy2,dz2,minmass,minmass1
    real tcol,v0,delta,kn,timestep,vol,rsum,mdr,daux,deltan,mass
    real deltavx,deltavy,deltavz,oo
    real, allocatable,dimension(:) :: fx,fy,fz

    integer ipart,jpart,ipartminmass,ntimesteps,i
    !This routine was originally written in C by Xavier Garcia, and has been
    !translated here to Fortran
    allocate(fx(nparticles),fy(nparticles),fz(nparticles))

    density=extract_scalar_field(state(1),"Density")
    fluid_density=0.0
    oo=1.0
    if (buoyancy) fluid_density=maxval(density%val)
    if (oneway) oo=0.0


    !//we need the minimum mass to calculate the collision time tcol=sqrt(mass/stiffness)
    minmass=huge(0.0)
    ipartminmass=1
    do ipart=1,nparticles
       minmass1=solid_density*4*pie*(radii(ipart)**3.0)/3
       if(minmass1<minmass) then
          minmass=minmass1
          ipartminmass=ipart
       end if
    end do

    !calculate Kn
    delta=0.001*radii(ipartminmass)
    v0=1.0  !max expected velocity during simulation.
    kn=minmass*(v0**2.0)/(delta**2.0)
    tcol=delta/v0

    timestep=tcol/250.0

    if(dt<timestep) then
       timestep=dt/10.0
       ntimesteps=10
    else
       ntimesteps=1+int(dt/timestep)
    end if
    !    ewrite(0,*) 'time : ',timestep,dt,ntimesteps
    !    ewrite(0,*) 'kn : ',kn,delta,v0,tcol


    do i=1,ntimesteps
       !Contact checking


       do ipart=1,nparticles
          vol=4*pie*(radii(ipart)**3.0)/3.0
          if(i==1.and.ipart==1) then
             ewrite(0,*)'drag forcedata: ',vol*solid_density,solid_density,forcex(ipart),forcey(ipart),forcez(ipart)
             ewrite(0,*)'buoyancy data : ',(solid_density-fluid_density)*vol*gravity_x,(solid_density-fluid_density)*vol*gravity_y,&
                  (solid_density-fluid_density)*vol*gravity_z
          end if
          !calculate total forces acting on particle (constant for this timestep)
          fx(ipart)=oo*forcex(ipart) + (solid_density-fluid_density)*vol*gravity_x
          fy(ipart)=oo*forcey(ipart) + (solid_density-fluid_density)*vol*gravity_y
          fz(ipart)=oo*forcez(ipart) + (solid_density-fluid_density)*vol*gravity_z
          !          ewrite(0,*) 'vol :',vol,fx(ipart),fy(ipart),fz(ipart)
       end do
       do ipart=1,nparticles
          !Check for overlaps with walls
          !walls in the x direction
          if (abs(centerpx(ipart)-xmin)<abs(xmax-centerpx(ipart))) then
             dxx=centerpx(ipart)-xmin
          else
             dxx=xmax-centerpx(ipart)
          end if
          dx2=dxx**2.0
          !          ewrite(0,*) 'dx2 :',dx2,centerpx(ipart),centerpy(ipart),centerpz(ipart)
          !          ewrite(0,*) 'radii :',radii(ipart)
          if(dx2<radii(ipart)**2.0) then
             deltan=radii(ipart)-dx2**0.5
             daux=kn*deltan**1.5/(dx2**0.5)
             fx(ipart) = fx(ipart) + daux*dxx
             !             ewrite(0,*) 'HIT THE WALL IN X!!!'
          end if

          !walls in the y direction
          if (abs(centerpy(ipart)-ymin)<abs(ymax-centerpy(ipart))) then
             dyy=centerpy(ipart)-ymin
          else
             dyy=ymax-centerpy(ipart)
          end if
          dy2=dyy**2.0

          if(dy2<radii(ipart)**2.0) then
             deltan=radii(ipart)-dy2**0.5
             daux=kn*deltan**1.5/(dy2**0.5)
             fy(ipart) = fy(ipart) + daux*dyy
             !             ewrite(0,*) 'HIT THE WALL IN Y!!!'
          end if

          !walls in the z direction
          if (abs(centerpz(ipart)-zmin)<abs(zmax-centerpz(ipart))) then
             dzz=centerpz(ipart)-zmin
          else
             dzz=zmax-centerpz(ipart)
          end if
          dz2=dzz**2.0
          if(dz2<radii(ipart)**2.0) then
             deltan=radii(ipart)-dz2**0.5
             daux=kn*deltan**1.5/(dz2**0.5)
             fz(ipart) = fz(ipart) + daux*dzz
             !             ewrite(0,*) 'HIT THE WALL IN Z!!!'
          end if
          !          ewrite(0,*) 'Walls dxx,dyy,dzz: ',dxx,dyy,dzz
          !          pause
          if(nparticles>1) then
             !Check for overlaps with other particles
             do jpart=ipart+1,nparticles
                rsum=radii(ipart) + radii(jpart)
                dxx=centerpx(ipart)-centerpx(jpart);dx2=dxx**2.0
                dyy=centerpy(ipart)-centerpy(jpart);dy2=dyy**2.0
                dzz=centerpz(ipart)-centerpz(jpart);dz2=dzz**2.0
                !                ewrite(0,*) 'Particles dxx :',dxx,dyy,dzz,rsum
                if (dx2+dy2+dz2.lt.rsum*rsum) then
                   mdr=(dx2+dy2+dz2)**0.5
                   deltan=(rsum - mdr)
                   daux=kn*deltan/mdr
                   fx(ipart) = fx(ipart) + daux*dxx
                   fy(ipart) = fy(ipart) + daux*dyy
                   fz(ipart) = fz(ipart) + daux*dzz
                   daux=-daux
                   fx(jpart) = fx(jpart) + daux*dxx
                   fy(jpart) = fy(jpart) + daux*dyy
                   fz(jpart) = fz(jpart) + daux*dzz
                   !                   ewrite(0,*) 'HIT ANOTHER PARTICLE!!!'
                end if
             end do
          end if
       end do

       do ipart=1, nparticles
          mass=solid_density*4*pie*(radii(ipart)**3.0)/3.0
          deltavx=fx(ipart)*timestep/mass
          deltavy=fy(ipart)*timestep/mass
          deltavz=fz(ipart)*timestep/mass
          if(i==1.and.ipart==1) then
             ewrite(0,*)'gravity data: ',mass*gravity_x,mass*gravity_y,mass*gravity_z
             ewrite(0,*)'part. data: ',mass,solid_density,fx(ipart),fy(ipart),fz(ipart)
             write(1231,*) fx(ipart)/mass,fy(ipart)/mass,fz(ipart)/mass
          end if
          centerpx(ipart)=centerpx(ipart) + (pvelx(ipart)+ deltavx/2.0)*timestep
          centerpy(ipart)=centerpy(ipart) + (pvely(ipart)+ deltavy/2.0)*timestep
          centerpz(ipart)=centerpz(ipart) + (pvelz(ipart)+ deltavz/2.0)*timestep
          pvelx(ipart) = pvelx(ipart) + deltavx
          pvely(ipart) = pvely(ipart) + deltavy
          pvelz(ipart) = pvelz(ipart) + deltavz
          !          ewrite(0,*) 'pvel: ',pvelx(ipart),pvely(ipart),pvelz(ipart)
          !          ewrite(0,*) 'centpx: ',centerpx(ipart),centerpy(ipart),centerpz(ipart)
          !          ewrite(0,*) 'deltav: ',deltavx,deltavy,deltavz,mass
       end do

    end do

  end subroutine hookforce

  function drag_coefficient(a,velocity,density,size,viscosity) result (cd)
    implicit none
    real, intent(in) :: velocity,density,size,viscosity,a
    real :: cd
    real Re

    !drag coefficient equation using Brown (2003) work
    !first calculate reynolds number
    Re=(1-a)*density*velocity*size/viscosity
    cd=(24.0/Re)*(1.0+0.15*Re**0.681)+0.407/(1+8710.0/Re)
!    ewrite(0,*) 'data3 :',a, Re,cd,viscosity,density,size,velocity
    return

  end function drag_coefficient

end module solidconfiguration
