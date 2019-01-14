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
!    but WITHOUT ANY WARRANTY; without even the implied arranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module ice_melt_interface
  use quadrature
  use elements
  use field_derivatives
  use fields
  use fields_base
  use field_options
  use fields_allocates
  use state_module
  use spud
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use state_fields_module
  use boundary_conditions
  use fields_manipulation
  use surface_integrals
  use fetools
  use vector_tools
  use FLDebug
  use integer_set_module
  use pickers_inquire
  use transform_elements
  use state_fields_module
  use global_parameters, only: PYTHON_FUNC_LEN
  use embed_python
implicit none

  private
  integer, save                   :: nnodes, dimen
  type(vector_field), save        :: surface_positions
  type(vector_field), save        :: funky_positions
  ! Normals to the iceshelf interface surface
  type(vector_field), save        :: unit_normal_vectors
  integer, dimension(:), pointer, save :: sf_nodes_point

  ! Fields and variables for the interface surface
  type(scalar_field), save        :: ice_surfaceT, ice_surfaceS! these are used to populate the bcs

  public :: melt_interface_initialisation, populate_iceshelf_boundary_conditions, melt_interface_calculate, melt_interface_boundary_condition, melt_interface_cleanup, melt_interface_allocate_surface

contains

  ! Calculation of the (ice) interface meltrate using the parameterisation described in Holland et al. 2008 doi:10.1175/2007JCLI1909.1

  !----------
  ! initialise parameters based on options
  !----------

  subroutine melt_interface_initialisation(state)

    type(state_type), intent(inout)     :: state
    type(scalar_field), pointer         :: T, S
    ! For the case when Dirichlet boundary conditions are used
    character(len=FIELD_NAME_LEN)       :: bc_type
    type(integer_set)                   :: surface_ids
    integer, dimension(:),allocatable   :: interface_surface_id
    integer, dimension(2)               :: shape_option
    type(mesh_type), pointer            :: mesh
    type(mesh_type)                     :: surface_mesh
    ! Allocated and returned by create_surface_mesh
    integer, dimension(:), pointer      :: surface_nodes
    integer, dimension(:), allocatable  :: surface_element_list
    integer                             :: i,the_node
    ! For input of the hydrostatic pressure
    character(len=PYTHON_FUNC_LEN)      :: subshelf_hydrostatic_pressure_function
    real                                :: current_time
    type(vector_field), pointer         :: positions
    type(scalar_field), pointer         :: subshelf_hydrostatic_pressure
    real :: c0, cI, L, TI, a, b, gammaT, gammaS, farfield_distance
    real                                :: T_steady, S_steady
    logical :: calculate_boundaries_T, calculate_boundaries_S
    logical :: calculate_boundary_temperature, calculate_boundary_salinity
    character(len=*), parameter         :: option_path = '/ocean_forcing/iceshelf_meltrate/Holland08'

    ewrite(1,*) 'Melt interface initialisation begins'

    call melt_interface_read_coefficients(c0=c0, cI=cI, L=L, TI=TI, a=a, b=b, gammaT=gammaT, gammaS=gammaS, farfield_distance=farfield_distance, T_steady=T_steady, S_steady=S_steady)

    calculate_boundary_temperature=have_option(trim(option_path)//'/calculate_boundaries/bc_value_temperature')
    calculate_boundary_salinity=have_option(trim(option_path)//'/calculate_boundaries/bc_value_salinity')
    if (calculate_boundary_temperature) write(3,*) "Melt interface, initialised T_steady", T_steady
    if (calculate_boundary_salinity) write(3,*) "Melt interface, initialised S_steady", S_steady

    ! bc= Dirichlet initialize T and S at the ice-ocean interface
    ! This change with bc type
    if (have_option(trim(option_path)//'/calculate_boundaries')) then
      call get_option(trim(option_path)//'/calculate_boundaries', bc_type)
      select case(bc_type)
      case('neumann')
        ewrite(3,*) 'Calculating steady Neumann boundary condition - nothing to be done here.'
      case('dirichlet')
        ewrite(3,*) 'Calculating steady Dirichlet boundary condition - nothing to be done here.'
        ! Define values of temperature and salinity at ice-ocean interface
        ! Get the surface_id of the ice-ocean interface
        ! Note this code needs working through - Neumann is definitely the best choice here!
        shape_option = option_shape(trim(option_path)//"/melt_surfaceID")
        allocate(interface_surface_id(1:shape_option(1)))
        call get_option(trim(option_path)//'/melt_surfaceID', interface_surface_id)
        call allocate(surface_ids)
        call insert(surface_ids,interface_surface_id)
        mesh => extract_mesh(state,"VelocityMesh")
        ! Obtain surface_mesh, surface_element_list, from mesh and surface_id
        call melt_surf_mesh(mesh,surface_ids,surface_mesh,surface_nodes,surface_element_list)

        T => extract_scalar_field(state,"Temperature")
        S => extract_scalar_field(state,"Salinity")
        do i=1, size(surface_nodes)
          the_node = surface_nodes(i)
          call set(T,the_node,T_steady)
          call set(S,the_node,S_steady)
        enddo
      case default
        FLAbort('Unknown boundary type for ice-ocean interface.')
      end select
    else

    endif

    ! Read the hydrostatic pressure
    if(have_option(trim(option_path)//'/subshelf_hydrostatic_pressure')) then
       call get_option(trim(option_path)//'/subshelf_hydrostatic_pressure/python', subshelf_hydrostatic_pressure_function)
       ! Get current time
         call get_option("/timestepping/current_time", current_time)
       ! Set initial condition from python function
       ! TODO: This field should be generated automatically and is not needed in the schema
       subshelf_hydrostatic_pressure => extract_scalar_field(state,"subshelf_hydrostatic_pressure")
       positions => extract_vector_field(state,"Coordinate")
       call set_from_python_function(subshelf_hydrostatic_pressure, trim(subshelf_hydrostatic_pressure_function), positions, current_time)
       ewrite(1,*) "Melt interface initialisation, found hydrostatic pressure field"
    endif

    ! Additional initialisation routines
    call melt_interface_allocate_surface(state)
    call melt_interface_calculate(state)
    ! Boundary conditions for the (ice) melt interface parameterisation
    if (have_option(trim(option_path)//'/calculate_boundaries')) then
      call melt_interface_boundary_condition(state)
    endif

    ewrite(1,*) "Melt interface initialisation end"

  end subroutine melt_interface_initialisation

  subroutine populate_iceshelf_boundary_conditions(state)
    type(state_type), intent(in)       :: state
    type(scalar_field), pointer        :: T,S
    character(len=FIELD_NAME_LEN)        :: bc_type
    integer, dimension(:), allocatable   :: surf_id
    integer, dimension(2)                :: shape_option
    type(integer_set)                    :: surface_ids
    !!
    type(mesh_type), pointer           :: surface_mesh
    type(scalar_field)                 :: scalar_surface_field

    ewrite(1,*) "-----*** Begin iceshelf BC-----"
    ! Get vector of surface ids
    shape_option=option_shape("/ocean_forcing/iceshelf_meltrate/Holland08/melt_surfaceID")
    allocate(surf_id(1:shape_option(1)))
    call get_option("/ocean_forcing/iceshelf_meltrate/Holland08/melt_surfaceID", surf_id)
    ewrite(1,*) "surf_id", surf_id
    call get_option("/ocean_forcing/iceshelf_meltrate/Holland08/calculate_boundaries", bc_type)
    call allocate(surface_ids)
    call insert(surface_ids,surf_id)
    ewrite(1,*) "set2vector(surface_ids)", set2vector(surface_ids)

    ewrite(1,*) "bc_type: ", bc_type
     ! Add boundary condition
    T => extract_scalar_field(state,"Temperature")
    S => extract_scalar_field(state,"Salinity")
    !    do i=1,node_count(T)
    !        ewrite(1,*) "line2100,i,node_val(T,i): ",i, node_val(T,i)
    !    enddo
    call add_boundary_condition(T, 'temperature_iceshelf_BC', bc_type, surf_id)
    call add_boundary_condition(S, 'salinity_iceshelf_BC', bc_type, surf_id)
    deallocate(surf_id)

    ! mesh of only the part of the surface where this b.c. applies for temperature
    call get_boundary_condition(T, 'temperature_iceshelf_BC', surface_mesh=surface_mesh)
    call allocate(scalar_surface_field, surface_mesh, name="value")
    call insert_surface_field(T, 'temperature_iceshelf_BC', scalar_surface_field)
    call deallocate(scalar_surface_field)

    ! mesh of only the part of the surface where this b.c. applies for salinity
    call get_boundary_condition(S, 'salinity_iceshelf_BC', surface_mesh=surface_mesh)
    call allocate(scalar_surface_field, surface_mesh, name="value")
    call insert_surface_field(S, 'salinity_iceshelf_BC', scalar_surface_field)
    call deallocate(scalar_surface_field)

    ewrite(1,*) "-----*** End iceshelf BC-----"
  end subroutine populate_iceshelf_boundary_conditions

  subroutine melt_interface_calculate(state)
    ! Calculate the melt rate

    type(state_type), intent(inout)     :: state
    type(scalar_field), pointer         :: Tb, Sb,MeltRate,Heat_flux,Salt_flux
    type(scalar_field), pointer         :: scalarfield
    type(vector_field), pointer         :: velocity,positions
    ! Debugging purpose pointer
    type(scalar_field), pointer         :: T_loc,S_loc,P_loc
    type(vector_field), pointer         :: V_loc,Location,Location_org
    real, dimension(:), allocatable     :: vel
    integer                             :: i
    ! Some internal variables
    real                                :: fv,fv_u,fv_k
    real                                :: T,S,P,Aa,Bb,Cc,topo,k
    real                                :: c0, cI, L, TI, a, b, gammaT, gammaS, farfield_distance
    real                                ::loc_Tb,loc_Sb,loc_meltrate,loc_heatflux,loc_saltflux, minimum_fv
    ! Aa*Sb^2+Bv*Sb+Cc
    ! Sink mesh part
    integer                             :: ele,stat,the_node
    real, dimension(:), allocatable     :: local_coord,coord
    integer, dimension(:), allocatable  :: surface_node_list
    type(scalar_field)                  :: re_pressure,re_DistanceToTop
    type(scalar_field), pointer         :: temperature, salinity, fric_vel
    character(len=*), parameter         :: option_path = '/ocean_forcing/iceshelf_meltrate/Holland08'

    
    call melt_interface_allocate_surface(state)
    ! Boundary conditions for the (ice) melt interface parameterisation
    ! if (have_option(trim(option_path)//'/calculate_boundaries')) then
    !   call melt_interface_boundary_condition(state)
    ! endif
    ! call melt_interface_initialisation(state)
    ! call melt_interface_boundary_condition(state)

    ewrite(1,*) "Melt interface calculation begins"

    call melt_interface_read_coefficients(c0=c0, cI=cI, L=L, TI=TI, a=a, b=b, gammaT=gammaT, gammaS=gammaS, farfield_distance=farfield_distance, minimum_fv=minimum_fv)

    !! All the variable under /ocean_forcing/iceshelf_meltrate/Holland08 should be in coordinate mesh
    !! coordinate mesh = continous mesh
    MeltRate => extract_scalar_field(state,"MeltRate")
    call set(MeltRate,0.0)
    Tb => extract_scalar_field(state,"Tb")
    Sb => extract_scalar_field(state,"Sb")
    ! call set(MeltRate,setnan(arg))
    call set(Tb,0.0)
    call set(Sb,-1.0)
    Heat_flux => extract_scalar_field(state,"Heat_flux")
    Salt_flux => extract_scalar_field(state,"Salt_flux")
    call set(Heat_flux,0.0)
    call set(Salt_flux,0.0)

    ! Local far field values
    T_loc => extract_scalar_field(state,"Tloc")
    S_loc => extract_scalar_field(state,"Sloc")
    P_loc => extract_scalar_field(state,"Ploc")
    V_loc => extract_vector_field(state,"Vloc")
    Location => extract_vector_field(state,"Location")
    Location_org => extract_vector_field(state,"Location_org")
    call set(T_loc,0.0)
    call set(S_loc,-1.0)
    call set(P_loc,-1.0)
    allocate(vel(V_loc%dim))
    vel = 0.0
    call set(V_loc,vel)
    call set(Location,vel)
    call set(Location_org,vel)

    positions => extract_vector_field(state,"Coordinate")

    ! Surface node list
    allocate(surface_node_list(size(sf_nodes_point)))
    ! sf_nodes is calculated in "melt_allocate_surface"
    surface_node_list=sf_nodes_point

    ! Pressure read the hydrostatic pressure, specified at schema
    scalarfield => extract_scalar_field(state,"subshelf_hydrostatic_pressure")
    call allocate(re_pressure,positions%mesh, name="RePressure")
    call remap_field(scalarfield,re_pressure,stat)

    if (have_option("/geometry/ocean_boundaries/scalar_field::DistanceToBottom")) then
        scalarfield => extract_scalar_field(state,"DistanceToBottom")
        call allocate(re_DistanceToTop,positions%mesh, name="ReDistanceToBottom")
        call remap_field(scalarfield,re_DistanceToTop,stat)
    endif

    allocate(local_coord(positions%dim+1))
    allocate(coord(positions%dim))

    temperature => extract_scalar_field(state,"Temperature")
    salinity    => extract_scalar_field(state,"Salinity")
    fric_vel => extract_scalar_field(state,"FrictionVelocity")
    velocity => extract_vector_field(state,"Velocity")

    ! Loop over the surface nodes to calculate melt rate
    do i=1,size(surface_node_list)
        the_node = surface_node_list(i)
        ! Interpolating
        coord = node_val(funky_positions,the_node)
        call picker_inquire(positions, coord, ele, local_coord,global=.false.)

        !! If sum(local_coord) is not equal to 1,
        !! we know that this coord (funky position) does not exist in the domain.
        if (sum(local_coord) .gt. 2.0) then
            ewrite(0,*) "Melt interface, funky coordinate: ", node_val(funky_positions,the_node)
            !! node_val(surface_positions,the_node) = node_val(positions,the_node)
            ewrite(0,*) "Melt interface, original element: ",ele
            ewrite(0,*) "Melt interface, sum of local_coord: ",  sum(local_coord)
            FLExit("Melt interface, your funky_positions is out of the domain. Change melt_LayerLength.")
        endif

        !Number of nodes per element for temperature
        ! Find T at a particular element given the field,element,
        !INPUT: field, element number
        !Output: real value

        call scalar_finder_ele(temperature,ele,positions%dim+1,local_coord,T)
        call scalar_finder_ele(salinity,ele,positions%dim+1,local_coord,S)
        call vector_finder_ele(velocity,ele,positions%dim+1,local_coord,vel)

        ! Pressure from the schema
        P = node_val(re_pressure,the_node)

        ! Get friction velocity based on constant y 
        ! (we're computing this in Python for now)
        call scalar_finder_ele(fric_vel,ele,positions%dim+1,local_coord,fv)
        fv = max(fv, minimum_fv)
        ! This computes fric vel based on constant y_plus
        ! fv_u = sqrt(sum(vel**2)) / 11.06
        ! fv_k = (0.09**0.25) * sqrt(k)
        ! fv = max(fv_k, fv_u, minimum_fv)

        ! constant = -7.53e-8 [C Pa^(-1)] comes from Holland and Jenkins Table 1
        ! TODO: Define as a constant explicitly, i.e. real, parameter :: topo = -7.53e-8
        topo = -7.53e-8*P

        ! Define Aa,Bb,Cc
        ! Aa*Sb**2 + Bb*Sb + Cc = 0.0
        Aa = -gammaS * fv * cI * a + a * c0 * gammaT * fv

        Bb = -gammaS * fv * L + gammaS * fv * S * cI * a
        Bb = Bb - gammaS * fv * cI * (b + topo) + gammaS * fv * cI * TI
        Bb = Bb - c0 * gammaT * fv * T + c0 * gammaT * fv * (b + topo)

        Cc = gammaS * fv * S * L + gammaS * fv * S * cI * (b + topo) + gammaS * fv * S * (-cI * TI)

        ! This could be a linear equation if Aa=0
        if (Aa .eq. 0.0) then
            loc_Sb = -Cc/Bb
        else
          ! Calculate for the 2nd oewrite(1,*) "size(surface_element_list)"rder polynomial.
          ! We have two solutions.
          loc_Sb = (-Bb + sqrt(Bb**2 - 4.0*Aa*Cc))/(2.0*Aa)
          ! loc_Sb has to be larger than 0; since the salinity in the ocean is positive definite.
          if (loc_Sb .lt. 0.0) then
            loc_Sb = (-Bb - sqrt(Bb**2 - 4.0*Aa*Cc))/(2.0*Aa)
          endif
          if (loc_Sb .lt. 0.0) then
            ewrite(0,*) "Melt interface, loc_Sb: ",  loc_Sb
            ewrite(0,*) "Melt interface, Aa: ",  Aa
            ewrite(0,*) "Melt interface, Bb: ",  Bb
            ewrite(0,*) "Melt interface, Cc: ",  Cc
            ewrite(0,*) "Melt interface, T: ",  T
            ewrite(0,*) "Melt interface, S: ",  S
            ewrite(0,*) "Melt interface, P: ",  P
            ewrite(0,*) "Melt interface, fv: ",  fv
            FLExit("Melt interface, Sb is negative. The range of Salinity is not right.")
          endif
        endif

        loc_Tb = a*loc_Sb + b + topo
        loc_meltrate = gammaS*fv*(S-loc_Sb)/loc_Sb
        !! Heat flux to the ocean
        loc_heatflux = (gammaT*fv + loc_meltrate)*(loc_Tb - T) ! or loc_meltrate*L + loc_meltrate*cI*(loc_Tb-TI)
        ! loc_heatflux = c0 * (gammaT*fv + loc_meltrate)*(loc_Tb - T) ! Added missing c0 coeff
        !! Salt flux to the ocean
        loc_saltflux = (gammaS*fv + loc_meltrate)*(loc_Sb - S)

        !! These are needed to implement BCs.
        call set(MeltRate, the_node, loc_meltrate)
        call set(Tb, the_node, loc_Tb)
        call set(Sb, the_node, loc_Sb)
        call set(Heat_flux, the_node,loc_heatflux)
        call set(Salt_flux, the_node,loc_saltflux)

        !!Far field values
        call set(T_loc,the_node,T)
        call set(S_loc,the_node,S)
        call set(P_loc,the_node,P)
        call set(V_loc,the_node,vel)
        call set(Location,the_node,node_val(funky_positions,the_node))
        call set(Location_org,the_node,node_val(positions,the_node))

    enddo

    deallocate(local_coord)
    deallocate(coord)
    call deallocate(re_pressure)
    ! TODO: check this below
    if (have_option("/geometry/ocean_boundaries/scalar_field::DistanceToBottom")) then
        call deallocate(re_DistanceToTop)
    endif

    ! Remap meltrate_p heat_flux_p salt_flux_p onto MeltRate Heat_flux Salt_flux
    ! Mappting continous mesh to discontinous mesh

    ewrite(1,*) "Melt interface calculation end"

    ! if (have_option(trim(option_path)//'/calculate_boundaries')) then
    !   call melt_interface_boundary_condition(state)
    ! endif

  end subroutine melt_interface_calculate

  subroutine melt_interface_boundary_condition(state)
    ! See preprocessor/Boundary_Conditions_From_options.F90 populate_iceshelf_boundary_conditions(states(1)) as well
    type(state_type), intent(inout)     :: state
    ! Boundary conditions
    type(scalar_field), pointer         :: TT,SS
    type(scalar_field), pointer         :: scalar_surface, scalar_surface2
    type(mesh_type), pointer            :: ice_mesh
    character(len=FIELD_NAME_LEN)       :: bc_type
    type(scalar_field)                  :: T_bc,S_bc
    type(scalar_field), pointer         :: Tflux,Sflux
    type(scalar_field), pointer         :: Tb,Sb
    integer                             :: i, the_node
    real                                :: Tz,Sz
    type(mesh_type), pointer            :: mesh
    type(mesh_type)                     :: surface_mesh
    integer, dimension(:), pointer      :: surface_nodes ! allocated and returned by create_surface_mesh
    integer, dimension(:), allocatable  :: surface_element_list
    type(integer_set)                   :: surface_ids
    integer, dimension(:),allocatable   :: surf_id
    integer, dimension(2) :: shape_option
    type(vector_field)                   :: unit_normal_vectors_vel
    type(scalar_field)                  :: heat_flux_vel,salt_flux_vel
    type(vector_field), pointer         :: velocity
    integer                             :: stat
    real, dimension(:), allocatable     :: vel
    real                                :: T_steady, S_steady
    logical :: calculate_boundary_temperature, calculate_boundary_salinity
    character(len=*), parameter         :: option_path = '/ocean_forcing/iceshelf_meltrate/Holland08'
    character(len=*), parameter         :: option_path_bc = '/ocean_forcing/iceshelf_meltrate/Holland08/calculate_boundaries'

    ewrite(1,*) "Melt interface boundary condition begins"
    call melt_interface_read_coefficients(T_steady=T_steady, S_steady=S_steady)
    calculate_boundary_temperature=have_option(trim(option_path)//'/calculate_boundaries/bc_value_temperature')
    calculate_boundary_salinity=have_option(trim(option_path)//'/calculate_boundaries/bc_value_salinity')

    !! See ./preprocessor/Boundary_Conditions_From_options.F90 populate_iceshelf_boundary_conditions(states(1)) as well
    ! Get the surface_id of the ice-ocean interface
    shape_option=option_shape(trim(option_path)//"/melt_surfaceID")
    allocate(surf_id(1:shape_option(1)))
    call get_option(trim(option_path)//'/melt_surfaceID',surf_id)
    call allocate(surface_ids)
    call insert(surface_ids,surf_id)

    TT=> extract_scalar_field(state,"Temperature")
    SS=> extract_scalar_field(state,"Salinity")
    velocity => extract_vector_field(state,"Velocity")

    if (node_count(TT) .eq. node_count(velocity)) then
      mesh => extract_mesh(state,"VelocityMesh")
    else
      mesh => extract_mesh(state,"CoordinateMesh")
    endif
    call melt_surf_mesh(mesh,surface_ids,surface_mesh,surface_nodes,surface_element_list)

    !! Remapt the contious mesh to discon, unit_normal vecotr

    call allocate(unit_normal_vectors_vel,velocity%dim,TT%mesh,"MyVelocitySurfaceMesh")

    allocate(vel(velocity%dim))
    vel = 0.0
    call set(unit_normal_vectors_vel,vel)
    deallocate(vel)

    call remap_field(unit_normal_vectors,unit_normal_vectors_vel,stat)

    call allocate(T_bc,TT%mesh, name="T_boundary")
    call allocate(S_bc,SS%mesh, name="S_boundary")

    ! Surface node list
    ! This change with bc type
    call get_option(trim(option_path_bc), bc_type)

    select case(bc_type)
    case("neumann")
      Tflux => extract_scalar_field(state,"Heat_flux")
      Sflux => extract_scalar_field(state,"Salt_flux")

      !! Remap the heat flux
      call allocate(heat_flux_vel,TT%mesh,"MyHeatFluxSurfaceMesh")
      call set(heat_flux_vel,0.0)
      call remap_field(Tflux,heat_flux_vel,stat)

      !! Remap the salt flux
      call allocate(salt_flux_vel,SS%mesh,"MySaltFluxSurfaceMesh")
      call set(salt_flux_vel,0.0)
      call remap_field(Sflux,salt_flux_vel,stat)

      ! create a surface mesh to place values onto. This is for the top surface
      call get_boundary_condition(TT, 'temperature_iceshelf_BC', surface_mesh=ice_mesh)
      call allocate(ice_surfaceT, ice_mesh, name="ice_surfaceT")
      !! Temperature
      scalar_surface => extract_surface_field(TT, 'temperature_iceshelf_BC', "value")

      do i=1,node_count(ice_surfaceT)
        the_node = surface_nodes(i)
        !                Kt = (sum(node_val(unit_normal_vectors_vel,the_node)*Kt_ar))
        !                Tz = node_val(heat_flux_vel,the_node)/(-Kt*c0)
        Tz = node_val(heat_flux_vel,the_node)
        call set(scalar_surface,i,Tz)
      enddo

      ! Salinity
      call get_boundary_condition(SS, 'salinity_iceshelf_BC', surface_mesh=ice_mesh)
      call allocate(ice_surfaceS, ice_mesh, name="ice_surfaceS")
      !! Salinity
      scalar_surface2 => extract_surface_field(SS, 'salinity_iceshelf_BC', "value")
      do i=1,node_count(ice_surfaceS)
        the_node = surface_nodes(i)
        !                Ks = abs(sum(node_val(unit_normal_vectors_vel,the_node)*Ks_ar))
        !                Sz = node_val(salt_flux_vel,the_node)/(-Ks)
        Sz = node_val(salt_flux_vel,the_node)
        !call set(ice_surfaceS,i,Sz)
        call set(scalar_surface2,i,Sz)
      enddo
      ! Deallocate the heat_flux and salt_flux on velocity mesh
      call deallocate(heat_flux_vel)
      call deallocate(salt_flux_vel)

      ! If the fluxes were constant
      if (calculate_boundary_temperature) then
        call set(T_bc,T_steady)
        ! create a surface mesh to place values onto. This is for the top surface
        call get_boundary_condition(TT, 'temperature_iceshelf_BC', surface_mesh=ice_mesh)
        call allocate(ice_surfaceT, ice_mesh, name="ice_surfaceT")
        do i=1,node_count(ice_surfaceT)
          the_node = surface_nodes(i)
          call set(ice_surfaceT,i,node_val(T_bc,the_node))
        enddo
        ! Temperature
        scalar_surface => extract_surface_field(TT, 'temperature_iceshelf_BC', "value")
        call remap_field(ice_surfaceT, scalar_surface)
      endif

      if (have_option(trim(option_path_bc)//'/bc_value_salinity')) then
        call set(S_bc,S_steady)
        ! Salinity
        call get_boundary_condition(SS, 'salinity_iceshelf_BC', surface_mesh=ice_mesh)
        call allocate(ice_surfaceS, ice_mesh, name="ice_surfaceS")
        do i=1,node_count(ice_surfaceS)
          the_node = surface_nodes(i)
          call set(ice_surfaceS,i,node_val(S_bc,the_node))
        enddo
        !! Salinity
          scalar_surface => extract_surface_field(SS, 'salinity_iceshelf_BC', "value")
          call remap_field(ice_surfaceS, scalar_surface)
      endif

    case("dirichlet")
      Tb => extract_scalar_field(state,"Tb")
      Sb => extract_scalar_field(state,"Sb")
      do i=1,node_count(T_bc)
        call set(T_bc,i,node_val(Tb,i))
        call set(S_bc,i,node_val(Sb,i))
      enddo
      ! When testing BC implementation
      if (calculate_boundary_temperature) then
        call set(T_bc,T_steady)
        ewrite(3,*) "Melt interface, initialised T_steady", T_steady
      endif
      if (calculate_boundary_salinity) then
        call set(S_bc,S_steady)
        ewrite(3,*) "Melt interface, initialised S_steady", S_steady
      endif

    case default
      FLAbort('Unknown boundary condition for ice-ocean interface')
    end select

    ! create a surface mesh to place values onto. This is for the top surface
    !call get_boundary_condition(TT, 'temperature_iceshelf_BC', surface_mesh=ice_mesh)
    !call allocate(ice_surfaceT, ice_mesh, name="ice_surfaceT")
    ! Define ice_surfaceT according to Heat_flux?

    !do i=1,node_count(ice_surfaceT)
    !    the_node = surface_node_list(i)
    !    call set(ice_surfaceT,i,node_val(T_bc,the_node))
    !enddo

    ! Salinity
    !call get_boundary_condition(SS, 'salinity_iceshelf_BC', surface_mesh=ice_mesh)
    !call allocate(ice_surfaceS, ice_mesh, name="ice_surfaceS")
    !do i=1,node_count(ice_surfaceS)
    !    the_node = surface_node_list(i)
    !    call set(ice_surfaceS,i,node_val(S_bc,the_node))
    !enddo

    !!! Temperature
    !scalar_surface => extract_surface_field(TT, 'temperature_iceshelf_BC', "value")
    !call remap_field(ice_surfaceT, scalar_surface)
    !!! Salinity
    !scalar_surface => extract_surface_field(SS, 'salinity_iceshelf_BC', "value")
    !call remap_field(ice_surfaceS, scalar_surface)

    call deallocate(T_bc)
    call deallocate(S_bc)
    call deallocate(ice_surfaceT)
    call deallocate(ice_surfaceS)
    call deallocate(unit_normal_vectors_vel)
    call deallocate(surface_mesh)

    ewrite(1,*) "Melt interface boundary condition end"

  end subroutine melt_interface_boundary_condition


  subroutine melt_interface_cleanup()
    ewrite(1,*) "Melt interface, cleaning up melt variables"
    call deallocate(surface_positions)
    call deallocate(funky_positions)
    call deallocate(unit_normal_vectors)

    deallocate(sf_nodes_point)
  end subroutine melt_interface_cleanup


  !------------------------------------------------------------------!
  !                       Private subroutines                        !
  !------------------------------------------------------------------!

  subroutine melt_interface_allocate_surface(state)
    type(state_type), intent(inout)     :: state
    type(vector_field), pointer         :: positions
    type(mesh_type), pointer            :: mesh
    type(mesh_type)                     :: surface_mesh
    ! Allocated and returned by create_surface_mesh
    integer, dimension(:), pointer      :: surface_nodes
    integer, dimension(:), allocatable  :: surface_element_list
    integer    :: face,i,j,k
    real, dimension(:), allocatable     :: coord
    ! For transform_facet_to_physical
    real, dimension(:,:), allocatable   :: normal
    real, dimension(:), allocatable     :: av_normal !Average of normal
    type(element_type), pointer     :: x_shape_f
    real, dimension(:), allocatable     :: xyz !New location of surface_mesh, farfield_distance away from the boundary
    type(integer_set)                   :: surface_ids
    integer, dimension(:),allocatable   :: interface_surface_id

    ! For the vector averaging schem, taking the area of the surface element etc
    real, dimension(:,:), allocatable    :: table
    integer, dimension(:,:), allocatable :: node_occupants
    integer                              :: st,en,node,dim_vec
    real                                 :: area_sum
    integer, dimension(2) :: shape_option
    real, dimension(:), allocatable     :: local_coord
    integer                             :: ele, counter,Nit
    real, dimension(:), allocatable     :: vel
    real :: c0, cI, L, TI, a, b, gammaT, gammaS, farfield_distance
    character(len=*), parameter         :: option_path = '/ocean_forcing/iceshelf_meltrate/Holland08/melt_surfaceID'

    ewrite(1,*) "Melt interface allocation of surface parameters"

    call melt_interface_read_coefficients(c0=c0, cI=cI, L=L, TI=TI, a=a, b=b, gammaT=gammaT, gammaS=gammaS, farfield_distance=farfield_distance)

    ! Get the surface_id of the ice-ocean interface
    shape_option=option_shape(trim(option_path))
    allocate(interface_surface_id(1:shape_option(1)))
    call get_option(trim(option_path), interface_surface_id)
    call allocate(surface_ids)
    call insert(surface_ids,interface_surface_id)

    mesh => extract_mesh(state,"CoordinateMesh")
    ! Input, mesh, surface_id
    ! Output surface_mesh, surface_element_list
    call melt_surf_mesh(mesh, surface_ids, surface_mesh, surface_nodes, surface_element_list)

    allocate(sf_nodes_point(size(surface_nodes)))
    sf_nodes_point = surface_nodes
    positions => extract_vector_field(state,"Coordinate")
    ! Sink the surface_mesh to the 'far field'
    call allocate(surface_positions,positions%dim,mesh,"MySurfacePosition")
    call allocate(funky_positions,surface_positions%dim, surface_positions%mesh,"MyFunkyPosition")
    call allocate(unit_normal_vectors,surface_positions%dim, surface_positions%mesh,"MyFunkyUnitVec")
    allocate(vel(positions%dim))
    vel = 0.0
    call set(unit_normal_vectors,vel)
    deallocate(vel)
    ! Check above

    ! Remap the positions vector to the surface mesh, surface_positions
    call remap_field_to_surface(positions, surface_positions, surface_element_list)

    ! Loop over # of surface elements
    !nf = size(surface_element_list)
    !2 comes because there are 2 nodes per a surface_element (assumption!)
    ! this number could be three for 3D problem, room for change
    ! dim_vec is the dimension of the positions, either 2 or 3.
    ! node_occupants(1,1:dim_vec*nf) = nodes that occupies the surface elements
    ! node_occupants(2,:) = surface elements
    ! table(1,:) = area of the face
    ! table(2,:) = normal_x
    ! table(3,:) = normal_y
    ! table(4,:) = normal_z
    dim_vec = positions%dim
    allocate(coord(dim_vec))
    allocate(node_occupants(2,dim_vec*size(surface_element_list)))
    allocate(table(dim_vec+1,dim_vec*size(surface_element_list)))
    allocate(av_normal(dim_vec))
    allocate(xyz(dim_vec))
    allocate(local_coord(positions%dim+1))

    do i=1,size(surface_element_list)
        st = 1+dim_vec*(i-1)
        en = dim_vec+dim_vec*(i-1)
        face = surface_element_list(i)
        node_occupants(2,st:en) = face
        node_occupants(1,st:en) = face_global_nodes(positions, face)
        ! Calculate the area of the surface element
        ! For 2D, the area is the length
        if (dim_vec .eq. 2) then
            table(1,st:en) = calc_area2(positions,node_occupants(1,st:en))
        endif

        ! For 3D, the area become the area of the triangle (assumption here).
        ! Subroutine calc_area3 uses Heron's formula.
        if (dim_vec .eq. 3) then
            table(1,st:en) = calc_area3(positions,node_occupants(1,st:en))
        endif
        ! Compute the normal vector at the surface element.
        x_shape_f=>face_shape(positions,face)
        allocate(normal(dim_vec,x_shape_f%ngi)) !(dim x x_shape_f%ngi)
        call transform_facet_to_physical(positions, face, normal=normal)
        av_normal = sum(normal,2)
        !"normal" should be the same; since the normal vector on the quadrature point does not change.
        av_normal = av_normal/(sum(av_normal*av_normal))**(0.5) ! normalize the vector
        deallocate(normal)
        ! Since av_normal is outward to the boundary, we will put -sign. We want the vector into the boundary
        av_normal = -av_normal
        ! Store av_normal, the normal vector averaged over quadrature points, in table
        do j=1,dim_vec
            table(j+1,st:en)=av_normal(j)
        enddo
    enddo
    ! Now loop over surface_nodes
    !        ewrite(1,*) "table(1,1:3): ", table(1,1:3)
    !        ewrite(1,*) "table(2,1:3),normal_x: ", table(2,1:3)
    !        ewrite(1,*) "table(3,1:3),normal_y: ", table(3,1:3)
    !        ewrite(1,*) "table(4,1:3),normal_z: ", table(4,1:3)
    ! In this loop, we will average adjacent normal vectors, using the area of the surface elements.


    do i=1,size(node_occupants(1,:))
        node = node_occupants(1,i)
        av_normal(:) = 0.0
        area_sum = 0.0
        ! This loop average the normal vector based on the areas of surface elements
        ! which share the "node". Using hash table may speed up this process
        do j=1,size(node_occupants(1,:))
            !Pick the surface elements that sheare the "node"
            if (node_occupants(1,j) .eq. node) then
                do k=1,dim_vec
                    !table(1,j) = area of the surface element
                    av_normal(k) = av_normal(k) + table(k+1,j)*table(1,j)
                enddo
                !The total areas of the surface elements that occupies the "node"
                area_sum = area_sum + table(1,j)
            endif
        enddo

        av_normal = av_normal / area_sum
        ! normalize
        av_normal = av_normal/(sum(av_normal*av_normal))**(0.5)

        farfield_distance = abs(farfield_distance)
        ! Shift the location of the surface nodes.
        !The coordinate of the surface node.

        coord = node_val(positions,node)
        ! farfield_distance = ||xyz - coord|| xyz and coord are vectors
        xyz = av_normal*farfield_distance + coord

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! have options /ocean_forcing/iceshelf_meltrate/Holland08/melt_LayerRelax
        if (have_option("/ocean_forcing/iceshelf_meltrate/Holland08/melt_LayerRelax")) then
          call get_option('/ocean_forcing/iceshelf_meltrate/Holland08/melt_LayerRelax/N',Nit)
          !!Check if xyz is within the domain
          call picker_inquire(positions, xyz, ele, local_coord,global=.false.)
          !! If sum(local_coord) is not equal to 1,
          !! we know that this coord (funky position) does not exist in the domain.

          counter = 1
          do while (sum(local_coord) .gt. 2.0 .and. counter .lt. Nit)
            xyz = av_normal*(farfield_distance-counter*(farfield_distance/Nit)) + coord
            call picker_inquire(positions, xyz, ele, local_coord,global=.false.)
            counter = counter+1
          enddo

          if (sum(local_coord) .gt. 2.0) then
            ewrite(1,*) "icehslef location, coord: ", coord
            ewrite(1,*) "funky position, xyz: ", xyz
            call picker_inquire(positions, coord, ele, local_coord,global=.false.)
            !! node_val(surface_positions,the_node) = node_val(positions,the_node)
            ewrite(1,*) "Original ele: ",ele
            ewrite(1,*) "sum of local_coord: ",  sum(local_coord)
            FLExit("In setting funky_position, your funky_positions is out of the domain. Change melt_LayerLength.")
          endif
        endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call set(funky_positions,node,xyz) ! Set the coordinate of sinked nodes, funky positions.
        call set(surface_positions,node,coord) !Original coordinate of the surface
        call set(unit_normal_vectors,node,av_normal) !save the unit normal vector for Kt and Ks calculation

        node = 0
    enddo

    ! Now deallocate junk
    deallocate(coord)
    deallocate(table)
    deallocate(node_occupants)
    deallocate(av_normal)
    deallocate(xyz)
    deallocate(local_coord)
    call deallocate(surface_ids)
    call deallocate(surface_mesh)
    ewrite(1,*) "Melt interface allocation of surface parameters end"

  end subroutine melt_interface_allocate_surface

  subroutine melt_surf_mesh(mesh,surface_ids,surface_mesh,surface_nodes,surface_element_list)

    type(mesh_type), intent(in)                       :: mesh
    type(integer_set), intent(in)                     :: surface_ids
    type(mesh_type), intent(out)                      :: surface_mesh
    integer, dimension(:), pointer, intent(out)       :: surface_nodes ! allocated and returned by create_surface_mesh
    integer, dimension(:), allocatable, intent(out)   :: surface_element_list
    !    integer, dimension(:), allocatable                :: surface_nodes_out
    type(integer_set)                                 :: surface_elements
    integer :: i
    ! create a set of surface elements that have surface id in 'surface_ids'

    call allocate(surface_elements)
    do i=1, surface_element_count(mesh)
      ! ewrite(1,*) "surf_normal surface_element_id(mesh, i)", surface_element_id(mesh, i)
      ! ewrite(1,*) "surf_normal surface_ids", set2vector(surface_ids)
      if (has_value(surface_ids, surface_element_id(mesh, i))) then
        call insert(surface_elements, i)
      end if
    end do

    allocate(surface_element_list(key_count(surface_elements)))
    surface_element_list=set2vector(surface_elements)
    call create_surface_mesh(surface_mesh, surface_nodes, mesh, surface_element_list, name=trim(mesh%name)//"ToshisMesh")
  end subroutine melt_surf_mesh

  subroutine melt_interface_read_coefficients(c0, cI, L, TI, a, b, gammaT, gammaS, farfield_distance, T_steady, S_steady, minimum_fv)
    real, intent(out), optional :: c0, cI, L, TI, a, b, gammaT, gammaS, farfield_distance, T_steady, S_steady, minimum_fv
    real :: Cd
    character(len=*), parameter :: option_path = '/ocean_forcing/iceshelf_meltrate/Holland08'

    ! Get the 6 model constants
    ! TODO: Check these exist first and error with a useful message if not - in preprocessor
    if (present(c0)) call get_option(trim(option_path)//'c0', c0, default = 3974.0)
    if (present(cI)) call get_option(trim(option_path)//'cI', cI, default = 2009.0)
    if (present(L))  call get_option(trim(option_path)//'/L', L, default = 3.35e5)
    if (present(TI)) call get_option(trim(option_path)//'/TI', TI, default = -25.0)
    if (present(a))  call get_option(trim(option_path)//'/a', a, default = -0.0573)
    if (present(b))  call get_option(trim(option_path)//'/b', b, default = 0.0832)
    if (present(farfield_distance)) call get_option(trim(option_path)//'/melt_LayerLength', farfield_distance)

    ! if (present(gammaT).or.present(gammaS)) call get_option(trim(option_path)//'/Cd', Cd, default = 1.5e-3)
    ! if (present(gammaT)) gammaT = sqrt(Cd)/(12.5*(7.0**(2.0/3.0))-9.0)
    ! if (present(gammaS)) gammaS = sqrt(Cd)/(12.5*(700.0**(2.0/3.0))-9.0)

    ! Use suggested values for gammaT, gammaS from Jenkins, Nicholls, Corr, 2010
    if (present(gammaT)) gammaT = 1.10e-2
    if (present(gammaS)) gammaS = 3.10e-4


    ! When steady boundary options are enabled
    if (present(T_steady)) call get_option(trim(option_path)//'/calculate_boundaries/bc_value_temperature',T_steady, default = 0.0)
    if (present(S_steady)) call get_option(trim(option_path)//'/calculate_boundaries/bc_value_salinity',S_steady, default = 34.0)
    
    if (present(minimum_fv)) call get_option(trim(option_path)//'/minimum_friction_velocity',minimum_fv, default = 1.0E-3)

  end subroutine melt_interface_read_coefficients

  subroutine scalar_finder_ele(scalar,ele,element_dim,local_coord,scalar_out)
    type(scalar_field), pointer, intent(in)           :: scalar
    integer, intent(in)                      :: ele
    integer, intent(in)                       :: element_dim
    real, dimension(element_dim), intent(in) :: local_coord
    real, intent(out)                        :: scalar_out
    integer                        :: j,node
    integer, dimension(:), allocatable  :: node_lists

    allocate(node_lists(size(ele_nodes(scalar,ele))))
    node_lists =  ele_nodes(scalar,ele) !Lists of nodes that occupies the element.
    scalar_out = 0.0
    do j=1,size(node_lists)
        node = node_lists(j)
        scalar_out = scalar_out + node_val(scalar,node)*local_coord(j)
    enddo
  end subroutine scalar_finder_ele

  subroutine vector_finder_ele(vector,ele,element_dim,local_coord,vector_out)
    type(vector_field), pointer, intent(in)           :: vector
    integer, intent(in)                      :: ele
    integer, intent(in)                      :: element_dim
    real, dimension(element_dim), intent(in) :: local_coord
    real, dimension(element_dim-1), intent(out) :: vector_out
    integer                        :: j,node
    integer, dimension(:), allocatable  :: node_lists

    allocate(node_lists(size(ele_nodes(vector,ele))))

    node_lists =  ele_nodes(vector,ele) !Lists of nodes that occupies the element.
    vector_out(:) = 0.0
    do j=1,size(node_lists)
        node = node_lists(j)
        vector_out = vector_out + node_val(vector,node)*local_coord(j)
    enddo
  end subroutine vector_finder_ele

  real function setnan(arg)
    real :: arg
    setnan = sqrt(arg)
  end function

  real function calc_area2(positions,nodes)
    type(vector_field), pointer        :: positions
    integer, dimension(2)              :: nodes
    real, dimension(2)                 :: coord1,coord2
    coord1 = node_val(positions,nodes(1))
    coord2 = node_val(positions,nodes(2))
    calc_area2 = (sum((coord2 - coord1)**2))**0.5
  end function

  real function calc_area3(positions,nodes)
    type(vector_field), pointer        :: positions
    integer, dimension(3)              :: nodes
    real, dimension(3)                 :: coord1,coord2,coord3
    real                               :: a,b,c,s
    coord1 = node_val(positions,nodes(1))
    coord2 = node_val(positions,nodes(2))
    coord3 = node_val(positions,nodes(3))
    ! Use Heron's formula to calculate the area of triangle
    a = (sum((coord2 - coord1)**2))**0.5
    b = (sum((coord3 - coord1)**2))**0.5
    c = (sum((coord2 - coord3)**2))**0.5
    s = 0.5*(a+b+c)
    calc_area3 = (s*(s-a)*(s-b)*(s-c))**0.5
  end function

end module ice_melt_interface
