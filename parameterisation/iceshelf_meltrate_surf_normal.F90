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

module iceshelf_meltrate_surf_normal

  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
  use fldebug
  use vector_tools
  use quadrature
  use elements
  use spud
  use integer_set_module
  use fields_allocates
  use fields_manipulation
  use transform_elements
  use fetools
  use fields
  use state_module
  use boundary_conditions
  use field_derivatives
  use field_options
  use state_fields_module
  use surface_integrals
  use pickers_inquire
  use state_fields_module

implicit none

  private
  real, save                      :: c0, cI, L, TI, &
                                     a, b, gammaT, gammaS,Cd, dist_meltrate
  integer, save                   :: nnodes,dimen
 
  type(vector_field), save        :: surface_positions
  type(vector_field), save        :: funky_positions 
  type(integer_set), save         :: sf_nodes !Nodes at the surface
  !BC
  integer, dimension(:), pointer, save :: ice_element_list
  ! these are the fields and variables for the surface values
  type(scalar_field), save             :: ice_surfaceT,ice_surfaceS! these are used to populate the bcs
  public :: melt_surf_init, melt_allocate_surface, melt_surf_calc, melt_bc



contains

!----------
! initialise parameters based on options
!----------

 subroutine melt_surf_init(state)

    type(state_type), intent(inout)     :: state
    character(len=OPTION_PATH_LEN)      :: melt_path
    ! hack
    type(scalar_field), pointer                  :: T,S
    ! When bc=Dirichlet 
     character(len=FIELD_NAME_LEN)       :: bc_type
    type(integer_set)                   :: surface_ids
    integer, dimension(:),allocatable   :: surf_id
    integer, dimension(:), allocatable  :: sf_nodes_ar
    integer, dimension(2) :: shape_option
    type(mesh_type), pointer            :: mesh
    type(mesh_type)                     :: surface_mesh
    integer, dimension(:), pointer      :: surface_nodes ! allocated and returned by create_surface_mesh
    integer, dimension(:), allocatable  :: surface_element_list
     integer                             :: i,the_node
      
    melt_path = "/ocean_forcing/iceshelf_meltrate/Holland08"

 
    ewrite(1,*) "--------Begin melt_init-------------"

    ! Get the 6 model constants
    call get_option(trim(melt_path)//'/c0', c0, default = 3974.0)
    call get_option(trim(melt_path)//'/cI', cI, default = 2009.0)
    call get_option(trim(melt_path)//'/L', L, default = 3.35e5)
    call get_option(trim(melt_path)//'/TI', TI, default = -25.0)
    call get_option(trim(melt_path)//'/a', a, default = -0.0573)
    call get_option(trim(melt_path)//'/b', b, default = 0.0832)
    call get_option(trim(melt_path)//'/Cd', Cd, default = 1.5e-3)
    call get_option(trim(melt_path)//'/melt_LayerLength', dist_meltrate)
 
    gammaT = sqrt(Cd)/(12.5*(7.0**(2.0/3.0))-9.0)
    gammaS = sqrt(Cd)/(12.5*(700.0**(2.0/3.0))-9.0)

    !! bc= Dirichlet initialize T and S at the ice-ocean interface
    !hack
    !This change with bc type
    call get_option("/ocean_forcing/iceshelf_meltrate/Holland08/calculate_boundaries", bc_type)   
       
        select case(bc_type)
        case("neumann")
         
            
        case("dirichlet") 
        !! Define the values at ice-ocean interface    
        ! Get the surface_id of the ice-ocean interface
            shape_option=option_shape(trim(melt_path)//"/melt_surfaceID")
            allocate(surf_id(1:shape_option(1)))
            call get_option(trim(melt_path)//'/melt_surfaceID',surf_id)
            call allocate(surface_ids)
            call insert(surface_ids,surf_id)
            mesh => extract_mesh(state,"VelocityMesh")
            ! Input, mesh, surface_id
            ! Output surface_mesh,surface_element_list
            call melt_surf_mesh(mesh,surface_ids,surface_mesh,surface_nodes,surface_element_list)

            T => extract_scalar_field(state,"Temperature")
            S => extract_scalar_field(state,"Salinity")
            do i=1,size(surface_nodes)
                the_node = surface_nodes(i)
                
                call set(T,the_node,0.0)
                call set(S,the_node,34.0)
               
            enddo
!            T_bc => extract_scalar_field(state,"Tb")
!            S_bc => extract_scalar_field(state,"Sb")
        case default
            FLAbort('Unknown BC for TKE')
        end select 
    

    ewrite(1,*) "---------End melt_surf_init---------------------------------"

 end subroutine melt_surf_init

subroutine melt_allocate_surface(state)

    type(state_type), intent(inout)     :: state
    type(vector_field), pointer         :: positions   
    type(mesh_type), pointer            :: mesh
    type(mesh_type)                     :: surface_mesh
    type(integer_set)                   :: surface_elements
    integer, dimension(:), pointer      :: surface_nodes ! allocated and returned by create_surface_mesh
    integer, dimension(:), allocatable  :: surface_element_list
    integer    :: ele, face,i,j,k
!! The local coordinates of the coordinate in the owning element
!!    real, dimension(:), allocatable    :: local_coord, local_coord_surf
    real, dimension(:), allocatable     :: coord
! for transform_facet_to_physical
    real, dimension(:,:), allocatable   :: normal
    real, dimension(:), allocatable     :: av_normal !Average of normal
!!    integer, dimension(:), pointer      :: ele_faces,face_neighs
    type(element_type), pointer     :: x_shape_f
    real                                :: al,be,ga
    real, dimension(:), allocatable     :: xyz !New location of surface_mesh, dist_meltrate away from the boundary
!integer, save         :: melt_surfaceID
    character(len=OPTION_PATH_LEN)      :: path
    type(integer_set)                   :: surface_ids
    integer, dimension(:),allocatable   :: surf_id

! For the vector averaging schem, taking the area of the surface element etc 
    real, dimension(:,:), allocatable   :: table
    integer, dimension(:,:), allocatable :: node_occupants
    integer                             :: st,en,node,dim_vec
    real                                :: area_sum
    integer, dimension(:), allocatable  :: sf_nodes_ar
    integer, dimension(2) :: shape_option

    ewrite(1,*) "-------Begin melt_allocate_surface---------"
    path = "/ocean_forcing/iceshelf_meltrate/Holland08"
    ! Get the surface_id of the ice-ocean interface
    shape_option=option_shape(trim(path)//"/melt_surfaceID")
    allocate(surf_id(1:shape_option(1)))
    call get_option(trim(path)//'/melt_surfaceID',surf_id)
    call allocate(surface_ids)
    call insert(surface_ids,surf_id)
!    deallocate(surf_id)
!    ewrite(1,*) "aft_surface_ids: ", set2vector(surface_ids)

    mesh => extract_mesh(state,"VelocityMesh")
! Input, mesh, surface_id
! Output surface_mesh,surface_element_list
    call melt_surf_mesh(mesh,surface_ids,surface_mesh,surface_nodes,surface_element_list)

    positions => extract_vector_field(state,"Coordinate")
!SINK THE surface_mesh!!
!   call get_option("/geometry/dimension/", dimension)
    call allocate(surface_positions,positions%dim,mesh,"MySurfacePosition")
    call allocate(funky_positions,surface_positions%dim, surface_positions%mesh,"MyFunkyPosition")
    

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
        !Since av_normal is outward to the boundary, we will put -sign. We want the vector into the boundary
        av_normal = -av_normal
        ! Store av_normal, the normal vector averaged over quadrature points, in table        
        do j=1,dim_vec
            table(j+1,st:en)=av_normal(j)
        enddo
    enddo
!! Now loop over surface_nodes   
!        ewrite(1,*) "table(1,1:3): ", table(1,1:3)
!        ewrite(1,*) "table(2,1:3),normal_x: ", table(2,1:3)
!        ewrite(1,*) "table(3,1:3),normal_y: ", table(3,1:3)
!        ewrite(1,*) "table(4,1:3),normal_z: ", table(4,1:3)
!!In this loop, we will average adjacent normal vectors, using the area of the surface elements.
    
    !allocate(sf_nodes_ar(size(node_occupants(1,:))))
    call allocate(sf_nodes)
!    call insert(sf_nodes,sf_nodes_ar)
    do i=1,size(node_occupants(1,:))
        node = node_occupants(1,i)
        av_normal(:) = 0.0
        area_sum = 0.0
        ! This loop average the normal vector based on the areas of surface elements
        !, which share the "node". Using hash table may speed up this process
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
         
        dist_meltrate = abs(dist_meltrate)
        ! Shift the location of the surface nodes.
        !The coordinate of the surface node.
        
        coord = node_val(positions,node) 
        
        ! dist_meltrate = ||xyz - coord|| xyz and coord are vectors
        xyz = av_normal*dist_meltrate + coord
        
        call set(funky_positions,node,xyz) ! Set the coordinate of sinked nodes, funky positions.
       
        call set(surface_positions,node,coord) !Original coordinate of the surface
!!        ewrite(1,*) "--------------------------------"
!!        ewrite(1,*) "node: ", node
!!        ewrite(1,*) "av_normal: ", av_normal
!!        ewrite(1,*) "funky: ", xyz
!!        ewrite(1,*) "coord: ", coord
!!        ! Save the corresponding node number.
!!        !sf_nodes_ar(i) = node
        call insert(sf_nodes,node)
       
        node = 0
    enddo
     

    deallocate(coord)
    deallocate(table)
    deallocate(node_occupants)    
    deallocate(av_normal)
    deallocate(xyz)
    call deallocate(surface_ids)
     ewrite(1,*) "-------End melt_allocate_surface---------"
end subroutine melt_allocate_surface



! Calculate the melt rate
 subroutine melt_surf_calc(state)

    type(state_type), intent(inout)     :: state
    type(scalar_field), pointer         :: Tb, Sb,MeltRate,Heat_flux,Salt_flux
    type(scalar_field), pointer         :: scalarfield
    type(vector_field), pointer         :: velocity,positions  
    !Debugging purpose pointer
    type(scalar_field), pointer         :: T_loc,S_loc,P_loc
    type(vector_field), pointer         :: V_loc,Location,Location_org 
    real, dimension(:), allocatable     :: vel
    integer                             :: i,j,k
    ! Some internal variables  
    real                                :: speed, T,S,P,Aa,Bb,Cc,topo
    real                                ::loc_Tb,loc_Sb,loc_meltrate,loc_heatflux,loc_saltflux
    ! Aa*Sb^2+Bv*Sb+Cc
    real :: arg = -1.0
    !Sink mesh part
    integer                             :: ele,face,node,stat,the_node
    real, dimension(:), allocatable     :: local_coord,coord
    type(element_type), pointer         :: x_shape_f
    integer, dimension(:), allocatable  :: surface_node_list,node_lists 
    type(scalar_field)                  :: re_temperature,re_salinity,re_pressure
    type(vector_field)                  :: re_velocity
  
    
    ewrite(1,*) "-------Begin melt_surf_calc------------"

    MeltRate => extract_scalar_field(state,"MeltRate")
    Tb => extract_scalar_field(state,"Tb")
    Sb => extract_scalar_field(state,"Sb")
    call set(MeltRate,setnan(arg))
    call set(Tb,setnan(arg))
    call set(Sb,setnan(arg))
    Heat_flux => extract_scalar_field(state,"Heat_flux")
    Salt_flux => extract_scalar_field(state,"Salt_flux")
    call set(Heat_flux,setnan(arg))
    call set(Salt_flux,setnan(arg))
    
! Debugging
    T_loc => extract_scalar_field(state,"Tloc")
    S_loc => extract_scalar_field(state,"Sloc")
    P_loc => extract_scalar_field(state,"Ploc")
    V_loc => extract_vector_field(state,"Vloc")
    Location => extract_vector_field(state,"Location")
    Location_org => extract_vector_field(state,"Location_org")
    call set(T_loc,setnan(arg))
    call set(S_loc,setnan(arg))
    call set(P_loc,setnan(arg))
    allocate(vel(V_loc%dim))
    vel = setnan(arg)
    call set(V_loc,vel)
    call set(Location,vel)
    call set(Location_org,vel)
    
    positions => extract_vector_field(state,"Coordinate")
    !my positions
  
      ! Surface node list
    allocate(surface_node_list(key_count(sf_nodes)))
    ! Make it to vector from integer_set.
    ! sf_nodes is calculated in "melt_allocate_surface"
    surface_node_list=set2vector(sf_nodes)
    
    ! Remap temperature, salinity, pressure, and velocity onto positions mesh
    scalarfield => extract_scalar_field(state,"Temperature")
    
    call allocate(re_temperature,positions%mesh,name="ReTemperature")
    
    call remap_field(scalarfield,re_temperature,stat)
   
    ! Salinity
    scalarfield=> extract_scalar_field(state,"Salinity")
 
    call allocate(re_salinity,positions%mesh, name="ReSalinity")
    call remap_field(scalarfield,re_salinity,stat)

    ! Pressure
    scalarfield => extract_scalar_field(state,"Pressure")
    call allocate(re_pressure,positions%mesh, name="RePressure")
    call remap_field(scalarfield,re_pressure,stat)
   
    ! Velocity    
    velocity => extract_vector_field(state,"Velocity")
    call allocate(re_velocity,velocity%dim,positions%mesh,name="ReVelocity")
    call remap_field(velocity, re_velocity, stat)
 
    allocate(local_coord(positions%dim+1))
    allocate(coord(positions%dim))
   
  


!!!!!!!!!!!!!!!!!!!!!!!!!
     !surface_positions

    !! Loope over the surface nodes to calculate melt rate etc.
    do i=1,size(surface_node_list)
        the_node = surface_node_list(i)
        !!! Interpolating   
        coord = node_val(funky_positions,the_node)
        call picker_inquire(positions, coord, ele, local_coord,global=.true.)

        !! If sum(local_coord) is not equal to 1,      
        !! we know that this coord (funky position) does not exist in the domain.
        
        if (sum(local_coord) .gt. 2.0) then 
            ewrite(1,*) "Funk coord: ", node_val(funky_positions,the_node)
            !! node_val(surface_positions,the_node) = node_val(positions,the_node)
            ewrite(1,*) "Original ele: ",ele   
            ewrite(1,*) "sum of local_coord: ",  sum(local_coord)    
            FLExit("Your funky_positions is out of the domain. Change melt_LayerLength.")
        endif
        !Number of nodes per element       
        allocate(node_lists(ele_and_faces_loc(positions, ele))) 
        node_lists =  ele_nodes(positions,ele) !Lists of nodes that occupies the element.
        ! This method of finding values at funky_positions works for P1, which are T,S, and velocity.
        coord = 0.0
        T = 0.0
        S = 0.0
        speed = 0.0
        vel = 0.0

        do j=1,size(node_lists)
            node = node_lists(j)
            coord = coord + node_val(positions,node)*local_coord(j)
            T = T + node_val(re_temperature,node)*local_coord(j)
            S = S + node_val(re_salinity,node)*local_coord(j)
            vel = vel + node_val(re_velocity,node)*local_coord(j)
        enddo
        deallocate(node_lists)
        ! Luckly P needs to be at the surface for the three equations
        P = node_val(re_pressure,the_node)
        speed = sqrt(sum(vel**2))
        
        if (speed .lt. 0.001) then
            speed = 0.001
             !ewrite(1,*) "----------iceshelf, speed less----", the_node,speed
        endif
        topo = -7.53e-8*P ! constant = -7.53e-8 [C Pa^(-1)] comes from Holland and Jenkins Table 1
        
        !! Define Aa,Bb,Cc
        !! Aa*Sb**2 + Bb*Sb + Cc = 0.0
        Aa = -gammaS*speed*cI*a + a*c0*gammaT*speed
        
        Bb = -gammaS*speed*L + gammaS*speed*S*cI*a
        Bb = Bb - gammaS*speed*cI*(b+topo) + gammaS*speed*cI*TI
        Bb = Bb - c0*gammaT*speed*T + c0*gammaT*speed*(b+topo)

        Cc = gammaS*speed*S*L +gammaS*speed*S*cI*(b+topo) + gammaS*speed*S*(-cI*TI)
        
        
        !! This could be a linear equation if Aa=0
        if (Aa .eq. 0.0) then
            loc_Sb = -Cc/Bb
        else
            !! Calculate for the 2nd oewrite(1,*) "size(surface_element_list)"rder polynomial.
            !! We have two solutions. 
            loc_Sb = (-Bb + sqrt(Bb**2 - 4.0*Aa*Cc))/(2.0*Aa)
            !! loc_Sb has to be larger than 0; since the salinity in the ocean is positive definite.
            if (loc_Sb .lt. 0.0) then
                loc_Sb = (-Bb - sqrt(Bb**2 - 4.0*Aa*Cc))/(2.0*Aa)
            endif
        endif
        !ewrite(1,*) "----------iceshelf, loc_Sb----", loc_Sb
        loc_Tb = a*loc_Sb + b + topo
        loc_meltrate = gammaS*speed*(S-loc_Sb)/loc_Sb
        !! Heat flux to the ocean
        loc_heatflux = c0*(gammaT*speed+ loc_meltrate)*(T-loc_Tb) ! or loc_meltrate*L + loc_meltrate*cI*(loc_Tb-TI)
        loc_saltflux = (gammaS*speed+loc_meltrate)*(S-loc_Sb)
        !! Some debugging
        !!ewrite(1,*) "melt_rate: ",loc_meltrate
        !!ewrite(1,*) "tLHS: ", c0*gammaT*speed*(T-loc_Tb)
        !!ewrite(1,*) "tRHS: ", loc_meltrate*L + loc_meltrate*cI*(loc_Tb-TI)
        !!ewrite(1,*) "sLHS: ", gammaS*speed*(S-loc_Sb)
        !!ewrite(1,*) "sRHS: ", loc_meltrate*loc_Sb
        !!ewrite(1,*) "bLHS: ", loc_Tb
        !!ewrite(1,*) "bRHS: ", a*loc_Sb + b + topo

        !! These are needed to implement BCs.
        call set(MeltRate, the_node, loc_meltrate)
        call set(Tb, the_node, loc_Tb)
        call set(Sb, the_node, loc_Sb)
        call set(Heat_flux, the_node,loc_heatflux)
        call set(Salt_flux, the_node,loc_saltflux)
        !! More or less for debugging purposes.
        call set(T_loc,the_node,T)
        call set(S_loc,the_node,S)
        call set(P_loc,the_node,P)
        call set(V_loc,the_node,vel)
        call set(Location,the_node,node_val(funky_positions,the_node))
        call set(Location_org,the_node,node_val(positions,the_node))
        !! BC test
        !!call set(TT,the_node,11.1)
!         ewrite(1,*) "----------iceshelf, loc_saltflux----", the_node,loc_saltflux
!        ewrite(1,*) "----------melt_surf_calc, loc_Tb----", the_node,loc_Tb
!        ewrite(1,*) "----------iceshelf, surfaceS----",node_val(re_salinity,the_node)
!        ewrite(1,*) "----------line 477 iceshelf, loc_meltrate----",loc_meltrate
!        ewrite(1,*) "----------iceshelf, node_val(TT,the_node)----",node_val(TT,the_node)
    enddo

    deallocate(local_coord)
    deallocate(coord)
    ewrite(1,*) "-----END melt_surf_calc-------"
end subroutine melt_surf_calc


subroutine melt_bc(state)
    type(state_type), intent(inout)     :: state
  !! BC
    type(scalar_field), pointer         :: TT,SS
    type(scalar_field), pointer         :: scalar_surface
    type(mesh_type), pointer            :: ice_mesh
    character(len=FIELD_NAME_LEN)       :: bc_type
    type(scalar_field), pointer         :: T_bc,S_bc
    integer, dimension(:), allocatable  :: surface_node_list
    integer                             :: i, the_node
!! Insert BC for temperature and salinity. This could be a separate subroutine

    
        TT=> extract_scalar_field(state,"Temperature")
        SS => extract_scalar_field(state,"Salinity")
        
        !This change with bc type
        call get_option("/ocean_forcing/iceshelf_meltrate/Holland08/calculate_boundaries", bc_type)   
       
        select case(bc_type)
        case("neumann")
            T_bc => extract_scalar_field(state,"Heat_flux")
            S_bc => extract_scalar_field(state,"Salt_flux")
            do i=1,node_count(T_bc)
                call set(T_bc,node_val(T_bc,i)*10.0**3)
                call set(S_bc,node_val(S_bc,i)*10.0**3)
            enddo
            
        case("dirichlet") 
            T_bc => extract_scalar_field(state,"Tb")
            S_bc => extract_scalar_field(state,"Sb")
        case default
            FLAbort('Unknown BC for TKE')
        end select 
    
    
        ! Surface node list
        allocate(surface_node_list(key_count(sf_nodes)))
        ! Make it to vector from integer_set.
        ! sf_nodes is calculated in "melt_allocate_surface"
        surface_node_list=set2vector(sf_nodes)
        
        ! create a surface mesh to place values onto. This is for the top surface
        call get_boundary_condition(TT, 'temperature_iceshelf_BC', surface_mesh=ice_mesh) 
        call allocate(ice_surfaceT, ice_mesh, name="ice_surfaceT")
        call get_boundary_condition(SS, 'salinity_iceshelf_BC', surface_mesh=ice_mesh) 
        call allocate(ice_surfaceS, ice_mesh, name="ice_surfaceS")
        ! Define ice_surfaceT according to Heat_flux?
        ewrite(1,*) "node_count(ice_surfaceT)",node_count(ice_surfaceT)
        do i=1,node_count(ice_surfaceT)
            the_node = surface_node_list(i)
            call set(ice_surfaceT,i,node_val(T_bc,the_node))
            call set(ice_surfaceS,i,node_val(S_bc,the_node))
        enddo
        !! Temperature
        scalar_surface => extract_surface_field(TT, 'temperature_iceshelf_BC', "value")    
        call remap_field(ice_surfaceT, scalar_surface)
        !! Salinity
        scalar_surface => extract_surface_field(SS, 'salinity_iceshelf_BC', "value")    
        call remap_field(ice_surfaceS, scalar_surface)
!        ewrite(1,*) "iceBC, ice_element_list",ice_element_list 
!        ewrite(1,*) "iceBC, node_count(scalar_surface)", node_count(scalar_surface)
!        ewrite(1,*) "node_val(scalar_surface,1): ", node_val(scalar_surface,1)

!! Insert BC for salinity
      

end subroutine melt_bc



!!------------------------------------------------------------------!
!!                       Private subroutines                        !
!!------------------------------------------------------------------!

subroutine melt_surf_mesh(mesh,surface_ids,surface_mesh,surface_nodes,surface_element_list)
    
    type(mesh_type), intent(in)                       :: mesh
    type(integer_set), intent(in)                     :: surface_ids
    type(mesh_type), intent(out)                      :: surface_mesh
    integer, dimension(:), pointer, intent(out)       :: surface_nodes ! allocated and returned by create_surface_mesh
    integer, dimension(:), allocatable, intent(out)   :: surface_element_list
    !    integer, dimension(:), allocatable                :: surface_nodes_out
    type(scalar_field)                                :: surface_field   
    type(integer_set)                                 :: surface_elements
    integer :: i
    ! create a set of surface elements that have surface id in 'surface_ids'

    call allocate(surface_elements)
    do i=1, surface_element_count(mesh)
!        ewrite(1,*) "surf_normal surface_element_id(mesh, i)", surface_element_id(mesh, i)
!        ewrite(1,*) "surf_normal surface_ids", set2vector(surface_ids)
        if (has_value(surface_ids, surface_element_id(mesh, i))) then
             call insert(surface_elements, i)
             
        end if
    end do
    
    allocate(surface_element_list(key_count(surface_elements)))
    surface_element_list=set2vector(surface_elements)
    call create_surface_mesh(surface_mesh, surface_nodes, mesh, surface_element_list, name=trim(mesh%name)//"ToshisMesh")
    

end subroutine melt_surf_mesh


real function setnan(arg)
     real :: arg
        setnan = sqrt(arg)
end function

real function calc_area2(positions,nodes)
     type(vector_field), pointer        :: positions
     integer, dimension(2)              :: nodes
     real, dimension(2)                 :: coord1,coord2
!     real                               :: area
     coord1 = node_val(positions,nodes(1))
     coord2 = node_val(positions,nodes(2))
    
     calc_area2 = (sum((coord2 - coord1)**2))**0.5
     
end function

real function calc_area3(positions,nodes)
     type(vector_field), pointer        :: positions
     integer, dimension(3)              :: nodes
     real, dimension(3)                 :: coord1,coord2,coord3
     real                               :: a,b,c,s
!     real                               :: area
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

end module iceshelf_meltrate_surf_normal
