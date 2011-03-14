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
  use quadrature
  use elements
  use field_derivatives
  use fields
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
implicit none

  private
  real, save                      :: c0, cI, L, TI, &
                                     a, b, gammaT, gammaS,Cd, dist_meltrate
  integer, save                   :: nnodes,dimen
 
  type(vector_field), save        :: surface_positions
  type(vector_field), save        :: funky_positions 
  
  
  public :: melt_surf_init, melt_surf_calc

!  public :: keps_init, keps_cleanup, keps_tke, keps_eps, keps_eddyvisc, keps_adapt_mesh, keps_check_options

  ! Outline:
  !  - Init in populate_state.
  !  - call keps_tke (which calculates production and sets source/absorption/diffusivity for solve).
  !  - call keps_eps (which sets source/absorption/diffusivity for solve).
  !  - After keps_eps solve, keps_eddyvisc recalculates the eddy viscosity and adds it to the viscosity field.
  !  - Wall functions are added to selected boundaries in keps_bcs and wall_functions.
  !  - keps_adapt_options repopulates the fields after an adapt.
  !  - When done, clean-up.

contains

!----------
! initialise parameters based on options
!----------

 subroutine melt_surf_init(state)

    type(state_type), intent(inout)     :: state
    character(len=OPTION_PATH_LEN)      :: melt_path
    integer                             :: i
    integer                             :: surf_id
   

 

    
    !    ewrite(1,*)'Now in k_epsilon turbulence model - keps_init'
    melt_path = "/ocean_forcing/iceshelf_meltrate/Holland08"

    ! Allocate the temporary, module-level variables
    call melt_allocate_fields(state)
    ewrite(1,*) "In melt_init, done melting allocate"
    !    call melt_countBC(state)

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

    call melt_allocate_surface(state,melt_path) !Return the surface, which we calculate the meltrate
 
    ewrite(1,*) "---------End melt_surf_init---------------------------------"

 end subroutine melt_surf_init

! Calculate the melt rate
 subroutine melt_surf_calc(state)

    type(state_type), intent(inout)     :: state
    type(scalar_field), pointer         :: Tb, Sb,MeltRate
    type(scalar_field), pointer         :: scalarfield
    type(vector_field), pointer         :: velocity,positions  
    !Debugging purpose pointer
    type(scalar_field), pointer         :: T_loc,S_loc,P_loc
    type(vector_field), pointer         :: V_loc,Location 
    real, dimension(:), allocatable     :: vel
    integer                             :: i,j,k
    ! Some internal variables  
    real                             :: speed, T,S,P,Aa,Bb,Cc,topo,loc_Tb,loc_Sb,loc_meltrate,dist_top
    ! Aa*Sb^2+Bv*Sb+Cc
    real :: arg = -1.0
    !Sink mesh part
    integer                             :: ele,face,node,stat
    !! The local coordinates of the coordinate in the owning element
    ! for transform_facet_to_physical
    real, dimension(:,:), allocatable   :: normal,av_normal
    real, dimension(:), allocatable     :: local_coord,coord,xyz
    integer, dimension(:), pointer      :: ele_faces
    type(element_type), pointer         :: x_shape_f
    real                                :: al,be,ga
    integer, dimension(:), allocatable  :: node_lists !Store the list of elements size(dim+1)
    type(scalar_field)                  :: re_temperature,re_salinity,re_pressure
    type(vector_field)                  :: re_velocity
    type(scalar_field), pointer         :: sc_nodes
    ewrite(1,*) "In melt_surf_calc, Hello world!"

    MeltRate => extract_scalar_field(state,"MeltRate")
    Tb => extract_scalar_field(state,"Tb")
    Sb => extract_scalar_field(state,"Sb")
    call set(MeltRate,setnan(arg))
    call set(Tb,setnan(arg))
    call set(Sb,setnan(arg))
! Debugging
    T_loc => extract_scalar_field(state,"Tloc")
    S_loc => extract_scalar_field(state,"Sloc")
    P_loc => extract_scalar_field(state,"Ploc")
    V_loc => extract_vector_field(state,"Vloc")
    Location => extract_vector_field(state,"Location")
    call set(T_loc,setnan(arg))
    call set(S_loc,setnan(arg))
    call set(P_loc,setnan(arg))
    allocate(vel(V_loc%dim))
    vel = setnan(arg)
    call set(V_loc,vel)
    call set(Location,vel)
    

    positions => extract_vector_field(state,"Coordinate")
    scalarfield => extract_scalar_field(state,"Temperature")
    ! Remapping business
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

    ewrite(1,*) "Start the shady stuff node_count(funky_positions): ", node_count(funky_positions)
    allocate(local_coord(positions%dim+1))
    allocate(coord(positions%dim))

   
    do i=1,node_count(funky_positions)

    !!! Interpolating
    ewrite(1,*) "Funk coord: ", node_val(funky_positions,i)
        coord = node_val(funky_positions,i)
        ewrite(1,*) "Funk coord, coord: ", coord
        call picker_inquire(positions, coord, ele, local_coord)
        ewrite(1,*) "Funk coord, picked, ele,local_coord: ", ele,local_coord
        ewrite(1,*) "Funk coord,ele_and_faces_loc(positions, ele): ", ele_and_faces_loc(positions, ele)
        allocate(node_lists(ele_and_faces_loc(positions, ele))) !Number of nodes per element
        node_lists =  ele_nodes(positions,ele)
        
        coord = 0.0
        T = 0.0
        S = 0.0
        P = 0.0
        speed = 0.0
        vel = 0.0
        ewrite(1,*) "Funk coord, node_lists: ", node_lists
        ewrite(1,*) "Funk coord, local_coord: ", local_coord  

        do j=1,size(node_lists)
            node = node_lists(j)
            coord = coord + node_val(positions,node)*local_coord(j)
            T = T + node_val(re_temperature,node)*local_coord(j)
            S = S + node_val(re_salinity,node)*local_coord(j)
            P = P + node_val(re_pressure,node)*local_coord(j)
            speed = speed + sqrt(sum(node_val(re_velocity,node)*node_val(re_velocity,node)))*local_coord(j)**2
            vel = vel + node_val(re_velocity,node)*local_coord(j)
        enddo
        speed = sqrt(sum(vel**2))
        ewrite(1,*) "Funk, speed: ", speed
        ewrite(1,*) "Funk, vel: ", vel
        ewrite(1,*) "Funk node_val(funky_positions,i): ", node_val(funky_positions,i)
        deallocate(node_lists)
        topo = -7.53e-8*P ! constant = -7.53e-8 [C Pa^(-1)] comes from Holland and Jenkins Table 1

        !! Define Aa,Bb,Cc
        Aa = -gammaS*speed*cI*a + a*c0*gammaT*speed
        Bb = -gammaS*speed*L + gammaS*speed*S*cI*a
        Bb = Bb - gammaS*speed*cI*(b+topo) + gammaS*speed*cI*TI
        Bb = Bb - c0*gammaT*speed*T + c0*gammaT*speed*(b+topo)

        Cc = gammaS*speed*S*L +gammaS*speed*S*cI*(b+topo) + gammaS*speed*S*(-cI*TI)
        !! This could be a linear if Aa=0
        if (Aa .eq. 0.0) then
            loc_Sb = -Cc/Bb
        else
        !! Calculate for 2nd order polynomial and have if statement to pick >0
        loc_Sb = (-Bb + sqrt(Bb**2 - 4.0*Aa*Cc))/(2.0*Aa)
        if (loc_Sb .lt. 0.0) then
            loc_Sb = (-Bb - sqrt(Bb**2 - 4.0*Aa*Cc))/(2.0*Aa)
            endif
        endif

        loc_Tb = a*loc_Sb + b + topo
        loc_meltrate = gammaS*speed*(S-loc_Sb)/loc_Sb


        ewrite(1,*) "melt_rate: ",loc_meltrate
        ewrite(1,*) "tLHS: ", c0*gammaT*speed*(T-loc_Tb)
        ewrite(1,*) "tRHS: ", loc_meltrate*L + loc_meltrate*cI*(loc_Tb-TI)
        ewrite(1,*) "sLHS: ", gammaS*speed*(S-loc_Sb)
        ewrite(1,*) "sRHS: ", loc_meltrate*loc_Sb
        ewrite(1,*) "bLHS: ", loc_Tb
        ewrite(1,*) "bRHS: ", a*loc_Sb + b + topo

        call set(MeltRate, i, loc_meltrate)
        call set(Tb, i, loc_Tb)
        call set(Sb, i, loc_Sb)
        call set(T_loc,i,T)
        call set(S_loc,i,S)
        call set(P_loc,i,P)
        call set(V_loc,i,vel)
        call set(Location,i,node_val(funky_positions,i))

    enddo

      
    deallocate(local_coord)
    deallocate(coord)
    ewrite(1,*) "melt_surf_calc_NORMAL, ------END------"
end subroutine melt_surf_calc





!!------------------------------------------------------------------!
!!                       Private subroutines                        !
!!------------------------------------------------------------------!
subroutine melt_allocate_fields(state)

    type(state_type), intent(inout) :: state
    type(vector_field), pointer     :: vectorField
    type(scalar_field), pointer     :: scalarField
    
    integer :: i
    ewrite(1,*) "In melt_allocate_fields"

    vectorField => extract_vector_field(state, "Velocity")
    nnodes = node_count(vectorField)
    
end subroutine melt_allocate_fields

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
    ewrite(1,*) "melt_surf_mesh, surface_ids:", surface_ids

    call allocate(surface_elements)
    do i=1, surface_element_count(mesh)
        ewrite(1,*) "In melt_surf_init_sub, i: ", i
        ewrite(1,*) "In melt_surf_init_sub, surface_element_id(mesh, i): ", surface_element_id(mesh, i)
        if (has_value(surface_ids, surface_element_id(mesh, i))) then
             call insert(surface_elements, i)
             
        end if
    end do
    
    allocate(surface_element_list(key_count(surface_elements)))
    surface_element_list=set2vector(surface_elements)
    ewrite(1,*) "In melt_surf_init_sub, set2vector done"
    call create_surface_mesh(surface_mesh, surface_nodes, mesh, surface_element_list, name=trim(mesh%name)//"ToshisMesh")
    

     ewrite(1,*) "In melt_surf_init_NORMAL, ------END------"
end subroutine melt_surf_mesh

subroutine melt_allocate_surface(state,path)

    type(state_type), intent(inout)     :: state
    type(vector_field), pointer         :: positions
    type(vector_field)                  :: surf_positions, surf_p


    type(scalar_field)                  :: sc_out
    type(scalar_field), pointer         :: sc_nodes
    type(mesh_type), pointer            :: mesh
    type(mesh_type)                     :: surface_mesh
    type(integer_set)                   :: surface_elements
    integer, dimension(:), pointer      :: surface_nodes ! allocated and returned by create_surface_mesh
    integer, dimension(:), allocatable  :: surface_element_list
    integer    :: ele, face,i,j,k
!! The local coordinates of the coordinate in the owning element
    real, dimension(:), allocatable    :: local_coord, local_coord_surf
    real, dimension(:), allocatable     :: coord
! for transform_facet_to_physical
    real, dimension(:,:), allocatable   :: normal
    real, dimension(:), allocatable     :: av_normal !Average of normal
    integer, dimension(:), pointer      :: ele_faces,face_neighs
    type(element_type), pointer     :: x_shape_f
    real                                :: al,be,ga
    real, dimension(:), allocatable     :: xyz !New location of surface_mesh, dist_meltrate away from the boundary
!integer, save         :: melt_surfaceID
    character(len=OPTION_PATH_LEN)      :: path
    type(integer_set)                   :: surface_ids
    integer                             :: surf_id
    real, dimension(:,:), allocatable   :: quad_val

!allocate(surface_ids)
    ewrite(1,*) "path: ", path
    call get_option(trim(path)//'/melt_surfaceID',surf_id)
    ewrite(1,*) "surf_id: ", surf_id
    ewrite(1,*) "surface_ids: ", surface_ids
    call allocate(surface_ids)
    call insert(surface_ids,surf_id)
    ewrite(1,*) "aft_surface_ids: ", surface_ids
    ewrite(1,*) "surf_id: ", surf_id


    mesh => extract_mesh(state,"VelocityMesh")
! Input, mesh, surface_id
! Output surface_mesh,surface_element_list
    call melt_surf_mesh(mesh,surface_ids,surface_mesh,surface_nodes,surface_element_list)

    positions => extract_vector_field(state,"Coordinate")
!SINK THE surface_mesh!!
!   call get_option("/geometry/dimension/", dimension)
    call allocate(surf_positions,positions%dim,surface_mesh,"MySurfacePosition")
    call allocate(funky_positions,surf_positions%dim, surf_positions%mesh,"MyFunkyPosition")
! Remap the positions vector to the surface mesh, surface_positions
    call remap_field_to_surface(positions, surf_positions, surface_element_list)  
!! something new here
call allocate(surf_p,positions%dim-1,surface_mesh,"surf_p")
call remap_field_to_surface(positions, surf_p, surface_element_list)  
!  call remap_field_to_surface(positions, funky_positions, surface_element_list)
    ewrite(1,*) "melt_allocate_surface, mesh_dim(positions%mesh): ",mesh_dim(positions%mesh)
    ewrite(1,*) "melt_allocate_surface, mesh_dim(surf_positions%mesh): ",mesh_dim(surf_positions%mesh)
    ewrite(1,*) "melt_allocate_surface, mesh_dim(funky_positions%mesh): ",mesh_dim(funky_positions%mesh)

   
    allocate(coord(positions%dim))
   
    ewrite(1,*) "melt_allocate_surface, surface_element_list: ", surface_element_list
    ewrite(1,*) "melt_allocate_surface, count(surface_element_list: ", size(surface_element_list)
    ewrite(1,*) "melt_allocate_surface, node_count(surf_positions): ", node_count(surf_positions)
    ewrite(1,*) "melt_allocate_surface, node_count(surface_nodes): ", size(surface_nodes)
!! Need to loop over i=1,node_count(surf_positions)
    allocate(av_normal(positions%dim))
    allocate(xyz(funky_positions%dim))
   
    call allocate(sc_out,surf_positions%mesh,name="scout")
    call insert(state,sc_out,name="SurfaceNodes")
    sc_nodes => extract_scalar_field(state,"SurfaceNodes")
    ewrite(1,*) "face_count(surf_p): ", face_count(surf_p)
    ewrite(1,*) "surface_element_count(mesh): ", surface_element_count(mesh)
    ewrite(1,*) "surface_element_count(positions): ", surface_element_count(positions)
    ewrite(1,*) "surface_element_list: ", surface_element_list
    k=1

    do i=1,size(surface_element_list)
        face = surface_element_list(i)
        ewrite(1,*) "--------------------------------------------- ", i
 
!!!!!!!!!!!!!!!!!!!!!!! Capability to calculate the normal of the neibors
!                face_neighs =>face_neigh_vector(positions, i)
        ewrite(1,*) "my face: ", face
        ewrite(1,*) "face that occupies the node: ", face_global_nodes(positions, face)
        ewrite(1,*) "face_neigns: ", face_neigh_vector(positions, face)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        x_shape_f=>face_shape(positions,face)
        allocate(normal(positions%dim,x_shape_f%ngi)) !(dim x x_shape_f%ngi)
        call transform_facet_to_physical(positions, face, normal=normal)
        ewrite(1,*) "normal: ", normal
        av_normal = sum(normal,2) !May need some clever averaging
        av_normal = av_normal/(sum(av_normal*av_normal))**(0.5) 
        ewrite(1,*) "av_normal: ", av_normal
        deallocate(normal)
!Since av_normal is outward to the boundary, we will put -sign. We want the vector into the boundary
        av_normal = -av_normal
!This is a distance away from the boundary where you calculate meltrate 
        dist_meltrate = abs(dist_meltrate)
        ewrite(1,*) "face_ngi(positions,face)", face_ngi(positions,face) 
        allocate(quad_val(positions%dim, face_ngi(positions, face)))
        ewrite(1,*) "quad allocated"
        quad_val = face_val_at_quad_vector(positions, face)
        ewrite(1,*) "quad_val", quad_val
                
        do j=1,2 !Little sketch 2 comes from the fact that there are 2 nodes per face

!                coord(1) = quad_val(1,j) !neighboring points
!                coord(2) = quad_val(2,j)
                coord(1) = quad_val(1,1+(j-1)*3) !neighboring points
                coord(2) = quad_val(2,1+(j-1)*3)

                ewrite(1,*) "melt_allocate_surface, coord(1): ", coord(1)
                ewrite(1,*) "melt_allocate_surface, coord(2): ", coord(2)
!For 2D
                if (positions%dim .eq. 2) then
                        al = av_normal(1)
                        be = av_normal(2)
                        xyz(1) = al*dist_meltrate + coord(1)
                        xyz(2) = be*dist_meltrate + coord(2)
                endif

                call set(funky_positions,k,xyz)
                k=k+1
        enddo

        ewrite(1,*) "hello548 size(surf_positions): ", size(surf_positions%val)
        ewrite(1,*) "hello548 size(surface_element_list): ", size(surface_element_list)
        deallocate(quad_val)


    enddo

        ewrite(1,*) "k: ", k
    ewrite(1,*) "node_count(funky_positions): ", node_count(funky_positions)
    ewrite(1,*) "node_count(positions): ", node_count(positions)
    ewrite(1,*) "size(surface_element_list): ", size(surface_element_list)
!ewrite(1,*) "melt_allocate_surface, node_count(funky_positions): ", node_count(funky_positions)

! Knowing the funcky_positions we can allocate some output
!
!        call insert(state,sc_out,name="MeltRate")
!call insert(state,sc_out,name="Tb")
!call insert(state,sc_out,name="Sb")
    call deallocate(surface_ids)
        

end subroutine melt_allocate_surface

real function setnan(arg)
     real :: arg
        setnan = sqrt(arg)
end function




end module iceshelf_meltrate_surf_normal
