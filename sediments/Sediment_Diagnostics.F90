!    Copyright (C) 2006-2009 Imperial College London and others.
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

module sediment_diagnostics
  use global_parameters, only:FIELD_NAME_LEN, OPTION_PATH_LEN
  use state_module
  use fields
  use vector_tools
  use spud
  use fetools
  use sediment
  use boundary_conditions

  implicit none
  
  private

  public calculate_sediment_flux, calculate_sinking_velocity

contains
  
  subroutine calculate_sediment_flux(state)
    !!< Calculate the advected flux of the sediment through the surfaces of
    !!< the domain.
    !!< Since we don't actually know what the advection scheme was, this is
    !!< only an estimate based on the value at the end of the timestep.
    type(state_type), intent(inout)                        :: state
    type(mesh_type), pointer                               :: surface_mesh
    type(vector_field), pointer                            :: X, U, gravity
    type(scalar_field_pointer), dimension(:), allocatable  :: sediment_field&
         &, deposition_field, sink_U
    type(scalar_field), dimension(:), allocatable          :: surface_field, erosion
    type(scalar_field), pointer                            :: erosion_flux
    type(scalar_field)                                     :: masslump
    integer                                                :: n_sediment_fields,&
         & i_field, i_bcs, node, ele, stat, n_bcs
    integer, dimension(2)                                  :: surface_id_count
    integer, dimension(:), allocatable                     :: surface_ids
    integer, dimension(:), pointer                         :: to_nodes, surface_element_list
    real                                                   :: dt
    real, dimension(:), allocatable                        :: values
    character(len=FIELD_NAME_LEN)                          :: field_name, bc_name,&
         & bc_type
    character(len=OPTION_PATH_LEN)                         :: bc_path

    ewrite(1,*) "In calculate_sediment_flux"

    n_sediment_fields = get_n_sediment_fields()
    call get_option("/timestepping/timestep", dt)
    
    allocate(sediment_field(n_sediment_fields))
    allocate(erosion(n_sediment_fields))
    allocate(surface_field(n_sediment_fields))
    allocate(deposition_field(n_sediment_fields))
    allocate(sink_U(n_sediment_fields))

    X => extract_vector_field(state, "Coordinate")
    U => extract_vector_field(state, "Velocity")
    gravity => extract_vector_field(state, "GravityDirection")

    do i_field=1, n_sediment_fields 
       
       field_name = get_sediment_field_name(i_field)
 
       sediment_field(i_field)%ptr => &
            extract_scalar_field(state,trim(field_name))

       deposition_field(i_field)%ptr => &
            extract_scalar_field(state,trim(field_name)//"SedimentDeposition")

       sink_U(i_field)%ptr => &
            extract_scalar_field(state,trim(field_name)//"SinkingVelocity",&
            & stat=stat)

       surface_mesh => deposition_field(i_field)%ptr%mesh%faces%surface_mesh
       
       ! allocate surface field that will contain calculated deposited sediment
       call allocate(surface_field(i_field), surface_mesh, "SurfaceDeposition")
       call zero(surface_field(i_field))

       ! extract the sediment_reentrainment BC which will contain the erosion flux
       ! out of the SurfaceDeposition field and store it in erosion
       call allocate(erosion(i_field), surface_mesh, "ErosionAmount")
       call zero(erosion(i_field))
       
       ! get boundary condition path and number of boundary conditions
       n_bcs = option_count(trim(sediment_field(i_field)%ptr%option_path)//'/prognostic/boundary_co&&
            && nditions')

       ! Loop over boundary conditions for field
       do i_bcs=0, n_bcs-1

          ! Get name and type of boundary condition
          call get_option(trim(bc_path)//"["//int2str(i_bcs)//"]"//"/name", bc_name)
          call get_option(trim(bc_path)//"["//int2str(i_bcs)//"]"//"/type[0]/name", bc_type)

          ! find reentrainment boundary condition (if there is one)
          if ((trim(bc_type) .eq. "sediment_reentrainment")) then

             ! get boundary condition info
             call get_boundary_condition(sediment_field(i_field)%ptr, name=bc_name, type=bc_type, &
               surface_element_list=surface_element_list)

             ! get erosion flux
             erosion_flux => extract_surface_field(sediment_field(i_field)%ptr, & 
                  bc_name=bc_name, name="value")
             
             ! set erosion field values
             allocate(values(erosion(i_field)%mesh%shape%loc))             
             do ele=1,ele_count(erosion_flux)
                to_nodes => ele_nodes(erosion(i_field), surface_element_list(ele))
                values = ele_val(erosion_flux,ele)
                do node=1,size(to_nodes)
                   call set(erosion(i_field),to_nodes(node),values(node))
                end do
             end do
             call scale(erosion(i_field),dt)
             deallocate(values)
             
          end if

       end do

    end do

    if(continuity(surface_mesh)>=0) then
       ! For continuous fields we need a global lumped mass. For dg we'll
       ! do the mass inversion on a per face basis inside the element loop.
       ! Continuity must be the same for all deposition meshes
       call allocate(masslump, surface_mesh, "SurfaceMassLump")
       call zero(masslump)
    end if
    
    sediment_fields: do i_field=1, n_sediment_fields

       ! get deposition surface ids
       surface_id_count=option_shape(trim(deposition_field(i_field)%ptr%option_path)//"/diagnos&
            &tic/surface_ids") 
       allocate(surface_ids(surface_id_count(1)))
       call get_option(trim(deposition_field(i_field)%ptr%option_path)//"/diagnostic/surface_id&&
            && s", surface_ids) 
       
       ! loop through elements in surface field
       elements: do ele=1,element_count(surface_field(i_field))

          ! check if element is on deposition surface
          if (.not.any(surface_element_id(deposition_field(i_field)%ptr,ele)&
               &==surface_ids)) then
             cycle elements
          end if

          ! assemble deposition element
          call assemble_sediment_flux_ele(ele, surface_field, sediment_field,&
               & X, U, sink_U(i_field)%ptr, gravity, masslump, dt, i_field)

       end do elements

    end do sediment_fields

    ! invert global lumped mass for continuous fields
    if(continuity(surface_mesh)>=0) then
       where (masslump%val/=0.0)
          masslump%val=1./masslump%val
       end where
    end if
    
    ! calculate net flux of sediment
    do i_field=1, n_sediment_fields

       if(continuity(surface_mesh)>=0) then
          call scale(surface_field(i_field), masslump)
       end if

       do node=1,node_count(surface_field(i_field))
          
          ! Add on sediment falling in and subtract sediment coming out
          call addto(deposition_field(i_field)%ptr, &
               deposition_field(i_field)%ptr%mesh%faces%surface_node_list(node), &
               node_val(surface_field(i_field), node) - node_val(erosion(i_field),node))
          
       end do

       call deallocate(surface_field(i_field))
       call deallocate(erosion(i_field))

    end do

    if(continuity(surface_mesh)>=0) then
       call deallocate(masslump)
    end if    

  end subroutine calculate_sediment_flux

  subroutine assemble_sediment_flux_ele(ele, surface_field, sediment_field,&
       & X, U, sink_U, gravity, masslump, dt, i_field)

    integer, intent(in)                                          :: ele, i_field
    type(scalar_field), dimension(:), intent(inout)              :: surface_field
    type(scalar_field_pointer), dimension(:), intent(inout)      :: sediment_field
    type(vector_field), intent(in)                               :: X, U, gravity
    type(scalar_field), pointer, intent(in)                      :: sink_U
    type(scalar_field), intent(inout)                            :: masslump
    real, intent(in)                                             :: dt
    integer, dimension(:), pointer                               :: s_ele
    real, dimension(ele_loc(surface_field(i_field), ele), &
         & ele_loc(surface_field(i_field), ele))                 :: invmass
    real, dimension(ele_loc(surface_field(i_field), ele))        :: flux
    real, dimension(ele_ngi(surface_field(i_field), ele))        :: detwei,&
         & U_normal_detwei, G_normal_detwei, U_sink_detwei
    real, dimension(U%dim, ele_ngi(surface_field(i_field), ele)) :: normal
    type(element_type), pointer                                  :: s_shape

    s_ele=>ele_nodes(surface_field(i_field), ele)
    s_shape=>ele_shape(surface_field(i_field), ele)
    
    call transform_facet_to_physical(X, ele, detwei, normal)

    if(continuity(surface_field(i_field))>=0) then
       call addto(masslump, s_ele, &
            sum(shape_shape(s_shape, s_shape, detwei), 1))
    else
       ! In the DG case we will apply the inverse mass locally.
       invmass=inverse(shape_shape(s_shape, s_shape, detwei))
    end if
    
    U_normal_detwei=sum(face_val_at_quad(U,ele)*normal,1)*detwei
    G_normal_detwei=sum(face_val_at_quad(gravity,ele)*normal,1)*detwei

    if(associated(sink_U)) then
       U_sink_detwei=G_normal_detwei&
            *face_val_at_quad(sink_U,ele)          
    else
       U_sink_detwei=0.0
    end if

    flux=dt*shape_rhs(s_shape, &
         face_val_at_quad(sediment_field(i_field)%ptr, ele)*&
         (U_normal_detwei+U_sink_detwei))

    if(continuity(surface_field(i_field))<0) then
       ! DG case.

       flux=matmul(invmass, flux)
    end if

    call addto(surface_field(i_field), s_ele, flux)

  end subroutine assemble_sediment_flux_ele

  subroutine calculate_sinking_velocity(state)
    
    type(state_type), intent(inout) :: state

    type(scalar_field_pointer), dimension(:), allocatable     :: sediment_concs
    type(scalar_field), pointer                               :: unhindered_sink_u, sink_u  
    type(scalar_field)                                        :: rhs
    integer                                                   :: n_sediment_fields,&
         & i_field, i_node
    character(len = OPTION_PATH_LEN)                          :: field_name

    ewrite(1,*) 'In calculate sediment sinking velocities'

    n_sediment_fields = get_n_sediment_fields()
    
    allocate(sediment_concs(n_sediment_fields))

    ! allocate storage for rhs and set all to 1
    field_name = get_sediment_field_name(1)
    sediment_concs(1)%ptr => extract_scalar_field(state, field_name)
    call allocate(rhs, sediment_concs(1)%ptr%mesh, name="Rhs")
    call set(rhs, 1.0)
       
    ! get sediment concentrations and remove from rhs
    do i_field=1, n_sediment_fields
       field_name = get_sediment_field_name(i_field)
       sediment_concs(i_field)%ptr => extract_scalar_field(state, field_name)
       call addto(rhs, sediment_concs(i_field)%ptr, scale=-1.0)
    end do
    
    ! raise rhs to power of 2.39
    do i_node = 1, node_count(rhs)
      call set(rhs, i_node, node_val(rhs, i_node)**2.39)
    end do 

    do i_field=1, n_sediment_fields
 
       ! check for diagnostic sinking velocity
       if (have_option(trim(sediment_concs(i_field)%ptr%option_path)// &
            &'/prognostic/scalar_field::SinkingVelocity/diagnostic')) then

          ewrite(2,*) 'Calculating diagnostic sink velocity for sediment field: '//&
               & trim(sediment_concs(i_field)%ptr%name)
          
          ! check for presence of unhindered sinking velocity value
          if (.not. have_option(trim(sediment_concs(i_field)%ptr%option_path)// &
            &'/prognostic/scalar_field::UnhinderedSinkingVelocity')) then
             FLExit("You must specify an unhindered sinking velocity field to be able &&
                && to calculate diagnostic sinking velocity field values for sediments")
          endif

          unhindered_sink_u => extract_scalar_field(state, &
               & trim(sediment_concs(i_field)%ptr%name)//'UnhinderedSinkingVelocity')
          ewrite_minmax(unhindered_sink_u)   

          sink_u => extract_scalar_field(state, &
               & trim(sediment_concs(i_field)%ptr%name)//'SinkingVelocity')
          
          ! calculate hindered sinking velocity
          call set(sink_u, unhindered_sink_u)
          call scale(sink_u, rhs)
          ewrite_minmax(sink_u)   

       endif

    end do

    call deallocate(rhs)
    deallocate(sediment_concs)

  end subroutine calculate_sinking_velocity

end module sediment_diagnostics
