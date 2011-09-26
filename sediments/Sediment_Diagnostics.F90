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
    !!<
    !!< Since we don't actually know what the advection scheme was, this is
    !!< only an estimate based on the value at the end of the timestep.
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: surface_mesh
    type(vector_field), pointer :: X, U, gravity
    type(scalar_field) :: masslump
    type(scalar_field_pointer), dimension(:), allocatable :: sediment_field&
         &, flux_field, sink_U
    type(scalar_field), dimension(:), allocatable :: surface_field, erosion
    type(scalar_field), pointer :: fluxUp

    integer :: sediment_classes, i, node, ele, stat
    integer, dimension(2) :: surface_id_count
    real :: dt
    character(len=FIELD_NAME_LEN) :: class_name
    integer, dimension(:), allocatable :: surface_ids
    integer :: id
    integer, dimension(:), pointer :: surface_element_list
    real, dimension(:), allocatable :: values
    integer, dimension(:), pointer :: to_nodes

    
    sediment_classes = get_nSediments()
    call get_option("/timestepping/timestep", dt)
    
    allocate(sediment_field(sediment_classes))
    allocate(erosion(sediment_classes))
    allocate(surface_field(sediment_classes))
    allocate(flux_field(sediment_classes))
    allocate(sink_U(sediment_classes))

    X=>extract_vector_field(state, "Coordinate")
    !! Slightly unsure whether this should be velocity or nonlinearvelocity
    !! at this stage of the timestep.
    U=>extract_vector_field(state, "Velocity")
    gravity=>extract_vector_field(state, "GravityDirection")

    surface_id_count=option_shape("/material_phase::"//trim(state%name)//& 
         "/sediment/scalar_field::SedimentDepositionTemplate/diagnostic/surface_id&
         &s") 
    allocate(surface_ids(surface_id_count(1)))
    call get_option("/material_phase::"//trim(state%name)//& 
         "/sediment/scalar_field::SedimentDepositionTemplate/diagnostic/surface_id&
         &s", surface_ids) 

    do i=1, sediment_classes
       class_name=get_sediment_name(i)
 
       sediment_field(i)%ptr=>&
            extract_scalar_field(state,trim(class_name))

       flux_field(i)%ptr=>&
            extract_scalar_field(state,"SedimentDeposition"//trim(class_name))

       sink_U(i)%ptr=>&
            extract_scalar_field(state,trim(class_name)//"SinkingVelocity",&
            & stat=stat)
       if (stat/=0) then
          nullify(sink_U(i)%ptr)
       end if

       surface_mesh=>flux_field(i)%ptr%mesh%faces%surface_mesh
       
       call allocate(surface_field(i), surface_mesh, "SurfaceDeposition")

       call zero(surface_field(i))

       ! extract the sediment_reentrainment BC which will contain the
       ! erosion flux out of the SurfaceDeposition and store in erosion
       call allocate(erosion(i), surface_mesh, "ErosionAmount")
       id = get_sediment_bc_id(i)
       call zero(erosion(i))
       if (id .eq. -1) then
           ! there isn't one...move along
           cycle
       end if
       allocate(values(erosion(i)%mesh%shape%loc))
       fluxUp => extract_surface_field(sediment_field(i)%ptr, id, "value")
       call get_boundary_condition(sediment_field(i)%ptr, id,surface_element_list=surface_element_list)
       do ele=1,ele_count(fluxUp)
           to_nodes => ele_nodes(erosion(i), surface_element_list(ele))
           values = ele_val(fluxUp,ele)
           do node=1,size(to_nodes)
               call set(erosion(i),to_nodes(node),values(node))
           end do
       end do
       call scale(erosion(i),dt)
       deallocate(values)


    end do

    if(continuity(surface_mesh)>=0) then
       ! For continuous fields we need a global lumped mass. For dg we'll
       ! do the mass inversion on a per face basis inside the element loop.
       call allocate(masslump, surface_mesh, "SurfaceMassLump")
       call zero(masslump)
    end if
    
    ! Note that this is a loop over the surface elements.
    do ele=1,element_count(surface_field(1))
       
       if (.not.any(surface_element_id(flux_field(1)%ptr,ele)&
            &==surface_ids)) then
          cycle
       end if

       call assemble_sediment_flux_ele(ele, surface_field, sediment_field,&
            & X, U, sink_U, gravity, masslump, dt)

    end do

    if(continuity(surface_mesh)>=0) then
       where (masslump%val/=0.0)
          masslump%val=1./masslump%val
       end where
    end if
    
    do i=1, sediment_classes
       if(continuity(surface_mesh)>=0) then
          call scale(surface_field(i), masslump)
       end if

       do node=1,node_count(surface_field(i))
          
          ! Add on sediment falling in and subtract sediment coming out
          call addto(flux_field(i)%ptr, &
               flux_field(i)%ptr%mesh%faces%surface_node_list(node), &
               node_val(surface_field(i), node) - node_val(erosion(i),node))
          
       end do

       call deallocate(surface_field(i))
       call deallocate(erosion(i))

    end do

    if(continuity(surface_mesh)>=0) then
       call deallocate(masslump)
    end if    

  end subroutine calculate_sediment_flux

  subroutine assemble_sediment_flux_ele(ele, surface_field, sediment_field,&
       & X, U, sink_U, gravity, masslump, dt)
    integer, intent(in) :: ele
    type(scalar_field), dimension(:), intent(inout) :: surface_field
    type(scalar_field_pointer), dimension(:), intent(inout) ::&
         & sediment_field
    type(vector_field), intent(in) :: X, U
    type(scalar_field_pointer), dimension(:), intent(in) ::&
         & sink_U
    type(vector_field), intent(in) :: gravity
    type(scalar_field), intent(inout) :: masslump
    real, intent(in) :: dt

    integer, dimension(:), pointer :: s_ele
    ! Note that all the surface_fields are on the same mesh.
    real, dimension(ele_loc(surface_field(1), ele), &
         ele_loc(surface_field(1), ele)) :: invmass
    real, dimension(ele_loc(surface_field(1), ele)) :: flux
    real, dimension(ele_ngi(surface_field(1), ele)) :: detwei,&
         & U_normal_detwei, G_normal_detwei, U_sink_detwei
    real, dimension(U%dim, ele_ngi(surface_field(1), ele)) :: normal
    type(element_type), pointer :: s_shape
    integer :: i

    s_ele=>ele_nodes(surface_field(1), ele)
    s_shape=>ele_shape(surface_field(1), ele)
    
    call transform_facet_to_physical(X, ele, detwei, normal)

    if(continuity(surface_field(1))>=0) then
       call addto(masslump, s_ele, &
            sum(shape_shape(s_shape, s_shape, detwei), 1))
    else
       ! In the DG case we will apply the inverse mass locally.
       invmass=inverse(shape_shape(s_shape, s_shape, detwei))
    end if
    
    U_normal_detwei=sum(face_val_at_quad(U,ele)*normal,1)*detwei
    G_normal_detwei=sum(face_val_at_quad(gravity,ele)*normal,1)*detwei

    do i=1,size(surface_field)
    
       if(associated(sink_U(i)%ptr)) then
          U_sink_detwei=G_normal_detwei&
               *face_val_at_quad(sink_U(i)%ptr,ele)          
       else
          U_sink_detwei=0.0
       end if

       flux=dt*shape_rhs(s_shape, &
            face_val_at_quad(sediment_field(i)%ptr, ele)*&
            (U_normal_detwei+U_sink_detwei))
            
       if(continuity(surface_field(1))<0) then
          ! DG case.
          flux=matmul(invmass, flux)
       end if

       call addto(surface_field(i), s_ele, flux)

    end do

  end subroutine assemble_sediment_flux_ele

  subroutine calculate_sinking_velocity(state)
    
    type(state_type), intent(inout) :: state

    type(scalar_field_pointer), dimension(:), allocatable :: sediment_concs
    type(scalar_field), pointer :: unhindered_sinkU, sinkU  
    type(scalar_field) :: rhs
    integer :: sediment_classes, i

    character(len = OPTION_PATH_LEN) :: class_name

    ewrite(1,*) 'In calculate sediment sinking velocities'

    sediment_classes = get_nSediments()
    
    allocate(sediment_concs(sediment_classes))

    class_name=get_sediment_name(1)
    sediment_concs(1)%ptr => extract_scalar_field(state,trim(class_name))

    call allocate(rhs, sediment_concs(1)%ptr%mesh, name="Rhs")
    call set(rhs, 1.0)
       
    ! get sediment concentrations and remove from rhs
    do i=1, sediment_classes
       class_name=get_sediment_name(i)
       sediment_concs(i)%ptr => extract_scalar_field(state,trim(class_name))
       call addto(rhs, sediment_concs(i)%ptr, scale=-1.0)
    end do
    
    ! raise rhs to power of 2.39
    do i = 1, node_count(rhs)
      call set(rhs, i, node_val(rhs, i)**2.39)
    end do 

    do i=1, sediment_classes
       class_name=get_sediment_name(i)
 
       ! check for diagnostic sinking velocity
       if (have_option(trim(sediment_concs(i)%ptr%option_path)// &
            &'/prognostic/scalar_field::SinkingVelocity/diagnostic')) then

          ewrite(2,*) 'Calculating diagnostic sink velocity for sediment field: '//&
               & trim(class_name)
          
          ! check for presence of unhindered sinking velocity value
          if (.not. have_option(trim(sediment_concs(i)%ptr%option_path)// &
            &'/prognostic/scalar_field::UnhinderedSinkingVelocity')) then
             FLExit("You must specify an unhindered sinking velocity field to be able &&
                && to calculate diagnostic sinking velocity field values for sediments")
          endif

          unhindered_sinkU => extract_scalar_field(state, trim(class_name)// &
               &'UnhinderedSinkingVelocity')
          ewrite_minmax(unhindered_sinkU)   

          sinkU => extract_scalar_field(state, trim(class_name)// &
               &'SinkingVelocity')
          
          ! calculate hindered sinking velocity
          call set(sinkU, unhindered_sinkU)
          call scale(sinkU, rhs)
          ewrite_minmax(sinkU)   

       endif

    end do

    call deallocate(rhs)
    deallocate(sediment_concs)

  end subroutine calculate_sinking_velocity

end module sediment_diagnostics
