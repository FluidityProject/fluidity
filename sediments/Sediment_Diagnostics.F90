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

  use fldebug
  use global_parameters, only:FIELD_NAME_LEN, OPTION_PATH_LEN, dt, timestep
  use futils, only: int2str
  use vector_tools
  use quadrature
  use elements
  use spud
  use sparse_tools
  use fetools
  use fields
  use state_module
  use fefields
  use boundary_conditions
  use sediment, only: get_n_sediment_fields, get_sediment_item

  implicit none
  
  private

  public calculate_sediment_flux, calculate_sediment_sinking_velocity,&
       & calculate_sediment_active_layer_d50, calculate_sediment_active_layer_sigma,&
       & calculate_sediment_active_layer_volume_fractions

contains
  
  subroutine calculate_sediment_flux(state)
    !!< Calculate the advected flux of the sediment through the surfaces of
    !!< the domain.
    !!< This is determined based upon a fixed theta of 0.5
    !!< Currently erosion bc is calculated explicitly and the value is consistent here
    type(state_type), intent(inout)                        :: state
    type(mesh_type), dimension(:), allocatable             :: surface_mesh
    type surface_nodes_array
       integer, dimension(:), pointer                      :: nodes
    end type surface_nodes_array
    type(surface_nodes_array), dimension(:), allocatable   :: surface_nodes
    type(vector_field), pointer                            :: X, old_U, new_U, gravity
    type(scalar_field), dimension(:), allocatable          :: deposited_sediment, erosion
    type(scalar_field), pointer                            :: erosion_flux, bedload_field&
         &, old_sediment_field, new_sediment_field, sink_U, diagnostic_field
    type(vector_field)                                     :: U
    type(scalar_field)                                     :: masslump, sediment_field
    integer                                                :: n_sediment_fields,&
         & i_field, i_bcs, i_node, ele, n_bcs, stat
    integer, dimension(2)                                  :: surface_id_count
    integer, dimension(:), allocatable                     :: surface_ids
    integer, dimension(:), pointer                         :: to_nodes, surface_element_list
    real, dimension(:), allocatable                        :: values
    character(len=FIELD_NAME_LEN)                          :: bc_name, bc_type
    character(len=OPTION_PATH_LEN)                         :: bc_path

    ewrite(1,*) "In calculate_sediment_bedload"

    ! obtain some required model variables
    n_sediment_fields = get_n_sediment_fields()
    if (n_sediment_fields == 0) return

    X => extract_vector_field(state, "Coordinate")
    gravity => extract_vector_field(state, "GravityDirection")

    new_U => extract_vector_field(state, "Velocity")
    old_U => extract_vector_field(state, "OldVelocity")
    call allocate(U, new_U%dim, new_U%mesh, name="CNVelocity")
    call zero(U)
    call addto(U, old_U, 0.5)
    call addto(U, new_U, 0.5)
    
    ! allocate space for erosion and deposit field arrays
    allocate(erosion(n_sediment_fields))
    allocate(deposited_sediment(n_sediment_fields))
    allocate(surface_mesh(n_sediment_fields))
    allocate(surface_nodes(n_sediment_fields))

    ! first loop obtains eroded sediment quantities from reentrainment bc's (calculated
    ! in sediment::set_sediment_reentrainment)
    erosion_fields_loop: do i_field=1, n_sediment_fields 

       ! obtain scalar fields for this sediment class
       call get_sediment_item(state, i_field, new_sediment_field)
       call get_sediment_item(state, i_field, "Bedload", bedload_field)

       ! generate a surface mesh for this field
       call create_surface_mesh(surface_mesh(i_field), surface_nodes(i_field)%nodes, &
            & mesh=bedload_field%mesh, name='SurfaceMesh')        

       ! allocate a field that will hold the quantity of sediment eroded from the bed in
       ! this timestep
       call allocate(erosion(i_field), surface_mesh(i_field), "ErosionAmount")
       call zero(erosion(i_field))
       
       ! get boundary condition path and number of boundary conditions
       bc_path = trim(new_sediment_field%option_path)//'/prognostic/boundary_conditions'
       n_bcs = option_count(bc_path)

       ! Loop over boundary conditions for field
       do i_bcs=0, n_bcs-1

          ! Get name and type of boundary condition
          call get_option(trim(bc_path)//"["//int2str(i_bcs)//"]"//"/name", bc_name)
          call get_option(trim(bc_path)//"["//int2str(i_bcs)//"]"//"/type[0]/name", bc_type)

          ! find reentrainment boundary condition (if there is one)
          if ((trim(bc_type) .eq. "sediment_reentrainment")) then

             ! get boundary condition info
             call get_boundary_condition(new_sediment_field, name=bc_name, type=bc_type, &
               surface_element_list=surface_element_list)

             ! get erosion flux
             erosion_flux => extract_surface_field(new_sediment_field, bc_name=bc_name,&
                  & name="value")
             
             ! set erosion field values
             allocate(values(erosion(i_field)%mesh%shape%loc))             
             do ele=1,ele_count(erosion_flux)
                to_nodes => ele_nodes(erosion(i_field), surface_element_list(ele))
                values = ele_val(erosion_flux,ele)
                do i_node=1,size(to_nodes)
                   call set(erosion(i_field),to_nodes(i_node),values(i_node))
                end do
             end do

             call scale(erosion(i_field),dt)
             deallocate(values)
             
          end if

       end do

    end do erosion_fields_loop

    ! second loop calculates the amount of sediment that has been deposited during the
    ! timestep 
    deposit_fields_loop: do i_field=1, n_sediment_fields

       ! obtain scalar fields for this sediment class
       call get_sediment_item(state, i_field, old_sediment_field, old = .true.)
       call get_sediment_item(state, i_field, new_sediment_field)
       call allocate(sediment_field, new_sediment_field%mesh, name="CNSedimentField")
       call zero(sediment_field)
       call addto(sediment_field, old_sediment_field, 0.5)
       call addto(sediment_field, new_sediment_field, 0.5)
       call get_sediment_item(state, i_field, "Bedload", bedload_field)
       call get_sediment_item(state, i_field, "SinkingVelocity", sink_U) 
       
       ! allocate surface field that will contain the calculated deposited sediment for
       ! this timestep
       call allocate(deposited_sediment(i_field), surface_mesh(i_field), "DepositedSediment")
       call zero(deposited_sediment(i_field))

       ! For continuous fields we need a global lumped mass. For dg we'll
       ! do the mass inversion on a per face basis inside the element loop.
       if(continuity(surface_mesh(i_field))>=0) then
          call allocate(masslump, surface_mesh(i_field), "SurfaceMassLump")
          call zero(masslump)
       end if
       
       ! obtain surface ids over which to record deposition
       surface_id_count=option_shape(trim(bedload_field%option_path)//&
            &"/prognostic/surface_ids")
       allocate(surface_ids(surface_id_count(1)))
       call get_option(trim(bedload_field%option_path)//"/prognostic/surface_ids", &
            & surface_ids)

       ! loop through elements in surface field
       elements: do ele=1,element_count(deposited_sediment(i_field))

          ! check if element is on bedload surface
          if (.not.any(surface_element_id(bedload_field, ele)&
               &==surface_ids)) then
             cycle elements
          end if

          ! assemble bedload element
          call assemble_sediment_flux_ele(ele, deposited_sediment,&
               & sediment_field, X, U, sink_U, gravity, masslump, i_field)

       end do elements

       deallocate(surface_ids)
       
       ! For continuous fields we divide by the inverse global lumped mass
       if(continuity(surface_mesh(i_field))>=0) then
          where (masslump%val/=0.0)
             masslump%val=1./masslump%val
          end where
          call scale(deposited_sediment(i_field), masslump)
          call deallocate(masslump)
       end if

       ! get erosion rate diagnostic field
       call get_sediment_item(state, i_field, "BedloadDepositRate", diagnostic_field, stat)
       if (stat == 0) then
          call zero(diagnostic_field)
          do i_node = 1, node_count(surface_mesh(i_field))
             call set(diagnostic_field, surface_nodes(i_field)%nodes(i_node), &
                  & node_val(deposited_sediment(i_field), i_node))
          end do
          call scale(diagnostic_field, 1./dt)
       end if

       call deallocate(sediment_field)

    end do deposit_fields_loop
    
    ! third loop to calculate net flux of sediment for this timestep
    net_flux_loop: do i_field=1, n_sediment_fields
       
       ! obtain scalar fields for this sediment class
       call get_sediment_item(state, i_field, "Bedload", bedload_field)
       
       if (.not. have_option(trim(bedload_field%option_path)//'/prognostic/disable_calculation')) then
          ! Add on sediment falling in and subtract sediment coming out
          do i_node = 1, node_count(surface_mesh(i_field))
             ! add deposited sediment
             call addto(bedload_field, surface_nodes(i_field)%nodes(i_node), &
                  & node_val(deposited_sediment(i_field), i_node))
             ! remove eroded sediment
             call addto(bedload_field, surface_nodes(i_field)%nodes(i_node), &
                  & -1.0 * node_val(erosion(i_field), i_node))
           end do
       end if

       ewrite_minmax(deposited_sediment(i_field)) 
       ewrite_minmax(erosion(i_field)) 
       ewrite_minmax(bedload_field) 

       call deallocate(deposited_sediment(i_field))
       call deallocate(erosion(i_field))
       call deallocate(surface_mesh(i_field))
       deallocate(surface_nodes(i_field)%nodes)

    end do net_flux_loop

    call deallocate(U)
    deallocate(deposited_sediment)
    deallocate(erosion)
    deallocate(surface_mesh)
    deallocate(surface_nodes)

  end subroutine calculate_sediment_flux

  subroutine assemble_sediment_flux_ele(ele, deposited_sediment, sediment_field,&
       & X, U, sink_U, gravity, masslump, i_field)

    integer, intent(in) :: ele, i_field
    type(scalar_field), dimension(:), intent(inout) :: deposited_sediment
    type(vector_field), intent(in) :: X, U, gravity
    type(scalar_field), intent(in) :: sink_U, sediment_field
    type(scalar_field), intent(inout) :: masslump

    integer, dimension(:), pointer :: s_ele
    real, dimension(ele_loc(deposited_sediment(i_field), ele), &
         & ele_loc(deposited_sediment(i_field), ele)) :: invmass
    real, dimension(ele_loc(deposited_sediment(i_field), ele)) :: flux
    real, dimension(ele_ngi(deposited_sediment(i_field), ele)) :: detwei,&
         & G_normal_detwei, U_sink_detwei
    real, dimension(U%dim, ele_ngi(deposited_sediment(i_field), ele)) :: normal
    type(element_type), pointer :: s_shape

    s_ele=>ele_nodes(deposited_sediment(i_field), ele)
    s_shape=>ele_shape(deposited_sediment(i_field), ele)
    
    call transform_facet_to_physical(X, ele, detwei, normal)

    if(continuity(deposited_sediment(i_field))>=0) then
       call addto(masslump, s_ele, &
            sum(shape_shape(s_shape, s_shape, detwei), 1))
    else
       ! In the DG case we will apply the inverse mass locally.
       invmass=inverse(shape_shape(s_shape, s_shape, detwei))
    end if
    
    G_normal_detwei=sum(face_val_at_quad(gravity,ele)*normal,1)*detwei
    U_sink_detwei=G_normal_detwei*face_val_at_quad(sink_U,ele)  

    flux=dt*shape_rhs(s_shape, &
         face_val_at_quad(sediment_field, ele)*U_sink_detwei)

    if(continuity(deposited_sediment(i_field))<0) then
       ! DG case.
       flux=matmul(invmass, flux)
    end if

    call addto(deposited_sediment(i_field), s_ele, flux)

  end subroutine assemble_sediment_flux_ele

  subroutine calculate_sediment_sinking_velocity(state)
    
    type(state_type), intent(inout) :: state

    type(scalar_field_pointer), dimension(:), allocatable     :: sediment_concs
    type(scalar_field), pointer                               :: unhindered_sink_u, sink_u 
    type(vector_field), pointer                               :: X
    type(scalar_field)                                        :: rhs, rhs_projection
    integer                                                   :: n_sediment_fields,&
         & i_field, i_node

    ewrite(1,*) 'In calculate sediment sinking velocities'

    n_sediment_fields = get_n_sediment_fields()
    
    allocate(sediment_concs(n_sediment_fields))

    ! allocate storage for rhs and set all to 1
    call get_sediment_item(state, 1, sediment_concs(1)%ptr)
    call allocate(rhs, sediment_concs(1)%ptr%mesh, name="Rhs")
    call set(rhs, 1.0)
       
    ! get sediment concentrations and remove from rhs
    do i_field=1, n_sediment_fields
       call get_sediment_item(state, i_field, sediment_concs(i_field)%ptr)
       call addto(rhs, sediment_concs(i_field)%ptr, scale=-1.0)
    end do
    
    ! raise rhs to power of 2.39
    do i_node = 1, node_count(rhs)
       if (node_val(rhs, i_node) > 1e-3) then
          call set(rhs, i_node, node_val(rhs, i_node)**2.39)
       else
          call set(rhs, i_node, 1e-3**2.39)
       end if
    end do 

    do i_field=1, n_sediment_fields
 
       ! check for diagnostic sinking velocity
       if (have_option(trim(sediment_concs(i_field)%ptr%option_path)// &
            &'/prognostic/scalar_field::SinkingVelocity/diagnostic')) then

          ewrite(2,*) 'Calculating diagnostic sink velocity for sediment field: ' //&
               & trim(sediment_concs(i_field)%ptr%name)
          
          ! check for presence of unhindered sinking velocity value
          if (.not. have_option(trim(sediment_concs(i_field)%ptr%option_path)// &
            &'/prognostic/scalar_field::UnhinderedSinkingVelocity')) then
             FLExit('You must specify an unhindered sinking velocity field to be able to calculate diagnostic sinking velocity field values for sediments')
          endif

          unhindered_sink_u => extract_scalar_field(state, &
               & trim(sediment_concs(i_field)%ptr%name)//'UnhinderedSinkingVelocity')
          ewrite_minmax(unhindered_sink_u)   

          sink_u => extract_scalar_field(state, &
               & trim(sediment_concs(i_field)%ptr%name)//'SinkingVelocity')

          ! calculate hindered sinking velocity
          call set(sink_u, unhindered_sink_u)
          if (rhs%mesh==sink_u%mesh) then
             call scale(sink_u, rhs)
          else
             call allocate(rhs_projection, sink_u%mesh, name="RhsProjection")
             X => extract_vector_field(state, 'Coordinate')
             call project_field(rhs, rhs_projection, X)
             call scale(sink_u, rhs_projection)
             call deallocate(rhs_projection)
          end if
          ewrite_minmax(sink_u) 
       endif

    end do

    call deallocate(rhs)
    deallocate(sediment_concs)

  end subroutine calculate_sediment_sinking_velocity

  subroutine calculate_sediment_active_layer_d50(state)

    type(state_type), intent(inout)             :: state
    type(scalar_field), pointer                 :: d50
    type(scalar_field)                          :: total_bedload
    type(scalar_field_pointer), dimension(:), allocatable :: sorted_bedload
    real, dimension(:), allocatable             :: sorted_diameter
    type(scalar_field_pointer)                  :: temp_bedload
    real                                        :: temp_diameter
    real                                        :: cumulative_bedload
    logical                                     :: sorted = .false.
    integer                                     :: i_field, n_fields, i_node, stat
    real                                        :: min_bedload = 1.0e-20

    ewrite(1,*) 'In calculate sediment_active_layer_d50'

    d50 => extract_scalar_field(state, 'SedimentBedActiveLayerD50', stat)
    if (stat /= 0) return

    n_fields = get_n_sediment_fields()

    allocate(sorted_bedload(n_fields))
    allocate(sorted_diameter(n_fields))

    do i_field = 1, n_fields
       call get_sediment_item(state, i_field, 'diameter', sorted_diameter(i_field), stat)
       if (stat /= 0) FLExit('All sediment fields must have a diameter to be able to calculate the SedimentBedActiveLayerD50')
       call get_sediment_item(state, i_field, 'Bedload', sorted_bedload(i_field)&
            &%ptr, stat)
    end do

    do while (.not. sorted)
       sorted = .true.
       do i_field = 2, n_fields
          if (sorted_diameter(i_field-1) > sorted_diameter(i_field)) then
             temp_diameter = sorted_diameter(i_field)
             sorted_diameter(i_field) = sorted_diameter(i_field-1)
             sorted_diameter(i_field-1) = temp_diameter
             temp_bedload = sorted_bedload(i_field)
             sorted_bedload(i_field) = sorted_bedload(i_field-1)
             sorted_bedload(i_field-1) = temp_bedload
             sorted = .false.
          end if
       end do
    end do

    call allocate(total_bedload, sorted_bedload(1)%ptr%mesh, 'TotalBedload')
    call zero(d50)
    call zero(total_bedload)

    do i_field = 1, n_fields
       call addto(total_bedload, sorted_bedload(i_field)%ptr)
    end do
 
    nodes: do i_node = 1, node_count(d50)

       if (node_val(total_bedload, i_node) > min_bedload) then 
          i_field = 0
          cumulative_bedload = 0.0
          do while (cumulative_bedload < 0.5*node_val(total_bedload, i_node))
             i_field = i_field + 1
             cumulative_bedload = cumulative_bedload + node_val(sorted_bedload(i_field)%ptr,&
                  & i_node)
          end do
          call set(d50, i_node, sorted_diameter(i_field))
       else 
          call set(d50, i_node, 0.0)
       end if

    end do nodes

    ewrite_minmax(d50)

    deallocate(sorted_diameter)
    deallocate(sorted_bedload)
    call deallocate(total_bedload)  

  end subroutine calculate_sediment_active_layer_d50

  subroutine calculate_sediment_active_layer_sigma(state)

    type(state_type), intent(inout)                       :: state
    type(scalar_field), pointer                           :: sigma
    type(mesh_type)                                       :: surface_mesh
    integer, dimension(:), pointer                        :: surface_node_list
    type(vector_field), pointer                           :: x
    type(scalar_field_pointer), dimension(:), allocatable :: bedload
    real, dimension(:), allocatable                       :: diameter
    type(scalar_field)                                    :: mean, &
         & masslump, sigma_surface
    integer                                               :: n_fields, i_field, i_ele,&
         & i_node, stat
    integer, dimension(2)                                 :: surface_id_count
    integer, dimension(:), allocatable                    :: surface_ids

    ewrite(1,*) 'In calculate_sediment_active_layer_sigma'
    
    sigma => extract_scalar_field(state, 'SedimentBedActiveLayerSigma', stat)
    if (stat /= 0) return
    x => extract_vector_field(state, 'Coordinate')

    n_fields = get_n_sediment_fields()
    allocate(bedload(n_fields))
    allocate(diameter(n_fields))    

    ! collect information required to calculate standard deviation
    data_collection_loop: do i_field = 1, n_fields
       call get_sediment_item(state, i_field, 'diameter', diameter(i_field), stat)
       if (stat /= 0) FLExit('All sediment fields must have a diameter to be able to calculate the SedimentBedActiveLayerSigma')
       call get_sediment_item(state, i_field, 'Bedload', bedload(i_field)%ptr)      
    end do data_collection_loop
       
    ! allocate surface field that will contain the calculated sigma values
    call create_surface_mesh(surface_mesh, surface_node_list, mesh=sigma%mesh, &
         &name='SurfaceMesh')
    call allocate(sigma_surface, surface_mesh, 'SigmaSurface')
    call zero(sigma_surface)

    ! For continuous fields we need a global lumped mass. For dg we'll
    ! do the mass inversion on a per face basis inside the element loop.
    if(continuity(sigma_surface)>=0) then
       call allocate(masslump, surface_mesh, 'SurfaceMassLump')
       call zero(masslump)
    end if

    ! obtain surface ids over which to calculate sigma
    surface_id_count=option_shape(trim(sigma%option_path)//'/diagnostic/surface_ids') 
    allocate(surface_ids(surface_id_count(1)))
    call get_option(trim(sigma%option_path)//'/diagnostic/surface_ids', surface_ids)

    ! loop through elements in surface field
    elements: do i_ele=1, element_count(sigma_surface)

       ! check if element is on prescribed surface
       if (.not.any(surface_element_id(sigma, i_ele) == surface_ids)) then
          cycle elements
       end if

       ! calculate sigma
       call calculate_sediment_active_layer_element_sigma(i_ele, sigma_surface, bedload,&
            & masslump, x, diameter, n_fields)

    end do elements

    ! For continuous fields we divide by the inverse global lumped mass
    if(continuity(surface_mesh)>=0) then
       where (masslump%val/=0.0)
          masslump%val=1./masslump%val
       end where
       call scale(sigma_surface, masslump)
       call deallocate(masslump)
    end if

    ! remap surface node values on to sigma field
    do i_node = 1, node_count(surface_mesh)
       call set(sigma, surface_node_list(i_node), node_val(sigma_surface, i_node))
    end do

    ewrite_minmax(sigma)

    deallocate(bedload)
    deallocate(diameter) 
    call deallocate(sigma_surface)
    call deallocate(surface_mesh)
    deallocate(surface_node_list)
    deallocate(surface_ids)
    
  end subroutine calculate_sediment_active_layer_sigma

  subroutine calculate_sediment_active_layer_element_sigma(i_ele, sigma_surface, bedload,&
       & masslump, x, diameter, n_fields)

    integer, intent(in)                                      :: i_ele
    type(scalar_field), intent(inout)                        :: sigma_surface
    type(scalar_field_pointer), dimension(:), intent(in)     :: bedload
    type(scalar_field), intent(inout)                        :: masslump
    type(vector_field), pointer, intent(in)                  :: x
    real, dimension(:), intent(in)                           :: diameter
    integer, intent(in)                                      :: n_fields
    integer, dimension(:), pointer                           :: ele
    type(element_type), pointer                              :: shape
    real, dimension(ele_loc(sigma_surface, i_ele), &
         & ele_loc(sigma_surface, i_ele))                    :: invmass
    real, dimension(ele_ngi(sigma_surface, i_ele))           :: detwei, total_bedload, &
         & mean_diameter, sigma_squared
    real, dimension(ele_loc(sigma_surface, i_ele))           :: sigma
    integer                                                  :: i_field, i_gi
    real                                                     :: min_bedload = 1.0e-20

    ele => ele_nodes(sigma_surface, i_ele)
    shape => ele_shape(sigma_surface, i_ele)

    call transform_facet_to_physical(x, i_ele, detwei)

    if(continuity(sigma_surface)>=0) then
       call addto(masslump, ele, &
            sum(shape_shape(shape, shape, detwei), 1))
    else
       ! In the DG case we will apply the inverse mass locally.
       invmass=inverse(shape_shape(shape, shape, detwei))
    end if

    do i_gi = 1, ele_ngi(sigma_surface, i_ele)
       total_bedload(i_gi) = 0.0
       mean_diameter(i_gi) = 0.0
       sigma_squared(i_gi) = 0.0
    end do
   
    mean_calculation_loop: do i_field = 1, n_fields
       total_bedload = total_bedload + face_val_at_quad(bedload(i_field)%ptr, i_ele)
       mean_diameter = mean_diameter + face_val_at_quad(bedload(i_field)%ptr, i_ele) &
            *diameter(i_field)
    end do mean_calculation_loop
    where ((total_bedload > min_bedload)) 
       mean_diameter = mean_diameter / total_bedload
    elsewhere
       mean_diameter = 0.0
    end where

    sigma_calculation_loop: do i_field = 1, n_fields
       sigma_squared = sigma_squared + face_val_at_quad(bedload(i_field)%ptr, i_ele) &
            *(diameter(i_field) - mean_diameter)**2.0       
    end do sigma_calculation_loop
    where ((total_bedload > min_bedload)) 
        sigma_squared = sigma_squared / total_bedload
    elsewhere
        sigma_squared = 0.0
    end where
    sigma = shape_rhs(shape, dsqrt(sigma_squared) * detwei)

    if(continuity(sigma_surface)<0) then
       ! DG case.
       sigma = matmul(invmass, sigma)
    end if

    call addto(sigma_surface, ele, sigma)
    
  end subroutine calculate_sediment_active_layer_element_sigma

  subroutine calculate_sediment_active_layer_volume_fractions(state)

    type(state_type), intent(inout)                       :: state
    type(scalar_field), pointer                           :: volume_fraction, bedload
    type(mesh_type)                                       :: surface_mesh
    integer, dimension(:), pointer                        :: surface_node_list
    type(vector_field), pointer                           :: x
    type(scalar_field)                                    :: total_bedload,&
         & volume_fraction_surface, masslump
    integer                                               :: n_fields, i_field, i_ele,&
         & i_node
    integer, dimension(2)                                 :: surface_id_count
    integer, dimension(:), allocatable                    :: surface_ids

    ewrite(1,*) 'In calculate_sediment_active_layer_volume_fractions'
    
    x => extract_vector_field(state, 'Coordinate')

    n_fields = get_n_sediment_fields()
    
    ! calculate combined bedload
    data_collection_loop: do i_field = 1, n_fields
       call get_sediment_item(state, i_field, 'Bedload', bedload) 
       if (i_field == 1) then
          call allocate(total_bedload, bedload%mesh, 'TotalBedload')
          call zero(total_bedload)
       end if
       call addto(total_bedload, bedload)
    end do data_collection_loop

    calculation_loop: do i_field = 1, n_fields

       ! get sediment bedload and volume fraction fields
       call get_sediment_item(state, i_field, 'Bedload', bedload) 
       call get_sediment_item(state, i_field, 'BedloadVolumeFraction', volume_fraction)  

       ! generate surface_mesh for calculation of volume fraction and create surface field
       call create_surface_mesh(surface_mesh, surface_node_list, & 
            & mesh=bedload%mesh, name='SurfaceMesh')
       call allocate(volume_fraction_surface, surface_mesh, 'VolumeFraction')
       call zero(volume_fraction_surface) 
       
       ! For continuous fields we need a global lumped mass. For dg we'll
       ! do the mass inversion on a per face basis inside the element loop.
       if(continuity(volume_fraction_surface)>=0) then
          call allocate(masslump, surface_mesh, 'SurfaceMassLump')
          call zero(masslump)
       end if

       ! obtain sediment bedload surface ids
       surface_id_count=option_shape(trim(bedload%option_path)//&
            &"/prognostic/surface_ids")
       allocate(surface_ids(surface_id_count(1)))
       call get_option(trim(bedload%option_path)//"/prognostic/surface_ids", &
            & surface_ids)
       
       ! loop through elements in surface field
       elements: do i_ele=1, element_count(volume_fraction_surface)

          ! check if element is on prescribed surface
          if (.not.any(surface_element_id(volume_fraction, i_ele) == surface_ids)) then
             cycle elements
          end if

          ! calculate volume_fraction
          call calculate_sediment_active_layer_element_volume_fractions(i_ele,&
               & volume_fraction_surface, bedload, total_bedload, masslump, x)

       end do elements

       ! For continuous fields we divide by the inverse global lumped mass
       if(continuity(volume_fraction_surface)>=0) then
          where (masslump%val/=0.0)
             masslump%val=1./masslump%val
          end where
          call scale(volume_fraction_surface, masslump)
          call deallocate(masslump)
       end if
       
       ! remap surface node values on to sigma field
       do i_node = 1, node_count(surface_mesh)
          call set(volume_fraction, surface_node_list(i_node),&
               & node_val(volume_fraction_surface, i_node))
       end do

       ewrite_minmax(volume_fraction)

       call deallocate(volume_fraction_surface)
       call deallocate(surface_mesh)
       deallocate(surface_node_list)
       deallocate(surface_ids)
       
    end do calculation_loop

    call deallocate(total_bedload)
    
  end subroutine calculate_sediment_active_layer_volume_fractions

  subroutine calculate_sediment_active_layer_element_volume_fractions(i_ele,&
       & volume_fraction_surface, bedload, total_bedload, masslump, x)

    integer, intent(in)                                      :: i_ele
    type(scalar_field), intent(inout)                        :: volume_fraction_surface
    type(scalar_field), intent(in)                           :: bedload, total_bedload
    type(scalar_field), intent(inout)                        :: masslump
    type(vector_field), pointer, intent(in)                  :: x
    integer, dimension(:), pointer                           :: ele
    type(element_type), pointer                              :: shape
    real, dimension(ele_loc(volume_fraction_surface, i_ele), &
         & ele_loc(volume_fraction_surface, i_ele))          :: invmass
    real, dimension(ele_ngi(volume_fraction_surface, i_ele)) :: detwei, &
         & volume_fraction_at_quad 
    real, dimension(ele_loc(volume_fraction_surface, i_ele)) :: volume_fraction
    real                                                     :: min_bedload = 1.0e-20

    ele => ele_nodes(volume_fraction_surface, i_ele)
    shape => ele_shape(volume_fraction_surface, i_ele)

    call transform_facet_to_physical(x, i_ele, detwei)

    if(continuity(volume_fraction_surface)>=0) then
       call addto(masslump, ele, &
            sum(shape_shape(shape, shape, detwei), 1))
    else
       ! In the DG case we will apply the inverse mass locally.
       invmass=inverse(shape_shape(shape, shape, detwei))
    end if

    where (face_val_at_quad(total_bedload, i_ele) > min_bedload) 
       volume_fraction_at_quad = face_val_at_quad(bedload, i_ele) / &
            & face_val_at_quad(total_bedload, i_ele) 
    elsewhere
       volume_fraction_at_quad = 0.0
    end where
    volume_fraction = shape_rhs(shape, volume_fraction_at_quad * detwei)

    if(continuity(volume_fraction_surface)<0) then
       ! DG case.
       volume_fraction = matmul(invmass, volume_fraction)
    end if

    call addto(volume_fraction_surface, ele, volume_fraction)
    
  end subroutine calculate_sediment_active_layer_element_volume_fractions

end module sediment_diagnostics
