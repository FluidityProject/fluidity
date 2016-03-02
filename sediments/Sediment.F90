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
  
module sediment

  use global_parameters, only:   OPTION_PATH_LEN, FIELD_NAME_LEN, dt, timestep
  use fldebug
  use futils, only: int2str
  use vector_tools
  use quadrature
  use elements
  use spud
  use fetools
  use fields
  use state_module
  use boundary_conditions
  use field_derivatives
  use sparse_matrices_fields
  use state_fields_module

  implicit none

  private 
  public set_sediment_reentrainment, sediment_check_options, get_n_sediment_fields, &
       & get_sediment_item, surface_horizontal_divergence

  interface get_sediment_item
     module procedure get_sediment_field, get_sediment_field_name,&
          & get_sediment_option_string, get_sediment_option_real,&
          & get_sediment_option_scalar_field
  end interface get_sediment_item

contains

  function get_n_sediment_fields() result (n_fields)

    ! Returns the number of sediment fields
    integer                      :: n_fields

    n_fields = option_count('/material_phase[0]/sediment/scalar_field')

    if (have_option('/material_phase[0]/sediment/scalar_field::SedimentBedActiveLayer&
         &D50')) n_fields = n_fields - 1
    if (have_option('/material_phase[0]/sediment/scalar_field::SedimentBedActiveLayer&
         &Sigma')) n_fields = n_fields - 1
    if (have_option('/material_phase[0]/sediment/scalar_field::ZeroSedimentConcentrat&
         &ionViscosity')) n_fields = n_fields - 1

  end function get_n_sediment_fields

  subroutine get_sediment_field(state, i_field, item, stat, old)

    ! Returns sediment field pointer
    type(state_type), intent(in)                :: state
    integer, intent(in)                         :: i_field
    type(scalar_field), pointer, intent(out)    :: item
    integer, intent(out), optional              :: stat
    character(len=FIELD_NAME_LEN)               :: name
    logical, intent(in), optional               :: old

    ! had to remove trim(state%option_path)// as this didn't work with flredecomp
    call get_option('/material_phase[0]/sediment/scalar_field['//int2str(i_field -&
         & 1)//']/name', name) 

    if (present(old) .and. old) then
      name = "Old"//trim(name)
    end if
    item => extract_scalar_field(state, trim(name), stat)

  end subroutine get_sediment_field

  subroutine get_sediment_field_name(state, i_field, item, stat)

    ! Returns sediment field string option
    type(state_type), intent(in)                :: state
    integer, intent(in)                         :: i_field
    character(len=FIELD_NAME_LEN), intent(out)  :: item
    integer, intent(out), optional              :: stat

    ! had to remove trim(state%option_path)// as this didn't work with flredecomp
    call get_option('/material_phase[0]/sediment/scalar_field['//int2str(i_field -&
         & 1)//']/name', item, stat=stat) 

  end subroutine get_sediment_field_name

  subroutine get_sediment_option_string(state, i_field, option, item, stat)

    ! Returns sediment field string option
    type(state_type), intent(in)                :: state
    integer, intent(in)                         :: i_field
    character(len=*), intent(in)                :: option
    character(len=FIELD_NAME_LEN), intent(out)  :: item
    integer, intent(out), optional              :: stat

    ! had to remove trim(state%option_path)// as this didn't work with flredecomp
    call get_option('/material_phase[0]/sediment/scalar_field['//int2str(i_field -&
         & 1)//']/prognostic/'//option, item, stat=stat) 

  end subroutine get_sediment_option_string

  subroutine get_sediment_option_real(state, i_field, option, item, stat, default)

    ! Returns sediment field real option
    type(state_type), intent(in)                :: state
    integer, intent(in)                         :: i_field
    character(len=*), intent(in)                :: option
    real, intent(out)                           :: item
    integer, intent(out), optional              :: stat
    real, intent(in), optional                  :: default
    
    ! had to remove trim(state%option_path)// as this didn't work with flredecomp
    call get_option('/material_phase[0]/sediment/scalar_field['//int2str(i_field -&
         & 1)//']/prognostic/'//option, item, stat = stat, default = default) 

  end subroutine get_sediment_option_real

  subroutine get_sediment_option_scalar_field(state, i_field, option, item, stat)

    ! Returns sediment field related scalar field
    type(state_type), intent(in)                :: state
    integer, intent(in)                         :: i_field
    type(scalar_field), pointer, intent(out)    :: item
    character(len=*), intent(in)                :: option
    integer, intent(out), optional              :: stat
    
    character(len=FIELD_NAME_LEN)               :: field_name

    ! had to remove trim(state%option_path)// as this didn't work with flredecomp
    call get_option('/material_phase[0]/sediment/scalar_field['//int2str(i_field -&
         & 1)//']/name', field_name) 
    item => extract_scalar_field(state, trim(field_name)//option, stat)

  end subroutine get_sediment_option_scalar_field

  subroutine set_sediment_reentrainment(state)

    type(state_type), intent(in)     :: state
    type(scalar_field), pointer      :: sediment_field
    integer                          :: i_field, i_bc, n_bc
    character(len = FIELD_NAME_LEN)  :: bc_name, bc_type
    character(len = OPTION_PATH_LEN) :: bc_path, bc_path_

    ewrite(1,*) "In set_sediment_reentrainment"

    sediment_fields: do i_field = 1, get_n_sediment_fields()

       ! extract sediment field from state
       call get_sediment_item(state, i_field, sediment_field)

       ! get boundary condition path and number of boundary conditions
       bc_path = trim(sediment_field%option_path)//'/prognostic/boundary_conditions'
       n_bc = option_count(trim(bc_path))

       ! Loop over boundary conditions for field
       boundary_conditions: do i_bc = 0, n_bc - 1
          
          ! Get name and type of boundary condition
          bc_path_ = trim(bc_path)//"["//int2str(i_bc)//"]"
          call get_option(trim(bc_path_)//"/name", bc_name)
          call get_option(trim(bc_path_)//"/type[0]/name", bc_type)

          ! skip if this is not a reentrainment boundary
          if (.not. (trim(bc_type) .eq. "sediment_reentrainment")) then
             cycle boundary_conditions
          end if

          ewrite(1,*) "Setting reentrainment boundary condition "//trim(bc_name)//" for fi&
               &eld: "//trim(sediment_field%name)

          call set_reentrainment_bc(state, sediment_field, bc_name, bc_path_, i_field)    

       end do boundary_conditions

    end do sediment_fields

  end subroutine set_sediment_reentrainment

  subroutine set_reentrainment_bc(state, sediment_field, bc_name, bc_path, i_field)

    type(state_type), intent(in)              :: state
    type(scalar_field), intent(in), pointer   :: sediment_field
    type(scalar_field), pointer               :: reentrainment, bedload, sink_U, d50,&
         & sigma, volume_fraction, diagnostic_field, old_diagnostic_field
    type(scalar_field)                        :: masslump, bedload_remap
    type(tensor_field), pointer               :: viscosity_pointer
    type(tensor_field), target                :: viscosity
    type(vector_field), pointer               :: x, shear_stress
    type(mesh_type), pointer                  :: surface_mesh
    character(len = FIELD_NAME_LEN)           :: bc_name, bc_path, algorithm
    integer                                   :: stat, i_ele, i_field, i_node, i_face, i, j
    integer, dimension(:), pointer            :: surface_element_list
    real, dimension(2,2)                      :: algorithm_viscosity

    ! get boundary condition field and zero
    reentrainment => extract_surface_field(sediment_field, bc_name, 'value')
    call set(reentrainment, 0.0)

    ! get boundary condition info
    call get_boundary_condition(sediment_field, name=bc_name,&
         & surface_mesh=surface_mesh, surface_element_list=surface_element_list)

    ! get bedload field
    call get_sediment_item(state, i_field, 'Bedload', bedload)
    
    ! get volume fraction
    call get_sediment_item(state, i_field, 'BedloadVolumeFraction', volume_fraction)

    ! get sinking velocity
    call get_sediment_item(state, i_field, 'SinkingVelocity', sink_U)

    ! get d50
    d50 => extract_scalar_field(state, 'SedimentBedActiveLayerD50', stat)

    ! get sigma
    sigma => extract_scalar_field(state, 'SedimentBedActiveLayerSigma', stat)

    call allocate(viscosity, sediment_field%mesh, "Viscosity")
    ! get viscosity
    call get_option(trim(bc_path)//"/type[0]/viscosity", algorithm_viscosity(1,1), stat&
         &=stat)
    if (stat == 0) then
       do j = 1, 2
          do i = 1, 2
             algorithm_viscosity(i,j) = algorithm_viscosity(1,1)
          end do
       end do
       call zero(viscosity)
       call set(viscosity, algorithm_viscosity)
       viscosity_pointer => viscosity
    else
       viscosity_pointer => extract_tensor_field(state, "Viscosity", stat)
       if (stat /= 0) then
          FLExit("A viscosity must be specified to calculate reentrainment")
       end if
    end if

    ! get shear stress
    shear_stress => extract_vector_field(state, "BedShearStress", stat)
    if (stat /= 0) then
       FLExit("A bed shear stress must be specified to calculate reentrainment")
    end if    

    ! get coordinate field
    x => extract_vector_field(state, "Coordinate")

    if (continuity(surface_mesh)>=0) then
       ! For continuous fields we need a global lumped mass. For dg we'll
       ! do the mass inversion on a per face basis inside the element loop.
       ! Continuity must be the same for all bedload meshes
       call allocate(masslump, surface_mesh, "SurfaceMassLump")
       call zero(masslump)
    end if

    call get_option(trim(bc_path)//"/type[0]/algorithm", algorithm) 
    ! loop through elements in surface field and calculate reentrainment
    elements: do i_ele = 1, element_count(reentrainment)

       select case(trim(algorithm))
       case("Generic")
          call assemble_generic_reentrainment_ele(state, i_field, i_ele, reentrainment,&
               & shear_stress, surface_element_list, x, masslump, sink_U, volume_fraction)
       case("Garcia_1991")
          call assemble_garcia_1991_reentrainment_ele(state, i_field, i_ele, reentrainment,&
               & x, masslump, surface_element_list, viscosity_pointer,&
               & shear_stress, d50, sink_U, sigma, volume_fraction)
       case default
          FLExit("A valid re-entrainment algorithm must be selected")
       end select   

    end do elements

    ! invert global lumped mass for continuous fields
    if(continuity(surface_mesh)>=0) then
       where (masslump%val/=0.0)
          masslump%val=1./masslump%val
       end where
       call scale(reentrainment, masslump)
       call deallocate(masslump)
    end if   

    ! check bound of entrainment so that it does not exceed the available sediment in the
    ! bed and is larger than zero.
    call allocate(bedload_remap, surface_mesh, name="bedload_remap")
    call remap_field_to_surface(bedload, bedload_remap, surface_element_list)
    nodes: do i_node = 1, node_count(reentrainment)
       if(dt/=0.0) then
          call set(reentrainment, i_node, min(max(node_val(reentrainment, i_node), 0.0),&
               & node_val(bedload_remap, i_node)/dt))
       else
          call set(reentrainment, i_node, max(node_val(reentrainment, i_node), 0.0))
       end if
    end do nodes
    
    ! store erosion rate in diagnositc field
    call get_sediment_item(state, i_field, "BedloadErosionRate", diagnostic_field, stat)
    if (stat == 0 .and. dt > 1e-15) then
       call zero(diagnostic_field)
       do i_ele = 1, ele_count(reentrainment)
          i_face=surface_element_list(i_ele)
          call set(diagnostic_field, face_global_nodes(diagnostic_field, i_face), &
               & ele_val(reentrainment, i_ele))
       end do
       ! I also need to set the old field value otherwise this get overwritten with zero straight away
       ! This is a bit messy but the important thing is that at the end of the timestep we have recorded 
       ! the erosion rate, this fools Fluidity such that we get that.
       old_diagnostic_field => extract_scalar_field(state, "Old"//trim(diagnostic_field%name), stat)
       if (stat == 0) then 
          call set(old_diagnostic_field, diagnostic_field)
       end if
    end if

    ! only for mms tests
    if (have_option(trim(bc_path)//"/type[0]/set_to_zero")) then
       ! zero reentrainment
       call zero(reentrainment)
    end if

    ewrite_minmax(bedload)  
    ewrite_minmax(reentrainment)  

    call deallocate(bedload_remap)
    call deallocate(viscosity)

  end subroutine set_reentrainment_bc
  
  subroutine assemble_garcia_1991_reentrainment_ele(state, i_field, i_ele, reentrainment,&
       & x, masslump, surface_element_list, viscosity, shear_stress, d50,&
       & sink_U, sigma, volume_fraction)

    type(state_type), intent(in)                     :: state
    integer, intent(in)                              :: i_ele, i_field
    type(tensor_field), pointer, intent(in)          :: viscosity
    type(vector_field), intent(in)                   :: x, shear_stress
    type(scalar_field), intent(inout)                :: masslump
    type(scalar_field), pointer, intent(inout)       :: reentrainment
    type(scalar_field), pointer, intent(in)          :: d50, sink_U, sigma,&
         & volume_fraction 
    integer, dimension(:), pointer, intent(in)       :: surface_element_list
    type(element_type), pointer                      :: shape
    integer, dimension(:), pointer                   :: ele
    real, dimension(ele_ngi(reentrainment, i_ele))   :: detwei
    real, dimension(ele_loc(reentrainment, i_ele), &
         & ele_loc(reentrainment, i_ele))            :: invmass
    real                                             :: A, R, d, g, density
    real, dimension(ele_ngi(reentrainment, i_ele))   :: R_p, u_star, Z
    real, dimension(ele_loc(reentrainment, i_ele))   :: E
    real, dimension(ele_ngi(reentrainment, i_ele))   :: shear, lambda_m
    real, dimension(shear_stress%dim, &
         & ele_ngi(reentrainment, i_ele))            :: shear_quad
    integer                                          :: i_gi, stat

    A = 1.3*10.0**(-7.0)
    
    ele => ele_nodes(reentrainment, i_ele)
    shape => ele_shape(reentrainment, i_ele)

    call transform_facet_to_physical(x, surface_element_list(i_ele), detwei)
    
    if(continuity(reentrainment)>=0) then
       call addto(masslump, ele, &
            sum(shape_shape(shape, shape, detwei), 1))
    else
       ! In the DG case we will apply the inverse mass locally.
       invmass=inverse(shape_shape(shape, shape, detwei))
    end if

    ! calculate particle Reynolds number
    call get_sediment_item(state, i_field, 'submerged_specific_gravity', R)
    call get_sediment_item(state, i_field, 'diameter', d)
    call get_option("/physical_parameters/gravity/magnitude", g)
    ! VISCOSITY ASSUMED TO BE ISOTROPIC - maybe should be in normal direction to surface
    where (face_val_at_quad(viscosity, surface_element_list(i_ele), 1, 1) > 0.0) 
       R_p = sqrt(R*g*d**3.0)/face_val_at_quad(viscosity, surface_element_list(i_ele), 1, 1)
    elsewhere 
       R_p = 0.0
    end where 
    
    ! calculate u_star (shear velocity)
    call get_option(trim(shear_stress%option_path)//"/diagnostic/density", density, stat)
    if (stat /= 0) density = 1.0
    ! calculate magnitude of shear stress at quadrature points
    shear_quad = face_val_at_quad(shear_stress, surface_element_list(i_ele))
    do i_gi = 1, ele_ngi(reentrainment, i_ele)
       shear(i_gi) = norm2(shear_quad(:, i_gi))
    end do
    u_star = sqrt(shear/density)

    ! calculate lambda_m
    lambda_m = 1.0 - 0.288 * face_val_at_quad(sigma, surface_element_list(i_ele))

    ! calculate Z
    where (face_val_at_quad(d50, surface_element_list(i_ele)) > 0.0)
       Z = lambda_m * u_star/face_val_at_quad(sink_U, surface_element_list(i_ele)) * R_p&
            &**0.6 * (d / face_val_at_quad(d50, surface_element_list(i_ele)))**0.2  
    elsewhere 
       Z = 0.0
    end where

    ! calculate reentrainment F*v_s*E
    E = shape_rhs(shape, face_val_at_quad(volume_fraction, surface_element_list(i_ele)) *&
         & face_val_at_quad(sink_U, surface_element_list(i_ele)) * A*Z**5 / (1 + A*Z**5&
         &/0.3) * detwei)  

    if(continuity(reentrainment)<0) then
       ! DG case.
       E = matmul(invmass, E)
    end if

    call addto(reentrainment, ele, E)
    
  end subroutine assemble_garcia_1991_reentrainment_ele

  subroutine assemble_generic_reentrainment_ele(state, i_field, i_ele, reentrainment,&
       & shear_stress, surface_element_list, x, masslump, sink_U, volume_fraction)

    type(state_type), intent(in)                     :: state
    integer, intent(in)                              :: i_ele, i_field
    type(scalar_field), pointer, intent(inout)       :: reentrainment
    type(vector_field), intent(in)                   :: x, shear_stress
    integer, dimension(:), pointer, intent(in)       :: surface_element_list
    type(scalar_field), intent(inout)                :: masslump
    type(scalar_field), pointer, intent(in)          :: sink_U, volume_fraction
    type(element_type), pointer                      :: shape
    integer, dimension(:), pointer                   :: ele
    real, dimension(ele_ngi(reentrainment, i_ele))   :: detwei
    real, dimension(ele_loc(reentrainment, i_ele), &
         & ele_loc(reentrainment, i_ele))            :: invmass
    integer, dimension(2)                            :: stat
    real                                             :: shear_crit, d, R, g, erod,&
         & poro, density, rho_0
    real, dimension(ele_loc(reentrainment, i_ele))   :: E
    real, dimension(ele_ngi(reentrainment, i_ele))   :: shear
    real, dimension(shear_stress%dim, &
         & ele_ngi(reentrainment, i_ele))            :: shear_quad

    integer                                          :: i_gi

    call get_sediment_item(state, i_field, 'critical_shear_stress', shear_crit, stat(1))
    call get_sediment_item(state, i_field, 'diameter', d, stat(2))
    ! non-dimensionalise shear stress
    call get_option('/material_phase::'//trim(state%name)//'/equation_of_state/fluids/line&
         &ar/reference_density', rho_0)
    ! get or calculate critical shear stress
    if (.not.any(stat .eq. 0)) then
       FLExit("You need to either specify a critical shear stress or a &
            &sediment diameter to use the generic formula for reentrainment")
    else if (stat(1) /= 0) then
       ! estimate of critical shear stress assuming grains larger than
       ! 10 microns and constant viscosity
       ! critical stress is either given by user (optional) or calculated
       ! using Shield's formula (depends on grain size and density and
       ! (vertical) viscosity) 
       call get_sediment_item(state, i_field, 'submerged_specific_gravity', R)
       call get_option("/physical_parameters/gravity/magnitude", g)
       shear_crit = 0.041 * R * rho_0 * g * d
    end if

    ! calculate eroded sediment flux and set reentrainment BC
    ! we only need to add to the source the reentrainment of sediment from the
    ! bedload into the neumann BC term
    !
    ! The source depends on the erodability of that sediment
    ! (usually 1 unless you want to do something like mud which sticks),
    ! the porosity in the bedload and the bed shear stress. 
    !
    ! Each sediment class has a critical shear stress, which if exceeded
    ! by the bed shear stress, sediment is placed into suspension
    !
    ! loop over nodes in bottom surface
    call get_sediment_item(state, i_field, 'erodability', erod, default=1.0)
    call get_sediment_item(state, i_field, 'bed_porosity', poro, default=0.3)
    
    ele => ele_nodes(reentrainment, i_ele)
    shape => ele_shape(reentrainment, i_ele) 

    call transform_facet_to_physical(x, surface_element_list(i_ele), detwei)  
    
    if(continuity(reentrainment)>=0) then
       call addto(masslump, ele, &
            sum(shape_shape(shape, shape, detwei), 1))
    else
       ! In the DG case we will apply the inverse mass locally.
       invmass=inverse(shape_shape(shape, shape, detwei))
    end if 

    ! calculate magnitude of shear stress at quadrature points
    shear_quad = face_val_at_quad(shear_stress, surface_element_list(i_ele))
    do i_gi = 1, ele_ngi(reentrainment, i_ele)
       shear(i_gi) = norm2(shear_quad(:, i_gi))
    end do
    ! non-dimensionalise shear stress
    call get_option(trim(shear_stress%option_path)//"/diagnostic/density", density)
    shear = shear / density

    ! calculate reentrainment F*vs*E
    E = shape_rhs(shape, face_val_at_quad(volume_fraction, surface_element_list(i_ele)) *&
         & face_val_at_quad(sink_U, surface_element_list(i_ele)) * erod * (1-poro) *&
         & (shear - shear_crit)/shear_crit * detwei) 

    if(continuity(reentrainment)<0) then
       ! DG case.
       E = matmul(invmass, E)
    end if

    call addto(reentrainment, ele, E)

  end subroutine assemble_generic_reentrainment_ele

 subroutine surface_horizontal_divergence(source, positions, output, surface_ids)
    !!< Return a field containing:
    !!<    div_HS source
    !!< where div_hs is a divergence operator restricted to the surface and
    !!< of spatial dimension one degree lower than the full mesh.
  
    type(vector_field), intent(in) :: source
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(inout) :: output
    integer, dimension(:), intent(in) :: surface_ids
    
    type(mesh_type) :: surface_mesh, source_surface_mesh, output_surface_mesh

    integer :: i, idx

    integer, dimension(surface_element_count(output)) :: surface_elements
    integer, dimension(:), pointer :: surface_nodes, source_surface_nodes,&
         output_surface_nodes
    real :: face_area, face_integral

    type(vector_field) :: surface_source, surface_positions
    type(scalar_field) :: surface_output
    real, dimension(mesh_dim(positions))  :: val
    
    call zero(output)

    !!! get the relevant surface. This wastes a bit of memory.
    idx=1
    do i = 1, surface_element_count(output)
       if (any(surface_ids == surface_element_id(positions, i))) then
          surface_elements(idx)=i
          idx=idx+1
       end if
    end do

    !!! make the surface meshes

    call create_surface_mesh(surface_mesh,  surface_nodes, &
         positions%mesh, surface_elements=surface_elements(:idx-1), &
         name='CoordinateSurfaceMesh')
    call create_surface_mesh(source_surface_mesh,  source_surface_nodes, &
         source%mesh, surface_elements=surface_elements(:idx-1), &
         name='SourceSurfaceMesh')
    call create_surface_mesh(output_surface_mesh,  output_surface_nodes, &
         output%mesh, surface_elements=surface_elements(:idx-1), &
         name='OutputSurfaceMesh')

    call allocate(surface_positions,mesh_dim(surface_mesh),surface_mesh,&
         "Coordinates")
    call allocate(surface_source,mesh_dim(surface_mesh),source_surface_mesh,&
         "Source")
    call allocate(surface_output,output_surface_mesh,"Divergence")

    do i=1,size(surface_nodes)
       val=node_val(positions,surface_nodes(i))
       call set(surface_positions,i,val(:mesh_dim(surface_mesh)))
    end do

    do i=1,size(source_surface_nodes)
       val=node_val(source,source_surface_nodes(i))
       call set(surface_source,i,val(:mesh_dim(surface_mesh)))
    end do

!!! now do the low dimensional divergence operation

    call div(surface_source,surface_positions,surface_output)

    do i=1,size(output_surface_nodes)
       call set(output,output_surface_nodes(i),node_val(surface_output,i))
    end do

!!! I *think* this surfices for parallel, due to the relatively locality of the 
!!! divergence operator

    call halo_update(output)

    call deallocate(surface_positions)
    call deallocate(surface_source)
    call deallocate(surface_output)

    call deallocate(surface_mesh)
    call deallocate(source_surface_mesh)
    call deallocate(output_surface_mesh)


    deallocate(surface_nodes, source_surface_nodes, output_surface_nodes)

  end subroutine surface_horizontal_divergence

  subroutine sediment_check_options

    character(len=FIELD_NAME_LEN)           :: field_mesh, sediment_mesh, bc_type
    character(len=OPTION_PATH_LEN)          :: field_option_path 
    integer                                 :: i_field, i_bc, i_bc_surf, i_bedload_surf,&
         & n_sediment_fields, nbcs
    integer, dimension(2)                   :: bc_surface_id_count, bedload_surface_id_count
    integer, dimension(:), allocatable      :: bc_surface_ids, bedload_surface_ids

    if (have_option('/material_phase[0]/sediment/')) then

       ewrite(1,*) 'Checking sediment model options'

       n_sediment_fields = get_n_sediment_fields()

       call get_option('/material_phase[0]/sediment/scalar_field[0]/prognostic/mesh[0]/name', sediment_mesh)

       sediment_fields: do i_field=1,n_sediment_fields

          field_option_path = '/material_phase[0]/sediment/scalar_field['//int2str(i_field - 1)//']/pro&
               &gnostic'

          ! check sinking velocity is specified for every sediment field
          if (.not.(have_option(trim(field_option_path)//'/scalar_field::SinkingVelocity')))&
               & then
             FLExit("You must specify a sinking velocity for every sediment field")
          end if

          ! check all sediment fields are on the same mesh
          call get_option(trim(field_option_path)//'/mesh[0]/name', field_mesh)
          if (.not.(trim(field_mesh) .eq. trim(sediment_mesh))) then
             FLExit("All sediment fields must be on the same mesh")
          end if

          ! check re-entrainment options
          ! get boundary condition path and number of boundary conditions
          nbcs=option_count(trim(field_option_path)//'/boundary_conditions')
          ! Loop over boundary conditions for field
          boundary_conditions: do i_bc=0, nbcs-1
          
             ! Get name and type of boundary condition
             call get_option(trim(field_option_path)//&
                  '/boundary_conditions['//int2str(i_bc)//&
                  ']/type[0]/name', bc_type)

             ! check whether this is a reentrainment boundary
             if (.not. (trim(bc_type) .eq. "sediment_reentrainment")) then
                cycle boundary_conditions
             end if

             ! check a 'BedShearStress' field exists
             if (.not.(have_option('/material_phase[0]/vector_field::BedShearStress'))) then
                FLExit("Reentrainment boundary condition requires a BedShearStress field")
             end if
             
             ! check boundary id's are the same for re-entrainment and bedload
             
             ! get bedload surface ids
             bedload_surface_id_count=option_shape(trim(field_option_path)// &
                  '/scalar_field::Bedload/prognostic/surface_ids')
             allocate(bedload_surface_ids(bedload_surface_id_count(1)))
             call get_option(trim(field_option_path)// &
                  '/scalar_field::Bedload/prognostic/surface_ids', bedload_surface_ids) 

             ! get reentrainment surface ids
             bc_surface_id_count=option_shape(trim(field_option_path)//'/boundary_conditions['&
                  &//int2str(i_bc)//']/surface_ids')
             allocate(bc_surface_ids(bc_surface_id_count(1)))
             call get_option(trim(field_option_path)// &
                  '/boundary_conditions['//int2str(i_bc)//']/surface_ids', bc_surface_ids) 

             bc_surface_id: do i_bc_surf=1, bc_surface_id_count(1)

                bedload_surface_id: do i_bedload_surf=1, bedload_surface_id_count(1)

                   if (bc_surface_ids(i_bc_surf) .eq. bedload_surface_ids(i_bedload_surf)) then
                      cycle bc_surface_id
                   end if

                end do bedload_surface_id

                FLExit("Reentrainment boundary condition is specified on a surface with no bedload")

             end do bc_surface_id

             deallocate(bc_surface_ids)
             deallocate(bedload_surface_ids)

          end do boundary_conditions

       end do sediment_fields

       ewrite(1,*) 'Sediment model options check complete'

    end if

  end subroutine sediment_check_options

end module sediment
