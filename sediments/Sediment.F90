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

  use quadrature
  use elements
  use field_derivatives
  use fields
  use sparse_matrices_fields
  use state_module
  use spud
  use global_parameters, only:   OPTION_PATH_LEN
  use state_fields_module
  use boundary_conditions
  use FLDebug

  implicit none

  public set_sediment_reentrainment, sediment_check_options, get_n_sediment_fields, &
       & get_sediment_item

  private 

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

  subroutine get_sediment_field(state, i_field, item, stat)

    ! Returns sediment field string option
    type(state_type), intent(in)                :: state
    integer, intent(in)                         :: i_field
    type(scalar_field), pointer, intent(out)    :: item
    integer, intent(out), optional              :: stat
    character(len=FIELD_NAME_LEN)               :: name

    call get_option(trim(state%option_path)//'/sediment/scalar_field['//int2str(i_field -&
         & 1)//']/name', name) 
    item => extract_scalar_field(state, trim(name), stat)

  end subroutine get_sediment_field

  subroutine get_sediment_field_name(state, i_field, item, stat)

    ! Returns sediment field string option
    type(state_type), intent(in)                :: state
    integer, intent(in)                         :: i_field
    character(len=FIELD_NAME_LEN), intent(out)  :: item
    integer, intent(out), optional              :: stat

    call get_option(trim(state%option_path)//'/sediment/scalar_field['//int2str(i_field -&
         & 1)//']/name', item) 

  end subroutine get_sediment_field_name

  subroutine get_sediment_option_string(state, i_field, option, item, stat)

    ! Returns sediment field string option
    type(state_type), intent(in)                :: state
    integer, intent(in)                         :: i_field
    character(len=*), intent(in)                :: option
    character(len=FIELD_NAME_LEN), intent(out)  :: item
    integer, intent(out), optional              :: stat

    call get_option(trim(state%option_path)//'/sediment/scalar_field['//int2str(i_field -&
         & 1)//']/prognostic/'//option, item) 

  end subroutine get_sediment_option_string

  subroutine get_sediment_option_real(state, i_field, option, item, stat)

    ! Returns sediment field real option
    type(state_type), intent(in)                :: state
    integer, intent(in)                         :: i_field
    character(len=*), intent(in)                :: option
    real, intent(out)                           :: item
    integer, intent(out), optional              :: stat
    
    call get_option(trim(state%option_path)//'/sediment/scalar_field['//int2str(i_field -&
         & 1)//']/prognostic/'//option, item, stat = stat) 

  end subroutine get_sediment_option_real

  subroutine get_sediment_option_scalar_field(state, i_field, option, item, stat)

    ! Returns sediment field related scalar field
    type(state_type), intent(in)                :: state
    integer, intent(in)                         :: i_field
    type(scalar_field), pointer, intent(out)    :: item
    character(len=*), intent(in)                :: option
    integer, intent(out), optional              :: stat
    
    character(len=FIELD_NAME_LEN)               :: field_name

    call get_option(trim(state%option_path)//'/sediment/scalar_field['//int2str(i_field -&
         & 1)//']/name', field_name) 
    item => extract_scalar_field(state, trim(field_name)//option, stat)

  end subroutine get_sediment_option_scalar_field

  subroutine set_sediment_reentrainment(state)

    type(state_type), intent(in)     :: state
    type(scalar_field), pointer      :: sediment_field
    integer                          :: i_field, i_bc, n_bc
    character(len = FIELD_NAME_LEN)  :: bc_name, bc_type, algorithm
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
    type(scalar_field), pointer               :: reentrainment, bedload, sink_U, d50
    type(scalar_field)                        :: masslump, bedload_remap
    type(tensor_field), pointer               :: viscosity_pointer
    type(tensor_field), target                :: viscosity
    type(vector_field), pointer               :: x, shear_stress
    type(mesh_type), pointer                  :: surface_mesh
    character(len = FIELD_NAME_LEN)           :: bc_name, bc_path, algorithm
    integer                                   :: stat, i_ele, i_field, i_node, i, j
    integer, dimension(:), pointer            :: surface_element_list
    real, dimension(2,2)                      :: algorithm_viscosity

    ! get boundary condition field and zero
    reentrainment => extract_surface_field(sediment_field, bc_name, 'value')
    call set(reentrainment, 0.0)

    ! get boundary condition info
    call get_boundary_condition(sediment_field, name=bc_name,&
         & surface_mesh=surface_mesh, surface_element_list=surface_element_list)

    ! get bedload field
    call get_sediment_item(state, i_field, 'SedimentBedload', bedload)

    ! get d50
    d50 => extract_scalar_field(state, 'SedimentBedActiveLayerD50', stat)
    if (stat /= 0) then
       FLExit("A SedimentBedActiveLayerD50 field must be activated to calculate reentrainm&
            &ent")
    end if

    ! get viscosity
    call get_option(trim(bc_path)//"/type[0]/viscosity", algorithm_viscosity(1,1), stat&
         &=stat)
    if (stat == 0) then
       do j = 1, 2
          do i = 1, 2
             algorithm_viscosity(i,j) = algorithm_viscosity(1,1)
          end do
       end do
       call allocate(viscosity, sediment_field%mesh, "Viscosity")
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
       
    ! loop through elements in surface field
    elements: do i_ele = 1, element_count(reentrainment)

       select case(trim(algorithm))
       case("Hill_2010")
          ! call hill_2010_reentrainment(state, sediment_field, bc_name)
       case("Garcia_1991")
          call assemble_garcia_1991_reentrainment_ele(state, i_field, i_ele, reentrainment,&
               & x, masslump, surface_mesh, surface_element_list, viscosity_pointer,&
               & shear_stress, bedload, d50)
       case default
          FLExit("A valid reentrainment algorithm must be selected")
       end select   

    end do elements

    ! invert global lumped mass for continuous fields
    if(continuity(surface_mesh)>=0) then
       where (masslump%val/=0.0)
          masslump%val=1./masslump%val
       end where
       call scale(reentrainment, masslump)
    end if

    ! check bound of entrainment so that it does not exceed the available sediment in the
    ! bed and is larger than zero.
    call allocate(bedload_remap, surface_mesh, name="bedload_remap")
    call remap_field_to_surface(bedload, bedload_remap, surface_element_list)
    nodes: do i_node = 1, node_count(reentrainment)
       call set(reentrainment, i_node, min(max(node_val(reentrainment, i_node), 0.0),&
            & node_val(bedload_remap, i_node)))
    end do nodes

    ewrite_minmax(bedload)  
    ewrite_minmax(reentrainment)  

    call deallocate(viscosity)

  end subroutine set_reentrainment_bc
  
  subroutine assemble_garcia_1991_reentrainment_ele(state, i_field, i_ele, reentrainment,&
       & x, masslump, surface_mesh, surface_element_list, viscosity, shear_stress, sink_U&
       &, d50)

    type(state_type), intent(in)                     :: state
    integer, intent(in)                              :: i_ele, i_field
    type(tensor_field), pointer, intent(in)          :: viscosity
    type(vector_field), intent(in)                   :: x, shear_stress
    type(scalar_field), intent(inout)                :: masslump
    type(scalar_field), pointer, intent(inout)       :: reentrainment
    type(scalar_field), pointer, intent(in)          :: sink_U, d50
    type(mesh_type), pointer, intent(in)             :: surface_mesh
    type(element_type), pointer                      :: shape
    integer, dimension(:), pointer, intent(in)       :: surface_element_list
    integer, dimension(:), pointer                   :: ele
    integer                                          :: i_node
    real, dimension(ele_ngi(reentrainment, i_ele))   :: detwei
    real, dimension(ele_loc(reentrainment, i_ele), &
         & ele_loc(reentrainment, i_ele))            :: invmass
    real                                             :: A, R, d, g, rho_0
    real                                             :: algorithm_viscosity
    real, dimension(ele_ngi(reentrainment, i_ele))   :: R_p, u_star, Z
    real, dimension(ele_loc(reentrainment, i_ele))   :: E

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
    R_p = sqrt(R*g*d**3)/face_val_at_quad(viscosity, surface_element_list(i_ele), 1, 1)
    
    ! calculate u_star (shear velocity)
    call get_option('/material_phase::'//trim(state%name)//'/equation_of_state/fluids/line&
         &ar/reference_density', rho_0)
    u_star = sqrt(norm2(face_val_at_quad(shear_stress, surface_element_list(i_ele)))&
         &/rho_0)

    ! calculate Z
    Z = (1 - 0.288) * u_star / face_val_at_quad(sink_U, surface_element_list(i_ele)) *&
         & R_p**0.6 * (d / face_val_at_quad(d50, surface_element_list(i_ele)))**0.2   

    ! calculate reentrainment
    E = shape_rhs(shape, A*Z**5 / (1 + A*Z**5/0.3) * detwei)

    if(continuity(reentrainment)<0) then
       ! DG case.
       E = matmul(invmass, E)
    end if

    call addto(reentrainment, ele, E)
    
  end subroutine assemble_garcia_1991_reentrainment_ele














  subroutine sediment_check_options

    character(len=FIELD_NAME_LEN)           :: field_mesh, sediment_mesh, bc_type
    character(len=OPTION_PATH_LEN)          :: field_option_path, bc_path, bc_path_i
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

          ! check reentrainment is set on a bedload surface and that a BedShearStress field is
          !  present if there is reentrainment

          ! get boundary condition path and number of boundary conditions
          nbcs=option_count(trim(field_option_path)//'/boundary_conditions')
          ! Loop over boundary conditions for field
          boundary_conditions: do i_bc=0, nbcs-1

             ! get name and type of boundary condition
             call get_option(trim(field_option_path)//'/boundary_conditions['//int2str(i_bc)//&
                  &']/type[0]/name', bc_type)

             ! check whether this is a reentrainment boundary
             if (.not. (trim(bc_type) .eq. "sediment_reentrainment")) then
                cycle boundary_conditions
             end if

             ! check a 'BedShearStress' field exists
             if (.not.(have_option('/material_phase[0]/vector_field::BedShearStress'))) then
                FLExit("Reentrainment boundary condition requires a BedShearStress field")
             end if

             ! get bedload surface ids
             bedload_surface_id_count=option_shape(trim(field_option_path)//'/scalar_field::SedimentB&
                  &edload/diagnostic/surface_ids')
             allocate(bedload_surface_ids(bedload_surface_id_count(1)))
             call get_option(trim(field_option_path)//'/scalar_field::SedimentBedload/diagnostic/surf&
                  &ace_ids', bedload_surface_ids) 

             ! get reentrainment surface ids
             bc_surface_id_count=option_shape(trim(field_option_path)//'/boundary_conditions['&
                  &//int2str(i_bc)//']/surface_ids')
             allocate(bc_surface_ids(bc_surface_id_count(1)))
             call get_option(trim(field_option_path)//'/boundary_conditions['//int2str(i_bc)//']/sur&
                  &face_ids', bc_surface_ids) 

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

!   subroutine set_sediment_reentrainment(state)

!     type(rtrmnt_info)                :: info
!     type(state_type), intent(in)     :: state
!     type(scalar_field), pointer      :: field, reentrainment
!     integer                          :: i_field, stat, i_bc, n_bc
!     character(len = FIELD_NAME_LEN)  :: bc_name, bc_type, algorithm
!     character(len = OPTION_PATH_LEN) :: bc_path, bc_path_

!     ewrite(1,*) "In set_sediment_reentrainment"

    

!     sediment_fields: do i_field = 1, info%n_fields

!        ! get boundary condition path and number of boundary conditions
!        bc_path = trim(info%fields(i_field)%base%option_path)//'/prognostic/boundary_conditions'
!        n_bc = option_count(trim(bc_path))

!        ! Loop over boundary conditions for field
!        boundary_conditions: do i_bc = 0, n_bc - 1
          
!           bc_path_ = trim(bc_path)//"["//int2str(i_bc)//"]"
!           ! Get name and type of boundary condition
!           call get_option(trim(bc_path_)//"/name", bc_name)
!           call get_option(trim(bc_path_)//"/type[0]/name", bc_type)

!           ! If this is a reentrainment boundary then store info
!           if (.not. (trim(bc_type) .eq. "sediment_reentrainment")) then
!              cycle boundary_conditions
!           end if

!           ewrite(1,*) "Setting reentrainment boundary conditions for field: ",trim(info&
!                &%fields(i_field)%base%name)

!           ! get boundary condition field and zero
!           reentrainment => extract_surface_field(info%fields(i_field)%base, bc_name, "value")
!           call set(reentrainment, 0.0)

!           call remap_info_to_boundary(info, i_field, bc_name)

!           call get_option(trim(bc_path_)//"/type[0]/algorithm", algorithm)     
!           select case(algorithm)
!           case("Hill_2010")
!              call hill_2010_reentrainment(info, i_field, reentrainment) 
!           case("Garcia_1991")
!              call garcia_1991_reentrainment(info, i_field, reentrainment)
!           case default
!              FLExit("A valid reentrainment algorithm must be selected")
!           end select          

!        end do boundary_conditions

!     end do sediment_fields  

!     call deallocate_info(info)

!   end subroutine set_sediment_reentrainment

!   subroutine garcia_1991_reentrainment(info, i_field, reentrainment)
    
!     type(rtrmnt_info), intent(in)         :: info
!     integer, intent(in)                   :: i_field
!     type(scalar_field), intent(inout)     :: reentrainment
!     logical                               :: have_d = .true.
!     real                                  :: R_p, u_star, mean_d, sigma_d, Z, d50, E, A
!     integer                               :: i_field_internal, i_node

!     do i_field_internal = 1, info%n_fields
!        if (.not. info%fields(i_field_internal)%have_d) have_d = .false.
!     end do

!     if (.not. (have_d .and. info%have_viscosity)) then
!        FLExit("All sediment fields must have a diameter, and a viscosity field must be&
!             & specified to calculate reentrainment using the Garcia_1991 formula")
!     end if

!     A = 1.3*10.0**(-7.0)

!     nodes: do i_node = 1, node_count(reentrainment)
       
!        ! calculate particle Reynolds number
!        R_p = sqrt(info%fields(i_field)%R*info%g*(info%fields(i_field)%d)**3)&
!             &/node_val(info%viscosity_iso_remap, i_node)

!        ! calculate u_star (shear velocity)
!        u_star = sqrt(norm2(node_val(info%shear_stress, i_node))/info%rho_0)

!        ! calculate d_50 (median grain size by volume)
!        call get_d50_sigma_d(info, i_node, d50, sigma_d)       

!        ! calculate Z
!        Z = (1 - 0.288 * sigma_d) * u_star / node_val(info%fields(i_field)%sink_U_remap,&
!             & i_node) * R_p**0.6 * (info%fields(i_field)%d / d50)**0.2   

!        ! calculate reentrainment (minimum is 0.0, maximum is bedload/dt)
!        E = max(min(&
!             & A*Z**5 / (1 + A*Z**5/0.3) &
!             &,0.0),node_val(info%fields(i_field)%bedload_remap, i_node)/info%dt)

!        call addto(reentrainment, i_node, E)

!     end do nodes

!   end subroutine garcia_1991_reentrainment

!   subroutine get_d50_sigma_d(info, i_node, d50, sigma_d)

!     type(rtrmnt_info), intent(in)  :: info
!     integer, intent(in)            :: i_node
!     real, intent(out)              :: d50, sigma_d
!     real, dimension(info%n_fields) :: sorted_bedload, sorted_diameter
!     logical                        :: sorted
!     real                           :: total_bedload, cumulative_bedload, temp, mean_d
!     integer                        :: i_field

!     do i_field = 1, info%n_fields
!        sorted_bedload(i_field) = node_val(info%fields(i_field)%bedload_remap, i_node)
!        sorted_diameter(i_field) = info%fields(i_field)%d
!     end do
!     sorted = .false.

!     do while (.not. sorted)
!        sorted = .true.
!        do i_field = 2, info%n_fields
!           if (sorted_diameter(i_field-1) > sorted_diameter(i_field)) then
!              temp = sorted_diameter(i_field)
!              sorted_diameter(i_field) = sorted_diameter(i_field-1)
!              sorted_diameter(i_field-1) = temp
!              temp = sorted_bedload(i_field)
!              sorted_bedload(i_field) = sorted_bedload(i_field-1)
!              sorted_bedload(i_field-1) = temp
!              sorted = .false.
!           end if
!        end do
!     end do

!     total_bedload = sum(sorted_bedload)

!     cumulative_bedload = 0.
!     i_field = 0
!     do while (cumulative_bedload < 0.5*total_bedload)
!        i_field = i_field + 1
!        cumulative_bedload = cumulative_bedload + sorted_bedload(i_field)
!     end do

!     d50 = sorted_diameter(i_field)   

!     ! calculate sigma (standard deviation of bed sediment by volume) 
!     mean_d = sum(sorted_diameter*sorted_bedload)/sum(sorted_bedload)
!     sigma_d = sqrt(sum((sorted_diameter-mean_d)**2))         

!   end subroutine get_d50_sigma_d

!   subroutine Hill_2010_reentrainment(info, i_field, reentrainment)
    
!     type(rtrmnt_info), intent(inout)           :: info
!     integer, intent(in)                        :: i_field
!     type(scalar_field), intent(inout)          :: reentrainment
!     real                                       :: E
!     integer                                    :: i_node

!     ! get or calculate critical shear stress
!     if ((.not. info%fields(i_field)%have_shear_crit) .and. info%fields(i_field)%have_d)&
!          & then
!        ! estimate of critical shear stress assuming grains larger than
!        ! 10 microns and constant viscosity
!        ! critical stress is either given by user (optional) or calculated
!        ! using Shield's formula (depends on grain size and density and
!        ! (vertical) viscosity) 
!        info%fields(i_field)%shear_crit = 0.041 * info%fields(i_field)%R * 1024. * info%g &
!             &* info%fields(i_field)%d
!     else
!        FLExit("You need to either specify a critical shear stress or a &&
!             && sediment diameter to use the j_hill formula for reentrainment")
!     end if

!     ! calculate eroded sediment flux and set reentrainment BC
!     ! we only need to add to the source the reentrainment of sediment from the
!     ! bedload into the neumann BC term
!     !
!     ! The source depends on the erodability of that sediment
!     ! (usually 1 unless you want to do something like mud which sticks),
!     ! the porosity in the bedload and the bed shear stress. 
!     !
!     ! Each sediment class has a critical shear stress, which if exceeded
!     ! by the bed shear stress, sediment is placed into suspension
!     !
!     ! loop over nodes in bottom surface
!     nodes: do i_node = 1, node_count(reentrainment)

!        ! calculate reentrainment (minimum is 0.0, maximum is bedload/dt)
!        E = max(min(&
!             &info%fields(i_field)%erod*(1-info%fields(i_field)%poro)&
!             &*((norm2(node_val(info%shear_stress, i_node)) - info%fields(i_field)&
!             &%shear_crit) / info%fields(i_field)%shear_crit)&
!             &,0.0),node_val(info%fields(i_field)%bedload_remap, i_node)/info%dt)

!        call addto(reentrainment, i_node, E)

!     end do nodes

!   end subroutine hill_2010_reentrainment

!   subroutine get_reentrainment_info(info, state)
    
!     type(rtrmnt_info), intent(inout) :: info
!     type(state_type), intent(in)     :: state
!     integer                          :: i_field, stat

!     ! find number of fields and allocate space in info object
!     info%n_fields = get_n_sediment_fields()
!     allocate(info%fields(info%n_fields))

!     ! get global information required for reentrainment algorithms
!     info%shear_stress => extract_vector_field(state, "BedShearStress")
!     info%viscosity => extract_tensor_field(state, "Viscosity", stat)
!     info%have_viscosity = .true.
!     if (stat /= 0) then
!        info%have_viscosity = .false.
!     end if
!     call get_option("/timestepping/timestep", info%dt)
!     call get_option("/physical_parameters/gravity/magnitude", info%g)
!     call get_option('/material_phase::'//trim(state%name)//'/equation_of_state/fluids/line&
!          &ar/reference_density', info%rho_0)   

!     ! get sediment field information required for reentrainment algorithms
!     sediment_fields: do i_field = 1, info%n_fields
!        call get_field_reentrainment_info(info%fields(i_field), state, i_field)
!     end do sediment_fields

!   end subroutine get_reentrainment_info

!   subroutine get_field_reentrainment_info(field_info, state, i_field)
    
!     type(rtrmnt_class_info), intent(inout) :: field_info
!     type(state_type), intent(in)           :: state
!     character(len=FIELD_NAME_LEN)          :: field_name
!     integer                                :: i_field, stat
    
!     ! field name
!     field_name = get_sediment_field_name(i_field)
!     field_info%base => extract_scalar_field(state, trim(field_name))
!     ! bedload
!     field_info%bedload => extract_scalar_field(state,trim(field_name)//"SedimentBedload",&
!          & stat=stat)
!     field_info%have_bedload = .true.
!     if (stat /= 0) then 
!        field_info%have_bedload = .false.
!     end if
!     ! diameter
!     if (have_option(trim(field_info%base%option_path)//"/prognostic/diameter")) then
!        call get_option(trim(field_info%base%option_path)//"/prognostic/diameter"&
!             &, field_info%d)
!        field_info%have_d = .true.
!     else
!        field_info%have_d = .false.
!     end if
!     ! erodability
!     call get_option(trim(field_info%base%option_path)//"/prognostic/erodability",&
!          & field_info%erod, default=1.0)
!     ! porosity
!     call get_option(trim(field_info%base%option_path)//"/prognostic/porosity", field_info&
!          &%poro, default=0.3)  
!     ! specific gravity
!     call get_option(trim(field_info%base%option_path)//"/prognostic/submerged_specific_gra&
!          &vity", field_info%R, default=1000.0)
!     ! sinking velocity
!     field_info%sink_U => extract_scalar_field(state, trim(field_name)//"SinkingVelocity")
!     ! critical shear stress
!     if (have_option(trim(field_info%base%option_path)//"/prognostic&
!          &/critical_shear_stress")) then
!        call get_option(trim(field_info%base%option_path)//"/prognostic/critical_shear_stre&
!             &ss", field_info%shear_crit)
!        field_info%have_shear_crit = .true.
!     else
!        field_info%have_shear_crit = .false.
!     end if

!   end subroutine get_field_reentrainment_info

!   subroutine remap_info_to_boundary(info, i_field, bc_name)

!     type(rtrmnt_info), intent(inout)            :: info
!     integer, intent(in)                         :: i_field
!     character(len = FIELD_NAME_LEN), intent(in) :: bc_name
!     type(mesh_type), pointer                    :: bottom_mesh
!     type(scalar_field)                          :: viscosity_iso
!     integer, dimension(:), pointer              :: surface_element_list
!     integer                                     :: i_field_internal

!     ! get boundary condition info
!     call get_boundary_condition(info%fields(i_field)%base, name=bc_name,&
!          & surface_mesh=bottom_mesh, surface_element_list=surface_element_list)

!     call deallocate_remapped_info(info)

!     ! convert viscosity to scalar field (ASSUMED ISOTROPIC) and remap to boundary
!     ! condition surface
!     call allocate(viscosity_iso, info%viscosity%mesh, name="viscosity_iso")
!     call set(viscosity_iso, info%viscosity, 1, 1)
!     call allocate(info%viscosity_iso_remap, bottom_mesh, name="viscosity_iso_remap")
!     call remap_field_to_surface(viscosity_iso, info%viscosity_iso_remap,&
!          & surface_element_list)
!     call deallocate(viscosity_iso)

!     ! remap shear stress to boundary condition surface
!     call allocate(info%shear_stress_remap, info%shear_stress%dim, bottom_mesh, name="shear&
!          &_stress_remap")
!     call remap_field_to_surface(info%shear_stress, info%shear_stress_remap, &
!          surface_element_list)

!     sediment_fields: do i_field_internal = 1, info%n_fields
       
!        ! remap bedload to boundary condition surface
!        call allocate(info%fields(i_field_internal)%bedload_remap, bottom_mesh, name&
!             &=trim(info%fields(i_field_internal)%base%name)//"bedload_remap")
!        call remap_field_to_surface(info%fields(i_field_internal)%bedload, info&
!             &%fields(i_field_internal)%bedload_remap, surface_element_list)
       
!        ! remap sink_U to boundary condition surface
!        call allocate(info%fields(i_field_internal)%sink_U_remap, bottom_mesh, name&
!             &=trim(info%fields(i_field_internal)%base%name)//"sink_U_remap")
!        call remap_field_to_surface(info%fields(i_field_internal)%sink_U, info&
!             &%fields(i_field_internal)%sink_U_remap, surface_element_list)

!     end do sediment_fields
    
!     ! if this has been remapped allocated fields will need deallocating later
!     info%remapped = .true.

!   end subroutine remap_info_to_boundary

!   subroutine deallocate_remapped_info(info)

!     type(rtrmnt_info), intent(inout) :: info
!     integer                          :: i_field

!     if (info%remapped) then
!        call deallocate(info%viscosity_iso_remap)
!        call deallocate(info%shear_stress_remap)
!        sediment_fields : do i_field = 1, info%n_fields
!           call deallocate(info%fields(i_field)%bedload_remap)
!           call deallocate(info%fields(i_field)%sink_U_remap)
!        end do sediment_fields
!     end if

!   end subroutine deallocate_remapped_info

!   subroutine deallocate_info(info)

!     type(rtrmnt_info), intent(inout) :: info

!     call deallocate_remapped_info(info)
!     deallocate(info%fields)

!   end subroutine deallocate_info

end module sediment
