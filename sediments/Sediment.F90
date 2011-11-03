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
       & get_sediment_field_name

  private

contains

  function get_n_sediment_fields()

    ! Returns the number of sediment fields

    integer :: get_n_sediment_fields

    get_n_sediment_fields = option_count('/material_phase[0]/sediment/scalar_field')

  end function get_n_sediment_fields

  function get_sediment_field_name(i_field)

    ! Returns 1-based sediment field name

    integer, intent(in)            :: i_field
    character(len=FIELD_NAME_LEN)  :: get_sediment_field_name

    call get_option('/material_phase[0]/sediment/scalar_field['//int2str(i_field - 1)&
         &//']/name', get_sediment_field_name) 

  end function get_sediment_field_name

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

  subroutine set_sediment_reentrainment(state)

    type(state_type), intent(in)                           :: state
    type(scalar_field_pointer), dimension(:), allocatable  :: bedload
    type(scalar_field), pointer                            :: field, erosion
    type(scalar_field)                                     :: bedload_surface
    type(vector_field), pointer                            :: bed_shear_stress
    type(vector_field)                                     :: shear_stress_surface
    type(tensor_field)                                     :: viscosity_surface
    type(mesh_type), pointer                               :: bottom_mesh
    real, dimension(:), allocatable                        :: diameter
    real                                                   :: R, g, erosion_flux, shear,&
         & dt, rho_0, sinking_velocity
    real                                                   :: erodibility, porosity,&
         & critical_shear_stress
    logical                                                :: have_erodibility,&
         & have_porosity, have_critical_shear_stress, have_viscosity
    logical, dimension(:), allocatable                     :: have_bedload, have_diameter
    integer, dimension(:), pointer                         :: surface_element_list
    integer                                                :: i_field, i_bc, j,&
         & n_sediment_fields, nbcs, stat
    character(len=FIELD_NAME_LEN)                          :: field_name, bc_name, bc_type
    character(len=OPTION_PATH_LEN)                         :: bc_path, bc_path_i

    ewrite(1,*) "In set_sediment_reentrainment"

    ! without the bed shear stress we cannot calculate reentrainment
    if (has_vector_field(state, "BedShearStress")) then
       bed_shear_stress => extract_vector_field(state, "BedShearStress")

       call get_option("/timestepping/timestep", dt)
       call get_option("/physical_parameters/gravity/magnitude", g)

       n_sediment_fields = get_n_sediment_fields()
       
       allocate (bedload(n_sediment_fields))
       allocate (have_bedload(n_sediment_fields))
       allocate (diameter(n_sediment_fields))
       allocate (have_diameter(n_sediment_fields))

       ! get deposited sediment for each sediment field
       do i_field=1,n_sediment_fields

          ! get field and field name
          field_name = get_sediment_field_name(i_field)
          field => extract_scalar_field(state, trim(field_name))

          ! get field bedload
          bedload(i_field)%ptr => extract_scalar_field(state,trim(field_name)//"SedimentBe&
               &dload", stat=stat)
          if (stat .eq. 0) then 
             have_bedload(i_field) = .true.
          else
             have_bedload(i_field) = .false.
          end if

          ! get sediment diameter
          if (have_option(trim(field%option_path)//"/prognostic/diameter")) then
             call get_option(trim(field%option_path)//"/prognostic/diameter",diameter)
             have_diameter(i_field) = .true.
          else
             have_diameter(i_field) = .false.
          end if

       end do

       do i_field=1,n_sediment_fields

          ! get field and field name
          field_name = get_sediment_field_name(i_field)
          field => extract_scalar_field(state, trim(field_name))

          ! get boundary condition path and number of boundary conditions
          bc_path=trim(field%option_path)//'/prognostic/boundary_conditions'
          nbcs=option_count(trim(bc_path))

          ! Loop over boundary conditions for field
          boundary_conditions: do i_bc=0, nbcs-1

             bc_path_i=trim(bc_path)//"["//int2str(i_bc)//"]"

             ! Get name and type of boundary condition
             call get_option(trim(bc_path_i)//"/name", bc_name)
             call get_option(trim(bc_path_i)//"/type[0]/name", bc_type)

             ! Only reentrainment boundaries are handled here
             if (.not. (trim(bc_type) .eq. "sediment_reentrainment")) then
                cycle boundary_conditions
             end if

             ewrite(1,*) "Setting reentrainment bc for field: ",trim(field_name)

             call get_erosion_parameters()

             call j_hill_reentrainment()

             call deallocate(bedload_surface)
             call deallocate(shear_stress_surface)

          end do boundary_conditions

       end do

       ! deallocate
       deallocate (bedload)
       deallocate (have_bedload)
       deallocate (diameter)

    end if

    contains

      subroutine get_erosion_parameters()

        type(tensor_field), pointer        :: viscosity

        call get_boundary_condition(field, name=bc_name, type=bc_type, &
             surface_mesh=bottom_mesh,surface_element_list=surface_element_list)   

        ! allocate surface fields for deposited sediment, shear stress and viscosity
        call allocate(bedload_surface, bottom_mesh, name="bedload_surface")
        call allocate(shear_stress_surface, bed_shear_stress%dim, bottom_mesh, name="sh&
             &ear_stress_surface")

        ! remap deposited sediment and shear stress fields to boundary condition
        ! surface mesh
        call remap_field_to_surface(bed_shear_stress, shear_stress_surface, &
             surface_element_list)
        call remap_field_to_surface(bedload(i_field)%ptr, bedload_surface, &
             surface_element_list)
        
        call get_option(trim(field%option_path)//"/prognostic/submerged_specific_gravity"&
             &,R)

        erosion => extract_surface_field(field, bc_name, "value")
        
        call get_option(trim(field%option_path)//"/prognostic/erodability", erodibility, default&
             &=1.0)

        call get_option(trim(field%option_path)//"/prognostic/porosity", porosity,&
             & default=0.3)

        call get_option('/material_phase::'//trim(state%name)//'/equation_of_state/fluids/&
             &linear/reference_density', rho_0)       

        call get_option(trim(field%option_path)//"/prognostic/SinkingVelocity", sinking_velocity)       

        if (have_option(trim(field%option_path)//"/prognostic/critical_shear_stress")) then
           viscosity => extract_tensor_field(state, "Viscosity")
           ! call allocate(viscosity_surface, mesh=bottom_mesh, dim=viscosity%dim, name="viscosity_sh&
           !      &ear_stress")
           ! call remap_field_to_surface(viscosity, viscosity_surface, surface_element_list) 
           
           have_viscosity = .true.
        else
           have_viscosity = .false.
        end if

        if (have_option(trim(field%option_path)//"/prognostic/diameter")) then
           call get_option(trim(field%option_path)//"/prognostic/diameter",diameter)
           have_diameter = .true.
        else
           have_diameter = .false.
        end if

        if (have_option(trim(field%option_path)//"/prognostic/critical_shear_stress")) then
           call get_option(trim(field%option_path)//"/prognostic/critical_shear_stress",&
                & critical_shear_stress)
           have_critical_shear_stress = .true.
        else
           have_critical_shear_stress = .false.
        end if

      end subroutine get_erosion_parameters

      subroutine Garcia_1991_reentrainment()

        real, dimension(:,:), allocatable  :: viscosity_node_val
        real                               :: R_p, u_star, d_50, A, sigma, Z
        integer                            :: node, row, column

        if (.not. (all(have_diameter) .and. have_viscosity)) then
           FLExit("All sediment fields must have a diameter, and a viscosity field must be&
                & specified to calculate erosion using the Garcia_1991 formula")
        end if

        ! allocate (viscosity_node_val(viscosity_surface%dim, viscosity_surface%dim))

        A = 1.3*10**(-7)
        
        do node = 1, node_count(bedload_surface)

           ! check viscosity is isotropic (NOT DOING THIS AS IT MIGHT BE SLOW - not in schema)
           ! viscosity_node_val = node_val(viscosity_surface, node)
           ! do row = 1, viscosity_surface%dim
           !    do column = 1, viscosity_surface%dim
           !       if (row .eq. column) then
           !          if (.not. (viscosity_node_val(row, column) .eq. viscosity_node_val(1,&
           !               & 1))) then
           !             FLExit("Garcia_1991 entrainment is only valid for isotropic viscosi&
           !                  &ty fields")
           !          end if
           !       else
           !          if (.not. (viscosity_node_val(row, column) .eq. 0.0)) then
           !             FLExit("Garcia_1991 entrainment is only valid for isotropic viscosi&
           !                  &ty fields")
           !          end if
           !       end if
           !    end do
           ! end do

           ! calculate particle Reynolds number
           R_p = (R*g*(diameter(i_field)/1000.0)**3)**0.5/viscosity_node_val(1,1)

           ! calculate u_star (shear velocity)
           u_star = (node_val(bedload_surface, node)/rho_0)**0.5

           ! calculate d_50 (median grain size by volume)
           d_50 = 0

           ! calculate sigma (standard deviation of bed sediment)
           sigma = 0

           ! calculate Z
           Z = (1 - 0.288 * sigma) * u_star / sinking_velocity * R_p**0.6 * &
                & ( diameter(i_field) / d_50 * 1000.0 )**0.2

           ! calculate erosion
           erosion_flux = A*Z**5 / (1 + A*Z**5/0.3)
           ! A limit is placed depending on how much of that sediment is in the
           ! bedload
           if (erosion_flux*dt > node_val(bedload_surface, node)) then
              erosion_flux = node_val(bedload_surface, node)/dt
           end if
           call set(erosion, node, erosion_flux)

        end do

        deallocate (viscosity_node_val)

      end subroutine Garcia_1991_reentrainment

      subroutine j_hill_reentrainment()

        ! get or calculate critical shear stress
        if ((.not. have_critical_shear_stress) .and. have_diameter(i_field)) then
           ! estimate of critical shear stress assuming grains larger than
           ! 10 microns and constant viscosity - note the conversion to mm!
           critical_shear_stress = 0.041 * R * 1024. * g * (diameter(i_field)/1000.)
        else
           FLExit("You need to either specify a critical shear stress or a &&
                && sediment diameter to use the j_hill formula for erosion")
        end if

        ! calculate eroded sediment flux and set reentrainment BC
        ! we only need to add to the source the erosion of sediment from the
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
        do j=1,node_count(bedload_surface)

           shear = norm2(node_val(shear_stress_surface, j))
           ! critical stress is either given by user (optional) or calculated
           ! using Shield's formula (depends on grain size and density and
           ! (vertical) viscosity)
           erosion_flux = erodibility*(1-porosity)*((shear - critical_shear_stress) /&
                & critical_shear_stress)

           if (erosion_flux < 0) then 
              erosion_flux = 0.0
           end if
           ! A limit is placed depending on how much of that sediment is in the
           ! bedload
           if (erosion_flux*dt > node_val(bedload_surface,j)) then
              erosion_flux = node_val(bedload_surface,j)/dt
           end if

           call set(erosion,j,erosion_flux)

        end do

      end subroutine j_hill_reentrainment

  end subroutine set_sediment_reentrainment

end module sediment
