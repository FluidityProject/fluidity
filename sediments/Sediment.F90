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

  public set_sediment_reentrainment, sediment_init, get_n_sediment_fields, &
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

subroutine sediment_init()

    

    ! Options check
    ! Chesk sinking velocity is specified for every sediment
    ! Check continuity is the same for all deposition fields
    ! Deposition on NoFlux or Zero Dirichlet boundary
    ! Reentrainment on non-deposition boundary

end subroutine sediment_init

subroutine set_sediment_reentrainment(state)

    type(state_type), intent(in):: state

    type(scalar_field), pointer        :: field, erosion, bedload
    type(vector_field), pointer        :: bed_shear_stress
    real                               :: erodibility, porosity
    real                               :: critical_shear_stress, shear
    real                               :: erosion_flux, diameter, density, g, s
    real                               :: viscosity
    integer                            :: i_field, i_bc, j, n_sediment_fields, nbcs
    character(len=FIELD_NAME_LEN)      :: field_name, bc_name, bc_type
    character(len=OPTION_PATH_LEN)     :: bc_path, bc_path_i
    type(scalar_field)                 :: bedload_surface
    type(vector_field)                 :: shear_stress_surface
    integer, dimension(:), pointer     :: surface_element_list
    type(mesh_type), pointer           :: bottom_mesh
    real                               :: dt

    ewrite(1,*) "In set_sediment_bc"
             
    bed_shear_stress => extract_vector_field(state, "BedShearStress")

    call get_option("/timestepping/timestep", dt)
    call get_option("/physical_parameters/gravity/magnitude", g)
    viscosity = 1. / 1000.

    n_sediment_fields = get_n_sediment_fields()

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
          
          ewrite(2,*) "Checking field BC: ",trim(field%name),' ',trim(bc_name),' ',&
               trim(bc_type)

          ! get deposited sediment field
          bedload => extract_scalar_field(state,trim(field_name)//"SedimentDeposition")

          call get_boundary_condition(field, name=bc_name, type=bc_type, &
               surface_mesh=bottom_mesh,surface_element_list=surface_element_list)   
      
          ! allocate surface fields for deposited sediment and shear stress
          call allocate(bedload_surface, bottom_mesh, name="bedload_surface")
          call allocate(shear_stress_surface, bed_shear_stress%dim, bottom_mesh, name="sh&
               &ear_stress_surface")

          ! remap deposited sediment and shear stress fields to boundary condition
          !  surface mesh
          call remap_field_to_surface(bed_shear_stress, shear_stress_surface, &
               surface_element_list)
          call remap_field_to_surface(bedload, bedload_surface, &
               surface_element_list)
       
          ! get or calculate critical shear stress
          call get_option(trim(field%option_path)//"/prognostic/erodability", erodibility, default=1.0)
          if (have_option(trim(field%option_path)//"/prognostic/critical_shear_stress")) then
             call get_option(trim(field%option_path)//"/prognostic/critical_shear_stress", critical_shear_stress)
          else
             if (have_option(trim(field%option_path)//"/prognostic/diameter")) then
                call get_option(trim(field%option_path)//"/prognostic/diameter",diameter)
             else
                FLExit("You need to either specify a critical shear stress or a &&
                     && sediment diameter")
             end if
             call get_option(trim(field%option_path)//"/prognostic/density",density)
             ! calc critical shear stress
             s = density/1000.
             !S_star = sqrt((s-1)*g*diameter**3)/viscosity
             !critical_shear_stress = 0.105*S_star**(-0.13) + &
             !                        0.045*exp(-35*S_star**(-0.59))
             ! estimate of critical shear stress assuming grains larger than
             ! 10 microns and constant viscosity - note the conversion to mm!
             critical_shear_stress = 0.041 * (s-1) * 1024. * g * (diameter/1000.)
          end if

          ! calculate eroded sediment flux and set reentrainment BC
          call get_option(trim(field%option_path)//"/prognostic/porosity", porosity, default=0.3)
          erosion => extract_surface_field(field, bc_name, "value")
          ! we only need to add to the source the erosion of sediment from the
          ! bedload into the neumann BC term
          !
          ! The source depends on the erodability of that sediment
          ! (usually 1 unless you want to do something like mud which sticks),
          ! the porosity in the bedload and the bed shear stress. 
          !
          ! Each sediment class has a critical shear stress, which if exceeded
          ! by the bed shear stress, sediment is placed into suspension

          ! loop over nodes in bottom surface
          do j=1,node_count(bedload_surface)

             shear = norm2(node_val(shear_stress_surface, j))
             ! critical stress is either given by user (optional) or calculated
             ! using Shield's formula (depends on grain size and density and
             ! (vertical) viscosity)
             erosion_flux = erodibility*(1-porosity)*((shear - critical_shear_stress) / critical_shear_stress)

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
          
          call deallocate(bedload_surface)
          call deallocate(shear_stress_surface)

       end do boundary_conditions

    end do

end subroutine set_sediment_reentrainment

end module sediment
