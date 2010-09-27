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
!    C.Pain@Imperial.ac.uk
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
  use equation_of_state
  use state_fields_module
  use boundary_conditions
  use FLDebug

  implicit none

  integer, dimension(:), allocatable :: sediment_boundary_condition_ids 
  integer                            :: nSedimentClasses

  public set_sediment_reentrainment, sediment_init, set_sediment_bc_id
  public sediment_cleanup, get_sediment_bc_id

  private

contains

subroutine sediment_init()

    nSedimentClasses = option_count('/material_phase[0]/sediment/sediment_class')
    allocate(sediment_boundary_condition_ids(nSedimentClasses))
    sediment_boundary_condition_ids = -1

end subroutine sediment_init

subroutine sediment_cleanup

    deallocate(sediment_boundary_condition_ids)

end subroutine sediment_cleanup

subroutine set_sediment_bc_id(name, id)

    character(len=FIELD_NAME_LEN), intent(in)  :: name
    integer, intent(in)                        :: id

    integer                        :: i
    character(len=FIELD_NAME_LEN)  :: class_name
    character(len=OPTION_PATH_LEN) :: option_path

    do i=1,nSedimentClasses

        option_path='/material_phase[0]/sediment/sediment_class['//int2str(i-1)//"]"
        
        call get_option(trim(option_path)//"/name", class_name)

        if ("SedimentConcentration"//class_name .eq. name) then
            sediment_boundary_condition_ids(i) = id
            exit
        end if
    end do

end subroutine set_sediment_bc_id

function get_sediment_bc_id(i) result (id)
    
    integer, intent(in) :: i
    integer             :: id


    id = sediment_boundary_condition_ids(i)

end function get_sediment_bc_id

subroutine set_sediment_reentrainment(state)

    type(state_type), intent(in):: state

    type(scalar_field), pointer        :: SedConc, erosion, bedload
    type(vector_field), pointer        :: bedShearStress
    real                               :: erodibility, porosity, critical_shear_stress, shear
    real                               :: erosion_flux, diameter, density, g, s, S_star
    real                               :: viscosity
    integer                            :: nNodes, i, j
    character(len=FIELD_NAME_LEN)      :: class_name
    character(len=OPTION_PATH_LEN)     :: option_path
    type(scalar_field)                 :: bedLoadSurface
    type(vector_field)                 :: shearStressSurface
    integer, dimension(:), pointer     :: surface_element_list
    type(mesh_type), pointer           :: bottom_mesh
    logical                            :: alloced
    real                               :: dt

    ewrite(1,*) "In set_sediment_bc"

    bedShearStress => extract_vector_field(state, "BedShearStress")

    call get_option("/timestepping/timestep", dt)

    call get_option("/material_phase[0]/sediment/scalar_field::SedimentTemplate/porosity", porosity, default=0.3)
    call get_option("/physical_parameters/gravity/magnitude", g)
    viscosity = 1. / 1000.

    alloced = .false.

    do i=1,nSedimentClasses

        if (sediment_boundary_condition_ids(i) .eq. -1) then
            cycle
        end if

        option_path='/material_phase::'//trim(state%name)//&
               '/sediment/sediment_class['//int2str(i-1)//"]"
          
        call get_option(trim(option_path)//"/name", class_name)

        SedConc => extract_scalar_field(state, "SedimentConcentration"//trim(class_name))
        bedload => extract_scalar_field(state,"SedimentFlux"//trim(class_name))

        if (.not. alloced) then
            call get_boundary_condition(SedConc, sediment_boundary_condition_ids(i), &
                        surface_mesh=bottom_mesh,surface_element_list=surface_element_list)         
            !call create_surface_mesh(bottom_mesh, surface_nodes, bedShearStress%mesh, surface_element_list, 'ErosionBed')
            call allocate(bedLoadSurface,bottom_mesh, name="bedLoadSurface")
            call allocate(shearStressSurface,bedShearStress%dim,bottom_mesh, name="shearStressSurface")
            call remap_field_to_surface(bedShearStress, shearStressSurface, &
                                        surface_element_list)
            
            nNodes = node_count(bedLoadSurface)
            alloced = .true.
        end if
        call get_option(trim(option_path)//"/erodibility", erodibility, default=1.0)
        call get_option(trim(option_path)//"/density",density)
        call get_option(trim(option_path)//"/diameter",diameter)
        if (have_option(trim(option_path)//"/critical_shear_stress")) then
            call get_option(trim(option_path)//"/critical_shear_stress", critical_shear_stress)
        else
            ! calc critical shear stress
            s = density/1024.
            S_star = sqrt((s-1)*g*diameter**3)/viscosity
            critical_shear_stress = 0.105*S_star**(-0.13) + &
                                    0.045*exp(-35*S_star**(-0.59))
        end if
        
        
        call remap_field_to_surface(bedload, bedLoadSurface, &
                                        surface_element_list)


        erosion => extract_surface_field(SedConc,sediment_boundary_condition_ids(i),"value")
        
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
        do j=1,nNodes
          
            shear = norm2(node_val(ShearStressSurface, j))
            ! critical stress is either given by user (optional) or calculated
            ! using Shield's formula (depends on grain size and density and
            ! (vertical) viscosity)
            erosion_flux = erodibility*(1-porosity)*((shear - critical_shear_stress) / critical_shear_stress)
                 
            if (erosion_flux < 0) then 
                erosion_flux = 0.0
            end if
            ! A limit is placed depending on how much of that sediment is in the
            ! bedload
            if (erosion_flux*dt > node_val(bedLoadSurface,j)) then
                erosion_flux = node_val(bedLoadSurface,j)/dt
            end if

            call set(erosion,j,erosion_flux)
       
        end do

    end do

    ! leave it up to the standard solve to add this term in as a neumann BC.
    if (alloced) then
        call deallocate(bedLoadSurface)
        call deallocate(shearStressSurface)
    end if
end subroutine set_sediment_reentrainment

end module sediment
