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

  public set_sediment_reentrainment

  private

contains

subroutine set_sediment_reentrainment(state)

    type(state_type), intent(in):: state

    type(scalar_field), pointer        :: SedConc, erosion, bedload
    type(vector_field), pointer        :: bedShearStress
    real                               :: erodibility, porosity, critical_shear_stress, shear
    real                               :: erosion_flux, diameter, density, g, s, S_star
    real                               :: viscosity
    integer                            :: nSediments, nNodes, i, j, k, id
    type(scalar_field)                 :: shear_mag
    integer, dimension(:), allocatable :: faceglobalnodes
    integer                            :: snloc,ele,sele,globnod
    character(len=FIELD_NAME_LEN)      :: class_name
    character(len=OPTION_PATH_LEN)     :: option_path

    bedShearStress => extract_vector_field(state, "BedShearStress")


    nSediments = option_count('/material_phase::'//trim(state%name)&
            //'/sediment/sediment_class')
    call get_option('/material_phase::'//trim(state%name)&
           //"ScalarField::SedimentTemplate/porosity", porosity)
    call get_option("/physical_parameters/gravity/magnitude", g)
    viscosity = 1. / 1000.

    call allocate(shear_mag,bedShearStress%mesh,"ShearMagnitude") 

    do i=1,nSediments

        option_path='/material_phase::'//trim(state%name)//&
               '/sediment/sediment_class['//int2str(i-1)//"]"
          
        call get_option(trim(option_path)//"/name", class_name)
        call get_option(trim(option_path)//"/erodibility", erodibility)
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

        SedConc => extract_scalar_field(state, "SedimentConcentration"//trim(class_name))
        bedload => extract_scalar_field(state,"SedimentFlux"//trim(class_name))
        call sediment_get_boundary_condition_ids(option_path//"/prognostic/boundary_conditions",id)
        if (id .eq. -1) then
            cycle
        end if
        erosion => extract_surface_field(SedConc, id,"value")
        nNodes = node_count(erosion)
    
        ! we only need to add to the source the erosion of sediment from the
        ! bedload into the neumann BC term
        !
        ! The source depends on the erodability of that sediment
        ! (usually 1 unless you want to do something like mud which sticks),
        ! the porosity in the bedload and the bed shear stress. 
        !
        ! Each sediment class has a critical shear stress, which if exceeded
        ! by the bed shear stress, sediment is placed into suspension

        snloc = face_loc(erosion, 1)
        allocate( faceglobalnodes(1:snloc) )
        ! loop over nodes in bottom surface
        do j=1,ele_count(erosion)

            
            faceglobalnodes = face_global_nodes(erosion, sele)
            do k = 1,snloc
                globnod = faceglobalnodes(k)
                shear = norm2(node_val(bedShearStress, globnod))

                ! critical stress is either given by user (optional) or calculated
                ! using Shield's formula (depends on grain size and density and
                ! (vertical) viscosity)
                erosion_flux = erodibility*(1-porosity)*((shear - critical_shear_stress) / critical_shear_stress)

                ! A limit is placed depending on how much of that sediment is in the
                ! bedload
                if (erosion_flux > node_val(bedload,globnod)) then
                    erosion_flux = node_val(bedload,globnod)
                end if

                call set(erosion,globnod,erosion_flux)
            end do
        end do

    end do

    ! leave it up to the standard solve to add this term in as a neumann BC.

    call deallocate(shear_mag)


end subroutine set_sediment_reentrainment



!----------------------------!
!     PRIVATE ROUTINES       !
!----------------------------!

subroutine sediment_get_boundary_condition_ids(bc_path, id)

    character(len=OPTION_PATH_LEN), intent(in)  :: bc_path
    integer, intent(out)                        :: id

    integer                                     :: nBCs, i
    character(len=OPTION_PATH_LEN)              :: bc_path_i
    character(len=FIELD_NAME_LEN)               :: bc_type 
    
    id = -1
    nbcs=option_count(trim(bc_path))
    ! Loop over boundary conditions
    do i=0, nbcs-1
        bc_path_i=trim(bc_path)//"["//int2str(i)//"]"
        call get_option(trim(bc_path_i)//"/type[0]/name", bc_type)
        if (trim(bc_type) .eq. "sediment_reentrainment") then
            id = i+1
        end if
    end do
    
end subroutine sediment_get_boundary_condition_ids

end module sediment
