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

  integer, dimension(:), allocatable :: sediment_boundary_condition_ids

  public set_sediment_reentrainment, sediment_init, set_sediment_bc_id
  public sediment_cleanup, get_sediment_bc_id
  public get_sediment_name, get_nSediments

  private

contains

function get_sediment_name(i)
    !!< Return the name of sediment class from the options tree at position i-1.
    !!< Loops over sediment should therefore be 1,nSedimentClasses in order to
    !!< make this funciton work.

    integer, intent(in) :: i

    character(len=FIELD_NAME_LEN) :: get_sediment_name

    call get_option('/material_phase[0]/sediment/class['//int2str(i-1)//']/name',get_sediment_name)        

end function get_sediment_name


function get_nSediments()
    !!< Return the number of sediment classes
    
    integer :: get_nSediments

    get_nSediments = option_count('/material_phase[0]/sediment/class')

end function get_nSediments

subroutine sediment_init()

    integer                        :: ic, bc, i, stat, nSedimentClasses
    integer                        :: nBoundaryConditions, nInitialConditions
    character(len=OPTION_PATH_LEN) :: option_path, temp_path
    character(len=FIELD_NAME_LEN)  :: class_name

    ! Set up some constants and allocate an array of sediment boundary condition
    ! IDs. 
    nSedimentClasses = option_count('/material_phase[0]/sediment/sediment_class')
    allocate(sediment_boundary_condition_ids(nSedimentClasses))
    sediment_boundary_condition_ids = -1

    ewrite(2,*) "In sediment_init"
    ewrite(2,*) "Found ",nSedimentClasses, "sediment classes"
    ! Check if we have SedimentTemplate on - if so, this is a fresh run, so 
    ! construct the sediment classes fro mthe template.
    ! If it doesn't exist, it's a checkpoint run, so no need to do this.
    ! Note that we have to also check options more carefully, as we can't 
    ! carry out all possible checks in options_check (our fields haven't been
    ! constructed at that point).
    if (have_option('/material_phase[0]/sediment/scalar_field::SedimentTemplate')) then

        ewrite(2,*) "Initialising sediment classes from Template"

        ! For each sediment class
        sediment_classes: do i=1,nSedimentClasses
            ! copy the template into a temp option path
            temp_path = '/material_phase[0]/sediment/Sediment_'//int2str(i-1)//'_temp'
            call copy_option('/material_phase[0]/sediment/scalar_field::SedimentTemplate',&
                             trim(temp_path))

            ! now go over all the top level options and copy them into the 
            ! temp option tree
            option_path='/material_phase[0]/sediment/sediment_class['//int2str(i-1)//']'
            call get_option(trim(option_path)//'/name',class_name)
            option_path='/material_phase[0]/sediment/sediment_class::'//trim(class_name)
            call delete_option(trim(temp_path)//"/name")
            call copy_option(trim(option_path)//"/name",trim(temp_path)//"/name")

            ! loop over all boundary conditions
            nBoundaryConditions = option_count(trim(option_path)//'/boundary_condition')
            bc_loop: do bc=1,nBoundaryConditions
                call move_option(trim(option_path)//'/boundary_condition['//int2str(bc)//']',&
                                 trim(temp_path)//'/boundary_condition['//int2str(bc)//']')
            call delete_option(trim(option_path),stat)
            end do bc_loop
                          
            ! loop over all initial conditions
            nInitialConditions = option_count(trim(option_path)//'/initial_condition')
            ic_loop: do ic=1,nInitialConditions
                call move_option(trim(option_path)//'/initial_condition['//int2str(ic)//']',&
                                 trim(temp_path//'/initial_condition['//int2str(ic)//']'))
            end do ic_loop

            ! now do the one off fields and options
            if (have_option(trim(option_path)//'/tensor_field::Diffusivity')) then
                if (have_option(trim(temp_path)//'/prognostic/tensor_field::Diffusivity')) then
                    call delete_option(trim(temp_path)//'/prognostic/tensor_field::Diffusivity')
                end if
                call move_option(trim(option_path)//'/tensor_field::Diffusivity',&
                & trim(temp_path)//'/prognostic/tensor_field::Diffusivity')
            end if

            ! Check that we have a density somewhere
            if (.not. have_option(trim(option_path)//'/density') .and. &
                .not. have_option(trim(temp_path)//'/density')) then
                FLExit("You must specify a sediment density under the Template &&
                && or under each sediment class")
            end if
            ! OK, so let's make sure we have the one under the sediment class
            ! if it exists.
            if (have_option(trim(option_path)//'/density')) then
                if (have_option(trim(temp_path)//'/density')) then
                    call delete_option(trim(temp_path)//'/density')
                end if
                call move_option(trim(option_path)//'/density',&
                & trim(temp_path)//'/density')
            end if

            if (have_option(trim(option_path)//'/diameter')) then
                if (have_option(trim(temp_path)//'/diameter')) then
                    call delete_option(trim(temp_path)//'/diameter')
                end if
                call move_option(trim(option_path)//'/diameter',&
                & trim(temp_path)//'/diameter')
            end if
            
            if (have_option(trim(option_path)//'/erodability')) then
                if (have_option(trim(temp_path)//'/erodability')) then
                    call delete_option(trim(temp_path)//'/erodability')
                end if
                call move_option(trim(option_path)//'/erodability',&
                & trim(temp_path)//'/erodability')
            end if
            
            if (have_option(trim(option_path)//'/critical_shear_stress')) then
                if (have_option(trim(temp_path)//'/critical_shear_stress')) then
                    call delete_option(trim(temp_path)//'/critical_shear_stress')
                end if
                call move_option(trim(option_path)//'/critical_shear_stress',&
                & trim(temp_path)//'/critical_shear_stress')
            end if
            
            if (have_option(trim(option_path)//'/scalar_field::SinkingVelocity')) then
                if (have_option(trim(temp_path)//'/prognostic/scalar_field::SinkingVelocity')) then
                    call delete_option(trim(temp_path)//'/prognostic/scalar_field::SinkingVelocity')
                end if
                call move_option(trim(option_path)//'/scalar_field::SinkingVelocity',&
                & trim(temp_path)//'/prognostic/scalar_field::SinkingVelocity')
            end if

            ! All done, move over the temporary one to the main material phase
            ! Note that this is no longer under /sediment, so it acts like a
            ! normal scalar field from now on...
            call move_option(temp_path,'/material_phase[0]/scalar_field::'//trim(class_name))
            ! ...however, after a restrat from checkpoint, we won't know about
            ! the sediment fields any longer, so add them to the option tree
            call add_option('/material_phase[0]/sediment/class['//int2str(i-1)//']',stat)
            call add_option('/material_phase[0]/sediment/class['//int2str(i-1)//']/name',stat)
            call set_option('/material_phase[0]/sediment/class['//int2str(i-1)//']/name',trim(class_name))

        end do sediment_classes

        ! Remove the Sediment template from the options tree
        call delete_option('/material_phase[0]/sediment/scalar_field::SedimentTemplate')
        ! delete the old sediment_class options
        do i=1,nSedimentClasses
            ! only ever need to delete 0 and when we remove a field, the next
            ! one will become zero
            call delete_option('/material_phase[0]/sediment/sediment_class[0]')
        end do

    else
        ! restart from checkpoint, so check the modified options tree
        nSedimentClasses = option_count('/material_phase[0]/sediment/class/')
        ewrite(2,*) "Initialising sediment classes from pre-existing run"
        ewrite(2,*) "Found ",nSedimentClasses, "sediment classes"
        ewrite(2,*) "They live in: "
        do i=1,nSedimentClasses
            call get_option('/material_phase[0]/sediment/class['//int2str(i-1)//']/name',class_name)
            ewrite(2,*) "Sediment class ",i," ",trim(class_name)
        end do

    end if

end subroutine sediment_init

subroutine sediment_cleanup

    deallocate(sediment_boundary_condition_ids)

end subroutine sediment_cleanup

subroutine set_sediment_bc_id(name, id)

    character(len=FIELD_NAME_LEN), intent(in)  :: name
    integer, intent(in)                        :: id

    integer                        :: i, nSedimentClasses
    character(len=FIELD_NAME_LEN)  :: class_name

    nSedimentClasses = get_nSediments()

    do i=1,nSedimentClasses
        
        class_name = get_sediment_name(i)

        if (class_name .eq. name) then
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
    real                               :: erosion_flux, diameter, density, g, s
    real                               :: viscosity
    integer                            :: nNodes, i, j, nSedimentClasses
    character(len=FIELD_NAME_LEN)      :: class_name
    character(len=OPTION_PATH_LEN)     :: option_path
    type(scalar_field)                 :: bedLoadSurface
    type(vector_field)                 :: shearStressSurface
    integer, dimension(:), pointer     :: surface_element_list
    type(mesh_type), pointer           :: bottom_mesh
    logical                            :: alloced
    real                               :: dt

    ewrite(1,*) "In set_sediment_bc"

    call get_option("/timestepping/timestep", dt)

    
    call get_option("/physical_parameters/gravity/magnitude", g)
    viscosity = 1. / 1000.

    alloced = .false.

    nSedimentClasses = get_nSediments()

    do i=1,nSedimentClasses

        if (sediment_boundary_condition_ids(i) .eq. -1) then
            cycle
        end if

        class_name = get_sediment_name(i)
        SedConc => extract_scalar_field(state, trim(class_name))
        option_path = SedConc%option_path
        bedload => extract_scalar_field(state,"SedimentFlux"//trim(class_name))

        call get_option(trim(option_path)//"/porosity", porosity, default=0.3)
        call get_option(trim(option_path)//"/name", class_name)

        if (.not. alloced) then
            bedShearStress => extract_vector_field(state, "BedShearStress")
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
        call get_option(trim(option_path)//"/erodability", erodibility, default=1.0)
        if (have_option(trim(option_path)//"/critical_shear_stress")) then
            call get_option(trim(option_path)//"/critical_shear_stress", critical_shear_stress)
        else
            if (have_option(trim(option_path)//"/diameter")) then
                call get_option(trim(option_path)//"/diameter",diameter)
            else
                FLExit("You need to either specify a critical shear stress or a &&
                && sediment diameter")
            end if
            call get_option(trim(option_path)//"/density",density)
            ! calc critical shear stress
            s = density/1000.
            !S_star = sqrt((s-1)*g*diameter**3)/viscosity
            !critical_shear_stress = 0.105*S_star**(-0.13) + &
            !                        0.045*exp(-35*S_star**(-0.59))
            ! estimate of critical shear stress assuming grains larger than
            ! 10 microns and constant viscosity - note the conversion to mm!
            critical_shear_stress = 0.041 * (s-1) * 1024. * g * (diameter/1000.)
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
