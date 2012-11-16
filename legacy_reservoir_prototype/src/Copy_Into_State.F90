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
  !    but WITHOUT ANY WARRANTY; without seven the implied warranty of
  !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  !    Lesser General Public License for more details.
  !
  !    You should have received a copy of the GNU Lesser General Public
  !    License along with this library; if not, write to the Free Software
  !    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
  !    USA
  
#include "fdebug.h"
  
  module Copy_BackTo_State
    !! This module enables the multiphase prototype code to interact with state by 
    !! copying what is needed back to state for output 

    use fldebug
    use state_module
    use fields
    use field_options
    use spud
    use populate_state_module
    use diagnostic_variables
    use diagnostic_fields
    use diagnostic_fields_wrapper
    use global_parameters, only: option_path_len
    use diagnostic_fields_wrapper_new
    use element_numbering
    use shape_functions
    use fefields
    use boundary_conditions
    use futils, only: int2str

    use Copy_Outof_State, only: print_from_state

    implicit none

    private

    public :: copy_into_state

  contains

    subroutine copy_into_state( state, &
         proto_saturations, &
         proto_temperatures, &
         proto_pressure, &
         proto_velocity_u, &
         proto_velocity_v, &
         proto_velocity_w, &                               
         proto_densities, &
         proto_components, &
         ncomp, &
         nphase, &
         cv_ndgln, &
         p_ndgln, &
         u_ndgln, &
         ndim )

      !!< Copy prototype solution arrays into fluidity state array for output

      type(state_type), dimension(:), intent(inout) :: state
      real, dimension(:), intent(in) :: proto_saturations
      real, dimension(:), intent(in) :: proto_temperatures      
      real, dimension(:), intent(in) :: proto_pressure
      real, dimension(:), intent(in) :: proto_velocity_u
      real, dimension(:), intent(in) :: proto_velocity_v
      real, dimension(:), intent(in) :: proto_velocity_w      
      real, dimension(:), intent(in) :: proto_densities
      real, dimension(:), intent(in) :: proto_components
      integer, intent(in) :: ncomp      
      integer, intent(in) :: nphase
      integer, dimension(:), intent(in) :: cv_ndgln
      integer, dimension(:), intent(in) :: p_ndgln
      integer, dimension(:), intent(in) :: u_ndgln
      integer, intent(in) :: ndim

      ! local variables        
      integer :: stat
      integer :: i,j,k,p,ele,jloc
      integer :: number_nodes
      integer :: nstates
      integer :: nloc, nlev
      integer, dimension(:), pointer :: element_nodes => null()
      type(scalar_field), pointer :: phasevolumefraction => null()
      type(scalar_field), pointer :: phasetemperature => null()      
      type(scalar_field), pointer :: pressure => null()
      type(vector_field), pointer :: velocity => null()
      type(scalar_field), pointer :: density => null()
      type(scalar_field), pointer :: componentmassfraction => null()
      character(len=option_path_len) :: material_phase_name, vel_element_type
      logical :: is_overlapping
      real, dimension(:), allocatable :: proto_velocity_u_tmp, proto_velocity_v_tmp, proto_velocity_w_tmp
      ewrite(3,*) "In copy_into_state"

      assert(size(state) >= nphase)

      ! Deal with overlapping velocity...
      ! Get the vel element type.
      call get_option('/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           vel_element_type)
      is_overlapping = .false.
      if ( trim( vel_element_type ) == 'overlapping' ) is_overlapping = .true. 

      velocity => extract_vector_field(state(1), "Velocity", stat=stat)
      if (stat /= 0) then 
         FLAbort('Failed to extract phase velocity from state in copy_into_state')
      end if
      number_nodes = node_count( velocity )

      allocate( proto_velocity_u_tmp( number_nodes * nphase ) ) ; proto_velocity_u_tmp = 0.
      allocate( proto_velocity_v_tmp( number_nodes * nphase ) ) ; proto_velocity_v_tmp = 0.
      allocate( proto_velocity_w_tmp( number_nodes * nphase ) ) ; proto_velocity_w_tmp = 0.

      if ( is_overlapping ) then

         ! in case of overlapping elements just take
         ! the average of the various levels
         !
         ! this means that this field cannot be used 
         ! for any cv-fem calculations
         !
         ! ONLY for visualisation purposes

         pressure => extract_scalar_field(state(1), "Pressure")
         nlev = ele_loc( pressure, 1 )
         nloc = ele_loc( velocity, 1 )

         do p = 1, nphase
            do i = 1, element_count( velocity )
               do j = 1, nlev
                  proto_velocity_u_tmp ( (p-1)*number_nodes + (i-1)*nloc + 1 : (p-1)*number_nodes + i*nloc ) = &
                       proto_velocity_u_tmp ( (p-1)*number_nodes + (i-1)*nloc + 1 : (p-1)*number_nodes + i*nloc ) + &
                       proto_velocity_u ( (p-1)*number_nodes*nlev + (i-1)*nloc*nlev + (j-1)*nloc+ 1 : & 
                       (p-1)*number_nodes*nlev + (i-1)*nloc*nlev + j*nloc )

                  proto_velocity_v_tmp ( (p-1)*number_nodes + (i-1)*nloc + 1 : (p-1)*number_nodes + i*nloc ) = &
                       proto_velocity_v_tmp ( (p-1)*number_nodes + (i-1)*nloc + 1 : (p-1)*number_nodes + i*nloc ) + &
                       proto_velocity_v ( (p-1)*number_nodes*nlev + (i-1)*nloc*nlev + (j-1)*nloc+ 1 : & 
                       (p-1)*number_nodes*nlev + (i-1)*nloc*nlev + j*nloc )

                  proto_velocity_w_tmp ( (p-1)*number_nodes + (i-1)*nloc + 1 : (p-1)*number_nodes + i*nloc ) = &
                       proto_velocity_w_tmp ( (p-1)*number_nodes + (i-1)*nloc + 1 : (p-1)*number_nodes + i*nloc ) + &
                       proto_velocity_w ( (p-1)*number_nodes*nlev + (i-1)*nloc*nlev + (j-1)*nloc+ 1 : & 
                       (p-1)*number_nodes*nlev + (i-1)*nloc*nlev + j*nloc )
               end do
            end do
         end do

         proto_velocity_u_tmp = proto_velocity_u_tmp / nlev
         proto_velocity_v_tmp = proto_velocity_v_tmp / nlev
         proto_velocity_w_tmp = proto_velocity_w_tmp / nlev

      else
         proto_velocity_u_tmp = proto_velocity_u
         proto_velocity_v_tmp = proto_velocity_v
         proto_velocity_w_tmp = proto_velocity_w
      end if

      phase_loop: do p = 1,nphase

         phasevolumefraction => extract_scalar_field(state(p), "PhaseVolumeFraction", stat=stat)

         if (stat /= 0) then 

            ewrite(1,*) 'Issue in prototype interface for phase ',p

            FLAbort('Failed to extract phase volume fraction from state in copy_into_state')

         end if

         number_nodes = node_count(phasevolumefraction)

         volf_ele_loop: do i = 1,element_count(phasevolumefraction)

            element_nodes => ele_nodes(phasevolumefraction,i)

            nloc = size(element_nodes)

            volf_node_loop: do j = 1,nloc

               call set(phasevolumefraction, &
                    element_nodes(j), &
                    proto_saturations((cv_ndgln((i-1)*nloc+j)) + (p-1)*number_nodes))

            end do volf_node_loop

         end do volf_ele_loop

         phasetemperature => extract_scalar_field(state(p), "Temperature", stat=stat)

         found_temp: if (stat == 0) then 

            ewrite(1,*) 'In copy in to state and found temperature for phase ',p

            number_nodes = node_count(phasetemperature)

            temp_ele_loop: do i = 1,element_count(phasetemperature)

               element_nodes => ele_nodes(phasetemperature,i)

               nloc = size(element_nodes)

               temp_node_loop: do j = 1,nloc

                  call set(phasetemperature, &
                       element_nodes(j), &
                       proto_temperatures((cv_ndgln((i-1)*nloc+j)) + (p-1)*number_nodes))

               end do temp_node_loop

            end do temp_ele_loop

         end if found_temp

         density => extract_scalar_field(state(p), "Density", stat=stat)

         if (stat /= 0) then 

            ewrite(1,*) 'Issue in prototype interface for phase ',p

            FLAbort('Failed to extract phase density from state in copy_into_state')

         end if

         number_nodes = node_count(density)

         den_ele_loop: do i = 1,element_count(density)

            element_nodes => ele_nodes(density,i)

            nloc = size(element_nodes)

            den_node_loop: do j = 1,nloc

               call set(density, &
                    element_nodes(j), &
                    proto_densities((cv_ndgln((i-1)*nloc+j)) + (p-1)*number_nodes))

            end do den_node_loop

         end do den_ele_loop

         velocity => extract_vector_field(state(p), "Velocity", stat=stat)

         if (stat /= 0) then 

            ewrite(1,*) 'Issue in prototype interface for phase ',p

            FLAbort('Failed to extract phase velocity from state in copy_into_state')

         end if

         number_nodes = node_count(velocity)

         vel_ele_loop: do i = 1,element_count(velocity)

            element_nodes => ele_nodes(velocity,i)

            nloc = size(element_nodes)

            vel_node_loop: do j = 1,nloc

               ! set u
               call set(velocity, &
                    1, &
                    element_nodes(j), &
                    proto_velocity_u_tmp( (i-1)*nloc + j + (p-1)*number_nodes) )

               ! set v
               if (ndim > 1) call set(velocity, &
                    2, &
                    element_nodes(j), &
                    proto_velocity_v_tmp( (i-1)*nloc + j + (p-1)*number_nodes) )

               ! set w
               if (ndim > 2) call set(velocity, &
                    3, &
                    element_nodes(j), &
                    proto_velocity_w_tmp( (i-1)*nloc + j + (p-1)*number_nodes) )

            end do vel_node_loop

         end do vel_ele_loop

      end do phase_loop

      deallocate( proto_velocity_u_tmp, proto_velocity_v_tmp, proto_velocity_w_tmp )

      ! comp is stored in the order
      !   comp1 phase1
      !   comp1 phase2
      !   comp2 phase1
      !   comp2 phase2
      !   etc
      have_comp: if (ncomp > 0) then

         ! there is a state for each phase and each component
         nstates = size(state)

         ! Assume for now that components have been inserted in state AFTER all the phases
         comp_loop: do i = nstates-ncomp+1,nstates

            ! inspect each scalar field of this state   
            sfield_loop: do j = 1,option_count("/material_phase[" // int2str(i-1) //"]/scalar_field")

               ! extract each scalar field of this state
               componentmassfraction => extract_scalar_field(state(i), j)

               ! determine if this scalar field has the magic name for components   
               is_compfrac: if (componentmassfraction%name(1:21) == "ComponentMassFraction") then

                  ! find the phase this component fraction is associated with   
                  call get_option("/material_phase[" // int2str(i-1) //"]/scalar_field[" // int2str(j-1) //&
                       &"]/material_phase_name", material_phase_name)

                  ! find the phase index   
                  phase_index_loop: do k = 1,nphase

                     phase_name_check: if (trim(material_phase_name) == state(k)%name) then

                        number_nodes = node_count(componentmassfraction)

                        comp_ele_loop: do ele = 1,element_count(componentmassfraction)

                           element_nodes => ele_nodes(componentmassfraction,ele)

                           nloc = size(element_nodes)

                           comp_node_loop: do jloc= 1,nloc

                              call set(componentmassfraction, &
                                   element_nodes(jloc), &
                                   proto_components((cv_ndgln((ele-1)*nloc+jloc)) + ((i-(1+nphase))*nphase+(k-1))*number_nodes))

                           end do comp_node_loop

                        end do comp_ele_loop

                        cycle sfield_loop

                     end if phase_name_check

                  end do phase_index_loop

               end if is_compfrac

            end do sfield_loop

         end do comp_loop

      end if have_comp

      pressure => extract_scalar_field(state(1), "Pressure")    

      press_ele_loop: do i = 1,ele_count(pressure)

         element_nodes => ele_nodes(pressure,i)

         nloc = size(element_nodes)

         press_node_loop: do j = 1,nloc

            call set(pressure, &
                 element_nodes(j), &
                 proto_pressure(p_ndgln((i-1)*nloc+j)))

         end do press_node_loop

      end do press_ele_loop

      ewrite(3,*) "Leaving copy_into_state"

      return

    end subroutine copy_into_state



  end module Copy_BackTo_State
