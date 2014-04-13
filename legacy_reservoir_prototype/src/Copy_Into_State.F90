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
         proto_densities, &
         proto_components, &
         ncomp, &
         nphase, &
         cv_ndgln, &
         p_ndgln )

      !!< Copy prototype solution arrays into fluidity state array for output

      type(state_type), dimension(:), intent(inout) :: state
      real, dimension(:), intent(in) :: proto_saturations
      real, dimension(:), intent(in) :: proto_temperatures
      real, dimension(:), intent(in) :: proto_densities
      real, dimension(:), intent(in) :: proto_components
      integer, intent(in) :: ncomp
      integer, intent(in) :: nphase
      integer, dimension(:), intent(in) :: cv_ndgln
      integer, dimension(:), intent(in) :: p_ndgln

      ! local variables
      integer :: stat
      integer :: i,j,k,p,ele,jloc, pos
      integer :: number_nodes
      integer :: nstates
      integer :: nloc
      integer, dimension(:), pointer :: element_nodes => null()
      type(scalar_field), pointer :: phasevolumefraction => null()
      type(scalar_field), pointer :: phasetemperature => null()      
      type(scalar_field), pointer :: density => null()
      type(scalar_field), pointer :: componentmassfraction => null()
      character(len=option_path_len) :: material_phase_name

      ewrite(3,*) "In copy_into_state"

      assert(size(state) >= nphase)

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

      end do phase_loop

      ewrite(3,*) "Leaving copy_into_state"

      return

    end subroutine copy_into_state

    function extract_or_insert_velocity_tensor(states,phase,position,velocity,pressure) result(vfields)

      type(state_type), dimension(:), target :: states
      integer :: phase
      type(vector_field) :: position,velocity
      type(scalar_field) :: pressure

      type(vector_field_pointer), dimension(ele_loc(pressure,1))  :: vfields
      type(vector_field), dimension(ele_loc(pressure,1)) :: vfield_array
      type(mesh_type) :: fem_velocity_mesh
      type(element_type) :: shape
      type(mesh_type) :: mesh
      type(state_type), pointer :: state
      integer :: stat, i, ic, nics
      character(len=OPTION_PATH_LEN) :: name

      ewrite(3,*) "In extract_or_insert_velocity_tensor"

      state=>states(phase)

      ewrite(3,*) "ndim "// int2str(mesh_dim(velocity)) 
      ewrite(3,*) "nloc "// int2str(ele_loc(pressure,1))
      
      

      if (has_vector_field(state,"VelocityChunk1") )then
         do i=1,ele_loc(pressure,1)
            name="VelocityChunk"//int2str(i)
            vfields(i)%ptr=>extract_vector_field(state,trim(name))
         end do
      else

         call copy_option(trim(velocity%mesh%option_path),&
              "/geometry/mesh::FEM"//trim(velocity%mesh%name),stat)
         call set_option("/geometry/mesh::FEM"//trim(velocity%mesh%name)//&
              "/from_mesh/mesh_shape/element_type",&
              "lagrangian")

         if (has_mesh(state,"FEM"//trim(velocity%mesh%name))) then
            mesh=extract_mesh(states(1),"FEM"//trim(velocity%mesh%name))
         else
            shape=ele_shape(velocity,1)
            mesh=make_mesh(model=pressure%mesh,shape=shape,continuity=-1,&
            name="FEM"//trim(velocity%mesh%name))
            mesh%option_path="/geometry/mesh::FEM"//trim(velocity%mesh%name)
            call insert(state,mesh,"FEM"//trim(velocity%mesh%name))
         end if

         call add_option(trim(complete_field_path(velocity%option_path))&
              //"/exclude_from_checkpointing",stat)

         do i=1,ele_loc(pressure,1)
            name="VelocityChunk"//int2str(i)
            call allocate(vfield_array(i),dim=mesh_dim(velocity),&
              mesh=mesh,name=trim(name))
            call add_option(trim(state%option_path)&
                 //"/vector_field::"//trim(name)//"/prognostic/mesh",stat)
            call add_option(trim(state%option_path)&
               //"/vector_field::"//trim(name)//&
               "/prognostic/initial_condition"&
               ,stat)
            call add_option(trim(state%option_path)&
               //"/vector_field::"//trim(name)//&
               "/prognostic/output/exclude_from_vtu"&
               ,stat)
            call set_option(trim(state%option_path)&
              //"/vector_field::"//trim(name)//"/prognostic/mesh/name",&
              "FEM"//trim(velocity%mesh%name),stat)
         
            vfield_array(i)%option_path=trim(state%option_path)&
                 //"/vector_field::"//trim(name)

            call insert(state,vfield_array(i),trim(name))
            call deallocate(vfield_array(i))
            vfields(i)%ptr=>extract_vector_field(state,trim(name))
         end do
      end if

      ewrite(3,*) "Leaving extract_or_insert_velocity_tensor"
      

    end function extract_or_insert_velocity_tensor


  end module Copy_BackTo_State
