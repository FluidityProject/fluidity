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
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module additional_field_locations

    use global_parameters, only: OPTION_PATH_LEN
    
    implicit none
      

    public :: additional_fields_absolute, additional_fields_relative, grandchild_paths 
    public :: additional_diagnostic_paths
    
    !! A list of locations in which additional scalar/vector/tensor fields
    !! are to be found. These are absolute paths in the schema.
    character(len=OPTION_PATH_LEN), dimension(8) :: additional_fields_absolute=&
       (/ &
       "/ocean_biology/pznd                                                                                                   ", &
       "/ocean_biology/six_component                                                                                          ", &
       "/ocean_forcing/iceshelf_meltrate/Holland08                                                                            ", &
       "/ocean_forcing/bulk_formulae/output_fluxes_diagnostics                                                                ", &
       "/porous_media                                                                                                         ", &
       "/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/les_model/dynamic_les ", &
       "/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/les_model/second_order", &
       "/material_phase[0]/sediment/                                                                                          " &
       /)
       
    !! A list of relative paths under /material_phase[i]
    !! that are searched for additional fields to be added.
    character(len=OPTION_PATH_LEN), dimension(13) :: additional_fields_relative=&
       (/ &
       "/subgridscale_parameterisations/Mellor_Yamada                                                       ", &
       "/subgridscale_parameterisations/prescribed_diffusivity                                              ", &
       "/subgridscale_parameterisations/GLS                                                                 ", &
       "/subgridscale_parameterisations/k-epsilon                                                           ", &
       "/subgridscale_parameterisations/k-epsilon/debugging_options/source_term_output_fields               ", &
       "/subgridscale_parameterisations/k-epsilon/debugging_options/prescribed_source_terms                 ", &
       "/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/les_model/second_order", &
       "/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/les_model/fourth_order", &
       "/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/les_model/wale        ", &
       "/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/les_model/dynamic_les ", &
       "/vector_field::Velocity/prognostic/equation::ShallowWater                                           ", &
       "/vector_field::Velocity/prognostic/equation::ShallowWater/bottom_drag                               ", &
       "/vector_field::BedShearStress/diagnostic/calculation_method/velocity_gradient                       " &
       /)

    !! Relative paths under a field that are searched for grandchildren
    !! (moved here because of extremely obscure intel ICE -Stephan)
    character(len=OPTION_PATH_LEN), dimension(1):: &
         grandchild_paths = (/&
         &    "/spatial_discretisation/inner_element" &
         /)

    !! List of other fields that Populate_State handles under the hood, but are
    !! needed for diagnostic dependencies.
    character(len=OPTION_PATH_LEN), dimension(1):: additional_diagnostic_paths = &
        (/ &
         &    "/geometry/ocean_boundaries" &
        /)


end module additional_field_locations
