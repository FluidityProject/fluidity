module darcy_transport_types

#include "fdebug.h"

  use spud
  use fields
  use state_module
  use fldebug
  use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
  
  use darcy_impes_leaching_types
  
  implicit none
  private

  public :: immobile_sfield_type, &
       darcy_impes_MIM_options_type, &
       prog_sfield_MIM_type,&
       heap_property_type


  type immobile_sfield_type
      type(scalar_field),pointer  :: sfield
      type(scalar_field),pointer  :: old_sfield
  end type immobile_sfield_type
   
  
   !the Mobile-Immobile model options
  type darcy_impes_MIM_options_type
      type(scalar_field_pointer), dimension(:), pointer :: immobile_saturation
      type(scalar_field_pointer), dimension(:), pointer :: old_immobile_saturation
      type(scalar_field_pointer), dimension(:), pointer :: mobile_saturation
      type(scalar_field_pointer), dimension(:), pointer :: old_mobile_saturation
      type(scalar_field_pointer), dimension(:), pointer :: mass_trans_coef
      type(scalar_field_pointer), dimension(:), pointer :: old_mass_trans_coef
      ! *** Flag for Whether there is Mobile-Immobile model
      logical, dimension(:), pointer :: have_MIM
      ! *** Flag to check wether the MIM exist in at least one phase
      logical :: have_MIM_phase
      ! *** Flag for Whether there is Mass transfer coefficient
      logical, dimension(:), pointer :: have_mass_trans_coef
      logical :: Lima_immobile_sat=.false. !if the immobile saturation if calculated by internal algrithm which is based on Lima,2005
      logical :: Lima_mass_trans=.false. !if the mass transfer coefficient if calculated by internal algrithm which is based on Lima,2005
  end type darcy_impes_MIM_options_type
   
  type prog_sfield_MIM_type
     logical :: have_MIM_source
     type(immobile_sfield_type) :: immobile_sfield
     type (scalar_field), pointer :: C_a !the averaged concentration
     type (scalar_field), pointer :: old_C_a
     type (scalar_field), pointer :: Fd !the weight factor for the mobile field, which is ratio of mobile part to average mass
     type (scalar_field), pointer :: Fs !the weight factor for the immobile field, which is ratio of immobile part to average mass
     type(leach_chemical_MIM_type) :: chem
  end type prog_sfield_MIM_type

  type heap_property_type
     type(scalar_field), pointer :: rock_d  !diameter of rock particle
  end type heap_property_type
  
end module darcy_transport_types
  
