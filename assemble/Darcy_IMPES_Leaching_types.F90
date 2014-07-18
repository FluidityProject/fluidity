module darcy_impes_leaching_types

#include "fdebug.h"

  use spud
  use fields
  use state_module
  use fldebug
  use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
  
  
  implicit none
  private

  public :: leach_chemical_di, &
            leach_chemical_prog_sfield_src, &
            leach_chem_src_type

  !!----Leaching chemical options-----------------------------------------------------------------------
  !options under darcy_impes_type                                                       
  type leach_chemical_di                                                                
    logical:: have_leach_chem_model = .false.  !if have leaching chemical model of this material phase
  end type leach_chemical_di
  
  !options under darcy_impes_generic_prog_sfield_type
  type leach_chemical_prog_sfield_src
    logical :: have_sol_src = .false. !if have the chemical source terms from solution phase
    logical :: have_dis_src = .false. !if have the chemical source terms from mineral disslotion
    !the array for the source terms from solution phase
    type(leach_chem_src_type), dimension(:), pointer :: sfield_sol_src =>null()
    !the array for the source terms from mineral dissolution
    type(leach_chem_src_type), dimension(:), pointer :: sfield_dis_src =>null()
  end type leach_chemical_prog_sfield_src
  
  !options for leach_chem_src_type for generic prog leach chemical src
  type leach_chem_src_type
    !The source name of the prog sfield from chemical reactions
    character(len=FIELD_NAME_LEN) :: lc_name
    !the stoichiometric factor of the species in the reaction equation 
    real :: sto_factor
  end type leach_chem_src_type


end module
