module darcy_impes_leaching_types

#include "fdebug.h"

  use spud
  use fields
  use state_module
  use fldebug
  use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
  
  
  implicit none
  private

  public :: leach_chemical_type, &
            leach_chemical_prog_sfield_src, &
            leaching_semi_empirical_model_type, &
            leaching_arrhenius_reaction_type, &
            leach_prog_sfield_heat_transfer_src, &
            leach_chemical_prog_sfield_subcycling, &
            leach_chemical_MIM_type

  !!----Leaching chemical options-----------------------------------------------------------------------
    !options for the bulk fluid dependency of the rate constant
  type bulk_fluid_type
     character(len=FIELD_NAME_LEN) :: lc_name
     real :: order
     integer:: phase
  end type bulk_fluid_type
  
  !options for the rate constant, arrhenius type
  type arrhenius_rate_constant_type
     type(scalar_field), pointer :: A !prefactor of the reaction
     real :: ae !activation energy
     real :: second_ae  !the second ae according to the change of the dissolution mechanics
     real :: gc !gas constant
     type(bulk_fluid_type), dimension(:), allocatable :: bulk
  end type arrhenius_rate_constant_type
  
  type leaching_bio_ferrous_oxi
     logical :: ifbio = .false.
     type(scalar_field_pointer), dimension(:), pointer ::  phi_ore => null()
     type(scalar_field_pointer), dimension(:), pointer ::  phi_l => null()
     type(scalar_field_pointer), dimension(:), pointer ::  miu => null()
     type(scalar_field_pointer), dimension(:), pointer :: phi_l_src => null()
     type(scalar_field) :: phi_l_src_t
     type(scalar_field) :: dcdt
     type(scalar_field), pointer :: cl => null()
     type(scalar_field), pointer :: cfe2 => null()
     real, dimension(:), allocatable  :: phi_max
     real :: Y
     real,dimension(:), allocatable  :: k1
     real,dimension(:), allocatable  :: k2
     real,dimension(:), allocatable  :: k_death
     real,dimension(:), allocatable  :: miu_max
     real,dimension(:), allocatable  :: T_shift
     real :: kmo
     real :: kmfe2
  end type leaching_bio_ferrous_oxi
    
  type leaching_arrhenius_reaction_type
     type(scalar_field), pointer :: dcdt => null() !the change rate of concentration per volume of heap
     type(arrhenius_rate_constant_type):: ak
     type(leaching_bio_ferrous_oxi):: bio
  end type leaching_arrhenius_reaction_type
  
  type leaching_reaction_cap_type
     logical :: have_cap= .false. !the cap used to limit the max value of the concentration, useful for the singular corner of heap
     integer, dimension(:), allocatable :: field_index
     integer, dimension(:), allocatable :: field_phase
     real, dimension(:), allocatable  :: cap_val !the maximun value of the field
  end type leaching_reaction_cap_type
  
  type leaching_semi_empirical_model_type
     type(scalar_field), pointer :: dcdt => null() !the change rate of concentration per volume of heap
     type(scalar_field), pointer :: ex_r => null() !the extraction rate
     type(scalar_field), pointer :: ex => null() !the extraction
     type(scalar_field), pointer :: mc => null() !the molar concentration of the mineral species per volume of heap
     type(arrhenius_rate_constant_type):: ak
     real, dimension(:,:), pointer :: spline_coe !the cubic spline coefficient
     integer :: ndata !the number of the experiment data points used to do spline interpolation
     real, dimension(:), allocatable :: exp_ex, exp_exrk !experiment data
     type(leaching_reaction_cap_type):: cap
     logical :: check_Eh= .false.
     real :: max_Eh !the maximun Eh to limit the dissolution
  end type leaching_semi_empirical_model_type
  
  !jarocite precipitation
  type leaching_internal_algorithm_jaro
     !the name of the scalar field used to calculate pH and Fe3 concentration
     character(len=FIELD_NAME_LEN) :: H_name 
     character(len=FIELD_NAME_LEN) :: Fe3_name
     real :: rate
     type(scalar_field), pointer :: dcdt => null() !the change rate of concentration of Fe3+ per volume of heap
     type(scalar_field), pointer :: js => null() !the molar concentration of jarosite per volume of heap
  end type leaching_internal_algorithm_jaro

  type leaching_internal_algorithm_oxdi
     !the name of the scalar field used for gas phase oxygen and liquid phase oxygen 
     character(len=FIELD_NAME_LEN) :: og_name !gas phase in phase 1
     character(len=FIELD_NAME_LEN) :: o2_name !liquid phase in phase 2
     type(scalar_field), pointer :: dcdt => null() !the change rate of concentration per volume of heap
  end type leaching_internal_algorithm_oxdi

  type leaching_internal_algorithm_sulf
     !the name of the scalar field used to calculate H and liquid phase oxygen concentration
     character(len=FIELD_NAME_LEN) :: H_name 
     character(len=FIELD_NAME_LEN) :: o2_name
     type(scalar_field), pointer :: dcdt => null() !the change rate of concentration per volume of heap
     type(scalar_field), pointer :: S0 => null() !the molar concentration of jarosite per volume of heap
     logical :: bio = .false. !Wether the S0 dissolution is under bacteria activity
     real :: ps  !the persentage of the S0 to dissolve when the leaching is non-bio
  end type leaching_internal_algorithm_sulf
  
  

  type leaching_internal_algorithm_gangue
     !name of H+ used to dissolve gangue mineral
     character(len=FIELD_NAME_LEN) :: H_name
     type(scalar_field), pointer :: dcdt => null() !the change rate of concentration per volume of heap
     real :: u !the reaction rate constant of gangue mineral dissolution, 1/s
  end type leaching_internal_algorithm_gangue
  
  !options for leaching solution phase reaction
  type leaching_solution_phase_type
     type(leaching_internal_algorithm_jaro):: jaro  !jarocite precipitation
     type(leaching_internal_algorithm_oxdi)::oxdi !oxygen dissolution
     type(leaching_arrhenius_reaction_type)::feox !ferrous oxidation
  end type leaching_solution_phase_type
  
  !options for leaching mineral dissolution
  type leaching_mineral_dissolution_type
     type(leaching_internal_algorithm_sulf):: sulf  !elemental sulfur dissolution
     type(leaching_semi_empirical_model_type)::pyri !dissolution of pyrite
     type(leaching_semi_empirical_model_type)::chal !dissolution of chalcopyrite
     type(leaching_internal_algorithm_gangue):: gang !gang mineral dissolution
  end type leaching_mineral_dissolution_type
  
  !options for leach_chem_src_type for generic prog leach chemical src
  type leach_chem_src_type
    !The source name of the prog sfield from chemical reactions
    character(len=FIELD_NAME_LEN) :: lc_name
    !the stoichiometric factor of the species in the reaction equation 
    real :: sto_factor
    type(scalar_field), pointer :: sfield => null()
  end type leach_chem_src_type
  
  type leach_chemical_prog_sfield_src_linearization
     logical:: have

  end type leach_chemical_prog_sfield_src_linearization  

  !options under darcy_impes_generic_prog_sfield_type
  type leach_chemical_prog_sfield_src
    logical :: have_sol_src = .false. !if have the chemical source terms from solution phase
    logical :: have_dis_src = .false. !if have the chemical source terms from mineral dissolution
    logical :: have_chem_src = .false. !if have any type of chemical source terms
    !the array for the source terms from solution phase
    type(leach_chem_src_type), dimension(:), pointer :: sfield_sol_src =>null()
    !the array for the source terms from mineral dissolution
    type(leach_chem_src_type), dimension(:), pointer :: sfield_dis_src =>null()
    type(leach_chemical_prog_sfield_src_linearization) :: src_linear
  end type leach_chemical_prog_sfield_src

  type leach_prog_sfield_heat_transfer_src
    logical :: have_heat_src = .false. !if have the leach temperature heat transfer sources
    type(leach_chemical_prog_sfield_src_linearization) :: src_linear
  end type leach_prog_sfield_heat_transfer_src

  type mineral_dissolution_heat_src_type
    real :: Enthalpy !the enthalpy of the dissolution
    type(scalar_field), pointer :: md_src !scalar field of mineral dissolution rate (mole/m^3_heap)
  end type mineral_dissolution_heat_src_type

  type solution_phase_heat_src_type
     real :: Enthalpy !the enthalpy of the reaction
     type(scalar_field), pointer :: sr_src !scalar field of solution reaction rate (mole/m^3_heap)
  end type solution_phase_heat_src_type
  
  type integer_array_type
    integer,  dimension(:), allocatable :: list
  end type integer_array_type
   
  
  type rock_surface_temperature_src_type
      logical :: have_rtss=.false. !if have rock surface temperature source
      type(integer_array_type),  dimension(:), allocatable :: ids !the boundary ids
      integer,  dimension(:,:), allocatable :: sele !the surface numbering of the boundary id, the first column is the source number, and the second column is the surface number
      integer,  dimension(:), allocatable :: transfer_type !the heat transfer type, 1 is convection, 2 is radiation
      real, dimension(:), allocatable :: coefficient !the heat transfer coefficient or emissivity
      real, dimension(:), allocatable :: ref_t !the reference temperature
      real, dimension(:), allocatable :: A !the surface area 
      type(scalar_field_pointer), dimension(:), pointer :: src=> null() !the surface heat transfer source
   end type rock_surface_temperature_src_type
  
  type leach_heat_transfer_model
     logical:: have_ht = .false.
     logical:: heat_transfer_single = .false.
     logical:: heat_transfer_two = .false.
     logical:: heat_transfer_three = .false.
     logical:: prog_liquid_temperature = .false.
     type(scalar_field), pointer :: liquid_temperature => null()
     type(scalar_field), pointer :: air_temperature => null()
     type(scalar_field), pointer :: rock_temperature => null()
     type(scalar_field), pointer :: rock_temperature_src => null()
     logical:: have_rock_temperature_src = .false.
     type(scalar_field), pointer :: rock_cp => null()
     type(scalar_field), pointer :: rock_density => null()
     type(scalar_field), pointer :: liquid_cp => null()
     type(scalar_field), pointer :: K_eff_ls => null() !liquid-solid effective heat transfer coefficient
     type(scalar_field_pointer), dimension(:), pointer :: two_phase_src_solid => null() !rock temperature heat transfer term
     type(scalar_field_pointer), dimension(:), pointer :: two_phase_src_liquid => null() !liquid temperature heat transfer term
      type(mineral_dissolution_heat_src_type), dimension(:), allocatable :: rock_md_src !mineral dissolution heat transfer sources
     type(solution_phase_heat_src_type), dimension(:), allocatable :: liquid_sr_src !solution phase reactions heat transfer sources
     type(rock_surface_temperature_src_type) :: rtss
  end type leach_heat_transfer_model

  type leach_wetting_efficiency_type
     logical :: have_wet_eff = .false.
     type(scalar_field), pointer :: wet_eff
  end type leach_wetting_efficiency_type

  !options for leaching chemical model under darcy_impes_type                                                       
  type leach_chemical_type                                                                
    logical:: have_leach_chem_model = .false.  !if have leaching chemical model of this material phase
    logical:: have_dis = .false. !if have leaching minderal dissolution
    logical:: have_sol = .false. !if have leaching solution phase reactions
    logical:: if_MIM= .false. !!if solve the leaching chemical model with MIM
    type(leaching_mineral_dissolution_type):: dis 
    type(leaching_solution_phase_type) :: sol 
    type(leach_heat_transfer_model)::ht
    type(leach_wetting_efficiency_type) :: wet_eff
    type(scalar_field), pointer :: Eh
  end type leach_chemical_type
 
  type leach_chemical_prog_sfield_subcycling
    logical:: have_leach_subcycle
    logical:: if_dinamic_dt
    real,pointer :: sub_dt !the dt of subcycling
    real :: beta !the weighting number
    integer::max_numsub !the maximun number of subcycling
    integer:: number_subcycle !the  numbmer of subcycling
    real :: start_extraction !the minimun extraction of chalcopyrite to start the dynaimic dt
    type(scalar_field), dimension(:), allocatable :: sub_lht !the tentative liquid holdup during subcycling
    type(scalar_field), dimension(:), allocatable :: sub_lht_im !the tentative liquid holdup for immobile one
    type(scalar_field), dimension(:), allocatable :: sub_rhs
    type(csr_matrix), dimension(:), allocatable :: sub_adv_diff !the tentative absorption, advection and diffusion during subcycling
    type(scalar_field), dimension(:), allocatable :: old_sub_lht
    type(scalar_field), dimension(:), allocatable :: old_sub_lht_im
    type(scalar_field), dimension(:), allocatable :: old_sub_rhs
    type(csr_matrix), dimension(:), allocatable :: old_sub_adv_diff 
    type(scalar_field), dimension(:), allocatable :: iterated_sfield
    type(scalar_field), dimension(:), allocatable :: iterated_imsfield
    type(scalar_field), pointer :: dy_field !field used to decide the time step of the dynamic subcycling 
    type(scalar_field), pointer :: old_dy_field
    type(scalar_field), pointer :: dynamic_dt !the pointer field to write the dynamic dt into the stat file as a scalar field, due to the fact that output a constant value to stat file need to change the 'diagnoctic_variables' file, therefore to avoid changing that file, we choose output dynamic_dt as a scalar field.
 end type leach_chemical_prog_sfield_subcycling

 type leach_chemical_MIM_field_type
    type (scalar_field), pointer :: sfield
    type (scalar_field) :: p_src !the positive part of the chemical source term
    type (scalar_field) :: n_src !the negative part of the chemical source term    
 end type leach_chemical_MIM_field_type
 

  type  leach_chemical_MIM_type
     logical :: have_chem = .false.
     logical :: if_src_linear= .false.
     type (leach_chemical_MIM_field_type) :: im_src !the chemical reaction src of the immobile part
     type (leach_chemical_MIM_field_type) :: mo_src !the chemical reaction src of the mobile part
  end type leach_chemical_MIM_type  
  
end module
