!This file contains all the subroutines of leaching chemical model

#include "fdebug.h"

module darcy_impes_leaching_chemical_model

  use spud
  use fields
  use state_module
  use fldebug
  use field_options
  use spud
  use data_structures
  use initialise_fields_module
  use global_parameters, only: OPTION_PATH_LEN
  use darcy_impes_leaching_types
  use vtk_cache_module, only:vtk_cache_finalise
  use fefields, only: compute_cv_mass

  use darcy_impes_assemble_type, only: darcy_impes_type
  use darcy_transport_model


  implicit none
  private
  
  public :: initialize_leaching_chemical_model, &
            finalize_leaching_chemical_model, &
            allocate_leaching_chemical_prog_sfield_src, &
            calculate_leaching_chemical_model, &
            allocate_leach_heat_transfer_prog_Temperature_src, &
            calculate_leach_heat_transfer_src, &
            calculate_leach_rock_temperature
         
  
  
  contains
  
  
  subroutine initialize_leaching_chemical_model(di)
     !!initialize the leaching chemical model
      
     type(darcy_impes_type), intent(inout) :: di
     
  
     !local parameter
     type(vector_field), pointer :: position =>null()
     integer :: p,f,flc, ns, nd, nb,fb, stat,ndata,nshape(2),n_md,f_md,i,nc,isf, ff,nphi,ncount
     character(len=OPTION_PATH_LEN) :: option_path, reaction_name, path_l,path_md,md_name,name_temp,cap_name, f_name,o_name,fe2_name
        !---------------------allocate the fields in the chemical model-------------
        !Extract scalar field directly under the model
        if (have_option('/Leaching_chemical_model/scalar_field')) then
           ns= option_count('/Leaching_chemical_model/scalar_field')
           option_path=('/Leaching_chemical_model/scalar_field')
           do f= 1, ns
              call get_option(trim(option_path)//'['//int2str(f-1)//']/name', f_name)
              select case(trim(f_name)) 
                case("Eh")
                  di%lc%Eh => extract_scalar_field(di%state(1), 'Eh', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract Eh in leaching chemical model')
                  end if
                case default
                  FLAbort("scalar field " // trim(f_name) // " leaching chemical model is not found")
              end select
           end do
           
        end if 
         
   
        !for the solution phase reactions
        if (have_option('/Leaching_chemical_model/SolutionPhaseReactions')) then           
           di%lc%have_sol=.true.
           !allocate the solution phase reaction
           ns= option_count('/Leaching_chemical_model/SolutionPhaseReactions/reaction')
           option_path=('/Leaching_chemical_model/SolutionPhaseReactions/reaction')
           do f= 1, ns
           
              call get_option(trim(option_path)//'['//int2str(f-1)//']/name', reaction_name)
              select case(trim(reaction_name)) 
                
                case("Ferrous_Oxidation")
                  di%lc%sol%feox%dcdt => extract_scalar_field(di%state(2), 'feox_dFe2_dt', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named feox_dFe2_dt')
                  end if
                  
                  if (have_option(trim(option_path)//'::Ferrous_Oxidation/Dissolution_Algorithm::bio_leaching')) then
                     
                     path_l=trim(option_path)//'::Ferrous_Oxidation/Dissolution_Algorithm::bio_leaching'
                     di%lc%sol%feox%bio%ifbio = .true.
                     nphi=0
                     do ncount=1, size(di%generic_prog_sfield)
                        if (di%generic_prog_sfield(ncount)%sfield%name(1:5)=='phi_l') then
                          nphi=nphi+1 
                        end if
                     end do 
                      
                     allocate(di%lc%sol%feox%bio%miu(nphi))
                     allocate(di%lc%sol%feox%bio%phi_l(nphi))
                     allocate(di%lc%sol%feox%bio%phi_ore(nphi))
                     allocate(di%lc%sol%feox%bio%phi_l_src(nphi))
                     allocate(di%lc%sol%feox%bio%phi_max(nphi))
                     allocate(di%lc%sol%feox%bio%miu_max(nphi))
                     allocate(di%lc%sol%feox%bio%k1(nphi))
                     allocate(di%lc%sol%feox%bio%k2(nphi))
                     allocate(di%lc%sol%feox%bio%k_death(nphi))
                     allocate(di%lc%sol%feox%bio%T_shift(nphi))

                     call allocate(di%lc%sol%feox%bio%dcdt, di%pressure_mesh)
                     
                     do ncount=1,size(di%lc%sol%feox%bio%phi_l)
                       di%lc%sol%feox%bio%miu(ncount)%ptr => extract_scalar_field(di%state(1), 'miu'//int2str(ncount)//'_feox_', stat=stat)
                       if (.not. stat==0) then
                         FLAbort('failed to extract the scaler field miu under ferrous oxidation by bacteria,check whether the number is correct')
                       end if
                       
                       di%lc%sol%feox%bio%phi_ore(ncount)%ptr => extract_scalar_field(di%state(1), 'phi_ore'//int2str(ncount)//'_feox_', stat=stat)
                       if (.not. stat==0) then
                         FLAbort('failed to extract the scaler field phi_ore under ferrous oxidation by bacteria, check whether the number is correct')
                       end if
                       
                       if (di%MIM_options%have_MIM(2)) then
                       
                                  
                         di%lc%sol%feox%bio%phi_l(ncount)%ptr => extract_scalar_field(di%state(2), 'phi_l'//int2str(ncount)//'Average_mass', stat=stat)
                         if (.not. stat==0) then
                           FLAbort('failed to extract the scaler field Average phi_l for ferrous oxidation by bacteria,check whether the number is correct')
                         end if
                       else   
          
                         di%lc%sol%feox%bio%phi_l(ncount)%ptr => extract_scalar_field(di%state(2), trim('phi_l'//int2str(ncount)), stat=stat)
                         
                         if (.not. stat==0) then
                           FLAbort('failed to extract the scaler field phi_l for ferrous oxidation by bacteria,check whether the number is correct')
                         end if
                       end if
                       
                       di%lc%sol%feox%bio%phi_l_src(ncount)%ptr => extract_scalar_field(di%state(2), 'phi_l'//int2str(ncount)//'_Bacteria_Ferrous_oxidation', stat=stat)
                       if (.not. stat==0) then
                         FLAbort('failed to extract the scaler field  Bacteria_Ferrous_oxidation under phi_l, check whether the number is correct')
                       end if
                       !extract the constant
                       path_l=trim(option_path)//'::Ferrous_Oxidation/Dissolution_Algorithm::bio_leaching'
                       
                       call get_option(trim(path_l)//'/scalar_field::miu'//int2str(ncount)//'/k1', di%lc%sol%feox%bio%k1(ncount))
                       call get_option(trim(path_l)//'/scalar_field::miu'//int2str(ncount)//'/k2', di%lc%sol%feox%bio%k2(ncount))
                       call get_option(trim(path_l)//'/scalar_field::miu'//int2str(ncount)//'/k_death', di%lc%sol%feox%bio%k_death(ncount))
                       call get_option(trim(path_l)//'/scalar_field::miu'//int2str(ncount)//'/phi_max', di%lc%sol%feox%bio%phi_max(ncount))
                       call get_option(trim(path_l)//'/scalar_field::miu'//int2str(ncount)//'/miu_max', di%lc%sol%feox%bio%miu_max(ncount))

                       call get_option(trim(path_l)//'/scalar_field::miu'//int2str(ncount)//'/T_shift', di%lc%sol%feox%bio%T_shift(ncount))
                       
                     end do

                     !get the name of liquid oxygen
                     call get_option(trim(path_l)//'/Oxygen_name', o_name)
                     !get the name of Ferrous
                     call get_option(trim(path_l)//'/Ferrous_name', fe2_name)
                     
                     call get_option(trim(path_l)//'/kmo', di%lc%sol%feox%bio%kmo)
                     call get_option(trim(path_l)//'/kmfe2', di%lc%sol%feox%bio%kmfe2)
                     
                     call get_option(trim(path_l)//'/Y', di%lc%sol%feox%bio%Y)
                     
                     if (di%MIM_options%have_MIM(2)) then
                       
                       di%lc%sol%feox%bio%cl => extract_scalar_field(di%state(2), trim(o_name)//'Average_mass', stat=stat)
                       if (.not. stat==0) then
                         FLAbort('failed to extract the scaler field Average liquid oxygen for ferrous oxidation by bacteria')
                       end if
                       
                       
                       di%lc%sol%feox%bio%cfe2 => extract_scalar_field(di%state(2), trim(fe2_name)//'Average_mass', stat=stat)
                       if (.not. stat==0) then
                         FLAbort('failed to extract the scaler field Average Ferrous concentration for ferrous oxidation by bacteria')
                       end if

                     else
                     
                       di%lc%sol%feox%bio%cl => extract_scalar_field(di%state(2), trim(o_name), stat=stat)
                       if (.not. stat==0) then
                         FLAbort('failed to extract the scaler field  liquid oxygen for ferrous oxidation by bacteria')
                       end if
                       
                       
                       di%lc%sol%feox%bio%cfe2 => extract_scalar_field(di%state(2), trim(fe2_name), stat=stat)
                       if (.not. stat==0) then
                         FLAbort('failed to extract the scaler field  Ferrous concentration for ferrous oxidation by bacteria')
                       end if
                       
                     end if

                  else
                     !get the reaction prefactor
                     di%lc%sol%feox%ak%A => extract_scalar_field(di%state(1), 'feox_prefactor', stat=stat)
                     if (.not. stat==0) then
                        FLAbort('failed to extract the leaching reaction named feox_prefactor')
                     end if

                     !get activation energy
                     path_l = trim(option_path)//'::Ferrous_Oxidation/Dissolution_Algorithm::Non-bio_leaching'
                     call get_option(trim(path_l)//'/rate_constant_Arrhenius/activation_energy',& 
                                                                          di%lc%sol%feox%ak%ae) 
                     !get gas constant
                     call get_option(trim(path_l)//'/rate_constant_Arrhenius/gas_constant',& 
                                                                     di%lc%sol%feox%ak%gc)
                     !count the reacting bulk species
                     nb=option_count(trim(path_l)//'/bulk_fluid_conditions/bulk_concentration')
                     if (nb < 1)  FLAbort('the number of reacting bulk concentration for ferrous oxidation should be at least 1')
                     allocate(di%lc%sol%feox%ak%bulk(nb))
                  
                     do fb=1, nb
                       !get related bulk species name
                       call get_option(trim(path_l)//'/bulk_fluid_conditions/bulk_concentration['//int2str(fb-1)//']/name', di%lc%sol%feox%ak%bulk(fb)%lc_name)
                       !gri the order of reaction
                       call get_option(trim(path_l)//'/bulk_fluid_conditions/bulk_concentration::'//trim(di%lc%sol%feox%ak%bulk(fb)%lc_name)&
                                                                                          //'/order', di%lc%sol%feox%ak%bulk(fb)%order)
                       call get_option(trim(path_l)//'/bulk_fluid_conditions/bulk_concentration::'//trim(di%lc%sol%feox%ak%bulk(fb)%lc_name)&
                                                                                                              //'/phase', di%lc%sol%feox%ak%bulk(fb)%phase)
                     end do
                     
                  end if
                   
                  
                case('Jarosite_Precipitation')
                  di%lc%sol%jaro%dcdt => extract_scalar_field(di%state(2), 'Jarosite_dM_dt', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named Jarosite_dM_dt')
                  end if
                  di%lc%sol%jaro%js => extract_scalar_field(di%state(2), 'Jarosite_molar_concentration', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named Jarosite_molar_concentration')
                  end if
                  path_l = trim(option_path)//'::'//trim(reaction_name)
                  !get the name of the the scalar field used to calculate pH and Fe3 concentration
                  call get_option(trim(path_l)//'/H_name',di%lc%sol%jaro%H_name)
                  call get_option(trim(path_l)//'/Fe3_name',di%lc%sol%jaro%Fe3_name)
                  !get the rate constant to calculate the jarosite precipitation
                  call get_option(trim(path_l)//'/rate_constant',di%lc%sol%jaro%rate)



                case('Oxygen_dissolution')
                  di%lc%sol%oxdi%dcdt => extract_scalar_field(di%state(2), 'oxdi_dOg_dt', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named oxdi_dOg_dt')
                  end if
                  path_l = trim(option_path)//'::'//trim(reaction_name)
                  !get the name of the the scalar field used for gas phase oxygen and liquid phase oxygen
                  call get_option(trim(path_l)//'/og_name',di%lc%sol%oxdi%og_name)
                  call get_option(trim(path_l)//'/o2_name',di%lc%sol%oxdi%o2_name)
                  
                case default
                  FLAbort("Leaching chemical algorithm " // trim(reaction_name) // " not found")
              end select
           
           end do
        
        else
           di%lc%have_sol=.false.

        end if
        
        !for mineral dissolution
        if (have_option('/Leaching_chemical_model/MineralDissolution')) then
           di%lc%have_dis=.true.
           !allocate the mineral dissolution
           nd=option_count('/Leaching_chemical_model/MineralDissolution/reaction')
           option_path=('/Leaching_chemical_model/MineralDissolution/reaction')
           do f=1, nd
             
              call get_option(trim(option_path)//'['//int2str(f-1)//']/name', reaction_name)
              select case(trim(reaction_name))
              
                case('CuFeS2_oxidation_aqueous_ferric_sulfate')
                  if (have_option(trim(trim(option_path)//'::'//trim(reaction_name)//'/Optimum_Eh'))) then
                     di%lc%dis%chal%check_Eh=.true.
                     call get_option(trim(option_path)//'::'//trim(reaction_name)//'/Optimum_Eh/Eh_optimum_value', & 
                             di%lc%dis%chal%max_Eh)
                     call get_option(trim(option_path)//'::'//trim(reaction_name)//'/Optimum_Eh/higher_activation_energy', & 
                          di%lc%dis%chal%ak%second_ae)

                  else
                     di%lc%dis%chal%check_Eh=.false.
                  end if

                  
                  if (have_option(trim(trim(option_path)//'::'//trim(reaction_name)//'/Cap_field'))) then
                     di%lc%dis%chal%cap%have_cap =.true.
                     nc=option_count(trim(trim(option_path)//'::'//trim(reaction_name)//'/Cap_field'))
                     allocate(di%lc%dis%chal%cap%field_index(nc))
                     allocate(di%lc%dis%chal%cap%cap_val(nc))
                     
                     do i =1, nc
                        call get_option(trim(option_path)//'::'//trim(reaction_name)//'/Cap_field/name', & 
                             cap_name)
                        call get_option(trim(option_path)//'::'//trim(reaction_name)//'/Cap_field::'//trim(cap_name)//'/field_maximun_value', di%lc%dis%chal%cap%cap_val(i))
                        
                        do isf=1,size(di%generic_prog_sfield)
                           if (trim(cap_name)==trim(di%generic_prog_sfield(isf)%sfield%name)) then
                              di%lc%dis%chal%cap%field_index(i)=isf                      
                           end if
                           
                        end do
                        
                     end do
                     
                  else
                     di%lc%dis%chal%cap%have_cap =.false.
                  end if 
                  di%lc%dis%chal%dcdt => extract_scalar_field(di%state(2), 'chal_dCuFeS2_dt', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named chal_dCuFeS2_dt')
                  end if
                  !get extraction rate
                  di%lc%dis%chal%ex_r => extract_scalar_field(di%state(2), 'chal_extraction_rate', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named chal_extraction_rate')
                  end if
                  !get current extraction
                  di%lc%dis%chal%ex => extract_scalar_field(di%state(2), 'chal_current_extraction', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named chal_current_extraction')
                  end if
                  !get the molar concentraction of the mineral, mole per volume of heap
                  di%lc%dis%chal%mc => extract_scalar_field(di%state(2), 'chal_molar_concentration', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named chal_molar_concentration')
                  end if                 

                  !get the reaction frefactor
                  di%lc%dis%chal%ak%A => extract_scalar_field(di%state(2), 'chal_prefactor', stat=stat)
                  if (.not. stat==0) then
                     FLAbort('failed to extract the leaching reaction named chal_prefactor')
                  end if
 
                  !get activation energy
                  path_l = trim(option_path)//'::'//trim(reaction_name)
                  call get_option(trim(path_l)//'/rate_constant_Arrhenius/activation_energy',&         
                                                                          di%lc%dis%chal%ak%ae)
                  !get gas constant
                  call get_option(trim(path_l)//'/rate_constant_Arrhenius/gas_constant',&
                                                                     di%lc%dis%chal%ak%gc)
                  !count the reacting bulk species
                  nb=option_count(trim(path_l)//'/bulk_fluid_conditions/bulk_concentration')
                  if (nb < 1)  FLAbort('the number of reacting bulk concentration for chalcopyrite dissolution should be at least 1')
                  allocate(di%lc%dis%chal%ak%bulk(nb))
                  do fb=1, nb
                    !get related bulk species name
                    call get_option(trim(path_l)//'/bulk_fluid_conditions/bulk_concentration['//int2str(fb-1)//']/name', di%lc%dis%chal%ak%bulk(fb)%lc_name)
                    !gri the order of reaction
                    call get_option(trim(path_l)//'/bulk_fluid_conditions/bulk_concentration::'//trim(di%lc%dis%chal%ak%bulk(fb)%lc_name)&
                                                                                          //'/order', di%lc%dis%chal%ak%bulk(fb)%order)
                    call get_option(trim(path_l)//'/bulk_fluid_conditions/bulk_concentration::'//trim(di%lc%dis%chal%ak%bulk(fb)%lc_name)&
                                                                                                              //'/phase', di%lc%dis%chal%ak%bulk(fb)%phase)
                  end do

                  !------allocate and do cubic spline interpolation of the experiment data-----------------------
                  call get_option(trim(path_l)//'/experiment_data/number_of_data_points', di%lc%dis%chal%ndata)
                  ndata=di%lc%dis%chal%ndata
                  allocate(di%lc%dis%chal%spline_coe(4,ndata), di%lc%dis%chal%exp_ex(ndata), di%lc%dis%chal%exp_exrk(ndata))
                  
                  !get the size of the experiment data
                  nshape=option_shape(trim(path_l)//'/experiment_data/extraction')
                  if (nshape(1) /= ndata) then
                  !the size of the data points should be the same with the number of data points
                    FLExit('The dimension of experiment data points of extraction for chalcopyrite dissolution should be the same with number of&
                               data points used to do cubic spline interpolation')
                  end if
                  call get_option(trim(path_l)//'/experiment_data/extraction', di%lc%dis%chal%exp_ex)

                  nshape=option_shape(trim(path_l)//'/experiment_data/empirical_extraction_rate_over_k') 
                  if (nshape(1) /= ndata) then
                  !the size of the data points should be the same with the number of data points
                    FLExit('The dimension of experiment data points of empirical extraction rate over k for chalcopyrite dissolution 
                        &should be the same with number of data points used to do cubic spline interpolation')
                  end if
                  call get_option(trim(path_l)//'/experiment_data/empirical_extraction_rate_over_k',di%lc%dis%chal%exp_exrk)                  

                  !do cubic spline interpolation
                  call cubic_spline_coefficient(di%lc%dis%chal%exp_ex,di%lc%dis%chal%exp_exrk,di%lc%dis%chal%spline_coe)
                
                
                case('FeS2_oxidation_aqueous_ferric_sulfate')
                  di%lc%dis%pyri%check_Eh=.false.
                  if (have_option(trim(trim(option_path)//'::'//trim(reaction_name)//'/Cap_field'))) then
                     di%lc%dis%pyri%cap%have_cap =.true.
                     nc=option_count(trim(trim(option_path)//'::'//trim(reaction_name)//'/Cap_field'))
                     allocate(di%lc%dis%pyri%cap%field_index(nc))
                     allocate(di%lc%dis%pyri%cap%cap_val(nc))
                     
                     do i =1, nc
                        call get_option(trim(option_path)//'::'//trim(reaction_name)//'/Cap_field/name', & 
                             cap_name)
                        call get_option(trim(option_path)//'::'//trim(reaction_name)//'/Cap_field::'//trim(cap_name)//'/field_maximun_value', di%lc%dis%pyri%cap%cap_val(i))
                        
                        do isf=1,size(di%generic_prog_sfield)
                           if (trim(cap_name)==trim(di%generic_prog_sfield(isf)%sfield%name)) then
                              di%lc%dis%pyri%cap%field_index(i)=isf
                           end if
                           
                        end do
                        
                     end do
                     
                  else
                     di%lc%dis%pyri%cap%have_cap =.false.
                  end if 
                   
                  di%lc%dis%pyri%dcdt => extract_scalar_field(di%state(2), 'pyri_dFeS2_dt', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named pyri_dFeS2_dt')
                  end if
                  !get extraction rate
                  di%lc%dis%pyri%ex_r => extract_scalar_field(di%state(2), 'pyri_extraction_rate', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named pyri_extraction_rate')
                  end if
                  !get current extraction
                  di%lc%dis%pyri%ex => extract_scalar_field(di%state(2), 'pyri_current_extraction', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named pyri_current_extraction')
                  end if
                  !get the molar concentraction of the mineral, mole per volume of heap
                  di%lc%dis%pyri%mc => extract_scalar_field(di%state(2), 'pyri_molar_concentration', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named pyri_molar_concentration')
                  end if

                 !get the reaction frefactor
                  di%lc%dis%pyri%ak%A => extract_scalar_field(di%state(2), 'pyri_prefactor', stat=stat)
                  if (.not. stat==0) then
                     FLAbort('failed to extract the leaching reaction named pyri_prefactor')
                  end if

                  !get activation energy
                  path_l = trim(option_path)//'::'//trim(reaction_name)
                  call get_option(trim(path_l)//'/rate_constant_Arrhenius/activation_energy',&
                                                                          di%lc%dis%pyri%ak%ae)
                  !get gas constant
                  call get_option(trim(path_l)//'/rate_constant_Arrhenius/gas_constant',&
                                                                     di%lc%dis%pyri%ak%gc)
                  !count the reacting bulk species
                  nb=option_count(trim(path_l)//'/bulk_fluid_conditions/bulk_concentration')
                  if (nb < 1)  FLAbort('the number of reacting bulk concentration for pyrite dissolution should be at least 1')
                  allocate(di%lc%dis%pyri%ak%bulk(nb))
                  do fb=1, nb
                    !get related bulk species name
                    call get_option(trim(path_l)//'/bulk_fluid_conditions/bulk_concentration['//int2str(fb-1)//']/name', di%lc%dis%pyri%ak%bulk(fb)%lc_name)
                    !gri the order of reaction
                    call get_option(trim(path_l)//'/bulk_fluid_conditions/bulk_concentration::'//trim(di%lc%dis%pyri%ak%bulk(fb)%lc_name)&
                                                                                          //'/order', di%lc%dis%pyri%ak%bulk(fb)%order)
                    call get_option(trim(path_l)//'/bulk_fluid_conditions/bulk_concentration::'//trim(di%lc%dis%pyri%ak%bulk(fb)%lc_name)&
                                                                                                              //'/phase', di%lc%dis%pyri%ak%bulk(fb)%phase)
                  end do
                  
                  !------allocate and do cubic spline interpolation of the experiment data-----------------------
                  call get_option(trim(path_l)//'/experiment_data/number_of_data_points', di%lc%dis%pyri%ndata)
                  ndata=di%lc%dis%pyri%ndata
                  allocate(di%lc%dis%pyri%spline_coe(4,ndata),di%lc%dis%pyri%exp_ex(ndata), di%lc%dis%pyri%exp_exrk(ndata))
    
                  !get the size of the experiment data
                  nshape=option_shape(trim(path_l)//'/experiment_data/extraction')
                  if (nshape(1) /= ndata) then
                  !the size of the data points should be the same with the number of data points
                    FLExit('The dimension of experiment data points of extraction for Pyrite dissolution should be the same with number of&
                               data points used to do cubic spline interpolation')
                  end if
                  call get_option(trim(path_l)//'/experiment_data/extraction',di%lc%dis%pyri%exp_ex)

                  nshape=option_shape(trim(path_l)//'/experiment_data/empirical_extraction_rate_over_k') 
                  if (nshape(1) /= ndata) then
                  !the size of the data points should be the same with the number of data points
                    FLExit('The dimension of experiment data points of empirical extraction rate over k for Pyrite dissolution&
                             should be the same with number of data points used to do cubic spline interpolation')
                  end if
                  call get_option(trim(path_l)//'/experiment_data/empirical_extraction_rate_over_k',di%lc%dis%pyri%exp_exrk)  
                  
                  !do cubic spline interpolation
                  call cubic_spline_coefficient(di%lc%dis%pyri%exp_ex,di%lc%dis%pyri%exp_exrk,di%lc%dis%pyri%spline_coe)

                
                case('S0_dissolution')
                  path_l = trim(option_path)//'::'//trim(reaction_name)
                  !get the name of the the scalar field used for H+ and liquid phase oxygen
                  call get_option(trim(path_l)//'/o2_name',di%lc%dis%sulf%o2_name)
                  call get_option(trim(path_l)//'/H_name',di%lc%dis%sulf%H_name)

                  call get_option(trim(path_l)//'/Dissolution_Algorithm/name',name_temp)
                  !wether the S0 is dissolved with bacteria or not
                  select case(trim(name_temp))
                    case('Non-bio_leaching')
                      di%lc%dis%sulf%bio = .false.
                      !get the percentage of S0 dissolution
                      call get_option(trim(path_l)//'/Dissolution_Algorithm'//'::'//trim(name_temp)//&
                      '/percentage_of_dissolve', di%lc%dis%sulf%ps)
                    case('bio_leaching')
                     di%lc%dis%sulf%bio = .true.
                     !not finish yet

                    case default
                      FLAbort("S0 dissolution algorithm" // trim(reaction_name) // " not found")
                  end select
                   
                  
                  di%lc%dis%sulf%S0 => extract_scalar_field(di%state(2), 'sulf_S0', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction species named S0')
                  end if
                  di%lc%dis%sulf%dcdt => extract_scalar_field(di%state(2), 'sulf_dS0_dt', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named dS0_dt')
                  end if

                case('Gangue_mineral_acid_dissolution')
                   path_l = trim(option_path)//'::'//trim(reaction_name)
                   !get the name of the the scalar field used for H+
                   call get_option(trim(path_l)//'/H_name',di%lc%dis%gang%H_name)
                   call get_option(trim(path_l)//'/rate_constant',di%lc%dis%gang%u) !the reaction rate constant, 1/s
                   
                   di%lc%dis%gang%dcdt => extract_scalar_field(di%state(2), 'gang_dG_dt', stat=stat)
                   if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named dG_dt')
                   end if
                   
                case default
                  FLAbort("Leaching chemical algorithm " // trim(reaction_name) // " not found")

             end select
             
           end do   
        else 
           di%lc%have_dis=.false.      

        end if

        !-------------allocate the temperature-----------------------
        !check the heat transfer model
        if (have_option('/Leaching_chemical_model/heat_transfer_model')) then
          option_path='/Leaching_chemical_model/heat_transfer_model'
          di%lc%ht%have_ht=.true.
          if (have_option(trim(option_path)//'/single_phase_heat_transfer')) then
             di%lc%ht%heat_transfer_single = .true.
             di%lc%ht%heat_transfer_two = .false.
             di%lc%ht%heat_transfer_three = .false.
          elseif (have_option(trim(option_path)//'/two_phases_heat_transfer')) then
             di%lc%ht%heat_transfer_two = .true.
             di%lc%ht%heat_transfer_single = .false.
             di%lc%ht%heat_transfer_three = .false.
          else
             di%lc%ht%heat_transfer_three = .true.
             di%lc%ht%heat_transfer_two = .false.
             di%lc%ht%heat_transfer_single = .false.
          end if

        else
         !default to calculate liquid phase heat transfer only
          di%lc%ht%heat_transfer_single = .true.
          di%lc%ht%heat_transfer_two = .false.
          di%lc%ht%heat_transfer_three = .false.
        end if

        !----------check wether the liquid temperature is the prognostic field or not------
        di%lc%ht%prog_liquid_temperature=(have_option('/material_phase::Phase2/scalar_field::Temperature/prognostic/'))
        ! when the leaching model is turned on, the liquid temperature could be a non-prognostic field only if it is a single phase heat transfer model
        if (.not. di%lc%ht%prog_liquid_temperature) then
          
           if (.not. di%lc%ht%heat_transfer_single) then

              FLAbort('The liquid temperature should be prognostic field if it is not a single phase heat transfer model')
           end if
           
        end if
        
        
        if (di%MIM_options%have_MIM(2) .and. (di%lc%ht%prog_liquid_temperature)) then
           di%lc%ht%liquid_temperature => extract_scalar_field(di%state(2), 'TemperatureAverage_mass', stat=stat)
           if (.not. stat==0) then              
              FLAbort('failed to extract the Average liquid temperature in phase 2 for leaching chemical model')              
           end if
        else           
           di%lc%ht%liquid_temperature => extract_scalar_field(di%state(2), 'Temperature', stat=stat)           
           if (.not. stat==0) then              
              FLAbort('failed to extract the liquid temperature in phase 2 for leaching chemical model')              
           end if
        end if
        
           
     
       if (.not. di%lc%ht%heat_transfer_single) then
       
           ! the rock temperature exist and extract it
           ! also initialize the rock temperature with the initial condition
           di%lc%ht%rock_temperature => extract_scalar_field(di%state(1), 'Rock_Temperature', stat=stat)
           if (.not. stat==0) then
              FLAbort('failed to extract the rock temperature in phase 1 for leaching chemical model')
           end if

           !initialize the rock temperature
           position => extract_vector_field(di%state(1), "Coordinate")
           call zero(di%lc%ht%rock_temperature)
           call initialise_field_over_regions(di%lc%ht%rock_temperature, &
                    trim(di%lc%ht%rock_temperature%option_path)//'/diagnostic/initial_condition', &
               position)

           call vtk_cache_finalise()

           nullify(position)

           !extract the source of rock temperature if exists
           if (have_option('/Leaching_chemical_model/heat_transfer_model/two_phases_heat_transfer/scalar_field::Rock_Temperature/scalar_field::Rock_Temperature_Source')) then
             
             di%lc%ht%have_rock_temperature_src = .true.
             
             di%lc%ht%rock_temperature_src => extract_scalar_field(di%state(1), 'Rock_Temperature_Source', stat=stat)

             if (.not. stat==0) then
               FLAbort('failed to extract the rock temperature source in phase 1 for leaching chemical model')
             end if
           end if         

           !extract rock temperature surface heat transfer source
           if (have_option('/Leaching_chemical_model/heat_transfer_model/two_phases_heat_transfer/scalar_field::Rock_Temperature/scalar_field::Rock_Temperature_Surface_Src1')) then

              call initialize_rock_temperature_surface_source(di)
           end if
           
           !extract the rock density and heat capacity
           di%lc%ht%rock_cp => extract_scalar_field(di%state(1), 'Rock_Cp', stat=stat)
           if (.not. stat==0) then
             FLAbort('failed to extract the Rock heat capacity for leaching chemical model')
           end if
           di%lc%ht%rock_density => extract_scalar_field(di%state(1), 'Rock_density', stat=stat)
           if (.not. stat==0) then
             FLAbort('failed to extract the Rock density for leaching chemical model')
           end if
        end if   

        if (di%lc%ht%heat_transfer_two) then
           !extract the rock-liquid effective heat transfer coefficient
           di%lc%ht%K_eff_ls => extract_scalar_field(di%state(1), 'K_eff_sl', stat=stat)
           if (.not. stat==0) then
             FLAbort('failed to extract rock-liquid effective heat transfer coefficient')
           end if

           !extract the rock heat transfer source
           ns= option_count(trim(option_path)//'/two_phases_heat_transfer/heat_transfer_sources/scalar_field')
           if (ns==0) FLAbort('For multiple phases heat transfer, please turn on the source term of the rock phase temperature')
           path_l= trim(option_path)//'/two_phases_heat_transfer/heat_transfer_sources/scalar_field'
           allocate(di%lc%ht%two_phase_src_solid(ns))
           do f= 1, ns
              call get_option(trim(path_l)//'['//int2str(f-1)//']/name', reaction_name)
              select case(trim(reaction_name))

                case('solid_liquid_heat_transfer_rock_phase') 

                   di%lc%ht%two_phase_src_solid(f)%ptr => extract_scalar_field(di%state(1), 'solid_liquid_heat_transfer_rock_phase', stat=stat)
                   if (.not. stat==0) then
                      FLAbort('failed to extract solid_liquid_heat_transfer_rock_phase under phase 1')
                   end if

                case('mineral_dissolution_heat_sources')
                    di%lc%ht%two_phase_src_solid(f)%ptr => extract_scalar_field(di%state(1), 'mineral_dissolution_heat_sources', stat=stat)
                    if (.not. stat==0) then
                      FLAbort('failed to extract mineral_dissolution_heat_sources under phase 1')
                    end if

                    !count the number of mineral dissolution sources
                    n_md=option_count(trim(path_l)//'::mineral_dissolution_heat_sources/mineral_dissolutions')
                    path_md=trim(path_l)//'::mineral_dissolution_heat_sources/mineral_dissolutions'
                    allocate(di%lc%ht%rock_md_src(n_md))
                    do f_md=1,n_md
                      call get_option(trim(path_md)//'['//int2str(f_md-1)//']/name',md_name)
                      di%lc%ht%rock_md_src(f_md)%md_src => extract_scalar_field(di%state(1), md_name, stat=stat)
                      if (.not. stat==0) then
                        FLAbort('failed to extract field '//md_name//' under phase 1')
                      end if
                      call get_option(trim(path_md)//'::'//trim(md_name)//'/Enthalpy', di%lc%ht%rock_md_src(f_md)%Enthalpy)                      
                    end do

                case default
                   FLAbort("Heat transfer algorithm " // trim(reaction_name) // "under rock temperature not found")

              end select
           end do
           
           !******************initialize the leaching heat transfer source term under liquid temperature
           if (.not. have_option('/material_phase::Phase2/scalar_field::Temperature/prognostic/leaching_temperature_sources')) then
              FLAbort("please turn on the liquid temperature as the prognostic field and its leaching temperature sources terms for two phase heat transfer")
           end if
           
           path_l='/material_phase::Phase2/scalar_field::Temperature/prognostic/leaching_temperature_sources/heat_transfer_sources/scalar_field'
           di%lc%ht%liquid_cp => extract_scalar_field(di%state(2), 'Liquid_Cp',stat=stat)
           if (.not. stat==0) then
              FLAbort('failed to extract liquid heat capacity under phase 2')
           end if

           !extract the heat transfer sources
           ns=option_count(trim(path_l))
           if (ns==0) FLAbort('For multiple phases heat transfer, please turn on the source term of liquid phase temperature')
           allocate(di%lc%ht%two_phase_src_liquid(ns)) 
           do f= 1, ns
             call get_option(trim(path_l)//'['//int2str(f-1)//']/name', reaction_name)
             select case(trim(reaction_name))
                case('solid_liquid_heat_transfer_liquid_phase')
                  di%lc%ht%two_phase_src_liquid(f)%ptr => extract_scalar_field(di%state(2), 'solid_liquid_heat_transfer_liquid_phase', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract solid_liquid_heat_transfer_liquid_phase under phase 2')
                  end if

                case('solution_phase_heat_sources')
                  di%lc%ht%two_phase_src_liquid(f)%ptr => extract_scalar_field(di%state(2), 'solution_phase_heat_sources', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract solution_phase_heat_sources under phase 2')
                  end if

                  !count the number of solution_phase_heat_sources
                  n_md=option_count(trim(path_l)//'::solution_phase_heat_sources/solution_phase_reactions')
                  path_md=trim(path_l)//'::solution_phase_heat_sources/solution_phase_reactions'
                  allocate(di%lc%ht%liquid_sr_src(n_md))
                  do f_md=1,n_md
                    call get_option(trim(path_md)//'['//int2str(f_md-1)//']/name',md_name)
                    di%lc%ht%liquid_sr_src(f_md)%sr_src => extract_scalar_field(di%state(1), md_name, stat=stat)
                    if (.not. stat==0) then
                       FLAbort('failed to extract field'//md_name// 'under phase 1')
                    end if
                    call get_option(trim(path_md)//'::'//trim(md_name)//'/Enthalpy', di%lc%ht%liquid_sr_src(f_md)%Enthalpy)
                  end do

                  case default
                    FLAbort("Heat transfer algorithm " // trim(reaction_name) // "under rock temperature not found")
             end select
           end do

        end if

        if (di%lc%ht%heat_transfer_three) then
          !the three phase heat transfer between air, .iquid and rock
          di%lc%ht%air_temperature => extract_scalar_field(di%state(1), 'Temperature', stat=stat)
          if (.not. stat==0) then
             FLAbort('failed to extract the air temperature in heat transfer model for leaching chemical model')
          end if
        end if

        !----------allocate the generic prognostic leaching source terms-------------       
        !--loop over field
        do f=1, size(di%generic_prog_sfield)
              
          option_path=di%generic_prog_sfield(f)%sfield%option_path
          p=di%generic_prog_sfield(f)%phase
          !check the maximun value cap
          if (have_option(trim(option_path)//'/prognostic/LeachingChemicalSourceTerm/field_maximun_value')) then
             di%generic_prog_sfield(f)%have_cap=.true.   
             call get_option(trim(option_path)//'/prognostic/LeachingChemicalSourceTerm/field_maximun_value', di%generic_prog_sfield(f)%cap_val)
          else
             di%generic_prog_sfield(f)%have_cap=.false.   
             
          end if
          
          !---check for heat transfer source under temperature field
          !check the source linearization
          if (have_option(trim(option_path)//'/prognostic/leaching_temperature_sources/Source_Linearization')) then
             di%generic_prog_sfield(f)%lh_src%src_linear%have=.true.

             !for coupled with MIM model
             if (di%generic_prog_sfield(f)%MIM%chem%have_chem) then
                di%generic_prog_sfield(f)%MIM%chem%if_src_linear= .true.
                call allocate(di%generic_prog_sfield(f)%MIM%chem%im_src%p_src, di%pressure_mesh)
                call allocate(di%generic_prog_sfield(f)%MIM%chem%im_src%n_src, di%pressure_mesh)
                call allocate(di%generic_prog_sfield(f)%MIM%chem%mo_src%p_src, di%pressure_mesh)
                call allocate(di%generic_prog_sfield(f)%MIM%chem%mo_src%n_src, di%pressure_mesh)
             else
                di%generic_prog_sfield(f)%MIM%chem%if_src_linear= .false.
             end if
             
          else
             di%generic_prog_sfield(f)%lh_src%src_linear%have=.false.
          end if

          if (have_option(trim(option_path)//'/prognostic/LeachingChemicalSourceTerm/Source_Linearization')) then
             di%generic_prog_sfield(f)%lc_src%src_linear%have=.true.
             if (di%generic_prog_sfield(f)%MIM%chem%have_chem) then
                di%generic_prog_sfield(f)%MIM%chem%if_src_linear= .true.
                call allocate(di%generic_prog_sfield(f)%MIM%chem%im_src%p_src, di%pressure_mesh)
                call allocate(di%generic_prog_sfield(f)%MIM%chem%im_src%n_src, di%pressure_mesh)
                call allocate(di%generic_prog_sfield(f)%MIM%chem%mo_src%p_src, di%pressure_mesh)
                call allocate(di%generic_prog_sfield(f)%MIM%chem%mo_src%n_src, di%pressure_mesh)
             else
                di%generic_prog_sfield(f)%MIM%chem%if_src_linear= .false.
             end if
          else
            di%generic_prog_sfield(f)%lc_src%src_linear%have=.false.
          end if
                            

          !---check for heat transfer source under temperature field
          ns=option_count(trim(option_path)//'/prognostic/leaching_temperature_sources/heat_transfer_sources/scalar_field')
          if (.not. ns==0) then
            if (di%lc%ht%heat_transfer_single) then
              FLAbort('Cannot use the single heat transfer model to calculate the solid-liquid heat transfer term under the temperature of phase'// int2str(f) )
            end if
      
            di%generic_prog_sfield(f)%lh_src%have_heat_src=.true.
          else
            di%generic_prog_sfield(f)%lh_src%have_heat_src=.false.
          end if
          !---check for solution phase source-----------
          ns = option_count(trim(option_path)//'/prognostic/LeachingChemicalSourceTerm/SolutionPhaseSource/scalar_field')
          if (.not. ns==0) then
            di%generic_prog_sfield(f)%lc_src%have_sol_src = .true.
            allocate(di%generic_prog_sfield(f)%lc_src%sfield_sol_src(ns))

            !extract the name of the chemical reaction and stoichemistry factor
            do flc=1, ns 
              call get_option(trim(option_path)//'/prognostic/LeachingChemicalSourceTerm/&
                             &SolutionPhaseSource/scalar_field['//int2str(flc-1)//']/name', & 
                              di%generic_prog_sfield(f)%lc_src%sfield_sol_src(flc)%lc_name)
              
              call get_option(trim(option_path)//'/prognostic/LeachingChemicalSourceTerm/&
                              &SolutionPhaseSource/scalar_field['//int2str(flc-1)//']/&
                              &diagnostic/stoichiometric_factor', &
                              di%generic_prog_sfield(f)%lc_src%sfield_sol_src(flc)%sto_factor)
              di%generic_prog_sfield(f)%lc_src%sfield_sol_src(flc)%sfield => extract_scalar_field(di%state(p), &
                                                               trim(di%generic_prog_sfield(f)%sfield%name)//'_'&
                                           //trim(di%generic_prog_sfield(f)%lc_src%sfield_sol_src(flc)%lc_name))
            end do
     
          else
            di%generic_prog_sfield(f)%lc_src%have_sol_src = .false.

          end if
          
          !-----check for mineral dissolution phase source------
          nd = option_count(trim(option_path)//'/prognostic/LeachingChemicalSourceTerm/MineralDissolutionSource/scalar_field')
          if (.not. nd==0) then
            di%generic_prog_sfield(f)%lc_src%have_dis_src = .true.
            allocate(di%generic_prog_sfield(f)%lc_src%sfield_dis_src(nd))

            !extract the name of the chemical reaction and stoichemistry factori
            do flc=1, nd
              call get_option(trim(option_path)//'/prognostic/LeachingChemicalSourceTerm/&
                              &MineralDissolutionSource/scalar_field['//int2str(flc-1)//']/name', &
                              di%generic_prog_sfield(f)%lc_src%sfield_dis_src(flc)%lc_name)

              call get_option(trim(option_path)//'/prognostic/LeachingChemicalSourceTerm/&
                              &MineralDissolutionSource/scalar_field['//int2str(flc-1)//']/&
                              &diagnostic/stoichiometric_factor', &
                              di%generic_prog_sfield(f)%lc_src%sfield_dis_src(flc)%sto_factor)
              di%generic_prog_sfield(f)%lc_src%sfield_dis_src(flc)%sfield => extract_scalar_field(di%state(p), &
                                                               trim(di%generic_prog_sfield(f)%sfield%name)//"_"&
                                           //trim(di%generic_prog_sfield(f)%lc_src%sfield_dis_src(flc)%lc_name))
            end do


          else
            di%generic_prog_sfield(f)%lc_src%have_dis_src = .false.
          
          endif

          if (di%generic_prog_sfield(f)%lc_src%have_dis_src .or. di%generic_prog_sfield(f)%lc_src%have_sol_src) then

             di%generic_prog_sfield(f)%lc_src%have_chem_src = .true.
          end if
          
        end do

        !-------------------solid liquid wetting efficiency------------------------------
        if (have_option(trim('/Leaching_chemical_model/liquid_solid_wetting_efficiency'))) then
           di%lc%wet_eff%have_wet_eff =.true.
           di%lc%wet_eff%wet_eff => extract_scalar_field(di%state(1), trim('Wetting_efficiency'))
        end if

        
  end subroutine initialize_leaching_chemical_model

  !--------------------------------------------------------------------------------------------------------
  subroutine finalize_leaching_chemical_model(di)  
     
     !finalize terms from leaching_chemical_model

     type(darcy_impes_type), intent(inout) :: di
     
     !local variables
     integer :: f, ns, nd, flc
     
     character(len=OPTION_PATH_LEN) :: option_path, reaction_name
     
     di%lc%have_leach_chem_model= .False.

     if (di%lc%wet_eff%have_wet_eff) then
        di%lc%wet_eff%have_wet_eff =.false.
        nullify(di%lc%wet_eff%wet_eff)
     end if
     
     
     !deallocate leaching chemical model
     if (di%lc%have_sol) then
       
       di%lc%have_sol=.false.
       ns= option_count('/Leaching_chemical_model/SolutionPhaseReactions/reaction')
       option_path=('/Leaching_chemical_model/SolutionPhaseReactions/reaction')
       do f= 1, ns

         call get_option(trim(option_path)//'['//int2str(f-1)//']/name', reaction_name)
         select case(trim(reaction_name))

           case("Ferrous_Oxidation")
              nullify(di%lc%sol%feox%dcdt)
              if (di%lc%sol%feox%bio%ifbio) then
                 di%lc%sol%feox%bio%ifbio=.false.
                 do flc=1, size(di%lc%sol%feox%bio%phi_l)
                    nullify(di%lc%sol%feox%bio%phi_ore(flc)%ptr)
                    nullify(di%lc%sol%feox%bio%phi_l(flc)%ptr)            
                    nullify(di%lc%sol%feox%bio%miu(flc)%ptr)
                    nullify(di%lc%sol%feox%bio%phi_l_src(flc)%ptr)
                 end do
                 nullify(di%lc%sol%feox%bio%cl)
                 nullify(di%lc%sol%feox%bio%cfe2)
                 deallocate(di%lc%sol%feox%bio%miu)
                 deallocate(di%lc%sol%feox%bio%phi_l)
                 deallocate(di%lc%sol%feox%bio%phi_ore)
                 deallocate(di%lc%sol%feox%bio%phi_l_src)
                 deallocate(di%lc%sol%feox%bio%phi_max)
                 deallocate(di%lc%sol%feox%bio%miu_max)
                 deallocate(di%lc%sol%feox%bio%T_shift)
                 deallocate(di%lc%sol%feox%bio%k1)
                 deallocate(di%lc%sol%feox%bio%k2) 
                 call deallocate(di%lc%sol%feox%bio%dcdt)      
              else
                 nullify(di%lc%sol%feox%ak%A)
                 deallocate(di%lc%sol%feox%ak%bulk)
              end if
              
           case('Jarosite_Precipitation')
              nullify(di%lc%sol%jaro%dcdt)

           case('Oxygen_dissolution')
              nullify(di%lc%sol%oxdi%dcdt)
              
         end select
         
       end do

     end if
     
     
     if (di%lc%have_dis) then
       
       di%lc%have_dis=.false.
       nd=option_count('/Leaching_chemical_model/MineralDissolution/reaction')
       option_path=('/Leaching_chemical_model/MineralDissolution/reaction')
       do f= 1, nd

         call get_option(trim(option_path)//'['//int2str(f-1)//']/name', reaction_name)
         select case(trim(reaction_name))

           case('CuFeS2_oxidation_aqueous_ferric_sulfate')
              nullify(di%lc%dis%chal%dcdt)
              nullify(di%lc%dis%chal%ex_r)
              nullify(di%lc%dis%chal%ex)
              nullify(di%lc%dis%chal%mc)
              nullify(di%lc%dis%chal%ak%A)
              deallocate(di%lc%dis%chal%ak%bulk)
              deallocate(di%lc%dis%chal%spline_coe)
              deallocate(di%lc%dis%chal%exp_ex)
              deallocate(di%lc%dis%chal%exp_exrk)
              if (di%lc%dis%chal%cap%have_cap) then
                 di%lc%dis%chal%cap%have_cap=.false.
                 deallocate(di%lc%dis%chal%cap%field_index)
                 deallocate(di%lc%dis%chal%cap%cap_val)
              end if
              di%lc%dis%chal%check_Eh=.false.
              nullify(di%lc%Eh)
              
           case('FeS2_oxidation_aqueous_ferric_sulfate')
              nullify(di%lc%dis%pyri%dcdt)
              nullify(di%lc%dis%pyri%ex_r)
              nullify(di%lc%dis%pyri%ex)
              nullify(di%lc%dis%pyri%mc)
              nullify(di%lc%dis%pyri%ak%A)
              deallocate(di%lc%dis%pyri%ak%bulk)
              deallocate(di%lc%dis%pyri%spline_coe)
              deallocate(di%lc%dis%pyri%exp_ex)
              deallocate(di%lc%dis%pyri%exp_exrk)
              if (di%lc%dis%pyri%cap%have_cap) then
                 di%lc%dis%pyri%cap%have_cap=.false.
                 deallocate(di%lc%dis%pyri%cap%field_index)
                 deallocate(di%lc%dis%pyri%cap%cap_val)
              end if

           case('S0_dissolution')
              nullify(di%lc%dis%sulf%S0)
              nullify(di%lc%dis%sulf%dcdt)
              di%lc%dis%sulf%bio = .false.

           case('Gangue_mineral_acid_dissolution')
              nullify(di%lc%dis%gang%dcdt)
              
         end select
           
       end do

     end if

     !deallocate leaching prognostic source field 
     do f=1, size(di%generic_prog_sfield)
        if (di%generic_prog_sfield(f)%have_cap) then
           di%generic_prog_sfield(f)%have_cap=.false.
        end if
        
        if (di%generic_prog_sfield(f)%MIM%chem%if_src_linear) then
           di%generic_prog_sfield(f)%MIM%chem%if_src_linear= .false.
           call deallocate(di%generic_prog_sfield(f)%MIM%chem%im_src%p_src)
           call deallocate(di%generic_prog_sfield(f)%MIM%chem%im_src%n_src)
           call deallocate(di%generic_prog_sfield(f)%MIM%chem%mo_src%p_src)
           call deallocate(di%generic_prog_sfield(f)%MIM%chem%mo_src%n_src)
        end if
        
       di%generic_prog_sfield(f)%lh_src%src_linear%have=.false.
       di%generic_prog_sfield(f)%lc_src%src_linear%have=.false.
       if (di%generic_prog_sfield(f)%lh_src%have_heat_src) &
         di%generic_prog_sfield(f)%lh_src%have_heat_src=.false.

       if (di%generic_prog_sfield(f)%lc_src%have_sol_src) then
         di%generic_prog_sfield(f)%lc_src%have_sol_src = .false.
         do flc=1, size(di%generic_prog_sfield(f)%lc_src%sfield_sol_src)
           nullify(di%generic_prog_sfield(f)%lc_src%sfield_sol_src(flc)%sfield)
         end do
         deallocate(di%generic_prog_sfield(f)%lc_src%sfield_sol_src)
       end if
       
       if (di%generic_prog_sfield(f)%lc_src%have_dis_src) then
         di%generic_prog_sfield(f)%lc_src%have_dis_src = .false.
         do flc=1, size(di%generic_prog_sfield(f)%lc_src%sfield_dis_src)
           nullify(di%generic_prog_sfield(f)%lc_src%sfield_dis_src(flc)%sfield)
         end do
         deallocate(di%generic_prog_sfield(f)%lc_src%sfield_dis_src)
       end if

     end do
   
     !deallocate heat transfer model
     if (have_option('/Leaching_chemical_model/heat_transfer_model')) then
       di%lc%ht%have_ht=.false.
       nullify(di%lc%ht%liquid_temperature)
       if (.not. di%lc%ht%heat_transfer_single) then
          nullify(di%lc%ht%rock_temperature)
          nullify(di%lc%ht%rock_cp)
          nullify(di%lc%ht%rock_density)

          if (di%lc%ht%heat_transfer_two) then
             !for rock phase
             nullify(di%lc%ht%K_eff_ls)
             do flc=1,size(di%lc%ht%two_phase_src_solid)
               nullify(di%lc%ht%two_phase_src_solid(flc)%ptr)
             end do
             deallocate(di%lc%ht%two_phase_src_solid)
             
             if (allocated(di%lc%ht%rock_md_src)) then
               do flc=1,size(di%lc%ht%rock_md_src)
                 nullify(di%lc%ht%rock_md_src(flc)%md_src)
               end do
               deallocate(di%lc%ht%rock_md_src)
             end if

             !for liquid phase
             nullify(di%lc%ht%liquid_cp)
             do flc=1,size(di%lc%ht%two_phase_src_liquid)
               nullify(di%lc%ht%two_phase_src_liquid(flc)%ptr)
             end do
             deallocate(di%lc%ht%two_phase_src_liquid)

             if (allocated(di%lc%ht%liquid_sr_src)) then
               do flc=1,size(di%lc%ht%liquid_sr_src)
                 nullify(di%lc%ht%liquid_sr_src(flc)%sr_src)
               end do
               deallocate(di%lc%ht%liquid_sr_src)
             end if  
          end if

          if (di%lc%ht%heat_transfer_three) then
            nullify(di%lc%ht%rock_temperature)
            nullify(di%lc%ht%air_temperature)
          end if    
       end if

       di%lc%ht%heat_transfer_single = .false.
       di%lc%ht%heat_transfer_two = .false.
       di%lc%ht%heat_transfer_three = .false.

       if (di%lc%ht%have_rock_temperature_src) then
          di%lc%ht%have_rock_temperature_src=.false.
          nullify(di%lc%ht%rock_temperature_src)
       end if

       if (di%lc%ht%rtss%have_rtss) then
          call finalize_rock_temperature_surface_source(di)
          
       end if
       
             
    end if

   end subroutine finalize_leaching_chemical_model


   !********The following are the subroutines to calculate fields for the chemical model****************
   
   !-------------calculate leaching chemical model----------------------------------------------!
   subroutine calculate_leaching_chemical_model(di)
      type(darcy_impes_type), intent(inout) :: di
      
      integer :: i
      real :: T,ft, dt
      real, dimension(:), allocatable :: A !pre-factor of the arrhenius rate constant
      type(scalar_field), pointer :: rock_temperature
      
      if (di%lc%have_dis) then  
         if (.not. di%lc%ht%heat_transfer_single) then
            rock_temperature =>di%lc%ht%rock_temperature
         else
            rock_temperature =>di%lc%ht%liquid_temperature
         end if
         
         if (di%lcsub%have_leach_subcycle) then
           dt=di%lcsub%sub_dt
         else
           dt=di%dt
        end if

        !--------------for mineral dissolution-----------------
        ! Chalcopyrite oxidation 
         if (associated(di%lc%dis%chal%dcdt)) call calculate_mineral_dissolution_semi_empirical_model(di,di%state,&
                                          rock_temperature,di%number_pmesh_node,dt,di%lc%dis%chal,&
                                          di%saturation(2)%ptr)
         !pyrite dissolution
         if (associated(di%lc%dis%pyri%dcdt)) call calculate_mineral_dissolution_semi_empirical_model(di,di%state,&
                                          rock_temperature,di%number_pmesh_node,dt,di%lc%dis%pyri,&
                                          di%saturation(2)%ptr)
         !S0 dissolution
         if (associated(di%lc%dis%sulf%dcdt)) call calculate_S0_dissolution(di)
         !gangue mineral dissolution
         if (associated(di%lc%dis%gang%dcdt)) call calculate_gangue_mineral_dissolution(di)
         !finalize                                 
         nullify(rock_temperature)                                 
      end if
      !----------------for solution phase reactions-------------------
      if (di%lc%have_sol) then
         !ferrous oxidation
         if (associated(di%lc%sol%feox%dcdt)) then
           if (di%lc%sol%feox%bio%ifbio) then
            call calculate_ferrous_oxidation_bio_terms(di)
            call set(di%lc%sol%feox%dcdt,di%lc%sol%feox%bio%dcdt)
            call scale(di%lc%sol%feox%dcdt, 1.0/0.056)

          else
            !calcultae the prefactor A, since the partial pressure of oxygen in the reaction rate equation are
            !transformed to the molar concentraction of dissolved oxygen by empirical equilibrium equatiion (Tromans1998)
            allocate(A(di%number_pmesh_node))
            
            node_loop: do i=1,di%number_pmesh_node
               !reaction only start when liquid saturation is not zero
               if (di%saturation(2)%ptr%val(i)>=1.0e-4) then
                  !calculate the equilibrium constant ft
                  !the rate of reaction is in the form of
                  !A1*[ft^-1][DO][Ferrous^2][H^-0.25]exp(-Ea/(RT)) in mol/m^3
                  !A1=A_org/1000.0, where A_org is the oringical prefactor in
                  ! 'A_org*Po*[Ferrous^2][H^-0.25]exp(-Ea/(RT))'
                 T=di%lc%ht%liquid_temperature%val(i)
                 ft=(0.046*(T**2.0)+203.35*T*DLOG(T/298.0)-(299.378+0.092*T)*(T-298)-20.591*(10.0**3.0))/(T*8.3144) 
                 ft=EXP(ft)
                 A(i)=di%lc%sol%feox%ak%A%val(i)/ft
                 
               else
                 A(i)=0.0
               end if
            end do node_loop
            call calculate_solution_phase_arrhenius_type_reaction_rate(di%state,&
                 A,di%lc%ht%liquid_temperature,di%number_pmesh_node,di%lc%sol%feox)
            
            !change the mole/(m^3 solution)/s to mole/(m^3 heap)/s     
            call scale(di%lc%sol%feox%dcdt,di%porosity_pmesh)
            
            call scale(di%lc%sol%feox%dcdt,di%saturation(2)%ptr)
            
            deallocate(A)
            
           end if

         end if

         !jarosite oxidation
         if (associated(di%lc%sol%jaro%dcdt)) then
           call calculate_jarosite_precipitation(di) 
         end if

         !oxygen dissolution
         if (associated(di%lc%sol%oxdi%dcdt)) then
           call calculate_oxygen_dissolution(di)
         end if

        if (di%lc%dis%chal%check_Eh) then
           call calculate_Eh(di)
         end if
         
      end if
     
      contains 

        subroutine calculate_mineral_dissolution_semi_empirical_model(di,states,temperature,node_number,dt,mineral,sat)
              type(darcy_impes_type), intent(inout) :: di
              type(leaching_semi_empirical_model_type), intent(inout) :: mineral
              type(state_type),dimension(:), intent(in) :: states
              integer, intent(in) :: node_number
              type(scalar_field), intent(in) :: temperature
              type(scalar_field), intent(in) :: sat !liquid saturation
              real, intent(in) :: dt !the time step
            
              
              real :: cb_n,mc, k_rate, ext, ext_rk, ext_r, dcdt !bode value of bulk concentration, molar concentration of the mineral, 
                                                    !tate constant,current extraction,the extraction rate with k, the extraction rate, the concentration change rate
              type(scalar_field_pointer), dimension(:), allocatable :: cb !the reacting bulk species
              real, dimension(:),allocatable ::  m ! the single node val of reacting bulk species and order of reaction
              real, dimension(:), pointer :: a,b,c,d !spline coefficienti
              real, dimension(4) :: a_c !the constants used to calculate arrhenius reattion rate
                                        !pre-factor,activation energy,gas constant,temperature
              character(len=FIELD_NAME_LEN) :: cb_name !the name of the reacting species
              integer :: p, nspecies, isp, node, stat,icap

              real :: wet_eff
              
              !get the number of species which take part in reaction
              nspecies = size(mineral%ak%bulk)             
              allocate(cb(nspecies), m(nspecies))
              
              !spline coefficient
              a => mineral%spline_coe(1,:)
              b => mineral%spline_coe(2,:)
              c => mineral%spline_coe(3,:)
              d => mineral%spline_coe(4,:)

              do isp = 1, nspecies
                 m(isp) = mineral%ak%bulk(isp)%order
                 
                 p = mineral%ak%bulk(isp)%phase
                 
                 if (di%MIM_options%have_MIM(p)) then
                    cb_name = trim(mineral%ak%bulk(isp)%lc_name)//'Average_mass'
                 else
                    cb_name = mineral%ak%bulk(isp)%lc_name
                 end if
                 
                 cb(isp)%ptr => extract_scalar_field(states(p), trim(cb_name), stat)
                 if (.not. stat==0) then
                    FLAbort('failed to extract the reacting species')
                 end if
              end do
              
              !Calculate exctraction rate, extraction, and concentration change rate for each node
              !node loop
              node_loop: do node=1, node_number
                   !the current extraction of the element node
                   ext = node_val(mineral%ex,node)
                
                   if (sat%val(node)<=1.0e-8) then
                     call set(mineral%ex_r, node, 0.0)
                     call set(mineral%ex, node, ext)
                     call set(mineral%dcdt, node, 0.0)
                     cycle node_loop
                   end if

                   !check the cap of the reaction
                   !if the cap field reach the maximun value, stop reaction of that node
                   if (mineral%cap%have_cap) then
                      
                      do icap=1, size(mineral%cap%field_index)
                         if (di%generic_prog_sfield(mineral%cap%field_index(icap))%sfield%val(node)>=mineral%cap%cap_val(icap)) then
                            call set(mineral%ex_r, node, 0.0)
                            call set(mineral%ex, node, ext)
                            call set(mineral%dcdt, node, 0.0)
                            cycle node_loop
                         end if                       
                      end do                      
                   end if
                   

                   !the pre_factor is ingnored in the semi-imperical model, set to 1.0
                   a_c(1)=node_val(mineral%ak%A, node)
                   !the activation energy
                   if (mineral%check_Eh) then
                      if (di%lc%Eh%val(node) >=0.65) then
                         a_c(2)=mineral%ak%second_ae
                      
                      else
                         a_c(2)=mineral%ak%ae
                         
                      end if
                   else
                      a_c(2)=mineral%ak%ae
                   end if
                   !the gas constant
                   a_c(3)=mineral%ak%gc
                   !the reaction temperature is based on the liquid temperature
                   a_c(4)=node_val(temperature, node)

                   !the molar concentration of the mineral at the unit of mole per volumn of heap
                   mc= node_val(mineral%mc, node)

                   !if extraction of the mineral is nearly one, not reacted
                   if ((1.0-ext)<=1.0D-15) then
                      k_rate=0.0  !stop reaction
                      go to 1
                   end if
                   
                   call calculate_arrhenius_reaction_rate_constant(nspecies,node,cb,m,a_c,k_rate)
                   
                 1 if (k_rate==0.0) then
                       call set(mineral%ex_r, node, 0.0)
                       call set(mineral%ex, node, ext)
                       call set(mineral%dcdt, node, 0.0)
                   else                           
                       !do cubic spline interpolation for current extraction rate with k
                       call cubic_spline_interpolation(a,b,c,d,mineral%exp_ex,ext,ext_rk)

                       !calculate wetting efficiency
                       if (di%lc%wet_eff%have_wet_eff) then
                          call calculate_solid_liquid_wetting_efficiency(di, node, wet_eff)
                          k_rate = wet_eff*k_rate
                       end if
                       
                       ext_r = ext_rk * k_rate
                   
                       !update the new extraction rate, unit is extraction per day
                       call set(mineral%ex_r, node, ext_r)
                   
                       !update the new extraction of the node
                       ext=ext+(ext_r/86400.0)*dt
                       call set(mineral%ex, node, ext)

                       !calculate the concentration change rate of the mineral
                       !the molar concentration %mc is at the unit of mole/m^3_heap
                       !the concentration change rate is at the unit of mole/m^3_heap/s
                       dcdt= -(ext_r/86400.0)*mc
                       call set(mineral%dcdt, node, dcdt)
                   end if 
                    
              end do node_loop

              !finalize
              !nullify before deallocate
              nullify(a,b,c,d)
              do isp = 1, nspecies
                nullify(cb(isp)%ptr)
              end do

              deallocate(cb)
              deallocate(m)
          
        end subroutine calculate_mineral_dissolution_semi_empirical_model

        subroutine calculate_solution_phase_arrhenius_type_reaction_rate(states,A,temperature,node_number,reaction)
              type(leaching_arrhenius_reaction_type), intent(inout) :: reaction
              type(state_type),dimension(:), intent(in) :: states
              integer, intent(in) :: node_number
              real, dimension(:),intent(in) :: A !pre-factor of the arrhenius rate constant
              type(scalar_field), intent(in) :: temperature
                           
              real, dimension(4) :: a_c !the constants used to calculate arrhenius reattion rate
                                        !pre-factor,activation energy,gas constant,temperature
              character(len=FIELD_NAME_LEN) :: cb_name !the name of the reacting species
              type(scalar_field_pointer), dimension(:), allocatable :: cb !the reacting bulk species
              real, dimension(:),allocatable ::  m ! the single node val of reacting bulk species and order of reaction
              real :: k_rate
              
              integer :: p, nspecies, isp, node, stat
                

              !get the number of species which take part in reaction
              nspecies = size(reaction%ak%bulk)
              allocate(cb(nspecies), m(nspecies))
              
              do isp = 1, nspecies
                 m(isp) = reaction%ak%bulk(isp)%order
                 
                 p = reaction%ak%bulk(isp)%phase

                 if (di%MIM_options%have_MIM(p)) then
                    cb_name = trim(reaction%ak%bulk(isp)%lc_name)//'Average_mass'
                 else
                    cb_name = reaction%ak%bulk(isp)%lc_name
                 end if
                 
                 cb(isp)%ptr => extract_scalar_field(states(p), trim(cb_name), stat)

              
                 if (.not. stat==0) then
                    FLAbort('failed to extract the reacting species')
                 end if
              end do

              !Calculate exctraction rate, extraction, and concentration change rate for each node
              !node loop
              node_loop: do node=1, node_number
                   
                   if (abs(A(node))<=1.0e-16) then
                      
                      call set(reaction%dcdt, node, 0.0)
                      
                   else
                     !the pre_factor is ingnored in the semi-imperical model, set to 1.0
                     a_c(1)=A(node)
                     !the activation energy
                     a_c(2)=reaction%ak%ae
                     !the gas constant
                     a_c(3)=reaction%ak%gc
                     !the reaction temperature is based on the liquid temperature
                     a_c(4)=node_val(temperature, node)
                     
                     call calculate_arrhenius_reaction_rate_constant(nspecies,node,cb,m,a_c,k_rate)

                     !mole per volumn of solution per second
                     call set(reaction%dcdt, node, k_rate)
                   end if
              end do node_loop

              !finalize
              do isp = 1, nspecies
               nullify(cb(isp)%ptr)
              end do

              deallocate(cb)
              deallocate(m)
 
        end subroutine calculate_solution_phase_arrhenius_type_reaction_rate

        subroutine calculate_arrhenius_reaction_rate_constant(nspecies,node,cb,m,ac,k_rate)
           integer, intent(in) :: nspecies, node
           type(scalar_field_pointer), dimension(:), intent(in) :: cb !the reacting bulk species
           real, dimension(4), intent(in) :: ac 
           real, dimension(:),intent(in) :: m 
           real, intent(inout) :: k_rate
           
           integer :: isp
           real ::cb_n
           
           !calculate rate constant, which is k=A*e^(Ea/(R*T))*(ab1**m1)*(ab2**m2)....
           k_rate=ac(1)*(EXP(ac(2)/(ac(3)*ac(4))))
           do isp= 1, nspecies                 
              cb_n = node_val(cb(isp)%ptr, node)

              if (cb_n <=0.00000001) then
                 if (m(isp) > 0.0) then
                   k_rate=0.0  !the species with positive order is the reactant, stop reaction

                   return
                else
                   cb_n  = 1.0 !the species with negative order is the product
                                !let it equal to 1 and make the reaction independent of it
                                !this might not be true, but zero product concentration 
                end if
              end if                                    
              k_rate = k_rate*(cb_n**m(isp)) 
           end do  
                              
        end subroutine calculate_arrhenius_reaction_rate_constant 

   end subroutine calculate_leaching_chemical_model
   
   subroutine calculate_S0_dissolution(di)
     type(darcy_impes_type), intent(inout) :: di
     type(scalar_field), pointer :: o2
     real ::dS0dt, dS0 
     integer :: i, stat
     character(len=OPTION_PATH_LEN) :: Oname

      if (di%MIM_options%have_MIM(2)) then
         Oname=trim(di%lc%dis%sulf%o2_name)//'Average_mass'

      else
         Oname=trim(di%lc%dis%sulf%o2_name)
         
      end if
    
     !only dissolve when there are enough dissolved oxygen
     o2 => extract_scalar_field(di%state(2),Oname, stat=stat)
        
     if (.not. stat==0) then
         FLAbort('failed to extract the scalar field of liquid  phase oxygen to calculate S0 dissolution')
     end if

     if (di%lc%dis%sulf%bio) then
       !not finished yet
     else
       node_loop: do i=1,di%number_pmesh_node
         !the dissolved S0 from chalcopyrite dissolution
         dS0dt = 2.0*di%lc%dis%chal%dcdt%val(i)*di%lc%dis%sulf%ps
         !only dissolve when there are enough dissolved oxygen
         if ((o2%val(i)+dS0dt)>=0.0) then
           !calculate the S0 dissolution, dS0/dt  
           !Assumed percentage of S0 generated by Chalcopyrite dissolution is dissolved to SO4
           ! the new calculated dCuFeS2/dt in mole/m^3_heap/s * the percentage of dissolution
           di%lc%dis%sulf%dcdt%val(i) = dS0dt

           !calculate the current S0
           dS0 = -2.0*di%lc%dis%chal%dcdt%val(i)*(1.0-di%lc%dis%sulf%ps)
           di%lc%dis%sulf%S0%val(i) = di%lc%dis%sulf%S0%val(i)+dS0
           
         else
           dS0 = -2.0*di%lc%dis%chal%dcdt%val(i)
           di%lc%dis%sulf%S0%val(i) = di%lc%dis%sulf%S0%val(i)+dS0
           di%lc%dis%sulf%dcdt%val(i) = 0.0
         end if

       end do node_loop
     end if

     nullify(o2)
      
   end subroutine calculate_S0_dissolution

   subroutine calculate_gangue_mineral_dissolution(di)
     type(darcy_impes_type), intent(inout) :: di

     type(scalar_field), pointer :: H
     integer:: stat,i
     character(len=OPTION_PATH_LEN) :: Hname

     if (di%MIM_options%have_MIM(2)) then
        Hname=trim(di%lc%dis%gang%H_name)//'Average_mass'
     else
        Hname=trim(di%lc%dis%gang%H_name)
     end if

     H => extract_scalar_field(di%state(2),Hname, stat=stat)
     if (.not. stat==0) then
        FLAbort('failed to extract the scalar field of H+ to calculate gangue dissolution')
     end if

     node_loop: do i=1,di%number_pmesh_node
        !in mole/m3_heap/s
        if (H%val(i)<0.1) then
           di%lc%dis%gang%dcdt%val(i)=0.0
        else
           
           di%lc%dis%gang%dcdt%val(i)=-H%val(i)*di%lc%dis%gang%u*di%porosity_pmesh%val(i)*di%saturation(2)%ptr%val(i)
        end if
        
     end do node_loop

     nullify(H)
     
   end subroutine calculate_gangue_mineral_dissolution
   
   subroutine calculate_Eh(di)
      type(darcy_impes_type), intent(inout) :: di
      type(scalar_field), pointer :: Fe2, Fe3
      character(len=OPTION_PATH_LEN) :: Fe2name, Fe3name 
      integer :: i, stat
      
      if (di%MIM_options%have_MIM(2)) then
         Fe2name=trim('Fe2')//'Average_mass'
         Fe3name=trim('Fe3')//'Average_mass'
      else
         Fe2name=trim('Fe2')
         Fe3name=trim('Fe3')
      end if
     
      Fe2 => extract_scalar_field(di%state(2),Fe2name, stat=stat) 
      if (.not. stat==0) then
        FLAbort('failed to extract the scalar field of Fe2 to calculate Eh')
      end if
      
      Fe3 => extract_scalar_field(di%state(2),Fe3name, stat=stat)
      if (.not. stat==0) then
        FLAbort('failed to extract the scalar field of Fe3 to calculate Eh')
      end if
      
      do i=1,di%number_pmesh_node         
         !if reactant Fe3 is near zero, stop reaction

         if ((Fe2%val(i)<=0.00001) .or. (Fe3%val(i)<=0.00001)) then
           di%lc%Eh%val(i)=0.0
         else
            di%lc%Eh%val(i)=0.67+0.059*DLOG10(Fe3%val(i)/Fe2%val(i))
         end if
      end do

      nullify(Fe2)
      nullify(Fe3)
      
   end subroutine calculate_Eh

   subroutine calculate_jarosite_precipitation(di)
      type(darcy_impes_type), intent(inout) :: di
      
      type(scalar_field), pointer :: H, Fe3 
      real :: pH, F, rhs, dpre 
      integer :: stat, i      
      character(len=OPTION_PATH_LEN) :: Hname, Fename 

      if (di%MIM_options%have_MIM(2)) then
         Hname=trim(di%lc%sol%jaro%H_name)//'Average_mass'
         Fename=trim(di%lc%sol%jaro%Fe3_name)//'Average_mass'
      else
         Hname=trim(di%lc%sol%jaro%H_name)
         Fename=trim(di%lc%sol%jaro%Fe3_name)
      end if
      
      H => extract_scalar_field(di%state(2),Hname, stat=stat)
      if (.not. stat==0) then
        FLAbort('failed to extract the scalar field of H+ to calculate pH')
      end if
      
      Fe3 => extract_scalar_field(di%state(2),Fename, stat=stat)
      if (.not. stat==0) then
        FLAbort('failed to extract the scalar field of Fe3 to calculate jarosite precipitation')
      end if

      do i=1,di%number_pmesh_node         
         !if reactant Fe3 is near zero, stop reaction
         if ((Fe3%val(i)<=0.01) .or. (di%saturation(2)%ptr%val(i)<1.0e-8)) then
           di%lc%sol%jaro%dcdt%val(i)=0.0
           cycle
           
         end if

         pH=H%val(i)/1000.0
         pH=-DLOG10(pH)
         F=Fe3%val(i)*0.055845  !molar weight of Fe3 is 0.055845 kg/mole         
         F=DLOG10(F)
         rhs=-1.4319*pH+0.8679
         if (F>rhs) then
           dpre=(-di%lc%sol%jaro%rate)*Fe3%val(i)           
           !change the unit from (mole/m^3 solution/s) to ((mole/m^3 heap/s)
           dpre=di%porosity_pmesh%val(i)*di%saturation(2)%ptr%val(i)*dpre
           di%lc%sol%jaro%js%val(i)=di%lc%sol%jaro%js%val(i)+dpre*(-1.0/3.0)
         else
           dpre=0.0
         end if
         
         di%lc%sol%jaro%dcdt%val(i)=dpre

      end do

      nullify(H)
      nullify(Fe3)
      

   end subroutine calculate_jarosite_precipitation

   subroutine calculate_oxygen_dissolution(di)
        type(darcy_impes_type), intent(inout) :: di
        
        type(scalar_field), pointer ::og,o2,Tl,og_src,o2_src,o2_m,o2_im
        real :: theta1, theta2, T, ft,d_O, theta_m, theta_im, Fd, Fs,src
        integer :: stat, i,f
        character(len=OPTION_PATH_LEN) :: Oname

        if (di%MIM_options%have_MIM(2)) then
           Oname=trim(di%lc%sol%oxdi%o2_name)//'Average_mass'
           do f=1, size(di%generic_prog_sfield)
            if (di%generic_prog_sfield(f)%phase>1) then
               if (trim(di%generic_prog_sfield(f)%MIM%C_a%name)==Oname) then
                  stat=0
                  exit
               else 
                  stat=1
               end if
            end if
            
           end do
           if (.not. stat==0) then
              FLAbort('failed to extract the scalar field of mobile liquid phase oxygen to calculate oxygen dissolution')
           end if
           !The mobile, immobile field of O2
           o2_m => di%generic_prog_sfield(f)%sfield
           o2_im => di%generic_prog_sfield(f)%MIM%immobile_sfield%sfield
    
        else
           Oname=trim(di%lc%sol%oxdi%o2_name)
        end if
        
        og => extract_scalar_field(di%state(1),trim(di%lc%sol%oxdi%og_name), stat=stat)
        if (.not. stat==0) then
          FLAbort('failed to extract the scalar field of gas phase oxygen to calculate oxygen dissolution')
        end if        
        !extract the oxygen dissolution source term under og
        og_src  => extract_scalar_field(di%state(1),trim(di%lc%sol%oxdi%og_name)//'_Oxygen_dissolution_gas_phase', stat=stat)
        if (.not. stat==0) then
          FLAbort('failed to extract the scalar field of gas phase source term to calculate oxygen dissolution')
        end if
        o2 => extract_scalar_field(di%state(2),Oname, stat=stat)
        if (.not. stat==0) then
          FLAbort('failed to extract the scalar field of liquid  phase oxygen to calculate oxygen dissolution')
        end if
        !extract the oxygen dissolution source term under o2
        o2_src  => extract_scalar_field(di%state(2),trim(di%lc%sol%oxdi%o2_name)//'_Oxygen_dissolution_liquid_phase', stat=stat)
        if (.not. stat==0) then
          FLAbort('failed to extract the scalar field of liquid phase source term to calculate oxygen dissolution')
        end if
        
        Tl =>  di%lc%ht%liquid_temperature
        if (.not. stat==0) then
          FLAbort('failed to extract the scalar field of liquid  temperature to calculate oxygen dissolution')
        end if
        
        !loop over nodes
        ! if only single phase exist, no dissolution
        node_loop: do i=1,di%number_pmesh_node
          if ((di%saturation(1)%ptr%val(i)<=1.0e-8) .or. (di%saturation(2)%ptr%val(i)<=1.0e-8)) then
             di%lc%sol%oxdi%dcdt%val(i)=0.0
             og_src%val(i)=0.0
             o2_src%val(i)=0.0
             cycle
          end if

          !calculate fluid hold-up
          theta1=di%porosity_pmesh%val(i)*di%saturation(1)%ptr%val(i)
          theta2=di%porosity_pmesh%val(i)*di%saturation(2)%ptr%val(i)


          !calculate the equilibrium constant ft, o2=og*ft
          T=Tl%val(i)
          ft=(0.046*(T**2.0)+203.35*T*DLOG(T/298.0)-(299.378+0.092*T)*(T-298)-20.591*(10.0**3.0))/(T*8.3144)
          ft=0.08205746*T*EXP(ft) !0.08205746 is gas constant in the unit of L*atm*K^-1*mol^-1
          
          ! the transfer of oxygen from liquid phase to gas phase in unit of mole/m^3_heap
          d_O=(1.0/(ft/theta1+1.0/theta2))*(o2%val(i)-ft*og%val(i))

          di%lc%sol%oxdi%dcdt%val(i)=d_O
          !add the equilibrium term of oxygen dissolution to gas phase and liquid phase oxygen
          !for gas phase
          og_src%val(i)=d_O/theta1
          og%val(i)=og%val(i)+og_src%val(i)
          !for liquid phase, if MIM, this is the average mass of liquid oxygen, mole per total volume of solution
          o2_src%val(i)=-d_O/theta2

          o2%val(i)=o2%val(i)+o2_src%val(i)
          !for the mobile part of liquid oxygen
          if (di%MIM_options%have_MIM(2)) then
             Fd=di%generic_prog_sfield(f)%MIM%Fd%val(i)
             Fs=di%generic_prog_sfield(f)%MIM%Fs%val(i)
             theta_m=di%porosity_pmesh%val(i)*di%MIM_options%mobile_saturation(2)%ptr%val(i)
             theta_im=di%porosity_pmesh%val(i)*di%MIM_options%immobile_saturation(2)%ptr%val(i)
             !concentration in mole per volume of immobile liquid solution
             o2_im%val(i)=o2%val(i)*theta2*Fs/theta_im
             !concentration in mole per volume of mobile liquid solution
             o2_m%val(i)=o2%val(i)*theta2*Fd/theta_m

             !the source of mobile part in mole per volume of immobile liquid solution
             src=di%generic_prog_sfield(f)%MIM%chem%im_src%sfield%val(i)
             di%generic_prog_sfield(f)%MIM%chem%im_src%sfield%val(i)=src+o2_src%val(i)*theta2*Fs/theta_im
             !the source of mobile part in mole per volume of mobile liquid solution
             src=di%generic_prog_sfield(f)%MIM%chem%mo_src%sfield%val(i)
             di%generic_prog_sfield(f)%MIM%chem%mo_src%sfield%val(i)=src+o2_src%val(i)*theta2*Fd/theta_m
          end if
        end do node_loop
        nullify(og)
        nullify(og_src)
        nullify(o2)
        nullify(o2_src)
        nullify(Tl)
        if (di%MIM_options%have_MIM(2)) then
           nullify(o2_m)
           nullify(o2_im)
        end if
        
   end subroutine calculate_oxygen_dissolution
   
   subroutine calculate_ferrous_oxidation_bio_terms(di)
     type(darcy_impes_type), intent(inout) :: di
     integer :: i,f,j
     real :: Ta, ft,cl,cfe2,rho_h,lh,phi_l,phi_o,src_o,old_phi_o,Fd,lh_m,miu, Tas
     
     call zero(di%lc%sol%feox%bio%dcdt)

     node_loop: do i=1,di%number_pmesh_node
         
         !change the concentration fron mole/m3 to kg/m3 
         cl=di%lc%sol%feox%bio%cl%val(i)*0.016
         cfe2=di%lc%sol%feox%bio%cfe2%val(i)*0.056
         
         !-------------calculate f(t)-----------------
         !the temperature T is based on the average of liquid and ore
         if (di%lc%ht%heat_transfer_single) then
           Ta=di%lc%ht%liquid_temperature%val(i)
         else
           Ta=(di%lc%ht%rock_temperature%val(i)+di%lc%ht%liquid_temperature%val(i))/2.0
         end if
         
         phi_loop: do j=1,size(di%lc%sol%feox%bio%miu)
           
           Tas=Ta-di%lc%sol%feox%bio%T_shift(j)
           ft=21830090.0*Tas*exp(-7000.0/Tas)/(1.0+exp(236.0-74000.0/Tas))

           !---------------calculate miu----------------------
           miu=di%lc%sol%feox%bio%miu_max(j)*ft*(cl/(di%lc%sol%feox%bio%kmo+cl))*(cfe2/(di%lc%sol%feox%bio%kmfe2+cfe2))

           di%lc%sol%feox%bio%miu(j)%ptr%val(i)=miu
         
           !-----------------calculate bacteria source terms---------------
           rho_h=(1-di%porosity_pmesh%val(i))*di%lc%ht%rock_density%val(i)
           lh=di%porosity_pmesh%val(i)*di%saturation(2)%ptr%val(i)
           phi_l=di%lc%sol%feox%bio%phi_l(j)%ptr%val(i)
           phi_o=di%lc%sol%feox%bio%phi_ore(j)%ptr%val(i)

           !in unit of cell/m^3 heap
           di%lc%sol%feox%bio%phi_l_src(j)%ptr%val(i)=(miu-di%lc%sol%feox%bio%k_death(j))*phi_l*lh &
                                             -di%lc%sol%feox%bio%k1(j)*phi_l*lh*(1-(phi_o/di%lc%sol%feox%bio%phi_max(j))) &
                                             +di%lc%sol%feox%bio%k2(j)*rho_h*phi_o
          
           !in unit of cell/m^3 heap
           src_o=(miu-di%lc%sol%feox%bio%k_death(j))*rho_h*phi_o &
                                             +di%lc%sol%feox%bio%k1(j)*phi_l*lh*(1-(phi_o/di%lc%sol%feox%bio%phi_max(j))) &
                                             -di%lc%sol%feox%bio%k2(j)*rho_h*phi_o
          

           
            

           di%lc%sol%feox%bio%dcdt%val(i)=di%lc%sol%feox%bio%dcdt%val(i)-di%lc%sol%feox%bio%phi_l(j)%ptr%val(i)*miu/di%lc%sol%feox%bio%Y
           !-------------------------calculate new phi_ore-------------------------------------
           !initially, phi_ore=phi_max
           if (di%current_time<=1.0e-10) then
              old_phi_o=0.0!di%lc%sol%feox%bio%phi_max
           else
              old_phi_o=phi_o
           end if
         
           !in kg/m^3_ore
           if (src_o>=0.0) then
              di%lc%sol%feox%bio%phi_ore(j)%ptr%val(i)=old_phi_o+src_o*di%dt/rho_h
           else
              di%lc%sol%feox%bio%phi_ore(j)%ptr%val(i)=old_phi_o/(1.0-src_o*di%dt/(rho_h*old_phi_o))
           end if
           
         end do phi_loop
         
     end do node_loop
     
   end subroutine calculate_ferrous_oxidation_bio_terms

   subroutine calculate_leach_heat_transfer_src(di)
     type(darcy_impes_type), intent(inout) :: di

     !local variables
     integer :: i,n_s,n_l,nsrc_s,nsrc_l,nd,nr,n2_s,n2_l
     real :: Crs,Crl,prts,dTsl,ktp,Hs,Hr,miu

     !---------------for two phase heat transfer--------------------
     !count the number of sources term for solid temperature
     nsrc_s=size(di%lc%ht%two_phase_src_solid)
     nsrc_l=size(di%lc%ht%two_phase_src_liquid)

     node_loop: do i=1,di%number_pmesh_node
       prts=1.0-di%porosity_pmesh%val(i)
       Crs=di%lc%ht%rock_cp%val(i)*di%lc%ht%rock_density%val(i)
       Crl=di%lc%ht%liquid_cp%val(i)*di%density(2)%ptr%val(i)
       dTsl=di%lc%ht%rock_temperature%val(i)-di%lc%ht%liquid_temperature%val(i)
       ktp=di%lc%ht%K_eff_ls%val(i)*dTsl*prts
       
       !loop the liquid temperature source terms
       src_loop_l: do n_l=1,nsrc_l
         
         select case(di%lc%ht%two_phase_src_liquid(n_l)%ptr%name)

           !the solid-liquid heat transfer source
           case('solid_liquid_heat_transfer_liquid_phase')
             di%lc%ht%two_phase_src_liquid(n_l)%ptr%val(i)=ktp/Crl
           case('solution_phase_heat_sources')
             Hr=0.0
             nr=size(di%lc%ht%liquid_sr_src)
             solution_reaction_loop: do n2_l=1,nr
               Hr=Hr+di%lc%ht%liquid_sr_src(n2_l)%sr_src%val(i)*di%lc%ht%liquid_sr_src(n2_l)%Enthalpy            

             end do solution_reaction_loop
             di%lc%ht%two_phase_src_liquid(n_l)%ptr%val(i)=Hr/Crl  !in (k/s)(m^3 solution/m^3 heap)
           case default
             FLAbort("liquid_phase temperature heat transfer source " // di%lc%ht%two_phase_src_liquid(n_l)%ptr%name // " not found")

         end select    
       end do src_loop_l

       !loop the rock temperature source terms
       src_loop_r: do n_s=1,nsrc_s

         select case(di%lc%ht%two_phase_src_solid(n_s)%ptr%name)
           
           !the solid-liquid heat transfer source
           case('solid_liquid_heat_transfer_rock_phase')
             di%lc%ht%two_phase_src_solid(n_s)%ptr%val(i)=-ktp/Crs
             
           !the heat source from mineral dissolutiom  
           case('mineral_dissolution_heat_sources')
             Hs=0.0
             nd=size(di%lc%ht%rock_md_src)
             dissolution_loop: do n2_s=1,nd
               Hs=Hs+di%lc%ht%rock_md_src(n2_s)%md_src%val(i)*di%lc%ht%rock_md_src(n2_s)%Enthalpy
             end do dissolution_loop
             di%lc%ht%two_phase_src_solid(n_s)%ptr%val(i)=Hs/Crs
           case default
             FLAbort("Solid_phase temperature heat transfer source "// di%lc%ht%two_phase_src_solid(n_s)%ptr%name // " not found")
         end select

      end do src_loop_r

      if (di%lc%ht%rtss%have_rtss) then
        call add_rock_temperature_surface_heat_transfer_src(di,Crs)
      end if
      

     end do node_loop

   end subroutine calculate_leach_heat_transfer_src

   subroutine calculate_leach_rock_temperature(di)
     type(darcy_impes_type), intent(inout) :: di

     !local temperature
     integer ::i, nsrc_s,n_s
     real :: prts, dt
     type(scalar_field) :: src, src_cv_mass
     
     call allocate(src,di%pressure_mesh)
     call zero(src)   
     nsrc_s=size(di%lc%ht%two_phase_src_solid)

     src_loop: do n_s=1, nsrc_s
        call addto(src,di%lc%ht%two_phase_src_solid(n_s)%ptr)
     end do src_loop

     if (di%lc%ht%have_rock_temperature_src) call addto(src,di%lc%ht%rock_temperature_src)

     if (di%lc%ht%rtss%have_rtss) then
        do i=1,size(di%lc%ht%rtss%src)
           call addto(src,di%lc%ht%rtss%src(i)%ptr)
        end do
                
     end if     
     
     if (di%lcsub%have_leach_subcycle) then
           dt=di%lcsub%sub_dt
     else
           dt=di%dt
     end if
     node_loop: do i=1,di%number_pmesh_node

         prts=1.0-di%porosity_pmesh%val(i)

         src%val(i)=src%val(i)*dt/prts

     end do node_loop

     !calculate the new rock temperature
     call addto(di%lc%ht%rock_temperature,src)
     call deallocate(src)
   end subroutine calculate_leach_rock_temperature 

   subroutine allocate_leach_heat_transfer_prog_Temperature_src(di,f,shared_rhs,shared_lhs,isub)
       type(darcy_impes_type), intent(inout) :: di
       integer, intent(in) :: f
       type(scalar_field),optional, intent(inout) :: shared_rhs,shared_lhs
       integer,optional,intent(in) :: isub
       
       !local variables
       integer :: i,n,nsrc,p
       real :: isub_e,isub_s,theta_d,old_theta_d
       type(scalar_field) :: src, src_p,src_n, theta_m, theta_im, src_cv_mass
       
       !for the liquid phase temperature heat transfer source
       nsrc=size(di%lc%ht%two_phase_src_liquid)
       call allocate(src,di%pressure_mesh)
       call zero(src)

       if (di%generic_prog_sfield(f)%lh_src%src_linear%have) then
          call allocate(src_p,di%pressure_mesh)
          call allocate(src_n,di%pressure_mesh)
          call zero(src_p)
          call zero(src_n)        
       end if       
       
       p=di%generic_prog_sfield(f)%phase
       
       if (di%lcsub%have_leach_subcycle) then
          isub_e=real(isub)/real(di%lcsub%number_subcycle)
       else
          isub_e=1.0          
       end if
       
       src_loop: do n=1,nsrc
          if (di%generic_prog_sfield(f)%lh_src%src_linear%have) then
             if (minval(di%lc%ht%two_phase_src_liquid(n)%ptr%val)<0.0) then
                !add to nagtive source
                call addto(src_n,di%lc%ht%two_phase_src_liquid(n)%ptr) !in (k/s)(m^3 solution/m^3 heap)
             else                            
                !add to positive source
                call addto(src_p,di%lc%ht%two_phase_src_liquid(n)%ptr) !in (k/s)(m^3 solution/m^3 heap)
             end if
             
          end if
          
          call addto(src,di%lc%ht%two_phase_src_liquid(n)%ptr)!in (k/s)(m^3 solution/m^3 heap)
       end do src_loop
       !Add leaching chemical source term to rhs
       if (di%generic_prog_sfield(f)%MIM%have_MIM_source) then
          
          call allocate(theta_m,di%pressure_mesh)
          call allocate(theta_im,di%pressure_mesh)
         
          call zero(theta_m)
          call zero(theta_im)
          
          if (di%lcsub%have_leach_subcycle) then
             !mobile liquid hold up
             call set(theta_m,di%lcsub%sub_lht(p))
             call scale(theta_m,isub_e)
             call addto(theta_m,di%lcsub%old_sub_lht(p),scale=(1-isub_e))

             !immobile liquid hold up  
             call set(theta_im,di%lcsub%sub_lht_im(p))
             call scale(theta_im,isub_e)
             call addto(theta_im,di%lcsub%old_sub_lht_im(p),scale=(1-isub_e))
         
          else
             !mobile liquid hold up
             call set(theta_m, di%porosity_pmesh)
             call scale(theta_m, di%MIM_options%mobile_saturation(p)%ptr)
             !immobile liquid hold up
             call set(theta_im, di%porosity_pmesh)
             call scale(theta_im, di%MIM_options%immobile_saturation(p)%ptr)             
          end if
           
          !allocate the chemical src to the immobile part
          call invert(theta_im)
          call set(di%generic_prog_sfield(f)%MIM%chem%im_src%sfield, src)
          call scale(di%generic_prog_sfield(f)%MIM%chem%im_src%sfield,di%generic_prog_sfield(f)%MIM%Fs)
          call scale(di%generic_prog_sfield(f)%MIM%chem%im_src%sfield,theta_im) !in k/s
          
          !allocate the chemical src to the mobile part
          call invert(theta_m)
          call set(di%generic_prog_sfield(f)%MIM%chem%mo_src%sfield, src)
          call scale(di%generic_prog_sfield(f)%MIM%chem%mo_src%sfield,di%generic_prog_sfield(f)%MIM%Fd)
          call scale(di%generic_prog_sfield(f)%MIM%chem%mo_src%sfield,theta_m) !in k/s
          
          !for the positive and negative part when doing source linearization
          if (di%generic_prog_sfield(f)%MIM%chem%if_src_linear) then

             !immobile part
             call set(di%generic_prog_sfield(f)%MIM%chem%im_src%p_src,src_p)
             call set(di%generic_prog_sfield(f)%MIM%chem%im_src%n_src,src_n)
             call scale(di%generic_prog_sfield(f)%MIM%chem%im_src%p_src,di%generic_prog_sfield(f)%MIM%Fs)
             call scale(di%generic_prog_sfield(f)%MIM%chem%im_src%p_src, theta_im) !in k/s
             call scale(di%generic_prog_sfield(f)%MIM%chem%im_src%n_src,di%generic_prog_sfield(f)%MIM%Fs)
             call scale(di%generic_prog_sfield(f)%MIM%chem%im_src%n_src, theta_im) !in k/s

             !mobile part
             call set(di%generic_prog_sfield(f)%MIM%chem%mo_src%p_src,src_p)
             call set(di%generic_prog_sfield(f)%MIM%chem%mo_src%n_src,src_n)
             call scale(di%generic_prog_sfield(f)%MIM%chem%mo_src%p_src,di%generic_prog_sfield(f)%MIM%Fd)
             call scale(di%generic_prog_sfield(f)%MIM%chem%mo_src%p_src, theta_m)  !in k/s
             call scale(di%generic_prog_sfield(f)%MIM%chem%mo_src%n_src,di%generic_prog_sfield(f)%MIM%Fd)
             call scale(di%generic_prog_sfield(f)%MIM%chem%mo_src%n_src, theta_m)   !in k/s
          end if
        
          call deallocate(theta_m)
          call deallocate(theta_im)
                 
       end if

       !*************Add to RHS and LHS************************************************
       !If the immobile  field exist, solve the source term of Mobile-immobile model implicitly
       if (di%generic_prog_sfield(f)%MIM%have_MIM_source) then
          if (di%lcsub%have_leach_subcycle) then
             call darcy_trans_MIM_prog_sfield_allocate_rhs_lhs(di,p,f,shared_rhs, shared_lhs, isub)
          else
             call darcy_trans_MIM_prog_sfield_allocate_rhs_lhs(di,p,f,shared_rhs, shared_lhs)
          end if
          
       else
          !----------When the MIM is turned off----------------------
          !Add the source term to the RHS
          call allocate(src_cv_mass,di%pressure_mesh)
          !with src linearization
          if(di%generic_prog_sfield(f)%lh_src%src_linear%have) then
             call compute_cv_mass(di%positions,src_cv_mass,src_p)
             call addto(di%rhs,src_cv_mass) !S_positive

             node_loop3: do i=1,di%number_pmesh_node
                if (di%lcsub%have_leach_subcycle) then
                 
                   if (di%lcsub%iterated_sfield(f)%val(i)<=1.0D-15) then
                      src_n%val(i)=0.0
                   else
                      theta_d=di%lcsub%sub_lht(p)%val(i)*isub_e+di%lcsub%old_sub_lht(p)%val(i)*(1-isub_e)
                      old_theta_d=di%lcsub%sub_lht(p)%val(i)*isub_s+di%lcsub%old_sub_lht(p)%val(i)*(1-isub_s)
                      src_n%val(i)=-src_n%val(i)*theta_d/(di%lcsub%iterated_sfield(f)%val(i)*old_theta_d) !'-S_negative*theta_d/(C_old*theta_d_old)' 
                   end if

                else
                 
                   if (di%generic_prog_sfield(f)%sfield%val(i)<=1.0D-15) then
                      src_n%val(i)=0.0
                   else
                      theta_d=di%sat_ADE%val(i)*di%porosity_pmesh%val(i)
                      old_theta_d=di%old_sat_ADE%val(i)*di%old_porosity_pmesh%val(i)
                      src_n%val(i)=-src_n%val(i)*theta_d/(di%generic_prog_sfield(f)%old_sfield%val(i)*old_theta_d) !'-S_negative*thta_d/(C_old*theta_d_old' 
                   end if
                 
                end if
             end do node_loop3

             call compute_cv_mass(di%positions, src_cv_mass, src_n )
             call addto(di%lhs, src_cv_mass)
          else
             !without source linearizatiom
             call compute_cv_mass(di%positions,src_cv_mass,src)
             call addto(di%rhs,src_cv_mass)
       
          end if
          call deallocate(src_cv_mass)
          
       end if
     
       call deallocate(src)

       if (di%generic_prog_sfield(f)%lh_src%src_linear%have) then
          call deallocate(src_p)
          call deallocate(src_n)       
       end if   
       
   end subroutine allocate_leach_heat_transfer_prog_Temperature_src
   

     !-------------Add the chemical source terms to RHS for solving the prognostic fields----------------
   subroutine allocate_leaching_chemical_prog_sfield_src(di,f,shared_rhs,shared_lhs,isub)
     !!Add the chemical source terms to the rhs of the ADE
     !!If MIM, calculate the immobile term as well
      type(darcy_impes_type), intent(inout) :: di
      integer, intent(in) :: f
      type(scalar_field),optional, intent(inout) :: shared_rhs,shared_lhs
      integer,optional,intent(in) :: isub
      
      !local variables
      type(scalar_field) :: leach_src,leach_src_p,leach_src_n, single_src, theta_m, theta_im, src_cv_mass
      integer :: n,p,i,stat
      real :: s_factor !the stoichemistry factor
      real :: isub_e,isub_s, theta_d, old_theta_d
      character(len=FIELD_NAME_LEN) :: lc_name,phin
      type(scalar_field), pointer :: src => null()
      
      call allocate(leach_src,di%pressure_mesh)
      call zero(leach_src)
      
      call allocate(single_src,di%pressure_mesh)
      call zero(single_src)

      if (di%generic_prog_sfield(f)%lc_src%src_linear%have) then
         call allocate(leach_src_p,di%pressure_mesh)
         call allocate(leach_src_n,di%pressure_mesh)
         call zero(leach_src_p)
         call zero(leach_src_n)        
      end if
      
      p=di%generic_prog_sfield(f)%phase

      if (di%lcsub%have_leach_subcycle) then

         isub_e=real(isub)/real(di%lcsub%number_subcycle)
         isub_s=real(isub-1)/real(di%lcsub%number_subcycle) 
      else         
         isub_e=1.0
         isub_s=0.0
      end if

      
      !for the solution phase reactions
      if (di%generic_prog_sfield(f)%lc_src%have_sol_src) then

        do n=1, size(di%generic_prog_sfield(f)%lc_src%sfield_sol_src) 
          lc_name = di%generic_prog_sfield(f)%lc_src%sfield_sol_src(n)%lc_name
          s_factor = di%generic_prog_sfield(f)%lc_src%sfield_sol_src(n)%sto_factor
          call zero(single_src)

          select case(trim(lc_name))
             case("Ferrous_Oxidation")
               if (.not. associated(di%lc%sol%feox%dcdt)) &
               FLAbort('Ferrous_Oxidation is turned off in the leaching chemical model, while its source term is turned on under the prognostic scaler field') 
               src => di%lc%sol%feox%dcdt             
             
             case('Jarosite_Precipitation')
               if (.not. associated(di%lc%sol%jaro%dcdt)) &
               FLAbort('Jarosite_Precipitation is turned off in the leaching chemical model, while its source term is turned on under the prognostic scaler field')  
               src => di%lc%sol%jaro%dcdt
            
             case('Oxygen_dissolution_liquid_phase') 
               if (.not. associated(di%lc%sol%oxdi%dcdt)) &
               FLAbort('Oxygen_dissolution is turned off in the leaching chemical model, while its source term is turned on under the prognostic scaler field')
               cycle !calculate elsewhere, in the calculate oxygen dissolution subroutine
            
             case('Oxygen_dissolution_gas_phase')
               if (.not. associated(di%lc%sol%oxdi%dcdt)) &
               FLAbort('Oxygen_dissolution is turned off in the leaching chemical model, while its source term is turned on under the prognostic scaler field')
               cycle !calculate elsewhere, in the calculate oxygen dissolution subroutine
               
             case('Bacteria_Ferrous_oxidation')
                 phin=di%generic_prog_sfield(f)%sfield%name
                 src=> extract_scalar_field(di%state(2), trim(phin)//'_Bacteria_Ferrous_oxidation', stat=stat)
                 if (.not. stat==0) then
                     FLAbort('failed to extract the scaler field  Bacteria_Ferrous_oxidation source, check whether the number is correct')
                 end if
            case default
               FLAbort("Leaching chemical algorithm " // trim(lc_name) // " not found")
          end select
          
          
          call addto(single_src,src,s_factor)
          
          call addto(leach_src, single_src) !addto the chemical source term with scale of the stoichemistry factor

          if (di%generic_prog_sfield(f)%lc_src%src_linear%have) then
            if (minval(single_src%val)<0.0) then
              !add to nagtive source
              call addto(leach_src_n, single_src)
            else
              !add to positive source
              call addto(leach_src_p, single_src)
           end if

         end if
           
          node_loop1: do i=1,di%number_pmesh_node
            single_src%val(i)=single_src%val(i)/(di%porosity_pmesh%val(i)*di%saturation(p)%ptr%val(i))
            !this the reaction src based on the total averaged concentration
            di%generic_prog_sfield(f)%lc_src%sfield_sol_src(n)%sfield%val(i)=single_src%val(i) !mole/m^3_solution               
          end do node_loop1
          
        end do
        
      end if

      !for the mineral dissolution ractions
      if (di%generic_prog_sfield(f)%lc_src%have_dis_src) then
        do n=1, size(di%generic_prog_sfield(f)%lc_src%sfield_dis_src)
          lc_name = di%generic_prog_sfield(f)%lc_src%sfield_dis_src(n)%lc_name
          s_factor = di%generic_prog_sfield(f)%lc_src%sfield_dis_src(n)%sto_factor
          call zero(single_src)
          select case(trim(lc_name))
             case("CuFeS2_oxidation_aqueous_ferric_sulfate")
               if (.not. associated(di%lc%dis%chal%dcdt)) &
               FLAbort('CuFeS2_oxidation_aqueous_ferric_sulfate is turned off in the leaching chemical model,while its source term is turned on under the prognostic scaler field')
               src => di%lc%dis%chal%dcdt

             case('FeS2_oxidation_aqueous_ferric_sulfate')
               if (.not. associated(di%lc%dis%pyri%dcdt)) &
               FLAbort('FeS2_oxidation_aqueous_ferric_sulfate is turned off in the leaching chemical model, while its source term is turned on under the prognostic scaler field')
               src => di%lc%dis%pyri%dcdt 

             case('S0_dissolution')
               if (.not. associated(di%lc%dis%sulf%dcdt)) &
               FLAbort('S0_dissolution is turned off in the leaching chemical model,while its source term is turned on under the prognostic scaler field')
               src => di%lc%dis%sulf%dcdt

             case('Gangue_mineral_acid_dissolution')
               if (.not. associated(di%lc%dis%gang%dcdt)) &
               FLAbort('gangue mineral dissolution is turned off in the leaching chemical model,while its source term is turned on under the prognostic scaler field')
               src => di%lc%dis%gang%dcdt
               
             case default
               FLAbort("Leaching chemical algorithm " // trim(lc_name) // " not found")
          end select
          
          call addto(single_src,src,s_factor)
          call addto(leach_src, single_src) !addto the chemical source term with scale of the stoichemistry factor, mole per volume of heap/s

          if (di%generic_prog_sfield(f)%lc_src%src_linear%have) then
            if (minval(single_src%val)<0.0) then
              !add to nagtive source
              call addto(leach_src_n, single_src)
            else
              !add to positive source
              call addto(leach_src_p, single_src)
           end if

         end if
         
          node_loop2: do i=1,di%number_pmesh_node
            single_src%val(i)=single_src%val(i)/(di%porosity_pmesh%val(i)*di%saturation(p)%ptr%val(i))
            !this the reaction src based on the total averaged concentration
            di%generic_prog_sfield(f)%lc_src%sfield_dis_src(n)%sfield%val(i)=single_src%val(i) !mole/m^3_solution/s
            
          end do node_loop2 
          
        end do
        
     end if
     
      !Add leaching chemical source term to rhs
      if (di%generic_prog_sfield(f)%MIM%have_MIM_source) then
         call allocate(theta_m,di%pressure_mesh)
         call allocate(theta_im,di%pressure_mesh)
         call zero(theta_m)
         call zero(theta_im)

         if (di%lcsub%have_leach_subcycle) then
             !mobile liquid hold up
             call set(theta_m,di%lcsub%sub_lht(p))
             call scale(theta_m,isub_e)
             call addto(theta_m,di%lcsub%old_sub_lht(p),scale=(1-isub_e))

             !immobile liquid hold up  
             call set(theta_im,di%lcsub%sub_lht_im(p))
             call scale(theta_im,isub_e)
             call addto(theta_im,di%lcsub%old_sub_lht_im(p),scale=(1-isub_e))
          
         else
             !mobile liquid hold up
             call set(theta_m, di%porosity_pmesh)
             call scale(theta_m, di%MIM_options%mobile_saturation(p)%ptr)
             !immobile liquid hold up
             call set(theta_im, di%porosity_pmesh)
             call scale(theta_im, di%MIM_options%immobile_saturation(p)%ptr)             
         end if
       
         !allocate the chemical src to the immobile part
         call invert(theta_im)
         call set(di%generic_prog_sfield(f)%MIM%chem%im_src%sfield, leach_src)
         call scale(di%generic_prog_sfield(f)%MIM%chem%im_src%sfield,di%generic_prog_sfield(f)%MIM%Fs)
         call scale(di%generic_prog_sfield(f)%MIM%chem%im_src%sfield, theta_im) !in mole per volume of immobile liquid/s
         
         !allocate the chemical src to the mobile part
         call invert(theta_m)                         
         call set(di%generic_prog_sfield(f)%MIM%chem%mo_src%sfield, leach_src)!in mole per volume of heap/s
         call scale(di%generic_prog_sfield(f)%MIM%chem%mo_src%sfield,di%generic_prog_sfield(f)%MIM%Fd)
         call scale(di%generic_prog_sfield(f)%MIM%chem%mo_src%sfield, theta_m) !in mole per volume of mobile liquid/s
         !for the positive and negative part when doing source linearization
         if (di%generic_prog_sfield(f)%MIM%chem%if_src_linear) then
            call zero(di%generic_prog_sfield(f)%MIM%chem%im_src%p_src)
            call zero(di%generic_prog_sfield(f)%MIM%chem%im_src%n_src)
            call zero(di%generic_prog_sfield(f)%MIM%chem%mo_src%p_src)
            call zero(di%generic_prog_sfield(f)%MIM%chem%mo_src%n_src)

            !immobile part
            call set(di%generic_prog_sfield(f)%MIM%chem%im_src%p_src,leach_src_p)
            call set(di%generic_prog_sfield(f)%MIM%chem%im_src%n_src,leach_src_n)
            call scale(di%generic_prog_sfield(f)%MIM%chem%im_src%p_src,di%generic_prog_sfield(f)%MIM%Fs)
            call scale(di%generic_prog_sfield(f)%MIM%chem%im_src%p_src, theta_im) !in mole per volume of immobile liquid/s
            call scale(di%generic_prog_sfield(f)%MIM%chem%im_src%n_src,di%generic_prog_sfield(f)%MIM%Fs)
            call scale(di%generic_prog_sfield(f)%MIM%chem%im_src%n_src, theta_im) !in mole per volume of immobile liquid/s

            !mobile part
            call set(di%generic_prog_sfield(f)%MIM%chem%mo_src%p_src,leach_src_p)
            call set(di%generic_prog_sfield(f)%MIM%chem%mo_src%n_src,leach_src_n)
            call scale(di%generic_prog_sfield(f)%MIM%chem%mo_src%p_src,di%generic_prog_sfield(f)%MIM%Fd)
            call scale(di%generic_prog_sfield(f)%MIM%chem%mo_src%p_src,theta_m) !in mole per volume of mobile liquid/s
            call scale(di%generic_prog_sfield(f)%MIM%chem%mo_src%n_src,di%generic_prog_sfield(f)%MIM%Fd)
            call scale(di%generic_prog_sfield(f)%MIM%chem%mo_src%n_src,theta_m) !in mole per volume of mobile liquid/s
         end if
         
         call deallocate(theta_im)
         call deallocate(theta_m)
      end if
      
      !*************Add to RHS and LHS************************************************
      !If the immobile  field exist, solve the source term of Mobile-immobile model implicitly
      if (di%generic_prog_sfield(f)%MIM%have_MIM_source) then
         if (di%lcsub%have_leach_subcycle) then
            call darcy_trans_MIM_prog_sfield_allocate_rhs_lhs(di,p,f,shared_rhs, shared_lhs, isub)
         else
            call darcy_trans_MIM_prog_sfield_allocate_rhs_lhs(di,p,f,shared_rhs, shared_lhs)
         end if
      else
         
         !----------When the MIM is turned off----------------------
         !Add the source term to the RHS
         call allocate(src_cv_mass,di%pressure_mesh)
         !with src linearization
         if(di%generic_prog_sfield(f)%lc_src%src_linear%have) then
            call compute_cv_mass(di%positions,src_cv_mass,leach_src_p)
            call addto(di%rhs,src_cv_mass) !S_positive
            !linearize the negative source and add to lhs
            node_loop3: do i=1,di%number_pmesh_node
              if (di%lcsub%have_leach_subcycle) then
                 
                 if (di%lcsub%iterated_sfield(f)%val(i)<=1.0D-15) then
                    leach_src_n%val(i)=0.0
                 else
                    theta_d=di%lcsub%sub_lht(p)%val(i)*isub_e+di%lcsub%old_sub_lht(p)%val(i)*(1-isub_e)
                    old_theta_d=di%lcsub%sub_lht(p)%val(i)*isub_s+di%lcsub%old_sub_lht(p)%val(i)*(1-isub_s)
                    leach_src_n%val(i)=-leach_src_n%val(i)*theta_d/(di%lcsub%iterated_sfield(f)%val(i)*old_theta_d) !'-S_negative*thta_d/(C_old*theta_d_old' 
                 end if
 
              else
                 
                 
                 if (di%generic_prog_sfield(f)%sfield%val(i)<=1.0D-15) then
                    leach_src_n%val(i)=0.0
                 else
                    theta_d=di%sat_ADE%val(i)*di%porosity_pmesh%val(i)
                    old_theta_d=di%old_sat_ADE%val(i)*di%old_porosity_pmesh%val(i)
                    leach_src_n%val(i)=-leach_src_n%val(i)*theta_d/(di%generic_prog_sfield(f)%sfield%val(i)*old_theta_d) !'-S_negative*thta_d/(C_old*theta_d_old'
                    
                 end if
                 
              end if
            end do node_loop3                   
            call compute_cv_mass(di%positions, src_cv_mass, leach_src_n )
            call addto(di%lhs, src_cv_mass)
            call compute_cv_mass(di%positions,src_cv_mass,leach_src_p)
            call addto(di%rhs,src_cv_mass)
            
         else
            !without source linearizatiom
            call compute_cv_mass(di%positions,src_cv_mass,leach_src)
            call addto(di%rhs,src_cv_mass)
         
         end if
         call deallocate(src_cv_mass)
         
      end if
  
      call deallocate(leach_src)
      call deallocate(single_src)

      if (di%generic_prog_sfield(f)%lc_src%src_linear%have) then
         call deallocate(leach_src_p)
         call deallocate(leach_src_n)
      end if

      nullify(src)

    end subroutine allocate_leaching_chemical_prog_sfield_src

    subroutine calculate_solid_liquid_wetting_efficiency(di, node, wet_eff)
      type(darcy_impes_type), intent(inout) :: di
      integer, intent(in) :: node
      real, intent(out) :: wet_eff
      
      real:: Ga, Re,u,mu,rho,g,poro,d

      !the velocity magnitude of darcy flux
      u=norm2(node_val(di%darcy_velocity(2)%ptr,node))

      !the gravity magnitude
      g=di%gravity_magnitude
      if (g<=1.0e-10) then
         g=9.8
      end if
     
      !the liquid viscosity
      mu=node_val(di%lviscosity_pmesh,node)

      !liquid density
      rho=node_val(di%density(2)%ptr,node)

      !heap porosity
      poro=node_val(di%porosity_pmesh,node)

      !rock particle diameter
      d=node_val(di%heap%rock_d,node)

      !Galileo number
      Ga=((d*poro/(1.0-poro))**3.0)*((rho/mu)**2.0)*g

      !Reynolds number
      Re=u*rho*d/mu/(1.0-poro)
      
      wet_eff=1.104*(Re**(1.0/3.0))*((1/Ga)**(1.0/9.0))

      call set(di%lc%wet_eff%wet_eff,node,wet_eff)
    end subroutine calculate_solid_liquid_wetting_efficiency

    subroutine initialize_rock_temperature_surface_source(di)
      type(darcy_impes_type), intent(inout) :: di
      integer :: ns, f, stat
      character(len=OPTION_PATH_LEN) :: option_path, src_name,path2,type_name
      integer :: shape_option(2)

      option_path=('/Leaching_chemical_model/heat_transfer_model/two_phases_heat_transfer/scalar_field::Rock_Temperature/scalar_field')
      ns=0
      do f=1, option_count(trim(option_path))
         call get_option(trim(option_path)//'['//int2str(f-1)//']/name', src_name)
         if (src_name(1:28)=='Rock_Temperature_Surface_Src') then
            ns=ns+1
         end if
         
      end do   
      
      if (ns>0) then
        di%lc%ht%rtss%have_rtss=.true.
        allocate(di%lc%ht%rtss%ids(ns))
        allocate(di%lc%ht%rtss%transfer_type(ns))      
        allocate(di%lc%ht%rtss%coefficient(ns))
        allocate(di%lc%ht%rtss%ref_t(ns))
        allocate(di%lc%ht%rtss%A(ns))
        allocate(di%lc%ht%rtss%src(ns))

        ns=0
        do f= 1, option_count(trim(option_path))
           call get_option(trim(option_path)//'['//int2str(f-1)//']/name', src_name)

         if (src_name(1:28)=='Rock_Temperature_Surface_Src') then
            ns=ns+1
            path2=trim(option_path)//'::'//trim(src_name)
            di%lc%ht%rtss%src(ns)%ptr=> extract_scalar_field(di%state(1),trim(src_name), stat=stat)
            if (.not. stat==0) then
               FLAbort('failed to extract' //trim(src_name))
            end if
            
            ! Get vector of surface ids
            shape_option=option_shape(trim(path2)//'/Boundary/surface_ID')
            allocate(di%lc%ht%rtss%ids(ns)%list(1:shape_option(1)))
            call get_option(trim(path2)//'/Boundary/surface_ID',di%lc%ht%rtss%ids(ns)%list)
            call get_option(trim(path2)//'/Boundary/heat_transfer_type[0]/name',type_name)
            select case(trim(type_name))
               case('heat_convection')
                    di%lc%ht%rtss%transfer_type(ns)=1
                    call get_option(trim(path2)//'/Boundary/heat_transfer_type::heat_convection/heat_transfer_coefficient',di%lc%ht%rtss%coefficient(ns))
                    call get_option(trim(path2)//'/Boundary/heat_transfer_type::heat_convection/surface_area',di%lc%ht%rtss%A(ns))
                    call get_option(trim(path2)//'/Boundary/heat_transfer_type::heat_convection/reference_temperature',di%lc%ht%rtss%ref_t(ns))

                 case('radiation')
                    di%lc%ht%rtss%transfer_type(ns)=2
                    call get_option(trim(path2)//'/Boundary/heat_transfer_type::radiation/emissivity',di%lc%ht%rtss%coefficient(ns))
                    call get_option(trim(path2)//'/Boundary/heat_transfer_type::radiation/surface_area',di%lc%ht%rtss%A(ns))
                    call get_option(trim(path2)//'/Boundary/heat_transfer_type::radiation/reference_temperature',di%lc%ht%rtss%ref_t(ns))
              
             end select
                  
            end if
            
         end do
         call allocate_rock_temperature_surface_numbers(di)
      end if
         
               

    end subroutine initialize_rock_temperature_surface_source

    subroutine finalize_rock_temperature_surface_source(di)
      type(darcy_impes_type), intent(inout) :: di
      integer :: i, ns
      
      ns=size(di%lc%ht%rtss%src)
      do i=1,ns
         nullify(di%lc%ht%rtss%src(i)%ptr)
         deallocate(di%lc%ht%rtss%ids(i)%list)
      end do
      deallocate(di%lc%ht%rtss%src)
      deallocate(di%lc%ht%rtss%ids)
      deallocate(di%lc%ht%rtss%sele)
      deallocate(di%lc%ht%rtss%transfer_type)      
      deallocate(di%lc%ht%rtss%coefficient)
      deallocate(di%lc%ht%rtss%ref_t)
      deallocate(di%lc%ht%rtss%A)
      di%lc%ht%rtss%have_rtss=.false.
        
    end subroutine finalize_rock_temperature_surface_source

    subroutine allocate_rock_temperature_surface_numbers(di)
      type(darcy_impes_type), intent(inout) :: di
      integer, dimension(:,:), allocatable :: sele_list
      integer ::  i,j,s_count,k,ns

      ns=size(di%lc%ht%rtss%ids)
      allocate(sele_list(ns*di%number_sele,2))
      s_count=0
      sele_loop: do i=1, di%number_sele
         do j=1,ns
  
            do k=1,size(di%lc%ht%rtss%ids(j)%list)
               
              if (di%lc%ht%rtss%ids(j)%list(k)==surface_element_id(di%lc%ht%rock_temperature, i)) then
                s_count=s_count+1
                !the first column is source number, the second column is surface element number
                sele_list(s_count,1)=j
                sele_list(s_count,2)=i
              end if
            end do
          
         end do
                    
      end do sele_loop
      allocate(di%lc%ht%rtss%sele(s_count,2))
     
      di%lc%ht%rtss%sele=sele_list(1:s_count,:)

    end subroutine allocate_rock_temperature_surface_numbers
    
    subroutine add_rock_temperature_surface_heat_transfer_src(di,Crs)
      type(darcy_impes_type), intent(inout) :: di
      real, intent(in) :: Crs
      
      integer :: i, num, ns, sele
      integer, dimension(:), pointer :: p_nodes_bdy
      real,    dimension(:),  pointer :: ghost_temp, ghost_src
      type(scalar_field) ::  cv_mass_pressure_mesh_with_source
      
      allocate(ghost_temp(face_loc(di%pressure_mesh,1)))
      allocate(ghost_src(face_loc(di%pressure_mesh,1)))
      allocate(p_nodes_bdy(face_loc(di%pressure_mesh,1)))
      
      ns=size(di%lc%ht%rtss%sele(:,1))

      !initialize the source
      do i=1, size(di%lc%ht%rtss%src)
        call zero(di%lc%ht%rtss%src(i)%ptr)         
      end do
      
      do i=1, ns
         num=di%lc%ht%rtss%sele(i,1)
         sele=di%lc%ht%rtss%sele(i,2)
         p_nodes_bdy = face_global_nodes(di%pressure_mesh, sele)
         ghost_temp=face_val(di%lc%ht%rock_temperature, sele)
         select case(di%lc%ht%rtss%transfer_type(num))
            
              case(1)
                !heat convection                 
                ghost_src=-(ghost_temp-di%lc%ht%rtss%ref_t(num))*di%lc%ht%rtss%A(num)*di%lc%ht%rtss%coefficient(num)/Crs
                call addto(di%lc%ht%rtss%src(num)%ptr, p_nodes_bdy, ghost_src)
              
              case(2)
                !radiation
                ghost_src=-(ghost_temp**4.0-di%lc%ht%rtss%ref_t(num)**4.0)*di%lc%ht%rtss%A(num)*di%lc%ht%rtss%coefficient(num)*5.6703D-11/Crs
                call addto(di%lc%ht%rtss%src(num)%ptr, p_nodes_bdy, ghost_src) 
                
         end select
                          
      end do          
              
    end subroutine add_rock_temperature_surface_heat_transfer_src
    
    
    
   !---------------------some accessory subroutines----------------------------------------

   subroutine cubic_spline_interpolation(a,b,c,d,x_data,x,y)
      real, dimension(:), intent(in) :: a,b,c,d,x_data
      real, intent(in) :: x
      real, intent(out) :: y
      
      integer :: n,ic,jc,pc

      n=size(x_data)

      if (x<=x_data(1)) then
        y=b(1)*(x-x_data(1))+a(1)
      else if (x>x_data(n)) then
        y=b(n)*(x-x_data(n))+a(n)
      else
        ic=0
        jc=n
        do while (jc-ic>1)
         if (mod(jc-ic,2)==0) then
            pc=(jc-ic)/2
         else
            pc=(jc-ic+1)/2
         end if
       
         if (x>x_data(ic+pc)) then
            ic=ic+pc
         else
            jc=ic+pc    
         end if
        end do
        y=a(ic)+b(ic)*(x-x_data(ic))+c(ic)*((x-x_data(ic))**2.0)+d(ic)*((x-x_data(ic))**3.0)
      end if
   end subroutine  cubic_spline_interpolation
   
   subroutine cubic_spline_coefficient(e,dedt_k,spline_coefficient)

      integer :: i
      real, intent(in) :: e(:), dedt_k(:) !the extraction and dextraction_over_k from experiment data
      !the parameters for spline interpolation
      integer :: n,n1,n2
      real, intent(inout) :: spline_coefficient(:,:)!cubic spline coefficients and the last node of data
     
      !parameters pass to thomas algrithm
      real :: dh(size(e)-1), dy(size(e)-1) !the interval between each two nodes
      real :: c(size(e))
      
      n=size(e)
      n1=n-1
      n2=n-2
     
      do i=1, n-1 
        dh(i)=e(i+1)-e(i)
        dy(i)=dedt_k(i+1)-dedt_k(i)
      end do
      
     
      call solve_tridiagonal_matrix(n2,dh,dy,c) 
     
      do i=1, n1
        !the zero order term coefficient
        spline_coefficient(1,i)=dedt_k(i)
        !the 1st order term coefficient
        spline_coefficient(2,i)=dy(i)/dh(i)-c(i)*(dh(i)/2.0)-(dh(i)/6.0)*(c(i+1)-c(i))
        !the 2nd order term coefficient
        spline_coefficient(3,i)=c(i)/2.0
        !the 3rd order term coefficient
        spline_coefficient(4,i)=(c(i+1)-c(i))/(6.0*dh(i))
      end do
      
      spline_coefficient(1,n)=dedt_k(n)
      spline_coefficient(2,n)=(dedt_k(n)-dedt_k(n-1))/(e(n)-e(n-1))
      spline_coefficient(3:4,n)=0.0D0
        
   end subroutine cubic_spline_coefficient

   subroutine solve_tridiagonal_matrix(n2,dh,dy,c)
     !solve tridiagonal matrix by Thomas Algrithm
     integer, intent(in) :: n2
     real, intent(in) :: dh(:), dy(:)
     real, intent(out) :: c(n2+2)

     integer :: j
     real :: x(n2)
     real :: coefficient(4,n2) !coefficient for matrix and vector
                               !a_(i)*x_(i-1)+b_(i)*x_(i)+c_(i)*x(i+1)=d_(i)

     do j=1,n2
        !coefficient a,b,c,d
        coefficient(1,j)=dh(j)
        coefficient(2,j)=2.0*(dh(j)+dh(j+1))
        coefficient(3,j)=dh(j+1)
        coefficient(4,j)=6.0*(dy(j+1)/dh(j+1)-dy(j)/dh(j))
     end do

     do j=1,n2
       if (j == 1) then
         coefficient(3,j)=coefficient(3,j)/coefficient(2,j)
         coefficient(4,j)=coefficient(4,j)/coefficient(2,j)
       else
         coefficient(3,j)=coefficient(3,j)/(coefficient(2,j)-coefficient(3,j-1)*coefficient(1,j))
         coefficient(4,j)=(coefficient(4,j)-coefficient(4,j-1)*coefficient(1,j))/&
                           (coefficient(2,j)-coefficient(3,j-1)*coefficient(1,j))
       end if
     end do

     do j=n2,1,-1
      if (j == n2) then
        x(j)=coefficient(4,j)
      else
        x(j)=(coefficient(4,j)-coefficient(3,j)*x(j+1))
      end if
     end do

     c=[0.0D0,x,0.0D0]

  end subroutine solve_tridiagonal_matrix   

end  module
