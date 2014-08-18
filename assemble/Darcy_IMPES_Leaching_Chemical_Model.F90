!This file contains all the subroutines of leaching chemical model

#include "fdebug.h"

module darcy_impes_leaching_chemical_model

  use spud
  use fields
  use state_module
  use fldebug
  use field_options
  use spud
  use fields_manipulation
  use data_structures
  use initialise_fields_module
  use global_parameters, only: OPTION_PATH_LEN
  use darcy_impes_leaching_types
  use fefields, only: compute_cv_mass

  use darcy_impes_assemble_type, only: darcy_impes_type
  
  implicit none
  private
  
  public :: initialize_leaching_chemical_model, &
            finalize_leaching_chemical_model, &
            add_leach_chemical_prog_src_to_rhs, &
            calculate_leaching_chemical_model
  

  
  
  
  
  contains
  
  
  subroutine initialize_leaching_chemical_model(di)
     !!initialize the leaching chemical model
      
     type(darcy_impes_type), intent(inout) :: di
     
  
     !local parameter
     integer :: p,f,flc, ns, nd, nb,fb, stat,ndata,nshape(2)
     character(len=OPTION_PATH_LEN) :: option_path, reaction_name, path_l
        !---------------------allocate the fields in the chemical model-------------
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
                  di%lc%sol%feox%dcdt => extract_scalar_field(di%state(1), 'feox_dFe2_dt', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named feox_dFe2_dt')
                  end if
                  !get activation energy
                  path_l = trim(option_path)//'::'//trim(reaction_name)
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
                   
                  
                case('Jarosite_Precipitation')
                  di%lc%sol%jaro%dcdt => extract_scalar_field(di%state(1), 'jaro_dM_dt', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named jaro_dM_dt')
                  end if

                case('Oxygen_dissolution')
                  di%lc%sol%oxdi%dcdt => extract_scalar_field(di%state(1), 'oxdi_dOg_dt', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named oxdi_dOg_dt')
                  end if

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
                  di%lc%dis%chal%dcdt => extract_scalar_field(di%state(1), 'chal_dCuFeS2_dt', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named chal_dCuFeS2_dt')
                  end if
                  !get extraction rate
                  di%lc%dis%chal%ex_r => extract_scalar_field(di%state(1), 'chal_extraction_rate', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named chal_extraction_rate')
                  end if
                  !get current extraction
                  di%lc%dis%chal%ex => extract_scalar_field(di%state(1), 'chal_current_extraction', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named chal_current_extraction')
                  end if
                  !get the molar concentraction of the mineral, mole per volume of heap
                  di%lc%dis%chal%mc => extract_scalar_field(di%state(1), 'chal_molar_concentration', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named chal_molar_concentration')
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
                  di%lc%dis%pyri%dcdt => extract_scalar_field(di%state(1), 'pyri_dFeS2_dt', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named pyri_dFeS2_dt')
                  end if
                  !get extraction rate
                  di%lc%dis%pyri%ex_r => extract_scalar_field(di%state(1), 'pyri_extraction_rate', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named pyri_extraction_rate')
                  end if
                  !get current extraction
                  di%lc%dis%pyri%ex => extract_scalar_field(di%state(1), 'pyri_current_extraction', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named pyri_current_extraction')
                  end if
                  !get the molar concentraction of the mineral, mole per volume of heap
                  di%lc%dis%pyri%mc => extract_scalar_field(di%state(1), 'pyri_molar_concentration', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named pyri_molar_concentration')
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
                  di%lc%dis%sulf%dcdt => extract_scalar_field(di%state(1), 'sulf_dS0_dt', stat=stat)
                  if (.not. stat==0) then
                    FLAbort('failed to extract the leaching reaction named sulf_dS0_dt')
                  end if
                
                case default
                  FLAbort("Leaching chemical algorithm " // trim(reaction_name) // " not found")

             end select
             
           end do   
        else 
           di%lc%have_dis=.false.      

        end if

        !----------allocate the generic prognostic leaching source terms-------------       
        !--loop over phase
        do f=1, size(di%generic_prog_sfield)
          option_path=di%generic_prog_sfield(f)%sfield%option_path
          p=di%generic_prog_sfield(f)%phase
          !---check for solution phase source-----------
          ns = option_count(trim(option_path)//'/prognostic/LeachingChemicalSourceTerm/SolutionPhaseSource/scalar_field')
          if (.not. ns==0) then
            di%generic_prog_sfield(f)%lc_src%have_sol_src = .true.
            allocate(di%generic_prog_sfield(f)%lc_src%sfield_sol_src(ns))

            !extract the name of the chemical reaction and stoichemistry factor
            do flc=1, ns 
              call get_option(trim(option_path)//'/prognostic/LeachingChemicalSourceTerm/&
                              SolutionPhaseSource/scalar_field['//int2str(flc-1)//']/name', & 
                              di%generic_prog_sfield(f)%lc_src%sfield_sol_src(flc)%lc_name)
              
              call get_option(trim(option_path)//'/prognostic/LeachingChemicalSourceTerm/&
                              SolutionPhaseSource/scalar_field['//int2str(flc-1)//']/&
                              diagnostic/stoichiometric_factor', &
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
                              MineralDissolutionSource/scalar_field['//int2str(flc-1)//']/name', &
                              di%generic_prog_sfield(f)%lc_src%sfield_dis_src(flc)%lc_name)

              call get_option(trim(option_path)//'/prognostic/LeachingChemicalSourceTerm/&
                              MineralDissolutionSource/scalar_field['//int2str(flc-1)//']/&
                              diagnostic/stoichiometric_factor', &
                              di%generic_prog_sfield(f)%lc_src%sfield_dis_src(flc)%sto_factor)
              di%generic_prog_sfield(f)%lc_src%sfield_dis_src(flc)%sfield => extract_scalar_field(di%state(p), &
                                                               trim(di%generic_prog_sfield(f)%sfield%name)//"_"&
                                           //trim(di%generic_prog_sfield(f)%lc_src%sfield_dis_src(flc)%lc_name))
            end do


          else
            di%generic_prog_sfield(f)%lc_src%have_dis_src = .false.
          
          endif

        end do

        !-------------allocate the temperature-----------------------
        !check the heat transfer model
        if (have_option('/Leaching_chemical_model/heat_transfer_model')) then
          if (have_option('/Leaching_chemical_model/heat_transfer_model/single_phase_heat_trasfer')) then
            di%lc%ht%heat_transfer_single = .true.
            di%lc%ht%heat_transfer_two = .false.
            di%lc%ht%heat_transfer_three = .false.
          elseif (have_option('/Leaching_chemical_model/heat_transfer_model/two_phases_heat_trasfer')) then
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
        
        di%lc%ht%liquid_temperature => extract_scalar_field(di%state(2), 'Temperature', stat=stat)
        if (.not. stat==0) then
           FLAbort('failed to extract the liquid temperature in phase 2 for leaching chemical model')
        end if

        if (di%lc%ht%heat_transfer_two) then
          di%lc%ht%air_temperature => extract_scalar_field(di%state(1), 'Temperature', stat=stat)
          if (.not. stat==0) then
            FLAbort('failed to extract the air temperature in phase 1 for leaching chemical model')
          end if

        elseif (di%lc%ht%heat_transfer_three) then
          di%lc%ht%rock_temperature => extract_scalar_field(di%state(1), 'Rock_Temperature', stat=stat)
          if (.not. stat==0) then
             FLAbort('failed to extract the rock temperature in heat transfer model for leaching chemical model')
          end if
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
              deallocate(di%lc%sol%feox%ak%bulk)

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
              deallocate(di%lc%dis%chal%ak%bulk)
              deallocate(di%lc%dis%chal%spline_coe)
              deallocate(di%lc%dis%chal%exp_ex)
              deallocate(di%lc%dis%chal%exp_exrk)

           case('FeS2_oxidation_aqueous_ferric_sulfate')
              nullify(di%lc%dis%pyri%dcdt)
              nullify(di%lc%dis%pyri%ex_r)
              nullify(di%lc%dis%pyri%ex)
              nullify(di%lc%dis%pyri%mc)
              deallocate(di%lc%dis%pyri%ak%bulk)
              deallocate(di%lc%dis%pyri%spline_coe)
              deallocate(di%lc%dis%pyri%exp_ex)
              deallocate(di%lc%dis%pyri%exp_exrk)

           case('S0_dissolution')
              nullify(di%lc%dis%sulf%dcdt)
              
         end select
           
       end do

     end if

     !deallocate leaching prognostic source field 
     do f=1, size(di%generic_prog_sfield)
       
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
       nullify(di%lc%ht%liquid_temperature)
       
       if (di%lc%ht%heat_transfer_two) nullify(di%lc%ht%air_temperature)
       if (di%lc%ht%heat_transfer_three) nullify(di%lc%ht%rock_temperature)

       di%lc%ht%heat_transfer_single = .false.
       di%lc%ht%heat_transfer_two = .false.
       di%lc%ht%heat_transfer_three = .false.
 
     end if

   end subroutine finalize_leaching_chemical_model


   !********The following are the subroutines to calculate fields for the chemical model****************
   
   !-------------Add the chemical source terms to RHS for solving the prognostic fields----------------
   subroutine add_leach_chemical_prog_src_to_rhs(di,f)
      type(darcy_impes_type), intent(inout) :: di
      integer, intent(in) :: f
      
      !local variables
      type(scalar_field) :: leach_src, single_src
      integer :: n
      real :: s_factor !the stoichemistry factor
      character(len=FIELD_NAME_LEN) :: lc_name
      type(scalar_field), pointer :: src => null()
      
      call allocate(leach_src,di%pressure_mesh)
      call zero(leach_src)
      
      call allocate(single_src,di%pressure_mesh)

      !for the solution phase reactions
      if (di%generic_prog_sfield(f)%lc_src%have_sol_src) then

        do n=1, size(di%generic_prog_sfield(f)%lc_src%sfield_sol_src) 
          lc_name = di%generic_prog_sfield(f)%lc_src%sfield_sol_src(n)%lc_name
          s_factor = di%generic_prog_sfield(f)%lc_src%sfield_sol_src(n)%sto_factor
          call zero(single_src)
          
          select case(trim(lc_name))
             case("Ferrous_Oxidation")
               if (.not. associated(di%lc%sol%feox%dcdt)) &
               FLAbort('Ferrous_Oxidation is turned off in the leaching chemical model, &
               while its source term is turned on under the prognostic scaler field') 
               src => di%lc%sol%feox%dcdt
               call addto(single_src,src,s_factor)
               call addto(leach_src, single_src) !addto the chemical source term with scale of the stoichemistry factor
             
             case('Jarosite_Precipitation')
               if (.not. associated(di%lc%sol%jaro%dcdt)) &
               FLAbort('Jarosite_Precipitation is turned off in the leaching chemical model, &
               while its source term is turned on under the prognostic scaler field')  
               src => di%lc%sol%jaro%dcdt
               call addto(single_src,src,s_factor)
               call addto(leach_src, single_src) !addto the chemical source term with scale of the stoichemistry factor
            
             case('Oxygen_dissolution_liquid_phase') 
               if (.not. associated(di%lc%sol%oxdi%dcdt)) &
               FLAbort('Oxygen_dissolution is turned off in the leaching chemical model, &
               while its source term is turned on under the prognostic scaler field')
               src => di%lc%sol%oxdi%dcdt
               call addto(single_src,src,s_factor)
               call addto(leach_src, single_src) !addto the chemical source term with scale of the stoichemistry factor
            
             case('Oxygen_dissolution_gas_phase')
               if (.not. associated(di%lc%sol%oxdi%dcdt)) &
               FLAbort('Oxygen_dissolution is turned off in the leaching chemical model, &
               while its source term is turned on under the prognostic scaler field')
               src => di%lc%sol%oxdi%dcdt
               call addto(single_src,src,s_factor)
               call addto(leach_src, single_src) !addto the chemical source term with scale of the stoichemistry factor
            
             case default
               FLAbort("Leaching chemical algorithm " // trim(lc_name) // " not found")
          end select
          
          call set(di%generic_prog_sfield(f)%lc_src%sfield_sol_src(n)%sfield, single_src)
          
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
               FLAbort('CuFeS2_oxidation_aqueous_ferric_sulfate is turned off in the leaching chemical model, &
               while its source term is turned on under the prognostic scaler field')
               src => di%lc%dis%chal%dcdt
               call addto(single_src,src,s_factor)
               call addto(leach_src, single_src) !addto the chemical source term with scale of the stoichemistry factor
 

             case('FeS2_oxidation_aqueous_ferric_sulfate')
               if (.not. associated(di%lc%dis%pyri%dcdt)) &
               FLAbort('FeS2_oxidation_aqueous_ferric_sulfate is turned off in the leaching chemical model, &
               while its source term is turned on under the prognostic scaler field')
               src => di%lc%dis%pyri%dcdt 
               call addto(single_src,src,s_factor)
               call addto(leach_src, single_src) !addto the chemical source term with scale of the stoichemistry factor

             case('S0_dissolution')
               if (.not. associated(di%lc%dis%sulf%dcdt)) &
               FLAbort('S0_dissolution is turned off in the leaching chemical model, &
               while its source term is turned on under the prognostic scaler field')
               src => di%lc%dis%sulf%dcdt
               call addto(single_src,src,s_factor)
               call addto(leach_src, single_src) !addto the chemical source term with scale of the stoichemistry factor

             case default
               FLAbort("Leaching chemical algorithm " // trim(lc_name) // " not found")
          end select
          
          call set(di%generic_prog_sfield(f)%lc_src%sfield_dis_src(n)%sfield, single_src)
          
        end do
        
      end if
      
      !Add leaching chemical source term to rhs
      call compute_cv_mass(di%positions, di%cv_mass_pressure_mesh_with_source, leach_src)
      call addto(di%rhs, di%cv_mass_pressure_mesh_with_source)

      call deallocate(leach_src)
      call deallocate(single_src)

      nullify(src)
   end subroutine add_leach_chemical_prog_src_to_rhs
   

   !-------------calculate leaching chemical model----------------------------------------------!
   subroutine calculate_leaching_chemical_model(di)
      type(darcy_impes_type), intent(inout) :: di

      !for mineral dissolution
      if (di%lc%have_dis) then            
         !--------for the chalcopyrite dissolution-----------------------------------!
         ! Chalcopyrite oxidation 
         if (associated(di%lc%dis%chal%dcdt)) call calculate_mineral_dissolution_semi_empirical_model(di%state,&
                                          di%lc%ht%liquid_temperature,di%number_pmesh_node,di%dt,di%lc%dis%chal)
         !pyrite dissolution
         if (associated(di%lc%dis%pyri%dcdt)) call calculate_mineral_dissolution_semi_empirical_model(di%state,&
                                          di%lc%ht%liquid_temperature,di%number_pmesh_node,di%dt,di%lc%dis%pyri)        
      end if
     
      contains 

        subroutine calculate_mineral_dissolution_semi_empirical_model(states,temperature,node_number,dt,mineral)
              type(leaching_semi_empirical_model_type), intent(inout) :: mineral
              type(state_type),dimension(:), intent(in) :: states
              integer, intent(in) :: node_number
              type(scalar_field), intent(in) :: temperature
              real, intent(in) :: dt !the time step
              real :: T_l,mc, k_rate, ext, ext_rk, ext_r, dcdt !temperature, molar concentration of the mineral, tate constant,current extraction 
                                                      !the extraction rate with k, the extraction rate, the concentration change rate
              type(scalar_field_pointer), dimension(:), allocatable :: cb !the reacting bulk species
              real, dimension(:),allocatable ::  cb_n, m ! the single node val of reacting bulk species and order of reaction
              real, dimension(:), pointer :: a,b,c,d !spline coefficient
              integer, dimension(:),allocatable :: p !the phase of the reacting species
              character(len=FIELD_NAME_LEN), dimension(:),allocatable :: cb_name !the name of the reacting species
              integer :: nspecies, isp, node, stat
             
              !get the number of species which take part in reaction
              nspecies = size(mineral%ak%bulk)             
              allocate(cb(nspecies),cb_n(nspecies), m(nspecies), p(nspecies), cb_name(nspecies))
              
              !spline coefficient
              a => mineral%spline_coe(1,:)
              b => mineral%spline_coe(2,:)
              c => mineral%spline_coe(3,:)
              d => mineral%spline_coe(4,:)

              do isp = 1, nspecies
                m(isp) = mineral%ak%bulk(isp)%order
                cb_name(isp) = mineral%ak%bulk(isp)%lc_name
                p(isp) = mineral%ak%bulk(isp)%phase
                cb(isp)%ptr => extract_scalar_field(states(p(isp)), trim(cb_name(isp)), stat)
                if (.not. stat==0) then
                  FLAbort('failed to extract the reacting species')
                end if
              end do
              
              !Calculate exctraction rate, extraction, and concentration change rate for each node
              !node loop
              node_loop: do node=1, node_number
                   
                   !the reaction temperature is based on the liquid temperature
                   T_l=node_val(temperature, node)

                   !the molar concentration of the mineral at the unit of mole per volumn of heap
                   mc= node_val(mineral%mc, node)

                   !calculate rate constant, which is k=e^(Ea/(R*T))
                   k_rate=EXP(mineral%ak%ae/(mineral%ak%gc*T_l))

                   !the current extraction of the element node
                   ext = node_val(mineral%ex,node)

                   !if extraction of the mineral is nearly one, not reacted
                   if ((1.0-ext)<=1.0D-15) then
                     ext_r = 0.0 
                     cycle node_loop                  
                   end if

                   do isp= 1, nspecies
                      cb_n(isp) = node_val(cb(isp)%ptr, node)
                      !if either of the reacting species is near zero, not reated
                      if  (cb_n(isp)<=1.0D-15) then 
                        if (m(isp) > 0.0) then
                            ext_r = 0.0  ! the species with positive order is the reactant, stop reaction
                            cycle node_loop
                        else
                            cb_n(isp) = 1.0 !the species with negative order is the product
                                            !let it equal to 1 and make the reaction independent of it
                                            !this might not be true, but zero product concentration 
                                            !normally only happen at beginning of reaction, so might only 
                                            !result in inaccuracy of the initial reaction 
                                                                                      
                        end if
                       
                      end if
                      k_rate = k_rate*(cb_n(isp)**m(isp))
                   end do
                   
                   !do cubic spline interpolation for current extraction rate with k
                   call cubic_spline_interpolation(a,b,c,d,mineral%exp_ex,ext,ext_rk) 
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

              end do node_loop

              !finalize
              !nullify before deallocate
              nullify(a,b,c,d)
              do isp = 1, nspecies
                nullify(cb(isp)%ptr)
              end do

              deallocate(cb)
              deallocate(cb_n)
              deallocate(m)
              deallocate(p)
              deallocate(cb_name)
          

        end subroutine calculate_mineral_dissolution_semi_empirical_model

   end subroutine calculate_leaching_chemical_model
    
   !---------------------some accessory subroutines----------------------------------------

   subroutine cubic_spline_interpolation(a,b,c,d,x_data,x,y)
      real, dimension(:), intent(in) :: a,b,c,d,x_data
      real, intent(in) :: x
      real, intent(out) :: y
      
      integer :: n,ic,jc,pc

      n=size(x_data)

      if (x<x_data(1)) then
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
