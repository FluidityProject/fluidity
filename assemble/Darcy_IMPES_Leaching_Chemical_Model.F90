!This file contains all the subroutines of leaching chemical model

#include "fdebug.h"

module darcy_impes_leaching_chemical_model

  use spud
  use fields
  use state_module
  use fldebug
  use global_parameters, only: OPTION_PATH_LEN
  use darcy_impes_leaching_types
  use fefields, only: compute_cv_mass

  use darcy_impes_assemble_type, only: darcy_impes_type
  
  implicit none
  private
  
  public :: initialize_leaching_chemical_model, &
            finalize_leaching_chemical_model, &
            add_leach_chemical_prog_src_to_rhs
  

  
  
  
  
  contains
  
  subroutine initialize_leaching_chemical_model(di)
     !!initialize the leaching chemical model
      
     type(darcy_impes_type), intent(inout) :: di
     
  
     !local parameter
     integer :: f,flc, ns, nd    
     character(len=OPTION_PATH_LEN) :: option_path
     
     !if there is a leaching chemical model, allocate the generic prognostic chemical source terms
     if (have_option('/Leaching_Chemical_Model')) then
        di%lc%have_leach_chem_model= .true.
        
        !--loop over phase
        do f=1, size(di%generic_prog_sfield)
          option_path=di%generic_prog_sfield(f)%sfield%option_path
          !---check for solution phase source-----------
          ns = option_count(trim(option_path)//'/prognostic/LeachingChemicalSourceTerm/SolutionPhaseSource')
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
            end do
     
          else
            di%generic_prog_sfield(f)%lc_src%have_sol_src = .false.

          end if
          
          !-----check for mineral dissolution phase source------
          nd = option_count(trim(option_path)//'/prognostic/LeachingChemicalSourceTerm/MineralDissolutionSource')
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
            end do


          else
            di%generic_prog_sfield(f)%lc_src%have_dis_src = .false.
          
          endif

        end do

     else 
        di%lc%have_leach_chem_model= .false.
     
     end if
  end subroutine

  !--------------------------------------------------------------------------------------------------------
  subroutine finalize_leaching_chemical_model(di) 
     
     !finalize terms from leaching_chemical_model

     type(darcy_impes_type), intent(inout) :: di
     
     !local variables
     integer :: f


     di%lc%have_leach_chem_model= .False.
       
     do f=1, size(di%generic_prog_sfield)
       
       if (di%generic_prog_sfield(f)%lc_src%have_sol_src) then
         di%generic_prog_sfield(f)%lc_src%have_sol_src = .false.
         deallocate(di%generic_prog_sfield(f)%lc_src%sfield_sol_src)
       end if
       
       if (di%generic_prog_sfield(f)%lc_src%have_dis_src) then
         di%generic_prog_sfield(f)%lc_src%have_dis_src = .false.
         deallocate(di%generic_prog_sfield(f)%lc_src%sfield_dis_src)
       end if

     end do
   
   end subroutine


   !********The following are the subroutines to calculate fields for the chemical model****************
   
   !-------------Add the chemical source terms to RHS for solving the prognostic fields----------------
   subroutine add_leach_chemical_prog_src_to_rhs(di,f)
      type(darcy_impes_type), intent(inout) :: di
      integer, intent(in) :: f
      
      !local variables
      type(scalar_field) :: leach_src
      integer :: n
      real :: s_factor !the stoichemistry factor
      character(len=FIELD_NAME_LEN) :: lc_name
      type(scalar_field), pointer :: src => null()
      
      call allocate(leach_src,di%pressure_mesh)
      call zero(leach_src)

      !for the solution phase reactions
      if (di%generic_prog_sfield(f)%lc_src%have_sol_src) then

        do n=1, size(di%generic_prog_sfield(f)%lc_src%sfield_sol_src) 
          lc_name = di%generic_prog_sfield(f)%lc_src%sfield_sol_src(n)%lc_name
          s_factor = di%generic_prog_sfield(f)%lc_src%sfield_sol_src(n)%sto_factor 
          src => extract_scalar_field(di%state(1), trim(lc_name))
          call addto(leach_src, src, s_factor) !addto the chemical source term with scale of the stoichemistry factor
        end do

      end if

      !for the mineral dissolution ractions
      if (di%generic_prog_sfield(f)%lc_src%have_dis_src) then
         
        do n=1, size(di%generic_prog_sfield(f)%lc_src%sfield_dis_src)
          lc_name = di%generic_prog_sfield(f)%lc_src%sfield_dis_src(n)%lc_name
          s_factor = di%generic_prog_sfield(f)%lc_src%sfield_dis_src(n)%sto_factor
          src => extract_scalar_field(di%state(1), trim(lc_name))
          call addto(leach_src, src, s_factor) !addto the chemical source term with scale of the stoichemistry factor
        end do

      end if
      
      !Add leaching chemical source term to rhs
      call compute_cv_mass(di%positions, di%cv_mass_pressure_mesh_with_source, leach_src)
      call addto(di%rhs, di%cv_mass_pressure_mesh_with_source)

      call deallocate(leach_src)

      nullify(src)
   end subroutine

end  module
