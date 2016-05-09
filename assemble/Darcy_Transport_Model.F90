!This file contains the subroutines of darcy transport model for leaching

#include "fdebug.h"

module darcy_transport_model

  use spud
  use fields
  use state_module
  use fldebug
  use field_options
  use spud
  use data_structures
  use initialise_fields_module
  use global_parameters, only: OPTION_PATH_LEN
  use vtk_cache_module, only:vtk_cache_finalise
  use fefields, only: compute_cv_mass
  use sparsity_patterns
  use solvers
  
  use darcy_impes_assemble_type, only: darcy_impes_type
  use darcy_impes_leaching_types


  implicit none
  private
  
  public :: darcy_trans_solve_MIM_saturations_and_mass_transfer_coefficient, &
       leaching_MIM_calculate_fields_and_ratio, &
       darcy_trans_MIM_prog_sfield_allocate_rhs_lhs, &
       darcy_trans_assemble_and_solve_immobile_sfield, &
       initialize_MIM_model,&       
       finalize_MIM_model,&
       finalize_heap_property,&
       initialize_heap_property

contains

  subroutine initialize_MIM_model(di)
      type(darcy_impes_type), intent(inout) :: di

      integer :: p, stat, f_count, f
      character(len=OPTION_PATH_LEN) :: tmp_char_option, tmp_char_option1, lalgorithm, path
      type(vector_field), pointer :: position =>null()

      ! Allocate the MIM model
      allocate(di%MIM_options%immobile_saturation(di%number_phase))
      allocate(di%MIM_options%old_immobile_saturation(di%number_phase))
      allocate(di%MIM_options%mobile_saturation(di%number_phase))
      allocate(di%MIM_options%old_mobile_saturation(di%number_phase))
      allocate(di%MIM_options%mass_trans_coef(di%number_phase))
      allocate(di%MIM_options%old_mass_trans_coef(di%number_phase)) 
      !flag of MIM model and mass transfer coefficient
      allocate(di%MIM_options%have_MIM(di%number_phase))
      allocate(di%MIM_options%have_mass_trans_coef(di%number_phase))

      do p = 1,di%number_phase
         if (p > 1) then                              
            di%MIM_options%immobile_saturation(p)%ptr  => extract_scalar_field(di%state(p), "ImmobileSaturation", stat = stat)
            if (stat == 0) then               
               di%MIM_options%have_MIM(p) = .true.
               di%MIM_options%old_immobile_saturation(p)%ptr  => extract_scalar_field(di%state(p), "OldImmobileSaturation")
               di%MIM_options%mobile_saturation(p)%ptr        => extract_scalar_field(di%state(p), "MobileSaturation")
               di%MIM_options%old_mobile_saturation(p)%ptr    => extract_scalar_field(di%state(p), "OldMobileSaturation")
                
    
            else
               di%MIM_options%have_MIM(p) = .false.
               di%MIM_options%immobile_saturation(p)%ptr  => di%constant_zero_sfield_pmesh
               di%MIM_options%mobile_saturation(p)%ptr    => di%constant_zero_sfield_pmesh
               di%MIM_options%old_immobile_saturation(p)%ptr  => di%constant_zero_sfield_pmesh
               di%MIM_options%old_mobile_saturation(p)%ptr    => di%constant_zero_sfield_pmesh

            end if
            
         else                       
             ! Cannot have MIM for phase 1 and dual phase
             di%MIM_options%have_MIM(p) = .false.
             di%MIM_options%immobile_saturation(p)%ptr  => di%constant_zero_sfield_pmesh
             di%MIM_options%mobile_saturation(p)%ptr    => di%constant_zero_sfield_pmesh
             di%MIM_options%old_immobile_saturation(p)%ptr  => di%constant_zero_sfield_pmesh
             di%MIM_options%old_mobile_saturation(p)%ptr    => di%constant_zero_sfield_pmesh
                         
         end if

         !If the MIM exists, check whether there is mass transfer coefficient
         if (di%MIM_options%have_MIM(p) ) then
            di%MIM_options%mass_trans_coef(p)%ptr => extract_scalar_field(di%state(p), "MassTransferCoefficient", stat = stat)
   
            di%MIM_options%old_mass_trans_coef(p)%ptr   => extract_scalar_field(di%state(p), "OldMassTransferCoefficient", stat = stat)
            if (stat == 0) then
                di%MIM_options%have_mass_trans_coef(p) = .true.
            else
                di%MIM_options%have_mass_trans_coef(p) = .false.
                di%MIM_options%mass_trans_coef(p)%ptr  => di%constant_zero_sfield_pmesh
                di%MIM_options%old_mass_trans_coef(p)%ptr => di%constant_zero_sfield_pmesh
            end if
         else
            !Cannot have mass transfer coefficient without MIM model
             di%MIM_options%have_mass_trans_coef(p) = .false.
            di%MIM_options%mass_trans_coef(p)%ptr  => di%constant_zero_sfield_pmesh
            di%MIM_options%old_mass_trans_coef(p)%ptr  => di%constant_zero_sfield_pmesh  
         end if

         path=trim('/material_phase['//int2str(p-1)//']'//'/MobileImmobileModel/scalar_field')
         if (have_option(trim(path)//'::ImmobileSaturation/diagnostic')) then
            call get_option(trim(path)//'::ImmobileSaturation/diagnostic/algorithm/name', lalgorithm, default = "Internal")
            
            select case(trim(lalgorithm))
               case("Internal")
                  di%MIM_options%Lima_immobile_sat=.true.
                  
               case("scalar_python_diagnostic")
                  di%MIM_options%Lima_immobile_sat=.false.

               case default
                  FLAbort("The internal algorithm is now the only diagnostic algorithm available to immobile saturation")

            end select
                                
         end if
         
         if (have_option(trim(path)//'::MassTransferCoefficient/diagnostic')) then
            call get_option(trim(path)//'::MassTransferCoefficient/diagnostic/algorithm/name', lalgorithm, default = "Internal")
            
            select case(trim(lalgorithm))
               case("Internal")
                  di%MIM_options%Lima_mass_trans=.true.
                  
               case("scalar_python_diagnostic")
                  di%MIM_options%Lima_mass_trans=.false.

               case default
                  FLAbort("The internal algorithm is now the only diagnostic algorithm available to the mass transfer coefficient")

           end select
                                
         end if
         
      end do

      !the flag to check wether the MIM exist in at least one phase
      di%MIM_options%have_MIM_phase = .false.
      do p = 1, di%number_phase
         if (di%MIM_options%have_MIM(p)) then
           di%MIM_options%have_MIM_phase = .true.
           exit
        end if
      end do 
      ewrite(1,*) 'if have MIM ', di%MIM_options%have_MIM_phase

      if (size(di%generic_prog_sfield) > 0) then
         f_count = 0
         do p = 1, di%number_phase

            do f = 1, option_count('/material_phase['//int2str(p-1)//']/scalar_field')

               if (have_option('/material_phase['//int2str(p-1)//']/scalar_field['//int2str(f-1)//']/prognostic')) then

                  call get_option('/material_phase['//int2str(p-1)//']/scalar_field['//int2str(f-1)//']/name', &
                                  tmp_char_option)

                  if ((trim(tmp_char_option) /= 'Pressure') .and. &
                       (trim(tmp_char_option) /= 'Saturation')) then

                     f_count = f_count + 1

                     if (di%MIM_options%have_MIM(p)) then
                        if (.not.(have_option('/material_phase['//int2str(p-1)//']/scalar_field['//int2str(f-1)//']/prognostic/scalar_field::Immobile'))) then
                           FLExit('The immobile field of field'//tmp_char_option// 'should be turned on')
                        end if
                        
                     end if
                     
                     if (have_option('/material_phase['//int2str(p-1)//']/scalar_field['//int2str(f-1)//']/prognostic/scalar_field::Immobile'))  then
                        if (.not.(di%MIM_options%have_MIM(p))) then
                        FLExit('The Mobileimmobile model of phase'//int2str(p)//'should be turned on')
                           
                        end if

                        di%generic_prog_sfield(f_count)%MIM%have_MIM_source = .true.
                        ! allocate the MIM sorce terms for the matrix to solve the prognostic sfield
                        tmp_char_option1 = trim(tmp_char_option)//'Immobile'

                        di%generic_prog_sfield(f_count)%MIM%immobile_sfield%sfield => extract_scalar_field(di%state(p), &
                             tmp_char_option1, stat)

                        if (.not. stat==0) then
                           FLExit('failed to extract the reacting species immobile field:'//tmp_char_option1)
                        end if
                        di%generic_prog_sfield(f_count)%MIM%immobile_sfield%old_sfield => extract_scalar_field(di%state(p), &
                             'Old'//trim(tmp_char_option1))

                        !get the initial condition of the immobile field
                        position => extract_vector_field(di%state(1), "Coordinate")
                        call zero(di%generic_prog_sfield(f_count)%MIM%immobile_sfield%sfield)
                        call initialise_field_over_regions(di%generic_prog_sfield(f_count)%MIM%immobile_sfield%sfield, &
                        trim(di%generic_prog_sfield(f_count)%MIM%immobile_sfield%sfield%option_path)//'/diagnostic/initial_condition', position)

                        call vtk_cache_finalise()

                        nullify(position)

                        !The average mass
                        tmp_char_option1 = trim(tmp_char_option)//'Average_mass'
                        di%generic_prog_sfield(f_count)%MIM%C_a =>extract_scalar_field(di%state(p), &
                             tmp_char_option1, stat)
 
                        if (.not. stat==0) then
                           FLExit('failed to extract the average mass of the  field:'//tmp_char_option)
                        end if

                        di%generic_prog_sfield(f_count)%MIM%old_C_a => extract_scalar_field(di%state(p), &
                             'Old'//trim(tmp_char_option1))

                         !the ratio of mobile and immobile
                        tmp_char_option1 = trim(tmp_char_option)//'Mobile_ratio'
                        di%generic_prog_sfield(f_count)%MIM%Fd => extract_scalar_field(di%state(p), &
                             trim(tmp_char_option1), stat)
                        if (.not. stat==0) then
                           FLExit('failed to extract the mobile ratio of the field:'//tmp_char_option)
                        end if

                        tmp_char_option1 = trim(tmp_char_option)//'Immobile_ratio'
                        di%generic_prog_sfield(f_count)%MIM%Fs => extract_scalar_field(di%state(p), &
                             trim(tmp_char_option1), stat)
                        if (.not. stat==0) then
                           FLExit('failed to extract the Immobile ratio of the field:'//tmp_char_option)
                        end if

                         !-----If the chemical leaching source terms exist, the MIM under the chemical source should also be turned on
                        if (have_option('/material_phase['//int2str(p-1)//']/scalar_field['//int2str(f-1)//']/prognostic/leaching_temperature_sources')) then
                           if (.not. have_option('/material_phase['//int2str(p-1)//']/scalar_field['//int2str(f-1)//']/prognostic/leaching_temperature_sources/Mobile_Immobile_Model')) then
                              FLExit('please turn on the Mobile_Immobile_Model under leaching_temperature_sources of field:'// tmp_char_option)
                           end if
                           
                           di%generic_prog_sfield(f_count)%MIM%chem%have_chem = .true.
                           tmp_char_option1 = trim(tmp_char_option)//'Mobile_chemical_src'
                           di%generic_prog_sfield(f_count)%MIM%chem%mo_src%sfield => extract_scalar_field(di%state(p), &
                                trim(tmp_char_option1), stat)
                           
                           if (.not. stat==0) then
                              FLExit('failed to extract the field:'//tmp_char_option1)
                           end if  

                           tmp_char_option1 = trim(tmp_char_option)//'Immobile_chemical_src'
                           di%generic_prog_sfield(f_count)%MIM%chem%im_src%sfield => extract_scalar_field(di%state(p), &
                                trim(tmp_char_option1), stat)                       
                           if (.not. stat==0) then
                              FLExit('failed to extract the field:'//tmp_char_option1)
                           end if

                        elseif (have_option('/material_phase['//int2str(p-1)//']/scalar_field['//int2str(f-1)//']/prognostic/LeachingChemicalSourceTerm')) then
                           if (.not. have_option('/material_phase['//int2str(p-1)//']/scalar_field['//int2str(f-1)//']/prognostic/LeachingChemicalSourceTerm/Mobile_Immobile_Model')) then
                              FLExit('please turn on the Mobile_Immobile_Model under leaching_temperature_sources of field:'// tmp_char_option)
                           end if

                           di%generic_prog_sfield(f_count)%MIM%chem%have_chem = .true.
                           tmp_char_option1 = trim(tmp_char_option)//'Mobile_chemical_src'
                           di%generic_prog_sfield(f_count)%MIM%chem%mo_src%sfield => extract_scalar_field(di%state(p), &
                                trim(tmp_char_option1), stat)
                           
                           if (.not. stat==0) then
                              FLExit('failed to extract the field:'//tmp_char_option1)
                           end if  

                           tmp_char_option1 = trim(tmp_char_option)//'Immobile_chemical_src'
                           di%generic_prog_sfield(f_count)%MIM%chem%im_src%sfield => extract_scalar_field(di%state(p), &
                                trim(tmp_char_option1), stat)
                           
                           if (.not. stat==0) then
                              FLExit('failed to extract the field:'//tmp_char_option1)
                           end if

                        else
                           di%generic_prog_sfield(f_count)%MIM%chem%have_chem=.false.
                       
                        end if

                         else
                        di%generic_prog_sfield(f_count)%MIM%have_MIM_source = .false.    
                     end if
                                      

                  end if


               end if


            end do

         end do
         
      end if
      
      
 end subroutine initialize_MIM_model

 subroutine initialize_heap_property(di)
   type(darcy_impes_type), intent(inout) :: di
   integer :: stat
   
   if (have_option(trim('/Leaching_chemical_model/liquid_solid_wetting_efficiency')) .or. di%MIM_options%Lima_immobile_sat .or. di%MIM_options%Lima_mass_trans ) then
      di%heap%rock_d=>extract_scalar_field(di%state(1), trim('Rock_diameter'), stat=stat)
      if (.not. stat==0) then
         FLAbort('please specify the rock diameter under porous media if you want to calculate the wetting efficiency or immobile saturation and mass transfer coeffcient diagnostically')         
      end if
   end if
   
 end subroutine initialize_heap_property
 

 subroutine finalize_heap_property(di)
   type(darcy_impes_type), intent(inout) :: di

   if (associated(di%heap%rock_d)) then
     nullify(di%heap%rock_d)
   end if
      
 end subroutine finalize_heap_property
 
 subroutine finalize_MIM_model(di)

   type(darcy_impes_type), intent(inout) :: di

   ! local variables
   integer :: p, f

   do p = 1,di%number_phase

      nullify(di%MIM_options%old_immobile_saturation(p)%ptr)
      nullify(di%MIM_options%immobile_saturation(p)%ptr)
      nullify(di%MIM_options%mobile_saturation(p)%ptr)
      nullify(di%MIM_options%old_mobile_saturation(p)%ptr)
      nullify(di%MIM_options%mass_trans_coef(p)%ptr)
      nullify(di%MIM_options%old_mass_trans_coef(p)%ptr)
      di%MIM_options%have_MIM(p)=.false.
   end do
   
   ! Deallocate the MIM model
   di%MIM_options%Lima_immobile_sat=.false.
   di%MIM_options%Lima_mass_trans=.false.
   di%MIM_options%have_mass_trans_coef = .false.
   di%MIM_options%have_MIM = .false.
   di%MIM_options%have_MIM_phase = .false.
   deallocate(di%MIM_options%immobile_saturation)
   deallocate(di%MIM_options%old_immobile_saturation)
   deallocate(di%MIM_options%old_mobile_saturation)
   deallocate(di%MIM_options%mobile_saturation)
   deallocate(di%MIM_options%mass_trans_coef)
   deallocate(di%MIM_options%have_MIM)
   deallocate(di%MIM_options%have_mass_trans_coef)
   deallocate(di%MIM_options%old_mass_trans_coef)

   if (size(di%generic_prog_sfield) > 0) then
      do f = 1, size(di%generic_prog_sfield)
         if ( di%generic_prog_sfield(f)%MIM%have_MIM_source) then
           nullify(di%generic_prog_sfield(f)%MIM%immobile_sfield%sfield)
           nullify(di%generic_prog_sfield(f)%MIM%immobile_sfield%old_sfield)
           nullify(di%generic_prog_sfield(f)%MIM%C_a)
           nullify(di%generic_prog_sfield(f)%MIM%old_C_a)
           nullify(di%generic_prog_sfield(f)%MIM%Fs)
           nullify(di%generic_prog_sfield(f)%MIM%Fd)
           di%generic_prog_sfield(f)%MIM%have_MIM_source  = .false.
           if (di%generic_prog_sfield(f)%MIM%chem%have_chem) then
              nullify(di%generic_prog_sfield(f)%MIM%chem%mo_src%sfield)
              nullify(di%generic_prog_sfield(f)%MIM%chem%im_src%sfield)
              di%generic_prog_sfield(f)%MIM%chem%have_chem=.false.
           end if
        end if        
      end do
   end if
   
  
   
 end subroutine finalize_MIM_model


 !Slove the mobile saturation if MIM exist
 subroutine darcy_trans_solve_MIM_saturations_and_mass_transfer_coefficient(di)

        type(darcy_impes_type), intent(inout) :: di
        integer :: i,node
        type(scalar_field), pointer :: total_sat => null()  !total saturation 
        type(scalar_field), pointer :: immobile_sat  => null()  ! immobile saturation
        type(scalar_field)  :: mobile_sat   !mobile saturation
        real:: Re,u,mu,rho,g,d,alpha
        
        call allocate(mobile_sat, di%pressure_mesh)
     
        do i= 2, di%number_phase

          total_sat      => di%saturation(i)%ptr
          immobile_sat   => di%MIM_options%immobile_saturation(i)%ptr

          if (di%MIM_options%have_MIM(i)) then

             if (di%MIM_options%Lima_immobile_sat .or. di%MIM_options%Lima_mass_trans) then
                do node=1,di%number_pmesh_node
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
                   !rock particle diameter
                   d=node_val(di%heap%rock_d,node)
                   !Reynolds number
                   Re=u*rho*d/mu

                   if (di%MIM_options%Lima_mass_trans) then
                      !use the emperical correlation from Lima,2005 to calculate mass transfer coefficient 
                      di%MIM_options%mass_trans_coef(i)%ptr%val(node)=1.59*(Re**0.578)/3600.0
                   end if

                   if (di%MIM_options%Lima_immobile_sat) then
                      !use the emperical correlation from Lima,2005 to calculate MIM saturations
                      !where alpha=Sat_d/Sat_s, (1/alpha=0.173Re^-0.286 nfrom Lima, 2005)
                      alpha=5.78*(Re**0.286)
                      immobile_sat%val(node)=total_sat%val(node)/(1.0+alpha)
                   end if
                                   
                end do              
             end if 
             
             
             call set(mobile_sat, total_sat)
             call addto(mobile_sat, immobile_sat, scale=-1.0)
             call set(di%MIM_options%mobile_saturation(i)%ptr, mobile_sat)
          end if
          nullify(immobile_sat)
          nullify(total_sat)
          call zero(mobile_sat)
        end do
        
        call deallocate(mobile_sat)

  end subroutine darcy_trans_solve_MIM_saturations_and_mass_transfer_coefficient
      


  subroutine darcy_trans_MIM_prog_sfield_allocate_rhs_lhs(di,p,f,shared_rhs, shared_lhs,isub)

         type(darcy_impes_type), intent(inout) :: di
         integer, intent(in) :: p,f
         type(scalar_field), intent(inout) :: shared_rhs, shared_lhs
         integer,optional,intent(in) :: isub
         
         type(scalar_field) :: MIM_src,MIM_src_s,leach_im_src,theta_s,theta_d,theta_s_old, theta_d_old, src_cv_mass, old_cd, old_cs
         real :: dt
         real :: isub_e, isub_s
         integer :: i

         if (di%lcsub%have_leach_subcycle) then
            dt=di%lcsub%sub_dt
         else
            dt=di%dt
         end if

         if (di%lcsub%have_leach_subcycle) then
            isub_e=real(isub)/real(di%lcsub%number_subcycle)
            isub_s=real(isub-1)/real(di%lcsub%number_subcycle) 
         else
            isub_s=0
            isub_e=1.0
         end if

         call allocate(MIM_src, di%pressure_mesh)
         call allocate(MIM_src_s, di%pressure_mesh)
         call allocate(theta_s, di%pressure_mesh)
         call allocate(theta_s_old, di%pressure_mesh)
         call allocate(theta_d, di%pressure_mesh)
         call allocate(theta_d_old, di%pressure_mesh)
         call allocate(src_cv_mass, di%pressure_mesh)
         call allocate(old_cd, di%pressure_mesh)
         call allocate(old_cs, di%pressure_mesh)

         if (di%lcsub%have_leach_subcycle) then
             !mobile liquid hold up
             call set(theta_d,di%lcsub%sub_lht(p))
             call scale(theta_d,isub_e)
             call addto(theta_d,di%lcsub%old_sub_lht(p),scale=(1-isub_e))
             !immobile liquid hold up  
             call set(theta_s,di%lcsub%sub_lht_im(p))
             call scale(theta_s,isub_e)
             call addto(theta_s,di%lcsub%old_sub_lht_im(p),scale=(1-isub_e))
             
             !!for the old liquid hold of previous subcycle
             !old mobile liquid hold up
             call set(theta_d_old,di%lcsub%sub_lht(p))
             call scale(theta_d_old,isub_s)
             call addto(theta_d_old,di%lcsub%old_sub_lht(p),scale=(1-isub_s))
             !old immobile liquid hold up  
             call set(theta_s_old,di%lcsub%sub_lht_im(p))
             call scale(theta_s_old,isub_s)
             call addto(theta_s_old,di%lcsub%old_sub_lht_im(p),scale=(1-isub_s))

             call set(old_cd,di%lcsub%iterated_sfield(f))
             call set(old_cs,di%lcsub%iterated_imsfield(f))
             
          else
             !mobile liquid hold up
             call set(theta_d,di%MIM_options%mobile_saturation(p)%ptr)
             call scale(theta_d, di%porosity_pmesh)
            
             !immobile liquid hold up
             call set(theta_s,di%MIM_options%immobile_saturation(p)%ptr)
             call scale(theta_s,di%porosity_pmesh) 
             
             !!for the old liquid hold up of previous timestep
             !mobile liquid hold up
             call set(theta_d_old,di%MIM_options%old_mobile_saturation(p)%ptr)
             call scale(theta_d_old,di%old_porosity_pmesh)
             !immobile liquid hold up
             call set(theta_s_old, di%MIM_options%old_immobile_saturation(p)%ptr)
             call scale(theta_s_old, di%old_porosity_pmesh)

             call set(old_cd,di%generic_prog_sfield(f)%old_sfield)
             call set(old_cs,di%generic_prog_sfield(f)%MIM%immobile_sfield%old_sfield)

         end if
         

         !-------the component addto the lhs matrix
         !calculate (theta_s_new+alpha*dt)
         call set(shared_lhs, theta_s)

         call addto(shared_lhs, di%MIM_options%mass_trans_coef(p)%ptr, scale=dt)
         
         !-------the component added to the rhs 
         !calculate (theta_s_old*Cs_old)
         call set(shared_rhs,theta_s_old)

         call scale(shared_rhs,old_cs)
        
         if (.not.(di%generic_prog_sfield(f)%MIM%chem%if_src_linear)) then
            
            !the MIM without any leaching chemistry souces
            !which is shared terms with the model with leaching chemistry souces but without source linearization
            !-------the component add tp the rhs
            !!!!calculate shared part rhs '(theta_s_old*Cs_old)/(theta_s_new+alpha*dt)'
            call invert(shared_lhs) ! now shared lhs '1/(theta_s_new+alpha*dt)'
            
            call scale(shared_rhs,shared_lhs)
        
            !***add to rhs
            call set(MIM_src, shared_rhs)
            !'(alpha*theta_s_old*Cs_old)/(theta_s_new+alpha*dt)'
            call scale(MIM_src, di%MIM_options%mass_trans_coef(p)%ptr)
            call compute_cv_mass(di%positions,src_cv_mass, MIM_src)
            call addto (di%rhs, src_cv_mass)
                    
            if(di%generic_prog_sfield(f)%MIM%chem%have_chem) then
               call allocate(leach_im_src, di%pressure_mesh)
               !change to the unit in per volume of heap
               call set(leach_im_src,di%generic_prog_sfield(f)%MIM%chem%im_src%sfield)
               call scale(leach_im_src,theta_s)
                
              !!! calculate the shared rhs, which is
              !!! (theta_s_old*Cs_old)/(theta_s_new+alpha*dt)+(dt*S_immobile)/(theta_s_new+alpha*dt)
               call set(MIM_src, shared_lhs) !shared lhs is '1/(theta_s_new+alpha*dt)'
               call scale(MIM_src,leach_im_src)
               call scale(MIM_src,dt) !(dt*S_immobile)/(theta_s_new+alpha*dt)
               call addto(shared_rhs,MIM_src)

               !add the chemistry source part into rhs
               !***(alpha*dt*S_immobile)/(theta_s_new+alpha*dt)
               call scale(MIM_src,di%MIM_options%mass_trans_coef(p)%ptr)

               !change the mobile source into the unit of per volume of heap
               call set(leach_im_src,di%generic_prog_sfield(f)%MIM%chem%mo_src%sfield)
               call scale(leach_im_src,theta_d)
               
               !***(alpha*dt*S_immobile)/(theta_s_new+alpha*dt)+S_mobile
               !final rhs is ('theta_d_old*Cd_old/dt)*(alpha*theta_s_old*Cs_old)/(theta_s_new+alpha*dt)+(alpha*dt*S_immobile)/(theta_s_new+alpha*dt)+S_mobile'
               call addto(MIM_src,leach_im_src)
               
               call compute_cv_mass(di%positions,src_cv_mass,MIM_src)

               call addto(di%rhs, src_cv_mass)

               call deallocate(leach_im_src)
               
            end if

            
            !-------the component addto the lhs matrix
            !!!!calculate (alpha*dt)/(theta_s_new+alpha*dt), which is the shared part           
            call scale(shared_lhs,dt)
            call scale(shared_lhs,di%MIM_options%mass_trans_coef(p)%ptr)
            !***add to lhs
            call set(MIM_src, shared_lhs)
            !calculate (alpha^2*dt)/(theta_s_new+alpha*dt)
            call scale(MIM_src,di%MIM_options%mass_trans_coef(p)%ptr)
            !calculate '(alpha^2*dt)/(theta_s_new+alpha*dt)-alpha'
            call addto(MIM_src,di%MIM_options%mass_trans_coef(p)%ptr, scale=-1.0)
            call compute_cv_mass(di%positions,src_cv_mass, MIM_src)
            call addto(di%lhs,src_cv_mass, scale=-1.0)
            
         else
            !----------the MIM with leaching chemistry souces and source linearization---------

            !---------add to rhs
            call allocate(leach_im_src, di%pressure_mesh)
            !the negative part of leaching immobile source, change to the unit in per volume of heap
            call set(leach_im_src,di%generic_prog_sfield(f)%MIM%chem%im_src%n_src)
            call scale(leach_im_src,theta_s)

            node_loop3: do i=1,di%number_pmesh_node
               if (shared_rhs%val(i)<=1.0D-15) then
                  MIM_src%val(i)=0.0

               else
                  MIM_src%val(i)=theta_s%val(i)*dt*leach_im_src%val(i)/shared_rhs%val(i)!(theta_s*dt*S_NegativeImmobie)/(theta_s_old*Cs_old)
               end if
             
            end do node_loop3
            
            call addto(MIM_src, shared_lhs, scale=-1.0)!-((theta_s_new+alpha*dt))+(theta_s*dt*S_NegativeImmobie)/(theta_s_old*Cs_old)
            call scale(MIM_src, -1.0) !(theta_s_new+alpha*dt)-(theta_s*dt*S_NegativeImmobie)/(theta_s_old*Cs_old)
            call invert(MIM_src) !1/((theta_s_new+alpha*dt)-(theta_s*dt*S_NegativeImmobie)/(theta_s_old*Cs_old))

            !the positive part of leaching immobile source, change to the unit in per volume of heap
            call set(leach_im_src,di%generic_prog_sfield(f)%MIM%chem%im_src%p_src)
            call scale(leach_im_src,theta_s)

            call scale(leach_im_src,dt)
            
            !!!the shared rhs is '(theta_s_old*Cs_old+dt*S_PositiveImmobile)/(theta_s_new+alpha*dt-theta_s*dt*S_NegativeImmobie/theta_s_old*Cs_old)'
            call addto(shared_rhs, leach_im_src)
            call scale(shared_rhs,MIM_src)

            !'alpha*(theta_s_old*Cs_old+dt*S_PositiveImmobile)/(theta_s_new+alpha*dt-theta_s*dt*S_NegativeImmobie/theta_s_old*Cs_old)'
            call set(MIM_src_s,shared_rhs)
            call scale(MIM_src_s,di%MIM_options%mass_trans_coef(p)%ptr)

            !the positive part of leaching mobile source, change to the unit in per volume of heap
            call set(leach_im_src,di%generic_prog_sfield(f)%MIM%chem%mo_src%p_src)
            call scale(leach_im_src,theta_d)

            !**add to rsh, the total parts which are added to rhs are
            !'alpha*(theta_s_old*Cs_old+dt*S_PositiveImmobile)/(theta_s_new+alpha*dt-theta_s*dt*S_NegativeImmobie/theta_s_old*Cs_old)  + S_PositiveMobile)'
            call addto(MIM_src_s, leach_im_src)
            
            call compute_cv_mass(di%positions,src_cv_mass,MIM_src_s)
            call addto(di%rhs, src_cv_mass)

            !--------------------add to the lhs matrix
            !!!the shared lhs is
            !!!alpha*dt/(theta_s_new+alpha*dt-theta_s*dt*S_NegativeImmobie/theta_s_old*Cs_old)
            call set(shared_lhs,MIM_src)
            call scale(shared_lhs,dt)
            call scale(shared_lhs,di%MIM_options%mass_trans_coef(p)%ptr)
            
            !the negative part of leaching mobile source, change to the unit in per volume of heap
            call set(leach_im_src,di%generic_prog_sfield(f)%MIM%chem%mo_src%n_src)
            call scale(leach_im_src,theta_d)            
            call scale(leach_im_src,theta_d) !S_NegativeMobile*theta_d
            
            call set(MIM_src, old_cd)
            call scale(MIM_src,theta_d_old)
            
            node_loop4: do i=1,di%number_pmesh_node
               if (MIM_src%val(i)<=1.0D-15) then
                  MIM_src%val(i)=0.0

               else
                  MIM_src%val(i)= leach_im_src%val(i)/MIM_src%val(i)  !S_NegativeMobile*theta_d/theta_d_old*Cd_old
               end if
             
            end do node_loop4
         
            call set(MIM_src_s, shared_lhs)
            call scale(MIM_src_s,di%MIM_options%mass_trans_coef(p)%ptr) !alpha^2*dt/(theta_s_new+alpha*dt-theta_s*dt*S_NegativeImmobie/theta_s_old*Cs_old)

            !****add to lhs matrix, the total part is
            !-(S_NegativeMobile*theta_d/theta_d_old*Cd_old+alpha^2*dt/(theta_s_new+alpha*dt-theta_s*dt*S_NegativeImmobie/theta_s_old*Cs_old)-alpha)
            call addto(MIM_src, MIM_src_s)
            call addto(MIM_src, di%MIM_options%mass_trans_coef(p)%ptr, scale=-1.0)
            call compute_cv_mass(di%positions,src_cv_mass, MIM_src)
            call addto(di%lhs, src_cv_mass, scale=-1.0)
            
            call deallocate(leach_im_src)
         end if
        
        call deallocate(MIM_src)
        call deallocate(MIM_src_s)
        call deallocate(theta_s)
        call deallocate(theta_s_old)
        call deallocate(theta_d)
        call deallocate(theta_d_old)
        call deallocate(src_cv_mass)
        call deallocate(old_cd)
        call deallocate(old_cs)
  end subroutine darcy_trans_MIM_prog_sfield_allocate_rhs_lhs

  
  !solve the immobile  sfield 
  subroutine darcy_trans_assemble_and_solve_immobile_sfield(di, p, f, shared_rhs, shared_lhs)
       
       type(darcy_impes_type), intent(inout) :: di
       type(scalar_field), intent(in) :: shared_rhs, shared_lhs
       integer, intent(in) :: p
       integer, intent(in) :: f

        
       ! calculate '(alpha*dt*C_d)/(theta_s+alpha*dt)'
       !or '(alpha*dt*C_d)/(theta_s_new+alpha*dt-theta_s*dt*S_NegativeImmobie/theta_s_old*Cs_old)' with src linearization
       call set(di%generic_prog_sfield(f)%MIM%immobile_sfield%sfield, shared_lhs)
     
       call scale(di%generic_prog_sfield(f)%MIM%immobile_sfield%sfield, di%generic_prog_sfield(f)%sfield)
      
       !calculate '(alpha*dt*C_d+theta_s_old*C_s_old+dt*S_immobile)/(theta_s+alpha*dt)'
       !or '(alpha*dt*C_d+theta_s_old*C_s_old+dt*S_PositiveImmobile)/(theta_s_new+alpha*dt-theta_s*dt*S_NegativeImmobie/theta_s_old*Cs_old)' with src linearization
       call addto(di%generic_prog_sfield(f)%MIM%immobile_sfield%sfield,shared_rhs)
       
       ewrite(1,*) 'Finished assemble and solve immobile prog sfield ',trim(di%generic_prog_sfield(f)%MIM%immobile_sfield%sfield%name),' of phase ',p

   end subroutine darcy_trans_assemble_and_solve_immobile_sfield

  subroutine leaching_MIM_calculate_fields_and_ratio(di,p,f)
     !calculate the average concentration based on the mobile immboile concentration
     !calculate the weighting constant ratio of mobile and immobile concentration
     
      type(darcy_impes_type), intent(inout) :: di
      integer, intent(in) ::p,f

      !local variables
      type(scalar_field) :: sfield,dfield,thetab,tb,tfield
      integer :: i

      !!!----------calculate the average concentration--------
      call allocate(sfield,di%pressure_mesh)
      call allocate(dfield,di%pressure_mesh)
      call allocate(tb,di%pressure_mesh)
      call allocate(tfield,di%pressure_mesh)
      call zero(sfield)
      call zero(dfield)
      call zero(tb)
      call zero(tfield) 

      !the total liquid hold up
      call set(tb,di%porosity_pmesh)
      call scale(tb,di%saturation(p)%ptr)
      call invert(tb,tfield)
      
      !the immobile contribution
      call set(sfield,di%generic_prog_sfield(f)%MIM%immobile_sfield%sfield)
      call scale(sfield,di%porosity_pmesh)
      call scale(sfield,di%MIM_options%immobile_saturation(p)%ptr)
      
      !the mobile contribution
      call set(dfield,di%generic_prog_sfield(f)%sfield)
      call scale(dfield,di%porosity_pmesh)
      call scale(dfield,di%MIM_options%mobile_saturation(p)%ptr)
      
      !the average concentration
      call set(di%generic_prog_sfield(f)%MIM%C_a ,sfield)
      call addto(di%generic_prog_sfield(f)%MIM%C_a,dfield)
      call scale(di%generic_prog_sfield(f)%MIM%C_a,tfield)

     !------------calculate the weighting constant ratio--------

      do i=1, di%number_pmesh_node
         if (di%generic_prog_sfield(f)%MIM%C_a%val(i)<=1.0D-16) then
            di%generic_prog_sfield(f)%MIM%Fs%val(i)=0.5
            di%generic_prog_sfield(f)%MIM%Fd%val(i)=0.5
         else
             tfield%val(i)=1.0/(tb%val(i)*di%generic_prog_sfield(f)%MIM%C_a%val(i))
             di%generic_prog_sfield(f)%MIM%Fd%val(i)=tfield%val(i)*dfield%val(i)

             di%generic_prog_sfield(f)%MIM%Fs%val(i)=tfield%val(i)*sfield%val(i)
          end if
       end do
                
      call deallocate(sfield)
      call deallocate(dfield)
      call deallocate(tb)
      call deallocate(tfield)
    end subroutine leaching_MIM_calculate_fields_and_ratio
    
end module darcy_transport_model
