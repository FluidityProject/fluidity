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
  use darcy_impes_leaching_types
  use vtk_cache_module, only:vtk_cache_finalise
  use fefields, only: compute_cv_mass

  use darcy_impes_assemble_type, only: darcy_impes_type
  use darcy_impes_leaching_types


  implicit none
  private
  
  public :: darcy_trans_MIM_assemble_and_solve_mobile_saturation, &
       leaching_MIM_calculate_fields_and_ratio, &
       darcy_trans_MIM_prog_sfield_allocate_rhs_lhs, &
       darcy_trans_assemble_and_solve_immobile_sfield, &
       darcy_trans_assemble_galerkin_projection_elemesh_to_pmesh, &
       initialize_MIM_model,&       
       finalize_MIM_model




contains

  subroutine initialize_MIM_model(di)
      type(darcy_impes_type), intent(inout) :: di

      integer :: p, stat, f_count, f
      character(len=OPTION_PATH_LEN) :: tmp_char_option, tmp_char_option1
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

                     if (have_option('/material_phase['//int2str(p-1)//']/scalar_field['//int2str(f-1)//']/prognostic/scalar_field::Immobile'))  then
                        if (.not.(di%MIM_options%have_MIM(p))) then
                        FLExit('The Mobileimmobile model of phase'//int2str(p)//'should be turned on')
                           
                        end if

                        di%generic_prog_sfield(f_count)%MIM%have_MIM_source = .true.
                        ! allocate the MIM sorce terms for the matrix to solve the prognostic sfield
                        call allocate(di%MIM_options%MIM_src, di%pressure_mesh)
                        call allocate(di%MIM_options%MIM_src_s, di%pressure_mesh)
                        tmp_char_option1 = 'Immobile'//trim(tmp_char_option)

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
                           di%generic_prog_sfield(f_count)%MIM%chem%mo_src => extract_scalar_field(di%state(p), &
                                trim(tmp_char_option1), stat)
                           
                           if (.not. stat==0) then
                              FLExit('failed to extract the field:'//tmp_char_option1)
                           end if  

                           tmp_char_option1 = trim(tmp_char_option)//'Immobile_chemical_src'
                           di%generic_prog_sfield(f_count)%MIM%chem%im_src => extract_scalar_field(di%state(p), &
                                trim(tmp_char_option1), stat)
                           
                           if (.not. stat==0) then
                              FLExit('failed to extract the field:'//tmp_char_option1)
                           end if

                        elseif (have_option('/material_phase['//int2str(p-1)//']/scalar_field['//int2str(f-1)//']/prognostic/leaching_temperature_sources')) then
                           if (.not. have_option('/material_phase['//int2str(p-1)//']/scalar_field['//int2str(f-1)//']/prognostic/leaching_temperature_source/Mobile_Immobile_Model')) then
                              FLExit('please turn on the Mobile_Immobile_Model under leaching_temperature_sources of field:'// tmp_char_option)
                           end if

                           di%generic_prog_sfield(f_count)%MIM%chem%have_chem = .true.
                           tmp_char_option1 = trim(tmp_char_option)//'Mobile_chemical_src'
                           di%generic_prog_sfield(f_count)%MIM%chem%mo_src => extract_scalar_field(di%state(p), &
                                trim(tmp_char_option1), stat)
                           
                           if (.not. stat==0) then
                              FLExit('failed to extract the field:'//tmp_char_option1)
                           end if  

                           tmp_char_option1 = trim(tmp_char_option)//'Immobile_chemical_src'
                           di%generic_prog_sfield(f_count)%MIM%chem%im_src => extract_scalar_field(di%state(p), &
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

 subroutine finalize_MIM_model(di)

   type(darcy_impes_type), intent(inout) :: di

   ! local variables
   integer :: p, f

   do p = 1,di%number_phase
      nullify(di%capilliary_pressure(p)%ptr)
      nullify(di%saturation(p)%ptr)
      nullify(di%old_saturation(p)%ptr)
      nullify(di%saturation_source(p)%ptr)
      nullify(di%relative_permeability(p)%ptr)
      nullify(di%old_relative_permeability(p)%ptr)
      nullify(di%viscosity(p)%ptr)
      nullify(di%darcy_velocity(p)%ptr)
      nullify(di%cfl(p)%ptr)
      nullify(di%mobility(p)%ptr)
      nullify(di%fractional_flow(p)%ptr)
      nullify(di%density(p)%ptr)
      nullify(di%old_density(p)%ptr)
      nullify(di%MIM_options%old_immobile_saturation(p)%ptr)
      nullify(di%MIM_options%immobile_saturation(p)%ptr)
      nullify(di%MIM_options%mobile_saturation(p)%ptr)
      nullify(di%MIM_options%old_mobile_saturation(p)%ptr)
      nullify(di%MIM_options%mass_trans_coef(p)%ptr)
      nullify(di%MIM_options%old_mass_trans_coef(p)%ptr)
      di%MIM_options%have_MIM(p)=.false.
   end do

   ! Deallocate the MIM model
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
      call deallocate(di%MIM_options%MIM_src)      
      call deallocate(di%MIM_options%MIM_src_s)
      nullify(di%generic_prog_sfield(f)%MIM%immobile_sfield%sfield)
      nullify(di%generic_prog_sfield(f)%MIM%immobile_sfield%old_sfield)
      nullify(di%generic_prog_sfield(f)%MIM%C_a)
      nullify(di%generic_prog_sfield(f)%MIM%old_C_a)
      nullify(di%generic_prog_sfield(f)%MIM%Fs)
      nullify(di%generic_prog_sfield(f)%MIM%Fd)
      di%generic_prog_sfield(f)%MIM%have_MIM_source  = .false.
      if (di%generic_prog_sfield(f)%MIM%chem%have_chem) then
         nullify(di%generic_prog_sfield(f)%MIM%chem%mo_src)
         nullify(di%generic_prog_sfield(f)%MIM%chem%im_src)
         di%generic_prog_sfield(f)%MIM%chem%have_chem=.false.
      end if
   end if
   
  
   
 end subroutine finalize_MIM_model
 

 !Slove the mobile saturation if MIM exist
 subroutine darcy_trans_MIM_assemble_and_solve_mobile_saturation(di)

        type(darcy_impes_type), intent(inout) :: di
        integer :: i
        type(scalar_field), pointer :: total_sat => null()  !total saturation 
        type(scalar_field), pointer :: immobile_sat  => null()  ! immobile saturation
        type(scalar_field)  :: mobile_sat   !mobile saturation

        call allocate(mobile_sat, di%pressure_mesh)
        
        do i= 2, di%number_phase

          total_sat      => di%saturation(i)%ptr
          immobile_sat   => di%MIM_options%immobile_saturation(i)%ptr

          if (di%MIM_options%have_MIM(i)) then
             ewrite(1, *) "calculate the mobile saturation of phase: ", i
             call set(mobile_sat, total_sat)
             call addto(mobile_sat, immobile_sat, scale=-1.0)
             call set(di%MIM_options%mobile_saturation(i)%ptr, mobile_sat)
          end if
          nullify(immobile_sat, total_sat)
          call zero(mobile_sat)
        end do
        
        call deallocate(mobile_sat)

  end subroutine darcy_trans_MIM_assemble_and_solve_mobile_saturation


  subroutine darcy_trans_MIM_prog_sfield_allocate_rhs_lhs(di,p,f,temp_MIM_src)

         type(darcy_impes_type), intent(inout) :: di
         integer, intent(in) :: p,f
         type(scalar_field), intent(inout) :: temp_MIM_src
         
         type(scalar_field) :: leach_im_src
        

         call zero(di%MIM_options%MIM_src)
         call zero(di%MIM_options%MIM_src_s)
                  
         !Addto the lhs matrix
         !calculate 1/(theta_s+alpha*dt)
         call set(di%MIM_options%MIM_src, di%MIM_options%immobile_saturation(p)%ptr)
         
         !check wether to scale with the porosity as a constant of a scalar field

         call scale(di%MIM_options%MIM_src, di%porosity_pmesh)

         call addto(di%MIM_options%MIM_src, di%MIM_options%mass_trans_coef(p)%ptr, scale=di%dt)
         call invert(di%MIM_options%MIM_src, temp_MIM_src) !temp_MIM_src is '1/(theta_s+alpha*dt)'

         call set(di%MIM_options%MIM_src, temp_MIM_src)
         call scale(di%MIM_options%MIM_src, di%MIM_options%mass_trans_coef(p)%ptr) 
         call scale(di%MIM_options%MIM_src, di%MIM_options%mass_trans_coef(p)%ptr)! repeated to scale with alpha**2
         call scale(di%MIM_options%MIM_src, di%dt)
         call addto(di%MIM_options%MIM_src, di%MIM_options%mass_trans_coef(p)%ptr, scale=-1.0)
           
         call compute_cv_mass(di%positions, di%MIM_options%MIM_src_s, di%MIM_options%MIM_src)

         call addto(di%lhs, di%MIM_options%MIM_src_s, scale=-1.0)

         !Addto the rhs matrix
         !before start to compute rhs, zero the field to be used
         call zero(di%MIM_options%MIM_src)

         !Add the source term with old immobile concentration 
         call set(di%MIM_options%MIM_src,temp_MIM_src)
         call scale(di%MIM_options%MIM_src, di%MIM_options%mass_trans_coef(p)%ptr)
         call scale(di%MIM_options%MIM_src, di%MIM_options%old_immobile_saturation(p)%ptr)
         call scale(di%MIM_options%MIM_src, di%generic_prog_sfield(f)%MIM%immobile_sfield%old_sfield)
         call scale(di%MIM_options%MIM_src, di%cv_mass_pressure_mesh_with_old_porosity) ! this has already inlcuded the cv pmesh
         call addto (di%rhs, di%MIM_options%MIM_src)

         !chemical leaching src
         !add the immobile src part to the ADE
         if (di%generic_prog_sfield(f)%MIM%chem%have_chem) then
           call allocate(leach_im_src, di%pressure_mesh)
           call zero(leach_im_src)

           call set(leach_im_src,di%generic_prog_sfield(f)%MIM%chem%im_src)
           call scale(leach_im_src,di%dt)
           call scale(leach_im_src,di%MIM_options%mass_trans_coef(p)%ptr)
           call scale(leach_im_src,temp_MIM_src)
           call addto(di%rhs,leach_im_src)

           call deallocate(leach_im_src)
        end if

  end subroutine darcy_trans_MIM_prog_sfield_allocate_rhs_lhs

  
  !solve the immobile  sfield 
  subroutine darcy_trans_assemble_and_solve_immobile_sfield(di, p, f, temp_MIM_src)
       
       type(darcy_impes_type), intent(inout) :: di
       type(scalar_field), intent(in) :: temp_MIM_src
       integer, intent(in) :: p
       integer, intent(in) :: f

       ! local variable
       type(scalar_field) :: temp_rhs
       
       call allocate(temp_rhs, di%pressure_mesh)
       
       ! calculate '(old_theta_s*old_C_s)/(theta_s+alpha*dt)'
       call set(temp_rhs, temp_MIM_src)
       call scale(temp_rhs, di%MIM_options%old_immobile_saturation(p)%ptr)

       call scale(temp_rhs, di%old_porosity_pmesh)

       call scale(temp_rhs, di%generic_prog_sfield(f)%MIM%immobile_sfield%old_sfield)

       call set(di%generic_prog_sfield(f)%MIM%immobile_sfield%sfield, temp_rhs)
       
       !calculate '(alpha*dt*C_d)/(theta_s+alpha*dt)'
       call set(temp_rhs, temp_MIM_src)
       call scale(temp_rhs,di%MIM_options%mass_trans_coef(p)%ptr)
       call scale(temp_rhs, di%dt)
       call scale(temp_rhs, di%generic_prog_sfield(f)%sfield) ! this use the C_d value at most recent time step n+1
       
       call addto(di%generic_prog_sfield(f)%MIM%immobile_sfield%sfield, temp_rhs)
       
       !chemical leaching src
       !add the immobile src '(dt*S)/(theta_s+alpha*dt)'
       if (di%generic_prog_sfield(f)%MIM%chem%have_chem) then
          call set(temp_rhs,di%generic_prog_sfield(f)%MIM%chem%im_src)
          call scale(temp_rhs,di%dt)
          call scale(temp_rhs,temp_MIM_src)
          call addto(di%generic_prog_sfield(f)%MIM%immobile_sfield%sfield,temp_rhs)
       end if

       call deallocate(temp_rhs)

       ewrite(1,*) 'Finished assemble and solve immobile prog sfield ',trim(di%generic_prog_sfield(f)%MIM%immobile_sfield%sfield%name),' of phase ',p

   end subroutine darcy_trans_assemble_and_solve_immobile_sfield


   subroutine darcy_trans_assemble_galerkin_projection_elemesh_to_pmesh(field, projected_field, positions, ele)

        type(scalar_field), intent(inout) :: field
        type(scalar_field), intent(in) :: projected_field
        type(vector_field), intent(in) :: positions
        integer, intent(in) :: ele
        type(element_type), pointer :: field_shape, proj_field_shape
        real, dimension(ele_loc(field, ele)) :: little_rhs
        real, dimension(ele_loc(field, ele), ele_loc(field, ele)) :: little_mass
        real, dimension(ele_loc(field, ele), ele_loc(projected_field, ele)) :: little_mba
        real, dimension(ele_loc(field, ele), ele_loc(projected_field, ele)) :: little_mba_int
        real, dimension(ele_ngi(field, ele)) :: detwei
        real, dimension(ele_loc(projected_field, ele)) :: proj_field_val 
        
        integer :: i, j, k


          field_shape => ele_shape(field, ele)
          proj_field_shape => ele_shape(projected_field, ele)

          call transform_to_physical(positions, ele, detwei=detwei)

          little_mass = shape_shape(field_shape, field_shape, detwei)

          ! And compute the product of the basis functions
          little_mba = 0
          do i=1,ele_ngi(field, ele)
           forall(j=1:ele_loc(field, ele), k=1:ele_loc(projected_field, ele))
             little_mba_int(j, k) = field_shape%n(j, i) * proj_field_shape%n(k, i)
           end forall
           little_mba = little_mba + little_mba_int * detwei(i)
          end do

          proj_field_val = ele_val(projected_field, ele)
          little_rhs = matmul(little_mba, proj_field_val)

          call solve(little_mass, little_rhs)
          call set(field, ele_nodes(field, ele), little_rhs)
 

  end subroutine darcy_trans_assemble_galerkin_projection_elemesh_to_pmesh


      
  subroutine leaching_MIM_calculate_fields_and_ratio(di,p,f)
     !calculate the average concentration based on the mobile immboile concentration
     !calculate the weighting constant ratio of mobile and immobile concentration
     
      type(darcy_impes_type), intent(inout) :: di
      integer, intent(in) ::p,f

      !local variables
      type(scalar_field) :: sfield,dfield,thetab,tb,tfield

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

      !!!----------calculate the weighting constant ratio--------
      call scale(tb,di%generic_prog_sfield(f)%MIM%C_a )
      call invert(tb,tfield)
      
      call set(di%generic_prog_sfield(f)%MIM%Fd, tfield)
      call scale(di%generic_prog_sfield(f)%MIM%Fd,dfield)

      call set(di%generic_prog_sfield(f)%MIM%Fs, tfield)
      call scale(di%generic_prog_sfield(f)%MIM%Fs,sfield)
      
      
      call deallocate(sfield)
      call deallocate(dfield)
      call deallocate(tb)
      call deallocate(tfield)
    end subroutine leaching_MIM_calculate_fields_and_ratio
    
end module darcy_transport_model
