#include "fdebug.h"

module interpolation_metric
!!< This module implements the standard Hessian-based
!!< error metric, as described in
!!<  "Tetrahedral mesh optimisation and adaptivity for steady-state and transient
!!<   finite element calculations", Pain et. al,
!!<  Comput. Methods Apll. Mech Engrg. 190 (2001) 3771-3796

  use spud
  use fldebug
  use metric_tools
  use fields
  use parallel_fields
  use state_module
  use vtk_interfaces
  use merge_tensors
  use halos
  use field_derivatives
  use field_options
  use form_metric_field
  use edge_length_module
  use aspect_ratios_module
  use field_preprocessing_module
  use project_metric_to_surface_module
  
  implicit none

  logical :: use_interpolation_metric

  private
  public :: use_interpolation_metric, initialise_interpolation_metric,&
            form_interpolation_metric

  contains

  subroutine initialise_interpolation_metric
    !!< For now we always want to use this metric.
    use_interpolation_metric = .true.
  end subroutine initialise_interpolation_metric

  subroutine form_interpolation_metric(state, error_metric)
    type(state_type), intent(inout), dimension(:) :: state
    type(tensor_field), intent(inout) :: error_metric

    type(tensor_field) :: tmp_tensor
    integer :: i, j, k, l
    type(scalar_field) :: edgelen, aspect_ratios
    type(scalar_field) :: adweit_s
    type(vector_field) :: adweit_v
    type(tensor_field) :: adweit_t

    integer, save :: adaptcnt = 0
    character(len=20) :: buf
    type(vector_field), pointer :: positions

    type(scalar_field) :: field_s, preprocessed_field_s
    type(vector_field) :: field_v
    type(tensor_field) :: field_t

    type(state_type) :: fields_state, weights_state
    type(scalar_field), pointer, dimension(:) :: fields_list, weights_list
    integer :: dim
    logical :: debug_metric, align_metric_vertically

    ! First let's get the fields that are actually being used.
    ! If the error is set to the special value 0.0, it means ignore that field
    ! for the purposes of computing the metric for adaptivity.

    positions => extract_vector_field(state(1), "Coordinate")
    dim = error_metric%dim(1)
    debug_metric = have_option("/mesh_adaptivity/hr_adaptivity/debug/write_metric_stages")
    ! is this metric going to be collapsed in the vertical to do horizontal adaptivity with it?
    align_metric_vertically = have_option("/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity/vertically_align_metric")

    do i=1,size(state)
      do j=1,scalar_field_count(state(i))
        field_s = extract_scalar_field(state(i), j)
        if (aliased(field_s)) then
          cycle
        end if

        if (have_adapt_opt(trim(field_s%option_path), "/adaptivity_options") &
          & .and. .not. have_adapt_opt(trim(field_s%option_path), "/adaptivity_options/no_interpolation_measure") &
          & .and. .not. have_adapt_opt(trim(field_s%option_path), "/adaptivity_options/anisotropic_zienkiewicz_zhu")) then
          call preprocess_field(field_s, positions, preprocessed_field_s)
          call insert(fields_state, preprocessed_field_s, trim(state(i)%name)//"::"//trim(field_s%name))
          call deallocate(preprocessed_field_s)
          adweit_s = extract_scalar_field(state(i), trim(field_s%name) // "InterpolationErrorBound")
          call insert(weights_state, adweit_s, trim(state(i)%name)//"::"//trim(adweit_s%name))
        end if
      end do

      do j=1,vector_field_count(state(i))
        field_v = extract_vector_field(state(i), j)
        if (aliased(field_v)) then
          cycle
        end if

        if (have_adapt_opt(trim(field_v%option_path), "/adaptivity_options") &
          & .and. .not. have_adapt_opt(trim(field_s%option_path), "/adaptivity_options/no_interpolation_measure")) then
          adweit_v = extract_vector_field(state(i), trim(field_v%name) // "InterpolationErrorBound")
          do k=1,field_v%dim
            adweit_s = extract_scalar_field(adweit_v, k)
            if (minval(adweit_s) > 0.0) then
              field_s = extract_scalar_field(field_v, k)
              call preprocess_field(field_s, positions, preprocessed_field_s)
              call insert(fields_state, preprocessed_field_s, trim(state(i)%name)//"::"//trim(field_s%name))
              call insert(weights_state, adweit_s, trim(state(i)%name)//"::"//trim(adweit_s%name))
              call deallocate(preprocessed_field_s)
            end if
          end do
        end if
      end do

      do j=1,tensor_field_count(state(i))
        field_t = extract_tensor_field(state(i), j)
        if (aliased(field_t)) then
          cycle
        end if

        if (have_adapt_opt(trim(field_t%option_path), "/adaptivity_options") &
          & .and. .not. have_adapt_opt(trim(field_s%option_path), "/adaptivity_options/no_interpolation_measure")) then
          adweit_t = extract_tensor_field(state(i), trim(field_t%name) // "InterpolationErrorBound")
          do k=1,field_t%dim(1)
            do l=1,field_t%dim(2)
              adweit_s = extract_scalar_field(adweit_t, k, l)
              if (minval(adweit_s) > 0.0) then
                field_s = extract_scalar_field(field_t, k, l)
                call preprocess_field(field_s, positions, preprocessed_field_s)
                call insert(fields_state, preprocessed_field_s, trim(state(i)%name)//"::"//trim(field_s%name))
                call insert(weights_state, adweit_s, trim(state(i)%name)//"::"//trim(adweit_s%name))
                call deallocate(preprocessed_field_s)
              end if
            end do
          end do
        end if
      end do
    end do

    call collapse_state(fields_state, fields_list)
    call collapse_state(weights_state, weights_list)

    if (size(fields_list) == 0) then
      ewrite(1,*) "Interpolation metric doing nothing"
      return
    end if

    ! If we have more than one scalar field, we'll need to compute each individual
    ! metric and convolve the two. We need some memory to store the second metric.
    if (size(fields_list) > 1) then
      call allocate(tmp_tensor, error_metric%mesh, "Tmp")
    end if

    if (debug_metric) then
      call allocate(edgelen, positions%mesh, "Desired edge lengths")
      call allocate(aspect_ratios, positions%mesh, "Metric aspect ratio")
    endif

    ewrite(2,*) "++: Forming interpolation metric"

    call halo_update(fields_list(1))
    call compute_hessian(fields_list(1), positions, error_metric)
    
    if (align_metric_vertically) then
      ! state only used for "GravityDirection", so state(1) is fine
      call vertically_align_metric(state(1), error_metric)
    end if

    ewrite(2,*) "++: Hessian", 1, "formed"

    if (debug_metric) then
      call get_edge_lengths(error_metric, edgelen)
      call vtk_write_fields(trim("interpolation_metric_hessian_1"), adaptcnt, positions, positions%mesh, &
                            sfields=(/fields_list(1), edgelen/), tfields=(/error_metric/))
    endif

    call form_metric(error_metric, fields_list(1), weights_list(1), state(1))
    ewrite(2,*) "++: Metric", 1, "formed"

    if (debug_metric) then
      call get_edge_lengths(error_metric, edgelen)
      call get_aspect_ratios(error_metric, aspect_ratios)
      call vtk_write_fields(trim("interpolation_metric_metric_1"), adaptcnt, positions, positions%mesh, &
                            sfields=(/fields_list(1), edgelen, aspect_ratios/), tfields=(/error_metric/))
    endif

    ! Now loop through the rest, construct each individual metric and merge with the first.
    do i=2,size(fields_list)
      write(buf, '(i0)') i
      call halo_update(fields_list(i))
      call compute_hessian(fields_list(i), positions, tmp_tensor)
      ewrite(2,*) "++: Hessian", i, "formed"

      if (align_metric_vertically) then
        ! state only used for "GravityDirection", so state(1) is fine
        call vertically_align_metric(state(1), tmp_tensor)
      end if

      if (debug_metric) then
        call get_edge_lengths(tmp_tensor, edgelen)
        call vtk_write_fields(trim("interpolation_metric_hessian_") // trim(buf), adaptcnt, positions, positions%mesh, &
                              sfields=(/fields_list(i), edgelen/), tfields=(/tmp_tensor/))
      endif

      call form_metric(tmp_tensor, fields_list(i), weights_list(i), state(1))
      ewrite(2,*) "++: Metric", i, "formed"

      if (debug_metric) then
        call get_edge_lengths(tmp_tensor, edgelen)
        call get_aspect_ratios(tmp_tensor, aspect_ratios)
        call vtk_write_fields(trim("interpolation_metric_metric_") // trim(buf), adaptcnt, positions, positions%mesh, &
                              sfields=(/fields_list(i), edgelen, aspect_ratios/), tfields=(/tmp_tensor/))
      endif

      call merge_tensor_fields(error_metric, tmp_tensor)
      ewrite(2,*) "++: Metric", i, "merged"

      if (debug_metric) then
        call get_edge_lengths(error_metric, edgelen)
        call get_aspect_ratios(error_metric, aspect_ratios)
        call vtk_write_fields(trim("interpolation_metric_merge_") // trim(buf), adaptcnt, positions, positions%mesh, &
                              sfields=(/fields_list(i), edgelen, aspect_ratios/), tfields=(/error_metric/))
      endif
    end do

   call check_metric(error_metric)
   
   if (size(fields_list) > 1) then 
    call deallocate(tmp_tensor)
   end if

   deallocate(fields_list)
   deallocate(weights_list)
   call deallocate(fields_state)
   call deallocate(weights_state)

   if (debug_metric) then
     call get_edge_lengths(error_metric, edgelen)
     call get_aspect_ratios(error_metric, aspect_ratios)
     call vtk_write_fields(trim("interpolation_metric_final"), adaptcnt, positions, positions%mesh, &
                              sfields=(/edgelen, aspect_ratios/), tfields=(/error_metric/))
     call deallocate(edgelen)
     call deallocate(aspect_ratios)
   endif

   adaptcnt = adaptcnt + 1

  end subroutine form_interpolation_metric

end module interpolation_metric
