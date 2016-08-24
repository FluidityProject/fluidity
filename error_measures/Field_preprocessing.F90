#include "fdebug.h"

module field_preprocessing_module

  use spud
  use global_parameters
  use fields
  use smoothing_module
  use field_options

  implicit none

  private
  public :: preprocess_field

  contains

    subroutine preprocess_field(field_in, positions, field_out)
      !! Do whatever actions are specified in the flml.
      !! Currently, only a Helmholtz smoother is available.

      type(scalar_field), intent(inout) :: field_in
      type(vector_field), intent(in) :: positions
      type(scalar_field), intent(out) :: field_out

      logical :: have_helmholtz_smoother
      integer :: stat

      have_helmholtz_smoother = have_option(trim(complete_field_path(trim(field_in%option_path), stat=stat)) &
           // "/adaptivity_options/preprocessing/helmholtz_smoother")

      if (.not. have_helmholtz_smoother) then
        field_out = field_in
        call incref(field_out)
        return
      else if (have_helmholtz_smoother) then
        call apply_helmholtz_smoother(field_in, positions, field_out)
      end if

    end subroutine preprocess_field

    subroutine apply_helmholtz_smoother(field_in, positions, field_out)
      type(scalar_field), intent(inout) :: field_in
      type(vector_field), intent(in) :: positions
      type(scalar_field), intent(out) :: field_out
      real :: alpha
      character(len=OPTION_PATH_LEN) :: path

      call get_option(trim(complete_field_path(trim(field_in%option_path))) &
               & // "/adaptivity_options/preprocessing/helmholtz_smoother/smoothing_scale_factor", alpha)

      call allocate(field_out, field_in%mesh, "Smoothed" // trim(field_in%name))
      call zero(field_out)

      path = trim(complete_field_path(trim(field_in%option_path))) // "/adaptivity_options/preprocessing/helmholtz_smoother"
      call smooth_scalar(field_in, positions, field_out, alpha, path=path)
    end subroutine apply_helmholtz_smoother

end module field_preprocessing_module
