#include "fdebug.h"

module convective_metric
!!< Reimplement CJC's
!!< convective error metric for Lucy's OODC problem.
!!< If the dynamical system is unstable at a node
!!< (measured by the gradient of temperature being unstable),
!!< replace the XX,YY components of the error metric with the ZZ
!!< component.

  use fields
  use field_derivatives
  use parallel_fields
  !use ieee_arithmetic
  use state_module
  use vtk_interfaces, only: vtk_write_fields
  implicit none

  logical :: use_convective_metric
  real :: T_scale
  integer :: convection_metric_option

  contains

  subroutine initialise_convective_metric
    implicit none

    integer :: io1, have_flag
    CHARACTER(LEN=4096) :: data_file

    data_file = ' '
    data_file(1:12) = 'mixadapt.dat'

    inquire(file=data_file,exist=use_convective_metric)

    if (use_convective_metric) then
      open(unit=2502, file=data_file, status='old', &
              iostat=io1)

      if(io1.ne.0) then
        use_convective_metric = .false.
      else
        read(unit=2502,fmt='(I9)') convection_metric_option
        read(unit=2502,fmt='(f10.3)') T_scale
        close(unit=2502)
        if(convection_metric_option.eq.0) then
          use_convective_metric = .false.
        else
          use_convective_metric = .true.
        end if
      end if
    end if

  end subroutine initialise_convective_metric

  subroutine form_convective_metric(state, error_metric)
    type(state_type), intent(in) :: state
    type(tensor_field), intent(inout) :: error_metric
    !locals
    type(tensor_field), pointer :: viscosity,diffusivity
    type(scalar_field), pointer :: temp_field
    type(scalar_field) :: diffusivity_zz, viscosity_zz
    type(scalar_field), dimension(1) :: grad_temp
    type(tensor_field) :: temp_hess
    type(vector_field), pointer :: positions

    integer :: i
#ifdef METRIC_DEBUG
    integer, save :: adaptcnt = 0
    character(len=20) :: buf
#endif
    real :: Ray_Crit = 1400.0, Ray
    real :: convection_metric
    
    ! these need to be fixed to new options: !!!!!
    real:: gravty
    real, dimension(1):: dengam

    ewrite(1,*) 'Subroutine form_convective_metric'
    
    FLExit("Convective metric needs to be fixed to use with new options.")

    temp_field => extract_scalar_field(state,'Temperature')
    positions => extract_vector_field(state,'Coordinate')
    diffusivity => extract_tensor_field(state,'TemperatureDiffusivity')
    viscosity => extract_tensor_field(state,'Viscosity')
    diffusivity_zz  = extract_scalar_field(diffusivity,3,3)
    viscosity_zz  = extract_scalar_field(viscosity,3,3)

    call allocate(grad_temp(1), temp_field%mesh, 'Temperature Z gradient')
    call allocate(temp_hess, temp_field%mesh)

    call differentiate_field(temp_field, positions, &
         & (/.false., .false., .true./), grad_temp)
    call compute_hessian(temp_field, positions, temp_hess)

    do i=1,temp_field%mesh%nodes
       ewrite(3,*) 'size grad_temp', size(grad_temp)
       ewrite(3,*) 'size grad_temp(1)%val', size(grad_temp(1)%val)
       if(grad_temp(1)%val(i)<0.0) then

          ewrite(3,*) 'grad temp', grad_temp(1)%val(i)
          ewrite(3,*) 'T_scale', T_scale
          ewrite(3,*) 'diffusivity', diffusivity_zz%val(i) 
          ewrite(3,*) 'viscosity', viscosity_zz%val(i)

          select case(convection_metric_option)
          case (1)
             convection_metric = (T_scale**2)*(dengam(1)* &
                  & gravty*abs(temp_hess%val(3,3,i))/(  &
                  & diffusivity_zz%val(i)* &
                  & viscosity_zz%val(i)))**(1.0/6.0)
          case (2)
             convection_metric = T_scale*abs(temp_hess%val(3,3,i))
          case default
             FLAbort('bad convection metric option')
          end select
          ewrite(3,*) 'convection_metric = ',convection_metric
          error_metric%val(1,1,i) = max(convection_metric, &
               error_metric%val(1,1,i))
          error_metric%val(2,2,i) = max(convection_metric, &
               error_metric%val(2,2,i))

!!$          Ray = CalculateRayleighNumber(grad_temp(1)%val(i), &
!!$               temp_hess%val(3,3,i),T_scale, &
!!$               diffusivity_zz%val(i),viscosity_zz%val(i))
!!$          if(.not.ieee_is_nan(Ray)) then
!!$             if(Ray>Ray_Crit) then
!!$                ewrite(3,*) 'Rayleigh , min lengthscale', Ray, Minimum_Lengthscale(Ray)
!!$                error_metric%val(1,1,i) = max(error_metric%val(1,1,i), &
!!$                     1.0/(Minimum_Lengthscale(Ray)**2))
!!$                ewrite(3,*) 'node #, metric', i, error_metric%val(1,1,i)
 !            end if
          !       end if
       end if
    end do

#ifdef METRIC_DEBUG
    write(buf, '(i0)') adaptcnt
    call vtk_write_fields(trim("convective_metric_") // trim(buf), adaptcnt, positions, positions%mesh, &
                          sfields=(/temp_field, grad_temp(1)/), tfields=(/error_metric/))
    adaptcnt = adaptcnt + 1
#endif

    call deallocate(grad_temp(1))
    call deallocate(temp_hess)

    ewrite(1,*) 'END Subroutine form_convective_metric'

  end subroutine form_convective_metric

!!$  function CalculateRayleighNumber(temp_z,temp_zz,t_scale, &
!!$       muzz,kappazz) result(Ray)
!!$    real, intent(in) :: temp_z,temp_zz,t_scale,muzz,kappazz
!!$    real :: Ray
!!$    !Rayleigh number = (g dengam/mu kappa) * DT * L^3
!!$    !                = (g dengam/mu kappa) * DT/DZ * L^4
!!$    !locals
!!$    real :: L
!!$
!!$    L = t_scale/sqrt(temp_zz)
!!$    Ray = -gravty*dengam(1)/(muzz*kappazz) * temp_z * L**4
!!$
!!$  end function CalculateRayleighNumber
!!$
!!$  function Minimum_Lengthscale(Ray) result(ml)
!!$    !output minimum lengthscale as obtained from Zoe Robert's data
!!$    real, intent(in) :: Ray
!!$    real :: ml
!!$    real, parameter :: gamm = 0.0001, A = 0.9262
!!$
!!$    ml = A*exp(-gamm*Ray)
!!$
!  end function Minimum_Lengthscale

end module convective_metric
