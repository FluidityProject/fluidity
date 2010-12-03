subroutine test_fluxes_reader_wrapper
  use unittest_tools
  implicit none

  real :: correct, value, correct_d, correct_t, values(2)
  logical :: fail

  ! Run the first GetScalar(..) test from FluxesReader test via the Fortran wrappers
  fail = .true.
  call fluxes_registerdatafile("../../tests/data/global_fluxes.nc")
  call fluxes_addfieldofinterest("_2t")
  call fluxes_setsimulationtimeunits("seconds since 1960-01-01 06:00:0.0")
  call fluxes_settimeseconds(0.0)

  correct = 299.237;
  call fluxes_getscalar("_2t",210.0,10.0, value)
  if(abs(value-correct) < 0.001) then
    fail = .false.
  else
    write(0,*) "Expected ", correct, ", got ", value
  end if  
  call report_test("[test_fluxes_reader_wrapper: GetScalar single point 1]", fail, .false., &
                     "Got incorrect value")

  ! Run the first GetScalars(..) test from FluxesReader test via the Fortran wrappers
  fail = .true.
  call fluxes_addfieldofinterest("_2d")
  call fluxes_settimeseconds(0.0)

  correct_t = 299.237;
  correct_d = 296.316;
  call fluxes_getscalars(210.0,10.0, values)
  if(abs(values(1)-correct_t) < 0.001) then
    fail = .false.
  else
    write(0,*) "t expected ", correct, ", got ", value
  end if  
  call report_test("[test_fluxes_reader_wrapper: GetScalars single point: t]", fail, .false., &
                     "Got incorrect value")
  if(abs(values(2)-correct_d) < 0.001) then
    fail = .false.
  else
    write(0,*) "d expected ", correct, ", got ", value
  end if 
  call report_test("[test_fluxes_reader_wrapper: GetScalars single point: d]", fail, .false., &
                     "Got incorrect value")

end subroutine test_fluxes_reader_wrapper
