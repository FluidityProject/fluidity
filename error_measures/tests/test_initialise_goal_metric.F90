subroutine test_initialise_goal_metric

  use unittest_tools
  use goal_metric
  use populate_state_module
  use spud
  implicit none

  logical :: fail
  integer :: stat

  call add_option("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity", stat=stat)
  call add_option("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity/enstrophy_goal", stat=stat)
  call set_option("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity/enstrophy_goal/dependencies", &
       & "Velocity%1 Velocity%2 Velocity%3", stat=stat)
  call set_option("/mesh_adaptivity/hr_adaptivity/goal_based_adaptivity/relative_tolerance", 0.10, stat=stat)

  call initialise_goal_metric

  if (goal_rel_tolerance /= 0.1) then
    fail = .true.
  else
    fail = .false.
  end if
  call report_test("[initialise_goal_metric]", fail, .false., "Get the relative tolerance correctly")

  if (use_goal_metric .eqv. .false.) then
    fail = .true.
  else
    fail = .false.
  end if
  call report_test("[initialise_goal_metric]", fail, .false., "We should use goal-based adaptivity")

  if (goal_name /= "enstrophy_goal") then
    fail = .true.
  else
    fail = .false.
  end if
  call report_test("[initialise_goal_metric]", fail, .false., "We are using the enstrophy goal")

  if (size(goal_deps) /= 3) then
    fail = .true.
  else
    fail = .false.
  end if
  call report_test("[initialise_goal_metric]", fail, .false., "Enstrophy has 3 dependencies")

  if (trim(goal_deps(1)) /= "Velocity%1") then
    fail = .true.
  else
    fail = .false.
  end if
  call report_test("[initialise_goal_metric]", fail, .false., "The first dependency is Velocity%1.")

  if (trim(goal_deps(2)) /= "Velocity%2") then
    fail = .true.
  else
    fail = .false.
  end if
  call report_test("[initialise_goal_metric]", fail, .false., "The second dependency is Velocity%2.")

  if (trim(goal_deps(3)) /= "Velocity%3") then
    write(0,*) "trim(goal_deps(3)) == ", trim(goal_deps(3))
    fail = .true.
  else
    fail = .false.
  end if
  call report_test("[initialise_goal_metric]", fail, .false., "The third dependency is Velocity%3.")

end subroutine test_initialise_goal_metric
