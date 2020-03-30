# ACX_lib_automagic(library, function, [action-if-found], [action-if-not-found])
#
AC_DEFUN([ACX_lib_automagic], [

AC_LANG_PUSH([C])

AC_MSG_CHECKING([$1 automagic])
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([[
void $2();
                   ]],[[
$2();
                   ]])],
  [
    AC_MSG_RESULT([yes])
    $3
  ],
  [
    AC_MSG_RESULT([no])
    AC_CHECK_LIB(
      [$1],
      [$2],
      [$3
       LIBS="-l$1 $LIBS"],
      [$4])
  ])

AC_LANG_POP([C])

])
AC_DEFUN([ACX_lib_fautomagic], [

AC_LANG_PUSH([Fortran])
AC_FC_SRCEXT(f90)

AC_MSG_CHECKING([$1 automagic])
AC_LINK_IFELSE(
  [
   program test
   
   interface
      subroutine $2()
      end subroutine
   end interface	 

   call $2()

   end program test

  ],
  [
    AC_MSG_RESULT([yes])
    $3
  ],
  [
    AC_MSG_RESULT([no])
    AC_CHECK_LIB(
      [$1],
      [$2],
      [$3
       LIBS="-l$1 $LIBS"],
      [$4])
  ])

AC_FC_SRCEXT(f)
AC_LANG_POP([Fortran])

])dnl ACX_lib_automagic
