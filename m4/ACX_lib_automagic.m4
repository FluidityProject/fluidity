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
    found="yes"
  ],
  [
    AC_MSG_RESULT([no])
    if test "x$found" = "xno" ; then
      AC_CHECK_LIB(
        [$1],
        [$2],
        [found="yes"],
        [found="no"])
    fi
  ])

if test "x$found" = "xyes" ; then
  true
  $3
elif test "x$found" = "xno" ; then
  true
  $4
fi

AC_LANG_POP([C])

])dnl ACX_lib_automagic
