#!/bin/sh
# Run this to generate all the initial makefiles, etc.

ACLOCAL_FLAGS="-I m4 -I . $ACLOCAL_FLAGS"
LIBTOOLIZE_FLAGS="--force $LIBTOOLIZE_FLAGS"
AUTOMAKE_FLAGS="--add-missing --copy --foreign $AUTOMAKE_FLAGS"

LIBTOOLIZE=`which libtoolize`
if [ "$LIBTOOLIZE" = "" ]; then
        LIBTOOLIZE=`which glibtoolize`
fi
if [ "$LIBTOOLIZE" = "" ]; then
        echo "libtoolize not found" 1>&2
        exit 1
fi

echo "+ making misc files ..."
touch NEWS README AUTHORS ChangeLog
echo
$LIBTOOLIZE $LIBTOOLIZE_FLAGS || {
	echo "libtoolize failed"
	exit 1
}
echo
echo "+ running aclocal ..."
aclocal $ACLOCAL_FLAGS || {
	echo "aclocal failed - check that all needed development files are present on system"
	exit 1
}
echo
echo "+ running autoheader ... "
autoheader || {
	echo "autoheader failed"
	exit 1
}
echo
echo "+ running autoconf ... "
autoconf || {
	echo "autoconf failed"
	exit 1
}
echo
echo "+ running automake ... "
automake $AUTOMAKE_FLAGS || {
	echo "automake failed"
	exit 1
}

