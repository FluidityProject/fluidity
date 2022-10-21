#!/bin/bash
#
# Package script for Diamond.
#
# Thanks: http://stackoverflow.com/questions/1596945/building-osx-app-bundle

# Instructions:
#  bzr co lp:spud mac_bundle
#  cd mac_bundle
#  ./configure
#  cd diamond/macosx
#  ./mac_bundle
#  *make a coffee or tea*
#  Done
#


# Also fix $INSTALLDIR/MacOS/diamond in case this number changes
PYVER=2.7
APP=Diamond.app
INSTALLDIR=$APP/Contents
LIBDIR=$INSTALLDIR/lib
APPDIR=/Applications/$INSTALLDIR
LOCALDIR=/opt/gtk

# set up our virtual python environment
virtualenv --python=python$PYVER --no-site-packages $INSTALLDIR

cd ../
# install diamond in it
sed -i .bak 's/packaging=False/packaging=True/' setup.py
macosx/$INSTALLDIR/bin/python setup.py install

# install dxdiff
cd ../dxdiff;
../diamond/macosx/$INSTALLDIR/bin/python setup.py install

cd ../diamond/macosx;

mkdir $INSTALLDIR/MacOS
mkdir $INSTALLDIR/Resources

# Sort out the Resources
cp PkgInfo $INSTALLDIR/
cp Info.plist $INSTALLDIR/
cp pango_rc $INSTALLDIR/Resources/
cp diamond.icns $INSTALLDIR/Resources/
cp diamond $INSTALLDIR/MacOS

# Now we have to play silly buggers with some bits of the diamond file
# as the Mac app packages adds a command line argument, which we want to ignore
sed -i .bak 's/sys.argv\[1:\]/sys.argv\[2:\]/' $INSTALLDIR/bin/diamond



# Let's get the latest fluidity release schema
# NOTE: UPDATE URL AFTER A RELEASE
if [ ! -d fluidity ]; then
	bzr branch lp:fluidity/4.1 fluidity
else
	cd fluidity
	bzr up
	cd ../
fi
# Make the schemata description
rm -rf $INSTALLDIR/Resources/schemata/fluidity
mkdir -p $INSTALLDIR/Resources/schemata/fluidity
cat > $INSTALLDIR/Resources/schemata/flml << EOF
Fluidity Markup Language
/Applications/Diamond.app/Contents/Resources/schemata/fluidity/fluidity_options.rng
EOF
cp fluidity/schemas/*.rng $INSTALLDIR/Resources/schemata/fluidity/
cp ../../schema/*.rng $INSTALLDIR/Resources/schemata/fluidity/
# clean up
#rm -rf fluidity

# Do the above for any other schema we want to distribute
# Don't forget the spud-base!



# Let's get lxml installed
$INSTALLDIR/bin/easy_install --allow-hosts=lxml.de,*.python.org lxml

# Temp. solution - Just manually copy stuff we know we need
SITEPACKAGES=$LIBDIR/python$PYVER/site-packages

mkdir -p $SITEPACKAGES

# This locates pygtk.pyc. We want the source file
pygtk=`python -c "import pygtk; print pygtk.__file__[:-1]"`
oldsite=`dirname $pygtk`
gobject=`python -c "import gobject; print gobject.__file__[:-1]"`
glib=`python -c "import glib; print glib.__file__[:-1]"`

# Copy PyGtk and related libraries
cp $pygtk $SITEPACKAGES
cp -r `dirname $gobject` $SITEPACKAGES
cp -r `dirname $glib` $SITEPACKAGES
cp -r $oldsite/cairo $SITEPACKAGES
cp -r $oldsite/gtk-2.0 $SITEPACKAGES
cp $oldsite/pygtk.pth $SITEPACKAGES


# Modules, config, etc.
for dir in lib etc/pango etc/gtk-2.0 share/themes; do
  mkdir -p $INSTALLDIR/$dir
  cp -r $LOCALDIR/$dir/* $INSTALLDIR/$dir
done

# Resources, are processed on startup
for dir in etc/gtk-2.0 etc/pango lib/gdk-pixbuf-2.0/2.10.0; do
  mkdir -p $INSTALLDIR/Resources/$dir
  cp $LOCALDIR/$dir/* $INSTALLDIR/Resources/$dir
done

# Somehow files are writen with mode 444
find $INSTALLDIR -type f -exec chmod u+w {} \;

function log() {
  echo $* >&2
}

function resolve_deps() {
  local lib=$1
  local dep
  otool -L $lib | grep -e "^.$LOCALDIR/" |\
      while read dep _; do
    echo $dep
  done
}

function fix_paths() {
  local lib=$1
  log Fixing $lib
  for dep in `resolve_deps $lib`; do
    log Fixing `basename $lib`
    log "|  $dep"
    install_name_tool -change $dep @executable_path/../lib/`basename $dep` $lib
  done
}

binlibs=`find $INSTALLDIR -type f -name '*.so'`

for lib in $binlibs; do
  log Resolving $lib
  resolve_deps $lib
  fix_paths $lib
done | sort -u | while read lib; do
  log Copying $lib
  cp $lib $LIBDIR
  chmod u+w $LIBDIR/`basename $lib`
  fix_paths $LIBDIR/`basename $lib`
done

for lib in $dylibs; do
  log Resolving $lib
  resolve_deps $lib
  fix_paths $lib
done | sort -u | while read lib; do
  log Copying $lib
  cp $lib $LIBDIR
  chmod u+w $LIBDIR/`basename $lib`
  fix_paths $LIBDIR/`basename $lib`
done


# Fix config files
sed -i .bak 's#/opt/gtk/#'$APPDIR'/#' $INSTALLDIR/etc/pango/pango.modules 
sed -i .bak 's#/opt/gtk/#'$APPDIR'/#' $INSTALLDIR/lib/gdk-pixbuf-2.0/2.10.0/loaders.cache
sed -i .bak 's#/opt/gtk/#'$APPDIR'/#' $INSTALLDIR/Resources/etc/pango/pango.modules 
sed -i .bak 's#/opt/gtk/#'$APPDIR'/#' $INSTALLDIR/Resources/lib/gdk-pixbuf-2.0/2.10.0/loaders.cache



# Package!
VERSION=0.01
# we now need to fiddle with the Python run path on the diamond script
# COMMENT THESE OUT IF YOU WANT TO TEST YOUR APP WITHOUT INSTALLING
# EDIT AS REQUIRED
sed -i .bak 's|/Users/amcg/Software/spud/mac_port/diamond|/Applications|' $INSTALLDIR/bin/diamond
sed -i .bak 's|/Users/amcg/Software/spud/mac_port/diamond|/Applications|' $INSTALLDIR/MacOS/diamond

zip -rq Diamond-$VERSION-osx.zip $APP
hdiutil create -srcfolder $APP Diamond-$VERSION.dmg


# To test, temporarily move your /opt/gtk folder
# Move the Diamond.app folder to /Applications and go
# open -a Console 
# helps to debug