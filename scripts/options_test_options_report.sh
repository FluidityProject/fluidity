#!/bin/bash -x

FLUIDITY_DIR=$1

IS_FLUIDITY=`head -1 ${FLUIDITY_DIR}/README | grep FLUIDITY`

if [ "X${IS_FLUIDITY}" == "X" ] ; then
   echo "ERROR: Invalid fluidity base directory supplied as parameter"
   exit 0 ;
fi

# Find all the source files of interest
files=`find ${FLUIDITY_DIR} -name "*.[cF][9p][0p]"`

# Generate .gcov files
for i in $files ; do
    dirname=`dirname $i`
    basename=`basename $i`
    pushd $dirname
    gcov $basename
    popd
done

# Remove the $PWD/ from the paths.
path_pattern=`(echo $PWD/ | sed 's/\//\\\\\//g')`
tmpfile=`mktemp gcovXXXXXX`
find ${FLUIDITY_DIR} -name "*.gcov" | eval "sed 's/${path_pattern}//'" > $tmpfile

# Generate fluidity options test coverage report
${FLUIDITY_DIR}/scripts/gcovifs2.py $tmpfile > report.txt
rm $tmpfile

