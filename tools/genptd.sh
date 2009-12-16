#!/bin/sh

rm -rf ptd; mkdir ptd
find . -mindepth 2 -maxdepth 2 -name '*.F'   -exec cp {} ptd \;
find . -mindepth 2 -maxdepth 2 -name '*.F90' -exec cp {} ptd \;

FFLAGS=$(grep ^FFLAGS Makefile | sed 's/FFLAGS  =  //')

cd ptd
for file in *
do
  gfortran $FFLAGS -I. -c -fdump-parse-tree $file > $file.ptd
  rm $file
done


