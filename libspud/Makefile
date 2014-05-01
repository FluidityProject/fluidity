#    Copyright (C) 2007 Imperial College London and others.
#
#    Please see the AUTHORS file in the main source directory for a full list
#    of copyright holders.
#
#    Applied Modelling and Computation Group
#    Department of Earth Science and Engineering
#    Imperial College London
#
#    David.Ham@Imperial.ac.uk
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation,
#    version 2.1 of the License.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#    USA

SHELL = /bin/sh

FC      = gfortran
FCFLAGS = -Iinclude  -DUSING_GFORTRAN=1 -DNDEBUG=1 -O3 -ffast-math -frecord-marker=4 -I/usr/include -ffree-line-length-none -ffixed-line-length-none

CXX     = g++
CXXFLAGS= -Iinclude  -DUSING_GFORTRAN=1 -DNDEBUG=1 -O3 

MAKE    = make
AR      = ar
ARFLAGS = cr

LIB = libspud.la
LIBS = $(shell echo   -lpthread -lm -lstdc++   -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lquadmath  -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lquadmath | sed 's/-L /-L/g')

HEADERS = include/spud.h include/spud include/spud_enums.h include/tinystr.h include/tinyxml.h
MODS = spud.mod
OBJS = spud.lo spud_interfaces.lo tinystr.lo tinyxml.lo tinyxmlerror.lo tinyxmlparser.lo fspud.lo

VPATH = src/

.SUFFIXES: .f90 .F90 .c .cpp .lo .la

.f90.lo:
	./libtool --mode=compile --tag=FC $(FC) $(FCFLAGS) -c $<
.cpp.lo:
	./libtool --mode=compile --tag=CXX $(CXX) $(CXXFLAGS) -c $<

default: libspud.la build-diamond

libspud.la: $(OBJS)
	./libtool --mode=link --tag=FC $(FC) $(FCFLAGS) -o $(LIB) $(OBJS) $(LIBS) -rpath /data/dunhui/ROM_Non-intrusive-ann/libspud/../lib 
	if test -f .libs/libspud.a; then cp .libs/libspud.a .; fi
	(if test -f .libs/libspud.so; then cp .libs/libspud.so ./libspud.so.0; ln -s libspud.so.0 libspud.so; fi) || true

build-diamond:
	cd diamond; python setup.py build; cd ..

test: unittest

unittest: libspud.la
	@cd src/tests; $(MAKE)

.PHONY:doc

doc: 
	@cd doc; $(MAKE) spud_manual.pdf

install: install-libspud install-spudtools install-diamond install-pyspud install-dxdiff

install-libspud: libspud.la
	/usr/bin/install -c -d $(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/../lib
	/usr/bin/install -c -d $(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/../include
	/usr/bin/install -c libspud.a $(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/../lib
	(if test -f libspud.so.0; then /usr/bin/install -c libspud.so.0 $(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/../lib; fi) || true
	(if test -f libspud.so; then cd $(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/../lib/; ln -s libspud.so.0 libspud.so; cd -; fi) || true
	/usr/bin/install -c $(MODS) $(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/../include
	/usr/bin/install -c $(HEADERS) $(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/../include

install-spudtools: 
	/usr/bin/install -c -d $(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/../share/spud
	/usr/bin/install -c -d $(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/../bin
	/usr/bin/install -c -m755 bin/spud-preprocess $(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/../bin
	/usr/bin/install -c -m755 bin/spud-set $(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/../bin
	/usr/bin/install -c -m644 schema/spud_base.rnc $(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/../share/spud
	/usr/bin/install -c -m644 schema/spud_base.rng $(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/../share/spud

install-diamond:
	cd diamond;	python setup.py install --prefix=$(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/..; cd ..

install-pyspud:
ifeq ($(origin BUILDING_DEBIAN),undefined)
	cd python; python setup.py install --prefix=$(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/..; cd ..
else
	cd python; for python in $(shell pyversions -r); do $$python setup.py install --prefix=$(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/.. --install-layout=deb; done; cd ..
endif

install-dxdiff:
ifeq ($(origin DESTDIR),undefined)
	cd dxdiff; python setup.py install --prefix=$(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/..; cd ..
else
	cd dxdiff; for python in $(shell pyversions -r); do $$python setup.py install --prefix=$(DESTDIR)/data/dunhui/ROM_Non-intrusive-ann/libspud/.. --install-layout=deb; done; cd ..
endif

clean:
	@cd doc; $(MAKE) clean
	rm -f *.o libspud.a libspud.so* *.o *.la *.mod *.lo
	rm -rf .libs
	@cd src/tests; $(MAKE) clean
	@cd diamond; rm -rf build; cd ..
	@cd python; rm -rf build; cd ..
	@cd dxdiff; rm -rf build; cd ..


distclean: clean
	@cd src/tests; $(MAKE) distclean
	rm -fr config.log config.status libtool autom4te.cache Makefile
	rm -f Makefile src/tests/Makefile examples/Makefile
	rm -f bin/spud-preprocess
	rm -f diamond/setup.py diamond/diamond/plugins.py diamond/diamond/preprocess.py
