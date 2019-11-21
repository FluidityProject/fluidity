TESTPYTHONPATH=$(shell python3 -m fluidity.state_types 2>&1 2>/dev/null; echo $$?)
TESTPYTHONLOCAL=$(shell PYTHONPATH=${PYTHONPATH}:${PWD}/../../python python3 -m fluidity.state_types 2>&1 2>/dev/null; echo $$?)

TESTPYTHON=${TESTPYTHONPATH}${TESTPYTHONLOCAL}

TESTBINFLUIDITY=$(shell which fluidity >/dev/null ; echo $$?)
TESTBINFLREDECOMP=$(shell which flredecomp >/dev/null ; echo $$?)

TESTBINFLUIDITYLOCAL=$(shell test -x ../../bin/fluidity >/dev/null ; echo $$?)
TESTBINFLREDECOMPLOCAL=$(shell test -x ../../bin/flredecomp >/dev/null ; echo $$?)

BINPREFIX=$(shell test -x ../../bin/fluidity >/dev/null && echo ${PWD}/../../bin/)

TESTBIN=${TESTBINFLUIDITY}${TESTBINFLREDECOMP}
TESTBINLOCAL=${TESTBINFLUIDITYLOCAL}${TESTBINFLREDECOMPLOCAL}

TESTFLOPTIONSLOCAL=$(shell test -r ../../schemas/fluidity_options.rng >/dev/null ; echo $$? )
FLOPTIONSPREFIX=$(shell test -r ../../schemas/fluidity_options.rng >/dev/null && echo ${PWD}/../../schemas/)

envcheck:
# Do we have both fluidity and flredecomp on the existing path?
ifneq (${TESTBIN},00)
 # Do we only have fluidity (and not flredecomp) on the existing path?
 ifeq (${TESTBIN},01)
	@echo "*** ERROR ***\nThe 'fluidity' binary is available on your system, but the 'flredecomp' binary\nisn't. This indicates that you have a partial install; perhaps you copied the\n'fluidity' binary into a system directory by hand? If you installed from\nsource, please rebuild Fluidity using the '--prefix' argument to specify an\ninstall directory, and run 'make install'. If Fluidity was installed for youon\nthe system, please pass this message on to your system administrator and ask\nthem to resolve the problem."
	@exit 1
 else
  # fluidity isn't on the path; do we have both fluidity and flredecomp locally?
  ifneq (${TESTBINLOCAL},00)
   # Do we only have fluidity (and not flredecomp) locally?
   ifeq (${TESTBINLOCAL},01)
	@echo "*** ERROR ***\nYou appear to be running within a partially-built Fluidity source tree. The\n'fluidity' binary is present, but the 'flredecomp' binary was not found. Please\nensure that 'make fltools' has been run in this source tree before you proceed\nto run examples."
	@exit 1
   else
    # No fluidity binary on the system path or locally
	@echo "*** ERROR ***\nThe Fluidity binaries needed to run this example can't be found.\nTo fix this either:\n\n   (a) install the Fluidity binary package if you are on Ubuntu\n   (b) compile Fluidity from source and install it in a system-accessable\n         location\n   (c) set your PATH environment variable to point to a built version of\n         Fluidity.\n\nRefer to the Fluidity manual for more instructions on any of the above,\nor contact the fluidity@imperial.ac.uk mailing list.\n"
	@exit 1
   endif 
  endif
 endif
 # Do we have both system AND local builds available?
 ifeq (${TESTBINLOCAL},00)
	@echo "*** WARNING ***\nYou are running in a built Fluidity source tree, with the fluidity binary also\non your default path. This example will run with the version as built in the\nlocal source tree.\n\n"
 endif
endif
# Did we entirely fail to find fluidity python?
ifeq (${TESTPYTHON},11)
	@echo "*** ERROR ***\nFluidity python support, required by this example, was not found on your\nsystem. Please either (a) install the Fluidity binary package if you are on\nUbuntu, or (b) get hold of a Fluidity source tree and set your PYTHONPATH to\nrefer to this, as described in the Fluidity manual."
	@exit 1
endif
# Do we have the fluidity options locally?
ifneq (${TESTFLOPTIONSLOCAL},0)
	@echo "*** ERROR ***\nA Fluidity options Relax NG file, required by this example, was not found on your\nsystem. Please either (a) install the Fluidity binary package if you are on\nUbuntu, or (b) get hold of a Fluidity source tree and run this example there."
	@exit 1
endif
 

.phony: envcheck
