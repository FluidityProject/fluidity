TESTPYTHONPATH=$(shell python -m fluidity.state_types 2>&1 2>/dev/null; echo $$?)
TESTPYTHONLOCAL=$(shell PYTHONPATH=${PYTHONPATH}:${PWD}/../../python python -m fluidity.state_types 2>&1 2>/dev/null; echo $$?)

TESTPYTHON=${TESTPYTHONPATH}${TESTPYTHONLOCAL}

TESTBINFLUIDITY=$(shell which fluidity >/dev/null ; echo $$?)
TESTBINFLDECOMP=$(shell which fldecomp >/dev/null ; echo $$?)

TESTBINFLUIDITYLOCAL=$(shell test -x ../../bin/fluidity >/dev/null ; echo $$?)
TESTBINFLDECOMPLOCAL=$(shell test -x ../../bin/fldecomp >/dev/null ; echo $$?)

BINPREFIX=$(shell test -x ../../bin/fluidity >/dev/null && echo ${PWD}/../../bin/)

TESTBIN=${TESTBINFLUIDITY}${TESTBINFLDECOMP}
TESTBINLOCAL=${TESTBINFLUIDITYLOCAL}${TESTBINFLDECOMPLOCAL}

envcheck:
# Do we have both fluidity and fldecomp on the existing path?
ifneq (${TESTBIN},00)
 # Do we only have fluidity (and not fldecomp) on the existing path?
 ifeq (${TESTBIN},01)
	@echo "The 'fluidity' binary is available on your system, but the 'fldecomp' binary isn't. This indicates that you have a partial install; perhaps you copied the 'fluidity' binary into a system directory by hand? If you installed from source, please rebuild Fluidity using the '--prefix' argument to specify an install directory, and run 'make install'. If Fluidity was installed for you on the system, please pass this message on to your system administrator and ask them to resolve the problem."
	@exit 1
 else
  # fluidity isn't on the path; do we have both fluidity and fldecomp locally?
  ifneq (${TESTBINLOCAL},00)
   # Do we only have fluidity (and not fldecomp) locally?
   ifeq (${TESTBINLOCAL},01)
	@echo "You appear to be running within a partially-built Fluidity source tree. The 'fluidity' binary is present, but the 'fldecomp' binary was not found. Please ensure that 'make fltools' has been run in this source tree before you proceed to run examples."
	@exit 1
   else
    # No fluidity binary on the system path or locally
	@echo "The Fluidity binaries needed to run this example can't be found.  To fix this, either (a) install the Fluidity binary package if you are on Ubuntu, (b) compile Fluidity from source and install it in a system-accessable location, or (c) set your PATH environment variable to point to a built version of Fluidity. Refer to the Fluidity manual for more instructions on any of the above, or contact the fluidity@imperial.ac.uk mailing list."
	@exit 1
   endif 
  endif
 endif
 # Do we have both system AND local builds available?
 ifeq (${TESTBINLOCAL},00)
	@echo "CAUTION: You are running in a built Fluidity source tree, with the fluidity binary also on your default path. This example will run with the version as built in the local source tree.\n\n"
 endif
endif
# Did we entirely fail to find fluidity python?
ifeq (${TESTPYTHON},11)
	@echo "Fluidity python support, required by this example, was not found on your system. Please either (a) install the Fluidity binary package if you are on Ubuntu, or (b) get hold of a Fluidity source tree and set your PYTHONPATH to refer to this, as described in the Fluidity manual."
	@exit 1
endif

.phony: envcheck
