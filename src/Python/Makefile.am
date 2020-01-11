if ENABLE_PYTHON

BUILT_SOURCES = H5hut_wrap.c
SWIG_SOURCES = H5hut.i numpy.i

pkgpython_PYTHON = H5hut.py __init__.py
pkgpyexec_LTLIBRARIES = _H5hut.la

_H5hut_la_SOURCES = H5hut_wrap.c $(SWIG_SOURCES)
_H5hut_la_LDFLAGS = -module
_H5hut_la_LIBADD = $(top_builddir)/src/h5core/libH5hut.la 
_H5hut_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS)
_H5hut_la_CPPFLAGS += -I$(top_srcdir)/src/include
_H5hut_la_CPPFLAGS += -I$(PYTHON_PREFIX)/include/python$(PYTHON_VERSION)/
_H5hut_la_CPPFLAGS += -I$(NUMPY_INCLUDE)

if ENABLE_PARALLEL

_H5hut_la_CPPFLAGS += -DPARALLEL_IO  -I$(MPI4PY_INCLUDE)

H5hut_wrap.c : $(SWIG_SOURCES)
	$(SWIG) $(SWIG_PYTHON_OPT) -python -I$(top_srcdir)/src/include -DPARALLEL_IO -I$(MPI4PY_INCLUDE) -o $@ $<

else

H5hut_wrap.c : $(SWIG_SOURCES)
	$(SWIG) $(SWIG_PYTHON_OPT) -python -I$(top_srcdir)/src/include -o $@ $<

endif

else
$(info Python is disabled)
endif

clean: clean-am

clean-local:
	$(RM) -f *~ H5hut_wrap.c H5hut.py *.pyc
