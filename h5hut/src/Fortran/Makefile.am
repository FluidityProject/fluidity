
if ENABLE_FORTRAN
AM_CPPFLAGS += -I$(top_srcdir)/src/include

F90_FILES =						\
	$(top_srcdir)/src/include/H5hut.f90		\
	$(top_srcdir)/src/include/H5_const.f90		\
	$(top_srcdir)/src/include/H5_file.f90		\
	$(top_srcdir)/src/include/H5_model.f90		\
	$(top_srcdir)/src/include/H5_file_attribs.f90	\
	$(top_srcdir)/src/include/H5_step_attribs.f90	\
	$(top_srcdir)/src/include/H5_err.f90		\
	$(top_srcdir)/src/include/H5_log.f90		\
	$(top_srcdir)/src/include/H5Part_model.f90	\
	$(top_srcdir)/src/include/H5Part_io.f90		\
	$(top_srcdir)/src/include/H5Block_model.f90	\
	$(top_srcdir)/src/include/H5Block_attribs.f90	\
	$(top_srcdir)/src/include/H5Block_io.f90

EXTRA_HEADERS =

EXTRA_DIST =			\
	$(F90_FILES)		\
	h5_private.h

nodist_include_HEADERS = 	\
	$(top_srcdir)/src/include/H5hut.f

lib_LTLIBRARIES = libH5hutF.la

include_HEADERS =		\
	$(F90_FILES)		\
	$(top_srcdir)/src/include/H5hut.f

libH5hutF_la_SOURCES =		\
	H5.c			\
	H5_attribs.c		\
	H5_model.c		\
	H5Part_model.c		\
	H5Part_io.c		\
	H5Block_model.c		\
	H5Block_attribs.c	\
	H5Block_io.c

libH5hutF_la_DEPENDENCIES =	\
	$(top_srcdir)/src/include/H5hut.f

libH5hutF_la_LDFLAGS = -version-info 2:0:0 -rpath '$(libdir)'

$(top_srcdir)/src/include/H5hut.f: $(F90_FILES)
	awk '/INTEGER\*8 :: /{print "       " $$0}' $^ >$@
	awk '/PARAMETER /{print "       " $$0}' $^ >> $@
	awk '/INTEGER\*8 FUNCTION/{print "       " $$1 " " $$3}' $^ >>$@

all-local: $(top_srcdir)/src/include/H5hut.f

endif

clean: clean-am

clean-local:
	$(RM) -f *~
