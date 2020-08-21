AM_CPPFLAGS += -I$(top_srcdir)/src/include
AM_LDFLAGS += -L${abs_top_builddir}/src/h5core/.libs

if ENABLE_DEBUG
AM_LDFLAGS += -static
endif

LDADD =
LDADD += -lH5hut

noinst_PROGRAMS =

if ENABLE_C
noinst_PROGRAMS += h5u_test h5b_test
endif

TESTS_ENVIRONMENT = env LD_LIBRARY_PATH=@HDF5_PREFIX@/lib:$(LD_LIBRARY_PATH)

EXTRA_PROGRAMS = h5u_test h5b_test

h5u_test_SOURCES =	\
		h5u_test.c \
		h5u_write.c \
		h5u_read.c \
		testframe.c \
		testframe.h \
		params.h

h5b_test_SOURCES =	\
		h5b_test.c \
		h5b_write.c \
		h5b_read.c \
		testframe.c \
		testframe.h \
		params.h

clean-local:
	${RM} -f *.h5 *~
