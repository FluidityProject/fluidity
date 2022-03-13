AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS = -I$(top_srcdir)/src/include @AM_CPPFLAGS@

EXTRA_HEADERS =					\
	private/h5_const.h			\
	private/h5_file.h			\
	private/h5_hdf5.h			\
	private/h5_hsearch.h			\
	private/h5_init.h			\
	private/h5_lustre.h			\
	private/h5_maps.h			\
	private/h5_model.h			\
	private/h5_mpi.h			\
	private/h5_qsort.h			\
	private/h5_io.h				\
	private/h5_types.h			\
	private/h5_va_macros.h			\
	private/h5b_model.h			\
	private/h5b_io.h			\
	private/h5b_types.h			\
	private/h5t_access.h			\
	private/h5t_adjacencies.h		\
	private/h5t_core.h			\
	private/h5t_err.h			\
	private/h5t_map.h			\
	private/h5t_model.h			\
	private/h5t_octree.h			\
	private/h5t_io.h			\
	private/h5t_ref_elements.h		\
	private/h5t_retrieve.h			\
	private/h5t_store.h			\
	private/h5t_tags.h			\
	private/h5t_types.h			\
	private/h5u_io.h			\
	private/h5u_types.h			\
	private/h5_attribs.h			\
	private/h5_log.h			\
	private/h5_err.h			\
	private/h5_fcmp.h

# Extra files that I wish to include in the dist tar ball.
EXTRA_DIST = $(EXTRA_HEADERS)

# What to build... Will be determined by configure script.
lib_LTLIBRARIES = libH5hut.la

# Listing of sources
libH5hut_la_SOURCES =				\
	h5_attachments.c			\
	h5_attribs.c				\
	h5_err.c				\
	h5_log.c				\
	h5_file.c				\
	h5_model.c				\
	h5_syscall.c				\
	h5u_io.c		                \
	h5b_io.c 	                        \
	h5u_model.c				\
        h5b_model.c				\
        h5b_attribs.c				\
	private/h5_hdf5.c			\
	h5_init.c			\
	private/h5_hsearch.c			\
	private/h5_maps.c			\
	private/h5_fcmp.c			\
	private/h5_qsort.c			\
	private/h5_qsort_r.c			\
	private/h5_io.c				\
	private/h5_lustre.c

# grid stuff
libH5hut_la_SOURCES +=				\
	h5t_adjacencies.c			\
	h5t_map.c				\
	h5t_model.c				\
	h5t_octree.c				\
	h5t_io.c				\
	h5t_retrieve.c				\
	h5t_store.c				\
	h5t_tags.c				\
	private/h5t_core.c			\
	private/h5t_core_trim.c			\
	private/h5t_core_tetm.c			\
	private/h5t_access_tetm.c		\
	private/h5t_access_trim.c		\
	private/h5t_adjacencies_tetm.c		\
	private/h5t_adjacencies_trim.c		\
	private/h5t_model_tetm.c		\
	private/h5t_model_trim.c		\
	private/h5t_retrieve_tetm.c		\
	private/h5t_retrieve_trim.c		\
	private/h5t_io_trim.c			\
	private/h5t_io_tetm.c			\
	private/h5t_ref_elements.c		\
	private/h5t_store_tetm.c		\
	private/h5t_store_trim.c


libH5hut_la_DEPENDENCIES = $(EXTRA_HEADERS)

libH5hut_la_LDFLAGS = -version-info 2:0:0

install-exec-local:
	@$(INSTALL) -d $(DESTDIR)$(includedir)/h5core
	@$(INSTALL) -m644 $(top_srcdir)/src/include/h5core/*.h $(DESTDIR)$(includedir)/h5core/

clean-local:
	$(RM) -f *~
