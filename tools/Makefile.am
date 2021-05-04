## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)

DEPS = $(top_builddir)/src/libnbi.la

bin_PROGRAMS = \
	nbi-solve-laplace \
	nbi-assemble-laplace \
	nbi-surface

nbi_solve_laplace_LDADD = \
	$(top_builddir)/src/.libs/libnbi.a 

nbi_assemble_laplace_LDADD = \
	$(top_builddir)/src/.libs/libnbi.a 

nbi_surface_LDADD = \
	$(top_builddir)/src/.libs/libnbi.a 
