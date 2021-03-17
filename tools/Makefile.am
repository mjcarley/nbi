## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)

DEPS = $(top_builddir)/src/libnbi.la

bin_PROGRAMS = \
	nbi-solve \
	nbi-assemble

nbi_solve_LDADD = \
	$(top_builddir)/src/.libs/libnbi.a 

nbi_assemble_LDADD = \
	$(top_builddir)/src/.libs/libnbi.a 
