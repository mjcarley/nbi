## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4

DEPS = $(top_builddir)/src/libnbi.la
LDADD = $(top_builddir)/src/.libs/libnbi.a

if SUBMODULES
AM_CPPFLAGS = -I@abs_top_builddir@ -I@abs_top_builddir@/src -I$(includedir) \
	-I@abs_top_builddir@/$(sqt_include) \
	-I@abs_top_builddir@/$(wbfmm_include) \
	-I@abs_top_builddir@/$(agg_include) \
	-DG_LOG_DOMAIN=\"NBI\"
AM_LDFLAGS = \
	-L@abs_top_builddir@/$(sqt_submodule) \
	-L@abs_top_builddir@/$(wbfmm_submodule) \
	-L@abs_top_builddir@/$(agg_submodule)
else
AM_CPPFLAGS = -I@abs_top_builddir@ -I@abs_top_builddir@/src -I$(includedir) \
		-DG_LOG_DOMAIN=\"NBI\"
endif

bin_PROGRAMS = \
	nbi-assemble-laplace \
	nbi-solve-laplace \
	nbi-field-laplace \
	nbi-assemble-helmholtz \
	nbi-solve-helmholtz \
	nbi-field-helmholtz \
	nbi-surface \
	nbi-process

nbi_solve_laplace_SOURCES = nbi-solve-laplace.c petsc-laplace.c
nbi_solve_helmholtz_SOURCES = nbi-solve-helmholtz.c petsc-helmholtz.c
