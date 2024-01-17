## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4

if SUBMODULES
AM_CPPFLAGS = -I@abs_top_builddir@ -I@abs_top_builddir@/src -I$(includedir) \
	-I@abs_top_builddir@/$(sqt_include) \
	-I@abs_top_builddir@/$(wbfmm_include) \
	-I@abs_top_builddir@/$(agg_include) \
	-DG_LOG_DOMAIN=\"NBI\"
AM_LDFLAGS = \
	-L@abs_top_builddir@/$(sqt_submodule) $(sqt_libs) \
	-L@abs_top_builddir@/$(wbfmm_submodule) $(wbfmm_libs) \
	-L@abs_top_builddir@/$(agg_submodule) $(agg_libs)
else
AM_CPPFLAGS = -I@abs_top_builddir@ -I@abs_top_builddir@/src -I$(includedir) \
		-DG_LOG_DOMAIN=\"NBI\"
AM_LDFLAGS = $(sqt_libs) $(wbfmm_libs) $(agg_libs)
endif

DEPS = $(top_builddir)/src/libnbi.la
LDADD = $(top_builddir)/src/.libs/libnbi.a

bin_PROGRAMS = \
	nbi-assemble-laplace \
	nbi-solve-laplace \
	nbi-field-laplace \
	nbi-assemble-helmholtz \
	nbi-solve-helmholtz \
	nbi-field-helmholtz \
	nbi-surface \
	nbi-process

dist = \
	petsc-makefile \
	petsc-config

nbi_solve_laplace_SOURCES = nbi-solve-laplace.c petsc-laplace.c
nbi_solve_helmholtz_SOURCES = nbi-solve-helmholtz.c petsc-helmholtz.c

if PETSC
nbi_solve_laplace_CPPFLAGS = $(AM_CPPFLAGS) `./petsc-config include`
nbi_solve_laplace_LDFLAGS  = $(AM_LDFLAGS) `./petsc-config lib`

nbi_solve_helmholtz_CPPFLAGS = $(AM_CPPFLAGS) `./petsc-config include`
nbi_solve_helmholtz_LDFLAGS  = $(AM_LDFLAGS) `./petsc-config lib`
else
nbi_solve_laplace_CPPFLAGS = $(AM_CPPFLAGS)
nbi_solve_laplace_LDFLAGS  = $(AM_LDFLAGS) `./petsc-config lib`

nbi_solve_helmholtz_CPPFLAGS = $(AM_CPPFLAGS)
nbi_solve_helmholtz_LDFLAGS  = $(AM_LDFLAGS) 

endif

