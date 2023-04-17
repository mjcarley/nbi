## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)

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
nbi_solve_laplace_CPPFLAGS = \
	-I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir) \
	`./petsc-config include`
nbi_solve_laplace_DEPS = $(top_builddir)/src/libnbi.la
nbi_solve_laplace_LDFLAGS = $(top_builddir)/src/.libs/libnbi.a \
	`./petsc-config lib`

nbi_solve_helmholtz_CPPFLAGS = \
	-I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir) \
	`./petsc-config include`
nbi_solve_helmholtz_DEPS = $(top_builddir)/src/libnbi.la
nbi_solve_helmholtz_LDFLAGS = $(top_builddir)/src/.libs/libnbi.a \
	`./petsc-config lib`
else
nbi_solve_laplace_CPPFLAGS = \
	-I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)
nbi_solve_laplace_DEPS = $(top_builddir)/src/libnbi.la
nbi_solve_laplace_LDFLAGS = $(top_builddir)/src/.libs/libnbi.a

nbi_solve_helmholtz_CPPFLAGS = \
	-I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)
nbi_solve_helmholtz_DEPS = $(top_builddir)/src/libnbi.la
nbi_solve_helmholtz_LDFLAGS = $(top_builddir)/src/.libs/libnbi.a
endif

