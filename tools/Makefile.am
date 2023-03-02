## Process this file with automake to produce Makefile.in


AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)

DEPS = $(top_builddir)/src/libnbi.la
LDADD = $(top_builddir)/src/.libs/libnbi.a

bin_PROGRAMS = \
	nbi-solve-laplace \
	nbi-assemble-laplace \
	nbi-solve-helmholtz \
	nbi-assemble-helmholtz \
	nbi-surface \
	nbi-process

if PETSC
include ${PETSC_DIR}/lib/petsc/conf/variables
AM_CPPFLAGS += ${PETSC_CC_INCLUDES} ${PETSC_FLAGS}
LDADD += $(PETSC_LIB)
##nbi_solve_laplace_CPPFLAGS = \
##	${PETSC_CC_INCLUDES} ${PETSC_FLAGS} $(AM_CPPFLAGS)
##nbi_solve_laplace_LDADD =  $(LDADD) $(PETSC_LIB)
##nbi_solve_helmholtz_CPPFLAGS = \
##	${PETSC_CC_INCLUDES} ${PETSC_FLAGS} $(AM_CPPFLAGS)
endif
