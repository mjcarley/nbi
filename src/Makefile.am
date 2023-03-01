AM_CPPFLAGS = -I$(top_srcdir) -I$(includedir) -DG_LOG_DOMAIN=\"NBI\"

bin_SCRIPTS=nbi-config

BUILT_SOURCES= \
	nbi-config

nbi-config: nbi-config.in

lib_LTLIBRARIES = libnbi.la

libnbi_la_LDFLAGS = -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic

libnbi_la_SOURCES = \
	nbi.c \
	geometries.c \
	laplace.c \
	helmholtz.c \
	tinyexpr.c \
	expressions.c \
	agg.c \
	export.c \
	functions.c \
	assemble.c \
	nbi-private.h \
	tinyexpr.h

## only want PETSc flags to affect code which uses PETSc
## (confined to solvers.c)
libnbi_la_LIBADD = libsolvers.la
noinst_LTLIBRARIES = libsolvers.la
libsolvers_la_SOURCES = solvers.c
if PETSC
include ${PETSC_DIR}/lib/petsc/conf/variables
libsolvers_la_CPPFLAGS = ${PETSC_CC_INCLUDES} ${PETSC_FLAGS} $(AM_CPPFLAGS)
endif

include_HEADERS = \
	nbi.h

##	tinyexpr.h

CLEANFILES = $(BUILT_SOURCES)
