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
	solvers.c \
	nbi-private.h \
	tinyexpr.h

include_HEADERS = \
	nbi.h

CLEANFILES = $(BUILT_SOURCES)
