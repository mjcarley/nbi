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
	nbi-private.h \
	laplace.c \
	solvers.c \
	tinyexpr.c \
	expressions.c

include_HEADERS = \
	nbi.h \
	tinyexpr.h

CLEANFILES = $(BUILT_SOURCES)
