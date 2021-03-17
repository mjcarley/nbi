AC_INIT([nbi], 0.1)
AC_CONFIG_SRCDIR([src/nbi.h])
AC_CONFIG_MACRO_DIR([m4])

NBI_MAJOR_VERSION=0
NBI_MINOR_VERSION=1
NBI_MICRO_VERSION=0
NBI_INTERFACE_AGE=0
NBI_BINARY_AGE=0
NBI_VERSION=$NBI_MAJOR_VERSION.$NBI_MINOR_VERSION.$NBI_MICRO_VERSION

AC_SUBST(NBI_MAJOR_VERSION)
AC_SUBST(NBI_MINOR_VERSION)
AC_SUBST(NBI_MICRO_VERSION)
AC_SUBST(NBI_VERSION)

# libtool versioning
LT_RELEASE=$NBI_MAJOR_VERSION.$NBI_MINOR_VERSION
LT_CURRENT=`expr $NBI_MICRO_VERSION - $NBI_INTERFACE_AGE`
LT_REVISION=$NBI_INTERFACE_AGE
LT_AGE=`expr $NBI_BINARY_AGE - $NBI_INTERFACE_AGE`
AC_SUBST(LT_RELEASE)
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)

# For automake.
VERSION=$NBI_VERSION
PACKAGE=nbi

AC_SUBST(PACKAGE)
AC_SUBST(VERSION)

AM_INIT_AUTOMAKE(subdir-objects)

##$PACKAGE, $VERSION, no-define)

# Specify a configuration file
AM_CONFIG_HEADER(config.h)

AC_DEFINE_UNQUOTED(NBI_MAJOR_VERSION, $NBI_MAJOR_VERSION, [Major version])
AC_DEFINE_UNQUOTED(NBI_MINOR_VERSION, $NBI_MINOR_VERSION, [Minor version])
AC_DEFINE_UNQUOTED(NBI_MICRO_VERSION, $NBI_MICRO_VERSION, [Micro version])
AC_DEFINE_UNQUOTED(NBI_INTERFACE_AGE, $NBI_INTERFACE_AGE, [Interface age])
AC_DEFINE_UNQUOTED(NBI_BINARY_AGE, $NBI_BINARY_AGE, [Binary age])

dnl Initialize libtool
AM_PROG_LIBTOOL

dnl Initialize maintainer mode
#AM_MAINTAINER_MODE

AC_PROG_CC
AM_PROG_CC_C_O

if test x$GCC = xyes ; then
  CFLAGS="$CFLAGS -Wall -Werror-implicit-function-declaration -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations"
fi

AC_PROG_AWK
AC_ISC_POSIX
AC_HEADER_STDC


AC_CHECK_PROG(SQT, sqt-config, yes, no)
if test x$SQT = xyes ; then
   AC_DEFINE(HAVE_SQT, 1, [singular quadrature library available])
   sqt_cflags=`sqt-config --cflags`
   sqt_libs=`sqt-config --libs`
else
  echo \#################################################################
  echo \#
  echo \# Singular Quadrature on Triangles library not installed or
  echo \# not properly installed. You can get it from
  echo \#
  echo \# https://github.com/mjcarley/sqt
  echo \#
  echo \#################################################################
  exit
fi

AC_CHECK_PROG(WBFMM, wbfmm-config, yes, no)
if test x$WBFMM = xyes ; then
   AC_DEFINE(HAVE_WBFMM, 1, [wide band FMM library available])
   wbfmm_cflags=`wbfmm-config --cflags`
   wbfmm_libs=`wbfmm-config --libs`
else
  echo \#################################################################
  echo \#
  echo \# Wide Band Fast Multipole Method library not installed or
  echo \# not properly installed. You can get it from
  echo \#
  echo \# https://github.com/mjcarley/wbfmm
  echo \#
  echo \#################################################################
  exit
fi

AC_CHECK_PROG(BLASWRAP, blaswrap-config, yes, no)
if test x$BLASWRAP = xyes ; then
   AC_DEFINE(HAVE_BLASWRAP, 1, [BLAS wrapper header available])
   blaswrap_cflags=`blaswrap-config --cflags`
   blaswrap_libs=`blaswrap-config --libs`
else
  echo \#################################################################
  echo \#
  echo \# BLAS wrapper library not installed or not properly installed.
  echo \# You can get it from
  echo \#
  echo \# https://github.com/mjcarley/blaswrap
  echo \#
  echo \#################################################################
  exit
fi

AX_GCC_X86_CPU_SUPPORTS(avx, hasavx=yes) 
AX_GCC_X86_CPU_SUPPORTS(avx2, hasavx2=yes) 
AX_GCC_X86_CPU_SUPPORTS(fma, hasfma=yes) 

if test x$hasavx = xyes; then
   SIMD_FLAGS="$SIMD_FLAGS -mavx"
fi
if test x$hasfma = xyes; then
   SIMD_FLAGS="$SIMD_FLAGS -mfma"
fi

AC_CHECK_LIB(m, cos)
AC_CONFIG_MACRO_DIR([m4])

AC_OPENMP
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN
AC_F77_WRAPPERS
ACX_BLAS
ACX_LAPACK

blas_libs="$BLAS_LIBS"
lapack_libs="$LAPACK_LIBS"

#DX_DOXYGEN_FEATURE(ON)
##DX_DOT_FEATURE(ON)
#DX_HTML_FEATURE(ON)
##DX_CHM_FEATURE(OFF)
##DX_CHI_FEATURE(OFF)
##DX_MAN_FEATURE(OFF)
##DX_RTF_FEATURE(OFF)
##DX_XML_FEATURE(OFF)
##DX_PDF_FEATURE(OFF)
##DX_PS_FEATURE(OFF)
#DX_INIT_DOXYGEN(nbi, doc/nbi.dxy)

CFLAGS="$OPENMP_CFLAGS $SIMD_FLAGS $glib_cflags $wbfmm_cflags $sqt_cflags $blaswrap_cflags $CFLAGS "
LDFLAGS="$LDFLAGS $glib_libs $wbfmm_libs $sqt_libs"

AC_SUBST(blas_libs)
AC_SUBST(lapack_libs)
AC_SUBST(CFLAGS)
AC_SUBST(LDFLAGS)
#AC_SUBST(DX_RULES)

LIBS="$glib_libs $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

AC_CONFIG_FILES([
Makefile
src/Makefile
src/nbi-config
tools/Makefile])

#doc/Makefile])

#doc/html/Makefile])

AC_OUTPUT

