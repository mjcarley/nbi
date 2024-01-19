AC_INIT([nbi],[1.0.0])
AC_CONFIG_SRCDIR([src/nbi.h])
AM_INIT_AUTOMAKE([subdir-objects])

NBI_MAJOR_VERSION=1
NBI_MINOR_VERSION=0
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

# Specify a configuration file
AC_CONFIG_HEADERS(config.h)

AC_DEFINE_UNQUOTED(NBI_MAJOR_VERSION, $NBI_MAJOR_VERSION, [Major version])
AC_DEFINE_UNQUOTED(NBI_MINOR_VERSION, $NBI_MINOR_VERSION, [Minor version])
AC_DEFINE_UNQUOTED(NBI_MICRO_VERSION, $NBI_MICRO_VERSION, [Micro version])
AC_DEFINE_UNQUOTED(NBI_INTERFACE_AGE, $NBI_INTERFACE_AGE, [Interface age])
AC_DEFINE_UNQUOTED(NBI_BINARY_AGE, $NBI_BINARY_AGE, [Binary age])

dnl Initialize libtool
LT_INIT

AC_PROG_CC
AM_PROG_CC_C_O

if test x$GCC = xyes ; then
  CFLAGS="$CFLAGS -Wall -Werror-implicit-function-declaration -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations"
fi

AC_PROG_AWK
AC_SEARCH_LIBS([strerror],[cposix])

submodules=no
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

if test -d "externals/sqt" ; then
  echo "SQT included as submodule."
  submodules=yes
  AC_CONFIG_SUBDIRS([externals/sqt])
  AC_DEFINE(HAVE_SQT, 1, [singular quadrature library available])

  sqt_submodule=externals/sqt/src/.libs
  sqt_libs=-lsqt
  sqt_include=externals/sqt/src
  ##ln -s ../externals/sqt/src/sqt.h src/sqt.h 
else
  echo "SQT not included as submodule."
  sqt_submodule=
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
fi

if test -d "externals/wbfmm" ; then
  echo "WBFMM included as submodule."
  submodules=yes
  AC_CONFIG_SUBDIRS([externals/wbfmm])
  AC_DEFINE(HAVE_WBFMM, 1, [wide band FMM library available])

  wbfmm_submodule=externals/wbfmm/src/.libs
  wbfmm_libs=-lwbfmm
  wbfmm_include=externals/wbfmm/src
  ##wbfmm_cflags=-I../${wbfmm_include}
  ##ln -s ../externals/wbfmm/src/wbfmm.h src/wbfmm.h 
else
  echo "WBFMM not included as submodule."
  wbfmm_submodule=
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
fi

if test -d "externals/agg" ; then
  echo "AGG included as submodule."
  submodules=yes
  AC_CONFIG_SUBDIRS([externals/agg])
  AC_DEFINE(HAVE_AGG, 1, [aerodynamic geometry library available])

  agg_submodule=externals/agg/src/.libs
  agg_libs=-lagg
  agg_include=externals/agg/src
  ##agg_cflags=-I../${agg_include}
  ##ln -s ../externals/agg/src/agg.h src/agg.h 
else
  AC_CHECK_PROG(AGG, agg-config, yes, no)
  if test x$AGG = xyes ; then
     AC_DEFINE(HAVE_AGG, 1, [aerodynamic geometry library available])
     agg_cflags=`agg-config --cflags`
     agg_libs=`agg-config --libs`
  else
    echo \#################################################################
    echo \#
    echo \# Aerodynamic Geometry Generator library not installed or
    echo \# not properly installed. You can get it from
    echo \#
    echo \# https://github.com/mjcarley/agg
    echo \#
    echo \#################################################################
    exit
  fi
fi

AC_CHECK_HEADERS(gmshc.h) 
AC_CHECK_LIB(gmsh, gmshInitialize, [], [
  echo \#################################################################
  echo \#
  echo \# GMSH API library not installed or not in library path
  echo \# You can get it from
  echo \#
  echo \# https://www.gmsh.info/
  echo \#
  echo \#################################################################
  exit		   
])

if test "${PETSC_DIR+set}" = set; then
  cd tools
  TMPFLAGS=`./petsc-config include`
  cc $TMPFLAGS petsc-version.c
  petsc_version=`./a.out`
  rm a.out
  cd ..
  petsc_major=`echo ${petsc_version} | cut -d . -f 1`
  petsc_minor=`echo ${petsc_version} | cut -d . -f 2`
  if [[ ${petsc_major} -lt 3 -o ${petsc_minor} -lt 17 ]]; then
  PETSC=no
  echo \#################################################################
  echo \#
  echo \# Environment variable PETSC_DIR=${PETSC_DIR}
  echo \# PETSc version ${petsc_version} detected
  echo \# PETSc version greater than or equal to 3.17 required
  echo \# PETSc solvers will NOT be available
  echo \#
  echo \#################################################################
  else
  AC_DEFINE(HAVE_PETSC, 1, [PETSC solvers available])
  PETSC=yes
  echo \#################################################################
  echo \#
  echo \# Environment variable PETSC_DIR=${PETSC_DIR}
  echo \# PETSc version ${petsc_version} detected
  echo \# PETSc solvers will be available
  echo \#
  echo \#################################################################
  fi  
else
  PETSC=no
  echo \#################################################################
  echo \#
  echo \# Environment variable PETSC_DIR not set: PETSc solvers will not
  echo \# be available
  echo \#
  echo \# If PETSc is installed on your system, check that PETSC_DIR
  echo \# is set correctly and run configure again.
  echo \#
  echo \# If PETSc is not installed, you will be limited to the built-in
  echo \# GMRES solver
  echo \#
  echo \# You can install PETSc by downloading it from here:
  echo \#
  echo \# https://petsc.org/
  echo \#
  echo \#################################################################
fi 				    

PKG_CHECK_MODULES([GLIB], [glib-2.0])

AX_GCC_X86_CPU_SUPPORTS(avx, hasavx=yes) 
AX_GCC_X86_CPU_SUPPORTS(avx2, hasavx2=yes) 
AX_GCC_X86_CPU_SUPPORTS(fma, hasfma=yes) 

if test x$hasavx = xyes; then
 SIMD_FLAGS="$SIMD_FLAGS -mavx -DWBFMM_USE_AVX"
fi
if test x$hasfma = xyes; then
 SIMD_FLAGS="$SIMD_FLAGS -mfma -DWBFMM_USE_AVX"
fi

AC_CHECK_LIB(m, cos)
AC_CONFIG_MACRO_DIR([m4])

AC_OPENMP
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN
AC_F77_WRAPPERS

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

CFLAGS="$CFLAGS $OPENMP_CFLAGS $SIMD_FLAGS $GLIB_CFLAGS $blaswrap_cflags $wbfmm_cflags $sqt_cflags $agg_cflags"
LIBS="$GLIB_LIBS $blaswrap_libs $LAPACK_LIBS $LIBS $FLIBS"

AC_SUBST(petsc_real)
##AC_SUBST(CPPFLAGS)
AC_SUBST(sqt_submodule)
AC_SUBST(sqt_include)
AC_SUBST(sqt_libs)
AC_SUBST(wbfmm_submodule)
AC_SUBST(wbfmm_include)
AC_SUBST(wbfmm_libs)
AC_SUBST(agg_submodule)
AC_SUBST(agg_include)
AC_SUBST(agg_libs)
AC_SUBST(subdirs)

AC_SUBST(CFLAGS)
AC_SUBST(LIBS)
AC_SUBST(LDFLAGS)

AM_CONDITIONAL(PETSC, test x$PETSC = xyes)
AM_CONDITIONAL(SUBMODULES, test x$submodules = xyes)

AC_CONFIG_FILES([
Makefile
src/Makefile
src/nbi-config
tools/Makefile
doc/Makefile
doc/html/Makefile])

AC_OUTPUT

