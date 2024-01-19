/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(name,NAME) name ## _

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* aerodynamic geometry library available */
#define HAVE_AGG 1

/* Define if avx2 instructions are supported */
#define HAVE_AVX2_INSTRUCTIONS 1

/* Define if avx instructions are supported */
#define HAVE_AVX_INSTRUCTIONS 1

/* BLAS wrapper header available */
#define HAVE_BLASWRAP 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define if fma instructions are supported */
#define HAVE_FMA_INSTRUCTIONS 1

/* Define to 1 if you have the <gmshc.h> header file. */
#define HAVE_GMSHC_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `gmsh' library (-lgmsh). */
#define HAVE_LIBGMSH 1

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* PETSC solvers available */
#define HAVE_PETSC 1

/* singular quadrature library available */
#define HAVE_SQT 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* wide band FMM library available */
#define HAVE_WBFMM 1

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Binary age */
#define NBI_BINARY_AGE 0

/* Interface age */
#define NBI_INTERFACE_AGE 0

/* Major version */
#define NBI_MAJOR_VERSION 1

/* Micro version */
#define NBI_MICRO_VERSION 0

/* Minor version */
#define NBI_MINOR_VERSION 0

/* Name of package */
#define PACKAGE "nbi"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#define PACKAGE_NAME "nbi"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "nbi 1.0.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "nbi"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.0.0"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "1.0.0"
