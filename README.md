NBI is a library and codes for the Nystrom Boundary Integral method.

It contains the tinyexpr parser and evaluation library written by
Lewis Van Winkle and released under the zlib licence. It is maintained
here:

https://github.com/codeplea/tinyexpr

# Prerequisites

To install NBI, you will need to install some or all of the following:

- BLAS wrapper https://github.com/mjcarley/blaswrap
- Singular Quadrature for Triangles https://github.com/mjcarley/sqt
- Wide Band Fast Multipole Method https://github.com/mjcarley/wbfmm
- Aerodynamic Geometry Generator https://github.com/mjcarley/agg
- GMSH API https://www.gmsh.info/
- Triangle API https://github.com/wo80/Triangle/
- PETSc https://petsc.org/

NBI and a number of its dependencies use a wrapper to interface with
BLAS, which you need to install first. This is available from:

- https://github.com/mjcarley/blaswrap

You will also need to have installed the following libraries:

- https://github.com/mjcarley/sqt
- https://github.com/mjcarley/wbfmm
- https://github.com/mjcarley/agg

If you clone the github repository using

`git clone --recursive https://github.com/mjcarley/nbi`

these libraries will be downloaded automatically in the `externals`
subdirectory of `nbi`, and will be built as part of the NBI build
process. 

You will also need the GMSH API library, available from

- https://www.gmsh.info/

If you are installing the AGG geometry library from within this tree,
you will need to install Christian Woltering's API for Jonathan
Shewchuk's Triangle code, available from:

- https://github.com/wo80/Triangle/

You should install PETSc if you want to use it in place of the (very)
basic built-in iterative solver (recommended):

- https://petsc.org/

and make sure the environment variable PETSC_DIR is set and visible in
the shell.

# Installation

If you have downloaded the source from github or equivalent, you will
need the autotools suite. Generate the configure script with

`. autogen.sh`

To configure and install the code,

`./configure [OPTIONS]`

`make`

`make install`

For information on options to control configuration, including the
installation location and where to find any required libraries:

  `./configure --help`

# Using and testing NBI

The main documentation, including details of examples, is in the
`.../doc` subdirectory of the distribution. The `.../examples`
subdirectory contains a number of examples which can be run using the
scripts test-all-laplace and test-all-helmholtz. Further details are
given in `.../examples/README.md`
