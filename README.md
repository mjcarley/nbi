NBI is a library and codes for the Nystrom Boundary Integral method.

It contains the tinyexpr parser and evaluation library written by
Lewis Van Winkle and released under the zlib licence. It is maintained
here:

https://github.com/codeplea/tinyexpr

# Getting NBI

NBI can be downloaded from https://github.com/mjcarley/nbi, but it is
recommended (see under Installation) to clone it using

`git clone --recursive https://github.com/mjcarley/nbi`

in order to facilitate installation of dependencies as submodules.

# Prerequisites

NBI requires the following packages:
- BLAS wrapper https://github.com/mjcarley/blaswrap
- SQT (Singular Quadrature for Triangles) https://github.com/mjcarley/sqt
- WBFMM (Wide Band Fast Multipole Method) https://github.com/mjcarley/wbfmm
- GMSH API https://www.gmsh.info/

NBI can also make use of:
- AGG (Aerodynamic Geometry Generator) https://github.com/mjcarley/agg,
  which requires the Triangle API https://github.com/wo80/Triangle/
- the iterative solvers in PETSc https://petsc.org/

AGG provides a geometry specification interface based on the
methods of Brenda Kulfan https://www.brendakulfan.com/research, which
can be used to generate surfaces such as wings found in aeronautical
applications. 

It is recommended to install PETSc to provide a choice of
high-performance iterative solvers for the boundary integral
equations. NBI does have a built-in basic GMRES solver but PETSc
provides access to a range of methods which are not otherwise available.

# Installation

You need to install the BLAS wrapper separately before installing NBI
or its dependencies, as it is used by a number of different
dependencies of NBI.

You will also need to install the GMSH API library, available from

- https://www.gmsh.info/

If you are installing the AGG geometry library (this is done by
default if you clone the repository from github) you will need to
install Christian Woltering's API for Jonathan Shewchuk's Triangle
code, available from:

- https://github.com/wo80/Triangle/

You should install PETSc if you want to use its solvers in place of
the (very) basic built-in iterative solver:

- https://petsc.org/

and make sure the environment variable PETSC_DIR is set and visible in
the shell.

To install NBI, it is recommended to clone the NBI repository using

`git clone --recursive https://github.com/mjcarley/nbi`

This will download the NBI repository including SQT, WBFMM, and AGG as
submodules. In the top directory of the repository, generate the
configuration scripts using

`. autogen.sh`

To configure the code:

`./configure [OPTIONS]`

You may need to set some variables to tell `configure` how to set
certain flags, such as the location of libraries and headers.  You can
set these variables when you call `configure`. For example:

`CPPFLAGS=-I/.../include LDFLAGS=-L/.../lib ./configure`

will set additional preprocessor and linker flags for headers or
libraries which are not in standard locations.

To build and install the library and codes:

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
