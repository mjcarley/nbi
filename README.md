NBI is a library and codes for the Nystrom Boundary Integral method.
Boundary integral methods are a popular approach to solving potential
problems, in particular for the Laplace and Helmholtz problems in such
areas as fluid dynamics, acoustics, and electromagnetism. They are a
natural choice for the solution of problems in unbounded domains, such
as wave scattering, where the radiation boundary condition is
automatically satisfied by the nature of the formulation. The Nystrom
method is one approach to the solution of boundary integral problems,
which lends itself to acceleration using the Fast Multipole Method
(FMM). NBI is a library for the solution of boundary integral problems
based on the approach of [Greengard et
al.](https://dx.doi.org/10.1016/j.jcpx.2021.100092) and the FMM
methods of Gumerov and Duraiswami
([2003](https://dx.doi.org/10.1137/S1064827501399705),
[2004](https://dx.doi.org/10.1016/b978-0-08-044371-3.x5000-5),
[2005](https://www.academia.edu/2791247/Comparison_of_the_efficiency_of_translation_operators_used_in_the_fast_multipole_method_for_the_3D_Laplace_equation),
[2009](https://dx.doi.org/10.1121/1.3021297)). The code includes a
number of executables which can be used to set up and solve problems
on realistic geometries, with a number of examples provided for
testing of the solver. Results can be visualized using
[GMSH](https://gmsh.info), a standard free meshing program.

NBI contains the tinyexpr parser and evaluation library written by
Lewis Van Winkle and released under the zlib licence. It is maintained
here:

https://github.com/codeplea/tinyexpr

# Getting NBI

NBI can be downloaded from https://github.com/mjcarley/nbi, but it is
recommended (see under Installation) to clone it using

`git clone --recursive https://github.com/mjcarley/nbi`

in order to facilitate installation of dependencies as submodules.

Alternatively, a Docker image of the current version is available from
https://github.com/mjcarley/nbi/pkgs/container/nbi More details on
using the Docker image are given below.

# Prerequisites

NBI requires the following packages:
- BLAS wrapper https://github.com/mjcarley/blaswrap
- SQT (Singular Quadrature for Triangles) https://github.com/mjcarley/sqt
- WBFMM (Wide Band Fast Multipole Method) https://github.com/mjcarley/wbfmm
- GMSH API https://www.gmsh.info/

NBI can also make use of:
- AGG (Aerodynamic Geometry Generator) https://github.com/mjcarley/agg,
  which requires the Triangle API https://github.com/wo80/Triangle/
- the iterative solvers in PETSc https://petsc.org/ with a version
  number greater than or equal to 3.17

AGG provides a geometry specification interface based on the
methods of Brenda Kulfan https://www.brendakulfan.com/research, which
can be used to generate surfaces such as wings found in aeronautical
applications. 

It is recommended to install PETSc to provide a choice of
high-performance iterative solvers for the boundary integral
equations. NBI does have a built-in basic GMRES solver but PETSc
provides access to a range of methods which are not otherwise available.
PETSc version 3.17 or higher is required. 

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

and make sure that configuration files for pkg-config are installed on
a default search path. It is recommended to install a version of PETSc
*without* MPI support, using the `--with-mpi=0` option to the PETSc
configure. NBI is a single-processor (threaded) code and using a
sequential version of PETSc avoids problems with linking to MPI
libraries whose features will not be used. 

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
scripts `test-all-laplace` and `test-all-helmholtz`. Further details are
given in `.../examples/README.md`

# Running NBI using the Docker image

It is recommended to compile and install NBI on your system, to take
advantage of properly optimized versions of libraries such as BLAS and
LAPACK, but a Docker image is available to allow some testing and
experimentation. It can be installed from
https://github.com/mjcarley/nbi/pkgs/container/nbi and the Dockerfile
used to generate it is distributed with the NBI source code.

Once installed, the Docker image contains results for two sample
calculations, one Laplace and one Helmholtz, which are generated using
the scripts `docker-laplace` and `docker-helmholtz` respectively, in
the directory `.../examples`. 

To install the Docker image:

`docker pull ghcr.io/mjcarley/nbi:main`

and to start an interactive shell which will allow you to run the NBI
codes:

`docker run -it -d [image name]`

`docker exec -i [container name] bash`

where `[container name]` can be found using `docker ps`.

In the newly running interactive shell, the full set of examples can
be run using

`cd /nbi/examples`

`./test-all-laplace`

`./test-all-helmholtz`

and results can be extracted to the host using

`docker cp [container name]:/nbi/examples .`

which will copy the full `examples` directory tree with all results
and visualization files.
